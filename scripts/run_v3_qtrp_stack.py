"""v3: fit the designs on the QTRP persulfidation label, per the registered plan.

Everything this script does is fixed in advance by
`protocols/v3_qtrp_analysis_plan.json` (status ``registered``, sha256 pinned in
`results/v3_plan_registration.json`).  Nothing here may be changed after a number
is seen; that is the point of the registration.

**Cohorts.** The primary cohort `v3_pooled_ph5_raw` pools the two pH 5.0 arms of
Fu 2020: a site is positive when a persulfidated form was detected in either arm,
and negative when a free thiol was detected in a pH 5.0 arm and no persulfidated
form was detected in either pH 5.0 arm.  Sites called SSH only in OTHER human
arms stay negative - that retention is exactly what distinguishes this cohort
from the demoted `arm_clean` rule, which emptied depth strata 5 to 7 and pushed
39% of positives into single-class strata
(`reports/QTRP_DEPTH_CONFOUND.md`, `reports/QTRP_TRAINED_DEPTH_CONDITIONED.md`).
`v3_pooled_ph5_armclean` and `v3_pooled_ph5_fullclean` are sensitivity cohorts.
pH 7.6 arms and S3-S5 never enter training; they only exclude negatives.

**Folds are grouped by homology component, not by protein**
(`results/v3_human_homology_components.csv`, 7-mer containment >= 0.40).  Sites on
accessions with no component are dropped and counted.  Outer: 5 component-grouped
folds.  Inner: 5 component-grouped folds inside each outer training set, used to
select, from the predeclared grid below, the HistGradientBoosting configuration
and the dimension-reduction option.  The inner objective is plain AUC on the
inner test folds.  Outer test folds enter nothing.

Predeclared grid: three configurations - the project default used by every prior
round, a shallower and more regularised variant, and a slower-learning variant
with more iterations - crossed with dimension reduction in {full, mutual-information
top 64, top 128, top 256} for designs wider than 200 columns and {full} otherwise.
Class imbalance is handled by ONE declared option, the project default of no class
weighting, so that v3 stays comparable with every earlier round.  Platt calibration
is not applied because it is monotone and cannot change a rank statistic; that is a
declared omission, not an oversight.

**Primary estimand.** Depth-stratified concordance: within-stratum AUC weighted by
n_positive * n_negative, strata being the exact number of distinct human QTRP arms
in which the site appears in either adduct form.  Depth is counted over ARMS, never
over rows: every persulfide site is also detected as a free thiol, so a row count
encodes the label and reached AUC 0.9478 in the discarded first version of that
probe.  Both calibers required by the plan are computed - global, and within
homology component crossed with depth.  The within-component caliber is reported as
`not_identifiable` unless at least 30 components contribute an informative stratum,
because the human QTRP cohort has little component redundancy; that threshold is
declared here, before the run.

**Predeclared differences** (on the stratified statistic): chem1046 - struct53,
struct53_plus_chem1046 - chem1046, chem1046 - v2_chem_zeroshot, chem1046 - depth1,
struct53 - probe_visibility10, and v2_chem_zeroshot - v2_full_zeroshot - the last of
which the depth-conditioned round omitted and which has so far been measured only on
cohorts carrying the arm_clean mechanism.

**Depth control, all five items from the plan.** Depth never enters a fitted design.
The primary metric is stratified.  A depth-matched paired subcohort is drawn at 1:2
negatives per positive within stratum, seed 20260915, and scored with plain AUC.
The Spearman correlation between each score and depth is reported: if a v3 design
correlates with depth more strongly than the v2 chemistry zero-shot score does, the
model absorbed the confound and must not be deployed whatever its AUC.  A
within-component label permutation control is run for chem1046 and must return to
about 0.5.

**Intervals.** Cluster bootstrap over homology components, 5000 replicates, seed
20260915, statistic recomputed in every replicate.

Reading rules: lower bound above 0.5 (or above 0 for a difference) discriminates;
an interval containing it is undetermined at this sample size; an upper bound below
it is inverted.

Outputs are written per cohort as they finish, so a failure late in the run still
leaves the primary cohort's result on disk.
"""
from __future__ import annotations

import collections
import json
import platform
import sys
import time

import numpy as np
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.feature_selection import mutual_info_classif
from sklearn.model_selection import GroupKFold

from common import RESULTS, ROOT, read_csv, sha256, write_csv, write_json
import v2_apply
from run_qtrp_persulfidation_test import (
    SITES, load_proteome, probe_matrix, structural_matrix,
)
from run_qtrp_negative_cleaning_test import COMPILED
from run_qtrp_trained_model import human_matrix
from run_qtrp_depth_conditioned import fast_auc, stratified_auc, verdict
from run_cross_protease_detectability_probe import (
    FEATURE_NAMES, HGB_KWARGS, VISIBILITY_ONLY,
)
from v2_stack import FEATURES

SCRIPT = ROOT / "scripts" / "run_v3_qtrp_stack.py"
PLAN = ROOT / "protocols" / "v3_qtrp_analysis_plan.json"
REGISTRATION = RESULTS / "v3_plan_registration.json"
COMPONENTS = RESULTS / "v3_human_homology_components.csv"
PH5_ARMS = {("S1", "5.0"), ("S2", "5.0")}
REPLICATES = 5000
SEED = 20260915
OUTER_FOLDS = INNER_FOLDS = 5
WIDE = 200
MIN_COMPONENTS_FOR_WITHIN_CALIBER = 30
MATCH_RATIO = 2

CONFIGS = (
    ("project_default", dict(HGB_KWARGS)),
    ("shallower_regularised", {**HGB_KWARGS, "max_leaf_nodes": 15,
                               "min_samples_leaf": 40, "l2_regularization": 1.0}),
    ("slower_learning", {**HGB_KWARGS, "learning_rate": 0.03, "max_iter": 600}),
)
DIFFS = (
    ("chem1046", "struct53"),
    ("struct53_plus_chem1046", "chem1046"),
    ("chem1046", "v2_chem_zeroshot"),
    ("chem1046", "depth1"),
    ("struct53", "probe_visibility10"),
    ("v2_chem_zeroshot", "v2_full_zeroshot"),
)


def pooled_labels(rows, compiled_ssh):
    """Positives, the three negative definitions, and the all-arm depth column."""
    seen_arms = collections.defaultdict(set)
    ph5_ssh, ph5_sh, ssh_any = set(), set(), set()
    for row in rows:
        if row.get("species") != "human":
            continue
        label = row["class_from_adduct"] or row["class_from_sheet"]
        if label not in ("SH", "SSH"):
            continue
        key = (row["accession"], int(row["site"]))
        arm = (row["source_table"], row["ph"])
        seen_arms[key].add(arm)
        if label == "SSH":
            ssh_any.add(key)
            if arm in PH5_ARMS:
                ph5_ssh.add(key)
        elif arm in PH5_ARMS:
            ph5_sh.add(key)
    depth = {key: len(arms) for key, arms in seen_arms.items()}
    raw = ph5_sh - ph5_ssh
    return {
        "v3_pooled_ph5_raw": (ph5_ssh, raw),
        "v3_pooled_ph5_armclean": (ph5_ssh, raw - ssh_any),
        "v3_pooled_ph5_fullclean": (ph5_ssh, (raw - ssh_any) - compiled_ssh),
    }, depth


def dimension_options(width):
    if width <= WIDE:
        return [("full", None)]
    return [("full", None)] + [(f"mi_top{k}", k) for k in (64, 128, 256)]


def mi_ranking(x_train, y_train):
    """Feature order by mutual information, on a median-imputed copy.

    The chemistry block carries NaN by construction: five frozen structural columns
    have no value for human sequences, and every earlier round in this tree passed
    those NaN straight to HistGradientBoosting, which handles them natively.
    `mutual_info_classif` refuses them, so the RANKING is computed on a copy whose
    NaN are replaced by the column median of the same training subset (a column that
    is entirely NaN gets zeros and therefore no mutual information).  The imputation
    exists only to order the columns; the model itself is always fitted on the raw
    matrix, NaN included, so no imputed value ever reaches a prediction.
    """
    filled = np.array(x_train, dtype=np.float64, copy=True)
    missing = ~np.isfinite(filled)
    if missing.any():
        for column in np.flatnonzero(missing.any(axis=0)):
            values = filled[:, column]
            finite = values[np.isfinite(values)]
            values[~np.isfinite(values)] = float(np.median(finite)) if finite.size else 0.0
    return np.argsort(-mutual_info_classif(filled, y_train, random_state=0))


def fit_predict(x_train, y_train, x_test, kwargs):
    model = HistGradientBoostingClassifier(**kwargs)
    model.fit(x_train, y_train)
    return model.predict_proba(x_test)[:, 1]


def select_and_score(x, y, groups, design_name, log):
    """Component-grouped outer folds with component-grouped inner selection."""
    scores = np.full(len(y), np.nan)
    options = dimension_options(x.shape[1])
    chosen = []
    outer = GroupKFold(n_splits=OUTER_FOLDS)
    for fold, (train, test) in enumerate(outer.split(x, y, groups), start=1):
        best = None
        inner_groups = groups[train]
        n_inner = min(INNER_FOLDS, len(np.unique(inner_groups)))
        splits = list(GroupKFold(n_splits=n_inner).split(x[train], y[train], inner_groups))
        # the mutual-information ranking depends on the inner training subset only,
        # so it is computed once per split and reused across k and across configurations
        ranking = {}
        if any(top_k is not None for _, top_k in options):
            for position, (inner_train, _) in enumerate(splits):
                a = train[inner_train]
                ranking[position] = mi_ranking(x[a], y[a])
        for option_name, top_k in options:
            for config_name, kwargs in CONFIGS:
                predictions = np.full(len(train), np.nan)
                for position, (inner_train, inner_test) in enumerate(splits):
                    a, b = train[inner_train], train[inner_test]
                    if y[a].sum() == 0 or y[a].sum() == len(a):
                        continue
                    cols = slice(None) if top_k is None else ranking[position][:top_k]
                    predictions[inner_test] = fit_predict(
                        x[a][:, cols], y[a], x[b][:, cols], kwargs)
                usable = np.isfinite(predictions)
                if usable.sum() < 10 or y[train][usable].sum() == 0:
                    continue
                objective = fast_auc(y[train][usable], predictions[usable])
                if best is None or objective > best[0]:
                    best = (objective, option_name, top_k, config_name, kwargs)
        if best is None:
            raise RuntimeError(f"inner selection failed: {design_name} fold {fold}")
        _, option_name, top_k, config_name, kwargs = best
        cols = slice(None) if top_k is None else mi_ranking(x[train], y[train])[:top_k]
        scores[test] = fit_predict(x[train][:, cols], y[train], x[test][:, cols], kwargs)
        chosen.append({"fold": fold, "inner_objective": round(best[0], 4),
                       "dimension": option_name, "config": config_name})
        log(f"    {design_name} fold {fold}: {option_name} / {config_name} "
            f"(inner AUC {best[0]:.4f})")
    return scores, chosen


def within_component_stratified(y, scores, strata, components):
    """The second caliber: strata are (component, depth) pairs."""
    keys = np.asarray([f"{c}|{d}" for c, d in zip(components, strata)])
    contributing = set()
    total_weight = total = 0.0
    for value in np.unique(keys):
        mask = keys == value
        n1 = int(y[mask].sum())
        n0 = int(mask.sum()) - n1
        if n1 == 0 or n0 == 0:
            continue
        contributing.add(value.split("|")[0])
        weight = float(n1 * n0)
        total += weight * fast_auc(y[mask], scores[mask])
        total_weight += weight
    if total_weight == 0:
        return float("nan"), 0
    return total / total_weight, len(contributing)


def depth_matched_subcohort(y, strata, seed=SEED, ratio=MATCH_RATIO):
    rng = np.random.default_rng(seed)
    keep = []
    for value in np.unique(strata):
        mask = np.flatnonzero(strata == value)
        positives = mask[y[mask] == 1]
        negatives = mask[y[mask] == 0]
        if positives.size == 0 or negatives.size == 0:
            continue
        take = min(negatives.size, ratio * positives.size)
        keep.extend(positives.tolist())
        keep.extend(rng.choice(negatives, size=take, replace=False).tolist())
    return np.asarray(sorted(keep))


def cluster_bootstrap(y, components, strata, score_map, diffs):
    unique = np.unique(components)
    index_of = {c: np.flatnonzero(components == c) for c in unique}
    names = list(score_map)
    draws = {f"{n}|{k}": [] for n in names for k in ("plain", "stratified")}
    diff_draws = {f"{a}__minus__{b}": [] for a, b in diffs}
    rng = np.random.default_rng(SEED)
    skipped = 0
    for _ in range(REPLICATES):
        picked = rng.integers(0, len(unique), len(unique))
        index = np.concatenate([index_of[unique[i]] for i in picked])
        yy, ss = y[index], strata[index]
        if yy.sum() == 0 or yy.sum() == len(yy):
            skipped += 1
            continue
        values = {}
        for name in names:
            vector = score_map[name][index]
            draws[f"{name}|plain"].append(fast_auc(yy, vector))
            value = stratified_auc(yy, vector, ss)
            values[name] = value
            draws[f"{name}|stratified"].append(value)
        if any(not np.isfinite(v) for v in values.values()):
            skipped += 1
            continue
        for left, right in diffs:
            diff_draws[f"{left}__minus__{right}"].append(values[left] - values[right])

    def interval(sample):
        array = np.asarray([v for v in sample if np.isfinite(v)], dtype=float)
        if array.size == 0:
            return [float("nan"), float("nan")]
        return [float(np.quantile(array, 0.025)), float(np.quantile(array, 0.975))]

    return {
        "replicates": REPLICATES, "usable_replicates": REPLICATES - skipped, "seed": SEED,
        "unit": "homology_component",
        "estimand": "stratified concordance conditional on exact detection depth, fixed out-of-fold scores",
        "auc": {n: {"plain": fast_auc(y, score_map[n]),
                    "plain_interval": interval(draws[f"{n}|plain"]),
                    "stratified": stratified_auc(y, score_map[n], strata),
                    "stratified_interval": interval(draws[f"{n}|stratified"])} for n in names},
        "stratified_differences": {k: {"mean_difference": float(np.mean(v)) if v else float("nan"),
                                       "percentile_95_interval": interval(v)}
                                   for k, v in diff_draws.items()},
    }


def main():
    started = time.time()
    def log(message):
        print(message, flush=True)

    plan = json.loads(PLAN.read_text())
    registration = json.loads(REGISTRATION.read_text())
    if sha256(PLAN) != registration["plan_sha256"]:
        raise RuntimeError("the registered plan on disk does not match its registration record")
    if plan["status"] != "registered":
        raise RuntimeError(f"plan status is {plan['status']!r}, refusing to run")

    component_of = {r["accession"]: r["component"] for r in read_csv(COMPONENTS)}
    rows = read_csv(SITES)
    sequences = load_proteome()
    compiled_ssh = {(r["accession"], int(r["site"]))
                    for r in read_csv(COMPILED) if r["compiled_ssh"] == "1"}
    cohorts, depth = pooled_labels(rows, compiled_ssh)
    visibility_columns = [i for i, n in enumerate(FEATURE_NAMES) if n in VISIBILITY_ONLY]

    union_keys = sorted({k for pair in cohorts.values() for k in pair[0] | pair[1]})
    residue_ok = np.asarray([
        bool(sequences.get(a)) and 1 <= p <= len(sequences[a]) and sequences[a][p - 1] == "C"
        for a, p in union_keys])
    S_all, _, structure_stats = structural_matrix(union_keys, sequences)
    P_all, probe_ok = probe_matrix(union_keys, sequences)
    covered = np.isfinite(S_all).any(axis=1) & residue_ok & probe_ok
    keys = [k for k, keep in zip(union_keys, covered) if keep]
    no_component = [k for k in keys if k[0] not in component_of]
    keys = [k for k in keys if k[0] in component_of]
    index_map = {k: i for i, k in enumerate(union_keys)}
    rows_kept = np.asarray([index_map[k] for k in keys])
    S_all, P_all = S_all[rows_kept], P_all[rows_kept]
    V_all, columns, zero_shot, zero_shot_audit = human_matrix(keys, sequences)
    row_of = {k: i for i, k in enumerate(keys)}
    log(json.dumps({"union": len(union_keys), "usable": len(keys),
                    "dropped_no_component": len(no_component),
                    "v2_matrix": list(V_all.shape), **structure_stats,
                    "seconds": round(time.time() - started, 1)}))

    summary, strata_rows, intervals, audit_cohorts = [], [], {}, {}
    for cohort, (positives, negatives) in cohorts.items():
        selected = [k for k in keys if k in positives or k in negatives]
        y = np.asarray([1 if k in positives else 0 for k in selected])
        if y.sum() < 20 or (y == 0).sum() < 20:
            audit_cohorts[cohort] = {"skipped": "too few usable sites"}
            continue
        index = np.asarray([row_of[k] for k in selected])
        comps = np.asarray([component_of[a] for a, _ in selected])
        proteins = np.asarray([a for a, _ in selected])
        strata = np.asarray([depth[k] for k in selected])
        log(f"  {cohort}: {len(selected)} sites, {int(y.sum())} positive, "
            f"{len(set(comps))} components, {len(set(proteins))} proteins")

        designs = {
            "struct53": S_all[index],
            "chem1046": V_all[np.ix_(index, columns["chem"])],
            "struct53_plus_chem1046": np.hstack([S_all[index],
                                                 V_all[np.ix_(index, columns["chem"])]]),
            "struct53_plus_zeroshot": np.hstack([
                S_all[index],
                zero_shot["v2_chem_zeroshot"][index].reshape(-1, 1)]),
            "probe_digest25": P_all[index],
            "probe_visibility10": P_all[np.ix_(index, visibility_columns)],
        }
        score_map, selection = {}, {}
        for name, x in designs.items():
            x = np.asarray(x, dtype=np.float64)
            score_map[name], selection[name] = select_and_score(x, y, comps, name, log)
            if not np.isfinite(score_map[name]).all():
                raise RuntimeError(f"non-finite out-of-fold scores: {cohort} {name}")
        score_map["v2_chem_zeroshot"] = zero_shot["v2_chem_zeroshot"][index]
        score_map["v2_full_zeroshot"] = zero_shot["v2_full_zeroshot"][index]
        score_map["depth1"] = strata.astype(float)

        dropped_positives = dropped_negatives = 0
        for value in np.unique(strata):
            mask = strata == value
            n1 = int(y[mask].sum())
            n0 = int(mask.sum()) - n1
            strata_rows.append({"cohort": cohort, "depth": int(value), "n_sites": int(mask.sum()),
                                "n_positive": n1, "n_negative": n0,
                                "positive_rate": round(n1 / max(1, int(mask.sum())), 4),
                                "informative": bool(n1 and n0)})
            if not (n1 and n0):
                dropped_positives += n1
                dropped_negatives += n0

        boot = cluster_bootstrap(y, comps, strata, score_map, DIFFS)
        matched = depth_matched_subcohort(y, strata)
        permuted = None
        if cohort == "v3_pooled_ph5_raw":
            rng = np.random.default_rng(SEED)
            y_perm = y.copy()
            for component in np.unique(comps):
                mask = np.flatnonzero(comps == component)
                y_perm[mask] = rng.permutation(y[mask])
            x = np.asarray(designs["chem1046"], dtype=np.float64)
            perm_scores, _ = select_and_score(x, y_perm, comps, "chem1046_permuted", log)
            permuted = {"design": "chem1046", "unit": "homology component",
                        "plain": round(fast_auc(y_perm, perm_scores), 4),
                        "stratified": round(stratified_auc(y_perm, perm_scores, strata), 4)}

        for name in score_map:
            within, contributing = within_component_stratified(y, score_map[name], strata, comps)
            entry = boot["auc"][name]
            summary.append({
                "cohort": cohort, "design": name,
                "n_sites": len(selected), "n_positive": int(y.sum()),
                "n_components": len(set(comps)), "n_proteins": len(set(proteins)),
                "auc_plain": round(entry["plain"], 4),
                "auc_depth_stratified": round(entry["stratified"], 4),
                "ci_low": round(entry["stratified_interval"][0], 4),
                "ci_high": round(entry["stratified_interval"][1], 4),
                "verdict_stratified": verdict(entry["stratified_interval"]),
                "plain_minus_stratified": round(entry["plain"] - entry["stratified"], 4),
                "auc_within_component_x_depth": ("not_identifiable"
                                                 if contributing < MIN_COMPONENTS_FOR_WITHIN_CALIBER
                                                 else round(within, 4)),
                "components_contributing_to_within_caliber": contributing,
                "auc_depth_matched_1_to_2": round(fast_auc(y[matched], score_map[name][matched]), 4),
                "spearman_score_vs_depth": round(float(spearmanr(score_map[name], strata).statistic), 4),
                "selection": json.dumps(selection.get(name, "not_fitted")),
            })

        intervals[cohort] = boot
        audit_cohorts[cohort] = {
            "n_sites": len(selected), "n_positive": int(y.sum()),
            "n_components": len(set(comps)), "n_proteins": len(set(proteins)),
            "depth_values": [int(v) for v in np.unique(strata)],
            "informative_strata": sum(1 for r in strata_rows
                                      if r["cohort"] == cohort and r["informative"]),
            "positives_in_uninformative_strata": dropped_positives,
            "negatives_in_uninformative_strata": dropped_negatives,
            "mean_depth_positive": round(float(strata[y == 1].mean()), 3),
            "mean_depth_negative": round(float(strata[y == 0].mean()), 3),
            "depth_matched_subcohort_sites": int(matched.size),
            "depth_matched_subcohort_positives": int(y[matched].sum()),
            "permutation_control": permuted,
            "inner_selection": selection,
        }
        write_csv(RESULTS / "v3_qtrp_stack_summary.csv", summary)
        write_csv(RESULTS / "v3_qtrp_stack_depth_strata.csv", strata_rows)
        write_json(RESULTS / "v3_qtrp_stack_intervals.json", intervals)
        log(f"  {cohort} done at {round((time.time() - started) / 60, 1)} min")

    primary = "v3_pooled_ph5_raw"
    criteria = {}
    if primary in intervals:
        entry = intervals[primary]["auc"]["chem1046"]
        diffs = intervals[primary]["stratified_differences"]
        chem_rho = abs(next(r["spearman_score_vs_depth"] for r in summary
                            if r["cohort"] == primary and r["design"] == "chem1046"))
        zero_rho = abs(next(r["spearman_score_vs_depth"] for r in summary
                            if r["cohort"] == primary and r["design"] == "v2_chem_zeroshot"))
        criteria = {
            "1_stratified_auc_lower_bound_above_0_65": bool(entry["stratified_interval"][0] > 0.65),
            "1_value": round(entry["stratified_interval"][0], 4),
            "2_chem1046_minus_v2_chem_zeroshot_above_0": bool(
                diffs["chem1046__minus__v2_chem_zeroshot"]["percentile_95_interval"][0] > 0),
            "2_value": [round(v, 4) for v in
                        diffs["chem1046__minus__v2_chem_zeroshot"]["percentile_95_interval"]],
            "3_chem1046_minus_depth1_above_0_and_spearman_not_higher": bool(
                diffs["chem1046__minus__depth1"]["percentile_95_interval"][0] > 0
                and chem_rho <= zero_rho),
            "3_values": {"difference_interval": [round(v, 4) for v in
                                                 diffs["chem1046__minus__depth1"]["percentile_95_interval"]],
                         "chem1046_abs_spearman_vs_depth": round(chem_rho, 4),
                         "v2_chem_zeroshot_abs_spearman_vs_depth": round(zero_rho, 4)},
            "4_external_qpers_sid_tier_a": "not_evaluated_in_this_round",
            "all_four_met": False,
        }
        criteria["all_four_met"] = bool(
            criteria["1_stratified_auc_lower_bound_above_0_65"]
            and criteria["2_chem1046_minus_v2_chem_zeroshot_above_0"]
            and criteria["3_chem1046_minus_depth1_above_0_and_spearman_not_higher"]
            and criteria["4_external_qpers_sid_tier_a"] is True)

    write_json(RESULTS / "v3_qtrp_stack_run_audit.json", {
        "completed": True,
        "elapsed_minutes": round((time.time() - started) / 60, 2),
        "host": platform.node(),
        "registered_plan": {"file": "protocols/v3_qtrp_analysis_plan.json",
                            "sha256": sha256(PLAN), "status": plan["status"],
                            "registered_utc": plan["registered_utc"],
                            "verified_against": "results/v3_plan_registration.json"},
        "cohorts": audit_cohorts,
        "sites_dropped_without_a_homology_component": len(no_component),
        "structure_coverage": structure_stats,
        "zero_shot_members": zero_shot_audit,
        "grid": {"configs": [name for name, _ in CONFIGS],
                 "config_kwargs": {name: kwargs for name, kwargs in CONFIGS},
                 "dimension_options_for_wide_designs": ["full", "mi_top64", "mi_top128", "mi_top256"],
                 "wide_threshold_columns": WIDE,
                 "class_weighting": "none, the project default, declared as the single option",
                 "mutual_information_ranking_on_imputed_copy": (
                     "the chemistry block carries NaN by construction (five frozen structural columns have no "
                     "value for human sequences). mutual_info_classif refuses NaN, so the RANKING is computed on "
                     "a copy whose NaN are replaced by the column median of the same training subset; an "
                     "all-NaN column gets zeros and no mutual information. Models are always fitted on the raw "
                     "matrix with NaN intact, as every earlier round in this tree did, so no imputed value "
                     "reaches a prediction."),
                 "platt_calibration": "not applied; monotone and cannot change a rank statistic"},
        "folds": {"outer": OUTER_FOLDS, "inner": INNER_FOLDS, "unit": "homology component",
                  "grouping_file": "results/v3_human_homology_components.csv"},
        "within_component_caliber_rule": (
            f"reported only when at least {MIN_COMPONENTS_FOR_WITHIN_CALIBER} components contribute an "
            "informative stratum; otherwise not_identifiable"),
        "depth_matched_subcohort": {"ratio": f"1:{MATCH_RATIO}", "seed": SEED},
        "deployment_criteria": criteria,
        "deployment_decision": (
            "criterion 4 requires the qPerS-SID tier A external evaluation, which is a separate one-time "
            "unlock and was not run here, so no deployment decision is made and the deployed weights "
            "0.64/0.20/0.16 stay unchanged"),
        "limits": [
            "314-positive scale against 1,046 features; component-grouped folds measure overfitting but do not remove it",
            "depth stratification cannot separate chemistry from occupancy: a site repeatedly seen as a persulfide may really carry more persulfide",
            "all positives come from one paper, one cell system, NaHS-loaded lysate at pH 5.0, which is not physiological persulfidation",
            "negatives remain 'no persulfidated form seen in this paper', not chemically confirmed unmodified",
            "the homology grouping is a declared grouping rule",
        ],
        "input_hashes": {
            "results/qtrp_sites_normalised.csv": sha256(SITES),
            "results/frbm_compiled_sites.csv": sha256(COMPILED),
            "results/v3_human_homology_components.csv": sha256(COMPONENTS),
            "features/v2_features.npz": sha256(FEATURES / "v2_features.npz"),
            "protocols/v3_qtrp_analysis_plan.json": sha256(PLAN),
            "scripts/run_v3_qtrp_stack.py": sha256(SCRIPT),
        },
        "versions": {"python": sys.version, "numpy": np.__version__,
                     "lightgbm": v2_apply.members.lgb.__version__},
    })
    log(f"done in {round((time.time() - started) / 60, 2)} min")


if __name__ == "__main__":
    main()
