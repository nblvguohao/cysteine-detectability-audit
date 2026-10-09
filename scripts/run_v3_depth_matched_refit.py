"""Is the chem1046-over-struct53 advantage separable from detection depth at all?

`reports/V3_QTRP_TRAINED_MODEL.md` left one question decidable and undecided.  On
the defensible raw caliber the advantage crosses zero (+0.0338 [-0.0205, +0.0864]);
it is clearly above zero only on the two cleaned cohorts, whose negative rule
inflates the depth confound (+0.1009 and +0.1077, with the score-depth Spearman
rising from 0.1476 to 0.2044 and 0.2696, and depth-alone AUC from 0.7340 to 0.8366
and 0.8391).  Stratification conditions on depth after the fit.  This round removes
depth from the cohort BEFORE the fit and refits everything from scratch.

**The design, and why it needs three cohorts rather than one.** A depth-matched
subcohort is necessarily smaller, so a difference that crosses zero there could mean
either "the advantage was depth-borne" or "the cohort got too small to see it".  The
two are separated by running a SIZE-MATCHED RANDOM cohort alongside: identical site
count and identical positive count, negatives drawn at random from the same raw pool
while ignoring depth.  Reading, fixed here before the run:

* advantage absent in the depth-matched cohort AND present in the size-matched
  random cohort  ->  the advantage is depth-borne;
* absent in BOTH  ->  the depth-matched null is a power statement, not evidence that
  depth carried the advantage;
* present in the depth-matched cohort  ->  the raw-caliber null in the v3 round was
  itself a power statement and the advantage survives depth removal.

Cohorts, all drawn from `v3_pooled_ph5_raw` (the caliber that retains negatives called
SSH in other arms), seed 20260915:

1. `depth_matched_1to2`  - primary. Within each exact-depth stratum: every positive,
   plus min(2 * n_positive, n_negative) negatives sampled without replacement. This is
   the subcohort the registered plan already named as the secondary analysis; the plan
   scored fixed out-of-fold scores on it, whereas this round REFITS on it.
2. `size_matched_random` - the power control described above.
3. `depth_matched_1to1`  - sensitivity; balances depth harder at the cost of half the
   negatives.

**Primary statistic is the plain AUC**, because depth is balanced by construction here
rather than conditioned on afterwards; the depth-stratified statistic is reported next
to it as a check, and the depth reference column is reported per cohort as the balance
audit - in a perfectly matched cohort it must sit at about 0.5 without any
stratification.

Everything else is inherited unchanged from `scripts/run_v3_qtrp_stack.py` so the
numbers stay comparable: homology-component outer and inner folds, the same
configuration and dimension grid, the NaN-safe mutual-information ranking with models
fitted on the raw matrix, cluster bootstrap over components with 5000 replicates and
seed 20260915, and the same three-way reading rule.

Predeclared differences: chem1046 - struct53 (the decisive one),
struct53_plus_chem1046 - chem1046, chem1046 - v2_chem_zeroshot, chem1046 - depth1,
struct53 - probe_visibility10.

No deployment criterion depends on this round; v3 is already decided as not deployed
(criterion 3 failed on the score-depth Spearman clause). This is a mechanism round.
"""
from __future__ import annotations

import collections
import json
import platform
import sys
import time

import numpy as np
from scipy.stats import spearmanr

from common import RESULTS, ROOT, read_csv, sha256, write_csv, write_json
import v2_apply
from run_qtrp_persulfidation_test import SITES, load_proteome, probe_matrix, structural_matrix
from run_qtrp_trained_model import human_matrix
from run_qtrp_depth_conditioned import fast_auc, stratified_auc, verdict
from run_v3_qtrp_stack import (
    COMPONENTS, PLAN, SEED, cluster_bootstrap, pooled_labels, select_and_score,
)
from run_cross_protease_detectability_probe import FEATURE_NAMES, VISIBILITY_ONLY
from v2_stack import FEATURES

SCRIPT = ROOT / "scripts" / "run_v3_depth_matched_refit.py"
COMPILED = RESULTS / "frbm_compiled_sites.csv"
PRIMARY = "v3_pooled_ph5_raw"
DIFFS = (
    ("chem1046", "struct53"),
    ("struct53_plus_chem1046", "chem1046"),
    ("chem1046", "v2_chem_zeroshot"),
    ("chem1046", "depth1"),
    ("struct53", "probe_visibility10"),
)


def depth_matched(y, strata, ratio, seed=SEED):
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


def size_matched_random(y, reference, seed=SEED):
    """Same site count and same positive count as `reference`, negatives ignoring depth."""
    rng = np.random.default_rng(seed + 1)
    n_pos = int(y[reference].sum())
    n_neg = int(reference.size - n_pos)
    positives = np.flatnonzero(y == 1)
    negatives = np.flatnonzero(y == 0)
    keep = np.concatenate([
        rng.choice(positives, size=min(n_pos, positives.size), replace=False),
        rng.choice(negatives, size=min(n_neg, negatives.size), replace=False)])
    return np.asarray(sorted(keep.tolist()))


def main():
    started = time.time()
    def log(message):
        print(message, flush=True)

    component_of = {r["accession"]: r["component"] for r in read_csv(COMPONENTS)}
    sequences = load_proteome()
    compiled_ssh = {(r["accession"], int(r["site"]))
                    for r in read_csv(COMPILED) if r["compiled_ssh"] == "1"}
    cohorts, depth = pooled_labels(read_csv(SITES), compiled_ssh)
    positives, negatives = cohorts[PRIMARY]
    visibility_columns = [i for i, n in enumerate(FEATURE_NAMES) if n in VISIBILITY_ONLY]

    union_keys = sorted(positives | negatives)
    residue_ok = np.asarray([
        bool(sequences.get(a)) and 1 <= p <= len(sequences[a]) and sequences[a][p - 1] == "C"
        for a, p in union_keys])
    S_all, _, structure_stats = structural_matrix(union_keys, sequences)
    P_all, probe_ok = probe_matrix(union_keys, sequences)
    covered = np.isfinite(S_all).any(axis=1) & residue_ok & probe_ok
    keys = [k for k, keep in zip(union_keys, covered) if keep and k[0] in component_of]
    index_map = {k: i for i, k in enumerate(union_keys)}
    rows_kept = np.asarray([index_map[k] for k in keys])
    S_all, P_all = S_all[rows_kept], P_all[rows_kept]
    V_all, columns, zero_shot, zero_shot_audit = human_matrix(keys, sequences)

    y_full = np.asarray([1 if k in positives else 0 for k in keys])
    comps_full = np.asarray([component_of[a] for a, _ in keys])
    strata_full = np.asarray([depth[k] for k in keys])
    designs_full = {
        "struct53": np.asarray(S_all, dtype=np.float64),
        "chem1046": np.asarray(V_all[:, columns["chem"]], dtype=np.float64),
        "struct53_plus_chem1046": np.asarray(
            np.hstack([S_all, V_all[:, columns["chem"]]]), dtype=np.float64),
        "struct53_plus_zeroshot": np.asarray(
            np.hstack([S_all, zero_shot["v2_chem_zeroshot"].reshape(-1, 1)]), dtype=np.float64),
        "probe_digest25": np.asarray(P_all, dtype=np.float64),
        "probe_visibility10": np.asarray(P_all[:, visibility_columns], dtype=np.float64),
    }
    log(json.dumps({"parent_cohort": PRIMARY, "sites": len(keys), "positives": int(y_full.sum()),
                    "components": len(set(comps_full)), **structure_stats,
                    "seconds": round(time.time() - started, 1)}))

    matched_1to2 = depth_matched(y_full, strata_full, 2)
    selections = {
        "depth_matched_1to2": matched_1to2,
        "size_matched_random": size_matched_random(y_full, matched_1to2),
        "depth_matched_1to1": depth_matched(y_full, strata_full, 1),
    }

    summary, intervals, audit_cohorts = [], {}, {}
    for cohort, index in selections.items():
        y, comps, strata = y_full[index], comps_full[index], strata_full[index]
        log(f"  {cohort}: {index.size} sites, {int(y.sum())} positive, {len(set(comps))} components")
        score_map, selection = {}, {}
        for name, matrix in designs_full.items():
            x = matrix[index]
            score_map[name], selection[name] = select_and_score(x, y, comps, name, log)
            if not np.isfinite(score_map[name]).all():
                raise RuntimeError(f"non-finite out-of-fold scores: {cohort} {name}")
        score_map["v2_chem_zeroshot"] = zero_shot["v2_chem_zeroshot"][index]
        score_map["v2_full_zeroshot"] = zero_shot["v2_full_zeroshot"][index]
        score_map["depth1"] = strata.astype(float)

        boot = cluster_bootstrap(y, comps, strata, score_map, DIFFS)
        intervals[cohort] = boot
        rates = {}
        for value in np.unique(strata):
            mask = strata == value
            n1 = int(y[mask].sum())
            rates[int(value)] = {"n": int(mask.sum()), "positive": n1,
                                 "negative": int(mask.sum()) - n1,
                                 "positive_rate": round(n1 / max(1, int(mask.sum())), 4)}
        for name, entry in boot["auc"].items():
            summary.append({
                "cohort": cohort, "design": name,
                "n_sites": int(index.size), "n_positive": int(y.sum()),
                "n_components": len(set(comps)),
                "auc_plain": round(entry["plain"], 4),
                "plain_ci_low": round(entry["plain_interval"][0], 4),
                "plain_ci_high": round(entry["plain_interval"][1], 4),
                "verdict_plain": verdict(entry["plain_interval"]),
                "auc_depth_stratified": round(entry["stratified"], 4),
                "stratified_ci_low": round(entry["stratified_interval"][0], 4),
                "stratified_ci_high": round(entry["stratified_interval"][1], 4),
                "spearman_score_vs_depth": round(float(spearmanr(score_map[name], strata).statistic), 4),
                "selection": json.dumps(selection.get(name, "not_fitted")),
            })
        audit_cohorts[cohort] = {
            "n_sites": int(index.size), "n_positive": int(y.sum()),
            "n_negative": int((y == 0).sum()), "n_components": len(set(comps)),
            "mean_depth_positive": round(float(strata[y == 1].mean()), 3),
            "mean_depth_negative": round(float(strata[y == 0].mean()), 3),
            "per_stratum": rates,
            "depth_alone_plain_auc": round(boot["auc"]["depth1"]["plain"], 4),
            "inner_selection": selection,
        }
        write_csv(RESULTS / "v3_depth_matched_refit_summary.csv", summary)
        write_json(RESULTS / "v3_depth_matched_refit_intervals.json", intervals)
        log(f"  {cohort} done at {round((time.time() - started) / 60, 1)} min")

    decisive = {}
    for cohort in selections:
        entry = intervals[cohort]["stratified_differences"]["chem1046__minus__struct53"]
        decisive[cohort] = {"mean_difference": round(entry["mean_difference"], 4),
                            "interval": [round(v, 4) for v in entry["percentile_95_interval"]],
                            "above_zero": bool(entry["percentile_95_interval"][0] > 0)}
    matched_above = decisive["depth_matched_1to2"]["above_zero"]
    random_above = decisive["size_matched_random"]["above_zero"]
    if not matched_above and random_above:
        reading = ("the advantage is depth-borne: it disappears once depth is balanced before the fit and "
                   "survives in a cohort of identical size whose negatives ignore depth")
    elif not matched_above and not random_above:
        reading = ("the depth-matched null is a power statement: the advantage is absent in a size-matched "
                   "cohort too, so removing depth is not what removed it")
    else:
        reading = ("the advantage survives depth removal, so the raw-caliber null in the v3 round was itself "
                   "a power statement")

    write_json(RESULTS / "v3_depth_matched_refit_audit.json", {
        "completed": True,
        "elapsed_minutes": round((time.time() - started) / 60, 2),
        "host": platform.node(),
        "question": ("is the chem1046-over-struct53 advantage separable from detection depth, or does it only "
                     "appear where the negative rule inflates that confound"),
        "relationship_to_the_registered_plan": (
            "the 1:2 depth-matched paired subcohort is the secondary analysis the registered plan already "
            "names; the plan scored fixed out-of-fold scores on it, this round REFITS on it. The "
            "size-matched random cohort is a declared addition, needed because a null in a smaller cohort is "
            "ambiguous between depth removal and power loss. No deployment criterion depends on this round."),
        "reading_rule_fixed_before_the_run": {
            "absent_in_matched_present_in_size_matched": "the advantage is depth-borne",
            "absent_in_both": "the matched null is a power statement, not evidence about depth",
            "present_in_matched": "the advantage survives depth removal"},
        "decisive_difference_chem1046_minus_struct53": decisive,
        "reading": reading,
        "cohorts": audit_cohorts,
        "parent_cohort": {"name": PRIMARY, "sites": len(keys), "positives": int(y_full.sum()),
                          "components": len(set(comps_full))},
        "zero_shot_members": zero_shot_audit,
        "limits": [
            "matching discards negatives, so every cohort here is smaller than the parent and all intervals are wider",
            "one matching draw per cohort, not a distribution over draws; the interval is the component bootstrap around that draw",
            "balancing depth cannot separate chemistry from occupancy either: a site repeatedly seen as a persulfide may really carry more persulfide",
            "the size-matched control balances size and positive count but not the depth distribution, which is the point",
        ],
        "input_hashes": {
            "results/qtrp_sites_normalised.csv": sha256(SITES),
            "results/frbm_compiled_sites.csv": sha256(COMPILED),
            "results/v3_human_homology_components.csv": sha256(COMPONENTS),
            "features/v2_features.npz": sha256(FEATURES / "v2_features.npz"),
            "protocols/v3_qtrp_analysis_plan.json": sha256(PLAN),
            "scripts/run_v3_depth_matched_refit.py": sha256(SCRIPT),
        },
        "versions": {"python": sys.version, "numpy": np.__version__,
                     "lightgbm": v2_apply.members.lgb.__version__},
    })
    log(json.dumps(decisive, indent=1))
    log(reading)
    log(f"done in {round((time.time() - started) / 60, 2)} min")


if __name__ == "__main__":
    main()
