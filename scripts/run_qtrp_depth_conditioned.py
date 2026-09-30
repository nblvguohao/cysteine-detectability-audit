"""The QTRP label has its own detectability confound. This conditions on it.

`reports/QTRP_PERSULFIDATION_TEST.md` treated the same-run free-thiol negative
set as free of the detectability problem, because both classes were seen by the
same probe in the same run.  That is true at the level of *peptide visibility*
and it is why the ten pure-visibility features fail there.  It is not true at
the level of *how often the site was seen*:

    depth = the number of distinct human QTRP arms in which the site appears at
            all, in either adduct form

Measured on the cohorts that report used, this single column beats every model
in it:

| arm | negatives | depth-alone AUC | struct53 |
|---|---|---|---|
| S1 pH5.0 | raw | 0.6833 | 0.6164 |
| S1 pH5.0 | arm_clean | 0.8176 | 0.6681 |
| S2 pH5.0 | raw | 0.7880 | 0.6892 |
| S2 pH5.0 | arm_clean | 0.9120 | 0.7240 |

and cleaning makes it worse rather than better, which is mechanical: `arm_clean`
removes the negatives that were called SSH in some other arm, and those are
exactly the negatives with high depth, so mean negative depth falls from 2.81 to
2.02 while positive depth stays at 3.88.  Part of the reported gain from
cleaning is therefore a covariate shift in the negative class, not noise
removal.

This script asks the question that survives: **at a fixed detection depth, does
anything still separate persulfidated from free-thiol-only cysteines?**

Method, predeclared here before the run:

* Strata are exact depth values.  A stratum with no positives or no negatives
  carries no information about ranking and is dropped; the dropped counts are
  reported.
* The conditional estimand is the depth-stratified concordance: the
  within-stratum AUC averaged with weights n_positive * n_negative, which is
  the probability that a randomly chosen positive outranks a randomly chosen
  negative *of the same depth*.
* Models are refitted per cohort with the same protein-grouped five folds and
  the same HistGradientBoosting configuration as every other probe row.  Depth
  is never a feature of any model except the `depth1` reference column.
* Intervals are the project's protein cluster bootstrap, 5,000 replicates, seed
  20260915, recomputing the stratified statistic on each draw.
* Reading rule, unchanged: lower bound above 0.5 is discrimination, an interval
  containing 0.5 is undetermined, an upper bound below 0.5 is inverted.

Predeclared comparisons, all on the stratified statistic:

* struct53 - probe_visibility10          the headline contrast, conditioned
* chem1046 - struct53                    does the v2 chemistry block add
* struct53 - v2_chem_zeroshot            in-domain structure against the
                                         deployed zero-shot ranker
* unconditional minus stratified, per design: how much of each reported AUC was
  detection depth rather than site chemistry.
"""
from __future__ import annotations

import collections
import json
import platform
import sys
import time

import numpy as np
from scipy.stats import rankdata
from sklearn.metrics import average_precision_score

from common import RESULTS, ROOT, read_csv, sha256, write_csv, write_json
import v2_apply
from run_qtrp_persulfidation_test import (
    SITES, grouped_scores, load_proteome, probe_matrix, structural_matrix,
)
from run_qtrp_negative_cleaning_test import COMPILED
from run_qtrp_trained_model import human_matrix
from run_cross_protease_detectability_probe import (
    FEATURE_NAMES, HGB_KWARGS, VISIBILITY_ONLY,
)
from v2_stack import FEATURES

SCRIPT = ROOT / "scripts" / "run_qtrp_depth_conditioned.py"
ARMS = (("S1", "5.0"), ("S2", "5.0"))
REPLICATES = 5000
SEED = 20260915
DIFFS = (
    ("struct53", "probe_visibility10"),
    ("chem1046", "struct53"),
    ("struct53", "v2_chem_zeroshot"),
)


def fast_auc(y, scores):
    n1 = int(y.sum())
    n0 = int(len(y)) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    ranks = rankdata(scores, method="average")
    return float((ranks[y == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def stratified_auc(y, scores, strata):
    """Weighted within-stratum concordance; weights are n_pos * n_neg."""
    total_weight, total = 0.0, 0.0
    for value in np.unique(strata):
        mask = strata == value
        n1 = int(y[mask].sum())
        n0 = int(mask.sum()) - n1
        if n1 == 0 or n0 == 0:
            continue
        weight = float(n1 * n0)
        total += weight * fast_auc(y[mask], scores[mask])
        total_weight += weight
    if total_weight == 0:
        return float("nan")
    return total / total_weight


def load_labels(rows, compiled_ssh):
    """Per-arm positives, the two negative definitions, and the depth column."""
    seen_arms = collections.defaultdict(set)
    per_arm = collections.defaultdict(lambda: collections.defaultdict(set))
    ssh_any = set()
    for row in rows:
        if row.get("species") != "human":
            continue
        label = row["class_from_adduct"] or row["class_from_sheet"]
        if label not in ("SH", "SSH"):
            continue
        key = (row["accession"], int(row["site"]))
        arm = (row["source_table"], row["ph"])
        seen_arms[key].add(arm)
        per_arm[arm][label].add(key)
        if label == "SSH":
            ssh_any.add(key)
    depth = {key: len(arms) for key, arms in seen_arms.items()}
    cohorts = {}
    for arm in ARMS:
        positives = set(per_arm[arm]["SSH"])
        raw = set(per_arm[arm]["SH"]) - positives
        cohorts[f"{arm[0]}_pH{arm[1]}|raw"] = (positives, raw)
        cohorts[f"{arm[0]}_pH{arm[1]}|arm_clean"] = (positives, raw - ssh_any)
        cohorts[f"{arm[0]}_pH{arm[1]}|full_clean"] = (
            positives, (raw - ssh_any) - compiled_ssh)
    return cohorts, depth


def cluster_bootstrap(y, groups, strata, score_map, diffs,
                      replicates=REPLICATES, seed=SEED):
    unique = np.unique(groups)
    index_of = {protein: np.flatnonzero(groups == protein) for protein in unique}
    names = list(score_map)
    draws = {f"{name}|{kind}": [] for name in names for kind in ("plain", "stratified")}
    diff_draws = {f"{a}__minus__{b}": [] for a, b in diffs}
    rng = np.random.default_rng(seed)
    skipped = 0
    for _ in range(replicates):
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
        "replicates": int(replicates),
        "usable_replicates": int(replicates - skipped),
        "seed": int(seed),
        "unit": "protein",
        "estimand": "stratified concordance conditional on exact detection depth, fixed out-of-fold scores",
        "auc": {
            name: {
                "plain": fast_auc(y, score_map[name]),
                "plain_interval": interval(draws[f"{name}|plain"]),
                "stratified": stratified_auc(y, score_map[name], strata),
                "stratified_interval": interval(draws[f"{name}|stratified"]),
            }
            for name in names
        },
        "stratified_differences": {
            key: {
                "mean_difference": float(np.mean(sample)) if sample else float("nan"),
                "percentile_95_interval": interval(sample),
            }
            for key, sample in diff_draws.items()
        },
    }


def verdict(bounds):
    low, high = bounds
    if not np.isfinite(low) or not np.isfinite(high):
        return "undefined"
    if low > 0.5:
        return "discriminates"
    if high < 0.5:
        return "inverted"
    return "undetermined_at_this_sample_size"


def main():
    started = time.time()
    rows = read_csv(SITES)
    sequences = load_proteome()
    compiled_ssh = {
        (r["accession"], int(r["site"]))
        for r in read_csv(COMPILED) if r["compiled_ssh"] == "1"
    }
    cohorts, depth = load_labels(rows, compiled_ssh)
    visibility_columns = [i for i, n in enumerate(FEATURE_NAMES) if n in VISIBILITY_ONLY]

    union_keys = sorted({key for pair in cohorts.values() for key in pair[0] | pair[1]})
    residue_ok = np.asarray([
        bool(sequences.get(a)) and 1 <= p <= len(sequences[a]) and sequences[a][p - 1] == "C"
        for a, p in union_keys
    ])
    S_all, _, structure_stats = structural_matrix(union_keys, sequences)
    P_all, probe_ok = probe_matrix(union_keys, sequences)
    covered = np.isfinite(S_all).any(axis=1) & residue_ok & probe_ok
    keys = [key for key, keep in zip(union_keys, covered) if keep]
    S_all, P_all = S_all[covered], P_all[covered]
    V_all, columns, zero_shot, zero_shot_audit = human_matrix(keys, sequences)
    row_of = {key: i for i, key in enumerate(keys)}
    print(json.dumps({"union": len(union_keys), "usable": len(keys),
                      "v2_matrix": list(V_all.shape), **structure_stats,
                      "seconds": round(time.time() - started, 1)}), flush=True)

    summary, intervals, audit_cohorts = [], {}, {}
    for cohort, (positives, negatives) in cohorts.items():
        selected = [key for key in keys if key in positives or key in negatives]
        y = np.asarray([1 if key in positives else 0 for key in selected])
        if y.sum() < 20 or (y == 0).sum() < 20:
            audit_cohorts[cohort] = {"skipped": "too few usable sites"}
            continue
        index = np.asarray([row_of[key] for key in selected])
        groups = np.asarray([a for a, _ in selected])
        strata = np.asarray([depth[key] for key in selected])

        designs = {
            "struct53": S_all[index],
            "chem1046": V_all[np.ix_(index, columns["chem"])],
            "probe_visibility10": P_all[np.ix_(index, visibility_columns)],
            "probe_digest25": P_all[index],
        }
        score_map = {}
        for name, x in designs.items():
            score_map[name] = grouped_scores(np.asarray(x, dtype=np.float64), y, groups)
            if not np.isfinite(score_map[name]).all():
                raise RuntimeError(f"Non-finite out-of-fold scores: {cohort} {name}")
        score_map["v2_chem_zeroshot"] = zero_shot["v2_chem_zeroshot"][index]
        score_map["v2_full_zeroshot"] = zero_shot["v2_full_zeroshot"][index]
        score_map["depth1"] = strata.astype(float)

        stratum_table = []
        dropped_positives = dropped_negatives = 0
        for value in np.unique(strata):
            mask = strata == value
            n1 = int(y[mask].sum())
            n0 = int(mask.sum()) - n1
            stratum_table.append({
                "cohort": cohort, "depth": int(value), "n_sites": int(mask.sum()),
                "n_positive": n1, "n_negative": n0,
                "positive_rate": round(n1 / max(1, int(mask.sum())), 4),
                "informative": bool(n1 and n0),
            })
            if not (n1 and n0):
                dropped_positives += n1
                dropped_negatives += n0

        boot = cluster_bootstrap(y, groups, strata, score_map, DIFFS)
        intervals[cohort] = boot
        audit_cohorts[cohort] = {
            "n_sites": len(selected),
            "n_positive": int(y.sum()),
            "n_proteins": int(len(set(groups))),
            "depth_values": [int(v) for v in np.unique(strata)],
            "informative_strata": sum(1 for row in stratum_table if row["informative"]),
            "positives_in_uninformative_strata": dropped_positives,
            "negatives_in_uninformative_strata": dropped_negatives,
            "mean_depth_positive": round(float(strata[y == 1].mean()), 3),
            "mean_depth_negative": round(float(strata[y == 0].mean()), 3),
            "strata": stratum_table,
        }
        for name in score_map:
            entry = boot["auc"][name]
            summary.append({
                "cohort": cohort,
                "design": name,
                "fit_regime": "zero_shot_annotation_trained" if name.endswith("_zeroshot")
                              else ("reference_column" if name == "depth1"
                                    else "fitted_in_domain_protein_grouped_folds"),
                "n_sites": len(selected),
                "n_positive": int(y.sum()),
                "n_proteins": int(len(set(groups))),
                "auc_plain": round(entry["plain"], 4),
                "auc_plain_ci_low": round(entry["plain_interval"][0], 4),
                "auc_plain_ci_high": round(entry["plain_interval"][1], 4),
                "auc_depth_stratified": round(entry["stratified"], 4),
                "auc_depth_stratified_ci_low": round(entry["stratified_interval"][0], 4),
                "auc_depth_stratified_ci_high": round(entry["stratified_interval"][1], 4),
                "plain_minus_stratified": round(entry["plain"] - entry["stratified"], 4),
                "verdict_stratified": verdict(entry["stratified_interval"]),
                "average_precision_plain": round(
                    float(average_precision_score(y, score_map[name])), 4),
            })
            print(f"{cohort:24s} {name:20s} plain {entry['plain']:.4f} "
                  f"stratified {entry['stratified']:.4f} "
                  f"[{entry['stratified_interval'][0]:.4f}, {entry['stratified_interval'][1]:.4f}] "
                  f"{verdict(entry['stratified_interval'])}", flush=True)
        for left, right in DIFFS:
            entry = boot["stratified_differences"][f"{left}__minus__{right}"]
            print(f"  {left} - {right}: {entry['mean_difference']:+.4f} "
                  f"[{entry['percentile_95_interval'][0]:+.4f}, "
                  f"{entry['percentile_95_interval'][1]:+.4f}]", flush=True)

    strata_rows = [row for payload in audit_cohorts.values()
                   for row in payload.get("strata", [])]
    write_csv(RESULTS / "qtrp_depth_conditioned_summary.csv", summary)
    write_csv(RESULTS / "qtrp_depth_strata.csv", strata_rows)
    write_json(RESULTS / "qtrp_depth_conditioned_intervals.json", intervals)
    write_json(RESULTS / "qtrp_depth_conditioned_audit.json", {
        "completed": True,
        "elapsed_minutes": round((time.time() - started) / 60, 2),
        "host": platform.node(),
        "question": (
            "at a fixed detection depth, does anything still separate "
            "persulfidated from free-thiol-only cysteines in the QTRP cohorts"
        ),
        "depth_definition": (
            "the number of distinct human QTRP arms in which the site appears at "
            "all, in either adduct form; counting rows instead would encode the "
            "label, because every persulfide site is also detected as a free thiol"
        ),
        "why_this_matters": (
            "the same-run negative set removes the peptide-visibility confound, "
            "which is why the ten pure-visibility features fail on this label, but "
            "it does not remove the repeated-detection confound; depth alone "
            "reaches AUC 0.6833 to 0.9120 across these cohorts, above every model "
            "reported on them, and negative-class cleaning raises it further "
            "because it removes high-depth negatives"
        ),
        "estimand": (
            "weighted within-stratum concordance, weights n_positive * n_negative; "
            "the probability that a positive outranks a negative of the same depth"
        ),
        "bootstrap": {
            "unit": "protein", "replicates": REPLICATES, "seed": SEED,
            "reading_rule_predeclared_in_script": (
                "lower bound above 0.5 = discriminates; interval containing 0.5 = "
                "undetermined; upper bound below 0.5 = inverted"
            ),
        },
        "predeclared_stratified_differences": [f"{a} - {b}" for a, b in DIFFS],
        "classifier": HGB_KWARGS,
        "zero_shot_models": zero_shot_audit,
        "cohorts": audit_cohorts,
        "structure_stats": structure_stats,
        "limits": [
            "stratifying on depth cannot separate chemistry from occupancy: a site repeatedly seen as a persulfide may genuinely carry more persulfide",
            "strata with one class only carry no ranking information and are dropped, which changes the effective cohort",
            "all positives come from one paper and one cell line",
            "negatives remain 'no persulfide seen in this paper' rather than chemically unmodified",
            "1,046 features against a few hundred positives risks overfitting; protein-grouped folds measure it but do not remove it",
        ],
        "input_hashes": {
            "results/qtrp_sites_normalised.csv": sha256(SITES),
            "results/frbm_compiled_sites.csv": sha256(COMPILED),
            "features/v2_features.npz": sha256(FEATURES / "v2_features.npz"),
            "scripts/run_qtrp_depth_conditioned.py": sha256(SCRIPT),
        },
        "versions": {
            "python": sys.version, "numpy": np.__version__,
            "lightgbm": v2_apply.members.lgb.__version__,
        },
    })
    print(json.dumps({k: v["stratified_differences"] for k, v in intervals.items()},
                     indent=2)[:1800])


if __name__ == "__main__":
    main()
