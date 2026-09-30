"""Does the QTRP-trained round survive depth control? Its stored scores answer it.

`reports/QTRP_DEPTH_CONFOUND.md` moved the primary caliber to raw negatives plus
detection-depth stratification, and `reports/V3_PLAN_QTRP_PERSULFIDATION.md`
names this as the first thing v3 must fix: the trained round
(`scripts/run_qtrp_trained_model.py`) reports plain AUCs on cohorts whose
negative definition - "free thiol in a pH 5.0 arm and never called SSH in ANY
human arm" - is the same mechanism that demoted `arm_clean`, and its own
`detection_depth1` probe already reached AUC 0.6468 [0.6045, 0.6866] with
`struct53 - detection_depth1` crossing zero at +0.0361 [-0.0173, +0.0904].

That round stored per-site out-of-fold scores in
`results/qtrp_trained_model_site_scores.csv`, so the depth-conditioned question
can be answered on the EXISTING fits. **No model is fitted in this script** and
the cohorts are not redefined; changing the negative definition to
`v3_pooled_ph5_raw` requires refitting and stays with the locked v3 round.

**Two depth columns, predeclared, because they are not the same thing.**

* `full_arm_depth`   - the number of distinct human QTRP arms in which the site
                       appears at all, in either adduct form, recomputed here
                       from `results/qtrp_sites_normalised.csv`.  This is the
                       definition used by the depth-conditioned round, values
                       1-7, and it is the PRIMARY stratifier.
* `cohort_arm_depth` - the `arms_detected` column already in the score table,
                       which counts only the arms the cohort admits (1-2 for the
                       pH 5.0 cohorts, 1-4 for `pooled_all_arms`).  Secondary.

Counting arms rather than rows is load-bearing: every persulfide site is also
detected as a free thiol, so a row count would encode the label - the error that
reported AUC 0.9478 and is documented in `reports/QTRP_DEPTH_CONFOUND.md`.

**Three estimands, predeclared.** The third exists because
`reports/DETECTION_MATCHED_DOUBLE_CONDITIONED.md` showed that a within-protein
and a global caliber can disagree in SIGN, so a stratified result reported in
only one caliber can be read backwards:

1. `plain`                   - global AUC, the published caliber.  Serves as the
                               consistency gate: it must reproduce
                               `results/qtrp_trained_model_summary.csv` to
                               within 0.001, since the scores are the same.
2. `depth_stratified`        - global concordance conditional on exact depth:
                               within-stratum AUC weighted by n_pos * n_neg.
3. `protein_x_depth`         - the same statistic with strata (protein, depth),
                               i.e. conditional on depth AND protein identity.
                               Every stratum is nested in a protein, so this is
                               a ratio of protein-level sums and the cluster
                               bootstrap resamples precomputed per-protein
                               (U, weight) pairs exactly.

Single-class strata carry no ranking information and are dropped; the positives
and negatives lost at each estimand are reported and must travel with the
numbers.

**Predeclared differences, on estimands 2 and 3:**

* chem1046 - struct53                     does the chemistry block add
* struct53_plus_chem1046 - chem1046       does structure add on top of it
* chem1046 - v2_chem_zeroshot             is training on this label worth it
* chem1046 - full_arm_depth1              is the fitted model more than depth
* struct53 - probe_visibility10           the headline contrast
* v2_chem_zeroshot - v2_full_zeroshot     the reversal that the depth round left
                                          unmeasured, flagged in
                                          reports/QTRP_INTERVALS_AND_CROSS_LABEL.md

**Predeclared reading rule**, unchanged: for an AUC a lower bound above 0.5 is
discrimination, an interval containing 0.5 is undetermined, an upper bound below
0.5 is inverted; for a difference the threshold is 0.

Cohorts: `pooled_ph5` is primary; `pooled_ph5_full_clean` and `pooled_all_arms`
are sensitivity checks, as in the round that produced the scores.

Intervals: protein cluster bootstrap, 5,000 replicates, seed 20260915.
"""
from __future__ import annotations

import collections
import json
import platform
import sys
import time

import numpy as np
import scipy
from scipy.stats import rankdata

from common import RESULTS, ROOT, read_csv, sha256, write_csv, write_json

SCRIPT = ROOT / "scripts" / "run_qtrp_trained_model_depth_conditioned.py"
SITE_SCORES = RESULTS / "qtrp_trained_model_site_scores.csv"
PUBLISHED = RESULTS / "qtrp_trained_model_summary.csv"
SITES = RESULTS / "qtrp_sites_normalised.csv"
DESIGNS = ("struct53", "chem1046", "struct53_plus_chem1046", "struct53_plus_zeroshot",
           "probe_digest25", "probe_visibility10", "detection_depth1",
           "v2_chem_zeroshot", "v2_full_zeroshot")
COHORTS = ("pooled_ph5", "pooled_ph5_full_clean", "pooled_all_arms")
GATE = 0.001
REPLICATES = 5000
SEED = 20260915
DIFFS = (
    ("chem1046", "struct53"),
    ("struct53_plus_chem1046", "chem1046"),
    ("chem1046", "v2_chem_zeroshot"),
    ("chem1046", "full_arm_depth1"),
    ("struct53", "probe_visibility10"),
    ("v2_chem_zeroshot", "v2_full_zeroshot"),
)


def fast_auc(y, scores):
    n1 = int(y.sum())
    n0 = int(len(y)) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    ranks = rankdata(scores, method="average")
    return float((ranks[y == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def stratified_auc(y, scores, strata):
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
    return total / total_weight if total_weight else float("nan")


def cell_u_and_weight(y, scores):
    """Mann-Whitney U with 0.5 for ties, and the weight n_positive * n_negative."""
    pos, neg = scores[y == 1], scores[y == 0]
    if pos.size == 0 or neg.size == 0:
        return 0.0, 0.0
    order = np.sort(neg)
    less = np.searchsorted(order, pos, side="left")
    equal = np.searchsorted(order, pos, side="right") - less
    return float(less.sum() + 0.5 * equal.sum()), float(pos.size * neg.size)


def full_arm_depth():
    """Distinct human QTRP arms in which each site appears, in either adduct form."""
    seen = collections.defaultdict(set)
    for row in read_csv(SITES):
        if row.get("species") != "human":
            continue
        label = row["class_from_adduct"] or row["class_from_sheet"]
        if label not in ("SH", "SSH"):
            continue
        seen[(row["accession"], int(row["site"]))].add((row["source_table"], row["ph"]))
    return {key: len(arms) for key, arms in seen.items()}


def protein_terms(y, proteins, cells, score_map):
    by_protein = collections.defaultdict(
        lambda: {"u": collections.defaultdict(float), "w": 0.0, "pos": 0, "neg": 0})
    order = np.argsort(cells, kind="stable")
    cell_sorted = cells[order]
    start = 0
    while start < len(order):
        stop = start
        while stop < len(order) and cell_sorted[stop] == cell_sorted[start]:
            stop += 1
        block = order[start:stop]
        yy = y[block]
        n1 = int(yy.sum())
        n0 = len(yy) - n1
        if n1 and n0:
            entry = by_protein[proteins[block[0]]]
            for name, vector in score_map.items():
                u, _ = cell_u_and_weight(yy, vector[block])
                entry["u"][name] += u
            entry["w"] += float(n1 * n0)
            entry["pos"] += n1
            entry["neg"] += n0
        start = stop
    return by_protein


def run_cohort(y, proteins, depth_full, depth_cohort, score_map):
    names = list(score_map)
    px_cells = np.asarray([f"{p}|d{d}" for p, d in zip(proteins, depth_full)])
    terms = protein_terms(y, proteins, px_cells, score_map)
    px_proteins = sorted(terms)
    px_weight = np.asarray([terms[p]["w"] for p in px_proteins])
    px_u = {n: np.asarray([terms[p]["u"][n] for p in px_proteins]) for n in names}

    unique = np.unique(proteins)
    index_of = {p: np.flatnonzero(proteins == p) for p in unique}
    position_of = {p: i for i, p in enumerate(px_proteins)}
    px_pick_index = [position_of.get(p, -1) for p in unique]

    point = {
        n: {
            "plain": fast_auc(y, score_map[n]),
            "depth_stratified": stratified_auc(y, score_map[n], depth_full),
            "cohort_depth_stratified": stratified_auc(y, score_map[n], depth_cohort),
            "protein_x_depth": float(px_u[n].sum() / px_weight.sum()),
        }
        for n in names
    }

    rng = np.random.default_rng(SEED)
    keys = ("plain", "depth_stratified", "protein_x_depth")
    draws = {f"{n}|{k}": [] for n in names for k in keys}
    diff_draws = {f"{a}__minus__{b}|{k}": []
                  for a, b in DIFFS for k in ("depth_stratified", "protein_x_depth")}
    skipped = 0
    for _ in range(REPLICATES):
        pick = rng.integers(0, len(unique), len(unique))
        index = np.concatenate([index_of[unique[i]] for i in pick])
        yy, dd = y[index], depth_full[index]
        if yy.sum() == 0 or yy.sum() == len(yy):
            skipped += 1
            continue
        px_rows = [px_pick_index[i] for i in pick if px_pick_index[i] >= 0]
        px_w = px_weight[px_rows].sum() if px_rows else 0.0
        values = {k: {} for k in keys}
        for n in names:
            vector = score_map[n][index]
            values["plain"][n] = fast_auc(yy, vector)
            values["depth_stratified"][n] = stratified_auc(yy, vector, dd)
            values["protein_x_depth"][n] = (
                float(px_u[n][px_rows].sum() / px_w) if px_w else float("nan"))
        for k in keys:
            for n in names:
                draws[f"{n}|{k}"].append(values[k][n])
        for k in ("depth_stratified", "protein_x_depth"):
            for a, b in DIFFS:
                left, right = values[k][a], values[k][b]
                if np.isfinite(left) and np.isfinite(right):
                    diff_draws[f"{a}__minus__{b}|{k}"].append(left - right)

    def interval(sample):
        array = np.asarray([v for v in sample if np.isfinite(v)], dtype=float)
        if array.size == 0:
            return [float("nan"), float("nan")]
        return [float(np.quantile(array, 0.025)), float(np.quantile(array, 0.975))]

    return {
        "usable_replicates": int(REPLICATES - skipped),
        "n_proteins_protein_x_depth": len(px_proteins),
        "positives_used_protein_x_depth": int(sum(terms[p]["pos"] for p in px_proteins)),
        "negatives_used_protein_x_depth": int(sum(terms[p]["neg"] for p in px_proteins)),
        "auc": {
            n: {**point[n], **{f"{k}_interval": interval(draws[f"{n}|{k}"]) for k in keys}}
            for n in names
        },
        "differences": {
            key: {"mean_difference": float(np.mean(sample)) if sample else float("nan"),
                  "percentile_95_interval": interval(sample)}
            for key, sample in diff_draws.items()
        },
    }


def verdict(bounds, centre=0.5):
    low, high = bounds
    if not np.isfinite(low) or not np.isfinite(high):
        return "undefined"
    if low > centre:
        return "discriminates"
    if high < centre:
        return "inverted"
    return "undetermined_at_this_sample_size"


def main():
    started = time.time()
    scores_rows = read_csv(SITE_SCORES)
    published = {(r["cohort"], r["design"]): float(r["roc_auc"]) for r in read_csv(PUBLISHED)}
    depth_map = full_arm_depth()

    summary, diff_rows, intervals, audit_cohorts = [], [], {}, {}
    for cohort in COHORTS:
        rows = [r for r in scores_rows if r["cohort"] == cohort]
        keys = [(r["accession"], int(r["site"])) for r in rows]
        missing = [k for k in keys if k not in depth_map]
        if missing:
            raise RuntimeError(f"{len(missing)} sites have no arm depth, e.g. {missing[:3]}")
        y = np.asarray([int(r["label"]) for r in rows])
        proteins = np.asarray([r["accession"] for r in rows])
        depth_full = np.asarray([depth_map[k] for k in keys])
        depth_cohort = np.asarray([int(r["arms_detected"]) for r in rows])
        score_map = {d: np.asarray([float(r[d]) for r in rows]) for d in DESIGNS}
        score_map["full_arm_depth1"] = depth_full.astype(float)

        gate = {d: round(fast_auc(y, score_map[d]) - published[(cohort, d)], 4) for d in DESIGNS}
        if max(abs(v) for v in gate.values()) > GATE:
            raise RuntimeError(f"Stored scores do not reproduce the published AUCs: {cohort} {gate}")

        strata_rows = []
        for value in np.unique(depth_full):
            mask = depth_full == value
            n1 = int(y[mask].sum())
            n0 = int(mask.sum()) - n1
            strata_rows.append({
                "cohort": cohort, "full_arm_depth": int(value), "n_sites": int(mask.sum()),
                "n_positive": n1, "n_negative": n0,
                "positive_rate": round(n1 / max(1, int(mask.sum())), 4),
                "informative": bool(n1 and n0),
            })
        dropped_pos = sum(r["n_positive"] for r in strata_rows if not r["informative"])
        dropped_neg = sum(r["n_negative"] for r in strata_rows if not r["informative"])

        boot = run_cohort(y, proteins, depth_full, depth_cohort, score_map)
        intervals[cohort] = boot
        audit_cohorts[cohort] = {
            "n_sites": len(rows), "n_positive": int(y.sum()),
            "n_proteins": int(len(set(proteins))),
            "full_arm_depth_values": [int(v) for v in np.unique(depth_full)],
            "cohort_arm_depth_values": [int(v) for v in np.unique(depth_cohort)],
            "mean_full_arm_depth_positive": round(float(depth_full[y == 1].mean()), 3),
            "mean_full_arm_depth_negative": round(float(depth_full[y == 0].mean()), 3),
            "positives_in_uninformative_depth_strata": dropped_pos,
            "negatives_in_uninformative_depth_strata": dropped_neg,
            "protein_x_depth_coverage": {
                "n_proteins": boot["n_proteins_protein_x_depth"],
                "positives_used": boot["positives_used_protein_x_depth"],
                "negatives_used": boot["negatives_used_protein_x_depth"],
            },
            "consistency_gate_vs_published": gate,
            "strata": strata_rows,
        }
        print(f"\n=== {cohort}  n={len(rows)} pos={int(y.sum())} "
              f"proteins={len(set(proteins))} depth_pos={audit_cohorts[cohort]['mean_full_arm_depth_positive']} "
              f"depth_neg={audit_cohorts[cohort]['mean_full_arm_depth_negative']} "
              f"dropped_pos={dropped_pos}", flush=True)
        for name in boot["auc"]:
            entry = boot["auc"][name]
            summary.append({
                "cohort": cohort, "design": name,
                "fit_regime": "reference_column" if name == "full_arm_depth1"
                              else ("zero_shot_annotation_trained" if name.endswith("_zeroshot")
                                    else "fitted_in_domain_protein_grouped_folds"),
                "n_sites": len(rows), "n_positive": int(y.sum()),
                "n_proteins": int(len(set(proteins))),
                "published_roc_auc": published.get((cohort, name), ""),
                "auc_plain": round(entry["plain"], 4),
                "auc_depth_stratified": round(entry["depth_stratified"], 4),
                "auc_depth_stratified_ci_low": round(entry["depth_stratified_interval"][0], 4),
                "auc_depth_stratified_ci_high": round(entry["depth_stratified_interval"][1], 4),
                "auc_cohort_depth_stratified": round(entry["cohort_depth_stratified"], 4),
                "auc_protein_x_depth": round(entry["protein_x_depth"], 4),
                "auc_protein_x_depth_ci_low": round(entry["protein_x_depth_interval"][0], 4),
                "auc_protein_x_depth_ci_high": round(entry["protein_x_depth_interval"][1], 4),
                "plain_minus_depth_stratified": round(entry["plain"] - entry["depth_stratified"], 4),
                "verdict_depth_stratified": verdict(entry["depth_stratified_interval"]),
                "verdict_protein_x_depth": verdict(entry["protein_x_depth_interval"]),
            })
            print(f"  {name:24s} plain {entry['plain']:.4f} depth {entry['depth_stratified']:.4f} "
                  f"[{entry['depth_stratified_interval'][0]:.4f}, {entry['depth_stratified_interval'][1]:.4f}] "
                  f"protein_x_depth {entry['protein_x_depth']:.4f} "
                  f"[{entry['protein_x_depth_interval'][0]:.4f}, {entry['protein_x_depth_interval'][1]:.4f}]",
                  flush=True)
        for a, b in DIFFS:
            for caliber in ("depth_stratified", "protein_x_depth"):
                entry = boot["differences"][f"{a}__minus__{b}|{caliber}"]
                low, high = entry["percentile_95_interval"]
                diff_rows.append({
                    "cohort": cohort, "caliber": caliber, "comparison": f"{a} - {b}",
                    "mean_difference": round(entry["mean_difference"], 4),
                    "ci_low": round(low, 4), "ci_high": round(high, 4),
                    "verdict": verdict(entry["percentile_95_interval"], centre=0.0),
                })
                print(f"    [{caliber:17s}] {a} - {b}: {entry['mean_difference']:+.4f} "
                      f"[{low:+.4f}, {high:+.4f}] "
                      f"{verdict(entry['percentile_95_interval'], centre=0.0)}", flush=True)

    write_csv(RESULTS / "qtrp_trained_depth_conditioned_summary.csv", summary)
    write_csv(RESULTS / "qtrp_trained_depth_conditioned_differences.csv", diff_rows)
    write_csv(RESULTS / "qtrp_trained_depth_strata.csv",
              [row for payload in audit_cohorts.values() for row in payload["strata"]])
    write_json(RESULTS / "qtrp_trained_depth_conditioned_intervals.json", intervals)
    write_json(RESULTS / "qtrp_trained_depth_conditioned_audit.json", {
        "completed": True,
        "elapsed_minutes": round((time.time() - started) / 60, 2),
        "host": platform.node(),
        "question": ("do the QTRP-trained round's results survive conditioning on detection "
                     "depth, and on depth together with protein identity"),
        "no_model_fitted": True,
        "score_source": "results/qtrp_trained_model_site_scores.csv, stored out-of-fold scores",
        "what_this_does_not_do": (
            "it does not change the negative definition; the cohorts here still exclude every "
            "site called SSH in any human arm, which is the arm_clean mechanism. The "
            "v3_pooled_ph5_raw cohort of reports/V3_PLAN_QTRP_PERSULFIDATION.md requires "
            "refitting and stays with the locked v3 round."
        ),
        "depth_definitions": {
            "full_arm_depth": ("distinct human QTRP arms in which the site appears at all, in "
                               "either adduct form, recomputed from results/qtrp_sites_normalised.csv; "
                               "PRIMARY stratifier"),
            "cohort_arm_depth": ("the arms_detected column in the score table, counting only arms the "
                                 "cohort admits; secondary"),
            "why_not_rows": ("every persulfide site is also detected as a free thiol, so a row count "
                             "would encode the label - the error that reported AUC 0.9478, documented "
                             "in reports/QTRP_DEPTH_CONFOUND.md"),
        },
        "estimands": {
            "plain": "global AUC, the published caliber, used as the consistency gate",
            "depth_stratified": "global concordance conditional on exact full_arm_depth",
            "protein_x_depth": ("concordance conditional on depth AND protein identity; strata are "
                                "nested in proteins so the cluster bootstrap resamples per-protein "
                                "(U, weight) pairs exactly"),
            "both_calibers_reported_because": ("reports/DETECTION_MATCHED_DOUBLE_CONDITIONED.md showed "
                                               "a within-protein and a global caliber can disagree in sign"),
        },
        "consistency_gate": {"threshold": GATE, "enforced": True},
        "bootstrap": {
            "unit": "protein", "replicates": REPLICATES, "seed": SEED,
            "reading_rule_predeclared_in_script": (
                "lower bound above 0.5 = discriminates; interval containing 0.5 = undetermined; "
                "upper bound below 0.5 = inverted; for differences the threshold is 0"),
        },
        "predeclared_differences": [f"{a} - {b}" for a, b in DIFFS],
        "cohorts": audit_cohorts,
        "limits": [
            "the negative definition is unchanged and still removes every site called SSH in any human arm, so these cohorts inherit the arm_clean covariate shift",
            "stratifying on depth cannot separate chemistry from occupancy: a site repeatedly seen as a persulfide may genuinely carry more persulfide",
            "single-class strata are dropped, which changes the effective cohort; the dropped counts are reported per cohort",
            "the protein-crossed estimand additionally drops proteins that carry only one class at a given depth; its coverage is reported",
            "all positives come from one paper and one cell line",
            "1,046 features against a few hundred positives risks overfitting; protein-grouped folds measure it but do not remove it",
            "NaHS-treated lysate at pH 5.0 is a strong exogenous sulfide load, not physiological persulfidation",
        ],
        "input_hashes": {
            "results/qtrp_trained_model_site_scores.csv": sha256(SITE_SCORES),
            "results/qtrp_trained_model_summary.csv": sha256(PUBLISHED),
            "results/qtrp_sites_normalised.csv": sha256(SITES),
            "scripts/run_qtrp_trained_model_depth_conditioned.py": sha256(SCRIPT),
        },
        "versions": {
            "python": sys.version, "numpy": np.__version__, "scipy": scipy.__version__,
            "note": ("run outside the project venv (3.9.6 / numpy 2.0.2); no model is fitted and only "
                     "rank statistics over stored out-of-fold scores are used"),
        },
    })


if __name__ == "__main__":
    main()
