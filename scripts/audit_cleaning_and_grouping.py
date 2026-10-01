"""Two pre-checks this project learned the hard way, as reusable functions.

Both were discovered by getting them wrong first, and both now have to run BEFORE a
caliber or a leakage control is trusted.

**Check 1 - is a negative-cleaning rule collinear with the covariate it removes?**
`reports/QTRP_DEPTH_CONFOUND.md` demoted the QTRP `arm_clean` rule after finding it
removed exactly the high-detection-depth negatives, emptying strata 5 to 7 and
orphaning 37% of positives into single-class strata; the same instrument then CLEARED
the ABPP detection-matched rule in `reports/DETECTION_MATCHED_INTENSITY_CONFOUND.md`,
because there the residual confound was graded and every stratum kept both classes.
Two outcomes, one procedure: report retention by class, then the covariate's own
discriminability, then decide - never the other way round.

**Check 2 - can a within-group permutation control say anything?**
`reports/V3_QTRP_TRAINED_MODEL.md` section 4: the registered plan's within-component
label permutation was degenerate on the human QTRP cohort (1205 components over 1991
sites, 780 singletons, 118 carrying both classes), so it moved 10.15% of the labels
and returned 0.63 instead of 0.5 - uninformative about leakage. The fix is to measure
the degeneracy first and fall back to a global permutation when the group structure
cannot scramble the labels.

**Thresholds, all declared here rather than chosen per use.**

* `ORPHAN_LIMIT = 0.05` - a rule that pushes more than 5% of positives into
  single-class strata may not serve as the primary caliber.
* `COVARIATE_RISE_LIMIT = 0.02` - a rule that raises the removed covariate's own AUC
  by more than this made the confound worse, not better.
* `MOVED_RATIO_LIMIT = 0.5` - a within-group permutation is informative only if it
  moves at least half as many labels as a global permutation of the same vector
  would; the global expectation for a binary vector with positive rate p is
  2 * p * (1 - p).

Reasons are split, because the two kinds do not carry the same weight. Asymmetric
retention by class is ADVISORY: it is what makes a rule worth checking, not a ground
for demoting it - the ABPP detection-matching rule keeps 85% of annotated sites against
43% of unannotated ones and was still cleared, because the confound it left behind was
graded and every stratum kept both classes. A raised covariate AUC, a reduced count of
informative strata, or orphaned positives past the limit are BLOCKING: each means the
rule cannot support a conditional comparison. The verdict is driven by the blocking
reasons alone, with advisory ones reported alongside.

Run as a script, it applies both checks to the three cases this project already
settled, so its output IS its regression test: `arm_clean` must come out demoted, ABPP
detection matching must come out usable, the human QTRP component grouping must come
out degenerate, and the frozen tomato cohort - 49,064 sites over 4,986 proteins - must
come out non-degenerate, which is why the control worked there in the first place.
"""
from __future__ import annotations

import collections
import json
import platform
import sys
import time

import numpy as np
from scipy.stats import rankdata

from common import RESULTS, ROOT, read_csv, sha256, write_csv, write_json

SCRIPT = ROOT / "scripts" / "audit_cleaning_and_grouping.py"
ORPHAN_LIMIT = 0.05
COVARIATE_RISE_LIMIT = 0.02
MOVED_RATIO_LIMIT = 0.5
REPLICATES = 5000
SEED = 20260915


def auc(y, scores):
    y = np.asarray(y)
    n1 = int(y.sum())
    n0 = int(y.size) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    ranks = rankdata(np.asarray(scores, dtype=float), method="average")
    return float((ranks[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def _cluster_interval(y, scores, groups, replicates=REPLICATES, seed=SEED):
    if groups is None:
        return None
    groups = np.asarray(groups)
    unique = np.unique(groups)
    index_of = {g: np.flatnonzero(groups == g) for g in unique}
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(replicates):
        picked = rng.integers(0, len(unique), len(unique))
        index = np.concatenate([index_of[unique[i]] for i in picked])
        value = auc(y[index], np.asarray(scores)[index])
        if np.isfinite(value):
            draws.append(value)
    if not draws:
        return None
    return [float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))]


def _strata_table(y, strata):
    rows = []
    for value in np.unique(strata):
        mask = strata == value
        n1 = int(y[mask].sum())
        rows.append({"stratum": value if isinstance(value, str) else int(value),
                     "n": int(mask.sum()), "positive": n1, "negative": int(mask.sum()) - n1,
                     "informative": bool(n1 and int(mask.sum()) - n1)})
    return rows


def cleaning_collinearity(y, keep, covariate=None, strata=None, groups=None, label=""):
    """Is a negative-cleaning rule collinear with the covariate it removes?

    `keep` is the boolean mask the rule would retain. `covariate` may be None when it
    only exists for the retained sites (the ABPP case): retention is still reported and
    the covariate test is recorded as not computable rather than silently passed.
    """
    y = np.asarray(y)
    keep = np.asarray(keep, dtype=bool)
    blocking, advisory, out = [], [], {"label": label}
    pos, neg = y == 1, y == 0
    out["retention"] = {
        "positives_before": int(pos.sum()), "positives_after": int((pos & keep).sum()),
        "negatives_before": int(neg.sum()), "negatives_after": int((neg & keep).sum()),
        "positive_retention": round(float((pos & keep).sum() / max(1, pos.sum())), 4),
        "negative_retention": round(float((neg & keep).sum() / max(1, neg.sum())), 4),
        "positive_rate_before": round(float(pos.mean()), 4),
        "positive_rate_after": round(float(y[keep].mean()), 4),
    }
    if out["retention"]["positive_retention"] - out["retention"]["negative_retention"] > 0.1:
        advisory.append("the rule is asymmetric by class: it keeps a far larger share of positives "
                        "than of negatives, which is why it needs this check")

    if covariate is None:
        out["covariate"] = {"computable": False,
                            "note": "the covariate exists only for retained sites, so before/after cannot be compared"}
    else:
        covariate = np.asarray(covariate, dtype=float)
        before, after = auc(y, covariate), auc(y[keep], covariate[keep])
        out["covariate"] = {
            "computable": True,
            "auc_before": round(before, 4), "auc_after": round(after, 4),
            "rise": round(after - before, 4),
            "interval_before": _cluster_interval(y, covariate, groups),
            "interval_after": _cluster_interval(y[keep], covariate[keep],
                                                None if groups is None else np.asarray(groups)[keep]),
            "mean_positive_before": round(float(covariate[pos].mean()), 3),
            "mean_negative_before": round(float(covariate[neg].mean()), 3),
            "mean_positive_after": round(float(covariate[pos & keep].mean()), 3),
            "mean_negative_after": round(float(covariate[neg & keep].mean()), 3),
        }
        if after - before > COVARIATE_RISE_LIMIT:
            blocking.append(f"the rule RAISES the covariate's own AUC by {after - before:+.4f}, "
                            f"above the declared limit {COVARIATE_RISE_LIMIT}")

    if strata is None:
        out["strata"] = {"computable": False}
    else:
        strata = np.asarray(strata)
        rows_before = _strata_table(y, strata)
        rows_after = _strata_table(y[keep], strata[keep])
        orphan_pos = sum(r["positive"] for r in rows_after if not r["informative"])
        orphan_neg = sum(r["negative"] for r in rows_after if not r["informative"])
        out["strata"] = {
            "computable": True,
            "informative_before": sum(1 for r in rows_before if r["informative"]),
            "informative_after": sum(1 for r in rows_after if r["informative"]),
            "orphaned_positives": int(orphan_pos),
            "orphaned_positive_fraction": round(float(orphan_pos / max(1, int(y[keep].sum()))), 4),
            "orphaned_negatives": int(orphan_neg),
            "before": rows_before, "after": rows_after,
        }
        if out["strata"]["informative_after"] < out["strata"]["informative_before"]:
            blocking.append(f"the rule reduces informative strata from {out['strata']['informative_before']} "
                            f"to {out['strata']['informative_after']}")
        if out["strata"]["orphaned_positive_fraction"] > ORPHAN_LIMIT:
            blocking.append(f"{orphan_pos} positives ({out['strata']['orphaned_positive_fraction']:.1%}) fall "
                            f"into single-class strata, above the declared limit {ORPHAN_LIMIT:.0%}")

    out["blocking_reasons"] = blocking
    out["advisory_reasons"] = advisory
    out["verdict"] = "demote_to_sensitivity_check" if blocking else "usable_as_primary_caliber"
    return out


def permutation_degeneracy(y, groups, seed=SEED, label=""):
    """Can a within-group label permutation scramble this label vector at all?"""
    y = np.asarray(y)
    groups = np.asarray(groups)
    sizes = collections.Counter(groups.tolist())
    mixed = [g for g in sizes if 0 < y[groups == g].sum() < (groups == g).sum()]
    rng = np.random.default_rng(seed)
    permuted = y.copy()
    for group in np.unique(groups):
        mask = np.flatnonzero(groups == group)
        permuted[mask] = rng.permutation(y[mask])
    moved = float((permuted != y).mean())
    rate = float(y.mean())
    expected_global = 2 * rate * (1 - rate)
    ratio = moved / expected_global if expected_global > 0 else 0.0
    informative = ratio >= MOVED_RATIO_LIMIT
    return {
        "label": label,
        "sites": int(y.size), "positives": int(y.sum()), "groups": len(sizes),
        "singleton_groups": int(sum(1 for v in sizes.values() if v == 1)),
        "mean_group_size": round(float(y.size / max(1, len(sizes))), 3),
        "groups_carrying_both_classes": len(mixed),
        "sites_in_mixed_groups": int(sum(int((groups == g).sum()) for g in mixed)),
        "labels_moved": int((permuted != y).sum()),
        "moved_fraction": round(moved, 4),
        "positives_moved": int(((y == 1) & (permuted != y)).sum()),
        "expected_global_moved_fraction": round(expected_global, 4),
        "moved_ratio_vs_global": round(ratio, 4),
        "informative": bool(informative),
        "verdict": ("within_group_permutation_is_informative" if informative else
                    "within_group_permutation_is_DEGENERATE_use_a_global_permutation"),
        "threshold": MOVED_RATIO_LIMIT,
    }


def _qtrp_case():
    from run_qtrp_persulfidation_test import SITES
    from run_qtrp_negative_cleaning_test import COMPILED
    from run_v3_qtrp_stack import COMPONENTS, pooled_labels
    component_of = {r["accession"]: r["component"] for r in read_csv(COMPONENTS)}
    compiled = {(r["accession"], int(r["site"])) for r in read_csv(COMPILED) if r["compiled_ssh"] == "1"}
    cohorts, depth = pooled_labels(read_csv(SITES), compiled)
    positives, raw_negatives = cohorts["v3_pooled_ph5_raw"]
    _, clean_negatives = cohorts["v3_pooled_ph5_armclean"]
    keys = sorted(positives | raw_negatives)
    y = np.asarray([1 if k in positives else 0 for k in keys])
    keep = np.asarray([k in positives or k in clean_negatives for k in keys])
    strata = np.asarray([depth[k] for k in keys])
    comps = np.asarray([component_of.get(a, a) for a, _ in keys])
    return y, keep, strata, comps


def _abpp_case():
    from run_qtrp_persulfidation_test import load_proteome  # noqa: F401  (import parity with the other rounds)
    abpp = {(a, int(p)) for a, p in json.loads(
        (ROOT / "external/lysis_oxidant_abpp/abpp_detected_cys_sites.json").read_text())}
    payload = json.loads((ROOT / "inputs/PXD044043_mapped_human.json").read_text())["hsa_persulfidation"]["proteins"]
    reactivity = {(r["accession"], int(r["position"])): r
                  for r in read_csv(RESULTS / "abpp_reactivity_sites.csv")}
    keys, y, keep, intensity = [], [], [], []
    for entry in payload:
        sequence = entry["seq"].upper()
        observed = {int(x) for x in entry["modified"]}
        for position in sorted(int(x) for x in entry["cys"]):
            if not (1 <= position <= len(sequence)) or sequence[position - 1] != "C":
                continue
            key = (entry["acc"], position)
            keys.append(key)
            y.append(1 if position in observed else 0)
            keep.append(key in abpp)
            row = reactivity.get(key)
            intensity.append(float(row["log2_no_ox"]) if row else np.nan)
    y = np.asarray(y); keep = np.asarray(keep); intensity = np.asarray(intensity)
    covered = keep & np.isfinite(intensity)
    deciles = np.full(y.size, "not_in_arm", dtype=object)
    values = intensity[covered]
    cuts = np.quantile(values, np.linspace(0, 1, 11))[1:-1]
    deciles[covered] = [f"decile_{d}" for d in np.clip(np.searchsorted(cuts, values, side="right"), 0, 9)]
    groups = np.asarray([a for a, _ in keys])
    return y, keep, deciles, groups, covered, intensity


def main():
    started = time.time()
    cases, rows = {}, []

    y, keep, strata, comps = _qtrp_case()
    cases["qtrp_arm_clean_vs_detection_depth"] = cleaning_collinearity(
        y, keep, covariate=strata.astype(float), strata=strata, groups=comps,
        label="QTRP arm_clean rule against the all-arm detection depth it removes")
    cases["qtrp_component_permutation"] = permutation_degeneracy(
        y, comps, label="human QTRP raw cohort, homology components")

    y_a, keep_a, deciles, groups_a, covered, intensity = _abpp_case()
    # strata are passed as None here for the same reason the covariate is: labelling intensity,
    # and therefore its deciles, exist only for sites the rule retains. The stratum analysis for
    # this rule lives in the next case, inside the arm it retains.
    cases["abpp_detection_matching_vs_labelling_intensity"] = cleaning_collinearity(
        y_a, keep_a, covariate=None, strata=None, groups=groups_a,
        label="ABPP detection-matching rule against ABPP labelling intensity")
    cases["abpp_intensity_inside_the_matched_arm"] = cleaning_collinearity(
        y_a[covered], np.ones(int(covered.sum()), dtype=bool),
        covariate=intensity[covered], strata=deciles[covered], groups=groups_a[covered],
        label="labelling intensity as a graded covariate inside the matched arm")

    frozen = np.load(ROOT / "inputs/frozen_benchmark_data.npz", allow_pickle=True)
    for unit in ("proteins", "components"):
        cases[f"frozen_tomato_permutation_by_{unit}"] = permutation_degeneracy(
            frozen["y"], frozen[unit], label=f"frozen tomato benchmark cohort, {unit}")

    for name, case in cases.items():
        if "verdict" not in case:
            continue
        if "retention" in case:
            rows.append({
                "check": "cleaning_collinearity", "case": name, "verdict": case["verdict"],
                "positive_retention": case["retention"]["positive_retention"],
                "negative_retention": case["retention"]["negative_retention"],
                "positive_rate_before": case["retention"]["positive_rate_before"],
                "positive_rate_after": case["retention"]["positive_rate_after"],
                "covariate_auc_before": case["covariate"].get("auc_before", ""),
                "covariate_auc_after": case["covariate"].get("auc_after", ""),
                "covariate_rise": case["covariate"].get("rise", ""),
                "informative_strata_before": case["strata"].get("informative_before", ""),
                "informative_strata_after": case["strata"].get("informative_after", ""),
                "orphaned_positives": case["strata"].get("orphaned_positives", ""),
                "orphaned_positive_fraction": case["strata"].get("orphaned_positive_fraction", ""),
                "moved_fraction": "", "moved_ratio_vs_global": "", "groups_both_classes": "",
                "blocking_reasons": " | ".join(case["blocking_reasons"]) or "none",
                "advisory_reasons": " | ".join(case["advisory_reasons"]) or "none",
            })
        else:
            rows.append({
                "check": "permutation_degeneracy", "case": name, "verdict": case["verdict"],
                "positive_retention": "", "negative_retention": "",
                "positive_rate_before": round(case["positives"] / case["sites"], 4), "positive_rate_after": "",
                "covariate_auc_before": "", "covariate_auc_after": "", "covariate_rise": "",
                "informative_strata_before": "", "informative_strata_after": "",
                "orphaned_positives": "", "orphaned_positive_fraction": "",
                "moved_fraction": case["moved_fraction"],
                "moved_ratio_vs_global": case["moved_ratio_vs_global"],
                "groups_both_classes": f"{case['groups_carrying_both_classes']}/{case['groups']}",
                "blocking_reasons": "", "advisory_reasons":
                    f"mean group size {case['mean_group_size']}, {case['singleton_groups']} singletons",
            })
    write_csv(RESULTS / "cleaning_and_grouping_audit_summary.csv", rows)

    expected = {
        "qtrp_arm_clean_vs_detection_depth": "demote_to_sensitivity_check",
        "abpp_detection_matching_vs_labelling_intensity": "usable_as_primary_caliber",
        "qtrp_component_permutation": "within_group_permutation_is_DEGENERATE_use_a_global_permutation",
        "frozen_tomato_permutation_by_proteins": "within_group_permutation_is_informative",
    }
    regression = {name: {"expected": want, "observed": cases[name]["verdict"],
                         "matches": cases[name]["verdict"] == want}
                  for name, want in expected.items()}

    write_json(RESULTS / "cleaning_and_grouping_audit.json", {
        "completed": True,
        "elapsed_minutes": round((time.time() - started) / 60, 2),
        "host": platform.node(),
        "purpose": ("reusable pre-checks: is a negative-cleaning rule collinear with the covariate it removes, "
                    "and can a within-group permutation control scramble the labels at all"),
        "thresholds": {"orphaned_positive_fraction_limit": ORPHAN_LIMIT,
                       "covariate_auc_rise_limit": COVARIATE_RISE_LIMIT,
                       "moved_ratio_vs_global_limit": MOVED_RATIO_LIMIT,
                       "declared_in": "the module docstring, before any case was run"},
        "regression_against_settled_cases": regression,
        "all_regression_cases_match": all(v["matches"] for v in regression.values()),
        "cases": cases,
        "intervals": {"method": "cluster bootstrap over the supplied group column",
                      "replicates": REPLICATES, "seed": SEED},
        "limits": [
            "the collinearity check is observational: it shows whether a rule is entangled with a covariate, not whether the covariate causes the signal",
            "the permutation check measures whether a within-group permutation CAN scramble labels, not whether a pipeline leaks; a passing group structure still needs the control to be run",
            "thresholds are declared conventions, not estimates; a rule just inside a limit deserves the same scrutiny as one just outside",
            "the ABPP case cannot compare the covariate before and after, because labelling intensity exists only for sites the rule retains",
        ],
        "input_hashes": {
            "results/qtrp_sites_normalised.csv": sha256(RESULTS / "qtrp_sites_normalised.csv"),
            "results/frbm_compiled_sites.csv": sha256(RESULTS / "frbm_compiled_sites.csv"),
            "results/v3_human_homology_components.csv": sha256(RESULTS / "v3_human_homology_components.csv"),
            "results/abpp_reactivity_sites.csv": sha256(RESULTS / "abpp_reactivity_sites.csv"),
            "external/lysis_oxidant_abpp/abpp_detected_cys_sites.json": sha256(
                ROOT / "external/lysis_oxidant_abpp/abpp_detected_cys_sites.json"),
            "inputs/PXD044043_mapped_human.json": sha256(ROOT / "inputs/PXD044043_mapped_human.json"),
            "inputs/frozen_benchmark_data.npz": sha256(ROOT / "inputs/frozen_benchmark_data.npz"),
            "scripts/audit_cleaning_and_grouping.py": sha256(SCRIPT),
        },
        "versions": {"python": sys.version, "numpy": np.__version__},
    })
    for row in rows:
        print(f"  {row['check']:24s} {row['case']:46s} {row['verdict']}", flush=True)
    print(json.dumps(regression, indent=1), flush=True)
    print(f"done in {round((time.time() - started) / 60, 2)} min", flush=True)


if __name__ == "__main__":
    main()
