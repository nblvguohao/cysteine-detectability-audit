"""Diagnostics a site-level PTM dataset should pass before its site preferences are believed.

This project reached three conclusions the hard way, each after getting it wrong once, and
each is now a function here rather than a paragraph in a report:

  * `reports/PTM_DETECTABILITY_CENSUS.md` / `reports/PHASE1_DETECTABILITY_SYNTHESIS.md`:
    whether a site table is dominated by mass-spectrometric detectability is a property of
    how its NEGATIVES were constructed, not of the modification chemistry.  Across 9 usable
    datasets, 3 switch chemistries and 3 species, pure-visibility features take 0.5963-0.9488
    of the achievable discriminability when negatives are "not pulled down", and that share
    falls by a mean 0.4645 when negatives become "seen in the same run but unmodified".
  * `reports/QTRP_DEPTH_CONFOUND.md`: a cleaning rule that removes the negatives carrying a
    confounder does not remove the confounder - it empties strata and orphans positives.
  * `reports/PHASE2_CLAIMS_UNDER_DETECTABILITY_CONTROL.md` and
    `reports/PHASE2B_CLAIMS_BACKFILL.md`: which published claim collapses under a
    detectability control is NOT predictable from a summary statistic.  The main session
    hypothesised that effect loss scales with the detectability propensity AUC and tested it:
    n=9, Spearman rho -0.2510, exact permutation p 0.5123, counterexample SFE-008 (propensity
    AUC 0.6220, retention 0.086).  The hypothesis is refuted, and that is precisely why a
    per-claim instrument is needed instead of a rule of thumb.

Everything declared below is fixed BEFORE any case runs, and running this file runs its own
regression: it must reproduce decisions that were already settled elsewhere in the tree.

DECLARED DEFINITIONS AND THRESHOLDS
-----------------------------------
share = (AUC_visibility - 0.5) / (AUC_digestion - 0.5), the fraction of the achievable
    discriminability that pure-visibility features already reach.  UNDEFINED when the
    denominator is <= SHARE_DENOMINATOR_FLOOR (0.02): if the richer feature set barely beats
    chance, the ratio is numerically unstable and must not be reported.  Reported without
    clipping; a share above 1 means the visibility subset beats its own superset and is a
    finding about the cohort, not an error to hide.

label semantics tiers - what a dataset's negative class licenses:
    T1_direct_adduct      both forms of the same cysteine observed in the same run, channel
                          assignment mutually exclusive at SITE level -> chemistry claims OK
    T2_parallel_ambiguous same-run parallel labelling but assignment not mutually exclusive
                          at site level, or negatives are same-run "detected but unmodified"
                          drawn from an enrichment run -> chemistry claims only WITH a
                          detectability control reported alongside
    T3_capture_annotation negatives are "not observed" (whole-proteome background, all Cys in
                          identified proteins, unmodified Cys in the same protein, curated
                          negatives from another database, random or simulated) -> only
                          detectability-ranking claims; chemistry claims NOT supported
    The tier is a property of the NEGATIVE SET, not of the chemistry or the species.  A T2
    dataset is not "better chemistry" than T3; it is a different question being asked.

depth confound strength: depth-alone AUC computed tie-corrected from the stratum table,
    plus informative strata (both classes present) and orphaned positives (positives in
    single-class strata).  ORPHANED_LIMIT 0.05 - above it the stratified statistic is not
    trustworthy, which is what demoted the QTRP arm_clean caliber.

claim verdict: the rule predeclared in results/phase2_claim_retest_audit.json, with the
    direction guard added in results/phase2b_claim_retest_audit.json.  Implemented here as a
    pure function so it can be regression-tested and reused; ATTENUATION_FLOOR 0.5.

REGRESSION - all of these must reproduce, or this file fails
------------------------------------------------------------
R1 every usable share in results/ptm_detectability_share.csv, recomputed from the stored
   per-feature-set global AUCs.
R2 depth-alone AUC for QTRP S1_pH5.0|raw and S2_pH5.0|raw, recomputed from
   results/qtrp_depth_strata.csv, against design=depth1 in
   results/qtrp_depth_conditioned_summary.csv. This case earned its keep on the first run:
   the concordance direction was inverted and it returned 1 - AUC (0.3172 for 0.6828).
   Plus arm_clean on both arms must be reported as NOT trustworthy by the orphaning rule.
R3 every verdict in results/phase2b_claim_retest_combined.csv that has the required numbers,
   replayed through claim_verdict().
R4 two synthetic cleaning rules with analytically known answers: a rule that removes exactly
   the high-covariate negatives must be BLOCKING; a rule that removes negatives uniformly at
   random must not be.
R6 every one of claim_verdict's 13 return paths, each landing on a distinct reason string.
   Added after the independent verification round found the 22 stored claims exercise only 7,
   so six paths were reachable but untested.
R5 the settled precedents recorded in results/cleaning_and_grouping_audit.json (QTRP
   arm_clean demoted, ABPP detection-matching cleared) - read from that audit rather than
   recomputed, and labelled as such.

KNOWN LIMITS OF THESE CHECKS
----------------------------
* The share needs a visibility feature block and a superset; on a dataset with no digestion
  model it is not computable and returns undefined rather than a number.
* The tier classifier reads declared metadata. It cannot detect that a paper's stated
  mutually-exclusive channels are not mutually exclusive in its own supplementary table -
  that is what caught PXD005168 (112 of 187 direct-channel sites also reported as blocked
  free thiols) and it took reading the table.
* Depth confound strength assumes the depth covariate does not encode the label. A covariate
  defined as "number of arms in which the site CLEARS the criterion" does encode it; this
  project made that mistake twice (row-counted probe depth, tier B arm depth).
* claim_verdict is a verdict rule, not a pipeline: it takes effect estimates as given.

Terminology: sulfenylation = 次磺酰化 (-SOH); sulfinylation = 亚磺酰化 (-SO2H).
"""
import csv, hashlib, json, sys, time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
SCRIPT = Path(__file__).resolve()

SHARE_DENOMINATOR_FLOOR = 0.02
ORPHANED_LIMIT = 0.05
ATTENUATION_FLOOR = 0.5
SEED = 20260915

TIER_BY_NEGATIVE_CLASS = {
    "same_run_mutually_exclusive": "T1_direct_adduct",
    "detected_background_proteome": "T2_parallel_ambiguous",
    "whole_proteome_annotation_background": "T3_capture_annotation",
    "all_cys_in_identified_proteins": "T3_capture_annotation",
    "unmodified_cys_same_protein": "T3_capture_annotation",
    "curated_negatives_other_db": "T3_capture_annotation",
    "random_or_simulated": "T3_capture_annotation",
    "not_defined": "unclassifiable",
}
TIER_LICENCE = {
    "T1_direct_adduct": "chemistry claims permitted",
    "T2_parallel_ambiguous": "chemistry claims permitted only with a detectability control reported",
    "T3_capture_annotation": "detectability-ranking claims only; chemistry claims not supported",
    "unclassifiable": "no claim; the negative set must be defined first",
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def auc(y, scores):
    y = np.asarray(y, dtype=float)
    scores = np.asarray(scores, dtype=float)
    pos, neg = scores[y == 1], scores[y == 0]
    if pos.size == 0 or neg.size == 0:
        return float("nan")
    order = np.argsort(np.concatenate([pos, neg]), kind="mergesort")
    ranks = np.empty(order.size, dtype=float)
    ranks[order] = np.arange(1, order.size + 1, dtype=float)
    values = np.concatenate([pos, neg])[order]
    i = 0
    while i < values.size:
        j = i
        while j + 1 < values.size and values[j + 1] == values[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = ranks[order[i:j + 1]].mean()
        i = j + 1
    return (ranks[: pos.size].sum() - pos.size * (pos.size + 1) / 2) / (pos.size * neg.size)


# ---------------------------------------------------------------- check 1: share
def detectability_share(auc_visibility, auc_digestion):
    """Fraction of achievable discriminability already reached by pure-visibility features."""
    denominator = auc_digestion - 0.5
    if not np.isfinite(denominator) or denominator <= SHARE_DENOMINATOR_FLOOR:
        return {"share": None, "defined": False,
                "reason": f"denominator {denominator:.4f} <= floor {SHARE_DENOMINATOR_FLOOR}; "
                          "the superset barely beats chance so the ratio is unstable"}
    return {"share": (auc_visibility - 0.5) / denominator, "defined": True, "reason": "ok"}


# --------------------------------------------------- check 2: label semantics tier
def label_semantics_tier(negative_class=None, same_run=None, mutually_exclusive_at_site=None):
    """Which claims a dataset's negative set licenses. Declared metadata in, tier out."""
    path = []
    if negative_class is not None:
        tier = TIER_BY_NEGATIVE_CLASS.get(negative_class, "unclassifiable")
        path.append(f"negative_class={negative_class} -> {tier}")
    elif same_run is None:
        tier, = ("unclassifiable",)
        path.append("no negative_class and no same_run flag -> unclassifiable")
    elif not same_run:
        tier = "T3_capture_annotation"
        path.append("negatives not from the same run -> T3")
    elif mutually_exclusive_at_site:
        tier = "T1_direct_adduct"
        path.append("same run and site-level mutually exclusive channels -> T1")
    else:
        tier = "T2_parallel_ambiguous"
        path.append("same run but assignment not mutually exclusive at site level -> T2")
    return {"tier": tier, "licence": TIER_LICENCE[tier], "decision_path": path}


# ------------------------------------------- check 3: depth confound strength
def depth_confound_strength(strata):
    """strata: iterable of (depth_value, n_positive, n_negative). Tie-corrected depth-alone AUC."""
    rows = sorted(((float(d), int(p), int(n)) for d, p, n in strata), key=lambda r: r[0])
    total_pos = sum(r[1] for r in rows)
    total_neg = sum(r[2] for r in rows)
    # AUC = P(depth of a positive > depth of a negative) + 0.5 P(tie): for each positive
    # stratum, count the negatives in the strata BELOW it. Getting this direction wrong
    # returns 1 - AUC, which is what the R2 regression caught on its first run.
    concordant = ties = 0.0
    for i, (_, p_i, _) in enumerate(rows):
        concordant += p_i * sum(r[2] for r in rows[:i])
        ties += p_i * rows[i][2]
    depth_auc = (concordant + 0.5 * ties) / (total_pos * total_neg) if total_pos and total_neg else float("nan")
    informative = [r for r in rows if r[1] > 0 and r[2] > 0]
    orphaned = sum(r[1] for r in rows if r[1] > 0 and r[2] == 0)
    fraction = orphaned / total_pos if total_pos else float("nan")
    # AUC is defined with depth increasing; a cohort where positives are shallower gives < 0.5
    return {"depth_alone_auc": depth_auc, "strata": len(rows), "informative_strata": len(informative),
            "orphaned_positives": orphaned, "orphaned_fraction": fraction,
            "trustworthy": bool(fraction <= ORPHANED_LIMIT),
            "reason": ("ok" if fraction <= ORPHANED_LIMIT else
                       f"orphaned positive fraction {fraction:.4f} > limit {ORPHANED_LIMIT}: "
                       "the stratified statistic is not trustworthy on this caliber")}


# ------------------------------------------------- check 4: claim verdict rule
def claim_verdict(baseline, baseline_ci, controlled_ci, random_control_ci,
                  controlled_point, claim_direction="positive_preference",
                  baseline_reproduced=True, precheck_blocking=False):
    """The rule predeclared in phase 2, with the phase-2b direction guard. Pure function."""
    def excludes_zero(ci):
        return ci is not None and (ci[0] > 0 or ci[1] < 0)
    if precheck_blocking:
        return {"verdict": "undecidable", "reason": "a precheck is BLOCKING on this caliber"}
    if baseline is None or baseline_ci is None:
        return {"verdict": "out_of_instrument_scope", "reason": "no baseline estimate"}
    if claim_direction == "null_no_preference":
        if excludes_zero(baseline_ci):
            return {"verdict": "undecidable", "reason": "the author's null does not reproduce here"}
        if excludes_zero(controlled_ci):
            return {"verdict": "reverses", "reason": "null broken by the control"}
        return {"verdict": "survives", "reason": "null survives; quote with the null resolution"}
    if excludes_zero(baseline_ci) and claim_direction == "positive_preference" and baseline < 0:
        return {"verdict": "baseline_contradicts_claim",
                "reason": "baseline effect is opposite to the claim; a surviving effect is not the claim"}
    if not baseline_reproduced and not excludes_zero(baseline_ci):
        return {"verdict": "undecidable", "reason": "baseline not reproducible and not established here"}
    if not excludes_zero(baseline_ci):
        return {"verdict": "undecidable", "reason": "baseline interval crosses zero in our hands"}
    if excludes_zero(controlled_ci):
        if controlled_point is not None and np.sign(controlled_point) != np.sign(baseline):
            return {"verdict": "reverses", "reason": "sign flipped with the interval clear of zero"}
        if controlled_point is not None and abs(controlled_point) < ATTENUATION_FLOOR * abs(baseline):
            return {"verdict": "attenuated", "reason": "sign kept, interval clear of zero, under half the baseline"}
        return {"verdict": "survives", "reason": "sign kept, interval clear of zero, at least half the baseline"}
    if excludes_zero(random_control_ci):
        return {"verdict": "vanishes", "reason": "controlled interval crosses zero while the size-matched random control does not"}
    return {"verdict": "undecidable", "reason": "both the controlled and the size-matched random control cross zero: power-limited"}


# R6 every return path of claim_verdict, as (kwargs, expected verdict). Added after the
# independent verification round (results/phase3_independent_verification.csv) found that the
# 22 stored claims exercise only 7 of the 13 paths, so six were reachable but untested.
VERDICT_BRANCHES = [
    ({"baseline": 1.0, "baseline_ci": (0.5, 1.5), "controlled_ci": (0.4, 1.2), "random_control_ci": (0.5, 1.5),
      "controlled_point": 0.8, "precheck_blocking": True}, "undecidable"),
    ({"baseline": None, "baseline_ci": None, "controlled_ci": None, "random_control_ci": None,
      "controlled_point": None}, "out_of_instrument_scope"),
    ({"baseline": 1.0, "baseline_ci": (0.5, 1.5), "controlled_ci": (0.4, 1.2), "random_control_ci": (0.5, 1.5),
      "controlled_point": 0.8, "claim_direction": "null_no_preference"}, "undecidable"),
    ({"baseline": 0.1, "baseline_ci": (-0.3, 0.5), "controlled_ci": (0.2, 0.9), "random_control_ci": (-0.2, 0.4),
      "controlled_point": 0.55, "claim_direction": "null_no_preference"}, "reverses"),
    ({"baseline": 0.1, "baseline_ci": (-0.3, 0.5), "controlled_ci": (-0.2, 0.4), "random_control_ci": (-0.2, 0.4),
      "controlled_point": 0.1, "claim_direction": "null_no_preference"}, "survives"),
    ({"baseline": -1.0, "baseline_ci": (-1.5, -0.5), "controlled_ci": (-1.2, -0.4), "random_control_ci": (-1.5, -0.5),
      "controlled_point": -0.8}, "baseline_contradicts_claim"),
    ({"baseline": 0.2, "baseline_ci": (-0.3, 0.7), "controlled_ci": (-0.4, 0.6), "random_control_ci": (-0.4, 0.6),
      "controlled_point": 0.1, "baseline_reproduced": False}, "undecidable"),
    ({"baseline": 0.2, "baseline_ci": (-0.3, 0.7), "controlled_ci": (-0.4, 0.6), "random_control_ci": (-0.4, 0.6),
      "controlled_point": 0.1, "baseline_reproduced": True}, "undecidable"),
    ({"baseline": 1.0, "baseline_ci": (0.5, 1.5), "controlled_ci": (-1.2, -0.4), "random_control_ci": (0.5, 1.5),
      "controlled_point": -0.8}, "reverses"),
    ({"baseline": 1.0, "baseline_ci": (0.5, 1.5), "controlled_ci": (0.1, 0.4), "random_control_ci": (0.5, 1.5),
      "controlled_point": 0.3}, "attenuated"),
    ({"baseline": 1.0, "baseline_ci": (0.5, 1.5), "controlled_ci": (0.6, 1.3), "random_control_ci": (0.5, 1.5),
      "controlled_point": 0.9}, "survives"),
    ({"baseline": 1.0, "baseline_ci": (0.5, 1.5), "controlled_ci": (-0.1, 0.5), "random_control_ci": (0.6, 1.4),
      "controlled_point": 0.2}, "vanishes"),
    ({"baseline": 1.0, "baseline_ci": (0.5, 1.5), "controlled_ci": (-0.1, 0.5), "random_control_ci": (-0.2, 1.4),
      "controlled_point": 0.2}, "undecidable"),
]


# ------------------------------------------------------------------- regression
def _regression():
    out = {"cases": [], "failures": []}

    def record(name, expected, got, ok, note=""):
        out["cases"].append({"case": name, "expected": expected, "got": got, "ok": bool(ok), "note": note})
        if not ok:
            out["failures"].append(name)

    # R1 shares
    share_rows = read_csv(RESULTS / "ptm_detectability_share.csv")
    by_key = {}
    for r in share_rows:
        by_key.setdefault((r["dataset_id"], r["negative_construction"]), {})[r["feature_set"]] = r
    n_ok = n_tot = 0
    for (dataset, neg), sets in sorted(by_key.items()):
        vis, dig = sets.get("VIS10"), sets.get("DIG25")
        if not vis or not dig or vis["usable"] != "1" or not vis["share_visibility"]:
            continue
        got = detectability_share(float(vis["auc_global"]), float(dig["auc_global"]))
        n_tot += 1
        ok = got["defined"] and abs(got["share"] - float(vis["share_visibility"])) < 1e-3
        n_ok += ok
        if not ok:
            record(f"R1 share {dataset} {neg}", vis["share_visibility"],
                   None if not got["defined"] else round(got["share"], 4), False)
    record("R1 all stored shares recomputed", n_tot, n_ok, n_ok == n_tot and n_tot >= 9)

    # R2 depth-alone AUC from the stratum table
    strata = read_csv(RESULTS / "qtrp_depth_strata.csv")
    summary = read_csv(RESULTS / "qtrp_depth_conditioned_summary.csv")
    stored = {r["cohort"]: r for r in summary if r["design"] == "depth1"}
    for cohort in ("S1_pH5.0|raw", "S2_pH5.0|raw"):
        rows = [(r["depth"], r["n_positive"], r["n_negative"]) for r in strata if r["cohort"] == cohort]
        got = depth_confound_strength(rows)
        want = float(stored[cohort]["auc_plain"])
        record(f"R2 depth-alone AUC {cohort}", round(want, 4), round(got["depth_alone_auc"], 4),
               abs(got["depth_alone_auc"] - want) < 1.5e-3)
    for cohort in ("S1_pH5.0|arm_clean", "S2_pH5.0|arm_clean"):
        rows = [(r["depth"], r["n_positive"], r["n_negative"]) for r in strata if r["cohort"] == cohort]
        got = depth_confound_strength(rows)
        record(f"R2 arm_clean orphans positives beyond the limit ({cohort})", False, got["trustworthy"],
               got["trustworthy"] is False,
               f"orphaned {got['orphaned_positives']} = {got['orphaned_fraction']:.4f}")

    # R3 replay every claim verdict
    comb = read_csv(RESULTS / "phase2b_claim_retest_combined.csv")
    def f(x):
        try:
            return float(x)
        except (TypeError, ValueError):
            return None
    replayed = matched = 0
    for r in comb:
        base, bl, bh = f(r["baseline_log2_or"]), f(r["baseline_ci_low"]), f(r["baseline_ci_high"])
        ml, mh, mp = f(r["matched_ci_low"]), f(r["matched_ci_high"]), f(r["matched_log2_or"])
        rl, rh = f(r["random_control_ci_low"]), f(r["random_control_ci_high"])
        if None in (base, bl, bh, ml, mh, mp):
            continue
        replayed += 1
        got = claim_verdict(base, (bl, bh), (ml, mh), (rl, rh) if None not in (rl, rh) else None, mp,
                            claim_direction=r["claim_direction"] or "positive_preference",
                            baseline_reproduced=(r["baseline_reproduction"] == "reproduced"))
        ok = got["verdict"] == r["verdict"]
        matched += ok
        if not ok:
            record(f"R3 verdict {r['claim_id']}", r["verdict"], got["verdict"], False, got["reason"])
    record("R3 stored verdicts replayed by the rule", replayed, matched, matched == replayed and replayed >= 9)

    # R4 synthetic cleaning rules with known answers
    sys.path.insert(0, str(ROOT / "scripts"))
    import audit_cleaning_and_grouping as acg
    rng = np.random.default_rng(SEED)
    n = 4000
    covariate = rng.integers(1, 8, size=n).astype(float)
    y = (rng.random(n) < (covariate / 14.0)).astype(int)
    groups = np.array([f"P{i // 4}" for i in range(n)])
    strata_q = np.digitize(covariate, np.quantile(covariate, [0.2, 0.4, 0.6, 0.8]))
    keep_targeted = ~((y == 0) & (covariate >= 5))
    keep_random = np.ones(n, dtype=bool)
    keep_random[rng.choice(np.flatnonzero(y == 0), size=int(0.35 * (y == 0).sum()), replace=False)] = False
    targeted = acg.cleaning_collinearity(y, keep_targeted, covariate=covariate, strata=strata_q,
                                         groups=groups, label="synthetic: removes high-covariate negatives")
    uniform = acg.cleaning_collinearity(y, keep_random, covariate=covariate, strata=strata_q,
                                        groups=groups, label="synthetic: removes negatives uniformly")
    record("R4 targeted removal is blocked", True, bool(targeted["blocking_reasons"]),
           bool(targeted["blocking_reasons"]), str(targeted["blocking_reasons"])[:160])
    record("R4 uniform removal is not blocked", False, bool(uniform["blocking_reasons"]),
           not uniform["blocking_reasons"], str(uniform["blocking_reasons"])[:160])

    # R6 every return path of claim_verdict
    hit = set()
    for kwargs, want in VERDICT_BRANCHES:
        got = claim_verdict(**kwargs)
        hit.add(got["reason"])
        if got["verdict"] != want:
            record(f"R6 branch -> {want}", want, got["verdict"], False, got["reason"])
    record("R6 every claim_verdict return path covered", len(VERDICT_BRANCHES), len(hit),
           len(hit) == len(VERDICT_BRANCHES),
           "each case must land on a distinct reason string, so a missed path shows up as a collision")

    # R5 settled precedents, read from the recorded audit (not recomputed)
    prior = json.loads((RESULTS / "cleaning_and_grouping_audit.json").read_text(encoding="utf-8"))
    qtrp = prior["cases"]["qtrp_arm_clean_vs_detection_depth"]["verdict"]
    abpp = prior["cases"]["abpp_detection_matching_vs_labelling_intensity"]["verdict"]
    record("R5 recorded precedent: QTRP arm_clean demoted", "demote_to_sensitivity_check", qtrp,
           qtrp == "demote_to_sensitivity_check", "read from the stored audit, not recomputed")
    record("R5 recorded precedent: ABPP detection matching cleared", "usable_as_primary_caliber", abpp,
           abpp == "usable_as_primary_caliber", "read from the stored audit, not recomputed")
    return out


def main():
    started = time.time()
    reg = _regression()

    # apply the tier classifier to the claim corpus and to the census datasets
    claims = read_csv(RESULTS / "ptm_site_preference_claims.csv")
    tiers = {}
    for r in claims:
        t = label_semantics_tier(negative_class=r["negative_set_class"])["tier"]
        tiers[t] = tiers.get(t, 0) + 1

    standard_rows = []
    for tier, licence in TIER_LICENCE.items():
        classes = sorted(k for k, v in TIER_BY_NEGATIVE_CLASS.items() if v == tier)
        standard_rows.append({
            "tier": tier,
            "negative_semantics": "; ".join(classes) or "n/a",
            "what_the_negative_means": {
                "T1_direct_adduct": "this cysteine was seen in the same run carrying the other adduct form",
                "T2_parallel_ambiguous": "this cysteine was seen in the same run but the channel assignment is not exclusive at site level",
                "T3_capture_annotation": "this cysteine was not observed as modified",
                "unclassifiable": "undefined by the authors",
            }[tier],
            "claims_licensed": licence,
            "required_control": {
                "T1_direct_adduct": "none beyond the usual grouping and interval discipline",
                "T2_parallel_ambiguous": "a detectability control must be reported with the claim",
                "T3_capture_annotation": "no chemistry claim; report the detectability share instead",
                "unclassifiable": "define the negative set",
            }[tier],
            "worked_example_in_this_tree": {
                "T1_direct_adduct": "QTRP low-pH arms: positive = persulfide form detected, negative = free thiol only, same run",
                "T2_parallel_ambiguous": "qPerS-SID tier B donor arms; YAP1C sulfenylation negatives drawn from enrichment-run peptides",
                "T3_capture_annotation": "the v2 tomato training label: pulled down and alkylated versus every other cysteine",
                "unclassifiable": "6 of 51 extracted claims state no negative set at all",
            }[tier],
            "counterexample": {
                "T1_direct_adduct": "PXD005168 declares two channels but 112 of 187 direct-channel sites are also reported as blocked free thiols, so it is not T1",
                "T2_parallel_ambiguous": "a same-run negative set is not automatically safe: the YAP1C cohort's share went UP under NEG_B (0.8962 versus 0.7322)",
                "T3_capture_annotation": "a T3 dataset can still support a ranking claim; T3 is not 'unusable'",
                "unclassifiable": "n/a",
            }[tier],
            "claims_in_this_corpus": tiers.get(tier, 0),
        })
    with (RESULTS / "label_semantics_tiers.csv").open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(standard_rows[0]))
        w.writeheader()
        w.writerows(standard_rows)

    audit = {
        "script": "scripts/ptm_detectability_diagnostics.py",
        "script_sha256": sha256(SCRIPT),
        "precheck_module_sha256": sha256(ROOT / "scripts" / "audit_cleaning_and_grouping.py"),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "declared": {
            "SHARE_DENOMINATOR_FLOOR": SHARE_DENOMINATOR_FLOOR,
            "ORPHANED_LIMIT": ORPHANED_LIMIT,
            "ATTENUATION_FLOOR": ATTENUATION_FLOOR,
            "tier_map": TIER_BY_NEGATIVE_CLASS,
            "tier_licence": TIER_LICENCE,
        },
        "regression": {
            "cases": reg["cases"],
            "failures": reg["failures"],
            "all_pass": not reg["failures"],
            "share_interval_caliber": ("the 0.5963-0.9488 share range is the NEG_A (not-observed) caliber "
                                       "only; the unconditional range over both constructions runs to "
                                       "-0.9075 and must not be quoted as the share range"),
            "recomputed_versus_recorded": {
                "recomputed": ["R1 shares", "R2 depth-alone AUC and orphaning", "R3 verdict replay",
                               "R4 synthetic cleaning rules", "R6 verdict branch coverage"],
                "recorded_only": ["R5 QTRP arm_clean and ABPP detection-matching verdicts"],
            },
        },
        "tier_distribution_over_51_claims": tiers,
        "inputs": {p: sha256(RESULTS / p) for p in (
            "ptm_detectability_share.csv", "qtrp_depth_strata.csv", "qtrp_depth_conditioned_summary.csv",
            "phase2b_claim_retest_combined.csv", "ptm_site_preference_claims.csv",
            "cleaning_and_grouping_audit.json")},
        "outputs": {"results/label_semantics_tiers.csv": sha256(RESULTS / "label_semantics_tiers.csv")},
        "input_coupling": ("R3 reads results/phase2b_claim_retest_combined.csv, which a parallel round was "
                           "writing when the diagnostics were first run; that file is now final and is a "
                           "declared input, hashed below"),
        "known_limits": [
            "the share needs a visibility block and a superset; undefined otherwise rather than reported",
            "the tier classifier reads declared metadata and cannot detect a paper whose declared exclusive channels are not exclusive in its own table",
            "depth confound strength assumes the depth covariate does not encode the label",
            "claim_verdict is a verdict rule and takes effect estimates as given",
        ],
        "versions": {"python": sys.version, "numpy": np.__version__},
        "elapsed_seconds": round(time.time() - started, 2),
    }
    (RESULTS / "ptm_detectability_diagnostics_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    for c in reg["cases"]:
        print(f"  {'PASS' if c['ok'] else 'FAIL'}  {c['case']}: expected {c['expected']} got {c['got']} {c['note'][:90]}")
    print("all regression cases pass:", not reg["failures"])
    print("tier distribution over 51 claims:", tiers)
    print("elapsed", audit["elapsed_seconds"], "s")


if __name__ == "__main__":
    main()
