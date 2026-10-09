"""Independent re-count of every number in the phase-3 deliverables, from the bottom tables up.

This is a cross-check of work done by another session, not a re-run of it.  Rule for this
round, declared before anything is computed: a number is CONFIRMED only when it is recounted
from a table that is UPSTREAM of the artefact being checked, and then compared against BOTH
sides separately - the audit JSON and the prose report - because a report can agree with its
own audit and both still be wrong about the tree.

WHAT IS RECOUNTED AND FROM WHERE (nothing here reads the phase-3 audit for its truth value)
-------------------------------------------------------------------------------------------
V1  detectability share: recomputed from results/ptm_detectability_share.csv per-feature-set
    global AUCs by the declared formula (AUC_VIS10 - 0.5) / (AUC_DIG25 - 0.5).  Counts usable
    shares, the min/max range, and the NEG_A -> NEG_B drop on the cohorts that have both.
V2  QTRP detection depth: depth is re-derived from the NORMALISED SITE TABLE
    results/qtrp_sites_normalised.csv by the definition predeclared in
    scripts/run_qtrp_depth_conditioned.py - "the number of distinct human QTRP arms
    (source_table, ph) in which the site appears at all, in either adduct form", over rows
    whose species is human and whose label (class_from_adduct, else class_from_sheet) is SH or
    SSH.  Labels per cohort are read from the site-level score table
    results/qtrp_site_scores.csv (arm, negative_definition, label), NOT from the stratum table.
    The stratum table results/qtrp_depth_strata.csv is then checked cell by cell against this
    reconstruction, and the depth-alone AUC is recomputed with sklearn.metrics.roc_auc_score -
    a different implementation from both the original round (scipy rankdata) and the
    diagnostics package (hand-rolled tie-corrected ranks).
V3  arm_clean orphaning: positives sitting in single-class depth strata, recounted from the
    same reconstruction, as a count and as a fraction of that cohort's positives.
V4  verdict rule: replayed over results/phase2_claim_retest.csv - the SETTLED 8-claim phase-2
    product (survives 5 / vanishes 1 / undecidable 1 / out-of-instrument-scope 1) - and the
    four precedents this round was told to reproduce are checked numerically against that
    table: SFE-006 baseline 1.1559 [0.9601, 1.3591] -> matched 0.1761 [-0.0814, 0.4379] with a
    size-matched random control 1.1826 [0.8993, 1.4720]; PERS-005/011/012 matched retention.
V5  label semantics tiers: the 51-claim distribution is recounted from
    results/ptm_site_preference_claims.csv with a mapping typed out independently in this file
    (TIER_EXPECTED below), and the report's reconciliation sentence - that phase 1's "30
    claims" is the three not-observed classes, a narrower set than the 39 - is recounted too.
V6  reference benchmark: the four cohorts' sites / positives / grouping units are recounted
    from their own upstream sources (results/v3_qtrp_stack_run_audit.json for the pooled v3
    cohort, results/ptm_detectability_census_audit.json for the census cohorts) and, for the
    two QTRP arms, from results/qtrp_site_scores.csv as well, which is upstream of both.

READING RULE: exact integer equality for counts; 1e-3 absolute for AUCs and shares (the stored
values are rounded to 4 decimals); a mismatch is reported as MISMATCH and never reconciled by
adjusting the expectation.  Discrepancies of provenance - two different rounds' numbers carried
under one cohort label - are reported as PROVENANCE even when both numbers are individually
correct, because that is what a reader would misread.

Terminology: sulfenylation = 次磺酰化 (-SOH); sulfinylation = 亚磺酰化 (-SO2H).
"""
import collections
import csv
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
SCRIPT = Path(__file__).resolve()
COUNT_TOL = 0
NUMERIC_TOL = 1e-3

# typed out here independently of scripts/ptm_detectability_diagnostics.py
TIER_EXPECTED = {
    "same_run_mutually_exclusive": "T1_direct_adduct",
    "detected_background_proteome": "T2_parallel_ambiguous",
    "whole_proteome_annotation_background": "T3_capture_annotation",
    "all_cys_in_identified_proteins": "T3_capture_annotation",
    "unmodified_cys_same_protein": "T3_capture_annotation",
    "curated_negatives_other_db": "T3_capture_annotation",
    "random_or_simulated": "T3_capture_annotation",
    "not_defined": "unclassifiable",
}
NOT_OBSERVED_CLASSES = ("whole_proteome_annotation_background",
                        "all_cys_in_identified_proteins",
                        "unmodified_cys_same_protein")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


FINDINGS = []


def check(item, recounted, side_audit, side_report, kind="count", note=""):
    """Compare one recounted truth against the audit side and the report side separately."""
    def agree(a, b):
        if a is None or b is None:
            return None
        if kind == "count":
            return int(a) == int(b)
        return abs(float(a) - float(b)) <= NUMERIC_TOL
    rec = {"item": item, "recounted": recounted, "audit_side": side_audit,
           "report_side": side_report, "audit_agrees": agree(recounted, side_audit),
           "report_agrees": agree(recounted, side_report), "kind": kind, "note": note}
    bad = [s for s in ("audit_agrees", "report_agrees") if rec[s] is False]
    rec["status"] = "MISMATCH" if bad else ("CONFIRMED" if rec["audit_agrees"] or rec["report_agrees"] else "UNCHECKED")
    FINDINGS.append(rec)
    return rec


def main():
    started = time.time()
    phase3_audit = json.loads((RESULTS / "ptm_detectability_diagnostics_audit.json").read_text(encoding="utf-8"))
    bench_audit = json.loads((RESULTS / "reference_benchmark_spec_audit.json").read_text(encoding="utf-8"))
    cases = {c["case"]: c for c in phase3_audit["regression"]["cases"]}

    # ---------------------------------------------------------------- V1 shares
    share_rows = read_csv(RESULTS / "ptm_detectability_share.csv")
    by_key = collections.defaultdict(dict)
    for r in share_rows:
        by_key[(r["dataset_id"], r["negative_construction"])][r["feature_set"]] = r
    shares = {}
    for key, sets in by_key.items():
        vis, dig = sets.get("VIS10"), sets.get("DIG25")
        if not vis or not dig or vis["usable"] != "1" or not vis["share_visibility"]:
            continue
        denom = num(dig["auc_global"]) - 0.5
        if denom <= 0.02:
            continue
        shares[key] = {"recomputed": (num(vis["auc_global"]) - 0.5) / denom,
                       "stored": num(vis["share_visibility"])}
    n_share = len(shares)
    max_dev = max(abs(v["recomputed"] - v["stored"]) for v in shares.values())
    check("V1 usable shares recomputed from stored per-feature-set AUCs", n_share,
          num(cases["R1 all stored shares recomputed"]["expected"]), 16,
          note=f"max |recomputed - stored| = {max_dev:.2e}")
    vals = sorted(v["stored"] for v in shares.values())
    # The 0.5963-0.9488 range quoted in the diagnostics docstring is CONDITIONAL on
    # NEG_A ("when negatives are 'not pulled down'"), so it must be recounted on the NEG_A
    # cohorts only. The first pass of this file checked it unconditionally and reported a
    # MISMATCH at -0.9075; that was this file's scoping error, not an error in the artefact,
    # and it is recorded here rather than silently deleted.
    neg_a = sorted(v["stored"] for k, v in shares.items() if k[1] == "NEG_A_not_observed")
    neg_b = sorted(v["stored"] for k, v in shares.items() if k[1] != "NEG_A_not_observed")
    check("V1 NEG_A share range minimum", round(neg_a[0], 4), None, 0.5963, kind="numeric",
          note=f"n_NEG_A={len(neg_a)}")
    check("V1 NEG_A share range maximum", round(neg_a[-1], 4), None, 0.9488, kind="numeric",
          note=f"n_NEG_A={len(neg_a)}")
    FINDINGS.append({
        "item": "V1 unconditional share range over all 16 usable cohorts",
        "recounted": f"[{vals[0]:.4f}, {vals[-1]:.4f}]", "audit_side": None, "report_side": None,
        "audit_agrees": None, "report_agrees": None, "kind": "context", "status": "NOTE",
        "note": f"NEG_B range [{neg_b[0]:.4f}, {neg_b[-1]:.4f}]; the minimum is qtrp_S2_ph5 "
                "NEG_B at -0.9075, i.e. visibility features rank INVERTED there. No document "
                "states the unconditional range; the quoted 0.5963-0.9488 is NEG_A-only and correct."})
    paired = {}
    for (ds, neg), v in shares.items():
        paired.setdefault(ds, {})[neg] = v["stored"]
    drops = [p["NEG_A_not_observed"] - p["NEG_B_observed_unmodified"]
             for p in paired.values() if {"NEG_A_not_observed", "NEG_B_observed_unmodified"} <= set(p)]
    check("V1 mean NEG_A -> NEG_B share drop", round(float(np.mean(drops)), 4), None, 0.4645,
          kind="numeric", note=f"n_paired_cohorts={len(drops)}")

    # ------------------------------------------------- V2 depth, rebuilt from the site table
    seen_arms = collections.defaultdict(set)
    for row in read_csv(RESULTS / "qtrp_sites_normalised.csv"):
        if row.get("species") != "human":
            continue
        label = row["class_from_adduct"] or row["class_from_sheet"]
        if label not in ("SH", "SSH"):
            continue
        seen_arms[(row["accession"], int(row["site"]))].add((row["source_table"], row["ph"]))
    depth = {k: len(v) for k, v in seen_arms.items()}

    scores = read_csv(RESULTS / "qtrp_site_scores.csv")
    stored_strata = read_csv(RESULTS / "qtrp_depth_strata.csv")
    summary = {(r["cohort"], r["design"]): r for r in read_csv(RESULTS / "qtrp_depth_conditioned_summary.csv")}
    cohort_facts = {}
    for arm in ("S1_pH5.0", "S2_pH5.0"):
        for negdef in ("raw", "arm_clean"):
            cohort = f"{arm}|{negdef}"
            rows = [r for r in scores if r["arm"] == arm and r["negative_definition"] == negdef]
            keys = [(r["accession"], int(r["site"])) for r in rows]
            missing = [k for k in keys if k not in depth]
            y = np.asarray([int(r["label"]) for r in rows])
            d = np.asarray([depth.get(k, -1) for k in keys], dtype=float)
            tab = collections.Counter()
            for di, yi in zip(d, y):
                tab[(int(di), int(yi))] += 1
            recount = {int(v): {"n_positive": tab[(int(v), 1)], "n_negative": tab[(int(v), 0)]}
                       for v in sorted(set(d))}
            stored = {int(r["depth"]): {"n_positive": int(r["n_positive"]), "n_negative": int(r["n_negative"])}
                      for r in stored_strata if r["cohort"] == cohort}
            cell_mismatch = [(k, recount.get(k), stored.get(k))
                             for k in sorted(set(recount) | set(stored))
                             if recount.get(k) != stored.get(k)]
            auc_here = float(roc_auc_score(y, d)) if 0 < y.sum() < len(y) else float("nan")
            orph = sum(v["n_positive"] for v in recount.values() if v["n_positive"] and not v["n_negative"])
            cohort_facts[cohort] = {
                "n_sites": len(rows), "n_positive": int(y.sum()), "n_negative": int((y == 0).sum()),
                "sites_without_depth": len(missing), "strata": len(recount),
                "informative_strata": sum(1 for v in recount.values() if v["n_positive"] and v["n_negative"]),
                "depth_alone_auc_sklearn": round(auc_here, 4),
                "orphaned_positives": orph,
                "orphaned_fraction": round(orph / max(1, int(y.sum())), 4),
                "stratum_cell_mismatches": cell_mismatch,
            }
    for cohort in ("S1_pH5.0|raw", "S2_pH5.0|raw"):
        f = cohort_facts[cohort]
        stored_auc = num(summary[(cohort, "depth1")]["auc_plain"])
        reg = cases[f"R2 depth-alone AUC {cohort}"]
        check(f"V2 depth-alone AUC {cohort} (sklearn, depth rebuilt from the site table)",
              f["depth_alone_auc_sklearn"], num(reg["got"]), stored_auc, kind="numeric",
              note=f"stored round summary auc_plain={stored_auc}; regression expected {reg['expected']}")
        check(f"V2 stratum table cells {cohort}", 0, len(f["stratum_cell_mismatches"]), 0,
              note="number of (depth, class) cells where the rebuilt count differs from qtrp_depth_strata.csv")
        check(f"V2 cohort positives {cohort}", f["n_positive"], None,
              272 if cohort.startswith("S1") else 80)
    for cohort in ("S1_pH5.0|arm_clean", "S2_pH5.0|arm_clean"):
        f = cohort_facts[cohort]
        reg = cases[f"R2 arm_clean orphans positives beyond the limit ({cohort})"]
        reported = {"S1_pH5.0|arm_clean": (102, 0.3750), "S2_pH5.0|arm_clean": (48, 0.6000)}[cohort]
        check(f"V3 orphaned positives {cohort}", f["orphaned_positives"], None, reported[0],
              note=f"regression note: {reg['note']}")
        check(f"V3 orphaned fraction {cohort}", f["orphaned_fraction"], None, reported[1], kind="numeric")

    # ------------------------------------------------------- V4 verdict rule on the settled 8
    sys.path.insert(0, str(ROOT / "scripts"))
    from ptm_detectability_diagnostics import claim_verdict
    claims_meta = {r["claim_id"]: r for r in read_csv(RESULTS / "ptm_site_preference_claims.csv")}
    phase2 = read_csv(RESULTS / "phase2_claim_retest.csv")
    replay = {"n": 0, "match": 0, "rows": []}
    for r in phase2:
        base, bl, bh = num(r["baseline_log2_or"]), num(r["baseline_ci_low"]), num(r["baseline_ci_high"])
        ml, mh, mp = num(r["matched_ci_low"]), num(r["matched_ci_high"]), num(r["matched_log2_or"])
        rl, rh = num(r["random_control_ci_low"]), num(r["random_control_ci_high"])
        direction = (claims_meta.get(r["claim_id"], {}).get("claim_direction") or "positive_preference")
        if None in (base, bl, bh, ml, mh, mp):
            replay["rows"].append({"claim_id": r["claim_id"], "stored": r["verdict"], "replayed": "incomplete_numbers"})
            continue
        got = claim_verdict(base, (bl, bh), (ml, mh), (rl, rh) if None not in (rl, rh) else None, mp,
                            claim_direction=direction,
                            baseline_reproduced=(r["baseline_reproduced"] in ("1", "True", "true", "reproduced")))
        replay["n"] += 1
        replay["match"] += got["verdict"] == r["verdict"]
        replay["rows"].append({"claim_id": r["claim_id"], "stored": r["verdict"], "replayed": got["verdict"],
                               "direction": direction})
    verdict_counts = collections.Counter(r["verdict"] for r in phase2)
    check("V4 settled phase-2 claims replayed by claim_verdict", replay["n"], None, replay["match"],
          note=f"stored verdict counts {dict(verdict_counts)}")
    sfe = {r["claim_id"]: r for r in phase2}.get("SFE-006", {})
    for field, want in (("baseline_log2_or", 1.1559), ("baseline_ci_low", 0.9601), ("baseline_ci_high", 1.3591),
                        ("matched_log2_or", 0.1761), ("matched_ci_low", -0.0814), ("matched_ci_high", 0.4379),
                        ("random_control_log2_or", 1.1826), ("random_control_ci_low", 0.8993),
                        ("random_control_ci_high", 1.4720)):
        check(f"V4 SFE-006 {field}", num(sfe.get(field)), None, want, kind="numeric")
    for cid, want in (("PERS-005", 0.832), ("PERS-011", 0.922), ("PERS-012", 0.716)):
        row = {r["claim_id"]: r for r in phase2}.get(cid, {})
        base, mp = num(row.get("baseline_log2_or")), num(row.get("matched_log2_or"))
        check(f"V4 {cid} matched retention (matched / baseline)", round(mp / base, 3) if base else None,
              None, want, kind="numeric",
              note=f"verdict {row.get('verdict')} matched {mp} baseline {base} "
                   f"matched CI [{row.get('matched_ci_low')}, {row.get('matched_ci_high')}]")

    # --------------------------------------------------------------- V5 tier distribution
    claim_rows = list(claims_meta.values())
    tiers = collections.Counter(TIER_EXPECTED.get(r["negative_set_class"], "unclassifiable") for r in claim_rows)
    audit_tiers = phase3_audit["tier_distribution_over_51_claims"]
    check("V5 claims in the corpus", len(claim_rows), sum(audit_tiers.values()), 51)
    for tier, reported in (("T3_capture_annotation", 39), ("unclassifiable", 6),
                           ("T2_parallel_ambiguous", 3), ("T1_direct_adduct", 3)):
        check(f"V5 tier {tier}", tiers[tier], audit_tiers.get(tier), reported)
    not_observed = sum(1 for r in claim_rows if r["negative_set_class"] in NOT_OBSERVED_CLASSES)
    check("V5 phase-1 'not observed' claims (the narrower 30)", not_observed, None, 30,
          note="report section three says the 39 is the wider set and phase 1's 30 is the three not-observed classes")
    tier_csv = read_csv(RESULTS / "label_semantics_tiers.csv")
    check("V5 label_semantics_tiers.csv rows", len(tier_csv), None, 4)
    csv_counts = {r["tier"]: int(r["claims_in_this_corpus"]) for r in tier_csv}
    for tier in tiers:
        check(f"V5 tier csv column {tier}", tiers[tier], csv_counts.get(tier), None)

    # ------------------------------------------------------------- V6 reference benchmark
    v3 = json.loads((RESULTS / "v3_qtrp_stack_run_audit.json").read_text(encoding="utf-8"))
    census = json.loads((RESULTS / "ptm_detectability_census_audit.json").read_text(encoding="utf-8"))["datasets"]
    spec = {r["cohort"]: r for r in read_csv(RESULTS / "reference_benchmark_spec.csv")}
    pooled = v3["cohorts"]["v3_pooled_ph5_raw"]
    check("V6 v3_pooled_ph5_raw sites", pooled["n_sites"], num(spec["v3_pooled_ph5_raw"]["n_sites"]), 1991)
    check("V6 v3_pooled_ph5_raw positives", pooled["n_positive"], num(spec["v3_pooled_ph5_raw"]["n_positive"]), 314)
    check("V6 v3_pooled_ph5_raw grouping units (homology components)", pooled["n_components"],
          num(spec["v3_pooled_ph5_raw"]["n_grouping_units"]), 1205)
    for cohort, cen_key, rep in (("S1_pH5.0|raw", "qtrp_S1_ph5", (1435, 277, 203)),
                                 ("S2_pH5.0|raw", "qtrp_S2_ph5", (1185, 80, 57)),
                                 ("qpers_sid_tierB", "qpers_sid_tierB", (2018, 724, 510))):
        src = census[cen_key]
        sites = src["n_positive"] + src["n_neg_B_observed_unmodified"]
        check(f"V6 {cohort} sites (census NEG_B)", sites, num(spec[cohort]["n_sites"]), rep[0])
        check(f"V6 {cohort} positives (census)", src["n_positive"], num(spec[cohort]["n_positive"]), rep[1])
        check(f"V6 {cohort} grouping units (census proteins WITH a positive)",
              src["n_proteins_with_positive"], num(spec[cohort]["n_grouping_units"]), rep[2])
    comp_rows = read_csv(RESULTS / "v3_human_homology_components.csv")
    check("V6 homology component file rows", len(comp_rows), bench_audit["components_file_rows"], 6574)
    eligible_flags = {c: spec[c]["eligible"] for c in spec}
    check("V6 cohorts marked eligible", sum(1 for v in eligible_flags.values() if v == "True"),
          len(bench_audit["eligible_cohorts"]), 4)
    identifiability_filled = sum(1 for c in spec if spec[c]["identifiable_units_in_within_caliber"].strip())
    check("V6 cohorts with the within-caliber identifiability count actually filled in",
          identifiability_filled, None, None,
          note="declared eligibility has three criteria (tier, >=30 identifiable units, homology folds) "
               "but the eligible flag is computed from tier alone")

    # ------------------------------------------- V7 gaps that are not arithmetic errors
    def gap(item, status, recounted, note):
        FINDINGS.append({"item": item, "recounted": recounted, "audit_side": None,
                         "report_side": None, "audit_agrees": None, "report_agrees": None,
                         "kind": "structural", "status": status, "note": note})

    # V7a branch coverage of the verdict rule under the regression that tests it
    combined = read_csv(RESULTS / "phase2b_claim_retest_combined.csv")
    branches = collections.Counter()
    for r in combined:
        base, bl, bh = num(r["baseline_log2_or"]), num(r["baseline_ci_low"]), num(r["baseline_ci_high"])
        ml, mh, mp = num(r["matched_ci_low"]), num(r["matched_ci_high"]), num(r["matched_log2_or"])
        rl, rh = num(r["random_control_ci_low"]), num(r["random_control_ci_high"])
        if None in (base, bl, bh, ml, mh, mp):
            continue
        g = claim_verdict(base, (bl, bh), (ml, mh), (rl, rh) if None not in (rl, rh) else None, mp,
                          claim_direction=r["claim_direction"] or "positive_preference",
                          baseline_reproduced=(r["baseline_reproduction"] == "reproduced"))
        branches[(g["verdict"], g["reason"][:40])] += 1
    n_return_paths = 13  # counted by reading claim_verdict
    gap("V7a claim_verdict return paths exercised by the R3 replay", "GAP",
        f"{len(branches)} of {n_return_paths}",
        "never exercised: precheck_blocking -> undecidable; baseline None -> "
        "out_of_instrument_scope; the two sign-flip 'reverses' paths; 'null survives'; and "
        "'baseline interval crosses zero' with the baseline reproduced. The completeness "
        "filter in R3 removes exactly the rows that would exercise the out-of-scope path, so "
        "that branch is shipped untested. Fixable with hand-built cases, no new data needed. "
        f"branches fired: {sorted(str(k) for k in branches)}")
    n_settled = len(read_csv(RESULTS / "phase2_claim_retest.csv"))
    settled = {r["claim_id"]: r["verdict"] for r in read_csv(RESULTS / "phase2_claim_retest.csv")}
    comb_v = {r["claim_id"]: r["verdict"] for r in combined}
    disagree = {c: (settled[c], comb_v.get(c)) for c in settled if comb_v.get(c) != settled[c]}
    gap("V7b R3 reads results/phase2b_claim_retest_combined.csv, written by a parallel round",
        "NOTE" if not disagree else "MISMATCH", f"{len(combined)} rows, 14 with complete numbers",
        f"the settled {n_settled}-claim verdicts are identical in both files "
        f"(disagreements: {disagree or 'none'}), so the package's regression does not contradict "
        "the settled product - but the regression is now coupled to a file that was being "
        "written concurrently with it. Pinning results/phase2_claim_retest.csv as a second, "
        "frozen R3 input would decouple them.")

    # V7c-V7f reference-benchmark provenance and rule-implementation gaps
    depth_round = {}
    for arm in ("S1_pH5.0", "S2_pH5.0"):
        rows = [r for r in scores if r["arm"] == arm and r["negative_definition"] == "raw"]
        depth_round[f"{arm}|raw"] = {
            "n_sites": len(rows), "n_positive": sum(int(r["label"]) for r in rows),
            "n_proteins": len({r["accession"] for r in rows}),
            "n_proteins_with_positive": len({r["accession"] for r in rows if r["label"] == "1"})}
    gap("V7c cohort label collision in results/reference_benchmark_spec.csv", "PROVENANCE",
        json.dumps(depth_round, ensure_ascii=False),
        "the rows named S1_pH5.0|raw and S2_pH5.0|raw carry CENSUS-round counts "
        "(277 positives / 1435 sites; 80 / 1185, negatives = NEG_B observed-unmodified), but in "
        "the rest of this tree those exact labels denote the DEPTH-round cohorts "
        "(272 / 1404; 80 / 1168). Every number is individually correct and traceable; the two "
        "site sets are different and share one name. S2 positives happen to coincide at 80, "
        "which makes the collision easy to miss.")
    gap("V7d n_grouping_units mixes two different quantities", "PROVENANCE",
        "v3: 1205 homology components; S1/S2/qpers: 203/57/510 proteins WITH a positive",
        "the column - and the report table header 分组单元 - carries homology components for the "
        "v3 row and proteins-with-a-positive for the other three. For comparison the depth-round "
        f"S1|raw cohort has {depth_round['S1_pH5.0|raw']['n_proteins']} proteins in total and "
        f"{depth_round['S1_pH5.0|raw']['n_proteins_with_positive']} with a positive, so the "
        "census's 203 is a third quantity again. Two columns, or a units column, would fix it.")
    gap("V7e the eligible flag implements one of the three declared criteria", "GAP",
        f"eligible=True for 4 of 4 cohorts; identifiability count filled for {identifiability_filled} of 4",
        "declared_eligibility names tier, at least 30 identifiable within-caliber units, and "
        "homology-component folds. In the code eligible = tier in (T1, T2); the identifiability "
        "criterion is evaluated only for v3_pooled_ph5_raw (44 units) and the folds criterion is "
        "met only there too - the spec's own grouping column reads 'protein / homology component' "
        "for S1/S2 and 'protein' for qpers_sid_tierB. Three of four cohorts are therefore marked "
        "eligible on the tier criterion alone, which the spec does not say.")
    ident_source = read_csv(RESULTS / "v3_qtrp_stack_summary.csv")
    ident_vals = sorted({r["components_contributing_to_within_caliber"] for r in ident_source
                         if r["cohort"] == "v3_pooled_ph5_raw"})
    gap("V7f the 44 is correct but not hash-pinned by the script that prints it", "GAP",
        f"v3_qtrp_stack_summary.csv components_contributing_to_within_caliber = {ident_vals} "
        f"for v3_pooled_ph5_raw (all designs)",
        "reports/V3_QTRP_TRAINED_MODEL.md states the same 44 against the pre-declared 30 "
        "threshold, so the number is sourced - but build_reference_benchmark_spec.py hardcodes it "
        "and results/v3_qtrp_stack_summary.csv is absent from that audit's inputs hash list. "
        "One line to read it instead.")
    docstring_item = json.loads((RESULTS / "qtrp_depth_confound_report_audit.json").read_text(
        encoding="utf-8")).get("docstring_versus_artifact_discrepancy")
    gap("V7g depth-alone AUC 0.6833/0.7880 versus 0.6828/0.7885", "NOTE",
        "bottom-up recount lands on 0.6828 / 0.7885 (the artefact side)",
        "the pre-declared table in the run_qtrp_depth_conditioned.py docstring lists 0.6833 and "
        "0.7880; its own summary output lists 0.6828 and 0.7885. This round's depth rebuilt from "
        "qtrp_sites_normalised.csv reproduces the ARTEFACT numbers, so the phase-3 report's "
        "'consistent with the original round' is accurate against the artefact. The docstring gap "
        "was already recorded upstream in results/qtrp_depth_confound_report_audit.json: "
        f"{json.dumps(docstring_item, ensure_ascii=False)[:220] if docstring_item else 'key absent'}")

    out = {
        "script": "scripts/verify_phase3_diagnostic_and_standard.py",
        "script_sha256": sha256(SCRIPT),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "purpose": "independent bottom-table re-count of the phase-3 deliverables produced by another session",
        "reading_rule": {"counts": "exact integer equality", "numeric": f"absolute {NUMERIC_TOL}"},
        "checked_artefacts": {p: sha256(RESULTS / p) for p in (
            "ptm_detectability_diagnostics_audit.json", "label_semantics_tiers.csv",
            "reference_benchmark_spec.csv", "reference_benchmark_spec_audit.json")},
        "upstream_inputs": {p: sha256(RESULTS / p) for p in (
            "ptm_detectability_share.csv", "qtrp_sites_normalised.csv", "qtrp_site_scores.csv",
            "qtrp_depth_strata.csv", "qtrp_depth_conditioned_summary.csv", "phase2_claim_retest.csv",
            "ptm_site_preference_claims.csv", "v3_qtrp_stack_run_audit.json",
            "ptm_detectability_census_audit.json", "v3_human_homology_components.csv")},
        "checked_script_sha256": {
            "scripts/ptm_detectability_diagnostics.py": sha256(ROOT / "scripts" / "ptm_detectability_diagnostics.py"),
            "scripts/build_reference_benchmark_spec.py": sha256(ROOT / "scripts" / "build_reference_benchmark_spec.py"),
        },
        "qtrp_cohort_recount": cohort_facts,
        "phase2_verdict_replay": replay,
        "tier_recount": dict(tiers),
        "share_recount": {"n_usable": n_share, "min": round(vals[0], 4), "max": round(vals[-1], 4),
                          "max_abs_deviation_from_stored": max_dev,
                          "paired_cohorts_for_drop": len(drops),
                          "mean_drop": round(float(np.mean(drops)), 4)},
        "findings": FINDINGS,
        "counts": {
            "rows": len(FINDINGS),
            **{s.lower(): sum(1 for f in FINDINGS if f["status"] == s)
               for s in ("CONFIRMED", "MISMATCH", "UNCHECKED", "NOTE", "GAP", "PROVENANCE")},
            "numeric_checks_against_a_stated_value": sum(
                1 for f in FINDINGS if f["status"] in ("CONFIRMED", "MISMATCH", "UNCHECKED")),
        },
        "branch_coverage_of_claim_verdict": {"fired": len(branches), "return_paths": n_return_paths},
        "versions": {"python": sys.version, "numpy": np.__version__},
        "elapsed_seconds": round(time.time() - started, 2),
    }
    (RESULTS / "phase3_independent_verification_audit.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    with (RESULTS / "phase3_independent_verification.csv").open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["item", "kind", "recounted", "audit_side", "report_side",
                                           "audit_agrees", "report_agrees", "status", "note"])
        w.writeheader()
        for f in FINDINGS:
            w.writerow(f)
    for f in FINDINGS:
        print(f"  {f['status']:9s} {f['item']}: recounted={f['recounted']} audit={f['audit_side']} "
              f"report={f['report_side']} {f['note'][:80]}")
    print("counts:", out["counts"])


if __name__ == "__main__":
    main()
