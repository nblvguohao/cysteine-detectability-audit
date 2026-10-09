# -*- coding: utf-8 -*-
"""E_bootstrap_coverage -- REVISION AFTER VERIFICATION (round 1), 2026-09-30.
POST HOC revision analysis; not registered, not pre-specified.

1. Task 3 re-derived with the SINGLE-DRAW yardstick as the PRIMARY calibration of the random-control
   interval (verifier problem 1; reasoning in r1_rc_yardsticks.py). Percentile intervals (baseline,
   matched, stratified) are rescaled by their w0 exactly as in calibrated_verdicts.py
   (calibrated_intervals.csv); the random-control interval is rescaled by
     primary     w_single = w0_sub * width(sub) / width(rc)  (calibrated width of ONE same-size draw)
     secondary   w0_rc    (nominal coverage of its own 20-draw-average point; the first report's choice)
     also        w_uncal  = width(sub) / width(rc)           (the percentile interval of one draw)
   Source of w_single: claim-matched runs (SFE-008: r1_sfe008_yardsticks.json, primary variant
   re_nb077, which matches both mean cluster size and singleton fraction; others:
   r1_rc_yardsticks_claims.csv); otherwise the grid (r1_rc_decomposition_grid.csv: nearest K on the
   log scale, folded prevalence, null_iid for protein-level and null_re for site-level claims, mean over
   the two size families; range = bootstrap CIs over both null scenarios and both size families).
   Each verdict is re-derived with the pipeline's rule (claims_table.verdict_from_intervals) at the
   primary value and at both ends of the range.
2. SNO-004 (null claim on its authors' cohort, stored null_broken_by_control): the verdict under the
   registered bucket map (NULL_VERDICT_BUCKET, run_phase2b_claims_backfill.py) and under the released
   tool's claim_verdict, and the post hoc adjusted readings WITH their null resolution (half-width of
   the matched interval, the quantity the phase-2b rules require to be quoted with a surviving null),
   and the concordant stratified interval (verifier problem 2).
3. The list of verdicts whose deciding bound is near 0 (verifier problem 6b).
Writes r1_calibrated_verdicts.csv, r1_rc_intervals.csv, r1_sno004_sensitivity.json, r1_summary.json.
"""
from __future__ import annotations

import json
import math
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
from scipy import stats as sps

import simlib as S
from claims_table import INTERVALS, verdict_from_intervals

OUT = r"/path/to/local/_cys_repo_work/public/revision_2026-09-30/results/E_bootstrap_coverage"
TYPES = ["baseline", "matched", "stratified", "random"]
NEGA = {"SNO-001", "SNO-004", "SNO-012", "SNO-014"}
Z = S.Z975
SFE008_PRIMARY = "re_nb077"


def rescale(p, lo, hi, w):
    return p - w * (p - lo), p + w * (hi - p)


def rc_single_source(claim, unit, K, pi, grid, cm, sfe):
    if claim == "SFE-008":
        v = sfe["variants"]
        prim = v[SFE008_PRIMARY]
        los = [x["w_single_lo"] for x in v.values()]
        his = [x["w_single_hi"] for x in v.values()]
        first = cm[(cm.claim == "SFE-008")]
        los += list(first.w_single_lo)
        his += list(first.w_single_hi)
        return {"w_single": prim["w_single"], "w_single_lo": min(los), "w_single_hi": max(his),
                "w_uncal": prim["w_uncal"], "w0_rc": prim["w0_rc"],
                "source": f"claim-matched r1_sfe008 {SFE008_PRIMARY} (range: CIs of re_nb2, re_nb077, iid_nb077 and the first run)"}
    g = cm[cm.claim == claim]
    if len(g):
        var = "iid" if unit == "protein" else ("nega" if claim in NEGA and (g.variant == "nega").any() else "re")
        prim = g[g.variant == var].iloc[0] if (g.variant == var).any() else g.iloc[0]
        return {"w_single": float(prim.w_single), "w_single_lo": float(g.w_single_lo.min()),
                "w_single_hi": float(g.w_single_hi.max()), "w_uncal": float(prim.w_uncal), "w0_rc": float(prim.w0_rc),
                "source": f"claim-matched {prim.variant} (range over {','.join(sorted(g.variant))})"}
    scen = "null_iid" if unit == "protein" else "null_re"
    Ks = np.array(sorted(grid.K.unique()))
    kk = int(Ks[np.argmin(np.abs(np.log(Ks) - math.log(max(K, 1))))])
    pis = np.array(sorted(grid.pi_x.unique()))
    pp = float(pis[np.argmin(np.abs(pis - min(pi, 1 - pi)))])
    cell = grid[(grid.K == kk) & (grid.pi_x == pp)]
    prim = cell[cell.scenario == scen]
    return {"w_single": float(prim.w_single.mean()), "w_single_lo": float(cell.w_single_lo.min()),
            "w_single_hi": float(cell.w_single_hi.max()), "w_uncal": float(prim.w_uncal.mean()),
            "w0_rc": float(prim.w0_rc.mean()), "source": f"grid {scen} K={kk} pi={pp} (range over null_iid,null_re x nb,lognormal)"}


def main():
    P = S.pipeline()
    ci = pd.read_csv(os.path.join(OUT, "claims_intervals.csv"))
    cal = pd.read_csv(os.path.join(OUT, "calibrated_intervals.csv"))
    grid = pd.read_csv(os.path.join(OUT, "r1_rc_decomposition_grid.csv"))
    cm = pd.read_csv(os.path.join(OUT, "r1_rc_yardsticks_claims.csv"))
    sfe = json.load(open(os.path.join(OUT, "r1_sfe008_yardsticks.json"), encoding="utf-8"))
    rc_rows, vrows = [], []
    for _, r in ci.iterrows():
        pi = (r.n_attribute / r.n_observations) if pd.notna(r.n_attribute) else 0.35
        K = float(r.bootstrap_clusters)
        w0 = {t: cal[(cal.claim_id == r.claim_id) & (cal.type == t)].iloc[0] for t in TYPES}
        rs = rc_single_source(r.claim_id, r.unit, K, pi, grid, cm, sfe)
        p, lo, hi = r.random_point, r.random_lo, r.random_hi
        rc_rows.append({"claim_id": r.claim_id, "point": p, "lo": lo, "hi": hi, "wstar": r.random_wstar,
                        "decisive": bool(r.random_decisive), "verdict_if_toggled": r.random_verdict_if_toggled,
                        **rs, "w0_rc_first_report": float(w0["random"].w0),
                        "near_bound_at_w_single": (p - rs["w_single"] * (p - lo)) if p >= 0 else (p + rs["w_single"] * (hi - p)),
                        "near_bound_at_w_single_lo": (p - rs["w_single_lo"] * (p - lo)) if p >= 0 else (p + rs["w_single_lo"] * (hi - p)),
                        "near_bound_at_w_single_hi": (p - rs["w_single_hi"] * (p - lo)) if p >= 0 else (p + rs["w_single_hi"] * (hi - p)),
                        "near_bound_at_w_uncal": (p - rs["w_uncal"] * (p - lo)) if p >= 0 else (p + rs["w_uncal"] * (hi - p)),
                        "flips_at_w_single": bool((r.random_wstar > 1) != (r.random_wstar > rs["w_single"])),
                        "could_flip_in_range": bool(min(1.0, rs["w_single_lo"]) < r.random_wstar < max(1.0, rs["w_single_hi"]))})

        def verdict(w_pct, w_rc):
            row = {"claim_direction": r.claim_direction, "precheck_blocked": False}
            for t in TYPES:
                pk, lk, hk = INTERVALS[t]
                pt, l_, h_ = r[f"{t}_point"], r[f"{t}_lo"], r[f"{t}_hi"]
                w = w_rc if t == "random" else w_pct[t]
                row[pk] = pt
                row[lk], row[hk] = rescale(pt, l_, h_, w)
            return verdict_from_intervals(row, P)[0]

        one = {t: 1.0 for t in TYPES}
        pct_pt = {t: float(w0[t].w0) for t in TYPES}
        pct_lo = {t: float(w0[t].w0_lo) for t in TYPES}
        pct_hi = {t: float(w0[t].w0_hi) for t in TYPES}
        assert verdict(one, 1.0) == r.stored_verdict, r.claim_id
        vrows.append({"claim_id": r.claim_id, "unit": r.unit, "is_transfer": bool(r.is_transfer),
                      "baseline_reproduction": r.baseline_reproduction, "stored_verdict": r.stored_verdict,
                      # primary: percentile intervals at w0, random control at the single-draw yardstick
                      "verdict_primary": verdict(pct_pt, rs["w_single"]),
                      "verdict_primary_low_end": verdict(pct_lo, rs["w_single_lo"]),
                      "verdict_primary_high_end": verdict(pct_hi, rs["w_single_hi"]),
                      "verdict_rc_single_draw_only": verdict(one, rs["w_single"]),
                      "verdict_rc_single_draw_lo_only": verdict(one, rs["w_single_lo"]),
                      "verdict_rc_single_draw_hi_only": verdict(one, rs["w_single_hi"]),
                      "verdict_rc_uncalibrated_single_draw_only": verdict(one, rs["w_uncal"]),
                      "verdict_secondary_rc_nominal": verdict(pct_pt, rs["w0_rc"]),
                      "verdict_pct_only_rc_as_stored": verdict(pct_pt, 1.0),
                      "w_single": rs["w_single"], "w_single_lo": rs["w_single_lo"], "w_single_hi": rs["w_single_hi"],
                      "rc_source": rs["source"]})
    rcd = pd.DataFrame(rc_rows)
    vd = pd.DataFrame(vrows)
    vd["flips_primary"] = vd.verdict_primary != vd.stored_verdict
    vd["could_flip_primary_range"] = vd.flips_primary | (vd.verdict_primary_low_end != vd.stored_verdict) | (
        vd.verdict_primary_high_end != vd.stored_verdict)
    vd["flips_secondary_rc_nominal"] = vd.verdict_secondary_rc_nominal != vd.stored_verdict
    # reading of the calibration analysis for each verdict:
    #   stable    the verdict is the stored one at the central value and at both ends of the range
    #   boundary  the two ends of the range give different verdicts (the analysis cannot decide it)
    #   flips     the central value and both ends agree on a verdict other than the stored one
    ends_differ = vd.verdict_primary_low_end != vd.verdict_primary_high_end
    all_other = (vd.verdict_primary != vd.stored_verdict) & (vd.verdict_primary_low_end != vd.stored_verdict) & (
        vd.verdict_primary_high_end != vd.stored_verdict)
    vd["reading"] = np.where(ends_differ, "boundary", np.where(all_other, "flips", np.where(
        vd.could_flip_primary_range, "boundary", "stable")))
    rcd.to_csv(os.path.join(OUT, "r1_rc_intervals.csv"), index=False)
    vd.to_csv(os.path.join(OUT, "r1_calibrated_verdicts.csv"), index=False)

    # ---- SNO-004 sensitivity (verifier problem 2) -------------------------------------------------
    s = ci[ci.claim_id == "SNO-004"].iloc[0]
    msum = json.load(open(os.path.join(OUT, "multiplicity_summary.json"), encoding="utf-8"))
    lv = msum["levels"]

    def adj(t, level):
        p, lo, hi = s[f"{t}_point"], s[f"{t}_lo"], s[f"{t}_hi"]
        z = sps.norm.ppf(1 - (1 - level) / 2)
        a, b = p - z * (p - lo) / Z, p + z * (hi - p) / Z
        return {"interval": [float(a), float(b)], "half_width_null_resolution": float((b - a) / 2),
                "excludes0": bool(a > 0 or b < 0), "max_odds_ratio_in_interval": float(2 ** b)}

    def pval(t):
        p, lo, hi = s[f"{t}_point"], s[f"{t}_lo"], s[f"{t}_hi"]
        se = (p - lo) / Z if p >= 0 else (hi - p) / Z
        return float(2 * sps.norm.sf(abs(p) / se))

    from cys_audit.verdict import claim_verdict, VERDICT_TO_CLAIM
    tool = claim_verdict(s.baseline_point, [s.baseline_lo, s.baseline_hi], [s.matched_lo, s.matched_hi],
                         [s.random_lo, s.random_hi], s.matched_point, s.claim_direction)
    bucket = P["NULL_VERDICT_BUCKET"]
    sno = {"label": "POST HOC sensitivity only (not registered); the registered reading is the stored one",
           "claim": "SNO-004", "claim_direction": s.claim_direction, "baseline_reproduction": s.baseline_reproduction,
           "is_transfer": bool(s.is_transfer), "stored_verdict": s.stored_verdict,
           "registered_bucket_map": bucket, "registered_bucket_of_stored_verdict": bucket.get(s.stored_verdict),
           "released_tool_claim_verdict": tool, "released_tool_status": VERDICT_TO_CLAIM[tool["verdict"]],
           "stored_intervals": {t: [float(s[f"{t}_point"]), float(s[f"{t}_lo"]), float(s[f"{t}_hi"])] for t in TYPES},
           "stored_null_resolution_matched_half_width": float((s.matched_hi - s.matched_lo) / 2),
           "wstar": {t: float(s[f"{t}_wstar"]) for t in TYPES},
           "p_normal_approx": {t: pval(t) for t in ("baseline", "matched", "stratified")},
           "adjusted_normal_approx": {
               "bonferroni_F84": {t: adj(t, lv["F84"]["bonferroni_level"]) for t in ("matched", "stratified")},
               "bh_fcr_F84": {t: adj(t, lv["F84"]["fcr_level"]) for t in ("matched", "stratified")},
               "bonferroni_F112": {t: adj(t, lv["F112"]["bonferroni_level"]) for t in ("matched", "stratified")},
               "bh_fcr_F112": {t: adj(t, lv["F112"]["fcr_level"]) for t in ("matched", "stratified")}},
           "calibration_w0_matched_range": [float(cal[(cal.claim_id == "SNO-004") & (cal.type == "matched")].w0_lo.iloc[0]),
                                            float(cal[(cal.claim_id == "SNO-004") & (cal.type == "matched")].w0_hi.iloc[0])],
           "calibration_w0_stratified_range": [float(cal[(cal.claim_id == "SNO-004") & (cal.type == "stratified")].w0_lo.iloc[0]),
                                               float(cal[(cal.claim_id == "SNO-004") & (cal.type == "stratified")].w0_hi.iloc[0])]}
    with open(os.path.join(OUT, "r1_sno004_sensitivity.json"), "w", encoding="utf-8") as fh:
        json.dump(sno, fh, indent=1, default=str)

    # ---- near-zero list (verifier problem 6b) ------------------------------------------------------
    near = []
    for _, r in ci.iterrows():
        for t in TYPES:
            if not r[f"{t}_decisive"]:
                continue
            w = r[f"{t}_wstar"]
            c = cal[(cal.claim_id == r.claim_id) & (cal.type == t)].iloc[0]
            if t == "random":
                rr = rcd[rcd.claim_id == r.claim_id].iloc[0]
                lo_, hi_ = rr.w_single_lo, rr.w_single_hi
            else:
                lo_, hi_ = c.w0_lo, c.w0_hi
            within = min(1.0, lo_) <= w <= max(1.0, hi_)
            if within or c.z_mc < 3:
                near.append({"claim_id": r.claim_id, "interval": t, "stored": [r[f"{t}_point"], r[f"{t}_lo"], r[f"{t}_hi"]],
                             "wstar": w, "yardstick_range": [lo_, hi_], "z_mc": c.z_mc,
                             "verdict_if_bound_crossed": r[f"{t}_verdict_if_toggled"], "stored_verdict": r.stored_verdict})
    summ = {"label": "POST HOC revision analysis E_bootstrap_coverage, revision after verification (not registered)",
            "primary_yardstick_random_control": "single same-size draw (calibrated width)",
            "reading_counts": vd.reading.value_counts().to_dict(),
            "boundary_verdicts": vd.loc[vd.reading == "boundary", ["claim_id", "stored_verdict", "verdict_primary",
                                                                   "verdict_primary_low_end", "verdict_primary_high_end"]].to_dict("records"),
            "flips_all_readings": vd.loc[vd.reading == "flips", ["claim_id", "stored_verdict", "verdict_primary"]].to_dict("records"),
            "flips_primary_central_value": vd.loc[vd.flips_primary, ["claim_id", "stored_verdict", "verdict_primary"]].to_dict("records"),
            "could_flip_primary_range": vd.loc[vd.could_flip_primary_range, ["claim_id", "stored_verdict", "verdict_primary",
                                                                             "verdict_primary_low_end", "verdict_primary_high_end"]].to_dict("records"),
            "flips_secondary_rc_nominal": vd.loc[vd.flips_secondary_rc_nominal, ["claim_id", "stored_verdict",
                                                                                  "verdict_secondary_rc_nominal"]].to_dict("records"),
            "tally_stored": vd.stored_verdict.value_counts().to_dict(),
            "tally_primary": vd.verdict_primary.value_counts().to_dict(),
            "own_data_claims_changing_to_vanishes_any_yardstick": vd.loc[(~vd.is_transfer) & (
                (vd.verdict_primary == "vanishes") | (vd.verdict_primary_low_end == "vanishes") | (vd.verdict_primary_high_end == "vanishes")
                | (vd.verdict_secondary_rc_nominal == "vanishes")) & (vd.stored_verdict != "vanishes"), "claim_id"].tolist(),
            "near_zero_deciding_bounds": near}
    with open(os.path.join(OUT, "r1_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summ, fh, indent=1, default=str)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 80)
    print(rcd[rcd.decisive][["claim_id", "point", "lo", "hi", "wstar", "w_single", "w_single_lo", "w_single_hi", "w_uncal",
                             "w0_rc", "near_bound_at_w_single", "near_bound_at_w_single_lo", "near_bound_at_w_single_hi",
                             "near_bound_at_w_uncal", "flips_at_w_single", "could_flip_in_range", "source"]].round(4).to_string())
    print(vd[["claim_id", "stored_verdict", "verdict_primary", "verdict_primary_low_end", "verdict_primary_high_end",
              "verdict_rc_uncalibrated_single_draw_only", "verdict_secondary_rc_nominal"]].to_string())
    print(json.dumps({k: v for k, v in summ.items() if k != "near_zero_deciding_bounds"}, indent=1, default=str))
    print(json.dumps(summ["near_zero_deciding_bounds"], indent=1, default=str))
    print(json.dumps(sno, indent=1, default=str))


if __name__ == "__main__":
    main()
