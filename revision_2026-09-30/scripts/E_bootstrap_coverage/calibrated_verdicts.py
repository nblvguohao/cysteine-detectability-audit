# -*- coding: utf-8 -*-
"""E_bootstrap_coverage -- POST HOC revision analysis (2026-09-30). Not registered, not pre-specified.

Task 3: which stored verdicts could the simulated miscalibration flip?

For every interval a verdict reads (baseline, matched, stratified, random control) this script takes
a calibration multiplier w0 = 95th percentile, under H0, of |point| / half-width-toward-0 (so that a
calibrated interval excludes 0 in 5% of null datasets). w0 > 1: the percentile interval is too
narrow; w0 < 1: too wide. Source of w0, in order of preference:
  1. the claim-matched null simulation of that claim and design (claims_matched_null_summary.csv;
     NEG_A variant for the NEG_A SNO cohorts, 're' as sensitivity, iid for protein-level claims);
  2. the grid (sim_summary_long.csv): design (baseline->full, matched->sub, stratified->full_mh,
     random->rc), scenario null_re (site-level) or null_iid (protein-level), nearest cluster count
     on the log scale (largest grid K = 1475 for larger cohorts), nearest attribute prevalence after
     folding p > 0.5 to 1 - p (the log OR only changes sign), mean over the two size families.
A stored interval is rescaled about its point estimate, lo' = p - w0 (p - lo), hi' = p + w0 (hi - p),
and the verdict re-derived with the pipeline's verbatim rule. Sensitivity: the same with w0 at the
ends of its range (Monte Carlo 95% CI, variants and size families). A decisive bound whose flip
multiplier w* lies inside that range is reported as 'could flip'.

Also reported: Monte Carlo error of each stored bound from B = 5000 replicates,
SD(bound) ~ sqrt(0.025*0.975/B) / phi(1.96) * (hi - lo)/3.92  (normal bootstrap distribution),
and z_MC = |bound| / SD: bounds with z_MC < 3 could move across 0 with a different bootstrap seed.
Writes calibrated_verdicts.csv, calibrated_intervals.csv, calibrated_summary.json.
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

OUT = r"C:/Users/admin/Desktop/小论文/_cys_repo_work/public/revision_2026-09-30/results/E_bootstrap_coverage"
TYPES = ["baseline", "matched", "stratified", "random"]
DESIGN = {"baseline": "full", "matched": "sub", "stratified": "full_mh", "random": "rc"}
CM_DESIGN = {"baseline": "full", "matched": "sub", "stratified": "mh", "random": "rc"}
NEGA = {"SNO-001", "SNO-004", "SNO-012", "SNO-014"}
B = 5000
MC_FACTOR = math.sqrt(0.025 * 0.975 / B) / sps.norm.pdf(1.959963984540054)


def grid_w0(sm, design, unit, K, pi):
    scen = "null_iid" if unit == "protein" else "null_re"
    g = sm[(sm.design == design) & (sm.method == "pct") & (sm.scenario == scen)]
    Ks = np.array(sorted(g.K.unique()))
    kk = int(Ks[np.argmin(np.abs(np.log(Ks) - math.log(max(K, 1))))])
    pis = np.array(sorted(g.pi_x.unique()))
    pf = min(pi, 1 - pi)
    pp = float(pis[np.argmin(np.abs(pis - pf))])
    gg = g[(g.K == kk) & (g.pi_x == pp)]
    # plausible range across both size families and all non-confounded + mildly confounded null scenarios
    allg = sm[(sm.design == design) & (sm.method == "pct") & (sm.K == kk) & (sm.pi_x == pp)
              & (sm.scenario.isin(["null_iid", "null_re", "null_conf15"]))]
    return {"w0": float(gg.w0.mean()), "w0_lo": float(min(gg.w0_ci_lo.min(), allg.w0.min())),
            "w0_hi": float(max(gg.w0_ci_hi.max(), allg.w0.max())),
            "type1": float(gg.flag0.mean()), "source": f"grid {design} {scen} K={kk} pi={pp}"}


def cm_w0(cm, claim, design, unit):
    g = cm[(cm.claim == claim) & (cm.design == design)]
    if g.empty:
        return None
    if unit == "protein":
        prim = g[g.variant == "iid"]
    elif claim in NEGA:
        prim = g[g.variant == "nega"]
    else:
        prim = g[g.variant == "re"]
    if prim.empty:
        prim = g
    return {"w0": float(prim.w0.iloc[0]), "w0_lo": float(g.w0_ci_lo.min()), "w0_hi": float(g.w0_ci_hi.max()),
            "type1": float(prim.type1.iloc[0]),
            "source": f"claim-matched {design} {prim.variant.iloc[0]} (range over {','.join(sorted(g.variant))})"}


def rc_as_subsample(sm, cm, claim, unit, K, pi):
    """Alternative yardstick for the random-control interval: the width a calibrated interval of ONE
    same-size random subsample would have (sub design, scaled by its own w0), relative to the
    random-control interval's width: w_alt = (width_sub * w0_sub) / width_rc."""
    g = cm[cm.claim == claim] if len(cm) else cm
    if len(g) and {"sub", "rc"} <= set(g.design):
        var = "iid" if unit == "protein" else ("nega" if claim in NEGA and (g.variant == "nega").any() else "re")
        gg = g[g.variant == var] if (g.variant == var).any() else g
        sub = gg[gg.design == "sub"].iloc[0]
        rc = gg[gg.design == "rc"].iloc[0]
        return float(sub.mean_width * sub.w0 / rc.mean_width), f"claim-matched {var}"
    scen = "null_iid" if unit == "protein" else "null_re"
    base = sm[(sm.method == "pct") & (sm.scenario == scen)]
    Ks = np.array(sorted(base.K.unique()))
    kk = int(Ks[np.argmin(np.abs(np.log(Ks) - math.log(max(K, 1))))])
    pis = np.array(sorted(base.pi_x.unique()))
    pp = float(pis[np.argmin(np.abs(pis - min(pi, 1 - pi)))])
    b = base[(base.K == kk) & (base.pi_x == pp)]
    sub = b[b.design == "sub"]
    rc = b[b.design == "rc"]
    return float((sub.mean_width * sub.w0).mean() / rc.mean_width.mean()), f"grid {scen} K={kk} pi={pp}"


def rescale(p, lo, hi, w):
    return p - w * (p - lo), p + w * (hi - p)


def _without_rc(verdict_at, cal):
    """Calibrated verdict with every interval rescaled except the random control (kept as stored)."""
    saved = dict(cal["random"])
    cal["random"].update({"w0": 1.0, "w0_lo": 1.0, "w0_hi": 1.0})
    v = verdict_at("point")
    cal["random"].clear()
    cal["random"].update(saved)
    return v


def main():
    P = S.pipeline()
    ci = pd.read_csv(os.path.join(OUT, "claims_intervals.csv"))
    sm = pd.read_csv(os.path.join(OUT, "sim_summary_long.csv"))
    cm_path = os.path.join(OUT, "claims_matched_null_summary.csv")
    cm = pd.read_csv(cm_path) if os.path.exists(cm_path) else pd.DataFrame(columns=["claim"])
    ivrows, vrows = [], []
    for _, r in ci.iterrows():
        pi = (r.n_attribute / r.n_observations) if pd.notna(r.n_attribute) else 0.35
        K = float(r.bootstrap_clusters)
        cal = {}
        for t in TYPES:
            src = cm_w0(cm, r.claim_id, CM_DESIGN[t], r.unit) if len(cm) else None
            if src is None:
                src = grid_w0(sm, DESIGN[t], r.unit, K, pi)
            cal[t] = src
            p, lo, hi = r[f"{t}_point"], r[f"{t}_lo"], r[f"{t}_hi"]
            bound = lo if p >= 0 else hi
            mc_sd = MC_FACTOR * (hi - lo) / 3.919927969080108
            ivrows.append({"claim_id": r.claim_id, "type": t, "point": p, "lo": lo, "hi": hi,
                           "excludes0": bool(r[f"{t}_excludes0"]), "decisive": bool(r[f"{t}_decisive"]),
                           "verdict_if_toggled": r[f"{t}_verdict_if_toggled"], "wstar": r[f"{t}_wstar"],
                           "w0": src["w0"], "w0_lo": src["w0_lo"], "w0_hi": src["w0_hi"],
                           "simulated_type1": src["type1"], "w0_source": src["source"],
                           "flips_at_w0": bool((r[f"{t}_wstar"] > 1) != (r[f"{t}_wstar"] > src["w0"])),
                           "could_flip_in_range": bool(min(1.0, src["w0_lo"]) < r[f"{t}_wstar"] <= max(1.0, src["w0_hi"])
                                                       or min(1.0, src["w0_lo"]) <= r[f"{t}_wstar"] < max(1.0, src["w0_hi"])),
                           "near_bound": bound, "bound_mc_sd_B5000": mc_sd,
                           "z_mc": abs(bound) / mc_sd if mc_sd > 0 else np.inf})
        w_rc_alt, rc_alt_src = rc_as_subsample(sm, cm, r.claim_id, r.unit, K, pi)

        def verdict_at(which):
            row = {"claim_direction": r.claim_direction, "precheck_blocked": False}
            for t in TYPES:
                pk, lk, hk = INTERVALS[t]
                p, lo, hi = r[f"{t}_point"], r[f"{t}_lo"], r[f"{t}_hi"]
                if which == "rc_alt":
                    w = w_rc_alt if t == "random" else 1.0
                else:
                    w = {"point": cal[t]["w0"], "low": cal[t]["w0_lo"], "high": cal[t]["w0_hi"], "one": 1.0}[which]
                row[pk] = p
                row[lk], row[hk] = rescale(p, lo, hi, w)
            return verdict_from_intervals(row, P)[0]
        v1 = verdict_at("one")
        assert v1 == r.stored_verdict, (r.claim_id, v1, r.stored_verdict)
        vrows.append({"claim_id": r.claim_id, "unit": r.unit, "clusters": int(K), "stored_verdict": r.stored_verdict,
                      "calibrated_verdict_w0": verdict_at("point"),
                      "verdict_w0_low_end": verdict_at("low"), "verdict_w0_high_end": verdict_at("high"),
                      "calibrated_verdict_w0_random_control_unchanged": _without_rc(verdict_at, cal),
                      "w_rc_alt": w_rc_alt, "w_rc_alt_source": rc_alt_src,
                      "verdict_rc_as_subsample": verdict_at("rc_alt")})
    ivd = pd.DataFrame(ivrows)
    vd = pd.DataFrame(vrows)
    vd["flips_at_w0"] = vd.calibrated_verdict_w0 != vd.stored_verdict
    vd["could_flip_in_range"] = (vd.verdict_w0_low_end != vd.stored_verdict) | (vd.verdict_w0_high_end != vd.stored_verdict) | vd.flips_at_w0
    ivd.to_csv(os.path.join(OUT, "calibrated_intervals.csv"), index=False)
    vd.to_csv(os.path.join(OUT, "calibrated_verdicts.csv"), index=False)
    out = {"label": "POST HOC revision analysis E_bootstrap_coverage, task 3 (not registered)",
           "flips_at_w0": vd.loc[vd.flips_at_w0, ["claim_id", "stored_verdict", "calibrated_verdict_w0"]].to_dict("records"),
           "could_flip_in_range": vd.loc[vd.could_flip_in_range, ["claim_id", "stored_verdict", "verdict_w0_low_end",
                                                                  "verdict_w0_high_end"]].to_dict("records"),
           "decisive_bounds_within_3_mc_sd": ivd.loc[ivd.decisive & (ivd.z_mc < 3),
                                                     ["claim_id", "type", "near_bound", "bound_mc_sd_B5000", "z_mc"]].to_dict("records"),
           "tally_stored": vd.stored_verdict.value_counts().to_dict(),
           "tally_calibrated_w0": vd.calibrated_verdict_w0.value_counts().to_dict()}
    with open(os.path.join(OUT, "calibrated_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1, default=str)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 60)
    print(ivd[ivd.decisive][["claim_id", "type", "wstar", "w0", "w0_lo", "w0_hi", "simulated_type1", "flips_at_w0",
                             "could_flip_in_range", "z_mc", "w0_source"]].round(3).to_string())
    print(vd.to_string())
    print(json.dumps({k: v for k, v in out.items() if k.startswith("tally") or k.startswith("flips")}, indent=1, default=str))


if __name__ == "__main__":
    main()
