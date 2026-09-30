# -*- coding: utf-8 -*-
"""E_bootstrap_coverage -- POST HOC revision analysis (2026-09-30). Not registered, not pre-specified.

Task 4: multiplicity across the 28 re-tested claims.

The stored intervals are 95% percentile intervals; the underlying bootstrap replicates are not
stored, so adjusted intervals are APPROXIMATED by a normal approximation on the log2 scale, one
standard error per side:  SE_lo = (point - lo)/1.96,  SE_hi = (hi - point)/1.96,
adjusted interval = [point - z* SE_lo, point + z* SE_hi]. The two-sided P value for 0 uses the
side facing 0: p = 2(1 - Phi(|point|/SE_toward_0)) = 2(1 - Phi(1.96 w*)).

Families: F84 = 28 claims x (baseline, matched, random control) [the task's family];
F112 = F84 + the stratified (secondary-caliber) intervals; F28 = each interval type on its own.
Adjustments: Bonferroni (level 1 - 0.05/m) and Benjamini-Hochberg at q = 0.05 with
FCR-adjusted intervals (Benjamini & Yekutieli 2005: level 1 - R q / m for all intervals, R = number
of BH rejections; equivalent to 'excludes 0 iff BH-rejected'). The stratified interval, which the
verdict also reads, is adjusted at the level of the family it belongs to (in F84, where it is not a
member, it is adjusted at the same level as the matched interval of that family).
Each verdict is then re-derived with the pipeline's verbatim rule (claims_table.verdict_from_intervals).

Expected chance exclusions under a global null: sum over intervals of the per-interval type-I
error, nominal (0.05) and as simulated for its design and cluster count (sim_summary_long.csv).
Writes multiplicity_intervals.csv, multiplicity_verdicts.csv, multiplicity_summary.json.
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
Q = 0.05
Z = S.Z975
TYPES = ["baseline", "matched", "random", "stratified"]


def se_sides(p, lo, hi):
    return (p - lo) / Z, (hi - p) / Z


def pval(p, lo, hi):
    se_lo, se_hi = se_sides(p, lo, hi)
    se = se_lo if p >= 0 else se_hi
    if se <= 0:
        return 0.0
    return float(2 * (1 - sps.norm.cdf(abs(p) / se)))


def bh(pvals, q=Q):
    p = np.asarray(pvals, float)
    m = p.size
    order = np.argsort(p)
    ranked = p[order]
    below = ranked <= q * (np.arange(1, m + 1) / m)
    R = int(np.max(np.flatnonzero(below)) + 1) if below.any() else 0
    rej = np.zeros(m, dtype=bool)
    rej[order[:R]] = True
    return rej, R


def adjusted(p, lo, hi, level):
    zadj = sps.norm.ppf(1 - (1 - level) / 2)
    se_lo, se_hi = se_sides(p, lo, hi)
    return p - zadj * se_lo, p + zadj * se_hi


def row_with(r, levels):
    """Stored row -> row whose four intervals are rebuilt at the given per-type levels."""
    out = {"claim_direction": r["claim_direction"], "precheck_blocked": False}
    for t in TYPES:
        pk, lk, hk = INTERVALS[t]
        p, lo, hi = r[f"{t}_point"], r[f"{t}_lo"], r[f"{t}_hi"]
        out[pk] = p
        if levels[t] is None:
            out[lk], out[hk] = lo, hi
        else:
            out[lk], out[hk] = adjusted(p, lo, hi, levels[t])
    return out


def main():
    P = S.pipeline()
    ci = pd.read_csv(os.path.join(OUT, "claims_intervals.csv"))
    recs = []
    for _, r in ci.iterrows():
        for t in TYPES:
            recs.append({"claim_id": r.claim_id, "type": t, "point": r[f"{t}_point"], "lo": r[f"{t}_lo"],
                         "hi": r[f"{t}_hi"], "wstar": r[f"{t}_wstar"], "excludes0": bool(r[f"{t}_excludes0"]),
                         "decisive": bool(r[f"{t}_decisive"]),
                         "p_normal_approx": pval(r[f"{t}_point"], r[f"{t}_lo"], r[f"{t}_hi"]),
                         "unit": r.unit, "clusters": r.bootstrap_clusters})
    iv = pd.DataFrame(recs)
    fams = {"F84": ["baseline", "matched", "random"], "F112": TYPES}
    res_levels = {}
    for fam, types in fams.items():
        sub = iv[iv.type.isin(types)]
        m = len(sub)
        rej, R = bh(sub.p_normal_approx.to_numpy())
        iv.loc[sub.index, f"{fam}_bh_reject"] = rej
        iv.loc[sub.index, f"{fam}_bonf_reject"] = sub.p_normal_approx.to_numpy() <= Q / m
        res_levels[fam] = {"m": m, "bonferroni_level": 1 - Q / m, "bh_R": R,
                           "fcr_level": (1 - R * Q / m) if R > 0 else 1 - Q / m}
    for t in TYPES:
        sub = iv[iv.type == t]
        rej, R = bh(sub.p_normal_approx.to_numpy())
        iv.loc[sub.index, "F28_bh_reject"] = rej
        iv.loc[sub.index, "F28_bonf_reject"] = sub.p_normal_approx.to_numpy() <= Q / 28
        res_levels[f"F28_{t}"] = {"m": 28, "bonferroni_level": 1 - Q / 28, "bh_R": R,
                                  "fcr_level": (1 - R * Q / 28) if R > 0 else 1 - Q / 28}

    schemes = {
        "stored_95": {t: None for t in TYPES},
        "bonferroni_F84": {t: res_levels["F84"]["bonferroni_level"] for t in TYPES},
        "bh_fcr_F84": {t: res_levels["F84"]["fcr_level"] for t in TYPES},
        "bonferroni_F112": {t: res_levels["F112"]["bonferroni_level"] for t in TYPES},
        "bh_fcr_F112": {t: res_levels["F112"]["fcr_level"] for t in TYPES},
        "bonferroni_F28_per_type": {t: res_levels[f"F28_{t}"]["bonferroni_level"] for t in TYPES},
        "bh_fcr_F28_per_type": {t: res_levels[f"F28_{t}"]["fcr_level"] for t in TYPES},
    }
    vrows = []
    for _, r in ci.iterrows():
        rec = {"claim_id": r.claim_id, "unit": r.unit, "stored_verdict": r.stored_verdict}
        for name, lv in schemes.items():
            v, _, _ = verdict_from_intervals(row_with(r, lv), P)
            rec[name] = v
        vrows.append(rec)
    vd = pd.DataFrame(vrows)
    assert (vd.stored_95 == vd.stored_verdict).all()
    changes = {name: vd.loc[vd[name] != vd.stored_verdict, ["claim_id", "stored_verdict", name]].to_dict("records")
               for name in schemes if name != "stored_95"}
    tallies = {name: vd[name].value_counts().to_dict() for name in schemes}

    # expected chance exclusions under a global null
    summ_path = os.path.join(OUT, "sim_summary_long.csv")
    exp = {"nominal_F84": 84 * Q, "nominal_F112": 112 * Q,
           "p_at_least_one_nominal_F84_independent": 1 - (1 - Q) ** 84}
    if os.path.exists(summ_path):
        sm = pd.read_csv(summ_path)
        design_of = {"baseline": ("full", "pct"), "matched": ("sub", "pct"), "random": ("rc", "pct"),
                     "stratified": ("full_mh", "pct")}
        alphas = []
        for _, x in iv.iterrows():
            d, meth = design_of[x.type]
            scen = "null_iid" if (x.unit == "protein") else "null_re"
            g = sm[(sm.design == d) & (sm.method == meth) & (sm.scenario == scen)]
            Ks = np.array(sorted(g.K.unique()))
            kk = Ks[np.argmin(np.abs(np.log(Ks) - math.log(max(float(x.clusters), 1.0))))]
            alphas.append(float(g[g.K == kk].flag0.mean()))
        iv["simulated_type1"] = alphas
        a84 = iv[iv.type != "stratified"].simulated_type1.to_numpy()
        exp.update({"simulated_F84": float(a84.sum()), "simulated_F112": float(iv.simulated_type1.sum()),
                    "p_at_least_one_simulated_F84_independent": float(1 - np.prod(1 - a84))})
    obs = {"observed_exclusions_F84": int(iv[iv.type != "stratified"].excludes0.sum()),
           "observed_exclusions_F112": int(iv.excludes0.sum()),
           "observed_exclusions_by_type": iv.groupby("type").excludes0.sum().astype(int).to_dict()}
    iv.to_csv(os.path.join(OUT, "multiplicity_intervals.csv"), index=False)
    vd.to_csv(os.path.join(OUT, "multiplicity_verdicts.csv"), index=False)
    out = {"label": "POST HOC revision analysis E_bootstrap_coverage, task 4 (not registered)",
           "approximation": "normal approximation on the log2 scale from stored 95% percentile bounds, one SE per side",
           "levels": res_levels, "expected_chance_exclusions_global_null": exp, "observed": obs,
           "verdict_changes": changes, "tallies": tallies}
    with open(os.path.join(OUT, "multiplicity_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1, default=str)
    print(json.dumps({k: v for k, v in out.items() if k not in ("tallies",)}, indent=1, default=str))


if __name__ == "__main__":
    main()
