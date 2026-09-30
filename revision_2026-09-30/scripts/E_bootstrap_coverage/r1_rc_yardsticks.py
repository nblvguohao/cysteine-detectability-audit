# -*- coding: utf-8 -*-
"""E_bootstrap_coverage -- REVISION AFTER VERIFICATION (round 1), 2026-09-30.
POST HOC revision analysis; not registered, not pre-specified.

Why this script exists (verifier problems 1 and 3). The pipeline's same-size random control
(random_control_effect) takes 20 fixed-seed draws of as many random negatives as positives; its point
estimate is the MEAN of the 20 draw estimates, while each of its 5,000 bootstrap replicates picks ONE
draw at random and resamples that draw's proteins. Two yardsticks can be used to say whether that
interval is 'too wide':

  nominal      w0_rc   = 95th percentile under H0 of |rc point| / rc half-width toward 0: the factor
                          that gives the rc interval 95% coverage of ITS OWN 20-draw-average point
                          estimate. Because the 20-draw average has close to baseline precision, an
                          rc interval rescaled this way re-runs the baseline test (the degenerate case
                          the released tool guards against, verdict.py lines 206-209).
  single draw  w_single = w0_sub * mean width(sub) / mean width(rc): the factor that gives the rc
                          interval the width of a CALIBRATED interval of ONE same-size random draw
                          (sub design = all positives + as many random negatives), i.e. the precision
                          of 'a smaller sample', which is what the control is declared to represent
                          (Experimental Procedures: 'a loss of power in a smaller sample'). PRIMARY.
  (also)       w_uncal  = mean width(sub) / mean width(rc): the width of the percentile interval that
                          one same-size draw actually has (the released tool's control is one draw).

Decomposition of the over-width (verifier problem 3):
  mixing factor  = mean width(rc) / mean width(sub)   (draw mixing inside the replicates)
  point SD ratio = SD(rc point) / SD(sub point)        (averaging 20 draws in the point estimate)
  1/w0_rc        = over-width relative to the rc's own point (both effects together)
  1/w_single     = over-width relative to one calibrated same-size draw

Inputs (all this item's own outputs): sim_cells/cell_*.csv.gz, sim_cell_parameters.csv,
claims_matched_null_sim.csv.gz. Uncertainty: 95% percentile interval from 2,000 bootstrap resamples of
datasets (paired: sub and rc of the same dataset), seed 20260932.
Writes r1_rc_decomposition_grid.csv, r1_rc_decomposition_by_K.csv, r1_rc_yardsticks_claims.csv.
"""
from __future__ import annotations

import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd

import simlib as S

OUT = r"C:/Users/admin/Desktop/小论文/_cys_repo_work/public/revision_2026-09-30/results/E_bootstrap_coverage"
SEED = 20260932
NBOOT = 2000


def q95(r):
    return float(np.quantile(r, 0.95))


def yardsticks(sub_pt, sub_lo, sub_hi, rc_pt, rc_lo, rc_hi, rng, nboot=NBOOT, full_pt=None):
    sub_pt, sub_lo, sub_hi = (np.asarray(v, float) for v in (sub_pt, sub_lo, sub_hi))
    rc_pt, rc_lo, rc_hi = (np.asarray(v, float) for v in (rc_pt, rc_lo, rc_hi))
    ok = np.isfinite(sub_lo) & np.isfinite(sub_hi) & np.isfinite(rc_lo) & np.isfinite(rc_hi)
    sub_pt, sub_lo, sub_hi, rc_pt, rc_lo, rc_hi = (v[ok] for v in (sub_pt, sub_lo, sub_hi, rc_pt, rc_lo, rc_hi))
    r_sub = S.ratio_to_zero(sub_pt, sub_lo, sub_hi)
    r_rc = S.ratio_to_zero(rc_pt, rc_lo, rc_hi)
    w_sub, w_rc = sub_hi - sub_lo, rc_hi - rc_lo

    def stats(ix):
        ws, wr = w_sub[ix].mean(), w_rc[ix].mean()
        w0s, w0r = q95(r_sub[ix]), q95(r_rc[ix])
        return np.array([w0s, w0r, w0s * ws / wr, ws / wr, wr / ws,
                         np.std(rc_pt[ix], ddof=1) / np.std(sub_pt[ix], ddof=1)])

    n = len(sub_pt)
    est = stats(np.arange(n))
    boots = np.array([stats(rng.integers(0, n, n)) for _ in range(nboot)])
    lo, hi = np.percentile(boots, [2.5, 97.5], axis=0)
    names = ["w0_sub", "w0_rc", "w_single", "w_uncal", "mixing_factor", "point_sd_ratio"]
    out = {"n_datasets": int(n), "type1_sub": float(np.mean((sub_lo > 0) | (sub_hi < 0))),
           "type1_rc": float(np.mean((rc_lo > 0) | (rc_hi < 0))),
           "mean_width_sub": float(w_sub.mean()), "mean_width_rc": float(w_rc.mean()),
           "sd_point_sub": float(np.std(sub_pt, ddof=1)), "sd_point_rc": float(np.std(rc_pt, ddof=1))}
    if full_pt is not None:
        fp = np.asarray(full_pt, float)[ok]
        out["sd_point_full"] = float(np.std(fp, ddof=1))
        out["sd_full_over_sub"] = out["sd_point_full"] / out["sd_point_sub"]
    for i, nm in enumerate(names):
        out[nm], out[f"{nm}_lo"], out[f"{nm}_hi"] = float(est[i]), float(lo[i]), float(hi[i])
    out["overwidth_vs_own_point"] = 1.0 / out["w0_rc"]
    out["overwidth_vs_single_draw"] = 1.0 / out["w_single"]
    return out


def grid(rng):
    par = pd.read_csv(os.path.join(OUT, "sim_cell_parameters.csv"))
    rows = []
    for _, c in par[par.scenario.isin(["null_iid", "null_re"])].iterrows():
        d = pd.read_csv(os.path.join(OUT, "sim_cells", f"cell_{int(c.cell_id):03d}.csv.gz"),
                        usecols=["full_point", "sub_point", "sub_pct_lo", "sub_pct_hi", "rc_point", "rc_lo", "rc_hi"])
        y = yardsticks(d.sub_point, d.sub_pct_lo, d.sub_pct_hi, d.rc_point, d.rc_lo, d.rc_hi, rng,
                       full_pt=d.full_point)
        rows.append({"cell_id": int(c.cell_id), "scenario": c.scenario, "K": int(c.K), "pi_x": c.pi_x,
                     "size": c["size"], **y})
    g = pd.DataFrame(rows)
    g.to_csv(os.path.join(OUT, "r1_rc_decomposition_grid.csv"), index=False)
    cols = ["w0_sub", "w0_rc", "w_single", "w_uncal", "mixing_factor", "point_sd_ratio", "sd_full_over_sub",
            "overwidth_vs_own_point", "overwidth_vs_single_draw", "type1_rc"]
    byk = g.groupby("K")[cols].agg(["mean", "min", "max"])
    byk.columns = [f"{a}_{b}" for a, b in byk.columns]
    byk = byk.reset_index()
    byk.to_csv(os.path.join(OUT, "r1_rc_decomposition_by_K.csv"), index=False)
    return g, byk


def claims(rng):
    df = pd.read_csv(os.path.join(OUT, "claims_matched_null_sim.csv.gz"))
    rows = []
    for (claim, variant), g in df.groupby(["claim", "variant"]):
        if "rc_lo" not in g or g.rc_lo.isna().all() or g.sub_lo.isna().all():
            continue
        y = yardsticks(g.sub_point, g.sub_lo, g.sub_hi, g.rc_point, g.rc_lo, g.rc_hi, rng)
        rows.append({"claim": claim, "variant": variant, "source": "claims_matched_null_sim.csv.gz (seed 20260931)", **y})
    c = pd.DataFrame(rows)
    c.to_csv(os.path.join(OUT, "r1_rc_yardsticks_claims.csv"), index=False)
    return c


def main():
    rng = np.random.default_rng(SEED)
    g, byk = grid(rng)
    c = claims(rng)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 60)
    show = ["K", "w0_rc_mean", "w_single_mean", "w_single_min", "w_single_max", "w_uncal_mean", "mixing_factor_min",
            "mixing_factor_max", "point_sd_ratio_min", "point_sd_ratio_max", "sd_full_over_sub_mean",
            "overwidth_vs_own_point_mean", "overwidth_vs_single_draw_min", "overwidth_vs_single_draw_max"]
    print(byk[show].round(3).to_string())
    print(c[["claim", "variant", "n_datasets", "w0_sub", "w0_rc", "w_single", "w_single_lo", "w_single_hi", "w_uncal",
             "mixing_factor", "point_sd_ratio"]].round(3).to_string())


if __name__ == "__main__":
    main()
