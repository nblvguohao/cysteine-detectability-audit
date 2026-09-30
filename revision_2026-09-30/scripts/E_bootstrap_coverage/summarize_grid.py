# -*- coding: utf-8 -*-
"""E_bootstrap_coverage -- POST HOC revision analysis (2026-09-30). Not registered.

Summarises the simulation grid (results/E_bootstrap_coverage/sim_cells/*.csv.gz):
  flag0      P(interval excludes 0)           -> type-I error for the null scenarios
                                                  (for null_conf*: rate of 'excludes 0' when the
                                                  WITHIN-protein association is exactly 0)
  cover      P(interval contains theta), theta = population pooled (marginal) log2 OR
  miss_lo / miss_hi  P(theta below / above the interval)
  w0         95th percentile over datasets of |point - theta| / (half-width toward theta):
             the factor by which the interval's half-widths must be multiplied for 95% coverage
             (1 = calibrated; > 1 anti-conservative; < 1 conservative), with a 95% order-statistic CI
  width      mean interval width (log2 units)
Monte Carlo SE for proportions: sqrt(p(1-p)/n).
Writes sim_summary_long.csv.
"""
from __future__ import annotations

import glob
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
from scipy import stats as sps

import simlib as S

OUT = r"C:/Users/admin/Desktop/小论文/_cys_repo_work/public/revision_2026-09-30/results/E_bootstrap_coverage"

DESIGN_METHODS = [
    ("full", "pct", "full_pct_lo", "full_pct_hi", "full_point"),
    ("full", "bca", "full_bca_lo", "full_bca_hi", "full_point"),
    ("full", "cr0", "full_cr0_lo", "full_cr0_hi", "full_wald_point"),
    ("full", "cr1t", "full_cr1t_lo", "full_cr1t_hi", "full_wald_point"),
    ("full", "mdt", "full_mdt_lo", "full_mdt_hi", "full_wald_point"),
    ("full_mh", "pct", "full_mh_lo", "full_mh_hi", "full_mh_point"),
    ("sub", "pct", "sub_pct_lo", "sub_pct_hi", "sub_point"),
    ("sub", "bca", "sub_bca_lo", "sub_bca_hi", "sub_point"),
    ("sub", "cr0", "sub_cr0_lo", "sub_cr0_hi", "sub_wald_point"),
    ("sub", "cr1t", "sub_cr1t_lo", "sub_cr1t_hi", "sub_wald_point"),
    ("sub", "mdt", "sub_mdt_lo", "sub_mdt_hi", "sub_wald_point"),
    ("rc", "pct", "rc_lo", "rc_hi", "rc_point"),
]


def quantile_ci(x, q=0.95, conf=0.95):
    """Distribution-free CI for a quantile from order statistics."""
    x = np.sort(np.asarray(x, float))
    n = x.size
    if n == 0:
        return np.nan, np.nan, np.nan
    est = float(np.quantile(x, q))
    lo_k = int(sps.binom.ppf((1 - conf) / 2, n, q))
    hi_k = int(sps.binom.ppf(1 - (1 - conf) / 2, n, q))
    lo_k = max(0, min(n - 1, lo_k - 1))
    hi_k = max(0, min(n - 1, hi_k))
    return est, float(x[lo_k]), float(x[hi_k])


def summarise_cell(df, cell):
    theta = float(cell["theta_pooled_log2"])
    out = []
    n = len(df)
    for design, method, lo_c, hi_c, pt_c in DESIGN_METHODS:
        if lo_c not in df:
            continue
        lo, hi, pt = df[lo_c].to_numpy(float), df[hi_c].to_numpy(float), df[pt_c].to_numpy(float)
        ok = np.isfinite(lo) & np.isfinite(hi)
        lo, hi, pt = lo[ok], hi[ok], pt[ok]
        m = int(ok.sum())
        flag0 = float(np.mean((lo > 0) | (hi < 0)))
        cover = float(np.mean((lo <= theta) & (theta <= hi)))
        within = float(cell["within_log2or"])
        cover_within = float(np.mean((lo <= within) & (within <= hi)))
        miss_lo = float(np.mean(theta < lo))
        miss_hi = float(np.mean(theta > hi))
        r = S.ratio_to_target(pt, lo, hi, theta)
        w0, w0_lo, w0_hi = quantile_ci(r)
        r0 = S.ratio_to_target(pt, lo, hi, 0.0)
        w00, w00_lo, w00_hi = quantile_ci(r0)
        out.append(dict(cell_id=int(cell["cell_id"]), scenario=cell["scenario"], K=int(cell["K"]),
                        pi_x=float(cell["pi_x"]), size=cell["size"], rho=float(cell["rho"]),
                        theta=theta, within_log2or=float(cell["within_log2or"]), design=design, method=method,
                        n_datasets=n, n_finite=m, flag0=flag0, flag0_mcse=float(np.sqrt(flag0 * (1 - flag0) / m)),
                        cover=cover, cover_mcse=float(np.sqrt(cover * (1 - cover) / m)),
                        cover_within=cover_within,
                        miss_lo=miss_lo, miss_hi=miss_hi, w0=w0, w0_ci_lo=w0_lo, w0_ci_hi=w0_hi,
                        w0_vs_zero=w00, w0_vs_zero_ci_lo=w00_lo, w0_vs_zero_ci_hi=w00_hi,
                        mean_width=float(np.mean(hi - lo)), mean_point=float(np.mean(pt)),
                        bias=float(np.mean(pt) - theta),
                        mean_K_eff=float(df[f"{design.split('_')[0]}_K"].mean()) if f"{design.split('_')[0]}_K" in df else np.nan,
                        mean_n11=float(df[f"{design.split('_')[0]}_n11"].mean()) if f"{design.split('_')[0]}_n11" in df else np.nan))
    return out


def main():
    cells = pd.read_csv(os.path.join(OUT, "sim_cell_parameters.csv")).set_index("cell_id")
    rows = []
    for f in sorted(glob.glob(os.path.join(OUT, "sim_cells", "cell_*.csv.gz"))):
        cid = int(os.path.basename(f)[5:8])
        df = pd.read_csv(f)
        cell = cells.loc[cid].to_dict()
        cell["cell_id"] = cid
        rows.extend(summarise_cell(df, cell))
    long = pd.DataFrame(rows)
    long.to_csv(os.path.join(OUT, "sim_summary_long.csv"), index=False)
    print("cells summarised:", long.cell_id.nunique())
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 500)
    show = long[(long.method == "pct")][["scenario", "K", "pi_x", "size", "design", "flag0", "cover", "miss_lo",
                                          "miss_hi", "w0", "w0_ci_lo", "w0_ci_hi", "mean_width"]]
    print(show.round(4).to_string())


if __name__ == "__main__":
    main()
