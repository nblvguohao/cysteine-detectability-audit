# -*- coding: utf-8 -*-
"""E_bootstrap_coverage -- POST HOC revision analysis (2026-09-30). Not registered.

Fidelity checks run before the simulation grid:
  V1  pipeline log_odds_ratio (verbatim) == tool log2_or on random tables
  V2  vectorised MH statistic == pipeline make_mh_log_odds_ratio (verbatim), incl. bootstrap weights
  V3  closed-form sandwich (CR0, Mancl-DeRouen) == statsmodels GEE (independence, logit),
      cov_type 'robust' and 'bias_reduced' (MLE, no Haldane term)
  V4  vectorised random control vs the pipeline's random_control_effect (verbatim) on the same data
  V5  tool percentile bootstrap vs the pipeline's bootstrap_interval (verbatim) on the same data
Writes results/E_bootstrap_coverage/validation.json
"""
from __future__ import annotations

import json
import os
import sys
import time

sys.dont_write_bytecode = True
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import simlib as S

OUT = r"C:/Users/admin/Desktop/小论文/_cys_repo_work/public/revision_2026-09-30/results/E_bootstrap_coverage"
SEED = 20260930


def v1():
    P = S.pipeline()
    rng = np.random.default_rng(SEED)
    worst = 0.0
    for _ in range(2000):
        t = rng.integers(0, 60, size=4)
        a = P["log_odds_ratio"](t[0], t[2], t[1], t[3])  # n11 attr&pos, n10 attr&neg, n01 noattr&pos, n00
        b = float(S.tool_stats.log2_or(t[0], t[1], t[2], t[3]))
        worst = max(worst, abs(a - b))
    return {"max_abs_diff": worst, "pass": worst < 1e-12}


def v2():
    P = S.pipeline()
    rng = np.random.default_rng(SEED + 1)
    par = S.calibrate(0.15, 0.35, 1.0, 1.0, 0.0)
    worst, n = 0.0, 0
    for rep in range(20):
        d = S.generate(rng, 164, "nb", par, 1.0, 1.0, 0.0)
        y, x, s, cid, K = d["y"], d["x"], d["s"], d["cid"], d["K"]
        strata = P["quantile_strata"](s)
        f = P["make_mh_log_odds_ratio"](y.astype(int), x.astype(int), strata)
        cell_code = np.where(x & y, 0, np.where(x & ~y, 1, np.where(~x & y, 2, 3)))
        cell = np.zeros((K, 20))
        np.add.at(cell, (cid, strata * 4 + cell_code), 1)
        for w in [np.ones(K)] + [rng.multinomial(K, np.full(K, 1 / K)).astype(float) for _ in range(10)]:
            a = f(w[cid])
            b = float(S.mh_stat((w @ cell)[None, :])[0])
            if np.isfinite(a) or np.isfinite(b):
                worst = max(worst, abs(a - b))
            n += 1
    return {"comparisons": n, "max_abs_diff": worst, "pass": worst < 1e-9}


def v3():
    import statsmodels.api as sm
    rng = np.random.default_rng(SEED + 2)
    rows = []
    for K, pix in ((44, 0.15), (164, 0.35), (44, 0.35)):
        par = S.calibrate(0.15, pix, 1.0, 1.0, 0.0)
        for rep in range(5):
            d = S.generate(rng, K, "nb", par, 1.0, 1.0, 0.0)
            counts = S.tool_counts(d["cid"], d["y"], d["x"], K)
            if (counts.sum(axis=0) == 0).any():
                continue
            mine = S.wald_sandwich(counts, haldane=False)
            X = sm.add_constant(d["x"].astype(float))
            res = {}
            for cov in ("robust", "bias_reduced"):
                m = sm.GEE(d["y"].astype(float), X, groups=d["cid"], family=sm.families.Binomial(),
                           cov_struct=sm.cov_struct.Independence())
                r = m.fit(cov_type=cov)
                res[cov] = (float(r.params[1] / S.LN2), float(r.bse[1] / S.LN2))
            rows.append({"K": K, "pi_x": pix, "beta_mine": mine["wald_point"], "beta_gee": res["robust"][0],
                         "se_cr0_mine": mine["se_cr0"], "se_robust_gee": res["robust"][1],
                         "se_md_mine": mine["se_md"], "se_bias_reduced_gee": res["bias_reduced"][1]})
    worst_b = max(abs(r["beta_mine"] - r["beta_gee"]) for r in rows)
    worst_cr0 = max(abs(r["se_cr0_mine"] / r["se_robust_gee"] - 1) for r in rows)
    worst_md = max(abs(r["se_md_mine"] / r["se_bias_reduced_gee"] - 1) for r in rows)
    return {"n": len(rows), "max_abs_diff_beta": worst_b, "max_rel_diff_se_cr0": worst_cr0,
            "max_rel_diff_se_md": worst_md, "rows": rows,
            "pass": worst_b < 1e-6 and worst_cr0 < 1e-4 and worst_md < 1e-3}


def v4(n_data=200, K=44):
    P = S.pipeline()
    rng = np.random.default_rng(SEED + 3)
    par = S.calibrate(0.15, 0.35, 1.0, 1.0, 0.0)
    flags_v, flags_p, w_v, w_p, t0 = [], [], [], [], time.time()
    for rep in range(n_data):
        d = S.generate(rng, K, "nb", par, 1.0, 1.0, 0.0)
        y, x, cid = d["y"], d["x"], d["cid"]
        seed = 7_000_000 + rep * 50
        rv = S.random_control_vectorised(cid, y, x, 5000, seed)
        pos = np.flatnonzero(y)
        point, iv, n_neg, _ = P["random_control_effect"](y.astype(int), x.astype(int), cid, pos, len(pos),
                                                        replicates=5000, seed=seed + 25, draws=20)
        flags_v.append(not (rv["rc_lo"] <= 0 <= rv["rc_hi"]))
        flags_p.append(not (iv[0] <= 0 <= iv[1]))
        w_v.append(rv["rc_hi"] - rv["rc_lo"])
        w_p.append(iv[1] - iv[0])
    return {"datasets": n_data, "K": K, "flag_rate_vectorised": float(np.mean(flags_v)),
            "flag_rate_pipeline_verbatim": float(np.mean(flags_p)),
            "mean_width_vectorised": float(np.mean(w_v)), "mean_width_pipeline_verbatim": float(np.mean(w_p)),
            "seconds": round(time.time() - t0, 1),
            "note": "same data, independent random streams; agreement is distributional"}


def v5(n_data=300, K=164):
    P = S.pipeline()
    rng = np.random.default_rng(SEED + 4)
    par = S.calibrate(0.15, 0.35, 1.0, 1.0, 0.0)
    agree, dl, dh, ft, fp = 0, [], [], [], []
    for rep in range(n_data):
        d = S.generate(rng, K, "nb", par, 1.0, 1.0, 0.0)
        y, x, cid = d["y"].astype(int), d["x"].astype(int), d["cid"]
        counts = S.tool_counts(d["cid"], d["y"], d["x"], K)
        _, (lo_t, hi_t), _, _ = S.tool_stats.boot_log2_or(counts, 5000, 9_000_000 + rep, 0.95)
        table = P["make_weighted_table"](y, x)
        (lo_p, hi_p), _ = P["bootstrap_interval"](lambda w: P["log_odds_ratio"](*table(w)), cid,
                                                  replicates=5000, seed=9_500_000 + rep)
        f_t = not (lo_t <= 0 <= hi_t)
        f_p = not (lo_p <= 0 <= hi_p)
        agree += int(f_t == f_p)
        ft.append(f_t)
        fp.append(f_p)
        dl.append(lo_t - lo_p)
        dh.append(hi_t - hi_p)
    return {"datasets": n_data, "K": K, "flag_rate_tool": float(np.mean(ft)), "flag_rate_pipeline": float(np.mean(fp)),
            "flag_agreement": agree / n_data, "mean_abs_diff_lo": float(np.mean(np.abs(dl))),
            "mean_abs_diff_hi": float(np.mean(np.abs(dh))), "mean_diff_lo": float(np.mean(dl)),
            "mean_diff_hi": float(np.mean(dh))}


def main():
    out = {"label": "POST HOC revision analysis E_bootstrap_coverage: fidelity checks (not registered)",
           "pipeline_main_sha256": S.sha256(S.PIPE_MAIN), "pipeline_2b_sha256": S.sha256(S.PIPE_2B),
           "tool_stats_sha256": S.sha256(S.TOOL_STATS)}
    for name, fn in (("V1_log2or", v1), ("V2_mh", v2), ("V3_sandwich_vs_statsmodels_gee", v3),
                     ("V4_random_control_K44", lambda: v4(400, 44)),
                     ("V4_random_control_K164", lambda: v4(200, 164)),
                     ("V5_tool_vs_pipeline_bootstrap", v5)):
        t0 = time.time()
        out[name] = fn()
        print(name, round(time.time() - t0, 1), "s", json.dumps({k: v for k, v in out[name].items() if k != "rows"}))
    with open(os.path.join(OUT, "validation.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
