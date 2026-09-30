# -*- coding: utf-8 -*-
"""E_bootstrap_coverage -- REVISION AFTER VERIFICATION (round 1), 2026-09-30.
POST HOC revision analysis; not registered, not pre-specified.

Verifier problem 1: SFE-008 (transfer re-test; stored random control 0.4111 [-0.0698, 0.8696],
w* = 0.855) flips undecidable -> vanishes under the NOMINAL yardstick for the random-control interval,
but under the purpose-matched SINGLE-DRAW yardstick it sits on the boundary (this item's first run:
w_single 0.871 [0.831, 0.902]; the verifier's independent run: 0.842). This script pins the single-draw
yardstick down at SFE-008's cohort shape with more datasets and under three cluster models:

  re_nb2     protein random intercepts SD 1 (label and attribute, independent), cluster size
             1 + NB(size 2), mean 1.786 (the first run's model; new seeds)                8,000 datasets
  re_nb077   as re_nb2 but NB size 0.77, which matches BOTH the mean cluster size (1.786) and the
             stored singleton fraction (2,128 of 3,662 clusters = 58.1%)                  4,000 datasets
  iid_nb077  as re_nb077 without random effects                                            4,000 datasets

Cohort shape (claims_intervals.csv): 3,662 clusters, 6,542 sites, 718 positives, attribute on 5,284 sites.
No attribute-label association. Designs: sub = all positives + as many random negatives (one same-size
draw; percentile interval by the released tool, cys_audit.stats.boot_log2_or, 5,000 replicates); rc =
the pipeline's same-size random control (simlib.random_control_vectorised, validated against the
verbatim pipeline function in validate.py V4).
Seeds: data SeedSequence([20260932, variant, dataset]); bootstrap seeds fixed functions of (variant, dataset).
Writes r1_sfe008_sim.csv.gz (per dataset) and r1_sfe008_yardsticks.json.
"""
from __future__ import annotations

import json
import math
import os
import sys
import time

sys.dont_write_bytecode = True
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
from scipy import special

import simlib as S
from r1_rc_yardsticks import yardsticks

OUT = r"C:/Users/admin/Desktop/小论文/_cys_repo_work/public/revision_2026-09-30/results/E_bootstrap_coverage"
MASTER = 20260932
REPS = 5000
K, N_SITES, N_POS, N_ATTR = 3662, 6542, 718, 5284
STORED = {"point": 0.4111, "lo": -0.0698, "hi": 0.8696}
VARIANTS = [  # name, sigma, NB size, datasets
    ("re_nb2", 1.0, 2.0, 8000),
    ("re_nb077", 1.0, 0.77, 4000),
    ("iid_nb077", 0.0, 0.77, 4000),
]


def sizes(rng, nb_size, mean_size=N_SITES / K):
    extra = mean_size - 1.0
    return 1 + rng.negative_binomial(nb_size, nb_size / (nb_size + extra), size=K)


def one(args):
    vidx, name, sigma, nb_size, rep, par = args
    rng = np.random.default_rng(np.random.SeedSequence([MASTER, vidx, rep]))
    m = sizes(rng, nb_size)
    z1, z2 = rng.standard_normal(K), rng.standard_normal(K)
    cid = np.repeat(np.arange(K), m)
    x = rng.random(cid.size) < special.expit(par["beta0"] + sigma * z2[cid])
    y = rng.random(cid.size) < special.expit(par["alpha0"] + sigma * z1[cid])
    base = 3_000_000_000 + vidx * 10_000_000 + rep * 100
    full_pt = float(S.log2_or_rows(S.tool_counts(cid, y, x, K).sum(axis=0)))
    idx = S.subsample_1to1(rng, y)
    inv, Ks = S.reindex(cid[idx])
    pt, (lo, hi), _, _ = S.tool_stats.boot_log2_or(S.tool_counts(inv, y[idx], x[idx], Ks), REPS, base + 3, 0.95)
    rc = S.random_control_vectorised(cid, y, x, REPS, base + 50)
    return {"variant": name, "rep": rep, "n": int(cid.size), "n_pos": int(y.sum()), "n_attr": int(x.sum()),
            "singleton_fraction": float(np.mean(m == 1)), "full_point": full_pt,
            "sub_point": pt, "sub_lo": lo, "sub_hi": hi, "sub_K": Ks,
            "rc_point": rc["rc_point"], "rc_lo": rc["rc_lo"], "rc_hi": rc["rc_hi"]}


def work(chunk):
    return [one(a) for a in chunk]


def main():
    t0 = time.time()
    tasks = []
    for vidx, (name, sigma, nb_size, n) in enumerate(VARIANTS):
        par = S.calibrate(N_POS / N_SITES, N_ATTR / N_SITES, sigma, sigma, 0.0)
        args = [(vidx, name, sigma, nb_size, r, {"alpha0": par["alpha0"], "beta0": par["beta0"]}) for r in range(n)]
        tasks += [args[i:i + 100] for i in range(0, n, 100)]
    rows = []
    from multiprocessing import Pool
    with Pool(4) as pool:
        for res in pool.imap_unordered(work, tasks):
            rows.extend(res)
    df = pd.DataFrame(rows).sort_values(["variant", "rep"]).reset_index(drop=True)
    df.to_csv(os.path.join(OUT, "r1_sfe008_sim.csv.gz"), index=False)
    rng = np.random.default_rng(MASTER)
    p, lo = STORED["point"], STORED["lo"]
    res = {"label": "POST HOC revision analysis E_bootstrap_coverage, revision after verification (not registered)",
           "claim": "SFE-008 (transfer re-test)", "stored_random_control": STORED,
           "wstar": p / (p - lo), "seconds": round(time.time() - t0), "variants": {}}
    for name, g in df.groupby("variant"):
        y = yardsticks(g.sub_point, g.sub_lo, g.sub_hi, g.rc_point, g.rc_lo, g.rc_hi, rng, full_pt=g.full_point)
        y["mean_singleton_fraction"] = float(g.singleton_fraction.mean())
        for key in ("w0_rc", "w_single", "w_single_lo", "w_single_hi", "w_uncal"):
            y[f"SFE008_rc_lower_bound_at_{key}"] = p - y[key] * (p - lo)
        y["verdict_at_w_single"] = "vanishes" if p - y["w_single"] * (p - lo) > 0 else "undecidable"
        res["variants"][name] = y
    with open(os.path.join(OUT, "r1_sfe008_yardsticks.json"), "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
