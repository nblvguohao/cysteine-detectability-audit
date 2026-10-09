# -*- coding: utf-8 -*-
"""E_bootstrap_coverage -- POST HOC revision analysis (2026-09-30). Not registered, not pre-specified.

Simulation grid (tasks 1-2). Every dataset is analysed in three designs that mirror the three
intervals a claim verdict reads:
  full  the whole cohort (baseline interval): percentile [paper], BCa, cluster-robust Wald
        (CR0-normal, CR1-t(G-1), Mancl-DeRouen-t(G-1)); plus the Mantel-Haenszel log2 OR over
        quintiles of a propensity-like covariate (secondary caliber), percentile
  sub   all positives + an equal number of random negatives, clusters re-indexed to those present
        (the size and structure of a matched set / one random-control draw): same five methods
  rc    the pipeline's same-size random control (20 draws, each replicate picks one draw),
        percentile
Grid: K (proteins) {44,164,776,1475} x attribute prevalence {0.15,0.35} x cluster-size family
{nb, lognormal} x scenario {null_iid, null_re, null_conf15, null_conf50, alt05, alt10}.
Bootstrap replicates 5,000 (the paper's number) everywhere.

usage: python run_grid.py [--datasets N] [--only-scenarios a,b] [--only-K 44,164] [--workers 4]
"""
from __future__ import annotations

import argparse
import json
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

import simlib as S

OUT = r"/path/to/local/_cys_repo_work/public/revision_2026-09-30/results/E_bootstrap_coverage"
CELLS_DIR = os.path.join(OUT, "sim_cells")
MASTER = 20260930
REPS = 5000
PI_Y = 0.15
SIGMA = 1.0  # label and attribute random-intercept SD (logit scale), anchored (anchor_parameters.json)

SCENARIOS = {
    "null_iid":    dict(su=0.0, sv=0.0, rho=0.0, target=None),
    "null_re":     dict(su=SIGMA, sv=SIGMA, rho=0.0, target=None),
    "null_conf15": dict(su=SIGMA, sv=SIGMA, rho=0.15, target=None),
    "null_conf50": dict(su=SIGMA, sv=SIGMA, rho=0.50, target=None),
    "alt05":       dict(su=SIGMA, sv=SIGMA, rho=0.0, target=0.5),
    "alt10":       dict(su=SIGMA, sv=SIGMA, rho=0.0, target=1.0),
}
KS = [44, 164, 776, 1475]
PIXS = [0.15, 0.35]
SIZES = ["nb", "lognormal"]


def cell_table():
    rows, cid = [], 0
    for scen in SCENARIOS:
        for K in KS:
            for pix in PIXS:
                for size in SIZES:
                    sc = SCENARIOS[scen]
                    par = S.calibrate(PI_Y, pix, sc["su"], sc["sv"], sc["rho"], target_log2or=sc["target"])
                    rows.append(dict(cell_id=cid, scenario=scen, K=K, pi_x=pix, size=size, su=sc["su"],
                                     sv=sc["sv"], rho=sc["rho"], within_log2or=par["gamma_log2"],
                                     theta_pooled_log2=par["theta_pooled_log2"], alpha0=par["alpha0"],
                                     beta0=par["beta0"], gamma=par["gamma"], positive_rate=par["positive_rate"]))
                    cid += 1
    return pd.DataFrame(rows)


def one_dataset(cell, rep):
    ss = np.random.SeedSequence([MASTER, int(cell["cell_id"]), int(rep)])
    rng = np.random.default_rng(ss)
    par = {"alpha0": cell["alpha0"], "beta0": cell["beta0"], "gamma": cell["gamma"]}
    d = S.generate(rng, int(cell["K"]), cell["size"], par, cell["su"], cell["sv"], cell["rho"])
    y, x, s, cid, K = d["y"], d["x"], d["s"], d["cid"], d["K"]
    base = 1_000_000_000 + int(cell["cell_id"]) * 10_000_000 + int(rep) * 100
    row = {"cell_id": int(cell["cell_id"]), "rep": int(rep)}
    full = S.analyse_table(S.tool_counts(cid, y, x, K), REPS, base + 1)
    row.update({f"full_{k}": v for k, v in full.items()})
    mh = S.mh_percentile(cid, K, y, x, s, REPS, base + 2)
    row.update({f"full_{k}": v for k, v in mh.items()})
    idx = S.subsample_1to1(rng, y)
    inv, Ks = S.reindex(cid[idx])
    sub = S.analyse_table(S.tool_counts(inv, y[idx], x[idx], Ks), REPS, base + 3)
    row.update({f"sub_{k}": v for k, v in sub.items()})
    rc = S.random_control_vectorised(cid, y, x, REPS, base + 50)
    row.update(rc)
    return row


def work(args):
    cell, start, stop = args
    t0 = time.time()
    rows = [one_dataset(cell, r) for r in range(start, stop)]
    return int(cell["cell_id"]), start, rows, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", type=int, default=0, help="overrides --datasets-by-K when > 0")
    ap.add_argument("--datasets-by-K", default="44:4000,164:3000,776:1000,1475:1000")
    ap.add_argument("--only-scenarios", default="")
    ap.add_argument("--only-K", default="")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--chunk", type=int, default=50)
    ap.add_argument("--timing", action="store_true")
    a = ap.parse_args()
    os.makedirs(CELLS_DIR, exist_ok=True)
    cells = cell_table()
    cells.to_csv(os.path.join(OUT, "sim_cell_parameters.csv"), index=False)
    sel = cells
    if a.only_scenarios:
        sel = sel[sel.scenario.isin(a.only_scenarios.split(","))]
    if a.only_K:
        sel = sel[sel.K.isin([int(k) for k in a.only_K.split(",")])]
    if a.timing:
        for K in KS:
            c = cells[(cells.K == K) & (cells.scenario == "null_re")].iloc[0].to_dict()
            t0 = time.time()
            for r in range(5):
                one_dataset(c, 10_000 + r)
            print("K", K, "sec/dataset", round((time.time() - t0) / 5, 3), flush=True)
        return
    by_k = {int(k): int(n) for k, n in (kv.split(":") for kv in a.datasets_by_K.split(","))}

    def n_for(K):
        return a.datasets if a.datasets > 0 else by_k[int(K)]

    tasks = []
    for _, c in sel.iterrows():
        path = os.path.join(CELLS_DIR, f"cell_{int(c.cell_id):03d}.csv.gz")
        if os.path.exists(path):
            continue
        n = n_for(c.K)
        for st in range(0, n, a.chunk):
            tasks.append((c.to_dict(), st, min(n, st + a.chunk)))
    # small cells first (early results for checking), then larger ones
    tasks.sort(key=lambda t: (t[0]["K"], t[0]["cell_id"], t[1]))
    print("tasks", len(tasks), "cells", sel.shape[0], flush=True)
    buf = {}
    need = {int(c.cell_id): n_for(c.K) for _, c in sel.iterrows()}
    t0 = time.time()
    from multiprocessing import Pool
    with Pool(a.workers) as pool:
        for cid, st, rows, dt in pool.imap_unordered(work, tasks):
            buf.setdefault(cid, []).extend(rows)
            if len(buf[cid]) >= need[cid]:
                df = pd.DataFrame(buf.pop(cid)).sort_values("rep")
                df.to_csv(os.path.join(CELLS_DIR, f"cell_{cid:03d}.csv.gz"), index=False)
                print(f"cell {cid} done ({round(time.time() - t0)} s elapsed)", flush=True)
    print("all done", round(time.time() - t0), "s")


if __name__ == "__main__":
    main()
