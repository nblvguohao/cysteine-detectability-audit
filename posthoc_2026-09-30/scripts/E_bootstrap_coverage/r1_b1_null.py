# -*- coding: utf-8 -*-
"""E_bootstrap_coverage -- REVISION AFTER VERIFICATION (round 1), 2026-09-30.
POST HOC revision analysis; not registered, not pre-specified.

Verifier problem 5: the first report applied 'the same interval construction' to both extended
strength-zero rates of the registered planted-artifact benchmark. B2 (abundance) is a rank AUC over
3,000 single-row proteins and is NOT addressed by this item's simulation. B1 (cleavage) IS a Haldane
log2 OR with a protein-clustered percentile bootstrap (5,000 replicates), on 300 synthetic proteins.
This script simulates B1's own null directly, with B1's own generator, instead of interpolating the grid:

  * cut_sites, detectable_positions and seed_for are compiled VERBATIM from
    repo/scripts/run_phase2_synthetic_benchmark_2026-09-22.py (AST extraction);
  * the strength-0 path of run_b1 (protein lengths, cut sites, detectable cysteines, near-cut flag,
    positives drawn uniformly with replacement then de-duplicated, statistic with +0.5 per cell) and
    cluster_boot_ci (seed + 1, rng.integers draws of clusters, 2.5th/97.5th percentiles) are RE-TYPED
    here in vectorised form, and checked against the stored registered replicates: the 50 strength-0
    replicates in phase2_synthetic_benchmark_replicates_2026-09-22.csv (point, interval, n) and the
    registered extension count (21 of 300, results/phase2_calibration_precision_2026-09-22.csv), whose
    replicates come from the same seed sequence (seed_for(20260922, 0, 0, rep), rep = 0..299);
  * then 10,000 FRESH strength-0 datasets (seeds SeedSequence([20260932, 7, rep]) for data and
    SeedSequence([20260932, 8, rep]) for the bootstrap) give B1's type-I error and calibration
    multiplier w0 at its own shape.
Writes r1_b1_null.csv.gz (per dataset, registered and fresh) and r1_b1_null.json.
"""
from __future__ import annotations

import ast
import json
import os
import sys
import time

sys.dont_write_bytecode = True
os.environ["OPENBLAS_NUM_THREADS"] = "1"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
from scipy import stats as sps

import simlib as S
from summarize_grid import quantile_ci

OUT = r"/path/to/local/_cys_repo_work/public/revision_2026-09-30/results/E_bootstrap_coverage"
BENCH_PY = r"/path/to/local/_cys_repo_work/repo/scripts/run_phase2_synthetic_benchmark_2026-09-22.py"
PROTOCOL = r"/path/to/local/_cys_repo_work/repo/protocols/phase2_synthetic_benchmark_preregistration_2026-09-22.json"
STORED_REPL = r"/path/to/local/_cys_repo_work/public/revision_2026-09-30/inputs/repo_results/phase2_synthetic_benchmark_replicates_2026-09-22.csv"
STORED_EXT = r"/path/to/local/_cys_repo_work/repo/results/phase2_calibration_precision_2026-09-22.csv"
MASTER = 20260932
N_FRESH = 10000
N_REGISTERED = 300


def _verbatim():
    src = open(BENCH_PY, encoding="utf-8").read()
    keep = [n for n in ast.parse(src).body if isinstance(n, ast.FunctionDef)
            and n.name in ("cut_sites", "detectable_positions", "seed_for")]
    ns = {"np": np, "__name__": "b1_verbatim"}
    exec(compile(ast.Module(body=keep, type_ignores=[]), BENCH_PY, "exec"), ns)
    return ns


NS = _verbatim()
P = json.load(open(PROTOCOL, encoding="utf-8"))
SPEC = P["benchmarks"]["B1_cleavage"]
REPS = int(P["bootstrap_reps"])


def generate_s0(rng):
    """Strength-0 path of run_b1, re-typed (same calls to rng, in the same order)."""
    spec = SPEC
    plens = rng.integers(spec["protein_length_range"][0], spec["protein_length_range"][1] + 1, size=spec["n_proteins"])
    prot_ids, dists = [], []
    for pi, L in enumerate(plens):
        cuts = NS["cut_sites"](L, spec["cleavage_site_rate_per_residue"], rng)
        det, dist = NS["detectable_positions"](L, cuts, spec["detectable_peptide_length_range"][0],
                                               spec["detectable_peptide_length_range"][1], spec["max_missed_cleavages"])
        is_cys = rng.random(L) < spec["cys_rate_per_residue"]
        cys_det = np.flatnonzero(is_cys & det)
        if len(cys_det) == 0:
            continue
        prot_ids.extend([pi] * len(cys_det))
        dists.extend(dist[cys_det].tolist())
    prot_ids = np.array(prot_ids)
    near = np.array(dists) <= spec["near_cut_window_residues"]
    n_bg = len(prot_ids)
    n_pos = max(1, int(round(spec["true_positive_fraction_of_detectable"] * n_bg)))
    _use_near_pool = rng.random(n_pos) < 0.0  # consumed exactly as in run_b1 (all False at s = 0)
    pos_idx = np.empty(n_pos, dtype=np.int64)
    for k in range(n_pos):
        pos_idx[k] = rng.integers(0, n_bg)
    pos_idx = np.unique(pos_idx)
    is_pos = np.zeros(n_bg, dtype=bool)
    is_pos[pos_idx] = True
    return prot_ids, near, is_pos


def stat_from_counts(a, b, c, d):
    a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    return np.log2((a / b) / (c / d))


def cluster_boot(prot_ids, near, is_pos, rng):
    """cluster_boot_ci of the benchmark script, vectorised: draws = rng.integers(0, n_clusters,
    (reps, n_clusters)) over sorted unique clusters; each replicate's rows are the union (with
    repetition) of the chosen clusters' rows, i.e. multiplicity-weighted per-cluster 2x2 tables."""
    uniq, inv = np.unique(prot_ids, return_inverse=True)
    G = len(uniq)
    code = np.where(is_pos, np.where(near, 0, 1), np.where(near, 2, 3))
    T = np.zeros((G, 4))
    np.add.at(T, (inv, code), 1.0)
    draws = rng.integers(0, G, size=(REPS, G))
    W = np.zeros((REPS, G))
    np.add.at(W, (np.repeat(np.arange(REPS), G), draws.ravel()), 1.0)
    tab = W @ T
    vals = stat_from_counts(tab[:, 0], tab[:, 1], tab[:, 2], tab[:, 3])
    return float(np.quantile(vals, 0.025)), float(np.quantile(vals, 0.975)), G


def one(args):
    kind, rep = args
    if kind == "registered":
        seed = NS["seed_for"](P["master_seed"], 0, 0, rep)
        rng = np.random.default_rng(seed)
        brng = np.random.default_rng(seed + 1)
    else:
        rng = np.random.default_rng(np.random.SeedSequence([MASTER, 7, rep]))
        brng = np.random.default_rng(np.random.SeedSequence([MASTER, 8, rep]))
    prot_ids, near, is_pos = generate_s0(rng)
    a = float(np.sum(is_pos & near)); b = float(np.sum(is_pos & ~near))
    c = float(np.sum(~is_pos & near)); d = float(np.sum(~is_pos & ~near))
    point = float(stat_from_counts(a, b, c, d))
    lo, hi, G = cluster_boot(prot_ids, near, is_pos, brng)
    return {"kind": kind, "rep": rep, "statistic": point, "ci_lo": lo, "ci_hi": hi, "n_positions": int(len(prot_ids)),
            "n_positive": int(is_pos.sum()), "n_clusters": G, "near_prevalence": float(near.mean()), "n11": int(a)}


def work(chunk):
    return [one(a) for a in chunk]


def summarise(g):
    pt, lo, hi = g.statistic.to_numpy(), g.ci_lo.to_numpy(), g.ci_hi.to_numpy()
    flag = (lo > 0) | (hi < 0)
    k, n = int(flag.sum()), len(g)
    w0, w0l, w0h = quantile_ci(S.ratio_to_zero(pt, lo, hi))
    ci = sps.binomtest(k, n).proportion_ci(0.95, method="wilson")
    return {"n": n, "flagged": k, "rate": k / n, "mcse": float(np.sqrt(k / n * (1 - k / n) / n)),
            "wilson_lo": float(ci.low), "wilson_hi": float(ci.high), "below0": float(np.mean(hi < 0)),
            "above0": float(np.mean(lo > 0)), "w0": w0, "w0_ci": [w0l, w0h],
            "mean_positions": float(g.n_positions.mean()), "mean_positives": float(g.n_positive.mean()),
            "mean_clusters": float(g.n_clusters.mean()), "mean_near_prevalence": float(g.near_prevalence.mean()),
            "mean_n11": float(g.n11.mean())}


def main():
    t0 = time.time()
    args = [("registered", r) for r in range(N_REGISTERED)] + [("fresh", r) for r in range(N_FRESH)]
    chunks = [args[i:i + 100] for i in range(0, len(args), 100)]
    rows = []
    from multiprocessing import Pool
    with Pool(4) as pool:
        for res in pool.imap_unordered(work, chunks):
            rows.extend(res)
    df = pd.DataFrame(rows).sort_values(["kind", "rep"]).reset_index(drop=True)
    df.to_csv(os.path.join(OUT, "r1_b1_null.csv.gz"), index=False)
    reg = df[df.kind == "registered"].set_index("rep")
    stored = pd.read_csv(STORED_REPL)
    stored = stored[(stored.benchmark == "B1_cleavage") & (stored.strength == 0.0)].set_index("replicate")
    j = stored.join(reg, rsuffix="_re", how="inner")
    ext = pd.read_csv(STORED_EXT).set_index("benchmark").loc["B1_cleavage"]
    val = {"stored_replicates_compared": int(len(j)),
           "max_abs_diff_statistic": float(np.max(np.abs(j.statistic - j.statistic_re))),
           "max_abs_diff_ci_lo": float(np.max(np.abs(j.ci_lo - j.ci_lo_re))),
           "max_abs_diff_ci_hi": float(np.max(np.abs(j.ci_hi - j.ci_hi_re))),
           "n_positions_equal": bool((j.n_positions == j.n_positions_re).all()),
           "n_positive_equal": bool((j.n_positive == j.n_positive_re).all()),
           "detected_first50_stored": int(j.detected.sum()),
           "detected_first50_recomputed": int(((j.ci_lo_re > 0) | (j.ci_hi_re < 0)).sum()),
           "extension_stored_detected_of_300": int(ext.n_detected),
           "extension_recomputed_detected_of_300": int(((reg.ci_lo > 0) | (reg.ci_hi < 0)).sum())}
    out = {"label": "POST HOC revision analysis E_bootstrap_coverage, revision after verification (not registered)",
           "benchmark": "B1_cleavage (registered planted-artifact benchmark), strength 0",
           "statistic": SPEC["statistic"], "bootstrap_reps": REPS, "validation_against_registered": val,
           "registered_first50": summarise(df[(df.kind == "registered") & (df.rep < 50)]),
           "registered_300": summarise(df[df.kind == "registered"]),
           "fresh_10000": summarise(df[df.kind == "fresh"]),
           "seconds": round(time.time() - t0)}
    with open(os.path.join(OUT, "r1_b1_null.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
