"""Phase 2: run the synthetic ground-truth benchmark. Every parameter is read from the pre-registration;
none is redefined here. See scripts/preregister_phase2_synthetic_benchmark_2026-09-22.py for the full
specification of each benchmark's generative model, injection mechanism, statistic, null value and the
pass rule.

Gate A0 refuses to run if the protocol JSON's sha256 differs from what the pre-registration's own audit
recorded. Every random draw traces to the master seed via a deterministic offset scheme (documented at
`seed_for`), so a re-run with the same seed reproduces the summary CSV byte for byte (control PC2).

OUTPUTS (new names)
  results/phase2_synthetic_benchmark_replicates_2026-09-22.csv   one row per (benchmark, strength, replicate)
  results/phase2_synthetic_benchmark_summary_2026-09-22.csv      one row per benchmark: metrics + pass/fail
  results/phase2_synthetic_benchmark_2026-09-22_audit.json
Interpreter: project venv (numpy).
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROTOCOL = os.path.join(ROOT, "protocols", "phase2_synthetic_benchmark_preregistration_2026-09-22.json")
PREREG_AUDIT = os.path.join(ROOT, "results", "phase2_synthetic_benchmark_prereg_2026-09-22_audit.json")
OUT_REPL = os.path.join(ROOT, "results", "phase2_synthetic_benchmark_replicates_2026-09-22.csv")
OUT_SUM = os.path.join(ROOT, "results", "phase2_synthetic_benchmark_summary_2026-09-22.csv")
OUT_AUD = os.path.join(ROOT, "results", "phase2_synthetic_benchmark_2026-09-22_audit.json")


def sha(p):
    with open(p, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def seed_for(master, bench_idx, strength_idx, replicate_idx):
    """deterministic per-draw seed: distinct for every (benchmark, strength, replicate) triple"""
    return int(master) * 1_000_000 + bench_idx * 10_000 + strength_idx * 1000 + replicate_idx


def cluster_boot_ci(clusters, values, reps, seed, statistic):
    """protein-clustered percentile bootstrap CI of `statistic` (a callable over per-cluster aggregates
    supplied via `values`, which is (cluster_id, *payload) rows -- statistic receives the resampled rows)."""
    uniq = sorted(set(clusters))
    idx_by_cluster = {c: [] for c in uniq}
    for i, c in enumerate(clusters):
        idx_by_cluster[c].append(i)
    rng = np.random.default_rng(seed)
    n_clusters = len(uniq)
    draws = rng.integers(0, n_clusters, size=(reps, n_clusters))
    cluster_arr = np.array(uniq, dtype=object)
    out = np.empty(reps, dtype=float)
    for r in range(reps):
        chosen = cluster_arr[draws[r]]
        rows = []
        for c in chosen:
            rows.extend(idx_by_cluster[c])
        out[r] = statistic(rows)
    return float(np.quantile(out, 0.025)), float(np.quantile(out, 0.975))


def auc_rank(scores, labels):
    """rank-based AUC via Mann-Whitney U; labels in {0,1}, ties split evenly."""
    order = np.argsort(scores, kind="mergesort")
    ranks = np.empty(len(scores))
    s_sorted = scores[order]
    i = 0
    r = 1
    while i < len(s_sorted):
        j = i
        while j < len(s_sorted) and s_sorted[j] == s_sorted[i]:
            j += 1
        avg_rank = (r + (r + (j - i) - 1)) / 2
        ranks[order[i:j]] = avg_rank
        r += (j - i)
        i = j
    labels = np.asarray(labels)
    n_pos, n_neg = int(labels.sum()), int((1 - labels).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    sum_ranks_pos = ranks[labels == 1].sum()
    u = sum_ranks_pos - n_pos * (n_pos + 1) / 2
    return float(u / (n_pos * n_neg))


def average_precision(scores, labels):
    order = np.argsort(-scores, kind="mergesort")
    labels_sorted = np.asarray(labels)[order]
    n_pos = int(labels_sorted.sum())
    if n_pos == 0:
        return float("nan")
    tp = np.cumsum(labels_sorted)
    precision = tp / (np.arange(len(labels_sorted)) + 1)
    return float((precision * labels_sorted).sum() / n_pos)


def spearman(x, y):
    def rank(v):
        order = np.argsort(v, kind="mergesort")
        r = np.empty(len(v))
        r[order] = np.arange(1, len(v) + 1)
        return r
    rx, ry = rank(np.asarray(x, float)), rank(np.asarray(y, float))
    return float(np.corrcoef(rx, ry)[0, 1])


# ---------------------------------------------------------------- B1 cleavage
def cut_sites(length, rate, rng):
    return np.flatnonzero(rng.random(length) < rate)


def detectable_positions(length, cuts, lo, hi, max_missed):
    """returns a boolean array (length,) marking positions covered by >=1 peptide of length in [lo,hi]
    with <=max_missed missed cleavages, and an int array of each covered position's min distance to the
    nearer flanking cut boundary (0-indexed positions; boundaries are the cut-site index +1, and 0/length
    for the termini)."""
    bounds = np.concatenate(([0], cuts + 1, [length]))
    bounds = np.unique(bounds)
    det = np.zeros(length, dtype=bool)
    dist = np.full(length, 10 ** 9, dtype=np.int64)
    nb = len(bounds) - 1
    for i in range(nb):
        for k in range(0, max_missed + 1):
            j = i + k
            if j >= nb:
                break
            start, end = bounds[i], bounds[j + 1]
            plen = end - start
            if lo <= plen <= hi:
                seg = np.arange(start, end)
                d = np.minimum(seg - start, end - seg)
                det[start:end] = True
                dist[start:end] = np.minimum(dist[start:end], d)
    return det, dist


def run_b1(spec, strengths, rpr, master, bidx):
    rows = []
    for sidx, s in enumerate(strengths):
        for rep in range(rpr):
            seed = seed_for(master, bidx, sidx, rep)
            rng = np.random.default_rng(seed)
            plens = rng.integers(spec["protein_length_range"][0], spec["protein_length_range"][1] + 1,
                                 size=spec["n_proteins"])
            prot_ids, offsets, dists = [], [], []
            for pi, L in enumerate(plens):
                cuts = cut_sites(L, spec["cleavage_site_rate_per_residue"], rng)
                det, dist = detectable_positions(L, cuts, spec["detectable_peptide_length_range"][0],
                                                 spec["detectable_peptide_length_range"][1], spec["max_missed_cleavages"])
                is_cys = rng.random(L) < spec["cys_rate_per_residue"]
                cys_det = np.flatnonzero(is_cys & det)
                if len(cys_det) == 0:
                    continue
                prot_ids.extend([pi] * len(cys_det))
                offsets.extend(cys_det.tolist())
                dists.extend(dist[cys_det].tolist())
            prot_ids = np.array(prot_ids)
            dists = np.array(dists)
            near = dists <= spec["near_cut_window_residues"]
            n_bg = len(prot_ids)
            n_pos = max(1, int(round(spec["true_positive_fraction_of_detectable"] * n_bg)))
            pool_idx = np.arange(n_bg)
            use_near_pool = rng.random(n_pos) < s
            pos_idx = np.empty(n_pos, dtype=np.int64)
            near_pool = pool_idx[near]
            for k in range(n_pos):
                if use_near_pool[k] and len(near_pool) > 0:
                    pos_idx[k] = near_pool[rng.integers(0, len(near_pool))]
                else:
                    pos_idx[k] = pool_idx[rng.integers(0, n_bg)]
            pos_idx = np.unique(pos_idx)
            bg_mask = np.ones(n_bg, dtype=bool)
            bg_mask[pos_idx] = False

            def stat(row_idx):
                row_idx = np.asarray(row_idx)
                pos_rows = row_idx[~bg_mask[row_idx]]
                bg_rows = row_idx[bg_mask[row_idx]]
                a = int(near[pos_rows].sum()) + 0.5
                b = int((~near[pos_rows]).sum()) + 0.5
                c = int(near[bg_rows].sum()) + 0.5
                d = int((~near[bg_rows]).sum()) + 0.5
                return float(np.log2((a / b) / (c / d)))

            point = stat(np.arange(n_bg))
            lo, hi = cluster_boot_ci(prot_ids.tolist(), None, spec["_bootstrap_reps"], seed + 1, lambda ridx: stat(ridx))
            rows.append({"strength": s, "replicate": rep, "statistic": point, "ci_lo": lo, "ci_hi": hi,
                        "n_positions": n_bg, "n_positive": len(pos_idx)})
    return rows


# ---------------------------------------------------------------- B2 abundance
def run_b2(spec, strengths, rpr, master, bidx):
    rows = []
    n = spec["n_proteins"]
    for sidx, s in enumerate(strengths):
        for rep in range(rpr):
            seed = seed_for(master, bidx, sidx, rep)
            rng = np.random.default_rng(seed)
            abund = rng.lognormal(mean=0.0, sigma=1.2, size=n)
            rank_q = (np.argsort(np.argsort(abund)) + 1) / n
            p_det = (1 - s) * spec["base_detection_rate"] + s * rank_q
            detected = rng.random(n) < p_det
            prot_ids = np.arange(n)

            def stat(row_idx):
                row_idx = np.asarray(row_idx)
                sc, lb = abund[row_idx], detected[row_idx].astype(int)
                return auc_rank(sc, lb)

            point = stat(np.arange(n))
            lo, hi = cluster_boot_ci(prot_ids.tolist(), None, spec["_bootstrap_reps"], seed + 1, lambda ridx: stat(ridx))
            rows.append({"strength": s, "replicate": rep, "statistic": point, "ci_lo": lo, "ci_hi": hi,
                        "n_positions": n, "n_positive": int(detected.sum())})
    return rows


# ---------------------------------------------------------------- B3 multicys
def run_b3(spec, strengths, rpr, master, bidx):
    rows = []
    n_single, n_multi, n_prot = spec["n_peptides_single"], spec["n_peptides_multi"], spec["n_proteins"]
    for sidx, s in enumerate(strengths):
        for rep in range(rpr):
            seed = seed_for(master, bidx, sidx, rep)
            rng = np.random.default_rng(seed)
            # single-Cys peptides
            plen_s = rng.integers(spec["peptide_length_range"][0], spec["peptide_length_range"][1] + 1, n_single)
            true_pos_s = rng.integers(1, plen_s + 1)  # 1-indexed position within the peptide
            reported_s = true_pos_s.astype(float) / plen_s
            prot_s = rng.integers(0, n_prot, n_single)
            # multi-Cys peptides
            plen_m = rng.integers(spec["peptide_length_range"][0], spec["peptide_length_range"][1] + 1, n_multi)
            n_cys = rng.integers(spec["multicys_n_cys_range"][0], spec["multicys_n_cys_range"][1] + 1, n_multi)
            reported_m = np.empty(n_multi, dtype=float)
            use_nterm = rng.random(n_multi) < s
            for i in range(n_multi):
                L, k = plen_m[i], n_cys[i]
                positions = rng.choice(np.arange(1, L + 1), size=k, replace=False)
                true_site = positions[rng.integers(0, k)]
                reported_site = positions.min() if use_nterm[i] else true_site
                reported_m[i] = reported_site / L
            prot_m = rng.integers(0, n_prot, n_multi)

            reported = np.concatenate([reported_s, reported_m])
            is_multi = np.concatenate([np.zeros(n_single, dtype=bool), np.ones(n_multi, dtype=bool)])
            prot_ids = np.concatenate([prot_s, prot_m])

            def stat(row_idx):
                row_idx = np.asarray(row_idx)
                mask_multi = is_multi[row_idx]
                m1 = reported[row_idx][mask_multi]
                m0 = reported[row_idx][~mask_multi]
                if len(m1) == 0 or len(m0) == 0:
                    return 0.0
                return float(m1.mean() - m0.mean())

            point = stat(np.arange(len(reported)))
            lo, hi = cluster_boot_ci(prot_ids.tolist(), None, spec["_bootstrap_reps"], seed + 1, lambda ridx: stat(ridx))
            rows.append({"strength": s, "replicate": rep, "statistic": point, "ci_lo": lo, "ci_hi": hi,
                        "n_positions": len(reported), "n_positive": n_multi})
    return rows


# ---------------------------------------------------------------- B4 depletion
def run_b4(spec, strengths, rpr, master, bidx):
    rows = []
    n_sites, n_prot, p_mod = spec["n_sites"], spec["n_proteins"], spec["true_modified_prevalence"]
    for sidx, s in enumerate(strengths):
        for rep in range(rpr):
            seed = seed_for(master, bidx, sidx, rep)
            rng = np.random.default_rng(seed)
            is_mod = rng.random(n_sites) < p_mod
            survive = is_mod | (rng.random(n_sites) < (1 - s))
            prot_ids = rng.integers(0, n_prot, n_sites)
            rep_prot = prot_ids[survive]
            rep_mod = is_mod[survive]

            def stat(row_idx):
                row_idx = np.asarray(row_idx)
                m = rep_mod[row_idx]
                return float(m.sum()) / max(1, len(row_idx))

            point = stat(np.arange(len(rep_prot)))
            lo, hi = cluster_boot_ci(rep_prot.tolist(), None, spec["_bootstrap_reps"], seed + 1, lambda ridx: stat(ridx))
            rows.append({"strength": s, "replicate": rep, "statistic": point, "ci_lo": lo, "ci_hi": hi,
                        "n_positions": int(survive.sum()), "n_positive": int(rep_mod.sum())})
    return rows


BENCHMARKS = [("B1_cleavage", run_b1), ("B2_abundance", run_b2), ("B3_multicys", run_b3), ("B4_depletion", run_b4)]
NULL_EXCLUDES = "detected"


def main():
    for p in (OUT_REPL, OUT_SUM, OUT_AUD):
        if os.path.exists(p):
            sys.exit(f"REFUSE: {p} exists")
    with open(PREREG_AUDIT, encoding="utf-8") as fh:
        pa = json.load(fh)
    if pa["protocol_sha256"] != sha(PROTOCOL):
        sys.exit("REFUSE: pre-registration protocol changed since it was registered")
    with open(PROTOCOL, encoding="utf-8") as fh:
        P = json.load(fh)
    master = P["master_seed"]
    strengths = P["strength_grid"]
    rpr = P["replicates_per_strength"]
    boot_reps = P["bootstrap_reps"]
    pass_rule = P["pass_rule"]

    all_rows, summary_rows = [], []
    for bidx, (name, fn) in enumerate(BENCHMARKS):
        spec = dict(P["benchmarks"][name])
        spec["_bootstrap_reps"] = boot_reps
        null_v = spec["null_value"]
        sign = 1.0 if spec["score_sign"].startswith("+1") else -1.0
        rows = fn(spec, strengths, rpr, master, bidx)
        for r in rows:
            r["benchmark"] = name
            r["detected"] = int(not (r["ci_lo"] <= null_v <= r["ci_hi"]))
            r["score"] = sign * r["statistic"]
        all_rows.extend(rows)

        sens = {}
        for s in strengths:
            sub = [r for r in rows if r["strength"] == s]
            sens[s] = sum(r["detected"] for r in sub) / len(sub)
        fpr0 = sens[strengths[0]]
        sens_s1 = sens[strengths[-1]]
        scores = np.array([r["score"] for r in rows])
        labels = np.array([1 if r["strength"] > 0 else 0 for r in rows])
        auroc = auc_rank(scores, labels)
        auprc = average_precision(scores, labels)
        rho = spearman(strengths, [sens[s] for s in strengths])
        thr = next((s for s in strengths if sens[s] >= pass_rule["sensitivity_at_s1_min"]), None)
        passed = (fpr0 <= pass_rule["fpr_at_s0_max"] and sens_s1 >= pass_rule["sensitivity_at_s1_min"]
                  and auroc >= pass_rule["auroc_min"] and rho >= pass_rule["spearman_strength_vs_sensitivity_min"])
        summary_rows.append({"benchmark": name, "artefact": spec["artefact"], "fpr_at_s0": round(fpr0, 4),
                             "sensitivity_at_s1": round(sens_s1, 4), "auroc": round(auroc, 4),
                             "auprc": round(auprc, 4), "spearman_strength_vs_sensitivity": round(rho, 4),
                             "detection_threshold_strength": thr, "pass": int(passed),
                             **{f"sensitivity_at_{s}": round(sens[s], 4) for s in strengths}})

    with open(OUT_REPL, "w", newline="", encoding="utf-8") as fh:
        fields = ["benchmark", "strength", "replicate", "statistic", "score", "ci_lo", "ci_hi", "detected",
                  "n_positions", "n_positive"]
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in all_rows:
            w.writerow({k: r[k] for k in fields})
    with open(OUT_SUM, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary_rows[0].keys()))
        w.writeheader()
        w.writerows(summary_rows)

    audit = {"script": "scripts/run_phase2_synthetic_benchmark_2026-09-22.py", "script_sha256": sha(os.path.abspath(__file__)),
             "protocol_sha256": sha(PROTOCOL), "prereg_audit_sha256": sha(PREREG_AUDIT),
             "master_seed": master, "strength_grid": strengths, "replicates_per_strength": rpr,
             "bootstrap_reps": boot_reps, "pass_rule": pass_rule,
             "n_replicate_rows": len(all_rows), "summary": summary_rows,
             "outputs": {os.path.relpath(OUT_REPL, ROOT): sha(OUT_REPL), os.path.relpath(OUT_SUM, ROOT): sha(OUT_SUM)}}
    with open(OUT_AUD, "w", encoding="utf-8") as fh:
        json.dump(audit, fh, ensure_ascii=False, indent=1)
    print(json.dumps(summary_rows, indent=1))


if __name__ == "__main__":
    main()
