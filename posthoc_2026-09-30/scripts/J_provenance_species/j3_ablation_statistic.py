# -*- coding: utf-8 -*-
"""J3. A threshold-free statistic for the three-arm ablation (Table 2, Fig. 7b), or why it cannot be run.

POST HOC revision analysis (2026-09-30), written in response to a pre-submission review; not
registered, not pre-specified.

QUESTION
Table 2 and Fig. 7b rank the public-refit arms by the registered "top-100 distal log2 odds ratio":
the Haldane log2 OR of carrying a trypsin cleavage site 6-12 residues from the cysteine, for the 100
top-scored cysteines of a cohort against the rest. The reviewer asks for a threshold-free statistic.
That needs the per-cysteine member scores of every arm and each cysteine's distal flag.

PART 1  Inventory. The per-cysteine out-of-fold scores of the public refit (one CSV per cohort and
        arm, columns accession, position, component, fold, label, blend, lgb_rank, lgb_bin,
        fusion_mlp) were written to external/public_refit_out_2026-09-21/ in the internal working
        tree; the ranking-member reading that Table 2 prints was written to
        public_refit_ablation_lgb_rank_2026-09-24.csv. This part checks which of these, and of the
        inputs needed to compute distal flags, exist on this machine, and records size and sha256 from
        the internal artifact manifest where the manifest has them.

PART 2  What the stored summary values already show about the top-100 statistic (no per-cysteine
        data needed). The registered analysis stores, per arm, n, the distal rate among the top 100
        (2 decimals, i.e. an integer count k out of 100) and the distal rate among the rest (4
        decimals). Because every arm scores the same cysteines, the total number D of distal cysteines
        is common to all arms; it is recovered as the one integer consistent with every stored row,
        and cross-checked against the label-stratum rates of Supplemental Note 9. With N, D and k the
        2x2 table is determined, so the stored log2 OR must be reproduced exactly: this verifies that
        each Table 2 value is a function of a single integer k. The script then reports
          * the statistic's range: its value at k = 0..100, its ceiling (k = 100) and its step size;
          * the no-association reference: k ~ Hypergeometric(N, D, 100);
          * a label-inheritance reference: the value expected if the 100 top-scored cysteines were 100
            label-positive cysteines taken at random with respect to the distal flag,
            k ~ Hypergeometric(N_pos, D_pos, 100). A model that ranks labelled cysteines first, and
            knows nothing about cleavage geometry within a label class, produces a top-100 statistic in
            this range, because labelled cysteines carry the distal flag more often (Supplemental
            Note 9). This is a descriptive reference, not a test: a model whose top 100 includes
            unlabelled cysteines has a lower label-inheritance expectation.

PART 2b (revision after verification, round 1; post hoc) Every arm reading against every reference.
        The round-1 text read the label-inheritance reference selectively. This part places each of the
        12 primary readings (chemistry set, minus win7/win10, minus all twenty; ranking member and
        registered blend; two cohorts), and the control arms for completeness, against
          * no association, k ~ Hypergeometric(N, D, 100);
          * label inheritance with a top 100 made of 100f labelled and 100(1-f) unlabelled cysteines drawn
            without regard to cleavage geometry, f = 1.0 (round-1 reference), 0.9, 0.8, 0.7, 0.5: k is the
            sum of two independent hypergeometric counts (exact convolution);
        and reports the exact one-sided P(K >= k_obs), the 95% range of K, and whether k_obs lies below,
        at the lower limit of, within, at the upper limit of, or above that range. It also records whether
        each arm's stored protein-clustered bootstrap interval overlaps the reference's log2 OR range. All
        of this is descriptive: the label composition of each top 100 is unknown, the tail probabilities are
        unadjusted for the 12 readings, and one cysteine moves the statistic by 0.2-0.75 log2 units near
        the observed values.

PART 2c (revision after verification, round 1; post hoc) Protein clustering of the label-inheritance
        reference, rice only. The references above draw the top 100 cysteine by cysteine. If a model's
        top 100 is concentrated in few proteins, the reference is wider. Distal flags need the full protein
        sequence (the released 15-residue windows are too short), and no human proteome file is on this
        machine, so this is computed for rice only, on the cysteines of the rice deposit whose protein is
        in UniProt UP000059680 release 2026_03 (W/external; 88% of the analysed cohort; windows checked to
        match). Two references for 100 labelled cysteines: drawn one by one (hypergeometric), and drawn
        protein by protein (proteins sampled without replacement with probability proportional to their
        number of labelled cysteines, every labelled cysteine of a drawn protein taken, the last protein
        subsampled to reach exactly 100; 20,000 draws, seed 20260930). The second is the extreme of
        clustering for a top 100 made of labelled cysteines.

PART 3  The threshold-free statistic itself, implemented and self-tested on synthetic data, ready to run
        on the refit files once they are deposited (python j3_ablation_statistic.py --oof-dir <dir>
        --hsa <hsa.fasta.gz> --osa <osa.fasta.gz>):
          * AUC of each member's score for distal-band membership over all cysteines;
          * the same AUC within each label stratum (label constant, so not inherited from the label;
            the estimand of Fig. 7a);
          * within-protein AUC for distal-band membership (mean over proteins holding both distal and
            non-distal cysteines);
          * paired differences arm minus chemistry-set control;
        all with protein-clustered bootstrap intervals (5,000 replicates, seed 20260921, the seed of the
        registered refit analysis). Distal flags follow the registered instrument
        (scripts/run_self_audit_five_proteases.py band_flags, DISTAL = (6, 12), trypsin rule of
        run_cross_protease_detectability_probe.PROTEASES), re-implemented verbatim below.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import re
import sys

import numpy as np
from scipy.stats import hypergeom
from sklearn.metrics import roc_auc_score

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from jcommon import MCP, OUT, REPO, W, Inputs, write_csv, write_fragment, write_json  # noqa: E402

HALDANE = 0.5
TOP_K = 100
DISTAL = (6, 12)
PROXIMAL = 5
BOOT = 5000
SEED = 20260921          # seed of the registered refit analysis (analyse_public_refit_2026-09-21.py)
SELFTEST_SEED = 20260930
COHORTS = ("human_PXD044043", "rice_PXD072089")
ARMS = ("control", "drop_w710", "drop_all_basic", "label_permuted", "detect_only", "full")
MEMBERS = ("lgb_rank", "blend", "lgb_bin", "fusion_mlp")
TRYPSIN = {"residues": "KR", "side": "after", "block_proline": True}


# ----------------------------------------------------------------------------- instrument (verbatim)
def boundaries_for(sequence, rule):
    """run_cross_protease_detectability_probe.py lines 122-137."""
    n = len(sequence)
    cuts = [0]
    residues, side, block = rule["residues"], rule["side"], rule["block_proline"]
    for i, residue in enumerate(sequence):
        if residue not in residues:
            continue
        cut = i + 1 if side == "after" else i
        if cut <= 0 or cut >= n:
            continue
        if block and sequence[cut] == "P":
            continue
        cuts.append(cut)
    if cuts[-1] != n:
        cuts.append(n)
    return sorted(set(cuts))


def band_flags(sequence, position, rule):
    """run_self_audit_five_proteases.py lines 105-126: (distal 6-12, proximal <=5, peptide length)."""
    i = position - 1
    cuts = boundaries_for(sequence, rule)
    interior = np.asarray([c for c in cuts if 0 < c < len(sequence)], dtype=int)
    cuts_full = np.concatenate([[0], interior, [len(sequence)]])
    idx = np.searchsorted(cuts_full, i, side="right")
    length = int(cuts_full[idx] - cuts_full[idx - 1])
    if interior.size == 0:
        return 0, 0, length
    d = np.abs(interior - i)
    return int(((d >= DISTAL[0]) & (d <= DISTAL[1])).any()), int((d <= PROXIMAL).any()), length


def second_cys_trypsinP(sequence, i):
    """revision round 2: 1 when the minimal Trypsin/P segment (cleavage after every K/R) around residue i holds
    a second cysteine, i.e. the cysteine can never be the only cysteine of a tryptic peptide (as j5)."""
    s = 0
    for j in range(i - 1, -1, -1):
        if sequence[j] in "KR":
            s = j + 1
            break
    e = len(sequence)
    for j in range(i + 1, len(sequence)):
        if sequence[j] in "KR":
            e = j + 1
            break
    return int(sequence[s:e].count("C") >= 2)


def log2_or(n11, n10, n01, n00):
    """run_phase2_claims_under_detectability_control.py line 192 (Haldane 0.5)."""
    a, b, c, d = n11 + HALDANE, n10 + HALDANE, n01 + HALDANE, n00 + HALDANE
    return float(np.log2((a * d) / (b * c)))


def top100_lor(k, D, N, K=TOP_K):
    # table(y=in_top, attribute=distal): n11 distal&top, n10 distal&rest, n01 nondistal&top, n00 nondistal&rest
    return log2_or(k, D - k, K - k, (N - K) - (D - k))


# ----------------------------------------------------------------------------- statistics
def make_weighted_auc(y, scores):
    """Weighted AUC closure, as run_phase2_claims_under_detectability_control.make_weighted_auc."""
    order = np.argsort(np.asarray(scores, dtype=float), kind="stable")
    y_sorted = np.asarray(y)[order]
    s_sorted = np.asarray(scores, dtype=float)[order]
    starts = np.flatnonzero(np.concatenate([[True], s_sorted[1:] != s_sorted[:-1]]))
    pos = (y_sorted == 1).astype(float)
    neg = (y_sorted == 0).astype(float)

    def auc(weights):
        w = np.asarray(weights, dtype=float)[order]
        bp = np.add.reduceat(w * pos, starts)
        bn = np.add.reduceat(w * neg, starts)
        tp, tn = bp.sum(), bn.sum()
        if tp <= 0 or tn <= 0:
            return float("nan")
        before = np.concatenate([[0.0], np.cumsum(bn)[:-1]])
        return float(np.dot(bp, before + 0.5 * bn) / (tp * tn))

    return auc


def per_protein_auc(groups, y, s):
    out = {}
    for g in np.unique(groups):
        m = groups == g
        yy = y[m]
        if 0 < yy.sum() < yy.size:
            out[g] = float(roc_auc_score(yy, s[m]))
    return out


def cluster_weights(groups, reps, rng):
    uniq, inv = np.unique(groups, return_inverse=True)
    counts = rng.multinomial(len(uniq), np.full(len(uniq), 1.0 / len(uniq)), size=reps)
    return counts[:, inv], uniq, counts


def threshold_free_table(cohort, acc, y, dist, scores_by_arm, members, boot=BOOT, seed=SEED):
    """scores_by_arm[arm][member] -> score vector aligned with acc/y/dist."""
    rows = []
    rng = np.random.default_rng(seed)
    wmat, uniq, counts = cluster_weights(acc, boot, rng)   # shared resamples: paired comparisons
    prot_index = {p: i for i, p in enumerate(uniq)}
    strata = {"all": np.ones(len(y), bool), "label_positive": y == 1, "label_negative": y == 0}
    ref = {}
    for arm, byMember in scores_by_arm.items():
        for member in members:
            s = byMember[member]
            rec = {"cohort": cohort, "arm": arm, "member": member}
            for sname, m in strata.items():
                f = make_weighted_auc(dist[m], s[m])
                point = f(np.ones(int(m.sum())))
                draws = np.asarray([f(w[m]) for w in wmat])
                ref[(arm, member, sname)] = draws
                rec[f"auc_distal_{sname}"] = point
                rec[f"auc_distal_{sname}_ci_low"] = float(np.nanquantile(draws, 0.025))
                rec[f"auc_distal_{sname}_ci_high"] = float(np.nanquantile(draws, 0.975))
            pp = per_protein_auc(acc, dist, s)
            keys = list(pp)
            vals = np.asarray([pp[k] for k in keys])
            idx = np.asarray([prot_index[k] for k in keys])
            cw = counts[:, idx].astype(float)
            draws = (cw @ vals) / cw.sum(axis=1)
            ref[(arm, member, "within_protein")] = draws
            rec["auc_distal_within_protein_mean"] = float(vals.mean())
            rec["auc_distal_within_protein_ci_low"] = float(np.quantile(draws, 0.025))
            rec["auc_distal_within_protein_ci_high"] = float(np.quantile(draws, 0.975))
            rec["n_proteins_within_protein"] = len(keys)
            rows.append(rec)
    diffs = []
    for arm in scores_by_arm:
        if arm == "control":
            continue
        for member in members:
            for sname in ("all", "label_positive", "label_negative", "within_protein"):
                a, c = ref[(arm, member, sname)], ref[("control", member, sname)]
                d = a - c
                pa = [r for r in rows if r["arm"] == arm and r["member"] == member][0]
                pc = [r for r in rows if r["arm"] == "control" and r["member"] == member][0]
                key = f"auc_distal_{sname}" if sname != "within_protein" else "auc_distal_within_protein_mean"
                diffs.append({"cohort": cohort, "arm": arm, "member": member, "stratum": sname,
                              "diff_vs_control": pa[key] - pc[key],
                              "ci_low": float(np.nanquantile(d, 0.025)), "ci_high": float(np.nanquantile(d, 0.975))})
    return rows, diffs


# ----------------------------------------------------------------------------- PART 2 helpers
def parse_note11(text):
    """Return {'lgb_rank': rows, 'blend': rows} from the two arm tables of Supplemental Note 11."""
    blocks, cur = [], None
    for line in text.splitlines():
        if line.startswith("| cohort | arm | n | top100_distal_rate"):
            cur = []
            blocks.append(cur)
            continue
        if cur is not None:
            if line.startswith("|---"):
                continue
            if line.startswith("| ") and ("human_PXD044043" in line or "rice_PXD072089" in line):
                cells = [c.strip() for c in line.strip().strip("|").split("|")]
                cur.append(dict(zip(["cohort", "arm", "n", "top100_distal_rate", "rate_rest", "log2_or", "ci_low",
                                     "ci_high", "within_protein_auc", "n_proteins_scored"], cells)))
            else:
                cur = None
    assert len(blocks) == 2, len(blocks)
    return {"lgb_rank": blocks[0], "blend": blocks[1]}


def feasible_D(rows, N):
    lo, hi = -math.inf, math.inf
    for r in rows:
        k = int(round(float(r["top100_distal_rate"]) * 100))
        rr = float(r["rate_rest"])
        rest = N - TOP_K
        a = math.ceil((rr - 0.00005) * rest - 1e-9) + k
        b = math.floor((rr + 0.00005) * rest + 1e-9) + k
        lo, hi = max(lo, a), min(hi, b)
    return lo, hi


def hyper_q(M, n, N, q):
    return int(hypergeom.ppf(q, M, n, N))


# ----------------------------------------------------------------------------- PART 2b helpers
PRIMARY_ARMS = ("control", "drop_w710", "drop_all_basic")
MIX_F = (1.0, 0.9, 0.8, 0.7, 0.5)
RICE_FASTA_2026 = W + "/external/UP000059680_39947.fasta.gz"
SITE_TABLE = REPO + "/external/public_cohort_site_scores.csv"
CLUSTER_DRAWS = 20000
CLUSTER_SEED = 20260930


def mixture_pmf(n_pos, d_pos, n_neg, d_neg, f, K=TOP_K):
    """pmf of k = distal count among a top K made of round(K f) labelled and K - round(K f) unlabelled
    cysteines, each group drawn without replacement and without regard to the distal flag."""
    kp = int(round(K * f))
    kn = K - kp
    pp = hypergeom.pmf(np.arange(kp + 1), n_pos, d_pos, kp) if kp > 0 else np.array([1.0])
    pn = hypergeom.pmf(np.arange(kn + 1), n_neg, d_neg, kn) if kn > 0 else np.array([1.0])
    pmf = np.convolve(pp, pn)
    return pmf / pmf.sum()


def pmf_q(pmf, q):
    """smallest k with P(K <= k) >= q (the definition of scipy's ppf)."""
    return int(np.searchsorted(np.cumsum(pmf), q - 1e-12, side="left"))


def position_vs_range(k, lo, hi):
    if k < lo:
        return "below"
    if k == lo:
        return "at lower limit"
    if k < hi:
        return "within"
    if k == hi:
        return "at upper limit"
    return "above"


def arm_vs_references(summary, tabs):
    rows = []
    for c, s in summary.items():
        N, D = s["N"], s["D_distal_total"]
        refs = [("no association", None, hypergeom.pmf(np.arange(TOP_K + 1), N, D, TOP_K))]
        for f in MIX_F:
            refs.append((f"label inheritance, {int(round(100 * f))} of 100 labelled", f,
                         mixture_pmf(s["N_label_positive"], s["D_label_positive"], s["N_label_negative"],
                                     s["D_label_negative"], f)))
        for member in ("lgb_rank", "blend"):
            for r in tabs[member]:
                if r["cohort"] != c:
                    continue
                k = int(round(float(r["top100_distal_rate"]) * 100))
                for rname, f, pmf in refs:
                    lo, hi = pmf_q(pmf, 0.025), pmf_q(pmf, 0.975)
                    lor_lo, lor_hi = top100_lor(lo, D, N), top100_lor(hi, D, N)
                    rows.append({
                        "cohort": c, "member": member, "arm": r["arm"], "primary_reading": int(r["arm"] in PRIMARY_ARMS),
                        "k_distal_in_top100": k, "log2_or": float(r["log2_or"]), "ci_low": float(r["ci_low"]),
                        "ci_high": float(r["ci_high"]), "reference": rname, "f_labelled_in_top100": "" if f is None else f,
                        "ref_mean_k": float(np.dot(np.arange(TOP_K + 1), pmf)), "ref_k_q025": lo, "ref_k_q975": hi,
                        "ref_log2_or_q025": lor_lo, "ref_log2_or_q975": lor_hi,
                        "P_K_ge_k_obs": float(pmf[k:].sum()),
                        "position_vs_reference_range": position_vs_range(k, lo, hi),
                        "stored_ci_overlaps_reference_log2_or_range": int(float(r["ci_low"]) <= lor_hi and float(r["ci_high"]) >= lor_lo),
                    })
    return rows


def count_positions(rows, reference):
    sel = [r for r in rows if r["primary_reading"] and r["reference"] == reference]
    out = {"n_readings": len(sel)}
    for pos in ("below", "at lower limit", "within", "at upper limit", "above"):
        hit = [r for r in sel if r["position_vs_reference_range"] == pos]
        out[pos] = {"n": len(hit), "readings": [f"{r['cohort'].split('_')[0]} {r['member']} {r['arm']} k={r['k_distal_in_top100']} "
                                                f"P={r['P_K_ge_k_obs']:.4f}" for r in hit]}
    out["n_P_lt_0.05"] = sum(r["P_K_ge_k_obs"] < 0.05 for r in sel)
    out["n_stored_ci_overlapping_reference_range"] = sum(r["stored_ci_overlaps_reference_log2_or_range"] for r in sel)
    out["stored_ci_not_overlapping"] = [f"{r['cohort'].split('_')[0]} {r['member']} {r['arm']}" for r in sel
                                        if not r["stored_ci_overlaps_reference_log2_or_range"]]
    return out


# ----------------------------------------------------------------------------- PART 2c helpers
def pps_cluster_draws(n_lab, d_lab, draws, seed, K=TOP_K, chunk=1000):
    """k for a top K of labelled cysteines drawn protein by protein: proteins sampled without replacement
    with probability proportional to their labelled count (Efraimidis-Spirakis keys), all labelled
    cysteines of each drawn protein taken, the last protein subsampled (hypergeometric) to reach K."""
    rng = np.random.default_rng(seed)
    n_lab = np.asarray(n_lab, dtype=float)
    d_lab = np.asarray(d_lab, dtype=int)
    out = np.empty(draws, dtype=int)
    done = 0
    while done < draws:
        m = min(chunk, draws - done)
        keys = np.log(rng.random((m, n_lab.size))) / n_lab            # larger key = drawn earlier
        order = np.argsort(-keys, axis=1)
        ns = n_lab[order].astype(int)
        ds = d_lab[order]
        cum = np.cumsum(ns, axis=1)
        last = np.argmax(cum >= K, axis=1)                             # index of the protein that reaches K
        rows_ = np.arange(m)
        before = np.where(last > 0, cum[rows_, last - 1], 0)
        dcum = np.cumsum(ds, axis=1)
        d_before = np.where(last > 0, dcum[rows_, last - 1], 0)
        need = K - before
        partial = rng.hypergeometric(ds[rows_, last], ns[rows_, last] - ds[rows_, last], need)
        out[done:done + m] = d_before + partial
        done += m
    return out


def rice_cluster_sensitivity(inp, summary, tabs):
    seqs = read_fasta_gz(inp.use(RICE_FASTA_2026))
    rows = [r for r in inp.read_csv(SITE_TABLE) if r["dataset"] == "rice_PXD072089"]
    kept, win_mismatch = [], 0
    for r in rows:
        a, p, w = r["accession"], int(r["position"]), r["context_15aa_or_terminal_shorter"]
        s = seqs.get(a)
        if s is None:
            continue
        if not (0 < p <= len(s)) or s[max(0, p - 8):p + 7] != w:
            win_mismatch += 1
            continue
        kept.append((a, p, int(r["observed_in_retrospective_dataset"]), band_flags(s, p, TRYPSIN)[0]))
    acc = np.asarray([x[0] for x in kept])
    y = np.asarray([x[2] for x in kept])
    dist = np.asarray([x[3] for x in kept])
    lab_prot = {}
    for a, yy, dd in zip(acc, y, dist):
        if yy == 1:
            n, d = lab_prot.get(a, (0, 0))
            lab_prot[a] = (n + 1, d + int(dd))
    n_lab = np.asarray([v[0] for v in lab_prot.values()])
    d_lab = np.asarray([v[1] for v in lab_prot.values()])
    Npos, Dpos = int(n_lab.sum()), int(d_lab.sum())
    s_full = summary["rice_PXD072089"]
    ind = hypergeom.pmf(np.arange(TOP_K + 1), Npos, Dpos, TOP_K)
    ks = pps_cluster_draws(n_lab, d_lab, CLUSTER_DRAWS, CLUSTER_SEED)
    clu = np.bincount(ks, minlength=TOP_K + 1) / ks.size
    res = {
        "scope": "rice deposit cysteines whose protein is in UniProt UP000059680 release 2026_03 and whose released window "
                 "matches that sequence; a subset of the analysed cohort (the refit used an earlier rice proteome file)",
        "rice_rows_in_site_table": len(rows), "cysteines_with_sequence": len(kept), "window_mismatches": win_mismatch,
        "analysed_cohort_N": s_full["N"], "share_of_analysed_cohort": len(kept) / s_full["N"],
        "labelled": Npos, "labelled_proteins": int(n_lab.size),
        "labelled_per_protein_distribution": {str(k): int(v) for k, v in sorted(zip(*np.unique(n_lab, return_counts=True)))},
        "distal_rate_labelled_subset": Dpos / Npos, "distal_rate_unlabelled_subset": float(dist[y == 0].mean()),
        "distal_rate_labelled_cohort_note9": s_full["D_label_positive"] / s_full["N_label_positive"],
        "distal_rate_unlabelled_cohort_note9": s_full["D_label_negative"] / s_full["N_label_negative"],
        "reference_independent": {"mean_k": float(np.dot(np.arange(TOP_K + 1), ind)), "k_q025": pmf_q(ind, 0.025),
                                  "k_q975": pmf_q(ind, 0.975), "var_k": float(np.dot((np.arange(TOP_K + 1) - np.dot(np.arange(TOP_K + 1), ind)) ** 2, ind))},
        "reference_protein_clustered": {"mean_k": float(ks.mean()), "k_q025": int(np.quantile(ks, 0.025, method="inverted_cdf")),
                                        "k_q975": int(np.quantile(ks, 0.975, method="inverted_cdf")), "var_k": float(ks.var()),
                                        "draws": CLUSTER_DRAWS, "seed": CLUSTER_SEED},
    }
    res["variance_ratio_clustered_vs_independent"] = res["reference_protein_clustered"]["var_k"] / res["reference_independent"]["var_k"]
    # one-way ANOVA intraclass correlation of the distal flag among labelled cysteines, clusters = proteins
    p_ = Dpos / Npos
    kc = n_lab.size
    msb = float(np.sum(n_lab * (d_lab / n_lab - p_) ** 2) / (kc - 1))
    msw = float((d_lab.sum() - np.sum(d_lab ** 2 / n_lab)) / (Npos - kc))
    n0 = float((Npos - np.sum(n_lab ** 2) / Npos) / (kc - 1))
    res["icc_distal_among_labelled_within_protein_anova"] = (msb - msw) / (msb + (n0 - 1) * msw)
    res["proteins_with_one_labelled_cysteine"] = int((n_lab == 1).sum())
    obs = []
    for member in ("lgb_rank", "blend"):
        for r in tabs[member]:
            if r["cohort"] != "rice_PXD072089" or r["arm"] not in PRIMARY_ARMS:
                continue
            k = int(round(float(r["top100_distal_rate"]) * 100))
            obs.append({"member": member, "arm": r["arm"], "k": k,
                        "P_ge_independent_subset": float(ind[k:].sum()), "P_ge_protein_clustered_subset": float(clu[k:].sum()),
                        "position_independent": position_vs_range(k, res["reference_independent"]["k_q025"], res["reference_independent"]["k_q975"]),
                        "position_protein_clustered": position_vs_range(k, res["reference_protein_clustered"]["k_q025"], res["reference_protein_clustered"]["k_q975"])})
    res["rice_primary_readings"] = obs
    return res


# ----------------------------------------------------------------------------- self-test
def selftest():
    rng = np.random.default_rng(SELFTEST_SEED)
    out = {}
    # 1 instrument band flags on hand-made sequences
    seq = "AAAAACAAAAAKAAAAAAAAAAAAAAAAAR"   # C at 6, K at 12 -> cut 12, d = |12-5| = 7 -> distal
    out["band_flags_case1"] = band_flags(seq, 6, TRYPSIN)
    seq2 = "AAAAACAKPAAAAAAAAAAAAAAAAAAAAR"  # K followed by P: no cut; only the protein end remains
    out["band_flags_case2"] = band_flags(seq2, 6, TRYPSIN)
    seq3 = "AAAAACAKAAAAAAAAAAAAAAAAAAAAAR"  # K at 8 -> cut 8, d = 2 -> proximal only
    out["band_flags_case3"] = band_flags(seq3, 6, TRYPSIN)
    ok1 = out["band_flags_case1"][0] == 1 and out["band_flags_case2"][0] == 0 and out["band_flags_case3"][:2] == (0, 1)
    # 2 weighted AUC equals sklearn
    n = 3000
    y = rng.integers(0, 2, n)
    s = rng.normal(size=n) + 0.8 * y
    s[:200] = np.round(s[:200], 1)        # ties
    f = make_weighted_auc(y, s)
    ok2 = abs(f(np.ones(n)) - roc_auc_score(y, s)) < 1e-12
    # 3 synthetic cohort: label depends on distal; one arm adds distal signal within labels
    P, C = 400, 12
    acc = np.repeat(np.array([f"P{i:04d}" for i in range(P)]), C)
    dist = (rng.random(P * C) < 0.79).astype(int)
    lab = (rng.random(P * C) < np.where(dist == 1, 0.20, 0.10)).astype(int)
    base = 2.0 * lab + rng.normal(size=P * C)
    arms = {"control": {"lgb_rank": base + 0.0 * dist}, "drop_all_basic": {"lgb_rank": base + 0.8 * dist}}
    rows, diffs = threshold_free_table("synthetic", acc, lab, dist, arms, ("lgb_rank",), boot=400, seed=SELFTEST_SEED)
    r0 = [r for r in rows if r["arm"] == "control"][0]
    r1 = [r for r in rows if r["arm"] == "drop_all_basic"][0]
    # control: no within-stratum association (CI covers 0.5) but an overall association inherited from the label
    ok3 = (r0["auc_distal_label_positive_ci_low"] <= 0.5 <= r0["auc_distal_label_positive_ci_high"]
           and r0["auc_distal_label_negative_ci_low"] <= 0.5 <= r0["auc_distal_label_negative_ci_high"]
           and r0["auc_distal_all_ci_low"] > 0.5
           and r1["auc_distal_label_negative_ci_low"] > 0.5)
    # 4 within-protein mean equals mean of sklearn per-protein AUCs
    pp = per_protein_auc(acc, dist, base)
    ok4 = abs(np.mean(list(pp.values())) - r0["auc_distal_within_protein_mean"]) < 1e-12
    # 5 Haldane top-100 formula reproduces the registered instrument on a toy table
    ok5 = abs(top100_lor(93, 14670, 18567) - log2_or(93, 14670 - 93, 7, 18467 - 14577)) < 1e-12
    # 6 (revision) the mixture reference reduces to the single hypergeometric at f = 1 and f = 0
    m1 = mixture_pmf(3206, 2877, 15361, 11793, 1.0)
    m0 = mixture_pmf(3206, 2877, 15361, 11793, 0.0)
    ok6 = (np.allclose(m1, hypergeom.pmf(np.arange(101), 3206, 2877, 100), atol=1e-14)
           and np.allclose(m0, hypergeom.pmf(np.arange(101), 15361, 11793, 100), atol=1e-14))
    # 7 (revision) protein-by-protein draws reduce to simple random sampling when every protein has one
    #   labelled cysteine: mean and variance match the hypergeometric within Monte Carlo error
    n1 = np.ones(700, dtype=int)
    d1 = (rng.random(700) < 0.84).astype(int)
    kk = pps_cluster_draws(n1, d1, 20000, SELFTEST_SEED)
    hm, hv = hypergeom.stats(700, int(d1.sum()), 100, moments="mv")
    ok7 = abs(kk.mean() - hm) < 0.1 and abs(kk.var() / hv - 1) < 0.05
    out.update({"band_flags_ok": bool(ok1), "weighted_auc_equals_sklearn": bool(ok2),
                "synthetic_label_inheritance_detected_and_separated": bool(ok3),
                "within_protein_mean_equals_sklearn": bool(ok4), "haldane_formula_ok": bool(ok5),
                "mixture_reference_limits_ok": bool(ok6), "protein_draws_reduce_to_hypergeometric": bool(ok7),
                "synthetic_control_row": {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r0.items()},
                "synthetic_added_signal_row": {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r1.items()}})
    out["all_passed"] = bool(ok1 and ok2 and ok3 and ok4 and ok5 and ok6 and ok7)
    return out


# ----------------------------------------------------------------------------- PART 3 runner (blocked until deposit)
def load_threshold_free_inputs(inp, oof_dir, fasta_by_cohort):
    data = {}
    for c, fa in fasta_by_cohort.items():
        seqs = read_fasta_gz(inp.use(fa))
        arms = {a: inp.read_csv(os.path.join(oof_dir, c, f"{c}_{a}_oof.csv")) for a in ARMS}
        data[c] = (seqs, arms)
    return data


def run_threshold_free(data, boot=BOOT, seed=SEED):
    """data[cohort] = (sequences by accession, {arm: oof rows}). Returns the threshold-free tables on all cysteines
    with a sequence and, revision round 2 (post hoc), on the cysteines that can be the only cysteine of a tryptic
    peptide (symmetric single-cysteine filter: both label classes restricted to coverable cysteines)."""
    tf_rows, tf_diffs, tfc_rows, tfc_diffs, diag = [], [], [], [], {}
    for c, (seqs, arms) in data.items():
        base = arms[ARMS[0]]
        acc = np.asarray([r["accession"] for r in base])
        pos = np.asarray([int(r["position"]) for r in base])
        keep = np.asarray([x in seqs and 0 < q <= len(seqs[x]) for x, q in zip(acc, pos)])
        y = np.asarray([int(r["label"]) for r in base])
        dist = np.asarray([band_flags(seqs[x], q, TRYPSIN)[0] if k else -1 for x, q, k in zip(acc, pos, keep)])
        ncov = np.asarray([second_cys_trypsinP(seqs[x], q - 1) if k else -1 for x, q, k in zip(acc, pos, keep)])
        by_arm = {}
        for a in ARMS:
            rr = arms[a]
            assert [r["accession"] for r in rr] == list(acc) and [int(r["position"]) for r in rr] == list(pos), (c, a)
            by_arm[a] = {m: np.asarray([float(r[m]) for r in rr])[keep] for m in MEMBERS}
        # recorded, not asserted: the labels are built from single-cysteine peptides (j5), so this should be 0;
        # human was checked only within the 15-residue windows before deposit
        diag[c] = {"rows": int(len(y)), "rows_with_sequence": int(keep.sum()),
                   "labelled_not_coverable": int(((ncov == 1) & (y == 1)).sum()),
                   "unlabelled_not_coverable": int(((ncov == 1) & (y == 0)).sum()),
                   "rows_coverable": int((ncov == 0).sum())}
        rows, diffs = threshold_free_table(c, acc[keep], y[keep], dist[keep], by_arm, MEMBERS, boot=boot, seed=seed)
        tf_rows += rows
        tf_diffs += diffs
        cov = ncov[keep] == 0
        by_arm_c = {a: {m: v[cov] for m, v in d.items()} for a, d in by_arm.items()}
        rows_c, diffs_c = threshold_free_table(c + "|coverable_only", acc[keep][cov], y[keep][cov], dist[keep][cov],
                                               by_arm_c, MEMBERS, boot=boot, seed=seed)
        tfc_rows += rows_c
        tfc_diffs += diffs_c
    return tf_rows, tf_diffs, tfc_rows, tfc_diffs, diag


# ----------------------------------------------------------------------------- main
def read_fasta_gz(path):
    out, acc, buf = {}, None, []
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                if acc:
                    out[acc] = "".join(buf)
                m = re.match(r">\w+\|([^|]+)\|", line)
                acc = m.group(1) if m else line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if acc:
        out[acc] = "".join(buf)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--oof-dir", default=REPO + "/external/public_refit_out_2026-09-21")
    ap.add_argument("--hsa", default=REPO + "/external/proteomes/hsa.fasta.gz")
    ap.add_argument("--osa", default=REPO + "/external/proteomes/osa.fasta.gz")
    args = ap.parse_args()
    inp = Inputs()
    outputs = []

    # ------------------------------------------------------------------ PART 1 inventory
    man = {r["path"]: r for r in inp.read_csv(REPO + "/results/artifact_manifest.csv")}
    need = []
    for c in COHORTS:
        for a in ARMS:
            need.append((f"external/public_refit_out_2026-09-21/{c}/{c}_{a}_oof.csv", "per-cysteine out-of-fold scores (blend, lgb_rank, lgb_bin, fusion_mlp)"))
            need.append((f"external/public_refit_out_2026-09-21/{c}/{c}_{a}_folds.json", "per-fold member selections and blend weights"))
        need.append((f"external/public_refit_out_2026-09-21/{c}/{c}_manifest.json", "run manifest (script and input sha256, host, threads)"))
        need.append((f"inputs/public_refit_2026-09-21/{c}_refit_inputs.npz", "feature matrix, labels, folds, embeddings used by the refit"))
    need += [("external/proteomes/hsa.fasta.gz", "human sequences for distal flags"),
             ("external/proteomes/osa.fasta.gz", "rice sequences for distal flags"),
             ("results/public_refit_ablation_lgb_rank_2026-09-24.csv", "ranking-member reading printed in Table 2 / Fig. 7b-c (Source Data source_table)"),
             ("<script that produced public_refit_ablation_lgb_rank_2026-09-24.csv>", "re-read of the lgb_rank column by the registered analysis code")]
    inv = []
    for rel, role in need:
        p = os.path.join(REPO, rel)
        m = man.get(rel, {})
        inv.append({"path_in_internal_tree": rel, "role": role, "present_on_this_machine": int(os.path.exists(p)),
                    "bytes_per_manifest": m.get("bytes", ""), "sha256_per_manifest": m.get("sha256", ""),
                    "manifest_entry": "yes" if m else "no",
                    "public_release_decision": "denied (R2 deny-prefix)" if rel.startswith(("external/", "inputs/")) else "not in release"})
    outputs.append(write_csv(OUT + "/ablation_missing_files.csv", inv))
    oof_present = all(r["present_on_this_machine"] for r in inv if r["path_in_internal_tree"].endswith("_oof.csv"))

    # ------------------------------------------------------------------ PART 2 stored-value analysis
    note11 = inp.read_text(MCP + "/supplemental/Supplemental_Note_11_public_feature_ablation.md")
    tabs = parse_note11(note11)
    stored = inp.read_csv(REPO + "/results/public_refit_ablation_2026-09-21.csv")
    # blend table in Note 11 must equal the stored registered CSV
    blend_ok = all(
        any(s["cohort"] == r["cohort"] and s["arm"] == r["arm"] and s["log2_or"] == r["log2_or"] and s["top100_distal_rate"] == r["top100_distal_rate"]
            and s["rate_rest"] == r["rate_rest"] for s in stored) for r in tabs["blend"])
    tex = inp.read_text(MCP + "/01_manuscript.tex")
    t2 = re.findall(r"^(human|rice) & (.+?) & ([\d,]+) & \\textbf\{([-\d.]+) \[([-\d.]+), ([-\d.]+)\]\} & ([\d.]+) \\\\", tex, re.M)
    armmap = {"chemistry set": "control", "minus win7/win10 basic composition": "drop_w710", "minus all twenty basic-composition columns": "drop_all_basic"}
    table2 = [{"cohort": ("human_PXD044043" if c == "human" else "rice_PXD072089"), "arm": armmap[a], "log2_or": v, "ci_low": lo, "ci_high": hi, "wp_auc": w}
              for c, a, _, v, lo, hi, w in t2]
    t2_ok = all(any(r["cohort"] == t["cohort"] and r["arm"] == t["arm"] and float(r["log2_or"]) == float(t["log2_or"])
                    and float(r["ci_low"]) == float(t["ci_low"]) and float(r["within_protein_auc"]) == float(t["wp_auc"])
                    for r in tabs["lgb_rank"]) for t in table2)

    note9 = inp.read_text(MCP + "/supplemental/Supplemental_Note_9_public_cleavage_geometry_self_audit.md")
    strat = {}
    for line in note9.splitlines():
        m = re.match(r"\| (human_PXD044043|rice_PXD072089) \| (v2r_chem_trained_on_rice|v2h_chem_trained_on_human) \| distal_6_12 \| (all|label_positive|label_negative) \| (\d+) \| ([\d.]+) \|", line)
        if m:
            strat[(m.group(1), m.group(3))] = (int(m.group(4)), float(m.group(5)))

    recon, gran, bench = [], [], []
    summary = {}
    for c in COHORTS:
        rows_c = [dict(r, member=mem) for mem in ("lgb_rank", "blend") for r in tabs[mem] if r["cohort"] == c]
        N = int(rows_c[0]["n"])
        lo, hi = feasible_D(rows_c, N)
        n_all, r_all = strat[(c, "all")]
        n_pos, r_pos = strat[(c, "label_positive")]
        n_neg, r_neg = strat[(c, "label_negative")]
        d_pos = int(round(r_pos * n_pos))
        d_neg_lo, d_neg_hi = math.ceil((r_neg - 0.00005) * n_neg), math.floor((r_neg + 0.00005) * n_neg)
        cands = [D for D in range(lo, hi + 1) if d_neg_lo <= D - d_pos <= d_neg_hi]
        assert len(cands) == 1, (c, lo, hi, cands)
        D = cands[0]
        summary[c] = {"N": N, "D_distal_total": D, "distal_prevalence": D / N, "feasible_D_from_arms": [lo, hi],
                      "N_label_positive": n_pos, "D_label_positive": d_pos, "N_label_negative": n_neg, "D_label_negative": D - d_pos,
                      "top_k_share_of_cohort": TOP_K / N}
        for r in rows_c:
            k = int(round(float(r["top100_distal_rate"]) * 100))
            val = top100_lor(k, D, N)
            recon.append({"cohort": c, "member": r["member"], "arm": r["arm"], "k_distal_in_top100": k,
                          "D_distal_total": D, "N": N, "stored_log2_or": r["log2_or"], "recomputed_log2_or": val,
                          "reproduces_at_4dp": int(abs(round(val, 4) - float(r["log2_or"])) < 1e-9),
                          "stored_ci": f"[{r['ci_low']}, {r['ci_high']}]", "within_protein_auc": r["within_protein_auc"],
                          "log2_or_if_one_fewer": top100_lor(k - 1, D, N), "log2_or_if_one_more": top100_lor(min(k + 1, 100), D, N)})
        for k in range(70, 101):
            gran.append({"cohort": c, "k_distal_in_top100": k, "log2_or": top100_lor(k, D, N),
                         "step_from_k_minus_1": (top100_lor(k, D, N) - top100_lor(k - 1, D, N))})
        # references
        for ref, M, nn in (("no association: 100 cysteines at random", N, D),
                           ("label inheritance: 100 label-positive cysteines at random", n_pos, d_pos)):
            q025, q50, q975 = hyper_q(M, nn, TOP_K, 0.025), hyper_q(M, nn, TOP_K, 0.5), hyper_q(M, nn, TOP_K, 0.975)
            bench.append({"cohort": c, "reference": ref, "expected_k": TOP_K * nn / M, "k_q025": q025, "k_median": q50, "k_q975": q975,
                          "log2_or_at_expected_k": top100_lor(TOP_K * nn / M, D, N),
                          "log2_or_at_q025": top100_lor(q025, D, N), "log2_or_at_q975": top100_lor(q975, D, N),
                          "ceiling_log2_or_k100": top100_lor(100, D, N)})
        for r in rows_c:
            k = int(round(float(r["top100_distal_rate"]) * 100))
            bench.append({"cohort": c, "reference": f"observed {r['member']} {r['arm']}", "expected_k": k, "k_q025": "", "k_median": "", "k_q975": "",
                          "log2_or_at_expected_k": float(r["log2_or"]),
                          "P_k_or_more_under_label_inheritance": float(hypergeom.sf(k - 1, n_pos, d_pos, TOP_K)),
                          "P_k_or_more_under_no_association": float(hypergeom.sf(k - 1, N, D, TOP_K))})
    outputs.append(write_csv(OUT + "/ablation_top100_reconstruction.csv", recon))
    outputs.append(write_csv(OUT + "/ablation_top100_granularity.csv", gran))
    outputs.append(write_csv(OUT + "/ablation_top100_references.csv", bench))

    # ------------------------------------------------------------------ PART 2b every arm against every reference
    avr = arm_vs_references(summary, tabs)
    for c in COHORTS:   # the f = 1.0 mixture must reproduce the round-1 hypergeometric label-inheritance range
        b = [x for x in bench if x["cohort"] == c and x["reference"].startswith("label inheritance")][0]
        a = [x for x in avr if x["cohort"] == c and x["f_labelled_in_top100"] == 1.0][0]
        assert (a["ref_k_q025"], a["ref_k_q975"]) == (b["k_q025"], b["k_q975"]), (c, a, b)
    outputs.append(write_csv(OUT + "/ablation_top100_arm_vs_reference.csv", avr))
    position_counts = {ref: count_positions(avr, ref) for ref in
                       ("no association", "label inheritance, 100 of 100 labelled", "label inheritance, 90 of 100 labelled",
                        "label inheritance, 80 of 100 labelled", "label inheritance, 50 of 100 labelled")}

    # ------------------------------------------------------------------ PART 2c protein clustering (rice only)
    rice_clu = rice_cluster_sensitivity(inp, summary, tabs)
    outputs.append(write_json(OUT + "/ablation_top100_rice_cluster_sensitivity.json", rice_clu))

    # ------------------------------------------------------------------ PART 3 threshold-free (blocked unless files exist)
    st = selftest()
    status = {"label": "POST HOC revision analysis (2026-09-30); not registered",
              "oof_files_present": bool(oof_present), "selftest": st,
              "stored_value_checks": {"note11_blend_table_equals_stored_registered_csv": bool(blend_ok),
                                      "table2_equals_note11_ranking_member_rows": bool(t2_ok),
                                      "table2_rows_parsed": len(table2)},
              "cohorts": summary,
              "arm_positions_vs_references (12 primary readings; revision after verification)": position_counts,
              "rice_protein_clustering_sensitivity": {k: v for k, v in rice_clu.items() if k != "labelled_per_protein_distribution"}}
    if oof_present and os.path.exists(args.hsa) and os.path.exists(args.osa):
        data = load_threshold_free_inputs(inp, args.oof_dir, {"human_PXD044043": args.hsa, "rice_PXD072089": args.osa})
        tf_rows, tf_diffs, tfc_rows, tfc_diffs, diag = run_threshold_free(data)
        outputs.append(write_csv(OUT + "/ablation_threshold_free.csv", tf_rows))
        outputs.append(write_csv(OUT + "/ablation_threshold_free_paired_diffs.csv", tf_diffs))
        outputs.append(write_csv(OUT + "/ablation_threshold_free_coverable_only.csv", tfc_rows))
        outputs.append(write_csv(OUT + "/ablation_threshold_free_paired_diffs_coverable_only.csv", tfc_diffs))
        status["threshold_free"] = "computed (all cysteines, and coverable cysteines only)"
        status["threshold_free_diagnostics"] = diag
    else:
        status["threshold_free"] = ("BLOCKED: the per-cysteine refit scores are not on this machine; see ablation_missing_files.csv. "
                                    "Run: python j3_ablation_statistic.py --oof-dir <external/public_refit_out_2026-09-21> "
                                    "--hsa <external/proteomes/hsa.fasta.gz> --osa <external/proteomes/osa.fasta.gz>. "
                                    "Revision round 2: the run also writes the symmetric-filter variant (both label classes "
                                    "restricted to cysteines that can be the only cysteine of a tryptic peptide).")
    outputs.append(write_json(OUT + "/ablation_threshold_free_status.json", status))
    write_fragment(os.path.basename(__file__), inp, outputs,
                   {"bootstrap_seed": SEED, "bootstrap_replicates": BOOT, "selftest_seed": SELFTEST_SEED,
                    "rice_cluster_reference_seed": CLUSTER_SEED, "rice_cluster_reference_draws": CLUSTER_DRAWS},
                   {"oof_present": bool(oof_present)})
    print(json.dumps({k: v for k, v in status.items() if k != "selftest"}, indent=1, ensure_ascii=False)[:3000])
    print("selftest all_passed:", st["all_passed"])


if __name__ == "__main__":
    main()
