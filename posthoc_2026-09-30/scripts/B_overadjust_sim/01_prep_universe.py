# -*- coding: utf-8 -*-
"""B_overadjust_sim step 1 (POST HOC): build the cysteine universe with real detection.

Data: PXD063463 unenriched (global) arm, as tabulated by the internal repository's phase-3 script
(`phase3_candidate_peptides_2026-09-22.csv`: one row per in-silico tryptic Cys-containing
candidate peptide of identified proteins; empirical_detected = the exact (protein, start,
sequence) triple was identified in the global arm). Sequences: UniProt UP000000589 (mouse),
release 2026_03.

Steps
  1. Verify every candidate peptide against the FASTA at its stated start; drop every protein
     with any mismatching row or absent from the FASTA (match rates reported).
  2. A cysteine is detected if >= 1 candidate peptide containing it has empirical_detected == 1.
  3. Universe = every cysteine of every retained protein that carries >= 1 detected cysteine.
  4. VIS10 per cysteine with the repository's own `phase2_claim_cohorts.site_feature_matrix`
     (tryptic rule, imported verbatim) and the claim attributes with the repository's own
     attribute functions (imported verbatim).
  5. Descriptive diagnostics: how strongly each attribute is associated with detection and how
     well VIS10 predicts each attribute (the over-adjustment potential).
"""
from __future__ import annotations

import json
import os
import sys
import time

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common_b as cb  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.linear_model import LogisticRegression  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

import phase2_claim_cohorts as p2c  # noqa: E402  (repo, read-only)
import phase2b_claim_cohorts as p2b  # noqa: E402
import phase2d_claim_cohorts as p2d  # noqa: E402
from run_cross_protease_detectability_probe import FEATURE_NAMES  # noqa: E402

VIS10_NAMES = [FEATURE_NAMES[i] for i in p2c.VIS10_INDEX]

# attribute name -> (claim, builder, label). Builders are the repository's own functions.
ATTRIBUTES = {
    "a1_SFE006_KR": ("SFE-006", lambda s, k: p2c.flank_residue_flag(s, k, p2c.SFE006_KR_OFFSETS, {"K", "R"}),
                     "K or R at -10, -8..-6, -4, -2, +4..+8 (SFE-006)"),
    "a1b_SNO021_K": ("SNO-021", lambda s, k: p2c._flank_residue_flag(s, k, p2c.SNO021_K_OFFSETS, "K"),
                     "K at +10, -10, +6, +5, +1 (SNO-021)"),
    "a2_SNO016_KRH": ("SNO-016", lambda s, k: p2b.flank_count_flag(s, k, p2b.SNO016_RESIDUES),
                      ">=3 of the 10 flanking positions within +/-5 in {K,R,H} (SNO-016)"),
    "a2b_PERS009_AKRV": ("PERS-009", lambda s, k: p2b.flank_count_flag(s, k, p2b.PERS009_RESIDUES),
                         ">=3 of the 10 flanking positions within +/-5 in {A,K,R,V} (PERS-009)"),
    "a3_SNO006_DE3": ("SNO-006", lambda s, k: p2c.flank_residue_flag(s, k, p2d.SNO006_NEG_PLUS3, p2d.NEGATIVE_RESIDUES),
                      "D or E at +3 (acidic control; SNO-006)"),
    "a3b_SFE002_E": ("SFE-002", lambda s, k: p2c.flank_residue_flag(s, k, p2d.SFE002_GLU_OFFSETS, {"E"}),
                     "E at -4, -3, +1, +3, +4, +5 (acidic control; SFE-002)"),
}


def log2_or_table(y, a):
    n11 = float(((a == 1) & (y == 1)).sum()); n10 = float(((a == 1) & (y == 0)).sum())
    n01 = float(((a == 0) & (y == 1)).sum()); n00 = float(((a == 0) & (y == 0)).sum())
    return float(np.log2(((n11 + .5) * (n00 + .5)) / ((n10 + .5) * (n01 + .5))))


def main():
    t0 = time.time()
    cand = pd.read_csv(cb.CANDIDATES, dtype={"protein": str, "sequence": str})
    fasta = p2b.read_fasta_gz(cb.FASTA)  # verbatim repository reader: accession = 2nd |-field
    n_rows = len(cand)
    seqs = cand["protein"].map(fasta)
    present = seqs.notna()
    ok = np.zeros(n_rows, dtype=bool)
    starts = cand["start"].to_numpy(); ends = cand["end"].to_numpy(); peps = cand["sequence"].to_numpy()
    sv = seqs.to_numpy()
    for i in range(n_rows):
        s = sv[i]
        if isinstance(s, str):
            ok[i] = s[starts[i] - 1:ends[i]] == peps[i]
    cand["fasta_match"] = ok
    per_prot = cand.groupby("protein")["fasta_match"].agg(["size", "sum"])
    prot_all = per_prot.index.to_numpy()
    prot_missing = sorted(set(prot_all) - set(fasta))
    prot_mismatch = sorted(p for p, r in per_prot.iterrows() if p in fasta and r["sum"] < r["size"])
    keep_prot = sorted(p for p, r in per_prot.iterrows() if p in fasta and r["sum"] == r["size"])
    keep_set = set(keep_prot)
    ver = {
        "candidate_rows": int(n_rows),
        "rows_protein_in_fasta": int(present.sum()),
        "rows_sequence_match": int(ok.sum()),
        "row_match_rate": float(ok.mean()),
        "row_match_rate_among_present": float(ok[present.to_numpy()].mean()),
        "proteins": int(len(prot_all)),
        "proteins_absent_from_fasta": int(len(prot_missing)),
        "proteins_with_any_mismatch": int(len(prot_mismatch)),
        "proteins_retained_all_rows_match": int(len(keep_prot)),
        "protein_match_rate": float(len(keep_prot) / len(prot_all)),
        "detected_rows_total": int((cand["empirical_detected"] == 1).sum()),
        "detected_rows_in_retained_proteins": int(((cand["empirical_detected"] == 1)
                                                   & cand["protein"].isin(keep_set)).sum()),
    }
    print(json.dumps(ver, indent=1))

    # detected cysteines
    det_rows = cand[(cand["empirical_detected"] == 1) & cand["protein"].isin(keep_set)]
    detected = set()
    for acc, st, pep in zip(det_rows["protein"], det_rows["start"], det_rows["sequence"]):
        for i, ch in enumerate(pep):
            if ch == "C":
                detected.add((acc, int(st) + i))
    # candidate coverage (any generated candidate peptide contains the Cys) in retained proteins
    covered = set(); covered_7_30 = set()
    kept_rows = cand[cand["protein"].isin(keep_set)]
    for acc, st, pep, L in zip(kept_rows["protein"], kept_rows["start"], kept_rows["sequence"], kept_rows["length"]):
        for i, ch in enumerate(pep):
            if ch == "C":
                covered.add((acc, int(st) + i))
                if 7 <= L <= 30:
                    covered_7_30.add((acc, int(st) + i))
    abundance = cand.groupby("protein")["log_abundance"].first().to_dict()

    # Peptide-level empirical detection model (phase-3 M3 feature set without the three
    # protein-level columns), used ONLY to build the 'steep' sensitivity scenario in which
    # positives are a more detectability-selected subset of detected cysteines. Per cysteine:
    # logit of the highest predicted detection probability among candidates carrying it.
    m3_feats = ["length", "missed_cleavages", "mass", "gravy", "frac_basic", "frac_acidic",
                "dist_n_term", "dist_c_term", "n_cys"]
    Xc = StandardScaler().fit_transform(kept_rows[m3_feats].to_numpy(dtype=float))
    m3 = LogisticRegression(C=1.0, max_iter=2000, solver="lbfgs").fit(Xc, kept_rows["empirical_detected"].to_numpy())
    p_cand = m3.predict_proba(Xc)[:, 1]
    m3_auc = float(roc_auc_score(kept_rows["empirical_detected"].to_numpy(), p_cand))
    best_p = {}
    for acc, st, pep, pp in zip(kept_rows["protein"], kept_rows["start"], kept_rows["sequence"], p_cand):
        for i, ch in enumerate(pep):
            if ch == "C":
                k = (acc, int(st) + i)
                if pp > best_p.get(k, -1.0):
                    best_p[k] = float(pp)
    m3_info = {"features": m3_feats, "candidate_level_auc_in_sample": m3_auc,
               "coefficients_standardised": dict(zip(m3_feats, [float(c) for c in m3.coef_[0]]))}

    universe_prot = sorted({a for a, _ in detected})
    keys = []
    for acc in universe_prot:
        s = fasta[acc]
        for i, ch in enumerate(s):
            if ch == "C":
                keys.append((acc, i + 1))
    assert all(fasta[a][p - 1] == "C" for a, p in detected)
    seq_of = {a: fasta[a] for a in universe_prot}
    t1 = time.time()
    vis = p2c.site_feature_matrix(seq_of, keys, "VIS10")  # verbatim repository VIS10
    print("VIS10 computed for %d cysteines in %.1fs" % (len(keys), time.time() - t1))
    df = pd.DataFrame({"protein": [a for a, _ in keys], "position": [p for _, p in keys]})
    df["detected"] = [1 if k in detected else 0 for k in keys]
    df["in_candidate"] = [1 if k in covered else 0 for k in keys]
    df["in_candidate_7_30"] = [1 if k in covered_7_30 else 0 for k in keys]
    df["log_abundance"] = [abundance.get(a, np.nan) for a, _ in keys]
    bp = np.asarray([best_p.get(k, np.nan) for k in keys], dtype=float)
    bp = np.clip(bp, 1e-6, 1 - 1e-6)
    df["m3_best_logit"] = np.log(bp / (1 - bp))
    for j, name in enumerate(VIS10_NAMES):
        df[name] = vis[:, j]
    for name, (_, fn, _) in ATTRIBUTES.items():
        df[name] = fn(seq_of, keys)
    df.to_csv(cb.UNIVERSE, index=False, compression="gzip")

    # ---------------- descriptive diagnostics ----------------
    y_det = df["detected"].to_numpy()
    X = StandardScaler().fit_transform(df[VIS10_NAMES].to_numpy())
    groups = df["protein"].to_numpy()
    rng = np.random.default_rng(cb.MASTER_SEED)
    uprot = np.unique(groups)
    fold_of = dict(zip(uprot, rng.permutation(len(uprot)) % 5))
    fold = np.asarray([fold_of[g] for g in groups])

    def cv_auc(target):
        oof = np.zeros(len(target))
        for k in range(5):
            tr, te = fold != k, fold == k
            m = LogisticRegression(C=1.0, max_iter=1000, solver="lbfgs").fit(X[tr], target[tr])
            oof[te] = m.predict_proba(X[te])[:, 1]
        return float(roc_auc_score(target, oof))

    diag = {"vis10_predicts_detection_auc_cv5": cv_auc(y_det)}
    attr_rows = []
    for name, (claim, _, label) in ATTRIBUTES.items():
        a = df[name].to_numpy()
        attr_rows.append({
            "attribute": name, "claim": claim, "label": label,
            "prevalence_all": float(a.mean()),
            "prevalence_detected": float(a[y_det == 1].mean()),
            "prevalence_undetected": float(a[y_det == 0].mean()),
            "log2_or_attribute_vs_detection": log2_or_table(y_det, a),
            "vis10_predicts_attribute_auc_cv5": cv_auc(a),
        })
        print(attr_rows[-1])
    summary = {
        "item": cb.ITEM, "label": "POST HOC revision analysis (not registered, not pre-specified)",
        "verification": ver,
        "n_universe_proteins": len(universe_prot), "n_universe_cysteines": len(keys),
        "n_detected_cysteines": int(y_det.sum()),
        "detected_fraction": float(y_det.mean()),
        "n_cys_in_any_candidate": int(df["in_candidate"].sum()),
        "n_cys_in_candidate_7_30": int(df["in_candidate_7_30"].sum()),
        "detected_per_protein_mean": float(df.groupby("protein")["detected"].sum().mean()),
        "cys_per_protein_mean": float(df.groupby("protein").size().mean()),
        "vis10_names": VIS10_NAMES,
        "peptide_detection_model_for_steep_scenario": m3_info,
        "attributes": attr_rows, "diagnostics": diag,
        "vis10_variance_zero_columns": [n for n in VIS10_NAMES if df[n].std() == 0],
        "runtime_s": round(time.time() - t0, 1),
        "inputs_sha256": {"candidates": cb.sha256_of(cb.CANDIDATES), "fasta": cb.sha256_of(cb.FASTA)},
    }
    with open(cb.PREP_SUMMARY, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps({k: v for k, v in summary.items() if k not in ("attributes",)}, indent=1))


if __name__ == "__main__":
    main()
