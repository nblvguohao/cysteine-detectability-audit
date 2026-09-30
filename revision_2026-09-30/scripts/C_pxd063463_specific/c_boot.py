"""Protein-clustered bootstrap helpers for item C_pxd063463_specific (POST HOC revision analysis, 2026-09-30).

All intervals are percentile intervals over REPS replicates of protein multiplicities drawn by Cys-Audit's own
stats.multiplicities (the same generator the tool uses), so single-arm shares computed here with the tool's seed are
identical to the tool's positive_detected_overlap test on the same rows.

Two-arm contrasts use a PAIRED protein bootstrap: the resampling unit is a protein accession in the union of the two
arms, and one draw of multiplicities is applied to both arms, because the four digests were cut from one enrichment
and the same protein contributes to both arms.
"""
from __future__ import annotations

import numpy as np

import c_common as C  # noqa: F401  (puts the Cys-Audit source on sys.path)
from cys_audit import stats


def pct(vals, level=0.95):
    return stats.percentile_interval(vals, level)


def share_boot(proteins, success, reps, seed, level=0.95):
    """Clustered share (the tool's boot_proportion)."""
    success = np.asarray(success, dtype=bool)
    if success.size == 0:
        return np.nan, np.nan, np.nan
    codes, labels = stats.cluster_index(list(proteins))
    pt, (lo, hi), _ = stats.boot_proportion(codes, len(labels), success, reps, seed, level)
    return pt, lo, hi


def per_protein_bins(df, proteins_index, bin_col, n_bins, value_col):
    """K x n_bins matrices of counts (rows identified) and successes for one arm."""
    K = len(proteins_index)
    n = np.zeros((K, n_bins), dtype=np.int64)
    k = np.zeros((K, n_bins), dtype=np.int64)
    if len(df):
        ui = df.protein.map(proteins_index).values
        b = df[bin_col].values.astype(int)
        np.add.at(n, (ui, b), 1)
        np.add.at(k, (ui, b), df[value_col].values.astype(np.int64))
    return n, k


def paired_contrast(dfT, dfX, value_col, reps, seed, bin_col=None, n_bins=1, level=0.95, protein_equal=False):
    """Contrast of the share of `value_col` between arm X and arm T (reference, usually trypsin).

    Returns dict with cX, cT (raw), cT_std (T standardised to X's bin distribution when bin_col is given),
    diff = cX - cT_std, and paired bootstrap intervals. With protein_equal=True the share is the unweighted mean over
    proteins of each protein's share (bins ignored).
    """
    prots = sorted(set(dfT.protein) | set(dfX.protein))
    index = {p: i for i, p in enumerate(prots)}
    if bin_col is None:
        dfT = dfT.assign(_bin=0)
        dfX = dfX.assign(_bin=0)
        bin_col, n_bins = "_bin", 1
    nT, kT = per_protein_bins(dfT, index, bin_col, n_bins, value_col)
    nX, kX = per_protein_bins(dfX, index, bin_col, n_bins, value_col)

    def stat(NT, KT, NX, KX):
        # NT etc: (m x B) totals, or per-protein matrices when protein_equal
        with np.errstate(invalid="ignore", divide="ignore"):
            cX = KX.sum(axis=-1) / NX.sum(axis=-1)
            cT = KT.sum(axis=-1) / NT.sum(axis=-1)
            pX = NX / NX.sum(axis=-1, keepdims=True)
            rT = np.where(NT > 0, KT / np.where(NT > 0, NT, 1), 0.0)
            # standardise T to X's bin distribution; a bin with X rows but no T rows makes the replicate undefined
            cTs = (pX * rT).sum(axis=-1)
            bad = ((pX > 0) & (NT == 0)).any(axis=-1)
            cTs = np.where(bad, np.nan, cTs)
        return cX, cT, cTs

    if protein_equal:
        tT, sT = nT.sum(axis=1), kT.sum(axis=1)
        tX, sX = nX.sum(axis=1), kX.sum(axis=1)
        inT, inX = tT > 0, tX > 0
        with np.errstate(invalid="ignore", divide="ignore"):
            rT, rX = sT / np.where(inT, tT, 1), sX / np.where(inX, tX, 1)
        cX = float(rX[inX].mean())
        cT = float(rT[inT].mean())
        pt = {"cX": cX, "cT": cT, "cT_std": cT, "diff_raw": cX - cT, "diff_std": cX - cT}
        vals = {k: [] for k in pt}
        for w in stats.multiplicities(len(prots), reps, seed):
            wx = w * inX
            wt = w * inT
            vx = (wx @ np.where(inX, rX, 0.0)) / wx.sum(axis=1)
            vt = (wt @ np.where(inT, rT, 0.0)) / wt.sum(axis=1)
            vals["cX"].append(vx)
            vals["cT"].append(vt)
            vals["cT_std"].append(vt)
            vals["diff_raw"].append(vx - vt)
            vals["diff_std"].append(vx - vt)
    else:
        cX, cT, cTs = stat(nT.sum(0)[None, :], kT.sum(0)[None, :], nX.sum(0)[None, :], kX.sum(0)[None, :])
        pt = {"cX": float(cX[0]), "cT": float(cT[0]), "cT_std": float(cTs[0]), "diff_raw": float(cX[0] - cT[0]),
              "diff_std": float(cX[0] - cTs[0])}
        vals = {k: [] for k in pt}
        for w in stats.multiplicities(len(prots), reps, seed):
            NT, KT, NX, KX = w @ nT, w @ kT, w @ nX, w @ kX
            a, b, c = stat(NT, KT, NX, KX)
            vals["cX"].append(a)
            vals["cT"].append(b)
            vals["cT_std"].append(c)
            vals["diff_raw"].append(a - b)
            vals["diff_std"].append(a - c)
    out = dict(pt)
    for k, v in vals.items():
        v = np.concatenate(v)
        lo, hi = pct(v, level)
        out[f"{k}_ci_low"], out[f"{k}_ci_high"] = lo, hi
        out[f"{k}_n_undefined_reps"] = int((~np.isfinite(v)).sum())
    out["n_proteins_union"] = len(prots)
    out["n_T"], out["n_X"] = int(nT.sum()), int(nX.sum())
    out["k_T"], out["k_X"] = int(kT.sum()), int(kX.sum())
    out["n_proteins_T"], out["n_proteins_X"] = int((nT.sum(1) > 0).sum()), int((nX.sum(1) > 0).sum())
    return out


def paired_matched(df, colT, colX, reps, seed, level=0.95):
    """Same cysteines scored in two arms (one row per shared cysteine): share colX - share colT, paired over proteins."""
    prots = sorted(set(df.protein))
    index = {p: i for i, p in enumerate(prots)}
    ui = df.protein.map(index).values
    K = len(prots)
    n = np.bincount(ui, minlength=K).astype(np.int64)
    a = np.bincount(ui, weights=df[colT].values, minlength=K).astype(np.int64)
    c = np.bincount(ui, weights=df[colX].values, minlength=K).astype(np.int64)
    pt = {"cT": a.sum() / n.sum(), "cX": c.sum() / n.sum()}
    pt["diff"] = pt["cX"] - pt["cT"]
    vals = {"cT": [], "cX": [], "diff": []}
    for w in stats.multiplicities(K, reps, seed):
        N = (w @ n).astype(float)
        vt, vx = (w @ a) / N, (w @ c) / N
        vals["cT"].append(vt)
        vals["cX"].append(vx)
        vals["diff"].append(vx - vt)
    out = dict(pt)
    for k, v in vals.items():
        lo, hi = pct(np.concatenate(v), level)
        out[f"{k}_ci_low"], out[f"{k}_ci_high"] = lo, hi
    out["n_cys"] = int(n.sum())
    out["n_proteins"] = K
    out["discordant_T1_X0"] = int(((df[colT] == 1) & (df[colX] == 0)).sum())
    out["discordant_T0_X1"] = int(((df[colT] == 0) & (df[colX] == 1)).sum())
    out["concordant_11"] = int(((df[colT] == 1) & (df[colX] == 1)).sum())
    out["concordant_00"] = int(((df[colT] == 0) & (df[colX] == 0)).sum())
    return out
