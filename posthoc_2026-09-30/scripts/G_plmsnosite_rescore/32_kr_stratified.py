"""G_plmsnosite_rescore, step 3c (system python): is the pLMSNOSite score's association with peptide
detectability among negatives carried by the local lysine/arginine count?

POST HOC revision analysis (2026-09-30), added after adversarial verification round 1, which showed
that the claim "pLMSNOSite has partly learned detectability" rests on an association that largely
runs through pep_cleavage_sites_within_20 (the number of K and R within 20 residues of the cysteine,
the DIG25 column that VIS10 omits because it could also encode a chemical acid-base motif).
Not registered, not pre-specified; descriptive.

Covariate: kr = pep_cleavage_sites_within_20 (count of K/R in sequence[i-20 : i+21], proline not
considered; repo/scripts/run_cross_protease_detectability_probe.py, features()). Strata = each observed
value of kr. Indicators (DIG25 block): pep_detectable_any_missed_cleavage (peptide detectable with up to
two missed cleavages) and pep_detectable_both (fully cleaved peptide of 7-30 residues and 700-3,500 Da).

Statistics (one pass of 5,000 protein-clustered resamples shared by all statistics; seed 20260930;
all 278 test proteins resampled, as in 30_compare.py, so replicates coincide with the main analysis):
 N1  among negatives: Spearman(score, kr) for pLMSNOSite, its two arms, DIG25, VIS10; Spearman(indicator, kr)
 N2  among negatives: AUROC of each score for each indicator, crude and K/R-stratified (pair-weighted mean of
     within-stratum AUROCs = probability that a random within-stratum (indicator 1, indicator 0) pair is ordered
     correctly; strata lacking either indicator value do not contribute)
 N3  among negatives: Spearman(pLMSNOSite, DIG25) and Spearman(pLMSNOSite, VIS10), crude and pooled within
     K/R strata (within-stratum weighted ranks, centred within stratum, pooled weighted Pearson)
 N4  among negatives, within each level of the <=2-missed-cleavage indicator: Spearman(score, kr)
 L1  all sites: label AUROC of each score, crude and K/R-stratified; paired differences pLMSNOSite - VIS10 /
     - DIG25 and recovery ratios (AUROC_det - 0.5)/(AUROC_pLMSNOSite - 0.5), crude and K/R-stratified
Outputs: kr_stratified_summary.csv, kr_strata_detail.csv (point estimates per stratum), kr_stratified.json
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "4")

import json
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

SEED, REPS = 20260930, 5000
SCORES = ("pLMSNOSite", "Embedding_arm", "ProtT5_arm", "DIG25", "VIS10")


class WRank:
    """Weighted average ranks of a fixed vector (rows duplicated by integer weights); as in 30_compare.py."""

    def __init__(self, x):
        x = np.asarray(x, dtype=float)
        self.order = np.argsort(x, kind="mergesort")
        _, self.inv = np.unique(x[self.order], return_inverse=True)
        self.n = len(x)

    def ranks(self, w):
        wg = np.bincount(self.inv, weights=w[self.order])
        before = np.concatenate(([0.0], np.cumsum(wg)[:-1]))
        r = np.empty(self.n)
        r[self.order] = (before + (wg + 1.0) / 2.0)[self.inv]
        return r


def wpearson(a, b, w):
    W = w.sum()
    if W <= 1:
        return float("nan")
    ma, mb = (w * a).sum() / W, (w * b).sum() / W
    va, vb = (w * (a - ma) ** 2).sum(), (w * (b - mb) ** 2).sum()
    if va <= 0 or vb <= 0:
        return float("nan")
    return float((w * (a - ma) * (b - mb)).sum() / np.sqrt(va * vb))


class StratAUC:
    """Pair-weighted K/R-stratified AUROC of `score` for binary `lab` over rows `idx`."""

    def __init__(self, score, lab, kr, idx):
        self.parts = []
        for v in np.unique(kr[idx]):
            j = idx[kr[idx] == v]
            if lab[j].min() == lab[j].max():
                continue
            self.parts.append((j, lab[j].astype(float), C.RankedScore(score[j], lab[j])))

    def value(self, w):
        num = den = 0.0
        for j, lj, rs in self.parts:
            wj = w[j]
            pairs = (wj * lj).sum() * (wj * (1.0 - lj)).sum()
            if pairs > 0:
                num += rs.auroc(wj) * pairs
                den += pairs
        return num / den if den > 0 else float("nan")


class Spear:
    """Weighted Spearman of a and b over rows `idx` (crude)."""

    def __init__(self, a, b, idx):
        self.idx, self.ra, self.rb = idx, WRank(a[idx]), WRank(b[idx])

    def value(self, w):
        ww = w[self.idx]
        return wpearson(self.ra.ranks(ww), self.rb.ranks(ww), ww)


class SpearWithin:
    """Pooled within-stratum Spearman: ranks within stratum, centred within stratum, pooled Pearson."""

    def __init__(self, a, b, kr, idx):
        self.parts = []
        for v in np.unique(kr[idx]):
            j = idx[kr[idx] == v]
            if len(j) >= 2:
                self.parts.append((j, WRank(a[j]), WRank(b[j])))

    def value(self, w):
        sab = saa = sbb = 0.0
        for j, ra_, rb_ in self.parts:
            wj = w[j]
            W = wj.sum()
            if W < 2:
                continue
            ra, rb = ra_.ranks(wj), rb_.ranks(wj)
            ra = ra - (wj * ra).sum() / W
            rb = rb - (wj * rb).sum() / W
            sab += (wj * ra * rb).sum()
            saa += (wj * ra * ra).sum()
            sbb += (wj * rb * rb).sum()
        return sab / np.sqrt(saa * sbb) if saa > 0 and sbb > 0 else float("nan")


def main():
    t0 = time.time()
    feats = pd.read_csv(C.OUT / "detectability_test_features.csv")
    plm = pd.read_csv(C.OUT / "plmsnosite_test_scores.csv")
    det = pd.read_csv(C.OUT / "detectability_test_scores.csv")
    for d in (feats, det):
        assert (d.UniProt.values == plm.UniProt.values).all() and (d.Position.values == plm.Position.values).all()
    manifest = {f"results/G_plmsnosite_rescore/{f}": C.sha256_file(C.OUT / f)
                for f in ("detectability_test_features.csv", "plmsnosite_test_scores.csv", "detectability_test_scores.csv")}
    y = plm.Target.to_numpy().astype(int)
    kr = feats.pep_cleavage_sites_within_20.to_numpy()
    assert np.allclose(kr, np.round(kr))
    kr = np.round(kr).astype(int)
    IND = {"le2_missed_cleavages": feats.pep_detectable_any_missed_cleavage.to_numpy().astype(int),
           "fully_cleaved": feats.pep_detectable_both.to_numpy().astype(int)}
    S = {"pLMSNOSite": plm.plmsnosite_prob.to_numpy(), "Embedding_arm": plm.embedding_base_prob.to_numpy(),
         "ProtT5_arm": plm.prott5_base_prob.to_numpy(), "DIG25": det.dig25_score.to_numpy(),
         "VIS10": det.vis10_score.to_numpy()}
    codes, labels = C.protein_codes(plm.UniProt.tolist())
    K = len(labels)
    allrows = np.arange(len(y))
    neg = np.flatnonzero(y == 0)
    krf = kr.astype(float)

    # ---------------------------------------------------------------- statistic objects
    stat = {}
    for k in SCORES:
        stat[f"N1|spearman_score_kr|negatives|{k}"] = Spear(S[k], krf, neg)
    for i, z in IND.items():
        stat[f"N1|spearman_indicator_kr|negatives|{i}"] = Spear(z.astype(float), krf, neg)
        for k in SCORES:
            stat[f"N2|auroc_indicator_crude|negatives|{i}|{k}"] = C.RankedScore(S[k][neg], z[neg])
            stat[f"N2|auroc_indicator_KRstrat|negatives|{i}|{k}"] = StratAUC(S[k], z, kr, neg)
    for b in ("DIG25", "VIS10"):
        stat[f"N3|spearman_pLMSNOSite_{b}_crude|negatives"] = Spear(S["pLMSNOSite"], S[b], neg)
        stat[f"N3|spearman_pLMSNOSite_{b}_withinKR|negatives"] = SpearWithin(S["pLMSNOSite"], S[b], kr, neg)
    z2 = IND["le2_missed_cleavages"]
    for lev in (1, 0):
        idx = neg[z2[neg] == lev]
        for k in ("pLMSNOSite", "Embedding_arm", "ProtT5_arm"):
            stat[f"N4|spearman_score_kr|negatives_indicator{lev}|{k}"] = Spear(S[k], krf, idx)
    for k in SCORES:
        stat[f"L1|label_auroc_crude|all|{k}"] = C.RankedScore(S[k], y)
        stat[f"L1|label_auroc_KRstrat|all|{k}"] = StratAUC(S[k], y, kr, allrows)

    def evaluate(w):
        out = {}
        wneg = w[neg]
        for name, obj in stat.items():
            if name.startswith("N2|auroc_indicator_crude"):
                out[name] = obj.auroc(wneg)
            elif name.startswith("L1|label_auroc_crude"):
                out[name] = obj.auroc(w)
            else:
                out[name] = obj.value(w)
        return out

    # sanity: crude objects equal sklearn / scipy on the full data; weighted objects equal expanded rows
    ones = np.ones(len(y))
    pt = evaluate(ones)
    for k in SCORES:
        assert abs(pt[f"L1|label_auroc_crude|all|{k}"] - roc_auc_score(y, S[k])) < 1e-12
        rho = stats.spearmanr(S[k][neg], kr[neg])[0]
        assert abs(pt[f"N1|spearman_score_kr|negatives|{k}"] - rho) < 1e-10, (k, rho)
        for i, z in IND.items():
            assert abs(pt[f"N2|auroc_indicator_crude|negatives|{i}|{k}"] - roc_auc_score(z[neg], S[k][neg])) < 1e-12
            # stratified AUROC equals the pair-weighted mean of sklearn within-stratum AUROCs
            num = den = 0.0
            for v in np.unique(kr[neg]):
                j = neg[kr[neg] == v]
                if z[j].min() == z[j].max():
                    continue
                pairs = z[j].sum() * (1 - z[j]).sum()
                num += roc_auc_score(z[j], S[k][j]) * pairs
                den += pairs
            assert abs(pt[f"N2|auroc_indicator_KRstrat|negatives|{i}|{k}"] - num / den) < 1e-12
    mult = next(C.multiplicities(K, 1, 7))
    wv = mult[codes].astype(float)
    ev = evaluate(wv)
    ex = np.repeat(allrows, wv.astype(int))
    exn = ex[y[ex] == 0]
    assert abs(ev["N1|spearman_score_kr|negatives|pLMSNOSite"] - stats.spearmanr(S["pLMSNOSite"][exn], kr[exn])[0]) < 1e-10
    assert abs(ev["L1|label_auroc_crude|all|DIG25"] - roc_auc_score(y[ex], S["DIG25"][ex])) < 1e-12
    num = den = 0.0
    for v in np.unique(kr[exn]):
        j = exn[kr[exn] == v]
        zz = IND["le2_missed_cleavages"][j]
        if zz.min() == zz.max():
            continue
        pairs = zz.sum() * (1 - zz).sum()
        num += roc_auc_score(zz, S["pLMSNOSite"][j]) * pairs
        den += pairs
    assert abs(ev["N2|auroc_indicator_KRstrat|negatives|le2_missed_cleavages|pLMSNOSite"] - num / den) < 1e-12

    # ---------------------------------------------------------------- bootstrap
    reps = [evaluate(m[codes].astype(float)) for m in C.multiplicities(K, REPS, SEED)]
    R = pd.DataFrame(reps)
    # derived: paired differences and recovery ratios (crude and K/R-stratified label AUROC)
    derived = {}
    for kind in ("crude", "KRstrat"):
        P = R[f"L1|label_auroc_{kind}|all|pLMSNOSite"]
        p0 = pt[f"L1|label_auroc_{kind}|all|pLMSNOSite"]
        for b in ("VIS10", "DIG25"):
            B = R[f"L1|label_auroc_{kind}|all|{b}"]
            b0 = pt[f"L1|label_auroc_{kind}|all|{b}"]
            derived[f"L1|label_auroc_diff_pLMSNOSite_minus_{b}_{kind}|all"] = (p0 - b0, P - B)
            derived[f"L1|recovery_{b}_vs_pLMSNOSite_{kind}|all"] = ((b0 - 0.5) / (p0 - 0.5), (B - 0.5) / (P - 0.5))
    for i in IND:
        for k in ("pLMSNOSite", "Embedding_arm", "ProtT5_arm"):
            c_ = f"N2|auroc_indicator_crude|negatives|{i}|{k}"
            s_ = f"N2|auroc_indicator_KRstrat|negatives|{i}|{k}"
            derived[f"N2|auroc_indicator_change_KRstrat_minus_crude|negatives|{i}|{k}"] = (pt[s_] - pt[c_], R[s_] - R[c_])
    for b in ("DIG25", "VIS10"):
        c_ = f"N3|spearman_pLMSNOSite_{b}_crude|negatives"
        s_ = f"N3|spearman_pLMSNOSite_{b}_withinKR|negatives"
        derived[f"N3|spearman_pLMSNOSite_{b}_change_withinKR_minus_crude|negatives"] = (pt[s_] - pt[c_], R[s_] - R[c_])

    rows = []
    for name in stat:
        vals = R[name].to_numpy(float)
        lo, hi = C.pctl(vals)
        rows.append({"statistic": name, "value": pt[name], "ci_low": lo, "ci_high": hi,
                     "n_valid_replicates": int(np.isfinite(vals).sum())})
    for name, (v0, vals) in derived.items():
        vals = np.asarray(vals, float)
        lo, hi = C.pctl(vals)
        rows.append({"statistic": name, "value": v0, "ci_low": lo, "ci_high": hi,
                     "n_valid_replicates": int(np.isfinite(vals).sum()),
                     "share_replicates_le_0": C.boot_p_le(vals) if "diff" in name or "change" in name else None})
    summ = pd.DataFrame(rows)
    summ.insert(0, "block", summ.statistic.str.split("|").str[0])
    summ.to_csv(C.OUT / "kr_stratified_summary.csv", index=False, float_format="%.10g")

    # ---------------------------------------------------------------- per-stratum detail (point estimates)
    det_rows = []
    for v in np.unique(kr):
        j_all = allrows[kr == v]
        j_neg = neg[kr[neg] == v]
        row = {"kr_count": int(v), "n_sites": len(j_all), "n_pos": int(y[j_all].sum()), "n_neg": len(j_neg),
               "n_proteins": int(plm.UniProt.iloc[j_all].nunique())}
        for i, z in IND.items():
            row[f"neg_n_{i}_1"] = int(z[j_neg].sum())
            for k in SCORES:
                ok = len(j_neg) and z[j_neg].min() != z[j_neg].max()
                row[f"neg_auroc_{i}|{k}"] = roc_auc_score(z[j_neg], S[k][j_neg]) if ok else None
        for k in ("pLMSNOSite", "DIG25", "VIS10"):
            ok = y[j_all].min() != y[j_all].max()
            row[f"label_auroc|{k}"] = roc_auc_score(y[j_all], S[k][j_all]) if ok else None
        for k in ("pLMSNOSite", "Embedding_arm", "ProtT5_arm"):
            row[f"neg_mean_score|{k}"] = float(S[k][j_neg].mean()) if len(j_neg) else None
        det_rows.append(row)
    detail = pd.DataFrame(det_rows)
    detail.to_csv(C.OUT / "kr_strata_detail.csv", index=False, float_format="%.10g")

    info = {"label": "POST HOC revision analysis 2026-09-30 (G_plmsnosite_rescore), added after verification round 1; "
                     "not registered; descriptive",
            "bootstrap": {"replicates": REPS, "seed": SEED, "unit": "protein (UniProt accession), all 278 test proteins",
                          "draws": "numpy default_rng(seed); rng.integers(0, K, K) per replicate (shared with 30_compare.py)",
                          "interval": "95% percentile"},
            "covariate": "pep_cleavage_sites_within_20 = number of K and R in sequence[i-20:i+21] (41 residues incl. the "
                         "cysteine; proline not considered)",
            "n_strata_all": int(len(np.unique(kr))), "n_negatives": int(len(neg)),
            "inputs_sha256": manifest, "runtime_s": round(time.time() - t0, 1)}
    (C.OUT / "kr_stratified.json").write_text(json.dumps(info, indent=1), encoding="utf-8")
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 90)
    print(summ.round(4).to_string())
    print(detail.round(3).to_string())
    print("runtime", round(time.time() - t0, 1), "s")


if __name__ == "__main__":
    main()
