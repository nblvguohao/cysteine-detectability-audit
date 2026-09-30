"""G_plmsnosite_rescore, step 3 (system python): compare the re-scored pLMSNOSite predictor with
the detectability-only baselines on the pLMSNOSite independent test set.

POST HOC revision analysis (2026-09-30) in response to review criticism; not registered,
not pre-specified. Seeds: bootstrap 20260930 (main), 20260924 (Monte Carlo stability check).

Inputs (written by steps 1-2): detectability_test_scores.csv, plmsnosite_test_scores.csv;
git blobs of the released data (sequence_train.csv, sequence_test.csv, protT5_test.csv).

Analyses
 T2  pLMSNOSite score-based AUROC / AUPRC with protein-clustered bootstrap 95% CIs; threshold-0.5
     sensitivity, specificity, their geometric and arithmetic means (the latter equals the AUROC that
     evaluate_model.py computes from thresholded predictions), MCC.
 T3  paired protein-clustered bootstrap of AUROC and AUPRC differences, pLMSNOSite - VIS10 and
     pLMSNOSite - DIG25 (all scores share every resample); DeLong (site-level) as a sensitivity.
 T4  does pLMSNOSite track detectability? (a) its AUROC within DIG25-score tertiles;
     (b) Spearman(pLMSNOSite, DIG25) within negatives and within positives; (c) descriptive
     logistic regression label ~ logit(pLMSNOSite) + logit(DIG25), protein-cluster-robust CIs;
     (d) AUROC of the unfitted rank average of pLMSNOSite and DIG25.
 T5  recovery ratios (AUROC_det - 0.5)/(AUROC_pLMSNOSite - 0.5) with the re-scored AUROC as the
     denominator, paired bootstrap CIs; the manuscript's quoted-anchor ratios recomputed alongside.
 R   robustness: data-quality exclusions, sequence-identity clusters, within-protein AUROC,
     DeLong, second bootstrap seed.
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "4")

import hashlib
import json
import platform
import sys
import time

import numpy as np
import pandas as pd
import scipy
import sklearn
import statsmodels
import statsmodels.api as sm
from scipy import stats
from sklearn.metrics import average_precision_score, roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

SEED = 20260930
SEED_ALT = 20260924
REPS = 5000
MODE = "per_replicate"


# ------------------------------------------------------------------ helpers
class WRank:
    """Weighted average ranks for a fixed vector, so Spearman on a bootstrap replicate
    (rows duplicated by integer weights) is computed without expanding rows."""

    def __init__(self, x):
        self.order = np.argsort(x, kind="mergesort")
        _, self.inv = np.unique(np.asarray(x)[self.order], return_inverse=True)
        self.n = len(x)

    def ranks(self, w):
        ws = w[self.order]
        wg = np.bincount(self.inv, weights=ws)
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


def within_protein_pairs(scores, y, codes, K):
    """Per-protein concordant-pair count (ties 0.5) and pair count."""
    conc, pairs = np.zeros(K), np.zeros(K)
    for k in range(K):
        idx = np.flatnonzero(codes == k)
        yp = y[idx].astype(bool)
        if yp.all() or (~yp).all():
            continue
        sp, sn = scores[idx][yp], scores[idx][~yp]
        d = sp[:, None] - sn[None, :]
        conc[k] = (d > 0).sum() + 0.5 * (d == 0).sum()
        pairs[k] = d.size
    return conc, pairs


def delong(y, score_list):
    """Fast DeLong (Sun & Xu 2014) for k correlated AUROCs; returns aucs and covariance."""
    y = np.asarray(y).astype(bool)
    m, n = y.sum(), (~y).sum()
    v10, v01, aucs = [], [], []
    for s in score_list:
        s = np.asarray(s, float)
        tz = stats.rankdata(s)
        tx = stats.rankdata(s[y])
        ty = stats.rankdata(s[~y])
        a = (tz[y] - tx) / n
        b = 1.0 - (tz[~y] - ty) / m
        v10.append(a)
        v01.append(b)
        aucs.append(a.mean())
    s10 = np.cov(np.vstack(v10)) if len(score_list) > 1 else np.var(v10[0], ddof=1)
    s01 = np.cov(np.vstack(v01)) if len(score_list) > 1 else np.var(v01[0], ddof=1)
    return np.array(aucs), np.atleast_2d(s10 / m + s01 / n)


def r4(x):
    return None if x is None or not np.isfinite(x) else float(x)


# ------------------------------------------------------------------ bootstrap engine
def run_bootstrap(df, scores, seed, reps=REPS, tertile_col=None, spearman_pairs=(), wp_scores=(),
                  thr_score=None, cluster_col="UniProt"):
    """One pass over `reps` protein-clustered resamples shared by every statistic.
    Returns dict name -> array of replicate values."""
    y = df["Target"].to_numpy().astype(int)
    codes, labels = C.protein_codes(df[cluster_col].tolist())
    K = len(labels)
    ranked = {k: C.RankedScore(v, y) for k, v in scores.items()}
    out = {f"auroc|{k}": [] for k in scores} | {f"auprc|{k}": [] for k in scores}
    out["prevalence"] = []
    # tertiles
    tert = None
    if tertile_col is not None:
        tert = df[tertile_col].to_numpy()
        for t in np.unique(tert):
            for k in scores:
                out[f"tert{t}|auroc|{k}"] = []
    # spearman
    wr = {}
    for a, b in spearman_pairs:
        for cls in (0, 1):
            sel = y == cls
            wr[(a, cls)] = (sel, WRank(scores[a][sel]))
            wr[(b, cls)] = (sel, WRank(scores[b][sel]))
            out[f"spearman|{a}~{b}|class{cls}"] = []
    # within-protein pooled AUROC
    wp = {k: within_protein_pairs(np.asarray(scores[k], float), y, codes, K) for k in wp_scores}
    for k in wp_scores:
        out[f"withinprot|{k}"] = []
    if thr_score is not None:
        call = np.asarray(scores[thr_score]) > 0.5
        for q in ("sn", "sp", "ba", "gmean", "mcc"):
            out[f"thr|{q}"] = []
    for mult in C.multiplicities(K, reps, seed, MODE):
        w = mult[codes].astype(float)
        out["prevalence"].append((w * y).sum() / w.sum())
        for k, rs in ranked.items():
            out[f"auroc|{k}"].append(rs.auroc(w))
            out[f"auprc|{k}"].append(rs.auprc(w))
        if tert is not None:
            for t in np.unique(tert):
                wt = w * (tert == t)
                for k, rs in ranked.items():
                    out[f"tert{t}|auroc|{k}"].append(rs.auroc(wt))
        for a, b in spearman_pairs:
            for cls in (0, 1):
                sel, ra = wr[(a, cls)]
                _, rb = wr[(b, cls)]
                ww = w[sel]
                out[f"spearman|{a}~{b}|class{cls}"].append(wpearson(ra.ranks(ww), rb.ranks(ww), ww))
        for k in wp_scores:
            conc, pairs = wp[k]
            out[f"withinprot|{k}"].append(float((mult * conc).sum() / (mult * pairs).sum()))
        if thr_score is not None:
            tp = (w * (y == 1) * call).sum()
            fn = (w * (y == 1) * ~call).sum()
            tn = (w * (y == 0) * ~call).sum()
            fp = (w * (y == 0) * call).sum()
            sn_, sp_ = tp / (tp + fn), tn / (tn + fp)
            den = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
            out["thr|sn"].append(sn_)
            out["thr|sp"].append(sp_)
            out["thr|ba"].append((sn_ + sp_) / 2)
            out["thr|gmean"].append(np.sqrt(sn_ * sp_))
            out["thr|mcc"].append((tp * tn - fp * fn) / den if den > 0 else np.nan)
    return {k: np.asarray(v, float) for k, v in out.items()}, K


def point_stats(df, scores, tertile_col=None, spearman_pairs=(), wp_scores=(), thr_score=None,
                cluster_col="UniProt"):
    y = df["Target"].to_numpy().astype(int)
    codes, labels = C.protein_codes(df[cluster_col].tolist())
    w = np.ones(len(y))
    res = {"prevalence": y.mean()}
    for k, v in scores.items():
        rs = C.RankedScore(v, y)
        res[f"auroc|{k}"] = rs.auroc()
        res[f"auprc|{k}"] = rs.auprc()
    if tertile_col is not None:
        tert = df[tertile_col].to_numpy()
        for t in np.unique(tert):
            for k, v in scores.items():
                rs = C.RankedScore(v, y)
                res[f"tert{t}|auroc|{k}"] = rs.auroc(w * (tert == t))
    for a, b in spearman_pairs:
        for cls in (0, 1):
            sel = y == cls
            res[f"spearman|{a}~{b}|class{cls}"] = float(stats.spearmanr(scores[a][sel], scores[b][sel])[0])
    for k in wp_scores:
        conc, pairs = within_protein_pairs(np.asarray(scores[k], float), y, codes, len(labels))
        res[f"withinprot|{k}"] = float(conc.sum() / pairs.sum())
    if thr_score is not None:
        call = np.asarray(scores[thr_score]) > 0.5
        tp = int(((y == 1) & call).sum()); fn = int(((y == 1) & ~call).sum())
        tn = int(((y == 0) & ~call).sum()); fp = int(((y == 0) & call).sum())
        sn_, sp_ = tp / (tp + fn), tn / (tn + fp)
        res.update({"thr|tp": tp, "thr|fn": fn, "thr|tn": tn, "thr|fp": fp, "thr|sn": sn_, "thr|sp": sp_,
                    "thr|ba": (sn_ + sp_) / 2, "thr|gmean": float(np.sqrt(sn_ * sp_)),
                    "thr|mcc": (tp * tn - fp * fn) / np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)),
                    "thr|ppv": tp / (tp + fp)})
    return res


def ci(arr):
    lo, hi = C.pctl(arr)
    return lo, hi


# ------------------------------------------------------------------ main
def main():
    t0 = time.time()
    manifest = {}
    det = pd.read_csv(C.OUT / "detectability_test_scores.csv")
    plm = pd.read_csv(C.OUT / "plmsnosite_test_scores.csv")
    for f in ("detectability_test_scores.csv", "plmsnosite_test_scores.csv"):
        manifest[f"results/G_plmsnosite_rescore/{f}"] = C.sha256_file(C.OUT / f)
    assert (det.UniProt.values == plm.UniProt.values).all() and (det.Position.values == plm.Position.values).all()
    assert (det.Target.values == plm.Target.values).all()
    train = C.read_blob_csv("data/train/sequence_train.csv", manifest)
    test = C.read_blob_csv("data/test/sequence_test.csv", manifest)
    pt5 = C.read_blob_csv("data/test/protT5_test.csv", manifest)
    assert (test.UniProt.values == det.UniProt.values).all() and (test.Position.values == det.Position.values).all()

    df = det[["row", "UniProt", "Position", "Target"]].copy()
    df["pLMSNOSite"] = plm["plmsnosite_prob"].to_numpy()
    df["VIS10"] = det["vis10_score"].to_numpy()
    df["DIG25"] = det["dig25_score"].to_numpy()
    df["ProtT5_arm"] = plm["prott5_base_prob"].to_numpy()
    df["Embedding_arm"] = plm["embedding_base_prob"].to_numpy()
    df["pLMSNOSite_call_0.5"] = (df["pLMSNOSite"] > 0.5).astype(float)
    n = len(df)
    rk = lambda v: stats.rankdata(v) / n  # noqa: E731
    df["rankavg_pLMSNOSite_DIG25"] = (rk(df["pLMSNOSite"]) + rk(df["DIG25"])) / 2
    df["rankavg_pLMSNOSite_VIS10"] = (rk(df["pLMSNOSite"]) + rk(df["VIS10"])) / 2
    q1, q2 = np.quantile(df["DIG25"], [1 / 3, 2 / 3])
    df["DIG25_tertile"] = np.where(df["DIG25"] <= q1, 1, np.where(df["DIG25"] <= q2, 2, 3))
    v1, v2 = np.quantile(df["VIS10"], [1 / 3, 2 / 3])
    df["VIS10_tertile"] = np.where(df["VIS10"] <= v1, 1, np.where(df["VIS10"] <= v2, 2, 3))

    # data-quality flags (for robustness)
    def win(seq, site, w=37):
        h = (w - 1) // 2
        s = "-" * h + seq + "-" * h
        site = site + h
        return s[site - 1 - h:site + h]
    train_w = set(win(s, p) for s, p in zip(train.sequences, train.Position))
    test_w = [win(s, p) for s, p in zip(test.sequences, test.Position)]
    df["flag_window_in_train"] = [w_ in train_w for w_ in test_w]
    df["flag_shared_id"] = df.UniProt.isin(set(train.UniProt)).to_numpy()
    df["flag_seq_in_train"] = test.sequences.isin(set(train.sequences)).to_numpy()
    df["flag_prott5_cropped"] = (pt5.New_Position.values != pt5.Position.values)
    df["flag_duplicate_site"] = df.duplicated(["UniProt", "Position"], keep="first").to_numpy()
    seq_cluster = {s: f"seqcluster_{i:04d}" for i, s in enumerate(sorted(set(test.sequences)))}
    df["seq_cluster"] = test.sequences.map(seq_cluster).to_numpy()

    main_scores = ["pLMSNOSite", "VIS10", "DIG25", "ProtT5_arm", "Embedding_arm", "pLMSNOSite_call_0.5",
                   "rankavg_pLMSNOSite_DIG25", "rankavg_pLMSNOSite_VIS10"]
    scores = {k: df[k].to_numpy(dtype=float) for k in main_scores}
    sp_pairs = [("pLMSNOSite", "DIG25"), ("pLMSNOSite", "VIS10")]
    wp_list = ["pLMSNOSite", "VIS10", "DIG25"]

    # sanity: weighted implementations equal sklearn on the full data
    y = df.Target.to_numpy()
    for k in ("pLMSNOSite", "VIS10", "DIG25"):
        assert abs(C.RankedScore(scores[k], y).auroc() - roc_auc_score(y, scores[k])) < 1e-12
        assert abs(C.RankedScore(scores[k], y).auprc() - average_precision_score(y, scores[k])) < 1e-12
    # sanity: weighted Spearman equals scipy on an explicitly expanded replicate
    codes_, labs_ = C.protein_codes(df.UniProt.tolist())
    mult = next(C.multiplicities(len(labs_), 1, 1))
    wv = mult[codes_].astype(float)
    sel = (y == 0)
    exp_idx = np.repeat(np.flatnonzero(sel), wv[sel].astype(int))
    rho_exp = stats.spearmanr(scores["pLMSNOSite"][exp_idx], scores["DIG25"][exp_idx])[0]
    rho_w = wpearson(WRank(scores["pLMSNOSite"][sel]).ranks(wv[sel]), WRank(scores["DIG25"][sel]).ranks(wv[sel]), wv[sel])
    assert abs(rho_exp - rho_w) < 1e-10, (rho_exp, rho_w)
    exp_all = np.repeat(np.arange(n), wv.astype(int))
    for k in ("pLMSNOSite", "DIG25"):
        assert abs(C.RankedScore(scores[k], y).auroc(wv) - roc_auc_score(y[exp_all], scores[k][exp_all])) < 1e-12
        assert abs(C.RankedScore(scores[k], y).auprc(wv) - average_precision_score(y[exp_all], scores[k][exp_all])) < 1e-12

    # ------------------------------------------------ main bootstrap (seed 20260930)
    pt = point_stats(df, scores, tertile_col="DIG25_tertile", spearman_pairs=sp_pairs, wp_scores=wp_list,
                     thr_score="pLMSNOSite")
    bs, K = run_bootstrap(df, scores, SEED, tertile_col="DIG25_tertile", spearman_pairs=sp_pairs,
                          wp_scores=wp_list, thr_score="pLMSNOSite")
    # VIS10-tertile sensitivity (smaller pass: AUROC only)
    pt_v = point_stats(df, {"pLMSNOSite": scores["pLMSNOSite"], "VIS10": scores["VIS10"]}, tertile_col="VIS10_tertile")
    bs_v, _ = run_bootstrap(df, {"pLMSNOSite": scores["pLMSNOSite"], "VIS10": scores["VIS10"]}, SEED,
                            tertile_col="VIS10_tertile")
    # second seed (Monte Carlo stability)
    bs_alt, _ = run_bootstrap(df, {k: scores[k] for k in ("pLMSNOSite", "VIS10", "DIG25")}, SEED_ALT)

    # save replicate draws for audit
    rep = pd.DataFrame({k: v for k, v in bs.items()})
    rep.insert(0, "replicate", np.arange(REPS))
    rep.to_csv(C.OUT / "bootstrap_replicates_seed20260930.csv.gz", index=False, float_format="%.10g",
               compression={"method": "gzip", "mtime": 0})

    # ---------- T2 / main metrics table
    rows = []
    dl_aucs, dl_cov = delong(y, [scores[k] for k in main_scores])
    for i, k in enumerate(main_scores):
        a_lo, a_hi = ci(bs[f"auroc|{k}"])
        p_lo, p_hi = ci(bs[f"auprc|{k}"])
        se = float(np.sqrt(dl_cov[i, i]))
        rows.append({"score": k, "n_sites": n, "n_pos": int(y.sum()), "n_proteins": K,
                     "auroc": pt[f"auroc|{k}"], "auroc_ci_low": a_lo, "auroc_ci_high": a_hi,
                     "auroc_delong_ci_low": dl_aucs[i] - 1.959963984540054 * se,
                     "auroc_delong_ci_high": dl_aucs[i] + 1.959963984540054 * se,
                     "auprc": pt[f"auprc|{k}"], "auprc_ci_low": p_lo, "auprc_ci_high": p_hi,
                     "auprc_chance": pt["prevalence"]})
    metrics = pd.DataFrame(rows)
    metrics.to_csv(C.OUT / "metrics_main.csv", index=False, float_format="%.10g")

    thr = {q: pt[f"thr|{q}"] for q in ("tp", "fn", "tn", "fp", "sn", "sp", "ba", "gmean", "mcc", "ppv")}
    thr_rows = []
    pub = {"sn": C.MANUSCRIPT["published_sn_pLMSNOSite"], "sp": C.MANUSCRIPT["published_sp_pLMSNOSite"],
           "ba": C.MANUSCRIPT["published_auroc_pLMSNOSite"], "gmean": C.MANUSCRIPT["published_auroc_pLMSNOSite"]}
    for q in ("sn", "sp", "ba", "gmean", "mcc"):
        lo, hi = ci(bs[f"thr|{q}"])
        thr_rows.append({"quantity": q, "value": thr[q], "ci_low": lo, "ci_high": hi,
                         "published_or_quoted": pub.get(q)})
    for q in ("tp", "fn", "tn", "fp", "ppv"):
        thr_rows.append({"quantity": q, "value": thr[q], "ci_low": None, "ci_high": None, "published_or_quoted": None})
    thr_df = pd.DataFrame(thr_rows)
    thr_df.to_csv(C.OUT / "plmsnosite_threshold_metrics.csv", index=False, float_format="%.10g")

    # ---------- T3 paired differences
    diff_rows = []
    pairs = [("pLMSNOSite", "VIS10"), ("pLMSNOSite", "DIG25"), ("rankavg_pLMSNOSite_DIG25", "pLMSNOSite"),
             ("rankavg_pLMSNOSite_VIS10", "pLMSNOSite"), ("DIG25", "VIS10"), ("Embedding_arm", "DIG25"),
             ("ProtT5_arm", "DIG25")]
    for a, b in pairs:
        for metric in ("auroc", "auprc"):
            d = bs[f"{metric}|{a}"] - bs[f"{metric}|{b}"]
            lo, hi = ci(d)
            row = {"comparison": f"{a} - {b}", "metric": metric, "difference": pt[f"{metric}|{a}"] - pt[f"{metric}|{b}"],
                   "ci_low": lo, "ci_high": hi, "share_replicates_le_0": C.boot_p_le(d),
                   "share_replicates_ge_0": float(np.mean(d >= 0))}
            if metric == "auroc":
                ia, ib = main_scores.index(a), main_scores.index(b)
                dd = dl_aucs[ia] - dl_aucs[ib]
                var = dl_cov[ia, ia] + dl_cov[ib, ib] - 2 * dl_cov[ia, ib]
                z = dd / np.sqrt(var)
                row.update({"delong_ci_low": dd - 1.959963984540054 * np.sqrt(var),
                            "delong_ci_high": dd + 1.959963984540054 * np.sqrt(var),
                            "delong_p_two_sided": float(2 * stats.norm.sf(abs(z)))})
            diff_rows.append(row)
    diffs = pd.DataFrame(diff_rows)
    diffs.to_csv(C.OUT / "paired_differences.csv", index=False, float_format="%.10g")

    # ---------- T4 detectability learning
    tert_rows = []
    for tcol, P, B, cuts in (("DIG25_tertile", pt, bs, (q1, q2)), ("VIS10_tertile", pt_v, bs_v, (v1, v2))):
        for t in (1, 2, 3):
            sel = df[tcol] == t
            for k in (("pLMSNOSite", "VIS10", "DIG25") if tcol == "DIG25_tertile" else ("pLMSNOSite", "VIS10")):
                key = f"tert{t}|auroc|{k}"
                lo, hi = ci(B[key])
                tert_rows.append({"stratifier": tcol, "tertile": t, "cut_points": f"{cuts[0]:.6g};{cuts[1]:.6g}",
                                  "n_sites": int(sel.sum()), "n_pos": int(df.loc[sel, "Target"].sum()),
                                  "n_proteins": int(df.loc[sel, "UniProt"].nunique()),
                                  "prevalence": float(df.loc[sel, "Target"].mean()),
                                  "score": k, "auroc": P[key], "ci_low": lo, "ci_high": hi})
    tert_df = pd.DataFrame(tert_rows)
    tert_df.to_csv(C.OUT / "detectability_learning_tertiles.csv", index=False, float_format="%.10g")

    sp_rows = []
    for a, b in sp_pairs:
        for cls in (0, 1):
            key = f"spearman|{a}~{b}|class{cls}"
            lo, hi = ci(bs[key])
            sel = y == cls
            sp_rows.append({"pair": f"{a} ~ {b}", "class": "positives" if cls else "negatives",
                            "n": int(sel.sum()), "spearman_rho": pt[key], "ci_low": lo, "ci_high": hi,
                            "scipy_p_unclustered": float(stats.spearmanr(scores[a][sel], scores[b][sel])[1])})
        rho_all = stats.spearmanr(scores[a], scores[b])
        sp_rows.append({"pair": f"{a} ~ {b}", "class": "all", "n": n, "spearman_rho": float(rho_all[0]),
                        "ci_low": None, "ci_high": None, "scipy_p_unclustered": float(rho_all[1])})
    sp_df = pd.DataFrame(sp_rows)
    sp_df.to_csv(C.OUT / "detectability_learning_spearman.csv", index=False, float_format="%.10g")

    # logistic regression (descriptive; protein-cluster-robust covariance)
    Lp, Ld, Lv = C.logit(scores["pLMSNOSite"]), C.logit(scores["DIG25"]), C.logit(scores["VIS10"])
    sd = {"logit_pLMSNOSite": Lp.std(ddof=1), "logit_DIG25": Ld.std(ddof=1), "logit_VIS10": Lv.std(ddof=1)}
    X_all = pd.DataFrame({"logit_pLMSNOSite": Lp, "logit_DIG25": Ld, "logit_VIS10": Lv})
    lr_rows = []
    specs = {"M1: pLMSNOSite": ["logit_pLMSNOSite"], "M2: DIG25": ["logit_DIG25"],
             "M3: pLMSNOSite + DIG25": ["logit_pLMSNOSite", "logit_DIG25"],
             "M4: VIS10": ["logit_VIS10"], "M5: pLMSNOSite + VIS10": ["logit_pLMSNOSite", "logit_VIS10"]}
    for name, cols in specs.items():
        X = sm.add_constant(X_all[cols])
        fit = sm.Logit(y, X).fit(disp=0, cov_type="cluster", cov_kwds={"groups": codes_})
        lin = np.asarray(X.to_numpy() @ fit.params.to_numpy())
        auc_in = roc_auc_score(y, lin)
        conf = fit.conf_int()
        for c in cols:
            b_, lo_, hi_ = fit.params[c], conf.loc[c, 0], conf.loc[c, 1]
            lr_rows.append({"model": name, "term": c, "beta_per_logit_unit": b_, "beta_ci_low": lo_, "beta_ci_high": hi_,
                            "sd_of_term": sd[c], "or_per_sd": np.exp(b_ * sd[c]), "or_per_sd_ci_low": np.exp(lo_ * sd[c]),
                            "or_per_sd_ci_high": np.exp(hi_ * sd[c]), "p_cluster_robust": fit.pvalues[c],
                            "mcfadden_r2": fit.prsquared, "in_sample_auroc_linear_predictor": auc_in,
                            "n": n, "n_clusters": K})
    lr_df = pd.DataFrame(lr_rows)
    b1 = lr_df.query("model == 'M1: pLMSNOSite'").beta_per_logit_unit.iloc[0]
    b3 = lr_df.query("model == 'M3: pLMSNOSite + DIG25' and term == 'logit_pLMSNOSite'").beta_per_logit_unit.iloc[0]
    b2 = lr_df.query("model == 'M2: DIG25'").beta_per_logit_unit.iloc[0]
    b32 = lr_df.query("model == 'M3: pLMSNOSite + DIG25' and term == 'logit_DIG25'").beta_per_logit_unit.iloc[0]
    lr_df.to_csv(C.OUT / "detectability_learning_logistic.csv", index=False, float_format="%.10g")
    attenuation = {"pLMSNOSite_beta_change_M1_to_M3": float((b1 - b3) / b1),
                   "DIG25_beta_change_M2_to_M3": float((b2 - b32) / b2)}

    # ---------- T5 recovery ratios
    rec_rows = []
    for det_name in ("VIS10", "DIG25"):
        num_b = bs[f"auroc|{det_name}"] - 0.5
        den_b = bs["auroc|pLMSNOSite"] - 0.5
        r_b = num_b / den_b
        lo, hi = ci(r_b)
        rec_rows.append({"detectability_model": det_name, "denominator": "pLMSNOSite re-scored score-based AUROC",
                         "denominator_value": pt["auroc|pLMSNOSite"],
                         "recovery": (pt[f"auroc|{det_name}"] - 0.5) / (pt["auroc|pLMSNOSite"] - 0.5),
                         "ci_low": lo, "ci_high": hi, "basis": "AUROC", "bootstrap": "paired, protein-clustered"})
        den_c = bs["auroc|pLMSNOSite_call_0.5"] - 0.5
        lo, hi = ci(num_b / den_c)
        rec_rows.append({"detectability_model": det_name,
                         "denominator": "pLMSNOSite AUROC of thresholded calls (= (Sn+Sp)/2; evaluate_model.py)",
                         "denominator_value": pt["auroc|pLMSNOSite_call_0.5"],
                         "recovery": (pt[f"auroc|{det_name}"] - 0.5) / (pt["auroc|pLMSNOSite_call_0.5"] - 0.5),
                         "ci_low": lo, "ci_high": hi, "basis": "AUROC", "bootstrap": "paired, protein-clustered"})
        for pub_name in ("pLMSNOSite", "PreSNO", "DeepNitro"):
            val = C.MANUSCRIPT[f"published_auroc_{pub_name}"]
            lo, hi = ci(num_b / (val - 0.5))
            rec_rows.append({"detectability_model": det_name, "denominator": f"{pub_name} tabulated (quoted, fixed)",
                             "denominator_value": val, "recovery": (pt[f"auroc|{det_name}"] - 0.5) / (val - 0.5),
                             "ci_low": lo, "ci_high": hi, "basis": "AUROC",
                             "bootstrap": "numerator only (denominator is a fixed quoted constant)",
                             "manuscript_value": C.MANUSCRIPT[f"recovery_{det_name}"][pub_name]})
        prev_b = bs["prevalence"]
        rp_b = (bs[f"auprc|{det_name}"] - prev_b) / (bs["auprc|pLMSNOSite"] - prev_b)
        lo, hi = ci(rp_b)
        rec_rows.append({"detectability_model": det_name, "denominator": "pLMSNOSite re-scored AUPRC (above prevalence)",
                         "denominator_value": pt["auprc|pLMSNOSite"],
                         "recovery": (pt[f"auprc|{det_name}"] - pt["prevalence"]) / (pt["auprc|pLMSNOSite"] - pt["prevalence"]),
                         "ci_low": lo, "ci_high": hi, "basis": "AUPRC above prevalence",
                         "bootstrap": "paired, protein-clustered"})
    rec_df = pd.DataFrame(rec_rows)
    rec_df.to_csv(C.OUT / "recovery_ratios.csv", index=False, float_format="%.10g")

    # ---------- R robustness
    rob_rows = []
    subsets = {
        "all sites": np.ones(n, bool),
        "excl. 4 proteins shared by ID with training": ~df.flag_shared_id.to_numpy(),
        "excl. proteins whose sequence occurs in training": ~df.flag_seq_in_train.to_numpy(),
        "excl. sites whose 37-residue window occurs in training": ~df.flag_window_in_train.to_numpy(),
        "excl. sites with cropped ProtT5 features (New_Position != Position)": ~df.flag_prott5_cropped.to_numpy(),
        "excl. duplicated test entry": ~df.flag_duplicate_site.to_numpy(),
        "excl. all of the above": ~(df.flag_shared_id | df.flag_seq_in_train | df.flag_window_in_train
                                    | df.flag_prott5_cropped | df.flag_duplicate_site).to_numpy(),
    }
    three = ("pLMSNOSite", "VIS10", "DIG25")
    for sname, mask in list(subsets.items()) + [("all sites, clusters = identical sequences", None)]:
        if mask is None:
            sub = df
            ccol = "seq_cluster"
        else:
            sub = df.loc[mask].reset_index(drop=True)
            ccol = "UniProt"
        sc = {k: sub[k].to_numpy(dtype=float) for k in three}
        P = point_stats(sub, sc, cluster_col=ccol)
        B, Kc = run_bootstrap(sub, sc, SEED, cluster_col=ccol)
        row = {"subset": sname, "n_sites": len(sub), "n_pos": int(sub.Target.sum()), "n_clusters": Kc}
        for k in three:
            lo, hi = ci(B[f"auroc|{k}"])
            row.update({f"auroc_{k}": P[f"auroc|{k}"], f"auroc_{k}_ci_low": lo, f"auroc_{k}_ci_high": hi})
            lo, hi = ci(B[f"auprc|{k}"])
            row.update({f"auprc_{k}": P[f"auprc|{k}"], f"auprc_{k}_ci_low": lo, f"auprc_{k}_ci_high": hi})
        for k in ("VIS10", "DIG25"):
            d = B["auroc|pLMSNOSite"] - B[f"auroc|{k}"]
            lo, hi = ci(d)
            row.update({f"auroc_diff_pLMSNOSite_minus_{k}": P["auroc|pLMSNOSite"] - P[f"auroc|{k}"],
                        f"auroc_diff_{k}_ci_low": lo, f"auroc_diff_{k}_ci_high": hi})
            d = B["auprc|pLMSNOSite"] - B[f"auprc|{k}"]
            lo, hi = ci(d)
            row.update({f"auprc_diff_pLMSNOSite_minus_{k}": P["auprc|pLMSNOSite"] - P[f"auprc|{k}"],
                        f"auprc_diff_{k}_ci_low": lo, f"auprc_diff_{k}_ci_high": hi})
            r = (B[f"auroc|{k}"] - 0.5) / (B["auroc|pLMSNOSite"] - 0.5)
            lo, hi = ci(r)
            row.update({f"recovery_{k}": (P[f"auroc|{k}"] - 0.5) / (P["auroc|pLMSNOSite"] - 0.5),
                        f"recovery_{k}_ci_low": lo, f"recovery_{k}_ci_high": hi})
        rob_rows.append(row)
    # within-protein AUROC and second seed
    wp_rows = []
    for k in wp_list:
        lo, hi = ci(bs[f"withinprot|{k}"])
        wp_rows.append({"quantity": f"pooled within-protein AUROC, {k}", "value": pt[f"withinprot|{k}"],
                        "ci_low": lo, "ci_high": hi})
    for k in ("VIS10", "DIG25"):
        d = bs["withinprot|pLMSNOSite"] - bs[f"withinprot|{k}"]
        lo, hi = ci(d)
        wp_rows.append({"quantity": f"pooled within-protein AUROC difference, pLMSNOSite - {k}",
                        "value": pt["withinprot|pLMSNOSite"] - pt[f"withinprot|{k}"], "ci_low": lo, "ci_high": hi})
        r = (bs[f"withinprot|{k}"] - 0.5) / (bs["withinprot|pLMSNOSite"] - 0.5)
        lo, hi = ci(r)
        wp_rows.append({"quantity": f"pooled within-protein recovery, {k}",
                        "value": (pt[f"withinprot|{k}"] - 0.5) / (pt["withinprot|pLMSNOSite"] - 0.5),
                        "ci_low": lo, "ci_high": hi})
    for k in three:
        for metric in ("auroc", "auprc"):
            lo, hi = ci(bs_alt[f"{metric}|{k}"])
            wp_rows.append({"quantity": f"{metric} {k}, bootstrap seed {SEED_ALT}", "value": pt[f"{metric}|{k}"],
                            "ci_low": lo, "ci_high": hi})
    for k in ("VIS10", "DIG25"):
        r = (bs_alt[f"auroc|{k}"] - 0.5) / (bs_alt["auroc|pLMSNOSite"] - 0.5)
        lo, hi = ci(r)
        wp_rows.append({"quantity": f"recovery {k}, bootstrap seed {SEED_ALT}",
                        "value": (pt[f"auroc|{k}"] - 0.5) / (pt["auroc|pLMSNOSite"] - 0.5), "ci_low": lo, "ci_high": hi})
        d = bs_alt["auroc|pLMSNOSite"] - bs_alt[f"auroc|{k}"]
        lo, hi = ci(d)
        wp_rows.append({"quantity": f"auroc difference pLMSNOSite - {k}, bootstrap seed {SEED_ALT}",
                        "value": pt["auroc|pLMSNOSite"] - pt[f"auroc|{k}"], "ci_low": lo, "ci_high": hi})
    rob = pd.DataFrame(rob_rows)
    rob.to_csv(C.OUT / "robustness_subsets.csv", index=False, float_format="%.10g")
    rob2 = pd.DataFrame(wp_rows)
    rob2.to_csv(C.OUT / "robustness_other.csv", index=False, float_format="%.10g")
    flags = df[["row", "UniProt", "Position", "Target", "DIG25_tertile", "VIS10_tertile", "seq_cluster",
                "flag_shared_id", "flag_seq_in_train", "flag_window_in_train", "flag_prott5_cropped",
                "flag_duplicate_site", "rankavg_pLMSNOSite_DIG25", "rankavg_pLMSNOSite_VIS10"]]
    flags.to_csv(C.OUT / "site_flags_and_combined_scores.csv", index=False, float_format="%.17g")
    # window-label conflict table
    tr_w = pd.DataFrame({"win": [win(s, p) for s, p in zip(train.sequences, train.Position)],
                         "train_label": train.Target})
    te_w = pd.DataFrame({"win": test_w, "test_label": df.Target})
    mm = te_w.merge(tr_w.groupby("win").train_label.mean().rename("train_label_mean"), left_on="win", right_index=True)
    conflict = {"test_windows_in_train": int(len(mm)),
                "label_conflicts": int((mm.test_label != mm.train_label_mean.round()).sum()),
                "crosstab": pd.crosstab(mm.test_label, mm.train_label_mean).to_dict()}

    summary = {
        "label": "POST HOC revision analysis 2026-09-30 (G_plmsnosite_rescore); not registered, not pre-specified",
        "bootstrap": {"replicates": REPS, "seed": SEED, "unit": "protein (UniProt accession), sorted labels",
                      "draws": "numpy default_rng(seed); rng.integers(0, K, K) per replicate", "interval": "95% percentile",
                      "estimand": "conditional on the fixed trained models; models are not refit within replicates"},
        "tertile_cut_points_DIG25": [q1, q2], "tertile_cut_points_VIS10": [v1, v2],
        "point": {k: (float(v) if isinstance(v, (float, np.floating, int, np.integer)) else v) for k, v in pt.items()},
        "logistic_attenuation": attenuation,
        "window_label_conflict": conflict,
        "flags_counts": {c: int(df[c].sum()) for c in df.columns if c.startswith("flag_")},
        "runtime_s": round(time.time() - t0, 1),
    }
    (C.OUT / "comparison_summary.json").write_text(json.dumps(summary, indent=1, default=str), encoding="utf-8")

    # input hashes of this step (item-level provenance.json is written by 40_provenance.py)
    (C.OUT / "comparison_inputs_sha256.json").write_text(json.dumps(manifest, indent=1), encoding="utf-8")

    # console summary
    pd.set_option("display.width", 250)
    print(metrics.round(4).to_string())
    print(thr_df.round(4).to_string())
    print(diffs.round(4).to_string())
    print(tert_df.round(4).to_string())
    print(sp_df.round(4).to_string())
    print(lr_df[["model", "term", "beta_per_logit_unit", "beta_ci_low", "beta_ci_high", "or_per_sd", "or_per_sd_ci_low",
                 "or_per_sd_ci_high", "p_cluster_robust", "in_sample_auroc_linear_predictor"]].round(4).to_string())
    print(attenuation)
    print(rec_df.round(4).to_string())
    print(rob.round(4).T.to_string())
    print(rob2.round(4).to_string())
    print(json.dumps(conflict, indent=1, default=str))
    print("runtime", round(time.time() - t0, 1), "s")


if __name__ == "__main__":
    main()
