"""Helpers for steps 8/8b (POST HOC, revision after adversarial verification, 2026-09-30).

Array-based versions of the rank-biserial, stratified and matched statistics of s04/s04b, so that the
calibrated null simulations of s08b can run fast. The matching reproduces s04b.match_ab_opp exactly
(positives in ascending log10 iBAQ, stable order; nearest unused negative in the same stratum).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats


def rank_biserial(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    if len(x) < 2 or len(y) < 2:
        return float("nan")
    rk = stats.rankdata(np.concatenate([x, y]))
    U = rk[:len(x)].sum() - len(x) * (len(x) + 1) / 2.0
    return float(2.0 * U / (len(x) * len(y)) - 1.0)


def strata_index(codes):
    """dict code -> indices (in original order)."""
    codes = np.asarray(codes)
    order = np.argsort(codes, kind="mergesort")
    cs = codes[order]
    cuts = np.flatnonzero(np.diff(cs)) + 1
    return {int(g[0]): o for g, o in zip(np.split(cs, cuts), np.split(order, cuts))}


def strat_r(lab, vals, sidx, min_n=10):
    rs = []
    for _, idx in sidx.items():
        l = lab[idx]
        a = vals[idx[l == 1]]; b = vals[idx[l == 0]]
        if len(a) >= min_n and len(b) >= min_n:
            rs.append(rank_biserial(a, b))
    return (float(np.mean(rs)) if rs else float("nan")), len(rs)


def sorted_strata_index(codes, x):
    """dict code -> indices sorted by x ascending (stable, ties in original order)."""
    order = np.argsort(x, kind="mergesort")
    codes = np.asarray(codes)[order]
    out = {}
    for c in np.unique(codes):
        out[int(c)] = order[codes == c]
    return out


def match_within(lab, x, ssidx):
    """1:1 nearest-x matching without replacement inside strata. Returns (pos_idx, neg_idx)."""
    P, N = [], []
    for _, idx in ssidx.items():
        l = lab[idx]
        pos = idx[l == 1]; neg = idx[l == 0]
        if len(pos) == 0 or len(neg) == 0:
            continue
        nv = x[neg]; used = np.zeros(len(neg), bool)
        for i in pos:
            d = np.abs(nv - x[i]); d[used] = np.inf
            j = int(d.argmin())
            if not np.isfinite(d[j]):
                break
            used[j] = True; P.append(i); N.append(neg[j])
    return np.asarray(P, int), np.asarray(N, int)


def decile_codes(v):
    e = np.percentile(v, np.linspace(0, 100, 11)); e[0] -= 1e-9; e[-1] += 1e-9
    return np.clip(np.searchsorted(e, v, side="left") - 1, 0, 9)


def tertile_codes(n):
    q1, q2 = np.quantile(n, [1 / 3, 2 / 3])
    return np.where(n <= q1, 0, np.where(n <= q2, 1, 2)), (float(q1), float(q2))


class Design:
    """Fixed (label-independent) arrays and strata for one data set."""

    def __init__(self, df, key_ndet="n_det", cap=20, tert_cuts=(3.0, 6.0)):
        """tert_cuts: opportunity tertile cut points; (3, 6) are the s04 tertiles of all 7,692 groups
        (n_det <= 3 / <= 6 / > 6); None derives them from this data set."""
        self.x = df["log10_ibaq"].to_numpy(float)
        self.length = df["length"].to_numpy(float)
        self.ncys = df["n_cys"].to_numpy(float)
        self.ndet = df[key_ndet].to_numpy(float)
        self.nundet = self.ncys - self.ndet
        self.dec = decile_codes(self.x)
        if tert_cuts is None:
            self.tert, self.tert_cuts = tertile_codes(self.ndet)
        else:
            q1, q2 = tert_cuts
            self.tert = np.where(self.ndet <= q1, 0, np.where(self.ndet <= q2, 1, 2))
            self.tert_cuts = (float(q1), float(q2))
        self.opp = np.minimum(self.ndet, cap).astype(int)
        self.s_dec = strata_index(self.dec)
        self.s_dt = strata_index(self.dec * 10 + self.tert)
        self.s_dx = strata_index(self.dec * 1000 + self.ndet.astype(int))
        self.ss_opp = sorted_strata_index(self.opp, self.x)
        self.X_m3u = np.column_stack([np.ones(len(self.x)), self.x, np.log2(1 + self.ndet),
                                      np.log2(self.length), self.nundet])

    def stats(self, lab, with_logit=True):
        lab = np.asarray(lab, int)
        out = {"n_pos": int(lab.sum())}
        out["r_dec_len"] = strat_r(lab, self.length, self.s_dec)[0]
        out["r_dec_cys"] = strat_r(lab, self.ncys, self.s_dec)[0]
        out["r_dt_len"] = strat_r(lab, self.length, self.s_dt)[0]
        out["r_dt_cys"] = strat_r(lab, self.ncys, self.s_dt)[0]
        out["r_cell_undet"], out["cells_dx"] = strat_r(lab, self.nundet, self.s_dx)
        P, N = match_within(lab, self.x, self.ss_opp)
        out["pairs"] = int(len(P))
        out["r_m_len"] = rank_biserial(self.length[P], self.length[N])
        out["r_m_cys"] = rank_biserial(self.ncys[P], self.ncys[N])
        out["r_m_undet"] = rank_biserial(self.nundet[P], self.nundet[N])
        out["mean_undet_pos"] = float(self.nundet[P].mean()) if len(P) else float("nan")
        out["mean_undet_neg"] = float(self.nundet[N].mean()) if len(N) else float("nan")
        if with_logit:
            try:
                m = sm.Logit(lab, self.X_m3u).fit(disp=0, maxiter=100)
                out["b_len_M3u"] = float(m.params[3]); out["b_undet_M3u"] = float(m.params[4])
            except Exception:  # pragma: no cover - separation would be reported as NaN
                out["b_len_M3u"] = out["b_undet_M3u"] = float("nan")
        return out


def simulate_chunk(args):
    """Worker for s08b: args = (design_df, p_null, reps, seed_entropy, spawn_key)."""
    df, p, reps, entropy, key = args
    D = Design(df)
    rng = np.random.default_rng(np.random.SeedSequence(entropy, spawn_key=key))
    out = []
    for _ in range(reps):
        y = (rng.random(len(p)) < p).astype(int)
        out.append(D.stats(y))
    return out


def features_for_window(seqs, seq_accessions, site_positions, digest, enz, lo, hi, mc):
    """n_det and S_map for every group under one detectability window."""
    nd, smap = [], []
    for acc, st in zip(seq_accessions, site_positions):
        s = seqs[acc]
        det = np.zeros(len(s), bool)
        for a0, b0 in digest(s, enz, missed=mc, min_len=lo, max_len=hi):
            det[a0 - 1:b0] = True
        cpos = [i for i, ch in enumerate(s) if ch == "C"]
        nd.append(int(det[cpos].sum()) if cpos else 0)
        smap.append(int(any(det[q - 1] for q in st)))
    return np.asarray(nd), np.asarray(smap)


def or_rows(m, label, model, terms, extra=None):
    ci = m.conf_int()
    rows = []
    for t in terms:
        if t in m.params.index:
            rows.append({"label": label, "model": model, "term": t, "OR": float(np.exp(m.params[t])),
                         "ci_low": float(np.exp(ci.loc[t, 0])), "ci_high": float(np.exp(ci.loc[t, 1])),
                         "p": float(m.pvalues[t]), "aic": float(m.aic), "llf": float(m.llf),
                         "n": int(m.nobs), **(extra or {})})
    return rows
