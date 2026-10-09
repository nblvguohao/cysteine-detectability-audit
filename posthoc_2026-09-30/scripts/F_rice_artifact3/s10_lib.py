"""Helpers for steps 10/10b (POST HOC; revision round 2, 2026-09-30).

Why this exists. Verifier round 2 (major 1) showed that the evidence used in round 1 cannot tell the
number of THEORETICALLY DETECTABLE cysteines (n_det, 'opportunity' under the manuscript's detectability
model) from the TOTAL number of cysteines (n_cys; Spearman 0.94). Step 10 therefore compares the two
counts directly, and step 10b adds total-count nulls beside the detectable-count nulls of s08b.

Design2 extends s08_lib.Design (whose statistics are kept unchanged, so that label sets drawn with the
s08b seeds reproduce the s08b statistics exactly) with:
  mirror design: 1:1 nearest-log10-iBAQ matching WITHOUT replacement inside strata of exact TOTAL
      cysteine count (capped at 25), then rank-biserial r of n_det, n_undet and length
      (if detectability matters beyond total count, positives should still carry more detectable Cys);
  r_cellc_det: mean r of n_det over iBAQ-decile x exact-n_cys cells (>= 10 per arm);
  logistic comparisons on the same groups:
      D  y ~ log10 iBAQ + log2 length + log2(1 + n_det)
      T  y ~ log10 iBAQ + log2 length + log2(1 + n_cys)
      B  y ~ log10 iBAQ + log2 length + log2(1 + n_det) + log2(1 + n_cys)
      L  y ~ log10 iBAQ + log2 length + n_det + n_undet          (symmetric linear counts)
    dAIC_D_minus_T = AIC(D) - AIC(T) (> 0: total count fits better; D and T have equal df);
    LR of each count given the other (B vs T, B vs D); b_det_B, b_cys_B; b_ndet_L, b_nundet_L and their
    difference; b_len_T (length at fixed total count).
Everything here is post hoc and was written after the verifier's round-2 report.
"""
from __future__ import annotations

import numpy as np
import statsmodels.api as sm

from s08_lib import Design, match_within, rank_biserial, sorted_strata_index, strat_r, strata_index

CAP_NCYS = 25


def _logit_llf(y, X):
    """Return (params, llf) of a logistic fit, or (None, nan) on failure (separation etc.)."""
    try:
        m = sm.Logit(y, X).fit(disp=0, maxiter=200)
        if not np.isfinite(m.llf):
            return None, float("nan")
        return m.params, float(m.llf)
    except Exception:  # pragma: no cover - reported as NaN
        return None, float("nan")


class Design2(Design):
    """Design plus the detectable-versus-total statistics of round 2."""

    def __init__(self, df, key_ndet="n_det", cap=20, tert_cuts=(3.0, 6.0), cap_ncys=CAP_NCYS):
        super().__init__(df, key_ndet=key_ndet, cap=cap, tert_cuts=tert_cuts)
        self.ncys_c = np.minimum(self.ncys, cap_ncys).astype(int)
        self.ss_cys = sorted_strata_index(self.ncys_c, self.x)
        self.s_dcx = strata_index(self.dec * 1000 + self.ncys_c)
        one = np.ones(len(self.x))
        l2len = np.log2(self.length)
        l2d = np.log2(1 + self.ndet)
        l2c = np.log2(1 + self.ncys)
        self.X_D = np.column_stack([one, self.x, l2len, l2d])
        self.X_T = np.column_stack([one, self.x, l2len, l2c])
        self.X_B = np.column_stack([one, self.x, l2len, l2d, l2c])
        self.X_L = np.column_stack([one, self.x, l2len, self.ndet, self.nundet])

    def stats2(self, lab, with_logit=True):
        """Round-2 statistics only (the Design.stats statistics are computed by stats())."""
        lab = np.asarray(lab, int)
        out = {}
        P, N = match_within(lab, self.x, self.ss_cys)
        out["pairs_c"] = int(len(P))
        out["r_mc_det"] = rank_biserial(self.ndet[P], self.ndet[N])
        out["r_mc_undet"] = rank_biserial(self.nundet[P], self.nundet[N])
        out["r_mc_len"] = rank_biserial(self.length[P], self.length[N])
        out["mean_det_pos_c"] = float(self.ndet[P].mean()) if len(P) else float("nan")
        out["mean_det_neg_c"] = float(self.ndet[N].mean()) if len(N) else float("nan")
        out["frac_pairs_equal_ncys"] = float((self.ncys[P] == self.ncys[N]).mean()) if len(P) else float("nan")
        out["r_cellc_det"], out["cells_dc"] = strat_r(lab, self.ndet, self.s_dcx)
        if with_logit:
            pD, lD = _logit_llf(lab, self.X_D)
            pT, lT = _logit_llf(lab, self.X_T)
            pB, lB = _logit_llf(lab, self.X_B)
            pL, lL = _logit_llf(lab, self.X_L)
            out["dAIC_D_minus_T"] = float(2 * (lT - lD))
            out["LR_det_given_T"] = float(2 * (lB - lT))
            out["LR_cys_given_D"] = float(2 * (lB - lD))
            out["b_det_B"] = float(pB[3]) if pB is not None else float("nan")
            out["b_cys_B"] = float(pB[4]) if pB is not None else float("nan")
            out["b_len_T"] = float(pT[2]) if pT is not None else float("nan")
            out["b_ndet_L"] = float(pL[3]) if pL is not None else float("nan")
            out["b_nundet_L"] = float(pL[4]) if pL is not None else float("nan")
            out["d_det_minus_undet_L"] = out["b_ndet_L"] - out["b_nundet_L"]
        return out

    def stats_all(self, lab, with_logit=True):
        out = self.stats(lab, with_logit=with_logit)
        out.update(self.stats2(lab, with_logit=with_logit))
        return out


def simulate_chunk2(args):
    """Worker for s10b: args = (design_df, p_null, reps, seed_entropy, spawn_key).

    The random stream is consumed exactly as in s08_lib.simulate_chunk (one uniform vector per label
    set), so a combination drawn with an s08b entropy reproduces the s08b label sets.
    """
    df, p, reps, entropy, key = args
    D = Design2(df)
    rng = np.random.default_rng(np.random.SeedSequence(entropy, spawn_key=key))
    out = []
    for _ in range(reps):
        y = (rng.random(len(p)) < p).astype(int)
        out.append(D.stats_all(y))
    return out
