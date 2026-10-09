"""VERIFIER round 4 (POST HOC) helpers for item F_rice_artifact3: own binary-regression MLE (Newton for logit,
Fisher scoring for cloglog, step-halving), own rank-biserial / stratified / matched statistics.
Written independently of the item's s08_lib / s10_lib (only the definitions in the report were used)."""
from __future__ import annotations

import numpy as np
import patsy
from scipy import stats
from scipy.special import expit


def _ll(y, eta, link):
    if link == "logit":
        return float(np.sum(y * eta - np.logaddexp(0.0, eta)))
    mu = -np.expm1(-np.exp(np.minimum(eta, 30)))
    mu = np.clip(mu, 1e-300, 1.0)
    return float(np.sum(y * np.log(mu) - (1 - y) * np.exp(np.minimum(eta, 30))))


def fit_binary(y, X, link="logit", maxit=200, tol=1e-10):
    """MLE of P(y=1) = g^-1(X b). Returns dict(b, se, llf, k, converged, cov). X must have full column rank."""
    y = np.asarray(y, float); X = np.asarray(X, float)
    n, k = X.shape
    assert np.linalg.matrix_rank(X) == k, "design not full rank"
    b = np.zeros(k)
    if link == "cloglog":
        pbar = min(max(y.mean(), 1e-4), 1 - 1e-4)
        # start: intercept-like via least squares of the working response
        b = np.linalg.lstsq(X, np.full(n, np.log(-np.log(1 - pbar))), rcond=None)[0]
    ll = _ll(y, X @ b, link)
    conv = False
    for _ in range(maxit):
        eta = X @ b
        if link == "logit":
            mu = expit(eta)
            w = mu * (1 - mu)
            grad = X.T @ (y - mu)
        else:
            e = np.exp(np.minimum(eta, 30))
            mu = np.clip(-np.expm1(-e), 1e-12, 1 - 1e-12)
            dmu = e * np.exp(-e)
            w = dmu ** 2 / (mu * (1 - mu))
            grad = X.T @ ((y - mu) * dmu / (mu * (1 - mu)))
        Hm = X.T @ (X * w[:, None])
        try:
            step = np.linalg.solve(Hm, grad)
        except np.linalg.LinAlgError:          # quasi-complete separation: minimum-norm step
            step = np.linalg.lstsq(Hm, grad, rcond=None)[0]
        t = 1.0
        while True:
            b_new = b + t * step
            ll_new = _ll(y, X @ b_new, link)
            if ll_new >= ll - 1e-12 or t < 1e-8:
                break
            t /= 2
        if abs(ll_new - ll) < tol and np.max(np.abs(b_new - b)) < 1e-7:
            b, ll, conv = b_new, ll_new, True
            break
        b, ll = b_new, ll_new
    eta = X @ b
    if link == "logit":
        mu = expit(eta); w = mu * (1 - mu)
    else:
        e = np.exp(np.minimum(eta, 30)); mu = np.clip(-np.expm1(-e), 1e-12, 1 - 1e-12)
        w = (e * np.exp(-e)) ** 2 / (mu * (1 - mu))
    cov = np.linalg.pinv(X.T @ (X * w[:, None]))
    return {"b": b, "se": np.sqrt(np.abs(np.diag(cov))), "cov": cov, "llf": ll, "k": k, "converged": conv,
            "aic": 2 * k - 2 * ll, "fitted": mu}


def spline_basis(x, df=5):
    """patsy natural cubic regression spline basis cr(x, df) WITHOUT intercept (its columns sum to one, so it
    spans the constant): the same column space as the item's 'Intercept + cr(x, df=5)'."""
    B = np.asarray(patsy.dmatrix(f"cr(x, df={df}) - 1", {"x": np.asarray(x, float)}))
    return B


def dummies(v, drop_first=True):
    lev = np.unique(v)
    cols = [(v == L).astype(float) for L in (lev[1:] if drop_first else lev)]
    return np.column_stack(cols) if cols else np.zeros((len(v), 0)), lev


def rank_biserial(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    U = stats.mannwhitneyu(a, b, alternative="two-sided", method="asymptotic").statistic
    return float(2 * U / (len(a) * len(b)) - 1)


def deciles(x):
    e = np.percentile(x, np.linspace(0, 100, 11)); e[0] -= 1e-9; e[-1] += 1e-9
    return np.clip(np.searchsorted(e, x, side="left") - 1, 0, 9)


def strat_mean_r(y, v, codes, min_n=10):
    rs = []
    for c in np.unique(codes):
        m = codes == c
        a = v[m & (y == 1)]; b = v[m & (y == 0)]
        if len(a) >= min_n and len(b) >= min_n:
            rs.append(rank_biserial(a, b))
    return float(np.mean(rs)) if rs else float("nan")


def greedy_match(y, x, strata):
    """1:1 nearest-x matching without replacement within strata; positives in ascending x (stable),
    nearest unused negative (first index on ties). Returns (pos_idx, neg_idx)."""
    P, N = [], []
    order = np.argsort(x, kind="mergesort")
    for c in np.unique(strata):
        idx = order[strata[order] == c]
        pos = idx[y[idx] == 1]; neg = idx[y[idx] == 0]
        if len(pos) == 0 or len(neg) == 0:
            continue
        nx = x[neg].copy()
        free = np.ones(len(neg), bool)
        for i in pos:
            if not free.any():
                break
            d = np.where(free, np.abs(nx - x[i]), np.inf)
            j = int(np.argmin(d))
            free[j] = False
            P.append(i); N.append(neg[j])
    return np.asarray(P, int), np.asarray(N, int)
