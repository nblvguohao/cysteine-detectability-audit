"""Propensity model for the claim retest (v0.2.0).

Ridge logistic regression by Newton-Raphson, minimising 0.5*||w||^2 + C * sum(log-loss) over standardised
covariates (intercept unpenalised). Missing covariate values are replaced by the column median (declared, simple,
and recorded in the audit); constant columns are dropped. All products use numpy.einsum rather than the @ operator,
so no BLAS call is made (macOS Accelerate raises spurious matmul warnings; CLAUDE.md s.9.19.1).
The model is fitted once on the full data; bootstrap replicates reuse the fitted strata (conventional; it slightly
understates uncertainty, stated in the report).
"""
from __future__ import annotations

import numpy as np

from .constants import PROPENSITY_C, PROPENSITY_STRATA


def fit_logit(X, y, C=PROPENSITY_C, max_iter=100, tol=1e-9):
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    n, d = X.shape
    Z = np.column_stack([np.ones(n), X])
    w = np.zeros(d + 1)
    pen = np.full(d + 1, 1.0 / C)
    pen[0] = 0.0
    for _ in range(max_iter):
        eta = np.einsum("ij,j->i", Z, w)
        p = 1.0 / (1.0 + np.exp(-np.clip(eta, -35, 35)))
        g = np.einsum("ij,i->j", Z, p - y) + pen * w
        W = p * (1 - p)
        H = np.einsum("ij,i,ik->jk", Z, W, Z) + np.diag(pen)
        step = np.linalg.solve(H, g)
        w = w - step
        if np.max(np.abs(step)) < tol:
            break
    eta = np.einsum("ij,j->i", Z, w)
    return w, eta


def prepare(columns):
    """columns: dict name -> array (may contain NaN). Returns standardised matrix, kept names, notes."""
    kept, cols, notes = [], [], []
    for name, v in columns.items():
        v = np.asarray(v, dtype=float).copy()
        miss = ~np.isfinite(v)
        if miss.all():
            notes.append(f"covariate {name} is all missing; dropped")
            continue
        if miss.any():
            v[miss] = np.median(v[~miss])
            notes.append(f"covariate {name}: {int(miss.sum())} missing values set to the median")
        sd = v.std()
        if sd == 0:
            notes.append(f"covariate {name} is constant; dropped")
            continue
        cols.append((v - v.mean()) / sd)
        kept.append(name)
    X = np.column_stack(cols) if cols else np.zeros((len(next(iter(columns.values()))), 0))
    return X, kept, notes


def strata_from_logit(eta, k=PROPENSITY_STRATA):
    cuts = np.quantile(eta, np.linspace(0, 1, k + 1)[1:-1])
    return np.clip(np.searchsorted(cuts, eta, side="right"), 0, k - 1)
