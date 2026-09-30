"""Statistics shared by the artefact tests.

All intervals are cluster-bootstrap percentile intervals: whole clusters (proteins by default, or
homology components when the input supplies a `cluster` column) are resampled with replacement.
Every statistic here is written in terms of per-cluster sufficient statistics so that 5000 replicates
stay cheap: a replicate is a vector of cluster multiplicities `w`, and a statistic is a function of
`w @ per_cluster_counts`.

Determinism: clusters are always sorted before indexing, so the result does not depend on set/dict
iteration order or PYTHONHASHSEED (the failure recorded in the analysis tree's CLAUDE.md s.9.64.2).
"""
from __future__ import annotations

import math

import numpy as np


def cluster_index(clusters):
    """Map cluster labels to 0..K-1 in sorted label order. Returns (codes, sorted_labels)."""
    labels = sorted(set(clusters))
    pos = {c: i for i, c in enumerate(labels)}
    return np.fromiter((pos[c] for c in clusters), dtype=np.int64, count=len(clusters)), labels


def multiplicities(n_clusters, reps, seed, chunk=250):
    """Yield (reps_in_chunk x n_clusters) arrays of bootstrap cluster multiplicities."""
    rng = np.random.default_rng(seed)
    done = 0
    while done < reps:
        m = min(chunk, reps - done)
        draws = rng.integers(0, n_clusters, size=(m, n_clusters))
        flat = (np.arange(m)[:, None] * n_clusters + draws).ravel()
        # int64 multiplicities: exact integer arithmetic, and integer matmul does not go through BLAS
        # (macOS Accelerate raises spurious matmul RuntimeWarnings on float64; see CLAUDE.md s.9.19.1)
        yield np.bincount(flat, minlength=m * n_clusters).reshape(m, n_clusters).astype(np.int64)
        done += m


def percentile_interval(values, level):
    a = (1.0 - level) / 2.0
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return float("nan"), float("nan")
    return float(np.quantile(v, a)), float(np.quantile(v, 1.0 - a))


def bootstrap_p(values, null):
    """Two-sided bootstrap p for H0: statistic == null (reported only; never used for a status)."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return float("nan")
    lo = float(np.mean(v <= null))
    hi = float(np.mean(v >= null))
    return float(max(min(1.0, 2.0 * min(lo, hi)), 1.0 / (v.size + 1)))


def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return c - h, c + h


def log2_or(a, b, c, d):
    """Haldane-corrected log2 odds ratio; a,b = flagged/unflagged positives, c,d = background."""
    a, b, c, d = (np.asarray(x, dtype=float) + 0.5 for x in (a, b, c, d))
    return np.log2((a / b) / (c / d))


def two_by_two_counts(codes, n_clusters, is_pos, flag, mask=None):
    """Per-cluster counts (K x 4): positives flagged, positives unflagged, background flagged, background
    unflagged. `mask` selects the rows that belong to the comparison at all."""
    is_pos = np.asarray(is_pos, dtype=bool)
    flag = np.asarray(flag, dtype=bool)
    keep = np.ones(len(codes), dtype=bool) if mask is None else np.asarray(mask, dtype=bool)
    out = np.zeros((n_clusters, 4), dtype=np.int64)
    for j, sel in enumerate((is_pos & flag, is_pos & ~flag, ~is_pos & flag, ~is_pos & ~flag)):
        s = sel & keep
        out[:, j] = np.bincount(codes[s], minlength=n_clusters)
    return out


def boot_log2_or(counts, reps, seed, level):
    """counts: K x 4 per-cluster table. Returns point, (lo, hi), bootstrap p (vs 0), replicate array."""
    tot = counts.sum(axis=0)
    point = float(log2_or(*tot))
    vals = []
    for w in multiplicities(counts.shape[0], reps, seed):
        s = w @ counts
        vals.append(log2_or(s[:, 0], s[:, 1], s[:, 2], s[:, 3]))
    vals = np.concatenate(vals)
    return point, percentile_interval(vals, level), bootstrap_p(vals, 0.0), vals


def boot_proportion(codes, n_clusters, success, reps, seed, level):
    """Share of rows that are successes, clustered. Returns point, (lo, hi), replicate array."""
    success = np.asarray(success, dtype=float)
    k = np.bincount(codes[success.astype(bool)], minlength=n_clusters).astype(np.int64)
    n = np.bincount(codes, minlength=n_clusters).astype(np.int64)
    point = float(k.sum() / n.sum()) if n.sum() else float("nan")
    vals = []
    for w in multiplicities(n_clusters, reps, seed):
        den = (w @ n).astype(float)
        with np.errstate(invalid="ignore", divide="ignore"):
            vals.append((w @ k).astype(float) / den)
    vals = np.concatenate(vals)
    return point, percentile_interval(vals, level), vals


def _weighted_auc(score_group, pos_w, neg_w, n_groups):
    """AUC from per-row weights; rows grouped by tied score (score_group ascending, 0..G-1)."""
    pg = np.bincount(score_group, weights=pos_w, minlength=n_groups)
    ng = np.bincount(score_group, weights=neg_w, minlength=n_groups)
    below = np.concatenate(([0.0], np.cumsum(ng)[:-1]))
    wp, wn = pg.sum(), ng.sum()
    if wp == 0 or wn == 0:
        return float("nan")
    return float((pg * (below + 0.5 * ng)).sum() / (wp * wn))


def auc(scores, labels):
    """Tie-corrected rank AUC: P(score of a positive > score of a background row) + 0.5 P(tie)."""
    scores = np.asarray(scores, dtype=float)
    labels = np.asarray(labels, dtype=bool)
    _, grp = np.unique(scores, return_inverse=True)
    return _weighted_auc(grp, labels.astype(float), (~labels).astype(float), int(grp.max()) + 1 if grp.size else 0)


def boot_auc(scores, labels, codes, n_clusters, reps, seed, level):
    scores = np.asarray(scores, dtype=float)
    labels = np.asarray(labels, dtype=bool)
    _, grp = np.unique(scores, return_inverse=True)
    g = int(grp.max()) + 1
    point = _weighted_auc(grp, labels.astype(float), (~labels).astype(float), g)
    pos_i, neg_i = labels.astype(float), (~labels).astype(float)
    vals = []
    for w in multiplicities(n_clusters, reps, seed):
        for r in range(w.shape[0]):
            rw = w[r][codes].astype(float)
            vals.append(_weighted_auc(grp, rw * pos_i, rw * neg_i, g))
    vals = np.asarray(vals)
    return point, percentile_interval(vals, level), bootstrap_p(vals, 0.5), vals


def mh_log2_or(strata_counts):
    """Mantel-Haenszel pooled log2 OR over strata. strata_counts: S x 4 (a, b, c, d)."""
    s = np.asarray(strata_counts, dtype=float)
    n = s.sum(axis=1)
    ok = n > 0
    s, n = s[ok], n[ok]
    num = (s[:, 0] * s[:, 3] / n).sum()
    den = (s[:, 1] * s[:, 2] / n).sum()
    if num <= 0 or den <= 0:
        # add 0.5 per cell in the informative strata, the same convention as log2_or
        s = s + 0.5
        n = s.sum(axis=1)
        num = (s[:, 0] * s[:, 3] / n).sum()
        den = (s[:, 1] * s[:, 2] / n).sum()
    return float(math.log2(num / den))


def bh_qvalues(pvals):
    """Benjamini-Hochberg q-values for a list that may contain None (returned as None)."""
    idx = [i for i, p in enumerate(pvals) if p is not None and np.isfinite(p)]
    q = [None] * len(pvals)
    if not idx:
        return q
    ps = np.array([pvals[i] for i in idx])
    order = np.argsort(ps)
    m = len(ps)
    ranked = ps[order] * m / (np.arange(m) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    ranked = np.minimum(ranked, 1.0)
    for k, o in enumerate(order):
        q[idx[o]] = float(ranked[k])
    return q
