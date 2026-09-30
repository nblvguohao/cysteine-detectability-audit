"""Greedy 1:1 nearest-neighbour matching without replacement on a scalar score (POST HOC, item
D_empirical_background). Positives are processed in random order; ties within the pool are broken at random
(the pool is shuffled before a stable sort); an exact tie in distance between the nearest free control on the
left and on the right is broken at random."""
from __future__ import annotations

import numpy as np

EPS = 1e-6


def logit(p):
    p = np.clip(np.asarray(p, dtype=float), EPS, 1 - EPS)
    return np.log(p / (1 - p))


class NextFree:
    """Union-find over sorted positions: nearest unused index to the right (or left) of i."""

    def __init__(self, n):
        self.r = np.arange(n + 1)   # sentinel n = none to the right
        self.l = np.arange(n + 1)   # index shifted by +1; sentinel 0 = none to the left

    @staticmethod
    def _find(a, i):
        root = i
        while a[root] != root:
            root = a[root]
        while a[i] != root:
            a[i], i = root, a[i]
        return root

    def right(self, i):
        return self._find(self.r, i)

    def left(self, i):
        return self._find(self.l, i + 1) - 1

    def take(self, i):
        self.r[i] = i + 1
        self.l[i + 1] = i


def nn_match(pos_score, pool_score, rng, caliper=None):
    """Returns, per positive, the index into the pool of its control (or -1 if none within the caliper)."""
    perm = rng.permutation(len(pool_score))
    srt = perm[np.argsort(pool_score[perm], kind="stable")]
    vals = pool_score[srt]
    n = len(vals)
    nf = NextFree(n)
    out = np.full(len(pos_score), -1, dtype=np.int64)
    js = np.searchsorted(vals, pos_score, side="left")
    for k in rng.permutation(len(pos_score)):
        s = pos_score[k]
        j = int(js[k])
        r = nf.right(j)
        l = nf.left(j - 1) if j >= 1 else -1
        dr = abs(vals[r] - s) if r < n else np.inf
        dl = abs(vals[l] - s) if l >= 0 else np.inf
        if dr == np.inf and dl == np.inf:
            continue
        if dr < dl:
            pick, d = r, dr
        elif dl < dr:
            pick, d = l, dl
        else:
            pick, d = (r, dr) if rng.integers(0, 2) == 0 else (l, dl)
        if caliper is not None and d > caliper:
            continue
        nf.take(pick)
        out[k] = srt[pick]
    return out


def match_design(label, protein, pos_mask, score, rng, caliper=None, within_protein=False):
    """role array (1 positive, 0 matched control, -1 excluded) and a small info dict. The pool is label == 0."""
    role = np.full(len(label), -1, dtype=np.int64)
    pool_idx = np.flatnonzero(label == 0)
    pos_idx = np.flatnonzero(pos_mask)
    if not within_protein:
        m = nn_match(score[pos_idx], score[pool_idx], rng, caliper)
        ok = m >= 0
        role[pos_idx[ok]] = 1
        role[pool_idx[m[ok]]] = 0
    else:
        by_p, pos_by_p = {}, {}
        for i in pool_idx:
            by_p.setdefault(protein[i], []).append(i)
        for i in pos_idx:
            pos_by_p.setdefault(protein[i], []).append(i)
        for p in sorted(pos_by_p):
            pi = np.array(pos_by_p[p])
            ci = np.array(by_p.get(p, []), dtype=np.int64)
            if ci.size == 0:
                continue
            m = nn_match(score[pi], score[ci], rng, caliper)
            ok = m >= 0
            role[pi[ok]] = 1
            role[ci[m[ok]]] = 0
    info = {"n_positive_input": int(pos_mask.sum()), "n_positive_matched": int((role == 1).sum()),
            "n_control": int((role == 0).sum())}
    return role, info
