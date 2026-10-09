"""Step 4c (POST HOC): the abundance ratio under the site-level labels.

The stored ratio (10.36-fold) uses label A, which counts a group as persulfidated when ANY peptide
(with or without a cysteine) of the persulfidation table has one of its accessions as leading razor
protein. 447 of those 1,197 groups have no identified Cys-containing peptide (s04). Here the same
ratio of medians and the same bootstrap (5,000 resamples over groups; seed 20260930) are computed
for: S (>=1 identified Cys-containing peptide), S1 (single-Cys peptide), and A-not-S (identified
only through peptides without cysteine), each against the groups not identified at all.
Output: s04c_label_ratios.csv
"""
import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import *  # noqa

F = pd.read_csv(RES / "s04_group_features.csv")
neg = F[F.A == 0]


def boot(pos_v, neg_v, n=5000, seed=SEED):
    rng = np.random.default_rng(seed)
    v = np.concatenate([pos_v, neg_v]); lab = np.r_[np.ones(len(pos_v), int), np.zeros(len(neg_v), int)]
    out = []
    for _ in range(n):
        s = rng.integers(0, len(v), len(v))
        vs, ls = v[s], lab[s]
        if ls.sum() == 0 or (1 - ls).sum() == 0:
            continue
        out.append(np.median(vs[ls == 1]) / np.median(vs[ls == 0]))
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


rows = []
for nm, pos in (("A_stored", F[F.A == 1]), ("S_any_Cys_peptide", F[F.S == 1]),
                ("S1_single_Cys_peptide", F[F.S1 == 1]), ("A_not_S_no_Cys_peptide", F[(F.A == 1) & (F.S == 0)])):
    for meas in ("ibaq", "quantity"):
        pv = pos[meas].dropna().values.astype(float); nv = neg[meas].dropna().values.astype(float)
        pv, nv = pv[pv > 0], nv[nv > 0]
        lo, hi = boot(pv, nv)
        rows.append({"positive_set": nm, "measure": meas, "n_pos": len(pv), "n_neg_not_identified": len(nv),
                     "median_ratio": float(np.median(pv) / np.median(nv)), "ci_low": lo, "ci_high": hi,
                     "bootstrap": "5,000 resamples over groups, seed 20260930"})
R = pd.DataFrame(rows)
R.to_csv(RES / "s04c_label_ratios.csv", index=False)
print(R.round(4).to_string())
