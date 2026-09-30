"""Artefact 2 - multi-cysteine peptides.

Question: are positives and observed background cysteines represented by multi-cysteine peptides to
a different extent? When they are, any filter or localisation rule applied to multi-cysteine
peptides acts asymmetrically on the two classes, and a neighbouring-cysteine or position feature can
change sign depending on that rule (in the analysis tree: -2.61 under asymmetric single-Cys filtering
against +3.34 with the filter applied to both classes).

Statistic: Haldane log2 odds ratio of "observed only on peptides with >= 2 cysteines"
(n_cys_peptide >= 2) between positives and OBSERVED background, clustered bootstrap, 95%.
Null 0, margin MMD_LOG2_OR. Needs n_cys_peptide and detected; undetected rows carry no peptide, so
this test always uses the observed background regardless of --background (reported in details).
"""
from __future__ import annotations

import numpy as np

from .. import stats
from ..constants import BOOTSTRAP_REPS, CI_LEVEL, MMD_LOG2_OR
from ..status import status_from_interval
from .common import rnd, too_small, undecidable

TEST = "multi_cysteine"


def run(ds, cfg):
    if ds.n_cys_peptide is None:
        return undecidable(TEST, "needs an 'n_cys_peptide' column")
    if ds.detected is None:
        return undecidable(TEST, "needs a 'detected' column (only observed cysteines have peptides)")
    have = np.isfinite(ds.n_cys_peptide)
    pos = (ds.label == 1) & have
    bg = (ds.label == 0) & (ds.detected == 1) & have
    msg = too_small(int(pos.sum()), int(bg.sum()))
    if msg:
        return undecidable(TEST, msg, n_positive=int(pos.sum()), n_background=int(bg.sum()))
    keep = pos | bg
    k_idx = np.flatnonzero(keep)
    codes, labels = stats.cluster_index([ds.cluster[i] for i in k_idx])
    flag = ds.n_cys_peptide[k_idx] >= 2
    counts = stats.two_by_two_counts(codes, len(labels), pos[k_idx], flag)
    point, (lo, hi), p, _ = stats.boot_log2_or(counts, cfg.get("reps", BOOTSTRAP_REPS), cfg["seed"] + 201, CI_LEVEL)
    st, why = status_from_interval(point, lo, hi, 0.0, MMD_LOG2_OR)
    tot = counts.sum(axis=0)
    return {"test": TEST, "status": st, "reason": why,
            "statistic": "log2_odds_ratio(multi-Cys-only support; positive vs observed background)",
            "null": 0.0, "margin": MMD_LOG2_OR, "ci_level": CI_LEVEL, "estimate": rnd(point),
            "ci": [rnd(lo), rnd(hi)], "p_value": rnd(p, 8), "n_positive": int(pos.sum()),
            "n_background": int(bg.sum()), "n_clusters": len(labels),
            "details": {"background_used": "observed (forced: undetected rows have no peptide)",
                        "share_multi_only_positive": rnd(tot[0] / (tot[0] + tot[1])),
                        "share_multi_only_background": rnd(tot[2] / (tot[2] + tot[3])),
                        "positives_without_n_cys_peptide": int(((ds.label == 1) & ~have).sum())}}
