"""Artefact 3 - protein abundance.

Question: does protein abundance separate positives from background? Abundant proteins are sampled
more deeply, so their cysteines are more often seen and more often called; abundance then rides
along with any protein-level property it correlates with.

Statistic: tie-corrected rank AUC of log abundance (a protein-level value copied to its rows),
positives vs background under --background, clustered bootstrap, 95%. Null 0.5, margin MMD_AUC.
Rows without an abundance value are excluded; if fewer than ABUNDANCE_MIN_COVERAGE of the tested
rows carry one, the test is UNDECIDABLE (the missing rows are rarely missing at random: low-abundance
proteins are the ones an external quantification fails to see).
"""
from __future__ import annotations

import numpy as np

from .. import stats
from ..constants import ABUNDANCE_MIN_COVERAGE, BOOTSTRAP_REPS, CI_LEVEL, MMD_AUC
from ..io import background_mask
from ..status import status_from_interval
from .common import rnd, too_small, undecidable

TEST = "protein_abundance"


def run(ds, cfg):
    if ds.abundance is None:
        return undecidable(TEST, "needs an 'abundance' column")
    bg = background_mask(ds, cfg["background"])
    if bg is None:
        return undecidable(TEST, "background 'observed' needs a 'detected' column")
    tested = (ds.label == 1) | bg
    have = np.isfinite(ds.abundance)
    coverage = float((tested & have).sum() / tested.sum()) if tested.sum() else 0.0
    if coverage < ABUNDANCE_MIN_COVERAGE:
        return undecidable(TEST, f"abundance present on {coverage:.3f} of tested rows < {ABUNDANCE_MIN_COVERAGE}",
                           details={"coverage": rnd(coverage)})
    pos = (ds.label == 1) & have
    bgm = bg & have
    msg = too_small(int(pos.sum()), int(bgm.sum()))
    if msg:
        return undecidable(TEST, msg, n_positive=int(pos.sum()), n_background=int(bgm.sum()))
    k_idx = np.flatnonzero(pos | bgm)
    codes, labels = stats.cluster_index([ds.cluster[i] for i in k_idx])
    score = np.log(ds.abundance[k_idx])
    point, (lo, hi), p, _ = stats.boot_auc(score, pos[k_idx], codes, len(labels),
                                           cfg.get("reps", BOOTSTRAP_REPS), cfg["seed"] + 301, CI_LEVEL)
    st, why = status_from_interval(point, lo, hi, 0.5, MMD_AUC)
    return {"test": TEST, "status": st, "reason": why,
            "statistic": "rank AUC of log abundance (positive vs background)", "null": 0.5, "margin": MMD_AUC,
            "ci_level": CI_LEVEL, "estimate": rnd(point), "ci": [rnd(lo), rnd(hi)], "p_value": rnd(p, 8),
            "n_positive": int(pos.sum()), "n_background": int(bgm.sum()), "n_clusters": len(labels),
            "details": {"background": cfg["background"], "coverage": rnd(coverage),
                        "median_log_abundance_positive": rnd(float(np.median(np.log(ds.abundance[pos])))),
                        "median_log_abundance_background": rnd(float(np.median(np.log(ds.abundance[bgm]))))}}
