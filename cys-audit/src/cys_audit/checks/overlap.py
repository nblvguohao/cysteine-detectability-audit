"""Artefact 4 - positive set = detected set (negative-class depletion).

Question: among cysteines the run observed at all, what share are reported as modified? When an
enrichment is specific, nearly every observed cysteine is a positive, the "detected but unmodified"
class nearly vanishes, and site-preference claims cannot be tested against an observed background -
only against a proteome background, which reintroduces detectability.

Statistic: coincidence = #(positive & detected) / #(detected), clustered bootstrap 95% interval
(Wilson interval also reported). Status on the bootstrap interval:
  FAIL        lower bound >= OVERLAP_FAIL_LOWER (0.90)
  WARNING     lower bound >= OVERLAP_WARNING_LOWER (0.50)
  PASS        upper bound <  OVERLAP_PASS_UPPER (0.50)
  UNDECIDABLE otherwise, or no 'detected' column
The count of detected-unmodified cysteines is reported because every observed-background test in this
audit depends on it.
"""
from __future__ import annotations

import numpy as np

from .. import stats
from ..constants import (BOOTSTRAP_REPS, CI_LEVEL, FAIL, MIN_BACKGROUND, OVERLAP_FAIL_LOWER, OVERLAP_PASS_UPPER,
                         OVERLAP_WARNING_LOWER, PASS, UNDECIDABLE, WARNING)
from .common import rnd, undecidable

TEST = "positive_detected_overlap"


def run(ds, cfg):
    if ds.detected is None:
        return undecidable(TEST, "needs a 'detected' column")
    det = ds.detected == 1
    n = int(det.sum())
    if n == 0:
        return undecidable(TEST, "no detected rows")
    k_idx = np.flatnonzero(det)
    codes, labels = stats.cluster_index([ds.cluster[i] for i in k_idx])
    succ = ds.label[k_idx] == 1
    k = int(succ.sum())
    point, (lo, hi), _ = stats.boot_proportion(codes, len(labels), succ, cfg.get("reps", BOOTSTRAP_REPS),
                                               cfg["seed"] + 401, CI_LEVEL)
    wl, wh = stats.wilson(k, n)
    if lo >= OVERLAP_FAIL_LOWER:
        st, why = FAIL, f"lower bound {lo:.4f} >= {OVERLAP_FAIL_LOWER}: the observed-unmodified class is nearly absent"
    elif lo >= OVERLAP_WARNING_LOWER:
        st, why = WARNING, f"lower bound {lo:.4f} >= {OVERLAP_WARNING_LOWER}: positives dominate the observed cysteines"
    elif hi < OVERLAP_PASS_UPPER:
        st, why = PASS, f"upper bound {hi:.4f} < {OVERLAP_PASS_UPPER}: an observed-unmodified class exists"
    else:
        st, why = UNDECIDABLE, f"interval [{lo:.4f}, {hi:.4f}] straddles {OVERLAP_PASS_UPPER}"
    n_unmod = n - k
    return {"test": TEST, "status": st, "reason": why,
            "statistic": "coincidence = positives / detected cysteines", "null": None,
            "margin": {"fail_lower": OVERLAP_FAIL_LOWER, "warning_lower": OVERLAP_WARNING_LOWER,
                       "pass_upper": OVERLAP_PASS_UPPER},
            "ci_level": CI_LEVEL, "estimate": rnd(point), "ci": [rnd(lo), rnd(hi)], "p_value": None,
            "n_positive": k, "n_background": n_unmod, "n_clusters": len(labels),
            "details": {"detected": n, "positive_and_detected": k, "detected_unmodified": n_unmod,
                        "wilson": [rnd(wl), rnd(wh)],
                        "observed_background_usable": bool(n_unmod >= MIN_BACKGROUND),
                        "positives_not_detected": int(((ds.label == 1) & ~det).sum())}}
