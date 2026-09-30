"""The single rule that turns an effect estimate into PASS / WARNING / FAIL / UNDECIDABLE.

Given a point estimate, its interval [lo, hi], the null value and a minimal meaningful difference
(MMD, delta):

  UNDECIDABLE  the inputs the test needs are missing or a class is below the minimum size; or the
               interval covers the null AND reaches beyond +/- delta (a meaningful artefact can be
               neither shown nor excluded)
  PASS         the whole interval lies strictly inside (null - delta, null + delta): an artefact
               signature of meaningful size is excluded (equivalence, not merely "not significant")
  FAIL         the interval excludes the null AND lies entirely at or beyond delta on one side:
               an artefact signature of at least meaningful size is established
  WARNING      the interval excludes the null but does not establish a meaningful size (it
               straddles delta, or lies between the null and delta)

The rule is direction-agnostic: cleavage geometry, multi-cysteine composition and abundance can each
push an apparent preference either way, so an association of either sign is an artefact signature.
A p-value never changes a status.
"""
from __future__ import annotations

import math

from .constants import FAIL, PASS, UNDECIDABLE, WARNING


def status_from_interval(point, lo, hi, null, delta):
    if any(x is None or (isinstance(x, float) and not math.isfinite(x)) for x in (point, lo, hi)):
        return UNDECIDABLE, "estimate or interval not finite"
    excludes_null = lo > null or hi < null
    if not excludes_null:
        if lo > null - delta and hi < null + delta:
            return PASS, f"interval [{lo:.4g}, {hi:.4g}] lies inside the equivalence margin +/-{delta:g}"
        return UNDECIDABLE, (f"interval [{lo:.4g}, {hi:.4g}] covers the null {null:g} but reaches beyond "
                             f"+/-{delta:g}: a meaningful artefact can be neither shown nor excluded")
    if lo >= null + delta or hi <= null - delta:
        return FAIL, (f"interval [{lo:.4g}, {hi:.4g}] excludes the null and lies beyond the margin "
                      f"{delta:g}: artefact signature of meaningful size")
    return WARNING, (f"interval [{lo:.4g}, {hi:.4g}] excludes the null {null:g} but does not establish "
                     f"a size beyond the margin {delta:g}")
