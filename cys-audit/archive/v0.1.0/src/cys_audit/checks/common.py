from __future__ import annotations

from ..constants import MIN_BACKGROUND, MIN_POSITIVES, UNDECIDABLE


def undecidable(test, reason, **extra):
    out = {"test": test, "status": UNDECIDABLE, "reason": reason, "statistic": None, "null": None,
           "margin": None, "estimate": None, "ci": None, "p_value": None, "n_positive": None,
           "n_background": None, "details": {}}
    out.update(extra)
    return out


def too_small(n_pos, n_bg):
    if n_pos < MIN_POSITIVES:
        return f"{n_pos} positives < minimum {MIN_POSITIVES}"
    if n_bg < MIN_BACKGROUND:
        return f"{n_bg} background rows < minimum {MIN_BACKGROUND}"
    return None


def rnd(x, k=6):
    return None if x is None else float(round(x, k))
