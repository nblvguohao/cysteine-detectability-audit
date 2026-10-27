"""Protease cleavage rules.

Each rule says which residues are cleavage-competent, on which side of the residue the peptide bond
is cut, and whether a following proline blocks the cut. The cleavage-geometry test asks whether a
cleavage-competent residue sits within a band of residues from the cysteine, so only the residue
set and the proline rule matter there; the side matters for in-silico digestion (not used by the
tests, provided for users who want it).
"""
from __future__ import annotations

RULES = {
    # name: (residues, side, blocked_by_following_proline)
    "trypsin":       ("KR", "C", True),
    "trypsin/p":     ("KR", "C", False),
    "lysc":          ("K", "C", False),
    "argc":          ("R", "C", True),
    "aspn":          ("D", "N", False),
    "gluc":          ("E", "C", False),
    "gluc_de":       ("DE", "C", False),
    "chymotrypsin":  ("FWY", "C", True),
    "chymotrypsin+": ("FWYLM", "C", True),
}


def rule(name):
    key = name.lower()
    if key not in RULES:
        raise ValueError(f"unknown protease {name!r}; choose one of {sorted(RULES)}")
    return RULES[key]


def competent_mask(seq, name):
    """Boolean list: residue i (0-based) is a cleavage-competent residue under the rule."""
    residues, side, block_p = rule(name)
    n = len(seq)
    out = [False] * n
    for i, ch in enumerate(seq):
        if ch in residues:
            if block_p and side == "C" and i + 1 < n and seq[i + 1] == "P":
                continue
            out[i] = True
    return out


def band_flag(seq, pos1, band, name, mask=None):
    """True if a cleavage-competent residue lies at a distance in [band[0], band[1]] residues from
    the residue at 1-based position pos1, on either side."""
    if mask is None:
        mask = competent_mask(seq, name)
    i0 = pos1 - 1
    lo, hi = band
    n = len(seq)
    for d in range(lo, hi + 1):
        for j in (i0 - d, i0 + d):
            if 0 <= j < n and mask[j]:
                return True
    return False


def digest(seq, name, missed=2, min_len=7, max_len=30):
    """In-silico digestion -> list of (start1, end1) peptides (1-based inclusive)."""
    residues, side, _ = rule(name)
    mask = competent_mask(seq, name)
    cuts = [0]
    for i, m in enumerate(mask):
        if m:
            b = i + 1 if side == "C" else i
            if 0 < b < len(seq):
                cuts.append(b)
    cuts.append(len(seq))
    cuts = sorted(set(cuts))
    peps = []
    for a in range(len(cuts) - 1):
        for k in range(missed + 1):
            b = a + 1 + k
            if b >= len(cuts):
                break
            s, e = cuts[a], cuts[b]
            if min_len <= e - s <= max_len:
                peps.append((s + 1, e))
    return peps
