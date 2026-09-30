"""Artefact 1 - protease cleavage geometry.

Question: is being a reported site associated with having a cleavage-competent residue near the
cysteine? Detectability depends on where the protease cuts (peptide length, termini), so a
positive-vs-background difference in cleavage-site proximity is a signature that the positive set
carries digestion geometry. A basic-residue "motif" around trypsin sites is the classic symptom.

Statistic: for each band in constants.CLEAVAGE_BANDS, the Haldane log2 odds ratio of
"a cleavage-competent residue within the band" between positives and background, with a
cluster-bootstrap interval at the Bonferroni level (two bands). Null 0, margin MMD_LOG2_OR.
The test status is the most severe band status (ties broken by the larger |estimate|).
Requires --fasta. Run it with the protease the data were digested with; running it again with a
different protease's rule is how the tool asks whether a signal follows the protease (an
orthogonal check, not a proof).
"""
from __future__ import annotations

import numpy as np

from .. import proteases, stats
from ..constants import (BOOTSTRAP_REPS, CLEAVAGE_BANDS, CLEAVAGE_CI_LEVEL, MMD_LOG2_OR, SEVERITY)
from ..io import background_mask
from ..status import status_from_interval
from .common import rnd, too_small, undecidable

TEST = "cleavage_geometry"


def band_flags(ds, protease):
    masks = {}
    flags = {b: np.zeros(len(ds), dtype=bool) for b in CLEAVAGE_BANDS}
    has_seq = np.zeros(len(ds), dtype=bool)
    for i, (p, pos) in enumerate(zip(ds.protein, ds.position)):
        s = ds.fasta.get(p)
        if s is None or pos > len(s):
            continue
        has_seq[i] = True
        if p not in masks:
            masks[p] = proteases.competent_mask(s, protease)
        for b, band in CLEAVAGE_BANDS.items():
            flags[b][i] = proteases.band_flag(s, int(pos), band, protease, masks[p])
    return flags, has_seq


def run(ds, cfg):
    if ds.fasta is None:
        return undecidable(TEST, "needs --fasta to locate cleavage sites")
    bg = background_mask(ds, cfg["background"])
    if bg is None:
        return undecidable(TEST, "background 'observed' needs a 'detected' column")
    flags, has_seq = band_flags(ds, cfg["protease"])
    pos = (ds.label == 1) & has_seq
    bgm = bg & has_seq
    msg = too_small(int(pos.sum()), int(bgm.sum()))
    if msg:
        return undecidable(TEST, msg, n_positive=int(pos.sum()), n_background=int(bgm.sum()))
    keep = pos | bgm
    codes, labels = stats.cluster_index([c for c, k in zip(ds.cluster, keep) if k])
    k_idx = np.flatnonzero(keep)
    bands = {}
    for j, (b, band) in enumerate(CLEAVAGE_BANDS.items()):
        counts = stats.two_by_two_counts(codes, len(labels), pos[k_idx], flags[b][k_idx])
        point, (lo, hi), p, _ = stats.boot_log2_or(counts, cfg.get("reps", BOOTSTRAP_REPS),
                                                   cfg["seed"] + 101 + j, CLEAVAGE_CI_LEVEL)
        st, why = status_from_interval(point, lo, hi, 0.0, MMD_LOG2_OR)
        tot = counts.sum(axis=0)
        bands[b] = {"band_residues": list(band), "status": st, "reason": why, "estimate": rnd(point),
                    "ci": [rnd(lo), rnd(hi)], "p_value": rnd(p, 8),
                    "share_flagged_positive": rnd(tot[0] / (tot[0] + tot[1])),
                    "share_flagged_background": rnd(tot[2] / (tot[2] + tot[3]))}
    head = max(bands, key=lambda b: (SEVERITY[bands[b]["status"]], abs(bands[b]["estimate"])))
    h = bands[head]
    return {"test": TEST, "status": h["status"], "reason": f"band {head}: {h['reason']}",
            "statistic": "log2_odds_ratio(cleavage-competent residue in band; positive vs background)",
            "null": 0.0, "margin": MMD_LOG2_OR, "ci_level": CLEAVAGE_CI_LEVEL, "estimate": h["estimate"],
            "ci": h["ci"], "p_value": h["p_value"], "n_positive": int(pos.sum()), "n_background": int(bgm.sum()),
            "n_clusters": len(labels),
            "details": {"protease": cfg["protease"], "background": cfg["background"], "headline_band": head,
                        "bands": bands, "rows_without_sequence": int((~has_seq).sum())}}
