"""Synthetic inputs with a planted artefact of known strength (0 = none, 1 = maximal).

Used by the unit tests, the demo, and the calibration experiment. Proteins are random sequences with a
fixed residue composition; every cysteine becomes a row. Labels are drawn so that at strength 0 they are
independent of the artefact axis, and at strength s a share s of positives is drawn from the artefact
pool (the same mixture design as the analysis tree's Phase 2 benchmark).

  cleavage   pool = cysteines with a trypsin-competent residue in the proximal band (1-3); background
             = proteome; detected column is all 1 (not the tested axis)
  abundance  protein abundance ~ lognormal; the positive probability of a cysteine rises with its
             protein's abundance rank with weight s; background = proteome
  multicys   each cysteine gets n_cys_peptide in {1, 2, 3} (P(>=2) = 0.3); pool = n_cys_peptide >= 2;
             all rows detected; background = observed
  overlap    true modified prevalence 0.2; a share s of the unmodified cysteines is removed from the
             detected set (they become detected 0), as in Phase 2 benchmark B4
  none       labels independent of everything (a null dataset for every axis)
"""
from __future__ import annotations

import csv
import os

import numpy as np

from . import proteases
from .constants import CLEAVAGE_BANDS

AA = "ACDEFGHIKLMNPQRSTVWY"
# composition roughly vertebrate-like; C ~2.3%, K+R ~11%
FREQ = np.array([7.0, 2.3, 4.8, 7.0, 3.7, 6.6, 2.6, 4.4, 5.7, 9.9, 2.2, 3.6, 6.2, 4.7, 5.6, 8.2, 5.3, 6.1, 1.2, 2.7])
FREQ = FREQ / FREQ.sum()


def random_proteome(n, rng, lo=200, hi=800):
    seqs = {}
    for i in range(n):
        L = int(rng.integers(lo, hi + 1))
        seqs[f"SYN{i:05d}"] = "M" + "".join(rng.choice(list(AA), size=L - 1, p=FREQ))
    return seqs


def _cys_rows(seqs):
    prot, pos = [], []
    for p in sorted(seqs):
        for i, ch in enumerate(seqs[p], 1):
            if ch == "C":
                prot.append(p)
                pos.append(i)
    return prot, np.array(pos)


def _mixture_labels(n, pool_mask, frac_pos, s, rng):
    """Draw round(frac_pos*n) positives without replacement: each from the pool w.p. s, else from all."""
    n_pos = int(round(frac_pos * n))
    label = np.zeros(n, dtype=int)
    pool = list(np.flatnonzero(pool_mask))
    allidx = list(range(n))
    rng.shuffle(pool)
    rng.shuffle(allidx)
    chosen = set()
    pi = ai = 0
    for _ in range(n_pos):
        take_pool = rng.random() < s and pi < len(pool)
        if take_pool:
            while pi < len(pool) and pool[pi] in chosen:
                pi += 1
            if pi < len(pool):
                chosen.add(pool[pi])
                continue
        while ai < n and allidx[ai] in chosen:
            ai += 1
        if ai < n:
            chosen.add(allidx[ai])
    label[list(chosen)] = 1
    return label


def generate(artefact, strength, seed, n_proteins=600, frac_pos=0.2):
    rng = np.random.default_rng(seed)
    seqs = random_proteome(n_proteins, rng)
    prot, pos = _cys_rows(seqs)
    n = len(prot)
    cols = {"protein": prot, "position": pos}
    s = float(strength)
    if artefact == "cleavage":
        band = CLEAVAGE_BANDS["proximal_1_3"]
        masks = {p: proteases.competent_mask(seqs[p], "trypsin") for p in seqs}
        pool = np.array([proteases.band_flag(seqs[p], int(q), band, "trypsin", masks[p]) for p, q in zip(prot, pos)])
        cols["label"] = _mixture_labels(n, pool, frac_pos, s, rng)
        cols["detected"] = np.ones(n, dtype=int)
    elif artefact == "abundance":
        ab_prot = {p: float(rng.lognormal(0.0, 1.2)) for p in sorted(seqs)}
        ab = np.array([ab_prot[p] for p in prot])
        rank_q = (np.argsort(np.argsort(ab)) + 1) / n
        w = (1 - s) + s * 2 * rank_q
        prob = np.clip(frac_pos * w, 0, 1)
        cols["label"] = (rng.random(n) < prob).astype(int)
        cols["detected"] = np.ones(n, dtype=int)
        cols["abundance"] = ab
    elif artefact == "multicys":
        ncys = rng.choice([1, 2, 3], size=n, p=[0.7, 0.2, 0.1])
        cols["label"] = _mixture_labels(n, ncys >= 2, frac_pos, s, rng)
        cols["detected"] = np.ones(n, dtype=int)
        cols["n_cys_peptide"] = ncys
    elif artefact == "overlap":
        is_mod = rng.random(n) < frac_pos
        survive = is_mod | (rng.random(n) >= s)
        cols["label"] = is_mod.astype(int)
        cols["detected"] = survive.astype(int)
    elif artefact == "none":
        cols["label"] = (rng.random(n) < frac_pos).astype(int)
        cols["detected"] = np.ones(n, dtype=int)
        ab_prot = {p: float(rng.lognormal(0.0, 1.2)) for p in sorted(seqs)}
        cols["abundance"] = np.array([ab_prot[p] for p in prot])
        cols["n_cys_peptide"] = rng.choice([1, 2, 3], size=n, p=[0.7, 0.2, 0.1])
    else:
        raise ValueError(f"unknown artefact {artefact!r}")
    return seqs, cols


def write(seqs, cols, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "proteome.fasta"), "w", encoding="utf-8") as fh:
        for p in sorted(seqs):
            fh.write(f">{p}\n")
            s = seqs[p]
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")
    names = list(cols)
    n = len(cols["protein"])
    with open(os.path.join(out_dir, "sites.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(names)
        for i in range(n):
            row = []
            for c in names:
                v = cols[c][i]
                row.append(f"{v:.6g}" if isinstance(v, (float, np.floating)) else str(v))
            w.writerow(row)


def simulate(artefact, strength, seed, out_dir, n_proteins=600):
    seqs, cols = generate(artefact, strength, seed, n_proteins=n_proteins)
    write(seqs, cols, out_dir)
    return seqs, cols
