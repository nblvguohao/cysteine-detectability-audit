"""Input schema, readers and validation.

One row per cysteine that the analysis considers. Positive rows are the reported modified sites;
every other row is a candidate background cysteine. The same table carries both, so the choice of
background (observed vs proteome) is a filter, not a second file.
"""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass, field

import numpy as np

SCHEMA = {
    "protein":      ("required", "str",   "protein accession; must match a FASTA identifier if --fasta is given"),
    "position":     ("required", "int",   "1-based residue index of the cysteine in the protein sequence"),
    "label":        ("required", "0/1",   "1 = reported modified site; 0 = candidate background cysteine"),
    "detected":     ("optional", "0/1",   "1 = the cysteine was observed in the run in any form (modified or not)"),
    "n_cys_peptide": ("optional", "int",  "smallest number of cysteines on any peptide that observed this cysteine"),
    "abundance":    ("optional", "float", "protein abundance (>0), identical on all rows of a protein; blank = missing"),
    "peptide_mass": ("optional", "float", "monoisotopic mass (Da) of the peptide supporting a positive row"),
    "cluster":      ("optional", "str",   "resampling cluster (e.g. homology component); default = protein"),
}


class InputError(ValueError):
    pass


def read_fasta(path):
    """Returns {identifier: sequence}. Identifier = second |-field of a UniProt header, else the first
    whitespace token. Both the UniProt accession and the full first token are indexed."""
    seqs, cur, buf = {}, None, []

    def flush():
        if cur is not None:
            s = "".join(buf).upper().rstrip("*")
            for key in cur:
                seqs.setdefault(key, s)

    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                flush()
                tok = line[1:].split()[0] if line[1:].split() else ""
                keys = [tok]
                parts = tok.split("|")
                if len(parts) >= 3:
                    keys.insert(0, parts[1])
                cur, buf = keys, []
            elif line:
                buf.append(line.strip())
        flush()
    return seqs


def _to_int01(v, col, i):
    if v in ("0", "1"):
        return int(v)
    raise InputError(f"row {i}: column '{col}' must be 0 or 1, got {v!r}")


def _to_float(v):
    if v is None or v.strip() == "" or v.strip().lower() in ("nan", "na"):
        return float("nan")
    return float(v)


@dataclass
class Dataset:
    protein: list
    position: np.ndarray
    label: np.ndarray
    detected: np.ndarray | None
    n_cys_peptide: np.ndarray | None
    abundance: np.ndarray | None
    peptide_mass: np.ndarray | None
    cluster: list
    extra: dict = field(default_factory=dict)
    fasta: dict | None = None
    notes: list = field(default_factory=list)

    def __len__(self):
        return len(self.protein)

    def subset(self, mask):
        mask = np.asarray(mask, dtype=bool)
        idx = np.flatnonzero(mask)
        pick = lambda a: None if a is None else a[mask]
        return Dataset(protein=[self.protein[i] for i in idx], position=self.position[mask],
                       label=self.label[mask], detected=pick(self.detected),
                       n_cys_peptide=pick(self.n_cys_peptide), abundance=pick(self.abundance),
                       peptide_mass=pick(self.peptide_mass), cluster=[self.cluster[i] for i in idx],
                       extra={k: v[mask] for k, v in self.extra.items()}, fasta=self.fasta,
                       notes=list(self.notes))


def read_sites(path, fasta=None, expand_background=False, max_residue_mismatch=0.01):
    """Read and validate the site table. See SCHEMA. Extra columns are kept (as float if every value
    parses, else dropped with a note) so a claim feature can be named on the command line."""
    with open(path, encoding="utf-8-sig", newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = rd.fieldnames or []
        missing = [c for c, (req, _, _) in SCHEMA.items() if req == "required" and c not in cols]
        if missing:
            raise InputError(f"missing required column(s): {missing}")
        rows = list(rd)
    if not rows:
        raise InputError("the site table has no rows")
    seen = set()
    protein, position, label = [], [], []
    for i, r in enumerate(rows, 2):
        p = r["protein"].strip()
        try:
            pos = int(r["position"])
        except ValueError:
            raise InputError(f"row {i}: position {r['position']!r} is not an integer")
        if pos < 1:
            raise InputError(f"row {i}: position must be >= 1")
        if (p, pos) in seen:
            raise InputError(f"row {i}: duplicate (protein, position) = ({p}, {pos})")
        seen.add((p, pos))
        protein.append(p)
        position.append(pos)
        label.append(_to_int01(r["label"].strip(), "label", i))
    n = len(rows)
    detected = (np.array([_to_int01(r["detected"].strip(), "detected", i) for i, r in enumerate(rows, 2)])
                if "detected" in cols else None)
    ncys = None
    if "n_cys_peptide" in cols:
        vals = []
        for r in rows:
            v = r["n_cys_peptide"].strip()
            vals.append(float(v) if v else float("nan"))
        ncys = np.array(vals)
    abundance = np.array([_to_float(r["abundance"]) for r in rows]) if "abundance" in cols else None
    pmass = np.array([_to_float(r["peptide_mass"]) for r in rows]) if "peptide_mass" in cols else None
    cluster = [r["cluster"].strip() or protein[k] for k, r in enumerate(rows)] if "cluster" in cols else list(protein)
    notes = []
    if abundance is not None:
        by_prot = {}
        for p, a in zip(protein, abundance):
            if math.isfinite(a):
                if a <= 0:
                    raise InputError(f"abundance must be > 0 (protein {p})")
                by_prot.setdefault(p, set()).add(a)
        bad = [p for p, s in by_prot.items() if len(s) > 1]
        if bad:
            raise InputError(f"abundance differs between rows of the same protein, e.g. {bad[:3]}")
    extra = {}
    for c in cols:
        if c in SCHEMA:
            continue
        try:
            extra[c] = np.array([_to_float(r[c]) for r in rows])
        except ValueError:
            notes.append(f"extra column '{c}' is not numeric and was ignored")
    ds = Dataset(protein=protein, position=np.array(position), label=np.array(label), detected=detected,
                 n_cys_peptide=ncys, abundance=abundance, peptide_mass=pmass, cluster=cluster,
                 extra=extra, fasta=None, notes=notes)
    if detected is not None:
        undetected_pos = int(((ds.label == 1) & (ds.detected == 0)).sum())
        if undetected_pos:
            notes.append(f"{undetected_pos} positive rows are marked detected=0; they are kept as positives")
    if fasta is not None:
        seqs = read_fasta(fasta)
        ds.fasta = seqs
        absent = sum(1 for p in protein if p not in seqs)
        wrong = sum(1 for p, pos in zip(protein, position)
                    if p in seqs and (pos > len(seqs[p]) or seqs[p][pos - 1] != "C"))
        if absent:
            notes.append(f"{absent} rows name a protein that is not in the FASTA")
        checked = n - absent
        if checked and wrong / checked > max_residue_mismatch:
            raise InputError(f"{wrong} of {checked} rows do not point at a cysteine in the FASTA "
                             f"(limit {max_residue_mismatch:.0%}); wrong FASTA version or 0-based positions?")
        if wrong:
            notes.append(f"{wrong} rows do not point at a cysteine in the FASTA (within the {max_residue_mismatch:.0%} tolerance)")
        if expand_background:
            ds = _expand(ds, seqs)
    elif expand_background:
        raise InputError("--expand-background needs --fasta")
    return ds


def _expand(ds, seqs):
    """Add every cysteine of every protein already in the table as a label-0, detected-0 row."""
    have = set(zip(ds.protein, ds.position.tolist()))
    add_p, add_pos, add_cl, add_ab = [], [], [], []
    ab_by = {}
    if ds.abundance is not None:
        for p, a in zip(ds.protein, ds.abundance):
            if math.isfinite(a):
                ab_by[p] = a
    cl_by = dict(zip(ds.protein, ds.cluster))
    for p in sorted(set(ds.protein)):
        s = seqs.get(p)
        if s is None:
            continue
        for i, ch in enumerate(s, 1):
            if ch == "C" and (p, i) not in have:
                add_p.append(p)
                add_pos.append(i)
                add_cl.append(cl_by[p])
                add_ab.append(ab_by.get(p, float("nan")))
    k = len(add_p)
    if k == 0:
        return ds
    nanarr = np.full(k, np.nan)
    out = Dataset(
        protein=ds.protein + add_p,
        position=np.concatenate([ds.position, np.array(add_pos, dtype=ds.position.dtype)]),
        label=np.concatenate([ds.label, np.zeros(k, dtype=ds.label.dtype)]),
        detected=None if ds.detected is None else np.concatenate([ds.detected, np.zeros(k, dtype=ds.detected.dtype)]),
        n_cys_peptide=None if ds.n_cys_peptide is None else np.concatenate([ds.n_cys_peptide, nanarr]),
        abundance=None if ds.abundance is None else np.concatenate([ds.abundance, np.array(add_ab)]),
        peptide_mass=None if ds.peptide_mass is None else np.concatenate([ds.peptide_mass, nanarr]),
        cluster=ds.cluster + add_cl,
        extra={c: np.concatenate([v, nanarr]) for c, v in ds.extra.items()},
        fasta=ds.fasta, notes=ds.notes + [f"--expand-background added {k} unlisted cysteines as label 0, detected 0"])
    return out


def background_mask(ds, background):
    """Rows that act as background under the chosen definition ('observed' or 'proteome')."""
    if background == "proteome":
        return ds.label == 0
    if background == "observed":
        if ds.detected is None:
            return None
        return (ds.label == 0) & (ds.detected == 1)
    raise InputError(f"unknown background {background!r}; use 'observed' or 'proteome'")
