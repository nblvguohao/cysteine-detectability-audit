"""Shared helpers for item H_artifact5_docs (POST HOC revision analysis, 2026-09-30).

Everything produced by this item is a POST HOC revision analysis written in response to a
pre-submission review. Nothing here is registered or pre-specified.

The helpers only read inputs (never write outside results/H_artifact5_docs) and keep a registry of
every input file read, with its sha256, so that provenance.json can list them.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import re

W = r"/path/to/local/_cys_repo_work/public/revision_2026-09-30"
REPO = r"/path/to/local/_cys_repo_work/repo"
MCP = r"/path/to/local/巯基化/MCP"
CYS_AUDIT_SRC = r"/path/to/local/_cys_repo_work/public/cys-audit/src"
ITEM = "H_artifact5_docs"
OUT = os.path.join(W, "results", ITEM)
SCRIPTS = os.path.join(W, "scripts", ITEM)

POSTHOC_LABEL = "POST HOC revision analysis (2026-09-30), not registered, not pre-specified"

# monoisotopic masses (same residue table as the repository's audit script)
AA = {'G': 57.02146, 'A': 71.03711, 'S': 87.03203, 'P': 97.05276, 'V': 99.06841, 'T': 101.04768,
      'C': 103.00919, 'L': 113.08406, 'I': 113.08406, 'N': 114.04293, 'D': 115.02694,
      'Q': 128.05858, 'K': 128.09496, 'E': 129.04259, 'M': 131.04049, 'H': 137.05891,
      'F': 147.06841, 'R': 156.10111, 'Y': 163.06333, 'W': 186.07931}
H2O = 18.010565
SULFIDE = 31.972071
DIOXIDATION = 31.989829
SEPARATION_DA = 0.017758          # dioxidation minus sulfide, as used throughout the manuscript
IODOTMT = 329.226595
CAM = 57.021464
OX_M = 15.9949
NEM = 125.047679

_INPUTS: dict[str, dict] = {}


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def sha256_lf(path: str) -> str:
    """sha256 after CRLF -> LF normalisation (the Windows checkout converts line endings)."""
    with open(path, "rb") as fh:
        data = fh.read().replace(b"\r\n", b"\n")
    return hashlib.sha256(data).hexdigest()


def register(path: str, role: str, text: bool = False) -> str:
    """Record an input file (sha256 of the bytes on disk; for text files also the LF-normalised hash)."""
    path = os.path.normpath(path)
    if path not in _INPUTS:
        rec = {"path": path.replace("\\", "/"), "role": role, "bytes": os.path.getsize(path),
               "sha256": sha256_file(path)}
        if text:
            rec["sha256_lf_normalised"] = sha256_lf(path)
        _INPUTS[path] = rec
    return path


def inputs_registry() -> list[dict]:
    return sorted(_INPUTS.values(), key=lambda r: r["path"])


def read_json(path: str, role: str):
    register(path, role, text=True)
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def read_text(path: str, role: str) -> str:
    register(path, role, text=True)
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def write_json(obj, name: str) -> str:
    path = os.path.join(OUT, name)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, ensure_ascii=False, indent=2, default=str)
    os.replace(tmp, path)
    return path


def write_csv(df, name: str) -> str:
    path = os.path.join(OUT, name)
    tmp = path + ".tmp"
    df.to_csv(tmp, index=False, encoding="utf-8", lineterminator="\n")
    os.replace(tmp, path)
    return path


def peptide_mass(seq: str) -> float:
    return sum(AA[c] for c in seq) + H2O


def read_fasta_gz(path: str, role: str) -> list[tuple[str, str]]:
    register(path, role)
    out, head, cur = [], None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if head is not None:
                    out.append((head, "".join(cur)))
                head, cur = line[1:], []
            else:
                cur.append(line.strip())
    if head is not None:
        out.append((head, "".join(cur)))
    return out


def tryptic_peptides(seq: str, min_len=7, max_len=30, max_mc=2):
    """Fully tryptic peptides (cleave after K/R, not before P), 0..max_mc missed cleavages."""
    sites = [i + 1 for i, c in enumerate(seq[:-1]) if c in "KR" and seq[i + 1] != "P"]
    bounds = [0] + sites + [len(seq)]
    for i in range(len(bounds) - 1):
        for mc in range(max_mc + 1):
            j = i + 1 + mc
            if j >= len(bounds):
                break
            pep = seq[bounds[i]:bounds[j]]
            if min_len <= len(pep) <= max_len:
                yield pep, mc


def comet_decoy(peptide: str) -> str:
    """Comet internal decoy as inferred from the output here: reverse all but the C-terminal residue."""
    return peptide[:-1][::-1] + peptide[-1]


def strip_latex(s: str) -> str:
    s = re.sub(r"\\textbf\{([^}]*)\}", r"\1", s)
    s = re.sub(r"\\emph\{([^}]*)\}", r"\1", s)
    return s
