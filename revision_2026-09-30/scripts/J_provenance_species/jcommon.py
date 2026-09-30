# -*- coding: utf-8 -*-
"""Shared helpers for revision item J_provenance_species (POST HOC revision analysis, 2026-09-30).

Everything produced under this item is a post hoc revision analysis written in response to a
pre-submission review. Nothing here is registered or pre-specified.

Rules honoured by every script of this item
  * read-only on the manuscript folder, the internal repository, W/inputs and W/external;
  * writes only to W/scripts/J_provenance_species, W/results/J_provenance_species and W/reports;
  * never imports code from the internal repository (importing would write __pycache__ there);
    the few functions that must match the repository's instrument are re-implemented verbatim
    and cited by file and line;
  * every input that is read is hashed (sha256) and recorded in a per-script provenance fragment,
    which j9_provenance.py merges into results/J_provenance_species/provenance.json.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import platform
import sys

sys.dont_write_bytecode = True

W = "C:/Users/admin/Desktop/小论文/_cys_repo_work/public/revision_2026-09-30"
REPO = "C:/Users/admin/Desktop/小论文/_cys_repo_work/repo"
PURGED = "C:/Users/admin/Desktop/小论文/_cys_repo_work/PURGED_BACKUP"
MCP = "C:/Users/admin/Desktop/小论文/巯基化/MCP"
OUT = W + "/results/J_provenance_species"
SCRIPTS = W + "/scripts/J_provenance_species"
REPORTS = W + "/reports"

# the collaborator tree is not on this machine; this is the path recorded in the repository audits
COLLAB = "/Users/lyuguohao/Documents/zhanghua_huibao/07_共性规律_方向一"

POSTHOC_LABEL = ("POST HOC revision analysis (2026-09-30), written in response to a pre-submission "
                 "review; not registered, not pre-specified")


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


TEXT_EXT = (".py", ".csv", ".json", ".md", ".txt", ".tex", ".tsv")


def sha256_lf(path: str) -> str:
    """sha256 after converting CRLF to LF. The internal repository was checked out on Windows with
    line-ending conversion, so hashes recorded in the repository's own audits (made on macOS) match
    this value, not the raw working-tree hash."""
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read().replace(b"\r\n", b"\n")).hexdigest()


class Inputs:
    """Records every input file read, with sha256 (raw and, for text files, LF-normalised) and size."""

    def __init__(self):
        self.files = {}

    def use(self, path: str) -> str:
        p = os.path.normpath(path).replace("\\", "/")
        if p not in self.files:
            rec = {"sha256": sha256_file(p), "bytes": os.path.getsize(p)}
            if p.lower().endswith(TEXT_EXT):
                lf = sha256_lf(p)
                if lf != rec["sha256"]:
                    rec["sha256_lf_normalised"] = lf
            self.files[p] = rec
        return p

    def read_text(self, path: str, encoding: str = "utf-8") -> str:
        p = self.use(path)
        with open(p, encoding=encoding) as fh:
            return fh.read()

    def read_json(self, path: str):
        return json.loads(self.read_text(path))

    def read_csv(self, path: str, encoding: str = "utf-8-sig", delimiter: str = ","):
        p = self.use(path)
        with open(p, encoding=encoding, newline="") as fh:
            return list(csv.DictReader(fh, delimiter=delimiter))


def write_csv(path: str, rows: list, fieldnames: list | None = None) -> str:
    if fieldnames is None:
        fieldnames = []
        for r in rows:
            for k in r:
                if k not in fieldnames:
                    fieldnames.append(k)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="raise")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})
    return path


def write_json(path: str, obj) -> str:
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, ensure_ascii=False, indent=1, sort_keys=False)
        fh.write("\n")
    return path


def versions() -> dict:
    out = {"python": sys.version.split()[0], "platform": platform.platform()}
    for mod in ("numpy", "pandas", "scipy", "sklearn", "Bio"):
        try:
            m = __import__(mod)
            out[mod] = getattr(m, "__version__", "unknown")
        except Exception as exc:  # pragma: no cover
            out[mod] = f"not importable: {exc}"
    return out


def write_fragment(name: str, inputs: Inputs, outputs: list, seeds: dict, notes: dict | None = None):
    frag = {
        "script": name,
        "script_sha256": sha256_file(os.path.join(SCRIPTS, name)),
        "label": POSTHOC_LABEL,
        "inputs": inputs.files,
        "outputs": {os.path.relpath(p, W).replace("\\", "/"): sha256_file(p) for p in outputs},
        "seeds": seeds,
        "versions": versions(),
        "notes": notes or {},
    }
    write_json(os.path.join(OUT, f"_provenance_fragment_{os.path.splitext(name)[0]}.json"), frag)
    return frag
