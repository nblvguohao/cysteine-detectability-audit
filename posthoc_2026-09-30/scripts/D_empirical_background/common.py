"""Shared paths and helpers for item D_empirical_background (POST HOC revision analysis, 2026-09-30).

Everything in this item is a post hoc analysis written in response to review criticism; nothing here was
registered or pre-specified.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import sys

sys.dont_write_bytecode = True  # never write __pycache__ into the read-only Cys-Audit source tree

W = "/path/to/local/_cys_repo_work/public/revision_2026-09-30"
ITEM = "D_empirical_background"
SCRIPTS = f"{W}/scripts/{ITEM}"
RESULTS = f"{W}/results/{ITEM}"
REPORT = f"{W}/reports/{ITEM}.md"
CYS_AUDIT_SRC = "/path/to/local/_cys_repo_work/public/cys-audit/src"
REPO = "/path/to/local/_cys_repo_work/repo"

IN_HYDP = f"{W}/inputs/phase4_inputs/pxd063463_Trypsin_HydP.tsv"
IN_HYDN = f"{W}/inputs/phase4_inputs/pxd063463_Trypsin_HydN.tsv"
IN_CAND = f"{W}/inputs/repo_results/phase3_candidate_peptides_2026-09-22.csv"
IN_EMP_SUMMARY = f"{W}/inputs/repo_results/phase3_empirical_detectability_summary_2026-09-22.json"
IN_FASTA_GZ = f"{W}/external/UP000000589_10090.fasta.gz"
STORED_AUDIT_PROTEOME = f"{REPO}/results/phase4_reports/pxd063463_Trypsin_HydP_trypsin_proteome/audit.json"
STORED_AUDIT_OBSERVED = f"{REPO}/results/phase4_reports/pxd063463_Trypsin_HydP_trypsin_observed/audit.json"
ORIG_EMP_SCRIPT = f"{REPO}/scripts/run_phase3_empirical_detectability_2026-09-22.py"
ORIG_EMP_PROTOCOL = f"{REPO}/protocols/phase3_empirical_detectability_preregistration_2026-09-22.json"

MASTER_SEED = 20260930          # seed family for everything new in this item
TOOL_SEED = 20260922            # Cys-Audit default seed; used for the reproduction runs
BOOT_REPS = 5000
N_JOBS = 4

if CYS_AUDIT_SRC not in sys.path:
    sys.path.insert(0, CYS_AUDIT_SRC)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_lf(path):
    """sha256 after CRLF -> LF normalisation (the stored audits hashed LF files; the Windows copies are CRLF)."""
    with open(path, "rb") as fh:
        b = fh.read()
    return hashlib.sha256(b.replace(b"\r\n", b"\n")).hexdigest()


def read_fasta_gz(path):
    """{UniProt accession: sequence} from a gzipped UniProt FASTA (second |-field of the header)."""
    seqs, acc, buf = {}, None, []
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if acc is not None:
                    seqs[acc] = "".join(buf).upper().rstrip("*")
                tok = line[1:].split()[0]
                parts = tok.split("|")
                acc = parts[1] if len(parts) >= 3 else tok
                buf = []
            elif line:
                buf.append(line.strip())
    if acc is not None:
        seqs[acc] = "".join(buf).upper().rstrip("*")
    return seqs


def dump_json(obj, path):
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=1, ensure_ascii=False, default=float)
        fh.write("\n")
