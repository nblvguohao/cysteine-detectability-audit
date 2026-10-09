"""Shared paths and helpers for item C_pxd063463_specific.

POST HOC revision analysis (2026-09-30) in response to pre-submission review criticism. Nothing in this item was
registered or pre-specified; every definition below was written after the stored Figure 2 results had been seen.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import platform
import sys

import numpy as np
import pandas as pd

W = "/path/to/local/_cys_repo_work/public/revision_2026-09-30"
ITEM = "C_pxd063463_specific"
SCRIPTS = f"{W}/scripts/{ITEM}"
RES = f"{W}/results/{ITEM}"
WORK = f"{RES}/_work"
IN_DIR = f"{W}/inputs/phase4_inputs"
FASTA_GZ = f"{W}/external/UP000000589_10090.fasta.gz"
FASTA = f"{WORK}/UP000000589_10090_2026_03.fasta"          # decompressed copy (derived, byte-identical content)
TOOL_SRC = "/path/to/local/_cys_repo_work/public/cys-audit/src"
TOOL_V010 = "/path/to/local/_cys_repo_work/public/cys-audit/archive/v0.1.0/src"
REPO = "/path/to/local/_cys_repo_work/repo"
STORED_REPORTS = f"{REPO}/results/phase4_reports"
STORED_DEMO = f"{W}/inputs/repo_results/phase4_public_demo_2026-09-22.csv"
SOURCE_FIG2 = f"{W}/inputs/mcp_package/Source_Data_Fig2_public_four_protease.csv"
SUPP_NOTE4 = "/path/to/local/巯基化/MCP/supplemental/Supplemental_Note_4_four_protease_public_deposit.md"
MANUSCRIPT = "/path/to/local/巯基化/MCP/01_manuscript.tex"

ARMS = ["Trypsin", "AspN", "CT", "GluC"]
OWN_RULE = {"Trypsin": "trypsin", "AspN": "aspn", "CT": "chymotrypsin", "GluC": "gluc"}   # as chosen in Phase 4
ARM_LABEL = {"Trypsin": "trypsin", "AspN": "AspN", "CT": "chymotrypsin", "GluC": "GluC"}
TOOL_SEED = 20260922          # Cys-Audit default seed, used by the stored Phase 4 audits
NEW_SEED = 20260930           # seed for the new (non-tool) bootstrap analyses of this item
REPS = 5000

if TOOL_SRC not in sys.path:
    sys.path.insert(0, TOOL_SRC)

PROV_FILE = f"{RES}/_provenance_parts.jsonl"


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def record_inputs(step, paths, extra=None):
    """Append one provenance line (inputs read by a step, with sha256)."""
    os.makedirs(RES, exist_ok=True)
    rec = {"step": step, "inputs": {p: sha256_file(p) for p in paths}}
    if extra:
        rec.update(extra)
    with open(PROV_FILE, "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec, ensure_ascii=False) + "\n")


def ensure_fasta():
    """Decompress the UniProt 2026_03 mouse FASTA once (Cys-Audit reads plain-text FASTA)."""
    os.makedirs(WORK, exist_ok=True)
    if not os.path.exists(FASTA):
        tmp = FASTA + ".part"
        with gzip.open(FASTA_GZ, "rb") as src, open(tmp, "wb") as dst:
            for block in iter(lambda: src.read(1 << 20), b""):
                dst.write(block)
        os.replace(tmp, FASTA)
    return FASTA


def load_fasta():
    from cys_audit.io import read_fasta
    return read_fasta(ensure_fasta())


def raw_path(arm, hyd):
    return f"{IN_DIR}/pxd063463_{arm}_{hyd}.tsv"


def filtered_path(arm, hyd):
    return f"{RES}/inputs_filtered/pxd063463_{arm}_{hyd}.tsv"


# Base table set for the count / coincidence / depth steps (c02, c04, c05):
#   filtered   (default) rows whose protein is absent from UniProt 2026_03 dropped (the base of every cleavage statistic)
#   unfiltered the deposit-converted tables as stored (the base of the manuscript's current counts, e.g. 1,359 / 5,328)
BASE = os.environ.get("C_BASE", "filtered")
assert BASE in ("filtered", "unfiltered"), BASE
OUT = RES if BASE == "filtered" else f"{RES}/unfiltered_base"


def arm_path(arm, hyd):
    return filtered_path(arm, hyd) if BASE == "filtered" else raw_path(arm, hyd)


def read_arm(arm, hyd):
    return pd.read_csv(arm_path(arm, hyd), sep="\t", dtype={"protein": str})


def env_record():
    import scipy
    import sklearn
    import statsmodels
    return {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
            "scipy": scipy.__version__, "statsmodels": statsmodels.__version__, "scikit-learn": sklearn.__version__,
            "platform": platform.platform()}


def window(seq, pos, half=10):
    """Sequence window of +/-half residues around 1-based pos, padded with '-'."""
    if seq is None or pos > len(seq):
        return None
    i0 = pos - 1
    left = seq[max(0, i0 - half):i0].rjust(half, "-")
    right = seq[i0 + 1:i0 + 1 + half].ljust(half, "-")
    return left + seq[i0] + right
