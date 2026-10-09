"""Shared paths and helpers for item I_fig5a_estimand.

POST HOC revision analysis (2026-09-30) in response to pre-submission review criticism of
Figure 5a ("detectability share"). Nothing in this folder is registered or pre-specified.
All inputs are opened read-only; outputs go only to results/I_fig5a_estimand/.
"""
from __future__ import annotations

import hashlib
import json
import pathlib

ITEM = "I_fig5a_estimand"
W = pathlib.Path(r"/path/to/local/_cys_repo_work/public/revision_2026-09-30")
REPO = pathlib.Path(r"/path/to/local/_cys_repo_work/repo")
MCP = pathlib.Path(r"/path/to/local/巯基化/MCP")
OUT = W / "results" / ITEM
OUT.mkdir(parents=True, exist_ok=True)

# primary inputs (read-only)
SHARE_CSV = W / "inputs/repo_results/ptm_detectability_share.csv"
SOURCE_DATA_FIG5 = W / "inputs/mcp_package/Source_Data_Fig5_three_axes.csv"
SD7_CENSUS = MCP / "supplemental/Supplemental_Data_7_dataset_census.csv"
SNOTE7 = MCP / "supplemental/Supplemental_Note_7_artefact4_public_deposit.md"
MANUSCRIPT = MCP / "01_manuscript.tex"
CENSUS_AUDIT = REPO / "results/ptm_detectability_census_audit.json"
CENSUS_SITES_AUDIT = REPO / "results/ptm_census_sites_audit.json"
CENSUS_VERIFICATION = REPO / "results/ptm_detectability_census_verification.csv"
DIAGNOSTICS_AUDIT = REPO / "results/ptm_detectability_diagnostics_audit.json"
CENSUS_SCRIPT = REPO / "scripts/run_ptm_detectability_census.py"
PROBE_SCRIPT = REPO / "scripts/run_cross_protease_detectability_probe.py"
INGEST_SCRIPT = REPO / "scripts/ingest_ptm_census_sites.py"
DIAGNOSTICS_SCRIPT = REPO / "scripts/ptm_detectability_diagnostics.py"
MOUSE_FASTA = W / "external/UP000000589_10090.fasta.gz"

SEED = 20260930

# Figure 5a order (top to bottom) and labels exactly as printed in the figure / Supplemental Note 7
COHORTS = [
    ("qtrp_S1_ph5", "QTRP S1 pH5 (human)"),
    ("qtrp_S2_ph5", "QTRP S2 pH5 (human)"),
    ("qpers_sid_tierB", "qPerS-SID tier B (human)"),
    ("cysboost2019_human_sno_hela", "Cys-BOOST HeLa (human)"),
    ("cysboost2019_human_sno_shsy5y", "Cys-BOOST SH-SY5Y (human)"),
    ("natcomm2023_ath_sno", "FAT-switch (Arabidopsis)"),
    ("abiotech2025_ath_sno", "PAT-switch (Arabidopsis)"),
    ("fps2020_ath_sulfenyl", "YAP1C reporter (Arabidopsis)"),
]
LABEL = dict(COHORTS)
NEG = {"NEG_A_not_observed": "NEG_A", "NEG_B_observed_unmodified": "NEG_B"}
FEATURE_SETS = ("VIS10", "DIG25", "CPL15")
SHARE_DENOM_FLOOR = 0.02  # declared in the original census script


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_lf(path: pathlib.Path) -> str:
    """sha256 after CRLF->LF normalisation (the Windows checkout converts line endings)."""
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def record_inputs(step: str, paths: list, extra: dict | None = None) -> None:
    """Write the list of inputs a step read, with hashes, for the provenance merge."""
    rec = {"step": step, "inputs": {}}
    for p in paths:
        p = pathlib.Path(p)
        rec["inputs"][str(p)] = {"exists": p.exists(),
                                 "sha256": sha256(p) if p.exists() and p.is_file() else None,
                                 "bytes": p.stat().st_size if p.exists() and p.is_file() else None}
    if extra:
        rec.update(extra)
    (OUT / f"_inputs_{step}.json").write_text(json.dumps(rec, indent=1, ensure_ascii=False), encoding="utf-8")
