# -*- coding: utf-8 -*-
"""Shared paths and imports for item B_overadjust_sim (POST HOC revision analysis, 2026-09-30).

Nothing in this item is registered or pre-specified: it was written in response to a
pre-submission review criticism (does the VIS10 matched-background control remove genuine
chemistry?).

The re-test instrument is IMPORTED, not rewritten, from the read-only internal repository:
  repo/scripts/run_phase2_claims_under_detectability_control.py  (log_odds_ratio, weighted_table,
      make_weighted_table, make_mh_log_odds_ratio, propensity_score, quantile_strata,
      nearest_neighbour_match, random_control_effect, bootstrap_interval, classify, DOWNGRADE,
      CALIPER_SD, HALDANE, RANDOM_DRAWS, SEED, SURVIVE_RETENTION)
  repo/scripts/phase2_claim_cohorts.py      (site_feature_matrix = VIS10, flank_residue_flag,
      _flank_residue_flag, SFE006_KR_OFFSETS, SNO021_K_OFFSETS)
  repo/scripts/phase2b_claim_cohorts.py     (flank_count_flag, SNO016_RESIDUES, PERS009_RESIDUES,
      read_fasta_gz)
  repo/scripts/phase2d_claim_cohorts.py     (SNO006_NEG_PLUS3, NEGATIVE_RESIDUES, SFE002_GLU_OFFSETS)
  repo/scripts/audit_cleaning_and_grouping.py (cleaning_collinearity: the matching pre-check)
  cys-audit/src/cys_audit (released tool v0.2.2: verdict.claim_retest with --claim-covariates)

sys.dont_write_bytecode is set before any import so that no __pycache__ is written into the
read-only repository or tool trees.
"""
from __future__ import annotations

import hashlib
import os
import sys

sys.dont_write_bytecode = True

W = r"C:/Users/admin/Desktop/小论文/_cys_repo_work/public/revision_2026-09-30"
REPO = r"C:/Users/admin/Desktop/小论文/_cys_repo_work/repo"
REPO_SCRIPTS = os.path.join(REPO, "scripts")
TOOL_SRC = r"C:/Users/admin/Desktop/小论文/_cys_repo_work/public/cys-audit/src"
ITEM = "B_overadjust_sim"
SCRIPTS = os.path.join(W, "scripts", ITEM)
RESULTS = os.path.join(W, "results", ITEM)
REPORT = os.path.join(W, "reports", ITEM + ".md")

CANDIDATES = os.path.join(W, "inputs", "repo_results", "phase3_candidate_peptides_2026-09-22.csv")
PHASE3_SUMMARY = os.path.join(W, "inputs", "repo_results", "phase3_empirical_detectability_summary_2026-09-22.json")
FASTA = os.path.join(W, "external", "UP000000589_10090.fasta.gz")
MCP_PKG = os.path.join(W, "inputs", "mcp_package")
SD3 = os.path.join(MCP_PKG, "Supplemental_Data_3_retest_round_a.csv")
SD4 = os.path.join(MCP_PKG, "Supplemental_Data_4_retest_round_b.csv")
SD5 = os.path.join(MCP_PKG, "Supplemental_Data_5_retest_round_d.csv")
SD6 = os.path.join(MCP_PKG, "Supplemental_Data_6_retest_round_e.csv")
SFE006_JSON = os.path.join(MCP_PKG, "Source_Data_text_phase2c_sfe006_reproduction.json")

UNIVERSE = os.path.join(RESULTS, "universe_cysteines.csv.gz")
PREP_SUMMARY = os.path.join(RESULTS, "prep_summary.json")

for p in (REPO_SCRIPTS, TOOL_SRC):
    if p not in sys.path:
        sys.path.insert(0, p)

MASTER_SEED = 20260930


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()
