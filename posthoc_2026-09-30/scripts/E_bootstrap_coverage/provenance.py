# -*- coding: utf-8 -*-
"""E_bootstrap_coverage -- POST HOC revision analysis (2026-09-30). Not registered.

Writes results/E_bootstrap_coverage/provenance.json: sha256 of every input read, of every script
of this item and of every output file; library versions; seeds; and a check of the registered
benchmark's strength-zero replicates (B1_cleavage, B2_abundance) with the same calibration statistic.
"""
from __future__ import annotations

import glob
import json
import os
import platform
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
import scipy
import statsmodels

import simlib as S
from summarize_grid import quantile_ci

W = r"/path/to/local/_cys_repo_work/public/revision_2026-09-30"
OUT = os.path.join(W, "results", "E_bootstrap_coverage")
SCRIPTS = os.path.join(W, "scripts", "E_bootstrap_coverage")
REPO = r"/path/to/local/_cys_repo_work/repo"

INPUTS = [
    os.path.join(W, "inputs", "phase4_inputs", "pxd063463_Trypsin_HydP.tsv"),
    os.path.join(W, "external", "UP000000589_10090.fasta.gz"),
    os.path.join(W, "inputs", "mcp_package", "Source_Data_Fig6_claim_retests.csv"),
    os.path.join(W, "inputs", "mcp_package", "Supplemental_Data_2_claim_verdict_tally.csv"),
    os.path.join(W, "inputs", "mcp_package", "Supplemental_Data_4_retest_round_b.csv"),
    os.path.join(W, "inputs", "mcp_package", "Supplemental_Data_5_retest_round_d.csv"),
    os.path.join(W, "inputs", "mcp_package", "Supplemental_Data_6_retest_round_e.csv"),
    os.path.join(W, "inputs", "repo_results", "phase2_synthetic_benchmark_replicates_2026-09-22.csv"),
    os.path.join(REPO, "results", "phase2_claim_retest_audit.json"),
    os.path.join(REPO, "results", "phase2b_claim_retest_audit.json"),
    os.path.join(REPO, "results", "phase2d_claim_retest_audit.json"),
    os.path.join(REPO, "results", "phase2e_claim_retest_audit.json"),
    os.path.join(REPO, "scripts", "run_phase2_claims_under_detectability_control.py"),
    os.path.join(REPO, "scripts", "run_phase2b_claims_backfill.py"),
    # added in the revision after verification: drivers whose verdict flow claims_table.py re-types,
    # and the registered benchmark generator / protocol / extension used by r1_b1_null.py
    os.path.join(REPO, "scripts", "run_phase2d_claims.py"),
    os.path.join(REPO, "scripts", "run_phase2e_claims.py"),
    os.path.join(REPO, "scripts", "run_phase2_synthetic_benchmark_2026-09-22.py"),
    os.path.join(REPO, "protocols", "phase2_synthetic_benchmark_preregistration_2026-09-22.json"),
    os.path.join(REPO, "results", "phase2_calibration_precision_2026-09-22.csv"),
    os.path.join(S.TOOL_SRC, "cys_audit", "stats.py"),
    os.path.join(S.TOOL_SRC, "cys_audit", "constants.py"),
    os.path.join(S.TOOL_SRC, "cys_audit", "__init__.py"),
    os.path.join(S.TOOL_SRC, "cys_audit", "verdict.py"),
    os.path.join(S.TOOL_SRC, "cys_audit", "io.py"),
    r"/path/to/local/巯基化/MCP/01_manuscript.tex",
    r"/path/to/local/巯基化/MCP/supplemental/Supplemental_Note_3_planted_artefact_benchmark.md",
]


def benchmark_check():
    df = pd.read_csv(INPUTS[7])
    out = {}
    for b in ("B1_cleavage", "B2_abundance"):
        g = df[(df.benchmark == b) & (df.strength == 0.0)]
        null = 0.0 if b == "B1_cleavage" else 0.5
        r = S.ratio_to_target(g.statistic.to_numpy(), g.ci_lo.to_numpy(), g.ci_hi.to_numpy(), null)
        w0, lo, hi = quantile_ci(r)
        out[b] = {"replicates": int(len(g)), "flagged": int(g.detected.sum()),
                  "w0_95th_percentile": w0, "w0_ci": [lo, hi],
                  "note": "50 registered strength-zero replicates only; the 300/200-replicate extension is stored as counts"}
    return out


def main():
    prov = {"label": "POST HOC revision analysis E_bootstrap_coverage (2026-09-30); not registered, not pre-specified",
            "python": sys.version, "platform": platform.platform(),
            "libraries": {"numpy": np.__version__, "scipy": scipy.__version__, "pandas": pd.__version__,
                          "statsmodels": statsmodels.__version__},
            "seeds": {"grid_master (run_grid.py, SeedSequence([20260930, cell_id, rep]))": 20260930,
                      "grid bootstrap seeds": "1e9 + cell_id*1e7 + rep*100 + {1 full, 2 MH, 3 sub, 50..69 random control}",
                      "claim-matched master (claims_matched_sim.py)": 20260931,
                      "claim-matched bootstrap seeds": "2e9 + config*1e7 + rep*100 + {1,2,3,50..69}",
                      "NEG_A calibration (calibrate_selected)": 12345,
                      "validation (validate.py)": "20260930 + {0..4}",
                      "revision after verification (r1_*.py) master": 20260932,
                      "r1_rc_yardsticks.py dataset-bootstrap rng": "default_rng(20260932), 2,000 resamples",
                      "r1_sfe008_sim.py": "SeedSequence([20260932, variant, rep]); bootstrap 3e9 + variant*1e7 + rep*100 + {3 sub, 50..69 rc}",
                      "r1_b1_null.py": "registered replicates seed_for(20260922, 0, 0, rep) (+1 bootstrap); fresh SeedSequence([20260932, 7|8, rep])",
                      "r1_confounding_anchored.py": "SeedSequence([20260932, 100 + cell, rep]); bootstrap 4e9 + cell*1e7 + rep*100 + 1",
                      "r2_anchor_refit.py": "deterministic (no random numbers)",
                      "r2_confounding_anchored.py": "SeedSequence([20260933, 200 + cell, rep]); bootstrap 5e9 + cell*1e7 + rep*100 + 1"},
            "bootstrap_replicates": 5000,
            "revision_after_verification_round1": {
                "date": "2026-09-30", "status": "POST HOC; not registered",
                "new_scripts": ["r1_rc_yardsticks.py", "r1_sfe008_sim.py", "r1_b1_null.py", "r1_confounding_anchored.py",
                                "r1_verdicts.py"],
                "changed_scripts": {"claims_table.py": "docstring only (verbatim vs re-typed parts); outputs re-generated and identical",
                                    "provenance.py": "hashes run_phase2d_claims.py, run_phase2e_claims.py and the benchmark files; separates verifier scripts"},
                "run_order": ["anchor_parameters.py", "validate.py", "run_grid.py", "summarize_grid.py", "claims_table.py",
                              "claims_matched_sim.py", "calibrated_verdicts.py", "multiplicity.py", "rule_difference.py",
                              "r1_rc_yardsticks.py", "r1_sfe008_sim.py", "r1_b1_null.py", "r1_confounding_anchored.py",
                              "r1_verdicts.py", "tables.py", "provenance.py"]},
            "revision_round2": {
                "date": "2026-09-30", "status": "POST HOC; not registered",
                "problem": ("verifier round 2, problem 1 (major): the anchoring random-intercept model was fitted with a fixed "
                            "20 x 20 non-adaptive Gauss-Hermite rule, inaccurate for clusters of up to 332 cysteines"),
                "new_scripts": ["r2_anchor_refit.py", "r2_confounding_anchored.py"],
                "changed_scripts": {"tables.py": "adds revision2_tables(); labels R1-g and R1-h as superseded",
                                    "provenance.py": "adds this block, the round-2 seeds, and hashes verify_r2_*.py separately"},
                "superseded_outputs": {
                    "anchor_parameters.json": "the 'glmm' block of each attribute (rho, SEs); descriptive cohort values stand",
                    "r1_confounding_anchored_cells.csv": "replaced by r2_confounding_anchored_cells.csv",
                    "r1_confounding_anchored_targets.csv": "replaced by r2_confounding_anchored_targets.csv",
                    "r1_confounding_anchored.csv.gz": "replaced by r2_confounding_anchored.csv.gz"},
                "run_order": ["anchor_parameters.py", "validate.py", "run_grid.py", "summarize_grid.py", "claims_table.py",
                              "claims_matched_sim.py", "calibrated_verdicts.py", "multiplicity.py", "rule_difference.py",
                              "r1_rc_yardsticks.py", "r1_sfe008_sim.py", "r1_b1_null.py", "r1_verdicts.py",
                              "r2_anchor_refit.py", "r2_confounding_anchored.py", "tables.py", "provenance.py"],
                "not_rerun": "r1_confounding_anchored.py (superseded; its outputs are kept for the record)"},
            "inputs_sha256": {p: S.sha256(p) for p in INPUTS if os.path.exists(p)},
            "scripts_sha256": {os.path.basename(p): S.sha256(p) for p in sorted(glob.glob(os.path.join(SCRIPTS, "*.py")))
                               if not os.path.basename(p).startswith(("verify_r1_", "verify_r2_"))},
            "verifier_scripts_sha256_not_run_by_implementer": {
                os.path.basename(p): S.sha256(p) for p in sorted(glob.glob(os.path.join(SCRIPTS, "verify_r*_*.py")))},
            "verifier_outputs_compared_in_report_sha256": {
                os.path.relpath(p, OUT): S.sha256(p) for p in sorted(glob.glob(os.path.join(OUT, "verify_r2", "anchor_*"))
                                                                     + glob.glob(os.path.join(OUT, "verify_r2", "confounding_*")))},
            "outputs_sha256": {os.path.relpath(p, OUT): S.sha256(p) for p in sorted(
                glob.glob(os.path.join(OUT, "*")) + glob.glob(os.path.join(OUT, "sim_cells", "*")))
                if os.path.isfile(p) and not p.endswith("provenance.json")},
            "registered_benchmark_strength_zero_check": benchmark_check(),
            "report_sha256": {"reports/E_bootstrap_coverage.md": S.sha256(os.path.join(W, "reports", "E_bootstrap_coverage.md"))}}
    with open(os.path.join(OUT, "provenance.json"), "w", encoding="utf-8") as fh:
        json.dump(prov, fh, indent=1)
    print(json.dumps(prov["registered_benchmark_strength_zero_check"], indent=1))
    print("inputs", len(prov["inputs_sha256"]), "scripts", len(prov["scripts_sha256"]), "outputs", len(prov["outputs_sha256"]))


if __name__ == "__main__":
    main()
