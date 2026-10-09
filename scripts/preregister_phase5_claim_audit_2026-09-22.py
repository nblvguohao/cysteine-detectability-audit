"""Phase 5 pre-registration: the claim audit matrix and an independent-control cross-check with Cys-Audit.

Written before scripts/run_phase5_claim_audit_2026-09-22.py exists or has run. The stored verdict tables are
already known (Phase 0 re-ran them byte-identically), so M1 below is a recount, not a blind test. X1 (the Cys-Audit
cross-check) has not been run on any claim; no tool status for any claim is known at registration.

M1 recount (a failure stops Phase 5). Recompute one verdict per claim from the four stored tables
   (results/phase2_claim_retest.csv, phase2b_claim_retest_combined.csv, phase2d_claim_retest.csv,
   phase2e_claim_retest.csv) with the de-duplication rule written in scripts/recount_claim_verdicts_2026-09-17.py
   (primary caliber only; phase2b combined supersedes phase2 for the eight A-layer claims; latest round wins;
   out_of_instrument_scope and blocked_material_unreachable = no measured baseline), re-implemented here without
   importing that script. Must equal results/claim_verdict_tally_2026-09-17.csv claim by claim, and give
   35 claims, 28 with a measured baseline: survives 12, undecidable 11, null_broken_by_control 2, attenuated 1,
   vanishes 1, baseline_contradicts_claim 1.
M2 claim audit matrix results/claim_audit_2026-09-22.csv, one row per claim (35), columns as in the brief s.14
   plus provenance. Definitions:
     effect_original  baseline log2 OR [CI] of the counted round (primary caliber)
     effect_control   matched log2 OR [CI] (propensity matching on the round's covariate set)
     p_original       baseline Fisher p as stored; p_control left empty with the note "interval-based; no p"
                      (no p was computed for the matched estimate, and none is invented here)
     direction_change no_baseline_effect (baseline CI covers 0) | lost_under_control (matched CI covers 0) |
                      same_sign | reversed
     verdict_status   survives->PASS, attenuated->ATTENUATED, undecidable->UNDECIDABLE,
                      vanishes / baseline_contradicts_claim / null_broken_by_control->FAIL
                      (null_broken_by_control is the tree's name, for a claim of NO preference, of what the
                      shared rule calls "reverses": the control establishes the preference the claim denied)
     evidence_level   brief s.8 hierarchy: L2_reanalysis_author_data (retest on the authors' own deposited
                      data), L2_transfer_other_cohort (is_transfer), none (no measured baseline)
     species_tested   from the reference-proteome file(s) the cohort builder opens (key -> species below);
                      "undetermined" if it opens none
M3 concentration: claims and verdicts per publication, recomputed from the matrix (s.9.18 reported 15 papers,
   max 7 from one; recomputed, not copied).
X1 Cys-Audit cross-check (descriptive; never changes a stored verdict).
   For each of the 28 claims with a baseline, rebuild the cohort with the round's own builder (canonical
   interpreter 3.11.16). Gate G1: rebuilt n_observations and n_positive equal the stored row; a claim failing G1
   is excluded and reported. Site-level cohorts only (protein-level: "not_applicable"): export protein, position,
   label = y, claim_attr = attribute, cluster = the builder's groups (homology components), FASTA from the
   cohort's sequences; audit with the frozen Phase 4 tool: --protease trypsin --background proteome
   --claim-feature claim_attr, direction positive (positive_preference) or null (null_no_preference), default
   seed, 5000 replicates. The tool's control stratifies on trypsin cleavage-band flags only (these cohorts carry
   no abundance or peptide columns), a WEAKER control than the tree's propensity matching on the round's
   covariate set. Pre-declared reading:
     agree      same four-level status
     lenient    tool status less severe than the tree's (order FAIL > ATTENUATED > UNDECIDABLE > PASS is NOT
                used; instead: tool PASS where tree is not PASS) - expected from the weaker control, not
                evidence against the tree
     stricter   tool FAIL or ATTENUATED where tree is PASS - flagged for review
     other      any remaining disagreement
   Also reported per claim: the tool's baseline estimate vs the stored baseline (same cohort, Haldane log2 OR in
   both; any difference is reported, not reconciled), and the tool's dataset-level cleavage status.
Stop rule: M1 mismatch, or G1 failing on more than 3 of 28 claims, stops Phase 5.
Standard library only.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "protocols", "phase5_claim_audit_preregistration_2026-09-22.json")
AUD = os.path.join(ROOT, "results", "phase5_claim_audit_prereg_2026-09-22_audit.json")
RUNNER = os.path.join(ROOT, "scripts", "run_phase5_claim_audit_2026-09-22.py")
FUTURE = [os.path.join(ROOT, "results", f) for f in ("claim_audit_2026-09-22.csv", "phase5_claim_audit_2026-09-22_audit.json",
                                                     "phase5_cys_audit_crosscheck_2026-09-22.csv")]
INPUTS = ["results/phase2_claim_retest.csv", "results/phase2b_claim_retest_combined.csv",
          "results/phase2d_claim_retest.csv", "results/phase2e_claim_retest.csv",
          "results/claim_verdict_tally_2026-09-17.csv", "results/ptm_site_preference_claims.csv",
          "scripts/phase2_claim_cohorts.py", "scripts/phase2b_claim_cohorts.py", "scripts/phase2d_claim_cohorts.py",
          "scripts/phase2e_claim_cohorts.py", "protocols/phase4_cys_audit_preregistration_2026-09-22.json"]
SPECIES = {"ath": "Arabidopsis thaliana", "hsa": "Homo sapiens", "osa": "Oryza sativa", "mmu": "Mus musculus",
           "pvu": "Phaseolus vulgaris", "tgo": "Toxoplasma gondii"}


def sha(p):
    with open(p, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def main():
    for p in (OUT, AUD):
        if os.path.exists(p):
            sys.exit(f"REFUSE: {p} exists")
    if os.path.exists(RUNNER) or any(os.path.exists(p) for p in FUTURE):
        sys.exit("REFUSE: the Phase 5 runner or its outputs already exist")
    protocol = {
        "phase": 5, "date": "2026-09-22",
        "inputs_sha256": {p: sha(os.path.join(ROOT, p)) for p in INPUTS},
        "M1_expected": {"claims": 35, "with_baseline": 28,
                        "counts": {"survives": 12, "undecidable": 11, "null_broken_by_control": 2, "attenuated": 1,
                                   "vanishes": 1, "baseline_contradicts_claim": 1},
                        "no_baseline": ["out_of_instrument_scope", "blocked_material_unreachable"]},
        "verdict_status": {"survives": "PASS", "attenuated": "ATTENUATED", "undecidable": "UNDECIDABLE",
                           "vanishes": "FAIL", "baseline_contradicts_claim": "FAIL", "null_broken_by_control": "FAIL"},
        "species_keys": SPECIES,
        "X1": {"tool_protocol": "protocols/phase4_cys_audit_preregistration_2026-09-22.json", "protease": "trypsin",
               "background": "proteome", "bootstrap_reps": 5000, "interpreter": "3.11.16",
               "direction_map": {"positive_preference": "positive", "null_no_preference": "null"},
               "gate_G1": "rebuilt n_observations and n_positive equal the stored row",
               "stop_if_G1_fail_more_than": 3},
    }
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(protocol, fh, indent=1, ensure_ascii=False)
    audit = {"script": "scripts/preregister_phase5_claim_audit_2026-09-22.py", "script_sha256": sha(os.path.abspath(__file__)),
             "protocol": os.path.relpath(OUT, ROOT), "protocol_sha256": sha(OUT), "runner_absent_when_written": True}
    with open(AUD, "w", encoding="utf-8") as fh:
        json.dump(audit, fh, indent=1)
    print(json.dumps(audit, indent=1))


if __name__ == "__main__":
    main()
