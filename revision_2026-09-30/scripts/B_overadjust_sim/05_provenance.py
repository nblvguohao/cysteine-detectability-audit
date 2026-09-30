# -*- coding: utf-8 -*-
"""B_overadjust_sim step 5 (POST HOC): provenance.json with sha256 of every input read, of every
imported repository / tool module, of this item's scripts and outputs, library versions and seeds."""
from __future__ import annotations

import glob
import json
import os
import platform
import sys
import time

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common_b as cb  # noqa: E402


def main():
    import matplotlib
    import numpy
    import pandas
    import scipy
    import sklearn
    import sim_engine as se
    import run_phase2_claims_under_detectability_control as rt
    from cys_audit import constants as tc

    inputs = {
        "phase3_candidate_peptides_2026-09-22.csv": cb.CANDIDATES,
        "UP000000589_10090.fasta.gz": cb.FASTA,
        "Supplemental_Data_4_retest_round_b.csv": cb.SD4,
        "Supplemental_Data_5_retest_round_d.csv": cb.SD5,
        "Supplemental_Data_6_retest_round_e.csv": cb.SD6,
        "Source_Data_text_phase2c_sfe006_reproduction.json": cb.SFE006_JSON,
    }
    repo_modules = ["run_phase2_claims_under_detectability_control.py", "phase2_claim_cohorts.py",
                    "phase2b_claim_cohorts.py", "phase2d_claim_cohorts.py",
                    "run_cross_protease_detectability_probe.py", "audit_cleaning_and_grouping.py", "common.py",
                    "run_phase2b_claims_backfill.py"]
    tool_modules = ["cys_audit/verdict.py", "cys_audit/propensity.py", "cys_audit/stats.py", "cys_audit/io.py",
                    "cys_audit/constants.py"]
    prov = {
        "item": cb.ITEM,
        "label": "POST HOC revision analysis in response to pre-submission review criticism; not registered, "
                 "not pre-specified",
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "inputs_sha256": {k: cb.sha256_of(v) for k, v in inputs.items()},
        "inputs_paths": inputs,
        "inputs_read_but_not_hashed_elsewhere": {
            "note": "Supplemental_Data_3 (round a) was inspected interactively; its SNO-021 and SFE-006 rows are "
                    "identical to the phase-2 rows carried in Supplemental_Data_4, which is the file the scripts read",
            "Supplemental_Data_3_retest_round_a.csv": cb.sha256_of(cb.SD3)},
        "repository_modules_imported_sha256": {m: cb.sha256_of(os.path.join(cb.REPO_SCRIPTS, m)) for m in repo_modules},
        "tool_modules_imported_sha256": {m: cb.sha256_of(os.path.join(cb.TOOL_SRC, m)) for m in tool_modules},
        "tool_version": tc.VERSION,
        "item_scripts_sha256": {os.path.basename(p): cb.sha256_of(p)
                                for p in sorted(glob.glob(os.path.join(cb.SCRIPTS, "*.py")))},
        "outputs_sha256": {os.path.basename(p): cb.sha256_of(p)
                           for p in sorted(glob.glob(os.path.join(cb.RESULTS, "*")))
                           if os.path.isfile(p) and not p.endswith("provenance.json")},
        "seeds": {"master_seed": cb.MASTER_SEED,
                  "replicate_seed": "numpy SeedSequence([20260930, scenario_code, attribute_code, "
                                    "round(effect*100), round(rate*1000), rep]); round-3 mixture tasks append "
                                    "[attribute2_code, round(effect2*100)]; round-4 length-selection tasks append "
                                    "[15, round(gamma*1000)]",
                  "scenario_codes": se.SCEN_CODES,
                  "pipeline_bootstrap_and_random_control_seed": rt.SEED,
                  "random_control_draw_seeds": "%d..%d" % (rt.SEED, rt.SEED + rt.RANDOM_DRAWS - 1),
                  "tool_seed": tc.DEFAULT_SEED,
                  "cv_fold_seed_prep": cb.MASTER_SEED},
        "settings": {"n_boot_intervals": se.N_BOOT, "n_boot_check": 5000, "caliper_sd": rt.CALIPER_SD,
                     "haldane": rt.HALDANE, "survive_retention": rt.SURVIVE_RETENTION,
                     "random_draws": rt.RANDOM_DRAWS, "latent_prevalence": se.LATENT_PREVALENCE,
                     "positive_rate_among_detected_primary": 0.20, "steep_gamma": 2.0,
                     "revision_steep_vis_gamma": 0.25, "revision_steep_nocomp_gamma": 3.0,
                     "revision_steep_strong_gamma": 4.0,
                     "revision_gamma_rule": "01e_calibrate_gamma_alt.py: gamma whose null propensity AUC is "
                                            "closest to the original steep model's 0.722 (grid extended "
                                            "downwards for steep_vis after the first run; documented there)",
                     "worker_processes": 4, "threads_per_process": 1},
        "revision_after_verification": {
            "date": "2026-09-30",
            "new_scripts": ["01d_extra_universe.py", "01e_calibrate_gamma_alt.py"],
            "modified_scripts": ["sim_engine.py (backward compatible; 07_rerun_check.py)", "02_run_simulation.py",
                                 "03_summarise.py", "04_figures.py", "06_report_tables.py", "07_rerun_check.py",
                                 "05_provenance.py"],
            "new_groups": ["permgrid", "steepalt", "art1extra", "art1steep4"],
            "derived_inputs": {"universe_extra.csv.gz": "01d_extra_universe.py from the same candidate table, "
                                                        "FASTA and universe"},
            "note_on_read_only_trees": "sys.dont_write_bytecode is set by every item script. Two __pycache__ "
                                       "folders exist in read-only trees and were not created by this item: "
                                       "cys-audit/src/cys_audit/__pycache__ (01:19:41, before this item started) "
                                       "and repo/scripts/__pycache__ (03:39, reported by the verifier as its own)."},
        "revision_round2": {
            "date": "2026-09-30",
            "label": "POST HOC; implements the verifier's round-2 proxy-chemistry model (genuine chemistry that "
                     "depends on the local K/R count within +/-20 residues) and the count-defined claim",
            "new_scripts": ["01f_krcount_universe.py"],
            "modified_scripts": ["sim_engine.py (backward compatible; 07_rerun_check.py)", "02_run_simulation.py",
                                 "03_summarise.py", "04_figures.py", "06_report_tables.py", "07_rerun_check.py",
                                 "05_provenance.py"],
            "new_groups": ["countclaim", "proxygrid", "proxysens", "proxyci", "proxyart1"],
            "attribute_and_feature_codes_added": {"x1_KRcount20_hi": 8, "x2_KR5_ge3": 9, "F_KRcount20_z": 10,
                                                  "F_KR5_z": 11, "F_AKRV10_z": 12},
            "reused_seeds": "countclaim null tables re-use the S0/S0steep/S0perm seeds; proxyart1 (scenarios "
                            "P1gridA1/P1steepgridA1) re-uses the P1grid/P1steepgrid seeds via seed_scenario; "
                            "all verified identical (rerun_check_countclaim_tables.csv)",
            "stopped_run": "a first launch of proxygrid with the Artifact-1 estimators on all ten attributes "
                           "(~7 s per table) was stopped before any output was written (run_log_round2.txt, first "
                           "'[group] proxygrid' line); proxygrid was re-run without them and proxyart1 computes them "
                           "on identical tables",
            "derived_inputs": {"universe_krcount.csv.gz": "01f_krcount_universe.py: the repository's DIG25 column "
                                                          "pep_cleavage_sites_within_20 (site_feature_matrix) and "
                                                          "flank_count_flag({K,R}) on the same universe"},
            "proxy_model": "P(positive | detected) = expit(a + b ln2 F), F = standardised K/R count within +/-20 "
                           "(population SD over the universe); b = log2 OR per SD; a solved so that the mean over "
                           "detected cysteines equals the positive rate",
            "proxy_grid_b": [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0],
            "proxyci_rule": "grid b nearest to the b whose mean simulated SFE-006 baseline equals the observed 1.384 "
                            "(02_run_simulation.proxy_ci_effect), applied per detection model",
            "verifier_outputs_read_sha256": {
                "verify_r2/universe_altchem.csv.gz (cross-check of the K/R count only)":
                    cb.sha256_of(os.path.join(cb.RESULTS, "verify_r2", "universe_altchem.csv.gz"))}},
        "revision_round3": {
            "date": "2026-09-30",
            "label": "POST HOC; implements the round-3 verifier's sequence-only discriminator (the association of "
                     "the label with a local residue count relative to its association with the claimed SFE-006 "
                     "offsets, in the same site table) and adds the mixture and cohort-size checks",
            "new_scripts": ["01g_discrim_universe.py", "08_discriminator.py", "09_sfe006_count_test.py"],
            "modified_scripts": ["sim_engine.py (backward compatible; 07_rerun_check.py)", "02_run_simulation.py",
                                 "04_figures.py", "06_report_tables.py", "07_rerun_check.py", "05_provenance.py"],
            "new_groups": ["discrim", "proxysensext", "discrim2"],
            "attribute_codes_added": {"x3_AKRV10_hi": 13, "x4_KR5_hi": 14},
            "scenario_codes_added": {k: v for k, v in se.SCEN_CODES.items() if v >= 31},
            "statistic": "baseline log2 OR of a count-defined attribute divided by the baseline log2 OR of the "
                         "SFE-006 offsets, in the same table (NEG_A background, no matching)",
            "count_attributes": "count above the median of the simulation universe: K/R within +/-20 > 4 "
                                "(01f), A/K/R/V within +/-10 > 5 and K/R within +/-5 > 1 (01g)",
            "calibration_rule": "02_run_simulation.exact_sfe006_effect: the planted strength at which the mean "
                                "simulated SFE-006 baseline equals the observed 1.3840, by linear interpolation "
                                "along the stored grid of that scenario; exact value, not the nearest grid point. "
                                "For the mixtures, exact_mix_c gives the offsets strength c that restores that "
                                "baseline at each count strength b. Both rules were fixed before the runs they "
                                "select; exact_mix_c reads the DMIXgrid output of the first stage.",
            "mixture_model": "a second planted term on the logit among detected cysteines (task['attribute2'], "
                             "task['effect2']); only tasks carrying it add those two fields to their seed, so every "
                             "earlier table is unchanged",
            "reused_seeds": "the six detectability-only cells (D0real, D0steep, D0steep4, D0steepvis, D0steepnc, "
                            "D0perm) re-use the S0/S0steep/S0art1steep4/S0steepvis/S0steepnc/S0perm seeds via "
                            "seed_scenario; all 1,600 verified identical (rerun_check_discrim_tables.csv)",
            "reading_rule": "results/B_overadjust_sim/discriminator_reading_rule.json, written by 08_discriminator.py "
                            "from the simulation alone and hashed into sfe006_count_test.json before any value was "
                            "computed on the authors' cohort",
            "cohort_test_status": "not computed here: results/phase2c_yang2014_sulfenyl_sites.csv and "
                                  "external/proteomes/hsa.fasta.gz are absent from this repository copy and no "
                                  "network was used (sfe006_count_test.json). 09_sfe006_count_test.py --selftest "
                                  "reproduces the engine's baselines exactly on simulated tables "
                                  "(sfe006_count_test_selftest.csv).",
            "derived_inputs": {"universe_discrim.csv.gz": "01g_discrim_universe.py: the A/K/R/V +/-10 and K/R +/-5 "
                                                          "counts of universe_krcount.csv.gz, thresholded at their "
                                                          "universe medians"},
            "verifier_outputs_read_sha256": {
                p: cb.sha256_of(os.path.join(cb.RESULTS, "verify_r3", p))
                for p in ("discriminator_r3_summary.csv", "r3_G4x.csv.gz")
                if os.path.exists(os.path.join(cb.RESULTS, "verify_r3", p))}},
        "revision_round4": {
            "date": "2026-09-30",
            "label": "POST HOC; answers the round-4 verification, which refuted the round-3 additions. Re-implements "
                     "the verifier's detectability-only model (selection of positives on tryptic-peptide length) in "
                     "this item's engine, calibrates it to the observed SFE-006 baseline, and reports what it does to "
                     "every attribute, verdict and exploratory estimator; extends the K/R +/-5 chemistry grid until it "
                     "reaches that baseline; and states the overlaps of the count/offset ratio",
            "withdrawn_in_round4": [
                "the count/offset ratio as a discriminator: its simulated ranges overlap once a detectability-only "
                "model reaches the observed baseline (discriminator_overlap_groups.csv); it is kept only as an "
                "exploratory description in Supplemental Note 15 and never appears in the proposed main text",
                "the reading rule discriminator_reading_rule.json: withdrawn before any cohort value was computed; "
                "09_sfe006_count_test.py no longer applies it and is descriptive only (the file and its sha256 are "
                "kept unchanged as a record of round 3)",
                "the bounds on the artifact reading ('neither reading is cheap'; 'detectability alone gives that "
                "count 0.3-0.6'; 'no simulated detection process produced a detectability-only loss of the size seen "
                "here'; '5-10x'): the length-selection model reproduces the observed baseline (1.381 against 1.3840) "
                "and propensity AUC (0.765 against 0.762), removes more than observed (loss 1.24 against 0.79), and "
                "gives the K/R +/-20 attribute a baseline of 1.18 with no chemistry",
                "the claim that an observed-negative background would separate artifact from chemistry directly: "
                "against detected non-positives the length-selection artifact keeps 0.89 of its baseline and "
                "K/R-count chemistry 0.92 (negb_detmatch_comparison.csv)"],
            "new_scripts": ["10_round4_summary.py"],
            "modified_scripts": ["sim_engine.py (backward compatible; 07_rerun_check.py)", "02_run_simulation.py",
                                 "04_figures.py", "06_report_tables.py", "07_rerun_check.py",
                                 "09_sfe006_count_test.py (rule no longer applied)", "05_provenance.py"],
            "new_groups": ["lenselgrid", "lenselcal", "lenselci", "kr5ext", "kr5cal"],
            "detection_mode_added": {"lensel7": "detectability only, no chemistry: among detected cysteines "
                                                "P(positive) = expit(a + gamma * -ln max(pep_len, 7)), pep_len being "
                                                "the VIS10 column 'pep_len'; re-implementation of the round-4 "
                                                "verifier's verify_r4_02_engine.py mode of the same name"},
            "scenario_codes_added": {k: v for k, v in se.SCEN_CODES.items() if v >= 44},
            "calibration_rules": "02_run_simulation.exact_lensel_gamma (gamma whose mean simulated SFE-006 baseline "
                                 "equals the observed 1.3840, by linear interpolation along the L0grid output) and "
                                 "exact_sfe006_effect over the K/R +/-5 grids extended to b = 8; both fixed before "
                                 "the runs they select",
            "calibrated_values": {"lensel7_gamma": 3.4978345453601807, "kr5_b_log2_or_per_sd": 6.093335419434448},
            "backward_compatibility": "only 'lensel7' tasks add [15, round(gamma*1000)] to their seed; 07_rerun_check "
                                      "re-ran 205 attribute-tables from 22 replicate files with the round-4 engine "
                                      "(max absolute difference 1.8e-15, no mismatching column) and re-generated "
                                      "4,700 stored tables identically",
            "verifier_outputs_read_sha256": {
                p: cb.sha256_of(os.path.join(cb.RESULTS, "verify_r4", p))
                for p in ("lensel_calib_summary.csv", "cells_summary.csv", "calib_points.json",
                          "negb_detmatch_summary.csv", "stress_summary.csv", "stress_match_summary.csv")
                if os.path.exists(os.path.join(cb.RESULTS, "verify_r4", p))}},
        "environment": {"python": sys.version.split()[0], "platform": platform.platform(),
                        "numpy": numpy.__version__, "pandas": pandas.__version__, "scipy": scipy.__version__,
                        "scikit-learn": sklearn.__version__, "matplotlib": matplotlib.__version__},
        "network": "none used",
    }
    with open(os.path.join(cb.RESULTS, "provenance.json"), "w", encoding="utf-8") as fh:
        json.dump(prov, fh, indent=1, ensure_ascii=False)
    print("wrote provenance.json with", len(prov["outputs_sha256"]), "outputs")


if __name__ == "__main__":
    main()
