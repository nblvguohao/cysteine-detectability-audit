"""Step 7 (POST HOC revision analysis, item D_empirical_background): collate the main and sensitivity tables and
write provenance.json (sha256 of every input read, library versions, seeds, sha256 of every output).

Outputs: main_table_kr_de.csv, sensitivity_table.csv, matched_vs_theoretical_paired.csv, provenance.json
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import glob
import json
import os
import platform

import numpy as np
import pandas as pd

from common import (CYS_AUDIT_SRC, IN_CAND, IN_EMP_SUMMARY, IN_FASTA_GZ, IN_HYDN, IN_HYDP, MASTER_SEED, ORIG_EMP_PROTOCOL,
                    ORIG_EMP_SCRIPT, REPO, RESULTS, SCRIPTS, STORED_AUDIT_OBSERVED, STORED_AUDIT_PROTEOME, TOOL_SEED, W,
                    dump_json, sha256_file, sha256_lf)

FEATURES = ["kr_dist", "kr_prox", "de_dist", "de_prox"]
LABELS = {"proteome": "(a) proteome: all other cysteines of identified proteins (Fig. 2a)",
          "theoretical": "(b) theoretically detectable (tryptic 7-30 aa, <=2 missed cleavages)",
          "match_theo": "(b') 1:1 matched on the theoretical flag",
          "empirical": "(c) observed in the unenriched (global) arm",
          "match_emp_lr": "(d) 1:1 matched on the empirical score (logistic, OOF)",
          "match_emp_gbm": "(d') 1:1 matched on the empirical score (gradient boosting, OOF)",
          "observed": "(e) identified in the same ABE arm without a site (Fig. 2b)"}
MATCHED = {"match_theo", "match_emp_lr", "match_emp_gbm", "match_emp_peptide_only", "match_emp_lr_within_protein"}


def main():
    D = pd.read_csv(f"{RESULTS}/background_log2or_by_design.csv").set_index("design")
    P = pd.read_csv(f"{RESULTS}/paired_fraction_removed.csv")
    Rs = pd.read_csv(f"{RESULTS}/rematch_bootstrap_summary.csv")
    RR = pd.read_csv(f"{RESULTS}/rematch_bootstrap_replicates.csv.gz")
    CB = pd.read_csv(f"{RESULTS}/conditional_bootstrap_summary.csv")

    def row_for(design, ps="all"):
        name = design if ps == "all" else f"{design}|{ps}"
        base = D.loc[name]
        r = {"positive_set": ps, "design": design, "label": LABELS.get(design, design),
             "n_positive": int(base.n_positive), "n_background": int(base.n_background)}
        for f in FEATURES:
            pr = P[(P.positive_set == ps) & (P.design == name) & (P.feature == f)].iloc[0]
            if design in MATCHED or design == "match_emp_obs_flag":
                cb = CB[(CB.design == name) & (CB.feature == f)].iloc[0]
                rs = Rs[(Rs.design == name) & (Rs.feature == f)]
                est = cb.estimate_median_200_match_seeds
                r[f"{f}_estimate"] = est
                r[f"{f}_ci95_low"], r[f"{f}_ci95_high"] = cb.cond_ci95_low, cb.cond_ci95_high
                r[f"{f}_ci975_low"], r[f"{f}_ci975_high"] = cb.cond_ci975_low, cb.cond_ci975_high
                r[f"{f}_estimate_primary_match_seed"] = base[f"{f}_estimate"]
                r[f"{f}_match_seed_p2.5"] = cb["match_seed_p2.5"]
                r[f"{f}_match_seed_p97.5"] = cb["match_seed_p97.5"]
                r[f"{f}_fraction_removed"] = cb.fraction_removed_point
                r[f"{f}_fraction_removed_ci95_low"] = cb.fraction_removed_ci95_low
                r[f"{f}_fraction_removed_ci95_high"] = cb.fraction_removed_ci95_high
                r[f"{f}_excess_over_observed"] = cb.excess_over_observed_point
                r[f"{f}_excess_over_observed_ci95_low"] = cb.excess_over_observed_ci95_low
                r[f"{f}_excess_over_observed_ci95_high"] = cb.excess_over_observed_ci95_high
                if len(rs):
                    rs = rs.iloc[0]
                    r[f"{f}_rematch_ci95_low"], r[f"{f}_rematch_ci95_high"] = rs.rematch_ci95_low, rs.rematch_ci95_high
                    r[f"{f}_rematch_boot_median"] = rs.rematch_boot_median
                r[f"{f}_interval_type"] = ("matched-sample protein-clustered bootstrap (5,000) mixed over 200 matchings; "
                                           "point = median over 200 matchings; rematch interval in separate columns")
            else:
                r[f"{f}_estimate"] = base[f"{f}_estimate"]
                r[f"{f}_ci95_low"], r[f"{f}_ci95_high"] = base[f"{f}_ci95_low"], base[f"{f}_ci95_high"]
                r[f"{f}_ci975_low"], r[f"{f}_ci975_high"] = base[f"{f}_ci975_low"], base[f"{f}_ci975_high"]
                r[f"{f}_fraction_removed"] = pr.fraction_removed_point
                r[f"{f}_fraction_removed_ci95_low"] = pr.fraction_removed_ci95_low
                r[f"{f}_fraction_removed_ci95_high"] = pr.fraction_removed_ci95_high
                r[f"{f}_excess_over_observed"] = pr.excess_over_observed_point
                r[f"{f}_excess_over_observed_ci95_low"] = pr.excess_over_observed_ci95_low
                r[f"{f}_excess_over_observed_ci95_high"] = pr.excess_over_observed_ci95_high
                r[f"{f}_interval_type"] = "protein-clustered bootstrap (5,000; Cys-Audit convention); fraction/excess paired"
            r[f"{f}_share_flagged_pos"] = base[f"{f}_share_flagged_pos"]
            r[f"{f}_share_flagged_bg"] = base[f"{f}_share_flagged_bg"]
        r["balance_auc_empirical_score_pos_vs_bg"] = base.get("bal_auc_score_lr_full_pos_vs_bg", np.nan)
        r["balance_smd_logit_empirical_score"] = base.get("bal_smd_logit_score_lr_full", np.nan)
        r["bg_share_identified_in_abe_arm"] = base.get("bal_share_identified_in_arm_bg", np.nan)
        return r

    main_rows = [row_for(d) for d in ("proteome", "theoretical", "match_theo", "empirical", "match_emp_lr",
                                      "match_emp_gbm", "observed")]
    M = pd.DataFrame(main_rows)
    # Cys-Audit status vocabulary (tool rule: Bonferroni 97.5% interval for two bands, margin 0.5 log2 units)
    sys.path.insert(0, CYS_AUDIT_SRC)
    from cys_audit.status import status_from_interval
    for f in ("kr_dist", "kr_prox"):
        M[f"{f}_cys_audit_status_975"] = [status_from_interval(e, lo, hi, 0.0, 0.5)[0] for e, lo, hi in
                                          zip(M[f"{f}_estimate"], M[f"{f}_ci975_low"], M[f"{f}_ci975_high"])]
    M.to_csv(f"{RESULTS}/main_table_kr_de.csv", index=False)

    sens_rows = []
    for d in ("theoretical_symmetric", "theoretical_trypsinP", "empirical_symmetric", "match_emp_peptide_only",
              "match_emp_lr_caliper", "match_emp_obs_flag", "match_emp_lr_within_protein",
              "match_emp_lr_within_protein_caliper", "proteome_same_protein", "observed_same_protein",
              "theoretical_same_protein", "empirical_same_protein", "observed_other_protein", "empirical_other_protein"):
        sens_rows.append(row_for(d))
    for ps in ("ha_strict", "ha_lenient"):
        for d in ("proteome", "theoretical", "match_theo", "empirical", "match_emp_lr", "observed"):
            sens_rows.append(row_for(d, ps))
    S = pd.DataFrame(sens_rows)
    S.to_csv(f"{RESULTS}/sensitivity_table.csv", index=False)

    # paired contrasts between the empirical and theoretical designs of the same family
    pc = []
    for ps in ("all", "ha_strict", "ha_lenient"):
        tag = "" if ps == "all" else f"|{ps}"
        for f in FEATURES:
            cb = CB[(CB.design == f"match_emp_lr{tag} minus match_theo{tag}") & (CB.feature == f)].iloc[0]
            pc.append({"positive_set": ps, "feature": f,
                       "contrast": "match_emp_lr - match_theo (matched-sample bootstrap mixed over 200 matchings)",
                       "point": float(cb.estimate_median_200_match_seeds), "ci95_low": float(cb.cond_ci95_low),
                       "ci95_high": float(cb.cond_ci95_high)})
            v = RR[f"match_emp_lr{tag}|{f}"] - RR[f"match_theo{tag}|{f}"]
            pc.append({"positive_set": ps, "feature": f, "contrast": "match_emp_lr - match_theo (rematch bootstrap)",
                       "point": float(np.median(v)), "ci95_low": float(np.quantile(v, 0.025)),
                       "ci95_high": float(np.quantile(v, 0.975))})
            pr = P[(P.positive_set == ps) & (P.design == f"empirical{tag}") & (P.feature == f)].iloc[0]
            pc.append({"positive_set": ps, "feature": f, "contrast": "empirical - theoretical (restriction; paired bootstrap)",
                       "point": float(pr.diff_vs_theoretical_point), "ci95_low": float(pr.diff_vs_theoretical_ci95_low),
                       "ci95_high": float(pr.diff_vs_theoretical_ci95_high)})
        if ps == "all":
            for f in FEATURES:
                cb = CB[(CB.design == "match_emp_gbm minus match_theo") & (CB.feature == f)].iloc[0]
                pc.append({"positive_set": ps, "feature": f,
                           "contrast": "match_emp_gbm - match_theo (matched-sample bootstrap mixed over 200 matchings)",
                           "point": float(cb.estimate_median_200_match_seeds), "ci95_low": float(cb.cond_ci95_low),
                           "ci95_high": float(cb.cond_ci95_high)})
    PC = pd.DataFrame(pc)
    PC.to_csv(f"{RESULTS}/matched_vs_theoretical_paired.csv", index=False)

    with pd.option_context("display.width", 250, "display.max_columns", 60):
        print(M[["design", "n_positive", "n_background", "kr_dist_estimate", "kr_dist_ci95_low", "kr_dist_ci95_high",
                 "kr_dist_fraction_removed", "kr_dist_fraction_removed_ci95_low", "kr_dist_fraction_removed_ci95_high",
                 "kr_dist_excess_over_observed", "kr_dist_excess_over_observed_ci95_low",
                 "kr_dist_excess_over_observed_ci95_high"]].round(4).to_string())
        print(M[["design", "kr_prox_estimate", "kr_prox_ci95_low", "kr_prox_ci95_high", "de_dist_estimate",
                 "de_dist_ci95_low", "de_dist_ci95_high", "de_prox_estimate", "de_prox_ci95_low", "de_prox_ci95_high",
                 "de_dist_share_flagged_pos", "de_dist_share_flagged_bg"]].round(4).to_string())
        print(S[["positive_set", "design", "n_positive", "n_background", "kr_dist_estimate", "kr_dist_ci95_low",
                 "kr_dist_ci95_high", "kr_dist_fraction_removed", "de_dist_estimate", "de_dist_ci95_low",
                 "de_dist_ci95_high"]].round(4).to_string())
        print(PC.round(4).to_string())

    # ------------------------------------------------------------------ provenance
    import lightgbm, scipy, sklearn, statsmodels
    sys.path.insert(0, CYS_AUDIT_SRC)
    from cys_audit import constants as CAC
    inputs = {
        IN_HYDP: {"role": "trypsin ABE arm +hydroxylamine (site table)"},
        IN_HYDN: {"role": "trypsin ABE arm -hydroxylamine (HA-specific sensitivity)"},
        IN_CAND: {"role": "phase-3 candidate peptides with global-arm identification flag"},
        IN_EMP_SUMMARY: {"role": "stored empirical-check summary (reproduction target)"},
        IN_FASTA_GZ: {"role": "UniProt 2026_03 mouse reference proteome"},
        STORED_AUDIT_PROTEOME: {"role": "stored Fig. 2a trypsin audit (reproduction target)"},
        STORED_AUDIT_OBSERVED: {"role": "stored Fig. 2b trypsin audit (reproduction target)"},
        ORIG_EMP_SCRIPT: {"role": "original empirical-check script (read for the model specification)"},
        ORIG_EMP_PROTOCOL: {"role": "original empirical-check protocol (read for the model specification)"},
        f"{REPO}/scripts/run_phase4_validation_2026-09-22.py": {"role": "original Fig. 2 audit script (read)"},
        "C:/Users/admin/Desktop/小论文/巯基化/MCP/01_manuscript.tex": {"role": "manuscript (read-only; text quoted)"},
        "C:/Users/admin/Desktop/小论文/巯基化/MCP/supplemental/Supplemental_Note_4_four_protease_public_deposit.md":
            {"role": "Supplemental Note 4 (read-only)"},
        f"{W}/inputs/mcp_package/Source_Data_Fig2_public_four_protease.csv": {"role": "Fig. 2 source data (read)"},
        # added in revision round 2 (step 10)
        f"{REPO}/results/phase1_probe_findings_2026-09-22.json":
            {"role": "phase-1 probe record of the deposit's mqpar_Global.xml / mqpar_ABE.xml (search design of the "
                     "arms; read in step 10)"},
        f"{REPO}/scripts/preregister_phase3_empirical_detectability_2026-09-22.py":
            {"role": "phase-3 pre-registration (read for the global arm's search settings: Trypsin/P, min length 7, "
                     "max 4,600 Da; not run)"},
        "C:/Users/admin/Desktop/小论文/巯基化/MCP/supplemental/Supplemental_Data_12_analysis_provenance.csv":
            {"role": "Supplemental Data 12 (read-only; column vocabulary for the proposed row)"},
    }
    for p in inputs:
        inputs[p]["sha256"] = sha256_file(p)
        if p.endswith((".tsv", ".csv", ".json", ".py", ".md", ".tex")):
            inputs[p]["sha256_lf_normalised"] = sha256_lf(p)
    cys_src = {}
    for f in sorted(glob.glob(f"{CYS_AUDIT_SRC}/cys_audit/**/*.py", recursive=True)):
        cys_src[os.path.relpath(f, CYS_AUDIT_SRC).replace("\\", "/")] = sha256_file(f)
    cli = json.load(open(f"{RESULTS}/cli_reproduction.json", encoding="utf-8"))
    outputs = {}
    verifier_files = {}
    for f in sorted(glob.glob(f"{RESULTS}/**/*", recursive=True)):
        rel = os.path.relpath(f, RESULTS).replace("\\", "/")
        if not os.path.isfile(f) or os.path.basename(f).startswith("_") or os.path.basename(f) == "provenance.json":
            continue
        if rel.startswith(("verify_r1/", "verify_r2/", "verify_r3/")):
            # written by the round-1/2/3 verifiers, not by this item's scripts; read only for comparison
            verifier_files[rel] = sha256_file(f)
        else:
            outputs[rel] = sha256_file(f)
    # this item's scripts only (the verify_r1_*.py files in the same folder were written by the verifier)
    scripts = {os.path.basename(f): sha256_file(f) for f in sorted(glob.glob(f"{SCRIPTS}/*.py"))
               if not os.path.basename(f).startswith("verify_")}
    verifier_scripts = {os.path.basename(f): sha256_file(f) for f in sorted(glob.glob(f"{SCRIPTS}/verify_*.py"))}
    prov = {
        "item": "D_empirical_background",
        "status": "POST HOC revision analysis (2026-09-30) in response to review criticism; not registered, not pre-specified",
        "inputs": inputs,
        "fasta_2026_03_decompressed_sha256": cli["fasta_2026_03_decompressed_sha256"],
        "cys_audit_source": {"path": CYS_AUDIT_SRC, "version": CAC.VERSION, "sha256": cys_src,
                             "note": "checks/cleavage.py, stats.py, proteases.py, io.py identical (LF-normalised) to archive/v0.1.0 used for the stored audits"},
        "scripts": scripts,
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
                        "scipy": scipy.__version__, "scikit_learn": sklearn.__version__,
                        "statsmodels": statsmodels.__version__, "lightgbm": lightgbm.__version__,
                        "platform": platform.platform(), "max_worker_processes": 4},
        "seeds": {"cys_audit_tool_seed": TOOL_SEED,
                  "feature_bootstrap_seeds": {"kr_prox": TOOL_SEED + 101, "kr_dist": TOOL_SEED + 102,
                                              "de_prox": TOOL_SEED + 103, "de_dist": TOOL_SEED + 104},
                  "cv_fold_assignment": "numpy default_rng(20260922 + 1) permutation of sorted proteins, fold = rank mod 5 (original)",
                  "lightgbm_random_state": MASTER_SEED, "paired_bootstrap_seed": MASTER_SEED + 7,
                  "primary_matching_seed": MASTER_SEED + 11, "matching_seed_family": "20260941 + 1000 + k, k = 0..199",
                  "rematch_bootstrap_seed": "default_rng([20260941, 7, replicate])", "mh_bootstrap_seed": MASTER_SEED + 13,
                  "pythonhashseed_cli": 0,
                  "step9_after_verification": {
                      "primary_paired_seed": MASTER_SEED + 7,
                      "alternative_paired_seeds": [MASTER_SEED + 2000 + k for k in range(20)],
                      "big_run_seed": MASTER_SEED + 3000, "big_run_replicates": 50000,
                      "detection_link_bootstrap_seeds": [MASTER_SEED + 4000 + i for i in range(4)],
                      "per_matching_status_seeds": {"kr_dist": TOOL_SEED + 102, "kr_prox": TOOL_SEED + 101},
                      "standardised_bprime_seed": MASTER_SEED + 7},
                  "step10_after_verification_round2": {
                      "primary_paired_seed": MASTER_SEED + 7, "big_run_seed": MASTER_SEED + 3000,
                      "big_run_replicates": 50000,
                      "stored_matchings_bprime_d": "20260941 + 1000 + k, k = 0..199 (reproduced)",
                      "dtheo_matchings": f"{MASTER_SEED + 5000} + k, k = 0..199",
                      "per_design_and_per_matching_status_seeds": {"kr_dist": TOOL_SEED + 102,
                                                                   "kr_prox": TOOL_SEED + 101,
                                                                   "de_dist": TOOL_SEED + 104,
                                                                   "de_prox": TOOL_SEED + 103}},
                  "step11_after_verification_round3": {
                      "run_A_multiplicities_seed": MASTER_SEED + 7,
                      "run_A_matching_generators": "default_rng([20260941, 7, replicate]) (step-04 convention; one "
                                                   "fresh generator per matched design)",
                      "run_B_multiplicities_seed": MASTER_SEED + 6000,
                      "run_B_matching_generators": "default_rng([20266941, 7, replicate])",
                      "replicates_per_run": 5000}},
        "bootstrap_replicates": 5000,
        "outputs_sha256": outputs,
        "run_order": "00, 01, 02, 03, 04, 04b, 05, 06, 08, 09, 10, 10b, 11, 07 (07 collates and writes this file last)",
        "revision_after_verification": {
            "date": "2026-09-30",
            "note": ("step 9 (09_revision_after_verification.py) was added after adversarial verification round 1; "
                     "all its quantities were chosen after every earlier result of this item had been seen"),
            "verifier_files_read_for_comparison_sha256": verifier_files,
            "verifier_scripts_sha256_not_run_by_this_item": verifier_scripts},
        "revision_round2": {
            "date": "2026-09-30",
            "note": ("step 10 (10_round2_combined_designs.py) was added after adversarial verification round 2: the "
                     "combined designs (theoretical AND empirical restriction; theo x emp standardisation; score "
                     "matching within theoretical strata), the candidate-aware empirical standardisation, the "
                     "decomposition of the (c) residual and the global arm's search design. All are post hoc. The "
                     "verifier's scripts were not imported or run; its stored values (verify_r2/) were read for the "
                     "comparison in r2_verifier_comparison.csv. Step 10 reproduces the stored step-9 contrasts and "
                     "per-design intervals exactly (self-checks in r2_summary.json)."),
            "outputs": sorted(k for k in outputs if k.startswith("r2_"))},
        "revision_round3": {
            "date": "2026-09-30",
            "note": ("step 11 (11_round3_rematch_fprime.py) was added after adversarial verification round 3: the "
                     "re-matching bootstrap of step 04 extended to (f') (score matching within theoretical strata) and "
                     "to paired differences, in the step-04 multiplicities (run A, which reproduces the stored step-04 "
                     "replicates of (a), (e), (b') and (d)) and in an independent replication run (run B). Post hoc. "
                     "The verifier's scripts were not imported or run; its stored values (verify_r3/) were read for "
                     "the comparison in r3_summary.json."),
            "outputs": sorted(k for k in outputs if k.startswith("r3_"))},
    }
    dump_json(prov, f"{RESULTS}/provenance.json")
    print("provenance written;", len(outputs), "outputs hashed")


if __name__ == "__main__":
    main()
