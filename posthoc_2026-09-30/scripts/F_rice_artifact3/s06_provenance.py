"""Step 6: provenance.json for item F_rice_artifact3 (POST HOC revision analysis, 2026-09-30).

Records sha256 of every input read, of every script and output of this item, library versions,
seeds and the network use (UniProt REST only, from fetched/fetch_summary.json).
"""
import platform

import numpy, pandas, scipy, statsmodels, sklearn
import pandas as pd

import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import *  # noqa

inputs = [PEP, PRO, FASTA_REF, STORED_A3, STORED_A3_SENS, REPO_SITE_SCORES, REPO_EFFECT_SIZES,
          REPO_RECON_ARMS, REPO_A3_SCRIPT,
          REPO / "scripts" / "artefact3_effect_sizes_and_matching_2026-09-19.py",
          REPO / "scripts" / "artefact3_positive_set_sensitivity_2026-09-19.py",
          REPO / "scripts" / "artefact3_reconstructed_original_arms_2026-09-19.py",
          REPO / "results" / "artefact3_osa_abundance_recomputed_2026-09-19_audit.json",
          REPO / "results" / "artefact3_reconstructed_original_arms_2026-09-19_audit.json",
          REPO / "scripts" / "run_retrospective_external_corrected.py",
          REPO / "scripts" / "run_v2_external_retrospective.py",
          REPO / "scripts" / "preregister_self_audit_public_2026-09-20.py",
          REPO / "scripts" / "run_self_audit_public_cohorts_2026-09-20.py",
          REPO / "scripts" / "baseline_phase0_rerun_correction_2026-09-21.py",
          REPO / "baseline" / "number_registry_2026-09-21.json",
          REPO / "supplemental" / "Supplemental_Note_2_public_ranker_training.md",
          CYS_AUDIT_SRC / "cys_audit" / "proteases.py",
          FETCHED / "fetched_sequences.tsv",
          pathlib.Path("/path/to/local/巯基化/MCP/01_manuscript.tex"),
          pathlib.Path("/path/to/local/巯基化/MCP/supplemental/Supplemental_Data_7_dataset_census.csv"),
          # read in the revision after adversarial verification (s08, s08b, s09 and the report)
          REPO / "baseline" / "baseline_metrics.csv",
          REPO / "baseline" / "baseline_report.md",
          REPO / "baseline" / "rerun_correction_2026-09-21_audit.json",
          REPO / "baseline" / "rerun_comparison_2026-09-21.csv",
          REPO / "baseline" / "manuscript_text_findings_2026-09-21.csv",
          REPO / "results" / "self_audit_public_cohorts_2026-09-20_audit.json",
          REPO / "results" / "artifact_manifest.csv",
          REPO / "reports" / "LABEL_PROVENANCE_AND_SEMANTICS.md",
          REPO / "results" / "artefact3_effect_sizes_2026-09-19.csv",
          pathlib.Path("/path/to/local/巯基化/MCP/supplemental/Supplemental_Note_2_public_ranker_training.md"),
          pathlib.Path("/path/to/local/巯基化/MCP/supplemental/Supplemental_Note_9_public_cleavage_geometry_self_audit.md"),
          pathlib.Path("/path/to/local/巯基化/MCP/supplemental/Supplemental_Note_10_public_five_protease_self_audit.md"),
          pathlib.Path("/path/to/local/巯基化/MCP/supplemental/Supplemental_Note_11_public_feature_ablation.md"),
          pathlib.Path("/path/to/local/巯基化/MCP/supplemental/Supplemental_Note_12_public_detectability_only.md"),
          pathlib.Path("/path/to/local/巯基化/MCP/supplemental/Supplemental_Note_13_nonnull_permutation_check.md"),
          pathlib.Path("/path/to/local/巯基化/MCP/supplemental/Supplemental_Data_12_analysis_provenance.csv"),
          # read in revision round 3 (report text: current cover-letter bullet, Note numbering)
          pathlib.Path("/path/to/local/巯基化/MCP/02_cover_letter.md"),
          pathlib.Path("/path/to/local/巯基化/MCP/supplemental/Supplemental_Note_20_fig5a_absolute_discrimination.md")]
fs = json.loads((FETCHED / "fetch_summary.json").read_text(encoding="utf-8"))
prov = {
    "item": ITEM,
    "label": "POST HOC revision analysis in response to pre-submission review criticism; not registered, not pre-specified",
    "date": "2026-09-30",
    "python": sys.version, "platform": platform.platform(),
    "libraries": {"numpy": numpy.__version__, "pandas": pandas.__version__, "scipy": scipy.__version__,
                  "statsmodels": statsmodels.__version__, "scikit-learn": sklearn.__version__},
    "seeds": {"posthoc_resampling_and_simulation": SEED,
              "reproduction_of_stored_run": SEED_ORIGINAL,
              "s04_bootstrap_per_combination": "SEED + 1..12",
              "s04b_matching_bootstrap": "SEED + 7",
              "s08_matched_pairs_bootstrap": "SEED + 801 (A), 802 (S), 803 (S1), 804 (S_map); 400 resamples",
              "s08b_calibrated_nulls": "numpy SeedSequence(entropy = SEED + 810 + k, spawn_key = (chunk,)), "
                                       "k = 0..4 over (S,N1),(S,N2),(S_map,N1),(S_map,N2),(A,N2); 4 chunks x 250",
              "s10_mirror_matched_bootstrap": "SEED + 831 (A), 832 (S), 833 (S1), 834 (S_map); 400 resamples",
              "s10b_nulls": "numpy SeedSequence(entropy, spawn_key = (chunk,)), 4 chunks x 250; the five s08b "
                            "combinations keep SEED + 810..814 (label sets identical to s08b, asserted); new: "
                            "(A,N1) SEED+815, (S,N3) 816, (S,N4) 817, (S_map,N3) 818, (S_map,N4) 819, (A,N3) 820, (A,N4) 821",
              "s12_null_extension": "numpy SeedSequence(entropy = SEED + 840 + i, spawn_key = (chunk,)), i = index of "
                                    "the combination in s10b's COMBOS (0..11), 16 chunks x 250 = 4,000 further sets; "
                                    "pooled with the 1,000 s10b sets"},
    "revision_round_3": {
        "date": "2026-09-30",
        "reason": "verifier round 3: (major) proposed Results text dropped the manuscript's post hoc disclosure "
                  "(text only); (noted in the verifier's comparison) the count of contrasts outside the 1,000-set "
                  "null range differed by one from an independent re-simulation for (A,N4) and (A,N1); step 0 of the "
                  "rice chain skipped the 22 decoy rows without a start position",
        "new_scripts": ["s12_null_extension.py", "s13_indicator_test.py"],
        "changed_scripts": ["s11_rice_chain_start.py (decoy rows without a start position counted at step 0; new "
                            "step 0b; steps 1-5 unchanged, asserted; records that none of the 400 lost accessions' "
                            "saved UniProtKB entry records carries a sequence)", "s06_provenance.py", "s07_key_numbers.py"],
        "run_order": "s12_null_extension -> s13_indicator_test -> s11_rice_chain_start -> s07_key_numbers -> s06_provenance",
        "declared_before_s12_ran": "pooled 5,000-set summaries become the primary null results whatever the outcome; "
                                   "'edge' = pooled two-sided P in [0.04, 0.06]",
        "network": "none in round 3",
        "verifier_files_in_results": "verify_r3/ was written by the round-3 verifier, not by this item's scripts"},
    "revision_round_2": {
        "date": "2026-09-30",
        "reason": "verifier round 2, major 1: detectable (n_det) versus total (n_cys) cysteine count",
        "new_scripts": ["s10_lib.py", "s10_count_vs_detectability.py", "s10b_null_total.py", "s11_rice_chain_start.py"],
        "changed_scripts": ["s06_provenance.py", "s07_key_numbers.py"],
        "run_order": "s10_count_vs_detectability -> s10b_null_total -> s11_rice_chain_start -> s07_key_numbers -> s06_provenance",
        "determinism_check": "all 14 outputs of s10/s10b/s11 (csv/json) byte-identical when rerun with PYTHONHASHSEED 12345 and 999",
        "network": "none in round 2",
        "fitting_note": "s10: logistic GLMs that IRLS stops at its 100-iteration cap (the unconstrained cr() spline "
                        "basis with a categorical count term is collinear with the intercept) are refitted by "
                        "Newton-Raphson; the s08b/s10b categorical-spline nulls (N2, N4) keep the IRLS fit so that "
                        "s08b is reproduced exactly; their fitted probabilities differ from the converged Newton "
                        "fit by at most 0.0017 (expected positives within 1.8 of observed)",
        "verifier_files_in_results": "verify_r2/ was written by the round-2 verifier, not by this item's scripts"},
    "revision_after_verification": {
        "date": "2026-09-30",
        "new_scripts": ["s08_lib.py", "s08_residual.py", "s08b_null.py", "s08c_undet_designs.py",
                        "s09_trace_artifact2.py"],
        "run_order": "s08_residual -> s08b_null -> s08c_undet_designs -> s09_trace_artifact2 -> s07_key_numbers -> s06_provenance",
        "determinism_check": "all 16 outputs of s08/s08b/s09 byte-identical when rerun with PYTHONHASHSEED 12345/999/7",
        "changed_scripts": ["s06_provenance.py", "s07_key_numbers.py"],
        "verifier_files_in_results": "verify_r1/ was written by the round-1 verifier (not by this item's scripts); "
                                     "hashed below with the other files of the results folder"},
    "inputs_sha256": {str(p): (sha256_file(p) if p.exists() else "MISSING") for p in inputs},
    # the internal repository is a Windows checkout with CRLF line endings; hashing its scripts
    # after CRLF->LF normalisation reproduces the sha256 recorded in their own audit files
    "repo_scripts_sha256_LF_normalised": {
        str(p): hashlib.sha256(p.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
        for p in inputs if str(p).startswith(str(REPO)) and p.suffix == ".py" and p.exists()},
    "network": {"allowed": "rest.uniprot.org only", "requests": fs["requests"],
                "bytes_downloaded": fs["bytes_downloaded"],
                "n_requested": fs["n_requested"], "n_active_uniprotkb": fs["n_active_uniprotkb"],
                "n_uniparc": fs["n_uniparc"], "n_not_found": fs["n_not_found"],
                "unisave_spot_check_identical": sum(1 for s in fs["unisave_spot_check"]
                                                    if s.get("identical_to_uniparc_highest_seq_version")),
                "unisave_spot_check_n": len(fs["unisave_spot_check"]),
                "s01b_inactive_reasons_400": (json.loads((FETCHED / "inactive_reasons_400.json").read_text(encoding="utf-8"))
                                              if (FETCHED / "inactive_reasons_400.json").exists() else None),
                "s01c_unisave_last_release_sample": (
                    pd.read_csv(FETCHED / "unisave_last_release_sample.csv")["lastRelease"].value_counts().to_dict()
                    if (FETCHED / "unisave_last_release_sample.csv").exists() else None)},
    "scripts_sha256": {p.name: sha256_file(p) for p in sorted(SCRIPTS.glob("*.py"))},
    "outputs_sha256": {str(p.relative_to(RES)): sha256_file(p) for p in sorted(RES.rglob("*"))
                       if p.is_file() and p.name != "provenance.json"},
    "worker_processes": "1 (all steps) except s08b_null.py, s10b_null_total.py and s12_null_extension.py: 4 processes each",
}
write_json(RES / "provenance.json", prov)
print("provenance written:", len(prov["inputs_sha256"]), "inputs,", len(prov["outputs_sha256"]), "outputs")
