# Release notes - v3.1.1 (2 October 2026)

Release candidate accompanying the Molecular & Cellular Proteomics submission (manuscript v3.1.1).
Changes from v3.1.0 (commit b94aab5). Cys-Audit itself is unchanged (0.2.2; 55 tests).

## 1. Added

### `artifacts12/` - scripts and stored outputs behind Artifacts 1 and 2

Copied unchanged from the collaborator working tree on 2 October 2026. All twelve files whose SHA-256 was
recorded at staging on 23 September 2026 (five scripts, six outputs and `MANIFEST.json`) match on their raw
bytes; no line-ending normalisation was needed. `artifacts12/README.md` gives the SHA-256 of every file, the
manuscript numbers each supports, how to re-run, the inputs that are not redistributed (PlantPTMViewer
exports, dbPTM files, PRIDE tables, UniProt tables, AlphaFold models; names, dates and SHA-256), and the
known gaps: the Arabidopsis proteome table read by `09_map_large_scale.py` (not retained at its recorded
path), the producer of `35_symmetric_ncys_test.json` (not retained), and the producers of
`24_extra_sequences.json` and of the dbPTM/PXD044043 site tables (not retained).
`artifacts12/number_registry.csv` maps each printed number to its file and JSON path.

### Files the manuscript, Supplemental Data and Source Data cite that v3.1.0 lacked

- `protocols/detectability_baseline_plmsnosite_preregistration_2026-09-24.json` (Supplemental Note 5;
  Supplemental Data 12, row pLMSNOSite).
- Figure 7 score tables named in `source_data_submitted/Source_Data_Fig7_self_audit_public.csv`:
  `results/self_audit_human_trained_rice_2026-09-23.csv` (+ `_audit.json`),
  `results/self_audit_human_trained_rice_five_proteases_2026-09-23.csv`,
  `results/self_audit_rice_trained_human_2026-09-24.csv` (+ `_audit.json`),
  `results/self_audit_rice_trained_human_five_proteases_2026-09-24.csv`,
  `results/public_refit_ablation_lgb_rank_2026-09-24.csv` (+ `_audit.json`).
- `results/search_space_mass_accuracy.csv` (+ `_audit.json`), cited by the Figure 3 Source Data.
- The scripts that wrote them (SHA-256 equal to the `script_sha256` recorded in their audit files):
  `run_self_audit_human_trained_rice_2026-09-23.py`, `run_self_audit_rice_trained_human_2026-09-24.py`,
  `reanalyse_public_refit_lgb_member_2026-09-24.py`, `analyse_public_refit_2026-09-21.py`,
  `train_v2_public_human_2026-09-23.py`, `train_v2_public_rice_2026-09-24.py`,
  `run_detectability_baseline_plmsnosite_2026-09-24.py` (`audit_search_space_and_mass_accuracy.py` was
  already present).

### Local modules imported by released scripts (previously missing)

`common.py`, `v2_stack.py`, `v2_members.py`, `v2_apply.py`, `v2_features.py`,
`run_cross_protease_detectability_probe.py`, `phase2_claim_cohorts.py`, `phase2b_claim_cohorts.py`,
`phase2d_claim_cohorts.py`, `phase2e_claim_cohorts.py`, `audit_cleaning_and_grouping.py`,
`ptm_detectability_diagnostics.py`, `run_qtrp_trained_model.py`, `run_v3_qtrp_stack.py`,
`run_detection_matched_persulfidation_test.py`, `run_qtrp_intervals_and_cross_label.py`,
`build_structural_chemistry_features.py`, `run_phase2b_claims_backfill.py` (Python imports), and
`run_self_audit_cleavage_flank.py`, `recompute_artefact3_osa_abundance_2026-09-19.py`,
`artefact3_effect_sizes_and_matching_2026-09-19.py`, `ingest_ptm_census_sites.py`,
`nc_figure_style_2026-09-21.py`, `plot_manuscript_f1_three_axes.py`,
`check_figure_text_overlaps_2026-09-19.py` (loaded by path). All are the authors' own code, copied
unchanged from the analysis tree; every local import in `scripts/` now resolves. Where the internal
artifact manifest records a SHA-256 for these files, the copies match it. `v2_apply.py` keeps an
environment-variable default that points to an earlier project's `model_pipeline.py` (not released); see
README, Scope.

### Other

- `figures_r31/README.md` (what the R31 set is, why it has no Figure 5, where the submitted figures are).
- `archive/README.md`, `RELEASE_NOTES_v3.1.1.md`, `.gitattributes` (`* -text`: no end-of-line conversion on
  checkout, so the manifest holds on every platform).

## 2. Updated

- `supplemental/` and `source_data_submitted/` synchronised with the v3.1.1 manuscript package:
  Supplemental Note 6 (Table S6.1 placeholders filled), Supplemental Note 7 (YAP1C check),
  Supplemental Data 12 (pLMSNOSite row), Source Data Figure 3 (`repo/results/` -> `results/`) and Figure 6
  (`source_table` -> `revision_2026-09-30/results/figures_final/Fig6/fig6_final_2026-09-30.py`). Other files
  of the package differ from the repository copies only in line endings.
- `CITATION.cff`: version v3.1.1, released 2026-10-02, authors Guohao Lv and Lichuan Gu (School of
  Artificial Intelligence, Anhui Agricultural University). `cys-audit/CITATION.cff`: the same authors in
  place of the placeholder.
- `cys-audit/LICENSE`: the placeholder ("licence not yet chosen") replaced by the MIT licence the manuscript
  and the root `LICENSE` state; copyright holder in both licences: Guohao Lv and Lichuan Gu.
- `recoding_audit/apply_recoded_count.py` (m10): paths resolved from the repository root, Fig. 1 Source Data
  read from `source_data_submitted/`, idempotent (on the released tables it changes nothing), LF output.
- `recoding_audit/{A,B}/*.py`: absolute Windows paths replaced by repository-relative paths. Note: run
  against the released Supplemental Data 9, which already carries the three full-text re-codings, the
  agreement scripts report 48/74 (A) and 40/74 (B) instead of the stored 46/74 and 37/74; the stored reports
  were computed against the original classes (`background_class_original` where filled).
- `README.md`: contents table, licence line, and a Scope section that states what re-runs and what does not.
- `MANIFEST_SHA256.txt` regenerated.

## 3. Moved

`scripts/baseline_phase0_{script_integrity,rerun,rerun_correction,findings_addendum}_2026-09-21.py` ->
`archive/scripts/` (internal number audits that read the internal and collaborator trees by absolute path;
the addendum imports an unreleased module). Their findings on Artifacts 1-2 were re-checked from the
deposited files (`artifacts12/README.md`).

## 4. Removed (internal or unrelated material)

None of these files is cited by the manuscript, the Supplemental Notes or Data, the Source Data, Cys-Audit,
or any README; no remaining file imports them (checked by name across all of these before removal). They
remain in the history of tag v3.1.0.

| file | reason |
|---|---|
| `scripts/build_formal_answer_part2.py` | Chinese-language report for the collaborating tomato-persulfidation project (gene-family question for that group); unrelated to the manuscript |
| `scripts/build_supplement_zhang.py` | "Supplement for Prof. Zhang Hua": experimental plan and data requests for the collaborator (Boltz-2 / WRKY interaction); unrelated |
| `scripts/run_af3_result_analysis.py` | AlphaFold 3 analysis of a three-protein interaction question (WRKY project); unrelated |
| `scripts/run_interface_proximity_check.py` | WRKY x PP2A-B56 AF3 interface distances; unrelated |
| `scripts/run_user_interaction_prep.py` | preparation for the same interaction question; unrelated |
| `scripts/run_user_motif_analysis.py` | short-linear-motif scan for the interaction/phosphorylation questions of the same project; unrelated |
| `scripts/run_kiae271_external_audit.py` | zero-shot scoring of the collaborator's own persulfidome with the earlier deployed rankers; not used in the manuscript |
| `scripts/run_v2_user_protein_report.py` | per-protein report for a collaborator's 15 proteins; not used; needs an unreleased encoder |
| `scripts/run_v2_candidates.py` | re-scoring of wet-lab candidate proteins; not used; needs an unreleased encoder |
| `scripts/score_experimental_candidates.py` | scoring of supplied candidate proteins for wet-lab follow-up; not used |
| `scripts/audit_manuscript_self_serving_reading_2026-09-19.py` | internal review of how the draft reads; drafting material |
| `scripts/preregister_phase6a_manuscript_2026-09-22.py` | registration of manuscript text corrections; drafting material |
| `protocols/phase6a_manuscript_preregistration_2026-09-22.json` | its protocol (text corrections R24 -> R25); drafting material |
| `protocols/science_session_plan_2026-09-15.json` | internal work-session plan (in Chinese); not an analysis protocol |
| `scripts/plot_fig1_biorender_r35_2026-09-23.py` | superseded BioRender draft of Figure 1 |
| `scripts/plot_fig1_biorender_r36_2026-09-23.py` | superseded BioRender draft of Figure 1 (a comment in `revision_2026-09-30/results/figures_final/Fig1/fig1_final.py` still names it as the origin of a mathtext setting) |
| `scripts/assemble_figures_r35_2026-09-23.py` | assembled the superseded R35 figure folder (not released) |
| `scripts/assemble_figures_r36_2026-09-23.py` | assembled the superseded R36 figure folder (not released) |
| `scripts/verify_reference_coverage_2026-09-18.py` | checked the R2 draft's reference list; drafting material |
| `scripts/verify_phase2_report_numbers_2026-09-22.py` | "house rule" number back-check of an internal phase report that is not released |
| `scripts/verify_phase3_report_numbers_2026-09-22.py` | the same for another internal phase report |
| `scripts/write_manual_fetch_list.py` | internal to-do list of papers to fetch by hand (in Chinese) |
| `scripts/write_manual_fetch_intake.py` | internal intake log of those papers (in Chinese) |

## 5. Line endings and the manifest

All text files outside `artifacts12/` are committed with LF line endings (as in v3.1.0). Six CSV files added
here were written with CRLF by Python's `csv` module and were normalised to LF; where an audit file records
the SHA-256 of the original, that value refers to the CRLF bytes:

| file | SHA-256 as written (CRLF) | SHA-256 committed (LF) |
|---|---|---|
| `results/public_refit_ablation_lgb_rank_2026-09-24.csv` | `fbd556a8258f769d5b8599619dd0d534e50930ce78bb5d0f33c9151b3e169069` | `6c044e3c32db24b5fa8594b4b56f031cefa70f897a107e0bb1eb474db2ad1446` |
| `results/self_audit_rice_trained_human_2026-09-24.csv` | `508e5009d2376913f1fcce080f8db125b005ed1c003ed3fc35f28f373dd6aa6a` | `b90f39114ce457351a8ef060c51dfb7aff0931d4d33c1686c02beddab7660508` |
| `results/search_space_mass_accuracy.csv` | `5984425bd851b25b272e17376a691a65318728f5778d006c431405d8296613d7` | `e5fcbcb4af4e53481524cba4201c8064e7ed7325ebc921c5d8273cd3748b7756` |
| `results/self_audit_rice_trained_human_five_proteases_2026-09-24.csv` | `acaa9c697d9c9972fba4617ad64d4b051bac21d96da3fdc8be9f40702e1da2ba` | `b429e8b15024944e47f5cdf8660d38ea218d9e67370ff7ac4d261ed253b35920` |
| `results/self_audit_human_trained_rice_2026-09-23.csv` | `5b45a5e18251925d4997f9df990647adce68b472d6a13273e199c98c52f68887` | `3dcbb3421348a8e6dddf1b879cbc427777b4dfa6cda133d0adb890218242e439` |
| `results/self_audit_human_trained_rice_five_proteases_2026-09-23.csv` | `4e18a665d857893a7fa06ab410a306241eae24b33f05dedb96409264b84edb06` | `a3704531f643b86591671ee45a13f7a2a888a58b1617021f3ec2f712569e39f0` |

`artifacts12/` is deposited byte for byte, so its three CRLF CSV files keep their line endings (`21_crossspecies_summary.csv`,
whose staged SHA-256 is of those bytes, `36_ncys_decompose.csv`, and the positional-profile CSV). With
`.gitattributes` set to `* -text`, every checkout reproduces the committed bytes, and
`MANIFEST_SHA256.txt` hashes exactly those bytes.

## 6. Files that keep absolute paths

Left unchanged because their SHA-256 values are recorded in provenance or audit files (editing them would
break that record), except where noted above. Mapping for re-running: `W` / `REV` -> `revision_2026-09-30/`;
`REPO` -> repository root; `MCP` -> repository root (`source_data/` -> `source_data_submitted/`);
`ROOT = /Users/.../2026-09-12_persulfidation_upgrade` -> repository root.

Paths in executed code (the script needs the constant changed before it runs):

- `archive/scripts/baseline_phase0_rerun_2026-09-21.py`
- `artifacts12/scripts/09_map_large_scale.py`
- `revision_2026-09-30/results/figures_final/Fig1/fig1_final.py`
- `revision_2026-09-30/results/figures_final/Fig1/qa_fig1_final.py`
- `revision_2026-09-30/results/figures_final/Fig3/fig3_final.py`
- `revision_2026-09-30/results/figures_final/Fig3/verify_fig3_final.py`
- `revision_2026-09-30/results/figures_final/Fig5/variant_SNO-012_undecidable/recolour_fig5c_sno012.py`
- `revision_2026-09-30/results/figures_final/Fig6/fig6_final_2026-09-30.py`
- `revision_2026-09-30/results/figures_final/Fig6/verify_fig6_final.py`
- `revision_2026-09-30/results/final_notes/Supplemental_Note_15_check_numbers.py`
- `revision_2026-09-30/scripts/B_overadjust_sim/common_b.py`
- `revision_2026-09-30/scripts/C_pxd063463_specific/c00_explore.py`
- `revision_2026-09-30/scripts/C_pxd063463_specific/c_common.py`
- `revision_2026-09-30/scripts/D_empirical_background/07_summarize.py`
- `revision_2026-09-30/scripts/D_empirical_background/common.py`
- `revision_2026-09-30/scripts/E_bootstrap_coverage/calibrated_verdicts.py`
- `revision_2026-09-30/scripts/E_bootstrap_coverage/claims_matched_sim.py`
- `revision_2026-09-30/scripts/E_bootstrap_coverage/claims_table.py`
- `revision_2026-09-30/scripts/E_bootstrap_coverage/multiplicity.py`
- `revision_2026-09-30/scripts/E_bootstrap_coverage/provenance.py`
- `revision_2026-09-30/scripts/E_bootstrap_coverage/r1_b1_null.py`
- `revision_2026-09-30/scripts/E_bootstrap_coverage/r1_rc_yardsticks.py`
- `revision_2026-09-30/scripts/E_bootstrap_coverage/r1_sfe008_sim.py`
- `revision_2026-09-30/scripts/E_bootstrap_coverage/r1_verdicts.py`
- `revision_2026-09-30/scripts/E_bootstrap_coverage/rule_difference.py`
- `revision_2026-09-30/scripts/E_bootstrap_coverage/run_grid.py`
- `revision_2026-09-30/scripts/E_bootstrap_coverage/simlib.py`
- `revision_2026-09-30/scripts/E_bootstrap_coverage/summarize_grid.py`
- `revision_2026-09-30/scripts/E_bootstrap_coverage/tables.py`
- `revision_2026-09-30/scripts/E_bootstrap_coverage/validate.py`
- `revision_2026-09-30/scripts/F_rice_artifact3/common.py`
- `revision_2026-09-30/scripts/F_rice_artifact3/s06_provenance.py`
- `revision_2026-09-30/scripts/F_rice_artifact3/verify_r4_build.py`
- `revision_2026-09-30/scripts/F_rice_artifact3/verify_r4_compare.py`
- `revision_2026-09-30/scripts/F_rice_artifact3/verify_r4_extra.py`
- `revision_2026-09-30/scripts/F_rice_artifact3/verify_r4_models.py`
- `revision_2026-09-30/scripts/F_rice_artifact3/verify_r4_null_recount.py`
- `revision_2026-09-30/scripts/F_rice_artifact3/verify_r4_null_sim.py`
- `revision_2026-09-30/scripts/F_rice_artifact3/verify_r4_rice.py`
- `revision_2026-09-30/scripts/G_plmsnosite_rescore/00_inspect_models.py`
- `revision_2026-09-30/scripts/G_plmsnosite_rescore/01_inspect_lambda.py`
- `revision_2026-09-30/scripts/G_plmsnosite_rescore/10_detectability_baseline.py`
- `revision_2026-09-30/scripts/G_plmsnosite_rescore/11_repro_diagnostics.py`
- `revision_2026-09-30/scripts/G_plmsnosite_rescore/20_run_plmsnosite.py`
- `revision_2026-09-30/scripts/G_plmsnosite_rescore/33_within_protein_mean.py`
- `revision_2026-09-30/scripts/G_plmsnosite_rescore/common.py`
- `revision_2026-09-30/scripts/H_artifact5_docs/h5_lib.py`
- `revision_2026-09-30/scripts/I_fig5a_estimand/common.py`
- `revision_2026-09-30/scripts/I_fig5a_estimand/s03_input_availability.py`
- `revision_2026-09-30/scripts/I_fig5a_estimand/verify_r1_rerun/common.py`
- `revision_2026-09-30/scripts/I_fig5a_estimand/verify_r1_rerun/s03_input_availability.py`
- `revision_2026-09-30/scripts/J_provenance_species/j5_singlecys_construction.py`
- `revision_2026-09-30/scripts/J_provenance_species/jcommon.py`
- `revision_2026-09-30/scripts/J_provenance_species/verify_r4_extra.py`
- `scripts/build_release_additions_r31_2026-09-23.py`
- `scripts/run_pxd015307_research_mass_accuracy.py`
- `scripts/run_retrospective_external_corrected.py`
- `scripts/v2_apply.py`
- `scripts/write_phase2e_report.py`
- `scripts/write_phase3_report.py`
- `scripts/write_self_audit_five_proteases_report.py`

Paths only in comments or docstrings (interpreter used, usage line; execution unaffected):

- `artifacts12/positional_profile/scripts/positional_kr_profile_public_2026-09-17.py`
- `scripts/check_figure_text_overlaps_2026-09-19.py`
- `scripts/gate_pka_instrument_equivalence.py`
- `scripts/plot_fig1_overview_2026-09-21.py`
- `scripts/plot_figs_2to5_nc_style_2026-09-21.py`
- `scripts/run_detectability_baseline_plmsnosite_2026-09-24.py`
- `scripts/run_phase2b_claims_backfill.py`
- `scripts/run_pka_ensemble_propka.py`
- `scripts/run_pka_round_verdict.py`
- `scripts/train_v2_public_human_2026-09-23.py`
- `scripts/train_v2_public_rice_2026-09-24.py`

Paths recorded as data (provenance and audit JSON, logs, reports; historical records, left as written):

- `README.md`
- `artifacts12/MANIFEST.json`
- `artifacts12/README.md`
- `protocols/pxd015307_research_preregistration_2026-09-17.json`
- `reports/PXD015307_RESEARCH_SEARCH_SPACE.md`
- `reports/V3_PREREGISTRATION_AND_COMPONENTS.md`
- `revision_2026-09-30/inputs/repo_reports/PXD015307_RESEARCH_SEARCH_SPACE.md`
- `revision_2026-09-30/inputs/repo_reports/V3_PREREGISTRATION_AND_COMPONENTS.md`
- `revision_2026-09-30/inputs/repo_reports/pxd015307_research_preregistration_2026-09-17.json`
- `revision_2026-09-30/inputs/repo_results/pxd015307_research_audit.json`
- `revision_2026-09-30/results/B_overadjust_sim/provenance.json`
- `revision_2026-09-30/results/C_pxd063463_specific/provenance.json`
- `revision_2026-09-30/results/C_pxd063463_specific/verify_r3/provenance_verify_r3_d_rule_difference.json`
- `revision_2026-09-30/results/D_empirical_background/provenance.json`
- `revision_2026-09-30/results/D_empirical_background/r2_global_arm_search_design.json`
- `revision_2026-09-30/results/E_bootstrap_coverage/provenance.json`
- `revision_2026-09-30/results/F_rice_artifact3/provenance.json`
- `revision_2026-09-30/results/F_rice_artifact3/s04_stdout.txt`
- `revision_2026-09-30/results/G_plmsnosite_rescore/provenance.json`
- `revision_2026-09-30/results/H_artifact5_docs/provenance.json`
- `revision_2026-09-30/results/I_fig5a_estimand/feature_identity_check.json`
- `revision_2026-09-30/results/I_fig5a_estimand/provenance.json`
- `revision_2026-09-30/results/J_provenance_species/provenance.json`
- `revision_2026-09-30/results/figures_final/Fig6/build_log.txt`
