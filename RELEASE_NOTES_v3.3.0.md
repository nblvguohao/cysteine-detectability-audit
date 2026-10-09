# Release notes - v3.3.0 (9 October 2026)

Accompanies the Molecular & Cellular Proteomics submission (manuscript package v3.3, title unchanged: *Five
detectability artifacts decide what site-level comparisons can show in cysteine-modification proteomics*).
Changes from v3.2.0 (commit b6fd491). Cys-Audit is unchanged (0.2.2). No analysis was re-run and no stored
analysis output changed: everything under `cys-audit/`, `protocols/`, `scripts/`, `results/`, `artifacts12/`,
`posthoc_2026-09-30/`, `posthoc_2026-10-02/`, `recoding_audit/`, `analysis/`, `source_data_submitted/` and
`figures_submitted/` is byte-identical to v3.2.0.

## 1. Supplemental material renumbered for v3.3 (`supplemental/`)

- Added `Supplemental_Results.md` (Supplemental Results, sections SR1-SR12) and
  `Supplemental_Figures_legends.md` (legends of Supplemental Figures S1-S7).
- Cross-references updated to the v3.3 numbering: figure numbers (for example Figure 2 of v3.2 is now
  Supplemental Figure S1, Figure 5a is now main Figure 2), Results sections (R1-R7) and Supplemental Results
  sections (SR1-SR12), and reference numbers, in the Supplemental Methods, Supplemental Notes 1, 2, 4-8, 14, 16-18,
  20 and 21, and Supplemental Data 12 (column `where_reported`, explained in its README).
- Added text: in the Supplemental Methods a "Data provenance" section that complements the Data Availability
  statement of the main text; in Supplemental Note 1 the table "Checks not shown in Table 1"; in Supplemental
  Note 14 a paragraph on the tool version and the `SFE-006` random control. Minor wording changes accompany the
  renumbering in some of these files.
- The other files are unchanged in content. Text files are committed with LF line endings, as in every release.

## 2. Figures redrawn for the v3.3 numbering

| v3.3 | content (v3.2 source) |
|---|---|
| Figure 1 | workflow illustration of v3.2 Figure 1a re-used, five steps renumbered, negative-set schematic, v3.2 Figure 1c |
| Figure 2 | v3.2 Figure 5a redrawn (detectability-only AUC against the two negative sets) |
| Figure 3 | cleavage geometry: positional K/R profile, matched enrichment z, distal-band odds ratios of PXD063463 |
| Figure 4 | identified-site share (Supplemental Notes 4 and 7), v3.2 Figure 3a and 3d (search space) |
| Figure 5 | v3.2 Figure 6 redrawn (site-level claims under the matched control) |
| Figure 6 | pLMSNOSite benchmark (Supplemental Note 22) and v3.2 Figure 7c |
| Figure S1, S2, S3, S6 | v3.2 Figures 2, 3, 7 and 4, unchanged |
| Figure S4 | protein-level claims of v3.2 Figure 6 |
| Figure S5 | v3.2 Figure 5b and 5c, placed unchanged |
| Figure S7 | v3.2 Figure 1b redrawn |

## 3. New folders

- `figures_final/`: the v3.3 figures and graphical abstract as submitted (PDF and PNG).
- `source_data_final/`: the Source Data submitted with v3.3 (14 figure files and the five `Source_Data_text_*`
  files, which are identical to those in `source_data_submitted/`). Their provenance fields name tables in
  `source_data_submitted/`, which therefore keeps its name and content.
- `posthoc_2026-10-09_figures/`: the build scripts of the two folders above, which read only files of this
  repository, and the value checks (`qa/`). Re-run in a clean checkout, they reproduce every PNG byte for byte and
  every Source Data file after line-ending normalisation; nine PDF files differ only in the timestamp of the
  embedded font program (see its README).

`figures_submitted/` and `source_data_submitted/` are kept unchanged as the v3.2 set.

## 4. Documentation

`README.md` (contents table and a section "Figures and Source Data") and `CITATION.cff` (version v3.3.0) updated.
`MANIFEST_SHA256.txt` regenerated.

## 5. Edited files and recorded hashes

Some files were edited to remove data and text that cannot be redistributed, and some files were removed; their recorded hashes in the listed provenance files refer to the earlier versions. The provenance, audit and documentation files that record such hashes (in full or abbreviated) are:

- `RELEASE_NOTES_v3.1.2.md`
- `artifacts12/MANIFEST.json`
- `artifacts12/README.md`
- `artifacts12/positional_profile/results/positional_kr_profile_public_audit.json`
- `posthoc_2026-09-30/inputs/repo_results/pxd015307_research_audit.json`
- `posthoc_2026-09-30/results/B_overadjust_sim/provenance.json`
- `posthoc_2026-09-30/results/C_pxd063463_specific/provenance.json`
- `posthoc_2026-09-30/results/D_empirical_background/provenance.json`
- `posthoc_2026-09-30/results/E_bootstrap_coverage/provenance.json`
- `posthoc_2026-09-30/results/F_rice_artifact3/provenance.json`
- `posthoc_2026-09-30/results/G_plmsnosite_rescore/detectability_reproduction.json`
- `posthoc_2026-09-30/results/G_plmsnosite_rescore/provenance.json`
- `posthoc_2026-09-30/results/G_plmsnosite_rescore/withinprot_mean.json`
- `posthoc_2026-09-30/results/H_artifact5_docs/provenance.json`
- `posthoc_2026-09-30/results/I_fig5a_estimand/feature_identity_check.json`
- `posthoc_2026-09-30/results/I_fig5a_estimand/provenance.json`
- `posthoc_2026-09-30/results/J_provenance_species/provenance.json`
- `posthoc_2026-09-30/results/figures_final/Fig5/fig5_rebuild.py`
- `posthoc_2026-09-30/results/figures_final/Fig6/build_log.txt`
- `posthoc_2026-09-30/scripts/H_artifact5_docs/h5_run_all.py`
- `protocols/phase4_cys_audit_preregistration_2026-09-22.json`
- `protocols/phase5b_covariate_control_preregistration_2026-09-22.json`
- `results/public_refit_ablation_lgb_rank_2026-09-24_audit.json`
- `results/pxd015307_research_audit.json`
- `results/search_space_mass_accuracy_audit.json`
- `results/self_audit_human_trained_rice_2026-09-23_audit.json`
- `results/self_audit_public_cohorts_2026-09-20_audit.json`
- `results/self_audit_rice_trained_human_2026-09-24_audit.json`

`MANIFEST_SHA256.txt` was regenerated and lists the current files.
