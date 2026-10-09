# Post hoc analyses (2026-09-30)

Folder names: `posthoc_2026-09-30/` and `posthoc_2026-10-02/` were named `revision_2026-09-30/` and `revision_2026-10-02/` before release v3.2.0. The files inside are unchanged byte for byte, so the SHA-256 values recorded in their `provenance.json` and audit files still hold; paths recorded inside those files and the working-folder constants of their scripts use the former names.

This folder holds the analyses added to the manuscript after the primary analyses were complete: Supplemental Notes 14 to 22, the new Figures 1, 3, 5 and 6, and the amendments to Supplemental Notes 1, 2, 5, 6 and 8 to 13. All of them are post hoc and are labeled as such in the manuscript and in Supplemental Data 12 (rows R14 to R22). The Notes themselves are in `../supplemental/`; the Source Data are in `../source_data_submitted/`.

| Note | subject | outputs | code |
|---|---|---|---|
| 14 | final claim tally; what the matched control changes, by attribute class | `results/final/` | `results/final/check_final_tables.py` |
| 15 | whether the matched background removes genuine chemistry (simulation) | `results/B_overadjust_sim/`, `results/final_notes/Supplemental_Note_15_number_check.csv` | `scripts/B_overadjust_sim/`, `results/final_notes/Supplemental_Note_15_check_numbers.py` |
| 16 | how much cleavage geometry each background removes (PXD063463) | `results/D_empirical_background/` | `scripts/D_empirical_background/` |
| 17 | calibration of the protein-clustered percentile interval | `results/E_bootstrap_coverage/` | `scripts/E_bootstrap_coverage/` |
| 18 | hydroxylamine-dependent sites and identification-matched backgrounds (PXD063463) | `results/C_pxd063463_specific/` | `scripts/C_pxd063463_specific/` |
| 19 | abundance matching, cysteine count and the rice cohort | `results/F_rice_artifact3/` | `scripts/F_rice_artifact3/` |
| 20 | absolute discrimination behind Figure 5a | `results/I_fig5a_estimand/` | `scripts/I_fig5a_estimand/` |
| 21 | how the self-audit labels were assigned; reference level of the top-100 statistic | `results/J_provenance_species/` | `scripts/J_provenance_species/` |
| 22 | the published S-nitrosylation predictor re-scored on its released test set | `results/G_plmsnosite_rescore/` | `scripts/G_plmsnosite_rescore/` |
| 6 | re-search of PXD015307: precursor accuracy, fragment binning, the deposit's sulfide PSMs | `results/H_artifact5_docs/` | `scripts/H_artifact5_docs/` |

Each results folder carries the files its Note cites and a `provenance.json` that records the SHA-256 of the inputs, scripts and outputs of the analysis, the seeds and the software versions. `results/figures_final/` holds the build scripts, checks and outputs of the four rebuilt figures.

**Inputs.** `inputs/` holds the tables derived from the public deposits and copies of the released tables the analyses read. External inputs (UniProt release 2026_03 proteomes, two PRIDE tables of PXD072035 and PXD072089, and the pLMSNOSite repository at commit e9158af) are not redistributed; `DOWNLOADS.md` gives their URLs and SHA-256 values.

**Software.** Python 3.10.11 with numpy 2.2.6, pandas 2.3.3, scipy 1.15.3, scikit-learn 1.7.2, statsmodels 0.14.6, LightGBM 4.7.0 and biopython 1.88; the released pLMSNOSite models were run under TensorFlow 2.15.1.

**Checks.** `python results/final/check_final_tables.py` verifies every verdict, estimate and retained fraction of Supplemental Note 14 against Supplemental Data 2, 4, 5, 6 and 10 and the `SFE-006` Source Data file. `python results/final_notes/Supplemental_Note_15_check_numbers.py` recomputes every number of Supplemental Note 15 from the stored simulation outputs.
