# Display numbering update, 2026-10-26

Use DISPLAY_CROSSWALK_2026-10-26.md for the current four main figures and Fig. S7. The retained documentation below describes AV1 display numbering and original source filenames. Numeric tables are unchanged.

# AV1 figure numbering update (24 October 2026)

This dictionary uses current filenames. Figure2b,c now show the AV1 common-target application; complete older seven-cohort panels are S6. Older figure5a–d filenames still refer to S1, as documented below.

| Current source | Meaning |
|---|---|
| SourceData_Fig2b.csv | Three compound AUROC A/B values, observed-minus-all-other differences and paired component-bootstrap 95% bounds, plus equal-compound mean difference. Blank A/B fields for the mean row are not zero. |
| SourceData_Fig2c.csv | Supported versus positive-bearing discovery proteins and their fraction; training support, not a test exclusion rule. |
| SourceData_FigS5.csv | Same AV1 points and support values in the supplementary display. |
| SourceData_FigS6a–c.csv | Complete prior seven-cohort fixed-scorer and refitted-model display. S6a equals current Fig2a. |
| CXXC_sparse_counts.tsv | Full common-support CXXC four-cell and event-protein counts for TableS14.2; counts independently checked against the archived universes. |

AV1 intervals condition on fitted OOF predictions, the supplied abundance scale, and one technical replicate pair. They are not biological-replicate intervals or chemical-truth validation.

# Figure source-data dictionary

Figure 2a and Figure S6 contain current strict_mc_v1 outputs from the versioned seven-cohort reconstruction; its complete fitted models and exact matrices are supplied in Model_Reconstruction. Other CSVs retain their original analysis scopes. CSVs are presentation tables; the renderer does not change estimates. The only table restructuring in this revision is Fig5d: a reported aggregate of two outcomes is expanded into the two one-claim categories explicitly described in SI S2.

| CSV | Rows / plotted content | Fields and interpretation |
|---|---|---|
| SourceData_Fig1c | 4 counts | `quantity` states the unit: the first three values are **proteins**, the fourth is **site–reference pairs**. `value` is a count; empty CI fields mean no interval was plotted. `source_file`/`filter` identify the original support-summary entry. |
| SourceData_Fig2a | 7 fixed VIS10 differences | `point`, `lo95`, `hi95`: current fixed NEG_A scorer AUROC(P,NEG_B) minus AUROC(P,NEG_A), paired unit-bootstrap 95% interval. `mde80` is a resolution summary, not an equivalence bound. |
| SourceData_FigS6b | 42 background-specific global estimates | Seven cohorts × three feature sets × two backgrounds. `point`, `lo95`, `hi95`: current strict_mc_v1 AUROC and unit-bootstrap 95% interval. Both background models share the union-based held-out folds. |
| SourceData_FigS6c | 42 background-specific macro estimates | Equal-accession mean of within-protein AUROCs over proteins with both classes in each respective background. Unit-bootstrap 95% interval; support can differ across backgrounds. This is distinct from fixed-scorer common-unit macro estimates in SI S3. |
| SourceData_Fig3a | own rule, trypsin rule, paired contrasts | `est`, `lo975`, `hi975`: log2 OR or paired own-minus-trypsin difference with 97.5% bootstrap bounds. `rule_kind`, `contrast`, `arm`, `background`, `band` identify the estimand. The intervals adjust for two bands, not all analyses. |
| SourceData_Fig3b | observed backgrounds and trypsin ladder | Same estimate fields as Fig3a. `n_pos`, `n_bg`, `identified_cys` are different denominators. `mde80_a025` is a resolution measure, not a directly observed effect. |
| SourceData_Fig3c | 14 positional-profile offsets | `offset`: residue distance excluding 0. `fg_frac_KR`, `bg_proteome_frac_KR`, `bg_observed_frac_KR`: descriptive fractions. `n_fg_windows` can vary by offset near sequence ends; `n_sites` is the cohort total. No per-offset inferential CI is supplied. |
| SourceData_Fig3d | raw and theoretically restricted residue comparisons | `est`, `lo95`, `hi95`: presence log2 OR with 95% bounds. `background=theoretically_detectable` is a **restriction**, not propensity matching, and also changes positive eligibility. `n_positive_sites` records this population change. |
| SourceData_Fig4a | null and planted-strength calibration summaries | `PASS`, `WARNING`, `FAIL`, `UNDECIDABLE` are proportions (0–1), converted to percent for plotting. `level` is the interval level. Panel a selects `strength=0`. `flag_excludes_null` and `flag_wilson_lo/hi` describe the null-exclusion frequency; it is distinct from FAIL alone. |
| SourceData_Fig4bcd | ISO/NONISO confirmation cells | `flag_rate` and Wilson bounds are proportions. `theta_pop` is a **legacy field name**: nonzero values are same-batch medians of uncorrected estimates, i.e. empirical references, not independently known population truths. Both null references are zero. `coverage` is coverage relative to that reference. `mean_point_corr` is the matched mean, not an independently proven unbiased estimator. `residual_bias_corr` must not be treated as population bias without the manuscript's qualification. These cells concern the single 95% near-cut statistic, not direct calibration of the deployed combined two-band rule. |
| SourceData_Fig5a | mixed input and diagnostic summary | `quantity` determines the type of `value` (count, proportion, status or text). Only the assignment fraction has a plotted 95% CI. `observed_background_usable=True` records an input-size property, not read-out validity or a PASS verdict. Historical cleavage status is UNDECIDABLE because no FASTA was supplied. The 0.90 convention was selected on this deposit. |
| SourceData_Fig5b | 4 global benchmark estimates | `value`, `ci_low`, `ci_high`: global AUROC and 95% bounds. The label-permuted DIG25 row is a control. |
| SourceData_Fig5c | 2 macro values and 1 paired difference | First two rows: macro within-protein AUROC. Third row: the paired pLMSNOSite−VIS10 difference and its own CI; it must not be obtained by subtracting individual interval endpoints. Unresolved does not establish equivalence. |
| SourceData_Fig5d | 5 claim categories | Counts are Survives 4, Attenuated 2, Undecidable 10, Null broken 1, Transfer contradiction 1, totaling 18. The last category is a transfer-cohort baseline contradiction, not reversal on the original cohort. Ten protein-level claims are outside this site-level denominator; their historical estimates remain in archived re-test tables. See Claim_Unit_Interpretation_2026-10-24.md. `source_file` cites SI S2; the separate claim ledger provides finer provenance when available. |

## Shared fields and display annotations

- `panel` is the manuscript panel letter. `source_file`, `filter`, `overlay_source`, `overlay_filter` and `note` are provenance metadata, not numeric observations.
- NEG_A comprises other cysteines in proteins bearing assigned sites. NEG_B comprises experiment-observed cysteines without an assignment. NEG_B may include proteins with no assigned site; it is not globally required to be a subset of NEG_A.
Source `filled` fields are historical display choices. Current Fig2 uses † for different positive/reference selections, ‡ for the qPerS-SID mixed Beads/Elution reference, and § for QTRP parallel SH/SSH profiles. All seven NEG_B memberships were traced from original workbooks; none is certified as a chemical-negative class.
- Historical cohort aliases: `aBIOTECH (Ath)` → PAT-switch (Ath); `Nat Commun (Ath)` → FAT-switch (Ath). FAT-switch has no NEG_B and is excluded from current Fig2. Its historical fixed-difference source panel was Fig2c; current fixed differences are Fig2a.
- A blank numeric field means missing or inapplicable, not zero. CI limits use the specific confidence level encoded in their column name.


Final layout note (20 October 2026): former main Figure 5 is Supplementary Figure S1, exported in Supplementary_Figures/FigS1.*. Historical SourceData_Fig5a–d.csv filenames and render audit labels are retained for traceability. In the v9 layout, main figures are Fig1–5; the new main Fig5 is the external ACR/IAM case, with source files listed below. The current Fig2 legend distinguishes source-defined populations after original-workbook membership tracing.


## v9 additions (21 October 2026)

| Current figure | Source data | Interpretation |
|---|---|---|
| Fig5a | Grounding_2026-10-21/SourceData_Fig5_support.csv | Positive-bearing versus IAM-reference-supported protein counts. |
| Fig5b,c | Grounding_2026-10-21/SourceData_Fig5_effects.csv | Common-support motif estimates and conditional protein-bootstrap percentiles. |
| FigS2 | Grounding_2026-10-21/SourceData_FigS2.csv | Actual assay subsets, fixed-score contrasts and unit-mass changes. |
| FigS3 | Grounding_2026-10-21/SourceData_FigS3.csv | All-positive-protein external estimates, before common-support restriction. |
| FigS4 | Grounding_2026-10-21/SourceData_FigS4.csv | Mean/std historical sequence-unit replay; QTRP S2 budget360, others400. |

Historical SourceData_Fig5a–d.csv belong to FigS1, not the new main Figure5. All new outcomes, failures and sparse-cell sensitivities are retained in Grounding_2026-10-21 and the accompanying 06 analysis archive.
