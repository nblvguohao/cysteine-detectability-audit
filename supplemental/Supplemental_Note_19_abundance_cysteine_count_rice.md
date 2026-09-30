# Supplemental Note 19. What abundance matching exposes in the rice protein-level contrast, and how the rice cohort was built (post hoc)

**Status.** Added in revision, post hoc, in response to review. The original Artifact 3 analysis, the abundance ratio and Mann–Whitney tests of protein length and cysteine count before and within abundance deciles, followed an analysis plan written before it ran; its numbers are reproduced here, not re-specified. That plan provided for two readings, and the outcome was the second (abundance does not account for the protein-level signal); reading the outcome as an artifact of matching is a post hoc interpretation, as stated in the main text. The rank-biserial effect sizes and the one-to-one abundance matching come from a follow-up written after the original outcome was known. Everything else in this Note, including the reconstruction of the rice site list (Part B), was designed after review. All values are read from stored outputs (Data and code, below).

## A. What abundance matching exposes

**Units and sets.** Protein groups come from the total-proteome report of PXD072035 (Spectronaut; 7,825 groups, 7,692 with a usable intensity-based absolute quantification, iBAQ) and peptides from the persulfidation peptide table of PXD072089 (MaxQuant). Each group is represented by its first-listed accession, whose iBAQ, the median over runs, gives the group's abundance; the report lists accessions alphabetically in all 7,825 rows. The original positive set is the 1,197 groups that contain a leading razor protein of the peptide table. The comparison groups are the other 6,495, none of whose accessions is a leading razor protein of the table; 803 of them contain an accession that the table lists among a peptide's proteins.

Source: `revision_2026-09-30/results/F_rice_artifact3/verify_r4/build_r4.json`.

**Reproduction.** Every number of the original analysis that does not depend on protein sequences reproduces exactly from the deposited tables, to the four decimals stored (Table S19.1): 1,197 against 6,495 groups, a median iBAQ ratio of 10.4 (95% CI 8.95–12.3), a raw-intensity ratio of 4.88 and a shuffled-label control of 0.97. The sequence file of the original run was not retained. We took each first-listed accession's sequence from UniProtKB release 2026_03 (1,698 groups from the rice reference proteome, 515 from other active entries) or, for the 5,479 accessions not returned as active entries of that release, from the latest UniParc sequence carrying the accession (identical to UniSave in 25 of 25 spot checks). This covers all 7,692 groups, against 6,363 in the original run, and changes the six rank-biserial values of the follow-up by at most 0.022, none of them in sign.

**Table S19.1. Original Artifact 3 values and their reproduction.** Reproduced with the original functions and seed (bootstrap: 5,000 resamples over protein groups). The first six rows do not depend on protein sequences.

| quantity | original | reproduced |
|---|---|---|
| protein groups: in the report / without a usable iBAQ / analyzed | 7,825 / 133 / 7,692 | 7,825 / 133 / 7,692 |
| original positive / comparison groups | 1,197 / 6,495 | 1,197 / 6,495 |
| leading razor proteins in the peptide table | 2,073 | 2,073 |
| median iBAQ ratio [95% CI] | 10.3633 [8.9464, 12.2902] | 10.3633 [8.9464, 12.2902] |
| raw-intensity ratio | 4.8789 | 4.8789 |
| shuffled-label control [95% CI] | 0.9717 [0.8054, 1.2248] | 0.9717 [0.8054, 1.2248] |
| groups with a sequence (positive / comparison) | 6,363 (1,062 / 5,301) | 7,692 (1,197 / 6,495) |
| length, rank-biserial r: unmatched / within iBAQ deciles / one-to-one matched on iBAQ | −0.009 / +0.201 / +0.170 | −0.026 / +0.199 / +0.150 |
| cysteine count, rank-biserial r: same designs | +0.104 / +0.290 / +0.256 | +0.082 / +0.286 / +0.244 |
| length, Mann–Whitney z: unmatched / stratified by decile / one-to-one matched | −0.48 / 8.56 / 6.79 | −1.43 / 8.94 / 6.35 |
| cysteine count, Mann–Whitney z: same designs | 5.35 / 12.77 / 10.23 | 4.53 / 13.28 / 10.36 |

Source: `revision_2026-09-30/results/F_rice_artifact3/s02_reproduction.csv` and `.../fetched/fetch_summary.json` (sequence sources and spot checks). With the rice reference proteome alone (6,266 groups), the rank-biserial values differ from the originals by at most 0.016.

**Why matching enlarges the contrast.** Abundant groups are shorter and poorer in cysteine (Spearman correlation with log10 iBAQ, −0.28 for length and −0.24 for cysteine count), and the positives are concentrated among abundant groups. Without matching, positives are compared with longer, cysteine-richer groups of lower abundance, which masks the contrast; matching on abundance exposes it.

Source: `revision_2026-09-30/results/F_rice_artifact3/s04_stdout.txt`.

**Site-level set and cysteine counts.** Of the 1,197 original positives, 750 contain the leading razor protein of at least one identified cysteine-containing peptide, decoys and contaminants removed; we call them the site-level set. The other 447 are represented in the table only by peptides without a cysteine. At the median, the site-level set is 7.2-fold more abundant than the comparison groups (95% CI 5.9–8.7) and the 447 are 17.8-fold (14.7–22.8). A cysteine is detectable when at least one fully tryptic peptide of 7–30 residues with up to two missed cleavages covers it (the detectability model of Experimental Procedures) and undetectable otherwise. 81% of cysteines are detectable, the detectable and total counts of a group correlate at Spearman 0.94, and 6,990 groups have at least one detectable cysteine.

Source: `revision_2026-09-30/results/F_rice_artifact3/s04c_label_ratios.csv` and `.../s10_summary.json`.

**The exposed contrast is one of cysteine count (Table S19.2).** At fixed abundance and length, the odds that a group is in the original positive set rise by 1.08 (95% CI 1.06–1.09) per detectable and by 1.07 (1.03–1.11) per undetectable cysteine (P = 0.80 for equality), and by 1.08 (1.06–1.09) per cysteine when the total count enters alone; in the site-level set the values are 1.11, 1.08 (P = 0.22) and 1.10. The two per-cysteine odds ratios do not differ under any detectability window tested (6–35 residues with up to two missed cleavages; 7–40 and 5–50 residues with up to three; tryptic cleavage also allowed before proline, 7–30 residues; P ≥ 0.18 in both sets). Models with the total count fit better than models with the detectable count in every model form, by 8.9–19.8 AIC units in the original set and 19.5–27.4 in the site-level set. Given the total count, the detectable count adds nothing in the original set (likelihood-ratio P ≥ 0.09 in every form tested, with or without the length term), whereas the total count adds to the detectable count in every form (P < 0.001). In the site-level set, the detectable count adds only the requirement that a positive have a detectable cysteine: P = 0.005 over all groups, 0.15 once that requirement enters both models as an indicator, and 0.10–0.32 among the 6,990 groups that meet it. The AIC separation rests on the primary window, under which about a fifth of cysteines are undetectable: under the three broader windows (87–96% of cysteines detectable; Spearman 0.96–0.98 between the counts), the difference lies between −1.9 and +3.5 AIC units in both sets.

**Length adds nothing at fixed cysteine count.** The odds ratio per doubling of length at fixed total count is 1.07 (0.97–1.17; P = 0.18) in the original set and 0.92 (0.82–1.03; P = 0.13) in the site-level set, with P ≥ 0.10 in every model form, the linear-count form included. These models measure abundance as iBAQ; at fixed raw intensity (which itself rises with length) and fixed cysteine counts, longer proteins are less often positive (0.81 per doubling, 0.74–0.89, original set).

**Table S19.2. Detectable versus total cysteine count.** Regressions of membership in the positive set on log10 iBAQ and log2 length. Odds ratios per cysteine: detectable and undetectable counts as separate linear terms, all 7,692 groups. Length: logistic with log2(1 + total count), all groups. ΔAIC: AIC of the model with log2(1 + detectable count) minus that of the model with log2(1 + total count); positive values favor the total count. P: likelihood-ratio tests, in the order of the first five ΔAIC rows (the categorical form has no nested test).

| statistic | original positive set (1,197) | site-level set (750) |
|---|---|---|
| odds ratio per detectable cysteine [95% CI] | 1.08 [1.06, 1.09] | 1.11 [1.09, 1.13] |
| odds ratio per undetectable cysteine | 1.07 [1.03, 1.11] | 1.08 [1.04, 1.12] |
| P, equal odds ratios | 0.80 | 0.22 |
| odds ratio per cysteine, total count alone | 1.08 [1.06, 1.09] | 1.10 [1.09, 1.12] |
| odds ratio per doubling of length, at fixed total count | 1.07 [0.97, 1.17] | 0.92 [0.82, 1.03] |
| ΔAIC: logistic, all 7,692 groups | 10.6 | 19.5 |
| ΔAIC: same, with an indicator for at least one detectable cysteine | 8.9 | 27.4 |
| ΔAIC: logistic, 6,990 groups with a detectable cysteine | 19.0 | 22.9 |
| ΔAIC: complementary log–log with a free exponent on the count, 6,990 groups | 17.0 | 23.4 |
| ΔAIC: spline in abundance, quadratic count, 6,990 groups | 19.8 | 24.1 |
| ΔAIC: spline in abundance, categorical count capped at 20, 6,990 groups | 12.5 | 22.3 |
| P, detectable count given total count | 0.26, 0.11, 0.95, 0.82, 0.48 | 0.005, 0.15, 0.10, 0.11, 0.32 |
| P, total count given detectable count (largest of the five) | 7.0 × 10⁻⁴ | 1.9 × 10⁻⁶ |

Source: `revision_2026-09-30/results/F_rice_artifact3/s10_linear_symmetric.csv`, `.../s10_comparisons.csv` (also the forms without length), `.../s13_indicator_test.csv`, `.../s10_windows.csv` (detectability windows) and `.../s08_sensitivity_undet.csv` (raw intensity).

**Matched designs.** After one-to-one matching on abundance within strata of exact detectable count, original positives still carry more cysteines, namely undetectable ones (rank-biserial r = +0.066, 95% CI 0.043–0.12; site-level set +0.097, 0.054–0.15), and are no longer than their matches (r = −0.004, −0.036 to 0.030). After matching within strata of exact total count, original positives carry no more detectable cysteines (r = +0.010, −0.012 to 0.023; 6.20 against 6.06 per group, Wilcoxon P = 0.72) and are no longer (r = −0.016, −0.045 to 0.024); site-level positives carry slightly more detectable cysteines (+0.028, 0.003–0.046) and are slightly shorter (−0.071, −0.10 to −0.018). Pairs share the exact stratum in 95–97% of cases; intervals come from 400 protein-group bootstrap resamples with the matching redone.

Source: `revision_2026-09-30/results/F_rice_artifact3/s08_matched_pairs.csv` and `.../s10_mirror_matched.csv`.

**Calibrated null models (Table S19.3).** Four null models, in which a group's chance of being positive depends only on its abundance and one cysteine count, were fitted to the observed labels of the 6,990 groups with a detectable cysteine (1,113 original and 746 site-level positives): with the total or the detectable count, either flexible (logistic, natural spline in log10 iBAQ with 5 degrees of freedom, count as a categorical term capped at 20) or rigid (complementary log–log with a free exponent on the count). Each was simulated 5,000 times, and 14 contrasts were computed on every simulated label set. Two-sided P is twice the smaller tail proportion (smallest attainable 0.0004), and P between 0.04 and 0.06 is called the edge of the range; this band, and the use of all 5,000 sets, were fixed before the last 4,000 sets were drawn but after earlier simulations had placed one contrast near the edge.

For the original set, the flexible total-count null places 13 of the 14 contrasts well inside its simulated 95% range (P ≥ 0.46), including the length and cysteine-count contrasts within abundance deciles (P = 0.92 and 0.98) and every contrast computed at a fixed number of detectable cysteines (rows 5–9). The fourteenth, the cysteine-count contrast within abundance decile × detectable-count tertile, lies at the edge (P = 0.049; an independent re-simulation of 10,000 sets during verification placed it just inside, P = 0.057). The flexible detectable-count null places 8 of the 14 outside its range (each P < 0.01), each of them a contrast between detectable and total count (rows 4, 6–8, 10 and 12–14). Both flexible nulls reproduce the within-decile contrasts, which therefore do not discriminate between the two counts. In the site-level set, the flexible total-count null leaves all nine cysteine contrasts inside its range and misses two length contrasts, in which positives are slightly shorter than it predicts (rows 5 and 11: r = −0.043 and −0.077; P = 0.022 and 0.013); the flexible detectable-count null misses 9. The rigid nulls fit less well with either count. In the original set they place 3 (total count) and 11 (detectable count) contrasts outside the range, the rigid total-count null missing both within-decile contrasts (P = 0.045 and 0.024) and row 4 (P = 0.0012); in the site-level set they place 2 (both length) and 9 outside. Three contrasts lie at the edge, all in the original set: row 4 under the flexible total-count null (P = 0.049), row 1 under the rigid total-count null (0.045) and row 3 under the rigid detectable-count null (0.041).

**Table S19.3. Original positive set: observed contrasts and their simulated 95% ranges under the two flexible nulls.** 6,990 groups with at least one detectable cysteine (1,113 positives); 5,000 simulated label sets per null; two-sided P in parentheses. r, rank-biserial r, positive minus comparison groups; matching is one-to-one on log10 iBAQ within strata of the exact count named. Rows 8 and 9 are logistic coefficients given log10 iBAQ, log2(1 + detectable count), log2 length and the undetectable count; row 12 is the ΔAIC of Table S19.2 on these groups; row 13 comes from the model with both counts; row 14 compares the per-cysteine coefficients of detectable and undetectable cysteines.

| row | contrast | observed | total-count null | detectable-count null |
|---|---|---|---|---|
| 1 | length, r within abundance deciles | +0.17 | 0.12 to 0.22 (0.92) | 0.10 to 0.20 (0.38) |
| 2 | cysteine count, r within abundance deciles | +0.26 | 0.21 to 0.31 (0.98) | 0.19 to 0.28 (0.20) |
| 3 | length, r within abundance decile × detectable-count tertile | +0.068 | 0.016 to 0.11 (0.88) | −0.004 to 0.092 (0.33) |
| 4 | cysteine count, r within the same strata | +0.186 | 0.091 to 0.185 (0.049) | 0.047 to 0.14 (0.0008) |
| 5 | length, r matched at exact detectable count | −0.007 | −0.036 to 0.039 (0.70) | −0.052 to 0.023 (0.65) |
| 6 | cysteine count, r, same matching | +0.023 | 0.011 to 0.044 (0.55) | −0.016 to 0.017 (0.008) |
| 7 | undetectable cysteines, r, same matching | +0.070 | 0.021 to 0.10 (0.70) | −0.041 to 0.041 (0.0008) |
| 8 | undetectable cysteines, log odds per cysteine | +0.087 | 0.038 to 0.11 (0.46) | −0.033 to 0.040 (0.0004) |
| 9 | length, log odds per doubling, same model | +0.036 | −0.083 to 0.11 (0.69) | −0.082 to 0.11 (0.62) |
| 10 | detectable cysteines, r matched at exact total count | +0.0004 | −0.020 to 0.015 (0.76) | 0.006 to 0.040 (0.009) |
| 11 | length, r, same matching | −0.025 | −0.054 to 0.021 (0.67) | −0.049 to 0.024 (0.53) |
| 12 | ΔAIC, detectable minus total count | +19.0 | −1.8 to 32.5 (0.63) | −33.6 to −2.2 (0.0004) |
| 13 | log2(1 + detectable count), coefficient given log2(1 + total count) | −0.008 | −0.15 to 0.27 (0.58) | 0.26 to 0.72 (0.0004) |
| 14 | per-cysteine log odds, detectable minus undetectable | −0.007 | −0.044 to 0.042 (0.79) | 0.038 to 0.13 (0.0004) |

Source: `revision_2026-09-30/results/F_rice_artifact3/s12_contrast_table_pooled.csv` (all four nulls and both sets) and `.../s12_null_pooled_summary.csv`; the independent re-simulation is in `.../verify_r4/comparison_r4.csv`.

**Reading.** These analyses identify the property that the matched contrast tracks, the number of cysteines a protein carries rather than their theoretical detectability or the protein's length, but they cannot by themselves separate more chances of capture, one per cysteine, from an association of persulfidation with cysteine-rich proteins.

## B. How the rice self-audit cohort was built

The rice deposit (PXD072089) holds one processed file, a MaxQuant peptide table, and no cysteine site table; the table's only site-level modification column is for methionine oxidation. A site is therefore a cysteine of a peptide identified in the enriched arm, placed on that peptide's leading razor protein at the table's start position. The list used in the self-audit was built in this project on 12 September 2026 with a script that is not available; the chain in Table S19.4 reproduces it site by site from the public table.

**Table S19.4. Reconstruction of the rice site list from the deposited peptide table.**

| step | sites | proteins | cysteines in these proteins |
|---|---|---|---|
| every cysteine of every identified cysteine-containing peptide, on its leading razor protein; decoys and contaminants removed | 1,753 | 1,164 | 10,662 |
| single-cysteine peptides only | 1,377 | 1,077 | 9,733 |
| isoform accessions mapped to the canonical entry (23 peptides, each located uniquely) | 1,377 | 1,077 | 9,747 |
| proteins with at least two cysteines | 1,340 | 1,040 | 9,710 |
| UniProtKB entry active in release 2026_03 | 817 | 640 | 5,598 |

Source: `revision_2026-09-30/results/F_rice_artifact3/s11_rice_chain.csv` and `.../s11_rice_chain.json`; the deposit's file listing is in the code repository (`results/pride_coincidence_candidates_2026-09-19.csv`). The last row is identical, site by site, to the self-audit list.

The last step was not a stated filter. All 400 proteins it removes (523 sites; 4,112 cysteines) have deleted UniProtKB entries (400 of 400 entry records: deleted, "not part of a reference proteome"; for 20 sampled entries, the last release carrying them was 2026_01, of 28 January 2026), and a UniProtKB query for them with the sequence field returned none; we recovered their sequences from UniParc. The project's records indicate that the list-building script took the sequences of proteins outside its reference proteome from a file of sequences retrieved one accession at a time, so the most likely explanation is that no sequence was obtained for these deleted entries; because the script is not available, this is an inference.

Source: `revision_2026-09-30/results/F_rice_artifact3/fetched/inactive_reasons_400.csv`, `.../fetched/unisave_last_release_sample.csv` and `.../fetched/fetched_sequences.tsv`.

When the zero-shot ranker scores the rice cohort, 594 of the 5,598 cysteines (80 of the 817 sites), in 59 proteins absent from the reference-proteome file used, are excluded, leaving 5,004 cysteines, 737 sites and 581 proteins. This step is taken from the repository's audit of that run, as recorded in `revision_2026-09-30/results/F_rice_artifact3/s11_rice_chain.json`, because the reference-proteome file was not retained. The within-cohort refits were trained on all 5,598 cysteines (code repository, `results/public_refit_inputs_2026-09-21_audit.json`). The rice self-audit and, according to the project's records, the rice analysis of Artifact 2 use the 817-site list; neither was re-run on the complete 1,340-site list.

The source publication (Lin et al., 2026) reports 1,691 persulfidated sites. Of 72 exploratory counting rules applied to the peptide table, the closest is every cysteine of the identified cysteine-containing peptides that are unique to one protein group (on the leading razor protein, decoys and contaminants removed), which gives 1,692 sites. The authors' counting rule is not deposited.

Source: `revision_2026-09-30/results/F_rice_artifact3/s05_rule_grid.csv`.

**Limits.** Everything beyond the reproduction is post hoc, and the rank-biserial values of the original follow-up can be recomputed only approximately, because its sequence file was not retained. The distinction between detectable and total count rests on the fifth of cysteines that the primary window calls undetectable, and under broader windows the two counts nearly coincide. Bootstrap intervals resample protein groups and are not homology-clustered. The sequences of 5,479 first-listed accessions come from UniParc; 25 of 25 spot checks against UniSave were identical, but untested accessions could differ.

**Data and code.** All values are read from stored outputs under `revision_2026-09-30/results/F_rice_artifact3/`: `s02_reproduction.csv` (sha256 517720dd…), `s04_stdout.txt` (9fb203ad…), `s04c_label_ratios.csv` (4c93f26d…), `s08_matched_pairs.csv` (0773e8fe…), `s08_sensitivity_undet.csv` (06ba8541…), `s10_summary.json` (be18d709…), `s10_linear_symmetric.csv` (fea8fea8…), `s10_comparisons.csv` (da23df26…), `s10_windows.csv` (79fd8b05…), `s10_mirror_matched.csv` (c2f42fee…), `s12_contrast_table_pooled.csv` (8182f378…), `s12_null_pooled_summary.csv` (03da4508…), `s13_indicator_test.csv` (d2d27538…), `s11_rice_chain.csv` (32391f76…), `s11_rice_chain.json` (24b1cefc…), `s05_rule_grid.csv` (44b8bab7…), `fetched/fetch_summary.json` (d03ada98…), `fetched/inactive_reasons_400.csv` (a4fc9954…), `fetched/unisave_last_release_sample.csv` (c6a1914e…), `fetched/fetched_sequences.tsv` (bc5893b2…), and, from the independent verification, `verify_r4/build_r4.json` (e3d4b36c…) and `verify_r4/comparison_r4.csv` (68d0ca4f…). The code-repository files cited are `results/pride_coincidence_candidates_2026-09-19.csv` (4f135b86…), `results/self_audit_public_cohorts_2026-09-20_audit.json` (976b11fb…; the 594-cysteine exclusion) and `results/public_refit_inputs_2026-09-21_audit.json` (aa4186e8…). The scripts are in `revision_2026-09-30/scripts/F_rice_artifact3/`; `provenance.json` in the results folder records the sha256 of inputs, scripts and outputs, the seeds and the environment.
