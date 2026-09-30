# Supplemental Note 16. How much tryptic cleavage geometry each background removes in the PXD063463 trypsin arm (post hoc)

**Status.** Added in revision, post hoc, in response to review. The analysis was designed after the review. The backgrounds that combine the sequence rule with the unenriched arm, the variants reported with them and the Holm family were defined after earlier results of this analysis had been seen. All values are read from stored outputs (Data and code, below).

**Why this arm.** Artifact 1 corrects cleavage geometry with a sequence rule for theoretical detectability. PXD063463 is the deposit in which the empirical detectability check was run on an unenriched arm of the same material (Experimental Procedures), and its trypsin arm also carries an observed background, the identified cysteines without a site (Fig. 2b; Supplemental Note 4). In this arm the Fig. 2 statistic can therefore be computed against backgrounds that correct for detectability in different ways.

**Data.** The trypsin arm treated with hydroxylamine; a site is a cysteine carrying the carbamidomethyl readout of the exchange step with a localization probability of at least 0.75 (Supplemental Note 4). Sequences are from UniProt release 2026_03 (mouse reference proteome, 21,860 entries) rather than the database used for Fig. 2; 18 of the 2,009 identified proteins, carrying 35 identified cysteines and 5 sites, are absent from that release. The analysis set is 1,991 proteins with 29,502 cysteines, of which 5,293 were identified in the arm and 1,354 are sites. A cysteine counts as identified in the unenriched (global) arm when it lies in a candidate peptide of the empirical detectability check that this arm identified, matched on protein, start position and sequence.

**Reproduction.** Cys-Audit v0.2.2, whose cleavage check and statistics, protease and input modules are identical to those that produced the stored audits, was rerun on the stored input table with the settings and seed of the stored audits and the 2026_03 sequences. Point estimates agree with the stored Fig. 2 values within 0.0096 and interval bounds within 0.0115 (Table S16.1); the differences follow from the change of sequence release.

**Table S16.1. Trypsin arm, trypsin rule: stored and reproduced values.** Log2 odds ratio [97.5% interval, the Bonferroni level of the stored audits].

| background, band | stored (Fig. 2) | reproduced (UniProt 2026_03) |
|---|---|---|
| proteome, distal 6–12 | 0.9730 [0.7115, 1.2756] | 0.9801 [0.7190, 1.2662] |
| proteome, proximal 1–3 | −0.1086 [−0.2853, 0.0688] | −0.1101 [−0.2902, 0.0678] |
| observed, distal 6–12 | 0.0608 [−0.2353, 0.3808] | 0.0705 [−0.2336, 0.3867] |
| observed, proximal 1–3 | −0.6949 [−0.8951, −0.4948] | −0.7043 [−0.9042, −0.5064] |
| sites / proteome background / observed background | 1,359 / 28,543 / 3,969 | 1,354 / 28,148 / 3,939 |

Source: `revision_2026-09-30/results/D_empirical_background/cli_reproduction.json`.

**Backgrounds and statistic.** The statistic is that of Fig. 2: the Haldane log2 odds ratio of carrying a lysine or arginine not followed by proline 6 to 12 residues from the cysteine on either side (the distal band), sites against background, here with 95% protein-clustered percentile intervals from 5,000 bootstrap replicates. Every background is drawn from the other cysteines of the identified proteins, and, as in Artifact 1, restrictions apply to the background only:
- proteome: all of them (Fig. 2a);
- theoretical: those in at least one fully tryptic peptide of 7 to 30 residues with at most two missed cleavages, the rule used throughout the paper;
- combined: those that meet the rule and were also identified in the unenriched arm;
- observed: those identified in the enriched arm itself without a site (Fig. 2b).

Differences between backgrounds were computed in one set of 50,000 paired protein resamples over the 1,991 proteins. Their two-sided bootstrap P values are Holm-adjusted across the ten between-background differences examined in this analysis.

**Table S16.2. Distal band, all 1,354 sites.** The share of sites with a distal K or R is 0.866.

| background | background cysteines | background share with a distal K or R | K or R: log2 odds ratio [95% CI] | share of the proteome estimate removed [95% CI] | D or E in the same band: log2 odds ratio [95% CI] |
|---|---|---|---|---|---|
| proteome | 28,148 | 0.765 | 0.980 [0.754, 1.229] | — | −0.071 [−0.287, 0.152] |
| theoretical | 22,490 | 0.839 | 0.302 [0.069, 0.566] | 0.692 [0.546, 0.916] | −0.133 [−0.345, 0.093] |
| combined | 3,618 | 0.869 | −0.041 [−0.329, 0.239] | 1.042 [0.796, 1.412] | −0.336 [−0.582, −0.076] |
| observed | 3,939 | 0.860 | 0.070 [−0.191, 0.346] | 0.928 [0.708, 1.255] | −0.335 [−0.591, −0.087] |

Log2 odds ratio intervals resample the proteins of each comparison (5,000 replicates); shares removed and their intervals come from the 50,000 paired resamples. The share removed is 1 minus the ratio of the estimate to the proteome estimate; above 1, the estimate fell below zero. It refers to this pooled log2 odds ratio and is not comparable with the reductions of the within-protein z statistic reported for Artifact 1. Source: `revision_2026-09-30/results/D_empirical_background/r2_designs.csv` and `.../r2_paired_contrasts.csv`.

**Table S16.3. Differences between backgrounds, distal band (50,000 paired protein resamples).**

| difference | estimate [95% CI] | Holm-adjusted P |
|---|---|---|
| theoretical − observed | 0.231 [0.108, 0.362] | 0.00064 |
| combined − theoretical | −0.343 [−0.488, −0.208] | 0.00020 |
| combined − observed | −0.112 [−0.266, 0.039] | 0.748 |
| D or E: combined − theoretical | −0.203 [−0.330, −0.079] | not in the family |
| K or R shift minus D or E shift, combined − theoretical | −0.140 [−0.324, 0.040] | not in the family |

Source: `revision_2026-09-30/results/D_empirical_background/r2_paired_contrasts.csv`, `.../r2_holm_round2.csv` (the ten-difference family) and `.../r2_kr_specific_part.csv`.

**Reading.**
- *The theoretical control is a partial correction*, as the main text states. Restricting the background to theoretically detectable cysteines lowered the distal estimate from 0.98 [0.75, 1.23] to 0.30 [0.07, 0.57], but left it 0.23 [0.11, 0.36] above the estimate against the observed background, 0.07 [−0.19, 0.35] (Holm-adjusted P = 0.0006).
- *Restricting the background by both the rule and the unenriched arm reached the level of the observed background.* Against the combined background the estimate was −0.04 [−0.33, 0.24], 0.34 [0.21, 0.49] below the theoretical restriction (Holm-adjusted P = 0.0002) and not resolved from the observed background (difference −0.11 [−0.27, 0.04]; Holm-adjusted P = 0.75). Sites and combined background carry a distal K or R equally often (0.866 and 0.869). The estimate stayed near zero when the sites were restricted in the same way (913 sites; −0.01 [−0.31, 0.33]), and within the 766 proteins that carry a site the combined background was level with the observed one (difference −0.02 [−0.24, 0.20]), whereas the theoretical restriction stayed 0.30 [0.10, 0.51] above it. Reweighting the whole background to the sites' distribution over the same two criteria, instead of restricting it, left more (0.24 [−0.01, 0.51]).
- *Qualification.* Relative to the rule alone, the combined background also shifted D or E in the same band, residues after which trypsin does not cut, by −0.20 [−0.33, −0.08], and the part of the K or R shift beyond that shift, −0.14 [−0.32, 0.04], is not resolved; the additional reduction is therefore not shown to be specific to cleavage residues. D or E is more frequent around cysteines that the unenriched arm identified (log2 odds ratio 0.32 [0.19, 0.45]), and the combined background reproduces the observed background's D or E contrast (−0.34 against −0.33), as a detectability background would. The unenriched arm, moreover, was searched without the blocking readout: carbamidomethylation of cysteine was a fixed modification and N-ethylmaleimide was not searched (MaxQuant 2.4.14.0, Trypsin/P), whereas the enriched arms searched both as variable modifications. As searched, it supplies a detectability background, not an observed negative class of detected-but-unmodified cysteines.

Source: `revision_2026-09-30/results/D_empirical_background/r2_designs.csv` (combined background with the sites restricted, within the 766 site-carrying proteins, and reweighted), `.../r2_paired_contrasts.csv`, `.../r1_comparison_feature.csv` (D or E against identification in the unenriched arm) and `.../r2_global_arm_search_design.json` (search parameters of the arms, read from the deposit's MaxQuant parameter files).

**Limits.** One S-palmitoylation cohort, one protease arm and one unenriched arm of the same material; all comparisons reuse the same sites, and the shares removed are specific to this arm and statistic. The unenriched-arm identifications come from an analysis designed for another purpose and are reused here post hoc. Whether the unenriched aliquot was itself exposed to N-ethylmaleimide cannot be read from its search parameters.

**Data and code.** All values are read from stored outputs under `revision_2026-09-30/results/D_empirical_background/`: `cli_reproduction.json` (sha256 30c81745…), `prepare_checks.json` (cd6e57fd…; analysis set), `r2_designs.csv` (c7b2c2a9…), `r2_paired_contrasts.csv` (2c03b947…), `r2_holm_round2.csv` (913475a2…), `r2_kr_specific_part.csv` (efd181ed…), `r1_comparison_feature.csv` (88712508…) and `r2_global_arm_search_design.json` (19fb9e8f…). These files also hold the other backgrounds examined in this analysis. The scripts are in `revision_2026-09-30/scripts/D_empirical_background/` (`00_reproduce_cli.py`, `01_prepare.py`, `03_backgrounds.py`, `09_revision_after_verification.py`, `10_round2_combined_designs.py`); `provenance.json` in the results folder records the sha256 of inputs, scripts and outputs, the seeds and the environment.
