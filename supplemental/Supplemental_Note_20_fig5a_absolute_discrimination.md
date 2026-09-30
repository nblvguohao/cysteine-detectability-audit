# Supplemental Note 20. Absolute discrimination behind Figure 5a (post hoc)

**Status.** Added in revision, post hoc, in response to review. Every value is read from, or computed arithmetically from, the stored outputs of the analysis behind Figure 5a (`results/ptm_detectability_share.csv`); no model was refitted.

**Models.** For each cohort and negative set, one fixed HistGradientBoosting configuration (no tuning, reweighting or calibration) was fitted with five-fold cross-validation grouped by homology component where available (human) or protein (Arabidopsis), on three feature sets computed under the tryptic rule:
- VIS10: pep_len, pep_log_len, pep_mass, pep_gravy, and the flags pep_detectable_length, pep_detectable_mass, pep_detectable_both, pep_mc1_detectable, pep_mc2_detectable, pep_detectable_any_missed_cleavage;
- CPL15, the columns DIG25 adds: pep_offset_in_peptide, pep_relative_offset, pep_dist_to_peptide_n_term, pep_dist_to_peptide_c_term, pep_mc1/mc2_left/right/min_len (six), pep_is_protein_n/c_terminal, pep_n_peptides_in_protein, pep_cys_count, pep_cleavage_sites_within_20 (K or R within ±20 residues);
- DIG25 = VIS10 + CPL15.

**AUCs.** Global AUC is over all sites. Within-protein AUC is pairwise concordance pooled over the proteins or components that carry both classes, reported when at least ten do. Intervals are 95% cluster-bootstrap percentiles (5,000 replicates, out-of-fold scores reweighted, no refit).

**The ratio.** Figure 5a plots the absolute AUCs. The ratio discussed in Results, (AUC_VIS10 − 0.5)/(AUC_DIG25 − 0.5) on the global AUC, is tabulated here: it is the fraction of the digest model's above-chance discrimination that its visibility subset reaches. Because two of the columns DIG25 adds fix peptide length exactly, it is not a partition of discrimination between visibility and chemistry; neither model contains chemical or structural features.

**Table S20.1. Absolute out-of-fold AUCs, point [95% CI].** Units are proteins or homology components; the number in brackets carries both classes. CPL15 AUCs are in the Source Data.

| Cohort | Neg. set | sites / negatives | units (both classes) | VIS10 global | DIG25 global | VIS10 within protein | DIG25 within protein | ratio |
|---|---|---|---|---|---|---|---|---|
| QTRP S1 pH5 | NEG_A | 277 / 2,056 | 201 (195) | 0.855 [0.835, 0.874] | 0.928 [0.914, 0.941] | 0.882 [0.857, 0.902] | 0.943 [0.929, 0.956] | 0.83 [0.79, 0.87] |
| QTRP S1 pH5 | NEG_B | 277 / 1,158 | 922 (99) | 0.571 [0.528, 0.613] | 0.662 [0.630, 0.692] | 0.620 [0.528, 0.708] | 0.792 [0.729, 0.847] | 0.44 [0.17, 0.70] |
| QTRP S2 pH5 | NEG_A | 80 / 391 | 56 (53) | 0.639 [0.562, 0.716] | 0.732 [0.675, 0.788] | 0.665 [0.561, 0.757] | 0.725 [0.619, 0.819] | 0.60 [0.31, 0.88] |
| QTRP S2 pH5 | NEG_B | 80 / 1,105 | 765 (33) | 0.433 [0.367, 0.499] | 0.574 [0.507, 0.638] | 0.419 [0.284, 0.577] | 0.535 [0.400, 0.675] | −0.91 [−3.5, −0.0047]* |
| qPerS-SID tier B | NEG_A | 724 / 4,057 | 504 (481) | 0.791 [0.775, 0.805] | 0.882 [0.869, 0.893] | 0.805 [0.785, 0.827] | 0.871 [0.856, 0.885] | 0.76 [0.73, 0.80] |
| qPerS-SID tier B | NEG_B | 724 / 1,294 | 1,072 (256) | 0.618 [0.592, 0.644] | 0.701 [0.676, 0.724] | 0.648 [0.599, 0.700] | 0.732 [0.696, 0.767] | 0.59 [0.47, 0.71] |
| Cys-BOOST HeLa | NEG_A | 8,184 / 38,843 | 3,541 (3,340) | 0.738 [0.732, 0.744] | 0.815 [0.808, 0.821] | 0.762 [0.754, 0.769] | 0.808 [0.800, 0.815] | 0.76 [0.74, 0.78] |
| Cys-BOOST HeLa | NEG_B | 8,184 / 15,457 | 6,276 (2,640) | 0.530 [0.523, 0.538] | 0.608 [0.598, 0.618] | 0.546 [0.533, 0.560] | 0.634 [0.620, 0.647] | 0.28 [0.21, 0.35] |
| Cys-BOOST SH-SY5Y | NEG_A | 2,134 / 15,474 | 1,405 (1,358) | 0.711 [0.700, 0.722] | 0.796 [0.785, 0.806] | 0.728 [0.712, 0.743] | 0.748 [0.731, 0.763] | 0.71 [0.68, 0.75] |
| Cys-BOOST SH-SY5Y | NEG_B | 2,134 / 403 | 1,405 (324) | 0.575 [0.548, 0.603] | 0.913 [0.901, 0.924] | 0.542 [0.511, 0.577] | 0.737 [0.697, 0.770] | 0.18 [0.12, 0.25] |
| FAT-switch | NEG_A | 1,985 / 12,541 | 1,463 (1,416) | 0.839 [0.831, 0.848] | 0.887 [0.880, 0.894] | 0.852 [0.840, 0.862] | 0.880 [0.870, 0.889] | 0.88 [0.86, 0.89] |
| FAT-switch | NEG_B | no observed-unmodified class | | | | | | |
| PAT-switch | NEG_A | 5,736 / 17,182 | 2,940 (2,690) | 0.759 [0.751, 0.766] | 0.852 [0.845, 0.858] | 0.785 [0.775, 0.794] | 0.847 [0.839, 0.855] | 0.74 [0.72, 0.76] |
| PAT-switch | NEG_B | 5,736 / 589 | 3,170 (281) | 0.624 [0.602, 0.646] | 0.816 [0.794, 0.837] | 0.564 [0.521, 0.609] | 0.824 [0.797, 0.849] | 0.39 [0.32, 0.46] |
| YAP1C reporter | NEG_A | 1,745 / 6,441 | 1,090 (978) | 0.726 [0.713, 0.739] | 0.809 [0.794, 0.823] | 0.752 [0.735, 0.768] | 0.811 [0.796, 0.826] | 0.73 [0.69, 0.77] |
| YAP1C reporter | NEG_B | 1,745 / 411 | 1,277 (149) | 0.772 [0.746, 0.796] | 0.803 [0.778, 0.827] | 0.765 [0.710, 0.817] | 0.843 [0.794, 0.887] | 0.90 [0.83, 0.96] |

\*Computed over the replicates in which the DIG25 AUC was at least 0.52 (see Interval note).

**Reading.**
- Against NEG_A, the visibility features alone separate sites from the other cysteines of the same proteins in every cohort, globally and within proteins, without any chemical or structural feature.
- NEG_B negatives were detected, so visibility features are expected to lose discrimination; VIS10 global AUC fell by 0.13–0.28 in six cohorts. Because 90–97% of units carry both classes under NEG_A but only 4–42% under NEG_B, the within-protein AUCs are the closer comparison between the two negative sets.
- Among detected cysteines (NEG_B), DIG25 within-protein AUC is 0.63–0.84 in six cohorts, with every interval above 0.5 (QTRP S2: 0.53 [0.40, 0.68]). The columns DIG25 adds carry this: CPL15 within-protein AUC is 0.61–0.86, against VIS10 0.54–0.65 in five of these cohorts and 0.76 in YAP1C.

**Table S20.2. Composition and within-protein ratios (point estimates).**

| Cohort | Neg. set | units with both classes | DIG25 global − within | ratio within protein | CPL15 / DIG25 above chance |
|---|---|---|---|---|---|
| QTRP S1 pH5 | NEG_A | 195 / 201 | −0.015 | 0.86 | 0.96 |
| QTRP S1 pH5 | NEG_B | 99 / 922 | −0.130 | 0.41 | 0.82 |
| QTRP S2 pH5 | NEG_A | 53 / 56 | +0.008 | 0.73 | 1.00 |
| QTRP S2 pH5 | NEG_B | 33 / 765 | +0.039 | −2.33 (denominator 0.035) | 1.27 |
| qPerS-SID tier B | NEG_A | 481 / 504 | +0.011 | 0.82 | 0.95 |
| qPerS-SID tier B | NEG_B | 256 / 1,072 | −0.031 | 0.64 | 0.84 |
| Cys-BOOST HeLa | NEG_A | 3,340 / 3,541 | +0.007 | 0.85 | 0.91 |
| Cys-BOOST HeLa | NEG_B | 2,640 / 6,276 | −0.026 | 0.35 | 0.91 |
| Cys-BOOST SH-SY5Y | NEG_A | 1,358 / 1,405 | +0.048 | 0.92 | 0.92 |
| Cys-BOOST SH-SY5Y | NEG_B | 324 / 1,405 | +0.176 | 0.18 | 0.99 |
| FAT-switch | NEG_A | 1,416 / 1,463 | +0.007 | 0.93 | 0.92 |
| PAT-switch | NEG_A | 2,690 / 2,940 | +0.005 | 0.82 | 0.93 |
| PAT-switch | NEG_B | 281 / 3,170 | −0.008 | 0.20 | 0.93 |
| YAP1C reporter | NEG_A | 978 / 1,090 | −0.0025 | 0.81 | 0.95 |
| YAP1C reporter | NEG_B | 149 / 1,277 | −0.040 | 0.77 | 0.98 |

**Interval note.** For QTRP S2 under NEG_B, replicates whose DIG25 AUC fell below 0.52 were excluded from the ratio's interval (more than 2.5% of replicates), so that interval is conditional. The within-protein ratios and the CPL15/DIG25 ratios have no intervals because the per-unit bootstrap draws were not stored.

**Verification.**
- All 15 tabulated ratios recompute from the stored AUCs within rounding (largest deviation 7.4 × 10⁻⁴).
- All Source Data values and all values quoted in Results match the analysis table.

**Source.** `results/ptm_detectability_share.csv` (sha256 732fd67a…); feature code `scripts/run_cross_protease_detectability_probe.py` (sha256 f9b990dd…, LF-normalized).
