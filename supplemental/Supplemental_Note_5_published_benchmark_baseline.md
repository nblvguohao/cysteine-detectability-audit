# Supplemental Note 5. A detectability-only baseline on a published S-nitrosylation predictor benchmark

**Question (registered before any feature was computed).** How much of the published independent-test AUROC on the pLMSNOSite S-nitrosylation benchmark can a model reach that sees only in-silico tryptic-digest features of the site's peptide, with no residue identity, evolutionary, structural or language-model information?

**Registration.** The question, the data and their hashes, the features, the model, the primary
statistic, the three reading branches and the two controls were recorded before the analysis script
existed, and the script refuses to run if either released data file differs from the hashes given below.
The branches are stated in full below, so the registration can be checked against the analysis as
reported. The protocol fixes the data, features, model, primary statistic, three reading branches and two controls, and commits to reporting whichever branch results.

**Data.** https://github.com/KCLabMTU/pLMSNOSite, commit e9158af06418f3e50c9bbbbeca77ea692e6507ce, data/train/sequence_train.csv and data/test/sequence_test.csv (Apache-2.0). SHA-256: `sequence_train.csv` `a3b4a1fe3858ed3e4a0aa62d47e52c8e1c118cceb06c4edcf75c0b82ff96bd0f`; `sequence_test.csv` `d4498e6f2afa28b89d150dacf44b1de6bbfd791e772f1fefe120aa4aabcc5551`. Used as released: training 3,276 positives and 3,276 negatives (1,918 proteins); test 351 positives and 3,168 negatives (278 proteins); 4 proteins occur in both. Every listed position is a cysteine in the released sequence. Declared discrepancy: the paper states 3,383 training SNO sites; the released balanced training file has 3,276. The released files are used and the published test-set numbers are taken as published.

**Features.** DIG25, 25 columns computed from the full sequence under the tryptic rule (cleavage after K or R, not before P; detectable peptide length 7–30 residues and mass 700–3,500 Da; up to two missed cleavages): `pep_len`, `pep_mass`, `pep_gravy`, `pep_cys_count`, `pep_offset_in_peptide`, `pep_relative_offset`, `pep_is_protein_n_terminal`, `pep_is_protein_c_terminal`, `pep_detectable_length`, `pep_detectable_mass`, `pep_detectable_both`, `pep_log_len`, `pep_mc1_left_len`, `pep_mc1_right_len`, `pep_mc1_min_len`, `pep_mc1_detectable`, `pep_mc2_left_len`, `pep_mc2_right_len`, `pep_mc2_min_len`, `pep_mc2_detectable`, `pep_detectable_any_missed_cleavage`, `pep_n_peptides_in_protein`, `pep_cleavage_sites_within_20`, `pep_dist_to_peptide_n_term`, `pep_dist_to_peptide_c_term`. VIS10, the 10-column visibility-only subset: `pep_len`, `pep_mass`, `pep_gravy`, `pep_detectable_length`, `pep_detectable_mass`, `pep_detectable_both`, `pep_log_len`, `pep_mc1_detectable`, `pep_mc2_detectable`, `pep_detectable_any_missed_cleavage`. VIS10 omits the position of the cysteine in its peptide and `pep_cleavage_sites_within_20`, the number of K and R within 20 residues, because a local basic-residue count can also be a genuine chemical feature (an acid–base motif) and would let chemistry pass as detectability.

**Model.** A scikit-learn HistGradientBoostingClassifier with the fixed configuration of the cross-protease detectability probe, fitted once to the full training set, without tuning or calibration. Primary statistic: DIG25 test AUROC with a 95% protein-clustered bootstrap interval (5,000 resamples of test proteins, seed 20260924). Recovery ratio: (AUROC − 0.5)/(published AUROC − 0.5).

**Branches (registered).** R >= 0.75 and the AUROC interval's lower bound > 0.5: detectability-only features recover most of the published discrimination on this benchmark (B1); 0.25 <= R < 0.75 and the lower bound > 0.5: detectability-only features recover a substantial minority of it (B2); R < 0.25, or the interval includes 0.5: the published discrimination on this benchmark is largely not reproduced by detectability features (B3). **Observed: B1.**

**Results.**

| model | n_features | test_auroc | ci_low | ci_high | test_auprc | recovery_vs_pLMSNOSite | recovery_vs_PreSNO | recovery_vs_DeepNitro |
|---|---|---|---|---|---|---|---|---|
| VIS10 | 10 | 0.7095 | 0.6817 | 0.7367 | 0.1677 | 0.8248 | 0.8184 | 0.9069 |
| DIG25 | 25 | 0.7765 | 0.7231 | 0.8274 | 0.2195 | 1.0886 | 1.0801 | 1.197 |
| DIG25 label-permuted (NC) | 25 | 0.5124 | 0.471 | 0.5546 |  |  |  |  |
| pLMSNOSite |  | 0.754 |  |  |  |  |  |  |
| PreSNO |  | 0.756 |  |  |  |  |  |  |
| DeepNitro |  | 0.731 |  |  |  |  |  |  |

Test prevalence is 351/3,519, so the chance level of AUPRC is 0.0997. Published AUROCs are those tabulated for the independent test set by Pratyush et al. (2023), with the values for PreSNO and DeepNitro adopted there from the PreSNO study; they were quoted, not recomputed.

**Controls.** Negative control: the same DIG25 model fitted after permuting the training labels gave test AUROC 0.5124 [0.471, 0.5546], whose interval covers 0.5 as required. Positive control: five-fold cross-validation on the training set, grouped by protein, gave AUROC 0.7669, 0.8344, 0.7412, 0.7371, 0.7467 (mean 0.7653).

**Reading.** The registered primary statistic, the recovery ratio of the 25-feature model against pLMSNOSite's tabulated AUROC, is 1.0886, above one. We do not lead with it. The value tabulated for pLMSNOSite in the AUROC column of that paper's independent-test table, 0.754, coincides with the geometric mean of sensitivity (0.735) and specificity (0.773) reported in the same paper's text, whereas the values tabulated for PreSNO and DeepNitro do not coincide with theirs; we therefore cannot verify that 0.754 is a score-based AUROC. Anchored instead on the highest tabulated value, PreSNO's 0.756, the ten visibility features recover 0.8184 of the above-chance AUROC and the 25 features 1.0801. Changing the anchor after reading the source paper is a departure from the registration, declared here. The comparison is confined to AUROC: the published precision-recall performance appears only as a figure, and at a test prevalence of 0.0997 our AUPRC is 0.1677 (ten features) and 0.2195 (25 features). We do not claim that a detectability model is a better S-nitrosylation predictor, only that on this benchmark, as constructed, its ranking performance approaches the published one.

**Resolved in revision (post hoc; Supplemental Note 22).** The authors' evaluation script computes the tabulated 0.754 from predictions thresholded at 0.5, as the mean of sensitivity and specificity, so it is not a ranking statistic. Run on the released test set and scored on their probabilities, the released pLMSNOSite models reach a test AUROC of 0.815 [0.763, 0.858]; against it the 25 features recover 0.88 [0.69, 1.08] of the above-chance AUROC across test sites, and within proteins the ten visibility features recover 0.94 [0.71, 1.24]. The main text uses the re-scored value as the denominator.

**Caveats (declared in the protocol).** The benchmark split is by site, not by protein; the number of proteins shared by training and test sets is reported. Negatives are other cysteines of SNO-bearing proteins, the construction this study questions; the comparison is about what the benchmark rewards, not about S-nitrosylation chemistry. Published AUROCs are quoted, not recomputed. The count of K and R near the cysteine enters DIG25 but not VIS10, which is why VIS10 is the conservative figure.
