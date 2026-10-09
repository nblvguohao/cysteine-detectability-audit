# Supplemental Note 22. The published S-nitrosylation predictor re-scored on its released test set (post hoc)

**Status.** Post hoc: added after the primary analyses were complete. Supplemental Note 5 compared two detectability-only models with the test-set AUROCs tabulated in the pLMSNOSite paper (Pratyush et al., 2023). This Note reproduces that analysis, runs the released pLMSNOSite models on the same test set, and compares the detectability-only models with the predictor's own scores.

**Data.** The data come from the pLMSNOSite repository (github.com/KCLabMTU/pLMSNOSite, commit e9158af). Every file was read as committed, from the git object store, so the training and test files carry the SHA-256 values listed in Supplemental Note 5. The training set has 6,552 cysteines (3,276 positives) in 1,918 proteins. The test set has 3,519 cysteines (351 positives) in 278 proteins, 4 of which also occur in training; its prevalence is 0.0997.

**Models.**
- *Detectability-only models.* VIS10 (the ten peptide-visibility features) and DIG25 (all 25 theoretical-digest features) are those defined in Supplemental Note 5, computed from the full sequence under the tryptic rule. Each set was fitted once to the whole training set with the fixed configuration of that Note, without tuning or calibration: a scikit-learn 1.7.2 HistGradientBoostingClassifier with 300 iterations, learning rate 0.06, at most 31 leaves, at least 20 samples per leaf, no L2 penalty, 255 bins, no early stopping and random state 0.
- *Released predictor.* The released ProtT5-feature base model, embedding base model and stacking model were run on the released test inputs (the authors' precomputed per-site ProtT5 features and 37-residue windows) with the authors' evaluation script, evaluate_model.py, unmodified. They ran under TensorFlow 2.15.1 (Keras 2.15.0); the models had been saved with Keras 2.8.0. Plotting and progress-bar libraries that the script imports but never uses were supplied as empty stand-ins. The stacking model's probabilities, which the script thresholds at 0.5, were kept; at every site the script's call equals probability > 0.5.
- *Check without TensorFlow.* A numpy re-computation of the three released models from their weights agrees with the TensorFlow probabilities to within 4 × 10⁻⁷ and gives the same calls, AUROC and AUPRC.

Source: `posthoc_2026-09-30/results/G_plmsnosite_rescore/detectability_reproduction.json`, `.../plmsnosite_run.json`, `.../numpy_forward_pass_check.json`.

**Statistics.**
- AUROC counts ties as one half, and AUPRC is step-wise average precision; both equal scikit-learn's roc_auc_score and average_precision_score.
- Intervals are 95% percentile intervals from 5,000 resamples of the 278 test proteins (seed 20260930), drawn once and shared by every score; models were not refit within resamples. Differences and shares are computed within each resample, so their intervals are paired. A difference is called resolved when its interval excludes zero; an unresolved difference is not evidence that two models are equal.
- A detectability model's share of pLMSNOSite's above-chance value is (model − chance)/(pLMSNOSite − chance), where chance is 0.5 for AUROC and, for AUPRC, the prevalence of the resample.
- Within-protein AUROC is defined as in Supplemental Note 12: the AUROC within each test protein that carries both labels (220 proteins, 3,252 sites, 284 positives), averaged with equal weight per protein. The other 58 test proteins carry one label only and do not enter. A protein drawn more than once in a resample counts once per draw.
- As a sensitivity, within-protein AUROC is also computed by pooling positive–negative pairs across the same 220 proteins, the definition of Supplemental Note 20, which weights each protein by its number of pairs.

**Reproduction of Supplemental Note 5.** The analysis of Supplemental Note 5 was re-implemented, and every AUROC, AUPRC and interval that it reports for the test set was reproduced exactly, to the four decimals given there:
- VIS10: AUROC 0.7095 [0.6817, 0.7367], AUPRC 0.1677;
- DIG25: AUROC 0.7765 [0.7231, 0.8274], AUPRC 0.2195;
- DIG25 trained on permuted labels: AUROC 0.5124 [0.471, 0.5546].

The three intervals (seed 20260924) are reproduced when one bootstrap generator serves DIG25, VIS10 and the label-permuted model, in that order. The five training-set folds of the positive control were not reproduced; their fold assignment was not recorded, and an unshuffled protein-grouped split gives a mean AUROC of 0.763 against the 0.765 of Supplemental Note 5.

Source: `posthoc_2026-09-30/results/G_plmsnosite_rescore/detectability_reproduction.csv`, `.../detectability_reproduction.json`, `.../reproduction_diagnostics.json`.

**What the tabulated pLMSNOSite AUROC is.** Run as released, evaluate_model.py prints sensitivity 0.7350 (258 of 351 positives), specificity 0.7724 (2,447 of 3,168 negatives) and "AUROC" 0.7537. The script computes that AUROC with scikit-learn's roc_auc_score from the predictions after thresholding them at 0.5. The AUROC of 0/1 calls is the mean of sensitivity and specificity at that threshold, (0.7350 + 0.7724)/2 = 0.7537, which matches the 0.754 tabulated for pLMSNOSite; the paper's own sensitivity and specificity, 0.735 and 0.773, also average to 0.754. The tabulated value is therefore the AUROC of thresholded calls, not a ranking statistic. The paper's specificity implies one to three more true negatives than this run, so the published run was not bit-identical to this one. Supplemental Note 5 noted that 0.754 coincides with the geometric mean of sensitivity and specificity; the two means agree here (0.7535 and 0.7537) because sensitivity and specificity are similar.

Source: `posthoc_2026-09-30/results/G_plmsnosite_rescore/evaluate_model_verbatim_stdout.txt`, `.../plmsnosite_threshold_metrics.csv`.

**The released predictor on its probabilities.** Scored on the probabilities that the script thresholds away, pLMSNOSite reaches AUROC 0.815 [0.763, 0.858] and AUPRC 0.307 [0.265, 0.362]. Table S22.1 gives every AUROC and AUPRC value of this Note. AUROC was the primary statistic of Supplemental Note 5, so the text reports ranking by AUROC; the detectability models' shares of the predictor's AUPRC gain over prevalence are lower than their AUROC shares (Table S22.1).

**Table S22.1. Discrimination on the released test set** (3,519 cysteines, 351 positives, 278 proteins; prevalence 0.0997). Point estimate [95% interval].

| Model or comparison | AUROC, all test sites | AUPRC, all test sites | Within-protein AUROC, mean of per-protein AUROCsᵃ | Within-protein AUROC, pairs pooled across proteinsᵇ |
|---|---|---|---|---|
| pLMSNOSite, released models (probabilities) | 0.815 [0.763, 0.858] | 0.307 [0.265, 0.362] | 0.702 [0.659, 0.745] | 0.771 [0.731, 0.808] |
| pLMSNOSite, calls thresholded at 0.5ᶜ | 0.754 [0.714, 0.789] | 0.220 [0.192, 0.253] | not computed | not computed |
| VIS10 | 0.709 [0.682, 0.737]ᵈ | 0.168 [0.130, 0.221] | 0.690 [0.647, 0.731] | 0.716 [0.675, 0.754] |
| DIG25 | 0.776 [0.725, 0.829]ᵈ | 0.220 [0.162, 0.322] | 0.684 [0.641, 0.726] | 0.800 [0.753, 0.846] |
| VIS10, share of pLMSNOSite's above-chance valueᵉ | 0.66 [0.56, 0.81] | 0.33 [0.22, 0.51] | 0.94 [0.71, 1.24] | 0.80 [0.60, 1.03] |
| DIG25, share of pLMSNOSite's above-chance valueᵉ | 0.88 [0.69, 1.08] | 0.58 [0.35, 1.02] | 0.91 [0.70, 1.16] | 1.11 [0.91, 1.30] |
| pLMSNOSite − VIS10, paired difference | 0.106 [0.051, 0.153] | 0.139 [0.085, 0.195] | 0.013 [−0.041, 0.065] | 0.055 [−0.007, 0.117] |
| pLMSNOSite − DIG25, paired difference | 0.039 [−0.021, 0.105] | 0.087 [−0.003, 0.152] | 0.019 [−0.028, 0.064] | −0.029 [−0.076, 0.024] |

ᵃ Over the 220 test proteins that carry both labels (3,252 sites, 284 positives), each protein weighted equally; the definition of Supplemental Note 12.
ᵇ Sensitivity: concordant positive–negative pairs within proteins over all such pairs in the same 220 proteins, so that each protein is weighted by its number of pairs; the definition of Supplemental Note 20.
ᶜ The value that evaluate_model.py prints as AUROC; it matches the 0.754 tabulated for pLMSNOSite.
ᵈ From the resamples shared by all scores (seed 20260930). The intervals of Supplemental Note 5 (seed 20260924), reproduced above, are [0.6817, 0.7367] for VIS10 and [0.7231, 0.8274] for DIG25.
ᵉ (model − chance)/(pLMSNOSite − chance); chance is 0.5 for AUROC and the prevalence for AUPRC.

Source: `posthoc_2026-09-30/results/G_plmsnosite_rescore/metrics_main.csv`, `.../paired_differences.csv`, `.../recovery_ratios.csv`, `.../withinprot_mean_summary.csv`; per-site scores in `.../plmsnosite_test_scores.csv` and `.../detectability_test_scores.csv`; per-protein AUROCs in `.../withinprot_mean_per_protein.csv`.

**Detectability-only models against the re-scored predictor.**
- *Across test sites.* The ten visibility features reach 0.66 [0.56, 0.81] of pLMSNOSite's above-chance AUROC, and the 25 detectability features 0.88 [0.69, 1.08]. The difference between pLMSNOSite and the 25 features is not resolved (0.039 [−0.021, 0.105]).
- *Within proteins* (within-protein AUROC as the mean of per-protein AUROCs over the 220 test proteins that carry both labels). pLMSNOSite reaches 0.702 [0.659, 0.745], VIS10 0.690 [0.647, 0.731] and DIG25 0.684 [0.641, 0.726]. The visibility features thus reach 0.94 [0.71, 1.24], and the 25 features 0.91 [0.70, 1.16], of pLMSNOSite's above-chance within-protein AUROC. Neither paired difference from pLMSNOSite is resolved: 0.013 [−0.041, 0.065] for VIS10 and 0.019 [−0.028, 0.064] for DIG25. Resampling only the 220 proteins with both labels, as in Supplemental Note 12, or using a second seed (20260924) moves no interval bound of these within-protein values by more than 0.02.
- *Within proteins, pairs pooled across proteins* (sensitivity). The shares are 0.80 [0.60, 1.03] for VIS10 and 1.11 [0.91, 1.30] for DIG25, and neither difference from pLMSNOSite is resolved (0.055 [−0.007, 0.117] and −0.029 [−0.076, 0.024]). The two definitions weight proteins differently: the 78 proteins with at most four positive–negative pairs carry 35% of the weight of the mean but 5% of the pairs.

Source: `posthoc_2026-09-30/results/G_plmsnosite_rescore/metrics_main.csv`, `.../paired_differences.csv`, `.../recovery_ratios.csv`, `.../withinprot_mean_summary.csv`, `.../withinprot_mean.json`, `.../withinprot_mean_by_size.csv`.

**Relation to the ratios of Supplemental Note 5.** Supplemental Note 5 took the three AUROCs tabulated for this test set as denominators, which gave the ten visibility features 0.82 to 0.91 of the above-chance value. The pLMSNOSite entry among them is the AUROC of thresholded calls, and the other two, 0.756 for PreSNO and 0.731 for DeepNitro, were not re-scored. This Note therefore uses the re-scored AUROC as the denominator, a change made post hoc, after the primary analysis. Against it, the 25-feature ratio of 0.88 still meets the criterion of the first branch of Supplemental Note 5 (a ratio of at least 0.75 with an AUROC interval above 0.5), although the ratio's own interval, 0.69 to 1.08, extends below 0.75.

Source: `posthoc_2026-09-30/results/G_plmsnosite_rescore/recovery_ratios.csv`.

**Reading.** On this benchmark as released, the 25 detectability features reach 0.88 of the released predictor's above-chance AUROC across test sites, with no resolved difference from it, and the ten visibility features alone reach two-thirds. Within proteins, the two feature sets reach 0.91 to 0.94 of the predictor's above-chance AUROC (mean of per-protein AUROCs), and neither differs resolvably from it under either within-protein definition. The comparison concerns what the benchmark rewards, not S-nitrosylation chemistry: the benchmark's negatives are other cysteines of S-nitrosylated proteins, and we do not claim that a detectability model is a better S-nitrosylation predictor.

**Caveats.**
- The benchmark is split by site, not by protein, and its negatives are the unmatched background examined in the main text.
- The re-scoring uses the authors' precomputed ProtT5 test features. PreSNO and DeepNitro were not re-scored.
- Intervals are conditional on the fitted models, which were not refit within resamples.

**Provenance.** All values are read from stored outputs under `posthoc_2026-09-30/results/G_plmsnosite_rescore/`, whose `provenance.json` records the SHA-256 of every input, script and output.
