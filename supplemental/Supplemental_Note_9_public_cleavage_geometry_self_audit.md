# Supplemental Note 9. Self-audit of the ranker: cleavage geometry within label strata, both directions

This is the disclosure the manuscript makes about its own model rather than about the literature. The ranker is trained on one public cohort and scored once, zero-shot, on the other (Supplemental Note 2); nothing is refit on the audited cohort. The association between score and cleavage geometry is measured within each label stratum, where the label is constant, so the association cannot be inherited from it. Models named `v2h_*` are trained on human PXD044043 and audited on rice PXD072089; models named `v2r_*` are trained on rice and audited on human. Labels are our own mapping of each deposit's identified peptides, because neither deposit's MaxQuant output contains a cysteine site table: a cysteine is labeled when it is the only cysteine of a peptide identified in the enriched arm (NEG_A: an unlabeled cysteine is not known to be unmodified; Supplemental Notes 2, 19 and 21).

**Provenance.** Only public data are used in this note. The instrument files are byte-identical to the copies whose correctness was established before this study's self-audit was first run (checked by SHA-256).

**Association between score and cleavage band** (`attribute_rate` is the rate of the cleavage attribute in the stratum; intervals are 95% protein-clustered bootstrap intervals, 5,000 replicates).

| cohort | model | band | stratum | n | attribute_rate | auc | ci_low | ci_high |
|---|---|---|---|---|---|---|---|---|
| rice_PXD072089 | v2h_chem_trained_on_human | distal_6_12 | all | 5004 | 0.76 | 0.6578 | 0.6379 | 0.678 |
| rice_PXD072089 | v2h_chem_trained_on_human | distal_6_12 | label_positive | 737 | 0.8453 | 0.6298 | 0.5733 | 0.6858 |
| rice_PXD072089 | v2h_chem_trained_on_human | distal_6_12 | label_negative | 4267 | 0.7453 | 0.6474 | 0.6255 | 0.6693 |
| rice_PXD072089 | v2h_chem_trained_on_human | proximal_le_5 | all | 5004 | 0.6441 | 0.6022 | 0.5821 | 0.6218 |
| rice_PXD072089 | v2h_chem_trained_on_human | proximal_le_5 | label_positive | 737 | 0.7137 | 0.5169 | 0.4688 | 0.5665 |
| rice_PXD072089 | v2h_chem_trained_on_human | proximal_le_5 | label_negative | 4267 | 0.6321 | 0.6055 | 0.5841 | 0.6264 |
| rice_PXD072089 | v2h_chem_trained_on_human | kr_sfe006_offsets | all | 5004 | 0.7008 | 0.6632 | 0.6446 | 0.6825 |
| rice_PXD072089 | v2h_chem_trained_on_human | kr_sfe006_offsets | label_positive | 737 | 0.7856 | 0.615 | 0.5647 | 0.6656 |
| rice_PXD072089 | v2h_chem_trained_on_human | kr_sfe006_offsets | label_negative | 4267 | 0.6862 | 0.6609 | 0.64 | 0.6811 |
| rice_PXD072089 | v2h_full_trained_on_human | distal_6_12 | all | 5004 | 0.76 | 0.6564 | 0.6365 | 0.6768 |
| rice_PXD072089 | v2h_full_trained_on_human | distal_6_12 | label_positive | 737 | 0.8453 | 0.6313 | 0.5694 | 0.693 |
| rice_PXD072089 | v2h_full_trained_on_human | distal_6_12 | label_negative | 4267 | 0.7453 | 0.6457 | 0.6239 | 0.6669 |
| rice_PXD072089 | v2h_full_trained_on_human | proximal_le_5 | all | 5004 | 0.6441 | 0.642 | 0.6232 | 0.6605 |
| rice_PXD072089 | v2h_full_trained_on_human | proximal_le_5 | label_positive | 737 | 0.7137 | 0.5442 | 0.4977 | 0.5921 |
| rice_PXD072089 | v2h_full_trained_on_human | proximal_le_5 | label_negative | 4267 | 0.6321 | 0.6545 | 0.6347 | 0.6737 |
| rice_PXD072089 | v2h_full_trained_on_human | kr_sfe006_offsets | all | 5004 | 0.7008 | 0.6473 | 0.6278 | 0.6662 |
| rice_PXD072089 | v2h_full_trained_on_human | kr_sfe006_offsets | label_positive | 737 | 0.7856 | 0.5823 | 0.5291 | 0.6355 |
| rice_PXD072089 | v2h_full_trained_on_human | kr_sfe006_offsets | label_negative | 4267 | 0.6862 | 0.6458 | 0.6252 | 0.6658 |
| human_PXD044043 | v2r_chem_trained_on_rice | distal_6_12 | all | 18567 | 0.7901 | 0.609 | 0.5981 | 0.62 |
| human_PXD044043 | v2r_chem_trained_on_rice | distal_6_12 | label_positive | 3206 | 0.8974 | 0.5906 | 0.5559 | 0.6231 |
| human_PXD044043 | v2r_chem_trained_on_rice | distal_6_12 | label_negative | 15361 | 0.7677 | 0.5813 | 0.5693 | 0.5931 |
| human_PXD044043 | v2r_chem_trained_on_rice | proximal_le_5 | all | 18567 | 0.6683 | 0.5332 | 0.5221 | 0.5444 |
| human_PXD044043 | v2r_chem_trained_on_rice | proximal_le_5 | label_positive | 3206 | 0.757 | 0.4476 | 0.4252 | 0.4701 |
| human_PXD044043 | v2r_chem_trained_on_rice | proximal_le_5 | label_negative | 15361 | 0.6498 | 0.5255 | 0.5133 | 0.538 |
| human_PXD044043 | v2r_chem_trained_on_rice | kr_sfe006_offsets | all | 18567 | 0.7236 | 0.5989 | 0.5878 | 0.6097 |
| human_PXD044043 | v2r_chem_trained_on_rice | kr_sfe006_offsets | label_positive | 3206 | 0.8194 | 0.5833 | 0.5575 | 0.6087 |
| human_PXD044043 | v2r_chem_trained_on_rice | kr_sfe006_offsets | label_negative | 15361 | 0.7037 | 0.5795 | 0.5677 | 0.5918 |
| human_PXD044043 | v2r_full_trained_on_rice | distal_6_12 | all | 18567 | 0.7901 | 0.6425 | 0.6312 | 0.6539 |
| human_PXD044043 | v2r_full_trained_on_rice | distal_6_12 | label_positive | 3206 | 0.8974 | 0.5875 | 0.552 | 0.6233 |
| human_PXD044043 | v2r_full_trained_on_rice | distal_6_12 | label_negative | 15361 | 0.7677 | 0.6199 | 0.6076 | 0.6326 |
| human_PXD044043 | v2r_full_trained_on_rice | proximal_le_5 | all | 18567 | 0.6683 | 0.5866 | 0.5752 | 0.5976 |
| human_PXD044043 | v2r_full_trained_on_rice | proximal_le_5 | label_positive | 3206 | 0.757 | 0.4593 | 0.4376 | 0.4812 |
| human_PXD044043 | v2r_full_trained_on_rice | proximal_le_5 | label_negative | 15361 | 0.6498 | 0.5907 | 0.5786 | 0.6022 |
| human_PXD044043 | v2r_full_trained_on_rice | kr_sfe006_offsets | all | 18567 | 0.7236 | 0.6084 | 0.5974 | 0.6195 |
| human_PXD044043 | v2r_full_trained_on_rice | kr_sfe006_offsets | label_positive | 3206 | 0.8194 | 0.5608 | 0.5337 | 0.5884 |
| human_PXD044043 | v2r_full_trained_on_rice | kr_sfe006_offsets | label_negative | 15361 | 0.7037 | 0.596 | 0.5843 | 0.6084 |

**Label's own association with the distal band** (the analytic ceiling of a binary predictor at the observed class rates): rice 0.5344, human 0.5558.

**Null control (permutation distribution).** For each cell, scores are permuted 1,000 times within the stratum and the distal AUC recomputed; the 2.5th-97.5th percentile range of the null must cover 0.5, and the permutation p of the observed AUC is reported. This control was added by a pre-registered amendment after a single-draw permutation -- which excludes 0.5 by chance in about 5% of cells under a correctly behaved estimator -- gave one cell whose interval excluded 0.5 by 0.0004 (v2h_full, label-negative: 0.4803 [0.4610, 0.4996]); that value is kept here. The reverse direction was pre-registered with this control from the start.

| cohort | model | stratum | n | observed_auc | null_mean | null_p2_5 | null_p97_5 | null_covers_half | permutation_p |
|---|---|---|---|---|---|---|---|---|---|
| rice_PXD072089 | v2h_chem_trained_on_human | label_positive | 737 | 0.6298 | 0.4987 | 0.4445 | 0.5573 | True | 0.001 |
| rice_PXD072089 | v2h_chem_trained_on_human | label_negative | 4267 | 0.6474 | 0.5 | 0.4807 | 0.5181 | True | 0.001 |
| rice_PXD072089 | v2h_full_trained_on_human | label_positive | 737 | 0.6313 | 0.5002 | 0.4436 | 0.5552 | True | 0.001 |
| rice_PXD072089 | v2h_full_trained_on_human | label_negative | 4267 | 0.6457 | 0.4998 | 0.4811 | 0.5181 | True | 0.001 |
| human_PXD044043 | v2r_chem_trained_on_rice | label_positive | 3206 | 0.5906 | 0.5008 | 0.4702 | 0.5343 | True | 0.001 |
| human_PXD044043 | v2r_chem_trained_on_rice | label_negative | 15361 | 0.5813 | 0.4997 | 0.4886 | 0.5103 | True | 0.001 |
| human_PXD044043 | v2r_full_trained_on_rice | label_positive | 3206 | 0.5875 | 0.5001 | 0.4697 | 0.5345 | True | 0.001 |
| human_PXD044043 | v2r_full_trained_on_rice | label_negative | 15361 | 0.6199 | 0.5001 | 0.489 | 0.5106 | True | 0.001 |

**Within-protein permutation** (the between-protein floor, not gating): 0.4883 to 0.5336 across the eight cells (values in the released audit records).

**Verdict.** In both training directions the interval lies above 0.5 in both label strata.

**Declared deviation.** Bootstrap clusters are proteins, not homology components; a k-mer homology grouping exists only for the human cohort and was not used. Intervals may be narrower than component clustering would give and must be quoted with this note.
