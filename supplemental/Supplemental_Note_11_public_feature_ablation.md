# Supplemental Note 11. Removing the declared composition columns does not remove the ordering (public cohorts)

Three arms of the chemistry set were refit from scratch on each public cohort with identical folds, member grids and seeds, differing only in the column set (all 1,046 chemistry columns; minus the win7/win10 basic-composition columns; minus all twenty basic-composition columns). Two further arms serve Supplemental Note 12 (a 25-column detectability-only set and the full 1,107-column set) and one is a control (the chemistry columns with labels permuted across proteins). The verdict statistic and its reading rule were fixed before any arm was fitted.

Public cohorts: human persulfidation deposit PXD044043 and rice persulfidation deposit PXD072089. The arms are refit within each cohort's own folds, on all 18,567 human and all 5,598 rice cysteines, and do not read the ranker audited in Supplemental Notes 2, 9 and 10. In the arm tables, rice values are computed on the 5,004 cysteines whose protein is in the reference proteome file used (578 proteins carry both labeled and unlabeled cysteines), whereas the rice paired differences are computed over the 637 proteins that carry both classes among all 5,598; human values use all 18,567 cysteines (1,475 proteins). The labels are our own mapping of each deposit's identified peptides: a cysteine is labeled when it is the only cysteine of a peptide identified in the enriched arm (a NEG_A negative: an unlabeled cysteine is not known to be unmodified; Supplemental Notes 2, 19 and 21). Bootstrap clusters are proteins, not homology components; a k-mer homology grouping exists only for the human cohort and was not used, so intervals may be narrower than component clustering would give. The five structural columns are missing for external sequences, so absolute performance is not comparable across cohorts, only across arms within one.

**Which score is read.** Every arm fits three members on the same folds: a LightGBM ranking model, a LightGBM binary model and a small neural network, blended with weights chosen on the validation folds. The neural network receives, besides the arm's own columns, a 1,024-dimensional ProtBERT embedding of each cysteine (the frozen last-layer token of the cysteine in a 31-residue window centered on it, shortened at the protein termini), in raw and protein-centered form. The blend of any arm therefore carries language-model information that the arm's column list does not describe, and a "detectability-only" blend is not detectability-only. The registered analysis read the blend. This note reports the LightGBM ranking member, which sees only the arm's own columns, as the primary reading, and the registered blend as the record. The binary LightGBM member fails the registered label-permutation control on the human cohort (the label-permuted arm's top-100 distal log2 odds ratio excludes zero), so it is not read. The choice of member to report was made after the registered analysis had been read and is post hoc; nothing was refitted, and the member scores are those stored by the registered run, re-read by the registered analysis code (which reproduces the stored blend outputs byte for byte).

**Arms, verdict statistic (top-100 distal log2 odds ratio) and within-protein AUC: LightGBM ranking member (primary).**

| cohort | arm | n | top100_distal_rate | rate_rest | log2_or | ci_low | ci_high | within_protein_auc | n_proteins_scored |
|---|---|---|---|---|---|---|---|---|---|
| human_PXD044043 | control | 18567 | 0.93 | 0.7894 | 1.7343 | 0.8724 | 3.1774 | 0.7935 | 1475 |
| human_PXD044043 | drop_w710 | 18567 | 0.95 | 0.7892 | 2.2132 | 1.2516 | 4.1207 | 0.7953 | 1475 |
| human_PXD044043 | drop_all_basic | 18567 | 0.97 | 0.7891 | 2.8961 | 1.8266 | 5.6708 | 0.7845 | 1475 |
| human_PXD044043 | label_permuted | 18567 | 0.83 | 0.7899 | 0.344 | -0.3222 | 1.2049 | 0.4933 | 1071 |
| human_PXD044043 | detect_only | 18567 | 0.91 | 0.7895 | 1.3611 | 0.5385 | 2.5604 | 0.8819 | 1475 |
| human_PXD044043 | full | 18567 | 0.93 | 0.7894 | 1.7343 | 0.8707 | 3.1456 | 0.8884 | 1475 |
| rice_PXD072089 | control | 5004 | 0.95 | 0.7561 | 2.486 | 1.5014 | 4.2985 | 0.7663 | 578 |
| rice_PXD072089 | drop_w710 | 5004 | 0.93 | 0.7565 | 2.0048 | 1.0923 | 3.4496 | 0.7636 | 578 |
| rice_PXD072089 | drop_all_basic | 5004 | 0.91 | 0.7569 | 1.6294 | 0.7985 | 2.8293 | 0.7818 | 578 |
| rice_PXD072089 | label_permuted | 5004 | 0.79 | 0.7594 | 0.229 | -0.4839 | 1.085 | 0.5128 | 384 |
| rice_PXD072089 | detect_only | 5004 | 0.96 | 0.7559 | 2.7921 | 1.7423 | 4.6079 | 0.8023 | 578 |
| rice_PXD072089 | full | 5004 | 0.9 | 0.7571 | 1.4675 | 0.6684 | 2.5763 | 0.901 | 578 |

All six ablation intervals exclude zero in both cohorts, and the label-permuted control crosses zero in both. Zero was the registered reference; Supplemental Note 21 (post hoc) compares these values with the level that the labels alone produce.

**Performance cost of deleting columns (post hoc, not pre-registered), LightGBM ranking member.** Paired protein-clustered bootstrap of the within-protein AUC difference against the full chemistry arm. 'Negligible' is deliberately not defined: no smallest meaningful AUC difference was registered.

| cohort | arm | n_proteins | mean_paired_diff_vs_control | ci_low | ci_high | reading |
|---|---|---|---|---|---|---|
| human_PXD044043 | drop_w710 | 1475 | 0.00177 | -0.00587 | 0.00956 | none |
| human_PXD044043 | drop_all_basic | 1475 | -0.00901 | -0.01703 | -0.00047 | cost |
| rice_PXD072089 | drop_w710 | 637 | -0.00244 | -0.0184 | 0.01289 | none |
| rice_PXD072089 | drop_all_basic | 637 | 0.01108 | -0.00395 | 0.02632 | none |

**Registered reading (blend of three members, including the embedding-fed network), kept as the record.**

| cohort | arm | n | top100_distal_rate | rate_rest | log2_or | ci_low | ci_high | within_protein_auc | n_proteins_scored |
|---|---|---|---|---|---|---|---|---|---|
| human_PXD044043 | control | 18567 | 0.96 | 0.7892 | 2.5182 | 1.4693 | 4.2822 | 0.8007 | 1475 |
| human_PXD044043 | drop_w710 | 18567 | 0.95 | 0.7892 | 2.2132 | 1.2816 | 4.0502 | 0.7968 | 1475 |
| human_PXD044043 | drop_all_basic | 18567 | 0.99 | 0.789 | 4.1488 | 2.7769 | 5.9529 | 0.7889 | 1475 |
| human_PXD044043 | label_permuted | 18567 | 0.83 | 0.7899 | 0.344 | -0.3222 | 1.2049 | 0.5021 | 1071 |
| human_PXD044043 | detect_only | 18567 | 0.9 | 0.7895 | 1.2004 | 0.4183 | 2.3721 | 0.8851 | 1475 |
| human_PXD044043 | full | 18567 | 0.97 | 0.7891 | 2.8961 | 1.7674 | 5.6986 | 0.887 | 1475 |
| rice_PXD072089 | control | 5004 | 0.84 | 0.7584 | 0.7069 | 0.0199 | 1.5736 | 0.8157 | 578 |
| rice_PXD072089 | drop_w710 | 5004 | 0.86 | 0.758 | 0.9302 | 0.2157 | 1.8554 | 0.8066 | 578 |
| rice_PXD072089 | drop_all_basic | 5004 | 0.85 | 0.7582 | 0.8156 | 0.0775 | 1.8229 | 0.8116 | 578 |
| rice_PXD072089 | label_permuted | 5004 | 0.69 | 0.7614 | -0.5321 | -1.1497 | 0.1918 | 0.5023 | 384 |
| rice_PXD072089 | detect_only | 5004 | 0.94 | 0.7563 | 2.2282 | 1.3277 | 3.7541 | 0.8371 | 578 |
| rice_PXD072089 | full | 5004 | 0.91 | 0.7569 | 1.6294 | 0.8088 | 2.8457 | 0.8953 | 578 |

| cohort | arm | n_proteins | mean_paired_diff_vs_control | ci_low | ci_high | reading |
|---|---|---|---|---|---|---|
| human_PXD044043 | drop_w710 | 1475 | -0.00393 | -0.01068 | 0.00286 | none |
| human_PXD044043 | drop_all_basic | 1475 | -0.01189 | -0.01937 | -0.00442 | cost |
| rice_PXD072089 | drop_w710 | 637 | -0.00553 | -0.01841 | 0.0072 | none |
| rice_PXD072089 | drop_all_basic | 637 | -0.00482 | -0.01531 | 0.00604 | none |
