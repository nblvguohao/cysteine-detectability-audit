# Supplemental Note 12. Detectability is most of this task (public cohorts)

A 25-column set containing no chemistry, answering only whether a peptide is identifiable, was refit on the same folds as the full 1,107-column model. The recovery ratio is the detectability-only within-protein AUC above 0.5 over the full model's above 0.5, with a paired protein-clustered bootstrap interval. The reading threshold, 0.5, was fixed before fitting.

Public cohorts: human persulfidation deposit PXD044043 and rice persulfidation deposit PXD072089. The arms are refit within each cohort's own folds, on all 18,567 human and all 5,598 rice cysteines, and do not read the ranker audited in Supplemental Notes 2, 9 and 10; rice values below are computed on the 5,004 cysteines whose protein is in the reference proteome file used (578 proteins carry both labeled and unlabeled cysteines), and human values on all 18,567 (1,475 proteins). The labels are our own mapping of each deposit's identified peptides: a cysteine is labeled when it is the only cysteine of a peptide identified in the enriched arm (a NEG_A negative: an unlabeled cysteine is not known to be unmodified; Supplemental Notes 2, 19 and 21). Bootstrap clusters are proteins, not homology components; a k-mer homology grouping exists only for the human cohort and was not used, so intervals may be narrower than component clustering would give. The five structural columns are missing for external sequences, so absolute performance is not comparable across cohorts, only across arms within one.

**Which score is read.** Every arm fits three members on the same folds: a LightGBM ranking model, a LightGBM binary model and a small neural network, blended with weights chosen on the validation folds. The neural network receives, besides the arm's own columns, a 1,024-dimensional ProtBERT embedding of each cysteine (the frozen last-layer token of the cysteine in a 31-residue window centered on it, shortened at the protein termini), in raw and protein-centered form. The blend of any arm therefore carries language-model information that the arm's column list does not describe, and a "detectability-only" blend is not detectability-only. The registered analysis read the blend. This note reports the LightGBM ranking member, which sees only the arm's own columns, as the primary reading, and the registered blend as the record. The binary LightGBM member fails the registered label-permutation control on the human cohort (the label-permuted arm's top-100 distal log2 odds ratio excludes zero), so it is not read. The choice of member to report was made after the registered analysis had been read and is post hoc; nothing was refitted, and the member scores are those stored by the registered run, re-read by the registered analysis code (which reproduces the stored blend outputs byte for byte).

**LightGBM ranking member (primary).**

| cohort | n_proteins | wp_auc_full | wp_auc_detect_only | recovery_ratio | ci_low | ci_high |
|---|---|---|---|---|---|---|
| human_PXD044043 | 1475 | 0.8884 | 0.8819 | 0.9832 | 0.9655 | 1.0016 |
| rice_PXD072089 | 578 | 0.901 | 0.8023 | 0.7538 | 0.7036 | 0.8049 |

**Registered reading (blend, including the embedding-fed network), kept as the record.**

| cohort | n_proteins | wp_auc_full | wp_auc_detect_only | recovery_ratio | ci_low | ci_high |
|---|---|---|---|---|---|---|
| human_PXD044043 | 1475 | 0.887 | 0.8851 | 0.9952 | 0.9779 | 1.0133 |
| rice_PXD072089 | 578 | 0.8953 | 0.8371 | 0.8527 | 0.8069 | 0.9 |

The embedding member accounts for much of the rice figure (0.8527 in the blend against 0.7538 in the ranking member) and little of the human one (0.9952 against 0.9832). Because these labels are close to 'was this residue reported', a high recovery ratio is the expected result and is not evidence that no chemical information exists.
