# Supplemental Note 7. Artifact 4 on public deposit PXD048216

PXD048216 is an acyl-biotin exchange study of S-palmitoylation in *Ustilaginoidea virens*, two biological replicates, searched in MaxQuant. Free thiols are blocked with N-ethylmaleimide (NEM); after release the newly freed cysteine is alkylated with iodoacetamide, so a carbamidomethyl (CAM) site is the modification readout and an NEM site a detected cysteine that was free. Observed cysteines are the cysteines of `peptides.txt` peptides on their leading razor protein; a site is the (protein, position) of a row of the CAM or NEM site table with localization probability at least the registered threshold, on that row's leading protein.

**Provenance.** Only public data are used in this note. The criteria were registered before the analysis (`protocols/public_artefact4_preregistration_2026-09-21.json`); the headline coincidence had been computed in an earlier exploratory pass. That pass applied a criterion of a Wilson lower bound of at least 0.95, which the union of both alkylation readouts met (0.9576) and the carbamidomethyl-only readout did not (0.9336); because both values had been seen, it was not used. The registered analysis fixed the readout and set the decision level at 0.90 after those values were known, so neither the 0.90 level nor the corresponding decision line in Figure 2c is blind to this deposit. Bootstrap clusters are proteins (5,000 replicates, registered seed).

**Coincidence: the share of identified cysteines carrying a site assignment.** `PRIMARY` is the registered definition; `SECONDARY` is the earlier asymmetric definition kept to show that it does not move the headline; `SENSITIVITY` drops the localization filter; `PER_REPLICATE` uses each replicate's own localization column.

| definition | reading | observed_cysteines | matched | p | wilson_low | wilson_high | boot_low | boot_high |
|---|---|---|---|---|---|---|---|---|
| PRIMARY_leading_both_sides | CAM_only | 3554 | 3347 | 0.94176 | 0.93357 | 0.94899 | 0.93277 | 0.94999 |
| PRIMARY_leading_both_sides | NEM_only | 3554 | 333 | 0.0937 | 0.08455 | 0.10372 | 0.08053 | 0.10799 |
| PRIMARY_leading_both_sides | union | 3554 | 3427 | 0.96427 | 0.95764 | 0.96988 | 0.95746 | 0.97077 |
| SECONDARY_exploratory_asymmetric | CAM_only | 3554 | 3347 | 0.94176 | 0.93357 | 0.94899 | 0.93277 | 0.94999 |
| SECONDARY_exploratory_asymmetric | NEM_only | 3554 | 333 | 0.0937 | 0.08455 | 0.10372 | 0.08053 | 0.10799 |
| SECONDARY_exploratory_asymmetric | union | 3554 | 3427 | 0.96427 | 0.95764 | 0.96988 | 0.95746 | 0.97077 |
| SENSITIVITY_all_localisation | CAM_only | 3554 | 3362 | 0.94598 | 0.93805 | 0.95294 | 0.93769 | 0.95379 |
| SENSITIVITY_all_localisation | union | 3554 | 3439 | 0.96764 | 0.9613 | 0.97297 | 0.96104 | 0.97362 |
| PER_REPLICATE_WT-1 | CAM_only | 3076 | 2912 | 0.94668 | 0.93817 | 0.95408 | 0.93795 | 0.95494 |
| PER_REPLICATE_WT-2 | CAM_only | 3106 | 2959 | 0.95267 | 0.94463 | 0.9596 | 0.94404 | 0.96118 |

**Partition of the observed cysteines**, stored values:

| key | value |
|---|---|
| observed | 3554 |
| CAM_and_NEM | 253 |
| CAM_only | 3094 |
| NEM_only | 80 |
| neither | 127 |
| no_CAM_site | 207 |
| share_no_CAM_site | 0.05824 |

**Position-permutation control** (each protein's CAM sites re-placed at random among its own cysteines):

| key | value |
|---|---|
| proteins_permuted | 1972 |
| observed_cysteines_in_set | 3554 |
| observed_share | 0.94176 |
| observed_wilson_low | 0.93357 |
| permuted_mean | 0.33847 |
| permuted_min | 0.32302 |
| permuted_max | 0.35537 |
| permuted_q025 | 0.32695 |
| permuted_q975 | 0.34976 |
| permutations | 200 |
| seed | 20260921 |
| reading_rule | pass if the permuted 97.5th percentile < observed Wilson lower bound |
| passes | true |

**Sequence source.** every identified peptide of the protein must match the recovered sequence at its stated Start position; first satisfying version, newest to oldest, is kept. Stored values:

| key | value |
|---|---|
| proteins_identified | 2024 |
| resolved_in_current_release | 595 |
| retired_needing_recovery | 1429 |
| recovered | 1429 |
| dropped | 0 |
| share_resolved | 1.0 |

**Detectability calibration** (Q3; trypsin, up to two missed cleavages, proline-blocked):

| rule | proteins | predicted_observable | observed | both | sensitivity | ppv |
|---|---|---|---|---|---|---|
| window_7_30 | 2024 | 13254 | 3554 | 3523 | 0.99128 | 0.26581 |
| window_6_35 | 2024 | 14110 | 3554 | 3543 | 0.9969 | 0.2511 |

**Multi-cysteine peptides against an in-silico digest** (Q4, descriptive only):

| n_cys | observed_peptides | expected_peptides | observed_share | expected_share | ratio |
|---|---|---|---|---|---|
| 1 | 2996 | 34239 | 0.89513 | 0.83163 | 1.0764 |
| 2 | 288 | 5731 | 0.08605 | 0.1392 | 0.6182 |
| 3 | 58 | 936 | 0.01733 | 0.02273 | 0.7622 |
| 4 | 4 | 213 | 0.0012 | 0.00517 | 0.231 |

**Modification-state composition of multi-cysteine peptides** (Q5, a prediction registered before the analysis):

| key | value |
|---|---|
| multi_cys_forms | 652 |
| all_CAM | 241 |
| mixed_CAM_and_NEM | 401 |
| share_all_CAM | 0.36963 |
| share_mixed | 0.61503 |
| neither_assignment | 10 |
| prediction | mixed CAM+NEM multi-cysteine peptides EXIST in non-trivial number under ABE |
| prediction_holds | true |
| note | the Modifications column of evidence.txt was glanced at before this script was written; the prediction itself was registered first |

**Audit.** The analysis audit for this note, including its correction log, is kept with the analysis records and is available on request.

**Scope.**
- Q4/Q5 describe the multi-cysteine composition of this ABE workflow; they are descriptive, not a replication

**The same statistic on published cohorts** (post hoc, descriptive; added after the analysis above). The share of identified cysteines carrying a site is computed from the positive and observed-unmodified sets already ingested for the detectability census of Figure 5a, as unique (accession, site) pairs, observed-unmodified minus positives, with Wilson 95% intervals. It is the Artifact 4 statistic only when the observed-unmodified class is the identified cysteinome of the same readout that produced the site calls; each cohort is classified by that rule from its own ingest definition, and only `comparable` and `comparable_threshold` rows are quoted in the main text.

| label | comparability | observed_class_definition | n_positive | n_observed_unmodified | share | wilson_low | wilson_high |
|---|---|---|---|---|---|---|---|
| PAT-switch (Arabidopsis) | not_comparable | observed class limited to peptides carrying a modification (subset by construction) | 5736 | 589 | 0.9069 | 0.8995 | 0.9138 |
| Cys-BOOST HeLa (human) | not_comparable | observed class from a different supplementary table than the SNO calls | 8184 | 15457 | 0.3462 | 0.3401 | 0.3523 |
| Cys-BOOST SH-SY5Y (human) | not_comparable | observed class from a different supplementary table than the SNO calls | 2134 | 403 | 0.8412 | 0.8264 | 0.8549 |
| YAP1C reporter (Arabidopsis) | comparable | positives cross-linked peptides, observed regular peptides of the same runs | 1745 | 411 | 0.8094 | 0.7922 | 0.8254 |
| FAT-switch (Arabidopsis) | not_computable | capture only: no observed-unmodified class reported | 1985 | 0 |  |  |  |
| qPerS-SID tier B (human) | comparable_threshold | same elution readout; positives by the authors' ratio threshold >= 1.30 | 724 | 1294 | 0.3588 | 0.3381 | 0.3799 |
| QTRP S1 pH5 (human) | comparable | positives -SSH sheet, observed -SH sheet of the same arm | 277 | 1158 | 0.193 | 0.1734 | 0.2143 |
| QTRP S2 pH5 (human) | comparable | positives -SSH sheet, observed -SH sheet of the same arm | 80 | 1105 | 0.0675 | 0.0546 | 0.0832 |
