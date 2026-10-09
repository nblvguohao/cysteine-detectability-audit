# Supplemental Note 4. The four-protease public deposit PXD063463

PXD063463 is an acyl-biotin exchange study in which one enrichment was split across four proteases. It is the only proteome-scale deposit our PRIDE screen found in which the same enrichment was digested more than one way, and it carries main Fig. 3c, d, Fig. 4a and Fig. S1.

**Input construction.** Each arm was converted from the deposit's MaxQuant output by `cys-audit/examples/convert_maxquant_abe.py`. A cysteine is *identified* when it lies in a peptide with at least one evidence row in that arm, and is a *site* when it carries carbamidomethyl with a localization probability of at least 0.75, which in this chemistry is the readout of the exchange step (free thiols are blocked with N-ethylmaleimide before the exchange). Proteins are leading razor proteins; reverse hits and contaminants are removed.

**The arms without hydroxylamine are a specificity check, and they are not empty.** Every arm treated with hydroxylamine has more sites than its matched arm without it, which is the registered sanity condition, but the arms without hydroxylamine still carry carbamidomethylated cysteines, so the positive sets are not free of non-specific capture.

| arm | identified cysteines (+hydroxylamine) | sites (+hydroxylamine) | identified (no hydroxylamine) | sites (no hydroxylamine) |
|---|---|---|---|---|
| trypsin | 5328 | 1359 | 2024 | 660 |
| AspN | 187 | 115 | 22 | 3 |
| chymotrypsin | 728 | 411 | 61 | 23 |
| GluC | 884 | 519 | 68 | 34 |

**Cleavage geometry.** The statistic is the Haldane log2 odds ratio of carrying a cleavage-competent residue within a band of the cysteine, positives against background, with protein-clustered percentile bootstrap intervals over 5000 replicates at the Bonferroni level for the two bands, recorded as `ci_level` 0.975 in each stored audit. Two bands are tested, 1 to 3 residues (proximal) and 6 to 12 residues (distal); the registered comparison between rules uses the band with the larger absolute estimate, so both bands are reported here and in Fig. S1a.

| arm | rule | background | proximal 1-3 | distal 6-12 |
|---|---|---|---|---|
| trypsin | trypsin | proteome | -0.108607 [-0.285308, 0.068814] | 0.972971 [0.711516, 1.275616] |
| AspN | aspn | proteome | 1.110573 [0.420214, 1.787266] | 1.539092 [0.793042, 2.492719] |
| chymotrypsin | chymotrypsin | proteome | 0.50446 [0.174646, 0.849023] | 1.194606 [0.801952, 1.63019] |
| GluC | gluc | proteome | 0.523234 [0.21904, 0.823891] | 0.944881 [0.573888, 1.347499] |
| AspN | trypsin | proteome | 0.109145 [-0.55249, 0.700185] | 0.118877 [-0.554927, 0.863666] |
| chymotrypsin | trypsin | proteome | 0.401074 [0.057153, 0.739983] | 0.334937 [-0.065642, 0.775046] |
| GluC | trypsin | proteome | -0.008047 [-0.351532, 0.332582] | 0.383989 [-0.064028, 0.846555] |
| trypsin | trypsin | observed | -0.694908 [-0.8951, -0.494816] | 0.060834 [-0.235279, 0.380768] |
| AspN | aspn | observed | -0.603101 [-1.563753, 0.278747] | 0.009696 [-1.346075, 1.279257] |
| chymotrypsin | chymotrypsin | observed | -0.898437 [-1.40473, -0.410929] | -0.182265 [-0.847299, 0.430168] |
| GluC | gluc | observed | -1.776797 [-2.297182, -1.293036] | -0.550606 [-1.186751, 0.04365] |

**Backgrounds.** *Proteome* is every other cysteine of the identified proteins, whether or not it was seen. *Observed* is the identified cysteines without a site, which is the within-experiment negative class. The distal signature is present against the proteome background in every arm under its own rule and absent against the observed background in all four arms, which is what a detectability explanation predicts. A proximal depletion remains against the observed background in three arms; it is reported here without interpretation.

**Coincidence.** The share of identified cysteines that carry a site, with a protein-clustered 95 per cent interval (`results/phase4_public_demo_2026-09-22.csv`).

| arm | coincidence | 95 per cent interval |
|---|---|---|
| trypsin | 0.255068 | 0.238042 to 0.272426 |
| AspN | 0.614973 | 0.540533 to 0.685992 |
| chymotrypsin | 0.56456 | 0.522788 to 0.607241 |
| GluC | 0.587104 | 0.547189 to 0.625432 |

**Limits.** Three of the four arms are small, so their intervals are wide and the comparison between rules rests mostly on the trypsin arm. The modification is S-palmitoylation in mouse macrophages, a different chemistry and organism from the redox modifications the rest of the paper is about. The site readout is carbamidomethyl, so the positive set inherits whatever non-specific capture the exchange step leaves, as the arms without hydroxylamine show.

