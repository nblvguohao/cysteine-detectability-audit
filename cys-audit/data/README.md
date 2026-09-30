# Data

No data are shipped here.

## Input schema (one row per cysteine)

| column | required | type | meaning |
|---|---|---|---|
| protein | yes | str | protein accession; must match a FASTA identifier if `--fasta` is given |
| position | yes | int | 1-based residue index of the cysteine |
| label | yes | 0/1 | 1 = reported modified site; 0 = candidate background cysteine |
| detected | no | 0/1 | 1 = the cysteine was observed in the run in any form |
| n_cys_peptide | no | int | smallest number of cysteines on any peptide that observed this cysteine |
| abundance | no | float | protein abundance (> 0), identical on all rows of one protein |
| peptide_mass | no | float | monoisotopic mass (Da) of the peptide supporting a positive row |
| cluster | no | str | resampling cluster (e.g. homology component); default = protein |

Any extra numeric column can be named as a claim feature (`--claim-feature`).
`cys-audit schema` prints this table.

## Public demonstration dataset

PXD063463 (PRIDE; Kim et al., *Proteomics* 2026, doi 10.1002/pmic.70100). Files used:
`evidence_ABE.txt`, `peptides_ABE.txt`, `proteinGroups_Global.txt`, `Mouse_2025_UP000000589_10090.fasta`.
Put them in `data/PXD063463/` and run `make reproduce`.

## Non-redistributable data

The tomato persulfidation site table used in the authors' earlier analyses was produced by a commercial
service for a collaborating laboratory and may not be redistributed. The tool does not need it: every
test runs on the schema above, and `cys-audit simulate` produces synthetic substitutes with a planted
artefact of known strength.
