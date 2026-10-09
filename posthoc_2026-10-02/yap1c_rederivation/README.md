# YAP1C re-derivation (Supplemental Note 7)

Script and output behind the paragraph "Re-derivation check for the YAP1C cohort" of Supplemental Note 7
(post hoc, added on 2 October 2026). The census value in the Note 7 table (0.8094) comes from the census
ingest; this folder recomputes the coincidence share from the search reports deposited with the study
(YAP1C disulfide-linked peptide reporter, Arabidopsis; Front. Plant Sci. 2020, doi 10.3389/fpls.2020.00777,
PMID 32714340; PRIDE PXD016723).

| file | what it is |
|---|---|
| `compute_fps2020_coincidence.py` | the computation (Python 3; needs `openpyxl`) |
| `coincidence_output.txt` | its standard output from the run of 2 October 2026 |

## Values used in Supplemental Note 7

From `coincidence_output.txt`, block `B_proteinLevel_S_runs` (pLink regular-peptide reports of the
protein-level IgG-enrichment runs, `_S` runs) and the final `B` lines (Wilson 95% intervals):

- every matching protein counted: coincidence 0.8469 [0.8333, 0.8596] -> **0.847 [0.833, 0.860]**
- first listed protein only: 0.8304 [0.8138, 0.8458] -> **0.830 [0.814, 0.846]**
- cysteines observed free (carbamidomethylated) in the same runs that are also sites: 714/1,153 =
  0.6193 [0.5909, 0.6468] -> **0.619 [0.591, 0.647]**

The other blocks (peptide-level runs, the non-enriched test samples of Dataset S1b) are printed for
context and are not cited in the manuscript.

## Inputs (not redistributed)

The inputs are third-party files and are not included. Download them and place them in one folder
(default `inputs/` next to the script, or set `YAP1C_INPUT_DIR`, or pass `--input-dir`; `--fasta`,
`--si-dir` and `--plink-dir` override single locations). The pLink zips can be given as zips or as their
extracted folders of the same names.

| file | source | bytes | SHA-256 |
|---|---|---:|---|
| `pLink_reports_untreated_proteinLevel.zip` | PRIDE PXD016723 (`ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2020/12/PXD016723/`) | 5,590,959 | `15bdbc22e813e53678601340dc7a4019514d9d09297494b335da0e2936bfffc0` |
| `pLink_reports_H2O2_proteinLevel.zip` | PRIDE PXD016723 | 4,648,005 | `4e91624a7587c3b9e34f1a660faadb8d2ef27920bca432aa213dbc9ff462435d` |
| `pLink_reports_untreated_peptideLevel.zip` | PRIDE PXD016723 | 1,359,634 | `c99024a91cde9423b129743d3e0280cfa83faa51d020949a2b65e1443dab1d84` |
| `pLink_reports_H2O2_peptideLevel.zip` | PRIDE PXD016723 | 1,008,984 | `3b296c1f8648b0b73f0a953abc9e6227e309af07d08c35c073eeb4ec24b2a428` |
| `Data_2.xlsx` (sheet `DatasetS2`, the site list) | supplementary data of the study (Frontiers) | 6,038,709 | `7ed4eb91ce0611f2e3251fce290a0b91d18bc4e5ac4ac75150a59af4716daa03` |
| `Data_1.XLSX` (sheet `DatasetS1b_regular_peptide`; used only for the S1b block) | supplementary data of the study (Frontiers) | 3,465,965 | `2c01b2b51e7dddd5f33ae18a29cb1780fa5e8cb1f9a56b841397ec83aec4e2b0` |
| `Ath_TAIR10_pep.fa.gz` | Ensembl Plants, Arabidopsis thaliana TAIR10 `pep.all` (downloaded 2 October 2026) | 9,698,724 | `084775e9ad341b0a133ff2d51c68f19a697f2d4a325cd8751ef548c418a3643d` |

pLink searched Araport11 representative models; 27 peptide-protein matches do not map to the TAIR10
sequences and are skipped (reported as "unmapped protein-matches").

## Re-running

```bash
python compute_fps2020_coincidence.py --input-dir /path/to/inputs > out.txt
cmp out.txt coincidence_output.txt
```

Checked on 2 October 2026 (Python 3.9.6, openpyxl 3.1.5) from the four zips alone and from the extracted
folders: standard output identical byte for byte to `coincidence_output.txt` (openpyxl prints a
conditional-formatting warning on standard error). The script differs from the run of 2 October 2026 only
in how input locations are set (the original hard-coded a temporary working folder) and in having the
Wilson-interval helper inlined; the computation is unchanged.
