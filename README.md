# Five detectability artefacts in cysteine-modification proteomics - code and derived data

This repository accompanies the manuscript *Five detectability artefacts that inflate, reverse, or
preclude site-level conclusions in cysteine-modification proteomics*. It contains the audit tool, the
pre-registration protocols for every registered analysis, the analysis outputs behind the figures and
tables, and the supplementary material exactly as the manuscript cites it.

## Contents

| path | what it is |
|---|---|
| `cys-audit/` | the Cys-Audit command-line tool (version 0.2.2), with its 55-test suite under `cys-audit/tests/` |
| `protocols/` | the pre-registration protocols, including those named in Supplemental Notes 3-7 |
| `scripts/` | the analysis code that produced the released outputs |
| `results/` | the stored analysis outputs named in the Source Data |
| `supplemental/` | the Supplemental Notes 1-22 and Supplemental Data 1-12, numbered as the manuscript cites them |
| `figures_submitted/` | the figures and the graphical abstract as submitted with the revised manuscript |
| `revision_2026-09-30/` | the analyses added in revision (Supplemental Notes 14-22) with their scripts, inputs and outputs, and the build scripts of the revised Figures 1, 3, 5 and 6; see its README |
| `source_data_r31/`, `figures_r31/` | the released figure assets and their source data |
| `recoding_audit/` | the blinded re-coding audit behind Supplemental Data 11 |
| `analysis/` | the non-null check of the normal approximation behind the binary-feature statistic (Supplemental Note 13) |
| `reports/` | the analytical reports corresponding to the Supplemental Notes |
- `source_data_submitted/` - the Source Data exactly as submitted with the manuscript, which is the
  serialisation the submitted figures were drawn from. `source_data_r31/` is the earlier working
  serialisation used by the phase scripts; where the two differ, `source_data_submitted/` is
  authoritative and is what the figures and tables cite.
| `MANIFEST_SHA256.txt` | SHA-256 of every file in this repository |

## Reuse

Code is released under the MIT licence (`LICENSE`); derived tables are released under CC BY 4.0
(`LICENSE-DATA`). The tool is installed with `pip install -e cys-audit`; the regression suite runs with
`pytest cys-audit/tests` and must pass before the tool is applied to anything.

Every analysis in the manuscript states whether it was registered and which decisions were taken after
a result had been seen; `supplemental/Supplemental_Data_12_analysis_provenance.csv` collects that record
in one place.

## Scope

This repository carries the artefacts the manuscript cites. It is not the authors' working tree: internal
project records, drafting scripts and material from unrelated projects are not part of it. For anything not
found here, contact the corresponding author.
