# cysteine-detectability-audit

Code and result tables for the manuscript *Comparison Backgrounds Shape Site-Level Inference in Cysteine Proteomics* (Guohao Lv and Lichuan Gu). The article reference will be added on publication. The manuscript and its Supporting Information are not part of this repository.

## Contents

| Path | Contents |
| --- | --- |
| `ac_minimal_reproducibility/` | Current scoring reconstruction: feature construction (VIS10, CPL15, DIG25), fixed-fold model fitting for the two training backgrounds (NEG_A, NEG_B), held-out scoring and the paired unit-bootstrap. Code and pinned dependencies only; its README gives the input format. |
| `cys-audit/` | Cys-Audit command-line tool, version 0.2.2, used for the supplementary diagnostic checks. Earlier archived versions of the tool are not included. |
| `results/source_data/` | Figure-level source tables and result tables for the current manuscript: comparison-support counts, fixed-scorer contrasts, the TMT-ABPP application test, protease and flank analyses, the external ACR/IAM case, simulations, the assay-label manifest and the claim ledger. |

## Where each displayed result is

Several result files keep the internal names under which they were generated. This table maps the current displays to their files (paths relative to `results/source_data/`).

| Display | Files |
| --- | --- |
| Figure 1c | `SourceData_Fig1c.csv` |
| Figure 2a | `Grounding_2026-10-21/SourceData_Fig5_support.csv` |
| Figure 2b,c | `Grounding_2026-10-21/SourceData_Fig5_effects.csv` |
| Figure 3a | `SourceData_Fig2a.csv` |
| Figure 3b | `SourceData_Fig2b.csv` |
| Figure 3c | `SourceData_Fig2c.csv` |
| Figure 4a–d | `SourceData_Fig3a.csv` to `SourceData_Fig3d.csv` |
| Figure S1a–d | `SourceData_Fig5a.csv` to `SourceData_Fig5d.csv` |
| Figures S2–S4 | `Grounding_2026-10-21/SourceData_FigS2.csv` to `SourceData_FigS4.csv` |
| Figure S5 | `SourceData_FigS5.csv` |
| Figure S6a–c | `SourceData_FigS6a.csv` to `SourceData_FigS6c.csv` |
| Figure S7 | `SourceData_Fig4a.csv`, `SourceData_Fig4bcd.csv` |
| Table S14.2 | `CXXC_sparse_counts.tsv` |

Column definitions are in `FIGURE_DATA_DICTIONARY.md`, and `DISPLAY_CROSSWALK_2026-10-26.md` records how the internal names relate to the current numbering. The other files in `Grounding_2026-10-21/` hold the underlying external-case and assay-composition results.

## Running the code

Current scoring reconstruction (Python 3.10 recommended):

```sh
cd ac_minimal_reproducibility
python3.10 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python run.py --sites sites.tsv --fasta reference.fasta --folds folds.tsv.gz --output new_run
.venv/bin/python -m unittest discover -s tests
```

The runner needs separately supplied, version-matched inputs: prepared site membership, full protein sequences and fixed folds. Their format is described in `ac_minimal_reproducibility/README.md`. These prepared inputs, fitted models and held-out predictions are not distributed here, so this repository alone does not regenerate the article's numerical values.

Cys-Audit:

```sh
cd cys-audit
pip install -e .
cys-audit schema
python -m pytest tests
```

## Data sources

Public deposits and published site tables are identified by accession and DOI in the manuscript's Supporting Information. Supplementary tables from other groups' publications are not redistributed here; they are available from their publishers.

## License

Code is released under the MIT License (`LICENSE`). Result tables under `results/` are released under CC BY 4.0 (`LICENSE-DATA`). Third-party materials keep their original terms.
