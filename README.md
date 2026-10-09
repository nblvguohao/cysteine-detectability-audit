# Five detectability artefacts in cysteine-modification proteomics - code and derived data

This repository accompanies the manuscript *Five detectability artifacts decide what site-level
comparisons can show in cysteine-modification proteomics*. It contains the audit tool, the
pre-registration protocols for every registered analysis, the analysis outputs behind the figures and
tables, and the supplementary material exactly as the manuscript cites it.

## Contents

| path | what it is |
|---|---|
| `cys-audit/` | the Cys-Audit command-line tool (version 0.2.2), with its 55-test suite under `cys-audit/tests/` |
| `protocols/` | the pre-registration protocols, including those named in Supplemental Notes 3-7 |
| `scripts/` | the analysis code that produced the released outputs, with the local modules it imports (see Scope) |
| `results/` | the stored analysis outputs named in the Source Data, including the score tables of Figure 7 of v3.2 (Supplemental Figure S3 of v3.3) |
| `artifacts12/` | the scripts and stored outputs behind Artifacts 1 and 2, deposited unchanged with their SHA-256; see its README for what re-runs and the four disclosed gaps |
| `supplemental/` | the Supplemental Methods, Supplemental Results, Supplemental Figure legends, Supplemental Notes 1-22 and Supplemental Data 1-12, numbered as manuscript v3.3 cites them |
| `source_data_final/` | the Source Data submitted with manuscript v3.3, for Figures 1-6, Supplemental Figures S1-S7, the graphical abstract and the text-cited tables; see "Figures and Source Data" below |
| `figures_final/` | the figures of manuscript v3.3 as submitted (Figures 1-6, Supplemental Figures S1-S7 and the graphical abstract; PDF and PNG) |
| `posthoc_2026-10-09_figures/` | the build scripts of `figures_final/` and `source_data_final/`, with the reproduction and value checks; see its README |
| `source_data_submitted/` | the Source Data of the v3.2 submission (Figures 1-7, v3.2 numbering), which is the serialisation the v3.2 figures were drawn from; it is authoritative wherever it differs from `source_data_r31/`, and the v3.3 Source Data cite its tables |
| `figures_submitted/` | the figures and the graphical abstract of the v3.2 submission (v3.2 numbering); the v3.3 build re-uses some of them unchanged |
| `posthoc_2026-09-30/` | the post hoc analyses (Supplemental Notes 14-22) with their scripts, inputs and outputs, and the build scripts of Figures 1, 3, 5 and 6 of the v3.2 submission; see its README |
| `posthoc_2026-10-02/` | the YAP1C re-derivation of Supplemental Note 7 (script and output; third-party inputs listed with SHA-256, not redistributed); see `yap1c_rederivation/README.md` |
| `source_data_r31/`, `figures_r31/` | the figure set and source data of the earlier R31 working version, kept because the R31 build and plotting scripts in `scripts/` read or write them; superseded by the two `*_submitted/` folders (see `figures_r31/README.md`) |
| `recoding_audit/` | the blinded re-coding audit behind Supplemental Data 11 |
| `analysis/` | the non-null check of the normal approximation behind the binary-feature statistic (Supplemental Note 13) |
| `reports/` | the analytical working reports corresponding to the Supplemental Notes (in Chinese) |
| `archive/` | four internal baseline-audit scripts and an unused prospective blind-test plan, kept for the record; outside the reproducibility claim (see its README) |
| `MANIFEST_SHA256.txt` | SHA-256 of every file in this repository |
| `RELEASE_NOTES_v3.3.0.md` | what changed from v3.2.0 (v3.3 supplemental material, figures, Source Data and figure build) |
| `RELEASE_NOTES_v3.2.0.md` | what changed from v3.1.2 (folder renames, Supplemental Methods, wording of the Supplemental Notes) |
| `RELEASE_NOTES_v3.1.2.md` | what changed from v3.1.1 (files cited by the Supplemental Notes, with their recorded SHA-256) |
| `RELEASE_NOTES_v3.1.1.md` | what changed from v3.1.0, including every file removed and why |

## Reuse

Code is released under the MIT licence (`LICENSE`, also `cys-audit/LICENSE`); derived tables are
released under CC BY 4.0 (`LICENSE-DATA`). Cite with `CITATION.cff` (release v3.3.0). The tool is installed with `pip install -e cys-audit`; the regression suite runs with
`pytest cys-audit/tests` and must pass before the tool is applied to anything.

Every analysis in the manuscript states whether it was registered and which decisions were taken after
a result had been seen; `supplemental/Supplemental_Data_12_analysis_provenance.csv` collects that record
in one place.

## Figures and Source Data

Manuscript v3.3 renumbered and redrew the figures. `figures_final/` and `source_data_final/` are the set
submitted with v3.3; `figures_submitted/` and `source_data_submitted/` are the set of the v3.2 submission, kept
unchanged because the v3.3 Source Data name their tables as sources
(`cysteine-detectability-audit/source_data_submitted/...`) and because the v3.3 build reads them. No value was
re-estimated for v3.3: every value in `source_data_final/` is read from a file of this repository. Supplemental
Figures S1, S2, S3 and S6 are Figures 2, 3, 7 and 4 of v3.2, copied unchanged; Supplemental Figure S5 places
panels b and c of v3.2 Figure 5 unchanged (letters and titles redrawn). The five `Source_Data_text_*` files are
identical in both folders. `posthoc_2026-10-09_figures/README.md` names the build script of each v3.3 figure and
its inputs, and records that the build reproduces `figures_final/` and `source_data_final/`.

## Scope

This repository carries the artefacts the manuscript cites. It is not the authors' working tree: internal
project records, drafting scripts and material from unrelated projects were removed in v3.1.1
(`RELEASE_NOTES_v3.1.1.md` lists each file). What can be re-run from it, and what cannot:

- **Re-runnable from the repository** (public inputs fetched as documented): the Cys-Audit tool and its
  tests (`cys-audit/`); the Supplemental Note 13 check (`analysis/`); the post hoc analyses
  (`posthoc_2026-09-30/`, inputs included, external downloads listed with SHA-256 in its `DOWNLOADS.md`);
  Artifacts 1 and 2 from the deposited site tables (`artifacts12/`, scripts 10, 12, 24, 36 and 21, with the
  exceptions its README states); the re-coding audit scripts (`recoding_audit/`); the v3.3 figures and their
  Source Data (`posthoc_2026-10-09_figures/`, macOS only because it embeds the system Helvetica).
- **`scripts/`** is the code as it was run in the authors' analysis tree. Every local module these scripts
  import is now included, so they import; but most of them read that tree's `inputs/` and `external/`
  folders (feature matrices, per-cysteine model scores, PRIDE downloads), which are not
  redistributed, so they document the computation rather than re-run it here. Functions that call
  `v2_apply.load_source_encoder()` (in `run_qtrp_trained_model.py`,
  `run_detection_matched_persulfidation_test.py`, `run_qtrp_intervals_and_cross_label.py`) and
  `run_retrospective_external_corrected.py` also need a feature encoder from an earlier project
  (`model_pipeline.py`) that is not released. `ingest_qpers_sid_supplementary.py` needs Python 3.10 or later.
- **Absolute paths.** Many scripts under `posthoc_2026-09-30/scripts/` and a few under `scripts/` set
  their working folders with machine-specific paths, written with a `/path/to/` placeholder (for example
  `W = "/path/to/.../revision_2026-09-30"`, `REPO = ".../_cys_repo_work/repo"`); set those constants
  to your checkout (`W` -> `posthoc_2026-09-30/`, `REPO` -> the repository root, `MCP` -> the repository
  root with `source_data_submitted/` for `source_data/`). `RELEASE_NOTES_v3.1.1.md` lists the files.
- **`archive/`** holds scripts that cannot run outside the authors' machines; they are kept for the record.
- **`archive/protocols_not_used/`** holds the plan, templates and prediction lock of a prospective blind test
  of the authors' own ranker (`BLIND_TEST_PROTOCOL*.md`, `blind_*template*.csv`,
  `pka_blind_prediction_lock_2026-09-17.json`; moved from `protocols/` in v3.1.2). That test is not part of
  the manuscript, and nothing in the manuscript, the Supplemental Notes and Data or the Source Data relies on
  these files.
- **Internal shorthand in hashed protocols.** Some protocol JSON files in `protocols/` (and frozen copies under
  `cys-audit/archive/`) contain free-text fields with the authors' internal shorthand, for example references to
  numbered "house rules" or to sections of a project notes file (`notes s. ...`). These refer to internal
  working notes that are not part of the release. The registered analysis is fully stated in their structured fields and in the
  Supplemental Notes.
- **Line endings.** Text files are committed with LF, except files whose SHA-256 was recorded on CRLF bytes
  (three CSV files in `artifacts12/` and eight files in `results/` added in v3.1.2); these are committed byte for
  byte so the recorded values hold (see `RELEASE_NOTES_v3.1.2.md`).
- **Folder names.** `posthoc_2026-09-30/` and `posthoc_2026-10-02/` were named `revision_2026-09-30/` and `revision_2026-10-02/` before release v3.2.0. The files inside are unchanged byte for byte, so the SHA-256 values recorded in their `provenance.json` and audit files still hold; paths recorded inside those files and the working-folder constants of their scripts use the former names.

For anything not found here, contact the corresponding author.
