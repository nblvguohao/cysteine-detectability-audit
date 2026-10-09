# Release notes - v3.2.0 (8 October 2026)

Accompanies the Molecular & Cellular Proteomics submission (manuscript package v3.2, title *Five detectability
artifacts decide what site-level comparisons can show in cysteine-modification proteomics*). Changes from
v3.1.2 (commit 6517b23). Cys-Audit is unchanged (0.2.2; 55 tests pass). No analysis was re-run and no stored
output changed.

## 1. Folders renamed

| before | after |
|---|---|
| `revision_2026-09-30/` | `posthoc_2026-09-30/` |
| `revision_2026-10-02/` | `posthoc_2026-10-02/` |

The files inside are unchanged byte for byte, so every SHA-256 recorded in their `provenance.json` and audit
files still holds. Paths recorded inside those files, and the working-folder constants of their scripts, keep
the former names (see Scope in `README.md`). The two folders' contents are the post hoc analyses of
Supplemental Notes 7 and 14-22, labeled as such in the manuscript and in Supplemental Data 12.

## 2. Supplemental material updated to the v3.2 package

- Added `supplemental/Supplemental_Methods.md`: the full experimental procedures, moved out of the main text,
  with reference numbers matching the v3.2 reference list.
- Supplemental Notes 2, 5, 6, 7, 11 and 14-22 and Supplemental Data 12: status lines and provenance wording
  now read "post hoc, added after the primary analyses were complete"; folder paths updated to the new names;
  reference numbers in Notes 6 and 17 updated to the v3.2 reference list; minor spelling (US English) and one
  added file reference in Note 6, as in the submitted package.
- `source_data_submitted/`: Source Data for Figures 5 and 6 updated in the same way (provenance text only;
  no value changed).

## 3. Documentation

`README.md`, `posthoc_2026-09-30/README.md`, `figures_r31/README.md`, `archive/README.md` and `CITATION.cff`
(version v3.2.0) updated for the new folder names and title. `MANIFEST_SHA256.txt` regenerated.
