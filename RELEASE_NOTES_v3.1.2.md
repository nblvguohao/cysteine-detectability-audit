# Release notes - v3.1.2 (2 October 2026)

Accompanies the Molecular & Cellular Proteomics submission (manuscript v3.1.1 package).
Changes from v3.1.1 (commit 5b25963). Cys-Audit is unchanged apart from one code comment (0.2.2; 55 tests
pass).

## 1. Added: files cited by Supplemental Notes 3, 6 and 19 that v3.1.1 lacked

Copied from the authors' analysis tree. For every file a SHA-256 is recorded in the Supplemental Notes or in
a provenance/audit file of this repository, and **every committed file matches its recorded SHA-256**.
These values were computed on the CRLF bytes the files had when they were recorded (the analysis ran on
Windows), so the files are committed with CRLF line endings, byte for byte; `.gitattributes` (`* -text`)
keeps those bytes on checkout. The LF-normalised SHA-256 is given for comparison (for the four audit JSON
files and `pxd015307_research_psms.csv`, the provenance files record it as `sha256_lf_normalised`).

| file | cited in | SHA-256 recorded (where) | SHA-256 committed (CRLF) | SHA-256 if LF-normalised |
|---|---|---|---|---|
| `results/phase2_synthetic_benchmark_summary_2026-09-22.csv` | Note 3 | `00407aec…8030` (`outputs` of the benchmark run's audit `phase2_synthetic_benchmark_2026-09-22_audit.json`, analysis tree; that audit is not in this repository) | `00407aec02810512ef09c5c695a279161b0666db2840452c940fc014502f8030` (match) | `7986b3aa52c896b3f41e341d98901af6e9f123e98b2e7999314c4a7466fe4a8c` |
| `results/phase2_calibration_precision_2026-09-22.csv` | Note 3 | `c0789db0…aaa71` (`revision_2026-09-30/results/E_bootstrap_coverage/provenance.json`) | `c0789db0b04d87ff038e34a68b64592621079f4a40695d94ad945e220a0aaa71` (match) | `c23f7b37d18f00b909e7b360d2baa80cb9823da74f085f8b29e882d45f78b517` |
| `results/pride_coincidence_candidates_2026-09-19.csv` | Note 19 | `4f135b86…` (Note 19) | `4f135b860b20439d34ee1c83d51ebe8c622e897f2f6bed110ec9ece74001722f` (match) | `c2bfc738e30d045f3da5b2986fe1784e7cb1492d94c7228270d439f39e906913` |
| `results/public_refit_inputs_2026-09-21_audit.json` | Note 19 | `aa4186e8…` (Note 19; full value in `revision_2026-09-30/results/J_provenance_species/provenance.json`) | `aa4186e814248ad4d81df2f0d08e53c856caae9ef3e4739916810c5bfb8ff61f` (match) | `0cf6723f3ef999d70afc67e497f920420cd533ee66c9ff9ffc2db6d1f75b6a5b` |
| `results/self_audit_public_cohorts_2026-09-20_audit.json` | Note 19 | `976b11fb…` (Note 19; full value in `revision_2026-09-30/results/F_rice_artifact3/provenance.json` and `J_provenance_species/provenance.json`) | `976b11fb0c668f2e4a1d264829779fca007d91cd808d3667e0ceb276092039a4` (match) | `b3236017d5d2372303a514ad13564db1fb155dc8c5b371eb1d6a286a4afe603d` |
| `results/pxd015307_research_psms.csv` | Note 6 | `8369c070…22ce` (`revision_2026-09-30/results/H_artifact5_docs/provenance.json`) | `8369c070c4cc8124030588e7f6cb77518c82be94603245a2e240a136e5ef22ce` (match) | `b94c946de56d6bd2d8c67cc675ba039d35178378e26138d44b38075d2e1ce602` |
| `results/pxd015307_research_audit.json` | Note 6 | `a85c8f61…bec6` (`H_artifact5_docs/provenance.json`) | `a85c8f61afd85945b396e0780f6d1daffc8a303d1df3336be0c00f9cea93bec6` (match) | `d9086c9d7110be38019547c55cb2fc933c8b82a8b63b076aeb0b712d1bcbeb37` |
| `results/pxd015307_posthoc_score_ties_audit.json` | Note 6 ("and its audit") | `42400e14…579d` (`H_artifact5_docs/provenance.json`) | `42400e14de2cf337e1acbc108029232a052dd50f56f19d2e08b00bffcd9a579d` (match) | `e1f4f13aa5cad6016f2118f14695fd17940a6ec503dd464e6440ca95247ec497` |

Notes:

- `pxd015307_research_psms.csv` holds only its header (`arm`): no arm of the PXD015307 re-search has a
  1%-FDR PSM with a +32-class cysteine modification (`pxd015307_research_audit.json`, `n_plus32_cys_psms`).
- The LF-normalised copies of `pxd015307_research_psms.csv` and `pxd015307_research_audit.json` under
  `revision_2026-09-30/inputs/repo_results/` are unchanged (they are the inputs that folder's scripts read,
  and provenance records both hashes).
- `pxd015307_research_audit.json` names the remote host directory on which the Comet search ran
  (`search_ran_on`). It is kept because the file's SHA-256 is recorded. None of the eight files contains a
  path on a personal computer.
- The scripts of revision F (`revision_2026-09-30/scripts/F_rice_artifact3/`) that read
  `results/self_audit_public_cohorts_2026-09-20_audit.json` now find it in the repository.

## 2. Added: `revision_2026-10-02/yap1c_rederivation/`

Script and output behind the YAP1C re-derivation paragraph of Supplemental Note 7 (0.847 [0.833, 0.860];
0.830 [0.814, 0.846]; 714/1,153 = 0.619 [0.591, 0.647]). The inputs (four pLink report zips from PRIDE
PXD016723, two supplementary data sheets of the study, the TAIR10 protein FASTA) are not redistributed;
the README lists them with sizes and SHA-256. The script now takes its input locations from
`--input-dir` / `YAP1C_INPUT_DIR` (default `inputs/` next to it) instead of a temporary folder, and can read
the zips directly. Re-run on 2 October 2026 from the original inputs: standard output identical byte for
byte to the stored `coincidence_output.txt`.

## 3. Moved

`protocols/BLIND_TEST_PROTOCOL.md`, `protocols/BLIND_TEST_PROTOCOL_V2.md`, the five
`protocols/blind_*template*.csv` files and `protocols/pka_blind_prediction_lock_2026-09-17.json` ->
`archive/protocols_not_used/` (plan, templates and prediction lock of a prospective blind test of the
authors' own ranker that is not part of the manuscript). None is cited by the manuscript, the Supplemental
Notes and Data or the Source Data (checked by name before moving). `archive/README.md` lists them and the two
scripts that still name the former paths.

## 4. Updated

- `cys-audit/src/cys_audit/constants.py`: a comment that cited an internal "house rule" now reads "standard
  procedure" (comment only; no code change; the 55 tests pass). The frozen copies under `cys-audit/archive/`
  are unchanged.
- `artifacts12/README.md`: the Known gaps list now states its count (four gaps; the former item 5 is a note,
  not a gap) and explains that `MANIFEST.json` keeps the original `source_tree` path verbatim because its
  SHA-256 was recorded at staging. Root `README.md` "two disclosed gaps" -> "four disclosed gaps".
- `README.md`: contents table (`revision_2026-10-02/`, `RELEASE_NOTES_v3.1.2.md`); Scope notes on
  `archive/protocols_not_used/`, on internal shorthand in the free-text fields of hashed protocol JSON files
  (left unchanged because their SHA-256 values were recorded at registration), and on line endings.
- `CITATION.cff`: version v3.1.2, released 2026-10-02.
- `MANIFEST_SHA256.txt` regenerated (same format as v3.1.1).
