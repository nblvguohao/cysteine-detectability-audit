# archive/ - kept for the record, not reproducible from this repository

The scripts here are outside the repository's reproducibility claim. They are kept, unchanged, because
other files cite them; they cannot be run from this repository.

| file | why it cannot run here |
|---|---|
| `scripts/baseline_phase0_script_integrity_2026-09-21.py` | compares the internal analysis tree and the collaborator tree against their audit records; reads both trees by absolute path |
| `scripts/baseline_phase0_rerun_2026-09-21.py` | re-ran producer scripts of both trees in a scratch mirror (`PHASE0_SCRATCH`) with machine-specific interpreters (absolute paths) |
| `scripts/baseline_phase0_rerun_correction_2026-09-21.py` | correction run of the above for the collaborator-tree rows (cited in `revision_2026-09-30/results/F_rice_artifact3/provenance.json` and `revision_2026-09-30/scripts/F_rice_artifact3/s09_trace_artifact2.py`, which give its former path `scripts/...`) |
| `scripts/baseline_phase0_findings_addendum_2026-09-21.py` | imports `baseline_phase0_summarise_2026-09-21.py`, which is not released because it reads the collaborator tree by absolute path |

What these runs established for Artifacts 1 and 2 (re-execution of 10, 12, 24 and 36) was repeated on
2 October 2026 from the deposited files themselves; see `artifacts12/README.md`, "Re-running".

## `protocols_not_used/`

Moved from `protocols/` in v3.1.2, unchanged. These files belong to a planned prospective blind test of the
authors' own ranker; that test is not part of the manuscript, and none of these files is cited by the
manuscript, the Supplemental Notes and Data or the Source Data.

| file | what it is |
|---|---|
| `BLIND_TEST_PROTOCOL.md`, `BLIND_TEST_PROTOCOL_V2.md` | the blind-test plan, first and second versions (in Chinese) |
| `blind_input_manifest_template.csv`, `blind_predictions_template.csv`, `blind_predictions_template_v2.csv`, `blind_protein_occurrence_template.csv`, `blind_site_outcomes_template.csv` | empty templates for that test |
| `pka_blind_prediction_lock_2026-09-17.json` | a prediction lock for a pKa round of that test |

Two files in `scripts/` still name the former paths and are left as they were run (code as run,
not edited after the fact): `run_pka_round_verdict.py` reads `protocols/pka_blind_prediction_lock_2026-09-17.json` (now
`archive/protocols_not_used/`), and the docstring of `run_detectability_stratified_reverse_test.py` cites
`protocols/BLIND_TEST_PROTOCOL_V2.md`. The `blind_analysis_plan*.json` files stay in `protocols/`.
