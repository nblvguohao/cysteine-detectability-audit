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
