# Figure 6 (revision of 2026-09-30)

The 28 re-tested claims under the final tally (`../../final/final_tally.csv`): 11 survive, 2 are attenuated, 13 are undecidable, 1 null claim is broken by the control and 1 baseline runs against the claimed direction. `SFE-006` is drawn from its re-run on its authors' own data, on which it is tallied, with its transfer re-test drawn just below its row; `SNO-012` is undecidable because its baseline interval covers zero.

| file | content |
|---|---|
| `fig6_final_2026-09-30.py` | builds the figure and its Source Data by copying stored values verbatim from the released Source Data, the `SFE-006` re-run and the final tally |
| `build_log.txt` | output of the build gates |
| `verify_fig6_final.py` | reads the PDF back and compares every drawn value with the Source Data |
| `Fig6_claim_retests.pdf`, `.png` | the figure and a 200-dpi preview |
| `Source_Data_Fig6_claim_retests.csv` | Source Data, as in `source_data_submitted/` |
| `Source_Data_Fig2_public_four_protease.csv` | the Figure 2 Source Data, in which the proximal-band rows against the observed background are marked as not plotted |
