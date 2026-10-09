# Figure 1 (revision of 2026-09-30)

Panel d shows the final tally of the 28 re-tested claims (`../../final/final_tally.csv`): 11 survive, 2 are attenuated, 13 are undecidable, 1 null claim is broken by the control and 1 baseline runs against the claimed direction. Site-level claims are solid, protein-level claims lighter and hatched. Panels a–c are the submitted page, unchanged.

| file | content |
|---|---|
| `fig1_final.py` | builds the PDF, the PNG preview and the Source Data from the submitted page and the final tally; the build is byte-reproducible and refuses to run unless the final tally equals the `final_verdict` column of Supplemental Data 2 claim by claim |
| `qa_fig1_final.py` | read-only checks: page size and fonts, panels a–c unchanged, every drawn segment and printed number equal to a Source Data row, recount of the tally |
| `Fig1_overview.pdf`, `.png` | the figure (7.2 × 4.33 in) and a 200-dpi preview |
| `Source_Data_Fig1_overview.csv` | Source Data, as in `source_data_submitted/` |

The scripts read the submitted figure from the authors' working tree; the paths are set at the top of each script.
