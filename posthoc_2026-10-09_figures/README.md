# Figure build of manuscript v3.3 (2026-10-09)

These scripts drew the figures of manuscript v3.3 (`../figures_final/`: main Figures 1-6, Supplemental Figures
S1-S7 and the graphical abstract) and wrote their Source Data (`../source_data_final/`, except the five
`Source_Data_text_*` files, which are copied unchanged from `../source_data_submitted/`). The figures were
redrawn for the new numbering; **no analysis was re-run and no value was re-estimated**. Every plotted or printed
value is read from a file of this repository and written to the figure's Source Data with that file's name.

| script | output |
|---|---|
| `fig1_overview.py` | Figure 1 (panel a re-uses the workflow illustration embedded in `../figures_submitted/Fig1_overview.pdf`) |
| `fig2_negative_sets.py` | Figure 2 |
| `fig3_cleavage_geometry.py` | Figure 3 |
| `fig4_removed_comparisons.py` | Figure 4 |
| `fig5_claim_retests.py` | Figure 5 and Supplemental Figure S4 |
| `fig6_predictors.py` | Figure 6 |
| `supp_figures.py` | Supplemental Figures S1, S2, S3, S6 (the earlier Figures 2, 3, 7, 4, copied unchanged), S5 (panels b and c of the earlier Figure 5, placed unchanged) and S7 (the earlier Figure 1b, redrawn) |
| `graphical_abstract.py` | graphical abstract |
| `common.py` | shared style, input paths and the canvas checks used by every script |
| `verify_values.py` | independent re-read of the plotted values (writes the two tables in `qa/`) |

**Inputs** (all in this repository): `../source_data_submitted/` and `../figures_submitted/` (the earlier figures
and their Source Data, whose values are re-used), the tables of Supplemental Notes 4, 7, 12, 14, 20 and 22 in
`../supplemental/` (unchanged in v3.3), `../artifacts12/positional_profile/results/positional_kr_profile_public_2026-09-17.csv`,
`../artifacts12/results/12_matched_{sno,so,ox}.json` and `../posthoc_2026-09-30/results/G_plmsnosite_rescore/`.
The scripts locate the repository root from their own location; they contain no machine-specific paths.

**Running.** From the repository root, for example `python3 posthoc_2026-10-09_figures/fig1_overview.py`, then
the other figure scripts. Outputs go to `posthoc_2026-10-09_figures/build_output/` (not tracked; set `OUT_DIR` to
write elsewhere). A script refuses to overwrite an existing output, refuses to write into `figures_final/` or
`source_data_final/`, and refuses to save a figure whose canvas check reports a text overlap, text off the page,
text below 6 pt or a printed number not found in its Source Data. The text face is Helvetica, extracted at run time
from the macOS system font `/System/Library/Fonts/Helvetica.ttc`; the build stops if that file is missing, so it
runs on macOS only. Versions used: Python 3.9.6, matplotlib 3.9.4, PyMuPDF 1.26.5, fontTools 4.60.2, pandas 2.3.3,
numpy 2.0.2.

**Reproduction check** (run on 2026-10-09 in a clean checkout of this release, macOS 26.6):

- all 14 PNG files are byte-identical to those in `../figures_final/`;
- of the 14 PDF files, 5 are byte-identical (S1, S2, S3, S6 and S7); the other 9 differ only in the embedded
  Helvetica font program (the `FontFile2` stream), whose `head.modified` timestamp is rewritten when the font is
  extracted at run time; their pages rendered at 150 dpi and their text layers are identical;
- the 14 Source Data files written by the scripts are identical to those in `../source_data_final/` after line
  endings are normalised (the scripts write CRLF; the repository stores LF).

**Value checks** (`qa/`, written by
`python3 posthoc_2026-10-09_figures/verify_values.py figures_final source_data_final posthoc_2026-10-09_figures/qa`):

- `verify_values.csv`: 1,007 plotted values re-read with separate parsing code from the table the manuscript cites
  for them (Supplemental Notes 4, 7, 12, 14, 20, 22, the positional-profile CSV, the matched JSON files, the
  pLMSNOSite result tables and the earlier Source Data), compared at that table's precision; 0 mismatches;
- `printed_numbers.csv`: the 307 numeric tokens in the text layers of Figures 1-6, S4 and S7 with the Source Data
  rows they match; the 48 tokens without a match are axis tick labels, claim identifiers (for example `SNO-021`)
  and attribute descriptions in row labels (for example `within ±5`), not data values.
