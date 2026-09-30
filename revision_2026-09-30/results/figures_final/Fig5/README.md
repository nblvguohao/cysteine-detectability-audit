# Figure 5 (revision of 2026-09-30)

Panel a shows the absolute out-of-fold AUCs of the visibility (VIS10) and digest (DIG25) models, over all sites and within proteins, under the two negative-set constructions (Supplemental Note 20); every plotted value is a stored cell of `results/ptm_detectability_share.csv` (copied in `../../../inputs/repo_results/`). Panels b and c are carried over from the submitted figure. Colors in b and c give verdicts, as in Figure 6.

| file | content |
|---|---|
| `fig5_rebuild.py` | draws panel a and places panels b and c of the submitted page at their original coordinates; checks every plotted value against the stored table and the tables of `../../I_fig5a_estimand/` |
| `variant_SNO-012_undecidable/recolour_fig5c_sno012.py` | final step: recolors the `SNO-012` row of panel c to the undecidable gray (its baseline interval covers zero) and adds the two counts printed in panel b's block headers to the Source Data; nothing else changes |
| `Fig5_three_axes.pdf`, `.png` | the figure after that step (7.2 × 5.91 in) and a 200-dpi preview |
| `Source_Data_Fig5_three_axes.csv` | Source Data, as in `source_data_submitted/` |
