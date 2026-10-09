# figures_r31/ - figure set of the R31 working version (superseded)

This folder keeps the figure files of the R31 working version of the manuscript (23 September 2026), with
their source data in `../source_data_r31/`. They are **not** the submitted figures:

- the figures as submitted are in `../figures_submitted/` (Figures 1-7 and the graphical abstract), with
  their source data in `../source_data_submitted/`;
- Figures 1, 3, 5 and 6 were rebuilt after the primary analyses; their build scripts, checks and outputs are in
  `../posthoc_2026-09-30/results/figures_final/`.

Contents: Figures 1, 2, 3, 4 and 6 (PDF and PNG; Figure 2 also as SVG). **There is no Figure 5 here**:
the R31 Figure 5 and its source table were not released (`../source_data_r31/SOURCE_DATA_INDEX_R31.csv`
still lists `Source_Data_Fig5_three_axes.csv`,
which is not in `../source_data_r31/`). The submitted Figure 5 is `../figures_submitted/Fig5_three_axes.pdf`,
built by `../posthoc_2026-09-30/results/figures_final/Fig5/fig5_rebuild.py`, with source data
`../source_data_submitted/Source_Data_Fig5_three_axes.csv`.

Of the five files, only Figure 2 is byte-identical to the submitted version; Figures 1, 3, 4 and 6
differ from the submitted ones.
