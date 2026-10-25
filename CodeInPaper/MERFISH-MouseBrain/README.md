# MERFISH mouse brain: density-correlated genes in astrocytes

Reproduces the mouse-brain case study of the Results subsection
"Density-correlated genes in mouse-brain astrocytes complement cell-type markers
and are reproducible across sections". DenMark is fitted to the astrocytes of
the Vizgen MERSCOPE mouse brain data (three coronal positions in triplicate,
483 genes), on a 30 x 30 grid restricted to an alpha-shape tissue domain, with
section S3R1 (slice 3, replicate 1) shown in detail.

## Files

| File | Purpose | Figures and tables |
|---|---|---|
| `S3R1-DenMark-5genes.R` | Section S3R1: cell-type map and composition, the genome-wide summary of the per-gene fits, fits of Aqp4, Sox9, Fgfr3, Ddr2 and Celsr2 with the full and the nested (a21 = 0) model, comparison with marker specificity and with SPARK-X, concordance across the nine sections, GSEA | Figure 3A-K, Figures S7, S8 and S9, Data S1 (S3R1 columns) |
| `depth_diagnostic_MERFISH_S3R1.R` | Per-cell detection depth against local cell density on the interior grids of S3R1 | Figure S13A, B and D |
| `DenMark-M2-vararea-a21zero.stan` | Nested model with a21 = 0 (rho = 0); same data block as `inst/stan/DenMark-vararea.stan`, so the two fits can be compared by WAIC | WAIC comparison in `S3R1-DenMark-5genes.R` (not reported in the manuscript) |
| `grid_area/grid_area_alpha_S*.csv` | Tissue domain of each of the nine sections on the 30 x 30 grid: `grid_index`, `area` (part of the grid inside the alpha-shape boundary, in square micrometers), `centroid_x`, `centroid_y` (micrometers), `type` (`interior` or `boundary`; exterior grids are omitted), `i`, `j` | Input to both scripts (S3R1); the domains are drawn in Figure S5 |

`S3R1-DenMark-5genes.R` was run section by section. Its figures go to
`OUTPUT_DIR/Revision-R1/Reanalysis-MERFISH/reanalysis-S3R1/S3R1-5genes/figures/`,
and its tables and fits to the folder held by `OUT_DIR` at that point, which is
`.../cluster-for-8-rest-slices/summary-S3R1/` from line 290 onward.

| Lines | Section | Output | Panel |
|---|---|---|---|
| 96-243 | cell-type annotation map with the tissue boundary | `00_celltype_map_S3R1.png` | 3A, left |
| 244-283 | cell-type composition figure | `00b_celltype_composition_S3R1.png` | 3A, middle |
| 284-507 | Add genome-wide summary | `panelwide_denmark_S3R1.csv` and `.rds`, `panelwide_denmark_corr_S3R1.png` | 3B; Data S1 |
| 508-598 | Astrocyte enrichment map (contour map) | astrocyte kernel density (drawn, not saved) | 3A, right |
| 599-1047 | Load data; keep ASTROCYTES only; grid-level counts; shared Stan inputs | `01_cell_level_expression.png`, `02_grid_level_counts_expression.png` | 3C, 3D |
| 1048-1276 | FIT ONLY (full and nested model, five genes), with a refit of Aqp4 at `adapt_delta = 0.99` | `fit_S3R1_<gene>.rds`, `waic_aqp4_combined_vs_y2.png` | inputs to the panels below |
| 1277-1593 | compare + visualize the 5 fitted genes | `07_rho_prior_posterior_S3R1.png`; `comparison_S3R1.csv`, `ranges_S3R1.csv`, `waic_comparison_S3R1.csv`, `06_waic_comparison_S3R1.png`, `08_waic_arrows_S3R1.png` | 3E; WAIC comparison |
| 1594-1998 | Correlation decay of the expression and cell-density fields | `16_decay_overlay_S3R1.png` | 3F |
| 1999-2100 | Latent-field decomposition, per gene (`decomp_gene()`) | drawn for each gene; the `ggsave()` calls are commented out | 3G (Aqp4), S7 (the other four) |
| 2101-2226 | Astrocyte marker specificity via `FindMarkers(test.use = "roc")` | `09_rho_vs_specificity_S3R1.png` | 3H |
| 2227-2378 | DenMark vs. SPARK | `10_dcg_vs_svg_rank_S3R1.png`, `dcg_vs_svg_S3R1.csv` | 3I |
| 2379-2540 | Slice by Slice, reproducibility panel | `11a_slice_concordance_heatmap.png`; `stable_core_genes.csv` with the plots `p_core` and `p_core_genes` | 3J; S8, S9 |
| 2541-2611 | GSEA | `12_gsea_dotplot_S3R1.png` | 3K |
| 2612-2657 | Mertk expression by cell type | printed table | Mertk values in the text |

## Inputs

Public data: the Vizgen MERSCOPE mouse brain receptor map, nine sections
(slices 1 to 3, replicates 1 to 3) with a `cell_by_gene` and a `cell_metadata`
table each, in the public Google Cloud bucket `public-datasets-vizgen-merfish`
(<https://console.cloud.google.com/storage/browser/public-datasets-vizgen-merfish>;
Vizgen Data Release Program). The scripts start from derived tables, below
`DATA_DIR`:

- `Revision-R1/Reanalysis-MERFISH/tissue-boundary-9-slices/annotation_out/cell_by_gene_with_celltype_S3R1.csv`
  and `cell_metadata_with_celltype_S3R1.csv`: the S3R1 tables with cell types
  (astrocyte, microglia, neuron, other) assigned from canonical markers (STAR
  Methods). Cell IDs are read as character.
- `Revision-R1/Reanalysis-MERFISH/cluster-for-8-rest-slices/results-S3R1/denmark_summary_S3R1_g*.csv`:
  one posterior summary per gene from the cluster fits of all genes in S3R1.
- `Revision-R1/Reanalysis-MERFISH/cluster-for-8-rest-slices/astar_summaries/astar_summary_<section>.csv`:
  rho summaries of all genes for each of the nine sections (`slice`, `gene_id`,
  `gene_name`, `mean`, `q2.5`, `q97.5`).

From this repository: `inst/stan/DenMark-vararea.stan`,
`DenMark-M2-vararea-a21zero.stan` and `grid_area/grid_area_alpha_S3R1.csv`.

## Run order

1. Upstream, not in this folder (see Missing files): annotate the nine sections,
   estimate their alpha-shape domains (`grid_area/`), fit every gene in every
   section on a cluster and summarize the fits.
2. `S3R1-DenMark-5genes.R`, section by section. The five-gene fits run 4 chains
   of 1,000 warm-up and 1,000 sampling iterations, for the full and the nested
   model, with 25 x 25 basis functions and a boundary factor of 1.5.
3. `depth_diagnostic_MERFISH_S3R1.R`, independent of step 2. It writes
   `MERFISH_S3R1_depth_diagnostic.png` and
   `MERFISH_S3R1_astrocyte_depth_density_perquadrat.csv` to
   `OUTPUT_DIR/Revision-R1/seq-depth-diagnostics/`. Panel C of Figure S13 (depth
   against the density of all cells) is not drawn by this version of the script.

## Missing files

- The cell-type annotation script that writes `annotation_out/` for the nine
  sections, and the alpha-shape boundary script that wrote `grid_area/`.
- The cluster scripts that fit every gene in every section
  (`results-<section>/denmark_summary_<section>_g*.csv`) and the script that
  summarizes them into `astar_summaries/astar_summary_<section>.csv`.
