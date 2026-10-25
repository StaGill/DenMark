# Ovarian cancer: program-level correlations across Visium HD and Xenium

Reproduces the cross-platform comparison of the Results subsection "DenMark
recovers program-level density correlations that transfer across platforms".
Tumor-cell density is correlated with four expression programs (proliferation,
epithelial identity, hypoxia and an EMT transcription-factor program) on adjacent
sections of one ovarian tumor profiled by Visium HD, Xenium v1 and Xenium Prime
5K. In each comparison, the programs are restricted to the genes shared by the
two panels, the total transcript count over those genes in each tile is the
expression offset on both platforms, and every section is gridded into 30 x 30
tiles over an alpha-shape domain.

## Notebook

| File | Purpose | Figures and tables |
|---|---|---|
| `denmark_xplatform_replication_umi_30x30-cleaned.ipynb` | Loads the three platforms, intersects their panels, marks tumor cells, and fits `inst/stan/DenMark-vararea-umi.stan` for each program (and for several single genes) on Visium HD and on the Xenium section of each pair | Figure 5G and Figure S19A (`fig_forest_HD_5K_programs.png`); Figure S19C (`fig_forest_HD_v1_programs.png`); tumor cells on each section (`fig_qc_tumour_mask.png`, compare Figure 5F); program genes on each shared panel (`fig_program_gene_survival.png`, compare Figure 5B); panel size and cell numbers (`fig_dataset_summary.png`); Data S4, shared-panel columns (`denmark_rho_pairwise_umi.csv`, `res_HD_a.csv`, `res_v1.csv`, `res_HD_b.csv`, `res_5K.csv`) |

Outputs go to `OUTPUT_DIR/xplatform_out/` (with `rho_draws.pkl`, `aux_all.pkl`
and `session_meta.pkl` for replotting) and the CmdStan files to
`OUTPUT_DIR/stan_out/`.

## Inputs

Public data from 10x Genomics: Visium HD
(<https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-ovarian-cancer-discovery-fresh-frozen>)
and Xenium v1 and Xenium Prime 5K
(<https://www.10xgenomics.com/datasets/xenium-comparison-fresh-frozen-human-ovarian-cancer>).
Below `DATA_DIR`, the notebook reads:

- `visium_ovarian/results/tables/visium_hd_ovarian_cancer_annotated.h5ad`: the
  Visium HD cells (segmented by Space Ranger, annotated with CellTypist), with
  raw counts and cell centroids in micrometers. It is produced by the authors'
  Visium HD pipeline, which is not in this repository (see Missing files).
- `Revision-R1/VisiumHD-Breast_Cancer/Xenium V1/` and
  `Revision-R1/VisiumHD-Breast_Cancer/Xenium 5K/`: the Xenium Onboard Analysis
  output folders of the two Xenium sections, each holding
  `cell_feature_matrix.h5` and `cells.parquet` (or `cells.csv.gz`). The parent
  folder name is the authors' original; it holds the ovarian data.

## How to run

Use Python 3.11 with numpy, pandas, scipy, matplotlib, scanpy, shapely, pyarrow
and CmdStanPy, and CmdStan 2.39 (on Windows, CmdStanPy also needs the RTools
toolchain). Edit the CONFIGURE ME cell, the first code cell, and run the
notebook from the top. Each fit runs 4 chains of 1,500 warm-up and 1,000
sampling iterations (`adapt_delta` 0.95, maximum tree depth 12) with 4 x 4 basis
functions and a boundary factor of 1.5, and is kept when R-hat is below 1.1, the
bulk effective sample size above 150 and the number of divergent transitions at
most 20.

The stored output of the cell that runs the comparisons ends with a GEOS memory
error raised by a three-way comparison that the cleaned notebook no longer
contains; the two pairwise comparisons had finished and were saved before it.

## Missing files

- The annotated Visium HD object
  `visium_hd_ovarian_cancer_annotated.h5ad` and the notebook that builds it. The
  same pipeline gives Figures 5A, 5C-E and S18 and the full-panel fits of Figure
  5E and Data S4.
- The fits with the tumor-cell count per tile as the offset (Figures S19B and
  S19D).
