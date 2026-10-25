# Xenium breast cancer: tumor and immune compartments

Reproduces the breast-tumor case study of the Results subsection
"Density-correlated genes reveal the spatial antagonism of tumor and immune
compartments in breast cancer". DenMark is fitted to all cells in three 1 mm
by 1 mm regions of interest (ROIs) dominated by invasive tumor, DCIS-1 and DCIS-2
cells, and, across the whole section, separately to DCIS cells (DCIS-1 and
DCIS-2 merged), invasive tumor cells and CD8+ T cells. SPARK-X on the same cells
and genes gives the spatially variable genes for comparison. The section is
replicate 1 of the 10x Xenium human breast cancer data (313 genes), with the
published cell types.

## Scripts

| File | Purpose | Figures and tables |
|---|---|---|
| `DenMark-framework-real-analysis-BC.R` | The single-gene DenMark workflow on this section: ROI windows (section "1B. decide the window of interest"), gene QC, gridding, HSGP fit, latent-field decomposition and spatial correlation functions. Also the gene QC of the whole-section DCIS group (section "for the large range analysis"), the spatial correlation functions of CDH1, POSTN and TRAC (from "Plot: the correlation decay for each of the genes at three regions", line 975) and the Sankey diagrams (from "Sankey plot for large range analysis", line 1576) | Figure S14B; Figure 4F; `DCIS-after-QC.csv` |
| `Viz-for-cluster-genes-short-range.R` | Collects the per-gene ROI fits into one table per ROI; posterior mean of rho against the length of its 95% interval ("Figure for valcano plots", line 267); genes shared across the ROIs by the sign of rho | Figure 4C; the counts of Figure 4D (67 genes share their sign in all three ROIs); Data S2, ROI columns (`DCG-DCIS1.csv`, `DCG-DCIS2.csv`, `DCG-invasive.csv`) |
| `Viz-for-cluster-genes-long-range.R` | Collects the per-gene whole-section fits; DCIS and invasive tumor cells against CD8+ T cells (from line 333); rho of CDH1, POSTN and TRAC with the prior ("Gene of interest ranking", line 662); cell-level and grid-level maps of the three genes (from lines 832 and 1115) | Figure 4E, 4G; Figures S14A and S15; Data S2, whole-section columns (`DCG-DCIS.csv`, `DCG-Tumor.csv`, `DCG-Immune.csv`) |
| `DenMark+SPARK-breast-cancer.R` | SPARK-X on the whole-section populations and on the ROIs; DenMark ranks against SPARK-X ranks ("Redraw into three panels", line 316) and the counts quoted in the text (line 425) | Figure S17 (`rank_dcg_vs_svg_by_celltype_3panel.png`); Data S2, SPARK-X columns (`SPARKSVG-*.csv`) |
| `GSEA-cd8+t-whole.R` | Rank-based GSEA of Gene Ontology terms in CD8+ T cells, genes ranked by the posterior mean of rho | Figure 4H, left (`12_gsea_dotplot_Immune.png`) |
| `Figure4H.R` | rho of LTB, CCR7, TCF7 and CD27 in CD8+ T cells against the LKJ(1.2) prior | Figure 4H, right (`Fig4H_forest_LKJ.pdf`, `.png`) |
| `whole-section/qc-select-genes-whole.R` | Gene QC for the merged DCIS group and the merged tumor class (top 95% of variance, detected in at least 5% of cells) | `qc-whole/<group>_qc_genes.rds`, `.csv`, `<group>_ngenes.txt` |
| `whole-section/DenMark-BC-RA-cluster-whole.R` | One whole-section DenMark fit for one gene of one group (an array task) | `BC-whole/BC-whole-<group>/<GENE>_M2_<k>.csv` |
| `whole-section/cluster-long-whole-DCIS.sh`, `cluster-long-whole-Tumor.sh` | SLURM arrays for group 2 (DCIS-1 and DCIS-2) and group 4 (DCIS-1, DCIS-2 and invasive tumor) | |
| `whole-section/summarize_BCwhole_Astar21.R` | One table of rho per group | `BC-whole-DCIS_Astar21_summary.csv`; `BC-whole-tumor_Astar21_summary.csv` (Data S3) |
| `whole-section/DCG-merged-tumor-vs-fine.R` | The merged tumor class against DCIS and invasive tumor cells analyzed separately | Figure S16C (`fig_composition_merged_vs_fine.png`), S16D (`fig_concordance_merged_vs_fine.png`); `DCG_merged_vs_fine_concordance.csv`, `DCG_merged_vs_fine_long.csv`, `DCG_merged_vs_fine_verdict.txt` |

The scripts are organized in sections that were run in turn.
`DenMark-framework-real-analysis-BC.R` also uses objects that it does not create
(for example `draws_df_HS` and `summary_gene_22`) and is a record of the
workflow rather than a script to run from top to bottom.

## Inputs

Public data: the 10x Genomics Xenium human breast cancer data
(<https://www.10xgenomics.com/products/xenium-in-situ/human-breast-dataset-explorer>),
replicate 1, with the cell types of Janesick et al. (2023, Nature Communications
14, 8353). The scripts read two preprocessed tables built from it: a cell-by-gene
count matrix, `cell_by_gene_expression.csv` (first column the cell, then one
column per gene), and a cell table, `cell_spatial.csv` (cell, `x` and `y` in
micrometers, `celltype`). They look for these tables below `DATA_DIR` at:

- `Cell-report-summary-codes-results/cluster-BC/short-range/cell_by_gene_expression.csv`
  and `10X-breast-cancer/after-preprocess/cell_spatial.csv`
  (`DenMark-framework-real-analysis-BC.R` at the start,
  `Viz-for-cluster-genes-long-range.R`, `DenMark+SPARK-breast-cancer.R`);
- `Cell-report-summary-codes-results/cluster-BC/large-range/cell_by_gene_expression.csv`
  and `cell_spatial.csv` (the correlation-decay sections of
  `DenMark-framework-real-analysis-BC.R`);
- `Dataset/xenium-breast-cancer/rep1/` (both files; the whole-section cluster
  scripts, which expect 313 genes besides the control probes).

The index `k` in the per-gene result names (`<GENE>_M2_<k>.csv`) is the gene's
position among the gene columns of the count matrix used for the fit.

Other inputs below `DATA_DIR`, which no active code in this folder writes:

| Path | Content | Read by |
|---|---|---|
| `Cell-report-summary-codes-results/cluster-BC/short-range/{DCIS-I,DCIS-II,Invasive}/<GENE>_M2_<k>.csv` | cmdstanr summary of one ROI fit per gene (row `Astar[2,1]` holds rho) | `Viz-for-cluster-genes-short-range.R` |
| `Cell-report-summary-codes-results/cluster-BC/large-range/{results-dcis,results-invasive-tumor,results-immune}/<GENE>_M2_<k>.csv` | the same for the whole-section fits of DCIS, invasive tumor and CD8+ T cells | `Viz-for-cluster-genes-long-range.R` |
| `Cell-report-summary-codes-results/cluster-BC/large-range/saved-after-QC-genenames/{Tumor,Immune}-after-QC.csv` | genes kept by the QC of the whole-section invasive tumor and CD8+ T-cell fits | `Viz-for-cluster-genes-long-range.R`, `DenMark+SPARK-breast-cancer.R` |
| `Cell-report-summary-codes-results/cluster-BC/short-range/saved-after-QC-genenames/{DCIS-I,DCIS-II,Invasive}-after-QC.csv` | genes kept by the QC of each ROI; the DCIS-I list is written by line 186 of `DenMark-framework-real-analysis-BC.R` once uncommented (to `OUTPUT_DIR`, so copy it into `DATA_DIR` before running the readers) | `Viz-for-cluster-genes-short-range.R`, `DenMark+SPARK-breast-cancer.R` |
| `Cell-report-summary-codes-results/dataset-multiple-genes-Xenium/{dcis,invasive,immune}_{CDH1,POSTN,TRAC}.rds` | posterior draws of the whole-section fits of CDH1, POSTN and TRAC | `DenMark-framework-real-analysis-BC.R` (Figure S14B) |
| `Cell-report-summary-codes-results/dataset-multiple-genes-Xenium/Gene-names-large-range/{dcis,tumor,immune}-{pos,neg}.csv` | genes with a positive or negative posterior mean in each population; written by lines 320 to 329 of `Viz-for-cluster-genes-long-range.R` once uncommented (to `OUTPUT_DIR`, so copy them into `DATA_DIR` before running the reader) | `DenMark-framework-real-analysis-BC.R` (Figure 4F) |
| `final-summarize-codes/dataset-analysis-genes/analysis_Aqp4.rds` | a saved cmdstanr fit, read at line 463 in place of the fit made just before it | `DenMark-framework-real-analysis-BC.R` |
| `Revision-R1/figures/Figure4Hsheet.xlsx` | sheet `Figure4H`: gene, posterior mean, 2.5% and 97.5% quantiles for CD8+ T cells (apparently the `DCG-Immune.csv` table) | `Figure4H.R` |

Files written by one script and read by another are under `OUTPUT_DIR`:
`DCIS-after-QC.csv` (written by `DenMark-framework-real-analysis-BC.R`), the
`DCG-*.csv` tables (the two `Viz-for-cluster-genes-*` scripts), the
`SPARKSVG-*.csv` tables (`DenMark+SPARK-breast-cancer.R`), `qc-whole/`
(`qc-select-genes-whole.R`) and the two `*_Astar21_summary.csv` tables
(`summarize_BCwhole_Astar21.R`).

## Run order

Regions of interest and whole-section populations:

1. `DenMark-framework-real-analysis-BC.R`: the data-loading section at the top
   ("[Step 0]") and the section "for the large range analysis", which writes
   the DCIS gene list.
2. `Viz-for-cluster-genes-short-range.R` and `Viz-for-cluster-genes-long-range.R`,
   which need the per-gene fits.
3. `DenMark+SPARK-breast-cancer.R` (needs the `DCG-*.csv` tables of step 2),
   `GSEA-cd8+t-whole.R` (needs `DCG-Immune.csv`) and `Figure4H.R`.
4. The correlation-decay and Sankey sections of
   `DenMark-framework-real-analysis-BC.R` (need the posterior draws and the sign
   lists).

Merged classes (Figure S16, Data S3), on a SLURM cluster, from `whole-section/`:

1. Export `DATA_DIR`, `OUTPUT_DIR` and `REPO_ROOT` in the shell and run
   `Rscript qc-select-genes-whole.R` once.
2. Set the same three variables and the allocation in the two `.sh` scripts,
   create `logs/` (`mkdir -p logs`), and submit one array per group:
   `sbatch --array=1-$(cat $OUTPUT_DIR/qc-whole/DCIS_ngenes.txt) cluster-long-whole-DCIS.sh`,
   and the same with `Tumor`. Each task fits one gene with 4 chains of 1,000
   warm-up and 1,000 sampling iterations on a 50 x 50 grid with 25 x 25 basis
   functions (requested: 10 tasks, 4 GB per CPU, 5 hours).
3. Run `summarize_BCwhole_Astar21.R`, which reads the per-gene fits in
   `OUTPUT_DIR/BC-whole/BC-whole-DCIS/` and `OUTPUT_DIR/BC-whole/BC-whole-Tumor/`
   and writes the two summary tables to `OUTPUT_DIR/BC-whole/`.
4. Put `full-division.xlsx` in `DATA_DIR` and run `DCG-merged-tumor-vs-fine.R`,
   which reads the summary tables in `OUTPUT_DIR/BC-whole/` and writes its
   figures and tables there.

## Missing files

- `M2-latent-field-mod.stan` (read from `STAN_DIR` by
  `DenMark-framework-real-analysis-BC.R`) and `M2-latent-field-mod-genes.stan`
  (read from `STAN_DIR` by `whole-section/DenMark-BC-RA-cluster-whole.R`). Both
  are called with the data list of `inst/stan/DenMark.stan`.
- The preprocessing that builds `cell_by_gene_expression.csv` and
  `cell_spatial.csv` from the 10x outputs and the published cell types.
- The cluster scripts that fit each gene in the three ROIs and in the
  whole-section DCIS, invasive tumor and CD8+ T-cell populations, and the QC
  gene lists, posterior draws, `analysis_Aqp4.rds` and `Figure4Hsheet.xlsx`
  listed under Inputs. `DenMark-BC-RA-cluster-whole.R` fits the whole-section
  groups 1 to 3 as well, but needs `qc-whole/Immune_qc_genes.rds` and
  `qc-whole/Invasive_qc_genes.rds`, which `qc-select-genes-whole.R` does not
  write.
- `full-division.xlsx`, read from `DATA_DIR` by `DCG-merged-tumor-vs-fine.R`: one block of
  `Gene`, `Mean`, `Low95` and `High95` for each of DCIS (columns B to E), CD8+ T
  cells (G to J) and invasive tumor (L to O), apparently the three `DCG-*.csv`
  tables of `Viz-for-cluster-genes-long-range.R` placed side by side.
