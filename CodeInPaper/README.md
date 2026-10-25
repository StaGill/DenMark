# Reproducing the DenMark manuscript analyses

This folder holds the scripts and notebooks behind the analyses and figures of
Xu, Schmidt and Zhang, "DenMark: A Bayesian Hierarchical Model for Identifying
Cell-Density Correlated Genes from Spatial Transcriptomics" (Cell Reports
Methods, under review). Each subfolder covers one case study. The code is kept as
it was run for the paper; only file paths were changed (see [Paths](#paths)).

| Subfolder | Case study | Results subsection | Figures | Dataset access |
|---|---|---|---|---|
| [`Simulation/`](Simulation/) | Calling density-correlated genes (DCGs) among simulated genes of known status | Validation of DenMark on simulated data | S3 | Simulated; no external data |
| [`MERFISH-MouseBrain/`](MERFISH-MouseBrain/) | Astrocytes in nine MERFISH mouse-brain sections, with section S3R1 in detail | Density-correlated genes in mouse-brain astrocytes complement cell-type markers and are reproducible across sections | 3A-K; S7, S8, S9; S13A, B, D | Vizgen MERSCOPE mouse brain data, public bucket `public-datasets-vizgen-merfish`: <https://console.cloud.google.com/storage/browser/public-datasets-vizgen-merfish> |
| [`Xenium-BreastCancer/`](Xenium-BreastCancer/) | Three regions of interest and three whole-section cell populations of a Xenium breast tumor | Density-correlated genes reveal the spatial antagonism of tumor and immune compartments in breast cancer | 4C-H; S14, S15, S16C-D, S17; Data S2, S3 | 10x Genomics Xenium human breast cancer: <https://www.10xgenomics.com/products/xenium-in-situ/human-breast-dataset-explorer> |
| [`OvarianCancer-VisiumHD-Xenium5K/`](OvarianCancer-VisiumHD-Xenium5K/) | Density-program correlations on adjacent ovarian tumor sections profiled by Visium HD, Xenium v1 and Xenium Prime 5K | DenMark recovers program-level density correlations that transfer across platforms | 5G (with 5B, 5F); S19A, S19C; Data S4 (shared panels) | 10x Genomics Visium HD: <https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-ovarian-cancer-discovery-fresh-frozen>; Xenium v1 and Xenium Prime 5K: <https://www.10xgenomics.com/datasets/xenium-comparison-fresh-frozen-human-ovarian-cancer> |
| [`Visium-ColorectalCancer/`](Visium-ColorectalCancer/) | Visium-like spots built from the Xenium breast section, and a Visium colorectal cancer section with nuclei counted from H&E | DenMark extends to spot-level Visium data and recovers single-cell density correlations after aggregation | 6C, 6E, 6G (with the super-grids of 6F) | 10x Genomics Visium CytAssist human colorectal cancer, 11 mm capture area, FFPE: <https://www.10xgenomics.com/datasets/human-colorectal-cancer-11-mm-capture-area-ffpe-2-standard>; Xenium breast data as above |

Figure numbers refer to the main text and the Supplemental Information of the
revised manuscript. Each subfolder README maps scripts to panels, lists the
inputs of every script and gives the run order.

## How to run

1. Clone the repository and, once it is tagged, check out v1.0.0, the release
   that will accompany the manuscript: `git clone https://github.com/StaGill/DenMark.git`,
   then `git checkout v1.0.0` inside the clone. The scripts read the Stan models from
   `inst/stan/` of this checkout and call CmdStan through cmdstanr (R) or
   CmdStanPy (Python); install CmdStan once with `cmdstanr::install_cmdstan()`.
   Each script loads its other packages at the top.
2. Download the public data listed above. The subfolder READMEs list the files
   each script reads and where it expects them below `DATA_DIR`.
3. Open the script or notebook. Its first lines are a **CONFIGURE ME** block that
   sets `DATA_DIR` (inputs), `OUTPUT_DIR` (results and figures) and `REPO_ROOT`
   (the root of this checkout), plus `STAN_DIR` where a model is not yet in the
   repository. Edit them for your environment. The defaults assume that R runs
   from the repository root and that Jupyter starts in the notebook's folder;
   with them, every script reads from `data/` and writes to `output/` at the
   repository root.
4. Run the scripts of a case study in the order given in its README, with the
   same `OUTPUT_DIR`.

Several scripts are organized in sections that were run in turn, and some
sections use objects created in other sections or sessions; the subfolder
READMEs name the sections that produce each panel. Stan compiles each model into an
executable next to its `.stan` file, so the first fit writes an executable into
`inst/stan/` or into the subfolder that holds the model.

## Paths

The scripts were written on the Compute Canada (Digital Research Alliance of
Canada) clusters and on a Windows workstation. Their absolute path prefixes have
been replaced by `DATA_DIR` for inputs and `OUTPUT_DIR` for outputs, and the
sub-folders below each prefix are kept so that the original layout stays
visible, for example
`file.path(DATA_DIR, "Revision-R1/Reanalysis-MERFISH/tissue-boundary-9-slices/annotation_out")`.

- Output files that were written to the working directory are written to
  `OUTPUT_DIR`; in `DenMark_realVisium_CRC_celltype.R` they go to
  `OUTPUT_DIR/Revision-R1/spot-Visium/`, where the second part of that script
  reads them back.
- A file that one script writes and another script (or a later section of the
  same script) reads is read from `OUTPUT_DIR`. Every other input is read from
  `DATA_DIR`. Pointing `DATA_DIR` and `OUTPUT_DIR` at the same folder therefore
  reproduces the authors' layout.
- Stan models shipped with the package are read from `inst/stan/`; the nested
  model of the mouse-brain analysis and the alpha-shape grid files are read from
  `MERFISH-MouseBrain/`. Models that are not yet in the repository are read from
  `STAN_DIR` (see [Missing files](#missing-files)).
- The SLURM scripts export `DATA_DIR`, `OUTPUT_DIR` and `REPO_ROOT` to the R
  script they launch, write their logs to `logs/`, and name the allocation
  `def-<account>`.
- Where an original script read a file from another folder, or with other
  letter case, than the script that writes it, the read path now matches the
  write path: `Simulation/04_experiment_TPR.R` (`FPR+FDR/`),
  `Xenium-BreastCancer/DenMark+SPARK-breast-cancer.R` (`SPARKSVG-DCIS1.csv`,
  `SPARKSVG-DCIS2.csv`) and the two summary scripts of
  `Xenium-BreastCancer/whole-section/` (`OUTPUT_DIR/BC-whole/`).

Nothing else in the code was changed. The stored outputs of the two notebooks
are those of the authors' runs, with the user-specific path prefixes in the
printed messages replaced by `<DATA_DIR>`, `<OUTPUT_DIR>`, `<REPO_ROOT>` or `~`.

## Software

The analyses in the paper used R 4.3.2, Python 3.11.15, CmdStan 2.39.0,
cmdstanr 0.8.1, CmdStanPy 1.3.0, SPARK 1.1.1 (SPARK-X), clusterProfiler 4.10.1,
CellTypist 1.7.1, Cellpose 4.1.1 (Cellpose-SAM) and Space Ranger 4.0.1 (Key
resources table of the STAR Methods). The SLURM scripts load the cluster's R
4.3.1 module.

## Missing files

The scripts need the following files, which are not in this repository and are
not written by any script in this folder. The subfolder READMEs give the
expected location of each.

**Stan models**

- `DenMark-M2-vararea-libsize.stan`, read by
  `Visium-ColorectalCancer/DenMark_realVisium_CRC_celltype.R` (variable-area
  model with the per-unit library size as the expression offset, data
  `log_lib`).
- `M2-latent-field-mod.stan`, read by
  `Xenium-BreastCancer/DenMark-framework-real-analysis-BC.R`, and
  `M2-latent-field-mod-genes.stan`, read by
  `Xenium-BreastCancer/whole-section/DenMark-BC-RA-cluster-whole.R`. Both are
  called with the data list of `inst/stan/DenMark.stan`.

**Upstream scripts and the files they write**

- Mouse brain: the cell-type annotation of the nine sections
  (`annotation_out/cell_by_gene_with_celltype_<section>.csv` and
  `cell_metadata_with_celltype_<section>.csv`); the alpha-shape tissue boundary
  that produced `MERFISH-MouseBrain/grid_area/`; the cluster fits of all genes
  in each section (`cluster-for-8-rest-slices/results-S3R1/denmark_summary_S3R1_g*.csv`)
  and their per-section summaries (`astar_summaries/astar_summary_<section>.csv`).
- Breast tumor: the preprocessing that turns the 10x outputs and the published
  cell types into `cell_by_gene_expression.csv` and `cell_spatial.csv`; the
  cluster fits of each gene in the three regions of interest
  (`cluster-BC/short-range/{DCIS-I,DCIS-II,Invasive}/`) and in the whole-section
  DCIS, invasive tumor and CD8+ T-cell populations
  (`cluster-BC/large-range/results-{dcis,invasive-tumor,immune}/`); the gene
  lists `Tumor-after-QC.csv`, `Immune-after-QC.csv`, `DCIS-I-after-QC.csv`,
  `DCIS-II-after-QC.csv` and `Invasive-after-QC.csv`; the posterior draws
  `{dcis,invasive,immune}_{CDH1,POSTN,TRAC}.rds`; `analysis_Aqp4.rds`;
  `Figure4Hsheet.xlsx`; `full-division.xlsx`. The six sign lists in
  `Gene-names-large-range/` are written by lines 320 to 329 of
  `Viz-for-cluster-genes-long-range.R`, which are commented out.
- Ovarian tumor: the CellTypist-annotated Visium HD object
  `visium_hd_ovarian_cancer_annotated.h5ad` and the notebook that builds it.
- Spot-level Visium: `pixel-size-raw.txt` and `radius-raw.txt` (one number each)
  for the Visium-like spots, and the optional `xenium_singlecell_coupling.csv`.

**Figures without a script here**

- Figure 2 and Figures S1, S2 and S4 (accuracy of the Hilbert space
  approximation, grid resolution, abundance and density contrast, missed cells).
- Figures S5 and S6 (tissue domains and cell types of the nine mouse-brain
  sections), and the per-section fits behind Figure 3J and Figures S8 and S9.
- Figures 4A, 4B, S11, S12 and S16A-B (breast-tumor maps, region grids, density
  and detection-depth diagnostics).
- Figures 5A, 5C-E and S18, the full-panel fits of Figure 5E and Data S4, and the
  cell-count offset fits of Figures S19B and S19D.
- Figures 6A, 6B, 6D and S20 and the maps of Figure 6F (simulated Visium spots,
  and the spot-level fits of 252 genes compared with single-cell estimates).
