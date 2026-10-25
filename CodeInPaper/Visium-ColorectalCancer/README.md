# Spot-level Visium: Visium-like spots from Xenium and a colorectal cancer section

Reproduces parts of the Results subsection "DenMark extends to spot-level Visium
data and recovers single-cell density correlations after aggregation". Visium-like
spots are built on the Xenium breast section to see what aggregation does to
rho, and DenMark is applied to a Visium CytAssist section of human colorectal
cancer, with cells counted from nuclei segmented in the H&E image and the total
UMI per super-grid as the expression offset.

## Files

| File | Purpose | Figures and tables |
|---|---|---|
| `segment_nuclei_settled.ipynb` | Segments nuclei in the full-resolution H&E image with Cellpose-SAM on a GPU, block by block, and writes their centroids; clusters the spots (Leiden) and labels a cluster as tumor when its mean epithelial marker score (EPCAM, KRT8, KRT18, KRT20, CDX2, CEACAM5) exceeds its mean stromal score | Figure 6E (`nuclei_seg_cellpose_zoom.png`; overview `nuclei_seg_cellpose.png`); `nuclei.csv` (547,877 nuclei), `nuclei_seg.npz`; `spot_celltype.csv` (5,487 tumor and 3,593 non-tumor spots), `spot_celltype_map.png` |
| `DenMark_realVisium_CRC_celltype.R` | Assigns each nucleus to the nearest spot center within 27.5 micrometers, keeps tumor spots with at least one nucleus, pools them into seven-spot super-grids (units of fewer than three spots are merged into the nearest super-grid), and fits DenMark to MKI67, MYC and PPIA with the log total UMI per super-grid as the expression offset | Figure 6G (`Fig_rho_CRC_tumor_MKI67_MYC.pdf`, `.png`); the super-grids summarized in Figure 6F (`Fig_pooling_CRC_tumor.png`); overview `Fig_realVisium_CRC_tumor.pdf`; `realVisium_CRC_tumor_Astar21.csv` |
| `DenMark_pseudoVisium_celltype_Astar21.R` | For one population (`CELLTYPE`: invasive tumor, DCIS or CD8+ T cells), places Visium-like spots (56 micrometers wide, 100 micrometers apart, on a hexagonal lattice) on the Xenium breast section and counts cells and transcripts per spot; pools the spots into super-grids and fits four genes; counts the spots each population occupies | Figure 6C (panels a and b of `Fig_pseudoVisium_<celltype>.pdf`); spot occupancy behind the counts in the text, 837 to 1,203 of 4,358 spots (`pseudoVisium_spot_occupancy.csv`); `pseudoVisium_<celltype>_Astar21.csv` |

The notebook and `DenMark_realVisium_CRC_celltype.R` write to
`OUTPUT_DIR/Revision-R1/spot-Visium/`; `DenMark_pseudoVisium_celltype_Astar21.R`
writes to `OUTPUT_DIR`.

The pooled fits of `DenMark_pseudoVisium_celltype_Astar21.R` are not reported in
the manuscript. There, each population is fitted at the spot level without
pooling, for 252 genes (Figure 6D and Figure S20); those cluster scripts are
not yet in this folder.

## Inputs

Colorectal cancer: the 10x Genomics Visium CytAssist human colorectal cancer
section, 11 mm capture area, FFPE
(<https://www.10xgenomics.com/datasets/human-colorectal-cancer-11-mm-capture-area-ffpe-2-standard>),
in `DATA_DIR/Revision-R1/spot-Visium/`:

- `CytAssist_11mm_FFPE_Human_Colorectal_Cancer_filtered_feature_bc_matrix.h5`
- `CytAssist_11mm_FFPE_Human_Colorectal_Cancer_tissue_image.tif` (full-resolution H&E)
- `CytAssist_11mm_FFPE_Human_Colorectal_Cancer_spatial/spatial/` (the unpacked
  spatial archive, with `tissue_positions.csv` and `scalefactors_json.json`)

`DenMark_realVisium_CRC_celltype.R` reads `nuclei.csv` and `spot_celltype.csv`
from `OUTPUT_DIR/Revision-R1/spot-Visium/`, where the notebook writes them.

Visium-like spots: the preprocessed tables of the Xenium breast section (see
`../Xenium-BreastCancer/README.md`) in
`DATA_DIR/Revision-R1/Visium-Breast-Cancer/xenium-breast-cancer/rep1/`:

- `cell_spatial.csv` (cell, `x` and `y` in micrometers, cell type) and
  `cell_by_gene_expression.csv`;
- `pixel-size-raw.txt` (micrometers per pixel, 0.5) and `radius-raw.txt` (spot
  radius in pixels, which gives 28 micrometers);
- optionally `xenium_singlecell_coupling.csv` (columns `celltype`, `gene`,
  `sc_coupling`), the single-cell estimates for the recovery plot; without it
  that plot is skipped.

## Run order

1. `segment_nuclei_settled.ipynb` (GPU; Cellpose 4 with the Cellpose-SAM model,
   scanpy for the spot clusters). Edit the CONFIGURE ME cell and run the cells
   in order. The overview plot of Step 6 uses `step`, `W_plot` and `plt`, which
   are first set in later cells (the tumor-spot cell and the Panel E cell); set
   them before that cell or skip it, since `nuclei.csv` is written in Step 8.
2. `DenMark_realVisium_CRC_celltype.R`. Each fit runs 4 chains of 1,500 warm-up
   and 1,500 sampling iterations (`adapt_delta` 0.97) with 8 x 8 basis
   functions. The last section replots rho for MKI67 and MYC from the saved
   table without refitting.
3. `DenMark_pseudoVisium_celltype_Astar21.R`, independent of steps 1 and 2, once
   per `CELLTYPE`.

## Missing files

- `DenMark-M2-vararea-libsize.stan`, read from `STAN_DIR` by
  `DenMark_realVisium_CRC_celltype.R`: the variable-area model with the
  centered log library size of each super-grid as the expression offset (data
  `log_lib`).
- `pixel-size-raw.txt`, `radius-raw.txt` and the optional
  `xenium_singlecell_coupling.csv` for the Visium-like spots, and the
  preprocessed Xenium tables (see `../Xenium-BreastCancer/README.md`).
- The scripts behind Figures 6A and 6B (simulated Visium spots), Figure 6D and
  Figure S20 (spot-level fits of 252 genes against single-cell estimates) and the
  maps of Figure 6F.
