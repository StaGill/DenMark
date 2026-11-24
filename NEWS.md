# DenMark (development version)

* `offset = "umi"` in `DenMark()` and `denmark_prepare()` (Python: `DenMark()`, `prepare_data()`): the expression counts take the total count of the focal cells in each grid cell as their exposure, `log(T + 0.05)` with `T` summed from a per-cell column named by `umi =`, instead of the number of cells. This is the model of the revised manuscript for sequencing data (`DenMark-vararea-umi.stan`), used for Visium HD in the ovarian analysis; rho then describes the share of the transcriptome taken by the gene rather than expression per cell, so the two offsets are not interchangeable. `buildgridpp(umi =)` returns the summed total counts as `T`, and `denmark_fields()` (Python: `fit.fields()`) adds that exposure, not the cell count, to `log_lambda2`.
* `denmark_pool_visium()` (Python: `pool_visium()`) pools spot-level Visium data into hexagonal super-grids of seven spots, with Mingchi Xu's aperture-7 rule from the colorectal analysis (`CodeInPaper/Visium-ColorectalCancer`) kept verbatim; pools of fewer than three spots are folded into the nearest larger pool. `DenMark(units = TRUE)` fits such pre-aggregated units directly, with their area as offset, alone or with `offset = "umi"`. Maps of these fits (`plot_denmark_grid()`, `plot_denmark_fields()`; Python: `plot_grid()`, `plot_fields()`) draw one point per unit.
* `denmark_scan()` and `scan_genes()` pass the columns that `offset = "umi"` and `units = TRUE` need (`umi`, `cells`, `area`) to every fit and leave them out of the default list of genes, as they do `unit` and `n_spots` of pooled units; a missing column stops the scan before the first fit. Parallel scans compile the model that the fits use.
* The scale of the length-scale prior, `mrange`, is now the largest distance between grid centroids divided by 2 x 2.75, as in every analysis of the paper, so that the prior mean puts the practical range at half that distance. On a full rectangular grid this is the previous value, the diagonal of the grid's bounding box. With a tissue boundary, which removes the grid cells outside the tissue, the largest distance is shorter than that diagonal, so fits with a boundary change slightly and now match the paper. It is computed on the convex hull of the centroids, in both languages; R now imports `grDevices` for `chull()`.
* The default `grid_size` of `DenMark()` and `denmark_prepare()` (Python: `DenMark()`, `prepare_data()`) is now 30 instead of 10, the grid of the mouse-brain and ovarian analyses. A fit that leaves it unset uses a finer grid, takes longer and changes; `grid_size = 10` gives the previous behavior.
* The calls in `class` from `denmark_scan()` and `scan_genes()` are now "positive DCG", "negative DCG" and "uninformative", the terms of the manuscript, in place of "PCG", "NCG" and "not DCG". The numbers are unchanged; code that compares the labels needs the new ones.
* New article [Gene programs on Visium HD and Xenium](https://stagill.github.io/DenMark/articles/programs.html), in R and Python: a program-level analysis from the cell table to the estimates (gene sets, effective number of genes, choice of offset, grid, fit, reading, and comparison of two platforms). It runs on a simulated section with a known answer, `inst/extdata/denmark_programs.csv.gz`, written by `data-raw/make_program_data.R`.
* Parity tests cover the total-count offset and the pooled units: `data-raw/make_toy_data.R` writes two more fixtures, `toy_stan_data_umi.json` and `toy_visium_units.json`, and the parity fixture with the tissue boundary carries the new prior scale.

# DenMark 0.9.0

First release as installable packages. The R functions and the Python module that were sourced from `DenMark/` are now the R package `DenMark`, at the repository root, and the Python package `denmark`, in `python/`. The R functions are those of Mingchi Xu's revision, which adds the tissue boundary (see below), and the Python package carries the same revision. Apart from how the Stan files are found and where compiled models are cached, the functions and the Stan model are carried over unchanged, with one exception that follows the revised manuscript: only the DenMark model is kept. The simplified model with independent fields (formerly `model = "M1"`), the `model` argument and WAIC are removed, with `denmark_compare()` and `fit_simple` in R and `compare_models()`, `waic()` and `fit_simple` in Python, and the Stan files drop the M2 label. Fits of the model give the same draws as `model = "M2"` did for the same Stan data.

* The R package installs from the StaGILL R-universe (`install.packages("DenMark", repos = c("https://stagill.r-universe.dev", "https://stan-dev.r-universe.dev", "https://cloud.r-project.org"))`) and the Python package from GitHub (`pip install "git+https://github.com/StaGill/DenMark#subdirectory=python"`); both build the same Stan data from the same input (with a tissue boundary, the grid areas agree up to the last bit, since sf and shapely compute polygon areas with different formulas).
* Bundled Stan models in `inst/stan/`, with byte copies in the Python package: the DenMark model (`DenMark.stan`), its variable-area version, which the packages fit with a tissue boundary (`DenMark-vararea.stan`), and its total-UMI version used in the revision analyses (`DenMark-vararea-umi.stan`).
* Compiled models are cached in a user directory, since the installed package may be read-only.
* Toy data with known density-expression correlations (geneA 0.6, geneB 0, geneC -0.6): `denmark_toy_data()` in R and `denmark.datasets.load_toy()` in Python.
* Tests in both languages (testthat and pytest), including R/Python parity of the Stan data, with and without the tissue boundary, and runs of the README examples.
* Documentation website (pkgdown) with a Get started vignette and articles on the Python interface, the toy data, the case studies and citation.
* `CodeInPaper/`: the scripts behind the analyses of the revised manuscript, with configurable data and output paths.

## Tissue boundary (revision by Mingchi Xu)

Mingchi Xu's revision of the analysis functions (`DenMark/DenMark-functions.R`, 2026-10-24) restricts inference to the tissue, as the revised manuscript does. The Python package has the same functions, ported from his R code.

* `denmark_alpha_shape()` (Python: `alpha_shape()`) delineates the tissue as an alpha shape of the cell locations: the union of the Delaunay triangles whose circumradius is below alpha, by default three times the median Delaunay edge length. It is computed from all cells of the section, because the region where the analyzed cell type could have been observed is the whole tissue.
* `boundary =` in `buildgridpp()`, `denmark_prepare()`, `DenMark()` and `denmark_scan()` (Python: `buildgridpp()`, `prepare_data()`, `DenMark()`, `scan_genes()`) classifies the grid cells as interior, boundary or exterior and gives each its area inside the tissue; from `denmark_prepare()` on, the exterior grid cells are dropped. `log_grid_area` then becomes a vector, and `DenMark()` fits `DenMark-vararea.stan`. `min_area_frac` optionally drops boundary grid cells with less than that fraction of a full grid cell inside the tissue.
* Without `boundary`, the bounding box is gridded with a common grid area, as before: the Stan data differ from those of the earlier functions only by floating-point rounding of `coords` and `L` (at most 9e-16 on the toy data), and the starting values agree up to the last bit unless the grids hold fewer than 0.1 cells on average.
* `plot_denmark_boundary()` (Python: `plot_boundary()`) shows the cells, the boundary and the classification of the grid cells; `plot_denmark_grid()` and `plot_denmark_fields()` draw the tissue outline and leave out the exterior grid cells.
* `denmark_simulate()` (Python: `simulate()`) gains `tissue = "irregular"`, which places the cells in a lobed tissue section, and `n_other` and `celltype`, which add cells of other types drawn from the locations that the thinning rejected, so that they are densest where the analyzed cells are sparse. The fields are now standardized and orthogonalized on a lattice over the tissue (over space rather than over the cells), so a seed no longer reproduces the data of the earlier functions.
* New toy data: the simulated tissue section of `Tutorial_DenMark.Rmd`, with a `celltype` column, 3,449 Tumor cells that carry the three genes and 4,000 Other cells that only delineate the tissue. The examples compute the boundary from all cells and analyze the Tumor cells.
* Input checks: a non-numeric expression column, such as a cell-type label picked as the default gene, and a boundary made of several separate features now stop with an explanatory error, and `buildgridpp()` reports the cells that lie outside the extent of the boundary.
* `denmark_init()` (Python: `make_init()`) starts `beta0` at the log of the number of cells per unit of tissue area.
* The boundary needs sf in R, a new suggested package, and shapely 2.0 or later in Python, installed with the extra `boundary` (`pip install "denmark[boundary] @ git+https://github.com/StaGill/DenMark#subdirectory=python"`) or with `pip install shapely`.

## Seurat and AnnData

* `denmark_from_seurat()` reads the cell centroids, a cell-type column and the raw counts from a Seurat object with imaging-based data (Xenium, MERSCOPE, CosMx), and `denmark_annotate()` writes the scan results into its gene metadata; SeuratObject is a new suggested package. In Python, `from_anndata()` reads an AnnData object (scanpy, squidpy, or a SpatialData table from spatialdata-io) and `annotate()` writes the results into `adata.var`.

## Scans of whole panels

* `denmark_scan(workers = )` and `scan_genes(n_jobs = )` fit several genes at a time (R: separate sessions through future.apply, a new suggested package; Python: threads, each driving its own CmdStan processes). Every gene keeps its seed, so the results equal those of the sequential scan.
* `inst/slurm/` holds a template that fits one gene per SLURM array task, can be resubmitted after an interruption, and collects the results with `denmark_scan()`.

## Help for users

* A FAQ page (installation, input data, grid and basis settings, convergence, reading the results), a Google Colab notebook of the toy-data workflow, and issue templates for bug reports and questions about data.

## Planned for 1.0.0

* Total-UMI expression offset for sequencing-based platforms (Visium HD and Visium).
* Pooling of Visium spots into seven-spot super-grids.
