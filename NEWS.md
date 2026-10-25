# DenMark 0.9.0

First release as installable packages. The R functions and the Python module that were sourced from `DenMark/` are now the R package `DenMark`, at the repository root, and the Python package `denmark`, in `python/`. The R functions are those of Mingchi Xu's revision, which adds the tissue boundary (see below), and the Python package carries the same revision. Apart from how the Stan files are found and where compiled models are cached, the functions and the Stan model are carried over unchanged, with one exception that follows the revised manuscript: only the DenMark model is kept. The simplified model with independent fields (formerly `model = "M1"`), the `model` argument and WAIC are removed, with `denmark_compare()` and `fit_simple` in R and `compare_models()`, `waic()` and `fit_simple` in Python, and the Stan files drop the M2 label. Fits of the model give the same draws as `model = "M2"` did for the same Stan data.

* Both packages install from GitHub (`remotes::install_github("StaGill/DenMark")`; `pip install "git+https://github.com/StaGill/DenMark#subdirectory=python"`) and build the same Stan data from the same input (with a tissue boundary, the grid areas agree up to the last bit, since sf and shapely compute polygon areas with different formulas).
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

## Planned for 1.0.0

* Total-UMI expression offset for sequencing-based platforms (Visium HD and Visium).
* Pooling of Visium spots into seven-spot super-grids.
