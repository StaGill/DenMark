# AGENTS.md

Guidance for coding agents working in this repository. Users of the packages should start from `README.md` or https://stagill.github.io/DenMark/.

## What this is

DenMark estimates, for one cell population in a spatial transcriptomics section, how the expression per cell of each gene varies with the local density of that population, and calls density-correlated genes (DCGs). It grids the cells of the population over the tissue (an alpha shape of all cells of the section), models cell counts and summed expression jointly as a marked log-Gaussian Cox process whose two latent Gaussian fields are approximated by a Hilbert space Gaussian process (HSGP), and samples the posterior with CmdStan. The density-expression correlation is `rho = Astar[2,1]`.

Two packages fit the same Stan models to the same gridded data: the R package `DenMark` at the repository root (cmdstanr) and the Python package `denmark` in `python/` (cmdstanpy). The method is described in Xu, Schmidt and Zhang, "DenMark: A Bayesian Hierarchical Model for Identifying Cell-Density Correlated Genes from Spatial Transcriptomics" (bioRxiv, doi:10.64898/2026.04.02.713482); the revised manuscript is under review at *Cell Reports Methods*. Version 0.9.0 packages the scripts that were sourced from `DenMark/`: the R functions of Mingchi Xu's revision (`DenMark/DenMark-functions.R`, GitHub commit 4c68cc9, 2026-10-24), which added the alpha-shape tissue boundary, and the Python module `DenMark/denmark.py` of commit e8b4cc9, with the revision ported to it. Release v1.0.0 will add the total-UMI offset and Visium spot pooling, and will accompany the manuscript.

## Set up and test

R, from the repository root:

```r
cmdstanr::install_cmdstan()   # once per machine; cmdstanr comes from https://stan-dev.r-universe.dev
devtools::document()          # regenerates NAMESPACE and man/ from the roxygen comments
devtools::test()
devtools::check()
devtools::install()
pkgdown::build_site()
```

Python:

```bash
pip install -e "python[test,plot,boundary]"
python -m cmdstanpy.install_cmdstan   # once per machine
pytest -q python/tests
```

Tests that fit a model skip when CmdStan is not installed; gridding, Stan data, toy data, R/Python parity and the Stan-file copies are tested without it. The tissue boundary needs sf in R (a suggested package) and shapely 2.0 or later in Python (the `boundary` extra), and so do the README examples. Both test suites run the code blocks of `README.md` as written, on the toy data. Blocks that read the reader's own files (`read.csv("...")`, `pd.read_csv("...")`) are skipped, and blocks that fit a model need CmdStan. An example that stops working therefore fails the build, so keep install commands in `bash` blocks rather than `r` or `python` blocks.

## CI

- `.github/workflows/R-CMD-check.yaml`: R CMD check on Ubuntu with R release and oldrel-1. CmdStan is built once per CmdStan release and cached, so the fit tests run. The built package has no `README.md`, so the release leg also runs `tests/testthat/test-readme.R` from the source tree.
- `.github/workflows/python-tests.yml`: pytest on Python 3.9 and 3.13 with the latest releases of the dependencies. Only the 3.13 leg installs CmdStan (cached) and runs the fits. The floors declared in `python/pyproject.toml` (numpy 1.21, pandas 1.3, scipy 1.7, cmdstanpy 1.2) are not tested yet; once they are, pin them on the 3.9 leg and change them in both files together.
- `.github/workflows/pkgdown.yaml`: builds the website, including the example fits of the vignette and articles, and deploys it to GitHub Pages. Pages deploys only once the repository is public (or on a paid plan).
- `.github/workflows/publish.yml`: see Releases.

## Layout

- `R/`, split along the section headers of the former `DenMark/DenMark-functions.R`:
  - `utils.R`: `q2.5()` to `q97.5()`, `matern32()` and internal input helpers.
  - `prepare.R`: `denmark_filter_genes()`, `denmark_alpha_shape()`, `buildgridpp()`, `grid_pts_centroid()`, `denmark_hsgp_indices()`, `denmark_prepare()`.
  - `fit.R`: Stan-file lookup, `denmark_model()`, `denmark_init()`, `DenMark()`, `DenMarkmain()`.
  - `postprocess.R`: `denmark_summary()`, `denmark_draws()`, `denmark_fields()`, `denmark_decay()`, `denmark_check_hsgp()`, `print.denmark()`, `denmark_save()`, `denmark_load()`.
  - `scan.R`: `denmark_scan()`.
  - `plot.R`: `plot_denmark_boundary()`, `plot_denmark_grid()`, `plot_denmark_fields()`, `plot_denmark_decay()`, `plot_denmark_scan()`.
  - `simulate.R`: `denmark_simulate()`.
  - `toy.R`: `denmark_toy_data()`.
  - `zzz.R`: the package environment that caches compiled models.
  - `DenMark-package.R`: package documentation and imports.
- `inst/stan/`: `DenMark.stan` (the DenMark model), `DenMark-vararea.stan` (the same model with an area offset per grid, which `DenMark()` fits with a tissue boundary) and `DenMark-vararea-umi.stan` (area offset per grid and total-UMI expression offset, used by the revision scripts). As in the revised manuscript, there is one model: the simplified model with independent fields (formerly M1), the M1/M2 labels and WAIC were removed in 0.9.0.
- `inst/extdata/denmark_toy.csv.gz`: the toy data. `data-raw/make_toy_data.R` rebuilds it, together with the parity fixtures `python/tests/data/toy_stan_data.json` (Tumor cells, no boundary) and `python/tests/data/toy_stan_data_boundary.json` (with the alpha-shape boundary of all cells: alpha, tissue area, grid areas and types, Stan data); copy the new file to `python/denmark/data/` afterwards, since a Python test checks that the two are identical.
- `tests/testthat/`: R tests.
- `vignettes/DenMark.Rmd`: the Get started vignette. `vignettes/articles/`: articles for the website only (Python, toy data, case studies, citation), excluded from the R build.
- `python/denmark/core.py`: the Python implementation, a function-by-function mirror of the R code, including the tissue boundary (`alpha_shape()`, `TissueBoundary`, `plot_boundary()`). `datasets.py`: `load_toy()`. `stan/`: byte copies of `inst/stan/*.stan`. `data/`: the toy data.
- `python/tests/`: pytest suite.
- `CodeInPaper/`: the scripts behind the manuscript's analyses, one subfolder per case study. Each script starts with a CONFIGURE ME block that sets `DATA_DIR`, `OUTPUT_DIR` and `REPO_ROOT`.
- `Tutorial_DenMark.Rmd` / `.pdf`: Mingchi Xu's tutorial for the script version of the R code (`DenMark/DenMark-functions.R`, which this release replaces with `R/`; it remains in the git history, last at commit 4c68cc9). It also fits the simplified model and compares models by WAIC, which the packages removed, so it does not run with the packages as written; they reproduce its tissue-boundary numbers (see below).

## Conventions

- Grid order: x varies fastest. In R, `grid_index = i + (j - 1) * nx` (1-based), so the grid vectors equal `as.vector(Ymat)`, whose rows index x. In Python, `grid_index = i + j * nx` (0-based) and the vectors equal `Ymat.ravel(order="F")`. A cell on an interior grid edge goes to the upper or right grid; cells on the outer edge of the region are kept.
- Tissue boundary: compute it from all cells of the section (`denmark_alpha_shape(x, y)`, `alpha_shape(x, y)`), not from the cell type analyzed. The region where a cell type could have been observed is the whole tissue; a boundary drawn around the analyzed cells alone, or no boundary, turns the empty grids outside the tissue into apparent regions of low density. The default alpha is 3 times the median Delaunay edge length (`alpha_mult = 3`).
- With `boundary =` (`buildgridpp()`, `denmark_prepare()`, `DenMark()`, and `denmark_scan()` through `...`; in Python `buildgridpp()`, `prepare_data()`, `DenMark()`, `scan_genes()`), the grid spans the extent of the tissue polygon and every grid cell is `interior`, `boundary` or `exterior`. Exterior grid cells are dropped (with `min_area_frac > 0`, also boundary grid cells with less than that fraction of a full grid cell inside the tissue), and every remaining grid cell enters with its area inside the tissue.
- `log_grid_area` in the Stan data is a scalar without a boundary (the whole bounding box, one common grid area) and a vector of length `n` with one. The variable-area model `DenMark-vararea.stan` is then selected automatically (`vararea` in `denmark_prepare()$info` and `denmark_model(vararea = TRUE)` in R). A custom `stan_file` used with a boundary must declare `vector[n] log_grid_area`; R stops otherwise.
- HSGP coordinates are the centroids of the grid cells that enter the model, centered at the center of their range, and `L1` and `L2` are the half ranges of those grid cells. Without a boundary they span the bounding box (by default the range of the data); with one, the grid spans the extent of the tissue polygon and only the retained grid cells count. The HSGP domain is `L = boundfactor * c(L1, L2)`.
- `rho = Astar[2,1]`, an entry of the Cholesky factor of the latent correlation matrix. The coregionalization coefficients are `a11 = sigma1`, `a21 = sigma2 * rho` and `a22 = sigma2 * sqrt(1 - rho^2)` (`denmark_draws()`, `DenMarkResult.draws()`).
- Symbols of the manuscript and their names in the code:

  | Manuscript | Code |
  |---|---|
  | intercepts β1 (cell density) and β2 (expression) | `beta0`, `beta1` |
  | latent fields ω1, ω2 | `w1`, `w2` |
  | HSGP weights η | `betab1`, `betab2` |
  | length-scales ℓ1, ℓ2 | `ell[1]`, `ell[2]` |
  | marginal standard deviations σ1, σ2 | `sigma[1]`, `sigma[2]` |
  | density-expression correlation ρ | `Astar[2,1]` |
  | cell count N_i and summed expression M_i of grid i | `Y`, `M` (`y1`, `y2` in the Stan data) |
  | grid area \|G\| (with a boundary, its part inside the tissue) | `grid_area`, or the `area` column of the grid with a boundary (`log_grid_area` in the Stan data) |
  | boundary factor c | `boundfactor` |
  | B1 x B2 basis functions | `basis` (`mstar` = B1 B2) |

  The intercepts are numbered from 1 in the manuscript and from 0 in the code.
- The expression model uses the offset `log(Y + 0.05)`, so grids without cells stay in the likelihood; there is no indicator for empty grids. The same constant appears in `denmark_init()`, `denmark_fields()` and their Python counterparts.
- `beta0` starts at the log of the number of cells per unit of tissue area (`denmark_init()`, `make_init()`). Without a boundary this is the earlier starting value, `log(mean(Y)) - log_grid_area`, unless the grids hold fewer than 0.1 cells on average.
- The Stan data field `is_centerted_PHI` is misspelled. Keep the name: the Stan models and earlier fits use it.
- `inst/stan/` is the single source of the Stan models. `python/denmark/stan/` holds byte copies, and a Python test fails when they differ: edit the file in `inst/stan/`, then copy it.
- R and Python build the same Stan data: grid indices, types, counts and coordinates exactly, and, with a tissue boundary, the grid areas up to the last bit, because sf computes polygon areas with GDAL's formula and shapely with GEOS's. `python/tests` compares the Stan data that Python builds from the toy data with the R output in `python/tests/data/toy_stan_data.json`, and the alpha shape, grid areas, grid types and Stan data with the boundary with `python/tests/data/toy_stan_data_boundary.json`; regenerate both fixtures with `data-raw/make_toy_data.R` when the R side changes on purpose.
- The toy data are a simulated tissue section: 7,449 cells with a `celltype` column, 3,449 `Tumor` cells (the analyzed population; true rho 0.6, 0 and -0.6 for geneA, geneB and geneC) and 4,000 `Other` cells that only delineate the tissue. Examples compute the boundary from all cells and analyze `tumor <- subset(toy, celltype == "Tumor", select = -celltype)` (Python: `toy[toy["celltype"] == "Tumor"].drop(columns="celltype")`). `denmark_scan()` and `scan_genes()` treat every column other than `x` and `y` as a gene, so drop `celltype` before a scan.
- In R, `matern32(d, sigma, rho)` names the length-scale `rho` for backward compatibility; Python's `matern32(d, sigma, ell)` calls it `ell`. Neither argument is the density-expression correlation.
- Python's `simulate()` draws with numpy's generator and does not reproduce R's `denmark_simulate()` for the same seed, which is why the toy data are shared as a file.

## Do not change the results silently

Until the manuscript is accepted, change packaging, documentation and tests, not the numerical internals. The function bodies in `R/` come from Mingchi Xu's revision of `DenMark/DenMark-functions.R` (GitHub commit 4c68cc9, 2026-10-24), apart from the Stan-file lookup, the compile cache and the removal of the simplified model and WAIC; the Stan code of the model is unchanged. `python/denmark/core.py` carries over `DenMark/denmark.py` of commit e8b4cc9 with the same exceptions, and its tissue-boundary code is a port of his R code, checked by parity tests against the R fixture `python/tests/data/toy_stan_data_boundary.json`. The planned refactors wait for v1.1 or v2.0: merging the three Stan variants into one model with the area and expression offsets passed as data, having `CodeInPaper/` call the packages, and adding names aligned with the manuscript while keeping the current ones as aliases.

Without a boundary, the revision leaves the Stan data and the starting values as they were, up to floating-point rounding: the new centering arithmetic of `denmark_prepare()` moves `coords` and `L` in the last bits (at most 9e-16 on the toy data; the bound scales with the coordinates), and the new `beta0` start agrees with the old one up to the last bit unless the grids hold fewer than 0.1 cells on average. On the toy data, the Stan data file that cmdstanr writes (15 significant digits) is byte-identical before and after; cmdstanpy writes full precision, so its file can differ in the last digits. Simulated data did change: `denmark_simulate()` no longer reproduces the data of the earlier scripts for the same seed (see `NEWS.md`).

A change to `R/`, `python/denmark/` or `inst/stan/` that alters a numerical output changes what the manuscript reports. Before such a change, record these outputs on the toy data, in R and in Python, for the Tumor cells (`tumor`, above) and the tissue boundary of all cells (`tissue <- denmark_alpha_shape(toy$x, toy$y)`; `tissue = dm.alpha_shape(toy["x"], toy["y"])`); afterwards, show that they are identical, or explain each difference in the commit message:

1. the alpha and the tissue area: `attr(tissue, "alpha")` and `sf::st_area(tissue)`, `tissue.alpha` and `tissue.area` in Python;
2. `buildgridpp(10, tumor, gene = "geneA", boundary = tissue)`: the counts, areas and types of the grid cells;
3. `denmark_prepare(tumor, gene = "geneA", grid_size = 10, basis = 10, boundary = tissue)$stan_data`, and `prepare_data(tumor, gene="geneA", grid_size=10, basis=10, boundary=tissue)["stan_data"]` in Python, and the same without `boundary`;
4. a short seeded fit with `boundary = tissue, grid_size = 10, basis = 10, chains = 2, iter_warmup = 300, iter_sampling = 300, seed = 2024`: the draws of the key parameters, `Astar[2,1]` among them.

Grid counts, grid indices, grid types and the Stan data must match exactly. Draws must be identical for the same CmdStan version and seed; record the CmdStan version with them, since another CmdStan or compiler version can change the draws.

Reference values, recorded on 2026-10-24 with R 4.6.1 and CmdStan 2.40.0 on Windows 11 (Rtools45):

| Output | Value |
|---|---|
| alpha; tissue area | 0.1068794; 8.5385318 (the bounding box of all cells has area 12.996216) |
| 10 x 10 grid with the boundary | 46 interior, 36 boundary and 18 exterior grid cells; `n` = 82, with all 3,449 Tumor cells |
| `L1`, `L2` with the boundary | 2.550369, 1.2739545 |
| Fit of geneA (item 4) | posterior mean of rho 0.648, 95% interval (0.387, 0.836); maximum R-hat 1.026, no divergent transitions, 278 of 600 transitions at the maximum tree depth |
| `denmark_scan()` of the Tumor cells with the same settings, rho (95% interval) | geneA 0.648 (0.387, 0.836), PCG; geneB 0.113 (-0.191, 0.408), not DCG; geneC -0.599 (-0.821, -0.315), NCG |

The same code reproduces the boundary numbers printed in `Tutorial_DenMark.pdf`, which simulates the toy data without rounding the coordinates: alpha 0.1068781, tissue area 8.5385318 and, on the 15 x 15 grid, 114 interior, 56 boundary and 55 exterior grid cells, `n` = 170, `L1` = 2.550369 and `L2` = 1.273955.

## Releases

Bump `Version` in `DESCRIPTION`, `__version__` in `python/denmark/__init__.py`, `version` in `CITATION.cff` and the heading of `NEWS.md` together. Publishing a GitHub release runs `.github/workflows/publish.yml`, which checks that the tag equals `v` followed by the DESCRIPTION version and by `__version__`, builds the Python package from `python/` and uploads it to PyPI by trusted publishing. A version number can be uploaded only once. v1.0.0 is reserved for the release that accompanies the manuscript.

## Style

- R documentation is roxygen2 with markdown. After editing a roxygen block, run `devtools::document()` and keep the regenerated `NAMESPACE` and `man/` in the change; do not edit them by hand.
- Python docstrings use the `:param:` style.
- No em-dashes in prose, comments or documentation: use a colon, semicolon, comma or parentheses. American spelling in new text. Comments inside function bodies that were copied verbatim keep their original wording.
- Examples in the README and the docs run on the toy data: the tissue boundary of all cells, the Tumor cells without the `celltype` column, and `grid_size = 10, basis = 10, chains = 2, seed = 2024`. The README and the vignettes use `iter_warmup = 500, iter_sampling = 500`, at which all three toy genes converge (maximum R-hat below 1.01); the man-page examples, the tests and the recorded outputs above use 300 and 300 to stay fast.
