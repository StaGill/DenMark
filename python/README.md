# denmark: DenMark for Python

`denmark` is the Python implementation of DenMark, a density-dependent marked point process framework for single-cell resolution spatial transcriptomics (MERFISH, Xenium, seqFISH). DenMark quantifies how the expression of a gene varies with local cell density and identifies density-correlated genes (DCGs). The package fits the same Stan model as the DenMark R package through [cmdstanpy](https://mc-stan.org/cmdstanpy/) and builds the same Stan data from the same input, so R and Python results are directly comparable.

The method is described in Xu, Schmidt and Zhang, under review at *Cell Reports Methods*.

## Installation

```bash
pip install "git+https://github.com/StaGill/DenMark#subdirectory=python"
python -m cmdstanpy.install_cmdstan    # once; on Windows add --compiler for the C++ toolchain
```

The package needs Python 3.9 or later. The plotting functions also need matplotlib, and the tissue boundary needs shapely; the `plot` and `boundary` extras install them:

```bash
pip install "denmark[plot,boundary] @ git+https://github.com/StaGill/DenMark#subdirectory=python"
```

The model compiles once, on first use, into `~/.cache/denmark`; set the environment variable `DENMARK_CACHE_DIR` to use another folder.

## Quick start

The toy data are a simulated tissue section: 3,449 cells of the analyzed type ("Tumor") and 4,000 cells of other types ("Other"). The true density-expression correlations of the Tumor cells are rho = 0.6 for geneA, 0 for geneB and -0.6 for geneC. The tissue boundary is the alpha shape of all cells; the model is then fitted to the Tumor cells inside it.

```python
import denmark as dm

toy = dm.datasets.load_toy()       # one row per cell: x, y, celltype, geneA, geneB, geneC
tissue = dm.alpha_shape(toy["x"], toy["y"])
tumor = toy[toy["celltype"] == "Tumor"].drop(columns="celltype")
prep = dm.prepare_data(tumor, gene="geneA", grid_size=10, basis=10, boundary=tissue)
print(prep["grid"].head())         # cell count Y, summed expression M and area inside the tissue
```

Fit DenMark to geneA. These short settings suit a first look; an analysis should keep the defaults (4 chains, 1,000 warm-up and 1,000 sampling iterations).

```python
fit = dm.DenMark(tumor, gene="geneA", grid_size=10, basis=10, boundary=tissue,
                 chains=2, iter_warmup=500, iter_sampling=500, seed=2024)
print(fit)                         # rho with its 95% credible interval and sampler diagnostics
```

Your own data take the same form, one row per cell:

```python
import pandas as pd

cells = pd.read_csv("cells.csv")   # columns x, y, celltype and one count column per gene
tissue = dm.alpha_shape(cells["x"], cells["y"])
focal = cells[cells["celltype"] == "Tumor"].drop(columns="celltype")
fit = dm.DenMark(focal, gene="CDH1", grid_size=15, basis=25, boundary=tissue)
```

## Learn more

The [main README](https://github.com/StaGill/DenMark#readme) walks through the full workflow in R and Python: latent fields, correlation decay and gene scans. The documentation is at <https://stagill.github.io/DenMark/>.

## Citation and license

Please cite the DenMark paper when you use the package; the reference is in [CITATION.cff](https://github.com/StaGill/DenMark/blob/main/CITATION.cff). The package is released under the MIT license.
