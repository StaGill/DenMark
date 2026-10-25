"""Toy data bundled with denmark for examples and tests."""

from importlib import resources

import pandas as pd


def load_toy():
    """Return the toy dataset: a simulated tissue section with 3,449 cells of
    the analyzed type ("Tumor") and 4,000 cells of other types ("Other"), one
    row per cell, with coordinates x and y, the cell type and the counts of
    three genes.

    The cells were simulated from the DenMark data-generating process in an
    irregular, lobed tissue with ``denmark_simulate()`` of the DenMark R
    package; ``data-raw/make_toy_data.R`` in the repository rebuilds the file.
    The true density-expression correlations of the Tumor cells are rho = 0.6
    for geneA, 0 for geneB and -0.6 for geneC, so a DenMark scan should call
    geneA a positively density-correlated gene (PCG), geneC a negatively
    density-correlated gene (NCG) and geneB neither. The other cells only help
    delineate the tissue: compute the boundary with ``alpha_shape()`` from all
    cells, then analyze the Tumor cells. The R package ships the same file,
    read by ``denmark_toy_data()``. Numbers are parsed to the nearest double,
    as R parses them: with ``float_precision="round_trip"`` before pandas 3,
    which made that the default and deprecated the argument.

    :return: a pandas DataFrame with columns x, y, celltype, geneA, geneB and
        geneC.
    """
    kwargs = {} if int(pd.__version__.split(".")[0]) >= 3 else {"float_precision": "round_trip"}
    ref = resources.files(__package__).joinpath("data/denmark_toy.csv.gz")
    with resources.as_file(ref) as path:
        return pd.read_csv(path, **kwargs)
