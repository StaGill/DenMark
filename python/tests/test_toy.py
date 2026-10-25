"""The bundled toy data."""

from importlib import resources

import numpy as np
import pytest

import denmark as dm

from .conftest import REPO_ROOT


def test_load_toy(toy):
    assert list(toy.columns) == ["x", "y", "celltype", "geneA", "geneB", "geneC"]
    assert toy["celltype"].value_counts().to_dict() == {"Other": 4000, "Tumor": 3449}
    assert not toy.isna().any().any()
    assert toy["x"].dtype == np.float64 and toy["y"].dtype == np.float64
    assert toy["x"].between(0, 6).all() and toy["y"].between(0, 4).all()
    for gene in ("geneA", "geneB", "geneC"):
        assert np.issubdtype(toy[gene].dtype, np.integer) and (toy[gene] >= 0).all(), gene


def test_toy_grids_without_losing_cells(tumor):
    sd = dm.prepare_data(tumor, gene="geneA", grid_size=10, basis=10)["stan_data"]
    assert sd["n"] == 100 and sd["mstar"] == 100
    assert sd["y1"].sum() == len(tumor)
    assert sd["y2"].sum() == tumor["geneA"].sum()


def test_toy_tissue_holds_every_cell(toy, tumor, tissue):
    """The alpha shape of all cells covers every cell, so no Tumor cell is lost."""
    import shapely

    assert shapely.intersects_xy(tissue.polygon, toy["x"], toy["y"]).all()
    prep = dm.prepare_data(tumor, gene="geneA", grid_size=10, basis=10, boundary=tissue)
    assert prep["info"]["n_cells"] == len(tumor)
    assert prep["stan_data"]["y2"].sum() == tumor["geneA"].sum()


def test_toy_file_matches_the_r_package():
    py = resources.files("denmark").joinpath("data/denmark_toy.csv.gz")
    r = REPO_ROOT / "inst" / "extdata" / "denmark_toy.csv.gz"
    if not (py.is_file() and r.is_file()):
        pytest.skip("the toy data file is missing from denmark/data/ or inst/extdata/")
    assert py.read_bytes() == r.read_bytes()
