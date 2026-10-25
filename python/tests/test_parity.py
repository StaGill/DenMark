"""R and Python build the same Stan data from the toy data.

The fixtures in tests/data/ are written by data-raw/make_toy_data.R (matrices
row-major, about 15 significant digits), for the Tumor cells of the toy data:
toy_stan_data.json holds R's denmark_prepare(tumor, gene = "geneA",
grid_size = 10, basis = 10, boundfactor = 1.5)$stan_data, and
toy_stan_data_boundary.json the same call with the alpha-shape tissue boundary
of all cells, denmark_alpha_shape(toy$x, toy$y), together with its alpha, the
tissue area and the area, type and index of each grid cell kept.

The alpha shape, the grid types and the coordinates agree exactly. The areas
agree to about 1e-15: sf computes planar areas with GDAL and shapely with
GEOS, whose formulas differ in the last bits.
"""

import json
from pathlib import Path

import numpy as np
import pytest

import denmark as dm

DATA = Path(__file__).resolve().parent / "data"
INTEGERS = ("n", "y1", "y2", "d", "mstar", "indices", "is_centerted_PHI")
FLOATS = ("coords", "log_grid_area", "mrange", "L")


def _fixture(name):
    path = DATA / name
    if not path.exists():
        pytest.skip(f"tests/data/{name} is missing; run data-raw/make_toy_data.R")
    return json.loads(path.read_text(encoding="utf-8"))


def _assert_same_stan_data(py, r, rtol):
    assert set(py) == set(r) == set(INTEGERS + FLOATS)
    for key in INTEGERS:
        np.testing.assert_array_equal(np.asarray(py[key]), np.asarray(r[key]), err_msg=key)
    for key in FLOATS:
        np.testing.assert_allclose(np.asarray(py[key], float), np.asarray(r[key], float),
                                   rtol=rtol, err_msg=key)


def test_stan_data_matches_r(tumor):
    r = _fixture("toy_stan_data.json")
    py = dm.prepare_data(tumor, gene="geneA", grid_size=10, basis=10, boundfactor=1.5)
    _assert_same_stan_data(py["stan_data"], r, rtol=1e-12)


def test_stan_data_with_boundary_matches_r(tumor, tissue):
    r = _fixture("toy_stan_data_boundary.json")
    assert tissue.alpha == pytest.approx(r["alpha"], rel=1e-12)
    assert tissue.area == pytest.approx(r["tissue_area"], rel=1e-12)
    py = dm.prepare_data(tumor, gene="geneA", grid_size=10, basis=10, boundfactor=1.5,
                         boundary=tissue)
    grid = py["grid"]
    assert py["info"]["n_grid_total"] == r["n_grid_total"]
    np.testing.assert_array_equal(grid["grid_index"] + 1, r["grid_index"])   # R counts from 1
    np.testing.assert_array_equal(grid["type"], r["grid_type"])
    np.testing.assert_allclose(grid["area"], r["grid_area"], rtol=1e-9)
    _assert_same_stan_data(py["stan_data"], r["stan_data"], rtol=1e-9)
