"""R and Python build the same Stan data from the toy data.

The fixtures in tests/data/ are written by data-raw/make_toy_data.R (matrices
row-major, about 15 significant digits), for the Tumor cells of the toy data:
toy_stan_data.json holds R's denmark_prepare(tumor, gene = "geneA",
grid_size = 10, basis = 10, boundfactor = 1.5)$stan_data, and
toy_stan_data_boundary.json the same call with the alpha-shape tissue boundary
of all cells, denmark_alpha_shape(toy$x, toy$y), together with its alpha, the
tissue area and the area, type and index of each grid cell kept.
toy_stan_data_umi.json holds the first call with offset = "umi", for a
per-cell total count of 20 plus the three genes, and toy_visium_units.json a
set of Visium-like spots with R's seven-spot pools (denmark_pool_visium()) and
the Stan data of the pools, prepared with units = TRUE and offset = "umi".

The alpha shape, the grid types, the grid coordinates and the pools agree
exactly. The areas agree to about 1e-15: sf computes planar areas with GDAL
and shapely with GEOS, whose formulas differ in the last bits. The centroids
of the pools, in micrometres, agree to within 1e-12, the precision with which
the fixture stores them.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import denmark as dm

DATA = Path(__file__).resolve().parent / "data"
INTEGERS = ("n", "y1", "y2", "d", "mstar", "indices", "is_centerted_PHI")
FLOATS = ("coords", "log_grid_area", "mrange", "L")
UMI_FLOATS = FLOATS + ("log_umi",)


def _fixture(name):
    path = DATA / name
    if not path.exists():
        pytest.skip(f"tests/data/{name} is missing; run data-raw/make_toy_data.R")
    return json.loads(path.read_text(encoding="utf-8"))


def _assert_same_stan_data(py, r, rtol, floats=FLOATS, atol=0):
    assert set(py) == set(r) == set(INTEGERS + floats)
    for key in INTEGERS:
        np.testing.assert_array_equal(np.asarray(py[key]), np.asarray(r[key]), err_msg=key)
    for key in floats:
        np.testing.assert_allclose(np.asarray(py[key], float), np.asarray(r[key], float),
                                   rtol=rtol, atol=atol, err_msg=key)


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


def test_stan_data_with_umi_offset_matches_r(tumor):
    r = _fixture("toy_stan_data_umi.json")
    cells = tumor.assign(UMI=20 + tumor["geneA"] + tumor["geneB"] + tumor["geneC"])
    py = dm.prepare_data(cells, gene="geneA", grid_size=10, basis=10, boundfactor=1.5,
                         offset="umi", umi="UMI")
    _assert_same_stan_data(py["stan_data"], r, rtol=1e-12, floats=UMI_FLOATS)


def test_pooled_visium_units_match_r():
    r = _fixture("toy_visium_units.json")
    units = dm.pool_visium(pd.DataFrame(r["spots"]), genes="geneA")
    np.testing.assert_array_equal(units.attrs["spot_unit"], r["spot_unit"])
    np.testing.assert_array_equal(units["n_spots"], r["n_spots"])
    py = dm.prepare_data(units, gene="geneA", basis=6, boundfactor=1.5, units=True,
                         offset="umi", umi="UMI")
    _assert_same_stan_data(py["stan_data"], r["stan_data"], rtol=1e-12, floats=UMI_FLOATS,
                           atol=1e-9)
