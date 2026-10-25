"""The tissue boundary: alpha shape, grid cells inside the tissue and the Stan
data of the variable-area model."""

import dataclasses
import sys

import numpy as np
import pandas as pd
import pytest

import denmark as dm

shapely = pytest.importorskip("shapely")


def lattice(nx, ny, x0=0.0, y0=0.0, hole=None):
    """Points on a unit lattice, optionally without those inside the open
    square `hole` = (lo, hi)."""
    x, y = np.meshgrid(x0 + np.arange(nx, dtype=float), y0 + np.arange(ny, dtype=float))
    x, y = x.ravel(), y.ravel()
    if hole is not None:
        out = (x > hole[0]) & (x < hole[1]) & (y > hole[0]) & (y < hole[1])
        x, y = x[~out], y[~out]
    return x, y


@pytest.fixture
def cells():
    """The cells of test_grid.py on [0, 4] x [0, 2], and a triangular tissue
    whose long side runs through the shared corner (2, 1) of the 2 x 2 grid:
    grid areas 2, 1, 1 and 0."""
    return pd.DataFrame({"x": [0.0, 0.5, 2.5, 3.5, 2.0, 4.0, 1.0],
                         "y": [0.0, 0.2, 0.1, 0.5, 1.0, 2.0, 1.9],
                         "g": [1, 2, 0, 2, 3, 4, 5]})


@pytest.fixture
def triangle():
    return shapely.Polygon([(0, 0), (4, 0), (0, 2)])


def test_alpha_shape_of_the_toy_tissue(toy, tissue):
    assert isinstance(tissue, dm.TissueBoundary)
    poly = tissue.polygon
    assert poly.geom_type == "Polygon" and poly.is_valid and not poly.interiors
    assert tissue.area == poly.area
    xmin, ymin, xmax, ymax = poly.bounds
    assert toy["x"].min() <= xmin < xmax <= toy["x"].max()
    assert toy["y"].min() <= ymin < ymax <= toy["y"].max()
    # a given alpha reproduces the default, alpha_mult scales it, duplicates are ignored
    assert dm.alpha_shape(toy["x"], toy["y"], alpha=tissue.alpha).polygon.equals_exact(poly, 0)
    assert dm.alpha_shape(toy["x"], toy["y"], alpha_mult=6).alpha == pytest.approx(2 * tissue.alpha)
    dup = pd.concat([toy, toy.head(50)])
    assert dm.alpha_shape(dup["x"], dup["y"]).polygon.equals_exact(poly, 0)


def test_tissue_boundary_is_frozen(tissue):
    with pytest.raises(dataclasses.FrozenInstanceError):
        tissue.alpha = 1.0


def test_alpha_keeps_triangles_below_the_threshold():
    """On a unit lattice every Delaunay triangle has circumradius sqrt(2) / 2."""
    x, y = lattice(5, 4)
    assert dm.alpha_shape(x, y, alpha=0.8).area == pytest.approx(12.0)
    with pytest.raises(ValueError, match="increase `alpha`"):
        dm.alpha_shape(x, y, alpha=0.7)


def test_largest_piece_or_all_pieces():
    xa, ya = lattice(5, 5)                      # area 16
    xb, yb = lattice(3, 3, x0=10.0)             # area 4
    x, y = np.r_[xa, xb], np.r_[ya, yb]
    largest = dm.alpha_shape(x, y, alpha=1.0)
    assert largest.polygon.geom_type == "Polygon" and largest.area == pytest.approx(16.0)
    both = dm.alpha_shape(x, y, alpha=1.0, keep="all")
    assert both.polygon.geom_type == "MultiPolygon" and both.area == pytest.approx(20.0)
    with pytest.raises(ValueError, match="keep"):
        dm.alpha_shape(x, y, keep="first")


def test_holes_are_filled_unless_asked_otherwise():
    # a 9 x 9 square with a 5 x 5 hole, whose four corners keep a triangle of area 1/2
    x, y = lattice(10, 10, hole=(2, 7))
    filled = dm.alpha_shape(x, y, alpha=1.0)
    assert not filled.polygon.interiors and filled.area == pytest.approx(81.0)
    holed = dm.alpha_shape(x, y, alpha=1.0, fill_holes=False)
    assert len(holed.polygon.interiors) == 1 and holed.area == pytest.approx(81 - 25 + 2)


def test_alpha_shape_needs_three_distinct_cells():
    with pytest.raises(ValueError, match="3 distinct"):
        dm.alpha_shape([0, 1, 1, 0], [0, 1, 1, 0])


def test_shapely_is_needed(monkeypatch):
    monkeypatch.setitem(sys.modules, "shapely", None)     # import shapely now fails
    with pytest.raises(ImportError, match=r"alpha_shape\(\) needs the 'shapely' package"):
        dm.alpha_shape([0, 1, 0], [0, 0, 1])


def test_prepare_data_drops_the_grid_cells_outside_the_tissue(cells, triangle, capsys):
    out = dm.prepare_data(cells, gene="g", grid_size=2, basis=3, boundary=triangle)
    assert "2 cell(s) fall in grid cells outside the tissue" in capsys.readouterr().out
    sd, grid, info = out["stan_data"], out["grid"], out["info"]
    np.testing.assert_array_equal(grid["grid_index"], [0, 1, 2])
    np.testing.assert_array_equal(sd["y1"], [2, 2, 1])
    np.testing.assert_array_equal(sd["y2"], [3, 2, 5])
    assert sd["n"] == 3 and sd["log_grid_area"].shape == (3,)
    np.testing.assert_allclose(np.exp(sd["log_grid_area"]), [2.0, 1.0, 1.0])
    assert info["vararea"] is True and info["boundary"] is triangle
    assert (info["n_grid_total"], info["n_interior"], info["n_boundary"]) == (4, 1, 2)
    assert info["n_cells"] == 5
    # the retained grid sets the domain: centroids (1, 0.5), (3, 0.5), (1, 1.5)
    np.testing.assert_allclose(info["centre"], [2.0, 1.0])
    np.testing.assert_allclose(sd["coords"], [[-1.0, -0.5], [1.0, -0.5], [-1.0, 0.5]])
    np.testing.assert_allclose(sd["L"], [1.5 * 2.0, 1.5 * 1.0])
    np.testing.assert_allclose(grid[["cx", "cy"]].to_numpy(), sd["coords"])


def test_min_area_frac_drops_small_boundary_cells(cells, triangle, capsys):
    keep_all = dm.prepare_data(cells, gene="g", grid_size=2, basis=3, boundary=triangle,
                               min_area_frac=0.4)
    assert keep_all["stan_data"]["n"] == 3
    out = dm.prepare_data(cells, gene="g", grid_size=2, basis=3, boundary=triangle,
                          min_area_frac=0.6)
    assert "5 cell(s) fall in grid cells outside the tissue" in capsys.readouterr().out
    assert out["stan_data"]["n"] == 1 and list(out["grid"]["type"]) == ["interior"]
    # the domain follows the one grid cell left, not the bounding box
    np.testing.assert_allclose(out["info"]["centre"], [1.0, 0.5])
    np.testing.assert_allclose([out["info"]["L1"], out["info"]["L2"]], [1.0, 0.5])


def test_toy_grid_inside_the_tissue(tumor, tissue):
    out = dm.prepare_data(tumor, gene="geneA", grid_size=10, basis=10, boundary=tissue)
    grid, info = out["grid"], out["info"]
    assert (info["n_grid_total"], info["n_interior"], info["n_boundary"]) == (100, 46, 36)
    assert out["stan_data"]["n"] == len(grid) == 82
    full = info["grid_area"]
    assert np.allclose(grid.loc[grid["type"] == "interior", "area"], full)
    assert ((grid.loc[grid["type"] == "boundary", "area"] > 0)
            & (grid.loc[grid["type"] == "boundary", "area"] < full)).all()
    # the grid spans the tissue, so the grid areas add up to the tissue area
    assert grid["area"].sum() == pytest.approx(tissue.area, rel=1e-6)


def test_make_init_uses_the_tissue_area(tumor, tissue):
    sd = dm.prepare_data(tumor, gene="geneA", grid_size=10, basis=10, boundary=tissue)["stan_data"]
    init = dm.make_init(sd)
    assert init["beta0"] == pytest.approx(np.log(len(tumor) / np.exp(sd["log_grid_area"]).sum()))


def test_cells_outside_the_boundary_extent_are_reported(cells, triangle, capsys):
    small = shapely.Polygon([(0, 0), (2, 0), (0, 1)])         # extent [0, 2] x [0, 1]
    dm.buildgridpp(2, cells, gene="g", boundary=small)
    assert "4 cell(s) lie outside the extent of the tissue boundary" in capsys.readouterr().out
    far = cells.assign(x=cells["x"] + 100)
    with pytest.raises(ValueError, match="same units"):
        dm.buildgridpp(2, far, gene="g", boundary=triangle)


def test_a_non_numeric_column_is_refused_as_the_expression():
    df = pd.DataFrame({"x": [1.0, 2.0, 3.0], "y": [1.0, 2.0, 3.0], "celltype": "T", "g": [1, 2, 3]})
    with pytest.raises(ValueError, match="Column 'celltype' is not numeric"):
        dm.buildgridpp(3, df)
    assert dm.buildgridpp(3, df, gene="g")["grid"]["M"].sum() == 6
