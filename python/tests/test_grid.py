"""Gridding: cell counts, summed expression, grid order, the bounding box and
the grid areas inside a tissue."""

import numpy as np
import pandas as pd
import pytest

import denmark as dm


@pytest.fixture
def cells():
    """Seven cells on [0, 4] x [0, 2]. On a 2 x 2 grid each grid cell is 2 wide
    and 1 high; (2, 1) sits on two interior edges and (4, 2) on the outer corner."""
    return pd.DataFrame({"x": [0.0, 0.5, 2.5, 3.5, 2.0, 4.0, 1.0],
                         "y": [0.0, 0.2, 0.1, 0.5, 1.0, 2.0, 1.9],
                         "g": [1, 2, 0, 2, 3, 4, 5]})


def test_counts_and_order(cells):
    g = dm.buildgridpp(2, cells, gene="g")
    grid = g["grid"]
    assert list(grid.columns) == ["grid_index", "i", "j", "x", "y", "Y", "M", "area", "type"]
    # x varies fastest: grid_index = i + j * nx (0-based; i + (j - 1) * nx in R)
    np.testing.assert_array_equal(grid["i"], [0, 1, 0, 1])
    np.testing.assert_array_equal(grid["j"], [0, 0, 1, 1])
    np.testing.assert_array_equal(grid["grid_index"], grid["i"] + grid["j"] * 2)
    # a cell on an interior edge goes to the upper/right grid cell; the outer boundary is kept
    np.testing.assert_array_equal(grid["Y"], [2, 2, 1, 2])
    np.testing.assert_array_equal(grid["M"], [3, 2, 5, 7])
    np.testing.assert_array_equal(g["loc_mark_mat"], np.column_stack([grid["Y"], grid["M"]]))
    np.testing.assert_array_equal(grid["x"], [1.0, 3.0, 1.0, 3.0])
    np.testing.assert_array_equal(grid["y"], [0.5, 0.5, 1.5, 1.5])
    assert g["grid_area"] == 2.0 and g["step"] == (2.0, 1.0) and g["grid_size"] == (2, 2)
    assert g["gene"] == "g"
    # without a tissue boundary every grid cell is interior, with the full area
    assert (grid["area"] == 2.0).all() and (grid["type"] == "interior").all()
    assert g["boundary"] is None


def test_area_and_type_inside_a_tissue(cells):
    """The tissue is the triangle (0, 0), (4, 0), (0, 2), whose long side runs
    through the shared corner (2, 1) of the four grid cells: the lower-left
    cell lies inside, the upper-right one only touches the tissue, and the
    tissue covers half of the other two."""
    shapely = pytest.importorskip("shapely")
    tri = shapely.Polygon([(0, 0), (4, 0), (0, 2)])
    g = dm.buildgridpp(2, cells, gene="g", boundary=tri)
    np.testing.assert_array_equal(g["bbox"], [0, 4, 0, 2])          # the extent of the tissue
    np.testing.assert_allclose(g["grid"]["area"], [2.0, 1.0, 1.0, 0.0])
    assert list(g["grid"]["type"]) == ["interior", "boundary", "boundary", "exterior"]
    np.testing.assert_array_equal(g["grid"]["Y"], [2, 2, 1, 2])       # counts do not change
    assert g["boundary"] is tri
    # a TissueBoundary gives the same grid, and the tissue sets the bounding box
    tb = dm.TissueBoundary(polygon=shapely.box(0.5, 0.0, 3.5, 2.0), alpha=1.0)
    g = dm.buildgridpp(2, cells, gene="g", boundary=tb)
    np.testing.assert_array_equal(g["bbox"], [0.5, 3.5, 0, 2])
    assert (g["grid"]["type"] == "interior").all() and g["grid"]["Y"].sum() == 5


def test_boundary_must_be_a_polygon(cells):
    with pytest.raises(TypeError, match="boundary"):
        dm.buildgridpp(2, cells, gene="g", boundary=[(0, 0), (4, 0), (0, 2)])


def test_matrices_have_rows_x(cells):
    g = dm.buildgridpp((4, 2), cells, gene="g")
    assert len(g["grid"]) == 8 and g["Ymat"].shape == (4, 2)
    np.testing.assert_array_equal(g["Ymat"], [[2, 0], [0, 1], [1, 1], [1, 1]])
    np.testing.assert_array_equal(g["Ymat"].ravel(order="F"), g["grid"]["Y"])
    np.testing.assert_array_equal(g["Mmat"].ravel(order="F"), g["grid"]["M"])


def test_bbox_drops_cells_outside(cells):
    g = dm.buildgridpp(2, cells, gene="g", bbox=(0, 2, 0, 2))
    np.testing.assert_array_equal(g["bbox"], [0, 2, 0, 2])
    np.testing.assert_array_equal(g["grid"]["Y"], [2, 0, 0, 2])
    np.testing.assert_array_equal(g["grid"]["M"], [3, 0, 0, 8])
    np.testing.assert_array_equal(g["grid"]["x"], [0.5, 1.5, 0.5, 1.5])


def test_bbox_larger_than_the_data(cells):
    g = dm.buildgridpp(2, cells, gene="g", bbox=(-1, 5, -1, 3))
    assert g["grid"]["Y"].sum() == len(cells)
    np.testing.assert_array_equal(g["grid"]["x"], [0.5, 3.5, 0.5, 3.5])
    np.testing.assert_array_equal(g["grid"]["y"], [0.0, 0.0, 2.0, 2.0])


@pytest.mark.parametrize("bbox", [(1, 0, 0, 1), (0, 1, 1, 0), (0, 1, 0)])
def test_invalid_bbox(cells, bbox):
    with pytest.raises(ValueError, match="bbox"):
        dm.buildgridpp(2, cells, gene="g", bbox=bbox)


def test_bbox_without_cells(cells):
    with pytest.raises(ValueError, match="No cells"):
        dm.buildgridpp(2, cells, gene="g", bbox=(10, 11, 10, 11))


def test_input_checks(cells):
    with pytest.raises(ValueError, match="grid_size"):
        dm.buildgridpp(1, cells, gene="g")
    with pytest.raises(ValueError, match="not a column"):
        dm.buildgridpp(2, cells, gene="absent")
    with pytest.raises(ValueError, match="non-negative"):
        dm.buildgridpp(2, cells.assign(g=-cells["g"]), gene="g")
    assert dm.buildgridpp(2, cells)["gene"] == "g"


def test_non_integer_counts_are_rounded(cells):
    with pytest.warns(UserWarning, match="not integers"):
        g = dm.buildgridpp(2, cells.assign(g=[0.4, 0.4, 0, 0, 0, 0, 0]), gene="g")
    np.testing.assert_array_equal(g["grid"]["M"], [1, 0, 0, 0])


def test_grid_pts_centroid_matches_the_grid(cells):
    gp = dm.grid_pts_centroid((4, 2), cells)
    grid = dm.buildgridpp((4, 2), cells, gene="g")["grid"]
    assert gp["N"] == 8 and gp["grid_area"] == 1.0
    np.testing.assert_array_equal(gp["grid_points"]["x"], grid["x"])
    np.testing.assert_array_equal(gp["grid_points"]["y"], grid["y"])
    gp = dm.grid_pts_centroid(2, cells, bbox=(0, 2, 0, 2))
    np.testing.assert_array_equal(gp["grid_points"]["x"], [0.5, 1.5, 0.5, 1.5])
    np.testing.assert_array_equal(gp["grid_points"]["y"], [0.5, 0.5, 1.5, 1.5])


def test_hsgp_indices_first_index_fastest():
    S = dm.hsgp_indices(3)
    assert S.shape == (9, 2) and np.issubdtype(S.dtype, np.integer)
    np.testing.assert_array_equal(S[:4], [[1, 1], [2, 1], [3, 1], [1, 2]])
    np.testing.assert_array_equal(dm.hsgp_indices((2, 3)),
                                  [[1, 1], [2, 1], [1, 2], [2, 2], [1, 3], [2, 3]])
