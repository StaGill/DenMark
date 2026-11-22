"""Figures with a tissue boundary (matplotlib, Agg backend)."""

import numpy as np
import pandas as pd
import pytest

import denmark as dm

pytest.importorskip("matplotlib").use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


@pytest.fixture(scope="module")
def grid(tumor, tissue):
    return dm.buildgridpp(10, tumor, gene="geneA", boundary=tissue)


def test_plot_boundary(toy, tissue, grid):
    fig = dm.plot_boundary(toy, tissue, grid=grid, highlight=toy["celltype"] == "Tumor")
    left, right = fig.axes
    assert left.get_title() == f"Alpha-shape tissue boundary (alpha = {tissue.alpha:.3g})"
    assert right.get_title() == "Grid cells: 46 interior, 36 boundary, 18 exterior"
    assert [t.get_text() for t in left.get_legend().get_texts()] == ["other cells",
                                                                     "analyzed cell type"]
    assert len(left.lines) == 1 and len(right.lines) == 1      # the tissue outline
    fig.canvas.draw()


def test_plot_boundary_without_grid_or_alpha(toy, tissue):
    fig = dm.plot_boundary(toy, tissue.polygon)
    assert len(fig.axes) == 1
    assert fig.axes[0].get_title() == "Alpha-shape tissue boundary (alpha = NA)"
    fig.canvas.draw()


def test_tile_maps_skip_the_grid_cells_outside_the_tissue(grid):
    fig = dm.plot_grid(grid)
    for ax in fig.axes[:2]:
        z = np.ma.filled(ax.images[0].get_array().astype(float), np.nan)
        assert np.isfinite(z).sum() == (grid["grid"]["type"] != "exterior").sum() == 82
        assert len(ax.lines) == 1
    fig.canvas.draw()


def test_tile_maps_without_a_boundary(tumor):
    g = dm.buildgridpp(10, tumor, gene="geneA")
    fig = dm.plot_grid(g)
    z = fig.axes[0].images[0].get_array()
    np.testing.assert_array_equal(z, g["Ymat"].T)        # rows = y
    assert not fig.axes[0].lines


def test_maps_of_units_draw_one_point_per_unit():
    grid = pd.DataFrame({"grid_index": [0, 1, 2], "i": np.nan, "j": np.nan,
                         "x": [0.0, 2.5, 1.0], "y": [0.0, 0.9, 2.6], "Y": [3, 5, 2],
                         "M": [1, 2, 0], "type": "unit"})
    fig = dm.plot_grid({"grid": grid, "gene": "geneA"})
    for ax in fig.axes[:2]:
        assert not ax.images
        assert len(ax.collections[0].get_offsets()) == 3
    fig.canvas.draw()
    f = grid.assign(eta1=[0.1, 0.2, 0.3], a11w1=[0.0, 0.1, 0.2], eta2=[-1.0, -0.8, -0.9],
                    a21w1=[0.1, 0.0, -0.1], a22w2=[0.2, 0.1, 0.0])
    res = dm.DenMarkResult(fit=None, gene="geneA", stan_data={}, grid=grid, info={}, settings={})
    fig = dm.plot_fields(res, fields=f)
    for k in (0, 1, 3, 4, 5):                       # the five maps; axes 2 is blank
        ax = fig.axes[k]
        assert not ax.images
        assert len(ax.collections[0].get_offsets()) == 3
    fig.canvas.draw()
