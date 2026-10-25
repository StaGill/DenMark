"""Figures with a tissue boundary (matplotlib, Agg backend)."""

import numpy as np
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
