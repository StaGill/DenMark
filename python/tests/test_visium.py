"""Total-count offset, seven-spot Visium pooling and pre-aggregated units,
mirroring tests/testthat/test-visium.R and test-prior-scale.R."""

import numpy as np
import pandas as pd
import pytest
from scipy.spatial.distance import pdist

import denmark as dm
from denmark import core


def hex_spots(nx=14, ny=12, seed=1):
    """A hexagonal Visium-like lattice, pitch 100: rows alternate by half a pitch."""
    i, j = np.meshgrid(np.arange(nx), np.arange(ny), indexing="xy")
    i, j = i.ravel(), j.ravel()
    rng = np.random.default_rng(seed)
    return pd.DataFrame({"x": 100 * (i + 0.5 * (j % 2)), "y": 100 * np.sqrt(3) / 2 * j,
                         "Y": rng.poisson(8, i.size), "UMI": rng.poisson(5000, i.size),
                         "geneA": rng.poisson(3, i.size), "geneB": rng.poisson(1, i.size)})


def umi_grid_df():
    """3 x 2 grid on [0, 3] x [0, 2]: grid cell k holds k cells, each with 2
    counts of the gene and 10 total counts."""
    k = np.repeat(np.arange(1, 7), np.arange(1, 7))
    return pd.DataFrame({"x": (k - 1) % 3 + 0.5, "y": (k - 1) // 3 + 0.5, "g": 2.0, "u": 10.0})


def test_buildgridpp_sums_total_counts():
    g = dm.buildgridpp((3, 2), umi_grid_df(), gene="g", bbox=(0, 3, 0, 2), umi="u")
    np.testing.assert_allclose(g["grid"]["T"], 10 * np.arange(1, 7))
    np.testing.assert_array_equal(g["grid"]["M"], 2 * np.arange(1, 7))
    g0 = dm.buildgridpp((3, 2), umi_grid_df(), gene="g", bbox=(0, 3, 0, 2))
    assert "T" not in g0["grid"]


def test_umi_column_is_not_the_default_gene():
    df = umi_grid_df()[["x", "y", "u", "g"]]
    assert dm.buildgridpp((3, 2), df, bbox=(0, 3, 0, 2), umi="u")["gene"] == "g"


def test_offset_umi_stan_data():
    p = dm.prepare_data(umi_grid_df(), gene="g", grid_size=(3, 2), basis=3,
                        bbox=(0, 3, 0, 2), offset="umi", umi="u")
    sd = p["stan_data"]
    np.testing.assert_allclose(sd["log_umi"], np.log(10 * np.arange(1, 7) + 0.05))
    assert np.size(sd["log_grid_area"]) == sd["n"]
    np.testing.assert_allclose(sd["log_grid_area"], np.zeros(6))
    assert p["info"]["offset"] == "umi" and p["info"]["units"] is False
    p0 = dm.prepare_data(umi_grid_df(), gene="g", grid_size=(3, 2), basis=3, bbox=(0, 3, 0, 2))
    assert "log_umi" not in p0["stan_data"]
    assert np.ndim(p0["stan_data"]["log_grid_area"]) == 0
    assert p0["info"]["offset"] == "cells"


def test_offset_umi_needs_the_column():
    df = umi_grid_df()
    with pytest.raises(ValueError, match="umi="):
        dm.prepare_data(df, gene="g", grid_size=3, offset="umi")
    with pytest.raises(ValueError, match="not a column"):
        dm.prepare_data(df, gene="g", grid_size=3, offset="umi", umi="nope")
    with pytest.raises(ValueError, match="offset"):
        dm.prepare_data(df, gene="g", grid_size=3, offset="reads")
    df.loc[0, "u"] = -1
    with pytest.raises(ValueError, match="non-negative"):
        dm.prepare_data(df, gene="g", grid_size=3, offset="umi", umi="u")


def test_make_init_beta1_per_total_count():
    p = dm.prepare_data(umi_grid_df(), gene="g", grid_size=(3, 2), basis=3,
                        bbox=(0, 3, 0, 2), offset="umi", umi="u")
    init = dm.make_init(p["stan_data"])
    expected = np.log(p["stan_data"]["y2"].sum() / (10 * np.arange(1, 7) + 0.05).sum())
    assert init["beta1"] == pytest.approx(expected)


def test_umi_offset_selects_the_umi_model():
    assert core._stan_path(umi=True).name == "DenMark-vararea-umi.stan"
    assert core._stan_path(vararea=True).name == "DenMark-vararea.stan"
    assert core._stan_path().name == "DenMark.stan"


def test_denmark_refuses_a_stan_file_without_log_umi(tmp_path):
    f = tmp_path / "custom.stan"
    f.write_text("data {\n  int n;\n  vector[n] log_grid_area;\n}\n")
    with pytest.raises(ValueError, match="log_umi"):
        dm.DenMark(umi_grid_df(), gene="g", grid_size=(3, 2), basis=3, bbox=(0, 3, 0, 2),
                   offset="umi", umi="u", stan_file=f)


def test_pool_visium_keeps_every_count():
    spots = hex_spots()
    u = dm.pool_visium(spots)
    sid = u.attrs["spot_unit"]
    assert len(sid) == len(spots)
    assert set(np.unique(sid)) == set(u["unit"])
    for c in ("Y", "UMI", "geneA", "geneB"):
        assert u[c].sum() == pytest.approx(spots[c].sum())
    assert u["n_spots"].sum() == len(spots)
    np.testing.assert_allclose(u["area"], u["n_spots"] * np.pi * 27.5 ** 2)
    assert (u["n_spots"] >= 3).all()


def test_full_pools_are_a_spot_and_its_six_neighbours():
    spots = hex_spots()
    u = dm.pool_visium(spots, min_spots=1)
    sid = u.attrs["spot_unit"]
    assert (u["n_spots"] <= 7).all()
    full = u.loc[u["n_spots"] == 7, "unit"]
    assert len(full) >= 5
    for k in full:
        s = spots.loc[sid == k, ["x", "y"]].to_numpy()
        d = np.sqrt(((s - s.mean(axis=0)) ** 2).sum(axis=1))
        np.testing.assert_allclose(np.sort(d), [0] + [100] * 6, atol=1e-8)


def test_small_pools_are_folded():
    spots = hex_spots()
    u3 = dm.pool_visium(spots, min_spots=3)
    u1 = dm.pool_visium(spots, min_spots=1)
    assert (u1["n_spots"] < 3).any()
    assert len(u3) < len(u1)
    assert u3["n_spots"].sum() == u1["n_spots"].sum()


def test_pool_visium_checks_its_input():
    spots = hex_spots()
    with pytest.raises(ValueError, match="Missing column"):
        dm.pool_visium(spots.drop(columns="Y"))
    bad = spots.astype({"geneA": float})          # pandas 3 refuses 0.5 in an integer column
    bad.loc[0, "geneA"] = 0.5
    with pytest.raises(ValueError, match="integer counts"):
        dm.pool_visium(bad)
    with pytest.raises(ValueError, match="not in"):
        dm.pool_visium(spots, genes=["nope"])
    u = dm.pool_visium(spots[["x", "y", "Y", "geneA"]], umi=None)
    assert list(u.columns) == ["unit", "x", "y", "Y", "n_spots", "area", "geneA"]


def test_units_stan_data():
    u = dm.pool_visium(hex_spots())
    p = dm.prepare_data(u, gene="geneA", basis=4, boundfactor=1.5,
                        units=True, offset="umi", umi="UMI")
    sd = p["stan_data"]
    assert sd["n"] == len(u)
    np.testing.assert_array_equal(sd["y1"], u["Y"])
    np.testing.assert_array_equal(sd["y2"], u["geneA"])
    np.testing.assert_allclose(sd["log_grid_area"], np.log(u["area"]))
    np.testing.assert_allclose(sd["log_umi"], np.log(u["UMI"] + 0.05))
    cx = u["x"] - (u["x"].min() + u["x"].max()) / 2
    cy = u["y"] - (u["y"].min() + u["y"].max()) / 2
    np.testing.assert_allclose(sd["coords"], np.column_stack([cx, cy]))
    np.testing.assert_allclose(sd["L"], 1.5 * np.array([cx.abs().max(), cy.abs().max()]))
    assert sd["mrange"] == pytest.approx(pdist(np.column_stack([cx, cy])).max() / 5.5)
    assert sd["mstar"] == 16
    assert p["info"]["units"] and p["info"]["vararea"]
    assert p["grid"]["i"].isna().all()
    assert "log_umi" not in dm.prepare_data(u, gene="geneA", basis=4, units=True)["stan_data"]
    with pytest.raises(ValueError, match="gene="):
        dm.prepare_data(u, basis=4, units=True)
    with pytest.raises(ValueError, match="Missing column"):
        dm.prepare_data(u, gene="geneA", basis=4, units=True, area="nope")


def test_max_dist_is_the_largest_pairwise_distance():
    rng = np.random.default_rng(3)
    pts = np.column_stack([rng.uniform(-2, 5, 300), rng.normal(0, 3, 300)])
    assert core._max_dist(pts) == pytest.approx(pdist(pts).max())
    line = np.column_stack([np.arange(10.0), 2 * np.arange(10.0)])     # collinear
    assert core._max_dist(line) == pytest.approx(np.hypot(9, 18))
    assert core._max_dist(np.array([[1.0, 2.0]])) == 0.0
    th = np.linspace(0, 2 * np.pi, 361)[:-1]
    assert core._max_dist(np.column_stack([np.cos(th), np.sin(th)])) == pytest.approx(2)


def test_fields_use_the_total_count_offset():
    class _Fit:
        def stan_variable(self, name):
            return {"w1": np.array([[0.1, 0.2, 0.3], [-0.1, 0.0, 0.1]]),
                    "w2": np.array([[0.0, 0.1, 0.2], [0.2, 0.1, 0.0]])}[name]

    grid = pd.DataFrame({"grid_index": [0, 1, 2], "i": np.nan, "j": np.nan,
                         "x": [0.0, 1.0, 2.0], "y": [0.0, 1.0, 0.0], "Y": [3, 5, 2], "M": [1, 2, 0]})
    sd = {"n": 3, "log_grid_area": np.zeros(3),
          "log_umi": np.log(np.array([100.0, 200.0, 300.0]) + 0.05)}
    res = dm.DenMarkResult(fit=_Fit(), gene="g", stan_data=sd, grid=grid, info={}, settings={})
    res.draws = lambda: pd.DataFrame({"beta0": [0.0, 0.1], "beta1": [-1.0, -0.9],
                                      "a11": [1.0, 1.1], "a21": [0.3, 0.2], "a22": [0.5, 0.6]})
    f = res.fields()
    np.testing.assert_allclose(f["log_lambda2"] - f["eta2"], sd["log_umi"])
    del sd["log_umi"]
    f = res.fields()
    np.testing.assert_allclose(f["log_lambda2"] - f["eta2"], np.log(grid["Y"] + 0.05))
