"""Stan data: keys, types, centering and the HSGP domain."""

import warnings

import numpy as np
import pandas as pd
import pytest

import denmark as dm

KEYS = {"n", "y1", "y2", "d", "mstar", "coords", "log_grid_area", "mrange", "indices", "L",
        "is_centerted_PHI"}


@pytest.fixture
def cells():
    rng = np.random.default_rng(7)
    n = 400
    return pd.DataFrame({"x": rng.uniform(2.0, 8.0, n), "y": rng.uniform(-1.0, 3.0, n),
                         "g": rng.poisson(2.0, n)})


def test_stan_data_keys_and_types(cells):
    sd = dm.prepare_data(cells, gene="g", grid_size=(6, 4), basis=(5, 3))["stan_data"]
    assert set(sd) == KEYS
    assert (sd["n"], sd["d"], sd["mstar"], sd["is_centerted_PHI"]) == (24, 2, 15, 0)
    for key in ("n", "d", "mstar", "is_centerted_PHI"):
        assert isinstance(sd[key], int), key
    for key in ("log_grid_area", "mrange"):
        assert isinstance(sd[key], float), key
    for key in ("y1", "y2", "indices"):
        assert np.issubdtype(sd[key].dtype, np.integer), key
    assert sd["y1"].shape == (24,) and sd["y2"].shape == (24,)
    assert sd["coords"].shape == (24, 2) and sd["coords"].dtype == np.float64
    assert isinstance(sd["L"], list) and len(sd["L"]) == 2
    np.testing.assert_array_equal(sd["indices"], dm.hsgp_indices((5, 3)))
    assert sd["y1"].sum() == len(cells) and sd["y2"].sum() == cells["g"].sum()


def test_outcomes_follow_the_grid(cells):
    out = dm.prepare_data(cells, gene="g", grid_size=(6, 4), basis=5)
    np.testing.assert_array_equal(out["stan_data"]["y1"], out["grid"]["Y"])
    np.testing.assert_array_equal(out["stan_data"]["y2"], out["grid"]["M"])
    g = dm.buildgridpp((6, 4), cells, gene="g")
    np.testing.assert_array_equal(out["info"]["Ymat"], g["Ymat"])
    np.testing.assert_array_equal(out["info"]["Mmat"], g["Mmat"])


def test_coordinates_are_centered(cells):
    out = dm.prepare_data(cells, gene="g", grid_size=(6, 4), basis=5)
    sd, info, grid = out["stan_data"], out["info"], out["grid"]
    bb = info["bbox"]
    np.testing.assert_allclose(info["centre"], [(bb[0] + bb[1]) / 2, (bb[2] + bb[3]) / 2])
    np.testing.assert_allclose(sd["coords"], grid[["x", "y"]].to_numpy() - info["centre"])
    np.testing.assert_allclose(grid[["cx", "cy"]].to_numpy(), sd["coords"])
    np.testing.assert_allclose(sd["coords"].mean(axis=0), 0.0, atol=1e-12)
    # centroids sit half a grid step inside [-L1, L1] x [-L2, L2]
    dx, dy = info["step"]
    np.testing.assert_allclose(sd["coords"].max(axis=0), [info["L1"] - dx / 2, info["L2"] - dy / 2])
    np.testing.assert_allclose(sd["coords"].min(axis=0), [dx / 2 - info["L1"], dy / 2 - info["L2"]])


@pytest.mark.parametrize("boundfactor", [1.0, 1.5, 2.2])
def test_hsgp_domain_is_boundfactor_times_half_range(cells, boundfactor):
    out = dm.prepare_data(cells, gene="g", grid_size=5, basis=4, boundfactor=boundfactor)
    bb, info = out["info"]["bbox"], out["info"]
    half = np.array([(bb[1] - bb[0]) / 2, (bb[3] - bb[2]) / 2])
    np.testing.assert_allclose([info["L1"], info["L2"]], half)
    np.testing.assert_allclose(out["stan_data"]["L"], boundfactor * half)
    assert info["boundfactor"] == boundfactor


def test_bbox_sets_the_domain(cells):
    out = dm.prepare_data(cells, gene="g", grid_size=4, basis=4, bbox=(0, 10, -2, 4))
    np.testing.assert_allclose(out["info"]["centre"], [5.0, 1.0])
    np.testing.assert_allclose([out["info"]["L1"], out["info"]["L2"]], [5.0, 3.0])
    np.testing.assert_allclose(out["stan_data"]["L"], [7.5, 4.5])


def test_scale_terms(cells):
    out = dm.prepare_data(cells, gene="g", grid_size=(6, 4), basis=5)
    sd, info = out["stan_data"], out["info"]
    c = sd["coords"]
    dmax = np.hypot(c[:, 0].max() - c[:, 0].min(), c[:, 1].max() - c[:, 1].min())
    assert info["max_dist"] == pytest.approx(dmax)
    assert sd["mrange"] == pytest.approx(dmax / 5.5)
    assert sd["log_grid_area"] == pytest.approx(np.log(info["step"][0] * info["step"][1]))
    assert info["n_cells"] == len(cells)


def test_info_without_a_boundary(cells):
    out = dm.prepare_data(cells, gene="g", grid_size=(6, 4), basis=5)
    info = out["info"]
    assert info["vararea"] is False and info["boundary"] is None
    assert (info["n_grid_total"], info["n_interior"], info["n_boundary"]) == (24, 24, 0)
    assert isinstance(out["stan_data"]["log_grid_area"], float)


def test_make_init(cells):
    sd = dm.prepare_data(cells, gene="g", grid_size=(6, 4), basis=(5, 3))["stan_data"]
    init = dm.make_init(sd, rho0=0.5)
    # beta0 starts at the log number of cells per unit area, beta1 at the log count per cell
    area = 24 * np.exp(sd["log_grid_area"])
    assert init["beta0"] == pytest.approx(np.log(len(cells) / area))
    assert init["beta1"] == pytest.approx(np.log(sd["y2"].sum() / (sd["y1"] + 0.05).sum()))
    assert init["ell"] == [sd["mrange"]] * 2
    assert len(init["betab1"]) == len(init["betab2"]) == 15
    np.testing.assert_allclose(init["Astar"], [[1, 0], [0.5, np.sqrt(0.75)]])
    assert dm.make_init(dict(sd, y1=0 * sd["y1"]))["beta0"] == pytest.approx(np.log(1 / area))
    with pytest.raises(ValueError, match="rho0"):
        dm.make_init(sd, rho0=1.0)


def test_boundfactor_below_one(cells):
    with pytest.raises(ValueError, match="boundfactor"):
        dm.prepare_data(cells, gene="g", boundfactor=0.9)


def test_supplied_half_ranges_are_only_checked(cells):
    half = dm.prepare_data(cells, gene="g")["info"]["L1"]
    with pytest.warns(UserWarning, match="the data value is used"):
        out = dm.prepare_data(cells, gene="g", L1=2 * half)
    assert out["info"]["L1"] == half
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        dm.prepare_data(cells, gene="g", L1=half)
