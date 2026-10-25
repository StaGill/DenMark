"""Short fits on the toy data; skipped when CmdStan is not installed."""

import numpy as np
import pytest

import denmark as dm

from .conftest import cmdstan_available

pytestmark = pytest.mark.skipif(not cmdstan_available(), reason="CmdStan is not installed")

# The fast settings used across the documentation.
FAST = dict(grid_size=10, basis=10, chains=2, iter_warmup=300, iter_sampling=300, seed=2024)
KEY_VARS = ["beta0", "beta1", "ell[1]", "ell[2]", "Astar[2,1]", "sigma[1]", "sigma[2]"]


@pytest.fixture(scope="module")
def fit(tumor):
    return dm.DenMark(tumor, gene="geneA", quiet=True, **FAST)


@pytest.fixture(scope="module")
def fit_tissue(tumor, tissue):
    return dm.DenMark(tumor, gene="geneA", boundary=tissue, quiet=True, **FAST)


def test_fit(fit):
    assert fit.gene == "geneA"
    s = fit.summary()
    assert list(s["variable"]) == KEY_VARS
    assert {"mean", "sd", "q2.5", "q97.5", "rhat", "ess_bulk", "ess_tail"} <= set(s.columns)
    assert np.isfinite(s[["mean", "sd", "q2.5", "q97.5"]].to_numpy(float)).all()
    rho = fit.rho
    assert rho.shape == (FAST["chains"] * FAST["iter_sampling"],)
    assert np.isfinite(rho).all() and (np.abs(rho) <= 1).all()
    assert len(fit.fields()) == fit.stan_data["n"]
    assert list(fit.decay()["ranges"]["field"]) == ["density", "expression"]
    assert repr(fit).startswith("DenMark fit")


def test_save_and_load(fit, tmp_path):
    path = dm.save_result(fit, tmp_path / "fit.pkl")
    np.testing.assert_array_equal(dm.load_result(path).rho, fit.rho)


def test_fit_with_a_tissue_boundary(fit_tissue, tissue, tmp_path):
    """The variable-area model on the 82 grid cells inside the tissue; geneA
    has a true rho of 0.6."""
    fit = fit_tissue
    assert fit.info["vararea"] and fit.stan_data["n"] == 82
    assert np.shape(fit.stan_data["log_grid_area"]) == (82,)
    assert "vararea" in fit.fit.metadata.cmdstan_config["model"]    # DenMark_vararea_model
    assert fit.rho.mean() > 0
    fields = fit.fields()
    assert len(fields) == 82
    np.testing.assert_allclose(fields["log_lambda1"] - fields["eta1"], fit.stan_data["log_grid_area"])
    loaded = dm.load_result(dm.save_result(fit, tmp_path / "fit.pkl"))
    assert loaded.info["boundary"] == tissue
    np.testing.assert_array_equal(loaded.rho, fit.rho)


def test_fields_map_with_a_tissue_boundary(fit_tissue):
    pytest.importorskip("matplotlib").use("Agg")
    import matplotlib.pyplot as plt

    fig = dm.plot_fields(fit_tissue)
    assert all(len(ax.lines) == 1 for ax in fig.axes if ax.images)   # the tissue outline
    plt.close(fig)


def test_scan_two_genes(tumor, tissue):
    scan = dm.scan_genes(tumor, genes=["geneA", "geneC"], boundary=tissue, verbose=False, **FAST)
    assert list(scan["gene"]) == ["geneA", "geneC"]      # sorted by decreasing rho
    assert scan["error"].isna().all()
    assert np.isfinite(scan[["rho_mean", "rho_lower", "rho_upper"]].to_numpy(float)).all()
    assert (scan["rho_lower"] <= scan["rho_upper"]).all()
    assert set(scan["class"]) <= {"PCG", "NCG", "not DCG"}
