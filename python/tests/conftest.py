"""Shared fixtures and helpers for the denmark tests."""

from importlib import resources
from pathlib import Path

import pytest

import denmark as dm

PKG_ROOT = Path(__file__).resolve().parents[1]    # python/
REPO_ROOT = Path(__file__).resolve().parents[2]   # the repository, root of the R package


def cmdstan_available():
    """True when cmdstanpy is installed and finds a CmdStan installation."""
    try:
        import cmdstanpy

        cmdstanpy.cmdstan_path()
    except Exception:
        return False
    return True


@pytest.fixture(scope="session")
def toy():
    """The bundled toy data; tests that use it skip until the data file is in place."""
    if not resources.files("denmark").joinpath("data/denmark_toy.csv.gz").is_file():
        pytest.skip("denmark/data/denmark_toy.csv.gz is not in this checkout")
    return dm.datasets.load_toy()


@pytest.fixture(scope="session")
def tumor(toy):
    """The analyzed cell type of the toy data, without the celltype column."""
    return toy[toy["celltype"] == "Tumor"].drop(columns="celltype")


@pytest.fixture(scope="session")
def tissue(toy):
    """The alpha-shape tissue boundary of all toy cells; skips without shapely."""
    pytest.importorskip("shapely")
    return dm.alpha_shape(toy["x"], toy["y"])
