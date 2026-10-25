"""Import-level checks."""

import re
import subprocess
import sys

import pytest

import denmark as dm

from .conftest import REPO_ROOT

# Public names carried over from DenMark/denmark.py before packaging; these stay public.
API = [
    "filter_genes", "buildgridpp", "grid_pts_centroid", "hsgp_indices", "prepare_data",
    "compile_model", "make_init", "DenMark", "denmark", "DenMarkmain", "DenMarkResult",
    "scan_genes", "matern32", "simulate", "save_result",
    "load_result", "plot_grid", "plot_fields", "plot_decay", "plot_scan",
]
# The tissue boundary, added with Mingchi Xu's revision of the R functions.
BOUNDARY_API = ["alpha_shape", "TissueBoundary", "plot_boundary"]


def test_version_is_string():
    assert isinstance(dm.__version__, str)
    assert dm.__version__.count(".") >= 1


def test_version_matches_the_r_package():
    desc = REPO_ROOT / "DESCRIPTION"
    if not desc.exists():
        pytest.skip("DESCRIPTION is not in this checkout")
    m = re.search(r"^Version:\s*(\S+)", desc.read_text(encoding="utf-8"), flags=re.M)
    assert m is not None, "DESCRIPTION has no Version field"
    assert m.group(1) == dm.__version__


def test_public_api_present():
    assert set(API + BOUNDARY_API) <= set(dm.__all__)
    for name in dm.__all__:
        assert hasattr(dm, name), f"denmark is missing public name: {name}"
    assert callable(dm.datasets.load_toy)


def test_shapely_is_optional():
    """shapely is imported only by the functions that need it."""
    code = "import sys, denmark; sys.exit('shapely' in sys.modules)"
    assert subprocess.run([sys.executable, "-c", code]).returncode == 0
