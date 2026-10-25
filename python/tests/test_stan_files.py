"""Stan files: the copies in denmark/stan/, their lookup and the compile cache."""

import os
from pathlib import Path

import pandas as pd
import pytest

import denmark as dm
from denmark import core

from .conftest import REPO_ROOT

PKG_STAN = Path(core.__file__).resolve().parent / "stan"
R_STAN = REPO_ROOT / "inst" / "stan"
STAN_FILES = ["DenMark-vararea-umi.stan", "DenMark-vararea.stan", "DenMark.stan"]


def test_packaged_stan_files():
    """The DenMark model and its two variable-area variants, nothing else."""
    assert sorted(p.name for p in PKG_STAN.glob("*.stan")) == STAN_FILES


def test_stan_files_match_the_r_package():
    """inst/stan/ is the single source; denmark/stan/ holds byte-identical copies."""
    if not (PKG_STAN.is_dir() and R_STAN.is_dir()):
        pytest.skip("denmark/stan/ or inst/stan/ is not in this checkout")
    names = sorted(p.name for p in PKG_STAN.glob("*.stan"))
    assert names == sorted(p.name for p in R_STAN.glob("*.stan"))
    for name in names:
        assert (PKG_STAN / name).read_bytes() == (R_STAN / name).read_bytes(), name


def test_model_is_found_inside_the_package():
    assert core._stan_path() == (PKG_STAN / "DenMark.stan").resolve()
    assert core._stan_path(vararea=True) == (PKG_STAN / "DenMark-vararea.stan").resolve()


def test_tissue_boundary_needs_a_vararea_stan_file():
    """Checked before compiling: one log_grid_area per grid cell."""
    shapely = pytest.importorskip("shapely")
    cells = pd.DataFrame({"x": [0.0, 1.0, 2.0, 0.5], "y": [0.0, 1.0, 0.0, 0.5], "g": [1, 0, 2, 3]})
    with pytest.raises(ValueError, match=r"needs `vector\[n\] log_grid_area`"):
        dm.DenMark(cells, gene="g", grid_size=2, basis=2, stan_file=PKG_STAN / "DenMark.stan",
                   boundary=shapely.box(0, 0, 2, 1))


def test_stan_dir_overrides_the_package(tmp_path):
    (tmp_path / "DenMark.stan").write_text("// a user copy\n")
    old = dm.STAN_DIR
    try:
        dm.STAN_DIR = tmp_path
        assert core.STAN_DIR == tmp_path
        assert core._stan_path() == (tmp_path / "DenMark.stan").resolve()
    finally:
        dm.STAN_DIR = old
    assert core.STAN_DIR is old


def test_missing_stan_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        core._stan_path(stan_file=tmp_path / "absent.stan")


def test_denmarkmain_rejects_m1():
    cells = pd.DataFrame({"x": [0.0, 1.0, 2.0], "y": [0.0, 1.0, 2.0], "g": [1, 0, 2]})
    with pytest.raises(ValueError, match="was removed in denmark 0.9.0"):
        dm.DenMarkmain(2, cells, 2, 1.5, denmark_version="M1")


@pytest.mark.parametrize("kwargs", [{"model": "M2"}, {"compute_waic": True}])
def test_removed_arguments_fail_before_fitting(kwargs):
    cells = pd.DataFrame({"x": [0.0, 1.0, 2.0], "y": [0.0, 1.0, 2.0], "g": [1, 0, 2]})
    with pytest.raises(TypeError, match="was removed in denmark 0.9.0"):
        dm.DenMark(cells, gene="g", **kwargs)


def test_scan_genes_rejects_fit_simple():
    cells = pd.DataFrame({"x": [0.0, 1.0, 2.0], "y": [0.0, 1.0, 2.0], "g": [1, 0, 2]})
    with pytest.raises(TypeError, match="was removed in denmark 0.9.0"):
        dm.scan_genes(cells, fit_simple=True)


def test_cache_copy_is_rewritten_only_when_changed(tmp_path, monkeypatch):
    monkeypatch.setenv("DENMARK_CACHE_DIR", str(tmp_path / "cache"))
    src = tmp_path / "src" / "model.stan"
    src.parent.mkdir()
    src.write_text("parameters { real mu; }\n")
    copy = core._cached_stan_file(src)
    assert tmp_path / "cache" in copy.parents
    assert copy.name == src.name and copy.read_bytes() == src.read_bytes()
    os.utime(copy, (1, 1))   # an unchanged source leaves the copy, and its time stamp, alone
    assert core._cached_stan_file(src) == copy and copy.stat().st_mtime == 1
    src.write_text("parameters { real nu; }\n")
    assert core._cached_stan_file(src) == copy and copy.read_bytes() == src.read_bytes()


def test_same_named_files_are_cached_apart(tmp_path, monkeypatch):
    monkeypatch.setenv("DENMARK_CACHE_DIR", str(tmp_path / "cache"))
    copies = []
    for folder in ("a", "b"):
        src = tmp_path / folder / "DenMark.stan"
        src.parent.mkdir()
        src.write_text(f"// {folder}\n")
        copies.append(core._cached_stan_file(src))
    assert copies[0] != copies[1]
    assert [c.read_text() for c in copies] == ["// a\n", "// b\n"]
