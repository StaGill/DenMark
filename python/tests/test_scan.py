"""Gene scans: argument checks and the parallel scan, which must return the
same table as the sequential one."""

import time

import numpy as np
import pandas as pd
import pytest

import denmark as dm

GENES = ["geneA", "geneB", "geneC"]


class _FakeFit:
    """Stands in for a DenMark fit, so that the scan logic runs without CmdStan."""

    def __init__(self, dataset, gene):
        v = dataset[gene].to_numpy(float)
        self.rho = np.tanh(np.linspace(-0.5, 0.5, 200) + v.mean() - 2)
        self.stan_data = {"y2": np.array([v.sum()])}
        self.diagnostics = {"max_rhat": 1.0, "min_ess_bulk": 400.0, "num_divergent": 0}
        self.time = 0.01


@pytest.fixture
def fake_fits(monkeypatch):
    def fake(dataset, gene=None, **kwargs):
        time.sleep(0.05)   # lets the threads overlap
        return _FakeFit(dataset, gene)
    monkeypatch.setattr(dm.core, "DenMark", fake)
    monkeypatch.setattr(dm.core, "compile_model", lambda *args, **kwargs: None)


def test_n_jobs_must_be_a_positive_whole_number(tumor):
    for bad in (0, -1, 1.5, True):
        with pytest.raises(ValueError, match="positive whole number"):
            dm.scan_genes(tumor, genes=["geneA"], n_jobs=bad)


def test_parallel_scan_returns_the_sequential_table(fake_fits, tumor):
    seq = dm.scan_genes(tumor, genes=GENES, verbose=False)
    par = dm.scan_genes(tumor, genes=GENES, verbose=False, n_jobs=3)
    pd.testing.assert_frame_equal(par, seq)
    assert set(par["gene"]) == set(GENES)


def test_parallel_scan_resumes_from_save_dir(fake_fits, tumor, tmp_path):
    first = dm.scan_genes(tumor, genes=GENES[:2], verbose=False, n_jobs=2, save_dir=tmp_path)
    again = dm.scan_genes(tumor, genes=GENES, verbose=False, n_jobs=2, save_dir=tmp_path)
    assert sorted(p.name for p in tmp_path.iterdir()) == [f"denmark_{g}.csv" for g in GENES]
    pd.testing.assert_series_equal(
        again.set_index("gene").loc[GENES[:2], "rho_mean"],
        first.set_index("gene").loc[GENES[:2], "rho_mean"])


def test_scan_keeps_the_offset_and_unit_columns(monkeypatch):
    seen, compiled = {}, []

    def fake(dataset, gene=None, **kwargs):
        seen[gene] = list(dataset.columns)
        return _FakeFit(dataset, gene)

    monkeypatch.setattr(dm.core, "DenMark", fake)
    monkeypatch.setattr(dm.core, "compile_model", lambda *args, **kwargs: compiled.append(kwargs))
    cells = pd.DataFrame({"x": [0.0, 1.0, 2.0], "y": [0.0, 1.0, 0.0], "UMI": [10.0, 12.0, 9.0],
                          "geneA": [1, 2, 0], "geneB": [0, 1, 1]})
    res = dm.scan_genes(cells, offset="umi", umi="UMI", verbose=False, n_jobs=2)
    assert sorted(res["gene"]) == ["geneA", "geneB"]
    assert seen["geneA"] == ["x", "y", "UMI", "geneA"]
    assert compiled[-1]["umi"] is True
    units = pd.DataFrame({"unit": [1, 2, 3], "x": [0.0, 1.0, 2.0], "y": [0.0, 1.0, 0.0],
                          "Y": [3, 5, 2], "UMI": [10.0, 12.0, 9.0], "n_spots": [7, 7, 4],
                          "area": [1.0, 1.0, 0.6], "geneA": [1, 2, 0]})
    res = dm.scan_genes(units, units=True, offset="umi", umi="UMI", verbose=False)
    assert list(res["gene"]) == ["geneA"]
    assert seen["geneA"] == ["x", "y", "UMI", "Y", "area", "geneA"]
    with pytest.raises(ValueError, match="Missing column"):
        dm.scan_genes(cells, units=True, verbose=False)
