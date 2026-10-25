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
