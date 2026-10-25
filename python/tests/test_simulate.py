"""simulate(): cells in a rectangle or an irregular tissue section."""

import numpy as np
import pandas as pd
import pytest

import denmark as dm
from denmark import core

RHO = {"up": 0.5, "down": -0.5}


def test_simulate_returns_cells_counts_and_the_truth():
    sim = dm.simulate(rho=RHO, width=3, height=2, n_cells=300, seed=3)
    assert list(sim.columns) == ["x", "y", "up", "down"]
    assert 0 < len(sim) <= 300
    assert sim["x"].between(0, 3).all() and sim["y"].between(0, 2).all()
    for gene in RHO:
        assert np.issubdtype(sim[gene].dtype, np.integer) and (sim[gene] >= 0).all()
    truth = sim.attrs["truth"]
    assert truth["rho"] == RHO and truth["window"] == (3, 2) and truth["sigma2"] == 0.8
    assert truth["tissue"] == "rectangle" and truth["outline"] is None
    pd.testing.assert_frame_equal(dm.simulate(rho=RHO, width=3, height=2, n_cells=300, seed=3),
                                  sim)


def test_irregular_tissue_with_other_cells():
    sim = dm.simulate(n_cells=1000, seed=2, tissue="irregular", n_other=800, celltype="Tumor")
    assert list(sim.columns) == ["x", "y", "celltype", "geneA", "geneB", "geneC"]
    assert sim["celltype"].value_counts().to_dict() == {"Tumor": 1000, "Other": 800}
    # every cell lies in the lobed tissue section
    assert core._in_tissue(sim["x"].to_numpy(), sim["y"].to_numpy(), 6, 4).all()
    other = sim[sim["celltype"] == "Other"]
    assert other[["geneA", "geneB", "geneC"]].to_numpy().mean() == pytest.approx(0.3, abs=0.05)
    truth = sim.attrs["truth"]
    assert truth["tissue"] == "irregular"
    outline = pd.DataFrame(truth["outline"])
    assert list(outline.columns) == ["x", "y"] and len(outline) == 361
    assert outline.iloc[0].tolist() == pytest.approx(outline.iloc[-1].tolist())   # closed
    assert outline["x"].between(0, 6).all() and outline["y"].between(0, 4).all()
    # the truth is plain data, so frames built from the simulation still concatenate
    parts = [sim[sim["celltype"] == "Tumor"], sim[sim["celltype"] == "Other"]]
    assert len(pd.concat(parts)) == len(sim)


def test_unknown_tissue():
    with pytest.raises(ValueError, match="tissue"):
        dm.simulate(tissue="round")
