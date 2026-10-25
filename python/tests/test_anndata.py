"""AnnData objects: the cell table read from one, and the scan written back."""

import numpy as np
import pandas as pd
import pytest

import denmark as dm

anndata = pytest.importorskip("anndata")
sparse = pytest.importorskip("scipy.sparse")

GENES = ["geneA", "geneB", "geneC"]


@pytest.fixture
def toy_adata(toy):
    obs = pd.DataFrame({"celltype": toy["celltype"].to_numpy()},
                       index=[f"cell{k}" for k in range(len(toy))])
    adata = anndata.AnnData(X=sparse.csr_matrix(toy[GENES].to_numpy(float)), obs=obs,
                            var=pd.DataFrame(index=GENES))
    adata.obsm["spatial"] = toy[["x", "y"]].to_numpy()
    return adata


def test_from_anndata_returns_the_cell_table(toy, toy_adata):
    cells = dm.from_anndata(toy_adata, celltype_key="celltype")
    assert list(cells.columns) == ["x", "y", "celltype"] + GENES
    assert list(cells.index[:2]) == ["cell0", "cell1"]
    pd.testing.assert_frame_equal(cells.reset_index(drop=True), toy[list(cells.columns)],
                                  check_dtype=False)
    two = dm.from_anndata(toy_adata, "celltype", genes=["geneC", "geneA"])
    assert list(two.columns) == ["x", "y", "celltype", "geneC", "geneA"]


def test_from_anndata_checks_its_input(toy_adata):
    with pytest.raises(ValueError, match="adata.obs"):
        dm.from_anndata(toy_adata, "nope")
    with pytest.raises(ValueError, match="adata.obsm"):
        dm.from_anndata(toy_adata, "celltype", spatial_key="X_umap")
    with pytest.raises(ValueError, match="Not variables"):
        dm.from_anndata(toy_adata, "celltype", genes=["geneZ"])
    toy_adata.layers["lognorm"] = toy_adata.X.log1p()
    with pytest.warns(UserWarning, match="raw counts"):
        dm.from_anndata(toy_adata, "celltype", layer="lognorm")


def test_annotate_writes_the_scan_into_var(toy_adata):
    scan = pd.DataFrame({"gene": ["geneC", "geneA"], "rho_mean": [-0.6, 0.6],
                         "rho_lower": [-0.8, 0.3], "rho_upper": [-0.3, 0.8],
                         "class": ["NCG", "PCG"]})
    out = dm.annotate(toy_adata, scan)
    assert out is toy_adata
    v = toy_adata.var
    assert v.loc["geneA", "denmark_rho"] == 0.6 and v.loc["geneC", "denmark_rho"] == -0.6
    assert np.isnan(v.loc["geneB", "denmark_rho"])
    assert v.loc["geneC", "denmark_class"] == "NCG"
    assert {"denmark_rho_lower", "denmark_rho_upper"} <= set(v.columns)
