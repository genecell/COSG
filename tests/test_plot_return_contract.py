"""A plotting call that draws should not also print an object.

`plotMarkerDotplot` returned its axes on the default path, so every notebook
cell showed `<Axes: ...>` under the figure -- and implied the object was still
useful after `plt.show()`, which for an inline backend it is not. The scanpy
backend already returned None; the native renderer did not.

The contract, which the `show` docstring has always described: `show=True`
draws and returns None, `show=False` returns the axes so the caller can keep
drawing on them.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import cosg

ad = pytest.importorskip("anndata")


@pytest.fixture
def adata():
    rng = np.random.default_rng(0)
    n, g = 60, 12
    a = ad.AnnData(
        sp.csr_matrix(rng.poisson(2, (n, g)).astype(np.float32)),
        obs=pd.DataFrame({"leiden": rng.choice(list("abc"), n)},
                         index=[f"c{i}" for i in range(n)]),
        var=pd.DataFrame(index=[f"g{i}" for i in range(g)]))
    a.obsm["X_umap"] = rng.normal(size=(n, 2))
    cosg.cosg(a, groupby="leiden", key_added="cosg", n_genes_user=3)
    return a


def test_the_default_draws_and_returns_nothing(adata):
    out = cosg.plotMarkerDotplot(adata, groupby="leiden", top_n_genes=2,
                                 use_rep="X_umap", key_cosg="cosg")
    assert out is None, f"the default call returned {type(out).__name__}"


def test_show_false_hands_back_the_axes(adata):
    out = cosg.plotMarkerDotplot(adata, groupby="leiden", top_n_genes=2,
                                 use_rep="X_umap", key_cosg="cosg", show=False)
    assert out is not None and hasattr(out, "figure"), type(out)
    import matplotlib.pyplot as plt
    plt.close("all")
