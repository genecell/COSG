"""scanpy is an optional dependency: only `plotMarkerDotplot` needs it.

1.0.4 required scanpy for the whole package, because `_plotting.py` imported
it at module level. That cost a plain ``pip install cosg`` ~60 MB and 9
packages (statsmodels, seaborn, umap-learn, patsy, ...) to support one
plotting function.

The import is now lazy and scanpy is declared as the ``cosg[dotplot]`` extra.
The load-bearing claim is that ``_dendrogram_order`` reproduces
``scanpy.tl.dendrogram``, so ordering behaviour did not change when scanpy
left the import path — that is what ``test_matches_scanpy_dendrogram`` pins.
"""
from __future__ import annotations

import ast
import pathlib

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from cosg._plotting import _dendrogram_order, _require_scanpy

PLOTTING = pathlib.Path(__file__).resolve().parent.parent / "cosg" / "_plotting.py"


@pytest.fixture
def adata():
    """Five well-separated groups in a 6-d representation, fixed seed."""
    rng = np.random.default_rng(0)
    groups, blocks = [], []
    centres = rng.normal(size=(5, 6)) * 5
    for i, centre in enumerate(centres):
        blocks.append(centre + rng.normal(scale=0.3, size=(20, 6)))
        groups += [f"g{i}"] * 20
    X = np.vstack(blocks)
    ad = AnnData(X=X.astype(np.float32))
    ad.obs["group"] = pd.Categorical(groups)
    ad.obsm["X_pca"] = X
    return ad


def test_scanpy_not_imported_at_module_level():
    """The whole point: `import cosg` must not pull scanpy in."""
    offenders = []
    for node in ast.parse(PLOTTING.read_text()).body:  # module level only
        names = []
        if isinstance(node, ast.Import):
            names = [a.name.split(".")[0] for a in node.names]
        elif isinstance(node, ast.ImportFrom) and node.module:
            names = [node.module.split(".")[0]]
        offenders += [f"line {node.lineno}" for n in names if n == "scanpy"]
    assert not offenders, f"scanpy imported at module level: {offenders}"


def test_import_cosg_without_scanpy(monkeypatch):
    """Simulate a bare install: cosg and its plotting names still import."""
    import builtins

    real_import = builtins.__import__

    def no_scanpy(name, *args, **kwargs):
        if name == "scanpy" or name.startswith("scanpy."):
            raise ImportError("No module named 'scanpy'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", no_scanpy)
    import cosg

    for name in ("plotMarkerDendrogram", "plotMarkerDotplot", "plotMarkerStream"):
        assert hasattr(cosg, name)

    with pytest.raises(ImportError, match=r"cosg\[dotplot\]"):
        _require_scanpy("plotMarkerDotplot")


def test_dendrogram_order_writes_scanpy_layout(adata):
    """sc.pl.dotplot(dendrogram=True) reads adata.uns; the keys must match."""
    _dendrogram_order(adata, groupby="group", use_rep="X_pca")
    res = adata.uns["dendrogram_group"]
    for key in (
        "linkage", "groupby", "use_rep", "cor_method", "linkage_method",
        "categories_ordered", "categories_idx_ordered", "dendrogram_info",
        "correlation_matrix",
    ):
        assert key in res, f"missing key scanpy writes: {key}"
    cats = list(adata.obs["group"].cat.categories)
    assert sorted(res["categories_ordered"]) == sorted(cats)
    assert [cats[i] for i in res["categories_idx_ordered"]] == list(
        res["categories_ordered"]
    )


def test_matches_scanpy_dendrogram(adata):
    """The reimplementation must agree with scanpy's, or ordering silently moved."""
    sc = pytest.importorskip("scanpy")

    ours = _dendrogram_order(adata, groupby="group", use_rep="X_pca", inplace=False)
    sc.tl.dendrogram(adata, groupby="group", use_rep="X_pca")
    theirs = adata.uns["dendrogram_group"]

    assert list(ours["categories_ordered"]) == list(theirs["categories_ordered"])
    np.testing.assert_allclose(
        ours["correlation_matrix"], theirs["correlation_matrix"], atol=1e-10
    )
    np.testing.assert_allclose(ours["linkage"], theirs["linkage"], atol=1e-10)


def test_rejects_missing_rep_and_single_category(adata):
    with pytest.raises(KeyError, match="X_umap"):
        _dendrogram_order(adata, groupby="group", use_rep="X_umap")

    adata.obs["one"] = pd.Categorical(["only"] * adata.n_obs)
    with pytest.raises(ValueError, match="at least 2"):
        _dendrogram_order(adata, groupby="one", use_rep="X_pca")
