"""cytome is an optional dependency: the plain-AnnData path must not need it.

1.1.1 regression: ``_cytome_dataset_classes()`` guarded its first lookup
(``from cytome.core.dataset import CytomeDataset``) with try/except
ImportError, but not the second (``getattr(_cytome(), "Dataset", None)``),
and ``_cytome()`` raises ImportError when cytome is absent. Every
``cosg.cosg(adata)`` call dispatches through ``_is_cytome_input`` ->
``_is_cytome_dataset`` -> ``_cytome_dataset_classes``, so a default install
(no ``cosg[cytome]`` extra) crashed on any AnnData input — exactly the
install the docstring promises a no-op (empty tuple) fallback for.

cytome IS installed in the dev environment, so these tests block the import
themselves: drop the function-attribute cache, then poison ``sys.modules``
with None for ``cytome`` and every already-imported ``cytome.*`` submodule
(a None entry makes ``import`` raise ImportError; the submodule entries
matter because ``from cytome.core.dataset import ...`` would otherwise be
served straight from ``sys.modules`` without touching the parent).

Also here: the matplotlib >= 3.9 regression guard. ``matplotlib.cm.get_cmap``
was removed in 3.9; ``_plotting.py`` must only use ``plt.get_cmap``.
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import cosg
from cosg import _cytome_streaming


@pytest.fixture
def blocked_cytome(monkeypatch):
    """Make ``import cytome`` (and any submodule) raise ImportError.

    Clears the ``_cytome_dataset_classes`` cache before AND after: before,
    so the lookup actually runs under the blocked import; after, so the
    empty tuple cached during this test cannot leak into tests that use
    the real cytome.
    """
    _cytome_streaming._cytome_dataset_classes.__dict__.pop("_cache", None)
    monkeypatch.setitem(sys.modules, "cytome", None)
    for name in [k for k in list(sys.modules) if k.startswith("cytome.")]:
        monkeypatch.setitem(sys.modules, name, None)
    yield
    _cytome_streaming._cytome_dataset_classes.__dict__.pop("_cache", None)


@pytest.fixture
def adata():
    """Small synthetic AnnData: 100 cells x 50 genes, 3 separated clusters."""
    rng = np.random.default_rng(0)
    n_cells, n_genes = 100, 50
    clusters = np.array([f"c{i % 3}" for i in range(n_cells)])
    X = rng.poisson(1.0, size=(n_cells, n_genes)).astype(np.float32)
    for i, cluster in enumerate(("c0", "c1", "c2")):
        rows = clusters == cluster
        X[np.ix_(rows, np.arange(i * 10, i * 10 + 10))] += 5.0
    ad = AnnData(X=X)
    ad.obs["cluster"] = pd.Categorical(clusters)
    return ad


def test_cosg_on_anndata_without_cytome(blocked_cytome, adata):
    """The headline bug: plain AnnData marker detection on a bare install."""
    cosg.cosg(adata, groupby="cluster", n_genes_user=10)
    assert "cosg" in adata.uns
    names = adata.uns["cosg"]["names"]
    assert set(names.dtype.names) == {"c0", "c1", "c2"}


def test_dataset_classes_empty_tuple_without_cytome(blocked_cytome):
    """Contract: the docstring promises a no-op (empty tuple) fallback."""
    assert _cytome_streaming._cytome_dataset_classes() == ()


def test_import_cosg_without_cm_get_cmap(monkeypatch):
    """matplotlib removed ``cm.get_cmap`` in 3.9; cosg must not touch it."""
    monkeypatch.delattr("matplotlib.cm.get_cmap", raising=False)
    for name in [k for k in list(sys.modules) if k == "cosg" or k.startswith("cosg.")]:
        monkeypatch.delitem(sys.modules, name)
    import cosg  # noqa: F811 — fresh import with cm.get_cmap gone

    # No call site may reach cm.get_cmap: grep the installed source.
    src_dir = pathlib.Path(cosg.__file__).resolve().parent
    offenders = []
    for py in sorted(src_dir.glob("*.py")):
        for lineno, line in enumerate(py.read_text().splitlines(), start=1):
            if "cm.get_cmap(" in line and not line.lstrip().startswith("#"):
                offenders.append(f"{py.name}:{lineno}: {line.strip()}")
    assert not offenders, (
        "cm.get_cmap was removed in matplotlib 3.9; use plt.get_cmap:\n"
        + "\n".join(offenders)
    )
