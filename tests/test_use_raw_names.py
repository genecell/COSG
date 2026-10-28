"""``use_raw=True`` names the genes of the matrix it scores.

scanpy's usual layout keeps every gene in ``adata.raw`` and the variable genes in
``adata``. COSG scored ``adata.raw.X`` but took the names from ``adata.var_names``,
so a marker at column ``j`` of raw was reported as the gene at position ``j`` of
the narrower table (or raised IndexError past its end). It entered in COSG 1.0.4,
which made ``use_raw`` score ``adata.raw.X`` at all: in 1.0.1 to 1.0.3 the line
after the raw branch set the matrix back to ``adata.X``, so ``use_raw`` had no
effect and the names were right.

Every gene-level output must describe raw's genes: the names, and the fold
changes and p-values read alongside them.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import cosg
from test_per_group_floor import numpy_gpu  # noqa: F401  (fixture)

GROUPS = ("A", "B", "C")


def _full(n=300, f=60, seed=0, dense=False):
    """Every group has ten markers among the first thirty genes; the last
    thirty are background, which is what the narrow table keeps."""
    import anndata as ad
    rng = np.random.default_rng(seed)
    gi = rng.integers(0, len(GROUPS), n)
    rate = np.full((len(GROUPS), f), 0.1)
    for g in range(len(GROUPS)):
        rate[g, g * 10:(g + 1) * 10] = 0.7
    X = np.log1p(((rng.random((n, f)) < rate[gi]) * rng.integers(1, 5, (n, f))).astype(np.float32))
    obs = pd.DataFrame({"grp": pd.Categorical(np.array(GROUPS)[gi]),
                        "batch": np.where(np.arange(n) % 2 == 0, "b1", "b2")},
                       index=[f"c{i}" for i in range(n)])
    a = ad.AnnData(X if dense else sp.csr_matrix(X), obs=obs)
    a.var_names = [f"g{j}" for j in range(f)]
    return a


def _narrow_with_raw(full):
    """The scanpy layout: raw holds every gene, adata the last thirty."""
    a = full.copy()
    a.raw = a
    return a[:, [f"g{j}" for j in range(30, full.n_vars)]].copy()


PATHS = {
    "streaming": {},
    "streaming_chunks": dict(cpu_chunk_size=64),
    "legacy": dict(cpu_chunk_size=0),
    "batch_streaming": dict(batch_key="batch"),
    "batch_legacy": dict(batch_key="batch", cpu_chunk_size=0),
}


def _fields(uns):
    return {k: {g: np.asarray(uns[k][g]) for g in uns[k].dtype.names}
            for k in uns if k not in ("params", "COSG")}


def _same(r_raw, r_ref):
    a, b = _fields(r_raw), _fields(r_ref)
    assert set(a) == set(b)
    for k in a:
        for g in a[k]:
            if a[k][g].dtype.kind in "fc":
                np.testing.assert_array_equal(a[k][g], b[k][g], err_msg=f"{k} {g}")
            else:
                assert list(a[k][g]) == list(b[k][g]), (k, g)


@pytest.mark.parametrize("path", list(PATHS))
@pytest.mark.parametrize("dense", [False, True])
def test_use_raw_reports_raw_genes(path, dense):
    full = _full(dense=dense)
    narrow = _narrow_with_raw(full)
    kw = dict(groupby="grp", n_genes_user=8, **PATHS[path])
    cosg.cosg(narrow, use_raw=True, **kw)
    cosg.cosg(full, **kw)
    _same(narrow.uns["cosg"], full.uns["cosg"])
    # the markers are raw-only genes, which the narrow table does not hold
    assert set(narrow.uns["cosg"]["names"]["A"][:5]) <= {f"g{j}" for j in range(10)}


def test_use_raw_fold_changes_and_pvalues_describe_raw_genes():
    full = _full()
    narrow = _narrow_with_raw(full)
    kw = dict(groupby="grp", n_genes_user=8, calculate_logfoldchanges=True, calculate_pvalues=True)
    cosg.cosg(narrow, use_raw=True, **kw)
    cosg.cosg(full, **kw)
    _same(narrow.uns["cosg"], full.uns["cosg"])
    assert "logfoldchanges" in narrow.uns["cosg"] and "pvals" in narrow.uns["cosg"]


def test_use_raw_without_raw_scores_x():
    full = _full()
    a = full.copy()
    cosg.cosg(a, groupby="grp", n_genes_user=8, use_raw=True)
    cosg.cosg(full, groupby="grp", n_genes_user=8)
    _same(a.uns["cosg"], full.uns["cosg"])


@pytest.mark.parametrize("path", [dict(), dict(gpu_chunk_size=64), dict(batch_key="batch")])
def test_use_raw_on_the_gpu_code(numpy_gpu, path):  # noqa: F811
    full = _full()
    narrow = _narrow_with_raw(full)
    kw = dict(groupby="grp", n_genes_user=8, device="gpu", **path)
    cosg.cosg(narrow, use_raw=True, **kw)
    cosg.cosg(full, **kw)
    _same(narrow.uns["cosg"], full.uns["cosg"])
