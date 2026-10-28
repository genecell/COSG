"""COSG's scores are finite on every path, and a filtered feature scores -1.

- The chunked GPU path, picked automatically when the matrix exceeds 70 % of
  free GPU memory, marked filtered features with -inf; ``indexByGene`` keeps
  it and ``iqrLogNormalize`` then returns NaN, which breaks the dendrogram and
  stream plots.
- With ``cpu_chunk_size=0``, ``batch_key`` and a dense matrix, a group too small
  in every batch had its cosine set to NaN, and the penalty's row sums carried
  the NaN into every group's score; NaN sorts first, so the "markers" were
  arbitrary. Such a group has no evidence and scores 0, as on the streaming and
  GPU batch paths.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import cosg
from test_per_group_floor import numpy_gpu  # noqa: F401  (fixture)


def _adata(dense=False, rare=False, n=300, f=60, seed=0):
    import anndata as ad
    rng = np.random.default_rng(seed)
    gi = rng.integers(0, 3, n)
    rate = np.full((3, f), 0.1)
    for g in range(3):
        rate[g, g * 10:(g + 1) * 10] = 0.7
    X = np.log1p(((rng.random((n, f)) < rate[gi]) * rng.integers(1, 5, (n, f))).astype(np.float32))
    lab = np.array(["A", "B", "C"])[gi]
    if rare:
        lab = lab.copy()
        lab[:4] = "R"                     # 2 cells per batch, under the 3-cell threshold
    obs = pd.DataFrame({"grp": pd.Categorical(lab), "batch": np.where(np.arange(n) % 2 == 0, "b1", "b2")},
                       index=[f"c{i}" for i in range(n)])
    a = ad.AnnData(X if dense else sp.csr_matrix(X), obs=obs)
    a.var_names = [f"g{j}" for j in range(f)]
    return a


def _scores(a):
    s = a.uns["cosg"]["scores"]
    return {g: np.asarray(s[g], dtype=float) for g in s.dtype.names}


def _removed(a):
    c = a.uns["cosg"]
    return {g: set(np.asarray(c["names"][g])[np.asarray(c["scores"][g]) == -1]) for g in c["names"].dtype.names}


@pytest.mark.parametrize("dense", [False, True])
def test_the_chunked_gpu_path_writes_minus_one(numpy_gpu, dense):  # noqa: F811
    runs = {}
    for name, kw in (("cpu", {}), ("gpu", dict(device="gpu")), ("gpu_chunked", dict(device="gpu", gpu_chunk_size=64))):
        a = _adata(dense=dense)
        cosg.cosg(a, groupby="grp", n_genes_user=a.n_vars, **kw)
        runs[name] = a
        assert all(np.isfinite(v).all() for v in _scores(a).values()), name
    assert _removed(runs["gpu_chunked"]) == _removed(runs["gpu"]) == _removed(runs["cpu"])
    assert any(_removed(runs["cpu"]).values())                 # the filter removed something
    df = cosg.indexByGene(runs["gpu_chunked"].uns["cosg"]["COSG"], set_nan_to_zero=True)
    assert np.isfinite(cosg.iqrLogNormalize(df).to_numpy()).all()


@pytest.mark.parametrize("path", [dict(cpu_chunk_size=0, dense=True), dict(cpu_chunk_size=0, dense=False),
                                  dict(dense=True), dict(dense=False)])
def test_a_group_too_small_in_every_batch_has_finite_scores(path):
    path = dict(path)
    a = _adata(dense=path.pop("dense"), rare=True)
    cosg.cosg(a, groupby="grp", batch_key="batch", n_genes_user=a.n_vars, **path)
    s = _scores(a)
    assert all(np.isfinite(v).all() for v in s.values())
    assert set(np.unique(s["R"])) <= {0.0, -1.0}             # no evidence: 0, or filtered
    # the other groups are scored as if R were not there to be scored
    b = _adata(dense=False, rare=True)
    cosg.cosg(b, groupby="grp", batch_key="batch", n_genes_user=b.n_vars)
    for g in ("A", "B", "C"):
        got = dict(zip(a.uns["cosg"]["names"][g], s[g]))
        want = dict(zip(b.uns["cosg"]["names"][g], _scores(b)[g]))
        for gene in want:
            assert got[gene] == pytest.approx(want[gene], abs=1e-6), (g, gene)


def test_a_group_too_small_in_every_batch_on_the_gpu_code(numpy_gpu):  # noqa: F811
    a = _adata(rare=True)
    cosg.cosg(a, groupby="grp", batch_key="batch", n_genes_user=a.n_vars, device="gpu")
    s = _scores(a)
    assert all(np.isfinite(v).all() for v in s.values())
    assert set(np.unique(s["R"])) <= {0.0, -1.0}
