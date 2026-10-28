"""The GPU path on a store, with numpy standing in for CuPy (its code, not the
device), and the GPU details of ``cosg.cosg``.

- ``device='gpu'`` on a store returned names, scores and groups only:
  ``calculate_pvalues`` was dropped without a word and ``iqr_normalize`` was
  not applied; ``device='auto'`` did the same on any machine with CuPy. The GPU
  path now shares the CPU's post-processing (every output format, the scale
  options), and a request it cannot serve (p-values, feature batching, a
  column-major matrix) raises for ``device='gpu'`` and runs on the CPU for
  ``device='auto'``.
- ``gpu_chunk_size=0``, documented as the monolithic transfer, raised
  ZeroDivisionError.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import cosg
from conftest import requires_cytome
from test_per_group_floor import numpy_gpu  # noqa: F401  (fixture)

N, F = 300, 60


def _data(seed=0):
    rng = np.random.default_rng(seed)
    gi = rng.integers(0, 3, N)
    rate = np.full((3, F), 0.1)
    for g in range(3):
        rate[g, g * 10:(g + 1) * 10] = 0.7
    X = ((rng.random((N, F)) < rate[gi]) * rng.integers(1, 5, (N, F))).astype(np.float32)
    return X, np.array(["A", "B", "C"])[gi]


@pytest.fixture(scope="module")
def stores(tmp_path_factory):
    pytest.importorskip("cytome")
    import cytome
    X, lab = _data()
    out = {}
    d = tmp_path_factory.mktemp("gpu_store")
    for name, column_major in (("row", False), ("col", True)):
        p = str(d / f"{name}.cytome")
        ds = cytome.create(p)
        ds.set_entity("cells", pd.DataFrame({"cell_idx": np.arange(N), "barcode": [f"c{i}" for i in range(N)],
                                             "grp": lab, "batch": np.where(np.arange(N) % 2, "b1", "b2")}))
        ds.set_entity("genes", pd.DataFrame({"gene_idx": np.arange(F), "gene_id": [f"g{j}" for j in range(F)]}))
        ds.flush()
        w = ds.create_layer_writer("RNA_counts", N, F, dtype=np.float32, compression="zlib")
        w.write_chunk(sp.csr_matrix(X), 0)
        w.finalize()
        ds.flush()
        if column_major:
            from test_per_group_floor import _column_major
            _column_major(ds, X, "RNA_counts", width=13)
        ds.close()
        out[name] = p
    return out


def _run(p, **kw):
    base = dict(groupby="grp", layer="counts", verbose=False, n_genes_user=20)
    base.update(kw)
    return cosg.run_cosg_cytome(p, **base)


def _close(a, b):
    assert set(a) == set(b)
    for k in a:
        if k == "names":
            assert np.array_equal(a[k], b[k])
        elif isinstance(a[k], np.ndarray):
            np.testing.assert_allclose(a[k], b[k], rtol=1e-4, atol=1e-6, err_msg=k)
        elif isinstance(a[k], dict) and k == "scores_dict":
            assert set(a[k]) == set(b[k])
            for key in a[k]:
                assert a[k][key] == pytest.approx(b[k][key], rel=1e-4, abs=1e-6)
        elif isinstance(a[k], pd.DataFrame):
            pd.testing.assert_frame_equal(a[k], b[k], rtol=1e-4, atol=1e-6)
        else:
            assert a[k] == b[k], k


@requires_cytome
@pytest.mark.parametrize("kw", [dict(), dict(iqr_normalize=True), dict(output_format="dict"),
                                dict(output_format="long"), dict(output_format="dense"),
                                dict(batch_key="batch", output_format="dict"),
                                dict(scale_rule="positive_median", output_format="dict")])
def test_the_gpu_path_returns_what_the_cpu_returns(numpy_gpu, stores, kw):  # noqa: F811
    _close(_run(stores["row"], device="gpu", **kw), _run(stores["row"], device="cpu", **kw))


@requires_cytome
def test_iqr_normalize_is_applied_on_the_gpu(numpy_gpu, stores):  # noqa: F811
    raw = _run(stores["row"], device="gpu")
    norm = _run(stores["row"], device="gpu", iqr_normalize=True)
    assert not np.allclose(raw["scores"], norm["scores"])


@requires_cytome
@pytest.mark.parametrize("store,kw,why", [
    ("row", dict(calculate_pvalues=True), "calculate_pvalues"),
    ("row", dict(feature_batching="chromosome"), "feature_batching"),
    ("col", dict(), "column-major"),
])
def test_what_the_gpu_cannot_do_raises_on_gpu_and_runs_on_the_cpu_for_auto(numpy_gpu, stores, store, kw, why):  # noqa: F811
    with pytest.raises(NotImplementedError, match=why):
        _run(stores[store], device="gpu", **kw)
    auto = _run(stores[store], device="auto", **kw)
    cpu = _run(stores[store], device="cpu", **kw)
    _close(auto, cpu)
    if kw.get("calculate_pvalues"):
        assert "pvals" in auto


def _adata():
    import anndata as ad
    X, lab = _data()
    a = ad.AnnData(sp.csr_matrix(np.log1p(X)), obs=pd.DataFrame({"grp": pd.Categorical(lab)},
                                                                 index=[f"c{i}" for i in range(N)]))
    a.var_names = [f"g{j}" for j in range(F)]
    return a


def test_gpu_chunk_size_zero_is_the_monolithic_transfer(numpy_gpu):  # noqa: F811
    runs = {}
    for name, cs in (("zero", 0), ("mono", None), ("all", N)):
        a = _adata()
        cosg.cosg(a, groupby="grp", device="gpu", gpu_chunk_size=cs, n_genes_user=10)
        runs[name] = a.uns["cosg"]
    for g in runs["mono"]["names"].dtype.names:
        assert list(runs["zero"]["names"][g]) == list(runs["mono"]["names"][g]) == list(runs["all"]["names"][g])
    assert runs["zero"]["scores"].tobytes() == runs["mono"]["scores"].tobytes() == runs["all"]["scores"].tobytes()
    with pytest.raises(ValueError, match="gpu_chunk_size"):
        cosg.cosg(_adata(), groupby="grp", device="gpu", gpu_chunk_size=-1)


@requires_cytome
def test_min_total_count_is_not_dropped_on_the_gpu(numpy_gpu, stores):  # noqa: F811
    """The GPU path is unbatched and never applied min_total_count, which the
    feature-batched CPU paths apply: device='auto' now runs such a request on
    the CPU, and device='gpu' refuses it."""
    kw = dict(min_total_count=40, feature_batching="auto", output_format="dict")
    with pytest.raises(NotImplementedError, match="min_total_count"):
        _run(stores["row"], device="gpu", **kw)
    _close(_run(stores["row"], device="auto", **kw), _run(stores["row"], device="cpu", **kw))
