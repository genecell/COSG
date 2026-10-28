"""Three behaviours fixed together.

A floor of 0 has one meaning on every path. With ``remove_lowly_expressed=True``
a group's threshold is ``max(N_k * floor_k, expressed_min_num_cells_in_target_group)``
everywhere, so ``expressed_pct=0`` leaves the 3-cell minimum. Before, it switched
the filter off altogether on the streaming paths (the default chunked CPU path,
the GPU's chunked path, every cytome path) and kept the minimum on the others, so
a feature seen in one cell of a group could be one of its markers or not
depending on the path. ``remove_lowly_expressed=False`` is the switch.

``min_total_count`` on the feature-batched cytome path applies whenever it is
above 0, as its docstring says. Before, it needed the nonzero counts that only
the floor kept, so a floor of 0 (or ``remove_lowly_expressed=False``) dropped it.

A group label holding "/" (cortex layers: "L2/3 IT") survives ``adata.write``.
h5ad and zarr keys cannot hold "/" (anndata 0.13 refuses them), and the by-group
frame ``adata.uns[key]['COSG']`` named its columns "names::L2/3 IT". The frame
now writes "/" as "%2F"; the record arrays keep the labels exactly, and
``indexByGene`` reads "%2F" back as "/".
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import cosg
from conftest import requires_cytome
from test_per_group_floor import (AD_PATHS, CY_PATHS, GPU_PATHS, GROUPS, LABELS, X,  # noqa: F401
                                  _ad_equal, _ad_run, _adata, _cy_equal, _cy_run, _real_cupy,
                                  _removed, numpy_gpu, requires_gpu, stores)

EVERY = dict(n_genes_user=X.shape[1])            # every feature reported, so every -1 is seen
#: 1e-9 of any group here is far below 3 cells: the floor that is the minimum alone
TINY = 1e-9


def _nnz_by_group(labels=LABELS):
    """Per group, the number of its cells each feature is nonzero in."""
    return {g: np.asarray((X[labels == g] != 0).sum(axis=0)).ravel() for g in GROUPS}


def _under_minimum(min_cells=3):
    nnz = _nnz_by_group()
    return {g: {f"g{j}" for j in np.flatnonzero(nnz[g] < min_cells)} for g in GROUPS}


def _filtered(uns):
    """Per group, the features the expression filter removed: -1, or -inf on the
    GPU's chunked path, which has always marked them so."""
    from test_per_group_floor import _ad_columns
    c = _ad_columns(uns)
    return {g: set(c["names"][g][(c["scores"][g] == -1) | np.isneginf(c["scores"][g])]) for g in c["names"]}


def _positive(res):
    """The (group, feature) pairs a cytome result scores above 0."""
    if "names" in res:
        out = set()
        for k, g in enumerate(map(str, res["groups_order"])):
            s = np.asarray(res["scores"][:, k], dtype=np.float64)
            out |= {(g, str(n)) for n, v in zip(res["names"][:, k], s) if v > 0}
        return out
    if "scores_dict" in res:
        return {(str(g), str(f)) for (g, f), s in res["scores_dict"].items() if s > 0}
    df = res["scores_df"]
    col = "score" if "score" in df else df.columns[-1]
    return {(str(g), str(f)) for g, f, s in zip(df["group"], df["feature"], df[col]) if s > 0}


def _same(r1, r2):
    try:
        _cy_equal(r1, r2)
        return True
    except AssertionError:
        return False


# ── a floor of 0: the minimum of cells, on every path ─────────────────────────

@pytest.mark.parametrize("path", list(AD_PATHS))
def test_anndata_a_zero_floor_keeps_the_minimum_of_cells(path):
    zero = _ad_run(path, 0.0, **EVERY)
    _ad_equal(zero, _ad_run(path, TINY, **EVERY))
    _ad_equal(zero, _ad_run(path, [0.0] * len(GROUPS), **EVERY))
    removed = _removed(zero)
    assert removed == _under_minimum(), "a floor of 0 removes exactly the features under 3 cells"
    assert any(removed.values()), "nothing under 3 cells here; the test proves nothing"
    # the minimum is the filter's: without the filter nothing is removed (cosg.cosg
    # takes no minimum below 1 cell, so the switch is remove_lowly_expressed)
    off = _ad_run(path, 0.0, remove_lowly_expressed=False, **EVERY)
    assert not any(_removed(off).values())


@pytest.mark.parametrize("pair", [("csr_stream", "csr_legacy"), ("dense_stream", "dense_legacy"),
                                  ("csr_stream_chunks", "csc_legacy"), ("csr_batch", "csr_batch_legacy")])
def test_the_paths_agree_on_what_a_zero_floor_removes(pair):
    """Before, the streaming path of each pair removed nothing at a floor of 0 and
    the legacy path removed the features under 3 cells."""
    a, b = (_removed(_ad_run(p, 0.0, **EVERY)) for p in pair)
    assert a == b


def test_score_from_sums_a_zero_floor_keeps_the_minimum_of_cells():
    gidx = np.searchsorted(GROUPS, LABELS)
    lam = sp.csr_matrix((np.ones(len(LABELS)), (np.arange(len(LABELS)), gidx)), shape=(len(LABELS), len(GROUPS)))
    Xs = sp.csr_matrix(X)
    dot = np.asarray((Xs.T @ lam).todense())
    sq = np.asarray(Xs.multiply(Xs).sum(axis=0)).ravel()
    nnz = np.asarray(((Xs != 0).astype(np.float64).T @ lam).todense())
    sizes = np.bincount(gidx, minlength=len(GROUPS)).astype(np.float64)
    zero = cosg.score_from_sums(dot, sq, sizes, nnz, expressed_pct=0.0)
    np.testing.assert_array_equal(zero, cosg.score_from_sums(dot, sq, sizes, nnz, expressed_pct=TINY))
    assert ((zero == -1) == (nnz < 3)).all()
    off = cosg.score_from_sums(dot, sq, sizes, nnz, remove_lowly_expressed=False)
    np.testing.assert_array_equal(
        off, cosg.score_from_sums(dot, sq, sizes, nnz, expressed_pct=0.0,
                                  expressed_min_num_cells_in_target_group=0))


@requires_cytome
@pytest.mark.parametrize("path", list(CY_PATHS))
def test_cytome_a_zero_floor_keeps_the_minimum_of_cells(stores, path):
    zero = _cy_run(stores, path, 0.0, **EVERY)
    _cy_equal(zero, _cy_run(stores, path, TINY, **EVERY))
    off = _cy_run(stores, path, 0.0, remove_lowly_expressed=False, **EVERY)
    _cy_equal(off, _cy_run(stores, path, 0.0, expressed_min_num_cells_in_target_group=0, **EVERY))
    assert not _same(zero, off), "the 3-cell minimum changes nothing here; the test proves nothing"
    # what the minimum removes is never a marker
    under = {(g, f) for g, fs in _under_minimum().items() for f in fs}
    assert not (_positive(zero) & under)
    assert _positive(off) & under, "no feature under 3 cells scores; the test proves nothing"


@requires_cytome
def test_cytome_and_anndata_remove_the_same_features_at_a_zero_floor(stores):
    res = _cy_run(stores, "rna_ndarray", 0.0, **EVERY)
    cy = {str(g): {str(n) for n, s in zip(res["names"][:, k], res["scores"][:, k]) if s == -1}
          for k, g in enumerate(res["groups_order"])}
    assert cy == _removed(_ad_run("csr_stream", 0.0, **EVERY)) == _under_minimum()


@pytest.mark.parametrize("path", list(GPU_PATHS))
def test_gpu_code_a_zero_floor_keeps_the_minimum_of_cells_on_numpy(numpy_gpu, path):
    zero = _ad_run(GPU_PATHS[path], 0.0, **EVERY)
    _ad_equal(zero, _ad_run(GPU_PATHS[path], TINY, **EVERY))
    assert _filtered(zero) == _under_minimum()
    off = _ad_run(GPU_PATHS[path], 0.0, remove_lowly_expressed=False, **EVERY)
    assert not any(_filtered(off).values())


@requires_cytome
@pytest.mark.parametrize("batch_key", [None, "batch"])
def test_cytome_gpu_code_a_zero_floor_keeps_the_minimum_of_cells_on_numpy(numpy_gpu, stores, batch_key):
    run = lambda pct, **kw: cosg.run_cosg_cytome(stores["rna"], groupby="grp", layer="counts", modality="RNA",
                                                 verbose=False, device="gpu", batch_key=batch_key,
                                                 expressed_pct=pct, **EVERY, **kw)
    zero = run(0.0)
    _cy_equal(zero, run(TINY))
    under = {(g, f) for g, fs in _under_minimum().items() for f in fs}
    assert not (_positive(zero) & under)


@requires_gpu
@pytest.mark.parametrize("path", list(GPU_PATHS))
def test_gpu_a_zero_floor_keeps_the_minimum_of_cells_on_the_device(path):
    _real_cupy()
    zero = _ad_run(GPU_PATHS[path], 0.0, **EVERY)
    _ad_equal(zero, _ad_run(GPU_PATHS[path], TINY, **EVERY))
    assert _filtered(zero) == _under_minimum()


@requires_gpu
@requires_cytome
@pytest.mark.parametrize("batch_key", [None, "batch"])
def test_cytome_gpu_a_zero_floor_keeps_the_minimum_of_cells_on_the_device(stores, batch_key):
    _real_cupy()
    run = lambda pct: cosg.run_cosg_cytome(stores["rna"], groupby="grp", layer="counts", modality="RNA",
                                           verbose=False, device="gpu", batch_key=batch_key,
                                           expressed_pct=pct, **EVERY)
    zero = run(0.0)
    _cy_equal(zero, run(TINY))
    under = {(g, f) for g, fs in _under_minimum().items() for f in fs}
    assert not (_positive(zero) & under)


# ── min_total_count on the feature-batched path, whatever the floor ───────────

FEATURE_BATCHED = [p for p, kw in CY_PATHS.items() if kw.get("feature_batching") in ("chromosome", "disk")
                   and not kw.get("min_total_count")]


@requires_cytome
@pytest.mark.parametrize("path", FEATURE_BATCHED)
@pytest.mark.parametrize("filters", [dict(remove_lowly_expressed=False), dict(expressed_pct=0.0),
                                     dict(expressed_pct=0.0, expressed_min_num_cells_in_target_group=0),
                                     dict(expressed_pct=0.1)])
def test_min_total_count_applies_whatever_the_floor(stores, path, filters):
    m = 12
    total = np.asarray((X != 0).sum(axis=0)).ravel()
    low = {f"g{j}" for j in np.flatnonzero(total <= m)}
    kw = dict(filters)
    pct = kw.pop("expressed_pct", 0.1)
    with_m = _cy_run(stores, path, pct, min_total_count=m, **EVERY, **kw)
    assert not {f for _, f in _positive(with_m)} & low, "a feature at or under min_total_count scored"
    without = _cy_run(stores, path, pct, min_total_count=0, **EVERY, **kw)
    if not (kw.get("remove_lowly_expressed", True) and pct >= 0.1):
        assert {f for _, f in _positive(without)} & low, "min_total_count removes nothing; the test proves nothing"


# ── a label holding "/" survives adata.write ─────────────────────────────────

SLASHED = {"a": "L2/3 IT", "b": "L5/6 NP", "c": "Micro/PVM", "d": "d", "e": "L4/5/6 IT"}


@pytest.mark.parametrize("path", ["csr_stream", "dense_stream", "csr_legacy", "csr_batch", "csr_batch_legacy",
                                  "csr_pvalues"])
def test_labels_with_a_slash_survive_adata_write(path, tmp_path):
    import anndata as ad
    kw = dict(AD_PATHS[path])
    a = _adata(kw.pop("fmt"), np.array([SLASHED[g] for g in LABELS]))
    cosg.cosg(a, groupby="grp", n_genes_user=15, mu=1, **kw)
    frame = a.uns["cosg"]["COSG"]
    labels = list(a.uns["cosg"]["names"].dtype.names)
    assert sorted(labels) == sorted(SLASHED.values())                 # the record arrays: exact labels
    assert not any("/" in c for c in frame.columns)
    assert "names::L2%2F3 IT" in frame.columns and "scores::L4%2F5%2F6 IT" in frame.columns
    by_gene = cosg.indexByGene(frame, set_nan_to_zero=True)
    assert list(by_gene.columns) == labels                            # read back as the labels
    a.write_h5ad(tmp_path / "x.h5ad")                                 # anndata 0.13 refused '/' here
    b = ad.read_h5ad(tmp_path / "x.h5ad")
    assert list(b.uns["cosg"]["names"].dtype.names) == labels
    pd.testing.assert_frame_equal(cosg.indexByGene(b.uns["cosg"]["COSG"], set_nan_to_zero=True), by_gene,
                                  check_dtype=False)
    # a frame from an earlier version, with '/' in its column names, reads the same
    old = frame.copy()
    old.columns = [c.replace("%2F", "/") for c in frame.columns]
    pd.testing.assert_frame_equal(cosg.indexByGene(old, set_nan_to_zero=True), by_gene)


def test_labels_with_a_slash_survive_write_zarr(tmp_path):
    pytest.importorskip("zarr")
    import anndata as ad
    a = _adata("csr", np.array([SLASHED[g] for g in LABELS]))
    cosg.cosg(a, groupby="grp", n_genes_user=10, mu=1)
    a.write_zarr(tmp_path / "x.zarr")                                 # zarr never allows '/' in a key
    b = ad.read_zarr(tmp_path / "x.zarr")
    assert list(cosg.indexByGene(b.uns["cosg"]["COSG"]).columns) == list(a.uns["cosg"]["names"].dtype.names)


def test_the_plotting_reader_sees_the_labels():
    from cosg._plotting import _resolve_cosg_scores
    a = _adata("csr", np.array([SLASHED[g] for g in LABELS]))
    cosg.cosg(a, groupby="grp", n_genes_user=10, mu=1)
    scores = _resolve_cosg_scores(a, None, "cosg", "grp", None)
    assert set(scores.columns) == set(a.obs["grp"].cat.categories)
