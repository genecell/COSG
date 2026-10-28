"""A default call on a store writes nothing into it.

A default ``run_cosg_cytome`` (or ``cosg.cosg``) on a path opened the store
read-write and left caches in it: ``RNA_cell_depth`` (through cytome's
``modality_cell_depth``) for log1p, ``<modality>_tfidf_params`` (through
PIASO) for TF-IDF. The depth cache was trusted whenever its length matched, so
a counts matrix replaced by one of the same shape kept the old depths. Now:

- log1p takes each cell's depth from the data: the row sums of the chunk being
  normalized (one pass over the columns first on a column-major matrix,
  ``cells.n_fragments`` for ATAC), the numbers cytome's ``matrix_sums`` gave;
- a path is opened read-only; only ``write_to_cytome`` writes, and a store it
  cannot write fails before the computation, not after.
"""
from __future__ import annotations

import os
import warnings

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import cosg
from conftest import requires_cytome, requires_piaso

pytestmark = requires_cytome

N, F = 300, 60


def _data(seed=0):
    rng = np.random.default_rng(seed)
    gi = rng.integers(0, 3, N)
    rate = np.full((3, F), 0.1)
    for g in range(3):
        rate[g, g * 10:(g + 1) * 10] = 0.7
    X = ((rng.random((N, F)) < rate[gi]) * rng.integers(1, 9, (N, F))).astype(np.float32)
    return X, np.array(["A", "B", "C"])[gi]


def _store(path, modality="RNA", column_major=False, X=None, n_fragments=None):
    import cytome
    if X is None:
        X, lab = _data()
    else:
        lab = _data()[1]
    ds = cytome.create(str(path))
    cells = pd.DataFrame({"cell_idx": np.arange(N), "barcode": [f"c{i}" for i in range(N)], "grp": lab})
    if n_fragments is not None:
        cells["n_fragments"] = n_fragments
    ds.set_entity("cells", cells)
    if modality == "RNA":
        ds.set_entity("genes", pd.DataFrame({"gene_idx": np.arange(F), "gene_id": [f"g{j}" for j in range(F)]}))
    else:
        starts = np.arange(F) * 1000
        ds.set_entity("peaks", pd.DataFrame({"peak_idx": np.arange(F), "chr": ["chr1"] * F, "start": starts,
                                             "end_": starts + 500, "peak_id": [f"chr1:{s}-{s + 500}" for s in starts]}))
    ds.flush()
    w = ds.create_layer_writer(f"{modality}_counts", N, F, dtype=np.float32, compression="zlib")
    w.write_chunk(sp.csr_matrix(X), 0)
    w.finalize()
    ds.flush()
    if column_major:
        from test_per_group_floor import _column_major
        _column_major(ds, X, f"{modality}_counts", width=13)
    ds.close()
    return str(path)


def _state(p):
    """The store's metadata keys and its files' bytes."""
    import cytome
    with open(p, "rb") as fh:
        main = fh.read()
    wal = os.path.getsize(p + "-wal") if os.path.exists(p + "-wal") else 0
    ds = cytome.open(p, mode="r")
    try:
        keys = sorted(ds.metadata.keys())
    finally:
        ds.close()
    return keys, main, wal


@pytest.mark.parametrize("column_major", [False, True])
def test_a_default_log1p_call_writes_nothing(tmp_path, column_major):
    p = _store(tmp_path / "s.cytome", column_major=column_major)
    before = _state(p)
    r1 = cosg.run_cosg_cytome(p, groupby="grp", verbose=False)              # layer auto -> log1p
    r2 = cosg.cosg(p, groupby="grp", output_format="dict")
    after = _state(p)
    assert after[0] == before[0], set(after[0]) ^ set(before[0])
    assert after[1] == before[1] and after[2] == 0
    assert "RNA_cell_depth" not in after[0]
    assert r1["names"].shape[0] == 50 and r2["scores_dict"]


@requires_piaso
def test_a_default_tfidf_call_writes_nothing(tmp_path):
    p = _store(tmp_path / "a.cytome", modality="ATAC", n_fragments=np.full(N, 1000.0))
    before = _state(p)
    cosg.run_cosg_cytome(p, groupby="grp", modality="ATAC", verbose=False)    # layer auto -> tfidf
    after = _state(p)
    assert after[0] == before[0] and "ATAC_tfidf_params" not in after[0]
    assert after[1] == before[1]


def _expected_depth(X):
    return np.asarray(sp.csr_matrix(X).sum(axis=1), dtype=np.float64).ravel()


@pytest.mark.parametrize("column_major", [False, True])
def test_the_log1p_depth_is_cytomes(tmp_path, column_major):
    """The depth equals cytome's matrix_sums, from which RNA_cell_depth was
    computed, so results on a store without a stale cache do not move."""
    import cytome
    from cytome.utils.modality import matrix_sums
    from cosg._cytome_streaming import _log1p_depth, _resolve_chunk_normalizer
    X, _ = _data()
    p = _store(tmp_path / "d.cytome", column_major=column_major)
    ds = cytome.open(p, mode="r")
    try:
        ref = matrix_sums(ds, "RNA", "counts")[0]
        depth = _log1p_depth(ds, "RNA")
        if column_major:
            np.testing.assert_array_equal(depth, ref)
        else:
            assert depth is None                       # taken from each chunk
        norm = _resolve_chunk_normalizer("log1p", ds, "RNA", use_cached_stats=True)
        from cosg._normalize import normalize_chunk_log1p
        chunk = sp.csr_matrix(X[:50])
        got = norm(chunk, np.arange(50))
        want = normalize_chunk_log1p(chunk, ref[:50], 1e4)
        assert (got != want).nnz == 0
    finally:
        ds.close()


def test_atac_log1p_depth_is_n_fragments(tmp_path):
    import cytome
    from cosg._cytome_streaming import _log1p_depth
    frags = np.arange(N, dtype=np.float64) + 500
    p = _store(tmp_path / "f.cytome", modality="ATAC", n_fragments=frags)
    ds = cytome.open(p, mode="r")
    try:
        np.testing.assert_array_equal(_log1p_depth(ds, "ATAC"), frags)
    finally:
        ds.close()


def test_a_stale_depth_cache_is_not_used(tmp_path):
    """A matrix replaced by one of the same shape kept the cached depths of the
    old one: the result followed the cache, not the data."""
    import cytome
    p = _store(tmp_path / "fresh.cytome")
    q = _store(tmp_path / "stale.cytome")
    ds = cytome.open(q)
    ds.metadata["RNA_cell_depth"] = np.full(N, 7.0)          # right length, wrong values
    ds.flush()
    ds.close()
    kw = dict(groupby="grp", layer="log1p", verbose=False, output_format="dict", n_genes_user=F)
    assert cosg.run_cosg_cytome(q, **kw)["scores_dict"] == cosg.run_cosg_cytome(p, **kw)["scores_dict"]


def test_a_path_is_opened_read_only(tmp_path):
    """On a store that cannot be written, a read-write open falls back to
    reading with a ReadOnlyWarning; a read-only open has nothing to warn of."""
    import cytome
    p = _store(tmp_path / "ro.cytome")
    for f in (p, p + "-wal", p + "-shm"):
        if os.path.exists(f):
            os.chmod(f, 0o444)
    os.chmod(tmp_path, 0o555)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", getattr(cytome, "ReadOnlyWarning", UserWarning))
            r = cosg.run_cosg_cytome(p, groupby="grp", layer="counts", verbose=False)
        assert r["names"].shape == (50, 3)
    finally:
        os.chmod(tmp_path, 0o755)


def test_write_to_cytome_fails_before_the_computation(tmp_path, monkeypatch):
    import cytome
    import cosg._cytome_streaming as cs
    p = _store(tmp_path / "w.cytome")
    calls = []
    real = cs._get_labels_and_gene_names
    monkeypatch.setattr(cs, "_get_labels_and_gene_names",
                        lambda *a, **k: (calls.append(1), real(*a, **k))[1])
    ds = cytome.open(p, mode="r")
    try:
        with pytest.raises(Exception, match="read.only|write"):
            cosg.run_cosg_cytome(ds, groupby="grp", layer="counts", verbose=False, write_to_cytome=True)
    finally:
        ds.close()
    assert not calls                                         # nothing was computed
    # and on a writable path it writes the markers, nothing else
    before = _state(p)[0]
    cosg.run_cosg_cytome(p, groupby="grp", layer="counts", verbose=False, write_to_cytome=True,
                         key_added="cosg_test")
    assert sorted(set(_state(p)[0]) - set(before)) == ["cosg_test"]
