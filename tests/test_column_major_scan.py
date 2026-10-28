"""COSG on a column-major matrix scores each feature as it is read.

One matrix is stored row-major and column-major. The row-major store takes
the accumulating path (every chromosome's statistics in memory), the
column-major one the scan. On counts the two must return the same features
with the same scores, to the bit, whatever the block size.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from conftest import requires_cytome, requires_piaso

pytestmark = requires_cytome

N_CELLS, N_PEAKS = 240, 330
CHROMS = ["chr1"] * 120 + ["chr2"] * 130 + ["chrX"] * 80


def _column_major(ds, X, name, width):
    csc = sp.csc_matrix(X)
    dtype = np.dtype(ds._conn.execute(
        "SELECT dtype FROM matrix_meta WHERE matrix_name=?", (name,)).fetchone()[0])
    ds._conn.execute("DELETE FROM matrix_chunks WHERE matrix_name=?", (name,))
    ds._conn.execute("DELETE FROM matrix_csc_chunks WHERE matrix_name=?", (name,))
    k = 0
    for k, start in enumerate(range(0, X.shape[1], width)):
        end = min(start + width, X.shape[1])
        sub = csc[:, start:end].tocsc()
        sub.sort_indices()
        ds._conn.execute(
            "INSERT INTO matrix_csc_chunks (matrix_name, chunk_idx, col_start, col_end, n_nonzero,"
            " data_blob, indices_blob, indptr_blob, dtype, compression) VALUES (?,?,?,?,?,?,?,?,?,'none')",
            (name, k, start, end, int(sub.nnz), sub.data.astype(dtype).tobytes(),
             sub.indices.astype(np.int32).tobytes(), sub.indptr.astype(np.int32).tobytes(), str(dtype)))
    ds._conn.execute("UPDATE matrix_meta SET n_chunks=0, has_csc=1, csc_chunk_size=?, csc_n_chunks=? "
                     "WHERE matrix_name=?", (width, k + 1, name))
    ds._conn.commit()


@pytest.fixture
def stores(tmp_path):
    import cytome
    rng = np.random.default_rng(5)
    groups = np.array(["a", "b", "c", "d", "e"])[rng.integers(0, 5, N_CELLS)]
    # group-specific accessibility, so there are markers to find
    rate = np.full((5, N_PEAKS), 0.05)
    for g in range(5):
        rate[g, g * 60:g * 60 + 50] = 0.6
    gi = np.searchsorted(np.array(["a", "b", "c", "d", "e"]), groups)
    X = ((rng.random((N_CELLS, N_PEAKS)) < rate[gi]) * rng.integers(1, 5, (N_CELLS, N_PEAKS))).astype(np.float32)
    out = []
    for name in ("row", "col"):
        p = tmp_path / f"{name}.cytome"
        ds = cytome.create(str(p))
        ds.set_entity("cells", pd.DataFrame({
            "cell_idx": np.arange(N_CELLS), "barcode": [f"BC{i}" for i in range(N_CELLS)],
            "grp": groups, "batch": np.where(np.arange(N_CELLS) % 3 == 0, "B1", "B2")}))
        starts = np.concatenate([np.arange(CHROMS.count(c)) * 1000 for c in ("chr1", "chr2", "chrX")])
        ds.set_entity("peaks", pd.DataFrame({
            "peak_idx": np.arange(N_PEAKS),
            "peak_id": [f"{c}:{s}-{s + 500}" for c, s in zip(CHROMS, starts)],
            "chr": CHROMS, "start": starts, "end_": starts + 500}))
        ds.add_matrix("ATAC_counts", sp.csr_matrix(X))
        ds.flush()
        if name == "col":
            _column_major(ds, X, "ATAC_counts", width=37)     # chunks straddle chromosomes
        ds.close()
        out.append(str(p))
    return out[0], out[1], X, groups


def _run(path, **kw):
    from cosg import run_cosg_cytome
    base = dict(groupby="grp", n_genes_user=25, mu=100.0, modality="ATAC", layer="counts",
                output_format="dict", iqr_normalize=True, feature_batching="chromosome", verbose=False)
    base.update(kw)
    return run_cosg_cytome(path, **base)


def _scores(res):
    return res["scores_dict"]


@pytest.mark.parametrize("kw", [
    dict(),
    dict(expressed_pct=0.05, min_total_count=2),
    dict(remove_lowly_expressed=False),
    dict(iqr_normalize=False),
    dict(mu=1.0),
    dict(feature_batching="disk"),
    dict(feature_batching="none"),
    dict(n_genes_user=400),                       # more than any group has
], ids=["default", "filters", "no_filter", "raw_scores", "mu1", "disk", "unbatched_request", "all"])
def test_scan_matches_the_accumulating_path_on_counts(stores, kw):
    row, col, _, _ = stores
    ref_kw = dict(kw)
    if ref_kw.get("feature_batching") == "none":
        ref_kw["feature_batching"] = "chromosome"
    a, b = _scores(_run(row, **ref_kw)), _scores(_run(col, **kw))
    assert len(a) > 0
    assert set(a) == set(b)
    assert all(a[k] == b[k] for k in a)


def test_a_cell_mask_is_honoured(stores):
    row, col, _, groups = stores
    mask = np.arange(N_CELLS) % 4 != 0
    a, b = _scores(_run(row, cell_mask=mask)), _scores(_run(col, cell_mask=mask))
    assert a == b and len(a) > 0
    assert a != _scores(_run(col))


def test_batch_key_matches(stores):
    row, col, _, _ = stores
    a, b = _scores(_run(row, batch_key="batch")), _scores(_run(col, batch_key="batch"))
    assert set(a) == set(b) and len(a) > 0
    assert max(abs(a[k] - b[k]) for k in a) < 1e-5


def test_the_block_size_does_not_change_the_result(stores):
    import cytome
    from cosg import _cytome_streaming as C
    _, col, _, _ = stores
    ref = _scores(_run(col))
    real = C._block_features
    try:
        for n in (1, 7, 37, 100, 10_000):
            C._block_features = lambda G, B, block_features=None, n=n: n
            assert _scores(_run(col)) == ref, n
    finally:
        C._block_features = real


def test_ndarray_output_matches(stores):
    row, col, _, _ = stores
    kw = dict(output_format="ndarray", n_genes_user=15)
    a, b = _run(row, **kw), _run(col, **kw)
    assert np.array_equal(a["scores"], b["scores"])
    # equal scores may be listed in either order; the sets per group agree
    for k in range(a["names"].shape[1]):
        assert set(a["names"][:, k]) == set(b["names"][:, k])


def test_ties_are_broken_by_feature_index(tmp_path):
    """Two features with the same score at the cut: the lower index is kept,
    whatever the block size."""
    import cytome
    from cosg import _cytome_streaming as C
    X = np.zeros((40, 12), dtype=np.float32)
    X[:20, :6] = 1.0          # six identical markers of group a
    X[20:, 6:] = 1.0          # six identical markers of group b
    p = tmp_path / "ties.cytome"
    ds = cytome.create(str(p))
    ds.set_entity("cells", pd.DataFrame({"cell_idx": np.arange(40), "barcode": [f"c{i}" for i in range(40)],
                                         "grp": ["a"] * 20 + ["b"] * 20}))
    ds.set_entity("peaks", pd.DataFrame({"peak_idx": np.arange(12), "peak_id": [f"chr1:{i * 1000}-{i * 1000 + 500}" for i in range(12)],
                                         "chr": ["chr1"] * 12, "start": np.arange(12) * 1000, "end_": np.arange(12) * 1000 + 500}))
    ds.add_matrix("ATAC_counts", sp.csr_matrix(X)); ds.flush()
    _column_major(ds, X, "ATAC_counts", width=5); ds.close()
    real = C._block_features
    try:
        got = []
        for n in (1, 4, 100):
            C._block_features = lambda G, B, block_features=None, n=n: n
            got.append(_scores(_run(str(p), n_genes_user=3, iqr_normalize=False)))
    finally:
        C._block_features = real
    assert got[0] == got[1] == got[2]
    assert sorted(f for g, f in got[0] if g == "a") == [f"chr1:{i * 1000}-{i * 1000 + 500}" for i in range(3)]


def test_chunks_that_skip_columns_are_refused(stores):
    import sqlite3
    _, col, _, _ = stores
    c = sqlite3.connect(col)
    c.execute("DELETE FROM matrix_csc_chunks WHERE matrix_name='ATAC_counts' AND chunk_idx=2")
    c.commit(); c.close()
    with pytest.raises(ValueError, match="not contiguous"):
        _run(col)


@requires_piaso
@pytest.mark.parametrize("layer", ["tfidf", "log1p"])
def test_a_normalised_layer_agrees_to_rounding(stores, layer):
    row, col, _, _ = stores
    a, b = _scores(_run(row, layer=layer)), _scores(_run(col, layer=layer))
    assert set(a) == set(b) and len(a) > 0
    assert max(abs(a[k] - b[k]) / max(abs(a[k]), 1e-12) for k in a) < 1e-5


@pytest.mark.parametrize("kw", [dict(), dict(batch_key="batch"), dict(remove_lowly_expressed=False)],
                         ids=["default", "batch_key", "no_filter"])
def test_the_accumulator_budget_splits_passes_without_changing_the_result(stores, kw):
    # "chromosome" on the row-major store: one pass while the sums fit the budget; with a
    # budget of one byte every chromosome is a pass of its own, and the scores are the same
    from cosg import _cytome_streaming as C
    row, _, _, _ = stores
    one = _scores(_run(row, **kw))
    split = _scores(_run(row, max_accumulator_bytes=1, **kw))
    assert len(one) > 0 and set(one) == set(split)
    assert all(one[k] == split[k] for k in one)
    ranges = {"chr1": (0, 120), "chr2": (120, 250), "chrX": (250, 330)}
    assert C._chrom_passes(["chr1", "chr2", "chrX"], ranges, 5, 1, True, 1) == [["chr1"], ["chr2"], ["chrX"]]
    assert C._chrom_passes(["chr1", "chr2", "chrX"], ranges, 5, 1, True, 2 * 1024 ** 3) == [["chr1", "chr2", "chrX"]]
    # chr1 alone is 120 x (5 x 4 + 4) + 120 x 5 x 4 = 5,280 bytes: a budget of 6,000 keeps it alone
    assert C._chrom_passes(["chr1", "chr2", "chrX"], ranges, 5, 1, True, 6_000)[0] == ["chr1"]


@requires_piaso
@pytest.mark.parametrize("layer", ["tfidf", "log1p"])
def test_a_column_chunk_is_normalized_with_its_own_features_factors(stores, layer):
    """The scan normalizes column chunks; TF-IDF indexes its IDF by the chunk's
    column indices, which are local to the chunk, so every chunk took the first
    chunk's IDF (COSG's cosine ignores a positive factor per feature, so its
    scores were right; the values were not). Each column chunk now gets its own
    features' factors, and equals the same columns normalized in full rows."""
    import cytome
    from cosg._cytome_streaming import _resolve_chunk_normalizer
    row, col, X, _ = stores
    rows = np.arange(N_CELLS)
    ds_row, ds_col = cytome.open(row, mode="r"), cytome.open(col, mode="r")
    try:
        full = _resolve_chunk_normalizer(layer, ds_row, "ATAC", True)(sp.csr_matrix(X), rows).toarray()
        by_col = _resolve_chunk_normalizer(layer, ds_col, "ATAC", True)
        for start in (0, 37, 111, 296):
            end = min(start + 37, N_PEAKS)
            part = by_col(sp.csr_matrix(X[:, start:end]), rows, col_start=start).toarray()
            np.testing.assert_allclose(part, full[:, start:end], rtol=1e-6, err_msg=f"{layer} {start}")
    finally:
        ds_row.close()
        ds_col.close()
