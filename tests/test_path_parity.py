"""Every path reports the same markers: one parity test.

COSG has many paths to the same scores: on an AnnData, the chunked CPU path
(CSR, CSC, dense), the legacy path (``cpu_chunk_size=0``) and the GPU paths
(monolithic, chunked, batch); on a store, the unbatched path, chromosome and
disk batching on a row-major matrix, the column scan of a column-major one, and
the GPU path. Most of the audit's findings were places where they had drifted
apart: the sentinel of a filtered feature, empty groups, the choice among tied
scores, the tail of the names array, what counts as expressed, ``use_raw``,
``n_genes_user`` above the number of features.

On integer counts every sum is exact, so the CPU paths must agree to the bit
on names and to float32 rounding on scores; the GPU code (numpy standing in
for CuPy) works in float32 and agrees to 1e-5. The data hold:

- ten features with the same values (exact ties among a group's markers);
- features under the detection floor (the -1 sentinel) and features no cell
  of a group expresses (score 0), with explicit stored zeros among them;
- a group of six cells, two to three per batch;
- a batch that holds cells of one group only (the "one-group stratum");
- more features asked for (``n_genes_user``) than there are.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import cosg
from conftest import requires_cytome
from test_per_group_floor import numpy_gpu, _column_major  # noqa: F401  (fixture, helper)

GROUPS = np.array(["A", "B", "C", "D"])
N_TOP = 50                      # more than the 40 features
CHROMS = ["chr1"] * 15 + ["chr2"] * 15 + ["chrX"] * 10


def _data(seed=7, n=240, f=40):
    rng = np.random.default_rng(seed)
    gi = np.r_[np.repeat([0, 1, 2], (n - 6) // 3), np.full(6, 3)]
    gi = np.r_[gi, np.zeros(n - gi.size, dtype=int)]
    rate = np.full((4, f), 0.08)
    rate[1, 10:20] = 0.6
    rate[2, 20:30] = 0.5
    rate[3, 30:34] = 0.9
    X = ((rng.random((n, f)) < rate[gi]) * rng.integers(1, 4, (n, f))).astype(np.float32)
    # ten exact ties among A's markers: the same column, repeated
    tie = ((rng.random(n) < np.where(gi == 0, 0.7, 0.05)) * 2).astype(np.float32)
    X[:, :10] = tie[:, None]
    X[:, 37] = 0.0                                   # a feature no cell expresses
    batch = np.where(gi == 2, "b3", np.where(np.arange(n) % 2 == 0, "b1", "b2"))   # b3: group C only
    Xs = sp.csr_matrix(X)
    # explicit stored zeros in feature 35 for half of A's cells
    rows = np.flatnonzero(gi == 0)[::2]
    lil = Xs.tolil()
    for i in rows:
        lil[i, 35] = 1.0
    Xs = lil.tocsr()
    for i in rows:
        s, e = Xs.indptr[i], Xs.indptr[i + 1]
        Xs.data[s:e][Xs.indices[s:e] == 35] = 0.0
    return Xs, GROUPS[gi], batch


X, LABELS, BATCH = _data()
NAMES = [f"p{j}" for j in range(X.shape[1])]


# ── AnnData paths ──────────────────────────────────────────────────────────

def _adata(fmt="csr", raw_wider=False):
    import anndata as ad
    mat = {"csr": X.copy(), "csc": X.tocsc(), "dense": X.toarray()}[fmt]
    obs = pd.DataFrame({"grp": pd.Categorical(LABELS, categories=list(GROUPS) + ["unused"]),
                        "batch": BATCH}, index=[f"c{i}" for i in range(len(LABELS))])
    a = ad.AnnData(mat, obs=obs)
    a.var_names = NAMES
    if raw_wider:
        # raw holds every feature, adata the last ten: use_raw must report raw's
        a.raw = a
        a = a[:, NAMES[30:]].copy()
    return a


AD_PATHS = {
    "csr_stream": dict(fmt="csr"),
    "csr_stream_chunks": dict(fmt="csr", cpu_chunk_size=37),
    "csc_stream": dict(fmt="csc"),
    "dense_stream": dict(fmt="dense"),
    "csr_legacy": dict(fmt="csr", cpu_chunk_size=0),
    "dense_legacy": dict(fmt="dense", cpu_chunk_size=0),
    "csr_use_raw": dict(fmt="csr", raw_wider=True, use_raw=True),
    "csr_legacy_use_raw": dict(fmt="csr", raw_wider=True, use_raw=True, cpu_chunk_size=0),
}
AD_GPU_PATHS = {
    "csr_gpu": dict(fmt="csr", device="gpu"),
    "dense_gpu": dict(fmt="dense", device="gpu"),
    "csr_gpu_chunked": dict(fmt="csr", device="gpu", gpu_chunk_size=37),
    "dense_gpu_chunked": dict(fmt="dense", device="gpu", gpu_chunk_size=37),
}


def _ad_result(path, batch=False, **extra):
    kw = dict(path)
    a = _adata(kw.pop("fmt"), kw.pop("raw_wider", False))
    if batch:
        kw["batch_key"] = "batch"
    cosg.cosg(a, groupby="grp", n_genes_user=N_TOP, **kw, **extra)
    u = a.uns["cosg"]
    return {g: (list(u["names"][g]), np.asarray(u["scores"][g], dtype=np.float64))
            for g in u["names"].dtype.names}


# ── store paths ────────────────────────────────────────────────────────────

def _store(path, modality, column_major=False):
    import cytome
    n, f = X.shape
    ds = cytome.create(str(path))
    ds.set_entity("cells", pd.DataFrame({"cell_idx": np.arange(n), "barcode": [f"c{i}" for i in range(n)],
                                         "grp": LABELS, "batch": BATCH}))
    if modality == "RNA":
        ds.set_entity("genes", pd.DataFrame({"gene_idx": np.arange(f), "gene_id": NAMES}))
    else:
        starts = np.arange(f) * 1000
        ds.set_entity("peaks", pd.DataFrame({"peak_idx": np.arange(f), "chr": CHROMS, "start": starts,
                                             "end_": starts + 500, "peak_id": NAMES}))
    ds.flush()
    w = ds.create_layer_writer(f"{modality}_counts", n, f, dtype=np.float32, compression="zlib")
    w.write_chunk(X, 0)
    w.finalize()
    ds.flush()
    if column_major:
        _column_major(ds, X, f"{modality}_counts", width=11)     # chunks straddle chromosomes
    ds.close()
    return str(path)


@pytest.fixture(scope="module")
def stores(tmp_path_factory):
    pytest.importorskip("cytome")
    d = tmp_path_factory.mktemp("parity")
    return {"rna": _store(d / "rna.cytome", "RNA"), "atac": _store(d / "atac.cytome", "ATAC"),
            "atac_col": _store(d / "atac_col.cytome", "ATAC", column_major=True)}


CY_PATHS = {
    "rna_none": dict(store="rna", modality="RNA", feature_batching="none"),
    "atac_none": dict(store="atac", modality="ATAC", feature_batching="none"),
    "atac_chromosome": dict(store="atac", modality="ATAC", feature_batching="chromosome"),
    "atac_one_pass_per_chromosome": dict(store="atac", modality="ATAC", feature_batching="chromosome",
                                         max_accumulator_bytes=1),
    "atac_disk": dict(store="atac", modality="ATAC", feature_batching="disk"),
    "atac_column_major": dict(store="atac_col", modality="ATAC", feature_batching="chromosome"),
}


def _cy_result(stores, path, batch=False, **extra):
    kw = dict(CY_PATHS[path]) if isinstance(path, str) else dict(path)
    p = stores[kw.pop("store")]
    if batch:
        kw["batch_key"] = "batch"
    r = cosg.run_cosg_cytome(p, groupby="grp", n_genes_user=N_TOP, layer="counts", verbose=False,
                             **kw, **extra)
    return {str(g): (list(r["names"][:, k]), r["scores"][:, k].astype(np.float64))
            for k, g in enumerate(r["groups_order"])}


# ── the comparison ─────────────────────────────────────────────────────────

def _same(got, ref, rtol=1e-6, label=""):
    assert list(got) == list(ref), (label, list(got), list(ref))
    for g in ref:
        assert got[g][0] == ref[g][0], (label, g, got[g][0][:12], ref[g][0][:12])
        np.testing.assert_allclose(got[g][1], ref[g][1], rtol=rtol, atol=1e-7, err_msg=f"{label} {g}")


@pytest.fixture(scope="module", params=[(False, True), (True, True), (False, False), (True, False)],
                ids=["no_batch", "batch", "no_batch_no_filter", "batch_no_filter"])
def reference(request):
    batch, filt = request.param
    return (batch, filt), _ad_result(AD_PATHS["csr_stream"], batch=batch, remove_lowly_expressed=filt)


def test_the_reference_holds_what_the_test_is_about(reference):
    (batch, filt), ref = reference
    assert list(ref) == ["A", "B", "C", "D"]                  # the unused category is not scored
    names, scores = ref["A"]
    assert len(names) == X.shape[1]                            # n_genes_user above the features: all
    assert names[:10] == [f"p{j}" for j in range(10)]          # the ten tied markers, by index
    assert len(set(np.round(scores[:10], 7))) == 1             # ... which are exact ties
    if filt:
        assert (scores == -1).any()                            # the sentinel, last, by index
        assert scores[names.index("p35")] == -1                # stored zeros are not expressing
        tail = names[-int((scores == -1).sum()):]
        assert tail == sorted(tail, key=lambda x: int(x[1:]))
    else:
        assert scores[names.index("p37")] == 0                 # expressed nowhere: 0, not -1
        assert not (scores == -1).any()
    if batch:
        assert set(BATCH[LABELS == "C"]) == {"b3"} and set(LABELS[BATCH == "b3"]) == {"C"}


@pytest.mark.parametrize("path", list(AD_PATHS))
def test_anndata_paths(reference, path):
    (batch, filt), ref = reference
    _same(_ad_result(AD_PATHS[path], batch=batch, remove_lowly_expressed=filt), ref, label=path)


@pytest.mark.parametrize("path", list(AD_GPU_PATHS))
def test_anndata_gpu_code(numpy_gpu, reference, path):  # noqa: F811
    (batch, filt), ref = reference
    if batch and "chunked" in path:
        pytest.skip("batch_key on the GPU takes the batch path whatever gpu_chunk_size is")
    _same(_ad_result(AD_GPU_PATHS[path], batch=batch, remove_lowly_expressed=filt), ref, rtol=1e-5,
          label=path)


@requires_cytome
@pytest.mark.parametrize("path", list(CY_PATHS))
def test_store_paths(stores, reference, path):
    (batch, filt), ref = reference
    if batch and CY_PATHS[path].get("max_accumulator_bytes"):
        pytest.skip("batch_key switches chromosome batching to disk batching")
    _same(_cy_result(stores, path, batch=batch, remove_lowly_expressed=filt), ref, label=path)


@requires_cytome
def test_store_gpu_code(numpy_gpu, stores, reference):  # noqa: F811
    (batch, filt), ref = reference
    _same(_cy_result(stores, "rna_none", batch=batch, device="gpu", remove_lowly_expressed=filt), ref,
          rtol=1e-5, label="store gpu")


@requires_cytome
@pytest.mark.parametrize("path", list(CY_PATHS))
def test_store_dict_output_is_the_positive_part(stores, path):
    """The dict output holds a group's positive scores, the same entries on
    every store path."""
    ref = _cy_result(stores, "rna_none")
    kw = dict(CY_PATHS[path])
    p = stores[kw.pop("store")]
    d = cosg.run_cosg_cytome(p, groupby="grp", n_genes_user=N_TOP, layer="counts", verbose=False,
                             output_format="dict", iqr_normalize=False, **kw)["scores_dict"]
    want = {(g, nm): s for g, (names, scores) in ref.items() for nm, s in zip(names, scores) if s > 0}
    assert set(d) == set(want)
    for key, v in d.items():
        assert v == pytest.approx(want[key], rel=1e-6)


@requires_cytome
@pytest.mark.parametrize("column_major", [False, True])
@pytest.mark.parametrize("fb", ["chromosome", "disk"])
def test_features_not_grouped_by_chromosome(tmp_path, column_major, fb):
    """chr1 split around chr2: the GROUP BY ranges overlapped, 40 features were
    scored twice (duplicate names in a group's markers) and the IQR scale was
    padded to more features than there are. The blocks are now runs of one
    chromosome, and the result equals the unbatched path's."""
    import cytome
    n, f = X.shape
    chroms = ["chr1"] * 10 + ["chr2"] * 10 + ["chr1"] * 15 + ["chrX"] * 5
    p = str(tmp_path / "split.cytome")
    ds = cytome.create(p)
    ds.set_entity("cells", pd.DataFrame({"cell_idx": np.arange(n), "barcode": [f"c{i}" for i in range(n)],
                                         "grp": LABELS, "batch": BATCH}))
    starts = np.arange(f) * 1000
    ds.set_entity("peaks", pd.DataFrame({"peak_idx": np.arange(f), "chr": chroms, "start": starts,
                                         "end_": starts + 500, "peak_id": NAMES}))
    ds.flush()
    w = ds.create_layer_writer("ATAC_counts", n, f, dtype=np.float32, compression="zlib")
    w.write_chunk(X, 0)
    w.finalize()
    ds.flush()
    if column_major:
        _column_major(ds, X, "ATAC_counts", width=7)
    ds.close()
    from cosg._cytome_streaming import _get_chrom_col_ranges
    ds = cytome.open(p, mode="r")
    try:
        ranges, order = _get_chrom_col_ranges(ds, "peaks", "peak_idx")
    finally:
        ds.close()
    assert order == ["chr1", "chr2", "chr1~2", "chrX"]
    assert sorted(ranges.values()) == [(0, 10), (10, 20), (20, 35), (35, 40)]
    stores = {"s": p}
    ref = _cy_result(stores, dict(store="s", modality="ATAC", feature_batching="none"))
    got = _cy_result(stores, dict(store="s", modality="ATAC", feature_batching=fb))
    _same(got, ref, label=fb)
    for g, (names, _) in got.items():
        assert len(set(names)) == len(names), g
    kw = dict(groupby="grp", n_genes_user=N_TOP, layer="counts", verbose=False, modality="ATAC",
              output_format="dense")
    d_ref = cosg.run_cosg_cytome(p, feature_batching="none", **kw)["scores_df"]
    d_got = cosg.run_cosg_cytome(p, feature_batching=fb, **kw)["scores_df"]
    pd.testing.assert_frame_equal(d_got, d_ref, rtol=1e-6)


@requires_cytome
def test_a_column_block_whose_features_all_fail_the_filter(tmp_path):
    """The column scan dropped a block whose features all failed the floor
    before scoring it, and with them its rows of the ndarray output: with
    n_genes_user past the live features the scan returned "" and 0.0 where
    every other path lists those features at -1."""
    import cytome
    from cosg._cytome_streaming import _run_cosg_batched
    rng = np.random.default_rng(4)
    n, f = 120, 30
    lab = np.array(["A", "B"])[rng.integers(0, 2, n)]
    Xd = ((rng.random((n, f)) < 0.6) * rng.integers(1, 4, (n, f))).astype(np.float32)
    Xd[:, 10:25] = 0.0
    Xd[:2, 10:25] = 1.0                                   # under any floor: two cells
    names = [f"chr1:{j * 1000}-{j * 1000 + 500}" for j in range(f)]
    p = str(tmp_path / "dead.cytome")
    ds = cytome.create(p)
    ds.set_entity("cells", pd.DataFrame({"cell_idx": np.arange(n), "barcode": [f"c{i}" for i in range(n)],
                                         "grp": lab}))
    ds.set_entity("peaks", pd.DataFrame({"peak_idx": np.arange(f), "chr": ["chr1"] * f,
                                         "start": np.arange(f) * 1000, "end_": np.arange(f) * 1000 + 500,
                                         "peak_id": names}))
    ds.flush()
    w = ds.create_layer_writer("ATAC_counts", n, f, dtype=np.float32, compression="zlib")
    w.write_chunk(sp.csr_matrix(Xd), 0)
    w.finalize()
    ds.flush()
    _column_major(ds, Xd, "ATAC_counts", width=5)
    ds.close()
    ref = cosg.run_cosg_cytome(p, groupby="grp", modality="ATAC", layer="counts", n_genes_user=40,
                               feature_batching="none", verbose=False)
    got = _run_cosg_batched(p, "grp", 40, 1.0, True, 0.1, 3, None, 2048, False, "ATAC", False, 0.95, 0.75,
                            0, "chromosome", layer="counts", block_features=5)
    for k in range(2):
        assert list(got["names"][:, k]) == list(ref["names"][:, k])
        np.testing.assert_allclose(got["scores"][:, k], ref["scores"][:, k], rtol=1e-6)
    assert "" not in set(got["names"].ravel())


@requires_cytome
@pytest.mark.parametrize("fb", ["none", "chromosome", "column_major"])
def test_feature_name_col_names_features_on_every_store_path(stores, tmp_path, fb):
    """The feature-batched path read the default name column for the dict,
    long and ndarray outputs, and feature_name_col only for the dense one."""
    import cytome
    import shutil
    src = stores["atac_col" if fb == "column_major" else "atac"]
    p = str(tmp_path / "named.cytome")
    shutil.copyfile(src, p)
    ds = cytome.open(p)
    ds.peaks["nearest_gene"] = np.array([f"G{j}" for j in range(X.shape[1])], dtype=object)
    ds.flush()
    ds.close()
    kw = dict(groupby="grp", modality="ATAC", layer="counts", verbose=False, n_genes_user=N_TOP,
              feature_name_col="nearest_gene",
              feature_batching="none" if fb == "none" else "chromosome")
    d = cosg.run_cosg_cytome(p, output_format="dict", **kw)["scores_dict"]
    nd = cosg.run_cosg_cytome(p, **kw)["names"]
    assert {f for _, f in d} <= {f"G{j}" for j in range(X.shape[1])}
    assert set(nd.ravel()) <= {f"G{j}" for j in range(X.shape[1])}
