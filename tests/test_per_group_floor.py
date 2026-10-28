"""A detection floor per group: ``expressed_pct`` as a dict, a Series or a sequence.

The floor removes a feature from one group's markers only: it is applied to the
group's own column of scores after the cosine and COSG's penalty, which use every
group but not the floor. So, by design and on every path:

- a scalar and the same value for every group give identical results;
- the result of group g under a per-group floor equals g's result from a run with
  floor[g] as the scalar.

A floor of 0 means "the 3-cell minimum only" on every path (before, it switched
the filter off altogether on the streaming paths); ``remove_lowly_expressed=False``
is the switch. ``test_zero_floor_and_slash_labels.py`` tests that meaning itself.

One thing couples groups, not in the floor itself, and is tested as what it is
rather than as the equivalence: ``scale_rule='positive_median'``
(``iqr_normalize``): a group with fewer than 20 positive scores is scaled by the
median of the other groups' scales, which their floors move. Its raw scores and
its ranking still follow its own floor alone. (``min_total_count`` on the
feature-batched cytome path coupled them too while a scalar 0 dropped it; it now
applies whatever the floors.)

The GPU paths are run here with numpy standing in for CuPy (their arithmetic,
not the device's), the threshold vector they build is checked against the
expression the device evaluates for a scalar, and they run on the device when
CuPy and a GPU are present.
"""
from __future__ import annotations

import sys
import types

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import cosg
from cosg import _floor
from conftest import requires_cytome

GROUPS = np.array(["a", "b", "c", "d", "e"])
#: Floors that cut each group's weak markers at a different place, and 0 for d,
#: which leaves d the 3-cell minimum alone.
FLOORS = {"a": 0.3, "b": 0.05, "c": 0.2, "d": 0.0, "e": 0.12}
N_TOP = 40


def _data(seed=11, n=700, f=200):
    """Five groups of 35 % to 5 % of the cells, each with 20 strong markers
    (detected in half its cells) and 20 weak ones, detected in 3 % to 35 % of
    its cells, so that the floors above cut different numbers of them."""
    rng = np.random.default_rng(seed)
    gi = rng.choice(len(GROUPS), n, p=[0.35, 0.25, 0.2, 0.15, 0.05])
    rate = np.full((len(GROUPS), f), 0.02)
    for g in range(len(GROUPS)):
        rate[g, g * 40:g * 40 + 20] = 0.5
        rate[g, g * 40 + 20:g * 40 + 40] = np.linspace(0.03, 0.35, 20)
    X = ((rng.random((n, f)) < rate[gi]) * rng.integers(1, 6, (n, f))).astype(np.float32)
    batch = np.where(np.arange(n) % 3 == 0, "B1", "B2")
    return X, GROUPS[gi], batch


X, LABELS, BATCH = _data()


# ── the resolver and the thresholds ────────────────────────────────────────

def test_a_number_is_one_floor_for_every_group():
    fl, per = _floor.resolve_expressed_pct(0.05, ["a", "b"])
    assert not per and fl.tolist() == [0.05, 0.05]
    fl, per = _floor.resolve_expressed_pct(None, ["a", "b"])
    assert not per and fl.tolist() == [0.1, 0.1]
    fl, per = _floor.resolve_expressed_pct(np.float32(0.25), ["a"])
    assert not per and fl.tolist() == [0.25]


def test_dict_series_and_sequence_give_the_same_floors():
    want = [FLOORS[g] for g in GROUPS]
    for given in (FLOORS, pd.Series(FLOORS), pd.Series(FLOORS).iloc[::-1], want, tuple(want),
                  np.asarray(want), pd.Index(want)):
        fl, per = _floor.resolve_expressed_pct(given, list(GROUPS))
        assert per and fl.dtype == np.float64 and fl.tolist() == want, type(given)


def test_keys_match_labels_by_their_string_form():
    fl, _ = _floor.resolve_expressed_pct({"0": 0.2, 1: 0.3}, [0, 1])
    assert fl.tolist() == [0.2, 0.3]
    fl, _ = _floor.resolve_expressed_pct(pd.Series({0: 0.2, 1: 0.3}), ["0", "1"])
    assert fl.tolist() == [0.2, 0.3]


def test_a_label_not_scored_here_is_ignored_and_a_group_without_cells_needs_none():
    fl, _ = _floor.resolve_expressed_pct({"a": 0.2, "b": 0.3, "z": 0.5}, ["a", "b"],
                                         known_labels=["a", "b", "z"])
    assert fl.tolist() == [0.2, 0.3]
    fl, _ = _floor.resolve_expressed_pct({"a": 0.2}, ["a", "b"], group_sizes=[10, 0])
    assert fl.tolist() == [0.2, _floor.DEFAULT_EXPRESSED_PCT]


@pytest.mark.parametrize("given,error,match", [
    ({"a": 0.1, "b": 0.1, "zz": 0.1}, ValueError, r"not groups of groupby: \['zz'\]"),
    ({"a": 0.1}, ValueError, r"no fraction for 1 group\(s\): \['b'\]"),
    ({"a": 0.1, "b": 1.5}, ValueError, r"between 0 and 1 .*'b': 1.5"),
    ({"a": 0.1, "b": -0.01}, ValueError, "between 0 and 1"),
    ({"a": 0.1, "b": float("nan")}, ValueError, "between 0 and 1"),
    ({"a": 0.1, "b": "high"}, TypeError, "must be numbers"),
    ({"a": 0.1, "b": 0.1, "0": 0.2}, ValueError, "not groups of groupby"),
    (pd.Series([0.1, 0.2, 0.3], index=["a", "b", "a"]), ValueError, "more than once"),
    ([0.1], ValueError, "1 fractions for 2 groups"),
    ([0.1, 0.2, 0.3], ValueError, "3 fractions for 2 groups"),
    ([0.1, 2.0], ValueError, "between 0 and 1"),
    (1.5, ValueError, "between 0 and 1"),
    (-0.1, ValueError, "between 0 and 1"),
    (float("nan"), ValueError, "between 0 and 1"),
    ("0.1", TypeError, "got str"),
    ({0.1, 0.2}, TypeError, "set"),
    (pd.DataFrame({"f": [0.1, 0.2]}, index=["a", "b"]), TypeError, "DataFrame"),
    (np.full((2, 1), 0.1), TypeError, "2-D"),
])
def test_bad_floors_raise_clear_errors(given, error, match):
    with pytest.raises(error, match=match):
        _floor.resolve_expressed_pct(given, ["a", "b"], known_labels=["a", "b"])


def test_two_keys_for_one_group_raise():
    with pytest.raises(ValueError, match="two fractions"):
        _floor.resolve_expressed_pct({0: 0.1, "0": 0.2, 1: 0.1}, [0, 1])


def test_a_key_matching_two_labels_by_string_raises():
    with pytest.raises(ValueError, match="more than one group"):
        _floor.resolve_expressed_pct({"k": 0.1}, [_Str("k", 1), _Str("k", 2)])


class _Str:
    """Two distinct labels whose string forms coincide."""
    def __init__(self, s, i):
        self.s, self.i = s, i

    def __str__(self):
        return self.s

    def __hash__(self):
        return hash((self.s, self.i))

    def __eq__(self, other):
        return isinstance(other, _Str) and (self.s, self.i) == (other.s, other.i)


SIZES = np.arange(0, 3001, dtype=np.float64)
FRACTIONS = (0.1, 0.05, 0.02, 0.01, 0.2, 0.3, 1 / 3, 0.0165289, 0.00527704, 1.0, 0.0)


@pytest.mark.parametrize("pct", FRACTIONS)
def test_the_float64_vector_reproduces_the_cpu_scalar_arithmetic(pct):
    got = _floor.floor_thresholds(SIZES, np.full(SIZES.size, pct), 3)
    sizes_int = SIZES.astype(np.int64)
    # the legacy loop: max(np.int64 size * float, 3); the streaming paths: float64 sizes
    want = np.array([max(n * pct, 3) for n in sizes_int], dtype=np.float64)
    np.testing.assert_array_equal(got, want)
    np.testing.assert_array_equal(got, np.maximum(SIZES * pct, 3))


@pytest.mark.parametrize("pct", FRACTIONS)
def test_the_gpu_threshold_vector_reproduces_the_device_scalar_arithmetic(pct):
    """What the GPU paths build, without a GPU. For a scalar floor the device
    formed ``float32 sizes * python float``, which CuPy (like numpy on an array)
    evaluates in float32, then the maximum with the 3-cell minimum; the vector
    of one fraction per group must give those thresholds to the bit."""
    sizes32 = SIZES.astype(np.float32)
    device = np.maximum(sizes32 * pct, 3)                     # cupy's expression, on numpy
    assert device.dtype == np.float32
    mono = _floor.floor_thresholds(sizes32, np.full(SIZES.size, pct), 3, dtype=np.float32)
    assert mono.dtype == np.float32
    np.testing.assert_array_equal(mono, device)                # gpu_core_cosg
    # gpu_core_cosg_chunked and the cytome GPU path: one group at a time,
    # max(float(size * pct), 3.0) and nothing for an empty group; a floor of 0
    # keeps the 3-cell minimum there too (it switched the filter off before)
    chunked = _floor.floor_thresholds(sizes32, np.full(SIZES.size, pct), 3, dtype=np.float32,
                                      skip_empty=True)
    per_group = [0.0 if s <= 0 else max(float(s * np.float32(pct)), 3.0) for s in sizes32]
    np.testing.assert_array_equal(chunked.astype(np.float64), np.asarray(per_group))


def test_the_gpu_threshold_vector_holds_each_groups_own_floor():
    sizes = np.array([245, 175, 140, 105, 35, 0], dtype=np.float32)
    floors = np.array([0.3, 0.05, 0.2, 0.0, 0.12, 0.5])
    got = _floor.floor_thresholds(sizes, floors, 3, dtype=np.float32)
    for k in range(len(sizes)):
        one = _floor.floor_thresholds(sizes[k:k + 1], floors[k:k + 1], 3, dtype=np.float32)
        assert got[k] == one[0] == np.maximum(sizes[k] * np.float32(floors[k]), np.float32(3))
    gated = _floor.floor_thresholds(sizes, floors, 3, dtype=np.float32, skip_empty=True)
    assert gated[3] == 3 and gated[5] == 0                      # floor 0: the minimum; no cells
    np.testing.assert_array_equal(gated[:5], got[:5])


def test_floors_for_columns_takes_numbers_and_sequences_only():
    assert _floor.floors_for_columns(0.2, 3).tolist() == [0.2, 0.2, 0.2]
    assert _floor.floors_for_columns([0.1, 0.2, 0.3], 3).tolist() == [0.1, 0.2, 0.3]
    with pytest.raises(TypeError, match="keyed by group"):
        _floor.floors_for_columns({"a": 0.1}, 1)
    with pytest.raises(ValueError, match="2 fractions for 3 groups"):
        _floor.floors_for_columns([0.1, 0.2], 3)


# ── AnnData ────────────────────────────────────────────────────────────────

def _adata(fmt, labels=LABELS):
    import anndata as ad
    mat = {"csr": sp.csr_matrix, "csc": sp.csc_matrix, "dense": np.array}[fmt](X)
    a = ad.AnnData(mat, obs=pd.DataFrame({"grp": pd.Categorical(labels), "batch": BATCH},
                                         index=[f"c{i}" for i in range(len(labels))]))
    a.var_names = [f"g{j}" for j in range(X.shape[1])]
    return a


AD_PATHS = {
    "csr_stream": dict(fmt="csr"),
    "csc_stream": dict(fmt="csc"),
    "dense_stream": dict(fmt="dense"),
    "csr_stream_chunks": dict(fmt="csr", cpu_chunk_size=64),
    "csr_legacy": dict(fmt="csr", cpu_chunk_size=0),
    "csc_legacy": dict(fmt="csc", cpu_chunk_size=0),
    "dense_legacy": dict(fmt="dense", cpu_chunk_size=0),
    "csr_batch": dict(fmt="csr", batch_key="batch"),
    "csr_batch_legacy": dict(fmt="csr", batch_key="batch", cpu_chunk_size=0),
    "dense_batch_legacy": dict(fmt="dense", batch_key="batch", cpu_chunk_size=0),
    "csr_pvalues": dict(fmt="csr", calculate_pvalues=True),
}


def _ad_run(path, expressed_pct, **extra):
    kw = dict(AD_PATHS[path]) if isinstance(path, str) else dict(path)
    a = _adata(kw.pop("fmt"))
    base = dict(groupby="grp", n_genes_user=N_TOP, mu=1)
    base.update(kw)
    base.update(extra)
    cosg.cosg(a, expressed_pct=expressed_pct, **base)
    return a.uns["cosg"]


def _ad_columns(uns):
    """{field: {group: array}} for every per-group field of the result."""
    return {k: {g: np.asarray(uns[k][g]) for g in uns[k].dtype.names}
            for k in uns if k not in ("params", "COSG")}


def _ad_equal(r1, r2, group=None):
    c1, c2 = _ad_columns(r1), _ad_columns(r2)
    assert set(c1) == set(c2)
    for k in c1:
        for g in ([group] if group is not None else c1[k]):
            a, b = c1[k][g], c2[k][g]
            if a.dtype.kind in "fc":
                assert np.array_equal(a, b, equal_nan=True), (k, g)
            else:
                assert np.array_equal(a, b), (k, g)


@pytest.mark.parametrize("path", list(AD_PATHS))
@pytest.mark.parametrize("pct", [0.1, 0.0, 0.3])
def test_anndata_a_scalar_and_a_uniform_vector_are_identical(path, pct):
    ref = _ad_run(path, pct)
    for given in ([pct] * 5, {g: pct for g in GROUPS}, pd.Series(pct, index=GROUPS)):
        _ad_equal(ref, _ad_run(path, given))
    assert pd.DataFrame(ref["COSG"]).equals(pd.DataFrame(_ad_run(path, [pct] * 5)["COSG"]))


@pytest.mark.parametrize("path", list(AD_PATHS))
def test_anndata_each_group_gets_its_own_floor(path):
    per_group = _ad_run(path, FLOORS)
    uniform = _ad_columns(_ad_run(path, 0.1))
    changed = 0
    for g, f in FLOORS.items():
        _ad_equal(per_group, _ad_run(path, f), group=g)
        changed += not np.array_equal(_ad_columns(per_group)["scores"][g], uniform["scores"][g])
    assert changed >= 3, "the floors chosen do not change the result; the test proves nothing"


def test_anndata_records_the_floor_of_each_group(tmp_path):
    import anndata as ad
    uns = _ad_run("csr_stream", FLOORS)
    rec = uns["params"]["expressed_pct"]
    assert pd.DataFrame(rec).iloc[0].to_dict() == FLOORS
    assert _ad_run("csr_stream", 0.05)["params"]["expressed_pct"] == 0.05      # a number, as given
    # labels holding '/' (L2/3 IT) survive adata.write, which a dict keyed by them
    # would not (anndata splits the key into groups, or refuses it); the 'COSG'
    # frame writes '/' as '%2F' (test_zero_floor_and_slash_labels.py)
    labels = np.where(LABELS == "a", "L2/3 IT", LABELS)
    a = _adata("csr", labels)
    floors = {("L2/3 IT" if g == "a" else g): f for g, f in FLOORS.items()}
    cosg.cosg(a, groupby="grp", n_genes_user=10, expressed_pct=floors)
    a.write_h5ad(tmp_path / "x.h5ad")
    back = ad.read_h5ad(tmp_path / "x.h5ad").uns["cosg"]["params"]["expressed_pct"]
    assert pd.DataFrame(back).iloc[0].to_dict() == floors


def test_anndata_groups_subset_needs_floors_for_the_selected_groups_only():
    # with a categorical groupby the result holds the selected groups only; a
    # key for an unselected category is accepted, and needs no fraction
    sel = ["a", "c", "e"]
    got = _ad_run(dict(fmt="csr"), {"a": 0.3, "c": 0.2, "e": 0.12}, groups=sel)
    also = _ad_run(dict(fmt="csr"), {"a": 0.3, "c": 0.2, "e": 0.12, "b": 0.9}, groups=sel)
    for g in sel:
        _ad_equal(got, also, group=g)
        _ad_equal(got, _ad_run(dict(fmt="csr"), FLOORS[g], groups=sel), group=g)


def test_anndata_int_labels_take_string_keys():
    a = _adata("csr", np.searchsorted(GROUPS, LABELS))           # categories 0..4
    cosg.cosg(a, groupby="grp", n_genes_user=10, expressed_pct={str(i): FLOORS[g] for i, g in enumerate(GROUPS)})
    assert pd.DataFrame(a.uns["cosg"]["params"]["expressed_pct"]).iloc[0].tolist() == list(FLOORS.values())


@pytest.mark.parametrize("given,error,match", [
    ([0.1] * 4, ValueError, "4 fractions for 5 groups"),
    ({**FLOORS, "zz": 0.1}, ValueError, r"not groups of groupby: \['zz'\]"),
    ({"a": 0.1}, ValueError, "no fraction for 4 group"),
    ({**FLOORS, "c": 1.2}, ValueError, "between 0 and 1"),
    (1.2, ValueError, "between 0 and 1"),
    ({1, 2}, TypeError, "set"),
])
def test_anndata_bad_floors_raise_before_anything_is_written(given, error, match):
    a = _adata("csr")
    with pytest.raises(error, match=match):
        cosg.cosg(a, groupby="grp", expressed_pct=given)
    assert "cosg" not in a.uns


# ── score_from_sums ────────────────────────────────────────────────────────

def _sums():
    onehot = (LABELS[:, None] == GROUPS[None, :]).astype(np.float64)
    Xd = X.astype(np.float64)
    return Xd.T @ onehot, (Xd ** 2).sum(axis=0), ((Xd > 0).T @ onehot).astype(np.int64), onehot.sum(axis=0)


def test_score_from_sums_takes_a_floor_per_column():
    dot, sq, nnz, sizes = _sums()
    vec = [FLOORS[g] for g in GROUPS]
    for pct in (0.1, 0.0, 0.3):
        np.testing.assert_array_equal(cosg.score_from_sums(dot, sq, sizes, nnz, expressed_pct=pct),
                                      cosg.score_from_sums(dot, sq, sizes, nnz, expressed_pct=[pct] * 5))
    got = cosg.score_from_sums(dot, sq, sizes, nnz, expressed_pct=vec, min_total_count=10)
    for k, f in enumerate(vec):
        one = cosg.score_from_sums(dot, sq, sizes, nnz, expressed_pct=f, min_total_count=10)
        np.testing.assert_array_equal(got[:, k], one[:, k])
    # here a floor of 0 keeps the 3-cell minimum, as a scalar 0 always has
    thr = np.maximum(sizes * np.asarray(vec), 3)
    np.testing.assert_array_equal(cosg.score_from_sums(dot, sq, sizes, nnz, expressed_pct=vec) == -1,
                                  nnz < thr[None, :])
    with pytest.raises(TypeError, match="keyed by group"):
        cosg.score_from_sums(dot, sq, sizes, nnz, expressed_pct=FLOORS)
    with pytest.raises(ValueError, match="4 fractions for 5 groups"):
        cosg.score_from_sums(dot, sq, sizes, nnz, expressed_pct=vec[:4])


# ── cytome ─────────────────────────────────────────────────────────────────

CHROMS = ["chr1"] * 80 + ["chr2"] * 85 + ["chrX"] * 35


def _column_major(ds, Xm, name, width):
    """Rewrite a stored matrix as column chunks (as tests/test_column_major_scan.py does)."""
    csc = sp.csc_matrix(Xm)
    dtype = np.dtype(ds._conn.execute(
        "SELECT dtype FROM matrix_meta WHERE matrix_name=?", (name,)).fetchone()[0])
    ds._conn.execute("DELETE FROM matrix_chunks WHERE matrix_name=?", (name,))
    ds._conn.execute("DELETE FROM matrix_csc_chunks WHERE matrix_name=?", (name,))
    k = 0
    for k, start in enumerate(range(0, Xm.shape[1], width)):
        end = min(start + width, Xm.shape[1])
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


def _store(path, modality, column_major=False):
    import cytome
    n, f = X.shape
    ds = cytome.create(str(path))
    ds.set_entity("cells", pd.DataFrame({"cell_idx": np.arange(n), "barcode": [f"BC{i}" for i in range(n)],
                                         "grp": LABELS, "batch": BATCH}))
    if modality == "RNA":
        ds.set_entity("genes", pd.DataFrame({"gene_idx": np.arange(f), "gene_id": [f"g{j}" for j in range(f)]}))
    else:
        starts = np.concatenate([np.arange(CHROMS.count(c)) for c in ("chr1", "chr2", "chrX")]) * 1000
        ds.set_entity("peaks", pd.DataFrame({"peak_idx": np.arange(f), "chr": CHROMS, "start": starts,
                                             "end_": starts + 500, "peak_id": [f"g{j}" for j in range(f)]}))
    ds.flush()
    # zlib, which the standard library has, so the stores can be written where
    # zstandard is not installed (the GPU environment)
    w = ds.create_layer_writer(f"{modality}_counts", n, f, dtype=np.float32, compression="zlib")
    w.write_chunk(sp.csr_matrix(X), 0)
    w.finalize()
    ds.flush()
    if column_major:
        _column_major(ds, X, f"{modality}_counts", width=23)       # chunks straddle chromosomes
    ds.close()
    return str(path)


@pytest.fixture(scope="module")
def stores(tmp_path_factory):
    pytest.importorskip("cytome")
    d = tmp_path_factory.mktemp("per_group_floor")
    return {"rna": _store(d / "rna.cytome", "RNA"), "atac": _store(d / "atac.cytome", "ATAC"),
            "atac_col": _store(d / "atac_col.cytome", "ATAC", column_major=True)}


CY_PATHS = {
    "rna_ndarray": dict(store="rna", modality="RNA", feature_batching="none"),
    "rna_dict": dict(store="rna", modality="RNA", feature_batching="none", output_format="dict"),
    "rna_long": dict(store="rna", modality="RNA", feature_batching="none", output_format="long"),
    "rna_batch_key": dict(store="rna", modality="RNA", feature_batching="none", batch_key="batch"),
    "rna_pvalues": dict(store="rna", modality="RNA", feature_batching="none", calculate_pvalues=True),
    "rna_mask": dict(store="rna", modality="RNA", feature_batching="none", output_format="dict",
                     cell_mask=np.arange(len(LABELS)) % 5 != 0),
    "atac_chromosome": dict(store="atac", modality="ATAC", feature_batching="chromosome", output_format="dict"),
    "atac_chromosome_ndarray": dict(store="atac", modality="ATAC", feature_batching="chromosome"),
    "atac_one_pass_per_chromosome": dict(store="atac", modality="ATAC", feature_batching="chromosome",
                                         output_format="dict", max_accumulator_bytes=1),
    "atac_disk": dict(store="atac", modality="ATAC", feature_batching="disk", output_format="dict"),
    "atac_disk_batch_key": dict(store="atac", modality="ATAC", feature_batching="disk", output_format="dict",
                                batch_key="batch"),
    "atac_column_major": dict(store="atac_col", modality="ATAC", feature_batching="chromosome",
                              output_format="dict"),
    "atac_column_major_ndarray": dict(store="atac_col", modality="ATAC", feature_batching="chromosome"),
    "atac_column_major_batch_key": dict(store="atac_col", modality="ATAC", feature_batching="chromosome",
                                        output_format="dict", batch_key="batch"),
    "atac_column_major_min_total_count": dict(store="atac_col", modality="ATAC", feature_batching="chromosome",
                                              output_format="dict", min_total_count=12),
    "atac_chromosome_min_total_count": dict(store="atac", modality="ATAC", feature_batching="chromosome",
                                            output_format="dict", min_total_count=12),
}


def _cy_run(stores, path, expressed_pct, **extra):
    kw = dict(CY_PATHS[path])
    p = stores[kw.pop("store")]
    kw.update(extra)
    base = dict(groupby="grp", n_genes_user=N_TOP, mu=1.0, layer="counts", verbose=False)
    base.update(kw)
    return cosg.run_cosg_cytome(p, expressed_pct=expressed_pct, **base)


def _cy_group(res, g):
    """What the result says about group g, in a comparable form."""
    if "names" in res:
        k = [str(x) for x in res["groups_order"]].index(g)
        out = {"names": res["names"][:, k].tolist(), "scores": res["scores"][:, k].tolist()}
        for key in ("pvals", "pvals_adj", "zscores"):
            if key in res:
                out[key] = np.nan_to_num(res[key][:, k], nan=-7.0).tolist()
        return out
    if "scores_dict" in res:
        return {f: s for (h, f), s in res["scores_dict"].items() if h == g}
    df = res["scores_df"]
    return df[df["group"] == g].reset_index(drop=True).to_dict("list")


def _cy_equal(r1, r2):
    """The same result; a per-group run's 'expressed_pct' record aside."""
    a = {k: v for k, v in r1.items() if k != "expressed_pct"}
    b = {k: v for k, v in r2.items() if k != "expressed_pct"}
    assert set(a) == set(b)
    for k in a:
        if isinstance(a[k], pd.DataFrame):
            pd.testing.assert_frame_equal(a[k], b[k])
        elif isinstance(a[k], np.ndarray):
            assert np.array_equal(a[k], b[k], equal_nan=a[k].dtype.kind == "f"), k
        else:
            assert a[k] == b[k], k


@requires_cytome
@pytest.mark.parametrize("path", list(CY_PATHS))
@pytest.mark.parametrize("pct", [0.1, 0.0, 0.3])
def test_cytome_a_scalar_and_a_uniform_vector_are_identical(stores, path, pct):
    ref = _cy_run(stores, path, pct)
    assert "expressed_pct" not in ref                 # a scalar run's result is unchanged
    for given in ([pct] * 5, {g: pct for g in GROUPS}):
        got = _cy_run(stores, path, given)
        _cy_equal(ref, got)
        assert set(got["expressed_pct"]) == {str(g) for g in got["groups_order"]}


@requires_cytome
@pytest.mark.parametrize("path", list(CY_PATHS))
def test_cytome_each_group_gets_its_own_floor(stores, path):
    per_group = _cy_run(stores, path, FLOORS)
    assert {str(k): v for k, v in per_group["expressed_pct"].items()} == FLOORS
    uniform = _cy_run(stores, path, 0.1)
    changed = 0
    for g, f in FLOORS.items():
        assert _cy_group(per_group, g) == _cy_group(_cy_run(stores, path, f), g), g
        changed += _cy_group(per_group, g) != _cy_group(uniform, g)
    assert changed >= 3, "the floors chosen do not change the result; the test proves nothing"


@requires_cytome
@pytest.mark.parametrize("path", ["atac_chromosome_min_total_count", "atac_column_major_min_total_count"])
def test_cytome_a_zero_floor_with_min_total_count(stores, path):
    """On the feature-batched path a group's floor of 0 is the scalar 0 with
    ``min_total_count`` too. (Before, a scalar 0 kept no nonzero counts and so
    dropped ``min_total_count``, while a group's 0 among positive floors got it;
    ``min_total_count`` now applies whatever the floors.) Its result does not
    depend on the other groups' floors, and ``min_total_count`` reaches it."""
    a = _cy_group(_cy_run(stores, path, FLOORS), "d")
    b = _cy_group(_cy_run(stores, path, {**FLOORS, "a": 0.9, "b": 0.01, "c": 0.5, "e": 0.3}), "d")
    assert a == b and len(a) > 0
    assert a == _cy_group(_cy_run(stores, path, 0.0), "d")
    # min_total_count does reach d (it removes features from the scale of its
    # normalized scores), at a scalar 0 as well
    no_mtc = _cy_group(_cy_run(stores, path, FLOORS, min_total_count=0), "d")
    assert a != no_mtc
    assert no_mtc == _cy_group(_cy_run(stores, path, 0.0, min_total_count=0), "d")


@requires_cytome
def test_cytome_dense_output_scalar_and_uniform_vector(stores):
    a = _cy_run(stores, "rna_dict", 0.1, output_format="dense")
    b = _cy_run(stores, "rna_dict", pd.Series(0.1, index=GROUPS), output_format="dense")
    pd.testing.assert_frame_equal(a["scores_df"], b["scores_df"])


@requires_cytome
def test_cytome_a_table_of_every_type_serves_any_mask(stores):
    # cell_mask removes every cell of e; its floor in the table is ignored
    mask = LABELS != "e"
    got = cosg.run_cosg_cytome(stores["rna"], groupby="grp", layer="counts", modality="RNA",
                               feature_batching="none", output_format="dict", verbose=False,
                               cell_mask=mask, expressed_pct=FLOORS)
    assert "e" not in got["groups_order"] and set(got["expressed_pct"]) == {"a", "b", "c", "d"}


@requires_cytome
def test_cytome_records_the_floors_with_the_markers(tmp_path):
    import cytome
    p = _store(tmp_path / "w.cytome", "RNA")      # its own store: this test writes to it
    cosg.run_cosg_cytome(str(p), groupby="grp", layer="counts", modality="RNA", verbose=False,
                         expressed_pct=FLOORS, write_to_cytome=True, key_added="cosg_floor")
    ds = cytome.open(str(p))
    try:
        assert ds.metadata["cosg_floor"]["expressed_pct"] == FLOORS
    finally:
        ds.close()


@requires_cytome
@pytest.mark.parametrize("given,error,match", [
    ([0.1] * 4, ValueError, "4 fractions for 5 groups"),
    ({**FLOORS, "zz": 0.1}, ValueError, r"not groups of groupby: \['zz'\]"),
    ({"a": 0.1, "b": 0.1}, ValueError, "no fraction for 3 group"),
    ({**FLOORS, "c": 1.2}, ValueError, "between 0 and 1"),
    (-0.5, ValueError, "between 0 and 1"),
])
@pytest.mark.parametrize("path", ["rna_ndarray", "atac_chromosome", "atac_column_major"])
def test_cytome_bad_floors_raise(stores, given, error, match, path):
    with pytest.raises(error, match=match):
        _cy_run(stores, path, given)


@requires_cytome
def test_cytome_positive_median_scale_couples_groups_only_through_the_scale(stores):
    """``scale_rule='positive_median'`` gives a group with fewer than 20 positive
    scores the median of the other groups' scales, so its normalized scores
    depend on their floors: there the equivalence holds for its features and
    their order, and for its raw scores, not for the normalized values."""
    kw = dict(scale_rule="positive_median", n_genes_user=200)
    floors = {**FLOORS, "e": 0.5}                    # e keeps fewer than 20 positive scores
    per_group = _cy_run(stores, "atac_chromosome", floors, **kw)
    alone = _cy_run(stores, "atac_chromosome", 0.5, **kw)
    e_pg, e_alone = _cy_group(per_group, "e"), _cy_group(alone, "e")
    assert 0 < len(e_pg) < 20
    assert sorted(e_pg, key=e_pg.get) == sorted(e_alone, key=e_alone.get)
    assert e_pg != e_alone                           # the stand-in scale moved with the others' floors
    raw = dict(kw, iqr_normalize=False)
    assert _cy_group(_cy_run(stores, "atac_chromosome", floors, **raw), "e") == \
        _cy_group(_cy_run(stores, "atac_chromosome", 0.5, **raw), "e")
    # and the groups with 20 or more positive scores keep their own scale
    for g in ("a", "c"):
        assert _cy_group(per_group, g) == _cy_group(_cy_run(stores, "atac_chromosome", floors[g], **kw), g)


# ── GPU paths, with numpy standing in for CuPy ─────────────────────────────

def _numpy_as_cupy():
    """numpy and scipy.sparse dressed as cupy and cupyx.scipy.sparse: enough for
    cosg._gpu and the cytome GPU path to run their own code on the CPU."""
    cp = types.ModuleType("cupy")
    for name in dir(np):
        if not name.startswith("__"):
            setattr(cp, name, getattr(np, name))
    cp.asnumpy = lambda a: np.asarray(a)
    cp.get_default_memory_pool = lambda: types.SimpleNamespace(free_all_blocks=lambda: None)
    cp.cuda = types.SimpleNamespace(
        Stream=types.SimpleNamespace(null=types.SimpleNamespace(synchronize=lambda: None)),
        Device=lambda *a, **k: types.SimpleNamespace(mem_info=(8 << 30, 16 << 30)),
        runtime=types.SimpleNamespace(getDeviceCount=lambda: 1))
    cupyx = types.ModuleType("cupyx")
    cupyx_scipy = types.ModuleType("cupyx.scipy")
    cupyx.scipy = cupyx_scipy
    cupyx_scipy.sparse = sp
    return {"cupy": cp, "cupyx": cupyx, "cupyx.scipy": cupyx_scipy, "cupyx.scipy.sparse": sp}


@pytest.fixture
def numpy_gpu(monkeypatch):
    import cosg._backend as backend
    for name, mod in _numpy_as_cupy().items():
        monkeypatch.setitem(sys.modules, name, mod)
    monkeypatch.delitem(sys.modules, "cosg._gpu", raising=False)   # re-imported on the stand-in
    monkeypatch.setattr(cosg, "_gpu", None, raising=False)
    monkeypatch.setattr(backend, "HAS_CUPY", True)
    yield


GPU_PATHS = {
    "csr_monolithic": dict(fmt="csr", device="gpu"),
    "dense_monolithic": dict(fmt="dense", device="gpu"),
    "csr_chunked": dict(fmt="csr", device="gpu", gpu_chunk_size=128),
    "dense_chunked": dict(fmt="dense", device="gpu", gpu_chunk_size=128),
    "csr_batch": dict(fmt="csr", device="gpu", batch_key="batch"),
    "dense_batch": dict(fmt="dense", device="gpu", batch_key="batch"),
}


def _removed(uns):
    """Per group, the genes the floor removed (scored -1)."""
    c = _ad_columns(uns)
    return {g: set(c["names"][g][c["scores"][g] == -1]) for g in c["names"]}


def _gpu_batch_removes_what_the_cpu_batch_path_removes():
    """The GPU batch path counts no nonzeros; the floor is applied after it,
    as after the CPU's batch path, and it used to be skipped there: with
    batch_key on the GPU, no gene was ever removed."""
    every = dict(n_genes_user=X.shape[1])
    for fmt in ("csr", "dense"):
        gpu = _removed(_ad_run(dict(fmt=fmt, device="gpu", batch_key="batch"), FLOORS, **every))
        cpu = _removed(_ad_run(dict(fmt=fmt, batch_key="batch"), FLOORS, **every))
        assert gpu == cpu and all(gpu[g] for g in ("a", "b", "c", "e"))
        off = _removed(_ad_run(dict(fmt=fmt, device="gpu", batch_key="batch"), FLOORS,
                               remove_lowly_expressed=False, **every))
        assert not any(off.values())


@pytest.mark.parametrize("path", list(GPU_PATHS))
def test_gpu_code_scalar_and_uniform_vector_on_numpy(numpy_gpu, path):
    for pct in (0.1, 0.0, 0.3):
        ref = _ad_run(GPU_PATHS[path], pct)
        _ad_equal(ref, _ad_run(GPU_PATHS[path], [pct] * 5))
        _ad_equal(ref, _ad_run(GPU_PATHS[path], {g: pct for g in GROUPS}))
    assert sys.modules["cosg._gpu"].cp is sys.modules["cupy"]     # the GPU code ran, on the stand-in


@pytest.mark.parametrize("path", list(GPU_PATHS))
def test_gpu_code_each_group_gets_its_own_floor_on_numpy(numpy_gpu, path):
    per_group = _ad_run(GPU_PATHS[path], FLOORS)
    for g, f in FLOORS.items():
        _ad_equal(per_group, _ad_run(GPU_PATHS[path], f), group=g)


def test_gpu_batch_code_applies_the_floor_on_numpy(numpy_gpu):
    _gpu_batch_removes_what_the_cpu_batch_path_removes()


@requires_cytome
@pytest.mark.parametrize("batch_key", [None, "batch"])
def test_cytome_gpu_code_on_numpy(numpy_gpu, stores, batch_key):
    run = lambda pct: cosg.run_cosg_cytome(stores["rna"], groupby="grp", n_genes_user=N_TOP, layer="counts",
                                           modality="RNA", verbose=False, device="gpu", batch_key=batch_key,
                                           expressed_pct=pct)
    for pct in (0.1, 0.0):
        _cy_equal(run(pct), run([pct] * 5))
    per_group = run(FLOORS)
    assert per_group["expressed_pct"] == FLOORS
    for g, f in FLOORS.items():
        assert _cy_group(per_group, g) == _cy_group(run(f), g)


# ── GPU paths on the device ────────────────────────────────────────────────

def _has_gpu():
    try:
        import cupy
        return cupy.cuda.runtime.getDeviceCount() > 0
    except Exception:
        return False


requires_gpu = pytest.mark.skipif(not _has_gpu(), reason="no cupy / no CUDA device")


def _real_cupy():
    """The stand-in above must not outlive its tests."""
    assert getattr(sys.modules.get("cupy"), "__file__", None), "cupy is the numpy stand-in"


@requires_gpu
@pytest.mark.parametrize("path", list(GPU_PATHS))
def test_gpu_scalar_uniform_vector_and_each_groups_floor_on_the_device(path):
    _real_cupy()
    for pct in (0.1, 0.0):
        _ad_equal(_ad_run(GPU_PATHS[path], pct), _ad_run(GPU_PATHS[path], [pct] * 5))
    per_group = _ad_run(GPU_PATHS[path], FLOORS)
    for g, f in FLOORS.items():
        _ad_equal(per_group, _ad_run(GPU_PATHS[path], f), group=g)


@requires_gpu
def test_gpu_batch_path_applies_the_floor_on_the_device():
    _real_cupy()
    _gpu_batch_removes_what_the_cpu_batch_path_removes()


@requires_gpu
@requires_cytome
@pytest.mark.parametrize("batch_key", [None, "batch"])
def test_cytome_gpu_on_the_device(stores, batch_key):
    _real_cupy()
    run = lambda pct: cosg.run_cosg_cytome(stores["rna"], groupby="grp", n_genes_user=N_TOP, layer="counts",
                                           modality="RNA", verbose=False, device="gpu", batch_key=batch_key,
                                           expressed_pct=pct)
    _cy_equal(run(0.1), run([0.1] * 5))
    per_group = run(FLOORS)
    for g, f in FLOORS.items():
        assert _cy_group(per_group, g) == _cy_group(run(f), g)
