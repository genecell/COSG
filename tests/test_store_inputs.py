"""What ``run_cosg_cytome`` accepts, and what it refuses before it streams.

- ``cell_mask``: an integer mask with ``batch_key`` raised IndexError (it was
  cast to bool), and a boolean mask 10 short was accepted with the last 10
  cells silently left out. Both forms now go through one normalizer.
- A cell without a label (NULL) crashed with TypeError; a NaN or an infinity
  in the matrix went through unnoticed. Both are refused, the latter from the
  sums of squares the pass forms anyway, naming the features.
- ``mu < 0``, ``n_genes_user < 1`` and a minimum of 0 cells were accepted (the
  second failing later with a cryptic error), as the AnnData path does not.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import cosg
from conftest import requires_cytome

pytestmark = requires_cytome

N, F = 300, 60


def _data(seed=0):
    rng = np.random.default_rng(seed)
    gi = rng.integers(0, 3, N)
    rate = np.full((3, F), 0.1)
    for g in range(3):
        rate[g, g * 10:(g + 1) * 10] = 0.7
    X = ((rng.random((N, F)) < rate[gi]) * rng.integers(1, 5, (N, F))).astype(np.float32)
    return X, np.array(["A", "B", "C"])[gi], np.where(np.arange(N) % 2 == 0, "b1", "b2")


def _store(path, X, lab, batch, modality="RNA", column_major=False, labels=None):
    import cytome
    ds = cytome.create(str(path))
    ds.set_entity("cells", pd.DataFrame({"cell_idx": np.arange(N), "barcode": [f"c{i}" for i in range(N)],
                                         "grp": lab if labels is None else labels, "batch": batch}))
    if modality == "RNA":
        ds.set_entity("genes", pd.DataFrame({"gene_idx": np.arange(F), "gene_id": [f"g{j}" for j in range(F)]}))
    else:
        starts = np.arange(F) * 1000
        ds.set_entity("peaks", pd.DataFrame({"peak_idx": np.arange(F), "chr": ["chr1"] * 30 + ["chr2"] * 30,
                                             "start": starts, "end_": starts + 500,
                                             "peak_id": [f"g{j}" for j in range(F)]}))
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


@pytest.fixture(scope="module")
def store(tmp_path_factory):
    X, lab, batch = _data()
    return _store(tmp_path_factory.mktemp("inputs") / "s.cytome", X, lab, batch)


def _run(p, **kw):
    base = dict(groupby="grp", layer="counts", verbose=False, output_format="dict", n_genes_user=20)
    base.update(kw)
    return cosg.run_cosg_cytome(p, **base)


@pytest.mark.parametrize("batch_key", [None, "batch"])
@pytest.mark.parametrize("fb", ["none", "chromosome"])
def test_integer_and_boolean_masks_are_the_same_mask(store, batch_key, fb):
    keep = np.arange(N) % 3 != 0
    idx = np.flatnonzero(keep)
    rng = np.random.default_rng(1)
    kw = dict(batch_key=batch_key, feature_batching=fb)
    r_bool = _run(store, cell_mask=keep, **kw)
    r_int = _run(store, cell_mask=idx, **kw)
    r_shuffled = _run(store, cell_mask=rng.permutation(idx), **kw)
    assert r_bool["scores_dict"] == r_int["scores_dict"] == r_shuffled["scores_dict"]
    assert sum(r_bool["group_sizes"].values()) == keep.sum()


@pytest.mark.parametrize("mask,match", [
    (np.ones(N - 10, dtype=bool), "one per cell"),
    (np.ones(N + 1, dtype=bool), "one per cell"),
    (np.array([0, 1, 1, 2]), "more than once"),
    (np.array([0, N]), "outside"),
    (np.array([-1, 3]), "outside"),
    (np.zeros(N, dtype=bool), "selects no cell"),
    (np.array([], dtype=np.int64), "selects no cell"),
    (np.ones((N, 1), dtype=bool), "one-dimensional"),
])
def test_bad_masks_are_refused(store, mask, match):
    with pytest.raises(ValueError, match=match):
        _run(store, cell_mask=mask)


def test_a_float_mask_is_refused(store):
    with pytest.raises(TypeError, match="boolean mask"):
        _run(store, cell_mask=np.ones(N, dtype=float))


def test_unlabelled_cells_are_refused_unless_masked_out(tmp_path):
    X, lab, batch = _data()
    labels = lab.astype(object)
    labels[:7] = None
    p = _store(tmp_path / "u.cytome", X, lab, batch, labels=labels)
    with pytest.raises(ValueError, match="7 of the 300 cells scored have no label"):
        _run(p)
    keep = np.array([l is not None for l in labels])
    r = _run(p, cell_mask=keep)
    assert sum(r["group_sizes"].values()) == N - 7


@pytest.mark.parametrize("fb,column_major", [("none", False), ("chromosome", False), ("chromosome", True)])
def test_non_finite_values_are_refused_and_named(tmp_path, fb, column_major):
    X, lab, batch = _data()
    X[5, 40] = np.nan
    X[9, 7] = np.inf
    p = _store(tmp_path / "nan.cytome", X, lab, batch, modality="ATAC", column_major=column_major)
    with pytest.raises(ValueError, match=r"holds NaN or infinite values in .*'g(7|40)'"):
        cosg.run_cosg_cytome(p, groupby="grp", modality="ATAC", layer="counts", feature_batching=fb,
                             verbose=False)


@pytest.mark.parametrize("kw,match", [
    (dict(mu=-1.0), "mu"), (dict(mu=float("nan")), "mu"),
    (dict(n_genes_user=0), "n_genes_user"), (dict(n_genes_user=2.5), "n_genes_user"),
    (dict(expressed_min_num_cells_in_target_group=-1), "at least 0"),
])
def test_scoring_arguments_are_checked_as_on_anndata(store, kw, match):
    with pytest.raises(ValueError, match=match):
        _run(store, **kw)


def test_the_unbatched_sums_are_float64(tmp_path):
    """The unbatched path summed in float32: over a store of many chunks the
    sums lost digits chunk after chunk, about 1e-6 of a score here, against
    float32 rounding (6e-8) once the sums are float64."""
    import cytome
    rng = np.random.default_rng(5)
    n, f = 6000, 30
    lab = np.array(["A", "B", "C"])[rng.integers(0, 3, n)]
    X = np.log1p(rng.poisson(1.3, (n, f)) * rng.uniform(0.5, 1.5, (n, f))).astype(np.float32)
    p = str(tmp_path / "f64.cytome")
    ds = cytome.create(p)
    ds.set_entity("cells", pd.DataFrame({"cell_idx": np.arange(n), "barcode": [f"c{i}" for i in range(n)],
                                         "grp": lab}))
    ds.set_entity("genes", pd.DataFrame({"gene_idx": np.arange(f), "gene_id": [f"g{j}" for j in range(f)]}))
    ds.add_matrix("RNA_counts", sp.csr_matrix(X))
    ds.flush()
    ds.close()
    r = cosg.run_cosg_cytome(p, groupby="grp", layer="counts", n_genes_user=f, batch_size=8,
                             remove_lowly_expressed=False, verbose=False)
    onehot = (lab[:, None] == np.array(["A", "B", "C"])[None, :]).astype(np.float64)
    X64 = X.astype(np.float64)
    want = cosg.score_from_sums(X64.T @ onehot, (X64 ** 2).sum(axis=0), onehot.sum(axis=0),
                                remove_lowly_expressed=False)
    worst = 0.0
    for k, g in enumerate(r["groups_order"]):
        got = dict(zip(r["names"][:, k], r["scores"][:, k].astype(np.float64)))
        ref = want[:, ["A", "B", "C"].index(g)]
        for j in range(f):
            worst = max(worst, abs(got[f"g{j}"] - ref[j]) / max(abs(ref[j]), 1e-12))
    assert worst < 3e-7, worst


@pytest.mark.parametrize("fb", ["none", "chromosome"])
def test_a_feature_name_held_twice(tmp_path, fb):
    """A 10x GRCh38 reference holds 10 gene symbols twice (HSPA14, TBCE, ...).
    The dict, long and dense outputs are keyed by name: one of the two rows of
    the dense frame stayed 0 and the other took whichever feature's score came
    last. The dense frame is now filled by feature index, each row with its
    own score, and the dict keeps the better of the two."""
    import cytome
    X, lab, batch = _data()
    X[:, 41] = X[:, 3] * 0.5 + (np.random.default_rng(2).random(N) < 0.3)   # another, different feature
    names = [f"g{j}" for j in range(F)]
    names[41] = "g3"                                                        # ... named like g3
    p = str(tmp_path / "dup.cytome")
    ds = cytome.create(p)
    ds.set_entity("cells", pd.DataFrame({"cell_idx": np.arange(N), "barcode": [f"c{i}" for i in range(N)],
                                         "grp": lab}))
    starts = np.arange(F) * 1000
    ds.set_entity("peaks", pd.DataFrame({"peak_idx": np.arange(F), "chr": ["chr1"] * F, "start": starts,
                                         "end_": starts + 500, "peak_id": names}))
    ds.flush()
    w = ds.create_layer_writer("ATAC_counts", N, F, dtype=np.float32, compression="zlib")
    w.write_chunk(sp.csr_matrix(X.astype(np.float32)), 0)
    w.finalize()
    ds.flush()
    ds.close()
    kw = dict(groupby="grp", modality="ATAC", layer="counts", verbose=False, n_genes_user=F,
              feature_batching=fb, iqr_normalize=False)
    nd = cosg.run_cosg_cytome(p, **kw)
    dense = cosg.run_cosg_cytome(p, output_format="dense", **kw)["scores_df"]
    d = cosg.run_cosg_cytome(p, output_format="dict", **kw)["scores_dict"]
    k = list(nd["groups_order"]).index("A")
    # the two features' own scores, from the feature-indexed ndarray output
    by_name = {}
    for nm, sc in zip(nd["names"][:, k], nd["scores"][:, k]):
        by_name.setdefault(nm, []).append(float(sc))
    s3, s41 = by_name["g3"]                                   # best first
    assert s3 != s41 and s3 > 0 and s41 > 0
    rows = np.flatnonzero(dense.index == "g3")
    assert list(rows) == [3, 41]
    got = sorted(dense["A"].to_numpy()[rows].tolist(), reverse=True)
    assert got == pytest.approx([s3, s41], rel=1e-6)          # each row its own score
    assert d[("A", "g3")] == pytest.approx(s3, rel=1e-6)      # the dict keeps the better one


def test_a_feature_name_held_twice_warns_where_it_cannot_be_kept(tmp_path):
    """dict and long are keyed by name: they keep the better score and say so;
    dense and ndarray keep both features, by position, and say nothing."""
    import warnings
    import cytome
    X, lab, _ = _data()
    X[:, 41] = X[:, 3]
    names = [f"g{j}" for j in range(F)]
    names[41] = "g3"
    p = str(tmp_path / "dup2.cytome")
    ds = cytome.create(p)
    ds.set_entity("cells", pd.DataFrame({"cell_idx": np.arange(N), "barcode": [f"c{i}" for i in range(N)],
                                         "grp": lab}))
    ds.set_entity("genes", pd.DataFrame({"gene_idx": np.arange(F), "gene_id": [f"ENSG{j}" for j in range(F)],
                                         "symbol": names}))
    ds.add_matrix("RNA_counts", sp.csr_matrix(X))
    ds.flush()
    ds.close()
    kw = dict(groupby="grp", layer="counts", verbose=False, n_genes_user=F, feature_name_col="symbol")
    for fmt in ("dict", "long"):
        with pytest.warns(UserWarning, match=r"1 feature name\(s\) belong to more than one scored feature, e.g. \['g3'\]"):
            cosg.run_cosg_cytome(p, output_format=fmt, **kw)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        dense = cosg.run_cosg_cytome(p, output_format="dense", **kw)["scores_df"]
        cosg.run_cosg_cytome(p, **kw)
    rows = np.flatnonzero(dense.index == "g3")
    np.testing.assert_array_equal(dense.iloc[rows[0]].to_numpy(), dense.iloc[rows[1]].to_numpy())
    assert (dense.iloc[rows[0]].to_numpy() > 0).any()               # neither row is all zero
