"""COSG's scores from per-group sums, and the scale of its normalized scores.

``cosg.score_from_sums`` is the last step of every streaming path, so the
same sums must give the scores those paths report, and the scores
``cosg.cosg`` reports on an AnnData. Sums of parts add up to the sums of the
whole, which is what makes the function useful outside COSG's own readers.

The scale ``iqr_normalize`` divides by keeps every group on a scale of its
own, also where the group scores above zero at fewer than 5 % of the
features and the 95th percentile is zero.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import cosg
from cosg._from_sums import _group_scale
from conftest import requires_cytome

GROUPS = np.array(list("abcde"))


def _data(seed=3, n=500, f=400, sparse_group=False):
    rng = np.random.default_rng(seed)
    gi = rng.choice(len(GROUPS), n, p=[0.3, 0.25, 0.2, 0.15, 0.1])
    rate = np.full((len(GROUPS), f), 0.05)
    for g in range(len(GROUPS)):
        rate[g, g * 70:g * 70 + 60] = 0.5
    if sparse_group:                       # e: open at 8 of 400 features only
        rate[4] = 0.0
        rate[4, 380:388] = 0.8
    X = (rng.random((n, f)) < rate[gi]) * rng.lognormal(0.5, 0.6, (n, f))
    return X.astype(np.float32), GROUPS[gi], np.where(np.arange(n) % 3 == 0, "B1", "B2")


def _sums(X, labels, groups=GROUPS):
    onehot = (labels[:, None] == groups[None, :]).astype(np.float64)
    Xd = X.astype(np.float64)
    return (Xd.T @ onehot, (Xd ** 2).sum(axis=0), ((Xd > 0).T @ onehot).astype(np.int64),
            onehot.sum(axis=0))


def _cytome(path, X, labels, batch, modality):
    import cytome
    n, f = X.shape
    ds = cytome.create(str(path))
    ds.set_entity("cells", pd.DataFrame({"cell_idx": np.arange(n), "barcode": [f"BC{i}" for i in range(n)],
                                         "grp": labels, "batch": batch}))
    if modality == "RNA":
        ds.set_entity("genes", pd.DataFrame({"gene_idx": np.arange(f), "gene_id": [f"F{i}" for i in range(f)]}))
    else:
        chroms = ["chr1"] * (f // 2) + ["chr2"] * (f - f // 2)
        starts = np.concatenate([np.arange(f // 2), np.arange(f - f // 2)]) * 1000
        ds.set_entity("peaks", pd.DataFrame({"peak_idx": np.arange(f), "chr": chroms, "start": starts,
                                             "end_": starts + 500,
                                             "peak_id": [f"F{i}" for i in range(f)]}))
    ds.add_matrix(f"{modality}_counts", sp.csr_matrix(X))
    ds.flush()
    ds.close()
    return str(path)


def _dense_raw(path, modality, f, **kw):
    """Every positive raw score run_cosg_cytome reports, as features x groups."""
    res = cosg.run_cosg_cytome(path, groupby="grp", modality=modality, layer="counts",
                               n_genes_user=f, output_format="dense", iqr_normalize=False,
                               verbose=False, **kw)
    df = res["scores_df"]
    return df.loc[[f"F{i}" for i in range(f)], list(GROUPS)].to_numpy()


@requires_cytome
@pytest.mark.parametrize("modality,batching", [("RNA", "none"), ("ATAC", "chromosome")])
@pytest.mark.parametrize("mu", [1.0, 100.0])
def test_the_streaming_paths_report_the_scores_of_their_sums(tmp_path, modality, batching, mu):
    X, labels, batch = _data()
    path = _cytome(tmp_path / f"{modality}.cytome", X, labels, batch, modality)
    got = _dense_raw(path, modality, X.shape[1], mu=mu, expressed_pct=0.1,
                     feature_batching=batching)
    dot, sq, nnz, sizes = _sums(X, labels)
    want = cosg.score_from_sums(dot, sq, sizes, nnz, mu=mu, expressed_pct=0.1)
    want = np.where(want > 0, want, 0.0)
    assert (got > 0).sum() > 100
    np.testing.assert_array_equal(got > 0, want > 0)
    np.testing.assert_allclose(got, want, rtol=2e-5, atol=1e-9)


@requires_cytome
def test_batch_sums_give_the_batch_corrected_scores(tmp_path):
    X, labels, batch = _data()
    path = _cytome(tmp_path / "b.cytome", X, labels, batch, "RNA")
    got = _dense_raw(path, "RNA", X.shape[1], mu=1.0, batch_key="batch")
    parts = [_sums(X[batch == b], labels[batch == b]) for b in ("B1", "B2")]
    dot = np.stack([p[0] for p in parts]); sq = np.stack([p[1] for p in parts])
    nnz = sum(p[2] for p in parts); bsizes = np.stack([p[3] for p in parts])
    want = cosg.score_from_sums(dot, sq, bsizes.sum(axis=0), nnz, batch_group_sizes=bsizes)
    want = np.where(want > 0, want, 0.0)
    np.testing.assert_array_equal(got > 0, want > 0)
    np.testing.assert_allclose(got, want, rtol=2e-5, atol=1e-9)


def test_anndata_cosg_reports_the_same_markers():
    ad = pytest.importorskip("anndata")
    X, labels, _ = _data()
    adata = ad.AnnData(sp.csr_matrix(X), obs=pd.DataFrame({"grp": pd.Categorical(labels)},
                                                          index=[f"c{i}" for i in range(len(labels))]))
    adata.var_names = [f"F{i}" for i in range(X.shape[1])]
    cosg.cosg(adata, groupby="grp", mu=1.0, n_genes_user=20, expressed_pct=0.1, key_added="cosg")
    names = pd.DataFrame(adata.uns["cosg"]["names"])
    scores = pd.DataFrame(adata.uns["cosg"]["scores"])
    dot, sq, nnz, sizes = _sums(X, labels)
    S = cosg.score_from_sums(dot, sq, sizes, nnz, mu=1.0, expressed_pct=0.1)
    for k, g in enumerate(GROUPS):
        top = np.argsort(-S[:, k], kind="stable")[:20]
        assert list(names[g]) == [f"F{i}" for i in top]
        np.testing.assert_allclose(scores[g].to_numpy(dtype=float), S[top, k], rtol=1e-5)


def test_sums_of_parts_add_up_to_the_whole():
    """Two halves of the cells, scored from their added sums, give the scores
    of all cells: the use of count-split halves, shards and added cells."""
    X, labels, _ = _data()
    h = np.arange(len(labels)) % 2 == 0
    a, b, whole = _sums(X[h], labels[h]), _sums(X[~h], labels[~h]), _sums(X, labels)
    added = cosg.score_from_sums(a[0] + b[0], a[1] + b[1], a[3] + b[3], a[2] + b[2])
    np.testing.assert_allclose(added, cosg.score_from_sums(whole[0], whole[1], whole[3], whole[2]),
                               rtol=1e-12, atol=1e-15)


def test_the_filters():
    X, labels, _ = _data()
    dot, sq, nnz, sizes = _sums(X, labels)
    kept = cosg.score_from_sums(dot, sq, sizes, nnz, remove_lowly_expressed=False)
    assert (kept >= 0).all()
    np.testing.assert_array_equal(kept, cosg.score_from_sums(dot, sq, sizes))
    filt = cosg.score_from_sums(dot, sq, sizes, nnz, expressed_pct=0.3)
    thr = np.maximum(sizes * 0.3, 3)
    np.testing.assert_array_equal(filt == -1, nnz < thr[None, :])
    rare = nnz.sum(axis=1) <= 40
    total = cosg.score_from_sums(dot, sq, sizes, nnz, remove_lowly_expressed=False, min_total_count=40)
    assert rare.any() and (total[rare] == -1).all() and (total[~rare] >= 0).all()


def test_shapes_are_checked():
    with pytest.raises(ValueError):
        cosg.score_from_sums(np.ones((4, 2)), np.ones(3), [2, 2])
    with pytest.raises(ValueError):
        cosg.score_from_sums(np.ones((2, 4, 2)), np.ones((2, 4)), [2, 2], batch_group_sizes=np.ones((3, 2)))


# ---- the scale of the normalized scores ----

def test_a_nonzero_spread_is_the_scale():
    col = np.linspace(0, 1, 101)
    assert _group_scale(col) == pytest.approx(0.95 - 0.75)


def test_ties_at_the_top_take_the_upper_quantile():
    col = np.r_[np.zeros(70), np.full(30, 0.4)]
    assert _group_scale(col) == pytest.approx(0.4)


def test_a_group_positive_at_few_features_keeps_a_scale_of_its_own():
    """q95 is zero when fewer than 5 % of the features score: the scale then
    is the 5th percentile of the positive scores, not the constant 1e-6."""
    pos = np.linspace(0.01, 0.2, 20)
    col = np.r_[np.zeros(980), pos]
    want = np.quantile(pos, 0.05)
    assert _group_scale(col) == pytest.approx(want)
    assert _group_scale(np.r_[np.zeros(990), np.full(10, 0.3)]) == pytest.approx(0.3)
    assert _group_scale(np.zeros(100)) == 1e-6
    assert _group_scale(np.array([])) == 1e-6


def test_the_scale_has_no_jump_where_the_upper_quantile_reaches_zero():
    """Just over 5 % positive, the spread over all features is a low quantile
    of the positive scores; just under, the fallback must land near it, not
    jump (QC doc 46.7: the spread of the positive scores, the earlier
    fallback, put the scales of two such groups up to 9.7 apart in log)."""
    rng = np.random.default_rng(3)
    scales = {}
    for share in (0.049, 0.051):
        n_pos = int(round(share * 100000))
        col = np.r_[np.zeros(100000 - n_pos), rng.lognormal(0.0, 1.0, n_pos)]
        scales[share] = _group_scale(col)
    assert 0.5 < scales[0.049] / scales[0.051] < 2.0, scales


def test_the_positive_median_rule():
    pos = np.linspace(0.01, 0.2, 20)
    for zeros in (80, 980, 9980):     # the rule does not move when zeros are added
        col = np.r_[np.zeros(zeros), pos]
        assert _group_scale(col, scale_rule="positive_median") == pytest.approx(np.median(pos))
    dense = np.linspace(0, 1, 101)    # dense groups too: the median of the positive scores
    assert _group_scale(dense, scale_rule="positive_median") == pytest.approx(np.median(dense[1:]))
    assert _group_scale(np.zeros(10), scale_rule="positive_median") == 1e-6
    assert _group_scale(np.array([]), scale_rule="positive_median") == 1e-6
    rng = np.random.default_rng(0)
    df = pd.DataFrame({"dense": rng.random(1000), "sparse": np.r_[np.zeros(980), pos]})
    out = cosg.iqrLogNormalize(df, scale_rule="positive_median")
    for c in df:
        v = df[c].to_numpy()
        np.testing.assert_allclose(out[c], np.log1p(v / np.median(v[v > 0])))
    with pytest.raises(ValueError):
        cosg.iqrLogNormalize(df, scale_rule="median")


def test_iqrLogNormalize_gives_a_zero_iqr_column_its_own_scale():
    rng = np.random.default_rng(0)
    dense = rng.random(1000)
    sparse_col = np.r_[np.zeros(980), np.linspace(0.01, 0.2, 20)]
    df = pd.DataFrame({"dense": dense, "sparse": sparse_col})
    out = cosg.iqrLogNormalize(df)
    iqr_dense = np.quantile(dense, 0.95) - np.quantile(dense, 0.75)
    np.testing.assert_allclose(out["dense"], np.log1p(dense / iqr_dense))
    np.testing.assert_allclose(out["sparse"], np.log1p(sparse_col / _group_scale(sparse_col)))


@requires_cytome
@pytest.mark.parametrize("modality,batching", [("RNA", "none"), ("ATAC", "chromosome")])
def test_a_sparse_group_is_normalized_on_its_own_scale(tmp_path, modality, batching):
    """Group e opens 8 of 400 features. Its normalized scores used to be set
    by the constant 1e-6 (about log(1e6) = 13.8 above its own scale); now
    they are on the scale of the other groups'."""
    X, labels, batch = _data(sparse_group=True)
    path = _cytome(tmp_path / f"s_{modality}.cytome", X, labels, batch, modality)
    res = cosg.run_cosg_cytome(path, groupby="grp", modality=modality, layer="counts", n_genes_user=50,
                               output_format="dict", iqr_normalize=True, feature_batching=batching,
                               expressed_pct=0.1, verbose=False)
    best = {}
    for (g, _), v in res["scores_dict"].items():
        best[g] = max(best.get(g, 0.0), v)
    others = max(v for g, v in best.items() if g != "e")
    assert best["e"] < 3 * others, best
    raw = _dense_raw(path, modality, X.shape[1], expressed_pct=0.1, feature_batching=batching)
    e = raw[:, 4]
    scale = _group_scale(e)
    top = sorted((v for (g, _), v in res["scores_dict"].items() if g == "e"), reverse=True)
    np.testing.assert_allclose(top, np.sort(np.log1p(e[e > 0] / scale))[::-1][:len(top)], rtol=1e-5)


def test_the_earlier_rules_are_kept_to_reproduce_earlier_scores():
    col = np.r_[np.zeros(980), np.linspace(0.01, 0.2, 20)]
    assert _group_scale(col, iqr_fallback="constant") == 1e-6
    assert _group_scale(np.r_[np.zeros(70), np.full(30, 0.4)], iqr_fallback="constant") == pytest.approx(0.4)
    rng = np.random.default_rng(0)
    df = pd.DataFrame({"dense": rng.random(1000), "sparse": col})
    old = df.quantile(0.95) - df.quantile(0.75)
    old = old.replace(0, old[old > 0].min())
    pd.testing.assert_frame_equal(cosg.iqrLogNormalize(df, iqr_fallback="smallest"), np.log1p(df / old))
    with pytest.raises(ValueError):
        cosg.iqrLogNormalize(df, iqr_fallback="constant")


@requires_cytome
@pytest.mark.parametrize("modality,batching", [("RNA", "none"), ("ATAC", "chromosome")])
def test_the_streaming_paths_take_the_positive_median_rule(tmp_path, modality, batching):
    X, labels, batch = _data(sparse_group=True)
    path = _cytome(tmp_path / f"m_{modality}.cytome", X, labels, batch, modality)
    kw = dict(groupby="grp", modality=modality, layer="counts", n_genes_user=50, output_format="dict",
              iqr_normalize=True, feature_batching=batching, expressed_pct=0.1, verbose=False)
    res = cosg.run_cosg_cytome(path, scale_rule="positive_median", **kw)
    raw = _dense_raw(path, modality, X.shape[1], expressed_pct=0.1, feature_batching=batching)
    groups = list(GROUPS)
    # a group with fewer than 20 positive scores takes the median of the other groups' own
    # medians (QC doc 52); the sparse group of these data has 8
    own = {k: np.median(raw[:, k][raw[:, k] > 0]) for k in range(len(groups))}
    n_pos = {k: int((raw[:, k] > 0).sum()) for k in range(len(groups))}
    pooled = np.median([own[k] for k in own if n_pos[k] >= 20])
    assert min(n_pos.values()) < 20 <= max(n_pos.values())
    for k, g in enumerate(groups):
        e = raw[:, k]
        top = sorted((v for (h, _), v in res["scores_dict"].items() if h == g), reverse=True)
        scale = own[k] if n_pos[k] >= 20 else pooled
        want = np.sort(np.log1p(e[e > 0] / scale))[::-1][:len(top)]
        np.testing.assert_allclose(top, want, rtol=1e-5)
    # the rankings within a group are the scale's to keep
    plain = cosg.run_cosg_cytome(path, **kw)
    for g in groups:
        a = [f for (h, f), _ in sorted(res["scores_dict"].items(), key=lambda x: -x[1]) if h == g]
        b = [f for (h, f), _ in sorted(plain["scores_dict"].items(), key=lambda x: -x[1]) if h == g]
        assert a == b
    via = cosg.cosg(path, scale_rule="positive_median", **{k: v for k, v in kw.items() if k != "verbose"})
    assert via["scores_dict"] == res["scores_dict"]
    with pytest.raises(ValueError):
        cosg.run_cosg_cytome(path, scale_rule="median", **kw)


@requires_cytome
def test_the_constant_fallback_reproduces_the_earlier_scale(tmp_path):
    X, labels, batch = _data(sparse_group=True)
    path = _cytome(tmp_path / "c.cytome", X, labels, batch, "ATAC")
    res = cosg.run_cosg_cytome(path, groupby="grp", modality="ATAC", layer="counts", n_genes_user=50,
                               output_format="dict", iqr_normalize=True, feature_batching="chromosome",
                               expressed_pct=0.1, verbose=False, iqr_fallback="constant")
    raw = _dense_raw(path, "ATAC", X.shape[1], expressed_pct=0.1, feature_batching="chromosome")
    e = raw[:, 4]
    top = sorted((v for (g, _), v in res["scores_dict"].items() if g == "e"), reverse=True)
    np.testing.assert_allclose(top, np.sort(np.log1p(e[e > 0] / 1e-6))[::-1][:len(top)], rtol=1e-5)
    with pytest.raises(ValueError):
        cosg.run_cosg_cytome(path, groupby="grp", modality="ATAC", iqr_fallback="smallest", verbose=False)
    # and through the unified entry point
    via = cosg.cosg(path, groupby="grp", modality="ATAC", layer="counts", n_genes_user=50,
                    output_format="dict", iqr_normalize=True, feature_batching="chromosome",
                    expressed_pct=0.1, iqr_fallback="constant")
    assert via["scores_dict"] == res["scores_dict"]


def test_the_scale_is_continuous_where_the_upper_quantile_reaches_zero():
    """QC doc 50.2: the spread fell to the smallest positive score just over 5 %
    positive and the fallback jumped back up just under it (0.044 against 0.198,
    log-normal); the spread is now floored at the fallback while the q_lower
    quantile is zero, so both sides meet, and groups further from 5 % keep the
    spread."""
    rng = np.random.default_rng(5)
    base = rng.lognormal(0.0, 1.0, 400000)
    scales = {}
    for share in (0.0499, 0.0500, 0.0501, 0.06):
        n_pos = int(round(share * 400000))
        col = np.r_[np.zeros(400000 - n_pos), base[:n_pos]]
        scales[share] = _group_scale(col)
    near = [scales[0.0499], scales[0.0500], scales[0.0501]]
    assert max(near) / min(near) < 1.1, scales
    n_pos = int(round(0.06 * 400000))
    col = np.r_[np.zeros(400000 - n_pos), base[:n_pos]]
    assert scales[0.06] == pytest.approx(np.quantile(col, 0.95))
    # the constant fallback of COSG 1.2.0 is left as it was
    n_pos = int(round(0.0501 * 400000))
    col = np.r_[np.zeros(400000 - n_pos), base[:n_pos]]
    assert _group_scale(col, iqr_fallback="constant") == pytest.approx(np.quantile(col, 0.95))


def test_the_median_rule_lends_a_group_with_few_positives_the_others_scale():
    """A group with fewer than 20 positive scores would be scaled by its own few
    markers (a lone marker normalized to log1p(1) however strong); it takes the
    median of the other groups' scales instead, on every path."""
    rng = np.random.default_rng(6)
    n = 5000
    df = pd.DataFrame({f"g{i}": np.r_[np.zeros(n - 1000), rng.lognormal(0.0, 0.5, 1000)] for i in range(3)})
    few = np.zeros(n)
    few[:3] = 50.0
    df["few"] = few
    out = cosg.iqrLogNormalize(df, scale_rule="positive_median")
    pooled = np.median([np.median(df[c][df[c] > 0]) for c in ("g0", "g1", "g2")])
    assert out["few"].max() == pytest.approx(np.log1p(50.0 / pooled))
    # a group with enough positives keeps its own median
    assert out["g0"].max() == pytest.approx(np.log1p(df["g0"].max() / np.median(df["g0"][df["g0"] > 0])))
    # without other groups (or none with enough positives) the group's own median stands
    assert _group_scale(few, scale_rule="positive_median") == pytest.approx(50.0)
