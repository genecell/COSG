"""The groups ``cosg.cosg`` reports, and what it refuses before computing.

- ``groups=`` with a categorical ``groupby`` kept every category: groups not
  asked for came back with arbitrary genes, and ``calculate_logfoldchanges``
  then raised IndexError. A category without cells (an unused one) did the
  same with ``groups='all'``. The result now holds the groups scored that have
  cells, as on the cytome paths, and the fold changes are those of the groups
  scored (rest = the other cells scored).
- ``calculate_logfoldchanges`` needed a categorical ``groupby``.
- A call that failed wiped the earlier result under the same key.
- A NaN or None label crashed with KeyError or TypeError; ``mu=NaN`` made every
  score NaN. Both are refused up front.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import cosg


def _adata(groups=("A", "B", "C", "D"), n=320, f=60, seed=0, cat=True, extra_categories=()):
    import anndata as ad
    rng = np.random.default_rng(seed)
    gi = rng.integers(0, len(groups), n)
    rate = np.full((len(groups), f), 0.1)
    for g in range(len(groups)):
        rate[g, g * 10:(g + 1) * 10] = 0.7
    X = np.log1p(((rng.random((n, f)) < rate[gi]) * rng.integers(1, 5, (n, f))).astype(np.float32))
    lab = np.array(groups)[gi]
    col = pd.Categorical(lab, categories=list(groups) + list(extra_categories)) if cat else lab
    a = ad.AnnData(sp.csr_matrix(X), obs=pd.DataFrame({"grp": col}, index=[f"c{i}" for i in range(n)]))
    a.var_names = [f"g{j}" for j in range(f)]
    return a


def _fields(uns):
    return {k: {g: np.asarray(uns[k][g]) for g in uns[k].dtype.names}
            for k in uns if k not in ("params", "COSG")}


def _equal(u1, u2):
    a, b = _fields(u1), _fields(u2)
    assert set(a) == set(b)
    for k in a:
        assert list(a[k]) == list(b[k]), k
        for g in a[k]:
            if a[k][g].dtype.kind == "f":
                np.testing.assert_allclose(a[k][g], b[k][g], rtol=1e-6, atol=1e-9, err_msg=f"{k} {g}")
            else:
                assert list(a[k][g]) == list(b[k][g]), (k, g)


@pytest.mark.parametrize("ccs", [None, 0])
def test_groups_subset_reports_the_groups_asked_for(ccs):
    a = _adata()
    cosg.cosg(a, groupby="grp", groups=["A", "C"], n_genes_user=6, cpu_chunk_size=ccs,
              calculate_logfoldchanges=True)
    assert list(a.uns["cosg"]["names"].dtype.names) == ["A", "C"]
    # the same as scoring the subset itself
    b = _adata()
    b = b[b.obs["grp"].isin(["A", "C"])].copy()
    b.obs["grp"] = b.obs["grp"].cat.remove_unused_categories()
    cosg.cosg(b, groupby="grp", n_genes_user=6, cpu_chunk_size=ccs, calculate_logfoldchanges=True)
    _equal(a.uns["cosg"], b.uns["cosg"])


def test_a_category_without_cells_is_not_scored():
    a = _adata(extra_categories=("Z",))
    cosg.cosg(a, groupby="grp", n_genes_user=6, calculate_logfoldchanges=True)
    assert list(a.uns["cosg"]["names"].dtype.names) == ["A", "B", "C", "D"]
    b = _adata()
    cosg.cosg(b, groupby="grp", n_genes_user=6, calculate_logfoldchanges=True)
    _equal(a.uns["cosg"], b.uns["cosg"])
    # a floor per group needs no fraction for it, and a sequence has one per group scored
    cosg.cosg(a, groupby="grp", n_genes_user=6, expressed_pct={"A": 0.1, "B": 0.1, "C": 0.1, "D": 0.1})
    cosg.cosg(a, groupby="grp", n_genes_user=6, expressed_pct=[0.1] * 4)


def test_fold_changes_without_a_categorical_groupby():
    a = _adata()
    b = _adata(cat=False)
    for x in (a, b):
        cosg.cosg(x, groupby="grp", n_genes_user=6, calculate_logfoldchanges=True)
    _equal(a.uns["cosg"], b.uns["cosg"])


def test_fold_changes_against_a_reference_group():
    a = _adata()
    cosg.cosg(a, groupby="grp", n_genes_user=6, calculate_logfoldchanges=True, reference="B")
    lfc = np.asarray(a.uns["cosg"]["logfoldchanges"]["B"])
    assert np.allclose(lfc, 0.0, atol=1e-6)                     # B against itself
    # a reference outside `groups` is compared with, not scored
    b = _adata()
    cosg.cosg(b, groupby="grp", groups=["A", "C"], n_genes_user=6, calculate_logfoldchanges=True,
              reference="B")
    assert list(b.uns["cosg"]["names"].dtype.names) == ["A", "C"]
    full = dict(zip(a.uns["cosg"]["names"]["A"], np.asarray(a.uns["cosg"]["logfoldchanges"]["A"])))
    for gene, v in zip(b.uns["cosg"]["names"]["A"], np.asarray(b.uns["cosg"]["logfoldchanges"]["A"])):
        if gene in full:
            assert v == pytest.approx(full[gene], abs=1e-5)           # A against B either way
    with pytest.raises(ValueError, match="reference"):
        cosg.cosg(a, groupby="grp", calculate_logfoldchanges=True, reference="nope")


@pytest.mark.parametrize("bad", [
    dict(pvalue_method="bogus", calculate_pvalues=True),
    dict(expressed_pct={"A": 0.1}),
    dict(mu=float("nan")),
    dict(mu=-1),
    dict(layer="nope"),
])
def test_a_failed_call_keeps_the_earlier_result(bad):
    a = _adata()
    cosg.cosg(a, groupby="grp", n_genes_user=3)
    before = {k: (v.copy() if hasattr(v, "copy") else v) for k, v in a.uns["cosg"].items()}
    with pytest.raises((ValueError, TypeError)):
        cosg.cosg(a, groupby="grp", n_genes_user=3, **bad)
    assert set(a.uns["cosg"]) == set(before)
    for k in ("names", "scores"):
        assert a.uns["cosg"][k].tobytes() == before[k].tobytes()


@pytest.mark.parametrize("cat", [True, False])
def test_cells_without_a_label_are_refused(cat):
    a = _adata(cat=cat)
    lab = a.obs["grp"].astype(object).to_numpy().copy()
    lab[:5] = None if not cat else np.nan
    a.obs["grp"] = pd.Categorical(lab) if cat else lab
    with pytest.raises(ValueError, match="5 of the 320 cells scored have no label"):
        cosg.cosg(a, groupby="grp")
    assert "cosg" not in a.uns
    # leaving them out with groups= is fine
    cosg.cosg(a, groupby="grp", groups=["A", "B", "C", "D"], n_genes_user=3)


# ── non-finite values, on every path ───────────────────────────────────────

def _with_nan(dense=False):
    a = _adata()
    X = a.X.toarray()
    X[0, 40] = np.nan
    X[3, 41] = np.inf
    a.X = X if dense else sp.csr_matrix(X)
    return a


@pytest.mark.parametrize("path", [dict(), dict(cpu_chunk_size=0), dict(batch_key="b"), dict(calculate_pvalues=True)])
@pytest.mark.parametrize("dense", [False, True])
def test_non_finite_values_are_refused(path, dense):
    a = _with_nan(dense)
    a.obs["b"] = np.where(np.arange(a.n_obs) % 2, "x", "y")
    with pytest.raises(ValueError, match=r"adata.X holds NaN or infinite values in 2 feature\(s\), e.g. \['g40', 'g41'\]"):
        cosg.cosg(a, groupby="grp", **path)
    assert "cosg" not in a.uns


def test_non_finite_values_in_raw_and_layers_are_named_by_their_source():
    a = _with_nan()
    a.layers["L"] = a.X.copy()
    with pytest.raises(ValueError, match=r"adata.layers\['L'\] holds NaN"):
        cosg.cosg(a, groupby="grp", layer="L")
    b = _adata()
    b.raw = _with_nan()
    with pytest.raises(ValueError, match="adata.raw.X holds NaN"):
        cosg.cosg(b, groupby="grp", use_raw=True)


@pytest.mark.parametrize("path", [dict(), dict(cpu_chunk_size=0), dict(batch_key="b"),
                                  dict(batch_key="b", cpu_chunk_size=0), dict(calculate_pvalues=True)])
def test_n_genes_user_above_the_number_of_genes_returns_every_gene(path):
    """argpartition raised (kth out of bounds) when n_genes_user exceeded the
    genes; the result now holds every gene of each group, ranked (the same on
    every store path: tests/test_path_parity.py)."""
    a = _adata(f=30)
    a.obs["b"] = np.where(np.arange(a.n_obs) % 2, "x", "y")
    cosg.cosg(a, groupby="grp", n_genes_user=50, **path)
    for g in a.uns["cosg"]["names"].dtype.names:
        names = list(a.uns["cosg"]["names"][g])
        assert len(names) == 30 and sorted(names) == sorted(a.var_names)
        scores = np.asarray(a.uns["cosg"]["scores"][g])
        assert np.all(np.diff(scores) <= 0)
