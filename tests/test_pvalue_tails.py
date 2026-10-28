"""The tail of the p-value: refined where the normal approximation is wrong,
reported in log space on every path, and never an exact 0.

- The refinement used to stop at a normal p of 1e-30, on the view that the
  saddlepoint is noise that deep. But there it is the normal p that is wrong,
  in the anti-conservative direction: a gene in 12 cells, all of a 100-cell
  group of 5,000, had a normal log10 p of -129.7 against an exact -20.7, and
  dense genes deep in the tail were tens of orders too small. The saddlepoint
  on atoms solves in log space and follows the exact tail there.
- The cytome path reported p = 0.0 for strong markers, had no ``neg_log10``
  columns, and threw away the saddlepoint's log p.
- ``pvalue_method='exact'`` truncates its table at ``h_max`` expressing cells
  of the group; a strong marker's tail lies past it and came back as 0.
- ``compute_pvalues(method='permutation')`` passed an argument that does not
  exist, and dropped ``batch``.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import cosg
from cosg._pvalue_atoms import exact_tail_atoms, exact_tail_bounded
from cosg._pvalues import analytic_pvalues, compute_pvalues, permutation_pvalues
from conftest import requires_cytome

N, NC = 5000, 100


def _exact_log10(x, lab, group="A"):
    v, w = np.unique(x[x > 0], return_counts=True)
    inside = lab == group
    t = x[inside].sum()
    k = int((x > 0).sum())
    p = exact_tail_atoms(v, w.astype(float), len(x), int(inside.sum()), t, int(min(k, inside.sum())))
    return np.log10(p)


def _data(seed=0):
    """Two sparse markers of A (the review's), and three denser genes deep in
    the tail, against a 100-cell group of 5,000."""
    rng = np.random.default_rng(seed)
    lab = np.where(np.arange(N) < NC, "A", "B")
    cols = []
    x = np.zeros(N); x[:12] = 1.0; cols.append(x)                      # 12 cells, all in A
    x = np.zeros(N); x[:30] = 1.0; x[NC:NC + 30] = 1.0; cols.append(x)  # 30 of its 60 in A
    for k, shift in ((1500, 2), (2000, 1), (4500, 3)):
        x = np.zeros(N)
        x[:NC] = 1 + shift + rng.poisson(1.0, NC)
        rest = NC + rng.choice(N - NC, k - NC, replace=False)
        x[rest] = 1 + rng.poisson(1.0, k - NC)
        cols.append(x)
    cols.append(rng.poisson(0.5, N).astype(float))                      # a null gene
    return np.column_stack(cols), lab


def test_the_saddlepoint_follows_the_exact_tail_where_the_normal_does_not():
    """Sparse markers are refined however deep their normal p (it used to stop
    at 1e-30, which left exactly them at the normal value); dense ones below
    1e-30 keep the normal p unless spa_min_p=0 asks for every entry."""
    X, lab = _data()
    P, Z, info = analytic_pvalues(sp.csr_matrix(X), lab, ["A", "B"], method="spa")
    got = info["log_p"][:, 0] / np.log(10)
    for j in range(2):                              # the two sparse markers
        want = _exact_log10(X[:, j], lab)
        assert abs(got[j] - want) < 1.0, (j, got[j], want)
    assert got[0] > -25 and got[1] > -40            # not -129.7 and -156.8
    _, Z0, normal = analytic_pvalues(sp.csr_matrix(X), lab, ["A", "B"], method="normal")
    for j in range(2, 5):                           # dense, deep: the normal value
        assert got[j] == pytest.approx(normal["log_p"][j, 0] / np.log(10), abs=1e-9)
    _, _, every = analytic_pvalues(sp.csr_matrix(X), lab, ["A", "B"], method="spa", spa_min_p=0.0)
    for j in range(5):
        want = _exact_log10(X[:, j], lab)
        assert abs(every["log_p"][j, 0] / np.log(10) - want) < 1.0, (j, want)


def test_no_p_below_the_permutation_floor():
    """At the edge of the null's support (the group holds a gene's top values)
    the saddlepoint cannot converge, and the normal p it kept could be
    hundreds of orders below the smallest p the permutation null can give,
    one labeling in prod_b C(n_b, n_c,b). That floor is the exact p there."""
    import warnings
    from scipy.special import gammaln
    rng = np.random.default_rng(0)
    n, nc = 20000, 2000
    lab = np.where(np.arange(n) < nc, "A", "B")
    bat = np.where(np.arange(n) % 2 == 0, "b1", "b2")
    x0 = np.log1p(rng.poisson(2, n).astype(float))
    x0[:nc] = np.log1p(1.7 * rng.poisson(40, nc))           # A holds the top values
    x1 = np.zeros(n)
    x1[:nc] = np.log1p(rng.uniform(1, 9, nc))                # expressed in A only
    X = sp.csr_matrix(np.column_stack([x0, x1]))
    lc = lambda N, K: gammaln(N + 1) - gammaln(K + 1) - gammaln(N - K + 1)
    for batch, floor in ((None, -lc(n, nc)), (bat, -2 * lc(n // 2, nc // 2))):
        for min_p in (0.0, 1e-30):
            with warnings.catch_warnings():
                warnings.simplefilter("error")               # no "did not converge" warning
                P, Z, info = analytic_pvalues(X, lab, ["A", "B"], batch=batch, method="spa",
                                              spa_min_p=min_p)
            lp = info["log_p"][:, 0]
            assert (lp >= floor - 1e-9).all(), (lp, floor)   # never below the floor
            # the exact p of both is the floor (the group is the one subset
            # reaching the maximum); the saddlepoint, where it converges, is
            # above it by a few orders, and where it fails the floor stands
            assert (lp <= floor + 25.0).all(), (lp, floor)
    _, _, normal = analytic_pvalues(X, lab, ["A", "B"], batch=bat, method="normal")
    assert (normal["log_p"][:, 0] >= -2 * lc(n // 2, nc // 2) - 1e-9).all()   # the normal p is floored too


@requires_cytome
def test_every_path_reports_log_p_and_no_exact_zero(tmp_path):
    import anndata as ad
    import cytome
    X, lab = _data()
    perfect = np.where(lab == "A", 50.0, 0.0)        # a perfect marker: its p underflows
    X = np.column_stack([X, perfect])
    names = [f"g{j}" for j in range(X.shape[1])]
    a = ad.AnnData(sp.csr_matrix(X.astype(np.float32)),
                   obs=pd.DataFrame({"grp": pd.Categorical(lab)}, index=[f"c{i}" for i in range(N)]))
    a.var_names = names
    cosg.cosg(a, groupby="grp", n_genes_user=X.shape[1], calculate_pvalues=True)
    ref = {f: dict(zip(a.uns["cosg"]["names"]["A"], np.asarray(a.uns["cosg"][f]["A"])))
           for f in ("pvals", "pvals_adj", "neg_log10_pvals", "neg_log10_pvals_adj")}

    p = str(tmp_path / "pv.cytome")
    ds = cytome.create(p)
    ds.set_entity("cells", pd.DataFrame({"cell_idx": np.arange(N), "barcode": [f"c{i}" for i in range(N)],
                                         "grp": lab}))
    ds.set_entity("genes", pd.DataFrame({"gene_idx": np.arange(X.shape[1]), "gene_id": names}))
    ds.add_matrix("RNA_counts", sp.csr_matrix(X.astype(np.float32)))
    ds.flush(); ds.close()
    nd = cosg.run_cosg_cytome(p, groupby="grp", layer="counts", n_genes_user=X.shape[1],
                              calculate_pvalues=True, verbose=False)
    k = list(nd["groups_order"]).index("A")
    for f in ref:
        got = dict(zip(nd["names"][:, k], nd[f][:, k]))
        for g in names:
            np.testing.assert_allclose(got[g], ref[f][g], rtol=1e-6, err_msg=f"{f} {g}")
    assert (nd["pvals"] > 0).all()
    # the perfect marker: the group draws exactly the expressing cells, one
    # subset of C(5000, 100), about 1e-211 (a float64 p that the cytome path
    # used to report as 0)
    from scipy.special import gammaln
    log10_subsets = (gammaln(N + 1) - gammaln(NC + 1) - gammaln(N - NC + 1)) / np.log(10)
    got = nd["neg_log10_pvals"][list(nd["names"][:, k]).index("g6"), k]
    assert abs(got - log10_subsets) < 1.0, (got, log10_subsets)
    d = cosg.run_cosg_cytome(p, groupby="grp", layer="counts", n_genes_user=4, calculate_pvalues=True,
                             output_format="dict", verbose=False)
    assert set(d["neg_log10_pvals"]) == set(d["pvals"]) == set(d["scores_dict"])
    for (g, f), v in d["neg_log10_pvals"].items():
        if g == "A":
            assert v == pytest.approx(ref["neg_log10_pvals"][f], rel=1e-6)
    assert all(v > 0 for v in d["pvals"].values())


def test_exact_never_reports_the_tail_its_table_cut_off():
    X, lab = _data()
    X[:NC, 6 - 1] = 20.0                             # replace the null gene: a strong integer marker
    X[NC:, 6 - 1] = (np.random.default_rng(1).random(N - NC) < 0.2) * 1.0
    Xs = sp.csr_matrix(X)
    P_ex, _, info = analytic_pvalues(Xs, lab, ["A", "B"], method="exact")
    P_spa, _, info_spa = analytic_pvalues(Xs, lab, ["A", "B"], method="spa")
    assert (P_ex > 0).all()
    lex, lspa = info["log_p"][:, 0], info_spa["log_p"][:, 0]
    assert np.abs(lex - lspa).max() / np.log(10) < 1.0


def test_the_bounded_exact_tail():
    rng = np.random.default_rng(0)
    n, n_c = 2000, 100
    x = np.zeros(n); x[:60] = 1 + rng.integers(0, 3, 60); x[n_c:n_c + 40] = 1.0
    v, w = np.unique(x[x > 0], return_counts=True)
    t = x[:n_c].sum()
    full = exact_tail_atoms(v, w.astype(float), n, n_c, t, 100)
    assert exact_tail_bounded(v, w.astype(float), n, n_c, t, 100) == full      # no cut: exact
    assert exact_tail_bounded(v, w.astype(float), n, n_c, t, 30) is None      # the tail lies past 30
    # a null-sized statistic: the cut leaves out a negligible mass
    t_mid = x.sum() * n_c / n
    p_mid = exact_tail_bounded(v, w.astype(float), n, n_c, t_mid, 40)
    assert p_mid is not None and p_mid == pytest.approx(exact_tail_atoms(v, w.astype(float), n, n_c, t_mid, 100),
                                                        rel=1e-3)


@pytest.mark.parametrize("batch", [None, "two"])
def test_compute_pvalues_permutation(batch):
    rng = np.random.default_rng(3)
    n = 300
    lab = rng.choice(["A", "B"], n)
    bat = None if batch is None else np.where(np.arange(n) < 150, "b1", "b2")
    X = sp.csr_matrix(rng.poisson(1.0, (n, 10)).astype(float))
    out = compute_pvalues(X, lab, ["A", "B"], batch=bat, method="permutation", n_permutations=50)
    ref = permutation_pvalues(X, lab, ["A", "B"], n_permutations=50, batch=bat)
    np.testing.assert_array_equal(out["p_value"], ref)
    np.testing.assert_array_equal(out["info"]["log_p"], np.log(ref))      # its own log p
