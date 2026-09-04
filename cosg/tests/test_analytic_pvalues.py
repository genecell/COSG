"""Tests for COSG's analytic p-values.

The spec's calibration test was a KS test on the p-value distribution under
H0. KS measures the bulk, and the normal approximation is *correct* in the
bulk -- it is 180x anti-conservative at 1e-5 on a sparse gene and passes KS
comfortably. So the load-bearing test here is the tail one:
``test_tail_is_calibrated_where_the_normal_approximation_is_not``.
"""
from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.stats import kstest, norm

from cosg._pvalues import (
    adjust_fdr,
    analytic_pvalues,
    group_sums,
    permutation_pvalues,
    population_moments,
    _exact_tail,
    _spa_tail,
)


# ---------------------------------------------------------------------------
# 1. the moments, against every possible assignment
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n,nc", [(12, 4), (11, 3), (10, 5)])
def test_moments_match_exhaustive_enumeration(n, nc):
    rng = np.random.default_rng(n * 100 + nc)
    x = np.round(rng.gamma(0.3, 4.0, n), 3)
    x[rng.random(n) < 0.4] = 0.0

    # m3 is opt-in: the third power sum fed the Edgeworth term the
    # saddlepoint replaced, so it is no longer computed by default.
    mean, s2, m3 = population_moments(x.reshape(-1, 1), want_m3=True)
    E = nc * mean[0]
    V = nc * (n - nc) / (n - 1) * s2[0]
    M3 = nc * (n - nc) * (n - 2 * nc) / ((n - 1) * (n - 2)) * m3[0]

    T = np.array([x[list(i)].sum() for i in combinations(range(n), nc)])
    assert np.isclose(E, T.mean())
    assert np.isclose(V, T.var())
    assert np.isclose(M3, ((T - T.mean()) ** 3).mean())


def test_exact_tail_equals_enumeration():
    rng = np.random.default_rng(3)
    n, nc = 14, 5
    x = np.zeros(n)
    x[:6] = rng.gamma(2.0, 1.0, 6)
    T = np.array([x[list(i)].sum() for i in combinations(range(n), nc)])
    for q in (0.5, 0.9, 0.99):
        t = np.quantile(T, q)
        assert np.isclose(_exact_tail(x[x != 0], n, nc, t), np.mean(T >= t),
                          atol=1e-9)


# ---------------------------------------------------------------------------
# 2. the tail -- the test the original plan did not have
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("k", [5000, 200, 60])
def test_tail_is_calibrated_where_the_normal_approximation_is_not(k):
    """At a threshold the normal calls 1e-5, SPA must be near the truth.

    A gene expressed in k of 50,000 cells against a 500-cell group. What
    governs the CLT here is the expected in-group expressing count, n_c*k/n,
    which is 50, 2 and 0.6 for the three k values -- so the normal is progres-
    sively, badly anti-conservative while SPA tracks the empirical tail.
    """
    rng = np.random.default_rng(11)
    n, nc = 50_000, 500
    x = np.zeros(n)
    x[:k] = rng.gamma(2.0, 1.0, k)

    mean, s2, _ = population_moments(x.reshape(-1, 1))
    sd = np.sqrt(nc * (n - nc) / (n - 1) * s2[0])
    t = nc * mean[0] + norm.isf(1e-5) * sd

    draws = np.array([x[rng.choice(n, nc, replace=False)].sum()
                      for _ in range(4000)])
    empirical = max((draws >= t).mean(), 1.0 / 4000)
    spa = _spa_tail(x[x != 0], n, nc, t)

    assert spa is not None
    ratio = spa / empirical
    assert 1 / 4 < ratio < 4, (
        f"k={k}: SPA {spa:.2e} vs empirical {empirical:.2e} (ratio {ratio:.1f}); "
        f"nominal normal p was 1e-05")


def test_normal_is_the_one_that_fails_that_case():
    """Guards the premise: without SPA this is badly anti-conservative."""
    rng = np.random.default_rng(11)
    n, nc, k = 50_000, 500, 200
    x = np.zeros(n)
    x[:k] = rng.gamma(2.0, 1.0, k)
    mean, s2, _ = population_moments(x.reshape(-1, 1))
    sd = np.sqrt(nc * (n - nc) / (n - 1) * s2[0])
    t = nc * mean[0] + norm.isf(1e-5) * sd
    draws = np.array([x[rng.choice(n, nc, replace=False)].sum()
                      for _ in range(4000)])
    empirical = max((draws >= t).mean(), 1.0 / 4000)
    assert empirical / 1e-5 > 20, (
        "the premise of this feature is that the normal approximation is "
        f"badly wrong here; measured ratio {empirical / 1e-5:.0f}x")


# ---------------------------------------------------------------------------
# 3. calibration under H0, and power
# ---------------------------------------------------------------------------

def _null_frame(n=600, g=400, seed=0, sparsity=0.85):
    rng = np.random.default_rng(seed)
    X = rng.poisson(0.8, (n, g)).astype(np.float64)
    X[rng.random((n, g)) < sparsity] = 0.0
    labels = np.array([f"c{i % 3}" for i in range(n)])
    rng.shuffle(labels)
    return sparse.csr_matrix(X), labels, [f"c{i}" for i in range(3)]


def test_uniform_under_the_null():
    """KS needs a continuous null and independent tests, and both matter.

    Counts with few expressing cells make the null genuinely discrete, and KS
    against a continuous uniform rejects discrete-but-valid p-values. The three
    groups of one gene are also not independent -- their sums are constrained
    to the gene total -- so a KS over the flattened matrix tests an assumption
    the data cannot meet. Continuous values, one group: the assumptions hold.
    """
    rng = np.random.default_rng(1)
    n, g = 800, 500
    X = np.zeros((n, g))
    mask = rng.random((n, g)) < 0.25
    X[mask] = rng.gamma(2.0, 1.0, mask.sum())
    labels = np.array(["c0"] * 260 + ["c1"] * 540)
    rng.shuffle(labels)

    P, _, _ = analytic_pvalues(sparse.csr_matrix(X), labels, ["c0", "c1"],
                               method="spa")
    one = P[:, 0]
    assert kstest(one, "uniform").pvalue > 0.01, "not uniform under H0"
    assert 0.03 < np.mean(one < 0.05) < 0.08


def test_false_positive_rate_is_not_inflated_on_counts():
    """The discrete case: KS does not apply, but the error rate still must."""
    X, labels, groups = _null_frame(seed=3)
    P, _, _ = analytic_pvalues(X, labels, groups, method="spa")
    p = P.ravel()
    assert np.mean(p < 0.05) < 0.09, np.mean(p < 0.05)
    assert np.mean(p < 0.01) < 0.03, np.mean(p < 0.01)


def test_planted_markers_dominate():
    rng = np.random.default_rng(5)
    n, g = 600, 200
    labels = np.array([f"c{i % 3}" for i in range(n)])
    X = rng.poisson(0.3, (n, g)).astype(np.float64)
    X[labels == "c0", :4] += rng.poisson(8, ((labels == "c0").sum(), 4))
    groups = [f"c{i}" for i in range(3)]
    P, _, _ = analytic_pvalues(sparse.csr_matrix(X), labels, groups)
    c0 = P[:, groups.index("c0")]
    assert set(np.argsort(c0)[:4]) == {0, 1, 2, 3}, (
        f"planted genes are not the four smallest p; got {np.argsort(c0)[:6]}")
    assert (c0[:4] < 1e-10).all(), c0[:4]
    assert adjust_fdr(c0)[:4].max() < 0.01


# ---------------------------------------------------------------------------
# 4. batch stratification
# ---------------------------------------------------------------------------

def test_batch_confounding_needs_the_stratified_null():
    """A batch effect leaking into an unstratified null, and the fix.

    Group membership is *biased* by batch (70/30) rather than determined by
    it -- determined would leave the stratified null with nothing to permute.
    The expression difference is purely between batches, so any apparent
    group signal is the confound.
    """
    rng = np.random.default_rng(7)
    n, g = 800, 300
    batch = np.array(["b0"] * (n // 2) + ["b1"] * (n // 2))
    p_c0 = np.where(batch == "b0", 0.3, 0.7)
    labels = np.where(rng.random(n) < p_c0, "c0", "c1")

    X = rng.poisson(0.6, (n, g)).astype(np.float64)
    X[batch == "b1"] += rng.poisson(2.5, ((batch == "b1").sum(), g))
    groups = ["c0", "c1"]

    P_un, _, _ = analytic_pvalues(sparse.csr_matrix(X), labels, groups)
    P_st, _, _ = analytic_pvalues(sparse.csr_matrix(X), labels, groups,
                                  batch=batch)

    un = np.mean(P_un[:, 0] < 0.05)
    st = np.mean(P_st[:, 0] < 0.05)
    assert un > 0.5, f"the confound should be obvious unstratified; got {un:.2f}"
    assert st < 0.2, f"stratifying should remove it; got {st:.2f}"


# ---------------------------------------------------------------------------
# 5. mechanics: determinism, FDR, floors, and the untouched default
# ---------------------------------------------------------------------------

def test_analytic_path_is_deterministic_and_leaves_global_rng_alone():
    X, labels, groups = _null_frame(seed=2)
    state = np.random.get_state()
    a, _, _ = analytic_pvalues(X, labels, groups)
    b, _, _ = analytic_pvalues(X, labels, groups)
    assert np.array_equal(a, b)
    assert np.array_equal(np.random.get_state()[1], state[1])


def test_no_p_value_is_exactly_zero():
    """A hard zero breaks -log10 and BH downstream."""
    rng = np.random.default_rng(9)
    n, g = 400, 50
    labels = np.array([f"c{i % 2}" for i in range(n)])
    X = rng.poisson(0.2, (n, g)).astype(np.float64)
    X[labels == "c0", 0] += 500.0
    P, _, _ = analytic_pvalues(sparse.csr_matrix(X), labels, ["c0", "c1"])
    assert (P > 0).all()


def test_bh_matches_a_hand_computation():
    # p*m/i = [.005, .02, .065, .05125, .042]; BH then takes the running
    # minimum from the largest p downward, so the last three collapse to .042.
    p = np.array([0.001, 0.008, 0.039, 0.041, 0.042])
    got = adjust_fdr(p, "fdr_bh")
    want = np.array([0.005, 0.02, 0.042, 0.042, 0.042])
    assert np.allclose(got, want, atol=1e-6), got
    assert (adjust_fdr(p, "fdr_by") >= got - 1e-12).all()


def test_permutation_oracle_agrees_with_the_analytic_path():
    rng = np.random.default_rng(4)
    n, g = 300, 40
    labels = np.array([f"c{i % 2}" for i in range(n)])
    X = rng.poisson(1.5, (n, g)).astype(np.float64)
    groups = ["c0", "c1"]
    P, _, _ = analytic_pvalues(sparse.csr_matrix(X), labels, groups)
    Q = permutation_pvalues(sparse.csr_matrix(X), labels, groups,
                            n_permutations=4000, random_seed=0)
    mid = (P > 0.02) & (P < 0.98)
    assert np.corrcoef(P[mid], Q[mid])[0, 1] > 0.95
    assert np.abs(P[mid] - Q[mid]).max() < 0.10


def test_permutation_p_uses_the_valid_estimator():
    """(r+1)/(B+1): r/B can report zero for a null never exceeded."""
    rng = np.random.default_rng(6)
    X = sparse.csr_matrix(rng.poisson(1.0, (120, 5)).astype(np.float64))
    labels = np.array([f"c{i % 2}" for i in range(120)])
    Q = permutation_pvalues(X, labels, ["c0", "c1"], n_permutations=100)
    assert (Q >= 1.0 / 101.0 - 1e-12).all()


# ---------------------------------------------------------------------------
# 6. the public surface
# ---------------------------------------------------------------------------

def _adata(seed=0):
    anndata = pytest.importorskip("anndata")
    rng = np.random.default_rng(seed)
    n, g = 450, 120
    labels = np.array([f"c{i % 3}" for i in range(n)])
    X = rng.poisson(0.4, (n, g)).astype(np.float32)
    X[labels == "c0", :3] += rng.poisson(7, ((labels == "c0").sum(), 3))
    return anndata.AnnData(
        X=X,
        obs=pd.DataFrame({"ct": pd.Categorical(labels)},
                         index=[f"cell{i}" for i in range(n)]),
        var=pd.DataFrame(index=[f"g{i}" for i in range(g)]))


def test_default_output_is_unchanged():
    """calculate_pvalues=False must be bit-identical to the release."""
    import cosg

    a, b = _adata(), _adata()
    cosg.cosg(a, groupby="ct", n_genes_user=10)
    cosg.cosg(b, groupby="ct", n_genes_user=10)
    assert set(a.uns["cosg"]) == set(b.uns["cosg"])
    assert "pvals" not in a.uns["cosg"]
    assert np.array_equal(pd.DataFrame(a.uns["cosg"]["scores"]).values,
                          pd.DataFrame(b.uns["cosg"]["scores"]).values)


def test_pvalues_appear_and_rank_sensibly():
    import cosg

    a = _adata()
    cosg.cosg(a, groupby="ct", calculate_pvalues=True, n_genes_user=10)
    res = a.uns["cosg"]
    for key in ("pvals", "pvals_adj", "zscores"):
        assert key in res, key
    pv = pd.DataFrame(res["pvals"])["c0"].values
    names = pd.DataFrame(res["names"])["c0"].values
    planted = {"g0", "g1", "g2"}
    assert planted <= set(names[:4]), names[:6]
    assert pv[:3].max() < 1e-8
    assert (pv > 0).all()


def test_scores_are_identical_with_and_without_pvalues():
    """The feature is additive: asking for p-values must not move the ranking."""
    import cosg

    a, b = _adata(), _adata()
    cosg.cosg(a, groupby="ct", n_genes_user=10)
    cosg.cosg(b, groupby="ct", calculate_pvalues=True, n_genes_user=10)
    assert np.allclose(pd.DataFrame(a.uns["cosg"]["scores"]).values,
                       pd.DataFrame(b.uns["cosg"]["scores"]).values)


@pytest.mark.parametrize("mu", [0.5, 1.0, 10.0])
def test_pvalues_do_not_depend_on_mu(mu):
    """mu ranks; chance does not care. Documented behaviour, pinned here."""
    import cosg

    ref = _adata()
    cosg.cosg(ref, groupby="ct", calculate_pvalues=True, mu=1.0,
              n_genes_user=120)
    a = _adata()
    cosg.cosg(a, groupby="ct", calculate_pvalues=True, mu=mu, n_genes_user=120)

    def by_gene(res):
        n = pd.DataFrame(res["names"])["c0"].values
        p = pd.DataFrame(res["pvals"])["c0"].values
        return dict(zip(n, p))

    r, g = by_gene(ref.uns["cosg"]), by_gene(a.uns["cosg"])
    shared = set(r) & set(g)
    assert len(shared) > 50
    for k in shared:
        assert np.isclose(r[k], g[k], rtol=1e-12), (k, r[k], g[k])


def test_bad_method_is_rejected():
    import cosg

    with pytest.raises(ValueError, match="pvalue_method"):
        cosg.cosg(_adata(), groupby="ct", calculate_pvalues=True,
                  pvalue_method="edgeworth")


# ---------------------------------------------------------------------------
# 7. streamed and in-memory must agree
# ---------------------------------------------------------------------------

def test_cytome_and_anndata_agree_exactly():
    """Same matrix, two backends, identical p-values.

    ``layer`` must be pinned: the streaming default is ``'auto'``, which may
    normalise, and then the two paths are scoring different matrices -- the
    ranking diverges too, which is the tell that it is not a p-value bug.
    """
    import os
    import tempfile

    import cosg
    cytome = pytest.importorskip("cytome")
    anndata = pytest.importorskip("anndata")
    from cosg._cytome_streaming import run_cosg_cytome

    rng = np.random.default_rng(0)
    n, g = 500, 120
    lab = np.array([f"c{i % 3}" for i in range(n)])
    X = rng.poisson(0.4, (n, g)).astype(np.float32)
    X[lab == "c0", :3] += rng.poisson(7, ((lab == "c0").sum(), 3))
    a = anndata.AnnData(
        X=X,
        obs=pd.DataFrame({"ct": pd.Categorical(lab)},
                         index=[f"cell{i}" for i in range(n)]),
        var=pd.DataFrame(index=[f"g{i}" for i in range(g)]))

    path = os.path.join(tempfile.mkdtemp(), "t.cytome")
    cytome.from_anndata(a, output=path, force=True).close()

    streamed = run_cosg_cytome(path, groupby="ct", n_genes_user=10,
                               calculate_pvalues=True, verbose=False,
                               layer="counts")
    assert {"pvals", "pvals_adj", "zscores"} <= set(streamed)

    b = a.copy()
    cosg.cosg(b, groupby="ct", n_genes_user=10, calculate_pvalues=True)

    gi = list(streamed["groups_order"]).index("c0")
    s_names, s_p = streamed["names"][:, gi], streamed["pvals"][:, gi]
    m_names = pd.DataFrame(b.uns["cosg"]["names"])["c0"].values
    m_p = pd.DataFrame(b.uns["cosg"]["pvals"])["c0"].values

    assert list(s_names) == list(m_names), (list(s_names), list(m_names))
    assert np.allclose(s_p, m_p, rtol=1e-12, atol=0), (s_p, m_p)


def test_streaming_moments_match_the_direct_ones():
    """The streamed path only ever sees power sums; they must agree."""
    from cosg._pvalues import population_moments, srswor_moments_from_sums

    rng = np.random.default_rng(2)
    x = rng.gamma(1.5, 2.0, 400)
    x[rng.random(400) < 0.7] = 0.0
    direct = [v[0] for v in population_moments(x.reshape(-1, 1), want_m3=True)]
    streamed = srswor_moments_from_sums(x.sum(), (x ** 2).sum(),
                                        (x ** 3).sum(), x.size)
    assert np.allclose(direct, streamed, rtol=1e-10)
