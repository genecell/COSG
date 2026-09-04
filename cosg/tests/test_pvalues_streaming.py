"""The streaming path's share of the significance layer.

`test_analytic_pvalues.py` covers the statistics. This file covers the parts
that only exist because the cytome path never holds the matrix: moments
rebuilt from streamed power sums, the decision about which entries need a tail
refinement, and the second pass that supplies the values for it.
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse

from cosg._pvalues import (
    population_moments, srswor_moments_from_sums, normal_tail, refine_tail,
    _exact_tail,
)


def test_moments_from_streamed_sums_match_the_in_memory_pass():
    """Streaming may not change the numbers, only where they come from."""
    rng = np.random.default_rng(23)
    X = sparse.csr_matrix(rng.poisson(0.5, (300, 40)).astype(np.float64))
    mean_a, s2_a, m3_a = population_moments(X, want_m3=True)

    sq = X.copy(); sq.data = sq.data ** 2
    cu = X.copy(); cu.data = cu.data ** 3
    mean_b, s2_b, m3_b = srswor_moments_from_sums(
        np.asarray(X.sum(axis=0)).ravel(),
        np.asarray(sq.sum(axis=0)).ravel(),
        np.asarray(cu.sum(axis=0)).ravel(),
        X.shape[0])

    assert np.allclose(mean_a, mean_b)
    assert np.allclose(s2_a, s2_b)
    assert np.allclose(m3_a, m3_b, atol=1e-9)


def test_power_sums_add_across_chunks():
    """What the chunk loop actually does: accumulate, then divide once."""
    rng = np.random.default_rng(24)
    X = sparse.csr_matrix(rng.poisson(0.6, (250, 30)).astype(np.float64))
    whole = population_moments(X, want_m3=True)

    s1 = np.zeros(X.shape[1]); s2 = np.zeros(X.shape[1]); s3 = np.zeros(X.shape[1])
    for lo in range(0, X.shape[0], 37):
        chunk = X[lo:lo + 37]
        sq = chunk.copy(); sq.data = sq.data ** 2
        cu = chunk.copy(); cu.data = cu.data ** 3
        s1 += np.asarray(chunk.sum(axis=0)).ravel()
        s2 += np.asarray(sq.sum(axis=0)).ravel()
        s3 += np.asarray(cu.sum(axis=0)).ravel()
    streamed = srswor_moments_from_sums(s1, s2, s3, X.shape[0])

    for a, b in zip(whole, streamed):
        assert np.allclose(a, b, atol=1e-9)


def test_flagging_refines_what_could_still_change_a_decision():
    """The trigger, pinned in both directions.

    A sparse gene whose normal p is small must be refined: that is the regime
    the normal approximation gets wrong, and it is why the refinement exists.
    A sparse gene whose normal p is nowhere near a threshold must NOT be --
    for a group smaller than half the cells the sum is right-skewed, so the
    true p exceeds the normal p and no refinement can make it significant.
    Flagging it anyway was what put 137,954 of 724,608 entries through the
    solver on a 10,000-cell dataset.
    """
    n = 5000
    # both genes are sparse (k = 200 of 5,000, expecting 20 in a 500-cell
    # group); they differ only in how extreme the observed group sum is.
    k = np.array([200, 200])
    mean = np.array([0.04, 0.04])
    s2 = np.array([0.5, 0.5])
    T = np.array([[60.0], [21.0]])       # far out; barely above expectation

    P, _, flagged = normal_tail(T, np.array([500.0]), n, mean, s2, None,
                                k_nonzero=k, method="spa")
    genes = set(flagged[:, 0].tolist())
    assert P[0, 0] <= 0.05 and 0 in genes, "a sparse gene in the tail must be refined"
    assert P[1, 0] > 0.05 and 1 not in genes, (
        "a sparse gene the normal approximation already puts in the bulk cannot "
        "become significant after refinement, so refining it is wasted work")


def test_normal_method_flags_nothing():
    """`pvalue_method='normal'` is the escape hatch; it must not refine."""
    n = 5000
    _, _, flagged = normal_tail(
        np.array([[20.0]]), np.array([500.0]), n,
        np.array([0.04]), np.array([0.5]), np.array([0.1]),
        k_nonzero=np.array([200]), method="normal")
    assert flagged.shape == (0, 2)


def test_refine_tail_enumerates_the_sparsest_genes():
    """Below the enumeration ceiling the answer is exact, not approximated."""
    rng = np.random.default_rng(29)
    vals = rng.gamma(2.0, 1.0, 8)
    n, nc = 500, 50
    t = float(vals.sum() * 0.5)
    assert np.isclose(refine_tail(vals, n, nc, t), _exact_tail(vals, n, nc, t))


def test_refine_tail_declines_rather_than_guesses():
    """No values, or a degenerate group, returns None so the caller keeps the
    normal p rather than being handed a number with nothing behind it."""
    assert refine_tail(np.empty(0), 100, 10, 1.0) is None
    assert refine_tail(np.array([1.0, 2.0]), 100, 0, 1.0) is None
    assert refine_tail(np.array([1.0, 2.0]), 100, 100, 1.0) is None


@pytest.mark.parametrize("nc", [10, 250])
def test_streamed_and_in_memory_agree_end_to_end(nc):
    """Same matrix, both routes to the normal p: identical to float error."""
    from cosg._pvalues import analytic_pvalues

    rng = np.random.default_rng(31)
    n, g = 500, 25
    labels = np.array(["A"] * nc + ["B"] * (n - nc))
    X = sparse.csr_matrix(rng.poisson(0.5, (n, g)).astype(np.float64))

    P_mem, _, _ = analytic_pvalues(X, labels, ["A", "B"], method="normal")

    mean, s2, m3 = population_moments(X)
    T = np.zeros((g, 2))
    for j, grp in enumerate(["A", "B"]):
        T[:, j] = np.asarray(X[labels == grp].sum(axis=0)).ravel()
    nc_all = np.array([(labels == "A").sum(), (labels == "B").sum()], dtype=float)
    P_str, _, _ = normal_tail(T, nc_all, n, mean, s2, m3,
                              k_nonzero=np.diff(X.tocsc().indptr),
                              method="normal")
    assert np.allclose(P_mem, P_str, rtol=1e-10, atol=1e-300)
