"""Stratified p-values (``batch_key`` with ``calculate_pvalues``) leave out the
strata that carry no information about a group.

Under the stratified null, labels are permuted within each batch. A batch that
holds none of a group's cells, or nothing but them, has a fixed group sum (0, or
the batch's total), with zero variance. COSG added such a batch's sum to the
statistic ``T`` but not to its mean ``E``, so a group alone in a batch had
``T - E`` shifted by that batch's total: on null data every expressed gene came
out near p = 1e-80.

The fix leaves those strata out of ``T``, ``E`` and ``V`` together, and out of
the saddlepoint's blocks. It equals the review's fix (add ``nc * mean`` to ``E``
with zero variance when the group is the whole batch) in exact arithmetic: the
group's sum over a batch it fills is the batch's total, ``nb * mean``, so both
change ``T - E`` and ``V`` by exactly nothing; a point-mass block shifts the
saddlepoint's null and its target by the same constant. Leaving the strata out
also avoids the rounding residue of ``T_b - nb * mean_b`` and the degenerate
saddlepoint blocks (inclusion probability 1).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import cosg
from cosg._pvalues import analytic_pvalues, permutation_pvalues, population_moments, group_sums


def _null(seed=1, n1=600, n2=400, f=200, rate=0.3):
    """Batch b1 mixes A and B; batch b2 holds only A. Expression is independent
    of the labels within each batch, and b2 is deeper than b1."""
    rng = np.random.default_rng(seed)
    lab = np.r_[rng.choice(["A", "B"], n1), np.repeat("A", n2)]
    bat = np.r_[np.repeat("b1", n1), np.repeat("b2", n2)]
    lam = np.r_[np.full(n1, 1.0), np.full(n2, 3.0)]
    X = (rng.random((n1 + n2, f)) < rate) * rng.poisson(lam[:, None], (n1 + n2, f))
    return sp.csr_matrix(X.astype(np.float64)), lab, bat


@pytest.mark.parametrize("method", ["spa", "normal"])
def test_a_group_alone_in_a_batch_is_calibrated(method):
    X, lab, bat = _null()
    P, Z, info = analytic_pvalues(X, lab, ["A", "B"], batch=bat, method=method)
    for j in range(2):
        assert 0.3 < np.median(P[:, j]) < 0.7, (j, np.median(P[:, j]))
        assert (P[:, j] < 1e-3).mean() <= 0.01
    # against the stratified permutation reference, gene by gene
    ref = permutation_pvalues(X, lab, ["A", "B"], n_permutations=999, batch=bat, random_seed=3)
    assert np.abs(P - ref).max() < 0.07, np.abs(P - ref).max()


def test_the_public_api_calls_no_marker_on_null_data():
    import anndata as ad
    X, lab, bat = _null()
    obs = pd.DataFrame({"grp": pd.Categorical(lab), "batch": bat}, index=[f"c{i}" for i in range(len(lab))])
    a = ad.AnnData(X.astype(np.float32), obs=obs)
    a.var_names = [f"g{j}" for j in range(X.shape[1])]
    for ccs in (None, 0):
        cosg.cosg(a, groupby="grp", batch_key="batch", calculate_pvalues=True, n_genes_user=X.shape[1],
                  cpu_chunk_size=ccs)
        for g in ("A", "B"):
            padj = np.asarray(a.uns["cosg"]["pvals_adj"][g], dtype=float)
            assert np.nansum(padj < 0.05) <= 2, (ccs, g, np.nansum(padj < 0.05))


def test_equals_the_reviews_fix():
    """The review's formula -- keep the batch's sum in T, add nc * mean to E
    with zero variance when the group fills the batch -- gives the same z."""
    X, lab, bat = _null(f=60)
    order = ["A", "B"]
    _, Z, _ = analytic_pvalues(X, lab, order, batch=bat, method="normal")
    T = np.zeros((X.shape[1], 2)); E = np.zeros_like(T); V = np.zeros_like(T)
    for b in np.unique(bat):
        m = bat == b
        nb = m.sum()
        mean_b, s2_b, _ = population_moments(X[m])
        Tb = group_sums(X[m], lab[m], order)
        T += Tb
        for j, g in enumerate(order):
            nc = (lab[m] == g).sum()
            if nc == 0:
                continue
            E[:, j] += nc * mean_b                       # also when nc == nb
            if nc < nb:
                V[:, j] += nc * (nb - nc) / (nb - 1) * s2_b
    Z_review = np.where(V > 0, (T - E) / np.sqrt(np.where(V > 0, V, 1)), 0.0)
    np.testing.assert_allclose(Z, Z_review, rtol=1e-9, atol=1e-9)


def test_the_saddlepoint_leaves_out_the_same_strata():
    """A sparse marker of A inside the mixed batch is refined by the
    saddlepoint; the batch A fills must not enter its null. Checked against a
    long stratified permutation run."""
    rng = np.random.default_rng(7)
    n1, n2 = 600, 400
    lab = np.r_[rng.choice(["A", "B"], n1, p=[0.1, 0.9]), np.repeat("A", n2)]
    bat = np.r_[np.repeat("b1", n1), np.repeat("b2", n2)]
    X = np.zeros((n1 + n2, 3))
    inA = np.flatnonzero(lab[:n1] == "A")
    inB = np.flatnonzero(lab[:n1] == "B")
    X[inA[:4], 0] = 1.0                                 # 4 of A's b1 cells
    X[inB[:8], 0] = 1.0                                 # and 8 B cells
    X[n1:, 0] = rng.integers(0, 2, n2)                  # b2: expression of its own
    X[:, 1] = rng.integers(0, 2, n1 + n2)               # null genes
    X[:, 2] = rng.integers(0, 3, n1 + n2)
    X = sp.csr_matrix(X)
    P, _, info = analytic_pvalues(X, lab, ["A", "B"], batch=bat, method="spa")
    assert info["n_spa"] > 0
    ref = permutation_pvalues(X, lab, ["A", "B"], n_permutations=20000, batch=bat, random_seed=11)
    assert ref[0, 0] < 0.05                              # a real, moderate signal
    assert abs(np.log10(P[0, 0]) - np.log10(ref[0, 0])) < 0.25, (P[0, 0], ref[0, 0])


def test_a_two_cell_stratum_counts_as_the_permutation_null_does():
    """Strata of fewer than three cells used to be left out of the statistic
    whatever they held; now a stratum counts whenever it informs the group, as
    in the permutation null."""
    rng = np.random.default_rng(5)
    n = 402
    bat = np.r_[np.repeat("b1", 400), ["tiny", "tiny"]]
    lab = np.r_[rng.choice(["A", "B"], 400), ["A", "B"]]
    X = sp.csr_matrix(rng.poisson(1.0, (n, 40)).astype(float))
    P, _, _ = analytic_pvalues(X, lab, ["A", "B"], batch=bat, method="normal")
    ref = permutation_pvalues(X, lab, ["A", "B"], n_permutations=999, batch=bat, random_seed=2)
    assert np.abs(P - ref).max() < 0.07
