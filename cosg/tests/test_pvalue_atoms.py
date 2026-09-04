"""The atom summary, the lattice correction, and the solver that reads them.

Each test pins a claim that was measured before it was implemented, and
several pin the *premise* as well: if a fixture stops exercising the regime a
correction exists for, the test says so rather than passing quietly.
"""
import itertools

import numpy as np
import pytest
from scipy import sparse

from cosg._pvalue_atoms import (
    GRID_NBINS, AtomSet, StreamingAtomAccumulator, atoms_from_columns,
    exact_tail_atoms, spa_atoms_tail, value_bin_index, values_are_integral,
)
from cosg._pvalues import _spa_tail, analytic_pvalues


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def _atoms_of(v):
    b = value_bin_index(v)
    c = np.bincount(b, minlength=GRID_NBINS)
    s = np.bincount(b, weights=v, minlength=GRID_NBINS)
    nz = c > 0
    return s[nz] / c[nz], c[nz].astype(float)


def _one_tail(v, n, n_c, t, lattice):
    av, aw = _atoms_of(np.asarray(v, dtype=np.float64))
    p, lp, ok = spa_atoms_tail(
        av, aw, np.array([0, av.size]), np.array([0]),
        np.array([float(n - len(v))]), np.array([float(n_c)]),
        np.array([float(t)]), lattice=lattice)
    return float(p[0]), float(lp[0]), bool(ok[0])


def _deep_t(v, n, n_c, z):
    v = np.asarray(v, dtype=np.float64)
    mean = v.sum() / n
    var = (v ** 2).sum() / n - mean ** 2
    return mean * n_c + z * np.sqrt(n_c * (n - n_c) / (n - 1) * var)


# --------------------------------------------------------------------------
# the exact tail, which everything else is measured against
# --------------------------------------------------------------------------

def test_exact_tail_equals_brute_force_enumeration():
    """The oracle is exact, not merely close: every subset, counted."""
    rng = np.random.default_rng(11)
    n, k, n_c = 16, 6, 5
    x = np.zeros(n)
    idx = rng.choice(n, k, replace=False)
    x[idx] = rng.geometric(0.4, k)
    vals, cnts = np.unique(x[idx], return_counts=True)
    for t in (3, 6, 9, 12):
        brute = np.mean([x[list(c)].sum() >= t
                         for c in itertools.combinations(range(n), n_c)])
        got = exact_tail_atoms(vals, cnts.astype(float), n, n_c, t, k)
        assert abs(brute - got) < 1e-12, (t, brute, got)


def test_exact_tail_refuses_non_integer_values():
    """It is a lattice construction; on continuous values it must decline
    rather than silently round."""
    v = np.array([0.3, 1.7, 2.2])
    assert exact_tail_atoms(v, np.ones(3), 100, 10, 3.0, 3) is None


# --------------------------------------------------------------------------
# compression: the summary changes the answer by less than the method's error
# --------------------------------------------------------------------------

@pytest.mark.parametrize("k,n_c", [(60, 500), (200, 500), (2000, 500), (10000, 1000)])
def test_atoms_reproduce_the_per_value_tail_on_normalised_data(k, n_c):
    """Atoms are a compression of the value multiset, and the saddlepoint
    reads a gene only through sums over that multiset — so the compressed and
    uncompressed tails must agree far inside the method's own error."""
    rng = np.random.default_rng(4 + k)
    n = 20_000
    v = np.log1p(rng.geometric(0.4, k) / rng.uniform(0.5, 2.0, k) * 1.3)
    t = _deep_t(v, n, n_c, 4.5)
    p_atoms, _, ok = _one_tail(v, n, n_c, t, lattice=False)
    p_values = _spa_tail(v, n, n_c, t)
    assert ok
    assert abs(p_atoms / p_values - 1.0) < 0.01, (p_atoms, p_values)
    # and the compression is real: far fewer atoms than values
    assert _atoms_of(v)[0].size < max(80, k // 10)


def test_the_summary_keeps_each_bin_mean_exactly():
    """Storing count *and* sum makes the atom the bin's own mean, so the first
    moment is preserved within every bin and the error is second order."""
    rng = np.random.default_rng(2)
    v = rng.gamma(2.0, 1.0, 5000)
    av, aw = _atoms_of(v)
    assert abs(float(av @ aw) - float(v.sum())) < 1e-8 * abs(float(v.sum()))
    assert int(aw.sum()) == v.size


# --------------------------------------------------------------------------
# the lattice correction
# --------------------------------------------------------------------------

@pytest.mark.parametrize("k,n_c", [(60, 500), (200, 500), (600, 1000), (2000, 500)])
def test_lattice_correction_recovers_the_exact_count_tail(k, n_c):
    """On integer data the sum lives on a lattice, and the continuous formula
    approximates a step function from the wrong side. The `t - 1/2` evaluation
    is what closes the gap."""
    rng = np.random.default_rng(5 + k)
    n = 20_000
    v = rng.geometric(0.4, k).astype(float)
    t = np.ceil(_deep_t(v, n, n_c, 4.5))
    vals, cnts = np.unique(v, return_counts=True)
    e_h = k * n_c / n
    h_max = int(min(k, n_c, e_h + 10 * np.sqrt(e_h) + 25))
    exact = exact_tail_atoms(vals, cnts.astype(float), n, n_c, t, h_max)

    p_lattice, _, ok_l = _one_tail(v, n, n_c, t, lattice=True)
    p_plain, _, _ = _one_tail(v, n, n_c, t, lattice=False)
    assert ok_l
    assert 0.90 <= p_lattice / exact <= 1.10, (p_lattice, exact)
    # The premise: without the correction the same row is materially
    # conservative. If this ever fails the fixture stopped exercising the
    # regime the correction exists for, and the test above became vacuous.
    assert p_plain / exact < 0.96, (p_plain, exact)


def test_integrality_is_detected_from_the_data():
    assert values_are_integral(np.array([1.0, 4.0, 17.0]))
    assert not values_are_integral(np.array([1.0, 4.5]))
    assert values_are_integral(np.array([], dtype=float))


# --------------------------------------------------------------------------
# the solver: convergence, the fallback, and log space
# --------------------------------------------------------------------------

def test_no_refined_row_falls_back_to_the_normal_p():
    """The failure mode this replaced: the vectorised solver gave up on the
    sparsest rows and the caller silently kept a normal p-value that was
    wrong by up to nine orders of magnitude — in exactly the regime the
    refinement exists for."""
    rng = np.random.default_rng(1)
    n = 20_000
    n_rows = 0
    vals, wts, blk_ptr, blk_row, nz, nc, tt = [], [], [0], [], [], [], []
    pos = 0
    for _ in range(200):
        k = int(rng.integers(20, 120))
        v = rng.geometric(0.4, k).astype(float)
        n_c = 500
        t = np.ceil(_deep_t(v, n, n_c, rng.uniform(3.0, 6.0)))
        av, aw = _atoms_of(v)
        vals.append(av); wts.append(aw)
        pos += av.size
        blk_ptr.append(pos); blk_row.append(n_rows)
        nz.append(float(n - k)); nc.append(float(n_c)); tt.append(float(t))
        n_rows += 1
    p, lp, ok = spa_atoms_tail(
        np.concatenate(vals), np.concatenate(wts),
        np.asarray(blk_ptr), np.asarray(blk_row),
        np.asarray(nz), np.asarray(nc), np.asarray(tt), lattice=True)
    assert ok.all(), f"{(~ok).sum()} of {n_rows} rows failed to converge"
    assert np.all((p > 0) & (p <= 1))


def test_log_p_agrees_with_p_and_survives_underflow():
    """`p` saturates at the float64 floor; the log form must not."""
    rng = np.random.default_rng(7)
    n, k, n_c = 50_000, 3000, 2000
    v = rng.geometric(0.3, k).astype(float)
    for z in (4.0, 12.0, 40.0):
        t = np.ceil(_deep_t(v, n, n_c, z))
        p, lp, ok = _one_tail(v, n, n_c, t, lattice=True)
        assert ok
        if p > 1e-300:
            assert abs(lp - np.log(p)) < 1e-6 * abs(np.log(p))
    # Deep enough that the float64 p is exactly zero. The log form must keep
    # going — and keep *ranking*, which is the whole point of reporting it.
    deep = [_one_tail(v, n, n_c, np.ceil(_deep_t(v, n, n_c, z)), lattice=True)
            for z in (100.0, 150.0)]
    assert all(ok for _p, _lp, ok in deep)
    assert all(_p == 0.0 for _p, _lp, _ok in deep), "fixture no longer underflows"
    assert deep[0][1] < np.log(1e-300)
    assert deep[1][1] < deep[0][1], "log p must still discriminate below the floor"


def test_a_row_with_more_mass_is_more_significant():
    """Monotonicity — a sanity property no approximation may break."""
    rng = np.random.default_rng(3)
    n, k, n_c = 20_000, 300, 500
    v = rng.geometric(0.4, k).astype(float)
    base = _deep_t(v, n, n_c, 3.0)
    ps = [_one_tail(v, n, n_c, np.ceil(base * m), lattice=True)[0]
          for m in (1.0, 1.2, 1.5)]
    assert ps[0] > ps[1] > ps[2]


# --------------------------------------------------------------------------
# the summary is additive: streamed and in-memory are identical
# --------------------------------------------------------------------------

def test_streaming_atoms_are_identical_to_in_memory_atoms():
    """Not 'close': the grid is fixed, so counts and sums per (gene, bin) add
    across chunks and the two backends see the same summary bit for bit."""
    rng = np.random.default_rng(8)
    n, g = 2000, 60
    X = sparse.csr_matrix((rng.random((n, g)) < 0.2) * rng.gamma(2.0, 1.0, (n, g)))
    want = np.array([1, 5, 9, 20, 44])
    direct = atoms_from_columns(X, want)
    for batch in (137, 500, 2000):
        acc = StreamingAtomAccumulator(want)
        for lo in range(0, n, batch):
            acc.add_chunk(X[lo:lo + batch])
        streamed = acc.finish()
        assert np.array_equal(streamed.ptr, direct.ptr)
        assert np.array_equal(streamed.wts, direct.wts)
        assert np.allclose(streamed.vals, direct.vals, rtol=0, atol=1e-12)


# --------------------------------------------------------------------------
# the trigger
# --------------------------------------------------------------------------

def test_trigger_skips_only_entries_that_cannot_become_significant():
    """A group under half the cells makes the sum right-skewed, so the true p
    exceeds the normal p and an entry already above the trigger cannot cross a
    threshold after refinement. Nothing with a small normal p is skipped."""
    rng = np.random.default_rng(6)
    n, g = 4000, 400
    X = sparse.csr_matrix(rng.poisson(0.3, (n, g)).astype(float))
    labels = np.array(["a"] * 300 + ["b"] * 1200 + ["c"] * 2500)
    P, Z, info = analytic_pvalues(X, labels, ["a", "b", "c"], method="spa")
    Pn, _, _ = analytic_pvalues(X, labels, ["a", "b", "c"], method="normal")
    changed = P != Pn
    # every changed entry had a normal p at or below the trigger
    assert np.all(Pn[changed] <= 0.05 + 1e-12)
    # and refinement never made a non-significant entry significant
    assert not np.any((Pn > 0.05) & (P <= 0.01))
    assert info["n_spa_failed"] == 0


def test_analytic_pvalues_reports_log_p_alongside_p():
    rng = np.random.default_rng(12)
    n, g = 1500, 120
    X = sparse.csr_matrix(rng.poisson(0.5, (n, g)).astype(float))
    labels = np.array(["a"] * 400 + ["b"] * 1100)
    P, Z, info = analytic_pvalues(X, labels, ["a", "b"], method="spa")
    assert "log_p" in info and info["log_p"].shape == P.shape
    fin = P > 1e-290
    assert np.allclose(info["log_p"][fin], np.log(P[fin]), rtol=1e-6)
    assert info["lattice"] is True          # Poisson counts are integral
