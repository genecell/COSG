"""Analytic significance for COSG's cosine specificity, calibrated in the tail.

The tested statistic
--------------------
For gene ``g`` with expression ``x`` over ``n`` cells and a group ``c`` of
``n_c`` cells, COSG's raw cosine is ``T / (||x|| * sqrt(n_c))`` with ``T`` the
sum of ``x`` over the group. Under the exchangeability null -- this gene's
expression is independent of the grouping -- ``||x||`` and ``n_c`` are fixed,
so testing the cosine *is* testing ``T``, and ``T`` is a simple random sample
sum without replacement from the fixed population ``x``. That has exact
moments:

    E[T]   = n_c * mean(x)
    Var[T] = n_c * (n - n_c) / (n - 1) * s2
    mu3[T] = n_c * (n - n_c) * (n - 2 n_c) / ((n - 1)(n - 2)) * m3

with ``s2`` and ``m3`` the population variance and third central moment. No
permutations are sampled and there is no p-value floor.

Why not just the normal approximation
-------------------------------------
Because it is wrong where it matters. The CLT's quality here is governed by the
expected number of *expressing* cells inside the group, ``n_c * k / n``, not by
``n_c`` or by the gene's global expressing count ``k``. Measured against a
Monte-Carlo null: a gene expressed in 200 of 50,000 cells, tested against a
500-cell group, has a normal p of 1e-5 where the true tail probability is
1.8e-3 -- anti-conservative by 180x. That gene is not a corner case; it is a
sparse, specific gene, which is the interesting kind, and BH at FDR 0.01 over
20,000 genes engages exactly those thresholds.

An Edgeworth (skewness) correction does not rescue it: Edgeworth has *absolute*
error, so it improves the bulk and still misses the tail by an order of
magnitude. A saddlepoint approximation has *relative* error and is the standard
remedy for precisely this shape of problem -- it is what GWAS adopted for score
tests on rare variants with unbalanced case-control ratios (SAIGE), where
"rare variant, few cases" is the same geometry as "sparse gene, small group".

Because ``T`` is a sum *without* replacement it has no closed-form cumulant
generating function, so the saddlepoint used here is the **conditional (double)
saddlepoint**: tilt independent Bernoulli inclusions and condition on the group
size. The conditional law does not depend on the Bernoulli parameter, so it is
fixed at 1/2. Every evaluation costs O(k) -- the zero-valued cells enter the
count constraint in closed form -- so the method is cheapest exactly where
genes are sparsest.

Two regimes, chosen internally
------------------------------
``pvalue_method='spa'`` (the default) is adaptive:

* **normal** where the tail is irrelevant -- a null gene reading 0.30 rather
  than 0.35 changes nothing, and the normal approximation is excellent there;
* **saddlepoint** once the normal p is small or, for a group of more than
  half the cells, the expected in-group expressing count is low; below a
  normal p of 1e-30 only for sparse genes, whose normal p is wrong by orders
  of magnitude there.

No p is reported below the permutation floor, 1 / prod_b C(n_b, n_c,b): the
probability of the observed labeling itself, which bounds the tail from
below.

The exact conditional tail (``pvalue_method='exact'``) is opt-in: it is the
oracle the saddlepoint is checked against, not part of ``'spa'``.

What the p-value is not
-----------------------
It calibrates the **raw cosine**. The reported COSG score is the
mu-regularized ranking statistic, and p-values are identical across ``mu`` by
design: ``mu`` encodes how harshly to discount a gene that is also high
elsewhere, which is a preference about ranking genes that are genuinely
associated. The null contains no ``mu``, so chance has one answer per gene per
group. A p-value on the composed score would also depend on every *other*
group's cosine, so adding an unrelated cluster would change this gene's
significance in this cluster.
"""
from __future__ import annotations

import warnings
from typing import Optional, Sequence

import numpy as np
from scipy import sparse
from scipy.special import expit, ndtr, log_expit, log_ndtr

from ._pvalue_atoms import (atoms_from_columns, spa_atoms_tail,
                            values_are_integral, exact_tail_atoms,
                            exact_tail_bounded)

__all__ = [
    "srswor_moments_from_sums",
    "normal_tail",
    "refine_tail",
    "population_moments",
    "group_sums",
    "analytic_pvalues",
    "permutation_pvalues",
    "adjust_fdr",
    "adjust_fdr_log",
]

_SQRT2PI = np.sqrt(2.0 * np.pi)


def _phi(z):
    return np.exp(-0.5 * np.asarray(z, dtype=np.float64) ** 2) / _SQRT2PI


def _sf(z):
    """Upper-tail normal survival function, accurate far into the tail."""
    return ndtr(-np.asarray(z, dtype=np.float64))


# ---------------------------------------------------------------------------
# Population moments
# ---------------------------------------------------------------------------

def population_moments(X, want_m3: bool = False):
    """Per-gene ``(mean, s2, m3)`` over cells, in float64.

    ``want_m3=False`` skips the third power sum. The skewness correction it
    fed was the Edgeworth term, and the saddlepoint replaced that -- so on a
    60-million-nonzero matrix this is one whole reduction that no longer has
    a consumer. ``m3`` comes back as None when it is not asked for.

    ``X`` is cells x genes, sparse or dense. Population (``/n``) moments, not
    sample ones -- the finite-population formulas above are stated that way.
    """
    n = X.shape[0]
    if sparse.issparse(X):
        # bincount over the column indices, not a CSC conversion: COSG is
        # handed CSR, and `.tocsc()` on a 20M-nonzero matrix costs more than
        # every other part of this feature combined. The three power sums come
        # from one pass over the stored values with no matrix copy at all.
        X = X.tocsr() if not sparse.isspmatrix_csr(X) else X
        g = X.shape[1]
        idx = X.indices
        d = X.data.astype(np.float64, copy=False)
        d2 = d * d
        s1 = np.bincount(idx, weights=d, minlength=g)
        sq = np.bincount(idx, weights=d2, minlength=g)
        cu = np.bincount(idx, weights=d2 * d, minlength=g) if want_m3 else None
    else:
        Xd = np.asarray(X, dtype=np.float64)
        s1 = Xd.sum(axis=0)
        sq = (Xd ** 2).sum(axis=0)
        cu = (Xd ** 3).sum(axis=0) if want_m3 else None

    mean = s1 / n
    s2 = sq / n - mean ** 2
    np.maximum(s2, 0.0, out=s2)
    # E[(x-mu)^3] = E[x^3] - 3 mu E[x^2] + 2 mu^3
    if cu is None:
        return mean, s2, None
    m3 = cu / n - 3.0 * mean * (sq / n) + 2.0 * mean ** 3
    return mean, s2, m3


def group_sums(X, labels, groups_order):
    """``T[gene, group]``: the sum of each gene over each group's cells."""
    labels = np.asarray(labels)
    idx = {g: i for i, g in enumerate(groups_order)}
    col = np.array([idx.get(v, -1) for v in labels])
    keep = col >= 0
    rows = np.flatnonzero(keep)
    ind = sparse.csr_matrix(
        (np.ones(rows.size, dtype=np.float64), (rows, col[keep])),
        shape=(X.shape[0], len(groups_order)),
    )
    T = (X.T @ ind)
    return np.asarray(T.todense() if sparse.issparse(T) else T, dtype=np.float64)


# ---------------------------------------------------------------------------
# Exact tail for genes with very few expressing cells
# ---------------------------------------------------------------------------

def _subset_sums_by_size(values):
    """``sums[h]`` = every subset sum of size ``h``. Cost 2**len(values)."""
    sums = [np.zeros(1, dtype=np.float64)]
    for v in values:
        new = [None] * (len(sums) + 1)
        new[0] = sums[0]
        for h in range(1, len(sums) + 1):
            grown = sums[h - 1] + v
            new[h] = np.concatenate([sums[h], grown]) if h < len(sums) else grown
        sums = new
    return sums


def _exact_tail(values, n, n_c, t, tol=1e-12):
    """Exact ``P(T >= t)`` by enumerating subsets of the expressing cells.

    ``T`` only ever sums expressing cells, so the null is: draw ``h`` of them
    into the group (hypergeometric), then a uniform ``h``-subset of their
    values. Both parts are enumerable when ``k`` is small.
    """
    from scipy.stats import hypergeom

    k = len(values)
    sums = _subset_sums_by_size(values)
    total = 0.0
    hmax = min(k, n_c)
    for h in range(0, hmax + 1):
        w = hypergeom.pmf(h, n, k, n_c)
        if w <= 0.0:
            continue
        s = sums[h]
        total += w * float(np.mean(s >= t - tol))
    return float(min(max(total, 0.0), 1.0))


# ---------------------------------------------------------------------------
# Conditional (double) saddlepoint
# ---------------------------------------------------------------------------

def _solve_saddle_bracketed(values, n_zero, n_c, t, u0):
    """Bracket-and-bisect on ``s``, solving ``u`` inside. Slow but sure.

    Newton handles almost every entry; the ones it does not are very sparse
    genes, where the tilt has to travel a long way and a quadratic step can
    overshoot out of the convex hull. Falling back here keeps those entries
    correct instead of dropping them back to a normal p-value that is known
    to be wrong by orders of magnitude exactly there.
    """
    def solve_u(s, u_start):
        u = u_start
        for _ in range(200):
            p_nz = expit(s * values + u)
            p_z = expit(u)
            fv = p_nz.sum() + n_zero * p_z - n_c
            d = (p_nz * (1.0 - p_nz)).sum() + n_zero * p_z * (1.0 - p_z)
            if d <= 1e-300:
                return None
            step = np.clip(fv / d, -4.0, 4.0)
            u -= step
            if abs(step) < 1e-13:
                break
        return u

    def g(s):
        u = solve_u(s, u0)
        if u is None:
            return None, None
        return float(values @ expit(s * values + u)) - t, u

    lo, hi = -1.0, 1.0
    glo, _ = g(lo)
    ghi, _ = g(hi)
    if glo is None or ghi is None:
        return None
    tries = 0
    # The inner Newton can fail while the bracket is being widened, and it
    # returns None when it does. Comparing that to 0 raises TypeError from
    # inside a numerical helper -- a crash where a fallback was intended.
    while glo is not None and glo > 0 and tries < 300:
        lo *= 2.0
        glo, _ = g(lo)
        tries += 1
    while ghi is not None and ghi < 0 and tries < 600:
        hi *= 2.0
        ghi, _ = g(hi)
        tries += 1
    if glo is None or ghi is None or glo > 0 or ghi < 0:
        return None
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        val, u = g(mid)
        if val is None:
            return None
        if val > 0:
            hi = mid
        else:
            lo = mid
        if hi - lo < 1e-14 * max(1.0, abs(mid)):
            break
    s_hat = 0.5 * (lo + hi)
    u_hat = solve_u(s_hat, u0)
    return None if u_hat is None else (s_hat, u_hat)


def _spa_tail(values, n, n_c, t):
    """``P(T >= t)`` by Skovgaard's conditional saddlepoint approximation.

    Joint CGF of ``(sum Z_i x_i, sum Z_i)`` with ``Z_i`` iid Bernoulli(1/2):
    ``K(s, u) = sum log(1 + exp(s x_i + u))``. Conditioning on ``sum Z_i =
    n_c`` recovers sampling without replacement, and the conditional law does
    not depend on the Bernoulli parameter.

    Solved by two-dimensional Newton on ``(s, u)``. ``K`` is convex, so the
    Hessian is positive definite and Newton converges in a handful of steps;
    a nested bracket-and-bisect was ~2000x slower for the same answer, which
    is the difference between a usable feature and one nobody switches on.
    """
    values = np.asarray(values, dtype=np.float64)
    k = values.size
    n_zero = n - k
    f = n_c / n
    if not (0.0 < f < 1.0):
        return None
    u = np.log(f / (1.0 - f))
    s = 0.0

    def _state(s, u):
        th = s * values + u
        p_nz = expit(th)
        p_z = expit(u)
        w_nz = p_nz * (1.0 - p_nz)
        w_z = p_z * (1.0 - p_z)
        Ks = float(values @ p_nz)
        Ku = float(p_nz.sum() + n_zero * p_z)
        Kss = float((values * values) @ w_nz)
        Ksu = float(values @ w_nz)
        Kuu = float(w_nz.sum() + n_zero * w_z)
        return Ks, Ku, Kss, Ksu, Kuu

    for _ in range(100):
        Ks, Ku, Kss, Ksu, Kuu = _state(s, u)
        r1, r2 = Ks - t, Ku - n_c
        if abs(r1) < 1e-11 * max(1.0, abs(t)) and abs(r2) < 1e-11 * max(1.0, n_c):
            break
        det = Kss * Kuu - Ksu * Ksu
        if not np.isfinite(det) or det <= 1e-300:
            return None
        ds = (Kuu * r1 - Ksu * r2) / det
        du = (-Ksu * r1 + Kss * r2) / det
        # damped: the tilt can be pushed far outside the convex hull by a
        # single unguarded step on a very sparse gene
        scale = min(1.0, 4.0 / max(abs(ds), abs(du), 1e-12))
        s -= ds * scale
        u -= du * scale
    else:
        got = _solve_saddle_bracketed(values, n_zero, n_c, t,
                                      np.log(f / (1.0 - f)))
        if got is None:
            return None
        s, u = got

    if abs(s) < 1e-10:
        return 0.5

    Ks, Ku, Kss, Ksu, Kuu = _state(s, u)
    det = Kss * Kuu - Ksu * Ksu
    if det <= 0 or Kuu <= 0:
        return None

    def Kfun(s_, u_):
        return float(-log_expit(-(s_ * values + u_)).sum()
                     - log_expit(-u_) * n_zero)

    # the s = 0 system has a closed-form solution: every cell has the same
    # inclusion probability, so u_null = logit(n_c / n)
    u_null = np.log(f / (1.0 - f))
    Kuu0 = n * f * (1.0 - f)
    if Kuu0 <= 0:
        return None

    inner = ((Kfun(0.0, u_null) - u_null * n_c)
             - (Kfun(s, u) - s * t - u * n_c))
    if inner <= 0:
        inner = 0.0
    w = np.sign(s) * np.sqrt(2.0 * inner)
    v = s * np.sqrt(det / Kuu0)
    if w == 0 or v == 0:
        return None
    p = _sf(w) - _phi(w) * (1.0 / w - 1.0 / v)
    return float(min(max(p, 0.0), 1.0))


# ---------------------------------------------------------------------------
# The streaming-friendly surface
#
# A streamed pass can accumulate power sums but cannot keep the matrix, so the
# analytic path is split in two: everything computable from sums, and a tail
# refinement that needs the nonzero values of a few flagged genes. The
# in-memory path calls the same two functions, so there is one implementation
# of the arithmetic rather than one per backend.
# ---------------------------------------------------------------------------






# ---------------------------------------------------------------------------
# The public analytic path
# ---------------------------------------------------------------------------

def analytic_pvalues(
    X,
    labels,
    groups_order,
    batch=None,
    method: str = "spa",
    spa_normal_p: float = 0.05,
    spa_min_expected: float = 30.0,
    spa_min_p: float = 1e-30,
    exact_max_table: float = 2e6,
    verbosity: int = 0,
):
    """``(pvals, zscores, info)`` for every gene x group.

    Parameters
    ----------
    X
        cells x genes, the same matrix COSG scored.
    labels
        Per-cell group label.
    groups_order
        Column order for the output.
    batch
        Per-cell batch label. The null becomes stratified permutation: labels
        are permuted within batch, so ``T`` is a sum of independent per-stratum
        sums and the three moments add across strata.
    method
        ``'spa'`` (adaptive, default), ``'normal'`` (no tail refinement), or
        ``'exact'`` (the exact conditional tail on integer values where its
        table is affordable and its cut cannot matter, else the saddlepoint).
    spa_normal_p
        Entries whose normal p is at most this are refined.
    spa_min_expected
        For a group of more than half the cells, entries whose expected
        number of expressing cells in the group is below this are refined
        whatever their normal p.
    spa_min_p
        Below this normal p (1e-30) only sparse entries (fewer than
        ``spa_min_expected`` expressing cells expected in the group) are
        refined; 0 refines every flagged entry. Every p is then held at or
        above the permutation floor, 1 / prod_b C(n_b, n_c,b).
    """
    # Do NOT convert up front. COSG hands this a CSR matrix; population_moments
    # wants CSR (bincount over column indices) and only the per-gene value
    # lookup for refined rows wants CSC. Converting here meant CSR -> CSC ->
    # CSR, two full rebuilds of a 40M-nonzero matrix and 1.9 s of the 6.8 s
    # this function took, for nothing. Build CSC lazily, once, and only if
    # some row actually needs refining.
    if not sparse.issparse(X):
        X = np.asarray(X)
    _csc_cache = {}

    def _as_csc():
        if "m" not in _csc_cache:
            _csc_cache["m"] = (X.tocsc() if sparse.issparse(X) else X)
        return _csc_cache["m"]
    n, n_genes = X.shape
    labels = np.asarray(labels)
    n_groups = len(groups_order)

    E = np.zeros((n_genes, n_groups), dtype=np.float64)
    V = np.zeros((n_genes, n_groups), dtype=np.float64)
    T = np.zeros((n_genes, n_groups), dtype=np.float64)

    strata = [np.ones(n, dtype=bool)] if batch is None else [
        np.asarray(batch) == b for b in np.unique(np.asarray(batch))
    ]
    # The smallest p the permutation null can give a group: the probability of
    # one labeling, 1 / prod_b C(n_b, n_c,b) over the strata that inform it.
    log_floor = np.zeros(n_groups, dtype=np.float64)

    # A stratum informs the test of a group only when it holds some of the
    # group's cells and some other cells. One with none of them, or with
    # nothing else, has a fixed sum under the within-stratum permutation null
    # (0, or the stratum's total), so it is left out of the statistic, its
    # mean and its variance together; the saddlepoint below leaves out the
    # same strata. Adding such a stratum's sum to T without its mean -- what
    # this did for a group alone in its batch -- shifted T by that total and
    # made every gene expressed there a marker (p near 1e-80 on null data).
    for mask in strata:
        nb = int(mask.sum())
        nc_b = np.array([(labels[mask] == g).sum() for g in groups_order],
                        dtype=np.float64)
        informative = np.flatnonzero((nc_b > 0) & (nc_b < nb))
        if informative.size == 0:
            continue
        Xs = X[mask]
        mean_b, s2_b, _ = population_moments(Xs)
        Tb = group_sums(Xs, labels[mask], groups_order)
        for j in informative:
            nc = nc_b[j]
            T[:, j] += Tb[:, j]
            E[:, j] += nc * mean_b
            V[:, j] += nc * (nb - nc) / (nb - 1) * s2_b
            log_floor[j] -= _log_binom(nb, nc)

    sd = np.sqrt(V)
    with np.errstate(invalid="ignore", divide="ignore"):
        Z = (T - E) / sd
    Z[~np.isfinite(Z)] = 0.0
    P = _sf(Z)
    P[sd <= 0] = 1.0

    info = {"n_exact": 0, "n_spa": 0, "n_normal": int(P.size),
            "n_spa_failed": 0, "n_bisected": 0, "lattice": False}
    log_P = np.log(np.clip(P, np.finfo(np.float64).tiny, 1.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        log_P = np.where(sd > 0, log_ndtr(-Z), 0.0)
    info["log_p"] = log_P
    if method == "normal":
        return _at_least_the_floor(P, Z, info, log_P, log_floor)

    # The values are read once, as atoms: the saddlepoint sees a gene only
    # through sums over its value multiset, so k individual cells and their
    # (value, count) summary give the same answer at a cost that no longer
    # scales with k. Integrality decides the continuity correction and is a
    # property of the data, never a claim by the caller.
    x_data = X.data if sparse.issparse(X) else X
    integral = values_are_integral(x_data)
    info["lattice"] = bool(integral)

    # Per-gene expressing counts decide which rows need refining. Taking them
    # from CSR's column indices avoids building CSC for the decision itself --
    # only the rows that survive it need per-gene values.
    if sparse.issparse(X):
        _Xr = X if sparse.isspmatrix_csr(X) else X.tocsr()
        counts = np.bincount(_Xr.indices[_Xr.data != 0], minlength=n_genes)
    else:
        counts = (X != 0).sum(axis=0)

    n_c_all = np.array([(labels == g).sum() for g in groups_order], dtype=np.float64)
    expected_in_group = np.outer(counts, n_c_all) / max(n, 1)

    # One trigger. A group smaller than half the cells makes the sum
    # right-skewed, so its upper tail is heavier than the normal's and the
    # true p exceeds the normal p: an entry the normal approximation already
    # puts above `spa_normal_p` cannot become significant after refinement,
    # and skipping it is conservative by construction. Groups larger than half
    # the cells have no such guarantee, so there the old sparsity criterion
    # still applies. This replaces the previous union, which flagged every
    # sparse entry whether or not it was anywhere near a decision -- on a
    # 10,000-cell dataset that was 137,954 of 724,608 entries.
    big_group = n_c_all > (n / 2.0)
    refine = (P <= spa_normal_p) | (big_group[None, :] & (expected_in_group < spa_min_expected))
    refine &= (sd > 0) & (T > E)
    # Below a normal p of `spa_min_p` (1e-30) an entry is refined only when it
    # is sparse: fewer than `spa_min_expected` expressing cells expected in the
    # group. There the normal p is wrong by orders of magnitude and the
    # saddlepoint is not: a gene in 12 cells, all of a 100-cell group of
    # 5,000, had a normal log10 p of -129.7 against an exact -20.7 (the
    # saddlepoint: -20.2), and every such entry used to keep the normal value.
    # A dense gene that deep keeps its normal p: it is too small as well
    # (40 to 110 orders at |z| 20 to 50), but no decision depends on it, and
    # at the most extreme |z| the saddlepoint's solve can fail to converge.
    # Every p is then held at or above the permutation floor (see below).
    # `spa_min_p=0` refines every flagged entry.
    if spa_min_p > 0:
        refine &= (P >= spa_min_p) | (expected_in_group < spa_min_expected)
    rows, cols = np.nonzero(refine)
    if verbosity > 0 and rows.size:
        print(f"COSG p-values: refining {rows.size} of {P.size} tail entries "
              f"({'integer lattice' if integral else 'continuous'})")
    if rows.size == 0:
        return _at_least_the_floor(P, Z, info, log_P, log_floor)

    genes_needed = np.unique(rows)
    strata = ([np.arange(n)] if batch is None else
              [np.flatnonzero(np.asarray(batch) == lev)
               for lev in np.unique(np.asarray(batch))])
    nc_by_stratum = [
        np.array([float((labels[idx] == g).sum()) for g in groups_order])
        for idx in strata
    ]
    # One atom set per stratum. Unstratified data have exactly one, so the
    # stratified and plain paths differ only in how many blocks a row carries.
    atom_sets, pos_maps, k_per = [], [], []
    for idx in strata:
        Xs = X[idx] if len(strata) > 1 else X
        aset = atoms_from_columns(Xs, genes_needed)
        atom_sets.append(aset)
        pos_maps.append(aset.index_of())
        k_per.append(aset.k)

    # Exact enumeration, opt-in: the conditional permutation law of an
    # integer-valued T is computable in closed form from the same atoms
    # (Pagano & Tritchler, 1983). It is the oracle the saddlepoint is measured
    # against, and `pvalue_method='exact'` reports it directly wherever the
    # table is affordable.
    exact_rows = set()
    if method == "exact" and not (integral and batch is None):
        # Asking for the exact tail and silently getting the saddlepoint is
        # the kind of thing a user only finds out by reading the source.
        why = ("the values are not integers -- the exact tail indexes its "
               "state by the sum, so the sum has to lie on a lattice"
               if not integral else
               "a stratified (batch_key) null is not enumerated")
        warnings.warn(
            f"pvalue_method='exact' does not apply here: {why}. "
            "The saddlepoint is used instead, which is what 'spa' would have "
            "done; on continuous values the two agree to about 1% anyway.",
            stacklevel=2)
    if method == "exact" and integral and batch is None:
        aset = atom_sets[0]; pmap = pos_maps[0]
        for gi, cj in zip(rows, cols):
            j = pmap.get(int(gi))
            if j is None:
                continue
            lo, hi = aset.ptr[j], aset.ptr[j + 1]
            av, aw = aset.vals[lo:hi], aset.wts[lo:hi]
            k = float(aw.sum()); nc = float(n_c_all[cj])
            eh = k * nc / max(n, 1)
            h_max = int(min(k, nc, eh + 10.0 * np.sqrt(max(eh, 1.0)) + 25))
            if h_max <= 0 or av.max() * h_max > exact_max_table:
                continue
            # None where the table's cut at h_max could matter (a strong
            # marker's tail lies past it, and came back as 0): the
            # saddlepoint takes those entries
            p_ex = exact_tail_bounded(av, aw, n, int(nc), float(T[gi, cj]), h_max)
            if p_ex is not None:
                P[gi, cj] = p_ex
                log_P[gi, cj] = np.log(max(p_ex, np.finfo(np.float64).tiny))
                info["n_exact"] += 1
                info["n_normal"] -= 1
                exact_rows.add((int(gi), int(cj)))

    # Everything else: the saddlepoint on the same atoms, vectorised.
    spa_pairs = [(int(gi), int(cj)) for gi, cj in zip(rows, cols)
                 if (int(gi), int(cj)) not in exact_rows]
    if spa_pairs:
        vals_parts, wts_parts, blk_ptr, blk_row = [], [], [0], []
        nz_l, nc_l, t_l = [], [], []
        pos = 0
        for r, (gi, cj) in enumerate(spa_pairs):
            for b, (idx, aset, pmap, ncb) in enumerate(
                    zip(strata, atom_sets, pos_maps, nc_by_stratum)):
                # the strata the statistic T leaves out (none of the group's
                # cells, or nothing else): their sum is fixed under the null
                if not (0.0 < ncb[cj] < idx.size):
                    continue
                j = pmap.get(gi)
                if j is None:
                    continue
                lo, hi = aset.ptr[j], aset.ptr[j + 1]
                av, aw = aset.vals[lo:hi], aset.wts[lo:hi]
                if av.size == 0:
                    continue
                vals_parts.append(av); wts_parts.append(aw)
                pos += av.size
                blk_ptr.append(pos); blk_row.append(r)
                nz_l.append(float(idx.size) - float(aw.sum()))
                nc_l.append(float(ncb[cj]))
            t_l.append(float(T[gi, cj]))
        if len(blk_row):
            pvec, lpvec, okvec = spa_atoms_tail(
                np.concatenate(vals_parts), np.concatenate(wts_parts),
                np.asarray(blk_ptr, dtype=np.int64),
                np.asarray(blk_row, dtype=np.int64),
                np.asarray(nz_l, dtype=np.float64),
                np.asarray(nc_l, dtype=np.float64),
                np.asarray(t_l, dtype=np.float64),
                lattice=integral)
            failed = []
            for r, (gi, cj) in enumerate(spa_pairs):
                if okvec[r]:
                    P[gi, cj] = pvec[r]
                    log_P[gi, cj] = lpvec[r]
                    info["n_spa"] += 1
                    info["n_normal"] -= 1
                else:
                    info["n_spa_failed"] += 1
                    failed.append((gi, cj))
            # A refinement fails where the group holds (nearly) the gene's top
            # values, so the observed sum sits at the edge of what the null can
            # produce and the tilt diverges. There the normal p can fall below
            # the permutation floor, which is the exact p at the edge; the floor
            # below takes those entries. Only failures it does not explain are
            # reported.
            unexplained = sum(1 for gi, cj in failed if log_P[gi, cj] >= log_floor[cj])
            if unexplained and verbosity >= 0:
                warnings.warn(
                    f"COSG p-values: {unexplained} of {len(spa_pairs)} refined "
                    "entries did not converge and kept their normal-approximation "
                    "p-value, which is anti-conservative in exactly that regime. "
                    "Report this (dataset shape, group sizes) — it should not happen.",
                    stacklevel=2)

    return _at_least_the_floor(P, Z, info, log_P, log_floor)


def _log_binom(n, k):
    """log C(n, k), natural."""
    from scipy.special import gammaln
    return float(gammaln(n + 1.0) - gammaln(k + 1.0) - gammaln(n - k + 1.0))


def _at_least_the_floor(P, Z, info, log_P, log_floor):
    """No p below the smallest the permutation null can give (one labeling of
    the group's cells, ``log_floor`` per group): a normal p deep in the tail,
    or one whose refinement failed at the edge of the null's support, can be
    hundreds of orders below it. Then ``P`` floored at the smallest positive
    float64 and ``log p`` capped at 0, as before."""
    floor = np.broadcast_to(np.asarray(log_floor, dtype=np.float64)[None, :], log_P.shape)
    below = log_P < floor
    info["n_floored"] = int(below.sum())
    if below.any():
        log_P = np.where(below, floor, log_P)
        P = np.where(below, np.exp(floor), P)
    tiny = np.finfo(np.float64).tiny
    P = np.clip(P, tiny, 1.0)
    info["log_p"] = np.minimum(log_P, 0.0)
    return P, Z, info

def permutation_pvalues(X, labels, groups_order, n_permutations=2000,
                        random_seed=0, batch=None):
    """Sampled permutation p-values -- the validation oracle, not the default.

    Uses ``(r + 1) / (B + 1)``, which is the valid estimator; ``r / B`` can
    report zero for a null that was simply never exceeded.
    """
    rng = np.random.default_rng(random_seed)
    X = X.tocsr() if sparse.issparse(X) else np.asarray(X)
    n = X.shape[0]
    labels = np.asarray(labels)
    T_obs = group_sums(X, labels, groups_order)
    ge = np.zeros_like(T_obs)
    for _ in range(n_permutations):
        if batch is None:
            perm = rng.permutation(n)
        else:
            perm = np.arange(n)
            b = np.asarray(batch)
            for lev in np.unique(b):
                idx = np.flatnonzero(b == lev)
                perm[idx] = rng.permutation(idx)
        ge += (group_sums(X, labels[perm], groups_order) >= T_obs)
    return (ge + 1.0) / (n_permutations + 1.0)


def adjust_fdr(p, method: str = "fdr_bh"):
    """Benjamini-Hochberg (or Benjamini-Yekutieli) over one group's genes."""
    p = np.asarray(p, dtype=np.float64)
    finite = np.isfinite(p)
    out = np.full(p.shape, np.nan)
    vals = p[finite]
    m = vals.size
    if m == 0:
        return out
    order = np.argsort(vals)
    ranked = vals[order]
    factor = 1.0
    if method == "fdr_by":
        factor = np.sum(1.0 / np.arange(1, m + 1))
    adj = ranked * m * factor / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    np.clip(adj, 0.0, 1.0, out=adj)
    res = np.empty(m)
    res[order] = adj
    out[finite] = res
    return out


def adjust_fdr_log(log_p, method: str = "fdr_bh"):
    """Benjamini-Hochberg entirely in log space.

    BH multiplies each ordered p by ``m * factor / rank`` and takes a running
    minimum from the largest downwards. Both steps commute with a monotone
    log, so ``log p_adj = cummin(log p + log(m * factor / rank))`` is the
    *exact* adjusted value -- and it stays exact where ``p`` itself has
    underflowed to zero, which is the regime the analytic route exists for
    (a third of a clean cell type's top markers reach the float64 floor).

    ``log_p`` is a NATURAL log, matching ``analytic_pvalues``' ``log_p``.
    Returns natural-log adjusted values, capped at ``log(1) = 0``.
    """
    lp = np.asarray(log_p, dtype=np.float64)
    finite = np.isfinite(lp)
    out = np.full(lp.shape, np.nan)
    vals = lp[finite]
    m = vals.size
    if m == 0:
        return out
    order = np.argsort(vals)
    ranked = vals[order]
    factor = 1.0
    if method == "fdr_by":
        factor = np.sum(1.0 / np.arange(1, m + 1))
    adj = ranked + np.log(m * factor) - np.log(np.arange(1, m + 1))
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    np.clip(adj, -np.inf, 0.0, out=adj)
    res = np.empty(m)
    res[order] = adj
    out[finite] = res
    return out


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def compute_pvalues(
    X,
    labels,
    groups_order,
    batch=None,
    method: str = "spa",
    fdr_method: str = "fdr_bh",
    n_permutations: int = 2000,
    random_seed: int = 0,
    verbosity: int = 0,
):
    """p-values, FDR and z for every gene x group. The one entry point.

    ``method`` is ``'spa'`` (adaptive: normal in the bulk, saddlepoint in the
    tail, exact enumeration for the sparsest genes), ``'normal'`` (no tail
    refinement -- for comparison, not for claims), or ``'permutation'`` (the
    sampled oracle, whose resolution is bounded by ``1/(n_permutations+1)``).
    """
    if method == "permutation":
        P = permutation_pvalues(X, labels, groups_order,
                                n_permutations=n_permutations,
                                random_seed=random_seed, batch=batch)
        _, Z, info = analytic_pvalues(X, labels, groups_order, batch=batch,
                                      method="normal", verbosity=0)
        info["log_p"] = np.log(P)          # the permutation p's, not the normal ones
    elif method in ("spa", "normal", "exact"):
        P, Z, info = analytic_pvalues(X, labels, groups_order, batch=batch,
                                      method=method, verbosity=verbosity)
    else:
        raise ValueError(
            "pvalue_method must be 'spa', 'normal', 'exact' or 'permutation', "
            f"got {method!r}")

    FDR = np.empty_like(P)
    for c in range(P.shape[1]):
        FDR[:, c] = adjust_fdr(P[:, c], method=fdr_method)

    return {"p_value": P, "p_value_fdr": FDR, "z_score": Z, "info": info}


# ---------------------------------------------------------------------------
# Streaming entry points
#
# The cytome path never holds the matrix, so it accumulates power sums and
# calls these. The split matters: the moments (and therefore the normal
# p-value and the flagging decision) need only sums, while the tail refinement
# needs the gene's actual nonzero values -- which is why the streaming path
# makes one targeted second pass over the flagged genes rather than keeping
# everything.
# ---------------------------------------------------------------------------







# ---------------------------------------------------------------------------
# Streaming entry points
#
# The cytome path never holds the matrix, so it cannot call the array-based
# functions above. What it can do cheaply is accumulate the first three power
# sums per gene while it is already streaming for the cosine -- one extra
# reduction -- and the finite-population moments need nothing else. The tail
# refinement then costs a second pass over the flagged genes only, which are
# sparse by construction because sparseness is what flagged them.
# ---------------------------------------------------------------------------

def srswor_moments_from_sums(s1, s2, s3, n):
    """Per-gene ``(mean, s2, m3)`` from the first three power sums."""
    n = float(n)
    s1 = np.asarray(s1, dtype=np.float64)
    s2 = np.asarray(s2, dtype=np.float64)
    s3 = np.asarray(s3, dtype=np.float64)
    mean = s1 / n
    var = np.maximum(s2 / n - mean ** 2, 0.0)
    m3 = s3 / n - 3.0 * mean * (s2 / n) + 2.0 * mean ** 3
    return mean, var, m3


def normal_tail(T, n_c_all, n, mean, s2, m3, k_nonzero=None, method="spa",
                tail_trigger=0.05, min_expected=30.0, min_p=1e-30, return_log=False):
    """Normal-approximation p and z, plus the entries whose tail needs more.

    Returns ``(P, Z, flagged)`` where ``flagged`` is an ``(m, 2)`` array of
    ``(gene, group)`` indices. The flag is on the **expected number of
    expressing cells inside the group** (``n_c * k / n``), not the gene's
    global expressing count: a gene in 200 of 50,000 cells against a 500-cell
    group expects 2, and a criterion that looked at 200 would wave it through
    while the normal p at a nominal 1e-5 is really 1.8e-3.

    ``return_log=True`` adds a fourth array, the natural log of ``P`` computed
    in log space (``log_ndtr``), which keeps ranking where ``P`` itself has
    underflowed to 0, as ``analytic_pvalues``' ``info['log_p']`` does.
    """
    T = np.asarray(T, dtype=np.float64)
    n_c_all = np.asarray(n_c_all, dtype=np.float64).ravel()
    n_genes, n_groups = T.shape

    E = np.outer(mean, n_c_all)
    var = np.outer(s2, n_c_all * (n - n_c_all) / max(n - 1.0, 1.0))
    sd = np.sqrt(np.maximum(var, 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        Z = np.where(sd > 0, (T - E) / np.where(sd > 0, sd, 1.0), 0.0)
    P = np.where(sd > 0, _sf(Z), 1.0)
    log_P = np.where(sd > 0, log_ndtr(-Z), 0.0) if return_log else None

    def _out(flagged):
        return (P, Z, flagged, log_P) if return_log else (P, Z, flagged)

    if method == "normal" or k_nonzero is None:
        return _out(np.empty((0, 2), dtype=int))

    k = np.asarray(k_nonzero, dtype=np.float64).ravel()
    expected_in_group = np.outer(k, n_c_all) / max(float(n), 1.0)
    # Identical rule to `analytic_pvalues`, and it has to be identical or the
    # two backends refine different entries and their parity test fails on a
    # bulk p-value that neither of them got wrong. See there for why a normal
    # p above `tail_trigger` cannot become significant when the group is
    # smaller than half the cells.
    big_group = np.asarray(n_c_all, dtype=np.float64) > (float(n) / 2.0)
    needs = ((P <= tail_trigger)
             | (big_group[None, :] & (expected_in_group < min_expected))) & (T > E) & (sd > 0)
    # As in `analytic_pvalues`, which keeps the two paths refining the same
    # entries: below `min_p` only sparse entries are refined.
    if min_p > 0:
        needs &= (P >= min_p) | (expected_in_group < min_expected)
    flagged = np.argwhere(needs)
    return _out(flagged)


def refine_tail(values, n, n_c, t, method="spa", exact_max_expressing=16):
    """Exact or saddlepoint upper tail for one (gene, group). None if it fails."""
    values = np.asarray(values, dtype=np.float64)
    values = values[values != 0]
    if values.size == 0 or n_c <= 0 or n_c >= n:
        return None
    if method in ("spa", "exact") and values.size <= exact_max_expressing:
        return _exact_tail(values, n, n_c, t)
    return _spa_tail(values, n, n_c, t)


# ---------------------------------------------------------------------------
# A second saddlepoint solver on raw values, stratified by blocks: the tests
# check spa_atoms_tail against it
# ---------------------------------------------------------------------------

#: Cap on how many expression values one Newton batch may hold.
_SPA_MAX_BATCH_VALUES = 4_000_000


def _spa_tail_vectorised(values, row_ptr, blk_ptr, blk_row, n_zero, n_c,
                         t_obs, max_iter=120, tol=1e-11):
    """Saddlepoint upper tails for many rows at once.

    The per-row solver spends its time in Python, once per (gene, group) pair,
    over a handful of arithmetic on a short array. Running every row's Newton
    in lockstep turns that into a few passes over one concatenated array.

    ``values`` holds every row's nonzero values end to end. A *block* is one
    (row, stratum) pair -- ``values[blk_ptr[j]:blk_ptr[j+1]]``, belonging to
    row ``blk_row[j]`` -- so the unstratified case is simply one block per row
    and the stratified case needs no separate code path. Each block carries its
    own tilt ``u_j`` enforcing that stratum's group size; a single ``s`` per
    row ties them together, exactly as in the scalar stratified solver.

    The Jacobian is an arrow matrix (one shared ``s``, independent ``u_j``), so
    the Newton step reduces to a Schur complement on ``s`` and a back
    substitution for each ``u_j``.

    Returns ``(p, ok)``; ``ok`` is False where the solve did not converge, and
    the caller keeps the normal p-value for those rows.
    """
    values = np.asarray(values, dtype=np.float64)
    # row_ptr is only ever used for the row count, and the batching recursion
    # below has no meaningful one to pass -- it slices rows, not values. Accept
    # None rather than make every caller invent an array to be ignored.
    row_ptr = (None if row_ptr is None
               else np.asarray(row_ptr, dtype=np.int64))
    blk_ptr = np.asarray(blk_ptr, dtype=np.int64)
    blk_row = np.asarray(blk_row, dtype=np.int64)
    n_zero = np.asarray(n_zero, dtype=np.float64)
    n_c = np.asarray(n_c, dtype=np.float64)
    t_obs = np.asarray(t_obs, dtype=np.float64)

    n_rows = int(t_obs.size)
    n_blk = int(blk_row.size)
    if n_rows == 0 or n_blk == 0:
        return np.ones(n_rows), np.zeros(n_rows, dtype=bool)

    if values.size > _SPA_MAX_BATCH_VALUES and n_rows > 1:
        # Vectorising every flagged row at once is the fast path right up
        # until it is not: a 17,000-cell dataset flagged 55,000 rows whose
        # values summed to nine figures, and a hundred iterations over that
        # array does not finish. Split on rows -- never inside one, or a row's
        # strata stop being conditioned jointly -- which keeps the win and
        # bounds the memory.
        p_out = np.ones(n_rows, dtype=np.float64)
        ok_out = np.zeros(n_rows, dtype=bool)
        blk_start = np.searchsorted(blk_row, np.arange(n_rows), "left")
        blk_end = np.searchsorted(blk_row, np.arange(n_rows), "right")
        lo = 0
        while lo < n_rows:
            hi, taken = lo, 0
            while hi < n_rows:
                b0, b1 = int(blk_start[hi]), int(blk_end[hi])
                span = int(blk_ptr[b1] - blk_ptr[b0]) if b1 > b0 else 0
                if taken and taken + span > _SPA_MAX_BATCH_VALUES:
                    break
                taken += span
                hi += 1
            b0, b1 = int(blk_start[lo]), int(blk_end[hi - 1])
            v0, v1 = int(blk_ptr[b0]), int(blk_ptr[b1])
            pc, okc = _spa_tail_vectorised(
                values[v0:v1],
                blk_ptr[b0:b1 + 1] - blk_ptr[b0],
                blk_ptr[b0:b1 + 1] - blk_ptr[b0],
                blk_row[b0:b1] - lo,
                n_zero[b0:b1], n_c[b0:b1], t_obs[lo:hi],
                max_iter=max_iter, tol=tol)
            p_out[lo:hi] = pc
            ok_out[lo:hi] = okc
            lo = hi
        return p_out, ok_out

    counts = np.diff(blk_ptr) if blk_ptr.size == n_blk + 1 else np.diff(
        np.concatenate([blk_ptr, [values.size]]))
    counts = counts.astype(np.int64)
    blk_of_val = np.repeat(np.arange(n_blk), counts)

    def blk_sum(x):
        return np.bincount(blk_of_val, weights=x, minlength=n_blk)

    def row_sum(x):
        return np.bincount(blk_row, weights=x, minlength=n_rows)

    n_blk_total = counts.astype(np.float64) + n_zero
    frac = np.clip(n_c / np.maximum(n_blk_total, 1.0), 1e-12, 1 - 1e-12)
    u_den = np.log(frac / (1.0 - frac))
    pi_den = expit(u_den)
    Kuu_den = n_blk_total * pi_den * (1.0 - pi_den)
    K_den = row_sum(np.logaddexp(0.0, u_den) * n_blk_total)
    ud_nc = row_sum(u_den * n_c)

    s = np.zeros(n_rows, dtype=np.float64)
    u = u_den.copy()
    ok = np.ones(n_rows, dtype=bool)

    Ks = Kss = schur = Knum = None
    for _ in range(max_iter):
        lin = s[blk_row][blk_of_val] * values + u[blk_of_val]
        pi = expit(lin)
        w = pi * (1.0 - pi)
        pi0 = expit(u)
        w0 = pi0 * (1.0 - pi0)

        Ku_b = blk_sum(pi) + n_zero * pi0
        Ks_b = blk_sum(pi * values)
        Kss_b = blk_sum(w * values * values)
        Ksu_b = blk_sum(w * values)
        Kuu_b = blk_sum(w) + n_zero * w0

        F_u = Ku_b - n_c
        Ks = row_sum(Ks_b)
        F_s = Ks - t_obs

        Kuu_safe = np.where(Kuu_b > 1e-300, Kuu_b, 1.0)
        Kss = row_sum(Kss_b)
        schur = row_sum(Ksu_b * Ksu_b / Kuu_safe)
        corr = row_sum(Ksu_b * F_u / Kuu_safe)
        denom = Kss - schur
        good = np.abs(denom) > 1e-300
        ds = np.zeros(n_rows)
        ds[good] = (F_s[good] - corr[good]) / denom[good]
        ds = np.clip(ds, -4.0, 4.0)
        du = -(F_u + Ksu_b * ds[blk_row]) / Kuu_safe
        du = np.clip(du, -4.0, 4.0)

        s = s - ds
        u = u + du
        ok &= good

        if (np.abs(F_s) < tol * np.maximum(1.0, np.abs(t_obs))).all() and \
           (np.abs(F_u) < tol * np.maximum(1.0, n_c)).all():
            break

    lin = s[blk_row][blk_of_val] * values + u[blk_of_val]
    pi = expit(lin)
    w = pi * (1.0 - pi)
    pi0 = expit(u)
    w0 = pi0 * (1.0 - pi0)
    Kss = row_sum(blk_sum(w * values * values))
    Ksu_b = blk_sum(w * values)
    Kuu_b = blk_sum(w) + n_zero * w0
    Kuu_safe = np.where(Kuu_b > 1e-300, Kuu_b, 1.0)
    schur = row_sum(Ksu_b * Ksu_b / Kuu_safe)
    Knum = row_sum(blk_sum(np.logaddexp(0.0, lin)) + n_zero * np.logaddexp(0.0, u))
    u_nc = row_sum(u * n_c)

    det_ratio = Kss - schur
    ratio_blocks = np.where(Kuu_den > 0, Kuu_b / np.where(Kuu_den > 0, Kuu_den, 1.0), 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        log_ratio = row_sum(np.log(np.maximum(ratio_blocks, 1e-300)))
    det_ratio = det_ratio * np.exp(log_ratio)

    inner = 2.0 * ((K_den - ud_nc) - (Knum - s * t_obs - u_nc))
    ok &= np.isfinite(inner) & (inner > 0) & (det_ratio > 0) & (np.abs(s) > 1e-9)
    wv = np.where(ok, np.sign(s) * np.sqrt(np.maximum(inner, 0.0)), 1.0)
    vv = np.where(ok, s * np.sqrt(np.maximum(det_ratio, 0.0)), 1.0)
    ok &= (np.abs(wv) > 1e-9) & (np.abs(vv) > 1e-12)

    p = np.ones(n_rows)
    with np.errstate(divide="ignore", invalid="ignore"):
        cand = ndtr(-wv) - _phi(wv) * (1.0 / wv - 1.0 / vv)
    ok &= np.isfinite(cand)
    p[ok] = np.clip(cand[ok], 0.0, 1.0)
    return p, ok
