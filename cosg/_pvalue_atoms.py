"""A gene's value distribution, compressed once, and the solver that reads it.

Why a summary at all
--------------------
The conditional saddlepoint reads a gene's expression only through its
cumulant generating function,

    K(s, u) = sum_i log(1 + exp(s x_i + u)),

which is a sum over the *multiset* of values -- two cells with the same value
contribute the same term twice. So the k individual values can be replaced by
atoms ``(value, multiplicity)`` without changing what the solver computes, and
the cost of a row stops scaling with how many cells express the gene. On a
10,000-cell dataset the flagged rows carried 381 million values; as atoms they
carry a few tens each.

The grid
--------
One fixed log-spaced grid, the same for every gene, every dataset and both
backends. Fixed rather than data-derived for three reasons: a streamed pass
cannot know the range before it starts, two backends must agree bit for bit,
and a grid that moves with the data would make p-values depend on which cells
were in the file. Each occupied bin stores its **count and its sum**, so the
atom's value is the bin's own mean -- the first moment is preserved exactly
within every bin and the residual is second order in the bin width. Measured
against the uncompressed solver at a deep-tail point, the compressed tail
agrees to 0.1% or better on normalised values and exactly on counts.

The oracle
----------
For integer data the conditional permutation law of ``T`` is not merely
approximable but computable, by polynomial arithmetic over the same atoms
(``exact_tail_atoms``; Pagano & Tritchler, 1983). That is what the
saddlepoint is measured against in ``tests/test_pvalue_atoms.py`` -- and the
exact tail is itself checked against brute-force enumeration of every subset,
so the chain of evidence bottoms out in counting rather than in another
approximation. Anyone changing the solver should check against it rather than
inventing a fixture. It is reachable as ``pvalue_method='exact'`` and is not
the default: it agrees with the saddlepoint to a few percent at 100-1000x the
cost, and it does not exist for non-integer values, since the construction
indexes its state by the sum and so needs the sum to lie on a lattice.

The lattice
-----------
When the values are integers the sum ``T`` lives on a lattice and
``P(T >= t)`` is a step function, which the continuous Lugannani-Rice formula
approximates from the wrong side of the step -- 7 to 20% conservative in the
sparse regime, measured against an exact enumeration. Evaluating the tail at
``t - 1/2`` (Daniels, 1987) removes it: the same rows come back at 0.94 to
1.00 times the exact value. Integrality is detected from the data, not
declared by the caller.
"""
from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.special import expit, log_expit, ndtr, log_ndtr

__all__ = [
    "exact_tail_atoms",
    "GRID_VMIN", "GRID_RATIO", "GRID_NBINS",
    "value_bin_index", "atoms_from_columns", "AtomSet",
    "spa_atoms_tail", "values_are_integral",
]

#: Smallest positive value the grid resolves; anything below shares bin 0.
GRID_VMIN = 1e-6
#: Ratio between consecutive bin edges. 1.03 keeps the compressed tail within
#: 0.1% of the uncompressed one, an order of magnitude below the saddlepoint's
#: own error.
GRID_RATIO = 1.03
#: Number of bins. 1,200 bins at ratio 1.03 span 1e-6 to ~5e9.
GRID_NBINS = 1200

_LOG_RATIO = np.log(GRID_RATIO)
_SQRT2PI = np.sqrt(2.0 * np.pi)


def _phi(z):
    return np.exp(-0.5 * np.asarray(z, dtype=np.float64) ** 2) / _SQRT2PI


def values_are_integral(values, sample: int = 1_000_000) -> bool:
    """Whether the stored values are whole numbers (raw or split counts).

    Sampled above ``sample`` values: a matrix that is integral is integral
    everywhere, and one that is not shows a fractional value almost at once.
    """
    v = np.asarray(values)
    if v.size == 0:
        return True
    if v.size > sample:
        step = max(1, v.size // sample)
        v = v[::step]
    v = v.astype(np.float64, copy=False)
    return bool(np.all(v == np.rint(v)))


def value_bin_index(values) -> np.ndarray:
    """Bin index on the shared grid for each positive value."""
    v = np.asarray(values, dtype=np.float64)
    idx = np.floor(np.log(np.maximum(v, GRID_VMIN) / GRID_VMIN) / _LOG_RATIO)
    return np.clip(idx, 0, GRID_NBINS - 1).astype(np.int64)


class AtomSet:
    """Per-gene atoms in CSR form: ``vals``/``wts`` sliced by ``ptr``.

    ``gene_ids[j]`` is the feature index of row ``j``; ``k[j]`` is how many
    cells that gene expresses (the atoms' total weight).
    """

    __slots__ = ("vals", "wts", "ptr", "gene_ids", "k")

    def __init__(self, vals, wts, ptr, gene_ids, k):
        self.vals = vals
        self.wts = wts
        self.ptr = ptr
        self.gene_ids = np.asarray(gene_ids, dtype=np.int64)
        self.k = np.asarray(k, dtype=np.float64)

    def __len__(self) -> int:
        return int(self.gene_ids.size)

    def index_of(self) -> dict:
        """feature index -> row in this set."""
        return {int(g): j for j, g in enumerate(self.gene_ids)}

    def n_atoms(self) -> np.ndarray:
        return np.diff(self.ptr)


def _atoms_from_dense_table(counts, sums, gene_ids) -> AtomSet:
    """Ragged atoms from a dense ``(n_genes, GRID_NBINS)`` count/sum table."""
    rows, cols = np.nonzero(counts)
    order = np.lexsort((cols, rows))
    rows, cols = rows[order], cols[order]
    w = counts[rows, cols].astype(np.float64)
    v = sums[rows, cols] / w
    ptr = np.zeros(counts.shape[0] + 1, dtype=np.int64)
    np.cumsum(np.bincount(rows, minlength=counts.shape[0]), out=ptr[1:])
    return AtomSet(v, w, ptr, gene_ids, counts.sum(axis=1))


def atoms_from_columns(X, gene_ids, chunk_genes: int = 4096) -> AtomSet:
    """Atoms for the given features of a cells x genes matrix.

    Only the requested features are read; the caller flags a small subset, so
    this never touches the whole matrix. Genes are processed in chunks so the
    dense binning table stays bounded regardless of how many were flagged.
    """
    gene_ids = np.asarray(gene_ids, dtype=np.int64)
    if gene_ids.size == 0:
        return AtomSet(np.empty(0), np.empty(0), np.zeros(1, dtype=np.int64),
                       gene_ids, np.empty(0))
    csc = X.tocsc() if sparse.issparse(X) else None
    parts = []
    for lo in range(0, gene_ids.size, chunk_genes):
        block = gene_ids[lo:lo + chunk_genes]
        counts = np.zeros((block.size, GRID_NBINS), dtype=np.int64)
        sums = np.zeros((block.size, GRID_NBINS), dtype=np.float64)
        for j, g in enumerate(block):
            if csc is not None:
                seg = slice(csc.indptr[g], csc.indptr[g + 1])
                v = csc.data[seg].astype(np.float64, copy=False)
            else:
                col = np.asarray(X[:, g], dtype=np.float64).ravel()
                v = col[col != 0]
            v = v[v > 0]
            if v.size == 0:
                continue
            b = value_bin_index(v)
            counts[j] = np.bincount(b, minlength=GRID_NBINS)
            sums[j] = np.bincount(b, weights=v, minlength=GRID_NBINS)
        parts.append(_atoms_from_dense_table(counts, sums, block))
    if len(parts) == 1:
        return parts[0]
    vals = np.concatenate([p.vals for p in parts])
    wts = np.concatenate([p.wts for p in parts])
    ptr = np.concatenate([[0]] + [np.diff(p.ptr) for p in parts])
    return AtomSet(vals, wts, np.cumsum(ptr), gene_ids,
                   np.concatenate([p.k for p in parts]))


class StreamingAtomAccumulator:
    """Per-gene atom counts accumulated over row chunks.

    A streamed pass sees each cell once and cannot keep the matrix, but the
    grid is fixed, so counts and sums per (gene, bin) simply add across chunks.
    That makes the streamed summary identical to the in-memory one rather than
    merely close to it.
    """

    def __init__(self, gene_ids):
        self.gene_ids = np.asarray(gene_ids, dtype=np.int64)
        self._pos = {int(g): j for j, g in enumerate(self.gene_ids)}
        self.counts = np.zeros((self.gene_ids.size, GRID_NBINS), dtype=np.int64)
        self.sums = np.zeros((self.gene_ids.size, GRID_NBINS), dtype=np.float64)

    def add_chunk(self, chunk) -> None:
        """Accumulate one cells x genes CSR chunk (full width)."""
        if not sparse.issparse(chunk):
            chunk = sparse.csr_matrix(chunk)
        chunk = chunk.tocsr()
        keep = np.isin(chunk.indices, self.gene_ids)
        if not keep.any():
            return
        cols = chunk.indices[keep]
        vals = chunk.data[keep].astype(np.float64, copy=False)
        pos_of = np.full(int(self.gene_ids.max()) + 1, -1, dtype=np.int64)
        pos_of[self.gene_ids] = np.arange(self.gene_ids.size)
        rows = pos_of[cols]
        good = (rows >= 0) & (vals > 0)
        rows, vals = rows[good], vals[good]
        key = rows * GRID_NBINS + value_bin_index(vals)
        size = self.gene_ids.size * GRID_NBINS
        self.counts += np.bincount(key, minlength=size).reshape(
            self.gene_ids.size, GRID_NBINS)
        self.sums += np.bincount(key, weights=vals, minlength=size).reshape(
            self.gene_ids.size, GRID_NBINS)

    def finish(self) -> AtomSet:
        return _atoms_from_dense_table(self.counts, self.sums, self.gene_ids)


def _lugannani_rice(w, v, log: bool = False):
    """Upper tail from the saddlepoint's ``(w, v)``, optionally in log space.

    The log form matters because ``p`` saturates: a third of the top markers
    of a clean cell type reach the float64 floor (~1e-308) and every one of
    them then reads as the same number. ``log_ndtr`` has no floor, so
    ``-log10 p`` keeps discriminating where ``p`` cannot.
    """
    w = np.asarray(w, dtype=np.float64)
    v = np.asarray(v, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        corr = _phi(w) * (1.0 / w - 1.0 / v)
        plain = ndtr(-w) - corr
    if not log:
        return plain
    # log(lead - corr) = log(lead) + log1p(-corr/lead). Both `lead` and `corr`
    # underflow well before their *ratio* does — past w ~ 38 each is 0 in
    # double precision while the ratio tends to 1 - w/v — so the ratio is
    # formed in log space and never as a quotient of two zeros. That is what
    # lets -log10 p keep ranking markers a float64 p-value cannot separate.
    with np.errstate(divide="ignore", invalid="ignore"):
        log_lead = log_ndtr(-w)
        log_phi = -0.5 * w * w - 0.5 * np.log(2.0 * np.pi)
        ratio = np.exp(log_phi - log_lead) * (1.0 / w - 1.0 / v)
        safe = np.isfinite(ratio) & (ratio < 1.0 - 1e-15)
        out = np.where(safe, log_lead + np.log1p(-np.where(safe, ratio, 0.0)),
                       np.log(np.clip(plain, np.finfo(np.float64).tiny, 1.0)))
    return out


def spa_atoms_tail(vals, wts, blk_ptr, blk_row, n_zero, n_c, t_obs,
                   lattice: bool = False, max_iter: int = 80, tol: float = 1e-11,
                   max_batch_values: int = 2_000_000):
    """``P(T >= t)`` for many (gene, group) rows at once, from atoms.

    A **block** is one (row, stratum) pair -- ``vals[blk_ptr[j]:blk_ptr[j+1]]``
    with weights ``wts[...]``, belonging to row ``blk_row[j]`` -- so the
    unstratified case is one block per row and stratification needs no second
    code path. One tilt ``s`` per row ties the blocks together; each block
    carries a ``u_j`` enforcing that stratum's group size. The Jacobian is an
    arrow matrix, so a Newton step is a Schur complement on ``s`` and a back
    substitution for each ``u_j``.

    Rows converge at different depths, so converged rows leave the iteration
    (the active set) rather than every row waiting for the slowest. Rows that
    Newton cannot solve go to a bisection on ``s``, which is monotone and
    therefore cannot fail -- the caller never has to fall back to a normal
    p-value that is known to be wrong by orders of magnitude exactly there.

    Returns ``(p, log_p, ok)``. ``ok`` is False only if even bisection failed,
    which the caller must report rather than silently paper over.
    """
    vals = np.asarray(vals, dtype=np.float64)
    wts = np.asarray(wts, dtype=np.float64)
    blk_ptr = np.asarray(blk_ptr, dtype=np.int64)
    blk_row = np.asarray(blk_row, dtype=np.int64)
    n_zero = np.asarray(n_zero, dtype=np.float64)
    n_c = np.asarray(n_c, dtype=np.float64)
    t_obs = np.asarray(t_obs, dtype=np.float64)
    n_rows = int(t_obs.size)
    n_blk = int(blk_row.size)
    if n_rows == 0 or n_blk == 0:
        return np.ones(n_rows), np.zeros(n_rows), np.zeros(n_rows, dtype=bool)

    # Split on rows when the batch is large. Never inside a row: a row's
    # strata are conditioned jointly and slicing them apart changes the null.
    if vals.size > max_batch_values and n_rows > 1:
        p_out = np.ones(n_rows)
        lp_out = np.zeros(n_rows)
        ok_out = np.zeros(n_rows, dtype=bool)
        blk_start = np.searchsorted(blk_row, np.arange(n_rows), "left")
        blk_end = np.searchsorted(blk_row, np.arange(n_rows), "right")
        lo = 0
        while lo < n_rows:
            hi, taken = lo, 0
            while hi < n_rows:
                b0, b1 = int(blk_start[hi]), int(blk_end[hi])
                span = int(blk_ptr[b1] - blk_ptr[b0]) if b1 > b0 else 0
                if taken and taken + span > max_batch_values:
                    break
                taken += span
                hi += 1
            b0, b1 = int(blk_start[lo]), int(blk_end[hi - 1])
            v0, v1 = int(blk_ptr[b0]), int(blk_ptr[b1])
            pc, lpc, okc = spa_atoms_tail(
                vals[v0:v1], wts[v0:v1], blk_ptr[b0:b1 + 1] - blk_ptr[b0],
                blk_row[b0:b1] - lo, n_zero[b0:b1], n_c[b0:b1], t_obs[lo:hi],
                lattice=lattice, max_iter=max_iter, tol=tol,
                max_batch_values=max_batch_values)
            p_out[lo:hi], lp_out[lo:hi], ok_out[lo:hi] = pc, lpc, okc
            lo = hi
        return p_out, lp_out, ok_out

    target = t_obs - 0.5 if lattice else t_obs

    counts = np.diff(blk_ptr).astype(np.int64)
    blk_of_val = np.repeat(np.arange(n_blk), counts)
    wv = wts * vals

    def blk_sum(x):
        return np.bincount(blk_of_val, weights=x, minlength=n_blk)

    def row_sum(x):
        return np.bincount(blk_row, weights=x, minlength=n_rows)

    n_blk_total = blk_sum(wts) + n_zero
    frac = np.clip(n_c / np.maximum(n_blk_total, 1.0), 1e-12, 1 - 1e-12)
    u_den = np.log(frac / (1.0 - frac))
    pi_den = expit(u_den)
    Kuu_den = n_blk_total * pi_den * (1.0 - pi_den)
    K_den = row_sum(-log_expit(-u_den) * n_blk_total)
    ud_nc = row_sum(u_den * n_c)

    def state(s, u):
        """(K_s, K_u per block, K_ss, K_su per block, K_uu per block, K)."""
        lin = s[blk_row][blk_of_val] * vals + u[blk_of_val]
        pi = expit(lin)
        wq = wts * pi * (1.0 - pi)
        pi0 = expit(u)
        w0 = n_zero * pi0 * (1.0 - pi0)
        Ku_b = blk_sum(wts * pi) + n_zero * pi0
        Ks = row_sum(blk_sum(wv * pi))
        Kss = row_sum(blk_sum(wq * vals * vals))
        Ksu_b = blk_sum(wq * vals)
        Kuu_b = blk_sum(wq) + w0
        K = row_sum(blk_sum(-wts * log_expit(-lin)) - n_zero * log_expit(-u))
        return Ks, Ku_b, Kss, Ksu_b, Kuu_b, K

    # Newton, on a shrinking set. Most rows converge in about six iterations
    # and a handful of very sparse ones need dozens; masking the *update* of
    # converged rows still evaluated them, so every row paid the slowest row's
    # iteration count -- 80 passes over 3.5 million atom entries where 6 were
    # needed. Converged rows are physically removed instead.
    s = np.zeros(n_rows)
    u = u_den.copy()
    a_rows = np.arange(n_rows)              # local row -> original row
    a_blocks = np.arange(n_blk)             # local block -> original block
    l_vals, l_wts, l_wv = vals, wts, wv
    l_bov, l_brow = blk_of_val, blk_row
    l_nz, l_nc, l_tgt = n_zero, n_c, target
    l_s, l_u = s.copy(), u.copy()
    n_r, n_b = n_rows, n_blk

    def local_state(ss, uu):
        lin = ss[l_brow][l_bov] * l_vals + uu[l_bov]
        pi = expit(lin)
        wq = l_wts * pi * (1.0 - pi)
        pi0 = expit(uu)
        Ku_b = np.bincount(l_bov, weights=l_wts * pi, minlength=n_b) + l_nz * pi0
        Ks = np.bincount(l_brow, weights=np.bincount(l_bov, weights=l_wv * pi,
                                                     minlength=n_b), minlength=n_r)
        Kss = np.bincount(l_brow, weights=np.bincount(l_bov, weights=wq * l_vals * l_vals,
                                                      minlength=n_b), minlength=n_r)
        Ksu_b = np.bincount(l_bov, weights=wq * l_vals, minlength=n_b)
        Kuu_b = np.bincount(l_bov, weights=wq, minlength=n_b) + l_nz * pi0 * (1.0 - pi0)
        return Ks, Ku_b, Kss, Ksu_b, Kuu_b

    for _ in range(max_iter):
        if n_r == 0:
            break
        Ks, Ku_b, Kss, Ksu_b, Kuu_b = local_state(l_s, l_u)
        F_s = Ks - l_tgt
        F_u = Ku_b - l_nc
        row_res_u = np.bincount(l_brow, weights=np.abs(F_u), minlength=n_r)
        row_nc = np.bincount(l_brow, weights=l_nc, minlength=n_r)
        done = (np.abs(F_s) < tol * np.maximum(1.0, np.abs(l_tgt))) & \
               (row_res_u < tol * np.maximum(1.0, row_nc))
        if done.any():
            s[a_rows[done]] = l_s[done]
            keep_b = ~done[l_brow]
            u[a_blocks[keep_b == False]] = l_u[~keep_b]
            if done.all():
                n_r = 0
                break
            keep_r = ~done
            keep_v = keep_b[l_bov]
            new_row_id = np.cumsum(keep_r) - 1
            new_blk_id = np.cumsum(keep_b) - 1
            l_vals, l_wts, l_wv = l_vals[keep_v], l_wts[keep_v], l_wv[keep_v]
            l_bov = new_blk_id[l_bov[keep_v]]
            l_brow = new_row_id[l_brow[keep_b]]
            l_nz, l_nc = l_nz[keep_b], l_nc[keep_b]
            l_tgt = l_tgt[keep_r]
            l_s, l_u = l_s[keep_r], l_u[keep_b]
            a_rows, a_blocks = a_rows[keep_r], a_blocks[keep_b]
            n_r, n_b = int(keep_r.sum()), int(keep_b.sum())
            if n_r == 0:
                break
            Ks, Ku_b, Kss, Ksu_b, Kuu_b = local_state(l_s, l_u)
            F_s = Ks - l_tgt
            F_u = Ku_b - l_nc
        Kuu_safe = np.where(Kuu_b > 1e-300, Kuu_b, 1.0)
        schur = np.bincount(l_brow, weights=Ksu_b * Ksu_b / Kuu_safe, minlength=n_r)
        corr = np.bincount(l_brow, weights=Ksu_b * F_u / Kuu_safe, minlength=n_r)
        denom = Kss - schur
        good = np.abs(denom) > 1e-300
        ds = np.zeros(n_r)
        ds[good] = (F_s[good] - corr[good]) / denom[good]
        ds = np.clip(ds, -4.0, 4.0)
        du = np.clip(-(F_u + Ksu_b * ds[l_brow]) / Kuu_safe, -4.0, 4.0)
        l_s = l_s - ds
        l_u = l_u + du
    if n_r:
        s[a_rows] = l_s
        u[a_blocks] = l_u

    # Rows Newton could not finish go to a bisection on `s`, which is monotone
    # along the constrained path and therefore cannot fail. This is what keeps
    # a refined row from silently reverting to its normal-approximation
    # p-value -- wrong by orders of magnitude in exactly this regime, and
    # exactly these rows.
    Ks, _Ku, _Kss, _Ksu, _Kuu, _K = state(s, u)
    unconverged = np.abs(Ks - target) > 1e-6 * np.maximum(1.0, np.abs(target))
    unconverged |= ~np.isfinite(s)
    if unconverged.any():
        # On the unconverged rows only. Bisection is ~300 outer steps with an
        # inner Newton inside each, so running it across the whole batch
        # because a handful of rows need it costs more than the whole Newton
        # phase -- on a 10,000-cell dataset that was most of the wall time.
        keep_b = unconverged[blk_row]
        keep_v = keep_b[blk_of_val]
        new_row = np.cumsum(unconverged) - 1
        new_blk = np.cumsum(keep_b) - 1
        sub_rows = int(unconverged.sum())
        sub_blks = int(keep_b.sum())
        sub_bov = new_blk[blk_of_val[keep_v]]
        sub_brow = new_row[blk_row[keep_b]]

        def sub_blk_sum(x):
            return np.bincount(sub_bov, weights=x, minlength=sub_blks)

        def sub_row_sum(x):
            return np.bincount(sub_brow, weights=x, minlength=sub_rows)

        s_b, u_b = _bisect_s(
            vals[keep_v], wts[keep_v], sub_bov, sub_brow,
            sub_blk_sum, sub_row_sum, n_zero[keep_b], n_c[keep_b],
            target[unconverged], u_den[keep_b],
            np.ones(sub_rows, dtype=bool), sub_rows, sub_blks)
        s = s.copy(); u = u.copy()
        s[unconverged] = s_b
        u[keep_b] = u_b

    Ks, Ku_b, Kss, Ksu_b, Kuu_b, Knum = state(s, u)
    Kuu_safe = np.where(Kuu_b > 1e-300, Kuu_b, 1.0)
    schur = row_sum(Ksu_b * Ksu_b / Kuu_safe)
    u_nc = row_sum(u * n_c)

    det_ratio = Kss - schur
    ratio_blocks = np.where(Kuu_den > 0, Kuu_b / np.where(Kuu_den > 0, Kuu_den, 1.0), 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        det_ratio = det_ratio * np.exp(row_sum(np.log(np.maximum(ratio_blocks, 1e-300))))

    inner = 2.0 * ((K_den - ud_nc) - (Knum - s * target - u_nc))
    ok = (np.abs(Ks - target) < 1e-5 * np.maximum(1.0, np.abs(target)))
    ok &= np.isfinite(inner) & (inner > 0) & (det_ratio > 0) & (np.abs(s) > 1e-9)
    w = np.where(ok, np.sign(s) * np.sqrt(np.maximum(inner, 0.0)), 1.0)
    v = np.where(ok, s * np.sqrt(np.maximum(det_ratio, 0.0)), 1.0)
    ok &= (np.abs(w) > 1e-9) & (np.abs(v) > 1e-12)

    p = np.ones(n_rows)
    log_p = np.zeros(n_rows)
    cand = _lugannani_rice(w, v)
    log_cand = _lugannani_rice(w, v, log=True)
    ok &= np.isfinite(cand)
    p[ok] = np.clip(cand[ok], 0.0, 1.0)
    log_p[ok] = np.minimum(log_cand[ok], 0.0)
    return p, log_p, ok


def _bisect_s(vals, wts, blk_of_val, blk_row, blk_sum, row_sum,
              n_zero, n_c, target, u_den, mask, n_rows, n_blk,
              outer: int = 200, inner: int = 60):
    """Bisection on ``s`` with the per-block ``u`` solved inside.

    ``K_s`` increases in ``s`` along the constrained path, so bracketing and
    halving always converges. Slower than Newton and used only where Newton
    gave up -- which is the deep tail of the sparsest genes, precisely the
    rows whose p-values must not silently revert to the normal approximation.
    """
    def solve_u(s):
        u = u_den.copy()
        for _ in range(inner):
            lin = s[blk_row][blk_of_val] * vals + u[blk_of_val]
            pi = expit(lin)
            pi0 = expit(u)
            f = blk_sum(wts * pi) + n_zero * pi0 - n_c
            d = blk_sum(wts * pi * (1.0 - pi)) + n_zero * pi0 * (1.0 - pi0)
            step = np.clip(f / np.where(d > 1e-300, d, 1.0), -4.0, 4.0)
            u = u - step
            if np.max(np.abs(step)) < 1e-13:
                break
        return u

    def g(s):
        u = solve_u(s)
        lin = s[blk_row][blk_of_val] * vals + u[blk_of_val]
        return row_sum(blk_sum(wts * vals * expit(lin))) - target, u

    lo = np.full(n_rows, -1.0)
    hi = np.full(n_rows, 1.0)
    for _ in range(60):
        glo, _ = g(lo)
        need = mask & (glo > 0)
        if not need.any():
            break
        lo = np.where(need, lo * 2.0, lo)
    for _ in range(60):
        ghi, _ = g(hi)
        need = mask & (ghi < 0)
        if not need.any():
            break
        hi = np.where(need, hi * 2.0, hi)
    u_out = u_den.copy()
    for _ in range(outer):
        mid = 0.5 * (lo + hi)
        val, u_mid = g(mid)
        u_out = np.where(mask[blk_row], u_mid, u_out)
        hi = np.where(mask & (val > 0), mid, hi)
        lo = np.where(mask & (val <= 0), mid, lo)
        if np.max(np.where(mask, hi - lo, 0.0)) < 1e-13:
            break
    s_out = 0.5 * (lo + hi)
    return s_out, u_out


def exact_tail_atoms(vals, wts, n, n_c, t, h_max, pi=None):
    """Exact ``P(T >= t)`` for integer-valued data, from the same atoms.

    Include each of the ``k`` expressing cells with an independent
    Bernoulli(pi) coin. The joint law of (number included ``h``, sum ``T``) is
    the product over atoms of ``(1 - pi + pi*y*z^v)^w`` -- one small
    two-dimensional polynomial. Conditional on ``h`` the included set is a
    uniform ``h``-subset whatever ``pi`` is, so each row of that table,
    normalised, is the exact law of ``T | h``; the ``n - k`` zero cells enter
    only through the hypergeometric law of ``h``. One polynomial per gene
    serves every group (Pagano & Tritchler, JASA 1983).

    ``h_max`` truncates the table where ``P(h)`` is negligible. Returns None
    when the values are not integers or the table would be unaffordable.
    """
    from scipy.stats import hypergeom, binom

    vals = np.asarray(vals, dtype=np.float64)
    wts = np.asarray(wts, dtype=np.float64)
    if vals.size == 0 or not np.all(vals == np.rint(vals)):
        return None
    k = int(round(float(wts.sum())))
    h_max = int(min(h_max, k, n_c))
    if h_max <= 0:
        return None
    v_int = np.rint(vals).astype(np.int64)
    w_int = np.rint(wts).astype(np.int64)
    T_max = int(v_int.max()) * h_max
    if (h_max + 1) * (T_max + 1) > 5e7:
        return None
    pi = min(max(n_c / n, 1e-6), 0.5) if pi is None else pi
    Q = np.zeros((h_max + 1, T_max + 1))
    Q[0, 0] = 1.0
    for v, w in zip(v_int, w_int):
        jmax = min(int(w), h_max)
        coef = binom.pmf(np.arange(jmax + 1), int(w), pi)
        new = coef[0] * Q
        for j in range(1, jmax + 1):
            if coef[j] <= 0.0:
                continue
            new[j:, v * j:] += coef[j] * Q[:h_max + 1 - j, :T_max + 1 - v * j]
        Q = new
    ph = hypergeom.pmf(np.arange(h_max + 1), n, k, int(n_c))
    lo = int(np.ceil(t))
    if lo > T_max:
        tail = np.zeros(h_max + 1)
    else:
        tail = Q[:, max(lo, 0):].sum(axis=1)
    tot = Q.sum(axis=1)
    ok = tot > 0
    return float(np.clip(np.sum(ph[ok] * tail[ok] / tot[ok]), 0.0, 1.0))
