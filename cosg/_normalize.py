"""Library-size + log1p normalization, implemented natively.

COSG used to borrow this from ``piaso.preprocessing._normalize_log1p``, which
made the *default* ``run_cosg_cytome`` call import PIASO at call time -- and
PIASO requires ``cosg>=1.0.3``, so COSG's headline streaming feature depended
on a package that depends on COSG.

The import was lazy, so it worked. It was still the wrong shape: log1p is
``log1p(counts / cell_depth * scale)``, the most standard normalization in the
field and no part of PIASO's contribution. COSG owning it means the default
path has no PIASO dependency at all.

``infog`` and ``tfidf`` still come from PIASO, and should: those *are* PIASO
methods, and a user asking for them has PIASO in mind.

Numerically identical to the PIASO implementation this replaces --
tests/test_normalize_log1p.py asserts bit-level equality against it when PIASO
is installed.
"""
from __future__ import annotations

import numpy as np
from scipy import sparse as sp


def normalize_chunk_log1p(chunk, cell_depth_chunk, scale_factor):
    """Per-chunk ``log1p(chunk / cell_depth * scale_factor)``.

    Sparse in, sparse out. That matters: COSG chains the result straight into
    a sparse matmul (``chunk.T @ Lam_chunk``) and into a ``chunk.data ** 2``
    accumulator, both of which break on a dense ``memoryview``.

    Parameters
    ----------
    chunk : scipy.sparse matrix or np.ndarray
        ``(n_cells_chunk, n_features)`` raw counts.
    cell_depth_chunk : np.ndarray
        ``(n_cells_chunk,)`` total counts for this chunk's rows. Zeros become
        1.0 -- an empty cell normalizes to all-zero rather than to NaN.
    scale_factor : float
        Target row sum before the log (conventionally 1e4).

    Returns
    -------
    Same kind as the input: CSR if sparse, ndarray if dense.
    """
    cd = np.asarray(cell_depth_chunk, dtype=np.float64)
    cd = np.where(cd == 0, 1.0, cd)
    scale = float(scale_factor)

    if sp.issparse(chunk):
        out = chunk.copy().astype(np.float32)
        if not sp.isspmatrix_csr(out):
            out = out.tocsr()
        # Row-scale through the indptr rather than building a diagonal
        # matrix -- no intermediate the size of the chunk.
        for i in range(out.shape[0]):
            start, end = out.indptr[i], out.indptr[i + 1]
            if end > start:
                out.data[start:end] = out.data[start:end] / cd[i] * scale
        # log1p(0) == 0, so operating on .data leaves the sparsity structure
        # untouched.
        np.log1p(out.data, out=out.data)
        return out

    dense = np.asarray(chunk, dtype=np.float32)
    dense = dense / cd.reshape(-1, 1) * scale
    return np.log1p(dense, out=dense)
