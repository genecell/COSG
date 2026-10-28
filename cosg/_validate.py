"""Checks on the input that every path shares.

A NaN in the matrix made its gene the top marker of every group whose floor
kept it, with the score NaN (NaN sorts first), on the default path, while the
legacy path raised; on the cytome paths it went through unnoticed. COSG's scores
are defined for finite values only, so a non-finite value is refused:

- on an in-memory matrix, before anything is computed (one pass over the
  stored values, in blocks, so no copy of the matrix is made);
- on a store, from the per-feature sums of squares the streaming pass forms
  anyway: a NaN or an infinity in a feature makes its sum non-finite, so the
  check costs nothing and still names the features.
"""
from __future__ import annotations

import numpy as np
from scipy import sparse

_BLOCK = 1 << 24


def _shown(names, idx, limit=5):
    idx = np.asarray(idx)[:limit]
    if names is None:
        return [int(i) for i in idx]
    names = np.asarray(names, dtype=object)
    return [str(names[int(i)]) for i in idx]


def check_finite_matrix(X, feature_names=None, what="the matrix"):
    """Raise ValueError naming the features of ``X`` (cells x features) that
    hold NaN or an infinity."""
    if sparse.issparse(X):
        data = X.data
        if data.dtype.kind not in "fc":
            return
        bad_pos = []
        for lo in range(0, data.size, _BLOCK):
            blk = data[lo:lo + _BLOCK]
            if not np.isfinite(blk).all():
                bad_pos.append(lo + np.flatnonzero(~np.isfinite(blk)))
        if not bad_pos:
            return
        pos = np.concatenate(bad_pos)
        if sparse.isspmatrix_csr(X) or getattr(X, "format", None) == "csr":
            feats = np.unique(X.indices[pos])
        elif sparse.isspmatrix_csc(X) or getattr(X, "format", None) == "csc":
            feats = np.unique(np.searchsorted(X.indptr, pos, side="right") - 1)
        else:
            Xc = X.tocsr()
            feats = np.flatnonzero(~np.isfinite(np.asarray(abs(Xc).sum(axis=0)).ravel()))
    else:
        A = np.asarray(X)
        if A.dtype.kind not in "fc":
            return
        flat = A.reshape(A.shape[0], -1) if A.ndim == 2 else A.reshape(-1, 1)
        bad = np.zeros(flat.shape[1], dtype=bool)
        step = max(1, _BLOCK // max(flat.shape[1], 1))
        for lo in range(0, flat.shape[0], step):
            blk = flat[lo:lo + step]
            fin = np.isfinite(blk)
            if not fin.all():
                bad |= ~fin.all(axis=0)
        if not bad.any():
            return
        feats = np.flatnonzero(bad)
    raise ValueError(
        f"{what} holds NaN or infinite values in {len(feats)} feature(s), e.g. "
        f"{_shown(feature_names, feats)}. COSG's scores are defined for finite "
        f"values; replace them (np.nan_to_num) or leave those features out.")


def check_finite_sums(sq, feature_names=None, what="the matrix", offset=0, lookup=None):
    """Raise ValueError when a feature's sum of squares is not finite, which is
    what a NaN or an infinity among its values makes it.

    ``sq`` is ``(n_features,)`` or ``(n_batches, n_features)``, for the features
    ``offset`` onwards. ``feature_names`` (aligned with ``sq``) or
    ``lookup(global_indices)`` names them in the message."""
    s = np.asarray(sq)
    ok = np.isfinite(s.reshape(-1, s.shape[-1])).all(axis=0) if s.ndim > 1 else np.isfinite(s)
    if ok.all():
        return
    local = np.flatnonzero(~ok)
    if feature_names is not None:
        shown = _shown(feature_names, local)
    elif lookup is not None:
        try:
            shown = [str(x) for x in lookup(local[:5] + int(offset))]
        except Exception:
            shown = [int(i) for i in local[:5] + int(offset)]
    else:
        shown = [int(i) for i in local[:5] + int(offset)]
    raise ValueError(
        f"{what} holds NaN or infinite values in {len(local)} feature(s), e.g. "
        f"{shown}. COSG's scores are defined for finite values; replace them or "
        f"leave those features out.")
