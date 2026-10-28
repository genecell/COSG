"""One rule for a group's top features, on every path.

COSG reports each group's ``n_genes_user`` best features. The paths used to
choose them differently: ``np.argpartition`` (the AnnData paths and the
cytome's unbatched path), a stable sort of the candidates in the order a pass
produced them (the cytome's row-major feature batching), and a sort by score
then index (the column-major scan). Among tied scores, which happens on peaks
and tiles with the same few cells open, each picked different features, and the
choice of ``argpartition`` is an implementation detail of numpy, not a rule.

The rule here: by score, highest first, and among equal scores by feature
index, lowest first. A request for more features than there are returns them
all. NaN ranks last (the callers refuse non-finite input before scoring).
"""
from __future__ import annotations

import numpy as np

__all__ = ["top_n", "top_n_of"]


def top_n(scores, n):
    """Indices of the ``n`` best ``scores``, best first, ties by index.

    Costs O(m) for the selection plus O(n log n) for the order, so the full
    column is never sorted when ``n`` is small.
    """
    s = np.asarray(scores)
    m = int(s.shape[0])
    n = int(min(max(int(n), 0), m))
    if n == 0:
        return np.empty(0, dtype=np.intp)
    if s.dtype.kind == "f" and np.isnan(s).any():
        s = np.where(np.isnan(s), -np.inf, s)
    if n < m:
        # the n-th best value; everything above it is in, and of the values
        # equal to it, the lowest indices fill what is left
        kth = s[np.argpartition(-s, n - 1)[n - 1]]
        above = np.flatnonzero(s > kth)
        tied = np.flatnonzero(s == kth)[: n - above.size]
        sel = np.concatenate([above, tied])
    else:
        sel = np.arange(m)
    return sel[np.lexsort((sel, -s[sel]))]


def top_n_of(values, index, n):
    """Positions of the ``n`` best ``(value, index)`` pairs: by value, highest
    first, then by ``index``, lowest first. For merging the candidates that
    feature blocks or passes produced, each chosen by :func:`top_n`."""
    values = np.asarray(values)
    index = np.asarray(index)
    if values.dtype.kind == "f" and np.isnan(values).any():
        values = np.where(np.isnan(values), -np.inf, values)
    order = np.lexsort((index, -values))
    return order[: max(int(n), 0)]
