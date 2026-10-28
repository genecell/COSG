"""COSG's scores from per-group sums, and the scale of its normalized scores.

COSG's score of a feature for a group depends on the data only through three
sums: the feature's values summed over the group's cells, its squared values
summed over all cells, and the group's size (plus, for the expression filter,
the number of the group's cells where the feature is nonzero). Every streaming
path accumulates exactly these and then scores them; :func:`score_from_sums`
is that last step, shared by all of them and public, so that sums added up
elsewhere can be scored without the matrix:

- sums from one half of a count split, with clusters from the other half,
  give markers whose p-values stay valid after clustering;
- sums of shards or files add, so COSG over a collection needs no single
  matrix holding it;
- sums update when cells are added.

:func:`_group_scale` is the scale ``iqr_normalize`` divides a group's raw
scores by, before ``log1p``, in every path (``scale_rule`` chooses the rule).
"""
from __future__ import annotations

import numpy as np


def score_from_sums(dot, sq, group_sizes, nnz=None, *, mu=1.0,
                    remove_lowly_expressed=True, expressed_pct=0.1,
                    expressed_min_num_cells_in_target_group=3,
                    min_total_count=0, batch_group_sizes=None,
                    batch_cell_number_threshold=3):
    """COSG's raw scores from per-group sums.

    For feature ``f`` and group ``k`` with ``N_k`` cells, the cosine between
    the feature's values and the group's indicator is
    ``dot[f, k] / (sqrt(sq[f]) * sqrt(N_k))``, and the score is COSG's
    penalty of it against the feature's cosines in all groups, as in
    :func:`cosg.cosg`. The arithmetic is that of the streaming paths of
    :func:`cosg.run_cosg_cytome`, which call this function, so the scores are
    the ones those paths report before any normalization.

    Parameters
    ----------
    dot : array, shape ``(n_features, n_groups)``
        Per feature and group, the feature's values summed over the group's
        cells. With ``batch_group_sizes``, shape ``(n_batches, n_features,
        n_groups)``: the sums within each batch.
    sq : array, shape ``(n_features,)``
        Per feature, its squared values summed over all cells; with batches,
        ``(n_batches, n_features)``.
    group_sizes : array, shape ``(n_groups,)``
        Cells per group.
    nnz : array, shape ``(n_features, n_groups)``, optional
        Per feature and group, the number of the group's cells where the
        feature is nonzero (over all batches). Needed by the expression
        filter and by ``min_total_count``; ``None`` skips both.
    mu : float, default 1.0
        COSG's penalty on features specific to several groups.
    remove_lowly_expressed : bool, default True
        Score -1 where the feature is nonzero in fewer than
        ``max(expressed_pct_k * N_k, expressed_min_num_cells_in_target_group)``
        of the group's cells.
    expressed_pct : float or sequence of float, default 0.1
        The detection floor: one fraction in [0, 1] for every group, or a
        sequence with one fraction per group, in the order of the columns
        (``[0.1, 0.05, 0.1]`` for three groups). The sums carry no group
        labels, so a dict or a Series keyed by group is refused here; pass
        ``series.reindex(groups_order).to_numpy()``. A group's floor changes
        only that group's column.
    expressed_min_num_cells_in_target_group : int, default 3
    min_total_count : int, default 0
        Score -1 in every group where the feature is nonzero in at most this
        many cells in all.
    batch_group_sizes : array, shape ``(n_batches, n_groups)``, optional
        Cells per group within each batch, for batch-corrected scores: the
        cosine is formed within each batch with that batch's own norms, then
        averaged over the batches where the group has at least
        ``batch_cell_number_threshold`` cells, as ``batch_key`` does.
    batch_cell_number_threshold : int, default 3

    Returns
    -------
    numpy.ndarray, shape ``(n_features, n_groups)``
        The raw scores, -1 where a filter removed the feature for the group
        (COSG's convention). Scores are not normalized: a group's best
        features are ``np.argsort(-scores[:, k])``, and
        ``cosg.iqrLogNormalize`` on the scores with -1 set to 0 gives the
        normalized scores of ``iqr_normalize=True``.

    Examples
    --------
    Sums over two halves add up to the sums over all cells:

    >>> scores = cosg.score_from_sums(dot_a + dot_b, sq_a + sq_b, sizes, nnz_a + nnz_b)

    A floor of 5 % for the second of three groups and 10 % for the others:

    >>> scores = cosg.score_from_sums(dot, sq, sizes, nnz, expressed_pct=[0.1, 0.05, 0.1])
    """
    from ._floor import floor_thresholds, floors_for_columns

    n_groups = np.shape(dot)[-1] if np.ndim(dot) else 0
    floors = floors_for_columns(expressed_pct, n_groups)
    thresholds = None
    if remove_lowly_expressed and nnz is not None:
        thresholds = floor_thresholds(group_sizes, floors, expressed_min_num_cells_in_target_group,
                                      skip_empty=True)
    return _score_from_sums(dot, sq, group_sizes, nnz, mu=mu, thresholds=thresholds,
                            min_total_count=min_total_count, batch_group_sizes=batch_group_sizes,
                            batch_cell_number_threshold=batch_cell_number_threshold)


def _score_from_sums(dot, sq, group_sizes, nnz=None, *, mu=1.0, thresholds=None,
                     min_total_count=0, batch_group_sizes=None, batch_cell_number_threshold=3):
    """:func:`score_from_sums` with the expression filter given as its
    thresholds: per group, the fewest cells a feature must be nonzero in
    (``cosg._floor.floor_thresholds``), 0 for a group left unfiltered;
    ``None`` applies no expression filter. The streaming paths call this, so
    that each passes the thresholds of its own semantics."""
    from ._cpu_streaming import _apply_cosg_penalty_cpu

    dot = np.asarray(dot)
    sq = np.asarray(sq)
    if batch_group_sizes is None:
        if dot.ndim != 2 or sq.ndim != 1 or sq.shape[0] != dot.shape[0]:
            raise ValueError(
                f"dot must be (n_features, n_groups) and sq (n_features,), got "
                f"{dot.shape} and {sq.shape}; pass batch_group_sizes with per-batch sums")
        group_norms_safe = np.maximum(np.sqrt(np.asarray(group_sizes, dtype=np.float64)), 1e-12)
        gene_norms = np.maximum(np.sqrt(sq), 1e-12)
        cosine = dot / (gene_norms[:, None] * group_norms_safe[None, :])
    else:
        bgs = np.asarray(batch_group_sizes, dtype=np.float64)
        if dot.ndim != 3 or sq.ndim != 2 or bgs.shape != (dot.shape[0], dot.shape[2]):
            raise ValueError(
                f"with batch_group_sizes, dot must be (n_batches, n_features, n_groups), "
                f"sq (n_batches, n_features) and batch_group_sizes (n_batches, n_groups), "
                f"got {dot.shape}, {sq.shape} and {bgs.shape}")
        acc = np.zeros(dot.shape[1:], dtype=np.float64)
        valid = np.zeros(dot.shape[2], dtype=np.float64)
        for b in range(dot.shape[0]):
            ok = (bgs[b] >= batch_cell_number_threshold).astype(np.float64)
            if not ok.any():
                continue
            gn = np.maximum(np.sqrt(sq[b]), 1e-12)
            gp = np.maximum(np.sqrt(bgs[b]), 1e-12)
            acc += (dot[b] / (gn[:, None] * gp[None, :])) * ok[None, :]
            valid += ok
        with np.errstate(invalid="ignore", divide="ignore"):
            acc /= valid[None, :]
        # a group valid in no batch has no evidence: 0
        cosine = np.nan_to_num(acc, nan=0.0).astype(np.float32)
    scores = _apply_cosg_penalty_cpu(cosine, mu)
    del cosine
    if nnz is not None:
        nnz = np.asarray(nnz)
        if min_total_count > 0:
            scores[nnz.sum(axis=1) <= min_total_count, :] = -1
        if thresholds is not None:
            # each group's own threshold, on its own column only; a threshold
            # of 0 (an empty group, or a floor of 0 with a minimum of 0) leaves
            # the column alone
            for k in range(scores.shape[1]):
                if thresholds[k] > 0:
                    scores[nnz[:, k] < thresholds[k], k] = -1
    return scores


_IQR_FALLBACKS = ("positive", "constant")
_SCALE_RULES = ("iqr", "positive_median")
# The quantile of a group's positive scores that scales it when its q_upper
# quantile over all features is zero (``iqr_fallback='positive'``).
_FALLBACK_POSITIVE_QUANTILE = 0.05
# The median rule needs this many positive scores to set a group's own scale; a
# group with fewer takes the median of the other groups' scales (QC doc 52).
_MIN_POSITIVE_FOR_MEDIAN = 20


def _check_iqr_fallback(iqr_fallback, allowed=_IQR_FALLBACKS):
    if iqr_fallback not in allowed:
        raise ValueError(f"iqr_fallback must be one of {allowed}, got {iqr_fallback!r}")


def _check_scale_rule(scale_rule):
    if scale_rule not in _SCALE_RULES:
        raise ValueError(f"scale_rule must be one of {_SCALE_RULES}, got {scale_rule!r}")


def _pooled_positive_median(cols):
    """The median rule's stand-in scale: the median, over the groups with at
    least ``_MIN_POSITIVE_FOR_MEDIAN`` positive scores, of each group's median
    positive score (each group counted once); ``None`` when no group has that
    many. ``cols`` holds one array per group (scores, or positive scores only)."""
    meds = []
    for c in cols:
        c = np.asarray(c)
        pos = c[c > 0]
        if pos.size >= _MIN_POSITIVE_FOR_MEDIAN:
            meds.append(float(np.median(pos)))
    return float(np.median(meds)) if meds else None


def _group_scale(col, q_upper=0.95, q_lower=0.75, iqr_fallback="positive", scale_rule="iqr",
                 pooled_median=None):
    """The scale ``iqr_normalize`` divides a group's raw scores by.

    ``col`` holds the group's scores over all features, filtered ones as 0.

    ``scale_rule='iqr'`` (the default): the spread between its ``q_upper``
    and ``q_lower`` quantiles, zeros included, as ``cosg.iqrLogNormalize``
    computes it. When that spread is zero:

    - the ``q_upper`` quantile, when it is positive (ties at the top);
    - when the group scores above zero at fewer than ``1 - q_upper`` of the
      features, so that the ``q_upper`` quantile is zero too, the 5th
      percentile of the group's positive scores. This is the common case on
      peaks and tiles, where a cell type opens a few percent of them. Before,
      the scale there was the constant 1e-6, which set the group's normalized
      scores by itself and made them incomparable with those of other groups.
    - and, so that the scale has no jump where the ``q_upper`` quantile
      reaches zero: while the ``q_lower`` quantile is zero (a group positive at
      under ``1 - q_lower`` of the features), the spread is floored at that
      5th percentile. As the positive share falls to ``1 - q_upper`` the
      spread falls to the group's smallest positive score, and the fallback
      below it would otherwise jump up from there (a 4.5-fold jump at 5 % on
      log-normal scores, QC doc 50.2); the floor changes only groups positive
      at 5 to 5.3 % of the features.

    ``scale_rule='positive_median'``: the median of the group's positive
    scores, for sparse features (peaks, tiles). There most groups score
    above zero at 5 to 25 % of the features, where the spread over all
    features is set by how many of them score at all (depth, the expression
    filter, the size of the feature set) more than by the group's own
    scores. The median of its positive scores gives every group the same
    footing relative to its own scores, and does not move when features
    scoring zero are added. A group with fewer than
    ``_MIN_POSITIVE_FOR_MEDIAN`` positive scores would be scaled by its own
    few markers, so it takes ``pooled_median`` instead (the median of the
    other groups' scales, :func:`_pooled_positive_median`), when given.

    1e-6 is left only for a group with no positive score, which reports
    nothing. ``iqr_fallback='constant'`` restores the rule of COSG 1.2.0 and
    earlier (1e-6 whenever the ``q_upper`` quantile is zero), to reproduce
    scores computed with them; it applies to ``scale_rule='iqr'`` only.
    """
    col = np.asarray(col)
    if col.size == 0:
        return 1e-6
    if scale_rule == "positive_median":
        pos = col[col > 0]
        if pos.size == 0:
            return 1e-6
        if pos.size < _MIN_POSITIVE_FOR_MEDIAN and pooled_median is not None:
            return float(pooled_median)
        return float(np.median(pos))
    q_up = float(np.quantile(col, q_upper))
    q_lo = float(np.quantile(col, q_lower))
    if q_up - q_lo > 0:
        if q_lo == 0 and iqr_fallback == "positive":
            pos = col[col > 0]
            return max(q_up - q_lo, float(np.quantile(pos, _FALLBACK_POSITIVE_QUANTILE)))
        return q_up - q_lo
    if q_up > 0:
        return q_up
    if iqr_fallback == "constant":
        return 1e-6
    pos = col[col > 0]
    if pos.size == 0:
        return 1e-6
    return float(np.quantile(pos, _FALLBACK_POSITIVE_QUANTILE))
