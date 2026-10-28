"""What the default scale of ``iqr_normalize`` depends on (audit part 5,
improvement 4, documented here, not changed).

Under ``scale_rule='iqr'`` (the default) a group's scale is the spread between
the 95th and 75th percentiles of its scores over all features, zeros
included. On sparse features, where a group scores above zero at a few
percent of them, that spread is set by how many features score at all, not by
the group's own scores: as the positive share crosses 5 to 10 % the scale moves
from the 5th percentile of the positive scores to about their median, about
five-fold on log-normal scores, for the same positive scores. The share moves
with the feature universe and the detection floor (QC doc 46.7).
``scale_rule='positive_median'`` takes the median of the positive scores and
does not move. Making it the default for peaks and tiles is a method change
left to the next round; this test pins the behaviour the decision is about.
"""
from __future__ import annotations

import numpy as np

from cosg._from_sums import _group_scale


def _column(n_positive, n_features=10_000, seed=0):
    """The same positive scores, among more or fewer zero-scoring features."""
    rng = np.random.default_rng(seed)
    pos = np.sort(rng.lognormal(mean=-3.0, sigma=1.0, size=n_positive))
    col = np.zeros(n_features)
    col[:n_positive] = pos
    return col, pos


def test_the_iqr_scale_follows_the_positive_share():
    shares = (0.03, 0.05, 0.07, 0.10, 0.15)
    pos = _column(500)[1]                                  # one set of positive scores
    scales = {}
    for share in shares:
        n_features = int(round(len(pos) / share))
        col = np.zeros(n_features)
        col[:len(pos)] = pos
        scales[share] = (_group_scale(col, scale_rule="iqr"), _group_scale(col, scale_rule="positive_median"))
    iqr = np.array([scales[s][0] for s in shares])
    med = np.array([scales[s][1] for s in shares])
    # positive_median: one scale whatever the share
    assert np.allclose(med, med[0])
    assert np.isclose(med[0], np.median(pos))
    # iqr: at or under 5 % the 5th percentile of the positive scores ...
    assert np.isclose(iqr[0], np.quantile(pos, 0.05)) and np.isclose(iqr[1], np.quantile(pos, 0.05))
    # ... rising with the share past it, for the same scores
    assert np.all(np.diff(iqr[1:]) > 0)
    assert iqr[3] / iqr[1] > 3.0, iqr / iqr[1]
