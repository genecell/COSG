"""COSG's detection floor, ``expressed_pct``: the share of a group's cells in
which a feature must be detected for the feature to be one of its markers.

A feature is kept for group ``k`` (``N_k`` cells) when it is nonzero in at
least ``max(expressed_pct_k * N_k, expressed_min_num_cells_in_target_group)``
of the group's cells, and scores -1 for that group otherwise. The fraction is
one number for every group, or one per group:

- a number in [0, 1]; ``None`` means the default, 0.1;
- a dict ``{group label: fraction}``, or a pandas Series indexed by group;
- a sequence (list, tuple, 1-D array) with one fraction per group, in the
  order of the groups in the result (``groups_order``).

With a dict or a Series every group that has cells needs a fraction: a group
left out is an error rather than a silent 0.1, because a floor read from the
data for each type (``cosg_pertype_floor.py`` in the PIASO benchmarks) can be
an order of magnitude below 0.1 on peaks and tiles, where a type left at 0.1
would lose most of its markers without a word. A key that names no label of
``groupby`` is an error too (a typo would otherwise leave its group without
the floor meant for it); a key naming a label that is not scored in this call
(removed by ``groups=`` or ``cell_mask``) is ignored. Keys are matched to the
labels as they are, or by their string form, so ``{"3": 0.05}`` names
cluster ``3``.

The floor removes a feature from one group's markers only: it is applied to
that group's own column of scores, after the cosine and COSG's penalty, which
use every group but not the floor. So the markers a group gets under a
per-group floor are those it gets with its own fraction as the scalar floor.
One coupling remains, not in the floor itself:
``scale_rule='positive_median'`` scales a group with fewer than 20 positive
scores by the median of the other groups' scales, which their floors move.

A floor of 0 has one meaning on every path: the group keeps the minimum of
``expressed_min_num_cells_in_target_group`` cells (3 by default), so a feature
seen in one or two of its cells is never one of its markers. Only
``remove_lowly_expressed=False`` switches the filter off. (Before, a 0 switched
it off altogether, the minimum included, on the streaming paths, and kept the
minimum on the others.)
"""
from __future__ import annotations

import numbers
from collections.abc import Mapping, Sequence, Set

import numpy as np

#: The floor when none is given, and of a group without cells (where no floor
#: has any effect).
DEFAULT_EXPRESSED_PCT = 0.1


def _pandas():
    try:
        import pandas as pd
    except ImportError:  # pragma: no cover - pandas is a dependency
        return None
    return pd


def _form(expressed_pct):
    """``'scalar'``, ``'mapping'`` or ``'sequence'``; TypeError for any other value."""
    if expressed_pct is None:
        return "scalar"
    if isinstance(expressed_pct, np.ndarray):
        if expressed_pct.ndim == 0:
            return "scalar"
        if expressed_pct.ndim == 1:
            return "sequence"
        raise TypeError(
            f"expressed_pct must be a number, a dict or Series keyed by group, or a "
            f"1-D sequence with one fraction per group; got a {expressed_pct.ndim}-D array.")
    if isinstance(expressed_pct, (numbers.Real, np.number)):
        return "scalar"
    pd = _pandas()
    if pd is not None:
        if isinstance(expressed_pct, pd.DataFrame):
            raise TypeError(
                "expressed_pct got a DataFrame; pass one column of it as a Series indexed "
                "by group, e.g. table.set_index('cell_type')['fraction'].")
        if isinstance(expressed_pct, pd.Series):
            return "mapping"
        if isinstance(expressed_pct, pd.Index):
            return "sequence"
    if isinstance(expressed_pct, Mapping):
        return "mapping"
    if isinstance(expressed_pct, Set):
        raise TypeError(
            "expressed_pct got a set, which has no order; pass a dict or a Series keyed "
            "by group, or a list in the order of the groups.")
    if isinstance(expressed_pct, Sequence) and not isinstance(expressed_pct, (str, bytes)):
        return "sequence"
    raise TypeError(
        f"expressed_pct must be a number, a dict or Series keyed by group, or a sequence "
        f"with one fraction per group; got {type(expressed_pct).__name__}.")


def is_per_group(expressed_pct):
    """Whether ``expressed_pct`` gives a fraction per group (dict, Series or sequence)."""
    return _form(expressed_pct) != "scalar"


def _scalar(expressed_pct):
    value = DEFAULT_EXPRESSED_PCT if expressed_pct is None else float(expressed_pct)
    if not 0.0 <= value <= 1.0:          # NaN fails this too
        raise ValueError(f"expressed_pct must be between 0 and 1, got {expressed_pct!r}.")
    return value


def _values(raw):
    try:
        values = np.asarray(list(raw), dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"expressed_pct's fractions must be numbers: {exc}") from None
    if values.ndim != 1:
        raise TypeError("expressed_pct's fractions must be numbers, one per group.")
    return values


def _check_range(labels, values):
    bad = ~((values >= 0.0) & (values <= 1.0))      # NaN is out of range
    if bad.any():
        shown = {labels[i]: float(values[i]) for i in np.flatnonzero(bad)[:10]}
        raise ValueError(
            f"expressed_pct must be between 0 and 1 for every group; out of range: {shown}"
            + (f" and {int(bad.sum()) - 10} more" if bad.sum() > 10 else "") + ".")


def _items(expressed_pct):
    pd = _pandas()
    if pd is not None and isinstance(expressed_pct, pd.Series):
        if expressed_pct.index.has_duplicates:
            dup = expressed_pct.index[expressed_pct.index.duplicated()].unique().tolist()
            raise ValueError(f"expressed_pct's index names a group more than once: {dup[:10]}.")
        return list(expressed_pct.index), _values(expressed_pct.to_numpy())
    keys = list(expressed_pct.keys())
    return keys, _values([expressed_pct[k] for k in keys])


def check_expressed_pct(expressed_pct):
    """Check the form of ``expressed_pct`` and the range of its fractions before
    the groups are known (they are checked against the groups by
    :func:`resolve_expressed_pct`). Raises TypeError or ValueError."""
    form = _form(expressed_pct)
    if form == "scalar":
        _scalar(expressed_pct)
    elif form == "mapping":
        keys, values = _items(expressed_pct)
        _check_range(keys, values)
    else:
        values = _values(expressed_pct)
        _check_range(list(range(len(values))), values)


def resolve_expressed_pct(expressed_pct, groups_order, group_sizes=None, known_labels=None):
    """One floor per group, in the order of ``groups_order``.

    Parameters
    ----------
    expressed_pct
        A number, ``None`` (0.1), a dict or Series keyed by group label, or a
        sequence with one fraction per group in the order of ``groups_order``.
    groups_order : sequence
        The groups scored, in the order of the result.
    group_sizes : sequence, optional
        Cells per group. A group without cells needs no fraction (no floor
        changes its result); it gets the default.
    known_labels : iterable, optional
        Every label of ``groupby``, scored or not. A dict or Series key naming
        one of them that is not in ``groups_order`` is ignored; a key naming
        none of them is an error.

    Returns
    -------
    floors : numpy.ndarray of float64, shape ``(len(groups_order),)``
    per_group : bool
        True when ``expressed_pct`` gave fractions per group (and should be
        recorded per group in the output).
    """
    groups = list(groups_order)
    n = len(groups)
    form = _form(expressed_pct)
    if form == "scalar":
        return np.full(n, _scalar(expressed_pct), dtype=np.float64), False

    if form == "sequence":
        values = _values(expressed_pct)
        if len(values) != n:
            raise ValueError(
                f"expressed_pct has {len(values)} fractions for {n} groups. A sequence gives one "
                f"fraction per group, in the order of the result's groups: {_shown(groups)}. A "
                f"dict or a Series keyed by group does not depend on that order.")
        _check_range(groups, values)
        return values.copy(), True

    keys, values = _items(expressed_pct)
    _check_range(keys, values)
    exact = {}
    by_str = {}
    for i, g in enumerate(groups):
        try:
            exact.setdefault(g, i)
        except TypeError:              # an unhashable label: matched by its string only
            pass
        by_str.setdefault(str(g), []).append(i)
    known = {str(x) for x in known_labels} if known_labels is not None else set()

    floors = np.full(n, np.nan, dtype=np.float64)
    source = [None] * n
    unknown = []
    for key, value in zip(keys, values):
        i = None
        try:
            i = exact.get(key)
        except TypeError:
            pass
        if i is None:
            same = by_str.get(str(key), [])
            if len(same) > 1:
                raise ValueError(
                    f"expressed_pct key {key!r} matches more than one group by its string "
                    f"form: {[groups[j] for j in same]}. Use the labels themselves as keys.")
            if same:
                i = same[0]
        if i is None:
            if str(key) not in known:
                unknown.append(key)
            continue                    # a label of groupby that is not scored in this call
        if source[i] is not None:
            raise ValueError(
                f"expressed_pct gives group {groups[i]!r} two fractions "
                f"(keys {source[i]!r} and {key!r}).")
        floors[i] = value
        source[i] = key
    if unknown:
        raise ValueError(
            f"expressed_pct names {len(unknown)} label(s) that are not groups of groupby: "
            f"{_shown(unknown)}. The groups are: {_shown(groups)}.")

    sizes = (np.ones(n) if group_sizes is None
             else np.asarray(group_sizes, dtype=np.float64).reshape(-1))
    missing = [groups[i] for i in range(n) if np.isnan(floors[i]) and sizes[i] > 0]
    if missing:
        raise ValueError(
            f"expressed_pct has no fraction for {len(missing)} group(s): {_shown(missing)}. "
            f"Give every group one -- for example "
            f"{{**{{g: 0.1 for g in groups}}, **your_fractions}} keeps 0.1 for the rest -- "
            f"or pass one number for all of them.")
    floors[np.isnan(floors)] = DEFAULT_EXPRESSED_PCT     # groups without cells
    return floors, True


def floors_for_columns(expressed_pct, n_groups):
    """``expressed_pct`` as one floor per column, for the functions that see the
    groups only as the columns of their sums: a number (``None`` is 0.1), or a
    sequence of ``n_groups`` fractions in the order of the columns. A dict or a
    Series raises TypeError there, since those functions have no labels to
    match its keys to."""
    form = _form(expressed_pct)
    if form == "scalar":
        return np.full(int(n_groups), _scalar(expressed_pct), dtype=np.float64)
    if form == "mapping":
        raise TypeError(
            "expressed_pct is keyed by group here, but this function sees the groups only "
            "as columns; pass one fraction per column in their order, e.g. "
            "series.reindex(groups_order).to_numpy().")
    values = _values(expressed_pct)
    if len(values) != int(n_groups):
        raise ValueError(
            f"expressed_pct has {len(values)} fractions for {int(n_groups)} groups; give one "
            f"per group, in the order of the columns.")
    _check_range(list(range(len(values))), values)
    return values.copy()


def _shown(items, limit=20):
    items = list(items)
    if len(items) <= limit:
        return repr(items)
    return repr(items[:limit])[:-1] + f", ... ({len(items)} in all)]"


def floor_thresholds(group_sizes, floors, min_cells, *, dtype=np.float64,
                     skip_empty=False):
    """The fewest expressing cells a feature needs in each group:
    ``max(N_k * floors[k], min_cells)``, in ``dtype``.

    The paths of COSG compute this in different precisions and leave different
    groups alone, and the per-group vector reproduces each path's scalar
    arithmetic exactly, so a scalar and the same value for every group give
    identical results:

    - ``dtype=np.float32`` for the GPU paths, which form the product in float32
      on the device (a float32 size times a Python float, which CuPy evaluates
      in float32); the CPU paths use float64.
    - ``skip_empty``: 0 (no filter) for a group without cells, for the paths
      that skip such groups. A threshold of 0 filters nothing, since no count
      is below it.

    A floor of 0 gives ``min_cells`` on every path.
    """
    sizes = np.asarray(group_sizes, dtype=np.float64).reshape(-1)
    fl = np.asarray(floors, dtype=np.float64).reshape(-1)
    if np.dtype(dtype) == np.float32:
        thr = np.maximum(sizes.astype(np.float32) * fl.astype(np.float32), np.float32(min_cells))
    else:
        thr = np.maximum(sizes * fl, np.float64(min_cells))
    if skip_empty:
        thr[sizes <= 0] = 0
    return thr


def tracks_floor(remove_lowly_expressed, floors, min_cells):
    """Whether the expression filter needs the nonzero counts: it is on, and
    some group's threshold can be above 0 (a floor above 0, or a minimum of
    cells above 0). A floor of 0 with a minimum of 0 filters nothing."""
    if not remove_lowly_expressed:
        return False
    return bool(np.any(np.asarray(floors) > 0)) or float(min_cells or 0) > 0


def floor_record(groups_order, floors):
    """``{group: fraction}`` in the order of ``groups_order``."""
    return {g: float(f) for g, f in zip(groups_order, np.asarray(floors, dtype=np.float64))}


def floor_record_array(groups_order, floors):
    """The per-group floors as a one-record structured array with a field per
    group, the form of ``adata.uns['cosg']['names']``: it survives
    ``adata.write`` with any label (a dict keyed by a label containing ``/``
    does not), and ``pd.DataFrame(rec).iloc[0]`` reads it back as a Series."""
    import pandas as pd
    return pd.DataFrame([np.asarray(floors, dtype=np.float64)], columns=list(groups_order)).to_records(
        index=False, column_dtypes="float64")


def describe_floor(floors, per_group, expressed_pct):
    """A short description of the floor for the logs."""
    if not per_group:
        return f"{expressed_pct}"
    fl = np.asarray(floors, dtype=np.float64)
    if not len(fl):
        return "per group (no groups)"
    return (f"per group ({len(fl)} groups; min {fl.min():.4g}, median {np.median(fl):.4g}, "
            f"max {fl.max():.4g})")
