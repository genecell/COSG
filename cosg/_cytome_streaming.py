"""
Cytome-native streaming for COSG marker gene detection.

Reads expression chunks directly from on-disk cytome storage without
materializing the full expression matrix. Peak RAM is bounded by one
chunk plus the O(n_genes × n_groups) accumulators, regardless of
dataset size.

CPU path reuses the accumulation math from _cpu_streaming.py.
GPU path reuses the penalty formula from _gpu.py.
"""

import json
import numpy as np
import time
from scipy import sparse as sp_sparse

from ._floor import (check_expressed_pct, describe_floor, floor_record, floor_thresholds,
                     is_per_group, resolve_expressed_pct, tracks_floor)
from ._from_sums import (_check_iqr_fallback, _check_scale_rule, _group_scale, _pooled_positive_median,
                         _score_from_sums)
from ._select import top_n, top_n_of


def _resolve_floors(expressed_pct, groups_order, group_sizes, all_labels, categories):
    """The detection floor of each group in ``groups_order``. A dict or Series
    key that names a label of the column but no group scored here (one that
    ``cell_mask`` left out) is ignored; one that names no label at all raises."""
    known = None
    if is_per_group(expressed_pct):
        known = set(np.asarray(all_labels, dtype=object).tolist()) | set(categories or [])
    return resolve_expressed_pct(expressed_pct, groups_order, group_sizes=group_sizes,
                                 known_labels=known)


# cytome is an optional dependency, imported inside the functions that use it:
# `import cosg` should neither require it nor pay to load a dataset engine for
# the AnnData path, which is most of COSG's use. Previously this was a
# module-level import guarded by a try/except in __init__.py, which meant a
# missing cytome made run_cosg_cytome silently *disappear* -- the caller got
# `AttributeError: module 'cosg' has no attribute 'run_cosg_cytome'`, which
# reads like a typo rather than a missing package.


def _cytome():
    """Import and return the cytome package."""
    try:
        import cytome
    except ImportError as exc:  # pragma: no cover - only on a broken install
        raise ImportError(
            "Reading a .cytome dataset needs cytome, which is an optional "
            "dependency of COSG.\n"
            "    pip install 'cosg[cytome]'   (or: pip install cytome)\n"
            "Marker detection on an AnnData does not need it."
        ) from exc
    return cytome


# Lazily resolved at first call — cytome's Dataset class isn't re-exported
# at the package root, so we look it up dynamically rather than couple to
# cytome's internal module layout.
def _cytome_dataset_classes():
    """Return a tuple of cytome Dataset class(es) for isinstance checks.

    Cached on the function attribute after the first successful lookup.
    Falls back to a no-op (empty tuple) if cytome internals change —
    callers then route everything as a path string.
    """
    cached = getattr(_cytome_dataset_classes, "_cache", None)
    if cached is not None:
        return cached
    candidates = []
    try:
        from cytome.core.dataset import CytomeDataset
        candidates.append(CytomeDataset)
    except ImportError:
        pass
    # Also accept any future ``cytome.Dataset`` alias if/when added.
    try:
        cytome_Dataset = getattr(_cytome(), "Dataset", None)
    except ImportError:
        cytome_Dataset = None
    if cytome_Dataset is not None and cytome_Dataset not in candidates:
        candidates.append(cytome_Dataset)
    result = tuple(candidates)
    _cytome_dataset_classes._cache = result
    return result


def _is_cytome_dataset(obj):
    """True if ``obj`` is an open cytome Dataset (any known subclass)."""
    classes = _cytome_dataset_classes()
    return bool(classes) and isinstance(obj, classes)


def _open_or_use(cytome_path, mode="r"):
    """Return ``(ds, opened_here)`` for a cytome path-or-Dataset arg.

    - ``cytome_path`` is an open cytome Dataset → return as-is with
      ``opened_here=False`` (caller owns lifecycle; do NOT close).
    - ``cytome_path`` is ``str`` / ``Path`` → ``cytome.open(path, mode)``
      and return ``opened_here=True`` so the caller closes at exit.

    A path is opened read-only (``mode='r'``) unless the caller writes:
    scoring writes nothing into the store, so a default call cannot leave a
    cache behind in someone else's file (it left ``RNA_cell_depth`` and
    ``<modality>_tfidf_params`` there, keyed only by their length).
    ``write_to_cytome`` opens the store for writing when it persists. A
    cytome before 0.3.5 has no read-only mode and is opened as before.
    """
    if _is_cytome_dataset(cytome_path):
        return cytome_path, False
    cy = _cytome()
    try:
        return cy.open(cytome_path, mode=mode), True
    except TypeError:                     # cytome < 0.3.5: open(path) only
        return cy.open(cytome_path), True


def _emit_filter_log(tag, *, remove_lowly_expressed, expressed_pct,
                      expressed_min_num_cells_in_target_group, mu):
    """Print a one-line summary of which filters are active for the run.

    Surfaces silently-applied thresholds at the top of the run so users
    notice them BEFORE compute is spent. Mirrors what shows up in the
    function-body filter logic.
    """
    parts = [f"remove_lowly_expressed={remove_lowly_expressed}"]
    if remove_lowly_expressed:
        parts.append(f"expressed_pct={expressed_pct}")
        parts.append(
            f"expressed_min_num_cells_in_target_group="
            f"{expressed_min_num_cells_in_target_group}"
        )
    else:
        parts.append("(expressed_pct ignored)")
    parts.append(f"mu={mu}")
    print(
        f"[{tag}] Filters: " + ", ".join(parts) +
        " (entries with raw COSG score ≤ 0 dropped from output)"
    )


# =====================================================================
# API v2 (Round 7, 2026-05-22): K-HARD migration
#
# The old API exposed two parameters for the number of top features per
# group (``top_k`` and ``n_top_genes``) plus implicit mode switching
# (``top_k=None`` → standard ndarray output; ``top_k=N`` → sparse
# dict/long/dense). The new API uses a single parameter ``n_genes_user``
# (matches ``cosg.cosg`` on AnnData side) and an explicit
# ``output_format`` switch:
#
#   - output_format="ndarray" (default) → (names, scores) ndarrays
#   - output_format="dict"|"long"|"dense" → sparse outputs
#
# Old callers passing ``top_k`` or ``n_top_genes`` get an actionable
# TypeError with a migration hint. No silent aliasing — the meanings
# diverged too much to keep both names.
# =====================================================================

_VALID_OUTPUT_FORMATS = ("ndarray", "dict", "long", "dense")
_SPARSE_OUTPUT_FORMATS = {"dict", "long", "dense"}


def _validate_cosg_api_v2(deprecated_kwargs, output_format, iqr_normalize, fn_name):
    """Validate the public COSG cytome API.

    Returns the resolved ``(use_sparse_output, iqr_normalize)`` pair so
    callers don't repeat the logic.

    Raises TypeError for:
      - Old kwargs ``top_k`` / ``n_top_genes`` with a migration hint.
      - Any other unexpected kwargs (mimics Python's default kwarg error).

    Raises ValueError for invalid ``output_format``.
    """
    # K-HARD: explicit migration messages for the renamed kwargs.
    if "top_k" in deprecated_kwargs:
        old_val = deprecated_kwargs["top_k"]
        raise TypeError(
            f"{fn_name}() no longer accepts 'top_k'. "
            f"Use output_format='dict' (or 'long' / 'dense') and "
            f"n_genes_user={old_val!r} instead. "
            f"Migration: replace `top_k={old_val!r}` with "
            f"`n_genes_user={old_val!r}, output_format='dict'`."
        )
    if "n_top_genes" in deprecated_kwargs:
        old_val = deprecated_kwargs["n_top_genes"]
        raise TypeError(
            f"{fn_name}() no longer accepts 'n_top_genes' — renamed to "
            f"'n_genes_user' to match cosg.cosg on the AnnData side. "
            f"Migration: replace `n_top_genes={old_val!r}` with "
            f"`n_genes_user={old_val!r}`."
        )
    # Round 9 (2026-05-23): cytome_layer renamed to layer for symmetry
    # with the AnnData side's `layer=` kwarg.
    if "cytome_layer" in deprecated_kwargs:
        old_val = deprecated_kwargs["cytome_layer"]
        raise TypeError(
            f"{fn_name}() no longer accepts 'cytome_layer' — renamed to "
            f"'layer' to match cosg.cosg on the AnnData side. "
            f"Migration: replace `cytome_layer={old_val!r}` with "
            f"`layer={old_val!r}`."
        )
    # Forward any remaining stray kwargs as a regular TypeError for
    # Python-default-like behaviour.
    if deprecated_kwargs:
        bad = ", ".join(repr(k) for k in deprecated_kwargs)
        raise TypeError(
            f"{fn_name}() got unexpected keyword argument(s): {bad}"
        )

    if output_format not in _VALID_OUTPUT_FORMATS:
        raise ValueError(
            f"Invalid output_format={output_format!r}. Must be one of "
            f"{list(_VALID_OUTPUT_FORMATS)}."
        )

    use_sparse_output = output_format in _SPARSE_OUTPUT_FORMATS

    # SS-B: iqr_normalize default depends on output_format.
    # - Sparse outputs (dict/long/dense) always normalize (matches the
    #   AnnData reference chain indexByGene → iqrLogNormalize).
    # - ndarray output preserves raw COSG λ-values by default for
    #   backward compatibility with existing downstream plotting code.
    if iqr_normalize is None:
        iqr_normalize = use_sparse_output

    # When iqr_normalize is explicitly False on a sparse output, that's
    # legitimate (returns raw λ in the dict/long/dense). When True on
    # ndarray, that's also legitimate (applies IQR-log1p to the top-N).
    return use_sparse_output, bool(iqr_normalize)


def _format_top_k_output(scores_dict, groups_order_str, group_size_dict,
                          output_format, ds, modality, verbose, feature_name_col=None,
                          dense_entries=None):
    """Convert the top_k scores_dict into the requested return shape.

    Three output formats share the same underlying dict; only the
    post-processing shape differs:

    - ``'dict'`` (default): returns the dict as-is.
    - ``'long'``: pandas DataFrame, columns ``[group, feature, score, rank]``.
      Rank 1 = highest score within each group, descending.
    - ``'dense'``: pandas DataFrame indexed by ALL features in the cytome's
      feature table (rows), columns = groups, missing entries filled with 0.
    """
    import pandas as pd

    if output_format == "dict":
        return {
            "scores_dict": scores_dict,
            "groups_order": groups_order_str,
            "group_sizes": group_size_dict,
        }

    if output_format == "long":
        # Group entries by group name first so per-group sort is O(n log n_per_group)
        from collections import defaultdict
        per_group = defaultdict(list)
        for (group, feature), score in scores_dict.items():
            per_group[group].append((feature, score))

        rows = []
        for group in groups_order_str:
            entries = per_group.get(group, [])
            entries.sort(key=lambda x: -x[1])  # descending by score
            for rank, (feature, score) in enumerate(entries, start=1):
                rows.append({
                    "group": group,
                    "feature": feature,
                    "score": score,
                    "rank": rank,
                })
        scores_df = pd.DataFrame(rows, columns=["group", "feature", "score", "rank"])
        return {
            "scores_df": scores_df,
            "groups_order": groups_order_str,
            "group_sizes": group_size_dict,
        }

    if output_format == "dense":
        # Row index = ALL features in the cytome's feature table for this
        # modality. Columns = groups_order_str. Fill missing with 0.
        all_features = _get_all_feature_names(ds, modality, feature_name_col)
        n_features = len(all_features)
        n_groups = len(groups_order_str)

        if verbose and n_features * n_groups > 50_000_000:
            print(
                f"[cytome_cpu] WARNING: dense format will allocate "
                f"{n_features:,} × {n_groups} × 4 bytes "
                f"≈ {(n_features * n_groups * 4) / 1e9:.1f} GB. "
                f"Consider 'long' or 'dict' if RAM is tight."
            )

        # Allocate float32 to keep memory in check. The pivot fills a numpy
        # array, which then becomes the frame (loc-based assignment is slow on
        # large indices, and pandas 3's copy-on-write makes a frame's .values
        # read-only, so the frame cannot be filled in place).
        values = np.zeros((n_features, n_groups), dtype=np.float32)
        if dense_entries is not None:
            # by feature index: a name the table holds twice (10 gene symbols
            # of a 10x GRCh38 reference) gets each feature's own score; by
            # name, one of its rows stayed 0 and the other took either score
            g_idx, f_idx, vals = (np.asarray(a) for a in dense_entries)
            if len(f_idx):
                values[f_idx.astype(np.int64), g_idx.astype(np.int64)] = vals
        else:
            feature_to_row = {f: i for i, f in enumerate(all_features)}
            group_to_col = {g: i for i, g in enumerate(groups_order_str)}
            for (group, feature), score in scores_dict.items():
                row = feature_to_row.get(feature)
                if row is None:
                    continue  # feature not in current feature table — defensive
                col = group_to_col[group]
                values[row, col] = score
        scores_df = pd.DataFrame(
            values,
            index=pd.Index(all_features, name="feature"),
            columns=groups_order_str,
        )
        return {
            "scores_df": scores_df,
            "groups_order": groups_order_str,
            "group_sizes": group_size_dict,
        }

    # Should never reach here — output_format is validated upstream.
    raise ValueError(f"Invalid output_format='{output_format}'.")


def _warn_repeated_names(collided, output_format):
    """Say when the dict or long output, keyed by feature name, could not hold
    two features that share a name (a 10x GRCh38 reference holds 10 gene
    symbols twice): each such name keeps its better score. The dense and
    ndarray outputs keep every feature, by position."""
    if not collided or output_format == "dense":
        return
    import warnings
    shown = sorted(collided)[:5]
    warnings.warn(
        f"{len(collided)} feature name(s) belong to more than one scored feature, e.g. "
        f"{shown}; output_format={output_format!r} is keyed by name and keeps each "
        f"name's better score. output_format='dense' or 'ndarray' keeps every feature, "
        f"and feature_name_col= can name features by a unique column (e.g. 'gene_id').",
        UserWarning, stacklevel=3)


def _get_all_feature_names(ds, modality, feature_name_col=None):
    """Return all feature names from the cytome's modality-specific feature
    table — used as the row index for the ``output_format='dense'`` DataFrame.

    Resolves the table + name column via the cytome modality registry
    (same single source of truth used by the rest of the streaming code),
    so gene-name resolution is consistent with how features are named in
    the dict / long outputs.
    """
    from cytome import modality_feature_table_info
    feature_table, idx_col, name_col = modality_feature_table_info(ds, modality, feature_name_col)
    rows = ds._conn.execute(
        f'SELECT "{name_col}" FROM {feature_table} ORDER BY {idx_col}'
    ).fetchall()
    return [str(r[0]) for r in rows]


def _get_labels_and_gene_names(ds, groupby, cell_mask=None, modality="RNA", load_gene_names=True, feature_name_col=None):
    """Extract group labels, feature names, and categorical ordering from a cytome dataset.

    When load_gene_names=False, gene_names is returned as None (lazy loading).

    Modality dispatch (RNA/GA/ATAC/tiles) is delegated to
    ``cytome.modality_feature_table_info`` — single source
    of truth shared with cytome.io and piaso plotting.
    """
    from cytome import modality_feature_table_info

    # Check the column exists BEFORE interpolating it into the query.
    #
    # SQLite has a compatibility misfeature: a double-quoted token that does
    # not resolve to a column is treated as a *string literal* rather than
    # raising. So `SELECT cell_idx, "Leiden" FROM cells` on a cytome with no
    # `Leiden` column returns the string 'Leiden' for every row -- COSG then
    # sees exactly one group, named after the typo, and every downstream number
    # is internally consistent. A run reported `1 groups` and a plausible
    # marker table for a column that was never there.
    #
    # A confident wrong answer is worse than a crash, so: look first.
    _cols = [r[1] for r in ds._conn.execute("PRAGMA table_info(cells)")]
    if groupby not in _cols:
        # SQLite column names are case-insensitive, so a case variant resolves
        # to the real column and must not be rejected here -- the guard exists
        # for a name that resolves to nothing at all.
        _close = [c for c in _cols if c.lower() == str(groupby).lower()]
        if not _close:
            raise KeyError(
                f"groupby={groupby!r} is not a column of the cytome cells "
                f"table. Available: {sorted(_cols)}"
            )

    # Read group labels from SQLite cells table
    rows = ds._conn.execute(
        f'SELECT cell_idx, "{groupby}" FROM cells ORDER BY cell_idx'
    ).fetchall()
    all_labels = np.array([r[1] for r in rows], dtype=object)

    if cell_mask is not None:
        cell_mask = np.asarray(cell_mask)
        if cell_mask.dtype == bool:
            keep_idx = np.where(cell_mask)[0]
        else:
            keep_idx = np.sort(cell_mask)
        labels = all_labels[keep_idx]
    else:
        labels = all_labels

    # Every cell scored needs a label: a missing one (NULL, or NaN) crashed in
    # the sort of the groups with a TypeError.
    _missing = np.fromiter((l is None or (isinstance(l, float) and l != l) for l in labels),
                           dtype=bool, count=len(labels))
    if _missing.any():
        raise ValueError(
            f"{int(_missing.sum())} of the {len(labels)} cells scored have no label in "
            f"cells.{groupby} (NULL or NaN). Leave them out with cell_mask (a "
            f"boolean mask that is True for the labelled cells).")

    # Check _column_meta for categorical ordering (matches AnnData's .cat.categories)
    categories = None
    try:
        cat_row = ds._conn.execute(
            "SELECT categories FROM _column_meta "
            "WHERE table_name = 'cells' AND column_name = ? "
            "AND dtype IN ('categorical', 'ordered_categorical')",
            (groupby,),
        ).fetchone()
        if cat_row and cat_row[0]:
            categories = json.loads(cat_row[0])
    except Exception:
        pass

    # Modality → entity-table routing via the cytome registry. Closes the
    # silent-fall-through bug for modality='GA' (was reading 'genes').
    feature_table, idx_col, name_col = modality_feature_table_info(ds, modality, feature_name_col)

    if load_gene_names:
        feat_rows = ds._conn.execute(
            f'SELECT {idx_col}, "{name_col}" FROM {feature_table} ORDER BY {idx_col}'
        ).fetchall()
        gene_names = np.array([r[1] for r in feat_rows], dtype=object)
    else:
        gene_names = None

    return labels, gene_names, all_labels, categories


def _build_group_mapping(labels, categories=None):
    """Build group ordering and cell-to-group index array.

    Parameters
    ----------
    labels : np.ndarray
        Group label for each cell.
    categories : list of str, optional
        Categorical ordering from _column_meta.  When provided, this order
        is used instead of alphabetical sort so that the column layout
        matches standard COSG (which reads ``group_info.cat.categories``).

    Returns
    -------
    groups_order : list of str
        Group names in canonical order.
    cell_group_indices : np.ndarray of int32
        Group index for each cell in labels.
    group_sizes : np.ndarray of float64
        Number of cells per group.
    """
    unique_vals = set(labels)

    if categories is not None:
        # Use the stored categorical order, filtering to groups that
        # actually appear in the data (same logic as COSG's cosg.py).
        groups_order = [c for c in categories if c in unique_vals]
    else:
        unique_sorted = np.unique(labels)

        def _is_all_numeric(groups):
            try:
                [float(x) for x in groups]
                return True
            except (ValueError, TypeError):
                return False

        if _is_all_numeric(unique_sorted):
            groups_order = sorted(unique_sorted, key=lambda x: float(x))
        else:
            groups_order = sorted(unique_sorted)

    group_to_idx = {g: i for i, g in enumerate(groups_order)}
    cell_group_indices = np.array(
        [group_to_idx[l] for l in labels], dtype=np.int32
    )
    group_sizes = np.zeros(len(groups_order), dtype=np.float64)
    for g_idx in cell_group_indices:
        group_sizes[g_idx] += 1

    return groups_order, cell_group_indices, group_sizes


def _get_feature_table_info(ds, modality, feature_name_col=None):
    """Return (feature_table, idx_col, name_col) for a given modality.

    Thin wrapper over ``cytome.modality_feature_table_info``
    so this and ``_get_labels_and_gene_names`` use the single canonical
    registry. Closes the GA-modality silent-fall-through bug.
    """
    from cytome import modality_feature_table_info
    return modality_feature_table_info(ds, modality, feature_name_col)


# Modality-appropriate normalization when layer='auto' (the default). COSG scores
# a cosine specificity, so the layer must match the modality's proper transform:
# RNA/GA are (log-)count data → log1p; ATAC/tiles are binary-ish accessibility →
# TF-IDF. `infog` (enhanced RNA) and `counts` (piaso-free raw) stay explicit opt-ins.
_AUTO_LAYER_BY_MODALITY = {
    "RNA": "log1p",
    "GA": "log1p",
    "ATAC": "tfidf",
    "tiles": "tfidf",
}

#: Layers COSG can compute on the fly. Each needs a piaso normalizer and has a
#: branch in _resolve_chunk_normalizer; 'counts' is handled before this, since
#: it needs no transform at all.
_ON_THE_FLY_LAYERS = ("log1p", "infog", "tfidf")


def _unknown_layer_message(layer):
    return (
        f"Unknown layer='{layer}'. Known: "
        f"'counts' | 'log1p' | 'infog' | 'tfidf' (or any layer name "
        f"already materialised in the cytome — pass compute_on_fly=False "
        f"and matching layer to read it directly)."
    )


def _stored_matrix_is_integer(ds, modality, n_probe=20000):
    """Whether the stored base matrix looks like raw integer counts.

    Reads at most ``n_probe`` stored values -- cheap, and enough: a
    log-normalised matrix has non-integers within its first chunk.
    Returns ``None`` when it cannot tell, so callers can leave behaviour
    unchanged rather than guess.
    """
    try:
        import numpy as _np
        # cytome >= 0.3.0 records this at write time, when the data was in
        # hand and the answer was certain. Prefer it: re-deriving a fact the
        # writer already knew is how the probe below shipped as dead code.
        try:
            row = ds._conn.execute(
                "SELECT is_integer FROM matrix_meta WHERE matrix_name = ?",
                (f"{modality}_counts",),
            ).fetchone()
            if row is not None and row[0] is not None:
                return bool(row[0])
        except Exception:
            pass                    # older schema: fall through to the probe
        name = f"{modality}_counts"
        meta = ds.matrix_meta(name)
        if meta is None:
            return None
        # cytome >= 0.3.0 records this at write time, which is the only moment
        # it is both cheap and certain. Prefer it: re-deriving a fact the
        # writer already knew is how this probe came to be dead code once.
        recorded = meta.get("is_integer") if hasattr(meta, "get") else None
        if recorded is not None:
            return bool(recorded)
        for chunk in ds.iter_chunks(modality=modality, layer="counts",
                                    batch_size=512):
            # iter_chunks yields (matrix, row_indices), not a bare matrix.
            mat = chunk[0] if isinstance(chunk, tuple) else chunk
            if _np.issubdtype(getattr(mat, "dtype", _np.float64), _np.integer):
                return True                      # stored as an integer dtype
            data = getattr(mat, "data", None)
            if data is None:
                data = _np.asarray(mat).ravel()
            if len(data) == 0:
                continue
            probe = _np.asarray(data[:n_probe], dtype="float64")
            return bool(_np.allclose(probe, _np.round(probe)))
    except Exception:
        return None
    return None


def _resolve_auto_layer(layer, modality, ds=None):
    """Resolve ``layer='auto'`` to a modality-appropriate normalization.

    Idempotent: any explicit layer (``'counts'|'log1p'|'infog'|'tfidf'`` or a
    materialised layer name) is returned unchanged. Called at the top of both
    low-level interpreters so every read path (cpu/gpu, on-the-fly/materialised)
    resolves ``auto`` identically regardless of entry point.

    **``auto`` looks at the data, not only at the modality.** ``from_anndata``
    stores whatever ``adata.X`` held in ``{modality}_counts`` -- including an
    already log-normalised matrix -- so resolving RNA to ``log1p`` by table
    lookup alone would log1p it a second time. The result stays plausible
    (correlation ~0.99 against the in-memory reference) while being wrong by
    enough to read as biology, which is the dangerous kind of wrong. When the
    stored matrix is not integer, ``auto`` uses it as given and says so.
    """
    if layer != "auto":
        return layer
    resolved = _AUTO_LAYER_BY_MODALITY.get(modality)
    if resolved is None:
        raise ValueError(
            f"run_cosg_cytome(layer='auto', modality={modality!r}): no default "
            f"normalization is defined for this modality. Pass an explicit "
            f"layer ('log1p' | 'infog' | 'tfidf' | 'counts')."
        )
    if resolved == "log1p" and ds is not None:
        # cytome >= 0.3.0 does not write `{modality}_counts` unless the values
        # really are counts, so on those files the matrix simply is not there.
        # Fall back to whatever `adata.X` was, which is what makes `auto`
        # equivalent to `cosg.cosg(adata)` -- the documented contract.
        if ds.matrix_meta(f"{modality}_counts") is None:
            recorded = None
            try:
                recorded = ds.metadata.get("_anndata_X_layer")
            except Exception:
                recorded = None
            if recorded and recorded.startswith(f"{modality}_"):
                main = recorded[len(modality) + 1:]
                import warnings
                warnings.warn(
                    f"run_cosg_cytome(layer='auto'): this file has no "
                    f"'{modality}_counts' -- cytome >= 0.3.0 only writes that "
                    f"name for raw counts. Using '{main}', the matrix that "
                    f"was adata.X. Pass an explicit layer to silence this.",
                    UserWarning, stacklevel=3,
                )
                return main
        is_int = _stored_matrix_is_integer(ds, modality)
        if is_int is False:
            import warnings
            warnings.warn(
                f"run_cosg_cytome(layer='auto'): {modality}_counts does not "
                f"hold integer counts, so it is already normalised -- applying "
                f"log1p would normalise it twice. Using the stored values as "
                f"given (equivalent to layer='counts'). Pass an explicit layer "
                f"to silence this.",
                UserWarning, stacklevel=3,
            )
            return "counts"
    return resolved


def _resolve_chunk_normalizer(layer, ds, modality, use_cached_stats):
    """Return a per-chunk normalizer callable, or ``None`` for raw counts.

    ``counts`` needs no transform and ``log1p`` is implemented natively in
    :mod:`cosg._normalize`, so **neither touches PIASO** -- including the
    public default, which is ``log1p`` for RNA. That matters because PIASO
    requires ``cosg>=1.0.3``: when COSG's default path imported PIASO back,
    the two were mutually dependent in practice, and COSG's headline streaming
    feature could not be exercised without installing the package that depends
    on COSG.

    ``infog`` and ``tfidf`` do lazy-import PIASO, which is correct -- they are
    PIASO normalizations, and asking for one means you want PIASO. Deferring
    to call time also sidesteps the import-order risk in ``piaso/__init__.py``
    (which eagerly imports cosg in ``_runGDR.py:13``): by the time this
    function fires, PIASO is fully loaded.

    Supported on-the-fly layers: ``log1p``, ``infog``, ``tfidf``.
    Materialised non-counts layers (e.g. ``RNA_infog`` already on disk)
    are read directly via ``iter_chunks(layer=...)`` by the caller —
    this function only fires when the matrix is missing AND
    ``compute_on_fly=True``.

    Parameters
    ----------
    layer : str
        ``'counts'``, ``'log1p'``, ``'infog'``, or ``'tfidf'``.
    ds : cytome.Dataset
    modality : str
    use_cached_stats : bool

    Returns
    -------
    Callable[[csr_matrix, np.ndarray], csr_matrix or ndarray] or None
        ``f(chunk, row_indices) -> normalized_chunk``. None means
        no transform — caller uses the raw counts chunk as-is.
    """
    layer = _resolve_auto_layer(layer, modality, ds=ds)
    if layer == "counts":
        return None

    # Reject an unrecognised layer BEFORE trying to import piaso. Otherwise a
    # typo'd layer name on a machine without piaso reports "needs piaso for
    # on-the-fly 'totaly_made_up' normalization" — telling the user to install
    # a package that could never have provided that layer.
    if layer not in _ON_THE_FLY_LAYERS:
        raise ValueError(_unknown_layer_message(layer))

    # Every non-'counts' layer needs a piaso normalizer. Import piaso lazily
    # (only here, at call time) so `import cosg` stays piaso-free; if piaso is
    # absent, raise a clear, actionable error instead of a bare ImportError on
    # a deep submodule path.
    # Only infog and tfidf need PIASO -- they are PIASO methods. log1p is
    # handled natively (cosg/_normalize.py) so the *default* call has no PIASO
    # dependency; PIASO requires cosg, and having COSG's default path require
    # PIASO back made the two mutually dependent in practice.
    if layer in ("infog", "tfidf"):
        try:
            import piaso  # noqa: F401
        except ImportError as exc:  # pragma: no cover - exercised only w/o piaso
            raise ImportError(
                f"cosg.run_cosg_cytome(layer='{layer}') needs PIASO: "
                f"'{layer}' is a PIASO normalization.\n"
                f"    pip install piaso-tools\n"
                f"(the distribution is `piaso-tools`; the import is `piaso`.) "
                f"layer='log1p' and layer='counts' need no PIASO."
            ) from exc

    if layer == "infog":
        # piaso.tools._normalization is loaded eagerly by piaso/_runGDR.py:13
        # so this top-level-equivalent lazy import is safe regardless of
        # caller import order.
        from piaso.plotting._plotEmbedding import _ensure_infog_params
        from piaso.tools._normalization import _normalize_chunk_infog
        params = _ensure_infog_params(
            ds, modality, use_cached_stats=use_cached_stats,
        )
        if params is None:
            raise RuntimeError(
                f"infog params not available for modality '{modality}'. "
                f"Run piaso.tl.infog(ds, modality='{modality}', "
                f"save_layer=False) to populate the cache, then retry."
            )
        cd = np.asarray(params["cell_depth"], dtype=np.float64)
        ig = np.asarray(params["inv_gene_depth"], dtype=np.float64)
        scale = float(params["scale"])
        counts_sum = float(params["counts_sum"])
        thr = params.get("threshold")
        # a column chunk (col_start given) gets its own features' factors
        return lambda chunk, indices, col_start=None: _normalize_chunk_infog(
            chunk, cd[indices], _feature_slice(ig, col_start, chunk), scale, counts_sum, thr,
        )

    if layer == "tfidf":
        from piaso.plotting._plotEmbedding import _ensure_tfidf_params
        from piaso.tools._runTFIDF import _normalize_chunk_tfidf
        params = _ensure_tfidf_params(
            ds, modality, use_cached_stats=use_cached_stats,
        )
        if params is None:
            raise RuntimeError(
                f"tfidf params not available for modality '{modality}'. "
                f"Run piaso.tl.compute_tfidf_stats(ds, modality='{modality}') "
                f"to populate the cache, then retry."
            )
        cd = np.asarray(params["cell_depth"], dtype=np.float64)
        idf = np.asarray(params["idf"], dtype=np.float64)
        scale_factor = float(params.get("scale_factor", 1e4))
        # A column chunk (col_start given) gets its own features' IDF: the
        # function indexes idf by the chunk's column indices, which are local
        # to the chunk, so every chunk took the first chunk's IDF. COSG's
        # cosine ignores a positive factor per feature, so its scores were
        # right; the normalized values were not.
        return lambda chunk, indices, col_start=None: _normalize_chunk_tfidf(
            chunk, cd[indices], _feature_slice(idf, col_start, chunk), scale_factor,
        )

    if layer == "log1p":
        from ._normalize import normalize_chunk_log1p as _normalize_chunk_log1p
        # Scale: prefer one recorded on the file (piaso.pp.normalize_log1p
        # records it at write-back time, so a COSG run matches whatever PIASO
        # already wrote); else the conventional 1e4.
        scale = 1e4
        params = ds.metadata.get(f"{modality}_log1p_params")
        if params is not None:
            scale = float(params.get("scale_factor", 1e4))
        # Each cell's depth from the data, never from a cache: on a row-major
        # matrix it is the row sum of the counts chunk being normalized, so
        # nothing is read twice, nothing is written into the store, and
        # nothing can go stale (the cache, RNA_cell_depth, was trusted
        # whenever its length matched, so a matrix replaced by one of the same
        # shape kept the old depths). The same numbers as cytome's
        # matrix_sums, which computed that cache.
        cell_depth = _log1p_depth(ds, modality)
        if cell_depth is None:
            def _from_rows(chunk, indices, col_start=None):
                if col_start is not None:     # part of each row: its sums are not depths
                    raise ValueError("log1p on a column chunk needs every cell's depth, "
                                     "which a row-major counts matrix gives per row chunk")
                return _normalize_chunk_log1p(chunk, _row_sums(chunk), scale)
            return _from_rows
        return lambda chunk, indices, col_start=None: _normalize_chunk_log1p(
            chunk, cell_depth[indices], scale,
        )

    # Unreachable: the layer was validated against _ON_THE_FLY_LAYERS above,
    # and every member has a branch. Kept as a guard so adding a name to that
    # tuple without a branch fails loudly rather than returning None (which
    # the caller would read as "raw counts, no transform").
    raise AssertionError(
        f"layer='{layer}' is in _ON_THE_FLY_LAYERS but has no normalizer branch"
    )


def _feature_slice(values, col_start, chunk):
    """A per-feature vector cut to a column chunk's features (all of it for a
    row chunk, ``col_start`` None)."""
    if col_start is None:
        return values
    return values[int(col_start):int(col_start) + chunk.shape[1]]


def _row_sums(chunk):
    """Each row's sum, as cytome's matrix_sums takes it (in the chunk's dtype,
    then float64), so that the depth of log1p is the one it computed."""
    return np.asarray(chunk.sum(axis=1), dtype=np.float64).ravel()


def _log1p_depth(ds, modality):
    """Every cell's depth for log1p when it cannot come from the chunk being
    normalized, else None (a row-major counts matrix: its row sums).

    - ATAC: ``cells.n_fragments`` when it is filled, the depth cytome's
      ``modality_cell_depth`` uses for ATAC (the importer stores it).
    - A column-major counts matrix: its column chunks hold part of each row,
      so the row sums are added up over one pass of the columns first, as
      ``matrix_sums`` does.
    Read only: no cache is consulted or written.
    """
    n_cells = int(ds.n_cells)
    if modality == "ATAC":
        try:
            n_frags = np.asarray(ds.cells["n_fragments"], dtype=np.float64)
            if n_frags.size == n_cells and float(n_frags.sum()) > 0:
                return n_frags
        except Exception:
            pass
    layer = _column_major_layer(ds, modality, "counts")
    if layer is None:
        return None
    depth = np.zeros(n_cells, dtype=np.float64)
    for _cs, _ce, chunk in layer.iter_columns():
        chunk = chunk.tocsc()
        if chunk.nnz:
            depth += np.asarray(chunk.sum(axis=1), dtype=np.float64).ravel()
    return depth


def _resolve_layer_to_read(ds, modality, layer, compute_on_fly):
    """Decide which on-disk layer iter_chunks should request.

    - ``layer == 'counts'``: read ``counts`` (no transform).
    - Materialised ``{modality}_{layer}`` exists: read it
      directly (storage > recomputation), no normalizer.
    - Layer missing AND ``compute_on_fly`` True: read ``counts``,
      let the normalizer transform on the fly.
    - Layer missing AND ``compute_on_fly`` False: raise actionable.

    Returns ``(layer_to_read, use_normalizer)`` where ``use_normalizer``
    indicates whether ``_resolve_chunk_normalizer`` should produce a
    transform.
    """
    layer = _resolve_auto_layer(layer, modality, ds=ds)
    if layer == "counts":
        return "counts", False
    matrix_name = f"{modality}_{layer}"
    if ds.matrix_meta(matrix_name) is not None:
        return layer, False
    if not compute_on_fly:
        raise ValueError(
            f"Matrix '{matrix_name}' not found in cytome and "
            f"compute_on_fly=False. Either materialise it (e.g. "
            f"piaso.pp.normalize_log1p(ds, modality='{modality}', "
            f"save_layer=True)) or pass compute_on_fly=True to compute "
            f"per-chunk on the fly from {modality}_counts."
        )
    return "counts", True



def _matrix_label(modality, layer_read, layer_asked, normalizer):
    """How an error names the matrix scored."""
    name = f"the {modality}_{layer_read} matrix"
    if normalizer is not None:
        name += f" (with {layer_asked} applied on the fly)"
    return name


def _column_major_layer(ds, modality, layer):
    """The layer as a MeasurementLayer when it is stored column-major, else None.

    A matrix written by the importer or quantifier with ``--layout csc`` has
    CSC chunks and no row chunks. ``iter_chunks`` refuses it; the column walk
    below reads it as the layout it is.
    """
    row = ds._conn.execute(
        "SELECT COALESCE(n_chunks, 0), COALESCE(has_csc, 0) FROM matrix_meta WHERE matrix_name = ?",
        (f"{modality}_{layer}",),
    ).fetchone()
    if row is None or row[0] != 0 or not row[1]:
        return None
    from cytome.core.measurement import MeasurementLayer
    return MeasurementLayer(ds._conn, f"{modality}_{layer}")


# Working memory one block of features may take while it is scored.
_BLOCK_BYTES = 64 * 1024 * 1024
_BLOCK_MIN, _BLOCK_MAX = 2_000, 200_000


# The most the row-major path's per-feature sums may hold in one pass under
# feature_batching="chromosome" (``max_accumulator_bytes`` overrides it).
_ACCUM_BUDGET_BYTES = 2 * 1024 ** 3


def _chrom_passes(chrom_order, chrom_ranges, n_groups, n_batches, track_nnz, budget):
    """The chromosomes, in order, packed into passes whose sums fit ``budget``
    bytes: per chromosome ``dot`` (features x groups) and ``sq`` (features), each
    with the batch axis, and ``nnz`` (features x groups) when it is tracked, four
    bytes each. A chromosome larger than the budget is a pass of its own."""
    passes, cur, cur_bytes = [], [], 0
    for chrom in chrom_order:
        cs, ce = chrom_ranges[chrom]
        n_feat = ce - cs
        b = n_batches * (n_feat * n_groups * 4 + n_feat * 4) + (n_feat * n_groups * 4 if track_nnz else 0)
        if cur and cur_bytes + b > budget:
            passes.append(cur)
            cur, cur_bytes = [chrom], b
        else:
            cur.append(chrom)
            cur_bytes += b
    if cur:
        passes.append(cur)
    return passes


def _block_features(n_groups, n_batches, block_features=None):
    """Features scored together. Per feature the scan holds the sums, the
    counts and the scoring temporaries: about sixty bytes per group, more
    with a batch axis."""
    if block_features:
        return max(1, int(block_features))
    per_feature = max(n_groups, 1) * (60 + 16 * max(n_batches - 1, 0))
    return int(min(max(_BLOCK_BYTES // per_feature, _BLOCK_MIN), _BLOCK_MAX))


def _scan_column_major(layer, keep_rows, cell_group_indices, n_groups, group_sizes,
                       batch_index, n_batches, grp_count_accum, batch_cell_number_threshold,
                       chrom_ranges, chunk_normalizer, track_nnz, mu,
                       min_total_count, thresholds,
                       n_keep, collect_positive, verbose, block_features=None, check_sums=None,
                       candidates_all=False):
    """COSG over a column-major matrix, scoring each feature as it is read.

    ``thresholds`` holds, per group, the fewest cells a feature must be
    nonzero in (``cosg._floor.floor_thresholds``: each group's own floor;
    0 for a group the filter leaves alone), or ``None`` when the expression
    filter is off; the nonzero counts are kept (``track_nnz``) when the
    filter or ``min_total_count`` needs them.

    A column chunk holds every entry of its features, so their statistics
    are complete when the chunk has been read: nothing has to be kept for
    the whole feature axis. The scan holds one block of features, scores
    it, keeps each group's candidates and positive scores, and lets the
    block go. Memory follows the block, not the number of features.

    The arithmetic is that of the accumulating path, feature by feature:
    ``dot = X^T @ Lam``, ``sq = sum(X^2)``, ``nnz`` per feature and group,
    then the cosine, the penalty and the filters. The sums are taken in
    float64 and stored in float32, so on counts the scores are identical;
    on a normalised layer they can differ in the last digit.

    ``candidates_all``: every feature is a candidate (the ndarray output,
    whose rows past a group's positive scores hold its next features by
    score, as on the unbatched path); otherwise only positive scores are.

    Returns ``(group_candidates, positives, n_chunks, n_cells, seconds)``:
    per group, ``(feature_indices, scores)`` sorted by score then index
    (``cosg._select``'s rule), and a list of float32 arrays of its positive
    scores (or ``None``).
    """
    t0 = time.time()
    n_rows, n_cols = layer.shape
    G = int(n_groups)
    if keep_rows is None:
        keep_rows = np.arange(n_rows, dtype=np.int64)
    keep_rows = np.asarray(keep_rows, dtype=np.int64)
    # The group (and batch) of every stored row; -1 for a row not kept.
    group_of_row = np.full(n_rows, -1, dtype=np.int64)
    group_of_row[keep_rows] = np.asarray(cell_group_indices, dtype=np.int64)
    if batch_index is not None:
        batch_of_row = np.zeros(n_rows, dtype=np.int64)
        batch_of_row[keep_rows] = np.asarray(batch_index, dtype=np.int64)
    B = int(n_batches)
    # Columns that belong to a chromosome range are scored; others never are.
    covered = np.zeros(n_cols, dtype=bool)
    for cs, ce in chrom_ranges.values():
        covered[cs:ce] = True
    if thresholds is not None:
        thresholds = np.asarray(thresholds, dtype=np.float64)
    all_rows = np.arange(n_rows)

    cand_idx = [[] for _ in range(G)]
    cand_val = [[] for _ in range(G)]
    cand_n = [0] * G
    positives = [[] for _ in range(G)] if collect_positive else None

    def trim(k):
        idx = np.concatenate(cand_idx[k]); val = np.concatenate(cand_val[k])
        # by score, then by feature index: the same candidates whatever the
        # block size or the order they arrived in
        order = np.lexsort((idx, -val))[:n_keep]
        cand_idx[k], cand_val[k], cand_n[k] = [idx[order]], [val[order]], len(order)

    def finish(start, dots, sqs, nnzs):
        n = sum(d.shape[1] for d in dots)
        if check_sums is not None:
            check_sums(np.concatenate(sqs, axis=1), start)
        rows = None                       # the features scored; None for all
        nnz = None
        scores = None
        if nnzs is not None:
            # A feature that fails the filters in every group scores -1
            # everywhere, whatever its cosine. On peaks and tiles that is most
            # features, so they are dropped before the dense arithmetic. The
            # penalty is computed row by row, so the rows kept are unaffected.
            nnz = np.concatenate(nnzs, axis=0)
            alive = covered[start:start + n].copy()
            if min_total_count > 0:
                alive &= nnz.sum(axis=1) > min_total_count
            if thresholds is not None:
                alive &= (nnz >= thresholds[None, :])[:, np.asarray(group_sizes) > 0].any(axis=1)
            rows = np.flatnonzero(alive)
            if not len(rows):
                if not candidates_all:
                    return
                # every feature of the block failed the filters: -1 in every
                # group; the ndarray output lists them as the unbatched path
                # does (the block vanished, and its rows became "" and 0.0)
                scores, rows = np.full((n, G), -1.0), None
            else:
                nnz = nnz[rows]
        if scores is None:
            dot = np.concatenate(dots, axis=1)                          # (B, n, G)
            sq = np.concatenate(sqs, axis=1)                            # (B, n)
            if rows is not None:
                dot, sq = dot[:, rows], sq[:, rows]
            dot = dot.astype(np.float32)
            sq = sq.astype(np.float32)
            # min_total_count has already removed its rows (alive)
            scores = _score_from_sums(
                dot[0] if batch_index is None else dot,
                sq[0] if batch_index is None else sq,
                group_sizes, nnz, mu=mu, thresholds=thresholds,
                batch_group_sizes=None if batch_index is None else grp_count_accum,
                batch_cell_number_threshold=batch_cell_number_threshold)
            del dot, sq
            if nnz is None:
                scores[~covered[start:start + n], :] = -1
        if candidates_all and rows is not None:
            # the features the filters dropped before scoring score -1 in every
            # group, as on the unbatched path
            full = np.full((n, G), -1.0)
            full[rows] = scores
            scores, rows = full, None
        for k in range(G):
            col = scores[:, k]
            pos = np.flatnonzero(col > 0)
            if positives is not None and len(pos):
                positives[k].append(col[pos].astype(np.float32))
            if candidates_all:
                loc = top_n(col, min(n_keep, col.shape[0]))
                feat, val = loc, col[loc]
            else:
                if not len(pos):
                    continue
                feat = pos if rows is None else rows[pos]
                val = col[pos]
                if len(pos) > n_keep:
                    top = np.lexsort((feat, -val))[:n_keep]
                    feat, val = feat[top], val[top]
            if not len(feat):
                continue
            cand_idx[k].append(feat.astype(np.int64) + start)
            cand_val[k].append(val.astype(np.float64))
            cand_n[k] += len(feat)
            if cand_n[k] > 4 * n_keep:
                trim(k)

    block = _block_features(G, B, block_features)
    dots, sqs, nnzs = [], [], ([] if track_nnz else None)
    held, block_start, n_chunks, expected = 0, 0, 0, 0
    for col_start, col_end, chunk in layer.iter_columns():
        if col_start != expected:
            raise ValueError(
                f"column chunks are not contiguous: expected a chunk at column "
                f"{expected}, found {col_start}")
        expected = col_end
        c = col_end - col_start
        if chunk_normalizer is not None:
            chunk = chunk_normalizer(chunk.tocsr(), all_rows, col_start=col_start).tocsc()
            chunk.sort_indices()
        g = group_of_row[chunk.indices]
        ok = g >= 0
        g = g[ok]
        col = np.repeat(np.arange(c, dtype=np.int64), np.diff(chunk.indptr))[ok]
        v = chunk.data[ok].astype(np.float64)
        if batch_index is None:
            flat = col * G + g
            feat = col
        else:
            b = batch_of_row[chunk.indices[ok]]
            feat = b * c + col
            flat = feat * G + g
        dots.append(np.bincount(flat, weights=v, minlength=B * c * G).reshape(B, c, G))
        sqs.append(np.bincount(feat, weights=v * v, minlength=B * c).reshape(B, c))
        if nnzs is not None:
            # a stored zero is not an expressing cell
            nnzs.append(np.bincount(col * G + g, weights=(v != 0), minlength=c * G)
                        .reshape(c, G).astype(np.int32))
        held += c
        n_chunks += 1
        if held >= block:
            finish(block_start, dots, sqs, nnzs)
            block_start += held
            dots, sqs, nnzs, held = [], [], ([] if track_nnz else None), 0
        if verbose and n_chunks % 2000 == 0:
            print(f"    {n_chunks} column chunks")
    if held:
        finish(block_start, dots, sqs, nnzs)

    group_candidates = {}
    for k in range(G):
        if cand_n[k]:
            trim(k)
            group_candidates[k] = (cand_idx[k][0], cand_val[k][0])
        else:
            group_candidates[k] = (np.empty(0, dtype=np.int64), np.empty(0))
    iqr_collector = None if positives is None else {k: positives[k] for k in range(G)}
    return group_candidates, iqr_collector, n_chunks, int(len(keep_rows)), time.time() - t0


def _get_chrom_col_ranges(ds, feature_table, idx_col):
    """Get chromosome → (start_col, end_col) column ranges from entity table.

    Peaks and tiles are registered in chromosome order during quantification,
    so each chromosome's features normally form one contiguous column range,
    and one GROUP BY query finds them. That was assumed, never checked: on a
    store whose features are not grouped by chromosome (chr1 split around
    chr2), the ranges overlapped, features were scored twice (duplicate names
    in a group's markers), and the IQR scale was padded to more features than
    there are. Now the ranges are checked to be disjoint and to cover every
    feature once; when they are not, the features are cut into runs of one
    chromosome in index order, which are (``chr1``, ``chr2``, ``chr1~2``, ...).
    Either way every feature is in exactly one block.
    """
    rows = ds._conn.execute(
        f'SELECT chr, MIN({idx_col}) as cs, MAX({idx_col}) + 1 as ce, COUNT(*) '
        f'FROM {feature_table} GROUP BY chr ORDER BY MIN({idx_col})'
    ).fetchall()
    n_total, lo, hi = ds._conn.execute(
        f'SELECT COUNT(*), MIN({idx_col}), MAX({idx_col}) FROM {feature_table}').fetchone()
    contiguous = bool(rows) and int(lo) == 0 and int(hi) + 1 == int(n_total)
    expect = 0
    for chrom, cs, ce, cnt in rows:
        if not contiguous:
            break
        contiguous = int(cs) == expect and int(ce) - int(cs) == int(cnt)
        expect = int(ce)
    if contiguous and expect == int(n_total):
        chrom_ranges = {}
        chrom_order = []
        for chrom, cs, ce, _cnt in rows:
            chrom_ranges[str(chrom)] = (int(cs), int(ce))
            chrom_order.append(str(chrom))
        return chrom_ranges, chrom_order

    # Runs of one chromosome in index order
    idx_chr = ds._conn.execute(
        f'SELECT {idx_col}, chr FROM {feature_table} ORDER BY {idx_col}').fetchall()
    idx = np.fromiter((r[0] for r in idx_chr), dtype=np.int64, count=len(idx_chr))
    if len(idx) and not np.array_equal(idx, np.arange(len(idx))):
        raise ValueError(
            f"COSG: the {feature_table} table's {idx_col} values are not 0..{len(idx) - 1} "
            f"(found {len(idx):,} rows from {int(idx.min())} to {int(idx.max())}), so its rows "
            f"cannot be matched to the matrix columns. Rebuild the table.")
    chroms = [str(r[1]) for r in idx_chr]
    chrom_ranges, chrom_order, seen = {}, [], {}
    start = 0
    for i in range(1, len(chroms) + 1):
        if i == len(chroms) or chroms[i] != chroms[start]:
            c = chroms[start]
            seen[c] = seen.get(c, 0) + 1
            key = c if seen[c] == 1 else f"{c}~{seen[c]}"
            chrom_ranges[key] = (start, i)
            chrom_order.append(key)
            start = i
    return chrom_ranges, chrom_order


def _check_feature_count(ds, modality, n_features, feature_name_col=None):
    """The feature (var) table and the modality's matrix must agree on the
    feature count. A 0-row or out-of-sync table (a modality whose feature
    table an older cytome.merge dropped) against a populated matrix crashed
    deep in the accumulation with a cryptic broadcast error, or, on the
    feature-batched path, left columns unscored. Fail early with a fix hint."""
    _mat_name, _mat_cols = None, None
    for _mn, _nc in ds._conn.execute(
        "SELECT matrix_name, n_cols FROM matrix_meta"
    ).fetchall():
        if _mn.split("_", 1)[0] == modality:
            _mat_name, _mat_cols = _mn, int(_nc)
            if _mn == f"{modality}_counts":
                break  # prefer the counts matrix as the reference
    if _mat_cols is not None and _mat_cols != n_features:
        from cytome import modality_feature_table_info
        _ftab, _idxc, _ = modality_feature_table_info(ds, modality, feature_name_col)
        raise ValueError(
            f"COSG: modality '{modality}' feature table '{_ftab}' has "
            f"{n_features:,} rows, but matrix '{_mat_name}' has {_mat_cols:,} "
            f"columns — the feature/var table is missing or out of sync (e.g. an "
            f"older cytome.merge carried the matrix but not the table). Rebuild "
            f"it so the '{_idxc}' order matches the matrix columns, e.g. "
            f"ds.set_entity('{_ftab}', <source {_ftab} dataframe>); ds.flush(), "
            f"then re-run COSG."
        )


def _lookup_feature_names(ds, feature_table, idx_col, name_col, indices):
    """Lookup feature names for specific indices (lazy, on-demand).

    Used by chromosome-batched COSG to avoid loading all 5M+ feature names
    upfront.  Only the final top-K indices (~36K) are looked up.
    """
    if len(indices) == 0:
        return np.array([], dtype=object)
    placeholders = ",".join(str(int(i)) for i in indices)
    rows = ds._conn.execute(
        f'SELECT {idx_col}, "{name_col}" FROM {feature_table} '
        f'WHERE {idx_col} IN ({placeholders}) ORDER BY {idx_col}'
    ).fetchall()
    idx_to_name = {r[0]: r[1] for r in rows}
    return np.array(
        [idx_to_name.get(int(i), f"feature_{i}") for i in indices],
        dtype=object
    )


def _persist_cosg_markers(cytome_path, result, groupby, modality,
                          n_genes_user, key_added):
    """Persist a ``run_cosg_cytome`` result to ``ds.metadata[key_added]`` as a
    long-form record list, regardless of the result's ``output_format``.

    Stored payload::

        {'groupby': ..., 'modality': ..., 'n_genes_user': ...,
         'records': [{'group', 'feature', 'score', 'rank'?}, ...],
         'expressed_pct'?: {group: fraction}}   # when the floor was per group
    """
    ds, _opened = _open_or_use(cytome_path, mode="r+")
    try:
        records = []
        if isinstance(result, dict) and "scores_df" in result:
            df = result["scores_df"]
            lower = {str(c).lower(): c for c in df.columns}
            if "score" in lower and "feature" in lower:
                # 'long' output — already tidy
                gcol = lower.get("group")
                fcol = lower["feature"]
                scol = lower["score"]
                rcol = lower.get("rank")
                for row in df.itertuples(index=False):
                    rd = row._asdict()
                    rec = {
                        "group": str(rd[gcol]) if gcol else "",
                        "feature": str(rd[fcol]),
                        "score": float(rd[scol]),
                    }
                    if rcol:
                        rec["rank"] = int(rd[rcol])
                    records.append(rec)
            else:
                # 'dense' output — features as index, groups as columns;
                # keep only strictly positive entries to stay bounded.
                for group in df.columns:
                    col = df[group]
                    nz = col[col > 0]
                    for feature, score in nz.items():
                        records.append({
                            "group": str(group),
                            "feature": str(feature),
                            "score": float(score),
                        })
        elif isinstance(result, dict) and "scores_dict" in result:
            for (group, feature), score in result["scores_dict"].items():
                records.append({
                    "group": str(group),
                    "feature": str(feature),
                    "score": float(score),
                })
        elif isinstance(result, dict) and "names" in result:
            names = result["names"]
            scores = result["scores"]
            groups = result["groups_order"]
            n_rows, n_groups = names.shape
            for k in range(n_groups):
                group = str(groups[k])
                for r in range(n_rows):
                    nm = names[r, k]
                    if nm is None or str(nm) == "":
                        continue
                    records.append({
                        "group": group,
                        "feature": str(nm),
                        "score": float(scores[r, k]),
                        "rank": r + 1,
                    })

        payload = {
            "groupby": groupby,
            "modality": modality,
            "n_genes_user": int(n_genes_user),
            "records": records,
        }
        # the detection floor of each group, when it was given per group
        if isinstance(result, dict) and isinstance(result.get("expressed_pct"), dict):
            payload["expressed_pct"] = {str(g): float(f)
                                        for g, f in result["expressed_pct"].items()}
        ds.metadata[key_added] = payload
        try:
            ds.flush()
        except Exception:
            pass
    finally:
        if _opened:
            ds.close()


def _normalize_cell_mask(cell_mask, n_cells):
    """``cell_mask`` as one boolean per cell of the store, or None.

    A boolean mask must have one entry per cell: one 10 short was accepted and
    the last 10 cells were silently left out. Integer indices name cells, each
    once and in range; they were cast to bool on the batch_key path, which
    raised IndexError. Every path below takes the boolean form.
    """
    if cell_mask is None:
        return None
    m = np.asarray(cell_mask)
    if m.ndim != 1:
        raise ValueError(f"cell_mask must be one-dimensional, got shape {m.shape}.")
    if m.dtype == bool:
        if m.shape[0] != n_cells:
            raise ValueError(
                f"cell_mask has {m.shape[0]:,} entries; a boolean mask needs one per cell "
                f"of the store ({n_cells:,}). Pass integer cell indices to name a subset.")
        out = m.copy()
    elif m.dtype.kind in "iu" or m.size == 0:
        m = m.astype(np.int64)
        if m.size and (m.min() < 0 or m.max() >= n_cells):
            raise ValueError(
                f"cell_mask holds cell indices outside 0..{n_cells - 1:,} "
                f"(min {int(m.min())}, max {int(m.max())}).")
        u = np.unique(m)
        if u.size != m.size:
            raise ValueError(f"cell_mask names {m.size - u.size:,} cell(s) more than once.")
        out = np.zeros(n_cells, dtype=bool)
        out[u] = True
    else:
        raise TypeError(
            f"cell_mask must be a boolean mask over the cells or integer cell indices, "
            f"got dtype {m.dtype}.")
    if not out.any():
        raise ValueError("cell_mask selects no cell.")
    return out


def _require_writable(target):
    """Raise now if ``write_to_cytome`` could not write ``target`` (a path,
    or an open Dataset): it used to fail only after the whole computation."""
    what = "run_cosg_cytome(write_to_cytome=True)"
    if _is_cytome_dataset(target):
        check = getattr(target, "require_writable", None)
        if callable(check):
            check(what)
        elif getattr(target, "read_only", False):
            raise PermissionError(f"{what}: the Dataset was opened read-only.")
        return
    check = getattr(_cytome(), "require_writable", None)
    if callable(check):
        check(target, what)


def _check_scoring_args(mu, n_genes_user, expressed_min_num_cells_in_target_group):
    """The same checks as cosg.cosg: the cytome paths accepted a negative mu
    and n_genes_user < 1 (which failed later with a cryptic error). A minimum
    of 0 cells is valid on every path: with a floor of 0 it filters nothing."""
    try:
        ok = bool(np.isfinite(mu)) and mu >= 0
    except TypeError:
        ok = False
    if not ok:
        raise ValueError(f"Parameter mu must be a non-negative number, got {mu!r}.")
    try:
        ok = (not isinstance(n_genes_user, bool) and n_genes_user >= 1
              and float(n_genes_user).is_integer())
    except (TypeError, ValueError):
        ok = False
    if not ok:
        raise ValueError(f"n_genes_user must be a whole number of at least 1, got {n_genes_user!r}.")
    if not expressed_min_num_cells_in_target_group >= 0:
        raise ValueError("expressed_min_num_cells_in_target_group must be at least 0.")


def _read_batch_labels(ds, batch_key, cell_mask, n_expected):
    """Batch labels for `batch_key`, aligned to the cells COSG will score.

    One owner so the CPU and feature-batched paths cannot drift on how a
    missing key, a None value or a cell_mask is handled.
    """
    try:
        allb = np.asarray(ds.cells[batch_key])
    except Exception as exc:
        raise ValueError(
            f"batch_key '{batch_key}' could not be read from the cells table: {exc}"
        ) from exc
    lab = allb if cell_mask is None else allb[np.asarray(cell_mask, dtype=bool)]
    if len(lab) != n_expected:
        raise ValueError(
            f"batch_key '{batch_key}' gave {len(lab):,} labels for {n_expected:,} cells — "
            f"the cells column and the group labels must describe the same cell set.")
    return np.asarray([("<NA>" if b is None else str(b)) for b in lab])


def _result_from_scores(genexlambda, gene_names, groups_order, group_sizes, n_genes_user, *,
                        use_sparse_output, output_format, iqr_normalize, q_upper, q_lower,
                        iqr_fallback, scale_rule, pvals_matrix=None, zscores_matrix=None,
                        logp_matrix=None, pvalue_fdr_method="fdr_bh", ds=None, modality=None,
                        verbose=False, feature_name_col=None, tag="cytome_cpu", t0=None,
                        t_post=None):
    """The result of an unbatched path (CPU or GPU) from every feature's scores.

    ``genexlambda`` is ``(n_features, n_groups)``, -1 where a filter removed
    the feature. One selector (``cosg._select.top_n``: by score, ties by
    feature index) picks each group's features; ``n_genes_user`` beyond the
    number of features returns them all (it raised a broadcast error).
    """
    n_genes, n_groups = genexlambda.shape
    n_out = min(int(n_genes_user), n_genes)
    if t_post is None:
        t_post = time.time()

    # --- Sparse output (dict/long/dense): extract top-K per group ---
    if use_sparse_output:
        scores_dict = {}
        dense_g, dense_f, dense_v = [], [], []
        collided = set()
        # p-values for exactly the entries that get reported
        _pv_dict, _pv_adj_dict, _nl_dict, _nl_adj_dict = {}, {}, {}, {}
        # the median rule's stand-in for a group with few positive scores (QC doc 52)
        pooled_median = (_pooled_positive_median(genexlambda[:, j] for j in range(n_groups))
                         if iqr_normalize and scale_rule == "positive_median" else None)
        for k in range(n_groups):
            col = genexlambda[:, k].copy()
            col[col < 0] = 0  # COSG sets filtered features to -1

            n_positive = int((col > 0).sum())
            if n_positive == 0:
                continue

            # All-values IQR: includes zero entries (matches
            # cosg.iqrLogNormalize on AnnData side). Clusters with
            # sparser positive signal end up with lower q_lower / q_upper
            # → smaller IQR → larger normalized scores, which is the
            # cross-cluster comparability the AnnData side has.
            if iqr_normalize:
                iqr = _group_scale(col, q_upper, q_lower, iqr_fallback, scale_rule, pooled_median)
            else:
                iqr = None  # raw-λ mode

            # the group's best positive features, ties by index
            n_keep = min(int(n_genes_user), n_positive)
            top_idx = top_n(col, n_keep)

            # Normalize (or keep raw) only the top-K
            if iqr_normalize:
                top_scores = np.log1p(col[top_idx] / iqr).astype(np.float32)
            else:
                top_scores = col[top_idx].astype(np.float32)
            top_names = gene_names[top_idx]
            group_name = str(groups_order[k])

            if pvals_matrix is not None:
                # BH over every gene tested in this group, then read off the
                # reported ones. Adjusting after selection would divide by the
                # size of the top-K rather than the size of the screen. The
                # same in log space, which does not underflow.
                from ._pvalues import adjust_fdr, adjust_fdr_log
                _tested = genexlambda[:, k] > -1
                _adj_col = np.full(pvals_matrix.shape[0], np.nan)
                _adj_col[_tested] = adjust_fdr(
                    pvals_matrix[_tested, k], pvalue_fdr_method)
                _ladj_col = np.full(pvals_matrix.shape[0], np.nan)
                _ladj_col[_tested] = adjust_fdr_log(
                    logp_matrix[_tested, k], pvalue_fdr_method)

            for _j, (name, score) in enumerate(zip(top_names, top_scores)):
                if score > 0:
                    dense_g.append(k); dense_f.append(int(top_idx[_j])); dense_v.append(score)
                    # a name the table holds twice keeps its best score (the
                    # first by the selector's order), not whichever came last
                    if (group_name, str(name)) in scores_dict:
                        collided.add(str(name))
                        continue
                    scores_dict[(group_name, str(name))] = float(score)
                    if pvals_matrix is not None:
                        _gi = int(top_idx[_j])
                        _key = (group_name, str(name))
                        _pv_dict[_key] = float(pvals_matrix[_gi, k])
                        _pv_adj_dict[_key] = float(_adj_col[_gi])
                        _nl_dict[_key] = float(-logp_matrix[_gi, k] / np.log(10.0))
                        _nl_adj_dict[_key] = float(-_ladj_col[_gi] / np.log(10.0))

        groups_order_str = [str(g) for g in groups_order]
        group_size_dict = {
            groups_order_str[k]: int(group_sizes[k])
            for k in range(n_groups)
        }

        if verbose:
            n_groups_with_entries = len({g for (g, _) in scores_dict.keys()})
            avg_per_group = (
                len(scores_dict) / n_groups_with_entries
                if n_groups_with_entries else 0
            )
            print(
                f"[{tag}] n_genes_user={n_genes_user} per group "
                f"(output_format={output_format!r}, "
                f"iqr_normalize={iqr_normalize}): "
                f"{len(scores_dict):,} entries written "
                f"(avg {avg_per_group:,.0f}/group; "
                f"entries with raw COSG score ≤ 0 dropped — "
                f"to include more, lower expressed_pct or mu)"
            )
            print(f"[{tag}] Post-processing: {time.time() - t_post:.3f}s")
            if t0 is not None:
                print(f"[{tag}] TOTAL: {time.time() - t0:.3f}s")

        _warn_repeated_names(collided, output_format)
        # Format dispatch: dict | long | dense
        _result = _format_top_k_output(
            scores_dict, groups_order_str, group_size_dict,
            output_format=output_format,
            ds=ds, modality=modality, verbose=verbose,
            feature_name_col=feature_name_col,
            dense_entries=(dense_g, dense_f, dense_v),
        )
        if pvals_matrix is not None and isinstance(_result, dict):
            _result['pvals'] = _pv_dict
            _result['pvals_adj'] = _pv_adj_dict
            _result['neg_log10_pvals'] = _nl_dict
            _result['neg_log10_pvals_adj'] = _nl_adj_dict
        return _result

    # --- ndarray output: (names, scores) per-group top-N ---
    names_out = np.empty((n_out, n_groups), dtype=object)
    scores_out = np.empty((n_out, n_groups), dtype=np.float32)
    pvals_out = pvals_adj_out = zscores_out = nl_out = nl_adj_out = None
    if pvals_matrix is not None:
        pvals_out = np.empty((n_out, n_groups), dtype=np.float64)
        pvals_adj_out = np.empty((n_out, n_groups), dtype=np.float64)
        zscores_out = np.empty((n_out, n_groups), dtype=np.float64)
        nl_out = np.empty((n_out, n_groups), dtype=np.float64)
        nl_adj_out = np.empty((n_out, n_groups), dtype=np.float64)

    # Pre-compute per-group IQR once if iqr_normalize is requested. Same
    # all-values IQR semantics as the sparse path, so the two output
    # shapes produce consistent normalized scores for the same input.
    if iqr_normalize:
        iqr_per_group = np.empty(n_groups, dtype=np.float32)
        pooled_median = (_pooled_positive_median(genexlambda[:, j] for j in range(n_groups))
                         if scale_rule == "positive_median" else None)
        for k in range(n_groups):
            col_for_iqr = genexlambda[:, k].copy()
            col_for_iqr[col_for_iqr < 0] = 0
            iqr_per_group[k] = _group_scale(col_for_iqr, q_upper, q_lower, iqr_fallback, scale_rule,
                                            pooled_median)

    for k in range(n_groups):
        scores_k = genexlambda[:, k]
        top_idx = top_n(scores_k, n_out)
        names_out[:, k] = gene_names[top_idx]
        raw_top = scores_k[top_idx].astype(np.float32)
        if iqr_normalize:
            # Clamp negatives to 0 before log1p (consistent with sparse path)
            scores_out[:, k] = np.log1p(
                np.maximum(raw_top, 0.0) / iqr_per_group[k]
            )
        else:
            scores_out[:, k] = raw_top

        if pvals_matrix is not None:
            # BH over the genes tested in this group, then subset to the
            # reported top-N -- correcting after selection would divide by
            # the wrong m.
            from ._pvalues import adjust_fdr, adjust_fdr_log
            tested = scores_k > -1
            col = pvals_matrix[:, k]
            adj = np.full(col.shape, np.nan)
            adj[tested] = adjust_fdr(col[tested], pvalue_fdr_method)
            ladj = np.full(col.shape, np.nan)
            ladj[tested] = adjust_fdr_log(logp_matrix[tested, k], pvalue_fdr_method)
            pvals_out[:, k] = col[top_idx]
            pvals_adj_out[:, k] = adj[top_idx]
            zscores_out[:, k] = zscores_matrix[top_idx, k]
            nl_out[:, k] = -logp_matrix[top_idx, k] / np.log(10.0)
            nl_adj_out[:, k] = -ladj[top_idx] / np.log(10.0)

    if verbose:
        print(f"[{tag}] Post-processing: {time.time() - t_post:.3f}s")
        if t0 is not None:
            print(f"[{tag}] TOTAL: {time.time() - t0:.3f}s")

    _nd_result = {
        'names': names_out,
        'scores': scores_out,
        'groups_order': groups_order,
    }
    if pvals_out is not None:
        _nd_result['pvals'] = pvals_out
        _nd_result['pvals_adj'] = pvals_adj_out
        _nd_result['zscores'] = zscores_out
        _nd_result['neg_log10_pvals'] = nl_out
        _nd_result['neg_log10_pvals_adj'] = nl_adj_out
    return _nd_result


def run_cosg_cytome(
    cytome_path,
    groupby,
    n_genes_user=50,
    mu=1.0,
    remove_lowly_expressed=True,
    expressed_pct=0.1,
    expressed_min_num_cells_in_target_group=3,
    cell_mask=None,
    batch_size=2048,
    verbose=True,
    modality="RNA",
    output_format="ndarray",
    iqr_normalize=None,
    q_upper=0.95,
    q_lower=0.75,
    iqr_fallback="positive",
    scale_rule="iqr",
    feature_batching="auto",
    min_total_count=0,
    layer="auto",
    compute_on_fly=True,
    use_cached_stats=True,
    device="cpu",
    calculate_pvalues=False,
    pvalue_method="spa",
    pvalue_fdr_method="fdr_bh",
    write_to_cytome=False,
    key_added="cosg",
    feature_name_col=None,
    batch_key=None,
    batch_cell_number_threshold=3,
    max_accumulator_bytes=None,
    **_deprecated_kwargs,
):
    """
    Run COSG marker gene detection directly from a cytome file.

    Equivalent to ``cosg.cosg(adata, groupby=groupby, ...)`` but reads
    chunks from disk instead of slicing an in-memory AnnData. Peak RAM
    is bounded by one chunk regardless of dataset size.

    Accepts a path string or an open ``cytome.Dataset``. Supports
    ``device='cpu' | 'gpu' | 'auto'`` for in-function GPU dispatch
    (the convenience ``run_cosg_cytome_gpu`` is a thin shim around
    this function).

    Parameters
    ----------
    batch_key : str, optional
        A ``cells`` column to correct for. When set, the cosine specificity is
        computed **within each batch, with that batch's own gene and group
        norms**, then averaged across the batches where the group had at least
        ``batch_cell_number_threshold`` cells. That is what makes the score
        batch-corrected rather than merely batch-weighted, and it matches the
        AnnData path's ``batch_key`` semantics exactly.

        Costs ``n_batches x n_features x n_groups x 4`` bytes of accumulators;
        the function refuses above 16 GB and tells you which knob to turn.
        Groups too small in *every* batch score 0 — they have no evidence — and
        the count is reported when ``verbose``.
    batch_cell_number_threshold : int, default 3
        Minimum cells a group must have **within a batch** for that batch to
        contribute to the group's score.
    cytome_path : str
        Path to ``.cytome`` file.
    groupby : str
        Column name in cell metadata for grouping.
    n_genes_user : int, default 50
        Number of top marker features per group. Matches the
        ``n_genes_user`` parameter on ``cosg.cosg`` (AnnData side).
    mu : float
        COSG penalty parameter.
    remove_lowly_expressed : bool
        Filter genes with low expression percentage.
    expressed_pct : float, dict, pandas.Series or sequence of float, default 0.1
        Minimum fraction of a group's cells that must express a feature for
        it to be one of the group's markers: the threshold is
        ``max(expressed_pct_k * N_k, expressed_min_num_cells_in_target_group)``
        cells. One number for every group, or a fraction per group:

        - a dict ``{group: fraction}`` or a pandas Series indexed by group.
          Every group scored needs a fraction (one left out raises rather
          than silently taking 0.1); a key naming no label of ``groupby``
          raises; a key naming a label that ``cell_mask`` left out is
          ignored, so one table of every type serves any mask. Keys match
          labels as they are or by their string form.
        - a sequence with one fraction per group, in the order of the
          result's ``groups_order``.

        Every fraction lies in [0, 1]. A group's fraction removes features
        from that group's markers only. A fraction of 0 leaves its group the
        minimum of ``expressed_min_num_cells_in_target_group`` cells, as
        ``expressed_pct=0`` does for every group, on every path;
        ``remove_lowly_expressed=False`` switches the filter off. Given per group, the floors are returned under
        ``'expressed_pct'`` (``{group: fraction}``, keyed as
        ``groups_order``) and stored with the markers by ``write_to_cytome``.

        Example: ``expressed_pct={'Microglia': 0.05, 'Astro': 0.1, 'Oligo': 0.1}``.

        A per-type floor read from the data -- the PIASO benchmark harness
        in ``benchmark/piasomarkerdb_chromatin`` (``detection_null.py`` and
        ``capture_floor.py``, run through COSG by ``cosg_pertype_floor.py``)
        writes ``<prefix>_capture.tsv`` with each type's floor as a count of
        nuclei, ``combined_count``: the larger of a background test at FDR
        5 % and 10 % of the type's capture -- is passed as fractions::

            t = pd.read_csv("kidney_l2_capture.tsv", sep="\\t").set_index("cell_type")
            floor = (t["combined_count"] - 0.5) / t["nuclei"]
            run_cosg_cytome(store, groupby="subclass_l2", expressed_pct=floor, ...)

        with ``nuclei`` the cells COSG scores in each type (the same mask);
        the 0.5 makes the threshold "``combined_count`` nuclei or more"
        whatever the rounding.
    expressed_min_num_cells_in_target_group : int
        Absolute minimum nonzero count threshold.
    cell_mask : np.ndarray, optional
        Boolean mask or sorted integer indices for cell selection.
    batch_size : int or None
        Number of rows to vstack from on-disk chunks before yielding.
        Default 2048. ``None`` uses raw 16-row on-disk chunks.
    verbose : bool
        Print progress.
    modality : str
        Modality to use (``"RNA"``, ``"ATAC"``, ``"GA"``, ``"tiles"``).
        Routes via the cytome modality registry to the matching feature
        entity (RNA → ``genes``, ATAC → ``peaks``, GA → ``GA_genes``,
        tiles → ``tiles``) and matrix layer.
    layer : str, default ``"auto"``
        Measurement layer to read for COSG scoring. Resolved by dispatch:

        - ``"auto"`` (**default**): resolve a modality-appropriate normalization —
          **RNA/GA → ``log1p``**, **ATAC/tiles → ``tfidf``**. This is the correct
          per-modality default (COSG scores a cosine specificity, so RNA needs
          log-normalization and ATAC needs TF-IDF; a single fixed default would be
          wrong for one of them). Resolution is idempotent and applied in the
          low-level layer interpreters, so every entry point (cpu/gpu, on-the-fly
          / materialised) resolves ``auto`` identically.
        - ``"log1p"``: on-the-fly ``log1p(x / cell_depth * scale)`` with
          ``scale = 1e4`` (**CP10K → log1p**, matching scanpy's
          ``normalize_total(target_sum=1e4)`` + ``log1p`` — *not* CPM / 1e6). A
          recorded ``{modality}_log1p_params['scale_factor']`` overrides the 1e4
          default. The ``auto`` choice for RNA/GA. Computed by COSG itself: it
          needs no piaso. Each cell's depth comes from the data, never from a
          cache: the row sums of the counts chunk being normalized (on a
          column-major matrix, one pass over its columns first), or
          ``cells.n_fragments`` for ATAC when it is filled.
        - ``"counts"``: **raw counts, used as-is with no transform**
          (``_resolve_chunk_normalizer`` returns ``None``). Needs no piaso,
          like ``log1p``. COSG's cosine specificity is per-feature
          scale-robust so this still works, but it does not correct cell-depth
          bias — use it for pre-normalized inputs, or raw counts by choice.
          (ATAC/tiles callers that want raw tile counts pin this explicitly.)
        - ``"infog"``: INFOG-weighted log1p (an **enhanced RNA normalization**;
          used by the PIASO pipeline and ``inferGRN``). Reads a pre-computed
          ``{modality}_infog`` layer if present, else computes on the fly from
          cached params. Requires **piaso** installed + imported.
        - ``"tfidf"``: TF-IDF (ATAC / tiles only); if the layer is
          missing, TF-IDF stats are loaded from
          ``ds.metadata['{modality}_tfidf_params']`` or computed +
          cached on first use.
        - Any other materialized layer name is read as-is from
          ``ds.{modality}.layer(name)``.
    output_format : str, default ``'ndarray'``
        Output shape:

        - ``'ndarray'`` (default): returns
          ``{'names': np.ndarray, 'scores': np.ndarray, 'groups_order': list}``
          with shape ``(n_genes_user, n_groups)``. Scores are raw COSG
          λ-values unless ``iqr_normalize=True``.
        - ``'dict'``: returns
          ``{'scores_dict': {(group, feature) → score}, 'groups_order': list, 'group_sizes': dict}``.
          Only entries with strictly positive raw COSG λ are included.
        - ``'long'``: returns ``{'scores_df': pd.DataFrame, ...}`` where
          ``scores_df`` has columns ``[group, feature, score, rank]``
          (rank 1 = highest within group). Sorted by group, then rank.
        - ``'dense'``: returns ``{'scores_df': pd.DataFrame, ...}`` where
          ``scores_df`` is indexed by ALL features in the cytome's
          feature table; missing entries fill with ``0.0``. Mirrors
          ``cosg.indexByGene + cosg.iqrLogNormalize`` on AnnData side.

        NOTE on dict size: ``n_genes_user`` IS per-group, but the
        ``'dict'`` / ``'long'`` / ``'dense'`` outputs only include
        entries with strictly positive raw COSG λ. Total size is
        therefore ``≤ n_genes_user × n_groups`` and typically much
        smaller (features below ``expressed_pct`` / negative λ after
        ``mu`` penalty are dropped). To include more entries per
        group, lower ``expressed_pct`` (e.g. 0.0), set
        ``remove_lowly_expressed=False``, or lower ``mu`` (default 1.0).
    iqr_normalize : bool or None, default None
        Whether to apply IQR-log1p normalization (per-group, all-values
        IQR matching ``cosg.iqrLogNormalize`` on the AnnData side).

        - ``None`` (default): True for sparse outputs
          (``dict`` / ``long`` / ``dense``), False for ``ndarray``.
          This matches the historical behaviour (sparse outputs were
          always normalized; ndarray output returned raw λ).
        - ``True``: always normalize (also lets ``ndarray`` mode return
          IQR-log1p-normalized scores).
        - ``False``: always return raw λ (also lets ``dict`` / ``long``
          / ``dense`` skip normalization — useful for raw-score debugging).

        Each group's scale is the spread between its ``q_upper`` and
        ``q_lower`` quantiles over all features, zeros included. Where that
        spread is zero the group keeps a scale of its own: the ``q_upper``
        quantile when it is positive, and otherwise, when the group scores
        above zero at fewer than ``1 - q_upper`` of the features (common on
        peaks and tiles, where a cell type opens a few percent of them), the
        5th percentile of its positive scores, where the spread reaches as
        that share falls to ``1 - q_upper``. ``scale_rule`` chooses another
        rule for sparse features.

        .. versionchanged:: 1.2.1
            The last case used the constant 1e-6, so those groups' normalized
            scores were set by the constant rather than by their data, and
            were larger than other groups' by about ``log(scale / 1e-6)``.
            The scores of groups with a nonzero spread are unchanged.
    use_cached_stats : bool, default True
        For ``infog`` and ``tfidf``, use the statistics PIASO cached in the
        store when they are fresh. ``log1p`` takes each cell's depth from the
        data whatever this is.

        A path is opened read-only, so a call writes nothing into the store
        (statistics PIASO computes for ``infog`` or ``tfidf`` are then
        recomputed by each call, not cached); only ``write_to_cytome=True``
        writes, and only its markers. A Dataset you opened is used as it is:
        opened read-write, PIASO may keep its statistics in it.
    write_to_cytome : bool, default False
        If True, persist the markers to the cytome under
        ``ds.metadata[key_added]`` as a long-form record list
        (``{'groupby', 'modality', 'n_genes_user', 'records': [{group,
        feature, score, rank}, ...]}``). The function still **returns** the
        result dict as usual (persist-and-return) so interactive callers can
        inspect it. Works for any ``output_format``.
    key_added : str, default 'cosg'
        Metadata key under which markers are stored when
        ``write_to_cytome=True`` (mirrors ``cosg.cosg``'s AnnData
        ``key_added`` → ``adata.uns['cosg']``).
    q_upper : float
        Upper quantile for IQR calculation.
    q_lower : float
        Lower quantile for IQR calculation.
    iqr_fallback : {'positive', 'constant'}, default 'positive'
        The scale of a group whose ``q_upper`` quantile is zero (see
        ``iqr_normalize``): ``'positive'``, the 5th percentile of its
        positive scores; ``'constant'``, the 1e-6 of COSG 1.2.0 and earlier,
        to reproduce normalized scores computed with them.
    scale_rule : {'iqr', 'positive_median'}, default 'iqr'
        How each group's scale is set when ``iqr_normalize`` is on.
        ``'iqr'``: the spread over all features above, as
        ``cosg.iqrLogNormalize`` and every earlier release. ``'positive_median'``:
        the median of the group's positive scores, for sparse features
        (peaks, tiles), where most groups score above zero at 5 to 25 % of
        the features and the spread over all features is set by how many
        of them score at all (depth, the expression filter, the size of the
        feature set) more than by the group's own scores. A scale divides a
        group's scores, so it never changes a ranking within a group; it
        matters where scores are compared across groups.

        .. versionadded:: 1.2.1
    feature_batching : str
        Feature batching strategy for ATAC data with millions of features.
        ``"auto"`` selects ``"chromosome"`` when n_features > 100K, else ``"none"``.
        ``"none"`` uses the original single-pass approach (RNA default).
        ``"chromosome"`` batches features by chromosome in a single pass (~950 MB).
        ``"disk"`` does multi-pass processing, one chromosome batch per pass (~200 MB).
        On a row-major matrix ``"chromosome"`` holds every chromosome's sums at
        once only while they fit ``max_accumulator_bytes``; beyond it the
        chromosomes are packed into passes that fit, so memory stops growing with
        the number of features times groups. The result is the same either way.
    max_accumulator_bytes : int or None
        The most memory the per-feature sums of ``"chromosome"`` batching may hold
        in one pass over the cells (default ``None``: 2 GiB). A column-major
        matrix is scored a block of features at a time and never needs it.
    min_total_count : int
        Drop features that are nonzero in this many cells or fewer, all groups
        together, before scoring. Applied when feature_batching != "none",
        whatever ``expressed_pct`` and ``remove_lowly_expressed`` are.

    Returns
    -------
    dict
        Shape depends on ``output_format`` (see above).

    Migration notes
    ---------------
    The old ``top_k`` and ``n_top_genes`` kwargs were removed. Passing
    them raises ``TypeError`` with an actionable migration message:

        - ``top_k=N``        → ``n_genes_user=N, output_format='dict'``
        - ``n_top_genes=N``  → ``n_genes_user=N``

    """
    t0 = time.time()
    _check_iqr_fallback(iqr_fallback)
    _check_scale_rule(scale_rule)
    _check_scoring_args(mu, n_genes_user, expressed_min_num_cells_in_target_group)
    n_genes_user = int(n_genes_user)
    # the form and range of the floor now; a dict or Series is matched to the
    # groups once their labels are read, before any chunk is streamed
    check_expressed_pct(expressed_pct)

    use_sparse_output, iqr_normalize = _validate_cosg_api_v2(
        _deprecated_kwargs, output_format, iqr_normalize,
        fn_name="run_cosg_cytome",
    )

    # Round 8 (2026-05-23): accept cytome.Dataset OR path. We only close
    # when we opened it ourselves, so callers who pass an open Dataset
    # retain ownership of its lifecycle.
    ds, _opened_here = _open_or_use(cytome_path)
    try:
        cell_mask = _normalize_cell_mask(cell_mask, int(ds.n_cells))
        if write_to_cytome:
            # a store that cannot be written fails now, not after the run
            _require_writable(cytome_path if _opened_here else ds)
    except Exception:
        if _opened_here:
            ds.close()
        raise

    # device='gpu' / 'auto': the GPU accumulates and scores, and the CPU
    # path's own post-processing makes the result, so every output format,
    # iqr_normalize and the scale options work there. It cannot compute
    # p-values (they need power sums and a second pass), batch the features,
    # or read a column-major matrix. Those requests used to return a result
    # without the p-values or the normalization asked for; now device='auto'
    # runs them on the CPU and device='gpu' raises.
    from ._backend import get_device
    _resolved_device = get_device(device)
    if _resolved_device == "gpu":
        _why = None
        if calculate_pvalues:
            _why = "calculate_pvalues"
        elif feature_batching not in (None, "none", "auto"):
            _why = f"feature_batching={feature_batching!r}"
        elif min_total_count and min_total_count > 0:
            # applied by the feature-batched CPU paths only; the GPU path is
            # unbatched and would return other markers without a word
            _why = "min_total_count"
        _gpu_out = None
        if _why is None:
            try:
                _gpu_out = _run_cosg_cytome_gpu_impl(
                    ds, groupby, mu=mu,
                    remove_lowly_expressed=remove_lowly_expressed,
                    expressed_pct=expressed_pct,
                    expressed_min_num_cells_in_target_group=
                        expressed_min_num_cells_in_target_group,
                    cell_mask=cell_mask, batch_size=batch_size, verbose=verbose,
                    modality=modality, layer=layer,
                    compute_on_fly=compute_on_fly,
                    use_cached_stats=use_cached_stats,
                    feature_name_col=feature_name_col,
                    batch_key=batch_key,
                    batch_cell_number_threshold=batch_cell_number_threshold,
                )
            except _GPUUnsupported as exc:
                _why = str(exc)
            except BaseException:
                if _opened_here:
                    ds.close()
                raise
        if _gpu_out is not None:
            genexlambda, gene_names, groups_order, group_sizes, floors, _floor_per_group = _gpu_out
            try:
                _gpu_result = _result_from_scores(
                    genexlambda, gene_names, groups_order, group_sizes, n_genes_user,
                    use_sparse_output=use_sparse_output, output_format=output_format,
                    iqr_normalize=iqr_normalize, q_upper=q_upper, q_lower=q_lower,
                    iqr_fallback=iqr_fallback, scale_rule=scale_rule,
                    ds=ds, modality=modality, verbose=verbose,
                    feature_name_col=feature_name_col, tag="cytome_gpu", t0=t0)
                if _floor_per_group:
                    _gpu_result['expressed_pct'] = floor_record(
                        [str(g) for g in groups_order] if use_sparse_output else groups_order, floors)
            finally:
                if _opened_here:
                    ds.close()
            if write_to_cytome:
                _persist_cosg_markers(cytome_path, _gpu_result, groupby,
                                      modality, n_genes_user, key_added)
            return _gpu_result
        if device == "gpu":
            if _opened_here:
                ds.close()
            raise NotImplementedError(
                f"device='gpu' does not support {_why}. Use device='cpu' (or "
                f"device='auto', which runs such a request on the CPU).")
        if verbose:
            print(f"[cosg] device='auto': the GPU path does not support {_why}; "
                  f"running on the CPU.")

    # Resolve feature_batching="auto"
    if feature_batching == "auto":
        if modality in ("tiles", "ATAC"):
            table = "tiles" if modality == "tiles" else "peaks"
            try:
                count_row = ds._conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()
                n_feat_check = count_row[0] if count_row else 0
            except Exception:
                n_feat_check = 0
            feature_batching = "chromosome" if n_feat_check > 100_000 else "none"
        else:
            feature_batching = "none"

    # The unbatched path streams row chunks, and a column-major matrix -- the
    # importer's and the quantifier's default for peaks and tiles -- has none.
    # The feature-batched path walks columns (`_scan_column_major`),
    # and its result is identical; only the pass structure differs. So send
    # column-major input there rather than teaching a second reader the layout.
    # Resolve the layer to read ONCE, up front (read a materialised layer, or
    # read counts and transform on the fly). The layout check needs the
    # matrix actually read: with layer='auto' it looked for '<modality>_auto',
    # found nothing, and a column-major matrix then crashed in iter_chunks.
    _layer_to_read, _need_normalizer = _resolve_layer_to_read(
        ds, modality, layer, compute_on_fly,
    )
    if feature_batching == "none" and _column_major_layer(ds, modality, _layer_to_read) is not None:
        if calculate_pvalues:
            if _opened_here:
                ds.close()
            raise NotImplementedError(
                f"calculate_pvalues needs the unbatched path, which reads row "
                f"chunks, and '{modality}_{_layer_to_read}' is stored "
                f"column-major so it has none. Re-import or re-quantify with "
                f"--layout csr for a significance run, or run without "
                f"calculate_pvalues.")
        if verbose:
            print("[cosg] matrix is column-major; using the feature-batched "
                  "path, which reads columns. The result is identical, only "
                  "the pass structure differs.")
        feature_batching = "chromosome"

    # Dispatch to chromosome-batched path for ATAC.
    # Round 8: only close if we opened it ourselves — Dataset input keeps
    # its lifecycle. _run_cosg_batched accepts both Path and Dataset via
    # _open_or_use too.
    if feature_batching in ("chromosome", "disk"):
        _batch_labels_fb = None
        if batch_key is not None:
            _batch_labels_fb = _read_batch_labels(ds, batch_key, cell_mask,
                                                  n_expected=len(_get_labels_and_gene_names(
                                                      ds, groupby, cell_mask, modality=modality,
                                                      load_gene_names=False)[0]))
            # Every accumulator gains a batch axis, so 'chromosome' (single pass,
            # all chromosomes resident) would hold n_batches x its usual RAM.
            # 'disk' splits chromosomes to a byte budget that now accounts for the
            # batch axis, keeping the target the same.
            if feature_batching == "chromosome":
                if verbose:
                    print("[cosg] batch_key: switching feature_batching "
                          "'chromosome' -> 'disk' so the batch axis does not multiply "
                          "peak RAM; the result is identical, only the pass structure "
                          "differs.")
                feature_batching = "disk"
        if _opened_here:
            ds.close()
        if calculate_pvalues:
            # This path holds only a block of features at a time, so the
            # second pass that refines a tail cannot reach the values it
            # needs. Say so rather than returning a table without the
            # columns that were asked for.
            raise NotImplementedError(
                "calculate_pvalues is not available on the feature-batched "
                "path (feature_batching=%r). Pass feature_batching='none' "
                "for a significance run, or use cosg.cosg() on an AnnData."
                % (feature_batching,))
        _batched_result = _run_cosg_batched(
            cytome_path, groupby, n_genes_user, mu,
            remove_lowly_expressed, expressed_pct,
            expressed_min_num_cells_in_target_group,
            cell_mask, batch_size, verbose, modality,
            use_sparse_output, q_upper, q_lower,
            min_total_count, feature_batching,
            layer=layer,
            compute_on_fly=compute_on_fly,
            use_cached_stats=use_cached_stats,
            iqr_normalize=iqr_normalize,
            iqr_fallback=iqr_fallback,
            scale_rule=scale_rule,
            output_format=output_format,
            feature_name_col=feature_name_col,
            batch_labels=_batch_labels_fb,
            batch_cell_number_threshold=batch_cell_number_threshold,
            max_accumulator_bytes=max_accumulator_bytes,
        )
        if write_to_cytome:
            _persist_cosg_markers(cytome_path, _batched_result, groupby,
                                  modality, n_genes_user, key_added)
        return _batched_result

    # Original "none" path (RNA default). The normalizer is resolved here, not
    # before the dispatch above: the batched path resolves its own, and on a
    # store opened read-only PIASO's statistics would be computed twice.
    _chunk_normalizer = (
        _resolve_chunk_normalizer(layer, ds, modality, use_cached_stats)
        if _need_normalizer else None
    )
    labels, gene_names, all_labels, categories = _get_labels_and_gene_names(
        ds, groupby, cell_mask, modality=modality, feature_name_col=feature_name_col
    )

    # Batch labels, aligned to the same cells as `labels`.
    batch_labels = (None if batch_key is None
                    else _read_batch_labels(ds, batch_key, cell_mask, len(labels)))
    groups_order, cell_group_indices, group_sizes = _build_group_mapping(
        labels, categories
    )
    floors, _floor_per_group = _resolve_floors(expressed_pct, groups_order, group_sizes,
                                               all_labels, categories)

    n_genes = len(gene_names)
    n_groups = len(groups_order)

    if verbose:
        print(f"[cytome_cpu] {len(labels):,} cells × {n_genes:,} genes, "
              f"{n_groups} groups")
        _emit_filter_log(
            "cytome_cpu",
            remove_lowly_expressed=remove_lowly_expressed,
            expressed_pct=describe_floor(floors, _floor_per_group, expressed_pct),
            expressed_min_num_cells_in_target_group=
                expressed_min_num_cells_in_target_group,
            mu=mu,
        )

    _check_feature_count(ds, modality, n_genes, feature_name_col)

    # ── batch-aware accumulators ────────────────────────────────────────────
    # Without batch_key: one (genes x groups) dot accumulator, as before.
    # With batch_key: one per batch, so each batch's cosine similarity can be
    # formed with ITS OWN gene/group norms before averaging — that is what makes
    # the score batch-corrected rather than merely batch-weighted. Mirrors the
    # AnnData path (cosg.py, the `batch_key is not None and cpu_chunk_size != 0`
    # branch) exactly: mask groups below the per-batch cell threshold, sum, then
    # divide by the number of batches in which the group was valid.
    if batch_labels is None:
        n_batches = 1
        batch_index = None
        batches_order = None
    else:
        batches_order, batch_index = np.unique(batch_labels, return_inverse=True)
        n_batches = len(batches_order)
        need_gb = n_batches * n_genes * n_groups * 8 / 1e9
        if verbose:
            print(f"[cytome_cpu] batch_key: {n_batches} batches x {n_genes:,} features "
                  f"x {n_groups} groups -> {need_gb:.2f} GB of accumulators")
        if need_gb > 16.0:
            raise MemoryError(
                f"batch-aware COSG needs {need_gb:.1f} GB of accumulators "
                f"({n_batches} batches x {n_genes:,} features x {n_groups} groups x 8 B). "
                f"Reduce the feature space (min_total_count), coarsen `groupby`, or use "
                f"feature_batching='chrom' to process features in blocks.")

    # Sums in float64, as on the AnnData paths: genes x groups is small here
    # (the feature-batched paths, built for millions of peaks or tiles, keep
    # float32 sums to bound their memory). In float32 the sums over a large
    # store lost digits chunk after chunk.
    dot_accum = np.zeros((n_batches, n_genes, n_groups), dtype=np.float64)
    gene_sq_accum = np.zeros((n_batches, n_genes), dtype=np.float64)
    grp_count_accum = np.zeros((n_batches, n_groups), dtype=np.float64)

    # The analytic p-value needs the first three population moments per gene.
    # Two of the three quantities are already here: dot_accum IS the group sum
    # T, and gene_sq_accum is the second power sum. Only the first and third
    # power sums are new, so significance costs two reductions in a pass that
    # already runs -- which is why the streaming path can carry it at all.
    if calculate_pvalues:
        gene_sum_accum = np.zeros(n_genes, dtype=np.float64)
        gene_cube_accum = np.zeros(n_genes, dtype=np.float64)
        gene_nnz_accum = np.zeros(n_genes, dtype=np.int64)
    else:
        gene_sum_accum = gene_cube_accum = gene_nnz_accum = None

    if tracks_floor(remove_lowly_expressed, floors, expressed_min_num_cells_in_target_group):
        nnz_accum = np.zeros((n_genes, n_groups), dtype=np.int32)
    else:
        nnz_accum = None

    # Stream chunks from cytome disk storage
    n_chunks = 0
    n_cells_seen = 0
    t_compute = 0.0

    for chunk_csr, row_indices in ds.iter_chunks(
        modality, cell_mask=cell_mask, batch_size=batch_size,
        layer=_layer_to_read,
    ):
        # On-the-fly per-chunk normalization (e.g. infog / tfidf / log1p
        # when the materialised layer doesn't exist).
        if _chunk_normalizer is not None:
            chunk_csr = _chunk_normalizer(chunk_csr, row_indices)
        t_comp = time.time()

        actual_chunk = chunk_csr.shape[0]

        # Map row_indices to positions in the labels array
        if cell_mask is not None:
            chunk_labels_idx = np.arange(
                n_cells_seen, n_cells_seen + actual_chunk, dtype=np.int32
            )
        else:
            chunk_labels_idx = row_indices

        chunk_groups = cell_group_indices[chunk_labels_idx]

        # Build indicator matrix for this chunk
        Lam_chunk = sp_sparse.csr_matrix(
            (np.ones(actual_chunk, dtype=np.float64),
             (np.arange(actual_chunk, dtype=np.int32), chunk_groups)),
            shape=(actual_chunk, n_groups)
        )

        # The chunk's values in float64 for the sums (the products and the
        # squares as well as their accumulation)
        chunk_csr = chunk_csr.astype(np.float64)

        # Accumulate dot products: (genes × groups), per batch
        if batch_index is None:
            D_chunk = chunk_csr.T @ Lam_chunk
            if sp_sparse.issparse(D_chunk):
                D_chunk = D_chunk.toarray()
            dot_accum[0] += D_chunk

            Xsq = chunk_csr.copy()
            Xsq.data = Xsq.data ** 2
            gene_sq_accum[0] += np.asarray(Xsq.sum(axis=0)).ravel()
            if calculate_pvalues:
                gene_sum_accum += np.asarray(chunk_csr.sum(axis=0)).ravel()
                _cube = chunk_csr.copy()
                _cube.data = _cube.data ** 3
                gene_cube_accum += np.asarray(_cube.sum(axis=0)).ravel()
                gene_nnz_accum += np.bincount(chunk_csr.indices[chunk_csr.data != 0],
                                              minlength=n_genes)
            del Xsq
            np.add.at(grp_count_accum[0], chunk_groups, 1.0)
        else:
            chunk_batches = batch_index[chunk_labels_idx]
            D_chunk = None
            for b in np.unique(chunk_batches):
                sel = np.flatnonzero(chunk_batches == b)
                sub = chunk_csr[sel]
                Lam_sub = Lam_chunk[sel]
                Db = sub.T @ Lam_sub
                if sp_sparse.issparse(Db):
                    Db = Db.toarray()
                dot_accum[b] += Db
                Xsq = sub.copy()
                Xsq.data = Xsq.data ** 2
                gene_sq_accum[b] += np.asarray(Xsq.sum(axis=0)).ravel()
                np.add.at(grp_count_accum[b], chunk_groups[sel], 1.0)
                del Xsq, sub, Lam_sub, Db

        # Accumulate nonzero counts
        if nnz_accum is not None:
            # "expressed" means a nonzero value: a stored zero is not counted
            X_bin = chunk_csr.copy()
            X_bin.data = (X_bin.data != 0).astype(X_bin.dtype)
            nnz_chunk = X_bin.T @ Lam_chunk
            if sp_sparse.issparse(nnz_chunk):
                nnz_chunk = nnz_chunk.toarray()
            nnz_accum += nnz_chunk.astype(np.int32)
            del X_bin, nnz_chunk

        t_compute += time.time() - t_comp

        n_chunks += 1
        n_cells_seen += actual_chunk

        del chunk_csr, D_chunk, Lam_chunk

        if verbose and n_chunks % 500 == 0:
            print(f"  Processed {n_chunks} chunks, {n_cells_seen:,} cells")

    if verbose:
        print(f"[cytome_cpu] Streaming: {n_chunks} chunks, "
              f"{n_cells_seen:,} cells, {t_compute:.3f}s compute")

    # a NaN or an infinity in a feature makes its sum of squares non-finite
    from ._validate import check_finite_sums
    check_finite_sums(gene_sq_accum, gene_names, what=_matrix_label(modality, _layer_to_read, layer,
                                                                      _chunk_normalizer))

    t_post = time.time()
    if batch_index is not None:
        # With batches the cosine is formed per batch with that batch's own
        # norms and averaged over the batches where the group has enough
        # cells (score_from_sums); a group valid in no batch scores 0.
        n_dropped = int(((grp_count_accum >= batch_cell_number_threshold).sum(axis=0) == 0).sum())
        if verbose and n_dropped:
            print(f"[cytome_cpu] batch_key: {n_dropped}/{n_groups} group(s) had "
                  f"< {batch_cell_number_threshold} cells in every batch; scored 0")

    pvals_matrix = zscores_matrix = logp_matrix = None
    if calculate_pvalues:
        from ._pvalues import srswor_moments_from_sums, normal_tail, _log_binom

        if n_batches > 1:
            raise NotImplementedError(
                "calculate_pvalues with batch_key is not supported on the "
                "streaming path yet; the stratified null needs per-batch "
                "power sums. Use cosg.cosg() on an AnnData for a batched "
                "significance run, or drop batch_key here.")

        n_tot = int(n_cells_seen)
        T_all = dot_accum[0].astype(np.float64)          # genes x groups
        nc_all = grp_count_accum[0].astype(np.float64)
        mean_g, s2_g, m3_g = srswor_moments_from_sums(
            gene_sum_accum, gene_sq_accum[0].astype(np.float64),
            gene_cube_accum, n_tot)

        # the log p as well, from log_ndtr: a strong marker's p underflows to
        # 0, and -log10 p keeps ranking it (the AnnData path's neg_log10_pvals)
        pvals_matrix, zscores_matrix, flagged, logp_matrix = normal_tail(
            T_all, nc_all, n_tot, mean_g, s2_g, m3_g,
            k_nonzero=gene_nnz_accum,
            method=pvalue_method,
            return_log=True,
        )

        # Second pass, over the flagged genes only, accumulating their value
        # atoms rather than their values. The grid is fixed, so counts and
        # sums per (gene, bin) add across chunks: the streamed summary is
        # *identical* to the in-memory one, not merely close, and memory is
        # bounded by the number of flagged genes instead of by how many cells
        # express them.
        if flagged.size and pvalue_method in ("spa", "exact"):
            from ._pvalue_atoms import (StreamingAtomAccumulator, spa_atoms_tail,
                                        values_are_integral, exact_tail_bounded)

            want = np.unique(flagged[:, 0])
            if verbose:
                print(f"  p-values: refining {flagged.shape[0]} gene x group "
                      f"entries over {want.size} genes (second pass)")
            acc_atoms = StreamingAtomAccumulator(want)
            integral = True
            for chunk_csr, _ri in ds.iter_chunks(
                    modality, cell_mask=cell_mask, batch_size=batch_size,
                    layer=_layer_to_read):
                if _chunk_normalizer is not None:
                    chunk_csr = _chunk_normalizer(chunk_csr, _ri)
                if integral and not values_are_integral(chunk_csr.data):
                    integral = False
                acc_atoms.add_chunk(chunk_csr)
            aset = acc_atoms.finish()
            pmap = aset.index_of()

            exact_done = set()
            if pvalue_method == "exact" and not integral:
                import warnings as _w
                _w.warn(
                    "pvalue_method='exact' does not apply here: the values are not "
                    "integers -- the exact tail indexes its state by the sum, so the "
                    "sum has to lie on a lattice. The saddlepoint is used instead.",
                    stacklevel=2)
            if pvalue_method == "exact" and integral:
                for gi, gj in flagged:
                    gi, gj = int(gi), int(gj)
                    j = pmap.get(gi)
                    if j is None:
                        continue
                    lo, hi = aset.ptr[j], aset.ptr[j + 1]
                    av, aw = aset.vals[lo:hi], aset.wts[lo:hi]
                    k = float(aw.sum()); nc = float(nc_all[gj])
                    if k == 0 or nc <= 0 or nc >= n_tot:
                        continue
                    eh = k * nc / max(n_tot, 1)
                    h_max = int(min(k, nc, eh + 10.0 * np.sqrt(max(eh, 1.0)) + 25))
                    # None where the table's cut at h_max could matter (a
                    # strong marker came back as exactly 0): the saddlepoint
                    # takes those entries
                    got = exact_tail_bounded(av, aw, n_tot, int(nc),
                                             float(T_all[gi, gj]), h_max)
                    if got is not None:
                        pvals_matrix[gi, gj] = got
                        logp_matrix[gi, gj] = np.log(got)
                        exact_done.add((gi, gj))

            # One solver, shared with the in-memory path. Two implementations
            # of one tail drift, and the drift shows up as a parity failure
            # that reads as a numerical accident rather than as two different
            # fallback decisions.
            pairs, vals_parts, wts_parts = [], [], []
            blk_ptr, blk_row, n_zero_l, n_c_l, t_obs = [0], [], [], [], []
            pos = 0
            for gi, gj in flagged:
                gi, gj = int(gi), int(gj)
                if (gi, gj) in exact_done:
                    continue
                j = pmap.get(gi)
                if j is None:
                    continue
                lo, hi = aset.ptr[j], aset.ptr[j + 1]
                av, aw = aset.vals[lo:hi], aset.wts[lo:hi]
                nc = float(nc_all[gj])
                if av.size == 0 or nc <= 0 or nc >= n_tot:
                    continue
                pairs.append((gi, gj))
                vals_parts.append(av); wts_parts.append(aw)
                pos += av.size
                blk_ptr.append(pos); blk_row.append(len(pairs) - 1)
                n_zero_l.append(float(n_tot) - float(aw.sum()))
                n_c_l.append(nc)
                t_obs.append(float(T_all[gi, gj]))
            if pairs:
                pvec, lpvec, okvec = spa_atoms_tail(
                    np.concatenate(vals_parts), np.concatenate(wts_parts),
                    np.asarray(blk_ptr, dtype=np.int64),
                    np.asarray(blk_row, dtype=np.int64),
                    np.asarray(n_zero_l, dtype=np.float64),
                    np.asarray(n_c_l, dtype=np.float64),
                    np.asarray(t_obs, dtype=np.float64),
                    lattice=integral)
                failed = []
                for r, (gi, gj) in enumerate(pairs):
                    if okvec[r]:
                        pvals_matrix[gi, gj] = pvec[r]
                        logp_matrix[gi, gj] = lpvec[r]
                    else:
                        failed.append((gi, gj))
                # a failure at the edge of the null's support (the group holds
                # the gene's top values) is taken by the permutation floor
                # below; only the others are reported, as on the AnnData path
                unexplained = sum(1 for gi, gj in failed
                                  if logp_matrix[gi, gj] >= -_log_binom(n_tot, nc_all[gj]))
                if unexplained:
                    import warnings as _w
                    _w.warn(
                        f"COSG p-values: {unexplained} of {len(pairs)} refined "
                        "entries did not converge and kept their "
                        "normal-approximation p-value, which is "
                        "anti-conservative in exactly that regime.",
                        stacklevel=2)

        # as the AnnData path: no p below the permutation floor (one labeling
        # of the group's cells), p floored at the smallest positive float64 (a
        # strong marker's tail underflows to 0), log p capped at 0
        from ._pvalues import _at_least_the_floor, _log_binom as _lb
        _floor = np.array([-_lb(n_tot, c) if 0 < c < n_tot else 0.0 for c in nc_all])
        pvals_matrix, zscores_matrix, _pinfo = _at_least_the_floor(
            pvals_matrix, zscores_matrix, {}, logp_matrix, _floor)
        logp_matrix = _pinfo["log_p"]

    # The cosine, COSG's penalty and the expression filter, from the sums:
    # each group's own threshold on its own column (none for a group without
    # cells; a floor of 0 keeps the minimum of cells).
    genexlambda = _score_from_sums(
        dot_accum[0] if batch_index is None else dot_accum,
        gene_sq_accum[0] if batch_index is None else gene_sq_accum,
        group_sizes, nnz_accum, mu=mu,
        thresholds=(None if nnz_accum is None else floor_thresholds(
            group_sizes, floors, expressed_min_num_cells_in_target_group,
            skip_empty=True)),
        batch_group_sizes=None if batch_index is None else grp_count_accum,
        batch_cell_number_threshold=batch_cell_number_threshold)
    del dot_accum, gene_sq_accum, nnz_accum

    _result = _result_from_scores(
        genexlambda, gene_names, groups_order, group_sizes, n_genes_user,
        use_sparse_output=use_sparse_output, output_format=output_format,
        iqr_normalize=iqr_normalize, q_upper=q_upper, q_lower=q_lower,
        iqr_fallback=iqr_fallback, scale_rule=scale_rule,
        pvals_matrix=pvals_matrix, zscores_matrix=zscores_matrix, logp_matrix=logp_matrix,
        pvalue_fdr_method=pvalue_fdr_method, ds=ds, modality=modality, verbose=verbose,
        feature_name_col=feature_name_col, tag="cytome_cpu", t0=t0, t_post=t_post)
    del genexlambda
    if _floor_per_group:
        _result['expressed_pct'] = floor_record(
            [str(g) for g in groups_order] if use_sparse_output else groups_order, floors)
    if _opened_here:
        ds.close()
    if write_to_cytome:
        _persist_cosg_markers(cytome_path, _result, groupby,
                              modality, n_genes_user, key_added)
    return _result


# ===================================================================
#  K-HARD migration stub — Round 9 (2026-05-23)
# ===================================================================

def run_cosg_cytome_cpu(*args, **kwargs):
    """REMOVED — renamed to :func:`run_cosg_cytome`.

    Stub function: raises ``TypeError`` at call time with a copy/pasteable
    migration hint. Kept here (rather than removed entirely) so callers
    using ``from cosg._cytome_streaming import run_cosg_cytome_cpu`` or
    ``cosg.run_cosg_cytome_cpu`` see a precise error pointing at their
    call site, with the new function name in the message.

    Implementation note: a stub function failing at call time is friendlier
    than a module ``__getattr__`` that fails at import time (which would
    break linters, completion engines, and any tool that just imports
    names for introspection). The traceback points exactly at the user's
    call, with the migration message right there.
    """
    raise TypeError(
        "cosg.run_cosg_cytome_cpu has been renamed to cosg.run_cosg_cytome. "
        "The renamed function accepts device='cpu' | 'gpu' | 'auto'. "
        "Migration: replace `cosg.run_cosg_cytome_cpu(...)` with "
        "`cosg.run_cosg_cytome(...)`. The kwarg `cytome_layer` was also "
        "renamed to `layer`."
    )


# ===================================================================
#  Chromosome-batched COSG for ATAC (Strategy 4 + Strategy 3)
# ===================================================================

def _run_cosg_batched(
    cytome_path, groupby, n_genes_user, mu,
    remove_lowly_expressed, expressed_pct,
    expressed_min_num_cells_in_target_group,
    cell_mask, batch_size, verbose, modality,
    use_sparse_output, q_upper, q_lower,
    min_total_count, mode,
    layer="auto",
    compute_on_fly=True,
    use_cached_stats=True,
    iqr_normalize=False,
    output_format="ndarray",
    feature_name_col=None,
    batch_labels=None,
    batch_cell_number_threshold=3,
    block_features=None,
    iqr_fallback="positive",
    scale_rule="iqr",
    max_accumulator_bytes=None,
):
    """Chromosome-batched COSG for ATAC data with millions of features.

    A column-major matrix (the importer's and the quantifier's default for
    peaks and tiles) is not batched by chromosome at all: it is scored as it
    is read, a block of features at a time (``_scan_column_major``), so its
    memory does not grow with the number of features or chromosomes.
    ``block_features`` fixes the block for tests; it never changes the result.

    Processes features in per-chromosome batches instead of all at once,
    reducing peak RAM from ~3.7 GB to ~950 MB (chromosome mode) or ~200 MB
    (disk/multi-pass mode).  Gene names are loaded lazily — only the final
    top-K feature names are queried from SQLite.

    Internal helper; the public API is ``run_cosg_cytome_cpu``. Callers
    pass ``use_sparse_output`` (bool) and ``iqr_normalize`` (bool, already
    resolved from the public sentinel default).

    Parameters
    ----------
    mode : str
        ``"chromosome"`` — single pass, all chrom accumulators in RAM.
        ``"disk"`` — multi-pass, one chromosome batch per pass (lower RAM).
    """
    t0 = time.time()
    # Round 8: accept Dataset OR path. Caller owns lifecycle if it passed
    # an open Dataset; we only close at our exit if we opened it.
    ds, _opened_here = _open_or_use(cytome_path)

    # Get labels only (lazy gene_names — saves ~400 MB)
    labels, _, all_labels, categories = _get_labels_and_gene_names(
        ds, groupby, cell_mask, modality=modality, load_gene_names=False
    )
    groups_order, cell_group_indices, group_sizes = _build_group_mapping(
        labels, categories
    )
    n_groups = len(groups_order)
    floors, _floor_per_group = _resolve_floors(expressed_pct, groups_order, group_sizes,
                                               all_labels, categories)

    # Batch dimension. Each chromosome's accumulators gain a leading batch axis
    # so the cosine can be formed per batch with that batch's own norms — the
    # same construction as the non-feature-batched path, just applied within a
    # chromosome block. `grp_count_accum` is global (a cell's group and batch do
    # not depend on which features we happen to be holding).
    if batch_labels is None:
        n_batches, batch_index = 1, None
    else:
        _bo, batch_index = np.unique(np.asarray(batch_labels), return_inverse=True)
        n_batches = len(_bo)
    grp_count_accum = np.zeros((n_batches, n_groups), dtype=np.float64)
    if batch_index is not None:
        for gi, bi in zip(cell_group_indices, batch_index):
            grp_count_accum[bi, gi] += 1.0
    else:
        grp_count_accum[0] = group_sizes

    # Get chromosome → column ranges
    feature_table, idx_col, name_col = _get_feature_table_info(ds, modality, feature_name_col)
    chrom_ranges, chrom_order = _get_chrom_col_ranges(ds, feature_table, idx_col)
    n_genes = sum(ce - cs for cs, ce in chrom_ranges.values())
    _check_feature_count(ds, modality, n_genes, feature_name_col)

    # Resolve layer + normalizer once for the entire batched accumulation
    _layer_to_read, _need_normalizer = _resolve_layer_to_read(
        ds, modality, layer, compute_on_fly,
    )
    _chunk_normalizer = (
        _resolve_chunk_normalizer(layer, ds, modality, use_cached_stats)
        if _need_normalizer else None
    )

    # The expression filter: each group's own threshold, max(N_k * pct_k,
    # min_cells), none for a group without cells (a floor of 0 keeps the
    # minimum of cells), and no filter with remove_lowly_expressed=False.
    # The nonzero counts are kept when the filter or min_total_count needs
    # them: min_total_count applies whatever the floor and the filter.
    _filter_on = tracks_floor(remove_lowly_expressed, floors, expressed_min_num_cells_in_target_group)
    track_nnz = _filter_on or min_total_count > 0
    thresholds = (floor_thresholds(group_sizes, floors, expressed_min_num_cells_in_target_group,
                                   skip_empty=True) if _filter_on else None)

    if verbose:
        print(f"[cosg_{mode}] {len(labels):,} cells × {n_genes:,} features, "
              f"{n_groups} groups, {len(chrom_order)} chromosomes")
        _emit_filter_log(
            f"cosg_{mode}",
            remove_lowly_expressed=remove_lowly_expressed,
            expressed_pct=describe_floor(floors, _floor_per_group, expressed_pct),
            expressed_min_num_cells_in_target_group=
                expressed_min_num_cells_in_target_group,
            mu=mu,
        )

    # Determine chromosome batches based on mode: "disk" packs the chromosomes into
    # passes of ~200 MB; "chromosome" holds them all in one pass while their sums fit
    # the accumulator budget, and packs them into passes of that budget beyond it (a
    # whole-brain atlas's tiles by its subclasses would otherwise hold some 13 GB).
    # Each feature's sums are made the same way in either case, so the result is too.
    budget = (_ACCUM_BUDGET_BYTES if max_accumulator_bytes is None else int(max_accumulator_bytes))
    chrom_batches = _chrom_passes(chrom_order, chrom_ranges, n_groups, n_batches, track_nnz,
                                  200 * 1024 * 1024 if mode != "chromosome" else budget)
    if verbose and (mode != "chromosome" or len(chrom_batches) > 1):
        print(f"[cosg_{mode}] Split into {len(chrom_batches)} passes"
              + (f" (the sums of every chromosome would exceed {budget / 2**30:.1f} GiB)"
                 if mode == "chromosome" else ""))

    # A NaN or an infinity in a feature makes its sum of squares non-finite:
    # checked per block, named through the feature table
    from ._validate import check_finite_sums
    _what = _matrix_label(modality, _layer_to_read, layer, _chunk_normalizer)

    def _names_of(idx):
        return _lookup_feature_names(ds, feature_table, idx_col, name_col, idx)

    # Output collectors. n_keep is the per-group cap on candidate features,
    # n_genes_user (and never more than the features). The ndarray output takes
    # a group's best features of any score, as the unbatched path does (its
    # tail held "" and 0.0 here); the other outputs report positive scores
    # only, so only those are candidates. Candidates are chosen and merged by
    # one rule, cosg._select's: by score, ties by feature index.
    n_keep = min(int(n_genes_user), int(n_genes))
    candidates_all = not use_sparse_output
    cand_idx = {k: [] for k in range(n_groups)}
    cand_val = {k: [] for k in range(n_groups)}
    group_candidates = None

    def _trim(k):
        idx = np.concatenate(cand_idx[k]); val = np.concatenate(cand_val[k])
        keep = top_n_of(val, idx, n_keep)
        cand_idx[k], cand_val[k] = [idx[keep]], [val[keep]]
    # iqr_collector is needed whenever IQR normalization will be applied.
    iqr_collector = {k: [] for k in range(n_groups)} if iqr_normalize else None

    t_stream_total = 0.0

    # A column-major matrix is scored as it is read, a block of features at a
    # time: no accumulator for the whole feature axis, so no chromosome
    # batches either, in either mode.
    _cm_layer = _column_major_layer(ds, modality, _layer_to_read)
    if _cm_layer is not None:
        if cell_mask is None:
            _keep_rows = None
        else:
            _cm = np.asarray(cell_mask)
            _keep_rows = np.flatnonzero(_cm) if _cm.dtype == bool else np.sort(_cm.astype(np.int64))
        (group_candidates, iqr_collector, _n_chunks, _n_seen,
         t_stream_total) = _scan_column_major(
            _cm_layer, _keep_rows, cell_group_indices, n_groups, group_sizes,
            batch_index, n_batches, grp_count_accum, batch_cell_number_threshold,
            chrom_ranges, _chunk_normalizer, track_nnz, mu,
            min_total_count, thresholds,
            n_keep, iqr_normalize, verbose, block_features=block_features,
            check_sums=lambda sq, start: check_finite_sums(sq, what=_what, offset=start, lookup=_names_of),
            candidates_all=candidates_all)
        if verbose:
            print(f"  Streaming: {_n_chunks} column chunks, {_n_seen:,} cells, "
                  f"{t_stream_total:.2f}s")
        chrom_batches = []

    # Process each batch of chromosomes
    for batch_idx, batch_chroms in enumerate(chrom_batches):
        if verbose and len(chrom_batches) > 1:
            print(f"  Pass {batch_idx + 1}/{len(chrom_batches)}: "
                  f"{len(batch_chroms)} chromosomes")

        # Initialize accumulators for this batch only
        chrom_accum = {}
        for chrom in batch_chroms:
            cs, ce = chrom_ranges[chrom]
            n_feat = ce - cs
            chrom_accum[chrom] = {
                'dot': np.zeros((n_batches, n_feat, n_groups), dtype=np.float32),
                'sq': np.zeros((n_batches, n_feat), dtype=np.float32),
            }
            if track_nnz:
                chrom_accum[chrom]['nnz'] = np.zeros(
                    (n_feat, n_groups), dtype=np.int32
                )

        # Stream chunks
        n_chunks = 0
        n_cells_seen = 0
        t_stream = 0.0

        for chunk_csr, row_indices in ds.iter_chunks(
            modality, cell_mask=cell_mask, batch_size=batch_size,
            layer=_layer_to_read,
        ):
            if _chunk_normalizer is not None:
                chunk_csr = _chunk_normalizer(chunk_csr, row_indices)
            t_comp = time.time()
            actual_chunk = chunk_csr.shape[0]

            if cell_mask is not None:
                chunk_labels_idx = np.arange(
                    n_cells_seen, n_cells_seen + actual_chunk, dtype=np.int32
                )
            else:
                chunk_labels_idx = row_indices

            chunk_groups = cell_group_indices[chunk_labels_idx]
            chunk_batches = (None if batch_index is None
                             else batch_index[chunk_labels_idx])

            # Build indicator matrix
            Lam_chunk = sp_sparse.csr_matrix(
                (np.ones(actual_chunk, dtype=np.float32),
                 (np.arange(actual_chunk, dtype=np.int32), chunk_groups)),
                shape=(actual_chunk, n_groups)
            )

            if chunk_csr.dtype != np.float32 and chunk_csr.dtype != np.float64:
                chunk_csr = chunk_csr.astype(np.float32)

            # Convert to CSC for efficient column slicing
            chunk_csc = chunk_csr.tocsc()
            del chunk_csr

            # Process each chromosome in this batch
            for chrom in batch_chroms:
                cs, ce = chrom_ranges[chrom]
                chrom_slice = chunk_csc[:, cs:ce]

                if chrom_slice.nnz == 0:
                    continue

                # Dot products: (n_chrom_feat × n_groups), per batch
                if batch_index is None:
                    D_chrom = chrom_slice.T @ Lam_chunk
                    if sp_sparse.issparse(D_chrom):
                        D_chrom = D_chrom.toarray()
                    chrom_accum[chrom]['dot'][0] += D_chrom
                    del D_chrom
                    Xsq = chrom_slice.copy()
                    Xsq.data = Xsq.data ** 2
                    chrom_accum[chrom]['sq'][0] += np.asarray(Xsq.sum(axis=0)).ravel()
                    del Xsq
                else:
                    for b in np.unique(chunk_batches):
                        sel = np.flatnonzero(chunk_batches == b)
                        sub = chrom_slice[sel]
                        Db = sub.T @ Lam_chunk[sel]
                        if sp_sparse.issparse(Db):
                            Db = Db.toarray()
                        chrom_accum[chrom]['dot'][b] += Db
                        Xsq = sub.copy()
                        Xsq.data = Xsq.data ** 2
                        chrom_accum[chrom]['sq'][b] += np.asarray(Xsq.sum(axis=0)).ravel()
                        del Db, Xsq, sub

                # Nonzero counts
                if track_nnz:
                    X_bin = chrom_slice.copy()
                    X_bin.data = (X_bin.data != 0).astype(X_bin.dtype)    # a stored zero is not expressing
                    nnz_chrom = X_bin.T @ Lam_chunk
                    if sp_sparse.issparse(nnz_chrom):
                        nnz_chrom = nnz_chrom.toarray()
                    chrom_accum[chrom]['nnz'] += nnz_chrom.astype(np.int32)
                    del X_bin, nnz_chrom

            del chunk_csc, Lam_chunk
            t_stream += time.time() - t_comp
            n_chunks += 1
            n_cells_seen += actual_chunk

            if verbose and n_chunks % 500 == 0:
                print(f"    {n_chunks} chunks, {n_cells_seen:,} cells")
        t_stream_total += t_stream

        if verbose:
            print(f"  Streaming: {n_chunks} chunks, {n_cells_seen:,} cells, "
                  f"{t_stream:.2f}s")

        # Post-process this batch: per chromosome, freed progressively
        for chrom in batch_chroms:
            acc = chrom_accum.pop(chrom)  # frees RAM
            cs, ce = chrom_ranges[chrom]
            n_feat = ce - cs
            check_finite_sums(acc['sq'], what=_what, offset=cs, lookup=_names_of)

            # The cosine (per batch when batch_key is set, each with its own
            # feature and group norms, averaged over the batches where the
            # group met the cell threshold), COSG's penalty and the filters:
            # the same construction as the non-feature-batched path.
            genexlambda = _score_from_sums(
                acc['dot'][0] if batch_index is None else acc['dot'],
                acc['sq'][0] if batch_index is None else acc['sq'],
                group_sizes, acc.get('nnz'), mu=mu, thresholds=thresholds,
                min_total_count=min_total_count,
                batch_group_sizes=None if batch_index is None else grp_count_accum,
                batch_cell_number_threshold=batch_cell_number_threshold)
            del acc['dot'], acc['sq']
            acc.pop('nnz', None)

            # Extract per-group top-K candidates for this chromosome
            for k in range(n_groups):
                col = genexlambda[:, k]

                # Collect positive scores for global IQR
                n_positive = int((col > 0).sum())
                if iqr_collector is not None and n_positive:
                    iqr_collector[k].append(col[col > 0])

                # this chromosome's best features (positive ones only unless
                # every feature is a candidate), ties by index
                n_cand = min(n_keep, col.shape[0] if candidates_all else n_positive)
                if n_cand == 0:
                    continue
                loc = top_n(col, n_cand)
                cand_idx[k].append(loc.astype(np.int64) + cs)
                cand_val[k].append(col[loc].astype(np.float64))
                if sum(len(a) for a in cand_idx[k]) > 4 * n_keep:
                    _trim(k)

            del genexlambda

        del chrom_accum

    # per group: (feature indices, scores), best first, ties by index
    if group_candidates is None:
        group_candidates = {}
        for k in range(n_groups):
            if cand_idx[k]:
                _trim(k)
                group_candidates[k] = (cand_idx[k][0], cand_val[k][0])
            else:
                group_candidates[k] = (np.empty(0, dtype=np.int64), np.empty(0))

    # Global IQR + merge
    t_post = time.time()

    def _compute_group_iqr_from_collector(k):
        """Compute the all-values IQR for group k.

        We only stored positive scores in iqr_collector[k]; the full
        feature column is virtually padded with (n_genes - n_positive)
        zeros to compute the AnnData-equivalent IQR over all features.

        Materialises a (n_genes,) float32 array per group and frees it
        between groups — bounded transient memory.
        """
        all_positive = (
            np.concatenate(iqr_collector[k])
            if iqr_collector[k] else np.array([], dtype=np.float32)
        )
        n_zeros = n_genes - len(all_positive)
        if n_zeros < 0:
            # Should not happen — defensive
            n_zeros = 0
        full_col = np.concatenate(
            [np.zeros(n_zeros, dtype=np.float32),
             all_positive.astype(np.float32)]
        )
        return _group_scale(full_col, q_upper, q_lower, iqr_fallback, scale_rule, _collector_pooled_median[0])

    # the median rule's stand-in for a group with few positive scores (QC doc 52), from every
    # group's collected positive scores
    _collector_pooled_median = [
        _pooled_positive_median((np.concatenate(iqr_collector[j]) if iqr_collector[j]
                                 else np.array([], dtype=np.float32)) for j in range(n_groups))
        if iqr_normalize and scale_rule == "positive_median" else None]

    if use_sparse_output:
        scores_dict = {}
        dense_g, dense_f, dense_v = [], [], []
        collided = set()
        for k in range(n_groups):
            top_indices, top_scores_raw = group_candidates[k]
            if not len(top_indices):
                continue

            # All-values IQR (matches cosg.iqrLogNormalize on AnnData side)
            iqr = (
                _compute_group_iqr_from_collector(k)
                if iqr_normalize else None
            )

            # Lazy feature name lookup
            top_names = _lookup_feature_names(
                ds, feature_table, idx_col, name_col, top_indices
            )

            group_name = str(groups_order[k])
            for fi, name, score_raw in zip(top_indices, top_names, top_scores_raw):
                if iqr_normalize:
                    score = float(np.log1p(score_raw / iqr))
                else:
                    score = float(score_raw)
                if score > 0:
                    dense_g.append(k); dense_f.append(int(fi)); dense_v.append(score)
                    # a name held twice keeps its best score (the first)
                    if (group_name, str(name)) in scores_dict:
                        collided.add(str(name))
                    else:
                        scores_dict[(group_name, str(name))] = score

        del iqr_collector, group_candidates

        groups_order_str = [str(g) for g in groups_order]
        group_size_dict = {
            groups_order_str[k]: int(group_sizes[k])
            for k in range(n_groups)
        }

        if verbose:
            n_groups_with_entries = len({g for (g, _) in scores_dict.keys()})
            avg_per_group = (
                len(scores_dict) / n_groups_with_entries
                if n_groups_with_entries else 0
            )
            print(
                f"[cosg_{mode}] n_genes_user={n_genes_user} per group "
                f"(output_format={output_format!r}, "
                f"iqr_normalize={iqr_normalize}): "
                f"{len(scores_dict):,} entries written "
                f"(avg {avg_per_group:,.0f}/group; "
                f"entries with raw COSG score ≤ 0 dropped — "
                f"to include more, lower expressed_pct or mu)"
            )

        # Capture Post timing BEFORE _format_top_k_output runs so the
        # printed "Post" still represents the per-group IQR + scores_dict
        # construction phase (not the dense formatting phase).
        t_post_end = time.time()
        _warn_repeated_names(collided, output_format)
        result = _format_top_k_output(
            scores_dict, groups_order_str, group_size_dict,
            output_format=output_format,
            ds=ds, modality=modality, verbose=verbose,
            feature_name_col=feature_name_col,
            dense_entries=(dense_g, dense_f, dense_v),
        )
        if _floor_per_group:
            result['expressed_pct'] = floor_record(groups_order_str, floors)
        t_format = time.time() - t_post_end

        if verbose:
            print(
                f"[cosg_{mode}] Stream: {t_stream_total:.2f}s, "
                f"Post: {t_post_end - t_post:.2f}s, "
                f"Format: {t_format:.2f}s, "
                f"Total: {time.time() - t0:.2f}s"
            )
        if _opened_here:
            ds.close()
        return result

    else:
        # ndarray output: extract top-N per group, optionally IQR-normalize.
        names_out = np.empty((n_keep, n_groups), dtype=object)
        scores_out = np.empty((n_keep, n_groups), dtype=np.float32)

        for k in range(n_groups):
            top_indices, top_raw = group_candidates[k]

            if not len(top_indices):
                names_out[:, k] = ""
                scores_out[:, k] = 0.0
                continue

            top_scores_vals = np.asarray(top_raw, dtype=np.float32)
            if iqr_normalize:
                # Same all-values IQR formula as sparse path
                iqr = _compute_group_iqr_from_collector(k)
                top_scores_vals = np.log1p(
                    np.maximum(top_scores_vals, 0.0) / iqr
                ).astype(np.float32)

            top_names = _lookup_feature_names(
                ds, feature_table, idx_col, name_col, top_indices
            )

            n_fill = len(top_indices)
            names_out[:n_fill, k] = top_names
            scores_out[:n_fill, k] = top_scores_vals
            if n_fill < n_keep:
                names_out[n_fill:, k] = ""
                scores_out[n_fill:, k] = 0.0

        if verbose:
            print(f"[cosg_{mode}] Post: {time.time() - t_post:.2f}s, "
                  f"Total: {time.time() - t0:.2f}s")

        if _opened_here:
            ds.close()
        result = {
            'names': names_out,
            'scores': scores_out,
            'groups_order': groups_order,
        }
        if _floor_per_group:
            result['expressed_pct'] = floor_record(groups_order, floors)
        return result


class _GPUUnsupported(Exception):
    """A request the GPU path cannot serve; run_cosg_cytome falls back to the
    CPU for device='auto' and raises NotImplementedError for device='gpu'."""


def _run_cosg_cytome_gpu_impl(
    ds,                                 # open cytome.Dataset
    groupby,
    *,
    mu=1.0,
    remove_lowly_expressed=True,
    expressed_pct=0.1,
    expressed_min_num_cells_in_target_group=3,
    cell_mask=None,
    batch_size=4096,
    verbose=True,
    modality="RNA",
    layer="auto",
    compute_on_fly=True,
    use_cached_stats=True,
    feature_name_col=None,
    batch_key=None,
    batch_cell_number_threshold=3,
):
    """GPU accumulation and scoring on an open ``cytome.Dataset``.

    Returns ``(genexlambda, gene_names, groups_order, group_sizes, floors,
    floor_per_group)``: every feature's scores on the host, which
    run_cosg_cytome turns into its result with the CPU path's own
    post-processing (one selector, one scale, every output format). The caller
    owns the Dataset. Raises :class:`_GPUUnsupported` for a column-major
    matrix, whose row chunks this path reads do not exist.
    """
    import cupy as cp
    import cupyx.scipy.sparse as cp_sparse
    from ._gpu import _apply_cosg_penalty_gpu

    t0 = time.time()

    labels, gene_names, all_labels, categories = _get_labels_and_gene_names(
        ds, groupby, cell_mask, modality=modality, feature_name_col=feature_name_col
    )
    groups_order, cell_group_indices, group_sizes = _build_group_mapping(
        labels, categories
    )
    floors, _floor_per_group = _resolve_floors(expressed_pct, groups_order, group_sizes,
                                               all_labels, categories)

    n_genes = len(gene_names)
    n_groups = len(groups_order)

    # Resolve layer + normalizer once before streaming (same logic as CPU path)
    _layer_to_read, _need_normalizer = _resolve_layer_to_read(
        ds, modality, layer, compute_on_fly,
    )
    if _column_major_layer(ds, modality, _layer_to_read) is not None:
        raise _GPUUnsupported(
            f"'{modality}_{_layer_to_read}' is stored column-major, and the GPU path reads row chunks")
    _chunk_normalizer = (
        _resolve_chunk_normalizer(layer, ds, modality, use_cached_stats)
        if _need_normalizer else None
    )

    group_sizes_gpu = cp.array(group_sizes, dtype=cp.float32)
    group_norms_gpu = cp.sqrt(group_sizes_gpu)

    if verbose:
        print(f"[cytome_gpu] {len(labels):,} cells × {n_genes:,} genes, "
              f"{n_groups} groups")
        _emit_filter_log(
            "cytome_gpu",
            remove_lowly_expressed=remove_lowly_expressed,
            expressed_pct=describe_floor(floors, _floor_per_group, expressed_pct),
            expressed_min_num_cells_in_target_group=
                expressed_min_num_cells_in_target_group,
            mu=mu,
        )

    # Batch dimension, mirroring the CPU path exactly: per-batch cosine with
    # that batch's own norms, averaged over batches meeting the cell threshold.
    if batch_key is None:
        n_batches, batch_index = 1, None
    else:
        _blab = _read_batch_labels(ds, batch_key, cell_mask, len(labels))
        _bo, _bi = np.unique(_blab, return_inverse=True)
        n_batches, batch_index = len(_bo), _bi
        need_gb = n_batches * n_genes * n_groups * 4 / 1e9
        if verbose:
            print(f"[cytome_gpu] batch_key: {n_batches} batches x {n_genes:,} features "
                  f"x {n_groups} groups -> {need_gb:.2f} GB of DEVICE accumulators")
        # GPU RAM is far scarcer than host RAM, so the ceiling is lower here.
        if need_gb > 8.0:
            raise MemoryError(
                f"batch-aware COSG needs {need_gb:.1f} GB of device accumulators "
                f"({n_batches} batches x {n_genes:,} features x {n_groups} groups x 4 B), "
                f"more than this path is willing to pin on the GPU. Use device='cpu' "
                f"(16 GB host ceiling), reduce the feature space with min_total_count, "
                f"or coarsen `groupby`.")

    # Initialize GPU accumulators. The group sizes per batch are counted on
    # the host: a device add per cell was one kernel launch each.
    dot_accum = cp.zeros((n_batches, n_genes, n_groups), dtype=cp.float32)
    gene_sq_accum = cp.zeros((n_batches, n_genes), dtype=cp.float32)
    grp_count_host = np.zeros((n_batches, n_groups), dtype=np.float64)

    if tracks_floor(remove_lowly_expressed, floors, expressed_min_num_cells_in_target_group):
        nnz_accum = cp.zeros((n_genes, n_groups), dtype=cp.float32)
    else:
        nnz_accum = None

    n_chunks = 0
    n_cells_seen = 0
    t_transfer = 0.0
    t_compute = 0.0

    for chunk_csr, row_indices in ds.iter_chunks(
        modality, cell_mask=cell_mask, batch_size=batch_size,
        layer=_layer_to_read,
    ):
        if _chunk_normalizer is not None:
            chunk_csr = _chunk_normalizer(chunk_csr, row_indices)
        actual_chunk = chunk_csr.shape[0]

        # Map to label positions
        if cell_mask is not None:
            chunk_labels_idx = np.arange(
                n_cells_seen, n_cells_seen + actual_chunk, dtype=np.int32
            )
        else:
            chunk_labels_idx = row_indices

        chunk_groups = cell_group_indices[chunk_labels_idx]

        # Transfer chunk to GPU
        t_xfer = time.time()
        if chunk_csr.dtype != np.float32:
            chunk_csr = chunk_csr.astype(np.float32)
        X_gpu = cp_sparse.csr_matrix(
            (cp.array(chunk_csr.data),
             cp.array(chunk_csr.indices),
             cp.array(chunk_csr.indptr)),
            shape=chunk_csr.shape
        )

        # Build dense indicator on GPU: (n_groups × chunk)
        Lam_dense_T = cp.zeros((n_groups, actual_chunk), dtype=cp.float32)
        Lam_dense_T[cp.array(chunk_groups), cp.arange(actual_chunk)] = 1.0

        cp.cuda.Stream.null.synchronize()
        t_transfer += time.time() - t_xfer

        # Compute on GPU
        t_comp = time.time()

        # Dot products: (groups × chunk) @ (chunk × genes) → (groups × genes)
        if batch_index is None:
            D_chunk = Lam_dense_T @ X_gpu
            if cp_sparse.issparse(D_chunk):
                D_chunk = D_chunk.toarray()
            dot_accum[0] += D_chunk.T  # → (genes × groups)
            Xsq = X_gpu.copy()
            Xsq.data = Xsq.data ** 2
            gene_sq_accum[0] += cp.array(Xsq.sum(axis=0)).ravel()
            del Xsq, D_chunk
        else:
            chunk_batches = batch_index[chunk_labels_idx]
            np.add.at(grp_count_host, (chunk_batches, chunk_groups), 1.0)
            for b in np.unique(chunk_batches):
                sel = np.flatnonzero(chunk_batches == b)
                sel_gpu = cp.array(sel)
                Xb = X_gpu[sel_gpu]
                Lb = Lam_dense_T[:, sel_gpu]
                Db = Lb @ Xb
                if cp_sparse.issparse(Db):
                    Db = Db.toarray()
                dot_accum[int(b)] += Db.T
                Xsq = Xb.copy()
                Xsq.data = Xsq.data ** 2
                gene_sq_accum[int(b)] += cp.array(Xsq.sum(axis=0)).ravel()
                del Xb, Lb, Db, Xsq

        # Nonzero counts
        if nnz_accum is not None:
            X_bin = X_gpu.copy()
            X_bin.data = (X_bin.data != 0).astype(cp.float32)    # a stored zero is not expressing
            nnz_chunk = Lam_dense_T @ X_bin
            if cp_sparse.issparse(nnz_chunk):
                nnz_chunk = nnz_chunk.toarray()
            nnz_accum += nnz_chunk.T
            del X_bin, nnz_chunk

        cp.cuda.Stream.null.synchronize()
        t_compute += time.time() - t_comp

        del X_gpu, Lam_dense_T
        cp.get_default_memory_pool().free_all_blocks()

        n_chunks += 1
        n_cells_seen += actual_chunk

        if verbose and n_chunks % 500 == 0:
            print(f"  Processed {n_chunks} chunks, {n_cells_seen:,} cells")

    if verbose:
        print(f"[cytome_gpu] Transfer: {t_transfer:.3f}s, "
              f"Compute: {t_compute:.3f}s")

    # A NaN or an infinity in a feature makes its sum of squares non-finite
    from ._validate import check_finite_sums
    check_finite_sums(cp.asnumpy(gene_sq_accum), gene_names,
                      what=_matrix_label(modality, _layer_to_read, layer, _chunk_normalizer))

    group_norms_safe = cp.maximum(group_norms_gpu, 1e-12)

    if batch_index is None:
        gene_norms = cp.maximum(cp.sqrt(gene_sq_accum[0]), 1e-12)
        cosine_sim = dot_accum[0] / (gene_norms[:, None] * group_norms_safe[None, :])
    else:
        # Per batch, with that batch's own gene and group norms, averaged over the
        # batches where the group met the cell threshold — identical construction
        # to the CPU path, just on device.
        grp_count_accum = cp.asarray(grp_count_host)
        _acc = cp.zeros((n_genes, n_groups), dtype=cp.float64)
        _valid = cp.zeros(n_groups, dtype=cp.float64)
        for b in range(n_batches):
            ok = (grp_count_accum[b] >= batch_cell_number_threshold).astype(cp.float64)
            if not bool(ok.any()):
                continue
            gn = cp.maximum(cp.sqrt(gene_sq_accum[b]), 1e-12)
            gp = cp.maximum(cp.sqrt(grp_count_accum[b]), 1e-12)
            _acc += (dot_accum[b] / (gn[:, None] * gp[None, :])) * ok[None, :]
            _valid += ok
        _acc /= cp.where(_valid > 0, _valid, cp.float64(1.0))[None, :]
        cosine_sim = cp.where(cp.isfinite(_acc), _acc, 0.0).astype(cp.float32)
        del _acc, _valid
    del dot_accum, gene_sq_accum

    genexlambda_gpu = _apply_cosg_penalty_gpu(cosine_sim, mu)
    del cosine_sim

    # Expressed pct filter on GPU: each group's own threshold, max(N_k *
    # pct_k, min_cells), in float32 as the device forms it (a float32 size
    # times a Python float); none for a group without cells (a floor of 0
    # keeps the minimum of cells)
    if nnz_accum is not None:
        thresholds = floor_thresholds(group_sizes, floors,
                                      expressed_min_num_cells_in_target_group,
                                      dtype=np.float32, skip_empty=True)
        for k in range(n_groups):
            if thresholds[k] > 0:
                genexlambda_gpu[nnz_accum[:, k] < float(thresholds[k]), k] = -1.0
        del nnz_accum

    genexlambda = cp.asnumpy(genexlambda_gpu)
    del genexlambda_gpu
    cp.get_default_memory_pool().free_all_blocks()
    if verbose:
        print(f"[cytome_gpu] Scored: {time.time() - t0:.3f}s")
    return genexlambda, gene_names, groups_order, group_sizes, floors, _floor_per_group


def run_cosg_cytome_gpu(
    cytome_path,
    groupby,
    n_genes_user=50,
    mu=1.0,
    remove_lowly_expressed=True,
    expressed_pct=0.1,
    expressed_min_num_cells_in_target_group=3,
    cell_mask=None,
    batch_size=4096,
    verbose=True,
    modality="RNA",
    layer="auto",
    compute_on_fly=True,
    use_cached_stats=True,
    batch_key=None,
    batch_cell_number_threshold=3,
    **_deprecated_kwargs,
):
    """GPU version of cytome-native COSG — convenience shortcut for
    ``run_cosg_cytome(..., device='gpu')``.

    Same kwargs as ``run_cosg_cytome`` for the ones listed here, and its
    ``'ndarray'`` result (``run_cosg_cytome(..., device='gpu',
    output_format=...)`` gives the other formats, computed on the GPU too).
    The GPU path computes no p-values, does not batch features, and does not
    read column-major matrices; ``device='auto'`` runs those requests on the
    CPU.

    Parameters
    ----------
    expressed_pct : float, dict, pandas.Series or sequence of float, default 0.1
        The detection floor, one fraction for every group or one per group,
        as in :func:`run_cosg_cytome` (a dict or Series keyed by group, or a
        sequence in the order of ``groups_order``), e.g.
        ``expressed_pct={'Microglia': 0.05, 'Astro': 0.1, 'Oligo': 0.1}``.
        Per-group floors are returned under ``'expressed_pct'``.
    modality : str
        Modality to use (``"RNA"``, ``"ATAC"``, ``"GA"``, ``"tiles"``).
        Routes via the cytome modality registry to the matching feature
        entity (RNA → ``genes``, ATAC → ``peaks``, GA → ``GA_genes``,
        tiles → ``tiles``) and matrix layer.
    layer : str, default ``"auto"``
        Measurement layer to read, resolved as in :func:`run_cosg_cytome`:

        - ``"auto"``: ``log1p`` for RNA/GA, ``tfidf`` for ATAC/tiles (and
          the stored values as given when they are not integer counts).
        - ``"log1p"``: ``log1p(counts / depth * 1e4)`` on the fly (CP10K, as
          scanpy's ``normalize_total(target_sum=1e4)`` + ``log1p``), or a
          stored ``{modality}_log1p`` layer.
        - ``"counts"``: raw counts, used as they are.
        - ``"infog"`` / ``"tfidf"``: PIASO's normalizations, on the fly from
          their cached parameters or from a stored layer.
        - Any other materialized layer name is read as-is.

    Migration notes — passing the old ``top_k`` or ``n_top_genes``
    kwargs raises ``TypeError`` with a migration hint:

        - ``n_top_genes=N`` → ``n_genes_user=N``
        - ``top_k=N``       → ``run_cosg_cytome(..., n_genes_user=N,
          output_format='dict')``.
    """
    # Validate API BEFORE importing cupy so the K-HARD migration TypeError
    # fires even on machines without GPU/CuPy installed.
    _validate_cosg_api_v2(
        _deprecated_kwargs,
        output_format="ndarray",
        iqr_normalize=False,
        fn_name="run_cosg_cytome_gpu",
    )
    return run_cosg_cytome(
        cytome_path,
        groupby=groupby,
        n_genes_user=n_genes_user,
        mu=mu,
        remove_lowly_expressed=remove_lowly_expressed,
        expressed_pct=expressed_pct,
        expressed_min_num_cells_in_target_group=
            expressed_min_num_cells_in_target_group,
        cell_mask=cell_mask,
        batch_size=batch_size,
        verbose=verbose,
        modality=modality,
        layer=layer,
        compute_on_fly=compute_on_fly,
        use_cached_stats=use_cached_stats,
        batch_key=batch_key,
        batch_cell_number_threshold=batch_cell_number_threshold,
        device="gpu",
    )
