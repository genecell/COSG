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
    cytome_Dataset = getattr(_cytome(), "Dataset", None)
    if cytome_Dataset is not None and cytome_Dataset not in candidates:
        candidates.append(cytome_Dataset)
    result = tuple(candidates)
    _cytome_dataset_classes._cache = result
    return result


def _is_cytome_dataset(obj):
    """True if ``obj`` is an open cytome Dataset (any known subclass)."""
    classes = _cytome_dataset_classes()
    return bool(classes) and isinstance(obj, classes)


def _open_or_use(cytome_path):
    """Return ``(ds, opened_here)`` for a cytome path-or-Dataset arg.

    - ``cytome_path`` is an open cytome Dataset → return as-is with
      ``opened_here=False`` (caller owns lifecycle; do NOT close).
    - ``cytome_path`` is ``str`` / ``Path`` → ``cytome.open(...)`` and
      return ``opened_here=True`` so the caller closes at exit.

    Mirrors the pattern used by ``piaso.tools._runGDR._open_cytome`` so
    the cytome streaming COSG functions accept already-open Datasets
    without forcing a close/reopen cycle.
    """
    if _is_cytome_dataset(cytome_path):
        return cytome_path, False
    return _cytome().open(cytome_path), True


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
                          output_format, ds, modality, verbose, feature_name_col=None):
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

        # Allocate float32 to keep memory in check
        scores_df = pd.DataFrame(
            np.zeros((n_features, n_groups), dtype=np.float32),
            index=pd.Index(all_features, name="feature"),
            columns=groups_order_str,
        )
        # Pivot via direct .values mutation for speed (loc-based assignment is slow on large indices)
        feature_to_row = {f: i for i, f in enumerate(all_features)}
        group_to_col = {g: i for i, g in enumerate(groups_order_str)}
        values = scores_df.values
        for (group, feature), score in scores_dict.items():
            row = feature_to_row.get(feature)
            if row is None:
                continue  # feature not in current feature table — defensive
            col = group_to_col[group]
            values[row, col] = score
        return {
            "scores_df": scores_df,
            "groups_order": groups_order_str,
            "group_sizes": group_size_dict,
        }

    # Should never reach here — output_format is validated upstream.
    raise ValueError(f"Invalid output_format='{output_format}'.")


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


def _quantile_with_zero_pad(positive, q, n_total):
    """Compute np.quantile of `positive` virtually padded to length n_total
    with leading zeros. Used by the chromosome-batched IQR computation
    where only positive scores are kept across chromosomes but the IQR
    must be computed over the full feature range (zeros + positives).

    Equivalent to:
        np.quantile(np.concatenate([np.zeros(n_total - len(positive)),
                                    np.sort(positive)]), q)
    but without the (n_total × 4 bytes) allocation. Used as Option α's
    sibling — Option α materializes the array, this is the zero-alloc
    variant kept for future use.
    """
    positive = np.asarray(positive)
    n_pos = len(positive)
    n_zeros = n_total - n_pos
    if n_total <= 0:
        return 0.0
    if n_pos == 0:
        return 0.0
    sorted_pos = np.sort(positive)
    r = q * (n_total - 1)
    floor = int(np.floor(r))
    frac = r - floor

    def at(i):
        return 0.0 if i < n_zeros else float(sorted_pos[i - n_zeros])

    if floor + 1 >= n_total:
        return at(n_total - 1)
    return at(floor) * (1 - frac) + at(floor + 1) * frac


def _get_labels_and_gene_names(ds, groupby, cell_mask=None, modality="RNA", load_gene_names=True, feature_name_col=None):
    """Extract group labels, feature names, and categorical ordering from a cytome dataset.

    When load_gene_names=False, gene_names is returned as None (lazy loading).

    Modality dispatch (RNA/GA/ATAC/tiles) is delegated to
    ``cytome.modality_feature_table_info`` — single source
    of truth shared with cytome.io and piaso plotting.
    """
    from cytome import modality_feature_table_info

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


def _resolve_auto_layer(layer, modality):
    """Resolve ``layer='auto'`` to a modality-appropriate normalization.

    Idempotent: any explicit layer (``'counts'|'log1p'|'infog'|'tfidf'`` or a
    materialised layer name) is returned unchanged. Called at the top of both
    low-level interpreters so every read path (cpu/gpu, on-the-fly/materialised)
    resolves ``auto`` identically regardless of entry point.
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
    layer = _resolve_auto_layer(layer, modality)
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
        return lambda chunk, indices: _normalize_chunk_infog(
            chunk, cd[indices], ig, scale, counts_sum, thr,
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
        return lambda chunk, indices: _normalize_chunk_tfidf(
            chunk, cd[indices], idf, scale_factor,
        )

    if layer == "log1p":
        from cytome import modality_cell_depth

        from ._normalize import normalize_chunk_log1p as _normalize_chunk_log1p
        cell_depth = modality_cell_depth(
            ds, modality, use_cached_stats=use_cached_stats,
        )
        # Scale: prefer one cached on the file (piaso.pp.normalize_log1p
        # records it at write-back time, so a COSG run matches whatever PIASO
        # already wrote); else the conventional 1e4.
        scale = 1e4
        params = ds.metadata.get(f"{modality}_log1p_params")
        if params is not None:
            scale = float(params.get("scale_factor", 1e4))
        return lambda chunk, indices: _normalize_chunk_log1p(
            chunk, cell_depth[indices], scale,
        )

    # Unreachable: the layer was validated against _ON_THE_FLY_LAYERS above,
    # and every member has a branch. Kept as a guard so adding a name to that
    # tuple without a branch fails loudly rather than returning None (which
    # the caller would read as "raw counts, no transform").
    raise AssertionError(
        f"layer='{layer}' is in _ON_THE_FLY_LAYERS but has no normalizer branch"
    )


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
    layer = _resolve_auto_layer(layer, modality)
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


def _get_chrom_col_ranges(ds, feature_table, idx_col):
    """Get chromosome → (start_col, end_col) column ranges from entity table.

    Tiles/peaks are registered in chromosome order during quantification,
    so each chromosome's features form a contiguous column range.
    """
    rows = ds._conn.execute(
        f'SELECT chr, MIN({idx_col}) as cs, MAX({idx_col}) + 1 as ce '
        f'FROM {feature_table} GROUP BY chr ORDER BY MIN({idx_col})'
    ).fetchall()
    chrom_ranges = {}
    chrom_order = []
    for chrom, cs, ce in rows:
        chrom_ranges[str(chrom)] = (int(cs), int(ce))
        chrom_order.append(str(chrom))
    return chrom_ranges, chrom_order


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
         'records': [{'group', 'feature', 'score', 'rank'?}, ...]}
    """
    ds, _opened = _open_or_use(cytome_path)
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

        ds.metadata[key_added] = {
            "groupby": groupby,
            "modality": modality,
            "n_genes_user": int(n_genes_user),
            "records": records,
        }
        try:
            ds.flush()
        except Exception:
            pass
    finally:
        if _opened:
            ds.close()


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
    feature_batching="auto",
    min_total_count=0,
    layer="auto",
    compute_on_fly=True,
    use_cached_stats=True,
    device="cpu",
    write_to_cytome=False,
    key_added="cosg",
    feature_name_col=None,
    batch_key=None,
    batch_cell_number_threshold=3,
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
    expressed_pct : float
        Minimum fraction of cells expressing a gene.
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
          cached ``{modality}_log1p_params['scale_factor']`` overrides the 1e4
          default. The ``auto`` choice for RNA/GA. Like all non-``counts`` layers
          it needs **piaso** — imported *lazily at call time* only, so ``import
          cosg`` itself stays piaso-free; a missing piaso raises a clear error
          pointing at ``pip install piaso`` or ``layer='counts'``.
        - ``"counts"``: **raw counts, used as-is with no transform**
          (``_resolve_chunk_normalizer`` returns ``None``). The one path with
          **zero piaso dependency**. COSG's cosine specificity is per-feature
          scale-robust so this still works, but it does not correct cell-depth
          bias — use it when piaso is unavailable, or for pre-normalized inputs.
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
    feature_batching : str
        Feature batching strategy for ATAC data with millions of features.
        ``"auto"`` selects ``"chromosome"`` when n_features > 100K, else ``"none"``.
        ``"none"`` uses the original single-pass approach (RNA default).
        ``"chromosome"`` batches features by chromosome in a single pass (~950 MB).
        ``"disk"`` does multi-pass processing, one chromosome batch per pass (~200 MB).
    min_total_count : int
        Pre-filter features with total count across all groups <= this value.
        Only applied when feature_batching != "none".

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
    from ._cpu_streaming import _apply_cosg_penalty_cpu

    t0 = time.time()

    use_sparse_output, iqr_normalize = _validate_cosg_api_v2(
        _deprecated_kwargs, output_format, iqr_normalize,
        fn_name="run_cosg_cytome",
    )

    # Round 8 (2026-05-23): accept cytome.Dataset OR path. We only close
    # when we opened it ourselves, so callers who pass an open Dataset
    # retain ownership of its lifecycle.
    ds, _opened_here = _open_or_use(cytome_path)

    # Round 8: device='gpu' / 'auto' dispatch — re-route to the GPU
    # implementation BEFORE consuming any cytome chunks on CPU. Validated
    # against the GPU-supported subset of kwargs (no output_format != ndarray,
    # no feature_batching != none).
    from ._backend import get_device
    _resolved_device = get_device(device)
    if _resolved_device == "gpu":
        if use_sparse_output:
            raise NotImplementedError(
                "device='gpu' currently supports only output_format='ndarray'. "
                "For 'dict' / 'long' / 'dense', use device='cpu'."
            )
        if feature_batching not in (None, "none", "auto"):
            # 'auto' is fine — _run_cosg_cytome_gpu_impl falls back to 'none'
            # at its modality check (RNA always uses none anyway).
            raise NotImplementedError(
                f"device='gpu' does not support feature_batching="
                f"{feature_batching!r}. Use device='cpu' for "
                f"chromosome / disk batching, or feature_batching='none'."
            )
        try:
            _gpu_result = _run_cosg_cytome_gpu_impl(
                ds, groupby,
                n_genes_user=n_genes_user, mu=mu,
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
        finally:
            if _opened_here:
                ds.close()
        if write_to_cytome:
            _persist_cosg_markers(cytome_path, _gpu_result, groupby,
                                  modality, n_genes_user, key_added)
        return _gpu_result

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

    # Resolve layer + normalizer ONCE before the streaming loops so the
    # decision (read materialised layer vs read counts and transform on
    # the fly) is made up front.
    _layer_to_read, _need_normalizer = _resolve_layer_to_read(
        ds, modality, layer, compute_on_fly,
    )
    _chunk_normalizer = (
        _resolve_chunk_normalizer(layer, ds, modality, use_cached_stats)
        if _need_normalizer else None
    )

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
            output_format=output_format,
            feature_name_col=feature_name_col,
            batch_labels=_batch_labels_fb,
            batch_cell_number_threshold=batch_cell_number_threshold,
        )
        if write_to_cytome:
            _persist_cosg_markers(cytome_path, _batched_result, groupby,
                                  modality, n_genes_user, key_added)
        return _batched_result

    # Original "none" path (RNA default) — unchanged below
    labels, gene_names, all_labels, categories = _get_labels_and_gene_names(
        ds, groupby, cell_mask, modality=modality, feature_name_col=feature_name_col
    )

    # Batch labels, aligned to the same cells as `labels`.
    batch_labels = (None if batch_key is None
                    else _read_batch_labels(ds, batch_key, cell_mask, len(labels)))
    groups_order, cell_group_indices, group_sizes = _build_group_mapping(
        labels, categories
    )

    n_genes = len(gene_names)
    n_groups = len(groups_order)
    group_norms = np.sqrt(group_sizes)

    if verbose:
        print(f"[cytome_cpu] {len(labels):,} cells × {n_genes:,} genes, "
              f"{n_groups} groups")
        _emit_filter_log(
            "cytome_cpu",
            remove_lowly_expressed=remove_lowly_expressed,
            expressed_pct=expressed_pct,
            expressed_min_num_cells_in_target_group=
                expressed_min_num_cells_in_target_group,
            mu=mu,
        )

    # Guard: the feature (var) table and the modality's matrix must agree on the
    # feature count. A 0-row / out-of-sync var table — e.g. a modality whose
    # feature table was dropped by an older cytome.merge — against a populated
    # matrix would otherwise crash deep in the accumulation with a cryptic
    # "operands could not be broadcast" error. Fail early with a fix hint.
    _mat_name, _mat_cols = None, None
    for _mn, _nc in ds._conn.execute(
        "SELECT matrix_name, n_cols FROM matrix_meta"
    ).fetchall():
        if _mn.split("_", 1)[0] == modality:
            _mat_name, _mat_cols = _mn, int(_nc)
            if _mn == f"{modality}_counts":
                break  # prefer the counts matrix as the reference
    if _mat_cols is not None and _mat_cols != n_genes:
        from cytome import modality_feature_table_info
        _ftab, _idxc, _ = modality_feature_table_info(ds, modality, feature_name_col)
        raise ValueError(
            f"COSG: modality '{modality}' feature table '{_ftab}' has "
            f"{n_genes:,} rows, but matrix '{_mat_name}' has {_mat_cols:,} "
            f"columns — the feature/var table is missing or out of sync (e.g. an "
            f"older cytome.merge carried the matrix but not the table). Rebuild "
            f"it so the '{_idxc}' order matches the matrix columns, e.g. "
            f"ds.set_entity('{_ftab}', <source {_ftab} dataframe>); ds.flush(), "
            f"then re-run COSG."
        )

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
        need_gb = n_batches * n_genes * n_groups * 4 / 1e9
        if verbose:
            print(f"[cytome_cpu] batch_key: {n_batches} batches x {n_genes:,} features "
                  f"x {n_groups} groups -> {need_gb:.2f} GB of accumulators")
        if need_gb > 16.0:
            raise MemoryError(
                f"batch-aware COSG needs {need_gb:.1f} GB of accumulators "
                f"({n_batches} batches x {n_genes:,} features x {n_groups} groups x 4 B). "
                f"Reduce the feature space (min_total_count), coarsen `groupby`, or use "
                f"feature_batching='chrom' to process features in blocks.")

    dot_accum = np.zeros((n_batches, n_genes, n_groups), dtype=np.float32)
    gene_sq_accum = np.zeros((n_batches, n_genes), dtype=np.float32)
    grp_count_accum = np.zeros((n_batches, n_groups), dtype=np.float64)

    if remove_lowly_expressed and expressed_pct > 0:
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
            (np.ones(actual_chunk, dtype=np.float32),
             (np.arange(actual_chunk, dtype=np.int32), chunk_groups)),
            shape=(actual_chunk, n_groups)
        )

        # Ensure float dtype for sparse operations
        if chunk_csr.dtype != np.float32 and chunk_csr.dtype != np.float64:
            chunk_csr = chunk_csr.astype(np.float32)

        # Accumulate dot products: (genes × groups), per batch
        if batch_index is None:
            D_chunk = chunk_csr.T @ Lam_chunk
            if sp_sparse.issparse(D_chunk):
                D_chunk = D_chunk.toarray()
            dot_accum[0] += D_chunk

            Xsq = chunk_csr.copy()
            Xsq.data = Xsq.data ** 2
            gene_sq_accum[0] += np.asarray(Xsq.sum(axis=0)).ravel()
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
            X_bin = chunk_csr.copy()
            X_bin.data = np.ones_like(X_bin.data)
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

    # Compute cosine similarity from accumulated sums
    t_post = time.time()

    if batch_index is None:
        gene_norms = np.maximum(np.sqrt(gene_sq_accum[0]), 1e-12)
        group_norms_safe = np.maximum(group_norms, 1e-12)
        cosine_sim = dot_accum[0] / (gene_norms[:, None] * group_norms_safe[None, :])
    else:
        # Per batch: cosine with that batch's own norms, groups below the cell
        # threshold zeroed, then averaged over the batches where the group was
        # valid. A group valid in no batch is NaN -> 0 (it gets no evidence).
        acc = np.zeros((n_genes, n_groups), dtype=np.float64)
        valid_counts = np.zeros(n_groups, dtype=np.float64)
        for b in range(n_batches):
            valid = (grp_count_accum[b] >= batch_cell_number_threshold).astype(np.float64)
            if not valid.any():
                continue
            gn = np.maximum(np.sqrt(gene_sq_accum[b]), 1e-12)
            gp = np.maximum(np.sqrt(grp_count_accum[b]), 1e-12)
            sim_b = dot_accum[b] / (gn[:, None] * gp[None, :])
            acc += sim_b * valid[None, :]
            valid_counts += valid
        with np.errstate(invalid="ignore", divide="ignore"):
            acc /= valid_counts[None, :]
        n_dropped = int((valid_counts == 0).sum())
        if verbose and n_dropped:
            print(f"[cytome_cpu] batch_key: {n_dropped}/{n_groups} group(s) had "
                  f"< {batch_cell_number_threshold} cells in every batch; scored 0")
        cosine_sim = np.nan_to_num(acc, nan=0.0).astype(np.float32)

    del dot_accum, gene_sq_accum

    # Apply COSG penalty formula
    genexlambda = _apply_cosg_penalty_cpu(cosine_sim, mu)
    del cosine_sim

    # Apply expressed_pct filter
    if nnz_accum is not None:
        for k in range(n_groups):
            n_cells_k = group_sizes[k]
            if n_cells_k > 0:
                threshold = max(
                    n_cells_k * expressed_pct,
                    expressed_min_num_cells_in_target_group
                )
                genexlambda[nnz_accum[:, k] < threshold, k] = -1
        del nnz_accum

    # --- Sparse output (dict/long/dense): extract top-K per group ---
    if use_sparse_output:
        scores_dict = {}
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
                q_up = float(np.quantile(col, q_upper))
                q_lo = float(np.quantile(col, q_lower))
                iqr = q_up - q_lo
                if iqr <= 0:
                    iqr = q_up if q_up > 0 else 1e-6
            else:
                iqr = None  # raw-λ mode

            # Extract top-K via argpartition (O(n), no full sort)
            n_keep = min(n_genes_user, n_positive)
            top_idx = np.argpartition(col, -n_keep)[-n_keep:]

            # Normalize (or keep raw) only the top-K
            if iqr_normalize:
                top_scores = np.log1p(col[top_idx] / iqr).astype(np.float32)
            else:
                top_scores = col[top_idx].astype(np.float32)
            top_names = gene_names[top_idx]
            group_name = str(groups_order[k])

            for name, score in zip(top_names, top_scores):
                if score > 0:
                    scores_dict[(group_name, str(name))] = float(score)

        del genexlambda

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
                f"[cytome_cpu] n_genes_user={n_genes_user} per group "
                f"(output_format={output_format!r}, "
                f"iqr_normalize={iqr_normalize}): "
                f"{len(scores_dict):,} entries written "
                f"(avg {avg_per_group:,.0f}/group; "
                f"entries with raw COSG score ≤ 0 dropped — "
                f"to include more, lower expressed_pct or mu)"
            )
            print(f"[cytome_cpu] Post-processing: {time.time() - t_post:.3f}s")
            print(f"[cytome_cpu] TOTAL: {time.time() - t0:.3f}s")

        # Format dispatch: dict | long | dense
        _result = _format_top_k_output(
            scores_dict, groups_order_str, group_size_dict,
            output_format=output_format,
            ds=ds, modality=modality, verbose=verbose,
            feature_name_col=feature_name_col,
        )
        if _opened_here:
            ds.close()
        if write_to_cytome:
            _persist_cosg_markers(cytome_path, _result, groupby,
                                  modality, n_genes_user, key_added)
        return _result

    # --- ndarray output: (names, scores) per-group top-N ---
    names_out = np.empty((n_genes_user, n_groups), dtype=object)
    scores_out = np.empty((n_genes_user, n_groups), dtype=np.float32)

    # Pre-compute per-group IQR once if iqr_normalize is requested. Same
    # all-values IQR semantics as the sparse path, so the two output
    # shapes produce consistent normalized scores for the same input.
    if iqr_normalize:
        iqr_per_group = np.empty(n_groups, dtype=np.float32)
        for k in range(n_groups):
            col_for_iqr = genexlambda[:, k].copy()
            col_for_iqr[col_for_iqr < 0] = 0
            q_up = float(np.quantile(col_for_iqr, q_upper))
            q_lo = float(np.quantile(col_for_iqr, q_lower))
            iqr = q_up - q_lo
            if iqr <= 0:
                iqr = q_up if q_up > 0 else 1e-6
            iqr_per_group[k] = iqr

    for k in range(n_groups):
        scores_k = genexlambda[:, k]
        top_idx = _select_top_n(scores_k, n_genes_user)
        names_out[:, k] = gene_names[top_idx]
        raw_top = scores_k[top_idx].astype(np.float32)
        if iqr_normalize:
            # Clamp negatives to 0 before log1p (consistent with sparse path)
            scores_out[:, k] = np.log1p(
                np.maximum(raw_top, 0.0) / iqr_per_group[k]
            )
        else:
            scores_out[:, k] = raw_top

    if verbose:
        print(f"[cytome_cpu] Post-processing: {time.time() - t_post:.3f}s")
        print(f"[cytome_cpu] TOTAL: {time.time() - t0:.3f}s")

    if _opened_here:
        ds.close()
    _nd_result = {
        'names': names_out,
        'scores': scores_out,
        'groups_order': groups_order,
    }
    if write_to_cytome:
        _persist_cosg_markers(cytome_path, _nd_result, groupby,
                              modality, n_genes_user, key_added)
    return _nd_result


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
):
    """Chromosome-batched COSG for ATAC data with millions of features.

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
    from ._cpu_streaming import _apply_cosg_penalty_cpu

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
    group_norms = np.sqrt(group_sizes)
    group_norms_safe = np.maximum(group_norms, 1e-12)

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
    feature_table, idx_col, name_col = _get_feature_table_info(ds, modality)
    chrom_ranges, chrom_order = _get_chrom_col_ranges(ds, feature_table, idx_col)
    n_genes = sum(ce - cs for cs, ce in chrom_ranges.values())

    # Resolve layer + normalizer once for the entire batched accumulation
    _layer_to_read, _need_normalizer = _resolve_layer_to_read(
        ds, modality, layer, compute_on_fly,
    )
    _chunk_normalizer = (
        _resolve_chunk_normalizer(layer, ds, modality, use_cached_stats)
        if _need_normalizer else None
    )

    track_nnz = remove_lowly_expressed and expressed_pct > 0

    if verbose:
        print(f"[cosg_{mode}] {len(labels):,} cells × {n_genes:,} features, "
              f"{n_groups} groups, {len(chrom_order)} chromosomes")
        _emit_filter_log(
            f"cosg_{mode}",
            remove_lowly_expressed=remove_lowly_expressed,
            expressed_pct=expressed_pct,
            expressed_min_num_cells_in_target_group=
                expressed_min_num_cells_in_target_group,
            mu=mu,
        )

    # Determine chromosome batches based on mode
    if mode == "chromosome":
        # Single pass: all chromosomes in one batch
        chrom_batches = [chrom_order]
    else:
        # Multi-pass: split chromosomes to fit in ~200 MB RAM budget
        max_batch_bytes = 200 * 1024 * 1024
        chrom_batches = []
        current_batch = []
        current_bytes = 0
        for chrom in chrom_order:
            cs, ce = chrom_ranges[chrom]
            n_feat = ce - cs
            # RAM per chrom: dot(n×g×4) + sq(n×4) + nnz(n×g×4)
            # x n_batches: dot/sq carry a batch axis when batch_key is set, so the
            # pass budget has to shrink accordingly or "disk" mode silently uses
            # n_batches times its stated RAM target.
            chrom_bytes = n_batches * (n_feat * n_groups * 4 + n_feat * 4)
            if track_nnz:
                chrom_bytes += n_feat * n_groups * 4
            if current_bytes + chrom_bytes > max_batch_bytes and current_batch:
                chrom_batches.append(current_batch)
                current_batch = [chrom]
                current_bytes = chrom_bytes
            else:
                current_batch.append(chrom)
                current_bytes += chrom_bytes
        if current_batch:
            chrom_batches.append(current_batch)

        if verbose:
            print(f"[cosg_{mode}] Split into {len(chrom_batches)} passes")

    # Output collectors. n_keep is the per-group cap on candidate features;
    # for ndarray (standard) and sparse modes alike, it's n_genes_user.
    n_keep = n_genes_user
    group_candidates = {k: [] for k in range(n_groups)}
    # iqr_collector is needed whenever IQR normalization will be applied.
    iqr_collector = {k: [] for k in range(n_groups)} if iqr_normalize else None

    t_stream_total = 0.0

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
                    X_bin.data = np.ones_like(X_bin.data)
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

            # Compute cosine similarity — per batch when batch_key is set, each
            # with its own feature and group norms, averaged over the batches
            # where the group met the cell threshold. Same construction as the
            # non-feature-batched path.
            if batch_index is None:
                gene_norms = np.maximum(np.sqrt(acc['sq'][0]), 1e-12)
                cosine_sim = acc['dot'][0] / (
                    gene_norms[:, None] * group_norms_safe[None, :]
                )
            else:
                _acc = np.zeros((n_feat, n_groups), dtype=np.float64)
                _valid = np.zeros(n_groups, dtype=np.float64)
                for b in range(n_batches):
                    ok = (grp_count_accum[b] >= batch_cell_number_threshold).astype(np.float64)
                    if not ok.any():
                        continue
                    gn = np.maximum(np.sqrt(acc['sq'][b]), 1e-12)
                    gp = np.maximum(np.sqrt(grp_count_accum[b]), 1e-12)
                    _acc += (acc['dot'][b] / (gn[:, None] * gp[None, :])) * ok[None, :]
                    _valid += ok
                with np.errstate(invalid="ignore", divide="ignore"):
                    _acc /= _valid[None, :]
                cosine_sim = np.nan_to_num(_acc, nan=0.0).astype(np.float32)
            del acc['dot'], acc['sq']

            # Apply COSG penalty
            genexlambda = _apply_cosg_penalty_cpu(cosine_sim, mu)
            del cosine_sim

            # min_total_count filter
            if min_total_count > 0 and 'nnz' in acc:
                total_nnz = acc['nnz'].sum(axis=1)
                genexlambda[total_nnz <= min_total_count, :] = -1

            # expressed_pct filter
            if 'nnz' in acc:
                for k in range(n_groups):
                    n_cells_k = group_sizes[k]
                    if n_cells_k > 0:
                        threshold = max(
                            n_cells_k * expressed_pct,
                            expressed_min_num_cells_in_target_group
                        )
                        genexlambda[acc['nnz'][:, k] < threshold, k] = -1
                del acc['nnz']

            # Extract per-group top-K candidates for this chromosome
            for k in range(n_groups):
                col = genexlambda[:, k]

                # Collect positive scores for global IQR
                if iqr_collector is not None:
                    positive = col[col > 0]
                    if len(positive) > 0:
                        iqr_collector[k].append(positive)

                # Top candidates from this chromosome
                n_positive = int((col > 0).sum())
                n_cand = min(n_keep, n_positive)
                if n_cand == 0:
                    continue

                top_local = np.argpartition(col, -n_cand)[-n_cand:]
                for li in top_local:
                    if col[li] > 0:
                        global_idx = cs + int(li)
                        group_candidates[k].append(
                            (global_idx, float(col[li]))
                        )

            del genexlambda

            # Trim candidates to prevent unbounded growth
            for k in range(n_groups):
                if len(group_candidates[k]) > 2 * n_keep:
                    group_candidates[k].sort(
                        key=lambda x: x[1], reverse=True
                    )
                    group_candidates[k] = group_candidates[k][:n_keep]

        del chrom_accum

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
        q_up = float(np.quantile(full_col, q_upper))
        q_lo = float(np.quantile(full_col, q_lower))
        iqr = q_up - q_lo
        if iqr <= 0:
            iqr = q_up if q_up > 0 else 1e-6
        return iqr

    if use_sparse_output:
        scores_dict = {}
        for k in range(n_groups):
            if not group_candidates[k]:
                continue

            # All-values IQR (matches cosg.iqrLogNormalize on AnnData side)
            iqr = (
                _compute_group_iqr_from_collector(k)
                if iqr_normalize else None
            )

            # Merge candidates across chromosomes, keep top-K globally
            candidates = group_candidates[k]
            candidates.sort(key=lambda x: x[1], reverse=True)
            top_candidates = candidates[:n_keep]

            # Lazy feature name lookup
            top_indices = [c[0] for c in top_candidates]
            top_scores_raw = [c[1] for c in top_candidates]
            top_names = _lookup_feature_names(
                ds, feature_table, idx_col, name_col, top_indices
            )

            group_name = str(groups_order[k])
            for name, score_raw in zip(top_names, top_scores_raw):
                if iqr_normalize:
                    score = float(np.log1p(score_raw / iqr))
                else:
                    score = float(score_raw)
                if score > 0:
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
        result = _format_top_k_output(
            scores_dict, groups_order_str, group_size_dict,
            output_format=output_format,
            ds=ds, modality=modality, verbose=verbose,
            feature_name_col=feature_name_col,
        )
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
        names_out = np.empty((n_genes_user, n_groups), dtype=object)
        scores_out = np.empty((n_genes_user, n_groups), dtype=np.float32)

        for k in range(n_groups):
            candidates = group_candidates[k]
            candidates.sort(key=lambda x: x[1], reverse=True)
            top_candidates = candidates[:n_genes_user]

            if not top_candidates:
                names_out[:, k] = ""
                scores_out[:, k] = 0.0
                continue

            top_indices = [c[0] for c in top_candidates]
            top_scores_vals = np.array(
                [c[1] for c in top_candidates], dtype=np.float32
            )
            if iqr_normalize:
                # Same all-values IQR formula as sparse path
                iqr = _compute_group_iqr_from_collector(k)
                top_scores_vals = np.log1p(
                    np.maximum(top_scores_vals, 0.0) / iqr
                ).astype(np.float32)

            top_names = _lookup_feature_names(
                ds, feature_table, idx_col, name_col, top_indices
            )

            n_fill = len(top_candidates)
            names_out[:n_fill, k] = top_names
            scores_out[:n_fill, k] = top_scores_vals
            if n_fill < n_genes_user:
                names_out[n_fill:, k] = ""
                scores_out[n_fill:, k] = 0.0

        if verbose:
            print(f"[cosg_{mode}] Post: {time.time() - t_post:.2f}s, "
                  f"Total: {time.time() - t0:.2f}s")

        if _opened_here:
            ds.close()
        return {
            'names': names_out,
            'scores': scores_out,
            'groups_order': groups_order,
        }


def _run_cosg_cytome_gpu_impl(
    ds,                                 # open cytome.Dataset
    groupby,
    *,
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
    feature_name_col=None,
    batch_key=None,
    batch_cell_number_threshold=3,
):
    """GPU implementation body — accepts an open ``cytome.Dataset``.

    Internal helper so both ``run_cosg_cytome(..., device='gpu')`` and
    the back-compat shim ``run_cosg_cytome_gpu(...)`` share one code
    path. The caller is responsible for opening/closing the Dataset.
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

    n_genes = len(gene_names)
    n_groups = len(groups_order)

    group_sizes_gpu = cp.array(group_sizes, dtype=cp.float32)
    group_norms_gpu = cp.sqrt(group_sizes_gpu)

    if verbose:
        print(f"[cytome_gpu] {len(labels):,} cells × {n_genes:,} genes, "
              f"{n_groups} groups")

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

    # Initialize GPU accumulators
    dot_accum = cp.zeros((n_batches, n_genes, n_groups), dtype=cp.float32)
    gene_sq_accum = cp.zeros((n_batches, n_genes), dtype=cp.float32)
    grp_count_accum = cp.zeros((n_batches, n_groups), dtype=cp.float64)

    if remove_lowly_expressed and expressed_pct > 0:
        nnz_accum = cp.zeros((n_genes, n_groups), dtype=cp.float32)
    else:
        nnz_accum = None

    # Resolve layer + normalizer once before streaming (same logic as CPU path)
    _layer_to_read, _need_normalizer = _resolve_layer_to_read(
        ds, modality, layer, compute_on_fly,
    )
    _chunk_normalizer = (
        _resolve_chunk_normalizer(layer, ds, modality, use_cached_stats)
        if _need_normalizer else None
    )

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
            del Xsq
            for gi in np.asarray(chunk_groups):
                grp_count_accum[0, int(gi)] += 1.0
        else:
            chunk_batches = batch_index[chunk_labels_idx]
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
                for gi in np.asarray(chunk_groups)[sel]:
                    grp_count_accum[int(b), int(gi)] += 1.0
                del Xb, Lb, Db, Xsq

        # Nonzero counts
        if nnz_accum is not None:
            X_bin = X_gpu.copy()
            X_bin.data = cp.ones_like(X_bin.data)
            nnz_chunk = Lam_dense_T @ X_bin
            if cp_sparse.issparse(nnz_chunk):
                nnz_chunk = nnz_chunk.toarray()
            nnz_accum += nnz_chunk.T
            del X_bin, nnz_chunk

        cp.cuda.Stream.null.synchronize()
        t_compute += time.time() - t_comp

        # D_chunk only exists on the non-batch branch; the batch branch builds
        # and frees a per-batch Db inside its loop.
        if batch_index is None:
            del D_chunk
        del X_gpu, Lam_dense_T
        cp.get_default_memory_pool().free_all_blocks()

        n_chunks += 1
        n_cells_seen += actual_chunk

        if verbose and n_chunks % 500 == 0:
            print(f"  Processed {n_chunks} chunks, {n_cells_seen:,} cells")

    if verbose:
        print(f"[cytome_gpu] Transfer: {t_transfer:.3f}s, "
              f"Compute: {t_compute:.3f}s")

    # Finalize on GPU
    t_post = time.time()

    group_norms_safe = cp.maximum(group_norms_gpu, 1e-12)

    if batch_index is None:
        gene_norms = cp.maximum(cp.sqrt(gene_sq_accum[0]), 1e-12)
        cosine_sim = dot_accum[0] / (gene_norms[:, None] * group_norms_safe[None, :])
    else:
        # Per batch, with that batch's own gene and group norms, averaged over the
        # batches where the group met the cell threshold — identical construction
        # to the CPU path, just on device.
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

    # Expressed pct filter on GPU
    if nnz_accum is not None:
        for k in range(n_groups):
            n_cells_k = group_sizes_gpu[k]
            if n_cells_k > 0:
                threshold = max(
                    float(n_cells_k * expressed_pct),
                    float(expressed_min_num_cells_in_target_group)
                )
                genexlambda_gpu[nnz_accum[:, k] < threshold, k] = -1.0
        del nnz_accum

    genexlambda = cp.asnumpy(genexlambda_gpu)
    del genexlambda_gpu
    cp.get_default_memory_pool().free_all_blocks()

    # Extract top genes per group (CPU)
    names_out = np.empty((n_genes_user, n_groups), dtype=object)
    scores_out = np.empty((n_genes_user, n_groups), dtype=np.float32)

    for k in range(n_groups):
        scores_k = genexlambda[:, k]
        top_idx = _select_top_n(scores_k, n_genes_user)
        names_out[:, k] = gene_names[top_idx]
        scores_out[:, k] = scores_k[top_idx].astype(np.float32)

    if verbose:
        print(f"[cytome_gpu] Post-processing: {time.time() - t_post:.3f}s")
        print(f"[cytome_gpu] TOTAL: {time.time() - t0:.3f}s")

    return {
        'names': names_out,
        'scores': scores_out,
        'groups_order': groups_order,
    }


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

    Same kwargs and return shape as the GPU path inside
    ``run_cosg_cytome``. Currently supports only the ndarray output
    shape; for ``output_format='dict' / 'long' / 'dense'`` use the
    CPU path.

    This is a thin shim around the unified ``run_cosg_cytome`` entry
    point — the body lives in ``_run_cosg_cytome_gpu_impl`` and is
    shared with ``run_cosg_cytome(..., device='gpu')``.

    Parameters
    ----------
    modality : str
        Modality to use (``"RNA"``, ``"ATAC"``, ``"GA"``, ``"tiles"``).
        Routes via the cytome modality registry to the matching feature
        entity (RNA → ``genes``, ATAC → ``peaks``, GA → ``GA_genes``,
        tiles → ``tiles``) and matrix layer.
    layer : str, default ``"counts"``
        Measurement layer to read for COSG scoring. Resolved via the
        same 4-way dispatch as ``run_cosg_cytome``:

        - ``"counts"``: raw counts, normalized on-the-fly to log1p-CPM.
        - ``"log1p"``: pre-computed log1p layer; used as-is.
        - ``"infog"``: pre-computed INFOG-weighted log1p layer; used
          as-is.
        - ``"tfidf"``: TF-IDF (ATAC / tiles only); cached params
          loaded from ``ds.metadata`` or computed on first use.
        - Any other materialized layer name is read as-is from
          ``ds.{modality}.layer(name)``.

    Migration notes — passing the old ``top_k`` or ``n_top_genes``
    kwargs raises ``TypeError`` with a migration hint:

        - ``n_top_genes=N`` → ``n_genes_user=N``
        - ``top_k=N``       → not supported on GPU; use
          ``run_cosg_cytome(..., n_genes_user=N, output_format='dict')``.
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


def _select_top_n(scores, n_top):
    """Select top n gene indices by score (descending)."""
    reference_indices = np.arange(scores.shape[0], dtype=int)
    n_top = min(n_top, len(scores))
    partition = np.argpartition(scores, -n_top)[-n_top:]
    partial_indices = np.argsort(scores[partition])[::-1]
    global_indices = reference_indices[partition][partial_indices]
    return global_indices
