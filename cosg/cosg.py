from anndata import AnnData
import numpy as np
import pandas as pd
from scipy import sparse
from typing import Iterable, Union, Optional



### Refer to: https://github.com/theislab/scanpy/blob/5533b644e796379fd146bf8e659fd49f92f718cd/scanpy/_compat.py
try:
    from typing import Literal
except ImportError:
    try:
        from typing_extensions import Literal
    except ImportError:

        class LiteralMeta(type):
            def __getitem__(cls, values):
                if not isinstance(values, tuple):
                    values = (values,)
                return type('Literal_', (Literal,), dict(__args__=values))

        class Literal(metaclass=LiteralMeta):
            pass

### Refer to Scanpy       
def _select_top_n(scores, n_top):
    reference_indices = np.arange(scores.shape[0], dtype=int)
    partition = np.argpartition(scores, -n_top)[-n_top:]
    partial_indices = np.argsort(scores[partition])[::-1]
    global_indices = reference_indices[partition][partial_indices]
    return global_indices

from scipy.sparse import csr_matrix
from sklearn.metrics.pairwise import cosine_similarity


_CYTOME_DEFAULT = object()   # sentinel for cytome-only kwargs in cosg()


def _is_anndata(obj):
    """True if ``obj`` is an AnnData (lazy import; AnnData may not be
    installed in some minimal cytome-only environments)."""
    try:
        import anndata
    except ImportError:
        return False
    return isinstance(obj, anndata.AnnData)


def _is_cytome_input(obj):
    """True if ``obj`` looks like a cytome input (path or open Dataset).

    Defers to ``cosg._cytome_streaming._is_cytome_dataset`` for the
    Dataset check, plus accepts ``str`` and ``pathlib.Path`` for paths.
    """
    import os
    from pathlib import Path
    if isinstance(obj, (str, Path)):
        return True
    try:
        from ._cytome_streaming import _is_cytome_dataset
    except ImportError:
        return False
    return _is_cytome_dataset(obj)


# AnnData-only kwargs the polymorphic ``cosg()`` recognises. If any of
# these is set to a NON-default value while the input is a cytome,
# the dispatcher raises TypeError naming the kwarg.
#
# Round 9 (2026-05-23): ``layer`` removed from this list. It's now a
# SHARED kwarg with the cytome side — same name on both, but means
# "adata.layers key" on AnnData and "cytome layer name" on cytome.
# Round 30 (2026-07-31): ``batch_key`` and ``batch_cell_number_threshold`` moved
# out of this list — the cytome streaming path implements them now (per-batch
# cosine with that batch's own norms, averaged over batches meeting the cell
# threshold), so they are SHARED kwargs like ``layer``.
_ANNDATA_ONLY_KWARGS = (
    "groups", "key_added",
    "calculate_logfoldchanges", "use_raw", "reference",
    "return_by_group", "column_delimiter", "copy", "gpu_chunk_size",
    "cpu_chunk_size",
)

# Cytome-only kwargs the polymorphic ``cosg()`` recognises. If any of
# these is non-sentinel while the input is an AnnData, dispatcher raises.
#
# Round 9: ``cytome_layer`` removed — replaced by shared ``layer``.
_CYTOME_ONLY_KWARGS = (
    "modality", "cell_mask", "batch_size", "feature_batching",
    "min_total_count", "output_format", "compute_on_fly",
    "use_cached_stats", "iqr_normalize", "q_upper", "q_lower",
)


def _dispatch_to_cytome(
    cytome_input,
    *,
    groupby, n_genes_user, mu,
    remove_lowly_expressed, expressed_pct,
    expressed_min_num_cells_in_target_group,
    device, verbosity,
    batch_key, batch_cell_number_threshold,   # Round 30: shared with the cytome path
    layer,                              # Round 9: shared kwarg, forward as layer=
    _anndata_only_values,
    # cytome-only kwargs (each may be _CYTOME_DEFAULT):
    modality, cell_mask, batch_size, feature_batching, min_total_count,
    output_format, compute_on_fly, use_cached_stats,
    iqr_normalize, q_upper, q_lower,
):
    """Internal: route from ``cosg.cosg(cytome_input, ...)`` to
    ``run_cosg_cytome``.

    Strict kwarg discipline:

    * Rejects AnnData-only kwargs (key_added, copy, use_raw, batch_key,
      calculate_logfoldchanges, gpu/cpu_chunk_size, groups, reference,
      return_by_group, column_delimiter, batch_cell_number_threshold)
      with TypeError naming the offender.
    * Forwards cytome-only kwargs only when they were explicitly set
      (so cytome defaults stay authoritative).
    * Maps AnnData-style ``verbosity: int`` → cytome-style
      ``verbose: bool`` (`verbosity > 0`).
    * ``layer`` is shared between AnnData and cytome paths: forwarded
      only when the user set it explicitly (i.e. not the AnnData
      default of None) so the cytome path's own default ('auto')
      wins when omitted.
    """
    from ._cytome_streaming import run_cosg_cytome

    # Reject AnnData-only kwargs set to non-default values.
    _anndata_defaults = {
        "groups": "all", "key_added": None,
        "calculate_logfoldchanges": False, "use_raw": False,
        "reference": "rest", "return_by_group": True,
        "column_delimiter": "::", "copy": False, "gpu_chunk_size": None,
        "cpu_chunk_size": None,
    }
    _bad_anndata_kwargs = [
        k for k, v in _anndata_only_values.items()
        if v != _anndata_defaults[k]
    ]
    if _bad_anndata_kwargs:
        bad = ", ".join(sorted(_bad_anndata_kwargs))
        raise TypeError(
            f"cosg.cosg(cytome_input, ...) does not accept AnnData-only "
            f"kwarg(s): {bad}. Drop the kwarg(s), or pass an AnnData as "
            f"the first argument instead."
        )

    # Build the forward kwargs — only set those the user explicitly chose.
    _fwd = dict(
        groupby=groupby,
        n_genes_user=n_genes_user,
        mu=mu,
        remove_lowly_expressed=remove_lowly_expressed,
        expressed_pct=expressed_pct,
        expressed_min_num_cells_in_target_group=
            expressed_min_num_cells_in_target_group,
        device=device,
        verbose=(verbosity > 0),
    )
    # Shared ``layer`` kwarg: forward only when the user set it (i.e.
    # not the AnnData default of None) so the cytome path's own default
    # of 'auto' wins when omitted.
    if layer is not None:
        _fwd["layer"] = layer
    # Shared batch kwargs (Round 30): forward only when set, so the cytome
    # defaults stay authoritative.
    if batch_key is not None:
        _fwd["batch_key"] = batch_key
        _fwd["batch_cell_number_threshold"] = batch_cell_number_threshold
    for _name, _val in (
        ("modality", modality), ("cell_mask", cell_mask),
        ("batch_size", batch_size), ("feature_batching", feature_batching),
        ("min_total_count", min_total_count),
        ("output_format", output_format),
        ("compute_on_fly", compute_on_fly),
        ("use_cached_stats", use_cached_stats),
        ("iqr_normalize", iqr_normalize),
        ("q_upper", q_upper), ("q_lower", q_lower),
    ):
        if _val is not _CYTOME_DEFAULT:
            _fwd[_name] = _val

    return run_cosg_cytome(cytome_input, **_fwd)


def cosg(
    adata,
    groupby: str = 'CellTypes',
    groups: Union[Literal['all'], Iterable[str]] = 'all',
    mu: float = 1,
    remove_lowly_expressed: bool = True,
    expressed_pct: Optional[float] = 0.1,
    expressed_min_num_cells_in_target_group: int = 3,
    n_genes_user: int = 50,
    batch_key: str = None,
    batch_cell_number_threshold: int = 3,
    key_added: Optional[str] = None,
    calculate_logfoldchanges: bool = False,
    use_raw: bool = False,
    layer: Optional[str] = None,
    reference: str = 'rest',
    return_by_group : bool = True,
    column_delimiter: str = "::",
    verbosity: int = 0,
    copy: bool = False,
    device: str = 'cpu',
    gpu_chunk_size: Optional[int] = None,
    cpu_chunk_size: Optional[int] = None,
    # ---- Cytome-only kwargs (Round 8, 2026-05-23) ----
    # These are recognised when ``adata`` is a cytome path / Dataset.
    # Passing any of them with a non-sentinel value while the input is
    # an AnnData raises ``TypeError``.
    # Round 9 (2026-05-23): ``cytome_layer`` removed — use shared ``layer``.
    modality=_CYTOME_DEFAULT,
    cell_mask=_CYTOME_DEFAULT,
    batch_size=_CYTOME_DEFAULT,
    feature_batching=_CYTOME_DEFAULT,
    min_total_count=_CYTOME_DEFAULT,
    output_format=_CYTOME_DEFAULT,
    compute_on_fly=_CYTOME_DEFAULT,
    use_cached_stats=_CYTOME_DEFAULT,
    iqr_normalize=_CYTOME_DEFAULT,
    q_upper=_CYTOME_DEFAULT,
    q_lower=_CYTOME_DEFAULT,
):
    """\
    Marker gene identification for single-cell sequencing data using COSG.
    
    Parameters
    ----------
    adata
        AnnData object with cell group information.
    groupby : str, default 'CellTypes'
        Key in `adata.obs` that defines the cell group labels. Defaults to 'CellTypes'.
    groups : {'all'} or iterable of str, default 'all'
        Subset of cell groups, e.g. [`'g1'`, `'g2'`, `'g3'`], to which comparison shall be restricted. The default value is 'all', and all groups will be compared.
    mu : float, default 1
        The penalty parameter restricting marker genes expressing in non-target cell groups. Larger value represents more strict restrictions. mu should be non-negative, and by default, mu = 1.
    remove_lowly_expressed : bool, default True
        If True, genes that express a percentage of target cells smaller than a specific value (`expressed_pct`) are not considered as marker genes for the target cells.

        .. note::
           Default changed from ``False`` to ``True`` to align with the
           cytome-streaming variant (``run_cosg_cytome``), which has
           always defaulted to ``True``. The 10% expression filter is the
           recommended default for marker detection across both backends —
           pass ``remove_lowly_expressed=False`` explicitly if you need
           the previous unfiltered behaviour.
    expressed_pct : float, optional, default 0.1
        When `remove_lowly_expressed` is set to True, genes that express a percentage of target cells smaller than a specific value (`expressed_pct`) are not considered as marker genes for the target cells. The default value for `expressed_pct` is 0.1 (10%).
    expressed_min_num_cells_in_target_group : int, default 3
        When `remove_lowly_expressed` is set to True, this sets the minimum absolute 
        number of cells that must express a gene for it to be considered as a marker.
        The actual threshold used is the maximum of (`n_cells * expressed_pct`) and 
        this value. This prevents overly permissive thresholds in small clusters.
        For example, with a cluster of 10 cells and `expressed_pct=0.1`, the threshold
        would be max(1, 3) = 3 cells instead of just 1 cell.
    n_genes_user : int, default 50
        The number of genes that appear in the returned tables. The default value is 50.
    batch_key : str, optional, default None
        Used to indicate which batch info in `adata.obs` to use for calculate the cosine similarities separately for each batch
        and then average them. The default value is None.
    batch_cell_number_threshold : int, default 3
        Minimum number of cells required in a batch for a given cluster to be considered
        for computing the cosine similarity score when `batch_key` is not None. If a cluster has fewer than this number of cells
        in a batch, the cosine similarity from that batch for the cluster will be ignored in the average.
    key_added : str, optional, default None
        Key under which the COSG results will be stored in `adata.uns`. If None, the default key 'cosg' is used.
    calculate_logfoldchanges : bool, default False
        If True, computes log fold-changes for the marker genes.
    use_raw : bool, default False
        If True and `adata.raw` exists, raw UMI counts saved in `adata.raw.X` are used for calculations.
    layer : str, optional, default None
        Key in `adata.layers` to specify an alternative gene expression layer for COSG.
    reference : str, default 'rest'
        Specifies the reference group for comparison. If `'rest'`, compare each group to the union of the rest of the group.
        If a group identifier, compare with respect to this group.
    return_by_group : bool, default True
        If True, the COSG scores are also summarized and stored separately by group in adata.uns[`key_added`]['COSG']. This will output another extra copy of the results. Defaults to True.
    column_delimiter : str, default "::"
        Delimiter to use when combining multi-index column names for `adata.write()` compatibility. Default is "::".
    verbosity : int, default 0
        Controls the verbosity of logging messages, defaults to 0. Higher values produce more detailed messages.
    copy: bool, default False
        If True, returns a copy of `adata` with the computed embeddings instead of modifying in place. Defaults to False.
    device : str, default 'cpu'
        Computation backend. Options:

        - ``'cpu'``: Use CPU (NumPy/SciPy). Default behavior, no extra dependencies.
        - ``'gpu'``: Use GPU (CuPy/cuSPARSE). Requires ``pip install cosg[gpu]``.
        - ``'auto'``: Use GPU if CuPy is available, otherwise fall back to CPU.
    gpu_chunk_size : int or None, default None
        Controls GPU streaming behavior for memory-efficient GPU processing.

        - None (default): Auto-select based on available GPU memory. Uses
          monolithic transfer when VRAM is sufficient (~15 GB), falls back
          to chunked streaming (~1 GB) when VRAM is tight.
        - 0 or n_cells: Force monolithic GPU transfer (requires ~15+ GB VRAM).
        - N (positive int): Use N cells per GPU chunk.
    cpu_chunk_size : int or None, default None
        Controls CPU streaming behavior for memory-efficient processing.

        - None (default): Auto-select chunk size based on available RAM.
          Uses chunked cell-axis accumulation that processes cells in blocks,
          reducing peak memory by ~2x while improving speed by ~1.8x.
        - 0: Use legacy CPU path (non-streaming). Required for batch_key mode.
          Slower and uses more memory, but produces bit-identical results to
          COSG versions prior to v1.2.
        - N (positive int): Use N cells per chunk. Smaller chunks use less
          memory but may be slightly slower due to per-chunk overhead.

    Returns
    -------
    None or AnnData
        The function stores the COSG marker gene identification results in `adata.uns[key_added]` (if `copy` is False)
        or returns a modified copy of `adata` (if `copy` is True). The results include structured arrays containing 
        gene names and their corresponding COSG scores for each cell group:
        names : structured `np.ndarray` (`.uns['cosg']`)
            Structured array to be indexed by group id storing the gene names. Ordered according to scores.
        scores : structured `np.ndarray` (`.uns['cosg']`)
            Structured array to be indexed by group id storing COSG scores for each gene for each
                group. Ordered according to scores.
        The marker genes and their COSG scores are also summarized and stored separately by group in
        adata.uns[`key_added`]['COSG'] if `return_by_group = True`.
    
    Notes
    -----
    COSG does not modify the input AnnData object's expression data (adata.X).
    It is safe to call cosg() without copying adata first.

    The function expects log-normalized expression data in adata.X (the output
    of sc.pp.normalize_total + sc.pp.log1p). Raw counts can be used via
    use_raw=True if adata.raw is set.

    GPU acceleration requires CuPy: pip install cupy-cuda12x

    Examples
    --------
    >>> import cosg as cosg
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> cosg.cosg(adata, key_added='cosg', groupby='bulk_labels')
    >>> ### Visualize the marker genes in a dot plot
    >>> cosg.plotMarkerDotplot(
    ...     adata,
    ...     groupby='bulk_labels',
    ...     top_n_genes=3,
    ...     key_cosg='cosg',
    ...     use_rep=None,
    ...     swap_axes=False,
    ...     standard_scale='var',
    ...     cmap='Spectral_r'
    ... )
    """
    # Round 8 (2026-05-23): polymorphic dispatch on the first argument.
    # When ``adata`` is a cytome path (str / Path) or an open cytome
    # Dataset, delegate to ``run_cosg_cytome``. Otherwise fall through
    # to the existing AnnData implementation below — that body is
    # unchanged from pre-Round-8.
    # Round 9 (2026-05-23): forward shared ``layer`` kwarg; no
    # ``cytome_layer`` (collapsed into ``layer``).
    if _is_cytome_input(adata):
        return _dispatch_to_cytome(
            adata,
            groupby=groupby,
            n_genes_user=n_genes_user,
            mu=mu,
            remove_lowly_expressed=remove_lowly_expressed,
            expressed_pct=expressed_pct,
            expressed_min_num_cells_in_target_group=
                expressed_min_num_cells_in_target_group,
            device=device,
            verbosity=verbosity,
            layer=layer,
            batch_key=batch_key,
            batch_cell_number_threshold=batch_cell_number_threshold,
            # AnnData-only kwargs — reject if any non-default.
            _anndata_only_values=dict(
                groups=groups,
                key_added=key_added,
                calculate_logfoldchanges=calculate_logfoldchanges,
                use_raw=use_raw,
                reference=reference,
                return_by_group=return_by_group,
                column_delimiter=column_delimiter,
                copy=copy,
                gpu_chunk_size=gpu_chunk_size,
                cpu_chunk_size=cpu_chunk_size,
            ),
            # Cytome-only kwargs — forward when non-sentinel.
            modality=modality,
            cell_mask=cell_mask,
            batch_size=batch_size,
            feature_batching=feature_batching,
            min_total_count=min_total_count,
            output_format=output_format,
            compute_on_fly=compute_on_fly,
            use_cached_stats=use_cached_stats,
            iqr_normalize=iqr_normalize,
            q_upper=q_upper,
            q_lower=q_lower,
        )

    if not _is_anndata(adata):
        raise TypeError(
            f"cosg.cosg expected an AnnData, str, pathlib.Path, or "
            f"cytome Dataset as the first argument; got "
            f"{type(adata).__name__!r}."
        )

    # AnnData path — reject any cytome-only kwarg that was set.
    # Round 9 (2026-05-23): ``cytome_layer`` removed from this list —
    # the shared ``layer`` kwarg handles both backends now.
    _cytome_set = {
        k: v for k, v in {
            "modality": modality, "cell_mask": cell_mask,
            "batch_size": batch_size, "feature_batching": feature_batching,
            "min_total_count": min_total_count,
            "output_format": output_format,
            "compute_on_fly": compute_on_fly,
            "use_cached_stats": use_cached_stats,
            "iqr_normalize": iqr_normalize,
            "q_upper": q_upper, "q_lower": q_lower,
        }.items() if v is not _CYTOME_DEFAULT
    }
    if _cytome_set:
        bad = ", ".join(sorted(_cytome_set))
        raise TypeError(
            f"cosg.cosg(AnnData, ...) does not accept cytome-only "
            f"kwarg(s): {bad}. Drop the kwarg(s), or pass a cytome "
            f"path / open cytome Dataset as the first argument instead."
        )

    ### Validate the input parameter mu
    if mu < 0:
        raise ValueError("Parameter mu must be non-negative.")

    ### Resolve device backend
    from ._backend import get_device
    _device = get_device(device)
    if verbosity > 0 and _device == 'gpu':
        print("[cosg] Using GPU (CuPy) backend.")
        
        
    # Validate groupby column exists
    if groupby not in adata.obs.columns:
        raise ValueError(
            f"Column '{groupby}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )
        
    # Validate groups parameter
    if groups != 'all':
        if isinstance(groups, (str, int)):
            raise ValueError("Specify a sequence of groups, not a single value.")
        
        available_groups = set(adata.obs[groupby].unique())
        requested_groups = set(groups)
        invalid_groups = requested_groups - available_groups
        
        if invalid_groups:
            raise ValueError(
                f"Groups {invalid_groups} not found in adata.obs['{groupby}']. "
                f"Available groups: {available_groups}"
            )
            
        if batch_key is not None:
            raise NotImplementedError(
                "Using `batch_key` with a subset of `groups` is not currently supported. "
                "Please use `groups='all'` when using batch correction, or set `batch_key=None`."
            )
            
    # Validate n_genes_user
    if n_genes_user < 1:
        raise ValueError("n_genes_user must be at least 1.")
    
    # Validate expressed_pct
    if expressed_pct is not None and not (0 <= expressed_pct <= 1):
        raise ValueError("expressed_pct must be between 0 and 1.")
        
    # Validate expressed_min_num_cells_in_target_group
    if expressed_min_num_cells_in_target_group < 1:
        raise ValueError("expressed_min_num_cells_in_target_group must be at least 1.")
    
    # Validate layer exists if specified
    if layer is not None and layer not in adata.layers:
        raise ValueError(
            f"Layer '{layer}' not found in adata.layers. "
            f"Available layers: {list(adata.layers.keys())}"
        )
    
    adata = adata.copy() if copy else adata

    if layer is not None:
        if use_raw:
            raise ValueError("Cannot specify `layer` and have `use_raw = True`.")
        cellxgene = adata.layers[layer]
    else:
        if use_raw and adata.raw is not None:
             cellxgene = adata.raw.X
        else:
            cellxgene = adata.X
    
    
    ### Refer to scanpy's framework
    # https://github.com/theislab/scanpy/blob/5533b644e796379fd146bf8e659fd49f92f718cd/scanpy/tools/_rank_genes_groups.py#L559
    if key_added is None:
        key_added = 'cosg'
    adata.uns[key_added] = {}
    adata.uns[key_added]['params'] = dict(
        groupby=groupby,
        reference=reference,
        groups=groups,
        batch_key=batch_key,
        method='COSG',
        use_raw=use_raw,
        layer=layer,
        mu=mu,
        remove_lowly_expressed=remove_lowly_expressed,  
        expressed_pct=expressed_pct, 
        expressed_min_num_cells_in_target_group=expressed_min_num_cells_in_target_group, 
    )
    
    ### Refer to: https://github.com/theislab/scanpy/blob/5533b644e796379fd146bf8e659fd49f92f718cd/scanpy/tools/_rank_genes_groups.py#L543
    if groups == 'all':
        ### group lable for each cell
        group_info=adata.obs[groupby]
    elif isinstance(groups, (str, int)):
        raise ValueError('Specify a sequence of groups')
    else:
        cells_selected=adata.obs[groupby].isin(groups)
        cells_selected=cells_selected.values
        if sparse.issparse(cellxgene):
            cellxgene=cellxgene[cells_selected]
        else:
            cellxgene=cellxgene[cells_selected,:]
            
            

        ### group lable for each cell
        group_info=adata.obs[groupby].copy()
        group_info=group_info[cells_selected]
        

    
    ### Use the categorical orders if it is a categorical variable
    if hasattr(group_info, "cat"):
        groups_order = list(group_info.cat.categories)
    else:
        unique_values=group_info.unique()
        ### add a helper function here, if the cell clusters are "1", "2", ... , "10", "11", ...
        ### order them as "1", "2", ... , "10", "11", ..., instead of being "1", "10", "11", ..., "2", ...
        def _is_all_numeric(groups):
            try:
                [float(x) for x in groups]
                return True
            except ValueError:
                return False

        if _is_all_numeric(unique_values):
            groups_order = sorted(unique_values, key=lambda x: float(x))
        else:
            groups_order = sorted(unique_values)
        
    
    
    n_cluster=len(groups_order)

    n_cell=cellxgene.shape[0]

    ### Efficiently create a sparse matrix for the cluster_mat matrix
    group_to_row = {group: i for i, group in enumerate(groups_order)}
    row_indices = np.array([group_to_row[group] for group in group_info])
    col_indices = np.arange(n_cell)
    data = np.ones_like(col_indices, dtype=int)
    cluster_mat = csr_matrix((data, (row_indices, col_indices)), shape=(n_cluster, n_cell))

    # ── GPU path ─────────────────────────────────────────────────────────
    _gpu_full_path = False  # flag: True when gpu_core_cosg_full() is used
    _gpu_chunked_path = False  # flag: True when gpu_core_cosg_chunked() is used
    _cpu_streaming_used = False  # flag: True when cpu_core_cosg_chunked() is used
    if _device == 'gpu':
        from ._gpu import (gpu_core_cosg, gpu_core_cosg_batched,
                           gpu_core_cosg_chunked, _should_use_chunked)

        if batch_key is not None:
            # Validate batch_key exists
            if batch_key not in adata.obs.columns:
                raise ValueError(
                    f"batch_key '{batch_key}' not found in adata.obs. "
                    f"Available columns: {list(adata.obs.columns)}"
                )
            batch_info = adata.obs[batch_key].values
            unique_batches = np.unique(batch_info)

            genexlambda = gpu_core_cosg_batched(
                cellxgene, cluster_mat, mu,
                batch_info=batch_info,
                unique_batches=unique_batches,
                batch_cell_number_threshold=batch_cell_number_threshold,
                is_sparse_input=sparse.issparse(cellxgene),
                verbosity=verbosity
            )
        else:
            # Decide: chunked streaming vs monolithic
            use_chunked = (gpu_chunk_size is not None or
                           _should_use_chunked(n_cell, cellxgene.shape[1],
                                               sparse.issparse(cellxgene), cellxgene))
            if use_chunked:
                genexlambda = gpu_core_cosg_chunked(
                    cellxgene, cluster_mat, mu,
                    is_sparse_input=sparse.issparse(cellxgene),
                    remove_lowly_expressed=remove_lowly_expressed,
                    expressed_pct=expressed_pct if expressed_pct is not None else 0.1,
                    expressed_min_num_cells_in_target_group=expressed_min_num_cells_in_target_group,
                    chunk_size=gpu_chunk_size,
                    verbosity=verbosity
                )
                _gpu_chunked_path = True
            else:
                genexlambda = gpu_core_cosg(
                    cellxgene, cluster_mat, mu,
                    is_sparse_input=sparse.issparse(cellxgene),
                    remove_lowly_expressed=remove_lowly_expressed,
                    expressed_pct=expressed_pct if expressed_pct is not None else 0.1,
                    expressed_min_num_cells_in_target_group=expressed_min_num_cells_in_target_group,
                    verbosity=verbosity
                )
    
    # ── CPU path ──────────────────────────────────────────────────────────
    else:
        _cpu_streaming_used = False

        # CPU streaming path (non-batch)
        if cpu_chunk_size != 0 and batch_key is None:
            from ._cpu_streaming import cpu_core_cosg_chunked

            genexlambda = cpu_core_cosg_chunked(
                cellxgene, cluster_mat, mu,
                is_sparse_input=sparse.issparse(cellxgene),
                remove_lowly_expressed=remove_lowly_expressed,
                expressed_pct=expressed_pct if expressed_pct is not None else 0.1,
                expressed_min_num_cells_in_target_group=expressed_min_num_cells_in_target_group,
                chunk_size=cpu_chunk_size,  # None = auto
                verbosity=verbosity
            )
            _cpu_streaming_used = True

        # CPU streaming batch path
        elif batch_key is not None and cpu_chunk_size != 0:
            from ._cpu_streaming import (cpu_cosine_sim_chunked,
                                         _apply_cosg_penalty_cpu)

            if batch_key not in adata.obs.columns:
                raise ValueError(
                    f"batch_key '{batch_key}' not found in adata.obs. "
                    f"Available columns: {list(adata.obs.columns)}"
                )

            batch_info = adata.obs[batch_key].values
            unique_batches = np.unique(batch_info)
            n_gene = cellxgene.shape[1]

            accum = np.zeros((n_gene, n_cluster), dtype=np.float64)
            valid_counts = np.zeros(n_cluster, dtype=np.float64)

            for batch in unique_batches:
                indices = np.flatnonzero(batch_info == batch)
                if indices.size == 0:
                    continue

                cluster_counts = np.array(
                    cluster_mat[:, indices].sum(axis=1)
                ).flatten()
                valid_mask = (
                    cluster_counts >= batch_cell_number_threshold
                ).astype(float)
                if not np.any(valid_mask):
                    continue

                sim_batch = cpu_cosine_sim_chunked(
                    cellxgene[indices],
                    cluster_mat[:, indices],
                    is_sparse_input=sparse.issparse(cellxgene),
                    chunk_size=cpu_chunk_size,
                    verbosity=verbosity
                )

                sim_batch *= valid_mask[np.newaxis, :]
                accum += sim_batch
                valid_counts += valid_mask

            for j in range(n_cluster):
                if valid_counts[j] > 0:
                    accum[:, j] /= valid_counts[j]
                else:
                    accum[:, j] = np.nan

            cosine_sim = accum
            cosine_sim_clean = np.nan_to_num(cosine_sim, nan=0.0)
            genexlambda = _apply_cosg_penalty_cpu(cosine_sim_clean, mu)

            # Apply expressed_pct filter using chunked nnz computation
            if remove_lowly_expressed:
                from ._cpu_streaming import cpu_nnz_counts_chunked
                _ept = expressed_pct if expressed_pct is not None else 0.1

                nnz_counts = cpu_nnz_counts_chunked(
                    cellxgene, cluster_mat,
                    is_sparse_input=sparse.issparse(cellxgene),
                    chunk_size=cpu_chunk_size,
                    verbosity=verbosity
                )

                group_sizes = np.array(cluster_mat.sum(axis=1)).ravel()
                for k in range(n_cluster):
                    n_cells_k = group_sizes[k]
                    threshold = max(
                        n_cells_k * _ept,
                        expressed_min_num_cells_in_target_group
                    )
                    genexlambda[nnz_counts[:, k] < threshold, k] = -1
                del nnz_counts

            _cpu_streaming_used = True
            del cosine_sim, cosine_sim_clean

        # Original CPU path (batch mode legacy or streaming explicitly disabled)
        elif batch_key is not None:
            
            # Validate batch_key exists
            if batch_key not in adata.obs.columns:
                raise ValueError(
                    f"batch_key '{batch_key}' not found in adata.obs. "
                    f"Available columns: {list(adata.obs.columns)}"
                )
            
            
            batch_info = adata.obs[batch_key].values  
            unique_batches = np.unique(batch_info)

            n_gene = cellxgene.shape[1]  # number of genes
            
            if sparse.issparse(cellxgene):
                # Initialize a sparse accumulator and valid count vector
                accum = sparse.csr_matrix((n_gene, n_cluster), dtype=float)
                
                ### This vector is used to keep track of how many batches provide valid
                ### (i.e. above the batch_cell_number_threshold) data for each cluster. 
                valid_counts = np.zeros(n_cluster, dtype=float)

                for batch in unique_batches:
                    indices = np.flatnonzero(batch_info == batch)
                    if indices.size == 0:
                        continue

                    # Compute number of cells per cluster for this batch.
                    cluster_counts = np.array(cluster_mat[:, indices].sum(axis=1)).flatten()
                    # valid_mask: 1.0 for clusters with enough cells, 0.0 otherwise
                    valid_mask = (cluster_counts >= batch_cell_number_threshold).astype(float)
                    if not np.any(valid_mask):
                        continue

                    # Compute cosine similarity for the batch (sparse output)
                    sim_batch = cosine_similarity(
                        X=cellxgene[indices].T,
                        Y=cluster_mat[:, indices],
                        dense_output=False
                    )  # shape: (n_gene, n_cluster)
                    # Multiply each column by valid_mask using a sparse diagonal matrix
                    valid_diag = sparse.diags(valid_mask)
                    sim_batch_valid = sim_batch.dot(valid_diag)

                    accum = accum + sim_batch_valid
                    valid_counts += valid_mask

                # Create a diagonal matrix for column-wise division.
                # For clusters with valid_counts==0, we assign np.nan.
                divisor = np.array([1/v if v > 0 else np.nan for v in valid_counts])
                divisor_diag = sparse.diags(divisor)
                # Multiply accumulator by divisor_diag to scale each column appropriately.
                accum = accum.dot(divisor_diag)
                cosine_sim = accum
            else:
                # Dense branch: use dense arrays.
                accum = np.zeros((n_gene, n_cluster))
                valid_counts = np.zeros(n_cluster)
                for batch in unique_batches:
                    indices = np.flatnonzero(batch_info == batch)
                    if indices.size == 0:
                        continue
                    cluster_counts = np.array(cluster_mat[:, indices].sum(axis=1)).flatten()
                    valid_mask = (cluster_counts >= batch_cell_number_threshold).astype(float)
                    if not np.any(valid_mask):
                        continue
                    sim_batch = cosine_similarity(
                        X=cellxgene[indices].T,
                        Y=cluster_mat[:, indices],
                        dense_output=True
                    )
                    sim_batch = sim_batch * valid_mask  # broadcasting valid_mask along columns
                    accum += sim_batch
                    valid_counts += valid_mask
                for j in range(n_cluster):
                    if valid_counts[j] > 0:
                        accum[:, j] /= valid_counts[j]
                    else:
                        accum[:, j] = np.nan  # or zero if you prefer
                cosine_sim = accum

        else:
            # no batch splitting
            if sparse.issparse(cellxgene):
                ### the dimension is: Gene x lambda
                cosine_sim = cosine_similarity(
                    X=cellxgene.T, 
                    Y=cluster_mat, 
                    dense_output=False
                )
            else:
                ### Convert to dense matrix
                cluster_mat=cluster_mat.toarray()
                ## Not using sparse matrix    
                cosine_sim=cosine_similarity(X=cellxgene.T, Y=cluster_mat, dense_output=True)     
        
        if not _cpu_streaming_used:
            if sparse.issparse(cellxgene):
                ### the dimension is: Gene x lambda

                ### Instead of using cosine_sim.multiply(cosine_sim), the following commands could keep the nonzero values order the same, which would be very useful for the downstream analysis
                ### Because all the calculation would be performed in 1D then
                genexlambda=cosine_sim.copy()
                genexlambda.data=np.multiply(genexlambda.data, genexlambda.data)
                #cosine_sim_data = cosine_sim.data  # Direct access to non-zero elements
                e_power2_sum = np.array(genexlambda.sum(axis=1)).flatten()  # Row-wise sum as a dense array
                if mu==1:
                    ### The np.diff(cosine_sim.indptr) is cool because np.diff(genexlambda.indptr) gives the number of non-zero elements per row
                    ### this avoids generating a large dense matrix and subseting it
                    ### as the .data will list all the nonzero values row by row, so every values are in the same order
                    ### add this basically gives the number of times (in an order, from pos 0 to pos N) to repreat for each element in the array
                    genexlambda.data = genexlambda.data / np.repeat(e_power2_sum, np.diff(genexlambda.indptr))
                else:
                    genexlambda.data=genexlambda.data/((1 - mu) * genexlambda.data + mu * np.repeat(e_power2_sum, np.diff(genexlambda.indptr)))
                ### Because I use genexlambda.data=np.multiply(genexlambda.data, genexlambda.data), so the nonzero values order are the same
                genexlambda.data=np.multiply(genexlambda.data, cosine_sim.data)

            ### If the cellxgene is not a sparse matrix
            else:
                pos_nonzero=cosine_sim!=0
                e_power2=np.multiply(cosine_sim,cosine_sim)
                e_power2_sum=np.sum(e_power2,axis=1)
                e_power2[pos_nonzero]=np.true_divide(e_power2[pos_nonzero],(1-mu)*e_power2[pos_nonzero]+mu*(np.dot(e_power2_sum.reshape(e_power2_sum.shape[0],1),np.repeat(1,e_power2.shape[1]).reshape(1,e_power2.shape[1]))[pos_nonzero]))
                e_power2[pos_nonzero]=np.multiply(e_power2[pos_nonzero],cosine_sim[pos_nonzero])
                genexlambda=e_power2

    ### Refer to scanpy
    rank_stats=None

    # ── Per-group loop (CPU, batched-GPU, and chunked-GPU paths) ─────────
    if True:

        ### Whether to calculate logfoldchanges, because this is required in scanpy 1.8
        anndata_obj = None  # Initialize before conditional creation
        if calculate_logfoldchanges:
            ### Calculate basic stats
            ### Refer to Scanpy
            # for clarity, rename variable
            if groups == 'all':
                groups_order2 = 'all'
            elif isinstance(groups, (str, int)):
                raise ValueError('Specify a sequence of groups')
            else:
                groups_order2 = list(groups)
                if isinstance(groups_order2[0], int):
                    groups_order2 = [str(n) for n in groups_order2]
                if reference != 'rest' and reference not in set(groups_order2):
                    groups_order2 += [reference]
            if reference != 'rest' and reference not in adata.obs[groupby].cat.categories:
                cats = adata.obs[groupby].cat.categories.tolist()
                raise ValueError(
                    f'reference = {reference} needs to be one of groupby = {cats}.'
                )
            pts=False
            anndata_obj = _RankGenes(adata, groups_order2, groupby, reference, use_raw, layer, pts)
            anndata_obj._basic_stats()


        ### Refer to Scanpy
        # for correct getnnz calculation
        ### get non-zeros for columns
        _skip_downstream_filter = (_device == 'gpu') or (_device != 'gpu' and _cpu_streaming_used)
        if not _skip_downstream_filter:
            # Count true nonzeros without mutating the input matrix
            # This avoids eliminate_zeros() which would modify adata.X in-place
            if sparse.issparse(cellxgene):
                get_nonzeros = lambda X: np.asarray((X != 0).sum(axis=0)).ravel()
            else:
                get_nonzeros = lambda X: np.count_nonzero(X, axis=0)


        # Pre-allocate result containers BEFORE the loop
        results_names = {}
        results_scores = {}
        results_logfoldchanges = {} if calculate_logfoldchanges else None

        for order_i, group_i in enumerate(groups_order):
            idx_i=group_info==group_i
            ### Convert to numpy array
            idx_i=idx_i.values

            # REVERT to original column-by-column extraction
            if sparse.issparse(genexlambda):
                scores = genexlambda[:, order_i].toarray()[:, 0]
            else:
                scores = genexlambda[:, order_i]

            if remove_lowly_expressed and not _skip_downstream_filter:
                # GPU and CPU streaming paths already applied expressed_pct filter internally
                n_cells_expressed = get_nonzeros(cellxgene[idx_i])
                n_cells_i = np.sum(idx_i)

                # Warn if cluster is smaller than minimum threshold
                if n_cells_i < expressed_min_num_cells_in_target_group:
                    if verbosity > 0:
                        print(
                            f"Warning: Group '{group_i}' has only {n_cells_i} cells, which is fewer than "
                            f"expressed_min_num_cells_in_target_group={expressed_min_num_cells_in_target_group}. "
                            f"All genes will be filtered out for this group."
                        )

                # Use the maximum of percentage-based and absolute minimum threshold
                threshold = max(
                    n_cells_i * expressed_pct,
                    expressed_min_num_cells_in_target_group
                )
                scores[n_cells_expressed < threshold] = -1

            global_indices = _select_top_n(scores, n_genes_user)

            # Store results in dictionaries instead of concatenating DataFrames
            results_names[group_i] = adata.var_names.values[global_indices]
            results_scores[group_i] = scores[global_indices]

            if calculate_logfoldchanges and anndata_obj is not None:
                group_index = np.where(anndata_obj.groups_order == group_i)[0][0]
                if anndata_obj.means is not None:
                    mean_group = anndata_obj.means[group_index]
                    mean_rest = (
                        anndata_obj.means_rest[group_index]
                        if anndata_obj.ireference is None
                        else anndata_obj.means[anndata_obj.ireference]
                    )
                    foldchanges = (anndata_obj.expm1_func(mean_group) + 1e-9) / (
                        anndata_obj.expm1_func(mean_rest) + 1e-9
                    )
                    results_logfoldchanges[group_i] = np.log2(foldchanges[global_indices])

        # Build DataFrame once after the loop
        columns_data = {}
        for group_i in groups_order:
            columns_data[(group_i, 'names')] = results_names[group_i]
            columns_data[(group_i, 'scores')] = results_scores[group_i]
            if calculate_logfoldchanges and group_i in results_logfoldchanges:
                columns_data[(group_i, 'logfoldchanges')] = results_logfoldchanges[group_i]
        rank_stats = pd.DataFrame(columns_data)

    #### also return a copy of the results showing the results by group
    if return_by_group:
        # Need to swap levels first to ensure attribute comes first, then cell group
        rank_stats_swapped = rank_stats.copy()
        rank_stats_swapped.columns = rank_stats_swapped.columns.swaplevel()
        
        # Create a new DataFrame with flattened column names using the specified delimiter
        # Format: attribute::cell_group
        flattened_columns = [column_delimiter.join(map(str, col)) for col in rank_stats_swapped.columns]
        cosg_results = rank_stats_swapped.copy()
        cosg_results.columns = flattened_columns
        
        # Store with flattened column names for compatibility with adata.write
        adata.uns[key_added]['COSG'] = cosg_results
    
    
    ## Refer to scanpy
    dtypes = {
            'names': 'O',
            'scores': 'float32',
            'logfoldchanges': 'float32',
    }
    
    rank_stats.columns = rank_stats.columns.swaplevel()
    for col in rank_stats.columns.levels[0]:
        adata.uns[key_added][col]=rank_stats[col].to_records(
            index=False, column_dtypes=dtypes[col]
        )
        
    if verbosity>0:    
        print(f"Finished identifying marker genes by COSG, and the results are in adata.uns['{key_added}'].")
        
    ### Return the result
    return adata if copy else None



def indexByGene(
    df: pd.DataFrame,
    gene_key: str = "names",
    score_key: str = "scores",
    column_delimiter: str = "::",
    set_nan_to_zero: bool= False,
    convert_negative_one_to_zero: bool = True
):
    """
    Reshapes a DataFrame with flattened column names (converted from MultiIndex columns) where 
    gene names are under the specified key and scores are under the corresponding key.
    The resulting DataFrame will have gene names as the index and scores for each cell group in separate columns.

    Note
    ----
    This function is designed for reindexing COSG's output stored in `adata.uns['cosg']['COSG']`.
    It works with the flattened column format where attribute and cell group are joined by a delimiter
    in the format 'attribute{delimiter}cell_group'.
    It is recommended to set `n_genes_user=adata.n_vars` when calling the `cosg.cosg` function to ensure 
    that scores for all genes are returned.
    
    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with MultiIndex columns.
    gene_key : str, optional, default="names"
        The key used for gene names in the first level of MultiIndex.
    score_key : str, optional, default="scores"
        The key used for scores in the first level of MultiIndex.
    column_delimiter : str, optional, default="::"
        Delimiter used to separate attribute and cell group in the column names.
    set_nan_to_zero : bool, optional, default=False
        If True, replaces all NaN values in the resulting DataFrame with 0.
        Set this parameter to True if `n_genes_user=adata.n_vars` is not set when calling the `cosg.cosg` function.
    convert_negative_one_to_zero : bool, optional, default=True
        If True, replaces all occurrences of -1.0 in the resulting DataFrame with 0.


    Returns
    -------
    pd.DataFrame
        Reshaped DataFrame with genes as index and scores per cell type.
        

    Example
    -------
    >>> cosg_df = indexByGene(
    ...     adata.uns['cosg']['COSG'],
    ...     gene_key="names",
    ...     score_key="scores",
    ...     column_delimiter="::",
    ...     set_nan_to_zero=True,
    ...     convert_negative_one_to_zero=True
    ... )
    >>>  ### Check the reindexed data frame
    >>> cosg_df.head()
    """
    
        # Validate input type
    if not isinstance(df, pd.DataFrame):
        raise TypeError("The input df must be a pandas DataFrame.")
    
    # Extract cell groups from the column names while preserving original order
    cell_groups = []
    for col in df.columns:
        parts = col.split(column_delimiter)
        if len(parts) > 2:
            raise ValueError(
                f"Column name '{col}' contains multiple instances of the delimiter "
                f"'{column_delimiter}'. This may indicate that gene names or group names "
                f"contain the delimiter. Please use a different column_delimiter."
            )
        if len(parts) == 2 and parts[0] == gene_key:  # Only look at gene columns to maintain order
            if parts[1] not in cell_groups:
                cell_groups.append(parts[1])
    
    # Gather all genes from all 'names' columns
    all_genes = []
    for cell_group in cell_groups:
        gene_col = f"{gene_key}{column_delimiter}{cell_group}"
        if gene_col in df.columns:
            all_genes.extend(df[gene_col].tolist())
    all_genes = pd.Series(all_genes).unique()
    
    # Create results DataFrame
    df_scores = pd.DataFrame(index=all_genes)
    
    # Fill in scores for each cell group
    for cell_group in cell_groups:
        gene_col = f"{gene_key}{column_delimiter}{cell_group}"
        score_col = f"{score_key}{column_delimiter}{cell_group}"
        
        if gene_col in df.columns and score_col in df.columns:
            # Create a temporary gene-to-score mapping
            gene_to_score = dict(zip(df[gene_col], df[score_col]))
            # Populate the scores column
            df_scores[cell_group] = df_scores.index.map(lambda x: gene_to_score.get(x, float('nan')))
    
    # Apply optional transformations
    if set_nan_to_zero:
        df_scores.fillna(0, inplace=True)
    
    # Optionally replace -1.0 values with 0
    if convert_negative_one_to_zero:
        df_scores = df_scores.replace(-1.0, 0)
    
    return df_scores




def iqrLogNormalize(
    df: pd.DataFrame,
    q_upper: float = 0.95,
    q_lower: float = 0.75
):
    """
    Applies an IQR-based scaling to a DataFrame followed by a log1p transformation.
    
    The function computes the interquartile range (IQR) for each column as the difference 
    between the q_upper and q_lower quantiles. Each column is then divided by its respective IQR.
    
    If any IQR value is zero, it is replaced with the smallest nonzero IQR found among all columns.
    If all IQR values are zero, a replacement value of 1e-6 is used.
    
    Optionally, the log1p-transformed result can be further scaled using min–max normalization 
    to rescale the values to the [0, 1] range (applied per column).
    
    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame containing numerical values.
    q_upper : float, optional, default=0.95
        The upper quantile used for IQR calculation.
    q_lower : float, optional, default=0.75
        The lower quantile used for IQR calculation.
      
    Returns
    -------
    pd.DataFrame
        A DataFrame with the same shape as the input containing the IQR-scaled and log1p-transformed values.
    
    Raises
    ------
    TypeError
        If the input df is not a pandas DataFrame.
    ValueError
        If q_lower is not in the range [0, q_upper) or if q_upper is not in (q_lower, 1].
    
    Example
    -------
    >>> ## cosg_df is a DataFrame with gene scores, output from the cosg.indexByGene function
    >>> cosg_df_transformed = iqr_log_normalize(cosg_df, q_upper=0.95, q_lower=0.75)
    """
    # Validate input type
    if not isinstance(df, pd.DataFrame):
        raise TypeError("The input df must be a pandas DataFrame.")
    
    # Validate quantile parameters
    if not (0 <= q_lower < q_upper <= 1):
        raise ValueError("q_lower must be >= 0 and less than q_upper, and q_upper must be <= 1.")
    
    # Compute the IQR for each column
    iqr = df.quantile(q_upper) - df.quantile(q_lower)
    
    # Replace zero IQR values with the smallest nonzero IQR (or 1e-6 if all are zero)
    if (iqr == 0).any():
        nonzero_iqr = iqr[iqr > 0]
        replacement_value = nonzero_iqr.min() if not nonzero_iqr.empty else 1e-6
        iqr = iqr.replace(0, replacement_value)
    
    # Scale the DataFrame by dividing each column by its corresponding IQR
    df_iqr_scaled = df / iqr
    
    # Apply the log1p transformation
    df_log_transformed = np.log1p(df_iqr_scaled)

    return df_log_transformed




#############
########
####
###



### Import from Scanpy
def select_groups(adata, groups_order_subset='all', key='groups'):
    """Get subset of groups in adata.obs[key]."""
    groups_order = adata.obs[key].cat.categories
    if key + '_masks' in adata.uns:
        groups_masks = adata.uns[key + '_masks']
    else:
        groups_masks = np.zeros(
            (len(adata.obs[key].cat.categories), adata.obs[key].values.size), dtype=bool
        )
        for iname, name in enumerate(adata.obs[key].cat.categories):
            # if the name is not found, fallback to index retrieval
            if adata.obs[key].cat.categories[iname] in adata.obs[key].values:
                mask = adata.obs[key].cat.categories[iname] == adata.obs[key].values
            else:
                mask = str(iname) == adata.obs[key].values
            groups_masks[iname] = mask
    groups_ids = list(range(len(groups_order)))
    if groups_order_subset != 'all':
        groups_ids = []
        for name in groups_order_subset:
            groups_ids.append(
                np.where(adata.obs[key].cat.categories.values == name)[0][0]
            )
        if len(groups_ids) == 0:
            # fallback to index retrieval
            groups_ids = np.where(
                np.in1d(
                    np.arange(len(adata.obs[key].cat.categories)).astype(str),
                    np.array(groups_order_subset),
                )
            )[0]
        if len(groups_ids) == 0:
            raise ValueError(
                f'{np.array(groups_order_subset)} invalid! Specify valid '
                f'groups_order (or indices) from {adata.obs[key].cat.categories}'
            )
#             from sys import exit

#             exit(0)
        groups_masks = groups_masks[groups_ids]
        groups_order_subset = adata.obs[key].cat.categories[groups_ids].values
    else:
        groups_order_subset = groups_order.values
    return groups_order_subset, groups_masks


### Import from Scanpy   
import numba




def _get_mean_var(X, *, axis=0):
    if sparse.issparse(X):
        mean, var = sparse_mean_variance_axis(X, axis=axis)
    else:
        mean = np.mean(X, axis=axis, dtype=np.float64)
        mean_sq = np.multiply(X, X).mean(axis=axis, dtype=np.float64)
        var = mean_sq - mean ** 2
    # enforce R convention (unbiased estimator) for variance
    var *= X.shape[axis] / (X.shape[axis] - 1)
    return mean, var


def sparse_mean_variance_axis(mtx: sparse.spmatrix, axis: int):
    """
    This code and internal functions are based on sklearns
    `sparsefuncs.mean_variance_axis`.
    Modifications:
    * allow deciding on the output type, which can increase accuracy when calculating the mean and variance of 32bit floats.
    * This doesn't currently implement support for null values, but could.
    * Uses numba not cython
    """
    assert axis in (0, 1)
    if isinstance(mtx, sparse.csr_matrix):
        ax_minor = 1
        shape = mtx.shape
    elif isinstance(mtx, sparse.csc_matrix):
        ax_minor = 0
        shape = mtx.shape[::-1]
    else:
        raise ValueError("This function only works on sparse csr and csc matrices")
    if axis == ax_minor:
        return sparse_mean_var_major_axis(
            mtx.data, mtx.indices, mtx.indptr, *shape, np.float64
        )
    else:
        return sparse_mean_var_minor_axis(mtx.data, mtx.indices, *shape, np.float64)


@numba.njit(cache=True)
def sparse_mean_var_minor_axis(data, indices, major_len, minor_len, dtype):
    """
    Computes mean and variance for a sparse matrix for the minor axis.
    Given arrays for a csr matrix, returns the means and variances for each
    column back.
    """
    non_zero = indices.shape[0]

    means = np.zeros(minor_len, dtype=dtype)
    variances = np.zeros_like(means, dtype=dtype)

    counts = np.zeros(minor_len, dtype=np.int64)

    for i in range(non_zero):
        col_ind = indices[i]
        means[col_ind] += data[i]

    for i in range(minor_len):
        means[i] /= major_len

    for i in range(non_zero):
        col_ind = indices[i]
        diff = data[i] - means[col_ind]
        variances[col_ind] += diff * diff
        counts[col_ind] += 1

    for i in range(minor_len):
        variances[i] += (major_len - counts[i]) * means[i] ** 2
        variances[i] /= major_len

    return means, variances


@numba.njit(cache=True)
def sparse_mean_var_major_axis(data, indices, indptr, major_len, minor_len, dtype):
    """
    Computes mean and variance for a sparse array for the major axis.
    Given arrays for a csr matrix, returns the means and variances for each
    row back.
    """
    means = np.zeros(major_len, dtype=dtype)
    variances = np.zeros_like(means, dtype=dtype)

    for i in range(major_len):
        startptr = indptr[i]
        endptr = indptr[i + 1]
        counts = endptr - startptr

        for j in range(startptr, endptr):
            means[i] += data[j]
        means[i] /= minor_len

        for j in range(startptr, endptr):
            diff = data[j] - means[i]
            variances[i] += diff * diff

        variances[i] += (minor_len - counts) * means[i] ** 2
        variances[i] /= minor_len

    return means, variances  

### Import from Scanpy
class _RankGenes:
    def __init__(
        self,
        adata,
        groups,
        groupby,
        reference='rest',
        use_raw=True,
        layer=None,
        comp_pts=False,
    ):

        if 'log1p' in adata.uns_keys() and adata.uns['log1p']['base'] is not None:
            self.expm1_func = lambda x: np.expm1(x * np.log(adata.uns['log1p']['base']))
        else:
            self.expm1_func = np.expm1

        self.groups_order, self.groups_masks = select_groups(
            adata, groups, groupby
        )

        # Singlet groups cause division by zero errors
        invalid_groups_selected = set(self.groups_order) & set(
            adata.obs[groupby].value_counts().loc[lambda x: x < 2].index
        )

        if len(invalid_groups_selected) > 0:
            raise ValueError(
                "Could not calculate statistics for groups {} since they only "
                "contain one sample.".format(', '.join(invalid_groups_selected))
            )

        adata_comp = adata
        if layer is not None:
            if use_raw:
                raise ValueError("Cannot specify `layer` and have `use_raw=True`.")
            X = adata_comp.layers[layer]
        else:
            if use_raw and adata.raw is not None:
                adata_comp = adata.raw
            X = adata_comp.X

        # Store reference without mutating the original matrix
        self.X = X
        self.var_names = adata_comp.var_names

        self.ireference = None
        if reference != 'rest':
            self.ireference = np.where(self.groups_order == reference)[0][0]

        self.means = None
        self.vars = None

        self.means_rest = None
        self.vars_rest = None

        self.comp_pts = comp_pts
        self.pts = None
        self.pts_rest = None

        self.stats = None

        # for logreg only
        self.grouping_mask = adata.obs[groupby].isin(self.groups_order)
        self.grouping = adata.obs.loc[self.grouping_mask, groupby]

    def _basic_stats(self):
        n_genes = self.X.shape[1]
        n_groups = self.groups_masks.shape[0]

        self.means = np.zeros((n_groups, n_genes))
        self.vars = np.zeros((n_groups, n_genes))
        self.pts = np.zeros((n_groups, n_genes)) if self.comp_pts else None

        if self.ireference is None:
            self.means_rest = np.zeros((n_groups, n_genes))
            self.vars_rest = np.zeros((n_groups, n_genes))
            self.pts_rest = np.zeros((n_groups, n_genes)) if self.comp_pts else None
        else:
            mask_rest = self.groups_masks[self.ireference]
            X_rest = self.X[mask_rest]
            self.means[self.ireference], self.vars[self.ireference] = _get_mean_var(
                X_rest
            )
            # deleting the next line causes a memory leak for some reason
            del X_rest

        if sparse.issparse(self.X):
            get_nonzeros = lambda X: np.asarray((X != 0).sum(axis=0)).ravel()
        else:
            get_nonzeros = lambda X: np.count_nonzero(X, axis=0)

        for imask, mask in enumerate(self.groups_masks):
            X_mask = self.X[mask]

            if self.comp_pts:
                self.pts[imask] = get_nonzeros(X_mask) / X_mask.shape[0]

            if self.ireference is not None and imask == self.ireference:
                continue

            self.means[imask], self.vars[imask] = _get_mean_var(X_mask)

            if self.ireference is None:
                mask_rest = ~mask
                X_rest = self.X[mask_rest]
                self.means_rest[imask], self.vars_rest[imask] = _get_mean_var(X_rest)
                # this can be costly for sparse data
                if self.comp_pts:
                    self.pts_rest[imask] = get_nonzeros(X_rest) / X_rest.shape[0]
                # deleting the next line causes a memory leak for some reason
                del X_rest

