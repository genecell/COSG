from .cosg import indexByGene, iqrLogNormalize

# scanpy is NOT imported at module level. Only `plotMarkerDotplot` needs it
# (it wraps `sc.pl.dotplot`), so it is imported lazily there and declared as
# the optional `cosg[dotplot]` extra. Importing it here would make scanpy —
# and statsmodels, seaborn, umap-learn behind it — mandatory for everyone who
# writes `import cosg`, including users who only compute markers.
from anndata import AnnData
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import matplotlib.colors as mcolors
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram as _sch_dendrogram
from scipy.spatial.distance import pdist, squareform
import matplotlib.patheffects as PathEffects
from scipy.sparse import issparse


def _require_scanpy(feature: str):
    """Import scanpy on demand, with an error that says what to install."""
    try:
        import scanpy as sc
    except ImportError as exc:  # pragma: no cover - exercised in a bare venv
        raise ImportError(
            f"{feature} needs scanpy, which is an optional dependency of COSG.\n"
            "    pip install 'cosg[dotplot]'   (or: pip install scanpy)\n"
            "Marker detection and the other plotting functions do not need it."
        ) from exc
    return sc


def _dendrogram_order(
    adata,
    groupby,
    use_rep,
    *,
    cor_method="pearson",
    linkage_method="complete",
    key_added=None,
    inplace=True,
):
    """Order `groupby` categories by hierarchical clustering of `use_rep`.

    A pandas/scipy reimplementation of ``scanpy.tl.dendrogram``, so that
    ordering does not drag scanpy into the import path. Same algorithm:
    per-category means of the representation, correlation between categories,
    ``1 - corr`` as the distance, then the leaf order of the linkage.

    The result is written to ``adata.uns[key_added]`` using scanpy's own key
    names, so anything that *reads* that entry — including
    ``sc.pl.dotplot(..., dendrogram=True)`` — keeps working unchanged.
    """
    if use_rep not in adata.obsm:
        raise KeyError(
            f"use_rep='{use_rep}' not found in adata.obsm "
            f"(available: {list(adata.obsm)})"
        )

    categorical = adata.obs[groupby]
    if not hasattr(categorical, "cat"):
        categorical = categorical.astype("category")
    categories = list(categorical.cat.categories)
    if len(categories) < 2:
        raise ValueError(
            f"groupby='{groupby}' has {len(categories)} category; a dendrogram "
            "needs at least 2."
        )

    rep_df = pd.DataFrame(np.asarray(adata.obsm[use_rep]), index=adata.obs_names)
    mean_df = rep_df.groupby(categorical.values, observed=True).mean()
    mean_df = mean_df.reindex(categories)

    corr_matrix = mean_df.T.corr(method=cor_method).clip(-1, 1)
    # squareform needs an exactly-zero diagonal; clip leaves 1.0 on it, but
    # floating point can leave 1 - corr at ~1e-16 instead of 0.
    dist = 1 - corr_matrix.values
    np.fill_diagonal(dist, 0.0)
    z_var = linkage(squareform(dist, checks=False), method=linkage_method)
    dendro_info = _sch_dendrogram(z_var, labels=categories, no_plot=True)
    categories_ordered = dendro_info["ivl"]

    result = dict(
        linkage=z_var,
        groupby=[groupby],
        use_rep=use_rep,
        cor_method=cor_method,
        linkage_method=linkage_method,
        categories_ordered=categories_ordered,
        categories_idx_ordered=[categories.index(c) for c in categories_ordered],
        dendrogram_info=dendro_info,
        correlation_matrix=corr_matrix.values,
    )
    if inplace:
        adata.uns[key_added or f"dendrogram_{groupby}"] = result
    return result


## Helper function
def _compute_gene_expression_percentage(adata, group_by, cosg_score_df, layer=None):
    """
    Computes the percentage of cells expressing genes in `cosg_score_df` within each cell type group.

    This optimized function calculates expression only for the genes present in `cosg_score_df`,
    making it significantly more efficient than computing for all genes.

    Parameters
    ----------
    adata : AnnData
        The AnnData object containing expression data.
    group_by : str
        The observation column in `adata.obs` to group cells.
    cosg_score_df : pd.DataFrame
        A DataFrame containing COSG scores, where index corresponds to marker genes.
    layer : str, optional (default: None)
        If provided, uses `adata.layers[layer]` for expression data; otherwise, uses `adata.X`.

    Returns
    -------
    pd.DataFrame
        A DataFrame where:
        - Rows correspond to genes.
        - Columns correspond to cell types.
        - Values represent the percentage of cells expressing the gene in that cell type.
    """

    # Get only the relevant genes that exist in adata
    genes_to_use = cosg_score_df.index.intersection(adata.var_names)
    if len(genes_to_use) == 0:
        raise ValueError("No valid genes from cosg_score_df are found in adata.var_names.")

    # Extract the relevant expression data
    expr_data = adata[:, genes_to_use].X if layer is None else adata[:, genes_to_use].layers[layer]

    # Convert to binary presence/absence (1 if expressed, 0 otherwise)
    if issparse(expr_data):
        expr_data = expr_data.copy()  # Ensure no modification of original data
        expr_data.data[:] = 1  # Convert all nonzero values to 1
    else:
        expr_data = (expr_data > 0).astype(int)

    # Convert to DataFrame (cells as rows, genes as columns)
    expr_df = pd.DataFrame(expr_data.toarray() if issparse(expr_data) else expr_data,
                           index=adata.obs_names, columns=genes_to_use)

    # Compute the sum of expressing cells for each group (vectorized)
    expr_sums = expr_df.groupby(adata.obs[group_by], observed=True).sum()  # ✅ Fix applied: observed=True

    # Compute percentage of expressing cells per group
    group_sizes = adata.obs[group_by].value_counts().reindex(expr_sums.index, fill_value=0).values
    expr_percentages = (expr_sums.div(group_sizes, axis=0) * 100).T  # Transpose so genes are rows

    return expr_percentages





# Radial layout helper function
def _build_subtree_sizes(G, node, subtree_size, visited):
    """
    Recursively compute the number of leaf descendants for each node.
    """
    if node in visited:
        return subtree_size[node]
    visited.add(node)
    children = list(G.successors(node))
    if len(children) == 0:
        subtree_size[node] = 1
    else:
        total = 0
        for c in children:
            total += _build_subtree_sizes(G, c, subtree_size, visited)
        subtree_size[node] = total
    return subtree_size[node]

def _radial_dendrogram_layout(G, root, radius_step=1.5, start_angle=0, end_angle=2*np.pi):
    """
    Compute a radial layout for a tree (with the root at the center).
    Angles are distributed in proportion to the number of leaf nodes.
    """
    subtree_size = {}
    _build_subtree_sizes(G, root, subtree_size, visited=set())
    pos = {}

    def recurse(node, r, alpha_start, alpha_end):
        alpha_mid = 0.5 * (alpha_start + alpha_end)
        pos[node] = (r * np.cos(alpha_mid), r * np.sin(alpha_mid))
        children = list(G.successors(node))
        if len(children) == 0:
            return
        total_leaves = sum(subtree_size[ch] for ch in children)
        r_next = r + radius_step
        angle_offset = alpha_start
        for ch in children:
            frac = subtree_size[ch] / total_leaves
            ch_alpha_start = angle_offset
            ch_alpha_end = angle_offset + frac * (alpha_end - alpha_start)
            recurse(ch, r_next, ch_alpha_start, ch_alpha_end)
            angle_offset = ch_alpha_end

    recurse(root, 0.0, start_angle, end_angle)
    return pos






### Plot marker specificity with cell type dendrogram information

# Import matplotlib modules needed for curved edges
import matplotlib.path as mpath
import matplotlib.patches as mpatches

def plotMarkerDendrogram(
    adata: AnnData,
    group_by: str,
    use_rep: str = 'X_pca',
    calculate_dendrogram_on_cosg_scores: bool = True,
    top_n_genes: int = 3,
    cosg_key: str = 'cosg',
    radius_step: float = 1.5,
    cmap: str = "Purples",
    cell_type_label_offset: float = 0,
    gene_label_offset: float = 0.25,
    gene_label_color: str = None,
    linkage_method: str = "ward",
    distance_metric: str = "euclidean",
    hierarchy_merge_scale: float = None,
    collapse_scale: float = None,
    add_cluster_node_for_single_node_cluster : bool = True,
    palette=None,
    gene_color_min: float = 0,
    gene_color_max: float = None,
    font_outline: float = 2,
    figure_size: tuple = (10, 10),
    node_shape_cell_type: str = 'o',
    node_shape_gene: str = 's',
    node_shape_internal: str = 'o',
    colorbar_width: float = 0.01,
    layer: str = None,
    gene_size_scale: float = 300,
    map_cell_type_gene: dict = None,
    cell_type_selected: list = None,
    color_root_node: str = '#D6EFD5',
    color_internal_node: str = 'lightgray',
    color_edge: str = 'lightgray',
    edge_curved: float = 0.0,
    show_figure: bool = True,
    save: str = None,
):
    """
    Visualizes a radial dendrogram of cell types with attached top marker genes.
    
    Computes a dendrogram in two modes:
    - If `calculate_dendrogram_on_cosg_scores` is True, uses COSG scores from `adata.uns['cosg']['COSG']`, 
      processed with `indexByGene()` and `iqrLogNormalize()`, then computes the dendrogram on the transposed DataFrame 
      with the specified `distance_metric` and `linkage_method`.
    - If False, aggregates `adata.obsm[use_rep]` by `adata.obs[group_by]`, computing distances with `pdist` 
      and linkage with the given `distance_metric` and `linkage_method`.

    When `collapse_scale` (0 to 1) is set and yields multiple clusters, cell types are grouped with cluster nodes; 
    if only one cluster, they attach directly to the root unless `add_cluster_node_for_single_node_cluster` is True, 
    which adds a cluster node for single-member clusters. If `collapse_scale` is None, `hierarchy_merge_scale` 
    (0 to 1) controls merging binary nodes into multi-child nodes based on distance similarity, with no merging 
    if None. Top marker genes (from COSG data) are added as nodes to cell type leaves, with labels offset by 
    `gene_label_offset` and colored by `gene_label_color` if provided.

    Cell type node colors come from `palette`:
    - Dictionary: Maps cell types to colors.
    - List: Assigns colors by cell type order.
    - None: Uses `adata.uns[f"{group_by}_colors"]` if available, else defaults to "lightblue".
    Marker gene node colors are scaled between `gene_color_min` and `gene_color_max` (max defaults to the 
    highest score). Node sizes reflect expression percentage (fraction of cells with expression > 0) from 
    `adata.X` or `adata.layers[layer]`.
    
    
    Parameters
    ----------
    adata : AnnData
        An AnnData object.
    group_by : str
        The observation key in adata.obs to group cell types.
    use_rep : str, optional, default='X_pca'
        The representation to use when aggregating data from adata.obsm.
    calculate_dendrogram_on_cosg_scores : bool, optional, default=True
        If True, compute the dendrogram on COSG scores derived using cosg.cosg, cosg.indexByGene and cosg.iqrLogNormalize.
        If False, compute the dendrogram on the aggregated representation from adata.obsm[use_rep].
    top_n_genes : int, optional, default=3
        Number of top marker genes (per cell type) to attach.
    cosg_key : str, optional, default='cosg'
        The key used to access the COSG marker gene identification results. Defaults to "cosg".
    radius_step : float, optional, default=1.5
        Radial distance between successive levels in the layout.
    cmap : str, optional, default="Purples"
        The matplotlib colormap to use for gene nodes.
    cell_type_label_offset : float, optional, default=0
        Fractional radial offset for cell type labels from the cell type node.
    gene_label_offset : float, optional, default=0.25
        Fractional radial offset for gene labels from the marker node.
    gene_label_color : str, optional, default=None
        If provided, this color is used for gene labels; otherwise, the gene node's colormap color is used.
    linkage_method : str, optional, default="ward"
        Linkage method to use when computing the dendrogram.
    distance_metric : str, optional, default="euclidean"
        Distance metric to use when computing the dendrogram.
    hierarchy_merge_scale : float or None, optional, default=None
        Controls the merging of binary nodes into multi-child nodes to simulate a non-binary hierarchy when
        collapse_scale is None. If provided, must be a float between 0 and 1, scaling the threshold relative to
        the range of linkage distances in Z. Nodes with distance differences below this scaled threshold are
        merged with their parent, allowing nodes to have more than two children.
        - 0: No merging (retains binary structure).
        - 1: Maximal merging (merges nodes if their distances differ by less than the full distance range).
        If None, no merging is performed, preserving the default binary dendrogram structure from Z.
        Raises ValueError if not between 0 and 1 when provided.
    collapse_scale : float or None, optional, default=None
        Controls the level of clustering in the dendrogram. If None, builds a full hierarchical dendrogram where
        nodes may have more than two children based on distance similarity. If a float between 0 and 1, scales the
        threshold relative to the min and max linkage distances in Z, collapsing leaves and internal nodes with
        distances below this scaled threshold into cluster nodes. 
        - 0: Maximal clustering (collapses at the minimum distance).
        - 1: Minimal clustering (collapses at the maximum distance).
        If only one cluster is found, no extra cluster node is added between the root and leaves. 
        Raises ValueError if not between 0 and 1 when provided.
    add_cluster_node_for_single_node_cluster : bool, optional, default=True
        Determines whether to create a cluster node for clusters containing only a single cell type when
        collapse_scale is provided. If True, a cluster node is added between the root and the single cell type
        node, maintaining a consistent hierarchy. If False, the single cell type node is connected directly to
        the root without an intermediate cluster node. Only applies when collapse_scale is not None and clustering
        results in single-member clusters.
    palette : dict, list, or None, optional, default=None
        Colors for cell type nodes. If a dict, keys are cell type names and values are colors.
        If a list, colors are assigned in order of cell types.
        If None and if adata.uns contains f"{group_by}_colors", that palette is used.
        Otherwise, cell type nodes default to "lightblue".
    gene_color_min : float, optional, default=0
        Minimum value for normalizing marker gene node colors.
    gene_color_max : float or None, optional, default=None
        Maximum value for normalizing marker gene node colors. If None, the maximum among marker scores is used.
    font_outline : float, optional, default=2
        Outline width for text labels.
    figure_size : tuple, optional, default=(10, 10)
        Size of the figure.
    node_shape_cell_type : str, optional, default='d'
        Shape of the cell type nodes. Default is 'd' (diamond). Can be any valid NetworkX node shape.
        Specification is as matplotlib.scatter marker, one of 'so^>v<dph8'. In detail:
        - 'o' : Circle
        - 's' : Square
        - 'd' : Diamond
        - 'v' : Triangle Down
        - '^' : Triangle Up
        - '<' : Triangle Left
        - '>' : Triangle Right
        - 'p' : Pentagon
        - 'h' : Hexagon
        - '8' : Octagon
    node_shape_gene : str, optional, default='o'
        Shape of marker gene nodes. Default is 'o' (circle).
    node_shape_internal : str, optional, default='o'
        Shape of internal dendrogram nodes. Default is 'o' (circle).
    colorbar_width : float, optional, default=0.01
        Width (in normalized figure coordinates) for the colorbar.
    layer : str, optional, default=None
        If provided, use adata.layers[layer] to calculate expression; otherwise, use adata.X.
    gene_size_scale : float, optional, default=300
        Base size for marker gene nodes; final size = gene_size_scale * (expression_percentage / 100).
    map_cell_type_gene : dict, optional, default=None
        Custom mapping of cell types to marker genes. If provided, this will be used instead of the top marker genes.
        Should be a dictionary where keys are cell type names and values are lists of gene names.
        Only genes present in adata.var_names will be included. It's okay if some cell types are not in the dict.
    cell_type_selected : list, optional, default=None
        List of cell types to include in the visualization. If provided, only these cell types will be shown.
        If None, all cell types will be included. Raises ValueError if none of the provided cell types are valid.
    color_root_node : str, optional, default='#D6EFD5'
        Color for the root node. Default is a dark gray (#404040).
    color_internal_node : str, optional, default='lightgray'
        Color for internal nodes and cluster nodes. Default is lightgray.
    color_edge : str, optional, default='lightgray'
        Color for all edges in the dendrogram. Default is lightgray.
    edge_curved : float, optional, default=0.0
        Controls the curvature of edges. 0.0 means straight lines, positive values increase curvature.
        Recommended range: 0.0 to 0.3. Default is 0.0 (straight lines).
    show_figure : bool, optional (default=True)
        Whether to display the figure after plotting.
    save : str or None, optional (default=None)
        File path to save the resulting figure. If None, the figure will not be saved.
    
    Returns
    -------
    None
        Displays a matplotlib figure of the radial dendrogram if `show_figure=True`.
    
    Example
    -------
    >>> import cosg
    >>> cosg.plotMarkerDendrogram(
    ...     adata,
    ...     group_by="CellTypes",
    ...     use_rep="X_pca",
    ...     calculate_dendrogram_on_cosg_scores=False,
    ...     top_n_genes=3,
    ...     radius_step=1.5,
    ...     cmap="Purples",
    ...     gene_label_offset=0.25,
    ...     gene_label_color="black",
    ...     linkage_method="ward",
    ...     distance_metric="correlation",
    ...     collapse_threshold=0.3,
    ...     palette=None,
    ...     gene_color_min=0,
    ...     gene_color_max=None,
    ...     font_outline=2,
    ...     figure_size=(10,10),
    ...     colorbar_width=0.02,
    ...     layer=None,
    ...     gene_size_scale=300
    ... )
    """
    # Compute the transformed COSG scores
    cosg_df = indexByGene(
        adata.uns[cosg_key]['COSG'],
        set_nan_to_zero=True,
        convert_negative_one_to_zero=True
    )
    cosg_score_df = iqrLogNormalize(cosg_df)
    
    # Decide which dendrogram to use
    if calculate_dendrogram_on_cosg_scores:
        data = cosg_score_df.T.values  # rows: cell types, columns: genes
        D = pdist(data, metric=distance_metric)
        Z = linkage(D, method=linkage_method)
        all_cell_types = list(cosg_score_df.columns)
    else:
        rep = adata.obsm[use_rep]
        df_rep = pd.DataFrame(rep, index=adata.obs_names)
        df_rep[group_by] = adata.obs[group_by].values
        group_means = df_rep.groupby(group_by, observed=True).mean()
        all_cell_types = list(group_means.index)
        data = group_means.values
        D = pdist(data, metric=distance_metric)
        Z = linkage(D, method=linkage_method)
    
    # Filter cell types if cell_type_selected is provided
    if cell_type_selected is not None:
        # Check which selected cell types are valid
        valid_selected_cell_types = [ct for ct in cell_type_selected if ct in all_cell_types]
        
        # If no valid cell types, raise an error
        if not valid_selected_cell_types:
            raise ValueError(f"None of the provided cell types {cell_type_selected} are valid. Valid cell types are: {all_cell_types}")
        
        # If only one valid cell type, we can't perform hierarchical clustering
        if len(valid_selected_cell_types) == 1:
            print(f"Only one valid cell type selected ({valid_selected_cell_types[0]}). Skipping hierarchical clustering.")
            # Set up a simplified tree with just one cell type
            G = nx.DiGraph()
            root = "root"
            G.add_node(root, node_type='root')
            
            ct = valid_selected_cell_types[0]
            G.add_node(ct, node_type='cell_type')
            G.add_edge(root, ct)
            
            # Update cell_types
            cell_types = valid_selected_cell_types
        else:
            # Create a mask for the distance matrix and linkage
            cell_type_indices = [all_cell_types.index(ct) for ct in valid_selected_cell_types]
            
            # Filter the data matrix
            if calculate_dendrogram_on_cosg_scores:
                filtered_data = cosg_score_df.iloc[:, [all_cell_types.index(ct) for ct in valid_selected_cell_types]].T.values
            else:
                filtered_data = group_means.iloc[[all_cell_types.index(ct) for ct in valid_selected_cell_types]].values
            
            # Recompute distances and linkage with filtered data
            D = pdist(filtered_data, metric=distance_metric)
            Z = linkage(D, method=linkage_method)
            
            # Update cell_types to only include selected ones
            cell_types = valid_selected_cell_types
    else:
        # Use all cell types
        cell_types = all_cell_types
    
    N = len(cell_types)
    
    # Check if we have at least two cell types for hierarchical clustering
    if N < 2:
        # We already handled the special case with only one selected cell type above,
        # but we'll keep this as a safety check for other code paths
        if 'G' not in locals():  # Only create the graph if not already created
            G = nx.DiGraph()
            root = "root"
            G.add_node(root, node_type='root')
            
            ct = cell_types[0]
            G.add_node(ct, node_type='cell_type')
            G.add_edge(root, ct)
    else:
        ### Build the tree graph
        G = nx.DiGraph()

        # Validate collapse_scale if provided
        if collapse_scale is not None:
            if not (0 <= collapse_scale <= 1):
                raise ValueError("collapse_scale must be between 0 and 1")

        # Validate hierarchy_merge_scale if provided
        if collapse_scale is None and hierarchy_merge_scale is not None:
            if not (0 <= hierarchy_merge_scale <= 1):
                raise ValueError("hierarchy_merge_scale must be between 0 and 1")

        # Calculate the range of distances in Z for scaling
        distances = Z[:, 2]  # Third column of Z contains the distances
        min_dist = np.min(distances)
        max_dist = np.max(distances)
        dist_range = max_dist - min_dist
        if collapse_scale is not None:
            if collapse_scale==0:
                scaled_collapse_threshold = min_dist - 1e-6
            else:
                # Scale the collapse_scale (0 to 1) to the actual distance range
                scaled_collapse_threshold = min_dist + collapse_scale * dist_range if dist_range > 0 else min_dist
        

        if collapse_scale is None:
            # Full hierarchical structure (not strictly binary)
            from collections import defaultdict
            # Track nodes and their children
            node_children = defaultdict(list)
            node_types = {}

            # Add cell types as leaf nodes
            for ct in cell_types:
                G.add_node(ct, node_type='cell_type')
                node_types[ct] = 'cell_type'

            # Process the linkage matrix Z to build the hierarchy
            for i, row in enumerate(Z):
                left_idx, right_idx, distance, _ = row
                left_idx, right_idx = int(left_idx), int(right_idx)
                internal_node = f"internal_{i+N}"
                G.add_node(internal_node, node_type='internal')
                node_types[internal_node] = 'internal'

                # Identify children (could be leaf or internal nodes)
                left_node = cell_types[left_idx] if left_idx < N else f"internal_{left_idx}"
                right_node = cell_types[right_idx] if right_idx < N else f"internal_{right_idx}"

                # Add edges from parent to children
                G.add_edge(internal_node, left_node)
                G.add_edge(internal_node, right_node)

                # Store children for potential merging into multi-child nodes
                node_children[internal_node].extend([left_node, right_node])

            # Root is the last internal node
            root = f"internal_{2*N - 2}"
            
            # Mark this node specifically as root_internal to apply the root color
            G.nodes[root]['node_type'] = 'root_internal'

            # Optional: Collapse binary nodes into multi-child nodes (simulating non-binary hierarchy)
            if hierarchy_merge_scale is not None:
                # Scale the hierarchy_merge_scale to the distance range
                merge_threshold = hierarchy_merge_scale * dist_range if dist_range > 0 else 0
                distance_dict = {f"internal_{i+N}": row[2] for i, row in enumerate(Z)}
                for node in list(G.nodes()):
                    if node_types.get(node) == 'internal' and len(node_children[node]) == 2:
                        parent = next(iter(G.predecessors(node)), None)
                        if parent and abs(distance_dict.get(node, 0) - distance_dict.get(parent, 0)) < merge_threshold:
                            # Merge with parent if distances are within the scaled threshold
                            children = node_children[node]
                            G.remove_node(node)
                            for child in children:
                                G.add_edge(parent, child)
                            node_children[parent].extend(children)
                            node_children.pop(node)

        else:
            from collections import defaultdict
            # Use the scaled threshold in fcluster
            cluster_labels = fcluster(Z, t=scaled_collapse_threshold, criterion='distance')
            unique_clusters = np.unique(cluster_labels)
            root = "root"
            G.add_node(root, node_type='root')
            clusters = defaultdict(list)
            for ct, lbl in zip(cell_types, cluster_labels):
                clusters[lbl].append(ct)
            if len(unique_clusters) == 1:
                for ct in clusters[unique_clusters[0]]:
                    G.add_node(ct, node_type='cell_type')
                    G.add_edge(root, ct)
            else:
                for lbl, members in clusters.items():
                    if len(members) > 1:
                        cluster_node = f"cluster_{lbl}"
                        G.add_node(cluster_node, node_type='cluster')
                        G.add_edge(root, cluster_node)
                        for ct in members:
                            G.add_node(ct, node_type='cell_type')
                            G.add_edge(cluster_node, ct)
                    else:
                        if add_cluster_node_for_single_node_cluster:
                            cluster_node = f"cluster_{lbl}"
                            G.add_node(cluster_node, node_type='cluster')
                            G.add_edge(root, cluster_node)
                            ct = members[0]
                            G.add_node(ct, node_type='cell_type')
                            G.add_edge(cluster_node, ct)
                        else:
                            ct = members[0]
                            G.add_node(ct, node_type='cell_type')
                            G.add_edge(root, ct)
    
    # Check if a custom cell type to gene mapping is provided
    if map_cell_type_gene is not None:
        # Get a set of valid genes (those in adata.var_names)
        valid_genes = set(adata.var_names)
        
        # Dictionary to store selected genes for each cell type
        selected_genes_dict = {}
        all_selected_genes = set()
        
        # Process each cell type in the mapping
        mapped_cell_types = []
        for ct, genes in map_cell_type_gene.items():
            if ct not in cell_types:
                continue  # Skip cell types not in the data
            
            mapped_cell_types.append(ct)
            
            # Filter genes to keep only those in adata.var_names
            valid_ct_genes = [gene for gene in genes if gene in valid_genes]
            
            # Store the valid genes for this cell type
            selected_genes_dict[ct] = valid_ct_genes
            all_selected_genes.update(valid_ct_genes)
        
        # Check if there's any overlap between cell_types and mapped cell types
        if len(mapped_cell_types) == 0:
            if cell_type_selected is not None:
                raise ValueError(f"No overlap between cell types in map_cell_type_gene {list(map_cell_type_gene.keys())} and the selected cell types {cell_type_selected}")
            else:
                raise ValueError(f"None of the cell types in map_cell_type_gene {list(map_cell_type_gene.keys())} are present in the data")
        
        # Convert set to list for further processing
        selected_genes = list(all_selected_genes)
    else:
        # Extract top N marker genes for all cell types at once (original behavior)
        marker_genes_df = pd.DataFrame(adata.uns[cosg_key]['names']).iloc[:top_n_genes]  # Slice once for efficiency
        
        # Filter to only include selected cell types if needed
        if cell_type_selected is not None:
            cols_to_keep = [col for col in marker_genes_df.columns if col in valid_selected_cell_types]
            marker_genes_df = marker_genes_df[cols_to_keep]
        
        selected_genes = marker_genes_df.values.flatten()  # Flatten to get all genes as a 1D list
        selected_genes = pd.Index(selected_genes).dropna().unique()  # Remove NaNs & duplicates


    # Attach marker gene nodes to each cell type leaf
    gene_nodes = {}
    
    if map_cell_type_gene is not None:
        # Use the custom mapping
        for ct in cell_types:
            if ct not in G or ct not in selected_genes_dict or not selected_genes_dict[ct]:
                continue
            
            # Get the custom genes for this cell type
            ct_genes = selected_genes_dict[ct]
            
            for gene in ct_genes:
                marker_node = f"{ct}__gene__{gene}"
                # Get the COSG score if available, otherwise use 0
                score = cosg_score_df.loc[gene, ct] if gene in cosg_score_df.index and ct in cosg_score_df.columns else 0
                G.add_node(marker_node, node_type='gene', score=score, gene=gene)
                G.add_edge(ct, marker_node)
                gene_nodes[marker_node] = score
    else:
        # Original behavior using top marker genes
        marker_genes_df = pd.DataFrame(adata.uns[cosg_key]['names']).iloc[:top_n_genes]
        
        # Filter to only include selected cell types if needed
        if cell_type_selected is not None:
            cols_to_keep = [col for col in marker_genes_df.columns if col in valid_selected_cell_types]
            marker_genes_df = marker_genes_df[cols_to_keep]
        
        for ct in cell_types:
            if ct not in G or ct not in marker_genes_df.columns:
                continue

            # Get precomputed top N marker genes for this cell type
            top_genes = marker_genes_df[ct].dropna()  # Drop NaNs to avoid issues

            for gene in top_genes:
                marker_node = f"{ct}__gene__{gene}"
                score = cosg_score_df.loc[gene, ct] if gene in cosg_score_df.index else 0  # Fetch COSG score
                G.add_node(marker_node, node_type='gene', score=score, gene=gene)
                G.add_edge(ct, marker_node)
                gene_nodes[marker_node] = score
 

    
    ### Calculate the expression percentage
    filtered_cosg_score_df = cosg_score_df.loc[selected_genes] if len(selected_genes) > 0 else cosg_score_df  # Keep only selected marker genes
    gene_expr_percentage = _compute_gene_expression_percentage(adata, group_by, filtered_cosg_score_df, layer=layer)

    
    gene_node_sizes = {}
    for n, d in G.nodes(data=True):
        if d.get('node_type') == 'gene':
            ct = n.split('__gene__')[0]
            gene_name = d['gene']
            percentage = gene_expr_percentage.loc[gene_name, ct] if gene_name in gene_expr_percentage.index else 0
            gene_node_sizes[n] = gene_size_scale * (percentage / 100)
    
    # Set final node sizes
    node_sizes = {}
    for n, d in G.nodes(data=True):
        ntype = d.get('node_type', '')
        if ntype in ['internal', 'root', 'root_internal', 'cluster']:
            node_sizes[n] = 50
        elif ntype == 'cell_type':
            node_sizes[n] = 600
        else:
            # Default size for any other node types
            node_sizes[n] = 100
            
    # Add gene node sizes separately (which use the expression percentage)
    for n, d in G.nodes(data=True):
        if d.get('node_type') == 'gene':
            node_sizes[n] = gene_node_sizes.get(n, 300)
    
    # Compute radial layout
    pos = _radial_dendrogram_layout(G, root, radius_step=radius_step)
    
    # Set palette for cell type nodes if not provided
    if palette is None and f"{group_by}_colors" in adata.uns:
        # If cell types are filtered, we need to filter the palette too
        if cell_type_selected is not None and f"{group_by}_colors" in adata.uns:
            all_palette = adata.uns[f"{group_by}_colors"]
            if hasattr(adata.obs[group_by], 'cat'):
                all_categories = adata.obs[group_by].cat.categories

                if len(all_palette) == len(all_categories):
                    filtered_indices = [list(all_categories).index(ct) for ct in valid_selected_cell_types if ct in all_categories]
                    palette = [all_palette[i] for i in filtered_indices]
                else:
                    palette = all_palette
            else:
                palette = all_palette
            
        else:
            palette = adata.uns[f"{group_by}_colors"]
    
    # Drawing
    fig = plt.figure(figsize=figure_size)
    
    # Create a square main Axes for the network plot (slightly smaller to make room for legends)
    ax_main = fig.add_axes([0.05, 0.05, 0.75, 0.85])
    ax_main.set_aspect('equal')
    
    # Create legends on the right side
    legend_x = 0.82  # Common x position for legends
    
    # Create node type legend
    ax_ns = fig.add_axes([legend_x, 0.65, 0.15, 0.15])  # Node type legend at top
    
    # Create expression percentage legend
    ax_ds = fig.add_axes([legend_x, 0.525, 0.15, 0.15])  # Expression % legend in middle
    
    # Create colorbar Axes
    ax_cb = fig.add_axes([legend_x + 0.04, 0.15, colorbar_width, 0.25])  # Colorbar at bottom, shifted right

    ### Set up the color map for gene nodes
    cmap_obj = plt.get_cmap(cmap)
    if gene_nodes:
        scores_array = np.array(list(gene_nodes.values()))
        vmin = gene_color_min
        vmax = gene_color_max if gene_color_max is not None else scores_array.max()
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    else:
        norm = mcolors.Normalize(vmin=0, vmax=1)
    
    ### Setup the node colors
    node_colors = {}
    for n, d in G.nodes(data=True):
        ntype = d.get('node_type', '')
        if ntype == 'internal':
            node_colors[n] = color_internal_node
        elif ntype in ['root', 'root_internal']:  # Handle both explicit root and root_internal nodes
            node_colors[n] = color_root_node
        elif ntype == 'cluster':
            node_colors[n] = color_internal_node
        elif ntype == 'cell_type':
            if palette is not None:
                if isinstance(palette, dict):
                    node_colors[n] = palette.get(n, "lightblue")
                elif isinstance(palette, list):
                    try:
                        idx = cell_types.index(n)
                        node_colors[n] = palette[idx] if idx < len(palette) else "lightblue"
                    except ValueError:
                        node_colors[n] = "lightblue"
                else:
                    node_colors[n] = "lightblue"
            else:
                node_colors[n] = "lightblue"
        elif ntype == 'gene':
            node_colors[n] = cmap_obj(norm(d['score']))
        else:
            node_colors[n] = 'lightgrey'
            
    ### Setup the node shapes
    node_shapes = {}
    for n, d in G.nodes(data=True):
        ntype = d.get('node_type', '')
        if ntype == 'internal':
            node_shapes[n] = node_shape_internal
        elif ntype == 'cell_type':
            node_shapes[n] = node_shape_cell_type
        elif ntype == 'gene':
            node_shapes[n] = node_shape_gene
        else:
            node_shapes[n] = 'o'
    ### Draw nodes seprately, because they are using different shapes
    #### Draw cell type nodes with the specified shape
    node_list = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'cell_type']
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=node_list,
        node_shape=node_shape_cell_type,
        node_color=[node_colors[n] for n in node_list],
        node_size=[node_sizes[n] for n in node_list],
        ax=ax_main
    )
    
    #### Draw gene nodes with the specified shape
    node_list = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'gene']
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=node_list,
        node_shape=node_shape_gene,
        node_color=[node_colors[n] for n in node_list],
        node_size=[node_sizes[n] for n in node_list],
        ax=ax_main
    )
    
    #### Draw internal nodes with the specified shape
    node_list = [n for n, d in G.nodes(data=True) if d.get('node_type') not in ('gene', 'cell_type')]
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=node_list,
        node_shape=node_shape_internal,
        node_color=[node_colors[n] for n in node_list],
        node_size=[node_sizes[n] for n in node_list],
        ax=ax_main
    )
    
    ### Draw edges:
    if edge_curved == 0:
        # Use the default straight edges
        nx.draw_networkx_edges(G, pos, ax=ax_main, arrows=False, edge_color=color_edge)
    else:
        # Draw curved edges
        for (u, v) in G.edges():
            # Get the start and end positions
            x1, y1 = pos[u]
            x2, y2 = pos[v]
            
            # Create the curved connection
            # Calculate control point for quadratic bezier curve
            mid_x = (x1 + x2) / 2
            mid_y = (y1 + y2) / 2
            
            # Calculate normal vector to the line connecting the nodes
            dx = x2 - x1
            dy = y2 - y1
            length = np.sqrt(dx*dx + dy*dy)
            
            if length > 0:
                # Perpendicular vector with length scaled by edge_curved
                norm_x = -dy / length * edge_curved
                norm_y = dx / length * edge_curved
                
                # Control point by moving midpoint perpendicular to the line
                ctrl_x = mid_x + norm_x
                ctrl_y = mid_y + norm_y
                
                # Create a Path with a quadratic curve
                path = mpath.Path([(x1, y1), (ctrl_x, ctrl_y), (x2, y2)], [mpath.Path.MOVETO, mpath.Path.CURVE3, mpath.Path.CURVE3])
                patch = mpatches.PathPatch(path, fill=False, edgecolor=color_edge, lw=1)
                ax_main.add_patch(patch)
            else:
                # Fall back to straight line if points are the same
                ax_main.plot([x1, x2], [y1, y2], color=color_edge, lw=1)
    
    ### Plot the labels
    for n, d in G.nodes(data=True):
        if d.get('node_type') == 'gene':
            parents = list(G.predecessors(n))
            if not parents:
                continue
            parent = parents[0]
            x_parent, y_parent = pos[parent]
            x_gene, y_gene = pos[n]
            vec = np.array([x_gene - x_parent, y_gene - y_parent])
            norm_vec = vec / (np.linalg.norm(vec) + 1e-9)
            label_pos = (x_gene + gene_label_offset * norm_vec[0],
                         y_gene + gene_label_offset * norm_vec[1])
            angle = np.degrees(np.arctan2(norm_vec[1], norm_vec[0]))
            if angle > 90:
                angle -= 180
            elif angle < -90:
                angle += 180
            text_color = gene_label_color if gene_label_color is not None else node_colors[n]
            txt = ax_main.text(label_pos[0], label_pos[1], d['gene'],
                    fontsize=8, color=text_color,
                    rotation=angle,
                    horizontalalignment='center',
                    verticalalignment='center')
            txt.set_path_effects([PathEffects.withStroke(linewidth=font_outline, foreground='white')])
        ### Adjust the direction of cell type labels:
        elif d.get('node_type') == 'cell_type':
            parents = list(G.predecessors(n))
            if not parents:
                continue
            parent = parents[0]
            x_parent, y_parent = pos[parent]
            x_ct, y_ct = pos[n]
            vec = np.array([x_ct - x_parent, y_ct - y_parent])
            norm_vec = vec / (np.linalg.norm(vec) + 1e-9)
            label_pos = (x_ct + cell_type_label_offset * norm_vec[0],
                         y_ct + cell_type_label_offset * norm_vec[1])
            angle = np.degrees(np.arctan2(norm_vec[1], norm_vec[0]))
            if angle > 90:
                angle -= 180
            elif angle < -90:
                angle += 180
            txt = ax_main.text(label_pos[0], label_pos[1], n,
                               fontsize=10, color='black',
                               rotation=angle,
                               horizontalalignment='center',
                               verticalalignment='center')
            txt.set_path_effects([PathEffects.withStroke(linewidth=font_outline, foreground='white')])
    
    ### Set up the color bar
    sm = plt.cm.ScalarMappable(cmap=cmap_obj, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, cax=ax_cb, orientation='vertical')
    # Move the label to the left side of the colorbar
    cbar.set_label("COSG Score", fontsize=12, rotation=270, labelpad=15, ha='center')
    
    
    ### Calculate the dot size dynamically
    min_expr = gene_expr_percentage[gene_expr_percentage > 0].min().min() if not gene_expr_percentage.empty else 5  # Ignore 0% values
    max_expr = gene_expr_percentage.max().max() if not gene_expr_percentage.empty else 100  # Maximum expression percentage

    # Round min/max to nearest multiple of 5
    min_expr_rounded = np.floor(min_expr / 5) * 5
    max_expr_rounded = np.ceil(max_expr / 5) * 5

    # Generate up to 5 evenly spaced values, ignoring 0%
    num_circles = min(5, int((max_expr_rounded - min_expr_rounded) / 5) + 1)
    legend_percentages = np.linspace(min_expr_rounded, max_expr_rounded, num=num_circles)
    legend_percentages = np.unique(np.round(legend_percentages / 5) * 5).astype(int)  # Ensure multiples of 5

    # Generate dot size legend (skip 0%)
    legend_markers = [
        plt.Line2D([0], [0], marker=node_shape_gene, color='black', label=f' {p}%', 
                   markerfacecolor='white', markersize=np.sqrt(gene_size_scale * (p/100)))
        for p in legend_percentages if p > 0
    ]

    
    
    # Calculate the largest gene node size
    max_gene_size = max(gene_node_sizes.values()) if gene_node_sizes else gene_size_scale
    
    # Create node shape legend with just cell type and gene nodes, using appropriate sizes
    node_shape_markers = [
        plt.Line2D([0], [0], marker=node_shape_cell_type, color='black', label=' Cell type', 
                   markerfacecolor='white', markersize=np.sqrt(600/np.pi)),
        plt.Line2D([0], [0], marker=node_shape_gene, color='black', label=' Gene', 
                   markerfacecolor='white', markersize=np.sqrt(max_gene_size/np.pi))
    ]
    
    # Place the node shape legend in the upper-left corner of ax_ns
    ax_ns.legend(handles=node_shape_markers, title="Node type", loc='upper left',
                 frameon=False, fontsize=12, title_fontsize=12)
    ax_ns.axis('off')
    
    # Place the dot size legend in the upper-left corner of ax_ds
    ax_ds.legend(handles=legend_markers, title="Expression %", loc='upper left',
                 frameon=False, fontsize=12, title_fontsize=12)
    ax_ds.axis('off')
    
    # ax_main.set_title("Radial Dendrogram of Cell Types with Top Marker Genes", fontsize=12)
    ax_main.axis('off')

        
    ### Whether to show the figure or not
    if show_figure:
        plt.show()  # Explicitly display the figure

    ### Save the figure
    if save:
        fig.savefig(save, bbox_inches='tight')  # Save the figure to file
        print("Figure saved to: ", save)
        plt.close(fig)  # Close the figure to prevent display
    elif not show_figure:
        plt.close(fig)  # Close the figure if not showing or saving
        
 
    
### packaged dotplot function in COSG

def plotMarkerDotplot(
    adata: AnnData,
    groupby: str,
    top_n_genes: int = 3,
    use_rep: str = 'X_pca',
    layer: str = None,
    key_cosg: str = 'cosg',
    swap_axes: bool = False,
    standard_scale: str = 'var',
    cmap: str = 'Spectral_r',
    save: str = None,
    **dotplot_kwargs
):
    """
    Generate a dot plot of top marker genes identified by COSG.

    The function computes the cell cluster ordering using a dendrogram (if `use_rep` is provided)
    or derives it from `adata.obs[groupby]`, extracts the top marker genes identified by COSG, 
    and plots a dotplot using Scanpy's `sc.pl.dotplot`.

    Requires scanpy, which is an optional dependency of COSG::

        pip install 'cosg[dotplot]'

    It is the only function in COSG that needs it. `plotMarkerDendrogram` and
    `plotMarkerStream` do not.

    Parameters
    ----------
    adata
        Annotated data object that includes COSG results.
    groupby : str
        The cell group key in `adata.obs`, should match with the `groupby` parameter used in COSG.
    top_n_genes : int, optional (default: 3)
        The number of top marker genes to show for each group.
    use_rep : str, optional (default: 'X_pca')
        The cell low-dimensional representation key (e.g., PCA, UMAP) in `adata.obsm` used to compute the dendrogram.
    layer : str or None, optional
        The layer key to use for expression values in the dotplot (default: None).
    key_cosg : str, optional (default: 'cosg')
        The key in `adata.uns` where COSG results are stored.
    swap_axes : bool, optional (default: False)
        Whether to swap axes in the dot plot.
    standard_scale : str or None, optional
        Whether to standardize expression values across `'var'` (genes) or `'group'` (cell groups or clusters).
        Can be `'var'`, `'group'`, or `None` (default: `'var'`).
    cmap : str, optional (default: 'Spectral_r')
        The colormap used for the dot plot.
    save : str or None, optional
        If provided, saves the plot to a file. The filename should include an extension (e.g., `"cosg_markers.pdf"`).
    **dotplot_kwargs : dict
        Additional keyword arguments to pass to `sc.pl.dotplot`.

    Returns
    -------
    None
        Displays the dot plot.

    Raises
    ------
    ValueError
        If required COSG results or dendrogram information are missing, or if the provided `groupby`
        does not match the one stored in COSG parameters.

    Example
    -------
    >>> import scanpy as sc
    >>> import cosg  # Assuming plotDotPlot is part of the cosg package
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> # Using a specific low-dimensional representation for dendrogram computation and cell type ordering:
    >>> cosg.plotMarkerDotplot(
    ...     adata,
    ...     groupby='bulk_labels',
    ...     top_n_genes=3,
    ...     key_cosg='cosg',
    ...     use_rep='X_pca',
    ...     swap_axes=False,
    ...     standard_scale='var',
    ...     cmap='Spectral_r'
    ... )
    >>> # Deriving cell order from adata.obs when use_rep is None:
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
    
    # Check that COSG results are available in adata.uns using the specified key.
    if key_cosg not in adata.uns or 'names' not in adata.uns[key_cosg]:
        raise ValueError(f"COSG results not found in `adata.uns['{key_cosg}']['names']`.")
    
    # Check that the provided groupby matches the one stored in COSG parameters.
    if 'params' not in adata.uns[key_cosg] or 'groupby' not in adata.uns[key_cosg]['params']:
        raise ValueError(
            f"The COSG results in `adata.uns['{key_cosg}']` do not contain a 'groupby' parameter in 'params'."
        )
    if adata.uns[key_cosg]['params']['groupby'] != groupby:
        raise ValueError(
            f"Provided groupby '{groupby}' does not match the groupby used in COSG results "
            f"('{adata.uns[key_cosg]['params']['groupby']}')."
        )

    
    # Set the dendrogram key
    dendro_key = 'dendrogram_' + groupby

    # Compute dendrogram or derive ordering based on use_rep
    if use_rep is not None:
        # Temporarily suppress scanpy verbosity
        _dendrogram_order(adata, groupby=groupby, use_rep=use_rep)


        if dendro_key not in adata.uns or 'categories_ordered' not in adata.uns[dendro_key]:
            raise ValueError(
                f"Dendrogram results for groupby='{groupby}' not found in "
                f"`adata.uns['{dendro_key}']['categories_ordered']`."
            )
        ordering = adata.uns[dendro_key]['categories_ordered']
    else:
        # Derive ordering locally from adata.obs[groupby] without writing to adata.uns.
        if hasattr(adata.obs[groupby], "cat"):
            ordering = list(adata.obs[groupby].cat.categories)
        else:
            unique_values=adata.obs[groupby].unique()
            
            ### add a helper function here, if the cell clusters are "1", "2", ... , "10", "11", ...
            ### order them as "1", "2", ... , "10", "11", ..., instead of being "1", "10", "11", ..., "2", ...
            def _is_all_numeric(groups):
                try:
                    [float(x) for x in groups]
                    return True
                except ValueError:
                    return False

            if _is_all_numeric(unique_values):
                ordering = sorted(unique_values, key=lambda x: float(x))
            else:
                ordering = sorted(unique_values)
            
        
    # Extract the top_n_genes marker genes for each group from the COSG results.
    df_tmp = pd.DataFrame(adata.uns[key_cosg]['names'][:top_n_genes,]).T

    # Reorder rows based on the derived ordering.
    df_tmp = df_tmp.reindex(ordering)

    # Convert the DataFrame rows to a dictionary of marker genes per group.
    marker_genes_list = {idx: list(row.values) for idx, row in df_tmp.iterrows()}
    marker_genes_list = {k: v for k, v in marker_genes_list.items() if not any(isinstance(x, float) for x in v)}

    # Enable dendrogram only if use_rep is provided. When it is, the ordering
    # above already wrote adata.uns['dendrogram_<groupby>'] in scanpy's own
    # layout, so sc.pl.dotplot reads it rather than recomputing.
    use_dendrogram = use_rep is not None
    # Generate and display the dot plot with the provided parameters.
    sc = _require_scanpy("plotMarkerDotplot")
    sc.pl.dotplot(
        adata,
        marker_genes_list,
        groupby=groupby,
        layer=layer,
        dendrogram=use_dendrogram,
        swap_axes=swap_axes,
        standard_scale=standard_scale,
        cmap=cmap,
        save=save,
        **dotplot_kwargs
    )
    
    
#############Plot Marker Stream##########
"""
plotMarkerStream: Visualize COSG marker gene specificity scores as a streamgraph.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d, make_interp_spline
from scipy.signal import savgol_filter
from typing import Optional, Union, Tuple, List, Dict
import warnings


def _validate_inputs(adata, groupby, cosg_scores, cosg_key, plot_type, smooth_method,
                     smooth_edge_mode, bracket_position, gene_label_position,
                     gene_label_mode, legend_mode, legend_loc, genes, orientation):
    """Validate all input parameters."""
    if groupby not in adata.obs.columns:
        raise ValueError(f"groupby '{groupby}' not found in adata.obs.")
    
    valid_plot_types = ['proportional', 'stacked', 'mirror', 'ridge']
    if plot_type not in valid_plot_types:
        raise ValueError(f"Invalid plot_type. Must be one of: {valid_plot_types}")
    
    valid_smooth_methods = ['gaussian', 'spline', 'savgol']
    if smooth_method not in valid_smooth_methods:
        raise ValueError(f"Invalid smooth_method. Must be one of: {valid_smooth_methods}")
    
    valid_smooth_edge_modes = ['constant', 'taper']
    if smooth_edge_mode not in valid_smooth_edge_modes:
        raise ValueError(f"Invalid smooth_edge_mode. Must be one of: {valid_smooth_edge_modes}")
    
    valid_bracket_positions = ['top', 'bottom', 'left', 'right']
    if bracket_position not in valid_bracket_positions:
        raise ValueError(f"Invalid bracket_position. Must be one of: {valid_bracket_positions}")
    
    valid_gene_label_positions = ['top', 'bottom', 'left', 'right']
    if gene_label_position not in valid_gene_label_positions:
        raise ValueError(f"Invalid gene_label_position. Must be one of: {valid_gene_label_positions}")
    
    valid_gene_label_modes = ['all', 'first', 'none']
    if gene_label_mode not in valid_gene_label_modes:
        raise ValueError(f"Invalid gene_label_mode. Must be one of: {valid_gene_label_modes}")
    
    valid_legend_modes = ['horizontal', 'vertical', 'none']
    if legend_mode not in valid_legend_modes:
        raise ValueError(f"Invalid legend_mode. Must be one of: {valid_legend_modes}")
    
    valid_legend_locs = ['bottom', 'right', 'top', 'left']
    if legend_loc not in valid_legend_locs:
        raise ValueError(f"Invalid legend_loc. Must be one of: {valid_legend_locs}")
    
    valid_orientations = ['vertical', 'horizontal']
    if orientation not in valid_orientations:
        raise ValueError(f"Invalid orientation. Must be one of: {valid_orientations}")
    
    _validate_genes_parameter(genes)
    
    if cosg_scores is None and genes is None:
        if cosg_key not in adata.uns:
            raise ValueError(f"COSG results not found in adata.uns['{cosg_key}'].")
        if 'COSG' not in adata.uns[cosg_key]:
            raise ValueError(f"adata.uns['{cosg_key}']['COSG'] not found.")


def _validate_genes_parameter(genes):
    """Validate the genes parameter."""
    if genes is None:
        return
    
    if not isinstance(genes, (list, dict)):
        raise TypeError(
            f"'genes' must be None, list, or dict. Got {type(genes).__name__}.\n"
            f"Examples:\n"
            f"  genes=['Gene1', 'Gene2']  # No brackets\n"
            f"  genes={{'Group1': ['Gene1'], 'Group2': ['Gene2']}}  # With brackets"
        )
    
    if isinstance(genes, list):
        if len(genes) == 0:
            raise ValueError("'genes' list cannot be empty.")
        non_strings = [g for g in genes if not isinstance(g, str)]
        if non_strings:
            raise TypeError(f"All items in 'genes' list must be strings. Found: {non_strings[:5]}")
    
    if isinstance(genes, dict):
        if len(genes) == 0:
            raise ValueError("'genes' dictionary cannot be empty.")
        for key, value in genes.items():
            if not isinstance(key, str):
                raise TypeError(f"Dict keys must be strings. Found: {type(key).__name__}")
            if not isinstance(value, (list, tuple)):
                raise TypeError(f"Dict values must be lists. For '{key}', got {type(value).__name__}")
            if len(value) == 0:
                raise ValueError(f"Gene list for group '{key}' cannot be empty.")
            non_strings = [g for g in value if not isinstance(g, str)]
            if non_strings:
                raise TypeError(f"All genes in '{key}' must be strings. Found: {non_strings[:3]}")


def _process_genes_parameter(genes, adata, cosg_score_df, show_brackets):
    """Process the genes parameter and validate gene existence."""
    if genes is None:
        return None, None, False, show_brackets
    
    available_genes = set(cosg_score_df.index)
    if hasattr(adata, 'var_names'):
        available_genes = available_genes.union(set(adata.var_names))
    
    if isinstance(genes, list):
        gene_list, missing = [], []
        for gene in genes:
            if gene in available_genes:
                gene_list.append(gene)
            else:
                missing.append(gene)
        if missing:
            warnings.warn(f"Genes not found: {missing[:10]}{'...' if len(missing) > 10 else ''}")
        if len(gene_list) == 0:
            raise ValueError("No valid genes found.")
        return gene_list, None, True, False
    
    elif isinstance(genes, dict):
        gene_dict, all_missing = {}, []
        for group, group_genes in genes.items():
            valid = [g for g in group_genes if g in available_genes]
            missing = [g for g in group_genes if g not in available_genes]
            all_missing.extend([(group, g) for g in missing])
            if valid:
                gene_dict[group] = valid
            else:
                warnings.warn(f"Group '{group}' has no valid genes.")
        if all_missing:
            warnings.warn(f"Genes not found: {all_missing[:10]}{'...' if len(all_missing) > 10 else ''}")
        if len(gene_dict) == 0:
            raise ValueError("No valid groups found.")
        return None, gene_dict, True, show_brackets
    
    return None, None, False, show_brackets


def _resolve_cosg_scores(adata, cosg_scores, cosg_key, groupby, genes):
    """Resolve COSG scores from user input or compute from adata.uns."""
    if cosg_scores is not None:
        cosg_score_df = cosg_scores.copy()
    else:
        if cosg_key in adata.uns and 'COSG' in adata.uns[cosg_key]:
            cosg_df = indexByGene(adata.uns[cosg_key]['COSG'], set_nan_to_zero=True, convert_negative_one_to_zero=True)
            cosg_score_df = iqrLogNormalize(cosg_df)
        elif genes is not None:
            warnings.warn("No COSG results. Creating uniform scores.")
            all_genes = genes if isinstance(genes, list) else [g for gl in genes.values() for g in gl]
            celltypes = list(adata.obs[groupby].cat.categories) if hasattr(adata.obs[groupby], 'cat') else list(adata.obs[groupby].unique())
            return pd.DataFrame(np.ones((len(all_genes), len(celltypes))), index=all_genes, columns=celltypes)
        else:
            raise ValueError(f"COSG results not found in adata.uns['{cosg_key}'].")
    
    if cosg_score_df.isna().any().any():
        cosg_score_df = cosg_score_df.fillna(0)
    if (cosg_score_df < 0).any().any():
        cosg_score_df = cosg_score_df.abs()
    return cosg_score_df


def _order_by_specificity_sum(cosg_score_df, n_top_genes, celltypes):
    """Order cell types by sum of their top N marker genes' scores."""
    sums = {}
    for ct in celltypes:
        if ct not in cosg_score_df.columns:
            sums[ct] = 0
        else:
            sums[ct] = cosg_score_df[ct].nlargest(n_top_genes).sum()
    return sorted(celltypes, key=lambda x: sums[x], reverse=True)


def _resolve_celltype_order(adata, cosg_score_df, groupby, celltype_order, use_dendrogram,
                            dendrogram_key, use_rep, calc_dendro_on_cosg, linkage_method,
                            distance_metric, order_by_specificity, n_top_genes, gene_dict):
    """Resolve cell type ordering."""
    if gene_dict is not None:
        if celltype_order is not None:
            ordering = [ct for ct in celltype_order if ct in gene_dict]
            if not ordering:
                ordering = list(gene_dict.keys())
        else:
            ordering = list(gene_dict.keys())
        return ordering
    
    if celltype_order is not None:
        ordering = list(celltype_order)
    elif order_by_specificity:
        celltypes = list(adata.obs[groupby].cat.categories) if hasattr(adata.obs[groupby], 'cat') else list(cosg_score_df.columns)
        celltypes = [ct for ct in celltypes if ct in cosg_score_df.columns]
        ordering = _order_by_specificity_sum(cosg_score_df, n_top_genes, celltypes)
    elif use_dendrogram:
        dendro_key = dendrogram_key or f'dendrogram_{groupby}'
        if dendro_key in adata.uns and 'categories_ordered' in adata.uns[dendro_key]:
            ordering = list(adata.uns[dendro_key]['categories_ordered'])
        elif calc_dendro_on_cosg:
            data = cosg_score_df.T.values
            if data.shape[0] < 2:
                ordering = list(cosg_score_df.columns)
            else:
                D = pdist(data, metric=distance_metric)
                Z = linkage(D, method=linkage_method)
                ordering = [cosg_score_df.columns[i] for i in leaves_list(Z)]
        else:
            info = _dendrogram_order(adata, groupby=groupby, use_rep=use_rep)
            ordering = list(info['categories_ordered'])
    else:
        ordering = list(adata.obs[groupby].cat.categories) if hasattr(adata.obs[groupby], 'cat') else list(cosg_score_df.columns)
    
    valid = [ct for ct in ordering if ct in cosg_score_df.columns]
    if not valid:
        raise ValueError("No valid cell types found.")
    return valid


def _resolve_colors(adata, groupby, ordering, celltype_colors, cmap, gene_dict):
    """Resolve colors for each cell type or group."""
    colors = {}
    if celltype_colors is not None:
        colors = dict(celltype_colors)
    elif gene_dict is None and f'{groupby}_colors' in adata.uns:
        cats = list(adata.obs[groupby].cat.categories)
        palette = adata.uns[f'{groupby}_colors']
        for i, cat in enumerate(cats):
            if i < len(palette):
                colors[cat] = palette[i]
    
    missing = [ct for ct in ordering if ct not in colors]
    if missing:
        cmap_obj = plt.get_cmap(cmap)
        for i, ct in enumerate(missing):
            colors[ct] = cmap_obj(i / max(1, len(missing) - 1))
    
    for ct in ordering:
        if ct not in colors:
            colors[ct] = 'gray'
    return colors


def _select_genes_and_build_matrix(adata, cosg_score_df, ordering, cosg_key, n_top_genes,
                                    cosg_scores_provided, gene_list=None, gene_dict=None):
    """Select genes and build score matrix."""
    if gene_list is not None:
        selected_genes = gene_list
        gene_celltype_map = ['custom'] * len(gene_list)
        gene_positions = list(range(len(gene_list)))
        celltype_boundaries = [0, len(gene_list)]
        score_columns = list(cosg_score_df.columns)
        score_matrix = np.zeros((len(selected_genes), len(score_columns)))
        for i, gene in enumerate(selected_genes):
            if gene in cosg_score_df.index:
                score_matrix[i, :] = cosg_score_df.loc[gene, score_columns].values
        return selected_genes, gene_celltype_map, gene_positions, celltype_boundaries, score_matrix, score_columns
    
    elif gene_dict is not None:
        selected_genes, gene_celltype_map, gene_positions = [], [], []
        position = 0
        celltype_boundaries = [0]
        for group in ordering:
            if group not in gene_dict:
                continue
            for gene in gene_dict[group]:
                selected_genes.append(gene)
                gene_celltype_map.append(group)
                gene_positions.append(position)
                position += 1
            celltype_boundaries.append(position)
        score_columns = ordering
        score_matrix = np.zeros((len(selected_genes), len(score_columns)))
        for i, gene in enumerate(selected_genes):
            if gene in cosg_score_df.index:
                for j, ct in enumerate(score_columns):
                    if ct in cosg_score_df.columns:
                        score_matrix[i, j] = cosg_score_df.loc[gene, ct]
                    elif ct == gene_celltype_map[i]:
                        score_matrix[i, j] = 1.0
        return selected_genes, gene_celltype_map, gene_positions, celltype_boundaries, score_matrix, score_columns
    
    else:
        marker_names = None
        if cosg_key in adata.uns and 'names' in adata.uns[cosg_key]:
            raw = adata.uns[cosg_key]['names']
            if hasattr(raw, 'dtype') and raw.dtype.names:
                marker_names = pd.DataFrame({n: raw[n] for n in raw.dtype.names})
            else:
                marker_names = pd.DataFrame(raw)
            marker_names = marker_names.iloc[:n_top_genes]
        
        if marker_names is None:
            marker_names = pd.DataFrame({ct: cosg_score_df[ct].nlargest(n_top_genes).index.tolist() for ct in ordering})
        
        selected_genes, gene_celltype_map, gene_positions = [], [], []
        position = 0
        celltype_boundaries = [0]
        
        for ct in ordering:
            if ct in marker_names.columns:
                col = marker_names[ct]
                top_genes = col.iloc[:n_top_genes].dropna().tolist() if hasattr(col, 'iloc') else [g for g in list(col[:n_top_genes]) if pd.notna(g)]
            else:
                top_genes = cosg_score_df[ct].nlargest(n_top_genes).index.tolist()
            
            for gene in top_genes:
                selected_genes.append(gene)
                gene_celltype_map.append(ct)
                gene_positions.append(position)
                position += 1
            celltype_boundaries.append(position)
        
        if not selected_genes:
            raise ValueError("No genes selected.")
        
        score_matrix = np.zeros((len(selected_genes), len(ordering)))
        for i, gene in enumerate(selected_genes):
            for j, ct in enumerate(ordering):
                if gene in cosg_score_df.index and ct in cosg_score_df.columns:
                    score_matrix[i, j] = cosg_score_df.loc[gene, ct]
        
        return selected_genes, gene_celltype_map, gene_positions, celltype_boundaries, score_matrix, ordering


def _normalize_scores(score_matrix, plot_type):
    """Normalize scores based on plot type."""
    if plot_type == 'stacked':
        return score_matrix.copy()
    row_sums = score_matrix.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    return score_matrix / row_sums


def _create_step_data(x_original, y_data):
    """Convert point data to step function."""
    n_genes = len(x_original)
    n_celltypes = y_data.shape[1]
    x_step = np.zeros(n_genes * 2)
    y_step = np.zeros((n_genes * 2, n_celltypes))
    for i in range(n_genes):
        x_step[2*i] = i - 0.5
        x_step[2*i + 1] = i + 0.5
        y_step[2*i, :] = y_data[i, :]
        y_step[2*i + 1, :] = y_data[i, :]
    return x_step, y_step


def _round_discrete_edges(x_step, y_step, sigma):
    """Apply smoothing to soften column edges."""
    n_points = len(x_step)
    n_celltypes = y_step.shape[1]
    n_smooth = n_points * 5
    x_smooth = np.linspace(x_step[0], x_step[-1], n_smooth)
    y_smooth = np.zeros((n_smooth, n_celltypes))
    for j in range(n_celltypes):
        f = interp1d(x_step, y_step[:, j], kind='linear', fill_value='extrapolate')
        y_smooth[:, j] = gaussian_filter1d(f(x_smooth), sigma=sigma * (n_smooth / n_points))
    return x_smooth, np.maximum(y_smooth, 0)


def _apply_continuous_smoothing(x_original, y_data, smooth_method, smooth_sigma, plot_type, smooth_edge_mode):
    """Apply continuous smoothing for flowing curves."""
    n_genes = len(x_original)
    n_celltypes = y_data.shape[1]
    
    if n_genes < 2:
        x_ext = np.array([-0.5, 0, 0.5])
        if smooth_edge_mode == 'constant':
            y_ext = np.vstack([y_data[0:1], y_data, y_data[-1:]])
        else:
            y_ext = np.vstack([np.zeros((1, n_celltypes)), y_data, np.zeros((1, n_celltypes))])
        return x_ext, y_ext
    
    n_smooth = max(n_genes * 5, 100)
    x_smooth = np.linspace(-0.5, n_genes - 0.5, n_smooth)
    x_ext = np.concatenate([[-0.5], x_original.astype(float), [n_genes - 0.5]])
    
    if smooth_edge_mode == 'constant':
        y_ext = np.vstack([y_data[0:1], y_data, y_data[-1:]])
    else:
        y_ext = np.vstack([np.zeros((1, n_celltypes)), y_data, np.zeros((1, n_celltypes))])
    
    smoothed = np.zeros((n_smooth, n_celltypes))
    for j in range(n_celltypes):
        y = y_ext[:, j]
        if smooth_method == 'gaussian':
            f = interp1d(x_ext, y, kind='linear', fill_value='extrapolate')
            smoothed[:, j] = gaussian_filter1d(f(x_smooth), sigma=smooth_sigma * (n_smooth / n_genes))
        elif smooth_method == 'spline':
            try:
                spline = make_interp_spline(x_ext, y, k=min(3, len(x_ext) - 1))
                smoothed[:, j] = spline(x_smooth)
            except:
                f = interp1d(x_ext, y, kind='linear', fill_value='extrapolate')
                smoothed[:, j] = f(x_smooth)
        else:
            f = interp1d(x_ext, y, kind='linear', fill_value='extrapolate')
            y_interp = f(x_smooth)
            window = min(max(5, int(n_smooth * 0.1) | 1), n_smooth)
            if window % 2 == 0:
                window -= 1
            smoothed[:, j] = savgol_filter(y_interp, window, min(3, window - 2))
    
    smoothed = np.maximum(smoothed, 0)
    if plot_type in ['proportional', 'mirror', 'ridge']:
        row_sums = smoothed.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1
        smoothed = smoothed / row_sums
    return x_smooth, smoothed


def _process_coordinates(x_original, y_data, discrete, discrete_rounded, discrete_round_sigma,
                         smooth, smooth_method, smooth_sigma, smooth_edge_mode, plot_type):
    """Process coordinates based on visualization mode."""
    n_genes = len(x_original)
    n_celltypes = y_data.shape[1]
    
    if discrete:
        x_step, y_step = _create_step_data(x_original, y_data)
        if discrete_rounded:
            x_plot, y_plot = _round_discrete_edges(x_step, y_step, discrete_round_sigma)
            if plot_type in ['proportional', 'mirror', 'ridge']:
                row_sums = y_plot.sum(axis=1, keepdims=True)
                row_sums[row_sums == 0] = 1
                y_plot = y_plot / row_sums
        else:
            x_plot, y_plot = x_step, y_step
    else:
        if smooth:
            x_plot, y_plot = _apply_continuous_smoothing(x_original, y_data, smooth_method, smooth_sigma, plot_type, smooth_edge_mode)
        else:
            x_ext = np.concatenate([[-0.5], x_original.astype(float), [n_genes - 0.5]])
            if smooth_edge_mode == 'constant':
                y_ext = np.vstack([y_data[0:1], y_data, y_data[-1:]])
            else:
                y_ext = np.vstack([np.zeros((1, n_celltypes)), y_data, np.zeros((1, n_celltypes))])
            x_plot, y_plot = x_ext, y_ext
    return x_plot, y_plot


def _draw_vertical_stream(ax, x_plot, y_plot, ordering, colors, plot_type):
    """Draw vertical stackplot."""
    y_arrays = [y_plot[:, j] for j in range(len(ordering))]
    color_list = [colors.get(ct, 'gray') for ct in ordering]
    baseline_map = {'proportional': 'zero', 'stacked': 'zero', 'mirror': 'sym', 'ridge': 'wiggle'}
    ax.stackplot(x_plot, *y_arrays, labels=ordering, colors=color_list,
                 baseline=baseline_map[plot_type], edgecolor='white', linewidth=0.3)


def _draw_vertical_brackets(ax, ordering, celltype_boundaries, bracket_position, bracket_height,
                            bracket_text_rotation, bracket_linewidth, celltype_fontsize, y_min, y_max):
    """Draw vertical brackets."""
    y_range = y_max - y_min if y_max != y_min else 1.0
    if bracket_position == 'top':
        bracket_base, arm_dir, text_va = y_max, 1, 'bottom'
    else:
        bracket_base, arm_dir, text_va = y_min, -1, 'top'
    
    arm_length = bracket_height * y_range * 0.3
    text_offset = bracket_height * y_range * 0.5
    
    for i in range(len(ordering)):
        if i + 1 >= len(celltype_boundaries):
            continue
        start, end = celltype_boundaries[i], celltype_boundaries[i + 1]
        if start >= end:
            continue
        
        left_x, right_x = start - 0.5, end - 0.5
        mid_x = (left_x + right_x) / 2
        bracket_top = bracket_base + arm_dir * arm_length
        
        ax.plot([left_x, left_x], [bracket_base, bracket_top], color='black', lw=bracket_linewidth, solid_capstyle='butt', clip_on=False)
        ax.plot([left_x, right_x], [bracket_top, bracket_top], color='black', lw=bracket_linewidth, solid_capstyle='butt', clip_on=False)
        ax.plot([right_x, right_x], [bracket_base, bracket_top], color='black', lw=bracket_linewidth, solid_capstyle='butt', clip_on=False)
        
        text_y = bracket_base + arm_dir * (arm_length + text_offset)
        rotation = bracket_text_rotation if len(ordering[i]) > 8 else 0
        ax.text(mid_x, text_y, ordering[i], ha='center', va=text_va, rotation=rotation,
                fontsize=celltype_fontsize, fontweight='medium', clip_on=False)
    
    return arm_length + text_offset * 1.5


def _draw_gene_labels_with_offset(ax, selected_genes, n_genes, gene_label_offset, gene_label_mode,
                                   celltype_boundaries, effective_gene_position, gene_tick_rotation,
                                   gene_fontsize, gene_label_ha, gene_label_va):
    """Draw gene labels with offset (labels only, not ticks)."""
    ax.set_xticks(list(range(n_genes)))
    ax.set_xticklabels([])
    
    if gene_label_mode == 'all':
        genes_to_label = [(i, selected_genes[i]) for i in range(n_genes)]
    elif gene_label_mode == 'first':
        genes_to_label = [(b, selected_genes[b]) for b in celltype_boundaries[:-1] if b < n_genes]
    else:
        return
    
    if effective_gene_position == 'top':
        actual_rotation = -gene_tick_rotation
        default_ha = 'left' if abs(gene_tick_rotation) >= 45 else 'center'
        default_va = 'bottom'
        y_pos = 1.03  # v4.4: More space above ticks
    else:
        actual_rotation = gene_tick_rotation
        default_ha = 'right' if gene_tick_rotation >= 45 else 'center'
        default_va = 'top'
        y_pos = -0.03  # v4.4: More space below ticks
    
    ha = gene_label_ha if gene_label_ha else default_ha
    va = gene_label_va if gene_label_va else default_va
    
    for i, gene in genes_to_label:
        x_pos = i + gene_label_offset
        ax.text(x_pos, y_pos, gene, ha=ha, va=va, rotation=actual_rotation,
                fontsize=gene_fontsize, transform=ax.get_xaxis_transform(), clip_on=False)


def _draw_horizontal_bars(ax, score_matrix, ordering, colors, n_genes):
    """Draw horizontal stacked bars."""
    y_positions = np.arange(n_genes)
    left = np.zeros(n_genes)
    for j, ct in enumerate(ordering):
        widths = score_matrix[:, j]
        ax.barh(y_positions, widths, left=left, height=0.8, color=colors.get(ct, 'gray'),
                edgecolor='white', linewidth=0.3, label=ct)
        left += widths


def _draw_horizontal_fill(ax, y_coords, x_data, ordering, colors):
    """Draw horizontal filled areas."""
    cumsum = np.zeros(len(y_coords))
    for j, ct in enumerate(ordering):
        x_vals = x_data[:, j]
        ax.fill_betweenx(y_coords, cumsum, cumsum + x_vals, color=colors.get(ct, 'gray'),
                         edgecolor='white', linewidth=0.3, label=ct)
        cumsum += x_vals


def _draw_horizontal_brackets(ax, ordering, celltype_boundaries, bracket_position, bracket_height,
                              bracket_text_rotation, bracket_linewidth, celltype_fontsize, x_min, x_max):
    """Draw horizontal brackets."""
    x_range = x_max - x_min if x_max != x_min else 1.0
    if bracket_position == 'right':
        bracket_base, arm_dir, text_ha = x_max, 1, 'left'
    else:
        bracket_base, arm_dir, text_ha = x_min, -1, 'right'
    
    arm_length = bracket_height * x_range * 0.3
    text_offset = bracket_height * x_range * 0.5
    
    for i in range(len(ordering)):
        if i + 1 >= len(celltype_boundaries):
            continue
        start, end = celltype_boundaries[i], celltype_boundaries[i + 1]
        if start >= end:
            continue
        
        top_y, bottom_y = start - 0.5, end - 0.5
        mid_y = (top_y + bottom_y) / 2
        bracket_edge = bracket_base + arm_dir * arm_length
        
        ax.plot([bracket_base, bracket_edge], [top_y, top_y], color='black', lw=bracket_linewidth, solid_capstyle='butt', clip_on=False)
        ax.plot([bracket_edge, bracket_edge], [top_y, bottom_y], color='black', lw=bracket_linewidth, solid_capstyle='butt', clip_on=False)
        ax.plot([bracket_base, bracket_edge], [bottom_y, bottom_y], color='black', lw=bracket_linewidth, solid_capstyle='butt', clip_on=False)
        
        text_x = bracket_base + arm_dir * (arm_length + text_offset)
        rotation = bracket_text_rotation if len(ordering[i]) > 8 else 0
        ax.text(text_x, mid_y, ordering[i], ha=text_ha, va='center', rotation=rotation,
                fontsize=celltype_fontsize, fontweight='medium', clip_on=False)
    
    return arm_length + text_offset * 1.5


def _calculate_legend_offset(gene_tick_rotation, gene_fontsize, effective_gene_position, legend_offset):
    """Calculate total offset for legend placement."""
    if effective_gene_position in ['top', 'left']:
        return legend_offset
    if abs(gene_tick_rotation) >= 60:
        label_space = 0.18 + (gene_fontsize - 8) * 0.01
    elif abs(gene_tick_rotation) >= 30:
        label_space = 0.12 + (gene_fontsize - 8) * 0.01
    else:
        label_space = 0.06 + (gene_fontsize - 8) * 0.01
    return label_space + legend_offset


def _setup_legend(ax, fig, ordering, colors, legend_mode, legend_loc, legend_ncol,
                  legend_offset, legend_fontsize, effective_gene_position,
                  gene_tick_rotation, gene_fontsize, orientation):
    """Set up the legend."""
    if legend_mode == 'none':
        return
    
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=colors.get(ct, 'gray'), edgecolor='white', linewidth=0.5) for ct in ordering]
    
    if legend_ncol is None:
        legend_ncol = min(8, max(4, len(ordering) // 3)) if legend_loc == 'bottom' else 1
    
    total_offset = _calculate_legend_offset(gene_tick_rotation, gene_fontsize, effective_gene_position, legend_offset)
    
    if legend_loc == 'bottom':
        ax.legend(handles, ordering, loc='upper center', bbox_to_anchor=(0.5, -total_offset),
                  ncol=legend_ncol, frameon=False, fontsize=legend_fontsize,
                  handlelength=1.5, handleheight=1, columnspacing=1)
    elif legend_loc == 'right':
        ax.legend(handles, ordering, loc='center left', bbox_to_anchor=(1.02, 0.5),
                  ncol=1, frameon=False, fontsize=legend_fontsize)
    elif legend_loc == 'top':
        ax.legend(handles, ordering, loc='lower center', bbox_to_anchor=(0.5, 1.02),
                  ncol=legend_ncol, frameon=False, fontsize=legend_fontsize)
    elif legend_loc == 'left':
        ax.legend(handles, ordering, loc='center right', bbox_to_anchor=(-0.02, 0.5),
                  ncol=1, frameon=False, fontsize=legend_fontsize)


def plotMarkerStream(
    adata: AnnData,
    groupby: str,
    cosg_scores: Optional[pd.DataFrame] = None,
    cosg_key: str = 'cosg',
    n_top_genes: int = 3,
    genes: Optional[Union[List[str], Dict[str, List[str]]]] = None,
    celltype_order: Optional[List[str]] = None,
    use_dendrogram: bool = True,
    dendrogram_key: Optional[str] = None,
    use_rep: str = 'X_pca',
    calculate_dendrogram_on_cosg_scores: bool = True,
    linkage_method: str = 'ward',
    distance_metric: str = 'euclidean',
    order_by_specificity: bool = False,
    celltype_colors: Optional[Dict[str, str]] = None,
    cmap: str = 'Spectral_r',
    plot_type: str = 'stacked',
    orientation: str = 'vertical',
    discrete: bool = True,
    discrete_rounded: bool = False,
    discrete_round_sigma: float = 0.3,
    smooth: bool = False,
    smooth_method: str = 'gaussian',
    smooth_sigma: float = 0.8,
    smooth_edge_mode: str = 'constant',
    show_brackets: bool = True,
    bracket_position: str = 'top',
    bracket_height: float = 0.08,
    bracket_text_rotation: float = 45,
    bracket_linewidth: float = 1.0,
    show_boundaries: bool = False,
    boundary_color: str = 'white',
    boundary_linewidth: float = 0.5,
    show_gene_labels: bool = True,
    gene_label_position: str = 'bottom',
    gene_label_mode: str = 'all',
    gene_tick_rotation: float = 90,
    gene_fontsize: float = 9,
    gene_label_ha: Optional[str] = None,
    gene_label_va: Optional[str] = None,
    gene_label_offset: float = 0.1,
    celltype_fontsize: float = 10,
    swap_axes_labels: bool = False,
    show_grid: bool = False,
    legend_mode: str = 'horizontal',
    legend_loc: str = 'bottom',
    legend_ncol: Optional[int] = None,
    legend_offset: float = 0.25,
    legend_fontsize: float = 8,
    figsize: Optional[Tuple[float, float]] = None,
    ax: Optional[plt.Axes] = None,
    title: Optional[str] = None,
    show: bool = True,
    save: Optional[str] = None,
    return_fig: bool = False,
    return_data: bool = False,
) -> Optional[Union[plt.Figure, Tuple[plt.Figure, pd.DataFrame]]]:
    """
    Visualize COSG marker gene specificity scores as a streamgraph.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data object.
    groupby : str
        Cell type column in adata.obs.
    cosg_scores : DataFrame, optional
        Pre-computed scores.
    cosg_key : str
        Key in adata.uns.
    n_top_genes : int
        Top markers per cell type.
    genes : list or dict, optional
        Custom genes: list (no brackets) or dict (with brackets).
    orientation : str
        'vertical' (default) or 'horizontal'.
    gene_label_offset : float
        Horizontal offset for labels (only moves labels, not ticks).
    
    Returns
    -------
    None, Figure, DataFrame, or tuple
    """
    # Validation
    _validate_inputs(adata, groupby, cosg_scores, cosg_key, plot_type, smooth_method,
                     smooth_edge_mode, bracket_position, gene_label_position,
                     gene_label_mode, legend_mode, legend_loc, genes, orientation)
    
    if not hasattr(adata.obs[groupby], 'cat'):
        adata.obs[groupby] = adata.obs[groupby].astype('category')
    
    # Resolve scores
    cosg_score_df = _resolve_cosg_scores(adata, cosg_scores, cosg_key, groupby, genes)
    
    # Process genes
    gene_list, gene_dict, use_custom_genes, effective_show_brackets = _process_genes_parameter(
        genes, adata, cosg_score_df, show_brackets
    )
    
    # Ordering
    ordering = _resolve_celltype_order(adata, cosg_score_df, groupby, celltype_order, use_dendrogram,
                                       dendrogram_key, use_rep, calculate_dendrogram_on_cosg_scores,
                                       linkage_method, distance_metric, order_by_specificity,
                                       n_top_genes, gene_dict)
    
    # Colors
    colors = _resolve_colors(adata, groupby, ordering, celltype_colors, cmap, gene_dict)
    
    # Gene selection
    (selected_genes, gene_celltype_map, gene_positions, celltype_boundaries,
     score_matrix, score_columns) = _select_genes_and_build_matrix(
        adata, cosg_score_df, ordering, cosg_key, n_top_genes,
        cosg_scores is not None, gene_list, gene_dict
    )
    n_genes = len(selected_genes)
    
    # Normalization
    normalized = _normalize_scores(score_matrix, plot_type)
    
    # Coordinates
    x_original = np.arange(n_genes)
    x_plot, y_plot = _process_coordinates(x_original, normalized, discrete, discrete_rounded,
                                           discrete_round_sigma, smooth, smooth_method, smooth_sigma,
                                           smooth_edge_mode, plot_type)
    
    # Effective positions
    if orientation == 'horizontal':
        effective_gene_position = 'left'
        effective_bracket_position = 'right'
        effective_legend_loc = legend_loc if legend_loc != 'bottom' else 'right'
    else:
        if swap_axes_labels:
            effective_gene_position = 'top'
            effective_bracket_position = 'bottom'
        else:
            effective_gene_position = gene_label_position
            effective_bracket_position = bracket_position
        effective_legend_loc = legend_loc
    
    # Figure setup
    if figsize is None:
        if orientation == 'vertical':
            width = max(10, n_genes * 0.3)
            height = 4 + (1.5 if effective_show_brackets else 0) + (1.0 if effective_legend_loc == 'bottom' and legend_mode != 'none' else 0)
        else:
            height = max(6, n_genes * 0.25)
            width = 8 + (2.0 if effective_show_brackets else 0) + (2.0 if effective_legend_loc == 'right' and legend_mode != 'none' else 0)
        figsize = (width, height)
    
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()
    
    # Phase 10: Grid control (v4.3 fix: avoid warning when show_grid=False)
    if show_grid:
        ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
        ax.set_axisbelow(True)
    else:
        ax.grid(False)
    
    # Drawing
    if orientation == 'vertical':
        _draw_vertical_stream(ax, x_plot, y_plot, score_columns, colors, plot_type)
    else:
        if discrete and not discrete_rounded:
            _draw_horizontal_bars(ax, normalized, score_columns, colors, n_genes)
        else:
            _draw_horizontal_fill(ax, x_plot, y_plot, score_columns, colors)
    
    # Boundaries
    if show_boundaries and len(celltype_boundaries) > 2:
        for boundary in celltype_boundaries[1:-1]:
            if orientation == 'vertical':
                ax.axvline(x=boundary - 0.5, color=boundary_color, lw=boundary_linewidth, ls='--', alpha=0.7, zorder=10)
            else:
                ax.axhline(y=boundary - 0.5, color=boundary_color, lw=boundary_linewidth, ls='--', alpha=0.7, zorder=10)
    
    # Axis limits
    if orientation == 'vertical':
        ax.set_xlim(-0.5, n_genes - 0.5)
        y_min, y_max = ax.get_ylim()
    else:
        ax.set_ylim(-0.5, n_genes - 0.5)
        ax.invert_yaxis()
        # v4.3 fix: Force x-axis to start at 0 (no gap on left)
        _, x_max_auto = ax.get_xlim()
        ax.set_xlim(0, x_max_auto)
        x_min, x_max = ax.get_xlim()
    
    # Brackets
    if effective_show_brackets and len(celltype_boundaries) > 1:
        if orientation == 'vertical':
            bracket_space = _draw_vertical_brackets(
                ax, ordering, celltype_boundaries, effective_bracket_position,
                bracket_height, bracket_text_rotation, bracket_linewidth, celltype_fontsize, y_min, y_max
            )
            if effective_bracket_position == 'top':
                ax.set_ylim(y_min, y_max + bracket_space)
            else:
                ax.set_ylim(y_min - bracket_space, y_max)
        else:
            bracket_space = _draw_horizontal_brackets(
                ax, ordering, celltype_boundaries, effective_bracket_position,
                bracket_height, bracket_text_rotation, bracket_linewidth, celltype_fontsize, x_min, x_max
            )
            if effective_bracket_position == 'right':
                ax.set_xlim(0, x_max + bracket_space)  # v4.3: Keep x starting at 0
            else:
                ax.set_xlim(x_min - bracket_space, x_max)
    
    # v4.3 fix: Ensure horizontal x-axis starts at 0 after all adjustments
    if orientation == 'horizontal':
        x_min_final, x_max_final = ax.get_xlim()
        if x_min_final != 0:
            ax.set_xlim(0, x_max_final)
    
    # Gene labels
    if show_gene_labels and gene_label_mode != 'none':
        if orientation == 'vertical':
            _draw_gene_labels_with_offset(
                ax, selected_genes, n_genes, gene_label_offset, gene_label_mode,
                celltype_boundaries, effective_gene_position, gene_tick_rotation,
                gene_fontsize, gene_label_ha, gene_label_va
            )
            if effective_gene_position == 'top':
                ax.xaxis.set_ticks_position('top')
                ax.tick_params(axis='x', top=True, bottom=False)
        else:
            if gene_label_mode == 'all':
                ax.set_yticks(list(range(n_genes)))
                ax.set_yticklabels(selected_genes, fontsize=gene_fontsize)
            elif gene_label_mode == 'first':
                tick_pos = [b for b in celltype_boundaries[:-1] if b < n_genes]
                ax.set_yticks(tick_pos)
                ax.set_yticklabels([selected_genes[p] for p in tick_pos], fontsize=gene_fontsize)
    else:
        if orientation == 'vertical':
            ax.set_xticks([])
        else:
            ax.set_yticks([])
    
    # Axis formatting
    if orientation == 'vertical':
        ax.set_ylabel('Specificity Score' if plot_type == 'stacked' else 'Relative Specificity', fontsize=10)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if effective_gene_position == 'top':
            ax.spines['top'].set_visible(True)
            ax.spines['bottom'].set_visible(False)
    else:
        ax.set_xlabel('Specificity Score' if plot_type == 'stacked' else 'Relative Specificity', fontsize=10)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    if title:
        ax.set_title(title, fontsize=12, fontweight='bold', pad=20)
    
    # Legend
    _setup_legend(ax, fig, score_columns, colors, legend_mode, effective_legend_loc, legend_ncol,
                  legend_offset, legend_fontsize, effective_gene_position, gene_tick_rotation,
                  gene_fontsize, orientation)
    
    # Layout
    plt.tight_layout()
    if orientation == 'vertical':
        if effective_legend_loc == 'bottom' and legend_mode != 'none':
            total_offset = _calculate_legend_offset(gene_tick_rotation, gene_fontsize, effective_gene_position, legend_offset)
            plt.subplots_adjust(bottom=min(0.15 + total_offset, 0.4))
        elif effective_legend_loc == 'right' and legend_mode != 'none':
            plt.subplots_adjust(right=0.85)
    else:
        if effective_legend_loc == 'right' and legend_mode != 'none':
            plt.subplots_adjust(right=0.82)
        elif effective_legend_loc == 'bottom' and legend_mode != 'none':
            plt.subplots_adjust(bottom=0.2)
    
    # Save
    if save:
        fig.savefig(save, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {save}")
    
    # Return data
    processed_df = None
    if return_data:
        processed_df = pd.DataFrame({
            'gene': selected_genes,
            'group': gene_celltype_map,
            'x_position': gene_positions,
        })
        for i, ct in enumerate(score_columns):
            processed_df[f'score_{ct}'] = score_matrix[:, i]
        if plot_type != 'stacked':
            for i, ct in enumerate(score_columns):
                processed_df[f'normalized_{ct}'] = normalized[:, i]
    
    # Show/Return
    if show:
        plt.show()
    
    if return_data and return_fig:
        return fig, processed_df
    elif return_fig:
        return fig
    elif return_data:
        if not show:
            plt.close(fig)
        return processed_df
    else:
        if not show:
            plt.close(fig)
        return None