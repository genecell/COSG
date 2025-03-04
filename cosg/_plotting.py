import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import matplotlib.colors as mcolors
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
import matplotlib.patheffects as PathEffects
from scipy.sparse import issparse

from .cosg import indexByGene, iqrLogNormalize


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
def plotMarkerDendrogram(
    adata,
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
        If provided, this color is used for gene labels; otherwise, the gene node’s colormap color is used.
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
        Specification is as matplotlib.scatter marker, one of ‘so^>v<dph8’. In detail:
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
        cell_types = list(cosg_score_df.columns)
    else:
        rep = adata.obsm[use_rep]
        df_rep = pd.DataFrame(rep, index=adata.obs_names)
        df_rep[group_by] = adata.obs[group_by].values
        group_means = df_rep.groupby(group_by, observed=True).mean()
        cell_types = list(group_means.index)
        data = group_means.values
        D = pdist(data, metric=distance_metric)
        Z = linkage(D, method=linkage_method)
    
    N = len(cell_types)
    
    
    
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
    
    
    
    # Extract top N marker genes for all cell types at once
    marker_genes_df = adata.uns[cosg_key]['COSG']['names'].iloc[:top_n_genes]  # Slice once for efficiency
    selected_genes = marker_genes_df.values.flatten()  # Flatten to get all genes as a 1D list
    selected_genes = pd.Index(selected_genes).dropna().unique()  # Remove NaNs & duplicates


    # Attach top marker gene nodes to each cell type leaf
    gene_nodes = {}
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
    filtered_cosg_score_df = cosg_score_df.loc[selected_genes]  # Keep only top marker genes
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
        if ntype in ['internal', 'root', 'cluster']:
            node_sizes[n] = 50
        elif ntype == 'cell_type':
            node_sizes[n] = 600
    for n, d in G.nodes(data=True):
        if d.get('node_type') == 'gene':
            node_sizes[n] = gene_node_sizes.get(n, 300)
    
    # Compute radial layout
    pos = _radial_dendrogram_layout(G, root, radius_step=radius_step)
    
    # Set palette for cell type nodes if not provided
    if palette is None and f"{group_by}_colors" in adata.uns:
        palette = adata.uns[f"{group_by}_colors"]
    
    # Drawing
    fig = plt.figure(figsize=figure_size)
    # Create a square main Axes for the network plot
    ax_main = fig.add_axes([0.1, 0.1, 0.75, 0.75])
    ax_main.set_aspect('equal')
    # Create separate Axes for the legends on the right side, side by side.
    # Create dot size legend Axes in the UPPER RIGHT corner
    ax_ds = fig.add_axes([0.8, 0.6, 0.05, 0.2])  # (left, bottom, width, height)
    # Create colorbar Axes in the LOWER RIGHT corner
    ax_cb = fig.add_axes([0.85, 0.1, colorbar_width, 0.2])  # (left, bottom, width, height)

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
            node_colors[n] = 'gray'
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
    nx.draw_networkx_edges(G, pos, ax=ax_main, arrows=False)
    

    
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
    min_expr = gene_expr_percentage[gene_expr_percentage > 0].min().min()  # Ignore 0% values
    max_expr = gene_expr_percentage.max().max()  # Maximum expression percentage

    # Round min/max to nearest multiple of 5
    min_expr_rounded = np.floor(min_expr / 5) * 5
    max_expr_rounded = np.ceil(max_expr / 5) * 5

    # Generate up to 5 evenly spaced values, ignoring 0%
    num_circles = min(5, int((max_expr_rounded - min_expr_rounded) / 5) + 1)
    legend_percentages = np.linspace(min_expr_rounded, max_expr_rounded, num=num_circles)
    legend_percentages = np.unique(np.round(legend_percentages / 5) * 5).astype(int)  # Ensure multiples of 5

    # Generate dot size legend (skip 0%)
    legend_markers = [
        plt.Line2D([0], [0], marker=node_shape_gene, color='black', label=f'{p}%', 
                   markerfacecolor='white', markersize=np.sqrt(gene_size_scale * (p/100)))
        for p in legend_percentages if p > 0
    ]

    
    
    # Place the dot size legend in the upper-left corner of ax_ds.
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
   
    
