from anndata import AnnData
import scanpy as sc
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

import numpy as np
from scipy.sparse import csr_matrix
from sklearn.metrics.pairwise import cosine_similarity


def cosg(
    adata,
    groupby = 'CellTypes',
    groups: Union[Literal['all'], Iterable[str]] = 'all',
    mu:float=1,
    remove_lowly_expressed: bool= False,
    expressed_pct: Optional[float] = 0.1,
    n_genes_user: int = 50,
    key_added: Optional[str] = None,
    calculate_logfoldchanges: bool = False,
    use_raw: bool = True,
    layer: Optional[str] = None,
    reference: str = 'rest',
    return_by_group : bool = True,
    verbosity: int = 0,
    copy: bool = False
):
    """\
    Marker gene identification for single-cell sequencing data using COSG.
    
    Parameters
    ----------
    adata
        Annotated data matrix. Note: input paramters are simliar to the parameters used for scanpy's rank_genes_groups() function.
    groupby
        The key of the cell groups in .obs, the default value is set to 'CellTypes'.
    groups
        Subset of cell groups, e.g. [`'g1'`, `'g2'`, `'g3'`], to which comparison shall be restricted. The default value is 'all', and all groups will be compared.
    mu
        The penalty restricting marker genes expressing in non-target cell groups. Larger value represents more strict restrictions. mu should be >= 0, and by default, mu = 1.
    remove_lowly_expressed
        If True, genes that express a percentage of target cells smaller than a specific value (`expressed_pct`) are not considered as marker genes for the target cells. The default value is False.
    expressed_pct
        When `remove_lowly_expressed` is set to True, genes that express a percentage of target cells smaller than a specific value (`expressed_pct`) are not considered as marker genes for the target cells. The default value for `expressed_pct`
        is 0.1 (10%).
    n_genes_user
        The number of genes that appear in the returned tables. The default value is 50.
    key_added
        The key in `adata.uns` information is saved to.
    calculate_logfoldchanges
        Calculate logfoldchanges.
    use_raw
        Use `raw` attribute of `adata` if present.
    layer
        Key from `adata.layers` whose value will be used to perform tests on.
    reference
        If `'rest'`, compare each group to the union of the rest of the group.
        If a group identifier, compare with respect to this group.
    return_by_group
        Whether return the COSG scores summarized by each group. This will output another extra copy of the results. Defaults to True.
    verbosity
        Controls the verbosity of logging messages, defaults to 0.
    copy
        If True, returns a copy of `adata` with the computed embeddings instead of modifying in place. Defaults to False.
 
    Returns
    -------
    names : structured `np.ndarray` (`.uns['rank_genes_groups']`)
        Structured array to be indexed by group id storing the gene names. Ordered according to scores.
    scores : structured `np.ndarray` (`.uns['rank_genes_groups']`)
        Structured array to be indexed by group id storing COSG scores for each gene for each
        group. Ordered according to scores.
    
    Examples
    --------
    >>> import cosg as cosg
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> cosg.cosg(adata, key_added='cosg', groupby='bulk_labels')
    >>> sc.pl.rank_genes_groups(adata, key='cosg')
    """
    
    adata = adata.copy() if copy else adata
        
    if layer is not None:
        if use_raw:
            raise ValueError("Cannot specify `layer` and have `use_raw = True`.")
        cellxgene = adata.layers[layer]
    else:
        if use_raw and adata.raw is not None:
             cellxgene = adata.raw.X
        cellxgene = adata.X
    
    
    ### Refer to scanpy's framework
    # https://github.com/theislab/scanpy/blob/5533b644e796379fd146bf8e659fd49f92f718cd/scanpy/tools/_rank_genes_groups.py#L559
    if key_added is None:
        key_added = 'rank_genes_groups'
    adata.uns[key_added] = {}
    adata.uns[key_added]['params'] = dict(
        groupby=groupby,
        reference=reference,
        groups=groups,
        method='COSG',
        use_raw=use_raw,
        layer=layer,
        mu=mu,
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
        

    
    groups_order=np.unique(group_info)
    n_cluster=len(groups_order)
    
    n_cell=cellxgene.shape[0]
    
    ### Efficiently create a sparse matrix for the cluster_mat matrix
    group_to_row = {group: i for i, group in enumerate(groups_order)}
    row_indices = np.array([group_to_row[group] for group in group_info])
    col_indices = np.arange(n_cell)
    data = np.ones_like(col_indices, dtype=int)
    cluster_mat = csr_matrix((data, (row_indices, col_indices)), shape=(n_cluster, n_cell))
    
    if sparse.issparse(cellxgene):
        ### the dimension is: Gene x lambda
        cosine_sim=cosine_similarity(X=cellxgene.T,Y=cluster_mat,dense_output=False) 
        
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
        genexlambda.data=np.multiply(genexlambda.data,cosine_sim.data)
    
    ### If the cellxgene is not a sparse matrix
    else:
        ### Convert to dense matrix
        cluster_mat=cluster_mat.toarray()
        
        ## Not using sparse matrix    
        cosine_sim=cosine_similarity(X=cellxgene.T, Y=cluster_mat, dense_output=True) 

        pos_nonzero=cosine_sim!=0
        e_power2=np.multiply(cosine_sim,cosine_sim)
        e_power2_sum=np.sum(e_power2,axis=1)
        e_power2[pos_nonzero]=np.true_divide(e_power2[pos_nonzero],(1-mu)*e_power2[pos_nonzero]+mu*(np.dot(e_power2_sum.reshape(e_power2_sum.shape[0],1),np.repeat(1,e_power2.shape[1]).reshape(1,e_power2.shape[1]))[pos_nonzero]))
        e_power2[pos_nonzero]=np.multiply(e_power2[pos_nonzero],cosine_sim[pos_nonzero])
        genexlambda=e_power2

    ### Refer to scanpy
    rank_stats=None
    
    ### Whether to calculate logfoldchanges, because this is required in scanpy 1.8
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
    if sparse.issparse(cellxgene):
        cellxgene.eliminate_zeros()
    if sparse.issparse(cellxgene):
        get_nonzeros = lambda X: X.getnnz(axis=0)
    else:
        get_nonzeros = lambda X: np.count_nonzero(X, axis=0)
    
    
    order_i=0
    for group_i in groups_order:    
        idx_i=group_info==group_i 

        ### Convert to numpy array
        idx_i=idx_i.values

        ## Compare the most ideal case to the worst case
        if sparse.issparse(cellxgene):
            scores=genexlambda[:,order_i].toarray()[:,0] ###Changed .A to .toarray() for future support
        else:
            scores=genexlambda[:,order_i]

        
        ### Mask these genes expressed in less than 3 cells in the cluster of interest
        if remove_lowly_expressed:
            n_cells_expressed=get_nonzeros(cellxgene[idx_i])
            n_cells_i=np.sum(idx_i)
            scores[n_cells_expressed<n_cells_i*expressed_pct]= -1

        global_indices = _select_top_n(scores, n_genes_user)

        # if rank_stats is None:
        #     idx = pd.MultiIndex.from_tuples([(group_i,'names')])
        #     rank_stats = pd.DataFrame(columns=idx)
        # rank_stats[group_i, 'names'] = adata.var_names[global_indices]
        # rank_stats[group_i, 'scores'] = scores[global_indices]
            
            
        # Initialize the DataFrame if rank_stats is None
        if rank_stats is None:
            rank_stats = pd.DataFrame()

        # Prepare data for new columns
        columns_data = {
            (group_i, 'names'): adata.var_names.values[global_indices],
            (group_i, 'scores'): scores[global_indices]
        }
        
        if calculate_logfoldchanges:
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
                )  # Avoid division by zero
                columns_data[(group_i, 'logfoldchanges')] = np.log2(foldchanges[global_indices])


        # Create a new DataFrame for the current group
        new_data = pd.DataFrame(columns_data)

        # Concatenate with rank_stats
        rank_stats = pd.concat([rank_stats, new_data], axis=1)
    
        order_i=order_i+1
    
    #### also return a copy of the results showing the results by group
    if return_by_group:
        adata.uns[key_added]['COSG']=rank_stats
    
    
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


import pandas as pd

def indexByGene(
    df: pd.DataFrame,
    gene_key: str = "names",
    score_key: str = "scores",
    set_nan_to_zero: bool= False,
    convert_negative_one_to_zero: bool = True
):
    """
    Reshapes a DataFrame with MultiIndex columns where gene names are under the specified key 
    and scores are under the corresponding key. The resulting DataFrame will have gene names as the index 
    and scores for each cell type in separate columns.

    Note
    ----
    This function is designed for reindexing COSG's output stored in `adata.uns['cosg']['COSG']`. 
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
    set_nan_to_zero : bool, optional, default=False
        If True, replaces all NaN values in the resulting DataFrame with 0.
        Set this parameter to True if `n_genes_user=adata.n_vars` is not set when calling the `cosg.cosg` function.
    convert_negative_one_to_zero : bool, optional, default=True
        If True, replaces all occurrences of -1.0 in the resulting DataFrame with 0.


    Returns
    -------
    pd.DataFrame
        Reshaped DataFrame with genes as index and scores per cell type.
        
    
    Raises
    ------
    TypeError
        If the input df is not a pandas DataFrame.
    ValueError
        If the input DataFrame does not have MultiIndex columns with at least two levels,
        or if the specified gene_key or score_key is not found in the first level of the columns.
        
    Example
    -------
    >>> cosg_df = indexByGene(
    ...     adata.uns['cosg']['COSG'],
    ...     gene_key="names",
    ...     score_key="scores",
    ...     set_nan_to_zero=True,
    ...     convert_negative_one_to_zero=True
    ... )
    >>>  ### Check the reindexed data frame
    >>> cosg_df.head()
    """
    
    # Validate input type
    if not isinstance(df, pd.DataFrame):
        raise TypeError("The input df must be a pandas DataFrame.")

    # Validate that columns are a MultiIndex with at least two levels
    if not isinstance(df.columns, pd.MultiIndex) or len(df.columns.levels) < 2:
        raise ValueError("The input DataFrame must have MultiIndex columns with at least two levels.")

    # Validate that gene_key and score_key are present in the first level of the MultiIndex
    if gene_key not in df.columns.get_level_values(0):
        raise ValueError(f"The gene_key '{gene_key}' is not present in the first level of the DataFrame columns.")
    if score_key not in df.columns.get_level_values(0):
        raise ValueError(f"The score_key '{score_key}' is not present in the first level of the DataFrame columns.")
    
    # Extract unique cell types from the second level of the MultiIndex
    cell_types = df.columns.get_level_values(1).unique()

    # Get all unique gene names across cell types using the gene_key column
    all_genes = pd.concat([df[(gene_key, ct)] for ct in cell_types]).unique()
    
    df_scores = pd.DataFrame(index=all_genes)

    # Assign scores to the corresponding gene names
    for ct in cell_types:
        scores_series = df.set_index((gene_key, ct))[(score_key, ct)]
        # Reindex the scores_series to align with all_genes
        df_scores[ct] = scores_series.reindex(df_scores.index)
    

    # Optionally replace NaN values with 0 if set_nan_to_zero is True
    if set_nan_to_zero and df_scores.isnull().values.any():
        df_scores.fillna(0, inplace=True)
        
    # Optionally replace -1.0 values with 0
    if convert_negative_one_to_zero:
        df_scores = df_scores.replace(-1.0, 0)

    return df_scores


import numpy as np
import pandas as pd

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
from scipy.sparse import issparse

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
            logg.debug(
                f'{np.array(groups_order_subset)} invalid! specify valid '
                f'groups_order (or indices) from {adata.obs[key].cat.categories}',
            )
            from sys import exit

            exit(0)
        groups_masks = groups_masks[groups_ids]
        groups_order_subset = adata.obs[key].cat.categories[groups_ids].values
    else:
        groups_order_subset = groups_order.values
    return groups_order_subset, groups_masks


### Import from Scanpy   
import numpy as np
from scipy import sparse
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

        # for correct getnnz calculation
        if issparse(X):
            X.eliminate_zeros()

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

        if issparse(self.X):
            get_nonzeros = lambda X: X.getnnz(axis=0)
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

