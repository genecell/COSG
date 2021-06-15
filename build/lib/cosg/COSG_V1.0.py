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

def _select_top_n(scores, n_top):
    reference_indices = np.arange(scores.shape[0], dtype=int)
    partition = np.argpartition(scores, -n_top)[-n_top:]
    partial_indices = np.argsort(scores[partition])[::-1]
    global_indices = reference_indices[partition][partial_indices]
    return global_indices
        
def cosg(
    adata,
    groupby='CellTypes',
    groups: Union[Literal['all'], Iterable[str]] = 'all',

    mu=1,
    remove_lowly_expressed:bool=False,
    expressed_pct:Optional[float] = 0.1,

    n_genes_user:int =50,
    key_added: Optional[str] = None,
    use_raw: bool = True,
    layer: Optional[str] = None,
    reference: str = 'rest',    

    copy:bool=False
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
    use_raw
        Use `raw` attribute of `adata` if present.
    layer
        Key from `adata.layers` whose value will be used to perform tests on.
    reference
        If `'rest'`, compare each group to the union of the rest of the group.
        If a group identifier, compare with respect to this group.
 
    Returns
    -------
    **names** : structured `np.ndarray` (`.uns['rank_genes_groups']`)
        Structured array to be indexed by group id storing the gene names. Ordered according to scores.
    **scores** : structured `np.ndarray` (`.uns['rank_genes_groups']`)
        Structured array to be indexed by group id storing COSG scores for each gene for each
        group. Ordered according to scores.
    Notes
    -----
    Contact: daimin@zju.edu.cn
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
            raise ValueError("Cannot specify `layer` and have `use_raw=True`.")
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
    method='cosg',
    use_raw=use_raw,
    layer=layer,
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
    cluster_mat=np.zeros(shape=(n_cluster,n_cell))

    ### To further restrict expression in other clusters, can think about a better way, such as taking the cluster similarities into consideration
    order_i=0
    for group_i in groups_order:    
        idx_i=group_info==group_i 
        cluster_mat[order_i,:][idx_i]=1
        order_i=order_i+1
    
    if sparse.issparse(cellxgene):
        ### Convert to sparse matrix
        from scipy.sparse import csr_matrix
        cluster_mat=csr_matrix(cluster_mat)

        from sklearn.metrics.pairwise import cosine_similarity

        ### the dimension is: Gene x lambda
        cosine_sim=cosine_similarity(X=cellxgene.T,Y=cluster_mat,dense_output=False) 

        pos_nonzero=cosine_sim.nonzero()
        genexlambda=cosine_sim.multiply(cosine_sim)

        e_power2_sum=genexlambda.sum(axis=1)


        if mu==1:
            genexlambda[pos_nonzero]=genexlambda[pos_nonzero]/(np.repeat(e_power2_sum,genexlambda.shape[1],axis=1)[pos_nonzero])
                                                    
        else:
            genexlambda[pos_nonzero]=genexlambda[pos_nonzero]/((1-mu)*genexlambda[pos_nonzero]+
                                                     mu*(
                                                        np.repeat(e_power2_sum,genexlambda.shape[1],axis=1)[pos_nonzero])
                                                    )
        genexlambda[pos_nonzero]=np.multiply(genexlambda[pos_nonzero],cosine_sim[pos_nonzero])
    
    ### If the cellxgene is not a sparse matrix
    else:
         ## Not using sparse matrix
        from sklearn.metrics.pairwise import cosine_similarity    
        cosine_sim=cosine_similarity(X=cellxgene.T,Y=cluster_mat,dense_output=True) 

        pos_nonzero=cosine_sim!=0
        e_power2=np.multiply(cosine_sim,cosine_sim)
        e_power2_sum=np.sum(e_power2,axis=1)
        e_power2[pos_nonzero]=np.true_divide(e_power2[pos_nonzero],(1-mu)*e_power2[pos_nonzero]+mu*(np.dot(e_power2_sum.reshape(e_power2_sum.shape[0],1),np.repeat(1,e_power2.shape[1]).reshape(1,e_power2.shape[1]))[pos_nonzero]))
        e_power2[pos_nonzero]=np.multiply(e_power2[pos_nonzero],cosine_sim[pos_nonzero])
        genexlambda=e_power2

    ### Refer to scanpy
    rank_stats=None
    
    
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
            scores=genexlambda[:,order_i].A[:,0]
        else:
            scores=genexlambda[:,order_i]

        
        ### Mask these genes expressed in less than 3 cells in the cluster of interest
        if remove_lowly_expressed:
            n_cells_expressed=get_nonzeros(cellxgene[idx_i])
            n_cells_i=np.sum(idx_i)
            scores[n_cells_expressed<n_cells_i*expressed_pct]= -1

        global_indices = _select_top_n(scores, n_genes_user)

        if rank_stats is None:
            idx = pd.MultiIndex.from_tuples([(group_i,'names')])
            rank_stats = pd.DataFrame(columns=idx)
        rank_stats[group_i, 'names'] = adata.var_names[global_indices]
        rank_stats[group_i, 'scores'] = scores[global_indices]

        order_i=order_i+1
            
    ## Refer to scanpy
    dtypes = {
    'names': 'O',
    'scores': 'float32',
    # 'logfoldchanges': 'float32',
    # 'pvals': 'float64',
    # 'pvals_adj': 'float64',
    }
    ###
    rank_stats.columns = rank_stats.columns.swaplevel()
    for col in rank_stats.columns.levels[0]:
        adata.uns[key_added][col]=rank_stats[col].to_records(
    index=False, column_dtypes=dtypes[col]
    )
        
    print('**finished identifying marker genes by COSG**')
        
    ### Return the result
    return adata if copy else None