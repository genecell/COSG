|Stars| |PyPI| |Docs| |Total downloads| |Monthly downloads|

.. |Stars| image:: https://img.shields.io/github/stars/genecell/COSG?logo=GitHub&color=yellow
   :target: https://github.com/genecell/COSG/stargazers
.. |PyPI| image:: https://img.shields.io/pypi/v/cosg?logo=PyPI
   :target: https://pypi.org/project/cosg
.. |Docs| image:: https://readthedocs.org/projects/cosg/badge/?version=latest
   :target: https://cosg.readthedocs.io
.. |Total downloads| image:: https://static.pepy.tech/personalized-badge/cosg?period=total&units=international_system&left_color=black&right_color=orange&left_text=downloads
   :target: https://pepy.tech/project/cosg
.. |Monthly downloads| image:: https://static.pepy.tech/personalized-badge/cosg?period=month&units=international_system&left_color=black&right_color=orange&left_text=downloads/month
 :target: https://pepy.tech/project/cosg

Accurate and fast cell marker gene identification with COSG
=======================================================================================================

Overview
---------
COSG is a cosine similarity-based method for more accurate and scalable marker gene identification.

- COSG is a general method for cell marker gene identification across different data modalities, e.g., scRNA-seq, scATAC-seq and spatially resolved transcriptome data.
- Marker genes or genomic regions identified by COSG are more indicative and with greater cell-type specificity.
- COSG is ultrafast for large-scale datasets, and is capable of identifying marker genes for one million cells in less than two minutes.

The method and benchmarking results are described in `Dai et al., (2022)`_. 

Documentation
--------------
The documentation for COSG is available `here <https://cosg.readthedocs.io/en/latest/>`_.

Tutorial
---------

The `COSG tutorial <https://nbviewer.jupyter.org/github/genecell/COSG/blob/main/tutorials/COSG-tutorial.ipynb>`_ provides a quick-start guide for using COSG and demonstrates the superior performance of COSG as compared with other methods, and the `Jupyter notebook <https://github.com/genecell/COSG/blob/main/tutorials/COSG-tutorial.ipynb>`_ is also available.

Question
---------
For questions about the code and tutorial, please contact Min Dai, dai@broadinstitute.org.

Example
---------
Run COSG:

.. code-block:: python
   
   import cosg
   n_gene=30
   groupby='CellTypes'
   cosg.cosg(adata,
       key_added='cosg',
       # use_raw=False, layer='log1p', ## e.g., if you want to use the log1p layer in adata
       mu=100,
       expressed_pct=0.1,
       remove_lowly_expressed=True,
        n_genes_user=100,
                  groupby=groupby)

Draw the dot plot:

.. code-block:: python
   
   sc.tl.dendrogram(adata,groupby=groupby)
   df_tmp=pd.DataFrame(adata.uns['cosg']['names'][:3,]).T
   df_tmp.reindex(adata.obs[groupby].cat.categories)
   marker_genes_list=np.ravel(df_tmp.reindex(adata.uns['dendrogram_'+groupby]['categories_ordered']))
   sc.pl.dotplot(adata, marker_genes_list,
                groupby=groupby,              
                dendrogram=True,
                 swap_axes=False,
                standard_scale='var',
                cmap='Spectral_r')


Output the marker list as pandas dataframe:

.. code-block:: python
   
   marker_gene=pd.DataFrame(adata.uns['cosg']['names'])
   marker_gene.head()

You could also check the COSG scores:

.. code-block:: python
   
   marker_gene_scores=pd.DataFrame(adata.uns['cosg']['scores'])
   marker_gene_scores.head()


Citation
---------
If COSG is useful for your research, please consider citing `Dai et al., (2022)`_.

.. _Dai et al., (2022): https://academic.oup.com/bib/advance-article-abstract/doi/10.1093/bib/bbab579/6511197?redirectedFrom=fulltext


