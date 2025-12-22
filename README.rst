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
=============================================================

Overview
---------

COSG is a cosine similarity-based method for more accurate and scalable marker gene identification.

- COSG is a general method for cell marker gene identification across different data modalities, e.g., scRNA-seq, scATAC-seq, and spatially resolved transcriptome data.

- Marker genes or genomic regions identified by COSG are more indicative and with greater cell-type specificity.

- COSG is ultrafast for large-scale datasets and is capable of identifying marker genes for one million cells in less than two minutes.

The method and benchmarking results are described in `Dai et al. (2022)`_.

Additionally, the R version of COSG is available `here <https://github.com/genecell/COSGR>`_.

Note I: we released our Python toolkit, `PIASO <https://piaso.org>`_, in which some methods were built upon COSG.

Note II: we have also recently released `PIASOmarkerDB <https://piaso.org/piasomarkerdb>`_ for beta testing! Besides cataloging comprehensive cell type-specific marker genes, PIASOmarkerDB also supports cell type inference based on your input gene list. The gene list could be the top marker genes identified by COSG in your cluster of interest, the differentially expressed genes in either snRNA-seq data or bulk RNA-seq, genes belonging to specific pathway or Gene Ontology term, or any gene set you are interested in!

Documentation
--------------

`COSG documentation <https://genecell.github.io/COSG/>`_.


Release notes
-------------

**Release v1.0.3** (March 11, 2025)


- Fixed the incompatibility with multiple index columns of ``adata.uns['cosg']['COSG']`` in ``adata.write`` function

- Enhanced ``plotMarkerDendrogram`` function with several new capabilities:

  - Implemented support for customized cell type-gene pairs
  - Added color control for nodes and edges
  - Added cell type filtering functionality
  - Integrated support for curved edges in visualization


**Release v1.0.2** (March 5, 2025)


- Added ``plotMarkerDotplot`` and ``plotMarkerDendrogram`` for enhanced marker gene visualization. 

- Introduced support for ``batch_key`` to compute cosine similarities separately across different batches.  

- Enabled calculation of normalized COSG scores for comparing gene expression specificity across cell types or datasets.  

- Resolved a SciPy version deprecation issue related to ``.A`` attribute usage.  

- Fixed a DataFrame manipulation warning.  

- Added verbosity control, allowing users to adjust log output levels.  

**Release v1.0.1** (June 15, 2021)


- First release in PyPI. 

Installation
------------
Stable version:

.. code-block:: bash

   pip install cosg

Development version:

.. code-block:: bash

   pip install git+https://github.com/genecell/COSG.git


Example
---------
Run COSG:

.. code-block:: python
   
   import cosg
   n_genes=30
   groupby='CellTypes'
   cosg.cosg(
      adata,
      key_added='cosg',
      # use_raw=False, layer='log1p', ## e.g., if you want to use the log1p layer in adata
      mu=100,
      expressed_pct=0.1,
      remove_lowly_expressed=True,
      n_genes_user=n_genes,
      groupby=groupby
   )

Draw the dot plot:

.. code-block:: python
   
   cosg.plotMarkerDotplot(
       adata,
       groupby=groupby,
       top_n_genes=3,
       key_cosg='cosg',
       use_rep='X_pca', ## Change use_rep to the cell embeddings key you'd like to use
       swap_axes=False,
       standard_scale='var',
       cmap='Spectral_r',
       # save='test.pdf'
   )



Output the marker list as pandas dataframe:

.. code-block:: python
   
   marker_gene=pd.DataFrame(adata.uns['cosg']['names'])
   marker_gene.head()

You could also check the COSG scores:

.. code-block:: python
   
   marker_gene_scores=pd.DataFrame(adata.uns['cosg']['scores'])
   marker_gene_scores.head()


Question
---------
For questions about the code and tutorial, please contact Min Dai, dai@broadinstitute.org.


Citation
---------
If COSG is useful for your research, please consider citing `Dai et al. (2022)`_.

.. _Dai et al. (2022): https://academic.oup.com/bib/advance-article-abstract/doi/10.1093/bib/bbab579/6511197?redirectedFrom=fulltext
