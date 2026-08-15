|Stars| |PyPI| |Bioconda| |Docs| |Total downloads| |Monthly downloads|

.. |Stars| image:: https://img.shields.io/github/stars/genecell/COSG?logo=GitHub&color=yellow
   :target: https://github.com/genecell/COSG/stargazers
.. |PyPI| image:: https://img.shields.io/pypi/v/cosg?logo=PyPI
   :target: https://pypi.org/project/cosg
.. |Bioconda| image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat
   :target: http://bioconda.github.io/recipes/cosg/README.html
   :alt: install with bioconda
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

Note II: we have also recently released `PIASOmarkerDB <https://piaso.org/piasomarkerdb>`_ for beta testing.

Note III: COSG is also available for online analysis via `Galaxy platform <https://usegalaxy.eu/root?tool_id=cosg>`_.

Documentation
--------------

`COSG documentation <https://genecell.github.io/COSG/>`_.


Installation
------------
Stable version (PyPI):

.. code-block:: bash

   pip install cosg

Stable version (bioconda):

.. code-block:: bash

   conda install -c conda-forge -c bioconda cosg

Development version:

.. code-block:: bash

   pip install git+https://github.com/genecell/COSG.git



Release notes
-------------

**Release v1.1.0** (August 12, 2026)


- Added a **streaming backend for** ``.cytome`` **datasets**: marker detection reads the file in chunks, so peak memory does not scale with the number of cells.

- The cytome streaming backend is an **optional extra**, on the same reasoning as scanpy above: ``pip install 'cosg[cytome]'``. Calling ``cosg.cosg()`` on a ``.cytome`` path without it raises an error naming the extra rather than failing obscurely.

- **The default streaming path no longer needs PIASO.** ``layer='log1p'`` (the RNA default) is computed by COSG itself and gives identical numbers to the implementation it replaced. Only ``layer='infog'`` and ``layer='tfidf'`` need PIASO, since those are PIASO normalizations; the error naming it now gives the correct install command, ``pip install piaso-tools``.

- Added a **GPU path** (CuPy), selected with ``device='cpu' | 'gpu' | 'auto'`` on both the in-memory and streaming paths.

- ``cosg.cosg()`` is now a **single polymorphic entry point**. It dispatches on its first argument:

  - an ``AnnData`` takes the in-memory path and writes ``adata.uns[key_added]``
  - a ``str`` or ``pathlib.Path`` to a ``.cytome`` file takes the streaming path and returns a dict
  - an **open** ``CytomeDataset`` (from ``cytome.open()``) also takes the streaming path; the caller keeps ownership of it and it is **not** closed

- Extended ``batch_key`` to the streaming, GPU and feature-batched paths.

- Added ``output_format`` for the streaming path: ``'dict'``, ``'long'`` or ``'dense'``.

- Added ``cosg.__version__``, which previously raised ``AttributeError``.

- **scanpy is now an optional extra rather than a required dependency.** Only ``plotMarkerDotplot`` uses it, so a plain ``pip install cosg`` is ~60 MB and 9 packages lighter. Install ``pip install 'cosg[dotplot]'`` to keep that function; calling it without scanpy raises an error naming the extra. Cell-type ordering that previously used ``scanpy.tl.dendrogram`` is now computed internally and reproduces it exactly.

- The streaming default layer is resolved from the modality: ``RNA``/``GA`` to ``log1p``, ``ATAC``/``tiles`` to ``tfidf``.

- **Behaviour change:** ``remove_lowly_expressed`` now defaults to ``True`` on the AnnData path, matching the streaming variant. Pass ``remove_lowly_expressed=False`` to restore the previous behaviour.

- **Behaviour change:** IQR normalisation is computed over all values per group, matching ``iqrLogNormalize``.


**Release v1.0.4** (March 5, 2026)


- Added ``plotMarkerStream`` for visualising marker gene specificity as a streamgraph.

- Added ``expressed_min_num_cells_in_target_group`` (default 3), which floors the expression threshold at ``max(n_cells * expressed_pct, 3)`` so small clusters do not get an overly permissive cutoff.

- Added input validation for ``groupby``, ``groups``, ``groups`` combined with ``batch_key``, and ``n_genes_user``.


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

Execution modes
---------------

COSG runs the same algorithm in several modes, selected by the available
hardware:

+------------------+----------+-----------+---------------------+
| Mode             | Speedup  | Memory    | Requirements        |
+==================+==========+===========+=====================+
| CPU chunked      | 1.8x     | ~19 GB    | Default, any system |
+------------------+----------+-----------+---------------------+
| CPU legacy       | 1.0x     | ~38 GB    | cpu_chunk_size=0    |
+------------------+----------+-----------+---------------------+
| GPU monolithic   | 6.7x     | 15 GB VRAM| 24+ GB GPU          |
+------------------+----------+-----------+---------------------+
| GPU chunked      | 4.7x     | 0.9 GB    | Any GPU (4+ GB)     |
+------------------+----------+-----------+---------------------+

Benchmarked on Allen Institute Human Neocortex (148K cells x 30K genes).

All four modes read an in-memory ``AnnData``, so peak memory still scales with
the size of the object. To keep memory bounded by the chunk size instead, read
straight from a file — see `Cytome datasets`_.

COSG does not modify the input ``adata.X``. It is safe to call
``cosg.cosg()`` without copying ``adata`` first.

Optional extras
---------------

``pip install cosg`` covers marker detection, ``plotMarkerDendrogram`` and
``plotMarkerStream``. Two features need extras::

    pip install 'cosg[dotplot]'   # plotMarkerDotplot (wraps scanpy.pl.dotplot)
    pip install 'cosg[gpu]'       # the CuPy GPU path, device='gpu'

scanpy was a hard dependency through 1.0.4. It is now optional because only
``plotMarkerDotplot`` uses it, and requiring it cost every install ~60 MB and
9 packages (statsmodels, seaborn, umap-learn, ...) for one plotting function.
Calling ``plotMarkerDotplot`` without it raises an error naming the extra.

Cytome datasets
---------------

Marker genes can be computed directly from a ``.cytome`` file, without loading
the matrix into memory and without going through AnnData:

.. code-block:: bash

    pip install cytome

.. code-block:: python

    import cosg

    markers = cosg.cosg(
        "atlas.cytome",
        groupby="cell_type",
        modality="RNA",
        layer="counts",
        n_genes_user=50,
    )

``cosg.cosg()`` dispatches on its first argument. An ``AnnData`` takes the
in-memory path and writes ``adata.uns[key_added]``; a path to a ``.cytome``
file takes the streaming path, reads the file in chunks and returns a dict of
``names``, ``scores`` and ``groups_order``. Peak memory is set by the chunk
size, not by the number of cells.

``layer=`` names a matrix stored in the file. Omit it and the default is
resolved from the modality — ``RNA``/``GA`` to ``log1p``, ``ATAC``/``tiles``
to ``tfidf`` — which normalizes on the fly and therefore also needs
``pip install piaso``. Pass a stored layer, as above, to run with cosg and
cytome alone. ``output_format=`` selects ``dict``, ``long`` or ``dense``.

cytome is a single-file, SQLite-backed format that holds matrices, cell and
feature metadata, embeddings and genomic fragments together:
https://github.com/genecell/cytome
