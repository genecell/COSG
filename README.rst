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

Tutorials — marker genes and their significance, on a cytome, across batches, on the GPU, and on
spatial data — are at
`piaso.org/tutorials/cosg-markers <https://piaso.org/tutorials/cosg-markers/>`_.

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

**Release v1.2.1** (October 28, 2026)

*Re-run if affected.* Four changes move results that may have been reported:

- **Fixed:** ``use_raw=True`` reported the wrong gene names in COSG 1.0.4 to 1.2.0. Affected: ``cosg.cosg(adata, use_raw=True)`` on an AnnData whose ``adata.raw`` holds more genes than ``adata``, which is scanpy's usual layout after selecting variable genes. COSG scored ``adata.raw.X`` but named its columns by ``adata.var_names``, so a marker at column *j* of raw was reported as the gene at position *j* of ``adata`` (a group's markers g2, g3 came back as g32, g33), or raised IndexError past its end. The names, fold changes and p-values now describe raw's genes. In 1.0.1 to 1.0.3 ``use_raw`` had no effect (``adata.X`` was scored, and named correctly). Re-run any ``use_raw=True`` analysis made with 1.0.4 or later.

- **Fixed:** stratified p-values with a group alone in a batch. Affected: ``cosg.cosg(..., batch_key=..., calculate_pvalues=True)`` where some batch holds cells of only one group (sorted populations per sample, one type per donor). That batch's sum entered the statistic but not its mean, so on null data the group's p-values were near 1e-80 and every gene it expresses was called. A batch that holds none of a group's cells, or nothing but them, carries no information about the group and is now left out of its test (statistic, mean, variance and saddlepoint together); the p-values match a stratified permutation reference. Strata of one or two cells now count when they inform the group, as in the permutation null. Runs without ``batch_key``, or whose every batch mixes groups, are unchanged; re-run the others.

- **Changed:** the normalized scores of a group that scores above zero at fewer than 5 % of the features. Affected: IQR-normalized scores, that is ``iqr_normalize`` on a cytome (on by default for the ``dict``, ``long`` and ``dense`` outputs) and ``cosg.iqrLogNormalize``, mostly on peaks and tiles, where a rare or shallowly sequenced cell type can score above zero at a few percent of them. The spread of such a group, between its ``q_upper`` and ``q_lower`` quantiles (the 95th and 75th percentiles by default), is zero: a cytome path then divided its scores by 1e-6, which put them about log(scale / 1e-6) above other groups', and ``iqrLogNormalize`` gave every column without a spread the smallest spread of another column. A group without a spread now keeps a scale of its own: its ``q_upper`` quantile when that is positive, else the 5th percentile of its positive scores (a group positive at 5 to 5.3 % of the features has its spread floored there, so the scale does not jump as the share crosses 5 %). The normalized scores of other groups are unchanged, and no ranking within a group changes. ``iqr_fallback='constant'`` on a cytome and ``iqr_fallback='smallest'`` in ``iqrLogNormalize`` reproduce 1.2.0. Re-run if you compared normalized scores across groups and one of them was this sparse.

- **Fixed:** ``cosg.cosg(adata, device='gpu', batch_key=...)`` applied no expression filter in 1.1.0 to 1.2.0, so ``remove_lowly_expressed`` and ``expressed_pct`` had no effect there and a gene detected in a handful of a group's cells could be one of its markers. The GPU batch path counts no nonzeros, and the filter after it was skipped as for the GPU paths that filter on the device. It is now applied after it, as after the CPU's batch path, and removes the same genes. Re-run such an analysis; its markers change wherever a gene was below the floor.

*Other changes to results.*

- **Changed:** p-values deep in the tail. The saddlepoint stopped at a normal p of 1e-30, which left the sparse markers, exactly the ones whose normal p is wrong by orders of magnitude, at their normal value: a gene in 12 cells, all of a 100-cell group of 5,000, had log10 p -129.7 against an exact -20.7. Sparse entries (fewer than 30 expressing cells expected in the group) are now refined however deep (-20.2 here; within 0.1 to 0.8 log10 of the exact tail); dense genes that deep keep the normal approximation's value, as before. No p is reported below the permutation floor, 1 / C(n, n_c) (per batch with ``batch_key``), the probability of the observed labeling itself. ``neg_log10_pvals`` of sparse and of the very strongest markers fall; decisions at any usual alpha do not change. ``pvalue_method='exact'`` no longer reports 0 for a strong marker whose tail lies past its table's cut (the saddlepoint takes those entries).

- **Changed:** the groups in the result are the groups scored that have cells. ``groups=`` with a categorical ``groupby`` returned every category, the unselected ones with arbitrary genes, and ``calculate_logfoldchanges`` then raised IndexError; an unused category did the same with ``groups='all'``. Such categories are now left out, as on the cytome paths; a sequence of per-group floors has one fraction per group in the result. Fold changes are those of the groups scored (``reference`` may be a label outside ``groups``) and no longer need a categorical ``groupby``. PIASO's supervised ``runGDR`` with ``batch_key`` loses, for a batch that lacks a type, the gene set of arbitrary genes that type used to contribute.

- **Changed:** one rule for a group's top features on every path: by score, ties by feature index (``argpartition`` and the passes of the feature-batched paths chose among tied scores differently, e.g. 10 tied peaks gave {p1, p4, p5} or {p0, p1, p2}). ``n_genes_user`` above the number of features returns them all (it raised). The ndarray output of the feature-batched cytome paths lists a group's next features by score, zeros and then -1, instead of "" and 0.0.

- **Changed:** "expressed" means a nonzero value on every path. The streaming, GPU and cytome paths counted stored entries, explicit zeros included, and dense streaming input counted values above 0. Results change only for matrices that store zeros or hold negative values.

- **Changed:** ``expressed_pct=0`` (or a group's 0) has one meaning on every path: the ``expressed_min_num_cells_in_target_group`` minimum alone, 3 cells by default. It used to switch the filter off altogether, the minimum included, on the default chunked CPU path, the GPU's chunked path and every path of ``run_cosg_cytome``, and to keep the minimum on the others, so a feature seen in one or two cells of a group could be among its markers or not depending on the path. Results change only for a 0 on those paths, where such features are no longer markers; ``remove_lowly_expressed=False`` switches the filter off.

- **Changed:** the cytome's unbatched path (the RNA default) sums in float64, as the AnnData paths do; its scores move at the level of float32 rounding of the sums (1.5e-6 over a store of many chunks). The feature-batched paths keep float32 sums. log1p normalization is vectorized (80 times faster on the column chunks of a column-major matrix), with the same numbers to the bit.

- **Changed:** a default call writes nothing into a cytome. A path is opened read-only (cytome 0.3.5 or later), and log1p takes each cell's depth from the data, the row sums of the chunk it normalizes, instead of the ``RNA_cell_depth`` cache, which was trusted whenever its length matched (a matrix replaced by one of the same shape kept the old depths). The depths are the ones the cache held when it was fresh, so results do not move; the first call on a store is faster (no separate pass for the depth). PIASO's infog and tfidf statistics are then recomputed by each call on a path rather than cached in the store, which costs a store without them one more pass over the matrix per call (about 25 % more time in our test); pass a Dataset opened read-write to keep caching them. ``write_to_cytome=True`` writes only its markers, and a store it cannot write fails before the computation.

- **Changed:** ``plotMarkerDotplot`` draws the plot itself. scanpy's ``sc.pl.dotplot``, which it called, cut the group labels beside the brackets to twice the number of genes per group (with ``swap_axes=True`` and ``top_n_genes=2``, 'L6 CT', 'L6 IT' and 'L6 IT Car3' all read 'L6 .') and drew the dendrogram on a scale of its own, out of register with the rows. The default renderer (``backend='native'``) needs no scanpy, keeps whole labels and draws the dendrogram on the rows; ``backend='scanpy'`` makes the earlier call. ``save`` is now a file path (scanpy pasted it after ``figures/dotplot_``, so an absolute path raised FileNotFoundError). The function returns None with ``show=True`` (the default), as before, and the dot grid's axes with ``show=False``. New arguments: ``show``, ``dpi``, ``expression_cutoff`` and ``backend``.

*New.*

- **A detection floor per group.** ``expressed_pct`` also takes one fraction per group: a dict ``{group: fraction}``, a pandas Series indexed by group, or a sequence in the order of the result's groups. The share of cells that suits an abundant, deeply sequenced type need not suit a rare or shallow one, and on peaks and tiles a floor read from the data can be an order of magnitude below 0.1. Every path honours it: ``cosg.cosg`` on CSR, CSC or dense input, chunked or not, with ``batch_key`` or p-values, on the CPU or the GPU; ``run_cosg_cytome`` unbatched, batched by chromosome or by disk, on row- or column-major matrices, and on the GPU; and ``score_from_sums``.

- A group's floor removes features from that group's markers only, after the cosine and the penalty, which do not involve it. So a group's markers under a per-group floor are those it gets with its own fraction as the single floor, and the same value for every group gives the result of that value as a number, to the bit. One thing couples groups, not the floor: under ``scale_rule='positive_median'`` a group with fewer than 20 positive scores takes the median of the other groups' scales, which their floors move.

- With a dict or a Series every group with cells needs a fraction: a group left out raises rather than silently taking 0.1. A key that names no label of ``groupby`` raises; a key for a label not scored in the call (left out by ``groups`` or ``cell_mask``) is ignored, so one table of every type serves any mask. A table of per-type cell counts is passed as fractions ``(count - 0.5) / n_cells``, which make the threshold "count cells or more" whatever the rounding.

- The floors used are recorded: in ``adata.uns[key_added]['params']['expressed_pct']`` as a one-record array with a field per group (the form of ``names``, which survives ``adata.write`` with labels such as ``L2/3 IT``), and in ``run_cosg_cytome``'s result under ``'expressed_pct'``, as ``{group: fraction}``, which ``write_to_cytome`` stores with the markers. A run with a single number returns exactly what it did before.

- ``expressed_pct`` is checked on the cytome paths too: a value outside [0, 1] raises there as it did on the AnnData path. ``expressed_pct=None`` means 0.1 on every path; the ``cpu_chunk_size=0`` path used to raise a TypeError for it.

- **Column-major cytome matrices.** A matrix stored as column (CSC) chunks only, which earlier versions could not read, is scored a block of features at a time: a column chunk holds every entry of its features, so memory follows the block, not the number of features times the number of groups. On counts the scores are those of the same matrix stored by rows.

- ``cosg.score_from_sums(dot, sq, group_sizes, nnz)`` returns COSG's raw scores from the sums every streaming path accumulates: per feature and group, the sum of the feature's values; per feature, the sum of its squares; the group sizes; and, for the expression filter, the nonzero counts. Sums of parts add up, so the sums of count-split halves, of shards or of added cells can be scored without the matrix, and per-batch sums (``batch_group_sizes``) give the batch-corrected scores of ``batch_key``. It is the last step of the streaming paths, whose results are unchanged by it to the bit.

- ``scale_rule='positive_median'`` (``run_cosg_cytome``, ``cosg.cosg`` on a cytome, ``iqrLogNormalize``) scales each group by the median of its positive scores, for sparse features (peaks, tiles), where most groups score above zero at 5 to 25 % of the features and the spread over all features is set by how many of them score at all (depth, the expression filter, the size of the feature set) more than by the group's own scores. A group with fewer than 20 positive scores takes the median of the other groups' scales. The default stays ``'iqr'``. A scale divides a group's scores, so no rule changes a ranking within a group.

- ``max_accumulator_bytes`` (``run_cosg_cytome``, 2 GiB by default): on a row-major matrix, chromosome batching holds every chromosome's per-feature sums in one pass only while they fit; beyond it the chromosomes are packed into passes that fit, so memory stops growing with the number of features times the number of groups (a whole-brain atlas's tiles by its subclasses would need about 13 GB). The scores are the same to the bit.

- ``neg_log10_pvals_adj``, beside ``neg_log10_pvals``: -log10 of the Benjamini-Hochberg adjusted p, computed in log space, so it stays finite where ``pvals_adj`` underflows to 0. ``neg_log10_pvals`` is of the raw p, so a line at 1.3 on it is raw 0.05, not FDR 0.05; this one is safe to threshold.

*Fixed.*

- **Fixed:** p-values on a cytome. Strong markers got exactly 0.0; p is now floored at the smallest positive float64, and ``neg_log10_pvals`` and ``neg_log10_pvals_adj`` are reported, as on the AnnData path (they agree to 1e-6 on the same counts).

- **Fixed:** ``device='gpu'`` (and ``device='auto'`` on a machine with CuPy) on a cytome dropped ``calculate_pvalues`` without a word and did not apply ``iqr_normalize``. The GPU path now shares the CPU's post-processing: every ``output_format`` and the scale options work there. What it cannot do (p-values, feature batching, a column-major matrix) raises for ``device='gpu'`` and runs on the CPU for ``device='auto'``.

- **Fixed:** the chunked GPU path (chosen automatically when the matrix exceeds 70 % of free GPU memory) marked filtered features with -inf instead of -1, which ``iqrLogNormalize`` turned into NaN. With ``cpu_chunk_size=0``, ``batch_key`` and a dense matrix, a group too small in every batch made every score NaN; it now scores 0, as on the other batch paths. ``gpu_chunk_size=0`` (the documented monolithic transfer) raised ZeroDivisionError.

- **Fixed:** a NaN or an infinity in the matrix made its gene the top marker of every group on the default path; non-finite values are now refused, naming the features (on a cytome, at no extra cost, from the sums of squares). A cell without a label raised KeyError or TypeError; it is refused with a count. ``mu=NaN`` made every score NaN. A call that failed wiped the earlier result under the same ``key_added``; it now leaves it in place.

- **Fixed:** ``cell_mask`` on a cytome. A boolean mask shorter than the cells was accepted and the last cells silently left out, and integer indices with ``batch_key`` raised IndexError. Both forms are validated (length, range, duplicates) and give the same result in any order.

- **Fixed:** a feature name held twice in a cytome's dense output. A 10x GRCh38 reference holds 10 gene symbols twice (HSPA14, TBCE, MATR3, ...). ``output_format='dense'`` mapped names to rows, so one row of each pair stayed 0 and the other took either feature's score; it is now filled by feature index, each row with its own score. The ``dict`` and ``long`` outputs, keyed by name, keep the better score of such a pair and warn.

- **Fixed:** chromosome batching assumed each chromosome's features contiguous; on a store not ordered by chromosome, features were scored twice and the IQR scale padded to more features than there are. The batches are now runs of one chromosome. A column-major matrix with the default ``layer='auto'`` crashed on the unbatched path. Column chunks are normalized with their own features' TF-IDF and INFOG factors.

- **Fixed:** a ``groupby`` column missing from a cytome's cells table was read by SQLite as a string, so every cell took that name as its label and the run reported one group named after it, with a plausible marker table. A missing column now raises, naming the columns there and a match that differs only in case; such a name is accepted, as SQLite resolves it.

- **Fixed:** ``min_total_count`` on the feature-batched path of ``run_cosg_cytome`` was dropped with ``expressed_pct=0`` or ``remove_lowly_expressed=False``, since it needed the nonzero counts only the floor kept. It now applies whenever it is above 0, as documented.

- **Fixed:** a group label holding "/" (``L2/3 IT``) left an AnnData that ``adata.write`` could not save: h5ad and zarr keys cannot hold "/", and anndata 0.13 refuses to write them by default. The by-group frame ``adata.uns[key]['COSG']`` now writes "/" as "%2F" in its column names (``names::L2%2F3 IT``); ``indexByGene``, and with it the plots that read the frame, reads it back as "/", for frames from earlier versions too, and the record arrays (``names``, ``scores``, ...) keep the labels exactly. Labels without "/" are written as before.

- **Fixed:** ``run_cosg_cytome(output_format='dense')`` raised "assignment destination is read-only" with pandas 3, whose copy-on-write makes a frame's ``.values`` read-only.

- ``run_cosg_cytome`` checks ``mu`` and ``n_genes_user`` as ``cosg.cosg`` does, and ``expressed_min_num_cells_in_target_group=0`` is accepted on both paths. COSGR, the R package, applies ``expressed_pct`` without the 3-cell minimum, so for groups under about 30 cells it can report other markers.

No dependency changes: the ``cosg[cytome]`` extra still accepts cytome 0.2.3 or later (before 0.3.5 a path is opened for writing, as in earlier releases).


**Release v1.2.0** (September 3, 2026)

- **New:** p-values for COSG's markers. ``cosg.cosg(..., calculate_pvalues=True)`` adds ``pvals``, ``pvals_adj`` (Benjamini-Hochberg within each group, or ``pvalue_fdr_method='fdr_by'``), ``zscores`` and ``neg_log10_pvals`` (computed in log space, with no floor) to the result. It is off by default, and the output is unchanged when it is off. Under the null that a gene's expression does not depend on the labels, the test is of the gene's sum over the group, drawn without replacement, with exact moments: ``pvalue_method='spa'`` (the default) uses the normal approximation in the bulk, a saddlepoint in the tail and exact enumeration for genes in very few cells; ``'normal'``, ``'exact'`` and ``'permutation'`` are the alternatives. ``batch_key`` stratifies the null on an AnnData; a cytome gives the same p-values as an AnnData and refuses ``batch_key`` with them. Closes `issue #6 <https://github.com/genecell/COSG/issues/6>`_.

- If the groups came from clustering the same matrix, these p-values are anti-conservative; labels from another source, or count splitting, avoid it. The `tutorial <https://piaso.org/tutorials/cosg-markers/>`_ shows both.

- Includes the 1.1.2 and 1.1.3 fixes.


**Release v1.1.3** (August 24, 2026)

- **Fixed:** ``run_cosg_cytome(layer='auto')`` normalized an already-normalized matrix a second time. ``auto`` chose ``log1p`` for RNA from a table keyed on the modality alone, without looking at the data. Because ``cytome.from_anndata`` stores whatever ``adata.X`` held under ``{modality}_counts`` — including a log-normalized matrix — any file written from normalized input had ``log1p`` applied twice, while the function is documented as equivalent to ``cosg.cosg(adata)``. Correlation with the in-memory reference stayed around 0.99 and specificity moved by up to 0.11: small enough to pass a glance, large enough to change which genes rank as markers. ``auto`` now probes the stored matrix and, when it is demonstrably not integer, uses the values as given and warns, naming the override. An explicit ``layer=`` is obeyed untouched.

- **Fixed:** the version had two sources of truth. ``pyproject.toml`` and ``cosg/__init__.py`` each carried a literal, and they agree only until a release bumps one of them. ``pyproject`` now declares the version dynamic and reads ``cosg.__version__``.

- **Interoperates with cytome 0.3.0.** That release stops ``from_anndata`` writing a non-integer ``adata.X`` to ``{modality}_counts``, so newer files carry ``{modality}_data`` (or a name you chose) instead. ``layer='auto'`` recognises this: it reads the recorded ``_anndata_X_layer``, uses that matrix, and says which one it picked. Where cytome records ``matrix_meta.is_integer``, that is preferred over probing the values. The minimum cytome version is unchanged at 0.2.3 — this release exists to read *older* files correctly, and requiring 0.3.0 would help none of them.

- If you have run ``run_cosg_cytome(layer='auto')`` on a cytome converted from a normalized AnnData, re-run it. Runs on raw-count files are unaffected, and any run that passed an explicit ``layer=`` was always correct.


**Release v1.1.2** (August 22, 2026)

- **Fixed:** ``cosg.cosg(adata)`` raised ``ImportError`` on a default install. The cytome ``Dataset`` alias lookup was not guarded like the first one, so a missing optional cytome dependency produced an error instead of the documented empty-tuple fallback. Plain AnnData no longer requires cytome to be installed.

- **Fixed:** ``plt.cm.get_cmap``, removed in matplotlib 3.9, replaced with ``plt.get_cmap``.


**Release v1.1.1** (August 15, 2026)

- Fixed the PyPI project page, which showed a single sentence instead of this README. The ``readme`` field pointed at an inline string rather than at ``README.rst``, so the long description was never packaged. Code is identical to v1.1.0.


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
