Usage
------

Import COSG as:
::
        import cosg as cosg

amd import Scanpy as:
::
        import scanpy as sc

Next, load the data via:
::
        adata = sc.datasets.pbmc68k_reduced()

then identify marker genes for each cell group by running:
::
        cosg.cosg(adata, key_added='cosg', groupby='bulk_labels')

and the top marker genes can be visualized via:
::
        sc.pl.rank_genes_groups(adata, key='cosg')
