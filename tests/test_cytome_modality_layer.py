"""COSG cytome streaming: modality (incl. GA) + layer support tests.

Two distinct concerns:

1. Modality dispatch — pre-fix, `_get_feature_table_info(ds, 'GA')` and
   `_get_labels_and_gene_names(ds, ..., modality='GA')` fell into the
   `else` branch and read the RNA `genes` table instead of `GA_genes`.
   Now both delegate to `cytome.utils.modality.modality_feature_table_info`
   so the GA branch is covered.

2. Layer support — pre-fix, the cytome COSG path always read `{modality}_counts`
   regardless of any user-specified `cytome_layer`. The new
   `cytome_layer` / `compute_on_fly` / `use_cached_stats` kwargs let
   COSG run on log1p / infog / tfidf — either reading the materialised
   matrix or computing per-chunk on the fly via piaso normalizers
   (lazy-imported only when needed).
"""
from __future__ import annotations

import sys
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp


def _build_multimodal_cytome(path, n_cells=20, seed=0):
    """RNA + GA + ATAC cytome with a deterministic per-modality value
    pattern so we can assert COSG saw the right matrix."""
    import cytome
    rng = np.random.default_rng(seed)
    leiden = np.array([f"g{i % 3}" for i in range(n_cells)])
    ds = cytome.create(path)
    ds.set_entity("cells", pd.DataFrame({
        "cell_idx": np.arange(n_cells),
        "barcode": [f"AAA-{i}" for i in range(n_cells)],
        "Leiden": leiden,
    }))
    # RNA: Sox2 high in g0
    rna_X = np.zeros((n_cells, 3), dtype=np.float32)
    rna_X[leiden == "g0", 0] = 5.0
    rna_X[leiden == "g0", 1] = 1.0
    rna_X[leiden == "g1", 1] = 5.0
    rna_X[leiden == "g2", 2] = 5.0
    ds.set_entity("genes", pd.DataFrame({
        "gene_idx": [0, 1, 2],
        "gene_id": ["Sox2", "Pax6", "Foxg1"],
    }))
    ds.add_matrix("RNA_counts", sp.csr_matrix(rna_X))
    # GA: deliberately DIFFERENT gene set + different per-cluster pattern
    ga_X = np.zeros((n_cells, 3), dtype=np.float32)
    ga_X[leiden == "g2", 0] = 7.0   # Olig2 marker for g2
    ga_X[leiden == "g0", 1] = 7.0   # Nestin marker for g0
    ga_X[leiden == "g1", 2] = 7.0   # Tubb3 marker for g1
    ds.set_entity("GA_genes", pd.DataFrame({
        "gene_idx": [0, 1, 2],
        "gene_id": ["Olig2", "Nestin", "Tubb3"],
    }))
    ds.add_matrix("GA_counts", sp.csr_matrix(ga_X))
    # ATAC stub
    ds.set_entity("peaks", pd.DataFrame({
        "peak_idx": [0],
        "peak_id": ["chr1:0-100"],
        "chr": ["chr1"], "start": [0], "end_": [100],
    }))
    ds.add_matrix(
        "ATAC_counts",
        sp.csr_matrix(np.ones((n_cells, 1), dtype=np.float32)),
    )
    ds.flush()
    return ds


# -----------------------------------------------------------------------
# Modality dispatch — the GA bug regression
# -----------------------------------------------------------------------

def test_get_feature_table_info_ga_routes_to_GA_genes(tmp_path):
    """Pre-fix: modality='GA' fell into the else branch and read 'genes'.
    Post-fix: routes to GA_genes via cytome.utils.modality."""
    import cytome
    from cosg._cytome_streaming import _get_feature_table_info
    ds = _build_multimodal_cytome(tmp_path / "x.cytome")
    assert _get_feature_table_info(ds, "GA") == ("GA_genes", "gene_idx", "gene_id")
    # And for sanity — RNA still routes correctly
    assert _get_feature_table_info(ds, "RNA") == ("genes", "gene_idx", "gene_id")
    ds.close()


def test_get_labels_and_gene_names_ga_returns_GA_gene_names(tmp_path):
    """Pre-fix: gene_names for GA modality came from the RNA `genes` table.
    Post-fix: comes from GA_genes."""
    import cytome
    from cosg._cytome_streaming import _get_labels_and_gene_names
    ds = _build_multimodal_cytome(tmp_path / "x.cytome")
    labels, gene_names, _, _ = _get_labels_and_gene_names(
        ds, groupby="Leiden", modality="GA",
    )
    assert list(gene_names) == ["Olig2", "Nestin", "Tubb3"]
    # If pre-fix code ran, gene_names would have been ["Sox2", "Pax6", "Foxg1"]
    ds.close()


def test_run_cosg_cytome_ga_modality_picks_GA_markers(tmp_path):
    """End-to-end on GA modality: COSG should rank Olig2 highly for g2,
    Nestin for g0, Tubb3 for g1 (per the deterministic GA_counts pattern).
    Pre-fix this would have read RNA_counts and ranked Sox2/Pax6/Foxg1."""
    from cosg._cytome_streaming import run_cosg_cytome
    ds = _build_multimodal_cytome(tmp_path / "x.cytome")
    ds.close()
    result = run_cosg_cytome(
        cytome_path=str(tmp_path / "x.cytome"),
        groupby="Leiden",
        modality="GA",
        n_genes_user=3,
        verbose=False,
    )
    names = result["names"]
    groups = list(result["groups_order"])
    # names is a structured array — the top marker for each group
    # should match the pattern (Olig2 for g2, Nestin for g0, Tubb3 for g1)
    g0_top = names[groups.index("g0")][0]
    g1_top = names[groups.index("g1")][0]
    g2_top = names[groups.index("g2")][0]
    assert g0_top == "Nestin", f"g0 top should be Nestin, got {g0_top}"
    assert g1_top == "Tubb3", f"g1 top should be Tubb3, got {g1_top}"
    assert g2_top == "Olig2", f"g2 top should be Olig2, got {g2_top}"


# -----------------------------------------------------------------------
# Layer support — counts default + materialized + on-the-fly
# -----------------------------------------------------------------------

def test_run_cosg_cytome_default_layer_is_counts(tmp_path):
    """Backward compat: default layer='counts' produces same
    markers as before."""
    from cosg._cytome_streaming import run_cosg_cytome
    ds = _build_multimodal_cytome(tmp_path / "x.cytome")
    ds.close()
    result = run_cosg_cytome(
        cytome_path=str(tmp_path / "x.cytome"),
        groupby="Leiden",
        modality="RNA",
        n_genes_user=3,
        verbose=False,
    )
    # Default layer='counts' — Sox2 is the g0 marker per the pattern
    g0_top = result["names"][list(result["groups_order"]).index("g0")][0]
    assert g0_top == "Sox2"


def test_run_cosg_cytome_materialized_log1p_layer(tmp_path):
    """When `RNA_log1p` is materialised, layer='log1p' reads it
    directly (no on-the-fly normalization)."""
    import cytome
    from cosg._cytome_streaming import run_cosg_cytome
    ds = _build_multimodal_cytome(tmp_path / "x.cytome")
    # Add a fake RNA_log1p with a deliberately-different pattern so we can
    # detect the layer was actually read.
    n_cells = ds.n_cells
    leiden = np.asarray(ds.cells["Leiden"]).astype(str)
    log1p_X = np.zeros((n_cells, 3), dtype=np.float32)
    # Reverse the pattern: Foxg1 high in g0 (NOT Sox2)
    log1p_X[leiden == "g0", 2] = 10.0
    ds.add_matrix("RNA_log1p", sp.csr_matrix(log1p_X))
    ds.flush()
    ds.close()

    result = run_cosg_cytome(
        cytome_path=str(tmp_path / "x.cytome"),
        groupby="Leiden",
        modality="RNA",
        layer="log1p",
        compute_on_fly=False,  # require materialised
        n_genes_user=3,
        verbose=False,
    )
    # If we read RNA_log1p, the g0 marker should be Foxg1 (per the fake
    # pattern); if we wrongly fell back to counts, it would be Sox2.
    g0_top = result["names"][list(result["groups_order"]).index("g0")][0]
    assert g0_top == "Foxg1", (
        f"Materialised-layer read failed: g0 top should be Foxg1 (from "
        f"RNA_log1p), got {g0_top} (looks like the counts matrix was read)."
    )


def test_run_cosg_cytome_strict_layer_missing_raises(tmp_path):
    """compute_on_fly=False and the layer isn't materialised → actionable
    ValueError that names compute_on_fly=True as the fix."""
    from cosg._cytome_streaming import run_cosg_cytome
    ds = _build_multimodal_cytome(tmp_path / "x.cytome")
    ds.close()
    with pytest.raises(ValueError, match="compute_on_fly"):
        run_cosg_cytome(
            cytome_path=str(tmp_path / "x.cytome"),
            groupby="Leiden",
            modality="RNA",
            layer="infog",
            compute_on_fly=False,    # require materialised, but nothing's there
            n_genes_user=3,
            verbose=False,
        )


def test_run_cosg_cytome_unknown_layer_raises(tmp_path):
    """Unknown layer with compute_on_fly=True still raises (no on-the-fly
    path for it)."""
    from cosg._cytome_streaming import run_cosg_cytome
    ds = _build_multimodal_cytome(tmp_path / "x.cytome")
    ds.close()
    # Round 9 (2026-05-23): cytome_layer renamed to layer, so the
    # error message no longer contains the old name. Match on "layer"
    # (more general) plus a recognisable bit of the new error text.
    with pytest.raises((ValueError, RuntimeError), match="layer"):
        run_cosg_cytome(
            cytome_path=str(tmp_path / "x.cytome"),
            groupby="Leiden",
            modality="RNA",
            layer="totally_made_up_layer",
            compute_on_fly=True,
            n_genes_user=3,
            verbose=False,
        )


def test_run_cosg_cytome_signature_has_new_kwargs():
    """Anti-doc-drift / contract test: the new kwargs must exist with
    the right defaults. Round 9 (2026-05-23): ``cytome_layer`` renamed
    to ``layer``."""
    import inspect
    from cosg._cytome_streaming import run_cosg_cytome
    sig = inspect.signature(run_cosg_cytome)
    # `layer` defaults to "auto" since Round 26: it resolves a
    # modality-appropriate normalization (log1p for RNA, TF-IDF for ATAC)
    # rather than one fixed default, because COSG scores a cosine
    # specificity and no single transform is right for both modalities.
    assert sig.parameters["layer"].default == "auto"
    assert sig.parameters["compute_on_fly"].default is True
    assert sig.parameters["use_cached_stats"].default is True


def test_run_cosg_cytome_gpu_signature_has_new_kwargs():
    """Same for the GPU variant. Round 9: ``cytome_layer`` → ``layer``."""
    import inspect
    from cosg._cytome_streaming import run_cosg_cytome_gpu
    sig = inspect.signature(run_cosg_cytome_gpu)
    # `layer` defaults to "auto" since Round 26: it resolves a
    # modality-appropriate normalization (log1p for RNA, TF-IDF for ATAC)
    # rather than one fixed default, because COSG scores a cosine
    # specificity and no single transform is right for both modalities.
    assert sig.parameters["layer"].default == "auto"
    assert sig.parameters["compute_on_fly"].default is True
    assert sig.parameters["use_cached_stats"].default is True
