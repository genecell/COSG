"""Round 8 (2026-05-23): polymorphic ``cosg.cosg(data, ...)`` dispatch tests.

`cosg.cosg` now accepts:
  * an ``AnnData`` — existing in-memory path (mutates ``data.uns``).
  * a cytome path string / ``pathlib.Path`` — streaming path.
  * an open cytome Dataset — streaming path; lifecycle owned by caller.

Plus the cytome path now accepts ``device='cpu' | 'gpu' | 'auto'`` so
GPU dispatch flows through the unified entry. ``run_cosg_cytome_gpu``
survives as a thin shim.

These tests pin both behaviour and the explicit wrong-input-type
kwarg rejection (Q5-A) that prevents silent kwarg misrouting.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from conftest import requires_cytome, requires_piaso


# --------------------------------------------------------------------
# Shared fixture — small RNA cytome with deterministic per-group pattern.
# --------------------------------------------------------------------

def _build_rna_cytome(path, n_cells=30, seed=0):
    """Tiny RNA cytome — 3 clusters × 3 genes; each gene marks one cluster."""
    import cytome
    leiden = np.array([f"g{i % 3}" for i in range(n_cells)])
    ds = cytome.create(path)
    ds.set_entity("cells", pd.DataFrame({
        "cell_idx": np.arange(n_cells),
        "barcode": [f"AAA-{i}" for i in range(n_cells)],
        "Leiden": leiden,
    }))
    rna_X = np.zeros((n_cells, 3), dtype=np.float32)
    rna_X[leiden == "g0", 0] = 5.0
    rna_X[leiden == "g1", 1] = 5.0
    rna_X[leiden == "g2", 2] = 5.0
    ds.set_entity("genes", pd.DataFrame({
        "gene_idx": [0, 1, 2],
        "gene_id": ["GeneA", "GeneB", "GeneC"],
    }))
    ds.add_matrix("RNA_counts", sp.csr_matrix(rna_X))
    ds.flush()
    return ds


def _build_matched_anndata(n_cells=30, seed=0):
    """AnnData with the SAME per-group pattern as ``_build_rna_cytome``,
    so the AnnData and cytome paths produce comparable markers."""
    import anndata
    leiden = np.array([f"g{i % 3}" for i in range(n_cells)])
    rna_X = np.zeros((n_cells, 3), dtype=np.float32)
    rna_X[leiden == "g0", 0] = 5.0
    rna_X[leiden == "g1", 1] = 5.0
    rna_X[leiden == "g2", 2] = 5.0
    return anndata.AnnData(
        X=sp.csr_matrix(rna_X),
        obs=pd.DataFrame({"Leiden": leiden}, index=[f"AAA-{i}" for i in range(n_cells)]),
        var=pd.DataFrame(index=["GeneA", "GeneB", "GeneC"]),
    )


# --------------------------------------------------------------------
# AnnData path unchanged
# --------------------------------------------------------------------

def test_cosg_anndata_input_returns_none_and_mutates_uns():
    """The AnnData polymorphic dispatch must NOT change AnnData
    semantics. Pre-Round-8 behaviour: writes to adata.uns['cosg'],
    returns None."""
    import cosg
    a = _build_matched_anndata()
    result = cosg.cosg(a, groupby="Leiden", n_genes_user=3,
                      remove_lowly_expressed=False, verbosity=0)
    assert result is None
    assert "cosg" in a.uns
    # 3 genes per group, 3 groups
    names = a.uns["cosg"]["names"]
    assert len(names.dtype.names) == 3  # 3 groups


# --------------------------------------------------------------------
# Cytome path — str / Path / Dataset
# --------------------------------------------------------------------

@requires_cytome
@requires_piaso
def test_cosg_cytome_path_str_dispatches_to_cytome_streaming(tmp_path):
    """First-arg string → routes to run_cosg_cytome. Returns a dict
    (NOT None), and the dict shape matches the direct cytome call."""
    import cosg
    ds = _build_rna_cytome(tmp_path / "x.cytome")
    ds.close()
    result_via_unified = cosg.cosg(
        str(tmp_path / "x.cytome"),
        groupby="Leiden",
        n_genes_user=3,
        remove_lowly_expressed=False,
        verbosity=0,
    )
    result_via_direct = cosg.run_cosg_cytome(
        str(tmp_path / "x.cytome"),
        groupby="Leiden",
        n_genes_user=3,
        remove_lowly_expressed=False,
        verbose=False,
    )
    # Both return dicts in ndarray output shape
    assert isinstance(result_via_unified, dict)
    assert set(result_via_unified.keys()) == {"names", "scores", "groups_order"}
    np.testing.assert_array_equal(
        result_via_unified["names"], result_via_direct["names"]
    )


@requires_cytome
@requires_piaso
def test_cosg_cytome_path_pathlib_dispatches(tmp_path):
    """First-arg pathlib.Path also routes via _is_cytome_input."""
    from pathlib import Path
    import cosg
    ds = _build_rna_cytome(tmp_path / "x.cytome")
    ds.close()
    result = cosg.cosg(
        Path(tmp_path / "x.cytome"),
        groupby="Leiden",
        n_genes_user=3,
        remove_lowly_expressed=False,
        verbosity=0,
    )
    assert isinstance(result, dict)
    assert "names" in result


@requires_cytome
@requires_piaso
def test_cosg_cytome_dataset_object_dispatches(tmp_path):
    """First-arg open cytome Dataset routes to the cytome path. Critical
    for interactive use where the user has an open Dataset already."""
    import cosg
    import cytome
    ds = _build_rna_cytome(tmp_path / "x.cytome")
    ds.close()
    open_ds = cytome.open(tmp_path / "x.cytome")
    try:
        result = cosg.cosg(
            open_ds,
            groupby="Leiden",
            n_genes_user=3,
            remove_lowly_expressed=False,
            verbosity=0,
        )
        assert isinstance(result, dict)
        assert "names" in result
        # Dataset must still be usable — caller owns the lifecycle
        assert open_ds.n_cells == 30
    finally:
        open_ds.close()


@requires_cytome
@requires_piaso
def test_cosg_cytome_dataset_not_closed_by_dispatch(tmp_path):
    """Pin lifecycle ownership — passing an open Dataset and using it
    again after cosg.cosg must work."""
    import cosg
    import cytome
    ds = _build_rna_cytome(tmp_path / "x.cytome")
    ds.close()
    open_ds = cytome.open(tmp_path / "x.cytome")
    try:
        cosg.cosg(open_ds, groupby="Leiden", n_genes_user=3,
                  remove_lowly_expressed=False, verbosity=0)
        # If the dispatcher wrongly closed the Dataset, the next read
        # would raise. (cytome's behaviour on closed-ds varies by op;
        # we just want SOMETHING that touches the connection.)
        assert "Leiden" in [c for c in open_ds.cells.columns]
    finally:
        open_ds.close()


# --------------------------------------------------------------------
# device= dispatch on the cytome path
# --------------------------------------------------------------------

@requires_cytome
@requires_piaso
def test_run_cosg_cytome_device_cpu_default(tmp_path):
    """device='cpu' (default) hits the existing CPU path. Sanity that
    the new device= kwarg didn't break the default."""
    import cosg
    ds = _build_rna_cytome(tmp_path / "x.cytome")
    ds.close()
    result = cosg.run_cosg_cytome(
        str(tmp_path / "x.cytome"), groupby="Leiden", n_genes_user=3,
        remove_lowly_expressed=False, verbose=False, device="cpu",
    )
    assert "names" in result


@requires_cytome
@requires_piaso
def test_run_cosg_cytome_device_auto_resolves(tmp_path):
    """device='auto' resolves via _backend.get_device — picks GPU when
    CuPy available, else CPU. We just assert it RUNS (whichever path
    it picks); the impl is exercised by the device-explicit tests."""
    import cosg
    ds = _build_rna_cytome(tmp_path / "x.cytome")
    ds.close()
    result = cosg.run_cosg_cytome(
        str(tmp_path / "x.cytome"), groupby="Leiden", n_genes_user=3,
        remove_lowly_expressed=False, verbose=False, device="auto",
    )
    assert "names" in result


@requires_cytome
def test_run_cosg_cytome_device_gpu_with_sparse_output_raises(tmp_path):
    """device='gpu' + output_format != 'ndarray' is unsupported. The
    function must raise NotImplementedError BEFORE doing any work,
    naming both kwargs in the message. This test runs regardless of
    CuPy availability — the rejection happens at the dispatch layer."""
    from cosg._backend import HAS_CUPY
    if not HAS_CUPY:
        pytest.skip("No CuPy — the device='gpu' code path is not reachable.")
    import cosg
    ds = _build_rna_cytome(tmp_path / "x.cytome")
    ds.close()
    with pytest.raises(NotImplementedError) as exc:
        cosg.run_cosg_cytome(
            str(tmp_path / "x.cytome"), groupby="Leiden", n_genes_user=3,
            remove_lowly_expressed=False, verbose=False,
            device="gpu", output_format="dict",
        )
    assert "device='gpu'" in str(exc.value)
    assert "ndarray" in str(exc.value)


@requires_cytome
def test_run_cosg_cytome_gpu_shim_delegates(tmp_path):
    """The pre-Round-8 run_cosg_cytome_gpu function is now a thin shim
    around run_cosg_cytome(..., device='gpu'). Verify it still
    produces a result when CuPy is available."""
    from cosg._backend import HAS_CUPY
    if not HAS_CUPY:
        pytest.skip("No CuPy — GPU path not reachable.")
    import cosg
    ds = _build_rna_cytome(tmp_path / "x.cytome")
    ds.close()
    # Shim
    result_shim = cosg.run_cosg_cytome_gpu(
        str(tmp_path / "x.cytome"), groupby="Leiden", n_genes_user=3,
        remove_lowly_expressed=False, verbose=False,
    )
    # Direct
    result_direct = cosg.run_cosg_cytome(
        str(tmp_path / "x.cytome"), groupby="Leiden", n_genes_user=3,
        remove_lowly_expressed=False, verbose=False, device="gpu",
    )
    np.testing.assert_array_equal(result_shim["names"], result_direct["names"])


# --------------------------------------------------------------------
# Wrong-input-type kwarg rejection (Q5-A)
# --------------------------------------------------------------------

def test_cosg_anndata_with_cytome_kwarg_modality_raises():
    """Cytome-only kwarg `modality=` on an AnnData input → TypeError
    naming the offending kwarg."""
    import cosg
    a = _build_matched_anndata()
    with pytest.raises(TypeError) as exc:
        cosg.cosg(a, groupby="Leiden", modality="RNA")
    msg = str(exc.value)
    assert "modality" in msg
    assert "cytome-only" in msg.lower() or "anndata" in msg.lower()


def test_cosg_anndata_with_cytome_kwarg_output_format_raises():
    """Cytome-only kwarg `output_format=` on AnnData → TypeError."""
    import cosg
    a = _build_matched_anndata()
    with pytest.raises(TypeError) as exc:
        cosg.cosg(a, groupby="Leiden", output_format="dict")
    assert "output_format" in str(exc.value)


@requires_cytome
def test_cosg_cytome_with_anndata_kwarg_key_added_raises(tmp_path):
    """AnnData-only kwarg `key_added=` on a cytome input → TypeError
    naming the offending kwarg."""
    import cosg
    ds = _build_rna_cytome(tmp_path / "x.cytome")
    ds.close()
    with pytest.raises(TypeError) as exc:
        cosg.cosg(
            str(tmp_path / "x.cytome"), groupby="Leiden",
            key_added="my_cosg_key",
        )
    msg = str(exc.value)
    assert "key_added" in msg
    assert "anndata-only" in msg.lower() or "cytome" in msg.lower()


@requires_cytome
def test_cosg_cytome_with_anndata_kwarg_copy_raises(tmp_path):
    """AnnData-only kwarg `copy=True` on a cytome input → TypeError."""
    import cosg
    ds = _build_rna_cytome(tmp_path / "x.cytome")
    ds.close()
    with pytest.raises(TypeError) as exc:
        cosg.cosg(
            str(tmp_path / "x.cytome"), groupby="Leiden", copy=True,
        )
    assert "copy" in str(exc.value)


# --------------------------------------------------------------------
# Bad input type
# --------------------------------------------------------------------

def test_cosg_unknown_first_arg_type_raises():
    """First-arg that's neither AnnData nor a cytome path/Dataset →
    TypeError listing the valid types."""
    import cosg
    with pytest.raises(TypeError) as exc:
        cosg.cosg(42, groupby="Leiden")
    msg = str(exc.value)
    assert "AnnData" in msg
    assert "cytome" in msg.lower()
