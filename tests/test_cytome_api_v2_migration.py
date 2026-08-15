"""K-HARD migration tests for the COSG cytome API.

Two breaking changes have been shipped:

* **Round 7** (2026-05-22): collapsed ``top_k`` + ``n_top_genes`` →
  ``n_genes_user`` + explicit ``output_format``.
* **Round 9** (2026-05-23): function rename
  ``run_cosg_cytome_cpu`` → ``run_cosg_cytome``, and kwarg rename
  ``cytome_layer`` → ``layer``.

K-HARD means: passing the old names raises ``TypeError`` immediately
(no silent aliasing). These tests pin both the error type AND the
actionable migration string so a future "let's be nice and accept the
old name as an alias" refactor doesn't accidentally regress the
explicit-failure contract.

Despite the file name (``_v2`` was the Round-7 validator name), this
file accumulates K-HARD migration tests across rounds — additions are
welcome here rather than fragmenting into ``_v3``, ``_v4``, …
"""
from __future__ import annotations

import pytest

from conftest import requires_piaso

pytest.importorskip(
    "cytome",
    reason="cytome not installed — this whole module tests the .cytome "
           "streaming backend (pip install -e '.[dev]')",
)


def _stub_cytome_path(tmp_path):
    """We don't actually need a real cytome — argument validation runs
    before ``cytome.open(path)``."""
    return str(tmp_path / "nonexistent.cytome")


# =====================================================================
# Round 7 K-HARD: top_k and n_top_genes removed
# =====================================================================

def test_top_k_kwarg_raises_typeerror_with_migration_hint(tmp_path):
    from cosg._cytome_streaming import run_cosg_cytome
    with pytest.raises(TypeError) as exc:
        run_cosg_cytome(
            _stub_cytome_path(tmp_path),
            groupby="leiden",
            top_k=2000,         # old kwarg — must fail loudly
        )
    msg = str(exc.value)
    assert "top_k" in msg, f"error msg should name 'top_k': {msg!r}"
    assert "n_genes_user" in msg, (
        f"error msg should mention the new kwarg name: {msg!r}"
    )
    assert "output_format" in msg, (
        f"error msg should mention the new output_format kwarg: {msg!r}"
    )
    # Migration hint should include the user's value so they can copy/paste
    assert "2000" in msg, (
        f"error msg should echo the user's value: {msg!r}"
    )


def test_n_top_genes_kwarg_raises_typeerror_with_migration_hint(tmp_path):
    from cosg._cytome_streaming import run_cosg_cytome
    with pytest.raises(TypeError) as exc:
        run_cosg_cytome(
            _stub_cytome_path(tmp_path),
            groupby="leiden",
            n_top_genes=50,     # renamed kwarg — must fail loudly
        )
    msg = str(exc.value)
    assert "n_top_genes" in msg
    assert "n_genes_user" in msg
    # Migration hint should include the user's value
    assert "50" in msg


def test_gpu_path_also_rejects_old_kwargs(tmp_path):
    """The GPU function uses the same validator. Even on machines without
    CuPy, the validation runs BEFORE the cupy import (because it's in the
    function body before the import), so the TypeError should fire
    regardless of GPU availability."""
    from cosg._cytome_streaming import run_cosg_cytome_gpu
    with pytest.raises(TypeError) as exc:
        run_cosg_cytome_gpu(
            _stub_cytome_path(tmp_path),
            groupby="leiden",
            n_top_genes=50,
        )
    msg = str(exc.value)
    assert "n_top_genes" in msg and "n_genes_user" in msg


def test_unexpected_kwarg_raises_typeerror(tmp_path):
    """Any other stray kwarg should also fail (defensive)."""
    from cosg._cytome_streaming import run_cosg_cytome
    with pytest.raises(TypeError) as exc:
        run_cosg_cytome(
            _stub_cytome_path(tmp_path),
            groupby="leiden",
            totally_made_up_kwarg=42,
        )
    assert "totally_made_up_kwarg" in str(exc.value)


# =====================================================================
# Round 9 K-HARD: cytome_layer renamed to layer
# =====================================================================

def test_cytome_layer_kwarg_raises_typeerror_with_migration_hint(tmp_path):
    """Round 9 (2026-05-23): cytome_layer was renamed to layer for
    symmetry with the AnnData side. Passing cytome_layer= must fail
    loudly with a copy/pasteable migration hint."""
    from cosg._cytome_streaming import run_cosg_cytome
    with pytest.raises(TypeError) as exc:
        run_cosg_cytome(
            _stub_cytome_path(tmp_path),
            groupby="leiden",
            cytome_layer="infog",
        )
    msg = str(exc.value)
    assert "cytome_layer" in msg, (
        f"error msg should name 'cytome_layer': {msg!r}"
    )
    assert "layer" in msg, (
        f"error msg should mention the new kwarg name 'layer': {msg!r}"
    )
    # Migration hint should echo the user's value
    assert "infog" in msg, (
        f"error msg should echo the user's value 'infog': {msg!r}"
    )


# =====================================================================
# Round 9 K-HARD: function rename run_cosg_cytome_cpu -> run_cosg_cytome
# =====================================================================

def test_run_cosg_cytome_cpu_old_name_raises_with_migration_hint(tmp_path):
    """Round 9: the old function name survives as a stub that raises
    TypeError at call time with a clear migration message. Pin both the
    error type AND the migration string."""
    import cosg
    with pytest.raises(TypeError) as exc:
        cosg.run_cosg_cytome_cpu(
            _stub_cytome_path(tmp_path), groupby="leiden",
        )
    msg = str(exc.value)
    assert "run_cosg_cytome_cpu" in msg
    assert "run_cosg_cytome" in msg
    assert "Round 9" in msg or "renamed" in msg.lower()


def test_run_cosg_cytome_cpu_deep_import_also_raises(tmp_path):
    """The stub must catch ``from cosg._cytome_streaming import
    run_cosg_cytome_cpu`` too — not just ``cosg.run_cosg_cytome_cpu``.
    The stub-function approach (vs module __getattr__) makes this
    automatic; this test pins it so a future refactor to module
    __getattr__ doesn't regress the deep-import case."""
    from cosg._cytome_streaming import run_cosg_cytome_cpu
    with pytest.raises(TypeError) as exc:
        run_cosg_cytome_cpu(
            _stub_cytome_path(tmp_path), groupby="leiden",
        )
    msg = str(exc.value)
    assert "renamed" in msg.lower() and "run_cosg_cytome" in msg


@requires_piaso
def test_run_cosg_cytome_new_name_is_callable_default_cpu(tmp_path):
    """Sanity: the new name actually works — opening a tiny cytome and
    running with defaults produces a result. This pins the rename
    didn't break the function body."""
    import cytome
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    import cosg

    p = str(tmp_path / "tiny.cytome")
    ds = cytome.create(p)
    n_cells = 12
    leiden = np.array([f"g{i % 3}" for i in range(n_cells)])
    ds.set_entity("cells", pd.DataFrame({
        "cell_idx": np.arange(n_cells),
        "barcode": [f"AAA-{i}" for i in range(n_cells)],
        "Leiden": leiden,
    }))
    X = np.zeros((n_cells, 3), dtype=np.float32)
    X[leiden == "g0", 0] = 5.0
    X[leiden == "g1", 1] = 5.0
    X[leiden == "g2", 2] = 5.0
    ds.set_entity("genes", pd.DataFrame({
        "gene_idx": [0, 1, 2],
        "gene_id": ["A", "B", "C"],
    }))
    ds.add_matrix("RNA_counts", sp.csr_matrix(X))
    ds.flush()
    ds.close()

    result = cosg.run_cosg_cytome(
        p, groupby="Leiden", n_genes_user=3,
        remove_lowly_expressed=False, verbose=False,
    )
    assert "names" in result
    assert result["names"].shape == (3, 3)


def test_layer_kwarg_works_as_replacement(tmp_path):
    """The new `layer=` kwarg should accept what `cytome_layer=` used
    to accept. Smoke-test on the default 'counts'."""
    import cytome
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    import cosg

    p = str(tmp_path / "tiny.cytome")
    ds = cytome.create(p)
    n_cells = 12
    leiden = np.array([f"g{i % 3}" for i in range(n_cells)])
    ds.set_entity("cells", pd.DataFrame({
        "cell_idx": np.arange(n_cells),
        "barcode": [f"AAA-{i}" for i in range(n_cells)],
        "Leiden": leiden,
    }))
    X = np.zeros((n_cells, 3), dtype=np.float32)
    X[leiden == "g0", 0] = 5.0
    X[leiden == "g1", 1] = 5.0
    X[leiden == "g2", 2] = 5.0
    ds.set_entity("genes", pd.DataFrame({
        "gene_idx": [0, 1, 2],
        "gene_id": ["A", "B", "C"],
    }))
    ds.add_matrix("RNA_counts", sp.csr_matrix(X))
    ds.flush()
    ds.close()

    result = cosg.run_cosg_cytome(
        p, groupby="Leiden", n_genes_user=3, layer="counts",
        remove_lowly_expressed=False, verbose=False,
    )
    assert "names" in result


# =====================================================================
# Signature contract (anti-doc-drift)
# =====================================================================

def test_cpu_signature_matches_v2_api():
    """Round 9: signature contract for ``run_cosg_cytome`` (was
    ``run_cosg_cytome_cpu`` pre-Round-9)."""
    import inspect
    from cosg._cytome_streaming import run_cosg_cytome
    sig = inspect.signature(run_cosg_cytome)
    params = sig.parameters
    assert "n_genes_user" in params and params["n_genes_user"].default == 50
    assert "output_format" in params and params["output_format"].default == "ndarray"
    assert "iqr_normalize" in params and params["iqr_normalize"].default is None
    # Old kwargs must NOT appear as positional/keyword parameters
    assert "top_k" not in params
    assert "n_top_genes" not in params
    assert "cytome_layer" not in params, (
        "Round 9: cytome_layer renamed to layer; old name must not "
        "appear as a positional/keyword parameter."
    )
    # `layer` defaults to "auto" since Round 26: it resolves a
    # modality-appropriate normalization (log1p for RNA, TF-IDF for ATAC)
    # rather than one fixed default, because COSG scores a cosine
    # specificity and no single transform is right for both modalities.
    assert "layer" in params and params["layer"].default == "auto"


def test_gpu_signature_matches_v2_api():
    import inspect
    from cosg._cytome_streaming import run_cosg_cytome_gpu
    sig = inspect.signature(run_cosg_cytome_gpu)
    params = sig.parameters
    assert "n_genes_user" in params and params["n_genes_user"].default == 50
    assert "top_k" not in params
    assert "n_top_genes" not in params
    assert "cytome_layer" not in params
    # `layer` defaults to "auto" since Round 26: it resolves a
    # modality-appropriate normalization (log1p for RNA, TF-IDF for ATAC)
    # rather than one fixed default, because COSG scores a cosine
    # specificity and no single transform is right for both modalities.
    assert "layer" in params and params["layer"].default == "auto"


def test_run_cosg_cytome_cpu_stub_is_callable_attribute():
    """The Round 9 K-HARD stub must be a callable attribute of the
    package (so `cosg.run_cosg_cytome_cpu` AND
    `from cosg._cytome_streaming import run_cosg_cytome_cpu` both
    resolve). Calling it raises; just looking it up must succeed."""
    import cosg
    assert callable(cosg.run_cosg_cytome_cpu)
    from cosg._cytome_streaming import run_cosg_cytome_cpu
    assert callable(run_cosg_cytome_cpu)


# =====================================================================
# output_format validation
# =====================================================================

def test_invalid_output_format_raises_valueerror(tmp_path):
    from cosg._cytome_streaming import run_cosg_cytome
    with pytest.raises(ValueError) as exc:
        run_cosg_cytome(
            _stub_cytome_path(tmp_path),
            groupby="leiden",
            output_format="ndarry",   # typo
        )
    msg = str(exc.value)
    assert "ndarry" in msg
    assert "ndarray" in msg  # valid options listed
    assert "dict" in msg


# =====================================================================
# SS-B: iqr_normalize default behaviour (Round 7)
# =====================================================================

def test_iqr_normalize_default_for_ndarray_is_false():
    """Validator should resolve iqr_normalize=None → False for ndarray."""
    from cosg._cytome_streaming import _validate_cosg_api_v2
    use_sparse, iqr = _validate_cosg_api_v2(
        {}, output_format="ndarray", iqr_normalize=None, fn_name="t",
    )
    assert use_sparse is False
    assert iqr is False


def test_iqr_normalize_default_for_sparse_is_true():
    """Validator should resolve iqr_normalize=None → True for sparse outputs."""
    from cosg._cytome_streaming import _validate_cosg_api_v2
    for fmt in ("dict", "long", "dense"):
        use_sparse, iqr = _validate_cosg_api_v2(
            {}, output_format=fmt, iqr_normalize=None, fn_name="t",
        )
        assert use_sparse is True, f"{fmt} should be a sparse output"
        assert iqr is True, f"{fmt} should default iqr_normalize=True"


def test_iqr_normalize_explicit_overrides_default():
    """User can override the default in either direction."""
    from cosg._cytome_streaming import _validate_cosg_api_v2
    # ndarray + explicit True
    _, iqr = _validate_cosg_api_v2(
        {}, output_format="ndarray", iqr_normalize=True, fn_name="t",
    )
    assert iqr is True
    # dict + explicit False
    _, iqr = _validate_cosg_api_v2(
        {}, output_format="dict", iqr_normalize=False, fn_name="t",
    )
    assert iqr is False
