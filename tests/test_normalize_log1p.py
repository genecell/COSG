"""COSG's native log1p must match the PIASO implementation it replaced.

COSG's default ``run_cosg_cytome(layer='log1p')`` used to import
``piaso.preprocessing._normalize_log1p``. PIASO requires ``cosg>=1.0.3``, so
COSG's headline streaming feature depended on a package that depends on COSG:
the import was lazy enough to work, but it meant the default path could not run
without installing PIASO, and it made the two mutually dependent in practice.

log1p is ``log1p(counts / cell_depth * scale)`` -- the field's most standard
normalization and no part of PIASO's contribution -- so COSG owns it now.

The point of this file is that "owns it now" must not mean "quietly computes
something slightly different". Existing results were produced through PIASO's
version; anyone re-running must get the same numbers.
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse as sp

from cosg._normalize import normalize_chunk_log1p


def _chunk(seed=0, n_cells=40, n_genes=25, density=0.3):
    rng = np.random.default_rng(seed)
    m = sp.random(n_cells, n_genes, density=density, format="csr",
                  random_state=rng, data_rvs=lambda k: rng.integers(1, 50, k))
    return m.astype(np.float32)


def test_matches_piaso_implementation():
    """Bit-level agreement with the function this replaced.

    Skipped when PIASO is absent -- which is the normal case in COSG's CI, and
    exactly why COSG should not have depended on it.
    """
    piaso_impl = pytest.importorskip(
        "piaso.preprocessing._normalize_log1p",
        reason="PIASO not installed; nothing to compare against",
    )._normalize_chunk_log1p

    for seed in range(5):
        chunk = _chunk(seed=seed)
        depth = np.asarray(chunk.sum(axis=1)).ravel()
        ours = normalize_chunk_log1p(chunk, depth, 1e4)
        theirs = piaso_impl(chunk, depth, 1e4)
        np.testing.assert_array_equal(
            ours.toarray(), theirs.toarray(),
            err_msg=f"seed={seed}: COSG's log1p diverged from PIASO's",
        )


def test_sparse_in_sparse_out():
    """COSG chains the result into a sparse matmul and ``.data ** 2``.

    Returning dense here does not fail loudly -- it fails later inside the
    scoring kernel with a confusing error about ``memoryview``. This has bitten
    the streaming path before.
    """
    chunk = _chunk()
    out = normalize_chunk_log1p(chunk, np.asarray(chunk.sum(axis=1)).ravel(), 1e4)
    assert sp.issparse(out), "log1p densified the chunk"
    assert sp.isspmatrix_csr(out)
    _ = out.data ** 2                      # the accumulator that breaks on dense
    assert out.nnz == chunk.nnz, "sparsity structure changed; log1p(0) must stay 0"


def test_rows_scale_to_the_target_before_the_log():
    chunk = _chunk()
    depth = np.asarray(chunk.sum(axis=1)).ravel()
    out = normalize_chunk_log1p(chunk, depth, 1e4)
    # Undo the log; every non-empty row should now sum to the scale factor.
    recovered = out.copy()
    recovered.data = np.expm1(recovered.data)
    sums = np.asarray(recovered.sum(axis=1)).ravel()
    np.testing.assert_allclose(sums[depth > 0], 1e4, rtol=1e-3)


def test_empty_cell_gives_zeros_not_nan():
    """A zero-depth cell divides by 1.0, not by 0.

    QC-filtered inputs should not contain empty cells, but a subset or a
    modality with no counts for some barcodes can produce one, and NaNs
    propagate silently into marker scores.
    """
    chunk = sp.csr_matrix(np.array([[0, 0, 0], [1, 2, 3]], dtype=np.float32))
    out = normalize_chunk_log1p(chunk, np.array([0.0, 6.0]), 1e4)
    assert not np.isnan(out.toarray()).any()
    assert out.toarray()[0].sum() == 0.0


def test_dense_input_stays_dense():
    chunk = _chunk().toarray()
    out = normalize_chunk_log1p(chunk, chunk.sum(axis=1), 1e4)
    assert isinstance(out, np.ndarray) and not sp.issparse(out)
    assert not np.isnan(out).any()


def test_default_layer_does_not_import_piaso():
    """The dependency-direction guarantee, asserted rather than assumed.

    If a future edit routes log1p back through PIASO, the mutual dependency
    returns silently -- every dev machine here has PIASO installed.
    """
    import ast
    import pathlib

    import cosg

    src = pathlib.Path(cosg.__file__).parent / "_cytome_streaming.py"
    tree = ast.parse(src.read_text())

    for node in ast.walk(tree):
        if not isinstance(node, ast.FunctionDef):
            continue
        if node.name != "_resolve_chunk_normalizer":
            continue
        # Collect the piaso imports and the layer names guarding them.
        for sub in ast.walk(node):
            mod = None
            if isinstance(sub, ast.ImportFrom) and sub.module:
                mod = sub.module
            elif isinstance(sub, ast.Import):
                mod = sub.names[0].name
            if mod and mod.split(".")[0] == "piaso":
                seg = ast.get_source_segment(src.read_text(), node) or ""
                idx = seg.find(mod)
                before = seg[:idx]
                assert "infog" in before or "tfidf" in before, (
                    f"a piaso import ({mod}) is reachable without first "
                    "branching on infog/tfidf -- the log1p and counts paths "
                    "must stay PIASO-free"
                )
        break
    else:
        pytest.fail("_resolve_chunk_normalizer not found")
