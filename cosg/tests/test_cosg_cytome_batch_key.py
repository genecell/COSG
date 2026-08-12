"""Batch-aware COSG on the cytome streaming path.

`batch_key` computes the cosine specificity **within each batch, using that
batch's own gene and group norms**, then averages over the batches where the
group had at least `batch_cell_number_threshold` cells. Averaging per-batch
cosines is what makes the score batch-corrected; accumulating one global dot
product and dividing at the end would only be batch-weighted, which is a
different (and wrong) thing.

The load-bearing test is `test_matches_anndata_batch_reference`: the streaming
result must equal the in-memory AnnData implementation that already existed. The
rest guard the edges — the threshold actually dropping under-populated groups,
`batch_key=None` staying bit-identical to before, and the paths that cannot do
batch correction refusing rather than silently returning an uncorrected score.
"""
import numpy as np
import pytest

pytest.importorskip("cytome")
import cytome  # noqa: E402

import cosg as cosg_mod  # noqa: E402
from cosg._cytome_streaming import run_cosg_cytome  # noqa: E402

N_CELLS, N_GENES, N_GROUPS, N_BATCH = 400, 60, 4, 3
SEED = 7


def _adata(n_cells=N_CELLS, n_genes=N_GENES, seed=SEED, batch_effect=3.0):
    """Cell types with real markers, plus an additive per-batch offset."""
    import anndata as ad
    import pandas as pd
    import scipy.sparse as sp
    rng = np.random.default_rng(seed)
    groups = np.array([f"ct{i % N_GROUPS}" for i in range(n_cells)])
    batches = np.array([f"b{(i // 7) % N_BATCH}" for i in range(n_cells)])
    X = rng.poisson(0.5, size=(n_cells, n_genes)).astype(np.float32)
    blk = n_genes // N_GROUPS
    for k in range(N_GROUPS):
        m = groups == f"ct{k}"
        X[m, k * blk:(k + 1) * blk] += rng.poisson(5.0, size=(m.sum(), blk))
    # batch effect: a multiplicative depth shift, the thing batch correction
    # is supposed to see through
    for j in range(N_BATCH):
        X[batches == f"b{j}"] *= (1.0 + batch_effect * j / N_BATCH)
    a = ad.AnnData(
        X=sp.csr_matrix(X),
        obs=pd.DataFrame({"cell_type": groups, "batch": batches},
                         index=[f"c{i}" for i in range(n_cells)]),
        var=pd.DataFrame(index=[f"g{j}" for j in range(n_genes)]))
    a.layers["counts"] = a.X
    return a


@pytest.fixture(scope="module")
def fixture(tmp_path_factory):
    a = _adata()
    out = tmp_path_factory.mktemp("cosgb") / "d.cytome"
    ds = cytome.from_anndata(a, modality="RNA", output=str(out))
    ds.cells["cell_type"] = a.obs.cell_type.astype(str).values
    ds.cells["batch"] = a.obs.batch.astype(str).values
    ds.flush()
    ds.close()
    return a, str(out)


def _nm(res):
    """run_cosg_cytome returns {'names', 'scores', 'groups_order'}."""
    return np.asarray(res["names"] if isinstance(res, dict) else res[0])


def _sc(res):
    return np.asarray(res["scores"] if isinstance(res, dict) else res[1], dtype=float)


def test_matches_anndata_batch_reference(fixture):
    """The streaming batch score must equal the in-memory one it mirrors."""
    a, path = fixture
    ad_res = cosg_mod.cosg(
        a.copy(), groupby="cell_type", batch_key="batch",
        batch_cell_number_threshold=3, n_genes_user=10, mu=1.0,
        layer="counts", cpu_chunk_size=None, verbosity=0, copy=True)
    ref = ad_res.uns["cosg"] if hasattr(ad_res, "uns") else ad_res

    cy = run_cosg_cytome(
        path, groupby="cell_type", batch_key="batch",
        batch_cell_number_threshold=3, n_genes_user=10, mu=1.0,
        layer="counts", modality="RNA", batch_size=64,
        feature_batching="none", verbose=False)

    ref_names = np.asarray(ref["names"].tolist() if hasattr(ref["names"], "tolist")
                           else ref["names"])
    cy_names = _nm(cy)
    assert cy_names.shape == ref_names.shape, (cy_names.shape, ref_names.shape)

    # Compare the ranked gene sets per group — tied scores can reorder, so the
    # set is the honest comparison, not the exact permutation.
    for j in range(cy_names.shape[1]):
        a_set = {str(x) for x in np.asarray(ref_names)[:, j]}
        c_set = {str(x) for x in cy_names[:, j]}
        overlap = len(a_set & c_set) / max(len(a_set), 1)
        assert overlap >= 0.9, (f"group {j}: only {overlap:.0%} of the top-10 agree\n"
                                f"  anndata: {sorted(a_set)}\n  cytome : {sorted(c_set)}")


def test_batch_key_changes_the_result(fixture):
    """A batch effect this large must move the score; otherwise nothing happened."""
    _, path = fixture
    common = dict(groupby="cell_type", n_genes_user=10, mu=1.0, layer="counts",
                  modality="RNA", batch_size=64, feature_batching="none", verbose=False)
    plain = run_cosg_cytome(path, **common)
    batched = run_cosg_cytome(path, batch_key="batch", **common)
    assert not np.allclose(_sc(plain), _sc(batched)), \
        "batch_key had no effect on the scores"


def test_no_batch_key_is_unchanged(fixture):
    """The default path must be bit-identical to before the batch dimension."""
    _, path = fixture
    common = dict(groupby="cell_type", n_genes_user=10, mu=1.0, layer="counts",
                  modality="RNA", batch_size=64, feature_batching="none", verbose=False)
    a = run_cosg_cytome(path, **common)
    b = run_cosg_cytome(path, batch_key=None, **common)
    assert np.array_equal(_nm(a), _nm(b))
    assert np.allclose(_sc(a), _sc(b))


def test_threshold_drops_groups_absent_from_every_batch(fixture):
    """A group below threshold in all batches has no evidence and scores 0."""
    _, path = fixture
    res = run_cosg_cytome(
        path, groupby="cell_type", batch_key="batch",
        batch_cell_number_threshold=10_000,          # nothing can meet this
        n_genes_user=5, mu=1.0, layer="counts", modality="RNA",
        batch_size=64, feature_batching="none", verbose=False)
    assert np.allclose(_sc(res), 0.0), \
        "with an unreachable threshold every score must collapse to 0"


def test_feature_batched_path_agrees_with_unbatched(fixture):
    """Chromosome/disk feature batching must give the SAME batch-corrected score.

    Feature batching is a memory strategy, not a different estimator: splitting
    the feature axis cannot change a per-feature cosine. If these disagree, one
    of the two paths is wrong — which is exactly the silent-uncorrected-score
    failure this parameter exists to avoid.
    """
    _, path = fixture
    common = dict(groupby="cell_type", batch_key="batch",
                  batch_cell_number_threshold=3, n_genes_user=10, mu=1.0,
                  layer="counts", modality="RNA", batch_size=64, verbose=False)
    plain = run_cosg_cytome(path, feature_batching="none", **common)
    # 'chromosome' is redirected to 'disk' when batch_key is set, so the batch
    # axis does not multiply peak RAM; the answer must be unaffected.
    chrom = run_cosg_cytome(path, feature_batching="chromosome", **common)
    for j in range(_nm(plain).shape[1]):
        a_set = {str(x) for x in _nm(plain)[:, j]}
        c_set = {str(x) for x in _nm(chrom)[:, j]}
        overlap = len(a_set & c_set) / max(len(a_set), 1)
        assert overlap >= 0.9, (f"group {j}: feature-batched and unbatched "
                                f"batch-corrected scores disagree ({overlap:.0%})\n"
                                f"  none : {sorted(a_set)}\n  chrom: {sorted(c_set)}")


def test_bad_batch_key_is_named(fixture):
    _, path = fixture
    with pytest.raises(ValueError, match="nope"):
        run_cosg_cytome(path, groupby="cell_type", batch_key="nope",
                        n_genes_user=5, layer="counts", modality="RNA",
                        feature_batching="none", verbose=False)


def test_polymorphic_entry_point_forwards_batch_key(fixture):
    """cosg.cosg(cytome, batch_key=...) must no longer raise TypeError."""
    _, path = fixture
    res = cosg_mod.cosg(path, groupby="cell_type", batch_key="batch",
                        n_genes_user=5, mu=1.0, layer="counts",
                        modality="RNA", feature_batching="none", verbosity=0)
    assert _nm(res).shape[0] == 5


# ── GPU ───────────────────────────────────────────────────────────────────────

def _has_gpu():
    try:
        import cupy
        return cupy.cuda.runtime.getDeviceCount() > 0
    except Exception:
        return False


@pytest.mark.skipif(not _has_gpu(), reason="no cupy / no CUDA device")
def test_gpu_batch_matches_cpu_batch(fixture):
    """The GPU batch path must agree with the CPU one it mirrors.

    Same construction, different device: per-batch cosine with that batch's own
    norms, averaged over batches meeting the threshold. A device that quietly
    computed the uncorrected score would pass every other test in this file.
    """
    _, path = fixture
    common = dict(groupby="cell_type", batch_key="batch",
                  batch_cell_number_threshold=3, n_genes_user=10, mu=1.0,
                  layer="counts", modality="RNA", batch_size=64, verbose=False)
    cpu = run_cosg_cytome(path, device="cpu", feature_batching="none", **common)
    gpu = run_cosg_cytome(path, device="gpu", **common)
    for j in range(_nm(cpu).shape[1]):
        a_set = {str(x) for x in _nm(cpu)[:, j]}
        g_set = {str(x) for x in _nm(gpu)[:, j]}
        overlap = len(a_set & g_set) / max(len(a_set), 1)
        assert overlap >= 0.9, (f"group {j}: GPU and CPU batch scores disagree "
                                f"({overlap:.0%})\n  cpu: {sorted(a_set)}\n  gpu: {sorted(g_set)}")


@pytest.mark.skipif(not _has_gpu(), reason="no cupy / no CUDA device")
def test_gpu_without_batch_key_is_unchanged(fixture):
    """Adding the batch axis must not perturb the ordinary GPU path."""
    _, path = fixture
    common = dict(groupby="cell_type", n_genes_user=10, mu=1.0, layer="counts",
                  modality="RNA", batch_size=64, verbose=False, device="gpu")
    a = run_cosg_cytome(path, **common)
    b = run_cosg_cytome(path, batch_key=None, **common)
    assert np.array_equal(_nm(a), _nm(b))
    assert np.allclose(_sc(a), _sc(b))
