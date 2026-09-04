"""The significance layer on a cytome, streamed.

The analytic route is what makes this possible: the moments it needs are power
sums, which add across chunks, so the only new cost in the streaming pass is
two reductions. Sampled permutations would have needed the matrix B times.

What has to hold: the streamed p-values must agree with the in-memory ones on
the same data, and the scores must not move because significance was asked for.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

pytest.importorskip("cytome")
pytest.importorskip("anndata")


def _adata(n_cells=600, n_genes=120, seed=0):
    import anndata as ad
    rng = np.random.default_rng(seed)
    groups = np.array([f"ct{k}" for k in range(3)]).repeat(n_cells // 3)
    X = rng.poisson(0.4, size=(n_cells, n_genes)).astype(np.float64)
    blk = 5
    for k in range(3):
        m = groups == f"ct{k}"
        X[m, k * blk:(k + 1) * blk] += rng.poisson(6.0, size=(m.sum(), blk))
    a = ad.AnnData(
        X=sparse.csr_matrix(X),
        obs=pd.DataFrame({"cell_type": groups},
                         index=[f"c{i}" for i in range(n_cells)]),
        var=pd.DataFrame(index=[f"g{j}" for j in range(n_genes)]))
    return a


@pytest.fixture(scope="module")
def cytome_fixture(tmp_path_factory):
    import cytome
    a = _adata()
    out = tmp_path_factory.mktemp("cosgp") / "d.cytome"
    ds = cytome.from_anndata(a, modality="RNA", output=str(out))
    ds.cells["cell_type"] = a.obs.cell_type.astype(str).values
    ds.flush()
    ds.close()
    return a, str(out)


def test_streamed_pvalues_match_the_in_memory_ones(cytome_fixture):
    """Same data, two engines, one answer."""
    from cosg._cytome_streaming import run_cosg_cytome
    from cosg._pvalues import analytic_pvalues

    a, path = cytome_fixture
    res = run_cosg_cytome(path, groupby="cell_type", n_genes_user=10,
                          layer="counts", calculate_pvalues=True,
                          pvalue_method="normal", verbose=False)
    assert "pvals" in res, f"streaming returned {sorted(res)}"

    order = list(res["groups_order"])
    P_mem, _, _ = analytic_pvalues(a.X, a.obs.cell_type.values, order,
                                   method="normal")
    names = np.asarray(res["names"])
    gene_ix = {g: i for i, g in enumerate(a.var_names)}
    streamed = np.asarray(res["pvals"], dtype=float)

    for j, grp in enumerate(order):
        for r in range(names.shape[0]):
            gi = gene_ix[names[r, j]]
            assert np.isclose(streamed[r, j], P_mem[gi, j], rtol=1e-6), (
                f"{grp}/{names[r, j]}: streamed {streamed[r, j]:.3e} "
                f"vs in-memory {P_mem[gi, j]:.3e}")


def test_scores_are_untouched_by_asking_for_pvalues(cytome_fixture):
    """The no-regression contract, on the streaming path."""
    from cosg._cytome_streaming import run_cosg_cytome

    _, path = cytome_fixture
    base = run_cosg_cytome(path, groupby="cell_type", n_genes_user=10,
                           layer="counts", verbose=False)
    withp = run_cosg_cytome(path, groupby="cell_type", n_genes_user=10,
                            layer="counts", calculate_pvalues=True,
                            pvalue_method="normal", verbose=False)
    assert np.array_equal(np.asarray(base["names"]), np.asarray(withp["names"]))
    assert np.allclose(np.asarray(base["scores"], dtype=float),
                       np.asarray(withp["scores"], dtype=float))


def test_planted_markers_are_significant_after_fdr(cytome_fixture):
    from cosg._cytome_streaming import run_cosg_cytome

    _, path = cytome_fixture
    res = run_cosg_cytome(path, groupby="cell_type", n_genes_user=8,
                          layer="counts", calculate_pvalues=True,
                          verbose=False)
    q = np.asarray(res["pvals_adj"], dtype=float)
    names = np.asarray(res["names"])
    order = list(res["groups_order"])
    for j, grp in enumerate(order):
        k = int(grp[-1])
        planted = {f"g{i}" for i in range(k * 5, (k + 1) * 5)}
        called = {names[r, j] for r in range(names.shape[0])
                  if np.isfinite(q[r, j]) and q[r, j] < 0.01}
        assert planted <= called, f"{grp}: missed {planted - called}"


def test_batch_key_with_pvalues_says_so_rather_than_being_wrong(cytome_fixture):
    """The stratified null is not implemented here; refusing beats guessing."""
    from cosg._cytome_streaming import run_cosg_cytome

    _, path = cytome_fixture
    with pytest.raises(NotImplementedError, match="batch"):
        run_cosg_cytome(path, groupby="cell_type", n_genes_user=5,
                        layer="counts", calculate_pvalues=True,
                        batch_key="cell_type", verbose=False)
