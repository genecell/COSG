"""Round-trip equivalence test for run_cosg_cytome(output_format='dense').

Validates that the cytome-streaming dense output matches the AnnData
reference chain (cosg.cosg → cosg.indexByGene → cosg.iqrLogNormalize)
on the same underlying counts matrix.

Important caveat: the two paths apply the all-values IQR computation to
slightly different inputs (cytome computes IQR over full COSG λ column
including zero-clamped negatives; AnnData chain computes IQR over the
top-N-per-group dense matrix with -1 → 0 fill). They should produce
strongly correlated dense matrices, but not bit-identity. We assert per
-column Pearson r ≥ 0.95 and shape match.

The 'long' output format is also covered: assert that pivoting it back
to dense via pandas matches the dict-based dense.

Round 7 (2026-05-22): updated to the new API — ``top_k`` is gone, use
``n_genes_user=N, output_format='dict'|'long'|'dense'``.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp


# Local-only fixture: a real MTG subset that is not distributed with the
# package. Point COSG_TEST_MTG_H5AD at your own .h5ad to run this comparison;
# without it the test skips. It used to hard-code one machine's absolute path,
# which leaked that path into the published repository and could never resolve
# for anyone else.
MTG_H5AD = os.environ.get("COSG_TEST_MTG_H5AD", "")


@pytest.fixture(scope="module")
def mtg_anndata_and_cytome(tmp_path_factory):
    """Build matched AnnData + cytome from the same MTG raw counts."""
    import anndata
    import cytome
    if not MTG_H5AD or not os.path.exists(MTG_H5AD):
        pytest.skip("set COSG_TEST_MTG_H5AD to a local MTG .h5ad to run this")
    src = anndata.read_h5ad(MTG_H5AD)
    raw = src.raw.X
    a = anndata.AnnData(
        X=raw.copy() if sp.issparse(raw) else raw.copy(),
        obs=src.obs.copy(),
        var=pd.DataFrame(index=list(src.raw.var_names)),
    )
    a.obs_names_make_unique()
    a.var_names_make_unique()
    cytome_path = tmp_path_factory.mktemp("dense_rt") / "mtg.cytome"
    ds = cytome.from_anndata(a, modality="RNA", output=str(cytome_path))
    ds.close()
    return a, str(cytome_path)


def _anndata_reference_dense(a, groupby, n_genes_user, mu, q_upper, q_lower):
    """Reproduce the AnnData reference chain:
        cosg.cosg → indexByGene(set_nan_to_zero, convert_negative_one_to_zero)
                 → iqrLogNormalize.
    Returns a (genes × groups) DataFrame.
    """
    import cosg
    a_local = a.copy()
    cosg.cosg(
        a_local, groupby=groupby, n_genes_user=n_genes_user, mu=mu,
        remove_lowly_expressed=True, expressed_pct=0.1,
    )
    ref = cosg.indexByGene(
        a_local.uns["cosg"]["COSG"],
        set_nan_to_zero=True,
        convert_negative_one_to_zero=True,
    )
    ref = cosg.iqrLogNormalize(ref, q_upper=q_upper, q_lower=q_lower)
    return ref


def test_cytome_dense_format_round_trip(mtg_anndata_and_cytome):
    """run_cosg_cytome(n_genes_user=n_features, output_format='dense')
    should match cosg.cosg(adata) + indexByGene + iqrLogNormalize on the
    same data.
    """
    from cosg._cytome_streaming import run_cosg_cytome

    a, cytome_path = mtg_anndata_and_cytome
    n_features = a.n_vars
    groupby = "CrossArea_subclass"
    mu = 1.0
    q_upper, q_lower = 0.95, 0.75

    # AnnData reference
    ref = _anndata_reference_dense(
        a, groupby, n_genes_user=n_features, mu=mu,
        q_upper=q_upper, q_lower=q_lower,
    )

    # Cytome dense output
    result = run_cosg_cytome(
        cytome_path,
        groupby=groupby,
        modality="RNA",
        n_genes_user=n_features,
        mu=mu,
        remove_lowly_expressed=True,
        expressed_pct=0.1,
        q_upper=q_upper,
        q_lower=q_lower,
        output_format="dense",
        verbose=False,
    )
    cy = result["scores_df"]

    # ---- Shape ----
    assert cy.shape[1] == ref.shape[1], (
        f"n_groups mismatch: cytome={cy.shape[1]}, anndata={ref.shape[1]}"
    )
    # cytome dense uses ALL features in the cytome's gene table; the AnnData
    # reference includes only genes that were top-N for at least one group.
    # So cy.shape[0] >= ref.shape[0] is the only guarantee.
    assert cy.shape[0] >= ref.shape[0], (
        f"cytome dense should index all features; cy={cy.shape[0]}, ref={ref.shape[0]}"
    )

    # ---- Column alignment ----
    common_groups = [g for g in cy.columns if g in ref.columns]
    assert len(common_groups) >= 1, "No overlapping groups between cy and ref"
    # Restrict to features present in both index sets
    common_features = ref.index.intersection(cy.index)
    assert len(common_features) >= 100, (
        f"Too few common features: {len(common_features)}"
    )

    cy_aligned = cy.loc[common_features, common_groups]
    ref_aligned = ref.loc[common_features, common_groups]

    # ---- Per-column Pearson r ----
    col_corrs = []
    for g in common_groups:
        cv = cy_aligned[g].values.astype(np.float64)
        rv = ref_aligned[g].values.astype(np.float64)
        # Handle constant columns gracefully
        if cv.std() == 0 or rv.std() == 0:
            continue
        r = np.corrcoef(cv, rv)[0, 1]
        col_corrs.append((g, r))
    if not col_corrs:
        pytest.fail("Could not compute any column correlations.")
    median_r = float(np.median([r for _, r in col_corrs]))
    min_r = float(np.min([r for _, r in col_corrs]))
    print(
        f"\n[round-trip] per-column Pearson r: "
        f"median={median_r:.4f}, min={min_r:.4f} "
        f"({len(col_corrs)} groups compared)"
    )
    if min_r < 0.85:
        worst = sorted(col_corrs, key=lambda x: x[1])[:3]
        print(f"  Worst groups: {worst}")
    assert median_r >= 0.95, (
        f"Per-column Pearson r too low: median={median_r:.4f}. "
        f"Cytome dense and AnnData reference should be strongly correlated."
    )


def test_cytome_long_format_pivots_to_dense(mtg_anndata_and_cytome):
    """The 'long' DataFrame should pivot back to the same numeric values
    as the 'dense' output (when restricted to the long-format's feature set).
    """
    from cosg._cytome_streaming import run_cosg_cytome

    _, cytome_path = mtg_anndata_and_cytome
    groupby = "CrossArea_subclass"

    long_result = run_cosg_cytome(
        cytome_path, groupby=groupby, modality="RNA",
        n_genes_user=1000, mu=1.0,
        remove_lowly_expressed=True, expressed_pct=0.1,
        output_format="long",
        verbose=False,
    )
    dense_result = run_cosg_cytome(
        cytome_path, groupby=groupby, modality="RNA",
        n_genes_user=1000, mu=1.0,
        remove_lowly_expressed=True, expressed_pct=0.1,
        output_format="dense",
        verbose=False,
    )

    long_df = long_result["scores_df"]
    dense_df = dense_result["scores_df"]

    # Pivot the long DataFrame back to dense
    pivoted = long_df.pivot_table(
        index="feature", columns="group", values="score", aggfunc="first",
    ).fillna(0.0)
    # Reindex to dense_df's feature/group order
    pivoted = pivoted.reindex(
        index=dense_df.index, columns=dense_df.columns, fill_value=0.0,
    )

    # Bit-equality (or near-bit) — both come from the same underlying dict
    assert pivoted.shape == dense_df.shape
    np.testing.assert_allclose(
        pivoted.values.astype(np.float32),
        dense_df.values.astype(np.float32),
        rtol=0, atol=1e-6,
        err_msg="long and dense formats should agree on per-(feature, group) score.",
    )


def test_cytome_iqr_normalize_standard_mode_no_op_check(mtg_anndata_and_cytome):
    """When iqr_normalize=True and output_format='ndarray', scores in the
    standard output should differ from iqr_normalize=False (i.e., the
    kwarg actually does something now). Pre-fix this was a silent no-op."""
    from cosg._cytome_streaming import run_cosg_cytome

    _, cytome_path = mtg_anndata_and_cytome
    groupby = "CrossArea_subclass"

    raw_result = run_cosg_cytome(
        cytome_path, groupby=groupby, modality="RNA",
        n_genes_user=20, mu=1.0,
        remove_lowly_expressed=True, expressed_pct=0.1,
        iqr_normalize=False, verbose=False,
    )
    norm_result = run_cosg_cytome(
        cytome_path, groupby=groupby, modality="RNA",
        n_genes_user=20, mu=1.0,
        remove_lowly_expressed=True, expressed_pct=0.1,
        iqr_normalize=True, verbose=False,
    )

    # Names should match (same top-N per group)
    np.testing.assert_array_equal(raw_result["names"], norm_result["names"])

    # Scores should differ: norm = log1p(raw / iqr) ≠ raw for non-trivial iqr
    raw_scores = np.asarray(raw_result["scores"])
    norm_scores = np.asarray(norm_result["scores"])
    assert raw_scores.shape == norm_scores.shape
    diff = float(np.linalg.norm(raw_scores - norm_scores))
    base = float(np.linalg.norm(raw_scores))
    rel = diff / base if base > 0 else float("inf")
    print(
        f"\n[iqr_normalize] ||raw - norm|| / ||raw|| = {rel:.4f} "
        f"(raw mean={raw_scores.mean():.4f}, norm mean={norm_scores.mean():.4f})"
    )
    assert rel > 0.05, (
        f"iqr_normalize=True should produce different scores than False, "
        f"but rel diff = {rel:.4f}. Suggests iqr_normalize is still a no-op."
    )
    # All normalized scores should be ≥ 0 (log1p(non-negative))
    assert (norm_scores >= 0).all(), "Normalized scores should all be ≥ 0."
