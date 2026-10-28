"""`neg_log10_pvals_adj`: a volcano axis that is safe to threshold.

`neg_log10_pvals` is -log10 of the RAW p. It exists because a third of a clean
cell type's top markers reach the float64 floor and `pvals` alone cannot rank
them -- but a line at 1.3 on that axis is raw 0.05, not FDR 0.05. The adjusted
twin is computed with BH in log space, so it neither underflows nor lies about
what it is.
"""

import numpy as np
import pandas as pd
import pytest
import anndata as ad

import cosg
from cosg._pvalues import adjust_fdr, adjust_fdr_log


def _fixture(n_cells=300, n_genes=60, seed=0):
    rng = np.random.default_rng(seed)
    X = rng.poisson(1.0, size=(n_cells, n_genes)).astype(np.float32)
    labels = np.array(["A", "B", "C"] * (n_cells // 3))
    # give each group a block of genuinely specific genes
    for gi, g in enumerate(["A", "B", "C"]):
        X[labels == g, gi * 5:(gi + 1) * 5] += 40
    a = ad.AnnData(X=X)
    a.var_names = [f"G{i}" for i in range(n_genes)]
    a.obs["ct"] = pd.Categorical(labels)
    return a


# ------------------------------------------------------- the BH identity --


def test_log_space_bh_equals_linear_bh_where_linear_works():
    rng = np.random.default_rng(0)
    p = rng.random(500) ** 3
    np.testing.assert_allclose(np.exp(adjust_fdr_log(np.log(p))),
                               adjust_fdr(p), rtol=1e-12)


def test_log_space_bh_still_ranks_where_linear_underflows():
    lp = np.array([-800.0, -750.0, -3.0, -0.1])       # exp() is 0, 0, ., .
    assert (adjust_fdr(np.exp(lp))[:2] == 0).all()     # linear cannot rank these
    ladj = adjust_fdr_log(lp)
    assert ladj[0] < ladj[1] < ladj[2] < ladj[3]      # log space still can


def test_log_space_bh_is_monotone_and_capped_at_one():
    rng = np.random.default_rng(3)
    lp = np.sort(-np.abs(rng.normal(size=200) * 50))
    ladj = adjust_fdr_log(lp, "fdr_bh")
    assert np.all(np.diff(ladj) >= -1e-12)            # ordered as the input is
    assert ladj.max() <= 0.0                          # log(1)
    assert np.isfinite(ladj).all()


def test_fdr_by_is_more_conservative_than_fdr_bh():
    rng = np.random.default_rng(4)
    lp = np.log(rng.random(100) ** 2)
    assert (adjust_fdr_log(lp, "fdr_by") >= adjust_fdr_log(lp, "fdr_bh") - 1e-12).all()


def test_nan_positions_are_preserved():
    lp = np.array([-5.0, np.nan, -1.0, np.nan])
    out = adjust_fdr_log(lp)
    assert np.isnan(out[[1, 3]]).all() and np.isfinite(out[[0, 2]]).all()


# ------------------------------------------------------ the COSG column --


def test_column_is_present_and_matches_pvals_adj():
    a = _fixture()
    cosg.cosg(a, groupby="ct", key_added="cosg", n_genes_user=15,
              calculate_pvalues=True)
    res = a.uns["cosg"]
    assert "neg_log10_pvals_adj" in res
    adj = np.asarray(pd.DataFrame(res["pvals_adj"]).values, dtype=float)
    nl = np.asarray(pd.DataFrame(res["neg_log10_pvals_adj"]).values, dtype=float)
    ok = np.isfinite(adj) & (adj > 1e-300) & np.isfinite(nl)
    assert ok.sum() > 0
    np.testing.assert_allclose(nl[ok], -np.log10(adj[ok]), rtol=1e-6)


def test_adjusted_axis_is_never_above_the_raw_one():
    """Adjusting can only make a p larger, so -log10 can only get smaller."""
    a = _fixture()
    cosg.cosg(a, groupby="ct", key_added="cosg", n_genes_user=20,
              calculate_pvalues=True)
    raw = np.asarray(pd.DataFrame(a.uns["cosg"]["neg_log10_pvals"]).values, dtype=float)
    adj = np.asarray(pd.DataFrame(a.uns["cosg"]["neg_log10_pvals_adj"]).values, dtype=float)
    ok = np.isfinite(raw) & np.isfinite(adj)
    assert (adj[ok] <= raw[ok] + 1e-9).all()


def test_absent_when_pvalues_are_not_requested():
    a = _fixture()
    cosg.cosg(a, groupby="ct", key_added="cosg", n_genes_user=10)
    assert "neg_log10_pvals_adj" not in a.uns["cosg"]
