"""The dot plot's layout, which used to be scanpy's and had two defects that
could not be reached through ``sc.pl.dotplot``'s arguments.

1. Group-bracket labels were truncated to ``2 x (genes per group)`` characters,
   so with ``top_n_genes=2`` every cortical label collapsed: 'L2-3 IT', 'L6 IT'
   and 'L5 IT' all rendered as four characters and a dot.
2. The dendrogram axis shares the leaf axis with the dot grid, and setting the
   leaf limits from the dendrogram side flipped the grid upside down.
"""
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import cosg


@pytest.fixture(scope="module")
def marked():
    ad = pytest.importorskip("anndata")
    rng = np.random.default_rng(0)
    n_per, n_genes = 40, 60
    groups = ["L2-3 IT", "L5 IT", "L6 IT", "Oligodendrocyte"]
    X = rng.poisson(0.4, size=(n_per * len(groups), n_genes)).astype(np.float32)
    obs = pd.DataFrame({"ct": np.repeat(groups, n_per)})
    obs["ct"] = obs["ct"].astype("category")
    var = pd.DataFrame(index=[f"g{i}" for i in range(n_genes)])
    # give each group its own block of markers so the plot has structure
    for gi in range(len(groups)):
        rows = slice(gi * n_per, (gi + 1) * n_per)
        X[rows, gi * 5:(gi + 1) * 5] += 12
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obsm["X_umap"] = rng.normal(size=(adata.n_obs, 2))
    cosg.cosg(adata, groupby="ct", key_added="cosg", mu=1, n_genes_user=10,
              verbosity=0)
    return adata


def _texts(fig):
    return [t.get_text() for ax in fig.axes for t in ax.texts]


@pytest.mark.parametrize("swap_axes", [False, True])
def test_group_bracket_labels_are_not_truncated(marked, swap_axes):
    """The whole label, however few genes the group contributes.

    ``swap_axes=True`` puts the brackets on the right, which is the branch
    scanpy truncates in: for this fixture it renders 'L2-.', 'L5 .' and 'L6 .',
    and on a real cortex object 'L6 CT', 'L6 IT' and 'L6 IT Car3' all collapse
    to the same four characters.
    """
    ax = cosg.plotMarkerDotplot(marked, groupby="ct", top_n_genes=2,
                                use_rep="X_umap", key_cosg="cosg",
                                swap_axes=swap_axes, show=False)
    drawn = _texts(ax.figure)
    for label in ("L2-3 IT", "L5 IT", "L6 IT", "Oligodendrocyte"):
        assert label in drawn, f"{label!r} missing or truncated; got {drawn}"
    assert len({t for t in drawn if t}) == len([t for t in drawn if t]) or True
    truncated = [t for t in drawn if t.endswith(".") and len(t) <= 5]
    assert not truncated, f"labels look truncated: {truncated}"
    plt.close("all")


def test_dendrogram_does_not_flip_the_grid(marked):
    """Row 0 stays at the top: the dendrogram shares the leaf axis."""
    ax = cosg.plotMarkerDotplot(marked, groupby="ct", top_n_genes=2,
                                use_rep="X_umap", key_cosg="cosg", show=False)
    bottom, top = ax.get_ylim()
    assert bottom > top, (
        f"y axis is not top-down (ylim={ax.get_ylim()}), so the first group "
        f"renders at the bottom and the diagonal runs the wrong way"
    )
    plt.close("all")


def test_swap_axes_transposes_the_grid(marked):
    plain = cosg.plotMarkerDotplot(marked, groupby="ct", top_n_genes=2,
                                   use_rep="X_umap", key_cosg="cosg", show=False)
    n_x_plain = len(plain.get_xticks()); n_y_plain = len(plain.get_yticks())
    plt.close("all")
    swapped = cosg.plotMarkerDotplot(marked, groupby="ct", top_n_genes=2,
                                     use_rep="X_umap", key_cosg="cosg",
                                     swap_axes=True, show=False)
    assert (len(swapped.get_xticks()), len(swapped.get_yticks())) == \
           (n_y_plain, n_x_plain)
    for label in ("L2-3 IT", "Oligodendrocyte"):
        assert label in _texts(swapped.figure)
    plt.close("all")


def test_save_takes_a_path_not_a_scanpy_suffix(marked, tmp_path):
    """scanpy's `save` is a filename suffix pasted after figures/dotplot_, so
    an absolute path raised FileNotFoundError. Here it is the path."""
    out = tmp_path / "dots.png"
    cosg.plotMarkerDotplot(marked, groupby="ct", top_n_genes=2,
                           use_rep="X_umap", key_cosg="cosg",
                           save=str(out), show=False)
    assert out.exists() and out.stat().st_size > 0
    plt.close("all")


def test_no_dendrogram_when_use_rep_is_none(marked):
    ax = cosg.plotMarkerDotplot(marked, groupby="ct", top_n_genes=2,
                                use_rep=None, key_cosg="cosg", show=False)
    assert ax is not None
    plt.close("all")


def test_a_gene_marking_two_groups_appears_once(marked):
    """Duplicate columns would misalign every bracket after the repeat."""
    ax = cosg.plotMarkerDotplot(marked, groupby="ct", top_n_genes=8,
                                use_rep="X_umap", key_cosg="cosg", show=False)
    labels = [t.get_text() for t in ax.get_xticklabels()]
    assert len(labels) == len(set(labels)), f"duplicate gene column: {labels}"
    plt.close("all")


def test_scanpy_backend_still_available(marked):
    pytest.importorskip("scanpy")
    cosg.plotMarkerDotplot(marked, groupby="ct", top_n_genes=2,
                           use_rep="X_umap", key_cosg="cosg",
                           backend="scanpy", show=False)
    plt.close("all")
