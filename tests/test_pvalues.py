"""The p-value layer, tested where it is actually load-bearing: the tail.

A marker screen acts on small p-values. BH at FDR 0.01 over 20,000 genes puts
the decision boundary near 1e-5, so a method that is correct in the bulk and
100x wrong at 1e-5 is wrong exactly where it is used -- and a KS test on the
p-value distribution will not notice, because KS measures the bulk. The tail
tests below are the ones that would have failed the original spec.
"""
from __future__ import annotations

from itertools import combinations

import numpy as np
import pytest
from scipy.stats import norm, kstest

from cosg._pvalues import (
    population_moments, group_sums, analytic_pvalues, permutation_pvalues,
    adjust_fdr, compute_pvalues, _spa_tail, _exact_tail,
)

sparse = pytest.importorskip("scipy.sparse")
anndata = pytest.importorskip("anndata")
pd = pytest.importorskip("pandas")
import cosg


# ---------------------------------------------------------------------------
# 1. The closed-form moments are exact, not approximate
# ---------------------------------------------------------------------------

def test_moments_match_complete_enumeration():
    """Enumerate every C(12, 4) assignment and compare all three moments."""
    rng = np.random.default_rng(0)
    n, n_c = 12, 4
    x = np.round(rng.gamma(0.3, 4.0, n), 3)
    x[rng.random(n) < 0.5] = 0.0          # zeros, as in real expression

    # want_m3: the third moment is opt-in now that the saddlepoint
    # replaced the Edgeworth term, but the formula is still part of the
    # streaming API's contract, so it stays under test.
    mean, s2, m3 = population_moments(x.reshape(-1, 1), want_m3=True)
    mean, s2, m3 = float(mean[0]), float(s2[0]), float(m3[0])

    T = np.array([x[list(c)].sum() for c in combinations(range(n), n_c)])
    E_cf = n_c * mean
    V_cf = n_c * (n - n_c) / (n - 1) * s2
    M3_cf = n_c * (n - n_c) * (n - 2 * n_c) / ((n - 1) * (n - 2)) * m3

    assert np.isclose(T.mean(), E_cf)
    assert np.isclose(T.var(), V_cf)
    assert np.isclose(((T - T.mean()) ** 3).mean(), M3_cf)


def test_exact_tail_matches_enumeration():
    """The enumerated tail is the truth for a gene in few cells."""
    rng = np.random.default_rng(3)
    n, k, n_c = 16, 6, 5          # C(16,5) = 4368 to enumerate
    x = np.zeros(n)
    x[:k] = rng.gamma(2.0, 1.0, k)
    t = np.quantile([x[list(c)].sum()
                     for c in combinations(range(n), n_c)], 0.97)
    brute = np.mean([x[list(c)].sum() >= t - 1e-12
                     for c in combinations(range(n), n_c)])
    assert np.isclose(_exact_tail(x[:k], n, n_c, t), brute, rtol=1e-9)


# ---------------------------------------------------------------------------
# 2. The tail: the reason the saddlepoint is here at all
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("k,expected_in_group", [(2000, 20.0), (200, 2.0)])
def test_saddlepoint_beats_normal_in_the_tail(k, expected_in_group):
    """At nominal 1e-5 the normal is off by 10-100x; the saddlepoint is not.

    This is the measurement that drove the design. The threshold is set from
    the normal approximation, the truth is sampled, and the saddlepoint is
    asked for the same tail.
    """
    rng = np.random.default_rng(1)
    n, n_c, B = 20_000, 200, 60_000
    x = np.zeros(n)
    x[:k] = rng.gamma(2.0, 1.0, k)
    assert np.isclose(n_c * k / n, expected_in_group)

    mean, s2, _ = population_moments(x.reshape(-1, 1))
    E = n_c * float(mean[0])
    sd = np.sqrt(n_c * (n - n_c) / (n - 1) * float(s2[0]))

    t = E + norm.isf(1e-4) * sd
    # Exact but O(k): only expressing cells contribute, so draw how many land
    # in the group (hypergeometric) and sum that many of their values without
    # replacement. Sampling whole permutations of 20,000 cells is the same
    # distribution and a hundred times slower.
    vals = x[:k]
    hits = rng.hypergeometric(k, n - k, n_c, size=B)
    exceed = 0
    for h in np.unique(hits):
        m = int((hits == h).sum())
        if h == 0:
            exceed += m * (0.0 >= t)
            continue
        # take the h smallest random keys per draw: a without-replacement
        # subset of size h, vectorised, in blocks so the key matrix stays small
        for lo in range(0, m, 2000):
            rows = min(2000, m - lo)
            keys = rng.random((rows, k))
            idx = np.argpartition(keys, h - 1, axis=1)[:, :h]
            exceed += int((vals[idx].sum(axis=1) >= t).sum())
    truth = exceed / B
    spa = _spa_tail(vals, n, n_c, t)

    assert spa is not None
    # the normal claims 1e-4 here; require the saddlepoint within 2x of truth
    assert truth / spa < 2.0 and spa / truth < 2.0, (
        f"k={k}: truth {truth:.2e}, saddlepoint {spa:.2e}, normal 1e-4")
    # and confirm the normal really is the worse one, so this test keeps its point
    assert truth / 1e-4 > 3.0, f"fixture no longer exercises the tail: {truth:.2e}"


# ---------------------------------------------------------------------------
# 3. Calibration and power
# ---------------------------------------------------------------------------

def _fixture(n=600, n_genes=200, n_groups=3, seed=0, planted=0):
    rng = np.random.default_rng(seed)
    lab = np.array([f"c{i % n_groups}" for i in range(n)])
    X = rng.poisson(0.5, (n, n_genes)).astype(np.float64)
    for j in range(planted):
        m = lab == f"c{j % n_groups}"
        X[m, j] += rng.poisson(10, m.sum())
    return X, lab, [f"c{i}" for i in range(n_groups)]


def test_null_pvalues_are_uniform():
    X, lab, order = _fixture(seed=5)
    P, _, _ = analytic_pvalues(X, lab, order, method="spa")
    assert kstest(P.ravel(), "uniform").pvalue > 0.01


def test_null_tail_is_not_anticonservative():
    """The check KS cannot make: how many p-values fall below a small cut."""
    X, lab, order = _fixture(n=900, n_genes=400, seed=7)
    P, _, _ = analytic_pvalues(X, lab, order, method="spa")
    for cut in (1e-2, 1e-3):
        rate = (P <= cut).mean()
        assert rate < 5 * cut, f"nominal {cut}, observed {rate:.2e}"


def test_planted_markers_clear_fdr():
    X, lab, order = _fixture(n=600, n_genes=200, planted=3, seed=1)
    out = compute_pvalues(X, lab, order, method="spa")
    fdr = out["p_value_fdr"]
    for j in range(3):
        assert fdr[j, j] < 1e-3, f"planted marker {j} not significant"


def test_analytic_agrees_with_sampled_permutation():
    """The oracle: sampled permutation p, on rows where sampling can resolve."""
    X, lab, order = _fixture(n=300, n_genes=60, planted=2, seed=2)
    P, _, _ = analytic_pvalues(X, lab, order, method="spa")
    Pp = permutation_pvalues(X, lab, order, n_permutations=4000, random_seed=0)
    mid = (Pp > 5e-3) & (Pp < 0.95)          # where 4000 draws have resolution
    assert mid.sum() > 30
    assert np.corrcoef(np.log10(P[mid]), np.log10(Pp[mid]))[0, 1] > 0.9


# ---------------------------------------------------------------------------
# 4. Batch: the null must be stratified or it is wrong
# ---------------------------------------------------------------------------

def test_batch_confounding_needs_stratification():
    """With batch structure the unstratified null is anti-conservative."""
    rng = np.random.default_rng(11)
    n, n_genes = 900, 150
    batch = np.array(["b0"] * (n // 2) + ["b1"] * (n - n // 2))
    # group is correlated with batch, and genes shift with batch: no real
    # group effect, only the confound
    lab = np.where(rng.random(n) < np.where(batch == "b0", 0.8, 0.2), "c0", "c1")
    X = rng.poisson(np.where(batch == "b0", 2.0, 0.5)[:, None], (n, n_genes)).astype(float)
    order = ["c0", "c1"]

    P_un, _, _ = analytic_pvalues(X, lab, order, method="spa")
    P_st, _, _ = analytic_pvalues(X, lab, order, batch=batch, method="spa")
    assert (P_un <= 1e-3).mean() > 10 * (P_st <= 1e-3).mean(), (
        "stratification should remove most of the confounded signal")


# ---------------------------------------------------------------------------
# 5. Contract: defaults, determinism, mu-invariance
# ---------------------------------------------------------------------------

def _adata(seed=0):
    X, lab, _ = _fixture(n=450, n_genes=120, planted=3, seed=seed)
    return anndata.AnnData(
        X=sparse.csr_matrix(X.astype(np.float32)),
        obs=pd.DataFrame({"ct": pd.Categorical(lab)},
                         index=[f"cell{i}" for i in range(X.shape[0])]),
        var=pd.DataFrame(index=[f"g{i}" for i in range(X.shape[1])]))


def test_default_is_unchanged():
    """calculate_pvalues=False must not add keys or alter the scores."""
    a, b = _adata(), _adata()
    cosg.cosg(a, groupby="ct", key_added="cosg", n_genes_user=5)
    cosg.cosg(b, groupby="ct", key_added="cosg", n_genes_user=5,
              calculate_pvalues=False)
    assert set(a.uns["cosg"]) == set(b.uns["cosg"])
    assert "pvals" not in a.uns["cosg"]
    np.testing.assert_array_equal(
        pd.DataFrame(a.uns["cosg"]["names"]).values,
        pd.DataFrame(b.uns["cosg"]["names"]).values)


def test_pvalues_do_not_depend_on_mu():
    """mu ranks; chance does not care.

    Compared per *gene*, not per row. The table lists the top-ranked genes and
    mu changes that ranking, so row i is a different gene at a different mu --
    the invariance is of the p-value attached to a gene, not of the table's
    layout. Worth pinning both halves: the same gene keeps its p-value, and
    the ranking is free to move.
    """
    per_gene = {}
    tables = {}
    for mu in (1.0, 100.0):
        a = _adata()
        cosg.cosg(a, groupby="ct", key_added="cosg", n_genes_user=8, mu=mu,
                  calculate_pvalues=True)
        names = pd.DataFrame(a.uns["cosg"]["names"])
        pvals = pd.DataFrame(a.uns["cosg"]["pvals"])
        tables[mu] = names
        per_gene[mu] = {(c, names[c][i]): float(pvals[c][i])
                        for c in names.columns for i in range(len(names))}

    shared = set(per_gene[1.0]) & set(per_gene[100.0])
    assert len(shared) >= 6, f"only {len(shared)} genes in common to compare"
    for key in shared:
        assert per_gene[1.0][key] == pytest.approx(per_gene[100.0][key],
                                                   rel=1e-12, abs=0), (
            f"p-value for {key} moved with mu")


def test_analytic_path_is_deterministic_and_leaves_global_rng_alone():
    X, lab, order = _fixture(planted=2, seed=4)
    before = np.random.get_state()[1][:8].copy()
    p1, _, _ = analytic_pvalues(X, lab, order, method="spa")
    p2, _, _ = analytic_pvalues(X, lab, order, method="spa")
    after = np.random.get_state()[1][:8]
    np.testing.assert_array_equal(p1, p2)
    np.testing.assert_array_equal(before, after)


def test_bad_method_is_rejected():
    X, lab, order = _fixture()
    with pytest.raises(ValueError, match="pvalue_method"):
        compute_pvalues(X, lab, order, method="magic")


def test_fdr_matches_a_known_case():
    p = np.array([0.001, 0.008, 0.039, 0.041, 0.042])
    got = adjust_fdr(p, method="fdr_bh")
    assert np.all(np.diff(got) >= -1e-12)          # monotone
    assert np.isclose(got[0], 0.005)               # 0.001 * 5 / 1
    assert np.all(got <= 1.0)


# ---------------------------------------------------------------------------
# 6. The cytome path must agree, and must actually receive the request
# ---------------------------------------------------------------------------

def test_cytome_path_forwards_the_pvalue_request(tmp_path):
    """`calculate_pvalues=True` was accepted and dropped on the way to cytome.

    The caller asked for p-values, got a result with no p-value keys, and no
    complaint -- the kwarg was simply absent from the forwarding list. Nothing
    downstream could notice, because a missing key looks like a caller that
    never asked.
    """
    cytome = pytest.importorskip("cytome")
    X, lab, _ = _fixture(n=450, n_genes=100, planted=3, seed=0)
    a = anndata.AnnData(
        X=sparse.csr_matrix(X.astype(np.float32)),
        obs=pd.DataFrame({"ct": pd.Categorical(lab)},
                         index=[f"cell{i}" for i in range(X.shape[0])]),
        var=pd.DataFrame(index=[f"g{i}" for i in range(X.shape[1])]))
    path = str(tmp_path / "t.cytome")
    cytome.from_anndata(a, output=path, force=True).close()

    out = cosg.cosg(path, groupby="ct", n_genes_user=5,
                    calculate_pvalues=True, layer="counts")
    assert any("pval" in k for k in out), f"no p-value keys returned: {list(out)}"


def test_cytome_and_in_memory_agree_on_the_same_layer(tmp_path):
    """Streamed and in-memory must give the same number for the same gene.

    On the same layer they agree exactly. They do *not* agree when the layer
    differs -- the cytome default is `layer='auto'`, which is normalised, so a
    parity check that omits `layer=` is comparing two different matrices and
    will look like a p-value bug.
    """
    cytome = pytest.importorskip("cytome")
    X, lab, _ = _fixture(n=450, n_genes=100, planted=3, seed=0)
    a = anndata.AnnData(
        X=sparse.csr_matrix(X.astype(np.float32)),
        obs=pd.DataFrame({"ct": pd.Categorical(lab)},
                         index=[f"cell{i}" for i in range(X.shape[0])]),
        var=pd.DataFrame(index=[f"g{i}" for i in range(X.shape[1])]))
    path = str(tmp_path / "t.cytome")
    cytome.from_anndata(a, output=path, force=True).close()

    cosg.cosg(a, groupby="ct", key_added="cosg", n_genes_user=5,
              calculate_pvalues=True)
    names = pd.DataFrame(a.uns["cosg"]["names"])
    pvals = pd.DataFrame(a.uns["cosg"]["pvals"])
    mem = {(c, names[c][i]): float(pvals[c][i])
           for c in names.columns for i in range(len(names))}

    out = cosg.cosg(path, groupby="ct", n_genes_user=5,
                    calculate_pvalues=True, layer="counts")
    nm = np.asarray(out["names"])
    pv = np.asarray(out["pvals"], dtype=float)
    order = list(out["groups_order"])

    compared = 0
    for j, c in enumerate(order):
        for i in range(nm.shape[0]):
            key = (c, nm[i, j])
            if key in mem:
                # 1e-9, not exact: both paths run the same Newton to the
                # same root, but they assemble their batches independently, so
                # the last couple of digits are convergence noise rather than
                # disagreement. A tolerance tighter than the solver's is a test
                # that fails for being precise about nothing.
                assert pv[i, j] == pytest.approx(mem[key], rel=1e-9, abs=0)
                compared += 1
    assert compared >= 9, f"only {compared} gene x group pairs compared"


# ---------------------------------------------------------------------------
# 6. The streamed path must agree with the in-memory one
# ---------------------------------------------------------------------------

def test_cytome_streaming_matches_in_memory():
    """Same data, same p-values, whether or not the matrix fits in memory.

    The streaming path never sees the matrix: it accumulates the first three
    power sums per gene during the pass it was making anyway, then refines the
    flagged tail entries in a second pass over those genes only.
    """
    cytome = pytest.importorskip("cytome")
    import tempfile, os

    a = _adata(seed=3)
    tmp = tempfile.mkdtemp()
    path = os.path.join(tmp, "t.cytome")
    cytome.from_anndata(a, output=path, force=True).close()

    cosg.cosg(a, groupby="ct", key_added="cosg", n_genes_user=5,
              calculate_pvalues=True, pvalue_method="spa")
    mem_names = pd.DataFrame(a.uns["cosg"]["names"])
    mem_p = pd.DataFrame(a.uns["cosg"]["pvals"]).astype(float)

    from cosg._cytome_streaming import run_cosg_cytome
    out = run_cosg_cytome(path, groupby="ct", modality="RNA", layer="counts",
                          n_genes_user=5, calculate_pvalues=True,
                          pvalue_method="spa")
    assert "pvals" in out and "pvals_adj" in out

    order = list(out["groups_order"])
    st_names = np.asarray(out["names"])
    st_p = np.asarray(out["pvals"], dtype=float)

    compared = 0
    for j, grp in enumerate(order):
        if grp not in mem_names.columns:
            continue
        lookup = dict(zip(mem_names[grp], mem_p[grp]))
        for i, gene in enumerate(st_names[:, j]):
            if gene in lookup:
                assert np.isclose(np.log10(st_p[i, j]), np.log10(lookup[gene]),
                                  rtol=1e-6, atol=1e-6), (
                    f"{grp}/{gene}: streamed {st_p[i, j]:.3e} vs "
                    f"in-memory {lookup[gene]:.3e}")
                compared += 1
    assert compared >= 6, f"only {compared} genes compared"


def test_streaming_rejects_batched_pvalues():
    """A stratified null needs per-batch power sums; say so rather than
    silently computing an unstratified p-value that would be wrong."""
    cytome = pytest.importorskip("cytome")
    import tempfile, os

    a = _adata(seed=6)
    a.obs["batch"] = pd.Categorical(
        ["b0" if i % 2 else "b1" for i in range(a.n_obs)])
    tmp = tempfile.mkdtemp()
    path = os.path.join(tmp, "t.cytome")
    cytome.from_anndata(a, output=path, force=True).close()

    from cosg._cytome_streaming import run_cosg_cytome
    with pytest.raises(NotImplementedError, match="batch_key"):
        run_cosg_cytome(path, groupby="ct", modality="RNA", layer="counts",
                        n_genes_user=5, calculate_pvalues=True,
                        batch_key="batch")


# ---------------------------------------------------------------------------
# 6. The streaming path must agree with the in-memory one
# ---------------------------------------------------------------------------

def test_cytome_streaming_matches_in_memory():
    """Same p-values from power sums as from the whole matrix.

    The streaming path accumulates the first three power sums per gene while
    it is already streaming for the cosine, then refines flagged genes in a
    second pass. It has no access to the matrix, so this is the check that the
    two routes compute the same quantity.
    """
    cytome = pytest.importorskip("cytome")
    import tempfile, os
    from cosg._cytome_streaming import run_cosg_cytome

    a = _adata(seed=0)
    d = tempfile.mkdtemp()
    path = os.path.join(d, "t.cytome")
    cytome.from_anndata(a, output=path, force=True).close()

    cosg.cosg(a, groupby="ct", key_added="cosg", n_genes_user=5,
              calculate_pvalues=True, verbosity=0)
    mem_p = pd.DataFrame(a.uns["cosg"]["pvals"]).astype(float)
    mem_n = pd.DataFrame(a.uns["cosg"]["names"])

    out = run_cosg_cytome(path, groupby="ct", n_genes_user=5, layer="counts",
                          calculate_pvalues=True, output_format="dict",
                          verbose=False)
    assert "pvals" in out, "streaming computed p-values and dropped them"
    pv = out["pvals"]

    checked = 0
    for c in mem_n.columns:
        top = str(mem_n[c][0])
        got = pv.get((str(c), top))
        assert got is not None, f"no streaming p-value for {c}/{top}"
        assert np.isclose(float(mem_p[c][0]), got, rtol=1e-9), (
            f"{c}/{top}: in-memory {mem_p[c][0]:.3e} vs streaming {got:.3e}")
        checked += 1
    assert checked == len(mem_n.columns)


def test_streaming_rejects_batched_pvalues_rather_than_faking_them():
    """A stratified null needs per-batch power sums; say so, do not approximate."""
    cytome = pytest.importorskip("cytome")
    import tempfile, os
    from cosg._cytome_streaming import run_cosg_cytome

    a = _adata(seed=1)
    a.obs["batch"] = pd.Categorical(
        ["b0"] * (a.n_obs // 2) + ["b1"] * (a.n_obs - a.n_obs // 2))
    d = tempfile.mkdtemp()
    path = os.path.join(d, "b.cytome")
    cytome.from_anndata(a, output=path, force=True).close()

    with pytest.raises(NotImplementedError, match="batch"):
        run_cosg_cytome(path, groupby="ct", n_genes_user=5, layer="counts",
                        batch_key="batch", calculate_pvalues=True,
                        output_format="dict", verbose=False)


# ---------------------------------------------------------------------------
# 7. Stratified tails, and the vectorised solver that makes them affordable
# ---------------------------------------------------------------------------

def test_vectorised_matches_the_scalar_solver():
    """One row at a time and many rows at once must give the same number."""
    from cosg._pvalues import _spa_tail, _spa_tail_vectorised

    rng = np.random.default_rng(0)
    n = 4000
    rows = []
    for k, nc in [(300, 400), (150, 200), (60, 400), (1200, 300), (40, 150)]:
        x = np.zeros(n)
        x[:k] = rng.gamma(2.0, 1.0, k)
        mean, s2, _ = population_moments(x.reshape(-1, 1))
        E = nc * float(mean[0])
        sd = np.sqrt(nc * (n - nc) / (n - 1) * float(s2[0]))
        rows.append((x[:k], nc, E + 3.5 * sd))

    scalar = [_spa_tail(v, n, nc, t) for v, nc, t in rows]
    vals = np.concatenate([r[0] for r in rows])
    row_ptr = np.cumsum([0] + [r[0].size for r in rows])
    p, ok = _spa_tail_vectorised(
        vals, row_ptr, row_ptr.copy(), np.arange(len(rows)),
        np.array([float(n - r[0].size) for r in rows]),
        np.array([float(r[1]) for r in rows]),
        np.array([r[2] for r in rows]))
    assert ok.all()
    for got, want in zip(p, scalar):
        assert got == pytest.approx(want, rel=1e-9)


def test_stratified_tail_matches_a_within_batch_permutation_null():
    """The reason the stratified tail had to exist before it could be used.

    Two batches of different depth, with the group taking 20% of one and 4% of
    the other. Permuting labels within batch is the correct null; permuting
    across batches is a different, much narrower one -- and using the
    unstratified tail on stratified data does not lose a little accuracy, it
    reports 1e-9 for something that happens once in a thousand.
    """
    from cosg._pvalues import _spa_tail, _spa_tail_vectorised

    rng = np.random.default_rng(3)
    nb, ncb = [1500, 1500], [300, 60]
    xs, st = [], []
    for b in (0, 1):
        xi = np.zeros(nb[b])
        xi[:250] = rng.gamma(2.0, 3.0 if b == 0 else 0.4, 250)
        xs.append(xi)
        st.append(np.full(nb[b], b))
    x = np.concatenate(xs)
    strat = np.concatenate(st)
    n = x.size

    idx_b = [np.flatnonzero(strat == b) for b in (0, 1)]
    T = np.array([sum(x[rng.choice(idx_b[b], ncb[b], replace=False)].sum()
                      for b in (0, 1)) for _ in range(60_000)])
    t = np.quantile(T, 1 - 1e-3)
    truth = (T >= t).mean()

    nz = x != 0
    order = np.argsort(strat[nz], kind="stable")
    vals, sb = x[nz][order], strat[nz][order]
    p, ok = _spa_tail_vectorised(
        vals,
        np.array([0, vals.size]),
        np.array([0, int((sb == 0).sum()), vals.size]),
        np.array([0, 0]),
        np.array([float((strat == 0).sum() - (sb == 0).sum()),
                  float((strat == 1).sum() - (sb == 1).sum())]),
        np.array([float(ncb[0]), float(ncb[1])]),
        np.array([t]))

    assert ok[0]
    assert 0.5 < p[0] / truth < 2.0, f"stratified {p[0]:.2e} vs truth {truth:.2e}"

    unstrat = _spa_tail(x[nz], n, sum(ncb), t)
    assert unstrat < truth / 100, (
        "the unstratified tail should be badly wrong here -- if it is not, "
        "this fixture no longer demonstrates why stratification matters")


def test_batched_run_refines_its_tail():
    """With batch_key the tail is refined, not skipped."""
    rng = np.random.default_rng(5)
    n, n_genes = 900, 120
    batch = np.array(["b0"] * (n // 2) + ["b1"] * (n - n // 2))
    lab = np.array([f"c{i % 3}" for i in range(n)])
    X = rng.poisson(np.where(batch == "b0", 2.0, 0.6)[:, None],
                    (n, n_genes)).astype(float)
    X[lab == "c0", 0] += rng.poisson(8, (lab == "c0").sum())
    order = ["c0", "c1", "c2"]

    _, _, info = analytic_pvalues(X, lab, order, batch=batch, method="spa")
    assert info["n_spa"] > 0, (
        f"batched run refined nothing: {info} -- the stratified tail is the "
        "whole point of this path")


# ---------------------------------------------------------------------------
# 7. The stratified tail: refined now, not skipped
# ---------------------------------------------------------------------------

def test_stratified_saddlepoint_matches_a_sampled_stratified_null():
    """With strata the tail is refined against the right null, not the pooled one.

    The unknowns become (s, u_1..u_B) -- one multiplier per stratum -- and the
    Jacobian is an arrow, dense in s and diagonal in the u_b, so a Schur
    complement solves it in O(B) per row. Checked against a sampled stratified
    null: draw how many expressing cells land in the group *within each
    stratum*, independently, and sum.
    """
    from cosg._pvalues import _spa_tail_vectorised

    rng = np.random.default_rng(4)
    nb, kb, ncb = [1200, 800], [120, 40], [300, 200]
    vals = [rng.gamma(2.0, 1.0, kb[0]), rng.gamma(2.0, 3.0, kb[1])]
    V = np.concatenate(vals)

    B = 60_000
    Ts = np.empty(B)
    for i in range(B):
        tot = 0.0
        for b in range(2):
            h = rng.hypergeometric(kb[b], nb[b] - kb[b], ncb[b])
            if h:
                tot += rng.choice(vals[b], h, replace=False).sum()
        Ts[i] = tot

    for q in (0.99, 0.999):
        t = np.quantile(Ts, q)
        empirical = (Ts >= t).mean()
        p, ok = _spa_tail_vectorised(
            V, np.array([0, V.size]), np.array([0, kb[0], V.size]),
            np.array([0, 0]),
            np.array([float(nb[0] - kb[0]), float(nb[1] - kb[1])]),
            np.array([float(ncb[0]), float(ncb[1])]),
            np.array([t]))
        assert ok[0], f"stratified solve did not converge at q={q}"
        assert 0.5 < empirical / p[0] < 2.0, (
            f"q={q}: empirical {empirical:.3e}, stratified SPA {p[0]:.3e}")


def test_batched_run_refines_its_tail(): 
    """A batched call must now report refined rows, not zero."""
    rng = np.random.default_rng(12)
    n, n_genes = 900, 200
    batch = np.array(["b0"] * (n // 2) + ["b1"] * (n - n // 2))
    lab = np.array([f"c{i % 3}" for i in range(n)])
    X = rng.poisson(np.where(batch == "b0", 2.0, 0.6)[:, None],
                    (n, n_genes)).astype(float)
    for j in range(3):
        m = lab == f"c{j}"
        X[m, j] += rng.poisson(12, m.sum())
    order = [f"c{i}" for i in range(3)]

    _, _, info = analytic_pvalues(X, lab, order, batch=batch, method="spa")
    assert info["n_spa"] > 0, (
        f"batched run refined nothing -- the stratified tail is not wired: {info}")


def test_vectorised_and_scalar_agree_where_both_apply():
    """One tail, two implementations, and they must not drift."""
    from cosg._pvalues import _spa_tail, _spa_tail_vectorised

    rng = np.random.default_rng(0)
    for (n, k, n_c) in [(2000, 200, 300), (5000, 400, 500), (3000, 900, 600)]:
        x = np.zeros(n)
        x[:k] = rng.gamma(2.0, 1.0, k)
        mean, s2, _ = population_moments(x.reshape(-1, 1))
        E = n_c * float(mean[0])
        sd = np.sqrt(n_c * (n - n_c) / (n - 1) * float(s2[0]))
        for zq in (3.0, 5.0, 8.0):
            t = E + zq * sd
            scalar = _spa_tail(x[:k], n, n_c, t)
            vec, ok = _spa_tail_vectorised(
                x[:k], np.array([0, k]), np.array([0, k]), np.array([0]),
                np.array([float(n - k)]), np.array([float(n_c)]), np.array([t]))
            if scalar is None or not ok[0]:
                continue
            assert vec[0] == pytest.approx(scalar, rel=1e-8), (
                f"n={n} k={k} n_c={n_c} z={zq}: {vec[0]:.4e} vs {scalar:.4e}")


# ---------------------------------------------------------------------------
# 8. The GPU path
# ---------------------------------------------------------------------------

def test_gpu_path_produces_the_same_pvalues():
    """device='gpu' changes where the score is computed, not the p-value.

    The p-value block consumes `cellxgene`, which is only ever sliced -- the
    GPU helpers do their own transfer -- so the tail arithmetic runs on the
    host either way and must agree exactly. Nothing device-specific needs
    writing; this test exists to keep that true.
    """
    pytest.importorskip("cupy")
    a_cpu, a_gpu = _adata(), _adata()
    cosg.cosg(a_cpu, groupby="ct", key_added="cosg", n_genes_user=5,
              calculate_pvalues=True, device="cpu")
    cosg.cosg(a_gpu, groupby="ct", key_added="cosg", n_genes_user=5,
              calculate_pvalues=True, device="gpu")

    def per_gene(a):
        nm = pd.DataFrame(a.uns["cosg"]["names"])
        pv = pd.DataFrame(a.uns["cosg"]["pvals"])
        return {(c, nm[c][i]): float(pv[c][i])
                for c in nm.columns for i in range(len(nm))}

    cpu, gpu = per_gene(a_cpu), per_gene(a_gpu)
    shared = set(cpu) & set(gpu)
    assert len(shared) >= 6, f"only {len(shared)} genes in common"
    for key in shared:
        assert gpu[key] == pytest.approx(cpu[key], rel=1e-9)


def test_pvalue_block_does_not_depend_on_the_device_argument():
    """Guard the above without a GPU: the matrix handed to the tail must be a
    host array, so no device branch can quietly change the arithmetic."""
    import inspect

    src = inspect.getsource(cosg.cosg)
    i = src.find("calculate_pvalues:")
    j = src.find("analytic_pvalues")
    assert j > 0, "the p-value block moved; re-check the device independence"
    # cellxgene is assigned from adata only, never from a device array
    assigns = [l.strip() for l in src.splitlines() if l.strip().startswith("cellxgene =")]
    assert assigns, "cellxgene assignment not found"
    for line in assigns:
        assert "adata" in line, (
            f"cellxgene is assigned from something other than adata ({line!r}); "
            "if it can become a device array the p-value path needs a transfer")


def test_pvalues_are_device_independent():
    """The GPU path gets p-values too, and must keep getting them.

    The p-value block sits after the device branch and reads the host matrix,
    so `device='gpu'` already returns p-values -- computed on the CPU. That is
    worth pinning rather than assuming: a future `if _device == 'gpu': skip`
    added for a good local reason would silently drop the columns on exactly
    the runs where the user is most likely to be at scale.

    Structural, because this machine has GPUs but no cupy, and a test that
    skips is a test that never fails.
    """
    import inspect
    import re

    src = inspect.getsource(cosg.cosg)
    start = src.index("if calculate_pvalues:")
    end = src.index("analytic_pvalues", start)
    head = src[start:end]
    assert not re.search(r"_device\s*==", head), (
        "a device guard appeared in the p-value branch; GPU runs would lose "
        "their p-values silently")


def test_pvalues_are_not_gated_on_the_device():
    """`device='gpu'` must still produce p-values.

    What makes that true is that the GPU kernels take the CPU matrix as an
    argument rather than replacing it, so the p-value block reads the same
    `cellxgene` on either path and needs no device-specific code. This test
    pins that property: the analytic call is fed `cellxgene`, and the block is
    not nested inside a device branch. A future refactor that moves the matrix
    onto the device has to come past this rather than silently drop p-values
    for GPU users.

    Structural rather than behavioural because cupy is not installed here --
    which is the honest reason, and the reason the check exists at all.
    """
    import pathlib as _p
    import cosg

    src = (_p.Path(cosg.__file__).parent / "cosg.py").read_text()
    call = src.index("analytic_pvalues(")
    window = src[call:call + 260]
    assert "cellxgene" in window, (
        "the analytic p-value call no longer reads the shared CPU matrix; on "
        f"the GPU path it may be reading a device array:\n{window[:200]}")

    # the enclosing block must be the calculate_pvalues gate, not a device one
    head = src[:call]
    gate = head.rindex("if calculate_pvalues:")
    assert "_device" not in head[gate:], (
        "a device test appears between `if calculate_pvalues:` and the call, "
        "so p-values may not run for device='gpu'")
