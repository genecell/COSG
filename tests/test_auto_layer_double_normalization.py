"""`layer='auto'` must not log1p an already-normalised matrix.

`from_anndata` stores whatever `adata.X` held in `{modality}_counts` -- the
name says counts, the contents are whatever you had. For a normalised AnnData
that is log-normalised values, and resolving RNA to `log1p` by modality lookup
alone applied log1p a second time. The damage was subtle: correlation ~0.99
against the in-memory reference, errors large enough to read as biology.
"""
import warnings

import numpy as np
import pytest
import scipy.sparse as sp


def _adata(normalised, n=200, g=60, seed=0):
    ad = pytest.importorskip("anndata")
    import pandas as pd
    rs = np.random.RandomState(seed)
    X = sp.csr_matrix(rs.poisson(2.0, (n, g)).astype(np.float32))
    a = ad.AnnData(X=X)
    a.var_names = [f"g{i}" for i in range(g)]
    a.obs["grp"] = pd.Categorical([f"c{i % 3}" for i in range(n)])
    if normalised:
        import scanpy as sc
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
    return a


def _write(a, tmp_path, name):
    cytome = pytest.importorskip("cytome")
    p = str(tmp_path / name)
    cytome.from_anndata(a, output=p).close()
    return p


def test_auto_keeps_log1p_on_real_counts(tmp_path):
    cytome = pytest.importorskip("cytome")
    from cosg._cytome_streaming import _resolve_auto_layer
    p = _write(_adata(normalised=False), tmp_path, "counts.cytome")
    ds = cytome.open(p)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error")        # must NOT warn
            assert _resolve_auto_layer("auto", "RNA", ds=ds) == "log1p"
    finally:
        ds.close()


def test_auto_refuses_to_normalise_twice(tmp_path):
    cytome = pytest.importorskip("cytome")
    pytest.importorskip("scanpy")
    from cosg._cytome_streaming import _resolve_auto_layer
    p = _write(_adata(normalised=True), tmp_path, "norm.cytome")
    ds = cytome.open(p)
    try:
        with pytest.warns(UserWarning):
            resolved = _resolve_auto_layer("auto", "RNA", ds=ds)
        # Never log1p. Which layer it lands on depends on the writer:
        # cytome < 0.3.0 stored the normalised matrix as RNA_counts, so the
        # probe catches it there; >= 0.3.0 does not write RNA_counts at all,
        # and auto falls back to whatever adata.X was.
        assert resolved != "log1p", resolved
        if ds.matrix_meta("RNA_counts") is not None:
            assert resolved == "counts", resolved
        else:
            assert ds.matrix_meta(f"RNA_{resolved}") is not None, resolved
    finally:
        ds.close()


def test_an_explicit_layer_is_never_second_guessed(tmp_path):
    cytome = pytest.importorskip("cytome")
    pytest.importorskip("scanpy")
    from cosg._cytome_streaming import _resolve_auto_layer
    p = _write(_adata(normalised=True), tmp_path, "norm2.cytome")
    ds = cytome.open(p)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            assert _resolve_auto_layer("log1p", "RNA", ds=ds) == "log1p"
            assert _resolve_auto_layer("counts", "RNA", ds=ds) == "counts"
    finally:
        ds.close()


def test_without_a_dataset_behaviour_is_unchanged():
    """Callers that cannot supply ds keep the old lookup rather than guessing."""
    from cosg._cytome_streaming import _resolve_auto_layer
    assert _resolve_auto_layer("auto", "RNA") == "log1p"
    assert _resolve_auto_layer("auto", "ATAC") == "tfidf"


def test_the_probe_reads_the_matrix_out_of_the_chunk_tuple(tmp_path):
    """Regression: `ds.iter_chunks` yields ``(matrix, row_indices)``.

    Reading ``.data`` off the tuple returned None for every chunk, so the
    loop skipped all of them and the probe always answered "cannot tell" --
    which made the double-normalisation guard dead code that still passed a
    smoke test. This asserts the probe reaches a verdict on real input, in
    both directions, which is what the tuple bug prevented.
    """
    import numpy as np
    import pytest
    import scipy.sparse as sp

    anndata = pytest.importorskip("anndata")
    cytome = pytest.importorskip("cytome")
    import pandas as pd
    from cosg._cytome_streaming import _stored_matrix_is_integer

    rs = np.random.RandomState(0)

    def build(name, normalize):
        a = anndata.AnnData(
            X=sp.csr_matrix(rs.poisson(2.0, (200, 30)).astype(np.float32)))
        a.obs["grp"] = pd.Categorical(["a", "b"] * 100)
        a.var_names = [f"g{i}" for i in range(30)]
        if normalize:
            a.X = sp.csr_matrix(np.log1p(
                np.asarray(a.X.todense()) / np.asarray(a.X.sum(1)) * 1e4))
        p = str(tmp_path / name)
        cytome.from_anndata(a, output=p).close()
        return p

    ds = cytome.open(build("raw.cytome", False))
    try:
        assert _stored_matrix_is_integer(ds, "RNA") is True
    finally:
        ds.close()

    ds = cytome.open(build("norm.cytome", True))
    try:
        # cytome >= 0.3.0 does not write RNA_counts for a normalized matrix,
        # so there is nothing to probe -- None is correct there. On older
        # files the probe must still reach a verdict rather than give up,
        # which is how the bug hid.
        verdict = _stored_matrix_is_integer(ds, "RNA")
        has_counts = ds.matrix_meta("RNA_counts") is not None
        assert verdict is (False if has_counts else None), (verdict, has_counts)
    finally:
        ds.close()


def test_recorded_is_integer_is_preferred_over_probing(tmp_path, monkeypatch):
    """cytome >= 0.3.0 records `is_integer` at write time; use it.

    The probe exists for older files. When the writer already recorded the
    answer, re-deriving it is wasted work and one more place to get wrong --
    this probe was dead code once precisely because every consumer had to
    implement it.
    """
    import numpy as np
    import pytest
    import scipy.sparse as sp

    anndata = pytest.importorskip("anndata")
    cytome = pytest.importorskip("cytome")
    from cosg._cytome_streaming import _stored_matrix_is_integer

    rs = np.random.RandomState(0)
    a = anndata.AnnData(X=sp.csr_matrix(rs.poisson(3.0, (40, 10)).astype(np.float32)))
    a.var_names = [f"g{i}" for i in range(10)]
    p = str(tmp_path / "r.cytome")
    cytome.from_anndata(a, output=p).close()

    ds = cytome.open(p)
    try:
        meta = ds.matrix_meta("RNA_counts")
        if meta is None or meta.get("is_integer") is None:
            pytest.skip("cytome older than 0.3.0 does not record is_integer")
        # make the probe path fail loudly: if it is used, the test errors
        monkeypatch.setattr(ds, "iter_chunks",
                            lambda *a, **k: (_ for _ in ()).throw(
                                AssertionError("probed despite a recorded flag")))
        assert _stored_matrix_is_integer(ds, "RNA") is True
    finally:
        ds.close()


def test_is_integer_from_matrix_meta_is_preferred_over_the_probe(tmp_path):
    """cytome >= 0.3.0 records the answer at write time; use it.

    Re-deriving a fact the writer already knew is how the probe in this module
    shipped as dead code -- it read `.data` off a `(matrix, row_indices)`
    tuple, answered "cannot tell" for every chunk, and the guard it feeds
    never fired while looking correct.
    """
    import sqlite3

    import numpy as np
    import pytest
    import scipy.sparse as sp

    anndata = pytest.importorskip("anndata")
    cytome = pytest.importorskip("cytome")
    from cosg._cytome_streaming import _stored_matrix_is_integer

    rs = np.random.RandomState(0)
    a = anndata.AnnData(
        X=sp.csr_matrix(rs.poisson(2.0, (80, 20)).astype(np.float32)))
    a.var_names = [f"g{i}" for i in range(20)]
    p = str(tmp_path / "raw.cytome")
    cytome.from_anndata(a, output=p).close()

    con = sqlite3.connect(p)
    try:
        cols = [r[1] for r in con.execute("PRAGMA table_info(matrix_meta)")]
        if "is_integer" not in cols:
            pytest.skip("cytome predates matrix_meta.is_integer")
        # contradict the data: the recorded flag must win, which is how we
        # know it is being read at all
        con.execute("UPDATE matrix_meta SET is_integer = 0 "
                    "WHERE matrix_name = 'RNA_counts'")
        con.commit()
    finally:
        con.close()

    ds = cytome.open(p)
    try:
        assert _stored_matrix_is_integer(ds, "RNA") is False
    finally:
        ds.close()
