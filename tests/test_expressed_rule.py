"""One meaning of "expressed" on every path: a nonzero value.

The streaming, GPU and cytome paths counted stored entries, explicit zeros
included; dense streaming counted values above 0; the legacy path counted
values different from 0. A gene with 21 true nonzeros (81 stored entries) in
a group of 89 cells passed ``expressed_pct=0.5`` on the default path and was
filtered on ``cpu_chunk_size=0``; a dense matrix with negative values had
fewer "expressing" cells on the default path than on the legacy one.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import cosg
from conftest import requires_cytome
from test_per_group_floor import numpy_gpu  # noqa: F401  (fixture)

N, F = 200, 20


def _matrix(seed=3):
    """Gene g5: 81 cells of A store an entry, 60 of them an explicit zero."""
    rng = np.random.default_rng(seed)
    lab = np.where(np.arange(N) < 89, "A", "B")
    X = sp.random(N, F, density=0.15, format="lil", random_state=rng, dtype=np.float32)
    for i in range(81):
        X[i, 5] = 1.0
    X = X.tocsr()
    X.data = np.round(X.data * 4 + 1).astype(np.float32)
    for i in range(21, 81):
        s, e = X.indptr[i], X.indptr[i + 1]
        X.data[s:e][X.indices[s:e] == 5] = 0.0
    assert X.nnz - sp.csr_matrix(X.toarray()).nnz == 60     # 60 stored zeros
    return X, lab


def _adata(X, lab):
    import anndata as ad
    a = ad.AnnData(X, obs=pd.DataFrame({"grp": pd.Categorical(lab), "batch": np.where(np.arange(N) % 2, "b1", "b2")},
                                       index=[f"c{i}" for i in range(N)]))
    a.var_names = [f"g{j}" for j in range(F)]
    return a


def _removed(a, g="A"):
    c = a.uns["cosg"]
    return set(np.asarray(c["names"][g])[np.asarray(c["scores"][g]) == -1])


PATHS = {"stream": {}, "stream_chunks": dict(cpu_chunk_size=50), "legacy": dict(cpu_chunk_size=0),
         "batch": dict(batch_key="batch"), "batch_legacy": dict(batch_key="batch", cpu_chunk_size=0)}


@pytest.mark.parametrize("path", list(PATHS))
def test_a_stored_zero_is_not_an_expressing_cell(path):
    X, lab = _matrix()
    a = _adata(X, lab)
    cosg.cosg(a, groupby="grp", n_genes_user=F, expressed_pct=0.5, **PATHS[path])
    assert "g5" in _removed(a)                    # 21 of 89 cells: under the floor
    ref = _adata(sp.csr_matrix(X.toarray()), lab)  # the same values without stored zeros
    cosg.cosg(ref, groupby="grp", n_genes_user=F, expressed_pct=0.5, **PATHS[path])
    assert _removed(a) == _removed(ref)


@pytest.mark.parametrize("path", [dict(device="gpu"), dict(device="gpu", gpu_chunk_size=50),
                                  dict(device="gpu", batch_key="batch")])
def test_a_stored_zero_is_not_an_expressing_cell_on_the_gpu_code(numpy_gpu, path):  # noqa: F811
    X, lab = _matrix()
    a = _adata(X, lab)
    cosg.cosg(a, groupby="grp", n_genes_user=F, expressed_pct=0.5, **path)
    assert "g5" in _removed(a)


@pytest.mark.parametrize("ccs", [None, 0, 50])
def test_negative_values_count_as_expressed_on_dense_input(ccs):
    rng = np.random.default_rng(0)
    lab = np.where(np.arange(N) < 100, "A", "B")
    X = np.zeros((N, F), dtype=np.float32)
    X[:100, :10] = rng.normal(0, 1, (100, 10))     # nonzero everywhere in A, half negative
    X[100:, 10:] = rng.normal(0, 1, (100, 10))
    a = _adata(X, lab)
    cosg.cosg(a, groupby="grp", n_genes_user=F, expressed_pct=0.6, cpu_chunk_size=ccs)
    # every one of A's genes is nonzero in all of A's cells: none removed
    assert not ({f"g{j}" for j in range(10)} & _removed(a))


@requires_cytome
@pytest.mark.parametrize("fb", ["none", "chromosome", "column_major"])
def test_a_stored_zero_is_not_an_expressing_cell_on_a_store(tmp_path, fb):
    import cytome
    X, lab = _matrix()
    p = str(tmp_path / "z.cytome")
    ds = cytome.create(p)
    ds.set_entity("cells", pd.DataFrame({"cell_idx": np.arange(N), "barcode": [f"c{i}" for i in range(N)],
                                         "grp": lab}))
    starts = np.arange(F) * 1000
    ds.set_entity("peaks", pd.DataFrame({"peak_idx": np.arange(F), "chr": ["chr1"] * F, "start": starts,
                                         "end_": starts + 500, "peak_id": [f"g{j}" for j in range(F)]}))
    ds.flush()
    w = ds.create_layer_writer("ATAC_counts", N, F, dtype=np.float32, compression="zlib")
    w.write_chunk(X, 0)
    w.finalize()
    ds.flush()
    if fb == "column_major":
        from test_per_group_floor import _column_major
        _column_major(ds, X, "ATAC_counts", width=7)
    ds.close()
    if fb != "column_major":
        ds = cytome.open(p, mode="r")
        try:
            zeros = sum(int((c.data == 0).sum()) for c, _ in ds.iter_chunks("ATAC", layer="counts"))
        finally:
            ds.close()
        assert zeros == 60                                     # the store keeps the stored zeros
    def run(pct):
        return cosg.run_cosg_cytome(p, groupby="grp", modality="ATAC", layer="counts", n_genes_user=F,
                                    expressed_pct=pct, feature_batching="none" if fb == "none" else "chromosome",
                                    output_format="dict", iqr_normalize=False, verbose=False)["scores_dict"]
    assert ("A", "g5") not in run(0.5)                          # 21 of 89 cells: filtered
    assert ("A", "g5") in run(0.2)                              # and kept under a floor it meets
