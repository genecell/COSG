"""A groupby that does not exist must fail, not invent a group.

SQLite treats a double-quoted token that does not resolve to a column as a
*string literal* rather than an error:

    SELECT cell_idx, "Leiden" FROM cells      -- no Leiden column
    -> [(0, 'Leiden'), (1, 'Leiden'), ...]

so every cell got the same label, COSG saw one group named after the typo, and
the marker table that came out looked entirely ordinary. The log said
``1 groups`` and nothing else complained. That is a confident wrong answer,
which is worse than a traceback.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

cytome = pytest.importorskip("cytome")
import cosg


def test_sqlite_really_does_this(tmp_path):
    """The premise, in isolation. If SQLite ever stops, the guard can go."""
    import sqlite3

    con = sqlite3.connect(":memory:")
    con.execute("CREATE TABLE cells (cell_idx INTEGER, ATAC_leiden TEXT)")
    con.executemany("INSERT INTO cells VALUES (?,?)", [(i, str(i % 3)) for i in range(6)])
    rows = con.execute('SELECT cell_idx, "Leiden" FROM cells').fetchall()
    assert {r[1] for r in rows} == {"Leiden"}, (
        "SQLite no longer returns an unresolved quoted identifier as a string "
        "literal; the guard in _get_labels_and_gene_names can be revisited")


@pytest.fixture
def ds(tmp_path):
    rng = np.random.default_rng(0)
    n, g = 60, 20
    a = pytest.importorskip("anndata").AnnData(
        sp.csr_matrix(rng.poisson(2, (n, g)).astype(np.float32)),
        obs=pd.DataFrame({"ATAC_leiden": rng.choice(list("012"), n),
                          "CellTypes": rng.choice(["A", "B"], n)},
                         index=[f"c{i}" for i in range(n)]),
        var=pd.DataFrame(index=[f"g{i}" for i in range(g)]))
    d = cytome.from_anndata(a, output=str(tmp_path / "t.cytome"))
    yield d
    d.close()


def test_a_missing_groupby_raises_and_names_the_columns(ds):
    with pytest.raises(KeyError) as e:
        cosg.run_cosg_cytome(ds, groupby="Leiden", n_genes_user=5, verbose=False)
    msg = str(e.value)
    assert "Leiden" in msg
    assert "ATAC_leiden" in msg, "the message should list what is actually there"


def test_a_case_variant_still_resolves(ds):
    """SQLite column names are case-insensitive, so `atac_leiden` really does
    name `ATAC_leiden` and must not be turned away by the guard."""
    out = cosg.run_cosg_cytome(ds, groupby="atac_leiden", n_genes_user=5,
                               output_format="long", verbose=False)
    groups = set(out["scores_df"]["group"])
    assert groups <= {"0", "1", "2"} and len(groups) > 1, groups
    assert "atac_leiden" not in groups


def test_the_real_column_still_works_and_gives_its_own_groups(ds):
    out = cosg.run_cosg_cytome(ds, groupby="ATAC_leiden", n_genes_user=5,
                               output_format="long", verbose=False)
    groups = set(out["scores_df"]["group"])
    assert groups <= {"0", "1", "2"} and len(groups) > 1, groups
    assert "ATAC_leiden" not in groups, "the column name is not a group label"
