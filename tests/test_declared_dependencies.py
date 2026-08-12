"""Every module-level third-party import must be a declared dependency.

COSG 1.0.4 depended on ``scanpy>=1.6.0``, and scanpy pulled in numba
transitively. A working tree without ``_plotting.py`` looks like it does not
need scanpy — it survives only in a docstring example — but ``cosg/cosg.py``
still does a module-level ``import numba`` to compile
``sparse_mean_var_{minor,major}_axis``. Dropping scanpy without declaring
numba made ``import cosg`` fail on a clean ``pip install cosg`` — the same
class of bug as the undeclared lz4 in cytome, and invisible on any developer
machine, where scanpy/umap have already dragged numba in.

Both are now declared. numba is declared *directly* rather than left to
scanpy, because that is exactly the transitive reliance that broke.

The parametrised check below is the general guard; the rest pin the specific
regressions.
"""
from __future__ import annotations

import ast
import importlib
import pathlib
import sys

import numpy as np
import pytest

PKG = pathlib.Path(__file__).resolve().parent.parent / "cosg"

#: Imported at module level and therefore required for the full public API.
#: networkx is needed by _plotting; matplotlib by both.
REQUIRED = [
    "anndata", "numpy", "pandas", "sklearn", "scipy", "numba",
    "matplotlib", "networkx",
]

#: Guarded by try/except or imported lazily inside a function, so a plain
#: install may legitimately lack them. scanpy is the ``cosg[dotplot]`` extra —
#: see test_scanpy_optional.py.
OPTIONAL = {"cupy", "cupyx", "cytome", "piaso", "scanpy"}


@pytest.mark.parametrize("mod", REQUIRED)
def test_required_dependency_is_importable(mod):
    """If this fails, `mod` has fallen out of [project] dependencies."""
    assert importlib.import_module(mod) is not None


def test_no_undeclared_module_level_imports():
    """Catch the next numba: a new top-level import of an undeclared package.

    Only module-level imports matter — a lazy import inside a function (the
    GPU path's ``import cupy``) fails at call time, which is recoverable and
    intended.
    """
    known = set(REQUIRED) | OPTIONAL | {"cosg"} | set(sys.stdlib_module_names)
    offenders = []
    for path in sorted(PKG.glob("*.py")):
        tree = ast.parse(path.read_text())
        for node in tree.body:  # module level only, not ast.walk
            names = []
            if isinstance(node, ast.Import):
                names = [a.name.split(".")[0] for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                names = [node.module.split(".")[0]]
            offenders += [
                f"{path.name}:{node.lineno} -> {n}" for n in names if n not in known
            ]
    assert not offenders, (
        "module-level imports of undeclared packages:\n  "
        + "\n  ".join(offenders)
        + "\nAdd them to [project] dependencies (and to REQUIRED here), or make "
        "the import lazy/guarded."
    )


def test_numba_kernels_actually_run():
    """numba is not merely importable: these kernels are jitted at import."""
    from scipy import sparse

    from cosg.cosg import sparse_mean_var_major_axis

    X = sparse.csr_matrix(np.array([[1.0, 0.0, 3.0], [0.0, 2.0, 0.0]]))
    means, variances = sparse_mean_var_major_axis(
        X.data, X.indices, X.indptr, X.shape[0], X.shape[1], np.dtype(np.float64)
    )
    np.testing.assert_allclose(means, X.toarray().mean(axis=1))
    # ddof=0: the kernel returns the population variance and leaves the
    # (n / (n - 1)) correction to its caller.
    np.testing.assert_allclose(variances, X.toarray().var(axis=1), rtol=1e-6)


def test_plotting_api_is_actually_exported():
    """The three plotMarker* names in __all__ must really be importable.

    __init__ imports _plotting inside a try/except ImportError, so anything
    that breaks that import — a missing _plotting.py, or dropping scanpy,
    which it imports at module level — removes the functions *silently* while
    __all__ still advertises them. `from cosg import *` then raises and
    `cosg.plotMarkerDotplot` is an AttributeError, for a package whose release
    notes promise this API is unchanged.
    """
    import cosg

    for name in ("plotMarkerDendrogram", "plotMarkerDotplot", "plotMarkerStream"):
        assert name in cosg.__all__
        assert hasattr(cosg, name), (
            f"{name} is in __all__ but not exported — _plotting failed to "
            "import. Check that cosg/_plotting.py is present and that its "
            "module-level imports (scanpy, networkx, matplotlib) are declared."
        )
