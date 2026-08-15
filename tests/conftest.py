"""Shared markers for optional-dependency test gating.

COSG's core (marker detection on an AnnData) has no optional dependencies.
Two features do, and their tests must *skip* rather than error when the
dependency is absent — otherwise a contributor who runs `pip install -e .`
and then pytest sees a wall of ModuleNotFoundError that looks like broken
code rather than a missing extra:

- ``cytome`` — the streaming backend (``cosg[dev]`` installs it)
- ``scanpy`` — ``plotMarkerDotplot`` only (``cosg[dotplot]``)

The release workflow asserts both are importable *before* running pytest, so
these skips cannot silently hide a broken extra in CI.
"""
from __future__ import annotations

import importlib.util

import pytest


def _missing(module: str) -> bool:
    return importlib.util.find_spec(module) is None


requires_cytome = pytest.mark.skipif(
    _missing("cytome"),
    reason="cytome not installed — needed for the .cytome streaming backend "
           "(pip install cytome, or pip install -e '.[dev]')",
)

#: piaso IS on PyPI — as the distribution `piaso-tools` (the import name and
#: the distribution name differ). It still cannot be used in CI, for two
#: reasons: it declares `cosg>=1.0.3`, so depending on it would be circular;
#: and the published piaso-tools 1.1.0 does not yet contain the streaming
#: normalizer symbols COSG's cytome path imports, so installing it would not
#: make these tests pass. They therefore run locally, where piaso is on
#: PYTHONPATH, and skip in CI. The layer='counts' path needs no piaso and is
#: covered everywhere.
requires_piaso = pytest.mark.skipif(
    _missing("piaso"),
    reason="piaso not installed — needed for on-the-fly log1p/infog/tfidf "
           "normalization (not on PyPI; install from source)",
)

requires_scanpy = pytest.mark.skipif(
    _missing("scanpy"),
    reason="scanpy not installed — needed for plotMarkerDotplot "
           "(pip install 'cosg[dotplot]')",
)
