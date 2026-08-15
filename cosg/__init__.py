"""COSG — COSine similarity-based marker Gene detection.

Public API
----------

``cosg.cosg(data, ...)`` is the recommended unified entry point.
It dispatches on the first argument:

- ``AnnData``                       → in-memory CPU/GPU path (writes
                                       ``data.uns[key_added]``).
- ``str`` / ``pathlib.Path``        → cytome streaming, bounded-RAM
                                       (returns a dict).
- ``CytomeDataset`` (open object,   → cytome streaming, same as above;
  as returned by ``cytome.open()``)    the caller owns the Dataset's
                                       lifecycle (it is NOT closed).

Power-user shortcuts (skip the polymorphic dispatch):

- ``cosg.run_cosg_cytome(path_or_dataset, ..., device='cpu')``
- ``cosg.run_cosg_cytome_gpu(path_or_dataset, ...)`` — thin shim for
  ``run_cosg_cytome(..., device='gpu')``.

GPU selection works on both paths via ``device='cpu' | 'gpu' | 'auto'``.
``device='auto'`` resolves to GPU when CuPy is available, else CPU.

Renames in the cytome streaming API:
  * ``run_cosg_cytome_cpu`` → ``run_cosg_cytome``. Old name imports
    succeed (a stub function is exported for clearer error messages),
    but **calling** it raises ``TypeError`` with a migration hint.
  * ``cytome_layer=`` kwarg → ``layer=``. Old kwarg raises ``TypeError``.
"""

from .cosg import cosg, indexByGene, iqrLogNormalize
from ._backend import HAS_CUPY
from ._backend import HAS_CUPY as _HAS_GPU

# Cytome-streaming variants — bounded-RAM marker calls reading directly from a
# .cytome file. Imported unconditionally: cytome is a required dependency, and
# _cytome_streaming no longer imports it at module scope, so this costs nothing
# at import time.
#
# This used to sit in a try/except. That made a missing cytome delete the names
# rather than explain itself -- callers got `AttributeError: module 'cosg' has
# no attribute 'run_cosg_cytome'`, which reads like a typo. Now the functions
# always exist and raise an ImportError naming the package if it is absent.
from ._cytome_streaming import (
    run_cosg_cytome,            # Round 9: canonical name
    run_cosg_cytome_cpu,        # Round 9: K-HARD migration stub
    run_cosg_cytome_gpu,        # shim → run_cosg_cytome(device='gpu')
)

try:
    from ._plotting import plotMarkerDendrogram, plotMarkerDotplot, plotMarkerStream
except ImportError:
    pass

__version__ = "1.1.1"

__all__ = [
    'cosg',
    'indexByGene',
    'iqrLogNormalize',
    # Published since 1.0.3 and imported above (guarded, because _plotting is
    # not in this working tree). They must stay in __all__: the public 1.0.4
    # __init__ defines no __all__ at all, so introducing one that omits them
    # would silently stop `from cosg import *` exporting them.
    'plotMarkerDendrogram',
    'plotMarkerDotplot',
    'plotMarkerStream',
    'run_cosg_cytome',
    'run_cosg_cytome_gpu',
    # run_cosg_cytome_cpu is intentionally NOT in __all__ — it's a
    # migration-aid stub that raises TypeError at call time. Listed
    # here as a comment so removal is a conscious decision.
]
