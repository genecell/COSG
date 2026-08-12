"""
Backend detection for COSG GPU support.

Lazily detects CuPy availability at import time.
No hard dependency — if CuPy is absent, GPU features gracefully degrade to CPU.
"""

HAS_CUPY = False
_cupy_import_error = None

try:
    import cupy as cp
    import cupyx.scipy.sparse as cp_sparse
    HAS_CUPY = True
except ImportError as e:
    _cupy_import_error = e


def require_gpu():
    """Raise ImportError with helpful message if CuPy is not available."""
    if not HAS_CUPY:
        raise ImportError(
            "GPU acceleration requires CuPy. Install with:\n"
            "  pip install cosg[gpu]        # if using pip extras\n"
            "  pip install cupy-cuda12x     # for CUDA 12.x\n"
            "  pip install cupy-cuda11x     # for CUDA 11.x\n"
            f"\nOriginal import error: {_cupy_import_error}"
        )


def get_device(device='auto'):
    """
    Resolve device string to actual backend.
    
    Parameters
    ----------
    device : str
        'auto' : use GPU if CuPy available, else CPU
        'gpu'  : require GPU, raise if unavailable
        'cpu'  : force CPU even if GPU is available
    
    Returns
    -------
    str : 'gpu' or 'cpu'
    """
    if device == 'auto':
        return 'gpu' if HAS_CUPY else 'cpu'
    elif device == 'gpu':
        require_gpu()
        return 'gpu'
    elif device == 'cpu':
        return 'cpu'
    else:
        raise ValueError(f"device must be 'auto', 'gpu', or 'cpu', got '{device}'")
