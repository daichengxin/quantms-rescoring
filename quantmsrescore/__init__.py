# =============================================================================
# Thread Control Configuration (MUST be done before importing heavy libraries)
# =============================================================================
# This module controls thread allocation for numerical libraries to prevent
# thread explosion in HPC/Slurm environments where multiple processes compete
# for resources.
#
# Problem: Libraries like PyTorch, TensorFlow, NumPy (via MKL/OpenBLAS) spawn
# threads equal to CPU count by default. With multiprocessing, this causes:
#   - processes × cpu_count threads competing for cpu_count cores
#   - Memory pressure from thread stacks
#   - Node OOM kills in Slurm when co-located with other jobs
#
# Solution: Set environment variables BEFORE importing heavy libraries.
# =============================================================================

import os
from typing import Optional

# Default threads per worker process
_DEFAULT_THREADS_PER_PROCESS = 1


def configure_threading(n_threads: Optional[int] = None, verbose: bool = False) -> None:
    """
    Configure thread counts for all numerical/ML libraries.

    This function MUST be called before importing numpy, torch, tensorflow,
    or any library that depends on them. It sets environment variables that
    control internal threading in these libraries.

    Parameters
    ----------
    n_threads : int, optional
        Number of threads per process. Defaults to 1 for HPC safety.
        Set higher only if you're sure about available resources.
    verbose : bool, optional
        If True, log the thread configuration.

    Notes
    -----
    For Slurm/HPC environments, the recommended approach is:
        - Set n_threads=1 (default)
        - Use --processes to control parallelism at the process level
        - Total parallelism = processes × n_threads

    Example
    -------
    For a 32-core node with --processes 8:
        - n_threads=1: 8 processes × 1 thread = 8 parallel units (safe)
        - n_threads=4: 8 processes × 4 threads = 32 parallel units (full utilization)
    """
    if n_threads is None:
        n_threads = _DEFAULT_THREADS_PER_PROCESS

    n_threads_str = str(n_threads)

    # OpenMP (used by many scientific libraries)
    os.environ["OMP_NUM_THREADS"] = n_threads_str

    # Intel MKL (NumPy/SciPy on Intel, PyTorch default backend)
    os.environ["MKL_NUM_THREADS"] = n_threads_str

    # NumExpr (used by pandas for some operations)
    os.environ["NUMEXPR_MAX_THREADS"] = n_threads_str

    # OpenBLAS (NumPy/SciPy alternative backend)
    os.environ["OPENBLAS_NUM_THREADS"] = n_threads_str

    # Apple Accelerate (macOS)
    os.environ["VECLIB_MAXIMUM_THREADS"] = n_threads_str

    # BLIS (another BLAS implementation)
    os.environ["BLIS_NUM_THREADS"] = n_threads_str

    # TensorFlow specific
    os.environ["TF_NUM_INTEROP_THREADS"] = n_threads_str
    os.environ["TF_NUM_INTRAOP_THREADS"] = n_threads_str

    # Disable TensorFlow GPU memory pre-allocation (helps with shared nodes)
    os.environ.setdefault("TF_FORCE_GPU_ALLOW_GROWTH", "true")

    # Reduce TensorFlow verbosity
    os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")

    if verbose:
        print(f"[quantms-rescoring] Thread configuration: {n_threads} threads per process")


def configure_torch_threads(n_threads: Optional[int] = None) -> None:
    """
    Configure PyTorch-specific thread settings.

    This should be called after PyTorch is imported, as it uses PyTorch's
    API directly rather than environment variables.

    Parameters
    ----------
    n_threads : int, optional
        Number of threads for PyTorch operations. Defaults to 1.
    """
    if n_threads is None:
        n_threads = _DEFAULT_THREADS_PER_PROCESS

    try:
        import torch
        torch.set_num_threads(n_threads)
        torch.set_num_interop_threads(n_threads)
    except ImportError:
        pass  # PyTorch not installed
    except RuntimeError:
        # Threads already set (can only be set once)
        pass


def get_safe_process_count(requested: int, memory_per_process_gb: float = 4.0) -> int:
    """
    Calculate a safe number of processes based on available resources.

    Parameters
    ----------
    requested : int
        Requested number of processes.
    memory_per_process_gb : float, optional
        Estimated memory per process in GB. Default is 4 GB.

    Returns
    -------
    int
        Safe number of processes (min of requested and calculated safe limit).
    """
    import multiprocessing

    cpu_count = multiprocessing.cpu_count()

    # Try to get available memory
    try:
        import psutil
        available_memory_gb = psutil.virtual_memory().available / (1024 ** 3)
        memory_based_limit = max(1, int(available_memory_gb / memory_per_process_gb))
    except ImportError:
        memory_based_limit = cpu_count  # Fall back to CPU count if psutil unavailable

    safe_count = min(requested, cpu_count, memory_based_limit)

    return max(1, safe_count)


# =============================================================================
# Apply default thread configuration immediately on import
# =============================================================================
# This ensures that any subsequent imports of numpy, torch, etc. respect
# the thread limits. Users can call configure_threading() again to adjust.
configure_threading(n_threads=_DEFAULT_THREADS_PER_PROCESS)


# =============================================================================
# Standard module initialization
# =============================================================================
from warnings import filterwarnings

# Suppress warnings about OPENMS_DATA_PATH
filterwarnings("ignore", message=".*OPENMS_DATA_PATH.*", category=UserWarning)

__version__ = "0.0.14"

__all__ = [
    "configure_threading",
    "configure_torch_threads",
    "get_safe_process_count",
    "__version__",
]
