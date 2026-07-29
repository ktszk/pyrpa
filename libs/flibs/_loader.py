from ctypes import byref, c_int64, c_double
from pathlib import Path
import numpy as np

# Resolve lib path relative to this file so import works from any current working directory.
_LIB_DIR = Path(__file__).resolve().parents[1]


def _cpu_has_avx512():
    """True if the running CPU supports AVX-512F (Zen 4, Intel AVX-512, ...)."""
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.startswith("flags"):
                    return "avx512f" in line.split()
    except OSError:
        pass
    return False


def _load_library():
    """Pick the fastest libfmod.so the current CPU can run.

    `make dispatch` builds per-arch variants (libfmod_znver4.so = AVX-512,
    libfmod_znver3.so = AVX2). Prefer the AVX-512 build only when the CPU
    supports it, so one checkout runs on both Zen 3 (e.g. 5700G) and Zen 4
    (e.g. 7900) yet uses AVX-512 only where it is safe. Fall back to the
    generic single-arch libfmod.so produced by a plain `make`.
    """
    candidates = []
    if _cpu_has_avx512():
        candidates.append("libfmod_znver4")
    candidates.append("libfmod_znver3")
    candidates.append("libfmod")  # generic fallback (plain `make`)
    for name in candidates:
        if (_LIB_DIR / f"{name}.so").exists():
            return np.ctypeslib.load_library(name, str(_LIB_DIR))
    # Nothing matched by name; let load_library raise a clear error for libfmod.
    return np.ctypeslib.load_library("libfmod", str(_LIB_DIR))


_lib = _load_library()


def i64(value):
    """byref-wrapped 64-bit integer scalar for ctypes Fortran calls."""
    return byref(c_int64(value))


def dbl(value):
    """byref-wrapped double scalar for ctypes Fortran calls."""
    return byref(c_double(value))
