from ctypes import byref, c_int64, c_double
from pathlib import Path
import numpy as np

# Resolve lib path relative to this file so import works from any current working directory.
_LIB_DIR = Path(__file__).resolve().parents[1]
_lib = np.ctypeslib.load_library("libfmod.so", str(_LIB_DIR))


def i64(value):
    """byref-wrapped 64-bit integer scalar for ctypes Fortran calls."""
    return byref(c_int64(value))


def dbl(value):
    """byref-wrapped double scalar for ctypes Fortran calls."""
    return byref(c_double(value))
