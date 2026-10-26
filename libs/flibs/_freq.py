import ctypes
import numpy as np
import libs.flibs._loader as _loader

# Public API of this submodule; keeps `from ._x import *` in __init__.py from
# re-exporting numpy, ctypes and the sibling packages imported above.
__all__ = ['fg_set_ir', 'fg_set_matsubara']

_loader._lib.fg_set_ir.argtypes = [
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),   # Nw, Ntau
    _loader.c_array(np.float64), _loader.c_array(np.float64),         # wF, wB [Nw]
    _loader.f_array(np.complex128), _loader.f_array(np.complex128),   # fwdF (Ntau,2Nw), bwdF (2Nw,Ntau)
    _loader.f_array(np.complex128), _loader.f_array(np.complex128),   # fwdB, bwdB
    _loader.c_array(np.int64),                                        # rev [Ntau] (1-based)
    _loader.c_array(np.complex128), _loader.c_array(np.complex128),   # stat, dens [2Nw]
]
_loader._lib.fg_set_ir.restype = None
_loader._lib.fg_set_matsubara.argtypes = []
_loader._lib.fg_set_matsubara.restype = None


def fg_set_ir(g: dict) -> None:
    """
    @brief Switch the Fortran frequency grid (module freq_grid) to the IR basis built by plibs.ir_grid.
    Every later FFT-convolution kernel call (chi0, Sigma) then works on its Nw stored sampling
    frequencies; Green's functions are built at those frequencies.
    @param g: the dict of plibs.ir_grid
    """
    _loader._lib.fg_set_ir(_loader.i64(g['Nw']), _loader.i64(g['Ntau']),
                           np.ascontiguousarray(g['wF']), np.ascontiguousarray(g['wB']),
                           g['fwdF'], g['bwdF'], g['fwdB'], g['bwdB'],
                           np.ascontiguousarray(g['rev'], dtype=np.int64),
                           np.ascontiguousarray(g['stat']), np.ascontiguousarray(g['dens']))


def fg_set_matsubara() -> None:
    """@brief Back to the uniform Matsubara grid (the default)."""
    _loader._lib.fg_set_matsubara()
