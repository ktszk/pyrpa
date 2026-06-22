from ctypes import POINTER, c_int64, c_double
import numpy as np
from ._loader import _lib, i64, dbl

# --- ctypes signature (set once at import); the Fortran entry is a subroutine.
_lib.riccati_chords.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),   # g   (out) [Ns,Nchord,Nw]
    np.ctypeslib.ndpointer(dtype=np.complex128),   # f   (out)
    np.ctypeslib.ndpointer(dtype=np.complex128),   # om  (in)
    np.ctypeslib.ndpointer(dtype=np.complex128),   # dd  (in)
    POINTER(c_double), POINTER(c_double),          # hvf, ds
    POINTER(c_int64), POINTER(c_int64), POINTER(c_int64),   # Ns, Nchord, Nw
]
_lib.riccati_chords.restype = None


def riccati_chords(om: np.ndarray, dd: np.ndarray, hvf: float, ds: float):
    """
    @fn riccati_chords
    @brief Scalar quasiclassical g, f along many chords via the stable tanh-step
    Riccati (Fortran).  Drop-in batched replacement for the per-chord Python
    _chord_gf / _chord_gf_pos loops in the surface / vortex / lattice solvers.
    @param om: (renormalized) frequency along each chord [Ns, Nchord, Nw] complex128
    @param dd: order parameter along each chord [Ns, Nchord, Nw] complex128
    @param hvf: hbar |v_F|;  ds: arc-length step
    @return (g, f): propagators [Ns, Nchord, Nw] complex128
    """
    om = np.ascontiguousarray(om, dtype=np.complex128)
    dd = np.ascontiguousarray(dd, dtype=np.complex128)
    Ns, Nchord, Nw = om.shape
    g = np.empty((Ns, Nchord, Nw), dtype=np.complex128)
    f = np.empty((Ns, Nchord, Nw), dtype=np.complex128)
    _lib.riccati_chords(g, f, om, dd, dbl(hvf), dbl(ds),
                        i64(Ns), i64(Nchord), i64(Nw))
    return g, f
