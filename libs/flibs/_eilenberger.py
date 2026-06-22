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


_lib.matrix_riccati_batch.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),   # g   (out) [Ns,Nw,2,2]
    np.ctypeslib.ndpointer(dtype=np.complex128),   # f   (out)
    np.ctypeslib.ndpointer(dtype=np.complex128),   # om  (in)  [Nw]
    np.ctypeslib.ndpointer(dtype=np.complex128),   # Dpath (in) [Ns,2,2]
    POINTER(c_double), POINTER(c_double), POINTER(c_double),   # hvf, ds, h
    POINTER(c_int64), POINTER(c_int64),            # Ns, Nw
]
_lib.matrix_riccati_batch.restype = None


def matrix_riccati_batch(om: np.ndarray, Dpath: np.ndarray, hvf: float, ds: float, h: float = 0.0):
    """
    @fn matrix_riccati_batch
    @brief 2x2 spin-matrix quasiclassical g, f along one inhomogeneous trajectory,
    batched over frequencies (Fortran).  Drop-in replacement for looping
    matrix_trajectory_gf over frequencies in the d-vector solvers.
    @param om: (renormalized) frequencies [Nw] complex128
    @param Dpath: 2x2 gap matrix along the path [Ns, 2, 2] complex128
    @param hvf, ds, h: hbar|v_F|, arc-length step, Zeeman energy
    @return (g, f): 2x2 propagators [Ns, Nw, 2, 2] complex128
    """
    om = np.ascontiguousarray(om, dtype=np.complex128)
    Dpath = np.ascontiguousarray(Dpath, dtype=np.complex128)
    Ns = Dpath.shape[0]
    Nw = om.shape[0]
    g = np.empty((Ns, Nw, 2, 2), dtype=np.complex128)
    f = np.empty((Ns, Nw, 2, 2), dtype=np.complex128)
    _lib.matrix_riccati_batch(g, f, om, Dpath, dbl(hvf), dbl(ds), dbl(h),
                              i64(Ns), i64(Nw))
    return g, f


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
