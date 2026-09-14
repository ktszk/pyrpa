import ctypes
import numpy as np
import libs.flibs._loader as _loader

# Public API of this submodule; keeps `from ._x import *` in __init__.py from
# re-exporting numpy, scipy, ctypes and the sibling packages imported above.
__all__ = [
    'matrix_riccati_batch', 'matrix_riccati_chords', 'riccati_chords', 'riccati_chords_bc'
]

# --- ctypes signature (set once at import); the Fortran entry is a subroutine.
_loader._lib.riccati_chords.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),           # g   (out) [Ns,Nchord,Nw]
    np.ctypeslib.ndpointer(dtype=np.complex128),           # f   (out)
    np.ctypeslib.ndpointer(dtype=np.complex128),           # om  (in)
    np.ctypeslib.ndpointer(dtype=np.complex128),           # dd  (in)
    ctypes.POINTER(ctypes.c_double),                                     # hvf
    np.ctypeslib.ndpointer(dtype=np.float64),              # ds  (in) [Nchord]
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),  # Ns, Nchord, Nw
]
_loader._lib.riccati_chords.restype = None


_loader._lib.matrix_riccati_batch.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),              # g   (out) [Ns,Nw,2,2]
    np.ctypeslib.ndpointer(dtype=np.complex128),              # f   (out)
    np.ctypeslib.ndpointer(dtype=np.complex128),              # om  (in)  [Nw]
    np.ctypeslib.ndpointer(dtype=np.complex128),              # Dpath (in) [Ns,2,2]
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),  # hvf, ds, h
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),                       # Ns, Nw
]
_loader._lib.matrix_riccati_batch.restype = None


def matrix_riccati_batch(om: np.ndarray, Dpath: np.ndarray, hvf: float, ds: float,
                         h: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    """
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
    _loader._lib.matrix_riccati_batch(g, f, om, Dpath, _loader.dbl(hvf), _loader.dbl(ds), _loader.dbl(h),
                                      _loader.i64(Ns), _loader.i64(Nw))
    return g, f


_loader._lib.matrix_riccati_chords.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),              # g   (out) [Ns,Nchord,Nw,2,2]
    np.ctypeslib.ndpointer(dtype=np.complex128),              # f   (out)
    np.ctypeslib.ndpointer(dtype=np.complex128),              # om  (in)  [Nw]
    np.ctypeslib.ndpointer(dtype=np.complex128),              # Dpath (in) [Ns,Nchord,2,2]
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),  # hvf, ds, h
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),     # Ns, Nchord, Nw
]
_loader._lib.matrix_riccati_chords.restype = None


def matrix_riccati_chords(om: np.ndarray, Dpath: np.ndarray, hvf: float, ds: float,
                          h: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    """
    @brief 2x2 spin-matrix g, f along many chords, batched over frequencies (Fortran,
    OpenMP over chord x frequency).  For the 2D vortex d-vector solver: one call per FS
    direction.  omega may be position dependent (Doppler v_F.Q on the vortex lattice).
    @param om: frequencies [Nw] (constant along the chord) OR [Ns, Nchord, Nw]
    @param Dpath: 2x2 gap matrix along each chord [Ns, Nchord, 2, 2] complex128
    @param hvf, ds, h: hbar|v_F|, arc-length step, Zeeman energy
    @return (g, f): 2x2 propagators [Ns, Nchord, Nw, 2, 2] complex128
    """
    Dpath = np.ascontiguousarray(Dpath, dtype=np.complex128)
    Ns, Nchord = Dpath.shape[0], Dpath.shape[1]
    om = np.asarray(om, dtype=np.complex128)
    Nw = om.shape[-1]
    if om.ndim == 1:                                   # broadcast constant-omega to the path
        om = np.broadcast_to(om, (Ns, Nchord, Nw))
    om = np.ascontiguousarray(om)
    g = np.empty((Ns, Nchord, Nw, 2, 2), dtype=np.complex128)
    f = np.empty((Ns, Nchord, Nw, 2, 2), dtype=np.complex128)
    _loader._lib.matrix_riccati_chords(g, f, om, Dpath, _loader.dbl(hvf), _loader.dbl(ds), _loader.dbl(h),
                                       _loader.i64(Ns), _loader.i64(Nchord), _loader.i64(Nw))
    return g, f


def riccati_chords(om: np.ndarray, dd: np.ndarray, hvf: float, ds: float) -> tuple[np.ndarray, np.ndarray]:
    """
    @brief Scalar quasiclassical g, f along many chords via the stable tanh-step
    Riccati (Fortran).  The single batched chord kernel for the surface / vortex /
    lattice solvers (forward gamma + backward gamma-tilde + g,f combine).
    @param om: (renormalized) frequency along each chord [Ns, Nchord, Nw] complex128
    @param dd: order parameter along each chord [Ns, Nchord, Nw] complex128
    @param hvf: hbar |v_F|
    @param ds: arc-length step, scalar or per-chord [Nchord] (e.g. dx/|cos beta|)
    @return (g, f): propagators [Ns, Nchord, Nw] complex128
    """
    om = np.ascontiguousarray(om, dtype=np.complex128)
    dd = np.ascontiguousarray(dd, dtype=np.complex128)
    Ns, Nchord, Nw = om.shape
    ds_arr = (np.full(Nchord, float(ds)) if np.isscalar(ds)
              else np.ascontiguousarray(ds, dtype=np.float64))
    g = np.empty((Ns, Nchord, Nw), dtype=np.complex128)
    f = np.empty((Ns, Nchord, Nw), dtype=np.complex128)
    _loader._lib.riccati_chords(g, f, om, dd, _loader.dbl(hvf), ds_arr,
                                _loader.i64(Ns), _loader.i64(Nchord), _loader.i64(Nw))
    return g, f


_loader._lib.riccati_chords_bc.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),           # g   (out)
    np.ctypeslib.ndpointer(dtype=np.complex128),           # f   (out)
    np.ctypeslib.ndpointer(dtype=np.complex128),           # om  (in)  [Nw]
    np.ctypeslib.ndpointer(dtype=np.complex128),           # dom (in)  [Ns,Nchord]
    np.ctypeslib.ndpointer(dtype=np.complex128),           # dd  (in)  [Ns,Nchord]
    ctypes.POINTER(ctypes.c_double),                                     # hvf
    np.ctypeslib.ndpointer(dtype=np.float64),              # ds  (in)  [Nchord]
    ctypes.POINTER(ctypes.c_int64),                                      # irow
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),  # Ns, Nchord, Nw
]
_loader._lib.riccati_chords_bc.restype = None


def riccati_chords_bc(om: np.ndarray, dd: np.ndarray, hvf: float, ds: float, dom: np.ndarray | None = None,
                      row: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """
    @brief Scalar Riccati chord integration with the frequency and the space
    dependence kept separate (Fortran): om_local = om[w] + dom[i,c] and a
    frequency-independent gap dd[i,c].  That is the structure of the clean vortex
    and vortex-lattice problem, so the caller avoids building two [Ns,Nchord,Nw]
    arrays per direction per iteration (tens of MB) only to hand them over --
    which is what dominated the solvers.  Use riccati_chords when the gap itself
    is frequency dependent (impurity self-energy Sigma_f).
    @param   om: frequency [Nw] complex128
    @param   dd: order parameter along each chord [Ns, Nchord] complex128
    @param  hvf: hbar |v_F|
    @param   ds: arc-length step, scalar or per-chord [Nchord]
    @param  dom: per-point frequency shift [Ns, Nchord] (Doppler / gauge), or None
    @param  row: None = return the full [Ns, Nchord, Nw] propagators;
                 an int i = return only that point of each chord, [Nchord, Nw]
                 (the anchor the gap equation and the lattice DOS need)
    @return (g, f)
    """
    om = np.ascontiguousarray(om, dtype=np.complex128).ravel()
    dd = np.ascontiguousarray(dd, dtype=np.complex128)
    Ns, Nchord = dd.shape
    Nw = om.size
    dom = (np.zeros((Ns, Nchord), dtype=np.complex128) if dom is None
           else np.ascontiguousarray(dom, dtype=np.complex128))
    ds_arr = (np.full(Nchord, float(ds)) if np.isscalar(ds)
              else np.ascontiguousarray(ds, dtype=np.float64))
    irow = 0 if row is None else int(row) + 1          # Fortran is 1-based
    shape = (Nchord, Nw) if row is not None else (Ns, Nchord, Nw)
    g = np.empty(shape, dtype=np.complex128)
    f = np.empty(shape, dtype=np.complex128)
    _loader._lib.riccati_chords_bc(g, f, om, dom, dd, _loader.dbl(hvf), ds_arr, _loader.i64(irow),
                                   _loader.i64(Ns), _loader.i64(Nchord), _loader.i64(Nw))
    return g, f
