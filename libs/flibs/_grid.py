"""
Grid kernels for the vortex / vortex-lattice solvers (Fortran, libs/src/fgrid.f90).

Bilinear sampling of small 2D fields at many scattered chord points, and the
Abrikosov lowest-Landau-level order parameter.  These replace the scipy
RegularGridInterpolator / numpy fancy-indexing versions that dominated the
runtime of the 2D vortex and periodic-lattice solvers, and are numerically
identical to them (same bilinear weights, same out-of-range fill).
"""
from ctypes import POINTER, c_int64, c_double
import numpy as np
from ._loader import _lib, i64, dbl

_c128 = np.ctypeslib.ndpointer(dtype=np.complex128)
_f64 = np.ctypeslib.ndpointer(dtype=np.float64)

_lib.bilinear_eval.argtypes = [
    _c128,                                   # out  [npt,nw]
    _c128,                                   # fld  [ng,ng,nw]
    POINTER(c_double), POINTER(c_double), POINTER(c_double),   # x0, dx, xhi
    POINTER(c_int64),                        # ng
    _f64, _f64,                              # px, py [npt]
    POINTER(c_int64), POINTER(c_int64),      # npt, nw
    _c128,                                   # fill (1-element array: Fortran takes it by reference)
]
_lib.bilinear_eval.restype = None

_lib.bilinear_cell_eval.argtypes = [
    _c128, _c128,                            # out [npt], fld [ng,ng]
    POINTER(c_int64),                        # ng
    _f64, _f64, POINTER(c_int64),            # px, py, npt
    POINTER(c_double), POINTER(c_double), POINTER(c_double),   # Dlx, Dly, zeta
]
_lib.bilinear_cell_eval.restype = None

_lib.abrikosov_z.argtypes = [
    _c128, _f64, _f64, POINTER(c_int64),     # out, px, py, npt
    POINTER(c_double), POINTER(c_double), POINTER(c_double),   # Dlx, Dly, zeta
    POINTER(c_int64),                        # nsum
]
_lib.abrikosov_z.restype = None


def _cplx(value):
    """complex128 scalar as a 1-element array (ctypes has no complex scalar; a
    Fortran scalar dummy takes the same by-reference pointer)."""
    return np.array([value], dtype=np.complex128)


def bilinear_eval(field: np.ndarray, xg: np.ndarray, px: np.ndarray, py: np.ndarray,
                  fill=0.0) -> np.ndarray:
    """
    @fn bilinear_eval
    @brief Bilinear interpolation of a field on a uniform square grid at scattered
    points, vectorized over a trailing (frequency) axis (Fortran).  Drop-in for
    scipy RegularGridInterpolator(('linear'), bounds_error=False, fill_value=fill)
    on a uniform grid, but without the per-call object and temporary overhead.
    @param field: [ng,ng] or [ng,ng,nw]; real input returns a real result
    @param    xg: the uniform axis shared by x and y (its first, last and spacing
                  are used; the exact endpoints keep boundary points in bounds)
    @param px,py: sample coordinates, any common shape
    @param  fill: value used outside the grid
    @return px.shape (+ (nw,)) array of interpolated values
    """
    field = np.asarray(field)
    real_in = not np.iscomplexobj(field)
    fld = np.ascontiguousarray(field, dtype=np.complex128)
    ng = fld.shape[0]
    if fld.shape[1] != ng:
        raise ValueError('bilinear_eval needs a square grid')
    nw = fld.shape[2] if fld.ndim == 3 else 1
    px = np.ascontiguousarray(px, dtype=np.float64)
    py = np.ascontiguousarray(py, dtype=np.float64)
    shape = px.shape
    npt = px.size
    out = np.empty((npt, nw), dtype=np.complex128)
    xg = np.asarray(xg, dtype=np.float64)
    _lib.bilinear_eval(out, fld, dbl(xg[0]), dbl(xg[1] - xg[0]), dbl(xg[-1]), i64(ng),
                       px.ravel(), py.ravel(), i64(npt), i64(nw), _cplx(fill))
    out = out.reshape(shape + ((nw,) if field.ndim == 3 else ()))
    return out.real.copy() if real_in else out


def bilinear_cell_eval(field: np.ndarray, px: np.ndarray, py: np.ndarray,
                       Dlx: float, Dly: float, zeta: float) -> np.ndarray:
    """
    @fn bilinear_cell_eval
    @brief Periodic bilinear interpolation of a cell-centred field on the (oblique)
    magnetic unit cell at scattered Cartesian points (Fortran).  Grid point k sits
    at fractional coordinate (k+0.5)/ng-0.5; the Cartesian -> fractional map is the
    inverse of x=(r1(1-zeta)+r2 zeta)Dlx, y=(-r1+r2)Dly.
    @param field: [ng,ng] real or complex cell field
    @param px,py: Cartesian sample coordinates, any common shape
    @param Dlx,Dly,zeta: cell geometry
    @return px.shape array of interpolated values (real in, real out)
    """
    field = np.asarray(field)
    real_in = not np.iscomplexobj(field)
    fld = np.ascontiguousarray(field, dtype=np.complex128)
    ng = fld.shape[0]
    px = np.ascontiguousarray(px, dtype=np.float64)
    py = np.ascontiguousarray(py, dtype=np.float64)
    out = np.empty(px.size, dtype=np.complex128)
    _lib.bilinear_cell_eval(out, fld, i64(ng), px.ravel(), py.ravel(), i64(px.size),
                            dbl(Dlx), dbl(Dly), dbl(zeta))
    out = out.reshape(px.shape)
    return out.real.copy() if real_in else out


def abrikosov_z(px: np.ndarray, py: np.ndarray, Dlx: float, Dly: float, zeta: float,
                nsum: int = 6) -> np.ndarray:
    """
    @fn abrikosov_z
    @brief Abrikosov (lowest-Landau-level) quasi-periodic order parameter of a
    magnetic unit cell, at scattered points (Fortran).  Its phase carries the
    vortex winding and its modulus vanishes at every core.
    @param px,py: Cartesian sample coordinates, any common shape
    @param Dlx,Dly,zeta: cell widths and obliqueness
    @param nsum: lattice-sum truncation (p = -nsum..nsum)
    @return px.shape complex array
    """
    px = np.ascontiguousarray(px, dtype=np.float64)
    py = np.ascontiguousarray(py, dtype=np.float64)
    out = np.empty(px.size, dtype=np.complex128)
    _lib.abrikosov_z(out, px.ravel(), py.ravel(), i64(px.size),
                     dbl(Dlx), dbl(Dly), dbl(zeta), i64(nsum))
    return out.reshape(px.shape)
