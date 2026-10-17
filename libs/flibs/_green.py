import ctypes
import numpy as np
import libs.flibs._loader as _loader

# Public API of this submodule; keeps `from ._x import *` in __init__.py from
# re-exporting numpy, scipy, ctypes and the sibling packages imported above.
__all__ = [
    'gen_Green0', 'gen_green', 'gen_green_from_eig', 'gen_tr_Greenw_0', 'gen_dos', 'get_a', 'get_QP',
    'pade_with_trace'
]

# --- ctypes signatures: set once at import.
# All Fortran entry points are subroutines, so restype is always None.
_loader._lib.gen_green0_.argtypes = [
    _loader.c_array(np.complex128),          # Gk (out)
    _loader.c_array(np.float64),             # eig
    _loader.c_array(np.complex128),          # uni
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),                 # mu, temp
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)  # Nk, Nw, Norb
]
_loader._lib.gen_green0_.restype = None
_loader._lib.gen_green_inv_.argtypes = [
    _loader.c_array(np.complex128),          # Gk (out)
    _loader.c_array(np.complex128),          # self
    _loader.c_array(np.complex128),          # hamk
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),                 # mu, temp
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)  # Nk, Nw, Norb
]
_loader._lib.gen_green_inv_.restype = None
_loader._lib.getinv_.argtypes = [
    _loader.c_array(np.complex128),          # Gk (inout)
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)  # Nk, Nw, Norb
]
_loader._lib.getinv_.restype = None
_loader._lib.gen_green_inv_from_eig.argtypes = [
    _loader.c_array(np.complex128),          # Gk (out)
    _loader.c_array(np.complex128),          # self
    _loader.c_array(np.complex128),          # uni
    _loader.c_array(np.float64),             # eig
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),                 # mu, temp
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)  # Nk, Nw, Norb
]
_loader._lib.gen_green_inv_from_eig.restype = None
_loader._lib.gen_tr_greenw_0.argtypes = [
    _loader.c_array(np.complex128),          # trGk (out)
    _loader.c_array(np.float64),             # wl
    _loader.c_array(np.float64),             # eig
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),                 # mu, delta
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)  # Nk, Nw, Norb
]
_loader._lib.gen_tr_greenw_0.restype = None
_loader._lib.gen_dos.argtypes = [
    _loader.c_array(np.complex128),          # Dos (out)
    _loader.c_array(np.float64),             # wl
    _loader.c_array(np.float64),             # eig
    _loader.c_array(np.complex128),          # uni
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),                 # mu, delta
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)  # Nk, Nw, Norb
]
_loader._lib.gen_dos.restype = None
_loader._lib.get_a_.argtypes = [
    _loader.c_array(np.complex128),  # a (out)
    _loader.c_array(np.complex128),  # xn
    _loader.c_array(np.complex128),  # inp_data
    ctypes.POINTER(ctypes.c_int64)                              # Np
]
_loader._lib.get_a_.restype = None
_loader._lib.get_qp_.argtypes = [
    _loader.c_array(np.complex128),  # P (out)
    _loader.c_array(np.complex128),  # Q (out)
    _loader.c_array(np.complex128),  # a
    _loader.c_array(np.complex128),  # xn
    _loader.c_array(np.complex128),  # wlist
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)            # Nw, Np
]
_loader._lib.get_qp_.restype = None
_loader._lib.pade_with_trace.argtypes = [
    _loader.c_array(np.complex128),  # arrayin
    _loader.c_array(np.complex128),  # arrayout (out)
    _loader.c_array(np.complex128),  # iwlist
    _loader.c_array(np.complex128),  # wlist
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nk, Niw
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)            # Nw, Norb
]
_loader._lib.pade_with_trace.restype = None


def gen_Green0(eig: np.ndarray, uni: np.ndarray,
               mu: float, temp: float, Nw: int) -> np.ndarray:
    """
    @brief Construct the non-interacting Matsubara Green's function G0(k, iω_n).
    @param   eig: Eigenvalues [Nk, Norb] float64
    @param   uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param    mu: Chemical potential in eV
    @param  temp: Temperature in eV
    @param    Nw: Number of Matsubara frequency points
    @return   Gk: Non-interacting Green's function [Norb, Norb, Nw, Nk] complex128
    """
    Nk = len(eig)
    Norb = eig.shape[1]
    Gk = np.zeros((Norb, Norb, Nw, Nk), dtype=np.complex128)
    _loader._lib.gen_green0_(Gk, eig, uni, _loader.dbl(mu), _loader.dbl(temp),
                             _loader.i64(Nk), _loader.i64(Nw), _loader.i64(Norb))
    return Gk


def gen_green(selfen: np.ndarray, hamk: np.ndarray,
              mu: float, temp: float) -> np.ndarray:
    """
    @brief Construct the interacting Matsubara Green's function G(k, iω_n) = [iω_n + mu - H(k) - Sigma(k,iω_n)]^-1.
    @param  selfen: Self-energy tensor [Norb, Norb, Nw, Nk] complex128
    @param    hamk: k-space Hamiltonian [Nk, Norb, Norb] complex128
    @param      mu: Chemical potential in eV
    @param    temp: Temperature in eV
    @return     Gk: Interacting Green's function [Norb, Norb, Nw, Nk] complex128
    """
    Nk = len(hamk)
    Norb = hamk.shape[1]
    Nw = selfen.shape[2]
    Gk = np.zeros((Norb, Norb, Nw, Nk), dtype=np.complex128)
    _loader._lib.gen_green_inv_(Gk, selfen, hamk, _loader.dbl(mu), _loader.dbl(temp),
                                _loader.i64(Nk), _loader.i64(Nw), _loader.i64(Norb))
    _loader._lib.getinv_(Gk, _loader.i64(Nk), _loader.i64(Nw), _loader.i64(Norb))
    return Gk


def gen_green_from_eig(selfen: np.ndarray, eig: np.ndarray,
                       uni: np.ndarray, mu: float, temp: float) -> np.ndarray:
    """
    @brief Construct the interacting Green's function using pre-computed eigenvalues and eigenvectors.
    @param  selfen: Self-energy tensor [Norb, Norb, Nw, Nk] complex128
    @param     eig: Eigenvalues [Nk, Norb] float64
    @param     uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param      mu: Chemical potential in eV
    @param    temp: Temperature in eV
    @return     Gk: Interacting Green's function [Norb, Norb, Nw, Nk] complex128
    """
    Nk = len(eig)
    Norb = eig.shape[1]
    Nw = selfen.shape[2]
    Gk = np.zeros((Norb, Norb, Nw, Nk), dtype=np.complex128)
    _loader._lib.gen_green_inv_from_eig(Gk, selfen, uni, eig, _loader.dbl(mu), _loader.dbl(temp),
                                        _loader.i64(Nk), _loader.i64(Nw), _loader.i64(Norb))
    _loader._lib.getinv_(Gk, _loader.i64(Nk), _loader.i64(Nw), _loader.i64(Norb))
    return Gk


def gen_tr_Greenw_0(eig: np.ndarray, mu: float,
                    wlist: np.ndarray, delta: float) -> np.ndarray:
    """
    @brief Compute the k-resolved spectral function (imaginary part of trace of G0) on the real-frequency axis.
    @param    eig: Eigenvalues [Nk, Norb] float64
    @param     mu: Chemical potential in eV
    @param  wlist: Real-frequency mesh [Nw] float64
    @param  delta: Lorentzian broadening parameter in eV
    @return  trGk: -Im[Tr G0(k,w+i*delta)] [Nk, Nw] float64
    """
    Nk, Nw = len(eig), len(wlist)
    Norb = eig.shape[1]
    trGk = np.zeros((Nk, Nw), dtype=np.complex128)
    _loader._lib.gen_tr_greenw_0(trGk, wlist, eig, _loader.dbl(mu), _loader.dbl(delta),
                                 _loader.i64(Nk), _loader.i64(Nw), _loader.i64(Norb))
    return -trGk.imag


def gen_dos(eig: np.ndarray, uni: np.ndarray, mu: float,
            wlist: np.ndarray, delta: float) -> np.ndarray:
    """
    @brief Compute the orbital-resolved partial density of states (pDOS) on the real-frequency axis.
    @param    eig: Eigenvalues [Nk, Norb] float64
    @param    uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param     mu: Chemical potential in eV
    @param  wlist: Real-frequency mesh [Nw] float64
    @param  delta: Lorentzian broadening parameter in eV
    @return  pDos: Orbital-resolved pDOS -Im[G(orb,w)] [Norb, Nw] float64
    """
    Nk, Nw = len(eig), len(wlist)
    Norb = eig.shape[1]
    pDos = np.zeros((Norb, Nw), dtype=np.complex128)
    _loader._lib.gen_dos(pDos, wlist, eig, uni, _loader.dbl(mu), _loader.dbl(delta),
                         _loader.i64(Nk), _loader.i64(Nw), _loader.i64(Norb))
    return -pDos.imag


def get_a(inp_data: np.ndarray, xlist: np.ndarray) -> np.ndarray:
    """
    @brief Compute Pade continued-fraction coefficients from input data on Matsubara grid.
    @param inp_data: Input data (e.g. Green's function) at Matsubara points [Np] complex128
    @param    xlist: Matsubara frequency points [Np] complex128
    @return       a: Pade coefficient array [Np] complex128
    """
    Np = len(inp_data)
    a = np.zeros(Np, dtype=np.complex128)
    _loader._lib.get_a_(a, xlist, inp_data, _loader.i64(Np))
    return a


def get_QP(a: np.ndarray, xlist: np.ndarray, wlist: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    @brief Evaluate Pade numerator Q and denominator P on the real-frequency axis from continued-fraction coefficients.
    @param      a: Pade coefficient array [Np] complex128 (from get_a)
    @param  xlist: Matsubara frequency points [Np] complex128
    @param  wlist: Real-frequency mesh [Nw] complex128
    @return     Q: Pade numerator evaluated on real axis [Nw] complex128
    @return     P: Pade denominator evaluated on real axis [Nw] complex128
    """
    Nw, Np = len(wlist), len(a)
    Q = np.zeros(Nw, dtype=np.complex128)
    P = np.zeros(Nw, dtype=np.complex128)
    _loader._lib.get_qp_(P, Q, a, xlist, wlist, _loader.i64(Nw), _loader.i64(Np))
    return Q, P


def pade_with_trace(A: np.ndarray, iwlist: np.ndarray, wlist: np.ndarray) -> np.ndarray:
    """
    @brief Perform Pade analytic continuation and trace over orbital indices to obtain the real-axis spectrum.
    @param      A: Matrix on Matsubara axis [Norb, Norb, Niw, Nk] complex128
    @param iwlist: Matsubara frequency points [Niw] complex128
    @param  wlist: Real-frequency mesh [Nw] complex128
    @return     B: Analytically-continued trace [Nk, Nw] complex128
    """
    # Callers typically pass a Matsubara slice G[:, :, :Niw, :], which is NOT contiguous;
    # ndpointer would hand Fortran the raw strided buffer, so copy to C order first.
    A = np.ascontiguousarray(A, dtype=np.complex128)
    iwlist = np.ascontiguousarray(iwlist, dtype=np.complex128)
    wlist = np.ascontiguousarray(wlist, dtype=np.complex128)
    Nk, Nw, Niw = len(A.T), len(wlist), len(iwlist)
    Norb = A.shape[0]
    B = np.zeros((Nk, Nw), dtype=np.complex128)
    _loader._lib.pade_with_trace(A, B, iwlist, wlist, _loader.i64(Nk), _loader.i64(Niw),
                                 _loader.i64(Nw), _loader.i64(Norb))
    return B
