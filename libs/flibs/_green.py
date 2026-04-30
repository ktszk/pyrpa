from ctypes import *
import numpy as np
from ._loader import _lib

def gen_Green0(eig: np.ndarray, uni: np.ndarray,
               mu: float, temp: float, Nw: int) -> np.ndarray:
    """
    @fn gen_Green0
    @brief Construct the non-interacting Matsubara Green's function G0(k, iω_n).
    @param   eig: Eigenvalues [Nk, Norb] float64
    @param   uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param    mu: Chemical potential in eV
    @param  temp: Temperature in eV
    @param    Nw: Number of Matsubara frequency points
    @return   Gk: Non-interacting Green's function [Norb, Norb, Nw, Nk] complex128
    """
    Nk = len(eig)
    Norb = int(eig.size / Nk)
    Gk = np.zeros((Norb, Norb, Nw, Nk), dtype=np.complex128)
    _lib.gen_green0_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        POINTER(c_double), POINTER(c_double),
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.gen_green0_.restype = c_void_p
    _lib.gen_green0_(Gk, eig, uni, byref(c_double(mu)), byref(c_double(temp)),
                     byref(c_int64(Nk)), byref(c_int64(Nw)), byref(c_int64(Norb)))
    return Gk

def gen_green(selfen: np.ndarray, hamk: np.ndarray,
              mu: float, temp: float) -> np.ndarray:
    """
    @fn gen_green
    @brief Construct the interacting Matsubara Green's function G(k, iω_n) = [iω_n + mu - H(k) - Sigma(k,iω_n)]^-1.
    @param  selfen: Self-energy tensor [Norb, Norb, Nw, Nk] complex128
    @param    hamk: k-space Hamiltonian [Nk, Norb, Norb] complex128
    @param      mu: Chemical potential in eV
    @param    temp: Temperature in eV
    @return     Gk: Interacting Green's function [Norb, Norb, Nw, Nk] complex128
    """
    Nk = len(hamk)
    Norb = int(np.sqrt(hamk.size / Nk))
    Nw = int(selfen.size / (Nk * Norb * Norb))
    Gk = np.zeros((Norb, Norb, Nw, Nk), dtype=np.complex128)
    _lib.gen_green_inv_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        POINTER(c_double), POINTER(c_double),
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.gen_green_inv_.restype = c_void_p
    _lib.gen_green_inv_(Gk, selfen, hamk, byref(c_double(mu)), byref(c_double(temp)),
                        byref(c_int64(Nk)), byref(c_int64(Nw)), byref(c_int64(Norb)))
    _lib.getinv_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.getinv_.restype = c_void_p
    _lib.getinv_(Gk, byref(c_int64(Nk)), byref(c_int64(Nw)), byref(c_int64(Norb)))
    return Gk

def gen_green_from_eig(selfen: np.ndarray, eig: np.ndarray,
                       uni: np.ndarray, mu: float, temp: float) -> np.ndarray:
    """
    @fn gen_green_from_eig
    @brief Construct the interacting Green's function using pre-computed eigenvalues and eigenvectors.
    @param  selfen: Self-energy tensor [Norb, Norb, Nw, Nk] complex128
    @param     eig: Eigenvalues [Nk, Norb] float64
    @param     uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param      mu: Chemical potential in eV
    @param    temp: Temperature in eV
    @return     Gk: Interacting Green's function [Norb, Norb, Nw, Nk] complex128
    """
    Nk = len(eig)
    Norb = int(eig.size / Nk)
    Nw = int(selfen.size / (Nk * Norb * Norb))
    Gk = np.zeros((Norb, Norb, Nw, Nk), dtype=np.complex128)
    _lib.gen_green_inv_from_eig.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        POINTER(c_double), POINTER(c_double),
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.gen_green_inv_from_eig.restype = c_void_p
    _lib.gen_green_inv_from_eig(Gk, selfen, uni, eig, byref(c_double(mu)), byref(c_double(temp)),
                                byref(c_int64(Nk)), byref(c_int64(Nw)), byref(c_int64(Norb)))
    _lib.getinv.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.getinv.restype = c_void_p
    _lib.getinv(Gk, byref(c_int64(Nk)), byref(c_int64(Nw)), byref(c_int64(Norb)))
    return Gk

def gen_tr_Greenw_0(eig: np.ndarray, mu: float,
                    wlist: np.ndarray, delta: float) -> np.ndarray:
    """
    @fn gen_tr_Greenw_0
    @brief Compute the k-resolved spectral function (imaginary part of trace of G0) on the real-frequency axis.
    @param    eig: Eigenvalues [Nk, Norb] float64
    @param     mu: Chemical potential in eV
    @param  wlist: Real-frequency mesh [Nw] float64
    @param  delta: Lorentzian broadening parameter in eV
    @return  trGk: -Im[Tr G0(k,w+i*delta)] [Nk, Nw] float64
    """
    Nk, Nw = len(eig), len(wlist)
    Norb = int(eig.size / Nk)
    trGk = np.zeros((Nk, Nw), dtype=np.complex128)
    _lib.gen_tr_greenw_0.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        POINTER(c_double), POINTER(c_double),
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.gen_tr_greenw_0.restype = c_void_p
    _lib.gen_tr_greenw_0(trGk, wlist, eig, byref(c_double(mu)), byref(c_double(delta)),
                         byref(c_int64(Nk)), byref(c_int64(Nw)), byref(c_int64(Norb)))
    return -trGk.imag

def gen_dos(eig: np.ndarray, uni: np.ndarray, mu: float,
            wlist: np.ndarray, delta: float) -> np.ndarray:
    """
    @fn gen_dos
    @brief Compute the orbital-resolved partial density of states (pDOS) on the real-frequency axis.
    @param    eig: Eigenvalues [Nk, Norb] float64
    @param    uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param     mu: Chemical potential in eV
    @param  wlist: Real-frequency mesh [Nw] float64
    @param  delta: Lorentzian broadening parameter in eV
    @return  pDos: Orbital-resolved pDOS -Im[G(orb,w)] [Norb, Nw] float64
    """
    Nk, Nw = len(eig), len(wlist)
    Norb = int(eig.size / Nk)
    pDos = np.zeros((Norb, Nw), dtype=np.complex128)
    _lib.gen_dos.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        POINTER(c_double), POINTER(c_double),
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.gen_dos.restype = c_void_p
    _lib.gen_dos(pDos, wlist, eig, uni, byref(c_double(mu)), byref(c_double(delta)),
                 byref(c_int64(Nk)), byref(c_int64(Nw)), byref(c_int64(Norb)))
    return -pDos.imag

def get_a(inp_data: np.ndarray, xlist: np.ndarray) -> np.ndarray:
    """
    @fn get_a
    @brief Compute Pade continued-fraction coefficients from input data on Matsubara grid.
    @param inp_data: Input data (e.g. Green's function) at Matsubara points [Np] complex128
    @param    xlist: Matsubara frequency points [Np] complex128
    @return       a: Pade coefficient array [Np] complex128
    """
    Np = len(inp_data)
    a = np.zeros(Np, dtype=np.complex128)
    _lib.get_a_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        POINTER(c_int64)
    ]
    _lib.get_a_.restype = None
    _lib.get_a_(a, xlist, inp_data, byref(c_int64(Np)))
    return a

def get_QP(a: np.ndarray, xlist: np.ndarray, wlist: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    @fn get_QP
    @brief Evaluate Pade numerator Q and denominator P on the real-frequency axis from continued-fraction coefficients.
    @param      a: Pade coefficient array [Np] complex128 (from get_a)
    @param  xlist: Matsubara frequency points [Np] complex128
    @param  wlist: Real-frequency mesh [Nw] complex128
    @retval     Q: Pade numerator evaluated on real axis [Nw] complex128
    @retval     P: Pade denominator evaluated on real axis [Nw] complex128
    """
    Nw, Np = len(wlist), len(a)
    Q = np.zeros(Nw, dtype=np.complex128)
    P = np.zeros(Nw, dtype=np.complex128)
    _lib.get_qp_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.get_qp_.restype = None
    _lib.get_qp_(P, Q, a, xlist, wlist, byref(c_int64(Nw)), byref(c_int64(Np)))
    return Q, P

def pade_with_trace(A: np.ndarray, iwlist: np.ndarray, wlist: np.ndarray) -> np.ndarray:
    """
    @fn pade_with_trace
    @brief Perform Pade analytic continuation and trace over orbital indices to obtain the real-axis spectrum.
    @param      A: Matrix on Matsubara axis [Norb, Norb, Niw, Nk] complex128
    @param iwlist: Matsubara frequency points [Niw] complex128
    @param  wlist: Real-frequency mesh [Nw] complex128
    @return     B: Analytically-continued trace [Nk, Nw] complex128
    """
    Nk, Nw, Niw = len(A.T), len(wlist), len(iwlist)
    Norb = int(np.sqrt(A.size / (Nk * Niw)))
    B = np.zeros((Nk, Nw), dtype=np.complex128)
    _lib.pade_with_trace.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.pade_with_trace.restype = None
    _lib.pade_with_trace(A, B, iwlist, wlist, byref(c_int64(Nk)), byref(c_int64(Niw)),
                         byref(c_int64(Nw)), byref(c_int64(Norb)))
    return B
