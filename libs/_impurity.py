from ctypes import *
import numpy as np
from ._loader import _lib

def gen_imp_ham(rvec: np.ndarray, ham_r: np.ndarray, ham_i: np.ndarray,
                rlist: np.ndarray, imp_list: np.ndarray, eps: float = 1.0e-5) -> np.ndarray:
    """
    @fn gen_imp_ham
    @brief Construct the real-space impurity Hamiltonian from bulk hopping blocks and impurity positions.
    @param     rvec: Real-space displacement vectors [Nr, 3] float64
    @param    ham_r: Real-part hopping blocks [Nr, Norb, Norb] complex128
    @param    ham_i: Imaginary-part hopping blocks [Nr, Norb, Norb] complex128
    @param    rlist: Impurity site positions [Nsite] float64
    @param imp_list: Impurity orbital indices [Nimp] int64
    @param      eps: Small tolerance for assembly
    @return ham_imp: Impurity Hamiltonian [Norb*Nsite, Norb*Nsite] complex128
    """
    Nr, Nimp, Nsite = len(rvec), len(imp_list), len(rlist)
    Norb = int(np.sqrt(ham_r.size / Nr))
    ham_imp = np.zeros((Norb * Nsite, Norb * Nsite), dtype=np.complex128)
    _lib.gen_imp_ham.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        POINTER(c_double),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.gen_imp_ham.restype = c_void_p
    _lib.gen_imp_ham(ham_imp, ham_r, rvec, ham_i, imp_list, rlist, byref(c_double(eps)),
                     byref(c_int64(Nimp)), byref(c_int64(Nsite)),
                     byref(c_int64(Nr)), byref(c_int64(Norb)))
    return ham_imp

def dft_imp_ham(ham_imp: np.ndarray, klist: np.ndarray, rlist: np.ndarray) -> np.ndarray:
    """
    @fn dft_imp_ham
    @brief Fourier-transform the real-space impurity Hamiltonian to k-space.
    @param  ham_imp: Real-space impurity Hamiltonian [Norb*Nsite, Norb*Nsite] complex128
    @param    klist: k-point list [Nk, 3] float64
    @param    rlist: Impurity site positions [Nsite] float64
    @return   ham_k: k-space impurity Hamiltonian [Norb*Nk, Norb*Nk] complex128
    """
    Nk, Nsite = len(klist), len(rlist)
    Norb = int(len(ham_imp) / Nsite)
    ham_k = np.zeros((Norb * Nk, Norb * Nk), dtype=np.complex128)
    _lib.get_dft_imp_ham.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64)
    ]
    _lib.get_dft_imp_ham.restype = c_void_p
    _lib.get_dft_imp_ham(ham_k, ham_imp, klist, rlist, byref(c_int64(Nk)),
                         byref(c_int64(Nsite)), byref(c_int64(Norb)))
    return ham_k

def get_imp_spectrum(uni: np.ndarray, eigs: np.ndarray, mu: float, wlist: np.ndarray,
                     klist: np.ndarray, rlist: np.ndarray, eta: float = 1.0e-3) -> np.ndarray:
    """
    @fn get_imp_spectrum
    @brief Compute the k- and frequency-resolved impurity spectral function A(k, w).
    @param     uni: Eigenvectors of the impurity Hamiltonian [Norb*Nsite, Norb*Nsite] complex128
    @param    eigs: Eigenvalues of the impurity Hamiltonian [Norb*Nsite] float64
    @param      mu: Chemical potential in eV
    @param   wlist: Real-frequency mesh [Nw] float64
    @param   klist: k-point list [Nk, 3] float64
    @param   rlist: Impurity site positions [Nsite] float64
    @param     eta: Lorentzian broadening parameter in eV
    @return spectrum: Impurity spectral function [Nk, Nw] complex128
    """
    Nw, Nk, Nsite = len(wlist), len(klist), len(rlist)
    Norb = int(len(eigs) / Nsite)
    spectrum = np.zeros((Nk, Nw), dtype=np.complex128)
    _lib.get_spectrum_spagehtti.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_double), POINTER(c_double)
    ]
    _lib.get_spectrum_spagehtti.restype = c_void_p
    _lib.get_spectrum_spagehtti(spectrum, uni, eigs, klist, rlist, wlist,
                                byref(c_int64(Nw)), byref(c_int64(Nk)),
                                byref(c_int64(Nsite)), byref(c_int64(Norb)),
                                byref(c_double(mu)), byref(c_double(eta)))
    return spectrum

def solve_cpa(hamk: np.ndarray, VA: np.ndarray, VB: np.ndarray,
              x: float, zlist: np.ndarray, pp: float = 0.5,
              maxiter: int = 500, tol: float = 1.0e-10,
              sigma_init: np.ndarray | None = None) -> np.ndarray:
    """
    @fn solve_cpa
    @brief Run CPA self-consistent loop over an array of complex frequencies.
    @param   hamk: k-space Hamiltonian [Nk, Norb, Norb] complex128
    @param     VA: Onsite potential of species A [Norb, Norb] complex128
    @param     VB: Onsite potential of species B [Norb, Norb] complex128
    @param      x: Concentration of species A (0 < x < 1)
    @param  zlist: Complex frequency array [Nw] complex128
    @param     pp: Linear mixing parameter (default 0.5)
    @param maxiter: Maximum CPA iterations per frequency (default 500)
    @param    tol: Convergence tolerance (default 1e-10)
    @param sigma_init: Initial guess [Nw, Norb, Norb] complex128 (default: VCA = x*VA + (1-x)*VB)
    @return sigma_cpa: Converged CPA self-energy [Nw, Norb, Norb] complex128
    """
    Nk = len(hamk)
    Norb = hamk.shape[1]
    Nw = len(zlist)
    VA_c = np.ascontiguousarray(VA, dtype=np.complex128)
    VB_c = np.ascontiguousarray(VB, dtype=np.complex128)
    zlist_c = np.ascontiguousarray(zlist, dtype=np.complex128)

    if sigma_init is not None:
        sigma_cpa = np.ascontiguousarray(sigma_init, dtype=np.complex128)
    else:
        sigma_cpa = np.zeros((Nw, Norb, Norb), dtype=np.complex128)
        vca = x * VA_c + (1.0 - x) * VB_c
        for iw in range(Nw):
            sigma_cpa[iw, :, :] = vca

    _lib.solve_cpa_array.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        POINTER(c_double),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        POINTER(c_double),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_double)
    ]
    _lib.solve_cpa_array.restype = None
    _lib.solve_cpa_array(sigma_cpa, hamk, VA_c, VB_c,
                         byref(c_double(x)), zlist_c,
                         byref(c_double(pp)),
                         byref(c_int64(Nk)), byref(c_int64(Norb)),
                         byref(c_int64(Nw)), byref(c_int64(maxiter)),
                         byref(c_double(tol)))
    return sigma_cpa
