import ctypes
import numpy as np
import libs.flibs._loader as _loader

# Public API of this submodule; keeps `from ._x import *` in __init__.py from
# re-exporting numpy, scipy, ctypes and the sibling packages imported above.
__all__ = [
    'gen_imp_ham', 'dft_imp_ham', 'get_imp_spectrum', 'solve_cpa'
]

# --- ctypes signatures: set once at import.
# All Fortran entry points are subroutines, so restype is always None.
_loader._lib.gen_imp_ham.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),  # ham_imp (out)
    np.ctypeslib.ndpointer(dtype=np.complex128),  # ham_r
    np.ctypeslib.ndpointer(dtype=np.float64),     # rvec
    np.ctypeslib.ndpointer(dtype=np.complex128),  # ham_i
    np.ctypeslib.ndpointer(dtype=np.complex128),  # ham_ri
    np.ctypeslib.ndpointer(dtype=np.int64),       # imp_list
    np.ctypeslib.ndpointer(dtype=np.float64),     # rlist
    ctypes.POINTER(ctypes.c_double),                            # eps
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nimp, Nsite
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)            # Nr, Norb
]
_loader._lib.gen_imp_ham.restype = None
_loader._lib.get_dft_imp_ham.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),  # ham_k (out)
    np.ctypeslib.ndpointer(dtype=np.complex128),  # ham_imp
    np.ctypeslib.ndpointer(dtype=np.float64),     # klist
    np.ctypeslib.ndpointer(dtype=np.float64),     # rlist
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nk, Nsite
    ctypes.POINTER(ctypes.c_int64)                              # Norb
]
_loader._lib.get_dft_imp_ham.restype = None
_loader._lib.get_spectrum_spaghetti.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),  # spa (out)
    np.ctypeslib.ndpointer(dtype=np.complex128),  # uni
    np.ctypeslib.ndpointer(dtype=np.float64),     # eigs
    np.ctypeslib.ndpointer(dtype=np.float64),     # klist
    np.ctypeslib.ndpointer(dtype=np.float64),     # rlist
    np.ctypeslib.ndpointer(dtype=np.float64),     # wlist
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nw, Nk
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nsite, Norb
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)          # mu, eta
]
_loader._lib.get_spectrum_spaghetti.restype = None
_loader._lib.solve_cpa_array.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),  # sigma_cpa_w (inout)
    np.ctypeslib.ndpointer(dtype=np.complex128),  # hamk
    np.ctypeslib.ndpointer(dtype=np.complex128),  # VA
    np.ctypeslib.ndpointer(dtype=np.complex128),  # VB
    ctypes.POINTER(ctypes.c_double),                            # x
    np.ctypeslib.ndpointer(dtype=np.complex128),  # zlist
    ctypes.POINTER(ctypes.c_double),                            # pp
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nk, Norb
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nw, maxiter
    ctypes.POINTER(ctypes.c_double)                             # tol
]
_loader._lib.solve_cpa_array.restype = None


def gen_imp_ham(rvec: np.ndarray, ham_r: np.ndarray, ham_i: np.ndarray, ham_ri: np.ndarray,
                rlist: np.ndarray, imp_list: np.ndarray, eps: float = 1.0e-5) -> np.ndarray:
    """
    @brief Construct the real-space supercell Hamiltonian for a host crystal with embedded impurities.
    @param     rvec: Real-space displacement vectors [Nr, 3] float64
    @param    ham_r: Host-host hopping blocks H_rr(R) [Nr, Norb, Norb] complex128
    @param    ham_i: Impurity-impurity hopping blocks H_ii(R) [Nr, Norb, Norb] complex128
    @param   ham_ri: Cross-species hopping blocks H_ri(R) [Nr, Norb, Norb] complex128
    @param    rlist: Site positions in fractional coords [Nsite, 3] float64
    @param imp_list: 0-based impurity site indices [Nimp] int64
    @param      eps: Small tolerance for R-vector matching
    @return ham_imp: Impurity Hamiltonian [Norb*Nsite, Norb*Nsite] complex128
    """
    Nr, Nimp, Nsite = len(rvec), len(imp_list), len(rlist)
    Norb = ham_r.shape[1]
    ham_imp = np.zeros((Norb * Nsite, Norb * Nsite), dtype=np.complex128)
    _loader._lib.gen_imp_ham(ham_imp, ham_r, rvec, ham_i, ham_ri, imp_list, rlist, _loader.dbl(eps),
                             _loader.i64(Nimp), _loader.i64(Nsite),
                             _loader.i64(Nr), _loader.i64(Norb))
    return ham_imp


def dft_imp_ham(ham_imp: np.ndarray, klist: np.ndarray, rlist: np.ndarray) -> np.ndarray:
    """
    @brief Fourier-transform the real-space impurity Hamiltonian to k-space.
    @param  ham_imp: Real-space impurity Hamiltonian [Norb*Nsite, Norb*Nsite] complex128
    @param    klist: k-point list [Nk, 3] float64
    @param    rlist: Impurity site positions [Nsite] float64
    @return   ham_k: k-space impurity Hamiltonian [Norb*Nk, Norb*Nk] complex128
    """
    Nk, Nsite = len(klist), len(rlist)
    Norb = int(len(ham_imp) / Nsite)
    ham_k = np.zeros((Norb * Nk, Norb * Nk), dtype=np.complex128)
    _loader._lib.get_dft_imp_ham(ham_k, ham_imp, klist, rlist, _loader.i64(Nk),
                                 _loader.i64(Nsite), _loader.i64(Norb))
    return ham_k


def get_imp_spectrum(uni: np.ndarray, eigs: np.ndarray, mu: float, wlist: np.ndarray,
                     klist: np.ndarray, rlist: np.ndarray, eta: float = 1.0e-3) -> np.ndarray:
    """
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
    _loader._lib.get_spectrum_spaghetti(spectrum, uni, eigs, klist, rlist, wlist,
                                        _loader.i64(Nw), _loader.i64(Nk),
                                        _loader.i64(Nsite), _loader.i64(Norb),
                                        _loader.dbl(mu), _loader.dbl(eta))
    return spectrum


def solve_cpa(hamk: np.ndarray, VA: np.ndarray, VB: np.ndarray,
              x: float, zlist: np.ndarray, pp: float = 0.5,
              maxiter: int = 500, tol: float = 1.0e-10,
              sigma_init: np.ndarray | None = None) -> np.ndarray:
    """
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

    _loader._lib.solve_cpa_array(sigma_cpa, hamk, VA_c, VB_c,
                                 _loader.dbl(x), zlist_c,
                                 _loader.dbl(pp),
                                 _loader.i64(Nk), _loader.i64(Norb),
                                 _loader.i64(Nw), _loader.i64(maxiter),
                                 _loader.dbl(tol))
    return sigma_cpa
