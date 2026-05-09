from ctypes import *
import numpy as np
from ._loader import _lib

def linearized_eliashberg(chi: np.ndarray, Gk: np.ndarray, uni: np.ndarray, init_delta: np.ndarray,
                          Smat: np.ndarray, Cmat: np.ndarray, olist: np.ndarray, plist: np.ndarray,
                          kmap: np.ndarray, invk: np.ndarray, Nx: int, Ny: int, Nz: int,
                          temp: float, gap_sym: int, eps: float = 1.0e-5, itemax: int = 300,
                          arnoldi_m: int = 10) -> np.ndarray:
    """
    @fn linearized_eliashberg
    @brief Solve the linearized Eliashberg gap equation (without SOC) to obtain the superconducting gap function.
    @param       chi: Irreducible susceptibility [Nchi, Nchi, Nw, Nk] complex128
    @param        Gk: Green's function on full k-grid [Norb, Norb, Nw, Nkall] complex128
    @param       uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param init_delta: Initial gap function [Norb, Norb, Nw, Nkall] float64
    @param      Smat: Spin interaction matrix [Nchi, Nchi] float64
    @param      Cmat: Charge interaction matrix [Nchi, Nchi] float64
    @param     olist: Orbital index list [Nchi] int64
    @param     plist: Parity list for gap symmetry [Norb] float64
    @param      kmap: Full-to-reduced k-point mapping [Nkall] int64
    @param      invk: Inverse k-point mapping [Nkall] int64
    @param    Nx,Ny,Nz: k-grid dimensions
    @param      temp: Temperature in eV
    @param   gap_sym: Gap symmetry index (0=s, 1=dx2-y2, etc.)
    @param       eps: Convergence tolerance
    @param    itemax: Maximum number of iterations
    @return    delta: Linearized gap function [Norb, Norb, Nw, Nkall] complex128
    """
    Norb, Nchi = len(Gk), len(Smat)
    Nkall, Nk, Nw = len(kmap), len(Gk[0, 0, 0]), len(Gk[0, 0])
    delta = np.zeros((Norb, Norb, Nw, Nk), dtype=np.complex128)
    _lib.lin_eliash.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        POINTER(c_double), POINTER(c_double),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64)
    ]
    _lib.lin_eliash.restype = None
    _lib.lin_eliash(delta, chi, Gk, uni, init_delta, Smat, Cmat, olist, plist, kmap, invk,
                    byref(c_double(temp)), byref(c_double(eps)), byref(c_int64(Nkall)),
                    byref(c_int64(Nk)), byref(c_int64(Nw)), byref(c_int64(Nchi)),
                    byref(c_int64(Norb)), byref(c_int64(Nx)), byref(c_int64(Ny)),
                    byref(c_int64(Nz)), byref(c_int64(itemax)), byref(c_int64(gap_sym)),
                    byref(c_int64(arnoldi_m)))
    return delta

def linearized_eliashberg_soc(chi: np.ndarray, Gk: np.ndarray, uni: np.ndarray, init_delta: np.ndarray,
                              Vmat: np.ndarray, sgnsig: np.ndarray, sgnsig2: np.ndarray, plist: np.ndarray,
                              slist: np.ndarray, olist: np.ndarray, kmap: np.ndarray, invk: np.ndarray,
                              invs: np.ndarray, invschi: np.ndarray, Nx: int, Ny: int, Nz: int, temp: float,
                              gap_sym: int, eps: float = 1.0e-4, itemax: int = 300,
                              arnoldi_m: int = 10) -> np.ndarray:
    """
    @fn linearized_eliashberg_soc
    @brief Solve the linearized Eliashberg gap equation with spin-orbit coupling (SOC).
    @param        chi: Irreducible susceptibility [Nchi, Nchi, Nw, Nk] complex128
    @param         Gk: Green's function on full k-grid [Norb, Norb, Nw, Nkall] complex128
    @param        uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param init_delta: Initial gap function [Norb, Norb, Nw, Nkall] float64
    @param       Vmat: SOC interaction matrix [Nchi, Nchi] float64
    @param     sgnsig: Sign matrix for spin-spin product [Norb, Norb] float64
    @param    sgnsig2: Sign matrix for chi orbital mapping [Nchi, Nchi] float64
    @param      plist: Parity list for gap symmetry [Norb] float64
    @param      slist: Spin label list [Norb] int64
    @param      olist: Orbital index list [Nchi] int64
    @param       kmap: Full-to-reduced k-point mapping [Nkall] int64
    @param       invk: Inverse k-point mapping [Nkall] int64
    @param       invs: Inverse spin index mapping [Norb] int64
    @param    invschi: Inverse chi orbital mapping [Nchi] int64
    @param   Nx,Ny,Nz: k-grid dimensions
    @param       temp: Temperature in eV
    @param    gap_sym: Gap symmetry index
    @param        eps: Convergence tolerance
    @param     itemax: Maximum number of iterations
    @return     delta: Linearized gap function [Norb, Norb, Nw, Nkall] complex128
    """
    Norb, Nchi = len(slist), len(Vmat)
    Nkall, Nk, Nw = len(kmap), len(Gk[0, 0, 0]), len(Gk[0, 0])
    delta = np.zeros((Norb, Norb, Nw, Nkall), dtype=np.complex128)
    _lib.lin_eliash_soc.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        POINTER(c_double), POINTER(c_double),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64)
    ]
    _lib.lin_eliash_soc.restype = None
    _lib.lin_eliash_soc(delta, chi, Gk, uni, init_delta, Vmat, sgnsig, sgnsig2, plist, olist,
                        slist, kmap, invk, invs, invschi, byref(c_double(temp)),
                        byref(c_double(eps)), byref(c_int64(Nkall)), byref(c_int64(Nk)),
                        byref(c_int64(Nw)), byref(c_int64(Nchi)), byref(c_int64(Norb)),
                        byref(c_int64(Nx)), byref(c_int64(Ny)), byref(c_int64(Nz)),
                        byref(c_int64(itemax)), byref(c_int64(gap_sym)),
                        byref(c_int64(arnoldi_m)))
    return delta

def nonlinear_eliashberg(delta_init: np.ndarray, Gk: np.ndarray, hamk: np.ndarray,
                         Smat: np.ndarray, Cmat: np.ndarray,
                         olist: np.ndarray, plist: np.ndarray, kmap: np.ndarray,
                         invk: np.ndarray, mu: float, temp: float, gap_sym: int,
                         Nx: int, Ny: int, Nz: int,
                         sw_sigma: bool = False, sw_Vconst: bool = False, eps: float = 1.0e-4,
                         itemax: int = 100, m_diis: int = 5) -> tuple[np.ndarray, np.ndarray]:
    """
    @fn nonlinear_eliashberg
    @brief Solve the non-linear FLEX-Eliashberg self-consistency loop (no SOC).
    Iteratively updates the gap Δ, normal/anomalous Green's functions G/F,
    irreducible χ_0_S = χ^GG + χ^FF, χ_0_C = χ^GG - χ^FF, and the FLEX
    vertices V_σ, V_Δ. Pulay/DIIS extrapolation is used for Δ mixing
    (linear mixing pp=0.3 falls back when m_diis<=1 or for the first iteration).
    Convergence: |Δ_new - Δ|_∞ / |Δ|_∞ < eps.

    @param  delta_init: Initial gap function [Norb, Norb, Nw, Nk] complex128.
                        Caller is responsible for providing a non-zero initial Δ.
                        Recommended sources:
                          - flibs.get_initial_delta(...) for a symmetry-adapted seed
                          - a converged linearized_eliashberg() output (scaled to a
                            physical magnitude, e.g. ~1.764·T) for a fast start
                        If the array is C-contiguous complex128 it is modified in place;
                        otherwise a contiguous copy is created.
    @param          Gk: Initial normal Green's function [Norb, Norb, Nw, Nk] complex128
                        (typically from flibs.gen_Green0 or a converged normal-state solution)
    @param        hamk: k-space Hamiltonian [Nk, Norb, Norb] complex128
    @param        Smat: Spin interaction matrix [Nchi, Nchi] float64
    @param        Cmat: Charge interaction matrix [Nchi, Nchi] float64
    @param       olist: Orbital index pairs for chi [Nchi, 2] int64
    @param       plist: Parity factors per orbital [Norb] float64
    @param        kmap: Full k-grid to irreducible-grid mapping [Nkall, 3] int64
    @param        invk: Inverse k-point list [Nkall, 3] int64
    @param          mu: Chemical potential in eV
    @param        temp: Temperature in eV
    @param     gap_sym: Gap symmetry (>=0 singlet, <0 triplet)
    @param    Nx,Ny,Nz: k-grid dimensions
    @param    sw_sigma: Include FLEX self-energy correction in the dressed G_0 (default False)
    @param    sw_Vconst: Use constant pairing interaction instead of FLEX vertices (default False)
    @param         eps: Convergence threshold on |Δ_new-Δ|/|Δ|
    @param      itemax: Maximum SCF iterations
    @param      m_diis: DIIS history length (m_diis<=1 → linear mixing only; default 5)
    @retval      delta: Converged gap function [Norb, Norb, Nw, Nk] complex128
                        (same array as delta_init when contiguous; otherwise a new copy)
    @retval     sigmak: FLEX self-energy [Norb, Norb, Nw, Nk] complex128
                        (zero array when sw_sigma=False)
    """
    if delta_init.dtype != np.complex128 or not delta_init.flags['C_CONTIGUOUS']:
        delta = np.ascontiguousarray(delta_init, dtype=np.complex128)
    else:
        delta = delta_init
    Norb, _, Nw, Nk = delta.shape
    Nkall, Nchi = len(kmap), len(Smat)
    sigmak = np.zeros((Norb, Norb, Nw, Nk), dtype=np.complex128)
    _lib.eliashberg.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),  # delta (inout)
        np.ctypeslib.ndpointer(dtype=np.complex128),  # sigmak (out/inout)
        np.ctypeslib.ndpointer(dtype=np.complex128),  # Gk (inout)
        np.ctypeslib.ndpointer(dtype=np.complex128),  # hamk
        np.ctypeslib.ndpointer(dtype=np.float64),     # Smat
        np.ctypeslib.ndpointer(dtype=np.float64),     # Cmat
        np.ctypeslib.ndpointer(dtype=np.int64),       # olist
        np.ctypeslib.ndpointer(dtype=np.float64),     # plist (= prt in Fortran)
        np.ctypeslib.ndpointer(dtype=np.int64),       # kmap
        np.ctypeslib.ndpointer(dtype=np.int64),       # invk
        POINTER(c_double), POINTER(c_double), POINTER(c_double),  # mu, temp, eps
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64),     # Nkall, Nk, Nw
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64),     # Nchi, Norb, Nx
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64),     # Ny, Nz, itemax
        POINTER(c_int64),                                         # gap_sym
        POINTER(c_bool), POINTER(c_bool),                         # sw_sigma,sw_Vconst
        POINTER(c_int64),                                          # m_diis
    ]
    _lib.eliashberg.restype = None
    _lib.eliashberg(delta, sigmak, Gk, hamk, Smat, Cmat, olist, plist,
                    kmap, invk,
                    byref(c_double(mu)), byref(c_double(temp)), byref(c_double(eps)),
                    byref(c_int64(Nkall)), byref(c_int64(Nk)), byref(c_int64(Nw)),
                    byref(c_int64(Nchi)), byref(c_int64(Norb)), byref(c_int64(Nx)),
                    byref(c_int64(Ny)), byref(c_int64(Nz)), byref(c_int64(itemax)),
                    byref(c_int64(gap_sym)), byref(c_bool(sw_sigma)), byref(c_bool(sw_Vconst)),
                    byref(c_int64(m_diis)))
    return delta, sigmak

def conv_delta_orb_to_band(delta: np.ndarray, uni: np.ndarray, invk: np.ndarray,
                           plist: np.ndarray, gap_sym) -> np.ndarray:
    """
    @fn conv_delta_orb_to_band
    @brief Convert the superconducting gap function from orbital basis to band basis (without SOC).
    @param   delta: Gap function in orbital basis [Norb, Norb, Nw, Nkall] complex128
    @param     uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param    invk: Inverse k-point mapping [Nkall] int64
    @param   plist: Parity list for gap symmetry [Norb] float64
    @param gap_sym: Gap symmetry index
    @return deltab: Gap function in band basis [Norb, Norb, Nkall] complex128
    """
    Nkall, Nk, Nw, Norb = len(invk), len(uni), len(delta[0, 0]), len(delta)
    deltab = np.zeros((Norb, Norb, Nkall), dtype=np.complex128)
    _lib.conv_delta_orb_to_band.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.conv_delta_orb_to_band.restype = None
    _lib.conv_delta_orb_to_band(deltab, delta, uni, plist, invk, byref(c_int64(Norb)),
                                byref(c_int64(Nkall)), byref(c_int64(Nk)),
                                byref(c_int64(Nw)), byref(c_int64(gap_sym)))
    return deltab

def conv_delta_orb_to_band_soc(delta: np.ndarray, uni: np.ndarray, invk: np.ndarray,
                               invs: np.ndarray, slist: np.ndarray) -> np.ndarray:
    """
    @fn conv_delta_orb_to_band_soc
    @brief Convert the superconducting gap function from orbital basis to band basis with SOC.
    @param   delta: Gap function in orbital basis [Norb, Norb, Nw, Nkall] complex128
    @param     uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param    invk: Inverse k-point mapping [Nkall] int64
    @param    invs: Inverse spin index mapping [Norb] int64
    @param   slist: Spin-orbit label for each orbital [Norb] int64
    @return deltab: Gap function in band basis [Norb, Norb, Nkall] complex128
    """
    Nkall, Nk, Nw, Norb = len(invk), len(uni), len(delta[0, 0]), len(delta)
    deltab = np.zeros((Norb, Norb, Nkall), dtype=np.complex128)
    _lib.conv_delta_orb_to_band_soc.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.conv_delta_orb_to_band_soc.restype = None
    _lib.conv_delta_orb_to_band_soc(deltab, delta, uni, invk, invs, slist, byref(c_int64(Norb)),
                                    byref(c_int64(Nkall)), byref(c_int64(Nk)), byref(c_int64(Nw)))
    return deltab

def gen_Fk(Gk: np.ndarray, delta: np.ndarray, invk: np.ndarray) -> np.ndarray:
    """
    @fn gen_Fk
    @brief Generate the anomalous Green's function Fk from G and delta (without SOC).
    @param     Gk: Normal Green's function [Norb, Norb, Nw, Nk] complex128
    @param  delta: Anomalous self-energy [Norb, Norb, Nw, Nkall] complex128
    @param   invk: Inverse k-point mapping [Nkall] int64
    @return    Fk: Anomalous Green's function [Norb, Norb, Nw, Nk] complex128
    """
    Nkall, Nk, Nw, Norb = len(invk), len(Gk[0, 0, 0]), len(delta[0, 0]), len(delta)
    Fk = np.zeros((Norb, Norb, Nw, Nk), dtype=np.complex128)
    _lib.mkfk_trs_nsoc_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.mkfk_trs_nsoc_.restype = None
    _lib.mkfk_trs_nsoc_(Fk, Gk, delta, byref(c_int64(Nk)), byref(c_int64(Nw)), byref(c_int64(Norb)))
    return Fk

def gen_Fk_soc(Gk: np.ndarray, delta: np.ndarray, invk: np.ndarray, invs: np.ndarray,
               slist: np.ndarray, gap_sym: int) -> np.ndarray:
    """
    @fn gen_Fk_soc
    @brief Generate the anomalous Green's function Fk from G and delta with SOC.
    @param     Gk: Normal Green's function [Norb, Norb, Nw, Nk] complex128
    @param  delta: Anomalous self-energy [Norb, Norb, Nw, Nkall] complex128
    @param   invk: Inverse k-point mapping [Nkall] int64
    @param   invs: Inverse spin index mapping [Norb] int64
    @param  slist: Spin-orbit label for each orbital [Norb] int64
    @param gap_sym: Gap symmetry index
    @return    Fk: Anomalous Green's function [Norb, Norb, Nw, Nkall] complex128
    """
    Nkall, Nk, Nw, Norb = len(invk), len(Gk[0, 0, 0]), len(delta[0, 0]), len(delta)
    Fk = np.zeros((Norb, Norb, Nw, Nkall), dtype=np.complex128)
    sgnsig = np.array([slist]).T.dot(np.array([slist])).astype(np.float64)
    _lib.mkfk_trs_soc_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64)
    ]
    _lib.mkfk_trs_soc_.restype = None
    _lib.mkfk_trs_soc_(Fk, Gk, delta, sgnsig, slist, invk, invs, byref(c_int64(Nkall)),
                       byref(c_int64(Nk)), byref(c_int64(Nw)), byref(c_int64(Norb)),
                       byref(c_int64(gap_sym)))
    return Fk

def get_initial_delta(init_delta: np.ndarray, uni: np.ndarray, kmap: np.ndarray,
                      invk: np.ndarray, Nw: int, gap_sym: int) -> np.ndarray:
    """
    @fn get_initial_delta
    @brief Convert a band-basis initial gap function into an orbital-basis gap function,
    broadcast over Matsubara frequencies, and normalize.
    Applies U * diag(init_delta) * U† at each k-point, then normalises by the total norm.
    @param  init_delta: Band-basis initial gap amplitudes [Norb, Nk] float64
                        (matches get_initial_gap output; C-order [Norb,Nk] == Fortran dimension(Nk,Norb))
    @param         uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param        kmap: Full k-grid to irreducible k-point mapping [Nkall, 3] int64
    @param        invk: Inverse k-point list [Nkall, 3] int64
    @param          Nw: Number of Matsubara frequency points
    @param     gap_sym: Gap symmetry index (0=s-wave diagonal, >0=singlet, <0=triplet)
    @return      delta: Orbital-basis gap function [Norb, Norb, Nw, Nk] complex128
    """
    Norb, Nk = init_delta.shape
    Nkall = len(kmap)
    delta = np.zeros((Norb, Norb, Nw, Nk), dtype=np.complex128)
    _lib.get_initial_delta_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64),
    ]
    _lib.get_initial_delta_.restype = None
    _lib.get_initial_delta_(delta, init_delta, uni, kmap, invk, byref(c_int64(Nkall)), byref(c_int64(Nk)),
                            byref(c_int64(Nw)), byref(c_int64(Norb)), byref(c_int64(gap_sym)))
    return delta

def remap_gap(delta0, plist, invk, gap_sym):
    """
    @fn remap_gap
    @brief Remap the anomalous self-energy from the reduced k-grid to the full k-grid applying gap symmetry.
    @param  delta0: Anomalous self-energy on reduced k-grid [Norb, Norb, Nw, Nk] complex128
    @param   plist: Parity list for gap symmetry [Norb] float64
    @param    invk: Inverse k-point mapping [Nkall] int64
    @param gap_sym: Gap symmetry index
    @return  delta: Anomalous self-energy on full k-grid [Norb, Norb, Nw, Nkall] complex128
    """
    Nkall, Nk, Norb = len(invk), len(delta0.T), len(plist)
    Nw = delta0.shape[2]
    delta = np.zeros((Norb, Norb, Nw, Nkall), dtype=np.complex128)
    _lib.remap_delta_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64)
    ]
    _lib.remap_delta_.restype = None
    _lib.remap_delta_(delta, delta0, plist, invk, byref(c_int64(Nkall)), byref(c_int64(Nk)),
                      byref(c_int64(Nw)), byref(c_int64(Norb)), byref(c_int64(gap_sym)))
    return delta
