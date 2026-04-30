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
    Nw = int(delta0.size / (Nk * Norb * Norb))
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
