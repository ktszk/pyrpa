from ctypes import *
import numpy as np
from ._loader import _lib

def calc_Lij(eig: np.ndarray, vk: np.ndarray, ffermi: np.ndarray, mu: float, w: float,
             idelta: float, temp: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    @fn calc_Lij
    @brief Evaluate the generalized thermoelectric transport tensors L11, L12, L22 at a given frequency.
    @param    eig: Eigenvalues [Nk, Norb] float64
    @param     vk: Group velocities [Nk, Norb, 3] complex128
    @param ffermi: Fermi-Dirac occupations [Nk, Norb] float64
    @param     mu: Chemical potential in eV
    @param      w: Frequency in eV
    @param idelta: Lorentzian broadening parameter in eV
    @param   temp: Temperature in eV
    @retval   L11: Charge conductivity tensor [3, 3] complex128
    @retval   L12: Thermoelectric tensor [3, 3] complex128
    @retval   L22: Thermal conductivity tensor [3, 3] complex128
    """
    Nk = len(eig)
    Norb = int(eig.size / Nk)
    L11 = np.zeros((3, 3), dtype=np.complex128)
    L12 = np.zeros((3, 3), dtype=np.complex128)
    L22 = np.zeros((3, 3), dtype=np.complex128)
    eps = idelta * 1e-3
    _lib.calc_lij.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_double),
        POINTER(c_double), POINTER(c_double),
        POINTER(c_double), POINTER(c_double),
    ]
    _lib.calc_lij.restype = c_void_p
    _lib.calc_lij(L11, L22, L12, vk, eig, ffermi, byref(c_int64(Norb)), byref(c_int64(Nk)),
                  byref(c_double(mu)), byref(c_double(w)), byref(c_double(idelta)),
                  byref(c_double(eps)), byref(c_double(temp)))
    return L11, L12, L22

def calc_Lij_wl(eig: np.ndarray, vk: np.ndarray, ffermi: np.ndarray, mu: float,
                wl: np.ndarray, idelta: float, temp: float
                ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    @fn calc_Lij_wl
    @brief Evaluate L11, L12, L22 for all frequencies wl in a single Fortran call.
    @param    eig: Eigenvalues [Nk, Norb] float64
    @param     vk: Group velocities [Nk, Norb, Norb, 3] complex128
    @param ffermi: Fermi-Dirac occupations [Nk, Norb] float64
    @param     mu: Chemical potential in eV
    @param     wl: Frequency mesh [Nw] float64
    @param idelta: Lorentzian broadening parameter in eV
    @param   temp: Temperature in eV
    @retval   L11: [Nw, 3, 3] complex128
    @retval   L12: [Nw, 3, 3] complex128
    @retval   L22: [Nw, 3, 3] complex128
    """
    Nk = len(eig)
    Norb = int(eig.size / Nk)
    Nw = len(wl)
    eps = idelta * 1e-3
    wl = np.ascontiguousarray(wl, dtype=np.float64)
    L11 = np.zeros((Nw, 3, 3), dtype=np.complex128, order='F')
    L12 = np.zeros((Nw, 3, 3), dtype=np.complex128, order='F')
    L22 = np.zeros((Nw, 3, 3), dtype=np.complex128, order='F')
    _lib.calc_lij_wl.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64),
        POINTER(c_double),
        np.ctypeslib.ndpointer(dtype=np.float64),
        POINTER(c_double), POINTER(c_double),
        POINTER(c_double),
    ]
    _lib.calc_lij_wl.restype = c_void_p
    _lib.calc_lij_wl(L11, L22, L12, vk, eig, ffermi,
                     byref(c_int64(Norb)), byref(c_int64(Nk)), byref(c_int64(Nw)),
                     byref(c_double(mu)), wl,
                     byref(c_double(idelta)), byref(c_double(eps)), byref(c_double(temp)))
    return L11, L12, L22

def calc_Kn(eig: np.ndarray, veloc: np.ndarray, kweight: np.ndarray, temp: float,
            mu: float, tau: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    @fn calc_Kn
    @brief Compute Boltzmann transport coefficients K0, K1, K2 for conductivity and thermopower.
    @param    eig: Eigenvalues [Nk, Norb] float64
    @param  veloc: Group velocities [Nk, Norb, 3] float64
    @param kweight: k-point weights [Nk] float64
    @param   temp: Temperature in eV
    @param     mu: Chemical potential in eV
    @param    tau: Relaxation time [Nk, Norb] float64
    @retval    K0: Charge transport tensor [3, 3] float64
    @retval    K1: Thermoelectric transport tensor [3, 3] float64
    @retval    K2: Thermal transport tensor [3, 3] float64
    """
    Nk = len(eig)
    Norb = int(eig.size / Nk)
    K0 = np.zeros((3, 3), dtype=np.float64)
    K1 = np.zeros((3, 3), dtype=np.float64)
    K2 = np.zeros((3, 3), dtype=np.float64)
    _lib.calc_kn.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        POINTER(c_double), POINTER(c_double),
        POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.calc_kn.restype = c_void_p
    _lib.calc_kn(K0, K1, K2, eig, veloc, kweight, tau, byref(c_double(temp)),
                 byref(c_double(mu)), byref(c_int64(Nk)), byref(c_int64(Norb)))
    return K0, K1, K2

def calc_sigmahall(eig: np.ndarray, veloc: np.ndarray, imass: np.ndarray,
                   kweight: np.ndarray, tau: np.ndarray, temp: float, mu: float) -> float:
    """
    @fn calc_sigmahall
    @brief Compute the Hall conductivity sigma_Hall via band integrals using velocity and inverse mass.
    @param     eig: Eigenvalues [Nk, Norb] float64
    @param   veloc: Group velocities [Nk, Norb, 3] float64
    @param   imass: Inverse effective mass tensor [Nk, Norb, 3, 3] float64
    @param kweight: k-point weights [Nk] float64
    @param     tau: Relaxation time [Nk, Norb] float64
    @param    temp: Temperature in eV
    @param      mu: Chemical potential in eV
    @return sigma_hall: Hall conductivity (scalar) float64
    """
    Nk = len(eig)
    Norb = int(eig.size / Nk)
    sigma_hall = c_double(0.0)
    _lib.calc_sigma_hall.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        POINTER(c_double), POINTER(c_double),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_double)
    ]
    _lib.calc_sigma_hall.restype = c_void_p
    _lib.calc_sigma_hall(eig, veloc, imass, kweight, tau, byref(c_double(temp)), byref(c_double(mu)),
                         byref(c_int64(Nk)), byref(c_int64(Norb)), byref(sigma_hall))
    return sigma_hall.value

def calc_tdf(eig: np.ndarray, veloc: np.ndarray, kweight: np.ndarray,
             tau: np.ndarray, Nw: int) -> np.ndarray:
    """
    @fn calc_tdf
    @brief Compute the transport distribution function (TDF) tensor as a function of energy.
    @param    eig: Eigenvalues [Nk, Norb] float64
    @param  veloc: Group velocities [Nk, Norb, 3] float64
    @param kweight: k-point weights [Nk] float64
    @param    tau: Relaxation time [Nk, Norb] float64
    @param     Nw: Number of energy points
    @return   tdf: Transport distribution function [Nw, 3, 3] float64
    """
    Nk = len(eig)
    Norb = int(eig.size / Nk)
    tdf = np.zeros((Nw, 3, 3), dtype=np.float64)
    _lib.calc_tdf.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64)
    ]
    _lib.calc_tdf.restype = c_void_p
    _lib.calc_tdf(tdf, eig, veloc, kweight, tau, byref(c_int64(Nw)),
                  byref(c_int64(Nk)), byref(c_int64(Norb)))
    return tdf

def get_tau(tauw: np.ndarray, eig: np.ndarray, tau_max: float, tau_mode: int, eps: float = 1.0e-4) -> np.ndarray:
    """
    @fn get_tau
    @brief Estimate k-resolved relaxation time tau(k) from frequency-dependent scattering weights.
    @param    tauw: Frequency-dependent scattering weight [Nw] float64
    @param     eig: Eigenvalues [Nk, Norb] float64
    @param tau_max: Upper bound for relaxation time values
    @param tau_mode: Mode selector for the Fortran tau integral equation
    @param     eps: Convergence tolerance
    @return    tau: k-resolved relaxation time [Nk, Norb] float64
    """
    Nk, Nw = len(eig), len(tauw)
    Norb = int(eig.size / Nk)
    tau = np.zeros((Nk, Norb), dtype=np.float64)
    _lib.get_tau.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        POINTER(c_double), POINTER(c_double),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.get_tau.restype = c_void_p
    _lib.get_tau(tau, tauw, eig, byref(c_double(tau_max)), byref(c_double(eps)),
                 byref(c_int64(tau_mode)), byref(c_int64(Nk)), byref(c_int64(Nw)),
                 byref(c_int64(Norb)))
    return tau

def calc_tau_epa(eig: np.ndarray, gavg: np.ndarray, wavg: np.ndarray,
                 edge: np.ndarray, step: np.ndarray, nbin: np.ndarray,
                 mu: float, temp: float) -> np.ndarray:
    """
    @fn calc_tau_epa
    @brief Compute EPA relaxation time from epa.x (job='egrid') output.
    @param   eig: Eigenvalues [Nk, Norb] float64 (eV)
    @param  gavg: EPA averaged |g|^2 [ngrid, nbin_max, nbin_max, nmodes] float64 (eV^2)
    @param  wavg: Averaged phonon frequencies per mode [nmodes] float64 (eV)
    @param  edge: Grid edges [ngrid] float64 (eV)
    @param  step: Grid steps [ngrid] float64 (eV)
    @param  nbin: Number of bins per grid [ngrid] int64
    @param    mu: Chemical potential (eV)
    @param  temp: Temperature kB*T (eV)
    @return  tau: Relaxation time [Nk, Norb] float64
    """
    Nk = len(eig)
    Norb = int(eig.size / Nk)
    ngrid = len(edge)
    nmodes = len(wavg)
    nbin_max = int(np.max(nbin))
    tau = np.zeros((Nk, Norb), dtype=np.float64)
    nbin_i64 = np.ascontiguousarray(nbin, dtype=np.int64)
    _lib.calc_tau_epa.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        POINTER(c_double), POINTER(c_double),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.calc_tau_epa.restype = c_void_p
    _lib.calc_tau_epa(tau, gavg, wavg, eig, edge, step,
                      byref(c_double(mu)), byref(c_double(temp)),
                      byref(c_int64(Nk)), byref(c_int64(Norb)),
                      byref(c_int64(nmodes)), nbin_i64,
                      byref(c_int64(ngrid)), byref(c_int64(nbin_max)))
    return tau
