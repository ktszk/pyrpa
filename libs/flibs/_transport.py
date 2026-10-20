import ctypes
import numpy as np
import libs.flibs._loader as _loader

# Public API of this submodule; keeps `from ._x import *` in __init__.py from
# re-exporting numpy, scipy, ctypes and the sibling packages imported above.
__all__ = [
    'calc_Lij', 'calc_Lij_wl', 'calc_Kn', 'calc_sigma_hall', 'calc_k0_band', 'calc_tdf', 'calc_tau_epa', 'epa_dos'
]

# --- ctypes signatures: set once at import.
# All Fortran entry points are subroutines, so restype is always None.
_loader._lib.calc_lij.argtypes = [
    _loader.c_array(np.complex128),  # L11 (out)
    _loader.c_array(np.complex128),  # L22 (out)
    _loader.c_array(np.complex128),  # L12 (out)
    _loader.c_array(np.complex128),  # vk
    _loader.c_array(np.float64),     # eig
    _loader.c_array(np.float64),     # ffermi
    _loader.c_array(np.float64),     # kweight
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Norb, Nk
    ctypes.POINTER(ctypes.c_double),                            # mu
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),         # w, idelta
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),         # eps, temp
]
_loader._lib.calc_lij.restype = None
_loader._lib.calc_lij_wl.argtypes = [
    _loader.f_array(np.complex128),           # L11 (out, Fortran order [Nw,3,3])
    _loader.f_array(np.complex128),           # L22 (out, Fortran order [Nw,3,3])
    _loader.f_array(np.complex128),           # L12 (out, Fortran order [Nw,3,3])
    _loader.c_array(np.complex128),           # vk
    _loader.c_array(np.float64),              # eig
    _loader.c_array(np.float64),              # ffermi
    _loader.c_array(np.float64),              # kweight
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),  # Norb, Nk, Nw
    ctypes.POINTER(ctypes.c_double),                                     # mu
    _loader.c_array(np.float64),              # wl
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),                  # idelta, eps
    ctypes.POINTER(ctypes.c_double),                                     # temp
]
_loader._lib.calc_lij_wl.restype = None
_loader._lib.calc_kn.argtypes = [
    _loader.c_array(np.float64),  # K0 (out)
    _loader.c_array(np.float64),  # K1 (out)
    _loader.c_array(np.float64),  # K2 (out)
    _loader.c_array(np.float64),  # eig
    _loader.c_array(np.float64),  # veloc
    _loader.c_array(np.float64),  # kweight
    _loader.c_array(np.float64),  # tau
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),      # temp, mu
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)         # Nk, Norb
]
_loader._lib.calc_kn.restype = None
_loader._lib.calc_sigma_hall.argtypes = [
    _loader.c_array(np.float64),  # shall (out)
    _loader.c_array(np.float64),  # sabs (out)
    _loader.c_array(np.float64),  # eig
    _loader.c_array(np.float64),  # veloc
    _loader.c_array(np.float64),  # imass
    _loader.c_array(np.float64),  # kweight
    _loader.c_array(np.float64),  # tau
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),      # temp, mu
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)         # Nk, Norb
]
_loader._lib.calc_sigma_hall.restype = None
_loader._lib.calc_k0_band.argtypes = [
    _loader.c_array(np.float64),  # K0n (out)
    _loader.c_array(np.float64),  # eig
    _loader.c_array(np.float64),  # veloc
    _loader.c_array(np.float64),  # kweight
    _loader.c_array(np.float64),  # tau
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),      # temp, mu
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)         # Nk, Norb
]
_loader._lib.calc_k0_band.restype = None
_loader._lib.calc_tdf.argtypes = [
    _loader.c_array(np.float64),  # tdf (out)
    _loader.c_array(np.float64),  # eig
    _loader.c_array(np.float64),  # veloc
    _loader.c_array(np.float64),  # kweight
    _loader.c_array(np.float64),  # tau
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),        # Nw, Nk
    ctypes.POINTER(ctypes.c_int64)                           # Norb
]
_loader._lib.calc_tdf.restype = None
_loader._lib.calc_tau_epa.argtypes = [
    _loader.c_array(np.float64),  # tau (out)
    _loader.c_array(np.float64),  # gavg
    _loader.c_array(np.float64),  # wavg
    _loader.c_array(np.float64),  # eig
    _loader.c_array(np.float64),  # edge
    _loader.c_array(np.float64),  # step
    _loader.c_array(np.float64),  # dos
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),      # e0, de
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),      # mu, temp
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),        # Nk, Norb
    ctypes.POINTER(ctypes.c_int64),                          # nmodes
    _loader.c_array(np.int64),    # nbin
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),        # ngrid, nbin_max
    ctypes.POINTER(ctypes.c_int64)                           # Ndos
]
_loader._lib.calc_tau_epa.restype = None


def calc_Lij(eig: np.ndarray, vk: np.ndarray, ffermi: np.ndarray, mu: float, w: float,
             idelta: float, temp: float,
             kweight: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    @brief Evaluate the generalized thermoelectric transport tensors L11, L12, L22 at a given frequency.
    @param    eig: Eigenvalues [Nk, Norb] float64
    @param     vk: Group velocities [Nk, Norb, 3] complex128
    @param ffermi: Fermi-Dirac occupations [Nk, Norb] float64
    @param     mu: Chemical potential in eV
    @param      w: Frequency in eV
    @param idelta: Lorentzian broadening parameter in eV
    @param   temp: Temperature in eV
    @param kweight: k-point weights [Nk] (None: uniform mesh, all 1); the sum is normalized by sum(kweight)
    @return   L11: Charge conductivity tensor [3, 3] complex128
    @return   L12: Thermoelectric tensor [3, 3] complex128
    @return   L22: Thermal conductivity tensor [3, 3] complex128
    """
    Nk = len(eig)
    Norb = eig.shape[1]
    L11 = np.zeros((3, 3), dtype=np.complex128)
    L12 = np.zeros((3, 3), dtype=np.complex128)
    L22 = np.zeros((3, 3), dtype=np.complex128)
    eps = idelta * 1e-3
    kweight = np.ones(Nk) if kweight is None else np.ascontiguousarray(kweight, dtype=np.float64)
    _loader._lib.calc_lij(L11, L22, L12, vk, eig, ffermi, kweight, _loader.i64(Norb), _loader.i64(Nk),
                          _loader.dbl(mu), _loader.dbl(w), _loader.dbl(idelta),
                          _loader.dbl(eps), _loader.dbl(temp))
    return L11, L12, L22


def calc_Lij_wl(eig: np.ndarray, vk: np.ndarray, ffermi: np.ndarray, mu: float,
                wl: np.ndarray, idelta: float, temp: float, kweight: np.ndarray | None = None
                ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    @brief Evaluate L11, L12, L22 for all frequencies wl in a single Fortran call.
    @param    eig: Eigenvalues [Nk, Norb] float64
    @param     vk: Group velocities [Nk, Norb, Norb, 3] complex128
    @param ffermi: Fermi-Dirac occupations [Nk, Norb] float64
    @param     mu: Chemical potential in eV
    @param     wl: Frequency mesh [Nw] float64
    @param idelta: Lorentzian broadening parameter in eV
    @param   temp: Temperature in eV
    @param kweight: k-point weights [Nk] (None: uniform mesh, all 1); the sum is normalized by sum(kweight)
    @return   L11: [Nw, 3, 3] complex128
    @return   L12: [Nw, 3, 3] complex128
    @return   L22: [Nw, 3, 3] complex128
    """
    Nk = len(eig)
    Norb = eig.shape[1]
    Nw = len(wl)
    eps = idelta * 1e-3
    wl = np.ascontiguousarray(wl, dtype=np.float64)
    kweight = np.ones(Nk) if kweight is None else np.ascontiguousarray(kweight, dtype=np.float64)
    L11 = np.zeros((Nw, 3, 3), dtype=np.complex128, order='F')
    L12 = np.zeros((Nw, 3, 3), dtype=np.complex128, order='F')
    L22 = np.zeros((Nw, 3, 3), dtype=np.complex128, order='F')
    _loader._lib.calc_lij_wl(L11, L22, L12, vk, eig, ffermi, kweight,
                             _loader.i64(Norb), _loader.i64(Nk), _loader.i64(Nw),
                             _loader.dbl(mu), wl,
                             _loader.dbl(idelta), _loader.dbl(eps), _loader.dbl(temp))
    return L11, L12, L22


def calc_Kn(eig: np.ndarray, veloc: np.ndarray, kweight: np.ndarray, temp: float,
            mu: float, tau: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    @brief Compute Boltzmann transport coefficients K0, K1, K2 for conductivity and thermopower.
    @param    eig: Eigenvalues [Nk, Norb] float64
    @param  veloc: Group velocities [Nk, Norb, 3] float64
    @param kweight: k-point weights [Nk] float64
    @param   temp: Temperature in eV
    @param     mu: Chemical potential in eV
    @param    tau: Relaxation time [Nk, Norb] float64
    @return    K0: Charge transport tensor [3, 3] float64
    @return    K1: Thermoelectric transport tensor [3, 3] float64
    @return    K2: Thermal transport tensor [3, 3] float64
    """
    Nk = len(eig)
    Norb = eig.shape[1]
    K0 = np.zeros((3, 3), dtype=np.float64)
    K1 = np.zeros((3, 3), dtype=np.float64)
    K2 = np.zeros((3, 3), dtype=np.float64)
    _loader._lib.calc_kn(K0, K1, K2, eig, veloc, kweight, tau, _loader.dbl(temp),
                         _loader.dbl(mu), _loader.i64(Nk), _loader.i64(Norb))
    return K0, K1, K2


def calc_sigma_hall(eig: np.ndarray, veloc: np.ndarray, imass: np.ndarray,
                    kweight: np.ndarray, tau: np.ndarray, temp: float,
                    mu: float) -> tuple[np.ndarray, np.ndarray]:
    """
    @brief Weak-field Boltzmann Hall kernel for all three field directions, band
    resolved and antisymmetrized in the two current indices:

        S[n,a,b,g] = 0.5*(A[n,a,b,g] - A[n,b,a,g]),
        A[n,a,b,g] = sum_k w_k tau^2 (-df/de) eps_{g,u,v} v_a M^-1_{b,u} v_v

    sigma^(1)_{ab}(B||g)/B is S[...,a,b,g] times the same prefactor that turns
    calc_Kn's K0 into sigma, so the caller applies 1/(sum(kweight)*Vuc) once.
    Contracted at (a,b,g)=(x,y,z) this is minus the textbook one-sided scalar
    -0.5 w (vx^2/m_yy + vy^2/m_xx - 2 vx vy/m_xy) term by term.  The tensor form
    is what lets the caller build rho^(1) = -sigma^-1 sigma^(1) sigma^-1 for a
    cell whose axes are not orthogonal (sigma_xy != 0 already at B=0) and for an
    arbitrary field direction; the band resolution feeds the two-fluid
    decomposition needed in a compensated metal.

    Sabs is sum |contribution| with the same indexing; |S|/Sabs is the
    cancellation ratio of the sign-changing Hall integrand, which is what
    predicts an R_H whose sign still moves with the k-mesh.
    @param     eig: Eigenvalues [Nk, Norb] float64 (eV)
    @param   veloc: Group velocities [Nk, Norb, 3] float64 (m/s)
    @param   imass: Inverse effective mass tensor [Nk, Norb, 3, 3] float64 (1/kg)
    @param kweight: k-point weights [Nk] float64
    @param     tau: Relaxation time [Nk, Norb] float64
    @param    temp: Temperature in eV
    @param      mu: Chemical potential in eV
    @return Shall: Hall kernel [Norb, 3, 3, 3] float64, indexed [n, a, b, g]
    @return  Sabs: sum of |contribution|, same shape
    """
    Nk = len(eig)
    Norb = eig.shape[1]
    Shall = np.zeros((Norb, 3, 3, 3), dtype=np.float64)
    Sabs = np.zeros((Norb, 3, 3, 3), dtype=np.float64)
    _loader._lib.calc_sigma_hall(Shall, Sabs, eig, veloc, imass, kweight, tau,
                                 _loader.dbl(temp), _loader.dbl(mu), _loader.i64(Nk), _loader.i64(Norb))
    return Shall, Sabs


def calc_k0_band(eig: np.ndarray, veloc: np.ndarray, kweight: np.ndarray,
                 temp: float, mu: float, tau: np.ndarray) -> np.ndarray:
    """
    @brief Per-band charge transport kernel K0n[n,i,j] = sum_k v_i v_j tau (-df/de) w_k.

    calc_Kn returns only the band sum; the per-band split is what allows the
    two-fluid decomposition sigma = sum_n sigma_n used to read a Hall
    coefficient in a compensated metal.  Sums over n to calc_Kn's K0.
    @param     eig: Eigenvalues [Nk, Norb] float64 (eV)
    @param   veloc: Group velocities [Nk, Norb, 3] float64 (m/s)
    @param kweight: k-point weights [Nk] float64
    @param    temp: Temperature in eV
    @param      mu: Chemical potential in eV
    @param     tau: Relaxation time [Nk, Norb] float64
    @return    K0n: [Norb, 3, 3] float64
    """
    Nk = len(eig)
    Norb = eig.shape[1]
    K0n = np.zeros((Norb, 3, 3), dtype=np.float64)
    _loader._lib.calc_k0_band(K0n, eig, veloc, kweight, tau, _loader.dbl(temp), _loader.dbl(mu),
                              _loader.i64(Nk), _loader.i64(Norb))
    return K0n


def calc_tdf(eig: np.ndarray, veloc: np.ndarray, kweight: np.ndarray,
             tau: np.ndarray, Nw: int) -> np.ndarray:
    """
    @brief Compute the transport distribution function (TDF) tensor as a function of energy.

    The energy grid is the Nw bin CENTRES E_i = min(eig) + (i-1/2)*dw,
    dw = (max(eig)-min(eig))/Nw.  Like calc_Kn this returns the kweight-weighted
    SUM over k (no 1/Nk), so int dE (-df/dE) tdf(E) equals calc_Kn's K0 and the
    caller must apply 1/(Nk*Vuc) exactly once.
    @param    eig: Eigenvalues [Nk, Norb] float64
    @param  veloc: Group velocities [Nk, Norb, 3] float64
    @param kweight: k-point weights [Nk] float64
    @param    tau: Relaxation time [Nk, Norb] float64
    @param     Nw: Number of energy points
    @return   tdf: Transport distribution function [Nw, 3, 3] float64
    """
    Nk = len(eig)
    Norb = eig.shape[1]
    tdf = np.zeros((Nw, 3, 3), dtype=np.float64)
    _loader._lib.calc_tdf(tdf, eig, veloc, kweight, tau, _loader.i64(Nw),
                          _loader.i64(Nk), _loader.i64(Norb))
    return tdf


def epa_dos(eig: np.ndarray, step: np.ndarray, nsub: int = 8) -> tuple[np.ndarray, float, float]:
    """
    @brief Electronic DOS per spin per unit cell for the EPA rate, on a grid nsub times finer
    than the finest EPA bin (linear, cloud-in-cell deposition of the eigenvalues, so the
    integral is the number of bands).  The EPA bins themselves are far too coarse for rho:
    bin-averaging it removes the band-edge onset that makes tau long there.
    @param  eig: Eigenvalues [Nk, Norb] (eV) on a uniform k-mesh
    @param step: EPA grid steps [ngrid] (eV)
    @return (dos [Ndos] (1/eV), e0, de): values on e0 + i*de
    """
    de = float(np.abs(step).min())/nsub
    e0 = float(eig.min())-de
    x = (np.ravel(eig)-e0)/de
    i = np.floor(x).astype(np.int64)
    t = x-i
    dos = np.zeros(int(i.max())+3)
    np.add.at(dos, i, 1.0-t)
    np.add.at(dos, i+1, t)
    return dos/(len(eig)*de), e0, de


def calc_tau_epa(eig: np.ndarray, gavg: np.ndarray, wavg: np.ndarray,
                 edge: np.ndarray, step: np.ndarray, nbin: np.ndarray,
                 mu: float, temp: float) -> np.ndarray:
    """
    @brief Compute EPA relaxation time from epa.x (job='egrid') output, using the
    energy-conserving EPA rate (Samsonidze-Kozinsky) as their BoltzTraP epa_calc does:
    1/tau(e) = 2*pi*sum_nu { g2(e,e+w)[nB+f(e+w)]rho(e+w) + g2(e,e-w)[nB+1-f(e-w)]rho(e-w) },
    with rho the electronic DOS per spin per unit cell (epa_dos, interpolated at e +- w),
    g2 constant per EPA bin, bins clamped into the grids and tau <= 1 ns.
    @param   eig: Eigenvalues [Nk, Norb] float64 (eV), uniform k-mesh
    @param  gavg: EPA averaged |g|^2 [ngrid, initial bin, final bin, nmodes] float64 (eV^2)
    @param  wavg: Averaged phonon frequencies per mode [nmodes] float64 (eV)
    @param  edge: Grid edges [ngrid] float64 (eV)
    @param  step: Grid steps [ngrid] float64 (eV)
    @param  nbin: Number of bins per grid [ngrid] int64
    @param    mu: Chemical potential (eV)
    @param  temp: Temperature kB*T (eV)
    @return  tau: Relaxation time [Nk, Norb] float64 (hbar/eV units)
    """
    eig = np.ascontiguousarray(eig, dtype=np.float64)
    Nk = len(eig)
    Norb = eig.shape[1]
    ngrid = len(edge)
    nmodes = len(wavg)
    nbin_max = int(np.max(nbin))
    dos, e0, de = epa_dos(eig, step)
    tau = np.zeros((Nk, Norb), dtype=np.float64)
    nbin_i64 = np.ascontiguousarray(nbin, dtype=np.int64)
    _loader._lib.calc_tau_epa(tau, np.ascontiguousarray(gavg, dtype=np.float64),
                              np.ascontiguousarray(wavg, dtype=np.float64), eig,
                              np.ascontiguousarray(edge, dtype=np.float64),
                              np.ascontiguousarray(step, dtype=np.float64), dos,
                              _loader.dbl(e0), _loader.dbl(de),
                              _loader.dbl(mu), _loader.dbl(temp),
                              _loader.i64(Nk), _loader.i64(Norb),
                              _loader.i64(nmodes), nbin_i64,
                              _loader.i64(ngrid), _loader.i64(nbin_max), _loader.i64(len(dos)))
    return tau
