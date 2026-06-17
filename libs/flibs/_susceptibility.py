from ctypes import *
import numpy as np
from ._loader import _lib, i64, dbl

# --- ctypes signatures: set once at import.
# All Fortran entry points are subroutines, so restype is always None.
_lib.get_chi0_conv_.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),  # chi
    np.ctypeslib.ndpointer(dtype=np.complex128),  # Gk
    np.ctypeslib.ndpointer(dtype=np.int64),       # kmap
    np.ctypeslib.ndpointer(dtype=np.int64),       # invk
    np.ctypeslib.ndpointer(dtype=np.int32),       # irr_chi
    np.ctypeslib.ndpointer(dtype=np.int32),       # chi_map
    np.ctypeslib.ndpointer(dtype=np.int64),       # olist
    POINTER(c_double), POINTER(c_int64),
    POINTER(c_int64), POINTER(c_int64),
    POINTER(c_int64), POINTER(c_int64),
    POINTER(c_int64), POINTER(c_int64),
    POINTER(c_int64)
]
_lib.get_chi0_conv_.restype = None
_lib.get_chi_irr.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.int64),
    np.ctypeslib.ndpointer(dtype=np.int64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    POINTER(c_int64), POINTER(c_int64),
    POINTER(c_int64), POINTER(c_int64),
    POINTER(c_double),
    POINTER(c_double), POINTER(c_double),
]
_lib.get_chi_irr.restype = None
_lib.chiq_map.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.int64),
    POINTER(c_double), POINTER(c_double),
    POINTER(c_double), POINTER(c_double),
    POINTER(c_int64), POINTER(c_int64),
    POINTER(c_int64), POINTER(c_int64), POINTER(c_int64)
]
_lib.chiq_map.restype = None
_lib.get_tr_chi.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.int64),
    POINTER(c_int64), POINTER(c_int64),
    POINTER(c_int64)
]
_lib.get_tr_chi.restype = None
_lib.get_phi_irr.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.int64),
    np.ctypeslib.ndpointer(dtype=np.int64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    POINTER(c_int64), POINTER(c_int64),
    POINTER(c_int64), POINTER(c_int64),
    POINTER(c_double), POINTER(c_double),
    POINTER(c_double), POINTER(c_double),
]
_lib.get_phi_irr.restype = None
_lib.phiq_map.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.int64),
    POINTER(c_double), POINTER(c_double),
    POINTER(c_double), POINTER(c_double),
    POINTER(c_double),
    POINTER(c_int64), POINTER(c_int64),
    POINTER(c_int64), POINTER(c_int64),
    POINTER(c_int64), POINTER(c_bool)
]
_lib.phiq_map.restype = None
_lib.get_tr_phi.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.int64),
    POINTER(c_int64), POINTER(c_int64),
    POINTER(c_int64)
]
_lib.get_tr_phi.restype = None
_lib.get_chis.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.float64),
    POINTER(c_int64), POINTER(c_int64),
]
_lib.get_chis.restype = None
_lib.get_chi0.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.int64),
    np.ctypeslib.ndpointer(dtype=np.int64),
    np.ctypeslib.ndpointer(dtype=np.int64),
    POINTER(c_double), POINTER(c_int64),
    POINTER(c_int64), POINTER(c_int64),
    POINTER(c_int64), POINTER(c_int64),
    POINTER(c_int64), POINTER(c_int64),
    POINTER(c_int64),
    POINTER(c_double)  # maxchi0s_out (Stoner factor)
]
_lib.get_chi0.restype = None
_lib.get_chi0_soc.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.int64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.int64),
    np.ctypeslib.ndpointer(dtype=np.int64),
    np.ctypeslib.ndpointer(dtype=np.int64),
    np.ctypeslib.ndpointer(dtype=np.int64),
    np.ctypeslib.ndpointer(dtype=np.int64),
    POINTER(c_double), POINTER(c_int64),
    POINTER(c_int64), POINTER(c_int64),
    POINTER(c_int64), POINTER(c_int64),
    POINTER(c_int64), POINTER(c_int64),
    POINTER(c_int64)
]
_lib.get_chi0_soc.restype = None
_lib.get_chis_chic_soc.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.int64),
    np.ctypeslib.ndpointer(dtype=np.int64),
    np.ctypeslib.ndpointer(dtype=np.int64),
    np.ctypeslib.ndpointer(dtype=np.int64),
    POINTER(c_int64), POINTER(c_int64),
    POINTER(c_int64), POINTER(c_int64)
]
_lib.get_chis_chic_soc.restype = None
_lib.get_chis_chic.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    POINTER(c_int64), POINTER(c_int64), POINTER(c_int64)
]
_lib.get_chis_chic.restype = None
_lib.get_chi_irr_sc.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.int64),
    np.ctypeslib.ndpointer(dtype=np.int64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    POINTER(c_int64), POINTER(c_int64),
    POINTER(c_int64), POINTER(c_int64),
    POINTER(c_double), POINTER(c_double),
    POINTER(c_double), POINTER(c_bool),
]
_lib.get_chi_irr_sc.restype = None
_lib.get_eig_or_tr_chi.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.complex128),
    np.ctypeslib.ndpointer(dtype=np.int64),
    POINTER(c_int64), POINTER(c_int64),
    POINTER(c_int64), POINTER(c_bool)
]
_lib.get_eig_or_tr_chi.restype = None

def _get_chi_map(olist: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Python port of fchi.f90 get_chi_map (not bind(C), so not callable directly):
    index map for the exchange symmetry chi_{1234}(q,iw) = chi*_{4321}(q,iw).
    Returns (chi_map, irr_chi) as int32 arrays whose C memory layout matches the
    Fortran dummies chi_map(Nchi,Nchi,2) and irr_chi(Nchi*(Nchi+1)/2,2),
    with 1-based index values as the Fortran side expects.
    """
    ol = np.asarray(olist, dtype=np.int64)
    Nchi = len(ol)
    chi_map = np.zeros((2, Nchi, Nchi), dtype=np.int32)   # [c,l1,m1] == Fortran chi_map(m1,l1,c)
    chi_irr = np.zeros((Nchi, Nchi), dtype=np.int32)
    for l1 in range(Nchi):
        for m1 in range(Nchi):
            t1 = (ol[l1, 0], ol[l1, 1], ol[m1, 0], ol[m1, 1])
            found = False
            for l2 in range(Nchi):
                for m2 in range(Nchi):
                    t2 = (ol[m2, 1], ol[m2, 0], ol[l2, 1], ol[l2, 0])
                    if t1 == t2:
                        chi_map[0, l1, m1] = m2 + 1
                        chi_map[1, l1, m1] = l2 + 1
                        if chi_irr[l2, m2] == 0:
                            chi_irr[l1, m1] = 1
                        found = True
                        break
                if found:
                    break
    irr_chi = np.zeros((2, Nchi * (Nchi + 1) // 2), dtype=np.int32)  # [c,i] == Fortran irr_chi(i,c)
    it = 0
    for l1 in range(Nchi):
        for m1 in range(Nchi):
            if chi_irr[l1, m1] == 1:
                irr_chi[0, it] = l1 + 1
                irr_chi[1, it] = m1 + 1
                it += 1
    return chi_map, irr_chi

def get_chi0_conv(Gk: np.ndarray, kmap: np.ndarray, invk: np.ndarray, olist: np.ndarray,
                  temp: float, Nx: int, Ny: int, Nz: int) -> np.ndarray:
    """
    @fn get_chi0_conv
    @brief Compute the irreducible susceptibility chi0 by convolution of Green's functions over the full BZ.
    Same as get_chi0 but without the interaction matrices (no Stoner-factor check/print).
    @param     Gk: Green's function tensor [Norb, Norb, Nw, Nk] complex128
    @param   kmap: Full k-grid FFT indices [Nkall, 3] int64
    @param   invk: Inverse k-point mapping [Nkall, 3] int64
    @param  olist: Orbital index pairs for chi [Nchi, 2] int64
    @param   temp: Temperature in eV
    @param     Nx: Number of k-points along kx
    @param     Ny: Number of k-points along ky
    @param     Nz: Number of k-points along kz
    @return   chi: Irreducible susceptibility tensor [Nchi, Nchi, Nw, Nk] complex128
    """
    Nkall, Nk = len(kmap), len(Gk[0, 0, 0])
    Norb, Nchi = len(Gk), len(olist)
    Nw = Gk.shape[2]
    chi = np.zeros((Nchi, Nchi, Nw, Nk), dtype=np.complex128)
    chi_map, irr_chi = _get_chi_map(olist)
    _lib.get_chi0_conv_(chi, Gk, kmap, invk, irr_chi, chi_map, olist,
                        dbl(temp), i64(Nx),
                        i64(Ny), i64(Nz), i64(Nw), i64(Nk),
                        i64(Nkall), i64(Norb), i64(Nchi))
    return chi

def get_chi_irr(uni: np.ndarray, eig: np.ndarray, ffermi: np.ndarray, qshift: np.ndarray,
                olist: np.ndarray, wlist: np.ndarray, idelta: float, temp: float) -> np.ndarray:
    """
    @fn get_chi_irr
    @brief Compute the irreducible (bare bubble) susceptibility chi0(q, w) for a given momentum shift.
    @param     uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param     eig: Eigenvalues [Nk, Norb] float64
    @param  ffermi: Fermi-Dirac occupations [Nk, Norb] float64
    @param  qshift: Index of k+q for each k-point [Nk] int64 (from get_qshift)
    @param   olist: Orbital index list [Nchi] int64
    @param   wlist: Real-frequency mesh [Nw] float64
    @param  idelta: Lorentzian broadening parameter in eV
    @param    temp: Temperature in eV
    @return   chi0: Irreducible susceptibility [Nw, Nchi, Nchi] complex128
    """
    Nk, Nw = len(eig), len(wlist)
    Norb, Nchi = eig.shape[1], len(olist)
    chi = np.zeros((Nw, Nchi, Nchi), dtype=np.complex128)
    eps = idelta * 1e-3
    _lib.get_chi_irr(chi, uni, eig, ffermi, qshift, olist, wlist, i64(Nchi),
                     i64(Norb), i64(Nk), i64(Nw),
                     dbl(idelta), dbl(eps), dbl(temp))
    return chi

def chis_qmap(uni: np.ndarray, eig: np.ndarray, ffermi: np.ndarray, klist: np.ndarray,
              Smat: np.ndarray, olist: np.ndarray, Nx: int, Ny: int, temp: float,
              ecut: float, idelta: float) -> tuple[np.ndarray, np.ndarray]:
    """
    @fn chis_qmap
    @brief Compute spin susceptibility and bare susceptibility maps on the qx-qy plane.
    @param     uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param     eig: Eigenvalues [Nk, Norb] float64
    @param  ffermi: Fermi-Dirac occupations [Nk, Norb] float64
    @param   klist: k-point list [Nk, 3] float64
    @param    Smat: Spin interaction matrix [Nchi, Nchi] float64
    @param   olist: Orbital index list [Nchi] int64
    @param      Nx: Number of q-points along qx
    @param      Ny: Number of q-points along qy
    @param    temp: Temperature in eV
    @param    ecut: Energy cutoff for the BZ summation in eV
    @param  idelta: Lorentzian broadening parameter in eV
    @retval   chis: Spin susceptibility map [Nx, Ny] complex128
    @retval    chi: Irreducible susceptibility map [Nx, Ny] complex128
    """
    Nk = len(eig)
    Norb, Nchi = eig.shape[1], len(olist)
    chi = np.zeros((Nx, Ny), dtype=np.complex128)
    chis = np.zeros((Nx, Ny), dtype=np.complex128)
    eps = idelta * 1e-3
    _lib.chiq_map(chis, chi, uni, eig, ffermi, klist, Smat, olist, dbl(temp),
                  dbl(ecut), dbl(idelta), dbl(eps),
                  i64(Nx), i64(Ny), i64(Nk),
                  i64(Norb), i64(Nchi))
    return chis, chi

def get_tr_chi(chis: np.ndarray, chi0: np.ndarray,
               olist: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    @fn get_tr_chi
    @brief Compute orbital-resolved traces of the interacting and bare susceptibilities.
    @param    chis: Interacting spin susceptibility [Nw, Nchi, Nchi] complex128
    @param    chi0: Irreducible susceptibility [Nw, Nchi, Nchi] complex128
    @param   olist: Orbital index list [Nchi] int64
    @retval trchis: Trace of interacting chi_s [Nw] complex128
    @retval trchi0: Trace of bare chi0 [Nw] complex128
    @retval chis_orb: Orbital-indexed chi_s [Nw, Norb+2] complex128
    """
    Nchi, Nw, Norb = len(olist), len(chi0), olist.max()
    trchis = np.zeros(Nw, dtype=np.complex128)
    trchi0 = np.zeros(Nw, dtype=np.complex128)
    chis_orb = np.zeros((Nw, Norb + 2), dtype=np.complex128)
    _lib.get_tr_chi(trchis, trchi0, chis_orb, chis, chi0, olist,
                    i64(Nw), i64(Nchi), i64(Norb))
    return trchis, trchi0, chis_orb

def get_phi_irr(uni: np.ndarray, eig: np.ndarray, ffermi: np.ndarray, qshift: np.ndarray,
                olist: np.ndarray, wlist: np.ndarray, idelta: float, mu: float,
                temp: float) -> np.ndarray:
    """
    @fn get_phi_irr
    @brief Compute the irreducible pairing susceptibility phi0(q, w) in the particle-particle channel.
    @param     uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param     eig: Eigenvalues [Nk, Norb] float64
    @param  ffermi: Fermi-Dirac occupations [Nk, Norb] float64
    @param  qshift: Index of k+q for each k-point [Nk] int64 (from get_qshift)
    @param   olist: Orbital index list [Nchi] int64
    @param   wlist: Frequency mesh [Nw] float64
    @param  idelta: Lorentzian broadening parameter in eV
    @param      mu: Chemical potential in eV
    @param    temp: Temperature in eV
    @return   phi0: Irreducible pairing susceptibility [Nw, Nchi, Nchi] complex128
    """
    Nk, Nw = len(eig), len(wlist)
    Norb, Nchi = eig.shape[1], len(olist)
    phi = np.zeros((Nw, Nchi, Nchi), dtype=np.complex128)
    eps = idelta * 1e-3
    _lib.get_phi_irr(phi, uni, eig, ffermi, qshift, olist, wlist, i64(Nchi),
                     i64(Norb), i64(Nk), i64(Nw),
                     dbl(idelta), dbl(eps),
                     dbl(mu), dbl(temp))
    return phi

def phi_qmap(uni: np.ndarray, eig: np.ndarray, ffermi: np.ndarray, klist: np.ndarray,
             olist: np.ndarray, Nx: int, Ny: int, mu: float, temp: float, ecut: float,
             idelta: float, sw_omega: bool) -> np.ndarray:
    """
    @fn phi_qmap
    @brief Compute the pairing susceptibility map on the qx-qy plane at a given energy cutoff.
    @param      uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param      eig: Eigenvalues [Nk, Norb] float64
    @param   ffermi: Fermi-Dirac occupations [Nk, Norb] float64
    @param    klist: k-point list [Nk, 3] float64
    @param    olist: Orbital index list [Nchi] int64
    @param       Nx: Number of q-points along qx
    @param       Ny: Number of q-points along qy
    @param       mu: Chemical potential in eV
    @param     temp: Temperature in eV
    @param     ecut: Energy cutoff for the BZ summation in eV
    @param   idelta: Lorentzian broadening parameter in eV
    @param sw_omega: If True, include Matsubara frequency summation; if False, use static limit
    @return     phi: Pairing susceptibility map [Nx, Ny] complex128
    """
    Nk = len(eig)
    Norb, Nchi = eig.shape[1], len(olist)
    phi = np.zeros((Nx, Ny), dtype=np.complex128)
    eps = idelta * 1e-3
    _lib.phiq_map(phi, uni, eig, ffermi, klist, olist, dbl(mu), dbl(temp),
                  dbl(ecut), dbl(idelta), dbl(eps),
                  i64(Nx), i64(Ny), i64(Nk),
                  i64(Norb), i64(Nchi), byref(c_bool(sw_omega)))
    return phi

def get_tr_phi(phi: np.ndarray, olist: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    @fn get_tr_phi
    @brief Compute orbital-resolved trace of the pairing susceptibility phi.
    @param     phi: Pairing susceptibility [Nw, Nchi, Nchi] complex128
    @param   olist: Orbital index list [Nchi] int64
    @retval  trphi: Trace of phi [Nw] complex128
    @retval phi_orb: Orbital-indexed phi [Nw, Norb+2] complex128
    """
    Nchi, Nw, Norb = len(olist), len(phi), olist.max()
    trphi = np.zeros(Nw, dtype=np.complex128)
    phi_orb = np.zeros((Nw, Norb + 2), dtype=np.complex128)
    _lib.get_tr_phi(trphi, phi_orb, phi, olist, i64(Nw),
                    i64(Nchi), i64(Norb))
    return trphi, phi_orb

def get_chis(chi0: np.ndarray, Smat: np.ndarray) -> np.ndarray:
    """
    @fn get_chis
    @brief Compute the RPA-corrected spin susceptibility chi_s from chi0 and the Stoner matrix.
    @param  chi0: Irreducible susceptibility [Nw, Nchi, Nchi] complex128
    @param  Smat: Spin interaction matrix [Nchi, Nchi] float64
    @return chis: Interacting spin susceptibility [Nw, Nchi, Nchi] complex128
    """
    Nchi, Nw = len(Smat), len(chi0)
    chis = np.zeros((Nw, Nchi, Nchi), dtype=np.complex128)
    _lib.get_chis(chis, chi0, Smat, i64(Nchi), i64(Nw))
    return chis

def get_chi0(Smat: np.ndarray, Cmat: np.ndarray, Gk: np.ndarray, olist: np.ndarray,
             kmap: np.ndarray, invk: np.ndarray, temp: float, Nx: int, Ny: int, Nz: int) -> tuple[np.ndarray, float]:
    """
    @fn get_chi0
    @brief Compute the irreducible susceptibility chi0 on the irreducible q-grid using the full BZ Green's function.
    @param   Smat: Spin interaction matrix [Nchi, Nchi] float64
    @param   Cmat: Charge interaction matrix [Nchi, Nchi] float64
    @param     Gk: Green's function on full k-grid [Norb, Norb, Nw, Nkall] complex128
    @param  olist: Orbital index list [Nchi] int64
    @param   kmap: Full-to-reduced k-point mapping [Nkall] int64
    @param   invk: Inverse k-point mapping [Nkall] int64
    @param   temp: Temperature in eV
    @param     Nx: Number of k-points along kx
    @param     Ny: Number of k-points along ky
    @param     Nz: Number of k-points along kz
    @return   chi: Irreducible susceptibility [Nchi, Nchi, Nw, Nk] complex128
    @return stoner: Maximum Stoner factor max(eig(chi0(q,0)*S)) over q (SDW instability when >= 1)
    """
    Norb, Nchi = len(Gk), len(olist)
    Nkall, Nk, Nw = len(kmap), len(Gk[0, 0, 0]), len(Gk[0, 0])
    chi = np.zeros((Nchi, Nchi, Nw, Nk), dtype=np.complex128)
    stoner = c_double(0.0)
    _lib.get_chi0(chi, Smat, Cmat, Gk, kmap, invk, olist, dbl(temp),
                  i64(Nx), i64(Ny), i64(Nz),
                  i64(Nw), i64(Nk), i64(Nkall),
                  i64(Norb), i64(Nchi), byref(stoner))
    return chi, stoner.value

def get_chi0_soc(Vmat: np.ndarray, Gk: np.ndarray, olist: np.ndarray, slist: np.ndarray,
                 kmap: np.ndarray, invk: np.ndarray, invs: np.ndarray, temp: float,
                 Nx: int, Ny: int, Nz: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    @fn get_chi0_soc
    @brief Compute the SOC-aware irreducible susceptibility chi0 and auxiliary spin-sign matrices.
    @param   Vmat: SOC interaction matrix [Nchi, Nchi] float64
    @param     Gk: Green's function on full k-grid [Norb, Norb, Nw, Nkall] complex128
    @param  olist: Orbital index list [Nchi] int64
    @param  slist: Spin-orbit label for each orbital [Norb] int64
    @param   kmap: Full-to-reduced k-point mapping [Nkall] int64
    @param   invk: Inverse k-point mapping [Nkall] int64
    @param   invs: Inverse spin index mapping [Norb] int64
    @param   temp: Temperature in eV
    @param     Nx: Number of k-points along kx
    @param     Ny: Number of k-points along ky
    @param     Nz: Number of k-points along kz
    @retval    chi: Irreducible susceptibility [Nchi, Nchi, Nw, Nk] complex128
    @retval sgnsig: Sign matrix [Norb, Norb] float64
    @retval sgnsig2: Sign matrix for chi orbital mapping [Nchi, Nchi] float64
    @retval invschi: Inverse chi orbital mapping [Nchi] int64
    """
    Norb, Nchi = len(slist), len(olist)
    Nkall, Nk, Nw = len(kmap), len(Gk[0, 0, 0]), len(Gk[0, 0])
    chi = np.zeros((Nchi, Nchi, Nw, Nk), dtype=np.complex128)
    sgnsig = np.zeros((Norb, Norb), dtype=np.float64)
    sgnsig2 = np.zeros((Nchi, Nchi), dtype=np.float64)
    invschi = np.zeros(Nchi, dtype=np.int64)
    _lib.get_chi0_soc(chi, sgnsig, sgnsig2, invschi, Vmat, Gk, kmap, invk, invs, olist, slist,
                      dbl(temp), i64(Nx), i64(Ny),
                      i64(Nz), i64(Nw), i64(Nk),
                      i64(Nkall), i64(Nchi), i64(Norb))
    return chi, sgnsig, sgnsig2, invschi

def get_chis_chic_soc(chi: np.ndarray, Vmat: np.ndarray, olist: np.ndarray, slist: np.ndarray,
                      invs: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    @fn get_chis_chic_soc
    @brief Decompose the SOC susceptibility into charge (chic), longitudinal spin (chiszz), and transverse (chispm) channels.
    @param    chi: Full susceptibility [Nchi, Nchi, Nk, Nw] complex128
    @param   Vmat: SOC interaction matrix [Nchi, Nchi] float64
    @param  olist: Orbital index list [Nchi] int64
    @param  slist: Spin-orbit label for each orbital [Norb] int64
    @param   invs: Inverse spin index mapping [Norb] int64
    @retval  chic: Charge susceptibility [Nchi/4, Nchi/4, Nk] complex128
    @retval chiszz: Longitudinal spin susceptibility [Nchi/4, Nchi/4, Nk] complex128
    @retval chispm: Transverse spin susceptibility [Nchi/4, Nchi/4, Nk] complex128
    """
    def get_orb_list(olist, slist, invs):
        orb_list = np.zeros(Nchi, dtype=np.int64)
        orbs = np.where(slist == 1)
        o1, o2 = np.meshgrid(orbs, orbs)
        orb_chi = np.array([o1.flatten(), o2.flatten()]).T
        for i, ol in enumerate(olist):
            o1_val = (ol[0] if slist[ol[0] - 1] == 1 else invs[ol[0] - 1]) - 1
            o2_val = (ol[1] if slist[ol[1] - 1] == 1 else invs[ol[1] - 1]) - 1
            for j, ob2 in enumerate(orb_chi):
                if o1_val == ob2[0] and o2_val == ob2[1]:
                    orb_list[i] = j + 1
                    break
        return orb_list
    Norb, Nchi = len(slist), len(olist)
    Nk, Nw = len(chi[0, 0, 0]), len(chi[0, 0])
    orb_list = get_orb_list(olist, slist, invs)
    chic = np.zeros((int(Nchi / 4), int(Nchi / 4), Nk), dtype=np.complex128)
    chiszz = np.zeros((int(Nchi / 4), int(Nchi / 4), Nk), dtype=np.complex128)
    chispm = np.zeros((int(Nchi / 4), int(Nchi / 4), Nk), dtype=np.complex128)
    _lib.get_chis_chic_soc(chic, chiszz, chispm, chi, Vmat, orb_list, olist, slist, invs,
                           i64(Nk), i64(Nw),
                           i64(Nchi), i64(Norb))
    return chic, chiszz, chispm

def get_chis_chic(chi: np.ndarray, Smat: np.ndarray, Cmat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    @fn get_chis_chic
    @brief Compute RPA-corrected spin susceptibility chis and charge susceptibility chic from chi0.
    @param   chi: Irreducible susceptibility [Nchi, Nchi, Nk, Nw] complex128
    @param  Smat: Spin interaction matrix [Nchi, Nchi] float64
    @param  Cmat: Charge interaction matrix [Nchi, Nchi] float64
    @retval chis: Interacting spin susceptibility [Nchi, Nchi, Nk] complex128
    @retval chic: Interacting charge susceptibility [Nchi, Nchi, Nk] complex128
    """
    Nk, Nw, Nchi = len(chi[0, 0, 0]), len(chi[0, 0]), len(Smat)
    chis = np.zeros((Nchi, Nchi, Nk), dtype=np.complex128)
    chic = np.zeros((Nchi, Nchi, Nk), dtype=np.complex128)
    _lib.get_chis_chic(chis, chic, chi, Smat, Cmat, i64(Nk),
                       i64(Nw), i64(Nchi))
    return chis, chic

def get_chi_irr_sc(uni: np.ndarray, eig: np.ndarray, ffermi: np.ndarray,
                   qshift: np.ndarray, ol: np.ndarray, wlist: np.ndarray,
                   idelta: float, temp: float, sw_spsym: bool) -> np.ndarray:
    """
    @fn get_chi_irr_sc
    @brief Compute the irreducible susceptibility in the superconducting state using BdG eigenstates.
    @param      uni: BdG eigenvector matrices [Nk, 2*Norb, 2*Norb] complex128
    @param      eig: BdG eigenvalues [Nk, 2*Norb] float64
    @param   ffermi: Fermi-Dirac occupations for BdG bands [Nk, 2*Norb] float64
    @param   qshift: Index of k+q for each k-point [Nk] int64 (from get_qshift)
    @param       ol: Orbital index pairs for chi [Nchi, 2] int64 (Fortran-contiguous)
    @param    wlist: Real-frequency mesh [Nw] float64
    @param   idelta: Lorentzian broadening parameter in eV
    @param     temp: Temperature in eV
    @param sw_spsym: True for triplet (dz) symmetry, False for singlet
    @return     chi: Irreducible susceptibility in SC state [Nw, Nchi, Nchi] complex128
    """
    Nk, Nw = len(eig), len(wlist)
    Norb = eig.shape[1] // 2
    Nchi = len(ol)
    eps = idelta * 1e-3
    chi = np.zeros((Nw, Nchi, Nchi,2), dtype=np.complex128)
    ol_f = np.asfortranarray(ol)
    _lib.get_chi_irr_sc(chi, uni, eig, ffermi, qshift, ol_f, wlist,
                        i64(Nchi), i64(Norb),
                        i64(Nk), i64(Nw),
                        dbl(idelta), dbl(eps),
                        dbl(temp), byref(c_bool(sw_spsym)))
    return chi

def get_eig_or_tr_chi(chi: np.ndarray, invk: np.ndarray, sw_eig: bool) -> np.ndarray:
    """
    @fn get_eig_or_tr_chi
    @brief Compute either the leading eigenvalue or the trace of the susceptibility on the full k-grid.
    @param    chi: Susceptibility [Nchi, Nchi, Nk] complex128
    @param   invk: Inverse k-point mapping [Nkall] int64
    @param sw_eig: If True, return leading eigenvalue; if False, return trace
    @return  chiq: Eigenvalue or trace on the full k-grid [Nkall] complex128
    """
    Nkall, Nk = len(invk), len(chi.T)
    Nchi = chi.shape[0]
    chiq = np.zeros(Nkall, dtype=np.complex128)
    _lib.get_eig_or_tr_chi(chiq, chi, invk, i64(Nkall), i64(Nk),
                            i64(Nchi), byref(c_bool(sw_eig)))
    return chiq
