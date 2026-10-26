import ctypes
import numpy as np
import libs.flibs._loader as _loader

# Public API of this submodule; keeps `from ._x import *` in __init__.py from
# re-exporting numpy, scipy, ctypes and the sibling packages imported above.
__all__ = [
    'get_chi0_conv', 'get_chi_irr', 'chiq_map', 'get_tr_chi', 'get_phi_irr', 'phiq_map', 'get_tr_phi', 'get_chis',
    'get_chi0', 'get_chi0_tail', 'get_chi0_soc', 'get_chis_chic_soc', 'get_chis_chic', 'get_chi_irr_sc',
    'get_eig_or_tr_chi', 'get_chi0_nosym'
]

# --- ctypes signatures: set once at import.
# All Fortran entry points are subroutines, so restype is always None.
_loader._lib.get_chi0_conv_.argtypes = [
    _loader.c_array(np.complex128),  # chi
    _loader.c_array(np.complex128),  # Gk
    _loader.c_array(np.int64),       # kmap
    _loader.c_array(np.int64),       # invk
    _loader.c_array(np.int32),       # irr_chi
    _loader.c_array(np.int32),       # chi_map
    _loader.f_array(np.int64),       # olist
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int64),          # temp, Nx
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Ny, Nz
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nw, Nk
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nkall, Norb
    ctypes.POINTER(ctypes.c_int64)                              # Nchi
]
_loader._lib.get_chi0_conv_.restype = None
_loader._lib.get_chi_irr.argtypes = [
    _loader.c_array(np.complex128),  # chi (out)
    _loader.c_array(np.complex128),  # uni
    _loader.c_array(np.float64),     # eig
    _loader.c_array(np.float64),     # ffermi
    _loader.c_array(np.int64),       # qshift
    _loader.f_array(np.int64),       # ol
    _loader.c_array(np.float64),     # wl
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nchi, Norb
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nk, Nw
    ctypes.POINTER(ctypes.c_double),                            # idelta
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),         # eps, temp
]
_loader._lib.get_chi_irr.restype = None
_loader._lib.chiq_map.argtypes = [
    _loader.c_array(np.complex128),          # trchis (out)
    _loader.c_array(np.complex128),          # trchi (out)
    _loader.c_array(np.complex128),          # uni
    _loader.c_array(np.float64),             # eig
    _loader.c_array(np.float64),             # ffermi
    _loader.c_array(np.float64),             # klist
    _loader.c_array(np.float64),             # Smat
    _loader.f_array(np.int64),               # ol
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),                 # temp, ecut
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),                 # idelta, eps
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),                   # Nx, Ny
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)  # Nk, Norb, Nchi
]
_loader._lib.chiq_map.restype = None
_loader._lib.get_tr_chi.argtypes = [
    _loader.c_array(np.complex128),  # trchis (out)
    _loader.c_array(np.complex128),  # trchi0 (out)
    _loader.c_array(np.complex128),  # chis_orb (out)
    _loader.c_array(np.complex128),  # chis
    _loader.c_array(np.complex128),  # chi0
    _loader.f_array(np.int64),       # olist
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nw, Nchi
    ctypes.POINTER(ctypes.c_int64)                              # Norb
]
_loader._lib.get_tr_chi.restype = None
_loader._lib.get_phi_irr.argtypes = [
    _loader.c_array(np.complex128),  # phi (out)
    _loader.c_array(np.complex128),  # uni
    _loader.c_array(np.float64),     # eig
    _loader.c_array(np.float64),     # ffermi
    _loader.c_array(np.int64),       # qshift
    _loader.f_array(np.int64),       # ol
    _loader.c_array(np.float64),     # wl
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nchi, Norb
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nk, Nw
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),         # idelta, eps
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),         # mu, temp
]
_loader._lib.get_phi_irr.restype = None
_loader._lib.phiq_map.argtypes = [
    _loader.c_array(np.complex128),  # trphi (out)
    _loader.c_array(np.complex128),  # uni
    _loader.c_array(np.float64),     # eig
    _loader.c_array(np.float64),     # ffermi
    _loader.c_array(np.float64),     # klist
    _loader.f_array(np.int64),       # ol
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),         # mu, temp
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),         # ecut, idelta
    ctypes.POINTER(ctypes.c_double),                            # eps
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nx, Ny
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nk, Norb
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_bool)             # Nchi, sw_omega
]
_loader._lib.phiq_map.restype = None
_loader._lib.get_tr_phi.argtypes = [
    _loader.c_array(np.complex128),  # trphi (out)
    _loader.c_array(np.complex128),  # phi_orb (out)
    _loader.c_array(np.complex128),  # phi
    _loader.f_array(np.int64),       # olist
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nw, Nchi
    ctypes.POINTER(ctypes.c_int64)                              # Norb
]
_loader._lib.get_tr_phi.restype = None
_loader._lib.get_chis.argtypes = [
    _loader.c_array(np.complex128),  # chis (out)
    _loader.c_array(np.complex128),  # chi0
    _loader.c_array(np.float64),     # Smat
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nchi, Nw
]
_loader._lib.get_chis.restype = None
_loader._lib.get_chi0.argtypes = [
    _loader.c_array(np.complex128),  # chi (out)
    _loader.c_array(np.float64),     # Smat
    _loader.c_array(np.float64),     # Cmat
    _loader.c_array(np.complex128),  # Gk
    _loader.c_array(np.int64),       # kmap
    _loader.c_array(np.int64),       # invk
    _loader.f_array(np.int64),       # olist
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int64),          # temp, Nx
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Ny, Nz
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nw, Nk
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nkall, Norb
    ctypes.POINTER(ctypes.c_int64),                             # Nchi
    ctypes.POINTER(ctypes.c_double)                             # maxchi0s_out (Stoner factor)
]
_loader._lib.get_chi0.restype = None
_loader._lib.get_chi0_tail.argtypes = [
    _loader.c_array(np.complex128),  # chi
    _loader.c_array(np.float64),     # Smat
    _loader.c_array(np.float64),     # Cmat
    _loader.c_array(np.complex128),  # Gk
    _loader.c_array(np.float64),     # eig
    _loader.c_array(np.complex128),  # uni
    _loader.c_array(np.int64),       # kmap
    _loader.c_array(np.int64),       # invk
    _loader.f_array(np.int64),       # olist
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),         # mu, temp
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nx, Ny
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nz, Nw
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nk, Nkall
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Norb, Nchi
    ctypes.POINTER(ctypes.c_double)                             # maxchi0s_out (Stoner factor)
]
_loader._lib.get_chi0_tail.restype = None
_loader._lib.get_chi0_soc.argtypes = [
    _loader.c_array(np.complex128),  # chi (out)
    _loader.c_array(np.float64),     # sgnsig (out)
    _loader.c_array(np.float64),     # sgnsig2 (out)
    _loader.c_array(np.int64),       # invschi (out)
    _loader.c_array(np.float64),     # Vmat
    _loader.c_array(np.complex128),  # Gk
    _loader.c_array(np.int64),       # kmap
    _loader.c_array(np.int64),       # invk
    _loader.c_array(np.int64),       # invs
    _loader.f_array(np.int64),       # olist
    _loader.c_array(np.int64),       # slist
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int64),          # temp, Nx
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Ny, Nz
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nw, Nk
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nkall, Nchi
    ctypes.POINTER(ctypes.c_int64)                              # Norb
]
_loader._lib.get_chi0_soc.restype = None
_loader._lib.get_chi0_nosym.argtypes = [
    _loader.c_array(np.complex128),  # chi (out)
    _loader.c_array(np.complex128),  # Gk
    _loader.c_array(np.int64),       # kmap
    _loader.f_array(np.int64),       # olist
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int64),          # temp, Nx
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Ny, Nz
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nw, Nkall
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)            # Norb, Nchi
]
_loader._lib.get_chi0_nosym.restype = None
_loader._lib.get_chis_chic_soc.argtypes = [
    _loader.c_array(np.complex128),  # chic (out)
    _loader.c_array(np.complex128),  # chiszz (out)
    _loader.c_array(np.complex128),  # chispm (out)
    _loader.c_array(np.complex128),  # chisxx (out)
    _loader.c_array(np.complex128),  # chisyy (out)
    _loader.c_array(np.complex128),  # chi
    _loader.c_array(np.float64),     # Vmat
    _loader.c_array(np.int64),       # orb_list
    _loader.f_array(np.int64),       # olist
    _loader.c_array(np.int64),       # slist
    _loader.c_array(np.int64),       # invs
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nk, Nw
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)            # Nchi, Norb
]
_loader._lib.get_chis_chic_soc.restype = None
_loader._lib.get_chis_chic.argtypes = [
    _loader.c_array(np.complex128),          # chis (out)
    _loader.c_array(np.complex128),          # chic (out)
    _loader.c_array(np.complex128),          # chi
    _loader.c_array(np.float64),             # Smat
    _loader.c_array(np.float64),             # Cmat
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)  # Nk, Nw, Nchi
]
_loader._lib.get_chis_chic.restype = None
_loader._lib.get_chi_irr_sc.argtypes = [
    _loader.c_array(np.complex128),  # chi (out)
    _loader.c_array(np.complex128),  # uni
    _loader.c_array(np.float64),     # eig
    _loader.c_array(np.float64),     # ffermi
    _loader.c_array(np.int64),       # qshift
    _loader.f_array(np.int64),       # ol
    _loader.c_array(np.float64),     # wl
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nchi, Norb
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nk, Nw
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),         # idelta, eps
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_bool),           # temp, sw_spsym
]
_loader._lib.get_chi_irr_sc.restype = None
_loader._lib.get_eig_or_tr_chi.argtypes = [
    _loader.c_array(np.complex128),  # chiqout (out)
    _loader.c_array(np.complex128),  # chi
    _loader.c_array(np.int64),       # invk
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nkall, Nk
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_bool)             # Nchi, sw_eig
]
_loader._lib.get_eig_or_tr_chi.restype = None


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
    _loader._lib.get_chi0_conv_(chi, Gk, kmap, invk, irr_chi, chi_map, olist,
                                _loader.dbl(temp), _loader.i64(Nx),
                                _loader.i64(Ny), _loader.i64(Nz), _loader.i64(Nw), _loader.i64(Nk),
                                _loader.i64(Nkall), _loader.i64(Norb), _loader.i64(Nchi))
    return chi


def get_chi_irr(uni: np.ndarray, eig: np.ndarray, ffermi: np.ndarray, qshift: np.ndarray,
                olist: np.ndarray, wlist: np.ndarray, idelta: float, temp: float) -> np.ndarray:
    """
    @brief Compute the irreducible (bare bubble) susceptibility chi0(q, w) for a given momentum shift.
    @param     uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param     eig: Eigenvalues [Nk, Norb] float64
    @param  ffermi: Fermi-Dirac occupations [Nk, Norb] float64
    @param  qshift: Index of k+q for each k-point [Nk] int64 (from get_qshift)
    @param   olist: Orbital-pair list of chi [Nchi, 2] int64 (1-based)
    @param   wlist: Real-frequency mesh [Nw] float64
    @param  idelta: Lorentzian broadening parameter in eV
    @param    temp: Temperature in eV
    @return   chi0: Irreducible susceptibility [Nw, Nchi, Nchi] complex128
    """
    Nk, Nw = len(eig), len(wlist)
    Norb, Nchi = eig.shape[1], len(olist)
    chi = np.zeros((Nw, Nchi, Nchi), dtype=np.complex128)
    eps = idelta * 1e-3
    _loader._lib.get_chi_irr(chi, uni, eig, ffermi, qshift, olist, wlist, _loader.i64(Nchi),
                             _loader.i64(Norb), _loader.i64(Nk), _loader.i64(Nw),
                             _loader.dbl(idelta), _loader.dbl(eps), _loader.dbl(temp))
    return chi


def chiq_map(uni: np.ndarray, eig: np.ndarray, ffermi: np.ndarray, klist: np.ndarray,
              Smat: np.ndarray, olist: np.ndarray, Nx: int, Ny: int, temp: float,
              ecut: float, idelta: float) -> tuple[np.ndarray, np.ndarray]:
    """
    @brief Compute spin susceptibility and bare susceptibility maps on the qx-qy plane.
    @param     uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param     eig: Eigenvalues [Nk, Norb] float64
    @param  ffermi: Fermi-Dirac occupations [Nk, Norb] float64
    @param   klist: k-point list [Nk, 3] float64
    @param    Smat: Spin interaction matrix [Nchi, Nchi] float64
    @param   olist: Orbital-pair list of chi [Nchi, 2] int64 (1-based)
    @param      Nx: Number of q-points along qx
    @param      Ny: Number of q-points along qy
    @param    temp: Temperature in eV
    @param    ecut: Energy cutoff for the BZ summation in eV
    @param  idelta: Lorentzian broadening parameter in eV
    @return   chis: Spin susceptibility map [Nx, Ny] complex128
    @return    chi: Irreducible susceptibility map [Nx, Ny] complex128
    """
    Nk = len(eig)
    Norb, Nchi = eig.shape[1], len(olist)
    chi = np.zeros((Nx, Ny), dtype=np.complex128)
    chis = np.zeros((Nx, Ny), dtype=np.complex128)
    eps = idelta * 1e-3
    _loader._lib.chiq_map(chis, chi, uni, eig, ffermi, klist, Smat, olist, _loader.dbl(temp),
                          _loader.dbl(ecut), _loader.dbl(idelta), _loader.dbl(eps),
                          _loader.i64(Nx), _loader.i64(Ny), _loader.i64(Nk),
                          _loader.i64(Norb), _loader.i64(Nchi))
    return chis, chi


def get_tr_chi(chis: np.ndarray, chi0: np.ndarray,
               olist: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    @brief Compute orbital-resolved traces of the interacting and bare susceptibilities.
    @param    chis: Interacting spin susceptibility [Nw, Nchi, Nchi] complex128
    @param    chi0: Irreducible susceptibility [Nw, Nchi, Nchi] complex128
    @param   olist: Orbital-pair list of chi [Nchi, 2] int64 (1-based)
    @return trchis: Trace of interacting chi_s [Nw] complex128
    @return trchi0: Trace of bare chi0 [Nw] complex128
    @return chis_orb: Orbital-indexed chi_s [Nw, Norb+2] complex128
    """
    Nchi, Nw, Norb = len(olist), len(chi0), olist.max()
    trchis = np.zeros(Nw, dtype=np.complex128)
    trchi0 = np.zeros(Nw, dtype=np.complex128)
    chis_orb = np.zeros((Nw, Norb + 2), dtype=np.complex128)
    _loader._lib.get_tr_chi(trchis, trchi0, chis_orb, chis, chi0, olist,
                            _loader.i64(Nw), _loader.i64(Nchi), _loader.i64(Norb))
    return trchis, trchi0, chis_orb


def get_phi_irr(uni: np.ndarray, eig: np.ndarray, ffermi: np.ndarray, qshift: np.ndarray,
                olist: np.ndarray, wlist: np.ndarray, idelta: float, mu: float,
                temp: float) -> np.ndarray:
    """
    @brief Compute the irreducible pairing susceptibility phi0(q, w) in the particle-particle channel.
    @param     uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param     eig: Eigenvalues [Nk, Norb] float64
    @param  ffermi: Fermi-Dirac occupations [Nk, Norb] float64
    @param  qshift: Index of k+q for each k-point [Nk] int64 (from get_qshift)
    @param   olist: Orbital-pair list of phi [Nchi, 2] int64 (1-based)
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
    _loader._lib.get_phi_irr(phi, uni, eig, ffermi, qshift, olist, wlist, _loader.i64(Nchi),
                             _loader.i64(Norb), _loader.i64(Nk), _loader.i64(Nw),
                             _loader.dbl(idelta), _loader.dbl(eps),
                             _loader.dbl(mu), _loader.dbl(temp))
    return phi


def phiq_map(uni: np.ndarray, eig: np.ndarray, ffermi: np.ndarray, klist: np.ndarray,
             olist: np.ndarray, Nx: int, Ny: int, mu: float, temp: float, ecut: float,
             idelta: float, sw_omega: bool) -> np.ndarray:
    """
    @brief Compute the pairing susceptibility map on the qx-qy plane at a given energy cutoff.
    @param      uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param      eig: Eigenvalues [Nk, Norb] float64
    @param   ffermi: Fermi-Dirac occupations [Nk, Norb] float64
    @param    klist: k-point list [Nk, 3] float64
    @param    olist: Orbital-pair list of phi [Nchi, 2] int64 (1-based)
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
    _loader._lib.phiq_map(phi, uni, eig, ffermi, klist, olist, _loader.dbl(mu), _loader.dbl(temp),
                          _loader.dbl(ecut), _loader.dbl(idelta), _loader.dbl(eps),
                          _loader.i64(Nx), _loader.i64(Ny), _loader.i64(Nk),
                          _loader.i64(Norb), _loader.i64(Nchi), ctypes.byref(ctypes.c_bool(sw_omega)))
    return phi


def get_tr_phi(phi: np.ndarray, olist: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    @brief Compute orbital-resolved trace of the pairing susceptibility phi.
    @param     phi: Pairing susceptibility [Nw, Nchi, Nchi] complex128
    @param   olist: Orbital-pair list of phi [Nchi, 2] int64 (1-based)
    @return  trphi: Trace of phi [Nw] complex128
    @return phi_orb: Orbital-indexed phi [Nw, Norb+2] complex128
    """
    Nchi, Nw, Norb = len(olist), len(phi), olist.max()
    trphi = np.zeros(Nw, dtype=np.complex128)
    phi_orb = np.zeros((Nw, Norb + 2), dtype=np.complex128)
    _loader._lib.get_tr_phi(trphi, phi_orb, phi, olist, _loader.i64(Nw),
                            _loader.i64(Nchi), _loader.i64(Norb))
    return trphi, phi_orb


def get_chis(chi0: np.ndarray, Smat: np.ndarray) -> np.ndarray:
    """
    @brief Compute the RPA-corrected spin susceptibility chi_s from chi0 and the Stoner matrix.
    @param  chi0: Irreducible susceptibility [Nw, Nchi, Nchi] complex128
    @param  Smat: Spin interaction matrix [Nchi, Nchi] float64
    @return chis: Interacting spin susceptibility [Nw, Nchi, Nchi] complex128
    """
    Nchi, Nw = len(Smat), len(chi0)
    chis = np.zeros((Nw, Nchi, Nchi), dtype=np.complex128)
    _loader._lib.get_chis(chis, chi0, Smat, _loader.i64(Nchi), _loader.i64(Nw))
    return chis


def get_chi0(Smat: np.ndarray, Cmat: np.ndarray, Gk: np.ndarray, olist: np.ndarray,
             kmap: np.ndarray, invk: np.ndarray, temp: float, Nx: int, Ny: int, Nz: int) -> tuple[np.ndarray, float]:
    """
    @brief Compute the irreducible susceptibility chi0 on the irreducible q-grid using the full BZ Green's function.
    @param   Smat: Spin interaction matrix [Nchi, Nchi] float64
    @param   Cmat: Charge interaction matrix [Nchi, Nchi] float64
    @param     Gk: Green's function on full k-grid [Norb, Norb, Nw, Nkall] complex128
    @param  olist: Orbital-pair list of chi [Nchi, 2] int64 (1-based)
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
    stoner = ctypes.c_double(0.0)
    _loader._lib.get_chi0(chi, Smat, Cmat, Gk, kmap, invk, olist, _loader.dbl(temp),
                          _loader.i64(Nx), _loader.i64(Ny), _loader.i64(Nz),
                          _loader.i64(Nw), _loader.i64(Nk), _loader.i64(Nkall),
                          _loader.i64(Norb), _loader.i64(Nchi), ctypes.byref(stoner))
    return chi, stoner.value


def get_chi0_tail(Smat: np.ndarray, Cmat: np.ndarray, Gk: np.ndarray, eig: np.ndarray,
                  uni: np.ndarray, olist: np.ndarray, kmap: np.ndarray, invk: np.ndarray,
                  mu: float, temp: float, Nx: int, Ny: int, Nz: int) -> tuple[np.ndarray, float]:
    """
    @brief Tail-corrected variant of get_chi0: chi0 = conv[G] - conv[G0] + analytic
    reference bubble.  By bilinearity the sharp-cutoff FFT convolution only acts on
    the fast-decaying residual G-G0 (~1/w^2), while the slowly decaying 1/(iw)
    reference part is evaluated exactly in imaginary time, reducing the Matsubara
    truncation error from O(1/Nw) to O(1/Nw^2).  G0 is built internally from
    (eig, uni) at the SAME chemical potential mu as Gk.
    @param   Smat: Spin interaction matrix [Nchi, Nchi] float64
    @param   Cmat: Charge interaction matrix [Nchi, Nchi] float64
    @param     Gk: Green's function on the irreducible k-grid [Norb, Norb, Nw, Nk] complex128
    @param    eig: Eigenvalues [Nk, Norb] float64
    @param    uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param  olist: Orbital index pairs for chi [Nchi, 2] int64
    @param   kmap: Full-to-reduced k-point mapping [Nkall] int64
    @param   invk: Inverse k-point mapping [Nkall] int64
    @param     mu: Chemical potential of Gk in eV
    @param   temp: Temperature in eV
    @param Nx,Ny,Nz: k-grid dimensions
    @return   chi: Tail-corrected irreducible susceptibility [Nchi, Nchi, Nw, Nk] complex128
    @return stoner: Maximum Stoner factor max(eig(chi0(q,0)*S)) over q
    """
    Norb, Nchi = len(Gk), len(olist)
    Nkall, Nk, Nw = len(kmap), len(Gk[0, 0, 0]), len(Gk[0, 0])
    chi = np.zeros((Nchi, Nchi, Nw, Nk), dtype=np.complex128)
    stoner = ctypes.c_double(0.0)
    _loader._lib.get_chi0_tail(chi, Smat, Cmat, Gk, eig, uni, kmap, invk, olist,
                               _loader.dbl(mu), _loader.dbl(temp),
                               _loader.i64(Nx), _loader.i64(Ny), _loader.i64(Nz),
                               _loader.i64(Nw), _loader.i64(Nk), _loader.i64(Nkall),
                               _loader.i64(Norb), _loader.i64(Nchi), ctypes.byref(stoner))
    return chi, stoner.value


def get_chi0_soc(Vmat: np.ndarray, Gk: np.ndarray, olist: np.ndarray, slist: np.ndarray,
                 kmap: np.ndarray, invk: np.ndarray, invs: np.ndarray, temp: float,
                 Nx: int, Ny: int, Nz: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    @brief Compute the SOC-aware irreducible susceptibility chi0 and auxiliary spin-sign matrices.
    @param   Vmat: SOC interaction matrix [Nchi, Nchi] float64
    @param     Gk: Green's function on full k-grid [Norb, Norb, Nw, Nkall] complex128
    @param  olist: Orbital-pair list of chi [Nchi, 2] int64 (1-based)
    @param  slist: Spin-orbit label for each orbital [Norb] int64
    @param   kmap: Full-to-reduced k-point mapping [Nkall] int64
    @param   invk: Inverse k-point mapping [Nkall] int64
    @param   invs: Inverse spin index mapping [Norb] int64
    @param   temp: Temperature in eV
    @param     Nx: Number of k-points along kx
    @param     Ny: Number of k-points along ky
    @param     Nz: Number of k-points along kz
    @return    chi: Irreducible susceptibility [Nchi, Nchi, Nw, Nk] complex128
    @return sgnsig: Sign matrix [Norb, Norb] float64
    @return sgnsig2: Sign matrix for chi orbital mapping [Nchi, Nchi] float64
    @return invschi: Inverse chi orbital mapping [Nchi] int64
    """
    Norb, Nchi = len(slist), len(olist)
    Nkall, Nk, Nw = len(kmap), len(Gk[0, 0, 0]), len(Gk[0, 0])
    chi = np.zeros((Nchi, Nchi, Nw, Nk), dtype=np.complex128)
    sgnsig = np.zeros((Norb, Norb), dtype=np.float64)
    sgnsig2 = np.zeros((Nchi, Nchi), dtype=np.float64)
    invschi = np.zeros(Nchi, dtype=np.int64)
    _loader._lib.get_chi0_soc(chi, sgnsig, sgnsig2, invschi, Vmat, Gk, kmap, invk, invs, olist, slist,
                              _loader.dbl(temp), _loader.i64(Nx), _loader.i64(Ny),
                              _loader.i64(Nz), _loader.i64(Nw), _loader.i64(Nk),
                              _loader.i64(Nkall), _loader.i64(Nchi), _loader.i64(Norb))
    return chi, sgnsig, sgnsig2, invschi


def get_chi0_nosym(Gk: np.ndarray, olist: np.ndarray, kmap: np.ndarray, temp: float,
                   Nx: int, Ny: int, Nz: int) -> np.ndarray:
    """
    @brief Compute the irreducible susceptibility chi0 on the FULL BZ without any symmetry (no TRS/inversion).
    Only Hermiticity chi(q,-iv) = chi(q,iv)^+ is used, so it is valid with Zeeman fields / magnetic order.
    @param     Gk: Green's function on the full k-grid [Norb, Norb, Nw, Nkall] complex128 (spin-orbital basis)
    @param  olist: Orbital-pair list of chi [Nchi, 2] int64 (1-based)
    @param   kmap: Integer mesh coordinates of each k-point [Nkall, 3] int64
    @param   temp: Temperature in eV
    @param     Nx: Number of k-points along kx
    @param     Ny: Number of k-points along ky
    @param     Nz: Number of k-points along kz
    @return   chi: Irreducible susceptibility [Nchi, Nchi, Nw, Nkall] complex128 (v_n >= 0)
    """
    Norb, Nchi = len(Gk), len(olist)
    Nkall, Nw = len(kmap), len(Gk[0, 0])
    chi = np.zeros((Nchi, Nchi, Nw, Nkall), dtype=np.complex128)
    _loader._lib.get_chi0_nosym(chi, Gk, kmap, olist, _loader.dbl(temp),
                                _loader.i64(Nx), _loader.i64(Ny), _loader.i64(Nz),
                                _loader.i64(Nw), _loader.i64(Nkall), _loader.i64(Norb),
                                _loader.i64(Nchi))
    return chi


def get_chis_chic_soc(chi: np.ndarray, Vmat: np.ndarray, olist: np.ndarray, slist: np.ndarray,
                      invs: np.ndarray) -> tuple[np.ndarray, ...]:
    """
    @brief Decompose the static SOC RPA susceptibility into charge (chic), longitudinal
    spin (chiszz), transverse (chispm = <S^- S^+>) and Cartesian transverse (chisxx,
    chisyy) channels.  chiszz, chisxx and chisyy share the normalization of S = sigma/2,
    so with spin-rotation symmetry chisxx = chisyy = chiszz (and chispm = 2 chiszz);
    SOC splits them, and chisxx != chisyy needs <S^+S^+> != 0 (in-plane anisotropy).
    @param    chi: Full susceptibility [Nchi, Nchi, Nk, Nw] complex128
    @param   Vmat: SOC interaction matrix [Nchi, Nchi] float64
    @param  olist: Orbital-pair list of chi [Nchi, 2] int64 (1-based)
    @param  slist: Spin-orbit label for each orbital [Norb] int64
    @param   invs: Inverse spin index mapping [Norb] int64
    @return  chic: Charge susceptibility [Nchi/4, Nchi/4, Nk] complex128
    @return chiszz: Longitudinal spin susceptibility [Nchi/4, Nchi/4, Nk] complex128
    @return chispm: Transverse spin susceptibility [Nchi/4, Nchi/4, Nk] complex128
    @return chisxx: S^x S^x spin susceptibility [Nchi/4, Nchi/4, Nk] complex128
    @return chisyy: S^y S^y spin susceptibility [Nchi/4, Nchi/4, Nk] complex128
    """
    def get_orb_list(olist: np.ndarray, slist: np.ndarray, invs: np.ndarray) -> np.ndarray:
        # Spinless orbital-pair index of every chi index.  The orbitals are those chiolist
        # actually carries (up-spin representatives), not every up orbital of the model: a
        # chiolist restricted to the correlated orbitals (d-p, periodic Anderson) has Nchi/4 =
        # Ncorr^2 output pairs, and indexing them by all Norb/2 orbitals wrote past the end of the
        # [Nchi/4, Nchi/4] outputs whenever a correlated orbital was not among the first ones.
        up = lambda o: (o if slist[o - 1] == 1 else invs[o - 1]) - 1
        orbs = sorted({up(o) for ol in olist for o in ol})
        if len(orbs)**2 != Nchi // 4 or 4*(Nchi // 4) != Nchi:
            raise ValueError(f'get_chis_chic_soc: chiolist must hold every spin combination of its orbitals '
                             f'({len(orbs)} orbitals -> {4*len(orbs)**2} pairs), got Nchi = {Nchi}')
        o1, o2 = np.meshgrid(orbs, orbs)
        pos = {(a, b): j + 1 for j, (a, b) in enumerate(zip(o1.flatten(), o2.flatten()))}
        return np.array([pos[(up(ol[0]), up(ol[1]))] for ol in olist], dtype=np.int64)
    Norb, Nchi = len(slist), len(olist)
    Nk, Nw = len(chi[0, 0, 0]), len(chi[0, 0])
    orb_list = get_orb_list(olist, slist, invs)
    chic = np.zeros((int(Nchi / 4), int(Nchi / 4), Nk), dtype=np.complex128)
    chiszz = np.zeros((int(Nchi / 4), int(Nchi / 4), Nk), dtype=np.complex128)
    chispm = np.zeros((int(Nchi / 4), int(Nchi / 4), Nk), dtype=np.complex128)
    chisxx = np.zeros((int(Nchi / 4), int(Nchi / 4), Nk), dtype=np.complex128)
    chisyy = np.zeros((int(Nchi / 4), int(Nchi / 4), Nk), dtype=np.complex128)
    _loader._lib.get_chis_chic_soc(chic, chiszz, chispm, chisxx, chisyy, chi, Vmat, orb_list, olist, slist, invs,
                                   _loader.i64(Nk), _loader.i64(Nw),
                                   _loader.i64(Nchi), _loader.i64(Norb))
    return chic, chiszz, chispm, chisxx, chisyy


def get_chis_chic(chi: np.ndarray, Smat: np.ndarray, Cmat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    @brief Compute RPA-corrected spin susceptibility chis and charge susceptibility chic from chi0.
    @param   chi: Irreducible susceptibility [Nchi, Nchi, Nk, Nw] complex128
    @param  Smat: Spin interaction matrix [Nchi, Nchi] float64
    @param  Cmat: Charge interaction matrix [Nchi, Nchi] float64
    @return chis: Interacting spin susceptibility [Nchi, Nchi, Nk] complex128
    @return chic: Interacting charge susceptibility [Nchi, Nchi, Nk] complex128
    """
    Nk, Nw, Nchi = len(chi[0, 0, 0]), len(chi[0, 0]), len(Smat)
    chis = np.zeros((Nchi, Nchi, Nk), dtype=np.complex128)
    chic = np.zeros((Nchi, Nchi, Nk), dtype=np.complex128)
    _loader._lib.get_chis_chic(chis, chic, chi, Smat, Cmat, _loader.i64(Nk),
                               _loader.i64(Nw), _loader.i64(Nchi))
    return chis, chic


def get_chi_irr_sc(uni: np.ndarray, eig: np.ndarray, ffermi: np.ndarray,
                   qshift: np.ndarray, ol: np.ndarray, wlist: np.ndarray,
                   idelta: float, temp: float, sw_spsym: bool) -> np.ndarray:
    """
    @brief Compute the irreducible susceptibility in the superconducting state using BdG eigenstates.
    @param      uni: BdG eigenvector matrices [Nk, 2*Norb, 2*Norb] complex128
    @param      eig: BdG eigenvalues [Nk, 2*Norb] float64
    @param   ffermi: Fermi-Dirac occupations for BdG bands [Nk, 2*Norb] float64
    @param   qshift: Index of k+q for each k-point [Nk] int64 (from get_qshift)
    @param       ol: Orbital index pairs for chi [Nchi, 2] int64 (Fortran-contiguous)
    @param    wlist: Real-frequency mesh [Nw] float64
    @param   idelta: Lorentzian broadening parameter in eV
    @param     temp: Temperature in eV
    @param sw_spsym: spin channel. False = the Yosida-suppressed one (singlet, or a
                     field parallel to the triplet d-vector); True = the one preserved
                     at T=0 (field perpendicular to d, e.g. in-plane for d||z)
    @return     chi: Irreducible susceptibility in SC state [Nw, Nchi, Nchi] complex128
    """
    Nk, Nw = len(eig), len(wlist)
    Norb = eig.shape[1] // 2
    Nchi = len(ol)
    eps = idelta * 1e-3
    chi = np.zeros((Nw, Nchi, Nchi, 2), dtype=np.complex128)
    ol_f = np.asfortranarray(ol)
    _loader._lib.get_chi_irr_sc(chi, uni, eig, ffermi, qshift, ol_f, wlist,
                                _loader.i64(Nchi), _loader.i64(Norb),
                                _loader.i64(Nk), _loader.i64(Nw),
                                _loader.dbl(idelta), _loader.dbl(eps),
                                _loader.dbl(temp), ctypes.byref(ctypes.c_bool(sw_spsym)))
    return chi


def get_eig_or_tr_chi(chi: np.ndarray, invk: np.ndarray, sw_eig: bool) -> np.ndarray:
    """
    @brief Compute either the leading eigenvalue or the trace of the susceptibility on the full k-grid.
    @param    chi: Susceptibility [Nchi, Nchi, Nk] complex128
    @param   invk: Inverse k-point mapping [Nkall] int64
    @param sw_eig: If True, return leading eigenvalue; if False, return trace
    @return  chiq: Eigenvalue or trace on the full k-grid [Nkall] complex128
    """
    Nkall, Nk = len(invk), len(chi.T)
    Nchi = chi.shape[0]
    chiq = np.zeros(Nkall, dtype=np.complex128)
    _loader._lib.get_eig_or_tr_chi(chiq, chi, invk, _loader.i64(Nkall), _loader.i64(Nk),
                                   _loader.i64(Nchi), ctypes.byref(ctypes.c_bool(sw_eig)))
    return chiq
