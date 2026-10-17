#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
Response functions: spin susceptibility, pairing susceptibility, conductivity, and gap symmetries.
"""
import numpy as np
import libs.flibs as flibs
import libs.plibs._bands as _bands
import libs.plibs._lattice as _lattice
import collections.abc

# Public API of this submodule; keeps `from ._x import *` in __init__.py from
# re-exporting numpy, scipy, ctypes and the sibling packages imported above.
__all__ = [
    'get_conductivity', 'chis_spectrum', 'chis_q_point', 'chis_q_point_sc', 'chis_spectrum_sc', 'chis_qmap',
    'phi_spectrum', 'phi_q_point', 'phi_qmap', 'get_chi_orb_list', 'expand_UJ_sites', 'site_block_mask',
    'set_harmonic_lattice', 'lattice_harmonic', 'gap_symms', 'get_initial_gap', 'sc_gap_on_mesh', 'bdg_bands_path'
]


class _NullPathWriter:
    def __enter__(self) -> '_NullPathWriter':
        return self

    def __exit__(self, exc_type: type | None, exc: BaseException | None, tb: object) -> bool:
        return False

    def begin_q(self, i: int, q: np.ndarray) -> None:
        return

    def write_trace(self, i: int, wlist: np.ndarray, trace: np.ndarray) -> None:
        return


class _PathWriter:
    def __init__(self, q_path: str = 'writeq.dat', trace_path: str | None = None,
                 trace_formatter: collections.abc.Callable | None = None) -> None:
        self.q_path = q_path
        self.trace_path = trace_path
        self.trace_formatter = trace_formatter or self._default_trace_formatter
        self._qfile = None
        self._trace_file = None

    @staticmethod
    def _default_trace_formatter(i: int, w: float, trval: complex) -> str:
        return f'{i:8.4f} {w:8.4f} {trval.imag:9.4f} {trval.real:9.4f}\n'

    def __enter__(self) -> '_PathWriter':
        self._qfile = open(self.q_path, 'w')
        if self.trace_path is not None:
            self._trace_file = open(self.trace_path, 'w')
        return self

    def __exit__(self, exc_type: type | None, exc: BaseException | None, tb: object) -> bool:
        if self._trace_file is not None:
            self._trace_file.close()
        if self._qfile is not None:
            self._qfile.close()
        return False

    def begin_q(self, i: int, q: np.ndarray) -> None:
        self._qfile.write(f'{i:d} {q[0]:5.3f} {q[1]:5.3f} {q[2]:5.3f}\n')
        self._qfile.flush()

    def write_trace(self, i: int, wlist: np.ndarray, trace: np.ndarray) -> None:
        if self._trace_file is None:
            return
        for w, trval in zip(wlist, trace):
            self._trace_file.write(self.trace_formatter(i, w, trval))
        self._trace_file.write('\n')
        self._trace_file.flush()


def _solve_chis_path(qlist: np.ndarray, wlist: np.ndarray, olist: np.ndarray, Smat: np.ndarray,
                     chi0_builder: collections.abc.Callable,
                     writer: _PathWriter | _NullPathWriter | None = None) -> tuple[np.ndarray, np.ndarray]:
    chisq = []
    chis_orbq = []
    writer = _NullPathWriter() if writer is None else writer
    with writer:
        for i, q in enumerate(qlist):
            writer.begin_q(i, q)
            chi0 = chi0_builder(q)
            chis = flibs.get_chis(chi0, Smat)
            trchis, trchi0, chis_orb = flibs.get_tr_chi(chis, chi0, olist)
            chisq.append(trchis)
            chis_orbq.append(chis_orb)
            writer.write_trace(i, wlist, trchi0)
    return np.array(chisq), np.array(chis_orbq)


def _prepare_sc_chi0_builder(hamk: np.ndarray, delta_k: np.ndarray, mu: float, temp: float,
                             klist: np.ndarray, olist: np.ndarray, wlist: np.ndarray,
                             idelta: float, sw_spsym: bool) -> collections.abc.Callable:
    Norb = hamk.shape[1]
    eig_BdG, uni_BdG = flibs.get_eig(flibs.mkBdGhamk(hamk - mu * np.eye(Norb), delta_k))
    # mu=0: chemical potential is already absorbed into hamBdGk (hamk - mu*I), so BdG
    # quasi-particle energies are measured from zero and the Fermi level is at 0.
    ffermi_BdG = flibs.get_ffermi(eig_BdG, 0., temp)

    def chi0_builder(q: np.ndarray) -> np.ndarray:
        qshift = flibs.get_qshift(klist, q)
        chi0sc = flibs.get_chi_irr_sc(uni_BdG, eig_BdG, ffermi_BdG, qshift, olist, wlist, idelta, temp, sw_spsym)
        return chi0sc[:, :, :, 0]+chi0sc[:, :, :, 1]

    return chi0_builder


def _solve_phi_path(qlist: np.ndarray, wlist: np.ndarray, olist: np.ndarray, phi_builder: collections.abc.Callable,
                    writer: _PathWriter | _NullPathWriter | None = None) -> tuple[np.ndarray, np.ndarray]:
    phiq = []
    phi_orbq = []
    writer = _NullPathWriter() if writer is None else writer
    with writer:
        for i, q in enumerate(qlist):
            writer.begin_q(i, q)
            phi = phi_builder(q)
            trphi, phi_orb = flibs.get_tr_phi(phi, olist)
            phiq.append(trphi)
            phi_orbq.append(phi_orb)
            writer.write_trace(i, wlist, trphi)
    return np.array(phiq), np.array(phi_orbq)


def get_conductivity(mu: float, temp: float, eig: np.ndarray, vk: np.ndarray, Nw: int, Emax: float,
                     idelta: float = 1.e-3,
                     kweight: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    @brief Compute transport coefficients L11, L12, L22 as a function of frequency using the Kubo formula.
    @param      mu: Chemical potential in eV
    @param    temp: Temperature in eV
    @param     eig: Eigenvalues array [Nk, Nband]
    @param      vk: Group velocities array [Nk, Nband, 3]
    @param      Nw: Number of frequency points
    @param    Emax: Maximum frequency in eV
    @param  idelta: Broadening parameter (Lorentzian width) in eV
    @param kweight: k-point weights [Nk] (None: uniform mesh)
    @return    L11: Charge conductivity spectrum array [Nw, 3, 3]
    @return    L12: Thermoelectric coefficient spectrum array [Nw, 3, 3]
    @return    L22: Thermal conductivity spectrum array [Nw, 3, 3]
    @return  wlist: Frequency mesh array [Nw]
    """
    wlist = np.linspace(0, Emax, Nw)
    ffermi = flibs.get_ffermi(eig, mu, temp)
    L11, L12, L22 = flibs.calc_Lij_wl(eig, vk, ffermi, mu, wlist, idelta, temp, kweight)
    return L11, L12, L22, wlist


def chis_spectrum(mu: float, temp: float, Smat: np.ndarray, klist: np.ndarray, qlist: np.ndarray,
                  olist: np.ndarray, eig: np.ndarray, uni: np.ndarray, Nw: int, Emax: float,
                  idelta: float = 1.e-3) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    @brief Compute the spin susceptibility chi_s spectrum along a q-path and write results to 'chi0.dat'.
    @param      mu: Chemical potential in eV
    @param    temp: Temperature in eV
    @param    Smat: Stoner interaction matrix [Norb, Norb]
    @param   klist: k-point list [Nk, 3]
    @param   qlist: q-point path list [Nq, 3]
    @param   olist: Orbital index list for susceptibility calculation
    @param     eig: Eigenvalues array [Nk, Nband]
    @param     uni: Eigenvectors array [Nk, Norb, Nband]
    @param      Nw: Number of frequency points
    @param    Emax: Maximum frequency in eV
    @param  idelta: Broadening parameter (Lorentzian width) in eV
    @return   chisq: Spin susceptibility trace at each q-point [Nq, Nw]
    @return chis_orbq: Orbital-resolved spin susceptibility [Nq, ...]
    @return   wlist: Frequency mesh array [Nw]
    """
    ffermi = flibs.get_ffermi(eig, mu, temp)
    wlist = np.linspace(0, Emax, Nw)

    def chi0_builder(q: np.ndarray) -> np.ndarray:
        qshift = flibs.get_qshift(klist, q)
        return flibs.get_chi_irr(uni, eig, ffermi, qshift, olist, wlist, idelta, temp)

    chisq, chis_orbq = _solve_chis_path(qlist, wlist, olist, Smat, chi0_builder, _PathWriter(trace_path='chi0.dat'))
    return chisq, chis_orbq, wlist


def chis_q_point(q: np.ndarray, eig: np.ndarray, uni: np.ndarray, Emax: float,
                 Nw: int, mu: float, temp: float, Smat: np.ndarray, klist: np.ndarray,
                 olist: np.ndarray, idelta: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    @brief Compute the spin susceptibility chi_s at a single q-point across a frequency mesh.
    @param       q: Single q-vector in fractional coordinates [3]
    @param     eig: Eigenvalues array [Nk, Nband]
    @param     uni: Eigenvectors array [Nk, Norb, Nband]
    @param    Emax: Maximum frequency in eV
    @param      Nw: Number of frequency points
    @param      mu: Chemical potential in eV
    @param    temp: Temperature in eV
    @param    Smat: Stoner interaction matrix [Norb, Norb]
    @param   klist: k-point list [Nk, 3]
    @param   olist: Orbital index list for susceptibility calculation
    @param  idelta: Broadening parameter (Lorentzian width) in eV
    @return  trchis: Trace of spin susceptibility chi_s [Nw]
    @return chis_orb: Orbital-resolved spin susceptibility
    @return   wlist: Frequency mesh array [Nw]
    """
    qlist = np.asarray([q], dtype=np.float64)
    ffermi = flibs.get_ffermi(eig, mu, temp)
    wlist = np.linspace(0, Emax, Nw)

    def chi0_builder(qvec: np.ndarray) -> np.ndarray:
        qshift = flibs.get_qshift(klist, qvec)
        return flibs.get_chi_irr(uni, eig, ffermi, qshift, olist, wlist, idelta, temp)

    chisq, chis_orbq = _solve_chis_path(qlist, wlist, olist, Smat, chi0_builder)
    return chisq[0], chis_orbq[0], wlist


def chis_q_point_sc(q: np.ndarray, hamk: np.ndarray, delta_k: np.ndarray, mu: float,
                    Emax: float, Nw: int, temp: float, Smat: np.ndarray,
                    klist: np.ndarray, olist: np.ndarray, idelta: float,
                    sw_spsym: bool = False) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    @brief Compute the spin susceptibility chi_s at a single q-point in the superconducting state.
    Builds the BdG Hamiltonian from hamk and the gap function delta_k, diagonalizes it,
    and calls get_chi_irr_sc to evaluate the irreducible susceptibility.
    @param        q: Single q-vector in fractional coordinates [3]
    @param     hamk: Normal-state k-space Hamiltonian [Nk, Norb, Norb] complex128
    @param  delta_k: Gap function (anomalous potential) on k-mesh [Nk, Norb, Norb] complex128
    @param       mu: Chemical potential in eV
    @param     Emax: Maximum frequency in eV
    @param       Nw: Number of frequency points
    @param     temp: Temperature in eV
    @param     Smat: Stoner interaction matrix [Nchi, Nchi] float64
    @param    klist: k-point list [Nk, 3] float64
    @param    olist: Orbital index pairs for susceptibility [Nchi, 2] int64 (chiolist)
    @param   idelta: Lorentzian broadening in eV
    @param sw_spsym: spin channel. False = the Yosida-suppressed one (singlet, or a
                     field parallel to the triplet d-vector); True = the one preserved
                     at T=0 (field perpendicular to d, e.g. in-plane for d||z)
    @return  trchis: Trace of spin susceptibility chi_s [Nw] complex128
    @return chis_orb: Orbital-resolved spin susceptibility
    @return   wlist: Frequency mesh [Nw] float64
    """
    qlist = np.asarray([q], dtype=np.float64)
    wlist = np.linspace(0, Emax, Nw)
    chi0_builder = _prepare_sc_chi0_builder(hamk, delta_k, mu, temp, klist, olist, wlist, idelta, sw_spsym)
    chisq, chis_orbq = _solve_chis_path(qlist, wlist, olist, Smat, chi0_builder)
    return chisq[0], chis_orbq[0], wlist


def chis_spectrum_sc(mu: float, temp: float, Smat: np.ndarray, hamk: np.ndarray,
                     delta_k: np.ndarray, klist: np.ndarray, qlist: np.ndarray,
                     olist: np.ndarray, Nw: int, Emax: float, idelta: float = 1.e-3,
                     sw_spsym: bool = False) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    @brief Compute the spin susceptibility chi_s spectrum along a q-path in the superconducting state.
    The BdG Hamiltonian is built and diagonalized once before the q-loop.
    Results are also written to 'chi0_sc.dat'.
    @param      mu: Chemical potential in eV
    @param    temp: Temperature in eV
    @param    Smat: Stoner interaction matrix [Nchi, Nchi] float64
    @param    hamk: Normal-state k-space Hamiltonian [Nk, Norb, Norb] complex128
    @param delta_k: Gap function (anomalous potential) on k-mesh [Nk, Norb, Norb] complex128
    @param   klist: k-point list [Nk, 3] float64
    @param   qlist: q-point path list [Nq, 3] float64
    @param   olist: Orbital index pairs for susceptibility [Nchi, 2] int64 (chiolist)
    @param      Nw: Number of frequency points
    @param    Emax: Maximum frequency in eV
    @param  idelta: Lorentzian broadening in eV
    @param sw_spsym: spin channel. False = the Yosida-suppressed one (singlet, or a
                     field parallel to the triplet d-vector); True = the one preserved
                     at T=0 (field perpendicular to d, e.g. in-plane for d||z)
    @return   chisq: Spin susceptibility trace at each q-point [Nq, Nw] complex128
    @return chis_orbq: Orbital-resolved spin susceptibility [Nq, ...]
    @return   wlist: Frequency mesh [Nw] float64
    """
    wlist = np.linspace(0, Emax, Nw)
    chi0_builder = _prepare_sc_chi0_builder(hamk, delta_k, mu, temp, klist, olist, wlist, idelta, sw_spsym)
    chisq, chis_orbq = _solve_chis_path(qlist, wlist, olist, Smat, chi0_builder, _PathWriter(trace_path='chi0_sc.dat'))
    return chisq, chis_orbq, wlist


def chis_qmap(Nx: int, Ny: int, Ecut: float, mu: float, temp: float, Smat: np.ndarray,
              klist: np.ndarray, olist: np.ndarray, eig: np.ndarray, uni: np.ndarray,
              idelta: float = 1.e-3) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    @brief Compute the spin susceptibility chi_s map on the qx-qy plane at a fixed energy cutoff.
    @param     Nx: Number of q-points along qx
    @param     Ny: Number of q-points along qy
    @param   Ecut: Energy cutoff for the susceptibility integration in eV
    @param     mu: Chemical potential in eV
    @param   temp: Temperature in eV
    @param   Smat: Stoner interaction matrix [Nchi, Nchi]
    @param  klist: k-point list [Nk, 3]
    @param  olist: Orbital-pair list of chi [Nchi, 2] int64 (1-based, from get_chi_orb_list)
    @param    eig: Eigenvalues array [Nk, Nband]
    @param    uni: Eigenvectors array [Nk, Norb, Nband]
    @param idelta: Broadening parameter (Lorentzian width) in eV
    @return  chis: Spin susceptibility map [Nx, Ny]
    @return  chi0: Bare susceptibility map [Nx, Ny]
    @return    qx: qx coordinate mesh [Nx, Ny]
    @return    qy: qy coordinate mesh [Nx, Ny]
    """
    ffermi = flibs.get_ffermi(eig, mu, temp)
    chis, chi0 = flibs.chiq_map(uni, eig, ffermi, klist, Smat, olist, Nx, Ny, temp, Ecut, idelta)
    x0 = np.linspace(0, 1, Nx, False)
    y0 = np.linspace(0, 1, Ny, False)
    # indexing='ij' so qx,qy have shape [Nx,Ny] matching the chi maps (chi[ix,iy])
    qx, qy = np.meshgrid(x0, y0, indexing='ij')
    return chis, chi0, qx, qy


def phi_spectrum(mu: float, temp: float, klist: np.ndarray, qlist: np.ndarray, olist: np.ndarray,
                 eig: np.ndarray,   uni: np.ndarray, Nw: int, Emax: float,
                 idelta: float = 1.e-3) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
        @brief Compute the pairing susceptibility phi spectrum along a q-path (anomalous susceptibility).
            Results are also written to 'phi0.dat'.
    @param      mu: Chemical potential in eV
    @param    temp: Temperature in eV
    @param   klist: k-point list [Nk, 3]
    @param   qlist: q-point path list [Nq, 3]
    @param   olist: Orbital index list for phi calculation
    @param     eig: Eigenvalues array [Nk, Nband]
    @param     uni: Eigenvectors array [Nk, Norb, Nband]
    @param      Nw: Number of frequency points
    @param    Emax: Maximum frequency in eV
    @param  idelta: Broadening parameter (Lorentzian width) in eV
    @return    phiq: Pairing susceptibility trace at each q-point [Nq, Nw]
    @return phi_orbq: Orbital-resolved pairing susceptibility [Nq, ...]
    @return   wlist: Frequency mesh array [Nw]
    """
    ffermi = flibs.get_ffermi(eig, mu, temp)
    wlist = np.linspace(0, Emax, Nw)

    def phi_builder(q: np.ndarray) -> np.ndarray:
        qshift = flibs.get_iqshift(klist, q)
        return flibs.get_phi_irr(uni, eig, ffermi, qshift, olist, wlist, idelta, mu, temp)

    phiq, phi_orbq = _solve_phi_path(qlist, wlist, olist, phi_builder,
                                     _PathWriter(trace_path='phi0.dat'))
    return phiq, phi_orbq, wlist


def phi_q_point(q: np.ndarray, eig: np.ndarray, uni: np.ndarray, Emax: float,
                Nw: int, mu: float, temp: float, klist: np.ndarray,
                olist: np.ndarray, idelta: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    @brief Compute the pairing susceptibility phi at a single q-point across a frequency mesh.
    @param       q: Single q-vector in fractional coordinates [3]
    @param     eig: Eigenvalues array [Nk, Nband]
    @param     uni: Eigenvectors array [Nk, Norb, Nband]
    @param    Emax: Maximum frequency in eV
    @param      Nw: Number of frequency points
    @param      mu: Chemical potential in eV
    @param    temp: Temperature in eV
    @param   klist: k-point list [Nk, 3]
    @param   olist: Orbital index list for phi calculation
    @param  idelta: Broadening parameter (Lorentzian width) in eV
    @return   trphi: Trace of pairing susceptibility phi [Nw]
    @return phi_orb: Orbital-resolved pairing susceptibility
    @return   wlist: Frequency mesh array [Nw]
    """
    qlist = np.asarray([q], dtype=np.float64)
    ffermi = flibs.get_ffermi(eig, mu, temp)
    wlist = np.linspace(0, Emax, Nw)

    def phi_builder(qvec: np.ndarray) -> np.ndarray:
        qshift = flibs.get_iqshift(klist, qvec)
        return flibs.get_phi_irr(uni, eig, ffermi, qshift, olist, wlist, idelta, mu, temp)

    phiq, phi_orbq = _solve_phi_path(qlist, wlist, olist, phi_builder)
    return phiq[0], phi_orbq[0], wlist


def phi_qmap(Nx: int, Ny: int, Ecut: float, mu: float, temp: float, klist: np.ndarray,
             olist: np.ndarray, eig: np.ndarray, uni: np.ndarray, idelta: float = 1.e-3,
             sw_omega: bool = True) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    @brief Compute the pairing susceptibility phi map on the qx-qy plane at a fixed energy cutoff.
    @param      Nx: Number of q-points along qx
    @param      Ny: Number of q-points along qy
    @param    Ecut: Energy cutoff for the susceptibility integration in eV
    @param      mu: Chemical potential in eV
    @param    temp: Temperature in eV
    @param   klist: k-point list [Nk, 3]
    @param   olist: Orbital-pair list of phi [Nchi, 2] int64 (1-based, from get_chi_orb_list)
    @param     eig: Eigenvalues array [Nk, Nband]
    @param     uni: Eigenvectors array [Nk, Norb, Nband]
    @param  idelta: Broadening parameter (Lorentzian width) in eV
    @param sw_omega: If True, integrate over Matsubara frequencies; if False, use static limit
    @return    phi: Pairing susceptibility map [Nx, Ny]
    @return     qx: qx coordinate mesh [Nx, Ny]
    @return     qy: qy coordinate mesh [Nx, Ny]
    """
    ffermi = flibs.get_ffermi(eig, mu, temp)
    phi = flibs.phiq_map(uni, eig, ffermi, klist, olist, Nx, Ny, mu, temp, Ecut, idelta, sw_omega)
    x0 = np.linspace(0, 1, Nx, False)
    y0 = np.linspace(0, 1, Ny, False)
    # indexing='ij' so qx,qy have shape [Nx,Ny] matching the phi map (phi[ix,iy])
    qx, qy = np.meshgrid(x0, y0, indexing='ij')
    return phi, qx, qy


def get_chi_orb_list(Norb: int, site_prof: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    @brief Generate the orbital index pair list for susceptibility calculation from site profile.
    @param    Norb: Total number of orbitals in the Hamiltonian
    @param site_prof: Array of orbital counts per site (e.g. [3, 3] for two 3-orbital sites)
    @return chiolist: Array of orbital index pairs (1-based) [Npairs, 2] int64, FORTRAN-contiguous:
                      the Fortran vertex/chi routines declare it as ol(Nchi,2) and read it in
                      column-major order, so a C-contiguous list would be silently transposed
    @return    site: Array of site indices corresponding to each orbital pair [Npairs]
    """
    if (len(site_prof) == 1):
        tmp = np.arange(Norb)+1
        o1, o2 = np.meshgrid(tmp, tmp)
        chiolist = np.array([o1.flatten(), o2.flatten()]).T
        site = np.ones(len(chiolist), dtype=np.int64)
    else:
        if (Norb == sum(site_prof)):
            chiolist = []
            site = []
            N0 = 1
            for i_site in site_prof:
                tmp = np.arange(i_site)+N0
                o1, o2 = np.meshgrid(tmp, tmp)
                chiolist += list(np.array([o1.flatten(), o2.flatten()]).T)
                # One site label per orbital pair (i_site^2 entries), not per row (i_site).
                site += [N0]*o1.size
                N0 += i_site
            chiolist = np.array(chiolist)
            site = np.array(site, dtype=np.int64)
        else:
            raise ValueError(f"site_prof={list(site_prof)} does not add up to Norb={Norb}")
    # Fortran memory order: the single-site branch gets it from the .T view above, the
    # multi-site one is built row by row and has to be converted.
    return np.asfortranarray(chiolist, dtype=np.int64), site


def expand_UJ_sites(Umat: np.ndarray, Jmat: np.ndarray, Norb: int,
                    site_prof: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    @brief Bring an orbital-dependent U/J pair into the Norb x Norb form the vertex
    routines index with the GLOBAL orbital number, replicating a per-site matrix on
    every site when one site worth of orbitals is given.
    The vertices only read pairs sitting on one site (site(i)==site(j) in get_scmat_orb),
    so the inter-site blocks are never used and are left at zero.
    @param      Umat: [Norb,Norb], or [Norb/Nsite, Norb/Nsite] when every site carries the
                      same number of orbitals (then replicated on each site)
    @param      Jmat: same shape as Umat
    @param      Norb: total number of orbitals of the hamiltonian
    @param site_prof: orbital counts per site, e.g. [5,5]; [Norb] or [1] for a single site
    @return Umat,Jmat: [Norb,Norb] float64, FORTRAN-contiguous, i.e. Umat[a,b] is the
                      Umat(a+1,b+1) the Fortran vertex routines read (immaterial for the
                      symmetric U/J of a Kanamori interaction, but the declared convention)
    """
    mats = []
    for name, mat in (('Umat', Umat), ('Jmat', Jmat)):
        arr = np.asarray(mat)
        if np.iscomplexobj(arr):
            raise ValueError(f'{name} must be a real-valued matrix')
        if arr.ndim != 2 or arr.shape[0] != arr.shape[1]:
            raise ValueError(f'{name} must be a square matrix, got shape {arr.shape}')
        mats.append(np.asfortranarray(arr, dtype=np.float64))
    Umat, Jmat = mats
    if Umat.shape != Jmat.shape:
        raise ValueError(f'Umat and Jmat shapes differ: {Umat.shape} vs {Jmat.shape}')
    Nsite = len(site_prof)
    if Umat.shape[0] == Norb:      # already indexed by the global orbital number
        return Umat, Jmat
    Nsorb = Umat.shape[0]
    if Nsite < 2:
        raise ValueError(f'Umat/Jmat shape must be ({Norb}, {Norb}), got {Umat.shape}')
    if sum(site_prof) != Norb:
        raise ValueError(f'site_prof={list(site_prof)} does not add up to Norb={Norb}')
    if any(ns != Nsorb for ns in site_prof):
        raise ValueError(f'Umat/Jmat shape must be ({Norb}, {Norb}), got {Umat.shape}: a per-site '
                         f'matrix is replicated only when every site has the same number of '
                         f'orbitals, but site_prof={list(site_prof)}')
    # one site worth of orbitals: put the same matrix on every site (block diagonal)
    Uout = np.zeros((Norb, Norb), dtype=np.float64, order='F')
    Jout = np.zeros((Norb, Norb), dtype=np.float64, order='F')
    for i_site in range(Nsite):
        n0 = i_site*Nsorb
        Uout[n0:n0+Nsorb, n0:n0+Nsorb] = Umat
        Jout[n0:n0+Nsorb, n0:n0+Nsorb] = Jmat
    return Uout, Jout


def site_block_mask(Norb: int, site_prof: np.ndarray) -> np.ndarray:
    """
    @brief Boolean [Norb,Norb] mask of the orbital pairs living on one site, i.e. exactly
    the elements of Umat/Jmat the on-site interaction vertices read.
    @param      Norb: total number of orbitals
    @param site_prof: orbital counts per site
    @return     mask: [Norb,Norb] bool
    """
    mask = np.zeros((Norb, Norb), dtype=bool)
    if len(site_prof) < 2:
        mask[:, :] = True
        return mask
    n0 = 0
    for ns in site_prof:
        mask[n0:n0+ns, n0:n0+ns] = True
        n0 += ns
    return mask


# --- Cartesian (R-shell) lattice harmonics ------------------------------------
# The fractional harmonics in gap_symms (cos(2*pi*k_x) etc.) carry their labelled
# symmetry only on a tetragonal / orthogonal single-site lattice.  Fixing that does
# NOT require the axes to be decomposed (sw_dec_axis, which would break the k-mesh
# periodicity the BZ sums and FFT convolutions rely on): the phase
# k.R = 2*pi*k_frac.n is basis independent, so it is enough to attach the CARTESIAN
# symmetry to the shell COEFFICIENTS,
#
#     phi_G(k) = sum_R K_G(R_hat) exp(i k.R),        R_cart = n . avec,
#
# which is periodic in G by construction (R is a real lattice vector) and transforms
# as the labelled irrep in any Bravais lattice.  On an orthogonal single-site lattice
# it reproduces the fractional formulas EXACTLY (tests/test_rpa_flex.py).
#
# The R set is fixed by seed patterns rather than by a shell index: R_cut is the
# longest of the seeds ((1,0,0)/(0,1,0) for the first-neighbour kinds, (1,+-1,0) for
# d_xy and s+-, (1,0,+-1)/(0,1,+-1) for d_xz/d_yz), and every R within R_cut enters.
# That keeps the two inequivalent stars of an orthorhombic cell (|a1| != |a2|) in the
# same harmonic, while on a high-symmetry lattice it picks up the complete star (all
# 12 nearest neighbours of an fcc cell, say) that the irrep needs.  K annihilates the
# shells that do not belong; the A1g kind ('s'), which K cannot discriminate, uses the
# single shell at R_cut instead.
#
# CAVEAT: for a multi-site cell the pairing bond is R + tau_j - tau_i, not R, so the
# sublattice phase is still missing here -- use the orbital-basis route there
# (eil_gap_orbital / eil_gap_file + the Nagai-Nakamura band projection).
_CART_K = {'s': lambda n: np.ones(len(n)),
           'dx2-y2': lambda n: n[:, 0] ** 2 - n[:, 1] ** 2,
           'dxy': lambda n: n[:, 0] * n[:, 1],
           'dxz': lambda n: n[:, 0] * n[:, 2],
           'dyz': lambda n: n[:, 1] * n[:, 2],
           'px': lambda n: n[:, 0],
           'py': lambda n: n[:, 1]}
_CART_ODD = {'px', 'py'}
# kind -> (seed patterns fixing R_cut, single shell only?, in-plane only?)
_CART_SEED = {'s':      ([(1, 1, 0), (1, -1, 0)], True, True),
              'dx2-y2': ([(1, 0, 0), (0, 1, 0)], False, False),
              'dxy':    ([(1, 1, 0), (1, -1, 0)], False, False),
              'dxz':    ([(1, 0, 1), (1, 0, -1)], False, False),
              'dyz':    ([(0, 1, 1), (0, 1, -1)], False, False),
              'px':     ([(1, 0, 0), (0, 1, 0)], False, False),
              'py':     ([(1, 0, 0), (0, 1, 0)], False, False)}
# gap_sym -> (kind, sign, target rms over the BZ).  The rms fixes the amplitude to the
# fractional convention (cos kx - cos ky and 2 sin kx sin ky have rms 1, 2 sin kx has
# rms sqrt(2)), so per-band delta0 amplitudes keep their meaning; the sign keeps the
# overall sign of the old formula.
_CART_GAP = {1: ('dx2-y2', +1.0, 1.0), 2: ('s', +1.0, 1.0), 3: ('dxy', -1.0, 1.0),
             4: ('dxz', -1.0, 1.0), 5: ('dyz', -1.0, 1.0),
             -1: ('px', +1.0, np.sqrt(2.0)), -2: ('py', +1.0, np.sqrt(2.0))}
_CART_CHIRAL = {6: (1, 3), 7: (4, 5), -3: (-1, -2)}   # chiral = partner_1 + i * partner_2

_HARM_AVEC = None
_HARM_RVEC = None


def set_harmonic_lattice(avec: np.ndarray | None = None, rvec: np.ndarray | None = None) -> None:
    """
    @brief Register the lattice that gap_symms uses to build its harmonics in the
    CARTESIAN basis (see lattice_harmonic).  main.py calls this once; with nothing
    registered gap_symms falls back to the fractional-coordinate formulas, which carry
    the labelled symmetry only for a tetragonal / orthogonal single-site lattice.
    @param avec: primitive translation vectors as rows [3,3] in Angstrom (None clears)
    @param rvec: R-vectors of the model [Nr,3] as integer fractional indices
    """
    global _HARM_AVEC, _HARM_RVEC
    _HARM_AVEC = None if avec is None else np.asarray(avec, dtype=np.float64)
    _HARM_RVEC = None if rvec is None else np.asarray(rvec, dtype=np.float64)


def _seed_rcut(kind: str, avec: np.ndarray) -> float:
    """|R| cutoff of the seed patterns of `kind` for this lattice (see _CART_SEED)."""
    seeds = _CART_SEED[kind][0]
    return max(np.linalg.norm(np.asarray(sd, dtype=np.float64)
                              @ np.asarray(avec, dtype=np.float64)) for sd in seeds)


def lattice_harmonic(klist: np.ndarray, avec: np.ndarray, rvec: np.ndarray, kind: str,
                     rms: float = None, rcut: float = None) -> np.ndarray:
    """
    @brief Symmetrized plane-wave harmonic phi(k) = sum_R K(R_hat) exp(i k.R) with K a
    CARTESIAN harmonic of the direction of R_cart = n . avec.  Periodic in the
    reciprocal lattice by construction and correct in any Bravais lattice, so the gap
    symmetry needs no axis decomposition.
    @param klist: k-points in FRACTIONAL coordinates [Nk,3] (the phase 2*pi*k.n is basis
                  independent, so they need no transformation)
    @param  avec: primitive translation vectors as rows [3,3]
    @param  rvec: R-vectors [Nr,3] as integer fractional indices
    @param  kind: 's','dx2-y2','dxy','dxz','dyz','px','py'
    @param   rms: rescale so that the BZ root-mean-square of phi is this (None: raw sum)
    @param  rcut: |R| cutoff in Angstrom (None: from the seed patterns of this kind)
    @return  phi: real array [Nk]  (cosine sum for even K, sine sum for odd K)
    """
    if kind not in _CART_K:
        raise ValueError(f"lattice_harmonic: unknown kind {kind!r}")
    rv = np.asarray(rvec, dtype=np.float64)
    av = np.asarray(avec, dtype=np.float64)
    Rc = rv @ av
    d = np.linalg.norm(Rc, axis=1)
    seeds, one_shell, in_plane = _CART_SEED[kind]
    if rcut is None:
        rcut = _seed_rcut(kind, av)
    radii = np.unique(np.round(d[d > 1.0e-8], 6))
    for rc in radii[radii >= rcut * (1.0 - 1.0e-6)] if len(radii) else []:
        # The seed cutoff is the right shell on an orthogonal lattice.  Where the
        # harmonic happens to vanish on it -- d_x2-y2 on the <111> neighbours of a bcc
        # cell, say -- walk out to the next shell instead of returning zero.
        sel = (d <= rc * (1.0 + 1.0e-6)) & (d > 1.0e-8)
        if one_shell:
            sel &= np.isclose(d, rc, rtol=1.0e-6)
        if in_plane:
            # s+- is A1g, which every shell carries, so K cannot pick the shell out: it
            # is defined by the IN-PLANE nesting, hence the R_z = 0 restriction (this is
            # what makes it 2 cos kx cos ky rather than the full cubic (110) sum).
            sel &= np.abs(Rc[:, 2]) < 1.0e-8
        if not sel.any():
            continue
        K = _CART_K[kind](Rc[sel] / d[sel][:, None])
        if np.abs(K).max() > 1.0e-12:
            break
    else:
        raise ValueError(f"lattice_harmonic: {kind} vanishes on every R shell of this "
                         "model (or the R set is too small for this gap symmetry)")
    ph = np.exp(2.0j * np.pi * (np.asarray(klist, dtype=np.float64) @ rv[sel].T))
    v = (ph * K).sum(axis=1)
    phi = v.imag if kind in _CART_ODD else v.real
    if rms is not None:
        # BZ rms is analytic: exp(i k.R) are orthonormal over the BZ and the R set is
        # closed under inversion, so <phi^2> = sum_R K_R^2 for both parities.
        cur = np.sqrt((K ** 2).sum())
        phi = phi * (rms / cur) if cur > 0.0 else phi
    return phi


def _cart_gap_row(klist: np.ndarray, gap_sym: int, avec: np.ndarray,
                  rvec: np.ndarray, rcut: float = None) -> np.ndarray:
    """One k-row of the Cartesian harmonic for gap_sym, normalized to the amplitude
    convention of the fractional formulas in gap_symms."""
    if gap_sym in _CART_CHIRAL:
        # A chiral form factor is only a proper 2D-irrep partner pair if BOTH members are
        # built on the SAME R set -- otherwise |phi| is not invariant under the rotation
        # that mixes them (hexagonal E2g fails by 77% with per-kind shells, 3e-15 with a
        # common one).  The smaller of the two seed cutoffs is used, and the per-kind
        # walk-out below restores the old shell where the partner vanishes on it (a
        # tetragonal cell, where the two are separate 1D irreps anyway).
        g1, g2 = _CART_CHIRAL[gap_sym]
        rc = min(_seed_rcut(_CART_GAP[g][0], avec) for g in (g1, g2))
        return (_cart_gap_row(klist, g1, avec, rvec, rcut=rc)
                + 1.0j * _cart_gap_row(klist, g2, avec, rvec, rcut=rc))
    kind, sgn, rms = _CART_GAP[gap_sym]
    return sgn * lattice_harmonic(klist, avec, rvec, kind, rms=rms, rcut=rcut)


def gap_symms(klist: np.ndarray, Norb: int, gap_sym: int,
              avec: np.ndarray = None, rvec: np.ndarray = None) -> np.ndarray:
    """
    @brief Pairing form factor phi(k) for one gap symmetry index, tiled over orbitals.

    The form factor depends only on k, so a single row is built and broadcast to
    every orbital; the caller supplies the orbital structure separately.  Two
    routes exist and the Cartesian one is preferred when it is available:

      * CARTESIAN (avec and rvec given, or set once via set_harmonic_lattice) --
        the harmonic is summed over the real R-shell of the model, so it carries
        the actual point-group content of the lattice.  Used for every gap_sym
        listed in _CART_GAP / _CART_CHIRAL.  If the R set is too small for the
        requested symmetry a warning is printed and the fractional form is used.
      * FRACTIONAL (fallback) -- the textbook cos(2*pi*k_x) - cos(2*pi*k_y) style
        expressions.  CAUTION: these carry the labelled symmetry only for a
        tetragonal / orthogonal SINGLE-SITE cell.  On a non-orthogonal lattice or
        a multi-site cell (the 2-Fe cell of the iron-based systems, say) the label
        and the true irreducible representation no longer coincide; prefer the
        orbital-basis route there (gap_orbital + project_gap_to_band).

    Indices 6, 7 and -3 are chiral and return a COMPLEX row; an unknown index
    returns zeros rather than raising, so a mistyped gap_sym shows up as a
    vanishing gap rather than a crash.

    @param   klist: k-points in fractional coordinates [Nk, 3] float64
    @param    Norb: number of orbitals; the row is tiled to [Norb, Nk]
    @param gap_sym: 0=s, 1=dx2-y2, 2=s+-, 3=dxy, 4=dxz, 5=dyz, 6=d+id, 7=dxz+idyz,
                    -1=px, -2=py, -3=p+ip (negative = triplet)
    @param    avec: primitive translation vectors as rows [3, 3] in Angstrom
                    (None = the value stored by set_harmonic_lattice)
    @param    rvec: R-vectors of the model [Nr, 3] as integer fractional indices
                    (None = the value stored by set_harmonic_lattice)
    @return phi: form factor [Norb, Nk], complex for the chiral indices
    """
    av = _HARM_AVEC if avec is None else np.asarray(avec, dtype=np.float64)
    rv = _HARM_RVEC if rvec is None else np.asarray(rvec, dtype=np.float64)
    if av is not None and rv is not None and (gap_sym in _CART_GAP or gap_sym in _CART_CHIRAL):
        try:
            return np.tile(_cart_gap_row(klist, gap_sym, av, rv), (Norb, 1))
        except ValueError as e:
            # too small an R set for this symmetry: keep working on the fractional
            # formula rather than aborting, but say so -- it is only labelled
            # correctly on an orthogonal single-site lattice.
            print(f'Warning: Cartesian gap harmonic unavailable ({e}); '
                  'falling back to the fractional-coordinate form factor', flush=True)
    # FRACTIONAL fallback (see the CAUTION in the docstring)
    A = 2*np.pi
    kx, ky, kz = klist[:, 0], klist[:, 1], klist[:, 2]
    if gap_sym == 0:        # s
        row = np.ones(len(klist), dtype=np.float64)
    elif gap_sym == 1:      # dx2-y2
        row = np.cos(A*kx)-np.cos(A*ky)
    elif gap_sym == 2:      # spm
        row = 2*np.cos(A*kx)*np.cos(A*ky)
    elif gap_sym == 3:      # dxy
        row = 2*np.sin(A*kx)*np.sin(A*ky)
    elif gap_sym == 4:      # dxz
        row = 2*np.sin(A*kx)*np.sin(A*kz)
    elif gap_sym == 5:      # dyz
        row = 2*np.sin(A*ky)*np.sin(A*kz)
    elif gap_sym == -1:     # px
        row = 2*np.sin(A*kx)
    elif gap_sym == -2:     # py
        row = 2*np.sin(A*ky)
    elif gap_sym == -3:     # p+ip (chiral, complex): px + i py
        row = 2*np.sin(A*kx)+2j*np.sin(A*ky)
    elif gap_sym == 6:      # d+id (chiral, complex): dx2-y2 + i dxy  (indices 1 + i*3)
        row = (np.cos(A*kx)-np.cos(A*ky))+2j*np.sin(A*kx)*np.sin(A*ky)
    elif gap_sym == 7:      # dxz+idyz (chiral, complex, kz dependent): indices 4 + i*5.
        # |phi| vanishes on the whole kz=0 (and zone-boundary) plane -> HORIZONTAL line
        # node, unlike the vertical nodes of the in-plane harmonics.
        row = 2*np.sin(A*kx)*np.sin(A*kz)+2j*np.sin(A*ky)*np.sin(A*kz)
    else:
        row = np.zeros(len(klist), dtype=np.float64)
    return np.tile(row, (Norb, 1))


def get_initial_gap(klist: np.ndarray, Norb: int, gap_sym: int,
                    avec: np.ndarray = None, rvec: np.ndarray = None) -> np.ndarray:
    """
    @brief Generate an initial gap function with the specified pairing symmetry for gap equation iteration.
    @param   klist: k-point list in fractional coordinates [Nk, 3]
    @param    Norb: Number of orbitals in the Hamiltonian
    @param gap_sym: Symmetry index: 0=s, 1=dx2-y2, 2=s+-, 3=dxy, 4=dxz, 5=dyz, -1=px, -2=py, -3=p+ip
    @param    avec: primitive translation vectors as rows [AA]
    @param    rvec: R-vectors [Nr, 3]
    @return init_gap: Initial gap function array [Norb, Nk]
    """
    if gap_sym >= 0:
        gapsym = ['s', 'dx2-y2', 'spm', 'dxy', 'dxz', 'dyz', 'd+id', 'dxz+idyz']
        print('gap symmetry is '+gapsym[gap_sym], flush=True)
    else:
        gapsym = ['s', 'px', 'py', 'p+ip']
        print('gap symmetry is '+gapsym[-gap_sym], flush=True)
    row = gap_symms(klist, Norb, gap_sym, avec=avec, rvec=rvec)
    # The linearized Eliashberg seed is passed to Fortran as real(c_double) and the
    # kernel itself is real, so a CHIRAL harmonic cannot be used as a seed: the two
    # members of a degenerate pair (e.g. dxz and dyz) must be obtained from separate
    # real runs and combined as Delta_1 +- i Delta_2 afterwards, with the chirality
    # selected below Tc by the condensation energy (linear theory cannot select it).
    if np.iscomplexobj(row) and np.abs(row.imag).max() > 0.0:
        raise ValueError(
            f'gap_sym={gap_sym} is a chiral (complex) harmonic and cannot seed the '
            'linearized Eliashberg equation (real kernel, real seed). Run the two real '
            'partners separately (e.g. 1 and 3 for d+id, 4 and 5 for dxz+idyz), check '
            'that they are degenerate and orthogonal, and build Delta_1 +- i Delta_2 '
            'afterwards. The chiral harmonics are for the Eilenberger form factor.')
    return row


def sc_gap_on_mesh(Nx: int, Ny: int, Nz: int, ham_r: np.ndarray, S_r: np.ndarray, rvec: np.ndarray,
                   gap_sym: int, delta0: float | list, plist: np.ndarray,
                   iperm: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """
    @brief Model superconducting gap Delta(k) in the orbital basis on the full k-mesh, for the
           SC-state susceptibility / NMR modes (BdG with a fixed gap).
           delta0 a number: the gap_sym form factor is projected to the orbital basis
           (get_initial_delta) and rescaled so that max|Delta| = delta0.
           delta0 a list: one amplitude per band (signs give s+-), with bands 1,2 damped by 1/2
           inside the Gamma-centred AFM diamond |kx|+|ky| < 1/2 (a model choice for the 5-band
           iron pnictide, applied only to the bands that exist).
    @param Nx, Ny, Nz: k-mesh
    @param ham_r, S_r, rvec: Wannier Hamiltonian, overlap ([] if orthogonal) and R-vectors
    @param gap_sym: gap symmetry index
    @param  delta0: maximum gap [eV], or per-band amplitudes [Norb]
    @param   plist: inversion parity per orbital (scalar delta0 only)
    @param   iperm: orbital permutation of the inversion operator (None = identity)
    @return  klist: full k-mesh on [0,1) without endpoints [Nk, 3] (required by get_qshift)
    @return deltak: gap [Nk, Norb, Norb] complex128, paired with klist
    """
    Norb = ham_r.shape[1]
    if np.isscalar(delta0):
        klist, kmap, invk = flibs.gen_irr_k_TRS(Nx, Ny, Nz)
        eig, uni = _bands.get_eigs(klist, ham_r, S_r, rvec)
        deltaini = flibs.remap_gap(flibs.get_initial_delta(get_initial_gap(klist, len(eig.T), gap_sym),
                                                           uni, kmap, invk, 1, gap_sym), plist, invk, gap_sym,
                                   iperm=iperm)
        deltaini_static = np.ascontiguousarray(deltaini[:, :, 0, :].transpose(2, 0, 1))  # [Nkall,Norb,Norb]
        print("maximum gap = %7.4f meV" % (delta0*1e3), flush=True)
        gap_max = abs(deltaini_static).max()
        if gap_max <= 0.0:
            print("Warning: initial gap is zero. Using unscaled delta", flush=True)
            deltak = deltaini_static.copy()
        else:
            deltak = delta0*deltaini_static/gap_max
        # Full-grid klist in kmap (FFT) order: keeps deltak <-> klist index pairing
        # consistent, and gives the [0,1) endpoint-free mesh required by get_qshift
        klist = kmap/np.array([Nx, Ny, Nz], dtype=np.float64)
    else:
        # sw_pp=False: get_qshift requires the periodic [0,1) mesh without endpoints
        _, klist = _lattice.gen_klist(Nx, Ny, Nz, sw_pp=False)
        eig, uni = _bands.get_eigs(klist, ham_r, S_r, rvec)
        inigap = gap_symms(klist, Norb, gap_sym)
        for i, k in enumerate(klist):
            # AFM-zone test in centered fractional coords: keep the full gap inside the
            # M-centered diamond (|kx|+|ky|>=1/2) and damp bands 1,2 inside the
            # Gamma-centered diamond (|kx|+|ky|<1/2).
            kc = k-np.round(k)
            if abs(kc[0])+abs(kc[1]) < 0.5:
                inigap[1:3, i] *= .5   # bands 1,2 where they exist
        # delta0 is per-band [Norb]; broadcast over k: (inigap.T)[Nk,Norb] * delta0[Norb] -> [Norb,Nk]
        deltaini = (inigap.T*np.array(delta0, dtype=np.complex128)).T
        deltak = flibs.get_band_to_orb_delta(deltaini, uni)
    return klist, deltak


def bdg_bands_path(ham_r: np.ndarray, S_r: np.ndarray, rvec: np.ndarray, mu: float, klist: np.ndarray,
                   gap_sym: int, delta0: float | list, Nk_path: int = 200) -> tuple[np.ndarray, np.ndarray]:
    """
    @brief BdG bands along (0,0,0) -> (0,0.5,0) for the model gap of sc_gap_on_mesh.
           The band form factor is scaled like the mesh gap: its maximum over klist is delta0
           (number), or one amplitude per band (list).
    @param ham_r, S_r, rvec: Wannier Hamiltonian, overlap ([] if orthogonal) and R-vectors
    @param      mu: chemical potential [eV]
    @param   klist: the k-mesh of the SC calculation (sets the form-factor normalization)
    @param gap_sym: gap symmetry index
    @param  delta0: maximum gap [eV], or per-band amplitudes [Norb]
    @param Nk_path: points along the path
    @return      ky: path coordinate k_y [Nk_path]
    @return eig_bdg: BdG eigenvalues [Nk_path, 2*Norb] (eV, relative to mu)
    """
    Norb = ham_r.shape[1]
    kpath_bdg = np.zeros((Nk_path, 3))
    kpath_bdg[:, 1] = np.linspace(0, 0.5, Nk_path)
    hamk_path = flibs.gen_ham(kpath_bdg, ham_r, rvec)
    _, uni_path = _bands.get_eigs(kpath_bdg, ham_r, S_r, rvec)
    inigap_path = gap_symms(kpath_bdg, Norb, gap_sym)
    if np.isscalar(delta0):
        # scale the band form factor so that its maximum over the k-mesh is delta0,
        # the amplitude deltak was normalized to
        fmax = np.abs(gap_symms(klist, 1, gap_sym)).max()
        scale = delta0 / fmax if fmax > 0 else delta0
        deltak_path = flibs.get_band_to_orb_delta((inigap_path * scale).astype(np.complex128), uni_path)
    else:
        deltak_path = flibs.get_band_to_orb_delta(
            (inigap_path.T * np.array(delta0, dtype=np.complex128)).T, uni_path)
    eig_BdG_path, _ = flibs.get_eig(flibs.mkBdGhamk(hamk_path - mu * np.eye(Norb), deltak_path))
    return kpath_bdg[:, 1], eig_BdG_path
