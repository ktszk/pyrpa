#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
Impurities in a finite supercell: real-space supercell Hamiltonian with on-site impurity
potentials, exact diagonalization, and the Bloch-unfolded spectral function A(k,w) of the
primitive cell (the "spaghetti" plot of a disordered crystal).
"""
import numpy as np
import scipy.linalg as sclin
import libs.plibs._bands as _bands
import libs.plibs._lattice as _lattice

# Public API of this submodule; keeps `from ._x import *` in __init__.py from
# re-exporting numpy, scipy, ctypes and the sibling packages imported above.
__all__ = [
    'random_impurity_sites', 'impurity_supercell', 'unfolded_spectrum', 'calc_impurity_spectrum'
]


def random_impurity_sites(Nsite: int, conc: float, seed: int | None = None) -> np.ndarray:
    """
    @brief Draw round(conc*Nsite) distinct impurity sites at random (at least one).
    @param Nsite: number of sites of the supercell
    @param  conc: impurity concentration (fraction of sites)
    @param  seed: seed of the random generator (None: not reproducible)
    @return sites: 0-based sorted site indices [Nimp] int64
    """
    nimp = int(min(max(round(conc*Nsite), 1), Nsite))
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(Nsite, size=nimp, replace=False)).astype(np.int64)


def impurity_supercell(ham_r: np.ndarray, rvec: np.ndarray, cell: list, imp_sites: np.ndarray,
                       imp_V: float | np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    @brief Real-space Hamiltonian of a cell[0] x cell[1] x cell[2] supercell (periodic) in which
           the sites imp_sites carry the extra on-site potential imp_V.  Hoppings are those of
           the host for every pair (host-host, host-impurity and impurity-impurity), so the
           impurity is a pure on-site perturbation, the same model as the CPA species B.
           Every H(R) is folded into the supercell: the block (site i, site j) is the sum of
           H(R) over all R with R_i + R = R_j modulo the cell.  This is the exact Hamiltonian of
           the periodic supercell for any cell size: a pristine cell reproduces H(k) at every
           commensurate k, also along an axis of size 1 (kz = 0 for cell[2] = 1, with the
           interlayer hoppings then summed into the on-site block) and when a hopping is longer
           than half the cell.  A small cell only means that an impurity sees its own images.
    @param     ham_r: host Wannier Hamiltonian H(R) [Nr, Norb, Norb] (convention of gen_ham)
    @param      rvec: R-vectors [Nr, 3] (integer fractional)
    @param      cell: supercell size per axis
    @param imp_sites: 0-based impurity site indices into the supercell (order of gen_rlist), distinct
    @param     imp_V: on-site potential [eV]: a number (times the unit matrix), a list of Norb
                      per-orbital values (diagonal) or a Hermitian [Norb, Norb] matrix
    @return  ham_imp: supercell Hamiltonian [Norb*Nsite, Norb*Nsite] (orbital fast, site slow)
    @return    rlist: site positions [Nsite, 3] (fractional, integer valued)
    """
    Norb = ham_r.shape[1]
    cell = np.asarray(cell, dtype=np.int64)
    rlist = np.ascontiguousarray(_lattice.gen_rlist(*cell), dtype=np.float64)
    Nsite = len(rlist)
    rv = np.asarray(rvec, dtype=np.float64)
    Ri = np.rint(rv).astype(np.int64)
    if not np.allclose(rv, Ri):
        raise ValueError('impurity_supercell: rvec must be integer lattice vectors')
    sites = np.asarray(imp_sites, dtype=np.int64).ravel()
    if np.any((sites < 0) | (sites >= Nsite)):
        raise ValueError(f'impurity_supercell: imp_sites must lie in [0, {Nsite-1}], got {sites.tolist()}')
    if len(np.unique(sites)) != len(sites):
        raise ValueError(f'impurity_supercell: imp_sites has duplicates: {sites.tolist()}')
    V = np.asarray(imp_V, dtype=np.complex128)
    if V.ndim == 0:
        V = V*np.eye(Norb)
    elif V.ndim == 1:
        if len(V) != Norb:
            raise ValueError(f'impurity_supercell: imp_V has {len(V)} entries, expected Norb={Norb}')
        V = np.diag(V)
    elif V.shape != (Norb, Norb):
        raise ValueError(f'impurity_supercell: imp_V must be a number, [Norb] or [Norb, Norb], got {V.shape}')
    if not np.allclose(V, V.conj().T):
        raise ValueError('impurity_supercell: imp_V must be Hermitian')

    # site index from integer position (gen_rlist order)
    rint = np.rint(rlist).astype(np.int64)
    site_of = np.empty(cell, dtype=np.int64)
    site_of[rint[:, 0], rint[:, 1], rint[:, 2]] = np.arange(Nsite)
    # gen_ham convention H(k) = sum_R H(R) exp(+2 pi i k.R): the block (site i, site i+R) is H(R)
    H4 = np.zeros((Nsite, Norb, Nsite, Norb), dtype=np.complex128)   # [site i, orb, site i+R, orb]
    src = np.arange(Nsite)
    for R, hR in zip(Ri, np.asarray(ham_r)):
        tgt = (rint + R) % cell
        # site i -> site i+R is a bijection for a fixed R, so no index repeats within one +=
        H4[src, :, site_of[tgt[:, 0], tgt[:, 1], tgt[:, 2]], :] += hR
    H4[sites, :, sites, :] += V
    return H4.reshape(Norb*Nsite, Norb*Nsite), rlist


def unfolded_spectrum(eigs: np.ndarray, uni: np.ndarray, rlist: np.ndarray, Norb: int, klist: np.ndarray,
                      wlist: np.ndarray, mu: float, eta: float) -> np.ndarray:
    """
    @brief Spectral function of the supercell eigenstates unfolded onto the primitive-cell k:
               A(k,w) = sum_n W_n(k) (eta/pi) / ((w - e_n + mu)^2 + eta^2),
               W_n(k) = (1/Nsite) sum_l |sum_R exp(+2 pi i k.R) u_n(l,R)|^2.
           The sign is the one for which a pristine supercell reproduces the eigenvalues of
           gen_ham's H(k) at k (and not at -k): checked on a model with complex hoppings,
           whose bands are not even in k.
           sum_n W_n(k) = Norb for every k, so a pristine supercell reproduces the band
           structure exactly (one unit of weight per band at k commensurate with the cell).
    @param  eigs: supercell eigenvalues [N] (N = Norb*Nsite)
    @param   uni: supercell eigenvectors as columns [N, N] (orbital fast, site slow)
    @param rlist: site positions [Nsite, 3] (fractional)
    @param  Norb: orbitals per site
    @param klist: primitive-cell k-points [Nk, 3] (fractional)
    @param wlist: energies relative to mu [Nw]
    @param    mu: chemical potential [eV]
    @param   eta: Lorentzian broadening [eV]
    @return  Akw: spectral function [Nk, Nw]
    """
    Nsite = len(rlist)
    u = np.asarray(uni).reshape(Nsite, Norb, -1)                       # [R, l, n]
    phase = np.exp(2j*np.pi*np.asarray(klist).dot(np.asarray(rlist).T))   # [k, R]
    amp = np.einsum('kr,rln->kln', phase, u)
    W = (np.abs(amp)**2).sum(axis=1)/Nsite                             # [k, n]
    lor = (eta/np.pi)/((wlist[None, :]-(eigs[:, None]-mu))**2+eta**2)   # [n, w]
    return W.dot(lor)


def calc_impurity_spectrum(ham_r: np.ndarray, rvec: np.ndarray, cell: list, imp_sites: np.ndarray,
                           imp_V: float | np.ndarray, fill: float, temp: float, k_sets: list, kmesh: int,
                           bvec: np.ndarray, Emin: float, Emax: float, Nw: int, eta: float) -> dict:
    """
    @brief Impurity spectral function along a symmetry line: build the supercell with
           impurity_supercell, diagonalize it, fix mu at the filling fill (electrons per site,
           as in calc_mu), and unfold the eigenstates (unfolded_spectrum).
    @param ham_r, rvec: host Wannier Hamiltonian and R-vectors
    @param cell, imp_sites, imp_V: supercell, impurity sites and potential (impurity_supercell)
    @param fill, temp: filling per site and temperature [eV]
    @param k_sets, kmesh, bvec: symmetry line of the primitive cell
    @param Emin, Emax, Nw: energy window relative to mu [eV] and number of points
    @param        eta: broadening [eV]
    @return dict with spa_length, xticks, wlist, Akw [Nk, Nw], mu, eigs, Nsite, Nimp
    """
    Norb = ham_r.shape[1]
    ham_imp, rlist = impurity_supercell(ham_r, rvec, cell, imp_sites, imp_V)
    Nsite = len(rlist)
    print(f"impurity supercell {list(cell)}: {Nsite} sites, {len(imp_sites)} impurities "
          f"(concentration {len(imp_sites)/Nsite:.3f}), matrix size {Norb*Nsite}", flush=True)
    eigs, uni = sclin.eigh(ham_imp)
    mu = _bands.calc_mu_imp(eigs, Nsite, fill, temp)
    print(f'chem. pot. = {mu:7.4f}', flush=True)
    klist, spa_length, xticks = _lattice.mk_klist(k_sets, kmesh, bvec)
    wlist = np.linspace(Emin, Emax, Nw, True)
    Akw = unfolded_spectrum(eigs, uni, rlist, Norb, klist, wlist, mu, eta)
    return {'spa_length': spa_length, 'xticks': xticks, 'wlist': wlist, 'Akw': Akw, 'mu': mu,
            'eigs': eigs, 'Nsite': Nsite, 'Nimp': len(imp_sites)}
