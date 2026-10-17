#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
Impurities in a finite supercell: real-space supercell Hamiltonian with on-site impurity
potentials, exact diagonalization, and the Bloch-unfolded spectral function A(k,w) of the
primitive cell (the "spaghetti" plot of a disordered crystal).
"""
import numpy as np
import scipy.linalg as sclin
import libs.flibs as flibs
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
    @param     ham_r: host Wannier Hamiltonian H(R) [Nr, Norb, Norb]
    @param      rvec: R-vectors [Nr, 3] (integer fractional)
    @param      cell: supercell size per axis
    @param imp_sites: 0-based impurity site indices into the supercell (order of gen_rlist)
    @param     imp_V: on-site potential [eV]: a number (times the unit matrix) or [Norb, Norb]
    @return  ham_imp: supercell Hamiltonian [Norb*Nsite, Norb*Nsite] (orbital fast, site slow)
    @return    rlist: site positions [Nsite, 3] (fractional, integer valued)
    """
    Norb = ham_r.shape[1]
    rlist = np.ascontiguousarray(_lattice.gen_rlist(*cell), dtype=np.float64)
    # the supercell is periodic with a minimum-image convention: a hopping longer than half
    # the cell would be folded onto the wrong neighbour
    rmax = np.abs(np.asarray(rvec)[np.abs(ham_r).max(axis=(1, 2)) > 0]).max(axis=0)
    short = [ax for ax in range(3) if cell[ax] > 1 and 2*rmax[ax] >= cell[ax]]
    if short:
        print(f"Warning: supercell {list(cell)} is not larger than twice the hopping range "
              f"{rmax.astype(int).tolist()} along axis {short}; long hoppings wrap around", flush=True)
    i0 = np.where(np.all(np.isclose(rvec, 0.0), axis=1))[0]
    if len(i0) > 1:
        raise ValueError('impurity_supercell: rvec has several R=0 entries')
    if len(i0) == 0:
        # no on-site block in the hopping file (zero on-site energies): add one to carry imp_V
        rvec = np.vstack([rvec, np.zeros(3)])
        ham_r = np.concatenate([ham_r, np.zeros((1, Norb, Norb), dtype=ham_r.dtype)])
        i0 = [len(rvec)-1]
    V = np.asarray(imp_V, dtype=np.complex128)
    V = V*np.eye(Norb) if V.ndim == 0 else V
    ham_i = np.array(ham_r, dtype=np.complex128, copy=True)
    ham_i[i0[0]] += V
    ham_imp = flibs.gen_imp_ham(np.ascontiguousarray(rvec, dtype=np.float64), np.ascontiguousarray(ham_r),
                                ham_i, np.ascontiguousarray(ham_r), rlist,
                                np.asarray(imp_sites, dtype=np.int64))
    return ham_imp, rlist


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
