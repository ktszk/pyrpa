#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
Boltzmann / linear-response transport drivers: relaxation-time models, weak-field Hall
coefficient with carrier analysis, thermoelectric tensors, optical constants and colours.
"""
import os
import warnings
import numpy as np
import scipy.linalg as sclin
import scipy.constants as scconst
import libs.flibs as flibs
import libs.plibs._units as _units
import libs.plibs._bands as _bands
import libs.plibs._calc as _calc
import libs.plibs._lattice as _lattice
import libs.plibs._color as _color
import libs.plibs._optical as _optical
import libs.plibs._response as _response
import libs.plibs._wannier_io as _wannier_io

# Public API of this submodule; keeps `from ._x import *` in __init__.py from
# re-exporting numpy, scipy, ctypes and the sibling packages imported above.
__all__ = [
    'print_matrix', 'build_tau', 'flex_quasiparticles', 'print_tau_summary', 'get_hall_coe',
    'calc_conductivity_Boltzmann', 'calc_conductivity_lrt'
]


def print_matrix(label: str, mat: np.ndarray, width: float = 10, prec: int = 3, ndigits: int = 10) -> None:
    """Print a labelled N x 3 real matrix, one row per line (transport tensors)."""
    print(label, flush=True)
    for row in np.asarray(mat).round(ndigits):
        print(f" {row[0]:{width}.{prec}e} {row[1]:{width}.{prec}e} {row[2]:{width}.{prec}e}", flush=True)


# tau_modes whose tau is scaled by tau_const ('epa' and 'flex' derive tau on their own)
_TAU_CONST_MODES = ('const', 'e2', 'dos')


def _tau_const_ok(tau_const: float | None, tau_kw: dict) -> bool:
    """tau_const must be positive only where the tau_mode uses it; None is allowed otherwise."""
    if tau_kw.get('tau_mode') in _TAU_CONST_MODES and (tau_const is None or tau_const <= 0):
        print(f"Error: Relaxation time (tau_const) must be positive for tau_mode='{tau_kw.get('tau_mode')}' "
              f"({tau_const})", flush=True)
        return False
    return True


def _inv_conducting(M: np.ndarray, rtol: float = 1.0e-12) -> tuple[np.ndarray, str]:
    """
    Inverse of a 3x3 transport tensor (K0, sigma, L11) on the directions that conduct, and a
    label of those that do not ('' if none).  A strictly 2D model, or Nz = 1 (only kz = 0, where
    v_z = dE/dkz vanishes), has sigma_zz = 0 exactly: the full inverse does not exist, yet the
    in-plane S, Peltier, kappa and Lorenz are well defined, and the SVD pseudo-inverse with the
    singular values below rtol*max dropped gives exactly the inverse of that block (the
    ordinary inverse when M is regular).  Quantities along a dead direction come out 0.
    """
    u, sv, vh = np.linalg.svd(M)
    keep = sv > rtol*sv.max() if sv.max() > 0 else np.zeros(3, dtype=bool)
    inv = (vh[keep].conj().T/sv[keep]).dot(u[:, keep].conj().T)
    dead = []
    for v in vh[~keep]:
        a = int(np.argmax(np.abs(v)))
        dead.append('xyz'[a] if abs(v[a]) > 1.0-1.0e-8 else
                    '(' + ', '.join(f'{c:.2f}' for c in np.real_if_close(v)) + ')')
    return inv, ', '.join(dead)


def build_tau(tau_mode: str, units: _units.Units, eig: np.ndarray, mu: float, temp: float, tau_const: float,
              rvec: np.ndarray | None = None, ham_r: np.ndarray | None = None, S_r: np.ndarray | None = None,
              avec: np.ndarray | None = None, Nx: int | None = None, Ny: int | None = None,
              Nz: int | None = None, Nw: int = 300, epa_file: str | None = None,
              dos_delta: float = 1.0e-2) -> tuple:
    '''
    Relaxation time tau[Nk,Norb] in units of tau_unit, per the tau_mode switch.

    R_H is INDEPENDENT of a constant tau (it carries tau^2/tau^2), so the whole
    temperature dependence of the Hall coefficient -- and any difference between
    the mobilities of two Fermi-surface sheets -- lives in tau(k).  With
    tau_mode='const' the R_H printed by get_hall_coe is a pure band-structure
    number; only a k-dependent tau makes R_H(T) meaningful.

    tau_mode: 'const', 'e2' (Fermi liquid), 'dos' (DOS ratio) or 'epa' (needs epa_file);
    'flex' is not a tau alone and is handled by flex_quasiparticles.
    dos_delta is the Lorentzian broadening [eV] of the DOS used by 'dos'.
    Returns (tau, label).
    '''
    if tau_mode == 'const':
        return eig*0.+tau_const, f'constant tau = {tau_const} '+units.tau_label
    if tau_mode == 'epa':
        if not epa_file:
            raise ValueError("tau_mode='epa' needs epa_file (epa.x job='egrid' output)")
        ngrid, nmodes, edge, step, nbin, wavg, gavg = _wannier_io.read_epa_output(epa_file)
        # calc_tau_epa works with hbar=1 and energies in eV, so it returns hbar/eV;
        # the Boltzmann kernels expect tau in units of tau_unit.
        tau = flibs.calc_tau_epa(eig, gavg, wavg, edge, step, nbin, mu, temp)*units.hbar/units.tau_unit
        return tau, f'EPA tau(k) from {epa_file} ({nmodes} modes)'
    if tau_mode == 'e2':
        # Fermi-liquid electron-electron scattering: 1/tau ∝ (e-mu)^2 + (pi T)^2, normalized so
        # that tau = tau_const on the Fermi level and falls off as 1/(e-mu)^2 away from it
        pit2 = (np.pi*temp)**2
        return (tau_const*pit2/((eig-mu)**2+pit2),
                f'Fermi-liquid tau: 1/tau ~ (e-mu)^2+(pi T)^2, tau(mu) = {tau_const} '
                + units.tau_label)
    if tau_mode == 'dos':
        # golden-rule scattering rate ∝ final-state DOS, used only as a RATIO to the Fermi level:
        # tau(e) = tau_const * D(mu)/D(e), so tau = tau_const on the Fermi level
        _, _, eig2, uni, _ = _bands.get_emesh(Nx, Ny, Nz, ham_r, S_r, rvec, avec.T, sw_uni=True)
        wlist = np.linspace(eig.min()-mu, eig.max()-mu, Nw, True)   # energies relative to mu
        dos = flibs.gen_dos(eig2, uni, mu, wlist, dos_delta).sum(axis=0)
        dos_mu = np.interp(0.0, wlist, dos)
        if dos_mu <= 1.0e-8*dos.max():
            raise ValueError("tau_mode='dos': the DOS vanishes at mu, so tau(e) = tau_const*D(mu)/D(e) "
                             "is undefined (insulator?)")
        # the broadened DOS is positive at every band energy; the floor only guards round-off
        dos_e = np.maximum(np.interp(eig-mu, wlist, dos), 1.0e-12*dos.max())
        return (tau_const*dos_mu/dos_e,
                f'DOS-ratio tau: tau = tau_const*D(mu)/D(e), tau(mu) = {tau_const} '
                + units.tau_label)
    if tau_mode == 'flex':
        raise ValueError("tau_mode='flex' renormalizes the bands as well as tau; "
                         "use flex_quasiparticles (the Boltzmann drivers do this) instead of build_tau")
    raise ValueError(f"unknown tau_mode '{tau_mode}' (use 'const', 'e2', 'dos', 'epa' or 'flex')")


def _fold_index(rvec_s: np.ndarray, Nx: int, Ny: int, Nz: int) -> np.ndarray:
    '''Flat index of R mod (Nx, Ny, Nz) on the C-ordered mesh of gen_klist(sw_pp=False).'''
    return ((rvec_s[:, 0] % Nx)*Ny+rvec_s[:, 1] % Ny)*Nz+rvec_s[:, 2] % Nz


def _sigma_k(sig_w: np.ndarray, idx: np.ndarray, Nx: int, Ny: int, Nz: int) -> np.ndarray:
    '''
    Sigma(k) = sum_R Sigma(R) e^{-2 pi i k.R} [Nk, Norb, Norb] at the points of the (Nx, Ny, Nz)
    mesh, for ONE frequency sig_w [Norb, Norb, Nr]: at k = m/N the phase depends on R only through
    R mod N (idx, from _fold_index), so folding R onto the mesh and one FFT give that sum exactly.
    '''
    Norb, _, Nr = sig_w.shape
    buf = np.zeros((Nx*Ny*Nz, Norb*Norb), dtype=np.complex128)
    np.add.at(buf, idx, sig_w.reshape(Norb*Norb, Nr).T)
    return np.fft.fftn(buf.reshape(Nx, Ny, Nz, Norb*Norb), axes=(0, 1, 2)).reshape(-1, Norb, Norb)


def _band_diag_linear(sig_m: np.ndarray, mesh: tuple, klist: np.ndarray, uni: np.ndarray) -> np.ndarray:
    '''
    Sigma_nn(k, e) [Nk, Norb, Ne] at the points klist [Nk, 3] (fractional) from Sigma(k_m, e)
    [Nm, Norb, Norb, Ne] on the periodic mesh `mesh` (gen_klist order), multilinear in k and
    projected with the eigenvectors uni (as _sigma_band_diag).  A convex combination of the
    mesh values, so what holds at every mesh point -- a negative semidefinite anti-Hermitian
    part (causality), dRe Sigma/de < 1 (Z > 0) -- holds everywhere; a Fourier interpolation of
    the real-axis Sigma^R, much sharper in k than Sigma(iw_n), rings and breaks both.
    '''
    m = np.array(mesh)
    x = np.mod(klist, 1.0)*m
    i0 = np.floor(x).astype(np.int64) % m
    t = np.where(m == 1, 0.0, x-np.floor(x))   # a single-point direction is constant
    corners = []
    for c in np.ndindex(2, 2, 2):
        c = np.array(c)
        if np.any((c == 1) & (m == 1)):
            continue
        wgt = np.prod(np.where(c == 1, t, 1.0-t), axis=1)
        i = (i0+c) % m
        corners.append((wgt, (i[:, 0]*m[1]+i[:, 1])*m[2]+i[:, 2]))
    Nk, Norb = uni.shape[:2]
    snn = np.empty((Nk, Norb, sig_m.shape[3]), dtype=np.complex128)
    for e in range(sig_m.shape[3]):   # one energy at a time: Sigma(k) is Nk*Norb^2 per energy
        sk = sum(wgt[:, None, None]*sig_m[idx, :, :, e] for wgt, idx in corners)
        snn[:, :, e] = np.einsum('kna,kab,knb->kn', uni, sk, uni.conj(), optimize=True)
    return snn


def _sigma_band_diag(sig_r: np.ndarray, rvec_s: np.ndarray, uni: np.ndarray, Nx: int, Ny: int,
                     Nz: int) -> np.ndarray:
    '''
    Band-diagonal self-energy Sigma_nn(k, w) = <n k|Sigma(k, w)|n k> [Nk, Norb, Nw] on the
    periodic mesh of _lattice.gen_klist(Nx, Ny, Nz, sw_pp=False), from the Wannier-R form
    sig_r [Norb, Norb, Nw, Nr] (self_en_wannier.npz, or its real-axis continuation).  sig_r is
    in the same orbital layout as gen_ham's H(k) (it goes straight into gen_green with it), so
    the eigenvectors of that H(k) project it directly.  uni[k] holds them as conjugated rows
    (U = uni[k]^H).  Sigma(k) = sum_R Sigma(R) e^{-2 pi i k.R} (as spectrum_on_path).
    '''
    Norb, _, Nw, Nr = sig_r.shape
    idx = _fold_index(rvec_s, Nx, Ny, Nz)
    snn = np.empty((Nx*Ny*Nz, Norb, Nw), dtype=np.complex128)
    for n in range(Nw):   # one frequency at a time: Sigma(k) itself is Nk*Norb^2 per frequency
        sk = _sigma_k(sig_r[:, :, n, :], idx, Nx, Ny, Nz)
        snn[:, :, n] = np.einsum('kna,kab,knb->kn', uni, sk, uni.conj(), optimize=True)
    return snn


def _self_mesh(npz: np.lib.npyio.NpzFile, fname: str) -> tuple:
    '''
    The k-mesh FLEX solved Sigma on.  Exports carry it ('mesh'); for older files it is the
    extent of the stored R (exact unless zero_tol dropped the outermost shell, in which case
    the mesh comes out coarser and its points are no longer FLEX's own).
    '''
    if 'mesh' in npz.files:
        return tuple(int(n) for n in npz['mesh'])
    r = np.asarray(npz['rvec'], dtype=np.int64)
    mesh = tuple(int(n) for n in r.max(axis=0)-r.min(axis=0)+1)
    print(f"{fname}: no 'mesh' stamp, FLEX k-mesh taken as {mesh[0]}x{mesh[1]}x{mesh[2]} from the R "
          'extent; rerun FLEX (sw_out_self=True) to store it', flush=True)
    return mesh


def _repair_causal(s: np.ndarray, zl: np.ndarray, tol: float) -> np.ndarray:
    '''
    s [Norb, Norb, Ne] on the real energies zl (zl[:3] = -h, 0, h, then an ascending grid)
    made causal: on the grid, every run of energies whose anti-Hermitian part has a positive
    eigenvalue (> tol) is replaced by the linear interpolation between the causal energies
    bounding it (the nearest one at a grid edge); -h, 0, h are then re-read from the grid if
    any of them was acausal or sits in a repaired run.  The run is first widened over the flanks of the dip, while the
    margin -lambda_max grows outward: its zero crossings would pin Gamma ~ 0 into the
    interpolation.  No causal grid energy: eigenvalues clipped to <= 0 instead.
    '''
    lmax = np.linalg.eigvalsh(np.moveaxis((s-np.conj(s.swapaxes(0, 1)))/2j, 2, 0)).max(axis=1)
    e = zl.real
    good = lmax[3:] <= tol
    if not good.any():
        w, v = np.linalg.eigh(np.moveaxis((s-np.conj(s.swapaxes(0, 1)))/2j, 2, 0))
        a = np.einsum('eij,ej,ekj->eik', v, np.minimum(w, 0.0), v.conj())
        return (s+np.conj(s.swapaxes(0, 1)))/2.0+1j*np.moveaxis(a, 0, 2)
    g = s[:, :, 3:].copy()
    eg, mg = e[3:], -lmax[3:]
    bad = ~good
    for j in np.where(np.diff(bad.astype(np.int8)) == 1)[0]:     # good j, bad j+1: left flank
        while j > 0 and mg[j-1] > mg[j]:
            bad[j] = True
            j -= 1
    for j in np.where(np.diff(bad.astype(np.int8)) == -1)[0]+1:  # bad j-1, good j: right flank
        while j < len(mg)-1 and mg[j+1] > mg[j]:
            bad[j] = True
            j += 1
    good = ~bad
    if not good.any():                      # the flanks swallowed the grid: keep the causal points
        good = lmax[3:] <= tol
    ig = np.where(good)[0]
    for j in np.where(~good)[0]:
        lo, hi = ig[ig < j], ig[ig > j]
        if len(lo) and len(hi):
            t = (eg[j]-eg[lo[-1]])/(eg[hi[0]]-eg[lo[-1]])
            g[:, :, j] = (1.0-t)*s[:, :, 3+lo[-1]]+t*s[:, :, 3+hi[0]]
        else:
            g[:, :, j] = s[:, :, 3+(lo[-1] if len(lo) else hi[0])]
    out = np.concatenate([s[:, :, :3], g], axis=2)
    jb = np.clip(np.searchsorted(eg, e[:3])-1, 0, len(eg)-2)   # grid interval holding -h, 0, h
    if (lmax[:3] > tol).any() or not (good[jb].all() and good[jb+1].all()):
        t = (e[:3]-eg[jb])/(eg[jb+1]-eg[jb])
        out[:, :, :3] = (1.0-t)*g[:, :, jb]+t*g[:, :, jb+1]
    return out


def _causal_retarded(sig_c: np.ndarray, iwl: np.ndarray, zl: np.ndarray, uni_c: np.ndarray, near: np.ndarray,
                     h: float, orders: np.ndarray, block: int) -> tuple:
    '''
    Sigma^R(k, e) [Nk, Norb, Norb, Ne] on the real energies zl (zl[:3] = -h, 0, h) by a Pade
    ensemble of the orbital matrix elements, at the k of FLEX's own mesh (sig_c [Nk, Norb, Norb,
    Niw] is the solver's data there, not an interpolation, so the continuation sees no
    interpolation noise).  An order is accepted at k when it is CAUSAL: the anti-Hermitian part
    (Sigma - Sigma^H)/2i is negative semidefinite at every e, and Z > 0 on the bands near the
    Fermi level (near [Nk, Norb]; Z from the slope of Re Sigma_nn over +-h; Z > 1 is allowed:
    away from the Fermi surface -Im Sigma need not be smallest at e = 0).  The element-wise
    median of the accepted orders is used, their mean (a convex combination, causal again)
    when the median is not.  k-points where no order is accepted (orders agreeing on a narrow
    acausal dip are the typical case, far from the Fermi surface) take the plain median and
    REPAIR it: on each energy interval where it is acausal the matrix is replaced by the linear
    interpolation between the causal energies at the ends (a convex combination, so causal);
    with no causal energy at all, the eigenvalues of its anti-Hermitian part are clipped to <= 0.
    Low orders are not tried there: they are causal but too smooth, and would put a vanishing
    Gamma at large |e| into the k-interpolation.
    Returns (Sigma^R, number of accepted orders [Nk] (0: repaired), spread of
    -Im Sigma_nn(0) over the accepted orders [Nk]).
    '''
    Nk, Norb, _, _ = sig_c.shape
    Ne = len(zl)
    out = np.empty((Nk, Norb, Norb, Ne), dtype=np.complex128)
    nok = np.zeros(Nk, dtype=np.int64)
    spread = np.zeros(Nk)
    tol = 1.0e-8*np.abs(sig_c[..., 0]).max()
    antiherm = lambda s: (s-np.conj(s.swapaxes(1, 2)))/2j

    def causal(s, u, nr):   # s [nk, Norb, Norb, Ne] -> bool [nk]
        lmax = np.linalg.eigvalsh(np.moveaxis(antiherm(s), 3, 1)).max(axis=(1, 2))
        snn = np.einsum('kna,kabe,knb->kne', u, s[..., :3], u.conj())
        slope = (snn[..., 2].real-snn[..., 0].real)/(2.0*h)
        return np.isfinite(s).all(axis=(1, 2, 3)) & (lmax <= tol) & ((slope < 1.0) | ~nr).all(axis=1)

    def ensemble(kk, ords):   # continuations of the k-points kk with the orders ords
        nb = len(kk)
        data = sig_c[kk].reshape(nb*Norb*Norb, -1)
        zz = np.tile(zl, (nb*Norb*Norb, 1))
        trial = np.empty((len(ords), nb, Norb, Norb, Ne), dtype=np.complex128)
        ok = np.zeros((len(ords), nb), dtype=bool)
        for i, n in enumerate(ords):
            trial[i] = flibs.pade_batch(data[:, :n], iwl[:n], zz).reshape(nb, Norb, Norb, Ne)
            ok[i] = causal(trial[i], uni_c[kk], near[kk])
        with np.errstate(invalid='ignore', divide='ignore'), warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)   # all-NaN slices: k without an accepted order
            raw = np.nanmedian(np.where(np.isfinite(trial), trial.real, np.nan), axis=0) \
                + 1j*np.nanmedian(np.where(np.isfinite(trial), trial.imag, np.nan), axis=0)
            trial[~ok] = complex(np.nan, np.nan)   # a bare np.nan would leave Im = 0 in the medians
            med = np.nanmedian(trial.real, axis=0)+1j*np.nanmedian(trial.imag, axis=0)
            avg = np.nanmean(trial, axis=0)
            g0 = -np.einsum('kna,okab,knb->okn', uni_c[kk], trial[..., 1], uni_c[kk].conj()).imag
            sp = (np.nanmax(g0, axis=0)-np.nanmin(g0, axis=0))/np.maximum(np.abs(np.nanmedian(g0, axis=0)), 1e-300)
            sp = np.nan_to_num(np.nanmax(np.where(near[kk], sp, np.nan), axis=1))   # widest band near E_F
        res = np.where(causal(np.nan_to_num(med), uni_c[kk], near[kk])[:, None, None, None], med, avg)
        return res, ok.sum(axis=0), sp, raw

    allk = np.arange(Nk)
    raw = np.empty_like(out)
    for k0 in range(0, Nk, block):
        kk = allk[k0:k0+block]
        out[kk], nok[kk], spread[kk], raw[kk] = ensemble(kk, orders)
    for k in allk[nok == 0]:
        out[k] = _repair_causal(np.nan_to_num(raw[k]), zl, tol)
    return out, nok, spread


def flex_quasiparticles(eig: np.ndarray, Nx: int, Ny: int, Nz: int, ham_r: np.ndarray, S_r: np.ndarray,
                        rvec: np.ndarray, temp: float, units: _units.Units, gamma_imp: float = 0.0,
                        n_pade: int = 40, ewin: float = 20.0, de: float = 0.25,
                        fname: str = 'self_en_wannier.npz') -> dict:
    '''
    Quasiparticle Boltzmann parameters (tau_mode='flex') from the FLEX self-energy.

    Sigma(k, iw_n) is continued to the real axis ON FLEX'S OWN k-mesh, where it is the solver's
    data rather than a Fourier interpolation (interpolated data put spurious near-axis Pade
    poles on ~10% of the states, whatever the arithmetic precision): the orbital matrix
    Sigma^R(k, e) is built on an energy grid by a causal Pade ensemble (_causal_retarded),
    interpolated multilinearly in k onto the transport mesh (which keeps it causal, unlike a
    Fourier interpolation of the sharp real-axis Sigma^R) and projected on the bands.  Around
    the Fermi level the pole of G_nn = 1/(e - xi - Sigma_nn(e)), xi = E - mu_FLEX, gives
        Z = 1/(1 - dRe Sigma^R/de|_0),   e* = Z (xi + Re Sigma^R(0)),
        Gamma = -Im Sigma^R(e*) + gamma_imp,   tau* = hbar/(2 Z Gamma),
    and the bubble (bare current vertex v) in this limit equals Boltzmann transport with
    energies e*, velocities Z v, inverse mass Z m^-1 and lifetime tau*:
    sigma ~ (Z v)^2 tau* = Z v^2 hbar/(2 Gamma) for each state.  Gamma is taken at e*, not
    at 0, so its particle-hole asymmetry -- what the Seebeck coefficient feeds on -- is kept.
    No current-vertex correction (d Sigma/dk, backflow): that is the Kubo bubble, not the
    conserving response.

    eig: band energies of the caller's mesh [Nk, Norb] (checked against this mesh's own,
    so the returned arrays line up with the caller's velocities).
    gamma_imp: extra k-independent damping [eV] (impurity scattering, Matthiessen's rule);
    n_pade: highest Pade order (Matsubara points); orders n_pade/2 .. n_pade are combined.
    ewin, de: the real-energy grid is -ewin*T .. ewin*T in steps of de*T; e* is clipped to it
    (states farther out carry no -df/de weight) and Gamma(e*) is interpolated linearly on it.
    Returns a dict: eig (= mu + e*), Z, gamma, tau [tau_unit] (all [Nk, Norb]), mu (= the
    FLEX chemical potential, the Fermi level of e*) and a label.
    '''
    if len(S_r) != 0:
        raise ValueError("tau_mode='flex' needs an orthogonal (non-MLO) basis, as FLEX does")
    try:
        npz = np.load(fname)
    except FileNotFoundError:
        raise FileNotFoundError(f"tau_mode='flex' needs '{fname}' (run FLEX with sw_out_self=True)") from None
    sig_r, rvec_s, iw = npz['sigma'], np.asarray(npz['rvec'], dtype=np.int64), npz['iw']
    mu_self, temp_self = float(npz['mu']), float(npz['temp'])
    # Sigma lives on the Matsubara mesh of ITS temperature; transport at another T would
    # silently mix two temperatures (the Fermi window here, the scattering there)
    if not np.isclose(temp_self, temp, rtol=1.0e-6, atol=0.0):
        raise ValueError(f"'{fname}' was written at T = {temp_self:.6e} eV but transport runs at "
                         f"T = {temp:.6e} eV: rerun FLEX at this temperature")
    Norb = ham_r.shape[1]
    if sig_r.shape[0] != Norb:
        raise ValueError(f"'{fname}' has Norb = {sig_r.shape[0]}, the Hamiltonian {Norb}")
    Nk, klist = _lattice.gen_klist(Nx, Ny, Nz, sw_pp=False)
    eig0, uni = _bands.get_eigs(klist, ham_r, S_r, rvec)
    if eig0.shape != eig.shape or not np.allclose(eig0, eig):
        raise ValueError('flex_quasiparticles: the caller\'s band energies are not on the '
                         'gen_klist(Nx, Ny, Nz, sw_pp=False) mesh')
    Npd = min(n_pade, len(iw))
    if Npd < 2:
        raise ValueError(f"tau_mode='flex' needs at least 2 Matsubara points (n_pade = {n_pade})")
    iwl = 1j*iw[:Npd]
    # --- Sigma^R on FLEX's mesh ------------------------------------------------
    # energies: -h, 0, h (Z and Re Sigma(0) from a central difference, h << T) and the grid
    # Gamma(e*) is read from
    h = 0.1*temp
    Ng = int(round(ewin/de))
    egrid = de*temp*np.arange(-Ng, Ng+1)
    zl = np.concatenate([[-h, 0.0, h], egrid]).astype(np.complex128)
    mx, my, mz = _self_mesh(npz, fname)
    idx = _fold_index(rvec_s, mx, my, mz)
    sig_c = np.stack([_sigma_k(sig_r[:, :, n, :], idx, mx, my, mz) for n in range(Npd)], axis=-1)
    _, klist_c = _lattice.gen_klist(mx, my, mz, sw_pp=False)
    eig_c, uni_c = _bands.get_eigs(klist_c, ham_r, S_r, rvec)
    near = np.abs(eig_c-mu_self) < ewin*temp
    orders = np.unique(np.linspace(max(2, Npd//2), Npd, 8).astype(np.int64))
    block = max(1, int(2.5e8//(16*Norb*Norb*len(zl)*(len(orders)+3))))   # ~250 MB of trial continuations
    sig_ret, nok, spread = _causal_retarded(sig_c, iwl, zl, uni_c, near, h, orders, block)
    # --- interpolate onto the transport mesh, project on the bands ---------------
    snn = _band_diag_linear(sig_ret, (mx, my, mz), klist, uni).reshape(Nk*Norb, len(zl))
    xi = eig.ravel()-mu_self
    resig0 = snn[:, 1].real
    with np.errstate(divide='ignore'):
        zr = 1.0/(1.0-(snn[:, 2].real-snn[:, 0].real)/(2.0*h))
    # Z > 1 (allowed away from the Fermi surface, see _causal_retarded): clipped to (0, 1]
    zbad = ~((zr > 0.0) & (zr <= 1.0))
    Z = np.where(zbad, np.clip(np.nan_to_num(zr, nan=1.0, posinf=1.0, neginf=1.0), 1.0e-3, 1.0), zr)
    estar = Z*(xi+resig0)
    x = np.clip(estar, egrid[0], egrid[-1])
    j = np.clip(np.floor((x-egrid[0])/(de*temp)).astype(np.int64), 0, 2*Ng-1)
    t = (x-egrid[j])/(de*temp)
    rows = np.arange(Nk*Norb)
    gam = -((1.0-t)*snn[rows, 3+j].imag+t*snn[rows, 4+j].imag)
    gfloor = 1.0e-3*temp
    gtot = gam+gamma_imp
    gam = np.maximum(gtot, gfloor)
    tau = units.hbar/(2.0*Z*gam)/units.tau_unit
    # --- report, over the states that carry the Fermi window ----------------
    w = 0.25*(1.-np.tanh(0.5*estar/temp)**2)/temp
    sel = w > 0.01*w.max()
    wa = lambda x: (w*x).sum()/w.sum()
    print(f'FLEX quasiparticles from {fname}: mu_FLEX = {mu_self:.4f} eV, Pade on the FLEX mesh '
          f'{mx}x{my}x{mz}, orders {orders.min()}..{orders.max()} ({len(orders)}), real-energy grid '
          f'+-{ewin:g}T / {de:g}T, gamma_imp = {gamma_imp:.3e} eV', flush=True)
    print(f'  Fermi-window averages: <Z> = {wa(Z):.4f}  <Gamma> = {wa(gam):.4e} eV  '
          f'<Re Sigma(0)> = {wa(resig0):.4e} eV', flush=True)
    print(f'  Z range [{Z[sel].min():.4f}, {Z[sel].max():.4f}],  Gamma range '
          f'[{gam[sel].min():.4e}, {gam[sel].max():.4e}] eV over {sel.sum()} states', flush=True)
    nfb = (nok == 0).sum()
    print(f'  Pade causal in {nok.mean()/len(orders)*100:.1f}% of the orders (FLEX-mesh average); '
          f'<spread of Gamma(0) over orders / Gamma> = {spread[nok > 0].mean() if (nok > 0).any() else 0.0:.3e}',
          flush=True)
    nz, nfl = zbad[sel].sum(), (gtot < gfloor)[sel].sum()
    if nfb or nz or nfl:
        print(f'  WARNING: {nfb} of {len(nok)} FLEX k-points had no causal Pade order (acausal energies of '
              f'the median repaired by interpolation); in the Fermi window {nz} states had Z outside (0, 1] (clipped) '
              f'and {nfl} sit on the Gamma floor {gfloor:.1e} eV (add gamma_imp if these dominate)', flush=True)
    return {'eig': (mu_self+estar).reshape(Nk, Norb), 'Z': Z.reshape(Nk, Norb),
            'gamma': gam.reshape(Nk, Norb), 'tau': tau.reshape(Nk, Norb), 'mu': mu_self,
            'label': f'FLEX quasiparticles: tau* = hbar/(2 Z Gamma(e*)), Pade on the FLEX mesh, orders '
                     f'{orders.min()}..{orders.max()}'
                     + (f', gamma_imp = {gamma_imp:.3e} eV' if gamma_imp else '')}


def _carrier_state(eig: np.ndarray, vk: np.ndarray, imass: np.ndarray | None, Nk: int, kweight: np.ndarray,
                   fill: float, temp: float, tau_const: float, units: _units.Units, tau_kw: dict,
                   rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray, avec: np.ndarray, Nx: int, Ny: int,
                   Nz: int, Nw: int) -> tuple:
    '''
    (eig, vk, imass, mu, tau, tau_label) the Boltzmann kernels run on.  For every tau_mode
    but 'flex' these are the bare bands, mu fixed by fill, and build_tau; 'flex' replaces
    them with the FLEX quasiparticles (renormalized energies and velocities, the FLEX mu).
    tau_kw['flex_kw'] holds the keyword arguments of flex_quasiparticles.
    '''
    tau_kw = dict(tau_kw)
    flex_kw = tau_kw.pop('flex_kw', None) or {}
    if tau_kw.get('tau_mode') == 'flex':
        qp = flex_quasiparticles(eig, Nx, Ny, Nz, ham_r, S_r, rvec, temp, units, **flex_kw)
        # Luttinger check: the quasiparticle Fermi volume should hold the electrons FLEX
        # fixed mu_FLEX with; a large mismatch means e* = Z(xi + Re Sigma(0)) is not a good
        # description of this Fermi surface (Gamma ~ band scale, or a poor continuation)
        nqp = (flibs.get_ffermi(qp['eig'], qp['mu'], temp).sum(axis=1)*kweight).sum()/kweight.sum()
        print(f'  quasiparticle filling = {nqp:.4f}  (fill = {fill:.4f})', flush=True)
        Zv = qp['Z'][:, :, None]
        return (qp['eig'], vk*Zv, None if imass is None else imass*Zv[..., None], qp['mu'], qp['tau'],
                qp['label'])
    mu = _bands.calc_mu(eig, Nk, fill, temp, kweight)
    tau, tau_label = build_tau(units=units, eig=eig, mu=mu, temp=temp, tau_const=tau_const, rvec=rvec, ham_r=ham_r,
                               S_r=S_r, avec=avec, Nx=Nx, Ny=Ny, Nz=Nz, Nw=Nw, **tau_kw)
    return eig, vk, imass, mu, tau, tau_label


def print_tau_summary(tau: np.ndarray, eig: np.ndarray, mu: float, temp: float, label: str, indent: str = '  ',
                      tau_label: str = 'fs') -> None:
    '''Report tau where it matters: weighted by the Fermi window -df/de, not over all bands.'''
    print(f'{indent}tau: {label}', flush=True)
    w = 0.25*(1.-np.tanh(0.5*(eig-mu)/temp)**2)/temp
    if w.sum() <= 0:
        return
    tav = (w*tau).sum()/w.sum()
    sel = w > 0.01*w.max()          # states actually inside the Fermi window
    unit = tau_label
    print(f'{indent}     <tau>_FS = {tav:.4g} {unit}   range over the Fermi window '
          f'[{tau[sel].min():.4g}, {tau[sel].max():.4g}] {unit}', flush=True)
    if tau.std() > 0:
        # R_H picks up <tau^2>/<tau>^2, which is 1 only for a constant tau
        r = ((w*tau*tau).sum()/w.sum())/tav**2
        print(f'{indent}     <tau^2>/<tau>^2 = {r:.4f} '
              + ('(constant tau: R_H is tau-independent)' if abs(r-1.) < 1e-6
                 else '(>1: tau(k) shifts R_H away from the constant-tau value)'), flush=True)


def _hall_at_mesh(rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray, avec: np.ndarray, Nx: int, Ny: int, Nz: int,
                  fill: float, temp: float, tau_const: float, gsp: float, units: _units.Units, tau_kw: dict,
                  Nw: int = 300, detail: bool = False) -> dict:
    '''
    Weak-field (B-linear) Boltzmann Hall response on ONE k-mesh.
    detail=True additionally builds the band-resolved (two-fluid) decomposition
    and the Luttinger-volume carrier densities; that part is only worth its cost
    on the mesh actually being reported, not on every step of a convergence sweep.
    Returns a dict; all printing is done by get_hall_coe.
    '''
    Nk, eig, vk, imass, kweight = _bands.get_emesh(Nx, Ny, Nz, ham_r, S_r, rvec, avec.T*units.ihbar,
                                                   sw_veloc=True, sw_mass=True)
    Vuc = sclin.det(avec)*1e-30
    wsum = kweight.sum()          # NOT Nk: stays correct if a symmetry-reduced mesh is ever used
    eig, vk, imass, mu, tau, tau_label = _carrier_state(eig, vk, imass, Nk, kweight, fill, temp, tau_const, units,
                                                        tau_kw, rvec, ham_r, S_r, avec, Nx, Ny, Nz, Nw)
    K0, K1, K2 = flibs.calc_Kn(eig, vk, kweight, temp, mu, tau)
    # Full Hall kernel S[n,a,b,g] for the three field directions g, band resolved.
    # The tensor (rather than the single B||z scalar) is what allows the proper
    # rho^(1) = -sigma^-1 sigma^(1) sigma^-1 below instead of assuming sigma is
    # diagonal in xyz.
    Shall_n, Sabs_n = flibs.calc_sigma_hall(eig, vk, imass/units.eC, kweight, tau, temp, mu)
    Shall = Shall_n.sum(axis=0)
    Rh = _calc.hall_coefficient(K0, Shall, Vuc, wsum, gsp)          # m^3/C, per field axis
    nH, kind = _calc.hall_carrier_density(Rh, units.eC)                 # cm^-3, 'electron'/'hole'
    sigma = gsp*units.tau_unit*units.eC*K0/(wsum*Vuc)                       # S/m
    # sigma^(1)/B carries one more power of e and of tau than sigma  -> S/(m T)
    sigma1 = gsp*units.eC*units.eC*units.tau_unit*units.tau_unit*Shall/(wsum*Vuc)
    # Hall mobility mu_H = R_H * sigma_aa with the current in the plane normal to B.
    # omega_c*tau = mu_H[SI]*B, so mu_H directly measures how far the weak-field
    # expansion this whole routine is built on can be trusted.
    muH = np.array([Rh[g]*sigma[(g+1)%3, (g+1)%3] for g in range(3)])   # m^2/(V s)
    # --- two k-mesh adequacy diagnostics -----------------------------------
    dfermi = 0.25*(1.-np.tanh(0.5*(eig-mu)/temp)**2)/temp
    # (i) effective number of states carrying the Fermi window, as the
    #     participation ratio of the -df/de weights.  Everything here is an
    #     average over Neff states, so Neff sets the statistical noise floor.
    Neff = dfermi.sum()**2/(dfermi**2).sum()
    # (ii) cancellation ratio of the Hall integrand: it changes sign around the
    #     Fermi surface, and when |sum| << sum|.| the total is a small residue of
    #     large opposing parts -- exactly the regime where R_H keeps changing
    #     sign with the mesh.  sigma_xx has no such cancellation, which is why it
    #     looks converged long before sigma^(1) is.
    cancel = abs(Shall[0, 1, 2])/max(Sabs_n[:, 0, 1, 2].sum(), 1e-300)
    out = {'mesh': (Nx, Ny, Nz), 'Nk': Nk, 'mu': mu, 'Rh': Rh, 'nH': nH, 'kind': kind,
           'sigma': sigma, 'sigma1': sigma1, 'muH': muH, 'Neff': Neff, 'cancel': cancel, 'K0': K0,
           'tau_label': tau_label}
    if detail:
        # keep the per-state arrays only for the mesh that gets reported: holding
        # eig/tau for every step of a sweep is hundreds of MB on a fine mesh
        out['tau'], out['eig'] = tau, eig
        # --- two-fluid decomposition -----------------------------------------
        # sigma and sigma^(1) are both plain band sums, so each sheet can be given
        # its own conductivity, Hall coefficient and mobility.  Comparing a band's
        # OWN n_H against its Luttinger volume shows whether that sheet behaves as
        # one simple closed pocket; the mismatch between the total n_H and n_e/n_h
        # is then the compensation the experiment is really seeing.
        K0n = flibs.calc_k0_band(eig, vk, kweight, temp, mu, tau)                   # [Norb,3,3]
        Norb = len(K0n)
        trK0 = np.trace(K0)
        Rh_n = np.full(Norb, np.nan)
        for n in range(Norb):
            if np.trace(K0n[n]) > 1e-8*trK0:
                Rh_n[n] = _calc.hall_coefficient(K0n[n], Shall_n[n], Vuc, wsum, gsp)[2]
        out['car'] = _calc.fs_carrier_density(eig, mu, Vuc, gsp, mesh=(Nx, Ny, Nz))
        out['sigma_n'] = gsp*units.tau_unit*units.eC*K0n/(wsum*Vuc)                        # [Norb,3,3] S/m
        out['Rh_n'] = Rh_n
    return out


def get_hall_coe(rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray, avec: np.ndarray, Nx: int, Ny: int, Nz: int,
                 fill: float, temp: float, tau_const: float | None, units: _units.Units, tau_kw: dict, Nw: int = 300,
                 with_spin: bool = False, mesh_list: list | None = None) -> list | None:
    '''
    Weak-field Hall coefficient R_H and the Hall carrier density n_H=1/(e|R_H|).

    R_H is formed from the FULL tensors, rho^(1)=-sigma^-1 sigma^(1) sigma^-1,
    for all three field directions.  Give mesh_list (e.g. [20,30,40] or
    [(20,20,10),...]) to run a k-mesh convergence sweep: sigma^(1) converges far
    more slowly than sigma and can come out with the WRONG SIGN on an
    under-converged mesh, so a single-mesh number is not usable on its own.
    tau_kw: keyword arguments of build_tau selecting the relaxation-time model
    (tau_mode, epa_file, dos_delta; flex_kw for tau_mode='flex', see flex_quasiparticles).
    tau_const: [tau_unit], used by tau_mode 'const', 'e2' and 'dos'; may be None for 'epa'/'flex'.
    '''
    # Parameter validation
    if temp <= 0:
        print("Error: Temperature (temp) is non-positive", flush=True)
        return
    if not _tau_const_ok(tau_const, tau_kw):
        return
    if sclin.det(avec) <= 0:
        print("Error: Unit cell volume (Vuc) is non-positive", flush=True)
        return
    gsp = (1.0 if with_spin else 2.0)  # spin degeneracy factor
    if mesh_list:
        meshes = [(m, m, m) if np.isscalar(m) else tuple(m) for m in mesh_list]
    else:
        meshes = [(Nx, Ny, Nz)]
    print(f"T = {temp/units.kb:.3f} K", flush=True)
    print(f"spin degeneracy gsp = {gsp:.0f} "
          f"({'SOC: both spins already in H' if with_spin else 'spin-degenerate bands'})", flush=True)
    print('  mesh          Nk      mu[eV]     Neff   cancel      R_H^x        R_H^y        R_H^z    [m^3/C]', flush=True)
    results = []
    for i, (nx, ny, nz) in enumerate(meshes):
        r = _hall_at_mesh(rvec, ham_r, S_r, avec, nx, ny, nz, fill, temp, tau_const, gsp, units, tau_kw, Nw,
                          detail=(i == len(meshes)-1))
        results.append(r)
        print(f'  {nx:3d}x{ny:3d}x{nz:3d} {r["Nk"]:9d}  {r["mu"]:9.4f} {r["Neff"]:8.3g} {r["cancel"]:8.2e} '
              + ' '.join(f'{v:+12.4e}' for v in r['Rh']), flush=True)
    res = results[-1]
    # --- k-mesh convergence verdict -----------------------------------------
    if len(results) > 1:
        prev = results[-2]['Rh']; last = res['Rh']
        print('  k-mesh convergence (finest vs previous):', flush=True)
        for g, ax in enumerate('xyz'):
            if abs(last[g]) < 1e-30:
                continue
            if prev[g]*last[g] < 0:
                msg = 'SIGN FLIP - not converged'
            else:
                d = abs(last[g]/prev[g]-1.)*100
                msg = f'{d:6.2f}% change'+('  <- not converged' if d > 5. else '')
            print(f'    B||{ax}: {msg}', flush=True)
    else:
        print('  (single mesh: pass mesh_list=[...] to check k-mesh convergence)', flush=True)
    # --- detail on the finest mesh ------------------------------------------
    print(f'--- {res["mesh"][0]}x{res["mesh"][1]}x{res["mesh"][2]} mesh, mu = {res["mu"]:.4f} eV ---', flush=True)
    print_tau_summary(res['tau'], res['eig'], res['mu'], temp, res['tau_label'], tau_label=units.tau_label)
    print_matrix('sigma matrix (S/m)', res['sigma'])
    for g, ax in enumerate('xyz'):
        print(f'  B||{ax}: R_H = {res["Rh"][g]:+.5e} m^3/C   '
              f'n_H = {res["nH"][g]:.5e} cm^-3 ({res["kind"][g]})   '
              f'mu_H = {1e4*res["muH"][g]:+.4e} cm^2/Vs   '
              f'omega_c*tau = {abs(res["muH"][g]):.3e} /T', flush=True)
    print(f'  sigma^(1)_xy/B = {res["sigma1"][0, 1, 2]:+.5e} S/(m T)  (B||z)', flush=True)
    # --- validity warnings ---------------------------------------------------
    print(f'  Fermi-window sampling: Neff = {res["Neff"]:.4g} states, '
          f'Hall-integrand cancellation |sum|/sum|.| = {res["cancel"]:.3e}', flush=True)
    if res['Neff'] < 1.0e3:
        print(f'  WARNING: only ~{res["Neff"]:.4g} states carry the Fermi window. '
              f'Refine the k-mesh (or raise T) before trusting R_H.', flush=True)
    if res['cancel'] < 1.0e-2:
        print(f'  WARNING: the Hall integrand cancels to {res["cancel"]:.2e} of its absolute weight. '
              f'R_H is a small residue of large opposing contributions here (near-compensation), '
              f'so it needs a much finer mesh than sigma_xx and its SIGN may still be mesh-dependent.', flush=True)
    Bmax = 0.1/max(abs(res['muH']).max(), 1e-30)
    print(f'  weak-field limit omega_c*tau<0.1 holds up to B ~ {Bmax:.3g} T', flush=True)
    # --- band-resolved (two-fluid) decomposition vs the Luttinger volume ------
    car = res['car']
    print('--- carrier densities: Luttinger volume vs Hall, band by band ---', flush=True)
    print('  band  type      n_e[cm^-3]    n_h[cm^-3]   pk e/h   sigma_xx[S/m]   R_H(band)'
          '    n_H(band)   mu_H[cm^2/Vs]', flush=True)
    for n, k in enumerate(car['kind']):
        if k in ('full', 'empty'):
            continue
        rn = res['Rh_n'][n]
        if np.isnan(rn) or abs(rn) < 1e-30:
            hall = f'{"-":>12s} {"-":>12s} {"-":>13s}'
        else:
            hall = (f'{rn:+12.4e} {1./(abs(rn)*units.eC)/1e6:12.4e}'
                    f' {1e4*abs(rn*res["sigma_n"][n, 0, 0]):13.4e}')
        print(f'  {n+1:4d}  {k:8s}  {car["ne_band"][n]:.5e}   {car["nh_band"][n]:.5e}  '
              f'{car["npocket_e"][n]:2d}/{car["npocket_h"][n]:<2d}   '
              f'{res["sigma_n"][n, 0, 0]:.5e}  {hall}', flush=True)
    ne, nh = car['n_e'], car['n_h']
    print(f'  Luttinger total: n_e = {ne:.5e}   n_h = {nh:.5e} cm^-3'
          + (f'   n_e/n_h = {ne/nh:.4f}' if nh > 0 else ''), flush=True)
    print(f'  Hall (all bands, B||z): n_H = {res["nH"][2]:.5e} cm^-3 ({res["kind"][2]})'
          '   <- the number to compare with experiment', flush=True)
    if (car['kind'] == 'open').any():
        print('  note: an OPEN sheet percolates the zone and has no Luttinger volume, so the '
              'n_e/n_h totals above miss whatever it carries.', flush=True)
    if ne > 0 and nh > 0 and abs(ne/nh-1.) < 0.2:
        print('  note: nearly compensated (n_e ~ n_h). n_H is then a small difference of '
              'large terms, extremely sensitive to mu, and is NOT a carrier density -- '
              'compare mobilities/two-fluid fits with experiment instead.', flush=True)
    return results


def calc_conductivity_Boltzmann(rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray, avec: np.ndarray, Nx: int,
                                Ny: int, Nz: int, fill: float, temp: float, tau_const: float | None,
                                units: _units.Units, tau_kw: dict, Nw: int = 300, with_spin: bool = False,
                                sw_tdf: bool = False) -> dict | None:
    '''
    calculate conductivities using Boltzmann equations
    tau_kw: keyword arguments of build_tau (tau_mode, epa_file, dos_delta), plus flex_kw
    (keyword arguments of flex_quasiparticles) for tau_mode='flex'
    tau_const: [tau_unit], used by tau_mode 'const', 'e2' and 'dos'; may be None for 'epa'/'flex'
    sw_tdf: compute the transport distribution function first, then energy-integrate
    Returns a dict of the tensors (sigma, kappa, kappa2, sigmaS, Seebeck, Peltier, Lorenz,
    Lorenz2, PF, mu); with sw_tdf also the transport distribution 'tdf' [Nw,3,3] on the
    energies 'tdf_e' (relative to mu).  None if the input is invalid.
    '''
    # no dep. T and mu or filling
    # Parameter validation
    if temp <= 0:
        print("Error: Temperature (temp) is non-positive", flush=True)
        return
    if not _tau_const_ok(tau_const, tau_kw):
        return

    Nk, eig, vk, kweight = _bands.get_emesh(Nx, Ny, Nz, ham_r, S_r, rvec, avec.T*units.ihbar, sw_veloc=True)
    Vuc = sclin.det(avec)*1e-30
    if Vuc <= 0:
        print("Error: Unit cell volume (Vuc) is non-positive", flush=True)
        return
    gsp = (1.0 if with_spin else 2.0)  # spin weight
    # calc_Kn / calc_tdf return kweight-weighted sums, so normalize by the total weight (as the
    # Hall kernels do), not by Nk: identical on the uniform mesh, correct on a weighted one
    wsum = kweight.sum()
    iNV = 1./(wsum*Vuc)
    itemp = 1./temp
    eig, vk, _, mu, tau, tau_label = _carrier_state(eig, vk, None, Nk, kweight, fill, temp, tau_const, units,
                                                    tau_kw, rvec, ham_r, S_r, avec, Nx, Ny, Nz, Nw)
    print(f"T = {temp/units.kb:.3f} K", flush=True)
    print(f"mu = {mu:.4f} eV", flush=True)
    print_tau_summary(tau, eig, mu, temp, tau_label, indent='', tau_label=units.tau_label)
    if sw_tdf:
        tdf = flibs.calc_tdf(eig, vk, kweight, tau, Nw)
        # absolute-energy grid matching the calc_tdf bin CENTRES: E_i = emin + (i-1/2)*dw
        dw = (eig.max()-eig.min())/Nw
        wlist = eig.min()+dw*(np.arange(1, Nw+1)-0.5)
        dfermi = 0.25*(1.-np.tanh(0.5*(wlist-mu)/temp)**2)/temp
        K0 = (dfermi*tdf.T).T.sum(axis=0)*dw
        K1 = (dfermi*(wlist-mu)*tdf.T).T.sum(axis=0)*dw
        K2 = (dfermi*(wlist-mu)**2*tdf.T).T.sum(axis=0)*dw
    else:
        K0, K1, K2 = flibs.calc_Kn(eig, vk, kweight, temp, mu, tau)
    # sigma [S/m] = gsp * tau[s] * e[C] * K0[eV·m^-3] / (sum_k w_k * V[m^3])
    sigma = gsp*units.tau_unit*units.eC*K0*iNV  # e cancels eV->J so units work out to S/m
    kappa = gsp*units.tau_unit*units.eC*units.kb*K2*iNV*itemp

    # inverses on the conducting directions only (see _inv_conducting): a direction without
    # conduction gets S = Peltier = L = 0 instead of discarding every component
    K0_inv, dead = _inv_conducting(K0)
    if dead:
        print(f"Note: no conduction along {dead} (sigma = 0 there, e.g. a 2D model or Nz = 1); "
              "Seebeck, Peltier and Lorenz are given on the conducting directions, 0 along it", flush=True)
    kappa2 = gsp*units.tau_unit*units.eC*units.kb*(K2-K1.dot(K0_inv.dot(K1)))*iNV*itemp
    Seebeck = -units.kb*K0_inv.dot(K1)*itemp
    Peltier = K1.dot(K0_inv)
    sigma_inv = _inv_conducting(sigma*temp)[0]
    # Lorenz tensor L = kb * kappa . (sigma*T)^-1 (matrix product, not element-wise)
    Lorenz = units.kb*kappa.dot(sigma_inv)
    Lorenz2 = units.kb*kappa2.dot(sigma_inv)
    sigmaS = gsp*units.tau_unit*units.kb*units.eC*K1*iNV*itemp
    PF = sigma*Seebeck**2
    print_matrix('sigma matrix (S/m)', sigma)
    print_matrix('kappa matrix (K22 only) (W/m/K)', kappa)
    print_matrix('kappa matrix (full) (W/m/K)', kappa2)
    print_matrix('sigmaS matrix (A/m/K)', sigmaS)
    print_matrix('Seebeck matrix (V/K)', Seebeck)
    print_matrix('Peltier matrix (V)', Peltier)
    print_matrix('Lorenz matrix (K22 only) (Wohm/K^2) (fe 2.44e-8)', Lorenz, width=9, prec=2)
    print_matrix('Lorenz matrix? (full) (Wohm/K^2)', Lorenz2, width=9, prec=2)
    print_matrix('Power Factor (SA/m^2/K)', PF, ndigits=8)
    out = {'mu': mu, 'sigma': sigma, 'kappa': kappa, 'kappa2': kappa2, 'sigmaS': sigmaS, 'Seebeck': Seebeck,
           'Peltier': Peltier, 'Lorenz': Lorenz, 'Lorenz2': Lorenz2, 'PF': PF}
    if sw_tdf:
        out['tdf'], out['tdf_e'] = tdf, wlist-mu
    return out


def calc_conductivity_lrt(rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray, avec: np.ndarray, Nx: int, Ny: int,
                          Nz: int, fill: float, temp: float, Nw: int, delta: float, units: _units.Units, Emax: float,
                          color_thick: float = 1.0e-3, color_scatt: float = 1.0e-6,
                          with_spin: bool = False) -> dict | None:
    '''
    calculation of linear response theory
    electric conductivity of LRT corresponds to Boltzmann then delta~O(10-1) (tau~1fs) at 300K
    thermal conductivity and L12 of LRT correspond to Boltzmann then delta~O(10-3) (tau~100fs) at 300K
    Emax: upper end of the frequency window [eV]; color_thick / color_scatt: slab thickness and
    Kubelka-Munk scattering length [m] of the non-metal colours
    Returns a dict for the figures (libs.plots.plot_optical): wlist, sigma [Nw,3,3] (S/m),
    epsilon, refl [Nw,3], wp_scr (screened plasma frequency or None), and the perceived
    colours cols / ncols.  None if the input is invalid.
    '''
    # Parameter validation
    if temp <= 0:
        print("Error: Temperature (temp) is non-positive", flush=True)
        return
    if delta <= 0:
        print("Error: Broadening (delta) is non-positive", flush=True)
        return

    Nk, eig, vk, kweight = _bands.get_emesh(Nx, Ny, Nz, ham_r, S_r, rvec, avec.T*units.ihbar, True, True)
    Vuc = sclin.det(avec)
    if Vuc <= 0:
        print("Error: Unit cell volume (Vuc) is non-positive", flush=True)
        return
    gsp = (1.0 if with_spin else 2.0)  # spin weight
    mu = _bands.calc_mu(eig, Nk, fill, temp, kweight)
    print(f'chemical potential = {mu:.4f} eV', flush=True)
    print(f'temperature = {temp/units.kb:.3f} K', flush=True)
    print(f'about tau = {1.e15*units.hbar/delta:.3f} fs', flush=True)
    print(f'delta = {delta:9.3e}', flush=True)
    L11, L12, L22, wlist = _response.get_conductivity(mu, temp, eig, vk, Nw, Emax, delta, kweight)
    sigmaconst = gsp*units.hbar*units.eC/Vuc*1.0e30
    kappaSconst = sigmaconst*units.kb/temp
    sigma = sigmaconst*L11
    kappa = kappaSconst*L22
    sigmaS = kappaSconst*L12
    # full electronic thermal conductivity subtracts the thermoelectric backflow:
    # kappa = kappaSconst*(L22 - L12.L11^-1.L12)  (cf. Boltzmann K2-K1.K0^-1.K1).
    # The subtraction is done on the bare L tensors (all share the same prefactor),
    # so kappa2 carries the same kappaSconst as the L22-only kappa.
    # inverses on the conducting directions only (see _inv_conducting), frequency by frequency
    inv_w = [_inv_conducting(l11) for l11 in L11]
    dead = sorted({d for _, d in inv_w if d})
    if dead:
        print(f"Note: no conduction along {'; '.join(dead)} at some frequencies (e.g. a 2D model or Nz = 1); "
              "Seebeck, kappa backflow and Lorenz are given on the conducting directions", flush=True)
    kappa2 = kappaSconst*np.array([l22-l12.dot(li.dot(l12)) for (li, _), l12, l22 in zip(inv_w, L12, L22)])
    # sigma = sigmaconst*L11, so sigma^-1 = L11^-1/sigmaconst on the same directions
    Seebeck = np.array([-(li/sigmaconst).dot(sS) for (li, _), sS in zip(inv_w, sigmaS)])
    print_matrix('sigma matrix (S/m)', sigma[0].real)
    print_matrix('kappa matrix (L22 only) (W/m/K)', kappa[0].real)
    print_matrix('kappa matrix (full) (W/m/K)', kappa2[0].real)
    print_matrix('sigmaS matrix (A/m/K)', sigmaS[0].real)
    print('Lorenz number (Wohm/K^2) (fe 2.44e-8)', flush=True)
    # Lorenz tensor L = kb * kappa . (sigma*T)^-1 (matrix product, not element-wise)
    sigmaT_inv = _inv_conducting(sigma[0]*temp)[0]
    Lorenz_matrix = (units.kb*kappa[0].dot(sigmaT_inv)).real.round(10)
    Lorenz_matrix2 = (units.kb*kappa2[0].dot(sigmaT_inv)).real.round(10)
    print(' (L22 only):', flush=True)
    for lor in Lorenz_matrix:
        print(f" {lor[0]:9.2e} {lor[1]:9.2e} {lor[2]:9.2e}", flush=True)
    print(' (full):', flush=True)
    for lor in Lorenz_matrix2:
        print(f" {lor[0]:9.2e} {lor[1]:9.2e} {lor[2]:9.2e}", flush=True)
    print_matrix('Seebeck coefficient matrix (V/K)', Seebeck[0].real)
    # dielectric function (relative permittivity) from the optical conductivity,
    # SI convention: eps(w) = 1 + i*sigma(w)/(eps0*w_rad) = 1 + i*hbar*sigma/(eps0*w),
    # sigma in S/m, w in eV (hbar converts the eV grid to angular frequency).
    # Frequency-resolved like sigma; the w=0 point is skipped (eps ~ 1/w there).
    eps0 = scconst.epsilon_0
    epsilon = 1.+1j*units.hbar/eps0*(sigma.T/(wlist+1e-8)).T
    # console summary at the lowest positive frequency (static-limit estimate)
    print(f'dielectric function at w = {wlist[1]:9.3e} eV', flush=True)
    print_matrix('Re[eps] matrix', epsilon[1].real)
    print_matrix('Im[eps] matrix', epsilon[1].imag)
    # --- optical constants: N, R, alpha, and the non-metal spectra -----------
    # N(w) = n + i*kappa = sqrt(eps(w)) on the passive branch (Im N >= 0), then
    # the normal-incidence Fresnel reflectivity R = |(N-1)/(N+1)|^2 and the
    # Beer-Lambert absorption coefficient alpha = 2*w*kappa/(hbar*c). Only the
    # principal (diagonal) components are used, so these describe light
    # polarized along a crystal axis; w=0 is skipped as for eps.
    Nref = _optical.refractive_index(epsilon.diagonal(axis1=1, axis2=2))  # [Nw,3] xx,yy,zz
    refl = _optical.reflectivity(Nref)
    alpha = _optical.absorption_coefficient(wlist, Nref)                  # [1/m]
    # what leaves a NON-metal: the slab transmittance at thickness color_thick
    # and the Kubelka-Munk diffuse (body) reflectance at scattering length
    # color_scatt. A metal is opaque, so both come out black -- correctly so.
    trans = _optical.transmittance(refl, alpha, color_thick)
    rinf = _optical.kubelka_munk(alpha, color_scatt)
    absorp = np.clip(1.0-refl-trans, 0.0, 1.0)          # absorptance of the slab
    print(f'reflectivity (normal incidence) at w = {wlist[1]:9.3e} eV', flush=True)
    print(f'  R = ({refl[1, 0]:.4f}, {refl[1, 1]:.4f}, {refl[1, 2]:.4f})  (xx, yy, zz)', flush=True)
    print(f'absorption coefficient at w = {wlist[1]:9.3e} eV', flush=True)
    print(f'  alpha = ({alpha[1, 0]:9.3e}, {alpha[1, 1]:9.3e}, {alpha[1, 2]:9.3e}) 1/m', flush=True)
    # --- perceived colour (metallic lustre) from R(w) ------------------------
    # R(w) is folded with the CIE 1931 2-deg observer under illuminant D65 and
    # converted to sRGB (see plibs._color), the standard route for first-
    # principles metal colours (Prandini et al., npj Comput. Mater. 5, 129
    # (2019)). The polycrystalline entry uses the directional average of R,
    # which is what an unpolarised measurement on a polycrystal sees.
    cols = [(lab, _color.spectrum_color(wlist[1:], r))
            for lab, r in [('xx', refl[1:, 0]), ('yy', refl[1:, 1]), ('zz', refl[1:, 2]),
                           ('poly', refl[1:].mean(axis=1))]]
    # --- non-metal colours: the light that got through the material ----------
    # A dielectric reflects almost nothing, so its colour is carried by the
    # transmitted (or diffusely back-scattered) light, both computed from the
    # directional average. The absorptance entry is the complement of what
    # leaves the sample, so it is reported as its photographic negative.
    ncols = [('T(d)', _color.spectrum_color(wlist[1:], trans[1:].mean(axis=1)), 'hex'),
             ('KM',  _color.spectrum_color(wlist[1:], rinf[1:].mean(axis=1)), 'hex'),
             ('A neg', _color.spectrum_color(wlist[1:], absorp[1:].mean(axis=1)), 'hex_neg')]
    print('perceived colour (CIE 1931 2-deg observer, D65, sRGB):', flush=True)
    # show the colour itself as a 24-bit ANSI block, but only on a real terminal
    # so redirected logs stay free of escape sequences
    tty = os.isatty(1)
    for lab, c in cols:
        r, g, b = (int(v) for v in c['rgb255'])
        sq = f'\x1b[48;2;{r};{g};{b}m    \x1b[0m ' if tty else ''
        print(f'  {lab:4s} {sq}{c["hex"]}  RGB = ({r:3d},{g:3d},{b:3d})'
              f'  xy = ({c["xy"][0]:.4f}, {c["xy"][1]:.4f})  Y = {c["Y"]:.4f}'
              + ('  [out of sRGB gamut, clipped]' if c['clipped'] else ''), flush=True)
    print(f'non-metal colour (d = {1e3*color_thick:.4g} mm, '
          f'l_scatt = {1e6*color_scatt:.4g} um, directional average):', flush=True)
    for lab, c, key in ncols:
        r, g, b = (int(v) for v in (c['rgb255_neg'] if key == 'hex_neg' else c['rgb255']))
        sq = f'\x1b[48;2;{r};{g};{b}m    \x1b[0m ' if tty else ''
        note = (' (opaque at this thickness)' if lab == 'T(d)' and c['Y'] < 1.e-4 else
                ' (negated: A is the complement of the transmitted light)'
                if key == 'hex_neg' else '')
        # Y belongs to the spectrum that was integrated, so on the negated row
        # it is labelled Y(A): the luminance of the absorptance, not of the swatch
        print(f'  {lab:5s} {sq}{c[key]}  RGB = ({r:3d},{g:3d},{b:3d})'
              f'  {"Y(A)" if key == "hex_neg" else "Y"} = {c["Y"]:.4f}{note}', flush=True)
    if not cols[0][1]['full_range']:
        print(f'  Warning: the spectrum does not cover the visible window '
              f'(380-780 nm = 1.59-3.26 eV); R was held constant outside '
              f'[{_color.ev_nm/wlist[-1]:.0f}, {_color.ev_nm/wlist[1]:.0f}] nm, so the colour is '
              f'only indicative. Increase Emax (and Nw) to fix this.', flush=True)
    # --- plasma frequency ----------------------------------------------------
    # (1) effective (unscreened) plasma frequency from the optical f-sum rule
    #     (hbar*wp_a)^2 = (2*hbar/(pi*eps0)) * int_0^Emax Re sigma_aa(w) dw
    #     (sigma in S/m, w in eV -> hbar*wp in eV). It is the spectral weight up
    #     to the cutoff Emax and, unlike a Drude fit, needs no broadening input.
    hwp = np.sqrt([abs(2*units.hbar/(np.pi*eps0)*np.trapezoid(sigma[:, i, i].real, wlist)) for i in range(3)])
    h_eVs = scconst.physical_constants['Planck constant in eV s'][0]  # h = 2*pi*hbar
    fp = hwp/h_eVs/1e12   # plasma frequency nu_p = wp/(2*pi) in THz
    print(f'plasma frequency (f-sum rule, integrated to Emax={Emax} eV):', flush=True)
    print(f'  hbar*wp = ({hwp[0]:.4f}, {hwp[1]:.4f}, {hwp[2]:.4f}) eV  (xx, yy, zz)', flush=True)
    print(f'  nu_p    = ({fp[0]:.2f}, {fp[1]:.2f}, {fp[2]:.2f}) THz', flush=True)
    # (2) screened plasma frequency: energy where Re eps_xx crosses zero upward
    #     (the longitudinal plasmon / reflectivity edge), by linear interpolation.
    e1, wpos = epsilon[1:, 0, 0].real, wlist[1:]
    zc = np.where((e1[:-1] < 0.0)&(e1[1:] >= 0.0))[0]
    if zc.size:
        j = zc[0]
        wp_scr = wpos[j]-e1[j]*(wpos[j+1]-wpos[j])/(e1[j+1]-e1[j])
        print(f'  screened plasma freq (Re eps_xx=0): hbar*wp_scr = {wp_scr:.4f} eV'
              f'  (nu_scr = {wp_scr/h_eVs/1e12:.2f} THz)', flush=True)
    else:
        wp_scr = None
        print('  screened plasma freq (Re eps_xx=0): no upward zero crossing in [0, Emax]', flush=True)
    # write sigma(w) and eps(w) (diagonal components) to a data file
    try:
        header = ('perceived colour (CIE1931 2-deg, D65, sRGB) '
                  + '  '.join(f'{lab}: {c["hex"]} RGB({c["rgb255"][0]},{c["rgb255"][1]},'
                              f'{c["rgb255"][2]})' for lab, c in cols)+'\n'
                  + f'non-metal colour (d = {color_thick:.4g} m, l_scatt = {color_scatt:.4g} m) '
                  + '  '.join(f'{lab}: {c[key]}' for lab, c, key in ncols)+'\n'
                  + f'plasma frequency (f-sum rule) hbar*wp [eV] = '
                  f'{hwp[0]:.6f} {hwp[1]:.6f} {hwp[2]:.6f}  (xx yy zz)\n'
                  + (f'screened plasma frequency (Re eps_xx=0) hbar*wp_scr [eV] = {wp_scr:.6f}\n'
                     if wp_scr is not None else 'screened plasma frequency (Re eps_xx=0): none in window\n')
                  + 'w[eV] '
                  'Re_sigma_xx Im_sigma_xx Re_sigma_yy Im_sigma_yy Re_sigma_zz Im_sigma_zz '
                  'Re_eps_xx Im_eps_xx Re_eps_yy Im_eps_yy Re_eps_zz Im_eps_zz '
                  'R_xx R_yy R_zz alpha_xx alpha_yy alpha_zz[1/m] '
                  'T_xx T_yy T_zz Rinf_xx Rinf_yy Rinf_zz')
        d = np.arange(3)
        data = np.column_stack([wlist[1:],
                                *[f(sigma[1:, i, i]) for i in d for f in (np.real, np.imag)],
                                *[f(epsilon[1:, i, i]) for i in d for f in (np.real, np.imag)],
                                refl[1:], alpha[1:], trans[1:], rinf[1:]])
        np.savetxt('optical.dat', data, header=header)
        print("optical conductivity, dielectric function and reflectivity written to 'optical.dat'", flush=True)
    except OSError as e:
        print(f"Error: Failed to write 'optical.dat': {e}", flush=True)
    return {'wlist': wlist, 'sigma': sigma, 'epsilon': epsilon, 'refl': refl, 'wp_scr': wp_scr,
            'cols': cols, 'ncols': ncols}
