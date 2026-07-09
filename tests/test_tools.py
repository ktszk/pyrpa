#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Self-tests of the shared test utilities (tests/_tools.py) -- these double as
usage examples.  Every reference implemented in _tools is checked against
either a closed form or the Fortran library itself, so a silent breakage of
the tooling cannot invalidate the rest of the suite unnoticed.

Runs standalone (no pytest needed):  python tests/test_tools.py
Also works under pytest if installed:  pytest tests/test_tools.py
Requires the Fortran library libfmod.so (cd libs && make FC=ifx SL=MKL).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import _tools as T
import libs.flibs as F


# --------------------------------------------------------------------------- #
#  distributions
# --------------------------------------------------------------------------- #
def test_fermi_bose_identities():
    temp = 0.03
    e = np.array([-50.0, -1.0, 0.0, 1.0, 50.0])
    f = T.fermi(e, temp)
    assert np.isfinite(f).all() and (f >= 0).all() and (f <= 1).all()
    # particle-hole symmetry f(e) + f(-e) = 1 and the T->0 step
    assert abs(T.fermi(0.7, temp) + T.fermi(-0.7, temp) - 1.0) < 1e-14
    assert T.fermi(50.0, temp) < 1e-300 or T.fermi(50.0, temp) == 0.0
    assert abs(T.fermi(0.0, temp, mu=0.0) - 0.5) < 1e-14
    # Bose: n(w) + n(-w) = -1, and n(w) ~ T/w - 1/2 for w << T
    w = np.array([0.5, 1.0, 2.0])
    assert np.abs(T.bose(w, temp) + T.bose(-w, temp) + 1.0).max() < 1e-12
    assert abs(T.bose(1e-6, temp) - (temp / 1e-6 - 0.5)) < 1e-3
    # Matsubara grids
    wl = T.fermionic_matsubara(temp, 4)
    assert np.allclose(wl, np.pi * temp * np.array([1, 3, 5, 7]))
    assert np.allclose(T.bosonic_matsubara(temp, 3),
                       2.0 * np.pi * temp * np.array([0, 1, 2]))


# --------------------------------------------------------------------------- #
#  model builders
# --------------------------------------------------------------------------- #
def test_one_orbital_bundle_green_is_free():
    """Gk of the bundle must equal 1/(iw + mu - eps) exactly."""
    md = T.one_orbital_square(Nx=4, Ny=4, Nw=8, temp=0.1, mu=0.2)
    wl = T.fermionic_matsubara(md['temp'], md['Nw'])
    eps = md['eig'][:, 0]
    Gref = 1.0 / (1j * wl[:, None] + md['mu'] - eps[None, :])
    assert np.abs(md['Gk'][0, 0] - Gref).max() < 1e-12
    # dispersion helper must agree with the reduced-grid eigenvalues at k=0
    assert abs(md['eps_full'][0, 0, 0] - (-2 * 1.0 - 2 * 0.5)) < 1e-12
    T.assert_green_causal(md['Gk'])


def test_two_orbital_bundle_matches_closed_form_bands():
    """Eigenvalues must be (e1+e2)/2 +- sqrt(((e1-e2)/2)^2 + V^2)."""
    md = T.two_orbital_square(Nx=4, Ny=4, Nw=4, t1=1.0, t2=0.5, de=0.4, V=0.2)
    for ik, k in enumerate(md['klist']):
        cx, cy = np.cos(2 * np.pi * k[0]), np.cos(2 * np.pi * k[1])
        e1 = -2.0 * md['t1'] * (cx + cy)
        e2 = -2.0 * md['t2'] * (cx + cy) + md['de']
        rt = np.sqrt(0.25 * (e1 - e2) ** 2 + md['V'] ** 2)
        ref = np.array([0.5 * (e1 + e2) - rt, 0.5 * (e1 + e2) + rt])
        assert np.abs(np.sort(md['eig'][ik]) - ref).max() < 1e-10
    T.assert_hermitian(md['hamk'], msg='two-orbital hamk')
    T.assert_green_causal(md['Gk'])


def test_band_diagonal_model_layout():
    eig, vk, veloc = T.band_diagonal_model(Nk=16, Norb=3, seed=1)
    assert eig.shape == (16, 3) and (np.diff(eig, axis=1) >= 0).all()
    # vk strictly band-diagonal and consistent with veloc
    off = vk.copy()
    for n in range(3):
        off[:, n, n, :] = 0.0
    assert np.abs(off).max() == 0.0
    assert np.abs(np.einsum('knni->kni', vk).real - veloc).max() < 1e-15


def test_cylindrical_fs_normalization():
    for gs in ('s', 'd', 'dxy'):
        wf, phi = T.cylindrical_fs(Nb=180, gap_sym=gs)
        # FS average <phi^2> = 1 (the normalization the solvers assume)
        assert abs((wf * phi ** 2).sum() / wf.sum() - 1.0) < 1e-12
    T.assert_raises(ValueError, T.cylindrical_fs, gap_sym='f')


# --------------------------------------------------------------------------- #
#  exact chi0 reference vs the Fortran library
# --------------------------------------------------------------------------- #
def test_two_pole_green_reduces_to_free():
    md = T.one_orbital_square(Nx=4, Ny=4, Nw=16, temp=0.05, mu=0.3)
    G = T.two_pole_green(md['eig'], md['temp'], md['mu'], md['Nw'], a=0.0)
    assert np.abs(G - md['Gk']).max() < 1e-12
    # finite a: still causal, and the two pole weights sum to 1
    E, W = T.two_pole_poles_weights(md['eig'][:, 0], a=0.5, delta=0.2)
    assert np.abs(W.sum(axis=0) - 1.0).max() < 1e-12
    Gi = T.two_pole_green(md['eig'], md['temp'], md['mu'], md['Nw'],
                          a=0.5, delta=0.2)
    T.assert_green_causal(Gi)


def test_exact_chi0_q0_sum_rule_and_tail_kernel():
    """exact_chi0_1orb against two independent checks: the closed-form
    q=0, nu=0 value sum_k f(1-f)/T, and the tail-corrected Fortran kernel
    (machine-exact at q=0, O(1/Nw^2) close elsewhere)."""
    md = T.one_orbital_square(Nx=8, Ny=8, Nw=256, temp=0.04, mu=0.3)
    mlist = [0, 1, 2]
    chi_ex = T.exact_chi0_1orb(md['eps_full'], md['klist'], mlist,
                               md['temp'], md['mu'])
    # closed form at q=0, nu=0
    f = T.fermi(md['eps_full'], md['temp'], md['mu'])
    ref00 = (f * (1.0 - f)).sum() / md['temp'] / md['eps_full'].size
    iq0 = int(np.where((np.abs(md['klist']) < 1e-12).all(axis=1))[0][0])
    assert abs(chi_ex[iq0, 0] - ref00) < 1e-12 * max(1.0, abs(ref00))
    # against the Fortran tail-corrected chi0 (O(1/Nw^2) truncation error)
    Gk = T.two_pole_green(md['eig'], md['temp'], md['mu'], md['Nw'])
    with T.silence():
        chi, _ = F.get_chi0_tail(md['Smat'], md['Cmat'], Gk, md['eig'],
                                 md['uni'], md['olist'], md['kmap'],
                                 md['invk'], md['mu'], md['temp'],
                                 md['Nx'], md['Ny'], md['Nz'])
    err = max(np.abs(chi[0, 0, m, :] - chi_ex[:, jm]).max()
              for jm, m in enumerate(mlist))
    assert err < 5e-3 * np.abs(chi_ex).max()


def test_exact_chi0_two_pole_against_sharp_cutoff_extrapolation():
    """With the interacting two-pole G the exact reference must be the
    Nw -> infinity limit of the plain sharp-cutoff convolution: the error
    halves when Nw doubles (O(1/Nw)) and fit_power_order sees slope ~1."""
    md = T.one_orbital_square(Nx=8, Ny=8, temp=0.04, mu=0.3)
    a, delta, mlist = 0.5, 0.2, [0, 1]
    chi_ex = T.exact_chi0_1orb(md['eps_full'], md['klist'], mlist,
                               md['temp'], md['mu'], a=a, delta=delta)
    nws, errs = [64, 128, 256], []
    for Nw in nws:
        Gk = T.two_pole_green(md['eig'], md['temp'], md['mu'], Nw,
                              a=a, delta=delta)
        with T.silence():
            chi, _ = F.get_chi0(md['Smat'], md['Cmat'], Gk, md['olist'],
                                md['kmap'], md['invk'], md['temp'],
                                md['Nx'], md['Ny'], md['Nz'])
        errs.append(max(np.abs(chi[0, 0, m, :] - chi_ex[:, jm]).max()
                        for jm, m in enumerate(mlist)))
    order = T.fit_power_order(nws, errs)
    assert 0.8 < order < 1.3, f"sharp-cutoff order {order:.2f}, expected ~1"


# --------------------------------------------------------------------------- #
#  numerics helpers
# --------------------------------------------------------------------------- #
def test_fit_power_order_recovers_synthetic_slopes():
    ns = np.array([16, 32, 64, 128, 256])
    for p in (0.5, 1.0, 2.0):
        assert abs(T.fit_power_order(ns, 3.7 / ns ** p) - p) < 1e-10
    assert abs(T.rel_err(np.array([1.0, 2.02]), np.array([1.0, 2.0])) - 0.01) < 1e-12


# --------------------------------------------------------------------------- #
#  physics assertions
# --------------------------------------------------------------------------- #
def test_assertions_detect_violations():
    h = np.array([[1.0, 0.5j], [-0.5j, 2.0]])
    T.assert_hermitian(h)
    T.assert_positive_semidefinite(h)
    T.assert_raises(AssertionError, T.assert_hermitian,
                    np.array([[0.0, 1.0], [0.0, 0.0]]))
    T.assert_raises(AssertionError, T.assert_positive_semidefinite,
                    np.array([[1.0, 0.0], [0.0, -1.0]]))
    # causality: flipping the sign of Im G must be caught
    md = T.one_orbital_square(Nx=2, Ny=2, Nw=4)
    T.assert_green_causal(md['Gk'])
    T.assert_raises(AssertionError, T.assert_green_causal, np.conj(md['Gk']))


def test_static_chi_of_stable_state_is_psd():
    """Usage example: the static RPA susceptibility matrices chis/chic(q) of
    a weakly interacting stable state must be Hermitian positive
    semi-definite in the orbital-pair basis."""
    md = T.two_orbital_square(Nx=4, Ny=4, Nw=64, temp=0.1, mu=0.0,
                              U=0.5, J=0.1)
    with T.silence():
        chi0, stoner = F.get_chi0(md['Smat'], md['Cmat'], md['Gk'],
                                  md['olist'], md['kmap'], md['invk'],
                                  md['temp'], md['Nx'], md['Ny'], md['Nz'])
        chis, chic = F.get_chis_chic(chi0, md['Smat'], md['Cmat'])
    assert 0.0 < stoner < 1.0                     # stable normal state
    for name, chi in (('chis', chis), ('chic', chic)):
        static = np.transpose(chi, (2, 0, 1))     # [Nk, Nchi, Nchi]
        T.assert_positive_semidefinite(static, tol=1e-10,
                                       msg=f'static {name}')


# --------------------------------------------------------------------------- #
#  standalone runner (no pytest required)
# --------------------------------------------------------------------------- #
if __name__ == '__main__':
    sys.exit(T.run_standalone(globals()))
