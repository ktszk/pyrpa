#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
Tests for the impurity supercell and its Bloch-unfolded spectral function
(libs/plibs/_impurity.py: impurity_supercell, unfolded_spectrum).

  * a pristine supercell must reproduce the primitive-cell bands at every commensurate k,
    with unit weight per band -- checked on a model whose bands are NOT even in k
    (complex hopping), which fixes the sign of the unfolding phase;
  * the unfolded weights obey the sum rule sum_n W_n(k) = Norb at any k, also with
    impurities present.

Runs standalone (no pytest needed):  python tests/test_impurity.py
Also works under pytest if installed:  pytest tests/test_impurity.py
Requires the Fortran library libfmod.so (cd libs && make FC=ifx SL=MKL).
"""
import os
import sys
import numpy as np
import scipy.linalg as sclin

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import libs.flibs as F
import libs.plibs as P


def _chiral_square(t: float = 1.0, phi: float = 0.7) -> tuple:
    """1-orbital square lattice with a complex x hopping: e(k) = -2t cos(2 pi k_x - phi) - 2t cos(2 pi k_y)."""
    rvec = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]], float)
    ham_r = np.zeros((5, 1, 1), complex)
    ham_r[1, 0, 0] = -t*np.exp(1j*phi)
    ham_r[2, 0, 0] = -t*np.exp(-1j*phi)
    ham_r[3:, 0, 0] = -t
    return rvec, ham_r


def test_pristine_supercell_unfolds_to_the_bands():
    rvec, ham_r = _chiral_square()
    H, rlist = P.impurity_supercell(ham_r, rvec, [8, 8, 1], np.array([], dtype=np.int64), 0.0)
    eigs, uni = sclin.eigh(H)
    ks = np.array([[0.25, 0.0, 0], [-0.25, 0.0, 0], [0.375, 0.125, 0], [0.125, -0.375, 0]])
    ek = F.get_eig(F.gen_ham(ks, ham_r, rvec))[0][:, 0]
    assert np.abs(ek[0]-ek[1]) > 1.0                                  # e(k) != e(-k): the test can see the sign
    w = np.linspace(-5, 5, 4001)
    A = P.unfolded_spectrum(eigs, uni, rlist, 1, ks, w, 0.0, 0.01)
    for e, a in zip(ek, A):
        assert abs(w[a.argmax()]-e) < 2*(w[1]-w[0])
    # the whole weight sits at e(k): the eigenstates at that energy (a degenerate set, when
    # other folded k share it) carry weight 1 together, all others none
    phase = np.exp(2j*np.pi*ks.dot(rlist.T))
    W = np.abs(phase.dot(uni))**2/len(rlist)
    for e, Wk in zip(ek, W):
        at_e = np.abs(eigs-e) < 1e-8
        assert np.isclose(Wk[at_e].sum(), 1.0) and np.allclose(Wk[~at_e], 0.0, atol=1e-10)


def test_unfolded_weight_sum_rule_with_impurities():
    rvec, ham_r = _chiral_square()
    H, rlist = P.impurity_supercell(ham_r, rvec, [6, 6, 1], np.array([0, 14, 21]), 1.5)
    eigs, uni = sclin.eigh(H)
    ks = np.random.default_rng(0).random((5, 3))
    ks[:, 2] = 0.0
    W = np.abs(np.exp(2j*np.pi*ks.dot(rlist.T)).dot(uni))**2/len(rlist)
    assert np.allclose(W.sum(axis=1), 1.0, atol=1e-12)                # Norb = 1


def test_missing_onsite_block_is_added():
    """A hopping file without an R=0 entry (zero on-site energies) still takes an impurity."""
    rvec, ham_r = _chiral_square()
    H0, _ = P.impurity_supercell(ham_r[1:], rvec[1:], [4, 4, 1], np.array([3]), 0.7)
    H1, _ = P.impurity_supercell(ham_r, rvec, [4, 4, 1], np.array([3]), 0.7)
    assert np.allclose(H0, H1)


def _pristine_bands_match(ham_r, rvec, cell, ks):
    """eigenvalues of the pristine supercell contain e_n(k) for every commensurate k in ks"""
    H, _ = P.impurity_supercell(ham_r, rvec, cell, np.array([], dtype=np.int64), 0.0)
    assert np.allclose(H, H.conj().T)
    eigs = sclin.eigvalsh(H)
    ek = F.get_eig(F.gen_ham(ks, ham_r, rvec))[0]
    for e in ek.ravel():
        assert np.min(np.abs(eigs-e)) < 1e-10


def test_axis_of_size_one_keeps_the_interlayer_hopping():
    """cell[2] = 1 is the kz = 0 cut: the Rz = +-1 hoppings must land in the on-site block."""
    rvec, ham_r = _chiral_square()
    rvec = np.vstack([rvec, [[0, 0, 1], [0, 0, -1], [1, 0, 1], [-1, 0, -1]]])
    ham_r = np.concatenate([ham_r, np.array([-0.4, -0.4, -0.15j, 0.15j]).reshape(4, 1, 1)])
    ks = np.array([[0.0, 0.0, 0.0], [0.25, 0.0, 0.0], [-0.25, 0.5, 0.0], [0.5, 0.25, 0.0]])
    _pristine_bands_match(ham_r, rvec, [4, 4, 1], ks)


def test_hopping_longer_than_half_the_cell_is_folded():
    """a third-neighbour hopping on a 4-site axis wraps onto the first neighbour exactly"""
    rvec, ham_r = _chiral_square()
    rvec = np.vstack([rvec, [[3, 0, 0], [-3, 0, 0], [2, 0, 0], [-2, 0, 0]]])
    ham_r = np.concatenate([ham_r, np.array([-0.2j, 0.2j, -0.1, -0.1]).reshape(4, 1, 1)])
    ks = np.array([[0.25, 0.0, 0.0], [-0.25, 0.25, 0.0], [0.5, 0.5, 0.0]])
    _pristine_bands_match(ham_r, rvec, [4, 4, 1], ks)


def test_per_orbital_potential_is_diagonal():
    """imp_V given as a list of Norb values is diag(imp_V), not added row-wise"""
    rvec = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0]], float)
    ham_r = np.zeros((3, 2, 2), complex)
    ham_r[0] = [[0.0, 0.3], [0.3, 0.5]]
    ham_r[1] = -np.eye(2)
    ham_r[2] = -np.eye(2)
    Hl, _ = P.impurity_supercell(ham_r, rvec, [4, 1, 1], np.array([1]), [1.0, 0.25])
    Hm, _ = P.impurity_supercell(ham_r, rvec, [4, 1, 1], np.array([1]), np.diag([1.0, 0.25]))
    H0, _ = P.impurity_supercell(ham_r, rvec, [4, 1, 1], np.array([], dtype=np.int64), 0.0)
    assert np.allclose(Hl, Hm)
    assert np.allclose((Hl-H0)[2:4, 2:4], np.diag([1.0, 0.25])) and np.isclose(np.abs(Hl-H0).sum(), 1.25)


def test_invalid_impurity_sites_are_rejected():
    rvec, ham_r = _chiral_square()
    for bad in ([0, 16], [-1], [3, 3]):
        try:
            P.impurity_supercell(ham_r, rvec, [4, 4, 1], np.array(bad), 1.0)
        except ValueError:
            continue
        raise AssertionError(f'imp_sites={bad} was accepted')


def test_nonorthogonal_unfolding_and_cpa():
    rvec = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0]], float)
    ham_r = np.array([[[0.1, 0.2j], [-0.2j, 0.7]],
                      [[-0.4, 0.1], [0.03j, -0.2]],
                      [[-0.4, -0.03j], [0.1, -0.2]]], complex)
    S_r = np.array([[[1, 0.15], [0.15, 1]],
                    [[0.08, 0.02j], [0.01, 0.04]],
                    [[0.08, 0.01], [-0.02j, 0.04]]], complex)
    sites = np.array([], dtype=int)
    H, rlist = P.impurity_supercell(ham_r, rvec, [4, 1, 1], sites, 0)
    S, _ = P.impurity_supercell(S_r, rvec, [4, 1, 1], sites, 0)
    eigs, uni = sclin.eigh(H, S)
    ks = np.array([[0, 0, 0], [0.25, 0, 0], [0.5, 0, 0]])
    hk, sk = F.gen_ham(ks, ham_r, rvec, Ovl_r=S_r)
    w, eta = np.linspace(-2, 2, 121), 0.07
    spectrum = P.unfolded_spectrum(eigs, uni, rlist, 2, ks, w, 0, eta, S)
    for i in range(len(ks)):
        bands = sclin.eigvalsh(hk[i], sk[i])
        expected = ((eta/np.pi)/((w[:, None]-bands)**2+eta**2)).sum(axis=1)
        assert np.allclose(spectrum[i], expected, atol=1e-10)

    H_imp, _ = P.impurity_supercell(ham_r, rvec, [4, 1, 1], np.array([1]), 0.5)
    eig_imp, u_imp = sclin.eigh(H_imp, S)
    k_any = np.array([[0.13, 0, 0], [0.39, 0, 0]])
    a = P.unfolded_spectrum(eig_imp, u_imp, rlist, 2, k_any, np.array([0.2]), 0, 1e5, S)
    assert np.allclose(a[:, 0]*np.pi*1e5, 2, atol=1e-8)

    va, vb = np.zeros((2, 2), complex), np.array([[0.4, 0.08j], [-0.08j, -0.2]])
    z = np.array([0.3+0.4j, -0.1+0.8j])
    sigma = F.solve_cpa(hk, va, vb, 0.4, z, Ovlk=sk)
    eye = np.eye(2)
    for iw, omega in enumerate(z):
        gloc = np.linalg.inv(omega*sk-hk-sigma[iw]).mean(axis=0)
        da, db = va-sigma[iw], vb-sigma[iw]
        residual = 0.4*da@np.linalg.inv(eye-gloc@da)+0.6*db@np.linalg.inv(eye-gloc@db)
        assert np.linalg.norm(residual) < 1e-9
    identity = np.broadcast_to(eye, hk.shape).copy()
    native = F.solve_cpa(hk, va, vb, 0.4, z)
    python = F.solve_cpa(hk, va, vb, 0.4, z, Ovlk=identity)
    assert np.allclose(native, python, atol=1e-9)


def test_nonorthogonal_cpa_host_spectrum_and_saved_self_energy(tmp_path):
    r = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0]], float)
    h = np.array([[[0.2]], [[-0.3]], [[-0.3]]], complex)
    s = np.array([[[1.]], [[0.1]], [[0.1]]], complex)
    path = [[0, 0, 0], [0.5, 0, 0]]
    mu, eta = 0.1, 0.04
    cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        _, _, w, A = P.calc_cpa(h, r, mu, 0.1, 4, 1, 1, 5, -0.5, 0.5, eta,
                                 path, 4, np.eye(3), x_cpa=0.0, VB_shift=1.0, S_r=s)
        assert np.allclose(np.load('self_en.npz')['arr_0'], 0)
    finally:
        os.chdir(cwd)
    klist, _, _ = P.mk_klist(path, 4, np.eye(3))
    hk, sk = F.gen_ham(klist, h, r, Ovl_r=s)
    reference = -np.imag(sk[:, 0, 0, None] /
                         ((w[None, :]+mu+1j*eta)*sk[:, 0, 0, None]-hk[:, 0, 0, None]))/np.pi
    assert np.allclose(A, reference, atol=1e-11)


if __name__ == '__main__':
    import _tools
    sys.exit(_tools.run_standalone(globals()))
