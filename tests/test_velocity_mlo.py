#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Regression / physics tests for the group velocity in a NON-ORTHOGONAL (MLO) basis
(get_veloc_mlo / get_vnm_mlo in libs/src/fmod.f90, wrapped by
libs/flibs/_hamiltonian.py get_vk / get_vnm / get_veloc / get_vnmk).

The MLO hoppings (ftype 3 = *.input MLO file, ftype 4 = ecalj job_mlo binaries)
carry an overlap matrix S(R), so the bands come from the GENERALIZED problem

    H(k) C = eps S(k) C ,     C^H S(k) C = 1        (get_eig_mlo)

whose Hellmann-Feynman velocity has an extra overlap-derivative term

    v_n = (C^H dH/dk C)_nn - eps_n (C^H dS/dk C)_nn .

Dropping that term (the plain get_veloc used for an orthogonal basis) gives a
velocity that is simply wrong -- O(100%) off for a real ecalj model.  The tests
below lock the correct expression in by comparing against a central finite
difference of the generalized eigenvalues themselves.

The models here are synthetic (random H(R), S(R) with the Hermiticity of a real
tight-binding model) so the tests need no input files.

Runs standalone (no pytest needed):  python tests/test_velocity_mlo.py
Also works under pytest if installed:  pytest tests/test_velocity_mlo.py
Requires the Fortran library libfmod.so (cd libs && make FC=ifx SL=MKL).
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import libs.flibs as F
import _tools as T

_mlo_model = T.mlo_model


# --------------------------------------------------------------------------- #
#  helpers
# --------------------------------------------------------------------------- #
def _eig_uni(klist, ham_r, S_r, rvec):
    hamk, Sk = F.gen_ham(klist, ham_r, rvec, Ovl_r=S_r)
    return F.get_eig(hamk, Sk)


def _fd_velocity(klist, ham_r, S_r, rvec, h=1e-5):
    """d eps_n / dk from a central difference.  get_vlm0 differentiates the phase
    exp(-2 pi i k.R), i.e. it returns d/d(2 pi k_frac), so the fractional-k
    finite difference is divided by 2 pi to match."""
    Nk, Norb = len(klist), ham_r.shape[1]
    fd = np.zeros((Nk, Norb, 3), dtype=np.float64)
    for a in range(3):
        dk = np.zeros(3)
        dk[a] = h
        ep = _eig_uni(np.ascontiguousarray(klist + dk), ham_r, S_r, rvec)[0]
        em = _eig_uni(np.ascontiguousarray(klist - dk), ham_r, S_r, rvec)[0]
        fd[:, :, a] = (ep - em) / (2.0 * h) / (2.0 * np.pi)
    return fd


def _klist(Nk=12, seed=7):
    rng = np.random.default_rng(seed)
    return np.ascontiguousarray(rng.random((Nk, 3)))


# --------------------------------------------------------------------------- #
#  tests
# --------------------------------------------------------------------------- #
def test_generalized_eigenvectors_are_S_orthonormal():
    """get_eig_mlo returns C with C^H S C = 1 -- the normalization get_veloc_mlo
    assumes when it subtracts eps * (C^H dS/dk C).  numpy sees the Fortran arrays
    transposed, so the Fortran C is uni[i].T and the check reads
    uni S uni^H = 1 on the numpy side (row b of uni[i] = eigenvector of band b)."""
    rvec, ham_r, S_r, Norb, _ = _mlo_model(seed=1)
    klist = _klist()
    _, Sk = F.gen_ham(klist, ham_r, rvec, Ovl_r=S_r)
    _, uni = _eig_uni(klist, ham_r, S_r, rvec)
    for i in range(len(klist)):
        gram = uni[i] @ Sk[i] @ uni[i].conj().T
        assert np.abs(gram - np.eye(Norb)).max() < 1e-10


def test_mlo_velocity_matches_finite_difference():
    """v_n = (C^H dH/dk C)_nn - eps_n (C^H dS/dk C)_nn reproduces d eps_n/dk."""
    rvec, ham_r, S_r, _, _ = _mlo_model(seed=2)
    klist = _klist()
    mrot = np.eye(3)
    eig, uni = _eig_uni(klist, ham_r, S_r, rvec)

    v = F.get_veloc(klist, ham_r, rvec, mrot, uni, S_r, eig)
    fd = _fd_velocity(klist, ham_r, S_r, rvec)
    assert np.abs(v - fd).max() < 1e-5 * max(1.0, np.abs(fd).max())


def test_overlap_term_is_not_negligible():
    """Guard against a silent regression to the orthogonal formula: without the
    dS/dk term the velocity misses the finite difference by a large margin."""
    rvec, ham_r, S_r, _, _ = _mlo_model(seed=2)
    klist = _klist()
    mrot = np.eye(3)
    eig, uni = _eig_uni(klist, ham_r, S_r, rvec)

    fd = _fd_velocity(klist, ham_r, S_r, rvec)
    v_no_ovl = F.get_veloc(klist, ham_r, rvec, mrot, uni)     # orthogonal formula
    assert np.abs(v_no_ovl - fd).max() > 0.05 * np.abs(fd).max()


def test_vnm_mlo_diagonal_equals_band_velocity():
    """The interband routine must reduce to get_veloc_mlo on the band diagonal:
    v_nn = H~_nn - (eps_n+eps_n)/2 S~_nn."""
    rvec, ham_r, S_r, Norb, _ = _mlo_model(seed=3)
    klist = _klist()
    mrot = np.eye(3)
    eig, uni = _eig_uni(klist, ham_r, S_r, rvec)

    v = F.get_veloc(klist, ham_r, rvec, mrot, uni, S_r, eig)
    vnm = F.get_vnmk(klist, ham_r, rvec, mrot, uni, S_r, eig)
    diag = np.einsum('knni->kni', vnm).real
    assert np.abs(diag - v).max() < 1e-12


def test_vnm_mlo_is_hermitian():
    """H~ and S~ are Hermitian and the (eps_n+eps_m)/2 prefactor is symmetric, so
    the interband velocity matrix stays Hermitian in the band indices."""
    rvec, ham_r, S_r, _, _ = _mlo_model(seed=4)
    klist = _klist()
    eig, uni = _eig_uni(klist, ham_r, S_r, rvec)
    vnm = F.get_vnmk(klist, ham_r, rvec, np.eye(3), uni, S_r, eig)
    assert np.abs(vnm - vnm.conj().transpose(0, 2, 1, 3)).max() < 1e-12


def test_axis_rotation_matches_orthogonal_convention():
    """mrot must act on the Cartesian axes exactly as in the orthogonal get_veloc.
    Fortran reads the C-ordered [3,3] array transposed, so vk(m) = sum_l mrot[l,m]*v(l),
    i.e. v_rot = v @ mrot -- the same for both routines."""
    rvec, ham_r, S_r, _, _ = _mlo_model(seed=5)
    klist = _klist()
    mrot = np.ascontiguousarray([[0.0, 1.0, 0.0],
                                 [2.0, 0.0, 0.0],
                                 [0.0, 0.0, 3.0]])
    eig, uni = _eig_uni(klist, ham_r, S_r, rvec)
    v1 = F.get_veloc(klist, ham_r, rvec, np.eye(3), uni, S_r, eig)
    v2 = F.get_veloc(klist, ham_r, rvec, mrot, uni, S_r, eig)
    assert np.abs(v2 - v1 @ mrot).max() < 1e-12

    # same convention on the orthogonal path, so ftype 4 stays consistent with the rest
    hamk = F.gen_ham(klist, ham_r, rvec)
    eo, uo = F.get_eig(hamk)
    w1 = F.get_veloc(klist, ham_r, rvec, np.eye(3), uo)
    w2 = F.get_veloc(klist, ham_r, rvec, mrot, uo)
    assert np.abs(w2 - w1 @ mrot).max() < 1e-12


def test_orthogonal_path_untouched():
    """S_r=[]/None must dispatch to the original get_veloc/get_vnm bit-for-bit."""
    rvec, ham_r, _, _, _ = _mlo_model(seed=6)
    klist = _klist()
    mrot = np.eye(3)
    hamk = F.gen_ham(klist, ham_r, rvec)
    eig, uni = F.get_eig(hamk)

    assert np.array_equal(F.get_veloc(klist, ham_r, rvec, mrot, uni),
                          F.get_veloc(klist, ham_r, rvec, mrot, uni, [], eig))
    assert np.array_equal(F.get_vnmk(klist, ham_r, rvec, mrot, uni),
                          F.get_vnmk(klist, ham_r, rvec, mrot, uni, None, eig))


def test_eig_required_with_overlap():
    """Passing sk0 without eig is a programming error, not a silent wrong answer."""
    rvec, ham_r, S_r, Norb, _ = _mlo_model(seed=8)
    klist = _klist(4)
    vk0 = F.get_vlm0(klist, ham_r, rvec)
    sk0 = F.get_vlm0(klist, S_r, rvec)
    _, uni = _eig_uni(klist, ham_r, S_r, rvec)
    for fn in (F.get_vk, F.get_vnm):
        try:
            fn(vk0, np.eye(3), uni, sk0)
        except ValueError:
            continue
        raise AssertionError(f'{fn.__name__} accepted sk0 without eig')


if __name__ == '__main__':
    import _tools
    sys.exit(_tools.run_standalone(globals()))
