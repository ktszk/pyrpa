#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Regression / physics tests for the band-basis inverse effective mass tensor
(get_imassk / get_imassk_mlo in libs/src/fmod.f90, wrapped by
libs/flibs/_hamiltonian.py get_imassk / get_mass).

Both routines must reproduce the SECOND derivative of the band energy,

    (1/m*)_ab ~ d2 eps_n / dk_a dk_b ,

which has two pieces the old implementation did not have:

  * the interband (level repulsion) term of 2nd-order perturbation theory,
        2 sum_{m!=n} Re[H~^a_nm H~^b_mn] / (eps_n - eps_m) ,
    present for ANY basis and typically the same size as the diagonal term --
    the old get_imassk was diag(U^H d2H/dk2 U) alone;
  * for a NON-ORTHOGONAL (MLO) basis (ftype 3/4), the overlap terms coming from
    H C = eps S C with the k-dependent normalization C^H S C = 1.

The reference here is a central finite difference of eps_n(k) itself, which knows
about all of it.  Bands that are near-degenerate somewhere on the test k-points are
excluded: the band-resolved mass genuinely does not exist at a crossing, and both
the routine (deg_thresh) and the finite difference break down there.

Runs standalone (no pytest needed):  python tests/test_effmass.py
Also works under pytest if installed:  pytest tests/test_effmass.py
Requires the Fortran library libfmod.so (cd libs && make FC=ifx SL=MKL).
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import libs.flibs as F
import _tools as T

FD_STEP = 2.0e-4          # h for the finite-difference Hessian
DEG_TOL = 2.0e-2          # bands closer than this anywhere are dropped from the comparison


# --------------------------------------------------------------------------- #
#  helpers
# --------------------------------------------------------------------------- #
def _eig_uni(klist, ham_r, S_r, rvec):
    if S_r is None or len(S_r) == 0:
        return F.get_eig(F.gen_ham(klist, ham_r, rvec))
    hamk, Sk = F.gen_ham(klist, ham_r, rvec, Ovl_r=S_r)
    return F.get_eig(hamk, Sk)


def _fd_hessian(klist, ham_r, S_r, rvec, h=FD_STEP):
    """d2 eps_n/dk_a dk_b by central differences.  get_imass0/get_vlm0 differentiate
    the phase exp(-2 pi i k.R), i.e. d/d(2 pi k_frac), so the fractional-k result is
    divided by (2 pi)^2 to match."""
    Nk, Norb = len(klist), ham_r.shape[1]

    def E(dk):
        return _eig_uni(np.ascontiguousarray(klist + dk), ham_r, S_r, rvec)[0]

    hess = np.zeros((Nk, Norb, 3, 3))
    e0 = E(np.zeros(3))
    for a in range(3):
        da = np.zeros(3)
        da[a] = h
        hess[:, :, a, a] = (E(da) + E(-da) - 2.0 * e0) / h ** 2
        for b in range(a + 1, 3):
            db = np.zeros(3)
            db[b] = h
            mix = (E(da + db) - E(da - db) - E(-da + db) + E(-da - db)) / (4.0 * h * h)
            hess[:, :, a, b] = mix
            hess[:, :, b, a] = mix
    return hess / (2.0 * np.pi) ** 2


def _nondegenerate(eig, tol=DEG_TOL):
    """Mask [Nk,Norb] of bands well separated from every other band at that k."""
    gap = np.abs(eig[:, :, None] - eig[:, None, :]) + np.eye(eig.shape[1]) * 1e9
    return gap.min(axis=2) > tol


def _klist(Nk=6, seed=3):
    """Generic k-points, kept off the high-symmetry planes where bands stick together."""
    rng = np.random.default_rng(seed)
    return np.ascontiguousarray(rng.random((Nk, 3)) * 0.3 + 0.11)


def _imass_without_interband(klist, ham_r, rvec, mrot, uni, eig):
    """What get_imassk computed before: diag(U^H d2H/dk2 U) only.  Reproduced by
    handing the routine energies so far apart that every 1/(eps_n-eps_m) vanishes."""
    return F.get_imassk(F.get_imass0(klist, ham_r, rvec), mrot, uni,
                        F.get_vlm0(klist, ham_r, rvec),
                        np.ascontiguousarray(np.zeros_like(eig) + 1.0e9))


# --------------------------------------------------------------------------- #
#  orthogonal basis (get_imassk)
# --------------------------------------------------------------------------- #
def test_imass_matches_finite_difference():
    """get_imassk (with the interband term) reproduces d2 eps/dk2."""
    rvec, ham_r, _, _, _ = T.mlo_model(seed=11)          # H(R) part only; S is ignored below
    klist = _klist()
    eig, uni = _eig_uni(klist, ham_r, [], rvec)
    ok = _nondegenerate(eig)

    imass = F.get_mass(klist, ham_r, rvec, np.eye(3), uni, True, None, eig)
    fd = _fd_hessian(klist, ham_r, [], rvec)
    assert ok.sum() > 0
    assert np.abs(imass - fd)[ok].max() < 1e-3 * max(1.0, np.abs(fd)[ok].max())


def test_interband_term_is_not_negligible():
    """Guard against a regression to diag(U^H d2H/dk2 U): the interband term is the
    same order as the diagonal one, so dropping it misses the finite difference badly."""
    rvec, ham_r, _, _, _ = T.mlo_model(seed=11)
    klist = _klist()
    eig, uni = _eig_uni(klist, ham_r, [], rvec)
    ok = _nondegenerate(eig)

    fd = _fd_hessian(klist, ham_r, [], rvec)
    old = _imass_without_interband(klist, ham_r, rvec, np.eye(3), uni, eig)
    assert np.abs(old - fd)[ok].max() > 0.1 * np.abs(fd)[ok].max()


def test_imass_is_axis_symmetric():
    """d2/dk_a dk_b is symmetric, so the tensor must be too (both bases)."""
    for seed, use_S in ((12, False), (12, True)):
        rvec, ham_r, S_r, _, _ = T.mlo_model(seed=seed)
        S = S_r if use_S else None
        klist = _klist()
        eig, uni = _eig_uni(klist, ham_r, S if S is not None else [], rvec)
        imass = F.get_mass(klist, ham_r, rvec, np.eye(3), uni, True, S, eig)
        assert np.abs(imass - imass.transpose(0, 1, 3, 2)).max() < 1e-10


def test_axis_rotation_convention():
    """mrot acts on both axes: M_rot = mrot^T M mrot (numpy side), matching the
    v -> v @ mrot convention of get_veloc."""
    rvec, ham_r, _, _, _ = T.mlo_model(seed=13)
    klist = _klist()
    eig, uni = _eig_uni(klist, ham_r, [], rvec)
    mrot = np.ascontiguousarray([[1.0, 0.2, 0.0],
                                 [0.0, 1.1, 0.0],
                                 [0.0, 0.0, 0.9]])
    m1 = F.get_mass(klist, ham_r, rvec, np.eye(3), uni, True, None, eig)
    m2 = F.get_mass(klist, ham_r, rvec, mrot, uni, True, None, eig)
    assert np.abs(m2 - mrot.T @ m1 @ mrot).max() < 1e-10


def test_degenerate_bands_stay_finite():
    """At an exact degeneracy the band-resolved mass does not exist; the routine must
    drop those pairs (deg_thresh) instead of returning inf/NaN."""
    # two decoupled identical orbitals -> every band exactly two-fold degenerate
    rvec = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0],
                     [0, 1, 0], [0, -1, 0]], dtype=np.float64)
    ham_r = np.zeros((len(rvec), 2, 2), dtype=np.complex128)
    for j in range(1, len(rvec)):
        ham_r[j] = -0.5 * np.eye(2)
    klist = _klist()
    eig, uni = _eig_uni(klist, ham_r, [], rvec)
    assert np.abs(eig[:, 0] - eig[:, 1]).max() < 1e-12          # really degenerate
    imass = F.get_mass(klist, ham_r, rvec, np.eye(3), uni, True, None, eig)
    assert np.isfinite(imass).all()
    # a degenerate pair contributes nothing, so this reduces to the diagonal term,
    # which for decoupled bands is still the exact answer
    fd = _fd_hessian(klist, ham_r, [], rvec)
    assert np.abs(imass - fd).max() < 1e-4


# --------------------------------------------------------------------------- #
#  non-orthogonal (MLO) basis (get_imassk_mlo)
# --------------------------------------------------------------------------- #
def test_mlo_imass_matches_finite_difference():
    """The full MLO expression (overlap + interband) reproduces d2 eps/dk2."""
    rvec, ham_r, S_r, _, _ = T.mlo_model(seed=14)
    klist = _klist()
    eig, uni = _eig_uni(klist, ham_r, S_r, rvec)
    ok = _nondegenerate(eig)

    imass = F.get_mass(klist, ham_r, rvec, np.eye(3), uni, True, S_r, eig)
    fd = _fd_hessian(klist, ham_r, S_r, rvec)
    assert ok.sum() > 0
    assert np.abs(imass - fd)[ok].max() < 1e-3 * max(1.0, np.abs(fd)[ok].max())


def test_mlo_overlap_terms_are_not_negligible():
    """Using the orthogonal routine on an MLO model (i.e. ignoring S) is badly wrong."""
    rvec, ham_r, S_r, _, _ = T.mlo_model(seed=14)
    klist = _klist()
    eig, uni = _eig_uni(klist, ham_r, S_r, rvec)
    ok = _nondegenerate(eig)

    fd = _fd_hessian(klist, ham_r, S_r, rvec)
    no_ovl = F.get_mass(klist, ham_r, rvec, np.eye(3), uni, True, None, eig)
    assert np.abs(no_ovl - fd)[ok].max() > 0.1 * np.abs(fd)[ok].max()


def test_mlo_reduces_to_orthogonal_when_S_is_identity():
    """With S(R) = delta_{R,0} * 1 the MLO routine must return exactly the orthogonal one."""
    rvec, ham_r, _, Norb, Nr = T.mlo_model(seed=15)
    S_r = np.zeros((Nr, Norb, Norb), dtype=np.complex128)
    S_r[np.argmin(np.abs(rvec).sum(axis=1))] = np.eye(Norb)     # S(R=0)=1, else 0
    klist = _klist()
    eig, uni = _eig_uni(klist, ham_r, S_r, rvec)

    m_mlo = F.get_mass(klist, ham_r, rvec, np.eye(3), uni, True, S_r, eig)
    m_ort = F.get_mass(klist, ham_r, rvec, np.eye(3), uni, True, None, eig)
    assert np.abs(m_mlo - m_ort).max() < 1e-9


def test_mass_is_inverse_of_imass():
    """sw_imass=False returns the matrix inverse of the sw_imass=True tensor."""
    rvec, ham_r, S_r, _, _ = T.mlo_model(seed=16)
    klist = _klist(4)
    eig, uni = _eig_uni(klist, ham_r, S_r, rvec)
    imass = F.get_mass(klist, ham_r, rvec, np.eye(3), uni, True, S_r, eig)
    mass = F.get_mass(klist, ham_r, rvec, np.eye(3), uni, False, S_r, eig)
    prod = np.einsum('knab,knbc->knac', mass, imass)
    assert np.abs(prod - np.eye(3)).max() < 1e-8


def test_eig_required():
    """get_mass without eig is a programming error, not a silent wrong answer."""
    rvec, ham_r, S_r, _, _ = T.mlo_model(seed=17)
    klist = _klist(3)
    _, uni = _eig_uni(klist, ham_r, S_r, rvec)
    T.assert_raises(ValueError, F.get_mass, klist, ham_r, rvec, np.eye(3), uni, True)


if __name__ == '__main__':
    import _tools
    sys.exit(_tools.run_standalone(globals()))
