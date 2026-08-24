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
import libs.plibs as P
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


def test_tilted_cylinder_orbit_area_follows_inverse_cosine():
    """A quasi-2D cylinder has A_ext(theta) = A_ext(0) / cos(theta)."""
    rvec = np.array(
        [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0],
         [0.0, 1.0, 0.0], [0.0, -1.0, 0.0]],
        dtype=np.float64,
    )
    ham_r = np.full((4, 1, 1), -1.0, dtype=np.complex128)
    mesh, mu, abz = 160, -3.0, 4.0 * np.pi**2
    areas = []
    for theta in (0.0, 35.0, 45.0):
        rotmat = P.make_rotmat(theta, 0.0)
        scan, _ = P.scan_fs_area(mesh, rvec, ham_r, [], rotmat, mu, abz, meshkz=9)
        kz_arr, area_arr = np.asarray(scan[0]).T
        candidates = P.find_extremal_kz(
            kz_arr, area_arr, 0, mesh, rvec, ham_r, [], rotmat, mu, abz,
        )
        assert candidates == ([0.0, 0.5] if theta == 0.0 else [0.0])
        eig = P.get_eigs_2d(mesh, rvec, ham_r, [], rotmat, candidates[0])
        contours, bands = P.get_kf_points(eig, mesh, mu, candidates[0])
        areas.append(P.get_band_area(contours, bands, 0, abz))
    for theta, area in zip((35.0, 45.0), areas[1:]):
        assert np.isclose(area / areas[0], 1.0 / np.cos(np.deg2rad(theta)), rtol=3e-3)


def _square_cylinder():
    """One-orbital square lattice with in-plane hoppings only: a cylinder along c."""
    rvec = np.array(
        [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0],
         [0.0, 1.0, 0.0], [0.0, -1.0, 0.0]],
        dtype=np.float64,
    )
    return rvec, np.full((4, 1, 1), -1.0, dtype=np.complex128)


_HEX = np.array([[3.0, 0.0, 0.0], [-1.5, 1.5 * np.sqrt(3.0), 0.0], [0.0, 0.0, 5.0]])


def test_slice_area_factor_is_the_reciprocal_cell_area():
    """At rotmat=I the Jacobian must be |b1 x b2|, not 4*pi^2/(|a1||a2|)."""
    for avec in (np.diag([3.0, 4.0, 5.0]), _HEX):
        bvec = P.get_bvec(avec)
        exact = np.linalg.norm(np.cross(bvec[0], bvec[1]))
        assert np.isclose(P.slice_area_factor(np.eye(3), bvec), exact, rtol=1e-12)
    # orthorhombic is the case the old hard-coded formula got right
    bvec = P.get_bvec(np.diag([3.0, 4.0, 5.0]))
    assert np.isclose(P.slice_area_factor(np.eye(3), bvec), 4.0 * np.pi**2 / 12.0)
    # a hexagonal cell is NOT: 4*pi^2/(a*a) misses the 1/sin(gamma)
    bvec = P.get_bvec(_HEX)
    assert not np.isclose(P.slice_area_factor(np.eye(3), bvec), 4.0 * np.pi**2 / 9.0)
    # a cubic cell has an isotropic Jacobian, so tilting cannot change it
    bvec = P.get_bvec(np.eye(3))
    facs = [P.slice_area_factor(P.make_rotmat_bfield(t, 20.0, bvec), bvec)
            for t in (0.0, 30.0, 60.0)]
    assert np.allclose(facs, 4.0 * np.pi**2)


def test_slice_is_perpendicular_to_the_requested_field():
    """make_rotmat_bfield must realise the requested (theta,phi) for any cell shape."""
    for avec in (np.diag([1.0, 1.0, 3.0]), _HEX, np.eye(3)):
        bvec = P.get_bvec(avec)
        for theta, phi in ((0.0, 0.0), (25.0, 0.0), (45.0, 30.0), (70.0, 115.0)):
            th, ph = np.deg2rad(theta), np.deg2rad(phi)
            want = np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])
            got = P.field_direction(P.make_rotmat_bfield(theta, phi, bvec), avec)
            assert np.allclose(got, want, atol=1e-12)
    # the plain fractional rotation only gets it right for a cubic cell
    bvec, avec = P.get_bvec(np.diag([1.0, 1.0, 3.0])), np.diag([1.0, 1.0, 3.0])
    got = P.field_direction(P.make_rotmat(45.0, 0.0), avec)
    assert np.isclose(np.degrees(np.arccos(got[2])), np.degrees(np.arctan(1.0 / 3.0)))


def test_tetragonal_cylinder_orbit_area_follows_inverse_cosine():
    """Same 1/cos(theta) law as the cubic test, but for c/a=3 where the old
    fixed ABZ and fractional rotation both failed."""
    rvec, ham_r = _square_cylinder()
    avec = np.diag([1.0, 1.0, 3.0])
    bvec = P.get_bvec(avec)
    mesh, mu = 240, -3.0
    P.open_contour_report()
    areas = {}
    for theta in (0.0, 15.0, 30.0, 45.0):
        rotmat = P.make_rotmat_bfield(theta, 0.0, bvec)
        abz = P.slice_area_factor(rotmat, bvec)
        eig = P.get_eigs_2d(mesh, rvec, ham_r, [], rotmat, 0.0)
        areas[theta] = P.get_band_area(*P.get_kf_points(eig, mesh, mu, 0.0), 0, abz)
    for theta in (15.0, 30.0):
        assert np.isclose(areas[theta] / areas[0.0],
                          1.0 / np.cos(np.deg2rad(theta)), rtol=1e-3)
    # at 45 deg the elongated orbit no longer fits in the slice box: it must be
    # reported as unavailable rather than silently truncated
    assert areas[45.0] is None
    assert P.open_contour_report() == 1


def test_open_contour_is_rejected():
    """A pocket straddling the zone boundary comes back as open pieces whose
    shoelace sum is not the orbit area."""
    rvec, ham_r = _square_cylinder()
    mesh, mu, abz = 200, 1.0, 4.0 * np.pi**2   # hole pocket around the zone corner
    eig = P.get_eigs_2d(mesh, rvec, ham_r, [], np.eye(3), 0.0)
    contours, bands = P.get_kf_points(eig, mesh, mu, 0.0)
    assert not all(np.allclose(c[0], c[-1]) for c in contours[0])
    P.open_contour_report()
    assert P.get_band_area(contours, bands, 0, abz) is None
    assert P.open_contour_report() == 1
    truncated = P.get_band_area(contours, bands, 0, abz, sw_strict=False)
    exact = 0.30831 * abz          # fraction of the BZ with eps > mu
    assert truncated < 0.5 * exact  # ... and it is wrong by a factor of ~3.6


def test_zone_boundary_extremal_gate():
    """kz=0.5 is extremal iff rotmat[2] is a reciprocal lattice vector, and kz=0 is
    only a candidate when the scan actually reaches it."""
    rvec, ham_r = _square_cylinder()
    kz = np.linspace(0.0, 0.5, 6)
    S = 1.0 + 0.1 * kz                     # monotone: no interior extremum to refine
    for theta, want in ((0.0, [0.0, 0.5]), (90.0, [0.0, 0.5]), (45.0, [0.0])):
        got = P.find_extremal_kz(kz, S, 0, 4, rvec, ham_r, [],
                                 P.make_rotmat(theta, 0.0), 0.0, 1.0)
        assert [round(k, 6) for k in got] == want
    # a band whose FS only appears at finite kz: kz_arr[0] is an ordinary point
    kz2 = np.linspace(0.2, 0.5, 6)
    got = P.find_extremal_kz(kz2, 1.0 + 0.1 * kz2, 0, 4, rvec, ham_r, [],
                             np.eye(3), 0.0, 1.0)
    assert [round(k, 6) for k in got] == [0.5]


def _cubic_tb(tx=1.0, ty=1.0, tz=1.0):
    """One-orbital tight binding on a cubic cell: eps = -2 sum_i t_i cos(k_i a)."""
    rvec = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0],
                     [0, -1, 0], [0, 0, 1], [0, 0, -1]], dtype=np.float64)
    ham_r = np.array([-tx, -tx, -ty, -ty, -tz, -tz],
                     dtype=np.complex128).reshape(6, 1, 1)
    return rvec, ham_r


def _extremal_orbits(rvec, ham_r, avec, mu, theta, phi, meshkz=25, dk=2 * np.pi / 240,
                     egrid=None):
    """The full plane-slicing pipeline, exactly as main.get_dhva_band drives it."""
    bvec = P.get_bvec(avec)
    bhat = P.bfield_hat(theta, phi)
    kf = P.fs_cartesian_points(rvec, ham_r, [], avec, mu, mesh=40, egrid=egrid)
    radius, dmax = P.fs_extent_along(kf, bhat, avec)
    period = P.slice_period(bvec, bhat)
    d_end = min(dmax, 0.5 * period) if period is not None else dmax
    sym = d_end if (period is not None and 0.5 * period <= dmax * (1 + 1e-9)) else None
    d_arr = np.linspace(0.0, d_end, meshkz, True)
    _, cache, _, _ = P.scan_plane_area(rvec, ham_r, [], avec, mu, bhat, d_arr,
                                       2.6 * radius, dk, egrid=egrid)
    ext = [e for br in P.track_orbit_branches(cache, d_arr)
           for e in P.branch_extrema(br, d_end=sym)]
    return P.dedup_extremal_orbits(ext, bvec)


def test_slice_period_commensurate_and_incommensurate():
    """The scan range must come from a real commensuration, never from a near-miss."""
    cub = P.get_bvec(np.eye(3))
    assert np.isclose(P.slice_period(cub, P.bfield_hat(0.0, 0.0)), 2 * np.pi)
    assert np.isclose(P.slice_period(cub, P.bfield_hat(90.0, 0.0)), 2 * np.pi)
    assert P.slice_period(cub, P.bfield_hat(35.0, 17.0)) is None
    # fcc: B along [111] repeats with b1.B_hat
    fcc = P.get_bvec(4.0 * np.array([[0, .5, .5], [.5, 0, .5], [.5, .5, 0]]))
    n = np.ones(3) / np.sqrt(3.0)
    assert np.isclose(P.slice_period(fcc, n), abs(fcc[0].dot(n)))
    # a field a whisker off [111] is incommensurate, and must not report a tiny period
    off = np.array([1.0, 1.0, 1.0 + 1e-3]); off /= np.linalg.norm(off)
    assert P.slice_period(fcc, off) is None


def test_plane_slice_reproduces_the_tilted_cylinder_law():
    """A(theta) = A(0)/cos(theta) at a tilt where the orbit no longer fits in one cell."""
    rvec = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0],
                     [0.0, 1.0, 0.0], [0.0, -1.0, 0.0]], dtype=np.float64)
    ham_r = np.full((4, 1, 1), -1.0, dtype=np.complex128)
    avec, mu, dk = np.diag([1.0, 1.0, 3.0]), -3.0, 2 * np.pi / 240
    a0 = None
    for theta in (0.0, 60.0, 75.0):
        bhat = P.bfield_hat(theta, 0.0)
        orbits, n_open, _, _ = P.find_plane_orbits(rvec, ham_r, [], avec, mu, bhat,
                                                   0.0, 4.0, dk, ngrow=4)
        o = min(orbits, key=lambda o: np.linalg.norm(o['cen2d']))
        if a0 is None:
            a0 = o['area']
        assert np.isclose(o['area'] / a0, 1.0 / np.cos(np.deg2rad(theta)), rtol=1e-3)


def test_plane_slice_is_isotropic_for_a_spherical_pocket():
    """A nearly free-electron pocket has an angle-independent extremal orbit, at
    field directions with no commensurate relation to the lattice."""
    rvec, ham_r = _cubic_tb()
    avec, mu = np.eye(3), -5.9
    areas = [max(e['area'] for e in _extremal_orbits(rvec, ham_r, avec, mu, th, ph))
             for th, ph in ((0.0, 0.0), (40.0, 17.0), (54.7356, 45.0), (71.3, 133.0))]
    assert max(areas) / min(areas) - 1.0 < 5e-3


def test_corrugated_cylinder_gives_belly_and_neck():
    """Both extremal orbits of a warped cylinder, against a 1D quadrature reference."""
    from scipy.integrate import quad

    def exact(mu_eff):
        xm = np.arccos(np.clip(mu_eff - 1.0, -1.0, 1.0)) / (2 * np.pi)
        y0 = lambda x: np.arccos(np.clip(mu_eff - np.cos(2 * np.pi * x), -1.0, 1.0)) / (2 * np.pi)
        return 4 * quad(y0, 0.0, xm, limit=400)[0] * 4 * np.pi**2

    tz, mu = 0.1, -3.0
    rvec, ham_r = _cubic_tb(tz=tz)
    got = sorted(e['area'] for e in
                 _extremal_orbits(rvec, ham_r, np.eye(3), mu, 0.0, 0.0,
                                  dk=2 * np.pi / 480))
    want = sorted([exact(-(mu + 2 * tz) / 2.0), exact(-(mu - 2 * tz) / 2.0)])
    assert len(got) == 2, f'expected belly and neck, got {got}'
    for g, w in zip(got, want):
        assert np.isclose(g, w, rtol=2e-4)


def test_energy_grid_interpolation_matches_direct_diagonalisation():
    """The interpolated sweep must reproduce the exact one it replaces."""
    tz, mu = 0.1, -3.0
    rvec, ham_r = _cubic_tb(tz=tz)
    direct = sorted(e['area'] for e in
                    _extremal_orbits(rvec, ham_r, np.eye(3), mu, 0.0, 0.0))
    egrid = P.build_fs_energy_grid(rvec, ham_r, [], mu, mesh=80)
    assert egrid is not None and egrid['bands'] == [0]
    interp = sorted(e['area'] for e in
                    _extremal_orbits(rvec, ham_r, np.eye(3), mu, 0.0, 0.0, egrid=egrid))
    assert len(interp) == len(direct)
    for a, b in zip(interp, direct):
        assert np.isclose(a, b, rtol=1e-4)


def test_dedup_separates_zone_copies_from_neighbouring_extrema():
    """The two questions the reduction must not confuse: 'same orbit, other zone?'
    is settled by the centre, 'same frequency, other place?' by the area."""
    bvec = P.get_bvec(np.diag([2.8, 2.8, 6.5]))
    rec = lambda band, area, cen: {'band': band, 'area': area, 'd': 0.0,
                                   'cen3d': np.asarray(cen, dtype=float)}
    base = rec(0, 1.0000, [0.3, 0.0, 0.0])
    # (a) exact lattice copy: same orbit, however far away it sits
    copy = rec(0, 1.0004, base['cen3d'] + bvec[0] - 2 * bvec[2])
    assert len(P.dedup_extremal_orbits([base, copy], bvec, tol_k=0.05)) == 1
    # (b) neighbouring extremum of the same tube, area within 0.1%: must survive.
    #     The old area-only rule (rtol 2e-3) collapsed exactly this case.
    neigh = rec(0, 1.0010, [0.3, 0.0, 0.4])
    assert len(P.dedup_extremal_orbits([base, neigh], bvec, tol_k=0.05)) == 2
    # (c) point-group image: no lattice vector relates them, but one frequency
    image = rec(0, 1.0000, [0.0, 0.3, 0.0])
    assert len(P.dedup_extremal_orbits([base, image], bvec, tol_k=0.05)) == 1
    # (d) a different band is never merged
    other = rec(1, 1.0000, [0.3, 0.0, 0.0])
    assert len(P.dedup_extremal_orbits([base, other], bvec, tol_k=0.05)) == 2
    # the reduced centre is what stage 1 compares, and it must fold copies together
    inv = np.linalg.inv(bvec)
    assert np.allclose(P.orbit_key(base, inv), P.orbit_key(copy, inv), atol=1e-9)


def test_corrugated_cylinder_extrema_survive_the_reduction():
    """Belly and neck of one tube differ by 35% here, but the reduction must keep
    both while folding away the copies the scan sees in neighbouring zones."""
    tz, mu = 0.1, -3.0
    rvec, ham_r = _cubic_tb(tz=tz)
    ext = _extremal_orbits(rvec, ham_r, np.eye(3), mu, 0.0, 0.0)
    assert len(ext) == 2
    assert all('_key' not in e for e in ext)


def test_encloses_origin():
    """The test that decides whether a section through the centre is a closed orbit."""
    ring = lambda r, c=(0.0, 0.0): np.stack(
        [c[0] + r * np.cos(np.linspace(0, 2 * np.pi, 97)),
         c[1] + r * np.sin(np.linspace(0, 2 * np.pi, 97))], axis=1)
    assert P.encloses_origin(ring(1.0))
    assert not P.encloses_origin(ring(0.4, c=(2.0, 0.0)))     # off to one side
    assert not P.encloses_origin(ring(0.4, c=(-2.0, 0.0)))    # and to the other
    assert P.encloses_origin(ring(3.0, c=(1.0, 1.0)))         # big and off-centre
    assert not P.encloses_origin(ring(1.0, c=(1.5, 0.0)))     # just misses
    # a contour whose bounding box contains the origin but which does not
    horseshoe = np.array([[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [0.5, 1.0],
                          [0.5, -0.5], [-1.0, -0.5], [-1.0, -1.0]])
    assert not P.encloses_origin(horseshoe)


def test_open_orbit_is_recognised_not_merely_discarded():
    """A section that does not close is a physical statement, not a failure: it says
    the orbit is open and carries no dHvA oscillation. Cutting a cylinder lengthwise
    is the simplest case."""
    rvec = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0],
                     [0.0, 1.0, 0.0], [0.0, -1.0, 0.0]], dtype=np.float64)
    ham_r = np.full((4, 1, 1), -1.0, dtype=np.complex128)
    avec, mu, dk = np.diag([1.0, 1.0, 3.0]), -3.0, 2 * np.pi / 240

    def probe(theta):
        orbits, _, _, open_bands = P.find_plane_orbits(
            rvec, ham_r, [], avec, mu, P.bfield_hat(theta, 0.0), 0.0, 4.0, dk, ngrow=4)
        return len(orbits), open_bands

    n, open_bands = probe(0.0)          # B along the tube: a closed cross-section
    assert n == 1 and open_bands == set()
    n, open_bands = probe(60.0)         # still closed, only wider
    assert n == 1 and open_bands == set()
    n, open_bands = probe(90.0)         # B across the tube: the section runs away
    assert n == 0 and open_bands == {0}


def test_follow_orbit_carries_an_orbit_to_the_next_angle():
    """An orbit already known must not have to be rediscovered: following it from its
    own centre must land on the same orbit at the next angle, belly and neck alike."""
    tz, mu = 0.1, -3.0
    rvec, ham_r = _cubic_tb(tz=tz)
    avec, bvec, dk = np.eye(3), P.get_bvec(np.eye(3)), 2 * np.pi / 240
    egrid = P.build_fs_energy_grid(rvec, ham_r, [], mu, mesh=120)
    seeds = _extremal_orbits(rvec, ham_r, avec, mu, 0.0, 0.0, egrid=egrid)
    assert len(seeds) == 2                       # belly and neck of the warped tube
    bhat = P.bfield_hat(10.0, 0.0)
    kf = P.fs_cartesian_points(rvec, ham_r, [], avec, mu, mesh=40, egrid=egrid)
    radius, dmax = P.fs_extent_along(kf, bhat, avec)
    step = min(dmax, 0.5 * P.slice_period(bvec, P.bfield_hat(0.0, 0.0))) / 24
    got = []
    for sd in seeds:
        f = P.follow_orbit(sd, rvec, ham_r, [], avec, mu, bhat, 2.6 * radius, dk,
                           step, egrid=egrid)
        assert f is not None, f"lost the orbit at {sd['area']:.4f}"
        assert f['band'] == sd['band']
        got.append(f['area'])
    # a tube tilted by theta has A -> A/cos(theta) to leading order
    for a, s in zip(sorted(got), sorted(e['area'] for e in seeds)):
        assert np.isclose(a, s / np.cos(np.deg2rad(10.0)), rtol=5e-3)


def test_follow_orbit_returns_none_when_the_orbit_is_gone():
    """Following must fail cleanly rather than latch onto whatever else is nearby.
    Turning B across the tube leaves the section open, so the belly has no successor."""
    tz, mu = 0.1, -3.0
    rvec, ham_r = _cubic_tb(tz=tz)
    avec, dk = np.eye(3), 2 * np.pi / 240
    seeds = _extremal_orbits(rvec, ham_r, avec, mu, 0.0, 0.0)
    belly = max(seeds, key=lambda e: e['area'])
    assert P.follow_orbit(belly, rvec, ham_r, [], avec, mu, P.bfield_hat(90.0, 0.0),
                          8.0, dk, 0.1) is None
    # a periodic image of the same orbit is the same orbit, and must still be found
    image = dict(belly)
    image['cen3d'] = belly['cen3d'] + np.array([0.0, 0.0, 2 * np.pi])
    f = P.follow_orbit(image, rvec, ham_r, [], avec, mu, P.bfield_hat(0.0, 0.0),
                       4.0, dk, 0.1)
    assert f is not None and np.isclose(f['area'], belly['area'], rtol=1e-3)


if __name__ == '__main__':
    import _tools
    sys.exit(_tools.run_standalone(globals()))
