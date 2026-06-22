#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Regression tests for the quasiclassical Eilenberger / Riccati suite (libs/plibs/
_eilenberger*.py and the Fortran kernels in libs/flibs/_eilenberger.py).

These lock in the physics benchmarks that were validated by hand while the suite
was built: Matsubara cutoff scaling, the Anderson theorem vs Abrikosov-Gor'kov
pair breaking, the Fortran<->Python equivalence of both Riccati kernels, the
surface gap suppression and zero-energy bound state, the vortex-core CdGM peak,
the Volovik field dependence, the model-FS reproduction of the cylinder, Pauli
limiting vs triplet immunity, and the surface / vortex-core d-vector textures.

Runs standalone (no pytest needed):  python tests/test_eilenberger.py
Also works under pytest if installed:  pytest tests/test_eilenberger.py
Requires the Fortran library libfmod.so (cd libs && make FC=ifx SL=MKL).
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from libs.plibs import _eilenberger as E
from libs.plibs import _eilenberger_surface as S
from libs.plibs import _eilenberger_vortex as V
from libs.plibs import _eilenberger_spin as SP
import libs.flibs as F


# --------------------------------------------------------------------------- #
#  helpers
# --------------------------------------------------------------------------- #
def _cyl(Nb=120, gap_sym='d'):
    """Synthetic cylindrical FS: uniform weights and the normalized form factor."""
    beta = np.linspace(0.0, 2.0 * np.pi, Nb, endpoint=False)
    wf = np.ones(Nb)
    if gap_sym == 's':
        phi = np.ones(Nb)
    elif gap_sym == 'd':
        phi = np.sqrt(2.0) * np.cos(2.0 * beta)
    elif gap_sym == 'dxy':
        phi = np.sqrt(2.0) * np.sin(2.0 * beta)
    else:
        raise ValueError(gap_sym)
    return wf, phi


def _ref_chord(om, dd, hvf, ds):
    """Pure-Python reference for one scalar chord (the deleted _integrate_vec/_chord_gf
    algorithm), to guard the Fortran riccati_chords.  om, dd are [Ns, Nw]."""
    Ns, Nw = dd.shape
    t = ds / hvf

    def integ(omp, ddp, g0):
        out = np.empty((Ns, Nw), dtype=np.complex128)
        g = g0.copy()
        out[0] = g
        for i in range(Ns - 1):
            D = 0.5 * (ddp[i] + ddp[i + 1])
            o = 0.5 * (omp[i] + omp[i + 1])
            R = np.sqrt(o ** 2 + (D * np.conj(D)).real)
            T = np.tanh(R * t) / R
            g = (g + T * (D - o * g)) / (1.0 + T * (o + np.conj(D) * g))
            out[i + 1] = g
        return out

    R0 = np.sqrt(om[0] ** 2 + dd[0] * np.conj(dd[0]))
    ga0 = np.where(np.abs(dd[0]) > 0, (R0 - om[0]) / np.where(np.abs(dd[0]) > 0, np.conj(dd[0]), 1.0), 0.0)
    gamma = integ(om, dd, ga0)
    Rn = np.sqrt(om[-1] ** 2 + dd[-1] * np.conj(dd[-1]))
    gb0 = np.where(np.abs(dd[-1]) > 0, (Rn - om[-1]) / np.where(np.abs(dd[-1]) > 0, dd[-1], 1.0), 0.0)
    gammat = integ(om[::-1], np.conj(dd[::-1]), gb0)[::-1]
    den = 1.0 + gamma * gammat
    return (1.0 - gamma * gammat) / den, 2.0 * gamma / den   # g, f


def _ref_matrix_traj(omega, Dpath, hvf, ds, h):
    """Pure-Python reference for the 2x2 spin-matrix Mobius Riccati along one path
    (the deleted matrix_trajectory_gf), to guard the Fortran matrix kernel."""
    from scipy.linalg import expm
    sz = np.array([[1, 0], [0, -1]], dtype=np.complex128)
    I2 = np.eye(2, dtype=np.complex128)
    Ns = len(Dpath)
    t = ds / hvf

    def gen(D, sgn):
        N = np.empty((4, 4), dtype=np.complex128)
        N[:2, :2] = -2.0 * omega * I2 + 1j * sgn * h * sz
        N[:2, 2:] = D
        N[2:, :2] = np.conj(D).T
        N[2:, 2:] = 1j * sgn * h * sz
        return N

    def root(D, is_a):
        d2 = 0.5 * np.trace(D @ np.conj(D).T).real
        E = np.sqrt(omega ** 2 + d2)
        return (D if is_a else np.conj(D).T) / (omega + E)

    def step(a, D, sgn):
        M = expm(gen(D, sgn) * t)
        u = M[:2, :2] @ a + M[:2, 2:]
        w = M[2:, :2] @ a + M[2:, 2:]
        return u @ np.linalg.inv(w)

    a = np.empty((Ns, 2, 2), dtype=np.complex128)
    b = np.empty((Ns, 2, 2), dtype=np.complex128)
    a[0] = root(Dpath[0], True)
    for i in range(Ns - 1):
        a[i + 1] = step(a[i], 0.5 * (Dpath[i] + Dpath[i + 1]), 1.0)
    b[Ns - 1] = root(Dpath[-1], False)
    for i in range(Ns - 1, 0, -1):
        b[i - 1] = step(b[i], 0.5 * (Dpath[i] + Dpath[i - 1]), -1.0)
    g = np.empty((Ns, 2, 2), dtype=np.complex128)
    f = np.empty((Ns, 2, 2), dtype=np.complex128)
    for i in range(Ns):
        P = np.linalg.inv(I2 + a[i] @ b[i])
        g[i] = P @ (I2 - a[i] @ b[i])
        f[i] = P @ (2.0 * a[i])
    return g, f


# --------------------------------------------------------------------------- #
#  homogeneous
# --------------------------------------------------------------------------- #
def test_matsubara_cutoff_scaling():
    """Number of Matsubara frequencies below a fixed cutoff grows as ~1/T."""
    wc = 0.5
    for T in (1e-3, 5e-4, 2e-4):
        n = len(E.matsubara(T, wc))
        expect = int(np.floor((wc / (np.pi * T) - 1.0) / 2.0)) + 1
        assert n == expect
    assert len(E.matsubara(2e-4, wc)) > len(E.matsubara(1e-3, wc))


def test_anderson_theorem_vs_AG():
    """Non-magnetic impurities leave the s-wave gap (Anderson) almost unchanged but
    suppress the sign-changing d-wave gap (Abrikosov-Gor'kov pair breaking)."""
    wc, T = 0.5, 5e-4
    om = E.matsubara(T, wc)
    g_imp = 5.0e-4
    wf_s, phi_s = _cyl(gap_sym='s')
    wf_d, phi_d = _cyl(gap_sym='d')
    Ds0 = E.solve_gap(T, wf_s, phi_s, om, 0.0, 1e8, 0.5)
    Ds1 = E.solve_gap(T, wf_s, phi_s, om, g_imp, 1e8, 0.5)
    Dd0 = E.solve_gap(T, wf_d, phi_d, om, 0.0, 1e8, 0.5)
    Dd1 = E.solve_gap(T, wf_d, phi_d, om, g_imp, 1e8, 0.5)
    assert Ds0 > 0 and Dd0 > 0
    assert Ds1 / Ds0 > 0.97                  # s-wave protected
    assert Dd1 / Dd0 < 0.9                    # d-wave pair-broken
    assert Dd1 / Dd0 < Ds1 / Ds0


def test_tc_bcs_ratio():
    """Weak-coupling BCS gap ratio 2 Delta0 / kTc ~ 3.53 (Delta0/Tc ~ 1.76)."""
    wc = 0.3
    wf, phi = _cyl(gap_sym='s')
    Tc = E.find_tc(wf, phi, 0.0, 1e8, 0.4, wc, 1e-5, 0.05)
    assert Tc > 0
    D0 = E.solve_gap(0.02 * Tc, wf, phi, E.matsubara(0.02 * Tc, wc), 0.0, 1e8, 0.4)
    assert 1.55 < D0 / Tc < 1.95             # ~1.76


# --------------------------------------------------------------------------- #
#  Fortran <-> reference kernels
# --------------------------------------------------------------------------- #
def test_riccati_chords_matches_reference():
    """Fortran scalar tanh-step kernel == the pure-Python reference (to ~1e-13)."""
    rng = np.random.default_rng(1)
    Ns, Nw = 40, 8
    om = ((2 * np.arange(Nw) + 1) * np.pi * 1e-3)[None, :] * np.ones((Ns, 1))
    om = om + 1j * 3e-4 * rng.standard_normal((Ns, 1))         # position-dependent
    dd = (rng.standard_normal((Ns, Nw)) + 1j * rng.standard_normal((Ns, Nw))) * 1e-3
    hvf, ds = 1.3, 0.7
    gR, fR = _ref_chord(om, dd, hvf, ds)
    gF, fF = F.riccati_chords(om[:, None, :], dd[:, None, :], hvf, ds)
    assert np.abs(gF[:, 0] - gR).max() < 1e-12
    assert np.abs(fF[:, 0] - fR).max() < 1e-12


def test_matrix_riccati_matches_python():
    """Fortran 2x2 matrix Mobius kernel == the pure-Python reference (to ~1e-11),
    including a finite Zeeman field."""
    rng = np.random.default_rng(2)
    ch = SP._default_dvector_channels()
    Smats = [c[2] for c in ch]
    Ns = 30
    Dpath = np.zeros((Ns, 2, 2), dtype=np.complex128)
    for a in range(2):
        Dpath += (rng.standard_normal(Ns) * 1e-3)[:, None, None] * Smats[a]
    om = (2 * np.arange(5) + 1) * np.pi * 1.5e-3
    for h in (0.0, 2e-3):
        gF, fF = F.matrix_riccati_batch(om.astype(complex), Dpath, 1.0, 0.6, h)
        for iw, wn in enumerate(om):
            gp, fp = _ref_matrix_traj(complex(wn), Dpath, 1.0, 0.6, h)
            assert np.abs(gF[:, iw] - gp).max() < 1e-11
            assert np.abs(fF[:, iw] - fp).max() < 1e-11
    # the batched-chords kernel must agree with the single-trajectory kernel
    gB, fB = F.matrix_riccati_batch(om.astype(complex), Dpath, 1.0, 0.6, 0.0)
    gC, fC = F.matrix_riccati_chords(om.astype(complex), Dpath[:, None], 1.0, 0.6, 0.0)
    assert np.abs(gC[:, 0] - gB).max() < 1e-13


# --------------------------------------------------------------------------- #
#  surface
# --------------------------------------------------------------------------- #
def test_surface_sign_changing_suppression():
    """The sign-changing d[110] gap is suppressed at a specular surface and heals to
    the bulk; the non-sign-changing d[100] stays close to flat."""
    wc, T = 0.5, 5e-4
    om = E.matsubara(T, wc)
    x, D110, Db = S.solve_surface(0.5, T, om, 'd', np.pi / 4, Nbeta=24, Lxi=8, nper=10)
    x, D100, Db2 = S.solve_surface(0.5, T, om, 'd', 0.0, Nbeta=24, Lxi=8, nper=10)
    assert D110[0] / Db < 0.1                 # strongly suppressed at the surface
    assert D110[-1] / Db > 0.95               # heals to bulk
    assert D100[0] / Db2 > 0.8                 # [100] surface barely suppressed


def test_surface_zero_energy_bound_state():
    """d[110] has a zero-energy surface bound state (large N(0,0)); d[100] does not."""
    wc, T = 0.5, 5e-4
    om = E.matsubara(T, wc)
    x, D110, Db = S.solve_surface(0.5, T, om, 'd', np.pi / 4, Nbeta=24, Lxi=8, nper=10)
    w = np.linspace(-2 * Db, 2 * Db, 41)
    n110 = S.surface_ldos(D110, x, w, 'd', np.pi / 4, ix=0, Dbulk=Db, Nbeta=48)
    x, D100, Db2 = S.solve_surface(0.5, T, om, 'd', 0.0, Nbeta=24, Lxi=8, nper=10)
    n100 = S.surface_ldos(D100, x, w, 'd', 0.0, ix=0, Dbulk=Db2, Nbeta=48)
    i0 = np.argmin(np.abs(w))
    assert n110[i0] > 5.0                      # zero-bias peak
    assert n110[i0] > 3.0 * n100[i0]           # absent for [100]


# --------------------------------------------------------------------------- #
#  vortex / lattice
# --------------------------------------------------------------------------- #
def test_vortex_core_peak():
    """Order parameter vanishes at the vortex core and a CdGM zero-energy peak sits
    there (N(core,0) well above the bulk N=1)."""
    wc, T = 0.5, 8e-4
    om = E.matsubara(T, wc)
    xg, Psi, Db, xi = V.solve_vortex2d(0.6, T, om, 'd', Lxi=7, ngrid=31, nbeta=18)
    ic = len(xg) // 2
    assert abs(Psi[ic, ic]) / Db < 0.05       # |Psi|->0 at the core
    n0 = V.vortex_ldos2d(Psi, xg, xi, np.array([0.0]), 'd', Db, nbeta=36)[ic, ic, 0]
    assert n0 > 5.0                            # core bound-state peak


def test_lattice_volovik_field_dependence():
    """The spatially-averaged zero-energy DOS of a d-wave vortex lattice grows with
    field (Volovik)."""
    wc, T = 0.5, 8e-4
    om = E.matsubara(T, wc)
    n = []
    for b in (0.1, 0.3):
        st = V.solve_lattice(0.6, T, om, gap_sym='d', field=b, kappa=5.0, Ng=14, nbeta=12)
        n.append(float(V.lattice_dos(st, 'd', np.array([0.0]), nbeta=48, delta=0.02 * st['Dbulk'])[0]))
    assert 0.0 < n[0] < n[1]                   # rises with B


# --------------------------------------------------------------------------- #
#  condensation free energy / supercurrent (je observables)
# --------------------------------------------------------------------------- #
def test_condensation_energy_universal():
    """The homogeneous condensation free energy is negative and coupling-independent:
    dOmega(0)/Damp0^2 -> a universal constant (~ -1/4 per spin)."""
    vals = []
    for g in (0.45, 0.5, 0.55):
        dO, D = E.condensation_energy(g, 1e-5, 0.3, 's')
        assert D > 0 and dO < 0
        vals.append(dO / D ** 2)
    assert max(vals) - min(vals) < 5e-3            # coupling-independent
    assert abs(np.mean(vals) - (-0.25)) < 0.02     # per-spin BCS value


def test_condensation_energy_vanishes_at_tc():
    """dOmega -> 0 as T -> Tc and grows in magnitude as T decreases."""
    g, wc = 0.4, 0.3
    Tc = E.find_tc(np.ones(240), np.ones(240), 0.0, 1e8, g, wc, 1e-5, 0.05)
    dO_lo, _ = E.condensation_energy(g, 0.2 * Tc, wc, 's')
    dO_hi, _ = E.condensation_energy(g, 0.95 * Tc, wc, 's')
    assert dO_lo < dO_hi < 0                        # |dOmega| larger at low T, ->0 at Tc
    assert abs(dO_hi) < 0.1 * abs(dO_lo)


def test_vortex_supercurrent_circulates():
    """The vortex supercurrent vanishes at the core, peaks near a coherence length,
    circulates with a single sign, and decays outward."""
    wc, T = 0.5, 8e-4
    om = E.matsubara(T, wc)
    xg, Psi, Db, xi = V.solve_vortex2d(0.6, T, om, 'd', Lxi=7, ngrid=31, nbeta=18)
    jx, jy = V.vortex_current2d(Psi, xg, xi, om, T, 'd', nbeta=18)
    rho, jphi = V.vortex_current_profile(jx, jy, xg)
    ipk = np.argmax(np.abs(jphi))
    assert abs(jphi[0]) < 0.05 * abs(jphi[ipk])    # ~0 at the core
    assert 0.5 < rho[ipk] / xi < 2.5               # peak near xi
    s = np.sign(jphi[ipk])
    assert np.all(jphi[1:] * s > -0.02 * abs(jphi[ipk]))   # single-sign circulation
    assert abs(jphi[-1]) < abs(jphi[ipk])          # decays outward


# --------------------------------------------------------------------------- #
#  model Fermi surface
# --------------------------------------------------------------------------- #
def test_model_fs_basics_and_cylinder_limit():
    """build_model_fs is normalized (sum nf=1), the isotropic FS has |v_F|=const, and
    its bulk gap matches the cylinder."""
    wc, T = 0.5, 8e-4
    om = E.matsubara(T, wc)
    fi = E.build_model_fs('iso', 120)
    assert abs(fi['nf'].sum() - 1.0) < 1e-12
    assert fi['vabs'].std() / fi['vabs'].mean() < 1e-6     # |v_F| isotropic
    D_fs = E.bulk_gap_fs(0.6, T, om, fi, 's')
    wf, phi = _cyl(gap_sym='s')
    D_cyl = E.solve_gap(T, wf, phi, om, 0.0, 1e8, 0.6)
    assert abs(D_fs - D_cyl) / D_cyl < 0.02


# --------------------------------------------------------------------------- #
#  spin: Pauli limiting vs triplet immunity
# --------------------------------------------------------------------------- #
def test_pauli_singlet_vs_triplet_immunity():
    """A Zeeman field suppresses the singlet gap (Pauli) but not a triplet with the
    d-vector perpendicular to the field (equal-spin pairing)."""
    wc, T = 0.5, 6e-4
    om = E.matsubara(T, wc)
    wf, phi = _cyl(Nb=64, gap_sym='s')
    D0 = SP.solve_gap_spin(0.6, T, om, wf, phi, 'singlet', h=0.0, damp_init=2e-3)
    h = 0.7 * D0                               # below the Chandrasekhar-Clogston jump
    Ds = SP.solve_gap_spin(0.6, T, om, wf, phi, 'singlet', h=h, damp_init=D0)
    Dt = SP.solve_gap_spin(0.6, T, om, wf, phi, 'triplet', dvec=(1, 0, 0), h=h, damp_init=D0)
    assert Ds / D0 < 0.98                      # singlet Pauli-suppressed
    assert Dt / D0 > 0.99                      # triplet d_|_h immune
    assert Dt / D0 > Ds / D0 + 0.02            # triplet clearly less suppressed


# --------------------------------------------------------------------------- #
#  d-vector textures
# --------------------------------------------------------------------------- #
def test_dvector_texture_surface():
    """At a surface the sign-changing dominant component is pair-broken and the
    subdominant rotates the d-vector (theta_d larger at the surface than in bulk)."""
    wc, T = 0.2, 1.5e-3
    om = E.matsubara(T, wc)
    x, D, Db = SP.solve_surface_dvector((0.8, 0.8 * 0.9), T, om, Lxi=7, nper=6, Nbeta=10, itemax=40)
    th = lambda j: np.degrees(np.arctan2(abs(D[1, j]), abs(D[0, j])))
    assert th(0) > th(-1) + 20.0               # d-vector rotates toward the surface
    assert th(0) > 70.0                         # nearly pure subdominant at the surface


def test_dvector_texture_vortex_core():
    """Around a vortex the winding dominant vanishes in the core where the
    core-localized subdominant survives; the d-vector tilts to ~90 deg at the core."""
    wc, T = 0.2, 1.5e-3
    om = E.matsubara(T, wc)
    xg, A, Db, xi = SP.solve_vortex2d_dvector((0.8, 0.8 * 0.95), T, om, windings=(1, 0),
                                              Lxi=7, ngrid=25, nbeta=10, itemax=30)
    ic = len(xg) // 2
    th = lambda j: np.degrees(np.arctan2(abs(A[1, ic + j, ic]), abs(A[0, ic + j, ic])))
    assert abs(A[0, ic, ic]) / np.max(np.abs(Db)) < 0.1    # dominant ->0 in the core
    assert th(0) > 75.0                                     # d ~ pure subdominant at core
    assert th(0) > th(len(xg) // 2 - 1) + 20.0             # texture: core vs edge


# --------------------------------------------------------------------------- #
#  standalone runner (no pytest required)
# --------------------------------------------------------------------------- #
if __name__ == '__main__':
    import time
    tests = [v for k, v in sorted(globals().items()) if k.startswith('test_') and callable(v)]
    npass = 0
    for t in tests:
        t0 = time.time()
        try:
            t()
            print(f"  PASS  {t.__name__:42s} ({time.time()-t0:5.1f}s)", flush=True)
            npass += 1
        except Exception as e:
            print(f"  FAIL  {t.__name__:42s} ({time.time()-t0:5.1f}s)  -> {type(e).__name__}: {e}", flush=True)
    print(f"\n{npass}/{len(tests)} passed")
    sys.exit(0 if npass == len(tests) else 1)
