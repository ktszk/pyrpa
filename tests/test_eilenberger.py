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


def test_multivortex_supercell_reduces():
    """A regular n^2-flux supercell of the periodic lattice reproduces the single-
    vortex cell DOS (the structure factor selects the primitive reciprocal lattice)."""
    wc, T = 0.5, 8e-4
    om = E.matsubara(T, wc)
    n0 = []
    for nfl, Ng in ((1, 14), (4, 28)):
        st = V.solve_lattice(0.6, T, om, gap_sym='d', field=0.2, kappa=5.0, Ng=Ng, nbeta=12, nflux=nfl)
        assert abs(st['S'] / nfl - V.solve_lattice(0.6, T, om, gap_sym='d', field=0.2,
                   kappa=5.0, Ng=14, nbeta=12, nflux=1)['S']) < 1e-6 * st['S']   # area/vortex fixed
        n0.append(float(V.lattice_dos(st, 'd', np.array([0.0]), nbeta=48, delta=0.02 * st['Dbulk'])[0]))
    assert abs(n0[0] - n0[1]) < 0.03 * n0[0]            # same DOS as the primitive cell


def test_field_tilt_enhances_zeeman():
    """A field tilt away from the c-axis raises the effective Maki energy (h/cos theta),
    splitting the s-wave vortex-core zero-energy peak: N(core,0) decreases with tilt."""
    wc, T = 0.5, 8e-4
    om = E.matsubara(T, wc)
    xg, Psi, Db, xi = V.solve_vortex2d(0.6, T, om, 's', Lxi=7, ngrid=29, nbeta=14, field=0.1)
    ic = len(xg) // 2
    h0 = 1.0e-3
    n_lo = V.vortex_ldos2d(Psi, xg, xi, np.array([0.0]), 's', Db, nbeta=24, field=0.1, h=h0)[ic, ic, 0]
    n_hi = V.vortex_ldos2d(Psi, xg, xi, np.array([0.0]), 's', Db, nbeta=24, field=0.1,
                           h=h0 / np.cos(np.radians(60.0)))[ic, ic, 0]
    assert n_hi < 0.8 * n_lo                            # tilt (larger h_eff) splits the core peak


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


def test_lattice_sc_formulation_a_node_and_volovik():
    """je-style self-consistent periodic lattice (formulation A): the complex order
    parameter develops a TRUE node at every core (|Psi|->0 there), the amplitude stays
    bounded (no binning blow-up), and the zero-energy DOS grows with field, sub-linearly
    (d-wave Volovik: between sqrt(B) and linear over this moderate field range)."""
    wc, T = 0.5, 8e-4
    om = E.matsubara(T, wc)
    n0 = {}
    for b in (0.1, 0.2):
        st = V.solve_lattice_sc(0.6, T, om, gap_sym='d', field=b, Ng=14, nbeta=12,
                                itemax=60, mix=0.4, eps=3e-3)
        D, Db = st['absD'], st['Dbulk']
        assert D.min() / Db < 0.3                   # core node (sampled at the nearest grid point)
        assert D.max() / Db < 1.05                  # bounded amplitude (no scatter blow-up)
        n0[b] = float(V.lattice_dos_sc(st, 'd', np.array([0.0]), nbeta=36, delta=0.03 * Db)[0])
    assert 0.0 < n0[0.1] < n0[0.2]                  # Volovik: DOS grows with B
    assert 1.25 < n0[0.2] / n0[0.1] < 2.1           # sub-linear growth (sqrt(B)..linear)


def test_lattice_sc_triangular_and_giant_vortex():
    """The formulation-A lattice supports the triangular (Abrikosov ground-state) cell
    and a multiply-quantized (Vw>1) giant vortex: both keep the core node and a Volovik
    zero-energy DOS, and the Vw=2 cell holds two flux quanta (area doubled)."""
    wc, T = 0.5, 8e-4
    om = E.matsubara(T, wc)
    sq = V.solve_lattice_sc(0.6, T, om, gap_sym='d', field=0.2, lattice='triangular',
                            Ng=16, nbeta=12, itemax=70, mix=0.4, eps=2.5e-3)
    assert abs(sq['S'] - 2 * np.pi / 0.2 * sq['xi'] ** 2) / sq['S'] < 1e-6   # one quantum/cell
    assert sq['absD'].min() / sq['Dbulk'] < 0.3 and sq['absD'].max() / sq['Dbulk'] < 1.05
    n0_tri = float(V.lattice_dos_sc(sq, 'd', np.array([0.0]), nbeta=36, delta=0.03 * sq['Dbulk'])[0])
    assert 0.0 < n0_tri < 1.2                         # finite Volovik DOS
    gv = V.solve_lattice_sc(0.6, T, om, gap_sym='d', field=0.2, Vw=2, Ng=18, nbeta=12,
                            itemax=90, mix=0.4, eps=2.5e-3)
    assert abs(gv['S'] - 2 * (2 * np.pi / 0.2) * gv['xi'] ** 2) / gv['S'] < 1e-6  # two quanta/cell
    assert gv['absD'].min() / gv['Dbulk'] < 0.1      # deep (multiply-quantized) core node


def test_lattice_free_energy_triangular():
    """The Eilenberger free-energy functional (Ichioka-Machida, lattice_free_energy)
    selects the vortex-lattice symmetry: for an isotropic s-wave gap the triangular
    lattice has the lower free energy, F(square) - F(triangular) > 0 (the generic
    Abrikosov result).  [The gap-anisotropy-driven square transition is a small fourfold
    term not resolved by this real-space framework -- see calc_vortex_lattice_symmetry.]"""
    wc, T = 0.5, 8e-4
    om = E.matsubara(T, wc)
    F = {}
    for latt in ('square', 'triangular'):
        st = V.solve_lattice_sc(0.6, T, om, gap_sym='s', field=0.2, lattice=latt,
                                Ng=18, nbeta=24, kappa=None, itemax=150, mix=0.4, eps=1e-3)
        F[latt] = V.lattice_free_energy(st, 0.6, T, om, 's', nbeta=48)
    assert F['square'] < 0 and F['triangular'] < 0       # condensed (free energy lowered)
    assert F['square'] - F['triangular'] > 1e-4          # triangular favored (isotropic)


def test_lattice_symmetry_gap_anisotropy_resolved():
    """The free energy resolves the gap anisotropy that selects the vortex-lattice
    symmetry: for d-wave it depends on the orientation theta0 of the gap relative to the
    lattice (fourfold ~cos4 theta0), enabling the triangular<->square transition scan;
    for the isotropic s-wave it is theta0-independent.  (The full apex+theta0 scan,
    calc_vortex_lattice_symmetry, reproduces d-wave triangular at low field -> square
    near Hc2.)"""
    wc, T = 0.5, 8e-4
    om = E.matsubara(T, wc)
    def Fth(gs, th):
        st = V.solve_lattice_sc(0.6, T, om, gap_sym=gs, field=0.2, lattice='square',
                                Ng=18, nbeta=24, kappa=None, itemax=150, mix=0.4,
                                eps=1e-3, theta0=th)
        return V.lattice_free_energy(st, 0.6, T, om, gs, nbeta=48, theta0=th)
    ds = abs(Fth('s', 0.0) - Fth('s', np.pi / 8))        # isotropic: no theta0 dependence
    dd = abs(Fth('d', 0.0) - Fth('d', np.pi / 8))        # d-wave: fourfold theta0 dependence
    assert ds < 1e-5
    assert dd > 1e-4 and dd > 30 * ds                    # d-wave clearly anisotropic vs s-wave


def test_lattice_symmetry_wannier_fs_rotation():
    """theta0 rotates the WHOLE crystal (Wannier/model FS + gap) rigidly against the vortex
    lattice -- not just the model gap.  Decisive invariance check on the C4-symmetric 'tb'
    FS: the lattice free energy is exactly invariant under a C4 symmetry rotation
    (theta0 -> theta0 + pi/2) yet changes under the non-symmetry pi/4 rotation, for BOTH
    s- and d-wave (the FS velocity anisotropy itself also couples to the lattice -- e.g.
    s-wave square lattices in the borocarbides -- so theta0 feeds FS *and* gap anisotropy
    into the symmetry selection)."""
    wc, T = 0.5, 8e-4
    om = E.matsubara(T, wc)
    # Cost-tuned: the C4 invariance is exact by symmetry, so it survives a coarse
    # discretisation (dsym ~ 1e-15 against the 1e-10 threshold), and the pi/4 signal
    # barely moves with it (dnon = 5.0e-3 / 1.3e-3 for s / d at Nth=48 as at Nth=160).
    fm = E.build_model_fs('tb', Nth=48, mu=-1.0, params=1.0)
    def Fth(gs, th):
        st = V.solve_lattice_sc(0.6, T, om, gap_sym=gs, field=0.2, Ng=12, nbeta=8,
                                kappa=None, itemax=60, mix=0.4, eps=1e-3, theta0=th, fs=fm)
        return V.lattice_free_energy(st, 0.6, T, om, gs, nbeta=16, theta0=th, fs=fm)
    for gs in ('s', 'd'):
        dsym = abs(Fth(gs, 0.0) - Fth(gs, np.pi / 2))    # C4 rotation: must be invariant
        dnon = abs(Fth(gs, 0.0) - Fth(gs, np.pi / 4))    # non-symmetry: must change
        assert dsym < 1e-10
        assert dnon > 1e-4 and dnon > 1e5 * dsym


def test_lattice_sc_finite_kappa_screening():
    """Finite-kappa A(r) back-reaction (je #2 connection) on the formulation-A lattice:
    the London-screened smooth vector potential keeps the core node (the complex Delta
    still winds) yet reduces the Volovik DOS -- weaker screening (larger kappa) gives a
    higher N(0), and all finite-kappa values lie below the bare (uniform-A-overcounting)
    extreme limit."""
    wc, T = 0.5, 8e-4
    om = E.matsubara(T, wc)
    n0 = {}
    for kap in (None, 2.0, 10.0):
        st = V.solve_lattice_sc(0.6, T, om, gap_sym='d', field=0.2, Ng=16, nbeta=14,
                                kappa=kap, itemax=90, mix=0.35, eps=2.5e-3)
        assert st['absD'].min() / st['Dbulk'] < 0.35     # core still suppressed (node)
        n0[kap] = float(V.lattice_dos_sc(st, 'd', np.array([0.0]), nbeta=30,
                                         delta=0.03 * st['Dbulk'])[0])
    assert n0[2.0] < n0[10.0] < n0[None]                 # screening: smaller kappa -> lower N(0)


def test_lattice_sc_dvector_texture():
    """je-style d-vector triplet on the TRUE periodic lattice (2x2 spin-matrix Riccati,
    formulation A): the dominant p_x(e_x) winds with the Abrikosov lattice and nodes at
    every core, while the subdominant p_y(e_z) nucleates core-localized (core > bulk) --
    the vortex-core d-vector texture, now on the periodic cell."""
    wc, T = 0.3, 3e-3
    om = E.matsubara(T, wc)
    # Ng=14, not 10: the core is a single point, so the shallowest grid value is a
    # resolution artefact -- min|A_dom|/Dref converges 0.265 / 0.130 / 0.080 / 0.058 for
    # Ng = 10 / 14 / 20 / 28.  (Ng=10 used to squeak under the threshold only because the
    # matrix Riccati's backward sweep was using the FORWARD generator, which is wrong for
    # the complex Abrikosov-phase gap -- see riccati_gen in feilenberger.f90.)
    st = SP.solve_lattice_sc_dvector((0.7, 0.66), T, om, windings=(1, 0), field=0.1,
                                     lattice='square', Ng=14, nbeta=8, Lchord=5.0,
                                     itemax=35, mix=0.4, eps=5e-3)
    A, Db = st['A'], st['Dbulk']; Dref = max(Db)
    assert Db[0] > 0 and Db[1] < 1e-3 * Db[0]          # dominant supercritical, subdom subcritical
    assert np.abs(A[0]).min() / Dref < 0.15           # dominant: node at the core
    r = np.sqrt(st['X'] ** 2 + st['Y'] ** 2); xi = st['xi']
    core, bulk = r < 0.8 * xi, r > 2.5 * xi
    assert np.abs(A[1])[core].mean() > 1.3 * np.abs(A[1])[bulk].mean()   # subdom core-localized
    # spatially-averaged DOS: finite, recovers to the normal state at high energy
    Nw = SP.lattice_dos_sc_dvector(st, np.array([0.0, 2.0 * Dref]), nbeta=20, delta=0.05 * Dref)
    assert Nw[0] > 0.0 and abs(Nw[1] - 1.0) < 0.25     # N(0)>0 (core states), N(large w)->1


def test_lattice_sc_self_consistent_A():
    """Fully self-consistent vector potential (je A_renew): A(K)=j_{s,T}(K)/K^2 from the
    actual quasiclassical current j_s=<v_F Im g>, solved together with Delta.  It keeps
    the core node, converges, and at weak field (near-uniform |Psi|) reduces to the
    analytic-London A (the current there equals the London (1/lambda^2) p_s)."""
    wc, T = 0.5, 8e-4
    om = E.matsubara(T, wc)
    stL = V.solve_lattice_sc(0.6, T, om, gap_sym='d', field=0.1, Ng=14, nbeta=12,
                             kappa=3.0, itemax=120, mix=0.4, eps=2.5e-3)
    stA = V.solve_lattice_sc(0.6, T, om, gap_sym='d', field=0.1, Ng=14, nbeta=12,
                             kappa=3.0, itemax=120, mix=0.4, eps=2.5e-3, self_consistent_A=True)
    assert stA['Afield'] is not None and stA['iters'] < 120          # converged with an A field
    assert stA['absD'].min() / stA['Dbulk'] < 0.3                    # core node preserved
    nL = float(V.lattice_dos_sc(stL, 'd', np.array([0.0]), nbeta=30, delta=0.03 * stL['Dbulk'])[0])
    nA = float(V.lattice_dos_sc(stA, 'd', np.array([0.0]), nbeta=30, delta=0.03 * stA['Dbulk'])[0])
    assert nA > 0 and abs(nA - nL) / nL < 0.20                       # ~reduces to London at weak field


def test_lattice_sc_grid_convergent():
    """The formulation-A gap map is interpolation-based (per-grid-point anchored
    trajectories), so the converged amplitude is grid-convergent -- unlike a scatter/
    binning estimator whose max|Psi| diverges as the grid is refined."""
    wc, T = 0.4, 1e-3
    om = E.matsubara(T, wc)
    mx = []
    for Ng in (12, 20):
        st = V.solve_lattice_sc(0.6, T, om, gap_sym='d', field=0.2, Ng=Ng, nbeta=12,
                                itemax=60, mix=0.4, eps=3e-3)
        mx.append(st['absD'].max() / st['Dbulk'])
    assert abs(mx[0] - mx[1]) < 0.06                # max|Psi| stable under refinement


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


def test_wannier_fs_matches_model_tb():
    """The real Wannier FS of the 1-orbital square lattice (inputs/square.hop) matches
    the analytic 'tb' model FS: same velocity convention, C4 isotropy, FS shape, and
    (decisively) the same d-wave bulk gap."""
    import libs.plibs as p
    hop = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       'inputs', 'square.hop')
    if not os.path.exists(hop):
        return                                          # skip if the input is absent
    rvec, ham_r, Norb, Nr = p.import_hoppings(hop, 1)
    avec = np.eye(3)
    mu = -1.0
    fw = E.build_wannier_fs(rvec, ham_r, [], avec, mu, mesh=240)
    fm = E.build_model_fs('tb', Nth=240, mu=mu, params=1.0)
    assert abs(fw['nf'].sum() - 1.0) < 1e-12
    # velocity convention: square analytic v_cart = (2 sin kx, 2 sin ky)
    i = len(fw['kx']) // 5
    assert abs(fw['vx'][i] - 2 * np.sin(fw['kx'][i])) < 1e-6
    assert abs(fw['vy'][i] - 2 * np.sin(fw['ky'][i])) < 1e-6
    # FS shape cos kx + cos ky = 0.5, and C4 isotropy
    res = np.cos(fw['kx']) + np.cos(fw['ky'])
    assert abs(res.mean() - 0.5) < 1e-3 and res.std() < 1e-3
    vx2 = (fw['nf'] * fw['vhx'] ** 2).sum()
    vy2 = (fw['nf'] * fw['vhy'] ** 2).sum()
    assert abs(vx2 / vy2 - 1.0) < 1e-3
    # decisive: same d-wave bulk gap as the model FS (same band physics)
    om = E.matsubara(8e-4, 0.5)
    Dw = E.bulk_gap_fs(0.6, 8e-4, om, fw, 'd')
    Dm = E.bulk_gap_fs(0.6, 8e-4, om, fm, 'd')
    assert abs(Dw - Dm) / Dm < 0.01


def test_multiband_gap_projection():
    """Low-energy projection of an orbital-basis pair potential onto the FS bands
    (Nagai-Nakamura JPSJ 85, 074707 (2016), Eq. 43): build_wannier_fs stores the band
    eigenvectors u(k_F), u(-k_F) and project_gap_to_band forms phi=u^dag Delta_orb u*(-k).
    On a 2-orbital model: an intra-orbital equal pairing (Delta=I) gives a band-uniform
    gap; an orbital-selective pairing (diag(1,0)) makes the band gap track the orbital
    weight; diag(1,-1) yields an s+-like sign change driven by the orbital character.
    The projected gap then feeds the quasiclassical solvers (bulk_gap_fs)."""
    import libs.plibs as p
    hop = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       'inputs', 'hop2.input')
    if not os.path.exists(hop):
        return
    rvec, ham_r, Norb, Nr = p.import_hoppings(hop, 1)
    assert Norb == 2
    fw = E.build_wannier_fs(rvec, ham_r, [], np.eye(3), 0.0, mesh=140)
    assert fw['evec'].shape == (len(fw['kx']), 2)            # (a) eigenvectors stored
    assert np.allclose((np.abs(fw['evec']) ** 2).sum(1), 1.0)   # orbital-weight normalized
    # (b) Delta = I  -> band-uniform gap (no orbital selectivity)
    E.project_gap_to_band(fw, np.eye(2))
    assert fw['phi'].real.std() / abs(fw['phi'].real.mean()) < 1e-6
    # (b) Delta = diag(1,0) -> band gap == orbital-0 weight |u_0|^2
    E.project_gap_to_band(fw, np.diag([1.0, 0.0]), normalize=False)
    assert np.allclose(fw['phi'].real, np.abs(fw['evec'][:, 0]) ** 2, atol=1e-6)
    w0 = np.abs(fw['evec'][:, 0]) ** 2
    assert w0.max() - w0.min() > 0.3                          # genuinely orbital-resolved
    # (b) Delta = diag(1,-1) -> s+- sign change from the orbital character
    E.project_gap_to_band(fw, np.diag([1.0, -1.0]), normalize=False)
    assert fw['phi'].real.min() < 0 < fw['phi'].real.max()
    # the projected gap drives the solver: build with gap_orbital=I gives a finite Tc gap
    fwp = E.build_wannier_fs(rvec, ham_r, [], np.eye(3), 0.0, mesh=140, gap_orbital=np.eye(2))
    D = E.bulk_gap_fs(0.6, 8e-4, E.matsubara(8e-4, 0.5), fwp, 0)
    assert D > 0


def test_gap_orbital_from_wannier(tmp_path):
    """The bridge that uses an RPA/FLEX gap exported as Wannier "hopping" (output_gap_wannier)
    as the quasiclassical pairing form factor: gap_orbital_from_wannier loads Delta(R,iw)
    and returns the inverse-FT callable kfrac->Delta_orb(k), which build_wannier_fs projects
    onto the FS bands (Nagai Eq 43).  Validated on a synthetic 1-orbital d-wave gap
    Delta(k)=cos2pi kx - cos2pi ky exported as Delta(R), recovering exactly that on the
    real square FS; n_avg averages consecutive Matsubara slices."""
    import libs.plibs as p
    # 1-orbital d-wave gap Delta(k)=cos2pi kx - cos2pi ky  ->  Delta(R) at (+-1,0),(0,+-1)
    rvec = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]])
    gap = np.zeros((1, 1, 2, 4), dtype=complex)
    gap[0, 0, 0, :] = [0.5, 0.5, -0.5, -0.5]                  # iw_0 slice: d-wave
    gap[0, 0, 1, :] = [0.5, 0.5, 0.5, 0.5]                    # iw_1 slice: extended-s (cos+cos)
    base = str(tmp_path / 'gap_wannier')
    np.savez(base, gap=gap, rvec=rvec, iw=np.array([0.1, 0.3]), temp=1e-3)
    g0 = p.gap_orbital_from_wannier(base)                     # iw_0 (d-wave)
    assert abs(g0(np.array([0.25, 0.0, 0.0]))[0, 0] - (np.cos(2 * np.pi * 0.25) - 1.0)) < 1e-12
    assert abs(g0(np.array([0.0, 0.25, 0.0]))[0, 0] - (1.0 - np.cos(2 * np.pi * 0.25))) < 1e-12
    # n_avg=2 averages the d-wave and extended-s slices -> (cos2pi kx - cos2pi ky + cos+cos)/2 = cos2pi kx
    g_avg = p.gap_orbital_from_wannier(base, iw_index=0, n_avg=2)
    assert abs(g_avg(np.array([0.1, 0.3, 0.0]))[0, 0] - np.cos(2 * np.pi * 0.1)) < 1e-12
    # end-to-end on the real square FS: projected fs['phi'] == the analytic d-wave shape
    hop = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       'inputs', 'square.hop')
    if not os.path.exists(hop):
        return
    rv, ham_r, Norb, Nr = p.import_hoppings(hop, 1)
    fw = E.build_wannier_fs(rv, ham_r, [], np.eye(3), -1.0, mesh=200, gap_orbital=g0)
    assert abs((fw['nf'] * np.abs(fw['phi']) ** 2).sum() - 1.0) < 1e-9   # normalized
    assert fw['phi'].real.min() < -0.3 and fw['phi'].real.max() > 0.3    # d-wave sign change
    fd = E.build_wannier_fs(rv, ham_r, [], np.eye(3), -1.0, mesh=200, gap_sym=1)
    assert abs(np.corrcoef(fw['phi'].real, fd['phi'].real)[0, 1]) > 0.99 # matches analytic d


def test_gap_orbital_from_wannier_round_trips_the_exporter(tmp_path):
    """Pin the Fourier sign against output_gap_wannier itself.

    The d-wave fixture above cannot do this: its Delta(R) is real and inversion
    symmetric, so exp(+2 pi i k.R) and exp(-2 pi i k.R) give the same answer and a flipped
    sign (which would return Delta(-k)) passes unnoticed.  Export a RANDOM gap through
    output_gap_wannier and require the exact Delta(k) that went in -- transposed, because
    the exporter is handed the reversed ctypes view (gap[a,b] = Delta_ba) and stores the
    physical Delta_ab, which is what the reader must hand to project_gap.  The fixture is
    neither symmetric nor Hermitian in the orbital indices, so that transpose is pinned
    here too, not just the Fourier sign."""
    import libs.plibs as p
    import libs.flibs as F
    from libs.plibs._wannier_io import _irr_to_full_kgrid
    Nx = Ny = 8
    Nz, Norb, Nw, temp = 1, 2, 3, 0.02
    klist, kmap, invk = F.gen_irr_k_TRS(Nx, Ny, Nz)
    rng = np.random.default_rng(3)
    gap = (rng.standard_normal((Norb, Norb, Nw, len(klist)))
           + 1j * rng.standard_normal((Norb, Norb, Nw, len(klist))))
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        p.output_gap_wannier(gap, kmap, invk, Nx, Ny, Nz, Nw, temp, N_cut=Nw, zero_tol=0.0)
    finally:
        os.chdir(cwd)

    g0 = p.gap_orbital_from_wannier(str(tmp_path / 'gap_wannier'))
    full = _irr_to_full_kgrid(gap, invk, kmap, Nx, Ny, Nz)      # what the exporter saw
    kfull = kmap / np.array([Nx, Ny, Nz], dtype=np.float64)
    err = wrong = flipped = 0.0
    for i, kf in enumerate(kfull):
        want = full[:, :, 0, kmap[i, 0], kmap[i, 1], kmap[i, 2]].T    # Delta_ab
        got = g0(kf)
        err = max(err, np.abs(got - want).max())
        # what the flipped Fourier sign would have produced at this k: Delta(-k)
        j = np.argmin(np.abs(kfull - (-kf) % 1.0).sum(axis=1))
        wrong = max(wrong, np.abs(full[:, :, 0, kmap[j, 0], kmap[j, 1], kmap[j, 2]].T - want).max())
        flipped = max(flipped, np.abs(want.T - want).max())           # ... and Delta^T(k)
    assert err < 1.0e-12, f'inverse FT does not round-trip the exporter (max err {err:.3e})'
    assert wrong > 1.0e-3, 'fixture too symmetric to distinguish Delta(k) from Delta(-k)'
    assert flipped > 1.0e-3, 'fixture too symmetric to distinguish Delta(k) from Delta^T(k)'
    # an export written before the orbital order was stamped (pair transposed, no
    # 'orb_conv' key) must still read back as the same physical Delta_ab
    z = np.load(str(tmp_path / 'gap_wannier.npz'))
    np.savez(str(tmp_path / 'legacy'), gap=z['gap'].swapaxes(0, 1), rvec=z['rvec'],
             iw=z['iw'], temp=z['temp'])
    g_leg = p.gap_orbital_from_wannier(str(tmp_path / 'legacy'))
    assert max(np.abs(g_leg(kf) - g0(kf)).max() for kf in kfull) < 1.0e-14


def test_gap_color_3d():
    """gap_color_3d (used by main.py's FERMI_3D, color_option=GAP) must reproduce the
    EXACT same phi(k) as the Eilenberger FS pipeline at arbitrary k-points, for both the
    phenomenological gap_sym/delta0 path and the gap_orbital (Nagai projection) path.  On
    the 1-orbital square model the eigenvector is gauge-trivial (u(k)=1), so phi(k) should
    equal the bare form factor exactly."""
    import libs.plibs as p
    hop = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       'inputs', 'square.hop')
    if not os.path.exists(hop):
        return
    rv, ham_r, Norb, Nr = p.import_hoppings(hop, 1)
    centers = [np.array([[0.1, 0.2, 0.0], [0.3, -0.1, 0.0], [-0.2, 0.4, 0.0]])]
    blist = [0]
    # gap_sym (d-wave) path == gap_symms directly
    c_sym = p.gap_color_3d(centers, blist, rv, ham_r, [], gap_sym=1)
    expect = np.cos(2 * np.pi * centers[0][:, 0]) - np.cos(2 * np.pi * centers[0][:, 1])
    assert np.allclose(c_sym[0], expect, atol=1e-12)
    # delta0 rescales the s-wave (gap_sym=0) row uniformly
    c_delta = p.gap_color_3d(centers, blist, rv, ham_r, [], gap_sym=0, delta0=[2.5])
    assert np.allclose(c_delta[0], 2.5, atol=1e-12)
    # gap_orbital callable (1x1 matrix) path: u(k)=u(-k)=1 for a single orbital, so
    # phi(k) == the bare orbital gap, matching the same d-wave form factor as gap_sym=1
    gorb = lambda kfrac: np.array([[np.cos(2 * np.pi * kfrac[0]) - np.cos(2 * np.pi * kfrac[1])]])
    c_orb = p.gap_color_3d(centers, blist, rv, ham_r, [], gap_orbital=gorb)
    assert np.allclose(c_orb[0], expect, atol=1e-12)


def test_gap_sym_index_and_delta0():
    """The gap symmetry can be set by the pyrpa integer gap_sym index (lattice
    harmonics via gap_symms, baked into the FS), and delta0 sets per-band gap
    amplitudes and signs (s+- = opposite signs on different sheets)."""
    import libs.plibs as p
    hop = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       'inputs', 'square.hop')
    if not os.path.exists(hop):
        return
    rvec, ham_r, Norb, Nr = p.import_hoppings(hop, 1)
    avec = np.eye(3)
    mu = -1.0
    fs_s = E.build_wannier_fs(rvec, ham_r, [], avec, mu, mesh=200, gap_sym=0)   # s
    fs_d = E.build_wannier_fs(rvec, ham_r, [], avec, mu, mesh=200, gap_sym=1)   # dx2-y2
    # phi is baked, normalized to <|phi|^2>=1, and returned by fs_form_factor
    assert 'phi' in fs_s
    assert abs((fs_s['nf'] * abs(fs_s['phi']) ** 2).sum() - 1.0) < 1e-9
    assert np.allclose(E.fs_form_factor(fs_s, 0), fs_s['phi'])
    assert fs_s['phi'].real.min() > 0.5                       # s: nodeless (one sign)
    assert fs_d['phi'].real.min() < -0.3 and fs_d['phi'].real.max() > 0.3   # d: sign-changing
    # delta0 per-band ratio + sign on a synthetic two-band FS
    n = len(fs_s['kx'])
    fs2 = E.build_wannier_fs(rvec, ham_r, [], avec, mu, mesh=160)
    m = len(fs2['kx'])
    fs2['band'] = np.where(np.arange(m) < m // 2, 0, 1)
    E.set_fs_gap(fs2, 0, delta0=[1.0, -2.0])
    b0 = fs2['phi'].real[fs2['band'] == 0]
    b1 = fs2['phi'].real[fs2['band'] == 1]
    assert b0.mean() > 0 and b1.mean() < 0                    # opposite signs (s+-)
    assert abs(abs(b1.mean()) / abs(b0.mean()) - 2.0) < 0.05  # |delta0| ratio 2
    assert abs((fs2['nf'] * abs(fs2['phi']) ** 2).sum() - 1.0) < 1e-9


def test_wannier_surface_zebs():
    """The specular surface works on a real Wannier FS (general-FS reflection): on the
    square lattice [100] surface the sign-changing dxy gap is suppressed with a
    zero-energy bound state, while the even d_{x^2-y^2} gap stays flat with none."""
    import libs.plibs as p
    hop = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       'inputs', 'square.hop')
    if not os.path.exists(hop):
        return
    rvec, ham_r, Norb, Nr = p.import_hoppings(hop, 1)
    avec = np.eye(3)
    mu = -1.0
    om = E.matsubara(5e-4, 0.5)
    res = {}
    for gs in (3, 1):                                 # dxy (sign-changing), d (even)
        fw = E.build_wannier_fs(rvec, ham_r, [], avec, mu, mesh=100, gap_sym=gs)
        x, D, Db = S.solve_surface_fs(0.5, 5e-4, om, fw, gs, Lxi=7, nper=6)
        w = np.linspace(-Db, Db, 31)
        ld = S.surface_ldos_fs(fw, D, x, w, gs, ix=0, Dbulk=Db)
        res[gs] = (D[0] / Db, ld[np.argmin(np.abs(w))])
    assert res[3][0] < 0.1 and res[3][1] > 5.0        # dxy: suppressed + ZEBS
    assert res[1][0] > 0.8 and res[1][1] < 1.0        # d: flat + no ZEBS


def test_spin_vortex_lattice_doppler():
    """The self-consistent d-vector vortex lattice (field>0, circular cell): converges
    with the supercurrent Doppler, keeps the core d-vector texture (theta_d ~ 90 deg at
    the core), and the Doppler modifies the inter-vortex subdominant vs the isolated
    vortex."""
    wc, T = 0.2, 1.5e-3
    om = E.matsubara(T, wc)
    prof = {}
    for field in (0.0, 0.5):
        xg, A, Db, xi = SP.solve_vortex2d_dvector((0.8, 0.8 * 0.95), T, om, windings=(1, 0),
                                                  Lxi=7, ngrid=23, nbeta=10, itemax=25, field=field)
        ic = len(xg) // 2
        Dref = np.max(np.abs(Db))
        apx = np.abs(A[0, ic:, ic]) / Dref
        apz = np.abs(A[1, ic:, ic]) / Dref
        thc = np.degrees(np.arctan2(apz[0], max(apx[0], 1e-12)))
        assert apx[0] < 0.1                            # dominant pair-broken in the core
        assert thc > 75.0                              # d-vector ~ pure subdominant at core
        prof[field] = apz
    assert np.abs(prof[0.5] - prof[0.0]).max() > 1e-4   # the supercurrent Doppler modifies the gap


def test_wannier_dvector_vortex():
    """The self-consistent d-vector vortex runs on a real Wannier FS (FS directions,
    Fermi velocities and nf from the band, downsampled to nbeta directions): the
    dominant component is pair-broken in the core and the d-vector tilts to ~90 deg."""
    import libs.plibs as p
    hop = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       'inputs', 'square.hop')
    if not os.path.exists(hop):
        return
    rvec, ham_r, Norb, Nr = p.import_hoppings(hop, 1)
    fw = E.build_wannier_fs(rvec, ham_r, [], np.eye(3), -1.0, mesh=120)   # FS geometry only
    om = E.matsubara(1.5e-3, 0.2)
    xg, A, Db, xi = SP.solve_vortex2d_dvector((0.8, 0.8 * 0.95), 1.5e-3, om, windings=(1, 0),
                                              Lxi=7, ngrid=23, nbeta=16, itemax=22, fs=fw)
    ic = len(xg) // 2
    apx = np.abs(A[0, ic:, ic]); apz = np.abs(A[1, ic:, ic])
    Dref = np.max(np.abs(Db))
    assert apx[0] / Dref < 0.1                          # dominant pair-broken in the core
    assert np.degrees(np.arctan2(apz[0], max(apx[0], 1e-12))) > 75.0   # d ~ pure subdominant at core


def test_vortex_maxwell_selfA():
    """The self-consistent finite-kappa vector potential (Maxwell back-reaction): the
    A_theta(rho) loop converges; the enclosed-flux edge value is fixed (A_theta(Rc) ~
    the extreme type-II value); and as kappa->infinity it reduces to extreme type-II
    (same core DOS)."""
    rg, aphi, n0_sc, n0_ex = V.calc_vortex_maxwell(0.6, 8e-4, 0.5, gap_sym='d', field=0.3,
                                                   kappa=80.0, ngrid=27, nbeta=14, itemax_a=3)
    Rc = rg[-1]
    a_edge_ext = rg[-1] / (2.0 * Rc ** 2)
    assert abs(aphi[-1] - a_edge_ext) / a_edge_ext < 0.1   # edge A fixed by enclosed flux
    assert abs(n0_sc - n0_ex) / n0_ex < 0.05               # kappa->inf reduces to extreme type-II


def test_chiral_pip_gap_sym_minus3():
    """gap_sym = -3 is the chiral p+ip state: a complex form factor (px + i py) that
    is fully gapped (|phi| > 0 everywhere)."""
    from libs.plibs._response import gap_symms
    kl = np.array([[0.25, 0.0, 0.0], [0.0, 0.25, 0.0]])
    row = gap_symms(kl, 1, -3)[0]
    assert np.iscomplexobj(row)
    assert abs(row[0] - 2.0) < 1e-9 and abs(row[1] - 2.0j) < 1e-9   # px + i py
    import libs.plibs as p
    hop = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       'inputs', 'square.hop')
    if not os.path.exists(hop):
        return
    rvec, ham_r, Norb, Nr = p.import_hoppings(hop, 1)
    fw = E.build_wannier_fs(rvec, ham_r, [], np.eye(3), -1.0, mesh=200, gap_sym=-3)
    assert np.iscomplexobj(fw['phi'])
    assert abs((fw['nf'] * abs(fw['phi']) ** 2).sum() - 1.0) < 1e-9
    assert abs(fw['phi']).min() > 0.3                 # chiral: fully gapped (no nodes)


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


def test_scalar_vortex_rejects_chiral():
    """The scalar self-consistent vortex solvers project onto a REAL amplitude
    ansatz, which is exact only for real form factors; a chiral phi (p+ip, ...)
    must be rejected (the d-vector solvers handle multi-component complex
    amplitudes), while real odd-parity harmonics (px/py) stay allowed."""
    T = 2e-3
    om = E.matsubara(T, 0.2)
    for gs in ('p+ip', 'p-ip', -3):
        for call in (lambda: V.solve_vortex2d(0.6, T, om, gap_sym=gs,
                                              ngrid=9, nbeta=4, itemax=1),
                     lambda: V.solve_lattice_sc(0.6, T, om, gap_sym=gs,
                                                field=0.2, Ng=6, nbeta=4, itemax=1)):
            try:
                call()
                raise AssertionError(f"chiral gap_sym={gs!r} must be rejected "
                                     "by the scalar vortex solvers")
            except ValueError:
                pass
    xg, Psi, Db, xi = V.solve_vortex2d(0.6, T, om, gap_sym='px',
                                       ngrid=9, nbeta=4, itemax=1)
    assert np.isfinite(np.abs(Psi)).all()                   # real px still runs


def _tb2(lam):
    """2-orbital 3D model.  lam=0 gives real hoppings, i.e. a time-reversal SYMMETRIC
    normal state H(-k)=conj(H(k)).  lam>0 adds an imaginary on-site inter-orbital term,
    which BREAKS time reversal in the normal state: H(k) is complex with H(-k)=H(k), so
    LAPACK returns u(-k)=u(k) and the legacy 'diag' pair partner is not conj(u(k))."""
    rv, hm = [], []
    for R, M in [((1, 0, 0), [[-1, 0], [0, .6]]), ((-1, 0, 0), [[-1, 0], [0, .6]]),
                 ((0, 1, 0), [[-1, 0], [0, .6]]), ((0, -1, 0), [[-1, 0], [0, .6]]),
                 ((0, 0, 1), [[-.3, 0], [0, -.3]]), ((0, 0, -1), [[-.3, 0], [0, -.3]]),
                 ((0, 0, 0), [[0, 1j * lam], [-1j * lam, .5]])]:
        rv.append(R); hm.append(np.array(M, dtype=complex))
    return np.array(rv, dtype=np.float64), np.array(hm)


def test_pair_gauge_is_gauge_invariant():
    """The band projection pairs |k> with a partner state, and the partner must be given
    in the SAME gauge as u(k).  Two invariants pin this down:

      (a) Delta_orb = identity is orbital blind, so the band gap MUST be band uniform
          (phi = 1) for every k -- whatever the eigenvector phases are.
      (b) phi must not change when each u(k) is multiplied by an arbitrary phase.

    The gauge-fixed routes ('trs', 'soc') satisfy both; the legacy 'diag' route (an
    independently diagonalized u(-k)) fails (a) badly on a complex Hamiltonian -- it
    produces a spurious anisotropy with a sign change on a single sheet -- which is why
    it is not the default.  A model with REAL hoppings cannot detect this at all
    (u(-k)=conj(u(k)) then holds exactly), so this test deliberately uses a complex one."""
    rvec, ham_r = _tb2(0.25)
    fs = E.build_wannier_fs(rvec, ham_r, [], np.eye(3), 0.2, mesh=120)
    assert np.abs(fs['evec'].imag).max() > 0.1                 # genuinely complex u(k)
    # (a) Delta = I -> phi == 1
    E.project_gap_to_band(fs, np.eye(2), normalize=False, gauge='trs')
    assert np.allclose(fs['phi'], 1.0, atol=1e-10)
    E.project_gap_to_band(fs, np.eye(2), normalize=False, gauge='diag')
    assert not np.allclose(fs['phi'], 1.0, atol=1e-3)           # legacy route is wrong
    assert fs['phi'].real.min() < 0 < fs['phi'].real.max()      # ...and even flips sign
    # (b) invariance under an arbitrary per-k phase on u(k)
    D = np.diag([1.0, -1.0])
    ref = E.project_gap(fs['evec'], D, fs['kf'], gauge='trs')
    rng = np.random.default_rng(0)
    ph = np.exp(2j * np.pi * rng.random(len(fs['evec'])))[:, None]
    got = E.project_gap(fs['evec'] * ph, D, fs['kf'], gauge='trs')
    assert np.abs(got - ref).max() < 1e-12
    # 'soc': the time-reversed partner (i sigma_y) conj(u) is gauge invariant too
    sm = E.spin_pair_map(2, 'block')
    r2 = E.project_gap(fs['evec'], D, fs['kf'], gauge='soc', spin_map=sm)
    g2 = E.project_gap(fs['evec'] * ph, D, fs['kf'], gauge='soc', spin_map=sm)
    assert np.abs(g2 - r2).max() < 1e-12
    # the diagnostic sees the discrepancy that makes 'diag' unsafe here
    assert E.check_trs_gauge(fs['evec'], fs['evec_mk'], warn=False) > 1e-6


def test_pair_gauge_survives_a_chiral_order_parameter():
    """A CHIRAL order parameter breaks time reversal in the SUPERCONDUCTING state, but
    the gauge fixing rests on a symmetry of the NORMAL state H(k) -- Delta never enters
    the construction of the pair partner.  So on a time-reversal-symmetric normal state
    the projection stays valid for a chiral gap, and the winding of Delta_orb(k) is
    passed through to phi rather than being destroyed or faked."""
    rvec, ham_r = _tb2(0.0)                                    # real h(R): TRS normal state
    fs = E.build_wannier_fs(rvec, ham_r, [], np.eye(3), 0.2, mesh=360)
    # the partner really is the same band at -k, independently of Delta
    assert E.check_pair_partner(fs['evec'], fs['evec_mk'], warn=False) > 1.0 - 1e-10
    A, kf = 2.0 * np.pi, fs['kf']
    dpid = ((np.cos(A * kf[:, 0]) - np.cos(A * kf[:, 1]))
            + 2j * np.sin(A * kf[:, 0]) * np.sin(A * kf[:, 1]))     # d + i dxy
    gorb = lambda k: np.eye(2) * ((np.cos(A * k[0]) - np.cos(A * k[1]))
                                  + 2j * np.sin(A * k[0]) * np.sin(A * k[1]))
    phi = E.project_gap(fs['evec'], gorb, kf, gauge='trs')
    assert np.abs(phi - dpid).max() < 1e-10                    # orbital-blind -> phi = f(k)
    assert np.abs(phi).min() > 0.5                             # chiral state: no node
    ang = np.unwrap(np.angle(phi[np.argsort(fs['th'])]))
    assert abs(abs(ang[-1] - ang[0]) / (2 * np.pi) - 2.0) < 0.05   # winding 2 (d+id)
    # still gauge invariant with a chiral Delta
    rng = np.random.default_rng(1)
    ph = np.exp(2j * np.pi * rng.random(len(fs['evec'])))[:, None]
    assert np.abs(E.project_gap(fs['evec'] * ph, gorb, kf, gauge='trs') - phi).max() < 1e-10


def test_pair_partner_detects_a_trs_broken_normal_state():
    """The failure mode of the gauge fixing is a time-reversal-broken NORMAL state
    (ferromagnetic superconductor, Zeeman term inside H, Haldane-type hoppings), not a
    chiral gap.  There T|k,n> is not the band-n state at -k, so no gauge-fixed scalar
    partner exists and only |phi| is meaningful -- check_pair_partner sees it."""
    rv_ok, hr_ok = _tb2(0.0)
    rv_bad, hr_bad = _tb2(0.25)
    fs_ok = E.build_wannier_fs(rv_ok, hr_ok, [], np.eye(3), 0.2, mesh=120)
    fs_bad = E.build_wannier_fs(rv_bad, hr_bad, [], np.eye(3), 0.2, mesh=120)
    assert E.check_pair_partner(fs_ok['evec'], fs_ok['evec_mk'], warn=False) > 1.0 - 1e-10
    assert E.check_pair_partner(fs_bad['evec'], fs_bad['evec_mk'], warn=False) < 0.95
    # |phi| survives the ambiguity even there; the phase does not
    D = np.diag([1.0, -1.0])
    a = E.project_gap(fs_bad['evec'], D, fs_bad['kf'], gauge='trs')
    assert np.isfinite(a).all()


def test_pair_gauge_real_model_agrees_with_legacy():
    """For REAL hoppings H(-k)=conj(H(k)) and LAPACK returns u(-k)=conj(u(k)) exactly,
    so the gauge-fixed 'trs' partner and the legacy 'diag' partner coincide: the change
    of default cannot move any published real-hopping number."""
    import libs.plibs as p
    hop = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       'inputs', 'hop2.input')
    if not os.path.exists(hop):
        return
    rvec, ham_r, Norb, Nr = p.import_hoppings(hop, 1)
    fs = E.build_wannier_fs(rvec, ham_r, [], np.eye(3), 0.0, mesh=140)
    assert E.check_trs_gauge(fs['evec'], fs['evec_mk'], warn=False) < 1e-12
    D = np.diag([1.0, -1.0])
    a = E.project_gap(fs['evec'], D, fs['kf'], fs['evec_mk'], gauge='trs')
    b = E.project_gap(fs['evec'], D, fs['kf'], fs['evec_mk'], gauge='diag')
    assert np.abs(a - b).max() < 1e-12


def test_chiral_harmonics_and_eliashberg_seed_guard():
    """gap_symms gains the two chiral combinations: 6 = d+id (in plane) and
    7 = dxz+idyz, whose |phi| vanishes on the whole kz=0 plane -- a HORIZONTAL line
    node, unlike the vertical nodes of the in-plane harmonics.  They are Eilenberger
    form factors only: the linearized Eliashberg kernel and seed are real, so seeding
    it with a chiral harmonic is refused with an explanation instead of failing in
    ctypes."""
    import libs.plibs as p
    k = np.array([[0.2, 0.1, 0.3], [0.2, 0.1, 0.0], [0.11, -0.07, 0.25]])
    d, dxy = p.gap_symms(k, 1, 1)[0], p.gap_symms(k, 1, 3)[0]
    assert np.allclose(p.gap_symms(k, 1, 6)[0], d + 1j * dxy)
    dxz, dyz = p.gap_symms(k, 1, 4)[0], p.gap_symms(k, 1, 5)[0]
    assert np.allclose(p.gap_symms(k, 1, 7)[0], dxz + 1j * dyz)
    assert p.gap_symms(k, 1, 7)[0][1] == 0.0                   # kz=0 -> horizontal node
    assert abs(p.gap_symms(k, 1, 7)[0][0]) > 0.1
    for gs in (6, 7, -3):
        try:
            p.get_initial_gap(k, 1, gs)
        except ValueError as e:
            assert 'chiral' in str(e)
        else:
            raise AssertionError(f'chiral gap_sym={gs} was accepted as an Eliashberg seed')
    # continuum counterpart: d+id is e^{2 i beta}, already <|phi|^2>=1
    beta = np.linspace(0.0, 2.0 * np.pi, 32, endpoint=False)
    assert E._gap_sym_str(6) == 'd+id'
    assert abs(np.mean(np.abs(E.form_factor(beta, 'd+id')) ** 2) - 1.0) < 1e-12
    assert np.abs(E.form_factor(beta, 'd+id') - np.exp(2j * beta)).max() < 1e-12


def test_wannier_fs_kz_stack():
    """build_wannier_fs(nkz>1) stacks k_z slices into a true 3D Fermi surface.  The
    decisive check is the one that motivated it: a k_z-dependent gap (7 = dxz+idyz)
    is IDENTICALLY ZERO on the single k_z=0 slice the FS used to be pinned to, and
    only becomes finite once the k_z stacking is on."""
    import libs.plibs as p
    base = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        'inputs', 'FeS')
    if not os.path.exists(base + '_hr.dat'):
        return
    rvec, ham_r, Norb, Nr = p.import_hoppings(base, 2)
    from libs.plibs._bands import get_emesh
    _, _, eig, _, _ = get_emesh(12, 12, 8, ham_r, [], rvec, np.eye(3), sw_uni=True)
    mu = float(np.quantile(eig.ravel(), 0.6))
    f1 = E.build_wannier_fs(rvec, ham_r, [], np.eye(3), mu, mesh=60, nkz=1, gap_sym=7)
    f3 = E.build_wannier_fs(rvec, ham_r, [], np.eye(3), mu, mesh=60, nkz=12, gap_sym=7)
    assert f1['nkz'] == 1 and f3['nkz'] == 12
    assert np.abs(f1['phi']).max() == 0.0                      # the k_z=0 pathology
    assert np.abs(f3['phi']).max() > 0.1                       # 3D FS: finite gap
    assert len(f3['kx']) > 5 * len(f1['kx'])
    for key in ('kz', 'vz', 'vabs3', 'nkz'):
        assert key in f3
    assert abs(f3['nf'].sum() - 1.0) < 1e-12
    assert abs((f3['nf'] * np.abs(f3['phi']) ** 2).sum() - 1.0) < 1e-10
    assert len(np.unique(f3['kz'])) == 12                      # every slice is present
    # vabs is the IN-PLANE speed (chord velocity and DOS measure); vabs3 the full one
    assert np.allclose(f3['vabs'], np.hypot(f3['vx'], f3['vy']))
    assert (f3['vabs3'] >= f3['vabs'] - 1e-12).all()
    assert np.abs(f3['vz']).max() > 0.0
    # the specular reflection partner stays inside its own k_z slice
    j = S._reflection_index(f3)
    assert np.allclose(f3['kz'][j], f3['kz'])
    assert (f3['band'][j] == f3['band']).all()


def test_build_fs_band_projection_3d():
    """build_fs (the full 3D k-mesh route) can now project an orbital-basis pair
    potential band-resolved onto the Fermi surface, which is what an ACCIDENTAL node
    needs: no gap_symms harmonic carries the orbital character that makes the gap
    depend on k_z and change sign between sheets.  Invariants: Delta=I is band
    uniform; an orbital-selective Delta tracks the orbital weight; delta0 gives the
    phenomenological per-sheet s+- signs the harmonic route used to be unable to."""
    import libs.plibs as p
    base = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        'inputs', 'FeS')
    if not os.path.exists(base + '_hr.dat'):
        return
    rvec, ham_r, Norb, Nr = p.import_hoppings(base, 2)
    from libs.plibs._bands import get_emesh
    _, klist, eig, uni, _ = get_emesh(16, 16, 8, ham_r, [], rvec, np.eye(3), sw_uni=True)
    mu = float(np.quantile(eig.ravel(), 0.6))
    # Delta = I -> band uniform (and hence exactly the normalized constant)
    wf, phi = E.build_fs(eig, klist, mu, 0, 2e-2, uni=uni, gap_orbital=np.eye(Norb))
    assert np.abs(phi.imag).max() < 1e-12
    assert phi.real.std() < 1e-10 and abs(phi.real.mean() - 1.0) < 1e-10
    # orbital selective -> phi == the orbital weight (up to the <|phi|^2>=1 scale)
    sel = np.zeros((Norb, Norb)); sel[0, 0] = 1.0
    wf2, phi2, info = E.build_fs(eig, klist, mu, 0, 2e-2, uni=uni, gap_orbital=sel,
                                 sw_band=True)
    ik, ibnd = info['ik'], info['band']
    u = uni[ik, ibnd, :]
    w0 = np.abs(u[:, 0]) ** 2 / (np.abs(u) ** 2).sum(axis=1)   # normalized orbital weight
    assert np.corrcoef(phi2.real, w0)[0, 1] > 0.999
    assert abs((wf2 * np.abs(phi2) ** 2).sum() / wf2.sum() - 1.0) < 1e-10
    assert set(info) == {'band', 'ik', 'kz'} and len(info['kz']) == len(wf2)
    # delta0: per-band signs on the harmonic route (multiband s+-)
    d0 = np.where(np.arange(Norb) < Norb // 2, 1.0, -1.0)
    _, phid = E.build_fs(eig, klist, mu, 0, 2e-2, delta0=d0)
    assert phid.real.min() < 0 < phid.real.max()
    # gap_orbital without eigenvectors, and the mesh route with the legacy gauge, refuse
    for kw in (dict(gap_orbital=np.eye(Norb)),
               dict(gap_orbital=np.eye(Norb), uni=uni, gauge='diag')):
        try:
            E.build_fs(eig, klist, mu, 0, 2e-2, **kw)
        except ValueError:
            pass
        else:
            raise AssertionError(f'build_fs accepted {sorted(kw)} silently')


def test_model_fs_3d_stacking():
    """build_model_fs gains the same k_z stacking as build_wannier_fs, plus genuinely
    3D dispersions.  Invariants: an IN-PLANE kind stacked over k_z reproduces the
    single-slice cylinder exactly (v_z=0, identical FS averages); a corrugated cylinder
    disperses along k_z and carries the fractional k so the k_z-dependent gap_symms
    harmonics work on it; a closed sheet drops the slices past its pole."""
    f1 = E.build_model_fs('tb', Nth=90, params=1.0)
    f8 = E.build_model_fs('tb', Nth=90, params=1.0, nkz=8)
    assert f1['nkz'] == 1 and f8['nkz'] == 8 and len(f8['kx']) == 8 * len(f1['kx'])
    assert np.abs(f8['vz']).max() == 0.0                    # in-plane kind: no k_z dispersion
    for key in ('vabs', 'vx'):                              # stacking must not bias the average
        assert abs((f8['nf'] * f8[key] ** 2).sum() - (f1['nf'] * f1[key] ** 2).sum()) < 1e-12
    assert abs(f8['nf'].sum() - 1.0) < 1e-12
    # corrugated cylinder: real k_z dispersion, and the lattice harmonics apply
    fc = E.build_model_fs('cyl', Nth=90, nkz=16, params=(1.0, 0.3))
    assert np.abs(fc['vz']).max() > 0.5 and 'kf' in fc
    assert (fc['vabs3'] >= fc['vabs'] - 1e-12).all()
    assert np.abs(E.fs_form_factor(fc, 7)).max() > 0.1      # dxz+idyz needs k_z: finite here
    assert np.abs(E.fs_form_factor(E.build_model_fs('cyl', Nth=90, params=(1.0, 0.3)), 7)).max() == 0.0
    # closed sheet: slices beyond the pole carry no FS and are dropped, not crashed on
    fsph = E.build_model_fs('sphere', Nth=60, nkz=24, kz_max=np.pi / 2)
    assert 0 < fsph['nkz'] < 24
    assert abs(fsph['nf'].sum() - 1.0) < 1e-12
    kk = np.hypot(np.hypot(fsph['kx'], fsph['ky']), fsph['kz'])
    assert kk.std() / kk.mean() < 1e-6                      # it really is a sphere


def test_superfluid_density_zz():
    """The c-axis superfluid density is the direct measure of three-dimensionality.
    On a corrugated cylinder rho_zz is suppressed by (v_z/v_x)^2, giving the mass
    anisotropy lambda_c/lambda_ab; with no k_z dispersion v_z == 0 and rho_zz is
    reported as 0 instead of 0/0.  For an isotropic gap every component shares the
    same temperature dependence, which pins the normalization."""
    T, wc, lam = 3e-4, 0.5, 0.5
    flat = E.build_model_fs('tb', Nth=90, params=1.0, nkz=8)
    D, rxx, ryy, rzz = E.superfluid_density_fs(lam, T, wc, flat, 's')
    assert D > 0 and rxx > 0 and rzz == 0.0                 # no v_z -> no c-axis response
    for tz in (0.3, 0.6):
        fc = E.build_model_fs('cyl', Nth=90, nkz=16, params=(1.0, tz))
        D, rxx, ryy, rzz = E.superfluid_density_fs(lam, T, wc, fc, 's')
        assert abs(rxx - ryy) < 1e-9                        # C4
        assert 0.0 < rzz and abs(rzz - rxx) < 1e-6          # s wave: same T dependence
    # stronger warping -> less anisotropic
    a3 = E.build_model_fs('cyl', Nth=90, nkz=16, params=(1.0, 0.3))
    a6 = E.build_model_fs('cyl', Nth=90, nkz=16, params=(1.0, 0.6))
    aniso = lambda f: np.sqrt((f['nf'] * f['vx'] ** 2).sum() / (f['nf'] * f['vz'] ** 2).sum())
    assert aniso(a3) > aniso(a6) > 1.0


def test_free_energy_selects_the_chiral_partner():
    """The condensation energy is the only thing that can choose between the two members
    of a degenerate pair -- linear theory leaves them degenerate at Tc.  On a 3D FS the
    real partners dxz and dyz must come out EXACTLY degenerate (C4), and the chiral
    combination dxz+idyz must win, having fewer nodes.  Both were unreachable before:
    calc_free_energy converted gap_sym through the 2D continuum map and raised on 4/5/7."""
    fc = E.build_model_fs('cyl', Nth=120, nkz=16, params=(1.0, 0.3))
    T, wc, lam = 5e-5, 0.5, 0.45
    dO4, D4 = E.condensation_energy(lam, T, wc, 4, fs=fc)
    dO5, D5 = E.condensation_energy(lam, T, wc, 5, fs=fc)
    dO7, D7 = E.condensation_energy(lam, T, wc, 7, fs=fc)
    assert D4 > 0 and abs(D4 - D5) < 1e-12 and abs(dO4 - dO5) < 1e-18   # exactly degenerate
    assert dO7 < dO4 and D7 > D4                            # chiral wins, larger gap
    assert dO4 < 0.0
    # the sweep driver accepts the k_z-dependent harmonic and an explicit phi alike
    rows = E.calc_free_energy(lam, T, wc, 7, fs=fc, t_list=np.array([T]),
                              fname=os.devnull)
    assert len(rows) == 1 and abs(rows[0][2] - dO7) < 1e-18
    phi = E.fs_form_factor(fc, 7)
    r2 = E.calc_free_energy(lam, T, wc, phi=phi, wnf=fc['nf'], t_list=np.array([T]),
                            fname=os.devnull)
    assert abs(r2[0][2] - dO7) < 1e-18
    # chiral continuum harmonics reach the cylinder branch too (used to be a KeyError)
    dO, D = E.condensation_energy(lam, T, wc, 'd+id')
    assert D > 0 and dO < 0


def test_pauli_and_spin_accept_the_projected_gap():
    """calc_pauli_limit / calc_spin_pauli used to call build_fs positionally, so the
    band-projected 3D gap could not reach them: Pauli limiting and the singlet-vs-triplet
    Zeeman response were stuck on the analytic harmonics.  They now take the same
    gap_orbital/delta0/gauge arguments as calc_eilenberger, and must reproduce the gap
    that build_fs + solve_gap give for the same projection."""
    import inspect
    import libs.plibs as p
    from libs.plibs import _eilenberger_spin as SP
    for fn in (E.calc_pauli_limit, SP.calc_spin_pauli):
        got = set(inspect.signature(fn).parameters)
        assert {'gap_orbital', 'delta0', 'gauge', 'spin_map'} <= got, fn.__name__
    base = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        'inputs', 'FeS')
    if not os.path.exists(base + '_hr.dat'):
        return
    rvec, ham_r, Norb, Nr = p.import_hoppings(base, 2)
    from libs.plibs._bands import get_emesh
    _, klist, eig, uni, _ = get_emesh(12, 12, 8, ham_r, [], rvec, np.eye(3), sw_uni=True)
    mu = float(np.quantile(eig.ravel(), 0.6))
    D = np.eye(Norb); D[4, 4] = D[Norb - 1, Norb - 1] = -1.0
    T, om = 2e-4, E.matsubara(2e-4, 0.5)
    wf, phi = E.build_fs(eig, klist, mu, 0, 2e-2, uni=uni, gap_orbital=D)
    wf0, phi0 = E.build_fs(eig, klist, mu, 0, 2e-2)
    assert np.abs(phi - phi0).max() > 1e-3                  # the projection really differs
    dref = E.solve_gap(T, wf, phi, om, 0.0, 1e8, 0.45)
    assert dref > 0
    cwd = os.getcwd()
    tmp = os.path.join(cwd, '_pauli_tmp')
    os.makedirs(tmp, exist_ok=True)
    try:
        os.chdir(tmp)
        E.calc_pauli_limit(12, 12, 8, 0.5, ham_r, [], rvec, np.eye(3), mu, T, 0, 0.45,
                           h_list=np.array([0.0]), fs_width=2e-2, gap_orbital=D)
        d0 = float(np.loadtxt('pauli_gap.dat')[1])          # Delta(h=0)/Delta0 == 1
        assert abs(d0 - 1.0) < 1e-6
    finally:
        os.chdir(cwd)
        for f in ('pauli_gap.dat', 'pauli_dos.dat'):
            if os.path.exists(os.path.join(tmp, f)):
                os.remove(os.path.join(tmp, f))
        os.rmdir(tmp)


def test_surface_gap_heals_to_bulk():
    """The far edge of the surface slab MUST sit at the bulk gap: whatever the surface
    does to the order parameter, several coherence lengths away it has to be the
    homogeneous solution.  solve_surface_fs used to fail this on every Fermi surface --
    it took Dbulk from the FULL FS while the kernel summed over the CHORD set, which is
    a different ensemble (near-grazing chords are dropped and the nearest-k_par specular
    map is not a bijection, so the chord weights sum to ~0.96-0.99).  The gap equation is
    exponential in the effective coupling, so that few per cent of missing weight parked
    the whole profile at the wrong level -- flat, but 5-12% low -- instead of healing.

    Both ensembles are now the same one, so the far field is the homogeneous fixed point
    by construction.  The sign-changing case is the sharp one: it must be strongly
    suppressed AT the surface and still heal to 1 away from it."""
    om = E.matsubara(3e-4, 0.5)
    for kind, par in (('iso', None), ('tb', 1.0)):
        fs = E.build_model_fs(kind, Nth=180, params=par)
        for gs in ('s', 'd'):                              # [100] surface: not sign changing
            x, Damp, Db = S.solve_surface_fs(0.5, 3e-4, om, fs, gs, Lxi=8.0, nper=8)
            assert Db > 0
            assert abs(Damp[-1] / Db - 1.0) < 5e-3, (kind, gs, Damp[-1] / Db)
            assert Damp[0] / Db > 0.9                      # no suppression for these
        # sign changing on this surface: suppressed at x=0, still heals to the bulk
        x, Damp, Db = S.solve_surface_fs(0.5, 3e-4, om, fs, 'dxy', Lxi=8.0, nper=8)
        assert Damp[0] / Db < 0.15 and abs(Damp[-1] / Db - 1.0) < 1e-2
        assert (np.diff(Damp) > -1e-12 * Db).all()         # monotonic recovery
    # a 3D (k_z stacked) FS must heal just as well
    fc = E.build_model_fs('cyl', Nth=90, nkz=8, params=(1.0, 0.3))
    x, Damp, Db = S.solve_surface_fs(0.5, 3e-4, om, fc, 's', Lxi=8.0, nper=8)
    assert abs(Damp[-1] / Db - 1.0) < 5e-3


def test_fs_hvf_sets_the_coherence_length():
    """xi = hvf/(pi*Dbulk) fixes the physical size of the slab, so hvf must be the real
    Fermi speed of the FS, not 1.  It must NOT also be handed to the chord kernel: the
    step there is already dx/|v_x| with the true velocity, so passing both double counts
    the speed -- which is why surface_ldos_fs no longer takes an hvf at all."""
    import inspect
    om = E.matsubara(3e-4, 0.5)
    fs = E.build_model_fs('tb', Nth=180, params=1.0)
    hv = E.fs_hvf(fs)
    assert abs(hv - (fs['nf'] * fs['vabs']).sum()) < 1e-12
    assert hv > 1.5                                        # this FS is nowhere near 1
    assert E.fs_hvf(None, 2.5) == 2.5                      # no FS -> the cylinder value
    x, Damp, Db = S.solve_surface_fs(0.5, 3e-4, om, fs, 's', Lxi=8.0, nper=8)
    assert abs(x[-1] - 8.0 * hv / (np.pi * Db)) < 1e-9     # default: xi from <|v_par|>
    x1, _, Db1 = S.solve_surface_fs(0.5, 3e-4, om, fs, 's', Lxi=8.0, nper=8, hvf=1.0)
    assert abs(x1[-1] - 8.0 / (np.pi * Db1)) < 1e-9        # explicit override still honoured
    assert abs(Db1 - Db) < 1e-12                           # ...and it is only a length unit
    assert 'hvf' not in inspect.signature(S.surface_ldos_fs).parameters


def test_reduce_fs_trajectories():
    """The vortex/lattice solvers cost one chord integration per FS point, so a 3D FS
    (10^3-10^4 points against the 24 directions of the model cylinder) has to be reduced.
    They see a trajectory only through (beta, |v_par|, phi), so points are quantized on
    that grid and each cell becomes one weighted direction.  What has to be preserved:
    the weights and the <|phi|^2>=1 convention (or lambda changes meaning), the mean
    speed (it sets xi), and the bulk gap -- convergently, as the bins are refined."""
    fc = E.build_model_fs('cyl', Nth=180, nkz=24, params=(1.0, 0.3))
    om = E.matsubara(3e-4, 0.5)
    n0 = len(fc['vabs'])
    assert n0 > 4000
    for gs in (0, 1):
        Dfull = E.bulk_gap_fs(0.5, 3e-4, om, fc, gs)
        prev = None
        for nb, nv, nph in ((24, 3, 3), (32, 3, 6), (64, 4, 8)):
            r = E.reduce_fs_trajectories(fc, gs, nbeta=nb, nv=nv, nphi=nph, verbose=False)
            assert r['ntraj'] == len(r['nf']) < n0 and r['nfs_full'] == n0
            assert abs(r['nf'].sum() - 1.0) < 1e-12
            assert abs((r['nf'] * np.abs(r['phi']) ** 2).sum() - 1.0) < 1e-12
            assert abs((r['nf'] * r['vabs']).sum() / E.fs_hvf(fc) - 1.0) < 0.02
            assert np.allclose(r['vabs'], np.hypot(r['vx'], r['vy']))
            err = abs(E.bulk_gap_fs(0.5, 3e-4, om, r, gs) / Dfull - 1.0)
            assert err < 0.02
            if prev is not None:
                assert err <= prev + 5e-3                  # refining bins does not diverge
            prev = err
    # phi is part of the key: sheets of OPPOSITE sign that share a direction and a speed
    # must not be merged, or an s+- gap would cancel itself into nothing
    th = np.linspace(0.0, 2.0 * np.pi, 64, endpoint=False)
    two = dict(vx=np.concatenate([np.cos(th)] * 2), vy=np.concatenate([np.sin(th)] * 2),
               vabs=np.ones(2 * len(th)), nf=np.full(2 * len(th), 1.0 / (2 * len(th))),
               phi=np.concatenate([np.ones(len(th)), -np.ones(len(th))]).astype(complex))
    r = E.reduce_fs_trajectories(two, nbeta=16, nv=1, nphi=2, verbose=False)
    assert r['phi'].real.min() < -0.9 and r['phi'].real.max() > 0.9     # both signs survive
    assert abs((r['nf'] * np.abs(r['phi']) ** 2).sum() - 1.0) < 1e-12
    assert abs((r['nf'] * r['phi']).sum()) < 1e-12         # ...and still cancel in the mean


def test_matrix_riccati_complex_gap():
    """The matrix Riccati must be EXACT for a homogeneous gap of any phase: f = Delta/R.
    It was exact only for a REAL one, because the backward (gamma-tilde) sweep used the
    FORWARD generator -- gamma-tilde obeys the same equation with Delta and Delta^dag
    interchanged, and the two coincide only when Delta^dag = Delta.  The error was 15.3%
    at beta = 90 deg (Delta purely imaginary), step- and grid-independent, and it inflated
    the chiral homogeneous fixed point to 1.494 x its bulk value.  The scalar kernel had
    it right all along (it integrates the reverse path with conjg(dd))."""
    from libs import flibs
    sx, sy, sz = SP._pauli()
    S = sz @ (1j * sy); Sd = np.conj(S).T
    D, w, ns = 1.8e-3, 3.6e-3, 32
    R = np.sqrt(w ** 2 + D ** 2)
    for bdeg in (0.0, 30.0, 45.0, 90.0, 137.0, 180.0):
        b = np.radians(bdeg)
        Dp = np.zeros((ns, 1, 2, 2), dtype=np.complex128)
        Dp[:, :] = D * np.exp(1j * b) * S
        _, f = flibs.matrix_riccati_chords(np.full((ns, 1, 1), w, dtype=np.complex128),
                                           np.ascontiguousarray(Dp), 1.0, 300.0, 0.0)
        fa = np.trace(Sd @ f[ns // 2, 0, 0]) / np.trace(Sd @ S)
        assert abs(fa - D * np.exp(1j * b) / R) < 1e-9 * D / R, (bdeg, fa)


def test_chiral_vortex():
    """A CHIRAL vortex: the two components are e^{+-i l beta} on ONE spin matrix, so their
    relative phase IS the chirality -- exactly what the scalar real-amplitude solver
    cannot hold (_reject_chiral_ff) and what the multi-component complex-amplitude solver
    can.  Their windings must differ by 2l or the ansatz is not axisymmetric.

    Two things are checked: the pure chiral state is a FIXED POINT (start on it and it
    stays), and the vortex has the right structure -- the dominant vanishes at the core
    and heals, while the opposite chirality is induced, vanishing BOTH at the core
    (winding m+2l, so ~rho^3) and in the bulk."""
    assert SP.chiral_windings(1, 1) == (1, 3) and SP.chiral_windings(1, 2) == (1, 5)
    ch = SP.chiral_channels(1)
    assert np.allclose(ch[0][2], ch[1][2])                 # ONE spin matrix
    beta = np.linspace(0.0, 2.0 * np.pi, 16, endpoint=False)
    assert np.allclose(np.abs(ch[0][1](beta)), 1.0)
    T, wc, lam = 4e-4, 0.5, 0.5
    om = E.matsubara(T, wc)
    omf = np.concatenate([om, -om])
    Dchi = float(SP._coupled_bulk_gap(np.array([lam]), T, omf,
                                      np.ones((1, 180)), np.full(180, 1.0 / 180))[0])
    assert Dchi > 0
    # (a) pure chiral is a fixed point: seed it exactly, A_- must stay ~0
    ng = 21
    A0 = np.zeros((2, ng, ng), dtype=np.complex128); A0[0] = Dchi
    _, A, _, _ = SP.solve_vortex2d_dvector((lam, lam), T, om, channels=ch, windings=(0, 0),
                                           Dbulk=np.array([Dchi, 0.0]), Lxi=7.0, ngrid=ng,
                                           nbeta=16, itemax=20, A_init=A0)
    ic = ng // 2
    assert abs(abs(A[0, ic, ic]) / Dchi - 1.0) < 0.02      # holds the bulk gap
    assert abs(A[1, ic, ic]) / Dchi < 5e-3                 # no opposite chirality
    # (b) the vortex itself
    ng = 31
    xg, A, _, xi = SP.solve_vortex2d_dvector((lam, lam), T, om, channels=ch,
                                             windings=SP.chiral_windings(1, 1),
                                             Dbulk=np.array([Dchi, 0.0]), Lxi=8.0,
                                             ngrid=ng, nbeta=16, itemax=60)
    ic = ng // 2
    ap = np.abs(A[0, ic:, ic]) / Dchi
    am = np.abs(A[1, ic:, ic]) / Dchi
    assert ap[0] < 0.02 and ap[-1] > 0.9                   # dominant: node at core, heals
    assert (np.diff(ap[:len(ap) // 2]) > -1e-3).all()      # monotonic out of the core
    assert am[0] < 0.02                                    # induced also vanishes at the core
    assert 0.01 < am.max() < 0.3                           # ...but is induced, and stays small
    assert am.argmax() not in (0, len(am) - 1)             # peaks in between


def test_field_direction_frame():
    """A field that is not along c is still a 2D problem: the vortex lines run along B, so
    the order parameter varies in the plane PERPENDICULAR to it.  fs_field_frame supplies
    the trajectory set in that plane.  Invariants: for B || c it reproduces the historical
    (x, y) set exactly; on an isotropic Fermi surface nothing depends on the direction; on
    a quasi-2D one the two axes of the vortex plane differ by <|v_ab|>/<|v_c|>; and the FS
    points whose velocity is PARALLEL to B (zero speed in the plane, infinite chord step)
    are dropped rather than handed to LAPACK."""
    fc = E.build_model_fs('cyl', Nth=90, nkz=24, params=(1.0, 0.3))
    fr = E.fs_field_frame(fc, (0, 0, 1), gap_sym=0)
    # 1e-9, not machine epsilon: with C4 the two rms velocities are equal only to
    # rounding, so the (area-preserving) rescale factors are 1 to about 1e-11, not to 1
    assert np.abs(fr['dirs'] - np.arctan2(fc['vy'], fc['vx'])).max() < 1e-9
    assert np.abs(fr['vabs'] - np.hypot(fc['vx'], fc['vy'])).max() < 1e-9
    assert abs(fr['aniso_ratio'] - 1.0) < 1e-9 and fr['keep'].all()
    assert np.allclose(fr['e1'], [1, 0, 0]) and np.allclose(fr['e2'], [0, 1, 0])
    # quasi-2D with B in-plane: strongly elliptical, and v || B points are dropped
    fx = E.fs_field_frame(fc, (1, 0, 0), gap_sym=0)
    assert fx['aniso_ratio'] > 3.0
    assert not fx['keep'].all() and (~fx['keep']).sum() < 0.01 * len(fx['keep'])
    assert fx['vabs'].min() > 0.0 and abs(fx['nf'].sum() - 1.0) < 1e-12
    # weaker warping -> less elliptical
    f6 = E.fs_field_frame(E.build_model_fs('cyl', Nth=90, nkz=24, params=(1.0, 0.6)),
                          (1, 0, 0), gap_sym=0)
    assert 1.0 < f6['aniso_ratio'] < fx['aniso_ratio']
    # isotropic sphere: no direction dependence
    sp = E.build_model_fs('sphere', Nth=90, nkz=48, kz_max=np.pi / 2)
    r0 = E.fs_field_frame(sp, (0, 0, 1), gap_sym=0)
    for bd in ((1, 0, 0), (1, 1, 1)):
        r = E.fs_field_frame(sp, bd, gap_sym=0)
        assert abs(r['aniso_ratio'] - 1.0) < 0.05
        assert abs(r['hvf_eff'] / r0['hvf_eff'] - 1.0) < 0.02
    # a strictly 2D FS has no orbital response to an in-plane field: refuse, do not fudge
    try:
        E.fs_field_frame(E.build_model_fs('tb', Nth=90, params=1.0), (1, 0, 0), gap_sym=0)
    except ValueError as e:
        assert '2D Fermi' in str(e)
    else:
        raise AssertionError('an in-plane field on a 2D FS was accepted')


def test_vortex_in_plane_field():
    """The vortex itself with the field off the c axis.  The plane perpendicular to B has
    its two axes rescaled by their rms velocities, which is an affine (area-preserving) map,
    so the square grid and the rotate-interpolate machinery still apply and the core comes
    out circular in the SCALED plane -- i.e. elliptical in real space, with the ellipticity
    set by the Fermi-velocity anisotropy.  Checked for the scalar solver and for the chiral
    two-component one."""
    from libs.plibs import _eilenberger_vortex as V
    T, wc, lam = 4e-4, 0.5, 0.5
    om = E.matsubara(T, wc)
    # Cost-tuned FS/grid sizes: every check here is a ratio or a node depth, which
    # converge long before the absolute profiles do (measured margins in the asserts).
    fc = E.build_model_fs('cyl', Nth=36, nkz=8, params=(1.0, 0.3))

    def core(xg, amp, Db, xi, scale):
        ic = len(xg) // 2
        r = xg[ic:] / xi
        c1 = np.interp(0.5, np.abs(amp[ic:, ic]) / Db, r) * scale[0]
        c2 = np.interp(0.5, np.abs(amp[ic, ic:]) / Db, r) * scale[1]
        return c1, c2

    # explicit B || c must be bit-identical to the default
    a = V.solve_vortex2d(lam, T, om, gap_sym=0, fs=fc, ngrid=17, Lxi=6.0, itemax=12)
    b = V.solve_vortex2d(lam, T, om, gap_sym=0, fs=fc, ngrid=17, Lxi=6.0, itemax=12,
                         bdir=(0, 0, 1))
    assert np.abs(a[1] - b[1]).max() == 0.0 and a[3] == b[3]
    for bd, lo, hi in (((0, 0, 1), 0.95, 1.05), ((1, 0, 0), 3.2, 4.4)):
        fr = E.fs_field_frame(fc, bd, gap_sym=0)
        xg, Psi, Db, xi = V.solve_vortex2d(lam, T, om, gap_sym=0, fs=fc, ngrid=21,
                                           Lxi=6.0, itemax=25, bdir=bd)
        assert np.abs(Psi[len(xg) // 2, len(xg) // 2]) / Db < 0.05     # node at the core
        c1, c2 = core(xg, Psi, Db, xi, fr['scale'])
        assert lo < c1 / c2 < hi, (bd, c1 / c2)                        # core ellipticity
        assert abs(c1 / c2 / fr['aniso_ratio'] - 1.0) < 0.08           # ...tracks <v> ratio
    # A FINITE field needs no extra treatment off the c axis: the rescaling of the two
    # axes is area preserving, so v.p_s is form invariant in the scaled plane (phase
    # gradient: v.grad(phi) = v~.grad~(phi); vector potential: A~_i = A_i/c_i gives
    # v.A = v~.A~ with the same curl B).  The check is that an ISOTROPIC Fermi surface
    # gives a field-direction-independent answer at finite field.
    sp = E.build_model_fs('sphere', Nth=32, nkz=16, kz_max=np.pi / 2)
    prof = []
    for bd in ((0, 0, 1), (1, 0, 0)):
        xg, Psi, Db, xi = V.solve_vortex2d(lam, T, om, gap_sym=0, fs=sp, ngrid=17,
                                           itemax=16, field=0.3, bdir=bd)
        ic = len(xg) // 2
        assert np.abs(Psi[ic, ic]) / Db < 0.02                     # node at the core
        prof.append(np.abs(Psi[ic:, ic]) / Db)
    assert np.abs(prof[0] - prof[1]).max() < 0.03


# --------------------------------------------------------------------------- #
#  standalone runner (no pytest required)
# --------------------------------------------------------------------------- #
if __name__ == '__main__':
    import _tools
    sys.exit(_tools.run_standalone(globals()))
