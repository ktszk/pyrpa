#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Regression / physics tests for the thermoelectric transport coefficients
(libs/flibs/_transport.py wrapping libs/src/fcond.f90).

Two independent routes compute the generalized transport tensors:

  * Boltzmann / relaxation-time (calc_Kn):
        K_n[ij] = sum_{k,n} v_i v_j * tau * (-df/de) * (e-mu)^n * w_k
  * Kubo linear response (calc_Lij / calc_Lij_wl):
        L_n(w) = (i/Nk) sum_{k,l,m} v v * [occupation] / [w + (e_m-e_l) + i.delta]
                 with the heat-current vertex carrying (e-mu)^n.

In the dc limit (w->0) the intraband Kubo term reduces to i/(i.delta)=1/delta, so
with a uniform k-weight 1/Nk and a constant relaxation time tau = 1/delta the two
routes must agree term by term:  Re L_n(0) == K_n.  These tests lock that in,
plus the Sommerfeld/Wiedemann-Franz ratio L22/L11 -> (pi^2/3) T^2 shared by both
routes, and the symmetrized interband heat-current vertex (e_l+e_m)/2 - mu that
keeps L22 positive (the physics fix in fcond.f90).

Runs standalone (no pytest needed):  python tests/test_transport.py
Also works under pytest if installed:  pytest tests/test_transport.py
Requires the Fortran library libfmod.so (cd libs && make FC=ifx SL=MKL).
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import libs.flibs as F


# --------------------------------------------------------------------------- #
#  helpers
# --------------------------------------------------------------------------- #
def _band_diagonal_model(Nk, Norb, seed=0, Ewidth=0.6):
    """A random band-diagonal model: eigenvalues + real band-diagonal velocity.

    Returns eig [Nk,Norb], the Kubo velocity vk [Nk,Norb,Norb,3] (purely diagonal,
    so no interband channel) and the Boltzmann velocity veloc [Nk,Norb,3] holding
    the same diagonal entries.
    """
    rng = np.random.default_rng(seed)
    eig = np.sort(rng.uniform(-Ewidth, Ewidth, size=(Nk, Norb)), axis=1).astype(np.float64)
    vdiag = rng.standard_normal((Nk, Norb, 3))            # real band-diagonal v
    vk = np.zeros((Nk, Norb, Norb, 3), dtype=np.complex128)
    for n in range(Norb):
        vk[:, n, n, :] = vdiag[:, n, :]
    veloc = np.ascontiguousarray(vdiag, dtype=np.float64)
    return eig, vk, veloc


def _kubo_reference(eig, vk, ff, mu, w, idelta, temp):
    """Pure-Python mirror of calc_lij, with the symmetrized heat-current vertex
    ebar = (e_l+e_m)/2 - mu.  Index map (C-order python <-> column-major fortran):
    fortran vk(axis, m, l, i)  ==  python vk[i, m, l, axis].

    One expression for every band pair: the occupation factor
    p = (f_l-f_m)/(e_m-e_l) tends to -df/de at degeneracy, so the l=m term is
    the Drude weight and the intra/inter-band crossover is smooth (mirrors
    occ_factor in fmod.f90)."""
    Nk, Norb = eig.shape
    eps = idelta * 1.0e-3
    de_min = max(1.0e-6 * temp, 1.0e-12)
    L11 = np.zeros((3, 3), dtype=np.complex128)
    L12 = np.zeros((3, 3), dtype=np.complex128)
    L22 = np.zeros((3, 3), dtype=np.complex128)
    for i in range(Nk):
        for l in range(Norb):
            for m in range(Norb):
                de = eig[i, m] - eig[i, l]
                if abs(de) < de_min:
                    pocc = ff[i, m] * (1.0 - ff[i, m]) / temp
                else:
                    pocc = (ff[i, l] - ff[i, m]) / de
                ebar = 0.5 * (eig[i, l] + eig[i, m]) - mu
                if l != m and abs(pocc) * max(1.0, ebar * ebar) <= eps:
                    continue
                denom = w + de + 1j * idelta
                for a in range(3):          # fortran row index "k"
                    for b in range(3):      # fortran col index "j"
                        tmp = vk[i, m, l, a] * vk[i, l, m, b] * pocc / denom
                        L11[a, b] += tmp
                        L12[a, b] += tmp * ebar
                        L22[a, b] += tmp * ebar * ebar
    return 1j * L11 / Nk, 1j * L12 / Nk, 1j * L22 / Nk


# --------------------------------------------------------------------------- #
#  Kubo (dc) <-> Boltzmann equivalence
# --------------------------------------------------------------------------- #
def test_kubo_dc_equals_boltzmann_single_band():
    """Single band: in the dc limit Re L_n(0) must equal the Boltzmann K_n
    when tau = 1/delta (constant) and the k-weight is uniform 1/Nk."""
    Nk, Norb = 64, 1
    temp, mu, delta = 0.025, 0.05, 1.0e-3
    eig, vk, veloc = _band_diagonal_model(Nk, Norb, seed=1)
    ff = F.get_ffermi(eig, mu, temp)

    L11, L12, L22 = F.calc_Lij(eig, vk, ff, mu, 0.0, delta, temp)

    kweight = np.full(Nk, 1.0 / Nk, dtype=np.float64)
    tau = np.full((Nk, Norb), 1.0 / delta, dtype=np.float64)
    K0, K1, K2 = F.calc_Kn(eig, veloc, kweight, temp, mu, tau)

    assert np.allclose(L11.real, K0, rtol=1e-6, atol=1e-9)
    assert np.allclose(L12.real, K1, rtol=1e-6, atol=1e-9)
    assert np.allclose(L22.real, K2, rtol=1e-6, atol=1e-9)
    # dc intraband response is purely real (no dissipative imaginary part at w=0)
    assert np.allclose(L11.imag, 0.0, atol=1e-9)


def test_kubo_dc_equals_boltzmann_multiband():
    """Multi-band but band-diagonal velocity (no interband channel): the band loop
    in both routines must still reproduce the term-by-term equivalence."""
    Nk, Norb = 48, 3
    temp, mu, delta = 0.03, -0.02, 1.0e-3
    eig, vk, veloc = _band_diagonal_model(Nk, Norb, seed=2)
    ff = F.get_ffermi(eig, mu, temp)

    L11, L12, L22 = F.calc_Lij(eig, vk, ff, mu, 0.0, delta, temp)

    kweight = np.full(Nk, 1.0 / Nk, dtype=np.float64)
    tau = np.full((Nk, Norb), 1.0 / delta, dtype=np.float64)
    K0, K1, K2 = F.calc_Kn(eig, veloc, kweight, temp, mu, tau)

    assert np.allclose(L11.real, K0, rtol=1e-6, atol=1e-9)
    assert np.allclose(L12.real, K1, rtol=1e-6, atol=1e-9)
    assert np.allclose(L22.real, K2, rtol=1e-6, atol=1e-9)


def test_calc_lij_wl_matches_scalar():
    """The batched calc_Lij_wl must reproduce the single-frequency calc_Lij at
    each frequency on the mesh."""
    Nk, Norb = 32, 2
    temp, mu, delta = 0.02, 0.0, 2.0e-3
    eig, vk, _ = _band_diagonal_model(Nk, Norb, seed=3)
    ff = F.get_ffermi(eig, mu, temp)
    wl = np.linspace(0.0, 0.4, 6)

    L11w, L12w, L22w = F.calc_Lij_wl(eig, vk, ff, mu, wl, delta, temp)
    for iw, w in enumerate(wl):
        a, b, c = F.calc_Lij(eig, vk, ff, mu, float(w), delta, temp)
        assert np.allclose(a, L11w[iw], atol=1e-10)
        assert np.allclose(b, L12w[iw], atol=1e-10)
        assert np.allclose(c, L22w[iw], atol=1e-10)


# --------------------------------------------------------------------------- #
#  Sommerfeld / Wiedemann-Franz ratio (shared by both routes)
# --------------------------------------------------------------------------- #
def test_sommerfeld_wiedemann_franz():
    """For a constant transport function (uniform DOS, unit velocity) the second
    moment ratio L22/L11 -> (pi^2/3) T^2 in the degenerate limit, identically for
    the Kubo and Boltzmann kernels."""
    Nk, Norb = 4001, 1
    temp, mu, delta = 0.01, 0.0, 1.0e-3
    W = 1.0                                              # half band width >> temp
    eig = np.linspace(-W, W, Nk).reshape(Nk, Norb).astype(np.float64)
    vk = np.zeros((Nk, Norb, Norb, 3), dtype=np.complex128)
    vk[:, 0, 0, 0] = 1.0                                 # unit velocity along x
    veloc = np.zeros((Nk, Norb, 3), dtype=np.float64)
    veloc[:, 0, 0] = 1.0
    ff = F.get_ffermi(eig, mu, temp)

    L11, L12, L22 = F.calc_Lij(eig, vk, ff, mu, 0.0, delta, temp)
    kweight = np.full(Nk, 1.0 / Nk, dtype=np.float64)
    tau = np.full((Nk, Norb), 1.0 / delta, dtype=np.float64)
    K0, K1, K2 = F.calc_Kn(eig, veloc, kweight, temp, mu, tau)

    target = (np.pi**2 / 3.0) * temp**2
    assert abs(L22[0, 0].real / L11[0, 0].real - target) / target < 0.02
    assert abs(K2[0, 0] / K0[0, 0] - target) / target < 0.02
    # particle-hole symmetric setup -> vanishing thermopower kernel L12, K1
    assert abs(L12[0, 0].real) / L11[0, 0].real < 1e-3
    assert abs(K1[0, 0]) / K0[0, 0] < 1e-3


# --------------------------------------------------------------------------- #
#  symmetrized interband heat-current vertex (the fcond.f90 fix)
# --------------------------------------------------------------------------- #
def test_interband_vertex_matches_reference():
    """A genuine multi-band model with full (Hermitian) interband velocity matrices:
    the Fortran L11/L12/L22 must reproduce the symmetrized-vertex Python reference."""
    Nk, Norb = 20, 2
    temp, mu, delta = 0.03, 0.0, 5.0e-3
    rng = np.random.default_rng(7)
    eig = np.sort(rng.uniform(-0.3, 0.3, size=(Nk, Norb)), axis=1).astype(np.float64)
    vk = (rng.standard_normal((Nk, Norb, Norb, 3))
          + 1j * rng.standard_normal((Nk, Norb, Norb, 3))).astype(np.complex128)
    for c in range(3):                                  # Hermitian in band indices
        vk[:, :, :, c] = 0.5 * (vk[:, :, :, c] + vk[:, :, :, c].conj().transpose(0, 2, 1))
    ff = F.get_ffermi(eig, mu, temp)

    L11, L12, L22 = F.calc_Lij(eig, vk, ff, mu, 0.0, delta, temp)
    R11, R12, R22 = _kubo_reference(eig, vk, ff, mu, 0.0, delta, temp)

    assert np.allclose(L11, R11, rtol=1e-8, atol=1e-10)
    assert np.allclose(L12, R12, rtol=1e-8, atol=1e-10)
    assert np.allclose(L22, R22, rtol=1e-8, atol=1e-10)


def test_interband_L22_nonnegative():
    """The symmetrized vertex (e_l+e_m)/2 - mu keeps the L22 diagonal >= 0.
    With mu pinned between two well-separated bands the OLD asymmetric weight
    (e_m-mu)(e_l-mu) was negative there; the fixed vertex gives ebar=0 -> >=0."""
    Nk, Norb = 16, 2
    temp, mu, delta = 2.0e-3, 0.0, 1.0e-3
    eig = np.zeros((Nk, Norb), dtype=np.float64)
    eig[:, 0] = -0.1                                     # mu sits at the midpoint
    eig[:, 1] = +0.1
    vk = np.zeros((Nk, Norb, Norb, 3), dtype=np.complex128)
    vk[:, 0, 1, 0] = 1.0                                 # pure interband x-velocity
    vk[:, 1, 0, 0] = 1.0
    ff = F.get_ffermi(eig, mu, temp)

    _, _, L22 = F.calc_Lij(eig, vk, ff, mu, 0.0, delta, temp)
    assert L22[0, 0].real >= -1e-12                      # never spuriously negative
    assert abs(L22[0, 0].real) < 1e-6                    # ebar=0 at the midpoint


# --------------------------------------------------------------------------- #
#  band degeneracies: the intra/inter-band crossover must be smooth
# --------------------------------------------------------------------------- #
def test_lij_continuous_through_band_degeneracy():
    """Two bands split by 2*d: L11, L12, L22 must be continuous as d -> 0.

    The old kernel had a hard |e_m-e_l| < 1e-9 switch plus a |f_l-f_m| > eps
    skip, so a pair split between those two scales lost its whole -df/de weight.
    Sweeping d straight through that window locks the fix in."""
    Nk, Norb = 64, 2
    temp, mu, delta = 0.02, 0.0, 1.0e-2
    rng = np.random.default_rng(11)
    e0 = rng.uniform(-0.3, 0.3, size=Nk)
    vdiag = rng.standard_normal((Nk, 3))
    ref = None
    for d in (0.0, 1e-9, 1e-8, 1e-6, 1e-5, 1e-4):
        eig = np.stack([e0 - d, e0 + d], axis=1).astype(np.float64)
        vk = np.zeros((Nk, Norb, Norb, 3), dtype=np.complex128)
        for n in range(Norb):                              # same velocity in both bands
            vk[:, n, n, :] = vdiag
        ff = F.get_ffermi(eig, mu, temp)
        got = np.array([x[0, 0] for x in F.calc_Lij(eig, vk, ff, mu, 0.0, delta, temp)])
        if ref is None:
            ref = got
            assert abs(ref[0]) > 0.0
        assert np.abs(got - ref).max() < 1e-4 * np.abs(ref).max()


def test_lij_uses_the_pair_vertex_at_an_exact_degeneracy():
    """At an exact degeneracy the l/=m pair must still use the INTERBAND vertex
    v_ml v_lm.  The old code routed every |e_m-e_l| < 1e-9 pair (including
    l/=m) into the intra-band branch with the diagonal vertex v_mm v_mm, which
    double-counted the Drude weight on band-crossing k-points."""
    Nk, Norb = 8, 2
    temp, mu, delta = 0.02, 0.0, 1.0e-2
    eig = np.zeros((Nk, Norb), dtype=np.float64)           # both bands exactly at mu
    vk = np.zeros((Nk, Norb, Norb, 3), dtype=np.complex128)
    vk[:, 0, 0, 0] = 1.0                                   # diagonal velocity
    vk[:, 1, 1, 0] = 1.0
    # interband vertex is zero -> the l/=m pairs must contribute nothing
    ff = F.get_ffermi(eig, mu, temp)
    L11, _, _ = F.calc_Lij(eig, vk, ff, mu, 0.0, delta, temp)
    # only the two l=m Drude terms survive: 2 * 1 * f(1-f)/T / (i*delta) * i
    expect = 2.0 * 0.25 / temp / delta
    assert abs(L11[0, 0].real - expect) < 1e-8 * expect


# --------------------------------------------------------------------------- #
#  transport distribution function: same normalization as calc_Kn
# --------------------------------------------------------------------------- #
def test_tdf_energy_integral_reproduces_calc_kn():
    """int dE (-df/dE) tdf(E) must equal calc_Kn's K0 (same kweight-weighted sum,
    no hidden 1/Nk), independently of Nk.  calc_tdf used to divide by Nk a second
    time, so the tdf route under-counted sigma by exactly a factor Nk."""
    for Nk in (64, 256):
        eig, vk, veloc = _band_diagonal_model(Nk, 1, seed=5, Ewidth=1.0)
        kweight = np.ones(Nk, dtype=np.float64)
        tau = np.ones_like(eig)
        temp, mu = 0.15, 0.0
        K0, _, _ = F.calc_Kn(eig, veloc, kweight, temp, mu, tau)
        Nw = 20000
        tdf = F.calc_tdf(eig, veloc, kweight, tau, Nw)
        dw = (eig.max() - eig.min()) / Nw
        wl = eig.min() + dw * (np.arange(1, Nw + 1) - 0.5)     # bin centres
        df = 0.25 * (1.0 - np.tanh(0.5 * (wl - mu) / temp) ** 2) / temp
        K0t = (df * tdf.T).T.sum(axis=0) * dw
        assert abs(K0t[0, 0] / K0[0, 0] - 1.0) < 5e-3          # no Nk factor left


# --------------------------------------------------------------------------- #
#  inverse effective mass (get_imass0 / get_imassk in fmod.f90)
# --------------------------------------------------------------------------- #
def test_imass_orbital_tensor_matches_analytic_second_derivative():
    """imk0 must equal -sum_R R_a R_b H(R) e^{-i 2pi k.R} in EVERY component.

    The axis symmetry (d2/dk_a dk_b is symmetric) and the orbital Hermiticity
    (M_lm = conj(M_ml)) are independent; the old symmetrization conjugated the
    axis transpose, which corrupted all m/=l entries and left the (a<b, l<m)
    block at zero."""
    rvec = np.array([[1, 1, 0], [-1, -1, 0], [1, 0, 0], [-1, 0, 0],
                     [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], dtype=np.float64)
    ham_r = np.zeros((len(rvec), 2, 2), dtype=np.complex128)
    for n, r in enumerate(rvec):
        ham_r[n, 0, 0] = -1.0
        ham_r[n, 1, 1] = -0.7
    ham_r[0, 0, 1] = 0.3 + 0.2j                       # complex inter-orbital hopping
    ham_r[1, 1, 0] = np.conj(ham_r[0, 0, 1])
    klist = np.array([[0.13, 0.27, 0.0], [0.41, -0.19, 0.33]])
    imk = F.get_imass0(klist, ham_r, rvec)            # numpy [Nk, l, m, b, a]
    for ik, k in enumerate(klist):
        M = np.zeros((3, 3, 2, 2), dtype=complex)     # M[a,b] in fortran (m,l) order
        for r, h in zip(rvec, ham_r):
            ph = np.exp(-2j * np.pi * (k @ r))
            for a in range(3):
                for b in range(3):
                    M[a, b] += -r[a] * r[b] * h.T * ph
        for l in range(2):
            for m in range(2):
                for a in range(3):
                    for b in range(3):
                        assert abs(imk[ik, l, m, b, a] - M[a, b, m, l]) < 1e-12


# --------------------------------------------------------------------------- #
#  weak-field Hall response (fcond.f90 calc_sigma_hall / plibs.hall_coefficient)
# --------------------------------------------------------------------------- #
def _free_electron_cell(gamma_deg, a=4.0):
    """One s-band on a simple lattice whose a and b axes meet at gamma_deg.

    gamma != 90 makes the Cartesian axes non-orthogonal to the lattice, so
    sigma_xy is already nonzero at B=0 -- the case where the old
    R_H = sigma_xy/(sigma_xx*sigma_yy) shortcut breaks down.
    """
    g = np.radians(gamma_deg)
    avec = np.array([[a, 0., 0.], [a*np.cos(g), a*np.sin(g), 0.], [0., 0., a]])
    rvec = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0],
                     [0, -1, 0], [0, 0, 1], [0, 0, -1]], dtype=np.float64)
    ham_r = np.full((6, 1, 1), -0.5, dtype=np.complex128)
    return avec, rvec, ham_r


def _hall_setup(gamma_deg, N=60, fill=0.005, TK=100.):
    import libs.plibs as P
    hbar = 6.582119569e-16                     # eV s
    kb = 8.617333262e-5                        # eV/K
    eC = 1.602176634e-19                       # C
    avec, rvec, ham_r = _free_electron_cell(gamma_deg)
    temp = TK*kb
    Nk, eig, vk, imass, kw = P.get_emesh(N, N, N, ham_r, [], rvec, avec.T*(1e-10/hbar),
                                         sw_veloc=True, sw_mass=True)
    mu = P.calc_mu(eig, Nk, fill, temp)
    tau = eig*0. + 1.0
    return dict(P=P, eig=eig, vk=vk, imass=imass/eC, kw=kw, tau=tau, temp=temp, mu=mu,
                avec=avec, Vuc=abs(np.linalg.det(avec))*1e-30, fill=fill, eC=eC)


def _hall_kernel(s, band_resolved=False):
    """Band-summed (or band-resolved) Hall kernel from fcond.f90."""
    S, _ = F.calc_sigma_hall(s['eig'], s['vk'], s['imass'], s['kw'], s['tau'],
                             s['temp'], s['mu'])
    return S if band_resolved else S.sum(axis=0)


def _multiorbital_model(Norb=3, seed=3, t=0.2):
    """Hermitian multi-orbital model: H(-R) = H(R)^dagger, on-site H(0) = H(0)^dagger.

    Getting this wrong (independent random blocks on +R and -R) makes H(k)
    non-hermitian, and the band-basis inverse mass then comes out ASYMMETRIC --
    which silently breaks every identity below, since d2e/dk_a dk_b is symmetric
    by construction.  Built explicitly so the test cannot regress into that.
    """
    rng = np.random.default_rng(seed)
    nn = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64)
    rvec = np.vstack([[[0., 0., 0.]], nn, -nn])
    ham_r = np.zeros((7, Norb, Norb), dtype=np.complex128)
    h0 = rng.standard_normal((Norb, Norb)) + 1j*rng.standard_normal((Norb, Norb))
    ham_r[0] = 0.5*(h0 + h0.conj().T)
    for i in range(1, 4):
        h = (rng.standard_normal((Norb, Norb)) + 1j*rng.standard_normal((Norb, Norb)))*t
        ham_r[i], ham_r[i+3] = h, h.conj().T
    return rvec, ham_r


def _sigma_hall_scalar_reference(eig, veloc, imass, kweight, tau, temp, mu):
    """Pure-numpy mirror of the textbook one-sided B||z Hall scalar,

        sum_k 0.5 (vx^2 M_yy + vy^2 M_xx - 2 vx vy M_xy) (-df/de) tau^2 w_k,

    band resolved.  Deliberately written in the symmetrized scalar form rather
    than the Levi-Civita contraction fcond.f90 uses, so agreement pins the
    kernel instead of restating it.  Equals -S[...,x,y,z] when M is symmetric.
    """
    df = 0.25*(1. - np.tanh(0.5*(eig - mu)/temp)**2)/temp
    w = df*tau*tau*kweight[:, None]
    vx, vy = veloc[:, :, 0], veloc[:, :, 1]
    integ = 0.5*(vx*vx*imass[:, :, 1, 1] + vy*vy*imass[:, :, 0, 0]
                 - 2.*vx*vy*imass[:, :, 0, 1])
    return (w*integ).sum(axis=0)                              # [Norb]


def test_sigma_hall_matches_the_scalar_reference():
    """calc_sigma_hall's (x,y,z) component must equal minus the textbook scalar,
    band by band, for a k- and band-dependent tau on a non-orthogonal cell.

    Contracting the Levi-Civita form at (a,b,g)=(x,y,z) gives
    -0.5 w (vx^2 M_yy + vy^2 M_xx - 2 vx vy M_xy) term by term, so the two agree
    exactly -- PROVIDED the inverse mass tensor is symmetric.  That makes this a
    live check on get_imassk (and on H(k) hermiticity) as much as on the kernel.
    """
    import libs.plibs as P
    hbar, eC, kb = 6.582119569e-16, 1.602176634e-19, 8.617333262e-5
    Norb = 3
    rvec, ham_r = _multiorbital_model(Norb)
    for gamma in (90., 113.):
        g = np.radians(gamma)
        avec = np.array([[4., 0., 0.], [4*np.cos(g), 4*np.sin(g), 0.], [0., 0., 4.]])
        N, temp = 16, 200.*kb
        Nk, eig, vk, imass, kw = P.get_emesh(N, N, N, ham_r, [], rvec, avec.T*(1e-10/hbar),
                                             sw_veloc=True, sw_mass=True)
        im = imass/eC
        assert np.abs(im - im.transpose(0, 1, 3, 2)).max() <= 1e-10*np.abs(im).max()
        mu = P.calc_mu(eig, Nk, 1.3, temp)
        tau = 1. + 0.5*np.abs(np.sin(eig*7))          # k- and band-dependent
        S, Sabs = F.calc_sigma_hall(eig, vk, im, kw, tau, temp, mu)
        ref = _sigma_hall_scalar_reference(eig, vk, im, kw, tau, temp, mu)
        assert np.abs(S[:, 0, 1, 2] + ref).max() <= 1e-9*np.abs(ref).max()   # per band
        Stot = S.sum(axis=0)
        assert abs(Stot[0, 1, 2] + ref.sum()) <= 1e-9*abs(ref.sum())
        # a Hall conductivity is antisymmetric in its two current indices
        assert np.abs(Stot + Stot.swapaxes(0, 1)).max() <= 1e-10*np.abs(Stot).max()
        # the absolute-value companion bounds the signed sum term by term
        assert np.all(np.abs(S) <= Sabs + 1e-30)
        assert np.abs(Stot).max() > 0.0


def test_hall_coefficient_is_lattice_angle_independent():
    """A nearly empty band is free-electron-like, at every lattice angle.

    R_H must be -1/(n e) with n = gsp*fill/Vuc (the dilute limit), and it cannot
    depend on the angle between the a,b axes, since that is the same electron gas
    described in skewed coordinates.  The old sigma-is-diagonal formula drifted by
    -9/-22/-35% at gamma = 105/120/135 deg because it ignored sigma_xy(B=0); the
    rho^(1) = -sigma^-1 sigma^(1) sigma^-1 form is flat to the non-parabolicity
    residue, which is what the second assertion pins.
    """
    gsp, err = 2.0, []
    for gamma in (90., 105., 120., 135.):
        s = _hall_setup(gamma, N=60, fill=0.005)
        K0, _, _ = F.calc_Kn(s['eig'], s['vk'], s['kw'], s['temp'], s['mu'], s['tau'])
        S = _hall_kernel(s)
        Rh = s['P'].hall_coefficient(K0, S, s['Vuc'], s['kw'].sum(), gsp)
        if gamma == 90.:
            assert np.all(Rh < 0.0)                           # electron-like
            assert np.abs(Rh/Rh[2] - 1.0).max() < 1e-8        # cubic cell: isotropic
        else:
            assert abs(K0[0, 1]) > 0.1*abs(K0[0, 0])          # sigma_xy(B=0) is real here
        err.append(Rh[2]/(-1.0/(gsp*s['fill']/s['Vuc']*s['eC'])) - 1.0)
    err = np.array(err)
    assert np.abs(err).max() < 0.05                           # each angle accurate
    assert np.abs(err - err[0]).max() < 5e-3                  # and the SAME error at every angle


def test_hall_carrier_density_sign_and_magnitude():
    """n_H must be positive for both carrier types, with the type read off sign(R_H)."""
    import libs.plibs as P
    eC = 1.602176634e-19
    n_H, kind = P.hall_carrier_density(np.array([-2.0e-8, +2.0e-8, -1.0e-9]), eC)
    assert kind == ['electron', 'hole', 'electron']
    assert np.all(n_H > 0.0)
    assert abs(n_H[0]/(1.0/(2.0e-8*eC)/1e6) - 1.0) < 1e-12
    assert abs(n_H[0] - n_H[1]) < 1e-6*n_H[0]                 # |R_H| alone sets the magnitude


def _two_pocket_model(a=4.0, t=0.5, tp=0.25, delta=0.4):
    """Two non-crossing bands: an electron pocket at Gamma and a hole pocket at R.

    e1 = -2t*S + 6t   (minimum 0 at Gamma),  e2 = -2tp*S - 6tp + delta
    (maximum delta at R), with S = cos kx + cos ky + cos kz.  e2 < e1 everywhere,
    so band 0 is the hole sheet and band 1 the electron sheet.  fill = 1 puts the
    two pocket volumes in exact compensation.
    """
    avec = np.eye(3)*a
    nn = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0],
                   [0, -1, 0], [0, 0, 1], [0, 0, -1]], dtype=np.float64)
    rvec = np.vstack([[[0., 0., 0.]], nn])
    ham_r = np.zeros((7, 2, 2), dtype=np.complex128)
    ham_r[0] = np.diag([6*t, -6*tp+delta]).astype(np.complex128)
    for i in range(1, 7):
        ham_r[i] = np.diag([-t, -tp]).astype(np.complex128)
    return avec, rvec, ham_r


def test_luttinger_volume_finds_both_pocket_types():
    """fill=1 of the two-pocket model is exactly compensated: n_e must equal n_h.

    Both pockets must also be found as ONE closed pocket each -- a pocket sitting
    on Gamma is split across all eight corners of the [0,1) mesh, so this only
    works if the components are merged across the periodic faces.
    """
    import libs.plibs as P
    kb = 8.617333262e-5
    avec, rvec, ham_r = _two_pocket_model()
    N, temp = 40, 150.*kb
    Nk, eig, kw = P.get_emesh(N, N, N, ham_r, [], rvec, avec.T)
    mu = P.calc_mu(eig, Nk, 1.0, temp)
    Vuc = np.linalg.det(avec)*1e-30
    car = P.fs_carrier_density(eig, mu, Vuc, 2.0, mesh=(N, N, N))
    assert list(car['kind']) == ['hole', 'electron']
    assert car['npocket_h'][0] == 1 and car['npocket_e'][0] == 0
    assert car['npocket_e'][1] == 1 and car['npocket_h'][1] == 0
    assert not car['open_sheet'].any()
    assert car['n_e'] > 0.0
    # residual is the lattice-count discretization of the two pocket volumes
    assert abs(car['n_e']/car['n_h'] - 1.0) < 0.08          # compensation by construction
    # doping breaks it: the hole pocket must close up and only electrons remain
    mu2 = P.calc_mu(eig, Nk, 1.1, temp)
    car2 = P.fs_carrier_density(eig, mu2, Vuc, 2.0, mesh=(N, N, N))
    assert car2['n_h'] == 0.0 and car2['n_e'] > car['n_e']


def test_luttinger_volume_flags_an_open_sheet():
    """A band dispersing along kx only has a percolating sheet, not a pocket."""
    import libs.plibs as P
    N = 24
    kx = np.linspace(0, 1, N, endpoint=False)
    e = np.cos(2*np.pi*kx)[:, None, None]*np.ones((N, N, N))   # flat in ky,kz
    eig = e.reshape(-1, 1)
    car = P.fs_carrier_density(eig, 0.0, 1e-28, 2.0, mesh=(N, N, N))
    assert car['kind'][0] == 'open'
    assert car['open_sheet'][0]
    assert car['n_e'] == 0.0 and car['n_h'] == 0.0


def test_luttinger_volume_of_a_dilute_pocket_matches_the_filling():
    """One nearly empty band: the pocket volume IS the electron count gsp*fill/Vuc."""
    import libs.plibs as P
    kb, gsp, fill = 8.617333262e-5, 2.0, 0.005
    avec, rvec, ham_r = _free_electron_cell(90.)
    N, temp = 40, 100.*kb
    Nk, eig, kw = P.get_emesh(N, N, N, ham_r, [], rvec, avec.T)
    mu = P.calc_mu(eig, Nk, fill, temp)
    Vuc = np.linalg.det(avec)*1e-30
    car = P.fs_carrier_density(eig, mu, Vuc, gsp, mesh=(N, N, N))
    assert car['kind'][0] == 'electron' and car['npocket_e'][0] == 1
    assert abs(car['n_e']/(gsp*fill/(Vuc*1e6)) - 1.0) < 0.05


def test_calc_k0_band_sums_to_calc_Kn():
    """The two-fluid decomposition must be a decomposition, not a reformulation."""
    eig, vk, veloc = _band_diagonal_model(Nk=200, Norb=4, seed=7)
    kw = np.ones(len(eig))
    tau = np.abs(np.random.default_rng(1).standard_normal(eig.shape)) + 0.5
    K0, _, _ = F.calc_Kn(eig, veloc, kw, 0.02, 0.0, tau)
    K0n = F.calc_k0_band(eig, veloc, kw, 0.02, 0.0, tau)
    assert np.abs(K0n.sum(axis=0) - K0).max() < 1e-9*np.abs(K0).max()
    assert np.abs(K0n - K0n.transpose(0, 2, 1)).max() < 1e-12*np.abs(K0n).max()


def test_hall_coefficient_survives_a_singular_conductivity():
    """A strictly 2-D band has sigma_zz = 0, so the full 3x3 sigma cannot be
    inverted -- yet R_H for B||z is perfectly well defined.  The in-plane 2x2
    fallback must reproduce the plain sigma_xy/(sigma_xx sigma_yy) value there,
    and report R_H for the two in-plane field directions as undefined rather than
    as a number.
    """
    import libs.plibs as P
    hbar, eC, kb = 6.582119569e-16, 1.602176634e-19, 8.617333262e-5
    avec = np.eye(3)*4.0
    rvec = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]], dtype=np.float64)
    ham_r = np.full((4, 1, 1), -0.5, dtype=np.complex128)
    N, temp = 30, 100.*kb
    Nk, eig, vk, imass, kw = P.get_emesh(N, N, 4, ham_r, [], rvec, avec.T*(1e-10/hbar),
                                         sw_veloc=True, sw_mass=True)
    mu = P.calc_mu(eig, Nk, 0.02, temp)
    tau = eig*0. + 1.0
    K0, _, _ = F.calc_Kn(eig, vk, kw, temp, mu, tau)
    assert abs(K0[2, 2]) < 1e-6*abs(K0[0, 0])                  # no dispersion along z
    S, _ = F.calc_sigma_hall(eig, vk, imass/eC, kw, tau, temp, mu)
    S = S.sum(axis=0)
    Vuc = np.linalg.det(avec)*1e-30
    Rh = P.hall_coefficient(K0, S, Vuc, kw.sum(), 2.0)
    scal = _sigma_hall_scalar_reference(eig, vk, imass/eC, kw, tau, temp, mu).sum()
    ref = -Vuc*kw.sum()*scal/(2.0*K0[0, 0]*K0[1, 1])
    assert np.isfinite(Rh[2]) and abs(Rh[2]/ref - 1.0) < 1e-10
    assert not np.isfinite(Rh[0]) and not np.isfinite(Rh[1])
    n_H, kind = P.hall_carrier_density(Rh, eC)
    assert kind == ['-', '-', 'electron']


# --------------------------------------------------------------------------- #
#  standalone runner (no pytest required)
# --------------------------------------------------------------------------- #
if __name__ == '__main__':
    import _tools
    sys.exit(_tools.run_standalone(globals()))
