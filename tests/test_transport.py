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
    fortran vk(axis, m, l, i)  ==  python vk[i, m, l, axis]."""
    Nk, Norb = eig.shape
    eps = idelta * 1.0e-3
    L11 = np.zeros((3, 3), dtype=np.complex128)
    L12 = np.zeros((3, 3), dtype=np.complex128)
    L22 = np.zeros((3, 3), dtype=np.complex128)
    for i in range(Nk):
        for l in range(Norb):
            for m in range(Norb):
                de = eig[i, m] - eig[i, l]
                for a in range(3):          # fortran row index "k"
                    for b in range(3):      # fortran col index "j"
                        if abs(de) < 1.0e-9:
                            e1 = eig[i, m] - mu
                            tmp = (vk[i, m, m, a] * vk[i, m, m, b]
                                   * ff[i, m] * (1.0 - ff[i, m])
                                   / (temp * (w + 1j * idelta)))
                            L11[a, b] += tmp
                            L12[a, b] += tmp * e1
                            L22[a, b] += tmp * e1 * e1
                        elif abs(ff[i, l] - ff[i, m]) > eps:
                            ebar = 0.5 * (eig[i, l] + eig[i, m]) - mu
                            tmp = (vk[i, m, l, a] * vk[i, l, m, b]
                                   * (ff[i, l] - ff[i, m])
                                   / (de * (w + de + 1j * idelta)))
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
#  standalone runner (no pytest required)
# --------------------------------------------------------------------------- #
if __name__ == '__main__':
    import _tools
    sys.exit(_tools.run_standalone(globals()))
