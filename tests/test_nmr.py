#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
Regression tests for the SC-state NMR observables (libs/plibs/_nmr.py):
Knight shift K(T) and the spin-lattice relaxation rate 1/(T1 T).

Runs standalone (no pytest needed):  python tests/test_nmr.py
Also works under pytest if installed:  pytest tests/test_nmr.py
Requires the Fortran library libfmod.so (cd libs && make FC=ifx SL=MKL).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import libs.flibs as F
import libs.plibs as P


# --------------------------------------------------------------------------- #
#  helpers
# --------------------------------------------------------------------------- #
def _square_lattice(N=192, mu=-0.6):
    """One-orbital 2D square lattice with t=1, t'=0.4, used for every physics test.

    N=192 puts 9 k-states inside idelta=1e-3 (nmr_scale_check's n_delta), and is
    measurably the coarsest mesh that reproduces the s-wave Hebel-Slichter peak and
    the d-wave T^3 law -- at N=128 (n_delta=3.9) both are absent.  See the table in
    nmr_scale_check.  Norb=1 keeps a bubble at this mesh in the millisecond range.
    """
    Nk, klist = P.gen_klist(N, N, 1, sw_pp=False)
    kx, ky = klist[:, 0], klist[:, 1]
    eps = (-2.0 * (np.cos(2 * np.pi * kx) + np.cos(2 * np.pi * ky))
           - 0.8 * np.cos(2 * np.pi * kx) * np.cos(2 * np.pi * ky))
    hamk = eps.reshape(Nk, 1, 1).astype(np.complex128)
    olist, site = P.get_chi_orb_list(1, [1])
    # U=J=0: bare bubble, so the tests probe the coherence factors and not the RPA
    # denominator, whose Stoner enhancement would add its own mesh sensitivity.
    Smat, _ = F.gen_SCmatrix(olist, site, 0.0, 0.0)
    return hamk, klist, olist, Smat, mu


def _gap(shape, Nk, D0):
    g = np.asarray(shape, dtype=np.float64)
    g = D0 * g / np.abs(g).max()
    return g.reshape(Nk, 1, 1).astype(np.complex128)


def test_nonorthogonal_sc_bubble_and_normal_limit():
    import scipy.linalg as la
    from libs.plibs._response import _prepare_sc_chi0_builder

    _, ks = P.gen_klist(4, 1, 1, sw_pp=False)
    r = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0]], float)
    hr = np.array([[[0.1, 0.2], [0.2, 0.7]],
                   [[0.05, 0.2], [0.03, -0.05]],
                   [[0.05, 0.03], [0.2, -0.05]]], complex)
    sr = np.array([[[1., 0.2], [0.2, 1.]],
                   [[0.05, 0.08], [0.02, 0.04]],
                   [[0.05, 0.02], [0.08, 0.04]]], complex)
    h, s = F.gen_ham(ks, hr, r, Ovl_r=sr)
    eig, uni = F.get_eig(h, s)
    _, model_gap = P.sc_gap_on_mesh(4, 1, 1, hr, sr, r, 1, 0.2)
    band_gap = uni @ model_gap @ uni.conj().transpose(0, 2, 1)
    assert np.allclose(band_gap[:, 0, 1], 0, atol=1e-12)
    assert np.allclose(band_gap[:, 0, 0], band_gap[:, 1, 1], atol=1e-12)
    olist, site = P.get_chi_orb_list(2, [2])
    q, wl, mu, temp, eta = np.array([0.25, 0, 0]), np.array([0.3]), 0.06, 0.12, 0.02
    shift = F.get_qshift(ks, q)
    builder = _prepare_sc_chi0_builder(h, np.zeros_like(h), mu, temp, ks, olist, wl, eta, False, s)
    normal = F.get_chi_irr(uni, eig, F.get_ffermi(eig, mu, temp), shift, olist, wl, eta, temp)
    assert np.allclose(builder(q), normal, atol=1e-11)

    # Independent generalized BdG eigenvectors and Lehmann sum, with a finite gap.
    gap = 0.12*s
    bdg_e, bdg_v = [], []
    for hk, sk, dk in zip(h, s, gap):
        xi = hk-mu*sk
        bdg = np.block([[xi, dk.conj().T], [dk, -xi]])
        metric = la.block_diag(sk, sk)
        en, vec = la.eigh(bdg, metric)
        bdg_e.append(en)
        bdg_v.append(vec.conj())
    bdg_e, bdg_v = np.array(bdg_e), np.array(bdg_v)
    occ = 1/(np.exp(bdg_e/temp)+1)
    ref = np.zeros((len(olist), len(olist)), complex)
    for k in range(len(ks)):
        kp = int(shift[k])-1
        for l in range(4):
            for m in range(4):
                weight = (occ[kp, l]-occ[k, m])/(wl[0]+bdg_e[k, m]-bdg_e[kp, l]+1j*eta)
                b, p = [], []
                for a, bidx in olist-1:
                    b.append(bdg_v[kp, a, l].conjugate()*bdg_v[k, bidx, m])
                    p.append(bdg_v[kp, bidx+2, l].conjugate()*bdg_v[k, a+2, m])
                x = np.array(b)+np.array(p)
                ref += 0.5*weight*np.outer(x, x.conj())/len(ks)
    got = _prepare_sc_chi0_builder(h, gap, mu, temp, ks, olist, wl, eta, False, s)(q)
    assert np.allclose(got[0], ref.T, atol=1e-11)
    Smat, _ = F.gen_SCmatrix(olist, site, 0, 0)
    res = P.nmr_sweep(h, gap, ks, np.array([[0, 0, 0], q]), mu, np.array([0.04]),
                      0.1, Smat, olist, idelta=eta, overlap=s, verbose=False)
    assert np.isfinite(res['K_sc']).all() and np.isfinite(res['R_sc']).all()


def test_nonorthogonal_model_gap_and_fs_scale():
    r = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0]], float)
    h = np.array([[[0.1]], [[-0.3]], [[-0.3]]], complex)
    s = np.array([[[1]], [[0.1]], [[0.1]]], complex)
    k, gap = P.sc_gap_on_mesh(4, 1, 1, h, s, r, 1, [0.2])
    _, metric = F.gen_ham(k, h, r, Ovl_r=s)
    assert np.allclose((gap/metric)[:, 0, 0], 0.2*P.gap_symms(k, 1, 1)[0])
    eig = np.array([[-0.1], [0.05], [0.3], [0.05]])
    result = P.fermi_surface_gap(eig, 0, gap, window=0.2, overlap=metric)
    reference = P.fermi_surface_gap(eig, 0, gap/metric, window=0.2)
    assert np.allclose(result, reference)
    h2 = np.array([[[0.1, 0.2j], [-0.2j, 0.7]]])
    s2 = np.array([[[1., 0.3], [0.3, 1.]]], complex)
    _, u = F.get_eig(h2, s2)
    band_gap = np.array([[0.2], [-0.1]], complex)
    orbital_gap = F.get_band_to_orb_delta(band_gap, u, s2)
    assert np.allclose(u @ orbital_gap @ u.conj().transpose(0, 2, 1), np.diag(band_gap[:, 0]))


def test_scalar_model_gap_preserves_orthogonal_result():
    r = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]], float)
    h = np.array([[[0.1, 0.2], [0.2, 0.7]],
                  [[-0.4, 0], [0, -0.2]], [[-0.4, 0], [0, -0.2]],
                  [[-0.2, 0], [0, -0.3]], [[-0.2, 0], [0, -0.3]]], complex)
    ki, kmap, invk = F.gen_irr_k_TRS(4, 4, 1)
    _, uni = P.get_eigs(ki, h, [], r)
    for sym in (0, 1, -1):
        initial = F.get_initial_delta(P.get_initial_gap(ki, 2, sym), uni, kmap, invk, 1, sym)
        old = F.remap_gap(initial, np.ones(2), invk, sym)[:, :, 0, :].transpose(2, 0, 1)
        old = 0.2*old/np.abs(old).max()
        k, gap = P.sc_gap_on_mesh(4, 4, 1, h, [], r, sym, 0.2)
        assert np.allclose(k, kmap/np.array([4, 4, 1]))
        assert np.allclose(gap, old, atol=1e-12)


def test_nonorthogonal_bdg_path_matches_mesh_gap():
    import scipy.linalg as la
    r = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]], float)
    h = np.array([[[0.1, 0.2], [0.2, 0.7]],
                  [[-0.4, 0.04], [0.01, -0.2]], [[-0.4, 0.01], [0.04, -0.2]],
                  [[-0.2, 0], [0, -0.3]], [[-0.2, 0], [0, -0.3]]], complex)
    s = np.array([[[1, 0.15], [0.15, 1]],
                  [[0.04, 0.01], [0.02, 0.03]], [[0.04, 0.02], [0.01, 0.03]],
                  [[0.02, 0], [0, 0.01]], [[0.02, 0], [0, 0.01]]], complex)
    for sym, amplitude in ((0, 0.2), (1, 0.2), (1, [0.2, -0.1])):
        k, gap = P.sc_gap_on_mesh(4, 4, 1, h, s, r, sym, amplitude)
        ky, path_e = P.bdg_bands_path(h, s, r, 0.05, k, sym, amplitude, Nk_path=3)
        for i, y in enumerate(ky):
            idx = np.flatnonzero(np.all(np.isclose(k, [0, y, 0]), axis=1))[0]
            hk, sk = F.gen_ham(k[idx:idx+1].copy(), h, r, Ovl_r=s)
            xi, d = hk[0]-0.05*sk[0], gap[idx]
            matrix = np.block([[xi, d.conj().T], [d, -xi]])
            assert np.allclose(path_e[i], la.eigvalsh(matrix, la.block_diag(sk[0], sk[0])), atol=1e-11)


# --------------------------------------------------------------------------- #
#  BCS gap interpolation
# --------------------------------------------------------------------------- #
def test_bcs_gap_ratio_limits():
    assert P.bcs_gap_ratio(1.0e-8, 1.0) > 0.999          # T -> 0 saturates at Delta(0)
    assert P.bcs_gap_ratio(1.0, 1.0) == 0.0              # closes exactly at Tc
    assert P.bcs_gap_ratio(1.5, 1.0) == 0.0              # and stays closed above
    assert 0.93 < P.bcs_gap_ratio(0.5, 1.0) < 0.95       # tanh(1.74) at t=1/2
    # monotonically decreasing in T
    t = np.linspace(0.02, 0.99, 40)
    r = np.array([P.bcs_gap_ratio(ti, 1.0) for ti in t])
    assert np.all(np.diff(r) < 0.0)


def test_bcs_temp_c_roundtrip():
    tc = P.bcs_temp_c(0.05)
    assert abs(tc - 0.05 / 1.764) < 1.0e-12
    assert abs(P.bcs_gap_ratio(1.0e-6 * tc, tc) - 1.0) < 1.0e-6


# --------------------------------------------------------------------------- #
#  q-mesh construction and q -> -q folding
# --------------------------------------------------------------------------- #
def test_nmr_qmesh_gamma_first_and_commensurate():
    qlist, qmult = P.nmr_qmesh(32, 32, 4, 8, 8, 2)
    assert np.allclose(qlist[0], 0.0)                    # Knight shift reads qlist[0]
    # every q must land exactly on the k-mesh or get_qshift cannot resolve k+q
    for q, n in zip(qlist.T, (32, 32, 4)):
        assert np.allclose(q * n, np.rint(q * n))


def test_nmr_qmesh_folding_conserves_total_weight():
    for size in ([8, 8, 2], [4, 4, 1], [16, 16, 1]):
        qf, mf = P.nmr_qmesh(32, 32, 4, *size)
        qu, mu_ = P.nmr_qmesh(32, 32, 4, *size, sw_fold=False)
        assert len(qf) < len(qu)                          # folding really removes points
        assert abs(mf.sum() - mu_.sum()) < 1.0e-9         # but conserves the BZ weight
        assert set(np.unique(mf)) <= {1.0, 2.0}


def test_nmr_qmesh_rounds_to_divisor():
    # 8 does not divide 30; a non-divisor stride would leave a short last interval
    # (biasing the BZ average) and break the q -> -q pairing, so it rounds down to 6.
    qlist, qmult = P.nmr_qmesh(30, 30, 1, 8, 8, 1)
    assert abs(qmult.sum() - 36.0) < 1.0e-9
    assert np.allclose(qlist[:, 0] * 30, np.rint(qlist[:, 0] * 30))


def test_nmr_qmesh_folding_reproduces_unfolded_sum():
    """The folded, multiplicity-weighted q-sum must equal the unfolded one."""
    N = 48
    hamk, klist, olist, Smat, mu = _square_lattice(N)
    D0 = 0.05
    tc = P.bcs_temp_c(D0)
    dk = _gap(np.ones(len(klist)), len(klist), D0)
    res = {}
    for tag, fold in (('fold', True), ('raw', False)):
        q, m = P.nmr_qmesh(N, N, 1, 6, 6, 1, sw_fold=fold)
        r = P.nmr_sweep(hamk, dk, klist, q, mu, [0.6 * tc], tc, Smat, olist,
                        idelta=4.0e-3, hf_weight=m * P.hyperfine_form_factor(q),
                        verbose=False)
        res[tag] = float(r['R_ratio'][0])
    assert abs(res['fold'] - res['raw']) < 1.0e-9 * max(1.0, abs(res['raw']))


# --------------------------------------------------------------------------- #
#  loading a self-consistent Eliashberg gap
# --------------------------------------------------------------------------- #
def _write_gap_wannier(tmpdir, Nx=8, Ny=8, Nz=1, Norb=2, Nw=6, temp=0.02, decay=0.3):
    """Export a known Delta(k, iw_n) through output_gap_wannier and return what it
    should read back as: (path, exact full-grid Delta(k,iw_0), klist, Delta(0) target)."""
    import libs.plibs._wannier_io as _wannier_io
    klist, kmap, invk = F.gen_irr_k_TRS(Nx, Ny, Nz)
    rng = np.random.default_rng(7)
    g = (rng.standard_normal((Norb, Norb, len(klist)))
         + 1j * rng.standard_normal((Norb, Norb, len(klist))))
    gap = np.empty((Norb, Norb, Nw, len(klist)), dtype=np.complex128)
    for n in range(Nw):
        # even in w_n, so the w_n^2 fit has an exact answer to find: Delta(0) = g
        wn = (2 * n + 1) * np.pi * temp
        gap[:, :, n, :] = g * (1.0 - decay * wn ** 2)
    cwd = os.getcwd()
    os.chdir(tmpdir)
    try:
        P.output_gap_wannier(gap, kmap, invk, Nx, Ny, Nz, Nw, temp, N_cut=Nw, zero_tol=0.0)
    finally:
        os.chdir(cwd)
    full = _wannier_io._irr_to_full_kgrid(gap, invk, kmap, Nx, Ny, Nz)
    kfull = kmap / np.array([Nx, Ny, Nz], dtype=np.float64)
    exact_iw0 = np.transpose(full[:, :, 0, kmap[:, 0], kmap[:, 1], kmap[:, 2]], (2, 0, 1))
    # the exact w_n -> 0 limit on the full grid, in the same [Nk,Norb,Norb] layout
    g_full = _wannier_io._irr_to_full_kgrid(g[:, :, None, :], invk, kmap, Nx, Ny, Nz)
    exact_w0 = np.transpose(g_full[:, :, 0, kmap[:, 0], kmap[:, 1], kmap[:, 2]], (2, 0, 1))
    return str(tmpdir), exact_iw0, exact_w0, kfull


def test_gap_from_eliashberg_round_trips_the_exported_gap(tmp_path):
    """Reading back what output_gap_wannier wrote must reproduce Delta(k) exactly.
    This pins the Fourier sign; the same convention is checked from the Eilenberger side
    by test_gap_orbital_from_wannier_round_trips_the_exporter."""
    d, exact_iw0, exact_w0, kfull = _write_gap_wannier(tmp_path)
    # single lowest slice: must match Delta(k, iw_0) to machine precision
    got = P.gap_from_eliashberg(kfull, f'{d}/gap_wannier', sw_extrapolate=False,
                                iw_index=0, n_avg=1, verbose=False)
    assert got.shape == exact_iw0.shape
    assert np.abs(got - exact_iw0).max() < 1.0e-10
    # the wrong sign would give Delta(-k); confirm the test can actually tell them apart
    assert np.abs(exact_iw0 - exact_iw0[_minus_k_index(kfull)]).max() > 1.0e-3


def _minus_k_index(kfull):
    """Index permutation k -> -k on a periodic [0,1) mesh."""
    key = {tuple(np.round(k, 9)): i for i, k in enumerate(kfull % 1.0)}
    return np.array([key[tuple(np.round((-k) % 1.0, 9))] for k in kfull])


def test_gap_from_eliashberg_extrapolates_away_the_matsubara_bias(tmp_path):
    """The stored gap is Delta(0)*(1 - c w_n^2); the w_n->0 fit must recover Delta(0)
    exactly, while reading iw_0 keeps a visible O((pi T)^2) error."""
    d, exact_iw0, exact_w0, kfull = _write_gap_wannier(tmp_path)
    fit = P.gap_from_eliashberg(kfull, f'{d}/gap_wannier', sw_extrapolate=True,
                                n_points=4, order=1, verbose=False)
    raw = P.gap_from_eliashberg(kfull, f'{d}/gap_wannier', sw_extrapolate=False,
                                verbose=False)
    scale = np.abs(exact_w0).max()
    assert np.abs(fit - exact_w0).max() < 1.0e-9 * scale        # fit is exact for this model
    assert np.abs(raw - exact_w0).max() > 1.0e-3 * scale        # iw_0 is measurably biased


def test_gap_from_eliashberg_rescales_to_the_fermi_surface(tmp_path):
    """gap_max sets the maximum ON the Fermi surface when eig/mu are supplied, and the
    zone-wide maximum otherwise -- the two differ for an anisotropic gap."""
    d, _, _, kfull = _write_gap_wannier(tmp_path)
    eps = -2.0 * (np.cos(2 * np.pi * kfull[:, 0]) + np.cos(2 * np.pi * kfull[:, 1]))
    eig = np.repeat(eps[:, None], 2, axis=1)
    target = 0.03
    zone = P.gap_from_eliashberg(kfull, f'{d}/gap_wannier', gap_max=target, verbose=False)
    assert abs(np.abs(zone).max() - target) < 1.0e-12
    fs = P.gap_from_eliashberg(kfull, f'{d}/gap_wannier', eig=eig, mu=-0.5,
                               gap_max=target, verbose=False)
    _, gmax_fs, _ = P.fermi_surface_gap(eig, -0.5, fs)
    assert abs(gmax_fs - target) < 1.0e-12
    # unscaled keeps the stored amplitude
    raw = P.gap_from_eliashberg(kfull, f'{d}/gap_wannier', verbose=False)
    assert np.abs(raw).max() > 1.0e-6 and abs(np.abs(raw).max() - target) > 1.0e-6


def test_gap_from_eliashberg_drives_a_sweep(tmp_path):
    """The loaded gap must be a drop-in for the gap_sym one: same layout, and it must
    reproduce the normal state above Tc."""
    N = 48
    hamk, klist, olist, Smat, mu = _square_lattice(N)
    d, _, _, _ = _write_gap_wannier(tmp_path, Nx=8, Ny=8, Nz=1, Norb=1)
    dk = P.gap_from_eliashberg(klist, f'{d}/gap_wannier', gap_max=0.05,
                               eig=F.get_eig(hamk)[0], mu=mu, verbose=False)
    assert dk.shape == (len(klist), 1, 1)
    tc = P.bcs_temp_c(0.05)
    q, m = P.nmr_qmesh(N, N, 1, 4, 4, 1)
    r = P.nmr_sweep(hamk, dk, klist, q, mu, np.array([0.5, 1.2]) * tc, tc, Smat, olist,
                    idelta=4.0e-3, hf_weight=m * P.hyperfine_form_factor(q), verbose=False)
    assert np.isfinite(r['K_ratio']).all() and np.isfinite(r['R_ratio']).all()
    assert abs(r['K_ratio'][1] - 1.0) < 1.0e-12          # T > Tc -> normal state
    assert r['K_ratio'][0] < 0.99                        # T < Tc -> suppressed


# --------------------------------------------------------------------------- #
#  hyperfine form factor
# --------------------------------------------------------------------------- #
def test_hyperfine_form_factor():
    qlist = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.0]])
    # B=0 -> flat weight A^2
    assert np.allclose(P.hyperfine_form_factor(qlist, 2.0, 0.0), 4.0)
    # A = -4B filters out the antiferromagnetic q = (pi,pi) response exactly
    w = P.hyperfine_form_factor(qlist, -4.0, 1.0)
    assert w[0] < 1.0e-20                                 # A(0) = -4B + 4B = 0
    assert w[1] > 1.0                                     # (pi,pi) survives this choice
    # ... and the complementary choice kills (pi,pi) instead
    w2 = P.hyperfine_form_factor(qlist, 4.0, 1.0)
    assert w2[1] < 1.0e-20


# --------------------------------------------------------------------------- #
#  energy-scale diagnostic
# --------------------------------------------------------------------------- #
def test_fermi_surface_gap_sees_through_the_form_factor():
    """max|Delta| over the zone is not the physical scale: a d_x2-y2 form factor is
    nearly blind to a pocket centred at (pi,pi), and Tc must follow the FS value."""
    N = 96
    Nk, klist = P.gen_klist(N, N, 1, sw_pp=False)
    kx, ky = klist[:, 0], klist[:, 1]
    eps = -2.0 * (np.cos(2 * np.pi * kx) + np.cos(2 * np.pi * ky))
    eig = eps.reshape(Nk, 1)
    D0 = 0.05
    d_wave = _gap(np.cos(2 * np.pi * kx) - np.cos(2 * np.pi * ky), Nk, D0)
    s_wave = _gap(np.ones(Nk), Nk, D0)

    # small electron pocket at (pi,pi), where the d-wave form factor nearly vanishes
    rms, gmax, f_low = P.fermi_surface_gap(eig, 2.83, d_wave)
    assert gmax < 0.5 * D0                      # far below the zone-wide maximum
    assert P.bcs_temp_c(gmax) < 0.5 * P.bcs_temp_c(D0)
    # an isotropic gap on the same pocket keeps the full amplitude
    rms_s, gmax_s, _ = P.fermi_surface_gap(eig, 2.83, s_wave)
    assert abs(gmax_s - D0) < 1.0e-12
    assert abs(rms_s - D0) < 1.0e-12
    # a large Fermi surface does reach the d-wave antinodes
    _, gmax_big, _ = P.fermi_surface_gap(eig, -0.4, d_wave)
    assert gmax_big > 0.8 * D0


def test_nmr_scale_check_uses_the_fermi_surface_gap():
    N = 96
    Nk, klist = P.gen_klist(N, N, 1, sw_pp=False)
    kx, ky = klist[:, 0], klist[:, 1]
    eig = (-2.0 * (np.cos(2 * np.pi * kx) + np.cos(2 * np.pi * ky))).reshape(Nk, 1)
    d_wave = _gap(np.cos(2 * np.pi * kx) - np.cos(2 * np.pi * ky), Nk, 0.05)
    # Judged against max|Delta|=0.05 this broadening looks fine; judged against the gap
    # actually present on the (pi,pi) pocket it is over-smeared, and only the second
    # reading is physical.
    without = P.nmr_scale_check(eig, 2.83, 5.0e-3, 0.05, verbose=False)
    with_fs = P.nmr_scale_check(eig, 2.83, 5.0e-3, 0.05, d_wave, verbose=False)
    assert with_fs['gap_fs'] < without['gap_fs']
    assert without['idelta'] < 0.2 * without['gap_fs']       # passes on the zone-wide gap
    assert with_fs['idelta'] >= 0.2 * with_fs['gap_fs']      # fails on the real one
    assert not with_fs['ok']


def test_nmr_scale_check_flags_too_coarse_kmesh():
    hamk, klist, olist, Smat, mu = _square_lattice(16)
    eig, _ = F.get_eig(hamk)
    # idelta far below the level spacing of a 16x16 mesh: the whole 1/T1T ratio would
    # be Lorentzian tail, and reducing idelta further would not help.
    bad = P.nmr_scale_check(eig, mu, 1.0e-6, 0.05, verbose=False)
    assert bad['n_delta'] < P.N_DELTA_MIN
    assert not bad['ok']
    # idelta comparable to the gap: over-smeared, also rejected
    over = P.nmr_scale_check(eig, mu, 0.04, 0.05, verbose=False)
    assert not over['ok']


def test_nmr_scale_check_accepts_valid_hierarchy():
    hamk, klist, olist, Smat, mu = _square_lattice(192)
    eig, _ = F.get_eig(hamk)
    good = P.nmr_scale_check(eig, mu, 1.0e-3, 0.05, verbose=False)
    assert good['n_delta'] >= P.N_DELTA_MIN and good['idelta'] < 0.2 * good['gap0']
    assert good['ok']


# --------------------------------------------------------------------------- #
#  physics: Knight shift
# --------------------------------------------------------------------------- #
def test_knight_shift_singlet_vanishes_and_triplet_survives():
    """K_s/K_n -> 0 for a singlet (Yosida) but stays finite in the channel
    perpendicular to the d-vector (sw_spsym=True)."""
    N = 192
    hamk, klist, olist, Smat, mu = _square_lattice(N)
    Nk = len(klist)
    D0 = 0.05
    tc = P.bcs_temp_c(D0)
    q, m = P.nmr_qmesh(N, N, 1, 4, 4, 1)
    hf = m * P.hyperfine_form_factor(q)
    dk = _gap(np.ones(Nk), Nk, D0)
    temps = np.array([0.12, 1.05]) * tc

    sing = P.nmr_sweep(hamk, dk, klist, q, mu, temps, tc, Smat, olist,
                       idelta=1.0e-3, sw_spsym=False, hf_weight=hf, verbose=False)
    trip = P.nmr_sweep(hamk, dk, klist, q, mu, temps, tc, Smat, olist,
                       idelta=1.0e-3, sw_spsym=True, hf_weight=hf, verbose=False)

    assert sing['K_ratio'][0] < 0.02          # Yosida: gone at T/Tc = 0.12
    assert trip['K_ratio'][0] > 0.80          # transverse channel survives
    # above Tc the gap is identically zero, so both must return exactly the normal state
    assert abs(sing['K_ratio'][1] - 1.0) < 1.0e-12
    assert abs(trip['K_ratio'][1] - 1.0) < 1.0e-12
    assert sing['gap_ratio'][1] == 0.0


def test_knight_shift_line_nodes_stay_finite_at_low_T():
    """A d-wave gap leaves a Fermi surface at the nodes, so K ~ T does not
    collapse the way the fully gapped s-wave one does."""
    N = 192
    hamk, klist, olist, Smat, mu = _square_lattice(N)
    Nk = len(klist)
    kx, ky = klist[:, 0], klist[:, 1]
    D0 = 0.05
    tc = P.bcs_temp_c(D0)
    q, m = P.nmr_qmesh(N, N, 1, 4, 4, 1)
    hf = m * P.hyperfine_form_factor(q)
    temps = np.array([0.15]) * tc

    s = P.nmr_sweep(hamk, _gap(np.ones(Nk), Nk, D0), klist, q, mu, temps, tc, Smat,
                    olist, idelta=1.0e-3, hf_weight=hf, verbose=False)
    d = P.nmr_sweep(hamk, _gap(np.cos(2 * np.pi * kx) - np.cos(2 * np.pi * ky), Nk, D0),
                    klist, q, mu, temps, tc, Smat, olist, idelta=1.0e-3,
                    hf_weight=hf, verbose=False)
    assert d['K_ratio'][0] > 5.0 * s['K_ratio'][0]
    assert 0.05 < d['K_ratio'][0] < 0.6


# --------------------------------------------------------------------------- #
#  physics: 1/(T1 T)
# --------------------------------------------------------------------------- #
def test_hebel_slichter_peak_for_full_gap():
    """A fully gapped singlet must overshoot the normal state just below Tc
    (case-I coherence factor + the BCS density-of-states pile-up)."""
    N = 192
    hamk, klist, olist, Smat, mu = _square_lattice(N)
    Nk = len(klist)
    D0 = 0.05
    tc = P.bcs_temp_c(D0)
    q, m = P.nmr_qmesh(N, N, 1, 4, 4, 1)
    hf = m * P.hyperfine_form_factor(q)
    temps = np.linspace(0.15, 1.02, 10) * tc
    r = P.nmr_sweep(hamk, _gap(np.ones(Nk), Nk, D0), klist, q, mu, temps, tc, Smat,
                    olist, idelta=1.0e-3, hf_weight=hf, verbose=False)
    peak = np.nanmax(r['R_ratio'])
    t_peak = r['t_red'][np.nanargmax(r['R_ratio'])]
    assert peak > 1.3                          # a real coherence peak, not just noise
    assert 0.5 < t_peak < 1.0                  # located just below Tc
    assert r['R_ratio'][0] < 0.5 * peak        # and it collapses on cooling


def test_line_nodes_give_power_law_relaxation():
    """Line nodes -> N(E) ~ E -> 1/T1 ~ T^3, i.e. (1/T1T) ratio ~ T^2 at low T."""
    N = 192
    hamk, klist, olist, Smat, mu = _square_lattice(N)
    Nk = len(klist)
    kx, ky = klist[:, 0], klist[:, 1]
    D0 = 0.05
    tc = P.bcs_temp_c(D0)
    # 16 per axis: coarser meshes (6 per axis -> only 20 folded points) do not resolve
    # the nodal-nodal scattering vectors and flatten the exponent below T^2.
    q, m = P.nmr_qmesh(N, N, 1, 16, 16, 1)
    hf = m * P.hyperfine_form_factor(q)
    temps = np.linspace(0.15, 0.45, 6) * tc
    r = P.nmr_sweep(hamk, _gap(np.cos(2 * np.pi * kx) - np.cos(2 * np.pi * ky), Nk, D0),
                    klist, q, mu, temps, tc, Smat, olist, idelta=1.0e-3,
                    hf_weight=hf, verbose=False)
    slope = np.polyfit(np.log(r['t_red']), np.log(np.abs(r['R_ratio'])), 1)[0]
    # 1/T1 ~ T^(slope+1).  The exponent itself is only good to about +-0.5 here: at low
    # T the nodal quasi-particles occupy four points of the Fermi surface, so chi''
    # collects weight from the few discrete q that connect nodes and a uniform q-mesh
    # samples them erratically (measured on this model, 1/T1 ~ T^1.4 ... T^2.9 as Nq per
    # axis runs 4..24).  Assert the band that excludes a constant and an exponential.
    assert 1.0 < slope < 3.0, f'1/T1 ~ T^{slope+1:.2f}, expected ~T^3'
    # robust part: the ratio grows monotonically with T (a power law, no coherence peak
    # in this window) and dwarfs the fully gapped case at the same temperature.
    assert np.all(np.diff(r['R_ratio']) > 0.0)
    s = P.nmr_sweep(hamk, _gap(np.ones(Nk), Nk, D0), klist, q, mu, temps[:1], tc, Smat,
                    olist, idelta=1.0e-3, hf_weight=hf, verbose=False)
    assert r['R_ratio'][0] > 3.0 * s['R_ratio'][0]


def test_above_tc_reduces_to_the_normal_state():
    """With Delta(T)=0 the SC branch is skipped and both observables must be
    exactly the normal-state ones -- the anchor the whole normalisation rests on."""
    N = 48
    hamk, klist, olist, Smat, mu = _square_lattice(N)
    Nk = len(klist)
    D0 = 0.05
    tc = P.bcs_temp_c(D0)
    q, m = P.nmr_qmesh(N, N, 1, 4, 4, 1)
    r = P.nmr_sweep(hamk, _gap(np.ones(Nk), Nk, D0), klist, q, mu,
                    np.array([1.1, 2.0]) * tc, tc, Smat, olist, idelta=4.0e-3,
                    hf_weight=m * P.hyperfine_form_factor(q), verbose=False)
    assert np.allclose(r['K_sc'], r['K_n'])
    assert np.allclose(r['R_sc'], r['R_n'])
    assert np.allclose(r['K_ratio'], 1.0)
    assert np.allclose(r['R_ratio'], 1.0)


def test_zero_gap_matches_normal_state_below_tc():
    """Delta(0)=0 with a forced Tc: the BdG bubble must reproduce the normal-state
    bubble at every temperature, checking the GG/FF split has no leftover offset."""
    N = 48
    hamk, klist, olist, Smat, mu = _square_lattice(N)
    Nk = len(klist)
    q, m = P.nmr_qmesh(N, N, 1, 4, 4, 1)
    dk = np.zeros((Nk, 1, 1), dtype=np.complex128)
    tc = 0.02
    r = P.nmr_sweep(hamk, dk, klist, q, mu, np.array([0.3, 0.7]) * tc, tc, Smat, olist,
                    idelta=4.0e-3, hf_weight=m * P.hyperfine_form_factor(q), verbose=False)
    # gap_ratio is non-zero (T < Tc) so the SC branch really did run
    assert np.all(r['gap_ratio'] > 0.0)
    assert np.allclose(r['K_ratio'], 1.0, atol=1.0e-8)
    assert np.allclose(r['R_ratio'], 1.0, atol=1.0e-8)


# --------------------------------------------------------------------------- #
#  input validation and I/O
# --------------------------------------------------------------------------- #
def test_write_nmr_dat_roundtrip(tmp_path=None):
    import tempfile
    import pathlib
    N = 48
    hamk, klist, olist, Smat, mu = _square_lattice(N)
    Nk = len(klist)
    D0 = 0.05
    tc = P.bcs_temp_c(D0)
    q, m = P.nmr_qmesh(N, N, 1, 4, 4, 1)
    r = P.nmr_sweep(hamk, _gap(np.ones(Nk), Nk, D0), klist, q, mu,
                    np.linspace(0.3, 1.1, 4) * tc, tc, Smat, olist, idelta=4.0e-3,
                    hf_weight=m * P.hyperfine_form_factor(q), verbose=False)
    d = pathlib.Path(tmp_path) if tmp_path is not None else pathlib.Path(tempfile.mkdtemp())
    fn = d / 'nmr_sc.dat'
    P.write_nmr_dat(r, str(fn))
    data = np.loadtxt(fn)
    assert data.shape == (4, 9)
    assert np.allclose(data[:, 1], r['t_red'], rtol=1.0e-4)
    assert np.allclose(data[:, 5], r['K_ratio'], rtol=1.0e-4)
    assert np.allclose(data[:, 8], r['R_ratio'], rtol=1.0e-4)


# --------------------------------------------------------------------------- #
#  standalone runner (no pytest required)
# --------------------------------------------------------------------------- #
if __name__ == '__main__':
    import _tools
    sys.exit(_tools.run_standalone(globals()))
