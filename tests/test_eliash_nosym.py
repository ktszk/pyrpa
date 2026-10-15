#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
Regression tests for the symmetry-free linearized Eliashberg solver (f_eliash_nosym.f90,
plibs._eliash_nosym: Zeeman field / FFLO pair momentum) and for the SOC solver fixes that
came out of cross-checking it against the no-SOC solver:

  * bare pairing term counted once, V_Delta = -V/2 + V chi V   (lin_eliash_soc, nosym)
  * triplet seed of lin_eliash_soc is a d_z, odd-parity triplet + parity projection
  * conv_delta_orb_to_band_soc is gauge covariant for triplets (gap_sym argument)

  * lin_eliash (no SOC) imposes inversion parity at the k=-k (TRIM) points, and both solvers run
    Arnoldi with a real inner product (their kernels are only real-linear), so the multi-orbital
    triplet (odd parity with an odd-frequency, orbital-antisymmetric admixture) is the exact
    eigenpair of the full kernel

Note on triplets in multi-orbital models: frequency parity is not a symmetry of the kernel there
(inter-orbital terms mix even and odd frequency); the multi-orbital triplet is checked against
the exact kernel projected on its true symmetry sector (Sz=0 triplet, odd parity).

Runs standalone (no pytest needed):  python tests/test_eliash_nosym.py
Also works under pytest if installed:  pytest tests/test_eliash_nosym.py
Requires the Fortran library libfmod.so (cd libs && make FC=ifx SL=MKL).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import _tools as T
import libs.flibs as F
import libs.plibs as P


# --------------------------------------------------------------------------- #
#  helpers
# --------------------------------------------------------------------------- #
_SX = np.array([[0, 1], [1, 0]], dtype=np.complex128)
_SY = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
_SZ = np.diag([1.0, -1.0]).astype(np.complex128)


def _full_mesh(Nx, Ny, Nz=1):
    """kmap of the complete mesh (gen_allk layout) and its fractional k-points."""
    _, kmap, _ = F.gen_irr_k_TRS(Nx, Ny, Nz)
    return kmap, kmap / np.array([Nx, Ny, Nz], dtype=np.float64)


def _rashba_zeeman_hk(klist, t=1.0, tp=-0.2, alpha=0.4, h=(0.0, 0.0, 0.15), shift=(0.0, 0.0)):
    """Spin-1/2 square lattice H(k+shift) = eps + Rashba + h.sigma: no TRS, no inversion."""
    kx = 2*np.pi*(klist[:, 0] + shift[0])
    ky = 2*np.pi*(klist[:, 1] + shift[1])
    eps = -2*t*(np.cos(kx) + np.cos(ky)) - 4*tp*np.cos(kx)*np.cos(ky)
    H = (eps[:, None, None]*np.eye(2) + alpha*(np.sin(ky)[:, None, None]*_SX - np.sin(kx)[:, None, None]*_SY)
         + h[0]*_SX + h[1]*_SY + h[2]*_SZ)
    return np.ascontiguousarray(H.astype(np.complex128))


def _spin_half(U):
    """Spin-orbital index lists and vertex of a single orbital with spin."""
    slist = np.array([1, -1], dtype=np.int64)
    invs = np.array([2, 1], dtype=np.int64)
    olist, site = P.get_chi_orb_list(2, [2])
    return slist, invs, olist, site, F.gen_Vmatrix(olist, slist, site, invs, U, 0.0)


def _nosym_lambda(Gk, olist, Vmat, kmap, Nx, Ny, temp, iQ=(0, 0, 0), m=60, seed=0):
    """chi0 -> V_Delta -> leading lambda of the nosym kernel (random complex seed)."""
    rng = np.random.default_rng(seed)
    Norb, Nkall = Gk.shape[0], Gk.shape[3]
    init = rng.standard_normal((Norb, Norb, Nkall)) + 1j*rng.standard_normal((Norb, Norb, Nkall))
    with T.silence():
        Vd = F.get_vdelta_nosym(F.get_chi0_nosym(Gk, olist, kmap, temp, Nx, Ny, 1), Vmat)
        return F.linearized_eliashberg_nosym(Vd, Gk, init, Vmat, olist, kmap, Nx, Ny, 1, temp, iQ=iQ,
                                             eps=1e-10, arnoldi_m=m)[1]


def _bruteforce_chi0_kernel(GF, VF, olist, temp, kmap, Nx, Ny, iQ):
    """Independent numpy implementation, in Fortran index space ([k, w, a, b] = X_ab(k, iw)),
    on the same circular 2Nw Matsubara grid as the FFT kernels:
      chi(q,v,m,l) = -(T/N) sum_{k,w} G_{l1 m1}(k+q,w+v) G_{m2 l2}(k,w)
      V_Delta = -V/2 + V (1 + chi V)^-1 chi V   (Nyquist slot -V/2)
      F(k,w)  = -G(k+Q,w) Delta(k,w) G^T(-k,-w),  Delta_ab(k,-w) = -Delta_ba(-k-Q,w)
      Delta'_ad(k,w) = -(T/N) sum_{q,v,b,c} V_Delta[q,v,m(c,d),l(a,b)] F_bc(k-q,w-v)
    Returns chi0 (v >= 0 half) and the kernel matrix on the stored w > 0 components."""
    Nk, Nw, Norb, _ = GF.shape
    N2, Nchi = 2*Nw, len(olist)
    ol = olist - 1
    idx = {(kmap[i, 0], kmap[i, 1]): i for i in range(Nk)}
    sh = lambda s, a: np.array([idx[((s*kmap[i, 0] + a[0]) % Nx, (s*kmap[i, 1] + a[1]) % Ny)] for i in range(Nk)])
    kadd = np.array([[idx[((kmap[a, 0] + kmap[b, 0]) % Nx, (kmap[a, 1] + kmap[b, 1]) % Ny)] for b in range(Nk)]
                     for a in range(Nk)])
    ksub = np.array([[idx[((kmap[a, 0] - kmap[b, 0]) % Nx, (kmap[a, 1] - kmap[b, 1]) % Ny)] for b in range(Nk)]
                     for a in range(Nk)])
    km, kq, kmq = sh(-1, (0, 0)), sh(1, iQ), sh(-1, (-iQ[0], -iQ[1]))
    Gf = np.zeros((Nk, N2, Norb, Norb), complex)
    Gf[:, :Nw] = GF
    for p in range(Nw, N2):
        Gf[:, p] = np.conj(np.swapaxes(GF[:, N2 - 1 - p], 1, 2))
    w = temp/Nk
    chi = np.zeros((Nk, N2, Nchi, Nchi), complex)
    pn = (np.arange(N2)[:, None] + np.arange(N2)[None, :]) % N2
    for q in range(Nk):
        Gkq = Gf[kadd[:, q]]
        for n in range(N2):
            Tm = np.einsum('kpac,kpdb->abcd', Gkq[:, pn[n]], Gf)
            chi[q, n] = -w*Tm[ol[None, :, 0], ol[None, :, 1], ol[:, None, 0], ol[:, None, 1]]
    Vd = np.empty_like(chi)
    eye = np.eye(Nchi)
    for q in range(Nk):
        for n in range(N2):
            if n == Nw:
                Vd[q, n] = -0.5*VF
                continue
            c1 = chi[q, n] @ VF
            Vd[q, n] = -0.5*VF + VF @ np.linalg.solve(eye + c1, c1)
    ichi = -np.ones((Norb, Norb), int)
    for l in range(Nchi):
        ichi[ol[l, 0], ol[l, 1]] = l
    Vt = np.zeros((Nk, N2, Norb, Norb, Norb, Norb), complex)            # [q, v, a, b, c, d]
    for a in range(Norb):
        for b in range(Norb):
            for c in range(Norb):
                for d in range(Norb):
                    Vt[:, :, a, b, c, d] = Vd[:, :, ichi[c, d], ichi[a, b]]
    pm = (np.arange(N2)[:, None] - np.arange(N2)[None, :]) % N2
    mw = N2 - 1 - np.arange(N2)
    dim = Nk*Nw*Norb*Norb
    K = np.zeros((dim, dim), complex)
    for col in range(dim):
        D = np.zeros(dim, complex)
        D[col] = 1.0
        D = D.reshape(Nk, Nw, Norb, Norb)
        Df = np.zeros((Nk, N2, Norb, Norb), complex)
        Df[:, :Nw] = D
        for p in range(Nw, N2):
            Df[:, p] = -np.swapaxes(D[kmq, N2 - 1 - p], 1, 2)
        Ff = -np.einsum('kpab,kpbc,kpdc->kpad', Gf[kq], Df, Gf[km][:, mw])
        new = np.empty((Nk, Nw, Norb, Norb), complex)
        for k in range(Nk):
            Fs = Ff[ksub[k]]
            for p in range(Nw):
                new[k, p] = -w*np.einsum('qnabcd,qnbc->ad', Vt, Fs[:, pm[p]])
        K[:, col] = new.reshape(-1)
    return chi[:, :Nw], K


def _two_orbital_dxz_dyz(t1=0.5, t2=0.05, t4=0.1):
    """dxz/dyz-like two-orbital model with inter-orbital hybridization -4 t4 sin kx sin ky."""
    H = {}

    def add(R, a, b, v):
        H.setdefault(R, np.zeros((2, 2), complex))[a, b] += v
    add((0, 0), 0, 0, 0.0)
    for s in (1, -1):
        add((s, 0), 0, 0, -t1); add((0, s), 0, 0, -t2)
        add((s, 0), 1, 1, -t2); add((0, s), 1, 1, -t1)
    for (rx, ry), v in (((1, 1), t4), ((-1, -1), t4), ((1, -1), -t4), ((-1, 1), -t4)):
        add((rx, ry), 0, 1, v); add((rx, ry), 1, 0, v)
    rvec = np.array([[*R, 0] for R in H], dtype=np.float64)
    return rvec, np.ascontiguousarray(np.array([H[R] for R in H]))


def _trs_solvers(rvec, ham_r, Nx, Ny, Nw, temp, mu, U, J, gs, eps=1e-10, m=80):
    """lambda (and gap) of lin_eliash (no SOC) and lin_eliash_soc (spin-doubled) for one gap_sym."""
    Norb = ham_r.shape[1]
    klist, kmap, invk = F.gen_irr_k_TRS(Nx, Ny, 1)
    olist, site = P.get_chi_orb_list(Norb, [Norb])
    eig, uni = P.get_eigs(klist, ham_r, [], rvec)
    Gk = F.gen_Green0(eig, uni, mu, temp, Nw)
    Smat, Cmat = F.gen_SCmatrix(olist, site, U, J)
    ham2, _, slist, invs = P.spin_double_model(ham_r, [])
    ol2, site2 = P.spin_double_chiolist(olist, site, Norb)
    eig2, uni2 = P.get_eigs(klist, ham2, [], rvec)
    Gk2 = F.gen_Green0(eig2, uni2, mu, temp, Nw)
    Vmat = F.gen_Vmatrix(ol2, slist, site2, invs, U, J)
    with T.silence():
        chi, _ = F.get_chi0(Smat, Cmat, Gk, olist, kmap, invk, temp, Nx, Ny, 1)
        _, lam_n = F.linearized_eliashberg(chi, Gk, uni, P.get_initial_gap(klist, Norb, gs), Smat, Cmat, olist,
                                           np.ones(Norb), kmap, invk, Nx, Ny, 1, temp, gs, eps=eps, arnoldi_m=m)
        chis, sg, sg2, ischi = F.get_chi0_soc(Vmat, Gk2, ol2, slist, kmap, invk, invs, temp, Nx, Ny, 1)
        d_s, lam_s = F.linearized_eliashberg_soc(chis, Gk2, uni2, P.get_initial_gap(klist, 2*Norb, gs), Vmat, sg,
                                                 sg2, np.ones(2*Norb), slist, ol2, kmap, invk, invs, ischi,
                                                 Nx, Ny, 1, temp, gs, eps=eps, arnoldi_m=m)
    return lam_n, lam_s, dict(delta=d_s, uni=uni2, invk=invk, invs=invs, slist=slist, kmap=kmap)


# --------------------------------------------------------------------------- #
#  tests
# --------------------------------------------------------------------------- #
def test_nosym_chi0_and_kernel_match_bruteforce_without_symmetry():
    """Rashba + Zeeman (no TRS, no inversion): chi0 and the leading eigenvalue of the
    nosym kernel equal an independent numpy implementation, also at finite Q."""
    Nx = Ny = 4
    Nw, temp, mu, U = 4, 0.1, -0.3, 2.0
    kmap, klist = _full_mesh(Nx, Ny)
    eig, uni = F.get_eig(_rashba_zeeman_hk(klist))
    Gk = F.gen_Green0(eig, uni, mu, temp, Nw)
    slist, invs, olist, site, Vmat = _spin_half(U)
    with T.silence():
        chi = F.get_chi0_nosym(Gk, olist, kmap, temp, Nx, Ny, 1)
    GF, VF = Gk.transpose(3, 2, 1, 0), Vmat.T
    for iQ in ((0, 0), (1, 2)):
        chib, K = _bruteforce_chi0_kernel(GF, VF, olist, temp, kmap, Nx, Ny, iQ)
        assert np.abs(chi.transpose(3, 2, 1, 0) - chib).max() < 1e-13
        ev = np.linalg.eigvals(K)
        lam_ref = ev[np.argmax(ev.real)]
        assert abs(lam_ref.imag) < 1e-12
        lam = _nosym_lambda(Gk, olist, Vmat, kmap, Nx, Ny, temp, iQ=(*iQ, 0), m=80)
        assert abs(lam - lam_ref.real) < 1e-11, (iQ, lam, lam_ref)


def test_singlet_lambda_agrees_between_nosoc_soc_and_nosym_on_anisotropic_band():
    """tx != ty mixes d and s, so the on-site pair amplitude is finite and a double-counted
    bare term (-V instead of -V/2) would split lin_eliash_soc from lin_eliash."""
    md = T.one_orbital_square(Nx=8, Ny=8, Nw=32, temp=0.05, mu=0.2, tx=1.0, ty=0.5)
    U = 2.0
    lam_n, lam_s, _ = _trs_solvers(md['rvec'], md['ham_r'], 8, 8, 32, 0.05, 0.2, U, 0.0, gs=1)
    assert abs(lam_n - lam_s) < 1e-9, (lam_n, lam_s)
    with T.silence():
        res = P.calc_lin_eliash_nosym(8, 8, 1, 32, md['ham_r'], [], md['rvec'], 0.2, 0.05, md['olist'],
                                      md['site'], False, U, 0.0, 1, False, arnoldi_m=80, eps=1e-10, sw_out=False)
    if res['character']['even_freq'] > 0.99:   # leading solution is the even-frequency singlet
        assert abs(res['best_lambda'] - lam_n) < 1e-8, (res['best_lambda'], lam_n)
    else:                                        # an odd-frequency state leads: it must lie above
        assert res['best_lambda'] > lam_n


def test_soc_triplet_matches_nosoc_triplet_one_orbital():
    """One orbital (frequency parity exact): a p_x seed must stay in the odd-parity Sz=0 triplet
    sector of lin_eliash_soc (d_z seed + parity projection) and reproduce lin_eliash, instead of
    drifting to an even-parity odd-frequency triplet."""
    md = T.one_orbital_square(Nx=8, Ny=8, Nw=32, temp=0.05, mu=0.2, tx=1.0, ty=0.5)
    lam_n, lam_s, info = _trs_solvers(md['rvec'], md['ham_r'], 8, 8, 32, 0.05, 0.2, 2.0, 0.0, gs=-1)
    assert abs(lam_n - lam_s) < 1e-10, (lam_n, lam_s)
    ch = P.gap_character_nosym(info['delta'], info['kmap'], info['slist'], 8, 8, 1)
    assert ch['even_freq'] > 0.999 and ch['equal_spin'] < 1e-12
    d, kmap = info['delta'], info['kmap']
    flat = {tuple(k): i for i, k in enumerate(kmap)}
    km = np.array([flat[tuple((-k) % np.array([8, 8, 1]))] for k in kmap])
    assert np.linalg.norm(d + d[..., km]) < 1e-8*np.linalg.norm(d)        # odd parity


def test_soc_multiorbital_triplet_matches_exact_symmetry_sector():
    """Two orbitals: lambda of lin_eliash_soc (p_x seed) equals the leading eigenvalue of the
    exact kernel restricted to its true symmetry sector, Sz=0 triplet with odd parity.  The
    solution mixes in a few % of odd frequency, which lin_eliash (no SOC) cannot represent."""
    N, Nw, temp, mu, U, J = 4, 4, 0.1, 0.35, 1.0, 0.3
    rvec, ham_r = _two_orbital_dxz_dyz()
    _, lam_s, _ = _trs_solvers(rvec, ham_r, N, N, Nw, temp, mu, U, J, gs=-1, m=200)
    ham2, _, slist, invs = P.spin_double_model(ham_r, [])
    olist, site = P.get_chi_orb_list(2, [2])
    ol2, site2 = P.spin_double_chiolist(olist, site, 2)
    kmap, klist = _full_mesh(N, N)
    eig, uni = P.get_eigs(klist, ham2, [], rvec)
    Gk = F.gen_Green0(eig, uni, mu, temp, Nw)
    Vmat = F.gen_Vmatrix(ol2, slist, site2, invs, U, J)
    _, K = _bruteforce_chi0_kernel(Gk.transpose(3, 2, 1, 0), Vmat.T, ol2, temp, kmap, N, N, (0, 0))
    Nk, No = len(kmap), 4
    flat = {tuple(k): i for i, k in enumerate(kmap)}
    km = np.array([flat[tuple((-k) % np.array([N, N, 1]))] for k in kmap])

    def proj(v):                                    # Sz=0, spin-swap symmetric, odd parity
        D = v.reshape(Nk, Nw, No, No).copy()
        D[:, :, :2, :2] = 0
        D[:, :, 2:, 2:] = 0
        sw = D.copy()
        sw[:, :, :2, 2:], sw[:, :, 2:, :2] = D[:, :, 2:, :2], D[:, :, :2, 2:]
        D = 0.5*(D + sw)
        return (0.5*(D - D[km])).reshape(-1)
    Pm = np.array([proj(e) for e in np.eye(K.shape[0])]).T
    assert np.linalg.norm(Pm @ K - K @ Pm) < 1e-12*np.linalg.norm(K)      # exact symmetry sector
    u, sv, _ = np.linalg.svd(Pm)
    B = u[:, sv > 0.5]
    ev = np.linalg.eigvals(B.conj().T @ K @ B)
    assert abs(lam_s - ev.real.max()) < 1e-8, (lam_s, ev.real.max())


def test_multiorbital_triplet_is_exact_eigenpair_nosoc_and_soc():
    """Two orbitals: lin_eliash and lin_eliash_soc (p_x seed) must both return an exact eigenpair
    of the full spin-orbital kernel in the Sz=0, odd-parity triplet sector, odd-frequency admixture
    included, with a Hermitian (TRS-gauge) orbital gap.  Guards the TRIM parity projection of
    lin_eliash (without it the triplet keeps a spurious component at k=-k and a wrong lambda) and
    the real-inner-product Arnoldi of both solvers (their kernels fill components by conjugation,
    so they are only real-linear; a complex Arnoldi gave errors up to 1e-4 in lambda)."""
    N, Nw, temp, mu, U = 4, 4, 0.1, 0.35, 1.0
    olist, site = P.get_chi_orb_list(2, [2])
    ol2, site2 = P.spin_double_chiolist(olist, site, 2)
    klist, kmap, invk = F.gen_irr_k_TRS(N, N, 1)
    flat = {tuple(k): i for i, k in enumerate(kmap.tolist())}
    km = np.array([flat[tuple(((-k) % np.array([N, N, 1])).tolist())] for k in kmap])
    trim = np.all((2*kmap) % np.array([N, N, 1]) == 0, axis=1)
    for t4, J in ((0.1, 0.3), (0.1, 0.0)):
        rvec, ham_r = _two_orbital_dxz_dyz(t4=t4)
        eig, uni = P.get_eigs(klist, ham_r, [], rvec)
        Gk = F.gen_Green0(eig, uni, mu, temp, Nw)
        Smat, Cmat = F.gen_SCmatrix(olist, site, U, J)
        with T.silence():
            chi, _ = F.get_chi0(Smat, Cmat, Gk, olist, kmap, invk, temp, N, N, 1)
            d, lam = F.linearized_eliashberg(chi, Gk, uni, P.get_initial_gap(klist, 2, -1), Smat, Cmat, olist,
                                             np.ones(2), kmap, invk, N, N, 1, temp, -1, eps=1e-12, arnoldi_m=200)
        A = np.array([(1 if s == 0 else -1)*d[:, :, :, ik - 1].transpose(2, 1, 0) for ik, s, _ in invk])
        assert np.abs(A[trim]).max() < 1e-12*np.abs(A).max()                 # odd parity at TRIM
        Ad = np.conj(np.swapaxes(A, 2, 3))
        assert np.linalg.norm(A - np.vdot(Ad, A)/abs(np.vdot(Ad, A))*Ad) < 1e-12*np.linalg.norm(A)   # Hermitian
        ham2, _, slist, invs = P.spin_double_model(ham_r, [])
        eig2, uni2 = P.get_eigs(kmap/np.array([N, N, 1]), ham2, [], rvec)
        Gk2 = F.gen_Green0(eig2, uni2, mu, temp, Nw)
        Vmat = F.gen_Vmatrix(ol2, slist, site2, invs, U, J)
        _, K = _bruteforce_chi0_kernel(Gk2.transpose(3, 2, 1, 0), Vmat.T, ol2, temp, kmap, N, N, (0, 0))
        X = np.zeros((len(kmap), Nw, 4, 4), complex)
        X[:, :, :2, 2:] = A
        X[:, :, 2:, :2] = A                                                   # Sz=0 triplet, d_z
        x = X.reshape(-1)
        assert np.linalg.norm(K @ x - lam*x) < 1e-10*np.linalg.norm(K @ x), (t4, J, lam)
        _, lam_s, info = _trs_solvers(rvec, ham_r, N, N, Nw, temp, mu, U, J, gs=-1, eps=1e-12, m=200)
        xs = np.asarray(info['delta']).transpose(3, 2, 1, 0).reshape(-1)
        assert abs(lam_s - lam) < 1e-10, (t4, J, lam_s, lam)
        assert np.linalg.norm(K @ xs - lam_s*xs) < 1e-10*np.linalg.norm(K @ xs), (t4, J, lam_s)

        def proj(v):                                # Sz=0, spin-swap symmetric, odd parity
            D = v.reshape(len(kmap), Nw, 4, 4).copy()
            D[:, :, :2, :2] = 0
            D[:, :, 2:, 2:] = 0
            sw = D.copy()
            sw[:, :, :2, 2:], sw[:, :, 2:, :2] = D[:, :, 2:, :2], D[:, :, :2, 2:]
            D = 0.5*(D + sw)
            return (0.5*(D - D[km])).reshape(-1)
        u, sv, _ = np.linalg.svd(np.array([proj(e) for e in np.eye(K.shape[0])]).T)
        B = u[:, sv > 0.5]
        assert abs(lam - np.linalg.eigvals(B.conj().T @ K @ B).real.max()) < 1e-9, (t4, J, lam)


def test_soc_with_spin_orbit_is_exact_eigenpair_in_parity_sector():
    """Genuine SOC (L_z s_z on dxz/dyz): lin_eliash_soc must give the leading eigenpair of the exact
    spin-orbital kernel in its parity sector -- even (gap_sym=+1) / odd (-1) parity with Delta on the
    irreducible k-points and no singlet/triplet spin structure imposed (with SOC they mix), and the
    whole kernel for parity mixing (-10, full BZ).  Imposing D_ud = -/+ D_du gave 0.054 instead of
    0.186 for +1; the complex Arnoldi gave 0.42 for -10."""
    N, Nw, temp, mu, U, J, lso = 4, 4, 0.1, 0.35, 1.0, 0.3, 0.2
    rvec, ham_r = _two_orbital_dxz_dyz()
    ham2, _, slist, invs = P.spin_double_model(ham_r, [])
    r0 = [i for i, r in enumerate(rvec) if not r.any()][0]
    ham2 = ham2.copy()
    ham2[r0] += 0.5*lso*np.kron(np.diag([1, -1]), np.array([[0, -1j], [1j, 0]]))
    olist, site = P.get_chi_orb_list(2, [2])
    ol2, site2 = P.spin_double_chiolist(olist, site, 2)
    klist, kmap, invk = F.gen_irr_k_TRS(N, N, 1)
    eig, uni = P.get_eigs(klist, ham2, [], rvec)
    Gk = F.gen_Green0(eig, uni, mu, temp, Nw)
    Vmat = F.gen_Vmatrix(ol2, slist, site2, invs, U, J)
    ef, uf = P.get_eigs(kmap/np.array([N, N, 1]), ham2, [], rvec)
    Gf = F.gen_Green0(ef, uf, mu, temp, Nw)
    _, K = _bruteforce_chi0_kernel(Gf.transpose(3, 2, 1, 0), Vmat.T, ol2, temp, kmap, N, N, (0, 0))
    Nkall = len(kmap)
    flat = {tuple(k): i for i, k in enumerate(kmap.tolist())}
    km = np.array([flat[tuple(((-k) % np.array([N, N, 1])).tolist())] for k in kmap])
    for gs in (1, -1, -10):
        with T.silence():
            chis, sg, sg2, isc = F.get_chi0_soc(Vmat, Gk, ol2, slist, kmap, invk, invs, temp, N, N, 1)
            d, lam = F.linearized_eliashberg_soc(chis, Gk, uni, P.get_initial_gap(klist, 4, 1 if gs > 0 else -1),
                                                 Vmat, sg, sg2, np.ones(4), slist, ol2, kmap, invk, invs, isc,
                                                 N, N, 1, temp, gs, eps=1e-12, arnoldi_m=200)
        x = np.asarray(d).transpose(3, 2, 1, 0).reshape(-1)
        if gs == -10:
            B = np.eye(K.shape[0])
        else:                                       # parity sector, no spin structure
            p = 1 if gs > 0 else -1
            Pm = np.array([(0.5*(e.reshape(Nkall, Nw, 4, 4) + p*e.reshape(Nkall, Nw, 4, 4)[km])).reshape(-1)
                           for e in np.eye(K.shape[0])]).T
            u, sv, _ = np.linalg.svd(Pm)
            B = u[:, sv > 0.5]
            X = x.reshape(Nkall, Nw, 4, 4)
            assert np.linalg.norm(X - p*X[km]) < 1e-12*np.linalg.norm(X)
        assert abs(lam - np.linalg.eigvals(B.conj().T @ K @ B).real.max()) < 1e-9, (gs, lam)
        assert np.linalg.norm(K @ x - lam*x) < 1e-10*np.linalg.norm(K @ x), gs


def test_soc_band_gap_is_gauge_covariant_with_physical_parity():
    """conv_delta_orb_to_band_soc must not depend on the basis chosen inside Kramers-degenerate
    band pairs, for singlets and triplets; the band gap is even (singlet) / odd (triplet) in k."""
    rvec, ham_r = _two_orbital_dxz_dyz()
    rng = np.random.default_rng(1)
    for gs, parity in ((1, 1.0), (-1, -1.0)):
        _, _, info = _trs_solvers(rvec, ham_r, 8, 8, 32, 0.03, 0.35, 1.0, 0.3, gs=gs)
        uni = info['uni']
        uni_b = uni.copy()
        for k in range(len(uni)):
            for j in range(0, uni.shape[1], 2):               # rows of uni[k] are the eigenvectors
                q, _ = np.linalg.qr(rng.standard_normal((2, 2)) + 1j*rng.standard_normal((2, 2)))
                uni_b[k][j:j+2, :] = q.T @ uni[k][j:j+2, :]
        bA = F.conv_delta_orb_to_band_soc(info['delta'], uni, info['invk'], info['invs'], info['slist'], gs)
        bB = F.conv_delta_orb_to_band_soc(info['delta'], uni_b, info['invk'], info['invs'], info['slist'], gs)
        dA = np.array([bA[n, n] for n in range(len(bA))])
        dB = np.array([bB[n, n] for n in range(len(bB))])
        assert np.abs(np.abs(dA) - np.abs(dB)).max() < 1e-12 * max(1.0, np.abs(dA).max()) + 1e-14
        kmap = info['kmap']
        flat = {tuple(k): i for i, k in enumerate(kmap)}
        km = np.array([flat[tuple((-k) % np.array([8, 8, 1]))] for k in kmap])
        assert np.linalg.norm(dA - parity*dA[:, km]) < 1e-10*np.linalg.norm(dA)


def test_fflo_pair_momentum_obeys_gauge_shift():
    """H_A(k) = H(k+A) shifts both pair partners by A: lambda_{H_A}(Q) = lambda_H(Q+2A)."""
    Nx = Ny = 8
    Nw, temp, mu, U = 16, 0.05, -0.6, 2.0
    kmap, klist = _full_mesh(Nx, Ny)
    slist, invs, olist, site, Vmat = _spin_half(U)
    A = (1, 2)
    G = {}
    for lab, shift in (('H', (0.0, 0.0)), ('HA', (A[0]/Nx, A[1]/Ny))):
        eig, uni = F.get_eig(_rashba_zeeman_hk(klist, tp=-0.25, alpha=0.3, h=(0, 0, 0.1), shift=shift))
        G[lab] = F.gen_Green0(eig, uni, mu, temp, Nw)
    for Q in ((0, 0), (1, 1)):
        la = _nosym_lambda(G['HA'], olist, Vmat, kmap, Nx, Ny, temp, iQ=(*Q, 0), m=80)
        lb = _nosym_lambda(G['H'], olist, Vmat, kmap, Nx, Ny, temp, iQ=(Q[0] + 2*A[0], Q[1] + 2*A[1], 0), m=80)
        assert abs(la - lb) < 1e-9, (Q, la, lb)


def test_calc_lin_eliash_nosym_zeeman_spin_rotation_and_pauli_suppression(tmp_path):
    """Without SOC the Zeeman direction is irrelevant (hx = hz), and the field suppresses the
    singlet lambda; also exercises the plibs flow (spin doubling, Zeeman term, Q loop, output)."""
    md = T.one_orbital_square(Nx=16, Ny=16, Nw=64, temp=0.04, mu=-0.6, tx=1.0, ty=1.0)
    rvec = np.vstack([md['rvec'], [[1, 1, 0], [1, -1, 0], [-1, 1, 0], [-1, -1, 0]]])
    ham_r = np.ascontiguousarray(np.concatenate([md['ham_r'], np.full((4, 1, 1), 0.25, dtype=np.complex128)]))
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        lam = {}
        with T.silence():
            for lab, h in (('0', (0, 0, 0)), ('z', (0, 0, 0.1)), ('x', (0.1, 0, 0))):
                r = P.calc_lin_eliash_nosym(16, 16, 1, 64, ham_r, [], rvec, -0.6, 0.04, md['olist'], md['site'],
                                            False, 2.0, 0.0, 1, False, hvec=h, qlist=[[0, 0, 0], [1, 0, 0]],
                                            arnoldi_m=60, eps=1e-10, sw_out=(lab == 'z'))
                lam[lab] = r['best_lambda']
                assert r['character']['even_freq'] > 0.99
        assert os.path.exists('lambda_Q.dat') and os.path.exists('gap_nosym.npz')
        assert len(np.loadtxt('lambda_Q.dat')) == 2
    finally:
        os.chdir(cwd)
    assert abs(lam['x'] - lam['z']) < 1e-9, lam
    assert lam['z'] < lam['0'], lam


# --------------------------------------------------------------------------- #
#  standalone runner (no pytest required)
# --------------------------------------------------------------------------- #
if __name__ == '__main__':
    sys.exit(T.run_standalone(globals()))
