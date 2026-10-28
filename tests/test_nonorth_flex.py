"""Non-orthogonal Green functions and FLEX basis-normalization invariance."""
import importlib.util
import os
import sys
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import libs.flibs as F  # noqa: E402
import libs.plibs as P  # noqa: E402
import _tools as T  # noqa: E402


def test_nonorthogonal_green_matches_matrix_dyson():
    ham = np.array([[[1.0, 0.2+0.3j], [0.2-0.3j, -0.4]]], dtype=np.complex128)
    overlap = np.array([[[1.3, 0.1+0.2j], [0.1-0.2j, 0.9]]], dtype=np.complex128)
    eig, uni = F.get_eig(ham, overlap)
    temp, mu, nw = 0.1, 0.2, 4
    sigma = np.zeros((2, 2, nw, 1), dtype=np.complex128)
    free = F.gen_Green0(eig, uni, mu, temp, nw)
    assert np.allclose(F.gen_green(sigma, ham, mu, temp, overlap), free, atol=1e-12)
    sigma[0, 0, :, 0] = 0.2-0.04j
    sigma[1, 1, :, 0] = -0.1-0.03j
    sigma[0, 1, :, 0] = 0.07+0.02j
    sigma[1, 0, :, 0] = 0.04-0.01j
    actual = F.gen_green(sigma, ham, mu, temp, overlap)
    assert np.allclose(F.gen_green_from_eig(sigma, eig, uni, mu, temp, overlap), actual, atol=1e-12)
    from libs.plibs._calc import _green_at, _metric_green
    wn = (2*np.arange(nw)+1)*np.pi*temp
    assert np.allclose(_green_at(sigma, ham, mu, wn, overlap), actual, atol=1e-12)
    weighted = _metric_green(actual, overlap)
    for j in range(nw):
        z = mu+1j*(2*j+1)*np.pi*temp
        expected = np.linalg.inv(z*overlap[0]-ham[0]-sigma[:, :, j, 0])
        assert np.allclose(actual[:, :, j, 0], expected, atol=1e-12)
        assert np.allclose(np.trace(weighted[:, :, j, 0]), np.trace(overlap[0] @ expected), atol=1e-12)


def test_nonorthogonal_flex_normal_and_soc_match_scaled_orthogonal_model():
    n, nw, temp, fill, scale, u = 4, 32, 0.12, 0.4, 1.7, 1.2
    m = T.one_orbital_square(Nx=n, Ny=n, Nw=nw, temp=temp, mu=-0.4, U=u)
    ham = np.ascontiguousarray(scale*m['hamk'])
    overlap = np.full_like(ham, scale)
    eig, uni = F.get_eig(ham, overlap)
    smat, cmat = F.gen_SCmatrix(m['olist'], m['site'], U=u*scale**2, J=0.0)
    ham2_r, _, slist, invs = P.spin_double_model(scale*m['ham_r'], [])
    ol2, site2 = P.spin_double_chiolist(m['olist'], m['site'], 1)
    ham2 = F.gen_ham(m['klist'], ham2_r, m['rvec'])
    ov2 = np.ascontiguousarray(np.broadcast_to(scale*np.eye(2), ham2.shape), dtype=np.complex128)
    eig2, uni2 = F.get_eig(ham2, ov2)
    vmat = F.gen_Vmatrix(ol2, slist, site2, invs, u*scale**2, 0.0)
    with T.silence():
        ref, mu_ref = F.mkself(m['Smat'], m['Cmat'], m['kmap'], m['invk'], m['olist'], m['hamk'],
                               m['eig'], m['uni'], m['mu'], fill, temp, nw, n, n, 1, False, False,
                               sub_sigma=0, eps=1e-9, scf_loop=200)
        actual, mu_actual = F.mkself(smat, cmat, m['kmap'], m['invk'], m['olist'], ham, eig, uni,
                                     m['mu'], fill, temp, nw, n, n, 1, False, False, sub_sigma=0,
                                     eps=1e-9, scf_loop=200, overlap=overlap)
        soc, mu_soc = F.mkself_soc(vmat, m['kmap'], m['invk'], invs, ol2, slist, ham2, eig2, uni2,
                                   m['mu'], 2*fill, temp, nw, n, n, 1, False, False, sub_sigma=0,
                                   eps=1e-9, scf_loop=200, overlap=ov2)
    assert abs(mu_actual-mu_ref) < 1e-8
    assert abs(mu_soc-mu_ref) < 1e-8
    assert T.rel_err(actual, scale*ref) < 1e-7
    assert T.rel_err(soc[:1, :1], actual) < 1e-7
    with T.silence():
        gref = F.gen_green(ref, m['hamk'], mu_ref, temp)
        g = F.gen_green(actual, ham, mu_actual, temp, overlap)
        gsoc = F.gen_green(soc, ham2, mu_soc, temp, ov2)
        chi_ref, _ = F.get_chi0(m['Smat'], m['Cmat'], gref, m['olist'], m['kmap'], m['invk'], temp, n, n, 1)
        chi, _ = F.get_chi0(smat, cmat, g, m['olist'], m['kmap'], m['invk'], temp, n, n, 1)
        init = P.get_initial_gap(m['klist'], 1, 1)
        dref, lam_ref = F.linearized_eliashberg(chi_ref, gref, m['uni'], init, m['Smat'], m['Cmat'],
                                                m['olist'], m['plist'], m['kmap'], m['invk'], n, n, 1,
                                                temp, 1, eps=1e-9)
        _, lam = F.linearized_eliashberg(chi, g, uni, init, smat, cmat, m['olist'], m['plist'],
                                         m['kmap'], m['invk'], n, n, 1, temp, 1, eps=1e-9, overlap=overlap)
        chis, sg, sg2, isc = F.get_chi0_soc(vmat, gsoc, ol2, slist, m['kmap'], m['invk'], invs,
                                            temp, n, n, 1)
        _, lam_soc = F.linearized_eliashberg_soc(chis, gsoc, uni2, P.get_initial_gap(m['klist'], 2, 1),
                                                 vmat, sg, sg2, np.ones(2), slist, ol2, m['kmap'], m['invk'],
                                                 invs, isc, n, n, 1, temp, 1, eps=1e-9, overlap=ov2)
        seed = np.ascontiguousarray(dref*0.01)
        d1, s1 = F.nonlinear_eliashberg(seed.copy(), gref.copy(), m['hamk'], m['Smat'], m['Cmat'],
                                        m['olist'], m['plist'], m['kmap'], m['invk'], mu_ref, temp, 1, n, n, 1,
                                        sw_sigma=True, sw_Vconst=True, itemax=2, m_diis=1,
                                        sw_amp_newton=False, gap_min=0.0, sub_sigma=0)
        d2, s2 = F.nonlinear_eliashberg(scale*seed.copy(), g.copy(), ham, smat, cmat, m['olist'], m['plist'],
                                        m['kmap'], m['invk'], mu_actual, temp, 1, n, n, 1,
                                        sw_sigma=True, sw_Vconst=True, itemax=2, m_diis=1,
                                        sw_amp_newton=False, gap_min=0.0, sub_sigma=0, overlap=overlap)
    assert abs(lam-lam_ref) < 1e-8
    assert abs(lam_soc-lam_ref) < 1e-8
    assert T.rel_err(d2, scale*d1) < 1e-7
    assert T.rel_err(s2, scale*s1) < 1e-7


def test_nonorthogonal_multiorbital_flex_soc_matches_nosoc():
    n, temp, fill = 4, 0.1, 0.8
    m = T.two_orbital_square(Nx=n, Ny=n, Nw=16, temp=temp, mu=-0.2, U=0.5, J=0.05)
    sr = np.zeros_like(m['ham_r'])
    sr[0] = [[1.2, 0.15], [0.15, 0.9]]
    sr[1] = sr[2] = [[0.08, 0.02], [0.02, -0.03]]
    ham, overlap = F.gen_ham(m['klist'], m['ham_r'], m['rvec'], Ovl_r=sr)
    eig, uni = F.get_eig(ham, overlap)
    ham2_r, sr2, slist, invs = P.spin_double_model(m['ham_r'], sr)
    ham2, overlap2 = F.gen_ham(m['klist'], ham2_r, m['rvec'], Ovl_r=sr2)
    eig2, uni2 = F.get_eig(ham2, overlap2)
    ol2, site2 = P.spin_double_chiolist(m['olist'], m['site'], 2)
    vmat = F.gen_Vmatrix(ol2, slist, site2, invs, 0.5, 0.05)
    grids = [None]
    if importlib.util.find_spec('sparse_ir') is not None:
        grids.append(P.ir_grid(temp, 10.0, eps=1e-8))
    for grid in grids:
        nw = 16 if grid is None else grid['Nw']
        try:
            P.set_freq_grid(grid)
            with T.silence():
                normal, mu1 = F.mkself(m['Smat'], m['Cmat'], m['kmap'], m['invk'], m['olist'], ham, eig, uni,
                                       m['mu'], fill, temp, nw, n, n, 1, False, False, sub_sigma=0,
                                       eps=1e-9, scf_loop=200, overlap=overlap)
                soc, mu2 = F.mkself_soc(vmat, m['kmap'], m['invk'], invs, ol2, slist, ham2, eig2, uni2,
                                        m['mu'], 2*fill, temp, nw, n, n, 1, False, False, sub_sigma=0,
                                        eps=1e-9, scf_loop=200, overlap=overlap2)
            assert abs(mu1-mu2) < 1e-8
            assert T.rel_err(normal, soc[:2, :2]) < 1e-7
        finally:
            P.set_freq_grid(None)


def test_nonorthogonal_reload_uses_overlap_for_both_spin_paths(tmp_path):
    from libs.plibs import _calc
    cwd = os.getcwd()
    os.chdir(tmp_path)   # _save_sigma / the reload use the working directory
    try:
        n, nw, temp, mu = 4, 8, 0.1, -0.2
        m = T.one_orbital_square(Nx=n, Ny=n, Nw=nw, temp=temp, mu=mu)
        sr = np.zeros_like(m['ham_r'])
        sr[0, 0, 0] = 1.4
        state = _calc._prepare_kspace_state(n, n, 1, m['ham_r'], sr, m['rvec'], need_ham=True)
        sigma = np.full((1, 1, nw, len(m['klist'])), 0.2-0.1j, dtype=np.complex128)
        _calc._save_sigma(sigma, mu, temp)
        with T.silence():
            result = _calc._prepare_green_state_normal(state, m['olist'], (m['Smat'], m['Cmat']), mu, 0.5,
                                                       temp, nw, n, n, 1, True, sw_from_file=True)
        assert np.allclose(result['Gk'], F.gen_green(sigma, state['ham_k'], mu, temp, state['overlap_k']))
        ham2_r, sr2, slist, invs = P.spin_double_model(m['ham_r'], sr)
        state2 = _calc._prepare_kspace_state(n, n, 1, ham2_r, sr2, m['rvec'], need_ham=True)
        ol2, site2 = P.spin_double_chiolist(m['olist'], m['site'], 1)
        vmat = F.gen_Vmatrix(ol2, slist, site2, invs, 0.0, 0.0)
        sigma2 = np.zeros((2, 2, nw, len(m['klist'])), dtype=np.complex128)
        sigma2[0, 0] = sigma2[1, 1] = sigma[0, 0]
        _calc._save_sigma(sigma2, mu, temp)
        with T.silence():
            result2 = _calc._prepare_green_state_soc(state2, ol2, slist, invs, vmat, mu, 1.0, temp, nw,
                                                     n, n, 1, True, sw_from_file=True)
        assert np.allclose(result2['Gk'], F.gen_green(sigma2, state2['ham_k'], mu, temp, state2['overlap_k']))
    finally:
        os.chdir(cwd)


def test_nonorthogonal_flex_with_soc_matches_scaled_orthogonal_model():
    from test_rpa_flex import _soc_model
    n, nw, temp, scale, u, j = 4, 12, 0.1, 1.5, 0.4, 0.05
    m = _soc_model(Nx=n, Ny=n, Nw=nw, temp=temp, U=u, J=j)
    ham = F.gen_ham(m['klist'], m['ham_r'], m['rvec'])
    ham_scaled = np.ascontiguousarray(scale*ham)
    overlap = np.ascontiguousarray(np.broadcast_to(scale*np.eye(m['Norb']), ham.shape), dtype=np.complex128)
    eig, uni = F.get_eig(ham_scaled, overlap)
    vmat = F.gen_Vmatrix(m['olist'], m['slist'], m['site'], m['invs'], u*scale**2, j*scale**2)
    common = (m['kmap'], m['invk'], m['invs'], m['olist'], m['slist'])
    with T.silence():
        ref, mu_ref = F.mkself_soc(m['Vmat'], *common, ham, m['eig'], m['uni'], 0.0, 2.0, temp, nw,
                                   n, n, 1, False, False, sub_sigma=0, eps=1e-9, scf_loop=200)
        actual, mu_actual = F.mkself_soc(vmat, *common, ham_scaled, eig, uni, 0.0, 2.0, temp, nw,
                                         n, n, 1, False, False, sub_sigma=0, eps=1e-9, scf_loop=200,
                                         overlap=overlap)
    assert abs(mu_actual-mu_ref) < 1e-8
    assert T.rel_err(actual, scale*ref) < 1e-7


def test_nonorthogonal_chi0_tail_correction_is_second_order():
    """chi0_ref_tau's tau=0 jump: G0 jumps by S(k)^-1, not 1, so J(q) is a k-convolution of S^-1 and rho.
    With the identity in its place the tail-corrected chi0 stays O(1/Nw) on an MLO basis."""
    n, temp, mu = 8, 0.1, -0.2
    m = T.two_orbital_square(Nx=n, Ny=n, Nw=16, temp=temp, mu=mu, U=0.5, J=0.05)
    sr = np.zeros_like(m['ham_r'])
    sr[0] = [[1.2, 0.15], [0.15, 0.9]]
    sr[1] = sr[2] = [[0.08, 0.02], [0.02, -0.03]]
    ham, overlap = F.gen_ham(m['klist'], m['ham_r'], m['rvec'], Ovl_r=sr)
    eig, uni = F.get_eig(ham, overlap)

    def chi(Nw, tail):
        G = F.gen_Green0(eig, uni, mu, temp, Nw)
        with T.silence():
            if tail:
                c = F.get_chi0_tail(m['Smat'], m['Cmat'], G, eig, uni, m['olist'], m['kmap'], m['invk'],
                                    mu, temp, n, n, 1)[0]
            else:
                c = F.get_chi0(m['Smat'], m['Cmat'], G, m['olist'], m['kmap'], m['invk'], temp, n, n, 1)[0]
        return c[:, :, :4]
    ref = 2*chi(8192, False)-chi(4096, False)    # Richardson: plain convolution is O(1/Nw)
    e64, e128 = T.rel_err(chi(64, True), ref), T.rel_err(chi(128, True), ref)
    assert e128 < 0.3*e64, (e64, e128)          # second order: ~1/4 per doubling (first order: 1/2)
    assert e128 < 1e-3


def test_mlo_orbital_normalization_keeps_bands_and_fixes_the_local_norm():
    """normalize_mlo_orbitals: S_aa(R=0) -> 1 with the bands unchanged (ecalj MLOs are not normalized)."""
    m = T.two_orbital_square(Nx=4, Ny=4, Nw=8, temp=0.1, mu=-0.2)
    sr = np.zeros_like(m['ham_r'])
    sr[0] = [[0.02, 0.003], [0.003, 0.5]]
    sr[1] = sr[2] = [[0.001, 0.0005], [0.0005, -0.02]]
    hr = m['ham_r']*np.sqrt(np.outer(sr[0].diagonal().real, sr[0].diagonal().real))
    h1, s1 = P.normalize_mlo_orbitals(m['rvec'], hr, sr)
    i0 = np.flatnonzero((np.abs(m['rvec']) < 1e-9).all(1))[0]
    assert np.allclose(s1[i0].diagonal(), 1.0)
    e0 = F.get_eig(*F.gen_ham(m['klist'], hr, m['rvec'], Ovl_r=sr), sw=False)
    e1 = F.get_eig(*F.gen_ham(m['klist'], h1, m['rvec'], Ovl_r=s1), sw=False)
    assert np.allclose(e0, e1, atol=1e-10)


def test_nonorthogonal_flex_is_covariant_under_an_orbital_rescaling():
    """Rescaling the orbitals by W = diag(w) (H, S -> W^-1 H W^-1, W^-1 S W^-1, the vertex by the four w's of its
    orbital pair indices) is the same physics: mu and the Stoner factor are unchanged and Sigma -> W^-1 Sigma W^-1.
    A FLEX that dropped S from G^-1 (the pre-MLO code) breaks this, unlike the scalar-overlap checks above."""
    n, temp, fill, nw = 4, 0.1, 0.8, 32
    m = T.two_orbital_square(Nx=n, Ny=n, Nw=16, temp=temp, mu=-0.2, U=1.0, J=0.1)
    sr = np.zeros_like(m['ham_r'])
    sr[0] = [[1.2, 0.15], [0.15, 0.9]]
    sr[1] = sr[2] = [[0.08, 0.02], [0.02, -0.03]]
    w = np.array([0.4, 1.6])
    fac = np.array([w[a-1]*w[b-1] for a, b in m['olist']])
    models = {'norm': (m['ham_r'], sr, m['Smat'], m['Cmat']),
              'raw': (m['ham_r']/np.outer(w, w), sr/np.outer(w, w), m['Smat']/np.outer(fac, fac),
                      m['Cmat']/np.outer(fac, fac))}
    res = {}
    for tag, (hr, s_r, smat, cmat) in models.items():
        ham, overlap = F.gen_ham(m['klist'], hr, m['rvec'], Ovl_r=s_r)
        eig, uni = F.get_eig(ham, overlap)
        with T.silence():
            sigma, mu = F.mkself(smat, cmat, m['kmap'], m['invk'], m['olist'], ham, eig, uni, m['mu'], fill, temp,
                                 nw, n, n, 1, False, False, sub_sigma=1, eps=1e-10, scf_loop=300, overlap=overlap)
            g = F.gen_green(sigma, ham, mu, temp, overlap=overlap)
            _, stoner = F.get_chi0(smat, cmat, g, m['olist'], m['kmap'], m['invk'], temp, n, n, 1)
        res[tag] = (sigma, mu, stoner)
    (s1, mu1, st1), (s2, mu2, st2) = res['norm'], res['raw']
    assert abs(mu1-mu2) < 1e-8 and abs(st1-st2) < 1e-8
    assert T.rel_err(s2, s1/np.outer(w, w)[:, :, None, None]) < 1e-7


if __name__ == '__main__':
    sys.exit(T.run_standalone(globals()))
