"""
IR (intermediate representation) frequency grid of the convolution kernels (plibs.ir_grid,
plibs.set_freq_grid).

  * the grid matrices: a fermionic function is reproduced from its sampled values, and the tau
    points are symmetric about beta/2 (the reflection of the bubble);
  * chi0 on the IR grid equals the exact Lindhard bubble at the sampled bosonic frequencies (no
    Matsubara cutoff error, unlike the uniform grid);
  * FLEX (mkself) on the IR grid gives the mu and Sigma that the uniform grid tends to as Nw grows.
The uniform grid is restored after each test.  Skipped without sparse-ir.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import libs.plibs as P  # noqa: E402
import libs.flibs as F  # noqa: E402
import _tools as T  # noqa: E402

try:
    import sparse_ir  # noqa: F401
    HAVE_IR = True
except ImportError:
    HAVE_IR = False


def test_ir_grid_reproduces_a_fermionic_function_and_reflects_tau():
    if not HAVE_IR:
        T.skip('sparse-ir not installed')
    temp = 0.05
    g = P.ir_grid(temp, 10.0, eps=1e-8)
    beta = 1.0/temp
    assert np.abs(g['tau'][g['rev']-1]-(beta-g['tau'])).max() < 1e-9*beta
    # G(iw) = 1/(iw - e) of poles inside the window, through its IR coefficients
    n = 2*np.arange(400)+1
    for e in (-3.0, 0.1, 2.5):
        x = 1.0/(1j*g['wF']-e)
        xl = P.ir_fit_fermion(g, x)
        assert T.rel_err(P.ir_eval_fermion(g, xl, n), 1.0/(1j*n*np.pi*temp-e)) < 1e-6, e
        # density weights: T sum_n G(iw_n) e^{iw 0+}, midpoint of tau = 0+-: f(e) - 1/2
        dens = (g['dens'][:g['Nw']]*x+g['dens'][g['Nw']:]*np.conj(x)).sum().real
        assert abs(dens-(T.fermi(e, temp)-0.5)) < 1e-6, e


def test_ir_chi0_is_the_exact_bubble_at_the_sampled_frequencies():
    if not HAVE_IR:
        T.skip('sparse-ir not installed')
    temp, N = 0.05, 8
    st = T.one_orbital_square(Nx=N, Ny=N, Nw=16, temp=temp, mu=-0.4)
    g = P.ir_grid(temp, 10.0, eps=1e-8)
    try:
        P.set_freq_grid(g)
        Gk = F.gen_Green0(st['eig'], st['uni'], st['mu'], temp, g['Nw'])
        with T.silence():
            chi, _ = F.get_chi0(st['Smat'], st['Cmat'], Gk, st['olist'], st['kmap'], st['invk'], temp, N, N, 1)
    finally:
        P.set_freq_grid(None)
    m = np.rint(g['wB']/(2*np.pi*temp)).astype(int)
    eps_full = T.full_grid_dispersion(st['ham_r'], st['rvec'], N, N)[..., 0]
    ex = T.exact_chi0_1orb(eps_full, st['klist'], m, temp, st['mu'])
    assert T.rel_err(chi[0, 0].T, ex) < 1e-6
    # the uniform grid at the same Nw is off by its O(1/Nw) cutoff error
    with T.silence():
        chiu, _ = F.get_chi0(st['Smat'], st['Cmat'], st['Gk'], st['olist'], st['kmap'], st['invk'], temp, N, N, 1)
    assert T.rel_err(chiu[0, 0, :1].T, ex[:, :1]) > 1e-3


def test_ir_ignores_the_chi0_tail_correction():
    if not HAVE_IR:
        T.skip('sparse-ir not installed')
    temp, N = 0.05, 8
    st = T.two_orbital_square(Nx=N, Ny=N, Nw=16, temp=temp, mu=-0.5, U=1.2, J=0.15)
    g = P.ir_grid(temp, 10.0, eps=1e-8)
    try:
        P.set_freq_grid(g)
        Gk = F.gen_Green0(st['eig'], st['uni'], st['mu'], temp, g['Nw'])
        with T.silence():
            chi, _ = F.get_chi0(st['Smat'], st['Cmat'], Gk, st['olist'], st['kmap'], st['invk'], temp, N, N, 1)
            chit, _ = F.get_chi0_tail(st['Smat'], st['Cmat'], Gk, st['eig'], st['uni'], st['olist'], st['kmap'],
                                      st['invk'], st['mu'], temp, N, N, 1)
            s = [F.mkself(st['Smat'], st['Cmat'], st['kmap'], st['invk'], st['olist'], st['hamk'], st['eig'],
                          st['uni'], st['mu'], 1.0, temp, g['Nw'], N, N, 1, False, False, sub_sigma=0, eps=1e-9,
                          scf_loop=200, sw_tail=tail) for tail in (False, True)]
    finally:
        P.set_freq_grid(None)
    assert np.array_equal(chit, chi)
    # two separate FLEX runs: equal up to the run-to-run roundoff of the threaded sums (gfortran + MKL)
    assert T.rel_err(s[1][0], s[0][0]) < 1e-10 and abs(s[1][1]-s[0][1]) < 1e-10


def test_ir_spectrum_is_continued_from_the_sampled_self_energy(tmp_path):
    if not HAVE_IR:
        T.skip('sparse-ir not installed')
    temp, N, mu, delta = 0.05, 8, -0.5, 0.02
    st = T.two_orbital_square(Nx=N, Ny=N, Nw=16, temp=temp, mu=mu)
    Norb = st['hamk'].shape[1]
    # local Sigma with two poles: the exact real-axis spectrum is known
    sig = lambda z: (0.15/(z-0.8)+0.15/(z+0.8))[..., None, None]*np.eye(Norb)
    g = P.ir_grid(temp, 10.0, eps=1e-8)
    sigmak = np.ascontiguousarray(np.broadcast_to(np.moveaxis(sig(1j*g['wF']), 0, 2)[..., None],
                                                  (Norb, Norb, g['Nw'], len(st['hamk']))))
    cwd = os.getcwd()
    os.chdir(tmp_path)   # the drivers write and read in the working directory
    try:
        try:
            P.set_freq_grid(g)
            with T.silence():
                P.output_self_wannier(sigmak, mu, st['kmap'], st['invk'], N, N, 1, g['Nw'], temp)
        finally:
            P.set_freq_grid(None)
        d = np.load('self_en_wannier.npz')
        assert np.allclose(d['iw_ir'], g['wF'])
        # the uniform-grid part keeps its format: Sigma at (2n+1) pi T, R = 0
        r0 = np.flatnonzero((d['rvec'] == 0).all(1))[0]
        assert T.rel_err(np.moveaxis(d['sigma'][:, :, :, r0], 2, 0), sig(1j*d['iw'])) < 1e-8
        with T.silence():
            sp, _, wl, A = P.spectrum_on_path([[0, 0, 0], [0.5, 0, 0]], 8, np.eye(3), mu, st['ham_r'], [], st['rvec'],
                                              -1.5, 1.5, delta, 301, temp)
        kl = np.linspace(0, 0.5, 9)[:, None]*[1, 0, 0]
        hk = F.gen_ham(kl, st['ham_r'], st['rvec'])
        z = wl+1j*delta
        Aex = -np.trace(np.linalg.inv((z+mu)[None, :, None, None]*np.eye(Norb)-hk[:, None]-sig(z)[None]),
                        axis1=2, axis2=3).imag
        # Pade on the sampled Sigma (a two-pole one is rational: exact); on its uniform-grid evaluation the
        # basis truncation would be continued into ~50 % spurious structure
        assert A.shape == Aex.shape
        assert np.abs(A-Aex).max() < 1e-4*Aex.max()
    finally:
        os.chdir(cwd)


def test_ir_flex_is_the_large_nw_limit_of_the_uniform_grid():
    if not HAVE_IR:
        T.skip('sparse-ir not installed')
    temp, N = 0.05, 8
    st = T.two_orbital_square(Nx=N, Ny=N, Nw=16, temp=temp, mu=-0.5, U=1.2, J=0.15)
    g = P.ir_grid(temp, 10.0, eps=1e-8)

    def run(Nw):
        with T.silence():
            return F.mkself(st['Smat'], st['Cmat'], st['kmap'], st['invk'], st['olist'], st['hamk'], st['eig'],
                            st['uni'], st['mu'], 1.0, temp, Nw, N, N, 1, False, False, sub_sigma=0, eps=1e-9,
                            scf_loop=200)
    try:
        P.set_freq_grid(g)
        sI, muI = run(g['Nw'])
    finally:
        P.set_freq_grid(None)
    # uniform grid: O(1/Nw) + O(1/Nw^2), removed by Richardson extrapolation over Nw = 1024, 2048, 4096
    n = np.rint((g['wF']/(np.pi*temp)-1)/2).astype(int)        # uniform index of each sampled w
    res = {Nw: run(Nw) for Nw in (1024, 2048, 4096)}
    ext = (8*res[4096][0][:, :, n]-6*res[2048][0][:, :, n]+res[1024][0][:, :, n])/3
    muext = (8*res[4096][1]-6*res[2048][1]+res[1024][1])/3
    assert abs(muI-muext) < 2e-6, (muI, muext)
    assert T.rel_err(sI, ext) < 2e-3
    assert T.rel_err(sI, ext) < 0.5*T.rel_err(res[4096][0][:, :, n], ext)


def test_ir_linearized_eliashberg_eigenvalue_is_the_large_nw_limit():
    if not HAVE_IR:
        T.skip('sparse-ir not installed')
    temp, N = 0.05, 8
    st = T.two_orbital_square(Nx=N, Ny=N, Nw=16, temp=temp, mu=-0.5, U=1.2, J=0.15)
    g = P.ir_grid(temp, 10.0, eps=1e-8)
    with T.silence():
        init = P.get_initial_gap(st['klist'], 2, 1)

    def lam(Nw):
        G = F.gen_Green0(st['eig'], st['uni'], st['mu'], temp, Nw)
        with T.silence():
            chi, _ = F.get_chi0(st['Smat'], st['Cmat'], G, st['olist'], st['kmap'], st['invk'], temp, N, N, 1)
            return F.linearized_eliashberg(chi, G, st['uni'], init, st['Smat'], st['Cmat'], st['olist'], st['plist'],
                                           st['kmap'], st['invk'], N, N, 1, temp, 1, eps=1e-10)[1]
    try:
        P.set_freq_grid(g)
        lI = lam(g['Nw'])
    finally:
        P.set_freq_grid(None)
    l512, l1024 = lam(512), lam(1024)
    ext = 2*l1024-l512                                       # the uniform grid is O(1/Nw)
    assert abs(lI-ext) < 1e-6*abs(ext), (lI, ext)
    assert abs(lI-ext) < 0.1*abs(l1024-ext)


def test_ir_soc_kernels_reproduce_the_no_soc_ones_on_a_spin_doubled_model():
    if not HAVE_IR:
        T.skip('sparse-ir not installed')
    temp, N, No = 0.05, 8, 2
    st = T.two_orbital_square(Nx=N, Ny=N, Nw=16, temp=temp, mu=-0.5, U=1.2, J=0.15)
    ham2 = np.zeros((len(st['rvec']), 2*No, 2*No), dtype=complex)
    ham2[:, :No, :No] = ham2[:, No:, No:] = st['ham_r']
    hk2 = F.gen_ham(st['klist'], ham2, st['rvec'])
    eig2, uni2 = F.get_eig(hk2)
    slist = np.array([1]*No+[-1]*No, dtype=np.int64)
    invs = np.concatenate([np.arange(No, 2*No), np.arange(No)]).astype(np.int64)+1
    olist2, site2 = P.get_chi_orb_list(2*No, [2*No])
    Vmat = F.gen_Vmatrix(olist2, slist, site2, invs, U=1.2, J=0.15)
    g = P.ir_grid(temp, 10.0, eps=1e-8)
    try:
        P.set_freq_grid(g)
        with T.silence():
            s1, mu1 = F.mkself(st['Smat'], st['Cmat'], st['kmap'], st['invk'], st['olist'], st['hamk'], st['eig'],
                               st['uni'], st['mu'], 1.0, temp, g['Nw'], N, N, 1, False, False, eps=1e-10,
                               scf_loop=300)
            s2, mu2 = F.mkself_soc(Vmat, st['kmap'], st['invk'], invs, olist2, slist, hk2, eig2, uni2, st['mu'], 2.0,
                                   temp, g['Nw'], N, N, 1, False, False, eps=1e-10, scf_loop=300)
            G1 = F.gen_Green0(st['eig'], st['uni'], st['mu'], temp, g['Nw'])
            G2 = F.gen_Green0(eig2, uni2, st['mu'], temp, g['Nw'])
            chi1, _ = F.get_chi0(st['Smat'], st['Cmat'], G1, st['olist'], st['kmap'], st['invk'], temp, N, N, 1)
            l1 = F.linearized_eliashberg(chi1, G1, st['uni'], P.get_initial_gap(st['klist'], No, 1), st['Smat'],
                                         st['Cmat'], st['olist'], st['plist'], st['kmap'], st['invk'], N, N, 1, temp,
                                         1, eps=1e-8)[1]
            chis, sg, sg2, isc = F.get_chi0_soc(Vmat, G2, olist2, slist, st['kmap'], st['invk'], invs, temp, N, N, 1)
            l2 = F.linearized_eliashberg_soc(chis, G2, uni2, P.get_initial_gap(st['klist'], 2*No, 1), Vmat, sg, sg2,
                                             np.ones(2*No), slist, olist2, st['kmap'], st['invk'], invs, isc, N, N, 1,
                                             temp, 1, eps=1e-8)[1]
    finally:
        P.set_freq_grid(None)
    assert abs(mu1-mu2) < 1e-9
    assert T.rel_err(s2[:No, :No], s1) < 1e-8 and T.rel_err(s2[No:, No:], s1) < 1e-8
    assert np.abs(s2[:No, No:]).max() < 1e-10
    assert abs(l1-l2) < 1e-8*abs(l1)


def test_ir_phonon_vertex_gives_the_isotropic_gap_equation_eigenvalue():
    """vph built on boson_frequencies of the IR grid: for one band at U = 0 the eigenvalue is that of
    lambda Delta_n = -T sum_m W(w_n - w_m) <1/(w_m^2 + xi^2)>_k Delta_m, here solved on a wide window."""
    if not HAVE_IR:
        T.skip('sparse-ir not installed')
    Nx = Ny = 16
    temp, mu, w0, g2 = 0.02, -0.1, 0.03, 0.004
    rvec = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]], dtype=np.float64)
    ham_r = np.zeros((5, 1, 1), dtype=np.complex128)
    ham_r[1:, 0, 0] = -0.25
    klist, kmap, invk = F.gen_irr_k_TRS(Nx, Ny, 1)
    eig, uni = F.get_eig(F.gen_ham(klist, ham_r, rvec))
    olist, site = P.get_chi_orb_list(1, [1])
    Smat = Cmat = np.zeros((1, 1))
    W = lambda v: -g2*2*w0/(v**2+w0**2)
    g = P.ir_grid(temp, 3.0, eps=1e-8)
    try:
        P.set_freq_grid(g)
        vj = P.boson_frequencies(temp, g['Nw'])
        assert np.allclose(vj, g['wB'])
        Gk = F.gen_Green0(eig, uni, mu=mu, temp=temp, Nw=g['Nw'])
        with T.silence():
            chi, _ = F.get_chi0(Smat, Cmat, Gk, olist, kmap, invk, temp, Nx, Ny, 1)
            _, lam = F.linearized_eliashberg(chi, Gk, uni, P.get_initial_gap(klist, 1, 0), Smat, Cmat, olist,
                                             F.get_plist(rvec, ham_r), kmap, invk, Nx, Ny, 1, temp, gap_sym=0,
                                             eps=1e-10, vph=W(vj)[None, None, :])
    finally:
        P.set_freq_grid(None)
    assert np.allclose(P.boson_frequencies(temp, 8), 2*np.pi*temp*np.arange(8))
    _, klf = P.gen_klist(Nx, Ny, 1, sw_pp=False)
    xi = P.get_eigs(klf, ham_r, [], rvec)[0][:, 0]-mu
    M = 1500                                           # |w_n| < 2 pi T M ~ 190 eV
    wn = (2*np.arange(-M, M)+1)*np.pi*temp
    chim = (1/(wn[None, :]**2+xi[:, None]**2)).mean(axis=0)
    ref = np.linalg.eigvals(-temp*W(wn[:, None]-wn[None, :])*chim[None, :]).real.max()
    assert lam > 0 and abs(lam-ref) < 1e-8*ref, (lam, ref)   # the uniform v_j on the IR slots: 1e-6


def test_ir_nosym_eliashberg_matches_the_symmetric_solver_and_the_large_nw_limit():
    """nosym kernels on IR: (a) a spin-doubled no-SOC model gives the no-SOC leading eigenvalue (singlet here);
    (b) Rashba + Zeeman with a finite pair momentum Q gives the large-Nw limit of the uniform grid."""
    if not HAVE_IR:
        T.skip('sparse-ir not installed')
    from test_eliash_nosym import _rashba_zeeman_hk, _spin_half, _nosym_lambda
    temp, N, No = 0.05, 8, 2
    st = T.two_orbital_square(Nx=N, Ny=N, Nw=16, temp=temp, mu=-0.5, U=1.2, J=0.15)
    _, kl = P.gen_klist(N, N, 1, sw_pp=False)
    kmap = np.ascontiguousarray(np.round(kl*np.array([N, N, 1])).astype(np.int64))
    ham2 = np.zeros((len(st['rvec']), 2*No, 2*No), dtype=complex)
    ham2[:, :No, :No] = ham2[:, No:, No:] = st['ham_r']
    e2, u2 = F.get_eig(F.gen_ham(kl, ham2, st['rvec']))
    sl2 = np.array([1]*No+[-1]*No, dtype=np.int64)
    iv2 = np.concatenate([np.arange(No, 2*No), np.arange(No)]).astype(np.int64)+1
    ol2, si2 = P.get_chi_orb_list(2*No, [2*No])
    V2 = F.gen_Vmatrix(ol2, sl2, si2, iv2, U=1.2, J=0.15)
    g = P.ir_grid(temp, 10.0, eps=1e-8)
    try:
        P.set_freq_grid(g)
        lN = _nosym_lambda(F.gen_Green0(e2, u2, st['mu'], temp, g['Nw']), ol2, V2, kmap, N, N, temp, sl2)
        G1 = F.gen_Green0(st['eig'], st['uni'], st['mu'], temp, g['Nw'])
        with T.silence():
            chi1, _ = F.get_chi0(st['Smat'], st['Cmat'], G1, st['olist'], st['kmap'], st['invk'], temp, N, N, 1)
            lS = F.linearized_eliashberg(chi1, G1, st['uni'], P.get_initial_gap(st['klist'], No, 0), st['Smat'],
                                         st['Cmat'], st['olist'], st['plist'], st['kmap'], st['invk'], N, N, 1, temp,
                                         0, eps=1e-10)[1]
    finally:
        P.set_freq_grid(None)
    assert abs(lN-lS) < 1e-8*lS, (lN, lS)
    # (b) Rashba + Zeeman, Q = (1/N, 0, 0)
    temp = 0.1
    e, u = F.get_eig(_rashba_zeeman_hk(kl))
    slist, invs, olist, site, Vmat = _spin_half(1.5)
    g = P.ir_grid(temp, 10.0, eps=1e-8)

    def lam(Nw):
        return _nosym_lambda(F.gen_Green0(e, u, 0.2, temp, Nw), olist, Vmat, kmap, N, N, temp, slist, iQ=(1, 0, 0))
    try:
        P.set_freq_grid(g)
        lI = lam(g['Nw'])
    finally:
        P.set_freq_grid(None)
    l1, l2, l3 = lam(128), lam(256), lam(512)
    ext = (8*l3-6*l2+l1)/3                             # O(1/Nw) and O(1/Nw^2) removed
    assert abs(lI-ext) < 1e-5*ext and abs(lI-ext) < 0.01*abs(l3-ext), (lI, ext)


def test_to_uniform_fermion_evaluates_a_self_energy_with_a_constant_part():
    """Sigma(iw) = S0 + sum_p A_p/(iw - e_p) (Hermitian S0, A_p): the constant S0 is not representable by the basis,
    to_uniform_fermion takes it off and the rest goes through the IR basis to the uniform w_n."""
    if not HAVE_IR:
        T.skip('sparse-ir not installed')
    temp, No = 0.05, 2
    rng = np.random.default_rng(1)
    S0 = rng.standard_normal((No, No))+1j*rng.standard_normal((No, No))
    S0 = S0+S0.conj().T
    poles = [(-1.3, 0.4), (0.2, 0.7), (2.1, 0.3)]
    vecs = [rng.standard_normal(No)+1j*rng.standard_normal(No) for _ in poles]

    def sigma(w):                                  # [No, No, len(w)]; array layout [a, b] = Sigma_ba, as flibs
        s = np.repeat(S0[:, :, None], len(w), axis=2).astype(complex)
        for (e, a), v in zip(poles, vecs):
            s += a*np.outer(v, v.conj())[:, :, None]/(1j*w[None, None, :]-e)
        return s.swapaxes(0, 1)
    g = P.ir_grid(temp, 5.0, eps=1e-10)
    try:
        P.set_freq_grid(g)
        x = sigma(g['wF'])[..., None]              # one k-point
        u = P.to_uniform_fermion(x, temp, 64)[..., 0]
        with np.testing.assert_raises(ValueError):
            P.to_uniform_fermion(x[:, :, :-1], temp, 64)   # not the grid's Nw
    finally:
        P.set_freq_grid(None)
    ref = sigma((2*np.arange(64)+1)*np.pi*temp)
    assert T.rel_err(u, ref) < 1e-7
    assert np.array_equal(P.to_uniform_fermion(ref[..., None], temp, 8)[..., 0], ref[:, :, :8])   # uniform: a slice


def test_ir_nonlinear_eliashberg_gap_is_the_large_nw_limit():
    """Converged nonlinear gap (one orbital, lambda_lin = 1.34, fixed RPA vertex) on IR = the uniform grid's
    large-Nw limit, which the uniform grid approaches slowly (3% too small at Nw = 4096)."""
    if not HAVE_IR:
        return
    N, temp, U, mu = 16, 0.02, 2.5, -0.6

    def solve(Nw):
        st = T.one_orbital_square(Nx=N, Ny=N, Nw=Nw, temp=temp, mu=mu, tx=1.0, ty=1.0, U=U)
        with T.silence():
            chi, _ = F.get_chi0(st['Smat'], st['Cmat'], st['Gk'], st['olist'], st['kmap'], st['invk'], temp, N, N, 1)
            d0, _ = F.linearized_eliashberg(chi, st['Gk'], st['uni'], P.get_initial_gap(st['klist'], 1, 1),
                                            st['Smat'], st['Cmat'], st['olist'], st['plist'], st['kmap'], st['invk'],
                                            N, N, 1, temp, 1, eps=1e-8)
            d0 *= 1.764*temp/np.abs(d0).max()
            d, _ = F.nonlinear_eliashberg(np.ascontiguousarray(d0), st['Gk'], st['hamk'], st['Smat'], st['Cmat'],
                                          st['olist'], st['plist'], st['kmap'], st['invk'], mu, temp, 1, N, N, 1,
                                          sw_Vconst=True, eps=1e-9, itemax=400)
        return d[0, 0]
    g = P.ir_grid(temp, 12.0, eps=1e-8)
    try:
        P.set_freq_grid(g)
        dI = solve(g['Nw'])
    finally:
        P.set_freq_grid(None)
    n = np.rint((g['wF']/(np.pi*temp)-1)/2).astype(int)
    keep = n < 1024
    u = {Nw: solve(Nw)[n[keep]] for Nw in (1024, 2048, 4096)}
    ext = (8*u[4096]-6*u[2048]+u[1024])/3
    assert T.rel_err(dI[keep], ext) < 5e-4
    assert T.rel_err(dI[keep], ext) < 0.05*T.rel_err(u[4096], ext)


def test_ir_nonlinear_map_with_the_vertex_rebuilt_from_the_sc_state():
    """One step of the nonlinear map with sw_Vconst=False (vertex from the SC-state chi0 = GG +- FF bubbles) for the
    same analytic seed Delta(k, iw) = D0 (cos kx - cos ky) w0^2/(w^2 + w0^2): IR = the uniform grid's large-Nw limit.
    With m_diis=1 and no amplitude Newton the first update is Delta = 0.3 M[seed] + 0.7 seed."""
    if not HAVE_IR:
        return
    N, temp, U, mu, D0, w0 = 16, 0.02, 2.5, -0.6, 0.01, 0.3

    def M(wn, Nw, vconst):
        st = T.one_orbital_square(Nx=N, Ny=N, Nw=Nw, temp=temp, mu=mu, tx=1.0, ty=1.0, U=U)
        kx, ky = 2*np.pi*st['klist'][:, 0], 2*np.pi*st['klist'][:, 1]
        seed = (D0*(np.cos(kx)-np.cos(ky))[None, :]*(w0**2/(wn**2+w0**2))[:, None]).astype(complex)[None, None]
        with T.silence():
            d, _ = F.nonlinear_eliashberg(np.ascontiguousarray(seed.copy()), st['Gk'], st['hamk'], st['Smat'],
                                          st['Cmat'], st['olist'], st['plist'], st['kmap'], st['invk'], mu, temp, 1,
                                          N, N, 1, sw_Vconst=vconst, eps=1e-14, itemax=1, m_diis=1, gap_min=0.0,
                                          sw_amp_newton=False)
        return ((d-0.7*seed)/0.3)[0, 0]
    g = P.ir_grid(temp, 12.0, eps=1e-8)
    try:
        P.set_freq_grid(g)
        mI, mIc = M(g['wF'], g['Nw'], False), M(g['wF'], g['Nw'], True)
    finally:
        P.set_freq_grid(None)
    n = np.rint((g['wF']/(np.pi*temp)-1)/2).astype(int)
    keep = n < 1024
    u = {Nw: M((2*np.arange(Nw)+1)*np.pi*temp, Nw, False)[n[keep]] for Nw in (1024, 2048, 4096)}
    ext = (8*u[4096]-6*u[2048]+u[1024])/3
    assert T.rel_err(mI[keep], ext) < 2e-4
    assert T.rel_err(mI, mIc) > 0.05                     # the FF bubble does enter (not a trivial check)


def test_ir_nonlinear_eliashberg_with_sigma_and_vertex_update_converges():
    """Strongly coupled one-band case (FLEX lambda = 1.23, Stoner 0.965, Delta ~ 8 pi T) with Sigma and the vertex
    rebuilt from the SC state: converges from a small seed, to the same gap for two mixing rates.  Needs the exact
    Nambu G (the former G = G0(1 + Delta F^+) sub-iteration diverged for Delta >> pi T), the static-Sigma
    subtraction (sub_sigma=1; mu is fixed in the loop) and the mixing of Sigma."""
    if not HAVE_IR:
        return
    N, temp, U, fill = 16, 0.01, 3.5, 0.41
    g = P.ir_grid(temp, 12.0, eps=1e-8)
    try:
        P.set_freq_grid(g)
        st = T.one_orbital_square(Nx=N, Ny=N, Nw=g['Nw'], temp=temp, mu=0.0, tx=1.0, ty=1.0, U=U)
        with T.silence():
            s, mu = F.mkself(st['Smat'], st['Cmat'], st['kmap'], st['invk'], st['olist'], st['hamk'], st['eig'],
                             st['uni'], -0.5, fill, temp, g['Nw'], N, N, 1, False, False, eps=1e-10, scf_loop=500)
            G = F.gen_green(s, st['hamk'], mu, temp)
            chi, _ = F.get_chi0(st['Smat'], st['Cmat'], G, st['olist'], st['kmap'], st['invk'], temp, N, N, 1)
            dl, lam = F.linearized_eliashberg(chi, G, st['uni'], P.get_initial_gap(st['klist'], 1, 1), st['Smat'],
                                              st['Cmat'], st['olist'], st['plist'], st['kmap'], st['invk'], N, N, 1,
                                              temp, 1, eps=1e-10)
            gaps = []
            for pp in (0.3, 0.15):
                d, _ = F.nonlinear_eliashberg(np.ascontiguousarray(dl*1e-4/np.abs(dl).max()), G.copy(),
                                              st['hamk'], st['Smat'], st['Cmat'], st['olist'], st['plist'],
                                              st['kmap'], st['invk'], mu, temp, 1, N, N, 1, sw_sigma=True,
                                              sw_Vconst=False, eps=1e-10, itemax=3000, gap_min=0.0,
                                              sw_amp_newton=False, m_diis=1, pp=pp)
                gaps.append(d[0, 0])
    finally:
        P.set_freq_grid(None)
    assert lam > 1.0
    assert np.all(np.isfinite(gaps[0])) and np.abs(gaps[0][0]).max() > 5*np.pi*temp
    assert T.rel_err(gaps[1], gaps[0]) < 1e-7


if __name__ == '__main__':
    sys.exit(T.run_standalone(globals()))
