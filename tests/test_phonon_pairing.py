"""
Phonon-mediated pairing from EPA data (plibs.epa_pairing_vertex) and its entry into the
linearized Eliashberg equation (flibs.linearized_eliashberg(..., vph=...)).

  * vph enters V_Delta with the sign and weight of the spin/charge vertex: for a single band at
    U = 0 the Eliashberg eigenvalue equals that of the isotropic gap equation
    lambda Delta_n = -T sum_m W(w_n - w_m) <1/(w_m^2 + xi^2)>_k Delta_m, solved directly;
  * epa_pairing_vertex builds W from the EPA table so that the Fermi-surface average of its
    band-basis amplitude is -sum_nu <g^2>_FS D_nu, i.e. lambda = 2 N_F sum_nu <g^2>/w_nu;
  * the overlap factor eta equals the explicit double Fermi-surface sum it stands for.
"""
import io
import os
import sys
import contextlib
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import libs.plibs as P  # noqa: E402
import libs.flibs as F  # noqa: E402


def _square(t=-0.25):
    rvec = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]], dtype=np.float64)
    ham_r = np.zeros((5, 1, 1), dtype=np.complex128)
    ham_r[1:, 0, 0] = t
    return rvec, ham_r


def _two_orbital(hyb=0.0):
    """Two square-lattice orbitals on one or two sites, coupled by an on-site hybridization."""
    rvec, h1 = _square()
    ham_r = np.zeros((5, 2, 2), dtype=np.complex128)
    ham_r[:, 0, 0] = h1[:, 0, 0]
    ham_r[:, 1, 1] = 0.8*h1[:, 0, 0]
    ham_r[0, 1, 1] = 0.05
    ham_r[0, 0, 1] = ham_r[0, 1, 0] = hyb
    return rvec, ham_r


def _write_epa(path, edge=(0.4, 0.45), step=(-0.1, 0.1), nbin=(8, 8), wcm=(200.0, 350.0), seed=0):
    rng = np.random.default_rng(seed)
    nm = len(wcm)
    lines = [f'{2:8d}{nm:8d}'] + [f'{e:14.8f}{s:14.8f}{n:8d}' for e, s, n in zip(edge, step, nbin)]
    lines.append(''.join(f'{w:18.8e}' for w in wcm))
    for ii in range(2):
        for jj in range(nbin[ii]):
            for kk in range(nbin[ii]):
                lines.append(f'{ii+1:8d}{jj+1:8d}{kk+1:8d}'
                             + ''.join(f'{x:18.8e}' for x in rng.uniform(1e-3, 5e-3, nm)))
    open(path, 'w').write('\n'.join(lines)+'\n')


def test_vph_enters_the_eliashberg_kernel_as_a_retarded_attraction():
    Nx = Ny = 16
    Nw, temp, mu = 64, 0.02, -0.1
    rvec, ham_r = _square()
    klist, kmap, invk = F.gen_irr_k_TRS(Nx, Ny, 1)
    eig, uni = F.get_eig(F.gen_ham(klist, ham_r, rvec))
    Gk = F.gen_Green0(eig, uni, mu=mu, temp=temp, Nw=Nw)
    olist, site = P.get_chi_orb_list(1, [1])
    Smat = Cmat = np.zeros((1, 1))
    vj = 2*np.pi*temp*np.arange(Nw)
    W = -0.004*2*0.03/(vj**2+0.03**2)
    with contextlib.redirect_stdout(io.StringIO()):
        chi, _ = F.get_chi0(Smat, Cmat, Gk, olist, kmap, invk, temp, Nx, Ny, 1)
        _, lam = F.linearized_eliashberg(chi, Gk, uni, P.get_initial_gap(klist, 1, 0), Smat, Cmat, olist,
                                         F.get_plist(rvec, ham_r), kmap, invk, Nx, Ny, 1, temp, gap_sym=0,
                                         eps=1e-9, vph=W[None, None, :])
    _, klf = P.gen_klist(Nx, Ny, 1, sw_pp=False)
    xi = P.get_eigs(klf, ham_r, [], rvec)[0][:, 0]-mu
    n = np.arange(-Nw, Nw)
    wn = (2*n+1)*np.pi*temp
    chim = (1/(wn[None, :]**2+xi[:, None]**2)).mean(axis=0)
    d = (n[:, None]-n[None, :]) % (2*Nw)          # the circular 2Nw convolution of the solver
    Wc = np.where(d < Nw, W[np.minimum(d, Nw-1)], np.where(d == Nw, 0.0, W[np.minimum(2*Nw-d, Nw-1)]))
    ref = np.linalg.eigvals(-temp*Wc*chim[None, :]).real.max()
    assert lam > 0 and abs(lam-ref) < 1e-8*ref


def test_epa_vertex_of_one_band_is_the_fermi_surface_coupling(tmp_path):
    _write_epa(str(tmp_path/'x.epa.e'))
    rvec, ham_r = _square()
    olist, site = P.get_chi_orb_list(1, [1])
    temp, mu, Nw = 0.01, 0.42, 16
    with contextlib.redirect_stdout(io.StringIO()):
        vph, info = P.epa_pairing_vertex(str(tmp_path/'x.epa.e'), ham_r+np.eye(1)*0.42*(np.arange(5) == 0)[:, None, None],
                                         [], rvec, 24, 24, 1, olist, site, mu, temp, Nw, gap_sym=0)
    assert abs(info['eta']-1.0) < 1e-12                       # one orbital: no overlap reduction
    _, _, edge, step, nbin, wavg, gavg = P.read_epa_output(str(tmp_path/'x.epa.e'))
    vj = 2*np.pi*temp*np.arange(Nw)
    W = -(info['g2_FS'][:, None]*2*wavg[:, None]/(vj**2+wavg[:, None]**2)).sum(axis=0)
    assert np.abs(vph[0, 0]-W).max() < 1e-14
    assert abs(info['lambda']-2*info['N_F']*np.sum(info['g2_FS']/wavg)) < 1e-12
    assert info['g2_FS'].min() >= gavg[gavg > 0].min() and info['g2_FS'].max() <= gavg.max()


def test_eta_is_the_fermi_surface_average_of_the_site_projected_overlaps(tmp_path):
    """eta = <sum_site |<u'|P_site|u>|^2>_FS,FS, the double sum written out, for one and two sites;
    two decoupled orbitals on one site give (Wa^2 + Wc^2)/(Wa + Wc)^2 with Wa, Wc their Fermi-window
    weights (no pair hopping between them: 1/2 for equal weights)."""
    _write_epa(str(tmp_path/'x.epa.e'))
    temp, mu, N = 0.02, 0.0, 12
    for hyb, prof in ((0.0, [2]), (0.07, [2]), (0.07, [1, 1])):
        rvec, ham_r = _two_orbital(hyb)
        olist, site = P.get_chi_orb_list(2, prof)
        with contextlib.redirect_stdout(io.StringIO()):
            _, info = P.epa_pairing_vertex(str(tmp_path/'x.epa.e'), ham_r, [], rvec, N, N, 1, olist, site,
                                           mu, temp, 4, gap_sym=0)
        _, _, eig, uni, _ = P.get_emesh(N, N, 1, ham_r, [], rvec, np.eye(3), sw_uni=True)
        w = (0.25/temp/np.cosh(0.5*(eig-mu)/temp)**2).ravel()
        u = np.conj(uni).reshape(-1, 2)                         # rows: <a|u_kn>
        blocks = [[0, 1]] if len(prof) == 1 else [[0], [1]]
        amp = sum(np.abs(u[:, b].conj().dot(u[:, b].T))**2 for b in blocks)
        ref = w.dot(amp).dot(w)/w.sum()**2
        assert abs(info['eta']-ref) < 1e-10*ref, (hyb, prof)
        if hyb == 0.0:                                          # decoupled: (Wa^2 + Wc^2)/(Wa + Wc)^2
            wb = (0.25/temp/np.cosh(0.5*(eig-mu)/temp)**2).sum(axis=0)
            assert abs(info['eta']-(wb**2).sum()/wb.sum()**2) < 1e-6  # eigenvectors mix at crossings


if __name__ == '__main__':
    import _tools
    sys.exit(_tools.run_standalone(globals()))
