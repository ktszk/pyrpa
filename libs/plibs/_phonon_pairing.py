"""
Phonon-mediated pairing interaction from the EPA (electron-phonon averaged) data of epa.x,
to be added to the spin/charge-fluctuation pairing of the linearized Eliashberg equation.
"""
import numpy as np
import libs.flibs as flibs
import libs.plibs._bands as _bands
import libs.plibs._freq_grid as _freq_grid
import libs.plibs._wannier_io as _wannier_io

# Public API of this submodule; keeps `from ._x import *` in __init__.py from
# re-exporting numpy and the sibling packages imported above.
__all__ = ['epa_pairing_vertex']


def _epa_bins(e: np.ndarray, edge: np.ndarray, step: np.ndarray, nbin: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """EPA grid (valence below the midgap, conduction above) and clamped bin of each energy, as calc_tau_epa."""
    ig = np.where(e < edge.sum()/len(edge), 0, len(edge)-1)
    x = np.clip((e-edge[ig])/step[ig], 1.0e-12, nbin[ig]-1.0e-12)
    return ig, x.astype(np.int64)


def epa_pairing_vertex(epa_file: str, ham_r: np.ndarray, S_r: np.ndarray, rvec: np.ndarray, Nx: int, Ny: int,
                       Nz: int, chiolist: np.ndarray, site: np.ndarray, mu: float, temp: float, Nw: int,
                       gap_sym: int, scale: float = 1.0) -> tuple[np.ndarray, dict]:
    """
    @brief Phonon-mediated pairing interaction for flibs.linearized_eliashberg (its vph) from epa.x output.

    EPA averages |g_nu|^2 over the states of energy bins and the phonon frequency over q, so what
    it fixes is the isotropic Fermi-surface coupling: the band-basis pair-scattering amplitude
    averaged over the Fermi surface,
        V_FS(iv) = -sum_nu <g_nu^2>_FS D_nu(iv),   D_nu(iv) = 2 w_nu/(v^2 + w_nu^2),
    i.e. lambda_nu = 2 N_F <g_nu^2>_FS / w_nu (McMillan's lambda in the Einstein form, N_F per spin
    per unit cell).  It is carried here by a retarded density-density attraction W(iv), local in
    each site and the same for every orbital pair of the site (U = U' = W, J = 0 in the Kanamori
    S/C convention of the spin/charge vertices), and W is normalized so that this orbital
    interaction reproduces V_FS on the Fermi surface: its band-basis amplitude between states u, u'
    is W sum_site |<u'|P_site|u>|^2, whose Fermi-surface average is W*eta with
        eta = sum_site Tr(P_site^2)/(sum w)^2,   P = sum_{k,n} w_kn |u_kn><u_kn|,  w = -df/de,
    so W = V_FS/eta.  Being q-independent it feeds the s-wave (k-even, local) channel only, in
    both singlet and orbital-antisymmetric triplet form.  No phonon self-energy is included: the
    1 + lambda mass renormalization of a full Migdal-Eliashberg treatment is absent.

    <g^2>_FS is the double Fermi-surface average of the EPA table over the bins of the states;
    pairs across the two EPA grids are not tabulated, so a Fermi surface split between them
    (E_F at the shared edge of a metal's valence and conduction grids) is averaged over the
    tabulated pairs only and flagged.
    @param  epa_file: epa.x (job='egrid') output
    @param ham_r,S_r,rvec: Hamiltonian in Wannier form (S_r = [] for an orthogonal basis)
    @param  Nx,Ny,Nz: k-mesh of the Fermi-surface averages (the Eliashberg mesh)
    @param  chiolist: orbital-pair basis of chi [Nchi, 2] (1-based), site: its site labels [Nchi]
    @param        mu: chemical potential [eV], on the same absolute scale as the EPA grid edges
    @param      temp: k_B T [eV]
    @param        Nw: number of stored bosonic frequencies: v_j = 2 pi j T, j = 0..Nw-1, on the uniform grid,
                      the sampling frequencies on the IR grid (plibs.set_freq_grid; set it before this call)
    @param   gap_sym: >= 0 singlet, < 0 triplet (the channel lin_eliash solves)
    @param     scale: factor on the phonon coupling (1 = EPA as is)
    @return vph [Nchi, Nchi, Nw] float64 and a dict (lambda_nu, lambda, omega_log, N_F, g2_FS, eta, W0)
    """
    ngrid, nmodes, edge, step, nbin, wavg, gavg = _wannier_io.read_epa_output(epa_file)
    Nk, _, eig, uni, _ = _bands.get_emesh(Nx, Ny, Nz, ham_r, S_r, rvec, np.eye(3), sw_uni=True)
    w = 0.25/temp/np.cosh(0.5*(eig-mu)/temp)**2                  # -df/de [Nk, Norb]
    wsum = w.sum()
    if wsum < 1.0e-12*w.size/temp:
        raise ValueError(f'epa_pairing_vertex: no states within the Fermi window at mu = {mu:.4f} eV')
    N_F = wsum/Nk
    # --- <g^2>_FS: double average of the EPA table with the Fermi-window weights of the bins
    ig, ib = _epa_bins(eig.ravel(), edge, step, nbin)
    h = np.zeros((ngrid, int(nbin.max())))
    np.add.at(h, (ig, ib), w.ravel())
    tab = sum(h[g].sum()**2 for g in range(ngrid))                 # weight of the tabulated pairs
    g2 = np.array([sum(h[g].dot(gavg[g, :, :, nu]).dot(h[g]) for g in range(ngrid)) for nu in range(nmodes)])/tab
    split = 1.0-tab/wsum**2
    # --- eta: Fermi-surface average of sum_site |<u'|P_site|u>|^2 for the site-local orbital W
    orb_site = {int(a): int(s) for (a, b), s in zip(np.asarray(chiolist), site) if a == b}
    c = np.conj(uni)                                               # c[k, n, a] = <a|u_kn>
    P = np.einsum('kn,kna,knb->ab', w, c, c.conj())
    eta = 0.0
    for s in set(orb_site.values()):
        idx = [a-1 for a, t in orb_site.items() if t == s]
        eta += np.sum(np.abs(P[np.ix_(idx, idx)])**2).real
    eta /= wsum**2
    # --- W(iv_j) and the vertex in the channel gap_sym selects
    ok = wavg > 1.0e-8
    vj = _freq_grid.boson_frequencies(temp, Nw)
    VFS = -scale*(g2[ok, None]*2.0*wavg[ok, None]/(vj[None, :]**2+wavg[ok, None]**2)).sum(axis=0)
    W = VFS/eta
    Norb = ham_r.shape[1]
    S1, C1 = flibs.gen_SCmatrix_orb(chiolist, site, np.ones((Norb, Norb)), np.zeros((Norb, Norb)))
    V1 = 0.5*(S1+C1) if gap_sym >= 0 else 0.5*(S1-C1)
    vph = np.ascontiguousarray(V1[:, :, None]*W[None, None, :])
    lam_nu = np.zeros(nmodes)
    lam_nu[ok] = scale*2.0*N_F*g2[ok]/wavg[ok]
    lam = lam_nu.sum()
    wlog = np.exp((lam_nu[ok]*np.log(wavg[ok])).sum()/lam) if lam > 0 else 0.0
    print(f'EPA phonon pairing from {epa_file}: N_F = {N_F:.4f} /eV/spin/cell, lambda = {lam:.4f}, '
          f'omega_log = {wlog*1e3:.2f} meV ({nmodes} modes, scale = {scale:g})', flush=True)
    print(f'  eta = {eta:.4f} (orbital/site overlap on the Fermi surface), W(iv=0) = {W[0]:.4e} eV '
          f'({"singlet" if gap_sym >= 0 else "triplet"} channel)', flush=True)
    if split > 1.0e-3:
        print(f'  WARNING: {split*100:.1f}% of the Fermi-surface pair weight straddles the two EPA grids and '
              'is not tabulated; <g^2> is averaged over the tabulated pairs (put E_F inside one grid)',
              flush=True)
    return vph, {'lambda_nu': lam_nu, 'lambda': lam, 'omega_log': wlog, 'N_F': N_F, 'g2_FS': g2,
                 'eta': eta, 'W0': W[0]}
