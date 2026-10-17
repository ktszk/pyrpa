#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
Boltzmann / linear-response transport drivers: relaxation-time models, weak-field Hall
coefficient with carrier analysis, thermoelectric tensors, optical constants and colours.
"""
import os
import numpy as np
import scipy.linalg as sclin
import scipy.constants as scconst
import libs.flibs as flibs
import libs.plibs._units as _units
import libs.plibs._bands as _bands
import libs.plibs._calc as _calc
import libs.plibs._color as _color
import libs.plibs._optical as _optical
import libs.plibs._response as _response
import libs.plibs._wannier_io as _wannier_io

# Public API of this submodule; keeps `from ._x import *` in __init__.py from
# re-exporting numpy, scipy, ctypes and the sibling packages imported above.
__all__ = [
    'print_matrix', 'build_tau', 'print_tau_summary', 'get_hall_coe', 'calc_conductivity_Boltzmann',
    'calc_conductivity_lrt'
]


def print_matrix(label: str, mat: np.ndarray, width: float = 10, prec: int = 3, ndigits: int = 10) -> None:
    """Print a labelled N x 3 real matrix, one row per line (transport tensors)."""
    print(label, flush=True)
    for row in np.asarray(mat).round(ndigits):
        print(f" {row[0]:{width}.{prec}e} {row[1]:{width}.{prec}e} {row[2]:{width}.{prec}e}", flush=True)


def build_tau(tau_mode: str, units: _units.Units, eig: np.ndarray, mu: float, temp: float, tau_const: float,
              rvec: np.ndarray | None = None, ham_r: np.ndarray | None = None, S_r: np.ndarray | None = None,
              avec: np.ndarray | None = None, Nx: int | None = None, Ny: int | None = None,
              Nz: int | None = None, Nw: int = 300, epa_file: str | None = None,
              dos_delta: float = 1.0e-2) -> tuple:
    '''
    Relaxation time tau[Nk,Norb] in units of tau_unit, per the tau_mode switch.

    R_H is INDEPENDENT of a constant tau (it carries tau^2/tau^2), so the whole
    temperature dependence of the Hall coefficient -- and any difference between
    the mobilities of two Fermi-surface sheets -- lives in tau(k).  With
    tau_mode='const' the R_H printed by get_hall_coe is a pure band-structure
    number; only a k-dependent tau makes R_H(T) meaningful.

    tau_mode: 'const', 'e2' (Fermi liquid), 'dos' (DOS ratio) or 'epa' (needs epa_file);
    dos_delta is the Lorentzian broadening [eV] of the DOS used by 'dos'.
    Returns (tau, label).
    '''
    if tau_mode == 'const':
        return eig*0.+tau_const, f'constant tau = {tau_const} '+units.tau_label
    if tau_mode == 'epa':
        if not epa_file:
            raise ValueError("tau_mode='epa' needs epa_file (epa.x job='egrid' output)")
        ngrid, nmodes, edge, step, nbin, wavg, gavg = _wannier_io.read_epa_output(epa_file)
        # calc_tau_epa works with hbar=1 and energies in eV, so it returns hbar/eV;
        # the Boltzmann kernels expect tau in units of tau_unit.
        tau = flibs.calc_tau_epa(eig, gavg, wavg, edge, step, nbin, mu, temp)*units.hbar/units.tau_unit
        return tau, f'EPA tau(k) from {epa_file} ({nmodes} modes)'
    if tau_mode == 'e2':
        # Fermi-liquid electron-electron scattering: 1/tau ∝ (e-mu)^2 + (pi T)^2, normalized so
        # that tau = tau_const on the Fermi level and falls off as 1/(e-mu)^2 away from it
        pit2 = (np.pi*temp)**2
        return (tau_const*pit2/((eig-mu)**2+pit2),
                f'Fermi-liquid tau: 1/tau ~ (e-mu)^2+(pi T)^2, tau(mu) = {tau_const} '
                + units.tau_label)
    if tau_mode == 'dos':
        # golden-rule scattering rate ∝ final-state DOS, used only as a RATIO to the Fermi level:
        # tau(e) = tau_const * D(mu)/D(e), so tau = tau_const on the Fermi level
        _, _, eig2, uni, _ = _bands.get_emesh(Nx, Ny, Nz, ham_r, S_r, rvec, avec.T, sw_uni=True)
        wlist = np.linspace(eig.min()-mu, eig.max()-mu, Nw, True)   # energies relative to mu
        dos = flibs.gen_dos(eig2, uni, mu, wlist, dos_delta).sum(axis=0)
        dos_mu = np.interp(0.0, wlist, dos)
        if dos_mu <= 1.0e-8*dos.max():
            raise ValueError("tau_mode='dos': the DOS vanishes at mu, so tau(e) = tau_const*D(mu)/D(e) "
                             "is undefined (insulator?)")
        # the broadened DOS is positive at every band energy; the floor only guards round-off
        dos_e = np.maximum(np.interp(eig-mu, wlist, dos), 1.0e-12*dos.max())
        return (tau_const*dos_mu/dos_e,
                f'DOS-ratio tau: tau = tau_const*D(mu)/D(e), tau(mu) = {tau_const} '
                + units.tau_label)
    raise ValueError(f"unknown tau_mode '{tau_mode}' (use 'const', 'e2', 'dos' or 'epa')")


def print_tau_summary(tau: np.ndarray, eig: np.ndarray, mu: float, temp: float, label: str, indent: str = '  ',
                      tau_label: str = 'fs') -> None:
    '''Report tau where it matters: weighted by the Fermi window -df/de, not over all bands.'''
    print(f'{indent}tau: {label}', flush=True)
    w = 0.25*(1.-np.tanh(0.5*(eig-mu)/temp)**2)/temp
    if w.sum() <= 0:
        return
    tav = (w*tau).sum()/w.sum()
    sel = w > 0.01*w.max()          # states actually inside the Fermi window
    unit = tau_label
    print(f'{indent}     <tau>_FS = {tav:.4g} {unit}   range over the Fermi window '
          f'[{tau[sel].min():.4g}, {tau[sel].max():.4g}] {unit}', flush=True)
    if tau.std() > 0:
        # R_H picks up <tau^2>/<tau>^2, which is 1 only for a constant tau
        r = ((w*tau*tau).sum()/w.sum())/tav**2
        print(f'{indent}     <tau^2>/<tau>^2 = {r:.4f} '
              + ('(constant tau: R_H is tau-independent)' if abs(r-1.) < 1e-6
                 else '(>1: tau(k) shifts R_H away from the constant-tau value)'), flush=True)


def _hall_at_mesh(rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray, avec: np.ndarray, Nx: int, Ny: int, Nz: int,
                  fill: float, temp: float, tau_const: float, gsp: float, units: _units.Units, tau_kw: dict,
                  Nw: int = 300, detail: bool = False) -> dict:
    '''
    Weak-field (B-linear) Boltzmann Hall response on ONE k-mesh.
    detail=True additionally builds the band-resolved (two-fluid) decomposition
    and the Luttinger-volume carrier densities; that part is only worth its cost
    on the mesh actually being reported, not on every step of a convergence sweep.
    Returns a dict; all printing is done by get_hall_coe.
    '''
    Nk, eig, vk, imass, kweight = _bands.get_emesh(Nx, Ny, Nz, ham_r, S_r, rvec, avec.T*units.ihbar,
                                                   sw_veloc=True, sw_mass=True)
    Vuc = sclin.det(avec)*1e-30
    wsum = kweight.sum()          # NOT Nk: stays correct if a symmetry-reduced mesh is ever used
    mu = _bands.calc_mu(eig, Nk, fill, temp, kweight)
    tau, tau_label = build_tau(units=units, eig=eig, mu=mu, temp=temp, tau_const=tau_const, rvec=rvec, ham_r=ham_r,
                               S_r=S_r, avec=avec, Nx=Nx, Ny=Ny, Nz=Nz, Nw=Nw, **tau_kw)
    K0, K1, K2 = flibs.calc_Kn(eig, vk, kweight, temp, mu, tau)
    # Full Hall kernel S[n,a,b,g] for the three field directions g, band resolved.
    # The tensor (rather than the single B||z scalar) is what allows the proper
    # rho^(1) = -sigma^-1 sigma^(1) sigma^-1 below instead of assuming sigma is
    # diagonal in xyz.
    Shall_n, Sabs_n = flibs.calc_sigma_hall(eig, vk, imass/units.eC, kweight, tau, temp, mu)
    Shall = Shall_n.sum(axis=0)
    Rh = _calc.hall_coefficient(K0, Shall, Vuc, wsum, gsp)          # m^3/C, per field axis
    nH, kind = _calc.hall_carrier_density(Rh, units.eC)                 # cm^-3, 'electron'/'hole'
    sigma = gsp*units.tau_unit*units.eC*K0/(wsum*Vuc)                       # S/m
    # sigma^(1)/B carries one more power of e and of tau than sigma  -> S/(m T)
    sigma1 = gsp*units.eC*units.eC*units.tau_unit*units.tau_unit*Shall/(wsum*Vuc)
    # Hall mobility mu_H = R_H * sigma_aa with the current in the plane normal to B.
    # omega_c*tau = mu_H[SI]*B, so mu_H directly measures how far the weak-field
    # expansion this whole routine is built on can be trusted.
    muH = np.array([Rh[g]*sigma[(g+1)%3, (g+1)%3] for g in range(3)])   # m^2/(V s)
    # --- two k-mesh adequacy diagnostics -----------------------------------
    dfermi = 0.25*(1.-np.tanh(0.5*(eig-mu)/temp)**2)/temp
    # (i) effective number of states carrying the Fermi window, as the
    #     participation ratio of the -df/de weights.  Everything here is an
    #     average over Neff states, so Neff sets the statistical noise floor.
    Neff = dfermi.sum()**2/(dfermi**2).sum()
    # (ii) cancellation ratio of the Hall integrand: it changes sign around the
    #     Fermi surface, and when |sum| << sum|.| the total is a small residue of
    #     large opposing parts -- exactly the regime where R_H keeps changing
    #     sign with the mesh.  sigma_xx has no such cancellation, which is why it
    #     looks converged long before sigma^(1) is.
    cancel = abs(Shall[0, 1, 2])/max(Sabs_n[:, 0, 1, 2].sum(), 1e-300)
    out = {'mesh': (Nx, Ny, Nz), 'Nk': Nk, 'mu': mu, 'Rh': Rh, 'nH': nH, 'kind': kind,
           'sigma': sigma, 'sigma1': sigma1, 'muH': muH, 'Neff': Neff, 'cancel': cancel, 'K0': K0,
           'tau_label': tau_label}
    if detail:
        # keep the per-state arrays only for the mesh that gets reported: holding
        # eig/tau for every step of a sweep is hundreds of MB on a fine mesh
        out['tau'], out['eig'] = tau, eig
        # --- two-fluid decomposition -----------------------------------------
        # sigma and sigma^(1) are both plain band sums, so each sheet can be given
        # its own conductivity, Hall coefficient and mobility.  Comparing a band's
        # OWN n_H against its Luttinger volume shows whether that sheet behaves as
        # one simple closed pocket; the mismatch between the total n_H and n_e/n_h
        # is then the compensation the experiment is really seeing.
        K0n = flibs.calc_k0_band(eig, vk, kweight, temp, mu, tau)                   # [Norb,3,3]
        Norb = len(K0n)
        trK0 = np.trace(K0)
        Rh_n = np.full(Norb, np.nan)
        for n in range(Norb):
            if np.trace(K0n[n]) > 1e-8*trK0:
                Rh_n[n] = _calc.hall_coefficient(K0n[n], Shall_n[n], Vuc, wsum, gsp)[2]
        out['car'] = _calc.fs_carrier_density(eig, mu, Vuc, gsp, mesh=(Nx, Ny, Nz))
        out['sigma_n'] = gsp*units.tau_unit*units.eC*K0n/(wsum*Vuc)                        # [Norb,3,3] S/m
        out['Rh_n'] = Rh_n
    return out


def get_hall_coe(rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray, avec: np.ndarray, Nx: int, Ny: int, Nz: int,
                 fill: float, temp: float, tau_const: float, units: _units.Units, tau_kw: dict, Nw: int = 300,
                 with_spin: bool = False, mesh_list: list | None = None) -> list | None:
    '''
    Weak-field Hall coefficient R_H and the Hall carrier density n_H=1/(e|R_H|).

    R_H is formed from the FULL tensors, rho^(1)=-sigma^-1 sigma^(1) sigma^-1,
    for all three field directions.  Give mesh_list (e.g. [20,30,40] or
    [(20,20,10),...]) to run a k-mesh convergence sweep: sigma^(1) converges far
    more slowly than sigma and can come out with the WRONG SIGN on an
    under-converged mesh, so a single-mesh number is not usable on its own.
    tau_kw: keyword arguments of build_tau selecting the relaxation-time model
    (tau_mode, epa_file, dos_delta).
    '''
    # Parameter validation
    if temp <= 0:
        print("Error: Temperature (temp) is non-positive", flush=True)
        return
    if tau_const <= 0:
        print("Error: Relaxation time (tau_const) is non-positive", flush=True)
        return
    if sclin.det(avec) <= 0:
        print("Error: Unit cell volume (Vuc) is non-positive", flush=True)
        return
    gsp = (1.0 if with_spin else 2.0)  # spin degeneracy factor
    if mesh_list:
        meshes = [(m, m, m) if np.isscalar(m) else tuple(m) for m in mesh_list]
    else:
        meshes = [(Nx, Ny, Nz)]
    print(f"T = {temp/units.kb:.3f} K", flush=True)
    print(f"spin degeneracy gsp = {gsp:.0f} "
          f"({'SOC: both spins already in H' if with_spin else 'spin-degenerate bands'})", flush=True)
    print('  mesh          Nk      mu[eV]     Neff   cancel      R_H^x        R_H^y        R_H^z    [m^3/C]', flush=True)
    results = []
    for i, (nx, ny, nz) in enumerate(meshes):
        r = _hall_at_mesh(rvec, ham_r, S_r, avec, nx, ny, nz, fill, temp, tau_const, gsp, units, tau_kw, Nw,
                          detail=(i == len(meshes)-1))
        results.append(r)
        print(f'  {nx:3d}x{ny:3d}x{nz:3d} {r["Nk"]:9d}  {r["mu"]:9.4f} {r["Neff"]:8.3g} {r["cancel"]:8.2e} '
              + ' '.join(f'{v:+12.4e}' for v in r['Rh']), flush=True)
    res = results[-1]
    # --- k-mesh convergence verdict -----------------------------------------
    if len(results) > 1:
        prev = results[-2]['Rh']; last = res['Rh']
        print('  k-mesh convergence (finest vs previous):', flush=True)
        for g, ax in enumerate('xyz'):
            if abs(last[g]) < 1e-30:
                continue
            if prev[g]*last[g] < 0:
                msg = 'SIGN FLIP - not converged'
            else:
                d = abs(last[g]/prev[g]-1.)*100
                msg = f'{d:6.2f}% change'+('  <- not converged' if d > 5. else '')
            print(f'    B||{ax}: {msg}', flush=True)
    else:
        print('  (single mesh: pass mesh_list=[...] to check k-mesh convergence)', flush=True)
    # --- detail on the finest mesh ------------------------------------------
    print(f'--- {res["mesh"][0]}x{res["mesh"][1]}x{res["mesh"][2]} mesh, mu = {res["mu"]:.4f} eV ---', flush=True)
    print_tau_summary(res['tau'], res['eig'], res['mu'], temp, res['tau_label'], tau_label=units.tau_label)
    print_matrix('sigma matrix (S/m)', res['sigma'])
    for g, ax in enumerate('xyz'):
        print(f'  B||{ax}: R_H = {res["Rh"][g]:+.5e} m^3/C   '
              f'n_H = {res["nH"][g]:.5e} cm^-3 ({res["kind"][g]})   '
              f'mu_H = {1e4*res["muH"][g]:+.4e} cm^2/Vs   '
              f'omega_c*tau = {abs(res["muH"][g]):.3e} /T', flush=True)
    print(f'  sigma^(1)_xy/B = {res["sigma1"][0, 1, 2]:+.5e} S/(m T)  (B||z)', flush=True)
    # --- validity warnings ---------------------------------------------------
    print(f'  Fermi-window sampling: Neff = {res["Neff"]:.4g} states, '
          f'Hall-integrand cancellation |sum|/sum|.| = {res["cancel"]:.3e}', flush=True)
    if res['Neff'] < 1.0e3:
        print(f'  WARNING: only ~{res["Neff"]:.4g} states carry the Fermi window. '
              f'Refine the k-mesh (or raise T) before trusting R_H.', flush=True)
    if res['cancel'] < 1.0e-2:
        print(f'  WARNING: the Hall integrand cancels to {res["cancel"]:.2e} of its absolute weight. '
              f'R_H is a small residue of large opposing contributions here (near-compensation), '
              f'so it needs a much finer mesh than sigma_xx and its SIGN may still be mesh-dependent.', flush=True)
    Bmax = 0.1/max(abs(res['muH']).max(), 1e-30)
    print(f'  weak-field limit omega_c*tau<0.1 holds up to B ~ {Bmax:.3g} T', flush=True)
    # --- band-resolved (two-fluid) decomposition vs the Luttinger volume ------
    car = res['car']
    print('--- carrier densities: Luttinger volume vs Hall, band by band ---', flush=True)
    print('  band  type      n_e[cm^-3]    n_h[cm^-3]   pk e/h   sigma_xx[S/m]   R_H(band)'
          '    n_H(band)   mu_H[cm^2/Vs]', flush=True)
    for n, k in enumerate(car['kind']):
        if k in ('full', 'empty'):
            continue
        rn = res['Rh_n'][n]
        if np.isnan(rn) or abs(rn) < 1e-30:
            hall = f'{"-":>12s} {"-":>12s} {"-":>13s}'
        else:
            hall = (f'{rn:+12.4e} {1./(abs(rn)*units.eC)/1e6:12.4e}'
                    f' {1e4*abs(rn*res["sigma_n"][n, 0, 0]):13.4e}')
        print(f'  {n+1:4d}  {k:8s}  {car["ne_band"][n]:.5e}   {car["nh_band"][n]:.5e}  '
              f'{car["npocket_e"][n]:2d}/{car["npocket_h"][n]:<2d}   '
              f'{res["sigma_n"][n, 0, 0]:.5e}  {hall}', flush=True)
    ne, nh = car['n_e'], car['n_h']
    print(f'  Luttinger total: n_e = {ne:.5e}   n_h = {nh:.5e} cm^-3'
          + (f'   n_e/n_h = {ne/nh:.4f}' if nh > 0 else ''), flush=True)
    print(f'  Hall (all bands, B||z): n_H = {res["nH"][2]:.5e} cm^-3 ({res["kind"][2]})'
          '   <- the number to compare with experiment', flush=True)
    if (car['kind'] == 'open').any():
        print('  note: an OPEN sheet percolates the zone and has no Luttinger volume, so the '
              'n_e/n_h totals above miss whatever it carries.', flush=True)
    if ne > 0 and nh > 0 and abs(ne/nh-1.) < 0.2:
        print('  note: nearly compensated (n_e ~ n_h). n_H is then a small difference of '
              'large terms, extremely sensitive to mu, and is NOT a carrier density -- '
              'compare mobilities/two-fluid fits with experiment instead.', flush=True)
    return results


def calc_conductivity_Boltzmann(rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray, avec: np.ndarray, Nx: int,
                                Ny: int, Nz: int, fill: float, temp: float, tau_const: float, units: _units.Units,
                                tau_kw: dict, Nw: int = 300, with_spin: bool = False,
                                sw_tdf: bool = False) -> dict | None:
    '''
    calculate conductivities using Boltzmann equations
    tau_kw: keyword arguments of build_tau (tau_mode, epa_file, dos_delta)
    sw_tdf: compute the transport distribution function first, then energy-integrate
    Returns a dict of the tensors (sigma, kappa, kappa2, sigmaS, Seebeck, Peltier, Lorenz,
    Lorenz2, PF, mu); with sw_tdf also the transport distribution 'tdf' [Nw,3,3] on the
    energies 'tdf_e' (relative to mu).  None if the input is invalid.
    '''
    # no dep. T and mu or filling
    # Parameter validation
    if temp <= 0:
        print("Error: Temperature (temp) is non-positive", flush=True)
        return
    if tau_const <= 0:
        print("Error: Relaxation time (tau_const) is non-positive", flush=True)
        return

    Nk, eig, vk, kweight = _bands.get_emesh(Nx, Ny, Nz, ham_r, S_r, rvec, avec.T*units.ihbar, sw_veloc=True)
    Vuc = sclin.det(avec)*1e-30
    if Vuc <= 0:
        print("Error: Unit cell volume (Vuc) is non-positive", flush=True)
        return
    gsp = (1.0 if with_spin else 2.0)  # spin weight
    # calc_Kn / calc_tdf return kweight-weighted sums, so normalize by the total weight (as the
    # Hall kernels do), not by Nk: identical on the uniform mesh, correct on a weighted one
    wsum = kweight.sum()
    iNV = 1./(wsum*Vuc)
    itemp = 1./temp
    mu = _bands.calc_mu(eig, Nk, fill, temp, kweight)
    tau, tau_label = build_tau(units=units, eig=eig, mu=mu, temp=temp, tau_const=tau_const, rvec=rvec, ham_r=ham_r,
                               S_r=S_r, avec=avec, Nx=Nx, Ny=Ny, Nz=Nz, Nw=Nw, **tau_kw)
    print(f"T = {temp/units.kb:.3f} K", flush=True)
    print(f"mu = {mu:.4f} eV", flush=True)
    print_tau_summary(tau, eig, mu, temp, tau_label, indent='', tau_label=units.tau_label)
    if sw_tdf:
        tdf = flibs.calc_tdf(eig, vk, kweight, tau, Nw)
        # absolute-energy grid matching the calc_tdf bin CENTRES: E_i = emin + (i-1/2)*dw
        dw = (eig.max()-eig.min())/Nw
        wlist = eig.min()+dw*(np.arange(1, Nw+1)-0.5)
        dfermi = 0.25*(1.-np.tanh(0.5*(wlist-mu)/temp)**2)/temp
        K0 = (dfermi*tdf.T).T.sum(axis=0)*dw
        K1 = (dfermi*(wlist-mu)*tdf.T).T.sum(axis=0)*dw
        K2 = (dfermi*(wlist-mu)**2*tdf.T).T.sum(axis=0)*dw
    else:
        K0, K1, K2 = flibs.calc_Kn(eig, vk, kweight, temp, mu, tau)
    # sigma [S/m] = gsp * tau[s] * e[C] * K0[eV·m^-3] / (sum_k w_k * V[m^3])
    sigma = gsp*units.tau_unit*units.eC*K0*iNV  # e cancels eV->J so units work out to S/m
    kappa = gsp*units.tau_unit*units.eC*units.kb*K2*iNV*itemp

    # Handle sclin.inv() failures
    try:
        K0_inv = sclin.inv(K0)
        kappa2 = gsp*units.tau_unit*units.eC*units.kb*(K2-K1.dot(K0_inv.dot(K1)))*iNV*itemp
        Seebeck = -units.kb*K0_inv.dot(K1)*itemp
        Peltier = K1.dot(K0_inv)
    except np.linalg.LinAlgError:
        print("Error: K0 is a singular matrix. Cannot compute Seebeck coefficient and Peltier", flush=True)
        kappa2 = kappa.copy()
        Seebeck = np.zeros_like(K0)
        Peltier = np.zeros_like(K0)

    try:
        sigma_inv = sclin.inv(sigma*temp)
        # Lorenz tensor L = kb * kappa . (sigma*T)^-1 (matrix product, not element-wise)
        Lorenz = units.kb*kappa.dot(sigma_inv)
        Lorenz2 = units.kb*kappa2.dot(sigma_inv)
    except np.linalg.LinAlgError:
        print("Warning: sigma*temp is a singular matrix. Cannot compute Lorenz coefficient", flush=True)
        Lorenz = np.zeros_like(sigma)
        Lorenz2 = np.zeros_like(sigma)
    sigmaS = gsp*units.tau_unit*units.kb*units.eC*K1*iNV*itemp
    PF = sigma*Seebeck**2
    print_matrix('sigma matrix (S/m)', sigma)
    print_matrix('kappa matrix (K22 only) (W/m/K)', kappa)
    print_matrix('kappa matrix (full) (W/m/K)', kappa2)
    print_matrix('sigmaS matrix (A/m/K)', sigmaS)
    print_matrix('Seebeck matrix (V/K)', Seebeck)
    print_matrix('Peltier matrix (V)', Peltier)
    print_matrix('Lorenz matrix (K22 only) (Wohm/K^2) (fe 2.44e-8)', Lorenz, width=9, prec=2)
    print_matrix('Lorenz matrix? (full) (Wohm/K^2)', Lorenz2, width=9, prec=2)
    print_matrix('Power Factor (SA/m^2/K)', PF, ndigits=8)
    out = {'mu': mu, 'sigma': sigma, 'kappa': kappa, 'kappa2': kappa2, 'sigmaS': sigmaS, 'Seebeck': Seebeck,
           'Peltier': Peltier, 'Lorenz': Lorenz, 'Lorenz2': Lorenz2, 'PF': PF}
    if sw_tdf:
        out['tdf'], out['tdf_e'] = tdf, wlist-mu
    return out


def calc_conductivity_lrt(rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray, avec: np.ndarray, Nx: int, Ny: int,
                          Nz: int, fill: float, temp: float, Nw: int, delta: float, units: _units.Units, Emax: float,
                          color_thick: float = 1.0e-3, color_scatt: float = 1.0e-6,
                          with_spin: bool = False) -> dict | None:
    '''
    calculation of linear response theory
    electric conductivity of LRT corresponds to Boltzmann then delta~O(10-1) (tau~1fs) at 300K
    thermal conductivity and L12 of LRT correspond to Boltzmann then delta~O(10-3) (tau~100fs) at 300K
    Emax: upper end of the frequency window [eV]; color_thick / color_scatt: slab thickness and
    Kubelka-Munk scattering length [m] of the non-metal colours
    Returns a dict for the figures (libs.plots.plot_optical): wlist, sigma [Nw,3,3] (S/m),
    epsilon, refl [Nw,3], wp_scr (screened plasma frequency or None), and the perceived
    colours cols / ncols.  None if the input is invalid.
    '''
    # Parameter validation
    if temp <= 0:
        print("Error: Temperature (temp) is non-positive", flush=True)
        return
    if delta <= 0:
        print("Error: Broadening (delta) is non-positive", flush=True)
        return

    Nk, eig, vk, kweight = _bands.get_emesh(Nx, Ny, Nz, ham_r, S_r, rvec, avec.T*units.ihbar, True, True)
    Vuc = sclin.det(avec)
    if Vuc <= 0:
        print("Error: Unit cell volume (Vuc) is non-positive", flush=True)
        return
    gsp = (1.0 if with_spin else 2.0)  # spin weight
    mu = _bands.calc_mu(eig, Nk, fill, temp, kweight)
    print(f'chemical potential = {mu:.4f} eV', flush=True)
    print(f'temperature = {temp/units.kb:.3f} K', flush=True)
    print(f'about tau = {1.e15*units.hbar/delta:.3f} fs', flush=True)
    print(f'delta = {delta:9.3e}', flush=True)
    L11, L12, L22, wlist = _response.get_conductivity(mu, temp, eig, vk, Nw, Emax, delta, kweight)
    sigmaconst = gsp*units.hbar*units.eC/Vuc*1.0e30
    kappaSconst = sigmaconst*units.kb/temp
    sigma = sigmaconst*L11
    kappa = kappaSconst*L22
    sigmaS = kappaSconst*L12
    # full electronic thermal conductivity subtracts the thermoelectric backflow:
    # kappa = kappaSconst*(L22 - L12.L11^-1.L12)  (cf. Boltzmann K2-K1.K0^-1.K1).
    # The subtraction is done on the bare L tensors (all share the same prefactor),
    # so kappa2 carries the same kappaSconst as the L22-only kappa.
    try:
        kappa2 = kappaSconst*np.array([l22-l12.dot(sclin.inv(l11).dot(l12))
                                       for l11, l12, l22 in zip(L11, L12, L22)])
    except np.linalg.LinAlgError:
        print("Warning: L11 is a singular matrix. Cannot subtract thermoelectric backflow from kappa", flush=True)
        kappa2 = kappa.copy()
    # Handle sclin.inv() failures
    try:
        Seebeck = np.array([-sclin.inv(s).dot(sS) for s, sS in zip(sigma, sigmaS)])
    except np.linalg.LinAlgError:
        print("Warning: sigma is a singular matrix. Cannot compute Seebeck coefficient", flush=True)
        Seebeck = np.zeros_like(sigma)
    print_matrix('sigma matrix (S/m)', sigma[0].real)
    print_matrix('kappa matrix (L22 only) (W/m/K)', kappa[0].real)
    print_matrix('kappa matrix (full) (W/m/K)', kappa2[0].real)
    print_matrix('sigmaS matrix (A/m/K)', sigmaS[0].real)
    print('Lorenz number (Wohm/K^2) (fe 2.44e-8)', flush=True)
    try:
        # Lorenz tensor L = kb * kappa . (sigma*T)^-1 (matrix product, not element-wise)
        sigmaT_inv = sclin.inv(sigma[0]*temp)
        Lorenz_matrix = (units.kb*kappa[0].dot(sigmaT_inv)).real.round(10)
        Lorenz_matrix2 = (units.kb*kappa2[0].dot(sigmaT_inv)).real.round(10)
        print(' (L22 only):', flush=True)
        for lor in Lorenz_matrix:
            print(f" {lor[0]:9.2e} {lor[1]:9.2e} {lor[2]:9.2e}", flush=True)
        print(' (full):', flush=True)
        for lor in Lorenz_matrix2:
            print(f" {lor[0]:9.2e} {lor[1]:9.2e} {lor[2]:9.2e}", flush=True)
    except np.linalg.LinAlgError:
        print("Warning: Failed to compute Lorenz coefficient (singular matrix)", flush=True)
    print_matrix('Seebeck coefficient matrix (V/K)', Seebeck[0].real)
    # dielectric function (relative permittivity) from the optical conductivity,
    # SI convention: eps(w) = 1 + i*sigma(w)/(eps0*w_rad) = 1 + i*hbar*sigma/(eps0*w),
    # sigma in S/m, w in eV (hbar converts the eV grid to angular frequency).
    # Frequency-resolved like sigma; the w=0 point is skipped (eps ~ 1/w there).
    eps0 = scconst.epsilon_0
    epsilon = 1.+1j*units.hbar/eps0*(sigma.T/(wlist+1e-8)).T
    # console summary at the lowest positive frequency (static-limit estimate)
    print(f'dielectric function at w = {wlist[1]:9.3e} eV', flush=True)
    print_matrix('Re[eps] matrix', epsilon[1].real)
    print_matrix('Im[eps] matrix', epsilon[1].imag)
    # --- optical constants: N, R, alpha, and the non-metal spectra -----------
    # N(w) = n + i*kappa = sqrt(eps(w)) on the passive branch (Im N >= 0), then
    # the normal-incidence Fresnel reflectivity R = |(N-1)/(N+1)|^2 and the
    # Beer-Lambert absorption coefficient alpha = 2*w*kappa/(hbar*c). Only the
    # principal (diagonal) components are used, so these describe light
    # polarized along a crystal axis; w=0 is skipped as for eps.
    Nref = _optical.refractive_index(epsilon.diagonal(axis1=1, axis2=2))  # [Nw,3] xx,yy,zz
    refl = _optical.reflectivity(Nref)
    alpha = _optical.absorption_coefficient(wlist, Nref)                  # [1/m]
    # what leaves a NON-metal: the slab transmittance at thickness color_thick
    # and the Kubelka-Munk diffuse (body) reflectance at scattering length
    # color_scatt. A metal is opaque, so both come out black -- correctly so.
    trans = _optical.transmittance(refl, alpha, color_thick)
    rinf = _optical.kubelka_munk(alpha, color_scatt)
    absorp = np.clip(1.0-refl-trans, 0.0, 1.0)          # absorptance of the slab
    print(f'reflectivity (normal incidence) at w = {wlist[1]:9.3e} eV', flush=True)
    print(f'  R = ({refl[1, 0]:.4f}, {refl[1, 1]:.4f}, {refl[1, 2]:.4f})  (xx, yy, zz)', flush=True)
    print(f'absorption coefficient at w = {wlist[1]:9.3e} eV', flush=True)
    print(f'  alpha = ({alpha[1, 0]:9.3e}, {alpha[1, 1]:9.3e}, {alpha[1, 2]:9.3e}) 1/m', flush=True)
    # --- perceived colour (metallic lustre) from R(w) ------------------------
    # R(w) is folded with the CIE 1931 2-deg observer under illuminant D65 and
    # converted to sRGB (see plibs._color), the standard route for first-
    # principles metal colours (Prandini et al., npj Comput. Mater. 5, 129
    # (2019)). The polycrystalline entry uses the directional average of R,
    # which is what an unpolarised measurement on a polycrystal sees.
    cols = [(lab, _color.spectrum_color(wlist[1:], r))
            for lab, r in [('xx', refl[1:, 0]), ('yy', refl[1:, 1]), ('zz', refl[1:, 2]),
                           ('poly', refl[1:].mean(axis=1))]]
    # --- non-metal colours: the light that got through the material ----------
    # A dielectric reflects almost nothing, so its colour is carried by the
    # transmitted (or diffusely back-scattered) light, both computed from the
    # directional average. The absorptance entry is the complement of what
    # leaves the sample, so it is reported as its photographic negative.
    ncols = [('T(d)', _color.spectrum_color(wlist[1:], trans[1:].mean(axis=1)), 'hex'),
             ('KM',  _color.spectrum_color(wlist[1:], rinf[1:].mean(axis=1)), 'hex'),
             ('A neg', _color.spectrum_color(wlist[1:], absorp[1:].mean(axis=1)), 'hex_neg')]
    print('perceived colour (CIE 1931 2-deg observer, D65, sRGB):', flush=True)
    # show the colour itself as a 24-bit ANSI block, but only on a real terminal
    # so redirected logs stay free of escape sequences
    tty = os.isatty(1)
    for lab, c in cols:
        r, g, b = (int(v) for v in c['rgb255'])
        sq = f'\x1b[48;2;{r};{g};{b}m    \x1b[0m ' if tty else ''
        print(f'  {lab:4s} {sq}{c["hex"]}  RGB = ({r:3d},{g:3d},{b:3d})'
              f'  xy = ({c["xy"][0]:.4f}, {c["xy"][1]:.4f})  Y = {c["Y"]:.4f}'
              + ('  [out of sRGB gamut, clipped]' if c['clipped'] else ''), flush=True)
    print(f'non-metal colour (d = {1e3*color_thick:.4g} mm, '
          f'l_scatt = {1e6*color_scatt:.4g} um, directional average):', flush=True)
    for lab, c, key in ncols:
        r, g, b = (int(v) for v in (c['rgb255_neg'] if key == 'hex_neg' else c['rgb255']))
        sq = f'\x1b[48;2;{r};{g};{b}m    \x1b[0m ' if tty else ''
        note = (' (opaque at this thickness)' if lab == 'T(d)' and c['Y'] < 1.e-4 else
                ' (negated: A is the complement of the transmitted light)'
                if key == 'hex_neg' else '')
        # Y belongs to the spectrum that was integrated, so on the negated row
        # it is labelled Y(A): the luminance of the absorptance, not of the swatch
        print(f'  {lab:5s} {sq}{c[key]}  RGB = ({r:3d},{g:3d},{b:3d})'
              f'  {"Y(A)" if key == "hex_neg" else "Y"} = {c["Y"]:.4f}{note}', flush=True)
    if not cols[0][1]['full_range']:
        print(f'  Warning: the spectrum does not cover the visible window '
              f'(380-780 nm = 1.59-3.26 eV); R was held constant outside '
              f'[{_color.ev_nm/wlist[-1]:.0f}, {_color.ev_nm/wlist[1]:.0f}] nm, so the colour is '
              f'only indicative. Increase Emax (and Nw) to fix this.', flush=True)
    # --- plasma frequency ----------------------------------------------------
    # (1) effective (unscreened) plasma frequency from the optical f-sum rule
    #     (hbar*wp_a)^2 = (2*hbar/(pi*eps0)) * int_0^Emax Re sigma_aa(w) dw
    #     (sigma in S/m, w in eV -> hbar*wp in eV). It is the spectral weight up
    #     to the cutoff Emax and, unlike a Drude fit, needs no broadening input.
    hwp = np.sqrt([abs(2*units.hbar/(np.pi*eps0)*np.trapezoid(sigma[:, i, i].real, wlist)) for i in range(3)])
    h_eVs = scconst.physical_constants['Planck constant in eV s'][0]  # h = 2*pi*hbar
    fp = hwp/h_eVs/1e12   # plasma frequency nu_p = wp/(2*pi) in THz
    print(f'plasma frequency (f-sum rule, integrated to Emax={Emax} eV):', flush=True)
    print(f'  hbar*wp = ({hwp[0]:.4f}, {hwp[1]:.4f}, {hwp[2]:.4f}) eV  (xx, yy, zz)', flush=True)
    print(f'  nu_p    = ({fp[0]:.2f}, {fp[1]:.2f}, {fp[2]:.2f}) THz', flush=True)
    # (2) screened plasma frequency: energy where Re eps_xx crosses zero upward
    #     (the longitudinal plasmon / reflectivity edge), by linear interpolation.
    e1, wpos = epsilon[1:, 0, 0].real, wlist[1:]
    zc = np.where((e1[:-1] < 0.0)&(e1[1:] >= 0.0))[0]
    if zc.size:
        j = zc[0]
        wp_scr = wpos[j]-e1[j]*(wpos[j+1]-wpos[j])/(e1[j+1]-e1[j])
        print(f'  screened plasma freq (Re eps_xx=0): hbar*wp_scr = {wp_scr:.4f} eV'
              f'  (nu_scr = {wp_scr/h_eVs/1e12:.2f} THz)', flush=True)
    else:
        wp_scr = None
        print('  screened plasma freq (Re eps_xx=0): no upward zero crossing in [0, Emax]', flush=True)
    # write sigma(w) and eps(w) (diagonal components) to a data file
    try:
        header = ('perceived colour (CIE1931 2-deg, D65, sRGB) '
                  + '  '.join(f'{lab}: {c["hex"]} RGB({c["rgb255"][0]},{c["rgb255"][1]},'
                              f'{c["rgb255"][2]})' for lab, c in cols)+'\n'
                  + f'non-metal colour (d = {color_thick:.4g} m, l_scatt = {color_scatt:.4g} m) '
                  + '  '.join(f'{lab}: {c[key]}' for lab, c, key in ncols)+'\n'
                  + f'plasma frequency (f-sum rule) hbar*wp [eV] = '
                  f'{hwp[0]:.6f} {hwp[1]:.6f} {hwp[2]:.6f}  (xx yy zz)\n'
                  + (f'screened plasma frequency (Re eps_xx=0) hbar*wp_scr [eV] = {wp_scr:.6f}\n'
                     if wp_scr is not None else 'screened plasma frequency (Re eps_xx=0): none in window\n')
                  + 'w[eV] '
                  'Re_sigma_xx Im_sigma_xx Re_sigma_yy Im_sigma_yy Re_sigma_zz Im_sigma_zz '
                  'Re_eps_xx Im_eps_xx Re_eps_yy Im_eps_yy Re_eps_zz Im_eps_zz '
                  'R_xx R_yy R_zz alpha_xx alpha_yy alpha_zz[1/m] '
                  'T_xx T_yy T_zz Rinf_xx Rinf_yy Rinf_zz')
        d = np.arange(3)
        data = np.column_stack([wlist[1:],
                                *[f(sigma[1:, i, i]) for i in d for f in (np.real, np.imag)],
                                *[f(epsilon[1:, i, i]) for i in d for f in (np.real, np.imag)],
                                refl[1:], alpha[1:], trans[1:], rinf[1:]])
        np.savetxt('optical.dat', data, header=header)
        print("optical conductivity, dielectric function and reflectivity written to 'optical.dat'", flush=True)
    except OSError as e:
        print(f"Error: Failed to write 'optical.dat': {e}", flush=True)
    return {'wlist': wlist, 'sigma': sigma, 'epsilon': epsilon, 'refl': refl, 'wp_scr': wp_scr,
            'cols': cols, 'ncols': ncols}
