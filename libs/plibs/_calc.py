#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Higher-level calculation routines: FLEX, Eliashberg, susceptibility spectra, carrier analysis.
"""
import numpy as np
import libs.flibs as flibs
from ._bands import get_eigs, get_emesh, calc_mu
from ._eilenberger import BCS_RATIO
from ._response import get_initial_gap, chis_spectrum, phi_spectrum
from ._wannier_io import output_self_wannier, output_gap_wannier


def _prepare_kspace_state(Nx:int, Ny:int, Nz:int, ham_r, S_r, rvec, need_ham:bool=False):
    klist, kmap, invk = flibs.gen_irr_k_TRS(Nx, Ny, Nz)
    eig, uni = get_eigs(klist, ham_r, S_r, rvec)
    state = {
        'klist': klist,
        'kmap': kmap,
        'invk': invk,
        'eig': eig,
        'uni': uni,
    }
    if need_ham:
        state['ham_k'] = flibs.gen_ham(klist, ham_r, rvec)
    return state


def _prepare_normal_interaction(olist, site, orb_dep:bool, U:float, J:float, Umat=None, Jmat=None):
    if orb_dep:
        return flibs.gen_SCmatrix_orb(olist, site, Umat, Jmat)
    return flibs.gen_SCmatrix(olist, site, U, J)


def _prepare_soc_interaction(olist, slist, site, invs, orb_dep:bool, U:float, J:float, Umat=None, Jmat=None):
    if orb_dep:
        return flibs.gen_Vmatrix_orb(olist, slist, site, invs, Umat, Jmat)
    return flibs.gen_Vmatrix(olist, slist, site, invs, U, J)


def _load_sigma_from_file():
    """Load (sigmak, mu_self) from self_en.npz, or None if the file is missing."""
    try:
        npz = np.load('self_en.npz')
        print("Import sigma from self_en.npz")
        return npz['arr_0'], npz['arr_1']
    except FileNotFoundError:
        print("Error: 'self_en.npz' not found", flush=True)
        return None


def _prepare_green_state_normal(state, olist, interaction, mu:float, fill:float, temp:float,
                                Nw:int, Nx:int, Ny:int, Nz:int, sw_self:bool,
                                sw_from_file:bool=False, sw_out_self:bool=False, sw_in_self:bool=False,
                                eps:float=1.0e-4, pp:float=0.5, m_diis:int=5, sw_rescale:bool=True,
                                sw_tail:bool=False):
    Smat, Cmat = interaction
    if sw_self:
        if sw_from_file:
            loaded = _load_sigma_from_file()
            if loaded is None:
                return None
            sigmak, mu_self = loaded
        else:
            sigmak, mu_self = flibs.mkself(Smat, Cmat, state['kmap'], state['invk'], olist,
                                           state['ham_k'], state['eig'], state['uni'],
                                           mu, fill, temp, Nw, Nx, Ny, Nz, sw_out_self, sw_in_self,
                                           eps=eps, pp=pp, m_diis=m_diis, sw_rescale=sw_rescale,
                                           sw_tail=sw_tail)
        print(f'chem. pot. with self= {mu_self:.4f} eV', flush=True)
        Gk = flibs.gen_green(sigmak, state['ham_k'], mu_self, temp)
    else:
        sigmak = None
        mu_self = mu
        Gk = flibs.gen_Green0(state['eig'], state['uni'], mu, temp, Nw)
    return {'Gk': Gk, 'mu_self': mu_self, 'sigmak': sigmak}


def _prepare_green_state_soc(state, olist, slist, invs, interaction, mu:float, fill:float, temp:float,
                             Nw:int, Nx:int, Ny:int, Nz:int, sw_self:bool,
                             sw_from_file:bool=False, sw_out_self:bool=False, sw_in_self:bool=False,
                             eps:float=1.0e-4, pp:float=0.5, m_diis:int=5, sw_rescale:bool=True):
    Vmat = interaction
    if sw_self:
        if sw_from_file:
            loaded = _load_sigma_from_file()
            if loaded is None:
                return None
            sigmak, mu_self = loaded
        else:
            sigmak, mu_self = flibs.mkself_soc(Vmat, state['kmap'], state['invk'], invs, olist, slist,
                                               state['ham_k'], state['eig'], state['uni'],
                                               mu, fill, temp, Nw, Nx, Ny, Nz, sw_out_self, sw_in_self,
                                               eps=eps, pp=pp, m_diis=m_diis, sw_rescale=sw_rescale)
        print(f'chem. pot. with self= {mu_self:.4f} eV', flush=True)
        Gk = flibs.gen_green(sigmak, state['ham_k'], mu_self, temp)
    else:
        sigmak = None
        mu_self = mu
        Gk = flibs.gen_Green0(state['eig'], state['uni'], mu, temp, Nw)
    return {'Gk': Gk, 'mu_self': mu_self, 'sigmak': sigmak}

def _get_chi0_normal(state, Gk, chiolist, Smat, Cmat, mu_self, temp, Nx, Ny, Nz, sw_tail):
    """Irreducible chi0 on the irreducible q-grid; sw_tail=True uses the
    tail-corrected convolution (conv[G]-conv[G0]+analytic reference, needs the
    chemical potential mu_self of Gk)."""
    if sw_tail:
        return flibs.get_chi0_tail(Smat, Cmat, Gk, state['eig'], state['uni'], chiolist,
                                   state['kmap'], state['invk'], mu_self, temp, Nx, Ny, Nz)
    return flibs.get_chi0(Smat, Cmat, Gk, chiolist, state['kmap'], state['invk'], temp, Nx, Ny, Nz)


def calc_path_spectrum(kind:str, mu:float, temp:float, klist, qlist, chiolist, eig, uni,
                       spa_length, Nw:int, Emax:float, delta:float, Smat=None):
    if kind == 'chis':
        print("calculate spin susceptibility", flush=True)
        spec, spec_orb, wlist = chis_spectrum(mu, temp, Smat, klist, qlist, chiolist, eig, uni, Nw, Emax, delta)
        outname = 'chis.dat'
    elif kind == 'phi':
        print("calculate sc susceptibility", flush=True)
        spec, spec_orb, wlist = phi_spectrum(mu, temp, klist, qlist, chiolist, eig, uni, Nw, Emax, delta)
        outname = 'phi.dat'
    else:
        raise ValueError(f"Unknown spectrum kind: {kind}")

    w, sp = np.meshgrid(wlist, spa_length)
    try:
        with open(outname, 'w') as f:
            for ww, ssp, spec_row in zip(w, sp, spec):
                for www, sssp, spec_val in zip(ww, ssp, spec_row):
                    f.write(f'{sssp:8.4f} {www:8.4f} {spec_val.imag:9.4f}\n')
                f.write('\n')
    except IOError as e:
        print(f"Error: Failed to write '{outname}': {e}", flush=True)

    if kind == 'chis':
        for i, spec_orb_row in enumerate(spec_orb.T):
            try:
                with open(f'chis_{i}.dat', 'w') as f:
                    for ww, ssp, spec_band in zip(w, sp, spec_orb_row.T):
                        for www, sssp, spec_val in zip(ww, ssp, spec_band):
                            f.write(f'{sssp:8.4f} {www:8.4f} {spec_val.imag:9.4f}\n')
                        f.write('\n')
            except IOError as e:
                print(f"Error: Failed to write 'chis_{i}.dat': {e}", flush=True)
                continue

    return w, sp, spec

def get_carrier_num(kmesh, rvec, ham_r, S_r, mu:float, Arot):
    Nk, eig, kweight = get_emesh(kmesh, kmesh, kmesh, ham_r, S_r, rvec, Arot)
    if Nk <= 0:
        print("Error: Number of k-points (Nk) is non-positive", flush=True)
        return
    fill = 0.0
    for i, en in enumerate(eig.T - mu):
        num_hole = float(np.where(en > 0)[0].size) / Nk
        num_particle = float(np.where(en <= 0)[0].size) / Nk
        print(i+1, round(num_hole, 4), round(num_particle, 4), flush=True)
        fill += num_particle
    print(f'sum of electrons is {fill:.4f}', flush=True)

def get_mu(ham_r, S_r, rvec, Arot, temp:float, fill:float, kmesh=40) -> float:
    if temp <= 0:
        print("Error: Temperature (temp) is non-positive", flush=True)
        return None
    if kmesh <= 0:
        print("Error: k-mesh size (kmesh) is non-positive", flush=True)
        return None
    print("calc chem. pot.", flush=True)
    print(f"band filling = {fill:.4f}", flush=True)
    Nk, eig, kweight = get_emesh(kmesh, kmesh, kmesh, ham_r, S_r, rvec, Arot)
    mu = calc_mu(eig, Nk, fill, temp)
    return mu

def output_gap_function(invk, kmap, gap, uni, plist, gap_sym, Nx:int,
                        soc=False, invs=None, slist=None, sw_orb=False):
    if sw_orb:
        if soc:
            gapb = gap[:, :, 0, :]
        else:
            gapb = flibs.remap_gap(gap[:, :, 0, :], plist, invk, gap_sym)
    else:
        if soc:
            gapb = flibs.conv_delta_orb_to_band_soc(gap, uni, invk, invs, slist)
        else:
            gapb = flibs.conv_delta_orb_to_band(gap, uni, invk, plist, gap_sym)
    print('output gap function')
    for iorb in range(len(gapb)):
        for jorb in range(len(gapb)):
            try:
                with open(f'gap_{iorb+1}{jorb+1}.dat', 'w') as f:
                    for i, km in enumerate(kmap):
                        if km[2] == 0:
                            f.write(f'{km[0]:3} {km[1]:3} {gapb[iorb,jorb,i].real:15.8e} {gapb[iorb,jorb,i].imag:15.8e}\n')
                            if km[0] == Nx - 1:
                                f.write('\n')
            except IOError as e:
                print(f"Error: Failed to write 'gap_{iorb+1}{jorb+1}.dat': {e}", flush=True)
                continue
    return 0

def calc_flex(Nx:int, Ny:int, Nz:int, Nw:int, ham_r, S_r, rvec, mu:float, temp:float,
              olist, site, orb_dep:bool, U:float, J:float, fill:float, sw_out_self:bool, sw_in_self:bool, 
              Umat=None, Jmat=None, eps=1.0e-4, pp=0.5, m_diis=5, sw_rescale:bool=True, sw_tail:bool=False):
    state = _prepare_kspace_state(Nx, Ny, Nz, ham_r, S_r, rvec, need_ham=True)
    # FLEX uses the irreducible k-mesh Green's function together with the same orbital-pair
    # interaction basis used in the response routines above.
    Smat, Cmat = _prepare_normal_interaction(olist, site, orb_dep, U, J, Umat, Jmat)
    sigmak, mu_self = flibs.mkself(Smat, Cmat, state['kmap'], state['invk'], olist, state['ham_k'], state['eig'], state['uni'],
                                   mu, fill, temp, Nw, Nx, Ny, Nz, sw_out_self, sw_in_self,
                                   eps=eps, pp=pp, m_diis=m_diis, sw_rescale=sw_rescale, sw_tail=sw_tail)
    if sw_out_self:
        np.savez('self_en', sigmak, mu_self)
        output_self_wannier(sigmak, mu_self, state['kmap'], state['invk'], Nx, Ny, Nz, Nw, temp)

def calc_flex_soc(Nx:int, Ny:int, Nz:int, Nw:int, ham_r, S_r, rvec, mu:float, temp:float,
                  olist, slist, invs, site,
                  orb_dep:bool, U:float, J:float, fill:float,
                  sw_out_self:bool, sw_in_self:bool,
                  Umat=None, Jmat=None, eps=1.0e-4, pp=0.5, m_diis=5, sw_rescale:bool=True):
    state = _prepare_kspace_state(Nx, Ny, Nz, ham_r, S_r, rvec, need_ham=True)
    Vmat = _prepare_soc_interaction(olist, slist, site, invs, orb_dep, U, J, Umat, Jmat)
    sigmak, mu_self = flibs.mkself_soc(Vmat, state['kmap'], state['invk'], invs, olist, slist, state['ham_k'], state['eig'], state['uni'],
                                       mu, fill, temp, Nw, Nx, Ny, Nz, sw_out_self, sw_in_self,
                                       eps=eps, pp=pp, m_diis=m_diis, sw_rescale=sw_rescale)
    if sw_out_self:
        np.savez('self_en', sigmak, mu_self)
        output_self_wannier(sigmak, mu_self, state['kmap'], state['invk'], Nx, Ny, Nz, Nw, temp)

def calc_lin_eliashberg_eq(Nx:int, Ny:int, Nz:int, Nw:int, ham_r, S_r, rvec, chiolist, site, plist,
                           mu:float, temp:float, gap_sym:int, sw_self:bool, orb_dep:bool, U:float, J:float,
                           fill:float, sw_from_file:bool, sw_out_self:bool, sw_in_self:bool,
                           Umat=None, Jmat=None, eps=1.0e-4, pp=0.5, m_diis=5, sw_rescale:bool=True,
                           sw_tail:bool=False):
    state = _prepare_kspace_state(Nx, Ny, Nz, ham_r, S_r, rvec, need_ham=sw_self)
    Smat, Cmat = _prepare_normal_interaction(chiolist, site, orb_dep, U, J, Umat, Jmat)
    green_state = _prepare_green_state_normal(state, chiolist, (Smat, Cmat), mu, fill, temp, Nw,
                                              Nx, Ny, Nz, sw_self, sw_from_file, sw_out_self,
                                              sw_in_self, eps, pp, m_diis, sw_rescale, sw_tail)
    if green_state is None:
        return
    Gk = green_state['Gk']
    init_delta = get_initial_gap(state['klist'], len(state['eig'].T), gap_sym)
    sw_eig = True  # True: leading eigenvalue (most-unstable mode); False: trace
    chi, stoner = _get_chi0_normal(state, Gk, chiolist, Smat, Cmat, green_state['mu_self'],
                                   temp, Nx, Ny, Nz, sw_tail)
    chis, chic = flibs.get_chis_chic(chi, Smat, Cmat)
    chisq = flibs.get_eig_or_tr_chi(chis, state['invk'], sw_eig)
    chicq = flibs.get_eig_or_tr_chi(chic, state['invk'], sw_eig)
    try:
        with open('chis.dat', 'w') as f, open('chic.dat', 'w') as f2:
            for i, k in enumerate(state['kmap']):
                if k[2] == 0.0:
                    f.write(f'{k[0]:6.4f} {k[1]:6.4f} {chisq[i].real:11.4e}\n')
                    f2.write(f'{k[0]:6.4f} {k[1]:6.4f} {chicq[i].real:11.4e}\n')
    except IOError as e:
        print(f"Error: Failed to write 'chis.dat' or 'chic.dat': {e}", flush=True)
    gap, lambda_eliash = flibs.linearized_eliashberg(chi, Gk, state['uni'], init_delta, Smat, Cmat, chiolist, plist,
                                      state['kmap'], state['invk'], Nx, Ny, Nz, temp, gap_sym)
    print(f'Stoner factor = {stoner:.6f}, lambda_eliash = {lambda_eliash:.6f}', flush=True)
    if sw_out_self:
        np.save('gap', gap)
        output_gap_wannier(gap, state['kmap'], state['invk'], Nx, Ny, Nz, Nw, temp)
    output_gap_function(state['invk'], state['kmap'], gap, state['uni'], plist, gap_sym, Nx)

def calc_lin_eliash_soc(Nx:int, Ny:int, Nz:int, Nw:int, ham_r, S_r, rvec, mu:float, temp:float,
                        chiolist, slist, plist, invs, site, orb_dep:bool, U:float, J:float, fill:float,
                        gap_sym:int, sw_self:bool, sw_from_file:bool, sw_out_self:bool, sw_in_self:bool,
                        Umat=None, Jmat=None, eps=1.0e-4, pp=0.5, m_diis=5, sw_rescale:bool=True):
    state = _prepare_kspace_state(Nx, Ny, Nz, ham_r, S_r, rvec, need_ham=sw_self)
    Vmat = _prepare_soc_interaction(chiolist, slist, site, invs, orb_dep, U, J, Umat, Jmat)
    green_state = _prepare_green_state_soc(state, chiolist, slist, invs, Vmat, mu, fill, temp, Nw,
                                           Nx, Ny, Nz, sw_self, sw_from_file, sw_out_self,
                                           sw_in_self, eps, pp, m_diis, sw_rescale)
    if green_state is None:
        return
    Gk = green_state['Gk']
    sw_eig = True  # True: leading eigenvalue (most-unstable mode); False: trace
    chi, sgnsig, sgnsig2, invschi = flibs.get_chi0_soc(Vmat, Gk, chiolist, slist, state['kmap'], state['invk'], invs, temp, Nx, Ny, Nz)
    chic, chiszz, chispm = flibs.get_chis_chic_soc(chi, Vmat, chiolist, slist, invs)
    chiszzq = flibs.get_eig_or_tr_chi(chiszz, state['invk'], sw_eig)
    chispmq = flibs.get_eig_or_tr_chi(chispm, state['invk'], sw_eig)
    chicq = flibs.get_eig_or_tr_chi(chic, state['invk'], sw_eig)
    try:
        with open('chis.dat', 'w') as f, open('chic.dat', 'w') as f2:
            for i, k in enumerate(state['kmap']):
                if k[2] == 0.0:
                    f.write(f'{k[0]:6.4f} {k[1]:6.4f} {chiszzq[i].real:11.4e} {chispmq[i].real:11.4e}\n')
                    f2.write(f'{k[0]:6.4f} {k[1]:6.4f} {chicq[i].real:11.4e}\n')
    except IOError as e:
        print(f"Error: Failed to write 'chis.dat' or 'chic.dat': {e}", flush=True)
    init_delta = get_initial_gap(state['klist'], len(slist), gap_sym)
    gap = flibs.linearized_eliashberg_soc(chi, Gk, state['uni'], init_delta, Vmat, sgnsig, sgnsig2, plist, slist, chiolist,
                                          state['kmap'], state['invk'], invs, invschi, Nx, Ny, Nz, temp, gap_sym)
    if sw_out_self:
        np.save('gap', gap)
        output_gap_wannier(gap, state['kmap'], state['invk'], Nx, Ny, Nz, Nw, temp)
    output_gap_function(state['invk'], state['kmap'], gap, state['uni'], plist, gap_sym, Nx, True, invs, slist)

def calc_eliashberg_eq(Nx:int, Ny:int, Nz:int, Nw:int, ham_r, S_r, rvec,
                       chiolist, site, plist, mu:float, temp:float, gap_sym:int, sw_self:bool,
                       orb_dep:bool, U:float, J:float, fill:float, sw_from_file:bool, sw_out_self:bool,
                       sw_in_self:bool, Umat=None, Jmat=None, eps=1.0e-4, pp=0.5, m_diis=5, sw_rescale:bool=True,
                       sw_check_only:bool=False, sw_tail:bool=False):
    """
    @param sw_check_only: If True, stop after the linearized Eliashberg solve (before the
                          nonlinear loop) and report the Stoner factor and lambda_eliash.
                          The calculation also stops early (regardless of this flag) when
                          the Stoner factor >= 1 (SDW/CDW instability) or lambda_eliash < 1
                          (T >= Tc; no superconducting instability at this temperature).
    """
    state = _prepare_kspace_state(Nx, Ny, Nz, ham_r, S_r, rvec, need_ham=True)
    Smat, Cmat = _prepare_normal_interaction(chiolist, site, orb_dep, U, J, Umat, Jmat)
    green_state = _prepare_green_state_normal(state, chiolist, (Smat, Cmat), mu, fill, temp, Nw,
                                              Nx, Ny, Nz, sw_self, sw_from_file, sw_out_self,
                                              sw_in_self, eps, pp, m_diis, sw_rescale, sw_tail)
    if green_state is None:
        return
    Gk = green_state['Gk']
    mu_self = green_state['mu_self']
    delta_init_band = get_initial_gap(state['klist'], len(state['eig'].T), gap_sym)
    # Use linearized Eliashberg to obtain a symmetry-correct initial gap.
    # Alternative: flibs.get_initial_delta (plain orbital projection, kept for reference).
    if True:
        # Reuse the linearized problem as a cheap physical gate before entering the much more
        # expensive nonlinear loop. This filters out SDW/CDW instabilities and T >= Tc cases.
        chi, stoner = _get_chi0_normal(state, Gk, chiolist, Smat, Cmat, mu_self,
                                       temp, Nx, Ny, Nz, sw_tail)
        if stoner >= 1.0:
            print(f"Stoner factor = {stoner:.6f} >= 1: SDW/CDW instability, stopping before "
                  f"nonlinear Eliashberg.", flush=True)
            return
        delta_init, lambda_eliash = flibs.linearized_eliashberg(chi, Gk, state['uni'], delta_init_band, Smat, Cmat,
                                                 chiolist, plist, state['kmap'], state['invk'], Nx, Ny, Nz, temp, gap_sym)
        print(f'Stoner factor = {stoner:.6f}, lambda_eliash = {lambda_eliash:.6f}', flush=True)
        if lambda_eliash < 1.0:
            print(f"lambda_eliash = {lambda_eliash:.6f} < 1: no superconducting instability "
                  f"at this temperature (T >= Tc), stopping before nonlinear Eliashberg.", flush=True)
            return
        if sw_check_only:
            print("sw_check_only=True: stopping after linearized Eliashberg.", flush=True)
            return
    else:
        delta_init_band = get_initial_gap(state['klist'], len(state['eig'].T), gap_sym)
        delta_init = flibs.get_initial_delta(delta_init_band, state['uni'], state['kmap'], state['invk'], Nw, gap_sym)
    # BCS weak-coupling ratio: Δ₀ = 1.764 kB Tc; scale initial gap to this amplitude.
    target_gap = BCS_RATIO * temp
    delta_max = np.abs(delta_init).max()
    if delta_max > 0.0:
        # Only the shape from the linearized eigenvector matters here; the nonlinear solver
        # converges more reliably if the starting amplitude is normalized to the target scale.
        delta_init *= target_gap / delta_max
    else:
        print("Warning: initial gap is zero; skip Tc-based scaling", flush=True)
    delta, sigmak = flibs.nonlinear_eliashberg(delta_init, Gk, state['ham_k'], Smat, Cmat, chiolist, plist,
                                               state['kmap'], state['invk'], mu_self, temp, gap_sym, Nx, Ny, Nz,
                                               sw_sigma=sw_self, sw_Vconst=True, eps=eps, m_diis=m_diis)
    if sw_out_self:
        np.save('gap', delta)
        output_gap_wannier(delta, state['kmap'], state['invk'], Nx, Ny, Nz, Nw, temp)
    output_gap_function(state['invk'], state['kmap'], delta, state['uni'], plist, gap_sym, Nx)
