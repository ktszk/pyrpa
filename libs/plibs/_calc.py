#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Higher-level calculation routines: FLEX, Eliashberg, susceptibility spectra, carrier analysis.
"""
import numpy as np
import libs.flibs as flibs
from ._bands import get_eigs, get_emesh, calc_mu
from ._response import get_initial_gap, chis_spectrum, phi_spectrum
from ._wannier_io import output_self_wannier, output_gap_wannier

def calc_chis_spectrum(mu:float, temp:float, Smat, klist, qlist, chiolist, eig, uni, spa_length,
                       Nw:int, Emax:float, delta:float):
    print("calculate spin susceptibility", flush=True)
    chisw, chisw_orb, wlist = chis_spectrum(mu, temp, Smat, klist, qlist, chiolist, eig, uni, Nw, Emax, delta)
    w, sp = np.meshgrid(wlist, spa_length)
    try:
        with open('chis.dat', 'w') as f:
            for ww, ssp, chis in zip(w, sp, chisw):
                for www, sssp, chi in zip(ww, ssp, chis):
                    f.write(f'{sssp:8.4f} {www:8.4f} {chi.imag:9.4f}\n')
                f.write('\n')
    except IOError as e:
        print(f"Error: Failed to write 'chis.dat': {e}", flush=True)
    for i, chiso in enumerate(chisw_orb.T):
        try:
            with open(f'chis_{i}.dat', 'w') as f:
                for ww, ssp, chis in zip(w, sp, chiso.T):
                    for www, sssp, chi in zip(ww, ssp, chis):
                        f.write(f'{sssp:8.4f} {www:8.4f} {chi.imag:9.4f}\n')
                    f.write('\n')
        except IOError as e:
            print(f"Error: Failed to write 'chis_{i}.dat': {e}", flush=True)
            continue
    return w, sp, chisw

def calc_phi_spectrum(mu:float, temp:float, klist, qlist, chiolist, eig, uni, spa_length,
                      Nw:int, Emax:float, delta:float):
    print("calculate sc susceptibility", flush=True)
    phiw, phiw_orb, wlist = phi_spectrum(mu, temp, klist, qlist, chiolist, eig, uni, Nw, Emax, delta)
    w, sp = np.meshgrid(wlist, spa_length)
    try:
        with open('phi.dat', 'w') as f:
            for ww, ssp, phi_val in zip(w, sp, phiw):
                for www, sssp, ph in zip(ww, ssp, phi_val):
                    f.write(f'{sssp:8.4f} {www:8.4f} {ph.imag:9.4f}\n')
                f.write('\n')
    except IOError as e:
        print(f"Error: Failed to write 'phi.dat': {e}", flush=True)
    return w, sp, phiw

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
              Umat=None, Jmat=None, eps=1.0e-4, pp=0.5, m_diis=5, sw_rescale:bool=True):
    klist, kmap, invk = flibs.gen_irr_k_TRS(Nx, Ny, Nz)
    eig, uni = get_eigs(klist, ham_r, S_r, rvec)
    ham_k = flibs.gen_ham(klist, ham_r, rvec)
    if orb_dep:
        Smat, Cmat = flibs.gen_SCmatrix_orb(olist, site, Umat, Jmat)
    else:
        Smat, Cmat = flibs.gen_SCmatrix(olist, site, U, J)
    sigmak, mu_self = flibs.mkself(Smat, Cmat, kmap, invk, olist, ham_k, eig, uni,
                                   mu, fill, temp, Nw, Nx, Ny, Nz, sw_out_self, sw_in_self,
                                   eps=eps, pp=pp, m_diis=m_diis, sw_rescale=sw_rescale)
    if sw_out_self:
        np.savez('self_en', sigmak, mu_self)
        output_self_wannier(sigmak, mu_self, kmap, invk, Nx, Ny, Nz, Nw, temp)

def calc_flex_soc(Nx:int, Ny:int, Nz:int, Nw:int, ham_r, S_r, rvec, mu:float, temp:float,
                  olist, slist, invs, site,
                  orb_dep:bool, U:float, J:float, fill:float,
                  sw_out_self:bool, sw_in_self:bool,
                  Umat=None, Jmat=None, eps=1.0e-4, pp=0.5, m_diis=5, sw_rescale:bool=True):
    klist, kmap, invk = flibs.gen_irr_k_TRS(Nx, Ny, Nz)
    eig, uni = get_eigs(klist, ham_r, S_r, rvec)
    ham_k = flibs.gen_ham(klist, ham_r, rvec)
    if orb_dep:
        Vmat = flibs.gen_Vmatrix_orb(olist, slist, site, invs, Umat, Jmat)
    else:
        Vmat = flibs.gen_Vmatrix(olist, slist, site, invs, U, J)
    sigmak, mu_self = flibs.mkself_soc(Vmat, kmap, invk, invs, olist, slist, ham_k, eig, uni,
                                       mu, fill, temp, Nw, Nx, Ny, Nz, sw_out_self, sw_in_self,
                                       eps=eps, pp=pp, m_diis=m_diis, sw_rescale=sw_rescale)
    if sw_out_self:
        np.savez('self_en', sigmak, mu_self)
        output_self_wannier(sigmak, mu_self, kmap, invk, Nx, Ny, Nz, Nw, temp)

def calc_lin_eliashberg_eq(Nx:int, Ny:int, Nz:int, Nw:int, ham_r, S_r, rvec, chiolist, site, plist,
                           mu:float, temp:float, gap_sym:int, sw_self:bool, orb_dep:bool, U:float, J:float,
                           fill:float, sw_from_file:bool, sw_out_self:bool, sw_in_self:bool,
                           Umat=None, Jmat=None, eps=1.0e-4, pp=0.5, m_diis=5, sw_rescale:bool=True):
    klist, kmap, invk = flibs.gen_irr_k_TRS(Nx, Ny, Nz)
    eig, uni = get_eigs(klist, ham_r, S_r, rvec)
    if orb_dep:
        Smat, Cmat = flibs.gen_SCmatrix_orb(chiolist, site, Umat, Jmat)
    else:
        Smat, Cmat = flibs.gen_SCmatrix(chiolist, site, U, J)
    if sw_self:
        ham_k = flibs.gen_ham(klist, ham_r, rvec)
        if sw_from_file:
            try:
                npz = np.load('self_en.npz')
                sigmak, mu_self = npz['arr_0'], npz['arr_1']
                print("Import sigma from self_en.npz")
            except FileNotFoundError:
                print("Error: 'self_en.npz' not found", flush=True)
                return
        else:
            sigmak, mu_self = flibs.mkself(Smat, Cmat, kmap, invk, chiolist, ham_k, eig, uni,
                                           mu, fill, temp, Nw, Nx, Ny, Nz, sw_out_self, sw_in_self,
                                           eps=eps, pp=pp, m_diis=m_diis, sw_rescale=sw_rescale)
        print(f'chem. pot. with self= {mu:.4f} eV', flush=True)
        Gk = flibs.gen_green(sigmak, ham_k, mu_self, temp)
    else:
        Gk = flibs.gen_Green0(eig, uni, mu, temp, Nw)
    init_delta = get_initial_gap(klist, len(eig.T), gap_sym)
    sw_eig = True  # True: leading eigenvalue (most-unstable mode); False: trace
    chi = flibs.get_chi0(Smat, Cmat, Gk, chiolist, kmap, invk, temp, Nx, Ny, Nz)
    chis, chic = flibs.get_chis_chic(chi, Smat, Cmat)
    chisq = flibs.get_eig_or_tr_chi(chis, invk, sw_eig)
    chicq = flibs.get_eig_or_tr_chi(chic, invk, sw_eig)
    try:
        with open('chis.dat', 'w') as f, open('chic.dat', 'w') as f2:
            for i, k in enumerate(kmap):
                if k[2] == 0.0:
                    f.write(f'{k[0]:6.4f} {k[1]:6.4f} {chisq[i].real:11.4e}\n')
                    f2.write(f'{k[0]:6.4f} {k[1]:6.4f} {chicq[i].real:11.4e}\n')
    except IOError as e:
        print(f"Error: Failed to write 'chis.dat' or 'chic.dat': {e}", flush=True)
    gap = flibs.linearized_eliashberg(chi, Gk, uni, init_delta, Smat, Cmat, chiolist, plist,
                                      kmap, invk, Nx, Ny, Nz, temp, gap_sym)
    if sw_out_self:
        np.save('gap', gap)
        output_gap_wannier(gap, kmap, invk, Nx, Ny, Nz, Nw, temp)
    output_gap_function(invk, kmap, gap, uni, plist, gap_sym, Nx)

def calc_lin_eliash_soc(Nx:int, Ny:int, Nz:int, Nw:int, ham_r, S_r, rvec, mu:float, temp:float,
                        chiolist, slist, plist, invs, site, orb_dep:bool, U:float, J:float, fill:float,
                        gap_sym:int, sw_self:bool, sw_from_file:bool, sw_out_self:bool, sw_in_self:bool,
                        Umat=None, Jmat=None, eps=1.0e-4, pp=0.5, m_diis=5, sw_rescale:bool=True):
    klist, kmap, invk = flibs.gen_irr_k_TRS(Nx, Ny, Nz)
    eig, uni = get_eigs(klist, ham_r, S_r, rvec)
    if orb_dep:
        Vmat = flibs.gen_Vmatrix_orb(chiolist, slist, site, invs, Umat, Jmat)
    else:
        Vmat = flibs.gen_Vmatrix(chiolist, slist, site, invs, U, J)
    if sw_self:
        ham_k = flibs.gen_ham(klist, ham_r, rvec)
        if sw_from_file:
            try:
                npz = np.load('self_en.npz')
                sigmak, mu_self = npz['arr_0'], npz['arr_1']
            except FileNotFoundError:
                print("Error: 'self_en.npz' not found", flush=True)
                return
        else:
            sigmak, mu_self = flibs.mkself_soc(Vmat, kmap, invk, invs, chiolist, slist, ham_k, eig, uni,
                                               mu, fill, temp, Nw, Nx, Ny, Nz, sw_out_self, sw_in_self,
                                               eps=eps, pp=pp, m_diis=m_diis, sw_rescale=sw_rescale)
        print(f'chem. pot. with self= {mu:.4f} eV', flush=True)
        Gk = flibs.gen_green(sigmak, ham_k, mu_self, temp)
    else:
        Gk = flibs.gen_Green0(eig, uni, mu, temp, Nw)
    sw_eig = True  # True: leading eigenvalue (most-unstable mode); False: trace
    chi, sgnsig, sgnsig2, invschi = flibs.get_chi0_soc(Vmat, Gk, chiolist, slist, kmap, invk, invs, temp, Nx, Ny, Nz)
    chic, chiszz, chispm = flibs.get_chis_chic_soc(chi, Vmat, chiolist, slist, invs)
    chiszzq = flibs.get_eig_or_tr_chi(chiszz, invk, sw_eig)
    chispmq = flibs.get_eig_or_tr_chi(chispm, invk, sw_eig)
    chicq = flibs.get_eig_or_tr_chi(chic, invk, sw_eig)
    try:
        with open('chis.dat', 'w') as f, open('chic.dat', 'w') as f2:
            for i, k in enumerate(kmap):
                if k[2] == 0.0:
                    f.write(f'{k[0]:6.4f} {k[1]:6.4f} {chiszzq[i].real:11.4e} {chispmq[i].real:11.4e}\n')
                    f2.write(f'{k[0]:6.4f} {k[1]:6.4f} {chicq[i].real:11.4e}\n')
    except IOError as e:
        print(f"Error: Failed to write 'chis.dat' or 'chic.dat': {e}", flush=True)
    init_delta = get_initial_gap(klist, len(slist), gap_sym)
    gap = flibs.linearized_eliashberg_soc(chi, Gk, uni, init_delta, Vmat, sgnsig, sgnsig2, plist, slist, chiolist,
                                          kmap, invk, invs, invschi, Nx, Ny, Nz, temp, gap_sym)
    if sw_out_self:
        np.save('gap', gap)
        output_gap_wannier(gap, kmap, invk, Nx, Ny, Nz, Nw, temp)
    output_gap_function(invk, kmap, gap, uni, plist, gap_sym, Nx, True, invs, slist)

def calc_eliashberg_eq(Nx:int, Ny:int, Nz:int, Nw:int, ham_r, S_r, rvec,
                       chiolist, site, plist, mu:float, temp:float, gap_sym:int, sw_self:bool,
                       orb_dep:bool, U:float, J:float, fill:float, sw_from_file:bool, sw_out_self:bool,
                       sw_in_self:bool, Umat=None, Jmat=None, eps=1.0e-4, pp=0.5, m_diis=5, sw_rescale:bool=True):
    klist, kmap, invk = flibs.gen_irr_k_TRS(Nx, Ny, Nz)
    ham_k = flibs.gen_ham(klist, ham_r, rvec)
    eig, uni = get_eigs(klist, ham_r, S_r, rvec)
    if orb_dep:
        Smat, Cmat = flibs.gen_SCmatrix_orb(chiolist, site, Umat, Jmat)
    else:
        Smat, Cmat = flibs.gen_SCmatrix(chiolist, site, U, J)
    if sw_self:
        if sw_from_file:
            try:
                npz = np.load('self_en.npz')
                sigmak, mu_self = npz['arr_0'], npz['arr_1']
                print("Import sigma from self_en.npz")
            except FileNotFoundError:
                print("Error: 'self_en.npz' not found", flush=True)
                return
        else:
            sigmak, mu_self = flibs.mkself(Smat, Cmat, kmap, invk, chiolist, ham_k, eig, uni,
                                           mu, fill, temp, Nw, Nx, Ny, Nz, sw_out_self, sw_in_self,
                                           eps=eps, pp=pp, m_diis=m_diis, sw_rescale=sw_rescale)
        print(f'chem. pot. with self= {mu_self:.4f} eV', flush=True)
        Gk = flibs.gen_green(sigmak, ham_k, mu_self, temp)
    else:
        mu_self = mu
        Gk = flibs.gen_Green0(eig, uni, mu, temp, Nw)
    delta_init_band = get_initial_gap(klist, len(eig.T), gap_sym)
    # Use linearized Eliashberg to obtain a symmetry-correct initial gap.
    # Alternative: flibs.get_initial_delta (plain orbital projection, kept for reference).
    if True:
        chi = flibs.get_chi0(Smat, Cmat, Gk, chiolist, kmap, invk, temp, Nx, Ny, Nz)
        delta_init = flibs.linearized_eliashberg(chi, Gk, uni, delta_init_band, Smat, Cmat, chiolist, plist,
                                                 kmap, invk, Nx, Ny, Nz, temp, gap_sym)
    else:
        delta_init_band = get_initial_gap(klist, len(eig.T), gap_sym)
        delta_init = flibs.get_initial_delta(delta_init_band, uni, kmap, invk, Nw, gap_sym)
    # BCS weak-coupling ratio: Δ₀ = 1.764 kB Tc; scale initial gap to this amplitude.
    target_gap = 1.764 * temp
    delta_max = np.abs(delta_init).max()
    if delta_max > 0.0:
        delta_init *= target_gap / delta_max
    else:
        print("Warning: initial gap is zero; skip Tc-based scaling", flush=True)
    delta, sigmak = flibs.nonlinear_eliashberg(delta_init, Gk, ham_k, Smat, Cmat, chiolist, plist,
                                               kmap, invk, mu_self, temp, gap_sym, Nx, Ny, Nz,
                                               sw_sigma=sw_self, sw_Vconst=True, eps=eps, m_diis=m_diis)
    if sw_out_self:
        np.save('gap', delta)
        output_gap_wannier(delta, kmap, invk, Nx, Ny, Nz, Nw, temp)
    output_gap_function(invk, kmap, delta, uni, plist, gap_sym, Nx)
