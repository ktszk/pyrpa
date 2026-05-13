#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Response functions: spin susceptibility, pairing susceptibility, conductivity, and gap symmetries.
"""
import numpy as np
import libs.flibs as flibs


def get_conductivity(mu: float, temp: float, eig: np.ndarray, vk: np.ndarray, Nw: int, Emax: float,
                     idelta: float = 1.e-3) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    @fn get_conductivity
    @brief Compute transport coefficients L11, L12, L22 as a function of frequency using the Kubo formula.
    @param      mu: Chemical potential in eV
    @param    temp: Temperature in eV
    @param     eig: Eigenvalues array [Nk, Nband]
    @param      vk: Group velocities array [Nk, Nband, 3]
    @param      Nw: Number of frequency points
    @param    Emax: Maximum frequency in eV
    @param  idelta: Broadening parameter (Lorentzian width) in eV
    @retval    L11: Charge conductivity spectrum array [Nw, 3, 3]
    @retval    L12: Thermoelectric coefficient spectrum array [Nw, 3, 3]
    @retval    L22: Thermal conductivity spectrum array [Nw, 3, 3]
    @retval  wlist: Frequency mesh array [Nw]
    """
    wlist=np.linspace(0,Emax,Nw)
    ffermi=flibs.get_ffermi(eig,mu,temp)
    L11,L12,L22=flibs.calc_Lij_wl(eig,vk,ffermi,mu,wlist,idelta,temp)
    return L11,L12,L22,wlist

def chis_spectrum(mu: float, temp: float, Smat: np.ndarray, klist: np.ndarray, qlist: np.ndarray,
                  olist: np.ndarray, eig: np.ndarray, uni: np.ndarray, Nw: int, Emax: float,
                  idelta: float = 1.e-3) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    @fn chis_spectrum
    @brief Compute the spin susceptibility chi_s spectrum along a q-path and write results to 'chi0.dat'.
    @param      mu: Chemical potential in eV
    @param    temp: Temperature in eV
    @param    Smat: Stoner interaction matrix [Norb, Norb]
    @param   klist: k-point list [Nk, 3]
    @param   qlist: q-point path list [Nq, 3]
    @param   olist: Orbital index list for susceptibility calculation
    @param     eig: Eigenvalues array [Nk, Nband]
    @param     uni: Eigenvectors array [Nk, Norb, Nband]
    @param      Nw: Number of frequency points
    @param    Emax: Maximum frequency in eV
    @param  idelta: Broadening parameter (Lorentzian width) in eV
    @retval   chisq: Spin susceptibility trace at each q-point [Nq, Nw]
    @retval chis_orbq: Orbital-resolved spin susceptibility [Nq, ...]
    @retval   wlist: Frequency mesh array [Nw]
    """
    ffermi=flibs.get_ffermi(eig,mu,temp)
    wlist=np.linspace(0,Emax,Nw)
    chisq=[]
    chis_orbq=[]
    with open('chi0.dat','w') as f, open('writeq.dat','w') as fq:
        for i,q in enumerate(qlist):
            fq.write(f'{i:d} {q[0]:5.3f} {q[1]:5.3f} {q[2]:5.3f}\n')
            fq.flush()
            qshift=flibs.get_qshift(klist,q)
            chi0=flibs.get_chi_irr(uni,eig,ffermi,qshift,olist,wlist,idelta,temp)
            chis=flibs.get_chis(chi0,Smat)
            trchis,trchi0,chis_orb=flibs.get_tr_chi(chis,chi0,olist)
            chisq.append(trchis)
            chis_orbq.append(chis_orb)
            for w,trchi in zip(wlist,trchi0):
                f.write(f'{i:8.4f} {w:8.4f} {trchi.imag:9.4f} {trchi.real:9.4f}\n')
            f.write('\n')
            f.flush()
    return np.array(chisq),np.array(chis_orbq),wlist

def chis_q_point(q: np.ndarray, eig: np.ndarray, uni: np.ndarray, Emax: float,
                 Nw: int, mu: float, temp: float, Smat: np.ndarray, klist: np.ndarray,
                 olist: np.ndarray, idelta: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    @fn chis_q_point
    @brief Compute the spin susceptibility chi_s at a single q-point across a frequency mesh.
    @param       q: Single q-vector in fractional coordinates [3]
    @param     eig: Eigenvalues array [Nk, Nband]
    @param     uni: Eigenvectors array [Nk, Norb, Nband]
    @param    Emax: Maximum frequency in eV
    @param      Nw: Number of frequency points
    @param      mu: Chemical potential in eV
    @param    temp: Temperature in eV
    @param    Smat: Stoner interaction matrix [Norb, Norb]
    @param   klist: k-point list [Nk, 3]
    @param   olist: Orbital index list for susceptibility calculation
    @param  idelta: Broadening parameter (Lorentzian width) in eV
    @retval  trchis: Trace of spin susceptibility chi_s [Nw]
    @retval chis_orb: Orbital-resolved spin susceptibility
    @retval   wlist: Frequency mesh array [Nw]
    """
    ffermi=flibs.get_ffermi(eig,mu,temp)
    wlist=np.linspace(0,Emax,Nw)
    qshift=flibs.get_qshift(klist,q)
    chi0=flibs.get_chi_irr(uni,eig,ffermi,qshift,olist,wlist,idelta,temp)
    chis=flibs.get_chis(chi0,Smat)
    trchis,trchi0,chis_orb=flibs.get_tr_chi(chis,chi0,olist)
    return trchis,chis_orb,wlist

def chis_q_point_sc(q: np.ndarray, hamk: np.ndarray, delta_k: np.ndarray, mu: float,
                    Emax: float, Nw: int, temp: float, Smat: np.ndarray,
                    klist: np.ndarray, olist: np.ndarray, idelta: float,
                    sw_spsym: bool = False) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    @fn chis_q_point_sc
    @brief Compute the spin susceptibility chi_s at a single q-point in the superconducting state.
    Builds the BdG Hamiltonian from hamk and the gap function delta_k, diagonalizes it,
    and calls get_chi_irr_sc to evaluate the irreducible susceptibility.
    @param        q: Single q-vector in fractional coordinates [3]
    @param     hamk: Normal-state k-space Hamiltonian [Nk, Norb, Norb] complex128
    @param  delta_k: Gap function (anomalous potential) on k-mesh [Nk, Norb, Norb] complex128
    @param       mu: Chemical potential in eV
    @param     Emax: Maximum frequency in eV
    @param       Nw: Number of frequency points
    @param     temp: Temperature in eV
    @param     Smat: Stoner interaction matrix [Nchi, Nchi] float64
    @param    klist: k-point list [Nk, 3] float64
    @param    olist: Orbital index pairs for susceptibility [Nchi, 2] int64 (chiolist)
    @param   idelta: Lorentzian broadening in eV
    @param sw_spsym: True for triplet (dz) symmetry, False for singlet
    @retval  trchis: Trace of spin susceptibility chi_s [Nw] complex128
    @retval chis_orb: Orbital-resolved spin susceptibility
    @retval   wlist: Frequency mesh [Nw] float64
    """
    Norb = hamk.shape[1]
    hamBdGk = flibs.mkBdGhamk(hamk - mu * np.eye(Norb), delta_k)
    eig_BdG, uni_BdG = flibs.get_eig(hamBdGk)
    # mu=0: chemical potential is already absorbed into hamBdGk (hamk - mu*I), so BdG
    # quasi-particle energies are measured from zero and the Fermi level is at 0.
    ffermi_BdG = flibs.get_ffermi(eig_BdG, 0., temp)
    wlist = np.linspace(0, Emax, Nw)
    qshift = flibs.get_qshift(klist, q)
    chi0 = flibs.get_chi_irr_sc(uni_BdG, eig_BdG, ffermi_BdG, qshift, olist, wlist, idelta, temp, sw_spsym)
    chis = flibs.get_chis(chi0, Smat)
    trchis, trchi0, chis_orb = flibs.get_tr_chi(chis, chi0, olist)
    return trchis, chis_orb, wlist

def chis_spectrum_sc(mu: float, temp: float, Smat: np.ndarray, hamk: np.ndarray,
                     delta_k: np.ndarray, klist: np.ndarray, qlist: np.ndarray,
                     olist: np.ndarray, Nw: int, Emax: float, idelta: float = 1.e-3,
                     sw_spsym: bool = False) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    @fn chis_spectrum_sc
    @brief Compute the spin susceptibility chi_s spectrum along a q-path in the superconducting state.
    The BdG Hamiltonian is built and diagonalized once before the q-loop.
    Results are also written to 'chi0_sc.dat'.
    @param      mu: Chemical potential in eV
    @param    temp: Temperature in eV
    @param    Smat: Stoner interaction matrix [Nchi, Nchi] float64
    @param    hamk: Normal-state k-space Hamiltonian [Nk, Norb, Norb] complex128
    @param delta_k: Gap function (anomalous potential) on k-mesh [Nk, Norb, Norb] complex128
    @param   klist: k-point list [Nk, 3] float64
    @param   qlist: q-point path list [Nq, 3] float64
    @param   olist: Orbital index pairs for susceptibility [Nchi, 2] int64 (chiolist)
    @param      Nw: Number of frequency points
    @param    Emax: Maximum frequency in eV
    @param  idelta: Lorentzian broadening in eV
    @param sw_spsym: True for triplet (dz) symmetry, False for singlet
    @retval   chisq: Spin susceptibility trace at each q-point [Nq, Nw] complex128
    @retval chis_orbq: Orbital-resolved spin susceptibility [Nq, ...]
    @retval   wlist: Frequency mesh [Nw] float64
    """
    Norb = hamk.shape[1]
    hamBdGk = flibs.mkBdGhamk(hamk - mu * np.eye(Norb), delta_k)
    eig_BdG, uni_BdG = flibs.get_eig(hamBdGk)
    # mu=0: chemical potential already absorbed into hamBdGk, BdG energies measured from 0.
    ffermi_BdG = flibs.get_ffermi(eig_BdG, 0., temp)
    wlist = np.linspace(0, Emax, Nw)
    chisq = []
    chis_orbq = []
    with open('chi0_sc.dat', 'w') as f, open('writeq.dat', 'w') as fq:
        for i, q in enumerate(qlist):
            fq.write(f'{i:d} {q[0]:5.3f} {q[1]:5.3f} {q[2]:5.3f}\n')
            fq.flush()
            qshift = flibs.get_qshift(klist, q)
            chi0 = flibs.get_chi_irr_sc(uni_BdG, eig_BdG, ffermi_BdG, qshift, olist, wlist, idelta, temp, sw_spsym)
            chis = flibs.get_chis(chi0, Smat)
            trchis, trchi0, chis_orb = flibs.get_tr_chi(chis, chi0, olist)
            chisq.append(trchis)
            chis_orbq.append(chis_orb)
            for w, trchi in zip(wlist, trchi0):
                f.write(f'{i:8.4f} {w:8.4f} {trchi.imag:9.4f} {trchi.real:9.4f}\n')
            f.write('\n')
            f.flush()
    return np.array(chisq), np.array(chis_orbq), wlist

def chis_qmap(Nx: int, Ny: int, Ecut: float, mu: float, temp: float, Smat: np.ndarray,
              klist: np.ndarray, olist: np.ndarray, eig: np.ndarray, uni: np.ndarray,
              idelta: float = 1.e-3) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    @fn chis_qmap
    @brief Compute the spin susceptibility chi_s map on the qx-qy plane at a fixed energy cutoff.
    @param     Nx: Number of q-points along qx
    @param     Ny: Number of q-points along qy
    @param   Ecut: Energy cutoff for the susceptibility integration in eV
    @param     mu: Chemical potential in eV
    @param   temp: Temperature in eV
    @param   Smat: Stoner interaction matrix [Norb, Norb]
    @param  klist: k-point list [Nk, 3]
    @param  olist: Orbital index list for susceptibility calculation
    @param    eig: Eigenvalues array [Nk, Nband]
    @param    uni: Eigenvectors array [Nk, Norb, Nband]
    @param idelta: Broadening parameter (Lorentzian width) in eV
    @retval  chis: Spin susceptibility map [Ny, Nx]
    @retval  chi0: Bare susceptibility map [Ny, Nx]
    @retval    qx: qx coordinate mesh [Ny, Nx]
    @retval    qy: qy coordinate mesh [Ny, Nx]
    """
    ffermi=flibs.get_ffermi(eig,mu,temp)
    chis,chi0=flibs.chis_qmap(uni,eig,ffermi,klist,Smat,olist,Nx,Ny,temp,Ecut,idelta)
    x0=np.linspace(0,1,Nx,False)
    y0=np.linspace(0,1,Ny,False)
    qx,qy=np.meshgrid(x0,y0)
    return chis,chi0,qx,qy

def phi_spectrum(mu: float, temp: float, klist: np.ndarray, qlist: np.ndarray, olist: np.ndarray,
                 eig: np.ndarray,   uni: np.ndarray, Nw: int, Emax: float,
                 idelta: float = 1.e-3) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    @fn phi_spectrum
    @brief Compute the pairing susceptibility phi spectrum along a q-path (anomalous susceptibility).
    @param      mu: Chemical potential in eV
    @param    temp: Temperature in eV
    @param   klist: k-point list [Nk, 3]
    @param   qlist: q-point path list [Nq, 3]
    @param   olist: Orbital index list for phi calculation
    @param     eig: Eigenvalues array [Nk, Nband]
    @param     uni: Eigenvectors array [Nk, Norb, Nband]
    @param      Nw: Number of frequency points
    @param    Emax: Maximum frequency in eV
    @param  idelta: Broadening parameter (Lorentzian width) in eV
    @retval    phiq: Pairing susceptibility trace at each q-point [Nq, Nw]
    @retval phi_orbq: Orbital-resolved pairing susceptibility [Nq, ...]
    @retval   wlist: Frequency mesh array [Nw]
    """
    ffermi=flibs.get_ffermi(eig,mu,temp)
    wlist=np.linspace(0,Emax,Nw)
    phiq=[]
    phi_orbq=[]
    with open('writeq.dat','w') as fq:
        for i,q in enumerate(qlist):
            fq.write(f'{i:d} {q[0]:5.3f} {q[1]:5.3f} {q[2]:5.3f}\n')
            fq.flush()
            qshift=flibs.get_iqshift(klist,q)
            phi=flibs.get_phi_irr(uni,eig,ffermi,qshift,olist,wlist,idelta,mu,temp)
            trphi,phi_orb=flibs.get_tr_phi(phi,olist)
            phiq.append(trphi)
            phi_orbq.append(phi_orb)
    return np.array(phiq),np.array(phi_orbq),wlist

def phi_qmap(Nx: int, Ny: int, Ecut: float, mu: float, temp: float, klist: np.ndarray,
             olist: np.ndarray, eig: np.ndarray, uni: np.ndarray, idelta: float = 1.e-3,
             sw_omega: bool = True) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    @fn phi_qmap
    @brief Compute the pairing susceptibility phi map on the qx-qy plane at a fixed energy cutoff.
    @param      Nx: Number of q-points along qx
    @param      Ny: Number of q-points along qy
    @param    Ecut: Energy cutoff for the susceptibility integration in eV
    @param      mu: Chemical potential in eV
    @param    temp: Temperature in eV
    @param   klist: k-point list [Nk, 3]
    @param   olist: Orbital index list for phi calculation
    @param     eig: Eigenvalues array [Nk, Nband]
    @param     uni: Eigenvectors array [Nk, Norb, Nband]
    @param  idelta: Broadening parameter (Lorentzian width) in eV
    @param sw_omega: If True, integrate over Matsubara frequencies; if False, use static limit
    @retval    phi: Pairing susceptibility map [Ny, Nx]
    @retval     qx: qx coordinate mesh [Ny, Nx]
    @retval     qy: qy coordinate mesh [Ny, Nx]
    """
    ffermi=flibs.get_ffermi(eig,mu,temp)
    phi=flibs.phi_qmap(uni,eig,ffermi,klist,olist,Nx,Ny,mu,temp,Ecut,idelta,sw_omega)
    x0=np.linspace(0,1,Nx,False)
    y0=np.linspace(0,1,Ny,False)
    qx,qy=np.meshgrid(x0,y0)
    return phi,qx,qy

def get_chi_orb_list(Norb: int, site_prof: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    @fn get_chi_orb_list
    @brief Generate the orbital index pair list for susceptibility calculation from site profile.
    @param    Norb: Total number of orbitals in the Hamiltonian
    @param site_prof: Array of orbital counts per site (e.g. [3, 3] for two 3-orbital sites)
    @retval chiolist: Array of orbital index pairs (1-based) [Npairs, 2]
    @retval    site: Array of site indices corresponding to each orbital pair [Npairs]
    """
    if(len(site_prof)==1):
        tmp=np.arange(Norb)+1
        o1,o2=np.meshgrid(tmp,tmp)
        chiolist=np.array([o1.flatten(),o2.flatten()]).T
        site=np.ones(len(chiolist),dtype=np.int64)
    else:
        if(Norb==sum(site_prof)):
            chiolist=[]
            site=[]
            N0=1
            for i_site in site_prof:
                tmp=np.arange(i_site)+N0
                o1,o2=np.meshgrid(tmp,tmp)
                chiolist+=list(np.array([o1.flatten(),o2.flatten()]).T)
                # One site label per orbital pair (i_site^2 entries), not per row (i_site).
                site+=[N0]*o1.size
                N0+=i_site
            chiolist=np.array(chiolist)
            site=np.array(site,dtype=np.int64)
        else:
            print("site_prof doesn't correspond to Hamiltonian")
            exit()
    return chiolist,site

def gap_symms(klist: np.ndarray, Norb: int, gap_sym: int):
    if gap_sym==0: #s
        init_gap=np.ones((Norb,len(klist)),dtype=np.float64)
    else:
        A=2*np.pi
        init_gap=np.zeros((Norb,len(klist)),dtype=np.float64)
        for i in range(Norb):
            if gap_sym==1: #dx2-y2
                init_gap[i,:]=np.cos(A*klist[:,0])-np.cos(A*klist[:,1])
            elif gap_sym==2: #spm
                init_gap[i,:]=2*np.cos(A*klist[:,0])*np.cos(A*klist[:,1])
            elif gap_sym==3: #dxy
                init_gap[i,:]=2*np.sin(A*klist[:,0])*np.sin(A*klist[:,1])
            elif gap_sym==4: #dxz
                init_gap[i,:]=2*np.sin(A*klist[:,0])*np.sin(A*klist[:,2])
            elif gap_sym==5: #dyz
                init_gap[i,:]=2*np.sin(A*klist[:,1])*np.sin(A*klist[:,2])
            elif gap_sym==-1: #px
                init_gap[i,:]=2*np.sin(A*klist[:,0])
            elif gap_sym==-2: #py
                init_gap[i,:]=2*np.sin(A*klist[:,1])
    return init_gap

def get_initial_gap(klist: np.ndarray, Norb: int, gap_sym: int) -> np.ndarray:
    """
    @fn get_initial_gap
    @brief Generate an initial gap function with the specified pairing symmetry for gap equation iteration.
    @param   klist: k-point list in fractional coordinates [Nk, 3]
    @param    Norb: Number of orbitals in the Hamiltonian
    @param gap_sym: Symmetry index: 0=s, 1=dx2-y2, 2=s+-, 3=dxy, -1=px, -2=py
    @return init_gap: Initial gap function array [Norb, Nk]
    """
    if gap_sym>=0:
        gapsym=['s','dx2-y2','spm','dxy']
        print('gap symmetry is '+gapsym[gap_sym])
    else:
        gapsym=['s','px','py']
        print('gap symmetry is '+gapsym[-gap_sym])
    return gap_symms(klist,Norb,gap_sym)
