#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Band structure, eigenvalues, chemical potential, and Fermi surface utilities.
"""
import numpy as np
import scipy.optimize as scopt
import scipy.linalg as sclin
import libs.flibs as flibs
from ._lattice import gen_klist


def get_eigs(klist: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray, rvec: np.ndarray,
             sw_uni: bool = False, sw_std: bool = False) -> tuple[np.ndarray, np.ndarray] | np.ndarray:
    """
    @fn get_eigs
    @brief This function generates eigenvalues of Hamiltonian
    @param  klist: list of k-points
    @param  ham_r: hopping parameters
    @param    S_r: overlap integrals
    @param   rvec: r vector of hoppings
    @param sw_uni: switch of output only unitary matrix or not
    @param sw_std: switch standardization of unitary matrix at MLO hoppings
    @retval   eig: eigenvalues of Hamiltonian
    @retval   uni: unitary matrix of eigenfunctions
    """
    if len(S_r)==0:
        ham_k=flibs.gen_ham(klist,ham_r,rvec)
        eig,uni=flibs.get_eig(ham_k)
        if sw_uni:
            return uni
        else:
            return eig,uni
    else:
        ham_k,S_k=flibs.gen_ham(klist,ham_r,rvec,Ovl_r=S_r)
        eig,uni=flibs.get_eig(ham_k,S_k)
        if sw_std:
            # Normalize each eigenvector column so orbital weights sum to 1 (MLO basis is non-orthogonal)
            norm=np.sqrt((abs(uni)**2).sum(axis=2))
            if sw_uni:
                return np.array([[u/nm for u,nm in zip(un,nor)] for un,nor in zip(uni,norm)])
            else:
                return eig,np.array([[u/nm for u,nm in zip(un,nor)] for un,nor in zip(uni,norm)])
        else:
            if sw_uni:
                return uni
            else:
                return eig,uni

def _solve_mu(func, emin: float, emax: float, label: str) -> float:
    """Bracket the chemical potential with brentq; if filling lies outside the band
    range, clamp to the closer boundary (shared by calc_mu and calc_mu_imp)."""
    try:
        return scopt.brentq(func, emin, emax)
    except ValueError:
        # filling is outside the computed band range; clamp to the closer boundary
        if func(emin)*func(emax) > 0:
            mu = emin if abs(func(emin)) < abs(func(emax)) else emax
            print(f"Warning: {label} could not bracket the chemical potential. Clamped to mu={mu:.4f} eV", flush=True)
            return mu
        raise

def calc_mu(eig,Nk,fill:float,temp:float)-> float:
    """
    @fn calc_mu()
    @brief This function obtains chemical potential mu
    @param   eig: Eigenvalues array
    @param    Nk: Number of k-points
    @param  fill: band filling
    @param  temp: Temperature
    @return   mu: chemical potential
    """
    if temp <= 0:
        raise ValueError(f"calc_mu requires temp > 0, got {temp}")
    def func(mu):
        sum_fermi=flibs.get_ffermi(eig,mu,temp).sum()
        return(fill*Nk-sum_fermi)  # zero when total electrons equals target filling * Nk
    return _solve_mu(func, eig.min(), eig.max(), 'calc_mu')

def calc_mu_imp(eigs,Nsite,fill:float,temp:float)-> float:
    """
    @fn calc_mu_imp
    @brief Calculate the chemical potential for an impurity system using the bisection method (brentq).
    @param   eigs: Eigenvalues array (all k-points and bands)
    @param  Nsite: Number of sites in the unit cell
    @param   fill: Target electron filling (electrons per site)
    @param   temp: Temperature in eV
    @return    mu: Chemical potential in eV
    """
    if temp <= 0:
        raise ValueError(f"calc_mu_imp requires temp > 0, got {temp}")
    itemp=1./temp
    def func(mu):
        return(fill*Nsite-0.5*(1.0-np.tanh(0.5*(eigs-mu)*itemp)).sum())
    return _solve_mu(func, eigs.min(), eigs.max(), 'calc_mu_imp')

def get_kf_points(eig: np.ndarray, mesh: int, mu: float, kz: float) -> tuple[list, list]:
    """
    @fn get_kf_points()
    @brief This function obtains fermi wave-number points from precomputed eigenvalues
    @param   eig: eigenvalues of Hamiltonian shape (Nk, Norb)
    @param  mesh: k-mesh
    @param    mu: chemical potential
    @param    kz: kz value (used to tag contour coordinates)
    """
    import skimage.measure as sk
    kf_points=[]
    fsband=[]
    for i,e in enumerate(eig.T-mu):
        if(e.max()*e.min() < 0. ):  # band crosses the Fermi level → has a Fermi surface sheet
            cont=sk.find_contours(e.reshape(mesh+1,mesh+1),0)
            # Rescale contour indices from [0, mesh] pixel space to [-0.5, 0.5] fractional BZ coords
            ct=[np.array([list(c)+[kz] for c in (cc-mesh/2)/mesh]) for cc in cont]
            kf_points.append(ct)
            fsband.append(i)
    return kf_points,fsband

def get_eigs_2d(mesh: int, rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray,
                RotMat: np.ndarray, kz: float) -> np.ndarray:
    """
    @fn get_eigs_2d()
    @brief Compute eigenvalues on 2D k-mesh at fixed kz (no contour finding)
    @param   mesh: k-mesh
    @param   rvec: r vector of hoppings
    @param  ham_r: hopping parameters
    @param    S_r: overlap integrals
    @param RotMat: rotation matrix
    @param     kz: kz value
    @return   eig: eigenvalues shape (Nk, Norb)
    """
    Nk,klist=gen_klist(mesh+1,mesh+1,kz=kz)
    # Rotate real-space vectors so that kz-slices are perpendicular to the magnetic field direction
    rvec1=RotMat.dot(rvec.T).T.copy()
    eig,_=get_eigs(klist,ham_r,S_r,rvec1)
    return eig

def shoelace_area(ct: np.ndarray) -> float:
    """
    @fn shoelace_area()
    @brief Compute area of a closed contour in normalized BZ coordinates via shoelace formula
    @param    ct: contour point array shape (N, 3), columns [kx, ky, kz] in normalized coords
    @return     : contour area in normalized BZ units
    """
    x,y=ct[:,0],ct[:,1]
    return 0.5*np.abs(np.dot(x,np.roll(y,-1))-np.dot(y,np.roll(x,-1)))

def get_band_area(v2: list, blist: list, band_idx: int, ABZ: float) -> float | None:
    """
    @fn get_band_area()
    @brief Compute total Fermi surface cross-section area for one band in AA^-2
           Sums over disconnected contour pieces. Returns None if band has no FS at this kz.
    @param      v2: list of contour lists, one entry per FS band
    @param   blist: list of band indices corresponding to entries in v2
    @param band_idx: target band index
    @param     ABZ: Brillouin zone area in AA^-2
    @return       : total area in AA^-2, or None if band_idx not in blist
    """
    if band_idx not in blist:
        return None
    j=blist.index(band_idx)
    return sum(shoelace_area(ct) for ct in v2[j])*ABZ

def make_rotmat(theta_deg: float, phi_deg: float) -> np.ndarray:
    """
    @fn make_rotmat()
    @brief Build rotation matrix R s.t. R @ B_hat = z_hat,
           so that kz-slices in the rotated frame are perpendicular to B.
           R = Ry(-theta) @ Rz(-phi)
    @param theta_deg: polar angle of B from z-axis [deg]
    @param  phi_deg: azimuthal angle of B [deg]
    @return  rotmat: 3x3 rotation matrix
    """
    th,ph=np.deg2rad(theta_deg),np.deg2rad(phi_deg)
    cp,sp=np.cos(ph),np.sin(ph)
    ct,st=np.cos(th),np.sin(th)
    # First rotate about z by -phi (align azimuth), then about y by -theta (align polar)
    # Result: R @ B_hat = z_hat, so dHvA cross-sections are perpendicular to B
    Rz=np.array([[ cp, sp, 0.],[-sp, cp, 0.],[0., 0., 1.]])   # Rz(-phi)
    Ry=np.array([[ ct, 0.,-st],[ 0., 1.,  0.],[ st, 0., ct]]) # Ry(-theta)
    return Ry@Rz

def scan_fs_area(mesh: int, rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray,
                 rotmat: np.ndarray, mu: float, ABZ: float,
                 meshkz: int=20) -> tuple[dict, dict]:
    """
    @fn scan_fs_area()
    @brief Phase-1 scan: compute cross-section area S(kz) at mu for all FS bands
    @param    mesh: k-mesh for 2D kx-ky grid
    @param    rvec: r vector of hoppings
    @param   ham_r: hopping parameters
    @param     S_r: overlap integrals
    @param  rotmat: rotation matrix mapping B_hat -> z_hat
    @param      mu: chemical potential
    @param     ABZ: Brillouin zone area in AA^-2
    @param  meshkz: number of kz scan points in [0, 0.5]
    @return S_scan: dict band_idx -> [(kz, S[AA^-2]), ...]
    @return eig_cache: dict kz -> eig array, reusable in Phase 2
    """
    kz0=np.linspace(0.,.5,meshkz,True)
    S_scan={}
    eig_cache={}
    for kz in kz0:
        eig=get_eigs_2d(mesh,rvec,ham_r,S_r,rotmat,kz)
        eig_cache[kz]=eig
        v2,blist=get_kf_points(eig,mesh,mu,kz)
        for band_idx in blist:
            S=get_band_area(v2,blist,band_idx,ABZ)
            S_scan.setdefault(band_idx,[]).append((kz,S))
    return S_scan,eig_cache

def find_extremal_kz(kz_arr: np.ndarray, S_arr: np.ndarray, band_idx: int,
                     mesh: int, rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray,
                     rotmat: np.ndarray, mu: float, ABZ: float) -> list:
    """
    @fn find_extremal_kz()
    @brief Phase-2 extremal search: find kz where dS/dkz=0 by sign-change detection,
           then refine each bracket with minimize_scalar. BZ boundaries always included.
    @param  kz_arr: kz values from Phase-1 scan
    @param   S_arr: S(kz) values from Phase-1 scan
    @param band_idx: target band index
    @param    mesh: k-mesh for 2D kx-ky grid
    @param    rvec: r vector of hoppings
    @param   ham_r: hopping parameters
    @param     S_r: overlap integrals
    @param  rotmat: rotation matrix mapping B_hat -> z_hat
    @param      mu: chemical potential
    @param     ABZ: Brillouin zone area in AA^-2
    @return cand_kz: list of extremal kz values
    """
    from scipy.optimize import minimize_scalar
    def S_at_kz(kz):
        eig=get_eigs_2d(mesh,rvec,ham_r,S_r,rotmat,kz)
        v2,blist=get_kf_points(eig,mesh,mu,kz)
        s=get_band_area(v2,blist,band_idx,ABZ)
        return s if s is not None else 0.
    dS=np.diff(S_arr)
    # BZ boundaries (kz=0 and kz=pi/2) are always extremal candidates
    cand_kz=[kz_arr[0],kz_arr[-1]]
    for i in range(len(dS)-1):
        if dS[i]*dS[i+1]<0:  # sign change in dS/dkz → local extremum in this bracket
            sign=-1. if dS[i]>0 else 1.
            res=minimize_scalar(lambda kz,s=sign: s*S_at_kz(kz),
                                bounds=(kz_arr[i],kz_arr[i+2]),method='bounded',
                                options={'xatol':1e-4,'maxiter':10})
            cand_kz.append(res.x)
    return cand_kz

def gen_3d_surf_points(mesh: int, rvec: np.ndarray, ham_r: np.ndarray,
                       S_r: np.ndarray, mu: float, kscale: float = 1.0) -> tuple[list, list, list]:
    """
    @fn gen_3d_surf_points()
    @brief This function obtains 3d fermi wave-numbers
    @param   mesh: k-mesh
    @param   rvec: r vector of hoppings
    @param  ham_r: hopping parameters
    @param    S_r: overlap integrals
    @param     mu: chemical potential
    @param kscale: change considering k-space area 1.0 is only 1st BZ
    """
    import skimage.measure as ski
    Nk,klist=gen_klist(mesh+1,mesh+1,mesh+1)
    klist=klist*kscale
    eig,uni=get_eigs(klist,ham_r,S_r,rvec)
    if isinstance(kscale,int):
        ks=kscale*np.array([1.,1.,1.])
    else:
        ks=np.array(kscale)
    fspolys=[]
    fscenters=[]
    fsband=[]
    for i,e in enumerate(eig.T-mu):
        if(e.max()*e.min()<0.):
            verts,faces, _, _=ski.marching_cubes(e.reshape(mesh+1,mesh+1,mesh+1),0,
                                                 spacing=(ks[0]*2*np.pi/mesh,ks[1]*2*np.pi/mesh,
                                                          ks[2]*2*np.pi/mesh))
            verts=verts-ks*np.pi
            fspolys.append(verts[faces])
            fscenters.append(verts[faces].mean(axis=1)*.5/np.pi)
            fsband.append(i)
    return fspolys,fscenters,fsband

def get_colors(klist: np.ndarray, blist, mrot, rvec: np.ndarray, ham_r: np.ndarray,
               S_r: np.ndarray, ol, color_option: int, sw_2d: bool = False):
    """
    @fn get_colors
    @brief Compute color weights for Fermi surface visualization (orbital weight or velocity magnitude).
    @param        klist: List of k-points on the Fermi surface contours
    @param        blist: List of band indices corresponding to each contour in klist
    @param         mrot: Rotation matrix for velocity projection
    @param         rvec: Real-space lattice vectors (Wannier R-vectors)
    @param        ham_r: Hamiltonian in real space (Wannier representation)
    @param          S_r: Overlap matrix in real space
    @param           ol: Orbital index or list of orbital indices for weight coloring
    @param color_option: 0 = no color, 1 = orbital weight, 2 = velocity magnitude
    @param        sw_2d: If True, treat klist as a nested list (2D Fermi surface mode)
    @return       clist: List of color arrays for each contour segment, or [] if color_option==0
    """
    def get_col(cl,ol):
        col=(np.abs(cl[:,ol])**2 if isinstance(ol,int)
             else (np.abs(cl[:,ol])**2).sum(axis=1)).round(4)
        return col
    if color_option==0:
        return []
    elif color_option==1: #orbital weight color
        if sw_2d:
            uni=[[get_eigs(k,ham_r,S_r,rvec,True,True)[:,b] for k in kk] for kk,b in zip(klist,blist)]
            clist=[[np.array([get_col(cl,ol[0]),get_col(cl,ol[1]),get_col(cl,ol[2])]).T for cl in clst]
                   for clst in uni]
        else:
            uni=[get_eigs(k,ham_r,S_r,rvec,True,True)[:,b] for k,b in zip(klist,blist)]
            clist=[np.array([get_col(clst,ol[0]),get_col(clst,ol[1]),get_col(clst,ol[2])]).T for clst in uni]
        return clist
    elif color_option==2: #velocity size color
        if sw_2d:
            uni=[[get_eigs(k,ham_r,S_r,rvec,True,False) for k in kk] for kk in klist]
            vk=[[flibs.get_veloc(k,ham_r,rvec,mrot,unkk)[:,b,:] for k,unkk in zip(kk,unk)]
                for kk,unk,b in zip(klist,uni,blist)]
            clist=[[np.sqrt((abs(vkkk)*abs(vkkk)).sum(axis=1)) for vkkk in vkk] for vkk in vk]
        else:
            uni=[get_eigs(k,ham_r,S_r,rvec,True,False) for k in klist]
            vk=[flibs.get_veloc(kk,ham_r,rvec,mrot,unk)[:,b,:] for kk,unk,b in zip(klist,uni,blist)]
            clist=[np.sqrt((abs(vkk)*abs(vkk)).sum(axis=1)) for vkk in vk]
        return clist

def get_emesh(Nx: int, Ny: int, Nz: int, ham_r: np.ndarray, S_r: np.ndarray, rvec: np.ndarray,
              avec: np.ndarray, sw_uni: bool = False, sw_veloc: bool = False, sw_mass: bool = False):
    """
    @fn get_emesh
    @brief Generate the band energy mesh on a full 3D k-grid with optional eigenvectors, velocities, and inverse mass.
    @param     Nx: Number of k-points along kx
    @param     Ny: Number of k-points along ky
    @param     Nz: Number of k-points along kz
    @param  ham_r: Hamiltonian in real space (Wannier representation)
    @param    S_r: Overlap matrix in real space
    @param   rvec: Real-space lattice vectors (Wannier R-vectors)
    @param   avec: Lattice vectors (used for velocity/mass calculation)
    @param sw_uni: If True, also return eigenvectors and k-point list
    @param sw_veloc: If True, compute and return group velocities
    @param sw_mass: If True (requires sw_veloc=True), also compute and return inverse effective mass tensor
    @return (Nk, eig, [vk], [imass], kweight): Tuple of k-point count, eigenvalues, and optionally velocities/mass
    """
    Nk, klist = gen_klist(Nx, Ny, Nz, sw_pp=False)
    eig, uni = get_eigs(klist, ham_r, S_r, rvec)
    kweight = np.ones(len(eig), dtype=np.float64)
    if sw_veloc:
        if sw_uni:
            vk=flibs.get_vnmk(klist,ham_r,rvec,avec,uni)
            return Nk,eig,vk,kweight
        else:
            vk=flibs.get_veloc(klist,ham_r,rvec,avec,uni)
            if sw_mass:
                imass=flibs.get_mass(klist,ham_r,rvec,avec,uni,True)
                return Nk,eig,vk,imass,kweight
            else:
                return Nk,eig,vk,kweight
    else:
        if sw_uni:
            return Nk,klist,eig,uni,kweight
        else:
            return Nk,eig,kweight

def calc_carrier(rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray, avec: np.ndarray,
                 Nx: int, Ny: int, Nz: int, fill: float, temp: float,
                 with_spin: bool = False) -> np.ndarray:
    """
    @fn calc_carrier
    @brief Calculate the carrier density (electrons/cm³) from the Fermi-Dirac distribution derivative.
    @param      rvec: Real-space lattice vectors (Wannier R-vectors)
    @param     ham_r: Hamiltonian in real space (Wannier representation)
    @param       S_r: Overlap matrix in real space
    @param      avec: Lattice vectors (rows are primitive vectors) in Angstrom
    @param        Nx: Number of k-points along kx
    @param        Ny: Number of k-points along ky
    @param        Nz: Number of k-points along kz
    @param      fill: Target electron filling (electrons per unit cell)
    @param      temp: Temperature in eV
    @param with_spin: If True, include spin degeneracy factor (not yet used)
    @return   n_carr: Carrier density array [Nband] in cm⁻³
    """
    Nk,eig,kweight=get_emesh(Nx,Ny,Nz,ham_r,S_r,rvec,avec.T)
    Vuc=sclin.det(avec)*1e-24
    mu=calc_mu(eig,Nk,fill,temp)
    dfermi=0.25*(1.-np.tanh(0.5*(eig-mu)/temp)**2)/temp
    n_carr=2*dfermi.sum(axis=0)/(Vuc*Nk)
    return n_carr
