#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
@file plibs.py
@package plibs
@brief python script for model calculations
"""
import libs.flibs as flibs
import numpy as np, scipy.optimize as scopt, scipy.linalg as sclin

def import_hoppings(fname:str,ftype:int) -> tuple[np.ndarray,np.ndarray,int,int]:
    """
    @fn import_hoppings()
    @brief This function imports hopping parameters from files
    @param fname: Name of import files
    @param ftype: File format of import files
    @retval  rvec: the array of r-vector, size:(nr,3)
    @retval ham_r: the array of hopping integrals, size:(nr,no,no)
    @retval    no: the number of orbitals
    @retval    nr: the number of r-vectors
    """
    def import_hop(name: str):
        rvec=np.loadtxt(f'{name}/irvec.txt')
        nr=rvec[:,0].size
        ndegen=np.loadtxt(f'{name}/ndegen.txt')
        tmp=np.array([complex(float(tp[0]),float(tp[1])) for tp in
                      [f.strip(' ()\n').split(',') for f in open(f'{name}/ham_r.txt','r')]])
        no=int(np.sqrt(tmp.size/nr))
        ham_r=(tmp.reshape(nr,no,no).T/ndegen).T.round(6).copy()
        return(rvec,ham_r,no,nr)

    def import_out(name:str):
        data=np.loadtxt(name)
        con=(data[:,:3]==data[0,:3]).prod(axis=1).sum()
        no,nr =int(np.sqrt(con)),data[:,0].size//con
        rvec=np.array(data[:nr,:3])
        ham_r=(data[:,3]+1j*data[:,4]).reshape(no*no,nr).T.reshape(nr,no,no).round(6).copy()
        return(rvec,ham_r,no,nr)

    def import_hr(name:str):
        tmp=[f.split() for f in open(f'{name}_hr.dat','r')]
        no, nr=int(tmp[1][0]), int(tmp[2][0])
        c2,tmp1=3,[]
        while not len(tmp1)==nr:
            tmp1.extend(tmp[c2])
            c2=c2+1
        ndegen=np.array([float(t) for t in tmp1])
        tmp1=[[float(t) for t in tp] for tp in tmp[c2:]]
        tmp=np.array([complex(tp[5],tp[6]) for tp in tmp1])
        rvec=np.array([tmp1[no*no*i][:3] for i in range(nr)])
        ham_r=(tmp.reshape((nr,no,no)).T/ndegen).T.round(6).copy()
        return(rvec,ham_r,no,nr)

    def import_Hopping(name:str,sw_axis=False):
        tmp=[f.split() for f in open(f'{name}/Hopping.dat','r')]
        axis=np.array([[float(tp) for tp in tpp] for tpp in tmp[1:4]])
        no,nr=int(tmp[4][0]),int(tmp[4][1])
        tmp1=np.array([[float(t) for t in tp] for tp in tmp[7+no:]])
        rvec=np.array([tmp1[no*no*i][:3] for i in range(nr)])
        tmp=np.array([complex(tp[8],tp[9]) for tp in tmp1])
        ham_r=tmp.reshape(nr,no,no)
        if sw_axis:
            return(rvec,ham_r,no,nr,axis)
        else:
            return(rvec,ham_r,no,nr)
    if ftype==0:
        rvec,ham_r,no,nr=import_hop(fname)
    elif ftype==1:
        rvec,ham_r,no,nr=import_out(fname)
    elif ftype==2:
        rvec,ham_r,no,nr=import_hr(fname)
    else:
        rvec,ham_r,no,nr=import_Hopping(fname)
    return(rvec,ham_r,no,nr)

def import_MLO_hoppings(name:str) -> tuple[np.ndarray,np.ndarray,np.ndarray,int,int]:
    """
    @fn import_MLO_hoppings()
    @brief This function imports MLO hopping parameters from files
    @param name: File name
    @retval  rvec: the array of r-vector, size:(nr,3)
    @retval ham_r: the array of hopping integrals, size:(nr,no,no)
    @retval   S_r: the array of overlap integrals, size:(nr,no,no)
    @retval    no: the number of orbitals
    @retval    nr: the number of r-vectors
    """
    tmp=[f.split() for f in open(f'{name}','r')]
    tmp1=np.array([[float(t) for t in tp] for tp in tmp])
    no=int(tmp1[:,0].max())
    nr=int(len(tmp1)/(no*no))
    tmp=np.array([complex(tp[5],tp[6]) for tp in tmp1])
    tmpS=np.array([complex(tp[7],tp[8]) for tp in tmp1])
    rvec=np.array([tmp1[i][2:5] for i in range(nr)])
    ham_r=tmp.reshape((no*no,nr)).T.reshape((nr,no,no)).round(6).copy()*13.6
    S_r=tmpS.reshape((no*no,nr)).T.reshape((nr,no,no)).round(6).copy()
    return rvec,ham_r,S_r,no,nr

def get_bvec(avec: np.ndarray) -> np.ndarray:
    """
    @fn get_bvec()
    @brief This function generates reciprocal lattice vector from primitive translation vector
    @param  avec: primitive translation vector
    @return bvec: reciprocal lattice vector
    """
    bvec=2*np.pi*sclin.inv(avec).T
    return bvec

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
    no=int(eig.size/len(eig))
    def func(mu):
        sum_fermi=flibs.get_ffermi(eig,mu,temp).sum()
        return(fill*Nk-sum_fermi)
    emax=eig.max()
    emin=eig.min()
    try:
        mu=scopt.brentq(func,emin,emax)
    except ValueError:
        # filling is outside the computed band range; clamp to the closer boundary
        if func(emin)*func(emax) > 0:
            mu=emin if abs(func(emin)) < abs(func(emax)) else emax
            print(f"Warning: calc_mu could not bracket the chemical potential. Clamped to mu={mu:.4f} eV",flush=True)
        else:
            raise
    return mu

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
    itemp=1./temp
    def func(mu):
        return(fill*Nsite-0.5*(1.0-np.tanh(0.5*(eigs-mu)*itemp)).sum())
    emax=eigs.max()
    emin=eigs.min()
    mu=scopt.brentq(func,emin,emax)
    return mu

def gen_rlist(Nx: int, Ny: int, Nz: int) -> np.ndarray:
    """
    @fn gen_rlist
    @brief This function generates the list of r-vectors
    @param     Nx: Number of sites along x-axis
    @param     Ny: Number of sites along y-axis
    @param     Nz: Number of sites along z-axis
    @return rlist: The array of r-vectors
    """
    x0=np.linspace(0,Nx,Nx,False)
    y0=np.linspace(0,Ny,Ny,False)
    z0=np.linspace(0,Nz,Nz,False)
    x,y,z=np.meshgrid(x0,y0,z0)
    rlist=np.array([x.ravel(),y.ravel(),z.ravel()]).T.copy()
    
    return rlist

def gen_klist_with_kmap(Nx: int, Ny: int, Nz: int) -> tuple[np.ndarray, np.ndarray]:
    """
    @fn gen_klist_with_kmap
    @brief This function generates the list of k-vectors
    @param     Nx: Number of kx mesh
    @param     Ny: Number of ky mesh
    @param     Nz: Number of kz mesh
    @retval klist: The array of k-vectors
    @retval  kmap: The array of property of k-points
    """
    x0=np.linspace(0,Nx,Nx,False,dtype=int)
    y0=np.linspace(0,Ny,Ny,False,dtype=int)
    z0=np.linspace(0,Nz,Nz,False,dtype=int)
    x,y,z=np.meshgrid(x0,y0,z0)
    kmap=np.array([x.ravel(),y.ravel(),z.ravel()]).T.copy()
    klist=np.array([x.ravel()/Nx,y.ravel()/Ny,z.ravel()/Nz]).T.copy()
    return klist,kmap

def gen_klist(Nx: int, Ny: int, Nz: int | None = None, sw_pp: bool = True, kz: float = 0.0) -> tuple[int, np.ndarray]:
    """
    @fn gen_klist()
    @brief This function generates a k-point list meshed in Nx,Ny,Nz
    @param     Nx: Number of axis 1 mesh (usually kx mesh)
    @param     Ny: Number of axis 2 mesh (usually ky mesh)
    @param     Nz: Number of axis 3 mesh (usually kz mesh)
    @param  sw_pp: switch output 2D mesh or 3D mesh
    @param     kz: value of axis 3 at 2D mesh (use only if sw_pp is True)
    @retval    Nk: number of k-points
    @retval klist: list of k-points
    """
    if sw_pp:
        kx=np.linspace(-0.5,0.5,Nx,True)
        ky=np.linspace(-0.5,0.5,Ny,True)
        if Nz is None:
            kz=np.array([kz])
        else:
            kz=np.linspace(-0.5,0.5,Nz,True)
    else:
        kx=np.linspace(0,1,Nx,False)
        ky=np.linspace(0,1,Ny,False)
        kz=np.linspace(0,1,Nz,False)
    x,y,z=np.meshgrid(kx,ky,kz)
    klist=np.array([x.ravel(),y.ravel(),z.ravel()]).T.copy()
    return len(klist),klist

def mk_klist(k_list: np.ndarray | list, N: int, bvec: np.ndarray) -> tuple[np.ndarray, np.ndarray, list]:
    """
    @fn mk_klist()
    @brief This function generates k-point list for symmetry line
    @param  k_list: the list of symmetry points
    @param       N: Number of mesh between symmetry points
    @param    bvec: reciprocal lattice vector
    @retval  klist: list of k-points at band plot
    @retval  splen: length of symmetry points
    @retval xticks: footnotes of klist at ticks of symmetry points
    """
    klist=[]
    splen=[]
    xticks=[]
    maxsplen=0
    for ks,ke in zip(k_list,k_list[1:]):
        dkv=np.array(ke)-np.array(ks)
        dkv_length=abs(dkv.dot(bvec)).sum()
        tmp=np.linspace(ks,ke,N,False)
        tmp2=np.linspace(0,dkv_length,N,False)+maxsplen
        maxsplen=tmp2.max()
        xticks+=[tmp2[0]]
        klist+=tmp.tolist()
        splen+=tmp2.tolist()
    klist+=[k_list[-1]]
    splen+=[maxsplen+dkv_length/N]
    xticks+=[splen[-1]]
    return np.array(klist),np.array(splen),xticks

def mk_qlist(k_set: np.ndarray | list, Nx: int, Ny: int, Nz: int, bvec: np.ndarray) -> tuple[np.ndarray, np.ndarray, list]:
    """
    @fn mk_qlist()
    @brief This function generates q-point list for symmetry line
    @param   k_set: the list of symmetry points
    @param      Nx: Number of x-mesh
    @param      Ny: Number of y-mesh
    @param      Nz: Number of z-mesh
    @param    bvec: reciprocal lattice vector
    @retval  qlist: list of q-points at band plot
    @retval  splen: length of symmetry points
    @retval xticks: footnotes of qlist at ticks of symmetry points
    """
    qlist=[]
    splen=[]
    xticks=[]
    maxsplen=0
    Narray=np.array([Nx,Ny,Nz])
    for ks,ke in zip(k_set,k_set[1:]):
        dk=np.array(ke)-np.array(ks)
        dk_length=abs(dk.dot(bvec)).sum()
        # ensure at least one division along non-zero component; use ceil to avoid zero due to truncation
        dN=np.asarray(np.ceil(abs(dk)*Narray),dtype=int)
        nonzero=dN[dN>0]
        if len(nonzero)==0:
            continue  # identical consecutive k-points; skip this segment
        N=nonzero.min()
        tmp=np.linspace(ks,ke,N,False)
        tmp2=np.linspace(0,dk_length,N,False)+maxsplen
        maxsplen+=dk_length
        xticks+=[tmp2[0]]
        qlist+=tmp.tolist()
        splen+=tmp2.tolist()
    qlist+=[k_set[-1]]
    splen+=[maxsplen+dk_length/N]
    xticks+=[splen[-1]]
    return np.array(qlist),np.array(splen),xticks

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
        if(e.max()*e.min() < 0. ):
            cont=sk.find_contours(e.reshape(mesh+1,mesh+1),0)
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
    cand_kz=[kz_arr[0],kz_arr[-1]]
    for i in range(len(dS)-1):
        if dS[i]*dS[i+1]<0:
            sign=-1. if dS[i]>0 else 1.
            res=minimize_scalar(lambda kz,s=sign: s*S_at_kz(kz),
                                bounds=(kz_arr[i],kz_arr[i+2]),method='bounded',
                                options={'xatol':1e-4,'maxiter':10})
            cand_kz.append(res.x)
    return cand_kz

def mk_kf(mesh: int, rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray, RotMat: np.ndarray, mu: float, kz: float) -> tuple[list, list]:
    """
    @fn mk_kf()
    @brief This function obtains 2d fermi wave-number
    @param   mesh: k-mesh
    @param   rvec: r vector of hoppings
    @param  ham_r: hopping parameters
    @param    S_r: overlap integrals
    @param RotMat: rotation matrix
    @param     mu: chemical potential
    @param     kz: kz value
    """
    eig=get_eigs_2d(mesh,rvec,ham_r,S_r,RotMat,kz)
    return get_kf_points(eig,mesh,mu,kz)

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

def get_ptv(alatt: np.ndarray, deg: np.ndarray, brav: int) -> tuple[np.ndarray, np.ndarray]:
    """
    @fn get_ptv
    @brief Generate primitive translation vectors (avec) and rotation matrix (Arot) for a given Bravais lattice type.
    @param  alatt: Lattice constants [a, b, c] in Angstrom
    @param    deg: Lattice angles [alpha, beta, gamma] in degrees (used for brav>=8)
    @param   brav: Bravais lattice type index (0=simple cubic, 1=FCC, 2=BCC, 3=hexagonal,
                   4=trigonal, 5=base-centered, 6=FCC2, 7=BCC2, >=8=general triclinic)
    @retval  avec: Primitive translation vectors as rows (shape [3,3])
    @retval  Arot: Rotation matrix defining the primitive cell in terms of conventional cell
    """
    if brav==0: #simple
        Arot=np.array([[ 1., 0., 0.],[ 0., 1., 0.],[ 0.,0., 1.]])
    elif brav==1: #face center
        Arot=np.array([[-.5, 0., .5],[0., .5, .5],[-.5,.5, 0.]])
    elif brav==2: #body center
        Arot=np.array([[.5,-.5, .5],[.5, .5, .5],[-.5,-.5, .5]])
    elif brav==3: #hexagonal
        Arot=np.array([[1., 0., 0.],[-.5,.5*np.sqrt(3.),0.],[0., 0., 1.]])
    elif brav==4: #trigonal
        cg=np.cos(np.pi*deg[2]/180.)
        tx,ty,tz=np.sqrt((1-cg)*0.5),np.sqrt((1.-cg)/6.),np.sqrt((1.+2*cg)/3.)
        Arot=np.array([[ tx,  -ty, tz],[ 0., 2*ty, tz],[-tx,  -ty, tz]])
    elif brav==5: #base center
        Arot=np.array([[ .5, .5, 0.],[-.5, .5, 0.],[ 0., 0., 1.]])
    elif brav==6: #face center 2
        Arot=np.array([[.5, 0., .5],[0., .5, .5],[.5,.5, 0.]])
    elif brav==7: #body center 2
        Arot=np.array([[-.5, .5, .5],[.5, -.5, .5],[.5, .5, -.5]])
    else:
        phase=[np.pi*deg[0]/180.,np.pi*deg[1]/180.,np.pi*deg[2]/180.]
        ax1,a2=alatt[1]/alatt[0],alatt[2]/alatt[0]
        r1,r2,r3=ax1*np.cos(phase[2]),ax1*np.sin(phase[2]),ax2*np.cos(phase[1])
        r4=ax2*(np.cos(phase[0])-np.cos(phase[1])*np.cos(phase[2]))/np.sin(phase[2])
        r5=ax2*np.sqrt(1+2*np.cos(phase[0])*np.cos(phase[1])*np.cos(phase[2])
                       -(np.cos(phase[0])**2+np.cos(phase[1])**2+np.cos(phase[2])**2))/np.sin(phase[2])
        Arot=np.array([[ 1., 0.,  0.],
                      [ r1, r2,  0.],
                      [ r3, r4, r5]])
    if brav in {0,1,2,3,4,5,6,7}:
        avec=alatt*Arot
    else:
        avec=alatt[0]*Arot
    return avec,Arot

def get_symm_line(brav:int)->tuple[list,list]:
    """
    @fn get_symm_line
    @brief Return the high-symmetry k-point path and corresponding axis labels for a given Bravais lattice.
    @param    brav: Bravais lattice type index (same convention as get_ptv)
    @retval k_list: List of high-symmetry k-points in fractional coordinates
    @retval xlabel: List of label strings (LaTeX) for each high-symmetry point
    """
    if brav==0: #simple
        k_list=[[0.,0.,.5],[0., 0., 0.],[.5, 0., 0.],[.5, .5, 0.],[0.,0.,0.]]
        xlabel=['Z',r'$\Gamma$','X','M',r'$\Gamma$']
    elif brav==1: #face center
        k_list=[[0.,0.,0.],[.5, 0., .5],[1., 0., 0.],[.5, .5, .5],[.5,.25,.75],[0.,0.,0.]]
        xlabel=[r'$\Gamma$','X',r'$\Gamma$','L','W',r'$\Gamma$']
    elif brav==2: #body center
        k_list=[[.5,.5,.5],[0., 0., 0.],[.5, 0., 0.],[.5, .5,-.5],[0.,0.,0.]]
        xlabel=['Z',r'$\Gamma$','X','M',r'$\Gamma$']
    elif brav==3: #hexagonal
        k_list=[[0.,0.,0.],[2./3.,-1./3., 0.],[.5, 0., 0.],[0., 0., 0.],[0.,0.,.5]]
        xlabel=[r'$\Gamma$','K','M',r'$\Gamma$','Z']
    elif brav==4: #trigonal
        k_list=[[0.,0.,0.],[.5,0.,.5],[.5,0.,0.],[0.,0.,0.],[.5,.5,.5]]
        xlabel=[r'$\Gamma$','K','M',r'$\Gamma$','Z']
    elif brav==5: #base center
        k_list=[[0.,0.,0.],[.5,0.,0.],[.5,.5,0.],[0.,0.,0.],[0.,0.,.5]]
        xlabel=[r'$\Gamma$','X','S',r'$\Gamma$','Z']
    elif brav==6:
        k_list=[[0.,0.,0.],[0., .5, .5],[.5, .5, .5],[.25,.75,.5],[0.,0.,0.]]
        xlabel=[r'$\Gamma$','X','L','W',r'$\Gamma$']
    elif brav==7: #body center (common)
        k_list=[[.5,.5,.5],[0., 0., 0.],[.5, 0., 0.],[.5, .5,-.5],[0.,0.,0.]]
        xlabel=['Z',r'$\Gamma$','X','M',r'$\Gamma$']
    else:
        k_list=[[0.,0.,0.],[.5,0.,0.],[.5,.5,0.],[0.,0.,0.]]
        xlabel=[r'$\Gamma$','X','M',r'$\Gamma$']
    return k_list,xlabel

def BZedge(brav: int) -> None:
    """
    @fn BZedge
    @brief Placeholder function for Brillouin zone edge drawing (not yet implemented).
    @param brav: Bravais lattice type index
    @return None
    """
    pass

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
    L11=[]
    L12=[]
    L22=[]
    for w in wlist:
        L11w,L12w,L22w=flibs.calc_Lij(eig,vk,ffermi,mu,w,idelta,temp)
        L11.append(L11w)
        L12.append(L12w)
        L22.append(L22w)
    L11=np.array(L11)
    L12=np.array(L12)
    L22=np.array(L22)
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
                site+=[N0]*len(o1)
                N0+=i_site
            chiolist=np.array(chiolist)
            site=np.array(site,dtype=np.int64)
        else:
            print("site_prof doesn't correspond to Hamiltonian")
            exit()
    return chiolist,site

def get_initial_gap(kmap: np.ndarray, klist: np.ndarray, Norb: int, gap_sym: int) -> np.ndarray:
    """
    @fn get_initial_gap
    @brief Generate an initial gap function with the specified pairing symmetry for gap equation iteration.
    @param    kmap: k-point index map [Nk, 3] (integer grid indices)
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
            elif gap_sym==-1: #px
                init_gap[i,:]=2*np.sin(A*klist[:,0])
            elif gap_sym==-2: #py
                init_gap[i,:]=2*np.sin(A*klist[:,1])
    return init_gap

def calc_carrier(rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray, avec: np.ndarray,
                 Nx: int, Ny: int, Nz: int, fill: float, temp: float, with_spin: bool = False) -> np.ndarray:
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

def read_epa_output(filename: str) -> tuple[int,int,np.ndarray,np.ndarray,np.ndarray,np.ndarray,np.ndarray]:
    """
    @fn read_epa_output
    @brief Read epa.x (job='egrid') output file and return EPA data arrays.
    @param filename: Path to the epa.x output file
    @retval   ngrid: Number of energy grids (typically 2: valence, conduction)
    @retval  nmodes: Number of phonon modes
    @retval    edge: Grid edges [ngrid] (eV)
    @retval    step: Grid steps [ngrid] (eV)
    @retval    nbin: Number of bins per grid [ngrid] int64
    @retval    wavg: Averaged phonon frequencies [nmodes] (eV)
    @retval    gavg: EPA averaged |g|^2 [ngrid, nbin_max, nbin_max, nmodes] (eV^2)
    """
    cm2ev = 1.23981e-4
    with open(filename) as f:
        ngrid, nmodes = map(int, f.readline().split())
        edge = np.zeros(ngrid)
        step = np.zeros(ngrid)
        nbin = np.zeros(ngrid, dtype=np.int64)
        for ii in range(ngrid):
            tokens = f.readline().split()
            edge[ii] = float(tokens[0])
            step[ii] = float(tokens[1])
            nbin[ii] = int(tokens[2])
        wavg_cm = np.array(f.readline().split(), dtype=np.float64)
        wavg = wavg_cm * cm2ev  # cm^-1 -> eV
        nbin_max = int(np.max(nbin))
        gavg = np.zeros((ngrid, nbin_max, nbin_max, nmodes), dtype=np.float64)
        for line in f:
            tokens = line.split()
            if len(tokens) < 3 + nmodes:
                continue
            ii = int(tokens[0]) - 1   # 0-based
            jj = int(tokens[1]) - 1
            kk = int(tokens[2]) - 1
            g2 = np.array(tokens[3:3+nmodes], dtype=np.float64)
            gavg[ii, kk, jj, :] = g2  # eV^2
    return ngrid, nmodes, edge, step, nbin, wavg, gavg
