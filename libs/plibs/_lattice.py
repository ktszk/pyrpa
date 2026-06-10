#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Lattice and reciprocal space utilities: hopping import, k/r-mesh generation, BZ visualization.
"""
import numpy as np
import scipy.linalg as sclin

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
        # Divide each R-block by its degeneracy weight (number of equivalent k-points)
        ham_r=(tmp.reshape(nr,no,no).T/ndegen).T.round(6).copy()
        return(rvec,ham_r,no,nr)

    def import_out(name:str):
        data=np.loadtxt(name)
        # Count rows that share the same R-vector as the first row → gives no^2 (orbital pairs per R)
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
    ham_r=tmp.reshape((no*no,nr)).T.reshape((nr,no,no)).round(6).copy()*13.6  # Rydberg -> eV
    S_r=tmpS.reshape((no*no,nr)).T.reshape((nr,no,no)).round(6).copy()
    return rvec,ham_r,S_r,no,nr

def get_bvec(avec: np.ndarray) -> np.ndarray:
    """
    @fn get_bvec()
    @brief This function generates reciprocal lattice vector from primitive translation vector
    @param  avec: primitive translation vector
    @return bvec: reciprocal lattice vector
    """
    # b_i = 2pi * (a^{-T})_i  (standard reciprocal lattice relation a_i . b_j = 2pi delta_ij)
    bvec=2*np.pi*sclin.inv(avec).T
    return bvec

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
        # Symmetric BZ: k in [-0.5, 0.5] for band/FS plots (includes both zone boundaries)
        kx=np.linspace(-0.5,0.5,Nx,True)
        ky=np.linspace(-0.5,0.5,Ny,True)
        if Nz is None:
            kz=np.array([kz])
        else:
            kz=np.linspace(-0.5,0.5,Nz,True)
    else:
        # Periodic BZ: k in [0, 1) without endpoint for FFT-compatible k-summations
        kx=np.linspace(0,1,Nx,False)
        ky=np.linspace(0,1,Ny,False)
        kz=np.linspace(0,1,Nz,False)
    # indexing='ij' so that eig.reshape(Nx,Ny[,Nz]) has axes ordered (kx,ky[,kz]);
    # the default 'xy' would swap the kx/ky axes in FS contour/marching-cubes plots
    x,y,z=np.meshgrid(kx,ky,kz,indexing='ij')
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
        # Arc length of segment in Cartesian reciprocal space (Angstrom^-1)
        dkv_length=np.linalg.norm(dkv.dot(bvec))
        tmp=np.linspace(ks,ke,N,False)       # N points, excluding endpoint (appended at end)
        tmp2=np.linspace(0,dkv_length,N,False)+maxsplen
        # advance by the full segment length so consecutive segments join seamlessly
        # (tmp2.max() would be one step short and shift every following xtick)
        maxsplen+=dkv_length
        xticks+=[tmp2[0]]
        klist+=tmp.tolist()
        splen+=tmp2.tolist()
    klist+=[k_list[-1]]        # append the final symmetry point
    splen+=[maxsplen]
    xticks+=[maxsplen]
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
    dk_length=0.0
    N=1
    for ks,ke in zip(k_set,k_set[1:]):
        dk=np.array(ke)-np.array(ks)
        # Arc length of segment in Cartesian reciprocal space (Angstrom^-1)
        dk_length=np.linalg.norm(dk.dot(bvec))
        # ensure at least one division along non-zero component; use ceil to avoid zero due to truncation
        # Compute per-axis step counts; use ceil so non-zero dk always yields ≥1 point
        dN=np.asarray(np.ceil(abs(dk)*Narray),dtype=int)
        nonzero=dN[dN>0]
        if len(nonzero)==0:
            dk_length=0.0  # identical consecutive k-points; skip this segment
            continue
        N=nonzero.min()   # use the tightest constraint to avoid oversampling any axis
        tmp=np.linspace(ks,ke,N,False)
        tmp2=np.linspace(0,dk_length,N,False)+maxsplen
        maxsplen+=dk_length
        xticks+=[tmp2[0]]
        qlist+=tmp.tolist()
        splen+=tmp2.tolist()
    qlist+=[k_set[-1]]
    splen+=[maxsplen]          # maxsplen already includes the full last segment
    xticks+=[maxsplen]
    return np.array(qlist),np.array(splen),xticks

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
        # Primitive vectors of rhombohedral cell expressed in orthogonal frame (Ashcroft & Mermin convention)
        tx,ty,tz=np.sqrt((1-cg)*0.5),np.sqrt((1.-cg)/6.),np.sqrt((1.+2*cg)/3.)
        Arot=np.array([[ tx,  -ty, tz],[ 0., 2*ty, tz],[-tx,  -ty, tz]])
    elif brav==5: #base center
        Arot=np.array([[ .5, .5, 0.],[-.5, .5, 0.],[ 0., 0., 1.]])
    elif brav==6: #face center 2
        Arot=np.array([[.5, 0., .5],[0., .5, .5],[.5,.5, 0.]])
    elif brav==7: #body center 2
        Arot=np.array([[-.5, .5, .5],[.5, -.5, .5],[.5, .5, -.5]])
    else:
        # General triclinic lattice: build orthogonal representation from lattice angles (alpha, beta, gamma)
        phase=[np.pi*deg[0]/180.,np.pi*deg[1]/180.,np.pi*deg[2]/180.]
        ax1,ax2=alatt[1]/alatt[0],alatt[2]/alatt[0]   # b/a, c/a ratios
        r1,r2,r3=ax1*np.cos(phase[2]),ax1*np.sin(phase[2]),ax2*np.cos(phase[1])
        sin_gamma=np.sin(phase[2])
        if abs(sin_gamma) < 1.0e-12:
            raise ValueError(f"Invalid gamma angle for triclinic cell: gamma={deg[2]} deg")
        r4=ax2*(np.cos(phase[0])-np.cos(phase[1])*np.cos(phase[2]))/sin_gamma
        # r5 = c * sqrt(1 - cos^2(alpha) - cos^2(beta) - cos^2(gamma) + 2cos(alpha)cos(beta)cos(gamma)) / sin(gamma)
        r5=ax2*np.sqrt(1+2*np.cos(phase[0])*np.cos(phase[1])*np.cos(phase[2])
                       -(np.cos(phase[0])**2+np.cos(phase[1])**2+np.cos(phase[2])**2))/sin_gamma
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

def BZedge(bvec: np.ndarray, ax=None) -> None:
    """
    @fn BZedge
    @brief Draw first Brillouin zone boundary edges using Wigner-Seitz construction in reciprocal space.
    @param bvec: reciprocal lattice vectors as rows, shape (3,3), in Angstrom^-1
    @param   ax: matplotlib Axes3D; uses plt.gca() if None
    """
    import itertools
    import matplotlib.pyplot as plt
    from scipy.spatial import Voronoi

    if ax is None:
        ax = plt.gca()

    # Generate reciprocal lattice points in Cartesian coordinates
    idx = range(-2, 3)
    pts = np.array([i*bvec[0] + j*bvec[1] + k*bvec[2]
                    for i, j, k in itertools.product(idx, repeat=3)])

    vor = Voronoi(pts)
    origin_idx = np.argmin(np.linalg.norm(pts, axis=1))

    # avec transforms Cartesian k to fractional*2pi plot coordinates: k_plot = avec @ k_cart
    avec = 2*np.pi * np.linalg.inv(bvec).T

    def _order_polygon(verts):
        c = verts.mean(axis=0)
        u = verts[0] - c
        u /= np.linalg.norm(u)
        v = None
        for i in range(1, len(verts)):
            n = np.cross(u, verts[i] - c)
            if np.linalg.norm(n) > 1e-10:
                v = np.cross(n / np.linalg.norm(n), u)
                break
        if v is None:
            return verts
        angles = np.arctan2([(p-c).dot(v) for p in verts], [(p-c).dot(u) for p in verts])
        return verts[np.argsort(angles)]

    for ridge_pts, ridge_verts in zip(vor.ridge_points, vor.ridge_vertices):
        if origin_idx not in ridge_pts or -1 in ridge_verts:
            continue
        verts_cart = vor.vertices[ridge_verts]
        verts_plot = (avec @ verts_cart.T).T
        ordered = _order_polygon(verts_plot)
        poly = np.vstack([ordered, ordered[0]])
        ax.plot(poly[:, 0], poly[:, 1], poly[:, 2], color='black', lw=0.5)
