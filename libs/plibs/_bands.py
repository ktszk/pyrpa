#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Band structure, eigenvalues, chemical potential, and Fermi surface utilities.
"""
import numpy as np
import scipy.optimize as scopt
import scipy.linalg as sclin
import libs.flibs as flibs
from ._lattice import gen_klist, get_bvec


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

# Orbits discarded because their contour ran into the edge of the 2D slice box.
# The rotated fractional box is not a Brillouin zone for a tilted field, so such a
# contour is a truncated piece of the real orbit and its shoelace area is meaningless.
_open_contour_state={'n':0}

def open_contour_report(reset: bool = True) -> int:
    """
    @fn open_contour_report()
    @brief Number of orbits get_band_area() discarded because the contour was open.
    @param reset: clear the counter after reading (default True)
    @return     : count since the last reset
    """
    n=_open_contour_state['n']
    if reset:
        _open_contour_state['n']=0
    return n

def get_band_area(v2: list, blist: list, band_idx: int, ABZ: float,
                  sw_strict: bool = True) -> float | None:
    """
    @fn get_band_area()
    @brief Compute total Fermi surface cross-section area for one band in AA^-2
           Sums over disconnected contour pieces. Returns None if band has no FS at this kz,
           or (with sw_strict) if any contour piece is open.
    @param      v2: list of contour lists, one entry per FS band
    @param   blist: list of band indices corresponding to entries in v2
    @param band_idx: target band index
    @param     ABZ: area Jacobian of the slice in AA^-2 (see slice_area_factor)
    @param sw_strict: reject orbits whose contour is not closed inside the slice box
    @return       : total area in AA^-2, or None if unavailable/unreliable
    """
    if band_idx not in blist:
        return None
    j=blist.index(band_idx)
    if sw_strict:
        for ct in v2[j]:
            # skimage repeats the first vertex on a closed contour; an open one ends
            # on the box edge, so the enclosed area is not the orbit area.
            if not np.allclose(ct[0],ct[-1]):
                _open_contour_state['n']+=1
                return None
    return sum(shoelace_area(ct) for ct in v2[j])*ABZ

def slice_area_factor(rotmat: np.ndarray, bvec: np.ndarray) -> float:
    """
    @fn slice_area_factor()
    @brief Area Jacobian [AA^-2] converting a contour area measured in the rotated
           fractional (k1,k2) box into a physical cross-section area.
           get_eigs_2d replaces R -> rotmat.R in the fractional basis, i.e.
           E'(k) = E(rotmat^T k), so a point of the kz-slice sits at Cartesian
           q = bvec^T rotmat^T k. The two in-plane tangents are therefore
           v_i = sum_k rotmat[i,k] b_k (i=1,2) and one unit of fractional area maps
           to |v1 x v2|. Reduces to 4*pi^2/(a1*a2) for an orthorhombic cell at rotmat=I;
           no orthogonality of rotmat is assumed.
    @param rotmat: rotation matrix applied to rvec (fractional basis)
    @param   bvec: reciprocal lattice vectors as rows [AA^-1]
    @return      : |v1 x v2| in AA^-2
    """
    v=np.asarray(rotmat,dtype=float)@np.asarray(bvec,dtype=float)
    fac=float(np.linalg.norm(np.cross(v[0],v[1])))
    if fac<=0.:
        raise ValueError('slice_area_factor: degenerate slice (rotmat rows 1,2 are collinear)')
    return fac

def field_direction(rotmat: np.ndarray, avec: np.ndarray) -> np.ndarray:
    """
    @fn field_direction()
    @brief Cartesian unit vector along which the kz-slices are actually stacked,
           i.e. the magnetic field direction the rotated Hamiltonian really describes.
           The slice normal is v1 x v2, which equals sum_i rotmat[2,i] a_i up to a factor.
    @param rotmat: rotation matrix applied to rvec (fractional basis)
    @param   avec: primitive translation vectors as rows [AA]
    @return      : unit vector B_hat in Cartesian coordinates
    """
    n=np.asarray(rotmat,dtype=float)[2].dot(np.asarray(avec,dtype=float))
    return n/np.linalg.norm(n)

def make_rotmat(theta_deg: float, phi_deg: float) -> np.ndarray:
    """
    @fn make_rotmat()
    @brief Build rotation matrix R = Ry(-theta) @ Rz(-phi), s.t. R @ B_hat = z_hat.
           get_eigs_2d applies it to the FRACTIONAL R-vectors, so the resulting
           kz-slices are perpendicular to B only for a cubic cell; for any other
           cell the realised field direction is sum_i R[2,i] a_i (see field_direction).
           Use make_rotmat_bfield() to hit a requested (theta,phi) exactly.
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

def make_rotmat_bfield(theta_deg: float, phi_deg: float, bvec: np.ndarray) -> np.ndarray:
    """
    @fn make_rotmat_bfield()
    @brief Build the fractional-basis matrix M whose kz-slices are EXACTLY perpendicular
           to B in Cartesian k-space, for an arbitrary (also non-orthogonal) lattice.
           The slice tangents are v_i = sum_k M[i,k] b_k, so v_i . B_hat = 0 requires
           rows 0,1 of M to be orthogonal to w_k = b_k . B_hat; taking M orthogonal with
           row 2 along w satisfies this and keeps the slice spacing uniform.
           For a cubic cell w is parallel to B_hat and this reduces to make_rotmat().
           Unlike make_rotmat(), which rotates the fractional R-vectors directly and
           therefore only realises the requested angle for a cubic lattice.
    @param theta_deg: polar angle of B from the Cartesian z-axis [deg]
    @param   phi_deg: azimuthal angle of B [deg]
    @param      bvec: reciprocal lattice vectors as rows [AA^-1]
    @return   rotmat: 3x3 matrix to hand to get_eigs_2d/scan_fs_area
    """
    th,ph=np.deg2rad(theta_deg),np.deg2rad(phi_deg)
    bhat=np.array([np.sin(th)*np.cos(ph),np.sin(th)*np.sin(ph),np.cos(th)])
    w=np.asarray(bvec,dtype=float).dot(bhat)   # w_k = b_k . B_hat
    n=w/np.linalg.norm(w)
    # complete {u1,u2,n} into a right-handed orthonormal triad; the in-plane choice is
    # arbitrary because |v1 x v2| is invariant under a rotation of u1,u2 within the plane
    t=np.array([1.,0.,0.]) if abs(n[0])<0.9 else np.array([0.,1.,0.])
    u1=np.cross(t,n); u1/=np.linalg.norm(u1)
    u2=np.cross(n,u1)
    return np.array([u1,u2,n])

# --------------------------------------------------------------------------- #
#  Cartesian-plane slicing (arbitrary real field direction)
#
#  The fractional-box slicing above cuts the [-1/2,1/2]^2 cell of the rotated
#  frame, which is a Brillouin zone only when the tilt is commensurate with the
#  lattice; anything larger than the cell leaves the box and is lost.
#  Here the plane perpendicular to B is sampled directly in Cartesian k, over a
#  window wide enough to enclose the orbit. H(k) is exactly periodic in the
#  fractional k it is fed, so the window may span any number of Brillouin zones
#  and no wrapping is needed -- which is what removes the commensurability
#  restriction on B_hat.
# --------------------------------------------------------------------------- #

def build_fs_energy_grid(rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray, mu: float,
                         mesh: int = 120, ewin: float = 0.0,
                         method: str = 'cubic') -> dict | None:
    """
    @fn build_fs_energy_grid()
    @brief Diagonalise once on a periodic fractional grid and wrap the result in an
           interpolator, so that an angle sweep costs lookups instead of eigen-solves.
           Only bands reaching the Fermi level (within ewin) are stored; sorted
           eigenvalues are continuous in k, so interpolating them is safe even across
           band crossings.
    @param   mesh: grid points per axis over the primitive cell
    @param   ewin: keep bands whose range comes within ewin of mu (0 = crossing bands only)
    @param method: interpolation order for scipy RegularGridInterpolator
    @return       : dict with 'interp' (callable on fractional k) and 'bands'
                    (original band indices), or None if no band reaches mu
    """
    from scipy.interpolate import RegularGridInterpolator
    Nk,klist=gen_klist(mesh,mesh,mesh,sw_pp=False)     # [0,1) grid, periodic
    eig,_=get_eigs(klist,ham_r,S_r,rvec)
    eig=eig.reshape(mesh,mesh,mesh,-1)
    lo,hi=eig.min(axis=(0,1,2)),eig.max(axis=(0,1,2))
    bands=[i for i in range(eig.shape[-1]) if lo[i]-ewin<=mu<=hi[i]+ewin]
    if not bands:
        return None
    npad=2 if method=='cubic' else 1
    grid=np.pad(eig[...,bands],((npad,npad),(npad,npad),(npad,npad),(0,0)),mode='wrap')
    ax=(np.arange(-npad,mesh+npad))/mesh
    interp=RegularGridInterpolator((ax,ax,ax),grid,method=method,
                                   bounds_error=False,fill_value=None)
    return {'interp':interp,'bands':bands,'mesh':mesh,'method':method}

def bfield_hat(theta_deg: float, phi_deg: float) -> np.ndarray:
    """
    @fn bfield_hat()
    @brief Cartesian unit vector of the field for polar/azimuthal angles.
    @param theta_deg: polar angle from the z-axis [deg]
    @param   phi_deg: azimuthal angle [deg]
    @return         : B_hat (3,)
    """
    th,ph=np.deg2rad(theta_deg),np.deg2rad(phi_deg)
    return np.array([np.sin(th)*np.cos(ph),np.sin(th)*np.sin(ph),np.cos(th)])

def plane_frame(bhat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    @fn plane_frame()
    @brief Right-handed orthonormal in-plane axes (e1,e2) with e1 x e2 = B_hat.
           The in-plane orientation is arbitrary; orbit areas do not depend on it.
    @param bhat: field direction (need not be normalised)
    @retval  e1: in-plane unit vector
    @retval  e2: in-plane unit vector
    """
    n=np.asarray(bhat,dtype=float)
    n=n/np.linalg.norm(n)
    t=np.array([1.,0.,0.]) if abs(n[0])<0.9 else np.array([0.,1.,0.])
    e1=np.cross(t,n); e1/=np.linalg.norm(e1)
    e2=np.cross(n,e1)
    return e1,e2

def get_eigs_plane(rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray, avec: np.ndarray,
                   bhat: np.ndarray, d: float, L: float, mesh: int,
                   egrid: dict | None = None) -> np.ndarray:
    """
    @fn get_eigs_plane()
    @brief Eigenvalues on a square Cartesian window of the plane k.B_hat = d:
           k = u*e1 + v*e2 + d*B_hat with u,v in [-L/2, L/2].
    @param   rvec: r vector of hoppings (fractional, NOT rotated -- the tilt lives
                   entirely in the plane parametrisation)
    @param  ham_r: hopping parameters
    @param    S_r: overlap integrals
    @param   avec: primitive translation vectors as rows [AA]
    @param   bhat: field direction (unit)
    @param      d: signed distance of the plane from the origin [AA^-1]
    @param      L: window edge length [AA^-1]
    @param   mesh: number of intervals per edge (grid is (mesh+1)^2)
    @param  egrid: optional interpolator from build_fs_energy_grid; when given the
                   band energies are looked up instead of diagonalised
    @return   eig: eigenvalues [Nk, Nb], reshapeable to (mesh+1, mesh+1, Nb).
                   With egrid the columns are the bands egrid['bands'], not all bands.
    """
    e1,e2=plane_frame(bhat)
    u=np.linspace(-.5*L,.5*L,mesh+1,True)
    U,V=np.meshgrid(u,u,indexing='ij')
    kcart=(U[...,None]*e1+V[...,None]*e2+d*np.asarray(bhat,dtype=float)).reshape(-1,3)
    # k_frac,i = a_i . k_cart / 2pi ; H(k) is periodic in k_frac, so no wrapping needed
    kfrac=kcart.dot(np.asarray(avec,dtype=float).T)/(2.*np.pi)
    if egrid is not None:
        return egrid['interp'](np.mod(kfrac,1.))
    eig,_=get_eigs(kfrac,ham_r,S_r,rvec)
    return eig

def encloses_origin(ct: np.ndarray) -> bool:
    """
    @fn encloses_origin()
    @brief Does this closed contour surround the field axis (the origin of the window)?
           Even-odd ray crossing along +u. This is the test that decides whether the
           section through the centre of the Fermi surface is a closed orbit, which a
           distance-to-the-axis threshold cannot do: a large Fermi surface always has
           some cut copy passing fairly close to the axis.
    @param ct: closed contour, (N,2+) array of in-plane coordinates
    @return  : True if the origin lies inside
    """
    x,y=np.asarray(ct)[:,0],np.asarray(ct)[:,1]
    x2,y2=np.roll(x,-1),np.roll(y,-1)
    straddle=(y>0.)!=(y2>0.)
    if not straddle.any():
        return False
    # u of the edge where it crosses y=0, for the straddling edges only
    xs,ys,xs2,ys2=x[straddle],y[straddle],x2[straddle],y2[straddle]
    xc=xs+(xs2-xs)*(-ys)/(ys2-ys)
    return bool((xc>0.).sum()%2==1)

def get_plane_orbits(eig: np.ndarray, mesh: int, L: float, mu: float,
                     bhat: np.ndarray, d: float,
                     band_map: list | None = None) -> tuple[list, int, set, dict, set]:
    """
    @fn get_plane_orbits()
    @brief Extract the closed Fermi-surface orbits of one plane window.
           Areas come out directly in AA^-2 (the grid is Cartesian, Jacobian 1).
    @param    eig: eigenvalues from get_eigs_plane
    @param   mesh: number of intervals per edge
    @param      L: window edge length [AA^-1]
    @param     mu: chemical potential
    @param   bhat: field direction (unit)
    @param      d: plane offset [AA^-1]
    @retval orbits: list of dicts with keys
                    band  : band index
                    area  : orbit area [AA^-2]
                    cen2d : (u,v) centroid in the window [AA^-1]
                    cen3d : Cartesian centroid [AA^-1]
                    d     : plane offset [AA^-1]
    @retval n_open: number of contour pieces that ran into the window edge
                    (open orbits, or a window too small to enclose them)
    @retval crossing: set of band indices that cross mu somewhere in the window
    @retval axial: bands having a closed orbit that encloses the field axis, i.e. whose
                    section through the centre is a genuine closed orbit
    @retval axis_open: bands whose section through the axis is an OPEN orbit -- nothing
                    closed surrounds the axis AND the contour nearest to it is an open
                    piece. The second half matters: a band whose pockets simply sit away
                    from the axis (electron pockets at X, say) also has nothing closed
                    around the axis, but its orbits are perfectly closed and it must not
                    be called open, nor drive the window to grow chasing it.
    """
    import skimage.measure as sk
    e1,e2=plane_frame(bhat)
    n=np.asarray(bhat,dtype=float)
    orbits=[]
    n_open=0
    crossing=set()
    open_near={}
    closed_near={}
    axial=set()
    for i,e in enumerate(eig.T-mu):
        if e.max()*e.min() >= 0.:
            continue
        i=i if band_map is None else band_map[i]
        crossing.add(i)
        for c in sk.find_contours(e.reshape(mesh+1,mesh+1),0.):
            ct=(c-mesh/2)/mesh*L        # pixel index -> (u,v) in AA^-1
            r=float(np.linalg.norm(ct[:,:2],axis=1).min())
            if not np.allclose(ct[0],ct[-1]):
                n_open+=1
                open_near[i]=min(open_near.get(i,np.inf),r)
                continue
            closed_near[i]=min(closed_near.get(i,np.inf),r)
            if encloses_origin(ct):
                axial.add(i)
            cen=ct[:-1].mean(axis=0)    # drop the repeated closing vertex
            orbits.append({'band':i,'area':shoelace_area(ct),
                           'cen2d':cen,'cen3d':cen[0]*e1+cen[1]*e2+d*n,'d':float(d)})
    axis_open={b for b,r in open_near.items()
               if b not in axial and r<closed_near.get(b,np.inf)}
    return orbits,n_open,crossing,axis_open,axial

def find_plane_orbits(rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray, avec: np.ndarray,
                      mu: float, bhat: np.ndarray, d: float, L0: float, dk: float,
                      ngrow: int = 3, mesh_max: int = 1200, egrid: dict | None = None,
                      ngrow_open: int = 2) -> tuple[list, int, float, set]:
    """
    @fn find_plane_orbits()
    @brief Orbits of one plane, widening the window until the section through the centre
           is resolved, and saying which bands are left with an open orbit there.

           The test for "resolved" is whether some closed contour ENCLOSES the field
           axis. Merely having a closed orbit somewhere in the window is not enough:
           a large Fermi surface always leaves small pockets closed out at the rim
           while the big central orbit is still cut in half by the window edge, which
           is what used to stop the growth too early.

           A section that still does not close after ngrow_open doublings (a few zone
           diameters) is an open orbit, not a window artefact: it carries no dHvA
           oscillation, and shows up in magnetoresistance instead. Growth stops there
           rather than chasing it to the mesh cap.
    @param     mu: chemical potential
    @param   bhat: field direction (unit)
    @param      d: plane offset [AA^-1]
    @param     L0: initial window edge length [AA^-1]
    @param     dk: target k-resolution [AA^-1] (mesh follows L to keep it fixed)
    @param  ngrow: how many times the window may be doubled
    @param mesh_max: cap on the grid size per edge
    @param ngrow_open: doublings spent on a central open section before calling it open
    @retval orbits: closed orbits (see get_plane_orbits)
    @retval n_open: pieces still open at the largest window tried
    @retval      L: window edge length actually used
    @retval open_bands: bands whose section through the centre is an open orbit
    """
    L=float(L0)
    orbits,n_open,central=[],0,set()
    for it in range(ngrow+1):
        mesh=int(min(max(np.ceil(L/dk),16),mesh_max))
        eig=get_eigs_plane(rvec,ham_r,S_r,avec,bhat,d,L,mesh,egrid)
        orbits,n_open,crossing,central,axial=get_plane_orbits(
            eig,mesh,L,mu,bhat,d,None if egrid is None else egrid['bands'])
        # grow only for a band whose axial section is unresolved: nothing closed around
        # the axis and an open piece is the thing nearest to it
        if mesh>=mesh_max or it>=ngrow or not central:
            break
        if it>=ngrow_open:
            break            # a few zone diameters and still no closed central orbit
        L*=2.
    return orbits,n_open,L,central

def fs_cartesian_points(rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray,
                        avec: np.ndarray, mu: float, mesh: int = 40,
                        egrid: dict | None = None) -> np.ndarray:
    """
    @fn fs_cartesian_points()
    @brief Coarse cloud of Cartesian k-points lying on the Fermi surface, computed once
           and reused for every field angle to size the slice window and the slice range.
    @param mesh: 3D sampling mesh per axis
    @return    : (N,3) Cartesian k-points [AA^-1]; empty if no band crosses mu
    """
    Nk,klist=gen_klist(mesh,mesh,mesh)
    eig=(get_eigs(klist,ham_r,S_r,rvec)[0] if egrid is None
         else egrid['interp'](np.mod(klist,1.)))
    tol=max(1.0e-3,(eig.max()-eig.min())/mesh)
    sel=(np.abs(eig-mu)<tol).any(axis=1)
    return klist[sel].dot(get_bvec(np.asarray(avec,dtype=float)))

def fs_extent_along(kf: np.ndarray, bhat: np.ndarray, avec: np.ndarray) -> tuple[float, float]:
    """
    @fn fs_extent_along()
    @brief Fermi-surface size relative to a field direction.
    @param     kf: Fermi-surface point cloud from fs_cartesian_points
    @param   bhat: field direction (unit)
    @param   avec: primitive translation vectors as rows, used as a fallback size
    @retval radius: largest in-plane distance from the B_hat axis [AA^-1]
    @retval   dmax: largest |k.B_hat| on the Fermi surface [AA^-1]
    """
    bv=get_bvec(np.asarray(avec,dtype=float))
    if len(kf)==0:   # no Fermi surface sampled: fall back to the cell size
        r=float(np.linalg.norm(np.abs(bv).sum(axis=0)))
        return r,r
    n=np.asarray(bhat,dtype=float)
    kpar=kf.dot(n)
    kperp=np.linalg.norm(kf-np.outer(kpar,n),axis=1)
    return float(kperp.max()),float(np.abs(kpar).max())

def scan_plane_area(rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray, avec: np.ndarray,
                    mu: float, bhat: np.ndarray, d_arr: np.ndarray, L0: float, dk: float,
                    ngrow: int = 3, egrid: dict | None = None) -> tuple[dict, dict, int, set]:
    """
    @fn scan_plane_area()
    @brief Scan the orbit area along the field direction: for every plane offset d,
           take for each band the orbit closest to the B_hat axis.
           (Step 1 reduction -- one representative branch per band; orbits further out
           are the same branches at other offsets, which the d-scan already covers.)
    @param  d_arr: plane offsets to scan [AA^-1]
    @param     L0: initial window edge length [AA^-1]
    @param     dk: target k-resolution [AA^-1]
    @retval S_scan: dict band_idx -> [(d, area[AA^-2]), ...]
    @retval  orbit_cache: dict d -> full orbit list (all bands, all copies)
    @retval n_open: number of contour pieces that never closed
    @retval open_bands: bands whose CENTRAL section (d=0) is an open orbit
    """
    S_scan={}
    orbit_cache={}
    n_open=0
    d_arr=np.asarray(d_arr,dtype=float)
    # Settle the window on a few probe offsets and then hold it fixed. Letting each
    # slice grow on its own changes how many orbit copies the window holds from slice
    # to slice, which breaks the continuity the branch tracking relies on.
    L=float(L0)
    open_bands=set()
    for j,dp in enumerate({d_arr[0],d_arr[len(d_arr)//2],d_arr[-1]} if len(d_arr) else {0.}):
        _,_,Lp,ob=find_plane_orbits(rvec,ham_r,S_r,avec,mu,bhat,float(dp),L0,dk,
                                    ngrow,egrid=egrid)
        L=max(L,Lp)
        # the reported verdict is about the CENTRAL section (d=0); at other offsets a
        # tube has simply shrunk away, which says nothing about the orbit topology
        if len(d_arr) and np.isclose(dp,d_arr[0]) and np.isclose(d_arr[0],0.):
            open_bands=ob
    for d in d_arr:
        orbits,nop,_,_=find_plane_orbits(rvec,ham_r,S_r,avec,mu,bhat,float(d),L,dk,
                                         0,egrid=egrid)
        n_open+=nop
        orbit_cache[float(d)]=orbits
        for b in {o['band'] for o in orbits}:
            o=min((o for o in orbits if o['band']==b),
                  key=lambda o: np.linalg.norm(o['cen2d']))
            S_scan.setdefault(b,[]).append((float(d),o['area']))
    return S_scan,orbit_cache,n_open,open_bands

def plane_area_at(rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray, avec: np.ndarray,
                  mu: float, bhat: np.ndarray, d: float, band_idx: int,
                  L0: float, dk: float, ngrow: int = 3,
                  egrid: dict | None = None) -> float:
    """
    @fn plane_area_at()
    @brief Area of the representative orbit of one band at one plane offset,
           0 when the band has no closed orbit there (so extremum refinement
           walks away from it rather than crashing).
    """
    orbits=find_plane_orbits(rvec,ham_r,S_r,avec,mu,bhat,float(d),L0,dk,ngrow,egrid=egrid)[0]
    cand=[o for o in orbits if o['band']==band_idx]
    if not cand:
        return 0.
    return min(cand,key=lambda o: np.linalg.norm(o['cen2d']))['area']

def find_extremal_d(d_arr: np.ndarray, S_arr: np.ndarray, band_idx: int,
                    rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray, avec: np.ndarray,
                    mu: float, bhat: np.ndarray, L0: float, dk: float,
                    ngrow: int = 3) -> list:
    """
    @fn find_extremal_d()
    @brief Plane offsets where dA/dd = 0, by sign-change detection on the coarse scan
           followed by a bounded refinement. d=0 is always extremal by time reversal
           (A(-d)=A(d)), provided the scan reaches it.
    @retval : list of extremal d values [AA^-1]
    """
    from scipy.optimize import minimize_scalar
    d_arr=np.asarray(d_arr,dtype=float)
    S_arr=np.asarray(S_arr,dtype=float)
    cand=[]
    if len(d_arr) and np.isclose(d_arr[0],0.):
        cand.append(0.)
    dS=np.diff(S_arr)
    ds_tol=1.0e-5*max(np.max(np.abs(S_arr)) if len(S_arr) else 1.,1.)
    for i in range(len(dS)-1):
        if (abs(dS[i])>ds_tol and abs(dS[i+1])>ds_tol and dS[i]*dS[i+1]<0):
            sign=-1. if dS[i]>0 else 1.
            res=minimize_scalar(
                lambda d,s=sign: s*plane_area_at(rvec,ham_r,S_r,avec,mu,bhat,d,
                                                 band_idx,L0,dk,ngrow),
                bounds=(d_arr[i],d_arr[i+2]),method='bounded',
                options={'xatol':1e-4,'maxiter':10})
            cand.append(float(res.x))
    step=(d_arr[-1]-d_arr[0])/max(len(d_arr)-1,1) if len(d_arr)>1 else 1.
    uniq=[]
    for d in sorted(cand):
        if not uniq or abs(d-uniq[-1])>0.05*step:
            uniq.append(d)
    return uniq

# --------------------------------------------------------------------------- #
#  Orbit branch tracking (step 2)
#
#  A wide window holds several copies of each orbit family. They are NOT the same
#  orbit repeated: a copy displaced by the in-plane part of a reciprocal vector G
#  is the same branch sampled at offset d - G.B_hat, so its area differs. Picking
#  "the one nearest the axis" therefore hops between branches as d advances.
#  Instead every orbit is followed through the scan, extrema are taken per branch,
#  and only at the end are branches that differ by a whole G collapsed.
# --------------------------------------------------------------------------- #

def slice_period(bvec: np.ndarray, bhat: np.ndarray, nshell: int = 4,
                 rtol: float = 1.0e-5) -> float | None:
    """
    @fn slice_period()
    @brief Repeat distance of the Fermi surface along B: the smallest positive value of
           G.B_hat over reciprocal lattice vectors G. Together with A(-d)=A(d) it makes
           [0, period/2] the fundamental domain of the slice scan, and both of its ends
           extremal planes.
           Returns None when B is incommensurate with the lattice -- then {G.B_hat} is
           dense, there is no repeat and no symmetry-plane extremum, and only interior
           dA/dd=0 orbits exist. That is physics, not a numerical failure.
    @param   bvec: reciprocal lattice vectors as rows [AA^-1]
    @param   bhat: field direction (unit)
    @param nshell: largest denominator tried is 8*nshell
    @param   rtol: relative tolerance for accepting the commensurate case
    @return      : period [AA^-1], or None if incommensurate
    """
    from math import gcd
    beta=np.asarray(bvec,dtype=float).dot(np.asarray(bhat,dtype=float))  # b_i . B_hat
    b0=float(np.abs(beta).max())
    if b0<=0.:
        return None
    # Commensurate means beta/b0 is rational with a small denominator. Testing that
    # directly is far safer than hunting for the smallest |sum h_i beta_i|: a field
    # a hair away from a symmetry direction makes some combination nearly cancel, and
    # that near-zero would masquerade as a very short period.
    for q in range(1,nshell*8+1):
        r=q*beta/b0
        if not np.allclose(r,np.rint(r),atol=rtol*max(1.,float(np.abs(r).max()))):
            continue
        m=[abs(int(v)) for v in np.rint(r) if abs(int(v))>0]
        if not m:
            continue
        g=m[0]
        for v in m[1:]:
            g=gcd(g,v)
        return b0*g/q     # smallest positive G.B_hat
    return None

def _orbit_link_cost(o1: dict, o2: dict) -> float:
    """Dimensionless distance between two orbits of neighbouring slices: centroid
    shift measured in orbit radii, plus relative area change."""
    r=np.sqrt(max(o1['area'],o2['area'])/np.pi)
    if r<=0.:
        return np.inf
    return (np.linalg.norm(o1['cen2d']-o2['cen2d'])/r
            +abs(o1['area']-o2['area'])/max(o1['area'],o2['area']))

def track_orbit_branches(orbit_cache: dict, d_arr: np.ndarray,
                         link_tol: float = 0.8) -> list:
    """
    @fn track_orbit_branches()
    @brief Follow each orbit through the slice scan, linking orbits of neighbouring
           offsets by centroid and area continuity (greedy nearest match per band).
    @param orbit_cache: dict d -> orbit list, from scan_plane_area
    @param      d_arr: slice offsets, in scan order
    @param   link_tol: maximum link cost; above it a new branch is started
    @return          : list of branches, each a list of orbit dicts ordered in d
    """
    branches=[]
    active=[]                       # (branch index, last orbit) still growing
    for d in np.asarray(d_arr,dtype=float):
        cur=orbit_cache.get(float(d),[])
        pairs=sorted(((_orbit_link_cost(o_prev,o),ia,io)
                      for ia,(_,o_prev) in enumerate(active)
                      for io,o in enumerate(cur)
                      if o_prev['band']==o['band']),key=lambda t:t[0])
        used_a,used_o,new_active=set(),set(),[]
        for cost,ia,io in pairs:
            if cost>link_tol or ia in used_a or io in used_o:
                continue
            used_a.add(ia); used_o.add(io)
            bi=active[ia][0]
            branches[bi].append(cur[io])
            new_active.append((bi,cur[io]))
        for io,o in enumerate(cur):  # unmatched orbits start their own branch
            if io not in used_o:
                branches.append([o])
                new_active.append((len(branches)-1,o))
        active=new_active
    return branches

def branch_extrema(branch: list, d0_is_extremal: bool = True,
                   d_end: float | None = None) -> list:
    """
    @fn branch_extrema()
    @brief Extremal orbits of one branch: dA/dd = 0, located by a sign change of the
           area difference and refined by a parabolic fit through the bracketing
           triple (A is quadratic near an extremum, so the vertex is second order).
           d=0 is extremal by time reversal, A(-d)=A(d).
    @param  d_end: end of the fundamental domain (half the slice period) when B is
                   commensurate with the lattice; A is symmetric about it too, so a
                   branch reaching it is extremal there. None for an incommensurate B.
    @return: list of dicts (area, d, cen3d, band) at the extrema
    """
    if len(branch)==1:
        # a lone orbit still counts when it sits on a symmetry plane
        o=branch[0]
        if (d0_is_extremal and np.isclose(o['d'],0.)) or \
           (d_end is not None and np.isclose(o['d'],d_end)):
            return [{'band':o['band'],'area':float(o['area']),'d':float(o['d']),
                     'cen3d':o['cen3d']}]
        return []
    d=np.array([o['d'] for o in branch])
    A=np.array([o['area'] for o in branch])
    out=[]
    if d0_is_extremal and np.isclose(d[0],0.):
        out.append({'band':branch[0]['band'],'area':float(A[0]),'d':0.,
                    'cen3d':branch[0]['cen3d']})
    if d_end is not None and np.isclose(d[-1],d_end):
        out.append({'band':branch[-1]['band'],'area':float(A[-1]),'d':float(d[-1]),
                    'cen3d':branch[-1]['cen3d']})
    dA=np.diff(A)
    tol=1.0e-5*max(A.max(),1.)
    for i in range(len(dA)-1):
        if abs(dA[i])<tol or abs(dA[i+1])<tol or dA[i]*dA[i+1]>=0.:
            continue
        x,y=d[i:i+3],A[i:i+3]
        den=(y[0]-2*y[1]+y[2])
        # vertex of the parabola through the three points (uniform or not)
        if abs(den)<1.0e-14:
            xv,yv=x[1],y[1]
        else:
            c=np.polyfit(x,y,2)
            xv=-0.5*c[1]/c[0]
            xv=float(np.clip(xv,x[0],x[2]))
            yv=float(np.polyval(c,xv))
        w=(xv-x[1])/max(x[2]-x[0],1e-30)
        cen=branch[i+1]['cen3d']+w*(branch[i+2]['cen3d']-branch[i]['cen3d'])
        out.append({'band':branch[i+1]['band'],'area':float(yv),'d':float(xv),'cen3d':cen})
    return out

def follow_orbit(seed: dict, rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray,
                 avec: np.ndarray, mu: float, bhat: np.ndarray, L: float, dk: float,
                 step: float, egrid: dict | None = None, nstep: int = 6,
                 reach: float = 1.5, widen: int = 2) -> dict | None:
    """
    @fn follow_orbit()
    @brief Re-find, at a new field direction, the orbit that was extremal at the
           previous one, and return its extremum there.

           The global branch search has to rediscover everything from scratch at every
           angle, and on a Fermi surface that fragments into hundreds of branches it
           drops orbits it found one angle earlier. An orbit already known is much
           easier to follow: it sits on the same Fermi surface, so its centre barely
           moves in k-space and the plane that carries it is simply the one through
           that centre, d = centre . B_hat. From there a short parabolic search along d
           finds dA/dd = 0 without any branch bookkeeping.
    @param   seed: extremal orbit dict from the previous angle (band, area, cen3d)
    @param      L: sampling window edge length [AA^-1]
    @param   step: offset step of the coarse scan, used as the search increment
    @param  nstep: maximum plane evaluations while bracketing the extremum
    @param  reach: how far from the predicted centre the orbit may have moved,
                   in units of its own radius
    @param  widen: how many times the search increment may be doubled when the whole
                   ladder runs monotonically and the extremum must lie further out
    @return      : extremal orbit dict, or None if the orbit is not there any more
    """
    c=np.asarray(seed['cen3d'],dtype=float)
    n=np.asarray(bhat,dtype=float)
    e1,e2=plane_frame(n)
    r_seed=np.sqrt(max(seed['area'],0.)/np.pi)
    tol=max(reach*r_seed,3.*dk)

    def at(d):
        orbits=find_plane_orbits(rvec,ham_r,S_r,avec,mu,n,float(d),L,dk,0,egrid=egrid)[0]
        p=np.array([(c-d*n).dot(e1),(c-d*n).dot(e2)])
        cand=[o for o in orbits if o['band']==seed['band']
              and np.linalg.norm(o['cen2d']-p)<tol]
        return min(cand,key=lambda o:np.linalg.norm(o['cen2d']-p)) if cand else None

    d0=float(c.dot(n))
    h=max(step,1.0e-4)
    o0=at(d0)
    if o0 is None:
        return None                    # the orbit is gone at this field direction
    if abs(d0)<0.5*h:
        # sits on the central plane, which time reversal makes extremal outright
        return {'band':o0['band'],'area':o0['area'],'d':0.,'cen3d':o0['cen3d']}
    # ladder out to both sides -- a maximum and a minimum are found the same way, and
    # walking uphill only would miss every neck orbit
    lad={0:o0}
    for sgn in (1,-1):
        for j in range(1,nstep+1):
            o=at(d0+sgn*j*h)
            if o is None:
                break
            lad[sgn*j]=o
    js=sorted(lad)
    A=np.array([lad[j]['area'] for j in js])
    dd=np.array([d0+j*h for j in js])
    if len(js)<3:
        return None
    dA=np.diff(A)
    best=None
    for i in range(len(dA)-1):
        if dA[i]*dA[i+1]>=0.:
            continue
        cand=0.5*(dd[i+1]+dd[i+1])
        if best is None or abs(dd[i+1]-d0)<abs(best[0]-d0):
            x,y=dd[i:i+3],A[i:i+3]
            co=np.polyfit(x,y,2)
            xv=float(np.clip(-0.5*co[1]/co[0],x[0],x[2])) if abs(co[0])>1e-30 else float(x[1])
            best=(xv,i+1)
    if best is None:
        if widen>0 and len(lad)>=2*nstep+1:
            # the whole ladder ran uphill or downhill: the extremum is further out
            return follow_orbit(seed,rvec,ham_r,S_r,avec,mu,bhat,L,dk,2.*h,egrid,
                                nstep,reach,widen-1)
        return None                    # monotonic: this orbit's extremum is elsewhere
    o=at(best[0])
    if o is None:
        o=lad[js[best[1]]]
        return {'band':o['band'],'area':o['area'],'d':float(dd[best[1]]),'cen3d':o['cen3d']}
    return {'band':o['band'],'area':o['area'],'d':float(best[0]),'cen3d':o['cen3d']}

def verify_orbits_window(records: list, rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray,
                         avec: np.ndarray, mu: float, bhat: np.ndarray, L: float, dk: float,
                         egrid: dict | None = None, grow: float = 1.6,
                         rtol: float = 5.0e-3) -> tuple[list, list]:
    """
    @fn verify_orbits_window()
    @brief Keep only orbits whose area does not move when the sampling window is widened.
           A window that cuts through a connected Fermi-surface region returns the
           boundary of whatever part of it fits, which is a closed curve but not an
           orbit: its area tracks the window. A real orbit does not.
    @param records: extremal orbit dicts (band, area, d, cen3d)
    @param       L: window edge length the records were found with [AA^-1]
    @param    grow: factor by which the window is widened for the check
    @param    rtol: relative area agreement required
    @retval  good: records confirmed window-independent
    @retval  bad : records rejected
    """
    good,bad=[],[]
    for r in records:
        orbits=find_plane_orbits(rvec,ham_r,S_r,avec,mu,bhat,r['d'],grow*L,dk,
                                 0,egrid=egrid)[0]
        same=[o for o in orbits
              if o['band']==r['band']
              and abs(o['area']-r['area'])<=rtol*max(r['area'],1e-30)]
        (good if same else bad).append(r)
    return good,bad

def orbit_key(record: dict, bvec_inv: np.ndarray) -> np.ndarray:
    """
    @fn orbit_key()
    @brief The orbit's 3D centre reduced into the first Brillouin zone, in fractional
           reciprocal coordinates. Two records of the SAME orbit seen in different
           zones differ by exactly a reciprocal lattice vector, so their keys coincide;
           two different extrema of one Fermi-surface tube sit at different offsets d
           and therefore keep different keys however close their areas are.
    @param   record: extremal orbit dict carrying 'cen3d'
    @param bvec_inv: inverse of the reciprocal lattice matrix (rows b_i)
    @return         : reduced fractional centre, each component in [-0.5, 0.5]
    """
    f=np.asarray(record['cen3d'],dtype=float).dot(bvec_inv)
    return f-np.rint(f)

def dedup_extremal_orbits(records: list, bvec: np.ndarray, tol_k: float = 0.02,
                          rtol_area: float = 2.0e-4) -> list:
    """
    @fn dedup_extremal_orbits()
    @brief Reduce the raw extremal orbits to one record per distinct dHvA frequency,
           in two stages that answer two different questions.

           Stage 1 -- is this the same orbit in another zone? The scan deliberately
           covers several repeats of the Fermi surface along B (that redundancy is what
           makes the orbit search robust), so one orbit turns up many times. At its own
           extremum a copy sits at centre + G exactly, so comparing reduced centres
           settles it with no reference to the area. Comparing areas here is what used
           to destroy genuinely distinct extrema of one tube whose areas happen to lie
           within a fraction of a percent of each other.

           Stage 2 -- are these different places on the Fermi surface with the same
           frequency? Point-group images (the X and Y electron pockets, the four
           corner pockets) are not related by any lattice vector, yet they oscillate at
           one frequency and must appear once. Their areas agree to numerical noise,
           which is why the tolerance here is tight.
    @param   records: extremal orbit dicts (band, area, d, cen3d)
    @param      bvec: reciprocal lattice vectors as rows [AA^-1]
    @param     tol_k: how close two reduced centres must be to be one orbit [AA^-1];
                      about one step of the offset scan is right -- above the scatter
                      between copies, below the spacing of distinct extrema
    @param rtol_area: relative area agreement for stage 2
    @return          : one record per distinct orbit, largest area first
    """
    inv=np.linalg.inv(np.asarray(bvec,dtype=float))
    keep=[]
    for r in sorted(records,key=lambda r:-r['area']):
        f=orbit_key(r,inv)
        dup=False
        for k in keep:
            if k['band']!=r['band']:
                continue
            g=f-k['_key']
            if np.linalg.norm((g-np.rint(g)).dot(bvec))<tol_k:
                dup=True
                break
        if not dup:
            r=dict(r); r['_key']=f
            keep.append(r)
    merged=[]
    for r in keep:
        if not any(k['band']==r['band']
                   and abs(k['area']-r['area'])<=rtol_area*max(k['area'],1e-30)
                   for k in merged):
            merged.append(r)
    for r in merged:
        r.pop('_key',None)
    return merged

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
            if S is None:   # open contour: orbit is truncated by the slice box, skip it
                continue
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
    # k_z=0 is fixed by time reversal (S(-kz)=S(kz)) and is therefore always extremal --
    # but only if the scan actually reaches it; a band whose FS appears at finite kz has
    # kz_arr[0]>0, which is an ordinary point.
    # k_z=0.5 is a symmetry plane iff rotmat[2] (=M^T e_z, the slice stacking vector in
    # the fractional basis) is a reciprocal lattice vector: then S(0.5+d)=S(0.5-d) follows
    # from time reversal plus periodicity. For a general tilt it is a spurious orbit.
    m3=np.asarray(rotmat,dtype=float)[2]
    cand_kz=[]
    if np.isclose(kz_arr[0],0.):
        cand_kz.append(kz_arr[0])
    if np.isclose(kz_arr[-1],0.5) and np.allclose(m3,np.rint(m3)):
        cand_kz.append(kz_arr[-1])
    ds_tol=1.0e-5*max(np.max(np.abs(S_arr)),1.0)
    for i in range(len(dS)-1):
        if (abs(dS[i]) > ds_tol and abs(dS[i+1]) > ds_tol
                and dS[i]*dS[i+1]<0):  # sign change in dS/dkz → local extremum in this bracket
            sign=-1. if dS[i]>0 else 1.
            res=minimize_scalar(lambda kz,s=sign: s*S_at_kz(kz),
                                bounds=(kz_arr[i],kz_arr[i+2]),method='bounded',
                                options={'xatol':1e-4,'maxiter':10})
            cand_kz.append(res.x)
    # drop duplicates (a refined extremum can land on the kz=0 / kz=0.5 candidate)
    uniq=[]
    for kz in sorted(cand_kz):
        if not uniq or abs(kz-uniq[-1])>1.0e-3:
            uniq.append(kz)
    return uniq

def gen_3d_surf_points(mesh: int, rvec: np.ndarray, ham_r: np.ndarray,
                       S_r: np.ndarray, mu: float, kscale: float = 1.0,
                       bvec: np.ndarray | None = None) -> tuple[list, list, list]:
    """
    @fn gen_3d_surf_points()
    @brief This function obtains 3d fermi wave-numbers
    @param   mesh: k-mesh
    @param   rvec: r vector of hoppings
    @param  ham_r: hopping parameters
    @param    S_r: overlap integrals
    @param     mu: chemical potential
    @param kscale: change considering k-space area 1.0 is only 1st BZ
    @param   bvec: reciprocal lattice vectors as rows [Angstrom^-1]. If given, the
                   surface vertices are returned in Cartesian k (same frame as BZedge);
                   otherwise they stay in 2pi*fractional units of the reciprocal cell.
    @retval fspolys: triangles of each Fermi sheet (Cartesian if bvec is given)
    @retval fscenters: triangle centers in FRACTIONAL coordinates (always), because
                   they are fed back to get_eigs/get_colors, which take fractional k
    @retval fsband: band index of each sheet
    """
    import skimage.measure as ski
    Nk,klist=gen_klist(mesh+1,mesh+1,mesh+1)
    klist=klist*kscale
    eig,uni=get_eigs(klist,ham_r,S_r,rvec)
    # accept a scalar (int or float) as well as a per-axis [kx,ky,kz] list
    ks=np.asarray(kscale,dtype=float)*np.ones(3)
    fspolys=[]
    fscenters=[]
    fsband=[]
    for i,e in enumerate(eig.T-mu):
        if(e.max()*e.min()<0.):
            verts,faces, _, _=ski.marching_cubes(e.reshape(mesh+1,mesh+1,mesh+1),0,
                                                 spacing=(ks[0]*2*np.pi/mesh,ks[1]*2*np.pi/mesh,
                                                          ks[2]*2*np.pi/mesh))
            verts=verts-ks*np.pi
            # k_frac in [-ks/2, ks/2]; klist fed to get_eigs uses the same fractional units
            kfrac=verts*.5/np.pi
            # Cartesian k = sum_i k_frac_i * b_i, so that the sheet and BZedge(bvec) share a frame
            fspolys.append((kfrac.dot(bvec) if bvec is not None else verts)[faces])
            fscenters.append(kfrac[faces].mean(axis=1))
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
        # eig is needed alongside uni: for a non-orthogonal (MLO) basis the velocity
        # carries the -eps*dS/dk term (get_veloc dispatches on S_r).
        if sw_2d:
            eu=[[get_eigs(k,ham_r,S_r,rvec) for k in kk] for kk in klist]
            vk=[[flibs.get_veloc(k,ham_r,rvec,mrot,unkk,S_r,eigk)[:,b,:]
                 for k,(eigk,unkk) in zip(kk,euk)]
                for kk,euk,b in zip(klist,eu,blist)]
            clist=[[np.sqrt((abs(vkkk)*abs(vkkk)).sum(axis=1)) for vkkk in vkk] for vkk in vk]
        else:
            eu=[get_eigs(k,ham_r,S_r,rvec) for k in klist]
            vk=[flibs.get_veloc(kk,ham_r,rvec,mrot,unk,S_r,eigk)[:,b,:]
                for kk,(eigk,unk),b in zip(klist,eu,blist)]
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
            vk=flibs.get_vnmk(klist,ham_r,rvec,avec,uni,S_r,eig)
            return Nk,eig,vk,kweight
        else:
            vk=flibs.get_veloc(klist,ham_r,rvec,avec,uni,S_r,eig)
            if sw_mass:
                imass=flibs.get_mass(klist,ham_r,rvec,avec,uni,True,S_r,eig)
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
