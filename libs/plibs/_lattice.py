#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
Lattice and reciprocal space utilities: hopping import, k/r-mesh generation, BZ visualization.
"""
import os
import numpy as np
import scipy.linalg as sclin

# Public API of this submodule; keeps `from ._x import *` in __init__.py from
# re-exporting numpy, scipy, ctypes and the sibling packages imported above.
__all__ = [
    'RY_TO_EV', 'import_hoppings', 'import_MLO_hoppings', 'read_HamRsMLO', 'normalize_mlo_orbitals', 'get_bvec', 'gen_rlist',
    'gen_klist_with_kmap', 'gen_klist', 'mk_klist', 'mk_qlist', 'get_ptv', 'get_symm_line',
    'bz_faces', 'bz_face_planes', 'fold_to_bz'
]

RY_TO_EV = 13.6          # Rydberg -> eV conversion for the MLO (ecalj) importers


def import_hoppings(fname: str, ftype: int) -> tuple[np.ndarray, np.ndarray, int, int]:
    """
    @brief This function imports hopping parameters from files
    @param fname: Name of import files
    @param ftype: File format of import files
    @return  rvec: the array of r-vector, size:(nr,3)
    @return ham_r: the array of hopping integrals, size:(nr,no,no)
    @return    no: the number of orbitals
    @return    nr: the number of r-vectors
    """
    def import_hop(name: str) -> tuple:
        rvec = np.loadtxt(f'{name}/irvec.txt')
        nr = rvec[:, 0].size
        ndegen = np.loadtxt(f'{name}/ndegen.txt')
        tmp = np.array([complex(float(tp[0]), float(tp[1])) for tp in
                        [f.strip(' ()\n').split(',') for f in open(f'{name}/ham_r.txt', 'r')]])
        no = int(np.sqrt(tmp.size/nr))
        # Divide each R-block by its degeneracy weight (number of equivalent k-points)
        ham_r = (tmp.reshape(nr, no, no).T/ndegen).T.round(6).copy()
        return (rvec, ham_r, no, nr)

    def import_out(name: str) -> tuple:
        data = np.loadtxt(name)
        # Count rows that share the same R-vector as the first row → gives no^2 (orbital pairs per R)
        con = (data[:, :3] == data[0, :3]).prod(axis=1).sum()
        no, nr = int(np.sqrt(con)), data[:, 0].size//con
        rvec = np.array(data[:nr, :3])
        ham_r = (data[:, 3]+1j*data[:, 4]).reshape(no*no, nr).T.reshape(nr, no, no).round(6).copy()
        return (rvec, ham_r, no, nr)

    def import_hr(name: str) -> tuple:
        tmp = [f.split() for f in open(f'{name}_hr.dat', 'r')]
        no, nr = int(tmp[1][0]), int(tmp[2][0])
        c2, tmp1 = 3, []
        while not len(tmp1) == nr:
            tmp1.extend(tmp[c2])
            c2 = c2+1
        ndegen = np.array([float(t) for t in tmp1])
        tmp1 = [[float(t) for t in tp] for tp in tmp[c2:]]
        tmp = np.array([complex(tp[5], tp[6]) for tp in tmp1])
        rvec = np.array([tmp1[no*no*i][:3] for i in range(nr)])
        ham_r = (tmp.reshape((nr, no, no)).T/ndegen).T.round(6).copy()
        return (rvec, ham_r, no, nr)

    def import_Hopping(name: str, sw_axis: bool = False) -> tuple:
        tmp = [f.split() for f in open(f'{name}/Hopping.dat', 'r')]
        axis = np.array([[float(tp) for tp in tpp] for tpp in tmp[1:4]])
        no, nr = int(tmp[4][0]), int(tmp[4][1])
        tmp1 = np.array([[float(t) for t in tp] for tp in tmp[7+no:]])
        rvec = np.array([tmp1[no*no*i][:3] for i in range(nr)])
        tmp = np.array([complex(tp[8], tp[9]) for tp in tmp1])
        ham_r = tmp.reshape(nr, no, no)
        if sw_axis:
            return (rvec, ham_r, no, nr, axis)
        else:
            return (rvec, ham_r, no, nr)
    if ftype == 0:
        rvec, ham_r, no, nr = import_hop(fname)
    elif ftype == 1:
        rvec, ham_r, no, nr = import_out(fname)
    elif ftype == 2:
        rvec, ham_r, no, nr = import_hr(fname)
    else:
        rvec, ham_r, no, nr = import_Hopping(fname)
    return (rvec, ham_r, no, nr)


def normalize_mlo_orbitals(rvec: np.ndarray, ham_r: np.ndarray,
                           S_r: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    @brief Scale every MLO orbital to unit norm: H(R) -> D^-1/2 H(R) D^-1/2, S(R) -> D^-1/2 S(R) D^-1/2 with
           D = diag S(R=0).  The bands are unchanged (the generalized eigenproblem H u = e S u is invariant
           under a diagonal rescaling of the basis), but every quantity in the orbital basis -- the local
           U, J, orbital weights, G, Sigma, the gap -- then refers to normalized orbitals, as for an orthogonal
           Wannier basis.  ecalj's MLOs are not normalized (e.g. Cu: S_aa(0) = 0.016 for s, 0.565 for d), and a
           U on the raw orbitals acts with a strength ~ 1/S_aa(0)^2 instead of U.
    @return ham_r, S_r: the rescaled arrays (new arrays; S_aa(R=0) = 1)
    """
    i0 = np.flatnonzero((np.abs(np.asarray(rvec)) < 1.0e-9).all(axis=1))
    if len(i0) != 1:
        raise ValueError('normalize_mlo_orbitals: R = 0 not found (once) among the R-vectors')
    d = S_r[i0[0]].diagonal().real
    if (d <= 0).any():
        raise ValueError(f'normalize_mlo_orbitals: non-positive orbital norm S_aa(R=0) = {d.min():.3e}')
    w = 1.0/np.sqrt(d)
    scale = w[:, None]*w[None, :]
    return ham_r*scale, S_r*scale


def import_MLO_hoppings(name: str, normalize: bool = True) -> tuple[np.ndarray, np.ndarray, np.ndarray, int, int]:
    """
    @brief This function imports MLO hopping parameters from files
    @param name: File name
    @param normalize: scale the orbitals to unit norm (normalize_mlo_orbitals); False keeps the raw MLOs
    @return  rvec: the array of r-vector, size:(nr,3)
    @return ham_r: the array of hopping integrals, size:(nr,no,no)
    @return   S_r: the array of overlap integrals, size:(nr,no,no)
    @return    no: the number of orbitals
    @return    nr: the number of r-vectors
    """
    tmp = [f.split() for f in open(f'{name}', 'r')]
    tmp1 = np.array([[float(t) for t in tp] for tp in tmp])
    no = int(tmp1[:, 0].max())
    nr = int(len(tmp1)/(no*no))
    tmp = np.array([complex(tp[5], tp[6]) for tp in tmp1])
    tmpS = np.array([complex(tp[7], tp[8]) for tp in tmp1])
    rvec = np.array([tmp1[i][2:5] for i in range(nr)])
    ham_r = tmp.reshape((no*no, nr)).T.reshape((nr, no, no)).round(6).copy()*RY_TO_EV
    S_r = tmpS.reshape((no*no, nr)).T.reshape((nr, no, no)).round(6).copy()
    if normalize:
        ham_r, S_r = normalize_mlo_orbitals(rvec, ham_r, S_r)
    return rvec, ham_r, S_r, no, nr


def read_HamRsMLO(directory: str = '.', normalize: bool = True) -> tuple[np.ndarray, np.ndarray, np.ndarray, int, int]:
    """
    @brief Read ecalj job_mlo binary output (HamRsMLO + HamiltonianPMTInfo).
           Returns the same format as import_MLO_hoppings.
    @param directory: path containing HamRsMLO and HamiltonianPMTInfo
    @param normalize: scale the orbitals to unit norm (normalize_mlo_orbitals); False keeps the raw MLOs
    @return  rvec: R-vectors (nr, 3) float64, in lattice-vector units
    @return ham_r: H(R) (nr, no, no) complex128, eV
    @return   S_r: S(R) (nr, no, no) complex128, dimensionless
    @return    no: number of MLO orbitals (ndimMTO)
    @return    nr: number of R-vectors
    """
    def _read_record(f: np.ndarray) -> bytes:
        raw = f.read(4)
        if len(raw) < 4:
            raise EOFError('Unexpected end of Fortran binary file')
        n = int(np.frombuffer(raw, dtype='<i4')[0])
        data = f.read(n)
        f.read(4)  # trailing record marker
        return data

    # ---- HamiltonianPMTInfo ----
    with open(os.path.join(directory, 'HamiltonianPMTInfo'), 'rb') as f:
        # Record 1: plat(3,3)[f64] nkk1,nkk2,nkk3,nbas[i32] qlat(3,3)[f64]
        rec = _read_record(f)
        nkk1, nkk2, nkk3, nbas = np.frombuffer(rec, dtype='<i4', count=4, offset=72)
        nbas = int(nbas)

        _read_record(f)  # Record 2: pos(3,nbas), alat  — not needed
        _read_record(f)  # Record 3: qplist(3,nkp)      — not needed

        # Record 4: npair(nbas,nbas)[i32] npairmx[i32]
        rec = _read_record(f)
        arr4 = np.frombuffer(rec, dtype='<i4')
        npair    = arr4[:nbas * nbas].reshape(nbas, nbas, order='F').copy()
        npairmx  = int(arr4[nbas * nbas])

        # Record 5: nlat(3,npairmx,nbas,nbas)[i32] nqwgt(npairmx,nbas,nbas)[i32]
        rec = _read_record(f)
        n_nlat = 3 * npairmx * nbas * nbas
        nlat = np.frombuffer(rec, dtype='<i4', count=n_nlat) \
                 .reshape(3, npairmx, nbas, nbas, order='F').copy()
        # nqwgt follows but is not needed

    # ---- HamRsMLO ----
    with open(os.path.join(directory, 'HamRsMLO'), 'rb') as f:
        # Record 1: ndimMTO, npairmx, nspx
        rec = _read_record(f)
        ndimMTO, _, nspx = np.frombuffer(rec, dtype='<i4')
        ndimMTO, nspx = int(ndimMTO), int(nspx)
        shape4 = (npairmx, ndimMTO, ndimMTO, nspx)
        expected_bytes = npairmx * ndimMTO * ndimMTO * nspx * 16  # complex128

        # Record 2: hammr(npairmx, ndimMTO, ndimMTO, nspx)
        rec = _read_record(f)
        hammr = np.frombuffer(rec, dtype='<c16').reshape(shape4, order='F').copy()

        # Record 3: ovlmr  — or hammhsor if SOC was enabled (shape uses 3 instead of nspx)
        rec = _read_record(f)
        if len(rec) != expected_bytes:
            # It is hammhsor (SOC); skip it and read ovlmr next
            rec = _read_record(f)
        ovlmr = np.frombuffer(rec, dtype='<c16').reshape(shape4, order='F').copy()

        # Record 4/5: ib_tableM(ndimMTO), k_tableM(ndimMTO), l_tableM(ndimMTO)
        rec = _read_record(f)
        tables   = np.frombuffer(rec, dtype='<i4').reshape(3, ndimMTO)  # C-order: row = table
        ib_tableM = tables[0] - 1  # Fortran 1-indexed -> 0-indexed atom index

    # ---- Build (rvec, ham_r, S_r) ----
    # hammr[it, i, j, sp] has R-vector nlat[:, it, ib_tableM[i], ib_tableM[j]].
    # Different (i,j) with different atom pairs use different nlat sub-arrays,
    # so we collect all unique R-vectors across all atom pairs.

    rvec_list: list[tuple[int, int, int]] = []
    rvec_idx:  dict[tuple[int, int, int], int] = {}
    rlookup:   dict[tuple[int, int], dict] = {}  # (ib1,ib2) -> {R_tuple -> it}

    for ib1 in range(nbas):
        for ib2 in range(nbas):
            nr_pair = int(npair[ib1, ib2])
            d: dict[tuple[int, int, int], int] = {}
            for it in range(nr_pair):
                R = (int(nlat[0, it, ib1, ib2]),
                     int(nlat[1, it, ib1, ib2]),
                     int(nlat[2, it, ib1, ib2]))
                d[R] = it
                if R not in rvec_idx:
                    rvec_idx[R] = len(rvec_list)
                    rvec_list.append(R)
            rlookup[(ib1, ib2)] = d

    nr = len(rvec_list)
    no = ndimMTO
    ham_r = np.zeros((nr, no, no), dtype=np.complex128)
    S_r   = np.zeros((nr, no, no), dtype=np.complex128)

    for i in range(no):
        ib1 = int(ib_tableM[i])
        for j in range(no):
            ib2 = int(ib_tableM[j])
            for R, it in rlookup[(ib1, ib2)].items():
                ir = rvec_idx[R]
                ham_r[ir, i, j] = hammr[it, i, j, 0]
                S_r[ir, i, j]   = ovlmr[it, i, j, 0]

    ham_r *= RY_TO_EV
    # float64 like the other importers: the Fortran kernels declare rvec as double
    rvec = np.array(rvec_list, dtype=np.float64)
    if normalize:
        ham_r, S_r = normalize_mlo_orbitals(rvec, ham_r, S_r)
    return rvec, ham_r, S_r, no, nr


def get_bvec(avec: np.ndarray) -> np.ndarray:
    """
    @brief This function generates reciprocal lattice vector from primitive translation vector
    @param  avec: primitive translation vector
    @return bvec: reciprocal lattice vector
    """
    # b_i = 2pi * (a^{-T})_i  (standard reciprocal lattice relation a_i . b_j = 2pi delta_ij)
    bvec = 2*np.pi*sclin.inv(avec).T
    return bvec


def gen_rlist(Nx: int, Ny: int, Nz: int) -> np.ndarray:
    """
    @brief This function generates the list of r-vectors
    @param     Nx: Number of sites along x-axis
    @param     Ny: Number of sites along y-axis
    @param     Nz: Number of sites along z-axis
    @return rlist: The array of r-vectors
    """
    x0 = np.linspace(0, Nx, Nx, False)
    y0 = np.linspace(0, Ny, Ny, False)
    z0 = np.linspace(0, Nz, Nz, False)
    x, y, z = np.meshgrid(x0, y0, z0)
    rlist = np.array([x.ravel(), y.ravel(), z.ravel()]).T.copy()

    return rlist


def gen_klist_with_kmap(Nx: int, Ny: int, Nz: int) -> tuple[np.ndarray, np.ndarray]:
    """
    @brief This function generates the list of k-vectors
    @param     Nx: Number of kx mesh
    @param     Ny: Number of ky mesh
    @param     Nz: Number of kz mesh
    @return klist: The array of k-vectors
    @return  kmap: The array of property of k-points
    """
    x0 = np.linspace(0, Nx, Nx, False, dtype=int)
    y0 = np.linspace(0, Ny, Ny, False, dtype=int)
    z0 = np.linspace(0, Nz, Nz, False, dtype=int)
    x, y, z = np.meshgrid(x0, y0, z0)
    kmap = np.array([x.ravel(), y.ravel(), z.ravel()]).T.copy()
    klist = np.array([x.ravel()/Nx, y.ravel()/Ny, z.ravel()/Nz]).T.copy()
    return klist, kmap


def gen_klist(Nx: int, Ny: int, Nz: int | None = None, sw_pp: bool = True, kz: float = 0.0) -> tuple[int, np.ndarray]:
    """
    @brief This function generates a k-point list meshed in Nx,Ny,Nz
    @param     Nx: Number of axis 1 mesh (usually kx mesh)
    @param     Ny: Number of axis 2 mesh (usually ky mesh)
    @param     Nz: Number of axis 3 mesh (usually kz mesh)
    @param  sw_pp: switch output 2D mesh or 3D mesh
    @param     kz: value of axis 3 at 2D mesh (use only if sw_pp is True)
    @return    Nk: number of k-points
    @return klist: list of k-points
    """
    if sw_pp:
        # Symmetric BZ: k in [-0.5, 0.5] for band/FS plots (includes both zone boundaries)
        kx = np.linspace(-0.5, 0.5, Nx, True)
        ky = np.linspace(-0.5, 0.5, Ny, True)
        if Nz is None:
            kz = np.array([kz])
        else:
            kz = np.linspace(-0.5, 0.5, Nz, True)
    else:
        # Periodic BZ: k in [0, 1) without endpoint for FFT-compatible k-summations
        kx = np.linspace(0, 1, Nx, False)
        ky = np.linspace(0, 1, Ny, False)
        kz = np.linspace(0, 1, Nz, False)
    # indexing='ij' so that eig.reshape(Nx,Ny[,Nz]) has axes ordered (kx,ky[,kz]);
    # the default 'xy' would swap the kx/ky axes in FS contour/marching-cubes plots
    x, y, z = np.meshgrid(kx, ky, kz, indexing='ij')
    klist = np.array([x.ravel(), y.ravel(), z.ravel()]).T.copy()
    return len(klist), klist


def mk_klist(k_list: np.ndarray | list, N: int, bvec: np.ndarray) -> tuple[np.ndarray, np.ndarray, list]:
    """
    @brief This function generates k-point list for symmetry line
    @param  k_list: the list of symmetry points
    @param       N: Number of mesh between symmetry points
    @param    bvec: reciprocal lattice vector
    @return  klist: list of k-points at band plot
    @return  splen: length of symmetry points
    @return xticks: footnotes of klist at ticks of symmetry points
    """
    klist = []
    splen = []
    xticks = []
    maxsplen = 0
    for ks, ke in zip(k_list, k_list[1:]):
        dkv = np.array(ke)-np.array(ks)
        # Arc length of segment in Cartesian reciprocal space (Angstrom^-1)
        dkv_length = np.linalg.norm(dkv.dot(bvec))
        tmp = np.linspace(ks, ke, N, False)       # N points, excluding endpoint (appended at end)
        tmp2 = np.linspace(0, dkv_length, N, False)+maxsplen
        # advance by the full segment length so consecutive segments join seamlessly
        # (tmp2.max() would be one step short and shift every following xtick)
        maxsplen += dkv_length
        xticks += [tmp2[0]]
        klist += tmp.tolist()
        splen += tmp2.tolist()
    klist += [k_list[-1]]        # append the final symmetry point
    splen += [maxsplen]
    xticks += [maxsplen]
    return np.array(klist), np.array(splen), xticks


def mk_qlist(k_set: np.ndarray | list, Nx: int, Ny: int, Nz: int,
             bvec: np.ndarray) -> tuple[np.ndarray, np.ndarray, list]:
    """
    @brief This function generates q-point list for symmetry line
    @param   k_set: the list of symmetry points
    @param      Nx: Number of x-mesh
    @param      Ny: Number of y-mesh
    @param      Nz: Number of z-mesh
    @param    bvec: reciprocal lattice vector
    @return  qlist: list of q-points at band plot
    @return  splen: length of symmetry points
    @return xticks: footnotes of qlist at ticks of symmetry points
    """
    qlist = []
    splen = []
    xticks = []
    maxsplen = 0
    Narray = np.array([Nx, Ny, Nz])
    dk_length = 0.0
    N = 1
    for ks, ke in zip(k_set, k_set[1:]):
        dk = np.array(ke)-np.array(ks)
        # Arc length of segment in Cartesian reciprocal space (Angstrom^-1)
        dk_length = np.linalg.norm(dk.dot(bvec))
        # ensure at least one division along non-zero component; use ceil to avoid zero due to truncation
        # Compute per-axis step counts; use ceil so non-zero dk always yields ≥1 point
        dN = np.asarray(np.ceil(abs(dk)*Narray), dtype=int)
        nonzero = dN[dN > 0]
        if len(nonzero) == 0:
            dk_length = 0.0  # identical consecutive k-points; skip this segment
            continue
        N = nonzero.min()   # use the tightest constraint to avoid oversampling any axis
        tmp = np.linspace(ks, ke, N, False)
        tmp2 = np.linspace(0, dk_length, N, False)+maxsplen
        maxsplen += dk_length
        xticks += [tmp2[0]]
        qlist += tmp.tolist()
        splen += tmp2.tolist()
    qlist += [k_set[-1]]
    splen += [maxsplen]          # maxsplen already includes the full last segment
    xticks += [maxsplen]
    # get_qshift resolves k+q by an exact lookup on the k-mesh, so every q has to be a mesh
    # point.  The step above is not always commensurate (a 1/3 coordinate such as hexagonal K,
    # or a k_z=1/2 end point on an Nz=1 mesh), so snap each q to the nearest mesh point;
    # the path length (splen) keeps the ideal positions for the plot axis.
    qlist = np.array(qlist, dtype=np.float64)
    qsnap = np.round(qlist*Narray)/Narray
    dev = np.abs(qsnap-qlist).max() if len(qlist) else 0.0
    if dev > 1.0e-8:
        print(f'mk_qlist: q-points snapped to the {Nx}x{Ny}x{Nz} k-mesh '
              f'(max shift {dev:.3g} in reciprocal-lattice units)', flush=True)
    return qsnap, np.array(splen), xticks


def get_ptv(alatt: np.ndarray, deg: np.ndarray, brav: int) -> tuple[np.ndarray, np.ndarray]:
    """
    @brief Generate primitive translation vectors (avec) and rotation matrix (Arot) for a given Bravais lattice type.
    @param  alatt: Lattice constants [a, b, c] in Angstrom
    @param    deg: Lattice angles [alpha, beta, gamma] in degrees (used for brav>=8)
    @param   brav: Bravais lattice type index (0=simple cubic, 1=FCC, 2=BCC, 3=hexagonal,
                   4=trigonal, 5=base-centered, 6=FCC2, 7=BCC2, >=8=general triclinic)
    @return  avec: Primitive translation vectors as rows (shape [3,3])
    @return  Arot: Rotation matrix defining the primitive cell in terms of conventional cell
    """
    if brav == 0:  # simple
        Arot = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    elif brav == 1:  # face center
        Arot = np.array([[-.5, 0., .5], [0., .5, .5], [-.5, .5, 0.]])
    elif brav == 2:  # body center
        Arot = np.array([[.5, -.5, .5], [.5, .5, .5], [-.5, -.5, .5]])
    elif brav == 3:  # hexagonal
        Arot = np.array([[1., 0., 0.], [-.5, .5*np.sqrt(3.), 0.], [0., 0., 1.]])
    elif brav == 4:  # trigonal
        cg = np.cos(np.pi*deg[2]/180.)
        # Primitive vectors of rhombohedral cell expressed in orthogonal frame (Ashcroft & Mermin convention)
        tx, ty, tz = np.sqrt((1-cg)*0.5), np.sqrt((1.-cg)/6.), np.sqrt((1.+2*cg)/3.)
        Arot = np.array([[tx,  -ty, tz], [0., 2*ty, tz], [-tx,  -ty, tz]])
    elif brav == 5:  # base center
        Arot = np.array([[.5, .5, 0.], [-.5, .5, 0.], [0., 0., 1.]])
    elif brav == 6:  # face center 2
        # a1,a2 in this order so that det(Arot)>0 (right-handed set)
        Arot = np.array([[0., .5, .5], [.5, 0., .5], [.5, .5, 0.]])
    elif brav == 7:  # body center 2
        Arot = np.array([[-.5, .5, .5], [.5, -.5, .5], [.5, .5, -.5]])
    else:
        # General triclinic lattice: build orthogonal representation from lattice angles (alpha, beta, gamma)
        phase = [np.pi*deg[0]/180., np.pi*deg[1]/180., np.pi*deg[2]/180.]
        ax1, ax2 = alatt[1]/alatt[0], alatt[2]/alatt[0]   # b/a, c/a ratios
        r1, r2, r3 = ax1*np.cos(phase[2]), ax1*np.sin(phase[2]), ax2*np.cos(phase[1])
        sin_gamma = np.sin(phase[2])
        if abs(sin_gamma) < 1.0e-12:
            raise ValueError(f"Invalid gamma angle for triclinic cell: gamma={deg[2]} deg")
        r4 = ax2*(np.cos(phase[0])-np.cos(phase[1])*np.cos(phase[2]))/sin_gamma
        # r5 = c * sqrt(1 - cos^2(alpha) - cos^2(beta) - cos^2(gamma) + 2cos(alpha)cos(beta)cos(gamma)) / sin(gamma)
        r5 = ax2*np.sqrt(1+2*np.cos(phase[0])*np.cos(phase[1])*np.cos(phase[2])
                         - (np.cos(phase[0])**2+np.cos(phase[1])**2+np.cos(phase[2])**2))/sin_gamma
        Arot = np.array([[1., 0.,  0.],
                        [r1, r2,  0.],
                         [r3, r4, r5]])
    if brav in {0, 1, 2, 3, 4, 5, 6, 7}:
        avec = alatt*Arot
    else:
        avec = alatt[0]*Arot
    # Cell volume must be positive: several routines use det(avec) directly as the
    # unit-cell volume (carrier density, conductivity), so a left-handed set would
    # silently flip their sign.
    vol = sclin.det(avec)
    if vol <= 0.:
        raise ValueError(f"primitive translation vectors for brav={brav} are not right-handed "
                         f"(cell volume = {vol})")
    return avec, Arot


def get_symm_line(brav: int) -> tuple[list, list]:
    """
    @brief Return the high-symmetry k-point path and corresponding axis labels for a given Bravais lattice.
    @param    brav: Bravais lattice type index (same convention as get_ptv)
    @return k_list: List of high-symmetry k-points in fractional coordinates
    @return xlabel: List of label strings (LaTeX) for each high-symmetry point
    """
    if brav == 0:  # simple
        k_list = [[0., 0., .5], [0., 0., 0.], [.5, 0., 0.], [.5, .5, 0.], [0., 0., 0.]]
        xlabel = ['Z', r'$\Gamma$', 'X', 'M', r'$\Gamma$']
    elif brav == 1:  # face center
        k_list = [[0., 0., 0.], [.5, 0., .5], [1., 0., 0.], [.5, .5, .5], [.5, .25, .75], [0., 0., 0.]]
        xlabel = [r'$\Gamma$', 'X', r'$\Gamma$', 'L', 'W', r'$\Gamma$']
    elif brav == 2:  # body center
        k_list = [[.5, .5, .5], [0., 0., 0.], [.5, 0., 0.], [.5, .5, -.5], [0., 0., 0.]]
        xlabel = ['Z', r'$\Gamma$', 'X', 'M', r'$\Gamma$']
    elif brav == 3:  # hexagonal
        k_list = [[0., 0., 0.], [2./3., -1./3., 0.], [.5, 0., 0.], [0., 0., 0.], [0., 0., .5]]
        xlabel = [r'$\Gamma$', 'K', 'M', r'$\Gamma$', 'Z']
    elif brav == 4:  # trigonal
        k_list = [[0., 0., 0.], [.5, 0., .5], [.5, 0., 0.], [0., 0., 0.], [.5, .5, .5]]
        xlabel = [r'$\Gamma$', 'K', 'M', r'$\Gamma$', 'Z']
    elif brav == 5:  # base center
        k_list = [[0., 0., 0.], [.5, 0., 0.], [.5, .5, 0.], [0., 0., 0.], [0., 0., .5]]
        xlabel = [r'$\Gamma$', 'X', 'S', r'$\Gamma$', 'Z']
    elif brav == 6:  # face center 2
        k_list = [[0., 0., 0.], [.5, 0., .5], [.5, .5, .5], [.75, .25, .5], [0., 0., 0.]]
        xlabel = [r'$\Gamma$', 'X', 'L', 'W', r'$\Gamma$']
    elif brav == 7:  # body center (common)
        k_list = [[.5, .5, .5], [0., 0., 0.], [.5, 0., 0.], [.5, .5, -.5], [0., 0., 0.]]
        xlabel = ['Z', r'$\Gamma$', 'X', 'M', r'$\Gamma$']
    else:
        k_list = [[0., 0., 0.], [.5, 0., 0.], [.5, .5, 0.], [0., 0., 0.]]
        xlabel = [r'$\Gamma$', 'X', 'M', r'$\Gamma$']
    return k_list, xlabel


def _order_bz_polygon(verts: np.ndarray) -> np.ndarray:
    """
    @brief Sort the vertices of a planar polygon by azimuth around its centre, so that
           consecutive vertices are the edges of the polygon.
    @param verts: polygon vertices [Nv,3] in arbitrary order
    @return     : the same vertices, ordered along the polygon
    """
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


def bz_faces(bvec: np.ndarray, nshell: int = 2) -> list:
    """
    @brief Faces of the 1st Brillouin zone, i.e. of the Wigner-Seitz cell of the
           reciprocal lattice, from the Voronoi cell of the origin.
    @param   bvec: reciprocal lattice vectors as rows, shape (3,3), in Angstrom^-1
    @param nshell: shell of reciprocal lattice points used for the Voronoi construction
    @return      : list of faces, each an ordered vertex array [Nv,3] in Cartesian k
    """
    import itertools
    import scipy.spatial as scspatial

    bv = np.asarray(bvec, dtype=float)
    idx = range(-nshell, nshell+1)
    pts = np.array([i*bv[0] + j*bv[1] + k*bv[2] for i, j, k in itertools.product(idx, repeat=3)])
    vor = scspatial.Voronoi(pts)
    origin_idx = np.argmin(np.linalg.norm(pts, axis=1))
    faces = []
    for ridge_pts, ridge_verts in zip(vor.ridge_points, vor.ridge_vertices):
        if origin_idx not in ridge_pts or -1 in ridge_verts:
            continue
        faces.append(_order_bz_polygon(vor.vertices[ridge_verts]))
    return faces


def bz_face_planes(bvec: np.ndarray, nshell: int = 2) -> tuple[np.ndarray, np.ndarray]:
    """
    @brief Half-space representation of the 1st BZ: k is inside iff k.n_i <= d_i for every
           face i. Face i is the perpendicular bisector of a reciprocal lattice vector G_i,
           so n_i = G_i/|G_i| and d_i = |G_i|/2; G is a face iff its own midpoint G/2 is not
           already cut off by the bisector of another G'.
    @param   bvec: reciprocal lattice vectors as rows, shape (3,3), in Angstrom^-1
    @param nshell: shell of reciprocal lattice vectors scanned for zone faces
    @return  nvec: outward unit normals of the zone faces [Nface,3]
    @return  dist: distance of each face from Gamma [Nface]
    """
    import itertools

    bv = np.asarray(bvec, dtype=float)
    gv = np.array([g for g in itertools.product(range(-nshell, nshell+1), repeat=3)
                   if g != (0, 0, 0)], dtype=float).dot(bv)
    half = .5*(gv*gv).sum(axis=1)                        # |G|^2/2 = G.(G/2)
    tol = 1.e-10*half.max()
    proj = (.5*gv).dot(gv.T)
    np.fill_diagonal(proj, -np.inf)                      # G.(G/2) = |G|^2/2 is the plane itself
    # strict inequality: a bisector that only touches the cell along an edge or at a corner
    # (G=(1,1,0) for a cubic lattice, say) cuts nothing off and is not a face
    isface = (proj < half[None, :]-tol).all(axis=1)
    faces = gv[isface]
    dist = .5*np.linalg.norm(faces, axis=1)
    return faces/(2.*dist[:, None]), dist


def _clip_half_space(poly: np.ndarray, nvec: np.ndarray, dist: float, eps: float) -> np.ndarray | None:
    """
    @brief Sutherland-Hodgman clip of a convex polygon by the half space k.n <= d.
    @param  poly: polygon vertices [Nv,3], ordered along the polygon
    @param  nvec: unit normal of the clipping plane
    @param  dist: distance of the clipping plane from the origin
    @param   eps: tolerance on the signed distance
    @return     : the clipped polygon, or None if nothing (or only a sliver) is left
    """
    s = poly.dot(nvec)-dist
    if (s <= eps).all():
        return poly
    if (s > eps).all():
        return None
    out = []
    for i in range(len(poly)):
        si, sj = s[i], s[(i+1) % len(poly)]
        if si <= eps:
            out.append(poly[i])
        if (si < -eps and sj > eps) or (si > eps and sj < -eps):   # the edge crosses the plane
            out.append(poly[i] + si/(si-sj)*(poly[(i+1) % len(poly)]-poly[i]))
    return np.array(out) if len(out) >= 3 else None


def fold_to_bz(polys: np.ndarray, bvec: np.ndarray, nshell: int = 2) -> list:
    """
    @brief Fold triangles sampled on the primitive reciprocal cell -- the parallelepiped
           k_frac in [-1/2,1/2]^3 -- into the 1st BZ, the Wigner-Seitz cell drawn by libs.plots.BZedge.
           For a primitive (P) lattice the two cells coincide and the triangles come back
           untouched, but for a centred lattice (body/face/base centred, and any skewed
           cell) they do not: the part of the sheet outside the zone has to be translated by
           a reciprocal lattice vector and cut at the zone faces, otherwise the Fermi surface
           sticks out of the BZ frame.
           The translates of the parallelepiped tile k-space, so every interior point of the
           zone is covered exactly once and no part of the sheet is duplicated or lost.
    @param  polys: triangles [Ntri,3,3] in Cartesian k [Angstrom^-1]
    @param   bvec: reciprocal lattice vectors as rows, shape (3,3), in Angstrom^-1
    @param nshell: shell of reciprocal lattice vectors scanned for zone faces
    @return      : list of convex polygons [Nv,3] (Nv>=3), all inside the zone
    """
    import itertools

    tris = np.asarray(polys, dtype=float)
    if tris.size == 0:
        return list(tris)
    bv = np.asarray(bvec, dtype=float)
    nvec, dist = bz_face_planes(bv, nshell)
    eps = 1.e-9*dist.max()
    if (tris.dot(nvec.T)-dist <= eps).all():
        # Nothing sticks out, and then nothing is waiting outside to be folded in either:
        # what would fill the part of the zone the parallelepiped misses is the translate
        # of the sheet in the part of the parallelepiped outside the zone, which is empty.
        return list(tris)
    # a point at fractional k is covered by the translate G = rint(k_frac), and the
    # per-component extremes of the (convex) zone are reached at its corners
    corners = np.vstack(bz_faces(bv, nshell))
    nmax = max(1, int(np.ceil(np.abs(corners.dot(sclin.inv(bv))).max())))
    tlo, thi = tris.reshape(-1, 3).min(axis=0), tris.reshape(-1, 3).max(axis=0)
    zlo, zhi = corners.min(axis=0), corners.max(axis=0)
    out = []
    for gfrac in itertools.product(range(-nmax, nmax+1), repeat=3):
        g = np.array(gfrac, dtype=float).dot(bv)
        if (tlo+g > zhi+eps).any() or (thi+g < zlo-eps).any():   # this image misses the zone
            continue
        s = (tris+g).dot(nvec.T)-dist                  # signed distances [Ntri,3,Nface]
        keep = ~(s > eps).all(axis=1).any(axis=1)      # drop what a single face cuts off entirely
        if not keep.any():
            continue
        full = (s <= eps).all(axis=(1, 2))
        out.extend(tris[full]+g)
        for tri in tris[keep & ~full]+g:               # only the triangles on the zone faces
            poly = tri
            for n, d in zip(nvec, dist):
                poly = _clip_half_space(poly, n, d, eps)
                if poly is None:
                    break
            if poly is not None:
                out.append(poly)
    return out
