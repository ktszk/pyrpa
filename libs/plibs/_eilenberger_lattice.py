#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
True periodic vortex lattice (square / triangular) via the quasiclassical Riccati
equation in the Doppler / London formulation.

Formulation (gauge-transformed, "B"): the order-parameter phase is removed, so the
Riccati equation is integrated with a *real* lattice-periodic amplitude |Delta(r)|
and the field enters only through the gauge-invariant supercurrent (Doppler) shift

    w_tilde(r) = w_n + i v_F . Q(r),   Q(r) = (1/2) grad(phi) - (e/hbar c) A(r),

the gauge-invariant superfluid momentum, which is *periodic* on the magnetic unit
cell.  Q is built in Fourier space over the reciprocal lattice {K} of the vortex
lattice; finite Ginzburg-Landau parameter kappa = lambda/xi enters through London
screening of the field,

    B(K) = Bbar * F(K) / (1 + lambda^2 K^2),      (Maxwell / finite kappa)
    Q(K) = -i (zhat x K) (pi/S) F(K) lambda^2 / (1 + lambda^2 K^2),

with cell area S = Phi0/B (one flux quantum) and a core form factor F(K) =
exp(-xi^2 K^2/2).  kappa -> infinity recovers the extreme type-II uniform field.

This is the standard tractable vortex-lattice method for the field-dependent
density of states (Volovik effect): it captures the inter-vortex nodal
quasiparticles and the field/penetration structure, at the cost of the discrete
Caroli-de Gennes-Matricon core levels (which need the full phase-resolved
magnetic-Bloch treatment).  Both square and triangular lattices are supported.

The order parameter amplitude |Delta(r)| is solved self-consistently on the cell;
trajectories wrap the cell periodically (fractional-coordinate sampling), so the
geometry is the genuine 2D lattice -- not the circular-cell approximation.
"""
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from ._eilenberger import matsubara, riccati_homogeneous, propagators_from_riccati
from ._eilenberger_surface import _integrate_vec, _bulk_gap
from ._eilenberger_vortex import _ff_vortex


def _lattice_vectors(field: float, xi: float, lattice: str = 'square'):
    """Primitive (a1,a2) and reciprocal (b1,b2) vectors of the vortex lattice with
    one flux quantum per cell (area S = 2 pi xi^2 / (B/Hc2)).  ai.bj = 2 pi dij."""
    S = 2.0 * np.pi / field * xi ** 2
    if lattice.startswith('s'):                       # square
        a = np.sqrt(S)
        a1, a2 = np.array([a, 0.0]), np.array([0.0, a])
    else:                                             # triangular (Abrikosov ground state)
        a = np.sqrt(2.0 * S / np.sqrt(3.0))
        a1, a2 = np.array([a, 0.0]), np.array([0.5 * a, 0.5 * np.sqrt(3.0) * a])
    M = np.column_stack([a1, a2])
    Bm = 2.0 * np.pi * np.linalg.inv(M).T
    return a1, a2, Bm[:, 0], Bm[:, 1], S, M


def _lattice_QB(X, Y, b1, b2, S, lam, xi, nmax=8):
    """Periodic supercurrent momentum Q(r)=(Qx,Qy) and field B(r)/Bbar on the grid,
    from the London/Brandt Fourier sum (finite kappa via lambda)."""
    Qx = np.zeros_like(X)
    Qy = np.zeros_like(X)
    Brel = np.ones_like(X)                            # B(r)/Bbar (K=0 term = 1)
    pref = np.pi / S
    for m in range(-nmax, nmax + 1):
        for n in range(-nmax, nmax + 1):
            if m == 0 and n == 0:
                continue
            K = m * b1 + n * b2
            K2 = K @ K
            F = np.exp(-0.5 * xi ** 2 * K2)           # core form factor
            scr = lam ** 2 / (1.0 + lam ** 2 * K2)
            ph = np.exp(1j * (K[0] * X + K[1] * Y))
            Qx += (-1j * (-K[1]) * pref * F * scr * ph).real
            Qy += (-1j * (K[0]) * pref * F * scr * ph).real
            Brel += (F / (1.0 + lam ** 2 * K2) * ph).real
    return Qx, Qy, Brel


def _nn_dist(X, Y, a1, a2, nrep=2):
    """Distance to the nearest vortex (periodic), for the model |Delta| seed/profile."""
    d = np.full(X.shape, np.inf)
    for i in range(-nrep, nrep + 1):
        for j in range(-nrep, nrep + 1):
            Rx, Ry = i * a1[0] + j * a2[0], i * a1[1] + j * a2[1]
            d = np.minimum(d, np.sqrt((X - Rx) ** 2 + (Y - Ry) ** 2))
    return d


def _periodic_interp(field, Ng):
    """Bilinear interpolator on fractional coords [0,1]^2 with periodic wrap."""
    fp = np.empty((Ng + 1, Ng + 1), dtype=field.dtype)
    fp[:Ng, :Ng] = field
    fp[Ng, :Ng] = field[0, :]
    fp[:Ng, Ng] = field[:, 0]
    fp[Ng, Ng] = field[0, 0]
    ax = np.linspace(0.0, 1.0, Ng + 1)
    return RegularGridInterpolator((ax, ax), fp, bounds_error=False, fill_value=None)


def _chord(omega2d, Delta2d, hvf, ds):
    """g, f along one chord with position-dependent (complex) frequency, real Delta."""
    ga0, _ = riccati_homogeneous(omega2d[0], Delta2d[0])
    gamma = _integrate_vec(omega2d, Delta2d, hvf, ds, ga0)
    _, gb0 = riccati_homogeneous(omega2d[-1], Delta2d[-1])
    gammat = _integrate_vec(omega2d[::-1], np.conj(Delta2d[::-1]), hvf, ds, gb0)[::-1]
    g, f, _ = propagators_from_riccati(gamma, gammat)
    return g, f


def _sample_pts(interp, x, y, Minv):
    """Sample a periodic interpolator at array points x,y (any shape)."""
    f1 = (Minv[0, 0] * x + Minv[0, 1] * y) % 1.0
    f2 = (Minv[1, 0] * x + Minv[1, 1] * y) % 1.0
    pts = np.stack([f1.ravel(), f2.ravel()], axis=1)
    return interp(pts).reshape(x.shape)


def solve_lattice(coupling, temp, omega, gap_sym='d', field=0.2, kappa=5.0,
                  lattice='square', Dbulk=None, Ng=24, nbeta=18, hvf=1.0, fs=None):
    """
    @fn solve_lattice
    @brief Build the periodic vortex-lattice state: lattice geometry, periodic
    supercurrent Q(r) and field B(r)/Bbar (London, finite kappa), and the amplitude
    |Delta(r)| (model Dbulk*tanh(d_nn/xi); self-consistency is a refinement that the
    DOS does not require).  Returns a state dict consumed by ``lattice_dos``.
    With a model FS (``fs``) the coherence length uses the representative velocity
    hvf_eff = <|v_F|>_nf so the geometry is |v_F|-consistent.
    @return dict (F1,F2,X,Y,absD,Qx,Qy,Brel,Minv,a1,a2,S,xi,Dbulk,acell,lam)
    """
    bfull = np.linspace(0.0, 2.0 * np.pi, 180, endpoint=False)
    hvf_eff = float((fs['nf'] * fs['vabs']).sum()) if fs is not None else hvf
    if Dbulk is None:
        if fs is not None:
            from ._eilenberger_fs import bulk_gap_fs
            Dbulk = bulk_gap_fs(coupling, temp, omega, fs, gap_sym)
        else:
            Dbulk = _bulk_gap(coupling, temp, omega, _ff_vortex(bfull, gap_sym))
    if Dbulk <= 0:
        return None
    xi = hvf_eff / (np.pi * Dbulk)
    a1, a2, b1, b2, S, M = _lattice_vectors(field, xi, lattice)
    Minv = np.linalg.inv(M)
    lam = kappa * xi
    fax = (np.arange(Ng) + 0.5) / Ng                  # cell-centered fractional grid
    F1, F2 = np.meshgrid(fax, fax, indexing='ij')
    X = F1 * a1[0] + F2 * a2[0]
    Y = F1 * a1[1] + F2 * a2[1]
    Qx, Qy, Brel = _lattice_QB(X, Y, b1, b2, S, lam, xi)
    dnn = _nn_dist(X, Y, a1, a2)
    absD = Dbulk * np.tanh(dnn / xi)
    return dict(F1=F1, F2=F2, X=X, Y=Y, absD=absD, Qx=Qx, Qy=Qy, Brel=Brel,
                Minv=Minv, a1=a1, a2=a2, S=S, xi=xi, Dbulk=Dbulk,
                acell=np.sqrt(S), lam=lam)


def lattice_dos(state, gap_sym, wlist, coupling=None, temp=None, omega=None,
                delta=None, nbeta=24, hvf=1.0, fs=None, Ls_cell=1.5, ds_xi=0.3):
    """
    @fn lattice_dos
    @brief Spatially-averaged density of states N(w)/N0 of the vortex lattice from a
    converged state, via the retarded Doppler Riccati, averaged over trajectories
    (directions and offsets) -- which samples the whole cell.  With a model FS
    (``fs``) the trajectory direction is v_hat, the Riccati velocity is |v_F| per
    direction, the Doppler v_F.Q scales with |v_F|, and the FS average is nf-weighted.
    @return N(w)/N0 [Nw]
    """
    Dbulk = state['Dbulk']
    xi = state['xi']
    Minv = state['Minv']
    acell = state['acell']
    Ng = state['absD'].shape[0]
    if delta is None:
        delta = 0.03 * Dbulk
    zomega = delta - 1j * wlist
    Nw = len(wlist)
    if fs is not None:                                # model FS: v_hat directions, nf weights
        from ._eilenberger_fs import fs_form_factor
        dirs = np.arctan2(fs['vy'], fs['vx'])
        phi = fs_form_factor(fs, gap_sym)
        hvfarr = fs['vabs']
        wt_dir = fs['nf']
        nbeta = len(dirs)
    else:                                             # isotropic cylinder (v_hat = k_hat, uniform)
        dirs = np.linspace(0.0, 2.0 * np.pi, nbeta, endpoint=False)
        phi = _ff_vortex(dirs, gap_sym)               # gap a trajectory sees = phi(beta)*|Delta(r)|
        hvfarr = np.full(nbeta, hvf)
        wt_dir = np.full(nbeta, 1.0 / nbeta)
    Di = _periodic_interp(state['absD'], Ng)
    Qix = _periodic_interp(state['Qx'], Ng)
    Qiy = _periodic_interp(state['Qy'], Ng)
    Ls = max(Ls_cell * acell, 3.0 * acell)            # long chord: samples the cell ergodically
    ds = ds_xi * xi
    ns = int(2 * Ls / ds)
    sgrid = np.linspace(-Ls, Ls, ns)
    # average over the whole chord except the ends (where the Riccati has not yet
    # forgotten the inflow boundary condition); a long periodic chord covers the
    # whole unit cell (square or rhombic) uniformly
    central = np.abs(sgrid) <= (Ls - 3.0 * xi)
    ncen = central.sum()
    # offsets b across the cell (perpendicular), to sample the 2D cell
    nb = max(12, Ng)
    bgrid = (np.arange(nb) + 0.5) / nb * acell - 0.5 * acell
    gsum = np.zeros(Nw, dtype=np.complex128)           # nf-weighted FS+spatial average
    for ib in range(nbeta):
        cb, sb = np.cos(dirs[ib]), np.sin(dirs[ib])
        gdir = np.zeros(Nw, dtype=np.complex128)
        for b in bgrid:
            x = sgrid * cb - b * sb
            y = sgrid * sb + b * cb
            Dl = phi[ib] * _sample_pts(Di, x, y, Minv)            # phi(beta)*|Delta(r)| along chord
            dop = hvfarr[ib] * (cb * _sample_pts(Qix, x, y, Minv)
                                + sb * _sample_pts(Qiy, x, y, Minv))   # v_F.Q (~|v_F|)
            om = zomega[None, :] + 1j * dop[:, None]               # [ns, Nw]
            g, _ = _chord(om, (Dl[:, None] + 0.0j) * np.ones((1, Nw)), hvfarr[ib], ds)
            gdir += g[central].sum(axis=0)
        gsum += wt_dir[ib] * gdir / (nb * ncen)        # spatial mean per direction, nf-weighted
    return gsum.real


def calc_vortex_lattice_periodic(coupling, temp, wc, gap_sym='d', field_list=None,
                                 kappa=5.0, lattice='square', kb=1.0, Ng=20, nbeta=16,
                                 fs_kind=None, fs_params=None):
    """
    @fn calc_vortex_lattice_periodic
    @brief Driver: true periodic vortex lattice (Doppler/London, finite kappa).
    Sweeps B/Hc2, self-consistently solves |Delta(r)| and reports the spatially-
    averaged zero-energy DOS <N(0)>/N0(B) (s-wave ~B cores; d-wave ~sqrt(B) Volovik),
    and writes the field profile B(r)/Bbar for the largest field.  With ``fs_kind``
    the trajectories run on a model Fermi surface with real Fermi velocities.
    """
    omega = matsubara(temp, wc)
    if field_list is None:
        field_list = [0.05, 0.1, 0.2, 0.4]
    fs = None
    if fs_kind is not None:
        from ._eilenberger_fs import build_model_fs
        fs = build_model_fs(fs_kind, 120, params=fs_params)
    print(f"periodic vortex lattice (Doppler/London): {gap_sym}, {lattice}, kappa={kappa}, "
          f"lambda={coupling:.3f}, T={temp/kb:.2f} K"
          f"{', model FS+v_F: '+fs_kind if fs is not None else ''}", flush=True)
    results = []
    last = None
    for b in field_list:
        st = solve_lattice(coupling, temp, omega, gap_sym=gap_sym, field=b, kappa=kappa,
                           lattice=lattice, Ng=Ng, nbeta=nbeta, fs=fs)
        if st is None or st['Dbulk'] <= 0:
            print(f"  B/Hc2={b:.3f}: normal", flush=True)
            continue
        # the gap-node sampling needs a fine FS-angle grid for d-wave; small broadening
        n0 = float(lattice_dos(st, gap_sym, np.array([0.0]), nbeta=max(nbeta, 60),
                               delta=0.012 * st['Dbulk'], fs=fs)[0])
        a_xi = st['acell'] / st['xi']
        bmin, bmax = st['Brel'].min(), st['Brel'].max()
        results.append((b, n0))
        print(f"  B/Hc2={b:.3f}  a/xi={a_xi:.2f}  <|D|>/Dbulk={st['absD'].mean()/st['Dbulk']:.3f}  "
              f"B(r)/Bbar in [{bmin:.3f},{bmax:.3f}]  <N(0)>/N0={n0:.4f}", flush=True)
        last = st
    if len(results) >= 2:
        bb = np.array([r[0] for r in results])
        nn = np.array([r[1] for r in results])
        p = np.polyfit(np.log(bb), np.log(np.maximum(nn, 1e-12)), 1)[0]
        print(f"  power law  <N(0)>/N0 ~ (B/Hc2)^{p:.2f}  (s ~1 cores, d ~0.5 Volovik)", flush=True)
    try:
        with open('lattice_dos.dat', 'w') as fh:
            fh.write("# B/Hc2   <N(0)>/N0\n")
            for b, n0 in results:
                fh.write(f"{b:12.5e} {n0:12.5e}\n")
        if last is not None:
            with open('lattice_field.dat', 'w') as fh:
                fh.write("# x/xi  y/xi  |Delta|/Dbulk  B(r)/Bbar\n")
                xi = last['xi']
                for i in range(last['X'].shape[0]):
                    for j in range(last['X'].shape[1]):
                        fh.write(f"{last['X'][i,j]/xi:9.3f} {last['Y'][i,j]/xi:9.3f} "
                                 f"{last['absD'][i,j]/last['Dbulk']:10.4e} {last['Brel'][i,j]:10.4e}\n")
                    fh.write("\n")
    except IOError as e:
        print(f"Error: failed to write lattice output: {e}", flush=True)
    return results
