#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Isolated-vortex quasiclassical Eilenberger solver via the Riccati equation
(first step toward the vortex lattice).

This adds the two ingredients that distinguish the vortex problem from the flat
surface: a genuine 2D real-space order parameter and the phase winding of a
vortex, Delta(r) = D(rho) e^{i theta}.  The quasiclassical propagators are
obtained by integrating the same stable Riccati step used for the surface, now
along straight chords through the 2D plane; the order parameter sampled along a
chord is complex and winds in phase, which is what binds the Caroli-de Gennes-
Matricon (CdGM) zero-energy core state.

Symmetry reduction (isotropic FS, axisymmetric s-wave vortex): the Fermi-surface
average over v_F directions equals an azimuthal average, because rotating the
trajectory direction is equivalent to rotating the observation point.  So it is
enough to integrate chords for a single v_F direction (taken along x, one chord
per impact parameter b = grid row), obtain the field on the (x, y) plane, and
average over the azimuth theta on circles of radius rho.

Scope of this first version: a single vortex, clean limit, s-wave, extreme
type-II (the vector potential / supercurrent Doppler shift is neglected; it
disperses the core state away from rho=0 but the rho=0 zero-energy peak remains).
The vortex *lattice* then adds: a periodic magnetic unit cell, the vector-
potential phase along trajectories, and quasi-periodic (magnetic-Bloch) boundary
conditions -- on top of the same chord integration and self-consistency here.
"""
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from ._eilenberger import riccati_homogeneous, propagators_from_riccati, matsubara
from ._eilenberger_surface import _integrate_vec, _bulk_gap, _bulk_gap_imp


def _chord_gf(omega: np.ndarray, Dchord: np.ndarray, hvf: float, ds: float):
    """Quasiclassical g, f along one straight chord (vectorized over frequency).
    gamma is integrated forward from the upstream bulk root, gamma-tilde backward
    from the downstream bulk root, both with the unconditionally stable step.
    @param  omega: frequencies [Nw] (Matsubara or retarded)
    @param Dchord: complex order parameter sampled along the chord [Ns]
    @param    hvf: hbar |v_F|
    @param     ds: step along the chord
    @return (g, f): propagators along the chord [Ns, Nw]
    """
    ga0, _ = riccati_homogeneous(omega, Dchord[0])
    gamma = _integrate_vec(omega, Dchord, hvf, ds, ga0)
    _, gb0 = riccati_homogeneous(omega, Dchord[-1])
    gammat = _integrate_vec(omega, np.conj(Dchord[::-1]), hvf, ds, gb0)[::-1]
    g, f, _ = propagators_from_riccati(gamma, gammat)
    return g, f


def _chord_gf_pos(omega2d: np.ndarray, Delta2d: np.ndarray, hvf: float, ds: float):
    """Like _chord_gf but with position-dependent renormalized fields along the
    chord: omega2d and Delta2d are [Ns, Nw] (needed once impurities make w_tilde(r)
    and Sigma_f(r) vary in space)."""
    ga0, _ = riccati_homogeneous(omega2d[0], Delta2d[0])
    gamma = _integrate_vec(omega2d, Delta2d, hvf, ds, ga0)
    _, gb0 = riccati_homogeneous(omega2d[-1], Delta2d[-1])
    gammat = _integrate_vec(omega2d[::-1], np.conj(Delta2d[::-1]), hvf, ds, gb0)[::-1]
    g, f, _ = propagators_from_riccati(gamma, gammat)
    return g, f


def _eval_field(field: np.ndarray, xg: np.ndarray, px: np.ndarray, py: np.ndarray,
                fill=0.0) -> np.ndarray:
    """Evaluate a complex grid field at scattered points (px, py).  ``field`` is
    [ng, ng] or [ng, ng, Nw]; returns px.shape (+ (Nw,)).  Real and imaginary
    parts are interpolated separately (bilinear)."""
    pts = np.stack([px.ravel(), py.ravel()], axis=1)
    re = RegularGridInterpolator((xg, xg), field.real, bounds_error=False, fill_value=fill)(pts)
    im = RegularGridInterpolator((xg, xg), field.imag, bounds_error=False, fill_value=fill)(pts)
    out = re + 1j * im
    shape = px.shape + ((field.shape[-1],) if field.ndim == 3 else ())
    return out.reshape(shape)


def _delta_grid(Dr: np.ndarray, rgrid: np.ndarray, Rg: np.ndarray, phase: np.ndarray,
                Dbulk: float) -> np.ndarray:
    """Build the complex order parameter on the 2D grid from the radial amplitude
    D(rho): Delta = D(rho) * e^{i theta}."""
    Damp = np.interp(Rg.ravel(), rgrid, Dr, right=Dbulk).reshape(Rg.shape)
    return Damp * phase


def _grid(Lxi: float, xi: float, ngrid: int):
    """Square (x, y) grid of half width Lxi*xi, the phase factor e^{i theta} and
    the radius on the grid, plus a radial grid for the amplitude profile."""
    R = Lxi * xi
    xg = np.linspace(-R, R, ngrid)
    X, Y = np.meshgrid(xg, xg, indexing='ij')        # chord = column j (fixed y), varies x
    Rg = np.sqrt(X ** 2 + Y ** 2)
    phase = np.where(Rg > 0, (X + 1j * Y) / np.where(Rg > 0, Rg, 1.0), 0.0)
    rgrid = np.linspace(0.0, R, ngrid // 2)
    return xg, Rg, phase, rgrid, R


def _azimuthal(field: np.ndarray, xg: np.ndarray, rgrid: np.ndarray,
               theta: np.ndarray, weight_phase: bool):
    """Azimuthal average of a complex (x, y, Nw) field on circles of radius rgrid.
    If weight_phase, multiply by e^{-i theta} first (gap equation: projects out the
    vortex winding); otherwise plain average (LDOS).
    @return averaged array [Nr, Nw]
    """
    fre = RegularGridInterpolator((xg, xg), field.real, bounds_error=False, fill_value=0.0)
    fim = RegularGridInterpolator((xg, xg), field.imag, bounds_error=False, fill_value=0.0)
    out = np.empty((len(rgrid), field.shape[-1]), dtype=np.complex128)
    emit = np.exp(-1j * theta)[:, None] if weight_phase else 1.0
    for k, rho in enumerate(rgrid):
        pts = np.stack([rho * np.cos(theta), rho * np.sin(theta)], axis=1)
        fc = fre(pts) + 1j * fim(pts)                 # [ntheta, Nw]
        out[k] = (fc * emit).mean(axis=0)
    return out


def solve_vortex(coupling: float, temp: float, omega: np.ndarray, Dbulk: float = None,
                 Lxi: float = 8.0, ngrid: int = 81, ntheta: int = 48, hvf: float = 1.0,
                 eps: float = 1.0e-3, itemax: int = 100, mix: float = 0.3):
    """
    @fn solve_vortex
    @brief Self-consistently solve the radial amplitude D(rho) of an isolated
    s-wave vortex (Matsubara, clean, extreme type-II).

    @param coupling: separable pairing coupling lambda
    @param     temp: temperature [eV]
    @param    omega: positive Matsubara frequencies [Nw]
    @param    Dbulk: bulk gap [eV]; computed from the homogeneous gap equation if None
    @param      Lxi: cell half width in coherence lengths xi = hvf/(pi*Dbulk)
    @param    ngrid: number of grid points per axis
    @param   ntheta: azimuthal samples for the FS / gap average
    @param      hvf: hbar |v_F| (length unit)
    @return (rgrid, Dr, Dbulk, xi): radial grid, amplitude profile, bulk gap, xi
    """
    if Dbulk is None:
        Dbulk = _bulk_gap(coupling, temp, omega, np.ones(16))   # s-wave, <phi^2>=1
    if Dbulk < 1.0e-6 * temp:
        xi = hvf / (np.pi * max(temp, 1e-12))
        rgrid = np.linspace(0.0, Lxi * xi, ngrid // 2)
        return rgrid, np.zeros_like(rgrid), 0.0, xi

    xi = hvf / (np.pi * Dbulk)
    xg, Rg, phase, rgrid, R = _grid(Lxi, xi, ngrid)
    dx = xg[1] - xg[0]
    theta = np.linspace(0.0, 2.0 * np.pi, ntheta, endpoint=False)
    Nw = len(omega)
    Dr = Dbulk * np.tanh(rgrid / xi)               # standard vortex ansatz
    for it in range(itemax):
        Delta = _delta_grid(Dr, rgrid, Rg, phase, Dbulk)
        f_grid = np.empty((ngrid, ngrid, Nw), dtype=np.complex128)
        for j in range(ngrid):                     # one chord per impact parameter (column)
            _, f = _chord_gf(omega, Delta[:, j], hvf, dx)
            f_grid[:, j] = f
        # gap equation: D(rho) = lambda 2T sum_n <f e^{-i theta}>_azimuth
        favg = _azimuthal(f_grid, xg, rgrid, theta, weight_phase=True)   # [Nr, Nw]
        Dr_new = coupling * 2.0 * temp * favg.sum(axis=1).real
        Dr_new[0] = 0.0                            # amplitude vanishes at the core
        err = np.abs(Dr_new - Dr).max() / Dbulk
        Dr = (1.0 - mix) * Dr + mix * Dr_new
        if err < eps:
            break
    return rgrid, Dr, Dbulk, xi


def vortex_ldos(Dr: np.ndarray, rgrid: np.ndarray, xi: float, wlist: np.ndarray,
                Dbulk: float, delta: float = None, Lxi: float = 8.0, ngrid: int = 81,
                ntheta: int = 48, hvf: float = 1.0, rho_list: np.ndarray = None):
    """
    @fn vortex_ldos
    @brief Local density of states N(rho, w)/N0 from the converged vortex profile,
    via the retarded Riccati equation (iw_n -> w + i delta).  N(0, w) shows the
    zero-energy Caroli-de Gennes-Matricon core bound state (the STM zero-bias peak).

    @param Dr,rgrid,xi: converged radial profile from solve_vortex
    @param   wlist: real-frequency grid [Nw]
    @param   Dbulk: bulk gap (for the default broadening)
    @param   delta: retarded broadening [eV] (default 0.02*Dbulk)
    @param rho_list: radii at which to return N (default = rgrid)
    @return (rho_out, ldos): radii [Nr2] and N(rho, w)/N0 [Nr2, Nw]
    """
    if delta is None:
        delta = 0.02 * Dbulk
    xg, Rg, phase, _, R = _grid(Lxi, xi, ngrid)
    dx = xg[1] - xg[0]
    theta = np.linspace(0.0, 2.0 * np.pi, ntheta, endpoint=False)
    zomega = delta - 1j * wlist
    Delta = _delta_grid(Dr, rgrid, Rg, phase, Dbulk)
    g_grid = np.empty((ngrid, ngrid, len(wlist)), dtype=np.complex128)
    for j in range(ngrid):
        g, _ = _chord_gf(zomega, Delta[:, j], hvf, dx)
        g_grid[:, j] = g
    rho_out = rgrid if rho_list is None else np.asarray(rho_list)
    gavg = _azimuthal(g_grid, xg, rho_out, theta, weight_phase=False)    # [Nr, Nw]
    return rho_out, gavg.real


def _ff_vortex(beta: np.ndarray, gap_sym: str) -> np.ndarray:
    """Pairing form factor phi_d(k_hat) vs the v_F direction beta, normalized so
    <phi^2> = 1 over the full circle.  For a d-wave vortex the gap a trajectory
    sees is phi(beta) * Psi(r), so the trajectory direction enters explicitly and
    the rotational symmetry of the s-wave vortex is reduced to fourfold."""
    if gap_sym == 's':
        return np.ones_like(beta)
    if gap_sym == 'd':
        return np.sqrt(2.0) * np.cos(2.0 * beta)
    if gap_sym == 'dxy':
        return np.sqrt(2.0) * np.sin(2.0 * beta)
    if gap_sym == 'px':
        return np.sqrt(2.0) * np.cos(beta)
    if gap_sym == 'py':
        return np.sqrt(2.0) * np.sin(beta)
    if gap_sym in ('p+ip', 'chiral'):
        return np.exp(1j * beta) * np.ones_like(beta, dtype=np.complex128)
    if gap_sym == 'p-ip':
        return np.exp(-1j * beta) * np.ones_like(beta, dtype=np.complex128)
    raise ValueError(f"unknown gap_sym: {gap_sym}")


def _doppler_chord(Lx, Ly, Rc, cb, sb, rho_min):
    """Field (vector-potential) Doppler shift along a chord for the circular cell.

    Here the order parameter keeps its vortex phase e^{i theta} (formulation A), so
    the Doppler term must be the *vector-potential* contribution ONLY -- the phase
    gradient grad(phi)/2 is already carried by Delta's phase and must NOT be added
    again (doing so double-counts the vortex supercurrent).  For the symmetric-gauge
    uniform field of a one-flux-quantum cell, (e/hbar c) A_theta = rho/(2 Rc^2), so

        delta(r) = - v_F . (e/hbar c) A = -(rho/(2 Rc^2)) (v_hat . theta_hat).

    (The full supervelocity p_s = grad(phi)/2 - (e/hbar c)A = (1/2rho)(1-rho^2/Rc^2)
    would be the right object only in the gauge-transformed real-|Delta| formulation
    B used by _eilenberger_lattice; mixing it with the kept phase is inconsistent.)
    """
    rho2 = Lx ** 2 + Ly ** 2
    rho2s = np.maximum(rho2, rho_min ** 2)
    vdoth = (-Ly * cb + Lx * sb)            # rho * (v_hat . theta_hat)
    delta = -(rho2 / Rc ** 2) * vdoth / (2.0 * rho2s)   # vector-potential part only
    return np.where(rho2 < Rc ** 2, delta, 0.0)


def solve_vortex2d(coupling: float, temp: float, omega: np.ndarray, gap_sym: str = 'd',
                   Dbulk: float = None, Lxi: float = 8.0, ngrid: int = 49, nbeta: int = 24,
                   hvf: float = 1.0, imp_gamma: float = 0.0, imp_c: float = 1.0e8,
                   field: float = 0.0, fs: dict = None, eps: float = 2.0e-3,
                   itemax: int = 80, mix: float = 0.3):
    """
    @fn solve_vortex2d
    @brief Self-consistently solve the complex 2D order-parameter field Psi(r) of a
    vortex for an arbitrary pairing symmetry.  Delta(r, k_hat) = phi(beta) Psi(r).

    For each Fermi-velocity direction beta the order parameter is sampled on a grid
    rotated to align trajectories with one axis, the Riccati equation is integrated
    along the chords, and the result is interpolated back to the (x, y) grid and
    accumulated with conj(phi(beta)).  This is the full vortex-lattice trajectory
    machinery (no symmetry reduction).

    Field (vortex lattice, circular-cell / Wigner-Seitz approximation): field =
    B/Hc2 > 0 sets a magnetic unit cell of radius Rc = sqrt(2/field)*xi (one flux
    quantum).  The order parameter keeps its vortex phase e^{i theta}, so the field
    enters only through the vector-potential Doppler shift w_tilde(r)=w_n - i v_F.(e/hbar c)A
    (see _doppler_chord; the phase gradient is already in Delta -- adding the full
    supervelocity would double-count it).  Extreme type-II, no Maxwell back-reaction.
    For the field-dependent DOS prefer the gauge-transformed periodic-lattice method
    in _eilenberger_lattice (square/triangular, finite kappa).

    With impurities (imp_gamma>0) the local non-magnetic T-matrix self-energy
    w_tilde(r,w_n), Sigma_f(r,w_n) is also relaxed together with Psi(r).

    @param  gap_sym: 's','d','dxy','px','py','p+ip'/'p-ip'
    @param      Lxi: cell half width in xi (isolated vortex, field=0)
    @param    ngrid,nbeta: grid points per axis / number of v_F directions
    @param imp_gamma,imp_c: normal-state scattering rate [eV] and T-matrix c (0=clean)
    @param    field: B/Hc2 (0 = isolated vortex; >0 = circular-cell vortex lattice)
    @return (xg, Psi, Dbulk, xi): grid axis, complex field [ngrid,ngrid], bulk gap, xi
    """
    no_imp = (imp_gamma == 0.0)
    use_fs = fs is not None
    if use_fs and (not no_imp or field > 0.0):
        raise NotImplementedError("model-FS vortex currently supports the clean, zero-field case")
    bfull = np.linspace(0.0, 2.0 * np.pi, 180, endpoint=False)
    if Dbulk is None:
        if use_fs:
            from ._eilenberger_fs import bulk_gap_fs
            Dbulk = bulk_gap_fs(coupling, temp, omega, fs, gap_sym)
        else:
            Dbulk = (_bulk_gap(coupling, temp, omega, _ff_vortex(bfull, gap_sym)) if no_imp
                     else _bulk_gap_imp(coupling, temp, omega, _ff_vortex(bfull, gap_sym), imp_gamma, imp_c))
    if Dbulk < 1.0e-6 * temp:
        xi = hvf / (np.pi * max(temp, 1e-12))
        xg = np.linspace(-Lxi * xi, Lxi * xi, ngrid)
        return xg, np.zeros((ngrid, ngrid), dtype=np.complex128), 0.0, xi

    xi = hvf / (np.pi * Dbulk)
    Rc = np.sqrt(2.0 / field) * xi if field > 0.0 else np.inf   # WS cell radius (one flux quantum)
    R = Rc if field > 0.0 else Lxi * xi
    use_pos = (not no_imp) or (field > 0.0)                     # position-dependent w_tilde
    xg = np.linspace(-R, R, ngrid)
    dx = xg[1] - xg[0]
    rho_min = 0.5 * dx
    X, Y = np.meshgrid(xg, xg, indexing='ij')
    Rg = np.sqrt(X ** 2 + Y ** 2)
    theta = np.arctan2(Y, X)                            # analytic vortex phase
    A = Dbulk * np.tanh(Rg / xi)                        # real amplitude field |Psi|
    if use_fs:                                          # model FS: trajectory along v_hat, nf weights
        from ._eilenberger_fs import fs_form_factor
        dirs = np.arctan2(fs['vy'], fs['vx'])           # v_F direction angles
        phi = fs_form_factor(fs, gap_sym)               # phi at the k-direction
        hvfarr = fs['vabs']                             # |v_F| per trajectory (Riccati velocity)
        wt_dir = fs['nf']                               # FS-average weight (sum 1)
        nbeta = len(dirs)
    else:                                               # isotropic cylinder (v_hat = k_hat, uniform)
        dirs = np.linspace(0.0, 2.0 * np.pi, nbeta, endpoint=False)
        phi = _ff_vortex(dirs, gap_sym)
        hvfarr = np.full(nbeta, hvf)
        wt_dir = np.full(nbeta, 1.0 / nbeta)
    SS, BB = np.meshgrid(xg, xg, indexing='ij')        # rotated-frame sample grid (s, b)
    Nw = len(omega)
    gpref = imp_gamma * (imp_c ** 2 + 1.0)
    wt = np.broadcast_to(omega, (ngrid, ngrid, Nw)).copy()   # impurity-renormalized w (no Doppler)
    sigf = np.zeros((ngrid, ngrid, Nw), dtype=np.complex128)
    for it in range(itemax):
        Ai = RegularGridInterpolator((xg, xg), A, bounds_error=False, fill_value=Dbulk)
        accf = np.zeros((ngrid, ngrid), dtype=np.complex128)        # sum_beta conj(phi) sum_n f
        gacc = np.zeros((ngrid, ngrid, Nw), dtype=np.complex128)    # sum_beta g (for <g>)
        facc = np.zeros((ngrid, ngrid, Nw), dtype=np.complex128)
        dwt = wt - omega                                # impurity deviation (->0 in bulk)
        for ib in range(nbeta):
            cb, sb = np.cos(dirs[ib]), np.sin(dirs[ib])
            hvf_i = hvfarr[ib]
            Lx = SS * cb - BB * sb
            Ly = SS * sb + BB * cb
            base = phi[ib] * Ai((Lx, Ly)) * np.exp(1j * np.arctan2(Ly, Lx))   # phi*Psi [ns,nb]
            sxy = X * cb + Y * sb
            bxy = -X * sb + Y * cb
            if not use_pos:
                sumf = np.empty((ngrid, ngrid), dtype=np.complex128)
                for j in range(ngrid):
                    _, f = _chord_gf(omega, base[:, j], hvf_i, dx)
                    sumf[:, j] = f.sum(axis=1)
                accf += wt_dir[ib] * np.conj(phi[ib]) * _eval_field(sumf, xg, sxy, bxy, fill=0.0)
            else:
                om_rot = np.broadcast_to(omega, (ngrid, ngrid, Nw)).astype(np.complex128).copy()
                if field > 0.0:                          # Doppler shift along the chord
                    om_rot = om_rot + 1j * _doppler_chord(Lx, Ly, Rc, cb, sb, rho_min)[:, :, None]
                Dtraj = np.broadcast_to(base[:, :, None], (ngrid, ngrid, Nw)).copy()
                if not no_imp:
                    om_rot = om_rot + _eval_field(dwt, xg, Lx, Ly, fill=0.0)
                    Dtraj = Dtraj + _eval_field(sigf, xg, Lx, Ly, fill=0.0)
                g_rot = np.empty((ngrid, ngrid, Nw), dtype=np.complex128)
                f_rot = np.empty((ngrid, ngrid, Nw), dtype=np.complex128)
                for j in range(ngrid):
                    g_rot[:, j], f_rot[:, j] = _chord_gf_pos(om_rot[:, j], Dtraj[:, j], hvf_i, dx)
                g_xy = _eval_field(g_rot, xg, sxy, bxy, fill=0.0)
                f_xy = _eval_field(f_rot, xg, sxy, bxy, fill=0.0)
                gacc += wt_dir[ib] * g_xy
                facc += wt_dir[ib] * f_xy
                accf += wt_dir[ib] * np.conj(phi[ib]) * f_xy.sum(axis=2)
        A_new = np.maximum(coupling * 2.0 * temp * (accf * np.exp(-1j * theta)).real, 0.0)
        err = np.abs(A_new - A).max() / Dbulk
        A = (1.0 - mix) * A + mix * A_new
        if not no_imp:
            avg_g = gacc
            avg_f = facc
            Dimp = imp_c ** 2 + avg_g ** 2 + avg_f * np.conj(avg_f)
            wt_new = omega + gpref * avg_g / Dimp
            sigf_new = gpref * avg_f / Dimp
            err = max(err, np.abs(wt_new - wt).max() / max(omega[0], 1e-12))
            wt = (1.0 - mix) * wt + mix * wt_new
            sigf = (1.0 - mix) * sigf + mix * sigf_new
        if err < eps:
            break
    return xg, A * np.exp(1j * theta), Dbulk, xi


def vortex_ldos2d(Psi: np.ndarray, xg: np.ndarray, xi: float, wlist: np.ndarray,
                  gap_sym: str, Dbulk: float, delta: float = None, nbeta: int = 36,
                  hvf: float = 1.0, imp_gamma: float = 0.0, imp_c: float = 1.0e8,
                  field: float = 0.0, h: float = 0.0, fs: dict = None, eps: float = 3.0e-3,
                  itemax: int = 60, mix: float = 0.5):
    """
    @fn vortex_ldos2d
    @brief 2D local density of states N(x,y,w)/N0 from a converged Psi field, via
    the retarded Riccati equation.  Reveals the fourfold d-wave core star; with
    impurities the retarded local T-matrix self-energy is solved self-consistently
    in space; with field>0 the supercurrent Doppler shift fills states between
    vortices (Volovik); with Zeeman h the spin sectors split w -> w -+ h (for a
    unitary fixed-d gap, N = (1/2) sum_sigma N(w - sigma h)), splitting the core
    bound states.
    @param field: B/Hc2 (0 = isolated vortex; >0 = circular-cell lattice Doppler shift)
    @param     h: Zeeman (Maki) energy mu*B (spin splitting of the core LDOS)
    @return ldos: N(x,y,w)/N0 [ngrid, ngrid, Nw]
    """
    if delta is None:
        delta = 0.02 * Dbulk
    ngrid = len(xg)
    dx = xg[1] - xg[0]
    rho_min = 0.5 * dx
    Rc = np.sqrt(2.0 / field) * xi if field > 0.0 else np.inf
    X, Y = np.meshgrid(xg, xg, indexing='ij')
    Nw = len(wlist)
    no_imp = (imp_gamma == 0.0)
    use_pos = (not no_imp) or (field > 0.0)
    use_fs = fs is not None
    if use_fs and use_pos:
        raise NotImplementedError("model-FS vortex LDOS currently supports the clean, zero-field case")
    if use_fs:
        from ._eilenberger_fs import fs_form_factor
        dirs = np.arctan2(fs['vy'], fs['vx'])
        phi = fs_form_factor(fs, gap_sym)
        hvfarr = fs['vabs']
        wt_dir = fs['nf']
        nbeta = len(dirs)
    else:
        dirs = np.linspace(0.0, 2.0 * np.pi, nbeta, endpoint=False)
        phi = _ff_vortex(dirs, gap_sym)
        hvfarr = np.full(nbeta, hvf)
        wt_dir = np.full(nbeta, 1.0 / nbeta)
    gpref = imp_gamma * (imp_c ** 2 + 1.0)
    SS, BB = np.meshgrid(xg, xg, indexing='ij')
    Ai = RegularGridInterpolator((xg, xg), np.abs(Psi), bounds_error=False, fill_value=Dbulk)
    base, sxyb, dopp = [], [], []
    for ib in range(nbeta):
        cb, sb = np.cos(dirs[ib]), np.sin(dirs[ib])
        Lx = SS * cb - BB * sb
        Ly = SS * sb + BB * cb
        base.append(phi[ib] * Ai((Lx, Ly)) * np.exp(1j * np.arctan2(Ly, Lx)))
        sxyb.append((X * cb + Y * sb, -X * sb + Y * cb, Lx, Ly))
        dopp.append(_doppler_chord(Lx, Ly, Rc, cb, sb, rho_min) if field > 0.0 else None)

    def _sector(zomega):
        """N(x,y,w) for a given (possibly spin-shifted) retarded frequency [Nw]."""
        wt = np.broadcast_to(zomega, (ngrid, ngrid, Nw)).copy()
        sigf = np.zeros((ngrid, ngrid, Nw), dtype=np.complex128)
        for it in range(1 if no_imp else itemax):
            gacc = np.zeros((ngrid, ngrid, Nw), dtype=np.complex128)
            facc = np.zeros((ngrid, ngrid, Nw), dtype=np.complex128)
            dwt = wt - zomega
            for ib in range(nbeta):
                sxy, bxy, Lx, Ly = sxyb[ib]
                if not use_pos:
                    g_rot = np.empty((ngrid, ngrid, Nw), dtype=np.complex128)
                    for j in range(ngrid):
                        g, _ = _chord_gf(zomega, base[ib][:, j], hvfarr[ib], dx)
                        g_rot[:, j] = g
                    gacc += wt_dir[ib] * _eval_field(g_rot, xg, sxy, bxy, fill=0.0)
                else:
                    om_rot = np.broadcast_to(zomega, (ngrid, ngrid, Nw)).astype(np.complex128).copy()
                    if field > 0.0:
                        om_rot = om_rot + 1j * dopp[ib][:, :, None]
                    Dtraj = np.broadcast_to(base[ib][:, :, None], (ngrid, ngrid, Nw)).copy()
                    if not no_imp:
                        om_rot = om_rot + _eval_field(dwt, xg, Lx, Ly, fill=0.0)
                        Dtraj = Dtraj + _eval_field(sigf, xg, Lx, Ly, fill=0.0)
                    g_rot = np.empty((ngrid, ngrid, Nw), dtype=np.complex128)
                    f_rot = np.empty((ngrid, ngrid, Nw), dtype=np.complex128)
                    for j in range(ngrid):
                        g_rot[:, j], f_rot[:, j] = _chord_gf_pos(om_rot[:, j], Dtraj[:, j], hvfarr[ib], dx)
                    gacc += wt_dir[ib] * _eval_field(g_rot, xg, sxy, bxy, fill=0.0)
                    if not no_imp:
                        facc += wt_dir[ib] * _eval_field(f_rot, xg, sxy, bxy, fill=0.0)
            if no_imp:
                return gacc.real
            avg_g = gacc
            avg_f = facc
            Dimp = imp_c ** 2 + avg_g ** 2 + avg_f * np.conj(avg_f)
            wt_new = zomega + gpref * avg_g / Dimp
            sigf_new = gpref * avg_f / Dimp
            err = np.abs(wt_new - wt).max() / max(abs(zomega[-1]), 1e-12)
            wt = (1.0 - mix) * wt + mix * wt_new
            sigf = (1.0 - mix) * sigf + mix * sigf_new
            if err < eps:
                break
        return avg_g.real

    base_z = delta - 1j * wlist
    if h == 0.0:
        return _sector(base_z)
    return 0.5 * (_sector(base_z + 1j * h) + _sector(base_z - 1j * h))


def vortex_field_profile(rgrid, Dr, Dbulk, xi, kappa, omega):
    """
    @fn vortex_field_profile
    @brief Self-consistent (finite-kappa Maxwell) magnetic field profile B(rho) of
    an isolated vortex from the quasiclassical local superfluid density.

    The local superfluid density (from the converged gap |Delta(rho)|),
        rho_s(rho)/rho_s_bulk = sum_n D(rho)^2/(w^2+D^2)^{3/2} / sum_n Dbulk^2/(...)
    sets a position-dependent penetration depth lambda(rho)=kappa*xi/sqrt(rho_s),
    which diverges at the core (where the condensate is depleted) and recovers to
    lambda_bulk=kappa*xi far away.  The radial London equation
        lambda(rho)^2 [-(1/rho) d/drho (rho dB/drho)] + B = const
    (Neumann BCs, cell-averaged flux fixed) is solved as one linear system; B(rho)
    is peaked at the core and screened over ~lambda away from it.
    @return (rho_s_ratio, lam, Brel): all on rgrid; Brel = B(rho)/<B> (mean 1)
    """
    Nr = len(rgrid)
    w2 = omega[None, :] ** 2
    sden = lambda D: (D ** 2 / (w2 + D ** 2) ** 1.5).sum(axis=1)        # sum_n per rho
    rs_bulk = float((Dbulk ** 2 / (omega ** 2 + Dbulk ** 2) ** 1.5).sum())
    rho_s = sden(Dr[:, None]) / max(rs_bulk, 1e-30)                     # rho_s(rho)/rho_s_bulk
    rho_s = np.clip(rho_s, 1e-4, None)
    lam = kappa * xi / np.sqrt(rho_s)                                   # local penetration depth
    # radial operator K[B] = -(1/rho) d/drho(rho dB/drho), Neumann BCs
    dr = rgrid[1] - rgrid[0]
    K = np.zeros((Nr, Nr))
    for i in range(Nr):
        rho = rgrid[i]
        if i == 0:
            K[0, 0] = 4.0 / dr ** 2          # -2 B'' at rho=0 (regular, B'(0)=0)
            K[0, 1] = -4.0 / dr ** 2
        elif i == Nr - 1:
            rm = 0.5 * (rgrid[i] + rgrid[i - 1])
            K[i, i] = rm / (rho * dr ** 2)    # B'(R)=0 (outer flux term drops)
            K[i, i - 1] = -rm / (rho * dr ** 2)
        else:
            rp = 0.5 * (rgrid[i] + rgrid[i + 1])
            rm = 0.5 * (rgrid[i] + rgrid[i - 1])
            K[i, i] = (rp + rm) / (rho * dr ** 2)
            K[i, i + 1] = -rp / (rho * dr ** 2)
            K[i, i - 1] = -rm / (rho * dr ** 2)
    # London eq  lambda^2 K[B] + B = S  with the vortex source S localized at the
    # core (the phase singularity, smeared over the coherence length xi); away from
    # the core B ~ K0(rho/lambda) -- peaked at the core, screened over lambda.
    src = np.exp(-0.5 * (rgrid / xi) ** 2)
    A = (lam ** 2)[:, None] * K + np.eye(Nr)
    u = np.linalg.solve(A, src)
    flux_mean = (u * rgrid).sum() / rgrid.sum()                        # cell flux average (2 pi rho dr)
    Brel = u / max(flux_mean, 1e-30)                                   # B(rho)/<B>, mean=1
    return rho_s, lam, Brel


def calc_vortex(coupling: float, temp: float, wc: float, gap_sym: str = 's', kb: float = 1.0,
                sw_ldos: bool = True, imp_gamma: float = 0.0, imp_c: float = 1.0e8,
                field: float = 0.0, h: float = 0.0, kappa: float = 0.0,
                fs_kind: str = None, fs_params=None, Lxi: float = 8.0, ngrid: int = 81):
    """
    @fn calc_vortex
    @brief Driver: self-consistent vortex profile and (optionally) LDOS.
    Clean s-wave at zero field/Zeeman uses the fast axisymmetric (radial) solver;
    any non-s pairing, finite impurity rate, finite field (field=B/Hc2>0,
    circular-cell vortex lattice), or finite Zeeman h uses the general 2D solver
    and writes the zero-energy LDOS map.

    @param  gap_sym: 's','d','dxy','px','py','p+ip'/'p-ip'
    @param  sw_ldos: also compute the real-frequency LDOS (core bound state)
    @param imp_gamma,imp_c: non-magnetic impurity scattering rate [eV] and T-matrix c
    @param    field: B/Hc2 (0 = isolated vortex; >0 = circular-cell vortex lattice)
    @param        h: Zeeman (Maki) energy mu*B (spin splitting of the core LDOS)
    """
    if gap_sym != 's' or imp_gamma != 0.0 or field > 0.0 or h != 0.0 or fs_kind is not None:
        fs = None
        if fs_kind is not None:
            from ._eilenberger_fs import build_model_fs
            fs = build_model_fs(fs_kind, 64, params=fs_params)
        return _calc_vortex_dwave(coupling, temp, wc, gap_sym, kb, sw_ldos, Lxi, ngrid,
                                  imp_gamma, imp_c, field, h, fs)
    omega = matsubara(temp, wc)
    print("isolated-vortex quasiclassical Eilenberger (Riccati, s-wave, clean)", flush=True)
    print(f"T={temp/kb:.2f} K, lambda={coupling:.3f}, {len(omega)} Matsubara freqs, "
          f"grid={ngrid}x{ngrid}, cell={Lxi} xi", flush=True)
    rgrid, Dr, Dbulk, xi = solve_vortex(coupling, temp, omega, Lxi=Lxi, ngrid=ngrid)
    if Dbulk <= 0.0:
        print("normal state (Dbulk=0); nothing to profile", flush=True)
        return rgrid, Dr
    print(f"Dbulk = {Dbulk:.6e} eV,  xi = {xi:.4g} (hvf=1 units)", flush=True)
    print(f"Delta(core)/Dbulk = {Dr[0]/Dbulk:.4f}   Delta(edge)/Dbulk = {Dr[-1]/Dbulk:.4f}", flush=True)
    try:
        with open('vortex_gap.dat', 'w') as fh:
            fh.write("# rho/xi   D(rho)/Dbulk   D(rho)[eV]\n")
            for rj, dj in zip(rgrid, Dr):
                fh.write(f"{rj/xi:12.5e} {dj/Dbulk:12.5e} {dj:14.6e}\n")
    except IOError as e:
        print(f"Error: failed to write 'vortex_gap.dat': {e}", flush=True)

    if kappa > 0.0:   # self-consistent finite-kappa Maxwell field profile B(rho)
        rho_s, lam, Brel = vortex_field_profile(rgrid, Dr, Dbulk, xi, kappa, omega)
        print(f"finite-kappa Maxwell: kappa={kappa:.1f}, lambda_bulk={kappa*xi:.4g} ({kappa*xi/xi:.1f} xi); "
              f"B(0)/<B>={Brel[0]:.3f}  rho_s(core)/rho_s_bulk={rho_s[0]:.3f}", flush=True)
        try:
            with open('vortex_field.dat', 'w') as fh:
                fh.write("# rho/xi  rho_s(rho)/rho_s_bulk  lambda(rho)/xi  B(rho)/<B>\n")
                for rj, rsj, lj, bj in zip(rgrid, rho_s, lam, Brel):
                    fh.write(f"{rj/xi:12.5e} {rsj:12.5e} {lj/xi:12.5e} {bj:12.5e}\n")
        except IOError as e:
            print(f"Error: failed to write 'vortex_field.dat': {e}", flush=True)

    if sw_ldos:
        wlist = np.linspace(-3.0 * Dbulk, 3.0 * Dbulk, 301)
        rho_probe = np.array([0.0, 0.5 * xi, 1.0 * xi, 2.0 * xi, 4.0 * xi])
        rho_out, ldos = vortex_ldos(Dr, rgrid, xi, wlist, Dbulk, Lxi=Lxi, ngrid=ngrid,
                                    rho_list=rho_probe)
        n0_core = ldos[0, np.argmin(np.abs(wlist))]
        print(f"core zero-bias LDOS  N(0,0)/N0 = {n0_core:.4f}", flush=True)
        try:
            with open('vortex_ldos.dat', 'w') as fh:
                fh.write("# w/Dbulk  " + "  ".join(f"N(rho={r/xi:.1f}xi)" for r in rho_out) + "\n")
                for iw, wj in enumerate(wlist):
                    fh.write(f"{wj/Dbulk:12.5e} " +
                             " ".join(f"{ldos[ir, iw]:11.5e}" for ir in range(len(rho_out))) + "\n")
        except IOError as e:
            print(f"Error: failed to write 'vortex_ldos.dat': {e}", flush=True)
    return rgrid, Dr


def _calc_vortex_dwave(coupling, temp, wc, gap_sym, kb, sw_ldos, Lxi, ngrid,
                       imp_gamma=0.0, imp_c=1.0e8, field=0.0, h=0.0, fs=None):
    """Driver for the general 2D vortex (non-s pairing, impurities, finite field =
    circular-cell vortex lattice, Zeeman h, and/or a model FS with Fermi velocities):
    full 2D self-consistent Psi(r), the azimuthally-averaged amplitude profile, the
    core spectrum, and the zero-energy LDOS map (the fourfold core star for d-wave;
    Zeeman-split core)."""
    # the 2D solver is far heavier per point than the radial one; use a coarser grid
    ng2d = min(ngrid, 49)
    omega = matsubara(temp, wc)
    imp_txt = 'clean' if imp_gamma == 0 else ('Born' if imp_c > 10 else 'unitary')
    cell_txt = f"lattice B/Hc2={field:.3f}" if field > 0 else f"isolated cell={Lxi} xi"
    if h != 0.0:
        cell_txt += f", Zeeman h={h:.3e} eV"
    if fs is not None:
        cell_txt += ", model FS+v_F"
    print(f"vortex quasiclassical Eilenberger (Riccati, {gap_sym}, {imp_txt}, 2D, {cell_txt})", flush=True)
    print(f"T={temp/kb:.2f} K, lambda={coupling:.3f}, {len(omega)} Matsubara freqs, "
          f"grid={ng2d}x{ng2d}, Gamma_N={imp_gamma:.3e} eV", flush=True)
    xg, Psi, Dbulk, xi = solve_vortex2d(coupling, temp, omega, gap_sym, Lxi=Lxi, ngrid=ng2d,
                                        imp_gamma=imp_gamma, imp_c=imp_c, field=field, fs=fs)
    if Dbulk <= 0.0:
        print("normal state (Dbulk=0); nothing to profile", flush=True)
        return xg, Psi
    X, Y = np.meshgrid(xg, xg, indexing='ij')
    Rg = np.sqrt(X ** 2 + Y ** 2)
    ic = ng2d // 2                                  # core index (x=y=0)
    Redge = xg[-1]
    print(f"Dbulk = {Dbulk:.6e} eV,  xi = {xi:.4g} (hvf=1 units), cell R/xi = {Redge/xi:.2f}", flush=True)
    print(f"|Psi(core)|/Dbulk = {np.abs(Psi[ic, ic])/Dbulk:.4f}   "
          f"|Psi(edge)|/Dbulk = {np.abs(Psi[ic, -1])/Dbulk:.4f}", flush=True)
    # azimuthally-averaged amplitude profile |Psi|(rho)
    rgrid = np.linspace(0.0, Redge, ng2d // 2)
    amp = np.abs(Psi)
    try:
        from scipy.interpolate import RegularGridInterpolator as _RGI
        ai = _RGI((xg, xg), amp, bounds_error=False, fill_value=Dbulk)
        th = np.linspace(0.0, 2.0 * np.pi, 48, endpoint=False)
        with open('vortex_gap.dat', 'w') as fh:
            fh.write("# rho/xi   |Psi|(rho)/Dbulk\n")
            for r in rgrid:
                val = ai(np.stack([r * np.cos(th), r * np.sin(th)], axis=1)).mean()
                fh.write(f"{r/xi:12.5e} {val/Dbulk:12.5e}\n")
    except IOError as e:
        print(f"Error: failed to write 'vortex_gap.dat': {e}", flush=True)

    if sw_ldos:
        # zero-energy map (shows the fourfold star) + core spectrum
        wmap = np.array([0.0])
        nmap = vortex_ldos2d(Psi, xg, xi, wmap, gap_sym, Dbulk, nbeta=48,
                             imp_gamma=imp_gamma, imp_c=imp_c, field=field, h=h, fs=fs)[:, :, 0]
        print(f"core zero-bias LDOS  N(0,0)/N0 = {nmap[ic, ic]:.4f}", flush=True)
        if field > 0.0:
            mask = Rg < Redge
            print(f"cell-averaged zero-energy DOS  <N(0)>/N0 = {nmap[mask].mean():.4f}", flush=True)
        # anisotropy: zero-energy LDOS at rho~xi along node [110] vs antinode [100]
        try:
            from scipy.interpolate import RegularGridInterpolator as _RGI
            ni = _RGI((xg, xg), nmap, bounds_error=False, fill_value=0.0)
            r1 = xi
            n_anti = float(ni([[r1, 0.0]])[0])                       # [100] antinode dir
            n_node = float(ni([[r1 / np.sqrt(2), r1 / np.sqrt(2)]])[0])  # [110] node dir
            print(f"zero-energy LDOS at rho=xi: node[110]/antinode[100] = "
                  f"{n_node:.3f}/{n_anti:.3f} = {n_node/max(n_anti,1e-9):.2f}", flush=True)
        except Exception:
            pass
        try:
            with open('vortex_ldos_map.dat', 'w') as fh:
                fh.write("# x/xi  y/xi  N(x,y,w=0)/N0\n")
                for i in range(ng2d):
                    for j in range(ng2d):
                        fh.write(f"{xg[i]/xi:10.4f} {xg[j]/xi:10.4f} {nmap[i, j]:12.5e}\n")
                    fh.write("\n")
        except IOError as e:
            print(f"Error: failed to write 'vortex_ldos_map.dat': {e}", flush=True)
        # core spectrum N(0,0,w)
        wlist = np.linspace(-3.0 * Dbulk, 3.0 * Dbulk, 121)
        spec = vortex_ldos2d(Psi, xg, xi, wlist, gap_sym, Dbulk, nbeta=36,
                             imp_gamma=imp_gamma, imp_c=imp_c, field=field, h=h, fs=fs)[ic, ic, :]
        try:
            with open('vortex_ldos.dat', 'w') as fh:
                fh.write("# w/Dbulk   N(core,w)/N0\n")
                for wj, nj in zip(wlist, spec):
                    fh.write(f"{wj/Dbulk:12.5e} {nj:12.5e}\n")
        except IOError as e:
            print(f"Error: failed to write 'vortex_ldos.dat': {e}", flush=True)
    return xg, Psi


def calc_vortex_lattice(coupling: float, temp: float, wc: float, gap_sym: str = 's',
                        field_list=None, kb: float = 1.0, imp_gamma: float = 0.0,
                        imp_c: float = 1.0e8, ngrid: int = 37, nbeta: int = 20):
    """
    @fn calc_vortex_lattice
    @brief Sweep the field B/Hc2 (circular-cell vortex lattice) and report the
    cell-averaged zero-energy density of states <N(0)>/N0(B), which distinguishes
    the pairing state: s-wave grows ~linearly with B (core-dominated), d-wave grows
    ~sqrt(B) (Volovik effect, nodal-quasiparticle Doppler shift).  Writes
    'vortex_lattice_dos.dat' and fits the power law <N(0)> ~ B^p.

    @param field_list: list of B/Hc2 values (default a log-spaced set 0.03..0.5)
    @return list of (field, <N(0)>/N0, Dbulk)
    """
    omega = matsubara(temp, wc)
    if field_list is None:
        field_list = [0.03, 0.06, 0.12, 0.25, 0.5]
    print(f"vortex lattice (circular-cell, Doppler) <N(0)>(B): {gap_sym}, lambda={coupling:.3f}, "
          f"T={temp/kb:.2f} K", flush=True)
    results = []
    for b in field_list:
        xg, Psi, Dbulk, xi = solve_vortex2d(coupling, temp, omega, gap_sym, ngrid=ngrid,
                                            nbeta=nbeta, field=b, imp_gamma=imp_gamma, imp_c=imp_c)
        if Dbulk <= 0.0:
            print(f"  B/Hc2={b:.3f}: normal", flush=True)
            continue
        Rc = np.sqrt(2.0 / b) * xi
        nmap = vortex_ldos2d(Psi, xg, xi, np.array([0.0]), gap_sym, Dbulk,
                             nbeta=max(nbeta, 36), field=b, imp_gamma=imp_gamma, imp_c=imp_c)[:, :, 0]
        X, Y = np.meshgrid(xg, xg, indexing='ij')
        mask = np.sqrt(X ** 2 + Y ** 2) < Rc
        n0 = float(nmap[mask].mean())
        results.append((b, n0, Dbulk))
        print(f"  B/Hc2={b:.3f}  Rc/xi={Rc/xi:.2f}  Dbulk={Dbulk:.3e}  <N(0)>/N0={n0:.4f}", flush=True)
    if len(results) >= 2:
        bb = np.array([r[0] for r in results])
        nn = np.array([r[1] for r in results])
        p = np.polyfit(np.log(bb), np.log(np.maximum(nn, 1e-12)), 1)[0]
        print(f"  power law  <N(0)>/N0 ~ (B/Hc2)^{p:.2f}   "
              f"(s-wave ~1 cores, d-wave ~0.5 Volovik)", flush=True)
    try:
        with open('vortex_lattice_dos.dat', 'w') as fh:
            fh.write("# B/Hc2   <N(0)>/N0   Dbulk[eV]\n")
            for b, n0, Db in results:
                fh.write(f"{b:12.5e} {n0:12.5e} {Db:14.6e}\n")
    except IOError as e:
        print(f"Error: failed to write 'vortex_lattice_dos.dat': {e}", flush=True)
    return results
