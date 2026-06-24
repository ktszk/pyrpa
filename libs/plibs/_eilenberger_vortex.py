#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Vortex quasiclassical Eilenberger solver via the Riccati equation: the isolated
vortex (axisymmetric s-wave and the general 2D solver) and, at the end of the
file, the true periodic vortex lattice (Doppler/London, square/triangular).

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
from ._eilenberger import matsubara
from ._eilenberger_surface import _bulk_gap, _bulk_gap_imp
from ..flibs import riccati_chords


def _chords_batch(om3: np.ndarray, dd3: np.ndarray, hvf: float, ds):
    """Batched scalar g, f over all chords [Ns, Nchord, Nw] via the Fortran
    riccati_chords kernel (forward gamma + backward gamma-tilde + combine).  ds is
    scalar or per-chord [Nchord]."""
    return riccati_chords(om3, dd3, hvf, ds)


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
        Delta = _delta_grid(Dr, rgrid, Rg, phase, Dbulk)   # [ngrid, ngrid] (s, impact param)
        om3 = np.broadcast_to(omega, (ngrid, ngrid, Nw)).astype(np.complex128)
        dd3 = np.broadcast_to(Delta[:, :, None], (ngrid, ngrid, Nw))
        _, f_grid = _chords_batch(om3, np.ascontiguousarray(dd3), hvf, dx)   # one call, all chords
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
    Nw = len(wlist)
    om3 = np.broadcast_to(zomega, (ngrid, ngrid, Nw)).astype(np.complex128)
    dd3 = np.broadcast_to(Delta[:, :, None], (ngrid, ngrid, Nw))
    g_grid, _ = _chords_batch(om3, np.ascontiguousarray(dd3), hvf, dx)
    rho_out = rgrid if rho_list is None else np.asarray(rho_list)
    gavg = _azimuthal(g_grid, xg, rho_out, theta, weight_phase=False)    # [Nr, Nw]
    return rho_out, gavg.real


_INT_GAP_STR = {0: 's', 1: 'd', 2: 's', 3: 'dxy', -1: 'px', -2: 'py', -3: 'p+ip'}   # gap_symms index -> continuum


def _ff_vortex(beta: np.ndarray, gap_sym) -> np.ndarray:
    """Pairing form factor phi_d(k_hat) vs the v_F direction beta, normalized so
    <phi^2> = 1 over the full circle.  For a d-wave vortex the gap a trajectory
    sees is phi(beta) * Psi(r), so the trajectory direction enters explicitly and
    the rotational symmetry of the s-wave vortex is reduced to fourfold.  ``gap_sym``
    may be a string or an integer (gap_symms index, mapped to the continuum harmonic)."""
    if isinstance(gap_sym, (int, np.integer)):
        gap_sym = _INT_GAP_STR.get(int(gap_sym), 's')
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
    B used by the periodic-lattice solver below; mixing it with the kept phase is inconsistent.)
    """
    rho2 = Lx ** 2 + Ly ** 2
    rho2s = np.maximum(rho2, rho_min ** 2)
    vdoth = (-Ly * cb + Lx * sb)            # rho * (v_hat . theta_hat)
    delta = -(rho2 / Rc ** 2) * vdoth / (2.0 * rho2s)   # vector-potential part only
    return np.where(rho2 < Rc ** 2, delta, 0.0)


def _doppler_aphi(Lx, Ly, Rc, cb, sb, rho_min, rgrid, aphi):
    """Field Doppler shift -v_F.A from a SELF-CONSISTENT radial vector potential
    A_theta(rho) (finite-kappa Maxwell back-reaction), replacing the uniform-field
    _doppler_chord.  delta = -A_theta(rho) (v_hat . theta_hat); outside the tabulated
    range A_theta ~ 1/rho (flux conservation).  aphi is A_theta on rgrid (0..Rc)."""
    rho = np.sqrt(np.maximum(Lx ** 2 + Ly ** 2, rho_min ** 2))
    vhth = (-Ly * cb + Lx * sb) / rho       # v_hat . theta_hat (unit)
    At = np.interp(rho, rgrid, aphi, right=aphi[-1] * rgrid[-1])  # A_theta(rho)
    far = rho > rgrid[-1]
    At = np.where(far, aphi[-1] * rgrid[-1] / rho, At)            # 1/rho tail beyond the cell
    return np.where(rho ** 2 < Rc ** 2, -At * vhth, 0.0)


def solve_vortex2d(coupling: float, temp: float, omega: np.ndarray, gap_sym: str = 'd',
                   Dbulk: float = None, Lxi: float = 8.0, ngrid: int = 49, nbeta: int = 24,
                   hvf: float = 1.0, imp_gamma: float = 0.0, imp_c: float = 1.0e8,
                   field: float = 0.0, fs: dict = None, qprof=None, eps: float = 2.0e-3,
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
    further below (solve_lattice / lattice_dos; square/triangular, finite kappa).

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
    # representative velocity scale (sets the coherence length xi); per-direction
    # |v_F| is used in each Riccati and Doppler so the field effect is consistent
    hvf_eff = float((fs['nf'] * fs['vabs']).sum()) if use_fs else hvf
    bfull = np.linspace(0.0, 2.0 * np.pi, 180, endpoint=False)
    if Dbulk is None:
        if use_fs:
            from ._eilenberger import bulk_gap_fs
            Dbulk = bulk_gap_fs(coupling, temp, omega, fs, gap_sym)
        else:
            Dbulk = (_bulk_gap(coupling, temp, omega, _ff_vortex(bfull, gap_sym)) if no_imp
                     else _bulk_gap_imp(coupling, temp, omega, _ff_vortex(bfull, gap_sym), imp_gamma, imp_c))
    if Dbulk < 1.0e-6 * temp:
        xi = hvf_eff / (np.pi * max(temp, 1e-12))
        xg = np.linspace(-Lxi * xi, Lxi * xi, ngrid)
        return xg, np.zeros((ngrid, ngrid), dtype=np.complex128), 0.0, xi

    xi = hvf_eff / (np.pi * Dbulk)
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
        from ._eilenberger import fs_form_factor
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
                om3 = np.broadcast_to(omega, (ngrid, ngrid, Nw)).astype(np.complex128)
                dd3 = np.broadcast_to(base[:, :, None], (ngrid, ngrid, Nw))
                _, f3 = _chords_batch(om3, np.ascontiguousarray(dd3), hvf_i, dx)
                sumf = f3.sum(axis=2)
                accf += wt_dir[ib] * np.conj(phi[ib]) * _eval_field(sumf, xg, sxy, bxy, fill=0.0)
            else:
                om_rot = np.broadcast_to(omega, (ngrid, ngrid, Nw)).astype(np.complex128).copy()
                if field > 0.0:                          # Doppler shift v_F.Q (scales with |v_F|)
                    dch = (_doppler_aphi(Lx, Ly, Rc, cb, sb, rho_min, qprof[0], qprof[1])
                           if qprof is not None else _doppler_chord(Lx, Ly, Rc, cb, sb, rho_min))
                    om_rot = om_rot + 1j * hvf_i * dch[:, :, None]
                Dtraj = np.broadcast_to(base[:, :, None], (ngrid, ngrid, Nw)).copy()
                if not no_imp:
                    om_rot = om_rot + _eval_field(dwt, xg, Lx, Ly, fill=0.0)
                    Dtraj = Dtraj + _eval_field(sigf, xg, Lx, Ly, fill=0.0)
                g_rot, f_rot = _chords_batch(om_rot, Dtraj, hvf_i, dx)
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
                  field: float = 0.0, h: float = 0.0, fs: dict = None, qprof=None,
                  eps: float = 3.0e-3, itemax: int = 60, mix: float = 0.5):
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
    if use_fs:
        from ._eilenberger import fs_form_factor
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
        if field <= 0.0:
            dopp.append(None)
        elif qprof is not None:
            dopp.append(_doppler_aphi(Lx, Ly, Rc, cb, sb, rho_min, qprof[0], qprof[1]))
        else:
            dopp.append(_doppler_chord(Lx, Ly, Rc, cb, sb, rho_min))

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
                    om3 = np.broadcast_to(zomega, (ngrid, ngrid, Nw)).astype(np.complex128)
                    dd3 = np.broadcast_to(base[ib][:, :, None], (ngrid, ngrid, Nw))
                    g_rot, _ = _chords_batch(om3, np.ascontiguousarray(dd3), hvfarr[ib], dx)
                    gacc += wt_dir[ib] * _eval_field(g_rot, xg, sxy, bxy, fill=0.0)
                else:
                    om_rot = np.broadcast_to(zomega, (ngrid, ngrid, Nw)).astype(np.complex128).copy()
                    if field > 0.0:
                        om_rot = om_rot + 1j * hvfarr[ib] * dopp[ib][:, :, None]
                    Dtraj = np.broadcast_to(base[ib][:, :, None], (ngrid, ngrid, Nw)).copy()
                    if not no_imp:
                        om_rot = om_rot + _eval_field(dwt, xg, Lx, Ly, fill=0.0)
                        Dtraj = Dtraj + _eval_field(sigf, xg, Lx, Ly, fill=0.0)
                    g_rot, f_rot = _chords_batch(om_rot, Dtraj, hvfarr[ib], dx)
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
                field: float = 0.0, h: float = 0.0, kappa: float = 0.0, tilt_deg: float = 0.0,
                fs_kind: str = None, fs_params=None, fs=None, Lxi: float = 8.0, ngrid: int = 81):
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
    @param tilt_deg: field tilt theta from the c-axis (quasi-2D).  The perpendicular
                     component B_z=B cos(theta) keeps setting the orbital vortex
                     (Doppler) ``field``, while the Zeeman couples to the total |B|,
                     so the effective Maki energy is h/cos(theta) -- a knob for the
                     Pauli/orbital ratio (Pauli-limiting / FFLO with tilt).
    """
    h_eff = h / np.cos(np.radians(tilt_deg)) if tilt_deg else h
    if gap_sym != 's' or imp_gamma != 0.0 or field > 0.0 or h_eff != 0.0 or fs_kind is not None or fs is not None:
        if fs is None and fs_kind is not None:        # build a model FS (else use the prebuilt fs)
            from ._eilenberger import build_model_fs
            fs = build_model_fs(fs_kind, 64, params=fs_params)
        return _calc_vortex_dwave(coupling, temp, wc, gap_sym, kb, sw_ldos, Lxi, ngrid,
                                  imp_gamma, imp_c, field, h_eff, fs)
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


# =========================================================================== #
#  True periodic vortex lattice (square / triangular): Doppler / London form
# --------------------------------------------------------------------------- #
# Gauge-transformed ("B") formulation: the order-parameter phase is removed, so
# the Riccati equation is integrated with a *real* lattice-periodic amplitude
# |Delta(r)| and the field enters only through the gauge-invariant supercurrent
# (Doppler) shift w_tilde(r) = w_n + i v_F . Q(r), with Q(r) periodic on the
# magnetic unit cell.  Q and B(r) are built in Fourier space over the reciprocal
# lattice {K}; finite GL kappa=lambda/xi enters via London screening
#   B(K) = Bbar F(K)/(1+lambda^2 K^2),  Q(K) = -i (zhat x K)(pi/S) F(K) lambda^2/(1+lambda^2 K^2),
# S = Phi0/B (one flux quantum), core form factor F(K)=exp(-xi^2 K^2/2).  This is
# the standard tractable vortex-lattice DOS method (Volovik effect); both square
# and triangular lattices are supported.  Trajectories wrap the cell periodically
# (fractional-coordinate sampling) -- the genuine 2D lattice, not a circular cell.


def _lattice_vectors(field: float, xi: float, lattice: str = 'square', nflux: int = 1):
    """Primitive (a1,a2) and reciprocal (b1,b2) vectors of the (super)cell holding
    ``nflux`` flux quanta (area S = nflux * 2 pi xi^2 / (B/Hc2), so the area per vortex
    -- hence the field -- is unchanged).  ai.bj = 2 pi dij."""
    S = nflux * 2.0 * np.pi / field * xi ** 2
    if lattice.startswith('s'):                       # square
        a = np.sqrt(S)
        a1, a2 = np.array([a, 0.0]), np.array([0.0, a])
    else:                                             # triangular (Abrikosov ground state)
        a = np.sqrt(2.0 * S / np.sqrt(3.0))
        a1, a2 = np.array([a, 0.0]), np.array([0.5 * a, 0.5 * np.sqrt(3.0) * a])
    M = np.column_stack([a1, a2])
    Bm = 2.0 * np.pi * np.linalg.inv(M).T
    return a1, a2, Bm[:, 0], Bm[:, 1], S, M


def _vortex_positions(nflux: int, a1: np.ndarray, a2: np.ndarray):
    """Cartesian positions of the ``nflux`` vortices in the (super)cell, on a regular
    n1 x n2 sub-grid (near-square factorization).  For nflux=1 -> the cell origin; for
    a perfect-square nflux the sub-grid coincides with the primitive vortex lattice so
    the structure factor reduces the supercell back to the single-vortex cell."""
    n1 = int(round(np.sqrt(nflux)))
    while n1 > 1 and nflux % n1 != 0:
        n1 -= 1
    n2 = nflux // n1
    return [(p / n1) * a1 + (q / n2) * a2 for p in range(n1) for q in range(n2)]


def _lattice_QB(X, Y, b1, b2, S, lam, xi, vpos, nmax=8):
    """Periodic supercurrent momentum Q(r)=(Qx,Qy) and field B(r)/Bbar on the grid,
    from the London/Brandt Fourier sum (finite kappa via lambda).  ``vpos`` are the
    Cartesian positions of the (nflux) vortices in the cell; their structure factor
    SF(K)=sum_v e^{-iK.r_v} places them (for one centred vortex SF=1)."""
    Qx = np.zeros_like(X)
    Qy = np.zeros_like(X)
    Brel = np.ones_like(X)                            # B(r)/Bbar (K=0 term = 1)
    nflux = len(vpos)
    pref = np.pi / S
    for m in range(-nmax, nmax + 1):
        for n in range(-nmax, nmax + 1):
            if m == 0 and n == 0:
                continue
            K = m * b1 + n * b2
            K2 = K @ K
            SF = sum(np.exp(-1j * (K[0] * rv[0] + K[1] * rv[1])) for rv in vpos)  # structure factor
            F = np.exp(-0.5 * xi ** 2 * K2)           # core form factor
            scr = lam ** 2 / (1.0 + lam ** 2 * K2)
            ph = np.exp(1j * (K[0] * X + K[1] * Y))
            Qx += (-1j * (-K[1]) * pref * SF * F * scr * ph).real
            Qy += (-1j * (K[0]) * pref * SF * F * scr * ph).real
            Brel += (F / (1.0 + lam ** 2 * K2) * (SF / nflux) * ph).real
    return Qx, Qy, Brel


def _nn_dist(X, Y, a1, a2, vpos, nrep=2):
    """Distance to the nearest vortex (periodic over the cell), for the model |Delta|
    seed/profile.  ``vpos`` are the in-cell vortex positions."""
    d = np.full(X.shape, np.inf)
    for rv in vpos:
        for i in range(-nrep, nrep + 1):
            for j in range(-nrep, nrep + 1):
                Rx = rv[0] + i * a1[0] + j * a2[0]
                Ry = rv[1] + i * a1[1] + j * a2[1]
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


def _sample_pts(interp, x, y, Minv):
    """Sample a periodic interpolator at array points x,y (any shape)."""
    f1 = (Minv[0, 0] * x + Minv[0, 1] * y) % 1.0
    f2 = (Minv[1, 0] * x + Minv[1, 1] * y) % 1.0
    pts = np.stack([f1.ravel(), f2.ravel()], axis=1)
    return interp(pts).reshape(x.shape)


def solve_lattice(coupling, temp, omega, gap_sym='d', field=0.2, kappa=5.0,
                  lattice='square', Dbulk=None, Ng=24, nbeta=18, hvf=1.0, fs=None, nflux=1):
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
            from ._eilenberger import bulk_gap_fs
            Dbulk = bulk_gap_fs(coupling, temp, omega, fs, gap_sym)
        else:
            Dbulk = _bulk_gap(coupling, temp, omega, _ff_vortex(bfull, gap_sym))
    if Dbulk <= 0:
        return None
    xi = hvf_eff / (np.pi * Dbulk)
    a1, a2, b1, b2, S, M = _lattice_vectors(field, xi, lattice, nflux)
    vpos = _vortex_positions(nflux, a1, a2)           # in-cell vortex positions
    Minv = np.linalg.inv(M)
    lam = kappa * xi
    fax = (np.arange(Ng) + 0.5) / Ng                  # cell-centered fractional grid
    F1, F2 = np.meshgrid(fax, fax, indexing='ij')
    X = F1 * a1[0] + F2 * a2[0]
    Y = F1 * a1[1] + F2 * a2[1]
    Qx, Qy, Brel = _lattice_QB(X, Y, b1, b2, S, lam, xi, vpos)
    dnn = _nn_dist(X, Y, a1, a2, vpos)
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
        from ._eilenberger import fs_form_factor
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
        # batch all offsets b for this direction: chord = axis ns, chords = nb offsets
        x = sgrid[:, None] * cb - bgrid[None, :] * sb           # [ns, nb]
        y = sgrid[:, None] * sb + bgrid[None, :] * cb
        Dl = phi[ib] * _sample_pts(Di, x, y, Minv)              # [ns, nb]
        dop = hvfarr[ib] * (cb * _sample_pts(Qix, x, y, Minv)
                            + sb * _sample_pts(Qiy, x, y, Minv))
        om3 = zomega[None, None, :] + 1j * dop[:, :, None]      # [ns, nb, Nw]
        dd3 = np.broadcast_to(Dl[:, :, None] + 0.0j, (ns, nb, Nw))
        g, _ = _chords_batch(np.ascontiguousarray(om3), np.ascontiguousarray(dd3), hvfarr[ib], ds)
        gdir = g[central].sum(axis=(0, 1))                      # sum over central s and offsets
        gsum += wt_dir[ib] * gdir / (nb * ncen)                # spatial mean per direction, nf-weighted
    return gsum.real


def calc_vortex_lattice_periodic(coupling, temp, wc, gap_sym='d', field_list=None,
                                 kappa=5.0, lattice='square', kb=1.0, Ng=20, nbeta=16,
                                 fs_kind=None, fs_params=None, fs=None, nflux=1):
    """
    @fn calc_vortex_lattice_periodic
    @brief Driver: true periodic vortex lattice (Doppler/London, finite kappa).
    Sweeps B/Hc2, self-consistently solves |Delta(r)| and reports the spatially-
    averaged zero-energy DOS <N(0)>/N0(B) (s-wave ~B cores; d-wave ~sqrt(B) Volovik),
    and writes the field profile B(r)/Bbar for the largest field.  With ``fs_kind``
    (model FS) or a prebuilt ``fs`` (e.g. build_wannier_fs) the trajectories run on a
    Fermi surface with real Fermi velocities.
    """
    omega = matsubara(temp, wc)
    if field_list is None:
        field_list = [0.05, 0.1, 0.2, 0.4]
    if fs is None and fs_kind is not None:
        from ._eilenberger import build_model_fs
        fs = build_model_fs(fs_kind, 120, params=fs_params)
    print(f"periodic vortex lattice (Doppler/London): {gap_sym}, {lattice}, kappa={kappa}, "
          f"lambda={coupling:.3f}, T={temp/kb:.2f} K, {nflux} vortices/cell"
          f"{', FS+v_F' if fs is not None else ''}", flush=True)
    results = []
    last = None
    for b in field_list:
        st = solve_lattice(coupling, temp, omega, gap_sym=gap_sym, field=b, kappa=kappa,
                           lattice=lattice, Ng=Ng * (1 if nflux == 1 else int(round(nflux ** 0.5))),
                           nbeta=nbeta, fs=fs, nflux=nflux)
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


def vortex_current2d(Psi, xg, xi, omega, temp, gap_sym='d', nbeta=24, hvf=1.0):
    """
    @fn vortex_current2d
    @brief Equilibrium charge supercurrent density j(r) of a converged vortex order
    parameter, j(r) ~ 2 pi T sum_{w_n>0} < v_hat Im g(r,k_hat,iw_n) >_FS (arb. units;
    prefactor 2 e N0 v_F0 dropped).  g is recomputed at the Matsubara frequencies on
    the rotated chords (same machinery as solve_vortex2d) from the |Psi| amplitude and
    the analytic vortex winding.
    @param Psi: converged complex order parameter field [ngrid,ngrid] (= A e^{i theta})
    @return (jx, jy): current-density components on the (x,y) grid [ngrid,ngrid]
    """
    ngrid = len(xg)
    dx = xg[1] - xg[0]
    X, Y = np.meshgrid(xg, xg, indexing='ij')
    SS, BB = np.meshgrid(xg, xg, indexing='ij')
    Nw = len(omega)
    AmpI = RegularGridInterpolator((xg, xg), np.abs(Psi), bounds_error=False,
                                   fill_value=float(np.abs(Psi).max()))
    dirs = np.linspace(0.0, 2.0 * np.pi, nbeta, endpoint=False)
    phi = _ff_vortex(dirs, gap_sym)
    jx = np.zeros((ngrid, ngrid))
    jy = np.zeros((ngrid, ngrid))
    for ib in range(nbeta):
        cb, sb = np.cos(dirs[ib]), np.sin(dirs[ib])
        Lx = SS * cb - BB * sb
        Ly = SS * sb + BB * cb
        base = phi[ib] * AmpI((Lx, Ly)) * np.exp(1j * np.arctan2(Ly, Lx))
        sxy = X * cb + Y * sb
        bxy = -X * sb + Y * cb
        om3 = np.broadcast_to(omega, (ngrid, ngrid, Nw)).astype(np.complex128)
        dd3 = np.broadcast_to(base[:, :, None], (ngrid, ngrid, Nw))
        g, _ = _chords_batch(om3, np.ascontiguousarray(dd3), hvf, dx)
        gsum = 2.0 * temp * g.sum(axis=2)                       # 2T sum_{n>0} g  [ns,nb]
        g_xy = _eval_field(gsum, xg, sxy, bxy, fill=0.0)        # -> (x,y)
        jx += (1.0 / nbeta) * cb * g_xy.imag
        jy += (1.0 / nbeta) * sb * g_xy.imag
    return jx, jy


def vortex_current_profile(jx, jy, xg, ntheta=72):
    """
    @fn vortex_current_profile
    @brief Azimuthal (circulating) supercurrent j_phi(rho) from the 2D current field,
    averaged on circles of radius rho: j_phi = < -jx sin(phi) + jy cos(phi) >_phi.
    @return (rho, jphi): radius grid and the azimuthal current [Nr]
    """
    R = 0.5 * (xg[-1] - xg[0])
    jxi = RegularGridInterpolator((xg, xg), jx, bounds_error=False, fill_value=0.0)
    jyi = RegularGridInterpolator((xg, xg), jy, bounds_error=False, fill_value=0.0)
    rho = np.linspace(0.0, 0.85 * R, len(xg) // 2)
    th = np.linspace(0.0, 2.0 * np.pi, ntheta, endpoint=False)
    jphi = np.empty(len(rho))
    for k, r in enumerate(rho):
        pts = np.stack([r * np.cos(th), r * np.sin(th)], axis=1)
        jphi[k] = np.mean(-jxi(pts) * np.sin(th) + jyi(pts) * np.cos(th))
    return rho, jphi


def calc_vortex_current(coupling, temp, wc, gap_sym='d', kb=1.0, Lxi=8.0, ngrid=49, nbeta=24):
    """
    @fn calc_vortex_current
    @brief Driver: converge an isolated vortex and report its circulating charge
    supercurrent j_phi(rho) (arb. units), which vanishes at the core, peaks near a
    coherence length, and decays outward.  Writes 'vortex_current.dat'.
    """
    omega = matsubara(temp, wc)
    print(f"vortex supercurrent (Riccati Eilenberger, {gap_sym}): T={temp/kb:.2f} K, "
          f"lambda={coupling:.3f}, grid={min(ngrid,49)}^2", flush=True)
    xg, Psi, Dbulk, xi = solve_vortex2d(coupling, temp, omega, gap_sym, Lxi=Lxi, ngrid=min(ngrid, 49), nbeta=nbeta)
    if Dbulk <= 0.0:
        print("normal state (Dbulk=0); no supercurrent", flush=True)
        return xg, None
    jx, jy = vortex_current2d(Psi, xg, xi, omega, temp, gap_sym, nbeta=nbeta)
    rho, jphi = vortex_current_profile(jx, jy, xg)
    ipk = int(np.argmax(np.abs(jphi)))
    print(f"Dbulk={Dbulk:.4e} eV, xi={xi:.4g}; peak |j_phi|={abs(jphi[ipk]):.3e} at rho/xi={rho[ipk]/xi:.2f}, "
          f"j_phi(core)={jphi[0]:.2e}", flush=True)
    try:
        with open('vortex_current.dat', 'w') as fh:
            fh.write("# rho/xi   j_phi[arb]\n")
            for r, jp in zip(rho, jphi):
                fh.write(f"{r/xi:10.4f} {jp:14.6e}\n")
    except IOError as e:
        print(f"Error writing vortex_current.dat: {e}", flush=True)
    return rho, jphi


def _radial_amp(Psi, xg, rgrid, ntheta=48):
    """Azimuthal average of |Psi| on circles of radius rgrid (radial gap profile)."""
    ai = RegularGridInterpolator((xg, xg), np.abs(Psi), bounds_error=False,
                                 fill_value=float(np.abs(Psi).max()))
    th = np.linspace(0.0, 2.0 * np.pi, ntheta, endpoint=False)
    return np.array([ai(np.stack([r * np.cos(th), r * np.sin(th)], axis=1)).mean() for r in rgrid])


def calc_vortex_maxwell(coupling, temp, wc, gap_sym='s', field=0.2, kappa=2.0, kb=1.0,
                        Lxi=8.0, ngrid=41, nbeta=24, itemax_a=6):
    """
    @fn calc_vortex_maxwell
    @brief Circular-cell vortex with the SELF-CONSISTENT (finite-kappa Maxwell)
    vector potential A_theta(rho) -- the je 'A_renew' back-reaction.  Outer loop:
    (i) solve the order parameter with the supercurrent Doppler from the current
    A_theta; (ii) from the converged |Delta|(rho) get the screened field B(rho)
    (finite-kappa London, vortex_field_profile) and integrate it to the new
    A_theta(rho)=(1/rho) int_0^rho B rho' drho'.  Finite kappa screens the
    supercurrent (A_theta below the extreme type-II rho/2Rc^2), reducing the
    field-induced zero-energy DOS.  Reports <N(0)> for the self-consistent A and,
    for reference, the extreme type-II (kappa->infinity) result.  Writes 'vortex_maxwell.dat'.
    """
    omega = matsubara(temp, wc)
    Dbulk = (_bulk_gap(coupling, temp, omega, _ff_vortex(np.linspace(0, 2 * np.pi, 180, endpoint=False), gap_sym)))
    if Dbulk <= 0:
        print("normal state (Dbulk=0)", flush=True)
        return None
    xi = 1.0 / (np.pi * Dbulk)
    Rc = np.sqrt(2.0 / field) * xi
    rgrid = np.linspace(0.0, Rc, ngrid // 2)
    Bbar = 1.0 / Rc ** 2
    aphi = rgrid / (2.0 * Rc ** 2)                      # seed: extreme type-II A_theta
    print(f"vortex Maxwell back-reaction (self-consistent A, finite kappa={kappa}): "
          f"{gap_sym}, B/Hc2={field}, T={temp/kb:.2f} K", flush=True)
    ng2 = min(ngrid, 41)
    for ita in range(itemax_a):
        xg, Psi, Db, xiv = solve_vortex2d(coupling, temp, omega, gap_sym, Lxi=Lxi, ngrid=ng2,
                                          nbeta=nbeta, field=field, qprof=(rgrid, aphi))
        Dr = _radial_amp(Psi, xg, rgrid)
        _, _, Brel = vortex_field_profile(rgrid, Dr, Db, xiv, kappa, omega)   # B(rho)/<B>, mean 1
        B = Brel * Bbar
        integ = np.concatenate([[0.0], np.cumsum(0.5 * (B[1:] * rgrid[1:] + B[:-1] * rgrid[:-1])
                                                 * np.diff(rgrid))])
        aphi_new = integ / np.maximum(rgrid, 0.5 * rgrid[1])           # A_theta=(1/rho)int B rho' drho'
        da = np.abs(aphi_new - aphi).max() / max(aphi.max(), 1e-30)
        aphi = 0.5 * aphi + 0.5 * aphi_new
        if da < 1e-3:
            break
    ic = len(xg) // 2
    n0_sc = vortex_ldos2d(Psi, xg, xiv, np.array([0.0]), gap_sym, Db, nbeta=48,
                          field=field, qprof=(rgrid, aphi))[ic, ic, 0]
    # reference: extreme type-II (no screening)
    xg2, Psi2, Db2, xiv2 = solve_vortex2d(coupling, temp, omega, gap_sym, Lxi=Lxi, ngrid=ng2,
                                          nbeta=nbeta, field=field)
    n0_ex = vortex_ldos2d(Psi2, xg2, xiv2, np.array([0.0]), gap_sym, Db2, nbeta=48,
                          field=field)[ic, ic, 0]
    a_ext = rgrid / (2.0 * Rc ** 2)                    # extreme type-II A_theta
    im = len(rgrid) // 2                                # mid-radius rho ~ Rc/2
    print(f"Dbulk={Db:.4e}, xi={xiv:.4g}, lambda/xi=kappa={kappa} (finite kappa concentrates B near the core)", flush=True)
    print(f"  A_theta(Rc/2): self-consistent={aphi[im]:.4e} vs extreme type-II={a_ext[im]:.4e} "
          f"(ratio {aphi[im]/max(a_ext[im],1e-30):.3f})", flush=True)
    print(f"  A_theta(edge): self={aphi[-1]:.4e} vs extreme={a_ext[-1]:.4e} (enclosed flux fixed)", flush=True)
    print(f"  N(core,0)/N0: self-consistent A = {n0_sc:.3f},  extreme type-II (kappa->inf) = {n0_ex:.3f}", flush=True)
    try:
        with open('vortex_maxwell.dat', 'w') as fh:
            fh.write("# rho/xi   A_theta(self)   A_theta(extreme)   B(rho)/<B>\n")
            for r, a in zip(rgrid, aphi):
                fh.write(f"{r/xiv:10.4f} {a:12.5e} {r/(2*Rc**2):12.5e} "
                         f"{np.interp(r,rgrid,Brel):12.5e}\n")
    except IOError as e:
        print(f"Error writing vortex_maxwell.dat: {e}", flush=True)
    return rgrid, aphi, n0_sc, n0_ex


# --------------------------------------------------------------------------- #
# je-style self-consistent PERIODIC vortex lattice (formulation A)
#
# This mirrors the reference Julia code (je): the order parameter keeps its
# vortex phase (complex Delta), the magnetic unit cell carries one flux quantum
# via the analytic Abrikosov quasi-periodic phase, and the gap equation is closed
# point by point -- for every grid point a trajectory is anchored AT that point
# and the Riccati f is read at the anchor (no scatter/binning, no rotated-frame
# remapping).  Only the real amplitude |Psi|(r) is interpolated (it is smooth and
# strictly periodic); the winding phase is evaluated analytically along the chord,
# so there is no magnetic-translation-phase bookkeeping in the accumulation.
# Extreme type-II (A=0): the supercurrent/Volovik shift emerges from the phase
# variation of the complex Delta sampled along each chord.
# --------------------------------------------------------------------------- #
def _cell_geom(field, xi, lattice='square', nflux=1):
    """Magnetic-cell geometry (je convention): rectangular generators (Dlx, Dly) with
    obliqueness zeta so the primitive translations are a1=((1-zeta)Dlx, -Dly),
    a2=(zeta Dlx, Dly) and the cell holds ``nflux`` flux quanta (area Dlx*Dly =
    nflux*2 pi xi^2 / field).  square -> (1,1,0); triangular (Abrikosov ground state)
    -> (1, sqrt3/2, 1/2).  Returns a dict with Dlx,Dly,zeta, a1,a2, reciprocal b1,b2, S."""
    if lattice.startswith('s'):
        Dlx0, Dly0, zeta = 1.0, 1.0, 0.0
    else:                                              # triangular / hexagonal
        Dlx0, Dly0, zeta = 1.0, np.sqrt(3.0) * 0.5, 0.5
    Lsc = xi * np.sqrt(nflux * 2.0 * np.pi / (field * Dlx0 * Dly0))
    Dlx, Dly = Dlx0 * Lsc, Dly0 * Lsc
    a1 = np.array([(1.0 - zeta) * Dlx, -Dly])
    a2 = np.array([zeta * Dlx, Dly])
    M = np.column_stack([a1, a2])
    B = 2.0 * np.pi * np.linalg.inv(M).T
    return dict(Dlx=Dlx, Dly=Dly, zeta=zeta, a1=a1, a2=a2,
                b1=B[:, 0], b2=B[:, 1], S=abs(Dlx * Dly), nflux=nflux)


def _to_frac(x, y, g):
    """Cartesian -> fractional cell coordinates (r1, r2) (inverse of the je map
    x=(r1(1-z)+r2 z)Dlx, y=(-r1+r2)Dly)."""
    return x / g['Dlx'] - g['zeta'] * y / g['Dly'], x / g['Dlx'] + (1.0 - g['zeta']) * y / g['Dly']


def _abrikosov_z(x, y, g, np_sum=6):
    """Abrikosov (lowest-Landau-level) quasi-periodic order parameter for the magnetic
    cell ``g`` (je's Sabrikosov; Dr=Dly/Dlx, obliqueness zeta, one zero per primitive
    cell at Dx0=-0.5(1+zeta), Dy0=-0.5).  Its phase is the vortex winding, its modulus
    ~1 between vortices and 0 at each core; works for square and triangular cells."""
    Dlx, Dly, zeta = g['Dlx'], g['Dly'], g['zeta']
    Dr = Dly / Dlx
    Dx0, Dy0 = -0.5 * (1.0 + zeta), -0.5
    xl, yl = x / Dlx, y / Dly
    Z = np.zeros(np.broadcast(x, y).shape, dtype=np.complex128)
    for p in range(-np_sum, np_sum + 1):
        Z += (np.exp(-np.pi * Dr * (yl + Dy0 + p) ** 2)
              * np.exp(-2j * np.pi * (p * (Dx0 + p * zeta * 0.5) + (Dy0 + p) * xl)))
    Z *= np.sqrt(np.sqrt(2.0 * Dr))
    return Z * np.exp(-1j * np.pi * xl * yl)


def _abrikosov_unit_phase(x, y, g, Vw=1):
    """Unit phase factor e^{i chi(r)} = (conj(Abrikosov)/|Abrikosov|)^Vw  (conjugate
    convention, je's conj(Zphase)^Vw); Vw>1 gives a multiply-quantized (giant) vortex."""
    Z = _abrikosov_z(x, y, g)
    az = np.abs(Z)
    u = np.where(az > 1e-12, np.conj(Z) / np.where(az > 1e-12, az, 1.0), 1.0 + 0.0j)
    return u if Vw == 1 else u ** Vw


def _periodic_eval(field, g, px, py):
    """Bilinear interpolation of a cell-centred field [Ng,Ng] (grid points at fractional
    (k+0.5)/Ng-0.5) at scattered Cartesian points (px,py), with periodic wrap on the
    (oblique) cell.  Real field only (the amplitude)."""
    Ng = field.shape[0]
    f1, f2 = _to_frac(px, py, g)
    u = ((f1 + 0.5) % 1.0) * Ng - 0.5                  # grid point k sits at u=k
    v = ((f2 + 0.5) % 1.0) * Ng - 0.5
    i0 = np.floor(u).astype(int); j0 = np.floor(v).astype(int)
    wu = u - i0; wv = v - j0
    i0m = i0 % Ng; i1 = (i0 + 1) % Ng
    j0m = j0 % Ng; j1 = (j0 + 1) % Ng
    return ((1 - wu) * (1 - wv) * field[i0m, j0m] + wu * (1 - wv) * field[i1, j0m]
            + (1 - wu) * wv * field[i0m, j1] + wu * wv * field[i1, j1])


def _london_A(g, lam, Vw=1, nfft=64):
    """Finite-kappa (London) vector potential A(r) on the magnetic cell ``g``, as the
    low-pass-filtered bare supervelocity:  A(K) = (grad(chi)/2)(K) / (1 + lambda^2 K^2),
    K = m b1 + n b2.  The filter makes A SMOOTH and regular at the core (so -v_F.A does
    NOT cancel the node-enforcing winding), keeps the uniform piece (K=0 -> removes the
    spurious bare uniform circulation) and screens the periodic supercurrent over lambda.
    Returned on the cell-centred fractional grid; sample with _periodic_eval.
    @return (Ax, Ay) [nfft,nfft]
    """
    fax = (np.arange(nfft) + 0.5) / nfft - 0.5
    F1, F2 = np.meshgrid(fax, fax, indexing='ij')
    X = F1 * g['a1'][0] + F2 * g['a2'][0]
    Y = F1 * g['a1'][1] + F2 * g['a2'][1]
    h = 1e-3 * np.sqrt(g['S'])
    U = _abrikosov_unit_phase(X, Y, g, Vw)
    gx = 0.5 * np.imag(np.conj(U) * (_abrikosov_unit_phase(X + h, Y, g, Vw)
                                     - _abrikosov_unit_phase(X - h, Y, g, Vw)) / (2 * h))
    gy = 0.5 * np.imag(np.conj(U) * (_abrikosov_unit_phase(X, Y + h, g, Vw)
                                     - _abrikosov_unit_phase(X, Y - h, g, Vw)) / (2 * h))
    m = np.fft.fftfreq(nfft) * nfft                    # integer Fourier indices
    M1, M2 = np.meshgrid(m, m, indexing='ij')
    KX = M1 * g['b1'][0] + M2 * g['b2'][0]
    KY = M1 * g['b1'][1] + M2 * g['b2'][1]
    H = 1.0 / (1.0 + lam ** 2 * (KX ** 2 + KY ** 2))
    Ax = np.fft.ifft2(np.fft.fft2(gx) * H).real
    Ay = np.fft.ifft2(np.fft.fft2(gy) * H).real
    return Ax, Ay


def _maxwell_A(jx, jy, g):
    """Coulomb-gauge vector potential from a 2D current density on the cell grid:
    A(K) = j_T(K) / K^2 with j_T the transverse (divergence-free) projection
    j_T(K) = j(K) - K (K.j)/K^2  (je Svector_potential).  K = m b1 + n b2.
    @return (Ax, Ay) on the same [Ng,Ng] grid."""
    Ng = jx.shape[0]
    m = np.fft.fftfreq(Ng) * Ng
    M1, M2 = np.meshgrid(m, m, indexing='ij')
    KX = M1 * g['b1'][0] + M2 * g['b2'][0]
    KY = M1 * g['b1'][1] + M2 * g['b2'][1]
    K2 = KX ** 2 + KY ** 2
    safe = K2 > 1e-12
    K2s = np.where(safe, K2, 1.0)
    Jx, Jy = np.fft.fft2(jx), np.fft.fft2(jy)
    Kdj = (KX * Jx + KY * Jy) / K2s
    Ax = np.where(safe, (Jx - KX * Kdj) / K2s, 0.0)
    Ay = np.where(safe, (Jy - KY * Kdj) / K2s, 0.0)
    return np.fft.ifft2(Ax).real, np.fft.ifft2(Ay).real


def solve_lattice_sc(coupling, temp, omega, gap_sym='d', field=0.2, lattice='square',
                     Ng=20, nbeta=16, hvf=1.0, fs=None, kappa=None, Vw=1, Lchord=6.0,
                     ds_xi=0.3, itemax=60, mix=0.4, eps=3.0e-3, anderson=True, m_and=4,
                     seed_profile=None, self_consistent_A=False, mixA=0.3):
    """
    @fn solve_lattice_sc
    @brief Self-consistent complex order parameter Psi(r)=|Psi|(r) e^{i chi(r)} of a
    true periodic vortex lattice, je-style: formulation A (phase kept in Delta) with
    the analytic Abrikosov winding chi(r), closed by per-grid-point anchored
    trajectories (the f-to-grid map is exact, not binned).
    @param lattice: 'square' or 'triangular' (Abrikosov ground state; oblique cell).
    @param Vw: vortex winding number / flux quanta per cell (1 = single vortex;
    Vw>1 = a multiply-quantized giant vortex, area scaled to hold Vw quanta).
    @param kappa: None -> bare extreme limit (Doppler = analytic grad(chi)/2 only, no
    uniform-A subtraction); finite -> the je finite-kappa back-reaction, where the
    smooth London vector potential A(r)=lowpass(grad chi/2) (_london_A) is subtracted
    as -v_F.A: the uniform vector potential is removed and the supercurrent is
    London-screened (lambda=kappa*xi), reducing the Volovik DOS.
    @param anderson,m_and: Anderson/Pulay acceleration of the self-consistency (m_and
    history vectors); falls back to linear ``mix`` for the first step / if disabled.
    @param seed_profile: optional [Ng,Ng] starting amplitude in units of the bulk gap
    (e.g. a converged absD/Dbulk from a nearby field) for a warm start.
    @param self_consistent_A: finite kappa only -- replace the analytic London A by the
    FULLY self-consistent vector potential (je A_renew): each step the actual
    quasiclassical supercurrent j_s(r)=<v_F Im g> is computed and A(K)=j_{s,T}(K)/K^2
    (transverse, Coulomb gauge, _maxwell_A) is updated together with Delta (mixA).  This
    captures the spatial superfluid-density variation (current suppressed at cores) that
    the bare-grad(chi) London A misses; reduces to it when |Psi| is uniform.
    @return state dict (X,Y,absD,Psi,chi,Dbulk,xi,a,...) for lattice_dos_sc
    """
    bfull = np.linspace(0.0, 2.0 * np.pi, 180, endpoint=False)
    hvf_eff = float((fs['nf'] * fs['vabs']).sum()) if fs is not None else hvf
    if fs is not None:
        from ._eilenberger import bulk_gap_fs, fs_form_factor
        Dbulk = bulk_gap_fs(coupling, temp, omega, fs, gap_sym)
    else:
        Dbulk = _bulk_gap(coupling, temp, omega, _ff_vortex(bfull, gap_sym))
    if Dbulk <= 0:
        return None
    xi = hvf_eff / (np.pi * Dbulk)
    g = _cell_geom(field, xi, lattice, nflux=Vw)        # oblique magnetic cell (Vw quanta)
    a = np.sqrt(g['S'])
    fax = (np.arange(Ng) + 0.5) / Ng - 0.5             # cell-centred fractional grid
    F1, F2 = np.meshgrid(fax, fax, indexing='ij')
    X = F1 * g['a1'][0] + F2 * g['a2'][0]
    Y = F1 * g['a1'][1] + F2 * g['a2'][1]
    chi_grid = np.angle(_abrikosov_unit_phase(X, Y, g, Vw))
    if fs is not None:
        dirs = np.arctan2(fs['vy'], fs['vx']); phi = fs_form_factor(fs, gap_sym)
        hvfarr = np.asarray(fs['vabs']); wt_dir = np.asarray(fs['nf']); nbd = len(dirs)
    else:
        dirs = np.linspace(0.0, 2.0 * np.pi, nbeta, endpoint=False)
        phi = _ff_vortex(dirs, gap_sym); hvfarr = np.full(nbeta, hvf)
        wt_dir = np.full(nbeta, 1.0 / nbeta); nbd = nbeta
    L = Lchord * xi; ds = ds_xi * xi
    ns = int(2 * L / ds) | 1                            # odd -> exact centre index
    ic = ns // 2
    s = np.linspace(-L, L, ns)
    ax = X.ravel(); ay = Y.ravel(); Nanch = ax.size
    lam = kappa * xi if kappa is not None else None
    Axg, Ayg = _london_A(g, lam, Vw) if lam is not None else (None, None)  # smooth screened A(r)
    # pre-compute the fixed per-direction chord geometry, winding phase, projector and
    # (finite kappa) the screening Doppler -v_hat.A(r)
    chord = []
    for ib in range(nbd):
        cb, sb = np.cos(dirs[ib]), np.sin(dirs[ib])
        px = ax[None, :] + s[:, None] * cb             # [ns, Nanch]
        py = ay[None, :] + s[:, None] * sb
        eichi = _abrikosov_unit_phase(px, py, g, Vw)   # e^{i chi} along the chord
        proj = wt_dir[ib] * np.conj(phi[ib]) * np.conj(eichi[ic])   # strip phase at anchor
        dch = (-(cb * _periodic_eval(Axg, g, px, py) + sb * _periodic_eval(Ayg, g, px, py))
               if lam is not None else None)
        chord.append((px, py, eichi, proj, dch))
    if seed_profile is not None:                       # warm start (units of Dbulk)
        A = Dbulk * np.asarray(seed_profile)
    else:
        A = Dbulk * np.tanh(np.abs(_abrikosov_z(X, Y, g)) / 0.5)   # node at each core
    Nw = len(omega)
    om0 = np.broadcast_to(omega, (ns, Nanch, Nw)).astype(np.complex128)
    scA = bool(self_consistent_A and lam is not None)
    if scA:                                            # bulk superfluid density (London calibration)
        q0 = 0.05 * Dbulk / max(hvf_eff, 1e-12)
        jb = 0.0
        for ib in range(nbd):
            wt = omega + 1j * hvfarr[ib] * np.cos(dirs[ib]) * q0
            Dk2 = (phi[ib] * np.conj(phi[ib])).real * Dbulk ** 2
            jb += wt_dir[ib] * hvfarr[ib] * np.cos(dirs[ib]) * np.imag(
                wt / np.sqrt(wt ** 2 + Dk2)).sum()
        rho_s = -2.0 * temp * jb / q0
        Cj = -1.0 / (lam ** 2 * rho_s) if rho_s != 0 else 0.0
        Axf, Ayf = _london_A(g, lam, Vw, nfft=Ng)      # seed A from the London solution (Ng grid)
        AuX, AuY = Axf.mean(), Ayf.mean()              # uniform (K=0) part: cancels the chi gauge term
                                                       # (gauge invariance); _maxwell_A gives only K!=0
    # iteration-independent per-direction frequency arrays (fixed Doppler; recomputed if scA)
    om_dir = [om0 if dch is None
              else np.ascontiguousarray(omega[None, None, :] + 1j * hvfarr[ib] * dch[:, :, None])
              for ib, (_, _, _, _, dch) in enumerate(chord)]

    def gap_map(Af, om_use):                            # one pass: returns A_new (+ current grids if scA)
        accf = np.zeros(Nanch, dtype=np.complex128)
        curx = np.zeros(Nanch); cury = np.zeros(Nanch)
        for ib in range(nbd):
            px, py, eichi, proj, _ = chord[ib]
            base = phi[ib] * _periodic_eval(Af, g, px, py) * eichi
            dd3 = np.broadcast_to(base[:, :, None], (ns, Nanch, Nw))
            gch, f = riccati_chords(om_use[ib], np.ascontiguousarray(dd3), hvfarr[ib], ds)
            accf += proj * f[ic].sum(axis=1)
            if scA:                                    # supercurrent <v_F Im g> at the anchor
                img = gch[ic].imag.sum(axis=1)
                cb, sb = np.cos(dirs[ib]), np.sin(dirs[ib])
                curx += wt_dir[ib] * hvfarr[ib] * cb * img
                cury += wt_dir[ib] * hvfarr[ib] * sb * img
        Anew = np.maximum(coupling * 2.0 * temp * accf.real, 0.0).reshape(Ng, Ng)
        if scA:
            return Anew, 2.0 * temp * curx.reshape(Ng, Ng), 2.0 * temp * cury.reshape(Ng, Ng)
        return Anew

    xs, fs = [], []                                    # Anderson/Pulay history (Walker-Ni) on Delta
    for it in range(itemax):
        if scA:                                        # Doppler from the current self-consistent A
            omd = [np.ascontiguousarray(omega[None, None, :] - 1j * hvfarr[ib]
                   * (np.cos(dirs[ib]) * _periodic_eval(Axf, g, chord[ib][0], chord[ib][1])
                      + np.sin(dirs[ib]) * _periodic_eval(Ayf, g, chord[ib][0], chord[ib][1]))[:, :, None])
                   for ib in range(nbd)]
            Anew, curx, cury = gap_map(A, omd)
            Axn, Ayn = _maxwell_A(Cj * curx, Cj * cury, g)   # K!=0 screening part
            Axn += AuX; Ayn += AuY                            # + uniform (gauge-cancelling) part
            errA = max(np.abs(Axn - Axf).max(), np.abs(Ayn - Ayf).max()) / max(np.abs(Axf).max(), 1e-12)
            Axf = (1.0 - mixA) * Axf + mixA * Axn
            Ayf = (1.0 - mixA) * Ayf + mixA * Ayn
        else:
            Anew = gap_map(A, om_dir); errA = 0.0
        res = Anew - A                                  # residual f_k = G(x_k) - x_k
        err = max(np.abs(res).max() / Dbulk, errA)
        if err < eps:
            A = Anew
            break
        xs.append(A.ravel().copy()); fs.append(res.ravel().copy())
        if len(fs) > m_and + 1:
            xs.pop(0); fs.pop(0)
        if not anderson or len(fs) == 1:
            A_next = A.ravel() + mix * res.ravel()
        else:                                          # min || f_k - dF.gamma ||
            dF = np.column_stack([fs[i + 1] - fs[i] for i in range(len(fs) - 1)])
            dX = np.column_stack([xs[i + 1] - xs[i] for i in range(len(xs) - 1)])
            gam, *_ = np.linalg.lstsq(dF, fs[-1], rcond=None)
            A_next = xs[-1] + mix * fs[-1] - (dX + mix * dF) @ gam
        A = np.maximum(A_next.reshape(Ng, Ng), 0.0)
    return dict(X=X, Y=Y, absD=A, Psi=A * np.exp(1j * chi_grid), chi=chi_grid,
                Dbulk=Dbulk, xi=xi, a=a, S=g['S'], geom=g, b1=g['b1'], b2=g['b2'],
                a1=g['a1'], a2=g['a2'], acell=a, kappa=kappa, Vw=Vw,
                Afield=((Axf, Ayf) if scA else None), iters=it + 1, err=err)


def lattice_dos_sc(state, gap_sym, wlist, delta=None, nbeta=24, hvf=1.0, fs=None,
                   Lchord=6.0, ds_xi=0.3):
    """
    @fn lattice_dos_sc
    @brief Spatially-averaged DOS N(w)/N0 of the je-style self-consistent periodic
    vortex lattice (formulation A): the same per-grid-point anchored trajectories as
    solve_lattice_sc, now on the retarded axis (z=delta-i*w).  The Volovik shift is
    carried by the winding phase of the complex Delta sampled along each chord, plus
    (finite kappa, from state['kappa']) the screening Doppler -v_F.A(r) from the smooth
    London vector potential (_london_A) -- consistent with solve_lattice_sc.
    N(w)/N0 = <Re g(anchor)>_{grid, FS}.
    @return N(w)/N0 [Nw]
    """
    A = state['absD']; xi = state['xi']
    Dbulk = state['Dbulk']; Ng = A.shape[0]
    X, Y = state['X'], state['Y']
    g = state['geom']; Vw = state.get('Vw', 1)
    kappa = state.get('kappa'); lam = kappa * xi if kappa is not None else None
    if state.get('Afield') is not None:                # self-consistent A (je A_renew)
        Axg, Ayg = state['Afield']
    else:
        Axg, Ayg = _london_A(g, lam, Vw) if lam is not None else (None, None)
    if delta is None:
        delta = 0.03 * Dbulk
    zomega = delta - 1j * np.asarray(wlist)
    Nw = zomega.size
    if fs is not None:
        from ._eilenberger import fs_form_factor
        dirs = np.arctan2(fs['vy'], fs['vx']); phi = fs_form_factor(fs, gap_sym)
        hvfarr = np.asarray(fs['vabs']); wt_dir = np.asarray(fs['nf']); nbd = len(dirs)
    else:
        dirs = np.linspace(0.0, 2.0 * np.pi, nbeta, endpoint=False)
        phi = _ff_vortex(dirs, gap_sym); hvfarr = np.full(nbeta, hvf)
        wt_dir = np.full(nbeta, 1.0 / nbeta); nbd = nbeta
    L = Lchord * xi; ds = ds_xi * xi
    ns = int(2 * L / ds) | 1; ic = ns // 2
    s = np.linspace(-L, L, ns)
    ax = X.ravel(); ay = Y.ravel(); Nanch = ax.size
    om0 = np.broadcast_to(zomega, (ns, Nanch, Nw)).astype(np.complex128)
    gsum = np.zeros(Nw)
    for ib in range(nbd):
        cb, sb = np.cos(dirs[ib]), np.sin(dirs[ib])
        px = ax[None, :] + s[:, None] * cb
        py = ay[None, :] + s[:, None] * sb
        eichi = _abrikosov_unit_phase(px, py, g, Vw)
        base = phi[ib] * _periodic_eval(A, g, px, py) * eichi
        dd3 = np.broadcast_to(base[:, :, None], (ns, Nanch, Nw))
        om3 = (om0 if lam is None else
               zomega[None, None, :] - 1j * hvfarr[ib]
               * (cb * _periodic_eval(Axg, g, px, py)
                  + sb * _periodic_eval(Ayg, g, px, py))[:, :, None])
        gg, _ = riccati_chords(np.ascontiguousarray(om3), np.ascontiguousarray(dd3), hvfarr[ib], ds)
        gsum += wt_dir[ib] * gg[ic].real.mean(axis=0)   # <Re g> over anchors
    return gsum


def calc_vortex_lattice_sc(coupling, temp, wc, gap_sym='d', field_list=None,
                           lattice='square', kb=1.0, Ng=20, nbeta=16, fs=None,
                           kappa=None, Vw=1, self_consistent_A=False):
    """
    @fn calc_vortex_lattice_sc
    @brief Driver: je-style self-consistent PERIODIC vortex lattice (formulation A).
    Sweeps B/Hc2, self-consistently solves the complex order parameter Psi(r) (true
    node at every core, full Abrikosov-lattice supercurrent) and reports the
    spatially-averaged zero-energy DOS <N(0)>/N0(B) (d-wave ~sqrt(B) Volovik).
    @param lattice: 'square' or 'triangular' (Abrikosov ground state).
    @param Vw: flux quanta per cell (1; Vw>1 = a multiply-quantized giant vortex).
    @param kappa: None -> bare extreme limit; finite -> the je finite-kappa A(r)
    back-reaction (London-screened supercurrent, smooth vector potential), which
    removes the spurious uniform-A overcount and screens the Volovik DOS.
    @param self_consistent_A: finite kappa only -- use the fully self-consistent A from
    the actual quasiclassical current (je A_renew) instead of the analytic London A.
    Writes the converged |Psi|(r) for the largest field to lattice_sc_op.dat.
    """
    omega = matsubara(temp, wc)
    if field_list is None:
        field_list = [0.05, 0.1, 0.2, 0.35]
    akind = ('extreme type-II' if kappa is None else
             f"finite kappa={kappa}{', self-consistent A' if self_consistent_A else ' (London A)'}")
    print(f"self-consistent periodic vortex lattice (formulation A, {akind}): "
          f"{gap_sym}, {lattice}, Vw={Vw}, lambda={coupling:.3f}, T={temp/kb:.2f} K"
          f"{', FS+v_F' if fs is not None else ''}", flush=True)
    results = []
    last = None
    seed = None                                        # warm start from the previous field
    for b in field_list:
        st = solve_lattice_sc(coupling, temp, omega, gap_sym=gap_sym, field=b,
                              lattice=lattice, Ng=Ng, nbeta=nbeta, fs=fs, kappa=kappa,
                              Vw=Vw, seed_profile=seed, self_consistent_A=self_consistent_A)
        if st is None or st['Dbulk'] <= 0:
            print(f"  B/Hc2={b:.3f}: normal", flush=True)
            continue
        seed = st['absD'] / st['Dbulk']                # dimensionless profile (same Ng/grid)
        n0 = float(lattice_dos_sc(st, gap_sym, np.array([0.0]), nbeta=max(nbeta, 36),
                                  delta=0.03 * st['Dbulk'], fs=fs)[0])
        results.append((b, n0))
        D = st['absD']; Db = st['Dbulk']
        print(f"  B/Hc2={b:.3f}  a/xi={st['a']/st['xi']:.2f}  iters={st['iters']}  "
              f"|D| min/mean/max/Db=[{D.min()/Db:.3f},{D.mean()/Db:.3f},{D.max()/Db:.3f}]  "
              f"<N(0)>/N0={n0:.4f}  N(0)/sqrt(B)={n0/np.sqrt(b):.3f}", flush=True)
        last = st
    if len(results) >= 2:
        bb = np.array([r[0] for r in results]); nn = np.array([r[1] for r in results])
        p = np.polyfit(np.log(bb), np.log(np.maximum(nn, 1e-12)), 1)[0]
        print(f"  power law  <N(0)>/N0 ~ (B/Hc2)^{p:.2f}  (d ~0.5 Volovik, saturating near Hc2)", flush=True)
    try:
        with open('lattice_sc.dat', 'w') as fh:
            fh.write("# B/Hc2   <N(0)>/N0\n")
            for b, n0 in results:
                fh.write(f"{b:12.5e} {n0:12.5e}\n")
        if last is not None:
            with open('lattice_sc_op.dat', 'w') as fh:
                fh.write("# x/xi  y/xi  |Psi|/Dbulk\n")
                xi = last['xi']; Db = last['Dbulk']
                for i in range(last['X'].shape[0]):
                    for j in range(last['X'].shape[1]):
                        fh.write(f"{last['X'][i,j]/xi:10.4f} {last['Y'][i,j]/xi:10.4f} "
                                 f"{last['absD'][i,j]/Db:10.5f}\n")
                    fh.write("\n")
    except OSError:
        pass
    return results
