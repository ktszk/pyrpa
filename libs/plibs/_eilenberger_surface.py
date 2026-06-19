#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Specular surface (1D) quasiclassical Eilenberger solver via the Riccati equation.

This is the first inhomogeneous application of the Riccati machinery in
``_eilenberger`` and the stepping stone toward the vortex lattice: it establishes
the real-space self-consistency of Delta(x), the trajectory<->grid sampling, the
boundary handling, and the observable (LDOS) extraction -- everything the vortex
solver needs except the magnetic field.

Geometry: a superconductor fills the half space x >= 0 with a flat, specularly
reflecting surface at x = 0.  A quasiclassical trajectory along the Fermi
velocity direction beta (measured from the surface normal x_hat) reflects at the
surface, k_x -> -k_x.  Using the *unfolding* trick the incoming branch
(beta_in = pi - beta, moving toward the surface) and the reflected outgoing
branch (beta, moving into the bulk) are concatenated into a single straight
Riccati trajectory whose order parameter form factor flips at the turning point
x = 0; no explicit boundary condition is then needed beyond the bulk inflow at
x = +infinity.

Pairing is the same separable model as the bulk solver, Delta(x, p) =
phi(p) * Delta_amp(x), with phi normalized to <phi^2>_FS = 1.  The Fermi surface
is an isotropic 2D cylinder (v_F || k); the surface orientation enters through
the form-factor angle ``beta_surf`` (0 = [100], pi/4 = [110] for d_{x^2-y^2}).
Feeding the real band's Fermi velocities is a later step; the model FS isolates
the new inhomogeneous physics for validation against known results
(Anderson-flat s-wave, surface-suppressed sign-changing gap, zero-energy
surface bound state).
"""
import numpy as np
from ._eilenberger import riccati_homogeneous, propagators_from_riccati


def form_factor(beta: np.ndarray, gap_sym: str, beta_surf: float = 0.0) -> np.ndarray:
    """
    @fn form_factor
    @brief Pairing form factor phi(beta) on the cylindrical FS, normalized so that
    the full-circle average <|phi|^2> = 1 (so the coupling lambda keeps its bulk
    BCS meaning).  May be complex (chiral / triplet states).
    @param     beta: FS angle(s) from the surface normal [rad]
    @param  gap_sym: singlet: 's', 'd' (d_{x^2-y^2}), 'dxy';
                     triplet (odd parity): 'px', 'py', 'p+ip' / 'p-ip' (chiral)
    @param beta_surf: rotation of the gap relative to the surface normal
                      (d-wave: 0 -> [100] surface, pi/4 -> [110] surface)
    @note Triplet states are treated in the single (pseudo-)spin sector with a
          fixed d-vector: the equal-spin chiral state p+ip = e^{i beta} is fully
          gapped in the bulk and carries topological edge / core states.
    """
    a = beta - beta_surf
    if gap_sym == 's':
        return np.ones_like(beta)
    if gap_sym == 'd':
        return np.sqrt(2.0) * np.cos(2.0 * a)
    if gap_sym == 'dxy':
        return np.sqrt(2.0) * np.sin(2.0 * a)
    if gap_sym == 'px':
        return np.sqrt(2.0) * np.cos(a)
    if gap_sym == 'py':
        return np.sqrt(2.0) * np.sin(a)
    if gap_sym in ('p+ip', 'chiral'):
        return np.exp(1j * a) * np.ones_like(beta, dtype=np.complex128)
    if gap_sym == 'p-ip':
        return np.exp(-1j * a) * np.ones_like(beta, dtype=np.complex128)
    raise ValueError(f"unknown gap_sym: {gap_sym}")


def _integrate_vec(omega: np.ndarray, Delta_s: np.ndarray, hvf: float, ds: float,
                   gamma0: np.ndarray) -> np.ndarray:
    """Integrate the Riccati coherence function along a trajectory, vectorized over
    the frequency axis, with an *unconditionally stable* analytic step.

    The Riccati equation is stiff: high Matsubara frequencies relax to the local
    bulk root over a length hbar*v_F/(2*omega_n) that is far smaller than the grid
    step on grazing trajectories, so explicit RK schemes overflow.  Writing
    gamma = u2/u1 with the linear spinor system
        hbar v_F d/ds [u1; u2] = [[omega, conj(D)], [D, -omega]] [u1; u2]
    the exact propagator over a step with piecewise-constant D maps
        gamma -> (gamma + T (D - omega gamma)) / (1 + T (omega + conj(D) gamma)),
        T = tanh(R t)/R,   R = sqrt(omega^2 + |D|^2),   t = ds/hbar v_F,
    which is bounded for any step size (tanh saturates) and reduces to the correct
    stable root, so the integration never blows up.

    omega and Delta_s may be constant along the path ([Nw]) or position dependent
    ([Ns, Nw]); the latter is needed once the impurity self-energy renormalizes the
    frequency and gap locally, omega -> w_tilde(x), Delta -> Delta_tilde(x).

    @param   omega: (renormalized/complex) frequencies [Nw] or [Ns, Nw]
    @param Delta_s: order parameter along the trajectory [Ns] or [Ns, Nw]
    @param     hvf: hbar |v_F|
    @param      ds: arc-length step
    @param  gamma0: initial gamma at s[0] [Nw]
    @return  gamma: gamma along the trajectory [Ns, Nw]
    """
    Ns = len(Delta_s)
    om2d = (np.ndim(omega) == 2)
    Nw = omega.shape[-1]
    out = np.empty((Ns, Nw), dtype=np.complex128)
    g = np.array(gamma0, dtype=np.complex128)
    out[0] = g
    t = ds / hvf
    for i in range(Ns - 1):
        D = 0.5 * (Delta_s[i] + Delta_s[i + 1])     # piecewise-constant over step
        om = 0.5 * (omega[i] + omega[i + 1]) if om2d else omega
        R = np.sqrt(om ** 2 + (D * np.conj(D)).real)
        T = np.tanh(R * t) / R
        g = (g + T * (D - om * g)) / (1.0 + T * (om + np.conj(D) * g))
        out[i + 1] = g
    return out


def _trajectory_gf(omega: np.ndarray, Damp: np.ndarray, phi_in: float, phi_out: float,
                   hvf: float, dx: float, cosb: float, sigf: np.ndarray = None):
    """
    @fn _trajectory_gf
    @brief Solve g, f on one unfolded surface trajectory for all frequencies.

    The unfolded path runs incoming (x: L->0, momentum beta_in, form factor
    phi_in) then outgoing (x: 0->L, momentum beta, form factor phi_out).  gamma
    is integrated forward from the incoming bulk; gamma-tilde backward from the
    outgoing bulk.  With impurities, the renormalized frequency ``omega`` and the
    additive (isotropic) self-energy ``sigf`` are position dependent on the x grid.

    @param  omega: frequencies [Nw] (clean) or w_tilde(x) [Ngrid, Nw] (with impurity)
    @param   Damp: gap amplitude profile on the x grid [Ngrid] (x=0..L)
    @param phi_in,phi_out: form factor on the incoming / outgoing branch
    @param    hvf: hbar |v_F|
    @param     dx: x-grid spacing
    @param   cosb: |cos(beta)| (perpendicular projection of v_F)
    @param   sigf: additive impurity self-energy Sigma_f(x) [Ngrid, Nw], or None (clean)
    @return (g_out, f_out, g_in, f_in): propagators [Ngrid, Nw] for the outgoing
            (direction beta) and incoming (direction beta_in) branches, each on
            the x=0..L grid.
    """
    Ng = len(Damp)
    ds = dx / cosb
    om2d = (np.ndim(omega) == 2)
    # unfolded order parameter Delta_tilde: incoming (x=L..0) then outgoing (x=0..L)
    if sigf is None:
        D_in = phi_in * Damp[::-1]
        D_out = phi_out * Damp
    else:
        D_in = phi_in * Damp[::-1, None] + sigf[::-1]
        D_out = phi_out * Damp[:, None] + sigf
    Dpath = np.concatenate([D_in, D_out])              # [2Ng] or [2Ng, Nw]
    om_fwd = np.concatenate([omega[::-1], omega]) if om2d else omega
    om_start = om_fwd[0] if om2d else omega
    om_end = om_fwd[-1] if om2d else omega
    # gamma: forward from incoming bulk (first path point)
    ga0, _ = riccati_homogeneous(om_start, Dpath[0])
    gamma = _integrate_vec(om_fwd, Dpath, hvf, ds, ga0)
    # gamma-tilde: backward from outgoing bulk (last path point)
    _, gb0 = riccati_homogeneous(om_end, Dpath[-1])
    om_bwd = om_fwd[::-1] if om2d else omega
    gammat_rev = _integrate_vec(om_bwd, np.conj(Dpath[::-1]), hvf, ds, gb0)
    gammat = gammat_rev[::-1]
    g, f, _ = propagators_from_riccati(gamma, gammat)
    # split: first half = incoming (x=L..0 -> reverse to x=0..L), second = outgoing
    g_in = g[:Ng][::-1]
    f_in = f[:Ng][::-1]
    g_out = g[Ng:]
    f_out = f[Ng:]
    return g_out, f_out, g_in, f_in


def _beta_grid(Nbeta: int, cos_min: float = 0.02):
    """Outgoing FS angles beta in (-pi/2, pi/2), excluding only a thin grazing
    sliver (|cos beta| < cos_min); the incoming partners are pi - beta.

    The exclusion must stay thin: over the full (-pi/2, pi/2) range 2*beta spans a
    whole period so <cos^2 2beta> = <sin^2 2beta> = 1/2 and the analytic sqrt(2)
    normalization of the form factor is orientation independent.  Truncating too
    much biases the half-range average and makes the bulk gap depend on the
    surface orientation.  The unconditionally stable analytic integrator tolerates
    the large step ds = dx/|cos beta| on near-grazing trajectories, so a thin
    sliver is all that is needed."""
    bmax = np.arccos(cos_min)
    beta = np.linspace(-bmax, bmax, Nbeta)
    return beta


def _surface_form_factors(beta: np.ndarray, gap_sym: str, beta_surf: float):
    """Outgoing/incoming form factors on the trajectory direction set, normalized
    so that the mean of |phi|^2 over the full set {beta} U {pi-beta} is 1.  This ties
    the surface kernel and the bulk reference to the *same* angular average, making
    the profile heal to Dbulk and the result independent of surface orientation."""
    phi_out = form_factor(beta, gap_sym, beta_surf)
    phi_in = form_factor(np.pi - beta, gap_sym, beta_surf)
    norm = np.sqrt(np.mean(np.abs(np.concatenate([phi_out, phi_in])) ** 2))
    if norm > 0:
        phi_out = phi_out / norm
        phi_in = phi_in / norm
    return phi_out, phi_in


def solve_surface(coupling: float, temp: float, omega: np.ndarray, gap_sym: str,
                  beta_surf: float = 0.0, Dbulk: float = None, Nbeta: int = 30,
                  Lxi: float = 8.0, nper: int = 16, hvf: float = 1.0,
                  imp_gamma: float = 0.0, imp_c: float = 1.0e8,
                  eps: float = 1.0e-4, itemax: int = 300, mix: float = 0.3):
    """
    @fn solve_surface
    @brief Self-consistently solve the gap profile Delta_amp(x) near a specular
    surface (Matsubara axis), optionally with non-magnetic T-matrix impurities.

    With impurities the self-energy is *local*: at each depth x the renormalized
    frequency w_tilde(x) and the additive pair self-energy Sigma_f(x) are built
    from the local Fermi-surface averages <g(x)>, <f(x)> and fed back into the
    Riccati trajectories.  The gap Delta(x), w_tilde(x) and Sigma_f(x) are relaxed
    together in one fixed-point loop.

    @param coupling: dimensionless separable pairing coupling lambda
    @param     temp: temperature [eV]
    @param    omega: positive Matsubara frequencies [Nw]
    @param  gap_sym: 's', 'd', or 'dxy'
    @param beta_surf: surface orientation (d-wave: 0=[100], pi/4=[110])
    @param    Dbulk: bulk gap amplitude [eV]; if None it is computed self-consistently
                     (with impurity if imp_gamma>0) from the homogeneous gap equation
    @param    Nbeta: number of outgoing trajectory angles
    @param      Lxi: domain length in units of the coherence length xi=hvf/(pi*Dbulk)
    @param     nper: x-grid points per xi
    @param      hvf: hbar |v_F| (length unit; physics depends only on x/xi)
    @param imp_gamma: normal-state impurity scattering rate Gamma_N [eV] (0=clean)
    @param    imp_c: T-matrix cot(delta0) (large=Born, 0=unitary)
    @return (x, Damp, Dbulk): x grid [Ngrid], gap profile [Ngrid], bulk gap
    """
    beta = _beta_grid(Nbeta)
    cosb = np.abs(np.cos(beta))
    phi_out, phi_in = _surface_form_factors(beta, gap_sym, beta_surf)
    phi_all = np.concatenate([phi_out, phi_in])     # full direction set, <phi^2>=1
    clean = (imp_gamma == 0.0)
    Nw = len(omega)

    if Dbulk is None:
        # bulk gap on the SAME (normalized) direction set used by the surface kernel,
        # so the profile heals to exactly Dbulk far from the surface and the result
        # is orientation independent (<phi^2>=1 enforced for every orientation)
        Dbulk = (_bulk_gap(coupling, temp, omega, phi_all, mix=mix) if clean
                 else _bulk_gap_imp(coupling, temp, omega, phi_all, imp_gamma, imp_c, mix=mix))
    # treat an exponentially small bulk gap (T >= Tc, e.g. AG-suppressed) as normal
    if Dbulk < 1.0e-6 * temp:
        xi = hvf / (np.pi * max(temp, 1e-12))
        x = np.linspace(0.0, Lxi * xi, int(Lxi * nper))
        return x, np.zeros_like(x), 0.0

    xi = hvf / (np.pi * Dbulk)
    Ngrid = int(Lxi * nper)
    x = np.linspace(0.0, Lxi * xi, Ngrid)
    dx = x[1] - x[0]
    ndir = 2 * Nbeta

    Damp = np.full(Ngrid, Dbulk, dtype=np.float64)
    gpref = imp_gamma * (imp_c ** 2 + 1.0)
    wt = np.tile(omega.astype(np.complex128), (Ngrid, 1))   # w_tilde(x,n) [Ngrid,Nw]
    sigf = np.zeros((Ngrid, Nw), dtype=np.complex128)       # Sigma_f(x,n)
    for it in range(itemax):
        acc = np.zeros(Ngrid, dtype=np.complex128)          # sum_n <phi f>(x)
        gsum = np.zeros((Ngrid, Nw), dtype=np.complex128)
        fsum = np.zeros((Ngrid, Nw), dtype=np.complex128)
        for ib in range(Nbeta):
            if clean:
                g_out, f_out, g_in, f_in = _trajectory_gf(omega, Damp, phi_in[ib],
                                                          phi_out[ib], hvf, dx, cosb[ib])
            else:
                g_out, f_out, g_in, f_in = _trajectory_gf(wt, Damp, phi_in[ib],
                                                          phi_out[ib], hvf, dx, cosb[ib], sigf=sigf)
            # gap projection uses conj(phi) (matters for complex/chiral form factors)
            acc += (np.conj(phi_out[ib]) * f_out + np.conj(phi_in[ib]) * f_in).sum(axis=1)
            if not clean:
                gsum += g_out + g_in
                fsum += f_out + f_in
        Damp_new = coupling * 2.0 * temp * (acc.real / ndir)
        err = np.abs(Damp_new - Damp).max() / max(Dbulk, 1e-12)
        Damp = (1.0 - mix) * Damp + mix * Damp_new
        if not clean:
            avg_g = gsum / ndir
            avg_f = fsum / ndir
            Dimp = imp_c ** 2 + avg_g ** 2 + avg_f * np.conj(avg_f)
            wt_new = omega[None, :] + gpref * avg_g / Dimp
            sigf_new = gpref * avg_f / Dimp
            err = max(err, np.abs(wt_new - wt).max() / max(omega[0], 1e-12))
            wt = (1.0 - mix) * wt + mix * wt_new
            sigf = (1.0 - mix) * sigf + mix * sigf_new
        if err < eps:
            break
    return x, Damp, Dbulk


def _bulk_gap(coupling: float, temp: float, omega: np.ndarray, phi: np.ndarray,
              eps: float = 1.0e-8, itemax: int = 500, mix: float = 0.5) -> float:
    """Homogeneous gap amplitude for the cylindrical FS (angle average over phi),
    consistent with the surface solver's normalization."""
    damp = 1.764 * temp
    Dk = phi[:, None]
    absphi2 = (phi * np.conj(phi)).real[:, None]            # |phi|^2
    for _ in range(itemax):
        R = np.sqrt(omega[None, :] ** 2 + absphi2 * damp ** 2)
        f = (Dk * damp) / R                                 # [Nbeta, Nw]
        phif = (np.conj(Dk) * f).mean(axis=0)               # angle average <phi* f> [Nw]
        damp_new = coupling * 2.0 * temp * phif.sum().real
        if damp_new <= 0.0:
            return 0.0
        if abs(damp_new - damp) < eps * max(damp_new, eps):
            return damp_new
        damp = (1.0 - mix) * damp + mix * damp_new
    return damp


def _bulk_gap_imp(coupling: float, temp: float, omega: np.ndarray, phi: np.ndarray,
                  gamma: float, cimp: float, eps: float = 1.0e-8, itemax: int = 500,
                  mix: float = 0.5) -> float:
    """Homogeneous gap amplitude on the cylindrical FS with the non-magnetic
    T-matrix impurity self-energy (the bulk reference for the surface solver).
    Nested loop: inner impurity self-consistency at fixed gap, outer gap update."""
    gpref = gamma * (cimp ** 2 + 1.0)
    damp = 1.764 * temp
    wt = omega.astype(np.complex128).copy()
    sigf = np.zeros(len(omega), dtype=np.complex128)
    for _ in range(itemax):
        for _inner in range(200):
            Dt = phi[:, None] * damp + sigf[None, :]
            R = np.sqrt(wt[None, :] ** 2 + Dt * np.conj(Dt))
            ag = (wt[None, :] / R).mean(axis=0)
            af = (Dt / R).mean(axis=0)
            Dimp = cimp ** 2 + ag ** 2 + af * np.conj(af)
            wt_n = omega + gpref * ag / Dimp
            sf_n = gpref * af / Dimp
            if max(np.abs(wt_n - wt).max(), np.abs(sf_n - sigf).max()) < 1e-12:
                wt, sigf = wt_n, sf_n
                break
            wt, sigf = wt_n, sf_n
        Dt = phi[:, None] * damp + sigf[None, :]
        R = np.sqrt(wt[None, :] ** 2 + Dt * np.conj(Dt))
        phif = (np.conj(phi[:, None]) * (Dt / R)).mean(axis=0)
        damp_new = coupling * 2.0 * temp * phif.sum().real
        if damp_new <= 0.0:
            return 0.0
        if abs(damp_new - damp) < eps * max(damp_new, eps):
            return damp_new
        damp = (1.0 - mix) * damp + mix * damp_new
    return damp


def surface_ldos(Damp: np.ndarray, x: np.ndarray, wlist: np.ndarray, gap_sym: str,
                 beta_surf: float = 0.0, ix: int = 0, delta: float = None,
                 Dbulk: float = 1.0, Nbeta: int = 60, hvf: float = 1.0,
                 imp_gamma: float = 0.0, imp_c: float = 1.0e8, h: float = 0.0,
                 eps: float = 1.0e-4, itemax: int = 200, mix: float = 0.5):
    """
    @fn surface_ldos
    @brief Local density of states N(x,w)/N0 = Re<g(x,w)> at grid index ix (default
    the surface x=0) from the converged gap profile, via the retarded Riccati
    equation (iw_n -> w + i*delta).

    With impurities the retarded T-matrix self-energy is elastic (frequency
    diagonal), so a single spatial self-consistency loop -- vectorized over the
    whole real-frequency grid -- gives w_tilde(x,w) and Sigma_f(x,w) at the fixed
    converged Delta(x).  Non-magnetic scattering broadens (and can split) the
    d-wave zero-energy surface bound state.

    Zeeman (Maki) field h splits the spin sectors w -> w -+ h: for a unitary gap
    with a fixed d-vector the spin-matrix Riccati block-diagonalizes into the two
    sigma=+-1 scalar sectors, so N = (1/2) sum_sigma N(w - sigma h).  This splits
    the d[110] zero-energy surface bound state into a peak at +-h.

    @param   Damp: converged gap profile [Ngrid] (from solve_surface)
    @param      x: x grid [Ngrid]
    @param  wlist: real-frequency grid [Nw]
    @param gap_sym,beta_surf: pairing / surface orientation (as in solve_surface)
    @param     ix: grid index at which to evaluate the LDOS (0 = surface)
    @param  delta: retarded broadening [eV] (default 0.02*Dbulk)
    @param  Dbulk: bulk gap (for the default broadening)
    @param  Nbeta: number of outgoing angles
    @param imp_gamma,imp_c: impurity scattering rate and T-matrix c (0=clean)
    @param      h: Zeeman (Maki) energy mu*B (spin splitting of the LDOS)
    @return  ldos: N(x,w)/N0 [Nw]
    """
    if delta is None:
        delta = 0.02 * Dbulk
    beta = _beta_grid(Nbeta)
    cosb = np.abs(np.cos(beta))
    phi_out, phi_in = _surface_form_factors(beta, gap_sym, beta_surf)
    dx = x[1] - x[0]
    Ng, Nw = len(Damp), len(wlist)
    ndir = 2 * Nbeta
    clean = (imp_gamma == 0.0)
    gpref = imp_gamma * (imp_c ** 2 + 1.0)

    def _sector(zomega):
        """LDOS N(ix) for a given (possibly spin-shifted) retarded frequency [Nw]."""
        wt = np.tile(zomega, (Ng, 1))
        sigf = np.zeros((Ng, Nw), dtype=np.complex128)
        for it in range(1 if clean else itemax):
            gsum_x = np.zeros((Ng, Nw), dtype=np.complex128)
            fsum_x = np.zeros((Ng, Nw), dtype=np.complex128)
            for ib in range(Nbeta):
                if clean:
                    g_out, _, g_in, _ = _trajectory_gf(zomega, Damp, phi_in[ib],
                                                       phi_out[ib], hvf, dx, cosb[ib])
                else:
                    g_out, f_out, g_in, f_in = _trajectory_gf(wt, Damp, phi_in[ib],
                                                              phi_out[ib], hvf, dx, cosb[ib], sigf=sigf)
                    fsum_x += f_out + f_in
                gsum_x += g_out + g_in
            if clean:
                return (gsum_x[ix] / ndir).real
            avg_g = gsum_x / ndir
            avg_f = fsum_x / ndir
            Dimp = imp_c ** 2 + avg_g ** 2 + avg_f * np.conj(avg_f)
            wt_new = zomega[None, :] + gpref * avg_g / Dimp
            sigf_new = gpref * avg_f / Dimp
            err = np.abs(wt_new - wt).max() / max(abs(zomega[-1]), 1e-12)
            wt = (1.0 - mix) * wt + mix * wt_new
            sigf = (1.0 - mix) * sigf + mix * sigf_new
            if err < eps:
                break
        return avg_g[ix].real

    base = delta - 1j * wlist                          # retarded frequency [Nw]
    if h == 0.0:
        return _sector(base)
    # spin-resolved: w -> w - sigma h  <=>  zomega -> base + i sigma h
    return 0.5 * (_sector(base + 1j * h) + _sector(base - 1j * h))


def calc_surface(coupling: float, temp: float, wc: float, gap_sym: str = 'd',
                 beta_surf: float = 0.0, kb: float = 1.0, sw_ldos: bool = True,
                 imp_gamma: float = 0.0, imp_c: float = 1.0e8, h: float = 0.0,
                 Nbeta: int = 30, Lxi: float = 8.0, nper: int = 16):
    """
    @fn calc_surface
    @brief Driver: solve the self-consistent gap profile Delta(x) at a specular
    surface and (optionally) the surface LDOS, writing results to files.

    @param coupling: separable pairing coupling lambda
    @param     temp: temperature [eV]
    @param       wc: Matsubara cutoff energy [eV]
    @param  gap_sym: 's', 'd' (d_{x^2-y^2}), or 'dxy'
    @param beta_surf: surface orientation (d-wave: 0=[100], pi/4=[110])
    @param       kb: Boltzmann constant [eV/K] for K-unit reporting
    @param  sw_ldos: also compute the real-frequency surface LDOS (the bound state)
    @param imp_gamma,imp_c: non-magnetic impurity scattering rate [eV] and T-matrix
                            c (0=clean, large=Born, 0=unitary)
    """
    from ._eilenberger import matsubara
    omega = matsubara(temp, wc)
    print("specular-surface quasiclassical Eilenberger (Riccati)", flush=True)
    print(f"gap_sym={gap_sym}, beta_surf={beta_surf:.4f} rad, "
          f"T={temp/kb:.2f} K, lambda={coupling:.3f}, {len(omega)} Matsubara freqs", flush=True)
    print(f"impurity: Gamma_N = {imp_gamma:.4e} eV, c = {imp_c:.3e} "
          f"({'clean' if imp_gamma == 0 else 'Born' if imp_c > 10 else 'unitary'})", flush=True)
    x, Damp, Dbulk = solve_surface(coupling, temp, omega, gap_sym, beta_surf,
                                   imp_gamma=imp_gamma, imp_c=imp_c, Nbeta=Nbeta, Lxi=Lxi, nper=nper)
    if Dbulk <= 0.0:
        print("normal state (Dbulk=0); nothing to profile", flush=True)
        return x, Damp
    xi = 1.0 / (np.pi * Dbulk)   # hvf=1
    print(f"Dbulk = {Dbulk:.6e} eV,  xi = {xi:.4g} (hvf=1 units)", flush=True)
    print(f"Delta(0)/Dbulk = {Damp[0]/Dbulk:.4f}   Delta(L)/Dbulk = {Damp[-1]/Dbulk:.4f}", flush=True)
    try:
        with open('surface_gap.dat', 'w') as fh:
            fh.write("# x/xi   Delta(x)/Dbulk   Delta(x)[eV]\n")
            for xj, dj in zip(x, Damp):
                fh.write(f"{xj/xi:12.5e} {dj/Dbulk:12.5e} {dj:14.6e}\n")
    except IOError as e:
        print(f"Error: failed to write 'surface_gap.dat': {e}", flush=True)

    if sw_ldos:
        wlist = np.linspace(-3.0 * Dbulk, 3.0 * Dbulk, 401)
        ldos = surface_ldos(Damp, x, wlist, gap_sym, beta_surf, ix=0, Dbulk=Dbulk, Nbeta=60,
                            imp_gamma=imp_gamma, imp_c=imp_c, h=h)
        n0 = ldos[np.argmin(np.abs(wlist))]
        if h != 0.0:
            print(f"Zeeman h = {h:.4e} eV (spin-split surface bound state)", flush=True)
        print(f"surface zero-energy LDOS  N(0,0)/N0 = {n0:.4f}", flush=True)
        try:
            with open('surface_ldos.dat', 'w') as fh:
                fh.write("# w/Dbulk   N(x=0,w)/N0\n")
                for wj, nj in zip(wlist, ldos):
                    fh.write(f"{wj/Dbulk:12.5e} {nj:12.5e}\n")
        except IOError as e:
            print(f"Error: failed to write 'surface_ldos.dat': {e}", flush=True)
    return x, Damp
