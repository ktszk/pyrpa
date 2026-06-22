#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Multi-orbital quasiclassical Eilenberger equation solver.

Phase 1 implements the *homogeneous* (spatially uniform) limit on the
Matsubara axis in pure Python.  In a uniform system the gradient term of the
Eilenberger transport equation vanishes, so the equation reduces to the
algebraic normalization condition for the quasiclassical propagators

    g = w_tilde / R,   f = D_tilde / R,   R = sqrt(w_tilde^2 + |D_tilde|^2)

evaluated band-resolved on the Fermi surface.  Non-magnetic impurity
scattering enters through a T-matrix self-energy (which contains the Born and
unitary limits as the parameters c -> infinity and c -> 0):

    w_tilde   = w_l + Gamma * <g> / Dimp
    D_tilde_k = Delta_k + Gamma * <f> / Dimp
    Dimp      = c^2 + <g>^2 + |<f>|^2

with <.> the (DOS-weighted) Fermi-surface average.  For a sign-changing gap
<f> -> 0, so the gap-renormalizing channel disappears and the impurity term
becomes pair breaking, reproducing the Anderson-theorem / Abrikosov-Gor'kov
physics for multiband superconductors.

The pairing interaction is currently a *separable* model,
V(k,k') = lambda * phi(k) phi(k') with phi the gap-symmetry form factor
(``gap_symms``), normalized to <phi^2>_FS = 1 so that ``lambda`` is the
dimensionless BCS-like coupling.  The Fermi-surface averaging and the gap
self-consistency are factored out from the propagator update so that:

  * the separable kernel can be replaced by an FS-projected RPA/FLEX vertex,
  * the homogeneous algebraic normalization can be replaced by a Riccati
    transport solver along trajectories (vortex lattice / inhomogeneous case),
  * the hot inner loops can be moved to Fortran,

without touching the high-level orchestration in ``calc_eilenberger``.
"""
import numpy as np
from scipy.optimize import brentq
from ._bands import get_emesh
from ._response import gap_symms


def _fs_average(values: np.ndarray, w: np.ndarray) -> np.ndarray:
    """DOS-weighted Fermi-surface average over the FS-point axis (axis 0).
    @param values: array [Nfs, ...] of the quantity to average
    @param      w: FS weights [Nfs]
    @return       averaged array [...] (axis 0 contracted)
    """
    return np.tensordot(w, values, axes=(0, 0)) / w.sum()


def build_fs(eig: np.ndarray, klist: np.ndarray, mu: float, gap_sym: int,
             width: float, w_cut: float = 1.0e-4):
    """
    @fn build_fs
    @brief Select Fermi-surface points from a band mesh and build their DOS
    weights and (normalized) pairing form factor.

    A broadened delta function delta(eps_{k,n} - mu) (Gaussian of width
    ``width``) provides the DOS weight, so no explicit FS triangulation is
    needed and all bands are handled uniformly.  Points with negligible weight
    are discarded because the quasiclassical propagators live on the FS.

    @param    eig: band energies on the mesh [Nk, Norb]
    @param  klist: k-points in fractional coordinates [Nk, 3]
    @param     mu: chemical potential [eV]
    @param gap_sym: gap-symmetry index passed to ``gap_symms``
    @param  width: Gaussian broadening of the FS delta function [eV]
    @param  w_cut: keep points with weight > w_cut * max(weight)
    @return (wf, phif): FS weights [Nfs] and normalized form factor [Nfs]
    """
    # Gaussian-broadened delta(eps - mu); the overall constant cancels in averages.
    de = eig - mu
    w = np.exp(-0.5 * (de / width) ** 2)            # [Nk, Norb]
    # form factor depends on k only; broadcast to every band
    phi_row = gap_symms(klist, 1, gap_sym)[0]        # [Nk]
    phi = np.repeat(phi_row[:, None], eig.shape[1], axis=1)  # [Nk, Norb]
    mask = w > w_cut * w.max()
    wf = w[mask]
    phif = phi[mask]
    # normalize so that <phi^2>_FS = 1  -> lambda is the dimensionless coupling
    norm = np.sqrt(_fs_average(phif ** 2, wf))
    if norm > 0:
        phif = phif / norm
    return wf, phif


def matsubara(temp: float, wc: float, Nw_max: int = 2000000) -> np.ndarray:
    """Positive fermionic Matsubara frequencies w_l = (2l+1) pi T below a *fixed*
    energy cutoff w_c.  The number of frequencies grows as 1/T, which keeps the
    pairing cutoff (and hence Tc) well defined; without a fixed cutoff the BCS
    log never saturates and Tc diverges.
    @param temp: temperature [eV]
    @param   wc: Matsubara cutoff energy [eV]
    @param Nw_max: hard cap on the number of frequencies (memory safety)
    """
    Nw = int(np.floor((wc / (np.pi * temp) - 1.0) / 2.0)) + 1
    Nw = max(1, min(Nw, Nw_max))
    return (2 * np.arange(Nw) + 1) * np.pi * temp


def _homogeneous_gf(wt: np.ndarray, Dt: np.ndarray, method: str = 'normalization'):
    """Homogeneous quasiclassical propagators (g, f) from the renormalized
    frequency w_tilde [Nw] and gap D_tilde [Nfs, Nw].

    Both routes are mathematically identical in the uniform limit; ``riccati``
    exercises the same coherence-function parametrization used by the trajectory
    solver (so it doubles as a consistency check of the Riccati machinery).
    @param method: 'normalization' (g=w/R, f=D/R) or 'riccati' (via gamma roots)
    """
    if method == 'riccati':
        gam, gamt = riccati_homogeneous(wt[None, :], Dt)
        g, f, _ = propagators_from_riccati(gam, gamt)
        return g, f
    R = np.sqrt(wt[None, :] ** 2 + Dt * np.conj(Dt))
    return wt[None, :] / R, Dt / R


def solve_quasiclassical(omega: np.ndarray, Delta: np.ndarray, wf: np.ndarray,
                         gamma: float, cimp: float, method: str = 'normalization',
                         eps: float = 1.0e-10, itemax: int = 2000):
    """
    @fn solve_quasiclassical
    @brief Solve the homogeneous Eilenberger equation with a self-consistent
    non-magnetic T-matrix impurity self-energy, for a *fixed* gap.

    @param  omega: positive Matsubara frequencies [Nw]
    @param  Delta: gap on the FS points (real or complex) [Nfs]
    @param     wf: FS DOS weights [Nfs]
    @param  gamma: normal-state scattering rate Gamma_N = n_imp/(pi N0) [eV]
                   (Gamma_N=0 -> clean).  This is the physical, c-independent rate:
                   in the normal state (g=1, f=0) it gives w_tilde = w + Gamma_N
                   in both the Born and unitary limits.
    @param   cimp: T-matrix cotangent parameter c = cot(delta0)
                   (large -> Born, 0 -> unitary).  The internal T-matrix prefactor
                   is Gamma = Gamma_N (c^2 + 1) so that the normal-state rate is
                   held fixed as c is varied.
    @param method: 'normalization' (algebraic g=w/R) or 'riccati' (coherence-
                   function route, shared with the vortex-lattice trajectory solver)
    @param    eps: convergence tolerance on (w_tilde, Sigma_f)
    @param itemax: maximum impurity iterations
    @return (g, f): quasiclassical propagators [Nfs, Nw] complex
    """
    Delta = Delta[:, None]                       # [Nfs, 1]
    if gamma == 0.0:
        wt = omega.astype(np.complex128)
        return _homogeneous_gf(wt, Delta, method)
    # T-matrix prefactor that holds the normal-state rate fixed across c (Born<->unitary)
    gpref = gamma * (cimp ** 2 + 1.0)
    # self-consistent impurity loop; w_tilde and Sigma_f are k-independent scalars(w_l)
    wt = omega.astype(np.complex128).copy()      # [Nw]
    sigf = np.zeros(omega.shape[0], dtype=np.complex128)
    for _ in range(itemax):
        Dt = Delta + sigf[None, :]               # [Nfs, Nw]
        g, f = _homogeneous_gf(wt, Dt, method)
        avg_g = _fs_average(g, wf)               # [Nw]
        avg_f = _fs_average(f, wf)
        Dimp = cimp ** 2 + avg_g ** 2 + avg_f * np.conj(avg_f)
        wt_new = omega + gpref * avg_g / Dimp
        sigf_new = gpref * avg_f / Dimp
        err = max(np.abs(wt_new - wt).max(), np.abs(sigf_new - sigf).max())
        wt, sigf = wt_new, sigf_new
        if err < eps:
            break
    Dt = Delta + sigf[None, :]
    return _homogeneous_gf(wt, Dt, method)


# --------------------------------------------------------------------------- #
#  Riccati formulation (vortex-lattice / inhomogeneous-ready)
# --------------------------------------------------------------------------- #
# The quasiclassical Green's function is parametrized by two coherence functions
# gamma (a) and gamma-tilde (b):
#       g  = (1 - gamma*gammat)/(1 + gamma*gammat)
#       f  = 2 gamma /(1 + gamma*gammat)
#       ft = 2 gammat/(1 + gamma*gammat)
# In a uniform system the Eilenberger transport equation reduces to the algebraic
# (stationary) Riccati root; in an inhomogeneous system (surface, SN junction,
# vortex lattice) gamma/gammat obey first-order ODEs along straight quasiclassical
# trajectories R(s) = R0 + s * v_F^hat, integrated with the renormalized fields
# w_tilde(R), Delta_tilde(R):
#       hbar v_F d(gamma)/ds = Delta_t - 2 w_t gamma - conj(Delta_t) gamma^2
# gamma is the stable (decaying) solution along +v_F; gammat along -v_F.  This is
# the building block the vortex-lattice solver will call per trajectory, sampling
# Delta(R) (with the magnetic phase) on a real-space grid under quasi-periodic
# boundary conditions.

def riccati_homogeneous(omega, Delta):
    """
    @fn riccati_homogeneous
    @brief Stable stationary root of the homogeneous Riccati equation,
    gamma = (R - w)/conj(Delta),  gammat = (R - w)/Delta,  R = sqrt(w^2+|Delta|^2).
    Reproduces g = w/R and f = Delta/R exactly (see propagators_from_riccati).
    @param omega: (renormalized) Matsubara frequency, broadcastable to Delta
    @param Delta: (renormalized) gap, complex array
    @return (gamma, gammat) with the shape of the broadcast inputs
    """
    omega = np.asarray(omega, dtype=np.complex128)
    Delta = np.asarray(Delta, dtype=np.complex128)
    R = np.sqrt(omega ** 2 + Delta * np.conj(Delta))
    absD = np.abs(Delta)
    safe = absD > 0.0
    gamma = np.zeros(np.broadcast(omega, Delta).shape, dtype=np.complex128)
    gammat = np.zeros_like(gamma)
    # (R - w)/conj(D); use where to avoid 0/0 in the normal (Delta=0) region
    gamma = np.where(safe, (R - omega) / np.where(safe, np.conj(Delta), 1.0), 0.0)
    gammat = np.where(safe, (R - omega) / np.where(safe, Delta, 1.0), 0.0)
    return gamma, gammat


def propagators_from_riccati(gamma, gammat):
    """
    @fn propagators_from_riccati
    @brief Quasiclassical propagators from the coherence functions.
    @return (g, f, ft): normal g, anomalous f and f-tilde, same shape as inputs
    """
    den = 1.0 + gamma * gammat
    return (1.0 - gamma * gammat) / den, 2.0 * gamma / den, 2.0 * gammat / den


def gap_kernel(f: np.ndarray, phif: np.ndarray, wf: np.ndarray, temp: float,
               coupling: float) -> complex:
    """
    @fn gap_kernel
    @brief Apply the separable pairing kernel to the anomalous propagator to get
    the updated gap amplitude.

    For V(k,k') = lambda phi(k) phi(k') the gap equation
    Delta_k = sum_{k'} V(k,k') T sum_l f_{k'l} collapses to
    Delta_k = phi_k * Damp with
        Damp = lambda * T sum_{all l} <phi f>_FS
             = lambda * 2 T sum_{l>=0} Re <phi f>_FS
    (the factor 2 accounts for negative Matsubara frequencies via f(-w)=f(w)*).

    @param      f: anomalous propagator on FS points [Nfs, Nw]
    @param   phif: normalized form factor on FS points [Nfs]
    @param     wf: FS DOS weights [Nfs]
    @param   temp: temperature [eV]
    @param coupling: dimensionless pairing coupling lambda
    @return  Damp: updated (scalar) gap amplitude
    """
    phif_w_avg = _fs_average(phif[:, None] * f, wf)   # [Nw]
    return coupling * 2.0 * temp * np.sum(phif_w_avg).real


def pairing_eigenvalue(temp: float, wf: np.ndarray, phif: np.ndarray,
                       omega: np.ndarray, gamma: float, cimp: float,
                       coupling: float, method: str = 'normalization',
                       h: float = 0.0, damp_probe: float = 1.0e-6) -> float:
    """
    @fn pairing_eigenvalue
    @brief Linearized Eilenberger pairing eigenvalue lambda_lin(T): the gap grows
    when lambda_lin > 1, so Tc is defined by lambda_lin(Tc) = 1.
    Computed as the ratio Damp_new / Damp in the Damp -> 0 limit.
    @param h: spin-singlet Zeeman (Maki) energy mu*B; shifts w_n -> w_n + i h so
              the two spin channels (sigma=+-1) are w_n +- i h.  Pauli depairing.
    """
    Delta = phif * damp_probe
    _, f = solve_quasiclassical(omega + 1j * h, Delta, wf, gamma, cimp, method)
    return gap_kernel(f, phif, wf, temp, coupling) / damp_probe


def solve_gap(temp: float, wf: np.ndarray, phif: np.ndarray, omega: np.ndarray,
              gamma: float, cimp: float, coupling: float, method: str = 'normalization',
              h: float = 0.0, damp_init: float = None, eps: float = 1.0e-8, itemax: int = 500,
              mix: float = 0.5) -> float:
    """
    @fn solve_gap
    @brief Self-consistently solve the non-linear homogeneous Eilenberger gap
    equation at temperature ``temp``; returns the gap amplitude (0 if normal).

    @param   damp_init: initial gap amplitude (default 1.764 * temp, BCS ratio)
    @param        mix: linear-mixing factor for the amplitude fixed point
    @param          h: spin-singlet Zeeman (Maki) energy mu*B (Pauli paramagnetic
                       depairing): the gap kernel uses Re f(w_n + i h), i.e. the
                       average of the sigma=+-1 spin channels w_n +- i h.  Gives the
                       Pauli-limited (Chandrasekhar-Clogston) suppression of the
                       singlet gap and the metastable branch up to the spinodal.
    @return     damp: converged gap amplitude [eV]
    """
    if damp_init is None:
        damp_init = 1.764 * temp
    omega_h = omega + 1j * h
    damp = damp_init
    for _ in range(itemax):
        Delta = phif * damp
        _, f = solve_quasiclassical(omega_h, Delta, wf, gamma, cimp, method)
        damp_new = gap_kernel(f, phif, wf, temp, coupling)
        if damp_new <= 0.0:
            return 0.0
        if abs(damp_new - damp) < eps * max(damp_new, eps):
            return damp_new
        damp = (1.0 - mix) * damp + mix * damp_new
    return damp


def find_tc(wf: np.ndarray, phif: np.ndarray, gamma: float, cimp: float,
            coupling: float, wc: float, t_lo: float, t_hi: float,
            method: str = 'normalization', tol: float = 1.0e-4, itemax: int = 60) -> float:
    """
    @fn find_tc
    @brief Locate Tc by bisection on lambda_lin(T) - 1 = 0 within [t_lo, t_hi].
    Returns 0.0 if no superconducting solution is bracketed (lambda_lin(t_lo) < 1).
    @param wc: fixed Matsubara cutoff energy [eV] (sets the pairing energy scale)
    """
    def lam(T):
        return pairing_eigenvalue(T, wf, phif, matsubara(T, wc), gamma, cimp, coupling, method)

    lo, hi = t_lo, t_hi
    if lam(lo) < 1.0:
        return 0.0
    if lam(hi) > 1.0:
        # widen upper bound a few times
        for _ in range(10):
            hi *= 1.5
            if lam(hi) < 1.0:
                break
        else:
            return hi
    for _ in range(itemax):
        mid = 0.5 * (lo + hi)
        if lam(mid) > 1.0:
            lo = mid
        else:
            hi = mid
        if (hi - lo) < tol * hi:
            break
    return 0.5 * (lo + hi)


def dos_zeeman(wlist: np.ndarray, Damp: float, wf: np.ndarray, phif: np.ndarray,
               h: float, delta_br: float) -> np.ndarray:
    """
    @fn dos_zeeman
    @brief Homogeneous (bulk) density of states N(w)/N0 of a singlet superconductor
    in a Zeeman field, via the retarded quasiclassical g of the two spin channels
    (sigma=+-1) with energies w -> w -+ h.  Shows the Zeeman splitting of the
    coherence peaks (Delta -> Delta +- h).
    @param  Damp: bulk gap amplitude [eV]
    @param wf,phif: Fermi-surface weights / form factor (|Delta(k)|=phif*Damp)
    @param     h: Zeeman (Maki) energy mu*B
    @param delta_br: retarded broadening [eV]
    @return N(w)/N0 [Nw]
    """
    Dk = (phif * Damp)[:, None]                         # [Nfs, 1] gap on FS
    Wn = wf.sum()
    out = np.zeros(len(wlist))
    for sigma in (+1.0, -1.0):
        z = (delta_br - 1j * wlist) + 1j * sigma * h    # spin channel [Nw]
        g = z[None, :] / np.sqrt(z[None, :] ** 2 + Dk * np.conj(Dk))   # [Nfs, Nw]
        out += 0.5 * (np.tensordot(wf, g.real, axes=(0, 0)) / Wn)
    return out


def calc_pauli_limit(Nx: int, Ny: int, Nz: int, wc: float, ham_r, S_r, rvec, avec,
                     mu: float, temp: float, gap_sym: int, coupling: float,
                     h_list=None, fs_width: float = 5.0e-3, kb: float = 1.0):
    """
    @fn calc_pauli_limit
    @brief Sweep the Zeeman (Maki) field h and report the singlet gap Delta(h) (the
    metastable SC branch) -- Pauli paramagnetic depairing.  At low T the SC solution
    survives up to the spinodal h* ~ Delta0, while the thermodynamic first-order
    (Chandrasekhar-Clogston) transition is at h_P = Delta0/sqrt(2).  Also writes the
    Zeeman-split bulk DOS at h = 0.5*Delta0 to 'pauli_dos.dat'.
    @param Nx,Ny,Nz,ham_r,S_r,rvec,avec,mu: Fermi-surface inputs (as in calc_eilenberger)
    @param h_list: list of Zeeman energies h [eV] (default fractions of Delta0)
    """
    omega = matsubara(temp, wc)
    Nk, klist, eig, uni, kweight = get_emesh(Nx, Ny, Nz, ham_r, S_r, rvec, avec, sw_uni=True)
    wf, phif = build_fs(eig, klist, mu, gap_sym, fs_width)
    D0 = solve_gap(temp, wf, phif, omega, 0.0, 1.0e8, coupling, h=0.0)
    if D0 <= 0:
        print("normal state at h=0; nothing to do", flush=True)
        return
    print(f"Pauli limiting: Delta0 = {D0:.4e} eV, T = {temp/kb:.2f} K", flush=True)
    print(f"  Chandrasekhar-Clogston h_P = Delta0/sqrt(2) = {D0/np.sqrt(2):.4e} eV", flush=True)
    if h_list is None:
        h_list = D0 * np.array([0.0, 0.3, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9, 1.0, 1.05])
    h_sp = 0.0
    try:
        fh = open('pauli_gap.dat', 'w'); fh.write("# h/Delta0   Delta(h)/Delta0\n")
        for h in h_list:
            # start from Delta0 to follow the metastable SC branch up to the spinodal
            Dh = solve_gap(temp, wf, phif, omega, 0.0, 1.0e8, coupling, h=float(h), damp_init=D0)
            fh.write(f"{h/D0:10.4f} {Dh/D0:10.4f}\n")
            print(f"  h/Delta0={h/D0:5.2f}: Delta/Delta0 = {Dh/D0:.4f}", flush=True)
            if Dh > 0.1 * D0:                 # robust SC threshold (ignore tiny residuals)
                h_sp = h
        fh.close()
    except IOError as e:
        print(f"Error writing pauli_gap.dat: {e}", flush=True)
    print(f"  spinodal h* (SC branch collapses) ~ {h_sp/D0:.2f} Delta0", flush=True)
    # Zeeman-split DOS at h = 0.5 Delta0
    wl = np.linspace(-3 * D0, 3 * D0, 301)
    Nz = dos_zeeman(wl, D0, wf, phif, 0.5 * D0, 0.03 * D0)
    try:
        with open('pauli_dos.dat', 'w') as f2:
            f2.write("# w/Delta0   N(w)/N0  (Zeeman split, h=0.5 Delta0)\n")
            for w, n in zip(wl, Nz):
                f2.write(f"{w/D0:10.4f} {n:12.5e}\n")
    except IOError as e:
        print(f"Error writing pauli_dos.dat: {e}", flush=True)


def calc_eilenberger(Nx: int, Ny: int, Nz: int, wc: float, ham_r, S_r, rvec, avec,
                     mu: float, temp: float, gap_sym: int, coupling: float,
                     imp_gamma: float = 0.0, imp_c: float = 1.0e8,
                     fs_width: float = 5.0e-3, kb: float = 1.0, method: str = 'normalization',
                     sw_find_tc: bool = False, sw_imp_sweep: bool = False,
                     imp_sweep: np.ndarray = None):
    """
    @fn calc_eilenberger
    @brief High-level driver for the homogeneous multi-orbital quasiclassical
    Eilenberger calculation (Matsubara, Python).

    Builds the Fermi surface from the band mesh, reports the linearized pairing
    eigenvalue lambda_lin(T), self-consistently solves the non-linear gap at the
    requested temperature, and optionally locates Tc and/or sweeps the impurity
    scattering rate to map Tc(Gamma) (pair-breaking curve).

    @param   Nx,Ny,Nz: k-mesh for the Fermi-surface average
    @param         wc: fixed Matsubara cutoff energy [eV] (sets the pairing scale;
                       the number of frequencies grows as 1/T so Tc stays finite)
    @param      ham_r,S_r,rvec,avec: tight-binding model / lattice
    @param         mu: chemical potential [eV]
    @param       temp: temperature [eV]
    @param    gap_sym: gap-symmetry index (see gap_symms)
    @param   coupling: dimensionless separable pairing coupling lambda
    @param  imp_gamma: impurity scattering strength Gamma [eV]
    @param      imp_c: T-matrix cot(delta0) (large -> Born, 0 -> unitary)
    @param   fs_width: Gaussian FS broadening [eV]
    @param         kb: Boltzmann constant [eV/K] for K-unit reporting
    @param     method: homogeneous (g,f) route: 'normalization' (default, fast)
                       or 'riccati' (coherence-function route shared with the
                       vortex-lattice trajectory solver)
    @param sw_find_tc: if True, bisect for Tc at the given impurity setting
    @param sw_imp_sweep: if True, sweep imp_sweep (Gamma values) and write Tc(Gamma)
    @param  imp_sweep: array of Gamma values [eV] for the sweep
    """
    print("calculate homogeneous quasiclassical Eilenberger equation", flush=True)
    Nk, klist, eig, uni, kweight = get_emesh(Nx, Ny, Nz, ham_r, S_r, rvec, avec, sw_uni=True)
    wf, phif = build_fs(eig, klist, mu, gap_sym, fs_width)
    print(f"Fermi-surface points kept: {len(wf)} (of {eig.size})", flush=True)
    if len(wf) == 0:
        print("Error: no Fermi-surface points found; check mu / fs_width", flush=True)
        return
    omega = matsubara(temp, wc)
    print(f"Matsubara cutoff wc = {wc:.4e} eV ({len(omega)} freqs at this T)", flush=True)
    print(f"pairing coupling lambda = {coupling:.4f}, gap_sym = {gap_sym}", flush=True)
    print(f"impurity: Gamma = {imp_gamma:.4e} eV, c = {imp_c:.3e} "
          f"({'clean' if imp_gamma == 0 else 'Born' if imp_c > 10 else 'unitary'})", flush=True)
    print(f"homogeneous solver method: {method}", flush=True)

    lam = pairing_eigenvalue(temp, wf, phif, omega, imp_gamma, imp_c, coupling, method)
    print(f"lambda_lin(T={temp/kb:.2f} K) = {lam:.6f} "
          f"({'SC' if lam > 1 else 'normal'})", flush=True)

    damp = solve_gap(temp, wf, phif, omega, imp_gamma, imp_c, coupling, method)
    print(f"gap amplitude Delta = {damp:.6e} eV "
          f"(Delta/kBT = {damp/temp:.4f})", flush=True)
    # frequency-resolved Fermi-surface-averaged anomalous propagator at the solved gap
    _, f = solve_quasiclassical(omega, phif * damp, wf, imp_gamma, imp_c, method)
    favg = _fs_average(f, wf)
    try:
        with open('eilenberger_f.dat', 'w') as fh:
            fh.write("# omega_n[eV]  Re<f>  Im<f>\n")
            for wl, fa in zip(omega, favg):
                fh.write(f"{wl:14.6e} {fa.real:14.6e} {fa.imag:14.6e}\n")
    except IOError as e:
        print(f"Error: failed to write 'eilenberger_f.dat': {e}", flush=True)

    if sw_find_tc:
        tc = find_tc(wf, phif, imp_gamma, imp_c, coupling, wc, 1.0e-4, 5.0 * temp, method)
        if tc > 0:
            print(f"Tc = {tc:.6e} eV ({tc/kb:.3f} K)", flush=True)
        else:
            print("no superconducting Tc found in the scanned range", flush=True)

    if sw_imp_sweep and imp_sweep is not None:
        print("sweep impurity scattering: Gamma[eV]  Tc[eV]  Tc[K]", flush=True)
        try:
            with open('eilenberger_tc.dat', 'w') as fh:
                fh.write("# Gamma[eV]  Tc[eV]  Tc[K]\n")
                for g_imp in imp_sweep:
                    tc = find_tc(wf, phif, float(g_imp), imp_c, coupling, wc, 1.0e-4, 5.0 * temp, method)
                    fh.write(f"{g_imp:14.6e} {tc:14.6e} {tc/kb:10.4f}\n")
                    print(f"  {g_imp:12.4e} {tc:12.4e} {tc/kb:10.3f}", flush=True)
        except IOError as e:
            print(f"Error: failed to write 'eilenberger_tc.dat': {e}", flush=True)
    return damp


# --------------------------------------------------------------------------- #
#  Superfluid density / penetration depth (quasiclassical electrodynamic response)
# --------------------------------------------------------------------------- #
# The London penetration depth that screens the field (the Maxwell constitutive
# relation, 1/lambda^2 ~ rho_s) follows from the quasiclassical superfluid density
#   rho_s(T)/rho_s(0) = pi <v_x^2 * 2T sum_{n>=0} Delta_k^2/(w_n^2+Delta_k^2)^{3/2}>_FS / <v_x^2>_FS
# (the bracket -> 1 as T->0).  Full-gap (s-wave) gives an exponentially flat
# low-T rho_s; line nodes (d-wave) give the linear-in-T rho_s (and lambda(T))
# that is the hallmark of nodal superconductors.

def superfluid_density(coupling: float, temp: float, wc: float, gap_sym: str = 's',
                       Nbeta: int = 360):
    """
    @fn superfluid_density
    @brief Self-consistent bulk gap Delta(T) and normalized superfluid density
    rho_s(T)/rho_s(0) on a model cylindrical Fermi surface (v_F || k).
    @param gap_sym: 's','d','dxy','px','py','p+ip'/'p-ip'
    @return (Delta, rho_ratio): bulk gap [eV] and rho_s(T)/rho_s(0) (0 if normal)
    """
    from ._eilenberger_surface import _bulk_gap
    from ._eilenberger_vortex import _ff_vortex
    omega = matsubara(temp, wc)
    beta = np.linspace(0.0, 2.0 * np.pi, Nbeta, endpoint=False)
    phi = _ff_vortex(beta, gap_sym)
    vx2 = np.cos(beta) ** 2
    Delta = _bulk_gap(coupling, temp, omega, phi)
    if Delta <= 1.0e-6 * temp:
        return 0.0, 0.0
    Dk = Delta * np.abs(phi)                                  # gap magnitude per direction
    fac = 2.0 * temp * (Dk[:, None] ** 2 /
                        (omega[None, :] ** 2 + Dk[:, None] ** 2) ** 1.5).sum(axis=1)
    rho_ratio = np.pi * (vx2 * fac).sum() / vx2.sum()        # -> 1 as T->0
    return Delta, rho_ratio


def calc_penetration_depth(coupling: float, temp: float, wc: float, gap_sym: str = 's',
                           t_list=None, kb: float = 1.0):
    """
    @fn calc_penetration_depth
    @brief Temperature sweep of the superfluid density rho_s(T)/rho_s(0) and the
    penetration depth lambda(T)/lambda(0)=1/sqrt(rho_s(T)/rho_s(0)).  Distinguishes
    full-gap (s-wave: exponentially flat low-T) from nodal (d-wave: linear-in-T)
    superconductors -- the classic penetration-depth fingerprint.  Writes
    'penetration_depth.dat'.
    @param t_list: temperatures [eV] (default a sweep from ~0.05 Tc to Tc)
    """
    # estimate Tc from the linearized eigenvalue scale: bracket via Delta(T)
    if t_list is None:
        # geometric-ish sweep up to where the gap vanishes
        t_hi = temp
        # ensure t_hi below Tc (gap nonzero); else scale down
        for _ in range(8):
            d, _r = superfluid_density(coupling, t_hi, wc, gap_sym)
            if d > 0:
                break
            t_hi *= 0.6
        t_list = np.linspace(0.05 * t_hi, 1.05 * t_hi, 16)
    print(f"penetration depth / superfluid density: {gap_sym}, lambda={coupling:.3f}", flush=True)
    print("  T[K]    Delta/D0   rho_s(T)/rho_s(0)   lambda(T)/lambda(0)", flush=True)
    rows = []
    D0 = None
    for T in t_list:
        D, rho = superfluid_density(coupling, T, wc, gap_sym)
        if D0 is None and D > 0:
            D0 = D
        lam = 1.0 / np.sqrt(rho) if rho > 1e-9 else np.inf
        rows.append((T, D, rho, lam))
        print(f"  {T/kb:6.2f}  {D/(D0 if D0 else 1):8.4f}   {rho:12.5f}     {lam:10.4f}", flush=True)
    try:
        with open('penetration_depth.dat', 'w') as fh:
            fh.write("# T[eV]   Delta[eV]   rho_s(T)/rho_s(0)   lambda(T)/lambda(0)\n")
            for T, D, rho, lam in rows:
                fh.write(f"{T:12.6e} {D:12.6e} {rho:12.6e} {lam:12.6e}\n")
    except IOError as e:
        print(f"Error writing penetration_depth.dat: {e}", flush=True)
    return rows


# =========================================================================== #
#  Model Fermi surface with Fermi velocities
# --------------------------------------------------------------------------- #
# Generalizes the isotropic v_F || k, |v_F|=const cylinder used by the bulk and
# trajectory solvers.  A 2D model dispersion eps(kx,ky) defines a (convex,
# single-sheet) Fermi surface; for each FS point we store k_F, the Fermi velocity
# v_F = grad eps (generally NOT parallel to k), |v_F|, the unit direction v_hat,
# and the DOS weight nf = dl/|v_F| (arc length over velocity).  This is the data
# the trajectory solvers need: the trajectory runs along v_hat, the chord step
# scales with the perpendicular velocity, and FS averages are nf-weighted (the je
# convention Dkx/Dky, Dvx/Dvy, Dnf).
#
#   'iso'     : eps = (kx^2+ky^2)/2 - mu          -> circular FS, v_F || k, |v_F| const
#   'ellipse' : eps = kx^2/(2 mx)+ky^2/(2 my)-mu  -> anisotropic mass (v_F not || k)
#   'tb'      : eps = -2 t (cos kx + cos ky) - mu -> tight-binding (below van Hove)


def _disp(kind, params):
    """Return (eps(kx,ky), grad eps -> (vx,vy), default mu, radial search rmax)."""
    if kind == 'iso':
        return (lambda kx, ky: 0.5 * (kx ** 2 + ky ** 2),
                lambda kx, ky: (kx, ky), 1.0, 10.0)
    if kind == 'ellipse':
        mx, my = (params or (1.0, 0.5))
        return (lambda kx, ky: 0.5 * (kx ** 2 / mx + ky ** 2 / my),
                lambda kx, ky: (kx / mx, ky / my), 1.0, 10.0)
    if kind == 'tb':
        t = (params or 1.0)
        t = t if np.isscalar(t) else t[0]
        return (lambda kx, ky: -2.0 * t * (np.cos(kx) + np.cos(ky)),
                lambda kx, ky: (2.0 * t * np.sin(kx), 2.0 * t * np.sin(ky)), -2.0 * t * 0.6, np.pi - 1e-6)
    raise ValueError(f"unknown FS kind: {kind}")


def build_model_fs(kind: str = 'iso', Nth: int = 360, mu: float = None, params=None):
    """
    @fn build_model_fs
    @brief Build a model Fermi surface (radial parametrization, convex single sheet)
    with Fermi velocities and DOS weights.
    @param kind: 'iso', 'ellipse' (params=(mx,my)), or 'tb' (params=t)
    @param Nth: number of FS points (angular samples)
    @param  mu: chemical potential (default per dispersion)
    @return dict: th,kx,ky,vx,vy,vabs,vhx,vhy,nf  (nf normalized to sum 1)
    """
    eps, grad, mu0, rmax = _disp(kind, params)
    if mu is None:
        mu = mu0
    th = np.linspace(0.0, 2.0 * np.pi, Nth, endpoint=False)
    kF = np.empty(Nth)
    for i, t in enumerate(th):
        c, s = np.cos(t), np.sin(t)
        f = lambda r: eps(r * c, r * s) - mu
        lo, hi = 1e-6, rmax
        try:
            kF[i] = brentq(f, lo, hi)
        except ValueError:
            rs = np.linspace(lo, hi, 200)
            fv = np.array([f(r) for r in rs])
            sgn = np.where(np.diff(np.sign(fv)) != 0)[0]
            kF[i] = brentq(f, rs[sgn[0]], rs[sgn[0] + 1]) if len(sgn) else np.nan
    kx, ky = kF * np.cos(th), kF * np.sin(th)
    vx, vy = grad(kx, ky)
    vabs = np.sqrt(vx ** 2 + vy ** 2)
    dth = th[1] - th[0]
    dkF = (np.roll(kF, -1) - np.roll(kF, 1)) / (2.0 * dth)
    dl = np.sqrt(kF ** 2 + dkF ** 2) * dth                # arc length element
    nf = dl / vabs
    nf = nf / nf.sum()
    return dict(th=th, kx=kx, ky=ky, vx=vx, vy=vy, vabs=vabs,
                vhx=vx / vabs, vhy=vy / vabs, nf=nf)


def fs_form_factor(fs: dict, gap_sym: str) -> np.ndarray:
    """Pairing form factor phi at each FS point (function of the k-direction),
    normalized to the nf-weighted <|phi|^2> = 1.  Matches the gap symmetries used
    elsewhere (s, d=d_{x^2-y^2}, dxy, px, py, chiral p+ip/p-ip)."""
    a = np.arctan2(fs['ky'], fs['kx'])
    if gap_sym == 's':
        phi = np.ones_like(a, dtype=np.complex128)
    elif gap_sym == 'd':
        phi = np.cos(2.0 * a).astype(np.complex128)
    elif gap_sym == 'dxy':
        phi = np.sin(2.0 * a).astype(np.complex128)
    elif gap_sym == 'px':
        phi = np.cos(a).astype(np.complex128)
    elif gap_sym == 'py':
        phi = np.sin(a).astype(np.complex128)
    elif gap_sym in ('p+ip', 'chiral'):
        phi = np.exp(1j * a)
    elif gap_sym == 'p-ip':
        phi = np.exp(-1j * a)
    else:
        raise ValueError(f"unknown gap_sym: {gap_sym}")
    norm = np.sqrt((fs['nf'] * np.abs(phi) ** 2).sum())
    return phi / norm if norm > 0 else phi


def bulk_gap_fs(coupling, temp, omega, fs, gap_sym, eps=1e-8, itemax=500, mix=0.5):
    """Homogeneous gap amplitude on a model FS (nf-weighted FS average)."""
    phi = fs_form_factor(fs, gap_sym)
    nf = fs['nf']
    a2 = (phi * np.conj(phi)).real
    damp = 1.764 * temp
    for _ in range(itemax):
        R = np.sqrt(omega[None, :] ** 2 + a2[:, None] * damp ** 2)
        f = (phi[:, None] * damp) / R
        phif = (nf[:, None] * np.conj(phi[:, None]) * f).sum(axis=0)
        damp_new = coupling * 2.0 * temp * phif.sum().real
        if damp_new <= 0.0:
            return 0.0
        if abs(damp_new - damp) < eps * max(damp_new, eps):
            return damp_new
        damp = (1.0 - mix) * damp + mix * damp_new
    return damp


def superfluid_density_fs(coupling, temp, wc, fs, gap_sym):
    """
    @fn superfluid_density_fs
    @brief Superfluid-density tensor (rho_xx, rho_yy)/rho_n(0) on a model FS with
    Fermi velocities.  rho_ii(T)/rho_ii(0) = pi <v_i^2 2T sum_n Dk^2/(w^2+Dk^2)^{3/2}>_nf
    / <v_i^2>_nf (the bracket -> 1 as T->0).  Anisotropic FS / nodes give an
    anisotropic penetration depth lambda_ii = 1/sqrt(rho_ii).
    @return (Delta, rho_xx, rho_yy)
    """
    omega = matsubara(temp, wc)
    Delta = bulk_gap_fs(coupling, temp, omega, fs, gap_sym)
    if Delta <= 1.0e-6 * temp:
        return 0.0, 0.0, 0.0
    phi = fs_form_factor(fs, gap_sym)
    Dk = Delta * np.abs(phi)
    fac = 2.0 * temp * (Dk[:, None] ** 2 /
                        (omega[None, :] ** 2 + Dk[:, None] ** 2) ** 1.5).sum(axis=1)
    nf = fs['nf']
    vx2, vy2 = fs['vhx'] ** 2, fs['vhy'] ** 2
    rho_xx = np.pi * (nf * vx2 * fac).sum() / (nf * vx2).sum()
    rho_yy = np.pi * (nf * vy2 * fac).sum() / (nf * vy2).sum()
    return Delta, rho_xx, rho_yy


def calc_fs_penetration(coupling, temp, wc, kind='ellipse', gap_sym='s', Nth=360,
                        mu=None, params=None, t_list=None, kb=1.0):
    """
    @fn calc_fs_penetration
    @brief Temperature sweep of the anisotropic superfluid density / penetration
    depth on a model Fermi surface with Fermi velocities.  An anisotropic FS gives
    lambda_xx != lambda_yy, and nodes give the linear-in-T penetration depth.
    Writes 'fs_penetration.dat'.
    """
    fs = build_model_fs(kind, Nth, mu, params)
    vx2 = (fs['nf'] * fs['vhx'] ** 2).sum()
    vy2 = (fs['nf'] * fs['vhy'] ** 2).sum()
    print(f"model FS '{kind}' (Nth={Nth}): <|v_F|>={fs['vabs'].mean():.4f}, "
          f"<v_x^2>/<v_y^2> = {vx2:.3f}/{vy2:.3f}, "
          f"absolute lambda_xx/lambda_yy(0) = {np.sqrt(vy2/vx2):.3f}", flush=True)
    print(f"penetration depth on model FS: {gap_sym}-wave, lambda={coupling:.3f}", flush=True)
    if t_list is None:
        t_hi = temp
        for _ in range(8):
            D, _x, _y = superfluid_density_fs(coupling, t_hi, wc, fs, gap_sym)
            if D > 0:
                break
            t_hi *= 0.6
        t_list = np.linspace(0.05 * t_hi, 1.02 * t_hi, 14)
    print("  T[K]    Delta     rho_xx   rho_yy   lam_xx/lam_yy", flush=True)
    rows = []
    for T in t_list:
        D, rxx, ryy = superfluid_density_fs(coupling, T, wc, fs, gap_sym)
        ratio = np.sqrt(ryy / rxx) if (rxx > 1e-9 and ryy > 1e-9) else np.nan
        rows.append((T, D, rxx, ryy, ratio))
        print(f"  {T/kb:6.2f}  {D:9.3e}  {rxx:7.4f}  {ryy:7.4f}   {ratio:7.4f}", flush=True)
    try:
        with open('fs_penetration.dat', 'w') as fh:
            fh.write("# T[eV]  Delta[eV]  rho_xx  rho_yy  lambda_xx/lambda_yy\n")
            for T, D, rxx, ryy, ratio in rows:
                fh.write(f"{T:12.6e} {D:12.6e} {rxx:12.6e} {ryy:12.6e} {ratio:12.6e}\n")
    except IOError as e:
        print(f"Error writing fs_penetration.dat: {e}", flush=True)
    return rows
