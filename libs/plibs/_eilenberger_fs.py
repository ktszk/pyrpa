#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Model Fermi surface with Fermi velocities for the quasiclassical Eilenberger
solvers (generalizes the isotropic v_F || k, |v_F|=const cylinder).

A 2D model dispersion eps(kx,ky) defines a (convex, single-sheet) Fermi surface;
for each FS point we store the position k_F, the Fermi velocity v_F = grad eps
(generally NOT parallel to k), its magnitude |v_F|, the unit direction v_hat, and
the density-of-states weight nf = dl/|v_F| (arc length over velocity).  This is
the data the trajectory solvers need: the trajectory runs along v_hat, the chord
step scales with the perpendicular velocity, and Fermi-surface averages are
nf-weighted (the je convention: Dkx/Dky, Dvx/Dvy, Dnf).

Built-in dispersions:
  'iso'     : eps = (kx^2+ky^2)/2 - mu          -> circular FS, v_F || k, |v_F| const
  'ellipse' : eps = kx^2/(2 mx)+ky^2/(2 my)-mu  -> anisotropic mass (v_F not || k)
  'tb'      : eps = -2 t (cos kx + cos ky) - mu -> tight-binding (below van Hove)
The real tight-binding band can be plugged in later by supplying its own contour.
"""
import numpy as np
from scipy.optimize import brentq


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
        # bracket
        try:
            kF[i] = brentq(f, lo, hi)
        except ValueError:
            # scan for a sign change
            rs = np.linspace(lo, hi, 200)
            fv = np.array([f(r) for r in rs])
            sgn = np.where(np.diff(np.sign(fv)) != 0)[0]
            kF[i] = brentq(f, rs[sgn[0]], rs[sgn[0] + 1]) if len(sgn) else np.nan
    kx, ky = kF * np.cos(th), kF * np.sin(th)
    vx, vy = grad(kx, ky)
    vabs = np.sqrt(vx ** 2 + vy ** 2)
    # arc length dl = sqrt(kF^2 + (dkF/dth)^2) dth (periodic derivative)
    dth = th[1] - th[0]
    dkF = (np.roll(kF, -1) - np.roll(kF, 1)) / (2.0 * dth)
    dl = np.sqrt(kF ** 2 + dkF ** 2) * dth
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
    from ._eilenberger import matsubara
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
    depth on a model Fermi surface with Fermi velocities.  Demonstrates the model
    FS + v_F capability: an anisotropic FS gives lambda_xx != lambda_yy, and nodes
    give the linear-in-T penetration depth.  Writes 'fs_penetration.dat'.
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
