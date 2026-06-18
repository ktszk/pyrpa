#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Spin-resolved (2x2) homogeneous quasiclassical Eilenberger solver: singlet and
triplet (d-vector) pairing with a Zeeman (Maki) field.

Where the scalar solver treats one pseudo-spin sector, this works in the full
Nambu (x) spin 4x4 space.  For a homogeneous superconductor the Eilenberger
equation reduces to the algebraic normalization

    g_hat = N_hat / sqrt(N_hat^2),   N_hat = [[ w_tilde , Delta_hat ],
                                              [ Delta_hat^dag , -w_tilde ]],

with the renormalized frequency carrying the Zeeman shift along z,
    w_tilde = w_n * I + i h sigma_z,   h = mu * B (Maki),
and the 2x2 spin gap matrix
    singlet : Delta_hat = i sigma_y                       (Delta = phi(k) Damp)
    triplet : Delta_hat = (d_hat . sigma) i sigma_y       (d-vector along d_hat).

The matrix square root is taken by eigendecomposition of N_hat^2 (batched over
Fermi-surface points and Matsubara frequencies).  The normal (g) and anomalous
(f) 2x2 blocks of g_hat give the DOS and the gap equation (projected onto the
pairing channel).  This reproduces singlet Pauli limiting and, for triplet,
the Pauli immunity when the d-vector is perpendicular to the field (equal-spin
pairing) versus depairing when d || h.
"""
import numpy as np
from ._eilenberger import matsubara, build_fs


def _pauli():
    sx = np.array([[0, 1], [1, 0]], dtype=np.complex128)
    sy = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
    sz = np.array([[1, 0], [0, -1]], dtype=np.complex128)
    return sx, sy, sz


def gap_matrix(channel: str, dvec=(0.0, 0.0, 1.0)) -> np.ndarray:
    """2x2 spin gap matrix Delta_hat (unit amplitude).
    @param channel: 'singlet' or 'triplet'
    @param   dvec: triplet d-vector direction (ignored for singlet)
    """
    sx, sy, sz = _pauli()
    isy = 1j * sy                                   # [[0,1],[-1,0]]
    if channel == 'singlet':
        return isy
    d = np.asarray(dvec, dtype=np.complex128)
    d = d / np.sqrt((np.abs(d) ** 2).sum())
    return (d[0] * sx + d[1] * sy + d[2] * sz) @ isy


def _ghat_homogeneous(omega: np.ndarray, Dk: np.ndarray, Shat: np.ndarray, h: float):
    """Nambu(x)spin g_hat = N/sqrt(N^2) for all (FS point, frequency).
    @param omega: Matsubara frequencies [Nw] (use the full +/- set)
    @param    Dk: complex gap amplitude phi(k)*Damp per FS point [Nfs]
    @param  Shat: 2x2 spin gap matrix (unit amplitude)
    @param     h: Zeeman energy along z
    @return (g, f): 2x2 normal / anomalous blocks [Nfs, Nw, 2, 2]
    """
    _, _, sz = _pauli()
    Nfs, Nw = len(Dk), len(omega)
    I2 = np.eye(2, dtype=np.complex128)
    wn = omega[None, :, None, None] * I2                        # frequency (tau3 in Nambu)
    zee = 1j * h * sz                                           # Zeeman (tau0: SAME sign in both blocks)
    upper = np.broadcast_to(wn + zee, (Nfs, Nw, 2, 2))
    lower = np.broadcast_to(-wn + zee, (Nfs, Nw, 2, 2))
    D = Dk[:, None, None, None] * Shat                          # [Nfs,1,2,2] gap matrix
    D = np.broadcast_to(D, (Nfs, Nw, 2, 2))
    Dd = np.conj(np.swapaxes(D, -1, -2))                        # Delta_hat^dag
    N = np.empty((Nfs, Nw, 4, 4), dtype=np.complex128)
    N[..., :2, :2] = upper
    N[..., :2, 2:] = D
    N[..., 2:, :2] = Dd
    N[..., 2:, 2:] = lower
    N2 = N @ N
    lam, V = np.linalg.eig(N2)                                  # batched [..,4],[..,4,4]
    inv_sqrt = 1.0 / np.sqrt(lam)                               # principal branch (Re>0)
    Vinv = np.linalg.inv(V)
    sq = V @ (inv_sqrt[..., None] * Vinv)                       # (N^2)^{-1/2}
    ghat = N @ sq
    return ghat[..., :2, :2], ghat[..., :2, 2:]


def solve_gap_spin(coupling, temp, omega, wf, phif, channel='singlet',
                   dvec=(0, 0, 1), h=0.0, damp_init=None, eps=1e-8, itemax=400, mix=0.5):
    """
    @fn solve_gap_spin
    @brief Self-consistent homogeneous gap amplitude for a singlet or triplet
    (d-vector) channel in a Zeeman field h (Maki), via the 4x4 Nambu(x)spin g_hat.
    Gap equation: Damp = lambda * T sum_n < phi* tr(Shat^dag f)/tr(Shat^dag Shat) >.
    @return Damp [eV]
    """
    Shat = gap_matrix(channel, dvec)
    trSS = np.trace(np.conj(Shat).T @ Shat).real               # tr(S^dag S)
    Sd = np.conj(Shat).T
    om = np.concatenate([omega, -omega])                       # full +/- Matsubara set
    if damp_init is None:
        damp_init = 1.764 * temp
    damp = damp_init
    Wn = wf.sum()
    for _ in range(itemax):
        Dk = phif * damp                                       # [Nfs] (real form factor)
        _, f = _ghat_homogeneous(om, Dk.astype(np.complex128), Shat, h)   # [Nfs,Nw,2,2]
        proj = np.einsum('ij,kwji->kw', Sd, f) / trSS          # tr(S^dag f) [Nfs,Nw]
        # gap channel amplitude: FS average of phi* proj, Matsubara sum
        amp = np.tensordot(wf, np.conj(phif)[:, None] * proj, axes=(0, 0)).sum() / Wn
        damp_new = coupling * temp * amp.real
        if damp_new <= 0.0:
            return 0.0
        if abs(damp_new - damp) < eps * max(damp_new, eps):
            return damp_new
        damp = (1.0 - mix) * damp + mix * damp_new
    return damp


def calc_spin_pauli(Nx, Ny, Nz, wc, ham_r, S_r, rvec, avec, mu, temp, coupling,
                    gap_sym=0, h_list=None, fs_width=5.0e-3, kb=1.0):
    """
    @fn calc_spin_pauli
    @brief Compare the Zeeman (Maki) field response of singlet vs triplet pairing
    using the 2x2 spin solver: singlet and triplet d||h are Pauli-limited, triplet
    d perpendicular to h is Pauli-immune (equal-spin pairing).  Sweeps h and writes
    'spin_pauli.dat' (Delta(h)/Delta0 for each channel).
    """
    from ._bands import get_emesh
    omega = matsubara(temp, wc)
    Nk, klist, eig, uni, kweight = get_emesh(Nx, Ny, Nz, ham_r, S_r, rvec, avec, sw_uni=True)
    wf, phif = build_fs(eig, klist, mu, gap_sym, fs_width)
    cases = [('singlet', (0, 0, 1), 'singlet'),
             ('triplet', (0, 0, 1), 'triplet d||h (z)'),
             ('triplet', (1, 0, 0), 'triplet d_|_h (x)')]
    D0 = {}
    for ch, dv, name in cases:
        D0[name] = solve_gap_spin(coupling, temp, omega, wf, phif, ch, dv, h=0.0, damp_init=2e-3)
    print(f"spin/Zeeman (Maki) response, T={temp/kb:.2f} K:", flush=True)
    for name in D0:
        print(f"  Delta0[{name}] = {D0[name]:.4e} eV", flush=True)
    Dref = max(D0.values())
    if Dref <= 0:
        print("normal state; nothing to sweep", flush=True)
        return
    if h_list is None:
        h_list = Dref * np.array([0.0, 0.3, 0.5, 0.7, 0.9, 1.1, 1.4, 1.8])
    try:
        with open('spin_pauli.dat', 'w') as fh:
            fh.write("# h/Delta0  " + "  ".join(n.replace(' ', '_') for _, _, n in cases) + "\n")
            print("  h/D0   " + "   ".join(f"{n[:12]:>12}" for _, _, n in cases), flush=True)
            for h in h_list:
                row = []
                for ch, dv, name in cases:
                    Dh = solve_gap_spin(coupling, temp, omega, wf, phif, ch, dv,
                                        h=float(h), damp_init=D0[name] if D0[name] > 0 else 2e-3)
                    row.append(Dh / Dref)
                fh.write(f"{h/Dref:8.4f} " + " ".join(f"{r:10.4f}" for r in row) + "\n")
                print(f"  {h/Dref:5.2f}  " + "   ".join(f"{r:12.4f}" for r in row), flush=True)
    except IOError as e:
        print(f"Error writing spin_pauli.dat: {e}", flush=True)
    return D0
