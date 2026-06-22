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
from scipy.linalg import expm
from ._eilenberger import matsubara, build_fs

try:                                       # optional Fortran acceleration of the 2x2 matrix Riccati
    from ..flibs import matrix_riccati_batch as _fort_mat
    _HAVE_FORT = True
except Exception:
    _HAVE_FORT = False


def _mat_batch(om, Dpath, hvf, ds, h):
    """g, f [Ns, Nw, 2, 2] along one trajectory for all frequencies: Fortran if
    available, else a Python loop over matrix_trajectory_gf."""
    if _HAVE_FORT:
        return _fort_mat(np.ascontiguousarray(om, dtype=np.complex128),
                         np.ascontiguousarray(Dpath, dtype=np.complex128), hvf, ds, h)
    Ns, Nw = len(Dpath), len(om)
    g = np.empty((Ns, Nw, 2, 2), dtype=np.complex128)
    f = np.empty((Ns, Nw, 2, 2), dtype=np.complex128)
    for iw, wn in enumerate(om):
        g[:, iw], f[:, iw] = matrix_trajectory_gf(complex(wn), Dpath, hvf, ds, h)
    return g, f


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


# --------------------------------------------------------------------------- #
#  Spin-matrix Riccati for INHOMOGENEOUS systems (surface / vortex / junction)
# --------------------------------------------------------------------------- #
# The 2x2 spin coherence matrix a(s) obeys, along a quasiclassical trajectory,
#     hbar v_F da/ds = Delta_hat - 2 w_n a - a Delta_hat^dag a + i h [sigma_z, a]
# (reduces to the scalar Riccati for a singlet without Zeeman).  The conjugate
# matrix b uses source Delta_hat^dag and -i h [sigma_z, b].  Each step is taken
# with the *exact* fractional-linear (Mobius) propagator of the associated linear
# spinor system M = exp(N ds/hbar v_F) -- the matrix generalization of the scalar
# tanh step -- which is unconditionally stable for the stiff high-frequency modes.

def _riccati_N(omega: complex, D: np.ndarray, h: float, is_a: bool) -> np.ndarray:
    """4x4 generator N for the spinor system whose ratio gives the matrix Riccati
    solution (a = u w^{-1}).  is_a: source D=Delta_hat; else source D=Delta_hat^dag."""
    _, _, sz = _pauli()
    src = D
    quad = np.conj(D).T                       # a-eq: Delta^dag (so -a Delta^dag a); b-eq: Delta
    sgn = 1.0 if is_a else -1.0
    N = np.empty((4, 4), dtype=np.complex128)
    N[:2, :2] = -2.0 * omega * np.eye(2) + 1j * sgn * h * sz   # P
    N[:2, 2:] = src                                            # Q
    N[2:, :2] = quad                                           # R
    N[2:, 2:] = 1j * sgn * h * sz                              # S
    return N


def _mobius_step(a: np.ndarray, N: np.ndarray, t: float) -> np.ndarray:
    """One exact fractional-linear step: a -> (M11 a + M12)(M21 a + M22)^{-1},
    M = exp(N t).  Unconditionally stable (contracts to the stable Riccati root)."""
    M = expm(N * t)
    u = M[:2, :2] @ a + M[:2, 2:]
    w = M[2:, :2] @ a + M[2:, 2:]
    return u @ np.linalg.inv(w)


def riccati_matrix_bulk(omega: complex, D: np.ndarray, is_a: bool) -> np.ndarray:
    """Homogeneous (bulk) matrix Riccati root for unitary gaps: a = D/(omega+E),
    E = sqrt(omega^2 + |D|^2) with |D|^2 = (D D^dag) assumed scalar*I (unitary)."""
    DDd = D @ np.conj(D).T
    d2 = 0.5 * np.trace(DDd).real            # |Delta|^2 (unitary: DDd = d2 * I)
    E = np.sqrt(omega ** 2 + d2)
    src = D if is_a else np.conj(D).T
    return src / (omega + E)


def integrate_riccati_matrix(omega: complex, Dpath, hvf: float, ds: float,
                             a0: np.ndarray, h: float, is_a: bool) -> np.ndarray:
    """Integrate the 2x2 matrix Riccati along a trajectory with the stable Mobius
    step.  Dpath is a list/array of 2x2 gap matrices along the path.
    @return a along the trajectory [Ns, 2, 2]
    """
    Ns = len(Dpath)
    out = np.empty((Ns, 2, 2), dtype=np.complex128)
    a = np.array(a0, dtype=np.complex128)
    out[0] = a
    t = ds / hvf
    for i in range(Ns - 1):
        Dmid = 0.5 * (Dpath[i] + Dpath[i + 1])
        a = _mobius_step(a, _riccati_N(omega, Dmid, h, is_a), t)
        out[i + 1] = a
    return out


def matrix_gf(a: np.ndarray, b: np.ndarray):
    """Quasiclassical 2x2 g, f from the coherence matrices a, b (je convention):
    g = (I + a b)^{-1}(I - a b),  f = (I + a b)^{-1} 2 a."""
    I = np.eye(2, dtype=np.complex128)
    P = np.linalg.inv(I + a @ b)
    return P @ (I - a @ b), P @ (2.0 * a)


def matrix_trajectory_gf(omega: complex, Dpath, hvf: float, ds: float, h: float = 0.0):
    """
    @fn matrix_trajectory_gf
    @brief Spin-matrix quasiclassical g, f along one inhomogeneous 1D trajectory
    (surface / junction / domain wall): a is integrated forward from the upstream
    bulk root, b backward (reversed path) from the downstream bulk root, then
    g=(I+ab)^{-1}(I-ab), f=(I+ab)^{-1}2a at every point.  Dpath = list of 2x2 gap
    matrices Delta_hat(x) along the path.
    @return (g, f): [Ns, 2, 2] each
    """
    Ns = len(Dpath)
    a0 = riccati_matrix_bulk(omega, Dpath[0], True)        # upstream bulk (h=0 seed; relaxes)
    a = integrate_riccati_matrix(omega, Dpath, hvf, ds, a0, h, True)
    b0 = riccati_matrix_bulk(omega, Dpath[-1], False)      # downstream bulk
    b_rev = integrate_riccati_matrix(omega, list(Dpath[::-1]), hvf, ds, b0, h, False)
    b = b_rev[::-1]
    g = np.empty((Ns, 2, 2), dtype=np.complex128)
    f = np.empty((Ns, 2, 2), dtype=np.complex128)
    for i in range(Ns):
        g[i], f[i] = matrix_gf(a[i], b[i])
    return g, f


# --------------------------------------------------------------------------- #
#  Self-consistent d-vector TEXTURE at a specular surface (triplet)
# --------------------------------------------------------------------------- #
def _default_dvector_channels():
    """Triplet d = cos(b) e_x (orbital p_x, sign-changing at an x-surface) +
    sin(b) e_z (orbital p_y, even).  Spin matrices Shat_x = sigma_x i sigma_y =
    diag(-1,1) and Shat_z = sigma_z i sigma_y = sigma_x are both traceless and give
    a unitary gap; the e_y component (Shat_y proportional to identity) is avoided
    because the identity piece is not handled by the unitary matrix Riccati bulk
    seed.  Each orbital factor is normalized to <|phi|^2>=1 over the circle."""
    sx, sy, sz = _pauli()
    isy = 1j * sy
    return [('px(e_x)', (lambda b: np.sqrt(2.0) * np.cos(b)), sx @ isy),
            ('py(e_z)', (lambda b: np.sqrt(2.0) * np.sin(b)), sz @ isy)]


def _bulk_dvector(couplings, temp, om, beta, w, phitil, eps=1e-10, itemax=2000, mix=0.5):
    """Coupled bulk amplitudes D_a of a multi-component unitary triplet with per-
    channel couplings.  R(beta) = sqrt(w^2 + sum_b (D_b phi_b)^2) couples the
    channels; gap_a: D_a = lambda_a T sum_n < phi_a^2 D_a / R > over the full
    direction set.  Returns the vector (D_a).  Subcritical channels relax to 0."""
    nc = len(phitil)
    lam = np.asarray(couplings, dtype=float)
    phi = np.array([phitil[a](beta) for a in range(nc)])               # [nc, Nb] outgoing
    phii = np.array([phitil[a](np.pi - beta) for a in range(nc)])      # reflected branch
    D = np.full(nc, 1.764 * temp)
    for _ in range(itemax):
        s2 = ((D[:, None] * phi) ** 2).sum(axis=0)                     # sum_b (D_b phi_b)^2 [Nb]
        s2i = ((D[:, None] * phii) ** 2).sum(axis=0)
        Ro = np.sqrt(om[None, :] ** 2 + s2[:, None])                   # [Nb, Nw]
        Ri = np.sqrt(om[None, :] ** 2 + s2i[:, None])
        Dnew = np.empty(nc)
        for a in range(nc):
            amp = (w * (phi[a][:, None] ** 2 / Ro).sum(axis=1)).sum() \
                + (w * (phii[a][:, None] ** 2 / Ri).sum(axis=1)).sum()
            Dnew[a] = max(lam[a] * temp * D[a] * amp, 0.0)
        if np.abs(Dnew - D).max() < eps * max(Dnew.max(), eps):
            return Dnew
        D = (1.0 - mix) * D + mix * Dnew
    return D


def solve_surface_dvector(couplings, temp, omega, channels=None, Dbulk=None,
                          Lxi=8.0, nper=8, Nbeta=16, hvf=1.0, cos_min=0.06,
                          eps=3e-3, itemax=60, mix=0.4):
    """
    @fn solve_surface_dvector
    @brief Self-consistent d-vector TEXTURE at a specular surface (x>=0) for a
    multi-component unitary triplet, via the 2x2 spin-matrix Riccati.  Each spin
    component a has an orbital form factor phi_a(beta), a spin matrix Shat_a, and its
    own coupling lambda_a; the (complex) spatial amplitude Delta_a(x) is solved
    self-consistently.  Components whose orbital factor is sign-changing under
    specular reflection (beta -> pi-beta) are suppressed at the surface, so a
    dominant + subdominant pair makes the net d-vector reorient in spin space as the
    surface is approached -- the texture.
    @param couplings: per-channel lambda_a (use a dominant + subdominant pair).
    @return (x, Damp [Ncomp, Ng] complex, Dbulk [Ncomp])
    """
    if channels is None:
        channels = _default_dvector_channels()
    nc = len(channels)
    lam = np.asarray(couplings, dtype=float)
    Smats = [S for _, _, S in channels]
    Sd = [np.conj(S).T for S in Smats]
    trSS = [np.trace(Sd[a] @ Smats[a]).real for a in range(nc)]
    phitil = [ch[1] for ch in channels]
    om = np.concatenate([omega, -omega])                  # full +/- Matsubara set
    bmax = 0.5 * np.pi * (1.0 - cos_min)
    beta = np.linspace(-bmax, bmax, Nbeta)
    cosb = np.cos(beta)
    w = 1.0 / (2.0 * Nbeta)                               # per branch (full circle uniform)
    if Dbulk is None:
        Dbulk = _bulk_dvector(lam, temp, om, beta, w, phitil)
    Dref = float(np.max(Dbulk))
    if Dref <= 1.0e-6 * temp:
        xi = hvf / (np.pi * max(temp, 1e-12))
        x = np.linspace(0.0, Lxi * xi, int(Lxi * nper))
        return x, np.zeros((nc, len(x)), dtype=np.complex128), Dbulk
    xi = hvf / (np.pi * Dref)
    Ng = int(Lxi * nper)
    x = np.linspace(0.0, Lxi * xi, Ng)
    dx = x[1] - x[0]
    # seed: each channel toward its bulk value (texture emerges from self-consistency)
    Damp = np.array([Dbulk[a] * np.tanh(x / xi) for a in range(nc)], dtype=np.complex128)
    fo = np.array([phitil[a](beta) for a in range(nc)])         # [nc, Nb] outgoing
    fi = np.array([phitil[a](np.pi - beta) for a in range(nc)]) # [nc, Nb] incoming
    for it in range(itemax):
        acc = np.zeros((nc, Ng), dtype=np.complex128)
        for ib in range(Nbeta):
            ds = dx / cosb[ib]
            # unfolded gap matrix path: incoming (x=L..0) then outgoing (x=0..L)
            Dpath = np.zeros((2 * Ng, 2, 2), dtype=np.complex128)
            for a in range(nc):
                amp_path = np.concatenate([fi[a, ib] * Damp[a][::-1], fo[a, ib] * Damp[a]])
                Dpath += amp_path[:, None, None] * Smats[a]
            _, fmat = _mat_batch(om, Dpath, hvf, ds, 0.0)          # [2Ng, Nw, 2, 2]
            f_out = fmat[Ng:]
            f_in = fmat[:Ng][::-1]
            for a in range(nc):
                po = np.einsum('ij,xwji->x', Sd[a], f_out) / trSS[a]   # summed over freq
                pi_ = np.einsum('ij,xwji->x', Sd[a], f_in) / trSS[a]
                acc[a] += w * (np.conj(fo[a, ib]) * po + np.conj(fi[a, ib]) * pi_)
        Damp_new = (lam[:, None] * temp) * acc            # complex (allows TRSB phase)
        err = np.abs(Damp_new - Damp).max() / Dref
        Damp = (1.0 - mix) * Damp + mix * Damp_new
        if err < eps:
            break
    return x, Damp, Dbulk


def surface_dvector_ldos(x, Damp, wlist, channels=None, Dbulk=None, delta=None,
                         Nbeta=16, hvf=1.0, cos_min=0.06):
    """
    @fn surface_dvector_ldos
    @brief Surface local density of states N(x=0, w)/N0 of a converged d-vector
    texture state, via the retarded 2x2 spin-matrix Riccati (w -> w + i delta).
    The sign-changing dominant component gives a zero-energy surface bound state
    (Andreev/ZEBS of the triplet edge); the subdominant / d-vector texture splits
    or shifts it.  Angle-averaged over the trajectory directions.
    @param x, Damp: grid and converged complex amplitudes [Ncomp, Ng] from solve_surface_dvector
    @return N(0, w)/N0 [Nw]
    """
    if channels is None:
        channels = _default_dvector_channels()
    nc = len(channels)
    Smats = [S for _, _, S in channels]
    phitil = [ch[1] for ch in channels]
    Ng = Damp.shape[1]
    Dref = float(np.max(np.abs(Damp))) if Dbulk is None else float(np.max(np.abs(Dbulk)))
    if delta is None:
        delta = 0.03 * Dref
    dx = x[1] - x[0]
    bmax = 0.5 * np.pi * (1.0 - cos_min)
    beta = np.linspace(-bmax, bmax, Nbeta)
    cosb = np.cos(beta)
    w = 1.0 / (2.0 * Nbeta)
    fo = np.array([phitil[a](beta) for a in range(nc)])
    fi = np.array([phitil[a](np.pi - beta) for a in range(nc)])
    Nw = len(wlist)
    zw = delta - 1j * wlist                                     # retarded continuation [Nw]
    ldos = np.zeros(Nw)
    for ib in range(Nbeta):                                     # batch all frequencies per direction
        ds = dx / cosb[ib]
        Dpath = np.zeros((2 * Ng, 2, 2), dtype=np.complex128)
        for a in range(nc):
            amp_path = np.concatenate([fi[a, ib] * Damp[a][::-1], fo[a, ib] * Damp[a]])
            Dpath += amp_path[:, None, None] * Smats[a]
        g, _ = _mat_batch(zw, Dpath, hvf, ds, 0.0)             # [2Ng, Nw, 2, 2]
        tr0 = np.einsum('wii->w', g[Ng]) + np.einsum('wii->w', g[Ng - 1])  # x=0 both branches
        ldos += w * tr0.real / 2.0
    return ldos


def calc_surface_dvector(coupling, temp, wc, kb=1.0, Lxi=8.0, nper=8, Nbeta=16,
                         sub_ratio=0.9, sw_ldos=True):
    """
    @fn calc_surface_dvector
    @brief Driver: self-consistent d-vector texture at a specular surface of a
    triplet superconductor.  The bulk orders in the dominant p_x channel (spin e_x,
    orbital cos(b), sign-changing at an x-surface, lambda); a subdominant p_y channel
    (spin e_z, orbital sin(b), even, lambda*sub_ratio) is also retained.  The dominant
    p_x is pair-broken at the surface, so the (relatively enhanced) p_y component
    rotates the net d-vector toward e_z -- the d-vector texture.  Reports the spatial
    component amplitudes and the rotation angle theta_d(x) = atan2(|Dpy|,|Dpx|), the
    relative phase (90 deg = time-reversal-broken p_x + i p_y surface state), and
    writes 'surface_dvector.dat'.
    @param sub_ratio: lambda_py / lambda_px (subdominant strength)
    """
    omega = matsubara(temp, wc)
    couplings = (coupling, coupling * sub_ratio)
    print("d-vector texture at a specular surface (triplet p_x + subdominant p_y, spin-matrix Riccati)", flush=True)
    print(f"T={temp/kb:.2f} K, lambda_px={couplings[0]:.3f}, lambda_py={couplings[1]:.3f}, "
          f"{len(omega)} Matsubara freqs, Nbeta={Nbeta}, grid={int(Lxi*nper)}", flush=True)
    x, Damp, Dbulk = solve_surface_dvector(couplings, temp, omega, Lxi=Lxi, nper=nper, Nbeta=Nbeta)
    Dref = float(np.max(np.abs(Dbulk)))
    if Dref <= 0.0:
        print("normal state (Dbulk=0); nothing to profile", flush=True)
        return x, Damp
    xi = 1.0 / (np.pi * Dref)
    apx, apy = np.abs(Damp[0]), np.abs(Damp[1])
    theta = np.degrees(np.arctan2(apy, np.maximum(apx, 1e-12)))            # d-vector tilt
    relph = np.degrees(np.angle(Damp[1] / np.where(np.abs(Damp[0]) > 1e-12 * Dref, Damp[0], 1.0)))
    print(f"Dbulk(px,py) = ({Dbulk[0].real:.4e}, {Dbulk[1].real:.4e}) eV,  xi = {xi:.4g}", flush=True)
    print(f"  surface: |Dpx|(0)/Db={apx[0]/Dref:.3f}  |Dpy|(0)/Db={apy[0]/Dref:.3f}  "
          f"theta_d(0)={theta[0]:.1f} deg  rel.phase={relph[0]:.1f} deg", flush=True)
    print(f"  bulk:    |Dpx|(L)/Db={apx[-1]/Dref:.3f}  |Dpy|(L)/Db={apy[-1]/Dref:.3f}  "
          f"theta_d(L)={theta[-1]:.1f} deg", flush=True)
    print(f"  TEXTURE: d-vector tilt theta_d {theta[-1]:.1f} deg (bulk) -> {theta[0]:.1f} deg (surface)", flush=True)
    try:
        with open('surface_dvector.dat', 'w') as fh:
            fh.write("# x/xi  |Dpx|/Db  |Dpy|/Db  theta_d[deg]  rel.phase[deg]\n")
            for j in range(len(x)):
                fh.write(f"{x[j]/xi:10.4f} {apx[j]/Dref:12.5e} {apy[j]/Dref:12.5e} "
                         f"{theta[j]:9.3f} {relph[j]:9.3f}\n")
    except IOError as e:
        print(f"Error writing surface_dvector.dat: {e}", flush=True)
    if sw_ldos:                                                # spectroscopic signature
        wlist = np.linspace(-2.0 * Dref, 2.0 * Dref, 121)
        ldos = surface_dvector_ldos(x, Damp, wlist, Dbulk=Dbulk, Nbeta=Nbeta)
        i0 = np.argmin(np.abs(wlist))
        print(f"  surface LDOS: N(0,0)/N0 = {ldos[i0]:.3f} "
              f"(sign-changing dominant -> zero-energy surface bound state)", flush=True)
        try:
            with open('surface_dvector_ldos.dat', 'w') as fh:
                fh.write("# w/Dbulk   N(0,w)/N0\n")
                for wr, n in zip(wlist, ldos):
                    fh.write(f"{wr/Dref:12.5e} {n:12.5e}\n")
        except IOError as e:
            print(f"Error writing surface_dvector_ldos.dat: {e}", flush=True)
    return x, Damp


# --------------------------------------------------------------------------- #
#  d-vector TEXTURE around a vortex (2D spin-matrix Riccati)
# --------------------------------------------------------------------------- #
def solve_vortex2d_dvector(couplings, temp, omega, channels=None, windings=(1, 0),
                           Dbulk=None, Lxi=7.0, ngrid=33, nbeta=16, hvf=1.0,
                           eps=4.0e-3, itemax=40, mix=0.4):
    """
    @fn solve_vortex2d_dvector
    @brief Self-consistent d-vector TEXTURE around an isolated vortex of a
    multi-component unitary triplet, via the 2x2 spin-matrix Riccati on the 2D plane
    (the vortex analogue of solve_surface_dvector).  Each spin component a has an
    orbital form factor phi_a(beta), a spin matrix Shat_a, a coupling lambda_a and a
    phase WINDING m_a; its complex amplitude field A_a(r) (winding removed) is solved
    self-consistently, Delta_hat(r,k) = sum_a phi_a(k) A_a(r) e^{i m_a theta} Shat_a.
    The dominant (m=1) component vanishes and is pair-broken in the core; a subdominant
    with a different winding (default m=0, core-localized) survives there, so the net
    d-vector reorients in spin space across the core -- the d-vector texture.  (For
    equal windings the axisymmetric profiles coincide and there is no radial texture.)
    @return (xg, A [Ncomp, ngrid, ngrid] complex, Dbulk [Ncomp], xi)
    """
    from scipy.interpolate import RegularGridInterpolator
    from ._eilenberger_vortex import _eval_field
    try:
        from ..flibs import matrix_riccati_chords as _fort_mc
        have_mc = True
    except Exception:
        have_mc = False
    if channels is None:
        channels = _default_dvector_channels()
    nc = len(channels)
    lam = np.asarray(couplings, dtype=float)
    Smats = [S for _, _, S in channels]
    Sd = [np.conj(S).T for S in Smats]
    trSS = [np.trace(Sd[a] @ Smats[a]).real for a in range(nc)]
    phitil = [ch[1] for ch in channels]
    om = np.concatenate([omega, -omega])
    if Dbulk is None:                                   # full-circle bulk amplitudes
        nb0 = 24
        bb = np.linspace(-0.5 * np.pi * 0.94, 0.5 * np.pi * 0.94, nb0)
        Dbulk = _bulk_dvector(lam, temp, om, bb, 1.0 / (2.0 * nb0), phitil)
    Dref = float(np.max(Dbulk))
    xi = hvf / (np.pi * Dref)
    R = Lxi * xi
    xg = np.linspace(-R, R, ngrid)
    dx = xg[1] - xg[0]
    X, Y = np.meshgrid(xg, xg, indexing='ij')
    Rg = np.sqrt(X ** 2 + Y ** 2)
    theta = np.arctan2(Y, X)
    SS, BB = np.meshgrid(xg, xg, indexing='ij')
    dirs = np.linspace(0.0, 2.0 * np.pi, nbeta, endpoint=False)
    wt_dir = 1.0 / nbeta
    phid = np.array([phitil[a](dirs) for a in range(nc)])      # [nc, nbeta]
    mwind = np.asarray(windings, dtype=int)
    # seed: winding-1 components vanish at the core (tanh); core-localized (m=0)
    # subcritical components get a small core-peaked seed so they can nucleate
    A = np.empty((nc, ngrid, ngrid), dtype=np.complex128)
    for a in range(nc):
        if Dbulk[a] > 1.0e-3 * Dref:
            A[a] = Dbulk[a] * (np.tanh(Rg / xi) if mwind[a] != 0 else 1.0)
        else:
            A[a] = 0.15 * Dref / np.cosh(Rg / xi)      # core-localized nucleation seed
    ewind = np.exp(1j * theta)
    for it in range(itemax):
        Ai = [RegularGridInterpolator((xg, xg), A[a], bounds_error=False,
                                      fill_value=complex(Dbulk[a])) for a in range(nc)]
        accf = np.zeros((nc, ngrid, ngrid), dtype=np.complex128)
        for ib in range(nbeta):
            cb, sb = np.cos(dirs[ib]), np.sin(dirs[ib])
            Lx = SS * cb - BB * sb
            Ly = SS * sb + BB * cb
            thr = np.exp(1j * np.arctan2(Ly, Lx))
            sxy = X * cb + Y * sb
            bxy = -X * sb + Y * cb
            Dpath = np.zeros((ngrid, ngrid, 2, 2), dtype=np.complex128)
            for a in range(nc):
                amp = phid[a, ib] * Ai[a]((Lx, Ly)) * thr ** mwind[a]   # [ns,nb]
                Dpath += amp[:, :, None, None] * Smats[a]
            if have_mc:
                _, fch = _fort_mc(om, Dpath, hvf, dx, 0.0)          # [ns,nb,Nw,2,2]
            else:
                fch = np.empty((ngrid, ngrid, len(om), 2, 2), dtype=np.complex128)
                for j in range(ngrid):
                    _, fch[:, j] = _mat_batch(om, Dpath[:, j], hvf, dx, 0.0)
            for a in range(nc):
                fa = np.einsum('ij,snwji->sn', Sd[a], fch) / trSS[a]   # summed over freq
                accf[a] += wt_dir * np.conj(phid[a, ib]) * _eval_field(fa, xg, sxy, bxy, fill=0.0)
        err = 0.0
        for a in range(nc):
            A_new = (lam[a] * temp) * (accf[a] * np.conj(ewind) ** mwind[a])   # remove winding
            err = max(err, np.abs(A_new - A[a]).max() / Dref)
            A[a] = (1.0 - mix) * A[a] + mix * A_new
        if err < eps:
            break
    return xg, A, Dbulk, xi


def calc_vortex_dvector(coupling, temp, wc, kb=1.0, Lxi=7.0, ngrid=33, nbeta=16,
                        sub_ratio=0.95):
    """
    @fn calc_vortex_dvector
    @brief Driver: self-consistent d-vector texture around an isolated vortex of a
    triplet superconductor (dominant p_x(e_x) + subdominant p_y(e_z), 2D spin-matrix
    Riccati).  The dominant component winds and is pair-broken in the core; the
    subdominant is relatively enhanced there, so the net d-vector tilt
    theta_d(r) = atan2(|A_pz|,|A_px|) reorients across the core.  Reports the radial
    texture (along +x) and writes 'vortex_dvector.dat'.
    """
    omega = matsubara(temp, wc)
    couplings = (coupling, coupling * sub_ratio)
    print("d-vector texture around a vortex (triplet p_x + subdominant p_y, 2D spin-matrix Riccati)", flush=True)
    print(f"T={temp/kb:.2f} K, lambda_px={couplings[0]:.3f}, lambda_py={couplings[1]:.3f}, "
          f"{len(omega)} Matsubara freqs, grid={ngrid}x{ngrid}, nbeta={nbeta}", flush=True)
    xg, A, Dbulk, xi = solve_vortex2d_dvector(couplings, temp, omega, Lxi=Lxi,
                                              ngrid=ngrid, nbeta=nbeta)
    Dref = float(np.max(np.abs(Dbulk)))
    if Dref <= 0.0:
        print("normal state (Dbulk=0); nothing to profile", flush=True)
        return xg, A
    ic = ngrid // 2
    r = xg[ic:]                                          # radial cut along +x (y=0)
    apx = np.abs(A[0, ic:, ic])
    apz = np.abs(A[1, ic:, ic])
    theta = np.degrees(np.arctan2(apz, np.maximum(apx, 1e-12)))
    print(f"Dbulk(px,pz) = ({Dbulk[0].real:.4e}, {Dbulk[1].real:.4e}) eV,  xi = {xi:.4g}", flush=True)
    print(f"  core r=0:  |A_px|/Db={apx[0]/Dref:.3f}  |A_pz|/Db={apz[0]/Dref:.3f}  theta_d={theta[0]:.1f} deg", flush=True)
    print(f"  bulk r=R:  |A_px|/Db={apx[-1]/Dref:.3f}  |A_pz|/Db={apz[-1]/Dref:.3f}  theta_d={theta[-1]:.1f} deg", flush=True)
    print(f"  TEXTURE: d-vector tilt theta_d {theta[-1]:.1f} deg (bulk) -> {theta[0]:.1f} deg (core)", flush=True)
    try:
        with open('vortex_dvector.dat', 'w') as fh:
            fh.write("# r/xi  |A_px|/Db  |A_pz|/Db  theta_d[deg]\n")
            for j in range(len(r)):
                fh.write(f"{r[j]/xi:10.4f} {apx[j]/Dref:12.5e} {apz[j]/Dref:12.5e} {theta[j]:9.3f}\n")
    except IOError as e:
        print(f"Error writing vortex_dvector.dat: {e}", flush=True)
    return xg, A
