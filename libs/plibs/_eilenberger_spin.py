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
from ._eilenberger import matsubara, build_fs, BCS_RATIO
from ..flibs import matrix_riccati_batch, matrix_riccati_chords


def _mat_batch(om, Dpath, hvf, ds, h):
    """g, f [Ns, Nw, 2, 2] along one trajectory, batched over frequencies, via the
    Fortran 2x2 spin-matrix Mobius Riccati kernel (the single matrix-Riccati path)."""
    return matrix_riccati_batch(np.ascontiguousarray(om, dtype=np.complex128),
                                np.ascontiguousarray(Dpath, dtype=np.complex128), hvf, ds, h)


def _beta_dirs(nbeta: int) -> np.ndarray:
    """Trajectory directions on the circle that provably never land on a gap node:
    the midpoint grid beta_i = (i + 1/2) * 2 pi / nbeta with an ODD number of points
    (an even nbeta is bumped by one).

    Why it matters.  The plain grid i * 2 pi / nbeta puts a direction exactly on
    beta = pi/2 whenever 4 | nbeta, which is where a p_x (or d_{x2-y2}) gap has a node,
    so that whole chord carries Delta == 0.  The matrix Riccati's bulk root
    a = D / (om + sqrt(om^2 + tr(D D^dag)/2)) is then 0/0 for the NEGATIVE Matsubara
    frequencies this module feeds it -- sqrt is the principal root, so for Re(om) < 0 the
    denominator cancels exactly -- and the NaN spreads along the chord and through the
    whole self-consistency.  Measured on a single p_x channel: nbeta=24 -> NaN,
    nbeta=25 -> the gap heals to 0.988 of its bulk value.

    Why odd + midpoint is enough for EVERY harmonic.  A node of cos(m beta) or
    sin(m beta) sits at beta = (2k+1) pi / (2m).  The midpoint grid is
    beta_i = (2i+1) pi / nbeta.  Equating them gives 2m(2i+1) = nbeta(2k+1): the left
    side is even, so a solution needs nbeta even.  With nbeta odd there is none, for any
    m and any k.  (A half-step offset alone does not do it -- nbeta=12 still samples the
    d-wave node at 45 degrees.)  The midpoint rule is also the better angular quadrature,
    so nothing is given up.
    """
    n = int(nbeta) | 1                                  # odd: never coincides with a node
    return (np.arange(n) + 0.5) * (2.0 * np.pi / n)


def _check_finite(x, where: str):
    """Fail loudly instead of returning NaN.  The known trigger is a trajectory whose
    gap vanishes identically (see _beta_dirs); a self-consistently developed node could
    do it too, and a silent NaN would otherwise propagate into the reported profile."""
    if not np.isfinite(x).all():
        raise FloatingPointError(
            f"{where}: the matrix Riccati returned non-finite values.  The usual cause is "
            "a trajectory on which the gap vanishes identically (an exactly nodal "
            "direction), where the bulk root is 0/0 for negative Matsubara frequencies -- "
            "see _beta_dirs.  Try a different nbeta.")
    return x


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
        damp_init = BCS_RATIO * temp
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
                    gap_sym=0, h_list=None, fs_width=5.0e-3, kb=1.0,
                    gap_orbital=None, delta0=None, gauge=None, spin_map=None):
    """
    @fn calc_spin_pauli
    @brief Compare the Zeeman (Maki) field response of singlet vs triplet pairing
    using the 2x2 spin solver: singlet and triplet d||h are Pauli-limited, triplet
    d perpendicular to h is Pauli-immune (equal-spin pairing).  Sweeps h and writes
    'spin_pauli.dat' (Delta(h)/Delta0 for each channel).
    @param gap_orbital,delta0,gauge,spin_map: gap specification on the 3D Fermi surface,
           exactly as in calc_eilenberger (band-projected orbital gap / per-sheet signs)
    """
    from ._bands import get_emesh
    from ._eilenberger import report_fs_gap
    omega = matsubara(temp, wc)
    Nk, klist, eig, uni, kweight = get_emesh(Nx, Ny, Nz, ham_r, S_r, rvec, avec, sw_uni=True)
    wf, phif, fsinfo = build_fs(eig, klist, mu, gap_sym, fs_width, uni=uni,
                                gap_orbital=gap_orbital, delta0=delta0,
                                gauge=gauge, spin_map=spin_map, sw_band=True)
    report_fs_gap(wf, phif, fsinfo, kb)
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


def chiral_channels(ell: int = 1, dvec: str = 'z'):
    """
    @fn chiral_channels
    @brief Two-component basis of a CHIRAL equal-spin triplet: phi_pm(beta) =
    e^{+-i ell beta} carried by ONE spin matrix (unlike _default_dvector_channels, whose
    two channels sit on DIFFERENT spin matrices and describe a d-vector texture in spin
    space).  Here the relative phase of the two amplitudes IS the chirality, which is
    exactly what the scalar real-amplitude vortex solver cannot represent and why it
    refuses chiral form factors (_reject_chiral_ff).

    ell=1 is p+ip / p-ip, ell=2 is d+id / d-id.  |phi| = 1 on the circle, so each channel
    already satisfies <|phi|^2> = 1 and the two are degenerate in the bulk -- they are
    time-reversed partners of the same irrep, so they take the SAME coupling and the
    chiral state is the one that spontaneously picks a single component.
    @param dvec: spin direction of the equal-spin pairing, 'z' (Shat=sigma_z i sigma_y)
                 or 'x'.  'y' is excluded: Shat_y is proportional to the identity, which
                 the unitary matrix-Riccati bulk seed does not handle.
    @return list of (name, phi(beta), Shat) channels, dominant chirality first
    """
    sx, sy, sz = _pauli()
    isy = 1j * sy
    if dvec == 'z':
        S = sz @ isy
    elif dvec == 'x':
        S = sx @ isy
    else:
        raise ValueError(f"chiral_channels: dvec={dvec!r} (use 'z' or 'x'; 'y' is "
                         "proportional to the identity and is not supported)")
    lab = {1: 'p', 2: 'd', 3: 'f'}.get(int(ell), f'L{int(ell)}')
    return [(f'{lab}+i{lab} (l=+{ell})', (lambda b, m=ell: np.exp(1j * m * b)), S),
            (f'{lab}-i{lab} (l=-{ell})', (lambda b, m=ell: np.exp(-1j * m * b)), S)]


def chiral_windings(m_dom: int = 1, ell: int = 1):
    """
    @fn chiral_windings
    @brief Phase windings of the two chiral_channels around an axisymmetric vortex.

    With Delta(r,k) = sum_pm A_pm(r) e^{+-i ell beta} and A_pm = a_pm(rho) e^{i m_pm theta},
    a rotation by alpha about z sends beta -> beta+alpha and theta -> theta+alpha, so the
    two channels pick up e^{i(+ell + m_+)alpha} and e^{i(-ell + m_-)alpha}.  The order
    parameter can only transform by ONE overall phase, so

        m_- = m_+ + 2 ell.

    A singly quantized vortex parallel to the chirality (m_+ = 1, ell = 1) therefore
    induces the opposite-chirality component with winding 3; the antiparallel vortex
    (m_+ = -1) induces it with winding 1.  Getting this wrong does not merely change the
    answer, it makes the ansatz non-axisymmetric and the self-consistency will not settle.

    @warning solve_vortex2d_dvector CANNOT yet run this.  Measured on the isotropic
    cylinder (lambda=0.5, T=4e-4, Lxi=7, ngrid=33, nbeta=24, itemax=60-80), the solver is
    sound only for the winding pattern it was built and validated for -- one winding-1
    dominant plus a winding-0 subdominant:

        channels (spin matrices)        windings   result
        cos/Sx + sin/Sz  (the default)   (1,0)     dominant heals to 0.991   OK
        cos/Sz + sin/Sz  (same Shat)     (1,0)     dominant heals to 0.986   OK
        cos/Sz + sin/Sz, degenerate      (1,0)     both heal to ~0.99        OK
        e^{+ib}/Sz + e^{-ib}/Sz          (1,3)     bulk (0.54,0.52): drifts to
                                                   the NON-chiral combination and
                                                   never heals to 1
        cos/Sz + sin/Sz, degenerate      (1,1)     NaN
        single channel cos/Sz            (1,)      NaN
        single channel e^{ib}/Sz         (1,)      bulk 0.32, never heals

    So sharing a spin matrix is fine and complex amplitudes are fine; what breaks is any
    configuration without exactly one winding-0 component.  The seed is the first
    suspect -- A[a] = Dbulk[a]*tanh(rho/xi) vanishes LINEARLY at the core, which is right
    only for |m|=1 (a winding-m component must vanish as rho^|m|, and a winding-0 one must
    not vanish at all) -- but that alone does not explain the NaN, so the iteration itself
    has to be looked at before a chiral driver can be trusted.
    @return (m_+, m_-)
    """
    return int(m_dom), int(m_dom) + 2 * int(ell)


def _bulk_dvector(couplings, temp, om, beta, w, phitil, eps=1e-10, itemax=2000, mix=0.5):
    """Coupled bulk amplitudes D_a of a multi-component unitary triplet with per-
    channel couplings.  R(beta) = sqrt(w^2 + sum_b (D_b phi_b)^2) couples the
    channels; gap_a: D_a = lambda_a T sum_n < phi_a^2 D_a / R > over the full
    direction set.  Returns the vector (D_a).  Subcritical channels relax to 0."""
    nc = len(phitil)
    lam = np.asarray(couplings, dtype=float)
    phi = np.array([phitil[a](beta) for a in range(nc)])               # [nc, Nb] outgoing
    phii = np.array([phitil[a](np.pi - beta) for a in range(nc)])      # reflected branch
    D = np.full(nc, BCS_RATIO * temp)
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


def _coupled_bulk_gap(lam, temp, om, phid_b, wfull):
    """Coupled multi-channel bulk amplitudes D_a on a dense (nf-weighted) direction
    set: same channel-coupled gap equation as _bulk_dvector but with a single
    branch and explicit FS weights (used by the vortex / lattice d-vector solvers).
    @param     lam: per-channel couplings [nc]
    @param      om: full +/- Matsubara set
    @param  phid_b: normalized form factors on the dense set [nc, ndense]
    @param   wfull: FS weights on the dense set [ndense] (sum 1)
    """
    nc = len(phid_b)
    a2 = np.abs(phid_b) ** 2                            # [nc, ndense]
    D = np.full(nc, BCS_RATIO * temp)
    for _ in range(2000):
        s2 = ((D[:, None] ** 2) * a2).sum(axis=0)       # sum_b (D_b|phi_b|)^2  [ndense]
        R = np.sqrt(om[None, :] ** 2 + s2[:, None])
        Dn = np.array([max(lam[a] * temp * D[a] * (wfull[:, None] * a2[a][:, None] / R).sum(), 0.0)
                       for a in range(nc)])
        if np.abs(Dn - D).max() < 1e-10 * max(Dn.max(), 1e-12):
            return Dn
        D = 0.5 * (D + Dn)
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
    # positive Matsubara only + factor 2 in the update: matrix_riccati_chords/_batch
    # return the TILDE propagator for Re(om) < 0 (f(-om) = conj(Delta)/R instead of
    # Delta/R), so the full +/- set delivers half the sum with the conjugate gap
    # phase.  See the note in solve_vortex2d_dvector.
    om = omega
    om_full = np.concatenate([omega, -omega])         # bulk gap equation (scalar)
    bmax = 0.5 * np.pi * (1.0 - cos_min)
    beta = np.linspace(-bmax, bmax, Nbeta)
    cosb = np.cos(beta)
    w = 1.0 / (2.0 * Nbeta)                               # per branch (full circle uniform)
    if Dbulk is None:
        Dbulk = _bulk_dvector(lam, temp, om_full, beta, w, phitil)
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
        Damp_new = (2.0 * lam[:, None] * temp) * acc      # complex (allows TRSB phase)
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
                           field=0.0, fs=None, eps=4.0e-3, itemax=40, mix=0.4,
                           A_init=None, bdir=None):
    """
    @fn solve_vortex2d_dvector
    @brief Self-consistent d-vector order parameter of a multi-component unitary
    triplet around a vortex (field=0: isolated vortex) or on the circular-cell vortex
    LATTICE (field=B/Hc2>0: Wigner-Seitz cell of radius Rc=sqrt(2/field)*xi with the
    supercurrent Doppler shift), via the 2x2 spin-matrix Riccati on the 2D plane.
    Each spin component a has an orbital form factor phi_a(beta), spin matrix Shat_a,
    coupling lambda_a and phase WINDING m_a; its complex amplitude field A_a(r) (winding
    removed) is solved self-consistently, Delta_hat(r,k)=sum_a phi_a(k) A_a(r) e^{i m_a theta} Shat_a.
    The dominant (m=1) component vanishes and is pair-broken in the core; a subdominant
    with a different winding (default m=0, core-localized) survives there, so the net
    d-vector reorients in spin space across the core -- the d-vector texture.  With
    field>0 the supercurrent v_F.Q fills the inter-vortex states (Volovik), entering as
    a position-dependent Doppler shift omega -> omega + i v_F.Q in the matrix Riccati.
    @param   bdir: field / vortex-line direction (default c axis).  The lines run along
           B, so the order parameter varies in the plane PERPENDICULAR to it and the
           problem stays 2D; fs_field_frame supplies the trajectory set there, with
           its two axes rescaled by their rms velocities.  Note the two frames are
           different objects: the TRAJECTORY direction is in the field frame, while
           the channel form factors phi(beta) are functions of the CRYSTAL in-plane
           angle and keep using it.  Needs a 3D FS for anything but B || c.
    @param A_init: explicit initial amplitudes [Ncomp, ngrid, ngrid], overriding the
           built-in seed.  Needed to ask whether a given state is a fixed point AT ALL
           (start exactly on it and see whether one iteration moves off it), which the
           built-in seed cannot express -- it always injects a nucleation amplitude into
           every subcritical channel.
    @return (xg, A [Ncomp, ngrid, ngrid] complex, Dbulk [Ncomp], xi)
    """
    from ._eilenberger_vortex import _eval_field, _doppler_chord
    if channels is None:
        channels = _default_dvector_channels()
    nc = len(channels)
    lam = np.asarray(couplings, dtype=float)
    Smats = [S for _, _, S in channels]
    Sd = [np.conj(S).T for S in Smats]
    trSS = [np.trace(Sd[a] @ Smats[a]).real for a in range(nc)]
    phitil = [ch[1] for ch in channels]
    # POSITIVE Matsubara frequencies only, with the factor 2 restored in the gap update.
    # matrix_riccati_chords returns the TILDE propagator for Re(om) < 0 -- measured on a
    # homogeneous chiral gap Delta = D e^{i beta} S: f(+om) = +Delta/R (correct) but
    # f(-om) = conj(Delta)/R, and g(-om) = -conj(g(+om)).  Feeding it the full +/- set
    # therefore delivers exactly HALF of the Matsubara sum with the CONJUGATE gap phase,
    # which for a chiral order parameter lands in the OPPOSITE chirality channel: starting
    # from a pure chiral state one iteration produced |A_-|/D = 0.211 = 0.4 (the mixing) x
    # 0.53, i.e. half the sum, misdirected.  It is invisible for a real gap (conj is the
    # identity), which is why every existing d-vector result is unaffected.  f is even in
    # om, so summing om > 0 and doubling is the same thing -- and it is the convention the
    # scalar solvers already use.
    om_full = np.concatenate([omega, -omega])          # bulk gap equation (scalar, exact)
    # ...and it is only valid while the frequencies stay REAL.  With the supercurrent
    # Doppler shift (field > 0) the argument becomes om + i v_F.Q, so +om and -om are no
    # longer related by conjugation and the reduction does not apply; there the +/- set is
    # kept and the chiral channel decomposition is still contaminated.  A proper fix
    # belongs in the kernel (return f, not f-tilde, for Re(om) < 0).
    zero_field = (field <= 0.0)
    om = omega if zero_field else om_full
    fac = 2.0 if zero_field else 1.0
    # trajectory directions, per-direction |v_F| and FS weight, and the orbital form
    # factor of each channel (evaluated at the k-direction); cylinder or a model/Wannier FS
    def _norm(ph, wt):                                  # normalize nf-weighted <|phi|^2>=1
        out = ph.copy()
        for a in range(nc):
            nrm = np.sqrt((wt * np.abs(out[a]) ** 2).sum())
            if nrm > 0:
                out[a] = out[a] / nrm
        return out
    # dense FS set for an accurate bulk gap, independent of the trajectory sampling
    if fs is not None:
        kfull = np.arctan2(fs['ky'], fs['kx'])
        wfull = np.asarray(fs['nf'], dtype=float)
    else:
        bfull = np.linspace(0.0, 2.0 * np.pi, 180, endpoint=False)
        kfull = bfull
        wfull = np.full(180, 1.0 / 180)
    phid_b = _norm(np.array([phitil[a](kfull) for a in range(nc)], dtype=np.complex128), wfull)
    # trajectory directions, per-direction |v_F| and FS weight, and the orbital form
    # factor of each channel (evaluated at the k-direction); cylinder or a model/Wannier FS
    if fs is not None:                                  # real FS: v_hat directions, nf weights
        from ._eilenberger import fs_field_frame
        frame = fs_field_frame(fs, bdir)
        if field > 0.0 and abs(frame['aniso_ratio'] - 1.0) > 0.05:
            raise NotImplementedError(
                f"solve_vortex2d_dvector: a finite field with an anisotropic vortex "
                f"plane (xi_1/xi_2 = {frame['aniso_ratio']:.3f}) is not supported yet; "
                "use field=0 (isolated vortex) for a tilted / in-plane field.")
        dirs = frame['dirs']                           # trajectory in the plane perp to B
        # fs_field_frame may drop FS points whose velocity is parallel to B, so the
        # crystal angles and the weights have to follow the same mask
        kang = kfull[frame['keep']]                    # CRYSTAL angle sets the form factor
        hvfarr = np.asarray(frame['vabs'], dtype=float)
        wt_dir = np.asarray(frame['nf'], dtype=float)
        if len(dirs) > nbeta:
            # Stratify by TRAJECTORY ANGLE and keep the heaviest member of each bin,
            # carrying that bin's total weight.  The old stride over the raw index
            # order is fine for an angle-ordered 2D contour but meaningless for a 3D
            # FS, whose points are ordered by (k_z slice, position along the contour):
            # 16 evenly spaced indices there sample a couple of k_z slices and miss the
            # rest of the Fermi surface entirely.
            # The key is (trajectory angle, CRYSTAL angle), not the trajectory angle
            # alone.  For B || c the two coincide and this is plain angular binning; for
            # a tilted or in-plane field they decouple completely, and a bin of fixed
            # trajectory direction then holds FS points spread over the whole crystal
            # angle -- i.e. over the whole range of the channel form factors phi(beta).
            # Keeping one representative of such a bin and giving it the bin's entire
            # weight is what made the in-plane chiral vortex diverge (A_+ ran to 2.3 and
            # the induced component to 5.3 times the bulk gap).
            nk = 8
            ib = np.minimum(((dirs + np.pi) / (2.0 * np.pi) * nbeta).astype(np.int64),
                            nbeta - 1)
            ik = np.minimum(((kang + np.pi) / (2.0 * np.pi) * nk).astype(np.int64), nk - 1)
            keys, inv = np.unique(ib * nk + ik, return_inverse=True)
            wsum = np.bincount(inv, weights=wt_dir)
            sel = np.array([np.where(inv == k)[0][np.argmax(wt_dir[inv == k])]
                            for k in range(len(keys))])
            dirs, kang = dirs[sel], kang[sel]
            hvfarr = hvfarr[sel]
            wt_dir = wsum / wsum.sum()                  # bin weights, normalized
        nbeta = len(dirs)
        phid = _norm(np.array([phitil[a](kang) for a in range(nc)], dtype=np.complex128), wt_dir)
        hvf_eff = frame['hvf_eff']
    else:                                              # isotropic cylinder (v_hat = k_hat)
        dirs = _beta_dirs(nbeta)
        nbeta = len(dirs)                              # _beta_dirs rounds up to odd
        hvfarr = np.full(nbeta, hvf)
        wt_dir = np.full(nbeta, 1.0 / nbeta)
        phid = np.array([phitil[a](dirs) for a in range(nc)], dtype=np.complex128)
        hvf_eff = hvf
    if Dbulk is None:                                   # coupled multi-channel bulk gap (dense set)
        Dbulk = _coupled_bulk_gap(lam, temp, om_full, phid_b, wfull)
    Dref = float(np.max(Dbulk))
    xi = hvf_eff / (np.pi * Dref)
    Rc = np.sqrt(2.0 / field) * xi if field > 0.0 else np.inf   # Wigner-Seitz cell radius
    R = Rc if field > 0.0 else Lxi * xi
    xg = np.linspace(-R, R, ngrid)
    dx = xg[1] - xg[0]
    rho_min = 0.5 * dx
    X, Y = np.meshgrid(xg, xg, indexing='ij')
    Rg = np.sqrt(X ** 2 + Y ** 2)
    theta = np.arctan2(Y, X)
    SS, BB = np.meshgrid(xg, xg, indexing='ij')
    mwind = np.asarray(windings, dtype=int)
    # seed: winding-1 components vanish at the core (tanh); core-localized (m=0)
    # subcritical components get a small core-peaked seed so they can nucleate
    A = np.empty((nc, ngrid, ngrid), dtype=np.complex128)
    tanh_r = np.tanh(Rg / xi)
    for a in range(nc):
        # A winding-m amplitude must vanish as rho^|m| at the core (only m=0 survives
        # there), so the seed carries tanh^|m|.  The old seed was tanh^1 for every m != 0
        # and core-PEAKED for any subcritical channel, which is right for the m=(1,0)
        # dominant+subdominant pair it was written for and wrong for anything else --
        # a winding-3 induced component was being started with a core antinode.
        prof = tanh_r ** abs(int(mwind[a])) if mwind[a] != 0 else 1.0
        if Dbulk[a] > 1.0e-3 * Dref:
            A[a] = Dbulk[a] * prof
        else:                                          # subcritical / induced: nucleation seed
            A[a] = 0.15 * Dref * prof / np.cosh(Rg / xi)
    if A_init is not None:
        A = np.array(A_init, dtype=np.complex128).reshape(nc, ngrid, ngrid).copy()
    ewind = np.exp(1j * theta)
    for it in range(itemax):
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
                amp = (phid[a, ib] * _eval_field(A[a], xg, Lx, Ly, fill=complex(Dbulk[a]))
                       * thr ** mwind[a])                               # [ns,nb]
                Dpath += amp[:, :, None, None] * Smats[a]
            if field > 0.0:           # supercurrent Doppler om -> om + i v_F.Q (position dependent)
                dop = hvfarr[ib] * _doppler_chord(Lx, Ly, Rc, cb, sb, rho_min)
                om_ch = om[None, None, :] + 1j * dop[:, :, None]    # [ns,nb,Nw]
            else:
                om_ch = om                                          # [Nw] (broadcast in wrapper)
            _, fch = matrix_riccati_chords(om_ch, np.ascontiguousarray(Dpath),
                                           hvfarr[ib], dx, 0.0)     # [ns,nb,Nw,2,2]
            _check_finite(fch, 'solve_vortex2d_dvector')
            for a in range(nc):
                fa = np.einsum('ij,snwji->sn', Sd[a], fch) / trSS[a]   # summed over freq
                accf[a] += wt_dir[ib] * np.conj(phid[a, ib]) * _eval_field(fa, xg, sxy, bxy, fill=0.0)
        err = 0.0
        for a in range(nc):
            # factor 2: only om > 0 is summed (f is even in om), see the note above
            A_new = (fac * lam[a] * temp) * (accf[a] * np.conj(ewind) ** mwind[a])  # remove winding
            err = max(err, np.abs(A_new - A[a]).max() / Dref)
            A[a] = (1.0 - mix) * A[a] + mix * A_new
        if err < eps:
            break
    return xg, A, Dbulk, xi


def solve_lattice_sc_dvector(couplings, temp, omega, channels=None, windings=(1, 0),
                             field=0.2, lattice='square', Ng=18, nbeta=16, hvf=1.0,
                             fs=None, kappa=None, Lchord=6.0, ds_xi=0.3,
                             eps=4.0e-3, itemax=60, mix=0.4, anderson=True, m_and=4):
    """
    @fn solve_lattice_sc_dvector
    @brief Self-consistent d-vector (multi-component unitary triplet) on the TRUE
    periodic vortex lattice, je-style (formulation A): the spin-matrix order parameter
    Delta_hat(r,k)=sum_a phi_a(k) A_a(r) e^{i m_a chi(r)} Shat_a keeps the analytic
    Abrikosov winding chi(r), and the gap equation is closed by per-grid-point anchored
    2x2 matrix-Riccati trajectories (exact f-to-grid map, no binning) on the magnetic
    cell (square/triangular, _cell_geom).  The dominant winding component (m=1) carries
    the vortex lattice and vanishes at every core; a subdominant with a different winding
    (default m=0) nucleates in the cores -> the d-vector reorients across each core.
    Finite ``kappa`` adds the smooth London Doppler -v_F.A(r) (screening, je #2).
    @return state dict (X,Y,A[nc,Ng,Ng],chi,Dbulk[nc],xi,geom,kappa,...) for lattice_dos_sc_dvector
    """
    from ._eilenberger_vortex import (_cell_geom, _abrikosov_unit_phase, _abrikosov_z,
                                      _periodic_eval, _london_A)
    if channels is None:
        channels = _default_dvector_channels()
    nc = len(channels)
    lam = np.asarray(couplings, dtype=float)
    Smats = [S for _, _, S in channels]
    Sd = [np.conj(S).T for S in Smats]
    trSS = [np.trace(Sd[a] @ Smats[a]).real for a in range(nc)]
    phitil = [ch[1] for ch in channels]
    # NOTE the +/- Matsubara set is kept here even though the kernel returns the TILDE
    # propagator for Re(om) < 0 (see solve_vortex2d_dvector): this solver always carries
    # the supercurrent Doppler shift om -> om + i v_F.Q, so +om and -om are not related by
    # conjugation and the positive-only reduction used there is not available.  Harmless
    # for a real order parameter, which is every configuration this solver is used for;
    # a CHIRAL lattice would need the kernel fixed first.
    om = om_full = np.concatenate([omega, -omega])
    mwind = np.asarray(windings, dtype=int)

    def _norm(ph, wt):                                  # normalize nf-weighted <|phi|^2>=1
        out = ph.copy()
        for a in range(nc):
            nrm = np.sqrt((wt * np.abs(out[a]) ** 2).sum())
            if nrm > 0:
                out[a] = out[a] / nrm
        return out
    if fs is not None:                                  # dense FS for the bulk gap
        kfull = np.arctan2(fs['ky'], fs['kx']); wfull = np.asarray(fs['nf'], float)
    else:
        kfull = np.linspace(0.0, 2.0 * np.pi, 180, endpoint=False); wfull = np.full(180, 1.0 / 180)
    phid_b = _norm(np.array([phitil[a](kfull) for a in range(nc)], dtype=np.complex128), wfull)
    if fs is not None:                                  # trajectory directions / weights
        dirs = np.arctan2(fs['vy'], fs['vx']); kang = kfull
        hvfarr = np.asarray(fs['vabs'], float); wt_dir = wfull
        if len(dirs) > nbeta:
            sel = np.linspace(0, len(dirs) - 1, nbeta).round().astype(int)
            dirs, kang, hvfarr = dirs[sel], kang[sel], hvfarr[sel]
            wt_dir = wt_dir[sel] / wt_dir[sel].sum()
        nbeta = len(dirs)
        phid = _norm(np.array([phitil[a](kang) for a in range(nc)], dtype=np.complex128), wt_dir)
        hvf_eff = float((np.asarray(fs['nf']) * np.asarray(fs['vabs'])).sum())
    else:
        dirs = _beta_dirs(nbeta)
        nbeta = len(dirs)                              # _beta_dirs rounds up to odd
        hvfarr = np.full(nbeta, hvf); wt_dir = np.full(nbeta, 1.0 / nbeta)
        phid = np.array([phitil[a](dirs) for a in range(nc)], dtype=np.complex128)
        hvf_eff = hvf
    Dbulk = _coupled_bulk_gap(lam, temp, om_full, phid_b, wfull)   # coupled multi-channel bulk gap
    Dref = float(np.max(Dbulk))
    if Dref <= 0:
        return None
    xi = hvf_eff / (np.pi * Dref)
    g = _cell_geom(field, xi, lattice, nflux=1)
    fax = (np.arange(Ng) + 0.5) / Ng - 0.5
    F1, F2 = np.meshgrid(fax, fax, indexing='ij')
    X = F1 * g['a1'][0] + F2 * g['a2'][0]
    Y = F1 * g['a1'][1] + F2 * g['a2'][1]
    chi1 = _abrikosov_unit_phase(X, Y, g, 1)           # winding-1 phase on the grid
    L = Lchord * xi; ds = ds_xi * xi
    ns = int(2 * L / ds) | 1; ic = ns // 2
    s = np.linspace(-L, L, ns)
    ax = X.ravel(); ay = Y.ravel(); Nanch = ax.size
    lamL = kappa * xi if kappa is not None else None
    Axg, Ayg = _london_A(g, lamL, 1) if lamL is not None else (None, None)
    chord = []                                          # fixed per-direction geometry
    for ib in range(nbeta):
        cb, sb = np.cos(dirs[ib]), np.sin(dirs[ib])
        px = ax[None, :] + s[:, None] * cb
        py = ay[None, :] + s[:, None] * sb
        eichi = _abrikosov_unit_phase(px, py, g, 1)
        dch = (-(cb * _periodic_eval(Axg, g, px, py) + sb * _periodic_eval(Ayg, g, px, py))
               if lamL is not None else None)
        chord.append((px, py, eichi, dch))
    Zg = np.abs(_abrikosov_z(X, Y, g))                  # ~0 at cores, ~1 between
    A = np.empty((nc, Ng, Ng), dtype=np.complex128)
    for a in range(nc):
        if Dbulk[a] > 1.0e-3 * Dref:
            A[a] = Dbulk[a] * (np.tanh(Zg / 0.5) if mwind[a] != 0 else 1.0)
        else:
            A[a] = 0.15 * Dref * (1.0 - np.tanh(Zg / 0.5))   # core-peaked nucleation seed
    Nw = len(om)
    om_dir = [om if dch is None
              else np.ascontiguousarray(om[None, None, :] + 1j * hvfarr[ib] * dch[:, :, None])
              for ib, (_, _, _, dch) in enumerate(chord)]
    def gap_map_dv(Acur):                              # one pass: A[nc] -> A_new[nc]
        accf = np.zeros((nc, Nanch), dtype=np.complex128)
        for ib in range(nbeta):
            px, py, eichi, _ = chord[ib]
            Dpath = np.zeros((ns, Nanch, 2, 2), dtype=np.complex128)
            for a in range(nc):
                amp = phid[a, ib] * _periodic_eval(Acur[a], g, px, py) * eichi ** mwind[a]
                Dpath += amp[:, :, None, None] * Smats[a]
            _, fch = matrix_riccati_chords(om_dir[ib], np.ascontiguousarray(Dpath),
                                           hvfarr[ib], ds, 0.0)
            fc = fch[ic]                                # [Nanch, Nw, 2, 2]
            for a in range(nc):
                tmp = np.einsum('ij,nwji->n', Sd[a], fc) / trSS[a]   # spin-project, sum freq
                accf[a] += wt_dir[ib] * np.conj(phid[a, ib]) * np.conj(eichi[ic] ** mwind[a]) * tmp
        return np.array([(lam[a] * temp) * accf[a].reshape(Ng, Ng) for a in range(nc)])

    xs, fs_h = [], []                                  # Anderson/Pulay history (stacked channels)
    sh = A.shape
    for it in range(itemax):
        res = gap_map_dv(A) - A
        err = np.abs(res).max() / Dref
        if err < eps:
            A = A + res
            break
        xs.append(A.ravel().copy()); fs_h.append(res.ravel().copy())
        if len(fs_h) > m_and + 1:
            xs.pop(0); fs_h.pop(0)
        if not anderson or len(fs_h) == 1:
            A = A + mix * res
        else:
            dF = np.column_stack([fs_h[i + 1] - fs_h[i] for i in range(len(fs_h) - 1)])
            dX = np.column_stack([xs[i + 1] - xs[i] for i in range(len(xs) - 1)])
            gam, *_ = np.linalg.lstsq(dF, fs_h[-1], rcond=None)
            A = (xs[-1] + mix * fs_h[-1] - (dX + mix * dF) @ gam).reshape(sh)
    return dict(X=X, Y=Y, A=A, chi=np.angle(chi1), Dbulk=Dbulk, xi=xi, geom=g,
                a=np.sqrt(g['S']), kappa=kappa, channels=channels, windings=mwind,
                iters=it + 1, err=err)


def lattice_dos_sc_dvector(state, wlist, delta=None, nbeta=16, hvf=1.0, fs=None,
                           Lchord=6.0, ds_xi=0.3):
    """
    @fn lattice_dos_sc_dvector
    @brief Spatially-averaged DOS N(w)/N0 of the self-consistent d-vector periodic
    lattice (state from solve_lattice_sc_dvector): per-grid-point anchored 2x2 matrix
    Riccati on the retarded axis (z=delta-i*w), N(w)/N0 = <(1/2) Tr Re g(anchor)>_{grid,FS},
    with the same Abrikosov winding and (finite kappa) London Doppler as the solver.
    @return N(w)/N0 [Nw]
    """
    from ._eilenberger_vortex import _abrikosov_unit_phase, _periodic_eval, _london_A
    A = state['A']; g = state['geom']; xi = state['xi']
    Dbulk = state['Dbulk']; Dref = float(np.max(Dbulk)); Ng = A.shape[1]
    X, Y = state['X'], state['Y']
    channels = state['channels']; mwind = state['windings']
    Smats = [S for _, _, S in channels]; phitil = [ch[1] for ch in channels]; nc = len(channels)
    kappa = state.get('kappa'); lam = kappa * xi if kappa is not None else None
    Axg, Ayg = _london_A(g, lam, 1) if lam is not None else (None, None)
    if delta is None:
        delta = 0.03 * Dref
    zomega = delta - 1j * np.asarray(wlist); Nw = zomega.size
    if fs is not None:
        dirs = np.arctan2(fs['vy'], fs['vx']); kang = np.arctan2(fs['ky'], fs['kx'])
        hvfarr = np.asarray(fs['vabs'], float); wt_dir = np.asarray(fs['nf'], float)
        if len(dirs) > nbeta:
            sel = np.linspace(0, len(dirs) - 1, nbeta).round().astype(int)
            dirs, kang, hvfarr = dirs[sel], kang[sel], hvfarr[sel]
            wt_dir = wt_dir[sel] / wt_dir[sel].sum()
        nbeta = len(dirs)
    else:
        dirs = _beta_dirs(nbeta); kang = dirs
        nbeta = len(dirs)                              # _beta_dirs rounds up to odd
        hvfarr = np.full(nbeta, hvf); wt_dir = np.full(nbeta, 1.0 / nbeta)
    phid = np.array([phitil[a](kang) for a in range(nc)], dtype=np.complex128)
    nrm = np.array([np.sqrt((wt_dir * np.abs(phid[a]) ** 2).sum()) for a in range(nc)])
    phid = phid / np.where(nrm[:, None] > 0, nrm[:, None], 1.0)
    L = Lchord * xi; ds = ds_xi * xi
    ns = int(2 * L / ds) | 1; ic = ns // 2
    s = np.linspace(-L, L, ns)
    ax = X.ravel(); ay = Y.ravel(); Nanch = ax.size
    gsum = np.zeros(Nw)
    for ib in range(nbeta):
        cb, sb = np.cos(dirs[ib]), np.sin(dirs[ib])
        px = ax[None, :] + s[:, None] * cb; py = ay[None, :] + s[:, None] * sb
        eichi = _abrikosov_unit_phase(px, py, g, 1)
        Dpath = np.zeros((ns, Nanch, 2, 2), dtype=np.complex128)
        for a in range(nc):
            amp = phid[a, ib] * _periodic_eval(A[a], g, px, py) * eichi ** mwind[a]
            Dpath += amp[:, :, None, None] * Smats[a]
        if lam is None:
            om_ch = zomega
        else:
            dch = -(cb * _periodic_eval(Axg, g, px, py) + sb * _periodic_eval(Ayg, g, px, py))
            om_ch = zomega[None, None, :] + 1j * hvfarr[ib] * dch[:, :, None]
        gch, _ = matrix_riccati_chords(om_ch, np.ascontiguousarray(Dpath), hvfarr[ib], ds, 0.0)
        # N(w) = (1/2) Tr Re g at the anchor, averaged over anchors
        trg = 0.5 * (gch[ic, :, :, 0, 0] + gch[ic, :, :, 1, 1]).real   # [Nanch, Nw]
        gsum += wt_dir[ib] * trg.mean(axis=0)
    return gsum


def calc_vortex_dvector(coupling, temp, wc, kb=1.0, Lxi=7.0, ngrid=33, nbeta=16,
                        sub_ratio=0.95, field=0.0, fs=None, bdir=None):
    """
    @fn calc_vortex_dvector
    @brief Driver: self-consistent d-vector order parameter of a triplet superconductor
    around an isolated vortex (field=0) or on the circular-cell vortex LATTICE
    (field=B/Hc2>0, with the supercurrent Doppler shift) -- dominant p_x(e_x) +
    subdominant p_y(e_z), 2D spin-matrix Riccati.  The dominant component winds and is
    pair-broken in the core; the subdominant is relatively enhanced there, so the net
    d-vector tilt theta_d(r) = atan2(|A_pz|,|A_px|) reorients across the core.  Reports
    the radial texture (along +x) and writes 'vortex_dvector.dat'.
    """
    omega = matsubara(temp, wc)
    couplings = (coupling, coupling * sub_ratio)
    cell = 'isolated vortex' if field <= 0 else f'circular-cell lattice B/Hc2={field:.3f}'
    fstxt = ', Wannier FS+v_F' if fs is not None else ''
    print(f"d-vector vortex ({cell}{fstxt}): triplet p_x + subdominant p_y, 2D spin-matrix Riccati", flush=True)
    print(f"T={temp/kb:.2f} K, lambda_px={couplings[0]:.3f}, lambda_py={couplings[1]:.3f}, "
          f"{len(omega)} Matsubara freqs, grid={ngrid}x{ngrid}", flush=True)
    xg, A, Dbulk, xi = solve_vortex2d_dvector(couplings, temp, omega, Lxi=Lxi,
                                              ngrid=ngrid, nbeta=nbeta, field=field,
                                              fs=fs, bdir=bdir)
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


def calc_vortex_chiral(coupling, temp, wc, ell=1, m_dom=1, dvec='z', kb=1.0, Lxi=8.0,
                       ngrid=41, nbeta=24, fs=None, eps=4.0e-3, itemax=100, bdir=None):
    """
    @fn calc_vortex_chiral
    @brief Driver: self-consistent isolated vortex of a CHIRAL equal-spin triplet
    (p+ip for ell=1, d+id for ell=2), via the multi-component spin-matrix Riccati.

    This is the case the scalar vortex solver refuses (_reject_chiral_ff): a chiral order
    parameter has a genuinely complex amplitude and the core induces the OPPOSITE
    chirality, neither of which a single real scalar field can carry.  The two components
    are e^{+-i ell beta} on ONE spin matrix (chiral_channels), their windings differ by
    2*ell (chiral_windings), and their complex amplitudes are solved self-consistently.

    The bulk reference is deliberately NOT the degenerate two-channel gap equation: the
    two chiralities are time-reversed partners of one irrep and share a coupling, so that
    equation returns them with equal amplitude -- the non-chiral (nematic) combination.
    The chiral state has one component at the isotropic-gap value and the other exactly
    zero, and it IS a fixed point of the solver (verified: starting on it, |A_-|/D stays
    at 4e-4 and A_+ holds 0.996 of the bulk).

    Measured for p+ip, m=(1,3), lambda=0.5, T=4e-4, Lxi=8, ngrid=41: A_+ rises from
    0.0003 at the core to 0.98 at the edge, and the induced opposite chirality peaks at
    |A_-|/Db = 0.073 around r = 3 xi, vanishing both at the core (winding 3 -> rho^3) and
    in the bulk.
    @param     ell: chirality (1 = p+ip, 2 = d+id)
    @param   m_dom: winding of the dominant component (+1 parallel to the chirality,
                    -1 antiparallel); the induced one follows chiral_windings
    @param    dvec: equal-spin direction, 'z' or 'x' (see chiral_channels)
    @param      fs: model / Wannier FS dict (None = isotropic cylinder)
    @return (xg, A): grid axis and the complex amplitudes [2, ngrid, ngrid]
    """
    omega = matsubara(temp, wc)
    om_full = np.concatenate([omega, -omega])
    channels = chiral_channels(ell, dvec)
    mw = chiral_windings(m_dom, ell)
    nd = 180
    # |phi_pm| = 1, so the dominant chirality obeys the isotropic-gap equation
    Dchi = float(_coupled_bulk_gap(np.array([coupling]), temp, om_full,
                                   np.ones((1, nd)), np.full(nd, 1.0 / nd))[0])
    fstxt = ', model/Wannier FS+v_F' if fs is not None else ', isotropic cylinder'
    print(f"chiral vortex (isolated{fstxt}): {channels[0][0]} dominant + induced "
          f"{channels[1][0]}, 2D spin-matrix Riccati", flush=True)
    print(f"T={temp/kb:.2f} K, lambda={coupling:.3f} (both chiralities, degenerate), "
          f"windings m=({mw[0]},{mw[1]}), {len(omega)} Matsubara freqs, "
          f"grid={ngrid}x{ngrid}", flush=True)
    if Dchi <= 0.0:
        print("normal state (Dbulk=0); nothing to profile", flush=True)
        return None, None
    xg, A, Dbulk, xi = solve_vortex2d_dvector((coupling, coupling), temp, omega,
                                              channels=channels, windings=mw,
                                              Dbulk=np.array([Dchi, 0.0]), Lxi=Lxi,
                                              ngrid=ngrid, nbeta=nbeta, field=0.0,
                                              fs=fs, eps=eps, itemax=itemax, bdir=bdir)
    ic = ngrid // 2
    r = xg[ic:]                                          # radial cut along +x (y=0)
    ap = np.abs(A[0, ic:, ic]) / Dchi
    am = np.abs(A[1, ic:, ic]) / Dchi
    jm = int(np.argmax(am))
    print(f"Dbulk = {Dchi:.4e} eV,  xi = {xi:.4g}", flush=True)
    print(f"  core  r=0:  |A_+|/Db={ap[0]:.4f}   |A_-|/Db={am[0]:.4f}", flush=True)
    print(f"  edge  r=R:  |A_+|/Db={ap[-1]:.4f}   |A_-|/Db={am[-1]:.4f}", flush=True)
    print(f"  induced opposite chirality: max |A_-|/Db = {am[jm]:.4f} at r/xi = "
          f"{r[jm]/xi:.2f}  (0 at the core and in the bulk = chiral state preserved)",
          flush=True)
    try:
        with open('vortex_chiral.dat', 'w') as fh:
            fh.write("# r/xi  |A_+|/Db  |A_-|/Db   (dominant / induced opposite chirality)\n")
            for j in range(len(r)):
                fh.write(f"{r[j]/xi:10.4f} {ap[j]:12.5e} {am[j]:12.5e}\n")
    except IOError as e:
        print(f"Error writing vortex_chiral.dat: {e}", flush=True)
    return xg, A


def calc_vortex_lattice_sc_dvector(coupling, temp, wc, kb=1.0, field=0.2, lattice='square',
                                   Ng=12, nbeta=10, sub_ratio=0.92, kappa=None, fs=None):
    """
    @fn calc_vortex_lattice_sc_dvector
    @brief Driver: self-consistent d-vector triplet on the TRUE periodic vortex lattice
    (je-style formulation A, square/triangular cell) -- dominant p_x(e_x) winding with
    the Abrikosov lattice (node at every core) + subdominant p_y(e_z) that nucleates,
    core-localized, where the dominant is suppressed (the d-vector texture).  Finite
    ``kappa`` adds the London screening Doppler.  Writes the converged fields to
    'lattice_dvector.dat'.
    """
    omega = matsubara(temp, wc)
    couplings = (coupling, coupling * sub_ratio)
    print(f"d-vector periodic vortex lattice (formulation A, {lattice}, B/Hc2={field:.3f}"
          f"{', finite kappa=%g' % kappa if kappa is not None else ', extreme'}"
          f"{', Wannier FS' if fs is not None else ''}): dominant p_x(e_x) + subdominant p_y(e_z)",
          flush=True)
    print(f"T={temp/kb:.2f} K, lambda=({couplings[0]:.3f},{couplings[1]:.3f}), "
          f"{len(omega)} Matsubara freqs, grid={Ng}x{Ng}", flush=True)
    st = solve_lattice_sc_dvector(couplings, temp, omega, windings=(1, 0), field=field,
                                  lattice=lattice, Ng=Ng, nbeta=nbeta, kappa=kappa, fs=fs)
    if st is None:
        print("normal state (Dbulk=0)", flush=True)
        return None
    A = st['A']; Db = st['Dbulk']; Dref = float(np.max(Db)); xi = st['xi']
    X, Y = st['X'], st['Y']; r = np.sqrt(X ** 2 + Y ** 2)
    core = r < 0.8 * xi; bulk = r > 2.0 * xi
    adom_core, adom_bulk = np.abs(A[0])[core].mean() / Dref, np.abs(A[0])[bulk].mean() / Dref
    asub_core = np.abs(A[1])[core].mean() / Dref
    asub_bulk = np.abs(A[1])[bulk].mean() / Dref if bulk.any() else 0.0
    print(f"Dbulk(dom,sub) = ({Db[0]:.4e}, {Db[1]:.4e}) eV,  xi = {xi:.4g},  iters={st['iters']}", flush=True)
    print(f"  dominant p_x(e_x):  min|A|/Db={np.abs(A[0]).min()/Dref:.3f} (core node)  "
          f"core={adom_core:.3f}  bulk={adom_bulk:.3f}", flush=True)
    print(f"  subdom  p_y(e_z):   core={asub_core:.4f}  bulk={asub_bulk:.4f}  "
          f"({'core-localized TEXTURE' if asub_core > 1.3 * max(asub_bulk, 1e-9) else 'uniform/absent'})",
          flush=True)
    n0 = float(lattice_dos_sc_dvector(st, np.array([0.0]), nbeta=max(nbeta, 20),
                                      delta=0.05 * Dref, fs=fs)[0])
    print(f"  spatially-averaged zero-energy DOS <N(0)>/N0 = {n0:.3f} "
          f"(triplet: core-bound-state + nodal contributions)", flush=True)
    try:
        with open('lattice_dvector.dat', 'w') as fh:
            fh.write("# x/xi  y/xi  |A_px|/Db  |A_pz|/Db\n")
            for i in range(Ng):
                for j in range(Ng):
                    fh.write(f"{X[i,j]/xi:10.4f} {Y[i,j]/xi:10.4f} "
                             f"{np.abs(A[0,i,j])/Dref:12.5e} {np.abs(A[1,i,j])/Dref:12.5e}\n")
                fh.write("\n")
    except IOError as e:
        print(f"Error writing lattice_dvector.dat: {e}", flush=True)
    return st
