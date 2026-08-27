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

BCS_RATIO = 1.764        # weak-coupling BCS gap ratio Delta0 / (kB Tc)


# --------------------------------------------------------------------------- #
#  Pair-partner gauge
# --------------------------------------------------------------------------- #
# The band-projected pair potential (Nagai-Nakamura, JPSJ 85, 074707 (2016) Eq. 43)
#
#     phi(k_F) = sum_ab conj(u_a(k_F)) Delta_orb,ab(k_F) conj(u_b(-k_F))
#
# pairs |k> with a PARTNER state, and the partner must be specified in the same
# gauge as u(k).  Feeding an independently diagonalized u(-k) does NOT do that:
# LAPACK fixes no phase convention, so u(k) -> e^{i theta(k)} u(k) sends
# phi -> e^{-i(theta(k)+theta(-k))} phi and the result is gauge DEPENDENT --
# a silently wrong, k-dependent phase (and, for real eigenvectors, sign) on phi.
#
# The partner is therefore always constructed from u(k) alone:
#
#   'trs'  spinless / single pseudo-spin sector with time-reversal symmetry,
#          H(-k) = conj(H(k)) (guaranteed when the hoppings h(R) are real).
#          NOTE this is a symmetry of the NORMAL-STATE H(k), not of Delta: a
#          CHIRAL order parameter breaks time reversal in the superconducting
#          state but leaves H untouched, so the gauge fixing stays valid and the
#          projected phi correctly inherits the winding of Delta_orb(k).  What
#          does invalidate it is a time-reversal-broken NORMAL state (a
#          ferromagnetic superconductor, a Zeeman term inside H, Haldane-type
#          complex hoppings) -- there the time-reversed state is not an
#          eigenstate at -k at all; check_pair_partner detects exactly that.
#          Then u(-k) = conj(u(k)) and conj(u(-k)) = u(k), so
#          phi = u^dag Delta u -- manifestly gauge invariant, real for a
#          Hermitian Delta, and equal to 1 for Delta = identity.  DEFAULT.
#   'soc'  spinful basis: the partner is the time-reversed state
#          T u = (i sigma_y) conj(u), whose conjugate has components
#          conj(T u)_{a,up} = +u_{a,dn},  conj(T u)_{a,dn} = -u_{a,up}.
#          Needs ``spin_map`` (see spin_pair_map) to know the spin partner of
#          each basis index.  Also gauge invariant.
#   'diag' legacy: the independently diagonalized u(-k).  Gauge DEPENDENT; kept
#          only for comparison / reproducing older numbers.  It happens to agree
#          with 'trs' exactly whenever LAPACK returns u(-k) = conj(u(k)), which
#          holds for real h(R) -- verify with ``check_trs_gauge``.

PAIR_GAUGE = 'trs'          # module default, see set_pair_gauge
PAIR_SPIN_MAP = None        # (jidx, sgn) for PAIR_GAUGE='soc'


def set_pair_gauge(gauge: str, spin_map=None):
    """
    @fn set_pair_gauge
    @brief Set the module-wide pair-partner gauge used by the band-projection
    routines (project_gap_to_band, gap_color_3d, build_fs).
    @param   gauge: 'trs' (default), 'soc' (needs spin_map) or 'diag' (legacy)
    @param spin_map: (jidx, sgn) arrays for 'soc', from ``spin_pair_map``
    """
    global PAIR_GAUGE, PAIR_SPIN_MAP
    if gauge not in ('trs', 'soc', 'diag'):
        raise ValueError(f"unknown pair gauge {gauge!r} (use 'trs', 'soc' or 'diag')")
    if gauge == 'soc' and spin_map is None and PAIR_SPIN_MAP is None:
        raise ValueError("pair gauge 'soc' needs spin_map=(jidx,sgn); see spin_pair_map")
    PAIR_GAUGE = gauge
    if spin_map is not None:
        PAIR_SPIN_MAP = (np.asarray(spin_map[0], dtype=np.int64),
                         np.asarray(spin_map[1], dtype=np.float64))
    return PAIR_GAUGE


def spin_pair_map(Nbasis: int, order: str = 'block'):
    """
    @fn spin_pair_map
    @brief Build the (jidx, sgn) spin-partner map required by the 'soc' pair gauge:
    conj(T u)_i = sgn[i] * u[jidx[i]], i.e. up -> +down, down -> -up.
    @param Nbasis: total number of basis functions (orbitals x spin, must be even)
    @param  order: 'block' = [orb_1..orb_N (up), orb_1..orb_N (down)] (pyrpa SOC
                   convention); 'interleave' = [orb_1 up, orb_1 dn, orb_2 up, ...]
    @return (jidx, sgn)
    """
    if Nbasis % 2:
        raise ValueError(f"spin_pair_map: Nbasis={Nbasis} is odd (not a spinful basis)")
    j = np.arange(Nbasis, dtype=np.int64)
    sgn = np.ones(Nbasis, dtype=np.float64)
    if order == 'block':
        N = Nbasis // 2
        j = np.concatenate([np.arange(N, Nbasis), np.arange(0, N)]).astype(np.int64)
        sgn[N:] = -1.0
    elif order == 'interleave':
        j = (j ^ 1).astype(np.int64)                 # 0<->1, 2<->3, ...
        sgn[1::2] = -1.0
    else:
        raise ValueError(f"unknown spin ordering {order!r} ('block' or 'interleave')")
    return j, sgn


def pair_partner_conj(ev: np.ndarray, evm: np.ndarray = None, gauge: str = None,
                      spin_map=None) -> np.ndarray:
    """
    @fn pair_partner_conj
    @brief The array that plays the role of ``conj(u_b(-k_F))`` in the band
    projection, built in the requested gauge (see the block comment above).
    @param  ev: band eigenvectors u(k_F) [Nfs, Nbasis]
    @param evm: independently diagonalized u(-k_F); only used by gauge='diag'
    @param gauge,spin_map: override the module defaults (set_pair_gauge)
    @return [Nfs, Nbasis] array to contract as conj(u(-k))
    """
    gauge = PAIR_GAUGE if gauge is None else gauge
    if gauge == 'trs':
        return ev                                    # conj(u(-k)) = conj(conj(u(k)))
    if gauge == 'soc':
        sm = PAIR_SPIN_MAP if spin_map is None else spin_map
        if sm is None:
            raise ValueError("pair gauge 'soc' needs spin_map=(jidx,sgn); see spin_pair_map")
        j, sgn = np.asarray(sm[0], dtype=np.int64), np.asarray(sm[1])
        if len(j) != ev.shape[1]:
            raise ValueError(f"spin_map length {len(j)} != basis size {ev.shape[1]}")
        return ev[:, j] * sgn[None, :]
    if gauge == 'diag':
        if evm is None:
            raise ValueError("pair gauge 'diag' needs the u(-k) eigenvectors (evm)")
        return np.conj(evm)
    raise ValueError(f"unknown pair gauge {gauge!r}")


def check_trs_gauge(ev: np.ndarray, evm: np.ndarray, tol: float = 1.0e-8,
                    warn: bool = True) -> float:
    """
    @fn check_trs_gauge
    @brief Diagnostic: max |u(-k) - conj(u(k))| over the Fermi surface.  It is 0 to
    machine precision for real hoppings (H(-k)=conj(H(k)) and LAPACK is conjugation
    consistent), in which case the legacy 'diag' route coincides with 'trs'.  A
    finite value means the two differ and only the gauge-fixed routes are correct.
    @return the max deviation (and warns above ``tol`` when ``warn``)
    """
    dev = float(np.abs(evm - np.conj(ev)).max()) if len(ev) else 0.0
    if warn and dev > tol:
        print(f"Warning: u(-k) != conj(u(k)) (max dev {dev:.3e}); the pair partner is "
              f"gauge fixed by PAIR_GAUGE={PAIR_GAUGE!r} -- 'diag' would be gauge dependent",
              flush=True)
    return dev


def check_pair_partner(ev: np.ndarray, evm: np.ndarray, gauge: str = None, spin_map=None,
                       tol: float = 1.0e-6, warn: bool = True) -> float:
    """
    @fn check_pair_partner
    @brief The decisive validity check of the gauge fixing: the time-reversed state
    T|k,n> must BE the band n state at -k, i.e. |<u(-k) | T u(k)>| = 1.  Unlike
    check_trs_gauge this is gauge invariant (insensitive to the arbitrary phase of
    u(-k)) and correct for both 'trs' and 'soc'.

    It is 1 whenever the NORMAL STATE has time-reversal symmetry -- including when
    the order parameter is chiral, since Delta plays no part in it.  It drops below
    1 only when H(k) itself breaks time reversal (ferromagnetic superconductor, a
    Zeeman/exchange term inside H, Haldane-type complex hoppings); no gauge-fixed
    scalar pair partner exists in that case and only |phi| stays meaningful, so the
    projection is reported as unreliable rather than silently used.
    @return the minimum overlap over the Fermi surface (warns below 1 - tol)
    """
    if not len(ev):
        return 1.0
    part = pair_partner_conj(ev, evm, gauge, spin_map)        # conj(T u)
    ov = float(np.abs(np.einsum('ia,ia->i', evm, part)).min())
    if warn and ov < 1.0 - tol:
        print(f"Warning: |<u(-k)|T u(k)>| drops to {ov:.4f} (should be 1): the NORMAL "
              f"state breaks time-reversal symmetry, so the band-projected pair "
              f"potential has no gauge-fixed phase -- only |phi| is meaningful "
              f"(a chiral Delta alone does NOT cause this)", flush=True)
    return ov


def _eval_gap_orbital(gap_orbital, kfrac: np.ndarray) -> np.ndarray:
    """Evaluate an orbital-basis pair potential on a batch of k-points [Nk,3].
    Accepts a constant matrix, a callable that already handles a batch (returning
    [Nk,N,N]), or a per-k callable (looped).  @return [Nk, N, N] complex."""
    kfrac = np.atleast_2d(np.asarray(kfrac, dtype=np.float64))
    if not callable(gap_orbital):
        D = np.asarray(gap_orbital, dtype=np.complex128)
        return np.broadcast_to(D, (len(kfrac),) + D.shape)
    try:                                             # batched callable (fast path)
        out = np.asarray(gap_orbital(kfrac), dtype=np.complex128)
        if out.ndim == 3 and out.shape[0] == len(kfrac):
            return out
    except Exception:
        pass
    return np.array([np.asarray(gap_orbital(k), dtype=np.complex128) for k in kfrac],
                    dtype=np.complex128)


def project_gap(ev: np.ndarray, gap_orbital, kfrac: np.ndarray, evm: np.ndarray = None,
                gauge: str = None, spin_map=None) -> np.ndarray:
    """
    @fn project_gap
    @brief Band-diagonal projection phi_i = sum_ab conj(u_a(k_i)) Delta_orb,ab(k_i) P_b(k_i)
    with the pair partner P = pair_partner_conj(...).  The single place where the
    Nagai-Nakamura Eq. 43 contraction is performed; shared by the FS-dict, 3D-mesh
    and 3D-plot routes so they cannot drift apart.
    @return phi [Nfs] complex (unnormalized)
    """
    part = pair_partner_conj(ev, evm, gauge, spin_map)
    kfrac = np.atleast_2d(np.asarray(kfrac, dtype=np.float64))
    phi = np.empty(len(ev), dtype=np.complex128)
    # Chunked: a 3D Fermi surface easily has 10^5 points, and materializing
    # Delta_orb(k) for all of them at once is Nfs*Norb^2 complex numbers.
    for i0 in range(0, len(ev), 8192):
        sl = slice(i0, i0 + 8192)
        D = _eval_gap_orbital(gap_orbital, kfrac[sl])            # [chunk, N, N]
        phi[sl] = np.einsum('ia,iab,ib->i', np.conj(ev[sl]), D, part[sl])
    return phi



def _fs_average(values: np.ndarray, w: np.ndarray) -> np.ndarray:
    """DOS-weighted Fermi-surface average over the FS-point axis (axis 0).
    @param values: array [Nfs, ...] of the quantity to average
    @param      w: FS weights [Nfs]
    @return       averaged array [...] (axis 0 contracted)
    """
    return np.tensordot(w, values, axes=(0, 0)) / w.sum()


def build_fs(eig: np.ndarray, klist: np.ndarray, mu: float, gap_sym: int,
             width: float, w_cut: float = 1.0e-4, uni: np.ndarray = None,
             gap_orbital=None, delta0=None, gauge: str = None, spin_map=None,
             sw_band: bool = False):
    """
    @fn build_fs
    @brief Select Fermi-surface points from a band mesh and build their DOS
    weights and (normalized) pairing form factor.

    A broadened delta function delta(eps_{k,n} - mu) (Gaussian of width
    ``width``) provides the DOS weight, so no explicit FS triangulation is
    needed and all bands are handled uniformly.  Points with negligible weight
    are discarded because the quasiclassical propagators live on the FS.
    Because the mesh is the FULL 3D k-grid, every sheet and every k_z is kept:
    this is the 3D Fermi-surface route.

    The form factor comes from one of three sources, in decreasing priority:

      * ``gap_orbital`` -- an orbital-basis pair potential (constant matrix or
        callable kfrac->NxN, e.g. an RPA/FLEX gap loaded by
        ``gap_orbital_from_wannier``) projected BAND-RESOLVED onto the Fermi
        surface (Nagai-Nakamura Eq. 43, ``project_gap``, gauge fixed).  This is
        the only route that can produce an ACCIDENTAL node: the k_z dependence
        and the multi-sheet sign structure come from the orbital character of
        each band, not from a chosen harmonic.  Needs ``uni``.
      * ``gap_sym`` -- an analytic lattice harmonic (``gap_symms``), identical on
        every band, optionally scaled per band by ``delta0`` (phenomenological
        multiband s+- signs / amplitude ratios).
      * ``gap_sym`` alone -- the historical behaviour.

    @param    eig: band energies on the mesh [Nk, Norb]
    @param  klist: k-points in fractional coordinates [Nk, 3]
    @param     mu: chemical potential [eV]
    @param gap_sym: gap-symmetry index passed to ``gap_symms``
    @param  width: Gaussian broadening of the FS delta function [eV]
    @param  w_cut: keep points with weight > w_cut * max(weight)
    @param    uni: band eigenvectors on the same mesh [Nk, Nband, Norb]
                   (``get_emesh(..., sw_uni=True)``); required by gap_orbital
    @param gap_orbital: orbital-basis pair potential to project (see above)
    @param delta0: per-band gap amplitudes/signs indexed by band; applied to the
                   harmonic route (and, if given, on top of the projection)
    @param gauge,spin_map: pair-partner gauge (see set_pair_gauge); the mesh route
                   has no independently diagonalized u(-k), so 'diag' is rejected
    @param sw_band: also return per-point provenance (band index, mesh index, k_z)
    @return (wf, phif) or (wf, phif, info): FS weights [Nfs], normalized form
            factor [Nfs] and optionally info = dict(band, ik, kz)
    """
    # Gaussian-broadened delta(eps - mu); the overall constant cancels in averages.
    de = eig - mu
    w = np.exp(-0.5 * (de / width) ** 2)            # [Nk, Norb]
    mask = w > w_cut * w.max()
    ik, ib = np.nonzero(mask)                       # mesh index / band index per FS point
    wf = w[mask]
    if gap_orbital is not None:
        if uni is None:
            raise ValueError("build_fs: gap_orbital needs the band eigenvectors "
                             "(get_emesh(..., sw_uni=True) -> uni)")
        if (PAIR_GAUGE if gauge is None else gauge) == 'diag':
            raise ValueError("build_fs: pair gauge 'diag' needs u(-k) from a second "
                             "diagonalization and is not available on the mesh route; "
                             "use 'trs' (real hoppings) or 'soc'")
        ev = uni[ik, ib, :]                          # u_band(k) on the FS  [Nfs, Norb]
        nrm = np.sqrt((np.abs(ev) ** 2).sum(axis=1))  # (MLO basis is non-orthogonal)
        ev = ev / np.where(nrm > 0, nrm, 1.0)[:, None]
        phif = project_gap(ev, gap_orbital, klist[ik], None, gauge, spin_map)
    else:
        phif = gap_symms(klist, 1, gap_sym)[0][ik].astype(np.complex128)   # k only
    if delta0 is not None:                           # per-band ratio / sign
        phif = phif * np.asarray(delta0, dtype=np.float64)[ib]
    # normalize so that <|phi|^2>_FS = 1 -> lambda is the dimensionless coupling
    # (|phi|^2, not phi^2, so complex/chiral form factors normalize correctly too)
    norm = np.sqrt(_fs_average(np.abs(phif) ** 2, wf))
    if norm > 0:
        phif = phif / norm
    if sw_band:
        return wf, phif, dict(band=ib, ik=ik, kz=klist[ik, 2])
    return wf, phif


def report_fs_gap(wf: np.ndarray, phif: np.ndarray, info: dict, kb: float = 1.0):
    """
    @fn report_fs_gap
    @brief Print the sheet-resolved structure of a Fermi-surface form factor: the
    weight and |phi| of every band, whether Re[phi] changes sign ON one sheet
    (accidental node) or only BETWEEN sheets (s+- type), and -- for a 3D FS -- how
    phi varies with |k_z|, which is where a HORIZONTAL line node shows up.
    @param wf,phif,info: the three outputs of build_fs(..., sw_band=True)
    """
    band, kz = info['band'], np.abs(info['kz'])
    wsum = wf.sum()
    # A chiral phi has no meaningful sign -- its phase winds -- so a node is where
    # |phi| vanishes, not where Re[phi] crosses zero (which would flag phantom nodes).
    chiral = np.abs(phif.imag).max() > 1.0e-10
    val = np.abs(phif) if chiral else phif.real
    lab = "|phi|" if chiral else "Re[phi]"
    small = 0.05 * np.abs(phif).mean()
    nodal = (lambda v: v.min() < small) if chiral else (lambda v: v.min() < 0.0 < v.max())
    print(f"  sheet   weight   <|phi|>    {lab} min    max   node"
          f"{'   (chiral: nodes are |phi|=0)' if chiral else ''}", flush=True)
    for b in np.unique(band):
        m = band == b
        print(f"  band {b:<3d} {wf[m].sum()/wsum:7.3f} {np.abs(phif[m]).mean():9.4f} "
              f"{val[m].min():+10.4f} {val[m].max():+8.4f}   "
              f"{'ON-SHEET' if nodal(val[m]) else '-':>8s}", flush=True)
    if len(np.unique(np.round(kz, 10))) > 2:          # 3D FS: k_z resolved profile
        edges = np.linspace(0.0, 0.5, 6)
        print(f"  {lab} vs |k_z| (0 -> 1/2), per sheet:", flush=True)
        for b in np.unique(band):
            m = band == b
            prof = [val[m & (kz >= a) & (kz < c)].mean()
                    if (m & (kz >= a) & (kz < c)).sum() > 3 else np.nan
                    for a, c in zip(edges[:-1], edges[1:])]
            row = " ".join(f"{v:+7.3f}" if np.isfinite(v) else "   --  " for v in prof)
            good = [v for v in prof if np.isfinite(v)]
            if chiral:      # |phi| dips to zero somewhere inside a k_z bin (binning the
                            # MEAN would wash a thin node plane out, so test the minima)
                hit = any(val[m & (kz >= a) & (kz < c)].min() < small
                          for a, c in zip(edges[:-1], edges[1:])
                          if (m & (kz >= a) & (kz < c)).sum() > 3)
            else:
                hit = good and min(good) < 0.0 < max(good)
            print(f"    band {b:<3d} {row}"
                  f"{' <- HORIZONTAL node along k_z' if hit else ''}", flush=True)


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
    phif_w_avg = _fs_average(np.conj(phif)[:, None] * f, wf)   # <phi* f> (conj for chiral)
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
        damp_init = BCS_RATIO * temp
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
                     h_list=None, fs_width: float = 5.0e-3, kb: float = 1.0,
                     gap_orbital=None, delta0=None, gauge: str = None, spin_map=None):
    """
    @fn calc_pauli_limit
    @brief Sweep the Zeeman (Maki) field h and report the singlet gap Delta(h) (the
    metastable SC branch) -- Pauli paramagnetic depairing.  At low T the SC solution
    survives up to the spinodal h* ~ Delta0, while the thermodynamic first-order
    (Chandrasekhar-Clogston) transition is at h_P = Delta0/sqrt(2).  Also writes the
    Zeeman-split bulk DOS at h = 0.5*Delta0 to 'pauli_dos.dat'.
    @param Nx,Ny,Nz,ham_r,S_r,rvec,avec,mu: Fermi-surface inputs (as in calc_eilenberger)
    @param h_list: list of Zeeman energies h [eV] (default fractions of Delta0)
    @param gap_orbital,delta0,gauge,spin_map: gap specification on the 3D Fermi surface,
           exactly as in calc_eilenberger (band-projected orbital gap / per-sheet signs)
    """
    omega = matsubara(temp, wc)
    Nk, klist, eig, uni, kweight = get_emesh(Nx, Ny, Nz, ham_r, S_r, rvec, avec, sw_uni=True)
    wf, phif, fsinfo = build_fs(eig, klist, mu, gap_sym, fs_width, uni=uni,
                                gap_orbital=gap_orbital, delta0=delta0,
                                gauge=gauge, spin_map=spin_map, sw_band=True)
    report_fs_gap(wf, phif, fsinfo, kb)
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
    ndos = dos_zeeman(wl, D0, wf, phif, 0.5 * D0, 0.03 * D0)
    try:
        with open('pauli_dos.dat', 'w') as f2:
            f2.write("# w/Delta0   N(w)/N0  (Zeeman split, h=0.5 Delta0)\n")
            for w, n in zip(wl, ndos):
                f2.write(f"{w/D0:10.4f} {n:12.5e}\n")
    except IOError as e:
        print(f"Error writing pauli_dos.dat: {e}", flush=True)


def calc_eilenberger(Nx: int, Ny: int, Nz: int, wc: float, ham_r, S_r, rvec, avec,
                     mu: float, temp: float, gap_sym: int, coupling: float,
                     imp_gamma: float = 0.0, imp_c: float = 1.0e8,
                     fs_width: float = 5.0e-3, kb: float = 1.0, method: str = 'normalization',
                     sw_find_tc: bool = False, sw_imp_sweep: bool = False,
                     imp_sweep: np.ndarray = None, gap_orbital=None, delta0=None,
                     gauge: str = None, spin_map=None):
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
    @param gap_orbital: orbital-basis pair potential (constant matrix or callable
                       kfrac->NxN, e.g. an RPA/FLEX gap from gap_orbital_from_wannier).
                       It is projected BAND-RESOLVED onto the 3D Fermi surface, so the
                       k_z dependence and the multi-sheet signs come from the orbital
                       character -- the only route to an ACCIDENTAL (e.g. horizontal)
                       node, which no gap_symms harmonic can represent.
    @param    delta0: per-band gap amplitudes/signs (phenomenological multiband s+-)
    @param gauge,spin_map: pair-partner gauge for the projection (set_pair_gauge)
    """
    print("calculate homogeneous quasiclassical Eilenberger equation", flush=True)
    Nk, klist, eig, uni, kweight = get_emesh(Nx, Ny, Nz, ham_r, S_r, rvec, avec, sw_uni=True)
    wf, phif, fsinfo = build_fs(eig, klist, mu, gap_sym, fs_width, uni=uni,
                                gap_orbital=gap_orbital, delta0=delta0,
                                gauge=gauge, spin_map=spin_map, sw_band=True)
    print(f"Fermi-surface points kept: {len(wf)} (of {eig.size})", flush=True)
    if len(wf) == 0:
        print("Error: no Fermi-surface points found; check mu / fs_width", flush=True)
        return
    if gap_orbital is not None:
        print(f"pairing form factor: band projection of an orbital gap "
              f"(pair gauge '{PAIR_GAUGE if gauge is None else gauge}')", flush=True)
    report_fs_gap(wf, phif, fsinfo, kb)
    omega = matsubara(temp, wc)
    print(f"Matsubara cutoff wc = {wc:.4e} eV ({len(omega)} freqs at this T)", flush=True)
    print(f"pairing coupling lambda = {coupling:.4f}, gap_sym = "
          f"{'projected orbital gap' if gap_orbital is not None else gap_sym}", flush=True)
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
    omega = matsubara(temp, wc)
    beta = np.linspace(0.0, 2.0 * np.pi, Nbeta, endpoint=False)
    phi = form_factor(beta, gap_sym)
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
    gap_sym = _gap_sym_str(gap_sym)                   # accept the global int index too
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
    """Return (eps(kx,ky,kz), grad eps -> (vx,vy,vz), default mu, radial search rmax).

    The in-plane kinds ('iso', 'ellipse', 'tb') ignore k_z and have v_z = 0, so
    stacking k_z slices reproduces the exact cylinder -- which makes them the
    control case for the genuinely 3D kinds:
      'cyl'      corrugated cylinder -2t(cos kx + cos ky) - 2 tz cos kz, params=(t,tz).
                 The canonical quasi-2D -> 3D testbed: tz/t tunes the warping, and a
                 k_z-dependent gap (dxz/dyz/dxz+idyz, horizontal nodes) lives on it.
      'sphere'   isotropic 3D (kx^2+ky^2+kz^2)/2 -- a CLOSED sheet, so slices beyond
                 the pole carry no Fermi surface and are dropped.
      'spheroid' anisotropic 3D, params=(mx,my,mz).
    """
    if kind == 'iso':
        return (lambda kx, ky, kz: 0.5 * (kx ** 2 + ky ** 2),
                lambda kx, ky, kz: (kx, ky, np.zeros_like(kx)), 1.0, 10.0)
    if kind == 'ellipse':
        mx, my = (params or (1.0, 0.5))
        return (lambda kx, ky, kz: 0.5 * (kx ** 2 / mx + ky ** 2 / my),
                lambda kx, ky, kz: (kx / mx, ky / my, np.zeros_like(kx)), 1.0, 10.0)
    if kind == 'tb':
        t = (params or 1.0)
        t = t if np.isscalar(t) else t[0]
        return (lambda kx, ky, kz: -2.0 * t * (np.cos(kx) + np.cos(ky)),
                lambda kx, ky, kz: (2.0 * t * np.sin(kx), 2.0 * t * np.sin(ky), np.zeros_like(kx)),
                -2.0 * t * 0.6, np.pi - 1e-6)
    if kind == 'cyl':
        t, tz = (params or (1.0, 0.2))[:2]
        return (lambda kx, ky, kz: -2.0 * t * (np.cos(kx) + np.cos(ky)) - 2.0 * tz * np.cos(kz),
                lambda kx, ky, kz: (2.0 * t * np.sin(kx), 2.0 * t * np.sin(ky),
                                    2.0 * tz * np.sin(kz) * np.ones_like(kx)),
                -2.0 * t * 0.6, np.pi - 1e-6)
    if kind == 'sphere':
        return (lambda kx, ky, kz: 0.5 * (kx ** 2 + ky ** 2 + kz ** 2),
                lambda kx, ky, kz: (kx, ky, kz * np.ones_like(kx)), 1.0, 10.0)
    if kind == 'spheroid':
        mx, my, mz = (params or (1.0, 1.0, 0.5))[:3]
        return (lambda kx, ky, kz: 0.5 * (kx ** 2 / mx + ky ** 2 / my + kz ** 2 / mz),
                lambda kx, ky, kz: (kx / mx, ky / my, kz / mz * np.ones_like(kx)), 1.0, 10.0)
    raise ValueError(f"unknown FS kind: {kind}")


def _radial_kf(eps, mu, th, kzv, rmax):
    """Radial Fermi wavenumber k_F(theta) of a convex sheet on one k_z slice, or None
    if the slice carries no Fermi surface (beyond the pole of a closed sheet)."""
    kF = np.empty(len(th))
    for i, t in enumerate(th):
        c, s = np.cos(t), np.sin(t)
        f = lambda r: eps(r * c, r * s, kzv) - mu
        lo, hi = 1e-6, rmax
        try:
            kF[i] = brentq(f, lo, hi)
        except ValueError:
            rs = np.linspace(lo, hi, 200)
            fv = np.array([f(r) for r in rs])
            sgn = np.where(np.diff(np.sign(fv)) != 0)[0]
            if not len(sgn):
                return None                            # no crossing on this slice
            kF[i] = brentq(f, rs[sgn[0]], rs[sgn[0] + 1])
    return kF


def build_model_fs(kind: str = 'iso', Nth: int = 360, mu: float = None, params=None,
                   nkz: int = 1, kz_max: float = np.pi):
    """
    @fn build_model_fs
    @brief Build a model Fermi surface (radial parametrization, convex single sheet)
    with Fermi velocities and DOS weights -- the analytic counterpart of
    build_wannier_fs, with the same k_z stacking and the same conventions.

    @param kind: in-plane 'iso', 'ellipse' (params=(mx,my)), 'tb' (params=t); 3D
                 'cyl' (corrugated cylinder, params=(t,tz)), 'sphere',
                 'spheroid' (params=(mx,my,mz)).  See _disp.
    @param Nth: FS points per k_z slice (angular samples)
    @param  mu: chemical potential (default per dispersion)
    @param nkz: number of k_z slices (1 = the historical single-slice, quasi-2D FS).
                Slices are uniform on k_z in [-kz_max, kz_max); those carrying no
                Fermi surface (a closed sheet past its pole) are dropped.
    @param kz_max: half width of the k_z range (pi for a lattice kind)
    @return dict: th,kx,ky,kz,vx,vy,vz,vabs,vabs3,vhx,vhy,nf,nkz [,kf]
            (nf normalized to sum 1)

    @note  ``vabs`` is the IN-PLANE speed, which is the correct chord velocity AND
           the correct DOS measure in 3D as well -- see the note on build_wannier_fs.
           The lattice kinds ('tb', 'cyl') also carry ``kf`` = k/(2 pi), the fractional
           coordinate, so the gap_symms lattice harmonics (including the k_z dependent
           4 dxz / 5 dyz / 7 dxz+idyz) work directly on the model FS.
    """
    eps, grad, mu0, rmax = _disp(kind, params)
    if mu is None:
        mu = mu0
    lattice = kind in ('tb', 'cyl')                    # Cartesian k with a = 1
    kzs = (np.array([0.0]) if nkz <= 1 else
           np.linspace(-kz_max, kz_max, int(nkz), endpoint=False))
    th = np.linspace(0.0, 2.0 * np.pi, Nth, endpoint=False)
    dth = th[1] - th[0]
    ths, kxs, kys, kzl, vxs, vys, vzs, nfs = ([] for _ in range(8))
    for kzv in kzs:
        kF = _radial_kf(eps, mu, th, kzv, rmax)
        if kF is None:
            continue                                   # closed sheet: past the pole
        kx, ky = kF * np.cos(th), kF * np.sin(th)
        vx, vy, vz = grad(kx, ky, kzv)
        vab = np.sqrt(vx ** 2 + vy ** 2)
        dkF = (np.roll(kF, -1) - np.roll(kF, 1)) / (2.0 * dth)
        dl = np.sqrt(kF ** 2 + dkF ** 2) * dth         # arc length element
        ths.append(th); kxs.append(kx); kys.append(ky)
        kzl.append(np.full(Nth, kzv))
        vxs.append(vx); vys.append(vy); vzs.append(np.broadcast_to(vz, kx.shape).copy())
        nfs.append(dl / np.maximum(vab, 1e-12))        # dk_z uniform -> cancels below
    if not kxs:
        raise ValueError(f"model FS '{kind}': no Fermi surface at mu={mu}")
    kx = np.concatenate(kxs); ky = np.concatenate(kys); kz = np.concatenate(kzl)
    vx = np.concatenate(vxs); vy = np.concatenate(vys); vz = np.concatenate(vzs)
    vabs = np.sqrt(vx ** 2 + vy ** 2)
    nf = np.concatenate(nfs)
    nf = nf / nf.sum()
    fs = dict(th=np.concatenate(ths), kx=kx, ky=ky, kz=kz, vx=vx, vy=vy, vz=vz,
              vabs=vabs, vabs3=np.sqrt(vx ** 2 + vy ** 2 + vz ** 2),
              vhx=vx / vabs, vhy=vy / vabs, nf=nf, nkz=len(kxs))
    if lattice:                                        # fractional k -> gap_symms harmonics
        fs['kf'] = np.stack([kx, ky, kz], axis=1) / (2.0 * np.pi)
    return fs


def fs_hvf(fs, hvf: float = 1.0) -> float:
    """
    @fn fs_hvf
    @brief The representative Fermi speed that sets the coherence length
    xi = hvf/(pi*Dbulk), i.e. the nf-weighted mean IN-PLANE speed <|v_par|> of a
    Fermi surface (``hvf`` itself when there is no FS, i.e. the isotropic cylinder).

    Single definition shared by every trajectory solver: xi fixes the physical size
    of the computational domain, so taking it as 1 while the chords run at the real
    |v_F| silently shrinks the box and the order parameter never heals to its bulk
    value.  Note this is only the GEOMETRY scale -- each chord still integrates with
    its own velocity.
    """
    if fs is None:
        return float(hvf)
    return float((np.asarray(fs['nf']) * np.asarray(fs['vabs'])).sum())


def fs_field_frame(fs, bdir=(0.0, 0.0, 1.0), gap_sym=None, aniso=True, vmin_frac=1.0e-3,
                   verbose=False):
    """
    @fn fs_field_frame
    @brief Trajectory set for vortex lines along an ARBITRARY field direction.

    The vortex lines run along B, so the order parameter varies only in the plane
    PERPENDICULAR to B and the transport term is v_perp . grad_perp: the problem stays
    TWO DIMENSIONAL, just in a different plane.  Everything the 2D chord solvers need is
    therefore still available -- a direction in the plane, a speed along it, a weight and
    a form factor -- and this returns exactly that, in an orthonormal frame (e1, e2, n)
    with n = B_hat.  For B along c the frame is (x, y) and the historical set is
    reproduced identically.

    Anisotropy.  With B in the plane of a quasi-2D Fermi surface the two directions of the
    vortex plane are wildly different: |v_e1| ~ v_F but |v_e2| ~ v_z << v_F, so the core is
    elliptical with xi_1/xi_2 = <|v_1|>/<|v_2|>, and a square grid would be mostly wasted.
    With aniso=True each axis is rescaled by its own rms velocity.  That is an AFFINE map,
    so trajectories stay straight and the square grid plus rotate-interpolate machinery
    apply verbatim; the vortex is a winding-1 point in the SCALED plane, which is exactly
    the anisotropic-GL solution.  Physical coordinates come back as
    r_1 = (s1/sbar) u,  r_2 = (s2/sbar) w  with the returned ``scale`` = (s1/sbar, s2/sbar).

    The Fermi-surface weight ``nf`` is untouched: it is the 3D DOS measure
    dS/|v_F| = dk_z dl/|v_par| and knows nothing about the field direction.  Only the
    chord VELOCITY changes, from |v_xy| to the in-plane part of v in the new frame.

    @param      fs: FS dict; needs vz for anything but B || c (build_wannier_fs(nkz>1) or
                    a 3D build_model_fs kind).  A 2D FS is treated as vz = 0, which is the
                    exact cylinder limit -- and then B in-plane has NO orbital effect at
                    all (v_perp -> |v_z| = 0), so it is rejected rather than silently
                    returning a zero-velocity trajectory set.
    @param    bdir: field direction (need not be normalized)
    @param   aniso: rescale the two axes by their rms velocities (see above)
    @param vmin_frac: drop trajectories whose perpendicular speed is below this fraction of
           the mean, and renormalize the weights.  A Fermi-surface point whose velocity is
           PARALLEL to B does not move in the vortex plane at all, so its chord step
           ds/|v_perp| is infinite -- the matrix exponential is then handed a non-finite
           generator and LAPACK fails deep inside with a stream of MKL parameter errors.
           Such points exist for any field direction that is not a symmetry axis of the FS
           (on a corrugated cylinder with B in-plane they are the points with k_y=k_z=0),
           but they are isolated, i.e. of measure zero, so dropping them is the right
           treatment rather than a fudge.
    @return dict(e1, e2, n, dirs, vabs, nf, phi, scale, hvf_eff, aniso_ratio); phi is None
            when neither gap_sym nor a baked-in fs['phi'] is available
    """
    n = np.asarray(bdir, dtype=np.float64)
    nn = np.sqrt((n ** 2).sum())
    if nn <= 0:
        raise ValueError("fs_field_frame: bdir must be a non-zero vector")
    n = n / nn
    # e1, e2 completing a right-handed frame.  e2 is taken as the most c-LIKE direction
    # available (the component of z_hat perpendicular to B), so that as B rotates the two
    # axes keep their meaning: e2 is the out-of-plane one and aniso_ratio is consistently
    # xi_(in plane)/xi_c.  Choosing e2 off whichever coordinate axis is least aligned with
    # B instead makes the pair SWAP as B turns -- for B rotating in the ab plane the ratio
    # jumped 0.667 -> 1.525 between phi = 0 and 15 degrees (the same physics with the two
    # labels exchanged), which is useless for an in-plane field-angle scan.  For B || c the
    # construction is degenerate and falls back to exactly (x, y).
    zc = np.array([0.0, 0.0, 1.0])
    e2 = zc - (zc @ n) * n
    if np.sqrt((e2 ** 2).sum()) < 1.0e-8:              # B || c: keep the historical frame
        e1, e2 = np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])
    else:
        e2 /= np.sqrt((e2 ** 2).sum())
        e1 = np.cross(e2, n)
    vz = fs.get('vz')
    v = np.stack([np.asarray(fs['vx']), np.asarray(fs['vy']),
                  np.zeros_like(fs['vx']) if vz is None else np.asarray(vz)], axis=1)
    v1, v2 = v @ e1, v @ e2
    nf = np.asarray(fs['nf'], dtype=np.float64)
    s1 = np.sqrt((nf * v1 ** 2).sum())
    s2 = np.sqrt((nf * v2 ** 2).sum())
    if min(s1, s2) <= 1.0e-12 * max(s1, s2, 1.0e-30):
        raise ValueError(
            f"fs_field_frame: the Fermi velocity has no component along e2 for this field "
            f"direction (<v_1^2>={s1**2:.3e}, <v_2^2>={s2**2:.3e}).  A strictly 2D Fermi "
            "surface has no orbital response to an in-plane field; use a 3D FS "
            "(build_wannier_fs(nkz>1) or a 3D build_model_fs kind).")
    if aniso:
        sbar = np.sqrt(s1 * s2)
        c1, c2 = sbar / s1, sbar / s2                  # scaled velocity = (c1 v1, c2 v2)
    else:
        sbar, c1, c2 = 1.0, 1.0, 1.0
    w1, w2 = c1 * v1, c2 * v2
    vperp = np.hypot(w1, w2)
    keep = vperp > vmin_frac * float((nf * vperp).sum())   # drop v || B (see vmin_frac)
    if not keep.all():
        w1, w2, vperp = w1[keep], w2[keep], vperp[keep]
        nf = nf[keep] / nf[keep].sum()
        fs = {k: (v[keep] if isinstance(v, np.ndarray) and len(v) == len(keep) else v)
              for k, v in fs.items()}
    # phi is optional: the d-vector solvers carry their own per-channel form factors and
    # pass no gap_sym, so asking fs_form_factor for one would fail on a bare FS dict
    have_phi = (gap_sym is not None) or ('phi' in fs)
    out = dict(e1=e1, e2=e2, n=n, dirs=np.arctan2(w2, w1), vabs=vperp, nf=nf,
               phi=fs_form_factor(fs, gap_sym).astype(np.complex128) if have_phi else None,
               scale=(1.0 / c1, 1.0 / c2), aniso_ratio=float(s1 / s2), keep=keep)
    out['hvf_eff'] = float((nf * out['vabs']).sum())
    if verbose:
        print(f"  field frame: B_hat=({n[0]:+.3f},{n[1]:+.3f},{n[2]:+.3f})  "
              f"e1=({e1[0]:+.3f},{e1[1]:+.3f},{e1[2]:+.3f})  "
              f"e2=({e2[0]:+.3f},{e2[1]:+.3f},{e2[2]:+.3f})", flush=True)
        print(f"  rms velocities <v_1^2>^1/2={s1:.4f} <v_2^2>^1/2={s2:.4f}  -> "
              f"xi_1/xi_2 = {s1/s2:.3f}"
              + ("  (axes rescaled)" if aniso else "  (NOT rescaled)"), flush=True)
        if not keep.all():
            print(f"  dropped {int((~keep).sum())} of {len(keep)} FS points whose velocity "
                  f"is nearly parallel to B", flush=True)
    return out


def reduce_fs_trajectories(fs, gap_sym=None, nbeta: int = 48, nv: int = 4, nphi: int = 8,
                           wmin: float = 1.0e-5, verbose: bool = True):
    """
    @fn reduce_fs_trajectories
    @brief Reduce a Fermi surface to a small set of WEIGHTED in-plane trajectories for
    the 2D-plane Riccati solvers (vortex, vortex lattice), whose cost is strictly linear
    in the number of trajectories: each one is a separate chord integration over the
    whole (x,y) grid, for every frequency, on every self-consistency iteration.  A 3D
    Fermi surface has 10^3-10^4 points, against the 24 directions of the model cylinder.

    Those solvers see a trajectory through exactly three numbers -- the in-plane
    direction beta, the in-plane speed |v_par| (the chord velocity) and the form factor
    phi -- so points are quantized on a (beta, |v_par|, Re phi, Im phi) grid and each
    occupied cell becomes ONE trajectory carrying the weighted mean of the three and the
    summed weight.

    phi is part of the key deliberately: two points can share a direction and a speed
    while sitting on sheets of OPPOSITE gap sign (s+-), and merging those would cancel a
    gap that is physically there.  |v_par| is binned by QUANTILE, so the bins follow the
    actual velocity distribution instead of being wasted on a sparse tail.  phi is
    rescaled back to the <|phi|^2> = 1 convention on the reduced set, so lambda keeps its
    meaning.

    Not for the surface solver: that one needs k_y, k_z and the band index to find the
    specular partner, which this reduction throws away.

    @param      fs: FS dict (build_model_fs / build_wannier_fs)
    @param gap_sym: pairing symmetry (ignored when fs already carries 'phi')
    @param   nbeta: direction bins over the full circle
    @param      nv: |v_par| quantile bins;  nphi: bins per phi component.
           Measured against the FULL 1920-direction d-wave vortex of a corrugated
           cylinder (ngrid=41, Lxi=6): the deviations shrink monotonically, so this is
           an ordinary discretization -- raise the bins until the answer stops moving.

             bins        dirs  speedup   dDbulk   rms d|Psi|   dN(core,0)  max_w dN
             (32,3,6)     108    x28.8    1.10%     0.013        0.063       0.230
             (48,4,8)     176    x17.7    0.54%     0.012        0.015       0.117   <- default
             (64,5,10)    284    x11.0    0.22%     0.004        0.017       0.082

           (N(core,0) itself is 7.05, so the default is ~0.2% on the core peak.)
    @param    wmin: drop cells below this fraction of the total weight
    @return a light FS dict (th,kx,ky,vx,vy,vabs,nf,phi,ntraj,nfs_full) for the solvers
    """
    beta = np.arctan2(np.asarray(fs['vy']), np.asarray(fs['vx']))
    v = np.asarray(fs['vabs'], dtype=np.float64)
    w = np.asarray(fs['nf'], dtype=np.float64)
    phi = np.asarray(fs_form_factor(fs, gap_sym), dtype=np.complex128)
    n0 = len(v)

    def _bin(x, nb):                                   # uniform bins over the data range
        if nb <= 1 or x.max() - x.min() < 1e-12 * max(1.0, abs(x).max()):
            return np.zeros(len(x), dtype=np.int64), 1
        e = np.linspace(x.min(), x.max(), nb + 1)[1:-1]
        return np.searchsorted(e, x).astype(np.int64), nb

    ib = np.minimum(((beta + np.pi) / (2.0 * np.pi) * nbeta).astype(np.int64), nbeta - 1)
    if nv > 1:                                         # quantile bins follow the |v| spread
        qe = np.unique(np.quantile(v, np.linspace(0.0, 1.0, nv + 1)))[1:-1]
        iv, nvb = np.searchsorted(qe, v).astype(np.int64), max(len(qe) + 1, 1)
    else:
        iv, nvb = np.zeros(n0, dtype=np.int64), 1
    ir, nrb = _bin(phi.real, nphi)
    ii, nib = _bin(phi.imag, nphi if np.abs(phi.imag).max() > 1e-12 else 1)

    key = ((ib * nvb + iv) * nrb + ir) * nib + ii
    _, inv = np.unique(key, return_inverse=True)
    ws = np.bincount(inv, weights=w)
    rep = lambda q: np.bincount(inv, weights=w * q) / ws
    br = np.arctan2(rep(np.sin(beta)), rep(np.cos(beta)))          # circular mean
    vr, pr = rep(v), rep(phi.real) + 1j * rep(phi.imag)

    keep = ws > wmin * ws.sum()
    br, vr, pr, ws = br[keep], vr[keep], pr[keep], ws[keep]
    ws = ws / ws.sum()
    nrm = np.sqrt((ws * np.abs(pr) ** 2).sum())                    # restore <|phi|^2> = 1
    if nrm > 0:
        pr = pr / nrm
    if verbose:
        print(f"  FS trajectories: {n0} points -> {len(ws)} weighted directions "
              f"(nbeta={nbeta}, nv={nv}, nphi={nphi}); <|v_par|> "
              f"{(w * v).sum():.4f} -> {(ws * vr).sum():.4f}", flush=True)
    out = dict(th=br, kx=np.cos(br), ky=np.sin(br), vx=vr * np.cos(br), vy=vr * np.sin(br),
               vabs=vr, vhx=np.cos(br), vhy=np.sin(br), nf=ws, phi=pr,
               ntraj=len(ws), nfs_full=n0)
    if 'nkz' in fs:
        out['nkz'] = fs['nkz']
    return out


_INT_GAP_STR = {0: 's', 1: 'd', 2: 's', 3: 'dxy', 6: 'd+id',
                -1: 'px', -2: 'py', -3: 'p+ip'}   # int -> continuum phi


def _gap_sym_str(gap_sym):
    """Map the global integer gap_sym index (gap_symms convention) to the continuum
    string used by the model-FS/cylinder routines; pass strings through unchanged.
    Indices without a 2D continuum harmonic (e.g. 4 dxz, 5 dyz) raise ValueError
    instead of silently falling back to 's'."""
    if isinstance(gap_sym, (int, np.integer)):
        try:
            return _INT_GAP_STR[int(gap_sym)]
        except KeyError:
            raise ValueError(
                f"gap_sym={int(gap_sym)} has no 2D continuum form factor "
                f"(supported indices: {sorted(_INT_GAP_STR)}); kz-dependent symmetries "
                "like 4 dxz / 5 dyz / 7 dxz+idyz need a 3D Fermi surface (the gap_symms "
                "lattice harmonics on a build_fs mesh or a build_wannier_fs(nkz>1) FS)") from None
    return gap_sym


def form_factor(beta: np.ndarray, gap_sym, beta_surf: float = 0.0) -> np.ndarray:
    """
    @fn form_factor
    @brief Continuum pairing form factor phi(beta) on the cylindrical FS,
    normalized so that the full-circle average <|phi|^2> = 1 (so the coupling
    lambda keeps its bulk BCS meaning).  May be complex (chiral / triplet
    states).  Shared by the surface and vortex solvers.
    @param     beta: FS angle(s) [rad]
    @param  gap_sym: singlet: 's', 'd' (d_{x^2-y^2}), 'dxy', chiral 'd+id'/'d-id';
                     triplet (odd parity): 'px', 'py', 'p+ip' / 'p-ip' (chiral);
                     or an integer (gap_symms index, mapped by _gap_sym_str)
    @param beta_surf: rotation of the gap relative to the reference axis
                      (d-wave: 0 -> [100] surface, pi/4 -> [110] surface)
    @note Triplet states are treated in the single (pseudo-)spin sector with a
          fixed d-vector: the equal-spin chiral state p+ip = e^{i beta} is fully
          gapped in the bulk and carries topological edge / core states.
    """
    gap_sym = _gap_sym_str(gap_sym)
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
    if gap_sym == 'd+id':          # d_{x2-y2} + i d_{xy} = e^{2 i beta} (|phi|=1 already)
        return np.exp(2j * a) * np.ones_like(beta, dtype=np.complex128)
    if gap_sym == 'd-id':
        return np.exp(-2j * a) * np.ones_like(beta, dtype=np.complex128)
    raise ValueError(f"unknown gap_sym: {gap_sym}")


def fs_form_factor(fs: dict, gap_sym) -> np.ndarray:
    """Pairing form factor phi at each FS point, normalized to nf-weighted <|phi|^2>=1.

    If the FS dict already carries a precomputed ``phi`` (e.g. a Wannier multiband FS
    built with a gap_sym/delta0; the per-band amplitudes and signs are baked in), it is
    returned directly.  Otherwise phi is built from ``gap_sym``: a STRING uses the
    continuum angular harmonics (s, d, dxy, px, py, chiral p+ip/p-ip); an INTEGER uses
    the same index convention as the rest of pyrpa (gap_symms: 0 s, 1 d_{x^2-y^2},
    2 s+-, 3 dxy, -1 px, -2 py) -- on a lattice FS (key 'kf') via the lattice harmonics,
    otherwise via the continuum angle."""
    if 'phi' in fs:
        return fs['phi']
    if isinstance(gap_sym, (int, np.integer)):
        if 'kf' in fs:                                   # lattice harmonic on the real FS
            phi = gap_symms(fs['kf'], 1, int(gap_sym))[0].astype(np.complex128)
            norm = np.sqrt((fs['nf'] * np.abs(phi) ** 2).sum())
            return phi / norm if norm > 0 else phi
        gap_sym = _gap_sym_str(gap_sym)                  # continuum fallback (model FS)
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
    elif gap_sym == 'd+id':
        phi = np.exp(2j * a)
    elif gap_sym == 'd-id':
        phi = np.exp(-2j * a)
    else:
        raise ValueError(f"unknown gap_sym: {gap_sym}")
    norm = np.sqrt((fs['nf'] * np.abs(phi) ** 2).sum())
    return phi / norm if norm > 0 else phi


def set_fs_gap(fs: dict, gap_sym, delta0=None):
    """Bake the pairing form factor into a (Wannier) FS dict: phi = gap_symms(k,gap_sym)
    with per-band amplitude and sign from ``delta0`` (indexed by fs['band']), normalized
    to nf-weighted <|phi|^2>=1.  Stores fs['phi'] and returns fs.
    @param gap_sym: integer symmetry index (gap_symms convention)
    @param  delta0: per-band gap amplitudes/signs (None = uniform 1); array indexed by band
    """
    phi = gap_symms(fs['kf'], 1, int(gap_sym))[0].astype(np.complex128)
    if delta0 is not None:
        d0 = np.asarray(delta0, dtype=np.float64)
        phi = phi * d0[fs['band']]                       # band-resolved ratio + sign
    norm = np.sqrt((fs['nf'] * np.abs(phi) ** 2).sum())
    if norm <= 0.0:
        # A k_z-dependent harmonic (4 dxz, 5 dyz, 7 dxz+idyz) vanishes identically on
        # the k_z = 0 plane -- exactly where a single-slice FS lives.  Silently
        # returning phi == 0 would make the gap zero for no visible reason.
        print(f"Warning: form factor gap_sym={gap_sym} vanishes on this Fermi surface "
              f"({fs.get('nkz', 1)} k_z slice(s)).  A k_z-dependent gap needs a 3D FS: "
              "build_wannier_fs(..., nkz>1) / eil_fs_nkz>1", flush=True)
    fs['phi'] = phi / norm if norm > 0 else phi
    return fs


def bulk_gap_fs(coupling, temp, omega, fs, gap_sym, eps=1e-8, itemax=500, mix=0.5):
    """Homogeneous gap amplitude on a model FS (nf-weighted FS average)."""
    phi = fs_form_factor(fs, gap_sym)
    nf = fs['nf']
    a2 = (phi * np.conj(phi)).real
    damp = BCS_RATIO * temp
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
    @brief Superfluid-density tensor (rho_xx, rho_yy, rho_zz)/rho_n(0) on a model FS
    with Fermi velocities.  rho_ii(T)/rho_ii(0) = pi <v_i^2 2T sum_n Dk^2/(w^2+Dk^2)^{3/2}>_nf
    / <v_i^2>_nf (the bracket -> 1 as T->0).  The weight uses the FULL velocity
    components v_i (rho_ii ~ integral dS/|v_F| v_i^2), not the unit direction:
    with nf = dl/|v_par| this weights fast FS sections by |v_F|^2, which sets both
    the anisotropy ratio and the T dependence.  Anisotropic FS / nodes give an
    anisotropic penetration depth lambda_ii = 1/sqrt(rho_ii).

    rho_zz is the out-of-plane (c-axis) response, so lambda_c/lambda_ab = sqrt(rho_xx/rho_zz)
    measures the mass anisotropy directly.  It needs a 3D Fermi surface: on a single
    k_z slice (or an in-plane model kind) v_z == 0 everywhere and rho_zz is returned as
    0 rather than 0/0.
    @return (Delta, rho_xx, rho_yy, rho_zz)
    """
    omega = matsubara(temp, wc)
    Delta = bulk_gap_fs(coupling, temp, omega, fs, gap_sym)
    if Delta <= 1.0e-6 * temp:
        return 0.0, 0.0, 0.0, 0.0
    phi = fs_form_factor(fs, gap_sym)
    Dk = Delta * np.abs(phi)
    fac = 2.0 * temp * (Dk[:, None] ** 2 /
                        (omega[None, :] ** 2 + Dk[:, None] ** 2) ** 1.5).sum(axis=1)
    nf = fs['nf']

    def rho(vi2):                                      # pi <v_i^2 fac> / <v_i^2>
        den = (nf * vi2).sum()
        return np.pi * (nf * vi2 * fac).sum() / den if den > 0.0 else 0.0

    vz = fs.get('vz')
    return (Delta, rho(fs['vx'] ** 2), rho(fs['vy'] ** 2),
            rho(vz ** 2) if vz is not None else 0.0)


def calc_fs_penetration(coupling, temp, wc, kind='ellipse', gap_sym='s', Nth=360,
                        mu=None, params=None, t_list=None, kb=1.0, nkz=1,
                        kz_max=np.pi, fs=None):
    """
    @fn calc_fs_penetration
    @brief Temperature sweep of the anisotropic superfluid density / penetration
    depth on a Fermi surface with Fermi velocities.  An anisotropic FS gives
    lambda_xx != lambda_yy, and nodes give the linear-in-T penetration depth.
    With a 3D FS (nkz>1, or a prebuilt ``fs``) the out-of-plane rho_zz and the mass
    anisotropy lambda_c/lambda_ab are reported as well.
    Writes 'fs_penetration.dat'.
    @param nkz,kz_max: k_z stacking of the model FS (see build_model_fs)
    @param     fs: prebuilt FS dict (build_model_fs / build_wannier_fs), overriding kind
    """
    fs = build_model_fs(kind, Nth, mu, params, nkz=nkz, kz_max=kz_max) if fs is None else fs
    if not (isinstance(gap_sym, (int, np.integer)) and 'kf' in fs):
        gap_sym = _gap_sym_str(gap_sym)               # continuum harmonic (no lattice FS)
    d3 = fs.get('nkz', 1) > 1
    # full-velocity weights <v_i^2>_nf (nf = dl/|v_par|), consistent with superfluid_density_fs
    vx2 = (fs['nf'] * fs['vx'] ** 2).sum()
    vy2 = (fs['nf'] * fs['vy'] ** 2).sum()
    vz2 = (fs['nf'] * fs.get('vz', np.zeros_like(fs['vx'])) ** 2).sum()
    print(f"model FS '{kind}' (Nth={Nth}, {fs.get('nkz', 1)} k_z slice(s), "
          f"{len(fs['kx'])} pts): <|v_par|>={fs['vabs'].mean():.4f}, "
          f"<v_x^2>/<v_y^2>/<v_z^2> = {vx2:.3f}/{vy2:.3f}/{vz2:.3f}, "
          f"absolute lambda_xx/lambda_yy(0) = {np.sqrt(vy2/vx2):.3f}"
          + (f", lambda_c/lambda_ab(0) = {np.sqrt(vx2/vz2):.3f}" if vz2 > 0 else ""),
          flush=True)
    print(f"penetration depth on model FS: gap {gap_sym}, lambda={coupling:.3f}", flush=True)
    if t_list is None:
        t_hi = temp
        for _ in range(8):
            if superfluid_density_fs(coupling, t_hi, wc, fs, gap_sym)[0] > 0:
                break
            t_hi *= 0.6
        t_list = np.linspace(0.05 * t_hi, 1.02 * t_hi, 14)
    print("  T[K]    Delta     rho_xx   rho_yy   rho_zz   lam_xx/lam_yy  lam_c/lam_ab",
          flush=True)
    rows = []
    for T in t_list:
        D, rxx, ryy, rzz = superfluid_density_fs(coupling, T, wc, fs, gap_sym)
        ratio = np.sqrt(ryy / rxx) if (rxx > 1e-9 and ryy > 1e-9) else np.nan
        rat_c = np.sqrt(rxx / rzz) if (rxx > 1e-9 and rzz > 1e-9) else np.nan
        rows.append((T, D, rxx, ryy, rzz, ratio, rat_c))
        print(f"  {T/kb:6.2f}  {D:9.3e}  {rxx:7.4f}  {ryy:7.4f}  {rzz:7.4f}   "
              f"{ratio:9.4f}     {rat_c:9.4f}", flush=True)
    try:
        with open('fs_penetration.dat', 'w') as fh:
            fh.write("# T[eV]  Delta[eV]  rho_xx  rho_yy  rho_zz  lambda_xx/lambda_yy"
                     "  lambda_c/lambda_ab\n")
            for T, D, rxx, ryy, rzz, ratio, rat_c in rows:
                fh.write(f"{T:12.6e} {D:12.6e} {rxx:12.6e} {ryy:12.6e} {rzz:12.6e} "
                         f"{ratio:12.6e} {rat_c:12.6e}\n")
    except IOError as e:
        print(f"Error writing fs_penetration.dat: {e}", flush=True)
    return rows


# =========================================================================== #
#  Condensation free energy
# --------------------------------------------------------------------------- #
# Homogeneous superconducting condensation free energy (Omega_s - Omega_n)/N0 in the
# standard quasiclassical / BCS Matsubara form,
#   dOmega/N0 = -2 pi T sum_{w_n>0} < E_k - w_n - |Delta_k|^2/(2 E_k) >_FS,
#   E_k = sqrt(w_n^2 + |Delta_k|^2),  |Delta_k| = Damp |phi(k)|.
# The summand decays as |Delta|^4/w_n^3, so the sum is cutoff-convergent and the
# result is coupling-independent: a universal ratio dOmega(0)/Damp0^2 and dOmega->0
# at Tc.  (The *inhomogeneous* case needs the Eilenberger free-energy functional /
# coupling-constant integration with the Riccati f -- gated on the self-consistent
# spatial order parameter.)


def free_energy_weights(gap_sym='s', fs=None, phi=None, wnf=None, Nb=240):
    """
    @fn free_energy_weights
    @brief Resolve the (phi, weights, label) triple the free-energy routines average
    over, from any of the three ways of specifying a gap.  Factored out so the
    condensation energy is not restricted to the handful of in-plane harmonics it
    used to hard-code -- the chiral-versus-real selection of a degenerate pair is
    decided by this quantity, and both members need to be expressible.

      * ``phi`` given   -- an explicit form factor with weights ``wnf`` (default
        uniform); this is what build_fs returns, so a gap on the full 3D k-mesh
        (band projected, k_z dependent, multi-sheet) can be fed straight in.
      * ``fs`` given    -- a model / Wannier FS dict; an INTEGER gap_sym then uses the
        gap_symms lattice harmonics on fs['kf'], so the k_z dependent 4/5/7 work.
      * neither         -- the isotropic cylinder, via the continuum ``form_factor``
        (every continuum harmonic, chiral 'd+id' / 'p+ip' included).
    @return (phi [Nfs] complex, wnf [Nfs] summing to 1, label)
    """
    if phi is not None:
        phi = np.asarray(phi, dtype=np.complex128)
        w = np.full(len(phi), 1.0 / len(phi)) if wnf is None else np.asarray(wnf, float)
        return phi, w / w.sum(), 'given phi'
    if fs is not None:
        return (fs_form_factor(fs, gap_sym).astype(np.complex128),
                np.asarray(fs['nf'], float), str(gap_sym))
    beta = np.linspace(0.0, 2.0 * np.pi, Nb, endpoint=False)
    lab = _gap_sym_str(gap_sym)
    return (form_factor(beta, lab).astype(np.complex128),
            np.full(Nb, 1.0 / Nb), lab)


def condensation_energy(coupling, temp, wc, gap_sym='s', fs=None, Damp=None,
                        phi=None, wnf=None):
    """
    @fn condensation_energy
    @brief Homogeneous superconducting condensation free energy (Omega_s-Omega_n)/N0
    [eV^2] at temperature ``temp`` (coupling-independent quasiclassical Matsubara form).
    @param gap_sym: pairing symmetry; with a model/Wannier FS pass ``fs`` (else cylinder)
    @param   Damp: bulk gap [eV] (computed self-consistently if None)
    @param phi,wnf: explicit form factor and FS weights (e.g. straight from build_fs),
                    overriding gap_sym/fs -- see free_energy_weights
    @return (dOmega, Damp): condensation energy per N0 [eV^2] (<=0) and the bulk gap
    """
    omega = matsubara(temp, wc)
    phi, wnf, _ = free_energy_weights(gap_sym, fs, phi, wnf)
    if Damp is None:
        Damp = solve_gap(temp, wnf, phi, omega, 0.0, 1.0e8, coupling)
    if Damp <= 1.0e-6 * temp:
        return 0.0, 0.0
    Dk2 = (Damp ** 2) * (phi * np.conj(phi)).real        # |Delta_k|^2 per FS point
    E = np.sqrt(omega[None, :] ** 2 + Dk2[:, None])      # [Nfs, Nw]
    summand = E - omega[None, :] - Dk2[:, None] / (2.0 * E)
    dOmega = -2.0 * np.pi * temp * (wnf[:, None] * summand).sum()   # FS-avg + 2T sum_{n>0}
    return dOmega, Damp


def calc_free_energy(coupling, temp, wc, gap_sym='s', fs=None, t_list=None, kb=1.0,
                     phi=None, wnf=None, fname='free_energy.dat'):
    """
    @fn calc_free_energy
    @brief Temperature sweep of the homogeneous condensation free energy
    (Omega_s-Omega_n)/N0 (coupling-constant integration).  Writes 'free_energy.dat'.
    @param fs,phi,wnf: how the gap is specified -- see free_energy_weights.  Any
           harmonic works, k_z dependent (4 dxz, 5 dyz, 7 dxz+idyz on a 3D FS) and
           chiral included, so two members of a degenerate pair can be compared:
           the one with the LOWER dOmega is the realized state, which is the only
           way to choose between them (linear theory leaves them degenerate at Tc).
    @param fname: output file (rename it when comparing several states)
    """
    phi, wnf, lab = free_energy_weights(gap_sym, fs, phi, wnf)
    if t_list is None:
        t_hi = temp
        for _ in range(8):
            if solve_gap(t_hi, wnf, phi, matsubara(t_hi, wc), 0.0, 1e8, coupling) > 0:
                break
            t_hi *= 0.6
        t_list = np.linspace(0.05 * t_hi, 1.05 * t_hi, 16)
    src = 'given phi' if lab == 'given phi' else ('FS dict' if fs is not None else 'cylinder')
    print(f"condensation free energy: gap {lab}, lambda={coupling:.3f}, {src} "
          f"({len(phi)} FS points)", flush=True)
    print("  T[K]      Delta[eV]     dOmega/N0 [eV^2]   dOmega/Delta0^2", flush=True)
    rows = []
    D0 = None
    for T in t_list:
        dO, D = condensation_energy(coupling, T, wc, phi=phi, wnf=wnf)
        if D0 is None and D > 0:
            D0 = D
        rows.append((T, D, dO))
    D0 = D0 or 1.0
    try:
        with open(fname, 'w') as fh:
            fh.write("# T[eV]  Delta[eV]  dOmega/N0[eV^2]  dOmega/Delta0^2\n")
            for T, D, dO in rows:
                print(f"  {T/kb:7.3f}  {D:12.5e}  {dO:14.6e}   {dO/D0**2:9.4f}", flush=True)
                fh.write(f"{T:12.6e} {D:12.6e} {dO:14.6e} {dO/D0**2:12.6e}\n")
    except IOError as e:
        print(f"Error writing {fname}: {e}", flush=True)
    return rows


def project_gap_to_band(fs, gap_orbital, normalize=True, gauge=None, spin_map=None):
    """
    @fn project_gap_to_band
    @brief Band-basis pairing form factor by the systematic low-energy projection of an
    orbital-basis pair potential onto the Fermi-surface bands (Nagai-Nakamura, JPSJ 85,
    074707 (2016), Eq. 43):  Delta_eff(k_F) = U(k_F)^dagger . Delta_orb(k_F) . U*(-k_F),
    where U(k_F)=fs['evec'] is the (orbital-content) eigenvector of the FS band and
    U*(-k_F)=conj(fs['evec_mk']).  For one band per sheet (M=1) this is the scalar
    phi(k_F)=sum_{ab} conj(u_a(k_F)) Delta_orb,ab(k_F) conj(u_b(-k_F)); the orbital
    character of the band (how u varies over/between sheets) is thus built into the gap
    -- e.g. orbital-selective or s+- structure emerges, rather than being imposed by delta0.
    @param fs: a Wannier FS dict carrying 'evec' and 'evec_mk' (build_wannier_fs)
    @param gap_orbital: the N x N orbital pair-potential matrix Delta_orb -- either a
           constant array (k-independent on-site/orbital pairing) or a callable
           kfrac(3,) -> N x N array (k-dependent, e.g. an extended/d-wave bond pairing)
    @param normalize: rescale to the nf-weighted <|phi|^2> = 1 convention
    @param gauge,spin_map: pair-partner gauge overriding the module default
           (see set_pair_gauge / pair_partner_conj); 'trs' (default) is gauge
           invariant, the legacy 'diag' route is not
    @return fs (with fs['phi'] set)
    """
    ev, evm = fs['evec'], fs['evec_mk']                # [Nfs, Norb]
    check_pair_partner(ev, evm, gauge, spin_map)       # warn if H breaks time reversal
    phi = project_gap(ev, gap_orbital, fs['kf'], evm, gauge, spin_map)
    if normalize:
        nrm = np.sqrt((fs['nf'] * np.abs(phi) ** 2).sum())
        phi = phi / nrm if nrm > 0 else phi
    fs['phi'] = phi
    return fs


def gap_orbital_from_wannier(fname='gap_wannier', iw_index=0, n_avg=1):
    """
    @fn gap_orbital_from_wannier
    @brief Load an RPA/FLEX gap Delta(R, iw_n) written by output_gap_wannier (the gap
    exported in Wannier-real-space "hopping" format) and return a callable
    kfrac(3,) -> Delta_orb(k) (the inverse Fourier transform on the orbital matrix),
    ready to pass as build_wannier_fs's ``gap_orbital`` -- which then projects it onto the
    Fermi-surface bands (Nagai-Nakamura Eq 43, project_gap_to_band) to set fs['phi'].
    This is the bridge that uses a self-consistent RPA/FLEX gap (e.g. the KFe2As2 gap of
    PRB 84, 144514) as the pairing form factor for the quasiclassical vortex(-lattice)
    solvers.

    The Eilenberger form factor is a STATIC phi(k_F); the lowest Matsubara slice iw_0
    (iw_index=0) carries the full symmetry / sign / node / anisotropy structure most
    sharply and matches the gap usually quoted on the Fermi surface.  n_avg>1 averages the
    first n_avg consecutive slices for noise robustness (smoother phi at the cost of
    slightly diluting the anisotropy, since Delta(k,iw_n) gets more isotropic with n).
    The absolute scale is irrelevant: project_gap_to_band renormalizes to <|phi|^2>=1.

    Fourier convention: output_gap_wannier builds Delta(R) with np.fft.ifftn over the
    k-axes, whose kernel is exp(+2 pi i k.R), so the inverse taken here carries the MINUS
    sign, matching the self-energy interpolation in main.py's plot_spectrum.  (This used
    to be written with a plus and therefore returned Delta(-k).  It changed nothing in
    practice -- the result is only ever used as a form factor that project_gap_to_band
    renormalizes, and Delta(-k)=Delta(k) for the even-parity gaps this path is used with,
    while an odd-parity one only picks up a global sign -- but the transform is now the
    true inverse, verified to round-trip an output_gap_wannier export to 3e-15.)

    @param fname: base name of the .npz written by output_gap_wannier (without extension)
    @param iw_index: starting Matsubara index (0 = lowest iw_0)
    @param n_avg: number of consecutive Matsubara slices to average (1 = single slice)
    @return callable kfrac(3,) -> Delta_orb(k) [Norb, Norb] complex
    """
    d = np.load(f'{fname}.npz')
    gap, rvec = d['gap'], d['rvec']                  # (Norb,Norb,N_cut,Nr), (Nr,3)
    n0, n1 = iw_index, min(iw_index + n_avg, gap.shape[2])
    if n0 >= gap.shape[2]:
        raise ValueError(f"iw_index={iw_index} out of range (N_iw={gap.shape[2]})")
    gR = gap[:, :, n0:n1, :].mean(axis=2)            # (Norb,Norb,Nr) over chosen slices
    rg = np.ascontiguousarray(rvec, dtype=np.float64)

    No, Nr = gR.shape[0], gR.shape[2]
    gR2 = gR.reshape(No * No, Nr)                    # BLAS-friendly (Norb^2, Nr)

    def gap_orbital(kfrac):                          # Delta_orb(k) = sum_R e^{-i 2pi k.R} Delta(R)
        k = np.asarray(kfrac, dtype=np.float64)
        if k.ndim == 1:                              # single k -> (Norb, Norb)
            return (gR * np.exp(-2j * np.pi * (rg @ k))).sum(axis=-1)
        out = np.empty((len(k), No, No), dtype=np.complex128)   # batch -> (Nk, Norb, Norb)
        for i0 in range(0, len(k), 4096):            # chunked: keeps the phase array small
            kc = k[i0:i0 + 4096]
            out[i0:i0 + len(kc)] = (gR2 @ np.exp(-2j * np.pi * (rg @ kc.T))).T.reshape(len(kc), No, No)
        return out
    return gap_orbital


def build_wannier_fs(rvec, ham_r, S_r, avec, mu, mesh=360, kz=0.0, nkz=1, kz_list=None,
                     band=None, RotMat=None, gap_sym=None, delta0=None, gap_orbital=None,
                     gauge=None, spin_map=None):
    """
    @fn build_wannier_fs
    @brief Build the Fermi surface (with Fermi velocities) of a REAL Wannier tight-
    binding model for the quasiclassical trajectory solvers -- the band-structure
    counterpart of build_model_fs.  The FS contour at mu is extracted from the 2D band
    mesh (get_kf_points), the group velocity v_F = grad_k eps is the Wannier band
    velocity (get_veloc), and the DOS weight is nf = dl/|v_F| (arc length over speed).
    @param rvec,ham_r,S_r: Wannier R-vectors, hopping, overlap (S_r=[] for orthogonal)
    @param avec: lattice vectors (velocity metric; Cartesian k = 2 pi * fractional)
    @param   mu: chemical potential [eV];  mesh: 2D k-mesh
    @param   kz: fixed k_z slice used when nkz == 1 and kz_list is None
    @param  nkz: number of k_z slices stacked into a full 3D Fermi surface
                 (1 = the historical single-slice, quasi-2D FS).  The slices are
                 uniform on k_z in [-1/2, 1/2), so every sheet and the whole k_z
                 dispersion enter -- required for k_z-dependent gaps (horizontal
                 line nodes, dxz/dyz/dxz+idyz), which vanish identically on the
                 k_z = 0 plane that the single-slice FS is pinned to.
    @param kz_list: explicit list of k_z slices (overrides nkz/kz)
    @param band: keep only this band index (None = all FS sheets)
    @param gap_sym: if given (integer index, gap_symms convention), bake in the pairing
                    form factor phi(k) via set_fs_gap (lattice harmonics on the real FS)
    @param  delta0: per-band gap amplitudes/signs (indexed by band) for a multiband gap
                    (e.g. s+- with opposite signs on different sheets); None = uniform
    @param gap_orbital: N x N orbital-basis pair potential (array or callable kfrac->NxN);
                    if given, the band gap phi(k) is the low-energy PROJECTION of it onto
                    the FS bands (Nagai-Nakamura Eq 43, project_gap_to_band) -- the orbital
                    character is built in, superseding the phenomenological gap_sym/delta0
    @param gauge,spin_map: pair-partner gauge for the projection (set_pair_gauge)
    @return dict th,kx,ky,kz,kf,band,vx,vy,vz,vabs,vabs3,vhx,vhy,nf,nkz,evec,evec_mk [,phi]

    @note  ``vabs`` (and vhx,vhy) is the IN-PLANE speed |v_par| = sqrt(vx^2+vy^2),
           and it is the right quantity in both places it is used, also in 3D:
           for a z-independent Delta(x,y) the transport term is v_F.grad =
           v_par.grad_par, so the chord velocity is |v_par|; and the Fermi-surface
           measure factorizes as
               int dS/|v_F| = int dk_z ( contour_int dl / |v_par| ),
           so the per-slice weight nf = dl/|v_par| with a uniform dk_z (which
           cancels in the normalization) is the exact 3D DOS weight.  The full
           3D velocity is kept as vz / vabs3 for quantities that need it
           (e.g. rho_zz, out-of-plane response).
    """
    from ._bands import get_eigs_2d, get_kf_points, get_eigs
    from .. import flibs
    if RotMat is None:
        RotMat = np.eye(3)
    if kz_list is not None:
        kzs = np.asarray(kz_list, dtype=np.float64)
    elif nkz > 1:
        kzs = np.linspace(-0.5, 0.5, int(nkz), endpoint=False)   # uniform dk_z
    else:
        kzs = np.array([float(kz)])
    kxs, kys, vxs, vys, vzs, nfs, kfs, bnd, evs, evms = ([] for _ in range(10))
    for kz_i in kzs:
        eig2d = get_eigs_2d(mesh, rvec, ham_r, S_r, RotMat, kz_i)
        kf, fsband = get_kf_points(eig2d, mesh, mu, kz_i)
        for contlist, b in zip(kf, fsband):
            if band is not None and b != band:
                continue
            for cont in contlist:                      # one connected FS piece
                pts = np.ascontiguousarray(cont, dtype=np.float64)   # [Np,3] fractional
                if len(pts) < 2:
                    continue
                eg, uni = get_eigs(pts, ham_r, S_r, rvec)
                # S_r/eg carry the -eps*dS/dk overlap term for a non-orthogonal (MLO) basis
                vk = flibs.get_veloc(pts, ham_r, rvec, avec.T, uni, S_r, eg).real  # [Np,Norb,3]
                v = vk[:, b, :]
                # band eigenvectors u(k_F) and u(-k_F) (orbital content) for the gap projection
                ub = uni[:, b, :].copy()               # [Np, Norb] (Eilenberger band a=b)
                ub /= np.sqrt((np.abs(ub) ** 2).sum(axis=1))[:, None]
                _, uni_m = get_eigs(np.ascontiguousarray(-pts), ham_r, S_r, rvec)
                umb = uni_m[:, b, :].copy()
                umb /= np.sqrt((np.abs(umb) ** 2).sum(axis=1))[:, None]
                kc = 2.0 * np.pi * pts[:, :2]          # Cartesian in-plane k (a=1)
                seg = np.sqrt((np.diff(kc, axis=0) ** 2).sum(axis=1))  # piece segments
                dl = np.zeros(len(pts))
                dl[:-1] += 0.5 * seg
                dl[1:] += 0.5 * seg
                vab = np.sqrt(v[:, 0] ** 2 + v[:, 1] ** 2)             # in-plane speed
                kxs.append(kc[:, 0]); kys.append(kc[:, 1])
                vxs.append(v[:, 0]); vys.append(v[:, 1]); vzs.append(v[:, 2])
                nfs.append(dl / np.maximum(vab, 1e-12))                # dk_z is uniform
                kfs.append(pts); bnd.append(np.full(len(pts), b, dtype=int))
                evs.append(ub); evms.append(umb)
    if not kxs:
        raise ValueError(f"no Fermi surface at this mu over {len(kzs)} k_z slice(s)")
    kx = np.concatenate(kxs); ky = np.concatenate(kys)
    vx = np.concatenate(vxs); vy = np.concatenate(vys); vz = np.concatenate(vzs)
    nf = np.concatenate(nfs)
    vabs = np.sqrt(vx ** 2 + vy ** 2)                  # in-plane |v_par| (see @note)
    nf = nf / nf.sum()
    kfrac = np.concatenate(kfs)
    fs = dict(th=np.arctan2(ky, kx), kx=kx, ky=ky, kz=kfrac[:, 2], kf=kfrac,
              band=np.concatenate(bnd), vx=vx, vy=vy, vz=vz, vabs=vabs,
              vabs3=np.sqrt(vx ** 2 + vy ** 2 + vz ** 2), nkz=len(kzs),
              vhx=vx / vabs, vhy=vy / vabs, nf=nf,
              evec=np.concatenate(evs), evec_mk=np.concatenate(evms))  # u(k_F), u(-k_F)
    if gap_orbital is not None:
        # band gap = gauge-fixed U^dag Delta_orb U  (Nagai-Nakamura Eq 43)
        project_gap_to_band(fs, gap_orbital, gauge=gauge, spin_map=spin_map)
    elif gap_sym is not None:
        set_fs_gap(fs, gap_sym, delta0)                # phenomenological per-band ratio/sign
    return fs


def gap_color_3d(centers, blist, rvec, ham_r, S_r, gap_sym=None, delta0=None, gap_orbital=None,
                 gauge=None, spin_map=None, sw_abs=None):
    """
    @fn gap_color_3d
    @brief Re[phi(k)] at each face center of the 3D Wannier Fermi surface
    (gen_3d_surf_points), using the EXACT SAME pairing-gap definition as the Eilenberger
    solvers (build_wannier_fs/set_fs_gap/project_gap_to_band) -- so FERMI_3D with
    color_option=ColorMode.GAP shows the actual gap driving the Eilenberger calculation,
    directly on the real 3D Fermi surface (all sheets/kz, not just the fixed-kz cut used
    by the vortex/surface solvers).  If gap_orbital is given, phi is the Nagai-Nakamura
    band projection U(k)^dagger . Delta_orb(k) . U(-k)^* (Eq.43); otherwise it is the
    phenomenological gap_sym/delta0 lattice harmonic (set_fs_gap).
    @param centers: list (per FS sheet) of face-center k-points [Nfaces,3] (fractional),
           as returned by gen_3d_surf_points
    @param   blist: band index per sheet (gen_3d_surf_points)
    @param rvec,ham_r,S_r: Wannier R-vectors, hopping, overlap (S_r=[] for orthogonal)
    @param gap_sym: phenomenological gap symmetry index (gap_symms convention)
    @param  delta0: per-band gap amplitude/sign (indexed by band); None = uniform 1
    @param gap_orbital: Norb x Norb orbital pair-potential (array or callable kfrac->NxN);
           supersedes gap_sym/delta0, same convention as build_wannier_fs
    @return clist: list (per FS sheet) of real arrays [Nfaces], Re[phi(k)]
    """
    from ._bands import get_eigs
    clist = []
    for k, b in zip(centers, blist):
        k = np.ascontiguousarray(k, dtype=np.float64)
        if gap_orbital is not None:
            _, uni = get_eigs(k, ham_r, S_r, rvec)
            ub = uni[:, b, :].copy()                              # u(k)  [Nfaces,Norb]
            ub /= np.sqrt((np.abs(ub) ** 2).sum(axis=1))[:, None]
            _, uni_m = get_eigs(-k, ham_r, S_r, rvec)
            umb = uni_m[:, b, :].copy()                           # u(-k) [Nfaces,Norb]
            umb /= np.sqrt((np.abs(umb) ** 2).sum(axis=1))[:, None]
            check_pair_partner(ub, umb, gauge, spin_map)
            phi = project_gap(ub, gap_orbital, k, umb, gauge, spin_map)
        else:
            phi = gap_symms(k, 1, int(gap_sym))[0].astype(np.complex128)
            if delta0 is not None:
                phi = phi * np.asarray(delta0, dtype=np.float64)[b]
        # A chiral phi has no meaningful sign: Re[phi] would draw spurious nodes where
        # only the phase winds, so colour such sheets by |phi| instead.
        use_abs = (np.abs(phi.imag).max() > 1e-12) if sw_abs is None else bool(sw_abs)
        clist.append(np.abs(phi) if use_abs else phi.real)
    return clist
