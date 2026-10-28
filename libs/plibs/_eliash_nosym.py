#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
Linearized Eliashberg equation without any symmetry (no TRS / inversion): Zeeman field, FFLO.

Everything lives on the full Brillouin zone in the spin-orbital basis of the SOC code.  A spinless
(non-SOC) model is spin-doubled here, so a Zeeman field can be applied to either kind of model.
The Fortran kernels are in libs/src/f_eliash_nosym.f90.

Flow of calc_lin_eliash_nosym:
  H(k) (+ Zeeman) on the full mesh -> G0 -> chi0 (once) -> V_Delta (once)
  -> for every pair momentum Q in qlist: lin_eliash_nosym -> lambda(Q)
"""
import numpy as np
import libs.flibs as flibs
import libs.plibs._bands as _bands
import libs.plibs._response as _response

# Public API of this submodule; keeps `from ._x import *` in __init__.py from
# re-exporting numpy and the sibling packages imported above.
__all__ = [
    'spin_double_model', 'spin_double_chiolist', 'add_zeeman', 'nosym_initial_gap', 'gap_character_nosym',
    'calc_lin_eliash_nosym'
]

_PAULI = (np.array([[0, 1], [1, 0]], dtype=np.complex128),
          np.array([[0, -1j], [1j, 0]], dtype=np.complex128),
          np.array([[1, 0], [0, -1]], dtype=np.complex128))


def spin_double_model(ham_r: np.ndarray, S_r: np.ndarray) -> tuple:
    """
    @brief Spin-double a spinless model: H(R) -> 1_2 (x) H(R), ordering (all up, all down) as main.py's SOC default.
    @param  ham_r: Real-space Hamiltonian blocks [Nr, Norb, Norb] complex128
    @param    S_r: Overlap blocks [Nr, Norb, Norb], or [] for an orthogonal basis
    @return ham_r2: Spin-doubled Hamiltonian [Nr, 2*Norb, 2*Norb] complex128
    @return   S_r2: Spin-doubled overlap (or [] unchanged)
    @return  slist: Spin label of each spin-orbital [2*Norb] int64 (+1 up, -1 down)
    @return   invs: 1-based index of the opposite-spin partner [2*Norb] int64
    """
    Norb = ham_r.shape[1]
    eye2 = np.eye(2)
    ham_r2 = np.ascontiguousarray(np.array([np.kron(eye2, h) for h in ham_r], dtype=np.complex128))
    S_r2 = S_r if len(S_r) == 0 else np.ascontiguousarray(np.array([np.kron(eye2, s) for s in S_r]))
    slist = np.concatenate([np.ones(Norb), -np.ones(Norb)]).astype(np.int64)
    invs = np.concatenate([np.arange(Norb, 2*Norb), np.arange(Norb)]).astype(np.int64) + 1
    return ham_r2, S_r2, slist, invs


def spin_double_chiolist(chiolist: np.ndarray, site: np.ndarray, Norb: int) -> tuple[np.ndarray, np.ndarray]:
    """
    @brief Spin-double a spinless orbital-pair list: each pair (a,b) on a site gives (a s, b s') for all s,s'.
    @param chiolist: Spinless orbital-pair list [Nchi, 2] int64 (1-based)
    @param     site: Site index of each pair [Nchi] int64
    @param     Norb: Number of spinless orbitals
    @return chiolist2: Spin-orbital pair list [4*Nchi, 2] int64 (1-based), FORTRAN-contiguous
    @return     site2: Site index of each spin-orbital pair [4*Nchi] int64
    """
    pairs, sites = [], []
    for s1 in (0, Norb):
        for s2 in (0, Norb):
            pairs.append(np.asarray(chiolist) + np.array([s1, s2]))
            sites.append(np.asarray(site))
    return (np.asfortranarray(np.concatenate(pairs).astype(np.int64)),
            np.concatenate(sites).astype(np.int64))


def add_zeeman(ham_r: np.ndarray, S_r: np.ndarray, rvec: np.ndarray, slist: np.ndarray, invs: np.ndarray,
               hvec: np.ndarray) -> np.ndarray:
    """
    @brief Add the Zeeman term H_Z = h.sigma (h in eV, i.e. g mu_B B/2 absorbed, spin only) to H(R).
    For an orthogonal basis it sits on R=0; with an overlap S(R) it is symmetrized as (Z S(R) + S(R) Z)/2,
    exact when S is spin diagonal with identical up/down blocks.
    @param  ham_r: Spin-orbital Hamiltonian [Nr, Norb, Norb] complex128 (not modified)
    @param    S_r: Overlap blocks [Nr, Norb, Norb], or [] for an orthogonal basis
    @param   rvec: R-vectors [Nr, 3]
    @param  slist: Spin label of each orbital [Norb] int64
    @param   invs: 1-based opposite-spin partner [Norb] int64
    @param   hvec: Zeeman field (hx, hy, hz) [eV]
    @return ham_r: New Hamiltonian with the Zeeman term [Nr, Norb, Norb] complex128
    """
    Norb = ham_r.shape[1]
    hs = sum(h*p for h, p in zip(hvec, _PAULI))   # 2x2 h.sigma, rows/cols (up, down)
    Z = np.zeros((Norb, Norb), dtype=np.complex128)
    for i in range(Norb):
        if slist[i] != 1:
            continue
        j = invs[i] - 1
        Z[i, i], Z[i, j], Z[j, i], Z[j, j] = hs[0, 0], hs[0, 1], hs[1, 0], hs[1, 1]
    out = np.array(ham_r, dtype=np.complex128, copy=True)
    if len(S_r) == 0:
        i0 = np.where(np.all(np.isclose(rvec, 0.0), axis=1))[0]
        if len(i0) != 1:
            raise ValueError('add_zeeman: rvec has no unique R=0 entry')
        out[i0[0]] += Z
    else:
        for ir, s in enumerate(S_r):
            out[ir] += 0.5*(Z @ s + s @ Z)
    return np.ascontiguousarray(out)


def nosym_initial_gap(klist: np.ndarray, slist: np.ndarray, invs: np.ndarray, gap_sym: int) -> np.ndarray:
    """
    @brief Orbital-diagonal spin-orbital seed for lin_eliash_nosym: Delta_{a up, a down}(k) = f(k),
    Delta_{a down, a up}(k) = -f(k) (gap_sym>=0, singlet) or +f(k) (gap_sym<0, Sz=0 triplet),
    with f the gap_sym form factor of _response.get_initial_gap.  The solver does not project onto
    this channel; the seed only sets the starting vector.
    @param   klist: Full-mesh k-points in fractional coordinates [Nkall, 3]
    @param   slist: Spin label of each orbital [Norb] int64
    @param    invs: 1-based opposite-spin partner [Norb] int64
    @param gap_sym: Symmetry index (see _response.get_initial_gap)
    @return   init: Initial gap [Norb, Norb, Nkall] complex128
    """
    Norb = len(slist)
    f = np.asarray(_response.get_initial_gap(klist, 1, gap_sym)[0], dtype=np.complex128)
    init = np.zeros((Norb, Norb, len(klist)), dtype=np.complex128)
    sgn = -1.0 if gap_sym >= 0 else 1.0
    for i in range(Norb):
        if slist[i] != 1:
            continue
        j = invs[i] - 1
        init[i, j] = f
        init[j, i] = sgn*f
    return init


def _stoner_nosym(chi0: np.ndarray, Vmat: np.ndarray, kmap: np.ndarray) -> tuple[float, np.ndarray]:
    """Max over q of Re eig(-chi0(q,0) V) (same quantity as ckchi in the SOC code) and its q."""
    mats = -np.einsum('nl,lmq->qnm', Vmat, chi0[:, :, 0, :])  # (V chi0)^T has the same spectrum
    st = np.linalg.eigvals(mats).real.max(axis=1)
    iq = int(st.argmax())
    return float(st[iq]), kmap[iq]


def gap_character_nosym(delta: np.ndarray, kmap: np.ndarray, slist: np.ndarray, Nx: int, Ny: int, Nz: int,
                        iQ: tuple = (0, 0, 0)) -> dict:
    """
    @brief Frequency parity and spin structure of a nosym gap (diagnostic; the solver projects onto neither).
    The w_n<0 half follows from Pauli, Delta_ab(k,-iw) = -Delta_ba(-k-Q,iw), so the even/odd-frequency
    parts are (Delta(iw) +- Delta(-iw))/2.  Without TRS/inversion, even and odd frequency can mix.
    An odd-frequency leading solution usually signals a channel the TRS solvers exclude by construction.
    @param delta: Gap function [Norb, Norb, Nw, Nkall] complex128 (python index [b,a] = Delta_ab)
    @param  kmap: Integer mesh coordinates [Nkall, 3] int64
    @param slist: Spin label of each orbital [Norb] int64
    @param Nx,Ny,Nz: k-grid dimensions
    @param    iQ: Pair momentum in integer mesh units
    @return  info: dict with relative weights 'even_freq', 'odd_freq', 'equal_spin', 'opposite_spin'
    """
    N = np.array([Nx, Ny, Nz])
    flat = {tuple(k): i for i, k in enumerate(kmap)}
    kmq = np.array([flat[tuple((-k - np.asarray(iQ)) % N)] for k in kmap])
    dneg = -np.swapaxes(delta, 0, 1)[:, :, :, kmq]    # Delta(-iw) in the same storage
    tot = np.sum(np.abs(delta)**2)
    even = np.sum(np.abs(0.5*(delta + dneg))**2)/tot
    odd = np.sum(np.abs(0.5*(delta - dneg))**2)/tot
    same = np.equal.outer(slist, slist)                 # [b,a]: same spin label
    eq = np.sum(np.abs(delta[same])**2)/tot
    return {'even_freq': float(even), 'odd_freq': float(odd), 'equal_spin': float(eq), 'opposite_spin': float(1.0 - eq)}


def calc_lin_eliash_nosym(Nx: int, Ny: int, Nz: int, Nw: int, ham_r: np.ndarray, S_r: np.ndarray,
                          rvec: np.ndarray, mu: float, temp: float, chiolist: np.ndarray, site: np.ndarray,
                          orb_dep: bool, U: float, J: float, gap_sym: int, sw_soc: bool,
                          slist: np.ndarray | None = None, invs: np.ndarray | None = None,
                          Umat: np.ndarray | None = None, Jmat: np.ndarray | None = None,
                          hvec: np.ndarray | list = (0.0, 0.0, 0.0), qlist: list | None = None,
                          init_delta: np.ndarray | None = None, eps: float = 1.0e-5, arnoldi_m: int = 60,
                          sw_out: bool = True, fill: float | None = None, mu_kmesh: int = 40) -> dict:
    """
    @brief Solve the linearized Eliashberg equation without symmetry (Zeeman field / FFLO) for a list of Q.

    The model is put on the full k-mesh in the spin-orbital basis (a spinless model is spin-doubled
    first), the Zeeman term h.sigma is added, and chi0 -> V_Delta is built ONCE (RPA, bare G).  The
    gap equation is then solved for each pair momentum Q in qlist; pairs are (k+Q,-k), Q is an integer
    mesh vector, Q=(iQ0/Nx, iQ1/Ny, iQ2/Nz).

    Side effects: 'lambda_Q.dat' gets one line per Q (iQ, Q in reciprocal-lattice units, lambda).
    With sw_out the gap of the Q with the largest lambda is saved to 'gap_nosym.npz' (orbital basis
    delta [Norb,Norb,Nw,Nkall], band basis deltab, kmap, iQ, slist, invs) and the kz=0 plane of its
    dominant orbital components at the lowest w_n is written to 'gap_nosym.dat'.

    The FLEX self-energy is not available without TRS: G is always the bare G0.

    Chemical potential: the Zeeman field moves the bands, so at a fixed electron count mu
    moves too.  With fill given and a non-zero field, mu is solved again for the Zeeman-split
    model, on the same mesh as get_mu (mu_kmesh^3) so that h=0 reproduces the mu passed in;
    without fill (or at h=0) the mu passed in is used as it is.

    @param Nx, Ny, Nz: k-mesh size
    @param      Nw: number of positive Matsubara frequencies
    @param   ham_r: hopping matrix H(R) [Nr, Norb, Norb] (spin-orbital if sw_soc, spinless otherwise)
    @param     S_r: overlap matrix S(R); [] for an orthogonal basis
    @param    rvec: R-vectors [Nr, 3]
    @param      mu: chemical potential [eV]
    @param    temp: temperature k_B*T [eV]
    @param chiolist: orbital-pair basis of chi [Nchi, 2] int64 (in the basis of ham_r)
    @param    site: site index of each pair [Nchi] int64
    @param orb_dep: True uses the orbital-dependent Umat/Jmat, False the scalar U/J
    @param    U, J: on-site Hubbard U and Hund J [eV]
    @param gap_sym: symmetry index of the seed (>=0 singlet, <0 Sz=0 triplet); no projection is applied
    @param  sw_soc: True if ham_r is already spin-orbital (then slist/invs are required)
    @param slist, invs: spin label and 1-based opposite-spin partner (sw_soc only)
    @param Umat, Jmat: orbital-resolved U and J in the basis of ham_r (orb_dep only)
    @param    hvec: Zeeman field (hx, hy, hz) [eV], H_Z = h.sigma
    @param   qlist: list of integer pair momenta [(iQx, iQy, iQz), ...]; None = [(0, 0, 0)]
    @param init_delta: explicit seed [Norb, Norb, Nkall] in the spin-orbital basis (overrides gap_sym)
    @param     eps: convergence tolerance
    @param arnoldi_m: Krylov dimension (0 = power method); small lambda may need >~60
    @param  sw_out: write the gap of the best Q
    @param    fill: band filling as in main.py (per spin for a spinless model, total for sw_soc);
                    None keeps mu fixed
    @param mu_kmesh: linear k-mesh of the mu solve (get_mu's default)
    @return results: dict with 'qlist', 'lambda', 'best_iQ', 'best_lambda', 'stoner', 'stoner_q',
                     'character' (gap_character_nosym of the best-Q gap), 'mu' (the one used)
    """
    if sw_soc:
        if slist is None or invs is None:
            raise ValueError('calc_lin_eliash_nosym: sw_soc=True needs slist and invs')
        ham_s, S_s = ham_r, S_r
        olist, site_s = chiolist, site
        Umat_s, Jmat_s = Umat, Jmat
    else:
        ham_s, S_s, slist, invs = spin_double_model(ham_r, S_r)
        olist, site_s = spin_double_chiolist(chiolist, site, ham_r.shape[1])
        Umat_s = None if Umat is None else np.tile(Umat, (2, 2))
        Jmat_s = None if Jmat is None else np.tile(Jmat, (2, 2))
    hvec = np.asarray(hvec, dtype=np.float64)
    print(f'nosym Eliashberg: Zeeman h = ({hvec[0]:.4e}, {hvec[1]:.4e}, {hvec[2]:.4e}) eV', flush=True)
    if np.any(hvec != 0.0):
        ham_s = add_zeeman(ham_s, S_s, rvec, slist, invs, hvec)
        if fill is not None:
            # same electron count in the Zeeman-split bands; the spin-doubled spinless model has
            # 2*fill electrons per cell (fill is per spin there), a spinful one fill
            nfill = fill if sw_soc else 2.0*fill
            Nkmu, eig_mu, _ = _bands.get_emesh(mu_kmesh, mu_kmesh, mu_kmesh, ham_s, S_s, rvec, np.eye(3))
            mu_h = _bands.calc_mu(eig_mu, Nkmu, nfill, temp)
            print(f'chem. pot. re-solved with the Zeeman field: {mu:.6f} -> {mu_h:.6f} eV '
                  f'(shift {mu_h-mu:+.3e} eV)', flush=True)
            mu = mu_h
    qlist = [(0, 0, 0)] if qlist is None else [tuple(int(x) for x in q) for q in qlist]

    # full-BZ mesh: kmap of gen_irr_k_TRS is the complete [0,N)^3 grid
    _, kmap, _ = flibs.gen_irr_k_TRS(Nx, Ny, Nz)
    klist = kmap / np.array([Nx, Ny, Nz], dtype=np.float64)
    eig, uni = _bands.get_eigs(klist, ham_s, S_s, rvec)
    Gk = flibs.gen_Green0(eig, uni, mu, temp, Nw)
    if orb_dep:
        Vmat = flibs.gen_Vmatrix_orb(olist, slist, site_s, invs, Umat_s, Jmat_s)
    else:
        Vmat = flibs.gen_Vmatrix(olist, slist, site_s, invs, U, J)

    chi = flibs.get_chi0_nosym(Gk, olist, kmap, temp, Nx, Ny, Nz)
    stoner, stoner_q = _stoner_nosym(chi, Vmat, kmap)
    print(f'Stoner factor = {stoner:.6f} at q = {stoner_q}', flush=True)
    if stoner >= 1.0:
        print('Warning: Stoner factor >= 1, the RPA normal state is magnetically unstable; '
              'lambda below is not meaningful', flush=True)
    Vdelta = flibs.get_vdelta_nosym(chi, Vmat, olist, slist)  # pairing channel; in place
    init = nosym_initial_gap(klist, slist, invs, gap_sym) if init_delta is None else init_delta

    lambdas, best = [], None
    for iQ in qlist:
        delta, lam = flibs.linearized_eliashberg_nosym(Vdelta, Gk, init, Vmat, olist, kmap, Nx, Ny, Nz, temp,
                                                       iQ=iQ, eps=eps, arnoldi_m=arnoldi_m)
        print(f'iQ = {iQ}: lambda_eliash = {lam:.6f}', flush=True)
        lambdas.append(lam)
        if best is None or lam > best[1]:
            best = (iQ, lam, delta)
    try:
        with open('lambda_Q.dat', 'w') as f:
            f.write(f'# iQx iQy iQz  Qx Qy Qz (r.l.u.)  lambda   h=({hvec[0]} {hvec[1]} {hvec[2]}) eV\n')
            for iQ, lam in zip(qlist, lambdas):
                f.write(f'{iQ[0]:4d} {iQ[1]:4d} {iQ[2]:4d} {iQ[0]/Nx:8.5f} {iQ[1]/Ny:8.5f} {iQ[2]/Nz:8.5f} '
                        f'{lam:14.8f}\n')
    except OSError as e:
        print(f"Error: Failed to write 'lambda_Q.dat': {e}", flush=True)
    print(f'best iQ = {best[0]}, lambda_eliash = {best[1]:.6f}', flush=True)
    char = gap_character_nosym(best[2], kmap, slist, Nx, Ny, Nz, best[0])
    print(f"gap character: even-freq {char['even_freq']:.4f}, odd-freq {char['odd_freq']:.4f}; "
          f"equal-spin {char['equal_spin']:.4f}, opposite-spin {char['opposite_spin']:.4f}", flush=True)
    if char['odd_freq'] > 0.5:
        print('Warning: the leading solution is mainly ODD-frequency (excluded by the TRS solvers); '
              'check the even-frequency candidates among the Ritz values', flush=True)
    if sw_out:
        _output_gap_nosym(best[2], uni, kmap, Nx, Ny, Nz, best[0], slist, invs)
    return {'qlist': qlist, 'lambda': np.array(lambdas), 'best_iQ': best[0], 'best_lambda': best[1],
            'stoner': stoner, 'stoner_q': stoner_q, 'character': char, 'mu': mu}


def _output_gap_nosym(delta: np.ndarray, uni: np.ndarray, kmap: np.ndarray, Nx: int, Ny: int, Nz: int,
                      iQ: tuple, slist: np.ndarray, invs: np.ndarray) -> None:
    """Save the gap of one Q (npz) and the kz=0 plane of its dominant orbital components (dat)."""
    deltab = flibs.conv_delta_orb_to_band_nosym(delta, uni, kmap, Nx, Ny, Nz, iQ)
    np.savez('gap_nosym', delta=delta, deltab=deltab, kmap=kmap, iQ=np.array(iQ), slist=slist, invs=invs)
    d0 = delta[:, :, 0, :]  # lowest w_n, [Norb, Norb, Nkall] (python index = Fortran transposed: d0[b,a]=D_ab)
    amp = np.abs(d0).max(axis=2)
    comps = [(a, b) for a in range(len(d0)) for b in range(len(d0)) if amp[b, a] > 1.0e-2*amp.max()]
    try:
        with open('gap_nosym.dat', 'w') as f:
            f.write(f'# iQ={tuple(iQ)}; kx ky then Re,Im of Delta_ab(k,iw_0) for (a,b) (1-based, spin-orbital): '
                    + ' '.join(f'({a+1},{b+1})' for a, b in comps) + '\n')
            for i, k in enumerate(kmap):
                if k[2] == 0:
                    vals = ' '.join(f'{d0[b, a, i].real:12.5e} {d0[b, a, i].imag:12.5e}' for a, b in comps)
                    f.write(f'{k[0]/Nx:7.4f} {k[1]/Ny:7.4f} {vals}\n')
    except OSError as e:
        print(f"Error: Failed to write 'gap_nosym.dat': {e}", flush=True)
