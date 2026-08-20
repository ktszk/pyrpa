#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
NMR observables in the superconducting state.

Everything here is a thin observable layer on top of the BdG spin susceptibility
already implemented in fchisc.f90 / _response.py.  Two quantities are produced:

  Knight shift          K(T)      ~  chi'_s(q=0, w=0)
  spin-lattice rate     1/(T1 T)  ~  <|A(q)|^2 chi''_s(q,w0)>_q / w0   (BZ average)

with w0 -> 0 (the NMR Larmor frequency is many orders below every electronic
scale, so only the slope of chi'' at the origin matters).

Both are also reported normalised to the NORMAL state evaluated at the same
temperature, k-mesh, q-mesh and broadening:

      K_s(T)/K_n(T)              -> Yosida function for a singlet gap
      (1/T1T)_s / (1/T1T)_n      -> Hebel-Slichter peak / power law

In that ratio the material prefactors (gamma_n, the hyperfine scale, g) cancel,
so what is left is purely the coherence-factor + gap-structure physics that the
BdG bubble encodes.  Normalising also buys a large amount of numerical robustness:
the two bubbles share their mesh errors, so the RATIO converges far faster than
either side alone.  On a small (pi,pi) electron pocket, refining 192 -> 1536 moved
the raw normal-state (1/T1T) over 1.5e-3 .. 1.2e-2 (a factor of 8) while the SC/normal
ratio settled at 0.24 .. 0.40.  Quote ratios; treat the raw R_sc, R_n columns of
nmr_sc.dat as diagnostics only.  The normal-state reference is computed through the
normal-state bubble (calc_chi in fchi.f90), which is the exact Delta->0 limit of
the SC bubble: at Delta=0 the BdG eigenvectors are block diagonal, the FF term
vanishes and GG reduces term by term to calc_chi.  Using it directly is 4x
cheaper (Norb instead of 2*Norb bands) and avoids the accidental
particle/hole eigenvector mixing that an exactly zero gap can produce in LAPACK.

Where the gap comes from
------------------------
Either source produces the same [Nk, Norb, Norb] static Delta(k) that nmr_sweep takes:

  * the phenomenological form factor (gap_symms / gap_sym + per-band delta0), built by
    main.py's SC-chi block -- quick symmetry scans;
  * a self-consistent RPA/FLEX/Eliashberg gap, via gap_from_eliashberg(), which reads the
    Delta(R, iw_n) written by output_gap_wannier.  Because that file is in real space it
    is mesh independent: the Eliashberg run can stay on the coarse k-mesh it needs to
    converge while the NMR sweep uses the far finer one that chi''(q,w->0) demands.

Either way the gap SHAPE is fixed once and only its amplitude is scaled with temperature
by the BCS interpolation.  Set the amplitude from the FERMI-SURFACE gap
(fermi_surface_gap), never from max|Delta| over the zone.

Sign / broadening notes
-----------------------
* The bubble uses weight = (f_l - f_m)/(w + E_m - E_l + i*idelta), i.e. the
  chi = -Pi_0 convention, so chi'' > 0 for w > 0 and no absolute value is taken.
* idelta is NOT a numerical detail here.  It plays the role of the
  quasi-particle lifetime that cuts off the BCS coherence peak in the DOS; with
  idelta -> 0 the Hebel-Slichter peak of an s-wave gap diverges logarithmically.
  Treat it as a physical (pair-breaking / smearing) parameter and check how the
  peak height moves with it.
* w0 must sit in the linear regime of chi''(w), i.e. w0 <~ idelta << Delta.
  The default w0 = 0.5*idelta satisfies this.
* The whole calculation only means anything inside the hierarchy
      dE_mesh << idelta << Delta_FS << W
  (k-mesh level spacing at E_F, broadening, the gap ON THE FERMI SURFACE --
  not max|Delta| over the zone -- and the bandwidth).  nmr_scale_check()
  is called automatically at the start of every sweep and warns when it is
  violated -- read its output before trusting a low-temperature number.  The
  usual failure is a k-mesh too coarse for the chosen idelta, which pins the
  T->0 value of (1/T1T)_s/(1/T1T)_n at a spurious plateau that does NOT go away
  when idelta is reduced.  See nmr_scale_check for the measured numbers.

Spin channel (sw_spsym)
-----------------------
sw_spsym flips the sign with which the anomalous (FF) term is added to the
normal (GG) term in fchisc.f90.  Its unambiguous operational meaning is:

  sw_spsym = False  (sgn=+1)  FF cancels GG as T->0  ->  Yosida suppression
  sw_spsym = True   (sgn=-1)  FF adds to GG          ->  response survives at T=0

Mapping that onto the standard triplet result chi_ij = chi_N[d_ij - dh_i dh_j (1-Y(T))]
(dh = unit d-vector, Y = Yosida function):

  * singlet: chi is isotropic and Yosida-suppressed in every direction, so
    sw_spsym=False is the only relevant setting and the same number drives K and 1/T1.
  * triplet: the component PARALLEL to d is Yosida-suppressed (sw_spsym=False), the
    components PERPENDICULAR to d are unaffected (sw_spsym=True).  For d || z that
    means field || z -> sw_spsym=False, field _|_ z -> sw_spsym=True.  Running the
    Knight shift with both settings gives exactly the K anisotropy that distinguishes
    a triplet (one axis stays flat through Tc) from a singlet (every axis drops).

NOTE: the comment in fchisc.f90 labels sgn=-1 as "triplet d||z: chi_zz preserved".
By the formula above the preserved component of a d||z triplet is the TRANSVERSE
one (chi_xx = chi_yy = chi_N), while chi_zz is Yosida-suppressed just like the
singlet.  The numerics are unaffected -- only that label is -- but keep it in mind
when reading the Fortran source.
"""
import numpy as np
import libs.flibs as flibs
from ._response import _prepare_sc_chi0_builder, _solve_chis_path

# BCS weak-coupling ratio 2*Delta(0)/(kB*Tc) = 3.528
BCS_GAP_RATIO = 1.764
# BCS interpolation coefficient in Delta(T)/Delta(0) = tanh(C*sqrt(Tc/T-1))
BCS_INTERP_C = 1.74
# Minimum number of k-mesh states that must lie within one broadening width at E_F for
# the Lorentzian to act as a delta function.  Calibrated in nmr_scale_check: at ~4 the
# Hebel-Slichter peak and the low-T power law are both absent, at ~9 both come out right.
N_DELTA_MIN = 8.0


def bcs_gap_ratio(temp: float, temp_c: float) -> float:
    """
    @fn bcs_gap_ratio
    @brief BCS interpolation formula for the gap amplitude Delta(T)/Delta(0).
    Delta(T)/Delta(0) = tanh(1.74*sqrt(Tc/T - 1)); exact at T->0 and reproduces
    the mean-field sqrt(1-T/Tc) onset at Tc to better than 1% over the whole range.
    The gap SHAPE (gap_sym / per-band delta0) is untouched: only the overall
    amplitude is scaled.
    @param   temp: temperature in eV
    @param temp_c: critical temperature in eV
    @return ratio: Delta(T)/Delta(0), 0 for temp >= temp_c
    """
    if temp >= temp_c or temp_c <= 0.0:
        return 0.0
    if temp <= 0.0:
        return 1.0
    return float(np.tanh(BCS_INTERP_C * np.sqrt(temp_c / temp - 1.0)))


def bcs_temp_c(delta_max: float) -> float:
    """
    @fn bcs_temp_c
    @brief Weak-coupling estimate of Tc from the zero-temperature gap: kB*Tc = Delta(0)/1.764.
    @param delta_max: maximum gap amplitude Delta(0) in eV
    @return   temp_c: critical temperature in eV
    """
    return abs(delta_max) / BCS_GAP_RATIO


def gap_from_eliashberg(klist: np.ndarray, fname: str = 'gap_wannier',
                        sw_extrapolate: bool = True, iw_index: int = 0, n_avg: int = 1,
                        n_points: int = 4, order: int = 1, eig: np.ndarray | None = None,
                        mu: float | None = None, gap_max: float | None = None,
                        verbose: bool = True) -> np.ndarray:
    """
    @fn gap_from_eliashberg
    @brief Load a self-consistent RPA/FLEX/Eliashberg gap and return Delta(k) for the BdG bubble.

    Reads the Delta(R, iw_n) written by output_gap_wannier (LIN_ELIASHBERG /
    NONLIN_ELIASHBERG with sw_out_self=True), reduces it to a STATIC orbital-basis
    Delta(k) and evaluates it on the supplied k-mesh -- the drop-in replacement for the
    phenomenological gap_symms/delta0 form factor used elsewhere in this module.

    Why the real-space file and not gap.npy: Delta(R) is mesh independent, so the
    Eliashberg run can stay on the coarse k-mesh it needs to converge while the NMR sweep
    runs on the much finer mesh that chi''(q,w->0) demands (see nmr_scale_check).  Reading
    gap.npy instead would force both to share a mesh.

    Matsubara -> static.  The BdG Hamiltonian takes a static gap, but Eliashberg produces
    Delta(k, iw_n) and w_n = (2n+1) pi T is never zero.  With sw_extrapolate (default) the
    w_n -> 0 limit is obtained by the same even-in-w_n polynomial fit used by
    output_gap_function (gap_extrapolate_w0), which removes the O((pi T)^2) bias of simply
    reading the lowest slice.  Doing the fit here in R-space is exactly equivalent to
    doing it in k-space: both the fit and the Fourier transform are linear in the data.
    Set sw_extrapolate=False to take slice iw_index (averaged over n_avg slices) instead,
    which is what the Eilenberger path does via eil_gap_iw / eil_gap_navg.

    Amplitude.  The LINEARIZED Eliashberg equation returns an eigenvector, so its overall
    scale is arbitrary and gap_max MUST be supplied to set it (the shape, signs, nodes and
    orbital content are what the calculation actually determines).  A gap from the
    NONLINEAR equation already carries a physical amplitude; pass gap_max=None to keep it.

    @param   klist: k-points to evaluate on, fractional [Nk, 3]
    @param   fname: base name of the .npz from output_gap_wannier (no extension)
    @param sw_extrapolate: fit Delta(iw_n) -> Delta(0) instead of reading one slice
    @param iw_index,n_avg: slice / slice-average to use when sw_extrapolate is False
    @param n_points,order: polynomial fit controls passed to gap_extrapolate_w0
    @param eig,mu: normal-state bands and chemical potential; when both are given the
                   rescaling target below is the maximum ON THE FERMI SURFACE rather than
                   over the whole zone, which is the only meaningful normalisation for a
                   gap whose form factor is anisotropic (see fermi_surface_gap)
    @param gap_max: rescale so that the largest gap equals this value in eV;
                    None = use the stored amplitude unchanged
    @param verbose: report what was loaded and how it was rescaled
    @return delta_k: static gap in orbital basis [Nk, Norb, Norb] complex128, laid out
                     exactly like the deltak built by the gap_sym path in main.py
    """
    from ._calc import gap_extrapolate_w0
    d = np.load(f'{fname}.npz')
    gap_R, rvec = d['gap'], d['rvec']                 # (Norb,Norb,N_iw,Nr), (Nr,3)
    temp_gap = float(d['temp']) if 'temp' in d else None
    Norb, N_iw, Nr = gap_R.shape[0], gap_R.shape[2], gap_R.shape[3]

    if sw_extrapolate:
        if temp_gap is None:
            raise ValueError(f"{fname}.npz has no 'temp'; cannot build the Matsubara mesh "
                             "for the w_n->0 fit.  Use sw_extrapolate=False.")
        if N_iw < n_points:
            raise ValueError(f"{fname}.npz holds only {N_iw} Matsubara slices, "
                             f"need >= n_points={n_points} for the w_n->0 fit")
        gR0, bias = gap_extrapolate_w0(gap_R, temp_gap, n_points, order)
        if verbose:
            print(f'  gap w_n->0 extrapolation: max correction over iw_0 ='
                  f' {bias.max():.3e} of the peak gap', flush=True)
    else:
        n1 = min(iw_index + n_avg, N_iw)
        if iw_index >= N_iw:
            raise ValueError(f"iw_index={iw_index} out of range (N_iw={N_iw})")
        gR0 = gap_R[:, :, iw_index:n1, :].mean(axis=2)                  # (Norb,Norb,Nr)

    # Delta(k) = sum_R Delta(R) exp(-2 pi i k.R).  output_gap_wannier builds Delta(R) with
    # np.fft.ifftn over the k-axes (kernel exp(+2 pi i k.R)), so this MINUS sign is its
    # exact inverse -- verified to round-trip a random gap to 3e-15.  Same convention as
    # gap_orbital_from_wannier (_eilenberger.py) and the self-energy interpolation in
    # main.py's plot_spectrum.
    kl = np.ascontiguousarray(klist, dtype=np.float64)
    rv = np.ascontiguousarray(rvec, dtype=np.float64)
    Nk = len(kl)
    delta_k = np.empty((Nk, Norb, Norb), dtype=np.complex128)
    # Cap the phase matrix at ~4M complex entries (~64 MB); Nk*Nr would otherwise be
    # gigabytes on the fine meshes this module needs.
    chunk = max(1, int(4e6 // max(Nr, 1)))
    for i0 in range(0, Nk, chunk):
        ph = np.exp(-2j * np.pi * (kl[i0:i0 + chunk] @ rv.T))            # [nk, Nr]
        # [k,a,b] = Delta_ab(k): the same layout as the gap_sym path's deltaini_static
        delta_k[i0:i0 + chunk] = np.einsum('kr,abr->kab', ph, gR0, optimize=True)

    raw_max = float(np.abs(delta_k).max())
    if verbose:
        print(f'  loaded {fname}.npz: Norb={Norb}, {Nr} R-vectors, {N_iw} Matsubara slices'
              + (f', T={temp_gap:.5f} eV' if temp_gap is not None else ''), flush=True)
    if gap_max is not None:
        if eig is not None and mu is not None:
            _, ref, _ = fermi_surface_gap(eig, mu, delta_k)
            where = 'on the Fermi surface'
        else:
            ref, where = raw_max, 'over the zone'
        if ref <= 0.0:
            raise ValueError('gap_from_eliashberg: the loaded gap vanishes '
                             f'{where}; cannot rescale to gap_max')
        delta_k *= gap_max / ref
        if verbose:
            print(f'  rescaled: max|Delta| {where} set to {gap_max*1e3:.4f} meV'
                  f' (raw zone max was {raw_max:.4e} eV)', flush=True)
    elif verbose:
        print(f'  amplitude used as stored (zone max {raw_max*1e3:.4f} meV).'
              '  A LINEARIZED-Eliashberg gap is an eigenvector with an arbitrary scale:'
              ' set gap_max unless this came from the nonlinear solver.', flush=True)
    return delta_k


def nmr_qmesh(Nx: int, Ny: int, Nz: int, Nqx: int, Nqy: int, Nqz: int = 1,
              sw_fold: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """
    @fn nmr_qmesh
    @brief Build the q-mesh for the 1/T1 Brillouin-zone sum, commensurate with the k-mesh.

    get_qshift resolves k+q by an exact hash lookup on the k-mesh, so every q must be
    a multiple of the k-mesh spacing.  The q-mesh is therefore a stride sub-sampling of
    the k-mesh: q = (jx*sx/Nx, jy*sy/Ny, jz*sz/Nz) with sa = Na//Nqa.  Gamma is index 0,
    which is what the Knight-shift extraction relies on.

    Cost note: one bubble costs O(Nk), so the whole sweep costs O(NT * Nq * Nk).  Setting
    the q-mesh equal to the k-mesh makes it O(Nk^2) and is not tractable beyond toy meshes
    (measured ~2.7e-4 s per k-point per bubble at Norb=5, Nchi=25, 24 threads).
    The two meshes have different jobs and need not match:
      * the k-sum is a RESOLUTION problem - it must resolve the Fermi surface, the gap
        nodes and the broadening, i.e. idelta >~ the k-mesh energy spacing v_F*dk;
      * the q-sum is an INTEGRATION of a function that is smooth in q, converging as
        O(Nq^-2).  8-16 points per axis is usually enough.  The exception is a system
        close to a magnetic instability (Stoner factor -> 1), where chi_s(q) develops a
        peak of width ~1/xi at q_AF and Nq per axis must exceed L/xi.
    Always confirm with nmr_qconvergence() rather than assuming.

    @param  Nx,Ny,Nz: k-mesh dimensions
    @param Nqx,Nqy,Nqz: requested q-mesh dimensions (rounded down to a divisor stride)
    @param   sw_fold: fold q and -q together using chi(q)=chi(-q).  Exact for the
                      centrosymmetric + time-reversal-symmetric systems that mkBdGhamk
                      already assumes, and halves the cost.
    @retval    qlist: q-points in fractional coordinates [Nq, 3] float64 (Gamma first)
    @retval    qmult: multiplicity of each retained q-point [Nq] float64
    """
    def axis(N, Nq):
        # The stride must DIVIDE N: otherwise the last interval of the sub-mesh is short
        # (biasing the BZ average) and the set is not closed under q -> -q, which would
        # break the folding below.  Round the request down to the nearest divisor of N.
        N, Nq = int(N), max(1, min(int(Nq), int(N)))
        n_use = max(d for d in range(1, N + 1) if N % d == 0 and d <= Nq)
        return np.arange(0, N, N // n_use, dtype=np.float64) / float(N)
    qx, qy, qz = axis(Nx, Nqx), axis(Ny, Nqy), axis(Nz, Nqz)
    gx, gy, gz = np.meshgrid(qx, qy, qz, indexing='ij')
    qlist = np.ascontiguousarray(np.array([gx.ravel(), gy.ravel(), gz.ravel()]).T, dtype=np.float64)
    if not sw_fold:
        return qlist, np.ones(len(qlist), dtype=np.float64)
    # Fold q against -q (mod 1).  Quantise to integers first so the pairing is exact
    # and independent of floating-point round-off in the stride division.
    N = np.array([Nx, Ny, Nz], dtype=np.int64)
    iq = np.rint(qlist * N).astype(np.int64) % N
    ineg = (-iq) % N
    keep, mult = [], []
    seen = set()
    for j in range(len(iq)):
        key, nkey = tuple(iq[j]), tuple(ineg[j])
        if key in seen:
            continue
        seen.add(key)
        if nkey == key:
            keep.append(j); mult.append(1.0)      # self-paired (Gamma, zone corners)
        else:
            seen.add(nkey)
            keep.append(j); mult.append(2.0)      # stands in for its -q partner too
    return (np.ascontiguousarray(qlist[keep], dtype=np.float64),
            np.array(mult, dtype=np.float64))


def hyperfine_form_factor(qlist: np.ndarray, hf_A: float = 1.0, hf_B: float = 0.0) -> np.ndarray:
    """
    @fn hyperfine_form_factor
    @brief Squared hyperfine form factor |A(q)|^2 entering the 1/T1 momentum sum.
    Uses the standard on-site + transferred (4 in-plane neighbours) model
        A(q) = hf_A + 2*hf_B*(cos(2 pi qx) + cos(2 pi qy)),
    which is the Mila-Rice form used for planar nuclei.  hf_B=0 gives the
    q-independent weight A(q)=hf_A (all q counted equally); hf_A = -4*hf_B
    filters out the antiferromagnetic q=(pi,pi) response, the usual way to
    contrast a plane site with an apical/ligand site.
    @param  qlist: q-points in fractional coordinates [Nq, 3]
    @param   hf_A: on-site hyperfine coupling (arbitrary units, only ratios matter)
    @param   hf_B: transferred hyperfine coupling to the 4 in-plane neighbours
    @return weight: |A(q)|^2 [Nq] float64
    """
    q = np.asarray(qlist, dtype=np.float64)
    Aq = hf_A + 2.0 * hf_B * (np.cos(2 * np.pi * q[:, 0]) + np.cos(2 * np.pi * q[:, 1]))
    return Aq * Aq


def _nmr_observables(chisq: np.ndarray, weight: np.ndarray, w0: float) -> tuple[float, float]:
    """
    Reduce the q-resolved susceptibility [Nq, 2] (w = 0 and w = w0) to the two
    NMR observables.  chisq[:,0] is the static response (real by construction in
    the w=0 branch of the bubble), chisq[:,1] the finite-w one.
    @return (K, R): K = chi'_s(q=0,w=0)
                    R = sum_q w(q) chi''_s(q,w0) / (w0 * sum_q w(q))   ( ~ 1/(T1 T) ),
                        w(q) = folding multiplicity * |A(q)|^2
    """
    K = float(chisq[0, 0].real)
    R = float(np.sum(weight * chisq[:, 1].imag) / (np.sum(weight) * w0))
    return K, R


def fermi_surface_gap(eig: np.ndarray, mu: float, delta_k: np.ndarray,
                      window: float | None = None) -> tuple[float, float]:
    """
    @fn fermi_surface_gap
    @brief Gap amplitude ON the Fermi surface, and the fraction of it that is left ungapped.

    max|Delta(k)| over the whole zone is the wrong scale for every NMR diagnostic: what
    controls the thermodynamics and the relaxation is the gap where there are states.
    A d_x2-y2 form factor on a pocket centred at (pi,pi), for instance, is ~0 all over
    that pocket while its zone-wide maximum is full size.
    @param     eig: normal-state eigenvalues [Nk, Norb] float64
    @param      mu: chemical potential in eV
    @param delta_k: gap in orbital basis [Nk, Norb, Norb] complex128
    @param  window: half-width in eV of the shell around E_F; default 2% of the bandwidth
    @retval gap_fs: r.m.s. of ||Delta(k)||_F/sqrt(Norb) over k-points with a band in the
                    shell -- the scale to compare against idelta and against T
    @retval gap_fs_max: largest gap on that shell (the antinodal value for a nodal gap).
                    THIS is the BCS scale that sets Tc, not the zone-wide max|Delta|
    @retval  f_low: fraction of those k-points whose gap is below gap_fs/10 (nodal weight)
    """
    e = np.asarray(eig, dtype=np.float64) - mu
    Norb = e.shape[1]
    if window is None:
        window = 0.02 * float(e.max() - e.min())
    for _ in range(8):
        on_fs = np.any(np.abs(e) < window, axis=1)
        if int(np.count_nonzero(on_fs)) >= 50:
            break
        window *= 2.0
    if not np.any(on_fs):
        return 0.0, 0.0, 1.0
    d = np.sqrt(np.sum(np.abs(np.asarray(delta_k)[on_fs]) ** 2, axis=(1, 2)) / Norb)
    gap_fs = float(np.sqrt(np.mean(d ** 2)))
    gap_fs_max = float(d.max())
    f_low = float(np.count_nonzero(d < 0.1 * gap_fs) / len(d)) if gap_fs > 0 else 1.0
    return gap_fs, gap_fs_max, f_low


def nmr_scale_check(eig: np.ndarray, mu: float, idelta: float, gap0: float,
                    delta_k: np.ndarray | None = None, verbose: bool = True) -> dict:
    """
    @fn nmr_scale_check
    @brief Verify the energy-scale hierarchy that 1/T1T from a real-frequency bubble requires.

    chi''(q,w0)/w0 is the most delicate quantity in this module, because the broadening
    idelta is a LORENTZIAN: its tails fall off only as 1/w^2, so every particle-hole pair
    in the band, however far from the Fermi level, leaks a finite amount of spectral
    weight down to w=0.  That leakage is what fills a fully gapped SC state in, and it
    only becomes negligible relative to the genuine quasi-elastic response when

        dE_mesh   <<   idelta   <<   Delta(0)   <<   W

    where dE_mesh is the mean level spacing of the k-mesh at E_F.

    The left inequality is the one that is usually violated, and violating it is silent.
    The sharp criterion is not dE_mesh < idelta but the NUMBER of k-mesh states inside one
    broadening width, n_delta = 2*idelta/dE_mesh: the Lorentzian only acts as a delta
    function once it contains ~10 levels.  Below that, both the normal and the SC bubble
    are pure Lorentzian tail, both scale linearly in idelta, and their RATIO saturates at
    a finite, meaningless value that does NOT improve when idelta is reduced.

    Measured on a single-orbital square lattice (t=1, t'=0.4, mu=-0.6, Delta(0)=0.05 eV),
    where the correct answers are a Hebel-Slichter peak of ~3.4 for the s-wave gap and
    1/T1 ~ T^3 for the d-wave one:

        N      dE_mesh   idelta    n_delta   HS peak   d-wave 1/T1 exponent
        96     9.0e-4    2.0e-3      4.4      1.05          0.9   (both wrong)
        128    5.2e-4    1.0e-3      3.9      1.12          1.1   (both wrong)
        192    2.2e-4    1.0e-3      9.0      3.39          2.9   (converged)
        256    1.2e-4    1.0e-3     16.0      3.44          3.1   (converged)

    At fixed N=128 the low-temperature ratio sticks at 0.085 for every idelta from 8e-3
    down to 5e-4: refining the K-MESH, not the broadening, is the fix.

    N_DELTA_MIN = 8 is a FLOOR calibrated on a large Fermi surface, not a guarantee.  A
    small Fermi pocket needs much more, because the normal-state chi'' then vanishes over
    most of the q-mesh (no states to scatter between at large q) while the SC-state
    Lorentzian tail does not, so the ratio is inflated by whatever tail survives.  Same
    single-orbital lattice, mu=2.83 (a small pocket at (pi,pi)), d-wave gap, T/Tc=0.3:

        N       idelta    n_delta   (1/T1T)_s/(1/T1T)_n
        192     1.0e-3      7.0            1.62
        384     1.0e-3     27.7            0.243
        768     5.0e-4     55.5            0.405
        1536    1.25e-4    55.4            0.381

    n_delta=7 overshoots the converged answer by a factor of four.  When the diagnostic
    only just passes, refine once more and confirm the ratio actually stops moving.

    The right inequality matters too: idelta >~ 0.2*Delta(0) smears the coherence peak of
    the BCS density of states and both flattens the Hebel-Slichter peak and leaks weight
    into the gap.

    @param     eig: normal-state eigenvalues [Nk, Norb] float64
    @param      mu: chemical potential in eV
    @param  idelta: Lorentzian broadening in eV
    @param    gap0: zero-temperature gap amplitude max|Delta(0)| in eV
    @param delta_k: gap in orbital basis [Nk, Norb, Norb]; when given, the comparison
                    against idelta uses the FERMI-SURFACE gap from fermi_surface_gap()
                    instead of the zone-wide maximum, which is the only version that
                    means anything for a nodal or strongly anisotropic gap
    @param verbose: print the hierarchy and any violation
    @retval  scales: dict with 'dE_mesh','n_delta','idelta','gap0','gap_fs','f_ungapped',
                     'bandwidth','ok'
    """
    e = np.asarray(eig, dtype=np.float64).ravel() - mu
    band_w = float(e.max() - e.min())
    # Mean level spacing at E_F: count states in a window wide enough to hold a decent
    # sample but still local, then divide the window by the count.
    win = max(10.0 * idelta, 0.02 * band_w)
    for _ in range(8):
        n_in = int(np.count_nonzero(np.abs(e) < win))
        if n_in >= 50:
            break
        win *= 2.0
    dE_mesh = 2.0 * win / max(n_in, 1)
    n_delta = 2.0 * idelta / dE_mesh              # k-mesh states inside one broadening
    if delta_k is not None:
        gap_fs, gap_fs_max, f_low = fermi_surface_gap(eig, mu, delta_k)
    else:
        gap_fs, gap_fs_max, f_low = abs(gap0), abs(gap0), float('nan')
    fine_enough = n_delta >= N_DELTA_MIN
    not_smeared = (idelta < 0.2 * gap_fs) if gap_fs > 0 else True
    if verbose:
        print(f'  energy scales: dE_mesh={dE_mesh:9.3e} < idelta={idelta:9.3e}'
              f' < Delta_FS={gap_fs:9.3e} (max {abs(gap0):9.3e}) < W={band_w:9.3e} eV'
              f'  (states within idelta: {n_delta:.1f})', flush=True)
        if not fine_enough:
            print(f'  WARNING: only {n_delta:.1f} k-mesh states lie within idelta'
                  f' (need >~{N_DELTA_MIN:.0f}).\n'
                  '           The Lorentzian is not yet acting as a delta function, so chi\'\' is'
                  ' mostly\n'
                  '           tail in BOTH the normal and the SC state and their ratio saturates at a'
                  '\n           spurious value: no Hebel-Slichter peak, no low-T power law.'
                  '\n           Refine the k-mesh (Nx,Ny,Nz) -- changing idelta does NOT fix this.',
                  flush=True)
        if not not_smeared:
            print(f'  WARNING: idelta ({idelta:.3e}) is not small against the Fermi-surface gap'
                  f' ({gap_fs:.3e}).'
                  '\n           The coherence peak is over-smeared: the Hebel-Slichter peak will be'
                  '\n           suppressed and weight leaks into the gap.', flush=True)
        if delta_k is not None and gap_fs_max < 0.5 * abs(gap0):
            print(f'  WARNING: the largest gap ON the Fermi surface ({gap_fs_max:.3e}) is well'
                  f' below the\n           zone-wide maximum ({abs(gap0):.3e}): the gap form'
                  ' factor has little weight where\n           the states are.  Tc is estimated'
                  ' from the Fermi-surface value, but check that the\n           gap symmetry is'
                  ' the one you meant for this Fermi surface.', flush=True)
    return {'dE_mesh': dE_mesh, 'n_delta': n_delta, 'idelta': idelta, 'gap0': abs(gap0),
            'gap_fs': gap_fs, 'gap_fs_max': gap_fs_max, 'f_ungapped': f_low,
            'bandwidth': band_w, 'ok': bool(fine_enough and not_smeared)}


def nmr_cost_estimate(Nk: int, Nq: int, NT: int, Norb: int, Nchi: int,
                      verbose: bool = True) -> float:
    """
    @fn nmr_cost_estimate
    @brief Rough wall-clock estimate of a full NMR temperature sweep, printed before the run.
    One SC bubble costs O(Nk * (2*Norb)^2 * Nchi^2) and was measured at 1.104 s for
    Nk=4096, Norb=5, Nchi=25 with 2 frequency points on 24 threads; the estimate scales
    that reference point.  The normal-state reference adds a further ~1/4 (Norb instead
    of 2*Norb bands).  Absolute numbers depend on the machine, but the SCALING is the
    point: the sweep is linear in Nq, so a q-mesh equal to the k-mesh is quadratic in Nk.
    @param    Nk: number of k-points
    @param    Nq: number of q-points (after folding)
    @param    NT: number of temperatures
    @param  Norb: number of orbitals
    @param  Nchi: number of orbital pairs
    @return t_est: estimated seconds
    """
    ref = 1.104 / (4096 * (2 * 5) ** 2 * 25 ** 2)     # s per (k-point * band-pair * Nchi^2)
    t_bubble = ref * Nk * (2 * Norb) ** 2 * Nchi ** 2
    t_est = NT * Nq * 1.25 * t_bubble
    if verbose:
        print(f'  estimated cost: {Nq} q-points x {NT} temperatures x {Nk} k-points'
              f'  ~ {t_est:.0f} s ({t_est/3600:.2f} h)', flush=True)
        if t_est > 6 * 3600:
            print('  WARNING: over 6 h.  Cost is linear in Nq; halve the q-mesh and check',
                  '\n           with nmr_qconvergence() whether the coarser sum is already converged.',
                  flush=True)
    return t_est


def nmr_qconvergence(hamk: np.ndarray, delta_k0: np.ndarray, klist: np.ndarray,
                     Nx: int, Ny: int, Nz: int, qsize, mu: float, temp: float,
                     temp_c: float, Smat: np.ndarray, olist: np.ndarray,
                     idelta: float = 1.e-3, w0: float | None = None,
                     sw_spsym: bool = False, hf_A: float = 1.0, hf_B: float = 0.0,
                     verbose: bool = True) -> dict:
    """
    @fn nmr_qconvergence
    @brief Check the q-mesh convergence of the 1/T1T sum at a single temperature.
    Evaluates the sum on the requested q-mesh and on the mesh with half the points per
    axis, and reports the relative change.  Run this once (deep in the SC state, where
    chi''(q) is most structured) before committing to a long temperature sweep: if the
    change is already at the percent level, the coarse mesh is enough and the sweep cost
    drops by 2^(number of refined axes).

    How much q-mesh you need depends strongly on the gap.  A full gap is forgiving: the
    Hebel-Slichter peak is a Fermi-surface average and 8 points per axis already give it.
    A NODAL gap is not: at low temperature the only quasi-particles sit at the nodes, so
    chi''(q,w->0) has weight on the few discrete q that connect one node to another, and
    a uniform mesh samples them erratically.  Measured on a single-orbital square lattice
    with a d_x2-y2 gap (192x192 k-mesh, idelta=1e-3), the fitted exponent of 1/T1 ~ T^n
    over T/Tc = 0.15..0.45 came out as

        Nq per axis   4     6     8    12    16    24
        n            1.4   1.9   2.9   2.3   2.9   2.5

    against the analytic n=3.  The qualitative power law is robust from ~8 per axis on,
    but the exponent itself is only good to about +-0.5 -- quote it with that caveat, or
    push Nq much higher for the low-temperature points only.
    @param ...: as nmr_sweep, plus the k-mesh dimensions and the requested q-mesh size
    @param  qsize: requested [Nqx, Nqy, Nqz]
    @retval result: dict with 'R_fine','R_coarse','rel_diff','Nq_fine','Nq_coarse'
    """
    out = {}
    for tag, size in (('fine', qsize),
                      ('coarse', [max(1, int(s) // 2) for s in qsize])):
        qlist, qmult = nmr_qmesh(Nx, Ny, Nz, *size)
        res = nmr_sweep(hamk, delta_k0, klist, qlist, mu, [temp], temp_c, Smat, olist,
                        idelta=idelta, w0=w0, sw_spsym=sw_spsym,
                        hf_weight=qmult * hyperfine_form_factor(qlist, hf_A, hf_B),
                        verbose=False)
        out['R_' + tag] = float(res['R_ratio'][0])
        out['Nq_' + tag] = len(qlist)
    denom = abs(out['R_fine']) if out['R_fine'] != 0 else 1.0
    out['rel_diff'] = abs(out['R_fine'] - out['R_coarse']) / denom
    if verbose:
        print(f'  q-mesh convergence at T={temp:.6f} eV (T/Tc={temp/temp_c:.3f}):', flush=True)
        print(f'    Nq={out["Nq_fine"]:5d} -> (1/T1T)ratio = {out["R_fine"]:.6f}', flush=True)
        print(f'    Nq={out["Nq_coarse"]:5d} -> (1/T1T)ratio = {out["R_coarse"]:.6f}', flush=True)
        print(f'    relative change = {out["rel_diff"]*100:.2f} %', flush=True)
    return out


def nmr_sweep(hamk: np.ndarray, delta_k0: np.ndarray, klist: np.ndarray, qlist: np.ndarray,
              mu: float, temp_list: np.ndarray, temp_c: float, Smat: np.ndarray,
              olist: np.ndarray, idelta: float = 1.e-3, w0: float | None = None,
              sw_spsym: bool = False, hf_weight: np.ndarray | None = None,
              verbose: bool = True, sw_scale_check: bool = True) -> dict:
    """
    @fn nmr_sweep
    @brief Temperature sweep of the SC-state NMR observables (Knight shift and 1/T1T).
    For every temperature the gap is rescaled by the BCS interpolation
    Delta(T)=Delta(0)*tanh(1.74*sqrt(Tc/T-1)), the BdG Hamiltonian is diagonalised
    ONCE, and the q-loop reuses that diagonalisation.  The normal-state reference
    is evaluated at the same temperature with the normal-state bubble.
    @param     hamk: normal-state k-space Hamiltonian [Nk, Norb, Norb] complex128
    @param delta_k0: zero-temperature gap in orbital basis [Nk, Norb, Norb] complex128
    @param    klist: k-point list, periodic [0,1) mesh without endpoints [Nk, 3]
    @param    qlist: q-mesh for the 1/T1 sum; qlist[0] MUST be Gamma [Nq, 3]
    @param       mu: chemical potential in eV
    @param temp_list: temperatures to evaluate, in eV [NT]
    @param   temp_c: critical temperature in eV (sets Delta(T))
    @param     Smat: Stoner interaction matrix [Nchi, Nchi] float64
    @param    olist: orbital index pairs for chi [Nchi, 2] int64 (chiolist)
    @param   idelta: Lorentzian broadening in eV; acts as the quasi-particle
                     lifetime that cuts off the s-wave coherence peak
    @param       w0: NMR probe frequency in eV; default 0.5*idelta (must be <~ idelta)
    @param sw_spsym: False -> Yosida-suppressed channel (singlet, or field || d for a
                     triplet); True -> channel preserved at T=0 (field _|_ d).
                     See the module docstring for the mapping.
    @param hf_weight: total q-weight [Nq]; normally qmult * |A(q)|^2, i.e. the folding
                      multiplicity from nmr_qmesh times hyperfine_form_factor.
                      None -> uniform weight
    @param  verbose: print per-temperature progress
    @param sw_scale_check: run nmr_scale_check on the normal-state bands first.
                     Set False when the caller has already reported it (main.py does,
                     so the warning appears before the q-convergence probe, not after)
    @retval  result: dict with keys
                     'temp'      temperatures [NT] (eV)
                     't_red'     reduced temperature T/Tc [NT]
                     'gap_ratio' Delta(T)/Delta(0) [NT]
                     'K_sc','K_n'        Knight shift, SC and normal [NT]
                     'R_sc','R_n'        1/(T1 T), SC and normal [NT]
                     'K_ratio','R_ratio' K_s/K_n and (1/T1T)_s/(1/T1T)_n [NT]
                     'temp_c','w0','idelta','sw_spsym'
    """
    qlist = np.ascontiguousarray(qlist, dtype=np.float64)
    if not np.allclose(qlist[0], 0.0):
        raise ValueError('nmr_sweep: qlist[0] must be Gamma (0,0,0) for the Knight shift')
    if w0 is None:
        w0 = 0.5 * idelta
    if w0 <= 0.0:
        raise ValueError('nmr_sweep: w0 must be positive (chi\'\'/w0 is the w->0 slope)')
    wlist = np.array([0.0, w0], dtype=np.float64)
    weight = np.ones(len(qlist)) if hf_weight is None else np.asarray(hf_weight, dtype=np.float64)

    # Normal-state bands are temperature independent: diagonalise once.
    eig_n, uni_n = flibs.get_eig(hamk)
    if verbose and sw_scale_check:
        nmr_scale_check(eig_n, mu, idelta, float(np.abs(delta_k0).max()), delta_k0)

    temp_list = np.asarray(temp_list, dtype=np.float64)
    out = {k: [] for k in ('gap_ratio', 'K_sc', 'K_n', 'R_sc', 'R_n')}
    for it, temp in enumerate(temp_list):
        ratio = bcs_gap_ratio(temp, temp_c)

        # ---- normal state (same mesh / broadening / w0) --------------------
        ffermi_n = flibs.get_ffermi(eig_n, mu, temp)

        def chi0_normal(q):
            qshift = flibs.get_qshift(klist, q)
            return flibs.get_chi_irr(uni_n, eig_n, ffermi_n, qshift, olist, wlist, idelta, temp)

        chisq_n, _ = _solve_chis_path(qlist, wlist, olist, Smat, chi0_normal)
        K_n, R_n = _nmr_observables(chisq_n, weight, w0)

        # ---- superconducting state ----------------------------------------
        if ratio <= 0.0:
            # T >= Tc: the gap is identically zero, so the SC bubble reduces to the
            # normal one exactly.  Skip it rather than diagonalising a gapless BdG
            # matrix, whose degenerate particle/hole pairs LAPACK may mix arbitrarily.
            K_sc, R_sc = K_n, R_n
        else:
            builder = _prepare_sc_chi0_builder(hamk, ratio * delta_k0, mu, temp,
                                               klist, olist, wlist, idelta, sw_spsym)
            chisq_sc, _ = _solve_chis_path(qlist, wlist, olist, Smat, builder)
            K_sc, R_sc = _nmr_observables(chisq_sc, weight, w0)

        for key, val in (('gap_ratio', ratio), ('K_sc', K_sc), ('K_n', K_n),
                         ('R_sc', R_sc), ('R_n', R_n)):
            out[key].append(val)
        if verbose:
            print(f'  [{it+1:3d}/{len(temp_list):3d}] T={temp:9.6f} eV  T/Tc={temp/temp_c:6.3f}'
                  f'  D/D0={ratio:5.3f}  K/Kn={K_sc/K_n if K_n!=0 else np.nan:8.5f}'
                  f'  (1/T1T)/(1/T1T)n={R_sc/R_n if R_n!=0 else np.nan:8.5f}', flush=True)

    res = {k: np.array(v, dtype=np.float64) for k, v in out.items()}
    res['temp'] = temp_list
    res['t_red'] = temp_list / temp_c
    with np.errstate(divide='ignore', invalid='ignore'):
        res['K_ratio'] = np.where(res['K_n'] != 0, res['K_sc'] / res['K_n'], np.nan)
        res['R_ratio'] = np.where(res['R_n'] != 0, res['R_sc'] / res['R_n'], np.nan)
    res.update({'temp_c': temp_c, 'w0': w0, 'idelta': idelta, 'sw_spsym': sw_spsym})
    return res


def write_nmr_dat(res: dict, fname: str = 'nmr_sc.dat') -> None:
    """
    @fn write_nmr_dat
    @brief Write the NMR temperature sweep to a text file.
    @param   res: dict returned by nmr_sweep
    @param fname: output file name
    """
    with open(fname, 'w') as f:
        f.write(f'# NMR observables in the SC state  (Tc={res["temp_c"]:.6e} eV, '
                f'w0={res["w0"]:.3e} eV, idelta={res["idelta"]:.3e} eV, '
                f'{"triplet d||z" if res["sw_spsym"] else "singlet"})\n')
        f.write('# T[eV]  T/Tc  D(T)/D(0)  K_sc  K_n  K_sc/K_n  '
                '(1/T1T)_sc  (1/T1T)_n  ratio\n')
        for i in range(len(res['temp'])):
            f.write(f'{res["temp"][i]:12.6e} {res["t_red"][i]:8.5f} {res["gap_ratio"][i]:8.5f} '
                    f'{res["K_sc"][i]:14.6e} {res["K_n"][i]:14.6e} {res["K_ratio"][i]:12.6f} '
                    f'{res["R_sc"][i]:14.6e} {res["R_n"][i]:14.6e} {res["R_ratio"][i]:12.6f}\n')


def plot_nmr(res: dict, fname: str = 'nmr_sc.png') -> None:
    """
    @fn plot_nmr
    @brief Two-panel summary plot of the NMR sweep.
    Left : K_s/K_n vs T/Tc (linear) -- Yosida curve for a singlet, finite residual
           for a triplet d||z.
    Right: (1/T1T)_s/(1/T1T)_n vs T/Tc (log-log) -- a Hebel-Slichter peak just below
           Tc for a full isotropic gap, its suppression for a sign-changing (s+-) gap,
           and the low-T power law (1/T1 ~ T^3 for line nodes, ~T^5 for point nodes,
           exponential for a full gap) read off as the slope.
    @param   res: dict returned by nmr_sweep
    @param fname: output png file name
    """
    import matplotlib.pyplot as plt
    t = res['t_red']
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))

    ax[0].plot(t, res['K_ratio'], 'o-', color='tab:blue', ms=4)
    ax[0].axhline(1.0, ls=':', lw=0.8, color='gray')
    ax[0].axvline(1.0, ls='--', lw=0.8, color='gray')
    ax[0].set_xlabel(r'$T/T_c$')
    ax[0].set_ylabel(r'$K_s/K_n$')
    ax[0].set_title('Knight shift')
    ax[0].set_xlim(0, max(1.05, t.max()))
    ax[0].set_ylim(bottom=0)

    ax[1].loglog(t, res['R_ratio'], 'o-', color='tab:red', ms=4)
    ax[1].axhline(1.0, ls=':', lw=0.8, color='gray')
    ax[1].axvline(1.0, ls='--', lw=0.8, color='gray')
    tt = t[t > 0]
    if len(tt):
        # The power law is quoted for 1/T1, but the axis shows 1/(T1 T): line nodes give
        # 1/T1 ~ T^3, hence a slope of 2 here.  Anchored at the highest temperature shown.
        ref = (tt / tt.max()) ** 2
        ax[1].loglog(tt, ref, 'k--', lw=0.8, label=r'$1/T_1\propto T^{3}$ (line nodes)')
        ax[1].legend(fontsize=8)
    ax[1].set_xlabel(r'$T/T_c$')
    ax[1].set_ylabel(r'$(1/T_1T)_s/(1/T_1T)_n$')
    ax[1].set_title(r'spin-lattice relaxation')

    fig.tight_layout()
    fig.savefig(fname, dpi=200)
    plt.close(fig)
