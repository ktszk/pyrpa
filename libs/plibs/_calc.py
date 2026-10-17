#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
Higher-level calculation routines: FLEX, Eliashberg, susceptibility spectra, carrier analysis.
"""
import numpy as np
import scipy.linalg as sclin
import libs.flibs as flibs
import libs.plibs._bands as _bands
import libs.plibs._eilenberger as _eilenberger
import libs.plibs._response as _response
import libs.plibs._wannier_io as _wannier_io

# Public API of this submodule; keeps `from ._x import *` in __init__.py from
# re-exporting numpy, scipy, ctypes and the sibling packages imported above.
__all__ = [
    'regrid_sigma_bin', 'calc_path_spectrum', 'hall_coefficient', 'hall_carrier_density', 'fs_carrier_density',
    'report_band_filling', 'get_mu', 'gap_extrapolate_w0', 'output_gap_function', 'calc_flex', 'calc_flex_soc',
    'calc_lin_eliashberg_eq', 'calc_lin_eliash_soc', 'calc_eliashberg_eq'
]


def _prepare_kspace_state(Nx: int, Ny: int, Nz: int, ham_r: np.ndarray, S_r: np.ndarray, rvec: np.ndarray,
                          need_ham: bool = False) -> dict:
    klist, kmap, invk = flibs.gen_irr_k_TRS(Nx, Ny, Nz)
    eig, uni = _bands.get_eigs(klist, ham_r, S_r, rvec)
    state = {
        'klist': klist,
        'kmap': kmap,
        'invk': invk,
        'eig': eig,
        'uni': uni,
    }
    if need_ham:
        state['ham_k'] = flibs.gen_ham(klist, ham_r, rvec)
    return state


def _prepare_normal_interaction(olist: np.ndarray, site: np.ndarray, orb_dep: bool, U: float, J: float,
                                Umat: np.ndarray | None = None,
                                Jmat: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    if orb_dep:
        return flibs.gen_SCmatrix_orb(olist, site, Umat, Jmat)
    return flibs.gen_SCmatrix(olist, site, U, J)


def _prepare_soc_interaction(olist: np.ndarray, slist: np.ndarray, site: np.ndarray, invs: np.ndarray, orb_dep: bool,
                             U: float, J: float, Umat: np.ndarray | None = None,
                             Jmat: np.ndarray | None = None) -> np.ndarray:
    if orb_dep:
        return flibs.gen_Vmatrix_orb(olist, slist, site, invs, Umat, Jmat)
    return flibs.gen_Vmatrix(olist, slist, site, invs, U, J)


def _load_sigma_from_file() -> tuple | None:
    """Load (sigmak, mu_self) from self_en.npz, or None if the file is missing."""
    try:
        npz = np.load('self_en.npz')
        print("Import sigma from self_en.npz", flush=True)
        return npz['arr_0'], npz['arr_1']
    except FileNotFoundError:
        print("Error: 'self_en.npz' not found", flush=True)
        return None


def _read_sigma_bin_meta(path: str = 'sigma.bin') -> tuple | None:
    """Metadata footer (temp, Nw, Norb, Nk) of a sigma.bin as a tuple, or None
    if the file is missing, malformed, or carries no footer (older io_sigma).
    Walks the sequential records by their 4-byte markers without loading the
    sigma payload; records larger than 2 GiB (compiler subrecords) are not
    handled and yield None."""
    import os, struct
    try:
        fsize = os.path.getsize(path)
        with open(path, 'rb') as f:
            last = None
            while f.tell() < fsize:
                head = f.read(4)
                if len(head) < 4:
                    return None
                n, = struct.unpack('=i', head)
                if n < 0 or f.tell() + n + 4 > fsize:
                    return None
                if n == 32:
                    last = f.read(32)
                else:
                    f.seek(n, 1)
                    last = None
                tail = f.read(4)
                if len(tail) < 4 or struct.unpack('=i', tail)[0] != n:
                    return None
    except OSError:
        return None
    if last is None:
        return None
    return struct.unpack('=4d', last)


def _auto_regrid_sigma_seed(temp: float, Nw: int, path: str = 'sigma.bin') -> None:
    """Bring the sigma.bin seed loaded via sw_in_self onto the current run's
    Matsubara mesh: if the metadata footer reports a different Nw, the file is
    re-gridded in place with regrid_sigma_bin before the Fortran side reads it.
    w_scale=temp_old/temp reproduces the frequency compression of the raw
    equal-Nw reuse, so T-annealing behaves the same whether or not Nw changed.
    Footer-less files (older io_sigma) are left untouched, and a Norb/Nk
    mismatch is left to the Fortran-side mesh check (such a seed is unusable
    for the run either way)."""
    meta = _read_sigma_bin_meta(path)
    if meta is None:
        return
    temp_old, Nw_old = meta[0], int(round(meta[1]))
    if Nw_old == Nw:
        return
    print(f"sw_in_self: sigma.bin seed has Nw={Nw_old}, run uses Nw={Nw}: "
          f"re-gridding in place", flush=True)
    regrid_sigma_bin(temp_new=temp, Nw_new=Nw, w_scale=temp_old / temp, path=path)


def regrid_sigma_bin(temp_old: float = None, temp_new: float = None, Norb: int = None,
                     Nw_old: int = None, Nw_new: int = None, w_scale: float = None,
                     path: str = 'sigma.bin', out_path: str = None) -> tuple:
    """
    @brief Re-grid the FLEX self-energy seed sigma.bin from the Matsubara mesh at
    temp_old onto the mesh at temp_new (T-annealing helper for sw_in_self).

    Sigma(k, i w_n) is interpolated along the imaginary axis by a cubic spline
    after extending the data to negative frequencies with the conjugate
    symmetry Sigma_lm(-iw) = Sigma_ml(iw)^* (TRS), so the low-frequency end
    needs no extrapolation.  The outermost stored frequency is discarded
    beforehand: it carries the wrap-around artifact of the circular Sigma
    convolution (its Im part can even flip sign) and would poison both the
    spline end and the tail fit.  New points beyond the remaining cutoff are
    filled with the two-parameter tail Sigma ~ c0 + c1/(iw) fitted to the two
    outermost kept points.  Interpolated diagonal elements are clipped to
    Im Sigma_ll <= 0 (causality).

    w_scale: the stored Sigma is evaluated at w_scale * w_n(temp_new) instead
    of w_n(temp_new).  The default 1.0 gives the faithful Sigma(T_old) on the
    new mesh.  w_scale = temp_old/temp_new reproduces the frequency
    compression of the raw index-mapped reuse (same file, no regrid) on an
    ARBITRARY Nw_new: empirically this seeds near-critical cooling runs much
    better, because the compression inflates the low-frequency Sigma by
    T_old/T_new, mimicking the physical growth of Sigma upon cooling.

    mu and mu_OLD stored in the file are passed through unchanged.  The file
    layout matches io_sigma in fself.f90 (sequential unformatted, 4-byte
    record markers; records: mu, mu_OLD, sigmak(Nk,Nw,Norb,Norb) in Fortran
    order, plus a trailing metadata record (temp, Nw, Norb, Nk) written by
    recent io_sigma versions).  When the metadata footer is present temp_old,
    Norb and Nw_old are read from it and may be omitted; values passed
    explicitly are checked against the footer and a mismatch raises
    ValueError.  Files written by older versions (no footer) require all
    three explicitly.  The re-gridded file always carries a footer with the
    NEW mesh, so chained annealing steps need no manual bookkeeping.
    The SOC path (mkself_soc) writes the same layout, so its files work here
    too; only files from the retired element-wise SOC format are unreadable.
    Records larger than 2 GiB (compiler subrecords) are not handled.

    @param temp_old: temperature [eV] the file was written at
                     (default: from the metadata footer)
    @param temp_new: target temperature [eV] (required)
    @param     Norb: number of orbitals (default: from the metadata footer)
    @param   Nw_old: Matsubara count of the stored sigma
                     (default: from the metadata footer)
    @param   Nw_new: Matsubara count of the target mesh (default: Nw_old)
    @param  w_scale: evaluation-frequency factor (default 1.0 = faithful;
                     temp_old/temp_new = index-map-like compression, see above)
    @param     path: input file (default 'sigma.bin')
    @param out_path: output file (default: overwrite path)
    @return (Nk, Nw_new): shape information of the written sigma
    """
    import scipy.io as scio
    import scipy.interpolate as scinterp
    if temp_new is None:
        raise ValueError("temp_new is required")
    if w_scale is None:
        w_scale = 1.0
    with scio.FortranFile(path, 'r') as f:
        mu = f.read_record(np.float64)
        mu_old = f.read_record(np.float64)
        flat = f.read_record(np.complex128)
        try:
            meta = f.read_record(np.float64)
            if meta.size != 4:
                meta = None
        except Exception:
            meta = None
    if meta is not None:
        stored = {'temp_old': meta[0], 'Nw_old': int(round(meta[1])),
                  'Norb': int(round(meta[2]))}
        for name, given in (('temp_old', temp_old), ('Nw_old', Nw_old),
                            ('Norb', Norb)):
            if given is not None and not np.isclose(given, stored[name],
                                                    rtol=1e-8, atol=0.0):
                raise ValueError(f"{name}={given} contradicts the sigma.bin "
                                 f"metadata footer ({stored[name]})")
        temp_old, Nw_old, Norb = stored['temp_old'], stored['Nw_old'], stored['Norb']
    elif temp_old is None or Nw_old is None or Norb is None:
        raise ValueError("sigma.bin carries no metadata footer (written by an "
                         "older io_sigma): temp_old, Norb and Nw_old must be "
                         "given explicitly")
    if Nw_new is None:
        Nw_new = Nw_old
    Nk, rem = divmod(flat.size, Nw_old * Norb * Norb)
    if rem != 0:
        raise ValueError(f"sigma.bin size {flat.size} is not divisible by "
                         f"Nw_old*Norb^2 = {Nw_old * Norb * Norb}: wrong Norb/Nw_old?")
    if meta is not None and Nk != int(round(meta[3])):
        raise ValueError(f"inferred Nk={Nk} contradicts the metadata footer "
                         f"({int(round(meta[3]))})")
    sig = flat.reshape((Nk, Nw_old, Norb, Norb), order='F')
    # drop the outermost frequency (circular-convolution wrap-around artifact)
    sig = sig[:, :Nw_old - 1]
    w_old = (2.0 * np.arange(Nw_old - 1) + 1.0) * np.pi * temp_old
    w_new = w_scale * (2.0 * np.arange(Nw_new) + 1.0) * np.pi * temp_new
    # conjugate-symmetric extension to negative frequencies
    mirror = np.conj(np.swapaxes(sig, 2, 3))[:, ::-1, :, :]
    wext = np.concatenate([-w_old[::-1], w_old])
    ext = np.concatenate([mirror, sig], axis=1)
    out = scinterp.CubicSpline(wext, ext, axis=1, extrapolate=False)(w_new)
    beyond = w_new > w_old[-1]
    if beyond.any():
        # tail Sigma ~ c0 + c1/(iw), least-squares over an interior window:
        # the last ~10-20% of the stored frequencies carry the growing
        # wrap-around contamination of the circular convolution, so fitting
        # at the very edge produces wildly wrong coefficients
        nk = len(w_old)
        i0, i1 = int(0.55 * nk), max(int(0.85 * nk), int(0.55 * nk) + 4)
        basis = np.stack([np.ones(i1 - i0), 1.0 / (1j * w_old[i0:i1])], axis=1)
        rhs = sig[:, i0:i1].reshape(Nk, i1 - i0, Norb * Norb).transpose(1, 0, 2)
        coef, *_ = np.linalg.lstsq(basis, rhs.reshape(i1 - i0, -1), rcond=None)
        c0 = coef[0].reshape(Nk, Norb, Norb)
        c1 = coef[1].reshape(Nk, Norb, Norb)
        out[:, beyond, :, :] = (c0[:, None] +
                                c1[:, None] / (1j * w_new[beyond])[None, :, None, None])
    di = np.arange(Norb)
    imdiag = out[:, :, di, di].imag
    nbad = int((imdiag > 1e-12).sum())
    if nbad:
        print(f"regrid_sigma_bin: clipped {nbad} causality-violating "
              f"diagonal points (Im Sigma_ll > 0)", flush=True)
    out[:, :, di, di] = out[:, :, di, di].real + 1j * np.minimum(imdiag, 0.0)
    with scio.FortranFile(out_path or path, 'w') as f:
        f.write_record(mu)
        f.write_record(mu_old)
        f.write_record(np.ascontiguousarray(out.ravel(order='F')))
        f.write_record(np.array([temp_new, float(Nw_new), float(Norb), float(Nk)]))
    print(f"regrid_sigma_bin: sigma.bin re-gridded T={temp_old:.6f} -> {temp_new:.6f} eV, "
          f"Nw={Nw_old} -> {Nw_new} (Nk={Nk}, Norb={Norb})", flush=True)
    return Nk, Nw_new


def _prepare_green_state_normal(state: dict, olist: np.ndarray, interaction: tuple, mu: float, fill: float, temp: float,
                                Nw: int, Nx: int, Ny: int, Nz: int, sw_self: bool,
                                sw_from_file: bool = False, sw_out_self: bool = False, sw_in_self: bool = False,
                                eps: float = 1.0e-4, pp: float = 0.5, m_diis: int = 5, sw_rescale: bool = True,
                                sw_tail: bool = False, sigma_scale: float = 1.0) -> dict | None:
    Smat, Cmat = interaction
    if sw_self:
        if sw_from_file:
            loaded = _load_sigma_from_file()
            if loaded is None:
                return None
            sigmak, mu_self = loaded
        else:
            if sw_in_self:
                _auto_regrid_sigma_seed(temp, Nw)
            sigmak, mu_self = flibs.mkself(Smat, Cmat, state['kmap'], state['invk'], olist,
                                           state['ham_k'], state['eig'], state['uni'],
                                           mu, fill, temp, Nw, Nx, Ny, Nz, sw_out_self, sw_in_self,
                                           eps=eps, pp=pp, m_diis=m_diis, sw_rescale=sw_rescale,
                                           sw_tail=sw_tail, sigma_scale=sigma_scale)
        print(f'chem. pot. with self= {mu_self:.4f} eV', flush=True)
        Gk = flibs.gen_green(sigmak, state['ham_k'], mu_self, temp)
    else:
        sigmak = None
        mu_self = mu
        Gk = flibs.gen_Green0(state['eig'], state['uni'], mu, temp, Nw)
    return {'Gk': Gk, 'mu_self': mu_self, 'sigmak': sigmak}


def _prepare_green_state_soc(state: dict, olist: np.ndarray, slist: np.ndarray, invs: np.ndarray, interaction: tuple,
                             mu: float, fill: float, temp: float, Nw: int, Nx: int, Ny: int, Nz: int, sw_self: bool,
                             sw_from_file: bool = False, sw_out_self: bool = False, sw_in_self: bool = False,
                             eps: float = 1.0e-4, pp: float = 0.5, m_diis: int = 5, sw_rescale: bool = True,
                             sigma_scale: float = 1.0) -> dict | None:
    Vmat = interaction
    if sw_self:
        if sw_from_file:
            loaded = _load_sigma_from_file()
            if loaded is None:
                return None
            sigmak, mu_self = loaded
        else:
            if sw_in_self:
                _auto_regrid_sigma_seed(temp, Nw)
            sigmak, mu_self = flibs.mkself_soc(Vmat, state['kmap'], state['invk'], invs, olist, slist,
                                               state['ham_k'], state['eig'], state['uni'],
                                               mu, fill, temp, Nw, Nx, Ny, Nz, sw_out_self, sw_in_self,
                                               eps=eps, pp=pp, m_diis=m_diis, sw_rescale=sw_rescale,
                                               sigma_scale=sigma_scale)
        print(f'chem. pot. with self= {mu_self:.4f} eV', flush=True)
        Gk = flibs.gen_green(sigmak, state['ham_k'], mu_self, temp)
    else:
        sigmak = None
        mu_self = mu
        Gk = flibs.gen_Green0(state['eig'], state['uni'], mu, temp, Nw)
    return {'Gk': Gk, 'mu_self': mu_self, 'sigmak': sigmak}


def _get_chi0_normal(state: dict, Gk: np.ndarray, chiolist: np.ndarray, Smat: np.ndarray, Cmat: np.ndarray,
                     mu_self: float, temp: float, Nx: int, Ny: int, Nz: int, sw_tail: bool) -> tuple[np.ndarray, float]:
    """Irreducible chi0 on the irreducible q-grid; sw_tail=True uses the
    tail-corrected convolution (conv[G]-conv[G0]+analytic reference, needs the
    chemical potential mu_self of Gk)."""
    if sw_tail:
        return flibs.get_chi0_tail(Smat, Cmat, Gk, state['eig'], state['uni'], chiolist,
                                   state['kmap'], state['invk'], mu_self, temp, Nx, Ny, Nz)
    return flibs.get_chi0(Smat, Cmat, Gk, chiolist, state['kmap'], state['invk'], temp, Nx, Ny, Nz)


def calc_path_spectrum(kind: str, mu: float, temp: float, klist: np.ndarray, qlist: np.ndarray, chiolist: np.ndarray,
                       eig: np.ndarray, uni: np.ndarray, spa_length: np.ndarray, Nw: int, Emax: float, delta: float,
                       Smat: np.ndarray | None = None) -> tuple:
    """
    @brief Compute a susceptibility spectrum along a symmetry line and write it out.

    Thin driver over _response.chis_spectrum / phi_spectrum: it runs the solver,
    writes the trace to '@<kind@>.dat' as (path length, omega, Im chi) blocks
    separated by a blank line (gnuplot pm3d layout), and for 'chis' additionally
    writes the per-orbital decomposition to 'chis_@<i@>.dat'.  A failed write is
    reported and skipped, not raised, so a long run is never lost to a full disk.

    @param    kind: 'chis' (spin susceptibility) or 'phi' (pairing susceptibility)
    @param      mu: chemical potential [eV]
    @param    temp: temperature k_B*T [eV]
    @param   klist: k-points of the BZ sum in fractional coordinates [Nk, 3]
    @param   qlist: q-points along the symmetry line [Nq, 3]
    @param chiolist: orbital-pair basis of chi [Nchi, 2] int64
    @param     eig: eigenvalues on klist [Nk, Norb]
    @param     uni: eigenvectors on klist [Nk, Norb, Norb]
    @param spa_length: cumulative path length of each q-point [Nq], the x axis of the plot
    @param      Nw: number of real frequencies
    @param    Emax: upper end of the frequency window [eV]
    @param   delta: spectral broadening eta [eV]
    @param    Smat: static spin vertex [Nchi, Nchi]; required for kind='chis', unused for 'phi'
    @return (w, sp, spec): frequency mesh, path-length mesh and the spectrum, all [Nq, Nw]
    """
    if kind == 'chis':
        print("calculate spin susceptibility", flush=True)
        spec, spec_orb, wlist = _response.chis_spectrum(mu, temp, Smat, klist, qlist, chiolist, eig, uni, Nw, Emax,
                                                        delta)
        outname = 'chis.dat'
    elif kind == 'phi':
        print("calculate sc susceptibility", flush=True)
        spec, spec_orb, wlist = _response.phi_spectrum(mu, temp, klist, qlist, chiolist, eig, uni, Nw, Emax, delta)
        outname = 'phi.dat'
    else:
        raise ValueError(f"Unknown spectrum kind: {kind}")

    w, sp = np.meshgrid(wlist, spa_length)
    try:
        with open(outname, 'w') as f:
            for ww, ssp, spec_row in zip(w, sp, spec):
                for www, sssp, spec_val in zip(ww, ssp, spec_row):
                    f.write(f'{sssp:8.4f} {www:8.4f} {spec_val.imag:9.4f}\n')
                f.write('\n')
    except OSError as e:
        print(f"Error: Failed to write '{outname}': {e}", flush=True)

    if kind == 'chis':
        for i, spec_orb_row in enumerate(spec_orb.T):
            try:
                with open(f'chis_{i}.dat', 'w') as f:
                    for ww, ssp, spec_band in zip(w, sp, spec_orb_row.T):
                        for www, sssp, spec_val in zip(ww, ssp, spec_band):
                            f.write(f'{sssp:8.4f} {www:8.4f} {spec_val.imag:9.4f}\n')
                        f.write('\n')
            except OSError as e:
                print(f"Error: Failed to write 'chis_{i}.dat': {e}", flush=True)
                continue

    return w, sp, spec

# --------------------------------------------------------------------------- #
#  Weak-field (B-linear) Hall response
# --------------------------------------------------------------------------- #
#  The k-space kernels themselves live in fortran (fcond.f90):
#    flibs.calc_sigma_hall  -> S[Norb,3,3,3], the B-linear Hall kernel for all
#                              three field directions plus its absolute-value sum
#    flibs.calc_k0_band     -> K0[Norb,3,3], the per-band charge kernel
#  Only the 3x3 algebra that turns them into R_H stays here.


def hall_coefficient(K0: np.ndarray, Shall: np.ndarray, Vuc: float, wsum: float,
                     gsp: float) -> np.ndarray:
    """
    @brief Hall coefficient for the three field directions from the full
    conductivity and Hall-kernel tensors.

        rho^(1) = -sigma^-1 sigma^(1) sigma^-1
        R_H^(g) = rho^(1)_{ba}/B = -(wsum*Vuc/gsp) [K0^-1 S^(g) K0^-1]_{b,a}

    with (a,b,g) cyclic.  All unit prefactors cancel between sigma^(1) and the
    two factors of sigma^-1 except (wsum*Vuc/gsp), so tau (and tau_unit) drop
    out of R_H for a k-independent tau -- which is why this is expressed in the
    bare calc_Kn / calc_sigma_hall sums.

    @param     K0: Boltzmann charge kernel [3,3] from calc_Kn (proportional to sigma)
    @param  Shall: Hall kernel [3,3,3] from flibs.calc_sigma_hall (band summed,
                  or one band's slice for a per-sheet R_H)
    @param    Vuc: Unit cell volume in m^3
    @param   wsum: Sum of the k-point weights (Nk for a uniform mesh)
    @param    gsp: Spin degeneracy factor (2 without SOC, 1 with)
    @return    Rh: Hall coefficient [3] in m^3/C, index = field direction
    """
    pref = -(wsum*Vuc/gsp)
    scale = (max(np.trace(K0), 0.0)/3.)**3
    Rh = np.full(3, np.nan)
    full_ok = abs(sclin.det(K0)) > 1.0e-12*scale
    K0inv = sclin.inv(K0) if full_ok else None
    for g in range(3):
        a, b = (g+1) % 3, (g+2) % 3
        if full_ok:
            Rh[g] = pref*K0inv.dot(Shall[:, :, g]).dot(K0inv)[b, a]
            continue
        # sigma singular: a band with no dispersion along some axis (a strictly 2D
        # model has sigma_zz = 0) makes the full tensor non-invertible while R_H
        # for B along that axis is still perfectly well defined.  Fall back to the
        # 2x2 block in the plane the current actually flows in, which reproduces
        # the diagonal expression exactly when sigma has no coupling to g.
        idx = np.ix_([a, b], [a, b])
        K2 = K0[idx]
        if abs(sclin.det(K2)) > 1.0e-12*(max(np.trace(K2), 0.0)/2.)**2:
            K2inv = sclin.inv(K2)
            Rh[g] = pref*K2inv.dot(Shall[:, :, g][idx]).dot(K2inv)[1, 0]
    return Rh


def hall_carrier_density(Rh: np.ndarray, eC: float) -> tuple[np.ndarray, list]:
    """
    @brief Convert R_H into a POSITIVE carrier density plus the carrier sign.

    n_H = 1/(e|R_H|).  The single-band relations are R_H = -1/(n_e e) for
    electrons and R_H = +1/(n_h e) for holes, so the sign of R_H selects the
    carrier type and only |R_H| sets the magnitude.  (The old
    nh = -1/(R_H e) returned a NEGATIVE number for a hole-like band and
    labelled it a hole density.)

    Valid as a true carrier count only for a single closed sheet; for a
    multi-band / compensated metal n_H is the two-fluid combination
    (p mu_h^2 - n mu_e^2)/... , not n or p.  Compare against the Luttinger
    volume from fs_carrier_density before reading it as a density.

    @param  Rh: Hall coefficient [3] in m^3/C
    @param  eC: Elementary charge in C
    @return (n_H, kind): density [3] in cm^-3 and 'electron'/'hole' per direction
    """
    Rh = np.atleast_1d(np.asarray(Rh, dtype=np.float64))
    with np.errstate(divide='ignore'):
        n_H = 1.0/(np.abs(Rh)*eC)/1.0e6                 # m^-3 -> cm^-3
    kind = ['-' if not np.isfinite(r) else ('electron' if r < 0 else 'hole') for r in Rh]
    return n_H, kind


def _periodic_components(mask: np.ndarray) -> list:
    """
    @brief Connected components of a 3-D boolean k-space mask under PERIODIC
    boundaries, each tagged with its volume fraction and whether it spans the zone.

    Both properties are essential here.  The periodic merge matters because a
    pocket centred on Gamma is split across the array edges of the [0,1) mesh
    get_emesh uses, and would be counted as up to 8 separate pockets.  The span
    flag separates a CLOSED pocket, whose volume is a Luttinger carrier count,
    from the percolating background of a partly filled band (and from a genuinely
    open Fermi sheet, which has no carrier count at all).

    @param mask: occupied (or unoccupied) region on the [Nx,Ny,Nz] mesh
    @return list of (fraction, spans) per component, largest fraction first
    """
    import scipy.ndimage as ndimage
    lab, n = ndimage.label(mask)
    if n == 0:
        return []
    parent = list(range(n+1))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for ax in range(3):                                  # glue the opposite faces
        lo, hi = np.take(lab, 0, axis=ax), np.take(lab, -1, axis=ax)
        sel = (lo > 0) & (hi > 0)
        if not sel.any():
            continue
        # collapse the face to its DISTINCT label pairs first: a face carries N^2
        # cells but only a handful of pairs, and the union loop is python
        pairs = np.unique(np.stack([lo[sel].ravel(), hi[sel].ravel()], axis=1), axis=0)
        for u, v in pairs:
            ru, rv = find(int(u)), find(int(v))
            if ru != rv:
                parent[rv] = ru
    labels = np.array([find(i) for i in range(n+1)])[lab]
    nlab = n+1
    size = np.bincount(labels.ravel(), minlength=nlab)
    # "spans axis ax" <=> the component appears in EVERY layer along ax.  Counting
    # occupied layers with one bincount per layer is O(Nk) for the whole axis;
    # testing each component with its own full-array mask would be O(Ncomp*Nk)
    # and dominates the carrier analysis on a fine mesh.
    spans = np.zeros(nlab, dtype=bool)
    for ax in range(3):
        n_ax = labels.shape[ax]
        cnt = np.zeros(nlab, dtype=np.int64)
        for i in range(n_ax):
            cnt += np.bincount(np.take(labels, i, axis=ax).ravel(), minlength=nlab) > 0
        spans |= (cnt == n_ax)
    total = labels.size
    out = [(size[u]/total, bool(spans[u])) for u in range(1, nlab) if size[u] > 0]
    return sorted(out, key=lambda p: -p[0])


def fs_carrier_density(eig: np.ndarray, mu: float, Vuc: float, gsp: float = 2.0,
                       mesh: tuple | None = None) -> dict:
    """
    @brief Carrier density from the LUTTINGER VOLUME of every closed Fermi-surface
    pocket -- the quantity an experimental carrier density is actually compared to.

    Each band's occupied region (e@<mu) and unoccupied region (e@>=mu) are split
    into connected components under periodic boundary conditions.  A CLOSED
    occupied component is an electron pocket holding gsp*frac/Vuc carriers; a
    CLOSED unoccupied component is a hole pocket holding the same of the opposite
    sign.  The component that percolates through the zone is the filled/empty
    background of that band and carries nothing.

    Working per pocket rather than per band matters: bands are sorted by energy at
    each k, so one band index routinely carries an electron pocket in one part of
    the zone and a hole pocket in another (wherever two bands cross).  A single
    occupied fraction per band would assign that band one type and silently drop
    the other pocket.

    An OPEN sheet (both regions percolate) has no Luttinger volume and is reported
    as such rather than as a number.

    This is NOT what the Hall coefficient returns.  n_H = 1/(e|R_H|) coincides
    with it only for a single closed sheet; comparing the two is how a
    discrepancy with experiment gets attributed to compensation rather than to a
    misplaced mu.

    @param   eig: Eigenvalues [Nk,Norb] float64 (eV), on the uniform [0,1) mesh
    @param    mu: Chemical potential in eV
    @param   Vuc: Unit cell volume in m^3
    @param   gsp: Spin degeneracy factor (2 without SOC, 1 with)
    @param  mesh: (Nx,Ny,Nz).  Required for the pocket decomposition; if None the
                  routine falls back to one occupied-fraction per band and sets
                  'approx' True.
    @return dict with per-band arrays 'ne_band','nh_band' [cm^-3], 'npocket_e',
            'npocket_h', 'open_sheet', 'focc', 'kind', and totals 'n_e','n_h'
    """
    Norb = eig.shape[1]
    scale = gsp/(Vuc*1.0e6)                              # BZ fraction -> cm^-3
    focc = np.array([(eig[:, n] < mu).mean() for n in range(Norb)])
    ne_band = np.zeros(Norb)
    nh_band = np.zeros(Norb)
    npocket_e = np.zeros(Norb, dtype=int)
    npocket_h = np.zeros(Norb, dtype=int)
    open_sheet = np.zeros(Norb, dtype=bool)
    kind = []
    for n in range(Norb):
        if not eig[:, n].min() < mu < eig[:, n].max():
            kind.append('full' if focc[n] > 0.5 else 'empty')
            continue
        if mesh is None:                                 # crude single-pocket fallback
            if focc[n] <= 0.5:
                ne_band[n], npocket_e[n] = scale*focc[n], 1
            else:
                nh_band[n], npocket_h[n] = scale*(1.-focc[n]), 1
        else:
            band = eig[:, n].reshape(mesh)
            for frac, spans in _periodic_components(band < mu):
                if not spans:
                    ne_band[n] += scale*frac
                    npocket_e[n] += 1
            for frac, spans in _periodic_components(band >= mu):
                if not spans:
                    nh_band[n] += scale*frac
                    npocket_h[n] += 1
            open_sheet[n] = (npocket_e[n] == 0 and npocket_h[n] == 0)
        if open_sheet[n]:
            kind.append('open')
        elif npocket_e[n] and npocket_h[n]:
            kind.append('both')
        elif npocket_h[n]:
            kind.append('hole')
        else:
            kind.append('electron')
    return {'ne_band': ne_band, 'nh_band': nh_band, 'npocket_e': npocket_e,
            'npocket_h': npocket_h, 'open_sheet': open_sheet, 'focc': focc,
            'kind': np.array(kind), 'approx': mesh is None,
            'n_e': ne_band.sum(), 'n_h': nh_band.sum()}


def report_band_filling(kmesh: int, rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray, mu: float,
                        Arot: np.ndarray) -> None:
    """
    @brief Print, for each band, the fraction of the mesh lying above and below mu.

    The unoccupied fraction is the hole-like side of that band and the occupied
    fraction the electron-like side, so the pair says at a glance whether a sheet
    sits near the bottom or the top of its band.  The occupied fractions sum to
    the electron count per unit cell, which is the consistency check against `fill`.

    Both fractions ARE Luttinger volumes, normalized per unit cell per spin
    rather than per cm^3: a band with no Fermi surface reads 1 (one electron per
    cell) or 0, and a band with a pocket reads that pocket's V_occ/V_BZ on the
    corresponding side.  fs_carrier_density reports the same volumes multiplied
    by gsp/V_uc, so n*V_uc/gsp returns the fraction printed here exactly.  This
    normalization is the one that pairs with `fill`; the cm^-3 one is what an
    experimental carrier concentration is compared against.
    @param  kmesh: k-mesh (used for all three axes)
    @param   rvec: Real-space lattice vectors (Wannier R-vectors)
    @param  ham_r: Hamiltonian in real space
    @param    S_r: Real-space overlap blocks ([] for an orthogonal basis)
    @param     mu: Chemical potential in eV
    @param   Arot: Rotation matrix passed through to get_emesh
    """
    Nk, eig, kweight = _bands.get_emesh(kmesh, kmesh, kmesh, ham_r, S_r, rvec, Arot)
    if Nk <= 0:
        print("Error: Number of k-points (Nk) is non-positive", flush=True)
        return
    fill = 0.0
    for i, en in enumerate(eig.T - mu):
        num_hole = float(np.where(en > 0)[0].size) / Nk
        num_particle = float(np.where(en <= 0)[0].size) / Nk
        print(i+1, round(num_hole, 4), round(num_particle, 4), flush=True)
        fill += num_particle
    print(f'sum of electrons is {fill:.4f}', flush=True)


def get_mu(ham_r: np.ndarray, S_r: np.ndarray, rvec: np.ndarray, Arot: np.ndarray, temp: float, fill: float,
           kmesh: int = 40) -> float:
    """
    @brief Chemical potential of the normal state on a cubic kmesh^3 mesh.

    Convenience wrapper around get_emesh + calc_mu for the callers that only need
    mu and have no use for the eigenvalues.  Returns None (after printing) rather
    than raising when temp or kmesh is non-positive, matching the other drivers
    in this module.

    @param   ham_r: hopping matrix H(R) [Nr, Norb, Norb] complex128
    @param     S_r: overlap matrix S(R) for a non-orthogonal (MLO) basis; [] if orthogonal
    @param    rvec: R-vectors [Nr, 3]
    @param    Arot: rotation matrix passed through to get_emesh (None = no rotation)
    @param    temp: temperature k_B*T [eV]
    @param    fill: band filling; mu solves sum f(eps-mu) = Nk*fill
    @param   kmesh: linear mesh size; the sum runs over kmesh^3 points
    @return mu: chemical potential [eV], or None when temp <= 0 or kmesh <= 0
    """
    if temp <= 0:
        print("Error: Temperature (temp) is non-positive", flush=True)
        return None
    if kmesh <= 0:
        print("Error: k-mesh size (kmesh) is non-positive", flush=True)
        return None
    print("calc chem. pot.", flush=True)
    print(f"band filling = {fill:.4f}", flush=True)
    Nk, eig, kweight = _bands.get_emesh(kmesh, kmesh, kmesh, ham_r, S_r, rvec, Arot)
    mu = _bands.calc_mu(eig, Nk, fill, temp)
    return mu


def gap_extrapolate_w0(gap: np.ndarray, temp: float, n_points: int = 4,
                       order: int = 1) -> tuple[np.ndarray, np.ndarray]:
    """
    @brief Extrapolate the gap function Delta(k,iw_n) to the w_n->0 limit along the
    imaginary axis, by a low-order polynomial fit in w_n^2 using the n_points lowest
    fermionic Matsubara points. w_n=(2n+1)*pi*temp is never exactly zero, so this
    limit is reached only by extrapolation, never by reading off a stored data point
    (in particular Delta(iw_0), the lowest point, differs from it by O((pi*temp)^2)).
    Delta(-iw_n)=Delta(iw_n)^* makes Delta(iw_n) even in w_n, so fitting in w_n^2
    (not w_n) removes the leading-order bias rather than just averaging noise.
    If Delta(z) is analytic in a neighborhood of z=0 (the generic case for a fully
    gapped or smoothly-nodal SC state), this imaginary-axis limit coincides with the
    real-axis retarded Delta(w=0).
    @param    gap: gap function [Norb, Norb, Nw, Nk] complex128 (Matsubara axis 3rd,
                  as produced by linearized_eliashberg/nonlinear_eliashberg)
    @param   temp: temperature in eV; sets w_n = (2n+1)*pi*temp
    @param n_points: number of lowest Matsubara points used in the fit (>= order+1)
    @param  order: polynomial order in w_n^2 (1: Delta0 + c*w_n^2, the leading-order
                  correction implied by the even-in-w_n symmetry above)
    @return  gap0: extrapolated Delta(k,w_n->0) [Norb, Norb, Nk] complex128
    @return  bias: |gap0 - Delta(iw_0)| / max|gap0|, the O((pi*temp)^2) correction
                  that reading Delta(iw_0) directly misses, per component but
                  normalized by the PEAK gap [Norb, Norb, Nk] float64. Normalizing
                  by the global peak (not the per-point value) keeps near-nodal /
                  off-diagonal components where gap0~0 from blowing up the reported
                  maximum; bias.max() is then "largest correction as a fraction of
                  the peak gap". It still grows near Tc, where the whole gap shrinks
                  toward the scale of pi*temp, flagging where iw_0 is unreliable.
    """
    if n_points < order + 1:
        raise ValueError(f"n_points ({n_points}) must be >= order+1 ({order+1})")
    Norb1, Norb2, Nw, Nk = gap.shape
    if Nw < n_points:
        raise ValueError(f"gap has only {Nw} Matsubara points, need >= {n_points}")
    wn2 = ((2*np.arange(n_points)+1)*np.pi*temp)**2                    # [n_points]
    data = np.transpose(gap[:, :, :n_points, :], (2, 0, 1, 3)).reshape(n_points, -1)  # [n_points, Norb^2*Nk]
    A = np.vander(wn2, order+1, increasing=True)                       # [n_points, order+1]
    coeffs, *_ = np.linalg.lstsq(A, data, rcond=None)                  # [order+1, Norb^2*Nk]
    gap0 = coeffs[0].reshape(Norb1, Norb2, Nk)
    iw0 = gap[:, :, 0, :]
    gap_scale = max(np.abs(gap0).max(), np.abs(iw0).max())
    if gap_scale == 0.0:
        bias = np.zeros((Norb1, Norb2, Nk), dtype=np.float64)
    else:
        bias = np.abs(gap0 - iw0) / gap_scale
    return gap0, bias


def output_gap_function(invk: np.ndarray, kmap: np.ndarray, gap: np.ndarray, uni: np.ndarray, plist: np.ndarray,
                        gap_sym: int, Nx: int, soc: bool = False, invs: np.ndarray | None = None,
                        slist: np.ndarray | None = None, sw_orb: bool = False, sw_extrapolate: bool = False,
                        temp: float | None = None, n_points: int = 4, order: int = 1,
                        iperm: np.ndarray | None = None) -> int:
    """
    @brief Convert the converged gap to the requested basis and write the kz=0 plane.

    Writes one file per orbital pair, 'gap_@<i@>@<j@>.dat', holding
    (kx index, ky index, Re Delta, Im Delta) for the kz = 0 slice, with a blank
    line after each kx row so gnuplot's pm3d reads it as a grid.  A failed write
    is reported and skipped rather than raised.

    Two conversions may happen first.  With sw_extrapolate the Matsubara axis is
    extrapolated to omega_n -> 0 (gap_extrapolate_w0), which removes the
    O((pi T)^2) bias of using the lowest slice alone; the maximum relative
    correction is printed.  With sw_orb=False the orbital-basis gap is rotated
    into the band basis, which needs the inversion operator (plist / iperm) to
    rebuild Delta(-k).

    @param    invk: index map k -> -k [Nkall, 3] int64
    @param    kmap: fractional index of every k-point [Nkall, 3] int64
    @param     gap: gap function [Norb, Norb, Nw, Nkall] complex128
    @param     uni: eigenvectors used for the band rotation [Nk, Norb, Norb]
    @param   plist: parity sign (+1/-1) per orbital [Norb] float64
    @param gap_sym: pairing symmetry index; its sign selects singlet/triplet
    @param      Nx: k-mesh size along x, used to break the output into rows
    @param     soc: True when the model carries spin-orbit coupling
    @param    invs: spin-flip index map [Norb] int64 (SOC only)
    @param   slist: spin label of each orbital [Norb] int64 (SOC only)
    @param  sw_orb: True keeps the orbital basis, False rotates to the band basis
    @param sw_extrapolate: extrapolate Delta(iw_n) -> Delta(0) before writing
    @param    temp: temperature k_B*T [eV], required when sw_extrapolate is True
    @param n_points: number of Matsubara slices used by the extrapolation fit
    @param   order: polynomial order of that fit
    @param   iperm: orbital permutation of the inversion operator (None = identity)
    @return 0
    """
    if sw_extrapolate:
        gap0, bias = gap_extrapolate_w0(gap, temp, n_points, order)
        print(f'gap w_n->0 extrapolation: max relative correction over iw_0 = {bias.max():.4e}', flush=True)
        gap = gap0[:, :, None, :]  # restore a length-1 Matsubara axis; downstream code reads index 0
    if sw_orb:
        if soc:
            gapb = gap[:, :, 0, :]
        else:
            gapb = flibs.remap_gap(gap[:, :, :1, :], plist, invk, gap_sym, iperm=iperm)[:, :, 0, :]
    else:
        if soc:
            gapb = flibs.conv_delta_orb_to_band_soc(gap, uni, invk, invs, slist, gap_sym)
        else:
            gapb = flibs.conv_delta_orb_to_band(gap, uni, invk, plist, gap_sym, iperm=iperm)
    print('output gap function', flush=True)
    for iorb in range(len(gapb)):
        for jorb in range(len(gapb)):
            try:
                with open(f'gap_{iorb+1}{jorb+1}.dat', 'w') as f:
                    for i, km in enumerate(kmap):
                        if km[2] == 0:
                            f.write(f'{km[0]:3} {km[1]:3} '
                                    f'{gapb[iorb, jorb, i].real:15.8e} '
                                    f'{gapb[iorb, jorb, i].imag:15.8e}\n')
                            if km[0] == Nx - 1:
                                f.write('\n')
            except OSError as e:
                print(f"Error: Failed to write 'gap_{iorb+1}{jorb+1}.dat': {e}", flush=True)
                continue
    return 0


def calc_flex(Nx: int, Ny: int, Nz: int, Nw: int, ham_r: np.ndarray, S_r: np.ndarray, rvec: np.ndarray, mu: float,
              temp: float, olist: np.ndarray, site: np.ndarray, orb_dep: bool, U: float, J: float, fill: float,
              sw_out_self: bool, sw_in_self: bool, Umat: np.ndarray | None = None, Jmat: np.ndarray | None = None,
              eps: float = 1.0e-4, pp: float = 0.5, m_diis: int = 5, sw_rescale: bool = True, sw_tail: bool = False,
              sigma_scale: float = 1.0) -> None:
    """
    @brief Solve the FLEX self-energy loop (no SOC) and optionally export it.

    Builds the irreducible-mesh Green's function and the static Coulomb vertex
    (S/C matrices) from the same orbital-pair basis the response routines use,
    then hands both to the Fortran FLEX driver, which iterates Sigma and mu to
    self-consistency.  With sw_out_self the converged self-energy is written
    twice: 'self_en.npz' (mesh bound, for a restart on the same grid) and the
    mesh-independent Wannier form via output_self_wannier.

    @param Nx, Ny, Nz: k-mesh size
    @param      Nw: number of positive Matsubara frequencies
    @param   ham_r: hopping matrix H(R) [Nr, Norb, Norb] complex128
    @param     S_r: overlap matrix S(R); [] for an orthogonal basis
    @param    rvec: R-vectors [Nr, 3]
    @param      mu: starting chemical potential [eV]
    @param    temp: temperature k_B*T [eV]
    @param   olist: orbital-pair basis of chi [Nchi, 2] int64
    @param    site: site index of each orbital [Norb] int64
    @param orb_dep: True uses the orbital-dependent Umat/Jmat, False the scalar U/J
    @param    U, J: on-site Hubbard U and Hund J [eV] (U' = U - 2J is used internally)
    @param    fill: band filling that fixes mu each iteration
    @param sw_out_self: write the converged self-energy to disk
    @param sw_in_self: start from the stored self-energy instead of Sigma = 0
    @param Umat, Jmat: orbital-resolved U and J [Norb, Norb] (orb_dep only)
    @param     eps: convergence tolerance of the self-consistency loop
    @param      pp: linear mixing weight of the new self-energy
    @param  m_diis: DIIS subspace size (0 disables DIIS)
    @param sw_rescale: rescale the interaction to keep the Stoner factor below 1
    @param sw_tail: add the high-frequency tail correction to chi0
    @param sigma_scale: scale factor applied to a self-energy read from file
    @return None
    """
    state = _prepare_kspace_state(Nx, Ny, Nz, ham_r, S_r, rvec, need_ham=True)
    # FLEX uses the irreducible k-mesh Green's function together with the same orbital-pair
    # interaction basis used in the response routines above.
    Smat, Cmat = _prepare_normal_interaction(olist, site, orb_dep, U, J, Umat, Jmat)
    if sw_in_self:
        _auto_regrid_sigma_seed(temp, Nw)
    sigmak, mu_self = flibs.mkself(Smat, Cmat, state['kmap'], state['invk'], olist, state['ham_k'], state['eig'],
                                   state['uni'], mu, fill, temp, Nw, Nx, Ny, Nz, sw_out_self, sw_in_self, eps=eps,
                                   pp=pp, m_diis=m_diis, sw_rescale=sw_rescale, sw_tail=sw_tail,
                                   sigma_scale=sigma_scale)
    if sw_out_self:
        np.savez('self_en', sigmak, mu_self)
        _wannier_io.output_self_wannier(sigmak, mu_self, state['kmap'], state['invk'], Nx, Ny, Nz, Nw, temp)


def calc_flex_soc(Nx: int, Ny: int, Nz: int, Nw: int, ham_r: np.ndarray, S_r: np.ndarray, rvec: np.ndarray, mu: float,
                  temp: float, olist: np.ndarray, slist: np.ndarray, invs: np.ndarray, site: np.ndarray, orb_dep: bool,
                  U: float, J: float, fill: float, sw_out_self: bool, sw_in_self: bool, Umat: np.ndarray | None = None,
                  Jmat: np.ndarray | None = None, eps: float = 1.0e-4, pp: float = 0.5, m_diis: int = 5,
                  sw_rescale: bool = True, sigma_scale: float = 1.0) -> None:
    """
    @brief Solve the FLEX self-energy loop WITH spin-orbit coupling.

    Same loop as calc_flex, but the spin-resolved basis makes the spin and charge
    channels mix, so the separate S/C matrices are replaced by the single SOC
    vertex Vmat and the Fortran entry point is mkself_soc.  Output handling is
    identical to calc_flex.

    @param Nx, Ny, Nz: k-mesh size
    @param      Nw: number of positive Matsubara frequencies
    @param   ham_r: hopping matrix H(R) [Nr, Norb, Norb] complex128 (spin-resolved basis)
    @param     S_r: overlap matrix S(R); [] for an orthogonal basis
    @param    rvec: R-vectors [Nr, 3]
    @param      mu: starting chemical potential [eV]
    @param    temp: temperature k_B*T [eV]
    @param   olist: orbital-pair basis of chi [Nchi, 2] int64
    @param   slist: spin label of each orbital [Norb] int64
    @param    invs: spin-flip index map [Norb] int64
    @param    site: site index of each orbital [Norb] int64
    @param orb_dep: True uses the orbital-dependent Umat/Jmat, False the scalar U/J
    @param    U, J: on-site Hubbard U and Hund J [eV]
    @param    fill: band filling that fixes mu each iteration (total, not per spin)
    @param sw_out_self: write the converged self-energy to disk
    @param sw_in_self: start from the stored self-energy instead of Sigma = 0
    @param Umat, Jmat: orbital-resolved U and J [Norb, Norb] (orb_dep only)
    @param     eps: convergence tolerance of the self-consistency loop
    @param      pp: linear mixing weight of the new self-energy
    @param  m_diis: DIIS subspace size (0 disables DIIS)
    @param sw_rescale: rescale the interaction to keep the Stoner factor below 1
    @param sigma_scale: scale factor applied to a self-energy read from file
    @return None
    """
    state = _prepare_kspace_state(Nx, Ny, Nz, ham_r, S_r, rvec, need_ham=True)
    Vmat = _prepare_soc_interaction(olist, slist, site, invs, orb_dep, U, J, Umat, Jmat)
    if sw_in_self:
        _auto_regrid_sigma_seed(temp, Nw)
    sigmak, mu_self = flibs.mkself_soc(Vmat, state['kmap'], state['invk'], invs, olist, slist, state['ham_k'],
                                       state['eig'], state['uni'], mu, fill, temp, Nw, Nx, Ny, Nz, sw_out_self,
                                       sw_in_self, eps=eps, pp=pp, m_diis=m_diis, sw_rescale=sw_rescale,
                                       sigma_scale=sigma_scale)
    if sw_out_self:
        np.savez('self_en', sigmak, mu_self)
        _wannier_io.output_self_wannier(sigmak, mu_self, state['kmap'], state['invk'], Nx, Ny, Nz, Nw, temp)


def calc_lin_eliashberg_eq(Nx: int, Ny: int, Nz: int, Nw: int, ham_r: np.ndarray, S_r: np.ndarray, rvec: np.ndarray,
                           chiolist: np.ndarray, site: np.ndarray, plist: np.ndarray, mu: float, temp: float,
                           gap_sym: int, sw_self: bool, orb_dep: bool, U: float, J: float, fill: float,
                           sw_from_file: bool, sw_out_self: bool, sw_in_self: bool, Umat: np.ndarray | None = None,
                           Jmat: np.ndarray | None = None, eps: float = 1.0e-4, pp: float = 0.5, m_diis: int = 5,
                           sw_rescale: bool = True, sw_tail: bool = False, sigma_scale: float = 1.0,
                           sw_gap_extrapolate: bool = False, iperm: np.ndarray | None = None,
                           gap_extrap_points: int = 4, gap_extrap_order: int = 1) -> None:
    """
    @brief Solve the linearized Eliashberg equation (no SOC) and report lambda.

    Runs the whole chain: k-space state -> Coulomb vertex -> Green's function
    (optionally through the FLEX self-energy when sw_self) -> chi0 -> chi_s/chi_c
    -> pairing kernel -> leading eigenvalue.  lambda_eliash = 1 marks Tc, so the
    printed value at a given temperature says how far the model is from ordering;
    the Stoner factor printed alongside it warns how close chi_s is to diverging.

    Side effects: 'chis.dat' and 'chic.dat' receive the kz=0 plane of the leading
    eigenvalue of chi_s and chi_c; the gap is written by output_gap_function, and
    with sw_out_self also as 'gap.npy' plus the mesh-independent Wannier form.
    Returns early (without writing) when the self-energy could not be prepared.

    @param Nx, Ny, Nz: k-mesh size
    @param      Nw: number of positive Matsubara frequencies
    @param   ham_r: hopping matrix H(R) [Nr, Norb, Norb] complex128
    @param     S_r: overlap matrix S(R); [] for an orthogonal basis
    @param    rvec: R-vectors [Nr, 3]
    @param chiolist: orbital-pair basis of chi [Nchi, 2] int64
    @param    site: site index of each orbital [Norb] int64
    @param   plist: parity sign (+1/-1) per orbital [Norb] float64
    @param      mu: chemical potential [eV]
    @param    temp: temperature k_B*T [eV]
    @param gap_sym: symmetry index of the initial gap; its sign selects singlet/triplet
    @param sw_self: True solves FLEX first, False uses the bare Green's function
    @param orb_dep: True uses the orbital-dependent Umat/Jmat, False the scalar U/J
    @param    U, J: on-site Hubbard U and Hund J [eV]
    @param    fill: band filling
    @param sw_from_file: read the self-energy from disk instead of solving FLEX
    @param sw_out_self: write the self-energy and the gap to disk
    @param sw_in_self: seed the FLEX loop from the stored self-energy
    @param Umat, Jmat: orbital-resolved U and J [Norb, Norb] (orb_dep only)
    @param eps, pp, m_diis, sw_rescale, sw_tail, sigma_scale: FLEX loop controls, see calc_flex
    @param sw_gap_extrapolate: extrapolate Delta(iw_n) -> Delta(0) before output
    @param   iperm: orbital permutation of the inversion operator (None = identity)
    @param gap_extrap_points: number of Matsubara slices used by that extrapolation
    @param gap_extrap_order: polynomial order of that fit
    @return None
    """
    state = _prepare_kspace_state(Nx, Ny, Nz, ham_r, S_r, rvec, need_ham=sw_self)
    Smat, Cmat = _prepare_normal_interaction(chiolist, site, orb_dep, U, J, Umat, Jmat)
    green_state = _prepare_green_state_normal(state, chiolist, (Smat, Cmat), mu, fill, temp, Nw,
                                              Nx, Ny, Nz, sw_self, sw_from_file, sw_out_self,
                                              sw_in_self, eps, pp, m_diis, sw_rescale, sw_tail,
                                              sigma_scale)
    if green_state is None:
        return
    Gk = green_state['Gk']
    init_delta = _response.get_initial_gap(state['klist'], len(state['eig'].T), gap_sym)
    sw_eig = True  # True: leading eigenvalue (most-unstable mode); False: trace
    chi, stoner = _get_chi0_normal(state, Gk, chiolist, Smat, Cmat, green_state['mu_self'],
                                   temp, Nx, Ny, Nz, sw_tail)
    chis, chic = flibs.get_chis_chic(chi, Smat, Cmat)
    chisq = flibs.get_eig_or_tr_chi(chis, state['invk'], sw_eig)
    chicq = flibs.get_eig_or_tr_chi(chic, state['invk'], sw_eig)
    try:
        with open('chis.dat', 'w') as f, open('chic.dat', 'w') as f2:
            for i, k in enumerate(state['kmap']):
                if k[2] == 0.0:
                    f.write(f'{k[0]:6.4f} {k[1]:6.4f} {chisq[i].real:11.4e}\n')
                    f2.write(f'{k[0]:6.4f} {k[1]:6.4f} {chicq[i].real:11.4e}\n')
    except OSError as e:
        print(f"Error: Failed to write 'chis.dat' or 'chic.dat': {e}", flush=True)
    gap, lambda_eliash = flibs.linearized_eliashberg(chi, Gk, state['uni'], init_delta, Smat, Cmat, chiolist, plist,
                                                     state['kmap'], state['invk'], Nx, Ny, Nz, temp, gap_sym,
                                                     iperm=iperm)
    print(f'Stoner factor = {stoner:.6f}, lambda_eliash = {lambda_eliash:.6f}', flush=True)
    if sw_out_self:
        np.save('gap', gap)
        _wannier_io.output_gap_wannier(gap, state['kmap'], state['invk'], Nx, Ny, Nz, Nw, temp)
    output_gap_function(state['invk'], state['kmap'], gap, state['uni'], plist, gap_sym, Nx, iperm=iperm,
                        sw_extrapolate=sw_gap_extrapolate, temp=temp,
                        n_points=gap_extrap_points, order=gap_extrap_order)


def calc_lin_eliash_soc(Nx: int, Ny: int, Nz: int, Nw: int, ham_r: np.ndarray, S_r: np.ndarray, rvec: np.ndarray,
                        mu: float, temp: float, chiolist: np.ndarray, slist: np.ndarray, plist: np.ndarray,
                        invs: np.ndarray, site: np.ndarray, orb_dep: bool, U: float, J: float, fill: float,
                        gap_sym: int, sw_self: bool, sw_from_file: bool, sw_out_self: bool, sw_in_self: bool,
                        Umat: np.ndarray | None = None, Jmat: np.ndarray | None = None, eps: float = 1.0e-4,
                        pp: float = 0.5, m_diis: int = 5, sw_rescale: bool = True, sigma_scale: float = 1.0,
                        sw_gap_extrapolate: bool = False, gap_extrap_points: int = 4,
                        gap_extrap_order: int = 1) -> None:
    """
    @brief Solve the linearized Eliashberg equation WITH spin-orbit coupling.

    Same chain as calc_lin_eliashberg_eq, but in the spin-resolved basis: one SOC
    vertex Vmat replaces the S/C pair, and chi0 comes back with the sign tables
    (sgnsig, sgnsig2) and the inverse orbital map invschi that the SOC pairing
    kernel needs.  Because spin is no longer a good quantum number, chi_s splits
    into a longitudinal (chi_szz) and a transverse (chi_spm) channel and both are
    written to 'chis.dat'; 'chic.dat' holds the charge channel.  Returns early
    when the self-energy could not be prepared.

    lambda_eliash is reported exactly as in the non-SOC path: the Fortran driver
    prints 'lambda_elsh=' every iteration and 'eliash=' on convergence, and hands
    the converged value back through lambda_out, which is printed here.  The only
    remaining difference is that no Stoner factor accompanies it: the SOC chi0
    entry point does not return one, and the quantity is defined differently in
    the spin-resolved basis, so it is deliberately not reported here.

    @param Nx, Ny, Nz: k-mesh size
    @param      Nw: number of positive Matsubara frequencies
    @param   ham_r: hopping matrix H(R) [Nr, Norb, Norb] complex128 (spin-resolved basis)
    @param     S_r: overlap matrix S(R); [] for an orthogonal basis
    @param    rvec: R-vectors [Nr, 3]
    @param      mu: chemical potential [eV]
    @param    temp: temperature k_B*T [eV]
    @param chiolist: orbital-pair basis of chi [Nchi, 2] int64
    @param   slist: spin label of each orbital [Norb] int64
    @param   plist: parity sign (+1/-1) per orbital [Norb] float64
    @param    invs: spin-flip index map [Norb] int64
    @param    site: site index of each orbital [Norb] int64
    @param orb_dep: True uses the orbital-dependent Umat/Jmat, False the scalar U/J
    @param    U, J: on-site Hubbard U and Hund J [eV]
    @param    fill: band filling (total, not per spin)
    @param gap_sym: symmetry index of the initial gap
    @param sw_self: True solves FLEX first, False uses the bare Green's function
    @param sw_from_file: read the self-energy from disk instead of solving FLEX
    @param sw_out_self: write the self-energy and the gap to disk
    @param sw_in_self: seed the FLEX loop from the stored self-energy
    @param Umat, Jmat: orbital-resolved U and J [Norb, Norb] (orb_dep only)
    @param eps, pp, m_diis, sw_rescale, sigma_scale: FLEX loop controls, see calc_flex_soc
    @param sw_gap_extrapolate: extrapolate Delta(iw_n) -> Delta(0) before output
    @param gap_extrap_points: number of Matsubara slices used by that extrapolation
    @param gap_extrap_order: polynomial order of that fit
    @return None
    """
    state = _prepare_kspace_state(Nx, Ny, Nz, ham_r, S_r, rvec, need_ham=sw_self)
    Vmat = _prepare_soc_interaction(chiolist, slist, site, invs, orb_dep, U, J, Umat, Jmat)
    green_state = _prepare_green_state_soc(state, chiolist, slist, invs, Vmat, mu, fill, temp, Nw,
                                           Nx, Ny, Nz, sw_self, sw_from_file, sw_out_self,
                                           sw_in_self, eps, pp, m_diis, sw_rescale,
                                           sigma_scale)
    if green_state is None:
        return
    Gk = green_state['Gk']
    sw_eig = True  # True: leading eigenvalue (most-unstable mode); False: trace
    chi, sgnsig, sgnsig2, invschi = flibs.get_chi0_soc(Vmat, Gk, chiolist, slist, state['kmap'], state['invk'], invs,
                                                       temp, Nx, Ny, Nz)
    chic, chiszz, chispm = flibs.get_chis_chic_soc(chi, Vmat, chiolist, slist, invs)
    chiszzq = flibs.get_eig_or_tr_chi(chiszz, state['invk'], sw_eig)
    chispmq = flibs.get_eig_or_tr_chi(chispm, state['invk'], sw_eig)
    chicq = flibs.get_eig_or_tr_chi(chic, state['invk'], sw_eig)
    try:
        with open('chis.dat', 'w') as f, open('chic.dat', 'w') as f2:
            for i, k in enumerate(state['kmap']):
                if k[2] == 0.0:
                    f.write(f'{k[0]:6.4f} {k[1]:6.4f} {chiszzq[i].real:11.4e} {chispmq[i].real:11.4e}\n')
                    f2.write(f'{k[0]:6.4f} {k[1]:6.4f} {chicq[i].real:11.4e}\n')
    except OSError as e:
        print(f"Error: Failed to write 'chis.dat' or 'chic.dat': {e}", flush=True)
    init_delta = _response.get_initial_gap(state['klist'], len(slist), gap_sym)
    gap, lambda_eliash = flibs.linearized_eliashberg_soc(chi, Gk, state['uni'], init_delta, Vmat, sgnsig, sgnsig2,
                                                         plist, slist, chiolist, state['kmap'], state['invk'],
                                                         invs, invschi, Nx, Ny, Nz, temp, gap_sym)
    print(f'lambda_eliash = {lambda_eliash:.6f}', flush=True)
    if sw_out_self:
        np.save('gap', gap)
        _wannier_io.output_gap_wannier(gap, state['kmap'], state['invk'], Nx, Ny, Nz, Nw, temp)
    output_gap_function(state['invk'], state['kmap'], gap, state['uni'], plist, gap_sym, Nx, True, invs, slist,
                        sw_extrapolate=sw_gap_extrapolate, temp=temp,
                        n_points=gap_extrap_points, order=gap_extrap_order)


def calc_eliashberg_eq(Nx: int, Ny: int, Nz: int, Nw: int, ham_r: np.ndarray, S_r: np.ndarray, rvec: np.ndarray,
                       chiolist: np.ndarray, site: np.ndarray, plist: np.ndarray, mu: float, temp: float, gap_sym: int,
                       sw_self: bool, orb_dep: bool, U: float, J: float, fill: float, sw_from_file: bool,
                       sw_out_self: bool, sw_in_self: bool, Umat: np.ndarray | None = None,
                       Jmat: np.ndarray | None = None, eps: float = 1.0e-4, pp: float = 0.5, m_diis: int = 5,
                       sw_rescale: bool = True, sw_check_only: bool = False, sw_tail: bool = False,
                       sigma_scale: float = 1.0, sw_gap_extrapolate: bool = False, gap_extrap_points: int = 4,
                       gap_extrap_order: int = 1, iperm: np.ndarray | None = None) -> None:
    """
    @brief Solve the nonlinear Eliashberg equation (no SOC), seeded by the linearized solution.

    The linearized problem is solved first, both as a cheap physical gate and as the
    source of a symmetry-correct initial gap.  The run stops before the expensive
    nonlinear loop when the Stoner factor >= 1 (SDW/CDW instability), when
    lambda_eliash < 1 (T >= Tc: no superconducting instability at this temperature),
    or when sw_check_only is set.  Otherwise only the SHAPE of the linearized
    eigenvector is kept: it is rescaled to the weak-coupling amplitude BCS_RATIO*temp
    and handed to flibs.nonlinear_eliashberg, which iterates the gap (and, with
    sw_self, the self-energy) to self-consistency.

    Side effects: the gap is written by output_gap_function, and with sw_out_self also
    as 'gap.npy' plus the mesh-independent Wannier form.  Returns early (without
    writing) when the self-energy could not be prepared.

    @param Nx, Ny, Nz: k-mesh size
    @param      Nw: number of positive Matsubara frequencies
    @param   ham_r: hopping matrix H(R) [Nr, Norb, Norb] complex128
    @param     S_r: overlap matrix S(R); [] for an orthogonal basis
    @param    rvec: R-vectors [Nr, 3]
    @param chiolist: orbital-pair basis of chi [Nchi, 2] int64
    @param    site: site index of each orbital [Norb] int64
    @param   plist: parity sign (+1/-1) per orbital [Norb] float64
    @param      mu: chemical potential [eV]
    @param    temp: temperature k_B*T [eV]
    @param gap_sym: symmetry index of the initial gap; its sign selects singlet/triplet
    @param sw_self: True includes the FLEX self-energy (also inside the nonlinear loop),
                    False uses the bare Green's function
    @param orb_dep: True uses the orbital-dependent Umat/Jmat, False the scalar U/J
    @param    U, J: on-site Hubbard U and Hund J [eV]
    @param    fill: band filling
    @param sw_from_file: read the self-energy from disk instead of solving FLEX
    @param sw_out_self: write the self-energy and the gap to disk
    @param sw_in_self: seed the FLEX loop from the stored self-energy
    @param Umat, Jmat: orbital-resolved U and J [Norb, Norb] (orb_dep only)
    @param eps, pp, m_diis, sw_rescale, sw_tail, sigma_scale: FLEX loop controls, see
                    calc_flex; eps and m_diis also control the nonlinear loop
    @param sw_check_only: If True, stop after the linearized Eliashberg solve (before the
                          nonlinear loop) and report the Stoner factor and lambda_eliash.
                          The calculation also stops early (regardless of this flag) when
                          the Stoner factor >= 1 (SDW/CDW instability) or lambda_eliash < 1
                          (T >= Tc; no superconducting instability at this temperature).
    @param sw_gap_extrapolate: extrapolate Delta(iw_n) -> Delta(0) before output
    @param gap_extrap_points: number of Matsubara slices used by that extrapolation
    @param gap_extrap_order: polynomial order of that fit
    @param   iperm: orbital permutation of the inversion operator (None = identity)
    """
    state = _prepare_kspace_state(Nx, Ny, Nz, ham_r, S_r, rvec, need_ham=True)
    Smat, Cmat = _prepare_normal_interaction(chiolist, site, orb_dep, U, J, Umat, Jmat)
    green_state = _prepare_green_state_normal(state, chiolist, (Smat, Cmat), mu, fill, temp, Nw,
                                              Nx, Ny, Nz, sw_self, sw_from_file, sw_out_self,
                                              sw_in_self, eps, pp, m_diis, sw_rescale, sw_tail,
                                              sigma_scale)
    if green_state is None:
        return
    Gk = green_state['Gk']
    mu_self = green_state['mu_self']
    delta_init_band = _response.get_initial_gap(state['klist'], len(state['eig'].T), gap_sym)
    # Use the linearized Eliashberg solution as a symmetry-correct initial gap (a plain orbital
    # projection, flibs.get_initial_delta, would also work as a seed).  Solving the linearized
    # problem first doubles as a cheap physical gate before the much more expensive nonlinear
    # loop: it filters out SDW/CDW instabilities and T >= Tc cases.
    chi, stoner = _get_chi0_normal(state, Gk, chiolist, Smat, Cmat, mu_self,
                                   temp, Nx, Ny, Nz, sw_tail)
    if stoner >= 1.0:
        print(f"Stoner factor = {stoner:.6f} >= 1: SDW/CDW instability, stopping before "
              f"nonlinear Eliashberg.", flush=True)
        return
    delta_init, lambda_eliash = flibs.linearized_eliashberg(chi, Gk, state['uni'], delta_init_band, Smat, Cmat,
                                                            chiolist, plist, state['kmap'], state['invk'], Nx, Ny,
                                                            Nz, temp, gap_sym, iperm=iperm)
    print(f'Stoner factor = {stoner:.6f}, lambda_eliash = {lambda_eliash:.6f}', flush=True)
    if lambda_eliash < 1.0:
        print(f"lambda_eliash = {lambda_eliash:.6f} < 1: no superconducting instability "
              f"at this temperature (T >= Tc), stopping before nonlinear Eliashberg.", flush=True)
        return
    if sw_check_only:
        print("sw_check_only=True: stopping after linearized Eliashberg.", flush=True)
        return
    # BCS weak-coupling ratio: Δ₀ = 1.764 kB Tc; scale initial gap to this amplitude.
    target_gap = _eilenberger.BCS_RATIO * temp
    delta_max = np.abs(delta_init).max()
    if delta_max > 0.0:
        # Only the shape from the linearized eigenvector matters here; the nonlinear solver
        # converges more reliably if the starting amplitude is normalized to the target scale.
        delta_init *= target_gap / delta_max
    else:
        print("Warning: initial gap is zero; skip Tc-based scaling", flush=True)
    delta, sigmak = flibs.nonlinear_eliashberg(delta_init, Gk, state['ham_k'], Smat, Cmat, chiolist, plist,
                                               state['kmap'], state['invk'], mu_self, temp, gap_sym, Nx, Ny, Nz,
                                               sw_sigma=sw_self, sw_Vconst=True, eps=eps, m_diis=m_diis, iperm=iperm)
    if sw_out_self:
        np.save('gap', delta)
        _wannier_io.output_gap_wannier(delta, state['kmap'], state['invk'], Nx, Ny, Nz, Nw, temp)
    output_gap_function(state['invk'], state['kmap'], delta, state['uni'], plist, gap_sym, Nx, iperm=iperm,
                        sw_extrapolate=sw_gap_extrapolate, temp=temp,
                        n_points=gap_extrap_points, order=gap_extrap_order)
