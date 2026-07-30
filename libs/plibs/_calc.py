#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Higher-level calculation routines: FLEX, Eliashberg, susceptibility spectra, carrier analysis.
"""
import numpy as np
import libs.flibs as flibs
from ._bands import get_eigs, get_emesh, calc_mu
from ._eilenberger import BCS_RATIO
from ._response import get_initial_gap, chis_spectrum, phi_spectrum
from ._wannier_io import output_self_wannier, output_gap_wannier


def _prepare_kspace_state(Nx:int, Ny:int, Nz:int, ham_r, S_r, rvec, need_ham:bool=False):
    klist, kmap, invk = flibs.gen_irr_k_TRS(Nx, Ny, Nz)
    eig, uni = get_eigs(klist, ham_r, S_r, rvec)
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


def _prepare_normal_interaction(olist, site, orb_dep:bool, U:float, J:float, Umat=None, Jmat=None):
    if orb_dep:
        return flibs.gen_SCmatrix_orb(olist, site, Umat, Jmat)
    return flibs.gen_SCmatrix(olist, site, U, J)


def _prepare_soc_interaction(olist, slist, site, invs, orb_dep:bool, U:float, J:float, Umat=None, Jmat=None):
    if orb_dep:
        return flibs.gen_Vmatrix_orb(olist, slist, site, invs, Umat, Jmat)
    return flibs.gen_Vmatrix(olist, slist, site, invs, U, J)


def _load_sigma_from_file():
    """Load (sigmak, mu_self) from self_en.npz, or None if the file is missing."""
    try:
        npz = np.load('self_en.npz')
        print("Import sigma from self_en.npz")
        return npz['arr_0'], npz['arr_1']
    except FileNotFoundError:
        print("Error: 'self_en.npz' not found", flush=True)
        return None


def _read_sigma_bin_meta(path:str='sigma.bin'):
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


def _auto_regrid_sigma_seed(temp:float, Nw:int, path:str='sigma.bin'):
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


def regrid_sigma_bin(temp_old:float=None, temp_new:float=None, Norb:int=None,
                     Nw_old:int=None, Nw_new:int=None, w_scale:float=None,
                     path:str='sigma.bin', out_path:str=None):
    """Re-grid the FLEX self-energy seed sigma.bin from the Matsubara mesh at
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
    from scipy.io import FortranFile
    from scipy.interpolate import CubicSpline
    if temp_new is None:
        raise ValueError("temp_new is required")
    if w_scale is None:
        w_scale = 1.0
    with FortranFile(path, 'r') as f:
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
    out = CubicSpline(wext, ext, axis=1, extrapolate=False)(w_new)
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
    with FortranFile(out_path or path, 'w') as f:
        f.write_record(mu)
        f.write_record(mu_old)
        f.write_record(np.ascontiguousarray(out.ravel(order='F')))
        f.write_record(np.array([temp_new, float(Nw_new), float(Norb), float(Nk)]))
    print(f"regrid_sigma_bin: sigma.bin re-gridded T={temp_old:.6f} -> {temp_new:.6f} eV, "
          f"Nw={Nw_old} -> {Nw_new} (Nk={Nk}, Norb={Norb})", flush=True)
    return Nk, Nw_new


def _prepare_green_state_normal(state, olist, interaction, mu:float, fill:float, temp:float,
                                Nw:int, Nx:int, Ny:int, Nz:int, sw_self:bool,
                                sw_from_file:bool=False, sw_out_self:bool=False, sw_in_self:bool=False,
                                eps:float=1.0e-4, pp:float=0.5, m_diis:int=5, sw_rescale:bool=True,
                                sw_tail:bool=False, sigma_scale:float=1.0):
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


def _prepare_green_state_soc(state, olist, slist, invs, interaction, mu:float, fill:float, temp:float,
                             Nw:int, Nx:int, Ny:int, Nz:int, sw_self:bool,
                             sw_from_file:bool=False, sw_out_self:bool=False, sw_in_self:bool=False,
                             eps:float=1.0e-4, pp:float=0.5, m_diis:int=5, sw_rescale:bool=True,
                             sigma_scale:float=1.0):
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

def _get_chi0_normal(state, Gk, chiolist, Smat, Cmat, mu_self, temp, Nx, Ny, Nz, sw_tail):
    """Irreducible chi0 on the irreducible q-grid; sw_tail=True uses the
    tail-corrected convolution (conv[G]-conv[G0]+analytic reference, needs the
    chemical potential mu_self of Gk)."""
    if sw_tail:
        return flibs.get_chi0_tail(Smat, Cmat, Gk, state['eig'], state['uni'], chiolist,
                                   state['kmap'], state['invk'], mu_self, temp, Nx, Ny, Nz)
    return flibs.get_chi0(Smat, Cmat, Gk, chiolist, state['kmap'], state['invk'], temp, Nx, Ny, Nz)


def calc_path_spectrum(kind:str, mu:float, temp:float, klist, qlist, chiolist, eig, uni,
                       spa_length, Nw:int, Emax:float, delta:float, Smat=None):
    if kind == 'chis':
        print("calculate spin susceptibility", flush=True)
        spec, spec_orb, wlist = chis_spectrum(mu, temp, Smat, klist, qlist, chiolist, eig, uni, Nw, Emax, delta)
        outname = 'chis.dat'
    elif kind == 'phi':
        print("calculate sc susceptibility", flush=True)
        spec, spec_orb, wlist = phi_spectrum(mu, temp, klist, qlist, chiolist, eig, uni, Nw, Emax, delta)
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
    except IOError as e:
        print(f"Error: Failed to write '{outname}': {e}", flush=True)

    if kind == 'chis':
        for i, spec_orb_row in enumerate(spec_orb.T):
            try:
                with open(f'chis_{i}.dat', 'w') as f:
                    for ww, ssp, spec_band in zip(w, sp, spec_orb_row.T):
                        for www, sssp, spec_val in zip(ww, ssp, spec_band):
                            f.write(f'{sssp:8.4f} {www:8.4f} {spec_val.imag:9.4f}\n')
                        f.write('\n')
            except IOError as e:
                print(f"Error: Failed to write 'chis_{i}.dat': {e}", flush=True)
                continue

    return w, sp, spec

def get_carrier_num(kmesh, rvec, ham_r, S_r, mu:float, Arot):
    Nk, eig, kweight = get_emesh(kmesh, kmesh, kmesh, ham_r, S_r, rvec, Arot)
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

def get_mu(ham_r, S_r, rvec, Arot, temp:float, fill:float, kmesh=40) -> float:
    if temp <= 0:
        print("Error: Temperature (temp) is non-positive", flush=True)
        return None
    if kmesh <= 0:
        print("Error: k-mesh size (kmesh) is non-positive", flush=True)
        return None
    print("calc chem. pot.", flush=True)
    print(f"band filling = {fill:.4f}", flush=True)
    Nk, eig, kweight = get_emesh(kmesh, kmesh, kmesh, ham_r, S_r, rvec, Arot)
    mu = calc_mu(eig, Nk, fill, temp)
    return mu

def gap_extrapolate_w0(gap: np.ndarray, temp: float, n_points: int = 4, order: int = 1):
    """
    @fn gap_extrapolate_w0
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
    @retval  gap0: extrapolated Delta(k,w_n->0) [Norb, Norb, Nk] complex128
    @retval  bias: |gap0 - Delta(iw_0)| / max|gap0|, the O((pi*temp)^2) correction
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
    data = np.transpose(gap[:, :, :n_points, :], (2,0,1,3)).reshape(n_points, -1)  # [n_points, Norb^2*Nk]
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

def output_gap_function(invk, kmap, gap, uni, plist, gap_sym, Nx:int,
                        soc=False, invs=None, slist=None, sw_orb=False,
                        sw_extrapolate=False, temp=None, n_points=4, order=1, iperm=None):
    if sw_extrapolate:
        gap0, bias = gap_extrapolate_w0(gap, temp, n_points, order)
        print(f'gap w_n->0 extrapolation: max relative correction over iw_0 = {bias.max():.4e}', flush=True)
        gap = gap0[:, :, None, :]  # restore a length-1 Matsubara axis; downstream code reads index 0
    if sw_orb:
        if soc:
            gapb = gap[:, :, 0, :]
        else:
            gapb = flibs.remap_gap(gap[:, :, 0, :], plist, invk, gap_sym, iperm=iperm)
    else:
        if soc:
            gapb = flibs.conv_delta_orb_to_band_soc(gap, uni, invk, invs, slist)
        else:
            gapb = flibs.conv_delta_orb_to_band(gap, uni, invk, plist, gap_sym, iperm=iperm)
    print('output gap function')
    for iorb in range(len(gapb)):
        for jorb in range(len(gapb)):
            try:
                with open(f'gap_{iorb+1}{jorb+1}.dat', 'w') as f:
                    for i, km in enumerate(kmap):
                        if km[2] == 0:
                            f.write(f'{km[0]:3} {km[1]:3} {gapb[iorb,jorb,i].real:15.8e} {gapb[iorb,jorb,i].imag:15.8e}\n')
                            if km[0] == Nx - 1:
                                f.write('\n')
            except IOError as e:
                print(f"Error: Failed to write 'gap_{iorb+1}{jorb+1}.dat': {e}", flush=True)
                continue
    return 0

def calc_flex(Nx:int, Ny:int, Nz:int, Nw:int, ham_r, S_r, rvec, mu:float, temp:float,
              olist, site, orb_dep:bool, U:float, J:float, fill:float, sw_out_self:bool, sw_in_self:bool, 
              Umat=None, Jmat=None, eps=1.0e-4, pp=0.5, m_diis=5, sw_rescale:bool=True, sw_tail:bool=False,
              sigma_scale:float=1.0):
    state = _prepare_kspace_state(Nx, Ny, Nz, ham_r, S_r, rvec, need_ham=True)
    # FLEX uses the irreducible k-mesh Green's function together with the same orbital-pair
    # interaction basis used in the response routines above.
    Smat, Cmat = _prepare_normal_interaction(olist, site, orb_dep, U, J, Umat, Jmat)
    if sw_in_self:
        _auto_regrid_sigma_seed(temp, Nw)
    sigmak, mu_self = flibs.mkself(Smat, Cmat, state['kmap'], state['invk'], olist, state['ham_k'], state['eig'], state['uni'],
                                   mu, fill, temp, Nw, Nx, Ny, Nz, sw_out_self, sw_in_self,
                                   eps=eps, pp=pp, m_diis=m_diis, sw_rescale=sw_rescale, sw_tail=sw_tail,
                                   sigma_scale=sigma_scale)
    if sw_out_self:
        np.savez('self_en', sigmak, mu_self)
        output_self_wannier(sigmak, mu_self, state['kmap'], state['invk'], Nx, Ny, Nz, Nw, temp)

def calc_flex_soc(Nx:int, Ny:int, Nz:int, Nw:int, ham_r, S_r, rvec, mu:float, temp:float,
                  olist, slist, invs, site,
                  orb_dep:bool, U:float, J:float, fill:float,
                  sw_out_self:bool, sw_in_self:bool,
                  Umat=None, Jmat=None, eps=1.0e-4, pp=0.5, m_diis=5, sw_rescale:bool=True,
                  sigma_scale:float=1.0):
    state = _prepare_kspace_state(Nx, Ny, Nz, ham_r, S_r, rvec, need_ham=True)
    Vmat = _prepare_soc_interaction(olist, slist, site, invs, orb_dep, U, J, Umat, Jmat)
    if sw_in_self:
        _auto_regrid_sigma_seed(temp, Nw)
    sigmak, mu_self = flibs.mkself_soc(Vmat, state['kmap'], state['invk'], invs, olist, slist, state['ham_k'], state['eig'], state['uni'],
                                       mu, fill, temp, Nw, Nx, Ny, Nz, sw_out_self, sw_in_self,
                                       eps=eps, pp=pp, m_diis=m_diis, sw_rescale=sw_rescale,
                                       sigma_scale=sigma_scale)
    if sw_out_self:
        np.savez('self_en', sigmak, mu_self)
        output_self_wannier(sigmak, mu_self, state['kmap'], state['invk'], Nx, Ny, Nz, Nw, temp)

def calc_lin_eliashberg_eq(Nx:int, Ny:int, Nz:int, Nw:int, ham_r, S_r, rvec, chiolist, site, plist,
                           mu:float, temp:float, gap_sym:int, sw_self:bool, orb_dep:bool, U:float, J:float,
                           fill:float, sw_from_file:bool, sw_out_self:bool, sw_in_self:bool,
                           Umat=None, Jmat=None, eps=1.0e-4, pp=0.5, m_diis=5, sw_rescale:bool=True,
                           sw_tail:bool=False, sigma_scale:float=1.0, sw_gap_extrapolate:bool=False, iperm=None,
                           gap_extrap_points=4, gap_extrap_order=1):
    state = _prepare_kspace_state(Nx, Ny, Nz, ham_r, S_r, rvec, need_ham=sw_self)
    Smat, Cmat = _prepare_normal_interaction(chiolist, site, orb_dep, U, J, Umat, Jmat)
    green_state = _prepare_green_state_normal(state, chiolist, (Smat, Cmat), mu, fill, temp, Nw,
                                              Nx, Ny, Nz, sw_self, sw_from_file, sw_out_self,
                                              sw_in_self, eps, pp, m_diis, sw_rescale, sw_tail,
                                              sigma_scale)
    if green_state is None:
        return
    Gk = green_state['Gk']
    init_delta = get_initial_gap(state['klist'], len(state['eig'].T), gap_sym)
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
    except IOError as e:
        print(f"Error: Failed to write 'chis.dat' or 'chic.dat': {e}", flush=True)
    gap, lambda_eliash = flibs.linearized_eliashberg(chi, Gk, state['uni'], init_delta, Smat, Cmat, chiolist, plist,
                                      state['kmap'], state['invk'], Nx, Ny, Nz, temp, gap_sym, iperm=iperm)
    print(f'Stoner factor = {stoner:.6f}, lambda_eliash = {lambda_eliash:.6f}', flush=True)
    if sw_out_self:
        np.save('gap', gap)
        output_gap_wannier(gap, state['kmap'], state['invk'], Nx, Ny, Nz, Nw, temp)
    output_gap_function(state['invk'], state['kmap'], gap, state['uni'], plist, gap_sym, Nx, iperm=iperm,
                        sw_extrapolate=sw_gap_extrapolate, temp=temp,
                        n_points=gap_extrap_points, order=gap_extrap_order)

def calc_lin_eliash_soc(Nx:int, Ny:int, Nz:int, Nw:int, ham_r, S_r, rvec, mu:float, temp:float,
                        chiolist, slist, plist, invs, site, orb_dep:bool, U:float, J:float, fill:float,
                        gap_sym:int, sw_self:bool, sw_from_file:bool, sw_out_self:bool, sw_in_self:bool,
                        Umat=None, Jmat=None, eps=1.0e-4, pp=0.5, m_diis=5, sw_rescale:bool=True,
                        sigma_scale:float=1.0, sw_gap_extrapolate:bool=False,
                        gap_extrap_points=4, gap_extrap_order=1):
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
    chi, sgnsig, sgnsig2, invschi = flibs.get_chi0_soc(Vmat, Gk, chiolist, slist, state['kmap'], state['invk'], invs, temp, Nx, Ny, Nz)
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
    except IOError as e:
        print(f"Error: Failed to write 'chis.dat' or 'chic.dat': {e}", flush=True)
    init_delta = get_initial_gap(state['klist'], len(slist), gap_sym)
    gap = flibs.linearized_eliashberg_soc(chi, Gk, state['uni'], init_delta, Vmat, sgnsig, sgnsig2, plist, slist, chiolist,
                                          state['kmap'], state['invk'], invs, invschi, Nx, Ny, Nz, temp, gap_sym)
    if sw_out_self:
        np.save('gap', gap)
        output_gap_wannier(gap, state['kmap'], state['invk'], Nx, Ny, Nz, Nw, temp)
    output_gap_function(state['invk'], state['kmap'], gap, state['uni'], plist, gap_sym, Nx, True, invs, slist,
                        sw_extrapolate=sw_gap_extrapolate, temp=temp,
                        n_points=gap_extrap_points, order=gap_extrap_order)

def calc_eliashberg_eq(Nx:int, Ny:int, Nz:int, Nw:int, ham_r, S_r, rvec,
                       chiolist, site, plist, mu:float, temp:float, gap_sym:int, sw_self:bool,
                       orb_dep:bool, U:float, J:float, fill:float, sw_from_file:bool, sw_out_self:bool,
                       sw_in_self:bool, Umat=None, Jmat=None, eps=1.0e-4, pp=0.5, m_diis=5, sw_rescale:bool=True,
                       sw_check_only:bool=False, sw_tail:bool=False, sigma_scale:float=1.0,
                       sw_gap_extrapolate:bool=False, gap_extrap_points=4, gap_extrap_order=1, iperm=None):
    """
    @param sw_check_only: If True, stop after the linearized Eliashberg solve (before the
                          nonlinear loop) and report the Stoner factor and lambda_eliash.
                          The calculation also stops early (regardless of this flag) when
                          the Stoner factor >= 1 (SDW/CDW instability) or lambda_eliash < 1
                          (T >= Tc; no superconducting instability at this temperature).
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
    delta_init_band = get_initial_gap(state['klist'], len(state['eig'].T), gap_sym)
    # Use linearized Eliashberg to obtain a symmetry-correct initial gap.
    # Alternative: flibs.get_initial_delta (plain orbital projection, kept for reference).
    if True:
        # Reuse the linearized problem as a cheap physical gate before entering the much more
        # expensive nonlinear loop. This filters out SDW/CDW instabilities and T >= Tc cases.
        chi, stoner = _get_chi0_normal(state, Gk, chiolist, Smat, Cmat, mu_self,
                                       temp, Nx, Ny, Nz, sw_tail)
        if stoner >= 1.0:
            print(f"Stoner factor = {stoner:.6f} >= 1: SDW/CDW instability, stopping before "
                  f"nonlinear Eliashberg.", flush=True)
            return
        delta_init, lambda_eliash = flibs.linearized_eliashberg(chi, Gk, state['uni'], delta_init_band, Smat, Cmat,
                                                 chiolist, plist, state['kmap'], state['invk'], Nx, Ny, Nz, temp, gap_sym,
                                                 iperm=iperm)
        print(f'Stoner factor = {stoner:.6f}, lambda_eliash = {lambda_eliash:.6f}', flush=True)
        if lambda_eliash < 1.0:
            print(f"lambda_eliash = {lambda_eliash:.6f} < 1: no superconducting instability "
                  f"at this temperature (T >= Tc), stopping before nonlinear Eliashberg.", flush=True)
            return
        if sw_check_only:
            print("sw_check_only=True: stopping after linearized Eliashberg.", flush=True)
            return
    else:
        delta_init_band = get_initial_gap(state['klist'], len(state['eig'].T), gap_sym)
        delta_init = flibs.get_initial_delta(delta_init_band, state['uni'], state['kmap'], state['invk'], Nw, gap_sym)
    # BCS weak-coupling ratio: Δ₀ = 1.764 kB Tc; scale initial gap to this amplitude.
    target_gap = BCS_RATIO * temp
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
        output_gap_wannier(delta, state['kmap'], state['invk'], Nx, Ny, Nz, Nw, temp)
    output_gap_function(state['invk'], state['kmap'], delta, state['uni'], plist, gap_sym, Nx, iperm=iperm,
                        sw_extrapolate=sw_gap_extrapolate, temp=temp,
                        n_points=gap_extrap_points, order=gap_extrap_order)
