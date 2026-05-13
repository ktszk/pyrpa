#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
EPA output reading and Wannier-R space file I/O for self-energy and gap functions.
"""
import numpy as np


def read_epa_output(filename: str) -> tuple[int,int,np.ndarray,np.ndarray,np.ndarray,np.ndarray,np.ndarray]:
    """
    @fn read_epa_output
    @brief Read epa.x (job='egrid') output file and return EPA data arrays.
    @param filename: Path to the epa.x output file
    @retval   ngrid: Number of energy grids (typically 2: valence, conduction)
    @retval  nmodes: Number of phonon modes
    @retval    edge: Grid edges [ngrid] (eV)
    @retval    step: Grid steps [ngrid] (eV)
    @retval    nbin: Number of bins per grid [ngrid] int64
    @retval    wavg: Averaged phonon frequencies [nmodes] (eV)
    @retval    gavg: EPA averaged |g|^2 [ngrid, nbin_max, nbin_max, nmodes] (eV^2)
    """
    cm2ev = 1.23981e-4
    with open(filename) as f:
        ngrid, nmodes = map(int, f.readline().split())
        edge = np.zeros(ngrid)
        step = np.zeros(ngrid)
        nbin = np.zeros(ngrid, dtype=np.int64)
        for ii in range(ngrid):
            tokens = f.readline().split()
            edge[ii] = float(tokens[0])
            step[ii] = float(tokens[1])
            nbin[ii] = int(tokens[2])
        wavg_cm = np.array(f.readline().split(), dtype=np.float64)
        wavg = wavg_cm * cm2ev  # cm^-1 -> eV
        nbin_max = int(np.max(nbin))
        gavg = np.zeros((ngrid, nbin_max, nbin_max, nmodes), dtype=np.float64)
        for line in f:
            tokens = line.split()
            if len(tokens) < 3 + nmodes:
                continue
            ii = int(tokens[0]) - 1   # 0-based
            jj = int(tokens[1]) - 1
            kk = int(tokens[2]) - 1
            g2 = np.array(tokens[3:3+nmodes], dtype=np.float64)
            gavg[ii, kk, jj, :] = g2  # eV^2
    return ngrid, nmodes, edge, step, nbin, wavg, gavg

def _irr_to_full_kgrid(data_irr: np.ndarray, invk: np.ndarray, kmap: np.ndarray,
                       Nx: int, Ny: int, Nz: int) -> np.ndarray:
    """
    @fn _irr_to_full_kgrid
    @brief Expand data from irreducible k-points to the full BZ using TRS symmetry.
    @param data_irr: (..., Nk_irr) complex128
    @param     invk: [Nkall, 3] int64 — col-0: Fortran 1-based irr-k index, col-1: 0=direct / !=0=TRS conjugate
    @param     kmap: [Nkall, 3] int64 — FFT grid indices (ix, iy, iz), 0-based
    @return    full: (..., Nx, Ny, Nz) complex128
    """
    ik_arr             = invk[:, 0] - 1        # Fortran 1-based → Python 0-based
    is_tr              = invk[:, 1] != 0       # True when this BZ point is reached via TRS: k -> -k
    ix_arr, iy_arr, iz_arr = kmap[:, 0], kmap[:, 1], kmap[:, 2]
    full = np.zeros(data_irr.shape[:-1] + (Nx, Ny, Nz), dtype=np.complex128)
    d = np.where(~is_tr)[0]
    t = np.where( is_tr)[0]
    # Direct k-points: copy irreducible value unchanged
    full[..., ix_arr[d], iy_arr[d], iz_arr[d]] = data_irr[..., ik_arr[d]]
    # TRS-related k-points: apply complex conjugation (TRS: f(k) = f*(-k))
    full[..., ix_arr[t], iy_arr[t], iz_arr[t]] = np.conj(data_irr[..., ik_arr[t]])
    return full

def _wannier_all_nonzero(data_full: np.ndarray, Nx: int, Ny: int, Nz: int,
                         N_cut: int, zero_tol: float) -> tuple[np.ndarray, np.ndarray]:
    """
    @fn _wannier_all_nonzero
    @brief IFFT over k-axes, apply zero_tol zeroing, return all R-vectors with non-zero components.
    @param data_full: (Norb, Norb, Nw, Nx, Ny, Nz) complex128
    @param     N_cut: number of Matsubara frequencies to retain
    @param  zero_tol: Re/Im below zero_tol * global_max are set to 0
    @retval data_out: (Norb, Norb, N_cut, Nr_nonzero) complex128
    @retval rvec_kept: (Nr_nonzero, 3) int64 — centered R-vector coordinates
    """
    N_cut_use = min(N_cut, data_full.shape[2])
    # IFFT over k-mesh (last 3 axes) to obtain Wannier R-space representation
    data_R    = np.fft.ifftn(data_full, axes=(-3, -2, -1))[:, :, :N_cut_use]
    abs_tol   = np.abs(data_R).max() * zero_tol
    # Zero out numerically negligible Re/Im parts to reduce output file size
    data_R.real[np.abs(data_R.real) < abs_tol] = 0.0
    data_R.imag[np.abs(data_R.imag) < abs_tol] = 0.0

    nonzero          = np.abs(data_R).max(axis=(0, 1, 2)) > 0   # (Nx, Ny, Nz)
    ix_nz, iy_nz, iz_nz = np.where(nonzero)

    # Wrap FFT indices to centered Wigner-Seitz convention: indices > N/2 map to negative R
    ix_c = np.where(ix_nz <= Nx // 2, ix_nz, ix_nz - Nx)
    iy_c = np.where(iy_nz <= Ny // 2, iy_nz, iy_nz - Ny)
    iz_c = np.where(iz_nz <= Nz // 2, iz_nz, iz_nz - Nz)

    rvec_kept = np.stack([ix_c, iy_c, iz_c], axis=1).astype(np.int64)
    data_out  = data_R[:, :, :, ix_nz, iy_nz, iz_nz]   # (Norb, Norb, N_cut, Nr_nonzero)
    return data_out, rvec_kept

def _write_wannier_dat(fname: str, data_out: np.ndarray, rvec_kept: np.ndarray, iw_grid: np.ndarray,
                       Nw_orig: int, label: str, mu: float | None = None, temp: float | None = None) -> None:
    """
    @fn _write_wannier_dat
    @brief Write Wannier-R Matsubara data as text in extended ham_r.dat format.
    @param    fname: output file path (without .dat)
    @param data_out: (Norb, Norb, N_cut, Nr) complex128
    @param rvec_kept: (Nr, 3) int64
    @param  iw_grid: (N_cut,) Matsubara frequencies in eV
    Format per line: Rx Ry Rz  io jo  iw   Re   Im
    """
    Norb, _, N_cut, Nr = data_out.shape
    hdr = f'# {label}: Norb={Norb} N_iw={N_cut} Nr={Nr}'
    if mu   is not None: hdr += f' mu={mu:.8f}'
    if temp is not None: hdr += f' temp={temp:.8e}'

    # Build flat index arrays — order: (ir, n, i, j) C-contiguous
    ir_v, n_v, i_v, j_v = np.meshgrid(
        np.arange(Nr), np.arange(N_cut), np.arange(Norb), np.arange(Norb), indexing='ij')
    vals = data_out[i_v, j_v, n_v, ir_v]   # (Nr, N_cut, Norb, Norb)

    Rx = rvec_kept[ir_v.ravel(), 0].astype(np.int64)
    Ry = rvec_kept[ir_v.ravel(), 1].astype(np.int64)
    Rz = rvec_kept[ir_v.ravel(), 2].astype(np.int64)
    io = (i_v.ravel() + 1).astype(np.int64)
    jo = (j_v.ravel() + 1).astype(np.int64)
    iw = (n_v.ravel() + 1).astype(np.int64)
    Re = vals.ravel().real
    Im = vals.ravel().imag

    try:
        lines = (np.char.mod('%4d', Rx) + np.char.mod('%4d', Ry) + np.char.mod('%4d', Rz)
                 + np.char.mod('%3d', io) + np.char.mod('%3d', jo) + np.char.mod('%5d', iw)
                 + np.char.mod('%16.8e', Re) + np.char.mod('%16.8e', Im))
        with open(f'{fname}.dat', 'w') as f:
            f.write(hdr + f'\n{Norb} {N_cut} {Nr}\n')
            f.write('\n'.join(lines.tolist()) + '\n')
    except IOError as e:
        print(f'Error: Failed to write {fname}.dat: {e}', flush=True)
        return

    print(f'{label}: Matsubara {Nw_orig} -> {N_cut} | Nr={Nr} | '
          f'{data_out.nbytes / 1e6:.1f} MB', flush=True)

def output_self_wannier(sigmak: np.ndarray, mu_self: float, kmap: np.ndarray, invk: np.ndarray,
                        Nx: int, Ny: int, Nz: int, Nw: int, temp: float, N_cut: int = 64,
                        zero_tol: float = 1e-5, fname: str = 'self_en_wannier') -> None:
    """
    @fn output_self_wannier
    @brief Convert self-energy Sigma(k, iw_n) -> Sigma(R, iw_n) in Wannier-real-space format.
    R-vectors are all grid points with at least one non-zero (orbital, Matsubara) component
    after zero_tol zeroing. Writes fname.npz (binary) and fname.dat (text).
    @param   sigmak: [Norb, Norb, Nw, Nk_irr] complex128 — output of mkself/mkself_soc
    @param  mu_self: chemical potential with self-energy correction (eV)
    @param     kmap: [Nkall, 3] int64 — FFT grid indices from gen_irr_k_TRS
    @param     invk: [Nkall, 3] int64 — irr-k mapping with TRS flag from gen_irr_k_TRS
    @param    N_cut: number of leading Matsubara frequencies to keep (default 64)
    @param  zero_tol: Re/Im components below zero_tol * max|Sigma| are set to 0 (default 1e-5)
    @param    fname: output file base name
    """
    sigma_full = _irr_to_full_kgrid(sigmak, invk, kmap, Nx, Ny, Nz)
    s11_iw0 = sigma_full[0, 0, 0]
    print(f'[Pre-FFT] Sigma_11(iw_0): Gamma={sigma_full[0,0,0,0,0,0]:.6e}, '
          f'max|..|={np.abs(s11_iw0).max():.6e}', flush=True)
    data_out, rvec_kept = _wannier_all_nonzero(sigma_full, Nx, Ny, Nz, N_cut, zero_tol)
    if data_out.size == 0:
        print('output_self_wannier: no non-zero R-vectors found', flush=True)
        return
    N_cut_used = data_out.shape[2]
    iw_grid    = (2 * np.arange(N_cut_used) + 1) * np.pi * temp
    try:
        np.savez(fname, sigma=data_out, rvec=rvec_kept, iw=iw_grid, mu=mu_self, temp=temp)
    except IOError as e:
        print(f'Error: Failed to write {fname}.npz: {e}', flush=True)
    _write_wannier_dat(fname, data_out, rvec_kept, iw_grid, Nw,
                       'Self-energy', mu=mu_self, temp=temp)

def output_gap_wannier(gap: np.ndarray, kmap: np.ndarray, invk: np.ndarray, Nx: int, Ny: int, Nz: int,
                       Nw: int, temp: float, N_cut: int = 64, zero_tol: float = 1e-5, fname: str = 'gap_wannier') -> None:
    """
    @fn output_gap_wannier
    @brief Convert gap function Delta(k, iw_n) -> Delta(R, iw_n) in Wannier-real-space format.
    R-vectors are all grid points with at least one non-zero (orbital, Matsubara) component
    after zero_tol zeroing. Writes fname.npz (binary) and fname.dat (text).
    Without SOC: gap is on irreducible k-points [Norb, Norb, Nw, Nk_irr]; expanded via invk+TRS.
    With SOC:    gap is on full BZ              [Norb, Norb, Nw, Nkall];   mapped directly via kmap.
    @param      gap: [Norb, Norb, Nw, Nk_irr or Nkall] complex128
    @param     kmap: [Nkall, 3] int64 — FFT grid indices from gen_irr_k_TRS
    @param     invk: [Nkall, 3] int64 — irr-k mapping with TRS flag from gen_irr_k_TRS
    @param    N_cut: number of leading Matsubara frequencies to keep (default 64)
    @param  zero_tol: Re/Im components below zero_tol * max|Delta| are set to 0 (default 1e-5)
    @param    fname: output file base name
    """
    Nkall = Nx * Ny * Nz
    if gap.shape[-1] == Nkall:   # SOC: already on full BZ
        Norb = gap.shape[0]
        ix_arr, iy_arr, iz_arr = kmap[:, 0], kmap[:, 1], kmap[:, 2]
        gap_full = np.zeros((Norb, Norb, Nw, Nx, Ny, Nz), dtype=np.complex128)
        gap_full[:, :, :, ix_arr, iy_arr, iz_arr] = gap
    else:                        # non-SOC: on irreducible k-points, expand with TRS
        gap_full = _irr_to_full_kgrid(gap, invk, kmap, Nx, Ny, Nz)
    g11_iw0 = gap_full[0, 0, 0]
    print(f'[Pre-FFT] Delta_11(iw_0): Gamma={gap_full[0,0,0,0,0,0]:.6e}, '
          f'max|..|={np.abs(g11_iw0).max():.6e}', flush=True)
    data_out, rvec_kept = _wannier_all_nonzero(gap_full, Nx, Ny, Nz, N_cut, zero_tol)
    if data_out.size == 0:
        print('output_gap_wannier: no non-zero R-vectors found', flush=True)
        return
    N_cut_used = data_out.shape[2]
    iw_grid    = (2 * np.arange(N_cut_used) + 1) * np.pi * temp
    try:
        np.savez(fname, gap=data_out, rvec=rvec_kept, iw=iw_grid, temp=temp)
    except IOError as e:
        print(f'Error: Failed to write {fname}.npz: {e}', flush=True)
    _write_wannier_dat(fname, data_out, rvec_kept, iw_grid, Nw,
                       'Gap function', temp=temp)
