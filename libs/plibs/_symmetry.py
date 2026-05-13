#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Point group symmetry operations, symmetry detection, and irreducible k-point generation.
"""
import numpy as np


def _standard_point_group_ops() -> list[np.ndarray]:
    """Return standard point group operations as integer 3×3 matrices (covers Oh + D4h)."""
    I = np.eye(3, dtype=int)
    ops = [I.copy()]
    ops.append(-I)                                                           # inversion
    for ax in range(3):                                                      # C2 about x, y, z
        R = -I.copy(); R[ax, ax] = 1; ops.append(R)
    ops.append(np.array([[ 0,-1, 0],[ 1, 0, 0],[ 0, 0, 1]]))               # C4z
    ops.append(np.array([[ 0, 1, 0],[-1, 0, 0],[ 0, 0, 1]]))               # C4z^{-1}
    ops.append(np.array([[ 1, 0, 0],[ 0, 0,-1],[ 0, 1, 0]]))               # C4x
    ops.append(np.array([[ 1, 0, 0],[ 0, 0, 1],[ 0,-1, 0]]))               # C4x^{-1}
    ops.append(np.array([[ 0, 0, 1],[ 0, 1, 0],[-1, 0, 0]]))               # C4y
    ops.append(np.array([[ 0, 0,-1],[ 0, 1, 0],[ 1, 0, 0]]))               # C4y^{-1}
    ops.append(np.array([[ 0,-1, 0],[ 1, 0, 0],[ 0, 0,-1]]))               # S4z
    ops.append(np.array([[ 0, 1, 0],[-1, 0, 0],[ 0, 0,-1]]))               # S4z^{-1}
    for ax in range(3):                                                      # mirror x, y, z
        R = I.copy(); R[ax, ax] = -1; ops.append(R)
    # C3 about all four body diagonals [±1,±1,±1] and their inverses
    ops.append(np.array([[ 0, 0, 1],[ 1, 0, 0],[ 0, 1, 0]]))               # C3 [111]
    ops.append(np.array([[ 0, 1, 0],[ 0, 0, 1],[ 1, 0, 0]]))               # C3^2 [111]
    ops.append(np.array([[ 0, 0,-1],[-1, 0, 0],[ 0, 1, 0]]))               # C3 [1-1-1]
    ops.append(np.array([[ 0,-1, 0],[ 0, 0, 1],[-1, 0, 0]]))               # C3^2 [1-1-1]
    ops.append(np.array([[ 0, 0, 1],[-1, 0, 0],[ 0,-1, 0]]))               # C3 [-11-1]
    ops.append(np.array([[ 0,-1, 0],[ 0, 0,-1],[ 1, 0, 0]]))               # C3^2 [-11-1]
    ops.append(np.array([[ 0, 0,-1],[ 1, 0, 0],[ 0,-1, 0]]))               # C3 [-1-11]
    ops.append(np.array([[ 0, 1, 0],[ 0, 0,-1],[-1, 0, 0]]))               # C3^2 [-1-11]
    # Face-diagonal mirrors σ_d and C2 rotations about face diagonals
    ops.append(np.array([[ 0, 1, 0],[ 1, 0, 0],[ 0, 0, 1]]))               # σ_d [110]
    ops.append(np.array([[ 0,-1, 0],[-1, 0, 0],[ 0, 0, 1]]))               # σ_d [1-10]
    ops.append(np.array([[ 0, 1, 0],[ 1, 0, 0],[ 0, 0,-1]]))               # C2 [110]
    ops.append(np.array([[ 0,-1, 0],[-1, 0, 0],[ 0, 0,-1]]))               # C2 [1-10]
    ops.append(np.array([[ 0, 0, 1],[ 0, 1, 0],[ 1, 0, 0]]))               # σ [101]
    ops.append(np.array([[ 0, 0,-1],[ 0, 1, 0],[-1, 0, 0]]))               # σ [10-1]
    ops.append(np.array([[ 0, 0, 1],[ 0,-1, 0],[ 1, 0, 0]]))               # C2 [101]
    ops.append(np.array([[ 0, 0,-1],[ 0,-1, 0],[-1, 0, 0]]))               # C2 [10-1]
    ops.append(np.array([[ 1, 0, 0],[ 0, 0, 1],[ 0, 1, 0]]))               # σ [011]
    ops.append(np.array([[ 1, 0, 0],[ 0, 0,-1],[ 0,-1, 0]]))               # σ [01-1]
    ops.append(np.array([[-1, 0, 0],[ 0, 0, 1],[ 0, 1, 0]]))               # C2 [011]
    ops.append(np.array([[-1, 0, 0],[ 0, 0,-1],[ 0,-1, 0]]))               # C2 [01-1]
    return ops


def _find_U_for_S(rvec_int: np.ndarray, ham_r: np.ndarray,
                  rvec_dict: dict, Si: np.ndarray,
                  tol: float, min_fraction: float) -> np.ndarray | None:
    """
    Numerically find the unitary U satisfying H(S·R) = U H(R) U† for available R-pairs.

    Strategy (tried in order):
      1. Full signed-permutation search (Norb ≤ 6):
         Enumerate all Norb! × 2^Norb monomial matrices.  Quick pre-filter on one pair.
      2. Block-structured signed-permutation search (Norb > 6):
         Factor Norb = n_sites × n_orb and assume U = P_site ⊗ U_orb (same orbital
         transformation for all equivalent sites).  Covers most physical multi-site models.
         Tries all (n_sites! site permutations) × (n_orb! × 2^n_orb orbital matrices).
      3. Null-space + rounding to nearest monomial (general fallback):
         Stack constraints H(SR)U = UH(R) as a linear system, find the null vector via SVD,
         project to unitary (polar decomposition), round to nearest monomial, and verify.

    Returns U (Norb×Norb complex128) if valid, None otherwise.
    """
    import math
    from itertools import permutations

    Norb = ham_r.shape[1]
    pairs: list[tuple[int, int]] = []
    for i, R in enumerate(rvec_int):
        SR = tuple((Si @ R).tolist())
        if SR in rvec_dict:
            pairs.append((i, rvec_dict[SR]))

    if len(pairs) < min_fraction * len(rvec_int):
        return None

    def _verify(U: np.ndarray) -> bool:
        return all(np.abs(ham_r[j] - U @ ham_r[i] @ U.conj().T).max() < tol
                   for i, j in pairs)

    i0, j0 = pairs[0]
    A0, B0 = ham_r[i0], ham_r[j0]

    def _try_monomial_block(n_orb: int, orb_perm: tuple, signs_int: int,
                             n_sites: int, site_perm: tuple) -> np.ndarray | None:
        """Build and quick-check one block-monomial candidate U."""
        U_orb = np.zeros((n_orb, n_orb), dtype=complex)
        for row in range(n_orb):
            U_orb[row, orb_perm[row]] = 1 if (signs_int >> row) & 1 else -1
        U_cand = np.zeros((Norb, Norb), dtype=complex)
        for s_out, s_in in enumerate(site_perm):
            sl_out = slice(s_out * n_orb, (s_out + 1) * n_orb)
            sl_in  = slice(s_in  * n_orb, (s_in  + 1) * n_orb)
            U_cand[sl_out, sl_in] = U_orb
        # Quick pre-filter on first pair
        if np.abs(B0 @ U_cand - U_cand @ A0).max() > tol * 20:
            return None
        return U_cand if _verify(U_cand) else None

    # --- Strategy 1: full signed-permutation search (small Norb) ---
    MAX_MONOMIAL = 100_000
    if math.factorial(Norb) * (2 ** Norb) <= MAX_MONOMIAL:
        for perm in permutations(range(Norb)):
            for signs_int in range(2 ** Norb):
                U = _try_monomial_block(Norb, perm, signs_int, 1, (0,))
                if U is not None:
                    return U
    else:
        # --- Strategy 2: null-space guided block-monomial search ---
        # Step 2a: null-space SVD on ALL available pairs (used as primary reject filter
        # AND to guide the block search).  A clear null space (s_min/s_max < tol_null)
        # is a necessary condition; without it we return None immediately.
        from scipy.optimize import linear_sum_assignment
        rows_ns = [np.kron(np.eye(Norb), ham_r[j]) - np.kron(ham_r[i].T, np.eye(Norb))
                   for i, j in pairs]
        M_ns = np.vstack(rows_ns)
        _, s_ns, Vh_ns = np.linalg.svd(M_ns, full_matrices=False)
        tol_null = max(tol * 0.5, 5e-3)              # tighten relative to verification tol
        if s_ns[-1] / max(s_ns[0], 1e-15) > tol_null:   # no clear null space → not a symmetry
            return None
        U_guide = Vh_ns[-1].reshape(Norb, Norb, order='F')   # raw (not unitary)

        # Step 2b: for each factorization Norb = n_sites × n_orb, use U_guide to infer
        # the block permutation structure, then enumerate only nearby monomial candidates.
        for n_orb in range(2, Norb):
            if Norb % n_orb != 0:
                continue
            if math.factorial(n_orb) * (2 ** n_orb) > MAX_MONOMIAL:
                continue
            n_sites = Norb // n_orb
            # Determine the most likely site permutation from U_guide block norms
            block_norm = np.array(
                [[np.linalg.norm(U_guide[r*n_orb:(r+1)*n_orb, c*n_orb:(c+1)*n_orb])
                  for c in range(n_sites)] for r in range(n_sites)]
            )
            _, site_col = linear_sum_assignment(-block_norm)
            site_perm = tuple(site_col.tolist())
            # Infer orbital permutation + signs from the dominant block
            dominant_block = U_guide[0:n_orb, site_perm[0]*n_orb:(site_perm[0]+1)*n_orb]
            _, orb_col = linear_sum_assignment(-np.abs(dominant_block))
            # Try the inferred candidate plus small perturbations (±signs)
            for signs_int in range(2 ** n_orb):
                U = _try_monomial_block(n_orb, tuple(orb_col.tolist()), signs_int,
                                        n_sites, site_perm)
                if U is not None:
                    return U
            # Fallback: enumerate all orbital permutations for this site permutation
            for orb_perm in permutations(range(n_orb)):
                for signs_int in range(2 ** n_orb):
                    U = _try_monomial_block(n_orb, orb_perm, signs_int, n_sites, site_perm)
                    if U is not None:
                        return U

    # --- Strategy 3: null-space + nearest-monomial (general fallback) ---
    from scipy.optimize import linear_sum_assignment
    subset = pairs[::max(1, len(pairs) // 200)][:200]
    rows = [np.kron(np.eye(Norb), ham_r[j]) - np.kron(ham_r[i].T, np.eye(Norb))
            for i, j in subset]
    M = np.vstack(rows)
    _, s, Vh = np.linalg.svd(M, full_matrices=False)
    if s[-1] / max(s[0], 1e-15) > 0.1:
        return None
    U_raw = Vh[-1].reshape(Norb, Norb, order='F')
    V, _, Wh = np.linalg.svd(U_raw)
    U_smooth = V @ Wh                       # nearest general unitary
    # Round to nearest monomial
    row_ind, col_ind = linear_sum_assignment(-np.abs(U_smooth))
    U_mono = np.zeros((Norb, Norb), dtype=complex)
    for r, c in zip(row_ind, col_ind):
        if abs(U_smooth[r, c]) < 0.3:
            break
        U_mono[r, c] = np.sign(U_smooth[r, c].real) if abs(U_smooth[r, c].imag) < 0.3 else 0
    else:
        if _verify(U_mono):
            return U_mono
    return U_smooth if _verify(U_smooth) else None


def detect_symm_ops_auto_U(rvec: np.ndarray, ham_r: np.ndarray,
                            candidate_S: list[np.ndarray] | None = None,
                            tol: float = 1e-3,
                            min_fraction: float = 0.9
                            ) -> list[tuple[np.ndarray, np.ndarray]]:
    """
    @fn detect_symm_ops_auto_U
    @brief Detect point group symmetry operations and automatically determine the orbital
    transformation matrix U for each, including multi-orbital systems with orbital mixing
    (e.g. E_g representations with dxz/dyz mixing under C4 rotation).
    For each candidate S, numerically solves H(S·R) = U H(R) U† for U using eigendecomposition
    of a random-weighted probe matrix, then verifies on all available R-pairs.
    Works with truncated R-vector sets (first-principles models) via min_fraction threshold.
    @param      rvec: R-vectors [Nr, 3] float64
    @param     ham_r: Hopping matrix [Nr, Norb, Norb] complex128
    @param candidate_S: Candidate 3×3 integer rotation matrices; uses standard ops if None
    @param       tol: Verification tolerance for H(SR) ≈ U H(R) U†
    @param min_fraction: Minimum fraction of R-vectors that must have their S-image present
    @retval  sym_ops: List of valid (S, U) pairs; identity is always first
    """
    rv_int = np.round(rvec).astype(int)
    rvec_dict = {(int(r[0]), int(r[1]), int(r[2])): i for i, r in enumerate(rv_int)}

    if candidate_S is None:
        candidate_S = _standard_point_group_ops()

    valid: list[tuple[np.ndarray, np.ndarray]] = []
    for S in candidate_S:
        Si = np.round(S).astype(int)
        U = _find_U_for_S(rv_int, ham_r, rvec_dict, Si, tol, min_fraction)
        if U is not None:
            valid.append((S, U))
    return valid


def detect_symm_ops(rvec: np.ndarray, ham_r: np.ndarray,
                    candidate_S: list[np.ndarray] | None = None,
                    U_list: list[np.ndarray] | None = None,
                    tol: float = 1e-4) -> list[tuple[np.ndarray, np.ndarray]]:
    """
    @fn detect_symm_ops
    @brief Detect point group symmetry operations that leave the Hamiltonian invariant.
    For each candidate (S, U) pair, checks H(S·R) ≈ U @ H(R) @ U† for all R in rvec.
    @param      rvec: R-vectors [Nr, 3] float64
    @param     ham_r: Hopping matrix [Nr, Norb, Norb] complex128
    @param candidate_S: Candidate 3×3 integer rotation matrices; uses standard ops if None
    @param    U_list: Norb×Norb orbital transformation matrices paired with candidate_S;
                      defaults to identity (no orbital mixing) if None
    @param       tol: Element-wise tolerance for comparison
    @retval  sym_ops: List of valid (S, U) pairs; identity is always first
    """
    rv_int = np.round(rvec).astype(int)
    rvec_dict = {(int(r[0]), int(r[1]), int(r[2])): i for i, r in enumerate(rv_int)}

    if candidate_S is None:
        candidate_S = _standard_point_group_ops()

    Norb = ham_r.shape[1]
    if U_list is None:
        U_list = [np.eye(Norb, dtype=complex)] * len(candidate_S)

    valid: list[tuple[np.ndarray, np.ndarray]] = []
    for S, U in zip(candidate_S, U_list):
        Si = np.round(S).astype(int)
        ok = True
        for i, R in enumerate(rv_int):
            SR = tuple((Si @ R).tolist())
            if SR not in rvec_dict:
                ok = False
                break
            j = rvec_dict[SR]
            if not np.allclose(ham_r[j], U @ ham_r[i] @ U.conj().T, atol=tol):
                ok = False
                break
        if ok:
            valid.append((S, U))
    return valid


def gen_irr_k_symm(Nx: int, Ny: int, Nz: int,
                   sym_ops: list[tuple[np.ndarray, np.ndarray]]
                   ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    @fn gen_irr_k_symm
    @brief Generate the irreducible k-point list using a set of point group symmetry operations.
    Each S in sym_ops acts on integer k-indices [ix, iy, iz] as Sk = S @ [ix, iy, iz] mod [Nx, Ny, Nz].
    Operations that do not map the mesh onto itself are silently skipped.
    @param    Nx: Number of k-points along kx
    @param    Ny: Number of k-points along ky
    @param    Nz: Number of k-points along kz
    @param sym_ops: List of (S, U) tuples from detect_symm_ops (identity must be first)
    @retval klist_irr: Irreducible k-coordinates [Nk_irr, 3] float64 in [0, 1)
    @retval    kmap: Full k-grid FFT indices [Nkall, 3] int64
    @retval    invk: Mapping [Nkall, 2] int64 — col0: irr-k index (0-based), col1: sym op index
    """
    N = np.array([Nx, Ny, Nz])
    Nkall = Nx * Ny * Nz

    ixs = np.arange(Nkall, dtype=np.int64) % Nx
    iys = (np.arange(Nkall, dtype=np.int64) // Nx) % Ny
    izs = np.arange(Nkall, dtype=np.int64) // (Nx * Ny)
    kmap = np.stack([ixs, iys, izs], axis=1)

    def to_flat(k_int: np.ndarray) -> int:
        km = k_int % N
        return int(km[0] + Nx * (km[1] + Ny * km[2]))

    irr_map  = np.full(Nkall, -1, dtype=np.int64)
    iop_map  = np.zeros(Nkall, dtype=np.int64)
    irr_flat: list[int] = []

    for flat_idx in range(Nkall):
        if irr_map[flat_idx] >= 0:
            continue
        irr_idx = len(irr_flat)
        irr_flat.append(flat_idx)
        irr_map[flat_idx] = irr_idx
        iop_map[flat_idx] = 0

        k_int = kmap[flat_idx]
        for iop, (S, _) in enumerate(sym_ops):
            if iop == 0:
                continue
            Sk = np.round(S @ k_int).astype(int)
            fi = to_flat(Sk)
            if irr_map[fi] < 0:
                irr_map[fi]  = irr_idx
                iop_map[fi]  = iop

    klist_irr = kmap[irr_flat].astype(np.float64) / N
    invk = np.stack([irr_map, iop_map], axis=1)
    return klist_irr, kmap, invk


def expand_irr_to_full_symm(data_irr: np.ndarray, invk: np.ndarray,
                              kmap: np.ndarray,
                              sym_ops: list[tuple[np.ndarray, np.ndarray]],
                              Nx: int, Ny: int, Nz: int,
                              orb_axes: tuple[int, int] | None = None) -> np.ndarray:
    """
    @fn expand_irr_to_full_symm
    @brief Expand data from irreducible k-points to the full BZ using symmetry operations.
    @param data_irr: Array whose last axis is Nk_irr.  Shape: [..., Nk_irr]
    @param     invk: [Nkall, 2] int64 — col0: irr-k index, col1: sym op index (from gen_irr_k_symm)
    @param     kmap: [Nkall, 3] int64 — FFT grid indices
    @param  sym_ops: List of (S, U) from detect_symm_ops
    @param  orb_axes: (ax1, ax2) in data shape that are orbital indices for U-rotation.
                      Pass None for scalar or already-diagonal data (e.g. eigenvalues).
    @retval    full: Array with shape [..., Nx, Ny, Nz]
    """
    full = np.zeros(data_irr.shape[:-1] + (Nx, Ny, Nz), dtype=data_irr.dtype)
    ik_arr  = invk[:, 0]
    iop_arr = invk[:, 1]
    ix_arr, iy_arr, iz_arr = kmap[:, 0], kmap[:, 1], kmap[:, 2]

    for flat_idx in range(Nx * Ny * Nz):
        d   = data_irr[..., ik_arr[flat_idx]]
        iop = int(iop_arr[flat_idx])
        if orb_axes is not None and iop > 0:
            _, U = sym_ops[iop]
            ax1, ax2 = orb_axes
            # Apply d → U d U† along orb_axes: contract U on ax1 (left), U* on ax2 (right)
            d = np.tensordot(U, d, axes=[[1], [ax1]])
            d = np.moveaxis(d, 0, ax1)
            d = np.tensordot(d, U.conj(), axes=[[ax2], [1]])
            d = np.moveaxis(d, -1, ax2)
        full[..., ix_arr[flat_idx], iy_arr[flat_idx], iz_arr[flat_idx]] = d

    return full
