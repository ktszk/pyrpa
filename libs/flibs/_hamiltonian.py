import ctypes
import numpy as np
import libs.flibs._loader as _loader

# Public API of this submodule; keeps `from ._x import *` in __init__.py from
# re-exporting numpy, scipy, ctypes and the sibling packages imported above.
__all__ = [
    'omp_params', 'gen_ham', 'get_eig', 'get_ffermi', 'get_vlm0', 'get_vk', 'get_vnm', 'get_vnmk', 'get_veloc',
    'get_imass0', 'get_imassk', 'get_mass', 'get_qshift', 'get_iqshift', 'gen_irr_k_TRS', 'gen_kpoint_weight',
    'mkBdGhamk', 'get_plist'
]

# --- ctypes signatures: set once at import.
# All Fortran entry points are subroutines, so restype is always None.
_loader._lib.openmp_params.argtypes = [
    ctypes.POINTER(ctypes.c_int64),                 # omp_num (out)
    ctypes.POINTER(ctypes.c_bool),                  # omp_check (out)
]
_loader._lib.openmp_params.restype = None
_loader._lib.gen_ham.argtypes = [
    _loader.c_array(np.complex128),          # ham_k (out)
    _loader.c_array(np.float64),             # klist
    _loader.c_array(np.complex128),          # ham_r
    _loader.c_array(np.float64),             # rvec
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)  # Nk, Nr, Norb
]
_loader._lib.gen_ham.restype = None
_loader._lib.get_eig.argtypes = [
    _loader.c_array(np.float64),     # eig (out)
    _loader.c_array(np.complex128),  # uni (out)
    _loader.c_array(np.complex128),  # ham_k
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)            # Nk, Norb
]
_loader._lib.get_eig.restype = None
_loader._lib.get_eig_mlo.argtypes = [
    _loader.c_array(np.float64),     # eig (out)
    _loader.c_array(np.complex128),  # uni (out)
    _loader.c_array(np.complex128),  # ham_k
    _loader.c_array(np.complex128),  # Ovlk
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)            # Nk, Norb
]
_loader._lib.get_eig_mlo.restype = None
_loader._lib.get_ffermi.argtypes = [
    _loader.c_array(np.float64),  # ffermi (out)
    _loader.c_array(np.float64),  # eig
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),      # mu, temp
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)         # Nk, Norb
]
_loader._lib.get_ffermi.restype = None
_loader._lib.get_vlm0.argtypes = [
    _loader.c_array(np.complex128),          # vk (out)
    _loader.c_array(np.float64),             # klist
    _loader.c_array(np.complex128),          # ham_r
    _loader.c_array(np.float64),             # rvec
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)  # Nk, Nr, Norb
]
_loader._lib.get_vlm0.restype = None
_loader._lib.get_veloc.argtypes = [
    _loader.c_array(np.float64),     # vk (out)
    _loader.c_array(np.complex128),  # vk0
    _loader.f_array(np.float64),     # mrot (Fortran order: v = mrot @ dE/dk)
    _loader.c_array(np.complex128),  # uni
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)            # Nk, Norb
]
_loader._lib.get_veloc.restype = None
_loader._lib.get_vnm.argtypes = [
    _loader.c_array(np.complex128),  # vk (out)
    _loader.c_array(np.complex128),  # vk0
    _loader.f_array(np.float64),     # mrot (Fortran order: v = mrot @ dE/dk)
    _loader.c_array(np.complex128),  # uni
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)            # Nk, Norb
]
_loader._lib.get_vnm.restype = None
_loader._lib.get_veloc_mlo.argtypes = [
    _loader.c_array(np.float64),     # vk (out)
    _loader.c_array(np.complex128),  # vk0
    _loader.c_array(np.complex128),  # sk0
    _loader.c_array(np.float64),     # eig
    _loader.f_array(np.float64),     # mrot (Fortran order: v = mrot @ dE/dk)
    _loader.c_array(np.complex128),  # uni
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)            # Nk, Norb
]
_loader._lib.get_veloc_mlo.restype = None
_loader._lib.get_vnm_mlo.argtypes = [
    _loader.c_array(np.complex128),  # vk (out)
    _loader.c_array(np.complex128),  # vk0 = dH/dk
    _loader.c_array(np.complex128),  # sk0 = dS/dk
    _loader.c_array(np.complex128),  # Ovlk = S(k)
    _loader.c_array(np.float64),     # eig
    _loader.f_array(np.float64),     # mrot (Fortran order: v = mrot @ dE/dk)
    _loader.c_array(np.complex128),  # uni
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)            # Nk, Norb
]
_loader._lib.get_vnm_mlo.restype = None
_loader._lib.get_imass0.argtypes = [
    _loader.c_array(np.complex128),          # imk (out)
    _loader.c_array(np.float64),             # klist
    _loader.c_array(np.complex128),          # ham_r
    _loader.c_array(np.float64),             # rvec
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)  # Nk, Nr, Norb
]
_loader._lib.get_imass0.restype = None
_loader._lib.get_imassk.argtypes = [
    _loader.c_array(np.float64),     # imk (out)
    _loader.c_array(np.complex128),  # imk0
    _loader.c_array(np.complex128),  # vk0
    _loader.c_array(np.float64),     # eig
    _loader.f_array(np.float64),     # mrot (Fortran order: v = mrot @ dE/dk)
    _loader.c_array(np.complex128),  # uni
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)            # Nk, Norb
]
_loader._lib.get_imassk.restype = None
_loader._lib.get_imassk_mlo.argtypes = [
    _loader.c_array(np.float64),     # imk (out)
    _loader.c_array(np.complex128),  # imk0
    _loader.c_array(np.complex128),  # smk0
    _loader.c_array(np.complex128),  # vk0
    _loader.c_array(np.complex128),  # sk0
    _loader.c_array(np.float64),     # eig
    _loader.f_array(np.float64),     # mrot (Fortran order: v = mrot @ dE/dk)
    _loader.c_array(np.complex128),  # uni
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)            # Nk, Norb
]
_loader._lib.get_imassk_mlo.restype = None
_loader._lib.get_qshift_.argtypes = [
    _loader.c_array(np.float64),  # qpoint
    _loader.c_array(np.float64),  # klist
    _loader.c_array(np.int64),    # qshift (out)
    ctypes.POINTER(ctypes.c_int64),                          # Nk
]
_loader._lib.get_qshift_.restype = None
_loader._lib.get_iqshift_.argtypes = [
    _loader.c_array(np.float64),  # qpoint
    _loader.c_array(np.float64),  # klist
    _loader.c_array(np.int64),    # qshift (out)
    ctypes.POINTER(ctypes.c_int64),                          # Nk
]
_loader._lib.get_iqshift_.restype = None
_loader._lib.generate_irr_kpoint_inv.argtypes = [
    _loader.c_array(np.float64),  # klist (out)
    _loader.c_array(np.int64),    # kmap (out)
    _loader.c_array(np.int64),    # invk_ft_list (out)
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),        # Nk, Nx
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)         # Ny, Nz
]
_loader._lib.generate_irr_kpoint_inv.restype = None
_loader._lib.get_kweight.argtypes = [
    _loader.c_array(np.float64),  # weight (out)
    _loader.c_array(np.int64),    # invk
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)         # Nk, Nkall
]
_loader._lib.get_kweight.restype = None
_loader._lib.mkbdghamk.argtypes = [
    _loader.c_array(np.complex128),  # hamBdGk (out)
    _loader.c_array(np.complex128),  # hamk
    _loader.c_array(np.complex128),  # delta
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nk, Norb
]
_loader._lib.mkbdghamk.restype = None
_loader._lib.get_parity_prop.argtypes = [
    _loader.c_array(np.float64),     # Pmn (inout)
    _loader.c_array(np.float64),     # rvec
    _loader.c_array(np.complex128),  # ham_r
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)            # Norb, Nr
]
_loader._lib.get_parity_prop.restype = None


def omp_params() -> tuple[int, bool]:
    """
    @brief Query the underlying Fortran library for OpenMP configuration.
    @return num_threads: Maximum number of OpenMP threads the Fortran code will use
    @return  valid_flag: True if OpenMP was initialized correctly
    """
    omp_num = ctypes.c_int64()
    omp_check = ctypes.c_bool()
    _loader._lib.openmp_params(ctypes.byref(omp_num), ctypes.byref(omp_check))
    return omp_num.value, omp_check.value


def gen_ham(klist: np.ndarray, ham_r: np.ndarray, rvec: np.ndarray,
            Ovl_r: np.ndarray | None = None) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    """
    @brief Compute the k-space Hamiltonian by Fourier-transforming real-space hopping blocks.
    @param  klist: k-point list in fractional coordinates [Nk, 3]
    @param  ham_r: Real-space Hamiltonian blocks [Nr, Norb, Norb] complex128
    @param   rvec: Real-space displacement vectors [Nr, 3]
    @param  Ovl_r: Optional overlap matrix blocks [Nr, Norb, Norb]; if given, also returns Ovlk
    @return  hamk: k-space Hamiltonian [Nk, Norb, Norb] complex128
    @return  Ovlk: k-space overlap matrix [Nk, Norb, Norb] complex128 (only when Ovl_r is given)
    """
    Nk, Nr = len(klist), len(rvec)
    assert ham_r.ndim == 3, "ham_r must be a 3‑D array"
    Norb = ham_r.shape[1]
    hamk = np.zeros((Nk, Norb, Norb), dtype=np.complex128)
    _loader._lib.gen_ham(hamk, klist, ham_r, rvec,
                         _loader.i64(Nk), _loader.i64(Nr), _loader.i64(Norb))
    if Ovl_r is not None:
        Ovlk = np.zeros((Nk, Norb, Norb), dtype=np.complex128)
        _loader._lib.gen_ham(Ovlk, klist, Ovl_r, rvec,
                             _loader.i64(Nk), _loader.i64(Nr), _loader.i64(Norb))
        return hamk, Ovlk
    else:
        return hamk


def get_eig(hamk: np.ndarray, Ovlk: np.ndarray | None = None,
            sw: bool = True) -> tuple[np.ndarray, np.ndarray] | np.ndarray:
    """
    @brief Diagonalize the k-space Hamiltonian, optionally solving the generalized eigenvalue problem.
    @param  hamk: k-space Hamiltonian [Nk, Norb, Norb] complex128
    @param  Ovlk: Overlap matrices [Nk, Norb, Norb] complex128; if given, solves generalized EVP
    @param    sw: If True, return both eigenvalues and eigenvectors; if False, return only eigenvalues
    @return  eig: Eigenvalues [Nk, Norb] float64
    @return  uni: Eigenvectors [Nk, Norb, Norb] complex128 (only returned when sw=True)
    """
    Nk = len(hamk)
    Norb = hamk.shape[1]
    eig = np.zeros((Nk, Norb), dtype=np.float64)
    uni = np.zeros((Nk, Norb, Norb), dtype=np.complex128)

    if Ovlk is None:
        _loader._lib.get_eig(eig, uni, hamk, _loader.i64(Nk), _loader.i64(Norb))
    else:
        _loader._lib.get_eig_mlo(eig, uni, hamk, Ovlk, _loader.i64(Nk), _loader.i64(Norb))

    if sw:
        return eig, uni
    else:
        return eig


def get_ffermi(eig: np.ndarray, mu: float, temp: float) -> np.ndarray:
    """
    @brief Compute Fermi-Dirac occupation numbers from eigenvalues.
    @param    eig: Eigenvalues [Nk, Norb] float64
    @param     mu: Chemical potential in eV
    @param   temp: Temperature in eV
    @return ffermi: Fermi-Dirac occupation numbers [Nk, Norb] float64
    """
    Nk = len(eig)
    Norb = eig.shape[1]
    ffermi = np.zeros((Nk, Norb), dtype=np.float64)
    _loader._lib.get_ffermi(ffermi, eig, _loader.dbl(mu), _loader.dbl(temp),
                            _loader.i64(Nk), _loader.i64(Norb))
    return ffermi


def get_vlm0(klist: np.ndarray, ham_r: np.ndarray, rvec: np.ndarray) -> np.ndarray:
    """
    @brief Compute unrotated band velocity matrices in orbital basis from real-space Hamiltonian.
    @param  klist: k-point list [Nk, 3] float64
    @param  ham_r: Real-space Hamiltonian blocks [Nr, Norb, Norb] complex128
    @param   rvec: Real-space displacement vectors [Nr, 3] float64
    @return   vk0: Velocity matrix elements [Nk, Norb, Norb, 3] complex128
    """
    Nk, Nr = len(klist), len(rvec)
    Norb = ham_r.shape[1]
    vk = np.zeros((Nk, Norb, Norb, 3), dtype=np.complex128)
    _loader._lib.get_vlm0(vk, klist, ham_r, rvec,
                          _loader.i64(Nk), _loader.i64(Nr), _loader.i64(Norb))
    return vk


def _fortran_mrot(mrot: np.ndarray) -> np.ndarray:
    """The Fortran kernels compute v_m = sum_l mrot(m,l) dE/dk_l, so mrot is handed over in
    FORTRAN order: whatever the memory layout of the caller's array (avec.T*ihbar is a
    transposed, Fortran-ordered view; np.eye(3) is C-ordered), the Python index order is
    the Fortran one and v = mrot @ dE/dk holds as written."""
    return np.asfortranarray(mrot, dtype=np.float64)


def get_vk(vk0: np.ndarray, mrot: np.ndarray, uni: np.ndarray,
           sk0: np.ndarray | None = None, eig: np.ndarray | None = None) -> np.ndarray:
    """
    @brief Transform orbital-basis velocity matrices to band-basis velocity expectation values.
    @param   vk0: Unrotated velocity matrices dH/dk [Nk, Norb, Norb, 3] complex128
    @param  mrot: Lattice-to-Cartesian matrix [3, 3] float64, v = mrot @ dE/dk (e.g. avec.T*ihbar)
    @param   uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param   sk0: Overlap derivative dS/dk [Nk, Norb, Norb, 3] complex128 for a non-orthogonal
                  (MLO) basis; None for an orthogonal basis
    @param   eig: Eigenvalues [Nk, Norb] float64, required together with sk0
    @return   vk: Band velocities [Nk, Norb, 3] float64
    """
    mrot = _fortran_mrot(mrot)
    Nk = len(uni)
    Norb = uni.shape[1]
    vk = np.zeros((Nk, Norb, 3), dtype=np.float64)
    if sk0 is None:
        _loader._lib.get_veloc(vk, vk0, mrot, uni, _loader.i64(Nk), _loader.i64(Norb))
    else:
        if eig is None:
            raise ValueError("get_vk: eig is required when sk0 (non-orthogonal basis) is given")
        _loader._lib.get_veloc_mlo(vk, vk0, sk0, eig, mrot, uni, _loader.i64(Nk), _loader.i64(Norb))
    return vk


def get_vnm(vk0: np.ndarray, mrot: np.ndarray, uni: np.ndarray,
            sk0: np.ndarray | None = None, eig: np.ndarray | None = None,
            Ovlk: np.ndarray | None = None) -> np.ndarray:
    """
    @brief Compute inter-band velocity matrix elements v_{nm}(k) for all band pairs and Cartesian directions.
    @param   vk0: Unrotated velocity matrices dH/dk [Nk, Norb, Norb, 3] complex128
    @param  mrot: Lattice-to-Cartesian matrix [3, 3] float64, v = mrot @ dE/dk (e.g. avec.T*ihbar)
    @param   uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param   sk0: Overlap derivative dS/dk [Nk, Norb, Norb, 3] complex128 for a non-orthogonal
                  (MLO) basis; None for an orthogonal basis
    @param   eig: Eigenvalues [Nk, Norb] float64, required together with sk0
    @param  Ovlk: k-space overlap S(k) [Nk, Norb, Norb] complex128, required together with sk0
                  (the exact Loewdin velocity needs S itself, not only its derivative)
    @return   vnm: Inter-band velocity matrix elements [Nk, Norb, Norb, 3] complex128
    """
    mrot = _fortran_mrot(mrot)
    Nk = len(uni)
    Norb = uni.shape[1]
    vk = np.zeros((Nk, Norb, Norb, 3), dtype=np.complex128)
    if sk0 is None:
        _loader._lib.get_vnm(vk, vk0, mrot, uni, _loader.i64(Nk), _loader.i64(Norb))
    else:
        if eig is None:
            raise ValueError("get_vnm: eig is required when sk0 (non-orthogonal basis) is given")
        if Ovlk is None:
            raise ValueError("get_vnm: Ovlk (S(k)) is required when sk0 (non-orthogonal basis) is given")
        _loader._lib.get_vnm_mlo(vk, vk0, sk0, Ovlk, eig, mrot, uni, _loader.i64(Nk), _loader.i64(Norb))
    return vk


def _get_sk0(klist: np.ndarray, S_r: np.ndarray, rvec: np.ndarray) -> np.ndarray | None:
    """dS/dk in the orbital basis, or None for an orthogonal basis (S_r empty/None).
    S(k) is built by the same Fourier sum as H(k), so get_vlm0 differentiates it too."""
    if S_r is None or len(S_r) == 0:
        return None
    return get_vlm0(klist, S_r, rvec)


def get_vnmk(klist: np.ndarray, ham_r: np.ndarray, rvec: np.ndarray,
             mrot: np.ndarray, uni: np.ndarray, S_r: np.ndarray | None = None,
             eig: np.ndarray | None = None) -> np.ndarray:
    """
    @brief Convenience wrapper: compute inter-band velocity matrix elements by composing get_vlm0 and get_vnm.
    @param  klist: k-point list [Nk, 3] float64
    @param  ham_r: Real-space Hamiltonian blocks [Nr, Norb, Norb] complex128
    @param   rvec: Real-space displacement vectors [Nr, 3] float64
    @param   mrot: Lattice-to-Cartesian matrix [3, 3] float64, v = mrot @ dE/dk (e.g. avec.T*ihbar)
    @param    uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param    S_r: Real-space overlap blocks [Nr, Norb, Norb] for a non-orthogonal (MLO) basis;
                   None or [] for an orthogonal basis
    @param    eig: Eigenvalues [Nk, Norb] float64, required together with S_r
    @return   vnm: Inter-band velocity matrix elements [Nk, Norb, Norb, 3] complex128
    """
    vk0 = get_vlm0(klist, ham_r, rvec)
    sk0 = _get_sk0(klist, S_r, rvec)
    Ovlk = None if sk0 is None else gen_ham(klist, S_r, rvec)
    return get_vnm(vk0, mrot, uni, sk0, eig, Ovlk)


def get_veloc(klist: np.ndarray, ham_r: np.ndarray, rvec: np.ndarray,
              mrot: np.ndarray, uni: np.ndarray, S_r: np.ndarray | None = None,
              eig: np.ndarray | None = None) -> np.ndarray:
    """
    @brief Convenience wrapper: compute band-diagonal velocity expectation values by composing get_vlm0 and get_vk.
    @param  klist: k-point list [Nk, 3] float64
    @param  ham_r: Real-space Hamiltonian blocks [Nr, Norb, Norb] complex128
    @param   rvec: Real-space displacement vectors [Nr, 3] float64
    @param   mrot: Lattice-to-Cartesian matrix [3, 3] float64, v = mrot @ dE/dk (e.g. avec.T*ihbar)
    @param    uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param    S_r: Real-space overlap blocks [Nr, Norb, Norb] for a non-orthogonal (MLO) basis;
                   None or [] for an orthogonal basis
    @param    eig: Eigenvalues [Nk, Norb] float64, required together with S_r
    @return    vk: Band velocities [Nk, Norb, 3] float64
    """
    vk0 = get_vlm0(klist, ham_r, rvec)
    return get_vk(vk0, mrot, uni, _get_sk0(klist, S_r, rvec), eig)


def get_imass0(klist: np.ndarray, ham_r: np.ndarray, rvec: np.ndarray) -> np.ndarray:
    """
    @brief Compute the unrotated inverse effective mass tensor in orbital basis.
    @param   klist: k-point list [Nk, 3] float64
    @param   ham_r: Real-space Hamiltonian blocks [Nr, Norb, Norb] complex128
    @param    rvec: Real-space displacement vectors [Nr, 3] float64
    @return imass0: Inverse mass tensor before band rotation [Nk, Norb, Norb, 3, 3] complex128
    """
    Nk, Nr = len(klist), len(rvec)
    Norb = ham_r.shape[1]
    imass0 = np.zeros((Nk, Norb, Norb, 3, 3), dtype=np.complex128)
    _loader._lib.get_imass0(imass0, klist, ham_r, rvec,
                            _loader.i64(Nk), _loader.i64(Nr), _loader.i64(Norb))
    return imass0


def get_imassk(imass0: np.ndarray, mrot: np.ndarray, uni: np.ndarray,
               vk0: np.ndarray, eig: np.ndarray,
               smass0: np.ndarray | None = None, sk0: np.ndarray | None = None) -> np.ndarray:
    """
    @brief Rotate the inverse effective mass tensor from orbital basis into band basis,
           including the interband (2nd-order perturbation) term
           2 sum_{m!=n} Re[H~^a_nm H~^b_mn]/(eps_n-eps_m), which needs dH/dk and eig.
           Band pairs degenerate to within 1e-6 eV are skipped (the mass diverges there).
    @param  imass0: Unrotated d2H/dk2 [Nk, Norb, Norb, 3, 3] complex128 (output of get_imass0)
    @param    mrot: Lattice-to-Cartesian matrix [3, 3] float64, v = mrot @ dE/dk (e.g. avec.T*ihbar)
    @param     uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param     vk0: Orbital-basis dH/dk [Nk, Norb, Norb, 3] complex128 (output of get_vlm0)
    @param     eig: Eigenvalues [Nk, Norb] float64
    @param  smass0: d2S/dk2 [Nk, Norb, Norb, 3, 3] for a non-orthogonal (MLO) basis; None otherwise
    @param     sk0: dS/dk [Nk, Norb, Norb, 3] for a non-orthogonal (MLO) basis; None otherwise
    @return  imass: Band-basis inverse effective mass tensor [Nk, Norb, 3, 3] float64
    """
    mrot = _fortran_mrot(mrot)
    Nk = len(uni)
    Norb = uni.shape[1]
    imass = np.zeros((Nk, Norb, 3, 3), dtype=np.float64)
    if smass0 is None and sk0 is None:
        _loader._lib.get_imassk(imass, imass0, vk0, eig, mrot, uni, _loader.i64(Nk), _loader.i64(Norb))
    else:
        if smass0 is None or sk0 is None:
            raise ValueError("get_imassk: smass0 and sk0 must be given together (non-orthogonal basis)")
        _loader._lib.get_imassk_mlo(imass, imass0, smass0, vk0, sk0, eig, mrot, uni, _loader.i64(Nk), _loader.i64(Norb))
    return imass


def get_mass(klist: np.ndarray, ham_r: np.ndarray, rvec: np.ndarray, mrot: np.ndarray,
             uni: np.ndarray, sw_imass: bool = False, S_r: np.ndarray | None = None,
             eig: np.ndarray | None = None) -> np.ndarray:
    """
    @brief Compute the effective mass tensor or its inverse for all bands at each k-point.
    @param    klist: k-point list [Nk, 3] float64
    @param    ham_r: Real-space Hamiltonian blocks [Nr, Norb, Norb] complex128
    @param     rvec: Real-space displacement vectors [Nr, 3] float64
    @param     mrot: Lattice-to-Cartesian matrix [3, 3] float64, v = mrot @ dE/dk (e.g. avec.T*ihbar)
    @param      uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param sw_imass: If True, return inverse mass (imass) [Nk, Norb, 3, 3]; if False, return mass
    @param      S_r: Real-space overlap blocks [Nr, Norb, Norb] for a non-orthogonal (MLO) basis;
                     None or [] for an orthogonal basis
    @param      eig: Eigenvalues [Nk, Norb] float64 -- REQUIRED (the interband term needs them)
    @return   mass: Effective mass or inverse mass tensor [Nk, Norb, 3, 3]
    """
    import scipy.linalg as sclin
    if eig is None:
        raise ValueError("get_mass: eig is required (the interband term of the mass needs it)")
    imass0 = get_imass0(klist, ham_r, rvec)
    vk0 = get_vlm0(klist, ham_r, rvec)
    if S_r is None or len(S_r) == 0:
        imass = get_imassk(imass0, mrot, uni, vk0, eig)
    else:
        imass = get_imassk(imass0, mrot, uni, vk0, eig,
                           get_imass0(klist, S_r, rvec), get_vlm0(klist, S_r, rvec))
    if sw_imass:
        return imass
    else:
        return np.array([[sclin.inv(im) for im in imas] for imas in imass])


def get_qshift(klist: np.ndarray, qpoint: np.ndarray) -> np.ndarray:
    """
    @brief For each k-point, find the index of k+q within the k-mesh (forward momentum shift).
    @param   klist: k-point list [Nk, 3] float64
    @param  qpoint: Momentum transfer vector [3] float64
    @return qshift: Index of k+q for each k-point [Nk] int64
    """
    Nk = len(klist)
    qshift = np.zeros(Nk, dtype=np.int64)
    _loader._lib.get_qshift_(qpoint, klist, qshift, _loader.i64(Nk))
    return qshift


def get_iqshift(klist: np.ndarray, qpoint: np.ndarray) -> np.ndarray:
    """
    @brief For each k-point, find the index of k-q within the k-mesh (backward momentum shift).
    @param   klist: k-point list [Nk, 3] float64
    @param  qpoint: Momentum transfer vector [3] float64
    @return qshift: Index of k-q for each k-point [Nk] int64
    """
    Nk = len(klist)
    qshift = np.zeros(Nk, dtype=np.int64)
    _loader._lib.get_iqshift_(qpoint, klist, qshift, _loader.i64(Nk))
    return qshift


def gen_irr_k_TRS(Nx: int, Ny: int, Nz: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    @brief Generate the irreducible k-point list for a time-reversal symmetric lattice.
    @param     Nx: Number of k-points along kx
    @param     Ny: Number of k-points along ky
    @param     Nz: Number of k-points along kz
    @return  klist: Irreducible k-point coordinates [Nk, 3] float64
    @return   kmap: Full k-grid to irreducible mapping [Nkall, 3] int64
    @return invk_ft_list: Inverse k-point list for Fortran FFT routines [Nkall, 3] int64
    """
    Nkall = Nx * Ny * Nz
    if Nkall % 2 == 0:
        if Nz % 2 == 0:
            Nk = int(Nkall / 2 + 4) if ((Nx % 2) == 0 and (Ny % 2) == 0) else (
                int(Nkall / 2 + 2) if ((Nx % 2) == 0 or (Ny % 2) == 0) else int(Nkall / 2 + 1)
            )
        else:
            Nk = int(Nkall / 2 + 2) if ((Nx % 2) == 0 and (Ny % 2) == 0) else int(Nkall / 2 + 1)
    else:
        Nk = int((Nkall + 1) / 2)
    klist = np.zeros((Nk, 3), dtype=np.float64)
    kmap = np.zeros((Nkall, 3), dtype=np.int64)
    invk_ft_list = np.zeros((Nkall, 3), dtype=np.int64)
    _loader._lib.generate_irr_kpoint_inv(klist, kmap, invk_ft_list, _loader.i64(Nk),
                                         _loader.i64(Nx), _loader.i64(Ny), _loader.i64(Nz))
    return klist, kmap, invk_ft_list


def gen_kpoint_weight(invk: np.ndarray, Nk: int) -> np.ndarray:
    """
    @brief Compute k-point weights for the irreducible BZ by counting degeneracy from the full k-grid mapping.
    @param  invk: Inverse k-point mapping [Nkall] int64
    @param    Nk: Number of irreducible k-points
    @return weight: k-point weights [Nk] float64
    """
    Nkall = len(invk)
    weight = np.zeros(Nk, dtype=np.float64)
    _loader._lib.get_kweight(weight, invk, _loader.i64(Nk), _loader.i64(Nkall))
    return weight


def mkBdGhamk(hamk: np.ndarray, delta: np.ndarray) -> np.ndarray:
    """
    @brief Build the Bogoliubov-de Gennes Hamiltonian from the normal-state Hamiltonian and gap function.
    @param   hamk: Normal-state k-space Hamiltonian [Nk, Norb, Norb] complex128
    @param  delta: Anomalous (pairing) potential [Nk, Norb, Norb] complex128
    @return hamBdGk: BdG Hamiltonian [Nk, 2*Norb, 2*Norb] complex128
    """
    Nk = len(hamk)
    Norb = hamk.shape[1]
    hamBdGk = np.zeros((Nk, 2 * Norb, 2 * Norb), dtype=np.complex128)
    _loader._lib.mkbdghamk(hamBdGk, hamk, delta, _loader.i64(Nk), _loader.i64(Norb))
    return hamBdGk


def get_plist(rvec: np.ndarray, ham_r: np.ndarray) -> np.ndarray:
    """
    @brief Determine the parity eigenvalues of each orbital from the Hamiltonian's inversion symmetry properties.
    @param  rvec: Real-space displacement vectors [Nr, 3] float64
    @param ham_r: Real-space Hamiltonian blocks [Nr, Norb, Norb] complex128
    @return plist: Parity sign (+1 or -1) for each orbital [Norb] float64
    """
    Nr = len(rvec)
    Norb = ham_r.shape[1]
    Pmn = np.zeros((Norb, Norb), dtype=np.float64)
    _loader._lib.get_parity_prop(Pmn, rvec, ham_r, _loader.i64(Norb), _loader.i64(Nr))
    return np.sign(Pmn[0, :])
