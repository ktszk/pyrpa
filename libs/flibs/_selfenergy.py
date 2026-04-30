from ctypes import *
import numpy as np
from ._loader import _lib

def get_Vsigma_nosoc_flex(chi: np.ndarray, Smat: np.ndarray, Cmat: np.ndarray) -> np.ndarray:
    """
    @fn get_Vsigma_nosoc_flex
    @brief Apply vertex corrections for FLEX (without SOC) to the susceptibility using spin and charge matrices.
    @param   chi: Susceptibility tensor [Nk, Nw, Nchi, Nchi] complex128 (modified in-place internally)
    @param  Smat: Spin interaction matrix [Nchi, Nchi] float64
    @param  Cmat: Charge interaction matrix [Nchi, Nchi] float64
    @return  chi: Copy of the vertex-corrected susceptibility tensor
    """
    Nk, Nw, Nchi = len(chi), len(chi[0]), len(Smat)
    _lib.get_vsigma_flex_nosoc_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.get_vsigma_flex_nosoc_.restype = c_void_p
    _lib.get_vsigma_flex_nosoc_(chi, Smat, Cmat, byref(c_int64(Nk)),
                                byref(c_int64(Nw)), byref(c_int64(Nchi)))
    return chi.copy()

def mkself(Smat: np.ndarray, Cmat: np.ndarray, kmap: np.ndarray, invk: np.ndarray,
           olist: np.ndarray, hamk: np.ndarray, eig: np.ndarray, uni: np.ndarray, mu: float,
           fill: float, temp: float, Nw: int, Nx: int, Ny: int, Nz: int, sw_out: bool,
           sw_in: bool, sub_sigma: int = 1, scf_loop: int = 300, eps: float = 1.0e-4,
           pp: float = 0.3, m_diis: int = 5, sw_rescale: bool = False) -> tuple[np.ndarray, float]:
    """
    @fn mkself
    @brief Iteratively compute the self-energy Sigma(k, iw_n) via FLEX self-consistency loop (without SOC).
    @param         Smat: Spin interaction matrix [Nchi, Nchi] float64
    @param         Cmat: Charge interaction matrix [Nchi, Nchi] float64
    @param         kmap: Full-to-reduced k-point mapping [Nkall] int64
    @param         invk: Inverse k-point mapping [Nkall] int64
    @param        olist: Orbital index list [Nchi] int64
    @param         hamk: k-space Hamiltonian [Nk, Norb, Norb] complex128
    @param          eig: Eigenvalues [Nk, Norb] float64
    @param          uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param           mu: Chemical potential in eV
    @param         fill: Target electron filling
    @param         temp: Temperature in eV
    @param           Nw: Number of Matsubara frequency points
    @param     Nx,Ny,Nz: k-grid dimensions
    @param       sw_out: If True, write intermediate sigma to file
    @param        sw_in: If True, read initial sigma from file
    @param   sub_sigma: HF double-counting subtraction mode (0: none, 1: lowest Matsubara, 2: frequency-average)
    @param     scf_loop: Maximum number of self-consistency iterations
    @param          eps: Convergence tolerance
    @param           pp: Linear mixing rate (0 < pp < 1)
    @param       m_diis: DIIS history length (m_diis=1: equivalent to linear mixing)
    @param   sw_rescale: If True, rescale chi0 when Stoner factor >= 1
    @retval      sigmak: Converged self-energy [Norb, Norb, Nw, Nk] complex128
    @retval     mu_self: Final chemical potential in eV
    """
    print("mixing rate: pp = %3.1f, DIIS m = %d" % (pp, m_diis))
    Nkall, Nk, Nchi = len(kmap), len(hamk), len(Smat)
    Norb = int(np.sqrt(hamk.size / Nk))
    mu_self = c_double()
    sigmak = np.zeros((Norb, Norb, Nw, Nk), dtype=np.complex128)
    _lib.mkself.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        POINTER(c_double),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        POINTER(c_double), POINTER(c_double),
        POINTER(c_double), POINTER(c_int64),
        POINTER(c_double), POINTER(c_double), POINTER(c_int64),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64), POINTER(c_bool), POINTER(c_bool),
        POINTER(c_int64),
        POINTER(c_bool)
    ]
    _lib.mkself.restype = c_void_p
    _lib.mkself(sigmak, byref(mu_self), Smat, Cmat, kmap, invk, olist, hamk, eig, uni,
        byref(c_double(mu)), byref(c_double(fill)), byref(c_double(temp)), byref(c_int64(scf_loop)),
        byref(c_double(pp)), byref(c_double(eps)), byref(c_int64(Nkall)), byref(c_int64(Nk)),
        byref(c_int64(Nw)), byref(c_int64(Norb)), byref(c_int64(Nchi)), byref(c_int64(Nx)),
        byref(c_int64(Ny)), byref(c_int64(Nz)), byref(c_int64(sub_sigma)), byref(c_bool(sw_out)),
        byref(c_bool(sw_in)), byref(c_int64(m_diis)), byref(c_bool(sw_rescale)))
    return sigmak, mu_self.value

def mkself_soc(Vmat: np.ndarray, kmap: np.ndarray, invk: np.ndarray, invs: np.ndarray,
               olist: np.ndarray, slist: np.ndarray, hamk: np.ndarray, eig: np.ndarray,
               uni: np.ndarray, mu: float, fill: float, temp: float, Nw: int,
               Nx: int, Ny: int, Nz: int, sw_out: bool, sw_in: bool,
               sub_sigma: int = 1, scf_loop: int = 300, eps: float = 1.0e-4,
               pp: float = 0.3, m_diis: int = 5, sw_rescale: bool = False) -> tuple[np.ndarray, float]:
    """
    @fn mkself_soc
    @brief Iteratively compute the self-energy Sigma(k, iw_n) via FLEX self-consistency loop with SOC.
    @param         Vmat: SOC interaction matrix [Nchi, Nchi] float64
    @param         kmap: Full-to-reduced k-point mapping [Nkall] int64
    @param         invk: Inverse k-point mapping [Nkall] int64
    @param         invs: Inverse spin index mapping [Norb] int64
    @param        olist: Orbital index list [Nchi] int64
    @param        slist: Spin label list for each orbital [Norb] int64
    @param         hamk: k-space Hamiltonian [Nk, Norb, Norb] complex128
    @param          eig: Eigenvalues [Nk, Norb] float64
    @param          uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param           mu: Chemical potential in eV
    @param         fill: Target electron filling
    @param         temp: Temperature in eV
    @param           Nw: Number of Matsubara frequency points
    @param     Nx,Ny,Nz: k-grid dimensions
    @param       sw_out: If True, write intermediate sigma to file
    @param        sw_in: If True, read initial sigma from file
    @param   sub_sigma: HF double-counting subtraction mode (0: none, 1: lowest Matsubara, 2: frequency-average)
    @param     scf_loop: Maximum number of self-consistency iterations
    @param          eps: Convergence tolerance
    @param           pp: Linear mixing rate (0 < pp < 1)
    @param       m_diis: DIIS history length (m_diis=1: equivalent to linear mixing)
    @retval      sigmak: Converged self-energy [Norb, Norb, Nw, Nk] complex128
    @retval     mu_self: Final chemical potential in eV
    """
    print("mixing rate: pp = %3.1f, DIIS m = %d" % (pp, m_diis))
    Nkall, Nk, Nchi = len(kmap), len(hamk), len(Vmat)
    Norb = int(np.sqrt(hamk.size / Nk))
    mu_self = c_double()
    sigmak = np.zeros((Norb, Norb, Nw, Nk), dtype=np.complex128)
    _lib.mkself_soc.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        POINTER(c_double),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        POINTER(c_double), POINTER(c_double),
        POINTER(c_double), POINTER(c_int64),
        POINTER(c_double), POINTER(c_double), POINTER(c_int64),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64), POINTER(c_bool), POINTER(c_bool),
        POINTER(c_int64), POINTER(c_bool)
    ]
    _lib.mkself_soc(sigmak, byref(mu_self), Vmat, kmap, invk, invs, olist, slist, hamk, eig, uni,
               byref(c_double(mu)), byref(c_double(fill)), byref(c_double(temp)), byref(c_int64(scf_loop)),
               byref(c_double(pp)), byref(c_double(eps)), byref(c_int64(Nkall)), byref(c_int64(Nk)),
               byref(c_int64(Nw)), byref(c_int64(Norb)), byref(c_int64(Nchi)), byref(c_int64(Nx)),
               byref(c_int64(Ny)), byref(c_int64(Nz)), byref(c_int64(sub_sigma)), byref(c_bool(sw_out)),
               byref(c_bool(sw_in)), byref(c_int64(m_diis)), byref(c_bool(sw_rescale)))
    return sigmak, mu_self.value
