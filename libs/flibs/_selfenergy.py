import ctypes
import numpy as np
import libs.flibs._loader as _loader
import libs.flibs._freq as _freq

# Public API of this submodule; keeps `from ._x import *` in __init__.py from
# re-exporting numpy, scipy, ctypes and the sibling packages imported above.
__all__ = [
    'get_Vsigma_nosoc_flex', 'mkself', 'mkself_soc'
]

# --- ctypes signatures: set once at import.
# All Fortran entry points are subroutines, so restype is always None.
_loader._lib.get_vsigma_flex_nosoc_.argtypes = [
    _loader.c_array(np.complex128),          # chi (inout)
    _loader.c_array(np.float64),             # Smat
    _loader.c_array(np.float64),             # Cmat
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)  # Nk, Nw, Nchi
]
_loader._lib.get_vsigma_flex_nosoc_.restype = None
_loader._lib.mkself.argtypes = [
    _loader.c_array(np.complex128),             # sigmak (out)
    ctypes.POINTER(ctypes.c_double),                                       # mu (out)
    _loader.c_array(np.float64),                # Smat
    _loader.c_array(np.float64),                # Cmat
    _loader.c_array(np.int64),                  # kmap
    _loader.c_array(np.int64),                  # invk
    _loader.f_array(np.int64),                  # olist
    _loader.c_array(np.complex128),             # hamk
    _loader.c_array(np.float64),                # eig
    _loader.c_array(np.complex128),             # uni
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),                    # mu_init, rfill
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int64),                     # temp, scf_loop
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int64),  # pp, eps, Nkall
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),                      # Nk, Nw
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),                      # Norb, Nchi
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),    # Nx, Ny, Nz
    ctypes.POINTER(ctypes.c_int64),                                        # sub_sigma
    ctypes.POINTER(ctypes.c_bool), ctypes.POINTER(ctypes.c_bool),          # sw_out, sw_in
    ctypes.POINTER(ctypes.c_int64),                                        # m_diis
    ctypes.POINTER(ctypes.c_bool),                                         # sw_rescale
    ctypes.POINTER(ctypes.c_bool),                                         # sw_tail
    ctypes.POINTER(ctypes.c_double)                                        # sigma_scale
]
_loader._lib.mkself.restype = None
_loader._lib.mkself_soc.argtypes = [
    _loader.c_array(np.complex128),             # sigmak (out)
    ctypes.POINTER(ctypes.c_double),                                       # mu (out)
    _loader.c_array(np.float64),                # Vmat
    _loader.c_array(np.int64),                  # kmap
    _loader.c_array(np.int64),                  # invk
    _loader.c_array(np.int64),                  # invs
    _loader.f_array(np.int64),                  # olist
    _loader.c_array(np.int64),                  # slist
    _loader.c_array(np.complex128),             # hamk
    _loader.c_array(np.float64),                # eig
    _loader.c_array(np.complex128),             # uni
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),                    # mu_init, rfill
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int64),                     # temp, scf_loop
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int64),  # pp, eps, Nkall
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),                      # Nk, Nw
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),                      # Norb, Nchi
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),    # Nx, Ny, Nz
    ctypes.POINTER(ctypes.c_int64),                                        # sub_sigma
    ctypes.POINTER(ctypes.c_bool), ctypes.POINTER(ctypes.c_bool),          # sw_out, sw_in
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_bool),                       # m_diis, sw_rescale
    ctypes.POINTER(ctypes.c_double)                                        # sigma_scale
]
_loader._lib.mkself_soc.restype = None


def get_Vsigma_nosoc_flex(chi: np.ndarray, Smat: np.ndarray, Cmat: np.ndarray) -> np.ndarray:
    """
    @brief Apply vertex corrections for FLEX (without SOC) to the susceptibility using spin and charge matrices.
    @param   chi: Susceptibility tensor [Nk, Nw, Nchi, Nchi] complex128 (modified in-place internally)
    @param  Smat: Spin interaction matrix [Nchi, Nchi] float64
    @param  Cmat: Charge interaction matrix [Nchi, Nchi] float64
    @return  chi: Copy of the vertex-corrected susceptibility tensor
    """
    Nk, Nw, Nchi = len(chi), len(chi[0]), len(Smat)
    _loader._lib.get_vsigma_flex_nosoc_(chi, Smat, Cmat, _loader.i64(Nk),
                                        _loader.i64(Nw), _loader.i64(Nchi))
    return chi.copy()


def mkself(Smat: np.ndarray, Cmat: np.ndarray, kmap: np.ndarray, invk: np.ndarray,
           olist: np.ndarray, hamk: np.ndarray, eig: np.ndarray, uni: np.ndarray, mu: float,
           fill: float, temp: float, Nw: int, Nx: int, Ny: int, Nz: int, sw_out: bool,
           sw_in: bool, sub_sigma: int = 1, scf_loop: int = 300, eps: float = 1.0e-4,
           pp: float = 0.3, m_diis: int = 5, sw_rescale: bool = False,
           sw_tail: bool = False, sigma_scale: float = 1.0) -> tuple[np.ndarray, float]:
    """
    @brief Iteratively compute the self-energy Sigma(k, iw_n) via FLEX self-consistency loop (without SOC).
    @param         Smat: Spin interaction matrix [Nchi, Nchi] float64
    @param         Cmat: Charge interaction matrix [Nchi, Nchi] float64
    @param         kmap: Full-to-reduced k-point mapping [Nkall] int64
    @param         invk: Inverse k-point mapping [Nkall] int64
    @param        olist: Orbital-pair list of chi [Nchi, 2] int64 (1-based)
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
    @param      sw_tail: If True, use the tail-corrected chi0 (conv[G]-conv[G0]+analytic
                         reference; O(1/Nw^2) Matsubara truncation error instead of O(1/Nw)).
                         Ignored on the IR grid, which has no truncation error to correct.
    @param  sigma_scale: Factor multiplying the sigma seed loaded via sw_in (only used when
                         sw_in=True). For U-annealing set ~(U_new/U_old)^2 since Sigma ~ U^2.
    @return      sigmak: Converged self-energy [Norb, Norb, Nw, Nk] complex128
    @return     mu_self: Final chemical potential in eV
    """
    print("mixing rate: pp = %3.1f, DIIS m = %d" % (pp, m_diis), flush=True)
    if sw_tail:
        if _freq.fg_is_uniform():
            print("chi0: tail-corrected convolution (analytic G0 reference)", flush=True)
        else:
            print("chi0: sw_tail ignored on the IR grid (no Matsubara truncation to correct)", flush=True)
    Nkall, Nk, Nchi = len(kmap), len(hamk), len(Smat)
    Norb = hamk.shape[1]
    mu_self = ctypes.c_double()
    sigmak = np.zeros((Norb, Norb, Nw, Nk), dtype=np.complex128)
    _loader._lib.mkself(sigmak, ctypes.byref(mu_self), Smat, Cmat, kmap, invk, olist, hamk, eig, uni, _loader.dbl(mu),
                        _loader.dbl(fill), _loader.dbl(temp), _loader.i64(scf_loop), _loader.dbl(pp), _loader.dbl(eps),
                        _loader.i64(Nkall), _loader.i64(Nk), _loader.i64(Nw), _loader.i64(Norb), _loader.i64(Nchi),
                        _loader.i64(Nx), _loader.i64(Ny), _loader.i64(Nz), _loader.i64(sub_sigma),
                        ctypes.byref(ctypes.c_bool(sw_out)), ctypes.byref(ctypes.c_bool(sw_in)), _loader.i64(m_diis),
                        ctypes.byref(ctypes.c_bool(sw_rescale)), ctypes.byref(ctypes.c_bool(sw_tail)),
                        _loader.dbl(sigma_scale))
    return sigmak, mu_self.value


def mkself_soc(Vmat: np.ndarray, kmap: np.ndarray, invk: np.ndarray, invs: np.ndarray,
               olist: np.ndarray, slist: np.ndarray, hamk: np.ndarray, eig: np.ndarray,
               uni: np.ndarray, mu: float, fill: float, temp: float, Nw: int,
               Nx: int, Ny: int, Nz: int, sw_out: bool, sw_in: bool,
               sub_sigma: int = 1, scf_loop: int = 300, eps: float = 1.0e-4,
               pp: float = 0.3, m_diis: int = 5, sw_rescale: bool = False,
               sigma_scale: float = 1.0) -> tuple[np.ndarray, float]:
    """
    @brief Iteratively compute the self-energy Sigma(k, iw_n) via FLEX self-consistency loop with SOC.
    @param        Vmat: SOC interaction matrix [Nchi, Nchi] float64
    @param        kmap: Full-to-reduced k-point mapping [Nkall] int64
    @param        invk: Inverse k-point mapping [Nkall] int64
    @param        invs: Inverse spin index mapping [Norb] int64
    @param       olist: Orbital-pair list of chi [Nchi, 2] int64 (1-based)
    @param       slist: Spin label list for each orbital [Norb] int64
    @param        hamk: k-space Hamiltonian [Nk, Norb, Norb] complex128
    @param         eig: Eigenvalues [Nk, Norb] float64
    @param         uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param          mu: Chemical potential in eV
    @param        fill: Target electron filling
    @param        temp: Temperature in eV
    @param          Nw: Number of Matsubara frequency points
    @param    Nx,Ny,Nz: k-grid dimensions
    @param      sw_out: If True, write intermediate sigma to file
    @param       sw_in: If True, read initial sigma from file
    @param   sub_sigma: HF double-counting subtraction mode (0: none, 1: lowest Matsubara, 2: frequency-average)
    @param    scf_loop: Maximum number of self-consistency iterations
    @param         eps: Convergence tolerance
    @param          pp: Linear mixing rate (0 < pp < 1)
    @param      m_diis: DIIS history length (m_diis=1: equivalent to linear mixing)
    @param  sw_rescale: rescale chi0 when the Stoner factor >= 1
    @param sigma_scale: Factor multiplying the sigma seed loaded via sw_in (only used when
                        sw_in=True). For U-annealing set ~(U_new/U_old)^2 since Sigma ~ U^2.
    @return      sigmak: Converged self-energy [Norb, Norb, Nw, Nk] complex128
    @return     mu_self: Final chemical potential in eV
    """
    print("mixing rate: pp = %3.1f, DIIS m = %d" % (pp, m_diis), flush=True)
    Nkall, Nk, Nchi = len(kmap), len(hamk), len(Vmat)
    Norb = hamk.shape[1]
    mu_self = ctypes.c_double()
    sigmak = np.zeros((Norb, Norb, Nw, Nk), dtype=np.complex128)
    _loader._lib.mkself_soc(sigmak, ctypes.byref(mu_self), Vmat, kmap, invk, invs, olist, slist, hamk, eig, uni,
                            _loader.dbl(mu), _loader.dbl(fill), _loader.dbl(temp), _loader.i64(scf_loop),
                            _loader.dbl(pp), _loader.dbl(eps), _loader.i64(Nkall), _loader.i64(Nk), _loader.i64(Nw),
                            _loader.i64(Norb), _loader.i64(Nchi), _loader.i64(Nx), _loader.i64(Ny), _loader.i64(Nz),
                            _loader.i64(sub_sigma), ctypes.byref(ctypes.c_bool(sw_out)),
                            ctypes.byref(ctypes.c_bool(sw_in)), _loader.i64(m_diis),
                            ctypes.byref(ctypes.c_bool(sw_rescale)), _loader.dbl(sigma_scale))
    return sigmak, mu_self.value
