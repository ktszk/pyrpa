import ctypes
import numpy as np
import libs.flibs._loader as _loader

# Public API of this submodule; keeps `from ._x import *` in __init__.py from
# re-exporting numpy, scipy, ctypes and the sibling packages imported above.
__all__ = [
    'linearized_eliashberg', 'linearized_eliashberg_soc', 'nonlinear_eliashberg', 'conv_delta_orb_to_band',
    'conv_delta_orb_to_band_soc', 'gen_Fk', 'gen_Fk_soc', 'get_initial_delta', 'get_band_to_orb_delta', 'remap_gap'
]

# --- ctypes signatures: set once at import.
# All Fortran entry points are subroutines, so restype is always None.
_loader._lib.lin_eliash.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),  # delta (out)
    np.ctypeslib.ndpointer(dtype=np.complex128),  # chi (inout)
    np.ctypeslib.ndpointer(dtype=np.complex128),  # Gk
    np.ctypeslib.ndpointer(dtype=np.complex128),  # uni
    np.ctypeslib.ndpointer(dtype=np.float64),     # init_delta
    np.ctypeslib.ndpointer(dtype=np.float64),     # Smat
    np.ctypeslib.ndpointer(dtype=np.float64),     # Cmat
    np.ctypeslib.ndpointer(dtype=np.int64),       # olist
    np.ctypeslib.ndpointer(dtype=np.float64),     # prt
    np.ctypeslib.ndpointer(dtype=np.int64),       # iprm (inversion orbital permutation)
    np.ctypeslib.ndpointer(dtype=np.int64),       # kmap
    np.ctypeslib.ndpointer(dtype=np.int64),       # invk
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),         # temp, eps
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nkall, Nk
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nw, Nchi
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Norb, Nx
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Ny, Nz
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # itemax, gap_sym
    ctypes.POINTER(ctypes.c_int64),                             # arnoldi_m
    ctypes.POINTER(ctypes.c_double)                             # lambda_out (physical Eliashberg eigenvalue)
]
_loader._lib.lin_eliash.restype = None
_loader._lib.lin_eliash_soc.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),  # delta (out)
    np.ctypeslib.ndpointer(dtype=np.complex128),  # chi (inout)
    np.ctypeslib.ndpointer(dtype=np.complex128),  # Gk
    np.ctypeslib.ndpointer(dtype=np.complex128),  # uni
    np.ctypeslib.ndpointer(dtype=np.float64),     # init_delta
    np.ctypeslib.ndpointer(dtype=np.float64),     # Vmat
    np.ctypeslib.ndpointer(dtype=np.float64),     # sgnsig
    np.ctypeslib.ndpointer(dtype=np.float64),     # sgnsig2
    np.ctypeslib.ndpointer(dtype=np.float64),     # prt
    np.ctypeslib.ndpointer(dtype=np.int64),       # olist
    np.ctypeslib.ndpointer(dtype=np.int64),       # slist
    np.ctypeslib.ndpointer(dtype=np.int64),       # kmap
    np.ctypeslib.ndpointer(dtype=np.int64),       # invk
    np.ctypeslib.ndpointer(dtype=np.int64),       # invs
    np.ctypeslib.ndpointer(dtype=np.int64),       # invschi
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),         # temp, eps
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nkall, Nk
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nw, Nchi
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Norb, Nx
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Ny, Nz
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # itemax, gap_sym
    ctypes.POINTER(ctypes.c_int64),                             # arnoldi_m
    ctypes.POINTER(ctypes.c_double)                             # lambda_out (out)
]
_loader._lib.lin_eliash_soc.restype = None
_loader._lib.eliashberg.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),              # delta (inout)
    np.ctypeslib.ndpointer(dtype=np.complex128),              # sigmak (out/inout)
    np.ctypeslib.ndpointer(dtype=np.complex128),              # Gk (inout)
    np.ctypeslib.ndpointer(dtype=np.complex128),              # hamk
    np.ctypeslib.ndpointer(dtype=np.float64),                 # Smat
    np.ctypeslib.ndpointer(dtype=np.float64),                 # Cmat
    np.ctypeslib.ndpointer(dtype=np.int64),                   # olist
    np.ctypeslib.ndpointer(dtype=np.float64),                 # plist (= prt in Fortran)
    np.ctypeslib.ndpointer(dtype=np.int64),                   # iprm (inversion orbital permutation)
    np.ctypeslib.ndpointer(dtype=np.int64),                   # kmap
    np.ctypeslib.ndpointer(dtype=np.int64),                   # invk
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),  # mu, temp, eps
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),     # Nkall, Nk, Nw
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),     # Nchi, Norb, Nx
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),     # Ny, Nz, itemax
    ctypes.POINTER(ctypes.c_int64),                                         # gap_sym
    ctypes.POINTER(ctypes.c_bool), ctypes.POINTER(ctypes.c_bool),                         # sw_sigma,sw_Vconst
    ctypes.POINTER(ctypes.c_int64),                                         # m_diis
    ctypes.POINTER(ctypes.c_double),                                        # gap_min
    ctypes.POINTER(ctypes.c_bool),                                          # sw_amp_newton
]
_loader._lib.eliashberg.restype = None
_loader._lib.conv_delta_orb_to_band.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),          # deltab (out)
    np.ctypeslib.ndpointer(dtype=np.complex128),          # delta
    np.ctypeslib.ndpointer(dtype=np.complex128),          # uni
    np.ctypeslib.ndpointer(dtype=np.float64),             # prt
    np.ctypeslib.ndpointer(dtype=np.int64),               # iprm (inversion orbital permutation)
    np.ctypeslib.ndpointer(dtype=np.int64),               # invk
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),                   # Norb, Nkall
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)  # Nk, Nw, gap_sym
]
_loader._lib.conv_delta_orb_to_band.restype = None
_loader._lib.conv_delta_orb_to_band_soc.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),  # deltab (out)
    np.ctypeslib.ndpointer(dtype=np.complex128),  # delta
    np.ctypeslib.ndpointer(dtype=np.complex128),  # uni
    np.ctypeslib.ndpointer(dtype=np.int64),       # invk
    np.ctypeslib.ndpointer(dtype=np.int64),       # invs
    np.ctypeslib.ndpointer(dtype=np.int64),       # slist
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Norb, Nkall
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)            # Nk, Nw
]
_loader._lib.conv_delta_orb_to_band_soc.restype = None
_loader._lib.mkfk_trs_nsoc_.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),          # fk (out)
    np.ctypeslib.ndpointer(dtype=np.complex128),          # Gk
    np.ctypeslib.ndpointer(dtype=np.complex128),          # delta
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)  # Nk, Nw, Norb
]
_loader._lib.mkfk_trs_nsoc_.restype = None
_loader._lib.mkfk_trs_soc_.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),  # fk (out)
    np.ctypeslib.ndpointer(dtype=np.complex128),  # Gk
    np.ctypeslib.ndpointer(dtype=np.complex128),  # delta
    np.ctypeslib.ndpointer(dtype=np.float64),     # sgnsig
    np.ctypeslib.ndpointer(dtype=np.int64),       # slist
    np.ctypeslib.ndpointer(dtype=np.int64),       # invk
    np.ctypeslib.ndpointer(dtype=np.int64),       # invs
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nkall, Nk
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nw, Norb
    ctypes.POINTER(ctypes.c_int64)                              # gap_sym
]
_loader._lib.mkfk_trs_soc_.restype = None
_loader._lib.get_initial_delta_.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),  # delta (out)
    np.ctypeslib.ndpointer(dtype=np.float64),     # init_delta
    np.ctypeslib.ndpointer(dtype=np.complex128),  # uni
    np.ctypeslib.ndpointer(dtype=np.int64),       # kmap
    np.ctypeslib.ndpointer(dtype=np.int64),       # invk
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nkall, Nk
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nw, Norb
    ctypes.POINTER(ctypes.c_int64),                             # gap_sym
]
_loader._lib.get_initial_delta_.restype = None
_loader._lib.get_band_to_orb_delta.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),  # delta (out)
    np.ctypeslib.ndpointer(dtype=np.complex128),  # init_delta
    np.ctypeslib.ndpointer(dtype=np.complex128),  # uni
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nk, Norb
]
_loader._lib.get_band_to_orb_delta.restype = None
_loader._lib.remap_delta.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.complex128),  # delta (out)
    np.ctypeslib.ndpointer(dtype=np.complex128),  # delta0
    np.ctypeslib.ndpointer(dtype=np.float64),     # prt
    np.ctypeslib.ndpointer(dtype=np.int64),       # iprm (inversion orbital permutation)
    np.ctypeslib.ndpointer(dtype=np.int64),       # invk
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nkall, Nk
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),           # Nw, Norb
    ctypes.POINTER(ctypes.c_int64)                              # gap_sym
]
_loader._lib.remap_delta.restype = None


def linearized_eliashberg(chi: np.ndarray, Gk: np.ndarray, uni: np.ndarray, init_delta: np.ndarray,
                          Smat: np.ndarray, Cmat: np.ndarray, olist: np.ndarray, plist: np.ndarray,
                          kmap: np.ndarray, invk: np.ndarray, Nx: int, Ny: int, Nz: int,
                          temp: float, gap_sym: int, eps: float = 1.0e-5, itemax: int = 300,
                          arnoldi_m: int = 10, iperm: np.ndarray | None = None) -> tuple[np.ndarray, float]:
    """
    @fn linearized_eliashberg
    @brief Solve the linearized Eliashberg gap equation (without SOC) to obtain the superconducting gap function.
    @param       chi: Irreducible susceptibility [Nchi, Nchi, Nw, Nk] complex128
    @param        Gk: Green's function on full k-grid [Norb, Norb, Nw, Nkall] complex128
    @param       uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param init_delta: Initial gap function [Norb, Norb, Nw, Nkall] float64
    @param      Smat: Spin interaction matrix [Nchi, Nchi] float64
    @param      Cmat: Charge interaction matrix [Nchi, Nchi] float64
    @param     olist: Orbital index list [Nchi] int64
    @param     plist: Parity list for gap symmetry [Norb] float64
    @param      kmap: Full-to-reduced k-point mapping [Nkall] int64
    @param      invk: Inverse k-point mapping [Nkall] int64
    @param    Nx,Ny,Nz: k-grid dimensions
    @param      temp: Temperature in eV
    @param   gap_sym: Gap symmetry index (0=s, 1=dx2-y2, etc.)
    @param       eps: Convergence tolerance
    @param    itemax: Maximum number of iterations
    @param     iperm: Orbital permutation of the inversion operator [Norb] int64.
                      None = identity, i.e. the diagonal-parity convention. Use
                      plibs.inversion_op() to obtain (plist, iperm) for multi-site cells.
    @return    delta: Linearized gap function [Norb, Norb, Nw, Nkall] complex128
    @return   lambda_eliash: Physical (largest) Eliashberg eigenvalue (>1 means T<Tc)
    """
    Norb, Nchi = len(Gk), len(Smat)
    Nkall, Nk, Nw = len(kmap), len(Gk[0, 0, 0]), len(Gk[0, 0])
    delta = np.zeros((Norb, Norb, Nw, Nk), dtype=np.complex128)
    lambda_out = ctypes.c_double(0.0)
    iprm = np.arange(Norb, dtype=np.int64) if iperm is None else np.asarray(iperm, dtype=np.int64)
    iprm = iprm + 1   # Fortran is 1-based
    _loader._lib.lin_eliash(delta, chi, Gk, uni, init_delta, Smat, Cmat, olist, plist, iprm, kmap, invk,
                            _loader.dbl(temp), _loader.dbl(eps), _loader.i64(Nkall),
                            _loader.i64(Nk), _loader.i64(Nw), _loader.i64(Nchi),
                            _loader.i64(Norb), _loader.i64(Nx), _loader.i64(Ny),
                            _loader.i64(Nz), _loader.i64(itemax), _loader.i64(gap_sym),
                            _loader.i64(arnoldi_m), ctypes.byref(lambda_out))
    return delta, lambda_out.value


def linearized_eliashberg_soc(chi: np.ndarray, Gk: np.ndarray, uni: np.ndarray, init_delta: np.ndarray,
                              Vmat: np.ndarray, sgnsig: np.ndarray, sgnsig2: np.ndarray, plist: np.ndarray,
                              slist: np.ndarray, olist: np.ndarray, kmap: np.ndarray, invk: np.ndarray,
                              invs: np.ndarray, invschi: np.ndarray, Nx: int, Ny: int, Nz: int, temp: float,
                              gap_sym: int, eps: float = 1.0e-4, itemax: int = 300,
                              arnoldi_m: int = 10) -> tuple[np.ndarray, float]:
    """
    @fn linearized_eliashberg_soc
    @brief Solve the linearized Eliashberg gap equation with spin-orbit coupling (SOC).
    @param        chi: Irreducible susceptibility [Nchi, Nchi, Nw, Nk] complex128
    @param         Gk: Green's function on full k-grid [Norb, Norb, Nw, Nkall] complex128
    @param        uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param init_delta: Initial gap function [Norb, Norb, Nw, Nkall] float64
    @param       Vmat: SOC interaction matrix [Nchi, Nchi] float64
    @param     sgnsig: Sign matrix for spin-spin product [Norb, Norb] float64
    @param    sgnsig2: Sign matrix for chi orbital mapping [Nchi, Nchi] float64
    @param      plist: Parity list for gap symmetry [Norb] float64
    @param      slist: Spin label list [Norb] int64
    @param      olist: Orbital index list [Nchi] int64
    @param       kmap: Full-to-reduced k-point mapping [Nkall] int64
    @param       invk: Inverse k-point mapping [Nkall] int64
    @param       invs: Inverse spin index mapping [Norb] int64
    @param    invschi: Inverse chi orbital mapping [Nchi] int64
    @param   Nx,Ny,Nz: k-grid dimensions
    @param       temp: Temperature in eV
    @param    gap_sym: Gap symmetry index
    @param        eps: Convergence tolerance
    @param     itemax: Maximum number of iterations
    @return     delta: Linearized gap function [Norb, Norb, Nw, Nkall] complex128
    @return   lambda_eliash: Physical (largest) Eliashberg eigenvalue (>1 means T<Tc)
    """
    Norb, Nchi = len(slist), len(Vmat)
    Nkall, Nk, Nw = len(kmap), len(Gk[0, 0, 0]), len(Gk[0, 0])
    delta = np.zeros((Norb, Norb, Nw, Nkall), dtype=np.complex128)
    lambda_out = ctypes.c_double(0.0)
    _loader._lib.lin_eliash_soc(delta, chi, Gk, uni, init_delta, Vmat, sgnsig, sgnsig2, plist, olist,
                                slist, kmap, invk, invs, invschi, _loader.dbl(temp),
                                _loader.dbl(eps), _loader.i64(Nkall), _loader.i64(Nk),
                                _loader.i64(Nw), _loader.i64(Nchi), _loader.i64(Norb),
                                _loader.i64(Nx), _loader.i64(Ny), _loader.i64(Nz),
                                _loader.i64(itemax), _loader.i64(gap_sym),
                                _loader.i64(arnoldi_m), ctypes.byref(lambda_out))
    return delta, lambda_out.value


def nonlinear_eliashberg(delta_init: np.ndarray, Gk: np.ndarray, hamk: np.ndarray,
                         Smat: np.ndarray, Cmat: np.ndarray,
                         olist: np.ndarray, plist: np.ndarray, kmap: np.ndarray,
                         invk: np.ndarray, mu: float, temp: float, gap_sym: int,
                         Nx: int, Ny: int, Nz: int,
                         sw_sigma: bool = False, sw_Vconst: bool = False, eps: float = 1.0e-4,
                         itemax: int = 100, m_diis: int = 5,
                         gap_min: float = 1.0e-2,
                         sw_amp_newton: bool = True, iperm: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """
    @fn nonlinear_eliashberg
    @brief Solve the non-linear FLEX-Eliashberg self-consistency loop (no SOC).
    Iteratively updates the gap Δ, normal/anomalous Green's functions G/F,
    irreducible χ_0_S = χ^GG + χ^FF, χ_0_C = χ^GG - χ^FF, and the FLEX
    vertices V_σ, V_Δ. Pulay/DIIS extrapolation is used for Δ mixing
    (linear mixing pp=0.3 falls back when m_diis<=1 or for the first iteration).
    Convergence: |Δ_new - Δ|_∞ / |Δ|_∞ < eps.

    @param  delta_init: Initial gap function [Norb, Norb, Nw, Nk] complex128.
                        Caller is responsible for providing a non-zero initial Δ.
                        Recommended sources:
                          - flibs.get_initial_delta(...) for a symmetry-adapted seed
                          - a converged linearized_eliashberg() output (scaled to a
                            physical magnitude, e.g. ~1.764·T) for a fast start
                        If the array is C-contiguous complex128 it is modified in place;
                        otherwise a contiguous copy is created.
    @param          Gk: Initial normal Green's function [Norb, Norb, Nw, Nk] complex128
                        (typically from flibs.gen_Green0 or a converged normal-state solution)
    @param        hamk: k-space Hamiltonian [Nk, Norb, Norb] complex128
    @param        Smat: Spin interaction matrix [Nchi, Nchi] float64
    @param        Cmat: Charge interaction matrix [Nchi, Nchi] float64
    @param       olist: Orbital index pairs for chi [Nchi, 2] int64
    @param       plist: Parity factors per orbital [Norb] float64
    @param        kmap: Full k-grid to irreducible-grid mapping [Nkall, 3] int64
    @param        invk: Inverse k-point list [Nkall, 3] int64
    @param          mu: Chemical potential in eV
    @param        temp: Temperature in eV
    @param     gap_sym: Gap symmetry (>=0 singlet, <0 triplet)
    @param    Nx,Ny,Nz: k-grid dimensions
    @param    sw_sigma: Include FLEX self-energy correction in the dressed G_0 (default False)
    @param    sw_Vconst: Use constant pairing interaction instead of FLEX vertices (default False)
    @param         eps: Convergence threshold on |Δ_new-Δ|/|Δ|
    @param      itemax: Maximum SCF iterations
    @param      m_diis: DIIS history length (m_diis<=1 → linear mixing only; default 5)
    @param     gap_min: Gap-collapse cutoff as a fraction of the initial |Δ|max (default 1e-2).
                        The loop stops early once |Δ|max falls below gap_min·|Δ_init|max,
                        i.e. the gap is decaying toward the trivial solution (T >= Tc).
                        Set gap_min<=0 to disable the early stop.
    @param sw_amp_newton: Enable the amplitude accelerator with adaptive anchoring (default True).
                        Near the bifurcated SC fixed point the slow mode is the gap amplitude.
                        The gap shape is anchored to the (normalized) seed eigenvector v0 (Δ=a·v0),
                        which suppresses drift onto the unphysical negative eigenmode; the v0-
                        projected gain g(a)=<v0,M[a v0]>/a follows the pitchfork form g(a)=λ0-c·a²,
                        and the amplitude is solved by the cubic root a*=sqrt((λ0-1)/c) (Phase 1).
                        Once |g-1|<2% it switches to DIIS to relax the shape (Phase 2), and after a
                        few polish steps re-anchors v0 to the relaxed shape and re-solves a
                        (adaptive anchoring), alternating until converged. Set False to use plain
                        DIIS/linear mixing throughout (uses m_diis).
    @retval      delta: Converged gap function [Norb, Norb, Nw, Nk] complex128
                        (same array as delta_init when contiguous; otherwise a new copy)
    @retval     sigmak: FLEX self-energy [Norb, Norb, Nw, Nk] complex128
                        (zero array when sw_sigma=False)
    """
    if delta_init.dtype != np.complex128 or not delta_init.flags['C_CONTIGUOUS']:
        delta = np.ascontiguousarray(delta_init, dtype=np.complex128)
    else:
        delta = delta_init
    Norb, _, Nw, Nk = delta.shape
    Nkall, Nchi = len(kmap), len(Smat)
    sigmak = np.zeros((Norb, Norb, Nw, Nk), dtype=np.complex128)
    iprm = np.arange(Norb, dtype=np.int64) if iperm is None else np.asarray(iperm, dtype=np.int64)
    iprm = iprm + 1   # Fortran is 1-based
    _loader._lib.eliashberg(delta, sigmak, Gk, hamk, Smat, Cmat, olist, plist, iprm, kmap, invk, _loader.dbl(mu),
                            _loader.dbl(temp), _loader.dbl(eps), _loader.i64(Nkall), _loader.i64(Nk), _loader.i64(Nw),
                            _loader.i64(Nchi), _loader.i64(Norb), _loader.i64(Nx), _loader.i64(Ny), _loader.i64(Nz),
                            _loader.i64(itemax), _loader.i64(gap_sym), ctypes.byref(ctypes.c_bool(sw_sigma)),
                            ctypes.byref(ctypes.c_bool(sw_Vconst)), _loader.i64(m_diis), _loader.dbl(gap_min),
                            ctypes.byref(ctypes.c_bool(sw_amp_newton)))
    return delta, sigmak


def conv_delta_orb_to_band(delta: np.ndarray, uni: np.ndarray, invk: np.ndarray,
                           plist: np.ndarray, gap_sym: int | str | None, iperm: np.ndarray | None = None) -> np.ndarray:
    """
    @fn conv_delta_orb_to_band
    @brief Convert the superconducting gap function from orbital basis to band basis (without SOC).
    @param   delta: Gap function in orbital basis [Norb, Norb, Nw, Nkall] complex128
    @param     uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param    invk: Inverse k-point mapping [Nkall] int64
    @param   plist: Parity list for gap symmetry [Norb] float64
    @param gap_sym: Gap symmetry index
    @return deltab: Gap function in band basis [Norb, Norb, Nkall] complex128
    """
    Nkall, Nk, Nw, Norb = len(invk), len(uni), len(delta[0, 0]), len(delta)
    deltab = np.zeros((Norb, Norb, Nkall), dtype=np.complex128)
    iprm = np.arange(Norb, dtype=np.int64) if iperm is None else np.asarray(iperm, dtype=np.int64)
    iprm = iprm + 1   # Fortran is 1-based
    _loader._lib.conv_delta_orb_to_band(deltab, delta, uni, plist, iprm, invk, _loader.i64(Norb),
                                        _loader.i64(Nkall), _loader.i64(Nk),
                                        _loader.i64(Nw), _loader.i64(gap_sym))
    return deltab


def conv_delta_orb_to_band_soc(delta: np.ndarray, uni: np.ndarray, invk: np.ndarray,
                               invs: np.ndarray, slist: np.ndarray) -> np.ndarray:
    """
    @fn conv_delta_orb_to_band_soc
    @brief Convert the superconducting gap function from orbital basis to band basis with SOC.
    @param   delta: Gap function in orbital basis [Norb, Norb, Nw, Nkall] complex128
    @param     uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param    invk: Inverse k-point mapping [Nkall] int64
    @param    invs: Inverse spin index mapping [Norb] int64
    @param   slist: Spin-orbit label for each orbital [Norb] int64
    @return deltab: Gap function in band basis [Norb, Norb, Nkall] complex128
    """
    Nkall, Nk, Nw, Norb = len(invk), len(uni), len(delta[0, 0]), len(delta)
    deltab = np.zeros((Norb, Norb, Nkall), dtype=np.complex128)
    _loader._lib.conv_delta_orb_to_band_soc(deltab, delta, uni, invk, invs, slist, _loader.i64(Norb),
                                            _loader.i64(Nkall), _loader.i64(Nk), _loader.i64(Nw))
    return deltab


def gen_Fk(Gk: np.ndarray, delta: np.ndarray, invk: np.ndarray) -> np.ndarray:
    """
    @fn gen_Fk
    @brief Generate the anomalous Green's function Fk from G and delta (without SOC).
    @param     Gk: Normal Green's function [Norb, Norb, Nw, Nk] complex128
    @param  delta: Anomalous self-energy [Norb, Norb, Nw, Nkall] complex128
    @param   invk: Inverse k-point mapping [Nkall] int64
    @return    Fk: Anomalous Green's function [Norb, Norb, Nw, Nk] complex128
    """
    Nkall, Nk, Nw, Norb = len(invk), len(Gk[0, 0, 0]), len(delta[0, 0]), len(delta)
    Fk = np.zeros((Norb, Norb, Nw, Nk), dtype=np.complex128)
    _loader._lib.mkfk_trs_nsoc_(Fk, Gk, delta, _loader.i64(Nk), _loader.i64(Nw), _loader.i64(Norb))
    return Fk


def gen_Fk_soc(Gk: np.ndarray, delta: np.ndarray, invk: np.ndarray, invs: np.ndarray,
               slist: np.ndarray, gap_sym: int) -> np.ndarray:
    """
    @fn gen_Fk_soc
    @brief Generate the anomalous Green's function Fk from G and delta with SOC.
    @param     Gk: Normal Green's function [Norb, Norb, Nw, Nk] complex128
    @param  delta: Anomalous self-energy [Norb, Norb, Nw, Nkall] complex128
    @param   invk: Inverse k-point mapping [Nkall] int64
    @param   invs: Inverse spin index mapping [Norb] int64
    @param  slist: Spin-orbit label for each orbital [Norb] int64
    @param gap_sym: Gap symmetry index
    @return    Fk: Anomalous Green's function [Norb, Norb, Nw, Nkall] complex128
    """
    Nkall, Nk, Nw, Norb = len(invk), len(Gk[0, 0, 0]), len(delta[0, 0]), len(delta)
    Fk = np.zeros((Norb, Norb, Nw, Nkall), dtype=np.complex128)
    sgnsig = np.array([slist]).T.dot(np.array([slist])).astype(np.float64)
    _loader._lib.mkfk_trs_soc_(Fk, Gk, delta, sgnsig, slist, invk, invs, _loader.i64(Nkall),
                               _loader.i64(Nk), _loader.i64(Nw), _loader.i64(Norb),
                               _loader.i64(gap_sym))
    return Fk


def get_initial_delta(init_delta: np.ndarray, uni: np.ndarray, kmap: np.ndarray,
                      invk: np.ndarray, Nw: int, gap_sym: int) -> np.ndarray:
    """
    @fn get_initial_delta
    @brief Convert a band-basis initial gap function into an orbital-basis gap function,
    broadcast over Matsubara frequencies, and normalize.
    Applies U * diag(init_delta) * U† at each k-point, then normalises by the total norm.
    @param  init_delta: Band-basis initial gap amplitudes [Norb, Nk] float64
                        (matches get_initial_gap output; C-order [Norb,Nk] == Fortran dimension(Nk,Norb))
    @param         uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param        kmap: Full k-grid to irreducible k-point mapping [Nkall, 3] int64
    @param        invk: Inverse k-point list [Nkall, 3] int64
    @param          Nw: Number of Matsubara frequency points
    @param     gap_sym: Gap symmetry index (0=s-wave diagonal, >0=singlet, <0=triplet)
    @return      delta: Orbital-basis gap function [Norb, Norb, Nw, Nk] complex128
    """
    Norb, Nk = init_delta.shape
    Nkall = len(kmap)
    delta = np.zeros((Norb, Norb, Nw, Nk), dtype=np.complex128)
    _loader._lib.get_initial_delta_(delta, init_delta, uni, kmap, invk, _loader.i64(Nkall), _loader.i64(Nk),
                                    _loader.i64(Nw), _loader.i64(Norb), _loader.i64(gap_sym))
    return delta


def get_band_to_orb_delta(init_delta: np.ndarray, uni: np.ndarray) -> np.ndarray:
    """
    @fn get_band_to_orb_delta
    @brief Transform the gap function from band-diagonal representation to orbital basis.
    The eigenvector matrices returned by get_eig/get_eigs are stored as row-eigenvectors,
    i.e. U(k) H(k) = E(k) U(k). Under this convention used with mkBdGhamk,
    the orbital-basis pairing matrix is
    delta_orb(k) = U(k)^T * diag(delta_n(k)) * U(k).
    This is the inverse of the band-projection in get_initial_delta / conv_delta_orb_to_band:
    given a band-diagonal gap it returns the full Norb×Norb orbital-basis matrix.
    @param init_delta: Band-diagonal gap amplitudes [Norb, Nk] complex128
                       (C-order [Norb, Nk] == Fortran dimension(Nk, Norb))
    @param       uni: Row-eigenvector matrices [Nk, Norb, Norb] complex128
    @return    delta: Gap function in orbital basis [Nk, Norb, Norb] complex128
                      (same layout as hamk from gen_ham; matches mkBdGhamk delta argument)
    """
    Norb, Nk = init_delta.shape
    # Fortran routine expects init_delta(Nk,Norb), uni(Norb,Norb,Nk),
    # and returns delta(Norb,Norb,Nk).
    init_delta_f = np.asfortranarray(init_delta.T)
    uni_f = np.asfortranarray(np.transpose(uni, (1, 2, 0)))
    delta_f = np.zeros((Norb, Norb, Nk), dtype=np.complex128, order='F')
    _loader._lib.get_band_to_orb_delta(delta_f, init_delta_f, uni_f, _loader.i64(Nk), _loader.i64(Norb))
    return np.ascontiguousarray(np.transpose(delta_f, (2, 0, 1)))


def remap_gap(delta0: np.ndarray, plist: np.ndarray, invk: np.ndarray, gap_sym: int | str | None,
              iperm: np.ndarray | None = None) -> np.ndarray:
    """
    @fn remap_gap
    @brief Remap the anomalous self-energy from the reduced k-grid to the full k-grid applying gap symmetry.
    @param  delta0: Anomalous self-energy on reduced k-grid [Norb, Norb, Nw, Nk] complex128
    @param   plist: Parity list for gap symmetry [Norb] float64
    @param    invk: Inverse k-point mapping [Nkall] int64
    @param gap_sym: Gap symmetry index
    @return  delta: Anomalous self-energy on full k-grid [Norb, Norb, Nw, Nkall] complex128
    """
    Nkall, Nk, Norb = len(invk), len(delta0.T), len(plist)
    Nw = delta0.shape[2]
    delta = np.zeros((Norb, Norb, Nw, Nkall), dtype=np.complex128)
    iprm = np.arange(Norb, dtype=np.int64) if iperm is None else np.asarray(iperm, dtype=np.int64)
    iprm = iprm + 1   # Fortran is 1-based
    _loader._lib.remap_delta(delta, delta0, plist, iprm, invk, _loader.i64(Nkall), _loader.i64(Nk),
                             _loader.i64(Nw), _loader.i64(Norb), _loader.i64(gap_sym))
    return delta
