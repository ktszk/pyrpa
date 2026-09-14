import ctypes
import numpy as np
import libs.flibs._loader as _loader

# Public API of this submodule; keeps `from ._x import *` in __init__.py from
# re-exporting numpy, scipy, ctypes and the sibling packages imported above.
__all__ = [
    'gen_SCmatrix', 'gen_SCmatrix_orb', 'gen_Vmatrix', 'gen_Vmatrix_orb'
]

# --- ctypes signatures: set once at import.
# All Fortran entry points are subroutines, so restype is always None.
_loader._lib.get_scmat.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64),               # Smat (out)
    np.ctypeslib.ndpointer(dtype=np.float64),               # Cmat (out)
    np.ctypeslib.ndpointer(dtype=np.int64),                 # ol
    np.ctypeslib.ndpointer(dtype=np.int64),                 # site
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int64)  # Uval, Jval, Nchi
]
_loader._lib.get_scmat.restype = None
_loader._lib.get_scmat_orb.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64),  # Smat (out)
    np.ctypeslib.ndpointer(dtype=np.float64),  # Cmat (out)
    np.ctypeslib.ndpointer(dtype=np.int64),    # ol
    np.ctypeslib.ndpointer(dtype=np.int64),    # site
    np.ctypeslib.ndpointer(dtype=np.float64),  # Umat
    np.ctypeslib.ndpointer(dtype=np.float64),  # Jmat
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)         # Nchi, Norb
]
_loader._lib.get_scmat_orb.restype = None
_loader._lib.get_vmat_soc.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64),  # Vmat (out)
    np.ctypeslib.ndpointer(dtype=np.int64),    # ol
    np.ctypeslib.ndpointer(dtype=np.int64),    # sl
    np.ctypeslib.ndpointer(dtype=np.int64),    # site
    np.ctypeslib.ndpointer(dtype=np.int64),    # invs
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),      # Uval, Jval
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)         # Nchi, Norb
]
_loader._lib.get_vmat_soc.restype = None
_loader._lib.get_vmat_soc_orb.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64),  # Vmat (out)
    np.ctypeslib.ndpointer(dtype=np.int64),    # ol
    np.ctypeslib.ndpointer(dtype=np.int64),    # sl
    np.ctypeslib.ndpointer(dtype=np.int64),    # site
    np.ctypeslib.ndpointer(dtype=np.int64),    # invs
    np.ctypeslib.ndpointer(dtype=np.float64),  # Umat
    np.ctypeslib.ndpointer(dtype=np.float64),  # Jmat
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)         # Nchi, Norb
]
_loader._lib.get_vmat_soc_orb.restype = None


def gen_SCmatrix(olist: np.ndarray, site: np.ndarray, U: float, J: float) -> tuple[np.ndarray, np.ndarray]:
    """
    @brief Generate spin (Smat) and charge (Cmat) vertex interaction matrices from uniform U and J parameters.
    @param  olist: Orbital-pair list of chi [Nchi, 2] int64 (1-based)
    @param   site: Site indices for each orbital pair [Nchi] int64
    @param      U: Hubbard Coulomb repulsion in eV
    @param      J: Hund's coupling in eV
    @return  Smat: Spin interaction matrix [Nchi, Nchi] float64
    @return  Cmat: Charge interaction matrix [Nchi, Nchi] float64
    """
    Nchi = len(olist)
    Smat = np.zeros((Nchi, Nchi), dtype=np.float64)
    Cmat = np.zeros((Nchi, Nchi), dtype=np.float64)
    _loader._lib.get_scmat(Smat, Cmat, olist, site, _loader.dbl(U), _loader.dbl(J), _loader.i64(Nchi))
    return Smat, Cmat


def gen_SCmatrix_orb(olist: np.ndarray, site: np.ndarray, Umat: np.ndarray,
                     Jmat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    @brief Generate spin (Smat) and charge (Cmat) interaction matrices using orbital-resolved U and J.
    @param  olist: Orbital-pair list of chi [Nchi, 2] int64 (1-based)
    @param   site: Site indices for each orbital pair [Nchi] int64
    @param   Umat: Orbital-resolved Coulomb matrix [Norb, Norb] float64
    @param   Jmat: Orbital-resolved Hund's coupling matrix [Norb, Norb] float64
    @return  Smat: Spin interaction matrix [Nchi, Nchi] float64
    @return  Cmat: Charge interaction matrix [Nchi, Nchi] float64
    """
    Nchi = len(olist)
    Norb = len(Umat)
    Smat = np.zeros((Nchi, Nchi), dtype=np.float64)
    Cmat = np.zeros((Nchi, Nchi), dtype=np.float64)
    # Keep ctypes argument order identical to Fortran ABI: Smat, Cmat, olist, site, Umat, Jmat, Nchi, Norb.
    _loader._lib.get_scmat_orb(Smat, Cmat, olist, site, Umat, Jmat, _loader.i64(Nchi), _loader.i64(Norb))
    return Smat, Cmat


def gen_Vmatrix(olist: np.ndarray, slist: np.ndarray, site: np.ndarray, invs: np.ndarray,
                U: float, J: float) -> np.ndarray:
    """
    @brief Generate the SOC interaction matrix Vmat from uniform U and J parameters.
    @param  olist: Orbital-pair list of chi [Nchi, 2] int64 (1-based)
    @param  slist: Spin-orbit label for each orbital [Norb] int64
    @param   site: Site indices for each orbital pair [Nchi] int64
    @param   invs: Inverse spin index mapping [Norb] int64
    @param      U: Hubbard Coulomb repulsion in eV
    @param      J: Hund's coupling in eV
    @return  Vmat: SOC interaction matrix [Nchi, Nchi] float64
    """
    Nchi, Norb = len(olist), len(slist)
    Vmat = np.zeros((Nchi, Nchi), dtype=np.float64)
    _loader._lib.get_vmat_soc(Vmat, olist, slist, site, invs, _loader.dbl(U),
                              _loader.dbl(J), _loader.i64(Nchi), _loader.i64(Norb))
    return Vmat


def gen_Vmatrix_orb(olist: np.ndarray, slist: np.ndarray, site: np.ndarray, invs: np.ndarray,
                    Umat: np.ndarray, Jmat: np.ndarray) -> np.ndarray:
    """
    @brief Generate the SOC interaction matrix Vmat using orbital-resolved U and J parameters.
    @param  olist: Orbital-pair list of chi [Nchi, 2] int64 (1-based)
    @param  slist: Spin-orbit label for each orbital [Norb] int64
    @param   site: Site indices for each orbital pair [Nchi] int64
    @param   invs: Inverse spin index mapping [Norb] int64
    @param   Umat: Orbital-resolved Coulomb matrix [Norb, Norb] float64
    @param   Jmat: Orbital-resolved Hund's coupling matrix [Norb, Norb] float64
    @return  Vmat: SOC interaction matrix [Nchi, Nchi] float64
    """
    Nchi, Norb = len(olist), len(Umat)
    Vmat = np.zeros((Nchi, Nchi), dtype=np.float64)
    _loader._lib.get_vmat_soc_orb(Vmat, olist, slist, site, invs, Umat, Jmat,
                                  _loader.i64(Nchi), _loader.i64(Norb))
    return Vmat
