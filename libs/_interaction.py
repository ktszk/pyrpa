from ctypes import *
import numpy as np
from ._loader import _lib

def gen_SCmatrix(olist: np.ndarray, site: np.ndarray, U: float, J: float) -> tuple[np.ndarray, np.ndarray]:
    """
    @fn gen_SCmatrix
    @brief Generate spin (Smat) and charge (Cmat) vertex interaction matrices from uniform U and J parameters.
    @param  olist: Orbital index list [Nchi] int64
    @param   site: Site indices for each orbital pair [Nchi] int64
    @param      U: Hubbard Coulomb repulsion in eV
    @param      J: Hund's coupling in eV
    @retval  Smat: Spin interaction matrix [Nchi, Nchi] float64
    @retval  Cmat: Charge interaction matrix [Nchi, Nchi] float64
    """
    Nchi = len(olist)
    Smat = np.zeros((Nchi, Nchi), dtype=np.float64)
    Cmat = np.zeros((Nchi, Nchi), dtype=np.float64)
    _lib.get_scmat.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        POINTER(c_double), POINTER(c_double), POINTER(c_int64)
    ]
    _lib.get_scmat.restype = c_void_p
    _lib.get_scmat(Smat, Cmat, olist, site, byref(c_double(U)), byref(c_double(J)), byref(c_int64(Nchi)))
    return Smat, Cmat

def gen_SCmatrix_orb(olist: np.ndarray, site: np.ndarray, Umat: np.ndarray, Jmat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    @fn gen_SCmatrix_orb
    @brief Generate spin (Smat) and charge (Cmat) interaction matrices using orbital-resolved U and J.
    @param  olist: Orbital index list [Nchi] int64
    @param   site: Site indices for each orbital pair [Nchi] int64
    @param   Umat: Orbital-resolved Coulomb matrix [Norb, Norb] float64
    @param   Jmat: Orbital-resolved Hund's coupling matrix [Norb, Norb] float64
    @retval  Smat: Spin interaction matrix [Nchi, Nchi] float64
    @retval  Cmat: Charge interaction matrix [Nchi, Nchi] float64
    """
    Nchi = len(olist)
    Norb = len(Umat)
    Smat = np.zeros((Nchi, Nchi), dtype=np.float64)
    Cmat = np.zeros((Nchi, Nchi), dtype=np.float64)
    _lib.get_scmat_orb.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.get_scmat_orb.restype = c_void_p
    _lib.get_scmat_orb(Smat, Cmat, olist, site, Umat, Jmat, byref(c_int64(Nchi)), byref(c_int64(Norb)))
    return Smat, Cmat

def gen_Vmatrix(olist: np.ndarray, slist: np.ndarray, site: np.ndarray, invs: np.ndarray,
                U: float, J: float) -> np.ndarray:
    """
    @fn gen_Vmatrix
    @brief Generate the SOC interaction matrix Vmat from uniform U and J parameters.
    @param  olist: Orbital index list [Nchi] int64
    @param  slist: Spin-orbit label for each orbital [Norb] int64
    @param   site: Site indices for each orbital pair [Nchi] int64
    @param   invs: Inverse spin index mapping [Norb] int64
    @param      U: Hubbard Coulomb repulsion in eV
    @param      J: Hund's coupling in eV
    @return  Vmat: SOC interaction matrix [Nchi, Nchi] float64
    """
    Nchi, Norb = len(olist), len(slist)
    Vmat = np.zeros((Nchi, Nchi), dtype=np.float64)
    _lib.get_vmat_soc.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        POINTER(c_double), POINTER(c_double),
        POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.get_vmat_soc.restype = None
    _lib.get_vmat_soc(Vmat, olist, slist, site, invs, byref(c_double(U)),
                      byref(c_double(J)), byref(c_int64(Nchi)), byref(c_int64(Norb)))
    return Vmat

def gen_Vmatrix_orb(olist: np.ndarray, slist: np.ndarray, site: np.ndarray, invs: np.ndarray,
                    Umat: np.ndarray, Jmat: np.ndarray) -> np.ndarray:
    """
    @fn gen_Vmatrix_orb
    @brief Generate the SOC interaction matrix Vmat using orbital-resolved U and J parameters.
    @param  olist: Orbital index list [Nchi] int64
    @param  slist: Spin-orbit label for each orbital [Norb] int64
    @param   site: Site indices for each orbital pair [Nchi] int64
    @param   invs: Inverse spin index mapping [Norb] int64
    @param   Umat: Orbital-resolved Coulomb matrix [Norb, Norb] float64
    @param   Jmat: Orbital-resolved Hund's coupling matrix [Norb, Norb] float64
    @return  Vmat: SOC interaction matrix [Nchi, Nchi] float64
    """
    Nchi, Norb = len(olist), len(Umat)
    Vmat = np.zeros((Nchi, Nchi), dtype=np.float64)
    _lib.get_vmat_soc_orb.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.get_vmat_soc_orb.restype = None
    _lib.get_vmat_soc_orb(Vmat, olist, slist, site, invs, Umat, Jmat,
                          byref(c_int64(Nchi)), byref(c_int64(Norb)))
    return Vmat
