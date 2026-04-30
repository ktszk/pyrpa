from ctypes import *
import numpy as np
from ._loader import _lib

def omp_params() ->tuple[int,bool]:
    """
    @fn omp_params
    @brief Query the underlying Fortran library for OpenMP configuration.
    @retval num_threads: Maximum number of OpenMP threads the Fortran code will use
    @retval  valid_flag: True if OpenMP was initialized correctly
    """
    omp_num=c_int64()
    omp_check=c_bool()
    _lib.openmp_params.argtypes=[POINTER(c_int64),POINTER(c_bool)]
    _lib.openmp_params.restype=None
    _lib.openmp_params(byref(omp_num),byref(omp_check))
    return omp_num.value,omp_check.value

def gen_ham(klist: np.ndarray, ham_r: np.ndarray, rvec: np.ndarray,
            Ovl_r: np.ndarray | None = None) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    """
    @fn gen_ham
    @brief Compute the k-space Hamiltonian by Fourier-transforming real-space hopping blocks.
    @param  klist: k-point list in fractional coordinates [Nk, 3]
    @param  ham_r: Real-space Hamiltonian blocks [Nr, Norb, Norb] complex128
    @param   rvec: Real-space displacement vectors [Nr, 3]
    @param  Ovl_r: Optional overlap matrix blocks [Nr, Norb, Norb]; if given, also returns Ovlk
    @retval  hamk: k-space Hamiltonian [Nk, Norb, Norb] complex128
    @retval  Ovlk: k-space overlap matrix [Nk, Norb, Norb] complex128 (only when Ovl_r is given)
    """
    Nk, Nr = len(klist), len(rvec)
    assert ham_r.ndim == 3, "ham_r must be a 3‑D array"
    Norb = int(np.sqrt(ham_r.size / Nr))
    hamk = np.zeros((Nk, Norb, Norb), dtype=np.complex128)
    _lib.gen_ham.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.gen_ham.restype = None
    _lib.gen_ham(hamk, klist, ham_r, rvec,
                 byref(c_int64(Nk)), byref(c_int64(Nr)), byref(c_int64(Norb)))
    if Ovl_r is not None:
        Ovlk = np.zeros((Nk, Norb, Norb), dtype=np.complex128)
        _lib.gen_ham(Ovlk, klist, Ovl_r, rvec,
                     byref(c_int64(Nk)), byref(c_int64(Nr)), byref(c_int64(Norb)))
        return hamk, Ovlk
    else:
        return hamk

def get_eig(hamk: np.ndarray, Ovlk: np.ndarray | None = None,
            sw: bool = True) -> tuple[np.ndarray, np.ndarray] | np.ndarray:
    """
    @fn get_eig
    @brief Diagonalize the k-space Hamiltonian, optionally solving the generalized eigenvalue problem.
    @param  hamk: k-space Hamiltonian [Nk, Norb, Norb] complex128
    @param  Ovlk: Overlap matrices [Nk, Norb, Norb] complex128; if given, solves generalized EVP
    @param    sw: If True, return both eigenvalues and eigenvectors; if False, return only eigenvalues
    @retval  eig: Eigenvalues [Nk, Norb] float64
    @retval  uni: Eigenvectors [Nk, Norb, Norb] complex128 (only returned when sw=True)
    """
    Nk = len(hamk)
    Norb = int(np.sqrt(hamk.size / Nk))
    eig = np.zeros((Nk, Norb), dtype=np.float64)
    uni = np.zeros((Nk, Norb, Norb), dtype=np.complex128)

    if Ovlk is None:
        _lib.get_eig.argtypes = [
            np.ctypeslib.ndpointer(dtype=np.float64),
            np.ctypeslib.ndpointer(dtype=np.complex128),
            np.ctypeslib.ndpointer(dtype=np.complex128),
            POINTER(c_int64), POINTER(c_int64)
        ]
        _lib.get_eig.restype = None
        _lib.get_eig(eig, uni, hamk, byref(c_int64(Nk)), byref(c_int64(Norb)))
    else:
        _lib.get_eig_mlo.argtypes = [
            np.ctypeslib.ndpointer(dtype=np.float64),
            np.ctypeslib.ndpointer(dtype=np.complex128),
            np.ctypeslib.ndpointer(dtype=np.complex128),
            np.ctypeslib.ndpointer(dtype=np.complex128),
            POINTER(c_int64), POINTER(c_int64)
        ]
        _lib.get_eig_mlo.restype = c_void_p
        _lib.get_eig_mlo(eig, uni, hamk, Ovlk, byref(c_int64(Nk)), byref(c_int64(Norb)))

    if sw:
        return eig, uni
    else:
        return eig

def get_ffermi(eig: np.ndarray, mu: float, temp: float) -> np.ndarray:
    """
    @fn get_ffermi
    @brief Compute Fermi-Dirac occupation numbers from eigenvalues.
    @param    eig: Eigenvalues [Nk, Norb] float64
    @param     mu: Chemical potential in eV
    @param   temp: Temperature in eV
    @return ffermi: Fermi-Dirac occupation numbers [Nk, Norb] float64
    """
    Nk = len(eig)
    Norb = int(eig.size / Nk)
    ffermi = np.zeros((Nk, Norb), dtype=np.float64)
    _lib.get_ffermi.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        POINTER(c_double), POINTER(c_double),
        POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.get_ffermi.restype = c_void_p
    _lib.get_ffermi(ffermi, eig, byref(c_double(mu)), byref(c_double(temp)),
                    byref(c_int64(Nk)), byref(c_int64(Norb)))
    return ffermi

def get_vlm0(klist: np.ndarray, ham_r: np.ndarray, rvec: np.ndarray) -> np.ndarray:
    """
    @fn get_vlm0
    @brief Compute unrotated band velocity matrices in orbital basis from real-space Hamiltonian.
    @param  klist: k-point list [Nk, 3] float64
    @param  ham_r: Real-space Hamiltonian blocks [Nr, Norb, Norb] complex128
    @param   rvec: Real-space displacement vectors [Nr, 3] float64
    @return   vk0: Velocity matrix elements [Nk, Norb, Norb, 3] complex128
    """
    Nk, Nr = len(klist), len(rvec)
    Norb = int(np.sqrt(ham_r.size / Nr))
    vk = np.zeros((Nk, Norb, Norb, 3), dtype=np.complex128)
    _lib.get_vlm0.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.get_vlm0.restype = c_void_p
    _lib.get_vlm0(vk, klist, ham_r, rvec,
                  byref(c_int64(Nk)), byref(c_int64(Nr)), byref(c_int64(Norb)))
    return vk

def get_vk(vk0: np.ndarray, mrot: np.ndarray, uni: np.ndarray) -> np.ndarray:
    """
    @fn get_vk
    @brief Transform orbital-basis velocity matrices to band-basis velocity expectation values.
    @param   vk0: Unrotated velocity matrices [Nk, Norb, Norb, 3] complex128
    @param  mrot: Rotation matrix (lattice-to-Cartesian) [3, 3] float64
    @param   uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @return   vk: Band velocities [Nk, Norb, 3] float64
    """
    Nk = len(uni)
    Norb = int(np.sqrt(uni.size / Nk))
    vk = np.zeros((Nk, Norb, 3), dtype=np.float64)
    _lib.get_veloc.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.get_veloc.restype = c_void_p
    _lib.get_veloc(vk, vk0, mrot, uni, byref(c_int64(Nk)), byref(c_int64(Norb)))
    return vk

def get_vnm(vk0: np.ndarray, mrot: np.ndarray, uni: np.ndarray) -> np.ndarray:
    """
    @fn get_vnm
    @brief Compute inter-band velocity matrix elements v_{nm}(k) for all band pairs and Cartesian directions.
    @param   vk0: Unrotated velocity matrices [Nk, Norb, Norb, 3] complex128
    @param  mrot: Rotation matrix (lattice-to-Cartesian) [3, 3] float64
    @param   uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @return   vnm: Inter-band velocity matrix elements [Nk, Norb, Norb, 3] complex128
    """
    Nk = len(uni)
    Norb = int(np.sqrt(uni.size / Nk))
    vk = np.zeros((Nk, Norb, Norb, 3), dtype=np.complex128)
    _lib.get_vnm.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.get_vnm.restype = c_void_p
    _lib.get_vnm(vk, vk0, mrot, uni, byref(c_int64(Nk)), byref(c_int64(Norb)))
    return vk

def get_vnmk(klist: np.ndarray, ham_r: np.ndarray, rvec: np.ndarray,
             mrot: np.ndarray, uni: np.ndarray) -> np.ndarray:
    """
    @fn get_vnmk
    @brief Convenience wrapper: compute inter-band velocity matrix elements by composing get_vlm0 and get_vnm.
    @param  klist: k-point list [Nk, 3] float64
    @param  ham_r: Real-space Hamiltonian blocks [Nr, Norb, Norb] complex128
    @param   rvec: Real-space displacement vectors [Nr, 3] float64
    @param   mrot: Rotation matrix (lattice-to-Cartesian) [3, 3] float64
    @param    uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @return   vnm: Inter-band velocity matrix elements [Nk, Norb, Norb, 3] complex128
    """
    vk0 = get_vlm0(klist, ham_r, rvec)
    return get_vnm(vk0, mrot, uni)

def get_veloc(klist: np.ndarray, ham_r: np.ndarray, rvec: np.ndarray,
              mrot: np.ndarray, uni: np.ndarray) -> np.ndarray:
    """
    @fn get_veloc
    @brief Convenience wrapper: compute band-diagonal velocity expectation values by composing get_vlm0 and get_vk.
    @param  klist: k-point list [Nk, 3] float64
    @param  ham_r: Real-space Hamiltonian blocks [Nr, Norb, Norb] complex128
    @param   rvec: Real-space displacement vectors [Nr, 3] float64
    @param   mrot: Rotation matrix (lattice-to-Cartesian) [3, 3] float64
    @param    uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @return    vk: Band velocities [Nk, Norb, 3] float64
    """
    vk0 = get_vlm0(klist, ham_r, rvec)
    return get_vk(vk0, mrot, uni)

def get_imass0(klist: np.ndarray, ham_r: np.ndarray, rvec: np.ndarray) -> np.ndarray:
    """
    @fn get_imass0
    @brief Compute the unrotated inverse effective mass tensor in orbital basis.
    @param   klist: k-point list [Nk, 3] float64
    @param   ham_r: Real-space Hamiltonian blocks [Nr, Norb, Norb] complex128
    @param    rvec: Real-space displacement vectors [Nr, 3] float64
    @return imass0: Inverse mass tensor before band rotation [Nk, Norb, Norb, 3, 3] complex128
    """
    Nk, Nr = len(klist), len(rvec)
    Norb = int(np.sqrt(ham_r.size / Nr))
    imass0 = np.zeros((Nk, Norb, Norb, 3, 3), dtype=np.complex128)
    _lib.get_imass0.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.get_imass0.restype = c_void_p
    _lib.get_imass0(imass0, klist, ham_r, rvec,
                    byref(c_int64(Nk)), byref(c_int64(Nr)), byref(c_int64(Norb)))
    return imass0

def get_imassk(imass0: np.ndarray, mrot: np.ndarray, uni: np.ndarray) -> np.ndarray:
    """
    @fn get_imassk
    @brief Rotate the inverse effective mass tensor from orbital basis into band basis.
    @param  imass0: Unrotated inverse mass tensor [Nk, Norb, Norb, 3, 3] complex128 (output of get_imass0)
    @param    mrot: Rotation matrix (lattice-to-Cartesian) [3, 3] float64
    @param     uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @return  imass: Band-basis inverse effective mass tensor [Nk, Norb, 3, 3] float64
    """
    Nk = len(uni)
    Norb = int(np.sqrt(uni.size / Nk))
    imass = np.zeros((Nk, Norb, 3, 3), dtype=np.float64)
    _lib.get_imassk.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.get_imassk.restype = c_void_p
    _lib.get_imassk(imass, imass0, mrot, uni, byref(c_int64(Nk)), byref(c_int64(Norb)))
    return imass

def get_mass(klist: np.ndarray, ham_r: np.ndarray, rvec: np.ndarray, mrot: np.ndarray,
             uni: np.ndarray, sw_imass: bool = False) -> np.ndarray:
    """
    @fn get_mass
    @brief Compute the effective mass tensor or its inverse for all bands at each k-point.
    @param    klist: k-point list [Nk, 3] float64
    @param    ham_r: Real-space Hamiltonian blocks [Nr, Norb, Norb] complex128
    @param     rvec: Real-space displacement vectors [Nr, 3] float64
    @param     mrot: Rotation matrix (lattice-to-Cartesian) [3, 3] float64
    @param      uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param sw_imass: If True, return inverse mass (imass) [Nk, Norb, 3, 3]; if False, return mass [Nk, Norb, 3, 3]
    @return   mass: Effective mass or inverse mass tensor [Nk, Norb, 3, 3]
    """
    import scipy.linalg as sclin
    imass0 = get_imass0(klist, ham_r, rvec)
    imass = get_imassk(imass0, mrot, uni)
    if sw_imass:
        return imass
    else:
        return np.array([[sclin.inv(im) for im in imas] for imas in imass])

def get_qshift(klist: np.ndarray, qpoint: np.ndarray) -> np.ndarray:
    """
    @fn get_qshift
    @brief For each k-point, find the index of k+q within the k-mesh (forward momentum shift).
    @param   klist: k-point list [Nk, 3] float64
    @param  qpoint: Momentum transfer vector [3] float64
    @return qshift: Index of k+q for each k-point [Nk] int64
    """
    Nk = len(klist)
    qshift = np.zeros(Nk, dtype=np.int64)
    _lib.get_qshift_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        POINTER(c_int64),
    ]
    _lib.get_qshift_.restype = c_void_p
    _lib.get_qshift_(qpoint, klist, qshift, byref(c_int64(Nk)))
    return qshift

def get_iqshift(klist: np.ndarray, qpoint: np.ndarray) -> np.ndarray:
    """
    @fn get_iqshift
    @brief For each k-point, find the index of k-q within the k-mesh (backward momentum shift).
    @param   klist: k-point list [Nk, 3] float64
    @param  qpoint: Momentum transfer vector [3] float64
    @return qshift: Index of k-q for each k-point [Nk] int64
    """
    Nk = len(klist)
    qshift = np.zeros(Nk, dtype=np.int64)
    _lib.get_iqshift_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        POINTER(c_int64),
    ]
    _lib.get_iqshift_.restype = c_void_p
    _lib.get_iqshift_(qpoint, klist, qshift, byref(c_int64(Nk)))
    return qshift

def gen_irr_k_TRS(Nx: int, Ny: int, Nz: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    @fn gen_irr_k_TRS
    @brief Generate the irreducible k-point list for a time-reversal symmetric lattice.
    @param     Nx: Number of k-points along kx
    @param     Ny: Number of k-points along ky
    @param     Nz: Number of k-points along kz
    @retval  klist: Irreducible k-point coordinates [Nk, 3] float64
    @retval   kmap: Full k-grid to irreducible mapping [Nkall, 3] int64
    @retval invk_ft_list: Inverse k-point list for Fortran FFT routines [Nkall, 3] int64
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
    _lib.generate_irr_kpoint_inv.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.generate_irr_kpoint_inv.restype = None
    _lib.generate_irr_kpoint_inv(klist, kmap, invk_ft_list, byref(c_int64(Nk)),
                                 byref(c_int64(Nx)), byref(c_int64(Ny)), byref(c_int64(Nz)))
    return klist, kmap, invk_ft_list

def gen_kpoint_weight(invk, Nk):
    """
    @fn gen_kpoint_weight
    @brief Compute k-point weights for the irreducible BZ by counting degeneracy from the full k-grid mapping.
    @param  invk: Inverse k-point mapping [Nkall] int64
    @param    Nk: Number of irreducible k-points
    @return weight: k-point weights [Nk] float64
    """
    Nkall = len(invk)
    weight = np.zeros(Nk, dtype=np.float64)
    _lib.get_kweight.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.int64),
        POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.get_kweight.restype = None
    _lib.get_kweight(weight, invk, byref(c_int64(Nk)), byref(c_int64(Nkall)))
    return weight

def get_plist(rvec, ham_r):
    """
    @fn get_plist
    @brief Determine the parity eigenvalues of each orbital from the Hamiltonian's inversion symmetry properties.
    @param  rvec: Real-space displacement vectors [Nr, 3] float64
    @param ham_r: Real-space Hamiltonian blocks [Nr, Norb, Norb] complex128
    @return plist: Parity sign (+1 or -1) for each orbital [Norb] float64
    """
    Nr = len(rvec)
    Norb = int(np.sqrt(ham_r.size / Nr))
    Pmn = np.zeros((Norb, Norb), dtype=np.float64)
    _lib.get_parity_prop.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        POINTER(c_int64), POINTER(c_int64)
    ]
    _lib.get_parity_prop.restype = None
    _lib.get_parity_prop(Pmn, rvec, ham_r, byref(c_int64(Norb)), byref(c_int64(Nr)))
    return np.sign(Pmn[0, :])
