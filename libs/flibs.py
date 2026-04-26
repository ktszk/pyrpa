from ctypes import *
import numpy as np
# Load the compiled Fortran shared library
flibs=np.ctypeslib.load_library("libs/libfmod.so",".")
# Each wrapper below sets argtypes/restype before calling to keep ctypes type-safe.
# All array arguments are passed as raw C pointers (ndpointer); scalars use byref(c_typeX).

def omp_params() ->tuple[int,bool]:
    """
    @fn omp_params
    @brief Query the underlying Fortran library for OpenMP configuration.
    @retval num_threads: Maximum number of OpenMP threads the Fortran code will use
    @retval  valid_flag: True if OpenMP was initialized correctly
    """
    omp_num=c_int64()
    omp_check=c_bool()
    flibs.openmp_params.argtypes=[POINTER(c_int64),POINTER(c_bool)]
    flibs.openmp_params.restype=None
    flibs.openmp_params(byref(omp_num),byref(omp_check))
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
    hamk = np.zeros((Nk, Norb, Norb), dtype=np.complex128)  # output buffer (C-contiguous)
    flibs.gen_ham.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64)
    ]
    flibs.gen_ham.restype = None
    flibs.gen_ham(hamk, klist, ham_r, rvec,
                   byref(c_int64(Nk)), byref(c_int64(Nr)), byref(c_int64(Norb)))
    if Ovl_r is not None:
        Ovlk = np.zeros((Nk, Norb, Norb), dtype=np.complex128)
        flibs.gen_ham(Ovlk, klist, Ovl_r, rvec,
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
        flibs.get_eig.argtypes = [
            np.ctypeslib.ndpointer(dtype=np.float64),    # eig
            np.ctypeslib.ndpointer(dtype=np.complex128), # uni
            np.ctypeslib.ndpointer(dtype=np.complex128), # hamk
            POINTER(c_int64), POINTER(c_int64)           # Nk, Norb
        ]
        flibs.get_eig.restype = None
        flibs.get_eig(eig, uni, hamk, byref(c_int64(Nk)), byref(c_int64(Norb)))
    else:
        flibs.get_eig_mlo.argtypes = [
            np.ctypeslib.ndpointer(dtype=np.float64),    # eig
            np.ctypeslib.ndpointer(dtype=np.complex128), # uni
            np.ctypeslib.ndpointer(dtype=np.complex128), # hamk
            np.ctypeslib.ndpointer(dtype=np.complex128), # Ovlk
            POINTER(c_int64), POINTER(c_int64)           # Nk, Norb
        ]
        flibs.get_eig_mlo.restype = c_void_p
        flibs.get_eig_mlo(eig, uni, hamk, Ovlk, byref(c_int64(Nk)), byref(c_int64(Norb)))

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

    flibs.get_ffermi.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64), # ffermi
        np.ctypeslib.ndpointer(dtype=np.float64), # eig
        POINTER(c_double), POINTER(c_double),      # mu, temp
        POINTER(c_int64), POINTER(c_int64)        # Nk, Norb
    ]
    flibs.get_ffermi.restype = c_void_p
    flibs.get_ffermi(ffermi, eig, byref(c_double(mu)), byref(c_double(temp)),
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
    flibs.get_vlm0.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64)
    ]
    flibs.get_vlm0.restype = c_void_p
    flibs.get_vlm0(vk, klist, ham_r, rvec,
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
    flibs.get_veloc.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        POINTER(c_int64), POINTER(c_int64)
    ]
    flibs.get_veloc.restype = c_void_p
    flibs.get_veloc(vk, vk0, mrot, uni, byref(c_int64(Nk)), byref(c_int64(Norb)))
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
    flibs.get_vnm.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        POINTER(c_int64), POINTER(c_int64)
    ]
    flibs.get_vnm.restype = c_void_p
    flibs.get_vnm(vk, vk0, mrot, uni, byref(c_int64(Nk)), byref(c_int64(Norb)))
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
    vk = get_vnm(vk0, mrot, uni)
    return vk

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
    vk = get_vk(vk0, mrot, uni)
    return vk

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
    Nk_p = byref(c_int64(Nk))
    Nr_p = byref(c_int64(Nr))
    Norb_p = byref(c_int64(Norb))
    flibs.get_imass0.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        POINTER(c_int64), POINTER(c_int64),
        POINTER(c_int64)
    ]
    flibs.get_imass0.restype = c_void_p
    flibs.get_imass0(imass0, klist, ham_r, rvec, Nk_p, Nr_p, Norb_p)
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
    flibs.get_imassk.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        POINTER(c_int64), POINTER(c_int64)
    ]
    flibs.get_imassk.restype = c_void_p
    flibs.get_imassk(imass, imass0, mrot, uni, byref(c_int64(Nk)), byref(c_int64(Norb)))
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
        # invert each 3x3 block separately
        return np.array([[sclin.inv(im) for im in imas] for imas in imass])

def gen_Green0(eig: np.ndarray, uni: np.ndarray,
               mu: float, temp: float, Nw: int) -> np.ndarray:
    """
    @fn gen_Green0
    @brief Construct the non-interacting Matsubara Green's function G0(k, iω_n).
    @param   eig: Eigenvalues [Nk, Norb] float64
    @param   uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param    mu: Chemical potential in eV
    @param  temp: Temperature in eV
    @param    Nw: Number of Matsubara frequency points
    @return   Gk: Non-interacting Green's function [Norb, Norb, Nw, Nk] complex128
    """
    Nk = len(eig)
    Norb = int(eig.size / Nk)
    Gk = np.zeros((Norb, Norb, Nw, Nk), dtype=np.complex128)
    flibs.gen_green0_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),        # Gk
        np.ctypeslib.ndpointer(dtype=np.float64),           # eig
        np.ctypeslib.ndpointer(dtype=np.complex128),        # uni
        POINTER(c_double), POINTER(c_double),                # mu, temp
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64)  # Nk, Nw, Norb
    ]
    flibs.gen_green0_.restype = c_void_p
    flibs.gen_green0_(Gk, eig, uni, byref(c_double(mu)), byref(c_double(temp)),
                      byref(c_int64(Nk)), byref(c_int64(Nw)), byref(c_int64(Norb)))
    return Gk

def gen_green(selfen: np.ndarray, hamk: np.ndarray,
              mu: float, temp: float) -> np.ndarray:
    """
    @fn gen_green
    @brief Construct the interacting Matsubara Green's function G(k, iω_n) = [iω_n + mu - H(k) - Sigma(k,iω_n)]^-1.
    @param  selfen: Self-energy tensor [Norb, Norb, Nw, Nk] complex128
    @param    hamk: k-space Hamiltonian [Nk, Norb, Norb] complex128
    @param      mu: Chemical potential in eV
    @param    temp: Temperature in eV
    @return     Gk: Interacting Green's function [Norb, Norb, Nw, Nk] complex128
    """
    Nk = len(hamk)
    Norb = int(np.sqrt(hamk.size / Nk))
    Nw = int(selfen.size / (Nk * Norb * Norb))
    Gk = np.zeros((Norb, Norb, Nw, Nk), dtype=np.complex128)
    flibs.gen_green_inv_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),        # Gk
        np.ctypeslib.ndpointer(dtype=np.complex128),        # selfen
        np.ctypeslib.ndpointer(dtype=np.complex128),        # hamk
        POINTER(c_double), POINTER(c_double),                # mu, temp
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64)  # Nk, Nw, Norb
    ]
    flibs.gen_green_inv_.restype = c_void_p
    flibs.gen_green_inv_(Gk, selfen, hamk, byref(c_double(mu)), byref(c_double(temp)),
                         byref(c_int64(Nk)), byref(c_int64(Nw)), byref(c_int64(Norb)))
    flibs.getinv_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),        # Gk
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64)  # Nk, Nw, Norb
    ]
    flibs.getinv_.restype = c_void_p
    flibs.getinv_(Gk, byref(c_int64(Nk)), byref(c_int64(Nw)), byref(c_int64(Norb)))
    return Gk

def gen_green_from_eig(selfen: np.ndarray, eig: np.ndarray,
                       uni: np.ndarray, mu: float, temp: float) -> np.ndarray:
    """
    @fn gen_green_from_eig
    @brief Construct the interacting Green's function using pre-computed eigenvalues and eigenvectors.
    @param  selfen: Self-energy tensor [Norb, Norb, Nw, Nk] complex128
    @param     eig: Eigenvalues [Nk, Norb] float64
    @param     uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param      mu: Chemical potential in eV
    @param    temp: Temperature in eV
    @return     Gk: Interacting Green's function [Norb, Norb, Nw, Nk] complex128
    """
    Nk = len(eig)
    Norb = int(eig.size / Nk)
    Nw = int(selfen.size / (Nk * Norb * Norb))
    Gk = np.zeros((Norb, Norb, Nw, Nk), dtype=np.complex128)
    flibs.gen_green_inv_from_eig.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),        # Gk
        np.ctypeslib.ndpointer(dtype=np.complex128),        # selfen
        np.ctypeslib.ndpointer(dtype=np.complex128),        # uni
        np.ctypeslib.ndpointer(dtype=np.float64),           # eig
        POINTER(c_double), POINTER(c_double),                # mu, temp
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64)  # Nk, Nw, Norb
    ]
    flibs.gen_green_inv_from_eig.restype = c_void_p
    flibs.gen_green_inv_from_eig(Gk, selfen, uni, eig, byref(c_double(mu)), byref(c_double(temp)),
                                 byref(c_int64(Nk)), byref(c_int64(Nw)), byref(c_int64(Norb)))
    flibs.getinv.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),        # Gk
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64)  # Nk, Nw, Norb
    ]
    flibs.getinv.restype = c_void_p
    flibs.getinv(Gk, byref(c_int64(Nk)), byref(c_int64(Nw)), byref(c_int64(Norb)))
    return Gk

def gen_tr_Greenw_0(eig: np.ndarray, mu: float,
                    wlist: np.ndarray, delta: float) -> np.ndarray:
    """
    @fn gen_tr_Greenw_0
    @brief Compute the k-resolved spectral function (imaginary part of trace of G0) on the real-frequency axis.
    @param    eig: Eigenvalues [Nk, Norb] float64
    @param     mu: Chemical potential in eV
    @param  wlist: Real-frequency mesh [Nw] float64
    @param  delta: Lorentzian broadening parameter in eV
    @return  trGk: -Im[Tr G0(k,w+i*delta)] [Nk, Nw] float64
    """
    Nk, Nw = len(eig), len(wlist)
    Norb = int(eig.size / Nk)
    trGk = np.zeros((Nk, Nw), dtype=np.complex128)
    flibs.gen_tr_greenw_0.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),  # trGk
        np.ctypeslib.ndpointer(dtype=np.float64),    # wlist
        np.ctypeslib.ndpointer(dtype=np.float64),    # eig
        POINTER(c_double), POINTER(c_double),        # mu, delta
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64)
    ]
    flibs.gen_tr_greenw_0.restype = c_void_p
    flibs.gen_tr_greenw_0(trGk, wlist, eig, byref(c_double(mu)), byref(c_double(delta)),
                          byref(c_int64(Nk)), byref(c_int64(Nw)), byref(c_int64(Norb)))
    return -trGk.imag

def gen_dos(eig: np.ndarray, uni: np.ndarray, mu: float,
            wlist: np.ndarray, delta: float) -> np.ndarray:
    """
    @fn gen_dos
    @brief Compute the orbital-resolved partial density of states (pDOS) on the real-frequency axis.
    @param    eig: Eigenvalues [Nk, Norb] float64
    @param    uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param     mu: Chemical potential in eV
    @param  wlist: Real-frequency mesh [Nw] float64
    @param  delta: Lorentzian broadening parameter in eV
    @return  pDos: Orbital-resolved pDOS -Im[G(orb,w)] [Norb, Nw] float64
    """
    Nk, Nw = len(eig), len(wlist)
    Norb = int(eig.size / Nk)
    pDos = np.zeros((Norb, Nw), dtype=np.complex128)
    flibs.gen_dos.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),  # Dos
        np.ctypeslib.ndpointer(dtype=np.float64),     # wlist
        np.ctypeslib.ndpointer(dtype=np.float64),     # eig
        np.ctypeslib.ndpointer(dtype=np.complex128),  # uni
        POINTER(c_double), POINTER(c_double),         # mu, delta
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64)
    ]
    flibs.gen_dos.restype = c_void_p
    flibs.gen_dos(pDos, wlist, eig, uni, byref(c_double(mu)), byref(c_double(delta)),
                  byref(c_int64(Nk)), byref(c_int64(Nw)), byref(c_int64(Norb)))
    return -pDos.imag

def get_chi0_conv(Gk: np.ndarray, kmap: np.ndarray, invk: np.ndarray, olist: np.ndarray,
                  temp: float, Nx: int, Ny: int, Nz: int) -> np.ndarray:
    """
    @fn get_chi0_conv
    @brief Compute the irreducible susceptibility chi0 by convolution of Green's functions over the full BZ.
    @param     Gk: Green's function tensor [Norb, Norb, Nw, Nk] complex128
    @param   kmap: Mapping from full k-mesh indices to reduced mesh [Nkall] int64
    @param   invk: Inverse k-point mapping [Nk] int64
    @param  olist: Orbital indices for susceptibility calculation [Nchi] int64
    @param   temp: Temperature in eV
    @param     Nx: Number of k-points along kx
    @param     Ny: Number of k-points along ky
    @param     Nz: Number of k-points along kz
    @return   chi: Irreducible susceptibility tensor [Nchi, Nchi, Nw, Nk] complex128
    """
    Nkall, Nk = len(kmap), len(Gk[0, 0, 0])
    Norb, Nchi = len(Gk), len(olist)
    Nw = int(Gk.size / (Norb * Norb * Nk))
    chi = np.zeros((Nchi, Nchi, Nw, Nk), dtype=np.complex128)
    flibs.get_chi0_conv_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),         # chi
        np.ctypeslib.ndpointer(dtype=np.complex128),         # Gk
        np.ctypeslib.ndpointer(dtype=np.int64),              # kmap
        np.ctypeslib.ndpointer(dtype=np.int64),              # invk
        np.ctypeslib.ndpointer(dtype=np.int64),              # olist
        POINTER(c_double), POINTER(c_int64), POINTER(c_int64),  # temp, Nx, Ny
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64),  # Nz, Nw, Nk
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64),  # Nkall, Norb, Nchi
    ]
    flibs.get_chi0_conv_.restype = c_void_p
    flibs.get_chi0_conv_(chi, Gk, kmap, invk, olist, byref(c_double(temp)), byref(c_int64(Nx)),
                         byref(c_int64(Ny)), byref(c_int64(Nz)), byref(c_int64(Nw)), byref(c_int64(Nk)),
                         byref(c_int64(Nkall)), byref(c_int64(Norb)), byref(c_int64(Nchi)))
    return chi

def get_chi0_sum(Gk: np.ndarray, invk: np.ndarray, klist: np.ndarray,
                 olist: np.ndarray, temp: float) -> np.ndarray:
    """
    @fn get_chi0_sum
    @brief Compute the irreducible susceptibility chi0 by direct BZ summation (alternative to convolution).
    @param     Gk: Green's function tensor [Norb, Norb, Nw, Nk] complex128
    @param   invk: Inverse k-point mapping [Nkall] int64
    @param  klist: k-point list [Nk, 3] float64
    @param  olist: Orbital indices for susceptibility calculation [Nchi] int64
    @param   temp: Temperature in eV
    @return   chi: Irreducible susceptibility tensor [Nchi, Nchi, Nw, Nk] complex128
    """
    Nkall, Nk = len(invk), len(klist)
    Norb, Nchi = len(Gk), len(olist)
    Nw = int(Gk.size / (Norb * Norb * Nk))
    chi = np.zeros((Nchi, Nchi, Nw, Nk), dtype=np.complex128)
    flibs.get_chi0_sum.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),   # chi
        np.ctypeslib.ndpointer(dtype=np.complex128),   # Gk
        np.ctypeslib.ndpointer(dtype=np.float64),      # klist
        np.ctypeslib.ndpointer(dtype=np.int64),        # olist
        POINTER(c_double), POINTER(c_int64),           # temp, Nw
        POINTER(c_int64), POINTER(c_int64),            # Nk, Nkall
        POINTER(c_int64), POINTER(c_int64)             # Norb, Nchi
    ]
    flibs.get_chi0_sum.restype = c_void_p
    flibs.get_chi0_sum(chi, Gk, klist, olist, byref(c_double(temp)), byref(c_int64(Nw)),
                       byref(c_int64(Nk)), byref(c_int64(Nkall)), byref(c_int64(Norb)),
                       byref(c_int64(Nchi)))
    return chi

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
    flibs.get_vsigma_flex_nosoc_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64)
    ]
    flibs.get_vsigma_flex_nosoc_.restype = c_void_p
    flibs.get_vsigma_flex_nosoc_(chi, Smat, Cmat, byref(c_int64(Nk)),
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
    flibs.mkself.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),          # sigmak
        POINTER(c_double),                                    # muself
        np.ctypeslib.ndpointer(dtype=np.float64),             # Smat
        np.ctypeslib.ndpointer(dtype=np.float64),             # Cmat
        np.ctypeslib.ndpointer(dtype=np.int64),               # kmap
        np.ctypeslib.ndpointer(dtype=np.int64),               # invk
        np.ctypeslib.ndpointer(dtype=np.int64),               # olist
        np.ctypeslib.ndpointer(dtype=np.complex128),          # hamk
        np.ctypeslib.ndpointer(dtype=np.float64),             # eig
        np.ctypeslib.ndpointer(dtype=np.complex128),          # uni
        POINTER(c_double), POINTER(c_double),                  # mu, fill
        POINTER(c_double), POINTER(c_int64),                   # temp, scf_loop
        POINTER(c_double), POINTER(c_double), POINTER(c_int64),  # pp, eps, Nkall
        POINTER(c_int64), POINTER(c_int64),                    # Nk, Nw
        POINTER(c_int64), POINTER(c_int64),                    # Nchi, Norb
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64),   # Nx, Ny, Nz
        POINTER(c_int64), POINTER(c_bool), POINTER(c_bool),      # sub_sigma, sw_out, sw_in
        POINTER(c_int64),                                        # m_diis
        POINTER(c_bool)                                          # sw_rescale
    ]
    flibs.mkself.restype = c_void_p
    flibs.mkself(sigmak, byref(mu_self), Smat, Cmat, kmap, invk, olist, hamk, eig, uni,
        byref(c_double(mu)), byref(c_double(fill)), byref(c_double(temp)), byref(c_int64(scf_loop)),
        byref(c_double(pp)), byref(c_double(eps)), byref(c_int64(Nkall)), byref(c_int64(Nk)),
        byref(c_int64(Nw)), byref(c_int64(Norb)), byref(c_int64(Nchi)), byref(c_int64(Nx)),
        byref(c_int64(Ny)), byref(c_int64(Nz)), byref(c_int64(sub_sigma)), byref(c_bool(sw_out)),
        byref(c_bool(sw_in)), byref(c_int64(m_diis)), byref(c_bool(sw_rescale)))
    return sigmak, mu_self.value

def mkself_soc(Vmat: np.ndarray, kmap: np.ndarray, invk: np.ndarray, invs: np.ndarray,
           olist: np.ndarray, slist: np.ndarray, hamk: np.ndarray, eig: np.ndarray, uni: np.ndarray, mu: float,
           fill: float, temp: float, Nw: int, Nx: int, Ny: int, Nz: int, sw_out: bool,
           sw_in: bool, sub_sigma: int = 1, scf_loop: int = 300, eps: float = 1.0e-4,
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
    flibs.mkself_soc.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),            # sigmak
        POINTER(c_double),                                      # muself
        np.ctypeslib.ndpointer(dtype=np.float64),               # Vmat
        np.ctypeslib.ndpointer(dtype=np.int64),                 # kmap
        np.ctypeslib.ndpointer(dtype=np.int64),                 # invk
        np.ctypeslib.ndpointer(dtype=np.int64),                 # invs
        np.ctypeslib.ndpointer(dtype=np.int64),                 # olist
        np.ctypeslib.ndpointer(dtype=np.int64),                 # slist
        np.ctypeslib.ndpointer(dtype=np.complex128),            # hamk
        np.ctypeslib.ndpointer(dtype=np.float64),               # eig
        np.ctypeslib.ndpointer(dtype=np.complex128),            # uni
        POINTER(c_double), POINTER(c_double),                   # mu, fill
        POINTER(c_double), POINTER(c_int64),                    # temp, scf_loop
        POINTER(c_double), POINTER(c_double), POINTER(c_int64), # pp, eps, Nkall
        POINTER(c_int64), POINTER(c_int64),                     # Nk, Nw
        POINTER(c_int64), POINTER(c_int64),                     # Nchi, Norb
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64),   # Nx, Ny, Nz
        POINTER(c_int64), POINTER(c_bool), POINTER(c_bool),       # sub_sigma, sw_out, sw_in
        POINTER(c_int64), POINTER(c_bool)                        # m_diis, sw_rescale
    ]
    flibs.mkself_soc(sigmak,byref(mu_self),Vmat,kmap,invk,invs,olist,slist,hamk,eig,uni,
               byref(c_double(mu)), byref(c_double(fill)), byref(c_double(temp)), byref(c_int64(scf_loop)),
               byref(c_double(pp)), byref(c_double(eps)), byref(c_int64(Nkall)), byref(c_int64(Nk)),
               byref(c_int64(Nw)), byref(c_int64(Norb)), byref(c_int64(Nchi)), byref(c_int64(Nx)),
               byref(c_int64(Ny)), byref(c_int64(Nz)), byref(c_int64(sub_sigma)), byref(c_bool(sw_out)), byref(c_bool(sw_in)), byref(c_int64(m_diis)), byref(c_bool(sw_rescale)))
    return sigmak, mu_self.value

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
    flibs.get_qshift_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),  # qpoint
        np.ctypeslib.ndpointer(dtype=np.float64),  # klist
        np.ctypeslib.ndpointer(dtype=np.int64),    # qshift
        POINTER(c_int64),                          # Nk
    ]
    flibs.get_qshift_.restype = c_void_p
    flibs.get_qshift_(qpoint, klist, qshift, byref(c_int64(Nk)))
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
    flibs.get_iqshift_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),  # qpoint
        np.ctypeslib.ndpointer(dtype=np.float64),  # klist
        np.ctypeslib.ndpointer(dtype=np.int64),    # qshift
        POINTER(c_int64),                          # Nk
    ]
    flibs.get_iqshift_.restype = c_void_p
    flibs.get_iqshift_(qpoint, klist, qshift, byref(c_int64(Nk)))
    return qshift

def get_chi_irr(uni: np.ndarray, eig: np.ndarray, ffermi: np.ndarray, qshift: np.ndarray,
                olist: np.ndarray, wlist: np.ndarray, idelta: float, temp: float) -> np.ndarray:
    """
    @fn get_chi_irr
    @brief Compute the irreducible (bare bubble) susceptibility chi0(q, w) for a given momentum shift.
    @param     uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param     eig: Eigenvalues [Nk, Norb] float64
    @param  ffermi: Fermi-Dirac occupations [Nk, Norb] float64
    @param  qshift: Index of k+q for each k-point [Nk] int64 (from get_qshift)
    @param   olist: Orbital index list [Nchi] int64
    @param   wlist: Real-frequency mesh [Nw] float64
    @param  idelta: Lorentzian broadening parameter in eV
    @param    temp: Temperature in eV
    @return   chi0: Irreducible susceptibility [Nw, Nchi, Nchi] complex128
    """
    Nk, Nw = len(eig), len(wlist)
    Norb, Nchi = int(eig.size / Nk), len(olist)
    chi = np.zeros((Nw, Nchi, Nchi), dtype=np.complex128)
    eps = idelta * 1e-3
    flibs.get_chi_irr.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128), # chi
        np.ctypeslib.ndpointer(dtype=np.complex128), # uni
        np.ctypeslib.ndpointer(dtype=np.float64),    # eig
        np.ctypeslib.ndpointer(dtype=np.float64),    # ffermi
        np.ctypeslib.ndpointer(dtype=np.int64),      # qshift
        np.ctypeslib.ndpointer(dtype=np.int64),      # olist
        np.ctypeslib.ndpointer(dtype=np.float64),    # wlist
        POINTER(c_int64), POINTER(c_int64),           # Nchi, Norb
        POINTER(c_int64), POINTER(c_int64),           # Nk, Nw
        POINTER(c_double),                            # idelta
        POINTER(c_double), POINTER(c_double),         # eps, temp
    ]
    flibs.get_chi_irr.restype = c_void_p
    flibs.get_chi_irr(chi, uni, eig, ffermi, qshift, olist, wlist, byref(c_int64(Nchi)),
                      byref(c_int64(Norb)), byref(c_int64(Nk)), byref(c_int64(Nw)),
                      byref(c_double(idelta)), byref(c_double(eps)), byref(c_double(temp)))
    return chi

def chis_qmap(uni: np.ndarray, eig: np.ndarray, ffermi: np.ndarray, klist: np.ndarray,
              Smat: np.ndarray, olist: np.ndarray, Nx: int, Ny: int, temp: float,
              ecut: float, idelta: float) -> tuple[np.ndarray, np.ndarray]:
    """
    @fn chis_qmap
    @brief Compute spin susceptibility and bare susceptibility maps on the qx-qy plane.
    @param     uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param     eig: Eigenvalues [Nk, Norb] float64
    @param  ffermi: Fermi-Dirac occupations [Nk, Norb] float64
    @param   klist: k-point list [Nk, 3] float64
    @param    Smat: Spin interaction matrix [Nchi, Nchi] float64
    @param   olist: Orbital index list [Nchi] int64
    @param      Nx: Number of q-points along qx
    @param      Ny: Number of q-points along qy
    @param    temp: Temperature in eV
    @param    ecut: Energy cutoff for the BZ summation in eV
    @param  idelta: Lorentzian broadening parameter in eV
    @retval   chis: Spin susceptibility map [Nx, Ny] complex128
    @retval    chi: Irreducible susceptibility map [Nx, Ny] complex128
    """
    Nk = len(eig)
    Norb, Nchi = int(eig.size / Nk), len(olist)
    chi = np.zeros((Nx, Ny), dtype=np.complex128)
    chis = np.zeros((Nx, Ny), dtype=np.complex128)
    eps = idelta * 1e-3
    flibs.chiq_map.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),        # chis
        np.ctypeslib.ndpointer(dtype=np.complex128),        # chi
        np.ctypeslib.ndpointer(dtype=np.complex128),        # uni
        np.ctypeslib.ndpointer(dtype=np.float64),           # eig
        np.ctypeslib.ndpointer(dtype=np.float64),           # ffermi
        np.ctypeslib.ndpointer(dtype=np.float64),           # klist
        np.ctypeslib.ndpointer(dtype=np.float64),           # Smat
        np.ctypeslib.ndpointer(dtype=np.int64),             # olist
        POINTER(c_double), POINTER(c_double),                # temp, ecut
        POINTER(c_double), POINTER(c_double),                # idelta, eps
        POINTER(c_int64), POINTER(c_int64),                  # Nx, Ny
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64) # Nk, Norb, Nchi
    ]
    flibs.chiq_map.restype = c_void_p
    flibs.chiq_map(chis, chi, uni, eig, ffermi, klist, Smat, olist, byref(c_double(temp)),
                   byref(c_double(ecut)), byref(c_double(idelta)), byref(c_double(eps)),
                   byref(c_int64(Nx)), byref(c_int64(Ny)), byref(c_int64(Nk)),
                   byref(c_int64(Norb)), byref(c_int64(Nchi)))
    return chis, chi

def get_tr_chi(chis: np.ndarray, chi0: np.ndarray,
               olist: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    @fn get_tr_chi
    @brief Compute orbital-resolved traces of the interacting and bare susceptibilities.
    @param    chis: Interacting spin susceptibility [Nw, Nchi, Nchi] complex128
    @param    chi0: Irreducible susceptibility [Nw, Nchi, Nchi] complex128
    @param   olist: Orbital index list [Nchi] int64
    @retval trchis: Trace of interacting chi_s [Nw] complex128
    @retval trchi0: Trace of bare chi0 [Nw] complex128
    @retval chis_orb: Orbital-indexed chi_s [Nw, Norb+2] complex128
    """
    Nchi, Nw, Norb = len(olist), len(chi0), olist.max()
    trchis = np.zeros(Nw, dtype=np.complex128)
    trchi0 = np.zeros(Nw, dtype=np.complex128)
    chis_orb = np.zeros((Nw, Norb + 2), dtype=np.complex128)
    flibs.get_tr_chi.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128), # trchis
        np.ctypeslib.ndpointer(dtype=np.complex128), # trchi0
        np.ctypeslib.ndpointer(dtype=np.complex128), # chis_orb
        np.ctypeslib.ndpointer(dtype=np.complex128), # chis
        np.ctypeslib.ndpointer(dtype=np.complex128), # chi0
        np.ctypeslib.ndpointer(dtype=np.int64),      # olist
        POINTER(c_int64), POINTER(c_int64),           # Nw, Nchi
        POINTER(c_int64)                              # Norb
    ]
    flibs.get_tr_chi.restype = c_void_p
    flibs.get_tr_chi(trchis, trchi0, chis_orb, chis, chi0, olist,
                     byref(c_int64(Nw)), byref(c_int64(Nchi)), byref(c_int64(Norb)))
    return trchis, trchi0, chis_orb

def get_phi_irr(uni: np.ndarray, eig: np.ndarray, ffermi: np.ndarray, qshift: np.ndarray,
                olist: np.ndarray, wlist: np.ndarray, idelta: float, mu: float,
                temp: float) -> np.ndarray:
    """
    @fn get_phi_irr
    @brief Compute the irreducible pairing susceptibility phi0(q, w) in the particle-particle channel.
    @param     uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param     eig: Eigenvalues [Nk, Norb] float64
    @param  ffermi: Fermi-Dirac occupations [Nk, Norb] float64
    @param  qshift: Index of k+q for each k-point [Nk] int64 (from get_qshift)
    @param   olist: Orbital index list [Nchi] int64
    @param   wlist: Frequency mesh [Nw] float64
    @param  idelta: Lorentzian broadening parameter in eV
    @param      mu: Chemical potential in eV
    @param    temp: Temperature in eV
    @return   phi0: Irreducible pairing susceptibility [Nw, Nchi, Nchi] complex128
    """
    Nk, Nw = len(eig), len(wlist)
    Norb, Nchi = int(eig.size / Nk), len(olist)
    phi = np.zeros((Nw, Nchi, Nchi), dtype=np.complex128)
    eps = idelta * 1e-3
    flibs.get_phi_irr.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128), # phi
        np.ctypeslib.ndpointer(dtype=np.complex128), # uni
        np.ctypeslib.ndpointer(dtype=np.float64),    # eig
        np.ctypeslib.ndpointer(dtype=np.float64),    # ffermi
        np.ctypeslib.ndpointer(dtype=np.int64),      # qshift
        np.ctypeslib.ndpointer(dtype=np.int64),      # olist
        np.ctypeslib.ndpointer(dtype=np.float64),    # wlist
        POINTER(c_int64), POINTER(c_int64),           # Nchi, Norb
        POINTER(c_int64), POINTER(c_int64),           # Nk, Nw
        POINTER(c_double), POINTER(c_double),         # eps, idelta
        POINTER(c_double), POINTER(c_double),         # mu, temp
    ]
    flibs.get_phi_irr.restype = c_void_p
    flibs.get_phi_irr(phi, uni, eig, ffermi, qshift, olist, wlist, byref(c_int64(Nchi)),
                      byref(c_int64(Norb)), byref(c_int64(Nk)), byref(c_int64(Nw)),
                      byref(c_double(idelta)), byref(c_double(eps)),
                      byref(c_double(mu)), byref(c_double(temp)))
    return phi

def phi_qmap(uni: np.ndarray, eig: np.ndarray, ffermi: np.ndarray, klist: np.ndarray,
             olist: np.ndarray, Nx: int, Ny: int, mu: float, temp: float, ecut: float,
             idelta: float, sw_omega: bool) -> np.ndarray:
    """
    @fn phi_qmap
    @brief Compute the pairing susceptibility map on the qx-qy plane at a given energy cutoff.
    @param      uni: Eigenvector matrices [Nk, Norb, Norb] complex128
    @param      eig: Eigenvalues [Nk, Norb] float64
    @param   ffermi: Fermi-Dirac occupations [Nk, Norb] float64
    @param    klist: k-point list [Nk, 3] float64
    @param    olist: Orbital index list [Nchi] int64
    @param       Nx: Number of q-points along qx
    @param       Ny: Number of q-points along qy
    @param       mu: Chemical potential in eV
    @param     temp: Temperature in eV
    @param     ecut: Energy cutoff for the BZ summation in eV
    @param   idelta: Lorentzian broadening parameter in eV
    @param sw_omega: If True, include Matsubara frequency summation; if False, use static limit
    @return     phi: Pairing susceptibility map [Nx, Ny] complex128
    """
    Nk = len(eig)
    Norb, Nchi = int(eig.size / Nk), len(olist)
    phi = np.zeros((Nx, Ny), dtype=np.complex128)
    eps = idelta * 1e-3
    flibs.phiq_map.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),        # phi
        np.ctypeslib.ndpointer(dtype=np.complex128),        # uni
        np.ctypeslib.ndpointer(dtype=np.float64),           # eig
        np.ctypeslib.ndpointer(dtype=np.float64),           # ffermi
        np.ctypeslib.ndpointer(dtype=np.float64),           # klist
        np.ctypeslib.ndpointer(dtype=np.int64),             # olist
        POINTER(c_double), POINTER(c_double),                # mu, temp
        POINTER(c_double), POINTER(c_double),                # ecut, idelta
        POINTER(c_double),                                   # eps
        POINTER(c_int64), POINTER(c_int64),                  # Nx, Ny
        POINTER(c_int64), POINTER(c_int64),                  # Nk, Norb
        POINTER(c_int64), POINTER(c_bool)                    # Nchi, sw_omega
    ]
    flibs.phiq_map.restype = c_void_p
    flibs.phiq_map(phi, uni, eig, ffermi, klist, olist, byref(c_double(mu)), byref(c_double(temp)),
                   byref(c_double(ecut)), byref(c_double(idelta)), byref(c_double(eps)),
                   byref(c_int64(Nx)), byref(c_int64(Ny)), byref(c_int64(Nk)),
                   byref(c_int64(Norb)), byref(c_int64(Nchi)), byref(c_bool(sw_omega)))
    return phi

def get_tr_phi(phi: np.ndarray, olist: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    @fn get_tr_phi
    @brief Compute orbital-resolved trace of the pairing susceptibility phi.
    @param     phi: Pairing susceptibility [Nw, Nchi, Nchi] complex128
    @param   olist: Orbital index list [Nchi] int64
    @retval  trphi: Trace of phi [Nw] complex128
    @retval phi_orb: Orbital-indexed phi [Nw, Norb+2] complex128
    """
    Nchi, Nw, Norb = len(olist), len(phi), olist.max()
    trphi = np.zeros(Nw, dtype=np.complex128)
    phi_orb = np.zeros((Nw, Norb + 2), dtype=np.complex128)
    flibs.get_tr_phi.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128), # trphi
        np.ctypeslib.ndpointer(dtype=np.complex128), # phi_orb
        np.ctypeslib.ndpointer(dtype=np.complex128), # phi
        np.ctypeslib.ndpointer(dtype=np.int64),      # olist
        POINTER(c_int64), POINTER(c_int64),           # Nw, Nchi
        POINTER(c_int64)                              # Norb
    ]
    flibs.get_tr_phi.restype = c_void_p
    flibs.get_tr_phi(trphi, phi_orb, phi, olist, byref(c_int64(Nw)),
                     byref(c_int64(Nchi)), byref(c_int64(Norb)))
    return trphi, phi_orb

def get_chis(chi0: np.ndarray, Smat: np.ndarray) -> np.ndarray:
    """
    @fn get_chis
    @brief Compute the RPA-corrected spin susceptibility chi_s from chi0 and the Stoner matrix.
    @param  chi0: Irreducible susceptibility [Nw, Nchi, Nchi] complex128
    @param  Smat: Spin interaction matrix [Nchi, Nchi] float64
    @return chis: Interacting spin susceptibility [Nw, Nchi, Nchi] complex128
    """
    Nchi, Nw = len(Smat), len(chi0)
    chis = np.zeros((Nw, Nchi, Nchi), dtype=np.complex128)
    flibs.get_chis.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128), # chis
        np.ctypeslib.ndpointer(dtype=np.complex128), # chi0
        np.ctypeslib.ndpointer(dtype=np.float64),    # Smat
        POINTER(c_int64), POINTER(c_int64),           # Nchi, Nw
    ]
    flibs.get_chis.restype = c_void_p
    flibs.get_chis(chis, chi0, Smat, byref(c_int64(Nchi)), byref(c_int64(Nw)))
    return chis

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
    flibs.get_scmat.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),             # Smat
        np.ctypeslib.ndpointer(dtype=np.float64),             # Cmat
        np.ctypeslib.ndpointer(dtype=np.int64),               # olist
        np.ctypeslib.ndpointer(dtype=np.int64),               # site
        POINTER(c_double), POINTER(c_double), POINTER(c_int64)  # U, J, Nchi
    ]
    flibs.get_scmat.restype = c_void_p
    flibs.get_scmat(Smat, Cmat, olist, site, byref(c_double(U)), byref(c_double(J)), byref(c_int64(Nchi)))
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
    flibs.get_scmat_orb.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64), # Smat
        np.ctypeslib.ndpointer(dtype=np.float64), # Cmat
        np.ctypeslib.ndpointer(dtype=np.float64), # Umat
        np.ctypeslib.ndpointer(dtype=np.float64), # Jmat
        np.ctypeslib.ndpointer(dtype=np.int64),   # olist
        np.ctypeslib.ndpointer(dtype=np.int64),   # site
        POINTER(c_int64), POINTER(c_int64)        # Nchi, Norb
    ]
    flibs.get_scmat_orb.restype = c_void_p
    flibs.get_scmat_orb(Smat, Cmat, olist, site, Umat, Jmat, byref(c_int64(Nchi)), byref(c_int64(Norb)))
    return Smat, Cmat

def calc_Lij(eig: np.ndarray, vk: np.ndarray, ffermi: np.ndarray, mu: float, w: float,
             idelta: float, temp: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    @fn calc_Lij
    @brief Evaluate the generalized thermoelectric transport tensors L11, L12, L22 at a given frequency.
    @param    eig: Eigenvalues [Nk, Norb] float64
    @param     vk: Group velocities [Nk, Norb, 3] complex128
    @param ffermi: Fermi-Dirac occupations [Nk, Norb] float64
    @param     mu: Chemical potential in eV
    @param      w: Frequency in eV
    @param idelta: Lorentzian broadening parameter in eV
    @param   temp: Temperature in eV
    @retval   L11: Charge conductivity tensor [3, 3] complex128
    @retval   L12: Thermoelectric tensor [3, 3] complex128
    @retval   L22: Thermal conductivity tensor [3, 3] complex128
    """
    Nk = len(eig)
    Norb = int(eig.size / Nk)
    L11 = np.zeros((3, 3), dtype=np.complex128)
    L12 = np.zeros((3, 3), dtype=np.complex128)
    L22 = np.zeros((3, 3), dtype=np.complex128)
    eps = idelta * 1e-3
    flibs.calc_lij.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128), # L11
        np.ctypeslib.ndpointer(dtype=np.complex128), # L12
        np.ctypeslib.ndpointer(dtype=np.complex128), # L22
        np.ctypeslib.ndpointer(dtype=np.complex128), # vk
        np.ctypeslib.ndpointer(dtype=np.float64),    # eig
        np.ctypeslib.ndpointer(dtype=np.float64),    # ffermi
        POINTER(c_int64), POINTER(c_int64),           # Norb, Nk
        POINTER(c_double),                           # mu
        POINTER(c_double), POINTER(c_double),         # w, idelta
        POINTER(c_double), POINTER(c_double),         # eps, temp
    ]
    flibs.calc_lij.restype = c_void_p
    flibs.calc_lij(L11, L22, L12, vk, eig, ffermi, byref(c_int64(Norb)), byref(c_int64(Nk)),
                   byref(c_double(mu)), byref(c_double(w)), byref(c_double(idelta)),
                   byref(c_double(eps)), byref(c_double(temp)))
    return L11, L12, L22

def calc_Kn(eig: np.ndarray, veloc: np.ndarray, kweight: np.ndarray, temp: float,
            mu: float, tau: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    @fn calc_Kn
    @brief Compute Boltzmann transport coefficients K0, K1, K2 for conductivity and thermopower.
    @param    eig: Eigenvalues [Nk, Norb] float64
    @param  veloc: Group velocities [Nk, Norb, 3] float64
    @param kweight: k-point weights [Nk] float64
    @param   temp: Temperature in eV
    @param     mu: Chemical potential in eV
    @param    tau: Relaxation time [Nk, Norb] float64
    @retval    K0: Charge transport tensor [3, 3] float64
    @retval    K1: Thermoelectric transport tensor [3, 3] float64
    @retval    K2: Thermal transport tensor [3, 3] float64
    """
    Nk = len(eig)
    Norb = int(eig.size / Nk)
    K0 = np.zeros((3, 3), dtype=np.float64)
    K1 = np.zeros((3, 3), dtype=np.float64)
    K2 = np.zeros((3, 3), dtype=np.float64)
    flibs.calc_kn.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64), # K0
        np.ctypeslib.ndpointer(dtype=np.float64), # K1
        np.ctypeslib.ndpointer(dtype=np.float64), # K2
        np.ctypeslib.ndpointer(dtype=np.float64), # eig
        np.ctypeslib.ndpointer(dtype=np.float64), # veloc
        np.ctypeslib.ndpointer(dtype=np.float64), # kweight
        np.ctypeslib.ndpointer(dtype=np.float64), # tau
        POINTER(c_double), POINTER(c_double),      # temp, mu
        POINTER(c_int64), POINTER(c_int64)        # Nk, Norb
    ]
    flibs.calc_kn.restype = c_void_p
    flibs.calc_kn(K0, K1, K2, eig, veloc, kweight, tau, byref(c_double(temp)),
                  byref(c_double(mu)), byref(c_int64(Nk)), byref(c_int64(Norb)))
    return K0, K1, K2

def calc_sigmahall(eig: np.ndarray, veloc: np.ndarray, imass: np.ndarray,
    kweight: np.ndarray, tau: np.ndarray, temp: float, mu: float) -> float:
    """
    @fn calc_sigmahall
    @brief Compute the Hall conductivity sigma_Hall via band integrals using velocity and inverse mass.
    @param     eig: Eigenvalues [Nk, Norb] float64
    @param   veloc: Group velocities [Nk, Norb, 3] float64
    @param   imass: Inverse effective mass tensor [Nk, Norb, 3, 3] float64
    @param kweight: k-point weights [Nk] float64
    @param     tau: Relaxation time [Nk, Norb] float64
    @param    temp: Temperature in eV
    @param      mu: Chemical potential in eV
    @return sigma_hall: Hall conductivity (scalar) float64
    """
    Nk = len(eig)
    Norb = int(eig.size / Nk)
    sigma_hall = c_double(0.0)
    flibs.calc_sigma_hall.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64), # eig
        np.ctypeslib.ndpointer(dtype=np.float64), # veloc
        np.ctypeslib.ndpointer(dtype=np.float64), # imass
        np.ctypeslib.ndpointer(dtype=np.float64), # kweight
        np.ctypeslib.ndpointer(dtype=np.float64), # tau
        POINTER(c_double), POINTER(c_double),      # temp, mu
        POINTER(c_int64), POINTER(c_int64),        # Nk, Norb
        POINTER(c_double)                         # sigma_hall
    ]
    flibs.calc_sigma_hall.restype = c_void_p
    flibs.calc_sigma_hall(eig, veloc, imass, kweight, tau, byref(c_double(temp)), byref(c_double(mu)),
                          byref(c_int64(Nk)), byref(c_int64(Norb)), byref(sigma_hall))
    return sigma_hall.value

def calc_tdf(eig: np.ndarray, veloc: np.ndarray, kweight: np.ndarray,
             tau: np.ndarray, Nw: int) -> np.ndarray:
    """
    @fn calc_tdf
    @brief Compute the transport distribution function (TDF) tensor as a function of energy.
    @param    eig: Eigenvalues [Nk, Norb] float64
    @param  veloc: Group velocities [Nk, Norb, 3] float64
    @param kweight: k-point weights [Nk] float64
    @param    tau: Relaxation time [Nk, Norb] float64
    @param     Nw: Number of energy points
    @return   tdf: Transport distribution function [Nw, 3, 3] float64
    """
    Nk = len(eig)
    Norb = int(eig.size / Nk)
    tdf = np.zeros((Nw, 3, 3), dtype=np.float64)
    flibs.calc_tdf.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64), # tdf
        np.ctypeslib.ndpointer(dtype=np.float64), # eig
        np.ctypeslib.ndpointer(dtype=np.float64), # veloc
        np.ctypeslib.ndpointer(dtype=np.float64), # kweight
        np.ctypeslib.ndpointer(dtype=np.float64), # tau
        POINTER(c_int64), POINTER(c_int64),        # Nw, Nk
        POINTER(c_int64)                           # Norb    
    ]
    flibs.calc_tdf.restype = c_void_p
    flibs.calc_tdf(tdf, eig, veloc, kweight, tau, byref(c_int64(Nw)),
                   byref(c_int64(Nk)), byref(c_int64(Norb)))
    return tdf

def get_tau(tauw: np.ndarray, eig: np.ndarray, tau_max: float, tau_mode: int, eps: float = 1.0e-4) -> np.ndarray:
    """
    @fn get_tau
    @brief Estimate k-resolved relaxation time tau(k) from frequency-dependent scattering weights.
    @param    tauw: Frequency-dependent scattering weight [Nw] float64
    @param     eig: Eigenvalues [Nk, Norb] float64
    @param tau_max: Upper bound for relaxation time values
    @param tau_mode: Mode selector for the Fortran tau integral equation
    @param     eps: Convergence tolerance
    @return    tau: k-resolved relaxation time [Nk, Norb] float64
    """
    Nk, Nw = len(eig), len(tauw)
    Norb = int(eig.size / Nk)
    tau = np.zeros((Nk, Norb), dtype=np.float64)
    flibs.get_tau.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64), # tau
        np.ctypeslib.ndpointer(dtype=np.float64), # tauw
        np.ctypeslib.ndpointer(dtype=np.float64), # eig
        POINTER(c_double), POINTER(c_double),      # tau_max, eps
        POINTER(c_int64), POINTER(c_int64),        # tau_mode, Nk
        POINTER(c_int64), POINTER(c_int64)         # Nw, Norb
    ]
    flibs.get_tau.restype = c_void_p
    flibs.get_tau(tau, tauw, eig, byref(c_double(tau_max)), byref(c_double(eps)),
                  byref(c_int64(tau_mode)), byref(c_int64(Nk)), byref(c_int64(Nw)),
                  byref(c_int64(Norb)))
    return tau

def calc_tau_epa(eig: np.ndarray, gavg: np.ndarray, wavg: np.ndarray,
                 edge: np.ndarray, step: np.ndarray, nbin: np.ndarray,
                 mu: float, temp: float) -> np.ndarray:
    """
    @fn calc_tau_epa
    @brief Compute EPA relaxation time from epa.x (job='egrid') output.
    @param   eig: Eigenvalues [Nk, Norb] float64 (eV)
    @param  gavg: EPA averaged |g|^2 [ngrid, nbin_max, nbin_max, nmodes] float64 (eV^2)
    @param  wavg: Averaged phonon frequencies per mode [nmodes] float64 (eV)
    @param  edge: Grid edges [ngrid] float64 (eV)
    @param  step: Grid steps [ngrid] float64 (eV)
    @param  nbin: Number of bins per grid [ngrid] int64
    @param    mu: Chemical potential (eV)
    @param  temp: Temperature kB*T (eV)
    @return  tau: Relaxation time [Nk, Norb] float64
    """
    Nk = len(eig)
    Norb = int(eig.size / Nk)
    ngrid = len(edge)
    nmodes = len(wavg)
    nbin_max = int(np.max(nbin))
    # Fortran expects (Norb,Nk) layout, (nmodes,nbin_max,nbin_max,ngrid) for gavg
    tau = np.zeros((Nk, Norb), dtype=np.float64)
    nbin_i64 = np.ascontiguousarray(nbin, dtype=np.int64)
    flibs.calc_tau_epa.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),  # tau
        np.ctypeslib.ndpointer(dtype=np.float64),  # gavg
        np.ctypeslib.ndpointer(dtype=np.float64),  # wavg
        np.ctypeslib.ndpointer(dtype=np.float64),  # eig
        np.ctypeslib.ndpointer(dtype=np.float64),  # edge
        np.ctypeslib.ndpointer(dtype=np.float64),  # step
        POINTER(c_double), POINTER(c_double),       # mu, temp
        POINTER(c_int64), POINTER(c_int64),         # Nk, Norb
        POINTER(c_int64),                           # nmodes
        np.ctypeslib.ndpointer(dtype=np.int64),     # nbin
        POINTER(c_int64), POINTER(c_int64)          # ngrid, nbin_max
    ]
    flibs.calc_tau_epa.restype = c_void_p
    flibs.calc_tau_epa(tau, gavg, wavg, eig, edge, step,
                       byref(c_double(mu)), byref(c_double(temp)),
                       byref(c_int64(Nk)), byref(c_int64(Norb)),
                       byref(c_int64(nmodes)), nbin_i64,
                       byref(c_int64(ngrid)), byref(c_int64(nbin_max)))
    return tau

def gen_imp_ham(rvec: np.ndarray, ham_r: np.ndarray, ham_i: np.ndarray,
                rlist: np.ndarray, imp_list: np.ndarray,eps: float = 1.0e-5) -> np.ndarray:
    """
    @fn gen_imp_ham
    @brief Construct the real-space impurity Hamiltonian from bulk hopping blocks and impurity positions.
    @param     rvec: Real-space displacement vectors [Nr, 3] float64
    @param    ham_r: Real-part hopping blocks [Nr, Norb, Norb] complex128
    @param    ham_i: Imaginary-part hopping blocks [Nr, Norb, Norb] complex128
    @param    rlist: Impurity site positions [Nsite] float64
    @param imp_list: Impurity orbital indices [Nimp] int64
    @param      eps: Small tolerance for assembly
    @return ham_imp: Impurity Hamiltonian [Norb*Nsite, Norb*Nsite] complex128
    """
    Nr, Nimp, Nsite = len(rvec), len(imp_list), len(rlist)
    Norb = int(np.sqrt(ham_r.size / Nr))
    ham_imp = np.zeros((Norb * Nsite, Norb * Nsite), dtype=np.complex128)
    flibs.gen_imp_ham.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128), # ham_imp
        np.ctypeslib.ndpointer(dtype=np.complex128), # ham_r
        np.ctypeslib.ndpointer(dtype=np.float64),    # rvec
        np.ctypeslib.ndpointer(dtype=np.complex128), # ham_i
        np.ctypeslib.ndpointer(dtype=np.int64),      # imp_list
        np.ctypeslib.ndpointer(dtype=np.float64),    # rlist
        POINTER(c_double),                           # eps
        POINTER(c_int64), POINTER(c_int64),           # Nimp, Nsite
        POINTER(c_int64), POINTER(c_int64)            # Nr, Norb
    ]
    flibs.gen_imp_ham.restype = c_void_p
    flibs.gen_imp_ham(ham_imp, ham_r, rvec, ham_i, imp_list, rlist, byref(c_double(eps)), byref(c_int64(Nimp)),
                      byref(c_int64(Nsite)), byref(c_int64(Nr)), byref(c_int64(Norb)))
    return ham_imp

def dft_imp_ham(ham_imp: np.ndarray, klist: np.ndarray, rlist: np.ndarray) -> np.ndarray:
    """
    @fn dft_imp_ham
    @brief Fourier-transform the real-space impurity Hamiltonian to k-space.
    @param  ham_imp: Real-space impurity Hamiltonian [Norb*Nsite, Norb*Nsite] complex128
    @param    klist: k-point list [Nk, 3] float64
    @param    rlist: Impurity site positions [Nsite] float64
    @return   ham_k: k-space impurity Hamiltonian [Norb*Nk, Norb*Nk] complex128
    """
    Nk, Nsite = len(klist), len(rlist)
    Norb = int(len(ham_imp) / Nsite)
    ham_k = np.zeros((Norb * Nk, Norb * Nk), dtype=np.complex128)
    flibs.get_dft_imp_ham.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128), # ham_k
        np.ctypeslib.ndpointer(dtype=np.complex128), # ham_imp
        np.ctypeslib.ndpointer(dtype=np.float64),    # klist
        np.ctypeslib.ndpointer(dtype=np.float64),    # rlist
        POINTER(c_int64), POINTER(c_int64),           # Nk, Nsite
        POINTER(c_int64)                             # Norb
    ]
    flibs.get_dft_imp_ham.restype = c_void_p
    flibs.get_dft_imp_ham(ham_k, ham_imp, klist, rlist, byref(c_int64(Nk)),
                          byref(c_int64(Nsite)), byref(c_int64(Norb)))
    return ham_k

def get_imp_spectrum(uni: np.ndarray, eigs: np.ndarray, mu: float, wlist: np.ndarray,
                     klist: np.ndarray, rlist: np.ndarray, eta: float = 1.0e-3) -> np.ndarray:
    """
    @fn get_imp_spectrum
    @brief Compute the k- and frequency-resolved impurity spectral function A(k, w).
    @param     uni: Eigenvectors of the impurity Hamiltonian [Norb*Nsite, Norb*Nsite] complex128
    @param    eigs: Eigenvalues of the impurity Hamiltonian [Norb*Nsite] float64
    @param      mu: Chemical potential in eV
    @param   wlist: Real-frequency mesh [Nw] float64
    @param   klist: k-point list [Nk, 3] float64
    @param   rlist: Impurity site positions [Nsite] float64
    @param     eta: Lorentzian broadening parameter in eV
    @return spectrum: Impurity spectral function [Nk, Nw] complex128
    """
    Nw, Nk, Nsite = len(wlist), len(klist), len(rlist)
    Norb = int(len(eigs) / Nsite)
    spectrum = np.zeros((Nk, Nw), dtype=np.complex128)
    flibs.get_spectrum_spagehtti.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128), # spectrum
        np.ctypeslib.ndpointer(dtype=np.complex128), # uni
        np.ctypeslib.ndpointer(dtype=np.float64),    # eigs
        np.ctypeslib.ndpointer(dtype=np.float64),    # klist
        np.ctypeslib.ndpointer(dtype=np.float64),    # rlist
        np.ctypeslib.ndpointer(dtype=np.float64),    # wlist
        POINTER(c_int64), POINTER(c_int64),           # Nw, Nk
        POINTER(c_int64), POINTER(c_int64),           # Nsite, Norb
        POINTER(c_double), POINTER(c_double)          # mu, eta
    ]
    flibs.get_spectrum_spagehtti.restype = c_void_p
    flibs.get_spectrum_spagehtti(spectrum, uni, eigs, klist, rlist, wlist, byref(c_int64(Nw)), byref(c_int64(Nk)),
                                 byref(c_int64(Nsite)), byref(c_int64(Norb)), byref(c_double(mu)), byref(c_double(eta)))
    return spectrum

def get_a(inp_data: np.ndarray, xlist: np.ndarray) -> np.ndarray:
    """
    @fn get_a
    @brief Compute Pade continued-fraction coefficients from input data on Matsubara grid.
    @param inp_data: Input data (e.g. Green's function) at Matsubara points [Np] complex128
    @param    xlist: Matsubara frequency points [Np] complex128
    @return       a: Pade coefficient array [Np] complex128
    """
    Np = len(inp_data)
    a = np.zeros(Np, dtype=np.complex128)
    flibs.get_a_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128), # a
        np.ctypeslib.ndpointer(dtype=np.complex128), # xlist
        np.ctypeslib.ndpointer(dtype=np.complex128), # inpdata
        POINTER(c_int64)                             # Np
    ]
    flibs.get_a_.restype = None
    flibs.get_a_(a, xlist, inp_data, byref(c_int64(Np)))
    return a

def get_QP(a: np.ndarray, xlist: np.ndarray, wlist: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    @fn get_QP
    @brief Evaluate Pade numerator Q and denominator P on the real-frequency axis from continued-fraction coefficients.
    @param      a: Pade coefficient array [Np] complex128 (from get_a)
    @param  xlist: Matsubara frequency points [Np] complex128
    @param  wlist: Real-frequency mesh [Nw] complex128
    @retval     Q: Pade numerator evaluated on real axis [Nw] complex128
    @retval     P: Pade denominator evaluated on real axis [Nw] complex128
    """
    Nw, Np = len(wlist), len(a)
    Q = np.zeros(Nw, dtype=np.complex128)
    P = np.zeros(Nw, dtype=np.complex128)
    flibs.get_qp_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128), # P
        np.ctypeslib.ndpointer(dtype=np.complex128), # Q
        np.ctypeslib.ndpointer(dtype=np.complex128), # a
        np.ctypeslib.ndpointer(dtype=np.complex128), # xlist
        np.ctypeslib.ndpointer(dtype=np.complex128), # wlist
        POINTER(c_int64), POINTER(c_int64)           # Nw, Np
    ]
    flibs.get_qp_.restype = None
    flibs.get_qp_(P, Q, a, xlist, wlist, byref(c_int64(Nw)), byref(c_int64(Np)))
    return Q, P

def pade_with_trace(A: np.ndarray, iwlist: np.ndarray, wlist: np.ndarray) -> np.ndarray:
    """
    @fn pade_with_trace
    @brief Perform Pade analytic continuation and trace over orbital indices to obtain the real-axis spectrum.
    @param      A: Matrix on Matsubara axis [Norb, Norb, Niw, Nk] complex128
    @param iwlist: Matsubara frequency points [Niw] complex128
    @param  wlist: Real-frequency mesh [Nw] complex128
    @return     B: Analytically-continued trace [Nk, Nw] complex128
    """
    Nk, Nw, Niw = len(A.T), len(wlist), len(iwlist)
    Norb = int(np.sqrt(A.size / (Nk * Niw)))
    B = np.zeros((Nk, Nw), dtype=np.complex128)
    flibs.pade_with_trace.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128), # A
        np.ctypeslib.ndpointer(dtype=np.complex128), # B
        np.ctypeslib.ndpointer(dtype=np.complex128), # iwlist
        np.ctypeslib.ndpointer(dtype=np.complex128), # wlist
        POINTER(c_int64), POINTER(c_int64),           # Nk, Niw
        POINTER(c_int64), POINTER(c_int64)            # Nw, Norb
    ]
    flibs.pade_with_trace.restype = None
    flibs.pade_with_trace(A, B, iwlist, wlist, byref(c_int64(Nk)), byref(c_int64(Niw)),
                          byref(c_int64(Nw)), byref(c_int64(Norb)))
    return B

def get_chi0(Smat: np.ndarray, Cmat: np.ndarray, Gk: np.ndarray, olist: np.ndarray,
             kmap: np.ndarray, invk: np.ndarray, temp: float, Nx: int, Ny: int, Nz: int) -> np.ndarray:
    """
    @fn get_chi0
    @brief Compute the irreducible susceptibility chi0 on the irreducible q-grid using the full BZ Green's function.
    @param   Smat: Spin interaction matrix [Nchi, Nchi] float64 (forwarded to Fortran for consistency)
    @param   Cmat: Charge interaction matrix [Nchi, Nchi] float64 (forwarded to Fortran for consistency)
    @param     Gk: Green's function on full k-grid [Norb, Norb, Nw, Nkall] complex128
    @param  olist: Orbital index list [Nchi] int64
    @param   kmap: Full-to-reduced k-point mapping [Nkall] int64
    @param   invk: Inverse k-point mapping [Nkall] int64
    @param   temp: Temperature in eV
    @param     Nx: Number of k-points along kx
    @param     Ny: Number of k-points along ky
    @param     Nz: Number of k-points along kz
    @return   chi: Irreducible susceptibility [Nchi, Nchi, Nw, Nk] complex128
    """
    Norb, Nchi = len(Gk), len(olist)
    Nkall, Nk, Nw = len(kmap), len(Gk[0, 0, 0]), len(Gk[0, 0])
    chi = np.zeros((Nchi, Nchi, Nw, Nk), dtype=np.complex128)
    flibs.get_chi0.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128), # chi
        np.ctypeslib.ndpointer(dtype=np.float64),    # Smat
        np.ctypeslib.ndpointer(dtype=np.float64),    # Cmat
        np.ctypeslib.ndpointer(dtype=np.complex128), # Gk
        np.ctypeslib.ndpointer(dtype=np.int64),      # kmap
        np.ctypeslib.ndpointer(dtype=np.int64),      # invk
        np.ctypeslib.ndpointer(dtype=np.int64),      # olist
        POINTER(c_double), POINTER(c_int64),          # temp, Nx
        POINTER(c_int64), POINTER(c_int64),           # Ny, Nz
        POINTER(c_int64), POINTER(c_int64),           # Nw, Nk
        POINTER(c_int64), POINTER(c_int64),           # Nkall, Nchi
        POINTER(c_int64)                              # Norb
    ]
    flibs.get_chi0.restype = None
    flibs.get_chi0(chi, Smat, Cmat, Gk, kmap, invk, olist, byref(c_double(temp)),
                   byref(c_int64(Nx)), byref(c_int64(Ny)), byref(c_int64(Nz)),
                   byref(c_int64(Nw)), byref(c_int64(Nk)), byref(c_int64(Nkall)),
                   byref(c_int64(Norb)), byref(c_int64(Nchi)))
    return chi

def linearized_eliashberg(chi: np.ndarray, Gk: np.ndarray, uni: np.ndarray, init_delta: np.ndarray,
                          Smat: np.ndarray, Cmat: np.ndarray, olist: np.ndarray, plist: np.ndarray,
                          kmap: np.ndarray,invk: np.ndarray, Nx: int, Ny: int, Nz: int,
                          temp: float, gap_sym: int, eps: float = 1.0e-5, itemax: int = 300,
                          arnoldi_m: int = 10) -> np.ndarray:
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
    @return    delta: Linearized gap function [Norb, Norb, Nw, Nkall] complex128
    """
    Norb, Nchi = len(Gk), len(Smat)
    Nkall, Nk, Nw = len(kmap), len(Gk[0, 0, 0]), len(Gk[0, 0])
    delta = np.zeros((Norb, Norb, Nw, Nk), dtype=np.complex128)
    flibs.lin_eliash.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128), # delta
        np.ctypeslib.ndpointer(dtype=np.complex128), # chi
        np.ctypeslib.ndpointer(dtype=np.complex128), # Gk
        np.ctypeslib.ndpointer(dtype=np.complex128), # uni
        np.ctypeslib.ndpointer(dtype=np.float64),    # init_delta
        np.ctypeslib.ndpointer(dtype=np.float64),    # Smat
        np.ctypeslib.ndpointer(dtype=np.float64),    # Cmat
        np.ctypeslib.ndpointer(dtype=np.int64),      # olist
        np.ctypeslib.ndpointer(dtype=np.float64),    # plist
        np.ctypeslib.ndpointer(dtype=np.int64),      # kmap
        np.ctypeslib.ndpointer(dtype=np.int64),      # invk
        POINTER(c_double), POINTER(c_double),         # temp, eps
        POINTER(c_int64), POINTER(c_int64),           # Nkall, Nk
        POINTER(c_int64), POINTER(c_int64),           # Nw, Nchi
        POINTER(c_int64), POINTER(c_int64),           # Norb, Nx
        POINTER(c_int64), POINTER(c_int64),           # Ny, Nz
        POINTER(c_int64), POINTER(c_int64),           # itemax, gapsym
        POINTER(c_int64)                              # arnoldi_m
    ]
    flibs.lin_eliash.restype = None
    flibs.lin_eliash(delta, chi, Gk, uni, init_delta, Smat, Cmat, olist, plist, kmap, invk,
                     byref(c_double(temp)), byref(c_double(eps)), byref(c_int64(Nkall)),
                     byref(c_int64(Nk)), byref(c_int64(Nw)), byref(c_int64(Nchi)),
                     byref(c_int64(Norb)), byref(c_int64(Nx)), byref(c_int64(Ny)),
                     byref(c_int64(Nz)), byref(c_int64(itemax)), byref(c_int64(gap_sym)),
                     byref(c_int64(arnoldi_m)))
    return delta

def linearized_eliashberg_soc(chi: np.ndarray, Gk: np.ndarray, uni: np.ndarray, init_delta: np.ndarray,
                              Vmat: np.ndarray, sgnsig: np.ndarray, sgnsig2: np.ndarray, plist: np.ndarray,
                              slist: np.ndarray, olist: np.ndarray, kmap: np.ndarray, invk: np.ndarray,
                              invs: np.ndarray, invschi: np.ndarray, Nx: int, Ny: int, Nz: int, temp: float,
                              gap_sym: int, eps: float = 1.0e-4, itemax: int = 300,
                              arnoldi_m: int = 10) -> np.ndarray:
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
    """
    Norb, Nchi = len(slist), len(Vmat)
    Nkall, Nk, Nw = len(kmap), len(Gk[0, 0, 0]), len(Gk[0, 0])
    delta = np.zeros((Norb, Norb, Nw, Nkall), dtype=np.complex128)
    flibs.lin_eliash_soc.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128), # delta
        np.ctypeslib.ndpointer(dtype=np.complex128), # chi
        np.ctypeslib.ndpointer(dtype=np.complex128), # Gk
        np.ctypeslib.ndpointer(dtype=np.complex128), # uni
        np.ctypeslib.ndpointer(dtype=np.float64),    # init_delta
        np.ctypeslib.ndpointer(dtype=np.float64),    # Vmat
        np.ctypeslib.ndpointer(dtype=np.float64),    # sgnsig
        np.ctypeslib.ndpointer(dtype=np.float64),    # sgnsig2
        np.ctypeslib.ndpointer(dtype=np.float64),    # plist
        np.ctypeslib.ndpointer(dtype=np.int64),      # olist
        np.ctypeslib.ndpointer(dtype=np.int64),      # slist
        np.ctypeslib.ndpointer(dtype=np.int64),      # kmap
        np.ctypeslib.ndpointer(dtype=np.int64),      # invk
        np.ctypeslib.ndpointer(dtype=np.int64),      # invs
        np.ctypeslib.ndpointer(dtype=np.int64),      # invschi
        POINTER(c_double), POINTER(c_double),         # temp, eps
        POINTER(c_int64), POINTER(c_int64),           # Nkall, Nk
        POINTER(c_int64), POINTER(c_int64),           # Nw, Nchi
        POINTER(c_int64), POINTER(c_int64),           # Norb, Nx
        POINTER(c_int64), POINTER(c_int64),           # Ny, Nz
        POINTER(c_int64), POINTER(c_int64),           # itemax, gapsym
        POINTER(c_int64)                              # arnoldi_m
    ]
    flibs.lin_eliash_soc.restype = None
    flibs.lin_eliash_soc(delta, chi, Gk, uni, init_delta, Vmat, sgnsig, sgnsig2, plist, olist,
                         slist, kmap, invk, invs, invschi, byref(c_double(temp)),
                         byref(c_double(eps)), byref(c_int64(Nkall)), byref(c_int64(Nk)),
                         byref(c_int64(Nw)), byref(c_int64(Nchi)), byref(c_int64(Norb)),
                         byref(c_int64(Nx)), byref(c_int64(Ny)), byref(c_int64(Nz)),
                         byref(c_int64(itemax)), byref(c_int64(gap_sym)),
                         byref(c_int64(arnoldi_m)))
    return delta

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
    flibs.get_vmat_soc.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64), # Vmat
        np.ctypeslib.ndpointer(dtype=np.int64),   # olist
        np.ctypeslib.ndpointer(dtype=np.int64),   # slist
        np.ctypeslib.ndpointer(dtype=np.int64),   # site
        np.ctypeslib.ndpointer(dtype=np.int64),   # invs
        POINTER(c_double), POINTER(c_double),      # U, J
        POINTER(c_int64), POINTER(c_int64)        # Nchi, Norb
    ]
    flibs.get_vmat_soc.restype = None
    flibs.get_vmat_soc(Vmat, olist, slist, site, invs, byref(c_double(U)),
                       byref(c_double(J)), byref(c_int64(Nchi)), byref(c_int64(Norb)))
    return Vmat

def gen_Vmatrix_orb(olist: np.ndarray, slist: np.ndarray, site: np.ndarray, invs: np.ndarray,
                    Umat: np.ndarray, Jmat: np.ndarray,) -> np.ndarray:
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
    flibs.get_vmat_soc_orb.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64), # Vmat
        np.ctypeslib.ndpointer(dtype=np.int64),   # olist
        np.ctypeslib.ndpointer(dtype=np.int64),   # slist
        np.ctypeslib.ndpointer(dtype=np.int64),   # site
        np.ctypeslib.ndpointer(dtype=np.int64),   # invs
        np.ctypeslib.ndpointer(dtype=np.float64), # Umat
        np.ctypeslib.ndpointer(dtype=np.float64), # Jmat
        POINTER(c_int64), POINTER(c_int64)        # Nchi, Norb
    ]
    flibs.get_vmat_soc_orb.restype = None
    flibs.get_vmat_soc_orb(Vmat, olist, slist, site, invs, Umat, Jmat,
                           byref(c_int64(Nchi)), byref(c_int64(Norb)))
    return Vmat

def get_chi0_soc(Vmat: np.ndarray, Gk: np.ndarray, olist: np.ndarray, slist: np.ndarray,
                 kmap: np.ndarray, invk: np.ndarray, invs: np.ndarray, temp: float,
                 Nx: int, Ny: int, Nz: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    @fn get_chi0_soc
    @brief Compute the SOC-aware irreducible susceptibility chi0 and auxiliary spin-sign matrices.
    @param   Vmat: SOC interaction matrix [Nchi, Nchi] float64
    @param     Gk: Green's function on full k-grid [Norb, Norb, Nw, Nkall] complex128
    @param  olist: Orbital index list [Nchi] int64
    @param  slist: Spin-orbit label for each orbital [Norb] int64
    @param   kmap: Full-to-reduced k-point mapping [Nkall] int64
    @param   invk: Inverse k-point mapping [Nkall] int64
    @param   invs: Inverse spin index mapping [Norb] int64
    @param   temp: Temperature in eV
    @param     Nx: Number of k-points along kx
    @param     Ny: Number of k-points along ky
    @param     Nz: Number of k-points along kz
    @retval    chi: Irreducible susceptibility [Nchi, Nchi, Nw, Nk] complex128
    @retval sgnsig: Sign matrix [Norb, Norb] float64
    @retval sgnsig2: Sign matrix for chi orbital mapping [Nchi, Nchi] float64
    @retval invschi: Inverse chi orbital mapping [Nchi] int64
    """
    Norb, Nchi = len(slist), len(olist)
    Nkall, Nk, Nw = len(kmap), len(Gk[0, 0, 0]), len(Gk[0, 0])
    chi = np.zeros((Nchi, Nchi, Nw, Nk), dtype=np.complex128)
    sgnsig = np.zeros((Norb, Norb), dtype=np.float64)
    sgnsig2 = np.zeros((Nchi, Nchi), dtype=np.float64)
    invschi = np.zeros(Nchi, dtype=np.int64)
    flibs.get_chi0_soc.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128), # chi
        np.ctypeslib.ndpointer(dtype=np.float64),    # sgnsig
        np.ctypeslib.ndpointer(dtype=np.float64),    # sgnsig2
        np.ctypeslib.ndpointer(dtype=np.int64),      # invschi
        np.ctypeslib.ndpointer(dtype=np.float64),    # Vmat
        np.ctypeslib.ndpointer(dtype=np.complex128), # Gk
        np.ctypeslib.ndpointer(dtype=np.int64),      # kmape
        np.ctypeslib.ndpointer(dtype=np.int64),      # invk
        np.ctypeslib.ndpointer(dtype=np.int64),      # invs
        np.ctypeslib.ndpointer(dtype=np.int64),      # olist
        np.ctypeslib.ndpointer(dtype=np.int64),      # slist
        POINTER(c_double), POINTER(c_int64),         # temp, Nx
        POINTER(c_int64), POINTER(c_int64),           # Ny, Nz
        POINTER(c_int64), POINTER(c_int64),           # Nw, Nk
        POINTER(c_int64), POINTER(c_int64),           # Nkall, Nchi
        POINTER(c_int64)                              # Norb
    ]
    flibs.get_chi0_soc.restype = None
    flibs.get_chi0_soc(chi, sgnsig, sgnsig2, invschi, Vmat, Gk, kmap, invk, invs, olist, slist,
                       byref(c_double(temp)), byref(c_int64(Nx)), byref(c_int64(Ny)),
                       byref(c_int64(Nz)), byref(c_int64(Nw)), byref(c_int64(Nk)),
                       byref(c_int64(Nkall)), byref(c_int64(Nchi)), byref(c_int64(Norb)))
    return chi, sgnsig, sgnsig2, invschi

def get_chis_chic_soc(chi: np.ndarray, Vmat: np.ndarray, olist: np.ndarray, slist: np.ndarray,
                      invs: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    @fn get_chis_chic_soc
    @brief Decompose the SOC susceptibility into charge (chic), longitudinal spin (chiszz), and transverse (chispm) channels.
    @param    chi: Full susceptibility [Nchi, Nchi, Nk, Nw] complex128
    @param   Vmat: SOC interaction matrix [Nchi, Nchi] float64
    @param  olist: Orbital index list [Nchi] int64
    @param  slist: Spin-orbit label for each orbital [Norb] int64
    @param   invs: Inverse spin index mapping [Norb] int64
    @retval  chic: Charge susceptibility [Nchi/4, Nchi/4, Nk] complex128
    @retval chiszz: Longitudinal spin susceptibility [Nchi/4, Nchi/4, Nk] complex128
    @retval chispm: Transverse spin susceptibility [Nchi/4, Nchi/4, Nk] complex128
    """
    def get_orb_list(olist, slist, invs):
        orb_list = np.zeros(Nchi, dtype=np.int64)
        orbs = np.where(slist == 1)
        o1, o2 = np.meshgrid(orbs, orbs)
        orb_chi = np.array([o1.flatten(), o2.flatten()]).T
        for i, ol in enumerate(olist):
            o1_val = (ol[0] if slist[ol[0] - 1] == 1 else invs[ol[0] - 1]) - 1
            o2_val = (ol[1] if slist[ol[1] - 1] == 1 else invs[ol[1] - 1]) - 1
            for j, ob2 in enumerate(orb_chi):
                if o1_val == ob2[0] and o2_val == ob2[1]:
                    orb_list[i] = j + 1
                    break
        return orb_list
    Norb, Nchi = len(slist), len(olist)
    Nk, Nw = len(chi[0, 0, 0]), len(chi[0, 0])
    orb_list = get_orb_list(olist, slist, invs)
    chic = np.zeros((int(Nchi / 4), int(Nchi / 4), Nk), dtype=np.complex128)
    chiszz = np.zeros((int(Nchi / 4), int(Nchi / 4), Nk), dtype=np.complex128)
    chispm = np.zeros((int(Nchi / 4), int(Nchi / 4), Nk), dtype=np.complex128)
    flibs.get_chis_chic_soc.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128), # chic
        np.ctypeslib.ndpointer(dtype=np.complex128), # chiszz
        np.ctypeslib.ndpointer(dtype=np.complex128), # chispm
        np.ctypeslib.ndpointer(dtype=np.complex128), # chi
        np.ctypeslib.ndpointer(dtype=np.float64),    # Vmat
        np.ctypeslib.ndpointer(dtype=np.int64),      # orb_list
        np.ctypeslib.ndpointer(dtype=np.int64),      # olist
        np.ctypeslib.ndpointer(dtype=np.int64),      # slist
        np.ctypeslib.ndpointer(dtype=np.int64),      # invs
        POINTER(c_int64), POINTER(c_int64),           # Nk, Nw
        POINTER(c_int64), POINTER(c_int64)            # Nchi, Norb
    ]
    flibs.get_chis_chic_soc.restype = None
    flibs.get_chis_chic_soc(chic, chiszz, chispm, chi, Vmat, orb_list, olist, slist, invs, byref(c_int64(Nk)),
                            byref(c_int64(Nw)), byref(c_int64(Nchi)), byref(c_int64(Norb)))
    return chic, chiszz, chispm

def get_chis_chic(chi: np.ndarray, Smat: np.ndarray, Cmat: np.ndarray,) -> tuple[np.ndarray, np.ndarray]:
    """
    @fn get_chis_chic
    @brief Compute RPA-corrected spin susceptibility chis and charge susceptibility chic from chi0.
    @param   chi: Irreducible susceptibility [Nchi, Nchi, Nk, Nw] complex128
    @param  Smat: Spin interaction matrix [Nchi, Nchi] float64
    @param  Cmat: Charge interaction matrix [Nchi, Nchi] float64
    @retval chis: Interacting spin susceptibility [Nchi, Nchi, Nk] complex128
    @retval chic: Interacting charge susceptibility [Nchi, Nchi, Nk] complex128
    """
    Nk, Nw, Nchi = len(chi[0, 0, 0]), len(chi[0, 0]), len(Smat)
    chis = np.zeros((Nchi, Nchi, Nk), dtype=np.complex128)
    chic = np.zeros((Nchi, Nchi, Nk), dtype=np.complex128)
    flibs.get_chis_chic.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),        # chis
        np.ctypeslib.ndpointer(dtype=np.complex128),        # chic
        np.ctypeslib.ndpointer(dtype=np.complex128),        # chi
        np.ctypeslib.ndpointer(dtype=np.float64),           # Smat
        np.ctypeslib.ndpointer(dtype=np.float64),           # Cmat
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64) # Nk, Nw, Nchi
    ]
    flibs.get_chis_chic.restype = None
    flibs.get_chis_chic(chis, chic, chi, Smat, Cmat, byref(c_int64(Nk)),
                        byref(c_int64(Nw)), byref(c_int64(Nchi)))
    return chis, chic

def conv_delta_orb_to_band(delta: np.ndarray, uni: np.ndarray, invk: np.ndarray, plist: np.ndarray,gap_sym) -> np.ndarray:
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
    flibs.conv_delta_orb_to_band.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),         # deltab
        np.ctypeslib.ndpointer(dtype=np.complex128),         # delta
        np.ctypeslib.ndpointer(dtype=np.complex128),         # uni
        np.ctypeslib.ndpointer(dtype=np.float64),            # plist
        np.ctypeslib.ndpointer(dtype=np.int64),              # invk
        POINTER(c_int64), POINTER(c_int64),                  # Norb, Nkall
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64) # Nk, Nw, gap_sym
    ]
    flibs.conv_delta_orb_to_band.restype = None
    flibs.conv_delta_orb_to_band(deltab, delta, uni, plist, invk, byref(c_int64(Norb)),
                                 byref(c_int64(Nkall)), byref(c_int64(Nk)), byref(c_int64(Nw)), byref(c_int64(gap_sym)))
    return deltab

def conv_delta_orb_to_band_soc(delta: np.ndarray, uni: np.ndarray, invk: np.ndarray,
                               invs: np.ndarray,slist: np.ndarray) -> np.ndarray:
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
    flibs.conv_delta_orb_to_band_soc.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128), # deltab
        np.ctypeslib.ndpointer(dtype=np.complex128), # delta
        np.ctypeslib.ndpointer(dtype=np.complex128), # uni
        np.ctypeslib.ndpointer(dtype=np.int64),      # invk
        np.ctypeslib.ndpointer(dtype=np.int64),      # invs
        np.ctypeslib.ndpointer(dtype=np.int64),      # slist
        POINTER(c_int64), POINTER(c_int64),           # Norb, Nkall
        POINTER(c_int64), POINTER(c_int64)            # Nk, Nw
    ]
    flibs.conv_delta_orb_to_band_soc.restype = None
    flibs.conv_delta_orb_to_band_soc(deltab, delta, uni, invk, invs, slist, byref(c_int64(Norb)),
                                     byref(c_int64(Nkall)), byref(c_int64(Nk)), byref(c_int64(Nw)))
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
    flibs.mkfk_trs_nsoc_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),         # Fk
        np.ctypeslib.ndpointer(dtype=np.complex128),         # Gk
        np.ctypeslib.ndpointer(dtype=np.complex128),         # delta
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64) # Nk, Nw, Norb
    ]
    flibs.mkfk_trs_nsoc_.restype = None
    flibs.mkfk_trs_nsoc_(Fk, Gk, delta, byref(c_int64(Nk)), byref(c_int64(Nw)), byref(c_int64(Norb)))

    return Fk

def gen_Fk_soc(Gk: np.ndarray, delta: np.ndarray, invk: np.ndarray, invs: np.ndarray,
               slist: np.ndarray,gap_sym: int) -> np.ndarray:
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
    sgnsig=np.array([slist]).T.dot(np.array([slist])).astype(np.float64)
    flibs.mkfk_trs_soc_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128), # Fk
        np.ctypeslib.ndpointer(dtype=np.complex128), # Gk
        np.ctypeslib.ndpointer(dtype=np.complex128), # delta
        np.ctypeslib.ndpointer(dtype=np.float64),    # sgnsig
        np.ctypeslib.ndpointer(dtype=np.int64),      # slist
        np.ctypeslib.ndpointer(dtype=np.int64),      # invk
        np.ctypeslib.ndpointer(dtype=np.int64),      # invs
        POINTER(c_int64), POINTER(c_int64),          # Nkall, Nk
        POINTER(c_int64), POINTER(c_int64),          # Nw, Norb
        POINTER(c_int64)                             # gap_sym
    ]
    flibs.mkfk_trs_soc_.restype = None
    flibs.mkfk_trs_soc_(Fk, Gk, delta, sgnsig, slist, invk, invs, byref(c_int64(Nkall)),
                        byref(c_int64(Nk)), byref(c_int64(Nw)), byref(c_int64(Norb)), byref(c_int64(gap_sym)))
    return Fk

def remap_gap(delta0,plist,invk,gap_sym):
    """
    @fn remap_gap
    @brief Remap the anomalous self-energy from the reduced k-grid to the full k-grid applying gap symmetry.
    @param  delta0: Anomalous self-energy on reduced k-grid [Norb, Norb, Nw, Nk] complex128
    @param   plist: Parity list for gap symmetry [Norb] float64
    @param    invk: Inverse k-point mapping [Nkall] int64
    @param gap_sym: Gap symmetry index
    @return  delta: Anomalous self-energy on full k-grid [Norb, Norb, Nw, Nkall] complex128
    """
    Nkall,Nk,Norb=len(invk),len(delta0.T),len(plist)
    Nw=int(delta0.size/(Nk*Norb*Norb))
    delta=np.zeros((Norb, Norb, Nw, Nkall), dtype=np.complex128)
    flibs.remap_delta_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128), # delta
        np.ctypeslib.ndpointer(dtype=np.complex128), # delta0
        np.ctypeslib.ndpointer(dtype=np.float64),    # plist
        np.ctypeslib.ndpointer(dtype=np.int64),      # invk
        POINTER(c_int64), POINTER(c_int64),          # Nkall, Nk
        POINTER(c_int64), POINTER(c_int64),          # Nw, Norb
        POINTER(c_int64)                             # gap_sym
        ]
    flibs.remap_delta_.restype = None
    flibs.remap_delta_(delta, delta0, plist, invk, byref(c_int64(Nkall)), byref(c_int64(Nk)),
                       byref(c_int64(Nw)), byref(c_int64(Norb)), byref(c_int64(gap_sym)))
    return delta

def get_eig_or_tr_chi(chi: np.ndarray, invk: np.ndarray, sw_eig:bool) -> np.ndarray:
    """
    @fn get_eig_or_tr_chi
    @brief Compute either the leading eigenvalue or the trace of the susceptibility on the full k-grid.
    @param    chi: Susceptibility [Nchi, Nchi, Nk] complex128
    @param   invk: Inverse k-point mapping [Nkall] int64
    @param sw_eig: If True, return leading eigenvalue; if False, return trace
    @return  chiq: Eigenvalue or trace on the full k-grid [Nkall] complex128
    """
    Nkall,Nk=len(invk),len(chi.T)
    Nchi=int(np.sqrt(chi.size/Nk))
    chiq=np.zeros(Nkall,dtype=np.complex128)
    flibs.get_eig_or_tr_chi.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128), # chiq
        np.ctypeslib.ndpointer(dtype=np.complex128), # chi
        np.ctypeslib.ndpointer(dtype=np.int64),      # invk
        POINTER(c_int64), POINTER(c_int64),          # Nkall, Nk
        POINTER(c_int64), POINTER(c_bool)            # Nchi, sw_eig
    ]
    flibs.get_eig_or_tr_chi.restype = None
    flibs.get_eig_or_tr_chi(chiq,chi, invk, byref(c_int64(Nkall)), byref(c_int64(Nk)),
                             byref(c_int64(Nchi)), byref(c_bool(sw_eig)))
    return chiq

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
        if Nz % 2 == 0:  # all even Nk=Nkall/2+4 else Nz even Nk=Nkall+2
            Nk = int(Nkall / 2 + 4) if ((Nx % 2) == 0 and (Ny % 2) == 0) else (
                int(Nkall / 2 + 2) if ((Nx % 2) == 0 or (Ny % 2) == 0) else int(Nkall / 2 + 1)
            )
        else:  # NxNy plane even Nk=Nkall/2+2 else Nk=Nkall/2+1
            Nk = int(Nkall / 2 + 2) if ((Nx % 2) == 0 and (Ny % 2) == 0) else int(Nkall / 2 + 1)
    else:  # all odd Nk=(Nkall+1)/2
        Nk = int((Nkall + 1) / 2)
    klist = np.zeros((Nk, 3), dtype=np.float64)
    kmap = np.zeros((Nkall, 3), dtype=np.int64)
    invk_ft_list = np.zeros((Nkall, 3), dtype=np.int64)
    flibs.generate_irr_kpoint_inv.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64), # klist
        np.ctypeslib.ndpointer(dtype=np.int64),   # kmap
        np.ctypeslib.ndpointer(dtype=np.int64),   # invk_ft_list
        POINTER(c_int64), POINTER(c_int64),        # Nk, Nx
        POINTER(c_int64), POINTER(c_int64)         # Ny, Nz
    ]
    flibs.generate_irr_kpoint_inv.restype = None
    flibs.generate_irr_kpoint_inv(klist, kmap, invk_ft_list, byref(c_int64(Nk)),
                                  byref(c_int64(Nx)), byref(c_int64(Ny)), byref(c_int64(Nz)))
    return klist, kmap, invk_ft_list

def gen_kpoint_weight(invk,Nk):
    """
    @fn gen_kpoint_weight
    @brief Compute k-point weights for the irreducible BZ by counting degeneracy from the full k-grid mapping.
    @param  invk: Inverse k-point mapping [Nkall] int64
    @param    Nk: Number of irreducible k-points
    @return weight: k-point weights [Nk] float64
    """
    Nkall=len(invk)
    weight=np.zeros(Nk,dtype=np.float64)
    flibs.get_kweight.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64), # invk
        np.ctypeslib.ndpointer(dtype=np.int64),   # weight
        POINTER(c_int64), POINTER(c_int64)        # Nk, Nkall
    ]
    flibs.get_kweight.restype = None
    flibs.get_kweight(weight,invk,byref(c_int64(Nk)),byref(c_int64(Nkall)))
    return weight

def get_plist(rvec,ham_r):
    """
    @fn get_plist
    @brief Determine the parity eigenvalues of each orbital from the Hamiltonian's inversion symmetry properties.
    @param  rvec: Real-space displacement vectors [Nr, 3] float64
    @param ham_r: Real-space Hamiltonian blocks [Nr, Norb, Norb] complex128
    @return plist: Parity sign (+1 or -1) for each orbital [Norb] float64
    """
    Nr=len(rvec)
    Norb=int(np.sqrt(ham_r.size/Nr))
    Pmn=np.zeros((Norb,Norb),dtype=np.float64)
    flibs.get_parity_prop.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),    #Pmn
        np.ctypeslib.ndpointer(dtype=np.float64),    #rvec
        np.ctypeslib.ndpointer(dtype=np.complex128), #ham_r
        POINTER(c_int64), POINTER(c_int64)           #Norb,Nr
    ]
    flibs.get_parity_prop.restype = None
    flibs.get_parity_prop(Pmn,rvec,ham_r,byref(c_int64(Norb)),byref(c_int64(Nr)))
    return np.sign(Pmn[0,:])

def solve_cpa(hamk: np.ndarray, VA: np.ndarray, VB: np.ndarray,
              x: float, zlist: np.ndarray, pp: float = 0.5,
              maxiter: int = 500, tol: float = 1.0e-10,
              sigma_init: np.ndarray | None = None) -> np.ndarray:
    """
    @fn solve_cpa
    @brief Run CPA self-consistent loop over an array of complex frequencies.
    @param   hamk: k-space Hamiltonian [Nk, Norb, Norb] complex128
    @param     VA: Onsite potential of species A [Norb, Norb] complex128
    @param     VB: Onsite potential of species B [Norb, Norb] complex128
    @param      x: Concentration of species A (0 < x < 1)
    @param  zlist: Complex frequency array [Nw] complex128 (松原 iω_n or 実軸 ω+iδ)
    @param     pp: Linear mixing parameter (default 0.5)
    @param maxiter: Maximum CPA iterations per frequency (default 500)
    @param    tol: Convergence tolerance (default 1e-10)
    @param sigma_init: Initial guess [Nw, Norb, Norb] complex128 (default: VCA = x*VA + (1-x)*VB)
    @return sigma_cpa: Converged CPA self-energy [Nw, Norb, Norb] complex128
    """
    Nk = len(hamk)
    Norb = hamk.shape[1]
    Nw = len(zlist)
    VA_c = np.ascontiguousarray(VA, dtype=np.complex128)
    VB_c = np.ascontiguousarray(VB, dtype=np.complex128)
    zlist_c = np.ascontiguousarray(zlist, dtype=np.complex128)

    # initial guess: VCA
    if sigma_init is not None:
        sigma_cpa = np.ascontiguousarray(sigma_init, dtype=np.complex128)
    else:
        sigma_cpa = np.zeros((Nw, Norb, Norb), dtype=np.complex128)
        vca = x * VA_c + (1.0 - x) * VB_c
        for iw in range(Nw):
            sigma_cpa[iw, :, :] = vca

    flibs.solve_cpa_array.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),  # sigma_cpa_w
        np.ctypeslib.ndpointer(dtype=np.complex128),  # hamk
        np.ctypeslib.ndpointer(dtype=np.complex128),  # VA
        np.ctypeslib.ndpointer(dtype=np.complex128),  # VB
        POINTER(c_double),                             # x
        np.ctypeslib.ndpointer(dtype=np.complex128),  # zlist
        POINTER(c_double),                             # pp
        POINTER(c_int64), POINTER(c_int64),            # Nk, Norb
        POINTER(c_int64), POINTER(c_int64),            # Nw, maxiter
        POINTER(c_double)                              # tol
    ]
    flibs.solve_cpa_array.restype = None
    flibs.solve_cpa_array(sigma_cpa, hamk, VA_c, VB_c,
                          byref(c_double(x)), zlist_c,
                          byref(c_double(pp)),
                          byref(c_int64(Nk)), byref(c_int64(Norb)),
                          byref(c_int64(Nw)), byref(c_int64(maxiter)),
                          byref(c_double(tol)))
    return sigma_cpa
