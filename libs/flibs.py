from ctypes import *
import numpy as np
#import fortran library
flibs=np.ctypeslib.load_library("libs/libfmod.so",".")
#interface for fmod subroutines

def omp_params() ->tuple[int,bool]:
    """Query the underlying library for OpenMP configuration.
    Returns a tuple (num_threads, valid_flag).  `num_threads` is the
    maximum number of threads the Fortran code will use; `valid_flag`
    indicates whether OpenMP was initialized correctly.
    """
    omp_num=c_int64()
    omp_check=c_bool()
    flibs.openmp_params.argtypes=[POINTER(c_int64),POINTER(c_bool)]
    flibs.openmp_params.restype=c_void_p
    flibs.openmp_params(byref(omp_num),byref(omp_check))
    return omp_num.value,omp_check.value

def gen_ham(klist: np.ndarray, ham_r: np.ndarray, rvec: np.ndarray,
            Ovl_r: np.ndarray | None = None) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    """Compute k‑space Hamiltonian from real‑space hopping blocks.
    Parameters
    ----------
    klist : (Nk,3) float64 list of reciprocal vectors.
    ham_r : (Nr, Norb, Norb) complex128 real‑space Hamiltonian blocks.
    rvec : (Nr,3) float64 corresponding real‑space vectors.
    Ovl_r : same as ``ham_r``, optional overlap blocks; if provided the function returns a tuple ``(hamk,Ovlk)``.

    Returns
    -------
    hamk : (Nk, Norb, Norb) complex128 Hamiltonian in k space.
    Ovlk : (Nk, Norb, Norb) complex128, optional Overlap matrix, only returned when ``Ovl_r`` is given.
    """
    Nk, Nr = len(klist), len(rvec)
    assert ham_r.ndim == 3, "ham_r must be a 3‑D array"
    Norb = int(np.sqrt(ham_r.size / Nr))
    hamk = np.zeros((Nk, Norb, Norb), dtype=np.complex128)
    flibs.gen_ham.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.complex128),
        np.ctypeslib.ndpointer(dtype=np.float64),
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64)
    ]
    flibs.gen_ham.restype = c_void_p
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
    """Diagonalize Hamiltonian (and overlap if provided).
    Parameters
    ----------
    hamk : (Nk, Norb, Norb) complex128 k‑space Hamiltonian.
    Ovlk : same shape as ``hamk`` or ``None`` overlap matrices; if supplied the generalized eigenvalue problem is solved.
    sw : bool if ``True`` return both eigenvalues and eigenvectors, otherwise just eigenvalues.

    Returns
    -------
    eig : (Nk, Norb) float64 eigenvalues for each k point.
    uni : (Nk, Norb, Norb) complex128 eigenvectors (only returned when ``sw`` is ``True``).
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
        flibs.get_eig.restype = c_void_p
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
    """Compute Fermi–Dirac occupations from eigenvalues.
    Parameters
    ----------
    eig : (Nk, Norb) float64 eigenvalues for each k point.
    mu : float chemical potential.
    temp : float temperature (same units as eigenvalues).

    Returns
    -------
    ffermi : (Nk, Norb) float64
        occupation numbers, one per state.
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
    """Compute band velocity matrices (unrotated) from real‑space data.
    Parameters
    ----------
    klist : (Nk,3) float64 list of k points.
    ham_r : (Nr, Norb, Norb) complex128 real‑space Hamiltonian blocks.
    rvec : (Nr,3) float64 real‑space displacement vectors.
    Returns
    -------
    vk : (Nk, Norb, Norb, 3) complex128 velocity matrix elements for each k point and Cartesian direction.
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
    """Transform orbital basis velocity data to band basis using rotation matrices.
    Parameters
    ----------
    vk0 : (Nk, Norb, Norb, 3) complex128 raw velocity matrices.
    mrot : (Nk, Norb, Norb) float64 rotation matrices for each k point.
    uni : (Nk, Norb, Norb) complex128 eigenvector matrices.

    Returns
    -------
    vk : (Nk, Norb, 3) float64 velocities along each Cartesian direction.
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
    """Compute velocity matrix elements with band indices.
    Produces a four‑dimensional array containing velocity
    elements between bands n and m for the three Cartesian
    directions.
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
    """Convenience wrapper returning band‑indexed velocities.
    This composes :func:`get_vlm0` and :func:`get_vnm` for a single call.
    """
    vk0 = get_vlm0(klist, ham_r, rvec)
    vk = get_vnm(vk0, mrot, uni)
    return vk

def get_veloc(klist: np.ndarray, ham_r: np.ndarray, rvec: np.ndarray,
              mrot: np.ndarray, uni: np.ndarray) -> np.ndarray:
    """Compute velocity expectation values in rotated band basis.
    Shortcut for :func:`get_vlm0` followed by :func:`get_vk`.
    """
    vk0 = get_vlm0(klist, ham_r, rvec)
    vk = get_vk(vk0, mrot, uni)
    return vk

def get_imass0(klist: np.ndarray, ham_r: np.ndarray, rvec: np.ndarray) -> np.ndarray:
    """Calculate the unrotated inverse of mass tensor at k points.
    Parameters
    ----------
    klist : (Nk,3) float64
    ham_r : (Nr, Norb, Norb) complex128
    rvec : (Nr,3) float64

    Returns
    -------
    imass0 : (Nk, Norb, Norb, 3, 3) complex128 Inverse of mass tensor before band rotation.
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
    """Rotate the inverse of mass tensor into band basis.
    Parameters are the output of :func:`get_imass0`, the rotation matrix ``mrot`` and eigenvectors ``uni``.
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
    """Return either the imaginary mass tensor or its inverse (mass).
    Parameters
    ----------
    sw_imass : bool if ``True`` return the rotated imaginary mass; otherwise return its matrix inverse using ``scipy.linalg.inv``.
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
    """Construct non‑interacting Green's function array G(k,iω)
    This routine wraps the Fortran subroutine ``gen_green0_`` which
    builds the frequency‑dependent Green's function from eigenvalues
    and eigenvectors of the Hamiltonian.  The output shape is
    ``(Norb, Norb, Nw, Nk)`` with complex entries.
    Parameters
    ----------
    eig : (Nk, Norb) float64 eigenvalues at each k‑point
    uni : (Nk, Norb, Norb) complex128 corresponding eigenvectors
    mu : float chemical potential
    temp : float temperature (same units as eigenvalues)
    Nw : int number of Matsubara frequencies (imaginary axis points)
    Returns
    -------
    Gk : ndarray[complex128] non‑interacting Green's function tensor
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
    """Generate interacting Green's function including self‑energy.
    Parameters
    ----------
    selfen : (Norb,Norb,Nw,Nk) complex128 self‑energy tensor in band basis
    hamk : (Nk, Norb, Norb) complex128 k‑space Hamiltonian
    mu : float chemical potential
    temp : float temperature
    Returns
    -------
    Gk : ndarray[complex128] full (interacting) Green's function array of shape ``(Norb,Norb,Nw,Nk)``
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
    """Build Green's function using eigvalues.
    This variant is similar to :func:`gen_green` but accepts the
    eigenvalues ``eig`` and eigenvectors ``uni`` separately.  It can
    be convenient when those quantities are already available and
    avoids recalculating them inside the Fortran routine.
    Parameters and return value are identical to :func:`gen_green`.
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
    """Compute imaginary part of trace of Green's function (corresponding to pDos).
    Parameters
    ----------
    eig : (Nk, Norb) float64 band energies
    mu : float chemical potential
    wlist : (Nw,) float64 frequency mesh on real axis
    delta : float small broadening parameter

    Returns
    -------
    trG : (Nk, Nw) complex128 imaginary part of the k‑resolved trace of the Green's function,
    returned with a negative sign to match physics conventions.
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
    """Partial density of states generated from Green's function.
    The Fortran routine ``gen_dos`` evaluates the k‑dependent density
    of states using eigenvectors ``uni`` and eigenvalues ``eig``.  The
    returned array has shape ``(Norb, Nw)`` and is complex; the
    physical DOS is the negative imaginary part.
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
    """Obtain the irreducible susceptibility chi0 using convolution over the full BZ.
    This wrapper calls the Fortran subroutine ``get_chi0_conv_`` which
    performs a convolution of the Green's functions stored in ``Gk``
    using the provided k‑point mapping arrays.

    Parameters
    ----------
    Gk : (Norb, Norb, Nw, Nk) complex128 Green's function tensor.
    kmap : (Nkall,) int64 Mapping from the full mesh to the reduced mesh.
    invk : (Nk,) int64 Inverse k‑point mapping for interpolation.
    olist : (Nchi,) int64 Orbital indices used in the susceptibility.
    temp : float temperature.
    Nx, Ny, Nz : int grid dimensions for real‑space convolution.

    Returns
    -------
    chi : (Nchi, Nchi, Nw, Nk) complex128 The computed susceptibility tensor.
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
    """Obtain irreducible susceptibility by summing over k‑points.
    The Fortran routine ``get_chi0_sum`` performs a simpler BZ
    summation compared to ``get_chi0_conv`` and is useful when
    convolution maps are not required.
    Parameters
    ----------
    Gk : (Norb, Norb, Nw, Nk) complex128 Green's function tensor.
    invk : (Nkall,) int64 inverse k mapping.
    klist : (Nk,3) float64 list of k‑points.
    olist : (Nchi,) int64 orbital indices.
    temp : float temperature.

    Returns
    -------
    chi : (Nchi, Nchi, Nw, Nk) complex128 irreducible susceptibility tensor.
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
    """Apply a flexible vertex correction without spin–orbit coupling.
    `chi` is modified in place by the Fortran routine but the
    original array is preserved by returning a copy.  The output is
    commonly used as an intermediate quantity in susceptibility
    calculations.
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
           sw_in: bool, sw_sub_sigma: bool = True, scf_loop: int = 300, eps: float = 1.0e-4,
           pp: float = 0.3) -> tuple[np.ndarray, float]:
    """Iteratively mix and compute the self‑energy sigma(k,ω).
    A high‑level wrapper around the Fortran ``mkself`` routine.  The
    mixing parameter ``pp`` is printed at the start to remind the user
    of its value.

    Parameters
    ----------
    Smat, Cmat : float64 matrices used for susceptibility corrections.
    kmap, invk : int64 k‑point mapping arrays.
    olist : int64 orbital list.
    hamk : complex128 k‑space Hamiltonian.
    eig, uni : eigenvalues and eigenvectors used in the calculation.
    mu, fill, temp : float chemical potential, filling, and temperature.
    Nw : int number of frequency points.
    Nx, Ny, Nz : int grid dimensions.
    sw_out, sw_in, sw_sub_sigma : bool switches controlling output, input, and sigma subspace.
    scf_loop : int maximum self‑consistency iterations.
    eps : float convergence tolerance.
    pp : float mixing rate (0<pp<1).
    Returns
    -------
    sigmak : ndarray[complex128] resulting self‑energy.
    mu_self : float final chemical potential computed by the routine.
    """
    print("mixing rate: pp = %3.1f" % pp)
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
        POINTER(c_bool), POINTER(c_bool), POINTER(c_bool)      # sw_sub_sigma, sw_out, sw_in
    ]
    flibs.mkself.restype = c_void_p
    flibs.mkself(sigmak, byref(mu_self), Smat, Cmat, kmap, invk, olist, hamk, eig, uni,
        byref(c_double(mu)), byref(c_double(fill)), byref(c_double(temp)), byref(c_int64(scf_loop)),
        byref(c_double(pp)), byref(c_double(eps)), byref(c_int64(Nkall)), byref(c_int64(Nk)), 
        byref(c_int64(Nw)), byref(c_int64(Norb)), byref(c_int64(Nchi)), byref(c_int64(Nx)), 
        byref(c_int64(Ny)), byref(c_int64(Nz)), byref(c_bool(sw_sub_sigma)), byref(c_bool(sw_out)), byref(c_bool(sw_in)))
    return sigmak, mu_self.value

def mkself_soc(Vmat: np.ndarray, kmap: np.ndarray, invk: np.ndarray, invs: np.ndarray,
           olist: np.ndarray, slist: np.ndarray, hamk: np.ndarray, eig: np.ndarray, uni: np.ndarray, mu: float,
           fill: float, temp: float, Nw: int, Nx: int, Ny: int, Nz: int, sw_out: bool,
           sw_in: bool, sw_sub_sigma: bool = True, scf_loop: int = 300, eps: float = 1.0e-4,
           pp: float = 0.3) -> tuple[np.ndarray, float]:
    print("mixing rate: pp = %3.1f" % pp)
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
        POINTER(c_bool), POINTER(c_bool), POINTER(c_bool)       # sw_sub_sigma, sw_out, sw_in
    ]
    flibs.mkself_soc(sigmak,byref(mu_self),Vmat,kmap,invk,invs,olist,slist,hamk,eig,uni,
               byref(c_double(mu)), byref(c_double(fill)), byref(c_double(temp)), byref(c_int64(scf_loop)),
               byref(c_double(pp)), byref(c_double(eps)), byref(c_int64(Nkall)), byref(c_int64(Nk)), 
               byref(c_int64(Nw)), byref(c_int64(Norb)), byref(c_int64(Nchi)), byref(c_int64(Nx)), 
               byref(c_int64(Ny)), byref(c_int64(Nz)), byref(c_bool(sw_sub_sigma)), byref(c_bool(sw_out)), byref(c_bool(sw_in)))
    return sigmak, mu_self.value

def get_qshift(klist: np.ndarray, qpoint: np.ndarray) -> np.ndarray:
    """Determine integer shifts of k‑mesh induced by a momentum transfer.
    The returned array has length ``Nk`` and indicates, for each k, the
    index of k+q within the mesh.  Useful for evaluating susceptibility
    sums with finite momentum transfer.
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
    """Compute the inverse of :func:`get_qshift`.
    For each k the returned value is the index of k-q; this is useful
    when working with backward momentum transfers.
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
    """Calculate irreducible susceptibility chi0 for given q‑shift.
    This routine wraps the Fortran ``get_chi_irr`` which computes the
    bubble diagram contribution at each Matsubara frequency.  The
    returned array has shape ``(Nw, Nchi, Nchi)``.
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
    """Generate q‑mapped susceptibility and its intermediate array.
    This helper wraps ``chiq_map`` which fills two arrays: ``chis``
    (final susceptibility on a q‑mesh) and ``chi`` (an intermediate storage used internally).
    Both have shape ``(Nx,Ny)``.
    Parameters
    ----------
    uni, eig, ffermi : arrays of eigenvectors, eigenvalues, and occupations used in the calculation.
    klist : (Nk,3) float64 list of k‑points.
    Smat : float64 matrix used in vertex corrections.
    olist : int64 orbital list.
    Nx, Ny : int dimensions of the q‑mesh.
    temp : float temperature.
    ecut : float energy cutoff used in the sum.
    idelta : float broadening parameter.

    Returns
    -------
    chis : (Nx,Ny) complex128 spin susceptibility on q‑mesh.
    chi  : (Nx,Ny) complex128 irreducible susceptibility on q-mesh.
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
    """Compute orbital‑resolved traces of susceptibilities.
    Returns three arrays: ``trchis`` (trace of interacting chi),
    ``trchi0`` (trace of irreducible chi0), and ``chis_orb`` which is an orbital‑indexed version with padding.
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
    """Evaluate irreducible pairing susceptibility (phi) for given q‑shift.
    Returns an array shape ``(Nw, Nchi, Nchi)`` similar to
    :func:`get_chi_irr` but for the particle–particle channel.
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
    """Map pairing susceptibility onto a q‑grid using frequency cutoff.
    This is the counterpart to ``chis_qmap`` but for the pairing channel;
    the returned matrix ``phi`` has dimension ``(Nx,Ny)``.
    The boolean ``sw_omega`` selects whether to include frequency dependence in the calculation.
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
    """Trace pairing susceptibility and produce orbital representation.
    Similar to :func:`get_tr_chi` but for the pairing channel.  The
    returned tuple is ``(trphi, phi_orb)`` where ``phi_orb`` has two
    extra columns of padding.
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
    """Compute corrected susceptibility matrix chis from chi0 and Smat.
    This wraps the Fortran ``get_chis`` which applies interaction
    corrections represented by ``Smat`` to the bare susceptibility
    ``chi0``.  The result has shape ``(Nw, Nchi, Nchi)``.
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
    """Generate spin and charge vertex interaction matrices.
    The underlying Fortran routine ``get_scmat`` fills two real
    matrices ``Smat`` and ``Cmat`` of dimension ``(Nchi, Nchi)`` using
    the orbital index list ``olist`` and site information.  The
    parameters ``U`` and ``J`` control the interaction strengths.

    Parameters
    ----------
    olist : (Nchi,) int64 orbital indices of chi.
    site : (Nchi,) int64 site indices corresponding to each orbital.
    U : float Coulomb repulsion parameter.
    J : float Hund's coupling parameter.

    Returns
    -------
    Smat, Cmat : (Nchi, Nchi) float64 spin and charge interaction matrices.
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
    """Construct spin/charge vertex interaction matrices using orbital‑resolved parameters.
    Each element of the input matrices ``Umat`` and ``Jmat`` corresponds
    to a pair of orbitals; this allows for orbital-dependent interactions.
    The result has the same shape as in :func:`gen_SCmatrix`.
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
    """Evaluate the generalized thermoelectric tensors L11, L12 and L22.
    These tensors appear in linear response formulas and are computed by
    the Fortran routine ``calc_lij`` using band energies, velocities, and Fermi occupations.

    Parameters
    ----------
    eig : (Nk, Norb) float64 eigenvalues.
    vk : (Nk, Norb, 3) complex128 velocities.
    ffermi : (Nk, Norb) float64 occupations.
    mu : float chemical potential.
    w : float frequency.
    idelta : float broadening parameter.
    temp : float temperature.

    Returns
    -------
    L11, L12, L22 : (3,3) complex128 response tensors.
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
    """Compute Kn transport coefficients K0, K1, and K2.
    Parameters mirror those of :func:`calc_Lij` but the outputs are
    real 3×3 matrices (float64) used in thermal and electrical
    conductivity formulae.
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
    """Compute the Hall coeeficient via band integrals."""
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
    """Compute transverse diffusion coefficient tensor `tdf`.
    Output shape is ``(Nw,3,3)`` and the values are real."""
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
    """Estimate relaxation time tau(k) from frequency‑dependent weights.
    The Fortran routine ``get_tau`` solves a simple integral equation
    to produce a k‑dependent relaxation time given input weights
    ``tauw``.  The output has shape ``(Nk, Norb)``.

    Parameters
    ----------
    tauw : (Nw,) float64 frequency‑dependent scattering weights.
    eig : (Nk, Norb) float64 eigenvalues used to determine Nk and Norb.
    tau_max : float upper bound for returned tau values.
    tau_mode : int mode selector passed to the Fortran routine.
    eps : float convergence tolerance.

    Returns
    -------
    tau : (Nk, Norb) float64 obtained relaxation times.
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

def gen_imp_ham(rvec: np.ndarray, ham_r: np.ndarray, ham_i: np.ndarray,
                rlist: np.ndarray, imp_list: np.ndarray,eps: float = 1.0e-5) -> np.ndarray:
    """Construct impurity Hamiltonian from real/imaginary components.
    Parameters
    ----------
    rvec : (Nr,3) float64 real‑space displacement vectors.
    ham_r : (Nr, Norb, Norb) complex128 real‑part hopping blocks.
    ham_i : (Nr, Norb, Norb) complex128 imaginary‑part blocks.
    rlist : (Nsite,) float64 list of impurity positions.
    imp_list : (Nimp,) int64 indices of impurity orbitals.
    eps : float small tolerance used in assembly.

    Returns
    -------
    ham_imp : (Norb*Nsite, Norb*Nsite) complex128 resulting impurity Hamiltonian matrix.
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
    """Perform Fourier transform of impurity Hamiltonian to k‑space.
    Parameters
    ----------
    ham_imp : (Norb*Nsite, Norb*Nsite) complex128 impurity Hamiltonian constructed by :func:`gen_imp_ham`.
    klist : (Nk,3) float64 list of k‑points.
    rlist : (Nsite,) float64 impurity positions.

    Returns
    -------
    ham_k : (Norb*Nk, Norb*Nk) complex128 k‑space representation of the impurity Hamiltonian.
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
    """Compute impurity spectral function over k and ω.
    The returned `spectrum` array has dimensions ``(Nk, Nw)`` and is
    complex; its imaginary part typically carries the physical signal.
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
    flibs.get_spectrum_spagehtti.retype = c_void_p
    flibs.get_spectrum_spagehtti(spectrum, uni, eigs, klist, rlist, wlist, byref(c_int64(Nw)), byref(c_int64(Nk)),
                                 byref(c_int64(Nsite)), byref(c_int64(Norb)), byref(c_double(mu)), byref(c_double(eta)))
    return spectrum

def get_a(inp_data: np.ndarray, xlist: np.ndarray) -> np.ndarray:
    """Compute auxiliary array `a` from input data and grid points.
    A small utility wrapper around the Fortran ``get_a_`` routine.
    """
    Np = len(inp_data)
    a = np.zeros(Np, dtype=np.complex128)
    flibs.get_a_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128), # a
        np.ctypeslib.ndpointer(dtype=np.complex128), # xlist
        np.ctypeslib.ndpointer(dtype=np.complex128), # inpdata
        POINTER(c_int64)                             # Np
    ]
    flibs.get_a_.retype = c_void_p
    flibs.get_a_(a, xlist, inp_data, byref(c_int64(Np)))
    return a

def get_QP(a: np.ndarray, xlist: np.ndarray, wlist: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Compute P and Q arrays from input `a`, grid `xlist`, and frequencies.
    This routine is another small helper wrapping a Fortran call;
    both returned arrays have length ``Nw``.
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
    flibs.get_qp_.retype = c_void_p
    flibs.get_qp_(P, Q, a, xlist, wlist, byref(c_int64(Nw)), byref(c_int64(Np)))
    return Q, P

def pade_with_trace(A: np.ndarray, iwlist: np.ndarray, wlist: np.ndarray) -> np.ndarray:
    """Perform Padé analytic continuation with trace operation.
    The input matrix ``A`` is transformed to produce output ``B`` of
    shape ``(Nk, Nw)``.  This is often used when continuing
    frequency‑dependent quantities from imaginary to real axis.
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
    flibs.pade_with_trace.retype = c_void_p
    flibs.pade_with_trace(A, B, iwlist, wlist, byref(c_int64(Nk)), byref(c_int64(Niw)),
                          byref(c_int64(Nw)), byref(c_int64(Norb)))
    return B

def get_chi0(Smat: np.ndarray, Cmat: np.ndarray, Gk: np.ndarray, olist: np.ndarray,
             kmap: np.ndarray, invk: np.ndarray, temp: float, Nx: int, Ny: int, Nz: int) -> np.ndarray:
    """Evaluate the irreducible susceptibility tensor ``chi0`` on an irreducible q-grid.
    This thin wrapper calls the Fortran ``get_chi0`` routine, which assembles
    a four‑dimensional array of shape ``(Nchi, Nchi, Nw, Nk)``
    representing orbital correlations at each frequency and irreducible k‑point.
    The input spin/charge matrices ``Smat`` and ``Cmat`` are forwarded unchanged
    and are useful for consistency with other routines that take these objects.

    Parameters
    ----------
    Smat, Cmat : (Nchi, Nchi) float64 Spin and charge interaction matrices.
    Gk : (Norb, Norb, Nw, Nkall) complex128 Green's function defined on the full (Nx×Ny×Nz) k grid.
    olist : (Nchi,) int64 Orbital index list for the susceptibility.
    kmap : (Nkall,3) int64 Mapping from full to irreducible k‑points.
    invk : (Nkall,3) int64 Inverse k‑point indices used for Fourier transforms.
    temp : float Temperature (energy units).
    Nx, Ny, Nz : int Dimensions of the full k grid.

    Returns
    -------
    chi : (Nchi, Nchi, Nw, Nk) complex128 Irreducible susceptibility on the reduced q-grid.
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
    flibs.get_chi0.retype = c_void_p
    flibs.get_chi0(chi, Smat, Cmat, Gk, kmap, invk, olist, byref(c_double(temp)),
                   byref(c_int64(Nx)), byref(c_int64(Ny)), byref(c_int64(Nz)),
                   byref(c_int64(Nw)), byref(c_int64(Nk)), byref(c_int64(Nkall)),
                   byref(c_int64(Norb)), byref(c_int64(Nchi)))
    return chi

def linearized_eliashberg(chi: np.ndarray, Gk: np.ndarray, uni: np.ndarray, init_delta: np.ndarray,
                          Smat: np.ndarray, Cmat: np.ndarray, olist: np.ndarray, plist: np.ndarray,
                          kmap: np.ndarray,invk: np.ndarray, Nx: int, Ny: int, Nz: int,
                          temp: float, gap_sym: int, eps: float = 1.0e-5, itemax: int = 300) -> np.ndarray:
    """Solve the linearized Eliashberg equation without SOC.
    The output ``delta`` is the superconducting gap function in orbital representation.
    It has shape (``Norb, Norb, Nw, Nkall``), where ``Nkall`` is the full number of k‑points.
    Convergence is attempted up to ``itemax`` iterations with tolerance ``eps``.

    Parameters
    ----------
    chi : (Nchi, Nchi, Nw, Nk) complex128 Irreducible susceptibility obtained from :func:`get_chi0`.
    Gk : (Norb, Norb, Nw, Nkall) complex128 Green's function on full k grid.
    uni : (Norb, Norb, Nk) complex128 Unitary transformation matrices.
    init_delta : (Norb, Norb, Nw, Nkall) float64 Initial guess for the gap function.
    Smat, Cmat : (Nchi, Nchi) float64 Interaction matrices.
    olist : (Nchi,) int64 Orbital indices.
    kmap, invk : arrays k‑point mapping information.
    Nx, Ny, Nz : int Grid dimensions for kmap/invk.
    temp : float Temperature.
    gap_sym : int Symmetry index for the gap function.
    eps : float Convergence tolerance.
    itemax : int Maximum iteration count.

    Returns
    -------
    delta : (Norb, Norb, Nw, Nkall) complex128 Linearized gap function solution.
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
        POINTER(c_int64), POINTER(c_int64)            # itemax, gapsym
    ]
    flibs.lin_eliash.retype = c_void_p
    flibs.lin_eliash(delta, chi, Gk, uni, init_delta, Smat, Cmat, olist, plist, kmap, invk,
                     byref(c_double(temp)), byref(c_double(eps)), byref(c_int64(Nkall)),
                     byref(c_int64(Nk)), byref(c_int64(Nw)), byref(c_int64(Nchi)),
                     byref(c_int64(Norb)), byref(c_int64(Nx)), byref(c_int64(Ny)),
                     byref(c_int64(Nz)), byref(c_int64(itemax)), byref(c_int64(gap_sym)))
    return delta

def linearized_eliashberg_soc(chi: np.ndarray, Gk: np.ndarray, uni: np.ndarray, init_delta: np.ndarray,
                              Vmat: np.ndarray, sgnsig: np.ndarray, sgnsig2: np.ndarray, plist: np.ndarray,
                              slist: np.ndarray, olist: np.ndarray, kmap: np.ndarray, invk: np.ndarray,
                              invs: np.ndarray, invschi: np.ndarray, Nx: int, Ny: int, Nz: int, temp: float,
                              gap_sym: int, eps: float = 1.0e-4, itemax: int = 300) -> np.ndarray:
    """Linearized Eliashberg solver including spin–orbit coupling effects.
    This more general version accepts additional sign and orbital mapping arrays associated 
    with SOC (`sgnsig`, `sgnsig2`, `invs`, `invschi`).
    The output ``delta`` has the same shape as in:func:`linearized_eliashberg`.
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
        POINTER(c_int64), POINTER(c_int64)            # itemax, gapsym
    ]
    flibs.lin_eliash_soc.retype = c_void_p
    flibs.lin_eliash_soc(delta, chi, Gk, uni, init_delta, Vmat, sgnsig, sgnsig2, plist, olist,
                         slist, kmap, invk, invs, invschi, byref(c_double(temp)),
                         byref(c_double(eps)), byref(c_int64(Nkall)), byref(c_int64(Nk)),
                         byref(c_int64(Nw)), byref(c_int64(Nchi)), byref(c_int64(Norb)),
                         byref(c_int64(Nx)), byref(c_int64(Ny)), byref(c_int64(Nz)),
                         byref(c_int64(itemax)), byref(c_int64(gap_sym)))
    return delta

def gen_Vmatrix(olist: np.ndarray, slist: np.ndarray, site: np.ndarray, invs: np.ndarray,
                U: float, J: float) -> np.ndarray:
    """Generate spin–orbit coupling interaction matrix Vmat in SOC calculation.
    The returned matrix has shape ``(Nchi, Nchi)`` and is constructed by the Fortran routine 
    ``get_vmat_soc`` using orbital lists and site/sign information.

    Parameters
    ----------
    olist : (Nchi,) int64 Orbital indices.
    slist : (Norb,) int64 Spin‑orbit label for each orbital.
    site : (Nchi,) int64 Site indices corresponding to each orbital pair.
    invs : (Norb,) int64 Inverse mapping for spin indices.
    U : float Hubbard interaction parameter.
    J : float Hund's coupling.

    Returns
    -------
    Vmat : (Nchi, Nchi) float64
        Interaction matrix suitable for SOC calculations.
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
    flibs.get_vmat_soc.restype = c_void_p
    flibs.get_vmat_soc(Vmat, olist, slist, site, invs, byref(c_double(U)),
                       byref(c_double(J)), byref(c_int64(Nchi)), byref(c_int64(Norb)))
    return Vmat

def gen_Vmatrix_orb(olist: np.ndarray, slist: np.ndarray, site: np.ndarray, invs: np.ndarray,
                    Umat: np.ndarray, Jmat: np.ndarray,) -> np.ndarray:
    """Construct SOC interaction `Vmat` using orbital‑resolved parameters.
    Each entry of ``Umat`` and ``Jmat`` corresponds to a pair of orbitals,
    allowing for orbital‑dependent interactions.
    The resulting matrix has the same shape as generated by:func:`gen_Vmatrix`.
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
    flibs.get_vmat_soc_orb.restype = c_void_p
    flibs.get_vmat_soc_orb(Vmat, olist, slist, site, invs, Umat, Jmat,
                           byref(c_int64(Nchi)), byref(c_int64(Norb)))
    return Vmat

def get_chi0_soc(Vmat: np.ndarray, Gk: np.ndarray, olist: np.ndarray, slist: np.ndarray,
                 kmap: np.ndarray, invk: np.ndarray, invs: np.ndarray, temp: float,
                 Nx: int, Ny: int, Nz: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Calculate SOC irreducible susceptibility and auxiliary sign matrices.
    The routine returns a tuple ``(chi, sgnsig, sgnsig2, invschi)`` where
    ``chi`` has shape ``(Nchi, Nchi, Nw, Nk)`` and the others are
    auxiliary arrays needed for later SOC-specific operations.
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
    flibs.get_chi0_soc.restype = c_void_p
    flibs.get_chi0_soc(chi, sgnsig, sgnsig2, invschi, Vmat, Gk, kmap, invk, invs, olist, slist,
                       byref(c_double(temp)), byref(c_int64(Nx)), byref(c_int64(Ny)),
                       byref(c_int64(Nz)), byref(c_int64(Nw)), byref(c_int64(Nk)),
                       byref(c_int64(Nkall)), byref(c_int64(Nchi)), byref(c_int64(Norb)))
    return chi, sgnsig, sgnsig2, invschi

def get_chis_chic_soc(chi: np.ndarray, Vmat: np.ndarray, olist: np.ndarray, slist: np.ndarray,
                      invs: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transform susceptibility to charge/chi channels in SOC case.
    The helper ``get_orb_list`` builds a mapping from orbital pairs to
    chi indices; the function then calls the Fortran routine to produce
    three reduced tensors ``chic``, ``chiszz`` and ``chispm``.
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
    flibs.get_chis_chic_soc.retype = c_void_p
    flibs.get_chis_chic_soc(chic, chiszz, chispm, chi, Vmat, orb_list, olist, slist, invs, byref(c_int64(Nk)),
                            byref(c_int64(Nw)), byref(c_int64(Nchi)), byref(c_int64(Norb)))
    return chic, chiszz, chispm

def get_chis_chic(chi: np.ndarray, Smat: np.ndarray, Cmat: np.ndarray,) -> tuple[np.ndarray, np.ndarray]:
    """Compute corrected susceptibilities ``chis`` and ``chic`` from ``chi``.
    The routine applies spin/charge vertex matrices ``Smat`` and ``Cmat`` to
    the irreducible susceptibility ``chi`` to produce two new tensors with dimensions ``(Nchi, Nchi, Nk)``.
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
    flibs.get_chis_chic.retype = c_void_p
    flibs.get_chis_chic(chis, chic, chi, Smat, Cmat, byref(c_int64(Nk)),
                        byref(c_int64(Nw)), byref(c_int64(Nchi)))
    return chis, chic

def conv_delta_orb_to_band(delta: np.ndarray, uni: np.ndarray, invk: np.ndarray, plist: np.ndarray,gap_sym) -> np.ndarray:
    """Convert orbital‑space gap function ``delta`` to band representation.
    Parameters
    ----------
    delta : (Norb, Norb, Nw, Nkall) complex128 Gap function in orbital basis.
    uni : (Norb, Norb, Nk) complex128 Unitary transformation matrices.
    invk : (Nkall,3) int64 Inverse k‑point mapping.

    Returns
    -------
    deltab : (Norb, Norb, Nkall) complex128 Gap function in band basis.
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
    flibs.conv_delta_orb_to_band.retype = c_void_p
    flibs.conv_delta_orb_to_band(deltab, delta, uni, plist, invk, byref(c_int64(Norb)),
                                 byref(c_int64(Nkall)), byref(c_int64(Nk)), byref(c_int64(Nw)), byref(c_int64(gap_sym)))
    return deltab

def conv_delta_orb_to_band_soc(delta: np.ndarray, uni: np.ndarray, invk: np.ndarray,
                               invs: np.ndarray,slist: np.ndarray) -> np.ndarray:
    """SOC-aware conversion of orbital gap to band representation.
    Same as :func:`conv_delta_orb_to_band` but includes additional
    spin indices ``invs`` and ``slist`` required by the SOC routine.
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
    flibs.conv_delta_orb_to_band_soc.retype = c_void_p
    flibs.conv_delta_orb_to_band_soc(deltab, delta, uni, invk, invs, slist, byref(c_int64(Norb)),
                                     byref(c_int64(Nkall)), byref(c_int64(Nk)), byref(c_int64(Nw)))
    return deltab

def gen_Fk(Gk: np.ndarray, delta: np.ndarray, invk: np.ndarray) -> np.ndarray:
    """Generate Fk tensor from Green's function and gap function.
    The result Fk has same dimensions as ``delta`` and is used in
    post‑processing of superconducting properties.
    """
    Nkall, Nk, Nw, Norb = len(invk), len(Gk[0, 0, 0]), len(delta[0, 0]), len(delta)
    Fk = np.zeros((Norb, Norb, Nw, Nk), dtype=np.complex128)
    flibs.mkfk_trs_nsoc_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128),         # Fk
        np.ctypeslib.ndpointer(dtype=np.complex128),         # Gk
        np.ctypeslib.ndpointer(dtype=np.complex128),         # delta
        POINTER(c_int64), POINTER(c_int64), POINTER(c_int64) # Nk, Nw, Norb
    ]
    flibs.mkfk_trs_nsoc_.retype = c_void_p
    flibs.mkfk_trs_nsoc_(Fk, Gk, delta, byref(c_int64(Nk)), byref(c_int64(Nw)), byref(c_int64(Norb)))

    return Fk

def gen_Fk_soc(Gk: np.ndarray, delta: np.ndarray, invk: np.ndarray, invs: np.ndarray,
               slist: np.ndarray,gap_sym: int) -> np.ndarray:
    """SOC-aware generation of Fk tensor from Green's function and gap function.
    Same as :func:`gen_Fk` but includes additional spin indices for SOC case.
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
    flibs.mkfk_trs_soc_.retype = c_void_p
    flibs.mkfk_trs_soc_(Fk, Gk, delta, sgnsig, slist, invk, invs, byref(c_int64(Nkall)),
                        byref(c_int64(Nk)), byref(c_int64(Nw)), byref(c_int64(Norb)), byref(c_int64(gap_sym)))
    return Fk

def remap_gap(delta0,plist,invk,gap_sym):
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
    flibs.remap_delta_.retype = c_void_p
    flibs.remap_delta_(delta, delta0, plist, invk, byref(c_int64(Nkall)), byref(c_int64(Nk)),
                       byref(c_int64(Nw)), byref(c_int64(Norb)), byref(c_int64(gap_sym)))
    return delta

def get_eig_or_tr_chi(chi: np.ndarray, invk: np.ndarray, sw_eig:bool) -> np.ndarray:
    Nkall,Nk=len(invk),len(chi.T)
    Nchi=int(np.sqrt(chi.size/Nk))
    chiq=np.zeros((Nchi,Nchi,Nkall),dtype=np.complex128)
    flibs.get_eig_or_tr_chi.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.complex128), # chiq
        np.ctypeslib.ndpointer(dtype=np.complex128), # chi
        np.ctypeslib.ndpointer(dtype=np.int64),      # invk
        POINTER(c_int64), POINTER(c_int64),          # Nkall, Nk
        POINTER(c_int64), POINTER(c_bool)            # Nchi, sw_eig
    ]
    flibs.get_eig_or_tr_chi.retype = c_void_p
    flibs.get_eig_or_tr_chi(chiq,chi, invk, byref(c_int64(Nkall)), byref(c_int64(Nk)),
                             byref(c_int64(Nchi)), byref(c_bool(sw_eig)))
    return chiq

def gen_irr_k_TRS(Nx: int, Ny: int, Nz: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Construct irreducible k‑point list for time‑reversal symmetric lattice.
    A fairly convoluted set of rules determine ``Nk`` (the number of irreducible k‑points)
    based on even/odd parity of the grid dimensions.
    The routine returns the k‑point coordinates, a full grid mapping ``kmap``,
    and an inverse list ``invk_ft_list`` used by Fortran routines.
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
    flibs.generate_irr_kpoint_inv.retype = c_void_p
    flibs.generate_irr_kpoint_inv(klist, kmap, invk_ft_list, byref(c_int64(Nk)),
                                  byref(c_int64(Nx)), byref(c_int64(Ny)), byref(c_int64(Nz)))
    return klist, kmap, invk_ft_list

def gen_kpoint_weight(invk,Nk):
    Nkall=len(invk)
    weight=np.zeros(Nk,dtype=np.float64)
    flibs.get_kweight.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64), # invk
        np.ctypeslib.ndpointer(dtype=np.int64),   # weight
        POINTER(c_int64), POINTER(c_int64)        # Nk, Nkall
    ]
    flibs.get_kweight.retype = c_void_p
    flibs.get_kweight(weight,invk,byref(c_int64(Nk)),byref(c_int64(Nkall)))
    return weight

def get_plist(rvec,ham_r):
    Nr=len(rvec)
    Norb=int(np.sqrt(ham_r.size/Nr))
    Pmn=np.zeros((Norb,Norb),dtype=np.float64)
    flibs.get_parity_prop.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),    #Pmn
        np.ctypeslib.ndpointer(dtype=np.float64),    #rvec
        np.ctypeslib.ndpointer(dtype=np.complex128), #ham_r
        POINTER(c_int64), POINTER(c_int64)           #Norb,Nr
    ]
    flibs.get_parity_prop.retype = c_void_p
    flibs.get_parity_prop(Pmn,rvec,ham_r,byref(c_int64(Norb)),byref(c_int64(Nr)))
    return np.sign(Pmn[0,:])
