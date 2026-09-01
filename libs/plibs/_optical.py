"""
Optical constants derived from the dielectric function, and the spectra that
carry the colour of a non-metal.

For a metal the colour lives in the specular reflectivity R(w). For a non-metal
R is small (R ~ 0.04-0.08 for n ~ 1.5) and almost colourless, so the colour has
to come from the light that has travelled through the material:

  * a transparent slab of thickness d transmits
        T = (1-R)^2 exp(-alpha d) / (1 - R^2 exp(-2 alpha d))
    (both surfaces, incoherent multiple reflections). Note that this colour is
    NOT an intrinsic material property: it depends on d.
  * a thick scattering powder (pigment) shows the Kubelka-Munk diffuse
    reflectance
        R_inf = 1 + K/S - sqrt((K/S)^2 + 2 K/S),   K/S = alpha * l
    where the effective scattering length l is the single free parameter (it
    absorbs the convention factor between the absorption coefficient alpha and
    the Kubelka-Munk absorption K, as well as the particle size).

The absorptance A = 1 - R - T is the complement of what leaves the sample, so a
colour computed from A comes out as the complementary colour and has to be
negated (see plibs._color, key 'hex_neg').
"""
import numpy as np
import scipy.constants as scconst

HBAR_EVS = scconst.physical_constants['reduced Planck constant in eV s'][0]
C_LIGHT = scconst.c             # speed of light [m/s]


def refractive_index(eps: np.ndarray) -> np.ndarray:
    """
    @fn refractive_index
    @brief Complex refractive index N = n + i*kappa from the dielectric function.
    @param eps: complex dielectric function, any shape
    @retval     N on the passive branch (Im N >= 0), same shape as eps
    """
    N = np.sqrt(np.asarray(eps, dtype=complex))
    return np.where(N.imag < 0.0, -N, N)


def reflectivity(N: np.ndarray) -> np.ndarray:
    """
    @fn reflectivity
    @brief Normal-incidence (vacuum/sample) Fresnel reflectivity from N.
    @param N: complex refractive index
    @retval    R = |(N-1)/(N+1)|^2, real, in [0, 1]
    """
    return (np.abs((N - 1.0)/(N + 1.0))**2).real


def absorption_coefficient(w: np.ndarray, N: np.ndarray) -> np.ndarray:
    """
    @fn absorption_coefficient
    @brief Beer-Lambert absorption coefficient alpha = 2*w*kappa/(hbar*c).
    @param w: photon energies [eV]; broadcast against the leading axes of N
    @param N: complex refractive index (kappa = Im N)
    @retval    alpha [1/m], same shape as N
    """
    N = np.asarray(N)
    w = np.asarray(w, dtype=float)
    while w.ndim < N.ndim:
        w = w[..., None]
    return 2.0*w*N.imag/(HBAR_EVS*C_LIGHT)


def transmittance(R: np.ndarray, alpha: np.ndarray, d: float) -> np.ndarray:
    """
    @fn transmittance
    @brief Slab transmittance including both surfaces and incoherent multiple
           internal reflections: T = (1-R)^2 e^-ad / (1 - R^2 e^-2ad).
    @param     R: normal-incidence reflectivity
    @param alpha: absorption coefficient [1/m]
    @param     d: sample thickness [m]
    @retval       T in [0, 1]; for alpha = 0 it reduces to (1-R)/(1+R)
    """
    e = np.exp(-np.clip(np.asarray(alpha)*d, 0.0, 700.0))   # no overflow to 0*inf
    den = 1.0 - (R*e)**2
    # den vanishes only for a lossless perfect mirror (R = 1, alpha = 0), where
    # the numerator vanishes too and the physical answer is T = 0
    return np.divide((1.0 - R)**2*e, den, out=np.zeros_like(den), where=den != 0.0)


def kubelka_munk(alpha: np.ndarray, ell: float) -> np.ndarray:
    """
    @fn kubelka_munk
    @brief Diffuse reflectance of an optically thick scattering layer (powder,
           pigment) in the Kubelka-Munk two-flux model, with K/S = alpha*ell.
    @param alpha: absorption coefficient [1/m]
    @param   ell: effective scattering length [m] (the model's free parameter)
    @retval       R_inf in [0, 1]: 1 for a non-absorbing powder, 0 for a black one
    """
    ks = np.clip(np.asarray(alpha)*ell, 0.0, None)
    return 1.0 + ks - np.sqrt(ks*(ks + 2.0))
