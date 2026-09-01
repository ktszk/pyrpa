"""
Colorimetry: the perceived colour of a metal from its computed reflectivity.

Follows the standard reflectance-colorimetry route used for first-principles
metal colours (G. Prandini, G.-M. Rignanese, N. Marzari, "Photorealistic
modelling of metals from first principles", npj Comput. Mater. 5, 129 (2019)):

    R(w)  ->  R(lambda)  ->  CIE XYZ tristimulus  ->  sRGB

    X = int S(l) R(l) xbar(l) dl / N ,  N = int S(l) ybar(l) dl   (Y, Z alike)

with S = CIE standard illuminant D65 and (xbar, ybar, zbar) the CIE 1931 2-deg
standard-observer colour matching functions. Normalising by N makes a perfect
mirror (R == 1) map exactly onto the D65 white point, i.e. the illuminant is
white-balanced away and what is left is the intrinsic colour of the metal --
the same convention as in the reference above. XYZ is then converted with the
sRGB (IEC 61966-2-1, D65) primaries and transfer function.

Reference data below: CIE 1931 2-deg colour matching functions and the D65
relative spectral power distribution, 380-780 nm in 5 nm steps (CVRL database).
Truncating to that window reproduces the tabulated D65 white point
x = 0.31272, y = 0.32903 to five decimals.
"""
import numpy as np
import scipy.constants as scconst

# lambda [nm], xbar, ybar, zbar (CIE 1931 2-deg), S_D65 (relative SPD)
_CIE = np.array([
    [380, 0.00136800, 0.00003900, 0.00645000,   49.9755],
    [385, 0.00223600, 0.00006400, 0.01054999,   52.3118],
    [390, 0.00424300, 0.00012000, 0.02005001,   54.6482],
    [395, 0.00765000, 0.00021700, 0.03621000,   68.7015],
    [400, 0.01431000, 0.00039600, 0.06785001,   82.7549],
    [405, 0.02319000, 0.00064000, 0.11020000,   87.1204],
    [410, 0.04351000, 0.00121000, 0.20740000,   91.4860],
    [415, 0.07763000, 0.00218000, 0.37130000,   92.4589],
    [420, 0.13438000, 0.00400000, 0.64560000,   93.4318],
    [425, 0.21477000, 0.00730000, 1.03905010,   90.0570],
    [430, 0.28390000, 0.01160000, 1.38560000,   86.6823],
    [435, 0.32850000, 0.01684000, 1.62296000,   95.7736],
    [440, 0.34828000, 0.02300000, 1.74706000,  104.8650],
    [445, 0.34806000, 0.02980000, 1.78260000,  110.9360],
    [450, 0.33620000, 0.03800000, 1.77211000,  117.0080],
    [455, 0.31870000, 0.04800000, 1.74410000,  117.4100],
    [460, 0.29080000, 0.06000000, 1.66920000,  117.8120],
    [465, 0.25110000, 0.07390000, 1.52810000,  116.3360],
    [470, 0.19536000, 0.09098000, 1.28764000,  114.8610],
    [475, 0.14210000, 0.11260000, 1.04190000,  115.3920],
    [480, 0.09564000, 0.13902000, 0.81295010,  115.9230],
    [485, 0.05795001, 0.16930000, 0.61620000,  112.3670],
    [490, 0.03201000, 0.20802000, 0.46518000,  108.8110],
    [495, 0.01470000, 0.25860000, 0.35330000,  109.0820],
    [500, 0.00490000, 0.32300000, 0.27200000,  109.3540],
    [505, 0.00240000, 0.40730000, 0.21230000,  108.5780],
    [510, 0.00930000, 0.50300000, 0.15820000,  107.8020],
    [515, 0.02910000, 0.60820000, 0.11170000,  106.2960],
    [520, 0.06327000, 0.71000000, 0.07824999,  104.7900],
    [525, 0.10960000, 0.79320000, 0.05725001,  106.2390],
    [530, 0.16550000, 0.86200000, 0.04216000,  107.6890],
    [535, 0.22574990, 0.91485010, 0.02984000,  106.0470],
    [540, 0.29040000, 0.95400000, 0.02030000,  104.4050],
    [545, 0.35970000, 0.98030000, 0.01340000,  104.2250],
    [550, 0.43344990, 0.99495010, 0.00875000,  104.0460],
    [555, 0.51205010, 1.00000000, 0.00575000,  102.0230],
    [560, 0.59450000, 0.99500000, 0.00390000,  100.0000],
    [565, 0.67840000, 0.97860000, 0.00275000,   98.1671],
    [570, 0.76210000, 0.95200000, 0.00210000,   96.3342],
    [575, 0.84250000, 0.91540000, 0.00180000,   96.0611],
    [580, 0.91630000, 0.87000000, 0.00165000,   95.7880],
    [585, 0.97860000, 0.81630000, 0.00140000,   92.2368],
    [590, 1.02630000, 0.75700000, 0.00110000,   88.6856],
    [595, 1.05670000, 0.69490000, 0.00100000,   89.3459],
    [600, 1.06220000, 0.63100000, 0.00080000,   90.0062],
    [605, 1.04560000, 0.56680000, 0.00060000,   89.8026],
    [610, 1.00260000, 0.50300000, 0.00034000,   89.5991],
    [615, 0.93840000, 0.44120000, 0.00024000,   88.6489],
    [620, 0.85444990, 0.38100000, 0.00019000,   87.6987],
    [625, 0.75140000, 0.32100000, 0.00010000,   85.4936],
    [630, 0.64240000, 0.26500000, 0.00005000,   83.2886],
    [635, 0.54190000, 0.21700000, 0.00003000,   83.4939],
    [640, 0.44790000, 0.17500000, 0.00002000,   83.6992],
    [645, 0.36080000, 0.13820000, 0.00001000,   81.8630],
    [650, 0.28350000, 0.10700000, 0.00000000,   80.0268],
    [655, 0.21870000, 0.08160000, 0.00000000,   80.1207],
    [660, 0.16490000, 0.06100000, 0.00000000,   80.2146],
    [665, 0.12120000, 0.04458000, 0.00000000,   81.2462],
    [670, 0.08740000, 0.03200000, 0.00000000,   82.2778],
    [675, 0.06360000, 0.02320000, 0.00000000,   80.2810],
    [680, 0.04677000, 0.01700000, 0.00000000,   78.2842],
    [685, 0.03290000, 0.01192000, 0.00000000,   74.0027],
    [690, 0.02270000, 0.00821000, 0.00000000,   69.7213],
    [695, 0.01584000, 0.00572300, 0.00000000,   70.6652],
    [700, 0.01135916, 0.00410200, 0.00000000,   71.6091],
    [705, 0.00811092, 0.00292900, 0.00000000,   72.9790],
    [710, 0.00579035, 0.00209100, 0.00000000,   74.3490],
    [715, 0.00410946, 0.00148400, 0.00000000,   67.9765],
    [720, 0.00289933, 0.00104700, 0.00000000,   61.6040],
    [725, 0.00204919, 0.00074000, 0.00000000,   65.7448],
    [730, 0.00143997, 0.00052000, 0.00000000,   69.8856],
    [735, 0.00099995, 0.00036110, 0.00000000,   72.4863],
    [740, 0.00069008, 0.00024920, 0.00000000,   75.0870],
    [745, 0.00047602, 0.00017190, 0.00000000,   69.3398],
    [750, 0.00033230, 0.00012000, 0.00000000,   63.5927],
    [755, 0.00023483, 0.00008480, 0.00000000,   55.0054],
    [760, 0.00016615, 0.00006000, 0.00000000,   46.4182],
    [765, 0.00011741, 0.00004240, 0.00000000,   56.6118],
    [770, 0.00008308, 0.00003000, 0.00000000,   66.8054],
    [775, 0.00005871, 0.00002120, 0.00000000,   65.0941],
    [780, 0.00004151, 0.00001499, 0.00000000,   63.3828],
])

CIE_LAM   = _CIE[:, 0]          # wavelength grid [nm]
CIE_CMF   = _CIE[:, 1:4]        # [81, 3] colour matching functions
CIE_D65   = _CIE[:, 4]          # [81] D65 relative spectral power

# eV <-> nm: lambda[nm] = (h/e)*c / E[eV]; h and e are exact SI definitions
EV_NM = scconst.h/scconst.e*scconst.c*1.e9

# XYZ -> linear sRGB (IEC 61966-2-1 primaries, D65 white point)
XYZ2RGB = np.array([[ 3.2404542, -1.5371385, -0.4985314],
                    [-0.9692660,  1.8760108,  0.0415560],
                    [ 0.0556434, -0.2040259,  1.0572252]])

# White point of the tabulated illuminant (i.e. the XYZ of a perfect mirror,
# R == 1) and its linear-RGB image. The 5 nm / 380-780 nm truncation shifts the
# white point from the nominal (0.95047, 1, 1.08883) in the fifth decimal, so
# the linear channels are divided by _WHITE_LIN: a perfect mirror then maps onto
# exactly (1, 1, 1) and the gamut-clipping flag stays free of round-off.
XYZ_WHITE = (CIE_D65[:, None]*CIE_CMF).sum(axis=0)/(CIE_D65*CIE_CMF[:, 1]).sum()
_WHITE_LIN = XYZ2RGB.dot(XYZ_WHITE)


def srgb_gamma(c: np.ndarray) -> np.ndarray:
    """
    @fn srgb_gamma
    @brief Apply the sRGB opto-electronic transfer function (linear -> encoded).
    @param  c: linear-light RGB in [0, 1]
    @retval    gamma-encoded sRGB in [0, 1]
    """
    c = np.clip(c, 0.0, 1.0)
    return np.where(c <= 0.0031308, 12.92*c, 1.055*np.power(c, 1.0/2.4) - 0.055)


def reflectivity_to_xyz(w: np.ndarray, refl: np.ndarray) -> tuple[np.ndarray, bool]:
    """
    @fn reflectivity_to_xyz
    @brief Integrate a reflectivity spectrum against the CIE observer under D65.
    @param     w: photon energies [eV], strictly positive, any order
    @param  refl: the spectral factor on w, values in [0, 1]
    @retval   XYZ: CIE 1931 tristimulus values (Y = 1 for a perfect mirror)
    @retval  full: True if the input spectrum spans the whole 380-780 nm window;
                   when False, R is held constant beyond the computed range
                   (edge clamping) and the colour is only indicative.
    """
    w = np.asarray(w, dtype=float)
    refl = np.asarray(refl, dtype=float)
    pos = w > 0.0
    lam = EV_NM/w[pos]                       # eV -> nm (descending in w)
    o = np.argsort(lam)                      # np.interp needs ascending x
    lam, r = lam[o], refl[pos][o]
    full = bool(lam[0] <= CIE_LAM[0] and lam[-1] >= CIE_LAM[-1])
    rl = np.interp(CIE_LAM, lam, r)          # clamps to edge values outside range
    sw = CIE_D65*rl
    XYZ = (sw[:, None]*CIE_CMF).sum(axis=0)/(CIE_D65*CIE_CMF[:, 1]).sum()
    return XYZ, full


def xyz_to_srgb(XYZ: np.ndarray) -> tuple[np.ndarray, bool]:
    """
    @fn xyz_to_srgb
    @brief Convert CIE XYZ to gamma-encoded sRGB, reporting gamut clipping.
    @param  XYZ: CIE 1931 tristimulus values
    @retval rgb: sRGB components in [0, 1]
    @retval clipped: True if the colour fell outside the sRGB gamut (or above
                     white) and had to be clipped
    """
    lin = XYZ2RGB.dot(XYZ)/_WHITE_LIN
    clipped = bool(np.any(lin < -1.e-6) or np.any(lin > 1.0 + 1.e-6))
    return srgb_gamma(lin), clipped


def spectrum_color(w: np.ndarray, spec: np.ndarray) -> dict:
    """
    @fn spectrum_color
    @brief Perceived colour of a spectrum that leaves a sample towards the eye.
    @param     w: photon energies [eV], the w = 0 point may be included
    @param  spec: the spectral factor on w, in [0, 1] -- the reflectivity R(w)
                  of a metal, the transmittance T(w) of a slab, or the diffuse
                  reflectance R_inf(w) of a powder
    @retval  dict with keys
             'XYZ'  : CIE 1931 tristimulus values [3]
             'xy'   : CIE 1931 chromaticity coordinates [2]
             'Y'    : luminance factor (1 = perfect mirror / clear sample)
             'rgb'  : gamma-encoded sRGB in [0, 1] [3]
             'rgb255': the same as 8-bit integers [3]
             'hex'  : '#rrggbb' string
             'rgb_neg', 'rgb255_neg', 'hex_neg' : the photographic negative
                  (1 - rgb in encoded sRGB). Feed an absorptance A(w) and read
                  these to get the colour of the light that survives -- what a
                  colour computed from A alone would otherwise report as its
                  own complement.
             'full_range': False if the 380-780 nm window is not fully covered
             'clipped'   : True if the colour was clipped into the sRGB gamut
    """
    XYZ, full = reflectivity_to_xyz(w, spec)
    rgb, clipped = xyz_to_srgb(XYZ)
    s = XYZ.sum()
    xy = XYZ[:2]/s if s > 0.0 else np.zeros(2)
    r255 = np.clip(np.rint(255*rgb), 0, 255).astype(int)
    neg = 1.0 - rgb
    n255 = np.clip(np.rint(255*neg), 0, 255).astype(int)
    return {'XYZ': XYZ, 'xy': xy, 'Y': float(XYZ[1]), 'rgb': rgb, 'rgb255': r255,
            'hex': '#%02x%02x%02x' % tuple(r255),
            'rgb_neg': neg, 'rgb255_neg': n255,
            'hex_neg': '#%02x%02x%02x' % tuple(n255),
            'full_range': full, 'clipped': clipped}


# the metal case, kept under its original name
reflectivity_color = spectrum_color
