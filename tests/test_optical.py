"""
Optical constants, the colour of metals (reflectivity -> CIE XYZ -> sRGB) and the
non-metal colour routes (libs/plibs/_optical.py, _color.py; option=6 of main.py).

Runs standalone (no pytest needed):  python tests/test_optical.py
Also works under pytest if installed:  pytest tests/test_optical.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import libs.plibs._color as _color
import libs.plibs._optical as _optical


# --------------------------------------------------------------------------- #
#  Reflectivity -> colour conversion
# --------------------------------------------------------------------------- #
# Physics / reference tests for the reflectivity -> colour conversion
# (libs/plibs/_color.py).
#
# The chain is  R(w) -> R(lambda) -> CIE XYZ (2-deg observer, illuminant D65)
# -> sRGB, the standard route for first-principles metal colours (Prandini et
# al., npj Comput. Mater. 5, 129 (2019)). It is pinned here by quantities whose
# values are fixed by the CIE definitions, so a corrupted colour-matching table
# or a wrong normalisation shows up immediately:
#
#   * a perfect mirror (R == 1) must land exactly on the D65 white point
#     x = 0.31272, y = 0.32903, Y = 1, i.e. sRGB #ffffff
#   * a grey reflector R = c keeps the white chromaticity and gives
#     Y = c and the sRGB transfer function value of c
#   * spectrally selective reflectors must produce the right hue (a 550 nm band
#     -> green near the spectral locus, a >600 nm long-pass -> red, ...)
#   * Lorentz-Drude models of real metals (Rakic et al., Appl. Opt. 37, 5271
#     (1998)) must come out in the right order: Au/Cu warm (x > x_D65), Ag/Al
#     nearly neutral, and every metal bright (Y > 0.5)


# a dense grid covering the whole visible window with room to spare
WGRID = np.linspace(0.4, 6.5, 2000)


def _lorentz_drude(w, wp, osc):
    """eps(w) of the Rakic Lorentz-Drude model: osc[0] is the Drude term."""
    f0, G0, _ = osc[0]
    eps = 1. - f0*wp**2/(w*(w + 1j*G0))
    for f, G, w0 in osc[1:]:
        eps += f*wp**2/((w0**2 - w**2) - 1j*G*w)
    return eps


def _reflectivity(eps):
    N = np.sqrt(eps)
    N = np.where(N.imag < 0., -N, N)
    return (np.abs((N - 1.)/(N + 1.))**2).real


def test_tables_are_intact():
    """The embedded CIE tables must be the 5 nm / 380-780 nm window and
    reproduce the tabulated D65 white point (0.95047, 1, 1.08883)."""
    assert _color.CIE_LAM[0] == 380. and _color.CIE_LAM[-1] == 780. and len(_color.CIE_LAM) == 81
    assert np.allclose(np.diff(_color.CIE_LAM), 5.)
    assert _color.CIE_CMF.shape == (81, 3) and _color.CIE_D65.shape == (81,)
    assert np.allclose(_color.XYZ_WHITE, [0.95047, 1.0, 1.08883], atol=5e-4)
    x, y = _color.XYZ_WHITE[:2]/_color.XYZ_WHITE.sum()
    assert abs(x - 0.31272) < 1e-4 and abs(y - 0.32903) < 1e-4


def test_perfect_mirror_is_white():
    """R == 1 defines the white point: pure white, unit luminance, no clipping."""
    c = _color.reflectivity_color(WGRID, np.ones_like(WGRID))
    assert c['hex'] == '#ffffff'
    assert abs(c['Y'] - 1.0) < 1e-9
    assert abs(c['xy'][0] - 0.31272) < 1e-4 and abs(c['xy'][1] - 0.32903) < 1e-4
    assert c['full_range'] and not c['clipped']


def test_grey_reflector_keeps_chromaticity():
    """A flat R = c is neutral: white chromaticity, Y = c, sRGB = gamma(c)."""
    for c0 in (0.2, 0.5, 0.85):
        c = _color.reflectivity_color(WGRID, c0*np.ones_like(WGRID))
        assert abs(c['Y'] - c0) < 1e-9
        assert np.allclose(c['xy'], _color.XYZ_WHITE[:2]/_color.XYZ_WHITE.sum(), atol=1e-9)
        assert np.allclose(c['rgb'], _color.srgb_gamma(np.full(3, c0)), atol=1e-9)


def test_selective_reflectors_have_the_right_hue():
    """Narrow-band and step reflectors must map onto the expected hues."""
    lam = _color.ev_nm/WGRID
    green = _color.reflectivity_color(WGRID, np.exp(-((lam - 550.)/12.)**2))
    red = _color.reflectivity_color(WGRID, (lam > 600.).astype(float))
    blue = _color.reflectivity_color(WGRID, (lam < 500.).astype(float))
    g, r, b = green['rgb255'], red['rgb255'], blue['rgb255']
    assert g[1] > g[0] and g[1] > g[2]        # 550 nm band -> green channel wins
    assert r[0] > r[1] and r[0] > r[2]        # long-pass    -> red channel wins
    assert b[2] > b[0] and b[2] > b[1]        # short-pass   -> blue channel wins
    assert green['xy'][1] > 0.6               # close to the spectral green locus
    assert red['xy'][0] > 0.6                 # close to the spectral red locus


def test_out_of_range_spectrum_is_flagged():
    """A spectrum that stops inside the visible window is edge-clamped and the
    'full_range' flag must say so (the colour is then only indicative)."""
    w = np.linspace(0.5, 2.0, 200)            # 620 nm and redder only
    c = _color.reflectivity_color(w, np.ones_like(w))
    assert not c['full_range']
    assert _color.reflectivity_color(WGRID, np.ones_like(WGRID))['full_range']


def test_real_metals_come_out_in_the_right_order():
    """Lorentz-Drude models (Rakic et al. 1998) of four metals: Au and Cu must
    be warm (chromaticity shifted to red/yellow with respect to D65), Ag and Al
    nearly neutral (Al is in fact marginally blue), and all four bright."""
    models = {
        'Au': (9.03, [(0.760, 0.053, 0.0), (0.024, 0.241, 0.415),
                      (0.010, 0.345, 0.830), (0.071, 0.870, 2.969),
                      (0.601, 2.494, 4.304), (4.384, 2.214, 13.32)]),
        'Ag': (9.01, [(0.845, 0.048, 0.0), (0.065, 3.886, 0.816),
                      (0.124, 0.452, 4.481), (0.011, 0.065, 8.185),
                      (0.840, 0.916, 9.083), (5.646, 2.419, 20.29)]),
        'Cu': (10.83, [(0.575, 0.030, 0.0), (0.061, 0.378, 0.291),
                       (0.104, 1.056, 2.957), (0.723, 3.213, 5.300),
                       (0.638, 4.305, 11.18)]),
        'Al': (14.98, [(0.523, 0.047, 0.0), (0.227, 0.333, 0.162),
                       (0.050, 0.312, 1.544), (0.166, 1.351, 1.808),
                       (0.030, 3.382, 3.473)]),
    }
    col = {m: _color.reflectivity_color(WGRID, _reflectivity(_lorentz_drude(WGRID, wp, osc)))
           for m, (wp, osc) in models.items()}
    xw = _color.XYZ_WHITE[0]/_color.XYZ_WHITE.sum()
    for m, c in col.items():
        assert c['Y'] > 0.5, f'{m} came out too dark: Y={c["Y"]}'
    for m in ('Au', 'Cu'):                       # red channel above the blue one
        assert col[m]['rgb255'][0] > col[m]['rgb255'][2], f'{m} is not warm'
    for m in ('Ag', 'Al'):                       # neutral to within a few counts
        assert abs(int(col[m]['rgb255'][0]) - int(col[m]['rgb255'][2])) <= 8, \
            f'{m} is not neutral'
    # gold and copper are visibly yellow/red, silver and aluminium are not
    assert col['Au']['xy'][0] - xw > 0.03
    assert col['Cu']['xy'][0] - xw > 0.02
    assert abs(col['Ag']['xy'][0] - xw) < 0.02
    assert abs(col['Al']['xy'][0] - xw) < 0.02
    # gold is the most saturated of the four (largest blue deficit)
    assert col['Au']['rgb255'][0] - col['Au']['rgb255'][2] > \
        col['Ag']['rgb255'][0] - col['Ag']['rgb255'][2]


def test_gamut_clipping_is_reported():
    """XYZ outside the sRGB gamut (or brighter than white) must be clipped and
    flagged rather than silently wrapped."""
    rgb, clipped = _color.xyz_to_srgb(2.0*_color.XYZ_WHITE)
    assert clipped and np.allclose(rgb, 1.0)
    rgb, clipped = _color.xyz_to_srgb(_color.XYZ_WHITE)
    assert not clipped and np.allclose(rgb, 1.0, atol=1e-6)


# --------------------------------------------------------------------------- #
#  standalone runner (no pytest required)
# --------------------------------------------------------------------------- #


# --------------------------------------------------------------------------- #
#  Optical constants and the non-metal colour routes
# --------------------------------------------------------------------------- #
# Physics tests for the optical constants and the non-metal colour routes
# (libs/plibs/_optical.py, used by option=6 in main.py).
#
#   N = sqrt(eps) (Im N >= 0),  R = |(N-1)/(N+1)|^2,  alpha = 2 w kappa/(hbar c)
#   slab:    T = (1-R)^2 e^-ad / (1 - R^2 e^-2ad)
#   powder:  R_inf = 1 + K/S - sqrt((K/S)^2 + 2 K/S),  K/S = alpha*l
#
# Pinned by closed-form limits that leave no room for a sign or factor slip:
# a lossless slab must transmit exactly (1-R)/(1+R), the Kubelka-Munk relation
# must invert to K/S = (1-R_inf)^2/(2 R_inf), and a model dielectric with an
# absorption band in the blue must transmit YELLOW -- while the colour computed
# from its absorptance must be the complementary blue whose negative is that
# same yellow.


W = np.linspace(0.5, 6.5, 4000)


def _model_dielectric(w):
    """n ~ 1.5 glass with one weak Lorentz band at 2.9 eV (428 nm, blue)."""
    return 2.25 + 0.02/(2.9**2 - w**2 - 1j*0.25*w)


def test_refractive_index_branch():
    """sqrt(eps) must land on the passive branch: kappa >= 0 and n >= 0, also
    for a metal below the plasma edge where eps is real and negative."""
    eps = np.array([2.25+0j, 4.0+1j, -9.0+0.001j, -0.5-0.0j])
    N = _optical.refractive_index(eps)
    assert np.all(N.imag >= 0.0) and np.all(N.real >= 0.0)
    assert np.allclose(N**2, eps, atol=1e-12)          # still a square root
    assert abs(N[0] - 1.5) < 1e-12                     # transparent glass


def test_reflectivity_closed_forms():
    """R for a lossless medium is the textbook ((n-1)/(n+1))^2."""
    for n in (1.5, 2.0, 3.4):
        R = _optical.reflectivity(_optical.refractive_index(np.array([n**2+0j])))[0]
        assert abs(R - ((n - 1.)/(n + 1.))**2) < 1e-12
    assert abs(_optical.reflectivity(np.array([1.0+0j]))[0]) < 1e-15   # index-matched


def test_absorption_coefficient_definition():
    """alpha = 2*w*kappa/(hbar*c), i.e. the intensity decays as exp(-alpha z)."""
    w = np.array([2.0, 3.0])
    N = np.array([1.5+0.01j, 1.5+0.02j])
    a = _optical.absorption_coefficient(w, N)
    assert np.allclose(a, 2.*w*N.imag/(_optical.hbar*_optical.c_light))
    # equivalent optics form alpha = 4*pi*kappa/lambda
    lam = 2.*np.pi*_optical.hbar*_optical.c_light/w                    # vacuum wavelength [m]
    assert np.allclose(a, 4.*np.pi*N.imag/lam, rtol=1e-12)
    # kappa = 0.01 at 2 eV (620 nm) penetrates a few micrometres
    assert abs(1./a[0] - 620e-9/(4.*np.pi*0.01))/(1./a[0]) < 1e-2


def test_lossless_slab_transmits_the_fresnel_limit():
    """With alpha = 0 the multiple-reflection sum collapses to (1-R)/(1+R)."""
    for n in (1.5, 2.0):
        R = _optical.reflectivity(_optical.refractive_index(np.array([n**2+0j])))
        T = _optical.transmittance(R, np.zeros_like(R), 1e-3)
        assert abs(T[0] - (1. - R[0])/(1. + R[0])) < 1e-12
        assert T[0] + 2*R[0]/(1. + R[0]) - 1.0 < 1e-12     # energy accounted for
    # opaque limit and the degenerate lossless perfect mirror
    R = np.array([0.5, 1.0])
    assert _optical.transmittance(R, np.array([1e9, 0.0]), 1e-3)[0] < 1e-12
    assert _optical.transmittance(R, np.array([1e9, 0.0]), 1e-3)[1] == 0.0


def test_transmittance_is_monotonic_in_thickness():
    """Beer-Lambert: thicker sample, less light, never more."""
    eps = _model_dielectric(W)
    N = _optical.refractive_index(eps)
    R, a = _optical.reflectivity(N), _optical.absorption_coefficient(W, N)
    Ts = [_optical.transmittance(R, a, d).mean() for d in (1e-6, 1e-5, 1e-4, 1e-3)]
    assert all(t1 > t2 for t1, t2 in zip(Ts, Ts[1:]))
    assert np.all(_optical.transmittance(R, a, 1e-5) <= 1.0)


def test_kubelka_munk_inverts():
    """R_inf must satisfy the Kubelka-Munk relation K/S = (1-R)^2/(2R)."""
    ks = np.array([1e-3, 0.1, 1.0, 10., 100.])
    R = _optical.kubelka_munk(ks, 1.0)
    assert np.allclose((1. - R)**2/(2.*R), ks, rtol=1e-10)
    assert abs(_optical.kubelka_munk(np.array([0.0]), 1.0)[0] - 1.0) < 1e-15
    assert _optical.kubelka_munk(np.array([1e12]), 1.0)[0] < 1e-6
    assert np.all(np.diff(_optical.kubelka_munk(ks, 1.0)) < 0.0)      # darker with alpha


def test_blue_absorber_transmits_yellow():
    """A dielectric with an absorption band in the blue must look yellow in
    transmission, and its absorptance colour must be the complementary blue
    whose NEGATIVE is that same yellow."""
    N = _optical.refractive_index(_model_dielectric(W))
    R, a = _optical.reflectivity(N), _optical.absorption_coefficient(W, N)
    T = _optical.transmittance(R, a, 50e-6)
    A = np.clip(1. - R - T, 0., 1.)
    ct, ca = _color.spectrum_color(W, T), _color.spectrum_color(W, A)
    r, g, b = ct['rgb255']
    assert r > b and g > b, f'transmission is not yellow: {ct["hex"]}'
    assert ct['xy'][0] > 0.35                       # warm chromaticity
    # the absorptance colour is the complement: blue-dominant ...
    ra, ga, ba = ca['rgb255']
    assert ba > ra and ba > ga, f'absorptance colour is not blue: {ca["hex"]}'
    # ... and its negative is yellow again (red/green above blue)
    rn, gn, bn = ca['rgb255_neg']
    assert rn > bn and gn > bn, f'negated absorptance is not yellow: {ca["hex_neg"]}'
    assert np.all(ca['rgb255'] + ca['rgb255_neg'] == 255)


def test_powder_and_slab_agree_on_the_hue():
    """The Kubelka-Munk body colour of the same material must share the hue of
    the transmission colour (both are 'what did not get absorbed')."""
    N = _optical.refractive_index(_model_dielectric(W))
    R, a = _optical.reflectivity(N), _optical.absorption_coefficient(W, N)
    ct = _color.spectrum_color(W, _optical.transmittance(R, a, 50e-6))
    ck = _color.spectrum_color(W, _optical.kubelka_munk(a, 10e-6))
    assert abs(ct['xy'][0] - ck['xy'][0]) < 0.08
    assert abs(ct['xy'][1] - ck['xy'][1]) < 0.08
    assert ck['rgb255'][0] > ck['rgb255'][2]        # yellow, not blue


def test_metal_is_opaque_in_both_non_metal_routes():
    """A Drude metal transmits nothing at 1 mm and its powder is dark: the
    non-metal routes must not invent a colour for it."""
    eps = 1. - 0.76*9.03**2/(W*(W + 1j*0.053))
    N = _optical.refractive_index(eps)
    R, a = _optical.reflectivity(N), _optical.absorption_coefficient(W, N)
    assert _color.spectrum_color(W, _optical.transmittance(R, a, 1e-3))['Y'] < 1e-6
    assert _color.spectrum_color(W, _optical.kubelka_munk(a, 1e-6))['Y'] < 0.05
    assert _color.spectrum_color(W, R)['Y'] > 0.5          # but it is a bright mirror


# --------------------------------------------------------------------------- #
#  standalone runner (no pytest required)
# --------------------------------------------------------------------------- #


if __name__ == '__main__':
    import _tools
    sys.exit(_tools.run_standalone(globals()))
