#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Physics / reference tests for the reflectivity -> colour conversion
(libs/plibs/_color.py).

The chain is  R(w) -> R(lambda) -> CIE XYZ (2-deg observer, illuminant D65)
-> sRGB, the standard route for first-principles metal colours (Prandini et
al., npj Comput. Mater. 5, 129 (2019)). It is pinned here by quantities whose
values are fixed by the CIE definitions, so a corrupted colour-matching table
or a wrong normalisation shows up immediately:

  * a perfect mirror (R == 1) must land exactly on the D65 white point
    x = 0.31272, y = 0.32903, Y = 1, i.e. sRGB #ffffff
  * a grey reflector R = c keeps the white chromaticity and gives
    Y = c and the sRGB transfer function value of c
  * spectrally selective reflectors must produce the right hue (a 550 nm band
    -> green near the spectral locus, a >600 nm long-pass -> red, ...)
  * Lorentz-Drude models of real metals (Rakic et al., Appl. Opt. 37, 5271
    (1998)) must come out in the right order: Au/Cu warm (x > x_D65), Ag/Al
    nearly neutral, and every metal bright (Y > 0.5)

Runs standalone (no pytest needed):  python tests/test_color.py
Also works under pytest if installed:  pytest tests/test_color.py
"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from libs.plibs._color import (CIE_LAM, CIE_CMF, CIE_D65, EV_NM, XYZ_WHITE,
                               srgb_gamma, reflectivity_color, xyz_to_srgb)

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
    assert CIE_LAM[0] == 380. and CIE_LAM[-1] == 780. and len(CIE_LAM) == 81
    assert np.allclose(np.diff(CIE_LAM), 5.)
    assert CIE_CMF.shape == (81, 3) and CIE_D65.shape == (81,)
    assert np.allclose(XYZ_WHITE, [0.95047, 1.0, 1.08883], atol=5e-4)
    x, y = XYZ_WHITE[:2]/XYZ_WHITE.sum()
    assert abs(x - 0.31272) < 1e-4 and abs(y - 0.32903) < 1e-4


def test_perfect_mirror_is_white():
    """R == 1 defines the white point: pure white, unit luminance, no clipping."""
    c = reflectivity_color(WGRID, np.ones_like(WGRID))
    assert c['hex'] == '#ffffff'
    assert abs(c['Y'] - 1.0) < 1e-9
    assert abs(c['xy'][0] - 0.31272) < 1e-4 and abs(c['xy'][1] - 0.32903) < 1e-4
    assert c['full_range'] and not c['clipped']


def test_grey_reflector_keeps_chromaticity():
    """A flat R = c is neutral: white chromaticity, Y = c, sRGB = gamma(c)."""
    for c0 in (0.2, 0.5, 0.85):
        c = reflectivity_color(WGRID, c0*np.ones_like(WGRID))
        assert abs(c['Y'] - c0) < 1e-9
        assert np.allclose(c['xy'], XYZ_WHITE[:2]/XYZ_WHITE.sum(), atol=1e-9)
        assert np.allclose(c['rgb'], srgb_gamma(np.full(3, c0)), atol=1e-9)


def test_selective_reflectors_have_the_right_hue():
    """Narrow-band and step reflectors must map onto the expected hues."""
    lam = EV_NM/WGRID
    green = reflectivity_color(WGRID, np.exp(-((lam - 550.)/12.)**2))
    red = reflectivity_color(WGRID, (lam > 600.).astype(float))
    blue = reflectivity_color(WGRID, (lam < 500.).astype(float))
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
    c = reflectivity_color(w, np.ones_like(w))
    assert not c['full_range']
    assert reflectivity_color(WGRID, np.ones_like(WGRID))['full_range']


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
    col = {m: reflectivity_color(WGRID, _reflectivity(_lorentz_drude(WGRID, wp, osc)))
           for m, (wp, osc) in models.items()}
    xw = XYZ_WHITE[0]/XYZ_WHITE.sum()
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
    rgb, clipped = xyz_to_srgb(2.0*XYZ_WHITE)
    assert clipped and np.allclose(rgb, 1.0)
    rgb, clipped = xyz_to_srgb(XYZ_WHITE)
    assert not clipped and np.allclose(rgb, 1.0, atol=1e-6)


# --------------------------------------------------------------------------- #
#  standalone runner (no pytest required)
# --------------------------------------------------------------------------- #
if __name__ == '__main__':
    import _tools
    sys.exit(_tools.run_standalone(globals()))
