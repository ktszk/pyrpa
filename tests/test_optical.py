#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Physics tests for the optical constants and the non-metal colour routes
(libs/plibs/_optical.py, used by option=6 in main.py).

  N = sqrt(eps) (Im N >= 0),  R = |(N-1)/(N+1)|^2,  alpha = 2 w kappa/(hbar c)
  slab:    T = (1-R)^2 e^-ad / (1 - R^2 e^-2ad)
  powder:  R_inf = 1 + K/S - sqrt((K/S)^2 + 2 K/S),  K/S = alpha*l

Pinned by closed-form limits that leave no room for a sign or factor slip:
a lossless slab must transmit exactly (1-R)/(1+R), the Kubelka-Munk relation
must invert to K/S = (1-R_inf)^2/(2 R_inf), and a model dielectric with an
absorption band in the blue must transmit YELLOW -- while the colour computed
from its absorptance must be the complementary blue whose negative is that
same yellow.

Runs standalone (no pytest needed):  python tests/test_optical.py
Also works under pytest if installed:  pytest tests/test_optical.py
"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from libs.plibs._optical import (refractive_index, reflectivity,
                                 absorption_coefficient, transmittance,
                                 kubelka_munk, HBAR_EVS, C_LIGHT)
from libs.plibs._color import spectrum_color

W = np.linspace(0.5, 6.5, 4000)


def _model_dielectric(w):
    """n ~ 1.5 glass with one weak Lorentz band at 2.9 eV (428 nm, blue)."""
    return 2.25 + 0.02/(2.9**2 - w**2 - 1j*0.25*w)


def test_refractive_index_branch():
    """sqrt(eps) must land on the passive branch: kappa >= 0 and n >= 0, also
    for a metal below the plasma edge where eps is real and negative."""
    eps = np.array([2.25+0j, 4.0+1j, -9.0+0.001j, -0.5-0.0j])
    N = refractive_index(eps)
    assert np.all(N.imag >= 0.0) and np.all(N.real >= 0.0)
    assert np.allclose(N**2, eps, atol=1e-12)          # still a square root
    assert abs(N[0] - 1.5) < 1e-12                     # transparent glass


def test_reflectivity_closed_forms():
    """R for a lossless medium is the textbook ((n-1)/(n+1))^2."""
    for n in (1.5, 2.0, 3.4):
        R = reflectivity(refractive_index(np.array([n**2+0j])))[0]
        assert abs(R - ((n - 1.)/(n + 1.))**2) < 1e-12
    assert abs(reflectivity(np.array([1.0+0j]))[0]) < 1e-15   # index-matched


def test_absorption_coefficient_definition():
    """alpha = 2*w*kappa/(hbar*c), i.e. the intensity decays as exp(-alpha z)."""
    w = np.array([2.0, 3.0])
    N = np.array([1.5+0.01j, 1.5+0.02j])
    a = absorption_coefficient(w, N)
    assert np.allclose(a, 2.*w*N.imag/(HBAR_EVS*C_LIGHT))
    # equivalent optics form alpha = 4*pi*kappa/lambda
    lam = 2.*np.pi*HBAR_EVS*C_LIGHT/w                    # vacuum wavelength [m]
    assert np.allclose(a, 4.*np.pi*N.imag/lam, rtol=1e-12)
    # kappa = 0.01 at 2 eV (620 nm) penetrates a few micrometres
    assert abs(1./a[0] - 620e-9/(4.*np.pi*0.01))/(1./a[0]) < 1e-2


def test_lossless_slab_transmits_the_fresnel_limit():
    """With alpha = 0 the multiple-reflection sum collapses to (1-R)/(1+R)."""
    for n in (1.5, 2.0):
        R = reflectivity(refractive_index(np.array([n**2+0j])))
        T = transmittance(R, np.zeros_like(R), 1e-3)
        assert abs(T[0] - (1. - R[0])/(1. + R[0])) < 1e-12
        assert T[0] + 2*R[0]/(1. + R[0]) - 1.0 < 1e-12     # energy accounted for
    # opaque limit and the degenerate lossless perfect mirror
    R = np.array([0.5, 1.0])
    assert transmittance(R, np.array([1e9, 0.0]), 1e-3)[0] < 1e-12
    assert transmittance(R, np.array([1e9, 0.0]), 1e-3)[1] == 0.0


def test_transmittance_is_monotonic_in_thickness():
    """Beer-Lambert: thicker sample, less light, never more."""
    eps = _model_dielectric(W)
    N = refractive_index(eps)
    R, a = reflectivity(N), absorption_coefficient(W, N)
    Ts = [transmittance(R, a, d).mean() for d in (1e-6, 1e-5, 1e-4, 1e-3)]
    assert all(t1 > t2 for t1, t2 in zip(Ts, Ts[1:]))
    assert np.all(transmittance(R, a, 1e-5) <= 1.0)


def test_kubelka_munk_inverts():
    """R_inf must satisfy the Kubelka-Munk relation K/S = (1-R)^2/(2R)."""
    ks = np.array([1e-3, 0.1, 1.0, 10., 100.])
    R = kubelka_munk(ks, 1.0)
    assert np.allclose((1. - R)**2/(2.*R), ks, rtol=1e-10)
    assert abs(kubelka_munk(np.array([0.0]), 1.0)[0] - 1.0) < 1e-15
    assert kubelka_munk(np.array([1e12]), 1.0)[0] < 1e-6
    assert np.all(np.diff(kubelka_munk(ks, 1.0)) < 0.0)      # darker with alpha


def test_blue_absorber_transmits_yellow():
    """A dielectric with an absorption band in the blue must look yellow in
    transmission, and its absorptance colour must be the complementary blue
    whose NEGATIVE is that same yellow."""
    N = refractive_index(_model_dielectric(W))
    R, a = reflectivity(N), absorption_coefficient(W, N)
    T = transmittance(R, a, 50e-6)
    A = np.clip(1. - R - T, 0., 1.)
    ct, ca = spectrum_color(W, T), spectrum_color(W, A)
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
    N = refractive_index(_model_dielectric(W))
    R, a = reflectivity(N), absorption_coefficient(W, N)
    ct = spectrum_color(W, transmittance(R, a, 50e-6))
    ck = spectrum_color(W, kubelka_munk(a, 10e-6))
    assert abs(ct['xy'][0] - ck['xy'][0]) < 0.08
    assert abs(ct['xy'][1] - ck['xy'][1]) < 0.08
    assert ck['rgb255'][0] > ck['rgb255'][2]        # yellow, not blue


def test_metal_is_opaque_in_both_non_metal_routes():
    """A Drude metal transmits nothing at 1 mm and its powder is dark: the
    non-metal routes must not invent a colour for it."""
    eps = 1. - 0.76*9.03**2/(W*(W + 1j*0.053))
    N = refractive_index(eps)
    R, a = reflectivity(N), absorption_coefficient(W, N)
    assert spectrum_color(W, transmittance(R, a, 1e-3))['Y'] < 1e-6
    assert spectrum_color(W, kubelka_munk(a, 1e-6))['Y'] < 0.05
    assert spectrum_color(W, R)['Y'] > 0.5          # but it is a bright mirror


# --------------------------------------------------------------------------- #
#  standalone runner (no pytest required)
# --------------------------------------------------------------------------- #
if __name__ == '__main__':
    import _tools
    sys.exit(_tools.run_standalone(globals()))
