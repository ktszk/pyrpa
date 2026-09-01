#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Tests for the Fortran grid kernels used by the vortex / vortex-lattice solvers
(libs/src/fgrid.f90 via libs/flibs/_grid.py) and for the broadcast-aware Riccati
chord kernel (riccati_chords_bc in libs/src/feilenberger.f90).

These replaced the scipy RegularGridInterpolator / numpy versions that dominated
the solvers' runtime, so what has to be pinned is that they are the SAME
numbers, not merely similar ones:

  * bilinear_eval must reproduce RegularGridInterpolator(method='linear',
    bounds_error=False, fill_value=fill) -- including points that sit exactly on
    the grid boundary, which is where a rebuilt x0+(ng-1)*dx endpoint silently
    differs from the linspace axis by an ulp and flips whole rows of chord points
    to the fill value
  * bilinear_cell_eval must reproduce the periodic cell-centred interpolation
  * abrikosov_z must reproduce the lowest-Landau-level lattice sum
  * riccati_chords_bc must reproduce riccati_chords fed the same fields
    broadcast to [Ns,Nchord,Nw], for the full output and for a single row

Runs standalone (no pytest needed):  python tests/test_grid.py
Also works under pytest if installed:  pytest tests/test_grid.py
Requires the Fortran library libfmod.so (cd libs && make).
"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from scipy.interpolate import RegularGridInterpolator as RGI

import libs.flibs as F
import libs.plibs._eilenberger_vortex as V


def _rgi(field, xg, px, py, fill):
    """The scipy reference the Fortran kernel replaces."""
    pts = np.stack([px.ravel(), py.ravel()], axis=1)
    out = (RGI((xg, xg), field.real, bounds_error=False, fill_value=fill)(pts)
           + 1j*RGI((xg, xg), field.imag, bounds_error=False, fill_value=fill)(pts))
    return out.reshape(px.shape + ((field.shape[-1],) if field.ndim == 3 else ()))


def test_bilinear_eval_matches_scipy():
    """Same values as RegularGridInterpolator, for a 3D (frequency) field and a
    plain 2D one, inside the grid and outside it."""
    rng = np.random.default_rng(0)
    ng, nw = 21, 4
    xg = np.linspace(-3.1, 4.7, ng)
    F3 = rng.normal(size=(ng, ng, nw)) + 1j*rng.normal(size=(ng, ng, nw))
    px = rng.uniform(-5, 6, size=(30, 7))
    py = rng.uniform(-5, 6, size=(30, 7))
    got = F.bilinear_eval(F3, xg, px, py, 0.0)
    assert got.shape == (30, 7, nw)
    assert np.allclose(got, _rgi(F3, xg, px, py, 0.0), rtol=0, atol=1e-12)
    F2 = F3[:, :, 0]
    assert np.allclose(F.bilinear_eval(F2, xg, px, py, 0.0),
                       _rgi(F2, xg, px, py, 0.0), rtol=0, atol=1e-12)
    # a real field stays real and keeps the real fill: this is how the solvers
    # sample the gap amplitude (fill = Dbulk outside the cell)
    out = F.bilinear_eval(F2.real, xg, px, py, 0.3)
    assert out.dtype == np.float64
    ref = RGI((xg, xg), F2.real, bounds_error=False,
              fill_value=0.3)(np.stack([px.ravel(), py.ravel()], axis=1)).reshape(px.shape)
    assert np.allclose(out, ref, atol=1e-12)
    # the fill is returned verbatim (the scipy pair-of-interpolators idiom used
    # to fill re and im with the SAME number, i.e. fill*(1+1j) -- never used with
    # a complex field in the solvers, and not reproduced here)
    outside = F.bilinear_eval(F2, xg, np.array([99.0]), np.array([99.0]), 0.3)
    assert outside[0] == 0.3 + 0j


def test_bilinear_eval_on_the_boundary():
    """Points exactly on the axis endpoints must be interpolated, not filled: a
    chord grid rotated by beta=0 lands whole rows on the boundary, and treating
    them as outside changed the vortex gap by tens of percent."""
    rng = np.random.default_rng(1)
    ng = 41
    xg = np.linspace(-7.3123, 7.3123, ng)          # x0+(ng-1)*dx != xg[-1] exactly
    fld = rng.normal(size=(ng, ng)) + 1j*rng.normal(size=(ng, ng))
    px = np.concatenate([xg, np.full(ng, xg[-1]), np.full(ng, xg[0])])
    py = np.concatenate([xg[::-1], xg, xg])
    got = F.bilinear_eval(fld, xg, px, py, 0.0)
    assert np.allclose(got, _rgi(fld, xg, px, py, 0.0), atol=1e-12)
    assert np.count_nonzero(got) == got.size      # nothing was filled away
    # a point genuinely outside is still filled
    out = F.bilinear_eval(fld, xg, np.array([xg[-1] + 1e-6]), np.array([0.0]), 7.0)
    assert out[0] == 7.0


def test_bilinear_cell_eval_matches_the_numpy_version():
    """Periodic cell-centred interpolation, real and complex fields."""
    rng = np.random.default_rng(2)
    ng = 16
    g = {'Dlx': 2.3, 'Dly': 3.1, 'zeta': 0.5}
    px = rng.uniform(-9, 9, size=(23, 5))
    py = rng.uniform(-9, 9, size=(23, 5))
    for fld in (rng.normal(size=(ng, ng)),
                rng.normal(size=(ng, ng)) + 1j*rng.normal(size=(ng, ng))):
        ref = _periodic_eval_reference(fld, g, px, py)
        got = F.bilinear_cell_eval(fld, px, py, g['Dlx'], g['Dly'], g['zeta'])
        assert got.shape == px.shape
        assert np.iscomplexobj(got) == np.iscomplexobj(fld)
        assert np.allclose(got, ref, atol=1e-12)


def _periodic_eval_reference(field, g, px, py):
    """The numpy fancy-indexing version the Fortran kernel replaced."""
    Ng = field.shape[0]
    f1 = px/g['Dlx'] - g['zeta']*py/g['Dly']
    f2 = px/g['Dlx'] + (1.0 - g['zeta'])*py/g['Dly']
    u = ((f1 + 0.5) % 1.0)*Ng - 0.5
    v = ((f2 + 0.5) % 1.0)*Ng - 0.5
    i0 = np.floor(u).astype(int); j0 = np.floor(v).astype(int)
    wu = u - i0; wv = v - j0
    i0m = i0 % Ng; i1 = (i0 + 1) % Ng
    j0m = j0 % Ng; j1 = (j0 + 1) % Ng
    return ((1 - wu)*(1 - wv)*field[i0m, j0m] + wu*(1 - wv)*field[i1, j0m]
            + (1 - wu)*wv*field[i0m, j1] + wu*wv*field[i1, j1])


def test_abrikosov_z_matches_the_lattice_sum():
    """The LLL order parameter, against a direct numpy evaluation of the sum."""
    rng = np.random.default_rng(3)
    g = {'Dlx': 1.7, 'Dly': 2.9, 'zeta': 0.5}
    x = rng.uniform(-6, 6, size=(19, 4))
    y = rng.uniform(-6, 6, size=(19, 4))
    Dr = g['Dly']/g['Dlx']
    Dx0, Dy0 = -0.5*(1.0 + g['zeta']), -0.5
    xl, yl = x/g['Dlx'], y/g['Dly']
    ref = np.zeros(x.shape, dtype=np.complex128)
    for p in range(-6, 7):
        ref += (np.exp(-np.pi*Dr*(yl + Dy0 + p)**2)
                * np.exp(-2j*np.pi*(p*(Dx0 + p*g['zeta']*0.5) + (Dy0 + p)*xl)))
    ref *= np.sqrt(np.sqrt(2.0*Dr))
    ref = ref*np.exp(-1j*np.pi*xl*yl)
    got = F.abrikosov_z(x, y, g['Dlx'], g['Dly'], g['zeta'], 6)
    assert np.allclose(got, ref, atol=1e-12)
    # the vortex cores are its zeros: |Z| dips at the cell corners/centre
    assert np.abs(F.abrikosov_z(np.array([0.0]), np.array([0.0]),
                                g['Dlx'], g['Dly'], g['zeta'], 6))[0] < 1e-6


def test_riccati_chords_bc_matches_the_full_kernel():
    """Feeding om[Nw]+dom[Ns,Nchord] and dd[Ns,Nchord] must give exactly what the
    original kernel gives when the same fields are broadcast to [Ns,Nchord,Nw]."""
    rng = np.random.default_rng(4)
    Ns, Nc, Nw = 31, 11, 8
    om = 0.01 + 1j*(0.1 + 0.3*np.arange(Nw))
    dd = rng.normal(size=(Ns, Nc)) + 1j*rng.normal(size=(Ns, Nc))
    dom = 0.05j*rng.normal(size=(Ns, Nc))
    ds = 0.3 + np.abs(rng.normal(size=Nc))
    om3 = np.ascontiguousarray(np.broadcast_to(om, (Ns, Nc, Nw)) + dom[:, :, None])
    dd3 = np.ascontiguousarray(np.broadcast_to(dd[:, :, None], (Ns, Nc, Nw)))
    g0, f0 = F.riccati_chords(om3, dd3, 1.3, ds)
    g1, f1 = F.riccati_chords_bc(om, dd, 1.3, ds, dom=dom)
    assert g1.shape == g0.shape
    assert np.allclose(g1, g0, atol=1e-11) and np.allclose(f1, f0, atol=1e-11)
    # single-row output: same values, only the anchor row
    irow = Ns//2
    g2, f2 = F.riccati_chords_bc(om, dd, 1.3, ds, dom=dom, row=irow)
    assert g2.shape == (Nc, Nw)
    assert np.allclose(g2, g0[irow], atol=1e-11) and np.allclose(f2, f0[irow], atol=1e-11)
    # no shift at all
    g3, _ = F.riccati_chords_bc(om, dd, 1.3, ds)
    g4, _ = F.riccati_chords(np.ascontiguousarray(np.broadcast_to(om, (Ns, Nc, Nw))), dd3, 1.3, ds)
    assert np.allclose(g3, g4, atol=1e-11)


def test_solver_helpers_use_the_kernels():
    """_eval_field / _periodic_eval / _abrikosov_z are the solvers' entry points to
    the kernels; check the wiring end to end (shapes, dtypes, fill)."""
    rng = np.random.default_rng(5)
    ng, nw = 17, 3
    xg = np.linspace(-2.0, 2.0, ng)
    fld = rng.normal(size=(ng, ng, nw)) + 1j*rng.normal(size=(ng, ng, nw))
    px = rng.uniform(-3, 3, size=(11, 6)); py = rng.uniform(-3, 3, size=(11, 6))
    assert np.allclose(V._eval_field(fld, xg, px, py, 0.0),
                       _rgi(fld, xg, px, py, 0.0), atol=1e-12)
    g = {'Dlx': 1.3, 'Dly': 1.9, 'zeta': 0.0}
    A = rng.normal(size=(ng, ng))
    assert np.allclose(V._periodic_eval(A, g, px, py),
                       _periodic_eval_reference(A, g, px, py), atol=1e-12)
    Z = V._abrikosov_z(px, py, g)
    assert Z.shape == px.shape and np.iscomplexobj(Z)
    # azimuthal average of a field with a pure e^{i theta} winding: weight_phase
    # projects out the winding and returns the radial amplitude
    X, Y = np.meshgrid(xg, xg, indexing='ij')
    fw = (np.sqrt(X**2 + Y**2)*np.exp(1j*np.arctan2(Y, X)))[:, :, None]
    th = np.linspace(0.0, 2*np.pi, 64, endpoint=False)
    rg = np.array([0.5, 1.0, 1.5])
    got = V._azimuthal(fw, xg, rg, th, True)[:, 0]
    assert np.allclose(got.real, rg, rtol=2e-3) and np.allclose(got.imag, 0.0, atol=1e-9)


# --------------------------------------------------------------------------- #
#  standalone runner (no pytest required)
# --------------------------------------------------------------------------- #
if __name__ == '__main__':
    import _tools
    sys.exit(_tools.run_standalone(globals()))
