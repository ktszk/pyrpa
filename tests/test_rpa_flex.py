#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Regression tests for the RPA / FLEX / Eliashberg building blocks outside the
quasiclassical Eilenberger suite.

Runs standalone (no pytest needed):  python tests/test_rpa_flex.py
Also works under pytest if installed:  pytest tests/test_rpa_flex.py
Requires the Fortran library libfmod.so (cd libs && make FC=ifx SL=MKL).
"""
import os
import sys
import contextlib
import io
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import libs.flibs as F
import libs.plibs as P


# --------------------------------------------------------------------------- #
#  helpers
# --------------------------------------------------------------------------- #
def _tiny_one_orbital_model(Nx=2, Ny=2, Nz=1, Nw=4, temp=0.1):
    """Small one-orbital tight-binding model for RPA/FLEX smoke tests."""
    rvec = np.array(
        [[0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]],
        dtype=np.float64,
    )
    ham_r = np.array(
        [[[0.0]], [[-1.0]], [[-1.0]], [[-0.5]], [[-0.5]]],
        dtype=np.complex128,
    )
    klist, kmap, invk = F.gen_irr_k_TRS(Nx, Ny, Nz)
    hamk = F.gen_ham(klist, ham_r, rvec)
    eig, uni = F.get_eig(hamk)
    Gk = F.gen_Green0(eig, uni, mu=0.0, temp=temp, Nw=Nw)
    olist, site = P.get_chi_orb_list(1, [1])
    Smat, Cmat = F.gen_SCmatrix(olist, site, U=0.0, J=0.0)
    plist = np.array([1.0], dtype=np.float64)
    return dict(
        Nx=Nx, Ny=Ny, Nz=Nz, Nw=Nw, temp=temp, rvec=rvec, ham_r=ham_r,
        klist=klist, kmap=kmap, invk=invk, hamk=hamk, eig=eig, uni=uni,
        Gk=Gk, olist=olist, site=site, Smat=Smat, Cmat=Cmat, plist=plist,
    )


def _silence_stdout():
    return contextlib.redirect_stdout(io.StringIO())


# --------------------------------------------------------------------------- #
#  interaction vertex / orbital-pair basis
# --------------------------------------------------------------------------- #
def test_chi_orb_list_multisite_site_length():
    """Multi-site chi orbital-pair lists must have one site label per pair."""
    chiolist, site = P.get_chi_orb_list(4, [2, 2])
    assert chiolist.shape == (8, 2)
    assert site.shape == (8,)
    assert np.array_equal(site[:4], np.ones(4, dtype=np.int64))
    assert np.array_equal(site[4:], np.full(4, 3, dtype=np.int64))
    assert chiolist.min() == 1 and chiolist.max() == 4


def test_uniform_sc_matrix_two_orbital_reference_values():
    """For a two-orbital single-site model, lock the U/J spin/charge vertices."""
    olist, site = P.get_chi_orb_list(2, [2])
    Smat, Cmat = F.gen_SCmatrix(olist, site, U=1.0, J=0.2)
    expected_s = np.array(
        [[1.0, 0.0, 0.0, 0.2],
         [0.0, 0.6, 0.2, 0.0],
         [0.0, 0.2, 0.6, 0.0],
         [0.2, 0.0, 0.0, 1.0]],
        dtype=np.float64,
    )
    expected_c = np.array(
        [[1.0, 0.0, 0.0, 1.0],
         [0.0, -0.2, 0.2, 0.0],
         [0.0, 0.2, -0.2, 0.0],
         [1.0, 0.0, 0.0, 1.0]],
        dtype=np.float64,
    )
    assert np.allclose(Smat, expected_s)
    assert np.allclose(Cmat, expected_c)


def test_orbital_dependent_vertices_reference_values():
    """Lock the orbital-dependent U/J vertex convention used by the Fortran wrapper."""
    olist, site = P.get_chi_orb_list(2, [2])
    Umat = np.full((2, 2), 1.0, dtype=np.float64)
    Jmat = np.full((2, 2), 0.2, dtype=np.float64)
    Smat_orb, Cmat_orb = F.gen_SCmatrix_orb(olist, site, Umat, Jmat)
    expected_s = np.array(
        [[1.0, 0.0, 0.0, 0.2],
         [0.0, 1.0, 0.2, 0.0],
         [0.0, 0.2, 1.0, 0.0],
         [0.2, 0.0, 0.0, 1.0]],
        dtype=np.float64,
    )
    expected_c = np.array(
        [[1.0, 0.0, 0.0, 1.8],
         [0.0, -0.6, 0.2, 0.0],
         [0.0, 0.2, -0.6, 0.0],
         [1.8, 0.0, 0.0, 1.0]],
        dtype=np.float64,
    )
    assert np.allclose(Smat_orb, expected_s)
    assert np.allclose(Cmat_orb, expected_c)


# --------------------------------------------------------------------------- #
#  RPA susceptibility algebra
# --------------------------------------------------------------------------- #
def test_scalar_rpa_chis_matches_closed_form():
    """For one orbital, chi_s = chi0 / (1 - U chi0)."""
    chi0 = np.array([[[0.20 + 0.10j]], [[0.05 + 0.02j]]], dtype=np.complex128)
    U = 1.2
    Smat = np.array([[U]], dtype=np.float64)
    chis = F.get_chis(chi0, Smat)
    expected = chi0[:, 0, 0] / (1.0 - U * chi0[:, 0, 0])
    assert chis.shape == chi0.shape
    assert np.allclose(chis[:, 0, 0], expected)


def test_chis_chic_reduce_to_bare_chi_when_vertices_are_zero():
    """With S=C=0, RPA spin/charge susceptibilities equal the bare chi0."""
    chi = np.zeros((2, 2, 1, 3), dtype=np.complex128)
    chi[:, :, 0, 0] = np.array([[0.10, 0.02j], [-0.02j, 0.20]])
    chi[:, :, 0, 1] = np.array([[0.05, 0.01], [0.03, 0.04]])
    Smat = np.zeros((2, 2), dtype=np.float64)
    Cmat = np.zeros((2, 2), dtype=np.float64)
    chis, chic = F.get_chis_chic(chi, Smat, Cmat)
    assert chis.shape == (2, 2, 3)
    assert chic.shape == (2, 2, 3)
    assert np.allclose(chis, chi[:, :, 0, :])
    assert np.allclose(chic, chi[:, :, 0, :])


# --------------------------------------------------------------------------- #
#  Green function / chi0 / FLEX vertex kernels
# --------------------------------------------------------------------------- #
def test_green0_noninteracting_inverse_for_one_orbital():
    """G0(k,iw) should be 1 / (iw_n + mu - eps_k) for a one-orbital model."""
    st = _tiny_one_orbital_model(Nw=5, temp=0.07)
    eig = st['eig'][:, 0]
    G = st['Gk'][0, 0]
    wn = (2 * np.arange(st['Nw']) + 1) * np.pi * st['temp']
    expected = 1.0 / (1j * wn[:, None] - eig[None, :])
    assert np.allclose(G, expected)


def test_chi0_convolution_matches_get_chi0_without_interaction():
    """The FLEX chi0 path and the standalone convolution helper should agree at U=J=0."""
    st = _tiny_one_orbital_model()
    with _silence_stdout():
        chi, stoner = F.get_chi0(
            st['Smat'], st['Cmat'], st['Gk'], st['olist'], st['kmap'], st['invk'],
            st['temp'], st['Nx'], st['Ny'], st['Nz'],
        )
    chi_conv = F.get_chi0_conv(
        st['Gk'], st['kmap'], st['invk'], st['olist'], st['temp'], st['Nx'], st['Ny'], st['Nz']
    )
    assert chi.shape == (1, 1, st['Nw'], len(st['klist']))
    assert np.isfinite(stoner)
    assert abs(stoner) < 1e-12
    assert np.allclose(chi, chi_conv)


def test_flex_vertex_zero_interaction_returns_zero():
    """FLEX self-energy vertex V_sigma must vanish when S=C=0."""
    chi = np.ones((3, 2, 1, 1), dtype=np.complex128) * (0.2 + 0.1j)
    Smat = np.zeros((1, 1), dtype=np.float64)
    Cmat = np.zeros((1, 1), dtype=np.float64)
    vsigma = F.get_Vsigma_nosoc_flex(chi.copy(), Smat, Cmat)
    assert vsigma.shape == chi.shape
    assert np.allclose(vsigma, 0.0)


# --------------------------------------------------------------------------- #
#  linearized / nonlinear Eliashberg smoke tests
# --------------------------------------------------------------------------- #
def test_linearized_eliashberg_zero_interaction_has_zero_lambda_and_finite_gap():
    """At S=C=0, the linearized Eliashberg eigenvalue should be zero and finite arrays returned."""
    st = _tiny_one_orbital_model()
    with _silence_stdout():
        chi, _ = F.get_chi0(
            st['Smat'], st['Cmat'], st['Gk'], st['olist'], st['kmap'], st['invk'],
            st['temp'], st['Nx'], st['Ny'], st['Nz'],
        )
        init_delta = P.get_initial_gap(st['klist'], 1, 0)
        gap, lam = F.linearized_eliashberg(
            chi, st['Gk'], st['uni'], init_delta, st['Smat'], st['Cmat'],
            st['olist'], st['plist'], st['kmap'], st['invk'],
            st['Nx'], st['Ny'], st['Nz'], st['temp'], gap_sym=0,
            itemax=3, arnoldi_m=2,
        )
    assert gap.shape == (1, 1, st['Nw'], len(st['klist']))
    assert np.isfinite(lam)
    assert abs(lam) < 1e-12
    assert np.isfinite(gap).all()


def test_nonlinear_eliashberg_zero_seed_and_zero_interaction_stays_zero():
    """The nonlinear FLEX-Eliashberg kernel should preserve the trivial zero-gap solution at S=C=0."""
    st = _tiny_one_orbital_model(Nw=3)
    delta0 = np.zeros((1, 1, st['Nw'], len(st['klist'])), dtype=np.complex128)
    with _silence_stdout():
        delta, sigmak = F.nonlinear_eliashberg(
            delta0, st['Gk'], st['hamk'], st['Smat'], st['Cmat'], st['olist'], st['plist'],
            st['kmap'], st['invk'], mu=0.0, temp=st['temp'], gap_sym=0,
            Nx=st['Nx'], Ny=st['Ny'], Nz=st['Nz'],
            sw_sigma=False, sw_Vconst=True, eps=1e-6, itemax=2, m_diis=1,
            gap_min=0.0, sw_amp_newton=False,
        )
    assert delta.shape == delta0.shape
    assert sigmak.shape == delta0.shape
    assert np.allclose(delta, 0.0)
    assert np.allclose(sigmak, 0.0)


# --------------------------------------------------------------------------- #
#  high-level calc helpers
# --------------------------------------------------------------------------- #
def test_load_sigma_from_file_missing_returns_none_in_temp_cwd():
    """_load_sigma_from_file should fail softly instead of crashing when self_en.npz is absent."""
    from libs.plibs._calc import _load_sigma_from_file

    old = os.getcwd()
    with tempfile.TemporaryDirectory() as td:
        os.chdir(td)
        try:
            with _silence_stdout():
                loaded = _load_sigma_from_file()
        finally:
            os.chdir(old)
    assert loaded is None


def test_output_gap_function_writes_expected_one_orbital_file(tmp_path):
    """The gap-output helper should write the one-orbital gap in kmap order."""
    st = _tiny_one_orbital_model(Nw=2)
    gap = np.ones((1, 1, st['Nw'], len(st['klist'])), dtype=np.complex128) * (0.3 + 0.1j)
    old = os.getcwd()
    os.chdir(tmp_path)
    try:
        with _silence_stdout():
            ret = P.output_gap_function(
                st['invk'], st['kmap'], gap, st['uni'], st['plist'], gap_sym=0, Nx=st['Nx'], sw_orb=False
            )
        out = Path('gap_11.dat')
        assert ret == 0
        assert out.exists()
        text = out.read_text()
    finally:
        os.chdir(old)
    assert '3.00000000e-01' in text
    assert '1.00000000e-01' in text


# --------------------------------------------------------------------------- #
#  standalone runner (no pytest required)
# --------------------------------------------------------------------------- #
if __name__ == '__main__':
    import inspect
    import time

    tests = [v for k, v in sorted(globals().items()) if k.startswith('test_') and callable(v)]
    npass = 0
    for t in tests:
        t0 = time.time()
        try:
            sig = inspect.signature(t)
            if 'tmp_path' in sig.parameters:
                with tempfile.TemporaryDirectory() as td:
                    t(Path(td))
            else:
                t()
            print(f"  PASS  {t.__name__:62s} ({time.time()-t0:5.1f}s)", flush=True)
            npass += 1
        except Exception as e:
            print(f"  FAIL  {t.__name__:62s} ({time.time()-t0:5.1f}s)  -> {type(e).__name__}: {e}", flush=True)
    print(f"\n{npass}/{len(tests)} passed")
    sys.exit(0 if npass == len(tests) else 1)
