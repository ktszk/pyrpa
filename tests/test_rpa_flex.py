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
#  tail-corrected chi0 (get_chi0_tail / chi0_tail_impl)
# --------------------------------------------------------------------------- #
def _lindhard_setup(Nx=8, Ny=8, temp=0.04, mu=0.3):
    """One-orbital model + exact Matsubara Lindhard chi0 on all irreducible q."""
    rvec = np.array([[0,0,0],[1,0,0],[-1,0,0],[0,1,0],[0,-1,0]], dtype=np.float64)
    ham_r = np.array([[[0.0]],[[-1.0]],[[-1.0]],[[-0.5]],[[-0.5]]], dtype=np.complex128)
    klist, kmap, invk = F.gen_irr_k_TRS(Nx, Ny, 1)
    hamk = F.gen_ham(klist, ham_r, rvec)
    eig, uni = F.get_eig(hamk)
    olist, site = P.get_chi_orb_list(1, [1])
    Smat, Cmat = F.gen_SCmatrix(olist, site, U=1.0, J=0.0)
    ix, iy = np.meshgrid(np.arange(Nx), np.arange(Ny), indexing='ij')
    kf = np.stack([ix.ravel()/Nx, iy.ravel()/Ny, np.zeros(Nx*Ny)], axis=1)
    epsg = F.gen_ham(kf, ham_r, rvec)[:, 0, 0].real.reshape(Nx, Ny)
    beta = 1.0/temp
    x = beta*(epsg-mu)
    fg = np.where(x > 0, np.exp(-x)/(1+np.exp(-x)), 1/(1+np.exp(x)))

    def chi_exact(mlist):
        out = np.zeros((len(klist), len(mlist)), dtype=complex)
        for iq, q in enumerate(klist):
            iqx = int(round(q[0]*Nx)) % Nx
            iqy = int(round(q[1]*Ny)) % Ny
            eq = np.roll(np.roll(epsg, -iqx, axis=0), -iqy, axis=1)
            fq = np.roll(np.roll(fg, -iqx, axis=0), -iqy, axis=1)
            for jm, m in enumerate(mlist):
                de = 2j*np.pi*m*temp + epsg - eq
                deg = (m == 0) & (np.abs(epsg-eq) < 1e-9)
                safe = np.where(np.abs(de) < 1e-30, 1.0, de)
                out[iq, jm] = np.where(deg, fg*(1-fg)/temp, (fq-fg)/safe).mean()
        return out

    def G0(Nw):
        wl = (2*np.arange(Nw)+1)*np.pi*temp
        G = 1.0/((mu+1j*wl)[None, :] - eig[:, 0][:, None])
        return np.ascontiguousarray(G.T[None, None, :, :])

    return dict(Nx=Nx, Ny=Ny, temp=temp, mu=mu, klist=klist, kmap=kmap, invk=invk,
                hamk=hamk, eig=eig, uni=uni, olist=olist, Smat=Smat, Cmat=Cmat,
                chi_exact=chi_exact, G0=G0)


def test_chi0_tail_matches_lindhard_and_converges_second_order():
    """chi0_tail = conv[G]-conv[G0]+analytic tau reference: against the exact
    Matsubara Lindhard function the sharp-cutoff conv error falls ~1/Nw while the
    tail-corrected error falls ~1/Nw^2 (and is much smaller at fixed Nw)."""
    st = _lindhard_setup()
    mlist = [0, 1, 2, 3]
    ex = st['chi_exact'](mlist)
    errs = {}
    for Nw in (32, 128):
        G = st['G0'](Nw)
        with _silence_stdout():
            co, _ = F.get_chi0(st['Smat'], st['Cmat'], G, st['olist'], st['kmap'],
                               st['invk'], st['temp'], st['Nx'], st['Ny'], 1)
            cn, _ = F.get_chi0_tail(st['Smat'], st['Cmat'], G, st['eig'], st['uni'],
                                    st['olist'], st['kmap'], st['invk'], st['mu'],
                                    st['temp'], st['Nx'], st['Ny'], 1)
        eo = max(np.abs(co[0, 0, m, :]-ex[:, jm]).max() for jm, m in enumerate(mlist))
        en = max(np.abs(cn[0, 0, m, :]-ex[:, jm]).max() for jm, m in enumerate(mlist))
        errs[Nw] = (eo, en)
    # first vs second order: quadrupling Nw gains ~4x (old) vs ~16x (new)
    gain_old = errs[32][0]/errs[128][0]
    gain_new = errs[32][1]/errs[128][1]
    assert 2.5 < gain_old < 7.0
    assert gain_new > 9.0
    # at Nw=128 the corrected chi0 is clearly more accurate
    assert errs[128][1] < 0.35*errs[128][0]


def test_chi0_tail_exact_at_q0_and_reduces_to_conv_api():
    """At q=0, nu=0 the reference bubble tau product is constant in tau, so the
    tail-corrected chi0(0,0) equals sum f(1-f)/T to machine precision; the plain
    get_chi0_conv wrapper (refactored over chi0_conv_acc) stays consistent with
    get_chi0."""
    st = _lindhard_setup()
    Nw = 32
    G = st['G0'](Nw)
    with _silence_stdout():
        cn, _ = F.get_chi0_tail(st['Smat'], st['Cmat'], G, st['eig'], st['uni'],
                                st['olist'], st['kmap'], st['invk'], st['mu'],
                                st['temp'], st['Nx'], st['Ny'], 1)
        co, _ = F.get_chi0(st['Smat'], st['Cmat'], G, st['olist'], st['kmap'],
                           st['invk'], st['temp'], st['Nx'], st['Ny'], 1)
        chi_conv = F.get_chi0_conv(G, st['kmap'], st['invk'], st['olist'],
                                   st['temp'], st['Nx'], st['Ny'], 1)
    ex0 = st['chi_exact']([0])[0, 0]              # q index 0 = Gamma
    assert abs(cn[0, 0, 0, 0] - ex0) < 1e-10
    assert np.allclose(co, chi_conv, atol=1e-12)


def test_mkself_sw_tail_smoke():
    """FLEX loop with the tail-corrected chi0 branch converges and returns a
    finite self-energy close to (but corrected from) the sharp-cutoff one."""
    m = _tiny_one_orbital_model(Nx=4, Ny=4, Nw=16, temp=0.1)
    Smat, Cmat = F.gen_SCmatrix(m['olist'], m['site'], U=1.0, J=0.0)
    with _silence_stdout():
        sg1, mu1 = F.mkself(Smat, Cmat, m['kmap'], m['invk'], m['olist'], m['hamk'],
                            m['eig'], m['uni'], 0.0, 0.8, m['temp'], m['Nw'],
                            m['Nx'], m['Ny'], m['Nz'], False, False, sw_tail=False)
        sg2, mu2 = F.mkself(Smat, Cmat, m['kmap'], m['invk'], m['olist'], m['hamk'],
                            m['eig'], m['uni'], 0.0, 0.8, m['temp'], m['Nw'],
                            m['Nx'], m['Ny'], m['Nz'], False, False, sw_tail=True)
    assert np.isfinite(sg2).all()
    smax = np.abs(sg1).max()
    assert smax > 0
    assert np.abs(sg1-sg2).max() < 0.6*smax    # a correction, not a different answer


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
