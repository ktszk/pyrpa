#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Regression tests for the RPA / FLEX / Eliashberg building blocks outside the
quasiclassical Eilenberger suite.

Runs standalone (no pytest needed):  python tests/test_rpa_flex.py
Also works under pytest if installed:  pytest tests/test_rpa_flex.py
Requires the Fortran library libfmod.so (cd libs && make FC=ifx SL=MKL).
"""
import json
import os
import sys
import contextlib
import io
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


def test_chi_orb_list_is_fortran_contiguous_for_the_vertex_routines():
    """ol is declared ol(Nchi,2) in Fortran: a C-contiguous pair list would be read transposed."""
    def ref_scmat(ol, site, U, J):  # get_scmat in python, indexing ol the way we read it here
        Nchi, Up = len(ol), U - 2 * J
        Smat, Cmat = np.zeros((Nchi, Nchi)), np.zeros((Nchi, Nchi))
        for i in range(Nchi):
            for j in range(Nchi):
                if site[i] != site[j]:
                    continue
                (a, b), (c, d) = ol[i], ol[j]
                if a == b and c == d:
                    Smat[j, i], Cmat[j, i] = (U, U) if a == c else (J, 2 * Up - J)
                elif a == c and b == d:
                    Smat[j, i], Cmat[j, i] = Up, -Up + 2 * J
                elif a == d and b == c:
                    Smat[j, i], Cmat[j, i] = J, J
        return Smat, Cmat

    for Norb, site_prof in ((4, [4]), (4, [2, 2]), (6, [3, 3]), (5, [3, 2])):
        olist, site = P.get_chi_orb_list(Norb, site_prof)
        assert olist.flags['F_CONTIGUOUS'] and olist.dtype == np.int64
        Sref, Cref = ref_scmat(olist, site, 1.0, 0.2)
        Smat, Cmat = F.gen_SCmatrix(olist, site, 1.0, 0.2)
        assert np.allclose(Smat, Sref) and np.allclose(Cmat, Cref), f'site_prof={site_prof}'
        # same vertex through the orbital-resolved route (off-diagonal U entries are U'=U-2J)
        Umat = np.full((Norb, Norb), 0.6) + 0.4 * np.eye(Norb)
        Jmat = np.full((Norb, Norb), 0.2)
        So, Co = F.gen_SCmatrix_orb(olist, site, Umat, Jmat)
        assert np.allclose(So, Sref) and np.allclose(Co, Cref), f'site_prof={site_prof}'


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


def test_read_UJ_stored_sets_feed_the_orbital_dependent_vertex(tmp_path):
    """read_UJ returns Norb x Norb matrices in the get_scmat_orb convention (diag U, off-diag U')."""
    root = Path(__file__).resolve().parents[1]
    Umat, Jmat = P.read_UJ('FeSe', str(root / 'UJ'))
    assert Umat.shape == Jmat.shape == (5, 5)
    assert Umat.dtype == np.float64 and Umat.flags['F_CONTIGUOUS']
    # yz/zx (orbitals 1,2) are degenerate, and both matrices are symmetric
    assert np.allclose(Umat, Umat.T) and np.allclose(Jmat, Jmat.T)
    assert Umat[1, 1] == Umat[2, 2]
    # the same set is reachable through a combined multi-material file
    combined_all = tmp_path / 'sets.json'
    combined_all.write_text(json.dumps({'FeSe': {'U': Umat.tolist(), 'J': Jmat.tolist()},
                                        'FeTe': {'U': [[0.0]], 'J': [[0.0]]}}))
    assert np.array_equal(Umat, P.read_UJ('FeSe', str(combined_all))[0])
    # a flat Norb^2 list and nested rows describe the same matrix
    flat = tmp_path / 'Flat.json'
    flat.write_text('{"U": [1.0, 0.6, 0.6, 1.0], "J": [0.0, 0.2, 0.2, 0.0]}')
    nest = tmp_path / 'Nest.json'
    nest.write_text('{"U": [[1.0, 0.6], [0.6, 1.0]], "J": [[0.0, 0.2], [0.2, 0.0]]}')
    Uf, Jf = P.read_UJ('Flat', str(tmp_path))
    assert np.array_equal(Uf, P.read_UJ('Nest', str(tmp_path))[0])
    # with U'=U-2J the orbital-resolved vertex must reproduce the scalar Kanamori one
    olist, site = P.get_chi_orb_list(2, [2])
    assert all(np.allclose(a, b) for a, b in
               zip(F.gen_SCmatrix_orb(olist, site, Uf, Jf), F.gen_SCmatrix(olist, site, 1.0, 0.2)))
    # an unknown name is an error, not a silent fallback to some other material
    combined = tmp_path / 'all.json'
    combined.write_text('{"A": {"U": [1.0], "J": [0.0]}}')
    for where, err in ((str(combined), KeyError), (str(tmp_path / 'missing'), FileNotFoundError)):
        try:
            P.read_UJ('Nonexistent', where)
        except err:
            continue
        raise AssertionError(f'read_UJ("Nonexistent", {where!r}) should raise {err.__name__}')


def test_expand_UJ_sites_replicates_one_site_on_every_site():
    """A stored one-site U/J set fills the on-site blocks of a multi-site hamiltonian."""
    root = Path(__file__).resolve().parents[1]
    Umat, Jmat = P.read_UJ('FeSe', str(root / 'UJ'))
    # 2-Fe (10-orbital) model: the 5-orbital set goes on both sites, inter-site blocks stay 0
    Uex, Jex = P.expand_UJ_sites(Umat, Jmat, 10, [5, 5])
    assert Uex.shape == Jex.shape == (10, 10)
    assert Uex.dtype == np.float64 and Uex.flags['F_CONTIGUOUS']
    for n0 in (0, 5):
        assert np.array_equal(Uex[n0:n0 + 5, n0:n0 + 5], Umat)
        assert np.array_equal(Jex[n0:n0 + 5, n0:n0 + 5], Jmat)
    mask = P.site_block_mask(10, [5, 5])
    assert not Uex[~mask].any() and not Jex[~mask].any()
    assert np.array_equal(mask, np.array([[i // 5 == j // 5 for j in range(10)] for i in range(10)]))
    # the vertex only reads the on-site pairs, so it is the 5-orbital one repeated per site
    olist, site = P.get_chi_orb_list(10, [5, 5])
    Smat, Cmat = F.gen_SCmatrix_orb(olist, site, Uex, Jex)
    olist1, site1 = P.get_chi_orb_list(5, [5])
    Smat1, Cmat1 = F.gen_SCmatrix_orb(olist1, site1, Umat, Jmat)
    assert np.allclose(Smat[:25, :25], Smat1) and np.allclose(Smat[25:, 25:], Smat1)
    assert np.allclose(Cmat[:25, :25], Cmat1) and np.allclose(Cmat[25:, 25:], Cmat1)
    # an already global-sized matrix is passed through, a mismatching one is an error
    assert np.array_equal(P.expand_UJ_sites(Uex, Jex, 10, [5, 5])[0], Uex)
    for args in (((Umat, Jmat, 10, [1]), ), ((Umat, Jmat, 11, [6, 5]), ), ((Umat[:, :4], Jmat[:, :4], 10, [5, 5]), )):
        try:
            P.expand_UJ_sites(*args[0])
        except ValueError:
            continue
        raise AssertionError(f'expand_UJ_sites{args[0][2:]} should raise ValueError')


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


def test_regrid_sigma_bin_roundtrip_and_reseed(tmp_path):
    """plibs.regrid_sigma_bin (T-annealing helper): the scipy reader/writer is
    binary-compatible with io_sigma, an identity re-grid reproduces the stored
    sigma exactly at every kept frequency (the dropped last point carries the
    wrap-around artifact), the result is causal, and a re-gridded seed with
    changed (T, Nw) is accepted by mkself via sw_in."""
    from scipy.io import FortranFile
    m = _tiny_one_orbital_model(Nx=4, Ny=4, Nw=32, temp=0.1)
    Smat, Cmat = F.gen_SCmatrix(m['olist'], m['site'], U=1.0, J=0.0)

    def run_mkself(temp, Nw, sw_out, sw_in):
        return F.mkself(Smat, Cmat, m['kmap'], m['invk'], m['olist'], m['hamk'],
                        m['eig'], m['uni'], 0.0, 0.8, temp, Nw,
                        m['Nx'], m['Ny'], m['Nz'], sw_out, sw_in, eps=1.0e-6)

    def read_sigma_bin(Nw):
        with FortranFile('sigma.bin', 'r') as f:
            f.read_record(np.float64)
            f.read_record(np.float64)
            return f.read_record(np.complex128).reshape((-1, Nw, 1, 1), order='F')

    olddir = os.getcwd()
    os.chdir(tmp_path)
    try:
        with _silence_stdout():
            run_mkself(0.1, 32, sw_out=True, sw_in=False)
            orig = read_sigma_bin(32)
            # temp_old/Norb/Nw_old are auto-detected from the metadata footer
            Nk, Nw2 = P.regrid_sigma_bin(temp_new=0.1)
            new = read_sigma_bin(32)
        assert (Nk, Nw2) == (orig.shape[0], 32)
        # spline is exact on the kept nodes; only the dropped last point changes
        assert np.abs(new[:, :31] - orig[:, :31]).max() < 1e-10
        assert (new[:, :, 0, 0].imag <= 1e-12).all()          # causality
        # explicit values contradicting the footer must be rejected
        try:
            with _silence_stdout():
                P.regrid_sigma_bin(0.1, 0.1, Norb=1, Nw_old=16)
            raise AssertionError("wrong Nw_old must be rejected by the footer check")
        except ValueError:
            pass
        # cool + enlarge Nw (impossible with raw reuse), index-map-like w_scale
        with _silence_stdout():
            run_mkself(0.1, 32, sw_out=True, sw_in=False)
            P.regrid_sigma_bin(temp_new=0.08, Nw_new=48, w_scale=0.1 / 0.08)
            sg, mu_s = run_mkself(0.08, 48, sw_out=False, sw_in=True)
        assert np.isfinite(sg).all() and np.isfinite(mu_s)
        assert np.abs(sg).max() > 0
    finally:
        os.chdir(olddir)


def test_auto_regrid_sigma_seed(tmp_path):
    """plibs._calc._auto_regrid_sigma_seed (sw_in_self helper): a sigma.bin
    seed whose footer Nw differs from the run's is re-gridded in place before
    mkself reads it; a matching Nw and a footer-less (old-format) file are
    left byte-identical."""
    from scipy.io import FortranFile
    from libs.plibs import _calc as PC
    m = _tiny_one_orbital_model(Nx=4, Ny=4, Nw=32, temp=0.1)
    Smat, Cmat = F.gen_SCmatrix(m['olist'], m['site'], U=1.0, J=0.0)

    def run_mkself(temp, Nw, sw_out, sw_in):
        return F.mkself(Smat, Cmat, m['kmap'], m['invk'], m['olist'], m['hamk'],
                        m['eig'], m['uni'], 0.0, 0.8, temp, Nw,
                        m['Nx'], m['Ny'], m['Nz'], sw_out, sw_in, eps=1.0e-6)

    olddir = os.getcwd()
    os.chdir(tmp_path)
    try:
        with _silence_stdout():
            run_mkself(0.1, 32, sw_out=True, sw_in=False)
        meta = PC._read_sigma_bin_meta()
        assert meta is not None
        temp0, Nw0, Norb0, Nk0 = meta
        assert (round(Nw0), round(Norb0)) == (32, 1) and abs(temp0 - 0.1) < 1e-12

        # matching Nw: the seed must not be touched
        before = Path('sigma.bin').read_bytes()
        with _silence_stdout():
            PC._auto_regrid_sigma_seed(0.1, 32)
        assert Path('sigma.bin').read_bytes() == before

        # Nw mismatch (cool + enlarge): re-gridded in place, footer updated,
        # and the seed is then accepted by mkself via sw_in
        with _silence_stdout():
            PC._auto_regrid_sigma_seed(0.08, 48)
        meta = PC._read_sigma_bin_meta()
        assert (round(meta[1]), round(meta[3])) == (48, round(Nk0))
        assert abs(meta[0] - 0.08) < 1e-12
        with _silence_stdout():
            sg, mu_s = run_mkself(0.08, 48, sw_out=False, sw_in=True)
        assert np.isfinite(sg).all() and np.isfinite(mu_s)
        assert np.abs(sg).max() > 0

        # footer-less (old io_sigma) file: meta is None and the file is left alone
        with FortranFile('sigma.bin', 'r') as f:
            mu = f.read_record(np.float64)
            mu_old = f.read_record(np.float64)
            flat = f.read_record(np.complex128)
        with FortranFile('sigma.bin', 'w') as f:
            f.write_record(mu)
            f.write_record(mu_old)
            f.write_record(flat)
        assert PC._read_sigma_bin_meta() is None
        before = Path('sigma.bin').read_bytes()
        with _silence_stdout():
            PC._auto_regrid_sigma_seed(0.1, 64)
        assert Path('sigma.bin').read_bytes() == before
    finally:
        os.chdir(olddir)


# --------------------------------------------------------------------------- #
#  real-frequency chi0 kernel: static limit and near-degenerate bands
#  (occ_factor in fmod.f90, used by calc_chi / irr_chi_sc)
# --------------------------------------------------------------------------- #
def _two_band_split_model(tp):
    """Two identical square-lattice bands hybridized on site by tp, so the pair is
    split by exactly 2*tp everywhere.  tp -> 0 must be a continuous limit."""
    rvec = np.array([[i, j, 0] for i in (-1, 0, 1) for j in (-1, 0, 1)], dtype=np.float64)
    ham_r = np.zeros((len(rvec), 2, 2), dtype=np.complex128)
    for n, r in enumerate(rvec):
        if abs(r[0]) + abs(r[1]) == 1:
            ham_r[n, 0, 0] = -1.0
            ham_r[n, 1, 1] = -1.0
        if abs(r[0]) + abs(r[1]) == 0:
            ham_r[n, 0, 1] = tp
            ham_r[n, 1, 0] = tp
    return rvec, ham_r


def test_static_chi0_is_real():
    """chi0(q, w=0) is a thermodynamic (real) response: idelta must not leak into
    it.  The old kernel divided by (e_m-e_l+i*idelta) even at w=0 and produced a
    spurious imaginary part."""
    rvec, ham_r = _two_band_split_model(0.3)
    Nk, klist = P.gen_klist(12, 12, 1, sw_pp=False)
    eig, uni = P.get_eigs(klist, ham_r, np.array([]), rvec)
    ff = F.get_ffermi(eig, -0.2, 0.03)
    olist = np.array([[1, 1], [1, 2], [2, 1], [2, 2]], dtype=np.int64)
    for q in ([0., 0., 0.], [0.5, 0.5, 0.], [0.25, 0.5, 0.]):
        qs = F.get_qshift(klist, np.array(q))
        chi = F.get_chi_irr(uni, eig, ff, qs, olist, np.array([0.0]), 1.0e-3, 0.03)[0]
        assert np.abs(chi.imag).max() < 1e-12


def test_static_chi0_continuous_across_band_splitting():
    """A band pair split by less than idelta must keep its full -df/de weight.

    The old kernel treated |de| < 1e-9 as degenerate and everything above it with
    a 1/(de + i*idelta) denominator, so for 1e-9 < |de| < idelta the intra-band
    weight vanished and chi0(q=0, w=0) collapsed to exactly half its correct
    value.  Here the splitting is swept straight through that window."""
    Nk, klist = P.gen_klist(16, 16, 1, sw_pp=False)
    olist = np.array([[1, 1], [1, 2], [2, 1], [2, 2]], dtype=np.int64)
    temp, mu, idelta = 0.02, 0.0, 1.0e-3
    q, wl = np.zeros(3), np.array([0.0])
    vals = []
    for tp in (0.0, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4):
        rvec, ham_r = _two_band_split_model(tp)
        eig, uni = P.get_eigs(klist, ham_r, np.array([]), rvec)
        ff = F.get_ffermi(eig, mu, temp)
        chi = F.get_chi_irr(uni, eig, ff, F.get_qshift(klist, q), olist, wl, idelta, temp)
        vals.append(chi[0, 0, 0].real)
    vals = np.array(vals)
    assert vals[0] > 0.0
    # every splitting far below T must reproduce the degenerate value
    assert np.abs(vals - vals[0]).max() < 1e-4 * vals[0]


def test_static_chi0_matches_exact_lindhard_one_orbital():
    """Single band: chi0(q,0) must equal the delta-free static Lindhard sum."""
    rvec = np.array([[i, j, 0] for i in (-1, 0, 1) for j in (-1, 0, 1)], dtype=np.float64)
    ham_r = np.zeros((len(rvec), 1, 1), dtype=np.complex128)
    for n, r in enumerate(rvec):
        if abs(r[0]) + abs(r[1]) == 1:
            ham_r[n, 0, 0] = -1.0
    Nk, klist = P.gen_klist(24, 24, 1, sw_pp=False)
    temp, mu = 0.03, -0.4
    eig, uni = P.get_eigs(klist, ham_r, np.array([]), rvec)
    ff = F.get_ffermi(eig, mu, temp)
    olist = np.array([[1, 1]], dtype=np.int64)
    for q in ([0., 0., 0.], [0.5, 0.5, 0.], [0.5, 0.0, 0.]):     # commensurate only
        qs = F.get_qshift(klist, np.array(q))
        got = F.get_chi_irr(uni, eig, ff, qs, olist, np.array([0.0]), 1.0e-3, temp)[0, 0, 0]
        kq = qs - 1
        de = eig[:, 0] - eig[kq, 0]
        big = np.abs(de) > 1e-9
        p = np.where(big, (ff[kq, 0] - ff[:, 0]) / np.where(big, de, 1.0),
                     ff[:, 0] * (1 - ff[:, 0]) / temp)
        assert abs(got.real - p.mean()) < 1e-7
        assert abs(got.imag) < 1e-12


def test_chi0_sc_reduces_to_normal_state_at_zero_gap():
    """The BdG susceptibility must go over into the normal-state chi0 as the gap
    is switched off, continuously (the +-E Bogoliubov branches are degenerate at
    the gap nodes, exactly the window the old kernel dropped)."""
    rvec = np.array([[i, j, 0] for i in (-1, 0, 1) for j in (-1, 0, 1)], dtype=np.float64)
    ham_r = np.zeros((len(rvec), 1, 1), dtype=np.complex128)
    for n, r in enumerate(rvec):
        if abs(r[0]) + abs(r[1]) == 1:
            ham_r[n, 0, 0] = -1.0
    Nk, klist = P.gen_klist(20, 20, 1, sw_pp=False)
    temp, mu, idelta = 0.02, -0.3, 1.0e-3
    olist = np.array([[1, 1]], dtype=np.int64)
    eig, uni = P.get_eigs(klist, ham_r, np.array([]), rvec)
    hamk = F.gen_ham(klist, ham_r, rvec)
    q, wl = np.array([0.5, 0.5, 0.]), np.array([0.0])
    qs = F.get_qshift(klist, q)
    normal = F.get_chi_irr(uni, eig, F.get_ffermi(eig, mu, temp),
                           qs, olist, wl, idelta, temp)[0, 0, 0].real
    ff_row = np.cos(2 * np.pi * klist[:, 0]) - np.cos(2 * np.pi * klist[:, 1])   # d-wave
    vals = []
    for d0 in (0.0, 1e-9, 1e-7, 1e-5, 1e-3, 1e-2):
        dk = (d0 * ff_row).astype(np.complex128).reshape(Nk, 1, 1)
        eb, ub = F.get_eig(F.mkBdGhamk(hamk - mu * np.eye(1), dk))
        c = F.get_chi_irr_sc(ub, eb, F.get_ffermi(eb, 0.0, temp), qs, olist,
                             wl, idelta, temp, False)
        vals.append((c[:, :, :, 0] + c[:, :, :, 1])[0, 0, 0].real)
    vals = np.array(vals)
    # gaps far below T (the first four) must be indistinguishable from Delta=0
    assert np.abs(vals[:4] - normal).max() < 1e-6 * abs(normal)
    # and the approach stays smooth once the gap becomes comparable to T
    assert np.abs(np.diff(vals)).max() < 0.05 * abs(normal)


# --------------------------------------------------------------------------- #
#  particle-particle (Cooper) bubble: pair-degenerate line must be kept
#  (pair_factor in fmod.f90, used by calc_phi)
# --------------------------------------------------------------------------- #
def _half_filled_square(mesh=32):
    """Square lattice at half filling: xi_k = 0 on a whole nested Fermi-surface
    line that passes exactly through mesh points, so xi_m + xi_l = 0 exactly for
    many pairs - the configuration that exposes the pair-degenerate handling."""
    rvec = np.array([[i, j, 0] for i in (-1, 0, 1) for j in (-1, 0, 1)], dtype=np.float64)
    ham_r = np.zeros((len(rvec), 1, 1), dtype=np.complex128)
    for n, r in enumerate(rvec):
        if abs(r[0]) + abs(r[1]) == 1:
            ham_r[n, 0, 0] = -1.0
    Nk, klist = P.gen_klist(mesh, mesh, 1, sw_pp=False)
    eig, uni = P.get_eigs(klist, ham_r, np.array([]), rvec)
    return Nk, klist, eig, uni


def test_static_phi0_matches_exact_pair_bubble():
    """phi0(q, w=0) must equal the delta-free particle-particle sum, including
    the xi_l = -xi_m pairs where numerator and denominator vanish together and
    the ratio tends to f(1-f)/T.  The old kernel either damped those with
    i*idelta or (idelta=0, phi_qmap) skipped them outright."""
    Nk, klist, eig, uni = _half_filled_square()
    mu, olist = 0.0, np.array([[1, 1]], dtype=np.int64)
    assert (np.abs(eig[:, 0] - mu) < 1e-12).sum() > 0        # FS hits mesh points
    for temp in (0.1, 0.02):
        ff = F.get_ffermi(eig, mu, temp)
        for q in ([0., 0., 0.], [0.5, 0.5, 0.]):
            qs = F.get_iqshift(klist, np.array(q))
            got = F.get_phi_irr(uni, eig, ff, qs, olist, np.array([0.0]),
                                1.0e-3, mu, temp)[0, 0, 0]
            kq = qs - 1
            s = (eig[:, 0] - mu) + (eig[kq, 0] - mu)
            big = np.abs(s) > 1e-12
            p = np.where(big, (1.0 - ff[kq, 0] - ff[:, 0]) / np.where(big, s, 1.0),
                         ff[:, 0] * (1 - ff[:, 0]) / temp)
            assert abs(got.real - p.mean()) < 1e-7
            assert abs(got.imag) < 1e-12


def test_static_phi0_keeps_growing_as_temperature_drops():
    """The Cooper channel must not saturate: phi0(0,0) has to keep growing as
    T falls (BCS log).  Dropping the Fermi-surface pairs made it level off."""
    Nk, klist, eig, uni = _half_filled_square()
    mu, olist = 0.0, np.array([[1, 1]], dtype=np.int64)
    qs = F.get_iqshift(klist, np.zeros(3))
    temps = np.array([0.1, 0.05, 0.02, 0.01])
    phi = []
    for temp in temps:
        ff = F.get_ffermi(eig, mu, temp)
        phi.append(F.get_phi_irr(uni, eig, ff, qs, olist, np.array([0.0]),
                                 1.0e-3, mu, temp)[0, 0, 0].real)
    phi = np.array(phi)
    assert (np.diff(phi) > 0).all()                          # monotonically growing
    # growth per decade of 1/T must not shrink (a saturating phi would)
    slopes = np.diff(phi) / np.diff(np.log(1.0 / temps))
    assert (np.diff(slopes) > 0).all()


# --------------------------------------------------------------------------- #
#  round-off threshold must be relative, not absolute
# --------------------------------------------------------------------------- #
def test_chi0_survives_when_its_natural_scale_is_tiny():
    """chi0 scales as 1/bandwidth.  With every energy blown up by 1e10 the whole
    bubble sits at ~1e-10, below the ABSOLUTE 1e-9 cut the code used to apply -
    which silently returned zero.  The cut is now relative to max|chi0|."""
    st = _tiny_one_orbital_model()
    chi_ref = F.get_chi0_conv(st['Gk'], st['kmap'], st['invk'], st['olist'],
                              st['temp'], st['Nx'], st['Ny'], st['Nz'])
    scale = 1.0e10                     # G ~ 1/E, so chi ~ G^2 * T ~ 1/scale
    chi_small = F.get_chi0_conv(st['Gk'] / scale, st['kmap'], st['invk'], st['olist'],
                                st['temp'] * scale, st['Nx'], st['Ny'], st['Nz'])
    assert np.abs(chi_small).max() > 0.0                     # not wiped out
    assert np.abs(chi_small).max() < 1.0e-9                  # genuinely below the old cut
    assert np.allclose(chi_small * scale, chi_ref, rtol=1e-10, atol=0.0)


# --------------------------------------------------------------------------- #
#  gap w_n -> 0 extrapolation (fermionic Matsubara frequencies are never
#  exactly zero, so reading Delta(iw_0) directly carries an O((pi*T)^2) bias)
# --------------------------------------------------------------------------- #
def test_gap_extrapolate_w0_recovers_delta0_and_flags_iw0_bias():
    """Delta(iw_n) = Delta0 + c*w_n^2 (the leading-order shape implied by
    Delta(-iw_n)=Delta(iw_n)^*): the w_n^2 fit must recover Delta0 far more
    accurately than reading off Delta(iw_0), and the reported bias must match
    the known analytic (pi*T)^2*|c|/Delta0 offset of that iw_0 shortcut."""
    from libs.plibs._calc import gap_extrapolate_w0
    temp = 0.02
    Nw, Norb, Nk = 200, 2, 3
    wn = (2 * np.arange(Nw) + 1) * np.pi * temp
    delta0_true = np.array([0.010, 0.006])
    c_true = -0.5
    gap = np.zeros((Norb, Norb, Nw, Nk), dtype=np.complex128)
    rng = np.random.default_rng(0)
    for k in range(Nk):
        scale = 1.0 + 0.1 * k
        for i in range(Norb):
            gap[i, i, :, k] = (delta0_true[i] * scale + c_true * wn**2
                              + 1.0e-6 * rng.standard_normal(Nw))

    gap0, bias = gap_extrapolate_w0(gap, temp, n_points=4, order=1)

    for i in range(Norb):
        assert abs(gap0[i, i, 0].real - delta0_true[i]) < 1.0e-3 * delta0_true[i]
    iw0 = gap[:, :, 0, :]
    naive_err = np.abs(iw0[0, 0, 0].real - delta0_true[0]) / delta0_true[0]
    assert naive_err > 0.15                                   # the iw_0 shortcut really is biased
    # bias is |gap0 - iw0| normalized by the PEAK gap (max over all orbitals/k);
    # here that peak is delta0_true[0]*(1+0.1*(Nk-1)) at k=Nk-1, orbital 0.
    gap_scale = delta0_true[0] * (1.0 + 0.1 * (Nk - 1))
    expected_bias_00 = abs(c_true) * wn[0]**2 / gap_scale
    assert abs(bias[0, 0, 0] - expected_bias_00) < 1.0e-2 * expected_bias_00
    # off-diagonal components are identically zero and must not blow up the max
    assert np.isfinite(bias).all()
    assert bias.max() < 1.0                                   # correction stays below the peak gap


# --------------------------------------------------------------------------- #
#  irreducible k-mesh under time reversal (finv.f90 generate_irr_kpoint_inv)
# --------------------------------------------------------------------------- #
def _irr_k_violations(Nx, Ny, Nz):
    """Count exact-integer violations of the three invariants invk must satisfy.

    Every full-grid point must (a) be its irreducible representative when the
    flag is 0, (b) be MINUS it mod 1 when the flag is 1, and (c) carry the
    correct full-grid index of its own -k.  Checking on integer mesh indices
    rather than on the float coordinates is the point: the routine used to
    match coordinates with a tolerance of 1/max(Nx,Ny,Nz) -- exactly one mesh
    spacing -- and silently matched a NEIGHBOUR whenever 1/N was not a binary
    fraction.  32^3 and 64^3 were fine; 12^3 mismapped 919 of 1728 points.
    """
    A = np.array([Nx, Ny, Nz], dtype=np.int64)
    klist, kmap, invk = F.gen_irr_k_TRS(Nx, Ny, Nz)
    assert invk[:, 0].min() >= 1 and invk[:, 0].max() <= len(klist)
    ki = np.rint(klist*A).astype(np.int64) % A                 # irreducible points as indices
    j = invk[:, 0] - 1
    rep = np.where(invk[:, 1][:, None] == 0, ki[j], (-ki[j]) % A)
    negk = (-kmap) % A
    flat = negk[:, 0] + 1 + negk[:, 1]*Nx + negk[:, 2]*Nx*Ny
    return int((~(rep == kmap).all(axis=1)).sum()) + int((invk[:, 2] != flat).sum())


def test_irr_k_mapping_is_exact_on_non_binary_meshes():
    """The TRS mapping must be exact for every mesh, not only for powers of two."""
    for mesh in [(2, 2, 2), (4, 4, 4), (5, 5, 5), (6, 6, 6), (9, 9, 9), (12, 12, 12),
                 (3, 3, 1), (6, 4, 2), (7, 7, 3), (10, 6, 5), (1, 1, 5), (8, 6, 4),
                 (5, 4, 3), (3, 2, 4), (7, 6, 5), (1, 2, 1)]:
        assert _irr_k_violations(*mesh) == 0, f"mesh {mesh}"


def test_irr_k_covers_the_full_grid_exactly_once():
    """Every full-grid point is claimed, and the weights sum back to Nkall."""
    for Nx, Ny, Nz in [(6, 6, 6), (4, 5, 3), (5, 4, 3), (12, 8, 6)]:
        klist, kmap, invk = F.gen_irr_k_TRS(Nx, Ny, Nz)
        Nkall = Nx*Ny*Nz
        assert len(kmap) == Nkall
        # kmap must be a permutation of the whole [0,Nx)x[0,Ny)x[0,Nz) grid
        flat = kmap[:, 0] + Nx*(kmap[:, 1] + Ny*kmap[:, 2])
        assert np.array_equal(np.sort(flat), np.arange(Nkall))
        # -k is an involution on the full grid
        assert np.array_equal(invk[invk[:, 2] - 1, 2], np.arange(1, Nkall + 1))
        # multiplicities of the irreducible points add up to the full grid
        assert np.bincount(invk[:, 0], minlength=len(klist) + 1).sum() == Nkall


def _Nk_from_trs_note(Nx, Ny, Nz):
    """N_k of the irreducible wedge, transcribed from TRS_irr.typ.

    k_z = 0 (and k_z = pi when N_z is even) are self-paired planes and take the
    2-D reduction; the remaining planes pair up with their -k_z partner and are
    kept whole.  The 2-D count depends on the parity of N_x and N_y, which is
    what leaves the +4 / +2 / +1 offsets on N/2.
    """
    ex, ey, ez = Nx % 2 == 0, Ny % 2 == 0, Nz % 2 == 0
    N = Nx*Ny*Nz
    if ez:
        return N//2 + 4 if (ex and ey) else (N//2 + 2 if (ex or ey) else N//2 + 1)
    return N//2 + 2 if (ex and ey) else (N//2 + 1 if (ex or ey) else (N + 1)//2)


def test_irr_wedge_matches_the_trs_note_for_every_parity():
    """gen_irr_k must build the wedge TRS_irr.typ describes, for all 8 parities.

    Four properties pin it down: the point count matches the note, the points
    are distinct, no k and -k are both present (that would double-count), and
    k plus -k covers the full grid (nothing is lost).

    Nx odd with Ny even used to fail all four: that branch of gen_irr_k was a
    bare `continue`, so the wedge was left unfilled and the routine wrote past
    the end of klist.  The note only derives the mixed-parity case for Nx even /
    Ny odd; the COUNT is symmetric under swapping the axes, but the
    CONSTRUCTION is not -- with Nx even both kx=0 and kx=pi are self-paired and
    lose half their ky, with Nx odd only kx=0 is.
    """
    meshes = [(4, 4, 4), (5, 5, 5), (4, 5, 5), (5, 4, 5),      # ee/oo/eo/oe with Nz odd
              (4, 4, 5), (5, 5, 4), (4, 5, 4), (5, 4, 4),      # ... and with Nz even
              (1, 2, 1), (3, 2, 4), (7, 6, 5), (1, 8, 6),      # more Nx-odd/Ny-even
              (6, 4, 2), (9, 9, 9), (12, 8, 6), (2, 2, 2)]
    for Nx, Ny, Nz in meshes:
        A = np.array([Nx, Ny, Nz], dtype=np.int64)
        klist, _, _ = F.gen_irr_k_TRS(Nx, Ny, Nz)
        ki = np.rint(klist*A).astype(np.int64) % A
        pts = {tuple(r) for r in ki}
        neg = {tuple(r) for r in ((-ki) % A)}
        mesh = (Nx, Ny, Nz)
        assert len(ki) == _Nk_from_trs_note(*mesh), f"{mesh}: count"
        assert len(pts) == len(ki), f"{mesh}: duplicate wedge points"
        assert not any(a != b and b in pts
                       for a, b in zip((tuple(r) for r in ki),
                                       (tuple(r) for r in ((-ki) % A)))), f"{mesh}: k and -k"
        assert len(pts | neg) == Nx*Ny*Nz, f"{mesh}: does not cover the grid"


def _multiorbital_trs_model(Norb=3, Nx=4, Ny=4, Nz=1, Nw=4, temp=0.05, mu=0.13, seed=3):
    """Multi-orbital model with H(R) real (TRS) and H(-R)=H(R)^T (hermitian H(k)), but
    H(R) deliberately NOT symmetric.  A one-orbital model, or one with symmetric H(k),
    cannot distinguish a transposed orbital index in the G(-k)=G(k)^T expansion."""
    rng = np.random.default_rng(seed)
    hop = {}
    for R in [(a, b, 0) for a in (-1, 0, 1) for b in (-1, 0, 1)]:
        if R == (0, 0, 0):
            A = rng.standard_normal((Norb, Norb))
            hop[R] = A + A.T
        elif R not in hop:
            A = rng.standard_normal((Norb, Norb))
            hop[R] = A
            hop[(-R[0], -R[1], -R[2])] = A.T
    rvec = np.array([[float(a) for a in R] for R in hop], dtype=np.float64)
    ham_r = np.array([hop[tuple(int(x) for x in R)] for R in rvec], dtype=np.complex128)
    assert np.abs(ham_r - np.transpose(ham_r, (0, 2, 1))).max() > 1e-3, "H(R) must be asymmetric"
    klist, kmap, invk = F.gen_irr_k_TRS(Nx, Ny, Nz)
    Nk = len(klist)
    Gk = F.gen_green(np.zeros((Norb, Norb, Nw, Nk), dtype=np.complex128),
                     F.gen_ham(klist, ham_r, rvec), mu, temp)
    olist, site = P.get_chi_orb_list(Norb, [1])
    Smat, Cmat = F.gen_SCmatrix(olist, site, U=0.0, J=0.0)
    return dict(Norb=Norb, Nx=Nx, Ny=Ny, Nz=Nz, Nw=Nw, Nk=Nk, Nkall=Nx*Ny*Nz, temp=temp,
                mu=mu, rvec=rvec, ham_r=ham_r, klist=klist, kmap=kmap, invk=invk,
                Gk=Gk, olist=olist, site=site, Smat=Smat, Cmat=Cmat)


def test_green_symmetry_relations_used_by_chi0_and_sigma():
    """The relations the k->-k / w->-w expansions rely on, for hermitian H with TRS:
    G(k,-iw) = G(k,iw)^dagger  and  G(-k,iw) = G(k,iw)^T."""
    st = _multiorbital_trs_model()
    rvec, ham_r, temp, mu, Norb = st['rvec'], st['ham_r'], st['temp'], st['mu'], st['Norb']

    def G(kf, n):
        w = (2 * n + 1) * np.pi * temp
        H = np.einsum('r,rab->ab', np.exp(-2j * np.pi * (np.asarray(kf) @ rvec.T)), ham_r)
        return np.linalg.inv((1j * w + mu) * np.eye(Norb) - H)

    rng = np.random.default_rng(11)
    ks = rng.random((8, 3))
    ks[:, 2] = 0.0
    for k in ks:
        for n in range(3):
            g = G(k, n)
            assert np.allclose(G(k, -n - 1), g.conj().T, atol=1e-12)   # hermiticity
            assert np.allclose(G(-k, n), g.T, atol=1e-12)              # TRS
            assert np.allclose(G(-k, -n - 1), g.conj(), atol=1e-12)    # both
    # the relations without the transpose must NOT hold, else the test is blind
    g = G(ks[0], 0)
    assert not np.allclose(G(ks[0], -1), g.conj(), atol=1e-6)
    assert not np.allclose(G(-ks[0], 0), g, atol=1e-6)


def test_chi0_symmetry_expansion_matches_unreduced_k_grid():
    """get_chi0 rebuilds half the BZ through G(-k,iw)=G(k,iw)^T.  Feeding it G on the full
    grid with every k its own representative bypasses that branch entirely; both must agree."""
    st = _multiorbital_trs_model()
    Nx, Ny, Nz, Nw, Norb = st['Nx'], st['Ny'], st['Nz'], st['Nw'], st['Norb']
    Nkall, kmap, invk = st['Nkall'], st['kmap'], st['invk']
    with _silence_stdout():
        chi_a, _ = F.get_chi0(st['Smat'], st['Cmat'], st['Gk'], st['olist'],
                              kmap, invk, st['temp'], Nx, Ny, Nz)
        kfull = np.ascontiguousarray(kmap / np.array([Nx, Ny, Nz]))
        invk_full = np.zeros_like(invk)
        invk_full[:, 0] = np.arange(1, Nkall + 1)
        invk_full[:, 2] = invk[:, 2]
        Gk_full = F.gen_green(np.zeros((Norb, Norb, Nw, Nkall), dtype=np.complex128),
                              F.gen_ham(kfull, st['ham_r'], st['rvec']), st['mu'], st['temp'])
        chi_b, _ = F.get_chi0(st['Smat'], st['Cmat'], Gk_full, st['olist'],
                              kmap, invk_full, st['temp'], Nx, Ny, Nz)
    sel = np.zeros(st['Nk'], dtype=int)
    for i in range(Nkall):
        if invk[i, 1] == 0:
            sel[invk[i, 0] - 1] = i
    assert np.allclose(chi_a, chi_b[:, :, :, sel], atol=1e-12)


def test_flex_filling_weights_irreducible_k_by_multiplicity():
    """The TRS wedge is not uniformly weighted: the time-reversal invariant momenta are
    their own partner (multiplicity 1) while every other k has 2.  An unweighted average
    over the wedge biases n(mu) -- and hence mu -- by O(1/Nk)."""
    for (Nx, Ny, Nz) in [(4, 4, 1), (8, 8, 1), (4, 4, 2)]:
        klist, kmap, invk = F.gen_irr_k_TRS(Nx, Ny, Nz)
        Nk, Nkall = len(klist), Nx * Ny * Nz
        mult = np.zeros(Nk)
        for i in range(Nkall):
            mult[invk[i, 0] - 1] += 1
        assert mult.sum() == Nkall
        assert set(np.unique(mult)) <= {1.0, 2.0}
        assert (mult == 1).sum() > 0, "there must be TRIM points, else the bias is invisible"
        temp = 0.1
        eig = -2 * (np.cos(2 * np.pi * klist[:, 0]) + np.cos(2 * np.pi * klist[:, 1]))
        occ = 0.5 * (1 - np.tanh(0.5 * (eig - 0.35) / temp))
        assert not np.isclose(occ.sum() / Nk, (mult * occ).sum() / Nkall, atol=1e-6)



# --------------------------------------------------------------------------- #
#  Cartesian (R-shell) gap harmonics
# --------------------------------------------------------------------------- #
_HARM_RVEC = np.array([[i, j, k] for i in (-2, -1, 0, 1, 2) for j in (-2, -1, 0, 1, 2)
                       for k in (-2, -1, 0, 1, 2)], dtype=np.float64)
_HARM_GAPS = (0, 1, 2, 3, 4, 5, -1, -2, 6, 7, -3)


def _cart_rot(kfrac, bvec, M):
    """Rotate k in CARTESIAN space by M and return the fractional coordinates."""
    return ((kfrac @ bvec) @ M.T) @ np.linalg.inv(bvec)


def test_cartesian_harmonic_reproduces_fractional_on_orthogonal_lattices():
    """The Cartesian R-shell route must be a drop-in replacement wherever the fractional
    formulas are already right: an orthogonal single-site lattice, every gap_sym
    (chiral ones included) and every axis ratio (a=b<c, a=b>c, a!=b!=c, cubic).
    Bit-level agreement is what makes switching it on by default safe for old inputs."""
    rng = np.random.default_rng(3)
    kf = rng.random((400, 3))
    deg = np.array([90., 90., 90.])
    for abc in ([4., 4., 6.], [4., 4., 3.], [3., 5., 7.], [4., 4., 4.]):
        avec, _ = P.get_ptv(np.array(abc), deg, 0)
        for gs in _HARM_GAPS:
            old = P.gap_symms(kf, 1, gs)[0]
            new = P.gap_symms(kf, 1, gs, avec=avec, rvec=_HARM_RVEC)[0]
            assert np.abs(old - new).max() < 1e-12, (abc, gs)


def test_cartesian_harmonic_is_periodic_and_carries_the_labelled_irrep():
    """On a centred lattice the fractional harmonic is periodic but is NOT the symmetry
    it is named after: cos(2 pi k_x) - cos(2 pi k_y) fails to change sign under the
    Cartesian C4z of an fcc/bcc cell.  The R-shell harmonic does, because the symmetry
    sits in the coefficients K(R_hat) while the phase 2 pi k.n stays basis independent
    -- so no axis decomposition (sw_dec_axis) is needed and periodicity is untouched."""
    rng = np.random.default_rng(7)
    kf = rng.random((300, 3))
    deg = np.array([90., 90., 90.])
    C4 = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    for brav in (1, 2, 6, 7):
        avec, _ = P.get_ptv(np.array([4., 4., 4.]), deg, brav)
        bvec = P.get_bvec(avec)
        kw = dict(avec=avec, rvec=_HARM_RVEC)
        k4 = _cart_rot(kf, bvec, C4)
        p0 = P.gap_symms(kf, 1, 1, **kw)[0].real
        p4 = P.gap_symms(k4, 1, 1, **kw)[0].real
        assert np.abs(p0).max() > 0.1
        assert np.abs(p4 + p0).max() / np.abs(p0).max() < 1e-10, brav   # B1g: phi(C4k) = -phi(k)
        for e in np.eye(3):                                             # periodic in every G
            assert np.abs(p0 - P.gap_symms(kf + e, 1, 1, **kw)[0].real).max() < 1e-10
        f0 = P.gap_symms(kf, 1, 1)[0].real                              # fractional: not B1g
        f4 = P.gap_symms(k4, 1, 1)[0].real
        assert np.abs(f4 + f0).max() / np.abs(f0).max() > 0.5, brav


def test_cartesian_chiral_partners_share_a_shell_on_a_hexagonal_lattice():
    """d+id is a 2D irrep (E2g) on a hexagonal lattice, so |phi| must be invariant under
    C3z -- which holds only if BOTH partners are built on the same R set (hence the
    common cutoff in _cart_gap_row).  The fractional harmonics miss this by O(1)."""
    rng = np.random.default_rng(11)
    kf = rng.random((300, 3))
    avec, _ = P.get_ptv(np.array([4., 4., 6.]), np.array([90., 90., 90.]), 3)
    bvec = P.get_bvec(avec)
    th = 2 * np.pi / 3
    C3 = np.array([[np.cos(th), -np.sin(th), 0.], [np.sin(th), np.cos(th), 0.], [0., 0., 1.]])
    k3 = _cart_rot(kf, bvec, C3)
    kw = dict(avec=avec, rvec=_HARM_RVEC)
    a0 = np.abs(P.gap_symms(kf, 1, 6, **kw)[0])
    a3 = np.abs(P.gap_symms(k3, 1, 6, **kw)[0])
    assert np.abs(a3 - a0).max() / a0.max() < 1e-10
    f0 = np.abs(P.gap_symms(kf, 1, 6)[0])
    f3 = np.abs(P.gap_symms(k3, 1, 6)[0])
    assert np.abs(f3 - f0).max() / f0.max() > 0.5


def test_set_harmonic_lattice_switches_the_route_and_survives_a_tiny_r_set():
    """The module default (what main.py registers) must reach gap_symms, be clearable,
    and degrade to the fractional formula with a warning -- not raise -- when the R set
    of the model cannot carry the requested symmetry."""
    rng = np.random.default_rng(5)
    kf = rng.random((64, 3))
    avec, _ = P.get_ptv(np.array([4., 4., 4.]), np.array([90., 90., 90.]), 1)
    frac = P.gap_symms(kf, 1, 1)[0]
    try:
        P.set_harmonic_lattice(avec, _HARM_RVEC)
        assert np.abs(P.gap_symms(kf, 1, 1)[0]
                      - P.gap_symms(kf, 1, 1, avec=avec, rvec=_HARM_RVEC)[0]).max() == 0.0
        assert np.abs(P.gap_symms(kf, 1, 1)[0] - frac).max() > 0.1
        P.set_harmonic_lattice(avec, np.zeros((1, 3)))          # only R = 0
        with contextlib.redirect_stdout(io.StringIO()) as out:
            fallback = P.gap_symms(kf, 1, 1)[0]
        assert 'falling back' in out.getvalue()
        assert np.abs(fallback - frac).max() == 0.0
    finally:
        P.set_harmonic_lattice(None, None)
    assert np.abs(P.gap_symms(kf, 1, 1)[0] - frac).max() == 0.0


# --------------------------------------------------------------------------- #
#  standalone runner (no pytest required)
# --------------------------------------------------------------------------- #
if __name__ == '__main__':
    import _tools
    sys.exit(_tools.run_standalone(globals()))
