"""
AiiDA sample interface (examples/aiida/pyrpa_aiida.py): every calcfunction returns what the direct
plibs / flibs calls give (RPA/FLEX: mu, Sigma, chi, lambda_eliash and the gap; Eilenberger: the
bulk, surface and vortex solutions), and the workflows record the provenance of every step.  Runs
on a temporary profile (nothing stored); skipped without aiida-core.
"""
import contextlib
import os
import sys
import tempfile
from unittest import mock
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'examples', 'aiida'))
import libs.plibs as P  # noqa: E402
import libs.flibs as F  # noqa: E402
import _tools as T  # noqa: E402

# importing aiida creates its configuration folder: keep it out of $HOME when there is none yet
if 'AIIDA_PATH' not in os.environ and not os.path.isdir(os.path.expanduser('~/.aiida')):
    os.environ['AIIDA_PATH'] = tempfile.mkdtemp(prefix='pyrpa_aiida_')
try:
    import aiida
    from aiida import orm
    HAVE_AIIDA = True
except ImportError:
    HAVE_AIIDA = False

HOP = os.path.join(ROOT, 'inputs', 'square.hop')
_profile_loaded = False
# small one-orbital square-lattice run below the magnetic instability (Stoner factor 0.72 in RPA, 0.59 in
# FLEX); U = 1 is far beyond it here (RPA 2.4), where FLEX settles on unphysical Stoner > 1 fixed points
FILL, TEMP, U, MESH, NW = 0.45, 0.02, 0.3, [8, 8, 1], 64
# Eilenberger: coupling and cutoff that put T well below Tc
EIL = dict(coupling=1.0, temp=0.005, wc=0.5)


def _setup():
    """Temporary profile (once per process) and the sample module."""
    global _profile_loaded
    if not HAVE_AIIDA:
        T.skip('aiida-core not installed')
    if not os.path.exists(HOP):
        T.skip('inputs/square.hop not found')
    if not _profile_loaded:
        from aiida.storage.sqlite_temp import SqliteTempBackend
        aiida.load_profile(SqliteTempBackend.create_profile('pyrpa-test'), allow_switch=True)
        _profile_loaded = True
    import pyrpa_aiida
    return pyrpa_aiida


@contextlib.contextmanager
def _workdir():
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as d:
        os.chdir(d)
        try:
            yield d
        finally:
            os.chdir(cwd)


def _direct_normal_state(mu):
    """The RPA/FLEX ingredients without AiiDA (as pyrpa_aiida._normal_state)."""
    rvec, ham_r, Norb, _ = P.import_hoppings(HOP, 1)
    klist, kmap, invk = F.gen_irr_k_TRS(*MESH)
    hamk = F.gen_ham(klist, ham_r, rvec)
    eig, uni = F.get_eig(hamk)
    olist, site = P.get_chi_orb_list(Norb, [Norb])
    Smat, Cmat = F.gen_SCmatrix(olist, site, U=U, J=0.0)
    return dict(rvec=rvec, ham_r=ham_r, Norb=Norb, klist=klist, kmap=kmap, invk=invk, hamk=hamk, eig=eig,
                uni=uni, olist=olist, Smat=Smat, Cmat=Cmat, mu=mu)


def _model_inputs():
    return dict(hop_file=orm.SinglefileData(HOP), ftype=orm.Int(1), fill=orm.Float(FILL), temp=orm.Float(TEMP),
                U=orm.Float(U), J=orm.Float(0.0), mesh=orm.List(MESH), Nw=orm.Int(NW))


def test_aiida_load_hoppings_and_mu_match_plibs():
    pa = _setup()
    hop = pa.load_hoppings(orm.SinglefileData(HOP), orm.Int(1))
    rvec, ham_r, _, _ = P.import_hoppings(HOP, 1)
    assert np.array_equal(hop.get_array('rvec'), rvec)
    assert np.array_equal(hop.get_array('ham_r'), ham_r)
    mesh = [16, 16, 1]
    mu = pa.chemical_potential(hop, orm.Float(FILL), orm.Float(TEMP), orm.List(mesh))
    Nk, klist = P.gen_klist(*mesh, sw_pp=False)
    ref = P.calc_mu(P.get_eigs(klist, ham_r, [], rvec)[0], Nk, FILL, TEMP)
    assert abs(mu.value-ref) < 1e-12


def test_aiida_load_wannier90_single_file(tmp_path):
    pa = _setup()
    path = tmp_path/'any-name.dat'
    path.write_text('Wannier90 test\n1\n3\n1 2 1\n'
                    '0 0 0 1 1 0.5 0.0\n'
                    '1 0 0 1 1 -2.0 0.0\n'
                    '-1 0 0 1 1 -1.0 0.0\n')
    with T.silence():
        hop = pa.load_hoppings(orm.SinglefileData(str(path)), orm.Int(2))
    assert np.array_equal(hop.get_array('rvec'), [[0, 0, 0], [1, 0, 0], [-1, 0, 0]])
    assert np.allclose(hop.get_array('ham_r')[:, 0, 0], [0.5, -1.0, -1.0])
    path.write_text('Wannier90 test\n2\n1\n1\n'
                    '0 0 0 1 1 0 0\n0 0 0 2 1 0 -1\n'
                    '0 0 0 1 2 0 1\n0 0 0 2 2 0 0\n')
    with T.silence():
        hop = pa.load_hoppings(orm.SinglefileData(str(path)), orm.Int(2))
    expected = np.array([[0, 1j], [-1j, 0]])
    assert np.array_equal(hop.get_array('ham_r')[0], expected)
    assert np.array_equal(F.gen_ham(np.zeros((1, 3)), hop.get_array('ham_r'),
                                    hop.get_array('rvec'))[0], expected)


def test_aiida_harmonics_use_only_recorded_lattice():
    pa = _setup()
    avec, _ = P.get_ptv(np.array([4., 4., 4.]), np.array([90., 90., 90.]), 1)
    rvec = np.array([[i, j, k] for i in range(-2, 3) for j in range(-2, 3)
                     for k in range(-2, 3)], dtype=float)
    square_r, square_h, _, _ = P.import_hoppings(HOP, 1)
    ham_r = np.zeros((len(rvec), 1, 1), dtype=complex)
    for r, h in zip(square_r, square_h):
        ham_r[np.all(rvec == r, axis=1)] = h
    hops = pa._arrays(rvec=rvec, ham_r=ham_r)
    lattice_hops = pa._arrays(rvec=rvec, ham_r=ham_r, avec=avec, s_rcut=np.array(4.0))
    params = dict(mu=orm.Float(-0.2), temp=orm.Float(EIL['temp']), wc=orm.Float(EIL['wc']),
                  coupling=orm.Float(EIL['coupling']), gap_sym=orm.Int(2), mesh=orm.List([16, 16, 1]),
                  fs_width=orm.Float(0.1))
    captured = []
    build_fs = P.build_fs

    def capture(*args, **kwargs):
        result = build_fs(*args, **kwargs)
        captured.append(result[1].copy())
        return result

    with mock.patch.object(P, 'build_fs', capture):
        try:
            with T.silence():
                P.set_harmonic_lattice(None, None)
                pa.eilenberger_bulk(hops, **params)
                pa.eilenberger_bulk(lattice_hops, **params)
                P.set_harmonic_lattice(np.diag([2., 3., 5.]), rvec, s_rcut=3.0)
                pa.eilenberger_bulk(hops, **params)
                pa.eilenberger_bulk(lattice_hops, **params)
        finally:
            P.set_harmonic_lattice(None, None)
    assert np.allclose(captured[0], captured[2])
    assert np.allclose(captured[1], captured[3])
    assert not np.allclose(captured[0], captured[1])
    _, klist, eig, _, _ = P.get_emesh(16, 16, 1, ham_r, [], rvec, np.eye(3), sw_uni=True)
    with T.silence():
        ref = build_fs(eig, klist, -0.2, 2, 0.1, harmonic_avec=avec, harmonic_rvec=rvec,
                       harmonic_s_rcut=4.0, use_global_harmonics=False)[1]
    assert np.allclose(captured[1], ref)


def test_aiida_lattice_provenance_and_eliashberg_seed():
    pa = _setup()
    avec, _ = P.get_ptv(np.array([4., 4., 4.]), np.array([90., 90., 90.]), 1)
    lattice = orm.Dict({'avec': avec.tolist(), 's_rcut': 4.0})
    seen = []
    initial_gap = P.get_initial_gap

    def capture(*args, **kwargs):
        seed = initial_gap(*args, **kwargs)
        seen.append((args, kwargs, seed.copy()))
        return seed

    with mock.patch.object(P, 'get_initial_gap', capture), T.silence():
        result, node = pa.eliashberg_workflow.run_get_node(
            **_model_inputs(), gap_sym=orm.Int(1), lattice=lattice)
    load = next(c for c in node.called if c.process_label == 'load_hoppings')
    assert load.inputs.lattice.uuid == lattice.uuid
    assert np.allclose(result['hoppings'].get_array('avec'), avec)
    assert float(result['hoppings'].get_array('s_rcut')) == 4.0
    args, kwargs, seed = seen[0]
    assert kwargs['use_global'] is False
    assert np.allclose(kwargs['avec'], avec)
    assert np.allclose(seed, initial_gap(*args, **kwargs))


def test_aiida_rejects_invalid_lattice():
    pa = _setup()
    for lattice in (
        {'avec': [[1, 0], [0, 1]]},
        {'avec': [[1, 0, 0], [0, 1, 0], [0, 0, 0]]},
        {'avec': [[1, 0, 0], [0, 1, 0], [0, 0, 1]], 's_rcut': 0},
        {'s_rcut': 1.0},
        {'avec': [[1, 0, 0], [0, 1, 0], [0, 0, 1]], 'unknown': 1},
    ):
        with T.silence():
            T.assert_raises(ValueError, pa.load_hoppings, orm.SinglefileData(HOP), orm.Int(1),
                            lattice=orm.Dict(lattice))


def test_aiida_frequency_grid_isolation_and_restoration():
    pa = _setup()
    try:
        import sparse_ir  # noqa: F401
    except ImportError:
        T.skip('sparse-ir not installed')
    opts = pa._opts(None)
    ir_opts = pa._opts(orm.Dict({'freq_grid': 'ir', 'ir_wmax': 8.0, 'ir_eps': 1e-6}))
    st = {'eig': np.array([[0.0, 1.0]])}
    old = P.active_freq_grid()
    grid = P.ir_grid(TEMP, 10.0, eps=1e-6)
    hop = pa.load_hoppings(orm.SinglefileData(HOP), orm.Int(1))
    args = (hop, orm.Float(-0.2), orm.Float(TEMP), orm.Float(U), orm.Float(0.0),
            orm.List(MESH), orm.Int(12))
    try:
        P.set_freq_grid(None)
        with T.silence():
            ref = pa.rpa_susceptibility(*args)
        P.set_freq_grid(grid)
        assert P.active_freq_grid() is grid and not F.fg_is_uniform()
        with pa._frequency_grid(opts, TEMP, -0.2, st) as selected:
            assert selected is None and P.active_freq_grid() is None and F.fg_is_uniform()
            eig = np.zeros((1, 1))
            uni = np.ones((1, 1, 1), complex)
            green = F.gen_Green0(eig, uni, 0.0, TEMP, 4)[0, 0, :, 0]
            assert np.allclose(green, 1/(1j*P.uniform_fermion_frequencies(TEMP, 4)))
            try:
                with pa._frequency_grid(ir_opts, TEMP, -0.2, st) as inner:
                    assert P.active_freq_grid() is inner and not F.fg_is_uniform()
                    raise RuntimeError('test restoration')
            except RuntimeError as err:
                assert str(err) == 'test restoration'
            assert P.active_freq_grid() is None and F.fg_is_uniform()
        assert P.active_freq_grid() is grid and not F.fg_is_uniform()
        with pa._frequency_grid(ir_opts, TEMP, -0.2, st) as inner:
            assert inner is not grid and P.active_freq_grid() is inner and not F.fg_is_uniform()
        assert P.active_freq_grid() is grid and not F.fg_is_uniform()
        try:
            with pa._frequency_grid(ir_opts, TEMP, -0.2, st):
                raise RuntimeError('test IR restoration')
        except RuntimeError as err:
            assert str(err) == 'test IR restoration'
        assert P.active_freq_grid() is grid and not F.fg_is_uniform()
        with T.silence():
            sigma = pa._arrays(sigma=np.zeros((1, 1, 1, 1), complex))
        def incompatible_grid():
            with pa._frequency_grid(ir_opts, TEMP, -0.2, st, sigma):
                raise AssertionError('incompatible sigma was accepted')
        T.assert_raises(ValueError, incompatible_grid)
        assert P.active_freq_grid() is grid and not F.fg_is_uniform()
        with T.silence():
            actual = pa.rpa_susceptibility(*args)
        assert np.allclose(actual['chi'].get_array('chis'), ref['chi'].get_array('chis'), atol=1e-12)
        assert abs(actual['stoner'].value-ref['stoner'].value) < 1e-12
        assert P.active_freq_grid() is grid and not F.fg_is_uniform()
        try:
            with pa._frequency_grid(opts, TEMP, -0.2, st):
                raise RuntimeError('test outer restoration')
        except RuntimeError as err:
            assert str(err) == 'test outer restoration'
        assert P.active_freq_grid() is grid and not F.fg_is_uniform()
    finally:
        P.set_freq_grid(old)


def test_aiida_flex_workflow_matches_mkself_and_records_provenance():
    pa = _setup()
    with T.silence():
        res, node = pa.flex_workflow.run_get_node(**_model_inputs())
    assert node.is_finished_ok
    st = _direct_normal_state(res['mu'].value)
    with T.silence():
        sigma, mu_self = F.mkself(st['Smat'], st['Cmat'], st['kmap'], st['invk'], st['olist'], st['hamk'],
                                  st['eig'], st['uni'], st['mu'], FILL, TEMP, NW, *MESH, False, False)
    assert abs(res['mu_self'].value-mu_self) < 1e-12
    assert np.allclose(res['sigma'].get_array('sigma'), sigma, rtol=0, atol=1e-12)
    # provenance: the FLEX step consumed the hoppings and the mu that the earlier steps produced
    calls = {c.process_label: c for c in node.called}
    assert set(calls) == {'load_hoppings', 'chemical_potential', 'flex_self_energy'}
    flex = calls['flex_self_energy']
    assert flex.inputs.hoppings.uuid == res['hoppings'].uuid == calls['load_hoppings'].outputs.result.uuid
    assert flex.inputs.mu.uuid == calls['chemical_potential'].outputs.result.uuid


def test_aiida_rpa_and_flex_eliashberg_match_flibs():
    pa = _setup()
    for use_flex in (False, True):
        with T.silence():
            res, node = pa.eliashberg_workflow.run_get_node(**_model_inputs(), gap_sym=orm.Int(1),
                                                            use_flex=orm.Bool(use_flex))
        assert node.is_finished_ok
        st = _direct_normal_state(res['mu'].value)
        with T.silence():
            if use_flex:
                sigma, mu = F.mkself(st['Smat'], st['Cmat'], st['kmap'], st['invk'], st['olist'], st['hamk'],
                                     st['eig'], st['uni'], st['mu'], FILL, TEMP, NW, *MESH, False, False)
                Gk = F.gen_green(sigma, st['hamk'], mu, TEMP)
            else:
                mu = st['mu']
                Gk = F.gen_Green0(st['eig'], st['uni'], mu, TEMP, NW)
            chi0, stoner = F.get_chi0(st['Smat'], st['Cmat'], Gk, st['olist'], st['kmap'], st['invk'], TEMP, *MESH)
            chis, _ = F.get_chis_chic(chi0, st['Smat'], st['Cmat'])
            gap, lam = F.linearized_eliashberg(chi0, Gk, st['uni'], P.get_initial_gap(st['klist'], st['Norb'], 1),
                                               st['Smat'], st['Cmat'], st['olist'],
                                               F.get_plist(st['rvec'], st['ham_r']), st['kmap'], st['invk'],
                                               *MESH, TEMP, 1)
        assert abs(res['stoner'].value-stoner) < 1e-12
        assert np.allclose(res['chi'].get_array('chis'), chis, rtol=0, atol=1e-12)
        assert abs(res['lambda_eliash'].value-lam) < 1e-8
        # the gap is defined up to a phase: compare the normalized |Delta|
        g = res['gap'].get_array('gap')
        assert np.allclose(np.abs(g)/np.abs(g).max(), np.abs(gap)/np.abs(gap).max(), atol=1e-6)
        assert 'gap' in res['gap_wannier'].get_arraynames() and 'rvec' in res['gap_wannier'].get_arraynames()
        labels = {c.process_label for c in node.called}
        assert ('flex_self_energy' in labels) == use_flex
        if use_flex:   # the FLEX susceptibility is built on the FLEX Sigma and mu
            calls = {c.process_label: c for c in node.called}
            assert calls['rpa_susceptibility'].inputs.sigma.uuid == res['sigma'].uuid
            assert calls['linearized_eliashberg'].inputs.mu.uuid == res['mu_self'].uuid


def test_aiida_eilenberger_bulk_matches_calc_eilenberger():
    pa = _setup()
    mesh = [64, 64, 1]
    hop = pa.load_hoppings(orm.SinglefileData(HOP), orm.Int(1))
    mu = pa.chemical_potential(hop, orm.Float(FILL), orm.Float(TEMP), orm.List([16, 16, 1]))
    with T.silence():
        res = pa.eilenberger_bulk(hop, mu, orm.Float(EIL['temp']), orm.Float(EIL['wc']),
                                  orm.Float(EIL['coupling']), orm.Int(1), orm.List(mesh), fs_width=orm.Float(0.05))
    rvec, ham_r, _, _ = P.import_hoppings(HOP, 1)
    with _workdir(), T.silence():   # calc_eilenberger writes eilenberger_f.dat
        damp = P.calc_eilenberger(*mesh, EIL['wc'], ham_r, [], rvec, np.eye(3), mu.value, EIL['temp'], 1,
                                  EIL['coupling'], fs_width=0.05)
    r = res['result'].get_dict()
    assert r['delta'] > 0 and r['lambda_lin'] > 1
    assert abs(r['delta']-damp) < 1e-10*damp
    assert len(res['f_avg'].get_array('omega')) == len(res['f_avg'].get_array('f_avg'))


def test_aiida_eliashberg_gap_feeds_eilenberger():
    pa = _setup()
    lattice = orm.Dict({'avec': np.eye(3).tolist()})
    with T.silence():
        res, node = pa.eliashberg_eilenberger_workflow.run_get_node(
            **_model_inputs(), gap_sym=orm.Int(1), wc=orm.Float(EIL['wc']), coupling=orm.Float(EIL['coupling']),
            eil_mesh=orm.List([64, 64, 1]), lattice=lattice)
    assert node.is_finished_ok
    r = res['eilenberger'].get_dict()
    assert r['form_factor'] == 'gap_wannier' and r['delta'] > 0
    # the Eilenberger step took the gap the Eliashberg step produced
    eil = next(c for c in node.called if c.process_label == 'eilenberger_bulk')
    el = next(c for c in node.called_descendants if c.process_label == 'linearized_eliashberg')
    assert eil.inputs.gap_wannier.uuid == el.outputs.gap_wannier.uuid == res['gap_wannier'].uuid
    load = next(c for c in node.called_descendants if c.process_label == 'load_hoppings')
    assert load.inputs.lattice.uuid == lattice.uuid
    assert np.array_equal(eil.inputs.hoppings.get_array('avec'), np.eye(3))


def test_aiida_eilenberger_surface_and_vortex_match_plibs():
    pa = _setup()
    kw = dict(coupling=orm.Float(EIL['coupling']), temp=orm.Float(EIL['temp']), wc=orm.Float(EIL['wc']))
    beta = np.pi/4   # d-wave [110]: suppressed surface gap, zero-energy bound state
    with T.silence():
        surf = pa.eilenberger_surface(**kw, gap_sym=orm.Str('d'), beta_surf=orm.Float(beta), Nbeta=orm.Int(12),
                                      Lxi=orm.Float(4.0), nper=orm.Int(8))
        vort = pa.eilenberger_vortex(**kw, gap_sym=orm.Str('s'), Lxi=orm.Float(4.0), ngrid=orm.Int(31))
        with _workdir():
            x, ds = P.calc_surface(EIL['coupling'], EIL['temp'], EIL['wc'], gap_sym='d', beta_surf=beta,
                                   Nbeta=12, Lxi=4.0, nper=8)
            rg, dv = P.calc_vortex(EIL['coupling'], EIL['temp'], EIL['wc'], gap_sym='s', Lxi=4.0, ngrid=31)
    assert np.allclose(surf['profile'].get_array('x'), x) and np.allclose(surf['profile'].get_array('delta'), ds)
    assert abs(ds[0]) < 0.1*abs(ds[-1])                    # [110]: the d-wave gap vanishes at the surface
    assert {'surface_gap.dat', 'surface_ldos.dat'} <= set(surf['files'].list_object_names())
    assert np.allclose(vort['profile'].get_array('grid'), rg) and np.allclose(vort['profile'].get_array('delta'), dv)
    assert abs(dv[0]) < 1e-6*abs(dv[-1])                   # the vortex core
    assert {'vortex_gap.dat', 'vortex_ldos.dat'} <= set(vort['files'].list_object_names())


def _write_mlo_file(path):
    """Unnormalized two-orbital MLO model in the ftype=3 text format (i j R ReH ImH ReS ImS, H in Ry)."""
    m = T.two_orbital_square(Nx=8, Ny=8, Nw=16, temp=0.05, mu=-0.2)
    sr = np.zeros_like(m['ham_r'])
    sr[0] = [[0.3, 0.05], [0.05, 0.8]]
    sr[1] = sr[2] = [[0.02, 0.01], [0.01, -0.03]]
    with open(path, 'w') as f:
        for i in range(2):
            for j in range(2):
                for r, R in enumerate(m['rvec']):
                    h, s = m['ham_r'][r, i, j]/P.RY_TO_EV, sr[r, i, j]
                    f.write(f'{i+1} {j+1} {R[0]:.0f} {R[1]:.0f} {R[2]:.0f} '
                            f'{h.real:.10f} {h.imag:.10f} {s.real:.10f} {s.imag:.10f}\n')


def test_aiida_mlo_flex_and_eliashberg_match_flibs(tmp_path):
    """MLO basis (ftype=3): the overlap of the normalized orbitals is carried through mu, FLEX and the
    linearized Eliashberg equation exactly as the direct calls with overlap= do."""
    pa = _setup()
    path = str(tmp_path/'model.mlo')
    _write_mlo_file(path)
    fill, temp, mesh, nw, u, j = 0.8, 0.05, [8, 8, 1], 32, 1.0, 0.1
    with T.silence():
        res, node = pa.eliashberg_workflow.run_get_node(
            hop_file=orm.SinglefileData(path), ftype=orm.Int(3), fill=orm.Float(fill), temp=orm.Float(temp),
            U=orm.Float(u), J=orm.Float(j), mesh=orm.List(mesh), Nw=orm.Int(nw), gap_sym=orm.Int(1),
            use_flex=orm.Bool(True))
    assert node.is_finished_ok
    rvec, ham_r, S_r, Norb, _ = P.import_MLO_hoppings(path)
    assert np.allclose(res['hoppings'].get_array('S_r'), S_r)
    Nk, kl = P.gen_klist(*mesh, sw_pp=False)
    assert abs(res['mu'].value-P.calc_mu(P.get_eigs(kl, ham_r, S_r, rvec)[0], Nk, fill, temp)) < 1e-12
    klist, kmap, invk = F.gen_irr_k_TRS(*mesh)
    hamk, ovk = F.gen_ham(klist, ham_r, rvec, Ovl_r=S_r)
    eig, uni = F.get_eig(hamk, ovk)
    olist, site = P.get_chi_orb_list(Norb, [Norb])
    Smat, Cmat = F.gen_SCmatrix(olist, site, U=u, J=j)
    with T.silence():
        sigma, mu = F.mkself(Smat, Cmat, kmap, invk, olist, hamk, eig, uni, res['mu'].value, fill, temp, nw,
                             *mesh, False, False, overlap=ovk)
        Gk = F.gen_green(sigma, hamk, mu, temp, overlap=ovk)
        chi0, stoner = F.get_chi0(Smat, Cmat, Gk, olist, kmap, invk, temp, *mesh)
        _, lam = F.linearized_eliashberg(chi0, Gk, uni, P.get_initial_gap(klist, Norb, 1), Smat, Cmat, olist,
                                         F.get_plist(rvec, ham_r), kmap, invk, *mesh, temp, 1, overlap=ovk)
    assert abs(res['mu_self'].value-mu) < 1e-12
    assert np.allclose(res['sigma'].get_array('sigma'), sigma, rtol=0, atol=1e-12)
    assert abs(res['stoner'].value-stoner) < 1e-12
    assert abs(res['lambda_eliash'].value-lam) < 1e-8


def test_aiida_options_flex_info_runtime_and_storage():
    """options: an unconverged FLEX loop fails (unless allowed) and flex_info reports it; store_sigma=False
    runs FLEX -> Eliashberg as one step with the same lambda and no Sigma node; the IR grid of the FLEX
    step is reused downstream; a uniform Umat/Jmat equals the uniform U, J; every step records runtime."""
    pa = _setup()
    inputs = dict(_model_inputs(), gap_sym=orm.Int(1), use_flex=orm.Bool(True))
    with T.silence():
        ref, node = pa.eliashberg_workflow.run_get_node(**inputs)
        lean = pa.eliashberg_workflow(**inputs, options=orm.Dict({'store_sigma': False, 'store_gap': False}))
    assert ref['flex_info']['converged'] and ref['flex_info']['iterations'] > 0
    assert abs(lean['lambda_eliash'].value-ref['lambda_eliash'].value) < 1e-12
    assert 'sigma' not in lean and 'gap' not in lean and 'gap_wannier' in lean
    assert all('runtime' in c.outputs for c in node.called if c.process_label != 'load_hoppings'
               and c.process_label != 'chemical_potential')
    um = [[U]]   # one orbital: U on the site, no U', J
    with T.silence():
        orb = pa.eliashberg_workflow(**inputs, options=orm.Dict({'Umat': um, 'Jmat': [[0.0]]}))
    assert abs(orb['lambda_eliash'].value-ref['lambda_eliash'].value) < 1e-12
    with T.silence():
        ir = pa.eliashberg_workflow(**inputs, options=orm.Dict({'freq_grid': 'ir', 'ir_wmax': 10.0}))
    assert ir['flex_info']['converged']
    assert 'ir_wmax' in ir['sigma'].get_arraynames() and ir['sigma'].get_array('sigma').shape[2] < NW
    assert abs(ir['mu_self'].value-ref['mu_self'].value) < 5e-3     # uniform Nw=64 vs IR: cutoff error
    failed = False   # AiiDA re-raises what the calcfunction raised (wrapped or not, depending on the version)
    try:
        with T.silence():
            pa.flex_workflow(**_model_inputs(), options=orm.Dict({'scf_loop': 2}))
    except Exception:
        failed = True
    assert failed, 'an unconverged FLEX loop was accepted'
    with T.silence():
        loose = pa.flex_workflow(**_model_inputs(), options=orm.Dict({'scf_loop': 2, 'require_converged': False}))
    assert not loose['flex_info']['converged'] and loose['flex_info']['iterations'] == 2
    try:
        pa._opts(orm.Dict({'no_such_option': 1}))
        raise AssertionError('an unknown option was accepted')
    except ValueError:
        pass


def test_aiida_soc_flex_and_eliashberg_match_flibs():
    """soc=True: the spinful square_2orb_soc model through FLEX (mkself_soc), the SOC susceptibility and
    the SOC linearized Eliashberg equation, exactly as the direct calls with main.py's spin labels; with
    store_sigma=False the one-step path gives the same lambda."""
    pa = _setup()
    hop = os.path.join(ROOT, 'inputs', 'square_2orb_soc.hop')
    if not os.path.exists(hop):
        T.skip('inputs/square_2orb_soc.hop not found')
    fill, temp, mesh, nw, u, j = 2.0, 0.05, [8, 8, 1], 32, 0.5, 0.05
    inputs = dict(hop_file=orm.SinglefileData(hop), ftype=orm.Int(1), fill=orm.Float(fill), temp=orm.Float(temp),
                  U=orm.Float(u), J=orm.Float(j), mesh=orm.List(mesh), Nw=orm.Int(nw), gap_sym=orm.Int(1),
                  use_flex=orm.Bool(True))
    with T.silence():
        res, node = pa.eliashberg_workflow.run_get_node(**inputs, options=orm.Dict({'soc': True}))
        lean = pa.eliashberg_workflow(**inputs, options=orm.Dict({'soc': True, 'store_sigma': False}))
    assert node.is_finished_ok and res['flex_info']['converged']
    assert 'stoner' not in res and 'chiszz' in res['chi'].get_arraynames()
    rvec, ham_r, Norb, _ = P.import_hoppings(hop, 1)
    klist, kmap, invk = F.gen_irr_k_TRS(*mesh)
    hamk = F.gen_ham(klist, ham_r, rvec)
    eig, uni = F.get_eig(hamk)
    slist = np.ones(Norb, dtype=np.int64)
    slist[(Norb+1)//2:] = -1
    invs = np.concatenate([np.arange((Norb+1)//2, Norb), np.arange((Norb+1)//2)])+1
    olist, site = P.get_chi_orb_list(Norb, [Norb])
    Vmat = F.gen_Vmatrix(olist, slist, site, invs, u, j)
    with T.silence():
        sigma, mu = F.mkself_soc(Vmat, kmap, invk, invs, olist, slist, hamk, eig, uni, res['mu'].value, fill,
                                 temp, nw, *mesh, False, False)
        Gk = F.gen_green(sigma, hamk, mu, temp)
        chi, sg, sg2, isc = F.get_chi0_soc(Vmat, Gk, olist, slist, kmap, invk, invs, temp, *mesh)
        _, chiszz, _, _, _ = F.get_chis_chic_soc(chi, Vmat, olist, slist, invs)
        _, lam = F.linearized_eliashberg_soc(chi, Gk, uni, P.get_initial_gap(klist, Norb, 1), Vmat, sg, sg2,
                                             F.get_plist(rvec, ham_r), slist, olist, kmap, invk, invs, isc,
                                             *mesh, temp, 1)
    assert abs(res['mu_self'].value-mu) < 1e-12
    assert np.allclose(res['sigma'].get_array('sigma'), sigma, rtol=0, atol=1e-12)
    assert np.allclose(res['chi'].get_array('chiszz'), chiszz, rtol=0, atol=1e-12)
    assert abs(res['lambda_eliash'].value-lam) < 1e-8
    assert abs(lean['lambda_eliash'].value-lam) < 1e-8 and 'sigma' not in lean


if __name__ == '__main__':
    sys.exit(T.run_standalone(globals()))
