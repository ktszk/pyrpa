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
# small one-orbital square-lattice run (Stoner factor ~0.8, below the magnetic instability)
FILL, TEMP, U, MESH, NW = 0.45, 0.02, 1.0, [8, 8, 1], 64
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
    with T.silence():
        res, node = pa.eliashberg_eilenberger_workflow.run_get_node(
            **_model_inputs(), gap_sym=orm.Int(1), wc=orm.Float(EIL['wc']), coupling=orm.Float(EIL['coupling']),
            eil_mesh=orm.List([64, 64, 1]))
    assert node.is_finished_ok
    r = res['eilenberger'].get_dict()
    assert r['form_factor'] == 'gap_wannier' and r['delta'] > 0
    # the Eilenberger step took the gap the Eliashberg step produced
    eil = next(c for c in node.called if c.process_label == 'eilenberger_bulk')
    el = next(c for c in node.called_descendants if c.process_label == 'linearized_eliashberg')
    assert eil.inputs.gap_wannier.uuid == el.outputs.gap_wannier.uuid == res['gap_wannier'].uuid


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


if __name__ == '__main__':
    sys.exit(T.run_standalone(globals()))
