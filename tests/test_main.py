#!/usr/bin/env python
import contextlib
import os
import sys
from unittest import mock

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import main as app


def _setup():
    return app.Setup(rvec=np.zeros((1, 3)), ham_r=np.zeros((1, 1, 1), complex), S_r=[],
                     Norb=1, plist=np.ones(1), avec=np.eye(3), Arot=np.eye(3))


def test_eilenberger_wannier_fs_uses_selected_chemical_potential():
    flags = dict(eil_fs_kind='wannier', eil_fs_traj=None, eil_surf_dvector=False,
                 eil_vort_maxwell=False, eil_vort_current=False, eil_vort_chiral=False,
                 eil_vort_lattice_sc=False, eil_field_list=None)
    for mode in ('surface', 'vortex', 'vortex_dvector'):
        for fixed in (0.0, -0.8, None):
            s = _setup()
            wanted = 0.37 if fixed is None else fixed
            with contextlib.ExitStack() as stack:
                stack.enter_context(mock.patch.multiple(app, **flags, mu0=fixed,
                                                        sw_calc_mu=fixed is None,
                                                        eil_vort_dvector=mode == 'vortex_dvector'))
                solve_mu = stack.enter_context(mock.patch.object(app.plibs, 'get_mu', return_value=0.37))
                build = stack.enter_context(mock.patch.object(app.plibs, 'build_wannier_fs', return_value={}))
                surface = stack.enter_context(mock.patch.object(app.plibs, 'calc_surface'))
                vortex = stack.enter_context(mock.patch.object(app.plibs, 'calc_vortex'))
                dvector = stack.enter_context(mock.patch.object(app.plibs, 'calc_vortex_dvector'))
                if mode == 'surface':
                    app.run_eilenberger_surface(s)
                    assert surface.call_count == 1
                else:
                    app.run_eilenberger_vortex(s)
                    assert (dvector if mode == 'vortex_dvector' else vortex).call_count == 1
                assert s.mu == wanted
                assert solve_mu.call_count == (1 if fixed is None else 0)
                assert len(build.call_args_list) == (2 if mode == 'vortex_dvector' else 1)
                assert all(call.args[4] == wanted for call in build.call_args_list)


def test_eilenberger_model_fs_does_not_solve_chemical_potential():
    s = _setup()
    with mock.patch.object(app, 'eil_fs_kind', 'ellipse'), \
         mock.patch.object(app.plibs, 'get_mu') as get_mu:
        assert app.get_eil_gap_fs(s) == (None, 'ellipse')
        assert s.mu is None
        get_mu.assert_not_called()


def test_eilenberger_wannier_fs_lies_on_fixed_mu():
    rvec, ham_r, norb, _ = app.plibs.import_hoppings(
        os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'inputs', 'square.hop'), 1)
    s = app.Setup(rvec=rvec, ham_r=ham_r, S_r=[], Norb=norb, plist=np.ones(norb),
                  avec=np.eye(3), Arot=np.eye(3))
    build = app.plibs.build_wannier_fs

    def small_mesh(*args, **kwargs):
        return build(*args, mesh=64, **kwargs)

    with mock.patch.multiple(app, eil_fs_kind='wannier', mu0=-0.4, sw_calc_mu=False,
                             gap_sym=0, delta0=None), \
         mock.patch.object(app.plibs, 'get_mu', side_effect=AssertionError('fixed mu was ignored')), \
         mock.patch.object(app.plibs, 'build_wannier_fs', side_effect=small_mesh):
        fs, _ = app.get_eil_gap_fs(s)
    eig, _ = app.plibs.get_eigs(fs['kf'], ham_r, [], rvec)
    assert np.max(np.abs(eig[:, 0]-(-0.4))) < 0.006


if __name__ == '__main__':
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import _tools
    sys.exit(_tools.run_standalone(globals()))
