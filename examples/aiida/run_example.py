#!/usr/bin/env python
"""
Run the pyrpa AiiDA sample (pyrpa_aiida.py) on inputs/square.hop.

  python examples/aiida/run_example.py [--calc CALC] [--profile NAME]

  CALC = flex         hopping file -> mu -> FLEX self-energy
         eliashberg   ... -> chi_s/chi_c, Stoner factor and the linearized Eliashberg lambda, gap
                      (RPA; --flex: with the FLEX self-energy)
         eilenberger  eliashberg, then the homogeneous Eilenberger gap with that gap as the form factor
         surface      Eilenberger surface state (model FS, d-wave [110]: zero-energy bound state)
         vortex       Eilenberger isolated vortex (model FS, s-wave)

The default temporary profile (SqliteTempBackend) needs no database, broker or `verdi presto` and
keeps nothing; --profile NAME stores the nodes in an existing profile.  The AiiDA configuration
folder is still created under $AIIDA_PATH (default ~/.aiida).
"""
import argparse
import os
import sys

import aiida
from aiida import orm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pyrpa_aiida as pa  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))


def model_inputs(hop: str) -> dict:
    """Inputs shared by the RPA/FLEX workflows: one-orbital square lattice near half filling."""
    return dict(hop_file=orm.SinglefileData(hop), ftype=orm.Int(1), fill=orm.Float(0.45), temp=orm.Float(0.02),
                U=orm.Float(1.0), J=orm.Float(0.0), mesh=orm.List([16, 16, 1]), Nw=orm.Int(128))


def print_provenance(node) -> None:
    print('provenance (processes called by the workflow):')
    for sub in node.called_descendants:
        ins = ', '.join(sorted(sub.inputs))
        outs = ', '.join(sorted(sub.outputs))
        print(f'  {sub.process_label:24s} pk={sub.pk:<4d} in: {ins}  out: {outs}')


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--calc', default='flex', choices=['flex', 'eliashberg', 'eilenberger', 'surface', 'vortex'])
    ap.add_argument('--flex', action='store_true', help='eliashberg/eilenberger: FLEX instead of RPA')
    ap.add_argument('--profile', default=None, help='existing AiiDA profile (default: a temporary one)')
    ap.add_argument('--hop', default=os.path.join(ROOT, 'inputs', 'square.hop'), help='hopping file (ftype 1)')
    args = ap.parse_args()
    if args.profile is None:
        from aiida.storage.sqlite_temp import SqliteTempBackend
        aiida.load_profile(SqliteTempBackend.create_profile('pyrpa-example'), allow_switch=True)
    else:
        aiida.load_profile(args.profile)

    if args.calc == 'flex':
        res, node = pa.flex_workflow.run_get_node(**model_inputs(args.hop))
        sigma = res['sigma'].get_array('sigma')
        print(f'bare mu = {res["mu"].value:.6f} eV, FLEX mu = {res["mu_self"].value:.6f} eV')
        print(f'Sigma: shape {sigma.shape}, Sigma_11(k=Gamma, iw_0) = {sigma[0, 0, 0, 0]:.6f}')
    elif args.calc in ('eliashberg', 'eilenberger'):
        kw = dict(model_inputs(args.hop), gap_sym=orm.Int(1), use_flex=orm.Bool(args.flex))   # d_x2-y2
        if args.calc == 'eliashberg':
            res, node = pa.eliashberg_workflow.run_get_node(**kw)
        else:
            res, node = pa.eliashberg_eilenberger_workflow.run_get_node(
                **kw, wc=orm.Float(1.0), coupling=orm.Float(1.0), eil_mesh=orm.List([96, 96, 1]))
        print(f'{"FLEX" if args.flex else "RPA"}: Stoner factor = {res["stoner"].value:.4f}, '
              f'lambda_eliash = {res["lambda_eliash"].value:.4f}')
        chis = res['chi'].get_array('chis_q').real
        print(f'max chi_s(q, v=0) = {chis.max():.4f} at q = {res["chi"].get_array("kmap")[chis.argmax()]}')
        if args.calc == 'eilenberger':
            print('Eilenberger with the Eliashberg gap:', res['eilenberger'].get_dict())
    elif args.calc == 'surface':
        res, node = pa.eilenberger_surface.run_get_node(
            coupling=orm.Float(1.0), temp=orm.Float(0.005), wc=orm.Float(0.5), gap_sym=orm.Str('d'),
            beta_surf=orm.Float(0.7853981633974483))
        d = res['profile'].get_array('delta')
        print(f'Delta(x): surface {abs(d[0]):.4e} eV -> bulk {abs(d[-1]):.4e} eV; files: '
              f'{res["files"].list_object_names()}')
    else:
        res, node = pa.eilenberger_vortex.run_get_node(
            coupling=orm.Float(1.0), temp=orm.Float(0.005), wc=orm.Float(0.5), gap_sym=orm.Str('s'))
        d = res['profile'].get_array('delta')
        print(f'Delta(rho): core {abs(d[0]):.4e} eV -> far {abs(d[-1]):.4e} eV; files: '
              f'{res["files"].list_object_names()}')
    print(f'{node.process_label} pk={node.pk}  state={node.process_state.value}  exit_status={node.exit_status}')
    if args.calc not in ('surface', 'vortex'):
        print_provenance(node)


if __name__ == '__main__':
    main()
