#!/usr/bin/env python
"""
Run the pyrpa AiiDA sample (pyrpa_aiida.py) on inputs/square.hop.

  python examples/aiida/run_example.py [--calc CALC] [--profile NAME] [--hop PATH --ftype N]
                                       [--options JSON] [--mlo-demo] [--soc]

  CALC = flex         hopping file -> mu -> FLEX self-energy
         eliashberg   ... -> chi_s/chi_c, Stoner factor and the linearized Eliashberg lambda, gap
                      (RPA; --flex: with the FLEX self-energy)
         eilenberger  eliashberg, then the homogeneous Eilenberger gap with that gap as the form factor
         surface      Eilenberger surface state (model FS, d-wave [110]: zero-energy bound state)
         vortex       Eilenberger isolated vortex (model FS, s-wave)

--hop / --ftype take any model main.py reads (ftype=3: an MLO file, ftype=4: the directory of the
ecalj job_mlo binaries, passed as FolderData).  --options is a JSON dict of pyrpa_aiida.OPTION_DEFAULTS,
e.g. '{"freq_grid": "ir", "Umat": [[3.0]], "Jmat": [[0.0]], "store_sigma": false}'.  --mlo-demo writes
a small non-orthogonal (MLO, ftype=3) two-orbital model and runs on it with the IR grid.  --soc runs
the spinful inputs/square_2orb_soc.hop (half filling, fill = 2) with soc=True; chi then holds the spin
channels zz, +-, xx, yy and no Stoner factor is given.

The default temporary profile (SqliteTempBackend) needs no database, broker or `verdi presto` and
keeps nothing; --profile NAME stores the nodes in an existing profile.  The AiiDA configuration
folder is still created under $AIIDA_PATH (default ~/.aiida).
"""
import argparse
import json
import os
import sys
import tempfile

import numpy as np

import aiida
from aiida import orm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pyrpa_aiida as pa  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))


def model_inputs(hop: str, ftype: int = 1, fill: float = 0.45, options: dict | None = None) -> dict:
    """Inputs shared by the RPA/FLEX workflows (default: one-orbital square lattice near half filling)."""
    hop_file = orm.FolderData(tree=hop) if os.path.isdir(hop) else orm.SinglefileData(hop)
    inputs = dict(hop_file=hop_file, ftype=orm.Int(ftype), fill=orm.Float(fill), temp=orm.Float(0.02),
                  U=orm.Float(0.3), J=orm.Float(0.0), mesh=orm.List([16, 16, 1]), Nw=orm.Int(128))
    if options:
        inputs['options'] = orm.Dict(options)
    return inputs


def write_mlo_demo(path: str) -> None:
    """A two-orbital square-lattice model in the MLO text format (ftype=3: i j R ReH ImH ReS ImS, H in Ry)
    with a non-orthogonal, unnormalized overlap -- the normalization of main.py is applied on reading."""
    sys.path.insert(0, ROOT)
    import libs.plibs as plibs
    rvec, ham_r, _, _ = plibs.import_hoppings(os.path.join(ROOT, 'inputs', 'square_2orb.hop'), 1)
    if not (np.abs(rvec) < 1e-9).all(axis=1).any():      # the file has no on-site block; S needs one
        rvec = np.vstack([rvec, np.zeros(3)])
        ham_r = np.concatenate([ham_r, np.zeros((1,)+ham_r.shape[1:], dtype=ham_r.dtype)])
    S_r = np.zeros_like(ham_r)
    i0 = int(np.flatnonzero((np.abs(rvec) < 1e-9).all(axis=1))[0])
    S_r[i0] = [[0.3, 0.05], [0.05, 0.8]]                 # on-site norms 0.3, 0.8 and an overlap between them
    for r, R in enumerate(rvec):
        if abs(R).sum() == 1:                            # nearest neighbours
            S_r[r] = [[0.02, 0.01], [0.01, -0.03]]
    with open(path, 'w') as f:
        for a in range(2):
            for b in range(2):
                for r, R in enumerate(rvec):
                    h, s = ham_r[r, a, b]/plibs.RY_TO_EV, S_r[r, a, b]
                    f.write(f'{a+1} {b+1} {R[0]:.0f} {R[1]:.0f} {R[2]:.0f} '
                            f'{h.real:.10f} {h.imag:.10f} {s.real:.10f} {s.imag:.10f}\n')


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
    ap.add_argument('--hop', default=os.path.join(ROOT, 'inputs', 'square.hop'),
                    help='hopping file, or for ftype=4 the job_mlo directory')
    ap.add_argument('--ftype', type=int, default=1, help='format of --hop, as in main.py')
    ap.add_argument('--fill', type=float, default=0.45, help='filling per spin')
    ap.add_argument('--options', default=None, help='JSON dict of pyrpa_aiida.OPTION_DEFAULTS')
    ap.add_argument('--mlo-demo', action='store_true',
                    help='run on a small two-orbital MLO model (written to a temporary file) with the IR grid')
    ap.add_argument('--soc', action='store_true',
                    help='spin-orbit coupling: inputs/square_2orb_soc.hop (unless --hop is given) at fill = 2')
    args = ap.parse_args()
    options = json.loads(args.options) if args.options else None
    if args.soc:
        if args.hop == ap.get_default('hop'):
            args.hop, args.fill = os.path.join(ROOT, 'inputs', 'square_2orb_soc.hop'), 2.0
        options = {'soc': True, **(options or {})}
    if args.mlo_demo:
        args.hop = os.path.join(tempfile.mkdtemp(prefix='pyrpa_mlo_'), 'demo.mlo')
        write_mlo_demo(args.hop)
        args.ftype, args.fill = 3, 0.8
        options = {'freq_grid': 'ir', **(options or {})}
    if args.profile is None:
        from aiida.storage.sqlite_temp import SqliteTempBackend
        aiida.load_profile(SqliteTempBackend.create_profile('pyrpa-example'), allow_switch=True)
    else:
        aiida.load_profile(args.profile)
    inputs = model_inputs(args.hop, args.ftype, args.fill, options)   # nodes need the profile loaded

    if args.calc == 'flex':
        res, node = pa.flex_workflow.run_get_node(**inputs)
        sigma = res['sigma'].get_array('sigma')
        print(f'bare mu = {res["mu"].value:.6f} eV, FLEX mu = {res["mu_self"].value:.6f} eV')
        print('FLEX loop:', res['flex_info'].get_dict())
        print(f'Sigma: shape {sigma.shape}, Sigma_11(k=Gamma, iw_0) = {sigma[0, 0, 0, 0]:.6f}')
    elif args.calc in ('eliashberg', 'eilenberger'):
        kw = dict(inputs, gap_sym=orm.Int(1), use_flex=orm.Bool(args.flex))   # d_x2-y2
        if args.calc == 'eliashberg':
            res, node = pa.eliashberg_workflow.run_get_node(**kw)
        else:
            res, node = pa.eliashberg_eilenberger_workflow.run_get_node(
                **kw, wc=orm.Float(1.0), coupling=orm.Float(1.0), eil_mesh=orm.List([96, 96, 1]))
        stoner = f'Stoner factor = {res["stoner"].value:.4f}, ' if 'stoner' in res else ''
        print(f'{"FLEX" if args.flex else "RPA"}: {stoner}lambda_eliash = {res["lambda_eliash"].value:.4f}')
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
