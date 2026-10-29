"""
Minimal AiiDA interface to pyrpa (test sample).

Every step is a calcfunction, so AiiDA records the inputs, the outputs and the links between them
(provenance) while the calculation itself runs in-process through libs.plibs / libs.flibs; no
computer, code or daemon has to be configured.  A CalcJob that submits main.py to a cluster is the
next step and is not here.

Model
  load_hoppings          SinglefileData (hopping file; ftype=4: FolderData of the ecalj job_mlo
                         binaries)                            -> ArrayData rvec, ham_r [, S_r]
                         An MLO basis (ftype=3, 4) carries its overlap S_r, normalized orbitals as in
                         main.py, and every step below then works with G^-1 = (iw+mu) S - H - Sigma.
  chemical_potential     hoppings + fill, temp, mesh          -> Float mu
RPA / FLEX, with or without SOC.  The optional options (Dict, keys and defaults in OPTION_DEFAULTS)
carry what main.py sets in its header: SOC (soc; slist, invs), the sites (site_prof) or a chi basis of
its own (chiolist, site), orbital-dependent U, J (Umat, Jmat), the Matsubara-axis grid (freq_grid
'matsubara' / 'ir', ir_wmax, ir_eps; sw_tail), the FLEX loop (eps, pp, m_diis, scf_loop, sub_sigma,
sw_rescale, require_converged) and what is stored (store_gap, store_sigma)
  flex_self_energy       hoppings + mu, fill, temp, U, J, mesh, Nw [, options]
                                                              -> sigma (ArrayData), mu_self, flex_info, runtime
  rpa_susceptibility     hoppings + mu, temp, U, J, mesh, Nw [, sigma, options]
                                                              -> chi (chi_s, chi_c at v=0), stoner, runtime
  linearized_eliashberg  hoppings + mu, temp, U, J, mesh, Nw, gap_sym [, sigma, options]
                                                              -> lambda_eliash, stoner, [gap,] gap_wannier, runtime
  flex_eliashberg        FLEX -> chi -> lambda in one step, Sigma kept in memory only (store_sigma=False)
                                                              -> mu_self, flex_info, chi, stoner, lambda_eliash, ...
  With SOC: chi holds the spin channels chiszz, chispm, chisxx, chisyy and chic; no Stoner factor is
  given (the SOC chi0 has none, as in main.py); the Eliashberg parity must be diagonal.
  flex_info (Dict): iterations, final esterr, converged, Stoner >= 1 warning of the FLEX loop;
  runtime (Dict): OpenMP threads and wait policy, libfmod build (OpenMP runtime), numpy, git commit
Quasiclassical Eilenberger
  eilenberger_bulk       hoppings + mu, temp, wc, coupling, gap_sym, mesh [, gap_wannier]
                         homogeneous multi-orbital gap on the Fermi surface of the model; with
                         gap_wannier the RPA/FLEX gap is the pairing form factor
                                                              -> result (Dict), f_avg (ArrayData)
  eilenberger_surface    coupling, temp, wc, gap_sym ('s', 'd', 'dxy'), beta_surf, ...  (model FS)
                                                              -> profile Delta(x), files
  eilenberger_vortex     coupling, temp, wc, gap_sym, field, ...  (model FS)
                                                              -> profile, files
Workflows (workfunctions chaining the above; all take options)
  flex_workflow                  hopping file -> mu -> FLEX Sigma
  eliashberg_workflow            hopping file -> mu -> [FLEX] -> chi, lambda, gap
                                 (FLEX with store_sigma=False: flex_eliashberg, Sigma not stored)
  eliashberg_eilenberger_workflow  eliashberg_workflow -> eilenberger_bulk with its gap

Arrays keep pyrpa's layout (Sigma, gap [Norb, Norb, Nw, Nk_irr] on the irreducible k-mesh of
gen_irr_k_TRS), complex128 as is.  Files a driver writes (LDOS, profiles) are kept as FolderData.
"""
import contextlib
import os
import pathlib
import re
import subprocess
import sys
import tempfile

import numpy as np
from aiida import orm
from aiida.engine import calcfunction, workfunction

# pyrpa is not an installed package: put the repository root on the path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
import libs.flibs as flibs  # noqa: E402
import libs.plibs as plibs  # noqa: E402


# --------------------------------------------------------------------------- #
#  helpers (plain Python, not recorded)
# --------------------------------------------------------------------------- #
def _model(hoppings: orm.ArrayData) -> tuple[np.ndarray, np.ndarray, np.ndarray | list]:
    """rvec, ham_r and S_r ([] for an orthogonal basis, as main.py passes it)."""
    S_r = hoppings.get_array('S_r') if 'S_r' in hoppings.get_arraynames() else []
    return hoppings.get_array('rvec'), hoppings.get_array('ham_r'), S_r


def _mesh(mesh: orm.List) -> tuple[int, int, int]:
    Nx, Ny, Nz = (int(n) for n in mesh.get_list())
    return Nx, Ny, Nz


def _arrays(**arrays: np.ndarray) -> orm.ArrayData:
    out = orm.ArrayData()
    for name, a in arrays.items():
        out.set_array(name, np.asarray(a))
    return out


def _npz_to_arraydata(path: str) -> orm.ArrayData:
    with np.load(path) as d:
        return _arrays(**{k: d[k] for k in d.files})


def _arraydata_to_npz(data: orm.ArrayData, path: str) -> None:
    np.savez(path, **{k: data.get_array(k) for k in data.get_arraynames()})


@contextlib.contextmanager
def _workdir():
    """Run a driver that writes files into the working directory inside a scratch folder."""
    with tempfile.TemporaryDirectory() as d, contextlib.chdir(d):
        yield d


def _folder(path: str) -> orm.FolderData:
    return orm.FolderData(tree=path)


OPTION_DEFAULTS = {
    # interaction (main.py: site_prof, orb_dep with Umat/Jmat, a user chiolist with its site list)
    'site_prof': None,          # orbital counts per site, e.g. [5, 5]; None: one site holding every orbital
    'Umat': None, 'Jmat': None,  # orbital-dependent U, J (Norb x Norb, or one site's worth); None: uniform U, J
    'chiolist': None, 'site': None,   # chi orbital-pair basis (1-based pairs) and its site labels
    # spin-orbit coupling (main.py: sw_soc, slist, invs): the model is spinful and U, J (or Umat, Jmat,
    # site_prof) refer to its spin-orbitals; spin labels default to main.py's (all up, then all down)
    'soc': False,
    'slist': None,              # spin label (+1 / -1) of each orbital
    'invs': None,               # 1-based opposite-spin partner of each orbital
    # Matsubara-axis grid (main.py: freq_grid, ir_wmax, ir_eps, sw_chi0_tail)
    'freq_grid': 'matsubara',   # 'ir': sparse sampling of the IR basis (Nw is then the basis's)
    'ir_wmax': None,            # IR spectral cutoff [eV]; None: 2 x the largest |e_k - mu|, as main.py
    'ir_eps': 1.0e-8,
    'sw_tail': False,           # tail-corrected chi0 (uniform grid; ignored on the IR grid)
    # FLEX loop (flibs.mkself)
    'eps': 1.0e-4, 'pp': 0.3, 'm_diis': 5, 'scf_loop': 300, 'sub_sigma': 1, 'sw_rescale': False,
    'require_converged': True,  # an unconverged FLEX loop fails the step instead of returning its Sigma
    # storage
    'store_gap': True,          # linearized_eliashberg: the full Delta(k, iw_n) besides gap_wannier
    'store_sigma': True,        # eliashberg_workflow with FLEX: store Sigma (False: flex_eliashberg)
}


def _opts(options: orm.Dict | None) -> dict:
    """The options with their defaults filled in; an unknown key is an error, not silently ignored."""
    opts = dict(OPTION_DEFAULTS)
    if options is not None:
        given = options.get_dict()
        unknown = set(given) - set(OPTION_DEFAULTS)
        if unknown:
            raise ValueError(f'unknown options {sorted(unknown)}; known: {sorted(OPTION_DEFAULTS)}')
        opts.update(given)
    if opts['freq_grid'] not in ('matsubara', 'ir'):
        raise ValueError(f"freq_grid={opts['freq_grid']!r} must be 'matsubara' or 'ir'")
    if (opts['Umat'] is None) != (opts['Jmat'] is None):
        raise ValueError('Umat and Jmat are given together')
    if (opts['chiolist'] is None) != (opts['site'] is None):
        raise ValueError('a user chiolist needs its site list (and vice versa)')
    if opts['soc'] and opts['sw_tail']:
        raise ValueError('sw_tail (tail-corrected chi0) exists for the no-SOC path only')
    return opts


def _normal_state(hoppings: orm.ArrayData, mesh: orm.List, U: orm.Float, J: orm.Float, opts: dict) -> dict:
    """Irreducible k-mesh, bands and the U, J vertex shared by the RPA/FLEX steps (as main.py builds them)."""
    rvec, ham_r, S_r = _model(hoppings)
    Nx, Ny, Nz = _mesh(mesh)
    Norb = ham_r.shape[1]
    klist, kmap, invk = flibs.gen_irr_k_TRS(Nx, Ny, Nz)
    if len(S_r):
        hamk, ovlk = flibs.gen_ham(klist, ham_r, rvec, Ovl_r=S_r)
        eig, uni = flibs.get_eig(hamk, ovlk)
    else:
        hamk, ovlk = flibs.gen_ham(klist, ham_r, rvec), None
        eig, uni = flibs.get_eig(hamk)
    site_prof = opts['site_prof'] if opts['site_prof'] is not None else [Norb]
    if opts['chiolist'] is not None:
        olist = np.asfortranarray(opts['chiolist'], dtype=np.int64)
        site = np.asarray(opts['site'], dtype=np.int64)
    else:
        olist, site = plibs.get_chi_orb_list(Norb, site_prof)
    st = dict(rvec=rvec, ham_r=ham_r, Nx=Nx, Ny=Ny, Nz=Nz, Norb=Norb, klist=klist, kmap=kmap, invk=invk,
              hamk=hamk, ovlk=ovlk, eig=eig, uni=uni, olist=olist, soc=bool(opts['soc']))
    UJ = None
    if opts['Umat'] is not None:
        UJ = plibs.expand_UJ_sites(np.asarray(opts['Umat'], dtype=np.float64),
                                   np.asarray(opts['Jmat'], dtype=np.float64), Norb, site_prof)
    if st['soc']:
        slist, invs = _spin_lists(Norb, opts)
        st.update(slist=slist, invs=invs,
                  Vmat=(flibs.gen_Vmatrix_orb(olist, slist, site, invs, *UJ) if UJ is not None else
                        flibs.gen_Vmatrix(olist, slist, site, invs, U.value, J.value)))
    else:
        st['Smat'], st['Cmat'] = (flibs.gen_SCmatrix_orb(olist, site, *UJ) if UJ is not None else
                                  flibs.gen_SCmatrix(olist, site, U=U.value, J=J.value))
    return st


def _spin_lists(Norb: int, opts: dict) -> tuple[np.ndarray, np.ndarray]:
    """slist, invs of a spinful basis: the options' if given, else main.py's default ordering
    (all up, then all down; the up block takes the extra orbital when Norb is odd)."""
    slist = opts['slist']
    if slist is None:
        slist = np.ones(Norb, dtype=np.int64)
        slist[(Norb+1)//2:] = -1
    invs = opts['invs']
    if invs is None:
        invs = np.concatenate([np.arange((Norb+1)//2, Norb), np.arange((Norb+1)//2)])+1
    slist, invs = np.asarray(slist, dtype=np.int64), np.asarray(invs, dtype=np.int64)
    if slist.shape != (Norb,) or invs.shape != (Norb,):
        raise ValueError(f'slist / invs need one entry per orbital (Norb = {Norb})')
    return slist, invs


@contextlib.contextmanager
def _frequency_grid(opts: dict, temp: float, mu: float, st: dict, sigma: orm.ArrayData | None = None):
    """
    Install the Matsubara-axis grid of a step and yield its Nw (None: the uniform grid, use the Nw
    input).  A FLEX sigma records its grid (ir_wmax, ir_eps), and a step built on it uses that same
    grid, so the auto cutoff is not re-derived from a shifted mu.
    """
    on_ir = sigma is not None and 'ir_wmax' in sigma.get_arraynames()
    if sigma is not None and on_ir != (opts['freq_grid'] == 'ir'):
        raise ValueError(f"sigma was computed on the {'IR' if on_ir else 'uniform'} grid, "
                         f"options ask for freq_grid={opts['freq_grid']!r}")
    if opts['freq_grid'] == 'matsubara':
        yield None
        return
    if on_ir:
        wmax, eps = float(sigma.get_array('ir_wmax')), float(sigma.get_array('ir_eps'))
    else:
        wmax = opts['ir_wmax'] if opts['ir_wmax'] is not None else 2.0*np.abs(st['eig']-mu).max()
        eps = opts['ir_eps']
    grid = plibs.ir_grid(temp, wmax, eps=eps)
    plibs.set_freq_grid(grid)
    try:
        grid['ir_wmax'], grid['ir_eps'] = wmax, eps
        yield grid
    finally:
        plibs.set_freq_grid(None)


def _grid_arrays(grid: dict | None) -> dict:
    """The arrays that let a later step rebuild the IR grid of this one."""
    return {} if grid is None else {'ir_wmax': np.array(grid['ir_wmax']), 'ir_eps': np.array(grid['ir_eps'])}


@contextlib.contextmanager
def _fortran_log():
    """Capture what the Fortran kernels print (file descriptor 1) and echo it afterwards: the FLEX
    loop reports its progress only there.  Yields a dict whose 'text' is filled on exit."""
    log = {}
    sys.stdout.flush()
    saved = os.dup(1)
    with tempfile.TemporaryFile(mode='w+b') as f:
        os.dup2(f.fileno(), 1)
        try:
            yield log
        finally:
            sys.stdout.flush()
            os.dup2(saved, 1)
            os.close(saved)
            f.seek(0)
            log['text'] = f.read().decode(errors='replace')
    sys.stdout.write(log['text'])
    sys.stdout.flush()


def _flex_info(text: str, opts: dict) -> dict:
    """Convergence of the FLEX loop from its log: the loop exits once esterr < eps.  A loop that ends with
    a Stoner factor >= 1 has settled on a magnetically unstable, unphysical fixed point and counts as not
    converged (near the instability the iteration can be thrown there by its first overshoots)."""
    esterr = [float(x) for x in re.findall(r'esterr=\s*([-+0-9.Ee]+)', text)]
    info = {'iterations': len(re.findall(r'iter=', text)), 'esterr': esterr[-1] if esterr else None,
            'eps': opts['eps'], 'stoner_unstable': 'Stoner factor >= 1' in text}
    info['converged'] = (info['esterr'] is not None and info['esterr'] < opts['eps']
                         and not info['stoner_unstable'])
    return info


def _runtime() -> orm.Dict:
    """How this step ran: OpenMP threads and wait policy, the libfmod build, numpy and the git commit."""
    lib = pathlib.Path(flibs._loader._lib._name)
    try:
        gomp = b'libgomp.so' in lib.read_bytes()
    except OSError:
        gomp = None
    try:
        commit = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd=lib.parent, capture_output=True,
                                text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        commit = None
    nthreads, _ = flibs.omp_params()
    return orm.Dict({'omp_threads': int(nthreads), 'omp_wait_policy': os.environ.get('OMP_WAIT_POLICY'),
                     'libfmod': lib.name, 'openmp_runtime': None if gomp is None else ('libgomp' if gomp else 'other'),
                     'numpy': np.__version__, 'pyrpa_commit': commit})


def _inversion(st: dict, tol: float = 3.0e-2) -> tuple[np.ndarray, np.ndarray | None]:
    """Parity of the Eliashberg solver as main.py sets it: the diagonal Wannier parity, or the
    signed orbital permutation when inversion exchanges sites (plibs.inversion_op, tol = main.py's
    inv_tol).  Returns (plist, iperm); iperm is None for the diagonal form."""
    plist = flibs.get_plist(st['rvec'], st['ham_r'])
    psgn, pperm, pstat, _ = plibs.inversion_op(st['rvec'], st['ham_r'], plist, tol=tol)
    if pstat == 'monomial':
        return psgn, pperm
    if pstat != 'diagonal':
        raise ValueError(f'inversion is not a signed permutation of the orbitals ({pstat}): '
                         'Delta(-k) cannot be rebuilt from Delta(k)')
    return plist, None


def _green(st: dict, mu: float, temp: float, Nw: int, sigma: orm.ArrayData | None) -> np.ndarray:
    """Bare G, or the FLEX one when sigma (output of flex_self_energy, with mu = its mu_self) is given."""
    if sigma is None:
        return flibs.gen_Green0(st['eig'], st['uni'], mu, temp, Nw)
    sigmak = sigma.get_array('sigma')
    if sigmak.shape[2] != Nw or sigmak.shape[3] != len(st['klist']):
        raise ValueError(f'sigma has shape {sigmak.shape}, not Nw = {Nw} on {len(st["klist"])} k-points')
    return flibs.gen_green(sigmak, st['hamk'], mu, temp, overlap=st['ovlk'])


# --------------------------------------------------------------------------- #
#  model
# --------------------------------------------------------------------------- #
@calcfunction
def load_hoppings(hop_file: orm.SinglefileData | orm.FolderData, ftype: orm.Int) -> orm.ArrayData:
    """
    Read the model as main.py does for this ftype: plibs.import_hoppings, or for an MLO basis
    import_MLO_hoppings (ftype=3, a file) / read_HamRsMLO (ftype=4, a FolderData holding HamRsMLO
    and HamiltonianPMTInfo), which also give the overlap S_r of the normalized orbitals.
    """
    # the readers take a path: copy the stored file(s) out of the repository
    with tempfile.TemporaryDirectory() as d:
        if isinstance(hop_file, orm.FolderData):
            hop_file.base.repository.copy_tree(d)
            path = d
        else:
            path = os.path.join(d, hop_file.filename)
            with hop_file.open(mode='rb') as src, open(path, 'wb') as dst:
                dst.write(src.read())
        S_r = None
        if ftype.value == 3:
            rvec, ham_r, S_r, _, _ = plibs.import_MLO_hoppings(path)
        elif ftype.value == 4:
            rvec, ham_r, S_r, _, _ = plibs.read_HamRsMLO(path)
        else:
            rvec, ham_r, _, _ = plibs.import_hoppings(path, ftype.value)
    arrays = dict(rvec=np.asarray(rvec, dtype=np.float64), ham_r=np.asarray(ham_r, dtype=np.complex128))
    if S_r is not None:
        arrays['S_r'] = np.asarray(S_r, dtype=np.complex128)
    return _arrays(**arrays)


@calcfunction
def chemical_potential(hoppings: orm.ArrayData, fill: orm.Float, temp: orm.Float, mesh: orm.List) -> orm.Float:
    """Chemical potential of the bare bands at filling fill on the full mesh (plibs.calc_mu)."""
    rvec, ham_r, S_r = _model(hoppings)
    Nx, Ny, Nz = _mesh(mesh)
    Nk, klist = plibs.gen_klist(Nx, Ny, Nz, sw_pp=False)
    eig = plibs.get_eigs(klist, ham_r, S_r, rvec)[0]
    return orm.Float(plibs.calc_mu(eig, Nk, fill.value, temp.value))


# --------------------------------------------------------------------------- #
#  RPA / FLEX
# --------------------------------------------------------------------------- #
def _run_flex(st: dict, mu: float, fill: float, temp: float, Nw: int, opts: dict) -> tuple:
    """FLEX loop (flibs.mkself, mkself_soc) on the installed grid; returns (Sigma, mu_self, flex_info)."""
    loop = dict(sub_sigma=opts['sub_sigma'], scf_loop=opts['scf_loop'], eps=opts['eps'], pp=opts['pp'],
                m_diis=opts['m_diis'], sw_rescale=opts['sw_rescale'], overlap=st['ovlk'])
    with _fortran_log() as log:
        if st['soc']:
            sigmak, mu_self = flibs.mkself_soc(st['Vmat'], st['kmap'], st['invk'], st['invs'], st['olist'],
                                               st['slist'], st['hamk'], st['eig'], st['uni'], mu, fill, temp, Nw,
                                               st['Nx'], st['Ny'], st['Nz'], False, False, **loop)
        else:
            sigmak, mu_self = flibs.mkself(st['Smat'], st['Cmat'], st['kmap'], st['invk'], st['olist'], st['hamk'],
                                           st['eig'], st['uni'], mu, fill, temp, Nw, st['Nx'], st['Ny'], st['Nz'],
                                           False, False, sw_tail=opts['sw_tail'], **loop)
    info = _flex_info(log['text'], opts)
    if opts['require_converged'] and not info['converged']:
        why = ('ended with a Stoner factor >= 1 (unphysical fixed point)' if info['stoner_unstable'] else
               f"did not converge in {info['iterations']} iterations (esterr = {info['esterr']}, eps = {opts['eps']})")
        raise RuntimeError(f'FLEX {why}; lower pp, raise scf_loop, tighten ir_eps on the IR grid, '
                           'or set require_converged=False')
    return sigmak, mu_self, info


def _chi0(st: dict, Gk: np.ndarray, mu: float, temp: float, opts: dict) -> tuple:
    """
    Irreducible chi0, the Stoner factor and, with SOC, the sign tables of the SOC pairing kernel
    (sgnsig, sgnsig2, invschi).  sw_tail: tail-corrected chi0 as main.py's sw_chi0_tail.  The SOC
    chi0 gives no Stoner factor (main.py does not report one either): it is None there.
    """
    if st['soc']:
        chi0, *tables = flibs.get_chi0_soc(st['Vmat'], Gk, st['olist'], st['slist'], st['kmap'], st['invk'],
                                           st['invs'], temp, st['Nx'], st['Ny'], st['Nz'])
        return chi0, None, tables
    if opts['sw_tail']:
        chi0, stoner = flibs.get_chi0_tail(st['Smat'], st['Cmat'], Gk, st['eig'], st['uni'], st['olist'],
                                           st['kmap'], st['invk'], mu, temp, st['Nx'], st['Ny'], st['Nz'])
    else:
        chi0, stoner = flibs.get_chi0(st['Smat'], st['Cmat'], Gk, st['olist'], st['kmap'], st['invk'], temp,
                                      st['Nx'], st['Ny'], st['Nz'])
    return chi0, stoner, None


def _stoner(stoner: float | None) -> dict:
    return {} if stoner is None else {'stoner': orm.Float(stoner)}


def _chi_arrays(st: dict, chi0: np.ndarray) -> orm.ArrayData:
    """chi_s, chi_c at v = 0 and their leading eigenvalues on the full mesh (order of kmap); with SOC the
    spin channels zz, +-, xx, yy (S = sigma/2) as calc_lin_eliash_soc writes them, chis_q being zz."""
    lead = lambda c: flibs.get_eig_or_tr_chi(c, st['invk'], True)
    if st['soc']:
        chic, chiszz, chispm, chisxx, chisyy = flibs.get_chis_chic_soc(chi0, st['Vmat'], st['olist'], st['slist'],
                                                                       st['invs'])
        return _arrays(qlist=st['klist'], kmap=st['kmap'], chic=chic, chiszz=chiszz, chispm=chispm,
                       chisxx=chisxx, chisyy=chisyy, chic_q=lead(chic), chis_q=lead(chiszz),
                       chispm_q=lead(chispm), chisxx_q=lead(chisxx), chisyy_q=lead(chisyy))
    chis, chic = flibs.get_chis_chic(chi0, st['Smat'], st['Cmat'])
    return _arrays(qlist=st['klist'], chis=chis, chic=chic, chis_q=lead(chis), chic_q=lead(chic), kmap=st['kmap'])


def _eliashberg(st: dict, chi0: np.ndarray, tables: list | None, Gk: np.ndarray, temp: float, Nw: int,
                gap_sym: int, opts: dict) -> dict:
    """Linearized Eliashberg equation on chi0 (overwritten); lambda, gap_wannier and, with store_gap, the gap."""
    plist, iperm = _inversion(st)
    init_delta = plibs.get_initial_gap(st['klist'], st['Norb'], gap_sym)
    if st['soc']:
        if iperm is not None:   # as main.py: the SOC kernels take a diagonal parity only
            raise ValueError('inversion exchanges orbitals (signed permutation), which the SOC Eliashberg '
                             'path does not support (diagonal parity only)')
        sgnsig, sgnsig2, invschi = tables
        gap, lam = flibs.linearized_eliashberg_soc(chi0, Gk, st['uni'], init_delta, st['Vmat'], sgnsig, sgnsig2,
                                                   plist, st['slist'], st['olist'], st['kmap'], st['invk'],
                                                   st['invs'], invschi, st['Nx'], st['Ny'], st['Nz'], temp, gap_sym,
                                                   overlap=st['ovlk'])
    else:
        gap, lam = flibs.linearized_eliashberg(chi0, Gk, st['uni'], init_delta, st['Smat'], st['Cmat'],
                                               st['olist'], plist, st['kmap'], st['invk'], st['Nx'], st['Ny'],
                                               st['Nz'], temp, gap_sym, iperm=iperm, overlap=st['ovlk'])
    with _workdir():
        plibs.output_gap_wannier(gap, st['kmap'], st['invk'], st['Nx'], st['Ny'], st['Nz'], Nw, temp)
        out = {'lambda_eliash': orm.Float(lam), 'gap_wannier': _npz_to_arraydata('gap_wannier.npz')}
    if opts['store_gap']:
        out['gap'] = _arrays(gap=gap, klist=st['klist'], iw=plibs.fermion_frequencies(temp, Nw))
    return out


@calcfunction
def flex_self_energy(hoppings: orm.ArrayData, mu: orm.Float, fill: orm.Float, temp: orm.Float,
                     U: orm.Float, J: orm.Float, mesh: orm.List, Nw: orm.Int, options: orm.Dict | None = None) -> dict:
    """FLEX self-energy (flibs.mkself) on the uniform Matsubara or the IR grid."""
    opts = _opts(options)
    st = _normal_state(hoppings, mesh, U, J, opts)
    with _frequency_grid(opts, temp.value, mu.value, st) as grid:
        nw = Nw.value if grid is None else grid['Nw']
        sigmak, mu_self, info = _run_flex(st, mu.value, fill.value, temp.value, nw, opts)
        iw = plibs.fermion_frequencies(temp.value, nw)
    sigma = _arrays(sigma=sigmak, klist=st['klist'], iw=iw, **_grid_arrays(grid))
    return {'sigma': sigma, 'mu_self': orm.Float(mu_self), 'flex_info': orm.Dict(info), 'runtime': _runtime()}


@calcfunction
def rpa_susceptibility(hoppings: orm.ArrayData, mu: orm.Float, temp: orm.Float, U: orm.Float, J: orm.Float,
                       mesh: orm.List, Nw: orm.Int, sigma: orm.ArrayData | None = None,
                       options: orm.Dict | None = None) -> dict:
    """
    RPA spin / charge susceptibility at v = 0 (leading eigenvalue of the orbital matrix) on the
    irreducible q-mesh; with sigma (and mu = mu_self) the FLEX one.  Also the Stoner factor.
    """
    opts = _opts(options)
    st = _normal_state(hoppings, mesh, U, J, opts)
    with _frequency_grid(opts, temp.value, mu.value, st, sigma) as grid:
        nw = Nw.value if grid is None else grid['Nw']
        Gk = _green(st, mu.value, temp.value, nw, sigma)
        chi0, stoner, _ = _chi0(st, Gk, mu.value, temp.value, opts)
    return {'chi': _chi_arrays(st, chi0), **_stoner(stoner), 'runtime': _runtime()}


@calcfunction
def linearized_eliashberg(hoppings: orm.ArrayData, mu: orm.Float, temp: orm.Float, U: orm.Float, J: orm.Float,
                          mesh: orm.List, Nw: orm.Int, gap_sym: orm.Int, sigma: orm.ArrayData | None = None,
                          options: orm.Dict | None = None) -> dict:
    """
    Linearized Eliashberg equation with the RPA (or, with sigma, FLEX) pairing interaction:
    the leading eigenvalue lambda_eliash, the gap Delta(k, iw_n) (store_gap) and its Wannier-R export
    (output_gap_wannier), which eilenberger_bulk takes as the pairing form factor.
    """
    opts = _opts(options)
    st = _normal_state(hoppings, mesh, U, J, opts)
    with _frequency_grid(opts, temp.value, mu.value, st, sigma) as grid:
        nw = Nw.value if grid is None else grid['Nw']
        Gk = _green(st, mu.value, temp.value, nw, sigma)
        chi0, stoner, tables = _chi0(st, Gk, mu.value, temp.value, opts)
        out = _eliashberg(st, chi0, tables, Gk, temp.value, nw, gap_sym.value, opts)
    return {**out, **_stoner(stoner), 'runtime': _runtime()}


@calcfunction
def flex_eliashberg(hoppings: orm.ArrayData, mu: orm.Float, fill: orm.Float, temp: orm.Float, U: orm.Float,
                    J: orm.Float, mesh: orm.List, Nw: orm.Int, gap_sym: orm.Int,
                    options: orm.Dict | None = None) -> dict:
    """
    FLEX self-energy, the FLEX susceptibility and the linearized Eliashberg equation in one step:
    Sigma (Norb^2 x Nw x Nk complex numbers) stays in memory and only the results are stored.
    """
    opts = _opts(options)
    st = _normal_state(hoppings, mesh, U, J, opts)
    with _frequency_grid(opts, temp.value, mu.value, st) as grid:
        nw = Nw.value if grid is None else grid['Nw']
        sigmak, mu_self, info = _run_flex(st, mu.value, fill.value, temp.value, nw, opts)
        Gk = flibs.gen_green(sigmak, st['hamk'], mu_self, temp.value, overlap=st['ovlk'])
        chi0, stoner, tables = _chi0(st, Gk, mu_self, temp.value, opts)
        chi = _chi_arrays(st, chi0)
        out = _eliashberg(st, chi0, tables, Gk, temp.value, nw, gap_sym.value, opts)
    return {**out, 'mu_self': orm.Float(mu_self), 'flex_info': orm.Dict(info), 'chi': chi,
            **_stoner(stoner), 'runtime': _runtime()}


# --------------------------------------------------------------------------- #
#  quasiclassical Eilenberger
# --------------------------------------------------------------------------- #
@calcfunction
def eilenberger_bulk(hoppings: orm.ArrayData, mu: orm.Float, temp: orm.Float, wc: orm.Float, coupling: orm.Float,
                     gap_sym: orm.Int, mesh: orm.List, fs_width: orm.Float = lambda: orm.Float(5.0e-3),
                     imp_gamma: orm.Float = lambda: orm.Float(0.0), imp_c: orm.Float = lambda: orm.Float(1.0e8),
                     find_tc: orm.Bool = lambda: orm.Bool(False),
                     gap_wannier: orm.ArrayData | None = None, gap_iw: orm.Int = lambda: orm.Int(0)) -> dict:
    """
    Homogeneous multi-orbital Eilenberger equation (the core of plibs.calc_eilenberger): the
    linearized pairing eigenvalue lambda_lin(T), the self-consistent gap amplitude Delta and the
    FS-averaged anomalous propagator <f>(w_n) at that gap; optionally Tc.  gap_wannier (from
    linearized_eliashberg) replaces the gap_sym harmonic by that gap at Matsubara slice gap_iw,
    projected band-resolved onto the Fermi surface.
    """
    rvec, ham_r, S_r = _model(hoppings)
    Nx, Ny, Nz = _mesh(mesh)
    _, klist, eig, uni, _ = plibs.get_emesh(Nx, Ny, Nz, ham_r, S_r, rvec, np.eye(3), sw_uni=True)
    gorb = None
    if gap_wannier is not None:
        with _workdir():
            _arraydata_to_npz(gap_wannier, 'gap_wannier.npz')
            gorb = plibs.gap_orbital_from_wannier('gap_wannier', gap_iw.value)
    wf, phif, _ = plibs.build_fs(eig, klist, mu.value, gap_sym.value, fs_width.value, uni=uni, gap_orbital=gorb,
                                 sw_band=True)
    if len(wf) == 0:
        raise ValueError('eilenberger_bulk: no Fermi-surface points (check mu / fs_width)')
    T, g, c, lam_c = temp.value, imp_gamma.value, imp_c.value, coupling.value
    omega = plibs.matsubara(T, wc.value)
    lam = plibs.pairing_eigenvalue(T, wf, phif, omega, g, c, lam_c)
    damp = plibs.solve_gap(T, wf, phif, omega, g, c, lam_c)
    _, f = plibs.solve_quasiclassical(omega, phif*damp, wf, g, c)
    result = {'lambda_lin': float(lam), 'delta': float(damp), 'n_fs': int(len(wf)),
              'form_factor': 'gap_wannier' if gorb is not None else f'gap_sym={gap_sym.value}'}
    if find_tc.value:
        result['tc'] = float(plibs.find_tc(wf, phif, g, c, lam_c, wc.value, 1.0e-4, 5.0*T))
    return {'result': orm.Dict(result), 'runtime': _runtime(),
            'f_avg': _arrays(omega=omega, f_avg=np.tensordot(wf, f, axes=(0, 0))/wf.sum())}


@calcfunction
def eilenberger_surface(coupling: orm.Float, temp: orm.Float, wc: orm.Float, gap_sym: orm.Str,
                        beta_surf: orm.Float = lambda: orm.Float(0.0), sw_ldos: orm.Bool = lambda: orm.Bool(True),
                        imp_gamma: orm.Float = lambda: orm.Float(0.0), imp_c: orm.Float = lambda: orm.Float(1.0e8),
                        h: orm.Float = lambda: orm.Float(0.0), Nbeta: orm.Int = lambda: orm.Int(30),
                        Lxi: orm.Float = lambda: orm.Float(8.0), nper: orm.Int = lambda: orm.Int(16)) -> dict:
    """Self-consistent gap profile Delta(x) at a specular surface (+ LDOS files), plibs.calc_surface, model FS."""
    with _workdir() as d:
        x, damp = plibs.calc_surface(coupling.value, temp.value, wc.value, gap_sym=gap_sym.value,
                                     beta_surf=beta_surf.value, sw_ldos=sw_ldos.value, imp_gamma=imp_gamma.value,
                                     imp_c=imp_c.value, h=h.value, Nbeta=Nbeta.value, Lxi=Lxi.value,
                                     nper=nper.value)
        files = _folder(d)
    return {'profile': _arrays(x=x, delta=damp), 'files': files}


@calcfunction
def eilenberger_vortex(coupling: orm.Float, temp: orm.Float, wc: orm.Float, gap_sym: orm.Str,
                       field: orm.Float = lambda: orm.Float(0.0), sw_ldos: orm.Bool = lambda: orm.Bool(True),
                       imp_gamma: orm.Float = lambda: orm.Float(0.0), imp_c: orm.Float = lambda: orm.Float(1.0e8),
                       h: orm.Float = lambda: orm.Float(0.0), Lxi: orm.Float = lambda: orm.Float(8.0),
                       ngrid: orm.Int = lambda: orm.Int(81)) -> dict:
    """
    Self-consistent vortex (field = 0: isolated, > 0: B/Hc2 of a circular-cell lattice) + LDOS
    files, plibs.calc_vortex on the model FS.  profile holds the solver's (grid, Delta): the radial
    grid and Delta(rho) of the clean s-wave axisymmetric solver, else the 2D grid and Delta(x, y).
    """
    with _workdir() as d:
        grid, prof = plibs.calc_vortex(coupling.value, temp.value, wc.value, gap_sym=gap_sym.value,
                                       sw_ldos=sw_ldos.value, imp_gamma=imp_gamma.value, imp_c=imp_c.value,
                                       field=field.value, h=h.value, Lxi=Lxi.value, ngrid=ngrid.value)
        files = _folder(d)
    return {'profile': _arrays(grid=grid, delta=prof), 'files': files}


# --------------------------------------------------------------------------- #
#  workflows
# --------------------------------------------------------------------------- #
@workfunction
def flex_workflow(hop_file: orm.SinglefileData | orm.FolderData, ftype: orm.Int, fill: orm.Float, temp: orm.Float,
                  U: orm.Float, J: orm.Float, mesh: orm.List, Nw: orm.Int, options: orm.Dict | None = None) -> dict:
    """hopping file -> bare mu -> FLEX self-energy, recorded as one provenance graph."""
    hoppings = load_hoppings(hop_file, ftype)
    mu = chemical_potential(hoppings, fill, temp, mesh)
    flex = flex_self_energy(hoppings, mu, fill, temp, U, J, mesh, Nw, **_given(options))
    return {'hoppings': hoppings, 'mu': mu, 'sigma': flex['sigma'], 'mu_self': flex['mu_self'],
            'flex_info': flex['flex_info']}


def _given(options: orm.Dict | None) -> dict:
    """options as a keyword of the steps, left out when not given (the steps then use the defaults)."""
    return {} if options is None else {'options': options}


@workfunction
def eliashberg_workflow(hop_file: orm.SinglefileData | orm.FolderData, ftype: orm.Int, fill: orm.Float,
                        temp: orm.Float, U: orm.Float, J: orm.Float, mesh: orm.List, Nw: orm.Int, gap_sym: orm.Int,
                        use_flex: orm.Bool = lambda: orm.Bool(False), options: orm.Dict | None = None) -> dict:
    """
    hopping file -> mu -> [FLEX Sigma] -> chi_s/chi_c and the linearized Eliashberg gap (RPA or FLEX).
    With FLEX and store_sigma=False the three steps run as one (flex_eliashberg) and Sigma is not stored.
    """
    hoppings = load_hoppings(hop_file, ftype)
    mu = chemical_potential(hoppings, fill, temp, mesh)
    out = {'hoppings': hoppings, 'mu': mu}
    opt = _given(options)
    if use_flex.value and not _opts(options)['store_sigma']:
        el = flex_eliashberg(hoppings, mu, fill, temp, U, J, mesh, Nw, gap_sym, **opt)
        out.update({k: el[k] for k in ('mu_self', 'flex_info', 'chi', 'stoner', 'lambda_eliash', 'gap_wannier', 'gap')
                    if k in el})
        return out
    extra = {}
    if use_flex.value:
        flex = flex_self_energy(hoppings, mu, fill, temp, U, J, mesh, Nw, **opt)
        out.update(sigma=flex['sigma'], mu_self=flex['mu_self'], flex_info=flex['flex_info'])
        mu, extra = flex['mu_self'], {'sigma': flex['sigma']}
    chi = rpa_susceptibility(hoppings, mu, temp, U, J, mesh, Nw, **extra, **opt)
    el = linearized_eliashberg(hoppings, mu, temp, U, J, mesh, Nw, gap_sym, **extra, **opt)
    out.update(chi=chi['chi'], lambda_eliash=el['lambda_eliash'], gap_wannier=el['gap_wannier'])
    out.update({k: v for k, v in (('stoner', chi.get('stoner')), ('gap', el.get('gap'))) if v is not None})
    return out


@workfunction
def eliashberg_eilenberger_workflow(hop_file: orm.SinglefileData | orm.FolderData, ftype: orm.Int, fill: orm.Float,
                                    temp: orm.Float, U: orm.Float, J: orm.Float, mesh: orm.List, Nw: orm.Int,
                                    gap_sym: orm.Int, wc: orm.Float, coupling: orm.Float, eil_mesh: orm.List,
                                    use_flex: orm.Bool = lambda: orm.Bool(False),
                                    options: orm.Dict | None = None) -> dict:
    """
    eliashberg_workflow, then the homogeneous Eilenberger gap with the RPA/FLEX gap as the pairing
    form factor (on the bare-band Fermi surface at the bare mu, on eil_mesh).
    """
    el = eliashberg_workflow(hop_file, ftype, fill, temp, U, J, mesh, Nw, gap_sym, use_flex, **_given(options))
    eil = eilenberger_bulk(el['hoppings'], el['mu'], temp, wc, coupling, gap_sym, eil_mesh,
                           gap_wannier=el['gap_wannier'])
    return {**el, 'eilenberger': eil['result'], 'f_avg': eil['f_avg']}
