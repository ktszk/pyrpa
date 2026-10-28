"""
Minimal AiiDA interface to pyrpa (test sample).

Every step is a calcfunction, so AiiDA records the inputs, the outputs and the links between them
(provenance) while the calculation itself runs in-process through libs.plibs / libs.flibs; no
computer, code or daemon has to be configured.  A CalcJob that submits main.py to a cluster is the
next step and is not here.

Model
  load_hoppings          SinglefileData (hopping file)        -> ArrayData rvec, ham_r
  chemical_potential     hoppings + fill, temp, mesh          -> Float mu
RPA / FLEX (no SOC, one site holding every orbital, uniform U, J)
  flex_self_energy       hoppings + mu, fill, temp, U, J, mesh, Nw
                                                              -> sigma (ArrayData), mu_self
  rpa_susceptibility     hoppings + mu, temp, U, J, mesh, Nw [, sigma]
                                                              -> chi (chi_s, chi_c at v=0), stoner
  linearized_eliashberg  hoppings + mu, temp, U, J, mesh, Nw, gap_sym [, sigma]
                                                              -> lambda_eliash, stoner, gap, gap_wannier
Quasiclassical Eilenberger
  eilenberger_bulk       hoppings + mu, temp, wc, coupling, gap_sym, mesh [, gap_wannier]
                         homogeneous multi-orbital gap on the Fermi surface of the model; with
                         gap_wannier the RPA/FLEX gap is the pairing form factor
                                                              -> result (Dict), f_avg (ArrayData)
  eilenberger_surface    coupling, temp, wc, gap_sym ('s', 'd', 'dxy'), beta_surf, ...  (model FS)
                                                              -> profile Delta(x), files
  eilenberger_vortex     coupling, temp, wc, gap_sym, field, ...  (model FS)
                                                              -> profile, files
Workflows (workfunctions chaining the above)
  flex_workflow                  hopping file -> mu -> FLEX Sigma
  eliashberg_workflow            hopping file -> mu -> [FLEX] -> chi, lambda, gap
  eliashberg_eilenberger_workflow  eliashberg_workflow -> eilenberger_bulk with its gap

Arrays keep pyrpa's layout (Sigma, gap [Norb, Norb, Nw, Nk_irr] on the irreducible k-mesh of
gen_irr_k_TRS), complex128 as is.  Files a driver writes (LDOS, profiles) are kept as FolderData.
"""
import contextlib
import os
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
def _model(hoppings: orm.ArrayData) -> tuple[np.ndarray, np.ndarray]:
    return hoppings.get_array('rvec'), hoppings.get_array('ham_r')


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


def _normal_state(hoppings: orm.ArrayData, mesh: orm.List, U: orm.Float, J: orm.Float) -> dict:
    """Irreducible k-mesh, bands and the uniform-U, J vertex shared by the RPA/FLEX steps."""
    rvec, ham_r = _model(hoppings)
    Nx, Ny, Nz = _mesh(mesh)
    Norb = ham_r.shape[1]
    klist, kmap, invk = flibs.gen_irr_k_TRS(Nx, Ny, Nz)
    hamk = flibs.gen_ham(klist, ham_r, rvec)
    eig, uni = flibs.get_eig(hamk)
    olist, site = plibs.get_chi_orb_list(Norb, [Norb])
    Smat, Cmat = flibs.gen_SCmatrix(olist, site, U=U.value, J=J.value)
    return dict(rvec=rvec, ham_r=ham_r, Nx=Nx, Ny=Ny, Nz=Nz, Norb=Norb, klist=klist, kmap=kmap, invk=invk,
                hamk=hamk, eig=eig, uni=uni, olist=olist, Smat=Smat, Cmat=Cmat)


def _green(st: dict, mu: float, temp: float, Nw: int, sigma: orm.ArrayData | None) -> np.ndarray:
    """Bare G, or the FLEX one when sigma (output of flex_self_energy, with mu = its mu_self) is given."""
    if sigma is None:
        return flibs.gen_Green0(st['eig'], st['uni'], mu, temp, Nw)
    sigmak = sigma.get_array('sigma')
    if sigmak.shape[2] != Nw or sigmak.shape[3] != len(st['klist']):
        raise ValueError(f'sigma has shape {sigmak.shape}, not Nw = {Nw} on {len(st["klist"])} k-points')
    return flibs.gen_green(sigmak, st['hamk'], mu, temp)


# --------------------------------------------------------------------------- #
#  model
# --------------------------------------------------------------------------- #
@calcfunction
def load_hoppings(hop_file: orm.SinglefileData, ftype: orm.Int) -> orm.ArrayData:
    """Read a hopping file with plibs.import_hoppings (ftype as in main.py)."""
    # import_hoppings takes a path: copy the stored file out of the repository
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, hop_file.filename)
        with hop_file.open(mode='rb') as src, open(path, 'wb') as dst:
            dst.write(src.read())
        rvec, ham_r, _, _ = plibs.import_hoppings(path, ftype.value)
    return _arrays(rvec=np.asarray(rvec, dtype=np.float64), ham_r=np.asarray(ham_r, dtype=np.complex128))


@calcfunction
def chemical_potential(hoppings: orm.ArrayData, fill: orm.Float, temp: orm.Float, mesh: orm.List) -> orm.Float:
    """Chemical potential of the bare bands at filling fill on the full mesh (plibs.calc_mu)."""
    rvec, ham_r = _model(hoppings)
    Nx, Ny, Nz = _mesh(mesh)
    Nk, klist = plibs.gen_klist(Nx, Ny, Nz, sw_pp=False)
    eig = plibs.get_eigs(klist, ham_r, [], rvec)[0]
    return orm.Float(plibs.calc_mu(eig, Nk, fill.value, temp.value))


# --------------------------------------------------------------------------- #
#  RPA / FLEX
# --------------------------------------------------------------------------- #
@calcfunction
def flex_self_energy(hoppings: orm.ArrayData, mu: orm.Float, fill: orm.Float, temp: orm.Float,
                     U: orm.Float, J: orm.Float, mesh: orm.List, Nw: orm.Int) -> dict:
    """FLEX self-energy on the uniform Matsubara grid (flibs.mkself)."""
    st = _normal_state(hoppings, mesh, U, J)
    sigmak, mu_self = flibs.mkself(st['Smat'], st['Cmat'], st['kmap'], st['invk'], st['olist'], st['hamk'],
                                   st['eig'], st['uni'], mu.value, fill.value, temp.value, Nw.value,
                                   st['Nx'], st['Ny'], st['Nz'], False, False)
    sigma = _arrays(sigma=sigmak, klist=st['klist'], iw=plibs.fermion_frequencies(temp.value, Nw.value))
    return {'sigma': sigma, 'mu_self': orm.Float(mu_self)}


@calcfunction
def rpa_susceptibility(hoppings: orm.ArrayData, mu: orm.Float, temp: orm.Float, U: orm.Float, J: orm.Float,
                       mesh: orm.List, Nw: orm.Int, sigma: orm.ArrayData | None = None) -> dict:
    """
    RPA spin / charge susceptibility at v = 0 (leading eigenvalue of the orbital matrix) on the
    irreducible q-mesh; with sigma (and mu = mu_self) the FLEX one.  Also the Stoner factor.
    """
    st = _normal_state(hoppings, mesh, U, J)
    Gk = _green(st, mu.value, temp.value, Nw.value, sigma)
    chi0, stoner = flibs.get_chi0(st['Smat'], st['Cmat'], Gk, st['olist'], st['kmap'], st['invk'], temp.value,
                                  st['Nx'], st['Ny'], st['Nz'])
    chis, chic = flibs.get_chis_chic(chi0, st['Smat'], st['Cmat'])
    chi = _arrays(qlist=st['klist'], chis=chis, chic=chic,
                  chis_q=flibs.get_eig_or_tr_chi(chis, st['invk'], True),   # full mesh, order of kmap
                  chic_q=flibs.get_eig_or_tr_chi(chic, st['invk'], True), kmap=st['kmap'])
    return {'chi': chi, 'stoner': orm.Float(stoner)}


@calcfunction
def linearized_eliashberg(hoppings: orm.ArrayData, mu: orm.Float, temp: orm.Float, U: orm.Float, J: orm.Float,
                          mesh: orm.List, Nw: orm.Int, gap_sym: orm.Int, sigma: orm.ArrayData | None = None) -> dict:
    """
    Linearized Eliashberg equation with the RPA (or, with sigma, FLEX) pairing interaction:
    the leading eigenvalue lambda_eliash, the gap Delta(k, iw_n) and its Wannier-R export
    (output_gap_wannier), which eilenberger_bulk takes as the pairing form factor.
    """
    st = _normal_state(hoppings, mesh, U, J)
    Gk = _green(st, mu.value, temp.value, Nw.value, sigma)
    chi0, stoner = flibs.get_chi0(st['Smat'], st['Cmat'], Gk, st['olist'], st['kmap'], st['invk'], temp.value,
                                  st['Nx'], st['Ny'], st['Nz'])
    plist = flibs.get_plist(st['rvec'], st['ham_r'])
    init_delta = plibs.get_initial_gap(st['klist'], st['Norb'], gap_sym.value)
    gap, lam = flibs.linearized_eliashberg(chi0, Gk, st['uni'], init_delta, st['Smat'], st['Cmat'], st['olist'],
                                           plist, st['kmap'], st['invk'], st['Nx'], st['Ny'], st['Nz'], temp.value,
                                           gap_sym.value)
    with _workdir():
        plibs.output_gap_wannier(gap, st['kmap'], st['invk'], st['Nx'], st['Ny'], st['Nz'], Nw.value, temp.value)
        gap_wannier = _npz_to_arraydata('gap_wannier.npz')
    return {'lambda_eliash': orm.Float(lam), 'stoner': orm.Float(stoner), 'gap_wannier': gap_wannier,
            'gap': _arrays(gap=gap, klist=st['klist'], iw=plibs.fermion_frequencies(temp.value, Nw.value))}


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
    rvec, ham_r = _model(hoppings)
    Nx, Ny, Nz = _mesh(mesh)
    _, klist, eig, uni, _ = plibs.get_emesh(Nx, Ny, Nz, ham_r, [], rvec, np.eye(3), sw_uni=True)
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
    return {'result': orm.Dict(result),
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
def flex_workflow(hop_file: orm.SinglefileData, ftype: orm.Int, fill: orm.Float, temp: orm.Float,
                  U: orm.Float, J: orm.Float, mesh: orm.List, Nw: orm.Int) -> dict:
    """hopping file -> bare mu -> FLEX self-energy, recorded as one provenance graph."""
    hoppings = load_hoppings(hop_file, ftype)
    mu = chemical_potential(hoppings, fill, temp, mesh)
    flex = flex_self_energy(hoppings, mu, fill, temp, U, J, mesh, Nw)
    return {'hoppings': hoppings, 'mu': mu, 'sigma': flex['sigma'], 'mu_self': flex['mu_self']}


@workfunction
def eliashberg_workflow(hop_file: orm.SinglefileData, ftype: orm.Int, fill: orm.Float, temp: orm.Float,
                        U: orm.Float, J: orm.Float, mesh: orm.List, Nw: orm.Int, gap_sym: orm.Int,
                        use_flex: orm.Bool = lambda: orm.Bool(False)) -> dict:
    """hopping file -> mu -> [FLEX Sigma] -> chi_s/chi_c and the linearized Eliashberg gap (RPA or FLEX)."""
    hoppings = load_hoppings(hop_file, ftype)
    mu = chemical_potential(hoppings, fill, temp, mesh)
    out = {'hoppings': hoppings, 'mu': mu}
    extra = {}
    if use_flex.value:
        flex = flex_self_energy(hoppings, mu, fill, temp, U, J, mesh, Nw)
        out.update(sigma=flex['sigma'], mu_self=flex['mu_self'])
        mu, extra = flex['mu_self'], {'sigma': flex['sigma']}
    chi = rpa_susceptibility(hoppings, mu, temp, U, J, mesh, Nw, **extra)
    el = linearized_eliashberg(hoppings, mu, temp, U, J, mesh, Nw, gap_sym, **extra)
    out.update(chi=chi['chi'], stoner=chi['stoner'], lambda_eliash=el['lambda_eliash'], gap=el['gap'],
               gap_wannier=el['gap_wannier'])
    return out


@workfunction
def eliashberg_eilenberger_workflow(hop_file: orm.SinglefileData, ftype: orm.Int, fill: orm.Float,
                                    temp: orm.Float, U: orm.Float, J: orm.Float, mesh: orm.List, Nw: orm.Int,
                                    gap_sym: orm.Int, wc: orm.Float, coupling: orm.Float, eil_mesh: orm.List,
                                    use_flex: orm.Bool = lambda: orm.Bool(False)) -> dict:
    """
    eliashberg_workflow, then the homogeneous Eilenberger gap with the RPA/FLEX gap as the pairing
    form factor (on the bare-band Fermi surface at the bare mu, on eil_mesh).
    """
    el = eliashberg_workflow(hop_file, ftype, fill, temp, U, J, mesh, Nw, gap_sym, use_flex)
    eil = eilenberger_bulk(el['hoppings'], el['mu'], temp, wc, coupling, gap_sym, eil_mesh,
                           gap_wannier=el['gap_wannier'])
    return {**el, 'eilenberger': eil['result'], 'f_avg': eil['f_avg']}
