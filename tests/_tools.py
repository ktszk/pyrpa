#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Shared utilities for the pyrpa test-suite.

Import from a test file living in tests/ as

    import _tools as T

(the directory of the test file is on sys.path both under pytest and when the
file is executed directly).  Contents:

  * output control ........ silence()
  * distributions ......... fermi, bose, fermionic_matsubara, bosonic_matsubara
  * model builders ........ one_orbital_square, two_orbital_square,
                            band_diagonal_model, cylindrical_fs,
                            full_grid_dispersion
  * exact references ...... two_pole_poles_weights, two_pole_green,
                            exact_chi0_1orb  (Lindhard bubble, optionally with
                            an exactly solvable two-pole self-energy)
  * numerics .............. rel_err, fit_power_order
  * physics assertions .... assert_hermitian, assert_positive_semidefinite,
                            assert_green_causal, assert_raises
  * runner ................ run_standalone (shared "no pytest needed" driver)

Requires the Fortran library libfmod.so (cd libs && make FC=ifx SL=MKL).
"""
import contextlib
import inspect
import io
import os
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import libs.flibs as F
import libs.plibs as P


# --------------------------------------------------------------------------- #
#  output control
# --------------------------------------------------------------------------- #
@contextlib.contextmanager
def silence():
    """Context manager that swallows stdout of the wrapped block at the
    file-descriptor level, so progress lines printed by the Fortran library
    (which bypass sys.stdout) are silenced too.  Python-level prints inside
    the block are also swallowed."""
    sys.stdout.flush()
    saved = os.dup(1)
    devnull = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull, 1)
        with contextlib.redirect_stdout(io.StringIO()):
            yield
    finally:
        sys.stdout.flush()
        os.dup2(saved, 1)
        os.close(saved)
        os.close(devnull)


# --------------------------------------------------------------------------- #
#  distributions / Matsubara grids
# --------------------------------------------------------------------------- #
def fermi(e, temp, mu=0.0):
    """Overflow-safe Fermi function f((e-mu)/T), element-wise."""
    x = np.atleast_1d(np.asarray(e - mu, dtype=np.float64)) / temp
    out = np.empty_like(x)
    pos = x > 0
    out[pos] = np.exp(-x[pos]) / (1.0 + np.exp(-x[pos]))
    out[~pos] = 1.0 / (1.0 + np.exp(x[~pos]))
    return out.reshape(np.shape(e - mu)) if np.ndim(e - mu) else out[0]


def bose(w, temp):
    """Overflow-safe Bose function n_B(w/T) for w != 0, element-wise."""
    x = np.asarray(w, dtype=np.float64) / temp
    ax = np.abs(x)
    n = np.exp(-ax) / (1.0 - np.exp(-ax))
    return np.where(x > 0, n, -(1.0 + n))


def fermionic_matsubara(temp, Nw):
    """Positive fermionic frequencies w_n = (2n+1) pi T, n = 0..Nw-1
    (the code stores only these; negative ones follow from TRS)."""
    return (2.0 * np.arange(Nw) + 1.0) * np.pi * temp


def bosonic_matsubara(temp, Nm):
    """Bosonic frequencies nu_m = 2 m pi T, m = 0..Nm-1 (code index j = m+1)."""
    return 2.0 * np.arange(Nm) * np.pi * temp


# --------------------------------------------------------------------------- #
#  model builders
# --------------------------------------------------------------------------- #
def full_grid_dispersion(ham_r, rvec, Nx, Ny, Nz=1):
    """Band energies on the full (unreduced) k-grid, shape [Nx,Ny,Nz,Norb].
    Index [ix,iy,iz] corresponds to k = (ix/Nx, iy/Ny, iz/Nz)."""
    ix, iy, iz = np.meshgrid(np.arange(Nx), np.arange(Ny), np.arange(Nz),
                             indexing='ij')
    kf = np.stack([ix.ravel() / Nx, iy.ravel() / Ny, iz.ravel() / Nz], axis=1)
    hk = F.gen_ham(kf, ham_r, rvec)
    ek = np.linalg.eigvalsh(hk)
    return ek.reshape(Nx, Ny, Nz, hk.shape[1])


def one_orbital_square(Nx=8, Ny=8, Nz=1, Nw=64, temp=0.1, mu=0.0,
                       tx=1.0, ty=0.5, U=0.0, J=0.0):
    """One-orbital square-lattice tight-binding bundle for RPA/FLEX tests.

    Nearest-neighbour hoppings -tx (x) and -ty (y); ty != tx keeps the C4
    breaking that exercises the k-map symmetry code.  Returns a dict with
    everything the flibs chi0/FLEX entry points need, plus eps_full
    [Nx,Ny,Nz] for exact references (exact_chi0_1orb)."""
    rvec = np.array(
        [[0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]],
        dtype=np.float64,
    )
    ham_r = np.array(
        [[[0.0]], [[-tx]], [[-tx]], [[-ty]], [[-ty]]],
        dtype=np.complex128,
    )
    klist, kmap, invk = F.gen_irr_k_TRS(Nx, Ny, Nz)
    hamk = F.gen_ham(klist, ham_r, rvec)
    eig, uni = F.get_eig(hamk)
    Gk = F.gen_Green0(eig, uni, mu=mu, temp=temp, Nw=Nw)
    olist, site = P.get_chi_orb_list(1, [1])
    Smat, Cmat = F.gen_SCmatrix(olist, site, U=U, J=J)
    eps_full = full_grid_dispersion(ham_r, rvec, Nx, Ny, Nz)[..., 0]
    plist = np.array([1.0], dtype=np.float64)
    return dict(
        Nx=Nx, Ny=Ny, Nz=Nz, Nw=Nw, temp=temp, mu=mu, rvec=rvec, ham_r=ham_r,
        klist=klist, kmap=kmap, invk=invk, hamk=hamk, eig=eig, uni=uni,
        Gk=Gk, olist=olist, site=site, Smat=Smat, Cmat=Cmat, plist=plist,
        eps_full=eps_full,
    )


def two_orbital_square(Nx=8, Ny=8, Nz=1, Nw=64, temp=0.1, mu=0.0,
                       t1=1.0, t2=0.5, de=0.4, V=0.2, U=0.0, J=0.0):
    """Two-orbital square-lattice bundle: diagonal hoppings -t1/-t2, crystal
    field de on orbital 2 and a constant on-site hybridization V.  The bands
    are exactly e_pm(k) = (e1+e2)/2 +- sqrt(((e1-e2)/2)^2 + V^2) with
    e1 = -2 t1 (cx+cy), e2 = -2 t2 (cx+cy) + de, so eigen-decomposition
    dependent code can be checked against a closed form."""
    rvec = np.array(
        [[0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]],
        dtype=np.float64,
    )
    onsite = np.array([[0.0, V], [V, de]])
    hop = lambda t1_, t2_: np.diag([-t1_, -t2_]).astype(np.complex128)
    ham_r = np.array(
        [onsite, hop(t1, t2), hop(t1, t2), hop(t1, t2), hop(t1, t2)],
        dtype=np.complex128,
    )
    klist, kmap, invk = F.gen_irr_k_TRS(Nx, Ny, Nz)
    hamk = F.gen_ham(klist, ham_r, rvec)
    eig, uni = F.get_eig(hamk)
    Gk = F.gen_Green0(eig, uni, mu=mu, temp=temp, Nw=Nw)
    olist, site = P.get_chi_orb_list(2, [2])
    Smat, Cmat = F.gen_SCmatrix(olist, site, U=U, J=J)
    plist = np.array([1.0, 1.0], dtype=np.float64)
    return dict(
        Nx=Nx, Ny=Ny, Nz=Nz, Nw=Nw, temp=temp, mu=mu, rvec=rvec, ham_r=ham_r,
        t1=t1, t2=t2, de=de, V=V,
        klist=klist, kmap=kmap, invk=invk, hamk=hamk, eig=eig, uni=uni,
        Gk=Gk, olist=olist, site=site, Smat=Smat, Cmat=Cmat, plist=plist,
    )


def band_diagonal_model(Nk, Norb, seed=0, Ewidth=0.6):
    """Random band-diagonal model for transport tests: sorted eigenvalues,
    the Kubo velocity vk [Nk,Norb,Norb,3] (purely band-diagonal, so no
    interband channel) and the Boltzmann velocity veloc [Nk,Norb,3] holding
    the same diagonal entries."""
    rng = np.random.default_rng(seed)
    eig = np.sort(rng.uniform(-Ewidth, Ewidth, size=(Nk, Norb)),
                  axis=1).astype(np.float64)
    vdiag = rng.standard_normal((Nk, Norb, 3))
    vk = np.zeros((Nk, Norb, Norb, 3), dtype=np.complex128)
    for n in range(Norb):
        vk[:, n, n, :] = vdiag[:, n, :]
    veloc = np.ascontiguousarray(vdiag, dtype=np.float64)
    return eig, vk, veloc


def cylindrical_fs(Nb=120, gap_sym='d'):
    """Synthetic cylindrical Fermi surface for the Eilenberger solvers:
    uniform weights and the FS-normalized form factor <phi^2> = 1.
    gap_sym in {'s', 'd', 'dxy'}."""
    beta = np.linspace(0.0, 2.0 * np.pi, Nb, endpoint=False)
    wf = np.ones(Nb)
    if gap_sym == 's':
        phi = np.ones(Nb)
    elif gap_sym == 'd':
        phi = np.sqrt(2.0) * np.cos(2.0 * beta)
    elif gap_sym == 'dxy':
        phi = np.sqrt(2.0) * np.sin(2.0 * beta)
    else:
        raise ValueError(gap_sym)
    return wf, phi


# --------------------------------------------------------------------------- #
#  exact references (single orbital)
# --------------------------------------------------------------------------- #
def two_pole_poles_weights(eps, a, delta):
    """Pole decomposition of G = 1/(x - eps - a^2/(x - delta)), x = iw + mu:
    returns E [2,...] and weights W [2,...] with G = sum_p W_p/(x - E_p).
    a = 0 reduces to the free propagator (both half-weight poles at eps)."""
    eps = np.asarray(eps, dtype=np.float64)
    if a == 0.0:
        E = np.stack([eps, eps])
        W = np.full_like(E, 0.5)
        return E, W
    tr = eps + delta
    det = eps * delta - a * a
    disc = np.sqrt(np.maximum(tr * tr - 4.0 * det, 1e-300))
    Ep = 0.5 * (tr + disc)
    Em = 0.5 * (tr - disc)
    wp = (Ep - delta) / (Ep - Em)
    wm = (Em - delta) / (Em - Ep)
    return np.stack([Ep, Em]), np.stack([wp, wm])


def two_pole_green(eig, temp, mu, Nw, a=0.0, delta=0.0):
    """Interacting single-orbital Green's function with the exactly solvable
    self-energy Sigma = a^2/(iw + mu - delta), in the code layout
    Gk [1,1,Nw,Nk] (positive fermionic frequencies only).  a = 0 gives the
    non-interacting G0."""
    wl = fermionic_matsubara(temp, Nw)
    x = mu + 1j * wl
    epsk = np.asarray(eig)[:, 0]
    G = 1.0 / (x[None, :] - epsk[:, None] - (a * a / (x[None, :] - delta) if a != 0.0 else 0.0))
    return np.ascontiguousarray(G.T[None, None, :, :])


def exact_chi0_1orb(eps_full, qlist, mlist, temp, mu, a=0.0, delta=0.0):
    """Exact Matsubara bubble chi0(q, i nu_m) of a single-orbital model, in
    the code sign convention chi = -T sum_w G(k+q, iw+inu) G(k, iw) > 0.

    eps_full : [Nx,Ny,Nz] dispersion on the full grid (full_grid_dispersion)
    qlist    : fractional q-vectors [Nq,3] (e.g. the irreducible klist)
    mlist    : bosonic indices m (nu_m = 2 pi m T); code index j = m+1
    a, delta : optional two-pole self-energy Sigma = a^2/(iw + mu - delta)

    Returns chi_ex [Nq, len(mlist)] complex.  The exact frequency sum is the
    Lindhard form sum_{p,p'} W_p(k+q) W_p'(k) (f_p - f_p')/(i nu + E_p' - E_p)
    with the degenerate m=0 limit f(1-f)/T, so this is free of any Matsubara
    cutoff and serves as the reference for chi0 tail-convergence tests."""
    eps_full = np.asarray(eps_full, dtype=np.float64)
    Nx, Ny, Nz = eps_full.shape
    E, W = two_pole_poles_weights(eps_full, a, delta)      # [2,Nx,Ny,Nz]
    out = np.zeros((len(qlist), len(mlist)), dtype=complex)
    for iq, q in enumerate(qlist):
        sh = tuple(-int(round(q[d] * n)) % n for d, n in enumerate((Nx, Ny, Nz)))
        Eq = np.roll(E, shift=sh, axis=(1, 2, 3))
        Wq = np.roll(W, shift=sh, axis=(1, 2, 3))
        for jm, m in enumerate(mlist):
            nu = 2.0 * np.pi * m * temp
            tot = 0.0 + 0.0j
            for p in range(2):
                for pp in range(2):
                    aE = Eq[p]                              # E_p(k+q)
                    bE = E[pp]                              # E_p'(k)
                    w2 = Wq[p] * W[pp]
                    fa = fermi(aE, temp, mu)
                    fb = fermi(bE, temp, mu)
                    de = 1j * nu + bE - aE
                    deg = (m == 0) & (np.abs(aE - bE) < 1e-9)
                    safe = np.where(np.abs(de) < 1e-30, 1.0, de)
                    val = np.where(deg, fb * (1.0 - fb) / temp,
                                   (fa - fb) / safe)
                    tot += (w2 * val).sum()
            out[iq, jm] = tot / (Nx * Ny * Nz)
    return out


# --------------------------------------------------------------------------- #
#  numerics
# --------------------------------------------------------------------------- #
def rel_err(a, b):
    """max|a-b| / max|b| over the whole arrays."""
    return np.abs(np.asarray(a) - np.asarray(b)).max() / np.abs(b).max()


def fit_power_order(ns, errs):
    """Least-squares power-law order p of err ~ C/n^p from paired samples
    (e.g. Nw values and truncation errors).  Use to lock in convergence
    orders: sharp Matsubara cutoff ~1, tail-corrected chi0 ~2."""
    ns = np.log(np.asarray(ns, dtype=np.float64))
    errs = np.log(np.asarray(errs, dtype=np.float64))
    return -np.polyfit(ns, errs, 1)[0]


# --------------------------------------------------------------------------- #
#  physics assertions
# --------------------------------------------------------------------------- #
def assert_hermitian(m, tol=1e-10, msg='matrix'):
    """Assert m[..., i, j] == conj(m[..., j, i]) within tol (absolute)."""
    m = np.asarray(m)
    dev = np.abs(m - np.conj(np.swapaxes(m, -1, -2))).max()
    assert dev < tol, f"{msg} is not Hermitian: max deviation {dev:.3e} >= {tol:.1e}"


def assert_positive_semidefinite(m, tol=1e-10, msg='matrix'):
    """Assert every [..., :, :] block is Hermitian positive semi-definite
    (eigenvalues > -tol).  Physical static susceptibilities chi(q, nu=0) of a
    stable normal state must pass this."""
    assert_hermitian(m, tol=max(tol, 1e-8), msg=msg)
    ev = np.linalg.eigvalsh(np.asarray(m))
    assert ev.min() > -tol, \
        f"{msg} is not positive semi-definite: min eigenvalue {ev.min():.3e}"
    return ev


def assert_green_causal(Gk, tol=1e-12, msg='Green function'):
    """Assert Im G_aa(k, iw_n) <= tol for all stored (positive) Matsubara
    frequencies -- the causality/Herglotz property of any physical G.
    Gk layout: [Norb, Norb, Nw, Nk]."""
    Gk = np.asarray(Gk)
    diag = np.einsum('aawk->awk', Gk)
    worst = diag.imag.max()
    assert worst <= tol, \
        f"{msg} violates causality: max Im G_diag = {worst:.3e} > {tol:.1e}"


def assert_raises(exc, fn, *args, **kwargs):
    """Minimal pytest.raises stand-in that also works in the standalone
    runner (no pytest import).  Returns the caught exception."""
    try:
        fn(*args, **kwargs)
    except exc as e:
        return e
    raise AssertionError(f"{fn} did not raise {exc.__name__}")


# --------------------------------------------------------------------------- #
#  standalone runner
# --------------------------------------------------------------------------- #
def run_standalone(namespace, argv=None):
    """Shared 'no pytest required' driver: run every test_* callable in
    namespace (a module's globals()), report PASS/FAIL with timings and
    return a shell exit code.  Supplies a temporary directory to tests with
    a tmp_path parameter.  Optional argv (default sys.argv[1:]): substring
    filters, only tests whose name contains one of them are run.

    Usage at the end of a test module:

        if __name__ == '__main__':
            import _tools
            sys.exit(_tools.run_standalone(globals()))
    """
    if argv is None:
        argv = sys.argv[1:]
    tests = [v for k, v in sorted(namespace.items())
             if k.startswith('test_') and callable(v)]
    if argv:
        tests = [t for t in tests if any(s in t.__name__ for s in argv)]
    width = max((len(t.__name__) for t in tests), default=10) + 2
    npass = 0
    for t in tests:
        t0 = time.time()
        try:
            kwargs = {}
            with contextlib.ExitStack() as st:
                if 'tmp_path' in inspect.signature(t).parameters:
                    td = st.enter_context(tempfile.TemporaryDirectory())
                    kwargs['tmp_path'] = Path(td)
                t(**kwargs)
            print(f"  PASS  {t.__name__:{width}s} ({time.time()-t0:5.1f}s)",
                  flush=True)
            npass += 1
        except Exception as e:
            print(f"  FAIL  {t.__name__:{width}s} ({time.time()-t0:5.1f}s)"
                  f"  -> {type(e).__name__}: {e}", flush=True)
    print(f"\n{npass}/{len(tests)} passed")
    return 0 if npass == len(tests) else 1
