# pyrpa test-suite

Prerequisite: the Fortran library is built (`cd libs && make`). Japanese version: [README_JP.md](README_JP.md).

## How to run

```sh
# recommended: pytest (from the repository root; pytest.ini picks up tests/)
pytest                        # all tests
pytest -k chi0                # filter by name
pytest --runslow              # also the tests marked slow (_tools.slow)

# without pytest (clusters etc.): every file runs standalone
python tests/test_rpa_flex.py           # one file
python tests/test_rpa_flex.py chi0      # only the tests whose name contains 'chi0'
python tests/test_eilenberger.py --runslow   # including the slow ones
python tests/run_all.py                 # all files + a summary table
python tests/run_all.py rpa trans -q    # file-name filter + quiet
python tests/run_all.py --runslow       # everything, slow tests included
```

## Files

| File | Covers |
|---|---|
| `test_bands.py` | Bands and Fermi-surface geometry: reading Wannier90 `_hr.dat` and the phase of `gen_ham`, the first Brillouin zone of the 3D Fermi surface, group velocity and inverse effective mass (orthogonal and MLO bases), cyclotron and dHvA extremal orbits |
| `test_transport.py` | Transport coefficients (Boltzmann / Kubo, Wiedemann-Franz, k-point weights, Hall response, relaxation-time models) |
| `test_optical.py` | Optical constants, reflectivity-to-colour conversion (CIE XYZ, D65, sRGB), the colour of non-metals (transmission, Kubelka-Munk) |
| `test_impurity.py` | Impurity supercell and the unfolded spectrum (phase sign, sum rule), CPA (also in a non-orthogonal basis) |
| `test_rpa_flex.py` | RPA / FLEX / Eliashberg (χ0, vertices, tail correction, mkself, mkself_soc and the double-counting cancellation), Green function, FLEX and Eliashberg in the non-orthogonal (MLO) basis, phonon-mediated pairing from EPA |
| `test_eliash_nosym.py` | Symmetry-free linearized Eliashberg (Zeeman, FFLO), consistency with the SOC version (bare term -V/2, triplet seed and parity projection, band-basis output), exact-kernel check of multi-orbital triplets, branch tracking under SOC |
| `test_freq_grid.py` | IR-basis frequency grid: FLEX, linearized / nonlinear Eliashberg, SOC and `sw_nosym` agree with the Nw→∞ limit of the uniform grid (SKIP without sparse-ir) |
| `test_nmr.py` | SC-state dynamic susceptibility (including the non-orthogonal BdG) and NMR (Knight shift, 1/T1T) |
| `test_eilenberger.py` | Quasiclassical Eilenberger (homogeneous, surface, vortex, vortex lattice) and the Fortran grid / chord kernels of the vortex solvers |
| `test_main.py` | Chemical potential of the Eilenberger modes in `main.py` |
| `test_aiida.py` | AiiDA sample (examples/aiida): the RPA/FLEX (μ, Σ, χ, linearized Eliashberg) and Eilenberger (bulk, surface, vortex) calcfunctions agree with direct plibs/flibs calls, and the provenance of the Eliashberg gap fed to Eilenberger is recorded (SKIP without aiida-core) |
| `test_tools.py` | Self-checks of the shared tools `_tools.py` (doubling as usage examples), CPU-based library selection, runtime settings of the gfortran build |
| `_tools.py` | Shared utilities (below) |
| `conftest.py` | pytest configuration (sys.path, `slow` marker, `silence` fixture) |
| `run_all.py` | Batch runner that needs no pytest |

## Shared tools `_tools.py`

Used from the test files (directly in tests/) as `import _tools as T`.

- **Model builders** — return, as a dict, the set handed to the Fortran entry points (klist/kmap/invk/hamk/eig/uni/Gk/olist/Smat/Cmat...)
  - `one_orbital_square(Nx, Ny, Nz, Nw, temp, mu, tx, ty, U, J)` — one-orbital square lattice; carries `eps_full` (dispersion on the full grid) for the exact references
  - `two_orbital_square(..., t1, t2, de, V)` — two orbitals + hybridization; bands checkable against the closed form `(e1+e2)/2 ± sqrt(((e1-e2)/2)² + V²)`
  - `band_diagonal_model(Nk, Norb, seed)` — random band-diagonal model for the transport tests
  - `cylindrical_fs(Nb, gap_sym)` — cylindrical FS for Eilenberger (normalized to `<φ²>=1`)
- **Exact references** — analytic results free of the Matsubara cutoff error
  - `exact_chi0_1orb(eps_full, qlist, mlist, temp, mu, a, delta)` — Lindhard bubble (the code's sign convention, χ>0); with `a≠0` also with the exactly solvable two-pole self-energy Σ = a²/(iω+μ−δ) → the reference of the Nw convergence-order tests of χ0
  - `two_pole_green(eig, temp, mu, Nw, a, delta)` — the corresponding G (code layout [1,1,Nw,Nk])
- **Numerical helpers**
  - `fit_power_order(ns, errs)` — log fit of p in err ∼ C/n^p (pins sharp cutoff ≈1, tail correction ≈2)
  - `rel_err(a, b)`
- **Physical assertions** — one-line checks of properties the physics must satisfy
  - `assert_hermitian(m)` / `assert_positive_semidefinite(m)` — e.g. the static χs/χc matrices
  - `assert_green_causal(Gk)` — Im G_aa(iω_n>0) ≤ 0 (causality)
  - `assert_raises(exc, fn, ...)` — a pytest.raises replacement that works without pytest
- **Other**
  - `fermi / bose / fermionic_matsubara / bosonic_matsubara` — overflow-safe distribution functions
  - `silence()` — context manager that suppresses the progress output of the Fortran wrappers
  - `run_standalone(globals())` — the shared `__main__` runner at the end of every file (tmp_path support, name filter from the arguments, `--runslow`)
  - `slow` — decorator for long tests (below)

## Template for a new test

```python
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _tools as T
import libs.flibs as F

def test_something():
    md = T.one_orbital_square(Nx=8, Ny=8, Nw=64, temp=0.05, mu=0.3)
    with T.silence():
        chi, stoner = F.get_chi0(md['Smat'], md['Cmat'], md['Gk'], md['olist'],
                                 md['kmap'], md['invk'], md['temp'],
                                 md['Nx'], md['Ny'], md['Nz'])
    ref = T.exact_chi0_1orb(md['eps_full'], md['klist'], [0], md['temp'], md['mu'])
    assert np.abs(chi[0, 0, 0, :] - ref[:, 0]).max() < ...

if __name__ == '__main__':
    sys.exit(T.run_standalone(globals()))
```

Mark long tests with `@T.slow` (`_tools.slow`). Under pytest, standalone and `run_all.py` alike they are
skipped by default and run only with `--runslow` (this works without pytest too).
At present the four marked tests of `test_eilenberger.py`; they alone take about 70% of that file's run
time, so the default run takes roughly half or less (measured: all 57 tests 854 s → the default 53 in 394 s).
