# pyrpa
pyrpa is a Python + Fortran code for model calculations on tight-binding Hamiltonians
(Wannier90, ecalj and other hopping files), from band structures to transport, magnetic
fluctuations and superconductivity. The heavy kernels are OpenMP-parallel Fortran in
`libs/libfmod.so`; everything is driven from the settings at the top of `main.py`.

## Features
- **Electronic structure**: bands, DOS, 2D/3D Fermi surfaces, spectral function $A(k,\omega)$
  (optionally with the FLEX self-energy), band filling, effective and cyclotron masses,
  dHvA frequency vs. field angle
- **Transport and optics**: Boltzmann and Kubo conductivities (electrical, thermal, Seebeck),
  weak-field Hall coefficient with band decomposition, relaxation-time models (constant, EPA,
  FLEX quasiparticles), optical conductivity, reflectivity and material colour
- **Fluctuations and superconductivity**: RPA spin/charge/pairing susceptibilities, FLEX
  self-energy, linearized and nonlinear Eliashberg equations, the symmetry-free Eliashberg
  equation (Zeeman field, FFLO), phonon-mediated pairing from EPA, SC-state spin
  susceptibility and NMR (Knight shift, $1/T_1T$)
- **Quasiclassical theory**: Eilenberger/Riccati solvers for $T_c$, penetration depth, surface
  Andreev bound states, vortices and vortex lattices
- **Other**: spin-orbit coupling, non-orthogonal MLO basis (ecalj), impurity supercells and CPA,
  the IR basis for the Matsubara axis, an AiiDA interface with provenance

## Requirements
- Fortran compiler: ifx/ifort, flang/amdflang or gfortran (tested: ifx, AMD Flang, gfortran 13.3)
- BLAS/LAPACK: MKL, AOCL (BLIS and FLAME) or OpenBLAS
- FFT: MKL's FFTW3 interface (ifx/ifort), AOCL-FFTZ (flang; `FFTZ_ROOT`, default
  `~/opt/aocl-fftz`), or the bundled pocketfft otherwise. FFTW / AMD-FFTW (GPL) is not used.
- Make (or CMake)
- Python 3.11 or later (tested: 3.11, 3.13, 3.14) with
  - numpy, scipy, matplotlib, scikit-image
  - sparse-ir (the IR basis, the default Matsubara grid of FLEX and Eliashberg; without it set
    `freq_grid = 'matsubara'` in `main.py`)
  - ana_cont (optional; MaxEnt analytic continuation from the library functions)
  - aiida-core (optional; only for `examples/aiida/`)
- Typst (optional; to build the manuals)

## Quick start
```bash
cd libs && make && cd ..     # compiler, BLAS/LAPACK and CPU target are detected automatically
python main.py               # runs the model and calculation mode set at the top of main.py
```
Edit the header of `main.py` to choose the model and the calculation, e.g. the band structure of
the bundled one-orbital square-lattice model:
```python
fname, ftype, brav, sw_soc = 'inputs/square.hop', 1, 0, False
option = 0            # 0 = band structure; see the CalcMode list in main.py
abc = [1., 1., 1.]    # lattice constants [Angstrom]
fill = 0.45           # filling (electrons per spin without SOC)
tempK = 200           # temperature [K]
```

## Building the Fortran kernel
`make` in `libs/` detects everything; values given on the command line take precedence.

| Variable | Meaning | Detected as |
|----------|---------|-------------|
| `FC` | compiler | the first of `ifx`, `ifort`, `flang`, `amdflang`, `gfortran` on `PATH` |
| `SL` | BLAS/LAPACK | `MKL` if `MKLROOT` is set, `AOCL` if `AOCL_ROOT` is set, else `OpenBLAS` |
| `ARCH` | CPU target | `znver3` on AMD, `x86-64-v3` on Intel (portable AVX2); `znver4` / `x86-64-v4` for AVX-512 only |

```bash
make FC=gfortran SL=OpenBLAS   # choose compiler and library
make ARCH=znver4               # choose the CPU target
make clean                     # remove build artifacts
```

On a cluster whose compute nodes differ from the login node (a library built for the wrong
CPU dies with `Illegal instruction`):
```bash
make auto        # one libfmod_<target>.so per CPU type found in the scheduler (Slurm/PBS)
make dispatch    # an AVX2 and an AVX-512 library without probing the nodes
```
At import, `libs/flibs/_loader.py` loads the best library the node can run.
`cmake -S libs -B build && cmake --build build` gives the same result as `make auto`.

## Repository layout
```
pyrpa/
├── main.py              # settings (top of the file) and one run_* function per calculation mode
├── libs/
│   ├── plibs/           # calculations on the Python side (model input, symmetry, drivers)
│   ├── flibs/           # ctypes wrappers of the Fortran kernel
│   ├── plots.py         # matplotlib figures
│   ├── src/             # Fortran sources (chi, FLEX, Eliashberg, transport, Eilenberger, ...)
│   ├── makefile, CMakeLists.txt, detect_arch.sh
│   └── libfmod*.so      # built by make
├── inputs/              # small test models (square.hop, square_2orb.hop, ...)
├── UJ/                  # stored cRPA U/J parameter sets (TOML)
├── examples/aiida/      # AiiDA interface (calcfunctions / workflows)
├── tests/               # regression and physics tests
└── Doc/                 # manuals and developer notes (Typst sources)
```

## Documentation
- User manual: `Doc/manual_EN.typ` (English) and `Doc/manual_JP.typ` (Japanese). Build the PDF with
  ```bash
  typst compile Doc/manual_EN.typ
  ```
  It covers the input formats, every calculation mode (`option`), all parameters,
  typical workflows and troubleshooting.
- Developer notes: `Doc/internals_EN.typ` / `Doc/internals_JP.typ` -- code layout, the functions
  each mode calls, call graphs, a function reference, implementation details, accuracy checks
  and the change history. The function reference and the call graphs are generated from the
  sources; rerun `python Doc/gen_internals.py` after changing the code.

## Tests
Build the library first, then from the repository root:
```bash
python tests/run_all.py            # all files, with a summary table (pytest not needed)
python tests/run_all.py --runslow  # also the slow tests
pytest tests                       # if pytest is installed
```
See [tests/README.md](tests/README.md) for what each file covers.

## License
pyrpa is released under the MIT License; see [LICENSE](LICENSE).
The bundled pocketfft header (`libs/src/pocketfft_hdronly.h`) is BSD-3-Clause; see
[THIRD_PARTY_NOTICES](THIRD_PARTY_NOTICES).
The stored U/J parameter sets in [UJ/](UJ/) are taken from published cRPA calculations;
each file records its reference (`source`, `doi`), which should be cited when the set is used.
