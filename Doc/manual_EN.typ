// Document format settings
#set text(font: "Times New Roman", size: 11pt)
#set page(paper: "a4", margin: (x: 2cm, y: 2cm))
#set math.equation(numbering: numbering.with("(1)"), supplement: "Eq.")
#set heading(numbering: "1.")
// External packages
#import "@preview/physica:0.9.3": *

#align(center)[
  #text(size: 18pt, weight: "bold")[pyrpa Manual]
  #v(0.5em)
  #text(size: 11pt)[2026]
]

#v(1em)

= About pyrpa

pyrpa is a Python code for computing and visualizing various physical properties from a model Hamiltonian (tight-binding model). It is primarily intended for post-processing results from first-principles calculations based on Wannier functions.

The main features currently supported are:

- Calculation and visualization of band structure and Fermi surfaces (2D / 3D)
- Density of states (DOS) calculation
- Spectral function $A(bold(k), omega)$ calculation
- Various transport coefficients (electrical conductivity, thermal conductivity, Seebeck coefficient, etc.) via Boltzmann theory and linear response theory
- Weak-field Hall coefficient $R_H$, Hall carrier density and Hall mobility, with a band-resolved (two-fluid) decomposition
- Spin susceptibility $chi_s$ and pairing susceptibility $phi$ based on RPA
- Dynamic spin susceptibility $chi_s^"SC"$ in the superconducting state
- Self-energy calculation via FLEX (Fluctuation EXchange approximation), on uniform Matsubara frequencies or the IR basis
- Linearized Eliashberg equation for the superconducting eigenvalue and gap function
- Nonlinear Eliashberg equation for self-consistent SC gap functions
- Band filling, cyclotron mass, and dHvA frequency calculations
- Quasiclassical Eilenberger/Riccati solvers: $T_c$, penetration depth, surface Andreev bound states, vortices and vortex lattices
- NMR in the superconducting state: Knight shift $K(T)$ and $1\/(T_1 T)$ from the BdG spin susceptibility
- Spectral function with impurities and via CPA
- Linearized Eliashberg equation without symmetry assumptions (Zeeman field, FFLO pair momentum $bold(Q)$)
- RPA / FLEX / Eliashberg in a non-orthogonal MLO basis (ecalj)
- Provenance-tracked runs through AiiDA (`examples/aiida`)

All settings in pyrpa are controlled by modifying variables at the top of `main.py`, above the library import lines.

= Notes on This Revision

The following changes alter results obtained with earlier versions; recompute the affected calculations (the notes of earlier revisions are in the git history).

- *FLEX and Eliashberg in a non-orthogonal MLO basis (ftype=3,4).* With a self-energy, $G^(-1)$ lacked the overlap $S$ ($G^(-1) = i omega + mu - H - Sigma$) and did not match the particle-number correction. For the Cu MLO model the FLEX $mu$ moved by about $+2.9$ eV from the bare value and the Stoner factor came out negative. Now $G^(-1) = (i omega + mu) S - H - Sigma$ and the particle number is $"tr"(S G)$. *Earlier MLO FLEX results, Eliashberg on top of them and self-energy-corrected spectra are all wrong.* RPA ($Sigma = 0$) built $G_0$ from the eigenvectors and is unchanged.
- *Normalized MLO orbitals.* ecalj's MLOs do not have unit norm (Cu: $S_(a a)(0) = 0.016$ for s, $0.565$ for d), so a U on the raw orbitals acted with a strength of about $U\/S_(a a)(0)^2$ (Cu RPA at U = 1: Stoner factor 85). The orbitals are now normalized on reading (the bands are unchanged). U, J, orbital weights and orbital-basis values of $Sigma$ and $Delta$ change, and files saved before do not follow the new convention. Interband velocity matrix elements (through the Löwdin orthogonalization) change by a few percent, and so do the optical conductivity and reflectivity (Cu: $Delta R lt.eq 0.02$). Group velocities and everything built only on them (Boltzmann transport, effective masses, dHvA) are unchanged.
- *`sw_chi0_tail` in an MLO basis.* The $tau = 0$ jump of $G_0$ was taken as the unit matrix, so in a non-orthogonal basis the truncation error stayed $O(1\/N_w)$ instead of $O(1\/N_w^2)$. The jump $S^(-1)(bold(k))$ is now convolved with the density matrix over $bold(k)$.
- *Nonlinear Eliashberg equation (option=23).* The normal Green's function of the SC state is now solved from the Nambu Dyson equation $G = [G_0^(-1) + Delta G_0^dagger Delta^dagger]^(-1)$ instead of the sub-iteration $G = G_0 (1 + Delta F^dagger)$ with the previous $F$ (which diverged once $Delta$ exceeded about $pi T$). With `sw_self=True` the static part of $Sigma$ is subtracted every iteration (`sub_sigma`; $mu$ is fixed in this loop) and $Sigma$ is mixed with the rate `pp`. A non-finite $Delta$ now stops with an error instead of being reported as converged. Converged gaps and $Sigma$ change.
- *Logical arguments of the Fortran kernels.* Some arguments passed from Python as a one-byte `c_bool` were received as a default four-byte `logical` (the spin channel `sw_spsym` of the SC-state $chi$, the eigenvalue/trace switch `sw_eig` of $chi$, ...). The other three bytes are undefined, so depending on the run they may have been misread; recompute option=12, 13, 27 results if in doubt.

New features (earlier results unchanged): the IR-basis frequency grid (`freq_grid='ir'`, see option=14), the symmetry-free linearized Eliashberg equation (`sw_nosym`, see option=15), the AiiDA interface, and the automatic setting for gfortran builds (see "Building the Fortran kernel").

= Requirements

The following Python packages are required:

- `numpy`
- `scipy`
- `matplotlib`
- `skimage` (used for Fermi surface plots option=2,3, cyclotron mass option=18, and dHvA orbits option=19)
- `sparse-ir` (optional; the IR basis, `freq_grid='ir'`; 1.x, 2.x, tested with 2.1.6)
- `ana_cont` (optional; MaxEnt analytic continuation, `sw_maxent=True`)
- `aiida-core` (optional; only for the AiiDA interface in `examples/aiida/`, see "AiiDA Interface")

Python 3.11 or later is required (the standard `tomllib` reads the stored U, J parameter sets). The Fortran kernel `libs/libfmod.so` must be compiled beforehand; `libs/flibs` and `libs/plibs` are pure-Python packages that wrap it. The FLEX, linear/nonlinear Eliashberg, and SC chi calculations are executed via OpenMP-parallel Fortran routines.

== Building the Fortran kernel

Run one of these in `libs`:

#table(
  columns: (auto, 1fr),
  [`make auto`], [*Recommended on a cluster.* Detects everything: the toolchain from what is installed and loaded, and the CPU targets from the compute nodes themselves. Builds one `libfmod_<target>.so` per distinct CPU type found, and `libs/flibs/_loader.py` picks the best one the running node supports at import.],
  [`make`], [One library for one target: `ARCH` (default: the portable AVX2 baseline of the build host's vendor, `znver3` on AMD and `x86-64-v3` on Intel), `FC` and `SL` as detected or as given on the command line.],
  [`make dispatch`], [A fixed AVX2 + AVX-512 pair (`znver3` + `znver4` on AMD, `x86-64-v3` + `x86-64-v4` on Intel; `DISPATCH_ARCHS` overrides it), i.e. `make auto` without the detection step.],
)

Explicit settings always win: `make FC=ifx SL=MKL ARCH=znver4`, or `make auto ARCH_LIST="znver3 znver4"`.

*What is detected, and where from.* `FC` is the first of `ifx`, `ifort`, `flang`, `amdflang`, `gfortran` on `PATH`; `SL` follows the modules actually loaded (`MKLROOT` $arrow$ MKL, `AOCL_ROOT` $arrow$ AOCL, otherwise OpenBLAS) — selecting a library whose root is not set is what used to yield a `libfmod.so` with unresolved BLAS symbols that failed only at import. The *CPU targets*, however, must not come from the build host: `-march=native` on a Zen 4 login node produces AVX-512 code that dies with SIGILL on a Zen 3 compute node. `libs/detect_arch.sh` therefore asks each node in the scheduler's inventory (Slurm `sinfo`, PBS `pbsnodes`, else just this machine) what `gcc -march=native` calls its CPU, and prints the distinct targets. Intel nodes are reported by their x86-64 level rather than their generation name (`x86-64-v4` with AVX-512 F/BW/CD/DQ/VL, `x86-64-v3` with AVX2/FMA/BMI2), so that per-generation extensions (AVX-VNNI, AMX, ...) never hand an older node a library built for a newer one; `x86-64-v3` also runs on Zen, so one library covers a mixed AMD/Intel cluster:

```
$ sh libs/detect_arch.sh --table
NODE         CPU                                      -march
salmon       AMD Ryzen 7 5700G with Radeon Graphics   znver3
unagi        AMD Ryzen 9 7900 12-Core Processor       znver4
yakisoba     AMD Ryzen 9 7900X 12-Core Processor      znver4
```

The result is cached in `libs/.arch_cache`; `sh libs/detect_arch.sh --refresh` re-probes, `--local` limits it to the current machine. A node that cannot be probed does not silently drop out: it contributes the portable `znver3` default and a warning. Setting `Features=znver3` / `znver4` in `slurm.conf` makes the probe a plain `sinfo` query with no job submission at all — worth doing if you administer the cluster.

At run time nothing has to be selected by hand: `_loader.py` reads the CPU flags and loads the highest-ISA `libfmod_*.so` the machine can actually execute, falling back to a plain `libfmod.so`. Each target is chosen only when the CPU has its whole set of flags (for `x86-64-v4`, all of AVX-512 F/BW/CD/DQ/VL). An Intel generation name can be built as well (`make ARCH=sapphirerapids`); it is then chosen only on nodes with that generation's extensions (AMX, AVX512-BF16, ...). MKL picks its instruction set by itself on Intel CPUs. For MKL on AMD also export `MKL_ENABLE_INSTRUCTIONS=AVX512` (Zen 4) or `AVX2` (Zen 3), and set `OMP_NUM_THREADS` from the scheduler (`$SLURM_CPUS_PER_TASK`) rather than hard-coding it.

*Building with gfortran.* gfortran rejects calls to one external procedure with different integer kinds within a file, so the integer arguments of BLAS/LAPACK are kept to one kind per file (`c_int64_t` in `fmod.f90`). With MKL (`SL=MKL`) the GNU threading layer is linked (`-lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core`). A gfortran build runs its OpenMP on libgomp, while the MKL build of numpy uses libiomp5: in one Python process the two runtimes keep their idle threads spinning and fight for the cores, which can slow a run down by orders of magnitude. `_loader.py` therefore sets `OMP_WAIT_POLICY=PASSIVE` before loading a `libfmod` linked against libgomp. libgfortran also holds Fortran output back until exit when stdout is not a terminal, so `GFORTRAN_UNBUFFERED_PRECONNECTED=y` is set as well (the FLEX iteration log then shows up as it runs; the AiiDA interface judges convergence from it). A value set by the user takes precedence in both cases. ifx/ifort builds share libiomp5 with numpy and are not affected.

*CMake.* `libs/CMakeLists.txt` builds the same thing with the same detection, for people who prefer it or need it for an IDE:

```bash
cmake -S libs -B build && cmake --build build -j     # = make auto
cmake -S libs -B build -DARCH=znver4                 # = make ARCH=znver4 (one libfmod.so)
cmake -S libs -B build -DFC_COMPILER=gfortran -DSL=MKL -DARCH_LIST="znver3;znver4"
```

Empty knobs mean auto: `FC_COMPILER` from what is on `PATH`, `SL` from `MKLROOT`/`AOCL_ROOT`, and the target list from `detect_arch.sh`. The libraries are written into `libs/`, not the build tree, because that is where `_loader.py` looks. The C++ compiler is only required on the pocketfft path — the MKL path does not enable the CXX language at all, so an oneAPI install with `ifx` but no `icpx` still configures. The makefile remains the reference build: with `flang` plus AOCC modules loaded it also sanitises `CPATH` and `LIBRARY_PATH` for its child processes, which CMake does not do.

= Input Files

== File Format Selection (`ftype`)

The format of the input Hamiltonian file is specified by the `ftype` variable. Assign the file name or directory path (without extension) as a string to `fname`.

#table(
  columns: (auto, 1fr),
  [*ftype*], [*File format*],
  [0], [Directory `{fname}/` containing `ham_r.txt`, `irvec.txt`, and `ndegen.txt`],
  [1], [The file given by `fname` (custom format; give the full name including its extension, e.g. `.input`, `.hop`)],
  [2], [`{fname}_hr.dat` file (Wannier90 default hopping file). *This is the standard choice when interfacing with first-principles codes.*],
  [3], [Text format of ecalj's (non-orthogonal) MLO basis, holding $H(R)$ and the overlap $S(R)$],
  [4], [Directory `{fname}/` containing `HamRsMLO` and `HamiltonianPMTInfo` (binary output of ecalj `job_mlo`). Same non-orthogonal MLO basis as ftype=3; the overlap $S(R)$ is read alongside $H(R)$. For ftype=3,4 every orbital is scaled to unit norm on reading ($H(R), S(R) -> D^(-1\/2) (dot) D^(-1\/2)$ with $D = "diag" S(R=0)$; the bands are unchanged; ecalj's MLOs are not normalized, e.g. Cu: $S_(a a)(0) = 0.016$ for s, $0.565$ for d)],
  [Other], [`Hopping.dat` file (ecalj hopping file)],
)

*Example for Wannier90 users:*

```python
fname = 'inputs/Sr2RuO4'   # path without the _hr.dat suffix
ftype = 2
```

== Spin-Orbit Coupling (`sw_soc`)

Setting `sw_soc = True` enables spin-orbit coupling (SOC). For SOC-active systems, the Hamiltonian must include spin degrees of freedom (i.e., the orbital dimension is doubled). FLEX (option=14) and the linearized Eliashberg equation (option=15) have SOC-aware implementations (`calc_flex_soc` / `calc_lin_eliash_soc` are dispatched internally). Both also work with the IR basis and a non-orthogonal MLO basis. The nonlinear Eliashberg solver (option=23) and SC chi (option=12,13) currently do not support SOC. The spin ordering is by default all orbitals spin up, then all spin down; for another ordering define `slist` (spin $plus.minus 1$ of each orbital) and `invs` (1-based opposite-spin partner).

= Bravais Lattice Setting (`brav`)

The `brav` variable specifies the type of Bravais lattice. This is used for computing reciprocal lattice vectors and auto-generating symmetry lines. Settings compatible with Quantum ESPRESSO (QE) and Wannier90 outputs are provided.

#table(
  columns: (auto, 1fr),
  [*brav*], [*Lattice type (corresponding code)*],
  [0], [Simple lattice (simple cubic / tetragonal / orthorhombic)],
  [1], [Face-centered cubic (QE default, equivalent to ibrav=2)],
  [2], [Body-centered cubic (QE default, equivalent to ibrav=3)],
  [3], [Hexagonal],
  [4], [Trigonal (QE ibrav=5)],
  [5], [Base-centered],
  [6], [Face-centered cubic (conventional orientation)],
  [7], [Body-centered cubic (conventional orientation)],
  [Other], [Monoclinic],
)

= Calculation Modes (`option`)

The type of calculation is selected by setting the integer variable `option`. The `CalcMode` IntEnum can also be used directly (e.g. `option = CalcMode.LIN_ELIASHBERG`). Each mode is described below.

== option=0: Band Plot (`BAND`)

Calculates and plots the band dispersion $E_n(bold(k))$ along a symmetry line. The symmetry path is specified via `k_sets` and `xlabel`, or auto-generated from the `brav` setting if not defined.

== option=1: Density of States (`DOS`)

$ "DOS"(omega) = - 1/pi sum_(bold(k), n) "Im" G^0_n (bold(k), omega + i delta) $

Plots both the total DOS and orbital-projected partial DOS.

== option=2: 2D Fermi Surface (`FERMI_2D`)

Draws the Fermi surface in the $k_x$-$k_y$ plane at a specified $k_z$ slice (default: $k_z = 0$, adjustable via `kz`). The chemical potential is determined self-consistently from the electron count `fill`, and the contour $E_n(bold(k)) = mu$ is extracted via skimage's `find_contours`. A rotation matrix `RotMat` can optionally be used to rotate the Fermi surface.

== option=3: 3D Fermi Surface (`FERMI_3D`)

Renders the three-dimensional Fermi surface as a polygon mesh using the Marching Cubes algorithm. The display scale along each axis can be adjusted with `kscale`.

The isosurface is extracted on the parallelepiped $k_("frac") in [-1\/2, 1\/2]^3$, i.e. on the *primitive cell* of the reciprocal lattice. For a primitive lattice that cell is the first Brillouin zone (the Wigner-Seitz cell drawn as the BZ frame), but for a centred (body/face/base) or skewed lattice it is not, and the Fermi surface would stick out of the zone. At `kscale=1` (the default) the part of the sheet outside the zone is therefore translated by a reciprocal lattice vector and cut at the zone faces, folding it back into the first BZ (the translates of the parallelepiped tile $bold(k)$-space, so no area is duplicated or lost). With `kscale` set to anything else you are deliberately looking at another $bold(k)$-window, and the sheet is left where the sampling put it.

With `color_option=ColorMode.GAP` (3), the surface is colored by $"Re"[phi(bold(k))]$, the *same* pairing form factor that drives the Eilenberger calculations (`gap_sym`/`delta0`, or `eil_gap_orbital`/`eil_gap_file` if set) — a quick way to check the actual gap (sign, nodes, anisotropy) on the real 3D Wannier Fermi surface (all sheets/$k_z$, not just the fixed-$k_z$ cut used by the vortex/surface solvers).

== option=4: Spectral Function (`SPECTRUM`)

Computes the electron spectral function from the imaginary part of the trace of the single-particle Green's function,

$ A(bold(k), omega) = -1/pi "Im" "Tr" [S(bold(k)) G(bold(k), omega + i delta)] $

and plots it along the symmetry line ($S = 1$ in an orthogonal basis). With `sw_self=True`, the self-energy in `self_en_wannier.npz` written by FLEX is Fourier-interpolated onto the $bold(k)$ path and continued to the real axis with Padé (MaxEnt with `sw_maxent=True`) to include interaction effects. A $Sigma$ solved on the IR basis is continued from its sampled values (`sigma_ir`; see "Frequency grid" under option=14).

== option=5: Boltzmann Transport (`CONDUCTIVITY_BT`)

Calculates transport coefficients from the Boltzmann equation in the relaxation-time approximation. The relaxation-time model is chosen with `tau_mode` (constant $tau$ by default).

Thermoelectric quantities:
- Electrical conductivity tensor $sigma_(i j)$ (unit: S/m)
- Thermal conductivity tensor $kappa_(i j)$ (unit: W/m/K)
- Seebeck coefficient tensor $S_(i j)$ (unit: V/K)
- Power factor $sigma S^2$ (unit: $upright(W \/ m \/K^2)$)
- Peltier coefficient, Lorenz number

The inverses in $S = -sigma^(-1) dots$ and the like are taken on the conducting directions only. A strictly two-dimensional model, or $N_z = 1$ (only $k_z = 0$, where $v_z = 0$), has $sigma_(z z) = 0$; the in-plane $S$, Peltier coefficient, $kappa$ and Lorenz number are still obtained, and the components along a non-conducting direction are reported as 0 (with a note).

=== Weak-field Hall response

The same mode also reports the $B$-linear Hall response. The Hall kernel is evaluated for *all three field directions* and resolved by band,

$ S_(a b)^((g)) = 1/2 sum_(bold(k),n) w_bold(k) tau^2 (-partial f \/ partial epsilon) epsilon_(g u v) [v_a (M^(-1))_(b u) v_v - (a <-> b)] $

and the Hall coefficient follows from the full tensors,

$ rho^((1)) = -sigma^(-1) sigma^((1)) sigma^(-1), quad R_H^((g)) = rho^((1))_(b a) \/ B quad (a,b,g " cyclic"). $

Using $sigma_(x y)\/(sigma_(x x) sigma_(y y))$ instead would assume $sigma$ is diagonal in $x y z$, which is false for any cell whose axes are not mutually orthogonal — a monoclinic cell has $sigma_(x y) eq.not 0$ already at $B=0$. If $sigma$ is singular (a strictly two-dimensional band has $sigma_(z z)=0$), the in-plane $2 times 2$ block is used for the field direction normal to it, and the other two directions are reported as undefined.

Output quantities:
- Hall coefficient $R_H$ for $bold(B) parallel x, y, z$ (unit: $upright(m^3 \/ C)$)
- Hall carrier density $n_H = 1\/(e |R_H|)$ (unit: $upright(c m^(-3))$), labelled `electron` or `hole` from $"sign"(R_H)$
- Hall mobility $mu_H = R_H sigma_(a a)$ (unit: $upright(c m^2 \/ V s)$) and $omega_c tau$ per tesla
- $sigma^((1))_(x y)\/B$ (unit: $upright(S \/ m T)$)
- The field range over which the weak-field expansion assumed by this mode ($omega_c tau < 0.1$) remains valid

#block(fill: luma(240), inset: 8pt, radius: 4pt, width: 100%)[
  $n_H = 1\/(e|R_H|)$ equals the true carrier density *only for a single closed Fermi-surface sheet*. In a multi-band metal it is a two-fluid average, $R_H = (p mu_h^2 - n mu_e^2)\/(e(p mu_h + n mu_e)^2)$, which can diverge or change sign near compensation. Compare it against the Luttinger volumes printed alongside (below) before reading it as a density.
]

=== $bold(k)$-mesh convergence

$sigma^((1))$ carries one more $bold(k)$-derivative than $sigma$ (through $M^(-1)$) and its integrand changes sign around the Fermi surface, so it converges far more slowly: on an under-converged mesh $R_H$ can come out with the *wrong sign*. Setting `hall_mesh_list` runs the whole calculation on a sequence of meshes and prints a convergence table, flagging a sign flip or a change above 5% between the two finest meshes.

Two diagnostics are printed per mesh:
- $N_"eff"$ — the participation ratio of the $-partial f\/partial epsilon$ weights, i.e. the effective number of states carrying the Fermi window. Everything is an average over $N_"eff"$ states, so this sets the noise floor; a warning is issued below $10^3$.
- The cancellation ratio $|sum| \/ sum |dot|$ of the Hall integrand. Where this is small the total is a tiny residue of large opposing contributions (near compensation), and the sign may still be mesh-dependent; a warning is issued below $10^(-2)$.

=== Band-resolved decomposition

$sigma$ and $sigma^((1))$ are both plain band sums, so each Fermi-surface sheet is also reported separately: its Luttinger-volume carrier density (the volume of each closed pocket, computed here), its $sigma_(x x)$, its own $R_H$ and $n_H$, and its mobility. Comparing a sheet's own $n_H$ against its Luttinger volume shows whether it behaves as one simple closed pocket; the gap between the total $n_H$ and $n_e\/n_h$ is then the compensation the experiment is actually seeing.

== option=6: Linear Response Transport (`CONDUCTIVITY_PT`)

Calculates electrical conductivity and related quantities based on the Kubo formula (linear response theory). The parameter `delta` corresponds to an effective relaxation time $tau approx hbar / delta$.

The optical conductivity $sigma(omega)$ is further post-processed into (diagonal components $x x, y y, z z$):

- the dielectric function $epsilon(omega) = 1 + i sigma(omega) \/ (epsilon_0 omega)$,
- the plasma frequency (unscreened, from the $f$-sum rule, and screened, from $"Re" epsilon_(x x) = 0$),
- the normal-incidence reflectivity $R(omega) = |(N - 1) \/ (N + 1)|^2$ with the complex refractive index $N = sqrt(epsilon)$,
- the perceived metallic colour: $R(omega)$ folded with the CIE 1931 2° standard observer under illuminant D65 gives the XYZ tristimulus values, converted to sRGB (hex, 8-bit RGB, chromaticity $x y$, luminance $Y$). This is the standard route for first-principles metal colours (Prandini _et al._, npj Comput. Mater. *5*, 129 (2019)); the `poly` entry uses the directional average of $R$, i.e. the polycrystalline colour.

A non-metal reflects almost nothing ($R approx 0.04$ for $n approx 1.5$) and its specular colour is nearly neutral, so its colour is carried by the light that travelled *through* the material. From the absorption coefficient $alpha = 2 omega kappa \/ (planck.reduce c)$ two further routes are reported:

- the *transmission colour* of a slab of thickness $d$, $T = (1-R)^2 e^(-alpha d) \/ (1 - R^2 e^(-2 alpha d))$ (both surfaces, incoherent multiple reflections; it reduces to $T = (1-R)\/(1+R)$ for $alpha = 0$). Note that this colour is *not* an intrinsic material property: it depends on $d$.
- the *body colour* of a powder or pigment, from the Kubelka--Munk two-flux diffuse reflectance $R_infinity = 1 + K\/S - sqrt((K\/S)^2 + 2 K\/S)$ with $K\/S = alpha l$, where the effective scattering length $l$ is the single free parameter (it absorbs the convention factor between $alpha$ and the Kubelka--Munk $K$, and the particle size).

A colour computed from the absorptance $A = 1 - R - T$ is also printed, but since $A$ is the complement of what leaves the sample it is reported as its *negative* (`hex_neg`, i.e. $1 - "rgb"$ in encoded sRGB). The two lengths are set by the parameters `color_thick` ($d$) and `color_scatt` ($l$). A metal is opaque, so its $T$ and $R_infinity$ colours come out black -- correctly so, and they are then dropped from the figure's swatches automatically.

Outputs are `optical.dat` ($sigma$, $epsilon$, $R$, $alpha$, $T$, $R_infinity$ spectra, with the colours recorded in the header), `dielectric.pdf/png` and `reflectivity.pdf/png` (with colour swatches). Set `Emax` above 3.3 eV so that the visible window 380--780 nm (1.59--3.26 eV) is covered; a warning is printed otherwise.

On accuracy: a metal's colour is dominated by the smooth Drude response, so moderate errors barely move the hue, whereas a non-metal's hue *is* the position of the absorption edge. An error in the model band gap translates directly into a hue shift, and excitonic effects and phonon-assisted indirect transitions are absent from this Kubo-bubble calculation; the `delta` broadening also blurs the edge. Check convergence, and the model itself, before quoting these colours.

== option=7: Spin Susceptibility Spectrum (`CHIS_SPECTRUM`)

Calculates the RPA spin susceptibility

$ chi_s (bold(q), omega) = frac(chi^0 (bold(q), omega), 1 - S chi^0 (bold(q), omega)) $

along the symmetry line and plots the imaginary part as a paramagnon spectrum. The interaction matrix $S$ is generated from the on-site parameters `U` and `J`. Results are written to `chis_spec.png`.

== option=8: Spin Susceptibility at a $bold(q)$ Point (`CHIS_QPOINT`)

Calculates and plots the frequency dependence $chi_s (bold(q), omega)$ at the wave vector $bold(q)$ specified by `at_point`.

== option=9: $bold(q)$-Space Spin Susceptibility Map (`CHIS_QMAP`)

Plots the spatial distribution of $chi_s (bold(q), omega_0)$ in the $k_x$-$k_y$ plane at the energy $omega_0$ specified by `Ecut`. Useful for visualizing nesting vectors. Output: `chi0map.png` and `chismap.png`.

== option=10: Pairing Susceptibility Spectrum (`PHI_SPECTRUM`)

Calculates the pairing susceptibility $phi(bold(q), omega)$ along the symmetry line and plots it (`phi_spec.png`).

== option=11: $bold(q)$-Space Pairing Susceptibility Map (`PHI_QMAP`)

Plots the $bold(q)$-space distribution of the pairing susceptibility at the energy specified by `Ecut` (`phimap.png`). Use `sw_omega` to switch between real and Matsubara frequency.

== option=12: SC-State Spin Susceptibility Spectrum (`CHIS_SPECTRUM_SC`)

Assumes a non-zero gap function $Delta(bold(k))$ and constructs the irreducible susceptibility in the superconducting state, which includes the anomalous bubble $F(bold(k))$:

$ chi_0^"SC" = chi^"GG" plus.minus chi^"FF" $

The dynamic spin susceptibility $chi_s^"SC" = chi_0^"SC" \/ (1 - S chi_0^"SC")$ is then obtained by RPA and plotted along the symmetry line. The initial gap symmetry is set by `gap_sym` (negative values for triplet). Output: `chis_sc_spec.png`.

The initial gap amplitude `delta0` can be specified in two ways (see also the parameter reference below):
- A single float (e.g. `delta0=1.e-2`): a single-band gap shape is generated internally and scaled to this maximum amplitude across all bands.
- A list of per-band values of length `Norb` (e.g. `delta0=[0.,0.2,0.3,-0.1,0.]`): enables *multi-gap* mode, where each band gets its own amplitude and sign. Mixing signs lets you represent sign-changing gaps such as $s^plus.minus$, as is common in multiband (e.g. Fe-based) superconductors.

== option=13: SC-State Spin Susceptibility at a $bold(q)$ Point (`CHIS_QPOINT_SC`)

Same SC framework as option=12, but evaluates $chi_s^"SC"(bold(q), omega)$ at the single $bold(q)$ specified by `at_point`. `delta0` is specified the same way as in option=12.

This mode additionally produces:

- The BdG (Bogoliubov–de Gennes) band dispersion along $bold(k) = (0,0,0) arrow.r (0,0.5,0)$ (`BdG_band.png`)
- The orbital-traced $chi_s^"SC"(omega)$ (`chisq.png`, `chis_sc.dat`)
- The orbital-resolved $chi_s^"SC"(omega)$ (`chisq_orb.png`, `chis_scorb.dat`)

== option=14: FLEX Self-Energy (`FLEX`)

Solves the FLEX equations self-consistently to obtain the electron self-energy $Sigma(bold(k), i omega_n)$ in Matsubara frequency space. With `sw_out_self=True` it is saved to `self_en.npz` (recording the frequencies `iw` and the temperature), `self_en_wannier.npz` (Wannier real-space form, read by option=4 and `tau_mode='flex'`) and, on the uniform grid, `sigma.bin` (seed for `sw_in_self`).

- For a non-orthogonal MLO basis (`ftype=3,4`), $G^(-1) = (i omega + mu) S - H - Sigma$ and the particle number uses $"tr"(S G)$. `U` and `J` are local interactions in the normalized MLO orbitals ($S_(a a)(R=0) = 1$, set when ftype=3,4 is read), not in an orthogonalized basis. Without the normalization U would act with a strength of about $U \/ S_(a a)(0)^2$ (some 4000 times U for the s orbital of Cu). Linear and nonlinear Eliashberg and self-energy-corrected spectra use the same convention. Self-energy files do not store the overlap, so use the same basis and model when loading them. Transport with `tau_mode='flex'` still requires an orthogonal basis.
- `sw_out_self=True`: write the self-energy to file
- `sw_in_self=True`: load the previous self-energy (`sigma.bin`) as the initial guess
- `sigma_in_scale`: factor multiplying the $Sigma$ seed loaded via `sw_in_self=True` (default 1.0). For $U$-annealing (continuation that raises $U$ stepwise, seeding each step from the `sigma.bin` converged at the smaller $U$) the weak-coupling relation $Sigma prop U^2$ suggests $(U_"new"\/U_"old")^2$. Unlike temperature annealing the Matsubara grid is unchanged, so no seed interpolation is needed (`sigma.bin` has no header — keep `Nw` identical between runs)
- Temperature annealing: for stepwise cooling at fixed `Nw` ($T_"new"\/T_"old" gt.eq 0.8$ or so) the raw reuse of `sigma.bin` is sufficient. Raw reuse amounts to a seed with the frequency axis compressed by $T_"old"\/T_"new"$, and this compression approximately mimics the growth of $Sigma$ upon cooling — near criticality it converges faster than a faithful interpolation. When `Nw` must change (the one case raw reuse cannot handle), the Python helper `plibs.regrid_sigma_bin(temp_new=..., Nw_new=..., w_scale=...)` re-interpolates `sigma.bin` onto the new mesh. The original temperature, `Nw` and orbital count are taken automatically from a metadata record appended at the end of `sigma.bin` by recent `io_sigma` versions (a trailing record, invisible to older readers); explicitly passed values that contradict the footer raise an error (files without the footer require explicit `temp_old`/`Norb`/`Nw_old`). The `sw_in_self` reader also checks the footer: a `sigma.bin` whose mesh (`Nw`/orbital count/k-points) disagrees with the run aborts with an error (a Fortran read silently accepts a record that is larger than requested, scrambling the layout — this check prevents that accident); a differing temperature alone is treated as T-annealing and only reported. The SOC path (`mkself_soc`) now writes the same single-record format, so the helper applies to its files as well (the retired element-wise SOC files are no longer readable) (cubic spline with the conjugate-symmetric extension $Sigma_(l m)(-i omega) = Sigma_(m l)(i omega)^*$ to negative frequencies; the outermost stored frequency, which carries the wrap-around artifact of the circular convolution, is discarded, and points beyond the cutoff are filled with a $c_0 + c_1\/(i omega)$ tail fitted on an interior window). For near-critical cooling use `w_scale=temp_old/temp_new` (the same compression as raw reuse); `w_scale=1` (default) gives the faithful interpolation
- `m_diis_num`: DIIS history length (default 5 if undefined)
- `sw_rescale_flex=True`: if the Stoner factor $S = max_bold(q) "eig"[chi_0(bold(q), 0) hat(S)]$ reaches 1 during the iteration, uniformly rescale $chi_0(bold(q), i nu_n)$ by $(1 - 10^(-4))\/S$ to avoid the magnetic-instability divergence
- `sw_chi0_tail=True`: evaluate $chi_0$ with a tail correction (default True; no-SOC path of option=14,15,23; also in a non-orthogonal MLO basis, and ignored on the IR basis, which has no truncation error).  The bubble is split bilinearly as $chi_0 = "conv"[G] - "conv"[G_0] + chi_0^"ref"$: the slowly decaying $1\/(i omega)$ reference part (the $G_0$ bubble) is evaluated exactly in imaginary time, and only the fast-decaying residual ($G - G_0 tilde 1\/omega^2$) passes through the FFT convolution.  The Matsubara truncation error improves from $O(1\/N_w)$ to $O(1\/N_w^2)$ (validated against the exact Lindhard function: convergence order 1.0 -> 2.0).  The $chi_0$ stage costs about 3x; effective for $N_w gt.eq 64$ (at very small $N_w$ the tau-discretization prefactor of the reference bubble can dominate).  The improved $chi_0$ propagates automatically into $chi_s$/$chi_c$, the FLEX vertex $V_sigma$ and the pairing vertex $V_Delta$.

*Notes on the built-in approximations* (all are deliberate choices made for computational cost and numerical stability):

- The `sw_rescale_flex` rescaling shrinks $chi_0$ uniformly at ALL momenta and frequencies, not just the diverging mode. It is, however, only a transient aid for the early iterations: the self-energy feedback normally pulls the Stoner factor below 1 by convergence, the rescaling becomes inactive, and the converged solution is then a genuine unrescaled FLEX fixed point (no bias). If the rescaling is still active at the final iteration, the normal state is magnetically unstable at this $T$/$U$ ($T < T_N$ within FLEX) and stdout prints `[FLEX] WARNING: Stoner factor >= 1 at the final iteration`; such results correspond to an artificially weakened interaction and should not be used quantitatively (raise $T$, or seed from a converged higher-$T$ $Sigma$ via `sw_in_self`). The clamp target just below 1 ($1 - 10^(-4)$) is deliberate: the resulting $Sigma$ overshoot is the mechanism that kicks the loop off the rescaled manifold. A softer clamp (e.g. 0.98) makes near-critical systems "converge" quickly to a spurious rescaled fixed point even when a genuine $S < 1$ solution exists — the extra iterations caused by the oscillation are the price of this safety.
- A static part of the self-energy is subtracted at every iteration (wrapper argument `sub_sigma`; `1`: Hermitian part of $Sigma(bold(k), i omega_0)$, the default; `2`: frequency average; `0`: no subtraction). This is a prescription that absorbs the static band shift into the chemical-potential adjustment and stabilizes the Fermi surface; it is NOT a strict Hartree--Fock-only subtraction (the value at $i omega_0$ also contains low-energy dynamical components). Use `sub_sigma=0` if absolute band shifts are of interest.
- The Matsubara-frequency convolutions ($chi_0$, $Sigma$, Eliashberg kernel) are circular FFT convolutions of length $2 N_w$ with a sharp cutoff and, by default, no high-frequency tail correction (only the farthest point of $V$ is approximated by the bare vertex). The effective cutoff $omega_c = (2 N_w - 1) pi T$ shrinks proportionally to $T$, so a temperature scan at fixed $N_w$ carries an $O(1\/N_w)$ systematic drift. Choose $N_w$ so that $omega_c$ safely exceeds the bandwidth and `U` even at the lowest temperature. `sw_chi0_tail=True` improves the $chi_0$ truncation error to $O(1\/N_w^2)$. The tail error of the $Sigma$ convolution itself is dominated by an $omega$-independent constant (an HF-like static shift) that is absorbed by `sub_sigma`/the chemical-potential adjustment, and the Eliashberg-kernel tail error only enters the uniform ($s$-wave) component of the gap and vanishes by symmetry for sign-changing gaps — hence no separate corrections are implemented for those ($V_Delta$ and $V_sigma$ still benefit through the corrected $chi_0$).

=== Frequency grid: uniform Matsubara and the IR basis (`freq_grid`)

With `freq_grid='ir'`, the frequency axis of option=14, 15, 16, 23 ($chi_0$, $Sigma$, the chemical-potential search and the linearized and nonlinear Eliashberg kernels; with and without SOC, and `sw_nosym`) is handled by the sparse sampling of the IR (intermediate representation) basis instead of $N_w$ uniform Matsubara frequencies. sparse-ir (1.x or 2.x; checked with 2.1.6) is required.

- The number of frequencies is fixed by the basis (it depends on $Lambda = omega_"max"\/T$ and the accuracy `ir_eps`; typically 15--30) and `Nw` is not used. There is no truncation error as on the uniform grid, so `sw_chi0_tail` is unnecessary and switched off automatically.
- `ir_wmax` (default None): spectral cutoff $omega_"max"$ [eV]; None uses $2 max_bold(k) |epsilon_bold(k) - mu|$. For the RPA $chi_0$ the band width suffices, but FLEX spreads spectral weight beyond the band, so about twice that is needed (000AsP, $T = 0.02$ eV: with $omega_"max" = 2.88$ eV $mu$ is off by 0.27 meV, with 6.0 eV it agrees with the $N_w -> infinity$ extrapolation of the uniform grid to $6 times 10^(-6)$ eV).
- `ir_eps` (default $10^(-8)$): accuracy of the basis.
- Accuracy (verified): $chi_0$ agrees with the exact Lindhard sum to $10^(-8)$, the FLEX $mu$ and Stoner factor with the extrapolated uniform grid to below $10^(-5)$, the linearized Eliashberg $lambda$ to $10^(-7)$--$10^(-8)$, the converged nonlinear Eliashberg gap to below $10^(-4)$ (one orbital, $16 times 16$, $T = 0.02$ eV, $lambda = 1.34$: the uniform grid gives a gap 3% too small even at $N_w = 4096$). On the uniform grid FLEX keeps an $O(1\/N_w)$ error from the truncation of $Sigma$ even with the $chi_0$ tail correction, and $mu$ can be off by about 1 meV at $N_w = 1024$. The run time drops by a factor of tens.
- Restrictions: `sigma.bin` (`sw_in_self`, the `sigma.bin` part of `sw_out_self`) and `sub_sigma=2` are uniform-grid only.
- Files: `self_en.npz` holds $Sigma$ on the sampling frequencies together with the frequencies (`iw`); a file from another grid is rejected (`sw_from_file`, option=16). `self_en_wannier.npz` keeps its format, the first 64 uniform $omega_n = (2n+1) pi T$, evaluated through the basis (the constant $Sigma(i omega -> infinity)$, which the basis cannot represent, is estimated from the Hermitian part at the three highest sampling frequencies together with the $1\/omega^2$, $1\/omega^4$ terms and taken off before the expansion). The gap (`gap.npy`, `gap_wannier.npz`) stays on the sampling frequencies, which are recorded in `gap_iw.npy` and in `iw` of `gap_wannier.npz`; the $omega_n -> 0$ extrapolation uses the recorded frequencies.
- *Padé continuation*: Padé is ill-conditioned: random errors of relative size $10^(-10)$ in the input change the spectrum by tens of percent. IR data evaluated through the basis at uniform $omega_n$ are more accurate on the Matsubara axis than the uniform grid (two-orbital square-lattice model (`square_2orb.hop`), $16 times 16$, $T = 0.05$ eV: the deviation of $G$ from the $N_w -> infinity$ extrapolation is $9 times 10^(-6)$, against $2 times 10^(-4)$ for $N_w = 4096$), but the basis-truncation error left in them is smooth and systematic in $omega_n$, and Padé continues it into spurious real-axis structure (negative spectra, lost peaks). Padé on the IR grid therefore takes *the sampled values themselves* (all $N_w$ points): `Gpade.dat` of option=16; with `sw_out_self`, `self_en_wannier.npz` stores the sampled values `sigma_ir`, `iw_ir` next to the uniform-grid `sigma`, and option=4 with `sw_self` and `tau_mode='flex'` use them. The accuracy is set by `ir_eps` (hardly by the FLEX convergence threshold); use `ir_eps=1e-14` or so when continuing with Padé (it adds only a few points). Relative L2 difference from the Richardson-extrapolated uniform-grid spectra (median / maximum over k): `square_2orb.hop` with `ir_eps=1e-8`: 12% / 45%, `1e-12`: 4% / 34%, `1e-14`: 3% / 6% (uniform $N_w = 1024$: 1.5% / 3%); 000AsP ($16 times 16 times 4$, $T = 0.05$ eV) with `1e-12`, `1e-14`: 4% / 13% ($N_w = 256$: 3% / 14%). The former route through the uniform-grid evaluation gave 25--30% / 50--70% at `1e-8`, with negative spectra. MaxEnt (`sw_maxent=True`) uses the uniform-grid evaluation and is stable for any `ir_eps`. Quantities on the Matsubara axis ($mu$, $lambda$, the gap, $chi$) are not affected.

== option=15: Linearized Eliashberg Equation (`LIN_ELIASHBERG`)

Solves the linearized Eliashberg equation as an eigenvalue problem using the effective pairing interaction from RPA or FLEX. The largest eigenvalue $lambda$ and the corresponding gap function $Delta(bold(k), i omega_n)$ are obtained. The superconducting transition temperature $T_c$ corresponds to $lambda = 1$. The default solver is a power method with shift + deflation; setting `arnoldi_m > 0` enables an Arnoldi solver.

- `sw_self=False`: no self-energy corrections (pure RPA)
- `sw_self=True`: uses the FLEX self-energy (FLEX is solved on the spot)
- `sw_from_file=True` (with `sw_self=True`): loads the self-energy from `self_en.npz` written by option=14 without re-running FLEX (rejected when the frequency grid or the temperature differ)
- `gap_sym`: sets the initial symmetry of the gap function (see below)
- *With SOC* (`sw_soc=True`) the sign of `gap_sym` fixes the *inversion parity* of the gap (positive: even, negative: odd), not the spin structure: spin singlet and triplet mix inside a parity sector once SOC is present. Such a sector also holds the opposite-spin, odd-frequency partner of the state you asked for, so the Arnoldi solver returns the branch whose spin character matches `gap_sym` instead of the largest eigenvalue of the sector; this keeps $lambda$ continuous as the SOC strength goes to zero. *Use the Arnoldi solver (`arnoldi_m > 0`) with SOC* — the power method has no branch tracking and converges to the dominant state of the sector, which may be the odd-frequency partner. Without SOC singlet and triplet are exact sectors and are projected out as before (detected at run time), so `gap_sym` keeps its usual singlet/triplet meaning there
- Gap functions are written to `gap.npy` and to one file per orbital pair — `gap_11.dat`, `gap_12.dat`, … (`gap_{i j}.dat`, 1-based orbital indices) (with `sw_out_self=True` also `gap_iw.npy` and the Wannier real-space `gap_wannier.npz`)

=== Linearized Eliashberg without symmetry assumptions (`sw_nosym`)

With `sw_nosym=True`, option=15 solves the gap equation on the full Brillouin zone without assuming time-reversal or inversion symmetry (bare $G$, RPA), for situations that break them, such as a Zeeman field or pairs of finite momentum (FFLO). A spinless model (`sw_soc=False`) is spin-doubled automatically. There is no singlet/triplet projection; `gap_sym` only seeds the iteration. The FLEX self-energy is not supported (with `sw_self` the bare $G$ is still used).

- `zeeman_h` (3 components, eV): Zeeman field $H_Z = bold(h) dot bold(sigma)$ ($g mu_B B\/2$ absorbed). The chemical potential is re-solved at the same `fill`.
- `fflo_q` (list of integer triples): the pair momenta $bold(Q) = (i_(Q x)\/N_x, i_(Q y)\/N_y, i_(Q z)\/N_z)$, e.g. `[[q, 0, 0] for q in range(5)]` scans along $k_x$. $chi_0$ and $V_Delta$ are computed once and reused for every $bold(Q)$; $lambda(bold(Q))$ goes to `lambda_Q.dat`, the gap to `gap_nosym.npz` and `gap_nosym.dat`.
- `nosym_arnoldi_m` (integer, default 60): dimension of the Krylov space; a small $lambda$ may need about 60 or more for $10^(-6)$ accuracy.

== option=16: Gap Function Post-Processing (`GAP_FUNCTION`)

Reads the gap function from `gap.npy` (produced by option=15) and computes:

- The anomalous Green's function $F(bold(k), i omega_n)$
- Analytic continuation to the real-frequency axis via Padé approximation (with `freq_grid='ir'` from the sampling points; `ir_eps=1e-14` or so is recommended, see "Frequency grid" under option=14)

== option=17: Band Filling (`BAND_FILLING`)

A quick diagnostic on the band occupations. For each band it prints the fraction of the $bold(k)$-mesh lying above and below $mu$: the unoccupied fraction is the hole-like side of that band, the occupied fraction the electron-like side, so the pair says at a glance whether a sheet sits near the bottom or the top of its band. The occupied fractions sum to the electron count per unit cell, which is the consistency check against `fill`. Both fractions are Luttinger volumes normalized per unit cell per spin rather than per $upright(c m^3)$: a band with no Fermi surface reads 1 (one electron per cell) or 0, and a band with a pocket reads that pocket's $V_"occ"\/V_"BZ"$ on the corresponding side. `fs_carrier_density`, used by option=5, reports the same volumes multiplied by $g_s\/V_"uc"$; this normalization is the one that pairs with `fill`, the $upright(c m^(-3))$ one is what an experimental carrier concentration is compared against.

== option=18: Cyclotron Mass (`CYCLOTRON_MASS`)

Computes the cyclotron mass via the Onsager relation $m^*_c = (hbar^2 \/ 2 pi)(partial S \/ partial E)$, where $S$ is the Fermi surface cross section. The field direction is set here by the rotation matrix `RotMat` (the run prints the implied $hat(B)$ and its $theta$, $phi$), unlike option=19, which scans `theta`/`phi` directly.

- Phase 1: Scan $S(k_z)$ across $k_z in [0, 0.5]$ (reduced coordinates) with `meshkz=20` points
- Phase 2: Refine extremal $k_z$, compute $partial S \/ partial E$ via central finite difference, and report $m^*_c$ in units of $m_e$

Internally uses skimage's `find_contours` and `marching_cubes`.

== option=19: dHvA Frequency vs Angle (`DHVA`)

Computes the geometric dHvA frequencies from extremal closed Fermi-surface cross sections using the Onsager relation $F = (hbar \/ (2 pi e)) A_"ext"$. Here $A_"ext"$ is an extremal area on a plane perpendicular to the magnetic field. `theta` is the polar angle measured from the Cartesian $z$ axis of the frame in which `avec` is given, and `phi` is the azimuth in the $x y$ plane; with the usual setting $c parallel z$, $theta=0$ therefore means $B parallel c$. The standard `main.py` path uses `phi=0` and samples `theta=0`–$90$ degrees at 40 equally spaced points.

The routine represents the Fermi surface in Cartesian reciprocal space and scans planes normal to the requested physical field direction. It identifies all closed contours on each plane, follows them as the plane is displaced, and retains the minimum and maximum cross sections of every continuous orbit branch. Reciprocal-zone copies are merged. This makes the calculation applicable to non-cubic cells and to tilted fields; `avec` must therefore contain the physical primitive lattice vectors in angstroms. The reported areas are in $angstrom^(-2)$ and frequencies are in tesla.

The result is saved as `dhva_band.png`, with frequency in tesla plotted against `theta`. The terminal also reports the Fermi-surface bands, scan window, number of orbit branches, and number of retained extremal orbits. This option calculates frequencies only: it does not calculate oscillation amplitudes, Dingle factors, Zeeman splitting, or cyclotron masses (use option=18 for the latter).

Open contours and contours whose area changes when the slice window is enlarged are discarded, because they do not define a closed dHvA orbit. A warning about discarded contour pieces is therefore not itself an error, but frequencies should be regarded as converged only after increasing the in-plane mesh `Nx`, the number of displacement planes `meshkz`, and the interpolation grid `grid_mesh` (default 120). This is especially important near $B parallel a b$, near a Lifshitz transition, or for very small pockets.

The geometry and unit conversion are covered by analytic-cylinder tests, including the $F(theta)=F(0)/cos(theta)$ law for a tetragonal lattice. As a material-level check, the 10-orbital FeS Wannier model with experimental lattice constants gives $0.489$–$3.123$ kT for $B parallel c$, consistent with the GGA frequencies of roughly $0.5$–$3$ kT reported for the same field direction in Fig. 3(b) of T. Terashima _et al._, Phys. Rev. B *99*, 134501 (2019). The measured frequencies in that work are smaller, $0.15$–$1.87$ kT. This gap is expected and is not a defect of the calculation, and it is not specific to the iron pnictides and chalcogenides. In a strongly correlated metal a DFT (LDA/GGA) band structure omits the band narrowing and the orbital-dependent band shifts caused by electron correlation, so it systematically overestimates the size of the Fermi-surface pockets while underestimating the mass enhancement and hence the density of states at $E_F$. The same bias carries over to the pairing calculations of this code: a Cooper instability evaluated on an unrenormalized DFT band structure comes out too weak, so the RPA/Eliashberg eigenvalue and $T_c$ are underestimated. Outside the iron-based family the discrepancy can be far larger. In the heavy-fermion antiferromagnet $"CeRhIn"_5$ a band calculation that treats the $4f$ electrons as itinerant gives the large Fermi surface, whereas the dHvA Fermi surface is essentially that of the $4f$-less reference $"LaRhIn"_5$, with cyclotron masses about an order of magnitude heavier (H. Shishido _et al._, J. Phys. Soc. Jpn. *71*, 162 (2002)). In $"Na"_x "CoO"_2$ the six small $e'_g$ hole pockets predicted by the LDA (D. J. Singh, Phys. Rev. B *61*, 13397 (2000)) are absent in photoemission, which sees only the $a_(1g)$ sheet (H.-B. Yang _et al._, Phys. Rev. Lett. *92*, 246403 (2004)); whether local correlations alone remove them is still debated. The benchmark for this option is therefore the frequencies calculated from the same class of band structure, not the raw experimental ones; reproducing the experiment requires a Wannier model built on a correlated band structure (QSGW, LDA+DMFT) or empirical band shifts.

== option=20: Electron Mass (`ELECTRON_MASS`)

Computes $m^* = hbar^2 (partial^2 E \/ partial bold(k)^2)^(-1)$ along the symmetry line (in units of $m_e$).

== option=21: Spectrum with Impurities (`SPECTRUM_IMPURITY`)

Diagonalizes the Hamiltonian of a periodic supercell carrying on-site impurity potentials and unfolds the spectral function $A(bold(k), omega)$ onto the $bold(k)$ of the original cell.

- `imp_cell` (3 components): supercell size per axis; it must exceed twice the hopping range (a matrix of dimension $N_"cell" times N_"orb"$ is diagonalized).
- `imp_sites` (`None` or list of integers): impurity sites in the supercell (0-based, `gen_rlist` order); `None` places them at random with concentration `imp_conc` (at least one).
- `imp_conc`, `imp_seed`: concentration and random seed of the random placement (`None`: a different placement every run).
- `imp_V` (number or $N_"orb" times N_"orb"$ matrix, eV): on-site impurity potential (a number multiplies the unit matrix).

== option=22: CPA Conductivity / Spectrum (`SIGMA_CPA`)

Solves the Coherent Potential Approximation (CPA) self-consistency for an alloy and outputs both the real-frequency spectral function (`cpa_spectrum.png`) and the Matsubara self-energy (`sigma.bin`, `self_en.npz`). The impurity concentration `x_cpa` and the on-site perturbations `VA`, `VB` are currently fixed inside the routine.

== option=23: Nonlinear Eliashberg Equation (`NONLIN_ELIASHBERG`)

Solves the fully nonlinear (self-consistent) SC FLEX-Eliashberg loop. Unlike the linearized solver, $Delta$ is allowed to grow to a finite amplitude below $T_c$, and the SC Dyson equations are iterated together with the FLEX self-energy and anomalous bubble.

The initial gap is now generated automatically inside the solver — there is no need to run option=15 beforehand to produce a `gap.npy` file.

1. The linearized Eliashberg equation is first solved internally to obtain the Stoner factor $S$ and the largest eigenvalue $lambda_"eliash"$.
   - If $S >= 1$ (SDW/CDW instability), the routine stops before entering the nonlinear loop.
   - If $lambda_"eliash" < 1$ ($T >= T_c$, no SC instability), the routine also stops.
   - If `sw_check_only=True`, the routine stops here and only reports $S$ and $lambda_"eliash"$ — useful for quickly scanning temperature to bracket $T_c$ without running the expensive nonlinear loop.
2. The symmetry-correct shape from this linear eigenvector is kept, and its amplitude is rescaled to the BCS weak-coupling value $Delta_0 = 1.764 k_upright(B) T_c$ before entering the nonlinear loop.
3. Each iteration: $Delta_"new" = T \/ N_k sum V_Delta dot F$ → amplitude-direction Newton (secant) acceleration + DIIS shape mixing (falls back to linear mixing pp=0.3 when bypassed) → update $Sigma$ (`sw_self=True`; static part subtracted, mixed with `pp`) → $G, F$ from the Nambu Dyson equation $G = [G_0^(-1) + Delta G_0^dagger Delta^dagger]^(-1)$, $F = -G_0 Delta G^dagger$ → recompute $V_sigma, V_Delta$ from $chi^"GG", chi^"FF"$ (skipped when $V$ is kept fixed)

- `m_diis_num`: DIIS history length ($>= 2$ enables Pulay extrapolation; otherwise linear mixing)
- `sw_self=True`: include FLEX self-energy via the $Sigma$-dressed Green's function
- `sw_from_file=True`: read the self-energy from `self_en.npz` written by option=14
- `sw_check_only`: see "Switch Variables" below
- `sw_out_self=True`: writes the gap to `gap.npy`, `gap_wannier.npz` and, with `sw_self=True`, the SC-state $Sigma$ to `self_en_sc.npz`, `self_en_sc_wannier.npz`
- A diverging iteration (non-finite $Delta$) stops with an error and writes no gap; reaching `itemax` without convergence prints a warning
- Works on the IR basis (`freq_grid='ir'`)
- Amplitude-direction Newton acceleration speeds up convergence of the gap magnitude (see the `sw_amp_newton` comments in `libs/src/ffeliash.f90`)
- See `libs/src/ffeliash.f90` for implementation details

== Quasiclassical Eilenberger Modes (option=24,25,26)

Options 24–26 solve the *quasiclassical Eilenberger equations* of superconductivity — the energy-integrated Gor'kov equations valid when the gap and disorder vary slowly on the scale of the Fermi wavelength ($Delta, T_c, hbar/tau << E_F$). This is the natural framework for $T_c$, the density of states, surface Andreev bound states, and the vortex/vortex-lattice state. The unknowns are the quasiclassical propagators $g(bold(k)_F, bold(r), omega_n)$ (normal) and $f$ (anomalous), parametrized by a single *Riccati amplitude* $a$ ($f = 2a\/(1+a a^*)$, $g = (1-a a^*)\/(1+a a^*)$), which is integrated along straight Fermi-velocity trajectories with a numerically stable Fortran kernel (`riccati_chords`; the $2 times 2$ spin version `matrix_riccati_chords` for the d-vector). The pairing is separable, $Delta(bold(k)_F, bold(r)) = Delta(bold(r)) phi(bold(k)_F)$, with the form factor $phi$ fixed by `gap_sym` and normalized to $⟨ |phi|^2 ⟩_"FS" = 1$, so the coupling `eil_coupling` is the dimensionless $lambda$.

The Fermi surface is shared by all three modes via `eil_fs_kind`: `None` is an isotropic cylinder (analytic angular average), `'iso'`/`'ellipse'`/`'tb'` are model FSs built from `eil_fs_params`, and `'wannier'` builds the real FS and Fermi velocities from the loaded Wannier band (the gap symmetry / multiband structure then comes from `gap_sym`, `delta0`, or `eil_gap_orbital`). The temperature is the global `tempK`/`temp`, and the gap symmetry is the global `gap_sym` (the model-FS routines map the integer index to its continuum harmonic, with $2$ ($s^plus.minus$) $-> s$).

=== option=24: Homogeneous Eilenberger (`EILENBERGER`)

The bulk (spatially uniform) solver. With all `eil_*` sub-mode flags off it self-consistently solves the gap $Delta(T)$ and reports $T_c$; `eil_find_tc=True` brackets $T_c$ by bisection. The sub-mode flags (mutually exclusive) select:

- `eil_imp_sweep=True`: sweep the non-magnetic impurity rate $Gamma$ (`eil_imp_list`) and write $T_c(Gamma)$ to `eilenberger_tc.dat` — the Abrikosov–Gor'kov pair-breaking curve (no suppression for an isotropic $s$-wave by Anderson's theorem; strong suppression for sign-changing gaps). `eil_imp_c` interpolates Born ($-> infinity$) to unitary ($-> 0$) scattering.
- `eil_pauli=True`: Zeeman (Maki) Pauli-limiting sweep — the singlet gap $Delta(h)$, the first-order spinodal/Chandrasekhar–Clogston transition, and the Zeeman-split DOS.
- `eil_spin=True`: the spin-resolved ($2 times 2$) Zeeman response, contrasting a singlet/parallel d-vector ($bold(d) parallel bold(h)$, Pauli-limited) with a perpendicular d-vector ($bold(d) perp bold(h)$, Zeeman-immune).
- `eil_lambda=True`: the superfluid density $rho_s(T)$ and penetration depth $lambda(T)$ (exponentially flat for a full gap, linear-in-$T$ for a nodal gap) → `penetration_depth.dat`.
- `eil_fs=True`: the same on a model FS with Fermi velocities, giving the anisotropic $lambda_(x x), lambda_(y y)$ → `fs_penetration.dat`.
- `eil_free_energy=True`: the condensation free energy $(Omega_s - Omega_n)\/N_0$ vs $T$ (coupling-constant integration) → `free_energy.dat`.

=== option=25: Surface Andreev Bound States (`EILENBERGER_SURFACE`)

Solves the self-consistent gap profile $Delta(x)$ near a specular surface by Riccati integration along reflected trajectories, and (with `eil_ldos=True`) the surface LDOS. The surface orientation is `eil_surf_beta` (for $d$-wave, $0 = [100]$ has no bound state; $pi\/4 = [110]$ produces the zero-energy Andreev bound state, the ZEBS, from the sign change felt on reflection). A Zeeman field `eil_zeeman` splits the ZEBS into $plus.minus h$. With `eil_surf_dvector=True` it instead solves the self-consistent triplet *d-vector texture* at the surface (a dominant + a subdominant component via the spin-matrix Riccati, coupling ratio `eil_dvec_subratio`).

=== option=26: Vortex and Vortex Lattice (`EILENBERGER_VORTEX`)

The inhomogeneous solver around a magnetic vortex. `eil_field` $= B\/H_(c 2)$ selects the geometry: `0` is an isolated vortex in a large circular cell (radius `eil_vort_lxi` in units of $xi$, grid `eil_vort_ngrid`); `>0` is a circular-cell vortex lattice with the Doppler shift. It writes the gap profile $Delta(rho)$ and (with `eil_ldos`) the zero-energy core LDOS (the Caroli–de Gennes–Matricon bound state; the Volovik $sqrt(B)$ DOS in the lattice). A Zeeman field `eil_zeeman` spin-splits the core states. Sub-modes:

- `eil_vort_current=True`: the circulating supercurrent $j_phi(rho)$ → `vortex_current.dat`.
- `eil_vort_field=True` / `eil_vort_maxwell=True`: the self-consistent finite-$kappa$ magnetic field $B(rho)$ / vector potential $bold(A)(bold(r))$ (Maxwell back-reaction; uses `eil_kappa`).
- `eil_vort_lattice_sc=True` with `eil_field_list`: the *true periodic* magnetic-Bloch vortex lattice (formulation A, extreme type-II): a complex order parameter $Psi(bold(r))$ with a real node at every core and the full Abrikosov supercurrent phase, swept over $B\/H_(c 2)$ to give $⟨ N(0) ⟩ (B)$ ($d$-wave $tilde sqrt(B)$ Volovik). `eil_lattice` is `'square'` or `'triangular'`, `eil_nvortex` sets the flux quanta per cell, finite `eil_kappa` adds London screening (and `eil_vort_scA=True` makes $bold(A)$ fully self-consistent from the quasiclassical current $bold(j)_s = ⟨ bold(v)_F "Im" g ⟩$, the `je` `A_renew` scheme).
- `eil_vort_dvector=True`: the self-consistent triplet d-vector vortex/lattice texture (dominant winding + core-localized subdominant; spin-matrix Riccati).
- `eil_gap_orbital`: an orbital-basis pair potential whose *low-energy projection* onto the FS bands sets the gap (Nagai–Nakamura multiband Eilenberger, JPSJ *85*, 074707 (2016), Eq. 43; needs a Wannier FS), superseding `gap_sym`/`delta0`.

The companion driver `calc_vortex_lattice_symmetry` (called from the library) minimizes the Ichioka–Machida lattice free energy over the cell apex angle and the gap-vs-lattice orientation $theta_0$ to determine the *stable vortex-lattice symmetry* and its field evolution (e.g. the $d$-wave triangular → square transition near $H_(c 2)$). When a Wannier FS is supplied, $theta_0$ rigidly rotates the whole crystal (FS + gap), so the Fermi-velocity anisotropy also enters the selection.

== option=27: NMR in the SC State (`NMR_SC`)

Computes the two standard NMR observables in the superconducting state from the same BdG spin susceptibility used by option=12,13: the Knight shift and the spin-lattice relaxation rate,

$ K(T) prop chi'_s (bold(q)=0, omega -> 0), quad 1\/(T_1 T) prop angle.l |A(bold(q))|^2 chi''_s (bold(q), omega_0) angle.r_bold(q) \/ omega_0 . $

Both are also evaluated in the *normal* state at the same temperature, $bold(k)$-mesh, $bold(q)$-mesh and broadening, and it is the ratios that carry the physics: $K_s\/K_n$ is the Yosida function for a singlet gap, and $(1\/T_1T)_s\/(1\/T_1T)_n$ shows the Hebel--Slichter peak of a full gap or the power law of a nodal one ($T^3$ for line nodes). In the ratio the material prefactors ($gamma_n$, the hyperfine scale, $g$) cancel and so does most of the mesh error, so the ratios converge far faster than either side alone — treat the raw `K_sc`, `K_n`, `(1/T1T)` columns as diagnostics and quote the ratios. The normal-state reference is computed with the normal-state bubble, which is the exact $Delta -> 0$ limit of the SC one and 4x cheaper.

The gap is either the `gap_sym`/`delta0` form factor (the same one option=12,13 build) or, with `nmr_gap_file`, a self-consistent RPA/FLEX/Eliashberg gap exported by option=15/23. Whichever it is, it is read as $Delta(0)$: at each temperature only its *amplitude* is rescaled by the BCS interpolation $Delta(T) = Delta(0) tanh(1.74 sqrt(T_c\/T - 1))$, the shape being held fixed. $T_c$ defaults to the weak-coupling estimate $max |Delta|_"FS" \/ 1.764$ taken from the maximum *on the Fermi surface* — the zone-wide maximum can be many times larger when the form factor has little weight on the sheets that carry states, which would put the whole sweep at the wrong reduced temperature — and can be overridden with `nmr_tc`.

The sweep covers `nmr_nt` temperatures over `nmr_trange` (in units of $T_c$) and writes `nmr_sc.dat` and `nmr_sc.png`. Before it starts, the run prints the spin channel in use, $Delta(0)$ and $T_c$, a scale check on the $bold(k)$-mesh/broadening hierarchy (the broadening must satisfy $delta gt.tilde v_F d k$), a cost estimate (the total scales as $N_T N_q N_k$), and — with `nmr_qconv=True` — the change in the $bold(q)$-sum between `nmr_qsize` and half that mesh.

= Color Plot Settings (`color_option`)

For option=0, 2, and 3, each point on the band or Fermi surface can be colored according to a physical quantity.

#table(
  columns: (auto, 1fr),
  [*color_option*], [*Meaning*],
  [0], [No color (black)],
  [1], [Orbital weights specified by `olist` are mapped to RGB (red/green/blue)],
  [2], [Group velocity magnitude $|bold(v)(bold(k))|$ is shown as a color gradient],
  [3], [(option=3, `FERMI_3D`, only) the Eilenberger pairing gap $"Re"[phi(bold(k))]$ — from `gap_sym`/`delta0`, or the Nagai–Nakamura projection of `eil_gap_orbital`/`eil_gap_file` if set],
)

For `color_option=1`, specify orbital indices in `olist` as `[R component, G component, B component]`. To assign multiple orbitals to the same color, use a nested list. Example:

```python
olist = [[0, 3], [1, 4], [2, 5]]
# Orbitals 0 and 3 → Red, orbitals 1 and 4 → Green, orbitals 2 and 5 → Blue
```

= Parameter Reference

This section explains all parameters in the upper part of `main.py`, including their physical meaning.

== Basic Settings

- `fname` (string): Path to the input file or directory, as `ftype` requires (ftype=1: the file name including its extension, ftype=2: the name without `_hr.dat`, ftype=0, 4: a directory).

- `ftype` (integer): Input Hamiltonian file format (see "Input Files").

- `brav` (integer): Bravais lattice type (see "Bravais Lattice Setting").

- `sw_soc` (bool): Switch to enable spin-orbit coupling.

- `inv_tol` (float, default $3 times 10^(-2)$): Relative residual below which the inversion operator is accepted as a symmetry of $H(bold(R))$ (option=14,15,16,23, which reconstruct $Delta(-bold(k))$ from $Delta(bold(k))$ with it). Fitted first-principles models carry a few percent of error, so the default is deliberately loose; raise it if a fitted model is rejected, and check the printed residual before trusting a run that needed a large value.

== Mesh Settings

- `Nx, Ny, Nz` (integer): Number of $bold(k)$-point mesh divisions in the first Brillouin zone along $x$, $y$, $z$. Used for 3D $bold(k)$-space integrations and FFT-based convolutions (FLEX, Eliashberg, chi). For 2D systems, set `Nz=1`. Powers of 2 (32, 64, ...) are preferred for FFT efficiency. Memory consumption scales as $tilde.equiv N_x N_y N_z dot N_w dot N_"orb"^2$.

- `Nw` (integer): Number of Matsubara frequencies for FLEX/Eliashberg/SC-chi calculations (option=12,13,14,15,16,23; not used with `freq_grid='ir'`), or number of real-frequency points for DOS/spectral function calculations (option=1,4, etc.). Matsubara frequencies are $omega_n = (2n+1) pi T$ for $n = 0, 1, \ldots, N_w - 1$. Lower temperatures require larger $N_w$ (typical: 256–1024). The effective frequency cutoff $omega_c = (2 N_w - 1) pi T$ shrinks proportionally to $T$, and the Matsubara convolutions use a sharp cutoff (no tail correction), so a temperature scan at fixed `Nw` accumulates a systematic drift; choose `Nw` so that $omega_c$ safely exceeds the bandwidth and `U` at the lowest temperature (see also the notes under option=14).

- `kmesh` (integer): Number of $bold(k)$ points along the symmetry line for band and spectral function plots. Larger values yield smoother plots (200–500 is typical).

== Lattice Constants

- `abc` (list, unit: Å): Lattice constants $a, b, c$. Used for computing group velocities $v = 1/ hbar (partial E)/(partial bold(k))$ and the physical length scale of symmetry paths. For Wannier90, use values consistent with the WIN/WOUT files.

- `alpha_beta_gamma` (list, unit: degrees): Lattice angles $alpha, beta, gamma$. For orthorhombic and cubic systems, set to `[90., 90., 90.]`.

== Temperature and Chemical Potential

- `tempK` (float, unit: K): Temperature in Kelvin. Converted internally to $T = k_B$ `tempK` (in eV).

- `temp` (float, unit: eV): Directly specifies $k_B T$ in eV. If both `tempK` and `temp` are defined, `temp` takes precedence.

- `fill` (float): Band filling. The chemical potential $mu$ is determined self-consistently from $sum_(bold(k), n) f(epsilon_n (bold(k)) - mu) = N_k dot f"ill"$. Allowed range is 0 to `no` (the number of orbitals in the Hamiltonian); `fill = no` corresponds to full filling. The mapping to the physical electron count per unit cell depends on the Hamiltonian convention:
  - Without SOC (`sw_soc=False`): the spin degree of freedom is not included explicitly, so `fill` represents the electron occupation per spin. For example, a 3-orbital model is full-filled at `fill = 3` (6 electrons including both spins) and half-filled at `fill = 1.5`.
  - With SOC (`sw_soc=True`): the Hamiltonian already includes spin (`no = 2 N_"orb"`), so `fill` directly equals the total electron count per unit cell. A 3-orbital SOC model (`no = 6`) is full-filled at `fill = 6` and half-filled at `fill = 3`.

- `mu0` (float, unit: eV): If defined, this value is used directly as the chemical potential, bypassing the self-consistent calculation from `fill`.

== Energy Range and Broadening

- `Emin, Emax` (float, unit: eV): Lower and upper bounds of the energy range for DOS and spectral function calculations.

- `delta` (float, unit: eV): Broadening parameter $delta$ for spectral calculations. This is the small imaginary part added to the Green's function, which broadens the Dirac delta function into a Lorentzian of finite width. Too large a value smears out features; too small a value introduces numerical noise (a typical range is 0.01–0.05 eV).

- `Ecut` (float, unit: eV): Fixed energy $omega_0$ for $bold(q)$-space susceptibility maps (option=9,11). Set near zero to probe the Fermi surface region.

- `delta0` (float, or a list of length `Norb`, unit: eV): Amplitude of the initial gap function for SC-chi calculations (option=12,13). Physically corresponds to the SC gap size (typical: $10^(-3)$–$10^(-2)$ eV ≈ 1–10 meV).
  - As a float: a single internally-generated gap shape (common to all bands) is scaled to this maximum amplitude. Setting it to 0 reduces to the normal-state calculation.
  - As a list (e.g. `delta0=[0.,0.2,0.3,-0.1,0.]`): enables multi-gap mode, where each band's amplitude and sign are set independently — use this to represent sign-changing gaps such as $s^plus.minus$.

== Transport Parameters

- `tau_const` (float, unit: fs): Constant relaxation time $tau$ used when `tau_mode='const'` (the Fermi-level $tau$ for `'e2'` and `'dos'`). `'epa'` and `'flex'` do not use it, so it need not be positive there, and the library functions also accept `None`. In the constant relaxation-time approximation this is a free parameter typically chosen by comparison with experiment (typical metals: 1–100 fs). Note that $R_H$ carries $tau^2\/tau^2$ and is therefore *independent* of a constant $tau$: with `tau_mode='const'` the Hall coefficient is a pure band-structure number with no temperature dependence beyond the Fermi window.

- `tau_mode` (str, default `'const'`): Relaxation-time model for option=5.
  - `'const'`: constant `tau_const`, $bold(k)$- and band-independent.
  - `'epa'`: $tau(bold(k))$ from the electron--phonon averaged coupling (Samsonidze--Kozinsky EPA), read from `epa_file`. This is what makes $R_H (T)$ and sheet-dependent mobilities meaningful.
  - `'e2'`: Fermi-liquid $tau$, $1\/tau prop (epsilon-mu)^2 + (pi k_B T)^2$ (electron--electron scattering), normalized to `tau_const` on the Fermi level.
  - `'dos'`: $tau prop 1\/D(epsilon)$ (golden-rule scattering rate), used as a ratio to the Fermi level: $tau(epsilon) = tau_"const" D(mu)\/D(epsilon)$.
  - `'flex'`: quasiparticle Boltzmann transport from the FLEX self-energy. Needs `self_en_wannier.npz` from a FLEX run with `sw_out_self=True` at the same temperature. The orbital matrix $Sigma(bold(k), i omega_n)$ is Padé-continued *on the k-mesh FLEX was solved on* (the solver's own data, not an interpolation, so no interpolation noise turns into spurious poles), and $Sigma^R (bold(k), epsilon)$ on a real-energy grid is interpolated multilinearly onto the transport mesh (a convex combination, so it stays causal) and projected on the bands; $Z = [1 - partial "Re" Sigma^R \/ partial epsilon]^(-1)$, $epsilon^* = Z(xi + "Re" Sigma^R (0))$ and $Gamma = -"Im" Sigma^R (epsilon^*)$ give the energies $epsilon^*$, velocities $Z v$ and lifetime $tau^* = hbar \/ (2 Z Gamma)$. Since $Gamma$ is taken at $epsilon^*$, the particle--hole asymmetry of the scattering that drives the Seebeck coefficient is kept. The chemical potential is the FLEX one (`fill` is not used) and `tau_const` is ignored. No current-vertex correction (bubble approximation).
  A summary is printed with the Fermi-surface average $angle.l tau angle.r$, its range over the Fermi window, and $angle.l tau^2 angle.r \/ angle.l tau angle.r^2$ — which is 1 exactly when $tau$ is constant and quantifies how far $tau(bold(k))$ shifts $R_H$ away from the constant-$tau$ value.

- `epa_file` (str or `None`): Path to the `epa.x` (`job='egrid'`) output read when `tau_mode='epa'`. The rate is evaluated as in the EPA authors' BoltzTraP implementation (`epa_calc`): $|g|^2$ constant per energy bin, the final-state DOS $rho(epsilon plus.minus omega)$ interpolated from a fine DOS of the transport k-mesh, states below the midgap assigned to the valence grid and above it to the conduction grid, states beyond a grid taking its outermost bin, and $tau lt.eq 1$ ns. The grid edges are absolute QE energies (VBM/CBM), so the Wannier Hamiltonian must use the same energy reference (a model shifted to put the Fermi level at 0 misaligns the bins).
- `epa_pair_file` (str or `None`, default `None`): Adds the phonon-mediated attraction built from this `epa.x` (`job='egrid'`) output to the spin/charge-fluctuation pairing of the linearized Eliashberg equation (option 15, no SOC). The interaction is a retarded density--density attraction $W(i nu) = -(1\/eta) sum_nu angle.l g_nu^2 angle.r_"FS" 2 omega_nu \/(nu^2 + omega_nu^2)$, local and equal for all orbitals of a site ($U = U' = W$, $J = 0$ in the Kanamori notation). $eta$ is the overlap factor of the Fermi-surface average of its band-basis amplitude, chosen so that the Fermi-surface-averaged isotropic coupling equals the EPA $lambda = 2 N_F sum_nu angle.l g_nu^2 angle.r \/ omega_nu$. Being $bold(q)$-independent it acts on the s-wave (isotropic) component only; the eigenvalue of a sign-changing gap such as d-wave is unchanged. The phonon self-energy ($1 + lambda$ mass renormalization) is not included, so the phonon contribution is somewhat overestimated. The EPA grid edges are absolute energies, so the Hamiltonian must use the same reference; a Fermi surface split between the two EPA grids is flagged (for a metal, let one grid contain the Fermi level).
- `epa_pair_scale` (float, default 1): Factor on that phonon coupling, for scanning the competition between spin fluctuations and phonons.

- `flex_gamma_imp` (float, unit: eV, default 0): Impurity damping $Gamma_"imp"$ added to $-"Im" Sigma^R$ for `tau_mode='flex'` (Matthiessen's rule), i.e. the residual resistivity once the FLEX $Gamma prop T^2$ becomes too small at low temperature.

- `flex_npade` (int, default 40): Largest number of Matsubara points in the Padé continuation for `tau_mode='flex'` (the median of the causal results of orders `flex_npade/2`..`flex_npade` is used). k-points where no order is causal are repaired by interpolating over their acausal energy intervals, and their number is reported as a warning. The FLEX k-mesh is stored in `self_en_wannier.npz`; for older files without it the mesh is inferred from the R extent, so rerun FLEX with `sw_out_self=True`. A $Sigma$ solved with `freq_grid='ir'` is continued from its sampled values (`sigma_ir`), all of them when there are fewer than `flex_npade`; `ir_eps=1e-14` or so is recommended (see "Frequency grid" under option=14).

- `hall_mesh_list` (list or `None`, default `None`): $bold(k)$-mesh convergence sweep for the Hall coefficient in option=5, e.g. `[20,30,40]` or `[(20,20,10),(30,30,15)]`. Each entry is either a single integer (cubic mesh) or an `(Nx,Ny,Nz)` tuple. `None` runs a single mesh `(Nx,Ny,Nz)` and gives no convergence information. Because $sigma^((1))$ can come out with the wrong sign on a coarse mesh, running a sweep at least once for a new material is strongly recommended.

- `sw_tdf` (bool): If `True`, the transport distribution function (TDF) is computed first, and transport coefficients are obtained by energy integration. Relevant when using an energy-dependent relaxation time.

Spin degeneracy in both transport modes (option=5,6) follows `sw_soc`: with `sw_soc=False` the bands are spin-degenerate and a factor $g_s=2$ is applied, with `sw_soc=True` both spins are already in the Hamiltonian and $g_s=1$ is used.

== Orbital and Interaction Parameters

- `olist` (list): Orbital indices for color plotting (`color_option=1`); see "Color Plot Settings".

- `U` (float, unit: eV): On-site Coulomb repulsion (Hubbard $U$). Used in FLEX/RPA calculations. This is a key parameter controlling magnetic and superconducting instabilities.

- `J` (float, unit: eV): On-site Hund's coupling constant. The screened interaction $U' = U - 2J$ (Kanamori screening) is used automatically.

- `chiolist`, `site` (integer arrays, undefined by default): The orbital-pair basis of $chi$ (pairs $(a, b)$ of 1-based orbital indices) and the site label of each pair. If undefined, all orbital pairs are built from `site_prof`. When only some orbitals interact -- d-p models, periodic Anderson models -- listing *only the correlated pairs* gives the RPA, FLEX and linearized Eliashberg results of the full-basis calculation with $U$ on those orbitals alone, in the much smaller basis of $N_"corr"^2$ pairs ($Sigma$ and $Delta$ appear on the correlated orbitals only). With SOC include every spin combination of the correlated orbitals ($4 n^2$ pairs for $n$ orbitals). The printed $chi$ is then the correlated block only, not the total susceptibility including the conduction electrons; and the static part of the FLEX $Sigma$ enters the correlated orbitals alone and effectively shifts $epsilon_d - epsilon_p$, so mind the `sub_sigma` (double-counting) choice.
- `orb_dep` (bool): If `True`, orbital-dependent interaction matrices `Umat` and `Jmat` are used. If `False` (default), the constant values `U` and `J` are applied uniformly to all orbitals. The diagonal of `Umat` is the intra-orbital $U_(i i)$ and its off-diagonal is the inter-orbital $U'_(i j)$ itself (the Kanamori relation $U' = U - 2J$ is *not* re-applied), `Jmat` holds $J_(i j)$; both must be real and $N_"orb" times N_"orb"$, or one site worth of orbitals when `site_prof` lists several equally sized sites -- the per-site matrix is then copied onto every site (e.g. a 5-orbital set on a 2-Fe, 10-orbital model). Only pairs living on the same site enter the vertices, so the inter-site blocks are never read.

- `UJ_material` (string, optional): With `orb_dep=True`, loads a stored parameter set (e.g. cRPA $U$/$J$ of the Fe-based superconductors) instead of writing `Umat`/`Jmat` out by hand. The set is read from `UJ_dir` (default `UJ`), which may be a directory of per-material files `{UJ_dir}/{UJ_material}.toml` or a single TOML file holding several materials as tables (`[FeSe]` ...); each entry stores `U` and `J` either flat (row-major $N_"orb"^2$ values) or as nested rows, and may carry `source`/`doi` (printed on loading) and `orbitals` (orbital names in file order, length-checked). JSON files (`.json`) are still read; a `.toml` of the same name takes precedence. Bundled sets: `FeSe`, `FeTe`, `LiFeAs`, `BaFe2As2`, `LaFeAsO`, `LaFeAsOnakamura` (5 Fe-$d$ orbitals, cRPA values; each file records its reference in `source`/`doi`, which should be cited when the set is used; they cover one site and are replicated over the sites of a multi-site model as described under `orb_dep`). The orbital order of the file must match that of the hamiltonian. `None` (default): use the `Umat`/`Jmat` written in the header.

- `m_diis_num` (integer, optional): DIIS (Pulay-accelerated mixing) history length for FLEX (option=14) and nonlinear Eliashberg (option=23). Values $>= 2$ enable Pulay extrapolation; $1$ falls back to linear mixing. Defaults to 5 if undefined. Larger values speed up convergence at the cost of memory ($N_x N_y N_z dot N_w dot N_"orb"^2$ per slot).

== Initial Gap Function Symmetry (`gap_sym`)

Specifies the symmetry of the initial gap function when solving the Eliashberg equation (option=15,23), or when generating the initial gap shape for the SC-chi calculations (option=12,13).

#table(
  columns: (auto, 1fr),
  [*gap_sym*], [*Symmetry*],
  [0], [$s$-wave (uniform positive sign for all $bold(k)$)],
  [1], [$d_(x^2-y^2)$-wave ($cos k_x - cos k_y$ type)],
  [2], [$s^plus.minus$-wave (sign changes across the nesting vector)],
  [3], [$d_(x y)$-wave ($sin k_x sin k_y$ type)],
  [4], [$d_(x z)$-wave ($sin k_x sin k_z$ type)],
  [5], [$d_(y z)$-wave ($sin k_y sin k_z$ type)],
  [6], [$d + i d$ (chiral, complex: index 1 $+ i dot$ index 3)],
  [7], [$d_(x z) + i d_(y z)$ (chiral, complex: index 4 $+ i dot$ index 5); $|phi|$ vanishes on the whole $k_z=0$ plane, i.e. a *horizontal* line node],
  [-1], [$p_x$-wave],
  [-2], [$p_y$-wave],
  [-3], [$p + i p$ (chiral, complex: $p_x + i p_y$)],
)

The chiral entries ($6$, $7$, $-3$) are complex and *cannot seed the linearized Eliashberg equation* (option=15), whose kernel and seed are real: run the two real partners separately (1 and 3 for $d+i d$, 4 and 5 for $d_(x z)+i d_(y z)$, $-1$ and $-2$ for $p+i p$), verify that they are degenerate and orthogonal, and form $Delta_1 plus.minus i Delta_2$ afterwards — linear theory cannot select the chirality, the condensation energy below $T_c$ does. They are meant for the Eilenberger form factor (option=24,25,26), where they are used directly.

Written as $cos(2 pi k_x)$ etc. these harmonics carry the symmetry they are named after only on a tetragonal / orthogonal single-site lattice: on an fcc, bcc or hexagonal cell $cos(2 pi k_x) - cos(2 pi k_y)$ is periodic but is *not* $B_(1g)$ — it fails to change sign under the Cartesian $C_(4z)$ of the crystal. The form factors are therefore built as *Cartesian lattice harmonics on the R shells of the model*,

$ phi_Gamma (bold(k)) = sum_(bold(R) in "shell") K_Gamma (hat(bold(R))) e^(i bold(k) dot bold(R)), quad bold(R)_"cart" = bold(n) dot "avec", $

with $K_Gamma$ the Cartesian harmonic ($hat(n)_x^2 - hat(n)_y^2$, $hat(n)_x hat(n)_y$, $hat(n)_x$, …) evaluated along each neighbour direction. Note that *no axis decomposition is involved*: the phase $bold(k) dot bold(R) = 2 pi bold(k)_"frac" dot bold(n)$ is basis independent, so only the shell *coefficients* need the Cartesian geometry — `sw_dec_axis` is not needed for the gap, and the $bold(k)$-mesh periodicity every BZ sum and FFT relies on is untouched ($phi(bold(k)+bold(G)) = phi(bold(k))$ holds by construction, since $bold(R)$ is a real lattice vector). On an orthogonal single-site lattice this reproduces the fractional formulas bit for bit, for every entry of the table above and any axis ratio, so existing inputs are unaffected.

`main.py` registers the lattice for this automatically (`plibs.set_harmonic_lattice(avec, rvec)`); `plibs.lattice_harmonic()` builds one harmonic directly, and `gap_symms(..., avec=, rvec=)` takes the lattice explicitly. If the model's $bold(R)$ set is too small to carry the requested symmetry, a warning is printed and the fractional form factor is used instead.

$s^plus.minus$ (`gap_sym=2`) is the one case the shell cannot be fixed by symmetry: every shell carries $A_(1g)$, so $K_Gamma$ cannot discriminate, and the default stands in for it with the length of the *fractional* seed pattern $(1, plus.minus 1, 0)$. That stand-in is basis dependent -- on a centred cell $(1,1,0)$ is not the Cartesian $[110]$, so the cutoff inflates and a far shell is taken ($"KFe"_2"As"_2$: $4a$ for the 2-Fe cell, $2 sqrt(2) a_"Fe"$ for the 1-Fe one; neither changes sign between the $Gamma$ and the zone-corner pocket). Which bond the sign change sits on is a physical choice, so the shell $|bold(R)|$ can be named explicitly in Angstrom with `gap_rcut` (`main.py`; `set_harmonic_lattice(..., s_rcut=)` internally). For the iron pnictides that shell is the crystal axis $a$ (3.825 A for $"KFe"_2"As"_2$): one number gives the same $2 cos(k'_x a_"Fe") cos(k'_y a_"Fe")$ in the 2-Fe cell (in-plane 1st shell) and in the 1-Fe cell (2nd, diagonal shell). `get_initial_gap` prints the $|bold(R)|$ of the shell actually used.

Two caveats remain. The shell of a chiral pair is shared between its two partners (otherwise $|phi|$ is not invariant under the rotation that mixes them — hexagonal $E_(2g)$ fails by 77% with per-partner shells), which the code does automatically. And for a *multi-site* cell the pairing bond is $bold(R) + bold(tau)_j - bold(tau)_i$ rather than $bold(R)$, so the sublattice phase is still missing: for the 2-Fe cell of the iron-based systems, or any model where inversion exchanges sites, use the orbital-basis route (`eil_gap_orbital` / `eil_gap_file` with the Nagai–Nakamura band projection).

== $bold(k)$-Space Settings

- `kz` (float, in reduced coordinates): The $k_z$ value for the 2D Fermi surface plot (option=2). Ranges from 0 to 0.5 ($k_z = 0$: $Gamma$-plane, $k_z = 0.5$: zone boundary plane).

- `kscale` (list or float): Display scale for each axis in the 3D Fermi surface plot (option=3). Example: `kscale=[1.0, 1.0, 0.5]` compresses the $k_z$ direction by half. The Fermi surface is folded into the first BZ (the Wigner-Seitz cell) only at the default `[1.0, 1.0, 1.0]`.

- `k_sets` (list of lists): Coordinates of symmetry line endpoints in reduced units (0 to 1). Define this to specify a custom symmetry path instead of the auto-generated one. Example: `k_sets=[[0,0,0],[0.5,0,0],[0.5,0.5,0]]`

- `xlabel` (list of strings): Labels for the points in `k_sets`. Example: `xlabel=[r'$\Gamma$','X','M']`

- `at_point` (list): Coordinates of the $bold(q)$ point (in reduced units) for the single-$bold(q)$ susceptibility calculation (option=8).

== Switch Variables

- `sw_unit` (bool): If `True` (default), physical constants in SI units are used and output is in physical units. If `False`, a dimensionless system with $hbar = k_upright(B) = e = 1$ is used.

- `sw_omega` (bool): Switch for option=11 — compute the pairing susceptibility on the real frequency axis (`True`) or Matsubara frequency axis (`False`).

- `sw_self` (bool): If `True`, the FLEX self-energy is incorporated into option=4 (spectrum), option=15 (linear Eliashberg), or option=23 (nonlinear Eliashberg) via the renormalized Green's function.

- `sw_out_self` (bool): If `True`, the FLEX self-energy is written to `self_en.npz`, `self_en_wannier.npz` (and `sigma.bin` on the uniform grid), and in option=15, 23 the gap to `gap.npy`, `gap_iw.npy`, `gap_wannier.npz`.

- `sw_in_self` (bool): If `True`, the previous FLEX self-energy is loaded from `sigma.bin` as the initial guess for the iterative self-consistent loop.

- `sigma_in_scale` (float, default 1.0): Factor multiplying the $Sigma$ seed loaded via `sw_in_self=True`. Intended for $U$-annealing; $(U_"new"\/U_"old")^2$ is a reasonable choice (see the notes under option=14).

- `sw_from_file` (bool): If `True` (with `sw_self=True`), the self-energy is read from `self_en.npz` and the FLEX calculation is skipped; the Eliashberg equation is then solved directly with it. A file written on another frequency grid or at another temperature is rejected.

- `sw_check_only` (bool): Used only by option=23 (nonlinear Eliashberg). If `True`, the routine stops right after the internal linearized-Eliashberg solve (reporting the Stoner factor $S$ and eigenvalue $lambda_"eliash"$) without running the nonlinear loop — handy for quickly bracketing $T_c$ via a temperature scan. Regardless of this flag, the nonlinear loop is also skipped automatically whenever $S >= 1$ or $lambda_"eliash" < 1$.

- `sw_rescale_flex` (bool): For FLEX (option=14), uniformly rescale $chi_0$ by $(1-10^(-4))\/S$ whenever the Stoner factor $S$ reaches 1 during the iteration, preventing the magnetic-instability divergence. This stabilizes the early iterations; the $Sigma$ feedback normally deactivates the rescaling by convergence, leaving the converged solution unbiased. If the rescaling is still active at the final iteration a warning (`[FLEX] WARNING`) is printed, meaning the normal state is magnetically unstable (see the notes under option=14).

- `freq_grid` (`'matsubara'` / `'ir'`, default `'matsubara'`), `ir_wmax` (float or None), `ir_eps` (float, default $10^(-8)$): frequency grid of option=14, 15, 16, 23; `'ir'` uses the sparse sampling of the IR basis (see "Frequency grid" under option=14; `ir_eps=1e-14` or so for a Padé continuation).

- `sw_chi0_tail` (bool): For option=14,15,23 (no-SOC path), evaluate $chi_0$ with the tail-corrected convolution (default True; ignored on the IR basis). The Matsubara truncation error improves from $O(1\/N_w)$ to $O(1\/N_w^2)$, so a smaller `Nw` reaches the same accuracy (the $chi_0$ stage costs about 3x; effective for $N_w gt.eq 64$). See the notes under option=14.

- `sw_nosym`, `zeeman_h`, `fflo_q`, `nosym_arnoldi_m`: the symmetry-free linearized Eliashberg equation (see option=15).

- `imp_cell`, `imp_sites`, `imp_conc`, `imp_seed`, `imp_V`: the impurity supercell (see option=21).

- `sw_dec_axis` (bool, *Fermi-surface plots only*): If `True`, the lattice vectors are re-expressed in the conventional orthogonal axes and $bold(R)$ is transformed with them. This is a *drawing* aid: it unfolds a centred cell — the body-centred 122 cell above all — into the familiar box, which is what makes a 10-orbital (2-Fe) Fermi surface directly comparable with a 5-orbital (1-Fe) one. $E(bold(k))$ at a given *physical* $bold(k)$ is preserved exactly (checked to $10^(-14)$ for every `brav`).

  For `brav=1,2,6,7` the transformed $bold(R)$ has half-integer components — that is precisely the unfolding — and the consequence is that $H(bold(k))$ is no longer periodic in the decomposed fractional $bold(k)$. Every BZ sum, FFT convolution, irreducible-wedge folding and cell volume ($det("avec")$, which is then not the volume per formula unit) is therefore invalid with this switch on, so it is accepted only for the plotting modes option=0 (band), option=2 and option=3 (Fermi surface); any other mode stops with an error. The band plot qualifies because it only evaluates $E(bold(k))$ along a path — the DOS (option=1) is a BZ sum and does not. With option=0 an auto-generated symmetry path is rebuilt for the decomposed (simple orthogonal) axes, since the path of the centred cell would mean a different physical $bold(k)$ there; a `k_sets` you supply yourself is used as given and must already be in the decomposed axes. Hexagonal and trigonal cells ($sqrt(3)\/2$ components) cannot be put on an integer orthogonal lattice at all. Gap symmetries no longer need this switch either — they are built from Cartesian $bold(R)$ shells (see `gap_sym` above).

=== Eilenberger Parameters (option=24,25,26)

These drive the quasiclassical Eilenberger solvers. The temperature is the global `tempK`/`temp` and the gap symmetry is the global `gap_sym`.

*Common (all three modes):*

- `eil_coupling` (float): the dimensionless separable pairing coupling $lambda$ (with $⟨ |phi|^2 ⟩_"FS" = 1$). Larger $lambda$ → higher $T_c$.
- `eil_wc` (float, unit: eV): the fixed Matsubara cutoff energy, which sets the pairing scale / $T_c$.
- `eil_fs_kind` (`None`/`'iso'`/`'ellipse'`/`'tb'`/`'cyl'`/`'sphere'`/`'spheroid'`/`'wannier'`): the Fermi surface. `None` = isotropic cylinder (the homogeneous penetration calc falls back to `'ellipse'`); the in-plane model FSs `'iso'`/`'ellipse'`/`'tb'` and the 3D model FSs `'cyl'` (corrugated cylinder), `'sphere'`, `'spheroid'` are built from `eil_fs_params` (the 3D kinds need `eil_fs_nkz` $> 1$); `'wannier'` = the real FS + Fermi velocities of the loaded band (gap symmetry/multiband from `gap_sym`, `delta0`, `eil_gap_orbital`).
- `eil_fs_params` (tuple): model-FS parameters — ellipse masses $(m_x, m_y)$, the `tb` hopping, $(t, t_z)$ for `'cyl'`, or $(m_x, m_y, m_z)$ for `'spheroid'`.
- `eil_fs_nkz` (int, default 1): number of $k_z$ slices stacked into the Fermi surface (`'wannier'` or a 3D model kind). $1$ = a single $k_z=0$ cut (quasi-2D); $> 1$ = a true 3D FS. *Required* ($gt.eq 16$) for a $k_z$-dependent gap: a horizontal line node lives on the $k_z=0$ plane, where a single-slice FS makes the gap identically zero.
- `eil_fs_traj` (`None`/int/tuple): trajectory reduction for the vortex/lattice solvers, whose cost is linear in the number of FS points (a 3D FS carries $10^3$–$10^4$ of them against the 24 directions of the model cylinder). `None` = every point; an integer = the number of $beta$ (direction) bins; `(n_beta, n_v, n_phi)` = full control over direction / $|v_parallel|$ quantile / $phi$ bins. As a reference point, $(48,4,8)$ reduces a 1920-direction FS to 176 with a $times 17.7$ speedup while moving the bulk gap by 0.5% and the vortex-core LDOS peak by 0.2%. Raise it until the answer stops moving. Not used by the surface solver, which needs $k_parallel$ and the band index to find the specular partner.
- `eil_pair_gauge` (`'trs'`/`'soc'`/`'diag'`): the pair-partner gauge of the band gap projection. `'trs'` (default) is for spinless/real hoppings, $phi = u^dagger Delta u$, gauge invariant; `'soc'` uses the time-reversed partner $(i sigma_y) u^*$ in a spinful basis; `'diag'` is the legacy independently-diagonalized $u(-bold(k))$ and is *gauge dependent* (comparison only). The requirement is time-reversal symmetry of the *normal* state $H(bold(k))$, not of $Delta$ — a chiral order parameter is fine (its winding is passed through to $phi$), a time-reversal-broken normal state (ferromagnetic SC, Zeeman inside $H$, Haldane-type hoppings) is not; the run warns when $|angle.l u(-bold(k))| T |u(bold(k)) angle.r| < 1$.
- `eil_spin_order` (`'block'`/`'interleave'`): basis ordering assumed by `eil_pair_gauge='soc'` — `'block'` = [$"orb"_1 dots "orb"_N$ up, $"orb"_1 dots "orb"_N$ down], `'interleave'` = [$"orb"_1$ up, $"orb"_1$ down, $dots$].
- `eil_imp_gamma` (float, unit: eV): the non-magnetic impurity scattering rate $Gamma$ ($0$ = clean).
- `eil_imp_c` (float): the T-matrix $cot delta_0$ — large = Born limit, $0$ = unitary limit.
- `eil_fs_width` (float, unit: eV): the Gaussian Fermi-surface broadening.
- `eil_zeeman` (float, unit: eV): the Zeeman (Maki) field for the LDOS (surface: splits the $d_[110]$ ZEBS into $plus.minus h$; vortex: spin-splits the core states; self-consistent periodic lattice `eil_vort_lattice_sc`: enters the gap equation as self-consistent Pauli depairing).

*Homogeneous (option=24):*

- `eil_method` (`'normalization'`/`'riccati'`): the $(g, f)$ route — `'normalization'` is fast; `'riccati'` matches the inhomogeneous kernel.
- `eil_find_tc` (bool): bisect for $T_c$ at the current impurity setting.
- `eil_imp_sweep` (bool), `eil_imp_list` (array): sweep $Gamma$ over `eil_imp_list` and write $T_c(Gamma)$ to `eilenberger_tc.dat`.
- `eil_pauli`, `eil_spin`, `eil_lambda`, `eil_fs`, `eil_free_energy` (bool): the mutually-exclusive sub-modes described under option=24 above.

*Surface (option=25):*

- `eil_surf_beta` (float, unit: rad): the surface orientation — $0 = [100]$, $pi\/4 approx 0.785 = [110]$ (the $d$-wave ZEBS).
- `eil_surf_dvector` (bool): self-consistent triplet d-vector surface texture.
- `eil_dvec_subratio` (float): the subdominant/dominant coupling ratio for the d-vector texture ($tilde 0.85$ is the bulk threshold).
- `eil_ldos` (bool): also compute the real-frequency surface/core LDOS.

*Vortex / vortex lattice (option=26):*

- `eil_field` (float): $B\/H_(c 2)$ — $0$ = isolated vortex, $>0$ = circular-cell lattice with the Doppler shift.
- `eil_field_list` (list): the $B\/H_(c 2)$ values to sweep on the *true periodic* lattice (e.g. `[0.04,0.08,0.16,0.32]`); `None` = single field.
- `eil_lattice` (`'square'`/`'triangular'`): the periodic-lattice geometry.
- `eil_kappa` (float): the GL parameter $kappa = lambda\/xi$ — large ($gt.eq 10^3$) = extreme type-II (no screening); finite = London screening / Maxwell back-reaction.
- `eil_nvortex` (int): flux quanta per computational cell (supercell).
- `eil_vort_lxi` (float), `eil_vort_ngrid` (int): the isolated-vortex cell half-width (in $xi$) and 2D grid size.
- `eil_vort_field`, `eil_vort_maxwell` (bool): the self-consistent finite-$kappa$ field $B(rho)$ / vector potential $bold(A)(bold(r))$.
- `eil_vort_current` (bool): the circulating supercurrent $j_phi(rho)$.
- `eil_vort_lattice_sc` (bool): the je-style fully self-consistent true periodic lattice; `eil_vort_scA=True` makes $bold(A)$ self-consistent from the quasiclassical current.
- `eil_vort_scA` (bool): with `eil_vort_lattice_sc` and a finite `eil_kappa`, make $bold(A)$ fully self-consistent from the quasiclassical current $bold(j)_s = angle.l bold(v)_F "Im" g angle.r$ (the `je` `A_renew` scheme) instead of using the analytic London $bold(A)$.
- `eil_vort_dvector` (bool): the self-consistent triplet d-vector vortex/lattice texture.
- `eil_vort_chiral` (bool): a self-consistent isolated *chiral* vortex ($p+i p$ for `eil_chiral_ell=1`, $d+i d$ for $2$) via the multi-component complex-amplitude spin-matrix solver — the case the scalar vortex solver refuses, since the core induces the *opposite* chirality and the amplitude is genuinely complex. Both chiralities share one coupling (degenerate partners of the same irrep), and the bulk is the single-channel isotropic gap, not the equal-amplitude (nematic) combination.
- `eil_chiral_ell` (int): the chirality of `eil_vort_chiral` — $1 = p+i p$, $2 = d+i d$.
- `eil_chiral_m` (int): the winding of the *dominant* chiral component ($+1$ parallel to the chirality, $-1$ antiparallel); the induced opposite component follows $m_- = m_+ + 2 ell$.
- `eil_chiral_dvec` (`'z'`/`'x'`): the equal-spin direction of the chiral triplet (`'y'` is proportional to the identity and is not supported by the unitary matrix-Riccati seed).
- `eil_field_dir` (`None` or 3-vector): the field / vortex-line direction for both the isolated vortex and the finite-field lattice, e.g. $(0,0,1)$ = $c$ axis (default), $(1,0,0)$ = in-plane. The vortex lines run along $bold(B)$, so the order parameter varies in the plane perpendicular to it and the problem stays 2D; that plane's two axes are rescaled by their rms Fermi velocities so that a square grid still fits an elliptical core ($xi_1\/xi_2$ is reported). Anything other than $bold(B) parallel c$ needs a 3D Fermi surface (`eil_fs_nkz` $> 1$ or a 3D `eil_fs_kind`). A finite `eil_field` combined with an anisotropic plane is not supported yet (the circular-cell Doppler shift assumes the unscaled plane) and raises an error.
- `eil_vort_tilt` (float, unit: deg): the field tilt from the $c$-axis (quasi-2D: orbital $B_z = B cos theta$, Zeeman $-> h\/cos theta$).
- `eil_gap_orbital` (`None` / $N_"orb" times N_"orb"$ matrix or callable): an orbital-basis pair potential whose low-energy projection onto the FS bands sets the gap (Nagai–Nakamura, JPSJ *85*, 074707 (2016), Eq. 43; needs a Wannier FS), superseding `gap_sym`/`delta0`.
- `eil_gap_file` (`None` / string): the base name (no extension) of a self-consistent RPA/FLEX gap exported as a Wannier-real-space "hopping" file by option=15/23 (`LIN_ELIASHBERG`/`NONLIN_ELIASHBERG`) with `sw_out_self=True` (`output_gap_wannier`, e.g. `'gap_wannier'`). When set, $Delta(bold(R), i omega_n)$ is loaded and used as `eil_gap_orbital` — its inverse Fourier transform $Delta_"orb"(bold(k)) = sum_bold(R) e^(-i 2 pi bold(k) dot bold(R)) Delta(bold(R))$ is projected onto the FS bands (the minus sign is the exact inverse of the `np.fft.ifftn` the exporter uses, verified by round-tripping an export). This is the route to use a *previously computed RPA gap* (e.g. for $"KFe"_2"As"_2$, PRB *84*, 144514) as the vortex pairing form factor. The exporting RPA run and the Eilenberger run must use the *same* Wannier Hamiltonian (same orbital basis, $bold(R)$/$bold(a)$ convention, and ideally $mu$/filling) so that the band eigenvectors and $Delta_"orb"$ share a basis. Supersedes `eil_gap_orbital`/`gap_sym`. The file stores the orbital pair in the physical order $Delta_(a b)$ (the exporter transposes the reversed ctypes view of the Fortran array once, on write, and stamps the npz with `orb_conv`); a gap exported before that stamp existed is still read correctly, with a one-line notice — re-export it to silence the notice. The order matters: $Delta$ is Hermitian but not symmetric, and its antisymmetric part changes sign under the transpose (63% of $max |Delta|$ for the 000AsP RPA solution, enough to make a fully gapped $s_(plus.minus)$ pocket look nodal).
- `eil_gap_iw` (int): the starting Matsubara index for `eil_gap_file` ($0$ = lowest $i omega_0$). The Eilenberger form factor is static; $i omega_0$ carries the symmetry / sign / node / anisotropy structure most sharply and matches the gap usually quoted on the Fermi surface.
- `eil_gap_navg` (int): number of consecutive Matsubara slices averaged for `eil_gap_file` ($1$ = single $i omega_("eil_gap_iw")$ slice). $> 1$ smooths slice noise at the cost of slightly diluting the anisotropy (since $Delta(bold(k), i omega_n)$ gets more isotropic with $n$). The absolute scale is irrelevant — the projected $phi$ is renormalized to $⟨ |phi|^2 ⟩ = 1$.

== NMR Parameters (option=27)

- `nmr_tc` (float or `None`): $T_c$ in eV. `None` uses the weak-coupling estimate $max |Delta|_"FS" \/ 1.764$ from the gap on the Fermi surface.
- `nmr_trange` (list): the reduced temperature range of the sweep, $[T\/T_c "min", T\/T_c "max"]$.
- `nmr_nt` (int): number of temperature points in the sweep.
- `nmr_qsize` (list): the $bold(q)$-mesh of the $1\/T_1$ BZ sum, sub-sampled from `Nx,Ny,Nz` (each entry rounded down to a divisor). The cost is linear in $N_q$ (total $tilde N_T N_q N_k$), so `nmr_qsize=[Nx,Ny,Nz]` is $O(N_k^2)$ and only tractable on toy meshes. The $bold(k)$-sum sets the *resolution* (it needs $delta gt.tilde v_F d k$) while the $bold(q)$-sum only *integrates* a smooth function (error $tilde N_q^(-2)$), so 8–16 per axis is normally plenty — except near a magnetic instability, where $chi_s(bold(q))$ sharpens at $bold(q)_"AF"$. Verify with `nmr_qconv`.
- `nmr_qfold` (bool): fold $bold(q)$ with $-bold(q)$ — exact for a centrosymmetric system with time-reversal symmetry, which the BdG basis already assumes — halving the cost.
- `nmr_qconv` (bool): before the sweep, compare the $bold(q)$-sum on `nmr_qsize` against half that mesh and print the relative change.
- `nmr_w0` (float or `None`, unit: eV): the NMR probe frequency $omega_0$ in $chi''(bold(q), omega_0)\/omega_0$. `None` uses $0.5 delta$; it must satisfy $omega_0 lt.tilde delta << Delta$.
- `nmr_hf_A`, `nmr_hf_B` (float): the Mila--Rice hyperfine form factor $A(bold(q)) = A + 2B(cos q_x + cos q_y)$. $B=0$ gives a flat weight; $A = -4B$ kills the $(pi,pi)$ response.
- `nmr_gap_file` (str or `None`): base name (no extension) of a gap exported by option=15/23 with `sw_out_self=True` (`output_gap_wannier`), used instead of the `gap_sym`/`delta0` form factor. $Delta(bold(R))$ is mesh independent, so the Eliashberg run may use a coarser $bold(k)$-mesh than this sweep. `gap_sym`/`delta0` must still be *valid* (the shared SC-chi block builds a gap from them before this one replaces it) but their values are then irrelevant — and `nmr_spsym` must be set, since `gap_sym` no longer describes the pairing channel. The stored orbital order and the handling of pre-`orb_conv` exports are the same as for `eil_gap_file`.
- `nmr_gap_extrapolate` (bool): fit $Delta(i omega_n) -> Delta(0)$, removing the $O((pi T)^2)$ bias of the lowest Matsubara slice. `False` uses slice `nmr_gap_iw` averaged over `nmr_gap_navg`.
- `nmr_gap_iw`, `nmr_gap_navg` (int): the Matsubara slice and the number of slices averaged when `nmr_gap_extrapolate=False`.
- `nmr_gap_max` (float or `None`, unit: eV): the overall $Delta(0)$ the loaded gap is rescaled to, measured as the maximum *on the Fermi surface*. Required for a linearized-Eliashberg gap (an eigenvector, whose scale is arbitrary); `None` keeps the stored amplitude of a nonlinear/self-consistent gap.
- `nmr_spsym` (bool or `None`): `None` takes the singlet/triplet channel from the sign of `gap_sym`; `True`/`False` overrides it (needed with `nmr_gap_file`). It selects the sign of the anomalous term, i.e. the Yosida-suppressed channel (singlet, or a triplet with $bold(h) parallel bold(d)$) versus the channel preserved at $T=0$ ($bold(h) perp bold(d)$).

= Typical Calculation Workflows

== Checking Band Structure and Fermi Surface

Start by verifying the band structure and Fermi surface.

```python
fname, ftype, brav, sw_soc = 'inputs/SomeMaterial', 2, 0, False
option = 0          # band plot
fill = 2.0          # electron count
abc = [4.0, 4.0, 4.0]
alpha_beta_gamma = [90., 90., 90.]
tempK = 300
```

Then inspect the Fermi surface:

```python
option = 2          # 2D Fermi surface
color_option = 1    # color by orbital weight
olist = [0, 1, 2]   # map each orbital to R/G/B
```

== Superconducting Gap Function (RPA)

```python
option = 15         # linearized Eliashberg equation
Nx, Ny, Nz, Nw = 32, 32, 1, 512
tempK = 50          # compute at relatively low temperature
fill = 2.0
U, J = 1.0, 0.1
gap_sym = 1         # initialize with d_{x^2-y^2}-wave symmetry
sw_self = False     # use RPA (no FLEX self-energy)
sw_out_self = True  # write gap function to file
```

When the calculation completes, the eigenvalue $lambda$ is printed to stdout and the gap function is written to `gap_{ij}.dat` and `gap.npy`.

== FLEX + Linearized Eliashberg

For a more refined calculation incorporating self-energy renormalization, first run option=14 for FLEX, then option=15 with `sw_self=True` and `sw_from_file=True`.

```python
# Step 1: FLEX self-energy calculation
option = 14
sw_out_self = True
m_diis_num = 5

# Step 2: Linearized Eliashberg with FLEX self-energy
option = 15
sw_self = True
sw_from_file = True
```

== Nonlinear Eliashberg (Self-Consistent SC Loop)

To grow $Delta$ to a finite amplitude below the temperature where the linearized Eliashberg eigenvalue reaches $lambda approx 1$, use option=23. The initial $Delta$ (symmetry shape and BCS-scaled amplitude) is generated automatically inside the solver, so there is no need to run option=15 first to produce a `gap.npy` file.

```python
# First, use sw_check_only=True to locate Tc without running the nonlinear loop
option = 23
gap_sym = 1
sw_self = True            # include FLEX self-energy
sw_from_file = True
sw_check_only = True
```

The Stoner factor $S$ and eigenvalue $lambda_"eliash"$ are printed to stdout. Once you have lowered the temperature enough that $lambda_"eliash" > 1$, run the nonlinear loop with `sw_check_only=False`.

```python
option = 23
sw_self = True            # include FLEX self-energy
sw_from_file = True
sw_check_only = False
m_diis_num = 5            # DIIS Pulay acceleration + amplitude-direction Newton acceleration
```

At low temperatures ($T tilde.equiv T_c \/ 5$), increasing the DIIS history to 5–10 typically accelerates convergence.

== Calculations on the IR Basis

At low temperatures, or for models whose bands reach high energies (an MLO basis, for example), the IR basis is accurate with far fewer frequencies than uniform Matsubara frequencies. Add to the FLEX + linearized Eliashberg flow above:

```python
freq_grid = 'ir'     # Nw is not used (the basis sets the number of frequencies)
ir_wmax = None       # None: 2 x max|e_k - mu|
ir_eps = 1.0e-8      # about 1e-14 for a Pade continuation to the real axis
```

== Dynamic Spin Susceptibility in the SC State

To probe the SC-gap dependence of spin excitations, compute $chi_s^"SC"(bold(q), omega)$ with a finite gap:

```python
option = 12               # SC chi spectrum along symmetry line
delta0 = 1.e-2            # initial gap amplitude (eV) ≈ 10 meV
gap_sym = 1               # d_{x^2-y^2}-wave
U, J = 0.8, 0.1
tempK = 50
```

Setting `option = 13` instead computes $chi_s^"SC"(omega)$ at the single $bold(q)$ point given by `at_point`.

= AiiDA Interface

`examples/aiida/pyrpa_aiida.py` runs pyrpa calculations as AiiDA `calcfunction`s / `workfunction`s, so the inputs and outputs of every step and the links between them (provenance) are recorded in the AiiDA database. The calculation itself runs in the same process through `libs.plibs` / `libs.flibs`, so no computer, code or daemon has to be configured. A `CalcJob` that submits `main.py` to a cluster is not implemented.

== Installation

```bash
pip install aiida-core            # tested with 2.9.4
```

It brings about 60 dependencies and may change the versions of packages such as `click` or `ipython`. To keep the working environment untouched, install it into a virtual environment that still sees the MKL builds of numpy/scipy:

```bash
python -m venv ~/venvs/pyrpa-aiida --system-site-packages
source ~/venvs/pyrpa-aiida/bin/activate
pip install aiida-core
```

The tests and `run_example.py` use a temporary profile (`SqliteTempBackend`) by default, so PostgreSQL, RabbitMQ and `verdi presto` are not needed (nothing is kept afterwards). To keep the results, create a profile (e.g. `verdi presto`) and pass it with `run_example.py --profile NAME`.

== Steps

#table(
  columns: (32%, 1fr, 26%),
  align: (left, left, left),
  table.header([*Step*], [*What it does*], [*Main outputs*]),
  [`load_hoppings`], [read the model (`ftype` as in `main.py`; for ftype=4 the directory is passed as `FolderData`)], [`rvec`, `ham_r` (, `S_r`)],
  [`chemical_potential`], [chemical potential of the bare bands (full $bold(k)$ mesh)], [`mu`],
  [`flex_self_energy`], [FLEX self-energy (`mkself`, with SOC `mkself_soc`)], [`sigma`, `mu_self`, `flex_info`],
  [`rpa_susceptibility`], [$chi_s$, $chi_c$ at $nu = 0$ (FLEX when `sigma` is given)], [`chi`, `stoner`],
  [`linearized_eliashberg`], [linearized Eliashberg equation (RPA / FLEX)], [`lambda_eliash`, `gap`, `gap_wannier`],
  [`flex_eliashberg`], [FLEX → $chi$ → $lambda$ in one step, $Sigma$ not stored], [the results of the three above],
  [`eilenberger_bulk`], [homogeneous Eilenberger equation (`gap_wannier` can be the form factor)], [`result`, `f_avg`],
  [`eilenberger_surface`, \ `eilenberger_vortex`], [surface and vortex (model Fermi surface)], [`profile`, `files`],
)

The workflows `flex_workflow` (model → $mu$ → FLEX), `eliashberg_workflow` (model → $mu$ → [FLEX] → $chi$, $lambda$) and `eliashberg_eilenberger_workflow` (then `eilenberger_bulk`) chain them. Every RPA/FLEX step also outputs `runtime` (OpenMP threads and wait policy, the `libfmod` build, the numpy version, the git commit).

Arrays keep pyrpa's layout ($Sigma$, $Delta$ as `[Norb, Norb, Nw, Nk]` on the irreducible $bold(k)$ mesh). An MLO basis (ftype=3,4) is normalized on reading as in `main.py`, and every later step uses $G^(-1) = (i omega + mu) S - H - Sigma$.

== Options (`options`)

The RPA/FLEX steps and the workflows take an optional `options` (`orm.Dict`) holding what `main.py` sets in its header; the defaults are in `OPTION_DEFAULTS`, and an unknown key is an error.

#table(
  columns: (24%, 14%, 1fr),
  align: (left, left, left),
  table.header([*Key*], [*Default*], [*Meaning (`main.py` counterpart)*]),
  [`soc`], [`False`], [spin-orbit coupling (`sw_soc`): the model is spinful, and U, J, `Umat` and `site_prof` refer to its spin-orbitals],
  [`slist`, `invs`], [`None`], [spin labels ($plus.minus 1$) and 1-based opposite-spin partners; `None`: `main.py`'s order, all up then all down],
  [`site_prof`], [`None`], [orbitals per site (e.g. `[5, 5]`); `None`: one site holding every orbital],
  [`chiolist`, `site`], [`None`], [orbital-pair basis of $chi$ (1-based) and its sites; restricting it to the d--d pairs, e.g., speeds up a U on the d orbitals only],
  [`Umat`, `Jmat`], [`None`], [orbital-dependent U, J (`orb_dep`); one site's worth is copied onto every site (`plibs.expand_UJ_sites`)],
  [`freq_grid`], [`'matsubara'`], [`'ir'`: the IR basis (`freq_grid`); the number of frequencies is set by the basis],
  [`ir_wmax`, `ir_eps`], [`None`, `1e-8`], [cutoff of the IR basis (`None`: $2 max|epsilon_bold(k) - mu|$) and its accuracy],
  [`sw_tail`], [`False`], [tail-corrected $chi_0$ (`sw_chi0_tail`); no SOC, uniform grid only],
  [`eps`, `pp`, `m_diis`, \ `scf_loop`, `sub_sigma`, \ `sw_rescale`], [`1e-4`, `0.3`, `5`, \ `300`, `1`, \ `False`], [FLEX loop controls (`flibs.mkself`)],
  [`require_converged`], [`True`], [fail the step when FLEX does not converge or ends with a Stoner factor $>= 1$],
  [`store_gap`], [`True`], [`False`: do not store $Delta$ on all frequencies, only `gap_wannier`],
  [`store_sigma`], [`True`], [`False`: run the FLEX `eliashberg_workflow` as the single step `flex_eliashberg`, without storing $Sigma$],
)

- *FLEX convergence (`flex_info`).* The Fortran FLEX loop does not return whether it converged, so its log (`iter=`, `esterr=`) gives the iteration count, the final esterr and the verdict. A loop that ends with a Stoner factor $>= 1$ has settled on a magnetically unstable, unphysical fixed point -- near the instability the first overshoots of the iteration can throw it there -- and does not count as converged.
- *The IR grid is passed on.* A FLEX $Sigma$ records the IR basis it was solved on (`ir_wmax`, `ir_eps`), and the steps built on it use the same basis (the cutoff is not re-derived from $mu_"self"$).
- *SOC.* $chi$ holds the spin channels `chiszz`, `chispm`, `chisxx`, `chisyy` and `chic` (`chis_q` is the zz leading eigenvalue). No Stoner factor is output (the SOC $chi_0$ has none, and `main.py` does not report one). The Eliashberg step needs a diagonal parity, and `sw_tail` is not available. Using an SOC gap as the form factor of `eilenberger_bulk` is not supported.
- *MLO and the frequency cutoff.* MLO bands reach high energies (about $+30$ eV for Cu), so on the uniform grid use `sw_tail=True` or a large enough `Nw`, or `freq_grid='ir'`. Without the tail correction at $N_w = 128$ (a 40 eV window), $lambda$ was off by several times in one example.

== Examples

`examples/aiida/run_example.py` runs the main steps:

```bash
python examples/aiida/run_example.py --calc flex                      # square.hop: model -> mu -> FLEX
python examples/aiida/run_example.py --calc eliashberg --flex         # FLEX + linearized Eliashberg
python examples/aiida/run_example.py --calc eliashberg --flex --soc   # SOC (square_2orb_soc.hop, fill = 2)
python examples/aiida/run_example.py --calc eliashberg --flex --mlo-demo   # two-orbital MLO model, IR basis
python examples/aiida/run_example.py --calc eliashberg --hop inputs/Cu --ftype 4 --fill 5.5 \
    --options '{"freq_grid": "ir", "require_converged": true}'   # directory of the ecalj job_mlo output
```

`--options` takes the keys of the table above as JSON. From Python:

```python
import aiida
from aiida import orm
from aiida.storage.sqlite_temp import SqliteTempBackend
aiida.load_profile(SqliteTempBackend.create_profile('pyrpa'), allow_switch=True)
import pyrpa_aiida as pa          # with examples/aiida on sys.path

res = pa.eliashberg_workflow(
    hop_file=orm.SinglefileData('inputs/square.hop'), ftype=orm.Int(1), fill=orm.Float(0.45),
    temp=orm.Float(0.02), U=orm.Float(0.3), J=orm.Float(0.0), mesh=orm.List([16, 16, 1]),
    Nw=orm.Int(128), gap_sym=orm.Int(1), use_flex=orm.Bool(True),
    options=orm.Dict({'freq_grid': 'ir', 'store_sigma': False}))
print(res['lambda_eliash'].value, res['flex_info'].get_dict())
```

= Test Suite

`tests/` holds regression tests of the numerical kernels and physics benchmarks: 17 files, 293 tests. Each test file can also be run directly as a Python script where `pytest` is not installed.

The Fortran shared library `libs/libfmod.so` must be compiled first; if it is not, run `make FC=<compiler> SL=<library>` in the `libs` directory.

== How to Run

Individually:

```bash
python tests/test_transport.py
python tests/test_rpa_flex.py
```

All at once with a summary table (no `pytest` needed):

```bash
python tests/run_all.py            # every file
python tests/run_all.py rpa        # only files whose name contains "rpa"
python tests/run_all.py --runslow  # also the tests marked slow
```

Where `pytest` is available, `pytest tests` works too. Tests that take a `tmp_path` get a temporary directory from the standalone runner as well.

== Coverage by File

#table(
  columns: (auto, auto, auto, 1fr),
  align: (left, right, right, left),
  table.header([*File*], [*Tests*], [*Time*], [*What is checked*]),
  [`test_eliash_nosym.py`], [10], [158 s],
    [the symmetry-free linearized Eliashberg equation (Zeeman field, FFLO), consistency with the SOC solver (bare term $-V\/2$, triplet seeds and parity projection, band-basis output), multi-orbital triplets against the exact kernel, branch tracking under SOC],
  [`test_eilenberger.py`], [59], [143 s],
    [Quasiclassical Eilenberger/Riccati solvers: homogeneous limits, Fortran kernels against Python references, surface Andreev states, vortices and self-consistent vortex lattices, condensation free energy, model and Wannier Fermi surfaces, Pauli limiting, triplet $d$-vector textures],
  [`test_freq_grid.py`], [13], [47 s],
    [the IR basis: reproduction by the basis, the exact $chi_0$, FLEX, linearized/nonlinear Eliashberg, SOC, `sw_nosym` and the phonon vertex against the $N_w -> infinity$ limit of the uniform grid, the tail correction being ignored on IR, Padé from the sampling points],
  [`test_nmr.py`], [22], [27 s],
    [NMR in the SC state: BCS gap interpolation, $bold(q)$-mesh folding, loading self-consistent Eliashberg gaps, Knight shift (Yosida and triplet), Hebel--Slichter peak, $T^3$ law of line nodes],
  [`test_rpa_flex.py`], [57], [18 s],
    [Building blocks of RPA / FLEX / Eliashberg: interaction vertices, RPA matrix algebra, $chi_0$ and $phi_0$ against exact Lindhard sums and pair bubbles, tail-corrected $chi_0$, `mkself` / `mkself_soc` and the double-counting cancellation, the Nambu solution of the nonlinear Eliashberg loop, self-energy regridding, the TRS irreducible $bold(k)$ mesh],
  [`test_aiida.py`], [9], [16 s],
    [The AiiDA interface (`examples/aiida`): every step and workflow matches the direct `plibs` / `flibs` calls (RPA / FLEX / Eliashberg, MLO, SOC, Eilenberger), the `options` (IR basis, `Umat`, `store_sigma`, an unconverged FLEX loop) and the recorded provenance. Skipped without `aiida-core`],
  [`test_nonorth_flex.py`], [8], [5 s],
    [Green's functions, FLEX and Eliashberg in a non-orthogonal (MLO) basis: the matrix Dyson equation, covariance under an orbital rescaling, agreement with the SOC solver, IR, reloading the self-energy, second-order convergence of the tail correction, orbital normalization],
  [`test_effmass.py`], [30], [3 s],
    [Inverse effective masses in orthogonal and MLO bases (interband term, overlap term, degeneracies), and dHvA extremal orbits: section geometry, rejection of open orbits, zone reduction],
  [`test_transport.py`], [32], [1 s],
    [Kubo$arrow.l.r$Boltzmann agreement in the DC limit, Wiedemann--Franz law, symmetrized interband heat vertex, band degeneracy, $bold(k)$-point weights, weak-field Hall response],
  [`test_bz.py`], [6], [< 1 s],
    [first-Brillouin-zone geometry of the 3D Fermi surface (faces, folding)],
  [`test_color.py`], [7], [< 1 s],
    [reflectivity to colour (CIE XYZ, D65, sRGB)],
  [`test_optical.py`], [9], [< 1 s],
    [optical constants and non-metal colours (refractive index, reflectivity, absorption, transmission, Kubelka--Munk)],
  [`test_grid.py`], [6], [< 1 s],
    [Fortran grid kernels of the vortex/vortex-lattice solvers and the Riccati chord kernel],
  [`test_impurity.py`], [7], [< 1 s],
    [impurity supercell and the unfolded spectrum (sign of the phase, sum rules)],
  [`test_phonon_pairing.py`], [3], [< 1 s],
    [phonon-mediated pairing from EPA and its entry into the linearized Eliashberg equation],
  [`test_velocity_mlo.py`], [7], [< 1 s],
    [Band velocities in a non-orthogonal (MLO) basis: the $-epsilon_n C^dagger (partial S\/partial k) C$ overlap term, hermiticity, exact dispatch to the orthogonal path],
  [`test_tools.py`], [8], [< 1 s],
    [the shared references the other files check against (Fermi/Bose identities, closed-form models, the exact $chi_0$), CPU-dependent library selection (AMD, Intel), the runtime environment of gfortran builds (`OMP_WAIT_POLICY`, stdout buffering)],
)

The whole suite takes about 7 minutes (ifx build, 24 threads), most of it in `test_eliash_nosym.py` and `test_eilenberger.py` (the symmetry-free solver works on the full Brillouin zone, and the quasiclassical solvers are self-consistent field problems on 2D grids). Every other file finishes within a minute. `test_aiida.py` is skipped without `aiida-core`, `test_freq_grid.py` without `sparse-ir`.

= Troubleshooting

- *Chemical potential does not converge*: Increase `Nx, Ny, Nz`, or slightly raise `tempK`. A coarse $bold(k)$-mesh can cause instability in the self-consistent $mu$ search.

- *Too much noise in the spectrum*: Increase `delta` slightly, or increase `Nw` to improve the frequency resolution.

- *Band structure looks wrong*: Check that `brav` is set correctly. In particular, when interfacing with QE, pay careful attention to the conventions for FCC (`brav=1`) and BCC (`brav=2`).

- *Linear Eliashberg does not converge*: Try reducing `U`, increasing `tempK`, or using a finer mesh (`Nx, Ny, Nz`). If the eigenvalue $lambda$ exceeds 1, the system is inside the superconducting phase at that temperature.

- *Nonlinear Eliashberg diverges*: Increase `m_diis_num` to 5–10, raise `tempK`, or increase `Nw`. Running with `sw_check_only=True` first to inspect the Stoner factor $S$ and eigenvalue $lambda_"eliash"$ is also useful. If $S$ exceeds 1, the system is magnetically unstable and `U` must be reduced; if $lambda_"eliash" < 1$, you are still above $T_c$ and need to lower the temperature (in either case the nonlinear loop is skipped automatically).

- *FLEX self-energy diverges*: Set `sw_rescale_flex=True` or reduce `U`.

- *`[FLEX] WARNING: Stoner factor >= 1 at the final iteration`*: the normal state is magnetically unstable at this temperature and interaction, and the $Sigma$ obtained is not a physical solution. Near the instability the first overshoots of the iteration can throw it onto an unphysical fixed point, and it may "converge" to different values for a different mixing rate or frequency grid. Lower `U` or raise the temperature.

- *A gfortran build is very slow*: two OpenMP runtimes compete (see "Building the Fortran kernel"). `_loader.py` sets `OMP_WAIT_POLICY=PASSIVE` automatically, but a value set in the environment takes precedence.

- *Huge Stoner factor or $lambda$ for an MLO model / gap or $Sigma$ disturbed at the frequency edge*: MLO bands reach high energies (about $+30$ eV for Cu), and when the uniform window $(2 N_w - 1) pi T$ is too small the truncation error of $chi_0$ grows large. Keep `sw_chi0_tail=True` and raise `Nw`, or use `freq_grid='ir'`. For a dp model with U on the d orbitals only, zero the s, p rows and columns of `Umat`/`Jmat` (`orb_dep`) and restrict `chiolist` to the d--d pairs, which is also faster.

- *Transport coefficients have wrong units*: Make sure `sw_unit=True`. Setting it to `False` switches to a dimensionless unit system.

- *The Hall coefficient changes sign with the $bold(k)$-mesh*: This is the normal failure mode of $sigma^((1))$, not a bug. Run a sweep with `hall_mesh_list` and read the two diagnostics: if $N_"eff"$ is below $10^3$ the Fermi window is carried by too few states, and if the cancellation ratio is below $10^(-2)$ the system is near compensation and needs a much finer mesh than $sigma_(x x)$ does. Raising `tempK` also widens the Fermi window.

- *The Hall carrier density disagrees with experiment*: Check it against the Luttinger volumes reported in the same output. If $n_e \/ n_h approx 1$ the metal is compensated and $n_H$ is not a carrier density at all — compare mobilities or fit a two-fluid model instead. Also verify that `sw_soc` is set correctly, since it selects the spin degeneracy factor, and that $omega_c tau$ times your experimental field is still $lt.tilde 0.1$.

- *$R_H$ shows no temperature dependence*: With `tau_mode='const'` this is expected — $R_H$ carries $tau^2\/tau^2$ and a constant $tau$ cancels exactly. Use `tau_mode='epa'` with an `epa_file` to get a $bold(k)$-dependent $tau$.
