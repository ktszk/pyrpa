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
  #text(size: 11pt)[2025]
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
- Self-energy calculation via FLEX (Fluctuation EXchange approximation)
- Linearized Eliashberg equation for the superconducting eigenvalue and gap function
- Nonlinear Eliashberg equation for self-consistent SC gap functions
- Band filling, cyclotron mass, and dHvA frequency calculations
- Quasiclassical Eilenberger/Riccati solvers: $T_c$, penetration depth, surface Andreev bound states, vortices and vortex lattices
- NMR in the superconducting state: Knight shift $K(T)$ and $1\/(T_1 T)$ from the BdG spin susceptibility
- Spectral function with impurities and via CPA

All settings in pyrpa are controlled by modifying variables at the top of `main.py`, above the library import lines.

= Notes on This Revision

Several changes alter results produced by earlier versions. The first two invalidate previously computed data outright, and the gap-harmonic change alters every result computed on a centred or hexagonal cell; the rest change what is reported or what is accepted. If you have output from before this revision, re-check it against the notes below.

- *Irreducible $bold(k)$-mesh under time reversal.* `generate_irr_kpoint_inv` matched $bold(k)$-point coordinates with a tolerance of $1\/max(N_x,N_y,N_z)$ — exactly one mesh spacing — and silently matched a *neighbouring* point whenever $1\/N$ was not a binary fraction. Meshes with $N$ a power of two ($4,8,16,32,64$) were unaffected; others were not, e.g. $12^3$ mismapped 919 of its 1728 points. The mapping is now built from integer mesh indices and is exact for every mesh. Any FLEX / Eliashberg / SC-$chi$ result obtained on a non-power-of-two mesh should be regarded as unreliable and recomputed. The same routine was $O(N_"kall"^2)$ and is now $O(N_"kall")$ (9.2 s $arrow$ 0.009 s at $64^3$).

- *Irreducible wedge for $N_x$ odd with $N_y$ even.* That parity combination was an unimplemented branch in `gen_irr_k`: the wedge was left unfilled and the routine wrote past the end of `klist`, crashing for $N_z >= 3$. It is now implemented, and the wedge has been verified against `TRS_irr.typ` for all 384 mesh/parity combinations up to $8 times 8 times 6$.

- *Spin degeneracy in transport.* option=5 and option=6 now take the spin factor from `sw_soc`. Previously $g_s=2$ was applied unconditionally, so results for a Hamiltonian that already carried both spins were a factor of two off.

- *Hall coefficient.* $R_H$ is now formed from the full conductivity tensors rather than from $sigma_(x y)\/(sigma_(x x) sigma_(y y))$, and the Hall carrier density is reported as a positive number with its carrier type. See option=5.

- *Gap form factors are Cartesian lattice harmonics.* `gap_sym` used to be evaluated as a function of the *fractional* coordinates ($cos(2 pi k_x) - cos(2 pi k_y)$ and so on), which carries the labelled symmetry only on a tetragonal / orthogonal single-site lattice: on an fcc, bcc or hexagonal cell that function is periodic but is not $B_(1g)$ at all — it does not change sign under the crystal's Cartesian $C_(4z)$ (checked: the violation is $O(1)$, and the hexagonal $E_(2g)$ doublet norm is off by 100%). The form factors are now built as $sum_bold(R) K_Gamma (hat(bold(R))) e^(i bold(k) dot bold(R))$ over the $bold(R)$ shells of the model, with $K_Gamma$ evaluated in Cartesian space; the labelled irrep and the $bold(k)$-space periodicity now both hold in any Bravais lattice (to $10^(-15)$). On an orthogonal single-site lattice the new harmonics reproduce the old ones bit for bit, for every `gap_sym` and any axis ratio, so results there are unchanged. Anything computed on a centred or hexagonal cell with a non-trivial `gap_sym` — SC-$chi$ (option=12,13), NMR (option=27), the Eilenberger form factor (option=24,25,26) and the Eliashberg seed (option=15,23) — changes, and the new value is the correct one. See the `gap_sym` section for the construction and its remaining multi-site caveat.

- *`sw_dec_axis` is restricted to the plotting modes.* The axis decomposition unfolds a centred cell into the conventional orthogonal axes, which is what makes a 10-orbital (2-Fe) Fermi surface or band structure comparable with a 5-orbital (1-Fe) one. It leaves $bold(R)$ with half-integer components, however, so $H(bold(k))$ is no longer periodic in the decomposed fractional $bold(k)$: BZ sums, FFT convolutions, irreducible-wedge folding and $det("avec")$ as a cell volume are all invalid. It is now accepted for option=0, 2 and 3 only, and any other mode stops with an error instead of returning a quietly wrong number. The gap symmetries that used to motivate it no longer need it (see the item above).

= Requirements

The following Python packages are required:

- `numpy`
- `scipy`
- `matplotlib`
- `skimage` (used for Fermi surface plots option=2,3, cyclotron mass option=18, and dHvA orbits option=19)

The Fortran kernel `libs/libfmod.so` must be compiled beforehand; `libs/flibs` and `libs/plibs` are pure-Python packages that wrap it. The FLEX, linear/nonlinear Eliashberg, and SC chi calculations are executed via OpenMP-parallel Fortran routines.

== Building the Fortran kernel

Run one of these in `libs`:

#table(
  columns: (auto, 1fr),
  [`make auto`], [*Recommended on a cluster.* Detects everything: the toolchain from what is installed and loaded, and the CPU targets from the compute nodes themselves. Builds one `libfmod_<target>.so` per distinct CPU type found, and `libs/flibs/_loader.py` picks the best one the running node supports at import.],
  [`make`], [One library for one target: `ARCH` (default `znver3`, the portable AVX2 baseline), `FC` and `SL` as detected or as given on the command line.],
  [`make dispatch`], [The fixed `znver3` + `znver4` pair, i.e. `make auto` without the detection step.],
)

Explicit settings always win: `make FC=ifx SL=MKL ARCH=znver4`, or `make auto ARCH_LIST="znver3 znver4"`.

*What is detected, and where from.* `FC` is the first of `ifx`, `ifort`, `flang`, `amdflang`, `gfortran` on `PATH`; `SL` follows the modules actually loaded (`MKLROOT` $arrow$ MKL, `AOCL_ROOT` $arrow$ AOCL, otherwise OpenBLAS) — selecting a library whose root is not set is what used to yield a `libfmod.so` with unresolved BLAS symbols that failed only at import. The *CPU targets*, however, must not come from the build host: `-march=native` on a Zen 4 login node produces AVX-512 code that dies with SIGILL on a Zen 3 compute node. `libs/detect_arch.sh` therefore asks each node in the scheduler's inventory (Slurm `sinfo`, PBS `pbsnodes`, else just this machine) what `gcc -march=native` calls its CPU, and prints the distinct targets:

```
$ sh libs/detect_arch.sh --table
NODE         CPU                                      -march
salmon       AMD Ryzen 7 5700G with Radeon Graphics   znver3
unagi        AMD Ryzen 9 7900 12-Core Processor       znver4
yakisoba     AMD Ryzen 9 7900X 12-Core Processor      znver4
```

The result is cached in `libs/.arch_cache`; `sh libs/detect_arch.sh --refresh` re-probes, `--local` limits it to the current machine. A node that cannot be probed does not silently drop out: it contributes the portable `znver3` default and a warning. Setting `Features=znver3` / `znver4` in `slurm.conf` makes the probe a plain `sinfo` query with no job submission at all — worth doing if you administer the cluster.

At run time nothing has to be selected by hand: `_loader.py` reads the CPU flags and loads the highest-ISA `libfmod_*.so` the machine can actually execute, falling back to a plain `libfmod.so`. For MKL on AMD also export `MKL_ENABLE_INSTRUCTIONS=AVX512` (Zen 4) or `AVX2` (Zen 3), and set `OMP_NUM_THREADS` from the scheduler (`$SLURM_CPUS_PER_TASK`) rather than hard-coding it.

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
  [1], [A single `{fname}.input` file (custom format)],
  [2], [`{fname}_hr.dat` file (Wannier90 default hopping file). *This is the standard choice when interfacing with first-principles codes.*],
  [3], [Non-orthogonal basis hopping in MLO (Maximally Localized Orbital) format],
  [4], [Directory `{fname}/` containing `HamRsMLO` and `HamiltonianPMTInfo` (binary output of ecalj `job_mlo`). Same non-orthogonal MLO basis as ftype=3; the overlap $S(R)$ is read alongside $H(R)$],
  [Other], [`Hopping.dat` file (ecalj hopping file)],
)

*Example for Wannier90 users:*

```python
fname = 'inputs/Sr2RuO4'   # path without the _hr.dat suffix
ftype = 2
```

== Spin-Orbit Coupling (`sw_soc`)

Setting `sw_soc = True` enables spin-orbit coupling (SOC). For SOC-active systems, the Hamiltonian must include spin degrees of freedom (i.e., the orbital dimension is doubled). FLEX (option=14) and the linearized Eliashberg equation (option=15) have SOC-aware implementations (`calc_flex_soc` / `calc_lin_eliash_soc` are dispatched internally). The nonlinear Eliashberg solver (option=23) and SC chi (option=12,13) currently do not support SOC.

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

#block(stroke: 1pt + red, inset: 8pt, radius: 4pt)[
  *Notice: Unverified modes*

  The following modes are tagged *(not implemented)* in the `CalcMode` definition in `main.py`, meaning that the correctness of their implementation has not been sufficiently verified. The physical validity and numerical accuracy of their output are *not guaranteed*; use them *at your own risk*. Independent cross-checks are strongly recommended before citing any results in publications.

  - option=20 (`ELECTRON_MASS`): electron mass along the symmetry line
  - option=21 (`SPECTRUM_IMPURITY`): spectral function with impurities
  - option=23 (`NONLIN_ELIASHBERG`): nonlinear Eliashberg equation
]

== option=0: Band Plot (`BAND`)

Calculates and plots the band dispersion $E_n(bold(k))$ along a symmetry line. The symmetry path is specified via `k_sets` and `xlabel`, or auto-generated from the `brav` setting if not defined.

== option=1: Density of States (`DOS`)

$ "DOS"(omega) = - 1/pi sum_(bold(k), n) "Im" G^0_n (bold(k), omega + i delta) $

Plots both the total DOS and orbital-projected partial DOS.

== option=2: 2D Fermi Surface (`FERMI_2D`)

Draws the Fermi surface in the $k_x$-$k_y$ plane at a specified $k_z$ slice (default: $k_z = 0$, adjustable via `kz`). The chemical potential is determined self-consistently from the electron count `fill`, and the contour $E_n(bold(k)) = mu$ is extracted via skimage's `find_contours`. A rotation matrix `RotMat` can optionally be used to rotate the Fermi surface.

== option=3: 3D Fermi Surface (`FERMI_3D`)

Renders the three-dimensional Fermi surface as a polygon mesh using the Marching Cubes algorithm. The display scale along each axis can be adjusted with `kscale`. With `color_option=ColorMode.GAP` (3), the surface is colored by $"Re"[phi(bold(k))]$, the *same* pairing form factor that drives the Eilenberger calculations (`gap_sym`/`delta0`, or `eil_gap_orbital`/`eil_gap_file` if set) — a quick way to check the actual gap (sign, nodes, anisotropy) on the real 3D Wannier Fermi surface (all sheets/$k_z$, not just the fixed-$k_z$ cut used by the vortex/surface solvers).

== option=4: Spectral Function (`SPECTRUM`)

Computes the electron spectral function from the imaginary part of the trace of the single-particle Green's function,

$ A(bold(k), omega) = -1/pi "Im" "Tr" G(bold(k), omega + i delta) $

and plots it along the symmetry line. With `sw_self=True`, the self-energy from a FLEX calculation is incorporated to include interaction effects.

== option=5: Boltzmann Transport (`CONDUCTIVITY_BT`)

Calculates transport coefficients from the Boltzmann equation in the relaxation-time approximation. The relaxation-time model is chosen with `tau_mode` (constant $tau$ by default).

Thermoelectric quantities:
- Electrical conductivity tensor $sigma_(i j)$ (unit: S/m)
- Thermal conductivity tensor $kappa_(i j)$ (unit: W/m/K)
- Seebeck coefficient tensor $S_(i j)$ (unit: V/K)
- Power factor $sigma S^2$ (unit: $upright(W \/ m \/K^2)$)
- Peltier coefficient, Lorenz number

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

Solves the FLEX equations self-consistently to obtain the electron self-energy $Sigma(bold(k), i omega_n)$ in Matsubara frequency space. The result is saved to `sigma.bin` and `self_en.npz`.

- `sw_out_self=True`: write the self-energy to file
- `sw_in_self=True`: load the previous self-energy (`sigma.bin`) as the initial guess
- `sigma_in_scale`: factor multiplying the $Sigma$ seed loaded via `sw_in_self=True` (default 1.0). For $U$-annealing (continuation that raises $U$ stepwise, seeding each step from the `sigma.bin` converged at the smaller $U$) the weak-coupling relation $Sigma prop U^2$ suggests $(U_"new"\/U_"old")^2$. Unlike temperature annealing the Matsubara grid is unchanged, so no seed interpolation is needed (`sigma.bin` has no header — keep `Nw` identical between runs)
- Temperature annealing: for stepwise cooling at fixed `Nw` ($T_"new"\/T_"old" gt.eq 0.8$ or so) the raw reuse of `sigma.bin` is sufficient. Raw reuse amounts to a seed with the frequency axis compressed by $T_"old"\/T_"new"$, and this compression approximately mimics the growth of $Sigma$ upon cooling — near criticality it converges faster than a faithful interpolation. When `Nw` must change (the one case raw reuse cannot handle), the Python helper `plibs.regrid_sigma_bin(temp_new=..., Nw_new=..., w_scale=...)` re-interpolates `sigma.bin` onto the new mesh. The original temperature, `Nw` and orbital count are taken automatically from a metadata record appended at the end of `sigma.bin` by recent `io_sigma` versions (a trailing record, invisible to older readers); explicitly passed values that contradict the footer raise an error (files without the footer require explicit `temp_old`/`Norb`/`Nw_old`). The `sw_in_self` reader also checks the footer: a `sigma.bin` whose mesh (`Nw`/orbital count/k-points) disagrees with the run aborts with an error (a Fortran read silently accepts a record that is larger than requested, scrambling the layout — this check prevents that accident); a differing temperature alone is treated as T-annealing and only reported. The SOC path (`mkself_soc`) now writes the same single-record format, so the helper applies to its files as well (the retired element-wise SOC files are no longer readable) (cubic spline with the conjugate-symmetric extension $Sigma_(l m)(-i omega) = Sigma_(m l)(i omega)^*$ to negative frequencies; the outermost stored frequency, which carries the wrap-around artifact of the circular convolution, is discarded, and points beyond the cutoff are filled with a $c_0 + c_1\/(i omega)$ tail fitted on an interior window). For near-critical cooling use `w_scale=temp_old/temp_new` (the same compression as raw reuse); `w_scale=1` (default) gives the faithful interpolation
- `m_diis_num`: DIIS history length (default 5 if undefined)
- `sw_rescale_flex=True`: if the Stoner factor $S = max_bold(q) "eig"[chi_0(bold(q), 0) hat(S)]$ reaches 1 during the iteration, uniformly rescale $chi_0(bold(q), i nu_n)$ by $(1 - 10^(-4))\/S$ to avoid the magnetic-instability divergence
- `sw_chi0_tail=True`: evaluate $chi_0$ with a tail correction (default False; no-SOC path of option=14,15,23).  The bubble is split bilinearly as $chi_0 = "conv"[G] - "conv"[G_0] + chi_0^"ref"$: the slowly decaying $1\/(i omega)$ reference part (the $G_0$ bubble) is evaluated exactly in imaginary time, and only the fast-decaying residual ($G - G_0 tilde 1\/omega^2$) passes through the FFT convolution.  The Matsubara truncation error improves from $O(1\/N_w)$ to $O(1\/N_w^2)$ (validated against the exact Lindhard function: convergence order 1.0 -> 2.0).  The $chi_0$ stage costs about 3x; effective for $N_w gt.eq 64$ (at very small $N_w$ the tau-discretization prefactor of the reference bubble can dominate).  The improved $chi_0$ propagates automatically into $chi_s$/$chi_c$, the FLEX vertex $V_sigma$ and the pairing vertex $V_Delta$.

*Notes on the built-in approximations* (all are deliberate choices made for computational cost and numerical stability):

- The `sw_rescale_flex` rescaling shrinks $chi_0$ uniformly at ALL momenta and frequencies, not just the diverging mode. It is, however, only a transient aid for the early iterations: the self-energy feedback normally pulls the Stoner factor below 1 by convergence, the rescaling becomes inactive, and the converged solution is then a genuine unrescaled FLEX fixed point (no bias). If the rescaling is still active at the final iteration, the normal state is magnetically unstable at this $T$/$U$ ($T < T_N$ within FLEX) and stdout prints `[FLEX] WARNING: Stoner factor >= 1 at the final iteration`; such results correspond to an artificially weakened interaction and should not be used quantitatively (raise $T$, or seed from a converged higher-$T$ $Sigma$ via `sw_in_self`). The clamp target just below 1 ($1 - 10^(-4)$) is deliberate: the resulting $Sigma$ overshoot is the mechanism that kicks the loop off the rescaled manifold. A softer clamp (e.g. 0.98) makes near-critical systems "converge" quickly to a spurious rescaled fixed point even when a genuine $S < 1$ solution exists — the extra iterations caused by the oscillation are the price of this safety.
- A static part of the self-energy is subtracted at every iteration (wrapper argument `sub_sigma`; `1`: Hermitian part of $Sigma(bold(k), i omega_0)$, the default; `2`: frequency average; `0`: no subtraction). This is a prescription that absorbs the static band shift into the chemical-potential adjustment and stabilizes the Fermi surface; it is NOT a strict Hartree--Fock-only subtraction (the value at $i omega_0$ also contains low-energy dynamical components). Use `sub_sigma=0` if absolute band shifts are of interest.
- The Matsubara-frequency convolutions ($chi_0$, $Sigma$, Eliashberg kernel) are circular FFT convolutions of length $2 N_w$ with a sharp cutoff and, by default, no high-frequency tail correction (only the farthest point of $V$ is approximated by the bare vertex). The effective cutoff $omega_c = (2 N_w - 1) pi T$ shrinks proportionally to $T$, so a temperature scan at fixed $N_w$ carries an $O(1\/N_w)$ systematic drift. Choose $N_w$ so that $omega_c$ safely exceeds the bandwidth and `U` even at the lowest temperature. `sw_chi0_tail=True` improves the $chi_0$ truncation error to $O(1\/N_w^2)$. The tail error of the $Sigma$ convolution itself is dominated by an $omega$-independent constant (an HF-like static shift) that is absorbed by `sub_sigma`/the chemical-potential adjustment, and the Eliashberg-kernel tail error only enters the uniform ($s$-wave) component of the gap and vanishes by symmetry for sign-changing gaps — hence no separate corrections are implemented for those ($V_Delta$ and $V_sigma$ still benefit through the corrected $chi_0$).

== option=15: Linearized Eliashberg Equation (`LIN_ELIASHBERG`)

Solves the linearized Eliashberg equation as an eigenvalue problem using the effective pairing interaction from RPA or FLEX. The largest eigenvalue $lambda$ and the corresponding gap function $Delta(bold(k), i omega_n)$ are obtained. The superconducting transition temperature $T_c$ corresponds to $lambda = 1$. The default solver is a power method with shift + deflation; setting `arnoldi_m > 0` enables an Arnoldi solver.

- `sw_self=False`: no self-energy corrections (pure RPA)
- `sw_self=True`: uses the FLEX self-energy (requires `sigma.bin`)
- `sw_from_file=True`: loads the self-energy from `sigma.bin` without re-running FLEX
- `gap_sym`: sets the initial symmetry of the gap function (see below)
- Gap functions are written to `gap.npy` and to one file per orbital pair — `gap_11.dat`, `gap_12.dat`, … (`gap_{i j}.dat`, 1-based orbital indices)

== option=16: Gap Function Post-Processing (`GAP_FUNCTION`)

Reads the gap function from `gap.npy` (produced by option=15) and computes:

- The anomalous Green's function $F(bold(k), i omega_n)$
- Analytic continuation to the real-frequency axis via Padé approximation

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

Builds a real-space supercell Hamiltonian containing impurities and computes the spectral function $A(bold(k), omega)$ via the impurity Green's function. The impurity sites are specified in `imp_list`.

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
3. Each iteration: $Delta_"new" = T \/ N_k sum V_Delta dot F$ → amplitude-direction Newton (secant) acceleration + DIIS shape mixing (falls back to linear mixing pp=0.3 when bypassed) → update $Sigma$ → SC Dyson updates of $G, F$ → recompute $V_sigma, V_Delta$ from $chi^"GG", chi^"FF"$

- `m_diis_num`: DIIS history length ($>= 2$ enables Pulay extrapolation; otherwise linear mixing)
- `sw_self=True`: include FLEX self-energy via the $Sigma$-dressed Green's function
- `sw_from_file=True`: read self-energy from `sigma.bin`
- `sw_check_only`: see "Switch Variables" below
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

- `fname` (string): Path to the input file. The format depends on `ftype`; the file extension is usually omitted.

- `ftype` (integer): Input Hamiltonian file format (see Section 4).

- `brav` (integer): Bravais lattice type (see Section 5).

- `sw_soc` (bool): Switch to enable spin-orbit coupling.

- `inv_tol` (float, default $3 times 10^(-2)$): Relative residual below which the inversion operator is accepted as a symmetry of $H(bold(R))$ (option=14,15,16,23, which reconstruct $Delta(-bold(k))$ from $Delta(bold(k))$ with it). Fitted first-principles models carry a few percent of error, so the default is deliberately loose; raise it if a fitted model is rejected, and check the printed residual before trusting a run that needed a large value.

== Mesh Settings

- `Nx, Ny, Nz` (integer): Number of $bold(k)$-point mesh divisions in the first Brillouin zone along $x$, $y$, $z$. Used for 3D $bold(k)$-space integrations and FFT-based convolutions (FLEX, Eliashberg, chi). For 2D systems, set `Nz=1`. Powers of 2 (32, 64, ...) are preferred for FFT efficiency. Memory consumption scales as $tilde.equiv N_x N_y N_z dot N_w dot N_"orb"^2$.

- `Nw` (integer): Number of Matsubara frequencies for FLEX/Eliashberg/SC-chi calculations (option=12,13,14,15,16,23), or number of real-frequency points for DOS/spectral function calculations (option=1,4, etc.). Matsubara frequencies are $omega_n = (2n+1) pi T$ for $n = 0, 1, \ldots, N_w - 1$. Lower temperatures require larger $N_w$ (typical: 256–1024). The effective frequency cutoff $omega_c = (2 N_w - 1) pi T$ shrinks proportionally to $T$, and the Matsubara convolutions use a sharp cutoff (no tail correction), so a temperature scan at fixed `Nw` accumulates a systematic drift; choose `Nw` so that $omega_c$ safely exceeds the bandwidth and `U` at the lowest temperature (see also the notes under option=14).

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

- `tau_const` (float, unit: fs): Constant relaxation time $tau$ used when `tau_mode='const'`. In the constant relaxation-time approximation this is a free parameter typically chosen by comparison with experiment (typical metals: 1–100 fs). Note that $R_H$ carries $tau^2\/tau^2$ and is therefore *independent* of a constant $tau$: with `tau_mode='const'` the Hall coefficient is a pure band-structure number with no temperature dependence beyond the Fermi window.

- `tau_mode` (str, default `'const'`): Relaxation-time model for option=5.
  - `'const'`: constant `tau_const`, $bold(k)$- and band-independent.
  - `'epa'`: $tau(bold(k))$ from the electron--phonon averaged coupling (Samsonidze--Kozinsky EPA), read from `epa_file`. This is what makes $R_H (T)$ and sheet-dependent mobilities meaningful.
  - `'dos1'` / `'dos2'`: DOS-based $tau$, scaled by `tau_const`.
  A summary is printed with the Fermi-surface average $angle.l tau angle.r$, its range over the Fermi window, and $angle.l tau^2 angle.r \/ angle.l tau angle.r^2$ — which is 1 exactly when $tau$ is constant and quantifies how far $tau(bold(k))$ shifts $R_H$ away from the constant-$tau$ value.

- `epa_file` (str or `None`): Path to the `epa.x` (`job='egrid'`) output read when `tau_mode='epa'`.

- `hall_mesh_list` (list or `None`, default `None`): $bold(k)$-mesh convergence sweep for the Hall coefficient in option=5, e.g. `[20,30,40]` or `[(20,20,10),(30,30,15)]`. Each entry is either a single integer (cubic mesh) or an `(Nx,Ny,Nz)` tuple. `None` runs a single mesh `(Nx,Ny,Nz)` and gives no convergence information. Because $sigma^((1))$ can come out with the wrong sign on a coarse mesh, running a sweep at least once for a new material is strongly recommended.

- `sw_tdf` (bool): If `True`, the transport distribution function (TDF) is computed first, and transport coefficients are obtained by energy integration. Relevant when using an energy-dependent relaxation time.

Spin degeneracy in both transport modes (option=5,6) follows `sw_soc`: with `sw_soc=False` the bands are spin-degenerate and a factor $g_s=2$ is applied, with `sw_soc=True` both spins are already in the Hamiltonian and $g_s=1$ is used.

== Orbital and Interaction Parameters

- `olist` (list): Orbital indices for color plotting (`color_option=1`); see Section 7.

- `U` (float, unit: eV): On-site Coulomb repulsion (Hubbard $U$). Used in FLEX/RPA calculations. This is a key parameter controlling magnetic and superconducting instabilities.

- `J` (float, unit: eV): On-site Hund's coupling constant. The screened interaction $U' = U - 2J$ (Kanamori screening) is used automatically.

- `orb_dep` (bool): If `True`, orbital-dependent interaction matrices `Umat` and `Jmat` are used. If `False` (default), the constant values `U` and `J` are applied uniformly to all orbitals. The diagonal of `Umat` is the intra-orbital $U_(i i)$ and its off-diagonal is the inter-orbital $U'_(i j)$ itself (the Kanamori relation $U' = U - 2J$ is *not* re-applied), `Jmat` holds $J_(i j)$; both must be real and $N_"orb" times N_"orb"$, or one site worth of orbitals when `site_prof` lists several equally sized sites -- the per-site matrix is then copied onto every site (e.g. a 5-orbital set on a 2-Fe, 10-orbital model). Only pairs living on the same site enter the vertices, so the inter-site blocks are never read.

- `UJ_material` (string, optional): With `orb_dep=True`, loads a stored parameter set (e.g. cRPA $U$/$J$ of the Fe-based superconductors) instead of writing `Umat`/`Jmat` out by hand. The set is read from `UJ_dir` (default `UJ`), which may be a directory of per-material files `{UJ_dir}/{UJ_material}.json` or a single json file holding several materials keyed by name; each entry stores `U` and `J` either flat (row-major $N_"orb"^2$ values) or as nested rows. Bundled sets: `FeSe`, `FeTe`, `LiFeAs`, `BaFe2As2`, `LaFeAsO`, `LaFeAsOnakamura` (5 Fe-$d$ orbitals; they cover one site and are replicated over the sites of a multi-site model as described under `orb_dep`). The orbital order of the file must match that of the hamiltonian. `None` (default): use the `Umat`/`Jmat` written in the header.

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

Two caveats remain. The shell of a chiral pair is shared between its two partners (otherwise $|phi|$ is not invariant under the rotation that mixes them — hexagonal $E_(2g)$ fails by 77% with per-partner shells), which the code does automatically. And for a *multi-site* cell the pairing bond is $bold(R) + bold(tau)_j - bold(tau)_i$ rather than $bold(R)$, so the sublattice phase is still missing: for the 2-Fe cell of the iron-based systems, or any model where inversion exchanges sites, use the orbital-basis route (`eil_gap_orbital` / `eil_gap_file` with the Nagai–Nakamura band projection).

== $bold(k)$-Space Settings

- `kz` (float, in reduced coordinates): The $k_z$ value for the 2D Fermi surface plot (option=2). Ranges from 0 to 0.5 ($k_z = 0$: $Gamma$-plane, $k_z = 0.5$: zone boundary plane).

- `kscale` (list or float): Display scale for each axis in the 3D Fermi surface plot (option=3). Example: `kscale=[1.0, 1.0, 0.5]` compresses the $k_z$ direction by half.

- `k_sets` (list of lists): Coordinates of symmetry line endpoints in reduced units (0 to 1). Define this to specify a custom symmetry path instead of the auto-generated one. Example: `k_sets=[[0,0,0],[0.5,0,0],[0.5,0.5,0]]`

- `xlabel` (list of strings): Labels for the points in `k_sets`. Example: `xlabel=[r'$\Gamma$','X','M']`

- `at_point` (list): Coordinates of the $bold(q)$ point (in reduced units) for the single-$bold(q)$ susceptibility calculation (option=8).

== Switch Variables

- `sw_unit` (bool): If `True` (default), physical constants in SI units are used and output is in physical units. If `False`, a dimensionless system with $hbar = k_upright(B) = e = 1$ is used.

- `sw_omega` (bool): Switch for option=11 — compute the pairing susceptibility on the real frequency axis (`True`) or Matsubara frequency axis (`False`).

- `sw_self` (bool): If `True`, the FLEX self-energy is incorporated into option=4 (spectrum), option=15 (linear Eliashberg), or option=23 (nonlinear Eliashberg) via the renormalized Green's function.

- `sw_out_self` (bool): If `True`, the FLEX self-energy is written to `sigma.bin` and `self_en.npz`. Also used to trigger writing the gap function in option=15.

- `sw_in_self` (bool): If `True`, the previous FLEX self-energy is loaded from `sigma.bin` as the initial guess for the iterative self-consistent loop.

- `sigma_in_scale` (float, default 1.0): Factor multiplying the $Sigma$ seed loaded via `sw_in_self=True`. Intended for $U$-annealing; $(U_"new"\/U_"old")^2$ is a reasonable choice (see the notes under option=14).

- `sw_from_file` (bool): If `True`, the self-energy is read from `sigma.bin` and the FLEX calculation is skipped. The Eliashberg equation is then solved directly with this pre-computed self-energy.

- `sw_check_only` (bool): Used only by option=23 (nonlinear Eliashberg). If `True`, the routine stops right after the internal linearized-Eliashberg solve (reporting the Stoner factor $S$ and eigenvalue $lambda_"eliash"$) without running the nonlinear loop — handy for quickly bracketing $T_c$ via a temperature scan. Regardless of this flag, the nonlinear loop is also skipped automatically whenever $S >= 1$ or $lambda_"eliash" < 1$.

- `sw_rescale_flex` (bool): For FLEX (option=14), uniformly rescale $chi_0$ by $(1-10^(-4))\/S$ whenever the Stoner factor $S$ reaches 1 during the iteration, preventing the magnetic-instability divergence. This stabilizes the early iterations; the $Sigma$ feedback normally deactivates the rescaling by convergence, leaving the converged solution unbiased. If the rescaling is still active at the final iteration a warning (`[FLEX] WARNING`) is printed, meaning the normal state is magnetically unstable (see the notes under option=14).

- `sw_chi0_tail` (bool): For option=14,15,23 (no-SOC path), evaluate $chi_0$ with the tail-corrected convolution (default False). The Matsubara truncation error improves from $O(1\/N_w)$ to $O(1\/N_w^2)$, so a smaller `Nw` reaches the same accuracy (the $chi_0$ stage costs about 3x; effective for $N_w gt.eq 64$). See the notes under option=14.

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
- `eil_zeeman` (float, unit: eV): the Zeeman (Maki) field for the LDOS (surface: splits the $d_[110]$ ZEBS into $plus.minus h$; vortex: spin-splits the core states).

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

= Test Suite

The `tests/` directory contains regression tests for the numerical kernels and physics benchmarks: 173 tests across seven files. Each file can be run directly as a Python script, so `pytest` is optional.

Before running the tests, the Fortran shared library `libs/libfmod.so` must be compiled. If it is missing, enter the `libs` directory and run `make FC=<compiler> SL=<library>`.

== How to Run

Run individual test files directly:

```bash
python tests/test_transport.py
python tests/test_rpa_flex.py
```

Run everything with a summary table (no `pytest` needed):

```bash
python tests/run_all.py          # everything
python tests/run_all.py rpa      # only files matching "rpa"
```

If `pytest` is available, the whole directory can also be run with `pytest tests`. A handful of tests take a `tmp_path` fixture and are reachable only under `pytest`; the standalone runner skips them.

== Coverage by file

#table(
  columns: (auto, auto, auto, 1fr),
  align: (left, right, right, left),
  table.header([*File*], [*Tests*], [*Time*], [*What it pins down*]),
  [`test_eilenberger.py`], [57], [437 s],
    [Quasiclassical Eilenberger/Riccati solvers: homogeneous limits, Fortran kernels against Python references, surface Andreev states, vortices and self-consistent vortex lattices, condensation free energy, model and Wannier Fermi surfaces, Pauli limiting, triplet $d$-vector textures],
  [`test_rpa_flex.py`], [34], [0.7 s],
    [RPA/FLEX/Eliashberg building blocks: interaction vertices, RPA algebra, $chi_0$ and $phi_0$ against exact Lindhard/pair sums, tail-corrected $chi_0$, self-energy regridding, and the irreducible $bold(k)$-mesh under time reversal],
  [`test_effmass.py`], [29], [3.8 s],
    [Inverse effective mass in the orthogonal and MLO bases (interband and overlap terms, degeneracies), and the dHvA extremal-orbit machinery: slice geometry, open-orbit rejection, and the zone reduction],
  [`test_nmr.py`], [22], [16 s],
    [NMR in the SC state: BCS gap interpolation, $bold(q)$-mesh folding, loading a self-consistent Eliashberg gap, Knight shift (Yosida vs triplet), Hebel--Slichter peak and the $T^3$ law of line nodes],
  [`test_transport.py`], [18], [0.8 s],
    [Kubo$arrow.l.r$Boltzmann equivalence in the dc limit, Wiedemann--Franz, the symmetrized interband heat-current vertex, band degeneracies, and the weak-field Hall response: $sigma^((1))$ against the scalar reference, $R_H arrow -1\/(n e)$ and its independence of the lattice angle, and Luttinger-volume carrier densities],
  [`test_velocity_mlo.py`], [7], [0.4 s],
    [Band velocity in a non-orthogonal (MLO) basis: the $-epsilon_n C^dagger (partial S\/partial k) C$ overlap term, Hermiticity, and bit-for-bit dispatch back to the orthogonal path],
  [`test_tools.py`], [6], [0.5 s],
    [The shared references the other files check against: Fermi/Bose identities, closed-form model bundles, and the exact $chi_0$ used to validate the Fortran convolution],
)

The whole suite takes about 7.5 minutes, and `test_eilenberger.py` is essentially all of it (96%) — the quasiclassical solvers are self-consistent field problems on a 2D grid. No single test dominates: the four heaviest (`test_vortex_in_plane_field`, `test_chiral_vortex`, `test_surface_gap_heals_to_bulk`, `test_free_energy_selects_the_chiral_partner`) take 57–77 s each. Reduce their `Ng`/`nbeta`/`ngrid`/`eil_fs_nkz` if a faster suite is needed; every other file finishes in under 20 s.

= Troubleshooting

- *Chemical potential does not converge*: Increase `Nx, Ny, Nz`, or slightly raise `tempK`. A coarse $bold(k)$-mesh can cause instability in the self-consistent $mu$ search.

- *Too much noise in the spectrum*: Increase `delta` slightly, or increase `Nw` to improve the frequency resolution.

- *Band structure looks wrong*: Check that `brav` is set correctly. In particular, when interfacing with QE, pay careful attention to the conventions for FCC (`brav=1`) and BCC (`brav=2`).

- *Linear Eliashberg does not converge*: Try reducing `U`, increasing `tempK`, or using a finer mesh (`Nx, Ny, Nz`). If the eigenvalue $lambda$ exceeds 1, the system is inside the superconducting phase at that temperature.

- *Nonlinear Eliashberg diverges*: Increase `m_diis_num` to 5–10, raise `tempK`, or increase `Nw`. Running with `sw_check_only=True` first to inspect the Stoner factor $S$ and eigenvalue $lambda_"eliash"$ is also useful. If $S$ exceeds 1, the system is magnetically unstable and `U` must be reduced; if $lambda_"eliash" < 1$, you are still above $T_c$ and need to lower the temperature (in either case the nonlinear loop is skipped automatically).

- *FLEX self-energy diverges*: Set `sw_rescale_flex=True` or reduce `U`.

- *Transport coefficients have wrong units*: Make sure `sw_unit=True`. Setting it to `False` switches to a dimensionless unit system.

- *The Hall coefficient changes sign with the $bold(k)$-mesh*: This is the normal failure mode of $sigma^((1))$, not a bug. Run a sweep with `hall_mesh_list` and read the two diagnostics: if $N_"eff"$ is below $10^3$ the Fermi window is carried by too few states, and if the cancellation ratio is below $10^(-2)$ the system is near compensation and needs a much finer mesh than $sigma_(x x)$ does. Raising `tempK` also widens the Fermi window.

- *The Hall carrier density disagrees with experiment*: Check it against the Luttinger volumes reported in the same output. If $n_e \/ n_h approx 1$ the metal is compensated and $n_H$ is not a carrier density at all — compare mobilities or fit a two-fluid model instead. Also verify that `sw_soc` is set correctly, since it selects the spin degeneracy factor, and that $omega_c tau$ times your experimental field is still $lt.tilde 0.1$.

- *$R_H$ shows no temperature dependence*: With `tau_mode='const'` this is expected — $R_H$ carries $tau^2\/tau^2$ and a constant $tau$ cancels exactly. Use `tau_mode='epa'` with an `epa_file` to get a $bold(k)$-dependent $tau$.
