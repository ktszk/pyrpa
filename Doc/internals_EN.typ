//document format
#set text(font: ("Times New Roman", "Liberation Serif"), size: 11pt)
#set page(paper: "a4", margin: (x: 2cm, y: 2cm))
#set math.equation(numbering: numbering.with("(1)"), supplement: "Eq.")
#set heading(numbering: "1.")
#import "callgraph.typ": module-graph-plibs, module-graph-kernels, option-graphs

#align(center)[
  #text(size: 18pt, weight: "bold")[pyrpa Developer Notes]
  #v(0.5em)
  #text(size: 11pt)[Internal functions, implementation, validation and change history]
]

#v(1em)

#outline(title: "Contents", depth: 2, indent: auto)

#pagebreak()

= About these notes

The user documentation is `Doc/manual_EN.typ`. These notes collect:

- arguments that `main.py` does not set, used when calling the functions of `libs.plibs` / `libs.flibs` directly
- internals of the build and of the implementation
- file-format compatibility
- accuracy checks
- the test coverage and the change history

Items refer to the manual by the name of the corresponding section.

= Code layout

- `main.py`: the settings (variables at the top) and the driver part: `setup()` (input checks, reading the model, lattice, interaction and chemical potential) and one `run_*` function per `option` (table `RUNNERS`).
- `libs/plibs`: the Python side of the calculations (model input, symmetry, drivers); no figures.
- `libs/flibs`: ctypes wrappers of the Fortran kernel `libfmod.so`.
- `libs/plots.py`: matplotlib figures.
- `libs/src/*.f90`: Fortran kernels ($chi$, FLEX, Eliashberg, transport, Eilenberger, ...).

#table(
  columns: (auto, 1fr),
  table.header([*Module*], [*Contents*]),
  [`plibs._lattice`], [lattice vectors, hopping input (all `ftype`), $bold(k)$/$bold(r)$ meshes, BZ],
  [`plibs._symmetry`], [point-group operations, symmetry detection, irreducible $bold(k)$ points],
  [`plibs._bands`], [eigenvalues, chemical potential, Fermi surface, cross-sections],
  [`plibs._response`], [spin / pairing susceptibility drivers, gap form factors, SC-state $chi$],
  [`plibs._calc`], [FLEX, Eliashberg, spectrum and CPA drivers; `sigma.bin` and gap I/O],
  [`plibs._eliash_nosym`], [symmetry-free linearized Eliashberg (Zeeman, FFLO)],
  [`plibs._freq_grid`], [IR-basis frequency grid],
  [`plibs._conductivity`], [Boltzmann and Kubo transport, Hall coefficient, relaxation-time models],
  [`plibs._optical`, `_color`], [optical constants and material colour],
  [`plibs._fermiology`], [cyclotron mass, dHvA, band masses on a path],
  [`plibs._impurity`], [impurity supercell and the unfolded spectrum],
  [`plibs._nmr`], [NMR in the SC state],
  [`plibs._phonon_pairing`], [phonon-mediated pairing from EPA],
  [`plibs._wannier_io`], [EPA input, U/J parameter sets, Wannier real-space I/O],
  [`plibs._maxent`], [MaxEnt analytic continuation (`ana_cont`)],
  [`plibs._eilenberger*`], [quasiclassical Eilenberger (homogeneous / model FS, surface, vortex and lattice, spin-matrix version)],
  [`plibs._units`], [physical constants (`sw_unit`)],
  [`flibs._loader`], [selecting and loading `libfmod*.so`, OpenMP settings],
  [`flibs._hamiltonian`, `_green`], [$H(bold(k))$, eigenvalues, Green's functions],
  [`flibs._susceptibility`, `_interaction`], [$chi_0$, RPA, interaction vertices],
  [`flibs._selfenergy`, `_eliashberg`], [FLEX self-energy, linearized and nonlinear Eliashberg],
  [`flibs._transport`], [transport kernels],
  [`flibs._freq`, `_grid`, `_impurity`, `_eilenberger`], [frequency transforms, vortex grids, impurities, Riccati kernels],
)

= `option` in `main.py` and the functions it calls

For every `option`, `setup()` first reads the model (`plibs.import_hoppings` / `import_MLO_hoppings` / `read_HamRsMLO`), sets up the lattice (`plibs.get_ptv`, `get_bvec`) and, where needed, the chemical potential (`plibs.get_mu`), the inversion operator (`plibs.inversion_op`), the orbital pairs of $chi$ (`plibs.get_chi_orb_list`) and the interaction matrices (`flibs.gen_SCmatrix`, `gen_SCmatrix_orb`). The main functions the `run_*` function of each `option` then calls are:

#table(
  columns: (auto, 1fr),
  table.header([*option*], [*main functions*]),
  [0 `BAND`], [`plibs.mk_klist`, `get_eigs` → `plots.plot_band`],
  [1 `DOS`], [`flibs.gen_dos` → `plots.plot_dos`],
  [2 `FERMI_2D`], [`plibs.get_eigs_2d`, `get_kf_points`, `get_colors` → `plots.plot_FS`],
  [3 `FERMI_3D`], [`plibs.gen_3d_surf_points` → `plots.plot_3d_surf`],
  [4 `SPECTRUM`], [`plibs.spectrum_on_path` → `plots.plot_path_map`],
  [5 `CONDUCTIVITY_BT`], [`plibs.calc_conductivity_Boltzmann`, `get_hall_coe`],
  [6 `CONDUCTIVITY_PT`], [`plibs.calc_conductivity_lrt` → `plots.plot_optical`],
  [7, 10], [`plibs.mk_qlist`, `calc_path_spectrum`],
  [8], [`plibs.chis_q_point` → `plots.plot_chis_qpoint`],
  [9, 11], [`plibs.chis_qmap`, `phi_qmap` → `plots.plot_qmap`],
  [12, 13], [`plibs.sc_gap_on_mesh`, `chis_spectrum_sc`, `chis_q_point_sc`, `bdg_bands_path`],
  [14 `FLEX`], [`plibs.calc_flex` / `calc_flex_soc`],
  [15 `LIN_ELIASHBERG`], [`plibs.calc_lin_eliashberg_eq` / `calc_lin_eliash_soc` / `calc_lin_eliash_nosym`],
  [16 `GAP_FUNCTION`], [`plibs.output_Fk`],
  [17 `BAND_FILLING`], [`plibs.report_band_filling`],
  [18 `CYCLOTRON_MASS`], [`plibs.get_mass`],
  [19 `DHVA`], [`plibs.get_dhva_band` → `plots.plot_dhva`],
  [20 `ELECTRON_MASS`], [`plibs.band_mass_path`, `report_band_mass`],
  [21 `SPECTRUM_IMPURITY`], [`plibs.random_impurity_sites`, `calc_impurity_spectrum`],
  [22 `SIGMA_CPA`], [`plibs.calc_cpa`],
  [23 `NONLIN_ELIASHBERG`], [`plibs.calc_eliashberg_eq`],
  [24 `EILENBERGER`], [`plibs.calc_eilenberger`, `calc_penetration_depth`, `calc_fs_penetration`, `calc_pauli_limit`, `calc_spin_pauli`, `calc_free_energy`],
  [25 `EILENBERGER_SURFACE`], [`plibs.calc_surface`, `calc_surface_dvector`],
  [26 `EILENBERGER_VORTEX`], [`plibs.calc_vortex`, `calc_vortex_lattice_periodic`, `calc_vortex_lattice_sc`, `calc_vortex_current`, `calc_vortex_maxwell`, `calc_vortex_chiral`, `calc_vortex_dvector`],
  [27 `NMR_SC`], [`plibs.nmr_qmesh`, `nmr_sweep`, `write_nmr_dat` → `plots.plot_nmr`],
)

= Call graphs

Call graphs (the counterpart of doxygen's call graphs) built by `Doc/gen_internals.py` from a static analysis of the sources. They are drawn with the Typst package diagraph (Graphviz compiled in), so no `dot` command is needed.

- Colours: #box(fill: rgb("#e8e8e8"), inset: 2pt)[`main.py`], #box(fill: rgb("#dbe9f7"), inset: 2pt)[`plibs`], #box(fill: rgb("#dcf0dc"), inset: 2pt)[`flibs`], #box(fill: rgb("#fbe3c8"), inset: 2pt, radius: 3pt)[Fortran (rounded; source file in brackets)]
- Arrows point from caller to callee. In the module graphs the line width counts the function-to-function calls.
- The graph of each `option` starts at the `run_*` function of `main.py` and follows 5 levels and up to 40 functions; a function with further callees ends in a dashed "…".
- The argument marshalling of `flibs._loader` (`i64`, `dbl`, `c_array`, ...) appears in every call and is left out.
- Being static, the analysis misses calls through variables holding functions and the like.

== Dependencies between the `plibs` modules

#module-graph-plibs

== `plibs` → `flibs` → Fortran

#module-graph-kernels

== Calls by `option`

#option-graphs

= Function reference

The public functions of `libs.plibs` and `libs.flibs` (the `__all__` of each submodule). Each summary is the first sentence of the function's docstring (its `@brief` if present), generated by `Doc/gen_internals.py`; see the docstrings for argument types and details. Rerun `python Doc/gen_internals.py` after changing the API.

#include "api_reference.typ"

= Library arguments (not set from `main.py`)

== Eigenvalue solver of the linearized Eliashberg equation (`arnoldi_m`)

`flibs.linearized_eliashberg` and `flibs.linearized_eliashberg_soc` take `arnoldi_m` (default 10): `arnoldi_m > 0` selects the Arnoldi method, `0` a power method with shift + deflation. `main.py` (`plibs.calc_lin_eliashberg_eq`, `calc_lin_eliash_soc`) uses the default and therefore always solves with Arnoldi.

With SOC the Arnoldi method is required. The power method has no branch tracking and converges to the dominant state of the parity sector, which may be the odd-frequency partner; the Arnoldi solver returns the branch whose spin character matches `gap_sym`.

The Krylov dimension of the symmetry-free version (`flibs.linearized_eliashberg_nosym`) is `nosym_arnoldi_m` in `main.py`.

== MaxEnt analytic continuation (`sw_maxent`)

`plibs.spectrum_on_path` (option=4) and `plibs.output_Fk` (option=16) use MaxEnt instead of Padé with `sw_maxent=True` (default False); the `ana_cont` package is required. `main.py` uses Padé only. For a $Sigma$ on the IR basis, MaxEnt is stable regardless of `ir_eps`.

== Static subtraction in FLEX (`sub_sigma`)

`flibs.mkself`, `mkself_soc` and the nonlinear Eliashberg wrapper take `sub_sigma`. Each iteration subtracts a static part of $Sigma$ so that the static band shift is absorbed into the chemical potential and the Fermi surface stays stable.

- `1`: the Hermitian part of $Sigma(bold(k), i omega_0)$ (default; used by `main.py`)
- `2`: the frequency average (uniform grid only; not available on the IR basis)
- `0`: none

This is not an exact Hartree--Fock subtraction; use `sub_sigma=0` when discussing absolute band shifts. In a model whose $chi$ pairs are restricted to the correlated orbitals (`chiolist`), the static part of $Sigma$ enters the correlated orbitals alone and effectively shifts $epsilon_d - epsilon_p$, so mind the double counting. In the AiiDA interface it is set with the `sub_sigma` key of `options`.

== Regridding `sigma.bin` (`plibs.regrid_sigma_bin`)

When the `sigma.bin` read by `sw_in_self=True` has a different `Nw` from the current run, `plibs._auto_regrid_sigma_seed` calls `regrid_sigma_bin(temp_new=temp, Nw_new=Nw, w_scale=temp_old/temp)` and rewrites the file in place. To call it by hand:

```python
plibs.regrid_sigma_bin(temp_new=..., Nw_new=..., w_scale=...)
```

- `w_scale=temp_old/temp_new` reproduces the frequency compression of reading an equal-`Nw` `sigma.bin` as it is (for cooling near criticality); `w_scale=1` (the function default) is a faithful interpolation.
- The interpolation is a cubic spline extended to negative frequencies with $Sigma_(l m)(-i omega) = Sigma_(m l)(i omega)^*$. The last point (contaminated by circular-convolution wrap-around) is dropped, and beyond the cutoff a $c_0 + c_1\/(i omega)$ tail fitted on an inner window is used.
- SOC runs (`mkself_soc`) write the same format and work with the same function.
- The source grid is taken from the metadata record (temperature, `Nw`, orbital count) at the end of `sigma.bin`; explicitly given values that contradict it raise an error. For files without the record see "File formats and compatibility".

== CPA concentration and potential

`plibs.calc_cpa(..., x_cpa=0.5, VB_shift=1.0)` is the CPA of a binary alloy: `x_cpa` is the concentration of species B and `VB_shift` its on-site level shift [eV]. `main.py` (option=22) passes the fixed values above.

== Gap form factors (`set_harmonic_lattice` and friends)

`main.py` registers the lattice for the form factors automatically (`plibs.set_harmonic_lattice(avec, rvec)`; the $s^plus.minus$ shell is set with `s_rcut=`, which receives `gap_rcut` of `main.py`). `plibs.lattice_harmonic()` builds one harmonic directly, and `gap_symms(..., avec=, rvec=)` takes the lattice explicitly. `get_initial_gap` prints the $|bold(R)|$ of the shell actually used.

On an orthogonal single-site lattice the construction reproduces the former fractional formulas ($cos(2 pi k_x) - cos(2 pi k_y)$ etc.) bit for bit, for every entry of the manual's table and any axis ratio.

== Vortex-lattice symmetry (`calc_vortex_lattice_symmetry`)

`plibs.calc_vortex_lattice_symmetry` is a driver not called from `main.py`. It minimizes the Ichioka–Machida lattice free energy over the cell apex angle and the gap-vs-lattice orientation $theta_0$ to determine the *stable vortex-lattice symmetry* and its field evolution (e.g. the $d$-wave triangular → square transition near $H_(c 2)$). When a Wannier FS is supplied, $theta_0$ rigidly rotates the whole crystal (FS + gap), so the Fermi-velocity anisotropy also enters the selection.

= Build internals

- *Automatic `SL` detection.* It follows the loaded modules (`MKLROOT`, `AOCL_ROOT`) because selecting a library whose root is not set yields a `libfmod.so` with unresolved BLAS symbols that fails only at import.
- *gfortran and BLAS integer kinds.* gfortran rejects calls of the same external procedure with different integer kinds in one file, so the integer arguments of BLAS/LAPACK have one kind per file (`c_int64_t` in `fmod.f90`).
- *CMake.* The makefile is the reference implementation. A C++ compiler is needed only for the pocketfft path; the MKL path does not enable the CXX language at all, so it configures on oneAPI installs that ship `ifx` but no `icpx`. With `flang` and the AOCC module the makefile also scrubs `CPATH` and `LIBRARY_PATH` for child processes, which the CMake build does not.

= FLEX and Eliashberg implementation

- *Clamp target of `sw_rescale_flex`.* Rescaling $chi_0$ by $(1 - 10^(-4))\/S$, just below 1, is deliberate: the overshoot of $Sigma$ right after clamping is what pulls the iteration out of the spurious rescaled state. Relaxing it to about 0.98 was found to make even near-critical systems that have a proper $S < 1$ fixed point "converge" quickly to the spurious one.
- *Temperature-annealing seed.* Reading an equal-`Nw` `sigma.bin` as it is gives a seed whose frequency axis is compressed by $T_"old"\/T_"new"$; near criticality this was found to converge faster than a faithful interpolation.
- *`sigma.bin` checks.* The reader checks `Nw`, the orbital count and the $bold(k)$-point count to guard against Fortran silently reading a larger record in a scrambled order.
- *Amplitude-direction Newton acceleration of the nonlinear Eliashberg loop.* See the `sw_amp_newton` comments in `libs/src/ffeliash.f90`.
- *Check of `sw_chi0_tail`.* Against the exact Lindhard function the convergence order of $chi_0$ rises from 1.0 to 2.0.
- *Writing $Sigma$ from the IR basis.* When `self_en_wannier.npz` is written in the uniform-$omega_n$ format, the constant $Sigma(i omega -> infinity)$, which the basis cannot represent, is estimated from the Hermitian part at the three highest sampling points together with $1\/omega^2$, $1\/omega^4$ terms and removed before the expansion.
- *Sign of the Fourier transform for `eil_gap_file`.* The minus sign in $Delta_"orb"(bold(k)) = sum_bold(R) e^(-i 2 pi bold(k) dot bold(R)) Delta(bold(R))$ is the exact inverse of the `np.fft.ifftn` the exporter uses, verified by round-tripping an export.

= File formats and compatibility

- *`sigma.bin`.* A metadata record (temperature, `Nw`, orbital count) is appended at the end, harmless to readers that predate it. Files without the record are not regridded automatically; pass `temp_old`/`Norb`/`Nw_old` to `regrid_sigma_bin`. The old per-element record format of the SOC version is not readable.
- *$bold(k)$ mesh in `self_en_wannier.npz`.* `tau_mode='flex'` reads the FLEX $bold(k)$ mesh from the file; for older files without it the mesh is inferred from the $bold(R)$ extent, so rerun FLEX with `sw_out_self=True`.
- *Orbital order of gap files (`orb_conv`).* `gap_wannier.npz` stores the orbital pair in the physical order $Delta_(a b)$ (the exporter transposes the reversed ctypes view of the Fortran array once, on write, and stamps the npz with `orb_conv`). A gap exported before the stamp existed is still read correctly, with a one-line notice; re-export it to silence the notice. The order matters: $Delta$ is Hermitian but not symmetric, and its antisymmetric part changes sign under the transpose (63% of $max |Delta|$ for the 000AsP RPA solution, enough to make a fully gapped $s_(plus.minus)$ pocket look nodal). `eil_gap_file` and `nmr_gap_file` handle it the same way.
- *U, J parameter sets.* JSON (`.json`) files are read besides TOML (`.toml`); a `.toml` of the same name takes precedence.

= Accuracy checks

== IR basis

- *Default `ir_wmax`.* For 000AsP at $T = 0.02$ eV, $omega_"max" = 2.88$ eV shifts $mu$ by 0.27 meV, while 6.0 eV agrees with the uniform-grid extrapolation to $6 times 10^(-6)$ eV; hence the default $2 times max|epsilon_bold(k) - mu|$.
- *Agreement with the extrapolated uniform grid.* $chi_0$ matches the exact Lindhard sum to $10^(-8)$, the FLEX $mu$ and Stoner factor to $10^(-5)$ or better, the linearized Eliashberg $lambda$ to $10^(-7)$--$10^(-8)$, and the converged nonlinear Eliashberg gap to $10^(-4)$ or better.
- *Uniform-grid error.* A one-orbital $16 times 16$ model at $T = 0.02$ eV with $lambda = 1.34$ still gives a 3% smaller gap at $N_w = 4096$. In FLEX the truncation of $Sigma$ leaves an $O(1\/N_w)$ error even with the $chi_0$ tail correction, so $mu$ can be off by about 1 meV at $N_w = 1024$.
- *Why Padé uses the sampling points.* Padé is ill-conditioned: random relative errors of $10^(-10)$ in its input change the spectrum by tens of percent. IR values re-evaluated on uniform $omega_n$ are more accurate on the Matsubara axis than the uniform grid (`square_2orb.hop`, $16 times 16$, $T = 0.05$ eV: $9 times 10^(-6)$ from the extrapolation for $G$, versus $2 times 10^(-4)$ at uniform $N_w = 4096$), but the remaining basis truncation error is smooth and systematic, and Padé continues it into spurious real-axis structure. The accuracy is set by `ir_eps`, almost independently of the FLEX convergence threshold.
- *Padé accuracy* (relative L2 difference from the Richardson-extrapolated uniform-grid spectrum, median / max over $bold(k)$):
  - `square_2orb.hop`: `ir_eps=1e-8` 12% / 45%, `1e-12` 4% / 34%, `1e-14` 3% / 6% (uniform $N_w = 1024$: 1.5% / 3%)
  - 000AsP ($16 times 16 times 4$, $T = 0.05$ eV): `1e-12`, `1e-14` 4% / 13% ($N_w = 256$: 3% / 14%)

== MLO basis

- Without the tail correction at $N_w = 128$ (a 40 eV window), $lambda$ of the Cu MLO model was off by several times (its bands reach about $+30$ eV).
- Without the normalization, ecalj's MLOs are far from unit norm: $1\/S_(a a)(0)^2$ is about 4000 for the s orbital of Cu.

== Others

- *`sw_dec_axis`.* $E(bold(k))$ at a given physical $bold(k)$ is preserved to $10^(-14)$ for every `brav`.
- *dHvA (option=19).* The geometry and unit conversion are covered by analytic-cylinder tests, including the $F(theta)=F(0)/cos(theta)$ law for a tetragonal lattice. As a material-level check, the 10-orbital FeS Wannier model with experimental lattice constants gives $0.489$–$3.123$ kT for $B parallel c$, consistent with the GGA frequencies of roughly $0.5$–$3$ kT in Fig. 3(b) of T. Terashima _et al._, Phys. Rev. B *99*, 134501 (2019).

== Background: DFT bands versus experiment

The measured frequencies in the FeS work above are smaller, $0.15$–$1.87$ kT. This gap is expected and is not a defect of the calculation, and it is not specific to the iron pnictides and chalcogenides. In a strongly correlated metal a DFT (LDA/GGA) band structure omits the band narrowing and the orbital-dependent band shifts caused by electron correlation, so it systematically overestimates the size of the Fermi-surface pockets while underestimating the mass enhancement and hence the density of states at $E_F$. The same bias carries over to the pairing calculations of this code: a Cooper instability evaluated on an unrenormalized DFT band structure comes out too weak, so the RPA/Eliashberg eigenvalue and $T_c$ are underestimated.

Outside the iron-based family the discrepancy can be far larger.

- In the heavy-fermion antiferromagnet $"CeRhIn"_5$ a band calculation that treats the $4f$ electrons as itinerant gives the large Fermi surface, whereas the dHvA Fermi surface is essentially that of the $4f$-less reference $"LaRhIn"_5$, with cyclotron masses about an order of magnitude heavier (H. Shishido _et al._, J. Phys. Soc. Jpn. *71*, 162 (2002)).
- In $"Na"_x "CoO"_2$ the six small $e'_g$ hole pockets predicted by the LDA (D. J. Singh, Phys. Rev. B *61*, 13397 (2000)) are absent in photoemission, which sees only the $a_(1g)$ sheet (H.-B. Yang _et al._, Phys. Rev. Lett. *92*, 246403 (2004)); whether local correlations alone remove them is still debated.

= Tests

How to run them is described in `tests/README.md` and in the manual's appendix. The whole suite takes about 7 minutes (ifx build, 24 threads), most of it in `test_eliash_nosym.py` (solved on the full Brillouin zone) and `test_eilenberger.py` (self-consistent field problems on 2D grids). `test_aiida.py` is skipped without `aiida-core` and `test_freq_grid.py` without `sparse-ir`.

#table(
  columns: (auto, 1fr),
  align: (left, left),
  table.header([*File*], [*What it checks*]),
  [`test_bands.py`], [Bands and Fermi-surface geometry: reading Wannier90 `_hr.dat` and the phase of `gen_ham`, the first Brillouin zone of the 3D Fermi surface (faces, folding), band velocity in the non-orthogonal (MLO) basis (overlap term, Hermiticity), inverse effective mass in orthogonal and MLO bases (interband and overlap terms, degeneracies), cyclotron and dHvA extremal orbits (section geometry, rejection of open orbits, zone reduction)],
  [`test_transport.py`], [Kubo$arrow.l.r$Boltzmann agreement in the DC limit, Wiedemann--Franz law, symmetrized interband heat-current vertex, band degeneracy, $bold(k)$-point weights, weak-field Hall response, relaxation-time models (EPA, FLEX)],
  [`test_optical.py`], [Optical constants and colour: reflectance-to-colour conversion (CIE XYZ, D65, sRGB), refractive index, reflectance, absorption, transmission, Kubelka--Munk],
  [`test_impurity.py`], [Impurity supercell and the unfolded spectrum (phase sign, sum rule), unfolding and CPA in a non-orthogonal basis],
  [`test_rpa_flex.py`], [Building blocks of RPA / FLEX / Eliashberg: interaction vertices, RPA matrix algebra, $chi_0$ and $phi_0$ against exact Lindhard and pair-bubble sums, the tail-corrected $chi_0$, `mkself` / `mkself_soc` and the double-counting cancellation, the Nambu solution of the nonlinear Eliashberg equation, self-energy regridding, the irreducible $bold(k)$ mesh under time reversal; Green function, FLEX and Eliashberg in the non-orthogonal (MLO) basis (matrix Dyson equation, covariance under scaling, orbital normalization); phonon-mediated pairing from EPA],
  [`test_eliash_nosym.py`], [Symmetry-free linearized Eliashberg (Zeeman field, FFLO), consistency with the SOC version (bare term $-V\/2$, triplet seed and parity projection, band-basis output), exact-kernel check of multi-orbital triplets, branch tracking under SOC],
  [`test_freq_grid.py`], [IR basis: reproducibility of the basis, exact $chi_0$, agreement of FLEX, linearized / nonlinear Eliashberg, SOC, `sw_nosym` and the phonon vertex with the $N_w -> infinity$ limit of the uniform grid, the tail correction being ignored, Padé from the sampling points],
  [`test_nmr.py`], [Dynamic susceptibility in the SC state (including the non-orthogonal BdG) and NMR: BCS gap interpolation, $bold(q)$-mesh folding, loading a self-consistent Eliashberg gap, Knight shift (Yosida and triplet), Hebel--Slichter peak, $T^3$ law of line nodes],
  [`test_eilenberger.py`], [Quasiclassical Eilenberger / Riccati solvers: homogeneous limits, Fortran kernels against Python references, surface Andreev bound states, vortices and the self-consistent vortex lattice, condensation free energy, model / Wannier Fermi surfaces, Pauli limiting, triplet $d$-vector textures; Fortran grid kernels of the vortex solvers and the Riccati chord kernel],
  [`test_main.py`], [Chemical potential of the Eilenberger modes in `main.py` (Wannier Fermi surface and fixed $mu$)],
  [`test_aiida.py`], [AiiDA interface (`examples/aiida`): steps and workflows agree with direct `plibs` / `flibs` calls (RPA / FLEX / Eliashberg, MLO, SOC, Eilenberger), `options`, the recorded lattice, provenance],
  [`test_tools.py`], [Reference implementations the other files check against (Fermi / Bose identities, closed-form models, exact $chi_0$), CPU-based library selection (AMD, Intel), runtime settings of the gfortran build (`OMP_WAIT_POLICY`, stdout buffering)],
)

= Change history

The following fixes change results obtained with earlier versions; rerun the affected calculations (for older changes see the git history).

- *FLEX and Eliashberg in the non-orthogonal MLO basis (ftype=3,4).* With a self-energy, $G^(-1)$ lacked the overlap $S$ ($G^(-1) = i omega + mu - H - Sigma$) and was inconsistent with the particle-number correction. For the Cu MLO model the FLEX $mu$ moved by about $+2.9$ eV from the bare value and the Stoner factor became negative. Now $G^(-1) = (i omega + mu) S - H - Sigma$ and the particle number is $"tr"(S G)$. *All earlier MLO FLEX results, Eliashberg results with FLEX, and spectra with a self-energy are wrong.* RPA ($Sigma = 0$) built $G_0$ from the eigenvectors and is unchanged.
- *Normalization of MLO orbitals.* ecalj's MLOs are not of unit norm (Cu: $S_(a a)(0) = 0.016$ for s and $0.565$ for d), so a U placed on the raw orbitals acted roughly as $U\/S_(a a)(0)^2$ (Cu RPA with U = 1 gave a Stoner factor of 85). Each orbital is now normalized on reading (bands unchanged). U, J, orbital weights and orbital-basis $Sigma$ and $Delta$ change, and files saved before the normalization use a different convention. Interband velocity matrix elements (through the Löwdin orthogonalization) change by a few percent, and so do the optical conductivity and reflectance (Cu: $Delta R lt.eq 0.02$). Group velocities and quantities that use only them (Boltzmann transport, effective masses, dHvA) are unchanged.
- *`sw_chi0_tail` in the MLO basis.* The jump of $G_0$ at $tau = 0$ was assumed to be the unit matrix, so in a non-orthogonal basis the truncation error stayed $O(1\/N_w)$ instead of $O(1\/N_w^2)$. The jump $S^(-1)(bold(k))$ and the $bold(k)$ convolution of the density matrix are now used.
- *Nonlinear Eliashberg equation (option=23).* The normal Green's function of the SC state is solved from the Nambu Dyson equation $G = [G_0^(-1) + Delta G_0^dagger Delta^dagger]^(-1)$ instead of the sub-iteration $G = G_0 (1 + Delta F^dagger)$ with the previous $F$ (which diverged once $Delta$ exceeded about $pi T$). With `sw_self=True` the static part of $Sigma$ is subtracted every iteration (`sub_sigma`; $mu$ is fixed in this loop) and $Sigma$ is mixed with the rate `pp`. A non-finite $Delta$ stops with an error instead of being reported as converged. Converged gaps and $Sigma$ change.
- *Logical kinds on the Fortran side.* Some arguments were passed from Python as 1-byte `c_bool` but received as default 4-byte `logical` (the spin channel `sw_spsym` of the SC-state $chi$, the leading-eigenvalue/trace switch `sw_eig` of $chi$, ...). The remaining 3 bytes are undefined, so they may have been read as wrong values on some platforms. Rerun option=12, 13, 27 if their results look suspicious.

Features added at the same time (earlier results unchanged): the IR frequency grid (`freq_grid='ir'`), the symmetry-free linearized Eliashberg equation (`sw_nosym`), the AiiDA interface, the automatic runtime settings of gfortran builds.
