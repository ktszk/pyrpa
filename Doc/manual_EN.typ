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
- Spin susceptibility $chi_s$ and pairing susceptibility $phi$ based on RPA
- Dynamic spin susceptibility $chi_s^"SC"$ in the superconducting state
- Self-energy calculation via FLEX (Fluctuation EXchange approximation)
- Linearized Eliashberg equation for the superconducting eigenvalue and gap function
- Nonlinear Eliashberg equation for self-consistent SC gap functions
- Carrier number, cyclotron mass, and dHvA frequency calculations
- Spectral function with impurities and via CPA

All settings in pyrpa are controlled by modifying variables at the top of `main.py`, above the library import lines.

= Requirements

The following Python packages are required:

- `numpy`
- `scipy`
- `matplotlib`
- `skimage` (used for Fermi surface plots option=2,3 and cyclotron mass option=18)

The internal Fortran libraries (`libs/flibs`, `libs/plibs`) must be compiled beforehand. The FLEX, linear/nonlinear Eliashberg, and SC chi calculations are executed via OpenMP-parallel Fortran routines.

= Input Files

=== File Format Selection (`ftype`)

The format of the input Hamiltonian file is specified by the `ftype` variable. Assign the file name or directory path (without extension) as a string to `fname`.

#table(
  columns: (auto, 1fr),
  [*ftype*], [*File format*],
  [0], [Directory `{fname}/` containing `ham_r.txt`, `irvec.txt`, and `ndegen.txt`],
  [1], [A single `{fname}.input` file (custom format)],
  [2], [`{fname}_hr.dat` file (Wannier90 default hopping file). *This is the standard choice when interfacing with first-principles codes.*],
  [3], [Non-orthogonal basis hopping in MLO (Maximally Localized Orbital) format],
  [Other], [`Hopping.dat` file (ecalj hopping file)],
)

*Example for Wannier90 users:*

```python
fname = 'inputs/Sr2RuO4'   # path without the _hr.dat suffix
ftype = 2
```

=== Spin-Orbit Coupling (`sw_soc`)

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
  [4], [Trigonal (QE sbrav=5)],
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

  - option=19 (`DHVA`): dHvA frequency vs angle
  - option=20 (`ELECTRON_MASS`): electron mass along the symmetry line
  - option=21 (`SPECTRUM_IMPURITY`): spectral function with impurities
  - option=23 (`NONLIN_ELIASHBERG`): nonlinear Eliashberg equation
]

=== option=0: Band Plot (`BAND`)

Calculates and plots the band dispersion $E_n(bold(k))$ along a symmetry line. The symmetry path is specified via `k_sets` and `xlabel`, or auto-generated from the `brav` setting if not defined.

=== option=1: Density of States (`DOS`)

$ "DOS"(omega) = - 1/pi sum_(bold(k), n) "Im" G^0_n (bold(k), omega + i delta) $

Plots both the total DOS and orbital-projected partial DOS.

=== option=2: 2D Fermi Surface (`FERMI_2D`)

Draws the Fermi surface in the $k_x$-$k_y$ plane at a specified $k_z$ slice (default: $k_z = 0$, adjustable via `kz`). The chemical potential is determined self-consistently from the electron count `fill`, and the contour $E_n(bold(k)) = mu$ is extracted via skimage's `find_contours`. A rotation matrix `RotMat` can optionally be used to rotate the Fermi surface.

=== option=3: 3D Fermi Surface (`FERMI_3D`)

Renders the three-dimensional Fermi surface as a polygon mesh using the Marching Cubes algorithm. The display scale along each axis can be adjusted with `kscale`.

=== option=4: Spectral Function (`SPECTRUM`)

Computes the electron spectral function from the imaginary part of the trace of the single-particle Green's function,

$ A(bold(k), omega) = -1/pi "Im" "Tr" G(bold(k), omega + i delta) $

and plots it along the symmetry line. With `sw_self=True`, the self-energy from a FLEX calculation is incorporated to include interaction effects.

=== option=5: Boltzmann Transport (`CONDUCTIVITY_BT`)

Calculates various transport coefficients using the Boltzmann transport equation under the relaxation time approximation ($tau = $ `tau_const` fs).

Output quantities:
- Electrical conductivity tensor $sigma_{i j}$ (unit: S/m)
- Thermal conductivity tensor $kappa_{i j}$ (unit: W/m/K)
- Seebeck coefficient tensor $S_{i j}$ (unit: V/K)
- Power factor $sigma S^2$ (unit: $upright(W \/ m \/K^2)$)
- Peltier coefficient, Lorenz number

=== option=6: Linear Response Transport (`CONDUCTIVITY_PT`)

Calculates electrical conductivity and related quantities based on the Kubo formula (linear response theory). The parameter `delta` corresponds to an effective relaxation time $tau approx hbar / delta$.

=== option=7: Spin Susceptibility Spectrum (`CHIS_SPECTRUM`)

Calculates the RPA spin susceptibility

$ chi_s (bold(q), omega) = frac(chi^0 (bold(q), omega), 1 - S chi^0 (bold(q), omega)) $

along the symmetry line and plots the imaginary part as a paramagnon spectrum. The interaction matrix $S$ is generated from the on-site parameters `U` and `J`. Results are written to `chis_spec.png`.

=== option=8: Spin Susceptibility at a $bold(q)$ Point (`CHIS_QPOINT`)

Calculates and plots the frequency dependence $chi_s (bold(q), omega)$ at the wave vector $bold(q)$ specified by `at_point`.

=== option=9: $bold(q)$-Space Spin Susceptibility Map (`CHIS_QMAP`)

Plots the spatial distribution of $chi_s (bold(q), omega_0)$ in the $k_x$-$k_y$ plane at the energy $omega_0$ specified by `Ecut`. Useful for visualizing nesting vectors. Output: `chi0map.png` and `chismap.png`.

=== option=10: Pairing Susceptibility Spectrum (`PHI_SPECTRUM`)

Calculates the pairing susceptibility $phi(bold(q), omega)$ along the symmetry line and plots it (`phi_spec.png`).

=== option=11: $bold(q)$-Space Pairing Susceptibility Map (`PHI_QMAP`)

Plots the $bold(q)$-space distribution of the pairing susceptibility at the energy specified by `Ecut` (`phimap.png`). Use `sw_omega` to switch between real and Matsubara frequency.

=== option=12: SC-State Spin Susceptibility Spectrum (`CHIS_SPECTRUM_SC`)

Assumes a non-zero gap function $Delta(bold(k))$ and constructs the irreducible susceptibility in the superconducting state, which includes the anomalous bubble $F(bold(k))$:

$ chi_0^"SC" = chi^"GG" plus.minus chi^"FF" $

The dynamic spin susceptibility $chi_s^"SC" = chi_0^"SC" \/ (1 - S chi_0^"SC")$ is then obtained by RPA and plotted along the symmetry line. The initial gap symmetry is set by `gap_sym` (negative values for triplet). Output: `chis_sc_spec.png`.

The initial gap amplitude `delta0` can be specified in two ways (see also the parameter reference below):
- A single float (e.g. `delta0=1.e-2`): a single-band gap shape is generated internally and scaled to this maximum amplitude across all bands.
- A list of per-band values of length `Norb` (e.g. `delta0=[0.,0.2,0.3,-0.1,0.]`): enables *multi-gap* mode, where each band gets its own amplitude and sign. Mixing signs lets you represent sign-changing gaps such as $s^plus.minus$, as is common in multiband (e.g. Fe-based) superconductors.

=== option=13: SC-State Spin Susceptibility at a $bold(q)$ Point (`CHIS_QPOINT_SC`)

Same SC framework as option=12, but evaluates $chi_s^"SC"(bold(q), omega)$ at the single $bold(q)$ specified by `at_point`. `delta0` is specified the same way as in option=12.

This mode additionally produces:

- The BdG (Bogoliubov–de Gennes) band dispersion along $bold(k) = (0,0,0) arrow.r (0,0.5,0)$ (`BdG_band.png`)
- The orbital-traced $chi_s^"SC"(omega)$ (`chisq.png`, `chis_sc.dat`)
- The orbital-resolved $chi_s^"SC"(omega)$ (`chisq_orb.png`, `chis_scorb.dat`)

=== option=14: FLEX Self-Energy (`FLEX`)

Solves the FLEX equations self-consistently to obtain the electron self-energy $Sigma(bold(k), i omega_n)$ in Matsubara frequency space. The result is saved to `sigma.bin` and `self_en.npz`.

- `sw_out_self=True`: write the self-energy to file
- `sw_in_self=True`: load the previous self-energy (`sigma.bin`) as the initial guess
- `m_diis_num`: DIIS history length (default 5 if undefined)
- `sw_rescale_flex=True`: dynamically rescale the self-energy so that max$|Sigma| approx U$, useful when Stoner factor is close to 1

=== option=15: Linearized Eliashberg Equation (`LIN_ELIASHBERG`)

Solves the linearized Eliashberg equation as an eigenvalue problem using the effective pairing interaction from RPA or FLEX. The largest eigenvalue $lambda$ and the corresponding gap function $Delta(bold(k), i omega_n)$ are obtained. The superconducting transition temperature $T_c$ corresponds to $lambda = 1$. The default solver is a power method with shift + deflation; setting `arnoldi_m > 0` enables an Arnoldi solver.

- `sw_self=False`: no self-energy corrections (pure RPA)
- `sw_self=True`: uses the FLEX self-energy (requires `sigma.bin`)
- `sw_from_file=True`: loads the self-energy from `sigma.bin` without re-running FLEX
- `gap_sym`: sets the initial symmetry of the gap function (see below)
- Gap functions are written to `gap_{ij}.dat` and `gap.npy`

=== option=16: Gap Function Post-Processing (`GAP_FUNCTION`)

Reads the gap function from `gap.npy` (produced by option=15) and computes:

- The anomalous Green's function $F(bold(k), i omega_n)$
- Analytic continuation to the real-frequency axis via Padé approximation

=== option=17: Carrier Number (`CARRIER_NUM`)

Calculates the electron (particle) and hole counts for each band from the volume of the Fermi surface. Useful for checking consistency with `fill`.

=== option=18: Cyclotron Mass (`CYCLOTRON_MASS`)

Computes the cyclotron mass via the Onsager relation $m^*_c = (planck^2 \/ 2 pi)(partial S \/ partial E)$, where $S$ is the Fermi surface cross section.

- Phase 1: Scan $S(k_z)$ across $k_z in [0, pi/2]$ with `meshkz=20` points
- Phase 2: Refine extremal $k_z$, compute $partial S \/ partial E$ via central finite difference, and report $m^*_c$ in units of $m_e$

Internally uses skimage's `find_contours` and `marching_cubes`.

=== option=19: dHvA Frequency vs Angle (`DHVA`)

Computes the dHvA frequency $F(theta)$ as a function of the polar angle $theta$ via the Onsager relation $F = planck \/ (2 pi e) dot A_"ext"$. The angle list scans 0–90° with 40 points by default.

=== option=20: Electron Mass (`ELECTRON_MASS`)

Computes $m^* = planck^2 (partial^2 E \/ partial bold(k)^2)^{-1}$ along the symmetry line (in units of $m_e$).

=== option=21: Spectrum with Impurities (`SPECTRUM_IMPURITY`)

Builds a real-space supercell Hamiltonian containing impurities and computes the spectral function $A(bold(k), omega)$ via the impurity Green's function. The impurity sites are specified in `imp_list`.

=== option=22: CPA Conductivity / Spectrum (`SIGMA_CPA`)

Solves the Coherent Potential Approximation (CPA) self-consistency for an alloy and outputs both the real-frequency spectral function (`cpa_spectrum.png`) and the Matsubara self-energy (`sigma.bin`, `self_en.npz`). The impurity concentration `x_cpa` and the on-site perturbations `VA`, `VB` are currently fixed inside the routine.

=== option=23: Nonlinear Eliashberg Equation (`NONLIN_ELIASHBERG`)

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

= Color Plot Settings (`color_option`)

For option=0, 2, and 3, each point on the band or Fermi surface can be colored according to a physical quantity.

#table(
  columns: (auto, 1fr),
  [*color_option*], [*Meaning*],
  [0], [No color (black)],
  [1], [Orbital weights specified by `olist` are mapped to RGB (red/green/blue)],
  [2], [Group velocity magnitude $|bold(v)(bold(k))|$ is shown as a color gradient],
)

For `color_option=1`, specify orbital indices in `olist` as `[R component, G component, B component]`. To assign multiple orbitals to the same color, use a nested list. Example:

```python
olist = [[0, 3], [1, 4], [2, 5]]
# Orbitals 0 and 3 → Red, orbitals 1 and 4 → Green, orbitals 2 and 5 → Blue
```

= Parameter Reference

This section explains all parameters in the upper part of `main.py`, including their physical meaning.

=== Basic Settings

- `fname` (string): Path to the input file. The format depends on `ftype`; the file extension is usually omitted.

- `ftype` (integer): Input Hamiltonian file format (see Section 3).

- `brav` (integer): Bravais lattice type (see Section 4).

- `sw_soc` (bool): Switch to enable spin-orbit coupling.

=== Mesh Settings

- `Nx, Ny, Nz` (integer): Number of $bold(k)$-point mesh divisions in the first Brillouin zone along $x$, $y$, $z$. Used for 3D $bold(k)$-space integrations and FFT-based convolutions (FLEX, Eliashberg, chi). For 2D systems, set `Nz=1`. Powers of 2 (32, 64, ...) are preferred for FFT efficiency. Memory consumption scales as $tilde.equiv N_x N_y N_z dot N_w dot N_"orb"^2$.

- `Nw` (integer): Number of Matsubara frequencies for FLEX/Eliashberg/SC-chi calculations (option=12,13,14,15,16,23), or number of real-frequency points for DOS/spectral function calculations (option=1,4, etc.). Matsubara frequencies are $omega_n = (2n+1) pi T$ for $n = 0, 1, \ldots, N_w - 1$. Lower temperatures require larger $N_w$ (typical: 256–1024).

- `kmesh` (integer): Number of $bold(k)$ points along the symmetry line for band and spectral function plots. Larger values yield smoother plots (200–500 is typical).

=== Lattice Constants

- `abc` (list, unit: Å): Lattice constants $a, b, c$. Used for computing group velocities $v = 1/ hbar (partial E)/(partial bold(k))$ and the physical length scale of symmetry paths. For Wannier90, use values consistent with the WIN/WOUT files.

- `alpha_beta_gamma` (list, unit: degrees): Lattice angles $alpha, beta, gamma$. For orthorhombic and cubic systems, set to `[90., 90., 90.]`.

=== Temperature and Chemical Potential

- `tempK` (float, unit: K): Temperature in Kelvin. Converted internally to $T = k_B$ `tempK` (in eV).

- `temp` (float, unit: eV): Directly specifies $k_B T$ in eV. If both `tempK` and `temp` are defined, `temp` takes precedence.

- `fill` (float): Band filling. The chemical potential $mu$ is determined self-consistently from $sum_(bold(k), n) f(epsilon_n (bold(k)) - mu) = N_k dot f"ill"$. Allowed range is 0 to `no` (the number of orbitals in the Hamiltonian); `fill = no` corresponds to full filling. The mapping to the physical electron count per unit cell depends on the Hamiltonian convention:
  - Without SOC (`sw_soc=False`): the spin degree of freedom is not included explicitly, so `fill` represents the electron occupation per spin. For example, a 3-orbital model is full-filled at `fill = 3` (6 electrons including both spins) and half-filled at `fill = 1.5`.
  - With SOC (`sw_soc=True`): the Hamiltonian already includes spin (`no = 2 N_"orb"`), so `fill` directly equals the total electron count per unit cell. A 3-orbital SOC model (`no = 6`) is full-filled at `fill = 6` and half-filled at `fill = 3`.

- `mu0` (float, unit: eV): If defined, this value is used directly as the chemical potential, bypassing the self-consistent calculation from `fill`.

=== Energy Range and Broadening

- `Emin, Emax` (float, unit: eV): Lower and upper bounds of the energy range for DOS and spectral function calculations.

- `delta` (float, unit: eV): Broadening parameter $delta$ for spectral calculations. This is the small imaginary part added to the Green's function, which broadens the Dirac delta function into a Lorentzian of finite width. Too large a value smears out features; too small a value introduces numerical noise (a typical range is 0.01–0.05 eV).

- `Ecut` (float, unit: eV): Fixed energy $omega_0$ for $bold(q)$-space susceptibility maps (option=9,11). Set near zero to probe the Fermi surface region.

- `delta0` (float, or a list of length `Norb`, unit: eV): Amplitude of the initial gap function for SC-chi calculations (option=12,13). Physically corresponds to the SC gap size (typical: $10^{-3}$–$10^{-2}$ eV ≈ 1–10 meV).
  - As a float: a single internally-generated gap shape (common to all bands) is scaled to this maximum amplitude. Setting it to 0 reduces to the normal-state calculation.
  - As a list (e.g. `delta0=[0.,0.2,0.3,-0.1,0.]`): enables multi-gap mode, where each band's amplitude and sign are set independently — use this to represent sign-changing gaps such as $s^plus.minus$.

=== Transport Parameters

- `tau_const` (float, unit: fs): Constant relaxation time $tau$ for Boltzmann theory (option=5). In the constant relaxation time approximation, this is a free parameter typically chosen by comparison with experiment (typical metals: 1–100 fs).

- `sw_tdf` (bool): If `True`, the transport distribution function (TDF) is computed first, and transport coefficients are obtained by energy integration. Relevant when using an energy-dependent relaxation time.

=== Orbital and Interaction Parameters

- `olist` (list): Orbital indices for color plotting (`color_option=1`); see Section 6.

- `U` (float, unit: eV): On-site Coulomb repulsion (Hubbard $U$). Used in FLEX/RPA calculations. This is a key parameter controlling magnetic and superconducting instabilities.

- `J` (float, unit: eV): On-site Hund's coupling constant. The screened interaction $U' = U - 2J$ (Kanamori screening) is used automatically.

- `orb_dep` (bool): If `True`, orbital-dependent interaction matrices `Umat` and `Jmat` are used. If `False` (default), the constant values `U` and `J` are applied uniformly to all orbitals.

- `m_diis_num` (integer, optional): DIIS (Pulay-accelerated mixing) history length for FLEX (option=14) and nonlinear Eliashberg (option=23). Values $>= 2$ enable Pulay extrapolation; $1$ falls back to linear mixing. Defaults to 5 if undefined. Larger values speed up convergence at the cost of memory ($N_x N_y N_z dot N_w dot N_"orb"^2$ per slot).

=== Initial Gap Function Symmetry (`gap_sym`)

Specifies the symmetry of the initial gap function when solving the Eliashberg equation (option=15,23), or when generating the initial gap shape for the SC-chi calculations (option=12,13).

#table(
  columns: (auto, 1fr),
  [*gap_sym*], [*Symmetry*],
  [0], [$s$-wave (uniform positive sign for all $bold(k)$)],
  [1], [$d_{x^2-y^2}$-wave ($cos k_x - cos k_y$ type)],
  [2], [$s^plus.minus$-wave (sign changes across the nesting vector)],
  [3], [$d_{x y}$-wave ($sin k_x sin k_y$ type)],
  [-1], [$p_x$-wave],
  [-2], [$p_y$-wave],
)

=== $bold(k)$-Space Settings

- `kz` (float, in reduced coordinates): The $k_z$ value for the 2D Fermi surface plot (option=2). Ranges from 0 to 0.5 ($k_z = 0$: $Gamma$-plane, $k_z = 0.5$: zone boundary plane).

- `kscale` (list or float): Display scale for each axis in the 3D Fermi surface plot (option=3). Example: `kscale=[1.0, 1.0, 0.5]` compresses the $k_z$ direction by half.

- `k_sets` (list of lists): Coordinates of symmetry line endpoints in reduced units (0 to 1). Define this to specify a custom symmetry path instead of the auto-generated one. Example: `k_sets=[[0,0,0],[0.5,0,0],[0.5,0.5,0]]`

- `xlabel` (list of strings): Labels for the points in `k_sets`. Example: `xlabel=[r'$\Gamma$','X','M']`

- `at_point` (list): Coordinates of the $bold(q)$ point (in reduced units) for the single-$bold(q)$ susceptibility calculation (option=8).

=== Switch Variables

- `sw_unit` (bool): If `True` (default), physical constants in SI units are used and output is in physical units. If `False`, a dimensionless system with $hbar = k_upright(B) = e = 1$ is used.

- `sw_omega` (bool): Switch for option=11 — compute the pairing susceptibility on the real frequency axis (`True`) or Matsubara frequency axis (`False`).

- `sw_self` (bool): If `True`, the FLEX self-energy is incorporated into option=4 (spectrum), option=15 (linear Eliashberg), or option=23 (nonlinear Eliashberg) via the renormalized Green's function.

- `sw_out_self` (bool): If `True`, the FLEX self-energy is written to `sigma.bin` and `self_en.npz`. Also used to trigger writing the gap function in option=15.

- `sw_in_self` (bool): If `True`, the previous FLEX self-energy is loaded from `sigma.bin` as the initial guess for the iterative self-consistent loop.

- `sw_from_file` (bool): If `True`, the self-energy is read from `sigma.bin` and the FLEX calculation is skipped. The Eliashberg equation is then solved directly with this pre-computed self-energy.

- `sw_check_only` (bool): Used only by option=23 (nonlinear Eliashberg). If `True`, the routine stops right after the internal linearized-Eliashberg solve (reporting the Stoner factor $S$ and eigenvalue $lambda_"eliash"$) without running the nonlinear loop — handy for quickly bracketing $T_c$ via a temperature scan. Regardless of this flag, the nonlinear loop is also skipped automatically whenever $S >= 1$ or $lambda_"eliash" < 1$.

- `sw_rescale_flex` (bool): For FLEX (option=14), dynamically rescale the self-energy when max$|Sigma|$ approaches `U` to prevent divergence. Useful when the Stoner factor is close to 1.

- `sw_dec_axis` (bool): If `True`, lattice vectors are decomposed appropriately to set up the reciprocal lattice vectors.

=== Eilenberger Parameters (option=24,25,26)

These drive the quasiclassical Eilenberger solvers. The temperature is the global `tempK`/`temp` and the gap symmetry is the global `gap_sym`.

*Common (all three modes):*

- `eil_coupling` (float): the dimensionless separable pairing coupling $lambda$ (with $⟨ |phi|^2 ⟩_"FS" = 1$). Larger $lambda$ → higher $T_c$.
- `eil_wc` (float, unit: eV): the fixed Matsubara cutoff energy, which sets the pairing scale / $T_c$.
- `eil_fs_kind` (`None`/`'iso'`/`'ellipse'`/`'tb'`/`'wannier'`): the Fermi surface. `None` = isotropic cylinder (the homogeneous penetration calc falls back to `'ellipse'`); `'iso'`/`'ellipse'`/`'tb'` = model FS from `eil_fs_params`; `'wannier'` = the real FS + Fermi velocities of the loaded band (gap symmetry/multiband from `gap_sym`, `delta0`, `eil_gap_orbital`).
- `eil_fs_params` (tuple): model-FS parameters — ellipse masses $(m_x, m_y)$ or the `tb` hopping.
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
- `eil_vort_dvector` (bool): the self-consistent triplet d-vector vortex/lattice texture.
- `eil_vort_tilt` (float, unit: deg): the field tilt from the $c$-axis (quasi-2D: orbital $B_z = B cos theta$, Zeeman $-> h\/cos theta$).
- `eil_gap_orbital` (`None` / $N_"orb" times N_"orb"$ matrix or callable): an orbital-basis pair potential whose low-energy projection onto the FS bands sets the gap (Nagai–Nakamura, JPSJ *85*, 074707 (2016), Eq. 43; needs a Wannier FS), superseding `gap_sym`/`delta0`.
- `eil_gap_file` (`None` / string): the base name (no extension) of a self-consistent RPA/FLEX gap exported as a Wannier-real-space "hopping" file by option=15/23 (`LIN_ELIASHBERG`/`NONLIN_ELIASHBERG`) with `sw_out_self=True` (`output_gap_wannier`, e.g. `'gap_wannier'`). When set, $Delta(bold(R), i omega_n)$ is loaded and used as `eil_gap_orbital` — its inverse Fourier transform $Delta_"orb"(bold(k)) = sum_bold(R) e^(i 2 pi bold(k) dot bold(R)) Delta(bold(R))$ is projected onto the FS bands. This is the route to use a *previously computed RPA gap* (e.g. for $"KFe"_2"As"_2$, PRB *84*, 144514) as the vortex pairing form factor. The exporting RPA run and the Eilenberger run must use the *same* Wannier Hamiltonian (same orbital basis, $bold(R)$/$bold(a)$ convention, and ideally $mu$/filling) so that the band eigenvectors and $Delta_"orb"$ share a basis. Supersedes `eil_gap_orbital`/`gap_sym`.
- `eil_gap_iw` (int): the starting Matsubara index for `eil_gap_file` ($0$ = lowest $i omega_0$). The Eilenberger form factor is static; $i omega_0$ carries the symmetry / sign / node / anisotropy structure most sharply and matches the gap usually quoted on the Fermi surface.
- `eil_gap_navg` (int): number of consecutive Matsubara slices averaged for `eil_gap_file` ($1$ = single $i omega_("eil_gap_iw")$ slice). $> 1$ smooths slice noise at the cost of slightly diluting the anisotropy (since $Delta(bold(k), i omega_n)$ gets more isotropic with $n$). The absolute scale is irrelevant — the projected $phi$ is renormalized to $⟨ |phi|^2 ⟩ = 1$.

= Typical Calculation Workflows

=== Checking Band Structure and Fermi Surface

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

=== Superconducting Gap Function (RPA)

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

=== FLEX + Linearized Eliashberg

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

=== Nonlinear Eliashberg (Self-Consistent SC Loop)

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

=== Dynamic Spin Susceptibility in the SC State

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

The `tests/` directory contains regression tests for the main numerical kernels and physics benchmarks. Each test file can be run directly as a Python script, so `pytest` is optional.

Before running the tests, the Fortran shared library `libs/libfmod.so` must be compiled. If it is missing, enter the `libs` directory and run `make FC=<compiler> SL=<library>`.

=== How to Run

Run individual test files directly:

```bash
python tests/test_eilenberger.py
python tests/test_rpa_flex.py
```

If `pytest` is available, the whole test directory can also be run with:

```bash
pytest tests
```

=== `tests/test_eilenberger.py`

This file tests the quasiclassical Eilenberger / Riccati solvers. It covers homogeneous systems, surfaces, vortices, vortex lattices, model Fermi surfaces, Pauli limiting, and triplet $d$-vector textures.

Main checks:

- Matsubara cutoff scaling with temperature
- Anderson theorem and Abrikosov--Gor'kov pair breaking
- Weak-coupling BCS ratio $2 Delta_0 / k_B T_c approx 3.53$
- Agreement between Fortran Riccati kernels and Python reference implementations
  - scalar `riccati_chords`
  - spin-matrix `matrix_riccati_batch`
  - batched chord `matrix_riccati_chords`
- Gap suppression and zero-energy bound states at a $d$-wave surface
- CdGM zero-energy peak at a vortex core
- Volovik-like field dependence in a vortex lattice
- Normalization of `build_model_fs` and the isotropic Fermi-surface limit
- Singlet Pauli suppression and robustness of triplet equal-spin pairing
- Triplet $d$-vector textures near surfaces and vortex cores

This is a physics-oriented benchmark suite and may take several tens of seconds depending on the machine.

=== `tests/test_rpa_flex.py`

This file provides lightweight regression tests for the RPA / FLEX / Eliashberg building blocks outside the Eilenberger suite. It targets Fortran wrappers, RPA matrix algebra, FLEX bubble / vertex kernels, and Eliashberg smoke tests.

Main checks:

- Consistency of multi-site orbital-pair and site-index generation in `get_chi_orb_list`
- Reference values of the two-orbital Kanamori-type vertex from `gen_SCmatrix`
- Reference values of orbital-dependent `Umat`, `Jmat` vertices from `gen_SCmatrix_orb`
- One-orbital RPA formula
  $ chi_s = chi^0 / (1 - U chi^0) $
- `get_chis_chic` reduces to bare $chi^0$ when `S=C=0`
- Analytic one-orbital Green's function check for `gen_Green0`:
  $ G^0(k,i omega_n) = 1 / (i omega_n + mu - epsilon_k) $
- Agreement between `get_chi0` and `get_chi0_conv`
- `get_Vsigma_nosoc_flex` returns a zero vertex at zero interaction
- `linearized_eliashberg` returns $lambda = 0$ and finite arrays at zero interaction
- `nonlinear_eliashberg` preserves the trivial solution $Delta=0$ for zero seed and zero interaction
- `_load_sigma_from_file` returns `None` without crashing when `self_en.npz` is absent
- One-orbital smoke test for `output_gap_function`

This test uses a very small one-orbital model and a small $k$ mesh, so it is intended to catch RPA/FLEX API regressions quickly.

= Troubleshooting

- *Chemical potential does not converge*: Increase `Nx, Ny, Nz`, or slightly raise `tempK`. A coarse $bold(k)$-mesh can cause instability in the self-consistent $mu$ search.

- *Too much noise in the spectrum*: Increase `delta` slightly, or increase `Nw` to improve the frequency resolution.

- *Band structure looks wrong*: Check that `brav` is set correctly. In particular, when interfacing with QE, pay careful attention to the conventions for FCC (`brav=1`) and BCC (`brav=2`).

- *Linear Eliashberg does not converge*: Try reducing `U`, increasing `tempK`, or using a finer mesh (`Nx, Ny, Nz`). If the eigenvalue $lambda$ exceeds 1, the system is inside the superconducting phase at that temperature.

- *Nonlinear Eliashberg diverges*: Increase `m_diis_num` to 5–10, raise `tempK`, or increase `Nw`. Running with `sw_check_only=True` first to inspect the Stoner factor $S$ and eigenvalue $lambda_"eliash"$ is also useful. If $S$ exceeds 1, the system is magnetically unstable and `U` must be reduced; if $lambda_"eliash" < 1$, you are still above $T_c$ and need to lower the temperature (in either case the nonlinear loop is skipped automatically).

- *FLEX self-energy diverges*: Set `sw_rescale_flex=True` or reduce `U`.

- *Transport coefficients have wrong units*: Make sure `sw_unit=True`. Setting it to `False` switches to a dimensionless unit system.
