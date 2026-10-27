#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
ftype: set input hamiltonian's format
0: ham_r.txt, irvec.txt, ndegen.txt in {fname} dir
1: .input file named {fname}
2: {fname}_hr.dat file (wannier90 default hopping file
3: for non-orthogonal basehoppings as MLO basis
4: HamRsMLO, HamiltonianPMTInfo binaries in {fname} dir (ecalj job_mlo output)
else: Hopping.dat file (ecalj hopping file)
brav: choose primitive translation vector S,FC,BC etc
0: simple 
1: face center (for QE default)
2: body center (for QE default)
3: hexagonal
4: trigonal    (for QE ibrav==5)
5: base center
6: face center (common)
7: body center (common)
else: monoclinic
"""
import enum

#fname, ftype, brav, sw_soc = 'inputs/Cu', 4, 6, False        # ecalj job_mlo binaries (abc=[3.597]*3, fill=5.5)
fname, ftype, brav, sw_soc = 'inputs/00010.input', 1, 2, False        # .input hopping file, body-centred cell
sw_dec_axis = False


class CalcMode(enum.IntEnum):
    """Calculation modes. The second tuple element is the human-readable label
    printed by main(); it is the single source of truth (replaces the old `opstr` list)."""
    def __new__(cls, value: int, description: str) -> 'CalcMode':
        obj = int.__new__(cls, value)
        obj._value_ = value
        obj.description = description
        return obj
    BAND               = (0,  "calculate band structure")                  # band structure plot
    DOS                = (1,  "calculate Dos")                              # density of states plot
    FERMI_2D           = (2,  "plot 2D Fermi surface")                     # 2D Fermi surface at kz plane (default: kz=0)
    FERMI_3D           = (3,  "plot 3D Fermi surface")                     # 3D Fermi surface plot
    SPECTRUM           = (4,  "calculate spectrum")                        # spectral function plot
    CONDUCTIVITY_BT    = (5,  "calculate conductivities using Boltzmann theory")  # conductivity via Boltzmann theory
    CONDUCTIVITY_PT    = (6,  "calculate conductivities with linear response")    # optical conductivity via linear response theory
    CHIS_SPECTRUM      = (7,  "calculate chis spectrum")                   # spin susceptibility spectrum along symmetry line
    CHIS_QPOINT        = (8,  "calculate chis at q-point")                 # spin susceptibility at specified q-point
    CHIS_QMAP          = (9,  "calculate chis qmap at Ecut")               # spin susceptibility q-map on Ecut plane
    PHI_SPECTRUM       = (10, "calc phi spectrum with symmetry line")      # pairing susceptibility spectrum along symmetry line
    PHI_QMAP           = (11, "calc phi on xy plane at Ecut")              # pairing susceptibility q-map on Ecut plane
    CHIS_SPECTRUM_SC   = (12, "calculate chis spectrum on sc")             # spin susceptibility spectrum in superconducting state
    CHIS_QPOINT_SC     = (13, "calculate chis at q-point on sc")           # spin susceptibility at specified q-point in superconducting state
    FLEX               = (14, "calc self energy")                          # self-energy calculation using FLEX
    LIN_ELIASHBERG     = (15, "solve linearized eliashberg equation")      # solve linearized Eliashberg equation
    GAP_FUNCTION       = (16, "gap_function")                              # post-process and output gap functions
    BAND_FILLING       = (17, "check band filling and electron/hole character")  # per-band filling & FS character
    CYCLOTRON_MASS     = (18, "calculate cyclotron mass")                  # cyclotron mass calculation
    DHVA               = (19, "plot dHvA frequency")
    ELECTRON_MASS      = (20, "calculate electron mass")                   # band (effective) mass tensor along the symmetry line
    SPECTRUM_IMPURITY  = (21, "spectrum with impurity")                    # unfolded spectral function of an impurity supercell
    SIGMA_CPA          = (22, "calculate sigma_cpa")                       # conductivity via CPA
    NONLIN_ELIASHBERG  = (23, "solve nonlinear eliashberg equation")       # solve nonlinear Eliashberg equation (no SOC)
    EILENBERGER        = (24, "solve quasiclassical Eilenberger equation")  # homogeneous multi-orbital Eilenberger (Matsubara)
    EILENBERGER_SURFACE = (25, "solve surface state via Riccati Eilenberger")  # specular surface gap profile + LDOS (model FS)
    EILENBERGER_VORTEX = (26, "solve isolated vortex via Riccati Eilenberger")  # vortex D(rho) + core LDOS (model FS)
    NMR_SC             = (27, "calculate NMR observables (K, 1/T1T) in the sc state")  # Knight shift & spin-lattice rate from the BdG chi_s


class ColorMode(enum.IntEnum):
    """Color modes for band/FS plots; second element is the printed label (replaces the old `cstr` list)."""
    def __new__(cls, value: int, description: str) -> 'ColorMode':
        obj = int.__new__(cls, value)
        obj._value_ = value
        obj.description = description
        return obj
    MONO     = (0, "no color")
    ORBITAL  = (1, "orbital weight")
    VELOCITY = (2, "velocity size")
    GAP      = (3, "Eilenberger pairing gap Re[phi(k)] (option 3 only)")


# ----- CalcMode capability groups: single source for the option-set checks in main() -----
M = CalcMode
# modes that need a symmetry-line k-path (k_sets / xlabel)
MODES_SYMLINE       = frozenset({M.BAND, M.SPECTRUM, M.CHIS_SPECTRUM, M.PHI_SPECTRUM, M.CHIS_SPECTRUM_SC,
                                 M.FLEX, M.ELECTRON_MASS, M.SPECTRUM_IMPURITY, M.SIGMA_CPA})
# modes that need a rotation matrix (RotMat) for the 2D / extremal-orbit cut
MODES_NEED_ROTMAT   = frozenset({M.FERMI_2D, M.CYCLOTRON_MASS})
# FLEX / Eliashberg family (self-energy, linearized gap, post gap)
MODES_FLEX_ELIASH   = frozenset({M.FLEX, M.LIN_ELIASHBERG, M.GAP_FUNCTION})
# chi at a single q-point (normal & SC state)
MODES_CHIS_QPOINT   = frozenset({M.CHIS_QPOINT, M.CHIS_QPOINT_SC})
# modes that color-code bands/FS by orbital or velocity
MODES_COLOR         = frozenset({M.BAND, M.FERMI_2D, M.FERMI_3D})
# susceptibility / Eliashberg modes that build the orbital-pair basis (chiolist)
MODES_NEED_CHI      = frozenset({M.CHIS_SPECTRUM, M.CHIS_QPOINT, M.CHIS_QMAP, M.PHI_SPECTRUM, M.PHI_QMAP,
                                 M.CHIS_SPECTRUM_SC, M.CHIS_QPOINT_SC, M.NMR_SC, M.FLEX, M.LIN_ELIASHBERG, M.NONLIN_ELIASHBERG})
# susceptibility modes that need the static Coulomb vertex (S/C matrices)
MODES_COULOMB_VERTEX = frozenset({M.CHIS_SPECTRUM, M.CHIS_QPOINT, M.CHIS_QMAP, M.CHIS_SPECTRUM_SC, M.CHIS_QPOINT_SC, M.NMR_SC})
# modes that print a constant U,J header
MODES_PRINT_UJ      = frozenset({M.CHIS_SPECTRUM, M.CHIS_QPOINT, M.CHIS_QMAP, M.FLEX, M.LIN_ELIASHBERG})
# modes that compute/define mu themselves (skip the common mu block)
MODES_SELF_MU       = frozenset({M.CONDUCTIVITY_BT, M.CONDUCTIVITY_PT, M.SPECTRUM_IMPURITY, M.SIGMA_CPA, M.EILENBERGER_SURFACE, M.EILENBERGER_VORTEX})
# k-mesh reporting groups
MODES_KMESH_SYMLINE = frozenset({M.BAND, M.SPECTRUM})
MODES_KMESH_SINGLE  = frozenset({M.FERMI_2D, M.FERMI_3D, M.CHIS_QMAP, M.PHI_QMAP, M.BAND_FILLING, M.CYCLOTRON_MASS, M.EILENBERGER})
# modes that may run with sw_dec_axis: the axis decomposition is a Fermi-surface DRAWING
# aid (unfolding a body-centred 122 cell into the conventional axes, so that the 10-orbital
# 2-Fe FS can be laid over the 5-orbital 1-Fe one).  The half-integer R components it
# leaves behind are exactly what that unfolding needs, but they also mean H(k) is no longer
# periodic in the decomposed fractional k, so anything that sums over the BZ, FFT-convolves,
# folds an irreducible wedge or divides by the cell volume is invalid -- see the guard in
# main().  The band plot is allowed for the same comparison: it only evaluates E(k) on a
# path, with no BZ sum (the DOS, option=1, IS a BZ sum and stays out).
MODES_DEC_AXIS      = frozenset({M.BAND, M.FERMI_2D, M.FERMI_3D})
MODES_MATSUBARA     = frozenset({M.FLEX, M.LIN_ELIASHBERG})
MODES_IR            = frozenset({M.FLEX, M.LIN_ELIASHBERG, M.GAP_FUNCTION, M.NONLIN_ELIASHBERG})  # may run on freq_grid='ir'
# dispatch groups
MODES_CHI_NORMAL    = frozenset({M.CHIS_SPECTRUM, M.CHIS_QPOINT, M.CHIS_QMAP, M.PHI_SPECTRUM, M.PHI_QMAP})
# SC-state chi family: all share the BdG gap/hamk construction block
MODES_CHIS_SC       = frozenset({M.CHIS_SPECTRUM_SC, M.CHIS_QPOINT_SC, M.NMR_SC})
# modes that reconstruct Delta(-k) from Delta(k) with the inversion operator (plist/iperm).
# For these the operator must really be a symmetry of H(R) -- see plibs.inversion_op.
# (FLEX alone never uses it; LIN_ELIASHBERG with sw_nosym is exempted in main().)
MODES_NEED_PARITY   = frozenset({M.LIN_ELIASHBERG, M.NONLIN_ELIASHBERG, M.GAP_FUNCTION})
del M

# option=CalcMode.CHIS_QPOINT_SC
option = CalcMode.FERMI_3D  # calculation mode to run (see the CalcMode enum above; 0-23 RPA/FLEX/transport, 24-26 Eilenberger superconductivity, 27 NMR in the SC state)
color_option = ColorMode.ORBITAL  # band/FS coloring (option 0,2,3): MONO=black, ORBITAL=olist weights->RGB, VELOCITY=|v(k)|, GAP=Re[phi(k)] from gap_sym/delta0/eil_gap_orbital/eil_gap_file (option 3 only) -- check the Eilenberger pairing gap on the real 3D FS

# Nx,Ny,Nz,Nw=256,256,4,200 #k and energy(or matsubara freq.) mesh size
Nx, Ny, Nz, Nw = 64, 64, 64, 2500
kmesh = 200  # number of k-points along the symmetry line for band/spectrum plots (larger=smoother)
kscale = [1.0, 1.0, 1.0]  # per-axis display scale for the 3D Fermi-surface plot (option 3); e.g. [1,1,0.5] compresses kz
kz = 0.0  # reduced kz of the 2D Fermi-surface cut (option 2): 0=Gamma plane, 0.5=zone-boundary plane
# RotMat=[[0,0,1],[0,1,0],[1,0,0]]

abc=[3.96,3.96,13.02]
#abc = [3.597, 3.597, 3.597]  # lattice constants a,b,c [Angstrom] (group velocities & symmetry-path lengths)
# alpha_beta_gamma=[90.,90.,90]  #lattice angles alpha,beta,gamma [deg] (default 90,90,90 if undefined)
# temp=2.0e-2 #2.59e-2   #directly set k_B*T [eV]; if defined it overrides tempK
tempK = 200  # Kelvin        #temperature [K] (converted internally to temp=k_B*tempK [eV])
fill = 5.5  # band filling; mu solved from sum f(eps-mu)=Nk*fill (no SOC: per spin, full=Norb; SOC: total, full=2*Norb)
# site_prof=[5]

inv_tol = 3.0e-2  # relative residual below which the inversion operator is accepted as a
# symmetry of H(R) (option 15,16,23). Fitted first-principles models
# carry a few % of error, so this is deliberately loose; see plibs.inversion_op
Emin, Emax = 0, 25.  # energy window [eV] for DOS / spectral-function plots (option 1,4)
delta = 2.72e-1  # spectral broadening eta [eV]: imaginary part added to G (Lorentzian width); too large smears, too small=noise
Ecut = 1.0e-2  # fixed energy omega_0 [eV] for the q-space susceptibility maps (option 9,11); ~0 probes the Fermi surface
tau_const = 100  # constant relaxation time tau [fs] for Boltzmann transport (option 5)
# --- colour of the optical response (option 6) ---
# The metallic lustre needs no parameter (it is R(w) itself). A non-metal is coloured by the
# light that got THROUGH the material, so it needs a length: color_thick is the slab thickness
# for the transmission colour, color_scatt the effective scattering length of the Kubelka-Munk
# body colour of a powder/pigment (its single free parameter, ~ the particle size).
color_thick = 1.0e-3  # sample thickness d [m] for the transmission colour (1e-3 = 1 mm)
color_scatt = 1.0e-6  # effective scattering length l [m] for the Kubelka-Munk body colour
# olist=[0,0,0]          #orbital indices mapped to R,G,B for orbital-weight coloring (color_option=1); nested lists group orbitals
olist = [0, [1, 2], 3]
# U,J=0.,0.
# U,J= 0.2, 0.025
# U,J= 0.4, 0.05
U, J = 0.6, 0.075
# U,J=1.2,0.15           #on-site Hubbard U and Hund J [eV] (FLEX/RPA); screened U'=U-2J used automatically
# U,J=1.8,0.225
# 0:s,1:dx2-y2,2:spm,3:dxy,-1:px,-2:py,-3:p+ip  (also drives ALL eilenberger modes; model FS/cylinder maps the int -> continuum harmonic, 2 spm -> s)
gap_sym = 2
# gap_rcut=3.825  # |R| [A] naming the shell of the s+- harmonic (gap_sym=2 only; None=seed-pattern
# cutoff of _CART_SEED).  Every shell carries A1g, so for s+- the shell is a physical
# choice and the seed-pattern stand-in is basis dependent (a centred cell inflates it
# and a far shell is picked).  For KFe2As2 use the crystal axis a=3.825: it gives the
# SAME harmonic 2cos(k'_x a_Fe)cos(k'_y a_Fe) in the 1-Fe (5-orbital) and 2-Fe
# (10-orbital) models -- the one that changes sign between Gamma and the corner pocket.

# use calculation of magnetic susceptibility at superconducting state
# delta0=1.e-2 #maximum gap size for calculating susceptibility in SC state; set to 0 for normal state
d0 = 1.e-1  # helper amplitude scale [eV] for building the per-band delta0 list below
delta0 = [0., d0*2., d0*3., -d0, 0.]  # initial SC gap amplitude/sign per band [eV] for SC-state chi (option 12,13,27); float=single shape, list=multi-gap (signs->s+-)

# --- NMR in the SC state (option 27): Knight shift K(T) and 1/(T1 T) from the BdG chi_s ---
# delta0 above is read as Delta(0); the gap SHAPE stays fixed and only its amplitude is
# scaled by the BCS interpolation Delta(T)=Delta(0)*tanh(1.74*sqrt(Tc/T-1)).
nmr_tc = None  # Tc [eV]; None -> weak-coupling estimate max|delta0|/1.764
nmr_trange = [0.05, 1.2]  # reduced temperature range [T/Tc min, T/Tc max] of the sweep
nmr_nt = 24  # number of temperature points in the sweep
# q-mesh for the 1/T1 BZ sum, sub-sampled from Nx,Ny,Nz (each entry rounded down to a divisor).
# COST IS LINEAR IN Nq: total ~ NT*Nq*Nk, so nmr_qsize=[Nx,Ny,Nz] is O(Nk^2) and only
# tractable on toy meshes.  The k-sum sets the RESOLUTION (needs delta >~ v_F*dk), the q-sum
# only INTEGRATES a smooth function (error ~Nq^-2), so 8-16/axis is normally plenty -- except
# near a magnetic instability, where chi_s(q) sharpens at q_AF.  Verify with nmr_qconv.
nmr_qsize = [8, 8, 1]
nmr_qfold = True  # fold q with -q (exact for centrosymmetric+TRS, which the BdG basis already assumes); halves the cost
nmr_qconv = True  # before the sweep, compare the q-sum on nmr_qsize and on half that mesh and print the relative change
nmr_w0 = None  # NMR probe frequency [eV] for chi''(q,w0)/w0; None -> 0.5*delta (needs w0 <~ delta << gap)
nmr_hf_A, nmr_hf_B = 1.0, 0.0  # hyperfine form factor A(q)=A+2B(cos qx+cos qy) (Mila-Rice); B=0 -> flat weight, A=-4B kills the (pi,pi) response
# --- use a self-consistent RPA/FLEX/Eliashberg gap instead of the gap_sym form factor ---
# Set to the base name (no extension) of the gap exported by LIN_ELIASHBERG/NONLIN_ELIASHBERG
# with sw_out_self=True (output_gap_wannier -> 'gap_wannier').  Delta(R) is mesh independent,
# so the Eliashberg run may use a coarser k-mesh than this sweep needs.  None -> gap_sym/delta0.
# NOTE: gap_sym/delta0 must still be VALID (the shared SC-chi block builds a gap from them
# before this one replaces it); their values are then irrelevant, but set nmr_spsym since
# gap_sym no longer describes the pairing channel.
nmr_gap_file = None
nmr_gap_extrapolate = True  # True: fit Delta(iw_n)->Delta(0) (removes the O((pi T)^2) bias of the lowest slice); False: use slice nmr_gap_iw averaged over nmr_gap_navg
nmr_gap_iw, nmr_gap_navg = 0, 1  # slice / slice-average used when nmr_gap_extrapolate=False
# Overall amplitude Delta(0) [eV] to rescale the loaded gap to, measured as the maximum ON THE
# FERMI SURFACE.  REQUIRED for a linearized-Eliashberg gap (an eigenvector: scale is arbitrary);
# set None to keep the stored amplitude of a nonlinear/self-consistent gap.
nmr_gap_max = None
nmr_spsym = None  # None: take the singlet/triplet channel from the sign of gap_sym; True/False: override (needed with nmr_gap_file, where gap_sym no longer describes the gap)
# mu0=9.85114560061123
# k_sets=[[0., 0., 0.],[.5, 0., 0.],[.5, .5, 0.]]
# xlabel=[r'$\Gamma$','X','M']
# m_diis_num=2
at_point = [0., .5, 0.]  # reduced q-point [0..1] for the single-q susceptibility (option 8)
# --- impurity supercell (option 21): on-site impurity potential in a periodic supercell ---
imp_cell = [8, 8, 1]  # supercell size per axis; must exceed twice the hopping range (size^3*Norb is diagonalized)
imp_sites = None  # 0-based impurity site indices in the supercell (gen_rlist order); None -> random at imp_conc
imp_conc = 0.05  # impurity concentration used when imp_sites is None (at least one impurity)
imp_seed = 0  # random seed of the impurity placement (None: different every run)
imp_V = 1.0  # on-site impurity potential [eV]: a number (times the unit matrix) or an Norb x Norb matrix
orb_dep = False  # use orbital dependence U,J (True: use Umat/Jmat matrices below; False: constant U,J)
UJ_material = None  # with orb_dep=True, name of a stored cRPA parameter set loaded from UJ_dir instead of
# writing Umat/Jmat out below, e.g. 'FeSe','FeTe','LiFeAs','BaFe2As2','LaFeAsO',
# 'LaFeAsOnakamura'. None: use the Umat/Jmat given here
UJ_dir = 'UJ'  # directory of per-material U/J files ({UJ_dir}/{UJ_material}.toml; .json also read), or a
# single multi-material TOML file (one [table] per material); each file records the reference
# of its set (source/doi). The orbital order must match that of the hamiltonian.
# The stored sets hold one site (5 Fe-d orbitals) and are copied onto every site
# of a multi-site model (2-Fe, 10 orbitals) following site_prof
Umat = None  # orbital-dependent U matrix (Norb x Norb, or one site worth of orbitals when
# site_prof has several equal sites: it is then copied onto every site);
# set when orb_dep=True and UJ_material=None
Jmat = None  # orbital-dependent J matrix, same shape as Umat
sw_unit = True  # True: use physical constants (SI/eV units), False: set all constants to 1 (dimensionless test mode)
sw_tdf = False  # True: compute the transport distribution function first, then energy-integrate (energy-dependent tau)
tau_mode = 'const'  # relaxation time model for the Boltzmann transport modes (option 5).
# 'const': tau_const [fs], k- and band-independent. R_H then carries tau^2/tau^2
#         and is tau-INDEPENDENT, i.e. a pure band-structure number with no T
#         dependence beyond the Fermi window.
# 'epa'  : tau(k) from the electron-phonon averaged coupling; needs epa_file.
#         This is what makes R_H(T) and sheet-dependent mobilities meaningful.
# 'e2'   : Fermi-liquid tau, 1/tau ∝ (e-mu)^2 + (pi k_B T)^2 (electron-electron scattering),
#         tau = tau_const on the Fermi level.
# 'dos'  : tau inversely proportional to the DOS at the state's energy (golden-rule scattering
#         rate), used as a ratio: tau(e) = tau_const*D(mu)/D(e), i.e. tau_const on the Fermi level.
# 'flex' : FLEX quasiparticles from self_en_wannier.npz (FLEX run with sw_out_self=True at the
#         SAME temperature). Sigma(k,iw_n) is Pade-continued on the FLEX k-mesh, interpolated
#         multilinearly onto the transport mesh and projected on the bands: Z, Re Sigma(0) and
#         Gamma(e*) = -Im Sigma^R at the quasiparticle energy e* give the energies e*, velocities
#         Z*v and tau* = hbar/(2 Z Gamma); mu is the FLEX one (fill is not used). tau_const is
#         ignored. Bubble only (no current-vertex correction).
epa_file = None  # path to the epa.x (job='egrid') output read when tau_mode='epa'
epa_pair_file = None  # epa.x (job='egrid') output: adds the phonon-mediated attraction built from it to the
# linearized Eliashberg pairing (option 15, no SOC), on top of the spin/charge fluctuations. None: off.
# Isotropic s-wave, normalized to the EPA lambda = 2 N_F sum_nu <g^2>/w_nu; no phonon self-energy.
epa_pair_scale = 1.0  # factor on that phonon coupling (e.g. to scan the phonon / spin-fluctuation balance)
flex_gamma_imp = 0.0  # tau_mode='flex': extra impurity damping Gamma_imp [eV] added to -Im Sigma^R (Matthiessen)
flex_npade = 40  # tau_mode='flex': number of Matsubara points used in the Pade continuation of Sigma
hall_mesh_list = None  # k-mesh convergence sweep for the Hall coefficient (option 5), e.g. [20,30,40]
# or [(20,20,10),(30,30,15)]. sigma^(1) converges much more slowly than sigma
# and can come out with the wrong SIGN on a coarse mesh, so check before use.
# None: single mesh (Nx,Ny,Nz), no convergence information.
sw_omega = False  # True: real freq, False: Matsubara freq.
sw_rescale_flex = True  # True: rescale self energy to make max|Sigma|~U, False: no rescaling
freq_grid = 'matsubara'  # Matsubara-axis grid of FLEX/Eliashberg (option 14,15,16,23): 'matsubara' = Nw uniform points,
# 'ir' = sparse sampling of the IR basis (needs sparse-ir; Nw is then fixed by the basis, ~15-30 points, and the
# truncation error of the uniform grid is gone). self_en_wannier.npz stays uniform-Matsubara
ir_wmax = None  # freq_grid='ir': spectral cutoff [eV]; None = 2 x the largest |e_k - mu| (FLEX spreads spectral weight
# beyond the bare band, so the band width alone is too small for FLEX)
ir_eps = 1.0e-8  # freq_grid='ir': accuracy of the IR basis
sw_chi0_tail = True  # True: tail-corrected chi0 in FLEX/Eliashberg (conv[G]-conv[G0]+analytic reference; O(1/Nw^2) Matsubara truncation error). no-SOC path only
sw_self = False  # True: use calculated self energy for spectrum band plot
sw_out_self = False  # True: write the FLEX self-energy to sigma.bin/self_en.npz (also triggers gap output in option 15,23; option 23 with sw_self also writes the SC-state self-energy to self_en_sc.npz)
sw_in_self = False  # True: load the previous self-energy from sigma.bin as the initial guess for the SC loop (a seed with a different Nw is auto re-gridded in place onto the current Nw)
sigma_in_scale = 1.0  # factor multiplying the sigma seed loaded via sw_in_self; for U-annealing set (U_new/U_old)**2 (Sigma ~ U^2)
sw_from_file = False  # True: read the self-energy from sigma.bin and skip FLEX (solve Eliashberg with it directly)
sw_check_only = False  # True: stop after linear Eliashberg (also stops if Stoner factor>=1 or lambda<1)
# ----- symmetry-free linearized Eliashberg (option 15 with sw_nosym): Zeeman field / FFLO -----
sw_nosym = False  # True: option 15 solves the gap equation on the full BZ with NO TRS/inversion (bare G, RPA);
# a spinless model (sw_soc=False) is spin-doubled automatically. No singlet/triplet projection: gap_sym only seeds
zeeman_h = [0.0, 0.0, 0.0]  # Zeeman field (hx,hy,hz) [eV] added as H_Z = h.sigma (g*mu_B*B/2 absorbed) (sw_nosym only); mu is re-solved at the same fill
fflo_q = [[0, 0, 0]]  # pair momenta Q in integer k-mesh units, Q=(iQx/Nx,iQy/Ny,iQz/Nz); lambda(Q) -> lambda_Q.dat
# e.g. [[q, 0, 0] for q in range(5)] scans along kx. chi0/V_Delta are computed once and reused for every Q
nosym_arnoldi_m = 60  # Krylov dimension of the nosym solver (small lambda may need >~60 for 1e-6 accuracy)
# ----- EILENBERGER (homogeneous quasiclassical) parameters -----
eil_coupling = 0.3  # dimensionless separable pairing coupling lambda (with <phi^2>_FS=1)
eil_wc = 0.5  # fixed Matsubara cutoff energy [eV] (sets the pairing scale / Tc)
eil_imp_gamma = 0.0  # non-magnetic impurity scattering strength Gamma [eV] (0=clean)
eil_imp_c = 1.0e8  # T-matrix cot(delta0): large=Born limit, 0=unitary limit
eil_fs_width = 5.0e-3  # Gaussian Fermi-surface broadening [eV]
eil_method = 'normalization'  # homogeneous (g,f) route: 'normalization' (fast) or 'riccati' (vortex-lattice-ready)
eil_find_tc = False  # True: bisect for Tc at the current impurity setting
eil_imp_sweep = False  # True: sweep Gamma and write Tc(Gamma) to eilenberger_tc.dat
eil_imp_list = None  # array of Gamma values [eV] for the sweep (e.g. np.linspace(0,0.05,11))
eil_pauli = False  # True: Zeeman/Maki Pauli-limiting sweep (singlet gap Delta(h), spinodal, Zeeman-split DOS)
eil_free_energy = False  # True: condensation free energy (Omega_s-Omega_n)/N0 vs T (coupling-independent; writes free_energy.dat)
eil_spin = False  # True: spin-2x2 Zeeman response, singlet vs triplet d-vector (d||h Pauli-limited, d_|_h immune)
eil_lambda = False  # True: superfluid density rho_s(T)/penetration depth lambda(T) sweep (s exp-flat, d linear-in-T)
eil_fs = False  # True: model-FS + Fermi-velocity penetration depth (anisotropic lambda_xx/lambda_yy)
eil_fs_kind = None  # Fermi surface for ALL eilenberger modes (homogeneous penetration / surface / vortex): None=isotropic (cylinder; homogeneous falls back to 'ellipse'), in-plane model FS 'iso'/'ellipse'(params=(mx,my))/'tb'(params=t), 3D model FS 'cyl'(corrugated cylinder, params=(t,tz))/'sphere'/'spheroid'(params=(mx,my,mz)) -- use with eil_fs_nkz>1, or 'wannier' (loaded band FS+v_F; surface/vortex symmetry/multiband from gap_sym,delta0)
eil_fs_params = (1.0, 0.4)  # model-FS parameters (ellipse masses or tb hopping)
eil_fs_nkz = 1  # k_z slices stacked into the FS (eil_fs_kind='wannier' or a 3D model kind): 1=single k_z=0 cut (quasi-2D), >1=true 3D Fermi surface. REQUIRED (>=16) for k_z-dependent gaps: a horizontal line node lives on the k_z=0 plane, where a single-slice FS makes the gap identically zero
eil_fs_traj = None  # trajectory reduction for the VORTEX/lattice solvers, whose cost is linear in the number of FS points (a 3D FS has 10^3-10^4 vs the 24 directions of the model cylinder): None=use every point, int=number of beta bins, (nbeta,nv,nphi)=full control (direction / |v_par| quantile / phi bins). default (48,4,8) reduces a 1920-direction FS to 176 with an x17.7 speedup, the bulk gap to 0.5% and the vortex-core LDOS peak to 0.2%. Raise until the answer stops moving. Not used by the surface solver (it needs k_par and the band index for the specular partner)
eil_pair_gauge = 'trs'  # pair-partner gauge of the band gap projection: 'trs' (default; spinless/real hoppings, phi=u^dag.Delta.u, gauge invariant), 'soc' (spinful basis, time-reversed partner (i sigma_y)u*), 'diag' (legacy independently diagonalized u(-k); GAUGE DEPENDENT, comparison only). The condition is time-reversal symmetry of the NORMAL state H(k), not of Delta: a CHIRAL order parameter is fine (its winding is passed through to phi), a time-reversal-broken NORMAL state (ferromagnetic SC, Zeeman inside H, Haldane-type hoppings) is not -- check_pair_partner warns when |<u(-k)|T u(k)>| < 1
eil_spin_order = 'block'  # basis ordering for eil_pair_gauge='soc': 'block'=[orb1..orbN up, orb1..orbN dn], 'interleave'=[orb1 up, orb1 dn, ...]
eil_zeeman = 0.0  # Zeeman (Maki) field [eV] for the LDOS (surface: splits the d[110] ZEBS into +-h; vortex: spin-splits the core bound states; periodic lattice (eil_vort_lattice_sc): self-consistent Pauli depairing)
# ----- EILENBERGER inhomogeneous (surface & vortex, Riccati; model cylindrical FS) -----
eil_ldos = True  # True: also compute the real-frequency LDOS (bound/core states) (surface & vortex)
eil_surf_beta = 0.785398  # surface orientation [rad]: 0=[100], pi/4(0.7854)=[110] (d-wave ZEBS)
eil_surf_dvector = False  # True: self-consistent triplet d-vector TEXTURE at the surface (dominant p_x(e_x) + subdominant p_y(e_z), spin-matrix Riccati)
eil_dvec_subratio = 0.9  # subdominant/dominant coupling ratio for the d-vector texture (~0.85 is the bulk threshold)
eil_vort_lxi = 8.0  # vortex cell half-width in coherence lengths xi (isolated vortex, field=0)
eil_vort_ngrid = 81  # vortex 2D grid points per axis
eil_vort_field = False  # True: also compute the self-consistent finite-kappa Maxwell field profile B(rho) of the vortex (uses eil_kappa)
eil_vort_chiral = False  # True: self-consistent isolated CHIRAL vortex (p+ip for eil_chiral_ell=1, d+id for 2) via the multi-component complex-amplitude spin-matrix solver -- the case the scalar vortex solver refuses, since the core induces the OPPOSITE chirality and the amplitude is genuinely complex. Both chiralities share one coupling (degenerate partners of the same irrep); the bulk is the single-channel isotropic gap, not the equal-amplitude (nematic) combination
eil_chiral_ell = 1  # chirality of eil_vort_chiral: 1=p+ip, 2=d+id
eil_chiral_m = 1  # winding of the DOMINANT chiral component (+1 parallel to the chirality, -1 antiparallel); the induced opposite component follows m_- = m_+ + 2*ell (chiral_windings)
eil_field_dir = None  # field / vortex-line direction (isolated vortex AND the finite-field lattice), e.g. (0,0,1)=c axis (default), (1,0,0)=in-plane. The vortex lines run along B, so the order parameter varies in the plane PERPENDICULAR to it and the problem stays 2D; that plane's two axes are rescaled by their rms Fermi velocities, so a square grid still fits the elliptical core (xi_1/xi_2 is reported). Needs a 3D Fermi surface (eil_fs_nkz>1 or a 3D eil_fs_kind) for anything but B||c. A finite eil_field with an anisotropic plane is NOT supported yet (the circular-cell Doppler assumes the unscaled plane) and raises
eil_chiral_dvec = 'z'  # equal-spin direction of the chiral triplet: 'z' or 'x' ('y' is proportional to the identity and is not supported by the unitary matrix-Riccati seed)
eil_vort_dvector = False  # True: self-consistent triplet d-vector TEXTURE around the vortex core (dominant p_x(e_x) winding + core-localized subdominant p_y(e_z), 2D spin-matrix Riccati; uses eil_dvec_subratio)
eil_vort_current = False  # True: circulating charge supercurrent j_phi(rho) of an isolated vortex (writes vortex_current.dat)
eil_vort_maxwell = False  # True: circular-cell vortex with the self-consistent finite-kappa vector potential A(r) (Maxwell back-reaction; needs eil_field>0, uses eil_kappa)
eil_vort_lattice_sc = False  # True: je-style self-consistent TRUE periodic lattice (formulation A, extreme type-II): complex Psi(r) with a real node at every core + full Abrikosov supercurrent; sweeps eil_field_list -> <N(0)>(B) (d~sqrt(B) Volovik)
eil_vort_scA = False  # True (lattice_sc, finite eil_kappa): fully self-consistent vector potential A from the quasiclassical current j_s=<v_F Im g> (je A_renew), instead of the analytic London A
eil_gap_orbital = None  # None, or an Norb x Norb orbital-basis pair-potential matrix (or callable kfrac->NxN): the band gap is its low-energy PROJECTION onto the FS bands (Nagai-Nakamura, JPSJ 85 074707 (2016) Eq.43; needs wannier FS); supersedes gap_sym/delta0
eil_gap_file = None  # None, or the base name (no extension) of an RPA/FLEX gap exported by LIN_ELIASHBERG/NONLIN_ELIASHBERG with sw_out_self=True (output_gap_wannier -> 'gap_wannier'): loads Delta(R,iw) and uses it as eil_gap_orbital (projected to the FS bands); supersedes eil_gap_orbital/gap_sym. Use a self-consistent RPA gap (e.g. KFe2As2) as the vortex pairing form factor
eil_gap_iw = 0  # starting Matsubara index for eil_gap_file (0=lowest iw_0; sharpest gap symmetry/anisotropy, matches the FS gap usually quoted)
eil_gap_navg = 1  # number of consecutive Matsubara slices to average for eil_gap_file (1=single iw_gap_iw slice; >1 smooths noise at the cost of slightly diluting the anisotropy)
eil_vort_tilt = 0.0  # field tilt theta [deg] from the c-axis (quasi-2D): orbital uses B_z=B cos(theta), Zeeman -> h/cos(theta) (Pauli/orbital ratio)
eil_nvortex = 1  # vortices (flux quanta) per computational cell of the periodic lattice (supercell; n^2 reduces to the primitive cell)
eil_field = 0.0  # vortex lattice field B/Hc2 (0=isolated vortex; >0=circular-cell lattice w/ Doppler)
eil_field_list = None  # list of B/Hc2 to sweep <N(0)>(B) on the TRUE periodic lattice (e.g. [0.04,0.08,0.16,0.32]); None=single field
eil_kappa = 100.0  # GL kappa=lambda/xi for the periodic lattice (large=extreme type-II; finite=London screening/Maxwell)
eil_lattice = 'square'  # periodic vortex lattice geometry: 'square' or 'triangular'
# ------------------------ initial parameters are above -------------------------------
# ----------------------------------main functions-------------------------------------
# -------------------------------- import packages ------------------------------------
import os
import numpy as np
import libs.flibs as flibs, libs.plibs as plibs
import libs.plots as plots   # figures (matplotlib); plibs holds the calculations only
import matplotlib.pyplot as plt
# ----------------------------- initial settings --------------------------------------
sw_gen_sym = False
try:
    kz
except NameError:
    kz = 0.0
try:
    mu0
    sw_calc_mu = False
except NameError:
    sw_calc_mu = True
try:
    kscale
except NameError:
    kscale = 1.0
try:
    gap_rcut
except NameError:
    gap_rcut = None
try:
    abc
except NameError:
    abc = [1., 1., 1.]
    print('abc not found, set to 1,1,1')
try:
    alpha_beta_gamma
except NameError:
    alpha_beta_gamma = [90., 90., 90]
    print('alpha_beta_gamma not found, set to 90,90,90')
alatt = np.array(abc)
deg = np.array(alpha_beta_gamma)
if option in MODES_SYMLINE:  # modes using k_sets/xlabel (symmetry line)
    try:
        k_sets
        xlabel
    except NameError:
        sw_gen_sym = True
        k_sets, xlabel2 = plibs.get_symm_line(brav)
        try:
            xlabel
            if len(xlabel) != len(k_sets):
                xlabel = xlabel2
        except NameError:
            xlabel = xlabel2
# physical constants in SI/eV units, or all 1 in the dimensionless test mode (sw_unit=False)
units = plibs.unit_constants(sw_unit)
kb, ihbar = units.kb, units.ihbar   # used directly below (the drivers take `units`)
try:
    tempK
    try:
        temp
        print('temp & tempK find, use temp')
    except NameError:
        temp = tempK*kb
except NameError:
    print('tempK not found')
    exit()
try:
    m_diis_num
except NameError:
    m_diis_num = 5
if option in MODES_NEED_ROTMAT:
    try:
        RotMat
        RotMat = np.array(RotMat)
    except NameError:
        print('No RotMat')
        RotMat = np.eye(3)
# ------------------------ define functions -------------------------------------------


def get_spin_lists(Norb: int) -> tuple[np.ndarray, np.ndarray]:
    """Spin label slist [Norb] (+1 up, -1 down) and 1-based opposite-spin partner invs [Norb]
    of a spinful (SOC) basis: the header's slist/invs if given, otherwise the default
    ordering (all up, then all down; the up block takes the extra orbital when Norb is odd)."""
    slist = globals().get('slist')
    if slist is None:
        slist = np.ones(Norb, dtype=np.int64)
        slist[(Norb+1)//2:] = -1
    invs = globals().get('invs')
    if invs is None:
        invs = np.concatenate([np.arange((Norb+1)//2, Norb), np.arange((Norb+1)//2)])+1
    return np.asarray(slist, dtype=np.int64), np.asarray(invs, dtype=np.int64)


def set_eil_gauge(Norb: int) -> dict:
    """Pair-partner gauge kwargs for the band-gap projection (eil_pair_gauge).
    'soc' additionally needs the spin-partner map of the spinful basis."""
    if eil_pair_gauge == 'soc':
        return dict(gauge='soc', spin_map=plibs.spin_pair_map(Norb, eil_spin_order))
    return dict(gauge=eil_pair_gauge)


def get_eil_gap_fs(rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray, avec: np.ndarray, Arot: np.ndarray,
                   temp: float, fill: float) -> tuple:
    """Shared FS/gap preparation for the Eilenberger modes: the RPA/FLEX gap
    (eil_gap_file) or the explicit orbital gap (eil_gap_orbital) as form factor,
    plus the real Wannier-band FS when eil_fs_kind='wannier'.
    Returns (fs, fs_kind): fs_kind=None when the gap is baked into fs['phi'],
    otherwise fs=None and fs_kind=eil_fs_kind (model FS/cylinder)."""
    gorb = (plibs.gap_orbital_from_wannier(eil_gap_file, eil_gap_iw, eil_gap_navg)  # RPA/FLEX gap as form factor
            if eil_gap_file else eil_gap_orbital)
    if eil_fs_kind == 'wannier':  # real Wannier-band FS + v_F (gap symmetry/multiband from gap_sym,delta0)
        fs = plibs.build_wannier_fs(rvec, ham_r, S_r, avec,
                                    plibs.get_mu(ham_r, S_r, rvec, Arot, temp, fill),
                                    nkz=eil_fs_nkz,  # >1: stack k_z slices into a true 3D FS
                                    gap_sym=gap_sym, delta0=delta0, gap_orbital=gorb,  # gap_orbital=projection (Nagai)
                                    **set_eil_gauge(len(ham_r[0])))
        return fs, None  # the (int) gap_sym is baked into fs['phi']
    return None, eil_fs_kind  # model FS/cylinder: the int gap_sym -> continuum harmonic


def main() -> None:
    omp_num, omp_check = flibs.omp_params()

    # ===== Input file and parameter validation =====
    # Validate ftype (any integer outside {0,1,2,3,4} uses the Hopping.dat branch)
    if not isinstance(ftype, int):
        print(f"Error: ftype={ftype!r} must be an integer", flush=True)
        return
    # Check the files each format actually reads: fname is a bare prefix for ftype 2
    # (_hr.dat is appended) and a directory for ftype 0, 4 and the Hopping.dat branch,
    # so testing fname itself would wrongly reject those.
    required = {0: ['irvec.txt', 'ndegen.txt', 'ham_r.txt'],
                1: [''],
                2: ['_hr.dat'],
                3: [''],
                4: ['HamRsMLO', 'HamiltonianPMTInfo']}.get(ftype, ['Hopping.dat'])
    sep = '' if ftype in (1, 2, 3) else os.sep       # ftype 2 appends a suffix, dirs join with /
    for f in required:
        if not os.path.exists(fname+sep+f):
            print(f"Error: Hamiltonian file '{fname+sep+f}' not found", flush=True)
            return
    # Validate brav (any integer outside {0..7} uses the monoclinic branch)
    if not isinstance(brav, int):
        print(f"Error: brav={brav!r} must be an integer", flush=True)
        return

    # ===== Computational parameter validation =====
    # Mesh size check
    if Nx <= 0 or Ny <= 0 or Nz <= 0 or Nw <= 0:
        print(f"Error: Invalid mesh size (Nx={Nx}, Ny={Ny}, Nz={Nz}, Nw={Nw})", flush=True)
        return
    # Temperature check
    # Many kernels use 1/temp or finite-T Fermi functions; T=0 is not supported in this implementation.
    if tempK <= 0:
        print(f"Error: Temperature (tempK) is non-positive ({tempK} K)", flush=True)
        return
    # U, J check
    if U < 0 or J < 0:
        print(f"Error: U, J are negative (U={U}, J={J})", flush=True)
        return
    # delta check
    if delta <= 0:
        print(f"Error: Broadening (delta) is non-positive ({delta})", flush=True)
        return
    # fill check
    if fill < 0:
        print(f"Error: Filling (fill) is negative ({fill})", flush=True)
        return
    # Emin, Emax check
    if Emin >= Emax:
        print(f"Error: Emin >= Emax ({Emin} >= {Emax})", flush=True)
        return
    # tau_const check (only the tau_modes scaled by it; 'epa' and 'flex' ignore it)
    if tau_mode in ('const', 'e2', 'dos') and (tau_const is None or tau_const <= 0):
        print(f"Error: Relaxation time (tau_const) must be positive for tau_mode='{tau_mode}' ({tau_const})",
              flush=True)
        return
    # =========================================

    # import hamiltonian
    if ftype == 3:
        rvec, ham_r, S_r, Norb, Nr = plibs.import_MLO_hoppings(fname)
    elif ftype == 4:
        # fname is the directory holding the ecalj job_mlo binaries
        rvec, ham_r, S_r, Norb, Nr = plibs.read_HamRsMLO(fname)
    else:
        rvec, ham_r, Norb, Nr = plibs.import_hoppings(fname, ftype)
        S_r = []
    plist = flibs.get_plist(rvec, ham_r)
    print("Effective parity of wannier functions:", plist, flush=True)
    # The Eliashberg solvers rebuild Delta(-k) from Delta(k) with this operator, so it has to
    # BE a symmetry of H(R). The diagonal parity above is exact for single-site cells but not
    # when inversion exchanges sites (2-Fe cells etc.), where the true operator is a signed
    # permutation. Validate, and switch to the permutation form when needed.
    iperm = None
    # the symmetry-free solver (sw_nosym) works on the full BZ and never rebuilds Delta(-k)
    if option in MODES_NEED_PARITY and not (option == CalcMode.LIN_ELIASHBERG and sw_nosym):
        psgn, pperm, pstat, pres = plibs.inversion_op(rvec, ham_r, plist, tol=inv_tol)
        res_str = 'n/a' if pres is None else f'{pres:.2e}'
        print(f"Inversion operator: {pstat} (diag(plist) residual {res_str})", flush=True)
        if pstat == 'monomial' and sw_soc:
            # the SOC Eliashberg kernels and output_Fk take a DIAGONAL parity only (no iperm)
            print("Error: inversion exchanges orbitals (signed permutation), which the SOC "
                  "Eliashberg path does not support (diagonal parity only).", flush=True)
            return
        if pstat == 'monomial':
            plist, iperm = psgn, pperm
            print(f"  -> using signed permutation: perm={pperm.tolist()} sgn={psgn.astype(int).tolist()}",
                  flush=True)
        elif pstat != 'diagonal':
            print(f"Error: inversion cannot be represented by a signed permutation ({pstat}).", flush=True)
            print("       Delta(-k) would be reconstructed with a non-symmetry, giving the leading", flush=True)
            print("       solution of a wrong eigenproblem. Raise inv_tol if the model is a fit, or", flush=True)
            print("       rebuild the Wannier functions in a symmetry-respecting gauge.", flush=True)
            return

    # ===== Parameter checks dependent on number of orbitals =====
    print(f"Number of orbital = {Norb}", flush=True)

    # Validate olist
    for i, ol in enumerate(olist):
        if isinstance(ol, int):
            if ol < 0 or ol >= Norb:
                print(f"Error: olist[{i}]={ol} is invalid. Valid range: 0-{Norb-1}", flush=True)
                return
        elif isinstance(ol, list):
            for j, olj in enumerate(ol):
                if olj < 0 or olj >= Norb:
                    print(f"Error: olist[{i}][{j}]={olj} is invalid. Valid range: 0-{Norb-1}", flush=True)
                    return

    # Check fill does not exceed number of bands
    if fill > Norb:
        print(f"Error: filling={fill} exceeds number of bands={Norb}. Valid range: 0 < fill <= {Norb}", flush=True)
        return

    if orb_dep:
        # Orbital-dependent interaction mode needs the U/J matrices, either written out in
        # the header (Umat/Jmat) or loaded from a stored set (UJ_material).
        if UJ_material is None:
            if globals().get('Umat') is None or globals().get('Jmat') is None:
                print("Error: orb_dep=True requires Umat and Jmat to be defined (or UJ_material to be set)", flush=True)
                return
            Umat_arr, Jmat_arr = Umat, Jmat
        else:
            try:
                Umat_arr, Jmat_arr = plibs.read_UJ(UJ_material, UJ_dir)
            except (OSError, KeyError, ValueError) as err:
                print(f"Error: {err}", flush=True)
                return
        # A stored set covers ONE site (e.g. the 5 Fe-d orbitals), so it is replicated on
        # every site of a multi-site hamiltonian (a 2-Fe, 10-orbital model); the checks on
        # shape/dtype live in expand_UJ_sites, which also returns contiguous float64 arrays
        # for stable ctypes calls to the Fortran wrappers.
        try:
            Umat_arr, Jmat_arr = plibs.expand_UJ_sites(Umat_arr, Jmat_arr, Norb,
                                                       globals().get('site_prof', [1]))
        except ValueError as err:
            print(f"Error: {err}", flush=True)
            return
        globals()['Umat'], globals()['Jmat'] = Umat_arr, Jmat_arr
    # =============================================

    # set lattice vector
    avec, Arot = plibs.get_ptv(alatt, deg, brav)
    # rotation axis
    if sw_dec_axis:
        # Fermi-surface drawing only.  Re-expressing R in the conventional (orthogonal)
        # axes is what lets a body-centred 122 Fermi surface be drawn in the familiar box
        # and compared between the 10-orbital (2-Fe) and 5-orbital (1-Fe) models; E(k) at
        # a given PHYSICAL k is preserved exactly (see the checks in the manual).  The
        # half-integer R components this leaves for brav=1,2,6,7 are intrinsic to that
        # unfolding -- and they are also why H(k) stops being periodic in the decomposed
        # fractional k, so BZ sums, FFT convolutions, the irreducible wedge and det(avec)
        # as a cell volume are all invalid here.  Hexagonal/trigonal cells cannot be put
        # on an integer orthogonal lattice at all (irrational components).  Hence the
        # whitelist: the transformation itself is left exactly as it was.
        if option not in MODES_DEC_AXIS:
            print(f"Error: sw_dec_axis=True is for the band and Fermi-surface plots only "
                  f"(option={int(CalcMode.BAND)},{int(CalcMode.FERMI_2D)},"
                  f"{int(CalcMode.FERMI_3D)}), not for "
                  f"option={int(option)} ({option.description}). The decomposed axes leave "
                  "H(k) non-periodic in the fractional k-mesh, which invalidates every BZ "
                  "sum, FFT convolution, irreducible wedge and cell volume. Set "
                  "sw_dec_axis=False; gap symmetries no longer need it either, since "
                  "plibs.lattice_harmonic builds them from Cartesian R shells.", flush=True)
            return
        print("sw_dec_axis: drawing in the decomposed (conventional) axes", flush=True)
        if option == CalcMode.BAND:
            # The auto-generated path came from the centred-cell brav and means a different
            # physical k in the decomposed axes.  Those axes ARE a simple orthogonal cell,
            # so the simple-lattice path is the meaningful one here -- and it is the path
            # the 1-Fe model is plotted along, which is the point of the comparison.
            if globals().get("sw_gen_sym", False):
                ks, xl = plibs.get_symm_line(0)
                globals()["k_sets"], globals()["xlabel"] = ks, xl
                print("sw_dec_axis: symmetry path regenerated for the decomposed "
                      "(simple orthogonal) axes", flush=True)
            else:
                print("sw_dec_axis: k_sets is used as given -- it must be in the "
                      "DECOMPOSED axes, not those of the centred cell", flush=True)
        rvec1 = Arot.T.dot(rvec.T).T
        rvec = rvec1.copy()
        if brav in {1, 2}:
            rvec[:, 2] *= 2.
            avec = alatt*np.eye(3)
            avec[:, 2] *= .5
        elif brav == 5:
            rvec *= 2.
            avec = (alatt*np.eye(3))*.5
        else:
            avec = alatt*np.eye(3)
    bvec = plibs.get_bvec(avec)  # set recp. lattice
    # Build the gap form factors from CARTESIAN harmonics on the R shells of this
    # lattice (plibs.lattice_harmonic).  Identical to the old fractional formulas on
    # an orthogonal single-site lattice, and correctly labelled on fcc/bcc/hexagonal
    # ones, where cos(2 pi k_x) etc. do not transform as the harmonic they name.
    plibs.set_harmonic_lattice(avec, rvec, s_rcut=gap_rcut)
    if omp_check:  # OMP properties
        print("OpenMP mode", flush=True)
        print(f"Number of OpenMP threads = {omp_num}", flush=True)
    print(f"calc mode {int(option)}: "+option.description, flush=True)

    # ===== Additional input validation (lower priority) =====
    # J > U warning
    if J > U:
        print(f"Warning: J={J} > U={U} is physically unusual. Please verify", flush=True)

    # Validate gap_sym during eliashberg/flex calculations
    if option in MODES_FLEX_ELIASH:  # eliashberg/flex calculations
        if gap_sym not in {-1, 0, 1, 2, 3}:
            print(f"Warning: gap_sym={gap_sym} is non-standard. Common values: -1,0,1,2,3", flush=True)
    if freq_grid not in ('matsubara', 'ir'):
        print(f"Error: freq_grid={freq_grid!r} must be 'matsubara' or 'ir'", flush=True)
        return
    if option == CalcMode.NONLIN_ELIASHBERG and sw_soc:
        print("Error: the nonlinear Eliashberg solver has no SOC version; set sw_soc=False", flush=True)
        return
    # the SC-state chi modes build their gap from delta0: a number (one shape) or one amplitude per band
    if option in MODES_CHIS_SC:
        if not (np.isscalar(delta0) or len(delta0) == Norb):
            print(f"Error: delta0 must be a number or a list of Norb={Norb} band amplitudes "
                  f"(got {len(delta0)})", flush=True)
            return

    # Validate at_point for chis at q-point calculation (normal and SC state)
    if option in MODES_CHIS_QPOINT:  # chis at q-point
        if len(at_point) != 3:
            print(f"Error: at_point has {len(at_point)} elements. Required: 3 elements", flush=True)
            return
        if any(p < 0 or p > 1 for p in at_point):
            print(f"Warning: at_point={at_point} values outside [0,1] range. May not be normalized by reciprocal lattice", flush=True)
    # ==========================================
    if option in MODES_COLOR:
        print("color mode: "+color_option.description, flush=True)
    print("Hamiltonian name is "+fname, flush=True)
    if option in MODES_PRINT_UJ:
        if orb_dep:  # orbital-dependent matrices: report their range instead of single numbers
            # Only the pairs sitting on one site enter the vertices; the inter-site blocks
            # are zeros that would otherwise be reported as the minimum of U' and J.
            offdiag = plibs.site_block_mask(Norb, globals().get('site_prof', [1])) & ~np.eye(Norb, dtype=bool)
            msg = (f'orbital dependent U,J'+(f' ({UJ_material} in {UJ_dir})' if UJ_material else '') +
                   f': U_ii= {np.diag(Umat).min():5.2f}-{np.diag(Umat).max():5.2f}')
            if offdiag.any():  # a one-orbital site has no inter-orbital U',J
                msg += (f", U'_ij= {Umat[offdiag].min():5.2f}-{Umat[offdiag].max():5.2f}, "
                        f'J_ij= {Jmat[offdiag].min():5.3f}-{Jmat[offdiag].max():5.3f}')
            print(msg)
        else:  # write constant U,J
            print(f'U= {U:5.2f} and J= {J:5.3f}')
    if option in MODES_NEED_CHI:
        # chiolist is the orbital-pair basis of chi, site the site label of each pair.
        # Header settings are read with globals().get: assigning a name anywhere in main()
        # makes it local for the whole function, so a bare `try: site_prof` would always
        # raise UnboundLocalError and silently discard the user's value.
        chiolist, site = globals().get('chiolist'), globals().get('site')
        if chiolist is None:
            # Build the orbital-pair basis once here so every response/Eliashberg branch
            # shares the same indexing convention when passing chi objects to Fortran.
            try:
                chiolist, site = plibs.get_chi_orb_list(len(ham_r[0]), globals().get('site_prof', [1]))
            except ValueError as err:
                print(f"Error: {err}", flush=True)
                return
        elif site is None:
            print("Error: a user-defined chiolist needs the matching site list (site)", flush=True)
            return
        else:
            chiolist = np.asfortranarray(chiolist, dtype=np.int64)   # Fortran reads ol(Nchi,2)
            site = np.asarray(site, dtype=np.int64)
        if option in MODES_COULOMB_VERTEX:
            print("generate coulomb vertex matrix S")
            # Susceptibility branches only need the static interaction vertex; FLEX/Eliashberg
            # rebuild their own matrices inside the dedicated solvers.
            if orb_dep:
                Smat, Cmat = flibs.gen_SCmatrix_orb(chiolist, site, Umat, Jmat)
            else:
                Smat, Cmat = flibs.gen_SCmatrix(chiolist, site, U, J)
    if option in MODES_KMESH_SYMLINE:
        if sw_gen_sym:
            print('generate symmetry line', flush=True)
        print(f'kmesh = {kmesh}', flush=True)
    elif option in MODES_KMESH_SINGLE:
        print(f'Number of k-mesh = {Nx}', flush=True)
    else:
        print(f'k-mesh is {Nx} {Ny} {Nz}', flush=True)
    if option in MODES_MATSUBARA and freq_grid == 'matsubara':
        print(f'Number of Matsubara freq. = {Nw}', flush=True)
    print("Lattice Vector (Angstrom)", flush=True)
    for i, a in enumerate(avec):
        print(f"a{i+1}: {a[0]:7.4f} {a[1]:7.4f} {a[2]:7.4f}", flush=True)
    print("Reciprocal Lattice Vector (Angstrom^-1)", flush=True)
    for i, b in enumerate(bvec):
        print(f"b{i+1}: %7.4f %7.4f %7.4f"%tuple(b), flush=True)
    if option in MODES_SELF_MU:  # conductivity and impurity/CPA functions calc or set mu themself
        pass
    else:
        if sw_calc_mu:
            mu = plibs.get_mu(ham_r, S_r, rvec, Arot, temp, fill)
        else:
            mu = mu0
            print('use fixed mu')
        print(f'Temperature = {temp:10.3e} eV ({temp/kb:.2f} K)', flush=True)
        print(f'chem. pot. = {mu:.4f} eV', flush=True)
    # Matsubara grid of the FLEX/Eliashberg kernels (Nw_mats, sw_tail_run: Nw and sw_chi0_tail as this run uses
    # them -- the globals themselves are not reassigned inside main)
    Nw_mats, sw_tail_run = Nw, sw_chi0_tail
    if freq_grid == 'ir' and option in MODES_IR:
        try:
            import sparse_ir  # noqa: F401
        except ImportError:
            print("Error: freq_grid='ir' needs the sparse-ir package (pip install sparse-ir)", flush=True)
            return
        wmax = ir_wmax
        if wmax is None:
            _, kl = plibs.gen_klist(min(Nx, 16), min(Ny, 16), min(Nz, 16), sw_pp=False)
            wmax = 2.0*np.abs(plibs.get_eigs(kl, ham_r, S_r, rvec)[0]-mu).max()
        grid = plibs.ir_grid(temp, wmax, eps=ir_eps)
        plibs.set_freq_grid(grid)
        Nw_mats = grid['Nw']
        print(f"frequency grid: IR basis, wmax = {wmax:.3f} eV (Lambda = {wmax/temp:.1f}), eps = {ir_eps:.1e}: "
              f"{Nw_mats} fermionic/bosonic points, {grid['Ntau']} tau points", flush=True)
        if sw_chi0_tail:
            print("  sw_chi0_tail is for the uniform grid (no truncation error on the IR grid): off", flush=True)
            sw_tail_run = False
        if sw_in_self or sw_out_self:
            print("  sigma.bin (sw_in_self / the sigma.bin part of sw_out_self) is uniform-grid only; "
                  "self_en.npz / self_en_wannier.npz are written as usual", flush=True)
    if option == CalcMode.BAND:  # plot band
        klist, spa_length, xticks = plibs.mk_klist(k_sets, kmesh, bvec)
        eig, uni = plibs.get_eigs(klist, ham_r, S_r, rvec)
        # plot_band loops over bands and reads orbital weights per k: uni(k,band,orb) -> (band,orb,k)
        plots.plot_band(eig.T-mu, spa_length, xlabel, xticks, uni.transpose(1, 2, 0), olist, color_option != 0)
    elif option == CalcMode.DOS:  # plot dos
        Nk, klist, eig, uni, kweight = plibs.get_emesh(Nx, Ny, Nz, ham_r, S_r, rvec, avec, sw_uni=True)
        wlist = np.linspace(Emin, Emax, Nw, True)
        Dos = flibs.gen_dos(eig, uni, mu, wlist, delta)
        plots.plot_dos(wlist, Dos, olist, color_option)
    elif option == CalcMode.FERMI_2D:  # 2D Fermi surface plot
        eig2d = plibs.get_eigs_2d(Nx, rvec, ham_r, S_r, RotMat, kz)
        klist, blist = plibs.get_kf_points(eig2d, Nx, mu, kz)
        clist = plibs.get_colors(klist, blist, ihbar*avec.T, rvec, ham_r, S_r, olist, color_option, True)
        plots.plot_FS(clist, klist, color_option)
    elif option == CalcMode.FERMI_3D:  # 3D Fermi surface plot (color_option=GAP: the Eilenberger pairing gap Re[phi(k)])
        polys, centers, blist = plibs.gen_3d_surf_points(Nx, rvec, ham_r, S_r, mu, kscale, bvec)
        # the gap is only needed for the GAP coloring: do not touch eil_gap_file otherwise
        if color_option == ColorMode.GAP:
            gorb = (plibs.gap_orbital_from_wannier(eil_gap_file, eil_gap_iw, eil_gap_navg)  # RPA/FLEX gap as form factor
                    if eil_gap_file else eil_gap_orbital)
            gauge_kw = set_eil_gauge(len(ham_r[0]))
        else:
            gorb, gauge_kw = None, None
        fspolys, fscenters, fscolors = plots.set_init_3dfsplot(color_option, polys, centers, blist, avec, rvec, ham_r,
                                                               S_r, olist, gap_sym, delta0, gorb,
                                                               gauge_kw, ihbar=ihbar)
        plots.plot_3d_surf(fspolys, fscenters, fscolors, color_option, bvec)
    elif option == CalcMode.SPECTRUM:  # plot spectrum
        spec = plibs.spectrum_on_path(k_sets, kmesh, bvec, mu, ham_r, S_r, rvec, Emin, Emax, delta, Nw, temp, sw_self)
        if spec is not None:
            spa_length, xticks, wlist, Akw = spec
            plots.plot_path_map(spa_length, xticks, xlabel, wlist, Akw)
    elif option == CalcMode.CONDUCTIVITY_BT:  # calc conductivity
        # relaxation-time model; dos_delta is the DOS broadening of tau_mode='dos'
        tau_kw = dict(tau_mode=tau_mode, epa_file=epa_file, dos_delta=delta,
                      flex_kw=dict(gamma_imp=flex_gamma_imp, n_pade=flex_npade))
        plibs.get_hall_coe(rvec, ham_r, S_r, avec, Nx, Ny, Nz, fill, temp, tau_const, units, tau_kw,
                           with_spin=sw_soc, mesh_list=hall_mesh_list)
        res = plibs.calc_conductivity_Boltzmann(rvec, ham_r, S_r, avec, Nx, Ny, Nz, fill, temp, tau_const, units,
                                                tau_kw, with_spin=sw_soc, sw_tdf=sw_tdf)
        if res is not None and sw_tdf:
            plots.plot_tdf(res['tdf_e'], res['tdf'])
    elif option == CalcMode.CONDUCTIVITY_PT:  # calc_optical conductivity
        res = plibs.calc_conductivity_lrt(rvec, ham_r, S_r, avec, Nx, Ny, Nz, fill, temp, Nw, delta, units, Emax,
                                          color_thick=color_thick, color_scatt=color_scatt, with_spin=sw_soc)
        if res is not None:
            plots.plot_optical(res)
    elif option in MODES_CHI_NORMAL:  # calc_chis_spectrum
        print("calculate electron energy", flush=True)
        Nk, klist, eig, uni, kweight = plibs.get_emesh(Nx, Ny, Nz, ham_r, S_r, rvec, avec, sw_uni=True)
        if option in {CalcMode.CHIS_SPECTRUM, CalcMode.PHI_SPECTRUM}:  # chis/phi spectrum with symmetry line
            print("generate qlist", flush=True)
            qlist, spa_length, xticks = plibs.mk_qlist(k_sets, Nx, Ny, Nz, bvec)
            if option == CalcMode.CHIS_SPECTRUM:
                w, sp, sus = plibs.calc_path_spectrum('chis', mu, temp, klist, qlist, chiolist, eig, uni, spa_length, Nw, Emax, delta, Smat)
                print("write chis spectrum in png file", flush=True)
                susfname = 'chis_spec.png'
            elif option == CalcMode.PHI_SPECTRUM:
                w, sp, sus = plibs.calc_path_spectrum('phi', mu, temp, klist, qlist, chiolist, eig, uni, spa_length, Nw, Emax, delta)
                print("write phi spectrum in png file", flush=True)
                susfname = 'phi_spec.png'
            plt.contourf(sp, w, abs(sus.imag), 100)
            plt.colorbar()
            plt.hot()
            plt.savefig(fname=susfname, dpi=300)
        elif option == CalcMode.CHIS_QPOINT:  # chis at q-point
            q_point = np.array(at_point)
            chis, chis_orb, wlist = plibs.chis_q_point(q_point, eig, uni, Emax, Nw, mu, temp, Smat, klist, chiolist, delta)
            plt.plot(wlist, chis.imag)
            plt.show()
            for cso in chis_orb.T:
                plt.plot(wlist, cso.imag)
            plt.show()
        else:  # chis/phi qmap at ecut plane
            if option == CalcMode.CHIS_QMAP:  # chis spectrum ecut plane
                sus, chi0, qx, qy = plibs.chis_qmap(Nx, Ny, Ecut, mu, temp, Smat, klist, chiolist, eig, uni, idelta=1.e-3)
                plt.contourf(qx, qy, abs(chi0.imag), 100)
                plt.colorbar()
                plt.jet()
                plt.savefig(fname='chi0map.png', dpi=300)
                susfname = 'chismap.png'
            elif option == CalcMode.PHI_QMAP:
                sus, qx, qy = plibs.phi_qmap(Nx, Ny, Ecut, mu, temp, klist, chiolist, eig, uni, idelta=1.e-3, sw_omega=sw_omega)
                susfname = 'phimap.png'
            if sw_omega:
                plt.contourf(qx, qy, abs(sus.imag), 100)
            else:
                plt.contourf(qx, qy, abs(sus.real), 100)
            plt.colorbar()
            plt.jet()
            # plt.show()
            plt.savefig(fname=susfname, dpi=300)
    elif option in MODES_CHIS_SC:  # calc chis at superconducting state
        sw_spsym = gap_sym < 0
        # model gap Delta(k) (gap_sym form factor, delta0 amplitude) on the full [0,1) k-mesh
        klist, deltak = plibs.sc_gap_on_mesh(Nx, Ny, Nz, ham_r, S_r, rvec, gap_sym, delta0, plist, iperm)
        hamk = flibs.gen_ham(klist, ham_r, rvec)
        if option == CalcMode.CHIS_SPECTRUM_SC:
            print("generate qlist", flush=True)
            qlist, spa_length, xticks = plibs.mk_qlist(k_sets, Nx, Ny, Nz, bvec)
            sus, _, wlist = plibs.chis_spectrum_sc(mu, temp, Smat, hamk, deltak, klist, qlist, chiolist, Nw, Emax,
                                                   delta, sw_spsym)
            w, sp = np.meshgrid(wlist, spa_length)   # [Nq,Nw], like sus
            print("write chis spectrum in png file", flush=True)
            plt.contourf(sp, w, abs(sus.imag), 100)
            plt.colorbar()
            plt.hot()
            plt.savefig(fname='chis_sc_spec.png', dpi=300)
            plt.close()
        elif option == CalcMode.CHIS_QPOINT_SC:
            # --- BdG band structure along (0,0,0) -> (0,0.5,0) ---
            plots.plot_bdg_path(*plibs.bdg_bands_path(ham_r, S_r, rvec, mu, klist, gap_sym, delta0))
            # --- chi calculation ---
            q_point = np.array(at_point)
            print(f"q-point is ({q_point[0]:.3f}, {q_point[1]:.3f}, {q_point[2]:.3f})", flush=True)
            chis, chis_orb, wlist = plibs.chis_q_point_sc(q_point, hamk, deltak, mu, Emax, Nw, temp,
                                                          Smat, klist, chiolist, delta, sw_spsym)
            plt.plot(wlist, chis.imag)
            plt.savefig('chisq.png', dpi=150)
            plt.clf()
            # plt.show()
            with open(f"chis_sc.dat", "w") as f:
                for iw, ic in zip(wlist, chis):
                    f.write(f"{iw:5.3f}, {ic.imag:12.8f}, {ic.real:12.8f}\n")
            for cso in chis_orb.T:
                plt.plot(wlist, cso.imag)
            # plt.show()
            plt.savefig('chisq_orb.png', dpi=150)
            plt.close()
            with open(f"chis_scorb.dat", 'w') as f:
                for iw, ic in zip(wlist, chis_orb):
                    f.write(f"{iw:5.3f}, ")
                    for cso in ic:
                        f.write(f"{cso.imag:12.8f}, ")
                    f.write("\n")
        elif option == CalcMode.NMR_SC:  # Knight shift K(T) and 1/(T1 T) from the SC-state chi_s
            eig_n = flibs.get_eig(hamk)[0]
            if nmr_gap_file is not None:
                # Replace the gap_sym/delta0 form factor by the self-consistent RPA/FLEX gap.
                # Delta(R) is mesh independent, so this klist need not match the Eliashberg one.
                print(f'using the Eliashberg gap from {nmr_gap_file}.npz', flush=True)
                deltak = plibs.gap_from_eliashberg(klist, nmr_gap_file, nmr_gap_extrapolate,
                                                   nmr_gap_iw, nmr_gap_navg, eig=eig_n, mu=mu,
                                                   gap_max=nmr_gap_max)
            if nmr_spsym is not None:
                sw_spsym = bool(nmr_spsym)
            # deltak is read as the T=0 gap; nmr_sweep rescales its AMPLITUDE
            # by the BCS interpolation at each temperature and leaves the shape alone.
            gap0 = float(abs(deltak).max())
            # Tc must come from the gap ON THE FERMI SURFACE: max|Delta| over the whole zone
            # can be many times larger when the form factor has little weight on the sheets
            # that carry states (e.g. dx2-y2 on a pocket centred at (pi,pi)), and a Tc taken
            # from it puts the whole sweep at the wrong reduced temperature.
            _, gap_fs_max, _ = plibs.fermi_surface_gap(eig_n, mu, deltak)
            temp_c = plibs.bcs_temp_c(gap_fs_max) if nmr_tc is None else float(nmr_tc)
            if temp_c <= 0.:
                print("Error: Tc is zero -- the gap vanishes on the Fermi surface."
                      " Set nmr_tc, raise delta0, or pick another gap_sym", flush=True)
                return
            # sw_spsym comes from the sign of gap_sym (or nmr_spsym): it selects the sign of
            # the anomalous term, i.e. the channel preserved at T=0 (field perpendicular to
            # d) vs the Yosida-suppressed one (singlet, or field parallel to d).  With
            # nmr_gap_file the gap no longer comes from gap_sym, so set nmr_spsym explicitly.
            print(f'spin channel: {"preserved at T=0 (triplet, field _|_ d)" if sw_spsym else "Yosida-suppressed (singlet, or field || d)"}', flush=True)
            print(f'Delta(0) = {gap0*1e3:7.4f} meV (max on Fermi surface {gap_fs_max*1e3:7.4f} meV),'
                  f' Tc = {temp_c:10.3e} eV ({temp_c/kb:.2f} K)'
                  f'{" (BCS estimate from the Fermi-surface gap)" if nmr_tc is None else ""}', flush=True)
            qlist, qmult = plibs.nmr_qmesh(Nx, Ny, Nz, *nmr_qsize, sw_fold=nmr_qfold)
            hf_weight = qmult*plibs.hyperfine_form_factor(qlist, nmr_hf_A, nmr_hf_B)
            print(f'q-mesh {nmr_qsize} -> {len(qlist)} points'
                  f'{" after q/-q folding" if nmr_qfold else ""}, '
                  f'hyperfine A={nmr_hf_A:.3f} B={nmr_hf_B:.3f}', flush=True)
            temps = np.linspace(nmr_trange[0], nmr_trange[1], nmr_nt)*temp_c
            # Report the k-mesh/broadening hierarchy BEFORE the q-convergence probe, so a
            # warning about an under-resolved k-mesh is not buried under the probe output.
            plibs.nmr_scale_check(eig_n, mu, delta, gap0, deltak)
            plibs.nmr_cost_estimate(len(klist), len(qlist), nmr_nt, Norb, len(chiolist))
            if nmr_qconv:
                # Probe deep in the SC state, where chi''(q) has the most q-structure.
                plibs.nmr_qconvergence(hamk, deltak, klist, Nx, Ny, Nz, nmr_qsize, mu,
                                       max(temps[0], 0.1*temp_c), temp_c, Smat, chiolist,
                                       delta, nmr_w0, sw_spsym, nmr_hf_A, nmr_hf_B)
            print('start NMR temperature sweep', flush=True)
            res = plibs.nmr_sweep(hamk, deltak, klist, qlist, mu, temps, temp_c, Smat, chiolist,
                                  delta, nmr_w0, sw_spsym, hf_weight, sw_scale_check=False)
            plibs.write_nmr_dat(res)
            plots.plot_nmr(res)
            print('wrote nmr_sc.dat and nmr_sc.png', flush=True)
    elif option == CalcMode.LIN_ELIASHBERG and sw_nosym:  # symmetry-free gap equation (Zeeman / FFLO)
        if sw_self or sw_from_file:
            print('Warning: sw_nosym has no FLEX self-energy (TRS-free FLEX not implemented); using bare G', flush=True)
        slist, invs = get_spin_lists(Norb) if sw_soc else (None, None)
        plibs.calc_lin_eliash_nosym(Nx, Ny, Nz, Nw_mats, ham_r, S_r, rvec, mu, temp, chiolist, site, orb_dep, U, J,
                                    gap_sym, sw_soc, slist, invs,
                                    Umat if orb_dep else None, Jmat if orb_dep else None,
                                    hvec=zeeman_h, qlist=fflo_q, arnoldi_m=nosym_arnoldi_m, fill=fill)
    elif option in MODES_FLEX_ELIASH:  # flex/eliashberg calculations
        if sw_soc:  # with soc
            slist, invs = get_spin_lists(Norb)
            if option == CalcMode.FLEX:  # calc self-energy using flex
                plibs.calc_flex_soc(Nx, Ny, Nz, Nw_mats, ham_r, S_r, rvec, mu, temp, chiolist, slist, invs, site,
                                    orb_dep, U, J, fill, sw_out_self, sw_in_self,
                                    Umat if orb_dep else None, Jmat if orb_dep else None,
                                    m_diis=m_diis_num, sw_rescale=sw_rescale_flex,
                                    sigma_scale=sigma_in_scale)
            elif option == CalcMode.LIN_ELIASHBERG:  # calc gap function
                if epa_pair_file:
                    print('Warning: epa_pair_file (phonon pairing) is implemented without SOC only; ignored here',
                          flush=True)
                plibs.calc_lin_eliash_soc(Nx, Ny, Nz, Nw_mats, ham_r, S_r, rvec, mu, temp, chiolist, slist, plist, invs, site,
                                          orb_dep, U, J, fill, gap_sym, sw_self, sw_from_file, sw_out_self, sw_in_self,
                                          Umat if orb_dep else None, Jmat if orb_dep else None,
                                          sigma_scale=sigma_in_scale)
            elif option == CalcMode.GAP_FUNCTION:  # post gap calculation, output gap function/anomalous green's function
                sig = plibs.output_Fk(Nx, Ny, Nz, Nw_mats, ham_r, S_r, rvec, plist, mu, temp, sw_self, sw_soc, invs, slist,
                                      gap_sym)
                if sig is not None:
                    plots.plot_sigma_iw(*sig)
        else:  # without soc
            if option == CalcMode.FLEX:  # calc self-energy using flex
                plibs.calc_flex(Nx, Ny, Nz, Nw_mats, ham_r, S_r, rvec, mu, temp, chiolist, site,
                                orb_dep, U, J, fill, sw_out_self, sw_in_self,
                                Umat if orb_dep else None, Jmat if orb_dep else None,
                                m_diis=m_diis_num, sw_rescale=sw_rescale_flex, sw_tail=sw_tail_run,
                                sigma_scale=sigma_in_scale)
            elif option == CalcMode.LIN_ELIASHBERG:  # calc gap function
                plibs.calc_lin_eliashberg_eq(Nx, Ny, Nz, Nw_mats, ham_r, S_r, rvec, chiolist, site, plist, mu, temp, gap_sym, sw_self,
                                             orb_dep, U, J, fill, sw_from_file, sw_out_self, sw_in_self,
                                             Umat if orb_dep else None, Jmat if orb_dep else None,
                                             sw_tail=sw_tail_run, sigma_scale=sigma_in_scale, iperm=iperm,
                                             epa_pair_file=epa_pair_file, epa_pair_scale=epa_pair_scale)
            elif option == CalcMode.GAP_FUNCTION:  # post gap calculation, output gap function/anomalous green's function
                sig = plibs.output_Fk(Nx, Ny, Nz, Nw_mats, ham_r, S_r, rvec, plist, mu, temp, sw_self, iperm=iperm)
                if sig is not None:
                    plots.plot_sigma_iw(*sig)
    elif option == CalcMode.EILENBERGER:  # solve homogeneous quasiclassical Eilenberger equation
        # the orbital gap (RPA/FLEX export or explicit matrix) is projected band-resolved
        # onto the FULL 3D Fermi surface -> k_z dependence and multi-sheet signs come from
        # the orbital character, which is what an accidental/horizontal node needs
        gorb = (plibs.gap_orbital_from_wannier(eil_gap_file, eil_gap_iw, eil_gap_navg)
                if eil_gap_file else eil_gap_orbital)
        gkw = dict(gap_orbital=gorb, delta0=delta0, **set_eil_gauge(len(ham_r[0])))
        if eil_fs:  # model FS + Fermi velocity: anisotropic penetration depth (+ rho_zz if 3D)
            plibs.calc_fs_penetration(eil_coupling, temp, eil_wc, kind=(eil_fs_kind or 'ellipse'),
                                      gap_sym=gap_sym, params=eil_fs_params, kb=kb,
                                      nkz=eil_fs_nkz)
        elif eil_lambda:  # superfluid density / penetration depth lambda(T) (s exp-flat, d linear-in-T)
            plibs.calc_penetration_depth(eil_coupling, temp, eil_wc, gap_sym=gap_sym, kb=kb)
        elif eil_spin:  # spin-2x2 Zeeman response: singlet vs triplet d-vector (d||h vs d_|_h)
            plibs.calc_spin_pauli(Nx, Ny, Nz, eil_wc, ham_r, S_r, rvec, avec, mu, temp, eil_coupling,
                                  gap_sym=gap_sym, fs_width=eil_fs_width, kb=kb, **gkw)
        elif eil_pauli:  # Zeeman (Maki) Pauli-limiting sweep: singlet gap Delta(h), spinodal, Zeeman-split DOS
            plibs.calc_pauli_limit(Nx, Ny, Nz, eil_wc, ham_r, S_r, rvec, avec, mu, temp, gap_sym, eil_coupling,
                                   fs_width=eil_fs_width, kb=kb, **gkw)
        elif eil_free_energy:  # condensation free energy (Omega_s-Omega_n)/N0 vs T (coupling-independent)
            # the 3D FS gap (any harmonic incl. kz-dependent/chiral, or the projected orbital
            # gap) -> the condensation energy that selects between degenerate partners
            _, klist, eig, uni, _ = plibs.get_emesh(Nx, Ny, Nz, ham_r, S_r, rvec, avec, sw_uni=True)
            wf, phif = plibs.build_fs(eig, klist, mu, gap_sym, eil_fs_width, uni=uni, **gkw)
            plibs.calc_free_energy(eil_coupling, temp, eil_wc, kb=kb, phi=phif, wnf=wf)
        else:
            # the orbital gap (RPA/FLEX export or explicit matrix) is projected band-resolved
            # onto the FULL 3D Fermi surface -> k_z dependence and multi-sheet signs come from
            # the orbital character, which is what an accidental/horizontal node needs
            plibs.calc_eilenberger(Nx, Ny, Nz, eil_wc, ham_r, S_r, rvec, avec, mu, temp, gap_sym, eil_coupling,
                                   imp_gamma=eil_imp_gamma, imp_c=eil_imp_c, fs_width=eil_fs_width, kb=kb,
                                   method=eil_method, sw_find_tc=eil_find_tc, sw_imp_sweep=eil_imp_sweep,
                                   imp_sweep=eil_imp_list, **gkw)
    elif option == CalcMode.EILENBERGER_SURFACE:  # specular surface state via Riccati Eilenberger (model FS)
        if eil_surf_dvector:  # self-consistent triplet d-vector texture (spin-matrix Riccati)
            plibs.calc_surface_dvector(eil_coupling, temp, eil_wc, kb=kb, sub_ratio=eil_dvec_subratio, sw_ldos=eil_ldos)
        else:
            sfs, sfk = get_eil_gap_fs(rvec, ham_r, S_r, avec, Arot, temp, fill)
            plibs.calc_surface(eil_coupling, temp, eil_wc, gap_sym=gap_sym, beta_surf=eil_surf_beta,
                               kb=kb, sw_ldos=eil_ldos, imp_gamma=eil_imp_gamma, imp_c=eil_imp_c, h=eil_zeeman,
                               fs_kind=sfk, fs_params=eil_fs_params, fs=sfs)
    elif option == CalcMode.EILENBERGER_VORTEX:  # vortex / vortex lattice via Riccati Eilenberger (model FS)
        eil_fs_obj, eil_fs_kw = get_eil_gap_fs(rvec, ham_r, S_r, avec, Arot, temp, fill)
        if eil_fs_obj is not None and eil_fs_traj is not None:  # one chord integration per FS point
            tj = [eil_fs_traj] if np.isscalar(eil_fs_traj) else list(eil_fs_traj)
            nbeta = tj[0]; nv = tj[1] if len(tj) > 1 else 4; nphi = tj[2] if len(tj) > 2 else 8
            # the reduction collapses onto the vortex plane, so it needs the field direction
            eil_fs_obj = plibs.reduce_fs_trajectories(eil_fs_obj, gap_sym, nbeta, nv, nphi,
                                                      bdir=eil_field_dir)
        if eil_vort_maxwell:  # self-consistent finite-kappa vector potential A(r) (Maxwell back-reaction)
            plibs.calc_vortex_maxwell(eil_coupling, temp, eil_wc, gap_sym=gap_sym, field=eil_field,
                                      kappa=eil_kappa, kb=kb, Lxi=eil_vort_lxi, ngrid=eil_vort_ngrid)
        elif eil_vort_current:  # circulating charge supercurrent j_phi(rho) of an isolated vortex
            plibs.calc_vortex_current(eil_coupling, temp, eil_wc, gap_sym=gap_sym, kb=kb,
                                      Lxi=eil_vort_lxi, ngrid=eil_vort_ngrid)
        elif eil_vort_chiral:  # chiral (p+ip / d+id) vortex: complex two-component amplitude, core-induced opposite chirality
            plibs.calc_vortex_chiral(eil_coupling, temp, eil_wc, ell=eil_chiral_ell, m_dom=eil_chiral_m,
                                     dvec=eil_chiral_dvec, kb=kb, Lxi=eil_vort_lxi,
                                     ngrid=eil_vort_ngrid, fs=eil_fs_obj, bdir=eil_field_dir)
        elif eil_vort_dvector:  # self-consistent triplet d-vector vortex/lattice (spin-matrix Riccati)
            dfs = (plibs.build_wannier_fs(rvec, ham_r, S_r, avec, plibs.get_mu(ham_r, S_r, rvec, Arot, temp, fill))
                   if eil_fs_kind == 'wannier' else None)  # FS geometry only (d-vector channels carry the gap)
            if eil_vort_lattice_sc:  # je-style TRUE periodic d-vector lattice (formulation A, square/triangular)
                plibs.calc_vortex_lattice_sc_dvector(eil_coupling, temp, eil_wc, kb=kb, field=eil_field,
                                                     lattice=eil_lattice, sub_ratio=eil_dvec_subratio,
                                                     kappa=(None if eil_kappa >= 1e3 else eil_kappa), fs=dfs)
            else:  # isolated vortex (eil_field=0) or circular-cell lattice (eil_field>0)
                plibs.calc_vortex_dvector(eil_coupling, temp, eil_wc, kb=kb, sub_ratio=eil_dvec_subratio,
                                          field=eil_field, fs=dfs, bdir=eil_field_dir)
        elif eil_vort_lattice_sc and eil_field_list is not None:  # je-style self-consistent periodic lattice (formulation A); eil_lattice square/triangular; eil_nvortex=Vw flux quanta/cell; finite eil_kappa = London A back-reaction, >=1e3 = bare extreme
            plibs.calc_vortex_lattice_sc(eil_coupling, temp, eil_wc, gap_sym=gap_sym,
                                         field_list=eil_field_list, lattice=eil_lattice, kb=kb, fs=eil_fs_obj,
                                         kappa=(None if eil_kappa >= 1e3 else eil_kappa), Vw=eil_nvortex,
                                         self_consistent_A=eil_vort_scA, bdir=eil_field_dir, h=eil_zeeman)
        elif eil_field_list is not None:  # sweep B/Hc2 on the TRUE periodic lattice -> <N(0)>(B) (d~sqrt(B) Volovik)
            plibs.calc_vortex_lattice_periodic(eil_coupling, temp, eil_wc, gap_sym=gap_sym,
                                               field_list=eil_field_list, kappa=eil_kappa, lattice=eil_lattice, kb=kb,
                                               fs_kind=eil_fs_kw, fs_params=eil_fs_params, fs=eil_fs_obj,
                                               nflux=eil_nvortex, bdir=eil_field_dir)
        else:  # single field (isolated vortex if eil_field=0, else circular-cell lattice)
            plibs.calc_vortex(eil_coupling, temp, eil_wc, gap_sym=gap_sym, kb=kb, sw_ldos=eil_ldos,
                              imp_gamma=eil_imp_gamma, imp_c=eil_imp_c, field=eil_field, h=eil_zeeman,
                              kappa=(eil_kappa if eil_vort_field else 0.0), tilt_deg=eil_vort_tilt,
                              fs_kind=eil_fs_kw, fs_params=eil_fs_params, fs=eil_fs_obj,
                              Lxi=eil_vort_lxi, ngrid=eil_vort_ngrid, bdir=eil_field_dir)
    elif option == CalcMode.BAND_FILLING:  # per-band filling and electron/hole character
        print('band, hole fraction, electron fraction:', flush=True)
        plibs.report_band_filling(Nx, rvec, ham_r, S_r, mu, Arot)
    elif option == CalcMode.CYCLOTRON_MASS:  # calc cyclotron mass
        plibs.get_mass(Nx, rvec, ham_r, S_r, avec, bvec, mu, RotMat, units)
    elif option == CalcMode.DHVA:  # plot dHvA frequency vs angle
        theta_list = np.linspace(0., 90., 40)
        dhva = plibs.get_dhva_band(Nx, rvec, ham_r, S_r, avec, bvec, mu, theta_list, units)
        if dhva:
            plots.plot_dhva(dhva, phi=0.)
    elif option == CalcMode.ELECTRON_MASS:  # band (effective) mass along the symmetry line
        res = plibs.band_mass_path(k_sets, kmesh, rvec, ham_r, S_r, avec, units)
        plibs.report_band_mass(res, xlabel, mu)
        plots.plot_band_mass(res, xlabel)
    elif option == CalcMode.SPECTRUM_IMPURITY:  # unfolded spectrum of an impurity supercell
        Nsite = int(np.prod(imp_cell))
        sites = (plibs.random_impurity_sites(Nsite, imp_conc, imp_seed) if imp_sites is None
                 else np.asarray(imp_sites, dtype=np.int64))
        res = plibs.calc_impurity_spectrum(ham_r, rvec, imp_cell, sites, imp_V, fill, temp, k_sets, kmesh, bvec,
                                           Emin, Emax, Nw, delta)
        plots.plot_impurity_spectrum(res, xlabel)
    elif option == CalcMode.SIGMA_CPA:
        # --- CPA calculation ---
        mu = plibs.get_mu(ham_r, S_r, rvec, Arot, temp, fill)
        print(f'Temperature = {temp:10.3e} eV ({temp/kb:.2f} K)', flush=True)
        print(f'chem. pot. = {mu:.4f} eV', flush=True)
        # binary alloy: species B shifts every on-site level by 1 eV, concentration 0.5
        spa_length, xticks, wlist, Akw = plibs.calc_cpa(ham_r, rvec, mu, temp, Nx, Ny, Nz, Nw, Emin, Emax, delta,
                                                        k_sets, kmesh, bvec, x_cpa=0.5, VB_shift=1.0)
        plots.plot_cpa_spectrum(spa_length, xticks, xlabel, wlist, Akw)
    elif option == CalcMode.NONLIN_ELIASHBERG:
        plibs.calc_eliashberg_eq(Nx, Ny, Nz, Nw_mats, ham_r, S_r, rvec, chiolist, site, plist,
                                 mu, temp, gap_sym, sw_self,
                                 orb_dep, U, J, fill, sw_from_file, sw_out_self, sw_in_self,
                                 Umat if orb_dep else None, Jmat if orb_dep else None,
                                 m_diis=m_diis_num, sw_rescale=sw_rescale_flex,
                                 sw_check_only=sw_check_only, sw_tail=sw_tail_run,
                                 sigma_scale=sigma_in_scale, iperm=iperm)


if __name__ == "__main__":
    main()
__license__ = "MIT"
