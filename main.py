#!/usr/bin/env python
#-*- coding:utf-8 -*-
from enum import IntEnum
"""
ftype: set input hamiltonian's format
0: ham_r.txt, irvec.txt, ndegen.txt in {fname} dir
1: .input file named {fname}
2: {fname}_hr.dat file (wannier90 default hopping file
3: for non-orthogonal basehoppings as MLO basis
else: Hopping.dat file (ecalj hopping file)
brav: choose primitive translation vector S,FC,BC etc
0: simple 
1: face center (for QE default)
2: body center (for QE default)
3: hexagonal
4: trigonal    (for QE sbrav==5)
5: base center
6: face center (common)
7: body center (common)
else: monoclinic
"""

#fname,ftype,brav,sw_soc='inputs/Sr2RuO4nso',0,7,False
#fname,ftype,brav,sw_soc='inputs/Sr2RuO4',2,2,True
#fname,ftype,brav,sw_soc='inputs/SiMLO.input',3,6,False
#fname,ftype,brav,sw_soc='inputs/NdFeAsO.input',1,0,False
#fname,ftype,brav,sw_soc='inputs/000AsP.input',1,0,False
#fname,ftype,brav,sw_soc='inputs/FeS',2,0,False
fname,ftype,brav,sw_soc='inputs/hop2.input',1,0,False
#fname,ftype,brav,sw_soc='inputs/hop2_soc.input',1,0,True
#fname,ftype,brav,sw_soc='inputs/square.hop',1,0,False
#fname,ftype,brav,sw_soc='inputs/square_soc.hop',1,0,True

sw_dec_axis=False

class CalcMode(IntEnum):
    """Calculation modes. The second tuple element is the human-readable label
    printed by main(); it is the single source of truth (replaces the old `opstr` list)."""
    def __new__(cls,value,description):
        obj=int.__new__(cls,value)
        obj._value_=value
        obj.description=description
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
    CARRIER_NUM        = (17, "calculate carrier number")                  # carrier number calculation
    CYCLOTRON_MASS     = (18, "calculate cyclotron mass")                  # cyclotron mass calculation
    DHVA               = (19, "plot dHvA frequency")                       # dHvA frequency vs angle plot (not implemented)
    ELECTRON_MASS      = (20, "calculate electron mass")                   # electron mass calculation (not implemented)
    SPECTRUM_IMPURITY  = (21, "spectrum with impurity")                    # spectral function with impurity (not implemented)
    SIGMA_CPA          = (22, "calculate sigma_cpa")                       # conductivity via CPA
    NONLIN_ELIASHBERG  = (23, "solve nonlinear eliashberg equation")       # solve nonlinear Eliashberg equation (not implemented)
    EILENBERGER        = (24, "solve quasiclassical Eilenberger equation")  # homogeneous multi-orbital Eilenberger (Matsubara)
    EILENBERGER_SURFACE= (25, "solve surface state via Riccati Eilenberger")  # specular surface gap profile + LDOS (model FS)
    EILENBERGER_VORTEX = (26, "solve isolated vortex via Riccati Eilenberger")  # vortex D(rho) + core LDOS (model FS)

class ColorMode(IntEnum):
    """Color modes for band/FS plots; second element is the printed label (replaces the old `cstr` list)."""
    def __new__(cls,value,description):
        obj=int.__new__(cls,value)
        obj._value_=value
        obj.description=description
        return obj
    MONO     = (0, "no color")
    ORBITAL  = (1, "orbital weight")
    VELOCITY = (2, "velocity size")

#----- CalcMode capability groups: single source for the option-set checks in main() -----
M=CalcMode
# modes that need a symmetry-line k-path (k_sets / xlabel)
MODES_SYMLINE       = frozenset({M.BAND,M.SPECTRUM,M.CHIS_SPECTRUM,M.PHI_SPECTRUM,M.CHIS_SPECTRUM_SC,
                                 M.FLEX,M.ELECTRON_MASS,M.SPECTRUM_IMPURITY,M.SIGMA_CPA})
# modes that need a rotation matrix (RotMat) for the 2D / extremal-orbit cut
MODES_NEED_ROTMAT   = frozenset({M.FERMI_2D,M.CYCLOTRON_MASS})
# FLEX / Eliashberg family (self-energy, linearized gap, post gap)
MODES_FLEX_ELIASH   = frozenset({M.FLEX,M.LIN_ELIASHBERG,M.GAP_FUNCTION})
# chi at a single q-point (normal & SC state)
MODES_CHIS_QPOINT   = frozenset({M.CHIS_QPOINT,M.CHIS_QPOINT_SC})
# modes that color-code bands/FS by orbital or velocity
MODES_COLOR         = frozenset({M.BAND,M.FERMI_2D,M.FERMI_3D})
# susceptibility / Eliashberg modes that build the orbital-pair basis (chiolist)
MODES_NEED_CHI      = frozenset({M.CHIS_SPECTRUM,M.CHIS_QPOINT,M.CHIS_QMAP,M.PHI_SPECTRUM,M.PHI_QMAP,
                                 M.CHIS_SPECTRUM_SC,M.CHIS_QPOINT_SC,M.FLEX,M.LIN_ELIASHBERG,M.NONLIN_ELIASHBERG})
# susceptibility modes that need the static Coulomb vertex (S/C matrices)
MODES_COULOMB_VERTEX= frozenset({M.CHIS_SPECTRUM,M.CHIS_QPOINT,M.CHIS_QMAP,M.CHIS_SPECTRUM_SC,M.CHIS_QPOINT_SC})
# modes that print a constant U,J header
MODES_PRINT_UJ      = frozenset({M.CHIS_SPECTRUM,M.CHIS_QPOINT,M.CHIS_QMAP,M.FLEX,M.LIN_ELIASHBERG})
# modes that compute/define mu themselves (skip the common mu block)
MODES_SELF_MU       = frozenset({M.CONDUCTIVITY_BT,M.CONDUCTIVITY_PT,M.SPECTRUM_IMPURITY,M.SIGMA_CPA,M.EILENBERGER_SURFACE,M.EILENBERGER_VORTEX})
# k-mesh reporting groups
MODES_KMESH_SYMLINE = frozenset({M.BAND,M.SPECTRUM})
MODES_KMESH_SINGLE  = frozenset({M.FERMI_2D,M.FERMI_3D,M.CHIS_QMAP,M.PHI_QMAP,M.CARRIER_NUM,M.CYCLOTRON_MASS,M.EILENBERGER})
MODES_MATSUBARA     = frozenset({M.FLEX,M.LIN_ELIASHBERG})
# dispatch groups
MODES_CHI_NORMAL    = frozenset({M.CHIS_SPECTRUM,M.CHIS_QPOINT,M.CHIS_QMAP,M.PHI_SPECTRUM,M.PHI_QMAP})
MODES_CHIS_SC       = frozenset({M.CHIS_SPECTRUM_SC,M.CHIS_QPOINT_SC})
del M

#option=CalcMode.CHIS_QPOINT_SC
option=CalcMode.LIN_ELIASHBERG  #calculation mode to run (see the CalcMode enum above; 0-23 RPA/FLEX/transport, 24-26 Eilenberger superconductivity)
color_option=ColorMode.VELOCITY  #band/FS coloring (option 0,2,3): MONO=black, ORBITAL=olist weights->RGB, VELOCITY=|v(k)|

#Nx,Ny,Nz,Nw=256,256,4,200 #k and energy(or matsubara freq.) mesh size
Nx,Ny,Nz,Nw=32,32,2,512
kmesh=200               #number of k-points along the symmetry line for band/spectrum plots (larger=smoother)
kscale=[1.0,1.0,1.0]    #per-axis display scale for the 3D Fermi-surface plot (option 3); e.g. [1,1,0.5] compresses kz
kz=0.0                  #reduced kz of the 2D Fermi-surface cut (option 2): 0=Gamma plane, 0.5=zone-boundary plane
#RotMat=[[0,0,1],[0,1,0],[1,0,0]]

#abc=[3.96*0.70711,3.96*0.70711,13.02*.5]
abc=[3.68,3.68,5.03]    #lattice constants a,b,c [Angstrom] (group velocities & symmetry-path lengths)
#alpha_beta_gamma=[90.,90.,90]  #lattice angles alpha,beta,gamma [deg] (default 90,90,90 if undefined)
#temp=2.0e-2 #2.59e-2   #directly set k_B*T [eV]; if defined it overrides tempK
tempK=85 #Kelvin        #temperature [K] (converted internally to temp=k_B*tempK [eV])
fill= 1.0 #2.9375       #band filling; mu solved from sum f(eps-mu)=Nk*fill (no SOC: per spin, full=Norb; SOC: total, full=2*Norb)
#site_prof=[5]

Emin,Emax=-3,1.         #energy window [eV] for DOS / spectral-function plots (option 1,4)
delta=5.0e-3            #spectral broadening eta [eV]: imaginary part added to G (Lorentzian width); too large smears, too small=noise
Ecut=1.0e-2            #fixed energy omega_0 [eV] for the q-space susceptibility maps (option 9,11); ~0 probes the Fermi surface
tau_const=100          #constant relaxation time tau [fs] for Boltzmann transport (option 5)
olist=[0,0,0]          #orbital indices mapped to R,G,B for orbital-weight coloring (color_option=1); nested lists group orbitals
#olist=[0,[1,2],3]
#U,J=0.,0.
#U,J= 0.2, 0.025
#U,J= 0.4, 0.05
#U,J= 0.6, 0.075
U,J=1.2,0.15           #on-site Hubbard U and Hund J [eV] (FLEX/RPA); screened U'=U-2J used automatically
#U,J=1.8,0.225
#0:s,1:dx2-y2,2:spm,3:dxy,-1:px,-2:py,-3:p+ip  (also drives ALL eilenberger modes; model FS/cylinder maps the int -> continuum harmonic, 2 spm -> s)
gap_sym=1

#use calculation of magnetic susceptibility at superconducting state
#delta0=1.e-2 #maximum gap size for calculating susceptibility in SC state; set to 0 for normal state
d0=1.e-1               #helper amplitude scale [eV] for building the per-band delta0 list below
delta0=[0.,d0*2.,d0*3.,-d0,0.]  #initial SC gap amplitude/sign per band [eV] for SC-state chi (option 12,13); float=single shape, list=multi-gap (signs->s+-)
#mu0=9.85114560061123
#k_sets=[[0., 0., 0.],[.5, 0., 0.],[.5, .5, 0.]]
#xlabel=[r'$\Gamma$','X','M']
#m_diis_num=2
at_point=[ 0., .5, 0.]  #reduced q-point [0..1] for the single-q susceptibility (option 8)
orb_dep=False  #use orbital dependence U,J (True: use Umat/Jmat matrices below; False: constant U,J)
Umat=None      #orbital-dependent U matrix (Norb x Norb); set when orb_dep=True
Jmat=None      #orbital-dependent J matrix (Norb x Norb); set when orb_dep=True
sw_unit=True    #True: use physical constants (SI/eV units), False: set all constants to 1 (dimensionless test mode)
sw_tdf=False   #True: compute the transport distribution function first, then energy-integrate (energy-dependent tau)
sw_omega=False #True: real freq, False: Matsubara freq.
sw_rescale_flex=True #True: rescale self energy to make max|Sigma|~U, False: no rescaling
sw_self=False  #True: use calculated self energy for spectrum band plot
sw_out_self=True #True: write the FLEX self-energy to sigma.bin/self_en.npz (also triggers gap output in option 15)
sw_in_self=False #True: load the previous self-energy from sigma.bin as the initial guess for the SC loop
sw_from_file=False #True: read the self-energy from sigma.bin and skip FLEX (solve Eliashberg with it directly)
sw_check_only=False #True: stop after linear Eliashberg (also stops if Stoner factor>=1 or lambda<1)
#----- EILENBERGER (homogeneous quasiclassical) parameters -----
eil_coupling=0.3     #dimensionless separable pairing coupling lambda (with <phi^2>_FS=1)
eil_wc=0.5           #fixed Matsubara cutoff energy [eV] (sets the pairing scale / Tc)
eil_imp_gamma=0.0    #non-magnetic impurity scattering strength Gamma [eV] (0=clean)
eil_imp_c=1.0e8      #T-matrix cot(delta0): large=Born limit, 0=unitary limit
eil_fs_width=5.0e-3  #Gaussian Fermi-surface broadening [eV]
eil_method='normalization' #homogeneous (g,f) route: 'normalization' (fast) or 'riccati' (vortex-lattice-ready)
eil_find_tc=False    #True: bisect for Tc at the current impurity setting
eil_imp_sweep=False  #True: sweep Gamma and write Tc(Gamma) to eilenberger_tc.dat
eil_imp_list=None    #array of Gamma values [eV] for the sweep (e.g. np.linspace(0,0.05,11))
eil_pauli=False      #True: Zeeman/Maki Pauli-limiting sweep (singlet gap Delta(h), spinodal, Zeeman-split DOS)
eil_free_energy=False #True: condensation free energy (Omega_s-Omega_n)/N0 vs T (coupling-independent; writes free_energy.dat)
eil_spin=False       #True: spin-2x2 Zeeman response, singlet vs triplet d-vector (d||h Pauli-limited, d_|_h immune)
eil_lambda=False     #True: superfluid density rho_s(T)/penetration depth lambda(T) sweep (s exp-flat, d linear-in-T)
eil_fs=False         #True: model-FS + Fermi-velocity penetration depth (anisotropic lambda_xx/lambda_yy)
eil_fs_kind=None     #Fermi surface for ALL eilenberger modes (homogeneous penetration / surface / vortex): None=isotropic (cylinder; homogeneous falls back to 'ellipse'), 'iso'/'ellipse'(params=(mx,my))/'tb'(params=t) model FS, or 'wannier' (loaded band FS+v_F; surface/vortex symmetry/multiband from gap_sym,delta0)
eil_fs_params=(1.0,0.4) #model-FS parameters (ellipse masses or tb hopping)
eil_zeeman=0.0       #Zeeman (Maki) field [eV] for the LDOS (surface: splits the d[110] ZEBS into +-h; vortex: spin-splits the core bound states)
#----- EILENBERGER inhomogeneous (surface & vortex, Riccati; model cylindrical FS) -----
eil_ldos=True            #True: also compute the real-frequency LDOS (bound/core states) (surface & vortex)
eil_surf_beta=0.785398   #surface orientation [rad]: 0=[100], pi/4(0.7854)=[110] (d-wave ZEBS)
eil_surf_dvector=False   #True: self-consistent triplet d-vector TEXTURE at the surface (dominant p_x(e_x) + subdominant p_y(e_z), spin-matrix Riccati)
eil_dvec_subratio=0.9    #subdominant/dominant coupling ratio for the d-vector texture (~0.85 is the bulk threshold)
eil_vort_lxi=8.0         #vortex cell half-width in coherence lengths xi (isolated vortex, field=0)
eil_vort_ngrid=81        #vortex 2D grid points per axis
eil_vort_field=False     #True: also compute the self-consistent finite-kappa Maxwell field profile B(rho) of the vortex (uses eil_kappa)
eil_vort_dvector=False   #True: self-consistent triplet d-vector TEXTURE around the vortex core (dominant p_x(e_x) winding + core-localized subdominant p_y(e_z), 2D spin-matrix Riccati; uses eil_dvec_subratio)
eil_vort_current=False   #True: circulating charge supercurrent j_phi(rho) of an isolated vortex (writes vortex_current.dat)
eil_vort_maxwell=False   #True: circular-cell vortex with the self-consistent finite-kappa vector potential A(r) (Maxwell back-reaction; needs eil_field>0, uses eil_kappa)
eil_vort_lattice_sc=False #True: je-style self-consistent TRUE periodic lattice (formulation A, extreme type-II): complex Psi(r) with a real node at every core + full Abrikosov supercurrent; sweeps eil_field_list -> <N(0)>(B) (d~sqrt(B) Volovik)
eil_vort_scA=False       #True (lattice_sc, finite eil_kappa): fully self-consistent vector potential A from the quasiclassical current j_s=<v_F Im g> (je A_renew), instead of the analytic London A
eil_gap_orbital=None     #None, or an Norb x Norb orbital-basis pair-potential matrix (or callable kfrac->NxN): the band gap is its low-energy PROJECTION onto the FS bands (Nagai-Nakamura, JPSJ 85 074707 (2016) Eq.43; needs wannier FS); supersedes gap_sym/delta0
eil_gap_file=None        #None, or the base name (no extension) of an RPA/FLEX gap exported by LIN_ELIASHBERG/NONLIN_ELIASHBERG with sw_out_self=True (output_gap_wannier -> 'gap_wannier'): loads Delta(R,iw) and uses it as eil_gap_orbital (projected to the FS bands); supersedes eil_gap_orbital/gap_sym. Use a self-consistent RPA gap (e.g. KFe2As2) as the vortex pairing form factor
eil_gap_iw=0             #starting Matsubara index for eil_gap_file (0=lowest iw_0; sharpest gap symmetry/anisotropy, matches the FS gap usually quoted)
eil_gap_navg=1           #number of consecutive Matsubara slices to average for eil_gap_file (1=single iw_gap_iw slice; >1 smooths noise at the cost of slightly diluting the anisotropy)
eil_vort_tilt=0.0        #field tilt theta [deg] from the c-axis (quasi-2D): orbital uses B_z=B cos(theta), Zeeman -> h/cos(theta) (Pauli/orbital ratio)
eil_nvortex=1            #vortices (flux quanta) per computational cell of the periodic lattice (supercell; n^2 reduces to the primitive cell)
eil_field=0.0            #vortex lattice field B/Hc2 (0=isolated vortex; >0=circular-cell lattice w/ Doppler)
eil_field_list=None      #list of B/Hc2 to sweep <N(0)>(B) on the TRUE periodic lattice (e.g. [0.04,0.08,0.16,0.32]); None=single field
eil_kappa=100.0          #GL kappa=lambda/xi for the periodic lattice (large=extreme type-II; finite=London screening/Maxwell)
eil_lattice='square'     #periodic vortex lattice geometry: 'square' or 'triangular'
#------------------------ initial parameters are above -------------------------------
#----------------------------------main functions-------------------------------------
#-------------------------------- import packages ------------------------------------
import os
import numpy as np
import libs.flibs as flibs, libs.plibs as plibs
import scipy.linalg as sclin, scipy.constants as scconst
import matplotlib.pyplot as plt, matplotlib.cm as cm
#----------------------------- initial settings --------------------------------------
sw_gen_sym=False
try:
    kz
except NameError:
    kz=0.0
try:
    mu0
    sw_calc_mu=False
except NameError:
    sw_calc_mu=True
try:
    kscale
except NameError:
    kscale=1.0
try:
    abc
except NameError:
    abc=[1.,1.,1.]
    print('abc not found, set to 1,1,1')
try:
    alpha_beta_gamma
except NameError:
    alpha_beta_gamma=[90.,90.,90]
    print('alpha_beta_gamma not found, set to 90,90,90')
alatt=np.array(abc)
deg=np.array(alpha_beta_gamma)
if option in MODES_SYMLINE: #modes using k_sets/xlabel (symmetry line)
    try:
        k_sets
        xlabel
    except NameError:
        sw_gen_sym=True
        k_sets,xlabel2=plibs.get_symm_line(brav)
        try:
            xlabel
            if len(xlabel)!=len(k_sets):
                xlabel=xlabel2
        except NameError:
            xlabel=xlabel2
if sw_unit:
    # Physical constants in SI/eV units; ihbar converts velocity: hbar/(eV·s) * 1e-10 (AA->m)^-1
    hbar=scconst.physical_constants['Planck constant over 2 pi in eV s'][0]
    ihbar=1.0e-10/hbar   # factor for v = (1/hbar) * dE/dk  [m/s per eV/AA^-1]
    kb=scconst.physical_constants['Boltzmann constant in eV/K'][0]
    eC=scconst.e         # elementary charge [C]
    tau_unit=1.e-15      # relaxation time unit: 1 fs = 1e-15 s
    emass=scconst.m_e    # electron rest mass [kg]
else:
    # Dimensionless units: all constants set to 1 for testing or code-level checks
    hbar=1.
    ihbar=1.
    kb=1.
    eC=1.
    tau_unit=1.
    emass=1.
try:
    tempK
    try:
        temp
        print('temp & tempK find, use temp')
    except NameError:
        temp=tempK*kb
except NameError:
    print('tempK not found')
    exit()
try:
    m_diis_num
except NameError:
    m_diis_num=5
if option in MODES_NEED_ROTMAT:
    try:
        RotMat
        RotMat=np.array(RotMat)
    except NameError:
        print('No RotMat')
        RotMat=np.eye(3)
#------------------------ define functions -------------------------------------------
def print_matrix(label,mat,width=10,prec=3,ndigits=10):
    """Print a labelled N x 3 real matrix, one row per line (transport tensors)."""
    print(label,flush=True)
    for row in np.asarray(mat).round(ndigits):
        print(f" {row[0]:{width}.{prec}e} {row[1]:{width}.{prec}e} {row[2]:{width}.{prec}e}",flush=True)

def flatten_orbs(seq):
    """Flatten an olist (mix of ints and int-lists) into a flat list of orbital indices."""
    out=[]
    for o in seq:
        out+=o if isinstance(o,list) else [o]
    return out

def plot_band(eig,spl,xlabel,xticks,uni,ol,color):
    def get_col(cl,ol):
        col=(np.abs(cl[ol])**2 if isinstance(ol,int)
             else (np.abs(cl[ol])**2).sum(axis=0)).round(4)
        return col
    # Check if olist has required 3 elements
    if len(ol) < 3:
        print(f"Warning: plot_band olist requires 3 elements (current: {len(ol)}). Color display disabled",flush=True)
        color=False
    fig=plt.figure()
    ax=plt.axes()
    for e,cl in zip(eig,uni): #band loop
        if color:
            # Normalize eigenvectors so that per-orbital weights (|u_n|^2) sum to 1
            norm=np.sqrt((abs(cl)**2).sum(axis=0))
            # Check for norm close to zero (use tolerance for floating-point comparison)
            if np.any(norm < 1e-14):
                print("Warning: plot_band found norm close to zero. Skipping",flush=True)
                continue
            cls=cl/norm
            # Map three orbital groups to RGB channels for color coding
            c1=get_col(cls,ol[0])
            c2=get_col(cls,ol[1])
            c3=get_col(cls,ol[2])
            clist=np.array([c1,c2,c3]).T
            for i in range(len(e)):
                plt.plot(spl[i:i+2],e[i:i+2],c=clist[i])
        else:
            for i in range(len(e)):
                plt.plot(spl[i:i+2],e[i:i+2],color='black')
    for x in xticks[1:-1]:
        plt.axvline(x,ls='-',lw=0.25,color='black')
    plt.ylim(eig.min()*1.1,eig.max()*1.1)
    plt.xlim(0,spl.max())
    plt.axhline(0.,ls='--',lw=0.25,color='black')
    plt.xticks(xticks,xlabel)
    plt.show()

def plot_FS(fscolors,klist,color_option:int):
    fig=plt.figure()
    ax=fig.add_subplot(111,aspect='equal')
    col=['r','g','b','c','m','y','k','w']
    if color_option==0:
        for kl in klist:
            for kk in kl:
                plt.plot(kk[:,0],kk[:,1],color='black',lw=2)
    else:
        if color_option==2:
            v=[]
            k=[]
            for vl,kl in zip(fscolors,klist):
                for vv,kk in zip(vl,kl):
                    v.extend(vv)
                    k.extend(kk)
            v=np.array(v)
            k=np.array(k)
            vmax=v.max()
            vmin=v.min()
            # Check for vmax=vmin case
            if vmax == vmin:
                print(f"Warning: plot_FS has vmax=vmin={vmax}. Color map disabled",flush=True)
                vmax = vmin + 1.0  # Set default value
        for kl,fscol,cb in zip(klist,fscolors,col):
            for kk,fcol in zip(kl,fscol):
                if color_option==2:
                    clist=cm.jet((fcol-vmin)/(vmax-vmin))
                else:
                    clist=fcol
                for k1,k2,clst in zip(kk,kk[1:],clist):
                    plt.plot([k1[0],k2[0]],[k1[1],k2[1]],c=clst,lw=2)
        if color_option==2:
            plt.scatter(k[:,0],k[:,1],s=0.1,c=v)
            plt.jet()
            plt.colorbar()
    plt.xlim(-0.5,0.5)
    plt.ylim(-0.5,0.5)
    plt.xticks([-0.5,0,0.5],[r'-$\pi$','0',r'$\pi$'])
    plt.yticks([-0.5,0,0.5],[r'-$\pi$','0',r'$\pi$'])
    plt.show()

def plot_3d_surf(fspolys,fscenters,fscolors,surface_opt,kscale,bvec):
    from mpl_toolkits.mplot3d import axes3d
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    import matplotlib.colors as colors
    fig=plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    clist=['r','g','b','c','m','y','k','w']
    if isinstance(kscale,int):
        ks=kscale*np.array([1.,1.,1.])
    else:
        ks=np.array(kscale)
    if surface_opt==0:
        for j,polys in enumerate(fspolys):
            ax.add_collection3d(Poly3DCollection(polys,facecolor=clist[j%6]))
    else:
        if surface_opt==1:
            tri=Poly3DCollection(fspolys,facecolors=fscolors,lw=0)
        else:
            clmax=fscolors.max()
            clmin=fscolors.min()
            # Check for clmax=clmin case
            if clmax == clmin:
                print(f"Warning: plot_3d_surf has clmax=clmin={clmax}. Using default colors",flush=True)
                nor_cols=np.ones_like(fscolors)*0.5
            else:
                nor_cols=(fscolors-clmin)/(clmax-clmin)
            tri=Poly3DCollection(fspolys,facecolors=cm.jet(nor_cols),lw=0)
            fs=ax.scatter(fscenters[:,0],fscenters[:,1],fscenters[:,2],
                          c=fscolors,cmap=cm.jet,s=0.1)
            plt.colorbar(fs,format='%.2e')
        ax.add_collection3d(tri)
    ax.grid(False)
    ax.set_xlim(-np.pi*ks[0], np.pi*ks[0])
    ax.set_ylim(-np.pi*ks[1], np.pi*ks[1])
    ax.set_zlim(-np.pi*ks[2], np.pi*ks[2])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    plt.tight_layout()
    plibs.BZedge(bvec, ax)
    plt.show()
    #plt.savefig(fname='3DFS.png',dpi=300)
    
def set_init_3dfsplot(color_option,polys,centers,blist,avec,rvec,ham_r,S_r,olist):
    if color_option==0:
        fspolys=polys
        fscenters=[]
        fscolors=[]
        return fspolys,fscenters,fscolors
    else:
        colors=plibs.get_colors(centers,blist,ihbar*avec.T,rvec,ham_r,S_r,olist,color_option)
        fspolys=[]
        fscenters=[]
        fscolors=[]
        for i in polys:
            fspolys.extend(i)
        for i in centers:
            fscenters.extend(i)
        for i in colors:
            fscolors.extend(i)
        fscenters=np.array(fscenters)*2.*np.pi
        fscolors=np.array(fscolors)
        return fspolys,fscenters,fscolors

def plot_spectrum(k_sets,xlabel,kmesh,bvec,mu:float,ham_r,S_r,rvec,Emin:float,Emax:float,
                  delta:float,Nw:int,sw_self=True):
    klist,spa_length,xticks=plibs.mk_klist(k_sets,kmesh,bvec)
    eig,uni=plibs.get_eigs(klist,ham_r,S_r,rvec)
    wlist=np.linspace(Emin,Emax,Nw)
    if sw_self:
        # Sigma is defined on the FLEX k-mesh, so evaluate it on the k-path by Fourier
        # interpolation of the Wannier-R self-energy written by output_self_wannier.
        try:
            npz=np.load('self_en_wannier.npz')
        except FileNotFoundError:
            print("Error: 'self_en_wannier.npz' not found (run FLEX with sw_out_self=True)",flush=True)
            return
        sig_r,rvec_s,iw_grid=npz['sigma'],npz['rvec'],npz['iw']
        mu_self=float(npz['mu'])
        print(f'chem. pot. with self = {mu_self:.4f} eV',flush=True)
        # Sigma(k) = sum_R Sigma(R) e^{-2pi i k.R}; Sigma is not Hermitian, so gen_ham cannot be reused
        phase=np.exp(-2j*np.pi*klist.dot(rvec_s.T))                          # [Nks,Nr]
        selfen=np.ascontiguousarray(np.einsum('mnwr,kr->mnwk',sig_r,phase))  # [Norb,Norb,Niw,Nks]
        ham_k=flibs.gen_ham(klist,ham_r,rvec)
        Gk_mats=flibs.gen_green(selfen,ham_k,mu_self,temp)
        Niw=min(40,len(iw_grid))  # Pade becomes unstable with too many Matsubara points
        # A(k,w) = -Im Tr G(k,w+i*delta)
        Gk=-flibs.pade_with_trace(Gk_mats[:,:,:Niw,:],1j*iw_grid[:Niw],wlist+1j*delta).imag
    else:
        Gk=flibs.gen_tr_Greenw_0(eig,mu,wlist,delta)
    w,x=np.meshgrid(wlist,spa_length)
    plt.contourf(x,w,Gk,cmap=plt.hot())
    for x in xticks[1:-1]:
        plt.axvline(x,ls='-',lw=0.25,color='black')
    plt.xlim(0,spa_length.max())
    plt.axhline(0.,ls='--',lw=0.25,color='black')
    plt.xticks(xticks,xlabel)
    plt.colorbar()
    plt.show()

def get_hall_coe(rvec,ham_r,S_r,avec,Nx:int,Ny:int,Nz:int,
                               fill:float,temp:float,tau_const,Nw=300,with_spin=False):
    # Parameter validation
    if temp <= 0:
        print("Error: Temperature (temp) is non-positive",flush=True)
        return
    if tau_const <= 0:
        print("Error: Relaxation time (tau_const) is non-positive",flush=True)
        return

    Nk,eig,vk,imass,kweight=plibs.get_emesh(Nx,Ny,Nz,ham_r,S_r,rvec,avec.T*ihbar,sw_veloc=True,sw_mass=True)
    Vuc=sclin.det(avec)*1e-30
    if Vuc <= 0:
        print("Error: Unit cell volume (Vuc) is non-positive",flush=True)
        return
    gsp=(1.0 if with_spin else 2.0) # spin degeneracy factor
    mu=plibs.calc_mu(eig,Nk,fill,temp)
    tau_mode=0
    if tau_mode==0:
        tau=eig*0.+tau_const  # constant relaxation time (k- and band-independent)
    else:
        # Energy-dependent tau from DOS (tau_mode != 0 path; not fully implemented)
        Nk,klist,eig,uni,kweight=plibs.get_emesh(Nx,Ny,Nz,ham_r,S_r,rvec,avec,sw_uni=True)
        wlist=np.linspace(eig.min()-mu,eig.max()-mu,Nw,True)
        Dos=flibs.gen_dos(eig,uni,mu,wlist,delta)
        tau=flibs.get_tau(Dos.sum(axis=0),eig,tau_const,tau_mode)
    print(f"T = {temp/kb:.3f} K",flush=True)
    print(f"mu = {mu:.4f} eV",flush=True)
    if tau_mode==0:
        print(f"tau = {tau_const} "+('fs' if sw_unit else ''),flush=True)
    else:
        print(f"max tau = {tau.max()} "+('fs' if sw_unit else ''),flush=True)
    sigma_hall=flibs.calc_sigmahall(eig,vk,imass/eC,kweight,tau,temp,mu)
    # K0 = charge transport kernel, K1 = thermoelectric kernel, K2 = thermal transport kernel
    K0,K1,K2=flibs.calc_Kn(eig,vk,kweight,temp,mu,tau)
    print(f"sigma_hall={sigma_hall:.6e}, K0[0,0]={K0[0,0]:.6e}, K0[1,1]={K0[1,1]:.6e}")
    # Check if K0 diagonal elements are non-zero (use tolerance for floating-point comparison)
    tol=1e-14
    if abs(K0[0,0]) < tol or abs(K0[1,1]) < tol:
        print("Error: K0 diagonal elements are too small. Cannot compute Hall coefficient",flush=True)
        return
    # Rh = -1/(n*e) from sigma_xy / (sigma_xx * sigma_yy); nh in cm^-3
    Rh=-Vuc*Nk*sigma_hall/(gsp*K0[0,0]*K0[1,1])
    nh=-1./(Rh*eC)/1e6
    print(f"Hall coefficient Rh = {Rh:.6e}",flush=True)
    print(f"Hole carrier density nh = {nh:.6e}",flush=True)

def calc_conductivity_Boltzmann(rvec,ham_r,S_r,avec,Nx:int,Ny:int,Nz:int,
                               fill:float,temp:float,tau_const,Nw=300,with_spin=False):
    '''
    calculate conductivities using Boltzmann equations
    '''
    #no dep. T and mu or filling
    # Parameter validation
    if temp <= 0:
        print("Error: Temperature (temp) is non-positive",flush=True)
        return
    if tau_const <= 0:
        print("Error: Relaxation time (tau_const) is non-positive",flush=True)
        return

    Nk,eig,vk,kweight=plibs.get_emesh(Nx,Ny,Nz,ham_r,S_r,rvec,avec.T*ihbar,sw_veloc=True)
    Vuc=sclin.det(avec)*1e-30
    if Vuc <= 0:
        print("Error: Unit cell volume (Vuc) is non-positive",flush=True)
        return
    gsp=(1.0 if with_spin else 2.0) #spin weight
    iNV=1./(Nk*Vuc)
    itemp=1./temp
    mu=plibs.calc_mu(eig,Nk,fill,temp)
    tau_mode=0
    if tau_mode==0:
        tau=eig*0.+tau_const
    else:
        Nk,klist,eig,uni,kweight=plibs.get_emesh(Nx,Ny,Nz,ham_r,S_r,rvec,avec,sw_uni=True)
        wlist=np.linspace(eig.min()-mu,eig.max()-mu,Nw,True)
        Dos=flibs.gen_dos(eig,uni,mu,wlist,delta)
        tau=flibs.get_tau(Dos.sum(axis=0),eig,tau_const,tau_mode)
    print(f"T = {temp/kb:.3f} K",flush=True)
    print(f"mu = {mu:.4f} eV",flush=True)
    if tau_mode==0:
        print(f"tau = {tau_const} "+('fs' if sw_unit else ''),flush=True)
    else:
        print(f"max tau = {tau.max()} "+('fs' if sw_unit else ''),flush=True)
    if sw_tdf:
        tdf=flibs.calc_tdf(eig,vk,kweight,tau,Nw)
        # absolute-energy grid matching the calc_tdf bins: E_i = emin + i*dw (i=1..Nw)
        dw=(eig.max()-eig.min())/Nw
        wlist=eig.min()+dw*np.arange(1,Nw+1)
        dfermi=0.25*(1.-np.tanh(0.5*(wlist-mu)/temp)**2)/temp
        K0=(dfermi*tdf.T).T.sum(axis=0)*dw
        K1=(dfermi*(wlist-mu)*tdf.T).T.sum(axis=0)*dw
        K2=(dfermi*(wlist-mu)**2*tdf.T).T.sum(axis=0)*dw
        plt.plot(wlist-mu,tdf[:,0,0])
        plt.show()
    else:
        K0,K1,K2=flibs.calc_Kn(eig,vk,kweight,temp,mu,tau)
    # sigma [S/m] = gsp * tau[s] * e[C] * K0[eV·m^-3] / (N*V[m^3])
    sigma=gsp*tau_unit*eC*K0*iNV  # e cancels eV->J so units work out to S/m
    kappa=gsp*tau_unit*eC*kb*K2*iNV*itemp

    # Handle sclin.inv() failures
    try:
        K0_inv = sclin.inv(K0)
        kappa2=gsp*tau_unit*eC*kb*(K2-K1.dot(K0_inv.dot(K1)))*iNV*itemp
        Seebeck=-kb*K0_inv.dot(K1)*itemp
        Peltier=K1.dot(K0_inv)
    except np.linalg.LinAlgError:
        print("Error: K0 is a singular matrix. Cannot compute Seebeck coefficient and Peltier",flush=True)
        kappa2=kappa.copy()
        Seebeck=np.zeros_like(K0)
        Peltier=np.zeros_like(K0)

    try:
        sigma_inv = sclin.inv(sigma*temp)
        # Lorenz tensor L = kb * kappa . (sigma*T)^-1 (matrix product, not element-wise)
        Lorenz=kb*kappa.dot(sigma_inv)
        Lorenz2=kb*kappa2.dot(sigma_inv)
    except np.linalg.LinAlgError:
        print("Warning: sigma*temp is a singular matrix. Cannot compute Lorenz coefficient",flush=True)
        Lorenz=np.zeros_like(sigma)
        Lorenz2=np.zeros_like(sigma)
    sigmaS=gsp*tau_unit*kb*eC*K1*iNV*itemp
    PF=sigma*Seebeck**2
    print_matrix('sigma matrix (S/m)',sigma)
    print_matrix('kappa matrix (K22 only) (W/m/K)',kappa)
    print_matrix('kappa matrix (full) (W/m/K)',kappa2)
    print_matrix('sigmaS matrix (A/m/K)',sigmaS)
    print_matrix('Seebeck matrix (V/K)',Seebeck)
    print_matrix('Peltier matrix (V)',Peltier)
    print_matrix('Lorenz matrix (K22 only) (Wohm/K^2) (fe 2.44e-8)',Lorenz,width=9,prec=2)
    print_matrix('Lorenz matrix? (full) (Wohm/K^2)',Lorenz2,width=9,prec=2)
    print_matrix('Power Factor (SA/m^2/K)',PF,ndigits=8)

def calc_conductivity_lrt(rvec,ham_r,S_r,avec,Nx:int,Ny:int,Nz:int,fill:float,
                         temp:float,Nw:int,delta,with_spin=False):
    '''
    calculation of linear response theory
    electric conductivity of LRT corresponds to Boltzmann then delta~O(10-1) (tau~1fs) at 300K
    thermal conductivity and L12 of LRT correspond to Boltzmann then delta~O(10-3) (tau~100fs) at 300K
    '''
    # Parameter validation
    if temp <= 0:
        print("Error: Temperature (temp) is non-positive",flush=True)
        return
    if delta <= 0:
        print("Error: Broadening (delta) is non-positive",flush=True)
        return

    Nk,eig,vk,kweight=plibs.get_emesh(Nx,Ny,Nz,ham_r,S_r,rvec,avec.T*ihbar,True,True)
    Vuc=sclin.det(avec)
    if Vuc <= 0:
        print("Error: Unit cell volume (Vuc) is non-positive",flush=True)
        return
    gsp=(1.0 if with_spin else 2.0) #spin weight
    mu=plibs.calc_mu(eig,Nk,fill,temp)
    print(f'chemical potential = {mu:.4f} eV',flush=True)
    print(f'temperature = {temp/kb:.3f} K',flush=True)
    print(f'about tau = {1.e15*hbar/delta:.3f} fs',flush=True)
    print(f'delta = {delta:9.3e}',flush=True)
    L11,L12,L22,wlist=plibs.get_conductivity(mu,temp,eig,vk,Nw,Emax,delta)
    sigmaconst=gsp*hbar*eC/Vuc*1.0e30
    kappaSconst=sigmaconst*kb/temp
    sigma=sigmaconst*L11
    kappa=kappaSconst*L22
    sigmaS=kappaSconst*L12
    # full electronic thermal conductivity subtracts the thermoelectric backflow:
    # kappa = kappaSconst*(L22 - L12.L11^-1.L12)  (cf. Boltzmann K2-K1.K0^-1.K1).
    # The subtraction is done on the bare L tensors (all share the same prefactor),
    # so kappa2 carries the same kappaSconst as the L22-only kappa.
    try:
        kappa2=kappaSconst*np.array([l22-l12.dot(sclin.inv(l11).dot(l12))
                                     for l11,l12,l22 in zip(L11,L12,L22)])
    except np.linalg.LinAlgError:
        print("Warning: L11 is a singular matrix. Cannot subtract thermoelectric backflow from kappa",flush=True)
        kappa2=kappa.copy()
    # Handle sclin.inv() failures
    try:
        Seebeck=np.array([-sclin.inv(s).dot(sS) for s,sS in zip(sigma,sigmaS)])
    except np.linalg.LinAlgError:
        print("Warning: sigma is a singular matrix. Cannot compute Seebeck coefficient",flush=True)
        Seebeck=np.zeros_like(sigma)
    print_matrix('sigma matrix (S/m)',sigma[0].real)
    print_matrix('kappa matrix (L22 only) (W/m/K)',kappa[0].real)
    print_matrix('kappa matrix (full) (W/m/K)',kappa2[0].real)
    print_matrix('sigmaS matrix (A/m/K)',sigmaS[0].real)
    print('Lorenz number (Wohm/K^2) (fe 2.44e-8)',flush=True)
    try:
        # Lorenz tensor L = kb * kappa . (sigma*T)^-1 (matrix product, not element-wise)
        sigmaT_inv=sclin.inv(sigma[0]*temp)
        Lorenz_matrix=(kb*kappa[0].dot(sigmaT_inv)).real.round(10)
        Lorenz_matrix2=(kb*kappa2[0].dot(sigmaT_inv)).real.round(10)
        print(' (L22 only):',flush=True)
        for lor in Lorenz_matrix:
            print(f" {lor[0]:9.2e} {lor[1]:9.2e} {lor[2]:9.2e}",flush=True)
        print(' (full):',flush=True)
        for lor in Lorenz_matrix2:
            print(f" {lor[0]:9.2e} {lor[1]:9.2e} {lor[2]:9.2e}",flush=True)
    except np.linalg.LinAlgError:
        print("Warning: Failed to compute Lorenz coefficient (singular matrix)",flush=True)
    print_matrix('Seebeck coefficient matrix (V/K)',Seebeck[0].real)
    fig=plt.figure()
    ax=fig.add_subplot(211)
    ax.plot(wlist,sigma[:,0,0].real)
    ax.plot(wlist,sigma[:,0,0].imag)
    ax2=fig.add_subplot(212)
    pol=1.+4*hbar*np.pi*1j*(sigma.T/(wlist+1e-8)).T
    ax2.plot(wlist[1:],pol[1:,0,0].real)
    ax2.plot(wlist[1:],pol[1:,0,0].imag)
    #ax2=fig.add_subplot(312)
    #ax2.plot(wlist,sigmaw.imag)
    #ax2.plot(wlist,(kappa[:,0,0]+kappa[:,1,1]+kappa[:,2,2]).real)
    #ax2.plot(wlist,(kappa[:,0,0]+kappa[:,1,1]+kappa[:,2,2]).imag)
    #ax3=fig.add_subplot(313)
    #ax3.plot(wlist,(sigmaS[:,0,0]+sigmaS[:,1,1]+sigmaS[:,2,2]).real)
    #ax3.plot(wlist,(sigmaS[:,0,0]+sigmaS[:,1,1]+sigmaS[:,2,2]).imag)
    plt.show()


def output_Fk(Nx:int,Ny:int,Nz:int,Nw:int,ham_r,S_r,rvec,plist,mu:float,temp:float,sw_self:bool,
              sw_soc=False,invs=None,slist=None,gap_sym=0):
    klist,kmap,invk=flibs.gen_irr_k_TRS(Nx,Ny,Nz)
    eig,uni=plibs.get_eigs(klist,ham_r,S_r,rvec)
    if sw_self:
        ham_k=flibs.gen_ham(klist,ham_r,rvec)
        npz=np.load('self_en.npz')
        sigmak,mu_self=npz['arr_0'],npz['arr_1']
        print(f'chem. pot. with self= {mu:.4f} eV',flush=True)
        Gk=flibs.gen_green(sigmak,ham_k,mu_self,temp)
    else:
        Gk=flibs.gen_Green0(eig,uni,mu,temp,Nw)
    try:
        gap=np.load('gap.npy')
    except FileNotFoundError:
        print("Error: 'gap.npy' not found",flush=True)
        return
    if sw_soc:
        Fk=flibs.gen_Fk_soc(Gk,gap,invk,invs,slist, gap_sym)
        Norb=int(len(slist)/2)
        if gap_sym>=0:
            Fks=Fk[:Norb,Norb:,0,:]-Fk[Norb:,:Norb,0,:]
            Fktr=np.array([f.diagonal().sum() for f in Fks.T])
        else:
            Fkt=Fk[:Norb,Norb:,0,:]+Fk[Norb:,:Norb,0,:]
            Fktr=np.array([f.diagonal().sum() for f in Fkt.T])
    else:
        Fk0=flibs.gen_Fk(Gk,gap,invk)
        Fk=flibs.remap_gap(Fk0[:,:,0,:],plist,invk,gap_sym)
        Fktr=np.array([f.diagonal().sum() for f in Fk.T])
    print('output anomalous green function')
    try:
        with open(f'Fk_tr.dat','w') as f:
            iwlist=(2*np.arange(Nw)+1)*np.pi*temp
            #wlistg=np.linspace(0,10,10)
            #Fkw=flibs.pade_with_trace(Fk[:,:,:30,:],iwlist[:30]*1j,wlistg+3e-2*1j)
            for i,km in enumerate(kmap):
                if km[2]==0:
                    #f.write(f'{km[0]:3} {km[1]:3} {Fks[i,0].real:12.8f} {Fks[i,0].imag:12.8f}\n')
                    f.write(f'{km[0]:3} {km[1]:3} {Fktr[i].real:12.8f} {Fktr[i].imag:12.8f}\n')
                    if km[0]==Nx-1:
                        f.write('\n')
    except IOError as e:
        print(f"Error: Failed to write 'Fk_tr.dat': {e}",flush=True)
    emax=(eig.max()-(mu_self if sw_self else mu))*1.2
    emin=(eig.min()-(mu_self if sw_self else mu))*1.2
    wlist=np.linspace(emin,emax,500)
    try:
        with open(f'Gpade.dat','w') as fp:
            idelta=1e-2
            Gkw=flibs.pade_with_trace(Gk[:,:,:40,:],iwlist[:40]*1j,wlist+idelta*1j)
            for i,km in enumerate(klist):
                if km[1]==0.0 and km[2]==0.0:
                    for j,w in enumerate(wlist):
                        fp.write(f'{km[0]:3} {w.real:12.8f} {-Gkw[i,j].imag:12.8f}\n')
                    fp.write('\n')
    except IOError as e:
        print(f"Error: Failed to write 'Gpade.dat': {e}",flush=True)
    #maxgap=abs(gap).max()
    #print((abs(gap)<maxgap*1.0e-6).sum())
    #print((abs(gap)<maxgap*1.0e-6).sum()/gap.size*100)
    #plt.plot(iwlist,gap[2,2,:,i].real)
    #plt.plot(iwlist,gap[2,2,:,100].imag)
    #plt.show()
    if sw_self:
        maxsigma=abs(sigmak).max()
        #print((abs(sigmak)<maxsigma*1.0e-6).sum())
        #print((abs(sigmak)<maxsigma*1.0e-6).sum()/sigmak.size*100)
        plt.plot(iwlist,sigmak[0,4,:,318].real,color='r')
        #plt.plot(-iwlist,Gk[4,0,:,318].real,color='r')
        plt.plot(iwlist,sigmak[0,4,:,318].imag,color='b')
        #plt.plot(-iwlist,-Gk[4,0,:,318].imag,color='b')
        plt.show()
    info=plibs.output_gap_function(invk,kmap,gap,uni,plist,gap_sym,Nx,sw_soc,invs,slist)


def get_mass(mesh,rvec,ham_r,S_r,mu:float,de=3.e-6,meshkz=20):
    """
    calculate cyclotron mass m*_c=hbar^2/2pi (dS/dE) where S is the area of Fermi surface cross section
    de: energy step for numerical derivative (eV)
    meshkz: number of kz points for coarse scan of S(kz) in [0, pi/2]
    """
    # Parameter validation
    if de <= 0:
        print("Error: Energy step (de) is non-positive",flush=True)
        return
    if meshkz <= 0:
        print("Error: kz mesh size (meshkz) is non-positive",flush=True)
        return

    al=alatt[:2]
    if al[0] <= 0 or al[1] <= 0:
        print("Error: Lattice constant (al) is non-positive",flush=True)
        return
    ABZ=4.*np.pi**2/(al[0]*al[1])

    # --- Phase 1: coarse scan S(kz) ---
    print("Phase 1: scanning S(kz)...",flush=True)
    S_scan,eig_cache=plibs.scan_fs_area(mesh,rvec,ham_r,S_r,RotMat,mu,ABZ,meshkz)

    # --- Phase 2: refine extremal kz, compute m* ---
    print("Phase 2: computing m* at extremal orbits...",flush=True)
    print("="*50,flush=True)
    for band_idx in sorted(S_scan.keys()):
        data=np.array(S_scan[band_idx])
        kz_arr,S_arr=data[:,0],data[:,1]
        if len(kz_arr)<3:
            continue
        cand_kz=plibs.find_extremal_kz(kz_arr,S_arr,band_idx,mesh,rvec,ham_r,S_r,RotMat,mu,ABZ)
        print(f"Band {band_idx+1}: extremal kz = {[f'{k:.4f}' for k in cand_kz]}",flush=True)

        results=[]
        for kz_ext in cand_kz:
            eig=eig_cache.get(kz_ext)
            if eig is None:
                eig=plibs.get_eigs_2d(mesh,rvec,ham_r,S_r,RotMat,kz_ext)
            v2,blist=plibs.get_kf_points(eig,mesh,mu,kz_ext)
            S0=plibs.get_band_area(v2,blist,band_idx,ABZ)
            if not S0:
                continue
            v2_p,blist_p=plibs.get_kf_points(eig,mesh,mu+de,kz_ext)
            v2_m,blist_m=plibs.get_kf_points(eig,mesh,mu-de,kz_ext)
            Sp=plibs.get_band_area(v2_p,blist_p,band_idx,ABZ)
            Sm=plibs.get_band_area(v2_m,blist_m,band_idx,ABZ)
            if Sp is None or Sm is None:
                continue
            # dS/dE via central finite difference; *1e20 converts AA^-2/eV -> m^-2/J
            dSdE_SI=(Sp-Sm)/(2.*de)*1.e20
            # Onsager: m*_c = (hbar^2 / 2pi) * dS/dE  [in units of electron mass m_e]
            mc=np.abs(dSdE_SI)*eC*hbar**2/(2*np.pi*emass)
            results.append((kz_ext,S0,mc))
            print(f"  kz={kz_ext:.4f}: S={S0:.4f} AA^-2, m*={mc:.4f} m_e",flush=True)

        if results:
            mc_vals=np.array([r[2] for r in results])
            kz_vals=np.array([r[0] for r in results])
            i_max,i_min=np.argmax(mc_vals),np.argmin(mc_vals)
            print(f"  >> Max m* = {mc_vals[i_max]:.4f} m_e  at kz = {kz_vals[i_max]:.4f}",flush=True)
            print(f"  >> Min m* = {mc_vals[i_min]:.4f} m_e  at kz = {kz_vals[i_min]:.4f}",flush=True)

def get_dhva_band(mesh,rvec,ham_r,S_r,mu:float,theta_list,phi=0.,meshkz=20):
    """
    Calculate dHvA frequency F vs magnetic field polar angle theta.
    F = hbar/(2*pi*e) * A_ext  (Onsager relation)
    theta_list: array of polar angles from z-axis [deg]
    phi       : azimuthal angle [deg] (fixes the rotation plane, default 0 = xz-plane)
    meshkz    : number of kz scan points per angle for the coarse S(kz) scan
    Returns   : dict  band_idx -> np.ndarray of shape (N, 2) columns=[theta, F[T]]
    """
    import matplotlib.pyplot as plt
    al=alatt[:2]
    ABZ=4.*np.pi**2/(al[0]*al[1])
    F_factor=hbar/(2.*np.pi)*1.e20   # AA^-2 -> T (Onsager): hbar[eV*s]/(2pi) * 1e20[AA^-2->m^-2]

    all_results={}  # band_idx -> [(theta, F), ...]

    for theta in theta_list:
        rotmat=plibs.make_rotmat(theta,phi)
        S_scan,eig_cache=plibs.scan_fs_area(mesh,rvec,ham_r,S_r,rotmat,mu,ABZ,meshkz)

        for band_idx in S_scan:
            data=np.array(S_scan[band_idx])
            kz_arr,S_arr=data[:,0],data[:,1]
            if len(kz_arr)<3:
                continue
            cand_kz=plibs.find_extremal_kz(kz_arr,S_arr,band_idx,mesh,rvec,ham_r,S_r,rotmat,mu,ABZ)

            for kz_ext in cand_kz:
                eig=eig_cache.get(kz_ext)
                if eig is None:
                    eig=plibs.get_eigs_2d(mesh,rvec,ham_r,S_r,rotmat,kz_ext)
                v2,blist_ext=plibs.get_kf_points(eig,mesh,mu,kz_ext)
                S0=plibs.get_band_area(v2,blist_ext,band_idx,ABZ)
                if not S0:
                    continue
                all_results.setdefault(band_idx,[]).append((theta,S0*F_factor))

        print(f"theta={theta:.1f} deg: FS bands={[b+1 for b in S_scan.keys()]}",flush=True)

    # convert to arrays
    all_results={k:np.array(v) for k,v in all_results.items()}

    # Group extremal F values by rank (ascending) at each theta and connect same-rank points
    # across angles to form continuous orbit branches (e.g., belly vs neck orbits)
    fig,ax=plt.subplots()
    prop_cycle=plt.rcParams['axes.prop_cycle'].by_key()['color']
    for ci,band_idx in enumerate(sorted(all_results.keys())):
        d=all_results[band_idx]
        color=prop_cycle[ci % len(prop_cycle)]
        thetas=np.unique(d[:,0])
        theta_F={th:sorted(d[d[:,0]==th,1]) for th in thetas}
        n_branches=max(len(v) for v in theta_F.values())
        for i in range(n_branches):
            pts=np.array([(th,flist[i])
                          for th,flist in sorted(theta_F.items()) if i<len(flist)])
            label=f'Band {band_idx+1}' if i==0 else '_nolegend_'
            ax.plot(pts[:,0],pts[:,1],'-o',color=color,markersize=3,label=label)
    ax.set_xlabel('theta (deg)')
    ax.set_ylabel('F (T)')
    ax.set_title(f'dHvA frequency vs field angle  phi={phi:.1f} deg')
    ax.legend()
    plt.tight_layout()
    plt.savefig('dhva_band.png',dpi=150)
    plt.show()
    return all_results

def calc_cpa_Akw(k_sets: list, kmesh: int, bvec: np.ndarray,
                 ham_r: np.ndarray, rvec: np.ndarray, mu: float,
                 sigma_cpa: np.ndarray, wlist: np.ndarray,
                 zlist_real: np.ndarray, xlabel: list) -> np.ndarray:
    """
    Compute and plot the CPA spectral function A(k,w) along a symmetry-line k-path.
    Returns Akw array of shape (Nks, Nw).
    """
    Norb = ham_r.shape[0]
    Nw = len(wlist)
    klist_s, spa_length, xticks = plibs.mk_klist(k_sets, kmesh, bvec)
    hamk_s = flibs.gen_ham(klist_s, ham_r, rvec)
    hamk_s[:, range(Norb), range(Norb)] -= mu
    Nks = len(klist_s)
    Akw = np.zeros((Nks, Nw))
    A_batch = (np.eye(Norb)[None, None, :, :] * zlist_real[:, None, None, None]
               - hamk_s[None, :, :, :]
               - sigma_cpa[:, None, :, :])  # (Nw, Nks, Norb, Norb)
    try:
        Gk_all = np.linalg.inv(A_batch)
        Akw = (-np.trace(Gk_all, axis1=2, axis2=3).imag / np.pi).T
    except np.linalg.LinAlgError:
        print("Warning: Batch matrix inversion failed, falling back to element-wise.", flush=True)
        for iw in range(Nw):
            zI = np.eye(Norb) * zlist_real[iw]
            sig = sigma_cpa[iw]
            for ik in range(Nks):
                try:
                    Gk = np.linalg.inv(zI - hamk_s[ik] - sig)
                    Akw[ik, iw] = -np.trace(Gk).imag / np.pi
                except np.linalg.LinAlgError:
                    print(f"Warning: Matrix inversion failed at k-point {ik}, frequency index {iw}. Setting to zero.", flush=True)
                    Akw[ik, iw] = 0.0
    w, x = np.meshgrid(wlist, spa_length)
    plt.contourf(x, w, Akw, 100, cmap=plt.get_cmap('hot'))
    for xt in xticks[1:-1]:
        plt.axvline(xt, ls='-', lw=0.25, color='white')
    plt.xlim(0, spa_length.max())
    plt.axhline(0., ls='--', lw=0.25, color='white')
    plt.xticks(xticks, xlabel)
    plt.ylabel('Energy (eV)')
    plt.title('CPA spectrum')
    plt.colorbar(label=r'$A(\mathbf{k},\omega)$')
    plt.tight_layout()
    plt.savefig('cpa_spectrum.png', dpi=150)
    plt.show()
    return Akw

def main():
    omp_num,omp_check=flibs.omp_params()

    # ===== Input file and parameter validation =====
    # Check if Hamiltonian file exists
    if not os.path.exists(fname):
        print(f"Error: Hamiltonian file '{fname}' not found",flush=True)
        return
    # Validate ftype (any integer outside {0,1,2,3} uses the Hopping.dat branch)
    if not isinstance(ftype, int):
        print(f"Error: ftype={ftype!r} must be an integer",flush=True)
        return
    # Validate brav (any integer outside {0..7} uses the monoclinic branch)
    if not isinstance(brav, int):
        print(f"Error: brav={brav!r} must be an integer",flush=True)
        return

    # ===== Computational parameter validation =====
    # Mesh size check
    if Nx <= 0 or Ny <= 0 or Nz <= 0 or Nw <= 0:
        print(f"Error: Invalid mesh size (Nx={Nx}, Ny={Ny}, Nz={Nz}, Nw={Nw})",flush=True)
        return
    # Temperature check
    # Many kernels use 1/temp or finite-T Fermi functions; T=0 is not supported in this implementation.
    if tempK <= 0:
        print(f"Error: Temperature (tempK) is non-positive ({tempK} K)",flush=True)
        return
    # U, J check
    if U < 0 or J < 0:
        print(f"Error: U, J are negative (U={U}, J={J})",flush=True)
        return
    # delta check
    if delta <= 0:
        print(f"Error: Broadening (delta) is non-positive ({delta})",flush=True)
        return
    # fill check
    if fill < 0:
        print(f"Error: Filling (fill) is negative ({fill})",flush=True)
        return
    # Emin, Emax check
    if Emin >= Emax:
        print(f"Error: Emin >= Emax ({Emin} >= {Emax})",flush=True)
        return
    # tau_const check
    if tau_const <= 0:
        print(f"Error: Relaxation time (tau_const) is non-positive ({tau_const})",flush=True)
        return
    # =========================================

    #import hamiltonian
    if ftype==3:
        rvec,ham_r,S_r,Norb,Nr=plibs.import_MLO_hoppings(fname)
    else:
        rvec,ham_r,Norb,Nr=plibs.import_hoppings(fname,ftype)
        S_r=[]
    plist=flibs.get_plist(rvec,ham_r)
    print("Effective parity of wannier functions:", plist, flush=True)

    # ===== Parameter checks dependent on number of orbitals =====
    print(f"Number of orbital = {Norb}",flush=True)

    # Validate olist
    for i, ol in enumerate(olist):
        if isinstance(ol, int):
            if ol < 0 or ol >= Norb:
                print(f"Error: olist[{i}]={ol} is invalid. Valid range: 0-{Norb-1}",flush=True)
                return
        elif isinstance(ol, list):
            for j, olj in enumerate(ol):
                if olj < 0 or olj >= Norb:
                    print(f"Error: olist[{i}][{j}]={olj} is invalid. Valid range: 0-{Norb-1}",flush=True)
                    return

    # Check fill does not exceed number of bands
    if fill > Norb:
        print(f"Error: filling={fill} exceeds number of bands={Norb}. Valid range: 0 < fill <= {Norb}",flush=True)
        return

    if orb_dep:
        # Orbital-dependent interaction mode requires full Norb x Norb U/J matrices.
        if 'Umat' not in globals() or 'Jmat' not in globals():
            print("Error: orb_dep=True requires Umat and Jmat to be defined",flush=True)
            return
        Umat_arr=np.asarray(Umat)
        Jmat_arr=np.asarray(Jmat)
        exp_shape=(Norb,Norb)
        if Umat_arr.shape != exp_shape or Jmat_arr.shape != exp_shape:
            print(f"Error: Umat/Jmat shape must be {exp_shape}, got {Umat_arr.shape}/{Jmat_arr.shape}",flush=True)
            return
        if np.iscomplexobj(Umat_arr) or np.iscomplexobj(Jmat_arr):
            print("Error: Umat/Jmat must be real-valued matrices",flush=True)
            return
        # Use contiguous float64 arrays for stable ctypes calls to Fortran wrappers.
        globals()['Umat']=np.ascontiguousarray(Umat_arr,dtype=np.float64)
        globals()['Jmat']=np.ascontiguousarray(Jmat_arr,dtype=np.float64)
    # =============================================

    #set lattice vector
    avec,Arot=plibs.get_ptv(alatt,deg,brav)
    #rotation axis
    if sw_dec_axis:
        rvec1=Arot.T.dot(rvec.T).T
        rvec=rvec1.copy()
        if brav in {1,2}:
            rvec[:,2]*=2.
            avec=alatt*np.eye(3)
            avec[:,2]*=.5
        elif brav==5:
            rvec*=2.
            avec=(alatt*np.eye(3))*.5
        else:
            avec=alatt*np.eye(3)
    bvec=plibs.get_bvec(avec) #set recp. lattice
    if omp_check: #OMP properties
        print("OpenMP mode",flush=True)
        print(f"Number of OpenMP threads = {omp_num}",flush=True)
    print(f"calc mode {int(option)}: "+option.description,flush=True)

    # ===== Additional input validation (lower priority) =====
    # J > U warning
    if J > U:
        print(f"Warning: J={J} > U={U} is physically unusual. Please verify",flush=True)

    # Validate gap_sym during eliashberg/flex calculations
    if option in MODES_FLEX_ELIASH:  # eliashberg/flex calculations
        if gap_sym not in {-1, 0, 1, 2, 3}:
            print(f"Warning: gap_sym={gap_sym} is non-standard. Common values: -1,0,1,2,3",flush=True)

    # Validate at_point for chis at q-point calculation (normal and SC state)
    if option in MODES_CHIS_QPOINT:  # chis at q-point
        if len(at_point) != 3:
            print(f"Error: at_point has {len(at_point)} elements. Required: 3 elements",flush=True)
            return
        if any(p < 0 or p > 1 for p in at_point):
            print(f"Warning: at_point={at_point} values outside [0,1] range. May not be normalized by reciprocal lattice",flush=True)
    # ==========================================
    if option in MODES_COLOR:
        print("color mode: "+color_option.description,flush=True)
    print("Hamiltonian name is "+fname,flush=True)
    print(f"Number of orbital = {Norb}",flush=True)
    if (orb_dep==False) and option in MODES_PRINT_UJ: #write constant U,J
        print(f'U= {U:5.2f} and J= {J:5.3f}')
    if option in MODES_NEED_CHI:
        """ chiolist is the list of orbital properties of index on chi """
        try:
            chiolist
        except NameError:
            try:
                site_prof
            except NameError:
                site_prof=[1] #one site (len(site_prof)=1)
            # Build the orbital-pair basis once here so every response/Eliashberg branch
            # shares the same indexing convention when passing chi objects to Fortran.
            chiolist,site=plibs.get_chi_orb_list(len(ham_r[0]),site_prof)
        if option in MODES_COULOMB_VERTEX:
            print("generate coulomb vertex matrix S")
            # Susceptibility branches only need the static interaction vertex; FLEX/Eliashberg
            # rebuild their own matrices inside the dedicated solvers.
            if orb_dep:
                Smat,Cmat=flibs.gen_SCmatrix_orb(chiolist,site,Umat,Jmat)
            else:
                Smat,Cmat=flibs.gen_SCmatrix(chiolist,site,U,J)
    if option in MODES_KMESH_SYMLINE:
        if sw_gen_sym:
            print('generate symmetry line',flush=True)
        print(f'kmesh = {kmesh}',flush=True)
    elif option in MODES_KMESH_SINGLE:
        print(f'Number of k-mesh = {Nx}',flush=True)
    else:
        print(f'k-mesh is {Nx} {Ny} {Nz}',flush=True)
    if option in MODES_MATSUBARA:
        print(f'Number of Matsubara freq. = {Nw}',flush=True)
    print("Lattice Vector (Angstrom)",flush=True)
    for i,a in enumerate(avec):
        print(f"a{i+1}: {a[0]:7.4f} {a[1]:7.4f} {a[2]:7.4f}",flush=True)
    print("Reciprocal Lattice Vector (Angstrom^-1)",flush=True)
    for i,b in enumerate(bvec):
        print(f"b{i+1}: %7.4f %7.4f %7.4f"%tuple(b),flush=True)
    if option in MODES_SELF_MU: #conductivity and impurity/CPA functions calc or set mu themself
        pass
    else:
        if sw_calc_mu:
            mu=plibs.get_mu(ham_r,S_r,rvec,Arot,temp,fill)
        else:
            mu=mu0
            print('use fixed mu')
        print(f'Temperature = {temp:10.3e} eV ({temp/kb:.2f} K)',flush=True)
        print(f'chem. pot. = {mu:.4f} eV',flush=True)
    if option==CalcMode.BAND: #plot band
        klist,spa_length,xticks=plibs.mk_klist(k_sets,kmesh,bvec)
        eig,uni0=plibs.get_eigs(klist,ham_r,S_r,rvec)
        uni=np.array([u.T for u in uni0]) #rotate uni(k,band,orb) to uni(k,orb,band)
        plot_band(eig.T-mu,spa_length,xlabel,xticks,uni.T,olist,(False if color_option==0 else True))
    elif option==CalcMode.DOS: #plot dos
        Nk,klist,eig,uni,kweight=plibs.get_emesh(Nx,Ny,Nz,ham_r,S_r,rvec,avec,sw_uni=True)
        wlist=np.linspace(Emin,Emax,Nw,True)
        Dos=flibs.gen_dos(eig,uni,mu,wlist,delta)
        if color_option==1:
            if len(olist)!=0:
                for ol,cl in zip(olist,['red','green','blue']):
                    if ol==[]:
                        pass
                    else:
                        if isinstance(ol,int):
                            plt.plot(wlist,Dos[ol],color=cl)
                        else:
                            plt.plot(wlist,Dos[ol,:].sum(axis=0),color=cl)
            plt.plot(wlist,Dos.sum(axis=0),color='black')
        elif color_option==0:
            clist=['k','r','g','b','c','m','y']
            tmp_ol_b=flatten_orbs(olist)
            plt.fill_between(wlist,Dos[tmp_ol_b].sum(axis=0),Dos.sum(axis=0),color=clist[0])
            for i in range(len(olist)-1):
                tmp_ol_u=flatten_orbs(olist[i:])
                tmp_ol_b=flatten_orbs(olist[i+1:])
                plt.fill_between(wlist,Dos[tmp_ol_b].sum(axis=0),Dos[tmp_ol_u].sum(axis=0),color=clist[i+1])
            plt.fill_between(wlist,0,Dos[tmp_ol_b].sum(axis=0),color=clist[len(olist)])
        plt.xlim(Emin,Emax)
        plt.ylim(0,max(Dos.sum(axis=0))*1.2)
        plt.show()
    elif option==CalcMode.FERMI_2D: #2D Fermi surface plot
        eig2d=plibs.get_eigs_2d(Nx,rvec,ham_r,S_r,RotMat,kz)
        klist,blist=plibs.get_kf_points(eig2d,Nx,mu,kz)
        clist=plibs.get_colors(klist,blist,ihbar*avec.T,rvec,ham_r,S_r,olist,color_option,True)
        plot_FS(clist,klist,color_option)
    elif option==CalcMode.FERMI_3D: #3D Fermi surface plot
        polys,centers,blist=plibs.gen_3d_surf_points(Nx,rvec,ham_r,S_r,mu,kscale)
        fspolys,fscenters,fscolors=set_init_3dfsplot(color_option,polys,centers,blist,avec,rvec,ham_r,S_r,olist)
        plot_3d_surf(fspolys,fscenters,fscolors,color_option,kscale,bvec)
    elif option==CalcMode.SPECTRUM: #plot spectrum
        plot_spectrum(k_sets,xlabel,kmesh,bvec,mu,ham_r,S_r,rvec,Emin,Emax,delta,Nw,sw_self)
    elif option==CalcMode.CONDUCTIVITY_BT: #calc conductivity
        get_hall_coe(rvec,ham_r,S_r,avec,Nx,Ny,Nz,fill,temp,tau_const)
        calc_conductivity_Boltzmann(rvec,ham_r,S_r,avec,Nx,Ny,Nz,fill,temp,tau_const)
    elif option==CalcMode.CONDUCTIVITY_PT: #calc_optical conductivity
        calc_conductivity_lrt(rvec,ham_r,S_r,avec,Nx,Ny,Nz,fill,temp,Nw,delta)
    elif option in MODES_CHI_NORMAL: #calc_chis_spectrum
        print("calculate electron energy",flush=True)
        Nk,klist,eig,uni,kweight=plibs.get_emesh(Nx,Ny,Nz,ham_r,S_r,rvec,avec,sw_uni=True)
        if option in {CalcMode.CHIS_SPECTRUM,CalcMode.PHI_SPECTRUM}: #chis/phi spectrum with symmetry line
            print("generate qlist",flush=True)
            qlist,spa_length,xticks=plibs.mk_qlist(k_sets,Nx,Ny,Nz,bvec)
            if option==CalcMode.CHIS_SPECTRUM:
                w,sp,sus=plibs.calc_path_spectrum('chis',mu,temp,klist,qlist,chiolist,eig,uni,spa_length,Nw,Emax,delta,Smat)
                print("write chis spectrum in png file",flush=True)
                susfname='chis_spec.png'
            elif option==CalcMode.PHI_SPECTRUM:
                w,sp,sus=plibs.calc_path_spectrum('phi',mu,temp,klist,qlist,chiolist,eig,uni,spa_length,Nw,Emax,delta)
                print("write phi spectrum in png file",flush=True)
                susfname='phi_spec.png'
            plt.contourf(sp,w,abs(sus.imag),100)
            plt.colorbar()
            plt.hot()
            plt.savefig(fname=susfname,dpi=300)
        elif option==CalcMode.CHIS_QPOINT: #chis at q-point
            q_point=np.array(at_point)
            chis,chis_orb,wlist=plibs.chis_q_point(q_point,eig,uni,Emax,Nw,mu,temp,Smat,klist,chiolist,delta)
            plt.plot(wlist,chis.imag)
            plt.show()
            for cso in chis_orb.T:
                plt.plot(wlist,cso.imag)
            plt.show()
        else: #chis/phi qmap at ecut plane
            if option==CalcMode.CHIS_QMAP: #chis spectrum ecut plane
                sus,chi0,qx,qy=plibs.chis_qmap(Nx,Ny,Ecut,mu,temp,Smat,klist,chiolist,eig,uni,idelta=1.e-3)
                plt.contourf(qx,qy,abs(chi0.imag),100)
                plt.colorbar()
                plt.jet()
                plt.savefig(fname='chi0map.png',dpi=300)
                susfname='chismap.png'
            elif option==CalcMode.PHI_QMAP:
                sus,qx,qy=plibs.phi_qmap(Nx,Ny,Ecut,mu,temp,klist,chiolist,eig,uni,idelta=1.e-3,sw_omega=sw_omega)
                susfname='phimap.png'
            if sw_omega:
                plt.contourf(qx,qy,abs(sus.imag),100)
            else:
                plt.contourf(qx,qy,abs(sus.real),100)
            plt.colorbar()
            plt.jet()
            #plt.show()
            plt.savefig(fname=susfname,dpi=300)
    elif option in MODES_CHIS_SC: #calc chis at superconducting state
        sw_spsym=True if gap_sym<0 else False
        if isinstance(delta0, float):
            klist,kmap,invk=flibs.gen_irr_k_TRS(Nx,Ny,Nz)
            eig,uni=plibs.get_eigs(klist,ham_r,S_r,rvec)
            deltaini=flibs.remap_gap(flibs.get_initial_delta(plibs.get_initial_gap(klist,len(eig.T),gap_sym),
                                                              uni, kmap, invk, 1, gap_sym),plist,invk,gap_sym)
            deltaini_static=np.ascontiguousarray(deltaini[:,:,0,:].transpose(2,0,1))  # [Nkall,Norb,Norb]
            print("maximum gap = %7.4f meV"%(delta0*1e3),flush=True)
            gap_max=abs(deltaini_static).max()
            if gap_max<=0.0:
                print("Warning: initial gap is zero. Using unscaled delta",flush=True)
                deltak=deltaini_static.copy()
            else:
                deltak=delta0*deltaini_static/gap_max
            # Full-grid klist in kmap (FFT) order: keeps deltak <-> klist index pairing
            # consistent, and gives the [0,1) endpoint-free mesh required by get_qshift
            klist=kmap/np.array([Nx,Ny,Nz],dtype=np.float64)
            Nk=len(klist)
        elif(isinstance(delta0, np.ndarray) or isinstance(delta0, list)):
            # sw_pp=False: get_qshift requires the periodic [0,1) mesh without endpoints
            Nk, klist = plibs.gen_klist(Nx, Ny, Nz, sw_pp=False)
            eig,uni=plibs.get_eigs(klist,ham_r,S_r,rvec)
            inigap=plibs.gap_symms(klist,Norb,gap_sym)
            for i,k in enumerate(klist):
                # AFM-zone test in centered fractional coords: keep the full gap inside the
                # M-centered diamond (|kx|+|ky|>=1/2) and damp bands 1,2 inside the
                # Gamma-centered diamond (|kx|+|ky|<1/2).
                kc=k-np.round(k)
                if(abs(kc[0])+abs(kc[1])>=0.5):
                    pass
                else:
                    inigap[1,i]=inigap[1,i]*.5
                    inigap[2,i]=inigap[2,i]*.5
            # delta0 is per-band [Norb]; broadcast over k: (inigap.T)[Nk,Norb] * delta0[Norb] → [Norb,Nk]
            deltaini=(inigap.T*np.array(delta0,dtype=np.complex128)).T
            deltak=flibs.get_band_to_orb_delta(deltaini,uni)
        hamk = flibs.gen_ham(klist, ham_r, rvec)
        if option==CalcMode.CHIS_SPECTRUM_SC:
            print("generate qlist",flush=True)
            qlist,spa_length,xticks=plibs.mk_qlist(k_sets,Nx,Ny,Nz,bvec)
            w,sp,sus=plibs.chis_spectrum_sc(mu, temp, Smat, hamk, deltak, klist, qlist, chiolist, Nw, Emax, delta, sw_spsym)
            print("write chis spectrum in png file",flush=True)
            plt.contourf(sp,w,abs(sus.imag),100)
            plt.colorbar()
            plt.hot()
            plt.savefig(fname='chis_sc_spec.png',dpi=300)
            plt.close()
        elif option==CalcMode.CHIS_QPOINT_SC:
            # --- BdG band structure along (0,0,0) -> (0,0.5,0) ---
            Nk_path = 200
            kpath_bdg = np.zeros((Nk_path, 3))
            kpath_bdg[:, 1] = np.linspace(0, 0.5, Nk_path)
            hamk_path = flibs.gen_ham(kpath_bdg, ham_r, rvec)
            _, uni_path = plibs.get_eigs(kpath_bdg, ham_r, S_r, rvec)
            inigap_path = plibs.gap_symms(kpath_bdg, Norb, gap_sym)
            if isinstance(delta0, float):
                scale = delta0 / gap_max if gap_max > 0 else delta0
                deltak_path = flibs.get_band_to_orb_delta(
                    (inigap_path * scale).astype(np.complex128), uni_path)
            else:
                deltak_path = flibs.get_band_to_orb_delta(
                    (inigap_path.T * np.array(delta0, dtype=np.complex128)).T, uni_path)
            eig_BdG_path, _ = flibs.get_eig(
                flibs.mkBdGhamk(hamk_path - mu * np.eye(Norb), deltak_path))
            ky = kpath_bdg[:, 1]
            for n in range(eig_BdG_path.shape[1]):
                plt.plot(ky, eig_BdG_path[:, n], 'b-', lw=0.5)
            plt.axhline(0, ls='--', lw=0.5, color='gray')
            plt.xlabel(r'$k_y \; (2\pi/a)$')
            plt.ylabel('Energy (eV)')
            plt.title('BdG bands  (0,0,0) → (0,0.5,0)')
            plt.tight_layout()
            plt.savefig('BdG_band.png', dpi=150)
            #plt.show()
            plt.clf()
            # --- chi calculation ---
            q_point=np.array(at_point)
            print(f"q-point is ({q_point[0]:.3f}, {q_point[1]:.3f}, {q_point[2]:.3f})",flush=True)
            chis,chis_orb,wlist=plibs.chis_q_point_sc(q_point, hamk, deltak, mu, Emax, Nw, temp,
                                                      Smat, klist, chiolist, delta, sw_spsym)
            plt.plot(wlist,chis.imag)
            plt.savefig('chisq.png', dpi=150)
            plt.clf()
            #plt.show()
            with open(f"chis_sc.dat","w") as f:
                for iw,ic in zip(wlist,chis):
                    f.write(f"{iw:5.3f}, {ic.imag:12.8f}, {ic.real:12.8f}\n")
            for cso in chis_orb.T:
                plt.plot(wlist,cso.imag)
            #plt.show()
            plt.savefig('chisq_orb.png', dpi=150)
            plt.close()
            with open(f"chis_scorb.dat",'w') as f:
                for iw,ic in zip(wlist,chis_orb):
                    f.write(f"{iw:5.3f}, ")
                    for cso in ic:
                        f.write(f"{cso.imag:12.8f}, ")
                    f.write("\n")
    elif option in MODES_FLEX_ELIASH: #flex/eliashberg calculations
        if sw_soc: #with soc
            try:
                slist
            except NameError: #default, up(~norb/2-1)>down(norb/2~)
                Norb=len(ham_r[0])
                # Correctly split even if Norb is odd (up-spin majority)
                slist=np.ones(Norb,dtype=np.int64)
                slist[(Norb+1)//2:]=-1  # Correct split: for odd case [+1,+1,+1,-1,-1] (Norb=5)
            try:
                invs
            except NameError:
                # Also fixed reverse index of split spins
                invs=np.concatenate([np.arange((Norb+1)//2,Norb),np.arange((Norb+1)//2)])+1
            if option==CalcMode.FLEX: #calc self-energy using flex
                plibs.calc_flex_soc(Nx,Ny,Nz,Nw,ham_r,S_r,rvec,mu,temp,chiolist,slist,invs,site,
                                    orb_dep,U,J,fill,sw_out_self,sw_in_self,
                                    Umat if orb_dep else None,Jmat if orb_dep else None,
                                    m_diis=m_diis_num,sw_rescale=sw_rescale_flex)
            elif option==CalcMode.LIN_ELIASHBERG: #calc gap function
                plibs.calc_lin_eliash_soc(Nx,Ny,Nz,Nw,ham_r,S_r,rvec,mu,temp,chiolist,slist,plist,invs,site,
                                          orb_dep,U,J,fill,gap_sym,sw_self,sw_from_file,sw_out_self,sw_in_self,
                                          Umat if orb_dep else None,Jmat if orb_dep else None)
            elif option==CalcMode.GAP_FUNCTION: #post gap calculation, output gap function/anomalous green's function
                output_Fk(Nx,Ny,Nz,Nw,ham_r,S_r,rvec,plist,mu,temp,sw_self,sw_soc,invs,slist,gap_sym)
        else: #without soc
            if option==CalcMode.FLEX: #calc self-energy using flex
                plibs.calc_flex(Nx,Ny,Nz,Nw,ham_r,S_r,rvec,mu,temp,chiolist,site,
                                orb_dep,U,J,fill,sw_out_self,sw_in_self,
                                Umat if orb_dep else None,Jmat if orb_dep else None,
                                m_diis=m_diis_num,sw_rescale=sw_rescale_flex)
            elif option==CalcMode.LIN_ELIASHBERG: #calc gap function
                plibs.calc_lin_eliashberg_eq(Nx,Ny,Nz,Nw,ham_r,S_r,rvec,chiolist,site,plist,mu,temp,gap_sym,sw_self,
                                             orb_dep,U,J,fill,sw_from_file,sw_out_self,sw_in_self,
                                             Umat if orb_dep else None,Jmat if orb_dep else None)
            elif option==CalcMode.GAP_FUNCTION: #post gap calculation, output gap function/anomalous green's function
                output_Fk(Nx,Ny,Nz,Nw,ham_r,S_r,rvec,plist,mu,temp,sw_self)
    elif option==CalcMode.EILENBERGER: #solve homogeneous quasiclassical Eilenberger equation
        if eil_fs: #model FS + Fermi velocity: anisotropic penetration depth lambda_xx/lambda_yy
            plibs.calc_fs_penetration(eil_coupling,temp,eil_wc,kind=(eil_fs_kind or 'ellipse'),
                                      gap_sym=gap_sym,params=eil_fs_params,kb=kb)
        elif eil_lambda: #superfluid density / penetration depth lambda(T) (s exp-flat, d linear-in-T)
            plibs.calc_penetration_depth(eil_coupling,temp,eil_wc,gap_sym=gap_sym,kb=kb)
        elif eil_spin: #spin-2x2 Zeeman response: singlet vs triplet d-vector (d||h vs d_|_h)
            plibs.calc_spin_pauli(Nx,Ny,Nz,eil_wc,ham_r,S_r,rvec,avec,mu,temp,eil_coupling,
                                  gap_sym=gap_sym,fs_width=eil_fs_width,kb=kb)
        elif eil_pauli: #Zeeman (Maki) Pauli-limiting sweep: singlet gap Delta(h), spinodal, Zeeman-split DOS
            plibs.calc_pauli_limit(Nx,Ny,Nz,eil_wc,ham_r,S_r,rvec,avec,mu,temp,gap_sym,eil_coupling,
                                   fs_width=eil_fs_width,kb=kb)
        elif eil_free_energy: #condensation free energy (Omega_s-Omega_n)/N0 vs T (coupling-independent)
            plibs.calc_free_energy(eil_coupling,temp,eil_wc,gap_sym=gap_sym,kb=kb)
        else:
            plibs.calc_eilenberger(Nx,Ny,Nz,eil_wc,ham_r,S_r,rvec,avec,mu,temp,gap_sym,eil_coupling,
                                   imp_gamma=eil_imp_gamma,imp_c=eil_imp_c,fs_width=eil_fs_width,kb=kb,
                                   method=eil_method,sw_find_tc=eil_find_tc,sw_imp_sweep=eil_imp_sweep,
                                   imp_sweep=eil_imp_list)
    elif option==CalcMode.EILENBERGER_SURFACE: #specular surface state via Riccati Eilenberger (model FS)
        if eil_surf_dvector: #self-consistent triplet d-vector texture (spin-matrix Riccati)
            plibs.calc_surface_dvector(eil_coupling,temp,eil_wc,kb=kb,sub_ratio=eil_dvec_subratio,sw_ldos=eil_ldos)
        else:
            gorb=(plibs.gap_orbital_from_wannier(eil_gap_file,eil_gap_iw,eil_gap_navg) #RPA/FLEX gap as form factor
                  if eil_gap_file else eil_gap_orbital)
            if eil_fs_kind=='wannier': #real Wannier-band FS + v_F (gap symmetry/multiband from gap_sym,delta0)
                sfs=plibs.build_wannier_fs(rvec,ham_r,S_r,avec,
                                           plibs.get_mu(ham_r,S_r,rvec,Arot,temp,fill),
                                           gap_sym=gap_sym,delta0=delta0,gap_orbital=gorb)
                sfk=None                    #the (int) gap_sym is baked into fs['phi']
            else:
                sfs,sfk=None,eil_fs_kind    #model FS/cylinder: the int gap_sym -> continuum harmonic
            plibs.calc_surface(eil_coupling,temp,eil_wc,gap_sym=gap_sym,beta_surf=eil_surf_beta,
                               kb=kb,sw_ldos=eil_ldos,imp_gamma=eil_imp_gamma,imp_c=eil_imp_c,h=eil_zeeman,
                               fs_kind=sfk,fs_params=eil_fs_params,fs=sfs)
    elif option==CalcMode.EILENBERGER_VORTEX: #vortex / vortex lattice via Riccati Eilenberger (model FS)
        gorb=(plibs.gap_orbital_from_wannier(eil_gap_file,eil_gap_iw,eil_gap_navg) #RPA/FLEX gap (e.g. KFe2As2) as form factor
              if eil_gap_file else eil_gap_orbital)
        if eil_fs_kind=='wannier': #real Wannier-band FS + Fermi velocities (mu from filling)
            eil_fs_obj=plibs.build_wannier_fs(rvec,ham_r,S_r,avec,
                                              plibs.get_mu(ham_r,S_r,rvec,Arot,temp,fill),
                                              gap_sym=gap_sym,delta0=delta0,gap_orbital=gorb) #gap_orbital=projection (Nagai)
            eil_fs_kw=None                  #the (int) gap_sym is baked into fs['phi']
        else:
            eil_fs_obj,eil_fs_kw=None,eil_fs_kind   #model FS/cylinder: the int gap_sym -> continuum harmonic
        if eil_vort_maxwell: #self-consistent finite-kappa vector potential A(r) (Maxwell back-reaction)
            plibs.calc_vortex_maxwell(eil_coupling,temp,eil_wc,gap_sym=gap_sym,field=eil_field,
                                      kappa=eil_kappa,kb=kb,Lxi=eil_vort_lxi,ngrid=eil_vort_ngrid)
        elif eil_vort_current: #circulating charge supercurrent j_phi(rho) of an isolated vortex
            plibs.calc_vortex_current(eil_coupling,temp,eil_wc,gap_sym=gap_sym,kb=kb,
                                      Lxi=eil_vort_lxi,ngrid=eil_vort_ngrid)
        elif eil_vort_dvector: #self-consistent triplet d-vector vortex/lattice (spin-matrix Riccati)
            dfs=(plibs.build_wannier_fs(rvec,ham_r,S_r,avec,plibs.get_mu(ham_r,S_r,rvec,Arot,temp,fill))
                 if eil_fs_kind=='wannier' else None)  #FS geometry only (d-vector channels carry the gap)
            if eil_vort_lattice_sc: #je-style TRUE periodic d-vector lattice (formulation A, square/triangular)
                plibs.calc_vortex_lattice_sc_dvector(eil_coupling,temp,eil_wc,kb=kb,field=eil_field,
                                                     lattice=eil_lattice,sub_ratio=eil_dvec_subratio,
                                                     kappa=(None if eil_kappa>=1e3 else eil_kappa),fs=dfs)
            else: #isolated vortex (eil_field=0) or circular-cell lattice (eil_field>0)
                plibs.calc_vortex_dvector(eil_coupling,temp,eil_wc,kb=kb,sub_ratio=eil_dvec_subratio,
                                          field=eil_field,fs=dfs)
        elif eil_vort_lattice_sc and eil_field_list is not None: #je-style self-consistent periodic lattice (formulation A); eil_lattice square/triangular; eil_nvortex=Vw flux quanta/cell; finite eil_kappa = London A back-reaction, >=1e3 = bare extreme
            plibs.calc_vortex_lattice_sc(eil_coupling,temp,eil_wc,gap_sym=gap_sym,
                                         field_list=eil_field_list,lattice=eil_lattice,kb=kb,fs=eil_fs_obj,
                                         kappa=(None if eil_kappa>=1e3 else eil_kappa),Vw=eil_nvortex,
                                         self_consistent_A=eil_vort_scA)
        elif eil_field_list is not None: #sweep B/Hc2 on the TRUE periodic lattice -> <N(0)>(B) (d~sqrt(B) Volovik)
            plibs.calc_vortex_lattice_periodic(eil_coupling,temp,eil_wc,gap_sym=gap_sym,
                                               field_list=eil_field_list,kappa=eil_kappa,lattice=eil_lattice,kb=kb,
                                               fs_kind=eil_fs_kw,fs_params=eil_fs_params,fs=eil_fs_obj,nflux=eil_nvortex)
        else: #single field (isolated vortex if eil_field=0, else circular-cell lattice)
            plibs.calc_vortex(eil_coupling,temp,eil_wc,gap_sym=gap_sym,kb=kb,sw_ldos=eil_ldos,
                              imp_gamma=eil_imp_gamma,imp_c=eil_imp_c,field=eil_field,h=eil_zeeman,
                              kappa=(eil_kappa if eil_vort_field else 0.0),tilt_deg=eil_vort_tilt,
                              fs_kind=eil_fs_kw,fs_params=eil_fs_params,fs=eil_fs_obj,
                              Lxi=eil_vort_lxi,ngrid=eil_vort_ngrid)
    elif option==CalcMode.CARRIER_NUM: #calc carrier number
        n_carr=plibs.calc_carrier(rvec,ham_r,S_r,avec,Nx,Ny,Nz,fill,temp)
        print(n_carr)
        print(n_carr.sum())
        plibs.get_carrier_num(Nx,rvec,ham_r,S_r,mu,Arot)
    elif option==CalcMode.CYCLOTRON_MASS: #calc cyclotron mass
        get_mass(Nx,rvec,ham_r,S_r,mu)
    elif option==CalcMode.DHVA: #plot dHvA frequency vs angle
        theta_list=np.linspace(0.,90.,40)
        get_dhva_band(Nx,rvec,ham_r,S_r,mu,theta_list)
    elif option==CalcMode.ELECTRON_MASS: #mass calc
        klist,spa_length,xticks=plibs.mk_klist(k_sets,kmesh,bvec)
        eig,uni=plibs.get_eigs(klist,ham_r,S_r,rvec)
        mass=flibs.get_mass(klist,ham_r,rvec,avec.T*ihbar,uni)*eC/emass
    elif option==CalcMode.SPECTRUM_IMPURITY: #calc spectrum with impurity
        klist,spa_length,xticks=plibs.mk_klist(k_sets,kmesh,bvec)
        rlist=plibs.gen_rlist(Nx,Ny,Nz)
        wlist=np.linspace(Emin,Emax,Nw,True)
        ham_i=ham_r
        ham_ri=ham_r  # placeholder: use host hopping for cross-species pairs
        imp_list=np.array([0])
        print("get imp Ham",flush=True)
        ham_imp=flibs.gen_imp_ham(rvec,ham_r,ham_i,ham_ri,rlist,imp_list)
        eigs,uni=sclin.eigh(ham_imp)
        print("get mu",flush=True)
        mu=plibs.calc_mu_imp(eigs,len(rlist),fill,temp)
        print('chem. pot. = %7.4f'%mu,flush=True)
        spectrum=flibs.get_imp_spectrum(uni,eigs,mu,wlist,klist,rlist)
        w,k=np.meshgrid(wlist,spa_length)
        plt.contourf(k,w,abs(spectrum.imag),100)
        plt.show()
    elif option==CalcMode.SIGMA_CPA:
        # --- CPA calculation ---
        mu=plibs.get_mu(ham_r,S_r,rvec,Arot,temp,fill)
        print(f'Temperature = {temp:10.3e} eV ({temp/kb:.2f} K)',flush=True)
        print(f'chem. pot. = {mu:.4f} eV',flush=True)
        # VA, VB are the onsite perturbation RELATIVE to the reference onsite in hamk
        # Species A (host): no perturbation -> VA = 0
        # Species B (impurity): diagonal shift
        VA = np.zeros((Norb, Norb), dtype=np.complex128)
        VB = 1.0 * np.eye(Norb, dtype=np.complex128)
        x_cpa = .5
        # sw_pp=False: periodic [0,1) mesh avoids double-counting the BZ boundary in the k-sum
        Nk, klist = plibs.gen_klist(Nx, Ny, Nz, sw_pp=False)
        hamk = flibs.gen_ham(klist, ham_r, rvec)
        # shift onsite by -mu so that w=0 corresponds to the Fermi level
        hamk[:, range(Norb), range(Norb)] -= mu
        # CPA on real axis for spectrum
        wlist = np.linspace(Emin, Emax, Nw)
        zlist_real = wlist + 1j*delta
        print(f"CPA: Norb={Norb}, Nk={Nk}, Nw={Nw}, x={x_cpa}",flush=True)
        sigma_cpa = flibs.solve_cpa(hamk, VA, VB, x_cpa, zlist_real, pp=0.5)
        print("CPA converged",flush=True)
        # --- spectrum along symmetry line ---
        Akw = calc_cpa_Akw(k_sets, kmesh, bvec, ham_r, rvec, mu,
                           sigma_cpa, wlist, zlist_real, xlabel)
        # --- CPA on Matsubara frequencies -> sigma.bin ---
        iwlist_mats = 1j*np.pi*temp*(2*np.arange(Nw, dtype=np.float64)+1)
        print(f"CPA Matsubara: Norb={Norb}, Nk={Nk}, Nw={Nw}, x={x_cpa}", flush=True)
        # reuse hamk (already shifted by -mu) computed above
        sigma_cpa_mats = flibs.solve_cpa(hamk, VA, VB, x_cpa, iwlist_mats, pp=0.5)
        print("CPA Matsubara converged", flush=True)
        # build (Norb, Norb, Nw, Nk_irr) array for sigma.bin
        # sigma_cpa_mats: (Nw, Norb, Norb) -> (Norb, Norb, Nw) -> broadcast over Nk_irr
        klist_irr, kmap, invk = flibs.gen_irr_k_TRS(Nx, Ny, Nz)
        Nk_irr = len(klist_irr)
        sigma_out = np.broadcast_to(
            sigma_cpa_mats.transpose(1, 2, 0)[:, :, :, np.newaxis],
            (Norb, Norb, Nw, Nk_irr)
        ).copy()
        # write sigma.bin in Fortran unformatted format
        # records: mu (float64), mu_OLD (float64), sigmak (Nk_irr,Nw,Norb,Norb) = (Norb,Norb,Nw,Nk_irr) C-order
        from scipy.io import FortranFile
        with FortranFile('sigma.bin', 'w') as f:
            f.write_record(np.float64(mu))
            f.write_record(np.float64(mu))   # mu_OLD = mu (no iteration history for CPA)
            f.write_record(sigma_out.flatten().astype(np.complex128))
        np.savez('self_en', sigma_out, np.float64(mu))
        print(f"CPA Matsubara self-energy written to sigma.bin and self_en.npz"
              f" (Nk_irr={Nk_irr}, Nw={Nw})", flush=True)
    elif option==CalcMode.NONLIN_ELIASHBERG:
        plibs.calc_eliashberg_eq(Nx,Ny,Nz,Nw,ham_r,S_r,rvec,chiolist,site,plist,
                                 mu,temp,gap_sym,sw_self,
                                 orb_dep,U,J,fill,sw_from_file,sw_out_self,sw_in_self,
                                 Umat if orb_dep else None,Jmat if orb_dep else None,
                                 m_diis=m_diis_num,sw_rescale=sw_rescale_flex,
                                 sw_check_only=sw_check_only)

if __name__=="__main__":
    main()
__license__="MIT"
