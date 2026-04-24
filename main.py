#!/usr/bin/env python
#-*- coding:utf-8 -*-
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
fname,ftype,brav,sw_soc='inputs/NdFeAsO.input',1,0,False
#fname,ftype,brav,sw_soc='inputs/FeS',2,0,False
#fname,ftype,brav,sw_soc='inputs/hop2.input',1,0,False
#fname,ftype,brav,sw_soc='inputs/hop2_soc.input',1,0,True
#fname,ftype,brav,sw_soc='inputs/square.hop',1,0,False
#fname,ftype,brav,sw_soc='inputs/square_soc.hop',1,0,True

sw_dec_axis=False

"""
option defines calculation modes
 0: band plot
 1: plot Dos
 2: write Fermi surface at kz (default: kz=0)
 3: write 3D Fermi surface
 4: write spectrum
 5: calc conductivity (Boltzmann Theory)
 6: calc conductivity (Purtubation Theory)
 7: calc chis spectrum with symmetry line
 8: calc chis at q-point
 9: calc chis on xy plane at Ecut
10: calc phi spectrum with symmetry line
11: calc phi on xy plane at Ecut
12: calc selfenergy using flex
13: solve linearized eliashberg equation
14: post process of gap functions
15: calc carrier num.
16: calc cycrotron mass
17: plot dHvA frequency vs angle (not implement)
18: mass calculation (not implement)
19: spectrum with impurity (not implement)
20: calc sigma_cpa (not implement)
color_option defines the meaning of color on Fermi surfaces
 0: band or mono color
 1: orbital weight settled by olist
 2: velocity size
"""
option=13
color_option=2

Nx,Ny,Nz,Nw=32,32,2,512 #k and energy(or matsubara freq.) mesh size
kmesh=200               #kmesh for spaghetti plot
kscale=[1.0,1.0,1.0]
kz=0.0
#RotMat=[[0,0,1],[0,1,0],[1,0,0]]

abc=[3.96*0.70711,3.96*0.70711,13.02*.5]
#abc=[3.68,3.68,5.03]
#alpha_beta_gamma=[90.,90.,90]
#temp=2.0e-2 #2.59e-2
tempK=500 #Kelvin
fill= 3.02
#site_prof=[5]

Emin,Emax=-3,3
delta=3.0e-2
Ecut=1.0e-2
tau_const=100
#olist=[0,0,0]
olist=[0,[1,2],3]
#olist=[[0,3],[1,4],[2,5]]
#olist=[[0,4],[1,2,5,6],[3,7]]
#U,J= 0.8, 0.1
U,J=1.2,0.15
#U,J=1.8,0.225
#0:s,1:dx2-y2,2:spm,3:dxy,-1:px,-2:py
gap_sym=2

#mu0=9.85114560061123
#k_sets=[[0., 0., 0.],[.5, 0., 0.],[.5, .5, 0.]]
#xlabel=[r'$\Gamma$','X','M']
#m_diis_num=2
at_point=[ 0., .5, 0.]
orb_dep=False  #use orbital dependence U,J
sw_unit=True    #set unit values unity (False) or not (True)
sw_tdf=False
sw_omega=False #True: real freq, False: Matsubara freq.
sw_rescale_flex=True #True: rescale self energy to make max|Sigma|~U, False: no rescaling
sw_self=False  #True: use calculated self energy for spectrum band plot
sw_out_self=True
sw_in_self=False
sw_from_file=False
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
try:
    alpha_beta_gamma
except NameError:
    alpha_beta_gamma=[90.,90.,90]
alatt=np.array(abc)
deg=np.array(alpha_beta_gamma)
if option in {0,4,7,12,20}:
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
    hbar=scconst.physical_constants['Planck constant over 2 pi in eV s'][0]
    ihbar=1.0e-10/hbar
    kb=scconst.physical_constants['Boltzmann constant in eV/K'][0]
    eC=scconst.e
    tau_unit=1.e-15
    emass=scconst.m_e
else:
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
    pass
try:
    m_diis_num
except NameError:
    m_diis_num=5
if option in {2,16}:
    try:
        RotMat
        RotMat=np.array(RotMat)
    except NameError:
        print('No RotMat')
        RotMat=np.eye(3)
#------------------------ define functions -------------------------------------------
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
            norm=np.sqrt((abs(cl)**2).sum(axis=0))
            # Check for norm close to zero (use tolerance for floating-point comparison)
            if np.any(norm < 1e-14):
                print("Warning: plot_band found norm close to zero. Skipping",flush=True)
                continue
            cls=cl/norm
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

def plot_3d_surf(fspolys,fscenters,fscolors,surface_opt,kscale):
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
    plibs.BZedge(0)
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
                  delta:float,Nw:int,sw_self=True,selfen=None):
    klist,spa_length,xticks=plibs.mk_klist(k_sets,kmesh,bvec)
    eig,uni=plibs.get_eigs(klist,ham_r,S_r,rvec)
    wlist=np.linspace(Emin,Emax,Nw)
    if sw_self:
        Gk0=flibs.gen_Green0(eig,uni,mu,temp,Nw)
        iwlist=np.pi*temp*(2*np.linspace(0,Nw,Nw,False)+1)*1j
        Gk=flibs.pade_with_trace(selfen,iwlist,wlist-1j*delta).imag
        w,x=np.meshgrid(wlist,spa_length)
        plt.contourf(x,w,Gk,cmap=plt.hot())
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
    gsp=(1.0 if with_spin else 2.0) #spin weight
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
    sigma_hall=flibs.calc_sigmahall(eig,vk,imass/eC,kweight,tau,temp,mu)
    K0,K1,K2=flibs.calc_Kn(eig,vk,kweight,temp,mu,tau)
    print(f"sigma_hall={sigma_hall:.6e}, K0[0,0]={K0[0,0]:.6e}, K0[1,1]={K0[1,1]:.6e}")
    # Check if K0 diagonal elements are non-zero (use tolerance for floating-point comparison)
    tol=1e-14
    if abs(K0[0,0]) < tol or abs(K0[1,1]) < tol:
        print("Error: K0 diagonal elements are too small. Cannot compute Hall coefficient",flush=True)
        return
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
        wlist=np.linspace(eig.min()-mu,eig.max()-mu,Nw)
        dw=(eig.max()-eig.min())/Nw
        dfermi=0.25*(1.-np.tanh(0.5*(wlist-mu)/temp)**2)/temp
        K0=(dfermi*tdf.T).T.sum(axis=0)*dw
        K1=(dfermi*(wlist-mu)*tdf.T).T.sum(axis=0)*dw
        K2=(dfermi*(wlist-mu)**2*tdf.T).T.sum(axis=0)*dw
        plt.plot(wlist-mu,tdf[:,0,0])
        plt.show()
    else:
        K0,K1,K2=flibs.calc_Kn(eig,vk,kweight,temp,mu,tau)
    sigma=gsp*tau_unit*eC*K0*iNV #a e is canceled effect of eV2J
    kappa=gsp*tau_unit*eC*kb*K2*iNV*itemp

    # Handle sclin.inv() failures
    try:
        K0_inv = sclin.inv(K0)
        kappa2=gsp*tau_unit*eC*kb*(K2-K1.dot(K0_inv.dot(K1)))*iNV*itemp
        Seebeck=-kb*K0_inv.dot(K1)*itemp
        Pertier=K1.dot(K0_inv)
    except np.linalg.LinAlgError:
        print("Error: K0 is a singular matrix. Cannot compute Seebeck coefficient and Pertier",flush=True)
        kappa2=kappa.copy()
        Seebeck=np.zeros_like(K0)
        Pertier=np.zeros_like(K0)

    try:
        sigma_inv = sclin.inv(sigma*temp)
        Lorenz=kb*kappa*sigma_inv
        Lorenz2=kb*kappa2*sigma_inv
    except np.linalg.LinAlgError:
        print("Warning: sigma*temp is a singular matrix. Cannot compute Lorenz coefficient",flush=True)
        Lorenz=np.zeros_like(sigma)
        Lorenz2=np.zeros_like(sigma)
    sigmaS=gsp*tau_unit*kb*eC*K1*iNV*itemp
    PF=sigma*Seebeck**2
    print('sigma matrix (S/m)',flush=True)
    for sig in sigma.round(10):
        print(f" {sig[0]:10.3e} {sig[1]:10.3e} {sig[2]:10.3e}",flush=True)
    print('kappa matrix (K22 only) (W/m/K)',flush=True)
    for kap in kappa.round(10):
        print(f" {kap[0]:10.3e} {kap[1]:10.3e} {kap[2]:10.3e}",flush=True)
    print('kappa matrix (full) (W/m/K)',flush=True)
    for kap in kappa2.round(10):
        print(f" {kap[0]:10.3e} {kap[1]:10.3e} {kap[2]:10.3e}",flush=True)
    print('sigmaS matrix (A/m/K)',flush=True)
    for sig in sigmaS.round(10):
        print(f" {sig[0]:10.3e} {sig[1]:10.3e} {sig[2]:10.3e}",flush=True)
    print('Seebeck matrix (V/K)',flush=True)
    for seeb in Seebeck.round(10):
        print(f" {seeb[0]:10.3e} {seeb[1]:10.3e} {seeb[2]:10.3e}",flush=True)
    print('Pertier matrix (V)',flush=True)
    for per in Pertier.round(10):
        print(f" {per[0]:10.3e} {per[1]:10.3e} {per[2]:10.3e}",flush=True)
    print('Lorenz matrix (K22 only) (Wohm/K^2) (fe 2.44e-8)',flush=True)
    for lor in Lorenz.round(10):
        print(f" {lor[0]:9.2e} {lor[1]:9.2e} {lor[2]:9.2e}",flush=True)
    print('Lorenz matrix? (full) (Wohm/K^2)',flush=True)
    for lor in Lorenz2.round(10):
        print(f" {lor[0]:9.2e} {lor[1]:9.2e} {lor[2]:9.2e}",flush=True)
    print('Power Factor (SA/m^2/K)',flush=True)
    for pofa in PF.round(8):
        print(f" {pofa[0]:10.3e} {pofa[1]:10.3e} {pofa[2]:10.3e}",flush=True)

def calc_conductivity_lrt(rvec,ham_r,S_r,avec,Nx:int,Ny:int,Nz:int,fill:float,
                         temp:float,Nw:int,delta,with_spin=False):
    '''
    calculation of linear response theory
    electric conductivity of LRT correponds to Boltzmann then delta~O(10-1) (tau~1fs) at 300K
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
    print(f'tempreture = {temp/kb:.3f} K',flush=True)
    print(f'about tau = {1.e15*hbar/delta:.3f} fs',flush=True)
    print(f'delta = {delta:9.3e}',flush=True)
    L11,L12,L22,wlist=plibs.get_conductivity(mu,temp,eig,vk,Nw,Emax,delta)
    sigmaconst=gsp*hbar*eC/Vuc*1.0e30
    kappaSconst=sigmaconst*kb/temp
    sigma=sigmaconst*L11
    kappa=kappaSconst*L22
    sigmaS=kappaSconst*L12
    # Handle sclin.inv() failures
    try:
        Seebeck=np.array([-sclin.inv(s).dot(sS) for s,sS in zip(sigma,sigmaS)])
    except np.linalg.LinAlgError:
        print("Warning: sigma is a singular matrix. Cannot compute Seebeck coefficient",flush=True)
        Seebeck=np.zeros_like(sigma)
    print('sigma matrix (S/m)',flush=True)
    for sig in sigma[0].real.round(10):
        print(f" {sig[0]:10.3e} {sig[1]:10.3e} {sig[2]:10.3e}",flush=True)
    print('kappa matrix (L22 only) (W/m/K)',flush=True)
    for kap in kappa[0].real.round(10):
        print(f" {kap[0]:10.3e} {kap[1]:10.3e} {kap[2]:10.3e}",flush=True)
    print('sigmaS matrix (A/m/K)',flush=True)
    for sig in sigmaS[0].real.round(10):
        print(f" {sig[0]:10.3e} {sig[1]:10.3e} {sig[2]:10.3e}",flush=True)
    print('Lorenz number (Wohm/K^2) (fe 2.44e-8)',flush=True)
    try:
        Lorenz_matrix=(kb*kappa[0]*sclin.inv(sigma[0]*temp)).real.round(10)
        for lor in Lorenz_matrix:
            print(f" {lor[0]:9.2e} {lor[1]:9.2e} {lor[2]:9.2e}",flush=True)
    except np.linalg.LinAlgError:
        print("Warning: Failed to compute Lorenz coefficient (singular matrix)",flush=True)
    print('Seebeck coefficient matrix (V/K)',flush=True)
    for seeb in Seebeck[0].real.round(10):
        print(f" {seeb[0]:10.3e} {seeb[1]:10.3e} {seeb[2]:10.3e}",flush=True)
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

def calc_chis_spectrum(mu:float,temp:float,Smat,klist,qlist,chiolist,eig,uni,spa_length,
                       Nw:int,Emax:float,delta:float):
    print("calculate spn susceptibility",flush=True)
    chisw,chisw_orb,wlist=plibs.chis_spectrum(mu,temp,Smat,klist,qlist,chiolist,eig,uni,Nw,Emax,delta)
    w,sp=np.meshgrid(wlist,spa_length)
    try:
        with open('chis.dat','w') as f:
            for ww,ssp,chis in zip(w,sp,chisw):
                for www,sssp,chi in zip(ww,ssp,chis):
                    f.write(f'{sssp:8.4f} {www:8.4f} {chi.imag:9.4f}\n')
                f.write('\n')
    except IOError as e:
        print(f"Error: Failed to write 'chis.dat': {e}",flush=True)
    for i,chiso in enumerate(chisw_orb.T):
        try:
            with open(f'chis_{i}.dat','w') as f:
                for ww,ssp,chis in zip(w,sp,chiso.T):
                    for www,sssp,chi in zip(ww,ssp,chis):
                        f.write(f'{sssp:8.4f} {www:8.4f} {chi.imag:9.4f}\n')
                    f.write('\n')
        except IOError as e:
            print(f"Error: Failed to write 'chis_{i}.dat': {e}",flush=True)
            continue
    return(w,sp,chisw)

def calc_phi_spectrum(mu:float,temp:float,klist,qlist,chiolist,eig,uni,spa_length,Nw:int,Emax:float,delta:float):
    print("calculate sc susceptibility",flush=True)
    phiw,phiw_orb,wlist=plibs.phi_spectrum(mu,temp,klist,qlist,chiolist,eig,uni,Nw,Emax,delta)
    w,sp=np.meshgrid(wlist,spa_length)
    try:
        with open('phi.dat','w') as f:
            for ww,ssp,phi in zip(w,sp,phiw):
                for www,sssp,ph in zip(ww,ssp,phi):
                    f.write(f'{sssp:8.4f} {www:8.4f} {ph.imag:9.4f}\n')
                f.write('\n')
    except IOError as e:
        print(f"Error: Failed to write 'phi.dat': {e}",flush=True)
    return(w,sp,phiw)

def calc_flex(Nx:int,Ny:int,Nz:int,Nw:int,ham_r,S_r,rvec,mu:float,temp:float,olist,site,eps=1.0e-4,pp=0.5,m_diis=5,sw_rescale:bool=True):
    klist,kmap,invk=flibs.gen_irr_k_TRS(Nx,Ny,Nz)
    eig,uni=plibs.get_eigs(klist,ham_r,S_r,rvec)
    ham_k=flibs.gen_ham(klist,ham_r,rvec)
    if orb_dep:
        Smat,Cmat=flibs.gen_SCmatrix_orb(olist,site,Umat,Jmat)
    else:
        Smat,Cmat=flibs.gen_SCmatrix(olist,site,U,J)
    sigmak,mu_self=flibs.mkself(Smat,Cmat,kmap,invk,olist,ham_k,eig,uni,mu,fill,temp,Nw,Nx,Ny,Nz,sw_out_self,sw_in_self,eps=eps,pp=pp,m_diis=m_diis,sw_rescale=sw_rescale)
    if sw_out_self:
        np.savez('self_en',sigmak,mu_self)
        plibs.output_self_wannier(sigmak,mu_self,kmap,invk,Nx,Ny,Nz,Nw,temp)

def calc_flex_soc(Nx:int,Ny:int,Nz:int,Nw:int,ham_r,S_r,rvec,mu:float,temp:float,olist,slist,invs,site,eps=1.0e-4,pp=0.5,m_diis=5,sw_rescale:bool=True):
    klist,kmap,invk=flibs.gen_irr_k_TRS(Nx,Ny,Nz)
    eig,uni=plibs.get_eigs(klist,ham_r,S_r,rvec)
    ham_k=flibs.gen_ham(klist,ham_r,rvec)
    if orb_dep:
        Vmat=flibs.gen_Vmatrix_orb(olist,slist,site,invs,Umat,Jmat)
    else:
        Vmat=flibs.gen_Vmatrix(olist,slist,site,invs,U,J)
    sigmak,mu_self=flibs.mkself_soc(Vmat,kmap,invk,invs,olist,slist,ham_k,eig,uni,mu,fill,temp,
                                    Nw,Nx,Ny,Nz,sw_out_self,sw_in_self,eps=eps,pp=pp,m_diis=m_diis,sw_rescale=sw_rescale)
    if sw_out_self:
        np.savez('self_en',sigmak,mu_self)
        plibs.output_self_wannier(sigmak,mu_self,kmap,invk,Nx,Ny,Nz,Nw,temp)


def output_gap_function(invk,kmap,gap,uni,plist,gap_sym,soc=False,invs=None,slist=None,sw_orb=False):
    if sw_orb:
        if soc:
            gapb=gap[:,:,0,:]
        else:
            gapb=flibs.remap_gap(gap[:,:,0,:],plist,invk,gap_sym)
    else:
        if soc:
            gapb=flibs.conv_delta_orb_to_band_soc(gap,uni,invk,invs,slist)
        else:
            gapb=flibs.conv_delta_orb_to_band(gap,uni,invk,plist,gap_sym)
    print('output gap function')
    for iorb in range(len(gapb)):
        for jorb in range(len(gapb)):
            try:
                with open(f'gap_{iorb+1}{jorb+1}.dat','w') as f:
                    for i,km in enumerate(kmap):
                        if km[2]==0:
                            f.write(f'{km[0]:3} {km[1]:3} {gapb[iorb,jorb,i].real:12.8f} {gapb[iorb,jorb,i].imag:12.8f}\n')
                            if km[0]==Nx-1:
                                f.write('\n')
            except IOError as e:
                print(f"Error: Failed to write 'gap_{iorb+1}{jorb+1}.dat': {e}",flush=True)
                continue
    return(0)

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
        no=int(len(slist)/2)
        if gap_sym>=0:
            Fks=Fk[:no,no:,0,:]-Fk[no:,:no,0,:]
            Fktr=np.array([f.diagonal().sum() for f in Fks.T])
        else:
            Fkt=Fk[:no,no:,0,:]+Fk[no:,:no,0,:]
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
    info=output_gap_function(invk,kmap,gap,uni,plist,gap_sym,sw_soc,invs,slist)

def calc_lin_eliashberg_eq(Nx:int,Ny:int,Nz:int,Nw:int,ham_r,S_r,rvec,chiolist,site,plist,
                           mu:float,temp:float,gap_sym:int,sw_self:bool,eps=1.0e-4,pp=0.5,m_diis=5,sw_rescale:bool=True):
    klist,kmap,invk=flibs.gen_irr_k_TRS(Nx,Ny,Nz)
    #weight=flibs.gen_kpoint_weight(invk,len(klist))
    eig,uni=plibs.get_eigs(klist,ham_r,S_r,rvec)
    if orb_dep:
        Smat,Cmat=flibs.gen_SCmatrix_orb(chiolist,site,Umat,Jmat)
    else:
        Smat,Cmat=flibs.gen_SCmatrix(chiolist,site,U,J)
    if sw_self:
        ham_k=flibs.gen_ham(klist,ham_r,rvec)
        if sw_from_file:
            try:
                npz=np.load('self_en.npz')
                sigmak,mu_self=npz['arr_0'],npz['arr_1']
            except FileNotFoundError:
                print("Error: 'self_en.npz' not found",flush=True)
                return
        else:
            sigmak,mu_self=flibs.mkself(Smat,Cmat,kmap,invk,chiolist,ham_k,eig,uni,
                                        mu,fill,temp,Nw,Nx,Ny,Nz,sw_out_self,sw_in_self,eps=eps,pp=pp,m_diis=m_diis,sw_rescale=sw_rescale)
        print(f'chem. pot. with self= {mu:.4f} eV',flush=True)
        Gk=flibs.gen_green(sigmak,ham_k,mu_self,temp)
    else:
        Gk=flibs.gen_Green0(eig,uni,mu,temp,Nw)
    init_delta=plibs.get_initial_gap(kmap,klist,len(eig.T),gap_sym)
    sw_eig=True
    chi=flibs.get_chi0(Smat,Cmat,Gk,chiolist,kmap,invk,temp,Nx,Ny,Nz)
    chis,chic=flibs.get_chis_chic(chi,Smat,Cmat)
    chisq=flibs.get_eig_or_tr_chi(chis,invk,sw_eig)
    chicq=flibs.get_eig_or_tr_chi(chic,invk,sw_eig)
    try:
        with open('chis.dat','w') as f, open('chic.dat','w') as f2:
            for i,k in enumerate(kmap):
                if k[2]==0.0:
                    f.write(f'{k[0]:6.4f} {k[1]:6.4f} {chisq[i].real:11.4e}\n')
                    f2.write(f'{k[0]:6.4f} {k[1]:6.4f} {chicq[i].real:11.4e}\n')
    except IOError as e:
        print(f"Error: Failed to write 'chis.dat' or 'chic.dat': {e}",flush=True)
    gap=flibs.linearized_eliashberg(chi,Gk,uni,init_delta,Smat,Cmat,chiolist,plist,kmap,invk,Nx,Ny,Nz,temp,gap_sym)
    if sw_out_self:
        np.save('gap',gap)
        plibs.output_gap_wannier(gap,kmap,invk,Nx,Ny,Nz,Nw,temp)
    info=output_gap_function(invk,kmap,gap,uni,plist,gap_sym)

def calc_lin_eliash_soc(Nx:int,Ny:int,Nz:int,Nw:int,ham_r,S_r,rvec,
                        mu:float,temp:float,chiolist,slist,plist,invs,site,eps=1.0e-4,pp=0.5,m_diis=5,sw_rescale:bool=True):
    klist,kmap,invk=flibs.gen_irr_k_TRS(Nx,Ny,Nz)
    eig,uni=plibs.get_eigs(klist,ham_r,S_r,rvec)
    if orb_dep:
        Vmat=flibs.gen_Vmatrix_orb(chiolist,slist,site,invs,Umat,Jmat)
    else:
        Vmat=flibs.gen_Vmatrix(chiolist,slist,site,invs,U,J)
    if sw_self:
        ham_k=flibs.gen_ham(klist,ham_r,rvec)
        if sw_from_file:
            try:
                npz=np.load('self_en.npz')
                sigmak,mu_self=npz['arr_0'],npz['arr_1']
            except FileNotFoundError:
                print("Error: 'self_en.npz' not found",flush=True)
                return
        else:
            sigmak,mu_self=flibs.mkself_soc(Vmat,kmap,invk,invs,chiolist,slist,ham_k,eig,uni,mu,fill,temp,
                                            Nw,Nx,Ny,Nz,sw_out_self,sw_in_self,eps=eps,pp=pp,m_diis=m_diis,sw_rescale=sw_rescale)
        print(f'chem. pot. with self= {mu:.4f} eV',flush=True)
        Gk=flibs.gen_green(sigmak,ham_k,mu_self,temp)
    else:
        Gk=flibs.gen_Green0(eig,uni,mu,temp,Nw)
    sw_eig=True
    chi,sgnsig,sgnsig2,invschi=flibs.get_chi0_soc(Vmat,Gk,chiolist,slist,kmap,invk,invs,temp,Nx,Ny,Nz)
    chic,chiszz,chispm=flibs.get_chis_chic_soc(chi,Vmat,chiolist,slist,invs)
    chiszzq=flibs.get_eig_or_tr_chi(chiszz,invk,sw_eig)
    chispmq=flibs.get_eig_or_tr_chi(chispm,invk,sw_eig)
    chicq=flibs.get_eig_or_tr_chi(chic,invk,sw_eig)
    try:
        with open('chis.dat','w') as f, open('chic.dat','w') as f2:
            for i,k in enumerate(kmap):
                if k[2]==0.0:
                    f.write(f'{k[0]:6.4f} {k[1]:6.4f} {chiszzq[i].real:11.4e} {chispmq[i].real:11.4e}\n')
                    f2.write(f'{k[0]:6.4f} {k[1]:6.4f} {chicq[i].real:11.4e}\n')
    except IOError as e:
        print(f"Error: Failed to write 'chis.dat' or 'chic.dat': {e}",flush=True)
    init_delta=plibs.get_initial_gap(kmap,klist,len(slist),gap_sym)
    gap=flibs.linearized_eliashberg_soc(chi,Gk,uni,init_delta,Vmat,sgnsig,sgnsig2,plist,slist,chiolist,
                                        kmap,invk,invs,invschi,Nx,Ny,Nz,temp,gap_sym)
    if sw_out_self:
        np.save('gap',gap)
        plibs.output_gap_wannier(gap,kmap,invk,Nx,Ny,Nz,Nw,temp)
    info=output_gap_function(invk,kmap,gap,uni,plist,gap_sym,True,invs,slist)

def get_carrier_num(kmesh,rvec,ham_r,S_r,mu:float,Arot):
    Nk,eig,kwieght=plibs.get_emesh(kmesh,kmesh,kmesh,ham_r,S_r,rvec,Arot)
    if Nk <= 0:
        print("Error: Number of k-points (Nk) is non-positive",flush=True)
        return
    fill=0.0
    for i,en in enumerate(eig.T-mu):
        num_hole=float(np.where(en>0)[0].size)/Nk
        num_particle=float(np.where(en<=0)[0].size)/Nk
        print(i+1,round(num_hole,4),round(num_particle,4),flush=True)
        fill+=num_particle
    print(f'sum of electrons is {fill:.4f}',flush=True)

def get_mu(ham_r,S_r,rvec,Arot,temp:float,kmesh=40)->float:
    # Parameter validation
    if temp < 0:
        print("Error: Temperature (temp) is negative",flush=True)
        return None
    if kmesh <= 0:
        print("Error: k-mesh size (kmesh) is non-positive",flush=True)
        return None

    print("calc chem. pot.",flush=True)
    print(f"band filling = {fill:.4f}",flush=True)
    Nk,eig,kweight=plibs.get_emesh(kmesh,kmesh,kmesh,ham_r,S_r,rvec,Arot)
    mu=plibs.calc_mu(eig,Nk,fill,temp)
    return mu

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
            dSdE_SI=(Sp-Sm)/(2.*de)*1.e20
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

    # split each band into orbit branches:
    # at each theta there are multiple extremal F values (kz=0, kz=pi/2, interior).
    # sort F ascending at each theta and connect same-rank points across theta -> branches.
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
    if tempK < 0:
        print(f"Error: Temperature (tempK) is negative ({tempK} K)",flush=True)
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
        rvec,ham_r,S_r,no,Nr=plibs.import_MLO_hoppings(fname)
    else:
        rvec,ham_r,no,Nr=plibs.import_hoppings(fname,ftype)
        S_r=[]
    plist=flibs.get_plist(rvec,ham_r)
    print("Effective parity of wanier functions:", plist)

    # ===== Parameter checks dependent on number of orbitals =====
    print(f"Number of orbital = {no}",flush=True)

    # Validate olist
    for i, ol in enumerate(olist):
        if isinstance(ol, int):
            if ol < 0 or ol >= no:
                print(f"Error: olist[{i}]={ol} is invalid. Valid range: 0-{no-1}",flush=True)
                return
        elif isinstance(ol, list):
            for j, olj in enumerate(ol):
                if olj < 0 or olj >= no:
                    print(f"Error: olist[{i}][{j}]={olj} is invalid. Valid range: 0-{no-1}",flush=True)
                    return

    # Check fill does not exceed number of bands
    if fill > no:
        print(f"Error: filling={fill} exceeds number of bands={no}. Valid range: 0 < fill <= {no}",flush=True)
        return
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
    opstr=["calculate band structure","calculate Dos","plot 2D Fermi surface",
           "plot 3D Fermi surface","calculate spectrum",
           "calculate conductivities using Boltzmann theory",
           "calculate conductivities with linear response","calculate chis spectrum",
           "calculate chis at q-point","calculate chis qmap at Ecut",
           "calc phi spectrum with symmetry line","calc phi on xy plane at Ecut",
           "calc self energy","solve linearized eliashberg equation",
           "gap_function","calculate carrier number","calculate cycrtron mass",
           "plot dHvA frequency","calculate electron mass","spectrum with impurity",
           "calculate sigma_cpa"]
    cstr=["no color",'orbital weight','velocity size']
    if omp_check: #OMP properties
        print("OpenMP mode",flush=True)
        print(f"Number of OpenMP threads = {omp_num}",flush=True)
    print(f"calc mode {option}: "+opstr[option],flush=True)

    # ===== Additional input validation (lower priority) =====
    # J > U warning
    if J > U:
        print(f"Warning: J={J} > U={U} is physically unusual. Please verify",flush=True)

    # Validate gap_sym during eliashberg/flex calculations
    if option in {12, 13, 14}:  # eliashberg/flex calculations
        if gap_sym not in {-1, 0, 1, 2, 3}:
            print(f"Warning: gap_sym={gap_sym} is non-standard. Common values: -1,0,1,2,3",flush=True)

    # Validate at_point for chis at q-point calculation
    if option in {8}:  # chis at q-point
        if len(at_point) != 3:
            print(f"Error: at_point has {len(at_point)} elements. Required: 3 elements",flush=True)
            return
        if any(p < 0 or p > 1 for p in at_point):
            print(f"Warning: at_point={at_point} values outside [0,1] range. May not be normalized by reciprocal lattice",flush=True)
    # ==========================================
    if option in {0,2,3}:
        print("color mode: "+cstr[color_option],flush=True)
    print("Hamiltonian name is "+fname,flush=True)
    print(f"Number of orbital = {no}",flush=True)
    if (orb_dep==False) and option in {7,8,9,12,13}: #write constant U,J
        print(f'U= {U:5.2f} and J= {J:5.3f}')
    if option in {7,8,9,10,11,12,13}:
        """ chiolist is the list of orbital properties of index on chi """
        try:
            chiolist
        except NameError:
            try:
                site_prof
            except NameError:
                site_prof=[1] #one site (len(site_prof)=1)
            chiolist,site=plibs.get_chi_orb_list(len(ham_r[0]),site_prof)        
    if option in {0,4}:
        if sw_gen_sym:
            print('generate symmetry line',flush=True)
        print(f'kmesh = {kmesh}',flush=True)
    elif option in {2,3,9,11,15,16}:
        print(f'Number of k-mesh = {Nx}',flush=True)
    else:
        print(f'k-mesh is {Nx} {Ny} {Nz}',flush=True)
    if option in {12,13}:
        print(f'Number of Matsubara freq. = {Nw}',flush=True)
    print("Lattice Vector (Angstrom)",flush=True)
    for i,a in enumerate(avec):
        print(f"a{i+1}: {a[0]:7.4f} {a[1]:7.4f} {a[2]:7.4f}",flush=True)
    print("Reciprocal Lattice Vector (Angstrom^-1)",flush=True)
    for i,b in enumerate(bvec):
        print(f"b{i+1}: %7.4f %7.4f %7.4f"%tuple(b),flush=True)
    if option in {5,6,19,20}: #conductivity (5,6) and impurity (18) functions calc or set mu themself
        pass
    else:
        if sw_calc_mu:
            mu=get_mu(ham_r,S_r,rvec,Arot,temp)
        else:
            mu=mu0
            print('use fixed mu')
        print(f'Temperature = {temp:10.3e} eV ({temp/kb:.2f} K)',flush=True)
        print(f'chem. pot. = {mu:.4f} eV',flush=True)
    if option==0: #plot band
        klist,spa_length,xticks=plibs.mk_klist(k_sets,kmesh,bvec)
        eig,uni0=plibs.get_eigs(klist,ham_r,S_r,rvec)
        uni=np.array([u.T for u in uni0]) #rotate uni(k,band,orb) to uni(k,orb,band)
        plot_band(eig.T-mu,spa_length,xlabel,xticks,uni.T,olist,(False if color_option==0 else True))
    elif option==1: #plot dos
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
            tmp_ol_b=[]
            for i in range(len(olist)):
                tmp_ol_b+=olist[i]if type(olist[i])==list else [olist[i]]
            plt.fill_between(wlist,Dos[tmp_ol_b].sum(axis=0),Dos.sum(axis=0),color=clist[0])
            for i in range(len(olist)-1):
                tmp_ol_u=[]
                tmp_ol_b=[]
                for j in range(i,len(olist)):
                    tmp_ol_u+=olist[j] if type(olist[j])==list else [olist[j]]
                for j in range(i+1,len(olist)):
                    tmp_ol_b+=olist[j] if type(olist[j])==list else [olist[j]]
                plt.fill_between(wlist,Dos[tmp_ol_b].sum(axis=0),Dos[tmp_ol_u].sum(axis=0),color=clist[i+1])
            plt.fill_between(wlist,0,Dos[tmp_ol_b].sum(axis=0),color=clist[len(olist)])
        plt.xlim(Emin,Emax)
        plt.ylim(0,max(Dos.sum(axis=0))*1.2)
        plt.show()
    elif option==2: #2D Fermi surface plot
        klist,blist=plibs.mk_kf(Nx,rvec,ham_r,S_r,RotMat,mu,kz)
        clist=plibs.get_colors(klist,blist,ihbar*avec.T,rvec,ham_r,S_r,olist,color_option,True)
        plot_FS(clist,klist,color_option)
    elif option==3: #3D Fermi surface plot
        polys,centers,blist=plibs.gen_3d_surf_points(Nx,rvec,ham_r,S_r,mu,kscale)
        fspolys,fscenters,fscolors=set_init_3dfsplot(color_option,polys,centers,blist,avec,rvec,ham_r,S_r,olist)
        plot_3d_surf(fspolys,fscenters,fscolors,color_option,kscale)
    elif option==4: #plot spectrum
        plot_spectrum(k_sets,xlabel,kmesh,bvec,mu,ham_r,S_r,rvec,Emin,Emax,delta,Nw,sw_self)
    elif option==5: #calc conductivity
        get_hall_coe(rvec,ham_r,S_r,avec,Nx,Ny,Nz,fill,temp,tau_const)
        calc_conductivity_Boltzmann(rvec,ham_r,S_r,avec,Nx,Ny,Nz,fill,temp,tau_const)
    elif option==6: #calc_optical conductivity
        calc_conductivity_lrt(rvec,ham_r,S_r,avec,Nx,Ny,Nz,fill,temp,Nw,delta)
    elif option in {7,8,9,10,11}: #calc_chis_spectrum
        print("calculate electron energy",flush=True)
        Nk,klist,eig,uni,kweight=plibs.get_emesh(Nx,Ny,Nz,ham_r,S_r,rvec,avec,sw_uni=True)
        if option in {7,8,9}:
            print("generate coulomb vertex matrix S")
            if orb_dep:
                Smat,Cmat=flibs.gen_SCmatrix_soc(chiolist,site,Umat,Jmat)
            else:
                Smat,Cmat=flibs.gen_SCmatrix(chiolist,site,U,J)
        if option in {7,10}: #chis/phi spectrum with symmetry line
            print("generate qlist",flush=True)
            qlist,spa_length,xticks=plibs.mk_qlist(k_sets,Nx,Ny,Nz,bvec)
            if option==7:
                w,sp,sus=calc_chis_spectrum(mu,temp,Smat,klist,qlist,chiolist,eig,uni,spa_length,Nw,Emax,delta)
                print("write chis spectrum in png file",flush=True)
                susfname='chis_spec.png'
            elif option==10:
                w,sp,sus=calc_phi_spectrum(mu,temp,klist,qlist,chiolist,eig,uni,spa_length,Nw,Emax,delta)
                print("write phi spectrum in png file",flush=True)
                susfname='phi_spec.png'
            plt.contourf(sp,w,abs(sus.imag),100)
            plt.colorbar()
            plt.hot()
            plt.savefig(fname=susfname,dpi=300)
        elif option==8:
            q_point=np.array(at_point)
            chis,chis_orb,wlist=plibs.chis_q_point(q_point,eig,uni,Emax,Nw,mu,temp,Smat,klist,chiolist,delta)
            plt.plot(wlist,chis.imag)
            plt.show()
        else:
            if option==9: #chis spectrum ecut plane
                sus,chi0,qx,qy=plibs.chis_qmap(Nx,Ny,Ecut,mu,temp,Smat,klist,chiolist,eig,uni,idelta=1.e-3)
                plt.contourf(qx,qy,abs(chi0.imag),100)
                plt.colorbar()
                plt.jet()
                plt.show()
                susfname='chismap.png'
            elif option==11:
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
    elif option in {12,13,14}:
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
            if option==12:
                calc_flex_soc(Nx,Ny,Nz,Nw,ham_r,S_r,rvec,mu,temp,chiolist,slist,invs,site,m_diis=m_diis_num,sw_rescale=sw_rescale_flex)
            elif option==13:
                calc_lin_eliash_soc(Nx,Ny,Nz,Nw,ham_r,S_r,rvec,mu,temp,chiolist,slist,plist,invs,site)
            elif option==14:
                output_Fk(Nx,Ny,Nz,Nw,ham_r,S_r,rvec,plist,mu,temp,sw_self,sw_soc,invs,slist,gap_sym)
        else: #without soc
            if option==12: #calc self-energy using flex
                calc_flex(Nx,Ny,Nz,Nw,ham_r,S_r,rvec,mu,temp,chiolist,site,m_diis=m_diis_num,sw_rescale=sw_rescale_flex)
            elif option==13: #calc gap function
                calc_lin_eliashberg_eq(Nx,Ny,Nz,Nw,ham_r,S_r,rvec,chiolist,site,plist,mu,temp,gap_sym,sw_self)
            elif option==14:
                output_Fk(Nx,Ny,Nz,Nw,ham_r,S_r,rvec,plist,mu,temp,sw_self)
    elif option==15: #calc carrier number
        n_carr=plibs.calc_carrier(rvec,ham_r,S_r,avec,Nx,Ny,Nz,fill,temp)
        print(n_carr)
        print(n_carr.sum())
        get_carrier_num(Nx,rvec,ham_r,S_r,mu,Arot)
    elif option==16: #calc cycrtron mass
        get_mass(Nx,rvec,ham_r,S_r,mu)
    elif option==17: #plot dHvA frequency vs angle
        theta_list=np.linspace(0.,90.,40)
        get_dhva_band(Nx,rvec,ham_r,S_r,mu,theta_list)
    elif option==18: #mass calc
        klist,spa_length,xticks=plibs.mk_klist(k_sets,kmesh,bvec)
        eig,uni=plibs.get_eigs(klist,ham_r,S_r,rvec)
        mass=flibs.get_mass(klist,ham_r,rvec,avec.T*ihbar,uni)*eC/emass
    elif option==19: #calc spectrum with impurity
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
    elif option==20:
        # --- CPA calculation ---
        mu=get_mu(ham_r,S_r,rvec,Arot,temp)
        print(f'Temperature = {temp:10.3e} eV ({temp/kb:.2f} K)',flush=True)
        print(f'chem. pot. = {mu:.4f} eV',flush=True)
        # VA, VB are the onsite perturbation RELATIVE to the reference onsite in hamk
        # Species A (host): no perturbation -> VA = 0
        # Species B (impurity): diagonal shift
        VA = np.zeros((no, no), dtype=np.complex128)
        VB = 1.0 * np.eye(no, dtype=np.complex128)
        x_cpa = .5
        Nk, klist = plibs.gen_klist(Nx, Ny, Nz)
        hamk = flibs.gen_ham(klist, ham_r, rvec)
        # shift onsite by -mu so that w=0 corresponds to the Fermi level
        hamk[:, range(no), range(no)] -= mu
        # CPA on real axis for spectrum
        wlist = np.linspace(Emin, Emax, Nw)
        zlist_real = wlist + 1j*delta
        print(f"CPA: Norb={no}, Nk={Nk}, Nw={Nw}, x={x_cpa}",flush=True)
        sigma_cpa = flibs.solve_cpa(hamk, VA, VB, x_cpa, zlist_real, pp=0.5)
        print("CPA converged",flush=True)
        # --- spectrum along symmetry line ---
        klist_s, spa_length, xticks = plibs.mk_klist(k_sets, kmesh, bvec)
        hamk_s = flibs.gen_ham(klist_s, ham_r, rvec)
        hamk_s[:, range(no), range(no)] -= mu
        Nks = len(klist_s)
        Akw = np.zeros((Nks, Nw))
        # A[iw,ik] = diag(zlist_real[iw]) - hamk_s[ik] - sigma_cpa[iw]
        A_batch = (np.eye(no)[None, None, :, :] * zlist_real[:, None, None, None]
                   - hamk_s[None, :, :, :]
                   - sigma_cpa[:, None, :, :])  # (Nw, Nks, no, no)
        try:
            Gk_all = np.linalg.inv(A_batch)
            Akw = (-np.trace(Gk_all, axis1=2, axis2=3).imag / np.pi).T
        except np.linalg.LinAlgError:
            print("Warning: Batch matrix inversion failed, falling back to element-wise.", flush=True)
            for iw in range(Nw):
                zI = np.eye(no) * zlist_real[iw]
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
if __name__=="__main__":
    main()
__license__="MIT"
