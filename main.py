#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
ftype: set input hamiltonian's format
0: ham_r.txt, irvec.txt, ndegen.txt in {fname} dir
1: .input file named {fname}
2: {fname}_hr.dat file (wannier90 default hopping file)
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
#fname,ftype,brav='Si',2,1
#fname,ftype,brav='00010.input',1,2
fname,ftype,brav='hop3.input',1,0

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
 8: calc chis on xy plane at Ecut
 9: calc carrier num.
10: calc cycrotron mass
11: calc selfenergy using flex
12: mass calculation
color_option defines the meaning of color on Fermi surfaces
 0: band or mono color
 1: orbital weight settled by olist
 2: velocity size
"""
option=0
color_option=1

Nx,Ny,Nz,Nw=40,40,1,500 #k and energy(or matsubara freq.) mesh size
kmesh=200               #kmesh for spaghetti plot
#kscale=[1.0,1.0,1.0]
#kz=0.5

abc=[3.96,3.96,13.02]
alpha_beta_gamma=[90.,90.,90]
temp=2.59e-3
fill=2.0

Emin,Emax=-10,5
delta=1.0e-1
Ecut=1.0e-3
tau_const=100
olist=[[0],[1],[2]]
U,J=0.8, 0.1 #1.2,0.15

#k_sets=[[0., 0., 0.],[.5, 0., 0.],[.5, .5, 0.]]
#xlabel=['$\Gamma$','X','M']

sw_calc_mu=True #calculate mu or not
sw_unit=True    #set unit values unity (False) or not (True)
sw_tdf=False
#----------------------------------main functions-------------------------------------
import numpy as np
import flibs, plibs
import scipy.linalg as sclin
import scipy.constants as scconst
import matplotlib.pyplot as plt
import matplotlib.cm as cm

alatt=np.array(abc)
deg=np.array(alpha_beta_gamma)

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
def plot_band(eig,spl,xlabel,xticks,uni,ol,color):
    def get_col(cl,ol):
        col=(np.abs(cl[ol])**2 if isinstance(ol,int)
             else (np.abs(cl[ol])**2).sum(axis=0)).round(4)
        return col
    fig=plt.figure()
    ax=plt.axes()
    for e,cl in zip(eig,uni):
        if color:
            c1=get_col(cl,ol[0])
            c2=get_col(cl,ol[1])
            c3=get_col(cl,ol[2])
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

def plot_FS(fscolors,klist,color_option):
    fig=plt.figure()
    ax=fig.add_subplot(111,aspect='equal')
    col=['r','g','b','c','m','y','k','w']
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
    for kl,fscol,cb in zip(klist,fscolors,col):
        for kk,fcol,in zip(kl,fscol):
            if color_option==0:
                for k1,k2 in zip(kk,kk[1:]):
                    plt.plot([k1[0],k2[0]],[k1[1],k2[1]],color='black')
            else:
                if color_option==2:
                    clist=cm.jet((fcol-vmin)/(vmax-vmin))
                else:
                    clist=fcol
                for k1,k2,clst in zip(kk,kk[1:],clist):
                    plt.plot([k1[0],k2[0]],[k1[1],k2[1]],c=clst,lw=3)
    if color_option==2:
        plt.scatter(k[:,0],k[:,1],s=0.1,c=v)
        plt.jet()
        plt.colorbar()
    plt.xlim(-0.5,0.5)
    plt.ylim(-0.5,0.5)
    plt.xticks([-0.5,0,0.5],['-$\pi$','0','$\pi$'])
    plt.yticks([-0.5,0,0.5],['-$\pi$','0','$\pi$'])
    plt.show()

def plot_3d_surf(fspolys,fscenters,fscolors,surface_opt):
    from mpl_toolkits.mplot3d import axes3d
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    import matplotlib.colors as colors
    fig=plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    clist=['r','g','b','c','m','y','k','w']
    if surface_opt==0:
        for j,polys in enumerate(fspolys):
            ax.add_collection3d(Poly3DCollection(polys,facecolor=clist[j%6]))
    else:
        if surface_opt==1:
            tri=Poly3DCollection(fspolys,facecolors=fscolors,lw=0)
        else:
            clmax=fscolors.max()
            clmin=fscolors.min()
            nor_cols=(fscolors-clmin)/(clmax-clmin)
            tri=Poly3DCollection(fspolys,facecolors=cm.jet(nor_cols),lw=0)
            fs=ax.scatter(fscenters[:,0],fscenters[:,1],fscenters[:,2],
                          c=fscolors,cmap=cm.jet,s=0.1)
            plt.colorbar(fs,format='%.2e')
        ax.add_collection3d(tri)
    ax.grid(False)
    ax.set_xlim(-np.pi, np.pi)
    ax.set_ylim(-np.pi, np.pi)
    ax.set_zlim(-np.pi, np.pi)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    plt.tight_layout()
    plibs.BZedge(0)
    plt.show()

def set_init_3dfsplot(color_option,polys,centers,blist,avec,rvec,ham_r,olist):
    if color_option==0:
        fspolys=polys
        fscenters=[]
        fscolors=[]
    else:
        colors=plibs.get_colors(centers,blist,ihbar*avec.T,rvec,ham_r,olist,color_option)
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

def plot_spectrum(k_sets,xlabel,kmesh,bvec,mu,ham_r,rvec,Emin,Emax,delta,Nw):
    klist,spa_length,xticks=plibs.mk_klist(k_sets,kmesh,bvec)
    ham_k=flibs.gen_ham(klist,ham_r,rvec)
    eig,uni=flibs.get_eig(ham_k)
    wlist=np.linspace(Emin,Emax,Nw)
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
    
def calc_conductivity_Bolzmann(rvec,ham_r,avec,fill,temp,tau_const,Nw=300,with_spin=False):
    #no dep. T and mu or filling
    Nk,eig,vk,kweight=plibs.get_emesh(Nx,Ny,Nz,ham_r,rvec,avec.T*ihbar,sw_veloc=True)
    Vuc=sclin.det(avec)*1e-30
    gsp=(1.0 if with_spin else 2.0) #spin weight
    iNV=1./(Nk*Vuc)
    tauconst=True
    tau=eig*0.+tau_const

    itemp=1./temp
    mu=plibs.calc_mu(eig,Nk,fill,temp)
    print("T = %6.3f K"%(temp/kb))
    print("mu = %6.3f eV"%mu)
    if tauconst:
        print("tau = %d"%tau_const+('fs' if sw_unit else ''))
    else:
        print("max tau = %d"%tau.max()+('fs' if sw_unit else ''))
    if sw_tdf:
        tdf=flibs.calc_tdf(eig,vk,kweight,tau,Nw)
        wlist=np.linspace(eig.min(),eig.max(),Nw)
        dw=(eig.max()-eig.min())/Nw
        dfermi=0.25*(1.-np.tanh(0.5*(wlist-mu)/temp)**2)/temp
        K0=(dfermi*tdf.T).T.sum(axis=0)*dw
        K1=(dfermi*(wlist-mu)*tdf.T).T.sum(axis=0)*dw
        K2=(dfermi*(wlist-mu)**2*tdf.T).T.sum(axis=0)*dw
        plt.plot(wlist-mu,tdf[:,0,0])
        plt.show()
    else:
        K0,K1,K2=flibs.calc_Kn(eig,vk,kweight,temp,mu,tau)
    sigma=gsp*tau_unit*eC*K0*iNV
    kappa=gsp*tau_unit*eC*kb*K2*iNV*itemp
    kappa2=gsp*tau_unit*eC*kb*(K2-K1.dot(sclin.inv(K0).dot(K1)))*iNV*itemp
    Seebeck=-kb*sclin.inv(K0).dot(K1)*itemp

    Lorenz=kb*kappa*sclin.inv(sigma*temp)
    Lorenz2=kb*kappa2*sclin.inv(sigma*temp)
    Pertier=K1.dot(sclin.inv(K0))
    sigmaS=gsp*tau_unit*kb*eC*K1*iNV*itemp
    PF=sigma*Seebeck**2
    print('sigma matrix (S/m)')
    print(sigma.round(10))
    print('kappa matrix (K22 only) (W/m/K')
    print(kappa.round(10))
    print('kappa matrix (full) (W/m/K)')
    print(kappa2.round(10))
    print('sigmaS matrix (A/m/K)')
    print(sigmaS.round(10))
    print('Seebeck matrix (V/K)')
    print(Seebeck.round(10))
    print('Pertier matrix (V)')
    print(Pertier.round(13))
    print('Lorenz matrix (K22 only) (Wohm/K^2)')
    print(Lorenz.round(10))
    print('Lorenz matrix? (full) (Wohm/K^2)')
    print(Lorenz2.round(10))
    print('Power Factor (SA/m^2/K)')
    print(PF.round(10))

def calc_conductivity_lr(rvec,ham_r,avec,Nx,Ny,Nz,fill,temp,Nw,delta,with_spin=False):
    '''
    calculation of linear response theory
    electric conductivity of LRT correponds to Boltzmann then delta~O(10-1) (tau~1fs) at 300K
    thermal conductivity and L12 of LRT correspond to Boltzmann then delta~O(10-3) (tau~100fs) at 300K
    '''
    Nk,eig,vk,kweight=plibs.get_emesh(Nx,Ny,Nz,ham_r,rvec,avec.T*ihbar,True,True)
    Vuc=sclin.det(avec)
    gsp=(1.0 if with_spin else 2.0) #spin weight
    mu=plibs.calc_mu(eig,Nk,fill,temp)
    #delta=hbar*1.e15/100
    print('chemical potential = %6.3f'%mu)
    print('tempreture = %6.3f'%(temp/kb))
    print('about tau = ',(hbar/delta))
    print('delta = %9.3e'%delta)
    L11,L12,L22,wlist=plibs.get_conductivity(mu,temp,eig,vk,Nw,Emax,delta)
    sigmaconst=gsp*hbar*eC/Vuc*1.0e30
    kappaSconst=sigmaconst*kb/temp
    sigma=sigmaconst*L11
    kappa=kappaSconst*L22
    sigmaS=kappaSconst*L12
    Seebeck=np.array([-sclin.inv(s).dot(sS) for s,sS in zip(sigma,sigmaS)])
    print('sigma matrix (S/m)')
    print(sigma[0].real.round(10))
    print('kappa matrix (L22 only) (W/m/K)')
    print(kappa[0].real.round(10))
    print('sigmaS matrix (A/m/K)')
    print(sigmaS[0].real.round(10))
    print('Lorenz number (Wohm/K^2)')
    print((kb*kappa[0]*sclin.inv(sigma[0]*temp)).real.round(10))
    print('Seebeck coefficient matrix (V/K)')
    print(Seebeck[0].real.round(10))
    print(hbar)
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

def calc_flex(Nx,Ny,Nz,Nw,ham_r,rvec,mu,temp,olist):
    klist,kmap=plibs.gen_klist_with_kmap(Nx,Ny,Nz)
    ham_k=flibs.gen_ham(klist,ham_r,rvec)
    eig,uni=flibs.get_eig(ham_k)
    Gk=flibs.gen_Green0(eig,uni,mu,temp,Nw)
    chi=flibs.get_chi0_comb(Gk,kmap,olist,Nx,Ny,Nz,Nw)
    print(chi)

def get_carrier_num(kmesh,rvec,ham_r,mu):
    Nk,eig,kwieght=plibs.get_emesh(kmesh,kmesh,kmesh,ham_r,rvec,Arot)
    for i,en in enumerate(eig.T-mu):
        num_hole=float(np.where(en>0)[0].size)/Nk
        num_particle=float(np.where(en<=0)[0].size)/Nk
        print(i+1,round(num_hole,4),round(num_particle,4))

def get_mu(ham_r,rvec,Arot,temp,kmesh=40):
    Nk,eig,kweight=plibs.get_emesh(kmesh,kmesh,kmesh,ham_r,rvec,Arot)
    mu=plibs.calc_mu(eig,Nk,fill,temp)
    return mu

def get_mass(mesh,rvec,ham_r,mu,de=3.e-4,meshkz=20):
    pass

def main():
    rvec,ham_r,no,Nr=plibs.import_hoppings(fname,ftype)
    avec,Arot=plibs.get_ptv(alatt,deg,brav)
    if sw_dec_axis:
        rvec1=Arot.T.dot(rvec.T).T
        rvec=rvec1
        if brav in {1,2}:
            rvec[:,2]*=2.
            avec=alatt*np.eye(3)
            avec[:,2]*=.5
        elif brav==5:
            rvec*=2.
            avec=(alatt*np.eye(3))*.5
        else:
            avec=alatt*np.eye(3)
    bvec=plibs.get_bvec(avec)
    print("Hamiltonian name is "+fname)
    print("Lattice Vector")
    print(avec)
    print("Reciprocal Lattice Vector")
    print(bvec)
    print("Number of orbital =",no)
    if option in {0}:
        print('kmesh = %d'%kmesh)
    else:
        print('k-mesh is %d %d %d'%(Nx,Ny,Nz))
    if option in {5,6}:
        pass
    else:
        if sw_calc_mu:
            mu=get_mu(ham_r,rvec,Arot,temp)
        else:
            try:
                mu=mu0
            except NameError:
                mu=get_mu(ham_r,rvec,Arot,temp)
        print('chem. pot. = %7.4f'%mu)

    if option in {0,4,7,12}:
        try:
            k_sets
            xlabel
        except NameError:
            print('generate symmetry line')
            k_sets,xlabel2=plibs.get_symm_line(brav)
            try:
                xlabel
                if len(xlabel)!=len(k_sets):
                    xlabel=xlabel2
            except NameError:
                xlabel=xlabel2
    if option==0: #band plot
        print("calculate band structure")
        klist,spa_length,xticks=plibs.mk_klist(k_sets,kmesh,bvec)
        ham_k=flibs.gen_ham(klist,ham_r,rvec)
        eig,uni0=flibs.get_eig(ham_k)
        uni=np.array([u.T for u in uni0]) #rotate uni(k,band,orb) to uni(k,orb,band)
        plot_band(eig.T-mu,spa_length,xlabel,xticks,uni.T,olist,(False if color_option==0 else True))
    elif option==1: #plot dos
        print("calculate Dos")
        Nk,eig,kweight=plibs.get_emesh(Nx,Ny,Nz,ham_r,rvec,avec)
        wlist=np.linspace(Emin,Emax,Nw,True)
        Dos=flibs.gen_tr_Greenw_0(eig,mu,wlist,delta).sum(axis=0)/Nk
        plt.plot(wlist,Dos,color='black')
        plt.show()
    elif option==2: #2D Fermi surface plot
        print("plot 2D Fermi surface")
        try:
            kz
        except NameError:
            kz=0.
        klist,blist=plibs.mk_kf(Nx,rvec,ham_r,mu,kz)
        clist=plibs.get_colors(klist,blist,ihbar*avec.T,rvec,ham_r,olist,color_option,True)
        plot_FS(clist,klist,color_option)
    elif option==3: #3D Fermi surface plot
        print("plot 3D Fermi surface")
        try:
            kscale
        except NameError:
            kscale=1.0
        polys,centers,blist=plibs.gen_3d_surf_points(Nx,rvec,ham_r,mu,kscale)
        fspolys,fscenters,fscolors=set_init_3dfsplot(color_option,polys,centers,blist,avec,rvec,ham_r,olist)
        plot_3d_surf(fspolys,fscenters,fscolors,color_option)
    elif option==4: #plot spectrum
        print("calculate spectrum")
        plot_spectrum(k_sets,xlabel,kmesh,bvec,mu,ham_r,rvec,Emin,Emax,delta,Nw)
    elif option==5: #calc conductivity
        print("calculate conductivities using Boltzmann theory")
        calc_conductivity_Bolzmann(rvec,ham_r,avec,fill,temp,tau_const)
    elif option==6: #calc_optical conductivity
        print("calculate conductivities")
        calc_conductivity_lr(rvec,ham_r,avec,Nx,Ny,Nz,fill,temp,Nw,delta)
    elif option in {7,8}: #calc_chis_spectrum
        print("calculate chis spectrum" if option==7 else "calculate chis qmap at Ecut")
        Nk,klist,eig,uni,kweight=plibs.get_emesh(Nx,Ny,Nz,ham_r,rvec,avec,sw_uni=True)
        try:
            chiolist
        except NameError:
            Norb=int(eig.size/Nk)
            tmp=np.arange(Norb)+1
            o1,o2=np.meshgrid(tmp,tmp)
            chiolist=np.array([o1.flatten(),o2.flatten()]).T
        Smat=flibs.gen_Smatrix(chiolist,U,J)
        if option==7: #chis spectrum with symmetry line
            qlist,spa_length,xticks=plibs.mk_qlist(k_sets,Nx,Ny,Nz,bvec)
            chisw,wlist=plibs.chis_spectrum(mu,temp,Smat,klist,qlist,chiolist,eig,uni,Nw,Emax,delta)
            w,sp=np.meshgrid(wlist,spa_length)
            f=open('chis.dat','w')
            for ww,ssp,chis in zip(w,sp,chisw):
                for www,sssp,chi in zip(ww,ssp,chis):
                    f.write(f'{www} {sssp} {chi.imag}\n')
                f.write('\n')
            f.close()
            plt.contourf(sp,w,abs(chisw.imag),100)
            plt.colorbar()
            #plt.jet()
            plt.hot()
            plt.show()
        else: #chis spectrum ecut plane
            chis,chi0,qx,qy=plibs.chis_qmap(Nx,Ny,Ecut,mu,temp,Smat,klist,chiolist,eig,uni,idelta=1.e-3)
            print('plot chis map')
            print(chis)
            plt.contourf(qx,qy,abs(chis.imag),100)
            plt.colorbar()
            plt.jet()
            plt.show()
            plt.contourf(qx,qy,abs(chi0.imag),100)
            plt.colorbar()
            plt.jet()
            plt.show()
    elif option==9: #calc carrier number
        print("calculate carrier number")
        get_carrier_num(Nx,rvec,ham_r,mu)
    elif option==10: #calc cycrtron mass
        print("calculate cycrtron mass")
        get_mass(Nx,rvec,ham_r,mu)
    elif option==11:
        print("calc self energy")
        Nk,klist,eig,uni,kweight=plibs.get_emesh(Nx,Ny,Nz,ham_r,rvec,avec,sw_uni=True)
        try:
            chiolist
        except NameError:
            Norb=int(eig.size/Nk)
            tmp=np.arange(Norb)+1
            o1,o2=np.meshgrid(tmp,tmp)
            chiolist=np.array([o1.flatten(),o2.flatten()]).T
        calc_flex(Nx,Ny,Nz,Nw,ham_r,rvec,mu,temp,chiolist)
    elif option==12: #mass calc
        klist,spa_length,xticks=plibs.mk_klist(k_sets,kmesh,bvec)
        ham_k=flibs.gen_ham(klist,ham_r,rvec)
        eig,uni=flibs.get_eig(ham_k)
        mass=flibs.get_mass(klist,ham_r,rvec,avec.T*ihbar,uni)*eC/emass
        print(mass[:,3,:,:])
        
if __name__=="__main__":
    main()
