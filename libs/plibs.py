#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
@file plibs.py
@package plibs
@brief python script for model calculations
"""
import libs.flibs as flibs
import numpy as np, scipy.optimize as scopt, scipy.linalg as sclin

def import_hoppings(fname:str,ftype:int):
    """
    @fn import_hoppings()
    @brief this function import hopping parameters from files
    @param fname: Name of import files
    @param ftype: File format of import files
    @retval  rvec: the array of r-vector, size:(nr,3)
    @retval ham_r: the array of hopping integrals, size:(nr,no,no)
    @retval    no: the number of orbitals
    @retval    nr: the number of r-vectors
    """
    def import_hop(name: str):
        rvec=np.loadtxt(f'{name}/irvec.txt')
        nr=rvec[:,0].size
        ndegen=np.loadtxt(f'{name}/ndegen.txt')
        tmp=np.array([complex(float(tp[0]),float(tp[1])) for tp in
                      [f.strip(' ()\n').split(',') for f in open(f'{name}/ham_r.txt','r')]])
        no=int(np.sqrt(tmp.size/nr))
        ham_r=(tmp.reshape(nr,no,no).T/ndegen).T
        return(rvec,ham_r,no,nr)

    def import_out(name:str):
        data=np.loadtxt(name)
        con=(data[:,:3]==data[0,:3]).prod(axis=1).sum()
        no,nr =int(np.sqrt(con)),data[:,0].size//con
        rvec=np.array(data[:nr,:3])
        ham_r=(data[:,3]+1j*data[:,4]).reshape(no*no,nr).T.reshape(nr,no,no).round(6).copy()
        return(rvec,ham_r,no,nr)

    def import_hr(name:str):
        tmp=[f.split() for f in open(f'{name}_hr.dat','r')]
        no, nr=int(tmp[1][0]), int(tmp[2][0])
        c2,tmp1=3,[]
        while not len(tmp1)==nr:
            tmp1.extend(tmp[c2])
            c2=c2+1
        ndegen=np.array([float(t) for t in tmp1])
        tmp1=[[float(t) for t in tp] for tp in tmp[c2:]]
        tmp=np.array([complex(tp[5],tp[6]) for tp in tmp1])
        rvec=np.array([tmp1[no*no*i][:3] for i in range(nr)])
        ham_r=(tmp.reshape((nr,no,no)).T/ndegen).T.round(6).copy()
        return(rvec,ham_r,no,nr)

    def import_Hopping(name:str,sw_axis=False):
        tmp=[f.split() for f in open(f'{name}/Hopping.dat','r')]
        axis=np.array([[float(tp) for tp in tpp] for tpp in tmp[1:4]])
        no,nr=int(tmp[4][0]),int(tmp[4][1])
        tmp1=np.array([[float(t) for t in tp] for tp in tmp[7+no:]])
        rvec=np.array([tmp1[no*no*i][:3] for i in range(nr)])
        tmp=np.array([complex(tp[8],tp[9]) for tp in tmp1])
        ham_r=tmp.reshape(nr,no,no)
        if sw_axis:
            return(rvec,ham_r,no,nr,axis)
        else:
            return(rvec,ham_r,no,nr)
    if ftype==0:
        rvec,ham_r,no,nr=import_hop(fname)
    elif ftype==1:
        rvec,ham_r,no,nr=import_out(fname)
    elif ftype==2:
        rvec,ham_r,no,nr=import_hr(fname)
    else:
        rvec,ham_r,no,nr=import_Hopping(fname)
    return(rvec,ham_r,no,nr)

def import_MLO_hoppings(name:str):
    """
    @fn import_MLO_hoppings()
    @brief this function import MLO hopping parameters from files
    @param name: File name
    @retval  rvec: the array of r-vector, size:(nr,3)
    @retval ham_r: the array of hopping integrals, size:(nr,no,no)
    @retval   S_r: the array of overlap integrals, size:(nr,no,no)
    @retval    no: the number of orbitals
    @retval    nr: the number of r-vectors
    """
    tmp=[f.split() for f in open(f'{name}','r')]
    tmp1=np.array([[float(t) for t in tp] for tp in tmp])
    no=int(tmp1[:,0].max())
    nr=int(len(tmp1)/(no*no))
    tmp=np.array([complex(tp[5],tp[6]) for tp in tmp1])
    tmpS=np.array([complex(tp[7],tp[8]) for tp in tmp1])
    rvec=np.array([tmp1[i][2:5] for i in range(nr)])
    ham_r=tmp.reshape((no*no,nr)).T.reshape((nr,no,no)).round(6).copy()*13.6
    S_r=tmpS.reshape((no*no,nr)).T.reshape((nr,no,no)).round(6).copy()*13.6
    return rvec,ham_r,S_r,no,nr

def get_bvec(avec):
    """
    @fn get_bvec()
    @brief THis function generate reciprocal lattice vector from primitive translation vector
    @param  avec: primitive translation vector
    @return bvec: reciprocal lattice vector
    """
    bvec=2*np.pi*sclin.inv(avec).T
    return bvec

def get_eigs(klist,ham_r,S_r,rvec,sw_uni=False,sw_std=False):
    """
    @fn get_eigs
    @brief This function generate eigenvalues of Hamiltonian
    @param  klist: list of k-points
    @param  ham_r: hopping parameters
    @param    S_r: overlap integrals
    @param   rvec: r vector of hoppings
    @param sw_uni: switch of output only unitary matrix or not
    @param sw_std: switch standardization of unitary matrix at MLO hoppings
    @retval   eig: eigenvalues of Hamiltonian
    @retval   uni: unitary matrix of eigenfunctions
    """
    if len(S_r)==0:
        ham_k=flibs.gen_ham(klist,ham_r,rvec)
        eig,uni=flibs.get_eig(ham_k)
        if sw_uni:
            return uni
        else:
            return eig,uni
    else:
        ham_k,S_k=flibs.gen_ham(klist,ham_r,rvec,Ovl_r=S_r)
        eig,uni=flibs.get_eig(ham_k,S_k)
        if sw_std:
            norm=np.sqrt((abs(uni)**2).sum(axis=2))
            if sw_uni:
                return np.array([[u/nm for u,nm in zip(un,nor)] for un,nor in zip(uni,norm)])
            else:
                return eig,np.array([[u/nor for u,nm in zip(un,nor)] for un,nor in zip(uni,norm)])
        else:
            if sw_uni:
                return uni
            else:
                return eig,uni

def calc_mu(eig,Nk,fill:float,temp:float)-> float:
    """
    @fn calc_mu()
    @brief This function obtains chemical potential mu
    @param   eig: Eigenvales array
    @param    Nk: Number of k-point
    @param  fill: band filling
    @param  temp: Temperature
    @return   mu: chemical potential
    """
    no=int(eig.size/len(eig))
    def func(mu):
        sum_fermi=flibs.get_ffermi(eig,mu,temp).sum()
        return(fill*Nk-sum_fermi)
    emax=eig.max()
    emin=eig.min()
    mu=scopt.brentq(func,emin,emax)
    return mu

def calc_mu_imp(eigs,Nsite,fill:float,temp:float)-> float:
    itemp=1./temp
    def func(mu):
        return(fill*Nsite-0.5*(1.0-np.tanh(0.5*(eigs-mu)*itemp)).sum())
    emax=eigs.max()
    emin=eigs.min()
    mu=scopt.brentq(func,emin,emax)
    return mu

def gen_rlist(Nx:int,Ny:int,Nz:int):
    """
    @fn gen_rlist
    @brief This function generate the list of r-vector
    @param     Nx: Number of site of x-axis
    @param     Ny: Number of site of y-axis
    @param     Nz: Number of site of z-axis
    @return rlist: The array of r-vectors
    """
    x0=np.linspace(0,Nx,Nx,False)
    y0=np.linspace(0,Ny,Ny,False)
    z0=np.linspace(0,Nz,Nz,False)
    x,y,z=np.meshgrid(x0,y0,z0)
    rlist=np.array([x.ravel(),y.ravel(),z.ravel()]).T.copy()
    
    return rlist

def gen_klist_with_kmap(Nx:int,Ny:int,Nz:int):
    """
    @fn gen_klist_with_kmap
    @brief This function generate the list of r-vector
    @param     Nx: Number of kx mesh
    @param     Ny: Number of ky mesh
    @param     Nz: Number of kz mesh
    @retval klist: The array of k-vectors
    @retval  kmap: The array of property of k-points
    """
    x0=np.linspace(0,Nx,Nx,False,dtype=int)
    y0=np.linspace(0,Ny,Ny,False,dtype=int)
    z0=np.linspace(0,Nz,Nz,False,dtype=int)
    x,y,z=np.meshgrid(x0,y0,z0)
    kmap=np.array([x.ravel(),y.ravel(),z.ravel()]).T.copy()
    klist=np.array([x.ravel()/Nx,y.ravel()/Ny,z.ravel()/Nz]).T.copy()
    return klist,kmap

def gen_klist(Nx:int,Ny:int,Nz=None,sw_pp=True,kz=0):
    """
    @fn gen_klist()
    @brief This function generate k-point list messhed Nx,Ny,Nz
    @param     Nx: Number of axis 1 mesh (usually kx mesh)
    @param     Ny: Number of axis 2 mesh (usually ky mesh)
    @param     Nx: Number of axis 3 mesh (usually kz mesh)
    @param  sw_pp: switch output 2D mesh or 3D mesh
    @param     kz: value of axis 3 at 2D mesh (use only if sw_pp is True)
    @retval    Nk: number of k-points
    @retval klist: list of k-points
    """
    if sw_pp:
        kx=np.linspace(-0.5,0.5,Nx,True)
        ky=np.linspace(-0.5,0.5,Ny,True)
        if Nz==None:
            kz=np.array([kz])
        else:
            kz=np.linspace(-0.5,0.5,Nz,True)
    else:
        kx=np.linspace(0,1,Nx,False)
        ky=np.linspace(0,1,Ny,False)
        kz=np.linspace(0,1,Nz,False)
    x,y,z=np.meshgrid(kx,ky,kz)
    klist=np.array([x.ravel(),y.ravel(),z.ravel()]).T.copy()
    return len(klist),klist

def mk_klist(k_list,N:int,bvec):
    """
    @fn mk_klist()
    @brief This function generate k-point list for symmetry line
    @param  k_list: the list of symmetry points
    @param       N: Number of mesh between symmetry points
    @param    bvec: reciprocal lattice vector
    @retval  klist: list of k-points at band plot
    @retval  splen: length of symmetry points
    @retval xticks: footnotes of klist at ticks of symmetry points
    """
    klist=[]
    splen=[]
    xticks=[]
    maxsplen=0
    for ks,ke in zip(k_list,k_list[1:]):
        dkv=np.array(ke)-np.array(ks)
        dkv_length=abs(dkv.dot(bvec)).sum()
        tmp=np.linspace(ks,ke,N,False)
        tmp2=np.linspace(0,dkv_length,N,False)+maxsplen
        maxsplen=tmp2.max()
        xticks+=[tmp2[0]]
        klist+=tmp.tolist()
        splen+=tmp2.tolist()
    klist+=[k_list[-1]]
    splen+=[maxsplen+dkv_length/N]
    xticks+=[splen[-1]]
    return np.array(klist),np.array(splen),xticks

def mk_qlist(k_set,Nx:int,Ny:int,Nz:int,bvec):
    """
    @fn mk_qlist()
    @brief This function generate q-point list for symmetry line
    @param   k_set: the list of symmetry points
    @param      Nx: Number of x-mesh
    @param      Ny: Number of y-mesh
    @param      Nz: Number of z-mesh
    @param    bvec: reciprocal lattice vector
    @retval  qlist: list of q-points at band plot
    @retval  splen: length of symmetry points
    @retval xticks: footnotes of qlist at ticks of symmetry points
    """
    qlist=[]
    splen=[]
    xticks=[]
    maxsplen=0
    Narray=np.array([Nx,Ny,Nz])
    for ks,ke in zip(k_set,k_set[1:]):
        dk=np.array(ke)-np.array(ks)
        dk_length=abs(dk.dot(bvec)).sum()
        dN=np.asarray(abs(dk)*Narray,dtype=int)
        N=dN[dN>0].min()
        tmp=np.linspace(ks,ke,N,False)
        tmp2=np.linspace(0,dk_length,N,False)+maxsplen
        maxsplen+=dk_length
        xticks+=[tmp2[0]]
        qlist+=tmp.tolist()
        splen+=tmp2.tolist()
    qlist+=[k_set[-1]]
    splen+=[maxsplen+dk_length/N]
    xticks+=[splen[-1]]
    return np.array(qlist),np.array(splen),xticks

def mk_kf(mesh,rvec,ham_r,S_r,RotMat,mu:float,kz):
    import skimage.measure as sk
    Nk,klist=gen_klist(mesh+1,mesh+1,kz=kz)
    rvec1=RotMat.dot(rvec.T).T.copy()
    eig,uni=get_eigs(klist,ham_r,S_r,rvec1)
    v2=[]
    fsband=[]
    for i,e in enumerate(eig.T-mu):
        if(e.max()*e.min() < 0. ):
            cont=sk.find_contours(e.reshape(mesh+1,mesh+1),0)
            ct=[np.array([list(c)+[kz] for c in (cc-mesh/2)/mesh]) for cc in cont]
            fsband.append(i)
            v2.append(ct)
    return v2,fsband

def gen_3d_surf_points(mesh,rvec,ham_r,S_r,mu,kscale=1.0):
    import skimage.measure as ski
    Nk,klist=gen_klist(mesh+1,mesh+1,mesh+1)
    klist=klist*kscale
    eig,uni=get_eigs(klist,ham_r,S_r,rvec)
    if isinstance(kscale,int):
        ks=kscale*np.array([1.,1.,1.])
    else:
        ks=np.array(kscale)
    fspolys=[]
    fscenters=[]
    fsband=[]
    for i,e in enumerate(eig.T-mu):
        if(e.max()*e.min()<0.):
            verts,faces, _, _=ski.marching_cubes(e.reshape(mesh+1,mesh+1,mesh+1),0,
                                                 spacing=(ks[0]*2*np.pi/mesh,ks[1]*2*np.pi/mesh,
                                                          ks[2]*2*np.pi/mesh))
            verts=verts-ks*np.pi
            fspolys.append(verts[faces])
            fscenters.append(verts[faces].mean(axis=1)*.5/np.pi)
            fsband.append(i)
    return fspolys,fscenters,fsband

def get_colors(klist,blist,mrot,rvec,ham_r,S_r,ol,color_option,sw_2d=False):
    def get_col(cl,ol):
        col=(np.abs(cl[:,ol])**2 if isinstance(ol,int)
             else (np.abs(cl[:,ol])**2).sum(axis=1)).round(4)
        return col
    if color_option==0:
        return []
    elif color_option==1: #orbital weight color
        if sw_2d:
            uni=[[get_eigs(k,ham_r,S_r,rvec,True,True)[:,b] for k in kk] for kk,b in zip(klist,blist)]
            clist=[[np.array([get_col(cl,ol[0]),get_col(cl,ol[1]),get_col(cl,ol[2])]).T for cl in clst]
                   for clst in uni]
        else:
            uni=[get_eigs(k,ham_r,S_r,rvec,True,True)[:,b] for k,b in zip(klist,blist)]
            clist=[np.array([get_col(clst,ol[0]),get_col(clst,ol[1]),get_col(clst,ol[2])]).T for clst in uni]
        return clist
    elif color_option==2: #velocity size color
        if sw_2d:
            uni=[[get_eigs(k,ham_r,S_r,rvec,True,False) for k in kk] for kk in klist]
            vk=[[flibs.get_veloc(k,ham_r,rvec,mrot,unkk)[:,b,:] for k,unkk in zip(kk,unk)]
                for kk,unk,b in zip(klist,uni,blist)]
            clist=[[np.sqrt((abs(vkkk)*abs(vkkk)).sum(axis=1)) for vkkk in vkk] for vkk in vk]
        else:
            uni=[get_eigs(k,ham_r,S_r,rvec,True,False) for k in klist]
            vk=[flibs.get_veloc(kk,ham_r,rvec,mrot,unk)[:,b,:] for kk,unk,b in zip(klist,uni,blist)]
            clist=[np.sqrt((abs(vkk)*abs(vkk)).sum(axis=1)) for vkk in vk]
        return clist

def get_emesh(Nx:int,Ny:int,Nz:int,ham_r,S_r,rvec,avec,sw_uni=False,sw_veloc=False):
    Nk,klist=gen_klist(Nx,Ny,Nz,sw_pp=False)
    eig,uni=get_eigs(klist,ham_r,S_r,rvec)
    kweight=np.ones(len(eig),dtype=np.float64)
    if sw_veloc:
        if sw_uni:
            vk=flibs.get_vnmk(klist,ham_r,rvec,avec,uni)
            return Nk,eig,vk,kweight
        else:
            vk=flibs.get_veloc(klist,ham_r,rvec,avec,uni)
            return Nk,eig,vk,kweight
    else:
        if sw_uni:
            return Nk,klist,eig,uni,kweight
        else:
            return Nk,eig,kweight

def get_ptv(alatt,deg,brav):
    if brav==0: #simple
        Arot=np.array([[ 1., 0., 0.],[ 0., 1., 0.],[ 0.,0., 1.]])
    elif brav==1: #face center
        Arot=np.array([[-.5, 0., .5],[0., .5, .5],[-.5,.5, 0.]])
    elif brav==2: #body center
        Arot=np.array([[.5,-.5, .5],[.5, .5, .5],[-.5,-.5, .5]])
    elif brav==3: #hexagonal
        Arot=np.array([[1., 0., 0.],[-.5,.5*np.sqrt(3.),0.],[0., 0., 1.]])
    elif brav==4: #trigonal
        cg=np.cos(np.pi*deg[2]/180.)
        tx,ty,tz=np.sqrt((1-cg)*0.5),np.sqrt((1.-cg)/6.),np.sqrt((1.+2*cg)/3.)
        Arot=np.array([[ tx,  -ty, tz],[ 0., 2*ty, tz],[-tx,  -ty, tz]])
    elif brav==5: #base center
        Arot=np.array([[ .5, .5, 0.],[-.5, .5, 0.],[ 0., 0., 1.]])
    elif brav==6: #face center 2
        Arot=np.array([[.5, 0., .5],[0., .5, .5],[.5,.5, 0.]])
    elif brav==7: #body center 2
        Arot=np.array([[-.5, .5, .5],[.5, -.5, .5],[.5, .5, -.5]])
    else:
        phase=[np.pi*deg[0]/180.,np.pi*deg[1]/180.,np.pi*deg[2]/180.]
        ax1,a2=alatt[1]/alatt[0],alatt[2]/alatt[0]
        r1,r2,r3=ax1*np.cos(phase[2]),ax1*np.sin(phase[2]),ax2*np.cos(phase[1])
        r4=ax2*(np.cos(phase[0])-np.cos(phase[1])*np.cos(phase[2]))/np.sin(phase[2])
        r5=ax2*np.sqrt(1+2*np.cos(phase[0])*np.cos(phase[1])*np.cos(phase[2])
                       -(np.cos(phase[0])**2+np.cos(phase[1])**2+np.cos(phase[2])**2))/np.sin(phase[2])
        Arot=np.array([[ 1., 0.,  0.],
                      [ r1, r2,  0.],
                      [ r3, r4, r5]])
    if brav in {0,1,2,3,4,5,6,7}:
        avec=alatt*Arot
    else:
        avec=alatt[0]*Arot
    return avec,Arot

def get_symm_line(brav:int)->tuple[list,list]:
    if brav==0: #simple
        k_list=[[0.,0.,.5],[0., 0., 0.],[.5, 0., 0.],[.5, .5, 0.],[0.,0.,0.]]
        xlabel=['Z',r'$\Gamma$','X','M',r'$\Gamma$']
    elif brav==1: #face center
        k_list=[[0.,0.,0.],[.5, 0., .5],[1., 0., 0.],[.5, .5, .5],[.5,.25,.75],[0.,0.,0.]]
        xlabel=[r'$\Gamma$','X',r'$\Gamma$','L','W',r'$\Gamma$']
    elif brav==2: #body center
        k_list=[[.5,.5,.5],[0., 0., 0.],[.5, 0., 0.],[.5, .5,-.5],[0.,0.,0.]]
        xlabel=['Z',r'$\Gamma$','X','M',r'$\Gamma$']
    elif brav==3: #hexagonal
        k_list=[[0.,0.,0.],[2./3.,-1./3., 0.],[.5, 0., 0.],[0., 0., 0.],[0.,0.,.5]]
        xlabel=[r'$\Gamma$','K','M',r'$\Gamma$','Z']
    elif brav==4: #trigonal
        k_list=[[0.,0.,0.],[.5,0.,.5],[.5,0.,0.],[0.,0.,0.],[.5,.5,.5]]
        xlabel=[r'$\Gamma$','K','M',r'$\Gamma$','Z']
    elif brav==6:
        k_list=[[0.,0.,0.],[0., .5, .5],[.5, .5, .5],[.25,.75,.5],[0.,0.,0.]]
        xlabel=[r'$\Gamma$','X','L','W',r'$\Gamma$']
    else:
        pass
    return k_list,xlabel

def BZedge(brav:int):
    pass

def get_conductivity(mu,temp,eig,vk,Nw,Emax,idelta=1.e-3):
    wlist=np.linspace(0,Emax,Nw)
    ffermi=flibs.get_ffermi(eig,mu,temp)
    L11=[]
    L12=[]
    L22=[]
    for w in wlist:
        L11w,L12w,L22w=flibs.calc_Lij(eig,vk,ffermi,mu,w,idelta,temp)
        L11.append(L11w)
        L12.append(L12w)
        L22.append(L22w)
    L11=np.array(L11)
    L12=np.array(L12)
    L22=np.array(L22)
    return L11,L12,L22,wlist

def chis_spectrum(mu,temp,Smat,klist,qlist,olist,eig,uni,Nw,Emax,idelta=1.e-3):
    ffermi=flibs.get_ffermi(eig,mu,temp)
    wlist=np.linspace(0,Emax,Nw)
    chisq=[]
    chis_orbq=[]
    f=open('chi0.dat','w')
    fq=open('writeq.dat','w')
    for i,q in enumerate(qlist):
        fq.write(f'{i:d} {q[0]:5.3f} {q[1]:5.3f} {q[2]:5.3f}\n')
        fq.flush()
        qshift=flibs.get_qshift(klist,q)
        chi0=flibs.get_chi_irr(uni,eig,ffermi,qshift,olist,wlist,idelta,temp)
        chis=flibs.get_chis(chi0,Smat)
        trchis,trchi0,chis_orb=flibs.get_tr_chi(chis,chi0,olist)
        chisq.append(trchis)
        chis_orbq.append(chis_orb)
        for w,trchi in zip(wlist,trchi0):
            f.write(f'{i:8.4f} {w:8.4f} {trchi.imag:9.4f} {trchi.real:9.4f}\n')
        f.write('\n')
        f.flush()
    f.close()
    fq.close()
    return np.array(chisq),np.array(chis_orbq),wlist

def chis_q_point(q,eig,uni,Emax,Nw,mu,temp,Smat,klist,olist,idelta):
    ffermi=flibs.get_ffermi(eig,mu,temp)
    wlist=np.linspace(0,Emax,Nw)
    qshift=flibs.get_qshift(klist,q)
    chi0=flibs.get_chi_irr(uni,eig,ffermi,qshift,olist,wlist,idelta,temp,0)
    chis=flibs.get_chis(chi0,Smat)
    trchis,trchi0,chis_orb=flibs.get_tr_chi(chis,chi0,olist)
    return trchis,chis_orb,wlist

def chis_qmap(Nx,Ny,Ecut,mu,temp,Smat,klist,olist,eig,uni,idelta=1.e-3):
    ffermi=flibs.get_ffermi(eig,mu,temp)
    chis,chi0=flibs.chis_qmap(uni,eig,ffermi,klist,Smat,olist,Nx,Ny,temp,Ecut,idelta)
    x0=np.linspace(0,1,Nx,False)
    y0=np.linspace(0,1,Ny,False)
    qx,qy=np.meshgrid(x0,y0)
    return chis,chi0,qx,qy

def phi_spectrum(mu,temp,klist,qlist,olist,eig,uni,Nw,Emax,idelta=1.e-3):
    ffermi=flibs.get_ffermi(eig,mu,temp)
    wlist=np.linspace(0,Emax,Nw)
    phiq=[]
    phi_orbq=[]
    fq=open('writeq.dat','w')
    for i,q in enumerate(qlist):
        fq.write(f'{i:d} {q[0]:5.3f} {q[1]:5.3f} {q[2]:5.3f}\n')
        fq.flush()
        qshift=flibs.get_iqshift(klist,q)
        phi=flibs.get_phi_irr(uni,eig,ffermi,qshift,olist,wlist,idelta,mu,temp)
        trphi,phi_orb=flibs.get_tr_phi(phi,olist)
        phiq.append(trphi)
        phi_orbq.append(phi_orb)
    fq.close()
    return np.array(phiq),np.array(phi_orbq),wlist

def phi_qmap(Nx,Ny,Ecut,mu,temp,klist,olist,eig,uni,idelta=1.e-3,sw_omega=True):
    ffermi=flibs.get_ffermi(eig,mu,temp)
    phi=flibs.phi_qmap(uni,eig,ffermi,klist,olist,Nx,Ny,mu,temp,Ecut,idelta,sw_omega)
    x0=np.linspace(0,1,Nx,False)
    y0=np.linspace(0,1,Ny,False)
    qx,qy=np.meshgrid(x0,y0)
    return phi,qx,qy

def get_chi_orb_list(Norb,site_prof):
    if(len(site_prof)==1):
        tmp=np.arange(Norb)+1
        o1,o2=np.meshgrid(tmp,tmp)
        chiolist=np.array([o1.flatten(),o2.flatten()]).T
    else:
        if(Norb==sum(site_prof)):
            chiolist=[]
            N0=1
            for i_site in site_prof:
                tmp=np.arange(i_site)+N0
                o1,o2=np.meshgrid(tmp,tmp)
                chiolist+=list(np.array([o1.flatten(),o2.flatten()]).T)
                N0+=i_site
            chiolist=np.array(chiolist)
        else:
            print("site_prof doesn't correspond to Hamiltonian")
            exit()
    return chiolist

def get_initial_gap(kmap,klist,Norb,gap_sym):
    if gap_sym>=0:
        gapsym=['s','dx2-y2','spm','dxy']
        print('gap symmetry is '+gapsym[gap_sym])
    else:
        gapsym=['s','px','py']
        print('gap symmetry is '+gapsym[-gap_sym])
    if gap_sym==0: #s
        init_gap=np.ones((Norb,len(kmap)),dtype=np.float64)
    else:
        A=2*np.pi/(kmap[:,0].max()+1)
        B=2*np.pi/(kmap[:,1].max()+1)
        init_gap=np.zeros((Norb,len(kmap)),dtype=np.float64)
        for i in range(Norb):
            if gap_sym==1: #dx2-y2
                init_gap[i,:]=np.cos(A*kmap[:,0])-np.cos(B*kmap[:,1])
            elif gap_sym==2: #spm
                init_gap[i,:]=2*np.cos(A*kmap[:,0])*np.cos(B*kmap[:,1])
            elif gap_sym==3: #dxy
                init_gap[i,:]=2*np.sin(A*kmap[:,0])*np.sin(B*kmap[:,1])
            elif gap_sym==-1: #px
                init_gap[i,:]=2*np.sin(A*kmap[:,0])
            elif gap_sym==-2: #py
                init_gap[i,:]=2*np.sin(B*kmap[:,1])
    return init_gap
