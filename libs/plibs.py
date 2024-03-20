#!/usr/bin/env python
#-*- coding:utf-8 -*-
from __future__ import print_function, division
import libs.flibs as flibs
import numpy as np, scipy.optimize as scopt, scipy.linalg as sclin

def import_hoppings(fname,ftype):
    def import_hop(name):
        rvec=np.loadtxt(name+'/irvec.txt')
        nr=rvec[:,0].size
        ndegen=np.loadtxt(name+'/ndegen.txt')
        tmp=np.array([complex(float(tp[0]),float(tp[1])) for tp in
                      [f.strip(' ()\n').split(',') for f in open(fname+'/ham_r.txt','r')]])
        no=int(np.sqrt(tmp.size/nr))
        ham_r=(tmp.reshape(nr,no,no).T/ndegen).T
        return(rvec,ham_r,no,nr)

    def import_out(name):
        data=np.loadtxt(name)
        con=(data[:,:3]==data[0,:3]).prod(axis=1).sum()
        no,nr =int(np.sqrt(con)),data[:,0].size//con
        rvec=np.array(data[:nr,:3])
        ham_r=(data[:,3]+1j*data[:,4]).reshape(no*no,nr).T.reshape(nr,no,no).round(6).copy()
        return(rvec,ham_r,no,nr)

    def import_hr(name):
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

    def import_Hopping(name,sw_axis=False):
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

def get_bvec(avec):
    bvec=2*np.pi*sclin.inv(avec).T
    return bvec

def calc_mu(eig,Nk,fill,temp):
    no=int(eig.size/len(eig))
    def func(mu):
        sum_tanh=np.tanh(0.5*(eig-mu)/temp).sum()
        return((2*fill-no)*Nk+sum_tanh)
    emax=eig.max()
    emin=eig.min()
    mu=scopt.brentq(func,emin,emax)
    return mu

def gen_klist_with_kmap(Nx,Ny,Nz):
    x0=np.linspace(0,Nx,Nx,False,dtype=int)
    y0=np.linspace(0,Ny,Ny,False,dtype=int)
    z0=np.linspace(0,Nz,Nz,False,dtype=int)
    x,y,z=np.meshgrid(x0,y0,z0)
    kmap=np.array([x.ravel(),y.ravel(),z.ravel()]).T.copy()
    klist=np.array([x.ravel()/Nx,y.ravel()/Ny,z.ravel()/Nz]).T.copy()
    return klist,kmap

def gen_klist(Nx,Ny,Nz=None,sw_pp=True,kz=0):
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
    Nk=len(klist)
    return Nk,klist

def mk_klist(k_list,N,bvec):
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

def mk_qlist(k_set,Nx,Ny,Nz,bvec):
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

def mk_kf(mesh,rvec,ham_r,mu,kz):
    import skimage.measure as sk
    Nk,klist=gen_klist(mesh+1,mesh+1,kz=kz)
    ham_k=flibs.gen_ham(klist,ham_r,rvec)
    eig,uni=flibs.get_eig(ham_k)
    v2=[]
    fsband=[]
    for i,e in enumerate(eig.T-mu):
        if(e.max()*e.min() < 0. ):
            cont=sk.find_contours(e.reshape(mesh+1,mesh+1),0)
            ct=[np.array([list(c)+[kz] for c in (cc-mesh/2)/mesh]) for cc in cont]
            fsband.append(i)
            v2.append(ct)
    return v2,fsband

def gen_3d_surf_points(mesh,rvec,ham_r,mu,kscale=1.0):
    import skimage.measure as ski
    Nk,klist=gen_klist(mesh+1,mesh+1,mesh+1)
    klist=klist*kscale
    ham_k=flibs.gen_ham(klist,ham_r,rvec)
    eig,uni=flibs.get_eig(ham_k)
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
                                                spacing=(ks[0]*2*np.pi/mesh,ks[1]*2*np.pi/mesh,ks[2]*2*np.pi/mesh))
            verts=verts-ks*np.pi
            fspolys.append(verts[faces])
            fscenters.append(verts[faces].mean(axis=1)*.5/np.pi)
            fsband.append(i)
    return fspolys,fscenters,fsband

def get_colors(klist,blist,mrot,rvec,ham_r,ol,color_option,sw_2d=False):
    def get_col(cl,ol):
        col=(np.abs(cl[:,ol])**2 if isinstance(ol,int)
             else (np.abs(cl[:,ol])**2).sum(axis=1)).round(4)
        return col
    if color_option==0:
        return []
    else:
        ham_k=[([flibs.gen_ham(k,ham_r,rvec) for k in kk] if sw_2d else flibs.gen_ham(kk,ham_r,rvec)) for kk in klist]
        if color_option==1:
            uni=[([flibs.get_uni(hk)[:,b] for hk in ham] if sw_2d else flibs.get_uni(ham)[:,b])
                 for ham,b in zip(ham_k,blist)]
            clist=[([np.array([get_col(cl,ol[0]),get_col(cl,ol[1]),get_col(cl,ol[2])]).T for cl in clst]
                    if sw_2d else np.array([get_col(clst,ol[0]),get_col(clst,ol[1]),get_col(clst,ol[2])]).T)
                   for clst in uni]
        elif color_option==2:
            uni=[([flibs.get_uni(hk) for hk in ham] if sw_2d else flibs.get_uni(ham)) for ham in ham_k]
            vk=[([flibs.get_veloc(k,ham_r,rvec,mrot,unkk)[:,b,:] for k,unkk in zip(kk,unk)]
                 if sw_2d else flibs.get_veloc(kk,ham_r,rvec,mrot,unk)[:,b,:])
                for kk,unk,b in zip(klist,uni,blist)]
            clist=[([np.sqrt((abs(vkkk)*abs(vkkk)).sum(axis=1)) for vkkk in vkk]
                    if sw_2d else np.sqrt((abs(vkk)*abs(vkk)).sum(axis=1))) for vkk in vk]
        return clist

def get_emesh(Nx,Ny,Nz,ham_r,rvec,avec,sw_uni=False,sw_veloc=False):
    Nk,klist=gen_klist(Nx,Ny,Nz,sw_pp=False)
    ham_k=flibs.gen_ham(klist,ham_r,rvec)
    eig,uni=flibs.get_eig(ham_k)
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

def get_symm_line(brav):
    if brav==0: #simple
        k_list=[[0.,0.,.5],[0., 0., 0.],[.5, 0., 0.],[.5, .5, 0.],[0.,0.,0.]]
        xlabel=['Z','$\Gamma$','X','M','$\Gamma$']
    elif brav==1: #face center
        k_list=[[0.,0.,0.],[.5, 0., .5],[1., 0., 0.],[.5, .5, .5],[.5,.25,.75],[0.,0.,0.]]
        xlabel=['$\Gamma$','X','$\Gamma$','L','W','$\Gamma$']
    elif brav==2: #body center
        k_list=[[.5,.5,.5],[0., 0., 0.],[.5, 0., 0.],[.5, .5,-.5],[0.,0.,0.]]
        xlabel=['Z','$\Gamma$','X','M','$\Gamma$']
    elif brav==3: #hexagonal
        k_list=[[0.,0.,0.],[2./3.,-1./3., 0.],[.5, 0., 0.],[0., 0., 0.],[0.,0.,.5]]
        xlabel=['$\Gamma$','K','M','$\Gamma$','Z']
    elif brav==4: #trigonal
        k_list=[[0.,0.,0.],[.5,0.,.5],[.5,0.,0.],[0.,0.,0.],[.5,.5,.5]]
        xlabel=['$\Gamma$','K','M','$\Gamma$','Z']
    else:
        pass
    return k_list,xlabel

def BZedge(brav):
    pass

def get_conductivity(mu,temp,eig,vk,Nw,Emax,idelta=1.e-3):
    wlist=np.linspace(0,Emax,Nw)
    ffermi=.5*(1.-np.tanh(.5*(eig-mu)/temp))
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
    ffermi=.5-.5*np.tanh(.5*(eig-mu)/temp)
    wlist=np.linspace(0,Emax,Nw)
    chisq=[]
    chis_orbq=[]
    f=open('chi0.dat','w')
    fq=open('writeq.dat','w')
    for i,q in enumerate(qlist):
        fq.write(f'{i:d} {q[0]:5.3f} {q[1]:5.3f} {q[2]:5.3f}\n')
        fq.flush()
        qshift=flibs.get_qshift(klist,q)
        #fkq=open(f'kq_{i:d}.dat','w')
        #for j,qs in enumerate(qshift):
        #    fkq.write(f'{qs:d}: {klist[qs-1][0]:5.3f}, {klist[qs-1][1]:5.3f}, {klist[qs-1][2]:5.3f}; {klist[j][0]:5.3f}, {klist[j][1]:5.3f}, {klist[j][2]:5.3f}\n')
        #fkq.close()
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
    ffermi=.5-.5*np.tanh(.5*(eig-mu)/temp)
    wlist=np.linspace(0,Emax,Nw)
    qshift=flibs.get_qshift(klist,q)
    chi0=flibs.get_chi_irr(uni,eig,ffermi,qshift,olist,wlist,idelta,temp)
    chis=flibs.get_chis(chi0,Smat)
    trchis,trchi0,chis_orb=flibs.get_tr_chi(chis,chi0,olist)
    return trchis,chis_orb,wlist

def chis_qmap(Nx,Ny,Ecut,mu,temp,Smat,klist,olist,eig,uni,idelta=1.e-3):
    ffermi=.5*(1.-np.tanh(.5*(eig-mu)/temp))
    chis,chi0=flibs.chis_qmap(uni,eig,ffermi,klist,Smat,olist,Nx,Ny,temp,Ecut,idelta)
    x0=np.linspace(0,1,Nx,False)
    y0=np.linspace(0,1,Ny,False)
    qx,qy=np.meshgrid(x0,y0)
    return chis,chi0,qx,qy
