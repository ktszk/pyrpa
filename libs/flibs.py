from ctypes import *
import numpy as np
#import fortran library
flibs=np.ctypeslib.load_library("libs/fmod.so",".")
#interface for fmod subroutines

def omp_params():
    omp_num=c_int64()
    omp_check=c_bool()
    flibs.openmp_params.argtypes=[POINTER(c_int64),POINTER(c_bool)]
    flibs.openmp_params.restype=c_void_p
    flibs.openmp_params(byref(omp_num),byref(omp_check))
    return omp_num.value,omp_check.value

def gen_ham(klist,ham_r,rvec,Ovl_r=[]):
    Nk,Nr=len(klist),len(rvec)
    Norb=int(np.sqrt(ham_r.size/Nr))
    hamk=np.zeros((Nk,Norb,Norb),dtype=np.complex128)
    flibs.gen_ham.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128),
                            np.ctypeslib.ndpointer(dtype=np.float64),
                            np.ctypeslib.ndpointer(dtype=np.complex128),
                            np.ctypeslib.ndpointer(dtype=np.float64),
                            POINTER(c_int64),POINTER(c_int64),POINTER(c_int64)]
    flibs.gen_ham.restype=c_void_p
    flibs.gen_ham(hamk,klist,ham_r,rvec,byref(c_int64(Nk)),byref(c_int64(Nr)),byref(c_int64(Norb)))
    if len(Ovl_r)!=0:
        Ovlk=np.zeros((Nk,Norb,Norb),dtype=np.complex128)
        flibs.gen_ham(Ovlk,klist,Ovl_r,rvec,byref(c_int64(Nk)),byref(c_int64(Nr)),byref(c_int64(Norb)))
        return hamk,Ovlk
    else:
        return hamk

def get_eig(hamk,Ovlk=[],sw=True):
    Nk=len(hamk)
    Norb=int(np.sqrt(hamk.size/Nk))
    eig=np.zeros((Nk,Norb),dtype=np.float64)
    uni=np.zeros((Nk,Norb,Norb),dtype=np.complex128)
    if len(Ovlk)==0:
        flibs.get_eig.argtypes=[np.ctypeslib.ndpointer(dtype=np.float64),    #eig
                                np.ctypeslib.ndpointer(dtype=np.complex128), #uni
                                np.ctypeslib.ndpointer(dtype=np.complex128), #hamk
                                POINTER(c_int64),POINTER(c_int64)]           #Nk,Norb
        flibs.get_eig.restype=c_void_p
        flibs.get_eig(eig,uni,hamk,byref(c_int64(Nk)),byref(c_int64(Norb)))
    else:
        flibs.get_eig_mlo.argtypes=[np.ctypeslib.ndpointer(dtype=np.float64),    #eig
                                    np.ctypeslib.ndpointer(dtype=np.complex128), #uni
                                    np.ctypeslib.ndpointer(dtype=np.complex128), #hamk
                                    np.ctypeslib.ndpointer(dtype=np.complex128), #Ovlk
                                    POINTER(c_int64),POINTER(c_int64)]           #Nk,Norb
        flibs.get_eig_mlo.restype=c_void_p
        flibs.get_eig_mlo(eig,uni,hamk,Ovlk,byref(c_int64(Nk)),byref(c_int64(Norb)))
    if sw:
        return eig,uni
    else:
        return eig

def get_uni(hamk,Ovlk=[]):
    eig,uni=get_eig(hamk,Ovlk)
    return uni
    
def get_vlm0(klist,ham_r,rvec):
    Nk,Nr=len(klist),len(rvec)
    Norb=int(np.sqrt(ham_r.size/Nr))
    vk=np.zeros((Nk,Norb,Norb,3),dtype=np.complex128)
    flibs.get_vlm0.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128),
                            np.ctypeslib.ndpointer(dtype=np.float64),
                            np.ctypeslib.ndpointer(dtype=np.complex128),
                            np.ctypeslib.ndpointer(dtype=np.float64),
                            POINTER(c_int64),POINTER(c_int64),POINTER(c_int64)]
    flibs.get_vlm0.restype=c_void_p
    flibs.get_vlm0(vk,klist,ham_r,rvec,byref(c_int64(Nk)),byref(c_int64(Nr)),byref(c_int64(Norb)))
    return vk

def get_vk(vk0,mrot,uni):
    Nk=len(uni)
    Norb=int(np.sqrt(uni.size/Nk))
    vk=np.zeros((Nk,Norb,3),dtype=np.float64)
    flibs.get_veloc.argtypes=[np.ctypeslib.ndpointer(dtype=np.float64),
                              np.ctypeslib.ndpointer(dtype=np.complex128),
                              np.ctypeslib.ndpointer(dtype=np.float64),
                              np.ctypeslib.ndpointer(dtype=np.complex128),
                              POINTER(c_int64),POINTER(c_int64)]
    flibs.get_veloc.restype=c_void_p
    flibs.get_veloc(vk,vk0,mrot,uni,byref(c_int64(Nk)),byref(c_int64(Norb)))
    return vk

def get_vnm(vk0,mrot,uni):
    Nk=len(uni)
    Norb=int(np.sqrt(uni.size/Nk))
    vk=np.zeros((Nk,Norb,Norb,3),dtype=np.complex128)
    flibs.get_vnm.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128),
                              np.ctypeslib.ndpointer(dtype=np.complex128),
                              np.ctypeslib.ndpointer(dtype=np.float64),
                              np.ctypeslib.ndpointer(dtype=np.complex128),
                              POINTER(c_int64),POINTER(c_int64)]
    flibs.get_vnm.restype=c_void_p
    flibs.get_vnm(vk,vk0,mrot,uni,byref(c_int64(Nk)),byref(c_int64(Norb)))
    return vk

def get_vnmk(klist,ham_r,rvec,mrot,uni):
    vk0=get_vlm0(klist,ham_r,rvec)
    vk=get_vnm(vk0,mrot,uni)
    return vk

def get_veloc(klist,ham_r,rvec,mrot,uni):
    vk0=get_vlm0(klist,ham_r,rvec)
    vk=get_vk(vk0,mrot,uni)
    return vk

def get_imass0(klist,ham_r,rvec):
    Nk,Nr=len(klist),len(rvec)
    Norb=int(np.sqrt(ham_r.size/Nr))
    imass0=np.zeros((Nk,Norb,Norb,3,3),dtype=np.complex128)
    Nk=byref(c_int64(Nk))
    Nr=byref(c_int64(Nr))
    Norb=byref(c_int64(Norb))
    flibs.get_imass0.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128),
                               np.ctypeslib.ndpointer(dtype=np.float64),
                               np.ctypeslib.ndpointer(dtype=np.complex128),
                               np.ctypeslib.ndpointer(dtype=np.float64),
                               POINTER(c_int64),POINTER(c_int64),
                               POINTER(c_int64)]
    flibs.get_imass0.restype=c_void_p
    flibs.get_imass0(imass0,klist,ham_r,rvec,Nk,Nr,Norb)
    return imass0

def get_imassk(imass0,mrot,uni):
    Nk=len(uni)
    Norb=int(np.sqrt(uni.size/Nk))
    imass=np.zeros((Nk,Norb,3,3),dtype=np.float64)
    flibs.get_imassk.argtypes=[np.ctypeslib.ndpointer(dtype=np.float64),
                               np.ctypeslib.ndpointer(dtype=np.complex128),
                               np.ctypeslib.ndpointer(dtype=np.float64),
                               np.ctypeslib.ndpointer(dtype=np.complex128),
                               POINTER(c_int64),POINTER(c_int64)]
    flibs.get_imassk.restype=c_void_p
    flibs.get_imassk(imass,imass0,mrot,uni,byref(c_int64(Nk)),byref(c_int64(Norb)))
    return imass

def get_mass(klist,ham_r,rvec,mrot,uni):
    import scipy.linalg as sclin
    imass0=get_imass0(klist,ham_r,rvec)
    imass=get_imassk(imass0,mrot,uni)
    mass=np.array([[sclin.inv(im) for im in imas] for imas in imass])
    return mass

def gen_Green0(eig,uni,mu,temp,Nw):
    Nk=len(eig)
    Norb=int(eig.size/Nk)
    Gk=np.zeros((Norb,Norb,Nw,Nk),dtype=np.complex128)
    flibs.gen_green0.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128),
                               np.ctypeslib.ndpointer(dtype=np.float64),
                               np.ctypeslib.ndpointer(dtype=np.complex128),
                               POINTER(c_double),POINTER(c_double),
                               POINTER(c_int64),POINTER(c_int64),POINTER(c_int64)]
    flibs.gen_green0.restype=c_void_p
    flibs.gen_green0(Gk,eig,uni,byref(c_double(mu)),byref(c_double(temp)),
                     byref(c_int64(Nk)),byref(c_int64(Nw)),byref(c_int64(Norb)))
    return Gk

def gen_tr_Greenw_0(eig,mu,wlist,delta):
    Nk,Nw=len(eig),len(wlist)
    Norb=int(eig.size/Nk)
    trGk=np.zeros((Nk,Nw),dtype=np.complex128)
    flibs.gen_tr_greenw_0.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128), #trGk
                                    np.ctypeslib.ndpointer(dtype=np.float64),    #wlist
                                    np.ctypeslib.ndpointer(dtype=np.float64),    #eig
                                    POINTER(c_double),POINTER(c_double),         #mu,delta
                                    POINTER(c_int64),POINTER(c_int64),POINTER(c_int64)]
    flibs.gen_tr_greenw_0.restype=c_void_p
    flibs.gen_tr_greenw_0(trGk,wlist,eig,byref(c_double(mu)),byref(c_double(delta)),
                          byref(c_int64(Nk)),byref(c_int64(Nw)),byref(c_int64(Norb)))
    return -trGk.imag

def get_chi0_comb(Gk,kmap,olist,Nx,Ny,Nz,Nw):
    Nk=len(Gk[0,0,0])
    Norb,Nchi=len(Gk),len(olist)
    chi=np.zeros((Nchi,Nchi,Nw,Nk),dtype=np.complex128)
    flibs.get_chi0_comb.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128),        #chi
                                  np.ctypeslib.ndpointer(dtype=np.complex128),        #Gk
                                  np.ctypeslib.ndpointer(dtype=np.int64),             #kmap
                                  np.ctypeslib.ndpointer(dtype=np.int64),             #olist
                                  POINTER(c_int64),POINTER(c_int64),POINTER(c_int64), #Nx,Ny,Nz
                                  POINTER(c_int64),POINTER(c_int64),                  #Nw,Nk
                                  POINTER(c_int64),POINTER(c_int64)]                  #Norb,Nchi
    flibs.get_chi0_comb.restype=c_void_p
    flibs.get_chi0_comb(chi,Gk,kmap,olist,byref(c_int64(Nx)),byref(c_int64(Ny)),
                        byref(c_int64(Nz)),byref(c_int64(Nw)),byref(c_int64(Nk)),
                        byref(c_int64(Norb)),byref(c_int64(Nchi)))
    #flibs.get_chi0_comb(chi,Gk,olist,byref(c_int64(Nw)),byref(c_int64(Nk)),byref(c_int64(Norb)),byref(c_int64(Nchi))) #,byref(c_int64(Nx)),byref(c_int64(Ny)),byref(c_int64(Nz)))
    return chi

def get_qshift(klist,qpoint):
    Nk=len(klist)
    qshift=np.zeros(Nk,dtype=np.int64)
    flibs.get_qshift.argtypes=[np.ctypeslib.ndpointer(dtype=np.float64), #qpoint
                               np.ctypeslib.ndpointer(dtype=np.float64), #klist
                               np.ctypeslib.ndpointer(dtype=np.int64),   #qshift
                               POINTER(c_int64)]                         #Nk
    flibs.get_qshift.restype=c_void_p
    flibs.get_qshift(qpoint,klist,qshift,byref(c_int64(Nk)))
    return qshift

def get_iqshift(klist,qpoint):
    Nk=len(klist)
    qshift=np.zeros(Nk,dtype=np.int64)
    flibs.get_iqshift.argtypes=[np.ctypeslib.ndpointer(dtype=np.float64), #qpoint
                               np.ctypeslib.ndpointer(dtype=np.float64), #klist
                               np.ctypeslib.ndpointer(dtype=np.int64),   #qshift
                               POINTER(c_int64)]                         #Nk
    flibs.get_iqshift.restype=c_void_p
    flibs.get_iqshift(qpoint,klist,qshift,byref(c_int64(Nk)))
    return qshift

def get_chi_irr(uni,eig,ffermi,qshift,olist,wlist,idelta,temp,i):
    Nk,Nw=len(eig),len(wlist)
    Norb,Nchi=int(eig.size/Nk),len(olist)
    chi=np.zeros((Nw,Nchi,Nchi),dtype=np.complex128)
    eps=idelta*1e-3
    flibs.get_chi_irr.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128), #chi
                                np.ctypeslib.ndpointer(dtype=np.complex128), #uni
                                np.ctypeslib.ndpointer(dtype=np.float64),    #eig
                                np.ctypeslib.ndpointer(dtype=np.float64),    #ffermi
                                np.ctypeslib.ndpointer(dtype=np.int64),      #qshift
                                np.ctypeslib.ndpointer(dtype=np.int64),      #olist
                                np.ctypeslib.ndpointer(dtype=np.float64),    #wlist
                                POINTER(c_int64),POINTER(c_int64),           #Nchi,Norb
                                POINTER(c_int64),POINTER(c_int64),           #Nk,Nw
                                POINTER(c_double),                           #idelta
                                POINTER(c_double),POINTER(c_double)
                                ,POINTER(c_int64)]         #eps,temp
    flibs.get_chi_irr.restype=c_void_p
    flibs.get_chi_irr(chi,uni,eig,ffermi,qshift,olist,wlist,byref(c_int64(Nchi)),
                      byref(c_int64(Norb)),byref(c_int64(Nk)),byref(c_int64(Nw)),
                      byref(c_double(idelta)),byref(c_double(eps)),byref(c_double(temp))
                      ,byref(c_int64(i)))
    return chi

def chis_qmap(uni,eig,ffermi,klist,Smat,olist,Nx,Ny,temp,ecut,idelta):
    Nk=len(eig)
    Norb,Nchi=int(eig.size/Nk),len(olist)
    chi=np.zeros((Nx,Ny),dtype=np.complex128)
    chis=np.zeros((Nx,Ny),dtype=np.complex128)
    eps=idelta*1e-3
    flibs.chiq_map.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128),        #chis
                             np.ctypeslib.ndpointer(dtype=np.complex128),        #chi
                             np.ctypeslib.ndpointer(dtype=np.complex128),        #uni
                             np.ctypeslib.ndpointer(dtype=np.float64),           #eig
                             np.ctypeslib.ndpointer(dtype=np.float64),           #ffermi
                             np.ctypeslib.ndpointer(dtype=np.float64),           #klist
                             np.ctypeslib.ndpointer(dtype=np.float64),           #Smat
                             np.ctypeslib.ndpointer(dtype=np.int64),             #olist
                             POINTER(c_double),POINTER(c_double),                #temp,ecut
                             POINTER(c_double),POINTER(c_double),                #idelta,eps
                             POINTER(c_int64),POINTER(c_int64),                  #Nx,Ny
                             POINTER(c_int64),POINTER(c_int64),POINTER(c_int64)] #Nk,Norb,Nchi
    flibs.chiq_map.restype=c_void_p
    flibs.chiq_map(chis,chi,uni,eig,ffermi,klist,Smat,olist,byref(c_double(temp)),
                   byref(c_double(ecut)),byref(c_double(idelta)),byref(c_double(eps)),
                   byref(c_int64(Nx)),byref(c_int64(Ny)),byref(c_int64(Nk)),
                   byref(c_int64(Norb)),byref(c_int64(Nchi)))
    return chis,chi

def phi_qmap(uni,eig,ffermi,klist,olist,Nx,Ny,temp,ecut,idelta):
    Nk=len(eig)
    Norb,Nchi=int(eig.size/Nk),len(olist)
    phi=np.zeros((Nx,Ny),dtype=np.complex128)
    eps=idelta*1e-3
    flibs.phiq_map.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128),        #phi
                             np.ctypeslib.ndpointer(dtype=np.complex128),        #uni
                             np.ctypeslib.ndpointer(dtype=np.float64),           #eig
                             np.ctypeslib.ndpointer(dtype=np.float64),           #ffermi
                             np.ctypeslib.ndpointer(dtype=np.float64),           #klist
                             np.ctypeslib.ndpointer(dtype=np.int64),             #olist
                             POINTER(c_double),POINTER(c_double),                #temp,ecut
                             POINTER(c_double),POINTER(c_double),                #idelta,eps
                             POINTER(c_int64),POINTER(c_int64),                  #Nx,Ny
                             POINTER(c_int64),POINTER(c_int64),POINTER(c_int64)] #Nk,Norb,Nchi
    flibs.phiq_map.restype=c_void_p
    flibs.phiq_map(phi,uni,eig,ffermi,klist,olist,byref(c_double(temp)),
                   byref(c_double(ecut)),byref(c_double(idelta)),byref(c_double(eps)),
                   byref(c_int64(Nx)),byref(c_int64(Ny)),byref(c_int64(Nk)),
                   byref(c_int64(Norb)),byref(c_int64(Nchi)))
    return phi

def get_tr_chi(chis,chi0,olist):
    Nchi,Nw,Norb=len(olist),len(chi0),olist.max()
    trchis=np.zeros(Nw,dtype=np.complex128)
    trchi0=np.zeros(Nw,dtype=np.complex128)
    chis_orb=np.zeros((Nw,Norb+2),dtype=np.complex128)
    flibs.get_tr_chi.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128), #trchis
                               np.ctypeslib.ndpointer(dtype=np.complex128), #trchi0
                               np.ctypeslib.ndpointer(dtype=np.complex128), #chis_orb
                               np.ctypeslib.ndpointer(dtype=np.complex128), #chis
                               np.ctypeslib.ndpointer(dtype=np.complex128), #chi0
                               np.ctypeslib.ndpointer(dtype=np.int64),      #olist
                               POINTER(c_int64),POINTER(c_int64),           #Nw,Nchi
                               POINTER(c_int64)]                            #Norb
    flibs.get_tr_chi.restype=c_void_p
    flibs.get_tr_chi(trchis,trchi0,chis_orb,chis,chi0,olist,byref(c_int64(Nw)),byref(c_int64(Nchi)),byref(c_int64(Norb)))
    return trchis,trchi0,chis_orb

def get_tr_phi(phi,olist):
    Nchi,Nw,Norb=len(olist),len(phi),olist.max()
    trphi=np.zeros(Nw,dtype=np.complex128)
    phi_orb=np.zeros((Nw,Norb+2),dtype=np.complex128)
    flibs.get_tr_phi.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128), #trphi
                               np.ctypeslib.ndpointer(dtype=np.complex128), #phi_orb
                               np.ctypeslib.ndpointer(dtype=np.complex128), #phi
                               np.ctypeslib.ndpointer(dtype=np.int64),      #olist
                               POINTER(c_int64),POINTER(c_int64),           #Nw,Nchi
                               POINTER(c_int64)]                            #Norb
    flibs.get_tr_phi.restype=c_void_p
    flibs.get_tr_phi(trphi,phi_orb,phi,olist,byref(c_int64(Nw)),
                     byref(c_int64(Nchi)),byref(c_int64(Norb)))
    return trphi,phi_orb

def get_chis(chi0,Smat):
    Nchi,Nw=len(Smat),len(chi0)
    chis=np.zeros((Nw,Nchi,Nchi),dtype=np.complex128)
    flibs.get_chis.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128), #chis
                             np.ctypeslib.ndpointer(dtype=np.complex128), #chi0
                             np.ctypeslib.ndpointer(dtype=np.float64),    #Smat
                             POINTER(c_int64),POINTER(c_int64)]           #Nchi,Nw
    flibs.get_chis.restype=c_void_p
    flibs.get_chis(chis,chi0,Smat,byref(c_int64(Nchi)),byref(c_int64(Nw)))
    return(chis)

def gen_Smatrix(olist,U,J):
    Nchi,Norb=len(olist),olist.max()
    Smat=np.zeros((Nchi,Nchi),dtype=np.float64)
    flibs.get_smat.argtypes=[np.ctypeslib.ndpointer(dtype=np.float64),  #Smat
                             np.ctypeslib.ndpointer(dtype=np.int64),    #olist
                             POINTER(c_double),POINTER(c_double),       #U,J
                             POINTER(c_int64),POINTER(c_int64)]         #Nchi,Norb
    flibs.get_smat.restype=c_void_p
    flibs.get_smat(Smat,olist,byref(c_double(U)),byref(c_double(J)),
                   byref(c_int64(Nchi)),byref(c_int64(Norb)))
    return Smat

def calc_Lij(eig,vk,ffermi,mu,w,idelta,temp):
    Nk=len(eig)
    Norb=int(eig.size/Nk)
    L11=np.zeros((3,3),dtype=np.complex128)
    L12=np.zeros((3,3),dtype=np.complex128)
    L22=np.zeros((3,3),dtype=np.complex128)
    eps=idelta*1e-3
    flibs.calc_lij.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128),    #L11
                                np.ctypeslib.ndpointer(dtype=np.complex128), #L12
                                np.ctypeslib.ndpointer(dtype=np.complex128), #L22
                                np.ctypeslib.ndpointer(dtype=np.complex128), #vk
                                np.ctypeslib.ndpointer(dtype=np.float64),    #eig
                                np.ctypeslib.ndpointer(dtype=np.float64),    #ffermi
                                POINTER(c_int64),POINTER(c_int64),           #Norb,Nk
                                POINTER(c_double),                           #mu
                                POINTER(c_double),POINTER(c_double),         #w,idelta
                                POINTER(c_double),POINTER(c_double)]         #eps,temp
    flibs.calc_lij.restype=c_void_p
    flibs.calc_lij(L11,L22,L12,vk,eig,ffermi,byref(c_int64(Norb)),byref(c_int64(Nk)),
                   byref(c_double(mu)),byref(c_double(w)),byref(c_double(idelta)),
                   byref(c_double(eps)),byref(c_double(temp)))
    return L11,L12,L22

def calc_Kn(eig,veloc,kweight,temp,mu,tau):
    Nk=len(eig)
    Norb=int(eig.size/Nk)
    K0=np.zeros((3,3),dtype=np.float64)
    K1=np.zeros((3,3),dtype=np.float64)
    K2=np.zeros((3,3),dtype=np.float64)
    flibs.calc_kn.argtypes=[np.ctypeslib.ndpointer(dtype=np.float64), #K0
                            np.ctypeslib.ndpointer(dtype=np.float64), #K1
                            np.ctypeslib.ndpointer(dtype=np.float64), #K2
                            np.ctypeslib.ndpointer(dtype=np.float64), #eig
                            np.ctypeslib.ndpointer(dtype=np.float64), #veloc
                            np.ctypeslib.ndpointer(dtype=np.float64), #kweight
                            np.ctypeslib.ndpointer(dtype=np.float64), #tau
                            POINTER(c_double),POINTER(c_double),      #temp,mu
                            POINTER(c_int64),POINTER(c_int64)]        #Nk,Norb
    flibs.calc_kn.restype=c_void_p
    flibs.calc_kn(K0,K1,K2,eig,veloc,kweight,tau,byref(c_double(temp)),
                  byref(c_double(mu)),byref(c_int64(Nk)),byref(c_int64(Norb)))
    return K0,K1,K2

def calc_tdf(eig,veloc,kweight,tau,Nw):
    Nk=len(eig)
    Norb=int(eig.size/Nk)
    tdf=np.zeros((Nw,3,3),dtype=np.float64)
    flibs.calc_tdf.argtypes=[np.ctypeslib.ndpointer(dtype=np.float64), #tdf
                            np.ctypeslib.ndpointer(dtype=np.float64), #eig
                            np.ctypeslib.ndpointer(dtype=np.float64), #veloc
                            np.ctypeslib.ndpointer(dtype=np.float64), #kweight
                            np.ctypeslib.ndpointer(dtype=np.float64), #tau
                            POINTER(c_int64),POINTER(c_int64),        #Nw,Nk
                            POINTER(c_int64)]                         #Norb    
    flibs.calc_tdf.restype=c_void_p
    flibs.calc_tdf(tdf,eig,veloc,kweight,tau,byref(c_int64(Nw)),
                   byref(c_int64(Nk)),byref(c_int64(Norb)))
    return tdf
