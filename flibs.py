from ctypes import *
import numpy as np
#import fortran library
flibs=np.ctypeslib.load_library("fmod.so",".")
#interface for fmod subroutines
def gen_ham(klist,ham_r,rvec):
    Nk,Nr=len(klist),len(rvec)
    Norb=int(np.sqrt(ham_r.size/Nr))
    hamk=np.zeros((Nk,Norb,Norb),dtype=np.complex128)
    Nk=byref(c_int64(Nk))
    Nr=byref(c_int64(Nr))
    Norb=byref(c_int64(Norb))
    flibs.gen_ham.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128),
                            np.ctypeslib.ndpointer(dtype=np.float64),
                            np.ctypeslib.ndpointer(dtype=np.complex128),
                            np.ctypeslib.ndpointer(dtype=np.float64),
                            POINTER(c_int64),POINTER(c_int64),POINTER(c_int64)]
    flibs.gen_ham.restype=c_void_p
    flibs.gen_ham(hamk,klist,ham_r,rvec,Nk,Nr,Norb)
    return hamk

def get_eig(hamk,sw=True):
    Nk=len(hamk)
    Norb=int(np.sqrt(hamk.size/Nk))
    eig=np.zeros((Nk,Norb),dtype=np.float64)
    uni=np.zeros((Nk,Norb,Norb),dtype=np.complex128)
    Nk=byref(c_int64(Nk))
    Norb=byref(c_int64(Norb))
    flibs.get_eig.argtypes=[np.ctypeslib.ndpointer(dtype=np.float64),
                            np.ctypeslib.ndpointer(dtype=np.complex128),
                            np.ctypeslib.ndpointer(dtype=np.complex128),
                            POINTER(c_int64),POINTER(c_int64)]
    flibs.get_eig.restype=c_void_p
    flibs.get_eig(eig,uni,hamk,Nk,Norb)
    if sw:
        return eig,uni
    else:
        return eig

def get_uni(hamk):
    eig,uni=get_eig(hamk)
    return uni
    
def get_vlm0(klist,ham_r,rvec):
    Nk,Nr=len(klist),len(rvec)
    Norb=int(np.sqrt(ham_r.size/Nr))
    vk=np.zeros((Nk,Norb,Norb,3),dtype=np.complex128)
    Nk=byref(c_int64(Nk))
    Nr=byref(c_int64(Nr))
    Norb=byref(c_int64(Norb))
    flibs.get_vlm0.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128),
                            np.ctypeslib.ndpointer(dtype=np.float64),
                            np.ctypeslib.ndpointer(dtype=np.complex128),
                            np.ctypeslib.ndpointer(dtype=np.float64),
                            POINTER(c_int64),POINTER(c_int64),POINTER(c_int64)]
    flibs.get_vlm0.restype=c_void_p
    flibs.get_vlm0(vk,klist,ham_r,rvec,Nk,Nr,Norb)
    return vk

def get_vk(vk0,mrot,uni):
    Nk=len(uni)
    Norb=int(np.sqrt(uni.size/Nk))
    vk=np.zeros((Nk,Norb,3),dtype=np.float64)
    Nk=byref(c_int64(Nk))
    Norb=byref(c_int64(Norb))
    flibs.get_veloc.argtypes=[np.ctypeslib.ndpointer(dtype=np.float64),
                              np.ctypeslib.ndpointer(dtype=np.complex128),
                              np.ctypeslib.ndpointer(dtype=np.float64),
                              np.ctypeslib.ndpointer(dtype=np.complex128),
                              POINTER(c_int64),POINTER(c_int64)]
    flibs.get_veloc.restype=c_void_p
    flibs.get_veloc(vk,vk0,mrot,uni,Nk,Norb)
    return vk

def get_vnm(vk0,mrot,uni):
    Nk=len(uni)
    Norb=int(np.sqrt(uni.size/Nk))
    vk=np.zeros((Nk,Norb,Norb,3),dtype=np.complex128)
    Nk=byref(c_int64(Nk))
    Norb=byref(c_int64(Norb))
    flibs.get_vnm.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128),
                              np.ctypeslib.ndpointer(dtype=np.complex128),
                              np.ctypeslib.ndpointer(dtype=np.float64),
                              np.ctypeslib.ndpointer(dtype=np.complex128),
                              POINTER(c_int64),POINTER(c_int64)]
    flibs.get_vnm.restype=c_void_p
    flibs.get_vnm(vk,vk0,mrot,uni,Nk,Norb)
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
    Nk=byref(c_int64(Nk))
    Norb=byref(c_int64(Norb))
    flibs.get_imassk.argtypes=[np.ctypeslib.ndpointer(dtype=np.float64),
                               np.ctypeslib.ndpointer(dtype=np.complex128),
                               np.ctypeslib.ndpointer(dtype=np.float64),
                               np.ctypeslib.ndpointer(dtype=np.complex128),
                               POINTER(c_int64),POINTER(c_int64)]
    flibs.get_imassk.restype=c_void_p
    flibs.get_imassk(imass,imass0,mrot,uni,Nk,Norb)
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
    Nk=byref(c_int64(Nk))
    Nw=byref(c_int64(Nw))
    Norb=byref(c_int64(Norb))
    mu=byref(c_double(mu))
    temp=byref(c_double(temp))
    flibs.gen_green0.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128),
                               np.ctypeslib.ndpointer(dtype=np.float64),
                               np.ctypeslib.ndpointer(dtype=np.complex128),
                               POINTER(c_double),POINTER(c_double),
                               POINTER(c_int64),POINTER(c_int64),POINTER(c_int64)]
    flibs.gen_green0.restype=c_void_p
    flibs.gen_green0(Gk,eig,uni,mu,temp,Nk,Nw,Norb)
    return Gk

def gen_tr_Greenw_0(eig,mu,wlist,delta):
    Nk=len(eig)
    Norb=int(eig.size/Nk)
    Nw=len(wlist)
    trGk=np.zeros((Nk,Nw),dtype=np.complex128)
    Nk=byref(c_int64(Nk))
    Nw=byref(c_int64(Nw))
    Norb=byref(c_int64(Norb))
    mu=byref(c_double(mu))
    delta=byref(c_double(delta))
    flibs.gen_tr_greenw_0.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128), #trGk
                                    np.ctypeslib.ndpointer(dtype=np.float64),    #wlist
                                    np.ctypeslib.ndpointer(dtype=np.float64),    #eig
                                    POINTER(c_double),POINTER(c_double),         #mu,delta
                                    POINTER(c_int64),POINTER(c_int64),POINTER(c_int64)]
    flibs.gen_tr_greenw_0.restype=c_void_p
    flibs.gen_tr_greenw_0(trGk,wlist,eig,mu,delta,Nk,Nw,Norb)
    return -trGk.imag

def get_chi0_comb(Gk,kmap,olist,Nx,Ny,Nz,Nw):
    Nk=len(Gk[0,0,0])
    Norb=len(Gk)
    Nchi=len(olist)
    chi=np.zeros((Nchi,Nchi,Nw,Nk),dtype=np.complex128)
    Nx=byref(c_int64(Nx))
    Ny=byref(c_int64(Ny))
    Nz=byref(c_int64(Nz))
    Nk=byref(c_int64(Nk))
    Nw=byref(c_int64(Nw))
    Norb=byref(c_int64(Norb))
    Nchi=byref(c_int64(Nchi))
    flibs.get_chi0_comb.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128),        #chi
                                  np.ctypeslib.ndpointer(dtype=np.complex128),        #Gk
                                  np.ctypeslib.ndpointer(dtype=np.int64),             #kmap
                                  np.ctypeslib.ndpointer(dtype=np.int64),             #olist
                                  POINTER(c_int64),POINTER(c_int64),POINTER(c_int64), #Nx,Ny,Nz
                                  POINTER(c_int64),POINTER(c_int64),                  #Nw,Nk
                                  POINTER(c_int64),POINTER(c_int64)]                  #Norb,Nchi
    flibs.get_chi0_comb.restype=c_void_p
    flibs.get_chi0_comb(chi,Gk,kmap,olist,Nx,Ny,Nz,Nw,Nk,Norb,Nchi)
    return chi

def get_qshift(klist,qpoint):
    Nk=len(klist)
    qshift=np.zeros(Nk,dtype=np.int64)
    Nk=byref(c_int64(Nk))
    flibs.get_qshift.argtypes=[np.ctypeslib.ndpointer(dtype=np.float64), #qpoint
                               np.ctypeslib.ndpointer(dtype=np.float64), #klist
                               np.ctypeslib.ndpointer(dtype=np.int64),   #qshift
                               POINTER(c_int64)]                         #Nk
    flibs.get_qshift.restype=c_void_p
    flibs.get_qshift(qpoint,klist,qshift,Nk)
    return qshift

def get_chi_irr(uni,eig,ffermi,qshift,olist,wlist,idelta,temp):
    Nk=len(eig)
    Nw=len(wlist)
    Norb=int(eig.size/Nk)
    Nchi=len(olist)
    chi=np.zeros((Nw,Nchi,Nchi),dtype=np.complex128)
    eps=idelta*1e-3
    Nk=byref(c_int64(Nk))
    Norb=byref(c_int64(Norb))
    Nchi=byref(c_int64(Nchi))
    Nw=byref(c_int64(Nw))
    idelta=byref(c_double(idelta))
    eps=byref(c_double(eps))
    temp=byref(c_double(temp))
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
                                POINTER(c_double),POINTER(c_double)]         #eps,temp
    flibs.get_chi_irr.restype=c_void_p
    flibs.get_chi_irr(chi,uni,eig,ffermi,qshift,olist,wlist,Nchi,Norb,Nk,Nw,idelta,eps,temp)
    return chi

def chis_qmap(uni,eig,ffermi,klist,Smat,olist,Nx,Ny,temp,ecut,idelta):
    Nk=len(eig)
    Norb=int(eig.size/Nk)
    Nchi=len(olist)
    chi=np.zeros((Nx,Ny),dtype=np.complex128)
    chis=np.zeros((Nx,Ny),dtype=np.complex128)
    eps=idelta*1e-3
    Nx=byref(c_int64(Nx))
    Ny=byref(c_int64(Ny))
    Nk=byref(c_int64(Nk))
    Norb=byref(c_int64(Norb))
    Nchi=byref(c_int64(Nchi))
    idelta=byref(c_double(idelta))
    eps=byref(c_double(eps))
    temp=byref(c_double(temp))
    ecut=byref(c_double(ecut))
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
    flibs.chiq_map(chis,chi,uni,eig,ffermi,klist,Smat,olist,temp,ecut,idelta,eps,Nx,Ny,Nk,Norb,Nchi)
    return chis,chi

def get_chis(chi,Smat):
    Nchi=len(Smat)
    Nw=len(chi)
    Nchi=byref(c_int64(Nchi))
    Nw=byref(c_int64(Nw))
    flibs.get_chis.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128), #chi0
                             np.ctypeslib.ndpointer(dtype=np.float64),    #Smat
                             POINTER(c_int64),POINTER(c_int64)]           #Nchi,Nw
    flibs.get_chis.restype=c_void_p
    flibs.get_chis(chi,Smat,Nchi,Nw)
    return(chi)

def gen_Smatrix(olist,U,J):
    Nchi=len(olist)
    Norb=olist.max()
    Smat=np.zeros((Nchi,Nchi),dtype=np.float64)
    Nchi=byref(c_int64(Nchi))
    Norb=byref(c_int64(Norb))
    U=byref(c_double(U))
    J=byref(c_double(J))
    flibs.get_smat.argtypes=[np.ctypeslib.ndpointer(dtype=np.float64),  #Smat
                             np.ctypeslib.ndpointer(dtype=np.int64),    #olist
                             POINTER(c_double),POINTER(c_double),       #U,J
                             POINTER(c_int64),POINTER(c_int64)]         #Nchi,Norb
    flibs.get_smat.restype=c_void_p
    flibs.get_smat(Smat,olist,U,J,Nchi,Norb)
    return Smat

def calc_Lij(eig,vk,ffermi,mu,w,idelta,temp):
    Nk=len(eig)
    Norb=int(eig.size/Nk)
    L11=np.zeros((3,3),dtype=np.complex128)
    L12=np.zeros((3,3),dtype=np.complex128)
    L22=np.zeros((3,3),dtype=np.complex128)
    eps=idelta*1e-3
    Nk=byref(c_int64(Nk))
    Norb=byref(c_int64(Norb))
    w=byref(c_double(w))
    idelta=byref(c_double(idelta))
    eps=byref(c_double(eps))
    mu=byref(c_double(mu))
    temp=byref(c_double(temp))
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
    flibs.calc_lij(L11,L22,L12,vk,eig,ffermi,Norb,Nk,mu,w,idelta,eps,temp)
    return L11,L12,L22

def calc_Kn(eig,veloc,kweight,temp,mu,tau):
    Nk=len(eig)
    Norb=int(eig.size/Nk)
    K0=np.zeros((3,3),dtype=np.float64)
    K1=np.zeros((3,3),dtype=np.float64)
    K2=np.zeros((3,3),dtype=np.float64)
    Nk=byref(c_int64(Nk))
    Norb=byref(c_int64(Norb))
    temp=byref(c_double(temp))
    mu=byref(c_double(mu))
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
    flibs.calc_kn(K0,K1,K2,eig,veloc,kweight,tau,temp,mu,Nk,Norb)
    return K0,K1,K2

def calc_tdf(eig,veloc,kweight,tau,Nw):
    Nk=len(eig)
    Norb=int(eig.size/Nk)
    tdf=np.zeros((Nw,3,3),dtype=np.float64)
    Nw=byref(c_int64(Nw))
    Nk=byref(c_int64(Nk))
    Norb=byref(c_int64(Norb))
    flibs.calc_tdf.argtypes=[np.ctypeslib.ndpointer(dtype=np.float64), #tdf
                            np.ctypeslib.ndpointer(dtype=np.float64), #eig
                            np.ctypeslib.ndpointer(dtype=np.float64), #veloc
                            np.ctypeslib.ndpointer(dtype=np.float64), #kweight
                            np.ctypeslib.ndpointer(dtype=np.float64), #tau
                            POINTER(c_int64),POINTER(c_int64),        #Nw,Nk
                            POINTER(c_int64)]                         #Norb    
    flibs.calc_tdf.restype=c_void_p
    flibs.calc_tdf(tdf,eig,veloc,kweight,tau,Nw,Nk,Norb)
    return tdf
