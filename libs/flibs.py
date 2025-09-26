from ctypes import *
import numpy as np
#import fortran library
flibs=np.ctypeslib.load_library("libs/libfmod.so",".")
#interface for fmod subroutines

def omp_params() ->tuple[int,bool]:
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

def get_ffermi(eig,mu:float,temp:float):
    Nk=len(eig)
    Norb=int(eig.size/Nk)
    ffermi=np.zeros((Nk,Norb),dtype=np.float64)
    flibs.get_ffermi.argtypes=[np.ctypeslib.ndpointer(dtype=np.float64), #ffermi
                               np.ctypeslib.ndpointer(dtype=np.float64), #eig
                               POINTER(c_double),POINTER(c_double),      #mu,temp
                               POINTER(c_int64),POINTER(c_int64)]        #Nk,Norb
    flibs.get_ffermi.restype=c_void_p
    flibs.get_ffermi(ffermi,eig,byref(c_double(mu)),byref(c_double(temp)),
                     byref(c_int64(Nk)),byref(c_int64(Norb)))
    return ffermi

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

def gen_Green0(eig,uni,mu:float,temp:float,Nw:int):
    Nk=len(eig)
    Norb=int(eig.size/Nk)
    Gk=np.zeros((Norb,Norb,Nw,Nk),dtype=np.complex128)
    flibs.gen_green0_.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128),        #Gk
                                np.ctypeslib.ndpointer(dtype=np.float64),           #eig
                                np.ctypeslib.ndpointer(dtype=np.complex128),        #uni
                                POINTER(c_double),POINTER(c_double),                #mu,temp
                                POINTER(c_int64),POINTER(c_int64),POINTER(c_int64)] #Nk,Nw,Norb
    flibs.gen_green0_.restype=c_void_p
    flibs.gen_green0_(Gk,eig,uni,byref(c_double(mu)),byref(c_double(temp)),
                     byref(c_int64(Nk)),byref(c_int64(Nw)),byref(c_int64(Norb)))
    return Gk

def gen_green(selfen,hamk,mu:float,temp:float):
    Nk=len(hamk)
    Norb=int(np.sqrt(hamk.size/Nk))
    Nw=int(selfen.size/(Nk*Norb*Norb))
    Gk=np.zeros((Norb,Norb,Nw,Nk),dtype=np.complex128)
    flibs.gen_green_inv_.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128),        #Gk
                                   np.ctypeslib.ndpointer(dtype=np.complex128),        #selfen
                                   np.ctypeslib.ndpointer(dtype=np.complex128),        #hamk
                                   POINTER(c_double),POINTER(c_double),                #mu,temp
                                   POINTER(c_int64),POINTER(c_int64),POINTER(c_int64)] #Nk,Nw,Norb
    flibs.gen_green_inv_.restype=c_void_p
    flibs.gen_green_inv_(Gk,selfen,hamk,byref(c_double(mu)),byref(c_double(temp)),
                        byref(c_int64(Nk)),byref(c_int64(Nw)),byref(c_int64(Norb)))
    flibs.getinv_.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128),        #Gk
                           POINTER(c_int64),POINTER(c_int64),POINTER(c_int64)] #Nk,Nw,Norb
    flibs.getinv_.restype=c_void_p
    flibs.getinv_(Gk,byref(c_int64(Nk)),byref(c_int64(Nw)),byref(c_int64(Norb)))
    return Gk

def gen_green_from_eig(selfen,eig,uni,mu:float,temp:float):
    Nk=len(eig)
    Norb=int(eig.size/Nk)
    Nw=int(selfen.size/(Nk*Norb*Norb))
    Gk=np.zeros((Norb,Norb,Nw,Nk),dtype=np.complex128)
    flibs.gen_green_inv_from_eig.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128),        #Gk
                                           np.ctypeslib.ndpointer(dtype=np.complex128),        #selfen
                                           np.ctypeslib.ndpointer(dtype=np.complex128),        #uni
                                           np.ctypeslib.ndpointer(dtype=np.float64),           #eig
                                           POINTER(c_double),POINTER(c_double),                #mu,temp
                                           POINTER(c_int64),POINTER(c_int64),POINTER(c_int64)] #Nk,Nw,Norb
    flibs.gen_green_inv_from_eig.restype=c_void_p
    flibs.gen_green_inv_from_eig(Gk,selfen,uni,eig,byref(c_double(mu)),byref(c_double(temp)),
                                 byref(c_int64(Nk)),byref(c_int64(Nw)),byref(c_int64(Norb)))
    flibs.getinv.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128),        #Gk
                           POINTER(c_int64),POINTER(c_int64),POINTER(c_int64)] #Nk,Nw,Norb
    flibs.getinv.restype=c_void_p
    flibs.getinv(Gk,byref(c_int64(Nk)),byref(c_int64(Nw)),byref(c_int64(Norb)))
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

def gen_dos(eig,uni,mu,wlist,delta):
    Nk,Nw=len(eig),len(wlist)
    Norb=int(eig.size/Nk)
    pDos=np.zeros((Norb,Nw),dtype=np.complex128)
    flibs.gen_dos.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128), #Dos
                            np.ctypeslib.ndpointer(dtype=np.float64),    #wlist
                            np.ctypeslib.ndpointer(dtype=np.float64),    #eig
                            np.ctypeslib.ndpointer(dtype=np.complex128), #uni
                            POINTER(c_double),POINTER(c_double),         #mu,delta
                            POINTER(c_int64),POINTER(c_int64),POINTER(c_int64)]
    flibs.gen_dos.restype=c_void_p
    flibs.gen_dos(pDos,wlist,eig,uni,byref(c_double(mu)),byref(c_double(delta)),
                          byref(c_int64(Nk)),byref(c_int64(Nw)),byref(c_int64(Norb)))
    return -pDos.imag

def get_chi0_conv(Gk,kmap,invk,olist,temp:float,Nx:int,Ny:int,Nz:int):
    Nkall,Nk=len(kmap),len(Gk[0,0,0])
    Norb,Nchi=len(Gk),len(olist)
    Nw=int(Gk.size/(Norb*Norb*Nk))
    chi=np.zeros((Nchi,Nchi,Nw,Nk),dtype=np.complex128)
    flibs.get_chi0_conv_.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128),         #chi
                                   np.ctypeslib.ndpointer(dtype=np.complex128),         #Gk
                                   np.ctypeslib.ndpointer(dtype=np.int64),              #kmap
                                   np.ctypeslib.ndpointer(dtype=np.int64),              #invk
                                   np.ctypeslib.ndpointer(dtype=np.int64),              #olist
                                   POINTER(c_double),POINTER(c_int64),POINTER(c_int64), #temp,Nx,Ny
                                   POINTER(c_int64),POINTER(c_int64),POINTER(c_int64),  #Nz,Nw,Nk
                                   POINTER(c_int64),POINTER(c_int64),POINTER(c_int64)]  #Nkall,Norb,Nchi
    flibs.get_chi0_conv_.restype=c_void_p
    flibs.get_chi0_conv_(chi,Gk,kmap,invk,olist,byref(c_double(temp)),byref(c_int64(Nx)),
                         byref(c_int64(Ny)),byref(c_int64(Nz)),byref(c_int64(Nw)),byref(c_int64(Nk)),
                         byref(c_int64(Nkall)),byref(c_int64(Norb)),byref(c_int64(Nchi)))
    return chi

def get_chi0_sum(Gk,invk,klist,olist,temp:float):
    Nkall,Nk=len(invk),len(klist)
    Norb,Nchi=len(Gk),len(olist)
    Nw=int(Gk.size/(Norb*Norb*Nk))
    chi=np.zeros((Nchi,Nchi,Nw,Nk),dtype=np.complex128)
    flibs.get_chi0_sum.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128),   #chi
                                 np.ctypeslib.ndpointer(dtype=np.complex128),   #Gk
                                 np.ctypeslib.ndpointer(dtype=np.float64),      #klist
                                 np.ctypeslib.ndpointer(dtype=np.int64),        #olist
                                 POINTER(c_double),POINTER(c_int64),            #temp,Nw
                                 POINTER(c_int64),POINTER(c_int64),             #Nk,Nkall
                                 POINTER(c_int64),POINTER(c_int64)]             #Norb,Nchi
    flibs.get_chi0_sum.restype=c_void_p
    flibs.get_chi0_sum(chi,Gk,klist,olist,byref(c_double(temp)),byref(c_int64(Nw)),
                       byref(c_int64(Nk)),byref(c_int64(Nkall)),byref(c_int64(Norb)),
                       byref(c_int64(Nchi)))
    return chi

def get_Vsigma_nosoc_flex(chi,Smat,Cmat):
    Nk,Nw,Nchi=len(chi),len(chi[0]),len(Smat)
    flibs.get_vsigma_flex_nosoc_.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128),
                                           np.ctypeslib.ndpointer(dtype=np.float64),
                                           np.ctypeslib.ndpointer(dtype=np.float64),
                                           POINTER(c_int64),POINTER(c_int64),POINTER(c_int64)]
    flibs.get_vsigma_flex_nosoc_.restype=c_void_p
    flibs.get_vsigma_flex_nosoc_(chi,Smat,Cmat,byref(c_int64(Nk)),
                                 byref(c_int64(Nw)),byref(c_int64(Nchi)))
    return chi.copy()

def mkself(Smat,Cmat,kmap,invk,olist,hamk,eig,uni,mu:float,fill:float,temp:float,
           Nw:int,Nx:int,Ny:int,Nz:int,sw_out:bool,sw_in:bool,sw_sub_sigma=True,scf_loop=300,eps=1.0e-4,pp=0.3):
    print('mixing rate: pp = %3.1f'%pp)
    Nkall,Nk,Nchi=len(kmap),len(hamk),len(Smat)
    Norb=int(np.sqrt(hamk.size/Nk))
    mu_self=c_double()
    sigmak=np.zeros((Norb,Norb,Nw,Nk),dtype=np.complex128)
    flibs.mkself.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128),          #sigmak
                           POINTER(c_double),                                    #muself
                           np.ctypeslib.ndpointer(dtype=np.float64),             #Smat
                           np.ctypeslib.ndpointer(dtype=np.float64),             #Cmat
                           np.ctypeslib.ndpointer(dtype=np.int64),               #kmap
                           np.ctypeslib.ndpointer(dtype=np.int64),               #invk
                           np.ctypeslib.ndpointer(dtype=np.int64),               #olist
                           np.ctypeslib.ndpointer(dtype=np.complex128),          #hamk
                           np.ctypeslib.ndpointer(dtype=np.float64),             #eig
                           np.ctypeslib.ndpointer(dtype=np.complex128),          #uni
                           POINTER(c_double),POINTER(c_double),                  #mu,fill
                           POINTER(c_double),POINTER(c_int64),                   #temp,scf_loop
                           POINTER(c_double),POINTER(c_double),POINTER(c_int64), #pp,eps,Nkall
                           POINTER(c_int64),POINTER(c_int64),                    #Nk,Nw
                           POINTER(c_int64),POINTER(c_int64),                    #Nchi,Norb
                           POINTER(c_int64),POINTER(c_int64),POINTER(c_int64),   #Nx,Ny,Nz
                           POINTER(c_bool),POINTER(c_bool),POINTER(c_bool)]      #sw_sub_sigma,sw_out,sw_in
    flibs.mkself.restype=c_void_p
    flibs.mkself(sigmak,byref(mu_self),Smat,Cmat,kmap,invk,olist,hamk,eig,uni,byref(c_double(mu)),
                 byref(c_double(fill)),byref(c_double(temp)),byref(c_int64(scf_loop)),
                 byref(c_double(pp)),byref(c_double(eps)),byref(c_int64(Nkall)),
                 byref(c_int64(Nk)),byref(c_int64(Nw)),byref(c_int64(Norb)),
                 byref(c_int64(Nchi)),byref(c_int64(Nx)),byref(c_int64(Ny)),byref(c_int64(Nz)),
                 byref(c_bool(sw_sub_sigma)),byref(c_bool(sw_out)),byref(c_bool(sw_in)))
    return sigmak,mu_self.value

def get_qshift(klist,qpoint):
    Nk=len(klist)
    qshift=np.zeros(Nk,dtype=np.int64)
    flibs.get_qshift_.argtypes=[np.ctypeslib.ndpointer(dtype=np.float64), #qpoint
                                np.ctypeslib.ndpointer(dtype=np.float64), #klist
                                np.ctypeslib.ndpointer(dtype=np.int64),   #qshift
                                POINTER(c_int64)]                         #Nk
    flibs.get_qshift_.restype=c_void_p
    flibs.get_qshift_(qpoint,klist,qshift,byref(c_int64(Nk)))
    return qshift

def get_iqshift(klist,qpoint):
    Nk=len(klist)
    qshift=np.zeros(Nk,dtype=np.int64)
    flibs.get_iqshift_.argtypes=[np.ctypeslib.ndpointer(dtype=np.float64), #qpoint
                                 np.ctypeslib.ndpointer(dtype=np.float64), #klist
                                 np.ctypeslib.ndpointer(dtype=np.int64),   #qshift
                                 POINTER(c_int64)]                         #Nk
    flibs.get_iqshift_.restype=c_void_p
    flibs.get_iqshift_(qpoint,klist,qshift,byref(c_int64(Nk)))
    return qshift

def get_chi_irr(uni,eig,ffermi,qshift,olist,wlist,idelta,temp):
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
                                POINTER(c_double),POINTER(c_double)]         #eps,temp
    flibs.get_chi_irr.restype=c_void_p
    flibs.get_chi_irr(chi,uni,eig,ffermi,qshift,olist,wlist,byref(c_int64(Nchi)),
                      byref(c_int64(Norb)),byref(c_int64(Nk)),byref(c_int64(Nw)),
                      byref(c_double(idelta)),byref(c_double(eps)),byref(c_double(temp)))
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

def get_phi_irr(uni,eig,ffermi,qshift,olist,wlist,idelta,mu,temp):
    Nk,Nw=len(eig),len(wlist)
    Norb,Nchi=int(eig.size/Nk),len(olist)
    phi=np.zeros((Nw,Nchi,Nchi),dtype=np.complex128)
    eps=idelta*1e-3
    flibs.get_phi_irr.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128), #phi
                                np.ctypeslib.ndpointer(dtype=np.complex128), #uni
                                np.ctypeslib.ndpointer(dtype=np.float64),    #eig
                                np.ctypeslib.ndpointer(dtype=np.float64),    #ffermi
                                np.ctypeslib.ndpointer(dtype=np.int64),      #qshift
                                np.ctypeslib.ndpointer(dtype=np.int64),      #olist
                                np.ctypeslib.ndpointer(dtype=np.float64),    #wlist
                                POINTER(c_int64),POINTER(c_int64),           #Nchi,Norb
                                POINTER(c_int64),POINTER(c_int64),           #Nk,Nw
                                POINTER(c_double),POINTER(c_double),         #eps,idelta
                                POINTER(c_double),POINTER(c_double)]         #mu,temp
    flibs.get_phi_irr.restype=c_void_p
    flibs.get_phi_irr(phi,uni,eig,ffermi,qshift,olist,wlist,byref(c_int64(Nchi)),
                      byref(c_int64(Norb)),byref(c_int64(Nk)),byref(c_int64(Nw)),
                      byref(c_double(idelta)),byref(c_double(eps)),
                      byref(c_double(mu)),byref(c_double(temp)))
    return phi

def phi_qmap(uni,eig,ffermi,klist,olist,Nx,Ny,mu,temp,ecut,idelta,sw_omega):
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
                             POINTER(c_double),POINTER(c_double),                #mu,temp
                             POINTER(c_double),POINTER(c_double),                #,ecut,idelta
                             POINTER(c_double),                                  #eps
                             POINTER(c_int64),POINTER(c_int64),                  #Nx,Ny
                             POINTER(c_int64),POINTER(c_int64),                  #Nk,Norb
                             POINTER(c_int64),POINTER(c_bool)] #Nchi
    flibs.phiq_map.restype=c_void_p
    flibs.phiq_map(phi,uni,eig,ffermi,klist,olist,byref(c_double(mu)),byref(c_double(temp)),
                   byref(c_double(ecut)),byref(c_double(idelta)),byref(c_double(eps)),
                   byref(c_int64(Nx)),byref(c_int64(Ny)),byref(c_int64(Nk)),
                   byref(c_int64(Norb)),byref(c_int64(Nchi)),byref(c_bool(sw_omega)))
    return phi

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

def gen_SCmatrix(olist,site,U,J):
    Nchi=len(olist)
    Smat=np.zeros((Nchi,Nchi),dtype=np.float64)
    Cmat=np.zeros((Nchi,Nchi),dtype=np.float64)
    flibs.get_scmat.argtypes=[np.ctypeslib.ndpointer(dtype=np.float64),             #Smat
                              np.ctypeslib.ndpointer(dtype=np.float64),             #Cmat
                              np.ctypeslib.ndpointer(dtype=np.int64),               #olist
                              np.ctypeslib.ndpointer(dtype=np.int64),               #site
                              POINTER(c_double),POINTER(c_double),POINTER(c_int64)] #U,J,Nchi
    flibs.get_scmat.restype=c_void_p
    flibs.get_scmat(Smat,Cmat,olist,site,byref(c_double(U)),byref(c_double(J)),byref(c_int64(Nchi)))
    return Smat,Cmat

def gen_SCmatrix_orb(olist,site,Umat,Jmat):
    Nchi=len(olist)
    Norb=len(Umat)
    Smat=np.zeros((Nchi,Nchi),dtype=np.float64)
    Cmat=np.zeros((Nchi,Nchi),dtype=np.float64)
    flibs.get_scmat_orb.argtypes=[np.ctypeslib.ndpointer(dtype=np.float64), #Smat
                                  np.ctypeslib.ndpointer(dtype=np.float64), #Cmat
                                  np.ctypeslib.ndpointer(dtype=np.float64), #Umat
                                  np.ctypeslib.ndpointer(dtype=np.float64), #Jmat
                                  np.ctypeslib.ndpointer(dtype=np.int64),   #olist
                                  np.ctypeslib.ndpointer(dtype=np.int64),   #site
                                  POINTER(c_int64),POINTER(c_int64)]        #Nchi,Norb
    flibs.get_scmat_orb.restype=c_void_p
    flibs.get_scmat_orb(Smat,Cmat,olist,site,Umat,Jmat,byref(c_int64(Nchi)),byref(c_int64(Norb)))
    return Smat,Cmat

def calc_Lij(eig,vk,ffermi,mu:float,w:float,idelta:float,temp:float):
    Nk=len(eig)
    Norb=int(eig.size/Nk)
    L11=np.zeros((3,3),dtype=np.complex128)
    L12=np.zeros((3,3),dtype=np.complex128)
    L22=np.zeros((3,3),dtype=np.complex128)
    eps=idelta*1e-3
    flibs.calc_lij.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128), #L11
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

def calc_Kn(eig,veloc,kweight,temp:float,mu:float,tau):
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

def calc_tdf(eig,veloc,kweight,tau,Nw:int):
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

def gen_imp_ham(rvec,ham_r,ham_i,rlist,imp_list,eps=1.0e-5):
    Nr,Nimp,Nsite=len(rvec),len(imp_list),len(rlist)
    Norb=int(np.sqrt(ham_r.size/Nr))
    ham_imp=np.zeros((Norb*Nsite,Norb*Nsite),dtype=np.complex128)
    flibs.gen_imp_ham.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128), #ham_imp
                                np.ctypeslib.ndpointer(dtype=np.complex128), #ham_r
                                np.ctypeslib.ndpointer(dtype=np.float64),    #rvec
                                np.ctypeslib.ndpointer(dtype=np.complex128), #ham_i
                                np.ctypeslib.ndpointer(dtype=np.int64),      #imp_list
                                np.ctypeslib.ndpointer(dtype=np.float64),    #rlist
                                POINTER(c_double),                           #eps
                                POINTER(c_int64),POINTER(c_int64),           #Nimp,Nsite
                                POINTER(c_int64),POINTER(c_int64)]           #Nr,Norb
    flibs.gen_imp_ham.restype=c_void_p
    flibs.gen_imp_ham(ham_imp,ham_r,rvec,ham_i,imp_list,rlist,byref(c_double(eps)),byref(c_int64(Nimp)),
                      byref(c_int64(Nsite)),byref(c_int64(Nr)),byref(c_int64(Norb)))
    return ham_imp

def dft_imp_ham(ham_imp,klist,rlist):
    Nk,Nsite=len(klist),len(rlist)
    Norb=int(len(ham_imp)/Nsite)
    ham_k=np.zeros((Norb*Nk,Norb*Nk),dtype=np.complex128)
    flibs.get_dft_imp_ham.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128), #ham_k
                                    np.ctypeslib.ndpointer(dtype=np.complex128), #ham_imp
                                    np.ctypeslib.ndpointer(dtype=np.float64),    #klist
                                    np.ctypeslib.ndpointer(dtype=np.float64),    #rlist
                                    POINTER(c_int64),POINTER(c_int64),           #Nk,Nsite
                                    POINTER(c_int64)]                            #Norb
    flibs.get_dft_imp_ham.restype=c_void_p
    flibs.get_dft_imp_ham(ham_k,ham_imp,klist,rlist,byref(c_int64(Nk)),
                          byref(c_int64(Nsite)),byref(c_int64(Norb)))
    return ham_k

def get_imp_spectrum(uni,eigs,mu:float,wlist,klist,rlist,eta=1.0e-3):
    Nw,Nk,Nsite=len(wlist),len(klist),len(rlist)
    Norb=int(len(eigs)/Nsite)
    spectrum=np.zeros((Nk,Nw),dtype=np.complex128)
    flibs.get_spectrum_spagehtti.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128), #spectrum
                                           np.ctypeslib.ndpointer(dtype=np.complex128), #uni
                                           np.ctypeslib.ndpointer(dtype=np.float64),    #eigs
                                           np.ctypeslib.ndpointer(dtype=np.float64),    #klist
                                           np.ctypeslib.ndpointer(dtype=np.float64),    #rlist
                                           np.ctypeslib.ndpointer(dtype=np.float64),    #wlist
                                           POINTER(c_int64),POINTER(c_int64),           #Nw,Nk
                                           POINTER(c_int64),POINTER(c_int64),           #Nsite,Norb
                                           POINTER(c_double),POINTER(c_double)]         #mu,eta
    flibs.get_spectrum_spagehtti.retype=c_void_p
    flibs.get_spectrum_spagehtti(spectrum,uni,eigs,klist,rlist,wlist,byref(c_int64(Nw)),byref(c_int64(Nk)),
                                 byref(c_int64(Nsite)),byref(c_int64(Norb)),byref(c_double(mu)),byref(c_double(eta)))
    return spectrum

def get_a(inp_data,xlist):
    Np=len(inp_data)
    a=np.zeros(Np,dtype=np.complex128)
    flibs.get_a_.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128), #a
                           np.ctypeslib.ndpointer(dtype=np.complex128), #xlist
                           np.ctypeslib.ndpointer(dtype=np.complex128), #inpdata
                           POINTER(c_int64)]                            #Np
    flibs.get_a_.retype=c_void_p
    flibs.get_a_(a,xlist,inp_data,byref(c_int64(Np)))
    return a

def get_QP(a,xlist,wlist):
    Nw,Np=len(wlist),len(a)
    Q=np.zeros(Nw,dtype=np.complex128)
    P=np.zeros(Nw,dtype=np.complex128)
    flibs.get_qp_.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128), #P
                            np.ctypeslib.ndpointer(dtype=np.complex128), #Q
                            np.ctypeslib.ndpointer(dtype=np.complex128), #a
                            np.ctypeslib.ndpointer(dtype=np.complex128), #xlist
                            np.ctypeslib.ndpointer(dtype=np.complex128), #wlist
                            POINTER(c_int64),POINTER(c_int64)]           #Nw,Np
    flibs.get_qp_.retype=c_void_p
    flibs.get_qp_(P,Q,a,xlist,wlist,byref(c_int64(Nw)),byref(c_int64(Np)))
    return Q,P

def pade_with_trace(A,iwlist,wlist):
    Nk,Nw,Niw=len(A.T),len(wlist),len(iwlist)
    Norb=int(np.sqrt(A.size/(Nk*Niw)))
    B=np.zeros((Nk,Nw),dtype=np.complex128)
    flibs.pade_with_trace.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128), #A
                                    np.ctypeslib.ndpointer(dtype=np.complex128), #B
                                    np.ctypeslib.ndpointer(dtype=np.complex128), #iwlist
                                    np.ctypeslib.ndpointer(dtype=np.complex128), #wlist
                                    POINTER(c_int64),POINTER(c_int64),           #Nk,Niw
                                    POINTER(c_int64),POINTER(c_int64)]           #Nw,Norb
    flibs.pade_with_trace.retype=c_void_p
    flibs.pade_with_trace(A,B,iwlist,wlist,byref(c_int64(Nk)),byref(c_int64(Niw)),
                          byref(c_int64(Nw)),byref(c_int64(Norb)))
    return B

def get_chi0(Smat,Cmat,Gk,olist,kmap,invk,temp:float,Nx:int,Ny:int,Nz:int):
    Norb,Nchi=len(Gk),len(olist)
    Nkall,Nk,Nw=len(kmap),len(Gk[0,0,0]),len(Gk[0,0])
    chi=np.zeros((Nchi,Nchi,Nw,Nk),dtype=np.complex128)
    flibs.get_chi0.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128), #chi
                             np.ctypeslib.ndpointer(dtype=np.float64),    #Smat
                             np.ctypeslib.ndpointer(dtype=np.float64),    #Cmat
                             np.ctypeslib.ndpointer(dtype=np.complex128), #Gk
                             np.ctypeslib.ndpointer(dtype=np.int64),      #kmap
                             np.ctypeslib.ndpointer(dtype=np.int64),      #invk
                             np.ctypeslib.ndpointer(dtype=np.int64),      #olist
                             POINTER(c_double),POINTER(c_int64),          #temp,Nx
                             POINTER(c_int64),POINTER(c_int64),           #Ny,Nz
                             POINTER(c_int64),POINTER(c_int64),           #Nw,Nk
                             POINTER(c_int64),POINTER(c_int64),           #Nkall,Nchi
                             POINTER(c_int64)]                            #Norb
    flibs.get_chi0.retype=c_void_p
    flibs.get_chi0(chi,Smat,Cmat,Gk,kmap,invk,olist,byref(c_double(temp)),
                   byref(c_int64(Nx)),byref(c_int64(Ny)),byref(c_int64(Nz)),
                   byref(c_int64(Nw)),byref(c_int64(Nk)),byref(c_int64(Nkall)),
                   byref(c_int64(Norb)),byref(c_int64(Nchi)))
    return chi

def linearized_eliashberg(chi,Gk,uni,init_delta,Smat,Cmat,olist,kmap,invk,Nx:int,Ny:int,Nz:int,temp:float,gap_sym:int,eps=1.0e-5,itemax=300):
    Norb,Nchi=len(Gk),len(Smat)
    Nkall,Nk,Nw=len(kmap),len(Gk[0,0,0]),len(Gk[0,0])
    delta=np.zeros((Norb,Norb,Nw,Nkall),dtype=np.complex128)
    flibs.lin_eliash.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128), #delta
                               np.ctypeslib.ndpointer(dtype=np.complex128), #chi
                               np.ctypeslib.ndpointer(dtype=np.complex128), #Gk
                               np.ctypeslib.ndpointer(dtype=np.complex128), #uni
                               np.ctypeslib.ndpointer(dtype=np.float64),    #init_delta
                               np.ctypeslib.ndpointer(dtype=np.float64),    #Smat
                               np.ctypeslib.ndpointer(dtype=np.float64),    #Cmat
                               np.ctypeslib.ndpointer(dtype=np.int64),      #olist
                               np.ctypeslib.ndpointer(dtype=np.int64),      #kmap
                               np.ctypeslib.ndpointer(dtype=np.int64),      #invk
                               POINTER(c_double),POINTER(c_double),         #temp,eps
                               POINTER(c_int64),POINTER(c_int64),           #Nkall,Nk
                               POINTER(c_int64),POINTER(c_int64),           #Nw,Nchi
                               POINTER(c_int64),POINTER(c_int64),           #Norb,Nx
                               POINTER(c_int64),POINTER(c_int64),           #Ny,Nz
                               POINTER(c_int64),POINTER(c_int64)]           #itemax,gapsym
    flibs.lin_eliash.retype=c_void_p
    flibs.lin_eliash(delta,chi,Gk,uni,init_delta,Smat,Cmat,olist,kmap,invk,
                     byref(c_double(temp)),byref(c_double(eps)),byref(c_int64(Nkall)),
                     byref(c_int64(Nk)),byref(c_int64(Nw)),byref(c_int64(Nchi)),
                     byref(c_int64(Norb)),byref(c_int64(Nx)),byref(c_int64(Ny)),
                     byref(c_int64(Nz)),byref(c_int64(itemax)),byref(c_int64(gap_sym)))
    return delta

def linearized_eliashberg_soc(chi,Gk,uni,init_delta,Vmat,sgnsig,sgnsig2,slist,olist,kmap,invk,invs,invschi,
                              Nx:int,Ny:int,Nz:int,temp:float,gap_sym:int,eps=1.0e-4,itemax=300):
    Norb,Nchi=len(slist),len(Vmat)
    Nkall,Nk,Nw=len(kmap),len(Gk[0,0,0]),len(Gk[0,0])
    delta=np.zeros((Norb,Norb,Nw,Nkall),dtype=np.complex128)
    flibs.lin_eliash_soc.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128), #delta
                                   np.ctypeslib.ndpointer(dtype=np.complex128), #chi
                                   np.ctypeslib.ndpointer(dtype=np.complex128), #Gk
                                   np.ctypeslib.ndpointer(dtype=np.complex128), #uni
                                   np.ctypeslib.ndpointer(dtype=np.float64),    #init_delta
                                   np.ctypeslib.ndpointer(dtype=np.float64),    #Vmat
                                   np.ctypeslib.ndpointer(dtype=np.float64),    #sgnsig
                                   np.ctypeslib.ndpointer(dtype=np.float64),    #sgnsig2
                                   np.ctypeslib.ndpointer(dtype=np.int64),      #olist
                                   np.ctypeslib.ndpointer(dtype=np.int64),      #slist
                                   np.ctypeslib.ndpointer(dtype=np.int64),      #kmap
                                   np.ctypeslib.ndpointer(dtype=np.int64),      #invk
                                   np.ctypeslib.ndpointer(dtype=np.int64),      #invs
                                   np.ctypeslib.ndpointer(dtype=np.int64),      #invschi
                                   POINTER(c_double),POINTER(c_double),         #temp,eps
                                   POINTER(c_int64),POINTER(c_int64),           #Nkall,Nk
                                   POINTER(c_int64),POINTER(c_int64),           #Nw,Nchi
                                   POINTER(c_int64),POINTER(c_int64),           #Norb,Nx
                                   POINTER(c_int64),POINTER(c_int64),           #Ny,Nz
                                   POINTER(c_int64),POINTER(c_int64)]           #itemax,gapsym
    flibs.lin_eliash_soc.retype=c_void_p
    flibs.lin_eliash_soc(delta,chi,Gk,uni,init_delta,Vmat,sgnsig,sgnsig2,olist,slist,
                         kmap,invk,invs,invschi,byref(c_double(temp)),byref(c_double(eps)),
                         byref(c_int64(Nkall)),byref(c_int64(Nk)),byref(c_int64(Nw)),
                         byref(c_int64(Nchi)),byref(c_int64(Norb)),byref(c_int64(Nx)),
                         byref(c_int64(Ny)),byref(c_int64(Nz)),byref(c_int64(itemax)),
                         byref(c_int64(gap_sym)))
    return delta

def gen_Vmatrix(olist,slist,site,invs,U,J):
    Nchi,Norb=len(olist),len(slist)
    Vmat=np.zeros((Nchi,Nchi),dtype=np.float64)
    flibs.get_vmat_soc.argtypes=[np.ctypeslib.ndpointer(dtype=np.float64), #Vmat
                                 np.ctypeslib.ndpointer(dtype=np.int64),   #olist
                                 np.ctypeslib.ndpointer(dtype=np.int64),   #slist
                                 np.ctypeslib.ndpointer(dtype=np.int64),   #site
                                 np.ctypeslib.ndpointer(dtype=np.int64),   #invs
                                 POINTER(c_double),POINTER(c_double),      #U,J
                                 POINTER(c_int64),POINTER(c_int64)]        #Nchi,Norb
    flibs.get_vmat_soc.restype=c_void_p
    flibs.get_vmat_soc(Vmat,olist,slist,site,invs,byref(c_double(U)),
                       byref(c_double(J)),byref(c_int64(Nchi)),byref(c_int64(Norb)))
    return Vmat

def gen_Vmatrix_orb(olist,slist,site,invs,Umat,Jmat):
    Nchi,Norb=len(olist),len(Umat)
    Vmat=np.zeros((Nchi,Nchi),dtype=np.float64)
    flibs.get_vmat_soc_orb.argtypes=[np.ctypeslib.ndpointer(dtype=np.float64), #Vmat
                                     np.ctypeslib.ndpointer(dtype=np.int64),   #olist
                                     np.ctypeslib.ndpointer(dtype=np.int64),   #slist
                                     np.ctypeslib.ndpointer(dtype=np.int64),   #site
                                     np.ctypeslib.ndpointer(dtype=np.int64),   #invs
                                     np.ctypeslib.ndpointer(dtype=np.float64), #Umat
                                     np.ctypeslib.ndpointer(dtype=np.float64), #Jmat
                                     POINTER(c_int64),POINTER(c_int64)]        #Nchi,Norb
    flibs.get_vmat_soc_orb.restype=c_void_p
    flibs.get_vmat_soc_orb(Vmat,olist,slist,site,invs,Umat,Jmat,
                           byref(c_int64(Nchi)),byref(c_int64(Norb)))
    return Vmat

def get_chi0_soc(Vmat,Gk,olist,slist,kmap,invk,invs,temp,Nx,Ny,Nz):
    Norb,Nchi=len(slist),len(olist)
    Nkall,Nk,Nw=len(kmap),len(Gk[0,0,0]),len(Gk[0,0])
    chi=np.zeros((Nchi,Nchi,Nw,Nk),dtype=np.complex128)
    sgnsig=np.zeros((Norb,Norb),dtype=np.float64)
    sgnsig2=np.zeros((Nchi,Nchi),dtype=np.float64)
    invschi=np.zeros(Nchi,dtype=np.int64)
    flibs.get_chi0_soc.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128), #chi
                                 np.ctypeslib.ndpointer(dtype=np.float64),    #sgnsig
                                 np.ctypeslib.ndpointer(dtype=np.float64),    #sgnsig2
                                 np.ctypeslib.ndpointer(dtype=np.int64),      #invschi
                                 np.ctypeslib.ndpointer(dtype=np.float64),    #Vmat
                                 np.ctypeslib.ndpointer(dtype=np.complex128), #Gk
                                 np.ctypeslib.ndpointer(dtype=np.int64),      #kmape
                                 np.ctypeslib.ndpointer(dtype=np.int64),      #invk
                                 np.ctypeslib.ndpointer(dtype=np.int64),      #invs
                                 np.ctypeslib.ndpointer(dtype=np.int64),      #olist
                                 np.ctypeslib.ndpointer(dtype=np.int64),      #slist
                                 POINTER(c_double),POINTER(c_int64),         #temp,Nx
                                 POINTER(c_int64),POINTER(c_int64),           #Ny,Nz
                                 POINTER(c_int64),POINTER(c_int64),           #Nw,Nk
                                 POINTER(c_int64),POINTER(c_int64),           #Nkall,Nchi
                                 POINTER(c_int64)]                            #Norb
    flibs.get_chi0_soc.restype=c_void_p
    flibs.get_chi0_soc(chi,sgnsig,sgnsig2,invschi,Vmat,Gk,kmap,invk,invs,olist,slist,
                       byref(c_double(temp)),byref(c_int64(Nx)),byref(c_int64(Ny)),
                       byref(c_int64(Nz)),byref(c_int64(Nw)),byref(c_int64(Nk)),
                       byref(c_int64(Nkall)),byref(c_int64(Nchi)),byref(c_int64(Norb)))
    return chi,sgnsig,sgnsig2,invschi

def get_chis_chic_soc(chi,Vmat,olist,slist,invs):
    def get_orb_list(olist,slist,invs):
        orb_list=np.zeros(Nchi,dtype=np.int64)
        orbs=np.where(slist==1)
        o1,o2=np.meshgrid(orbs,orbs)
        orb_chi=np.array([o1.flatten(),o2.flatten()]).T
        for i,ol in enumerate(olist):
            o1=(ol[0] if slist[ol[0]-1]==1 else invs[ol[0]-1])-1
            o2=(ol[1] if slist[ol[1]-1]==1 else invs[ol[1]-1])-1
            for j,ob2 in enumerate(orb_chi):
                if(o1==ob2[0] and o2==ob2[1]):
                    orb_list[i]=j+1
                    break
        return orb_list
    Norb,Nchi=len(slist),len(olist)
    Nk,Nw=len(chi[0,0,0]),len(chi[0,0])
    orb_list=get_orb_list(olist,slist,invs)
    chic=np.zeros((int(Nchi/4),int(Nchi/4),Nk),dtype=np.complex128)
    chiszz=np.zeros((int(Nchi/4),int(Nchi/4),Nk),dtype=np.complex128)
    chispm=np.zeros((int(Nchi/4),int(Nchi/4),Nk),dtype=np.complex128)
    flibs.get_chis_chic_soc.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128), #chic
                                      np.ctypeslib.ndpointer(dtype=np.complex128), #chiszz
                                      np.ctypeslib.ndpointer(dtype=np.complex128), #chispm
                                      np.ctypeslib.ndpointer(dtype=np.complex128), #chi
                                      np.ctypeslib.ndpointer(dtype=np.float64),    #Vmat
                                      np.ctypeslib.ndpointer(dtype=np.int64),      #orb_list
                                      np.ctypeslib.ndpointer(dtype=np.int64),      #olist
                                      np.ctypeslib.ndpointer(dtype=np.int64),      #slist
                                      np.ctypeslib.ndpointer(dtype=np.int64),      #invs
                                      POINTER(c_int64),POINTER(c_int64),           #Nk,Nw
                                      POINTER(c_int64),POINTER(c_int64)]           #Nchi,Norb
    flibs.get_chis_chic_soc.retype=c_void_p
    flibs.get_chis_chic_soc(chic,chiszz,chispm,chi,Vmat,orb_list,olist,slist,invs,byref(c_int64(Nk)),
                            byref(c_int64(Nw)),byref(c_int64(Nchi)),byref(c_int64(Norb)))
    return chic,chiszz,chispm

def get_chis_chic(chi,Smat,Cmat):
    Nk,Nw,Nchi=len(chi[0,0,0]),len(chi[0,0]),len(Smat)
    chis=np.zeros((Nchi,Nchi,Nk),dtype=np.complex128)
    chic=np.zeros((Nchi,Nchi,Nk),dtype=np.complex128)
    flibs.get_chis_chic.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128),        #chis
                                  np.ctypeslib.ndpointer(dtype=np.complex128),        #chic
                                  np.ctypeslib.ndpointer(dtype=np.complex128),        #chi
                                  np.ctypeslib.ndpointer(dtype=np.float64),           #Smat
                                  np.ctypeslib.ndpointer(dtype=np.float64),           #Cmat
                                  POINTER(c_int64),POINTER(c_int64),POINTER(c_int64)] #Nk,Nw,Nchi
    flibs.get_chis_chic.retype=c_void_p
    flibs.get_chis_chic(chis,chic,chi,Smat,Cmat,byref(c_int64(Nk)),
                        byref(c_int64(Nw)),byref(c_int64(Nchi)))
    return chis,chic

def conv_delta_orb_to_band(delta,uni,invk):
    Nkall,Nk,Nw,Norb=len(invk),len(uni),len(delta[0,0]),len(delta)
    deltab=np.zeros((Norb,Norb,Nkall),dtype=np.complex128)
    flibs.conv_delta_orb_to_band.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128), #deltab
                                           np.ctypeslib.ndpointer(dtype=np.complex128), #delta
                                           np.ctypeslib.ndpointer(dtype=np.complex128), #uni
                                           np.ctypeslib.ndpointer(dtype=np.int64),      #invk
                                           POINTER(c_int64),POINTER(c_int64),           #Norb,Nkall
                                           POINTER(c_int64),POINTER(c_int64)]           #Nk,Nw
    flibs.conv_delta_orb_to_band.retype=c_void_p
    flibs.conv_delta_orb_to_band(deltab,delta,uni,invk,byref(c_int64(Norb)),
                                 byref(c_int64(Nkall)),byref(c_int64(Nk)),byref(c_int64(Nw)))
    return deltab

def conv_delta_orb_to_band_soc(delta,uni,invk,invs,slist):
    Nkall,Nk,Nw,Norb=len(invk),len(uni),len(delta[0,0]),len(delta)
    deltab=np.zeros((Norb,Norb,Nkall),dtype=np.complex128)
    flibs.conv_delta_orb_to_band_soc.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128), #deltab
                                               np.ctypeslib.ndpointer(dtype=np.complex128), #delta
                                               np.ctypeslib.ndpointer(dtype=np.complex128), #uni
                                               np.ctypeslib.ndpointer(dtype=np.int64),      #invk
                                               np.ctypeslib.ndpointer(dtype=np.int64),      #invs
                                               np.ctypeslib.ndpointer(dtype=np.int64),      #slist
                                               POINTER(c_int64),POINTER(c_int64),           #Norb,Nkall
                                               POINTER(c_int64),POINTER(c_int64)]           #Nk,Nw
    flibs.conv_delta_orb_to_band_soc.retype=c_void_p
    flibs.conv_delta_orb_to_band_soc(deltab,delta,uni,invk,invs,slist,byref(c_int64(Norb)),
                                     byref(c_int64(Nkall)),byref(c_int64(Nk)),byref(c_int64(Nw)))
    return deltab

def gen_Fk(Gk,delta,invk):
    Nkall,Nk,Nw,Norb=len(invk),len(Gk[0,0,0]),len(delta[0,0]),len(delta)
    Fk=np.zeros((Norb,Norb,Nw,Nkall),dtype=np.complex128)
    flibs.mkfk_trs_nsoc_.argtypes=[np.ctypeslib.ndpointer(dtype=np.complex128), #Fk
                                   np.ctypeslib.ndpointer(dtype=np.complex128), #Gk
                                   np.ctypeslib.ndpointer(dtype=np.complex128), #delta
                                   np.ctypeslib.ndpointer(dtype=np.int64),      #invk
                                   POINTER(c_int64),POINTER(c_int64),           #Nkall,Nk
                                   POINTER(c_int64),POINTER(c_int64)]           #Nw,Norb
    flibs.mkfk_trs_nsoc_.retype=c_void_p
    flibs.mkfk_trs_nsoc_(Fk,Gk,delta,invk,byref(c_int64(Nkall)),byref(c_int64(Nk)),
                         byref(c_int64(Nw)),byref(c_int64(Norb)))

    return Fk

def gen_irr_k_TRS(Nx,Ny,Nz):
    Nkall=Nx*Ny*Nz
    if(Nkall%2==0):
        if(Nz%2==0): #all even Nk=Nkall/2+4 else Nz even Nk=Nkall+2
            Nk=int(Nkall/2+4) if((Nx%2)==0 and (Ny%2)==0) else (int(Nkall/2+2) if((Nx%2)==0 or (Ny%2)==0) else int(Nkall/2+1))
        else: #NxNy plane even Nk=Nkall/2+2 else Nk=Nkall/2+1
            Nk=int(Nkall/2+2) if((Nx%2)==0 and (Ny%2)==0) else int(Nkall/2+1)
    else: #all odd Nk=(Nkall+1)/2
        Nk=int((Nkall+1)/2)
    klist=np.zeros((Nk,3),dtype=np.float64)
    kmap=np.zeros((Nkall,3),dtype=np.int64)
    invk_ft_list=np.zeros((Nkall,3),dtype=np.int64)
    flibs.generate_irr_kpoint_inv.argtypes=[np.ctypeslib.ndpointer(dtype=np.float64), #klist
                                            np.ctypeslib.ndpointer(dtype=np.int64),   #kmap
                                            np.ctypeslib.ndpointer(dtype=np.int64),   #invk_ft_list
                                            POINTER(c_int64),POINTER(c_int64),        #Nk,Nx
                                            POINTER(c_int64),POINTER(c_int64)]        #Ny,Nz
    flibs.generate_irr_kpoint_inv.retype=c_void_p
    flibs.generate_irr_kpoint_inv(klist,kmap,invk_ft_list,byref(c_int64(Nk)),
                                  byref(c_int64(Nx)),byref(c_int64(Ny)),byref(c_int64(Nz)))
    return klist,kmap,invk_ft_list
