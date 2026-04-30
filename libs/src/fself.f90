! Caches forward/backward FFTW plans across calls within a loop.
! Call init_fft_plans before the orbital loop and destroy_fft_plans after.
module fftw_plan_cache
  use,intrinsic::iso_c_binding, only:c_int64_t
  implicit none
  integer(c_int64_t),save :: plan_fwd=0, plan_bwd=0
end module fftw_plan_cache

subroutine init_fft_plans(cmat,tmp,Nx,Ny,Nz,Nw)
  use fftw_plan_cache
  use,intrinsic::iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nx,Ny,Nz,Nw
  complex(c_double),intent(inout),dimension(Nx,Ny,Nz,Nw):: cmat,tmp
  integer(c_int32_t),dimension(4):: Nlist
  Nlist=(/int(Nx,c_int32_t),int(Ny,c_int32_t),int(Nz,c_int32_t),int(Nw,c_int32_t)/)
  call dfftw_plan_dft(plan_fwd,4,Nlist,cmat,tmp,-1,64)
  call dfftw_plan_dft(plan_bwd,4,Nlist,cmat,tmp, 1,64)
end subroutine init_fft_plans

subroutine destroy_fft_plans()
  use fftw_plan_cache
  implicit none
  call dfftw_destroy_plan(plan_fwd)
  call dfftw_destroy_plan(plan_bwd)
  plan_fwd=0
  plan_bwd=0
end subroutine destroy_fft_plans

subroutine FFT(cmat,tmp,Nx,Ny,Nz,Nw,SW)
  !> 4D Fast Fourier transform function
  !!@param cmat,inout: FFT source and results
  !!@param  tmp,inout: temporary
  !!@param      Nx,in: x mesh
  !!@param      Ny,in: y mesh
  !!@param      Nz,in: z mesh
  !!@param      Nw,in: w mesh
  !!@param      SW,in: FFT or IFFT switch (if true IFFT)
  use fftw_plan_cache
  use,intrinsic::iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nx,Ny,Nz,Nw
  logical,intent(in):: SW
  complex(c_double),intent(inout),dimension(Nx,Ny,Nz,Nw):: cmat,tmp

  integer(c_int64_t) plan
  integer(c_int32_t) Inv
  integer(c_int32_t),dimension(4):: Nlist

  Nlist=(/int(Nx,c_int32_t),int(Ny,c_int32_t),int(Nz,c_int32_t),int(Nw,c_int32_t)/)
  if(plan_fwd /= 0)then
     ! Reuse cached plans: avoids heap alloc/free on every call
     if(SW)then
        call dfftw_execute_dft(plan_fwd,cmat,tmp)
     else
        call dfftw_execute_dft(plan_bwd,cmat,tmp)
     end if
  else
     if(SW)then
        Inv=-1
     else
        Inv=1
     end if
     call dfftw_plan_dft(plan,4,Nlist,cmat,tmp,Inv,64)
     call dfftw_execute(plan)
     call dfftw_destroy_plan(plan)
  end if
  if(.not. SW)then
     cmat(:,:,:,:)=tmp(:,:,:,:)/product(Nlist)
  else
     cmat(:,:,:,:)=tmp(:,:,:,:)
  end if
end subroutine FFT

subroutine gen_green0(Gk,eig,uni,mu,temp,Nk,Nw,Norb) bind(C,name="gen_green0_")
  !> This function obtains green function G0
  !!@param  Gk,out: green function
  !!@param  eig,in: energies at k-points
  !!@param  uni,in: unitary matrix
  !!@param   mu,in: chemical potential
  !!@param temp,in: Temperature
  !!@param   Nk,in: The number of k-points
  !!@param   Nw,in: The number of energies mesh
  !!@param Norb,in: The number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  use constant
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,norb
  real(c_double),intent(in):: mu,temp
  real(c_double),intent(in),dimension(Norb,Nk):: eig
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(c_double),intent(out),dimension(Nk,Nw,Norb,Norb):: Gk
  
  integer(c_int32_t) i,j,l,m,n
  complex(c_double) iw

   Gk(:,:,:,:)=0.0d0

  !$omp parallel private(l,m)
  do l=1,Norb
     do m=1,Norb
        !$omp do private(iw,i,j,n)
        wloop: do j=1,Nw !ien=pi(2l+1)/beta l=0,1,... beta=(kBT)^-1
           iw=cmplx(mu,dble(2*j-1)*pi*temp,kind=c_double) !iω_j = μ + i(2j-1)πT (j=1: ω_0=πT)
           band_loop: do n=1,Norb
              kloop: do i=1,Nk
                 Gk(i,j,m,l)=Gk(i,j,m,l)+uni(m,n,i)*conjg(uni(l,n,i))/(iw-eig(n,i))
              end do kloop
           end do band_loop
        end do wloop
        !$omp end do
     end do
  end do
  !$omp end parallel
end subroutine gen_green0

subroutine gen_green_inv(Gk,self,hamk,mu,temp,Nk,Nw,Norb) bind(C,name="gen_green_inv_")
  !> This function obtains inverse of green function G^-1
  !!@param  Gk,out: inverse of green function
  !!@param self,in: self energies
  !!@param hamk,in: hamiltonian at k-point
  !!@param   mu,in: chemical potential
  !!@param temp,in: Temperature
  !!@param   Nk,in: The number of k-points
  !!@param   Nw,in: The number of energies mesh
  !!@param Norb,in: The number of orbitals
  use,intrinsic:: iso_c_binding, only: c_int64_t,c_double,c_int32_t
  use constant
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,Norb
  real(c_double),intent(in):: mu,temp
  complex(c_double),intent(in),dimension(norb,norb,Nk):: hamk
  complex(c_double),intent(in),dimension(Nk,Nw,norb,norb):: self
  complex(c_double),intent(out),dimension(Nk,Nw,norb,norb):: Gk

  integer(c_int32_t)i,j,l,m
  complex(c_double) iw

  !G^-1=G^-1_0-sigma (=iwI-Hk-sigma)
  !$omp parallel private(l,m)
  do l=1,Norb
     do m=1,Norb
        !$omp do private(iw,i)
        do j=1,Nw
           iw=cmplx(mu,dble(2*j-1)*pi*temp,kind=c_double)
           do i=1,Nk
              Gk(i,j,m,l)=-(hamk(m,l,i)+self(i,j,m,l))
              if(l==m)then
                 Gk(i,j,l,l)=Gk(i,j,l,l)+iw
              end if
           end do
        end do
        !$omp end do
     end do
  end do
  !$omp end parallel
end subroutine gen_green_inv

subroutine gen_green_inv_from_eig(Gk,self,uni,eig,mu,temp,Nk,Nw,Norb) bind(C)
  !> This function obtains inverse of green function G^-1 form energy and unitary matrix
  !!@param  Gk,out: inverse of green function
  !!@param self,in: self energies
  !!@param  uni,in: unitary matrix
  !!@param  eig,in: energies at k-points
  !!@param   mu,in: chemical potential
  !!@param temp,in: Temperature
  !!@param   Nk,in: The number of k-points
  !!@param   Nw,in: The number of energies mesh
  !!@param Norb,in: The number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  use constant
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,norb
  real(c_double),intent(in):: mu,temp
  real(c_double),intent(in),dimension(norb,Nk):: eig
  complex(c_double),intent(in),dimension(norb,norb,Nk):: uni
  complex(c_double),intent(in),dimension(Nk,Nw,norb,norb):: self
  complex(c_double),intent(out),dimension(Nk,Nw,norb,norb):: Gk

  integer(c_int32_t)i,j,l,m,n
  complex(c_double) iw

  !$omp parallel private(l,m)
  !$omp workshare
  Gk(:,:,:,:)=0.0d0
  !$omp end workshare
  do l=1,norb
     do m=1,norb
        !$omp do private(iw,i,n)
        wloop: do j=1,Nw
           iw=cmplx(mu,dble(2*j-1)*pi*temp,kind=c_double) !(2(j-1)+1)pi/beta
           kloop: do i=1,Nk
              band_loop: do n=1,norb
                 Gk(i,j,m,l)=Gk(i,j,m,l)+uni(m,n,i)*conjg(uni(l,n,i))*(iw-eig(n,i))
              end do band_loop
              Gk(i,j,m,l)=Gk(i,j,m,l)-self(i,j,m,l)
           end do kloop
        end do wloop
        !$omp end do
     end do
  end do
  !$omp end parallel
end subroutine gen_green_inv_from_eig

subroutine getinv(Gk,Nk,Nw,Norb) bind(C,name="getinv_")
  !> This function obtain green function form G^-1
  !!@param Gk,inout: green function
  !!@param    Nk,in: The number of k-points
  !!@param    Nw,in: The number of energies mesh
  !!@param  Norb,in: The number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,Norb
  complex(c_double),intent(inout),dimension(Nk,Nw,Norb,Norb):: Gk

  integer(c_int32_t) i,j,ipiv(Norb),info
  complex(c_double) tmp(Norb,Norb),work(2*Norb)

  !$omp parallel do private(i,j,tmp,work,ipiv,info)
  do i=1,Nw
     do j=1,Nk
        tmp(:,:)=Gk(j,i,:,:)
        call zgetrf(Norb,Norb,tmp,Norb,ipiv,info)
        if(info/=0)then; print*,'zgetrf failed in getinv: info=',info; stop; end if
        call zgetri(Norb,tmp,Norb,ipiv,work,2*Norb,info)
        if(info/=0)then; print*,'zgetri failed in getinv: info=',info; stop; end if
        Gk(j,i,:,:)=tmp(:,:)
     end do
  end do
  !$omp end parallel do
end subroutine getinv

subroutine ckchi_impl(chi,Smat,Cmat,kmap,invk,Nk,Nkall,Nchi,Nw,maxchi0s_out)
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nkall,Nchi,Nw
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
  real(c_double),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  real(c_double),intent(out):: maxchi0s_out
  complex(c_double),intent(in),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(c_int32_t) i,info,chisk,chick,chiskall,chickall
  real(c_double) maxchi0s,maxchi0c,maxchi0s2,maxchi0c2
  real(c_double),dimension(2*Nchi):: rwork
  complex(c_double),dimension(Nchi*Nchi*4+1):: work
  complex(c_double),dimension(Nchi):: eigs,eigc
  complex(c_double),dimension(Nchi,Nchi):: chi0s,chi0c,tmp1,tmp2,Smat_c,Cmat_c,chi_tmp

  Smat_c = cmplx(Smat, 0.0d0, kind=c_double)
  Cmat_c = cmplx(Cmat, 0.0d0, kind=c_double)
  maxchi0s2=-1.0d5
  maxchi0c2=-1.0d5
  do i=1,Nk
     chi_tmp = chi(i,1,:,:)
     call zgemm('N','N',Nchi,Nchi,Nchi,(1.0d0,0.0d0),chi_tmp,Nchi,Smat_c,Nchi,(0.0d0,0.0d0),chi0s,Nchi)
     call zgemm('N','N',Nchi,Nchi,Nchi,(-1.0d0,0.0d0),chi_tmp,Nchi,Cmat_c,Nchi,(0.0d0,0.0d0),chi0c,Nchi)
     call zgeev('N','N',Nchi,chi0s,Nchi,eigs,tmp1,Nchi,tmp2,Nchi,work,Nchi*Nchi*4+1,rwork,info)
     call zgeev('N','N',Nchi,chi0c,Nchi,eigc,tmp1,Nchi,tmp2,Nchi,work,Nchi*Nchi*4+1,rwork,info)
     maxchi0s=maxval(dble(eigs))
     maxchi0c=maxval(dble(eigc))
     if(maxchi0s>maxchi0s2)then
        chisk=i
        maxchi0s2=maxchi0s
     end if
     if(maxchi0c>maxchi0c2)then
        chick=i
        maxchi0c2=maxchi0c
     end if
  end do
  do i=1,Nkall !get kmap footnote
     if(invk(2,i)==0)then
        if(invk(1,i)==chisk)chiskall=i
        if(invk(1,i)==chick)chickall=i
     end if
  end do
  print'(A3,3I4,F12.8)','SDW',kmap(:,chiskall),maxchi0s2
  print'(A3,3I4,F12.8)','CDW',kmap(:,chickall),maxchi0c2
  maxchi0s_out=maxchi0s2
end subroutine ckchi_impl

subroutine get_chi0(chi,Smat,Cmat,Gk,kmap,invk,olist,temp,Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi) bind(C)
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
  real(c_double),intent(in):: temp
  real(c_double),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(c_double),intent(out),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(c_int32_t),dimension(Nchi,Nchi,2)::chi_map
  integer(c_int32_t),dimension(Nchi*(Nchi+1)/2,2)::irr_chi
  real(c_double) dummy_maxchi0s

  call get_chi_map(chi_map,irr_chi,olist,Nchi)
  call get_chi0_conv(chi,Gk,kmap,invk,irr_chi,chi_map,olist,temp,Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi)
  call ckchi_impl(chi,Smat,Cmat,kmap,invk,Nk,Nkall,Nchi,Nw,dummy_maxchi0s)
end subroutine get_chi0

subroutine get_chi_map(chi_map,irr_chi,olist,Nchi)
  !> This function generate index of exchange symmetry chi1234(q,iw)=chi*4321(q,iw).
  !> This symmetry can use system has no spin dependence and TRS.
  !!@param chi_map,out: mapping list of chi index
  !!@param irr_chi,out: irreducible index of chii
  !!@param    olist,in: orbital index list of chi index
  !!@param     Nchi,in: Number of chi index
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nchi
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int32_t),intent(out),dimension(Nchi,Nchi,2):: chi_map
  integer(c_int32_t),intent(out),dimension(Nchi*(Nchi+1)/2,2):: irr_chi

  integer(c_int32_t) l1,m1,l2,m2,iter
  integer(c_int32_t),dimension(Nchi,Nchi):: chi_irr
  integer(c_int32_t),dimension(4):: tmp1,tmp2
  logical ck

  chi_irr(:,:)=0
  chi_map(:,:,:)=0
  do l1=1,Nchi
     tmp1(1)=olist(l1,1) !1 of 1234
     tmp1(2)=olist(l1,2) !2 of 1234
     do m1=1,Nchi
        tmp1(3)=olist(m1,1) !3 of 1234
        tmp1(4)=olist(m1,2) !4 of 1234
        ck=.false.
        do l2=1,Nchi
           tmp2(3)=olist(l2,2) !2 of 4321
           tmp2(4)=olist(l2,1) !1 of 4321
           do m2=1,Nchi
              tmp2(1)=olist(m2,2) !4 of 4321
              tmp2(2)=olist(m2,1) !3 of 4321
              if(sum(abs(tmp1(:)-tmp2(:)))==0)then !4321 correspond to 1234
                 chi_map(m1,l1,1)=m2
                 chi_map(m1,l1,2)=l2
                 if(chi_irr(m2,l2)==0)then !4321 is not irreducible
                    chi_irr(m1,l1)=1 !1 is irreducible index
                 end if
                 ck=.true.
                 exit
              end if
           end do
           if(ck)exit
        end do
     end do
  end do

  !get list of irreducible chi index
  iter=1
  do l1=1,Nchi
     do m1=1,Nchi
        if(chi_irr(m1,l1)==1)then
           irr_chi(iter,1)=l1
           irr_chi(iter,2)=m1
           iter=iter+1
        end if
     end do
  end do
end subroutine get_chi_map

subroutine get_chi0_conv(chi,Gk,kmap,invk,irr_chi,chi_map,olist,temp,&
     Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi) bind(C,name='get_chi0_conv_')
  !> This function obtains chi_0 using convolution.
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nw,Norb,Nchi,Nkall,Nk,Nx,Ny,Nz
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
  integer(c_int32_t),intent(in),dimension(Nchi,Nchi,2):: chi_map
  integer(c_int32_t),intent(in),dimension(Nchi*(Nchi+1)/2,2):: irr_chi
  real(c_double),intent(in):: temp
  complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(c_double),intent(out),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(c_int32_t) i,j,k,l,m,n,iorb
  integer(c_int32_t) ii(0:Nx-1),ij(0:Ny-1),ik(0:Nz-1),iw(2*Nw)
  real(c_double) weight
  real(c_double),parameter:: eps=1.0d-9
  complex(c_double),dimension(0:Nx-1,0:Ny-1,0:Nz-1,2*Nw):: tmp,tmpgk13,tmpgk42
  
  weight=temp/dble(Nkall)
  ii(0)=0
  ij(0)=0
  ik(0)=0
  iw(1)=1
  !$omp parallel
  !$omp do
  do i=1,Nx-1
     ii(i)=Nx-i
  end do
  !$omp end do
  !$omp do
  do i=1,Ny-1
     ij(i)=Ny-i
  end do
  !$omp end do
  !$omp do
  do i=1,Nz-1
     ik(i)=Nz-i
  end do
  !$omp end do
  !$omp do
  do i=2,2*Nw
     iw(i)=2*Nw-i+2
  end do
  !$omp end do
  !$omp end parallel
  call init_fft_plans(tmpgk13,tmp,Nx,Ny,Nz,2*Nw)
  orb_loop:do iorb=1,Nchi*(Nchi+1)/2
     l=irr_chi(iorb,1)
     m=irr_chi(iorb,2)
     !use symmetry G^lm(k,iw)=G^ml(k,-iw) from Hermitian symmetry of Hamiltonian
     !$omp parallel do private(i)
     w_loop_Gk_to_tmp:do j=1,Nw
        k_loop_Gk_to_tmp:do i=1,Nkall
           if(invk(2,i)==0)then
              !iw
              tmpgk13(kmap(1,i),kmap(2,i),kmap(3,i),j)=Gk(invk(1,i),j,olist(l,1),olist(m,1)) !G13(k,iw)
              tmpgk42(kmap(1,i),kmap(2,i),kmap(3,i),j)=Gk(invk(1,i),j,olist(m,2),olist(l,2)) !G42(k,iw)
              !-iw
              tmpgk13(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=conjg(Gk(invk(1,i),j,olist(m,1),olist(l,1))) !G13(k,-iw)=G^*31(k,iw)
              tmpgk42(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=conjg(Gk(invk(1,i),j,olist(l,2),olist(m,2))) !G42(k,-iw)=G^*24(k,iw)
           else if(invk(2,i)==1)then
              !iw
              tmpgk13(kmap(1,i),kmap(2,i),kmap(3,i),j)=Gk(invk(1,i),j,olist(m,1),olist(l,1)) !G13(-k,iw)=G^31(k,iw)
              tmpgk42(kmap(1,i),kmap(2,i),kmap(3,i),j)=Gk(invk(1,i),j,olist(l,2),olist(m,2)) !G42(-k,iw)=G^24(k,iw)
              !-iw
              tmpgk13(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=conjg(Gk(invk(1,i),j,olist(l,1),olist(m,1))) !G13(-k,-iw)=G^*13(k,iw)
              tmpgk42(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=conjg(Gk(invk(1,i),j,olist(m,2),olist(l,2))) !G42(-k,-iw)=G^*42(k,iw)
           end if
        end do k_loop_Gk_to_tmp
     end do w_loop_Gk_to_tmp
     !$omp end parallel do
     call FFT(tmpgk13,tmp,Nx,Ny,Nz,2*Nw,.true.)
     call FFT(tmpgk42,tmp,Nx,Ny,Nz,2*Nw,.true.)
     !calculate G(r)G(-r)
     !$omp parallel do private(i,j,k)
     w_loop_conv:do n=1,2*Nw
        z_loop:do k=0,Nz-1
           y_loop:do j=0,Ny-1
              !$omp simd
              x_loop:do i=0,Nx-1
                 tmp(i,j,k,n)=-tmpgk13(i,j,k,n)*tmpgk42(ii(i),ij(j),ik(k),iw(n))
              end do x_loop
              !$omp end simd
           end do y_loop
        end do z_loop
     end do w_loop_conv
     !$omp end parallel do
     call FFT(tmp,tmpgk13,Nx,Ny,Nz,2*Nw,.false.)
     !$omp parallel do private(i,j)
     w_loop_tmp_to_chi:do j=1,Nw
        k_loop_tmp_to_chi:do i=1,Nkall
           if(invk(2,i)==0)then
              chi(invk(1,i),j,m,l)=tmp(kmap(1,i),kmap(2,i),kmap(3,i),j)*weight
              if(abs(dble(chi(invk(1,i),j,m,l)))<eps) chi(invk(1,i),j,m,l)=cmplx(0.0d0,imag(chi(invk(1,i),j,m,l)),kind=c_double)
              if(abs(imag(chi(invk(1,i),j,m,l)))<eps) chi(invk(1,i),j,m,l)=cmplx(dble(chi(invk(1,i),j,m,l)),0.0d0,kind=c_double)
           end if
        end do k_loop_tmp_to_chi
     end do w_loop_tmp_to_chi
     !$omp end parallel do
     chi(:,:,chi_map(m,l,1),chi_map(m,l,2))=conjg(chi(:,:,m,l))
  end do orb_loop
  call destroy_fft_plans()
end subroutine get_chi0_conv

subroutine get_chi0_sum(chi,Gk,klist,invk,irr_chi,chi_map,olist,temp,Nw,Nk,Nkall,Norb,Nchi) bind(C)
  !> It obtains chi_0 using summation. Its cost is O(Nk^2), so it is heavy. You should use get_chi0_conv.
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nw,Norb,Nchi,Nk,Nkall
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),dimension(3,Nkall):: invk
  integer(c_int32_t),intent(in),dimension(Nchi,Nchi,2):: chi_map
  integer(c_int32_t),intent(in),dimension(Nchi*(Nchi+1)/2,2):: irr_chi
  real(c_double),intent(in),dimension(3,Nkall):: klist
  real(c_double),intent(in):: temp
  complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(c_double),intent(out),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(c_int32_t) i,j,l,m,iq,iw,iorb
  integer(c_int32_t),dimension(2*Nw):: wshift
  integer(c_int64_t),dimension(Nkall):: qshift
  real(c_double) weight
  complex(c_double),dimension(Nkall,2*Nw):: tmpgk13,tmpgk42

  !$omp parallel workshare
  chi(:,:,:,:)=0.0d0
  !$omp end parallel workshare

  weight=temp/dble(Nkall)
  orb_loop:do iorb=1,Nchi*(Nchi+1)/2
     l=irr_chi(iorb,1)
     m=irr_chi(iorb,2)
     !$omp parallel do private(i)
     do j=1,Nw
        do i=1,Nkall
           if(invk(2,i)==0)then
              tmpgk13(i,j)=Gk(i,j,olist(m,1),olist(l,1))
              tmpgk42(i,j)=Gk(i,j,olist(l,2),olist(m,2))
              tmpgk13(i,2*Nw-j+1)=conjg(Gk(i,j,olist(l,1),olist(m,1)))
              tmpgk42(i,2*Nw-j+1)=conjg(Gk(i,j,olist(m,2),olist(l,2)))
           else if(invk(2,i)==1)then
              tmpgk13(i,j)=Gk(i,j,olist(l,1),olist(m,1))
              tmpgk42(i,j)=Gk(i,j,olist(m,2),olist(l,2))
              tmpgk13(i,2*Nw-j+1)=conjg(Gk(i,j,olist(m,1),olist(l,1)))
              tmpgk42(i,2*Nw-j+1)=conjg(Gk(i,j,olist(l,2),olist(m,2)))
           end if
        end do
     end do
     !$omp end parallel do
     wloop: do iw=1,Nw
        wl_loop: do j=1,2*Nw
           wshift(j)=mod(j+iw-1,2*Nw)
           if(wshift(j)==0)then
              wshift(j)=2*Nw
           end if
        end do wl_loop
        qloop: do iq=1,Nkall
           if(invk(2,iq)==0)then
              call get_qshift(klist(:,iq),klist,qshift,Nk)
              !$omp parallel do private(i) reduction(+: chi)
              wloop2: do j=1,2*Nw
                 !$omp simd
                 kloop: do i=1,Nkall
                    chi(invk(1,iq),iw,m,l)=chi(invk(1,iq),iw,m,l)-tmpgk13(i,j)*tmpgk42(qshift(i),wshift(j))
                 end do kloop
                 !$omp end simd
              end do wloop2
              !$omp end parallel do
           end if
        end do qloop
     end do wloop
     chi(:,:,chi_map(m,l,1),chi_map(m,l,2))=conjg(chi(:,:,m,l))
  end do orb_loop
  chi(:,:,:,:)=chi(:,:,:,:)*weight
end subroutine get_chi0_sum

subroutine get_vsigma_flex_nosoc(chi,Smat,Cmat,Nk,Nw,Nchi) bind(C,name='get_vsigma_flex_nosoc_')
  !> This function obtains interaction V_sigma without soc
  !!@param  chi,inout: irreducible susceptibility and pairing interaction
  !!@param    Smat,in: S-matrix
  !!@param    Cmat,in: C-matrix
  !!@param      Nk,in: Number of k-points
  !!@param      Nw,in: Number of Matsubara frequencies
  !!@param    Nchi,in: Number of footnote of chi
   use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,Nchi
  real(c_double),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  complex(c_double),intent(inout),dimension(Nk,Nw,Nchi,Nchi):: chi

   integer(c_int32_t) i,j,l,info
  integer(c_int32_t),dimension(Nchi):: ipiv
   complex(c_double),dimension(Nchi,Nchi):: cmat1,cmat2,cmat3,cmat4,cmat5,Smat_c,Cmat_c,SC_c

   Smat_c(:,:)=cmplx(Smat(:,:),0.0d0,kind=c_double)
   Cmat_c(:,:)=cmplx(Cmat(:,:),0.0d0,kind=c_double)
   SC_c(:,:)=Smat_c(:,:)+Cmat_c(:,:)

  !$omp parallel do collapse(2) private(i,cmat1,cmat2,cmat3,cmat4,cmat5,ipiv,info,l)
  do j=1,Nw
     do i=1,Nk
        call zgemm('N','N',Nchi,Nchi,Nchi,(-1.0d0,0.0d0),chi(i,j,:,:),Nchi,Smat_c,Nchi,(0.0d0,0.0d0),cmat1,Nchi) !-chi0S
        call zgemm('N','N',Nchi,Nchi,Nchi,(1.0d0,0.0d0),chi(i,j,:,:),Nchi,Cmat_c,Nchi,(0.0d0,0.0d0),cmat2,Nchi)  !chi0C
        cmat3(:,:)=-cmat1(:,:) !chi0S (RHS for first solve)
        cmat4(:,:)=cmat2(:,:)  !chi0C (RHS for second solve)
        cmat5(:,:)=cmat3(:,:)+cmat4(:,:) !chi0S+chi0C (save before zgesv overwrites)
        do l=1,Nchi
           cmat1(l,l)=cmat1(l,l)+1.0d0 !I-chi0S
           cmat2(l,l)=cmat2(l,l)+1.0d0 !I+chi0C
        end do
        call zgesv(Nchi,Nchi,cmat1,Nchi,ipiv,cmat3,Nchi,info) !(I-chi0S)X=chi0S -> cmat3=chiS
        if(info/=0)then; print*,'zgesv failed: info=',info; stop; end if
        call zgesv(Nchi,Nchi,cmat2,Nchi,ipiv,cmat4,Nchi,info) !(I+chi0C)X=chi0C -> cmat4=chiC
        if(info/=0)then; print*,'zgesv failed: info=',info; stop; end if
        ! FLEX self-energy kernel: V_σ = 3/2·S·χ_s + 1/2·C·χ_c - 1/4·(S+C)·(χ_s+χ_c) + static (3/2·S - 1/2·C)
        cmat1(:,:)=1.5d0*Smat(:,:)-0.5d0*Cmat(:,:)   !static (HF-like) bare vertex
        call zgemm('N','N',Nchi,Nchi,Nchi,(1.5d0,0.0d0),Smat_c,Nchi,cmat3,Nchi,(1.0d0,0.0d0),cmat1,Nchi)  !+3/2·S·χ_s
        call zgemm('N','N',Nchi,Nchi,Nchi,(0.5d0,0.0d0),Cmat_c,Nchi,cmat4,Nchi,(1.0d0,0.0d0),cmat1,Nchi)  !+1/2·C·χ_c
        call zgemm('N','N',Nchi,Nchi,Nchi,(-0.25d0,0.0d0),SC_c,Nchi,cmat5,Nchi,(1.0d0,0.0d0),cmat1,Nchi) !-1/4·(S+C)·(χ_s+χ_c), subtract double count
        chi(i,j,:,:)=cmat1(:,:)
     end do
  end do
  !$omp end parallel do
end subroutine get_vsigma_flex_nosoc

subroutine calc_sigma(sigmak,Gk,Vsigma,Smat,Cmat,kmap,invk,olist,temp,Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz)
  !> This function obtain new gap function Delta=sum VF
  !!@param sigmak,out: self energy
  !!@param      Gk,in: green function
  !!@param  Vsigma,in: irreducible susceptibility and interaction
  !!@param    Smat,in: S-matrix
  !!@param    Cmat,in: C-matrix
  !!@param    kmap,in: property of k-point
  !!@param    invk,in: list of reverse k-points
  !!@param   olist,in: property of chi footnote
  !!@param    temp,in: Temperature
  !!@param   Nkall,in: Number of all k-points
  !!@param      Nk,in: Number of k-points
  !!@param      Nw,in: Number of Matsubara frequencies
  !!@param    Nchi,in: Number of footnote of chi
  !!@param    Norb,in: Number of orbitals
  !!@param      Nx,in: Number of kx mesh
  !!@param      Ny,in: Number of ky mesh
  !!@param      Nz,in: Number of kz mesh
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  real(c_double),intent(in):: temp
  real(c_double),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  complex(c_double),intent(in),dimension(Nk,Nw,Nchi,Nchi):: Vsigma
  complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(c_double),intent(out),dimension(Nk,Nw,Norb,Norb):: sigmak

  integer(c_int32_t) i,j,k,n,l,m
  real(c_double) weight
  real(c_double),parameter:: eps=1.0d-9
  complex(c_double),dimension(0:Nx-1,0:Ny-1,0:Nz-1,2*Nw):: tmpVsigma,tmp,tmpgk

  weight=temp/dble(Nkall)
  !$omp parallel workshare
  sigmak(:,:,:,:)=0.0d0
  !$omp end parallel workshare
  call init_fft_plans(tmpVsigma,tmp,Nx,Ny,Nz,2*Nw)
  do l=1,Nchi
     do m=1,Nchi
        !$omp parallel
        !$omp do
        do i=1,Nkall !j=1>iw=0, j=Nw>iw=inf
           if(invk(2,i)==0)then
              tmpVsigma(kmap(1,i),kmap(2,i),kmap(3,i),1)=Vsigma(invk(1,i),1,m,l) !Vsigma(k,iw)
           else if(invk(2,i)==1)then
              tmpVsigma(kmap(1,i),kmap(2,i),kmap(3,i),1)=Vsigma(invk(1,i),1,l,m) !Vsigma(-k,iw)
           end if
           tmpVsigma(kmap(1,i),kmap(2,i),kmap(3,i),Nw+1)=1.5d0*Smat(m,l)-0.5d0*Cmat(m,l)
        end do
        !$omp end do
        !$omp do private(i)
        do j=2,Nw
           do i=1,Nkall
              if(invk(2,i)==0)then
                 tmpVsigma(kmap(1,i),kmap(2,i),kmap(3,i),j)=Vsigma(invk(1,i),j,m,l) !Vsigma(k,iw)
                 tmpVsigma(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+2)=conjg(Vsigma(invk(1,i),j,l,m)) !Vsigma(k,-iw)ml=Vsigma^*lm(k,iw)
              else if(invk(2,i)==1)then
                 tmpVsigma(kmap(1,i),kmap(2,i),kmap(3,i),j)=Vsigma(invk(1,i),j,l,m) !Vsigma(-k,iw)=Vsigma(-k,iw)^T
                 tmpVsigma(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+2)=conjg(Vsigma(invk(1,i),j,m,l)) !Vsigma(-k,-iw)=Vsigma^*(k,iw)
              end if
           end do
        end do
        !$omp end do
        !$omp do private(i)
        do j=1,Nw
           do i=1,Nkall
              if(invk(2,i)==0)then
                 tmpgk(kmap(1,i),kmap(2,i),kmap(3,i),j)=Gk(invk(1,i),j,olist(m,2),olist(l,2)) !G42(iw)
                 tmpgk(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=conjg(Gk(invk(1,i),j,olist(l,2),olist(m,2))) !G42(k,-iw)=G^*24(k,iw)
              else if(invk(2,i)==1)then
                 tmpgk(kmap(1,i),kmap(2,i),kmap(3,i),j)=Gk(invk(1,i),j,olist(l,2),olist(m,2)) !G42(-k,iw)=G24(k,iw)
                 tmpgk(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=conjg(Gk(invk(1,i),j,olist(m,2),olist(l,2))) !G42(-k,-iw)=G^*42(k,iw)
              end if
           end do
        end do
        !$omp end do
        !$omp end parallel
        call FFT(tmpVsigma,tmp,Nx,Ny,Nz,2*Nw,.true.)
        call FFT(tmpgk,tmp,Nx,Ny,Nz,2*Nw,.true.)
        !$omp parallel
        !$omp do private(i,j,k,n)
        do n=1,2*Nw
           do k=0,Nz-1
              do j=0,Ny-1
                 do i=0,Nx-1
                    tmp(i,j,k,n)=tmpVsigma(i,j,k,n)*tmpgk(i,j,k,n)
                 end do
              end do
           end do
        end do
        !$omp end do
        !$omp workshare
      tmpgk=cmplx(0.0d0,0.0d0,kind=c_double)
        !$omp end workshare
        !$omp end parallel
        call FFT(tmp,tmpgk,Nx,Ny,Nz,2*Nw,.false.)
        !$omp parallel do private(i,j)
        do j=1,Nw
           do i=1,Nkall
              if(invk(2,i)==0)then
                 sigmak(invk(1,i),j,olist(m,1),olist(l,1))=sigmak(invk(1,i),j,olist(m,1),olist(l,1))&
                      +tmp(kmap(1,i),kmap(2,i),kmap(3,i),j)
              end if
           end do
        end do
        !$omp end parallel do
     end do
  end do

  call destroy_fft_plans()
  do l=1,Norb
     do m=1,Norb
        !$omp parallel do private(i,j)
        do j=1,Nw
           do i=1,Nk
              sigmak(i,j,m,l)=sigmak(i,j,m,l)*weight
              if(abs(dble(sigmak(i,j,m,l)))<eps) sigmak(i,j,m,l)=cmplx(0.0d0,imag(sigmak(i,j,m,l)),kind=c_double)
              if(abs(imag(sigmak(i,j,m,l)))<eps) sigmak(i,j,m,l)=cmplx(dble(sigmak(i,j,m,l)),0.0d0,kind=c_double)
           end do
        end do
        !$omp end parallel do
     end do
  end do
end subroutine calc_sigma

subroutine mkself(sigmak,mu,Smat,Cmat,kmap,invk,olist,hamk,eig,uni,mu_init,rfill,temp,&
     scf_loop,pp,eps,Nkall,Nk,Nw,Norb,Nchi,Nx,Ny,Nz,sub_sigma,sw_out,sw_in,m_diis,sw_rescale) bind(C)
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nw,Norb,Nchi,Nkall,Nk,Nx,Ny,Nz,scf_loop,m_diis
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
  logical(1),intent(in):: sw_in,sw_out,sw_rescale
  integer(c_int64_t),intent(in):: sub_sigma
  real(c_double),intent(in):: temp,eps,pp,rfill,mu_init
  real(c_double),intent(in),dimension(Norb,Nk):: eig
  real(c_double),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  real(c_double),intent(out):: mu
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni,hamk
  complex(c_double),intent(out),dimension(Nk,Nw,Norb,Norb):: sigmak

  integer(c_int32_t) scf_i
  integer(c_int32_t),dimension(Nchi,Nchi,2)::chi_map
  integer(c_int32_t),dimension(Nchi*(Nchi+1)/2,2)::irr_chi
  real(c_double)esterr,mu_OLD,eps_sgm,maxchi0s_global
  complex(c_double),dimension(Nk,Nw,Norb,Norb):: Gk,sigmak0
  complex(c_double),dimension(Nk,Nw,Nchi,Nchi):: chi
  ! DIIS
  integer(c_int32_t):: n_hist,i_hist
  complex(c_double),allocatable:: xout_hist(:,:,:,:,:),res_hist(:,:,:,:,:)
  real(c_double),allocatable:: B_diis(:,:),rhs_diis(:)
  integer(c_int32_t),allocatable:: ipiv_diis(:)
  ! bracket cache for renew_mu
  real(c_double):: muS_cache,muL_cache
  logical:: bracket_valid
  ! eigenvalues of H(k)+Sigma(k,iw) precomputed before mu search
  complex(c_double),allocatable:: eig_gkinv(:,:,:)

  eps_sgm=1.0d-10
  mu=mu_init
  if(sw_in)then
     print*,"load self"
     call io_sigma(.false.)
     !$omp parallel
     !$omp workshare
     sigmak0(:,:,:,:)=sigmak(:,:,:,:)
     !$omp end workshare
     !$omp end parallel
     call gen_green_inv(Gk,sigmak,hamk,mu,temp,Nk,Nw,Norb)
     call getinv(Gk,Nk,Nw,Norb)
  else
     mu_old=mu*1.2
     !$omp parallel workshare
     sigmak0(:,:,:,:)=0.0d0
     Gk(:,:,:,:)=0.0d0 !gen_green0 need to initialization of Gk
     !$omp end parallel workshare
     call gen_green0(Gk,eig,uni,mu,temp,Nk,Nw,Norb)
  end if
  call get_chi_map(chi_map,irr_chi,olist,Nchi)
  allocate(xout_hist(Nk,Nw,Norb,Norb,m_diis))
  allocate(res_hist (Nk,Nw,Norb,Norb,m_diis))
  allocate(B_diis(m_diis+1,m_diis+1))
  allocate(rhs_diis(m_diis+1))
  allocate(ipiv_diis(m_diis+1))
  allocate(eig_gkinv(Norb,Nk,Nw))
  n_hist=0
  i_hist=0
  bracket_valid=.false.
  muS_cache=0.0d0
  muL_cache=0.0d0
  iter_loop: do scf_i=1,scf_loop
     print'(A5,I5)','iter=',scf_i
     call get_chi0_conv(chi,Gk,kmap,invk,irr_chi,chi_map,olist,temp,Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi)
     call ckchi_impl(chi,Smat,Cmat,kmap,invk,Nk,Nkall,Nchi,Nw,maxchi0s_global)
     if(sw_rescale .and. maxchi0s_global>=1.0d0)then
        print'(A,F10.6,A)','[FLEX] Stoner factor=',maxchi0s_global,'>= 1: rescaling chi0'
        chi(:,:,:,:)=chi(:,:,:,:)*(1.0d0-1.0d-4)/maxchi0s_global
     end if
     call get_Vsigma_flex_nosoc(chi,Smat,Cmat,Nk,Nw,Nchi)
     print'(A16,E12.4,A5,E12.4)','Re V_sigma: max:',maxval(dble(chi)),' min:',minval(dble(chi))
     call calc_sigma(sigmak,Gk,chi,Smat,Cmat,kmap,invk,olist,temp,Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz)
     if(sub_sigma>0)then
        sub_self:block
          integer(c_int32_t) iw,l,m
          complex(c_double),dimension(Nk,Norb,Norb):: sub_sigmak
          if(sub_sigma==1)then
             !$omp parallel
             do l=1,Norb
                do m=l,Norb
                   !$omp workshare
                   sub_sigmak(:,l,m)=(sigmak(:,1,l,m)+conjg(sigmak(:,1,m,l)))*0.5d0
                   !$omp end workshare
                   if(l.ne.m)sub_sigmak(:,m,l)=conjg(sub_sigmak(:,l,m))
                end do
             end do
             !$omp end parallel
          else  !sub_sigma==2: frequency-average HF subtraction
             !$omp parallel
             do l=1,Norb
                do m=l,Norb
                   !$omp workshare
                   sub_sigmak(:,l,m)=sum(sigmak(:,1:Nw,l,m)+conjg(sigmak(:,1:Nw,m,l)),dim=2)*(0.5d0/dble(Nw))
                   !$omp end workshare
                   if(l.ne.m)sub_sigmak(:,m,l)=conjg(sub_sigmak(:,l,m))
                end do
             end do
             !$omp end parallel
          end if
          !$omp parallel do private(iw)
          do iw=1,Nw
             sigmak(:,iw,:,:)=sigmak(:,iw,:,:)-sub_sigmak(:,:,:)
          end do
          !$omp end parallel do
        end block sub_self
     end if
     call compare_sigma()
     if(esterr<eps)then
        exit
     end if
     call renew_mu()
     call gen_green_inv(Gk,sigmak,hamk,mu,temp,Nk,Nw,Norb)
     call getinv(Gk,Nk,Nw,Norb)
     !$omp parallel
     !$omp workshare
     sigmak0(:,:,:,:)=sigmak(:,:,:,:)
     !$omp end workshare
     !$omp end parallel
  end do iter_loop
  deallocate(xout_hist,res_hist,B_diis,rhs_diis,ipiv_diis)
  if(sw_out)then
     call io_sigma(.true.)
  end if
  call renew_mu()
  deallocate(eig_gkinv)
contains
  subroutine compare_sigma()
    integer(c_int32_t) i,j,l,m,kerr,iwerr,lerr,merr,ih,jh,idx_i,idx_j,info,n_cur
   real(c_double) est, eps_reg
    complex(c_double) cksigm
    complex(c_double),dimension(Nk,Nw,Norb,Norb):: sigma_diis

    ! --- Store current output and residual in circular buffer ---
    i_hist=mod(i_hist,int(m_diis,c_int32_t))+1
    n_hist=min(n_hist+1,int(m_diis,c_int32_t))
    xout_hist(:,:,:,:,i_hist)=sigmak(:,:,:,:)
    res_hist(:,:,:,:,i_hist)=sigmak(:,:,:,:)-sigmak0(:,:,:,:)
    n_cur=n_hist

    ! --- Build Pulay matrix B of size (n_cur+1 x n_cur+1) ---
    B_diis(1:n_cur+1,1:n_cur+1)=0.0d0
    do ih=1,n_cur
       idx_i=modulo(i_hist-n_cur+ih-1,int(m_diis,c_int32_t))+1
       do jh=1,n_cur
          idx_j=modulo(i_hist-n_cur+jh-1,int(m_diis,c_int32_t))+1
          B_diis(ih,jh)=real(sum(conjg(res_hist(:,:,:,:,idx_i))*res_hist(:,:,:,:,idx_j)),kind=c_double)
       end do
       B_diis(ih,n_cur+1)=1.0d0
       B_diis(n_cur+1,ih)=1.0d0
    end do

    ! Regularize diagonal to improve numerical stability
    eps_reg=1.0d-12
    if(n_cur>0)then
       do ih=1,n_cur
          B_diis(ih,ih)=B_diis(ih,ih)+eps_reg
       end do
    end if

    ! --- Right-hand side vector (Lagrange constraint: sum c_i = 1) ---
    rhs_diis(1:n_cur)=0.0d0
    rhs_diis(n_cur+1)=1.0d0

    ! --- Solve B*c = rhs via LAPACK dgesv ---
    call dgesv(n_cur+1,1,B_diis(1:n_cur+1,1:n_cur+1),n_cur+1, &
               ipiv_diis(1:n_cur+1),rhs_diis(1:n_cur+1),n_cur+1,info)
    ! If dgesv fails (singular matrix), fall back to using most-recent vector
    if(info/=0)then
       print*,'DIIS: dgesv failed (info=',info,'), fallback to most-recent entry'
       rhs_diis(1:n_cur)=0.0d0
       if(n_cur>=1) then
          rhs_diis(n_cur)=1.0d0
       end if
    end if

    ! --- DIIS extrapolation: sigma_diis = sum_i c_i * xout_hist_i ---
    sigma_diis(:,:,:,:)=0.0d0
    do ih=1,n_cur
       idx_i=modulo(i_hist-n_cur+ih-1,int(m_diis,c_int32_t))+1
       sigma_diis=sigma_diis+rhs_diis(ih)*xout_hist(:,:,:,:,idx_i)
    end do

    ! --- Convergence check and mixing ---
    ! n_cur=1 (first step or m_diis=1): fall back to linear mixing
    ! n_cur>=2 (DIIS active): use sigma_diis directly to preserve optimal extrapolation
    esterr=0.0d0
    est=100.0d0
    do l=1,Norb
       do m=1,Norb
          do j=1,Nw
             do i=1,Nk
                cksigm=sigma_diis(i,j,m,l)
                if(abs(cksigm)>eps_sgm)then
                   est=abs((sigmak0(i,j,m,l)-cksigm)/cksigm)
                   if(est>esterr)then
                      esterr=est
                      kerr=i
                      iwerr=j
                      lerr=l
                      merr=m
                   end if
                end if
                if(n_cur>=2)then
                   sigmak(i,j,m,l)=cksigm
                else
                   sigmak(i,j,m,l)=pp*cksigm+(1.0d0-pp)*sigmak0(i,j,m,l)
                end if
             end do
          end do
       end do
    end do
    do i=1,Nkall
       if(invk(1,i)==kerr)then
          print '(A7,E9.2,A14,3(1X,I3),I5,2I3,I5)','esterr=',esterr,' at k,iw,m,l=', &
               & kmap(1,i),kmap(2,i),kmap(3,i),iwerr,merr,lerr,kerr
          exit
       end if
    end do
  end subroutine compare_sigma

  subroutine renew_mu()
    integer(c_int32_t) i_iter
    integer(c_int32_t),parameter:: itemax=100
    logical flag
    logical bracket_found
    real(c_double) rnS,rnL,rnc,rnM,muc,mud,muL,muS,muM,eps,dmu

    call precompute_eig_gkinv()
    if(esterr>1.0d-2)then
       eps= 1.0d-8
    else
       eps= esterr*1.0d-1
    end if
    dmu= abs(mu-mu_OLD)*2.0d0
    if (dmu<eps*4.0d0) dmu= eps*4.0d0
    mu_OLD=mu

    bracket_found=.false.
    if(bracket_valid)then
       call get_rn(rnS,muS_cache)
       call get_rn(rnL,muL_cache)
       if(rnS<rfill .and. rnL>rfill)then
          muS=muS_cache
          muL=muL_cache
          bracket_found=.true.
          print'(A)','[mu] cached bracket valid'
       else
          print'(A)','[mu] cached bracket invalid, re-searching'
       end if
    end if

    if(.not.bracket_found)then
       muL= mu+dmu
       muS= mu-dmu
       upper_lim: do i_iter=1,itemax
          mu=muL
          rnL=0.0d0
          call get_rn(rnL,mu)
          if(rnL<rfill)then
             if(abs(rfill-rnL)>0.5d0)then
                muL=muL+1.0d0
             else if(abs(rfill-rnL)>dmu)then
                muL=muL+0.5d0
             else
                muL= muL +dmu
             end if
          else
             exit
          end if
          if(i_iter==itemax)then
             print*,'Too many'
             stop
          end if
       end do upper_lim

       lower_lim: do i_iter=1,itemax
          mu=muS
          rnS=0.0d0
          call get_rn(rnS,mu)
          if(rnS>rfill)then
             if(abs(rnS-rfill)>0.5d0)then
                muS=muS-1.0d0
             else if(abs(rnS-rfill)>dmu)then
                muS=muS-0.5d0
             else
                muS=muS-dmu
             end if
          else
             exit
          end if
          if(i_iter==itemax)then
             print*,'Too many'
             stop
          end if
       end do lower_lim
    end if

    ! save wide bracket before Brent narrows it
    muS_cache=muS
    muL_cache=muL
    bracket_valid=.true.

    rnL=rnL-rfill
    rnS=rnS-rfill
    rnc=rnS
    muc=muS
    mud=0.0d0
    flag=.false.
    brent_loop: do i_iter=1,itemax
       if(rnc/=rnS .and. rnc/=rnL)then
          muM=(muL*rnS*rnc*(rnS-rnc)+muS*rnL*rnc*(rnc-rnL)+muc*rnL*rnS*(rnL-rnS))&
               /((rnL-rnS)*(rnc-rnL)*(rnc-rnS))
       else
          muM=muL-rnL*(muL-muS)/(rnL-rnS)
       end if
       if((0.25d0*(3.0d0*muS+muL) > muM .or. muM > muL) .or. &
            (flag .and. abs(muM-muL) >= abs(muL-muc)*0.5d0) .or. &
            (flag .eqv. .false. .and. abs(muM-muL) >= abs(muc-mud)*0.5d0) .or. &
            (flag .and. abs(muL-muc)<1.0d-8) .or. &
            (flag.eqv. .false. .and. abs(muc-mud)<1.0d-8))then
          muM= (muL+muS)*0.5d0
          flag=.true.
       else
          flag=.false.
       end if
       if(abs(muL-muM)<eps)exit
       mu=muM
       rnM=0
       call get_rn(rnM,muM)
       !print '(1x,a,2f22.16,l2)','muM,rnM=   ',muM,rnM,flag
       mud=muc
       muc=muL
       rnc=rnL
       if(rnS*(rnM-rfill)<0)then
          muL=muM
          rnL=rnM-rfill
       else
          muS=muM
          rnS=rnM-rfill
       end if
       if(abs(rnS)<abs(rnL))then
          muM=muL
          muL=muS
          muS=muM
          rnM=rnL
          rnL=rnS
          rnS=rnM
       end if
       if(i_iter==itemax)then
          print *,'Too many loop!'
          stop
       end if
    end do brent_loop
    if(rnL==rnS)then
       mu=(muS+muL)*0.5d0
    else
       mu= (muS*rnL-muL*rnS)/(rnL-rnS)
    end if
    call get_rn(rnM,mu)
    mu_OLD=mu
    print'(A4,F8.4,A10,F8.4)','mu  =',mu,' filling =',rnM
  end subroutine renew_mu
  
  subroutine get_rn(rn,rmu)
    use constant
    use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
    real(c_double),intent(in):: rmu
    real(c_double),intent(out):: rn

    integer(c_int32_t) i,j,n
    real(c_double) tmp,deltagk
    complex(c_double) iw_im

    ! n(μ) = (1/Nk){Σ_{k,n} f(ε_{kn}) + 2T·Σ_{k,j}[G(k,iω_j) - G_0(k,iω_j)]}
    ! non-interacting Fermi sum (unchanged)
    tmp=sum(0.5d0*(1.0d0-tanh(0.5d0*(eig(:,:)-rmu)/temp)))
    ! correction from self-energy: use precomputed eigenvalues of H(k)+Sigma(k,iw): no getinv needed
    deltagk=0.0d0
    !$omp parallel do reduction(+:deltagk) private(n,iw_im)
    do j=1,Nw
       iw_im=cmplx(0.0d0,dble(2*j-1)*pi*temp,kind=c_double)
       do i=1,Nk
          do n=1,Norb
             deltagk=deltagk+dble(1.0d0/(rmu+iw_im-eig_gkinv(n,i,j)) &
                                  -1.0d0/(rmu+iw_im-eig(n,i)))
          end do
       end do
    end do
    !$omp end parallel do
    rn=(tmp+2*temp*deltagk)/Nk
  end subroutine get_rn

  subroutine precompute_eig_gkinv()
    use constant
    use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
    integer(c_int32_t) i,j,info_eig
    complex(c_double),dimension(Norb,Norb):: A_tmp
    complex(c_double),dimension(Norb):: eig_tmp
    complex(c_double),dimension(Norb,Norb):: vl_dum,vr_dum
    real(c_double),dimension(2*Norb):: rwork_eig
    complex(c_double),dimension(2*Norb):: work_eig

    !$omp parallel do private(A_tmp,eig_tmp,vl_dum,vr_dum,rwork_eig,work_eig,info_eig)
    do j=1,Nw
       do i=1,Nk
          A_tmp(:,:)=hamk(:,:,i)+sigmak(i,j,:,:)
          call zgeev('N','N',Norb,A_tmp,Norb,eig_tmp,vl_dum,Norb,vr_dum,Norb, &
                     work_eig,2*Norb,rwork_eig,info_eig)
          if(info_eig/=0)then
             print*,'zgeev failed in precompute_eig_gkinv: info=',info_eig
             stop
          end if
          eig_gkinv(:,i,j)=eig_tmp(:)
       end do
    end do
    !$omp end parallel do
  end subroutine precompute_eig_gkinv
  
  subroutine io_sigma(sw)
    logical,intent(in):: sw !True: out, False: in
    open(55,file='sigma.bin',form='unformatted')
    if(sw)then
       write(55)mu
       write(55)mu_OLD
       write(55)sigmak
    else
       read(55)mu
       read(55)mu_OLD
       read(55)sigmak
    end if
    close(55)
  end subroutine io_sigma
end subroutine mkself

subroutine get_chis_chic(chis,chic,chi,Smat,Cmat,Nk,Nw,Nchi) bind(C)
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,Nchi
  real(c_double),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  complex(c_double),intent(in),dimension(Nk,Nw,Nchi,Nchi):: chi
  complex(c_double),intent(out),dimension(Nk,Nchi,Nchi):: chis,chic

  integer(c_int32_t) i,l,info
  integer(c_int32_t),dimension(Nchi):: ipiv
  complex(c_double),dimension(Nchi,Nchi):: cmat1,cmat2,cmat3,cmat4,Smat_c,Cmat_c

  Smat_c=cmplx(Smat,0.0d0,kind=c_double)
  Cmat_c=cmplx(Cmat,0.0d0,kind=c_double)

  qloop:do i=1,Nk
     call zgemm('N','N',Nchi,Nchi,Nchi,(-1.0d0,0.0d0),chi(i,1,:,:),Nchi,Smat_c,Nchi,(0.0d0,0.0d0),cmat1,Nchi) !-chi0S
     call zgemm('N','N',Nchi,Nchi,Nchi,(1.0d0,0.0d0),chi(i,1,:,:),Nchi,Cmat_c,Nchi,(0.0d0,0.0d0),cmat2,Nchi)  !chi0C
     do l=1,Nchi
        cmat1(l,l)=cmat1(l,l)+1.0d0 !I-chi0S
        cmat2(l,l)=cmat2(l,l)+1.0d0 !I+chi0C
     end do
     cmat3(:,:)=chi(i,1,:,:) !RHS copy for spin susceptibility
     cmat4(:,:)=chi(i,1,:,:) !RHS copy for charge susceptibility
     ! RPA spin susceptibility:   chi_s = (I - chi0·S)^{-1} chi0
     call zgesv(Nchi,Nchi,cmat1,Nchi,ipiv,cmat3,Nchi,info) !(I-chi0S)X=chi0 -> cmat3=chiS
     if(info/=0)then; print*,'zgesv failed: info=',info; stop; end if
     ! RPA charge susceptibility: chi_c = (I + chi0·C)^{-1} chi0
     call zgesv(Nchi,Nchi,cmat2,Nchi,ipiv,cmat4,Nchi,info) !(I+chi0C)X=chi0 -> cmat4=chiC
     if(info/=0)then; print*,'zgesv failed: info=',info; stop; end if
     chis(i,:,:)=cmat3(:,:)
     chic(i,:,:)=cmat4(:,:)
  end do qloop
end subroutine get_chis_chic

subroutine get_eig_or_tr_chi(chiqout,chi,invk,Nkall,Nk,Nchi,sw_eig) bind(C)
   use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
   implicit none
   integer(c_int64_t),intent(in):: Nkall,Nk,Nchi
   integer(c_int64_t),intent(in),dimension(3,Nkall):: invk
   logical,intent(in):: sw_eig
   complex(c_double),intent(in),dimension(Nk,Nchi,Nchi):: chi
   complex(c_double),intent(out),dimension(Nkall):: chiqout

   integer(c_int32_t) i,l,info
   real(c_double),dimension(2*Nchi):: rwork
   complex(c_double),dimension(Nchi*Nchi*4+1):: work
   complex(c_double),dimension(Nchi):: eig
   complex(c_double),dimension(Nchi,Nchi):: tmp,tmp1,tmp2

   !$omp parallel
   !$omp workshare
   chiqout(:)=0.0d0
   !$omp end workshare
   !$omp do private(tmp,tmp1,tmp2,rwork,work,eig,i,l,info)
   do i=1,Nkall
      if(sw_eig)then
         tmp(:,:)=chi(invk(1,i),:,:)
         call zgeev('N','N',Nchi,tmp,Nchi,eig,tmp1,Nchi,tmp2,Nchi,work,Nchi*Nchi*4+1,rwork,info)
         chiqout(i)=maxval(dble(eig))
      else
         do l=1,Nchi
            chiqout(i)=chiqout(i)+chi(invk(1,i),l,l)
         end do
      end if
   end do
   !$omp end do
   !$omp end parallel
end subroutine