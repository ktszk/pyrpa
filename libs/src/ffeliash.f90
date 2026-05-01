!non-linearied eliashberg equations solver
subroutine eliashberg(delta,sigmak,hamk,uni,eig,init_delta,Smat,Cmat,olist,prt,kmap,invk,mu,temp,eps,&
  Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,itemax,gap_sym,sw_sigma,m_diis) bind(C)
   use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t,c_bool
  implicit none
  integer(c_int64_t),intent(in):: Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,itemax,gap_sym,m_diis
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
   logical(c_bool),intent(in):: sw_sigma
  real(c_double),intent(in):: temp,eps,mu
  real(c_double),intent(in),dimension(Norb):: prt
  real(c_double),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  real(c_double),intent(in),dimension(Nk,Norb):: init_delta
  real(c_double),intent(in),dimension(Norb,Nk):: eig
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni,hamk
  complex(c_double),intent(out),dimension(Nk,Nw,Norb,Norb):: delta,sigmak

  integer(c_int32_t) i_iter,i_eig,i
  integer(c_int32_t),dimension(Nchi,Nchi,2)::chi_map
  integer(c_int32_t),dimension(Nchi*(Nchi+1)/2,2)::irr_chi
  logical(1) sw_pair
  real(c_double) norm,inorm,weight,lambda_rq,lambda_prev,lambda_phys,lambda1,vec_err,proj
  real(c_double)esterr,mu_OLD,eps_sgm,maxchi0s_global
  complex(c_double),dimension(Nk,Nw,Norb,Norb):: newdelta,fk,delta1,Gk,Gk0
  complex(c_double),allocatable,dimension(:,:,:,:):: Vdelta,Vsigma

  allocate(Vdelta(Nchi,Nchi,Nw,Nk),Vsigma(Nchi,Nchi,Nw,Nk))
  if(gap_sym>=0)then
     sw_pair=.true.
     print'(A7)','singlet'
  else
     sw_pair=.false.
     print'(A7)','triplet'
  end if
  weight=temp/dble(Nkall)
  call gen_green0(Gk,eig,uni,mu,temp,Nk,Nw,Norb)
  call get_initial_delta(delta,init_delta,uni,kmap,invk,Nkall,Nk,Nw,Norb,gap_sym)
  call mkfk_trs_nsoc(fk,Gk,delta,Nk,Nw,Norb)
  call get_chi0_conv(Vdelta,Gk,kmap,invk,irr_chi,chi_map,olist,temp,Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi)
  call ckchi_impl(Vdelta,Smat,Cmat,kmap,invk,Nk,Nkall,Nchi,Nw,maxchi0s_global)
  call get_V_delta_nsoc_flex(Vdelta,Smat,Cmat,Nk,Nw,Nchi,sw_pair)
  ! mkV_flex_nosoc expects both channels (chi0S/chi0C). Initialize spin channel from charge channel at startup.
  Vsigma(:,:,:,:)=Vdelta(:,:,:,:)
  if(sw_sigma)then
     call get_Vsigma_flex_nosoc(Vsigma,Smat,Cmat,Nk,Nw,Nchi)
  end if
  iter_loop: do i_iter=1,itemax
     print'(A5,I5)','iter=',i_iter
     call mkdelta_nsoc(newdelta,fk,Vdelta,Smat,Cmat,kmap,invk,prt,olist,Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,sw_pair)
     if(sw_sigma)then !Gk0 is (G^-1-sigma)^-1
        call calc_sigma(sigmak,Gk,Vsigma,Smat,Cmat,kmap,invk,olist,temp,Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz)
        call gen_green_inv(Gk0,sigmak,hamk,mu,temp,Nk,Nw,Norb)
        call getinv(Gk0,Nk,Nw,Norb)
     else
        call gen_green0(Gk0,eig,uni,mu,temp,Nk,Nw,Norb)
     end if
     !G=Gk0(I+DeltaFk)
     call mkgreliah_trs_nsoc(Gk,Gk0,fk,delta,Nk,Nw,Norb)
     !F=Gk0DeltaGk
     call mkfkeliash_trs_nsoc(fk,Gk0,Gk,delta,Nk,Nw,Norb) 
     call mkV_flex_nosoc(Vdelta,Vsigma,Smat,Cmat,Nk,Nw,Nchi,sw_pair)
  end do iter_loop
  deallocate(Vdelta,Vsigma)
end subroutine eliashberg

subroutine mkgreliah_trs_nsoc(Gk,Gk0,fk,delta,Nk,Nw,Norb)
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,Norb
   complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk0,fk,delta
   complex(c_double),intent(out),dimension(Nk,Nw,Norb,Norb):: Gk
 
  integer(c_int32_t) i,j,l,m
  complex(c_double),dimension(Norb,Norb):: cmat1

  !$omp parallel do collapse(2) private(i,j,l,m,cmat1)
  do j=1,Nw
     do i=1,Nk
        ! cmat1 = I+Delta^+(k,iw) * F^+(k,iw)
        call zgemm('N','C',Norb,Norb,Norb,(1.0d0,0.0d0),delta(i,j,:,:),Norb,fk(i,j,:,:),Norb,(0.0d0,0.0d0),cmat1,Norb)
        do l=1,Norb
            cmat1(l,l)=cmat1(l,l)+1.0d0
        end do
        ! Gk = Gk0* cmat1
        call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),Gk0(i,j,:,:),Norb,cmat1,Norb,(0.0d0,0.0d0),Gk(i,j,:,:),Norb)
     end do
  end do
  !$omp end parallel do
end subroutine mkgreliah_trs_nsoc

subroutine mkfkeliash_trs_nsoc(fk,Gk0,Gk,delta,Nk,Nw,Norb)
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,Norb
  complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk,Gk0,delta
  complex(c_double),intent(out),dimension(Nk,Nw,Norb,Norb):: fk
 
  integer(c_int32_t) i,j,l,m
  complex(c_double),dimension(Norb,Norb):: cmat1

  !$omp parallel do collapse(2) private(i,j,l,m,cmat1)
  do j=1,Nw
     do i=1,Nk
        ! cmat1 = Delta(k,iw) * G(k,iw)
        call zgemm('N','C',Norb,Norb,Norb,(1.0d0,0.0d0),delta(i,j,:,:),Norb,Gk(i,j,:,:),Norb,(0.0d0,0.0d0),cmat1,Norb)
                ! enforce Hermitian form used in the original implementation
                ! fk = Gk0* cmat1
        call zgemm('N','N',Norb,Norb,Norb,(-1.0d0,0.0d0),Gk0(i,j,:,:),Norb,cmat1,Norb,(0.0d0,0.0d0),fk(i,j,:,:),Norb)
        do l=1,Norb
           fk(i,j,l,l)=dble(fk(i,j,l,l))
           do m=l+1,Norb
              fk(i,j,l,m)=conjg(fk(i,j,m,l))
           end do
        end do
     end do
  end do
  !$omp end parallel do
end subroutine mkfkeliash_trs_nsoc

subroutine mkV_flex_nosoc(Vdelta,Vsigma,Smat,Cmat,Nk,Nw,Nchi,sw_pair)
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,Nchi
  logical(1),intent(in):: sw_pair
  real(c_double),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  complex(c_double),intent(inout),dimension(Nchi,Nchi,Nw,Nk):: Vdelta,Vsigma

  integer(c_int32_t) i,j,l,info
  integer(c_int32_t),dimension(Nchi):: ipiv
  complex(c_double),dimension(Nchi,Nchi):: cmat1,cmat2,cmat3,cmat4,cmat5,Smat_c,Cmat_c,V0_c,SC_c

  Smat_c(:,:)=cmplx(Smat(:,:),0.0d0,kind=c_double)
  Cmat_c(:,:)=cmplx(Cmat(:,:),0.0d0,kind=c_double)
  SC_c(:,:)=Smat_c(:,:)+Cmat_c(:,:)
  if(sw_pair)then
    V0_c(:,:)=cmplx(0.5d0*(Smat(:,:)+Cmat(:,:)),0.0d0,kind=c_double) !bare Vud=(C+S)/2
  else
    V0_c(:,:)=cmplx(0.5d0*(Smat(:,:)-Cmat(:,:)),0.0d0,kind=c_double) !bare Vuu=(S-C)/2
  end if
  !$omp parallel do collapse(2) private(i,cmat1,cmat2,cmat3,cmat4,cmat5,ipiv,info,l)
  wloop:do j=1,Nw
     qloop:do i=1,Nk
        call zgemm('N','N',Nchi,Nchi,Nchi,(-1.0d0,0.0d0),Vsigma(:,:,j,i),Nchi,Smat_c,Nchi,(0.0d0,0.0d0),cmat1,Nchi) !-chi0S
        call zgemm('N','N',Nchi,Nchi,Nchi,(1.0d0,0.0d0),Vdelta(:,:,j,i),Nchi,Cmat_c,Nchi,(0.0d0,0.0d0),cmat2,Nchi)  !chi0C
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
        cmat1(:,:)=cmat3(:,:) !cmat1=chiS
        cmat2(:,:)=cmat4(:,:) !cmat2=chiC
        cmat3(:,:)=1.5d0*Smat(:,:)-0.5d0*Cmat(:,:)   !static (HF-like) bare vertex
        cmat4(:,:)=V0_c(:,:)   !start from bare static vertex
        call zgemm('N','N',Nchi,Nchi,Nchi,(1.5d0,0.0d0),Smat_c,Nchi,cmat1,Nchi,(1.0d0,0.0d0),cmat3,Nchi)  !+3/2·S·χ_s
        call zgemm('N','N',Nchi,Nchi,Nchi,(0.5d0,0.0d0),Cmat_c,Nchi,cmat2,Nchi,(1.0d0,0.0d0),cmat3,Nchi)  !+1/2·C·χ_c
        call zgemm('N','N',Nchi,Nchi,Nchi,(-0.25d0,0.0d0),SC_c,Nchi,cmat5,Nchi,(1.0d0,0.0d0),cmat3,Nchi) !-1/4·(S+C)·(χ_s+χ_c), subtract double count
        if(sw_pair)then !singlet
           call zgemm('N','N',Nchi,Nchi,Nchi,(1.5d0,0.0d0),Smat_c,Nchi,cmat1,Nchi,(1.0d0,0.0d0),cmat4,Nchi)  !+3/2·S·χ_s
           call zgemm('N','N',Nchi,Nchi,Nchi,(-0.5d0,0.0d0),Cmat_c,Nchi,cmat2,Nchi,(1.0d0,0.0d0),cmat4,Nchi) !-1/2·C·χ_c
        else !triplet
           call zgemm('N','N',Nchi,Nchi,Nchi,(-0.5d0,0.0d0),Smat_c,Nchi,cmat1,Nchi,(1.0d0,0.0d0),cmat4,Nchi) !-1/2·S·χ_s
           call zgemm('N','N',Nchi,Nchi,Nchi,(-0.5d0,0.0d0),Cmat_c,Nchi,cmat2,Nchi,(1.0d0,0.0d0),cmat4,Nchi) !-1/2·C·χ_c
        end if
        Vsigma(i,j,:,:)=cmat3(:,:)
        Vdelta(i,j,:,:)=cmat4(:,:)
     end do qloop
  end do wloop
  !$omp end parallel do
end subroutine mkV_flex_nosoc

subroutine get_chi0sc(Vsigma,Vdelta,Gk,fk,kmap,invk,irr_chi,chi_map,olist,temp,Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi,sw_pair)
   use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
   implicit none
   integer(c_int64_t),intent(in):: Nx,Ny,Nz,Nk,Nw,Norb,Nkall,Nchi
   integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
   integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
   integer(c_int32_t),dimension(Nchi,Nchi,2)::chi_map
   integer(c_int32_t),dimension(Nchi*(Nchi+1)/2,2)::irr_chi
   logical(1),intent(in):: sw_pair
   real(c_double),intent(in):: temp
   complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk,fk
   complex(c_double),intent(out),dimension(Nchi,Nchi,Nw,Nk):: Vsigma,Vdelta
  
   integer(c_int32_t) i,j,l,m
   complex(c_double),dimension(Nchi,Nchi):: cmat1,cmat2
   

   call get_chi0_conv(Vsigma,Gk,kmap,invk,irr_chi,chi_map,olist,temp,Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi)
   call get_chi0_conv_ff(Vdelta,fk,kmap,invk,irr_chi,chi_map,olist,temp,Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi,sw_pair)
   !$omp parallel do collapse(2) private(i,j,l,m,cmat1)
   do j=1,Nw
       do i=1,Nk
         cmat1(:,:)=Vsigma(:,:,j,i)+Vdelta(:,:,j,i) !chi0S (RHS for first solve)
         cmat2(:,:)=Vsigma(:,:,j,i)-Vdelta(:,:,j,i) !chi0C (RHS for second solve)
         Vsigma(:,:,j,i)=cmat1(:,:) !cmat1=chi0S
         Vdelta(:,:,j,i)=cmat2(:,:) !cmat2=chi0
       end do
   end do
   !$omp end parallel do
end subroutine get_chi0sc

subroutine get_chi0_conv_ff(chi,Fk,kmap,invk,irr_chi,chi_map,olist,temp,Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi,sw_pair)
  !> This function obtains chi_0 using convolution.
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nw,Norb,Nchi,Nkall,Nk,Nx,Ny,Nz
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
  integer(c_int32_t),intent(in),dimension(Nchi,Nchi,2):: chi_map
  integer(c_int32_t),intent(in),dimension(Nchi*(Nchi+1)/2,2):: irr_chi
  logical(1),intent(in):: sw_pair
  real(c_double),intent(in):: temp
  complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: Fk
  complex(c_double),intent(out),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(c_int32_t) i,j,k,l,m,n,iorb
  integer(c_int32_t) ii(0:Nx-1),ij(0:Ny-1),ik(0:Nz-1),iw(2*Nw)
  real(c_double) weight,sgn
  real(c_double),parameter:: eps=1.0d-9
  complex(c_double),dimension(0:Nx-1,0:Ny-1,0:Nz-1,2*Nw):: tmp,tmpfk14,tmpfk23
  
   if(sw_pair)then
     sgn=-1.0d0 !triplet_dz
   else
     sgn=+1.0d0 !singlet
   end if
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
  call init_fft_plans(tmpfk14,tmp,Nx,Ny,Nz,2*Nw)
  orb_loop:do iorb=1,Nchi*(Nchi+1)/2
     l=irr_chi(iorb,1)
     m=irr_chi(iorb,2)
     !use symmetry G^lm(k,iw)=G^ml(k,-iw) from Hermitian symmetry of Hamiltonian
     !$omp parallel do private(i)
     w_loop_Gk_to_tmp:do j=1,Nw
        k_loop_Gk_to_tmp:do i=1,Nkall
           if(invk(2,i)==0)then
              !iw
              tmpfk14(kmap(1,i),kmap(2,i),kmap(3,i),j)=conjg(Fk(invk(1,i),j,olist(m,2),olist(l,1))) !F^+14(k,iw)
              tmpfk23(kmap(1,i),kmap(2,i),kmap(3,i),j)=Fk(invk(1,i),j,olist(l,2),olist(m,1)) !F23(k,iw)
              !-iw
              tmpfk14(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=sgn*conjg(Fk(invk(1,i),j,olist(l,1),olist(m,2))) !F^+14(k,-iw)=sgnF^*14(k,iw)
              tmpfk23(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=sgn*Fk(invk(1,i),j,olist(m,1),olist(l,2)) !F23(k,-iw)=sgnF32(k,iw)
           else if(invk(2,i)==1)then
              !iw
              tmpfk14(kmap(1,i),kmap(2,i),kmap(3,i),j)=conjg(Fk(invk(1,i),j,olist(m,1),olist(l,1))) !F^+14(-k,iw)=F31(k,iw)
              tmpfk23(kmap(1,i),kmap(2,i),kmap(3,i),j)=Fk(invk(1,i),j,olist(l,2),olist(m,2)) !F23(-k,iw)=F^24(k,iw)
              !-iw
              tmpfk14(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=sgn*conjg(Fk(invk(1,i),j,olist(l,1),olist(m,2))) !F^14(-k,-iw)=F^14(k,iw)
              tmpfk23(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=sgn*Fk(invk(1,i),j,olist(l,2),olist(m,1)) !G23(-k,-iw)=F32(k,iw)
           end if
        end do k_loop_Gk_to_tmp
     end do w_loop_Gk_to_tmp
     !$omp end parallel do
     call FFT(tmpfk14,tmp,Nx,Ny,Nz,2*Nw,.true.)
     call FFT(tmpfk23,tmp,Nx,Ny,Nz,2*Nw,.true.)
     !calculate G(r)G(-r)
     !$omp parallel do private(i,j,k)
     w_loop_conv:do n=1,2*Nw
        z_loop:do k=0,Nz-1
           y_loop:do j=0,Ny-1
              !$omp simd
              x_loop:do i=0,Nx-1
                 tmp(i,j,k,n)=-tmpfk14(i,j,k,n)*tmpfk23(ii(i),ij(j),ik(k),iw(n))
              end do x_loop
              !$omp end simd
           end do y_loop
        end do z_loop
     end do w_loop_conv
     !$omp end parallel do
     call FFT(tmp,tmpfk14,Nx,Ny,Nz,2*Nw,.false.)
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
end subroutine get_chi0_conv_ff