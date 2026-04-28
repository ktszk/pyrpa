module constant
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  real(c_double),parameter:: pi=3.141592653589793238462643383279d0
end module constant

subroutine openmp_params(omp_num,omp_check) bind(C)
  !$ use omp_lib
  use,intrinsic:: iso_c_binding, only:c_int64_t
  implicit none
  integer(c_int64_t),intent(out):: omp_num
  logical(1),intent(out):: omp_check
  !$ if(.true.)then
  !$   omp_check=.true.
  !$   omp_num=omp_get_max_threads()
  !$ else
       omp_check=.false.
       omp_num=0
  !$ end if
end subroutine openmp_params

subroutine gen_ham(ham_k,klist,ham_r,rvec,Nk,Nr,Norb) bind(C)
  !> This function generate model hamiltonian from hoppings
  !!@param ham_k,out: k-space Hamiltonian
  !!@param  klist,in: k-points list
  !!@param  ham_r,in: hopping integrals
  !!@param   rvec,in: The list of lattice site
  !!@param     Nk,in: The number of k-points
  !!@param     Nr,in: The number of r-vector
  !!@param   Norb,in: The number of orbitals
  use constant
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nr,Norb
  real(c_double),intent(in),dimension(3,Nk):: klist
  real(c_double),intent(in),dimension(3,Nr):: rvec
  complex(c_double),intent(in),dimension(Norb,Norb,Nr):: ham_r
  complex(c_double),intent(out),dimension(Norb,Norb,Nk):: ham_k

  integer(c_int32_t) i,j,l,m
  real(c_double) phase
  complex(c_double) cphase

   ham_k(:,:,:)=(0.0d0,0.0d0)

  ! H(k) = sum_R  H(R) * exp(-2pi*i * k.R)
  !$omp parallel do private(l,m,j,phase,cphase)
  klop: do i=1,Nk
     rloop: do j=1,Nr
        phase=2*pi*sum(klist(:,i)*rvec(:,j))
        cphase=cmplx(cos(phase),-sin(phase),kind=c_double)   ! exp(-i*phi)
        do l=1,Norb
           do m=l,Norb
              ham_k(m,l,i)=ham_k(m,l,i)+ham_r(m,l,j)*cphase
           end do
        end do
     end do rloop
     ! Enforce Hermitian symmetry: diagonal must be real; upper = conj(lower)
     do l=1,Norb
        ham_k(l,l,i)=dble(ham_k(l,l,i))
        do m=l+1,Norb
           ham_k(l,m,i)=conjg(ham_k(m,l,i))
        end do
     end do
  end do klop
  !$omp end parallel do
end subroutine gen_ham

subroutine get_eig(eig,uni,ham_k,Nk,Norb) bind(C)
  !> This function obtains energies
  !!@param  eig,out: energies at k-points
  !!@param  uni,out: unitary matrix at k-points
  !!@param ham_k,in: k-space Hamiltonian
  !!@param    Nk,in: The number of k-points
  !!@param  Norb,in: The number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,Norb
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: ham_k
  real(c_double),intent(out),dimension(Norb,Nk):: eig
  complex(c_double),intent(out),dimension(Norb,Norb,Nk):: uni

  integer(c_int32_t) i,info
  integer(c_int64_t) lwork
  real(c_double) rwork(3*Norb-2),eq(Norb)
  complex(c_double) work_query(1),en(Norb,Norb)
  complex(c_double),allocatable :: work(:)

  call zheev('V','U',Norb,en,Norb,eq,work_query,-1_c_int64_t,rwork,info)
  lwork = int(dble(work_query(1)), c_int64_t)

  !$omp parallel private(en,eq,work,rwork,info)
  allocate(work(lwork))
  !$omp do
  kloop: do i=1,Nk
     en(:,:)=ham_k(:,:,i)
     call zheev('V','U',Norb,en,Norb,eq,work,lwork,rwork,info)
     if(info/=0)then
        print*,'zheev failed in get_eig: info=',info
        stop
     end if
     uni(:,:,i)=en(:,:)
     eig(:,i)=eq(:)
  end do kloop
  !$omp end do
  deallocate(work)
  !$omp end parallel
end subroutine get_eig

subroutine get_eig_mlo(eig,uni,ham_k,Ovlk,Nk,Norb) bind(C)
  !> This function obtains energies
  !!@param  eig,out: energies at k-points
  !!@param  uni,out: unitary matrix at k-points
  !!@param ham_k,in: k-space Hamiltonian
  !!@param  Ovlk,in: k-space Overlap integrals
  !!@param    Nk,in: The number of k-points
  !!@param  Norb,in: The number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,norb
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: ham_k,Ovlk
  real(c_double),intent(out),dimension(Norb,Nk):: eig
  complex(c_double),intent(out),dimension(Norb,Norb,Nk):: uni

  integer(c_int32_t) i,j,k,l,m,info
  integer(c_int64_t) lwork
  real(c_double) rwork(3*Norb-2),eq(Norb),norm
  complex(c_double) work_query(1)
  complex(c_double),allocatable :: work(:)
  complex(c_double),dimension(Norb,Norb):: tmp,tmp2,tmp3
  real(c_double),parameter:: ovl_thresh=1.0d-8

  call zheev('V','U',norb,tmp,norb,eq,work_query,-1_c_int64_t,rwork,info)
  lwork = int(dble(work_query(1)), c_int64_t)

  !$omp parallel private(tmp,tmp2,tmp3,norm,eq,work,rwork,info)
  allocate(work(lwork))
  !$omp do
  kloop: do i=1,Nk
     tmp(:,:)=Ovlk(:,:,i)
     call zheev('V','U',norb,tmp,norb,eq,work,lwork,rwork,info)
     if(info/=0)then
        print*,'zheev failed in get_eig_mlo (overlap): info=',info
        stop
     end if
     ! Canonical orthogonalization: rescale eigenvectors by 1/sqrt(eigenvalue) of S
     do j=1,Norb
        if(eq(j)>ovl_thresh)then
           tmp2(:,j)=tmp(:,j)/sqrt(cmplx(eq(j),kind=c_double))
        else
           tmp2(:,j)=(0.0d0,0.0d0) !discard near-null basis vector (near-linear dependence)
        end if
     end do
     ! Transform H to orthonormal basis: H_orth = tmp2^H * H * tmp2
     call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),ham_k(:,:,i),Norb,tmp2,Norb,(0.0d0,0.0d0),tmp3,Norb)
     call zgemm('C','N',Norb,Norb,Norb,(1.0d0,0.0d0),tmp2,Norb,tmp3,Norb,(0.0d0,0.0d0),tmp,Norb)
     call zheev('V','U',norb,tmp,norb,eq,work,lwork,rwork,info)
     if(info/=0)then
        print*,'zheev failed in get_eig_mlo (ham): info=',info
        stop
     end if
     eig(:,i)=eq(:)
     ! uni = tmp2 * tmp  (compose transformations)
     call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),tmp2,Norb,tmp,Norb,(0.0d0,0.0d0),tmp3,Norb)
     uni(:,:,i)=tmp3(:,:)
  end do kloop
  !$omp end do
  deallocate(work)
  !$omp end parallel
end subroutine get_eig_mlo

subroutine get_ffermi(ffermi,eig,mu,temp,Nk,Norb) bind(C)
  !> This function obtains fermi functions
  !!@param ffermi,out: fermi distribute functions
  !!@param     eig,in: energies at k-points
  !!@param      mu,in: chemical potential
  !!@param    temp,in: Temperature
  !!@param      Nk,in: The number of k-points
  !!@param    Norb,in: The number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,Norb
  real(c_double),intent(in):: mu,temp
  real(c_double),intent(in),dimension(Norb,Nk):: eig
  real(c_double),intent(out),dimension(Norb,Nk):: ffermi

  integer(c_int32_t) i,j
   real(c_double) itemp,temp_safe

   temp_safe=max(temp,1.0d-12)
   itemp=0.5d0/temp_safe
  ! f(e) = 1/(exp((e-mu)/T)+1) = 0.5 - 0.5*tanh((e-mu)/(2T))  [numerically stable form]
  !$omp parallel do private(j)
  do i=1,Nk
     do j=1,Norb
        ffermi(j,i)=0.5d0-0.5d0*tanh((eig(j,i)-mu)*itemp)
     end do
  end do
  !$omp end parallel do
end subroutine get_ffermi

subroutine get_imass0(imk,klist,ham_r,rvec,Nk,Nr,Norb) bind(C)
  use constant
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nr,norb
  real(c_double),intent(in),dimension(3,Nk):: klist
  real(c_double),intent(in),dimension(3,Nr):: rvec
  complex(c_double),intent(in),dimension(Norb,Norb,Nr):: ham_r
  complex(c_double),intent(out),dimension(3,3,Norb,Norb,Nk):: imk

  integer(c_int32_t) i,j,k,l,m,n
  real(c_double) phase
  complex(c_double) cphase

   imk(:,:,:,:,:)=0.0d0

  !$omp parallel do private(l,m,j,k,n,phase,cphase)
  kloop: do i=1,Nk
     rloop: do j=1,Nr
        phase=2*pi*sum(klist(:,i)*rvec(:,j))
        cphase=cmplx(cos(phase),-sin(phase),kind=c_double)
        do l=1,Norb
           do m=l,Norb
              !$omp simd
              axis1: do k=1,3
                 axis2: do n=k,3
                    imk(n,k,m,l,i)=imk(n,k,m,l,i)-rvec(n,j)*rvec(k,j)*ham_r(m,l,j)*cphase
                 end do axis2
              end do axis1
              !$omp end simd
           end do
        end do
     end do rloop
     do l=1,Norb
        do m=l,Norb
           !$omp simd
           axis12: do k=1,3
              axis22: do n=k,3
                 imk(k,n,m,l,i)=imk(n,k,m,l,i)
                 imk(n,k,l,m,i)=conjg(imk(n,k,m,l,i))
                 imk(k,n,m,l,i)=imk(n,k,l,m,i)
              end do axis22
           end do axis12
           !$omp end simd
        end do
     end do
  end do kloop
  !$omp end parallel do
end subroutine get_imass0

subroutine get_imassk(imk,imk0,mrot,uni,Nk,Norb) bind(C)
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,Norb
  real(c_double),intent(in),dimension(3,3):: mrot
  complex(c_double),intent(in),dimension(3,3,Norb,Norb,Nk):: imk0
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni
  real(c_double),intent(out),dimension(3,3,Norb,Nk):: imk

  integer(c_int32_t) i,j_ax,k_ax,l,m,n
  complex(c_double) tmp(3,3,Norb)
  complex(c_double) Mtmp(Norb,Norb),Wtmp(Norb,Norb),res_mat(Norb,Norb)

  !$omp parallel do private(tmp,Mtmp,Wtmp,res_mat,j_ax,k_ax,l,m,n)
  kloop: do i=1,Nk
     !rotate orb to band via zgemm: tmp(k_ax,j_ax,n)=diag(U^H * imk0(k_ax,j_ax,:,:,i) * U)_n
     do k_ax=1,3
        do j_ax=1,3
           Mtmp(:,:)=imk0(k_ax,j_ax,:,:,i)
           call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),Mtmp,Norb,uni(:,:,i),Norb,(0.0d0,0.0d0),Wtmp,Norb)
           call zgemm('C','N',Norb,Norb,Norb,(1.0d0,0.0d0),uni(:,:,i),Norb,Wtmp,Norb,(0.0d0,0.0d0),res_mat,Norb)
           do n=1,norb
              tmp(k_ax,j_ax,n)=res_mat(n,n)
           end do
        end do
     end do
     !rotate axis
     imk(:,:,:,i)=0.0d0
     band_loop2: do n=1,norb
        do l=1,3
           do m=1,3
              do j_ax=1,3
                 do k_ax=1,3
                    imk(m,l,n,i)=imk(m,l,n,i)+mrot(m,k_ax)*mrot(l,j_ax)*dble(tmp(k_ax,j_ax,n))
                 end do
              end do
           end do
        end do
     end do band_loop2
  end do kloop
  !$omp end parallel do
end subroutine get_imassk

subroutine get_vlm0(vk,klist,ham_r,rvec,Nk,Nr,Norb) bind(C)
  !> This function obtain orbital basis velocity
  !!@param   vk,out: orbital-basis group velocities
  !!@param klist,in: k-points list
  !!@param ham_r,in: hopping integrals
  !!@param  rvec,in: The list of lattice site
  !!@param    Nk,in: The number of k-points
  !!@param    Nr,in: The number of r-vector
  !!@param  Norb,in: The number of orbitals
  use constant
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nr,Norb
  real(c_double),intent(in),dimension(3,Nk):: klist
  real(c_double),intent(in),dimension(3,Nr):: rvec
  complex(c_double),intent(in),dimension(Norb,Norb,Nr):: ham_r
  complex(c_double),intent(out),dimension(3,Norb,Norb,Nk):: vk

  integer(c_int32_t) i,j,k,l,m
  real(c_double) phase
  complex(c_double) cphase

   vk(:,:,:,:)=0.0d0

  !$omp parallel do private(l,m,j,k,phase,cphase)
  kloop: do i=1,Nk
     rloop: do j=1,Nr
        phase=2*pi*sum(klist(:,i)*rvec(:,j))
        ! d/dk_alpha exp(-i*phi) = -i*R_alpha * exp(-i*phi)
        ! Re/Im part: -(sin+i*cos) = -i*(cos-i*sin) = -i*exp(-i*phi)
        cphase=cmplx(sin(phase),cos(phase),kind=c_double)   ! = -i*exp(-i*phi) in real/imag components
        do l=1,Norb
           do m=l,Norb
              !$omp simd
              vaxis: do k=1,3
                 vk(k,m,l,i)=vk(k,m,l,i)-rvec(k,j)*ham_r(m,l,j)*cphase
              end do vaxis
              !$omp end simd
           end do
        end do
     end do rloop
     do l=1,Norb
        do m=l,Norb
           vaxis2: do k=1,3
              vk(k,l,m,i)=conjg(vk(k,m,l,i))
           end do vaxis2
        end do
     end do
  end do kloop
  !$omp end parallel do
end subroutine get_vlm0

subroutine get_veloc(vk,vk0,mrot,uni,Nk,Norb) bind(C)
  !> This function obtain group velocity at bands
  !!@param  vk,out: band-basis group velocities
  !!@param  vk0,in: orbital-basis group velocities
  !!@param mrot,in: rotation matrix for real space
  !!@param  uni,in: unitary matrix orbital to band
  !!@param   Nk,in: The number of k-points
  !!@param Norb,in: The number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,Norb
  real(c_double),intent(in),dimension(3,3):: mrot
  complex(c_double),intent(in),dimension(3,Norb,Norb,Nk):: vk0
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni
  real(c_double),intent(out),dimension(3,Norb,Nk):: vk

  integer(c_int32_t) i,j_dir,l,m,n
  complex(c_double) tmp(3,Norb)
  complex(c_double) vtmp(Norb,Norb),Wtmp(Norb,Norb),tmp_j(Norb,Norb)

  !$omp parallel do private(tmp,vtmp,Wtmp,tmp_j,j_dir,l,m,n)
  kloop: do i=1,Nk
     !rotate orb to band via zgemm: tmp_j = U^H * vk0(j,:,:,i) * U, take diagonal
     do j_dir=1,3
        vtmp(:,:)=vk0(j_dir,:,:,i)
        call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),vtmp,Norb,uni(:,:,i),Norb,(0.0d0,0.0d0),Wtmp,Norb)
        call zgemm('C','N',Norb,Norb,Norb,(1.0d0,0.0d0),uni(:,:,i),Norb,Wtmp,Norb,(0.0d0,0.0d0),tmp_j,Norb)
        do n=1,Norb
           tmp(j_dir,n)=tmp_j(n,n)
        end do
     end do
     !rotate axis
     band_loop2: do n=1,Norb
        do l=1,3
           do m=1,3
              vk(m,n,i)=vk(m,n,i)+mrot(m,l)*dble(tmp(l,n))
           end do
        end do
     end do band_loop2
  end do kloop
  !$omp end parallel do
end subroutine get_veloc

subroutine get_vnm(vk,vk0,mrot,uni,Nk,Norb) bind(C)
  !> This function obtain band-basis group velocity with inter-band
  !!@param  vk,out: band-basis group velocities
  !!@param  vk0,in: orbital-basis group velocities
  !!@param mrot,in: rotation matrix for real space
  !!@param  uni,in: unitary matrix orbital to band
  !!@param   Nk,in: The number of k-points
  !!@param Norb,in: The number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,Norb
  real(c_double),intent(in),dimension(3,3):: mrot
  complex(c_double),intent(in),dimension(3,Norb,Norb,Nk):: vk0
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(c_double),intent(out),dimension(3,Norb,Norb,Nk):: vk

  integer(c_int32_t) i,j_dir,l,m,n,k
  complex(c_double) tmp(3,Norb,Norb)
  complex(c_double),dimension(Norb,Norb):: vtmp,Wtmp,tmp_j

  !$omp parallel
  !$omp workshare
  vk(:,:,:,:)=0.0d0
  !$omp end workshare
  !$omp do private(tmp,vtmp,Wtmp,tmp_j,j_dir,l,m,n,k)
  kloop: do i=1,Nk
     ! tmp(j,:,:) = uni^H * vk0(j,:,:,i) * uni  for j=1..3
     do j_dir=1,3
        vtmp(:,:)=vk0(j_dir,:,:,i)
        call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),vtmp,Norb,uni(:,:,i),Norb,(0.0d0,0.0d0),Wtmp,Norb)
        call zgemm('C','N',Norb,Norb,Norb,(1.0d0,0.0d0),uni(:,:,i),Norb,Wtmp,Norb,(0.0d0,0.0d0),tmp_j,Norb)
        tmp(j_dir,:,:)=tmp_j(:,:)
     end do
     !rotate axis
     band_loop3: do n=1,Norb
        do k=1,Norb
           do l=1,3
              do m=1,3
                 vk(m,k,n,i)=vk(m,k,n,i)+mrot(m,l)*dble(tmp(l,k,n))
              end do
           end do
        end do
     end do band_loop3
  end do kloop
  !$omp end do
  !$omp end parallel
end subroutine get_vnm

subroutine gen_tr_greenw_0(trGk,wl,eig,mu,delta,Nk,Nw,Norb) bind(C)
  !> This function obtain orbital(band) trace of green function
  !!@param trGk,out: trace of green function
  !!@param    wl,in: list of energies
  !!@param   eig,in: energies at k-points
  !!@param    mu,in: chemical potential
  !!@param delta,in: dumping factor
  !!@param    Nk,in: The number of k-points
  !!@param    Nw,in: The number of energies mesh
  !!@param  Norb,in: The number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,Norb
  real(c_double),intent(in):: mu,delta
  real(c_double),intent(in),dimension(Norb,Nk):: eig
  real(c_double),intent(in),dimension(Nw):: wl
  complex(c_double),intent(out),dimension(Nw,Nk):: trGk

  integer(c_int32_t) i,j,n

  trGk(:,:) = (0.0d0, 0.0d0)

  !$omp parallel do private(i,n)
  kloop: do j=1,Nk
     !$omp simd
     wloop: do i=1,Nw
        band_loop: do n=1,Norb
           trGk(i,j)=trGk(i,j)+1./cmplx(wl(i)-eig(n,j)+mu,delta,kind=c_double)
        end do band_loop
     end do wloop
     !$omp end simd
  end do kloop
  !$omp end parallel do
end subroutine gen_tr_greenw_0

subroutine gen_dos(Dos,wl,eig,uni,mu,delta,Nk,Nw,Norb) bind(C)
  !> This function obtain partial Dos (and sum of them is total)
  !!@param  Dos,out: partial Density of states
  !!@param    wl,in: list of energies
  !!@param   eig,in: energies at k-points
  !!@param   uni,in: unitary matrix at k-points  
  !!@param    mu,in: chemical potential
  !!@param delta,in: dumping factor
  !!@param    Nk,in: The number of k-points
  !!@param    Nw,in: The number of energies mesh
  !!@param  Norb,in: The number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,Norb
  real(c_double),intent(in):: mu,delta
  real(c_double),intent(in),dimension(Norb,Nk):: eig
  real(c_double),intent(in),dimension(Nw):: wl
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(c_double),intent(out),dimension(Nw,Norb):: Dos

  integer(c_int32_t) i,j,k,n

  !$omp parallel
  !$omp workshare
  Dos(:,:)=0.0d0
  !$omp end workshare
  orb_loop: do j=1,Norb
     !$omp do private(i,k,n)
     wloop: do i=1,Nw
        kloop: do k=1,Nk
           bandloop: do n=1,Norb
              Dos(i,j)=Dos(i,j)+uni(j,n,k)*conjg(uni(j,n,k))/cmplx(wl(i)-eig(n,k)+mu,delta,kind=c_double)
           end do bandloop
        end do kloop
     end do wloop
     !$omp end do
  end do orb_loop
  !$omp workshare
  Dos(:,:)=Dos(:,:)/Nk
  !$omp end workshare
  !$omp end parallel
end subroutine gen_dos

subroutine get_parity_prop(Pmn,rvec,ham_r,Norb,Nr) bind(C)
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Norb,Nr
  real(c_double),intent(in),dimension(3,Nr):: rvec
  complex(c_double),intent(in),dimension(Norb,Norb,Nr):: ham_r
  real(c_double),intent(inout),dimension(Norb,Norb):: Pmn

  integer(c_int32_t) i,j,l,m
  real(c_double) diff_r
  real(c_double),parameter:: eps=1e-6
  real(c_double),dimension(Norb,Norb):: parity_den,parity_num
  
  !$omp parallel
  !$omp workshare
  parity_num(:,:)=0.0d0
  parity_den(:,:)=0.0d0
  !$omp end workshare
  !$omp do private(i,j,diff_r) reduction(+: parity_num,parity_den)
  ! Find inversion-partner pairs: R_j = -R_i (i.e. R_i + R_j = 0)
  ! Compare H(R) and H(-R): if Re(<H(R),H(-R)>) > 0 → even parity (+1), else odd (-1)
  do i=1,Nr
     do j=i+1,Nr
        diff_r=sum(abs(rvec(:,i)+rvec(:,j)))
        if(diff_r<eps)then      ! found partner: rvec(:,j) = -rvec(:,i)
           do l=1,Norb
              do m=1,Norb
                 parity_num(m,l)=parity_num(m,l)+dble(ham_r(m,l,j)*conjg(ham_r(m,l,i)))
                 parity_den(m,l)=parity_den(m,l)+abs(ham_r(m,l,i))**2
              end do
           end do
           exit
        end if
     end do
  end do
  !$omp end do
  !$omp do private(l,m)
  do l=1,Norb
     do m=1,Norb
        if(parity_den(m,l)>1e-10)then
           Pmn(m,l)=parity_num(m,l)/parity_den(m,l)
        else
           Pmn(m,l)=0.0d0
        end if
     end do
  end do
  !$omp end do
  !$omp do private(l,m)
  do l=1,Norb
     do m=1,Norb
        if(Pmn(m,l)>=0.0d0)then
           Pmn(m,l)=1.0d0
        else
           Pmn(m,l)=-1.0d0
        end if
     end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine get_parity_prop
