module constant
  use,intrinsic:: iso_fortran_env, only:real64
  implicit none
  real(real64),parameter:: pi=3.141592653589793238462643383279d0
end module constant

subroutine openmp_params(omp_num,omp_check) bind(C)
  !$ use omp_lib
  use,intrinsic:: iso_fortran_env, only:int64
  implicit none
  integer(int64),intent(out):: omp_num
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
  use constant
  use,intrinsic:: iso_fortran_env, only:int64,real64
  implicit none
  integer(int64),intent(in):: Nk,Nr,Norb
  real(real64),intent(in),dimension(3,Nk):: klist
  real(real64),intent(in),dimension(3,Nr):: rvec
  complex(real64),intent(in),dimension(Norb,Norb,Nr):: ham_r
  complex(real64),intent(out),dimension(Norb,Norb,Nk):: ham_k

  integer(int64) i,j,l,m
  real(real64) phase

  !$omp parallel do private(l,m,j,phase)
  klop: do i=1,Nk
     do l=1,Norb
        do m=l,Norb
           !$omp simd
           rloop: do j=1,Nr
              phase=2*pi*sum(klist(:,i)*rvec(:,j))
              ham_k(m,l,i)=ham_k(m,l,i)+ham_r(m,l,j)*cmplx(cos(phase),-sin(phase))
           end do rloop
           !$omp end simd
           ham_k(l,m,i)=conjg(ham_k(m,l,i)) !Hamiltonian is Hermite
        end do
     end do
  end do klop
  !$omp end parallel do
end subroutine gen_ham

subroutine get_eig(eig,uni,ham_k,Nk,Norb) bind(C)
  use,intrinsic:: iso_fortran_env, only:int64,real64
  implicit none
  integer(int64),intent(in):: Nk,Norb
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: ham_k
  real(real64),intent(out),dimension(Norb,Nk):: eig
  complex(real64),intent(out),dimension(Norb,Norb,Nk):: uni

  integer(int64) i,info
  real(real64) rwork(3*Norb-2),eq(Norb)
  complex(real64) work(2*Norb-1),en(Norb,Norb)

  !$omp parallel do private(en,eq,work,rwork,info)
  kloop: do i=1,Nk
     en(:,:)=ham_k(:,:,i)
     call zheev('V','U',Norb,en,Norb,eq,work,2*Norb-1,rwork,info)
     uni(:,:,i)=en(:,:)
     eig(:,i)=eq(:)
  end do kloop
  !$omp end parallel do
end subroutine get_eig

subroutine get_eig_mlo(eig,uni,ham_k,Ovlk,Nk,norb) bind(C)
  use,intrinsic:: iso_fortran_env, only:int64,real64
  implicit none
  integer(int64),intent(in):: Nk,norb
  complex(real64),intent(in),dimension(norb,norb,Nk):: ham_k,Ovlk
  real(real64),intent(out),dimension(norb,Nk):: eig
  complex(real64),intent(out),dimension(norb,norb,Nk):: uni

  integer(int64) i,j,k,l,m,info
  real(real64) rwork(3*norb-2),eq(norb),norm
  complex(real64) work(2*norb-1)
  complex(real64),dimension(norb,norb):: tmp,tmp2,tmp3

  !$omp parallel do private(tmp,tmp2,tmp3,norm,eq,work,rwork,info)
  kloop: do i=1,Nk
     tmp(:,:)=Ovlk(:,:,i)
     call zheev('V','U',norb,tmp,norb,eq,work,2*norb-1,rwork,info)
     do j=1,norb
        tmp2(:,j)=tmp(:,j)/sqrt(cmplx(eq(j)))
     end do
     tmp(:,:)=0.0d0
     do j=1,norb
        do k=1,norb
           do l=1,norb
              do m=1,norb
                 tmp(k,j)=tmp(k,j)+conjg(tmp2(m,k))*ham_k(m,l,i)*tmp2(l,j)
              end do
           end do
        end do
     end do
     call zheev('V','U',norb,tmp,norb,eq,work,2*norb-1,rwork,info)
     eig(:,i)=eq(:)
     tmp3(:,:)=0.0d0
     do j=1,norb
        do k=1,norb
           do l=1,norb
              tmp3(k,j)=tmp3(k,j)+tmp2(k,l)*tmp(l,j)
           end do
        end do
     end do
     uni(:,:,i)=tmp3(:,:)
  end do kloop
  !$omp end parallel do
end subroutine get_eig_mlo

subroutine get_ffermi(ffermi,eig,mu,temp,Nk,Norb) bind(C)
  use,intrinsic:: iso_fortran_env, only:int64,real64
  implicit none
  integer(int64),intent(in):: Nk,Norb
  real(real64),intent(in):: mu,temp
  real(real64),intent(in),dimension(Norb,Nk):: eig
  real(real64),intent(out),dimension(Norb,Nk):: ffermi

  integer(int64) i,j
  real(real64) itemp

  itemp=0.5d0/temp
  !$omp parallel do private(j)
  do i=1,Nk
     do j=1,Norb
        ffermi(j,i)=0.5d0-0.5d0*tanh((eig(j,i)-mu)*itemp)
     end do
  end do
  !$omp end parallel do
end subroutine get_ffermi

subroutine get_imass0(imk,klist,ham_r,rvec,Nk,Nr,Norb) bind(C)
  use,intrinsic:: iso_fortran_env, only:int64,real64
  use constant
  implicit none
  integer(int64),intent(in):: Nk,Nr,norb
  real(real64),intent(in),dimension(3,Nk):: klist
  real(real64),intent(in),dimension(3,Nr):: rvec
  complex(real64),intent(in),dimension(Norb,Norb,Nr):: ham_r
  complex(real64),intent(out),dimension(3,3,Norb,Norb,Nk):: imk

  integer(int64) i,j,k,l,m,n
  real(real64) phase

  !$omp parallel do private(l,m,j,phase)
  kloop: do i=1,Nk
     do l=1,Norb
        do m=l,Norb
           !$omp simd
           rloop: do j=1,Nr
              phase=2*pi*sum(klist(:,i)*rvec(:,j))
              axis1: do k=1,3
                 axis2: do n=k,3
                    imk(n,k,m,l,i)=imk(n,k,m,l,i)-rvec(n,j)*rvec(k,j)&
                         *ham_r(m,l,j)*cmplx(cos(phase),-sin(phase))
                 end do axis2
              end do axis1
           end do rloop
           !$omp end simd
           axis12: do k=1,3
              axis22: do n=k,3
                 imk(k,n,m,l,i)=imk(n,k,l,m,i)
                 imk(n,k,l,m,i)=conjg(imk(n,k,m,l,i))
                 imk(k,n,m,l,i)=imk(n,k,l,m,i)
              end do axis22
           end do axis12
        end do
     end do
  end do kloop
  !$omp end parallel do
end subroutine get_imass0

subroutine get_imassk(imk,imk0,mrot,uni,Nk,Norb) bind(C)
  use,intrinsic:: iso_fortran_env, only:int64,real64
  implicit none
  integer(int64),intent(in):: Nk,Norb
  real(real64),intent(in),dimension(3,3):: mrot
  complex(real64),intent(in),dimension(3,3,Norb,Norb,Nk):: imk0
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  real(real64),intent(out),dimension(3,3,Norb,Nk):: imk

  integer(int64) i,j,k,l,m,n
  complex(real64) tmp(3,3,Norb)

  !$omp parallel do private(tmp,l,m,n,j)
  kloop: do i=1,Nk
     !rotate orb to band
     tmp(:,:,:)=0.0d0
     do l=1,norb
        do m=1,norb
           band_loop: do n=1,norb
              do j=1,3
                 do k=1,3 
                    tmp(k,j,n)=tmp(k,j,n)+conjg(uni(l,n,i))*imk0(k,j,l,m,i)*uni(m,n,i)
                 end do
              end do
           end do band_loop
        end do
     end do
     !rotate axis
     band_loop2: do n=1,norb
        do l=1,3
           do m=1,3
              do j=1,3
                 do k=1,3
                    imk(m,l,n,i)=imk(m,l,n,i)+mrot(m,k)*mrot(l,j)*dble(tmp(k,j,n))
                 end do
              end do
           end do
        end do
     end do band_loop2
  end do kloop
  !$omp end parallel do
end subroutine get_imassk

subroutine get_vlm0(vk,klist,ham_r,rvec,Nk,Nr,Norb) bind(C)
  use,intrinsic:: iso_fortran_env, only:int64,real64
  use constant
  implicit none
  integer(int64),intent(in):: Nk,Nr,Norb
  real(real64),intent(in),dimension(3,Nk):: klist
  real(real64),intent(in),dimension(3,Nr):: rvec
  complex(real64),intent(in),dimension(Norb,Norb,Nr):: ham_r
  complex(real64),intent(out),dimension(3,Norb,Norb,Nk):: vk

  integer(int64) i,j,k,l,m
  real(real64) phase

  !$omp parallel do private(l,m,j,phase)
  kloop: do i=1,Nk
     do l=1,Norb
        do m=l,Norb
           !$omp simd
           rloop: do j=1,Nr
              phase=2*pi*sum(klist(:,i)*rvec(:,j))
              vaxis: do k=1,3
                 vk(k,m,l,i)=vk(k,m,l,i)-rvec(k,j)*ham_r(m,l,j)*cmplx(sin(phase),cos(phase))
              end do vaxis
           end do rloop
           !$omp end simd
           vaxis2: do k=1,3
              vk(k,l,m,i)=conjg(vk(k,m,l,i))
           end do vaxis2
        end do
     end do
  end do kloop
  !$omp end parallel do
end subroutine get_vlm0

subroutine get_veloc(vk,vk0,mrot,uni,Nk,Norb) bind(C)
  use,intrinsic:: iso_fortran_env, only:int64,real64
  implicit none
  integer(int64),intent(in):: Nk,Norb
  real(real64),intent(in),dimension(3,3):: mrot
  complex(real64),intent(in),dimension(3,Norb,Norb,Nk):: vk0
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  real(real64),intent(out),dimension(3,Norb,Nk):: vk

  integer(int64) i,j,l,m,n
  complex(real64) tmp(3,Norb)

  !$omp parallel do private(tmp,l,m,n,j)
  kloop: do i=1,Nk
     !rotate orb to band
     tmp(:,:)=0.0d0
     do l=1,Norb
        do m=1,Norb
           band_loop: do n=1,Norb
              do j=1,3
                 tmp(j,n)=tmp(j,n)+conjg(uni(l,n,i))*vk0(j,l,m,i)*uni(m,n,i)
              end do
           end do band_loop
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
  use,intrinsic:: iso_fortran_env, only:int64,real64
  implicit none
  integer(int64),intent(in):: Nk,Norb
  real(real64),intent(in),dimension(3,3):: mrot
  complex(real64),intent(in),dimension(3,Norb,Norb,Nk):: vk0
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(real64),intent(out),dimension(3,Norb,Norb,Nk):: vk

  integer(int64) i,j,l,m,n,k
  complex(real64) tmp(3,Norb,Norb)

  !$omp parallel
  !$omp workshare
  vk(:,:,:,:)=0.0d0
  !$omp end workshare
  !$omp do private(tmp,l,m,n,j,k)
  kloop: do i=1,Nk
     !rotate orb to band
     tmp(:,:,:)=0.0d0
     orb_loop: do l=1,Norb
        orb_loop2: do m=1,Norb
           band_loop: do n=1,Norb
              band_loop2: do k=1,Norb
                 do j=1,3
                    tmp(j,k,n)=tmp(j,k,n)+conjg(uni(m,k,i))*vk0(j,m,l,i)*uni(l,n,i)
                 end do
              end do band_loop2
           end do band_loop
        end do orb_loop2
     end do orb_loop
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
  use,intrinsic:: iso_fortran_env, only:int64,real64
  implicit none
  integer(int64),intent(in):: Nk,Nw,Norb
  real(real64),intent(in):: mu,delta
  real(real64),intent(in),dimension(Norb,Nk):: eig
  real(real64),intent(in),dimension(Nw):: wl
  complex(real64),intent(out),dimension(Nw,Nk):: trGk

  integer(int64) i,j,n

  !$omp parallel do private(i,n)
  kloop: do j=1,Nk
     !$omp simd
     wloop: do i=1,Nw
        band_loop: do n=1,Norb
           trGk(i,j)=trGk(i,j)+1./cmplx(wl(i)-eig(n,j)+mu,delta)
        end do band_loop
     end do wloop
     !$omp end simd
  end do kloop
  !$omp end parallel do
end subroutine gen_tr_greenw_0
