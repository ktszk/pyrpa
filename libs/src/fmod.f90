module constant
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
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
  integer(int64),intent(in):: Nk,Nr,Norb
  real(real64),intent(in),dimension(3,Nk):: klist
  real(real64),intent(in),dimension(3,Nr):: rvec
  complex(real64),intent(in),dimension(Norb,Norb,Nr):: ham_r
  complex(real64),intent(out),dimension(Norb,Norb,Nk):: ham_k

  integer(int32) i,j,l,m
  real(real64) phase

  !$omp parallel do private(l,m,j,phase)
  klop: do i=1,Nk
     do l=1,Norb
        do m=l,Norb
           rloop: do j=1,Nr
              phase=2*pi*sum(klist(:,i)*rvec(:,j))
              ham_k(m,l,i)=ham_k(m,l,i)+ham_r(m,l,j)*cmplx(cos(phase),-sin(phase))
           end do rloop
           if(l==m)then
              ham_k(l,l,i)=dble(ham_k(l,l,i)) !diagonal is real
           else
              ham_k(l,m,i)=conjg(ham_k(m,l,i)) !Hamiltonian is Hermite
           end if
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
  use,intrinsic:: iso_fortran_env, only:int32,int64,real64
  implicit none
  integer(int64),intent(in):: Nk,Norb
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: ham_k
  real(real64),intent(out),dimension(Norb,Nk):: eig
  complex(real64),intent(out),dimension(Norb,Norb,Nk):: uni

  integer(int32) i,info
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

subroutine get_eig_mlo(eig,uni,ham_k,Ovlk,Nk,Norb) bind(C)
  !> This function obtains energies
  !!@param  eig,out: energies at k-points
  !!@param  uni,out: unitary matrix at k-points
  !!@param ham_k,in: k-space Hamiltonian
  !!@param  Ovlk,in: k-space Overlap integrals
  !!@param    Nk,in: The number of k-points
  !!@param  Norb,in: The number of orbitals
  use,intrinsic:: iso_fortran_env, only:int32,int64,real64
  implicit none
  integer(int64),intent(in):: Nk,norb
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: ham_k,Ovlk
  real(real64),intent(out),dimension(Norb,Nk):: eig
  complex(real64),intent(out),dimension(Norb,Norb,Nk):: uni

  integer(int32) i,j,k,l,m,info
  real(real64) rwork(3*Norb-2),eq(Norb),norm
  complex(real64) work(2*Norb-1)
  complex(real64),dimension(Norb,Norb):: tmp,tmp2,tmp3

  !$omp parallel do private(tmp,tmp2,tmp3,norm,eq,work,rwork,info)
  kloop: do i=1,Nk
     tmp(:,:)=Ovlk(:,:,i)
     call zheev('V','U',norb,tmp,norb,eq,work,2*norb-1,rwork,info)
     !$omp simd
     do j=1,Norb
        tmp2(:,j)=tmp(:,j)/sqrt(cmplx(eq(j)))
     end do
     !$omp end simd
     tmp(:,:)=0.0d0
     !$omp simd
     do j=1,Norb
        do k=1,Norb
           do l=1,Norb
              do m=1,Norb
                 tmp(k,j)=tmp(k,j)+conjg(tmp2(m,k))*ham_k(m,l,i)*tmp2(l,j)
              end do
           end do
        end do
     end do
     !$omp end simd
     call zheev('V','U',norb,tmp,norb,eq,work,2*norb-1,rwork,info)
     eig(:,i)=eq(:)
     tmp3(:,:)=0.0d0
     do j=1,Norb
        do k=1,Norb
           do l=1,Norb
              tmp3(k,j)=tmp3(k,j)+tmp2(k,l)*tmp(l,j)
           end do
        end do
     end do
     uni(:,:,i)=tmp3(:,:)
  end do kloop
  !$omp end parallel do
end subroutine get_eig_mlo

subroutine get_ffermi(ffermi,eig,mu,temp,Nk,Norb) bind(C)
  !> This function obtains fermi functions
  !!@param ffermi,out: fermi distribute functions
  !!@param     eig,in: energies at k-points
  !!@param      mu,in: chemical potential
  !!@param    temp,in: Temperature
  !!@param      Nk,in: The number of k-points
  !!@param    Norb,in: The number of orbitals
  use,intrinsic:: iso_fortran_env, only:int32,int64,real64
  implicit none
  integer(int64),intent(in):: Nk,Norb
  real(real64),intent(in):: mu,temp
  real(real64),intent(in),dimension(Norb,Nk):: eig
  real(real64),intent(out),dimension(Norb,Nk):: ffermi

  integer(int32) i,j
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
  use constant
  implicit none
  integer(int64),intent(in):: Nk,Nr,norb
  real(real64),intent(in),dimension(3,Nk):: klist
  real(real64),intent(in),dimension(3,Nr):: rvec
  complex(real64),intent(in),dimension(Norb,Norb,Nr):: ham_r
  complex(real64),intent(out),dimension(3,3,Norb,Norb,Nk):: imk

  integer(int32) i,j,k,l,m,n
  real(real64) phase

  !$omp parallel do private(l,m,j,phase)
  kloop: do i=1,Nk
     do l=1,Norb
        do m=l,Norb
           rloop: do j=1,Nr
              phase=2*pi*sum(klist(:,i)*rvec(:,j))
              !$omp simd
              axis1: do k=1,3
                 axis2: do n=k,3
                    imk(n,k,m,l,i)=imk(n,k,m,l,i)-rvec(n,j)*rvec(k,j)&
                         *ham_r(m,l,j)*cmplx(cos(phase),-sin(phase))
                 end do axis2
              end do axis1
              !$omp end simd
           end do rloop
           !$omp simd
           axis12: do k=1,3
              axis22: do n=k,3
                 imk(k,n,m,l,i)=imk(n,k,l,m,i)
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
  use,intrinsic:: iso_fortran_env, only:int32,int64,real64
  implicit none
  integer(int64),intent(in):: Nk,Norb
  real(real64),intent(in),dimension(3,3):: mrot
  complex(real64),intent(in),dimension(3,3,Norb,Norb,Nk):: imk0
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  real(real64),intent(out),dimension(3,3,Norb,Nk):: imk

  integer(int32) i,j,k,l,m,n
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
  integer(int64),intent(in):: Nk,Nr,Norb
  real(real64),intent(in),dimension(3,Nk):: klist
  real(real64),intent(in),dimension(3,Nr):: rvec
  complex(real64),intent(in),dimension(Norb,Norb,Nr):: ham_r
  complex(real64),intent(out),dimension(3,Norb,Norb,Nk):: vk

  integer(int32) i,j,k,l,m
  real(real64) phase

  !$omp parallel do private(l,m,j,phase)
  kloop: do i=1,Nk
     do l=1,Norb
        do m=l,Norb
           rloop: do j=1,Nr
              phase=2*pi*sum(klist(:,i)*rvec(:,j))
              !$omp simd
              vaxis: do k=1,3
                 vk(k,m,l,i)=vk(k,m,l,i)-rvec(k,j)*ham_r(m,l,j)*cmplx(sin(phase),cos(phase))
              end do vaxis
              !$omp end simd
           end do rloop
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
  use,intrinsic:: iso_fortran_env, only:int32,int64,real64
  implicit none
  integer(int64),intent(in):: Nk,Norb
  real(real64),intent(in),dimension(3,3):: mrot
  complex(real64),intent(in),dimension(3,Norb,Norb,Nk):: vk0
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  real(real64),intent(out),dimension(3,Norb,Nk):: vk

  integer(int32) i,j,l,m,n
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
  !> This function obtain band-basis group velocity with inter-band
  !!@param  vk,out: band-basis group velocities
  !!@param  vk0,in: orbital-basis group velocities
  !!@param mrot,in: rotation matrix for real space
  !!@param  uni,in: unitary matrix orbital to band
  !!@param   Nk,in: The number of k-points
  !!@param Norb,in: The number of orbitals
  use,intrinsic:: iso_fortran_env, only:int32,int64,real64
  implicit none
  integer(int64),intent(in):: Nk,Norb
  real(real64),intent(in),dimension(3,3):: mrot
  complex(real64),intent(in),dimension(3,Norb,Norb,Nk):: vk0
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(real64),intent(out),dimension(3,Norb,Norb,Nk):: vk

  integer(int32) i,j,l,m,n,k
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
  !> This function obtain orbital(band) trace of green function
  !!@param trGk,out: trace of green function
  !!@param    wl,in: list of energies
  !!@param   eig,in: energies at k-points
  !!@param    mu,in: chemical potential
  !!@param delta,in: dumping factor
  !!@param    Nk,in: The number of k-points
  !!@param    Nw,in: The number of energies mesh
  !!@param  Norb,in: The number of orbitals
  use,intrinsic:: iso_fortran_env, only:int32,int64,real64
  implicit none
  integer(int64),intent(in):: Nk,Nw,Norb
  real(real64),intent(in):: mu,delta
  real(real64),intent(in),dimension(Norb,Nk):: eig
  real(real64),intent(in),dimension(Nw):: wl
  complex(real64),intent(out),dimension(Nw,Nk):: trGk

  integer(int32) i,j,n

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
  use,intrinsic:: iso_fortran_env, only:int32,int64,real64
  implicit none
  integer(int64),intent(in):: Nk,Nw,Norb
  real(real64),intent(in):: mu,delta
  real(real64),intent(in),dimension(Norb,Nk):: eig
  real(real64),intent(in),dimension(Nw):: wl
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(real64),intent(out),dimension(Nw,Norb):: Dos

  integer(int32) i,j,k,n

  !$omp parallel
  !$omp workshare
  Dos(:,:)=0.0d0
  !$omp end workshare
  orb_loop: do j=1,Norb
     !$omp do private(i,k,n)
     wloop: do i=1,Nw
        kloop: do k=1,Nk
           bandloop: do n=1,Norb
              Dos(i,j)=Dos(i,j)+uni(j,n,k)*conjg(uni(j,n,k))/cmplx(wl(i)-eig(n,k)+mu,delta)
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
