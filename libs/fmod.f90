module constant
  implicit none
  real(8),parameter:: pi=3.141592653589793238462643383279d0
end module constant

subroutine FFT(cmat,tmp,Nx,Ny,Nz,Nw,SW)
  implicit none
  integer(8),intent(in):: Nx,Ny,Nz,Nw
  logical(4),intent(in):: SW
  complex(8),intent(inout),dimension(Nx,Ny,Nz,Nw):: cmat,tmp

  integer(8) plan
  integer(4) Inv
  integer(4),dimension(4):: Nlist
  
  Nlist=(/Nx,Ny,Nx,Nw/)
  if(SW)then
     Inv=-1
  else
     Inv=1
  end If
  call dfftw_plan_dft(plan,4,Nlist,cmat,tmp,Inv,64)
  call dfftw_execute(plan)
  call dfftw_destroy_plan(plan)
  if(.not. SW) cmat=tmp/product(Nlist)
end subroutine FFT

subroutine openmp_params(omp_num,omp_check) bind(C)
  !$ use omp_lib
  implicit none
  integer(8),intent(out):: omp_num
  logical(1),intent(out):: omp_check
  !$ if(.true.)then
  !$   omp_check=.true.
  !$   omp_num=omp_get_max_threads()
  !$ else
       omp_check=.false.
       omp_num=0
  !$ end if
end subroutine openmp_params

subroutine gen_ham(ham_k,klist,ham_r,rvec,Nk,Nr,norb) bind(C)
  use constant
  implicit none
  integer(8),intent(in):: Nk,Nr,norb
  real(8),intent(in),dimension(3,Nk):: klist
  real(8),intent(in),dimension(3,Nr):: rvec
  complex(8),intent(in),dimension(norb,norb,Nr):: ham_r
  complex(8),intent(out),dimension(norb,norb,Nk):: ham_k

  integer(8) i,j,l,m
  real(8) phase

  !$omp parallel do private(l,m,j,phase)
  klop: do i=1,Nk
     do l=1,norb
        do m=l,norb
           !$omp simd
           rloop: do j=1,Nr
              phase=2*pi*sum(klist(:,i)*rvec(:,j))
              ham_k(m,l,i)=ham_k(m,l,i)+ham_r(m,l,j)*cmplx(cos(phase),-sin(phase))
           end do rloop
           !$omp end simd
           ham_k(l,m,i)=conjg(ham_k(m,l,i))
        end do
     end do
  end do klop
  !$omp end parallel do
end subroutine gen_ham

subroutine get_eig(eig,uni,ham_k,Nk,norb) bind(C)
  implicit none
  integer(8),intent(in):: Nk,norb
  complex(8),intent(in),dimension(norb,norb,Nk):: ham_k
  real(8),intent(out),dimension(norb,Nk):: eig
  complex(8),intent(out),dimension(norb,norb,Nk):: uni

  integer(8) i,info
  real(8) rwork(3*norb-2),eq(norb)
  complex(8) work(2*norb-1),en(norb,norb)

  !$omp parallel do private(en,eq,work,rwork,info)
  kloop: do i=1,Nk
     en(:,:)=ham_k(:,:,i)
     call zheev('V','U',norb,en,norb,eq,work,2*norb-1,rwork,info)
     uni(:,:,i)=en(:,:)
     eig(:,i)=eq(:)
  end do kloop
  !$omp end parallel do
end subroutine get_eig

subroutine get_imass0(imk,klist,ham_r,rvec,Nk,Nr,norb) bind(C)
  use constant
  implicit none
  integer(8),intent(in):: Nk,Nr,norb
  real(8),intent(in),dimension(3,Nk):: klist
  real(8),intent(in),dimension(3,Nr):: rvec
  complex(8),intent(in),dimension(norb,norb,Nr):: ham_r
  complex(8),intent(out),dimension(3,3,norb,norb,Nk):: imk

  integer(8) i,j,k,l,m,n
  real(8) phase

  !$omp parallel do private(l,m,j,phase)
  kloop: do i=1,Nk
     do l=1,norb
        do m=l,norb
           !$omp simd
           rloop: do j=1,Nr
              phase=2*pi*sum(klist(:,i)*rvec(:,j))
              axis1: do k=1,3
                 axis2: do n=k,3
                    imk(n,k,m,l,i)=imk(n,k,m,l,i)-rvec(n,j)*rvec(k,j)*ham_r(m,l,j)*cmplx(cos(phase),-sin(phase))
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

subroutine get_imassk(imk,imk0,mrot,uni,Nk,norb) bind(C)
  implicit none
  integer(8),intent(in):: Nk,norb
  real(8),intent(in),dimension(3,3):: mrot
  complex(8),intent(in),dimension(3,3,norb,norb,Nk):: imk0
  complex(8),intent(in),dimension(norb,norb,Nk):: uni
  real(8),intent(out),dimension(3,3,norb,Nk):: imk

  integer(8) i,j,k,l,m,n
  complex(8) tmp(3,3,norb)

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

subroutine get_vlm0(vk,klist,ham_r,rvec,Nk,Nr,norb) bind(C)
  use constant
  implicit none
  integer(8),intent(in):: Nk,Nr,norb
  real(8),intent(in),dimension(3,Nk):: klist
  real(8),intent(in),dimension(3,Nr):: rvec
  complex(8),intent(in),dimension(norb,norb,Nr):: ham_r  
  complex(8),intent(out),dimension(3,norb,norb,Nk):: vk

  integer(8) i,j,k,l,m
  real(8) phase

  !$omp parallel do private(l,m,j,phase)
  kloop: do i=1,Nk
     do l=1,norb
        do m=l,norb
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

subroutine get_veloc(vk,vk0,mrot,uni,Nk,norb) bind(C)
  implicit none
  integer(8),intent(in):: Nk,norb
  real(8),intent(in),dimension(3,3):: mrot
  complex(8),intent(in),dimension(3,norb,norb,Nk):: vk0
  complex(8),intent(in),dimension(norb,norb,Nk):: uni
  real(8),intent(out),dimension(3,norb,Nk):: vk

  integer(8) i,j,l,m,n
  complex(8) tmp(3,norb)

  !$omp parallel do private(tmp,l,m,n,j)
  kloop: do i=1,Nk
     !rotate orb to band
     tmp(:,:)=0.0d0
     do l=1,norb
        do m=1,norb
           band_loop: do n=1,norb
              do j=1,3
                 tmp(j,n)=tmp(j,n)+conjg(uni(l,n,i))*vk0(j,l,m,i)*uni(m,n,i)
              end do
           end do band_loop
        end do
     end do
     !rotate axis
     band_loop2: do n=1,norb
        do l=1,3
           do m=1,3
              vk(m,n,i)=vk(m,n,i)+mrot(m,l)*dble(tmp(l,n))
           end do
        end do
     end do band_loop2
  end do kloop
  !$omp end parallel do
end subroutine get_veloc

subroutine get_vnm(vk,vk0,mrot,uni,Nk,norb) bind(C)
  implicit none
  integer(8),intent(in):: Nk,norb
  real(8),intent(in),dimension(3,3):: mrot
  complex(8),intent(in),dimension(3,norb,norb,Nk):: vk0
  complex(8),intent(in),dimension(norb,norb,Nk):: uni
  complex(8),intent(out),dimension(3,norb,norb,Nk):: vk

  integer(8) i,j,l,m,n,k
  complex(8) tmp(3,norb,norb)

  !$omp parallel
  !$omp workshare
  vk(:,:,:,:)=0.0d0
  !$omp end workshare
  !$omp do private(tmp,l,m,n,j,k)
  kloop: do i=1,Nk
     !rotate orb to band
     tmp(:,:,:)=0.0d0
     orb_loop: do l=1,norb
        orb_loop2: do m=1,norb
           band_loop: do n=1,norb
              band_loop2: do k=1,norb
                 do j=1,3
                    tmp(j,k,n)=tmp(j,k,n)+conjg(uni(m,k,i))*vk0(j,m,l,i)*uni(l,n,i)
                 end do
              end do band_loop2
           end do band_loop
        end do orb_loop2
     end do orb_loop
     !rotate axis
     band_loop3: do n=1,norb
        do k=1,norb
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

subroutine gen_green0(Gk,eig,uni,mu,temp,Nk,Nw,norb) bind(C)
  use constant
  implicit none
  integer(8),intent(in):: Nk,Nw,norb
  real(8),intent(in):: mu,temp
  real(8),intent(in),dimension(norb,Nk):: eig
  complex(8),intent(in),dimension(norb,norb,Nk):: uni
  complex(8),intent(out),dimension(Nk,Nw,norb,norb):: Gk
  
  integer(8) i,j,l,m,n
  complex(8) iw

  do l=1,norb
     do m=1,norb
        band_loop: do n=1,norb
           !$omp parallel do private(iw,i,j)
           wloop: do j=1,Nw
              iw=cmplx(mu,dble(2*(j-1)+1)*pi*temp)
              !!$omp simd
              kloop: do i=1,Nk
                 Gk(i,j,m,l)=Gk(i,j,m,l)+uni(l,n,i)*conjg(uni(m,n,i))/(iw-eig(n,i))
              end do kloop
              !!$omp end simd
           end do wloop
           !$omp end parallel do
        end do band_loop
     end do
  end do
end subroutine gen_green0

subroutine get_chi0_comb(chi,Gk,kmap,olist,Nx,Ny,Nz,Nw,Nk,Norb,Nchi) bind(C)
  implicit none
  integer(8),intent(in):: Nx,Ny,Nz,Nw,Norb,Nchi,Nk
  integer(8),intent(in),dimension(2,Nchi):: olist
  integer(8),intent(in),dimension(3,Nk):: kmap
  complex(8),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(8),intent(out),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(8) i,j,k,l,m,n
  integer(8) ii(0:Nx-1),ij(0:Ny-1),ik(0:Nz-1)
  complex(8),dimension(0:Nx-1,0:Ny-1,0:Nz-1,2*Nw):: tmp1,tmp2,tmpchi,tmp

  do i=0,Nx-1
     ii(i)=mod(Nx-i,Nx)
  end do
  
  do j=0,Ny-1
     ij(j)=mod(Ny-j,Ny)
  end do
  
  do k=0,Nz-1
     ik(k)=mod(Nz-k,Nz)
  end do

  do l=1,Nchi
     do m=1,Nchi
        !use symmetry G^lm(k,iw)=G^ml(k,-iw)
        !$omp parallel do private(i)
        do j=1,Nw
           do i=1,Nk
              tmp1(kmap(1,i),kmap(2,i),kmap(3,i),j)=Gk(i,j,olist(1,l),olist(2,m)) !G13(iw)
              tmp2(kmap(1,i),kmap(2,i),kmap(3,i),j)=Gk(i,j,olist(2,l),olist(1,m)) !G42(iw)
              tmp1(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=conjg(Gk(i,j,olist(2,m),olist(1,l))) !G13(-iw)
              tmp2(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=conjg(Gk(i,j,olist(1,m),olist(2,l))) !G42(-iw)
           end do
        end do
        !$omp end parallel do
        call FFT(tmp1,tmp,Nx,Ny,Nz,2*Nw,.true.)
        call FFT(tmp2,tmp,Nx,Ny,Nz,2*Nw,.true.)
        !calculate G(r)G(-r)
        !$omp parallel do private(i,j,k)
        do n=1,2*Nw
           do k=0,Nz-1
              do j=0,Ny-1
                 !$omp simd
                 do i=0,Nx-1
                    tmpchi(i,j,k,n)=tmp1(i,j,k,n)*tmp2(ii(i),ij(j),ik(k),2*nw-n+1)
                 end do
                 !$omp end simd
              end do
           end do
        end do
        !$omp end parallel do
        call FFT(tmpchi,tmp,Nx,Ny,Nz,2*Nw,.false.)
        !$omp parallel do private(i,j)
        do j=1,Nw
           do i=1,Nk
              chi(i,j,m,l)=tmpchi(kmap(1,i),kmap(2,i),kmap(3,i),j)
           end do
        end do
        !$omp end parallel do
     end do
  end do
end subroutine get_chi0_comb

subroutine gen_tr_greenw_0(trGk,wl,eig,mu,delta,Nk,Nw,Norb) bind(C)
  implicit none
  integer(8),intent(in):: Nk,Nw,Norb
  real(8),intent(in):: mu,delta
  real(8),intent(in),dimension(Norb,Nk):: eig
  real(8),intent(in),dimension(Nw):: wl
  complex(8),intent(out),dimension(Nw,Nk):: trGk

  integer(8) i,j,n

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

subroutine get_smat(Smat,ol,Uval,Jval,Nchi,Norb) bind(C)
  implicit none
  integer(8),intent(in):: Nchi,Norb
  integer(8),intent(in),dimension(2,Nchi):: ol
  real(8),intent(in):: Uval,Jval
  real(8),intent(out),dimension(Nchi,Nchi):: Smat

  integer(8) i,j

  !$omp parallel
  !$omp workshare
  Smat(:,:)=0.0d0
  !$omp end workshare
  !$omp do private(j)
  do i=1,Nchi
     do j=1,Nchi
        if((ol(1,i)==ol(2,i)).and.(ol(1,j)==ol(2,j)))then
           if(ol(1,i)==ol(1,j))then
              Smat(j,i)=Uval
           else
              Smat(j,i)=Jval
           end if
        else if((ol(1,i)==ol(1,j)).and.(ol(2,i)==ol(2,j)))then
           Smat(j,i)=Uval-2*Jval
        else if((ol(1,i)==ol(2,j)).and.(ol(2,i)==ol(1,j)))then
           Smat(j,i)=Jval
        end if
     end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine get_smat

