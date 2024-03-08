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
  use constant
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
        !call FFT(tmp1,tmp,Nx,Ny,Nz,2*Nw,.true.)
        !call FFT(tmp2,tmp,Nx,Ny,Nz,2*Nw,.true.)
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
        !call FFT(tmpchi,tmp,Nx,Ny,Nz,2*Nw,.false.)
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
  use constant
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

subroutine get_qshift(qpoint,klist,qshift,Nk) bind(C)
  implicit none
  integer(8),intent(in):: Nk
  real(8),intent(in),dimension(3,Nk):: klist
  real(8),intent(in),dimension(3):: qpoint
  integer(8),intent(out),dimension(Nk):: qshift

  call set_qshift(qpoint,klist,qshift,Nk)
end subroutine get_qshift

subroutine set_qshift(qpoint,klist,qshift,Nk)
  implicit none
  integer(8),intent(in):: Nk
  real(8),intent(in),dimension(3,Nk):: klist
  real(8),intent(in),dimension(3):: qpoint
  integer(8),intent(out),dimension(Nk):: qshift

  integer(8) i,j,k,ck
  real(8) tmp
  real(8),dimension(3,Nk):: kqlist

  !$omp parallel
  !$omp workshare
  kqlist(:,:)=0.0d0
  !$omp end workshare
  !$omp do private(i,j)
  kloop:do i=1,Nk
     kqlist(:,i)=klist(:,i)+qpoint(:)
     do j=1,3
        if(kqlist(j,i)>=1.0d0)then
           kqlist(j,i)=kqlist(j,i)-1.0d0
        else if(kqlist(j,i)<0.0d0)then
           kqlist(j,i)=kqlist(j,i)+1.0d0           
        end if
     end do
  end do kloop
  !$omp end do
  !$omp workshare
  qshift(:)=1
  !$omp end workshare
  !$omp do private(i,j,tmp)
  kq_loop: do i=1,Nk
     k_loop: do j=1,Nk
        tmp=sum(abs(klist(:,j)-kqlist(:,i)))
        if(tmp<1.0d-9)then
           qshift(i)=j
           exit
        end if
     end do k_loop
  end do kq_loop
  !$omp end do
  !$omp end parallel
end subroutine set_qshift

module calc_irr_chi
  implicit none
contains
  function calc_chi(Nk,Norb,Nchi,uni,eig,ffermi,ol,temp,qshift,w,idelta,eps)
    implicit none
    integer(8),intent(in):: Nk,Norb,Nchi
    integer(8),intent(in),dimension(Nk):: qshift
    integer(8),intent(in),dimension(2,Nchi):: ol
    real(8),intent(in):: temp,eps,idelta,w
    real(8),intent(in),dimension(Norb,Nk):: eig,ffermi
    complex(8),intent(in),dimension(Norb,Norb,Nk):: uni

    integer(8) i,j,k,l,m
    complex(8) unitmp
    complex(8),dimension(Nchi,Nchi):: chi,calc_chi

    chi(:,:)=0.0d0
    kloop: do k=1,Nk
       band1_loop: do l=1,Norb
          band2_loop: do m=1,Norb
             chiorb1_loop: do j=1,Nchi
                chiorb2_loop:do i=1,Nchi
                   unitmp=uni(ol(1,j),l,qshift(k))*conjg(uni(ol(1,i),l,qshift(k)))&
                        *uni(ol(2,i),m,k)*conjg(uni(ol(2,j),m,k))
                   if(abs(w)==0.0d0 .and. abs(eig(m,k)-eig(l,qshift(k)))<1.0d-9)then
                      chi(i,j)=chi(i,j)+unitmp*ffermi(m,k)*(1.0d0-ffermi(m,k))/temp
                   else if(abs(ffermi(l,qshift(i))-ffermi(m,i))>eps)then
                      chi(i,j)=chi(i,j)+unitmp*(ffermi(l,qshift(k))-ffermi(m,k))&
                           /cmplx(w+eig(m,k)-eig(l,qshift(k)),idelta)
                   end if
                end do chiorb2_loop
             end do chiorb1_loop
          end do band2_loop
       end do band1_loop
    end do kloop
    calc_chi=chi(:,:)/Nk
  end function calc_chi
end module calc_irr_chi

subroutine get_tr_chi(trchis,trchi0,chis,chi0,olist,Nw,Nchi) bind(C)
  implicit none
  integer(8),intent(in):: Nchi,Nw
  integer(8),intent(in),dimension(2,Nchi):: olist
  complex(8),intent(in),dimension(Nchi,Nchi,Nw):: chis,chi0
  complex(8),intent(out),dimension(Nw):: trchis,trchi0

  integer(8) i,j,k
  !$omp parallel do private(j,k)
  wloop:do i=1,Nw
     orb_lop1:do j=1,Nchi
        if(olist(1,j)==olist(2,j))then
           orb_loop2:do k=1,Nchi
              if(olist(1,k)==olist(2,k))then
                 trchis(i)=trchis(i)+chis(k,j,i)
                 trchi0(i)=trchi0(i)+chi0(k,j,i)
              end if
           end do orb_loop2
        end if
     end do orb_lop1
  end do wloop
  !$omp end parallel do
end subroutine get_tr_chi

subroutine get_chi_irr(chi,uni,eig,ffermi,qshift,ol,wl,Nchi,Norb,Nk,Nw,idelta,eps,temp) bind(C)
  use calc_irr_chi
  implicit none
  integer(8),intent(in):: Nk,Norb,Nw,Nchi
  integer(8),intent(in),dimension(Nk):: qshift
  integer(8),intent(in),dimension(2,Nchi):: ol
  real(8),intent(in):: temp,eps,idelta
  real(8),intent(in),dimension(Norb,Nk):: eig,ffermi
  real(8),intent(in),dimension(Nw):: wl
  complex(8),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(8),intent(out),dimension(Nchi,Nchi,Nw):: chi

  integer(8) i

  !$omp parallel do private(i)
  wloop: do i=1,Nw
     chi(:,:,i)=calc_chi(Nk,Norb,Nchi,uni,eig,ffermi,ol,temp,qshift,wl(i),idelta,eps)
  end do wloop
  !$omp end parallel do
end subroutine get_chi_irr

subroutine chiq_map(trchis,trchi,uni,eig,ffermi,klist,Smat,ol,temp,ecut,idelta,eps,Nx,Ny,Nk,Norb,Nchi) bind(C)
  use calc_irr_chi
  implicit none
  integer(8),intent(in):: Nx,Ny,Nk,Norb,Nchi
  integer(8),intent(in),dimension(2,Nchi):: ol
  real(8),intent(in):: ecut,idelta,eps,temp
  real(8),intent(in),dimension(3,Nk):: klist
  real(8),intent(in),dimension(Norb,Nk):: eig,ffermi
  real(8),intent(in),dimension(Nchi,Nchi):: Smat
  complex(8),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(8),intent(out),dimension(Ny,Nx):: trchis,trchi

  integer(4) info
  integer(8) i,j,l,m,n
  integer(8),dimension(Nk):: qshift
  integer(4),dimension(Nchi):: ipiv
  real(8),dimension(3):: qpoint
  complex(8),dimension(Nchi,Nchi):: chi,tmp,tmp2
  complex(8),dimension(2*Nchi):: work

  qpoint(3)=0
  !$omp parallel
  !$omp workshare
  trchi(:,:)=0.0d0
  trchis(:,:)=0.0d0
  !$omp end workshare
  !$omp do private(i,j,l,m,n,chi,tmp,tmp2,qpoint,qshift,ipiv,work,info)
  do i=1,Nx
     do j=1,Ny
        qpoint(1)=dble(i-1)/Nx
        qpoint(2)=dble(j-1)/Ny
        call set_qshift(qpoint,klist,qshift,Nk)
        chi(:,:)=calc_chi(Nk,Norb,Nchi,uni,eig,ffermi,ol,temp,qshift,ecut,idelta,eps)
        tmp(:,:)=0.0d0
        do l=1,Nchi
           tmp(l,l)=1.0d0
           do m=1,Nchi
              !$omp simd
              do n=1,Nchi
                 tmp(m,l)=tmp(m,l)-Smat(m,n)*chi(n,l)
              end do
              !$omp end simd
           end do
        end do
        call zgetrf(Nchi,Nchi,tmp,Nchi,ipiv,info)
        call zgetri(Nchi,tmp,Nchi,ipiv,work,2*Nchi,info)
        tmp2(:,:)=0.0d0
        do l=1,Nchi
           do m=1,Nchi
              !$omp simd
              do n=1,Nchi
                 tmp2(m,l)=tmp2(m,l)+chi(m,n)*tmp(n,l)
              end do
              !$omp end simd
           end do
        end do
        !take chis_llmm
        do l=1,Nchi
           if(ol(1,l)==ol(2,l))then
              do m=1,Nchi
                 if(ol(1,m)==ol(2,m))then
                    trchis(j,i)=trchis(j,i)+tmp2(l,l)
                    trchi(j,i)=trchi(j,i)+chi(l,l)
                 end if
              end do
           end if
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine chiq_map

subroutine get_chis(chis,chi0,Smat,Nchi,Nw) bind(C)
  implicit none
  integer(8),intent(in):: Nchi,Nw
  real(8),dimension(Nchi,Nchi):: Smat
  complex(8),intent(in),dimension(Nchi,Nchi,Nw):: chi0
  complex(8),intent(out),dimension(Nchi,Nchi,Nw):: chis
  integer(8) i,l,m,n
  integer(4) info
  integer(4),dimension(Nchi):: ipiv
  complex(8),dimension(Nchi):: work
  complex(8),dimension(Nchi,Nchi):: tmp

  !$omp parallel do private(tmp,l,m,n,work,ipiv,info)
  do i=1,Nw
     tmp(:,:)=0.0d0
     do l=1,Nchi
        do n=1,Nchi
           !$omp simd
           do m=1,Nchi
              tmp(m,l)=tmp(m,l)-chi0(m,n,i)*Smat(n,l)
           end do
           !$omp end simd
        end do
        tmp(l,l)=tmp(l,l)+1.0d0
     end do
     call zgetrf(Nchi,Nchi,tmp,Nchi,ipiv,info)
     call zgetri(Nchi,tmp,Nchi,ipiv,work,Nchi,info)
     do l=1,Nchi
        do n=1,Nchi
           !$omp simd
           do m=1,Nchi
              chis(m,l,i)=chis(m,l,i)+tmp(m,n)*chi0(n,m,i)
           end do
           !$omp end simd
        end do
     end do
  end do
  !$omp end parallel do 
end subroutine get_chis

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

subroutine calc_lij(L11,L22,L12,vk,eig,ffermi,Norb,Nk,mu,w,idelta,eps,temp) bind(C)
  implicit none
  integer(8),intent(in):: Nk,Norb
  real(8),intent(in):: temp,eps,w,idelta,mu
  real(8),intent(in),dimension(Norb,Nk):: eig,ffermi
  complex(8),intent(in),dimension(3,Norb,Norb,Nk):: vk
  complex(8),intent(out),dimension(3,3):: L11,L12,L22
  
  integer(8) i,j,k,l,m
  complex(8) tmp
  complex(8),parameter::ii=(0.0d0,1.0d0)

  !$omp parallel
  !$omp workshare
  L11(:,:)=0.0d0
  L12(:,:)=0.0d0
  L22(:,:)=0.0d0
  !$omp end workshare
  !$omp do reduction(+: L11,L12,L22) private(i,l,m,j,k,tmp)
  k_loop: do i=1,Nk
     band_loop1: do l=1,Norb
        !$omp simd
        band_loop2: do m=1,Norb
           do j=1,3
              do k=1,3
                 if(abs(eig(m,i)-eig(l,i))<1.0d-9)then
                    tmp=vk(k,m,m,i)*vk(j,m,m,i)*ffermi(m,i)*(1.0d0-ffermi(m,i))/(temp*cmplx(w,idelta))
                    L11(k,j)=L11(k,j)+tmp
                    L12(k,j)=L12(k,j)+tmp*(eig(m,i)-mu)
                    L22(k,j)=L22(k,j)+tmp*(eig(m,i)-mu)*(eig(m,i)-mu)
                 else if(abs(ffermi(l,i)-ffermi(m,i))>eps)then
                    tmp=vk(k,m,l,i)*vk(j,l,m,i)*(ffermi(l,i)-ffermi(m,i))/((eig(m,i)-eig(l,i))&
                         *cmplx(w+eig(m,i)-eig(l,i),idelta))
                    L11(k,j)=L11(k,j)+tmp
                    L12(k,j)=L12(k,j)+tmp*(eig(l,i)-mu)
                    L22(k,j)=L22(k,j)+tmp*(eig(m,i)-mu)*(eig(l,i)-mu)
                 end if
              end do
           end do
        end do band_loop2
        !$omp end simd
     end do band_loop1
  end do k_loop
  !$omp end do
  !$omp workshare
  L11(:,:)=ii*L11(:,:)/Nk
  L12(:,:)=ii*L12(:,:)/Nk
  L22(:,:)=ii*L22(:,:)/Nk
  !$omp end workshare
  !$omp end parallel
end subroutine calc_lij

subroutine calc_kn(K0,K1,K2,eig,veloc,kweight,tau,temp,mu,Nk,Norb) bind(C)
  !
  ! calc_Kn
  ! Kn_ij=sum_k(v_ki*v_kj*(e_k-mu)^n*(-df(e_k)/de))
  !
  implicit none
  integer(8),intent(in):: Nk,Norb
  real(8),intent(in):: temp,mu
  real(8),intent(in),dimension(Norb,Nk):: eig,tau
  real(8),intent(in),dimension(Nk):: kweight
  real(8),intent(in),dimension(3,Norb,Nk):: veloc
  real(8),intent(out),dimension(3,3):: K0,K1,K2

  real(8),dimension(Norb,Nk):: dfermi
  integer(8) i,j,l,m
  real(8) tmp

  !$omp parallel
  !$omp do private(j)
  get_dfermi: do i=1,Nk
     do j=1,Norb
        dfermi(j,i)=0.25d0*(1.0d0-tanh(0.5d0*(eig(j,i)-mu)/temp)**2)/temp
     end do
  end do get_dfermi
  !$omp end do
  !$omp workshare
  K0(:,:)=0.0d0
  K1(:,:)=0.0d0
  K2(:,:)=0.0d0
  !$omp end workshare
  
  !$omp do private(j,l,m,tmp) reduction(+:K0,K1,K2)
  get_Kn: do i=1,Nk
     band_loop: do j=1,Norb
        axis1: do l=1,3
           !$omp simd
           axis2: do m=1,3
              tmp=veloc(m,j,i)*veloc(l,j,i)*dfermi(j,i)*tau(j,i)*kweight(i)
              K0(m,l)=K0(m,l)+tmp
              K1(m,l)=K1(m,l)+tmp*(eig(j,i)-mu)
              K2(m,l)=K2(m,l)+tmp*(eig(j,i)-mu)*(eig(j,i)-mu)
           end do axis2
           !$omp end simd
        end do axis1
     end do band_loop
  end do get_Kn
  !$omp end do
  !$omp end parallel
end subroutine calc_kn

subroutine calc_tdf(tdf,eig,veloc,kweight,tau,Nw,Nk,Norb) bind(C)
  !
  ! calc tdf function
  ! sum_k(v_ki*v_kj*tau)
  !
  implicit none
  integer(8),intent(in):: Nk,Norb,Nw
  real(8),intent(in),dimension(Norb,Nk):: eig,tau
  real(8),intent(in),dimension(Nk):: kweight
  real(8),intent(in),dimension(3,Norb,Nk):: veloc
  real(8),intent(out),dimension(3,3,Nw):: tdf

  integer(8) i,j,l,m,iw
  real(8) tmp,emax,emin,id,dw
  id=1.0d-3
  emax=maxval(eig)
  emin=minval(eig)
  dw=(emax-emin)/Nw
  !$omp parallel
  !$omp workshare
  tdf(:,:,:)=0.0d0
  !$omp end workshare
  omega_loop: do iw=1,Nw
     axis1: do l=1,3
        axis2: do m=l,3
           !$omp do private(i,j) reduction(+:tdf)
           k_loop: do i=1,Nk
              band_loop: do j=1,Norb
                 tdf(m,l,iw)=tdf(m,l,iw)+veloc(m,j,i)*veloc(l,j,i)*tau(j,i)*kweight(i)/((iw*dw+emin-eig(j,i))**2+id*id)
              end do band_loop
           end do k_loop
           !$omp end do
           tdf(l,m,iw)=tdf(m,l,iw)
        end do axis2
     end do axis1
  end do omega_loop
  !$omp workshare
  tdf(:,:,:)=tdf(:,:,:)*id/Nk
  !$omp end workshare
  !$omp end parallel
end subroutine calc_tdf
