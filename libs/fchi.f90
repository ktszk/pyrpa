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
