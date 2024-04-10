subroutine get_iqshift(qpoint,klist,qshift,Nk) bind(C)
  use,intrinsic:: iso_fortran_env, only:int64,real64
  implicit none
  integer(int64),intent(in):: Nk
  real(real64),intent(in),dimension(3,Nk):: klist
  real(real64),intent(in),dimension(3):: qpoint
  integer(int64),intent(out),dimension(Nk):: qshift
  
  call set_iqshift(qpoint,klist,qshift,Nk)
end subroutine get_iqshift
  
subroutine set_iqshift(qpoint,klist,qshift,Nk)
  !shift k to -k+q
  use,intrinsic:: iso_fortran_env, only:int64,real64
  implicit none
  integer(int64),intent(in):: Nk
  real(real64),intent(in),dimension(3,Nk):: klist
  real(real64),intent(in),dimension(3):: qpoint
  integer(int64),intent(out),dimension(Nk):: qshift
  
  integer(int64) i,j,k,ck
  real(real64) tmp
  real(real64),dimension(3,Nk):: kqlist
  
  !$omp parallel
  !$omp workshare
  kqlist(:,:)=0.0d0
  !$omp end workshare
  !$omp do private(j)
  kloop:do i=1,Nk
     kqlist(:,i)=1.0d0-klist(:,i)+qpoint(:)
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
  !$omp do private(j,tmp)
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
end subroutine set_iqshift

module calc_irr_phi
  use,intrinsic:: iso_fortran_env, only:int32,int64,real64
  implicit none
contains
  function calc_phi(Nk,Norb,Nchi,uni,eig,ffermi,ol,mu,temp,qshift,w,idelta,eps)
    integer(int64),intent(in):: Nk,Norb,Nchi
    integer(int64),intent(in),dimension(Nk):: qshift
    integer(int64),intent(in),dimension(Nchi,2):: ol
    real(real64),intent(in):: mu,temp,eps,idelta,w
    real(real64),intent(in),dimension(Norb,Nk):: eig,ffermi
    complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  
    integer(int64) i,j,k,l,m
    complex(real64) unitmp
    complex(real64),dimension(Nchi,Nchi):: phi,calc_phi
  
    phi(:,:)=0.0d0
    kloop: do k=1,Nk
       band1_loop: do l=1,Norb
          band2_loop: do m=1,Norb
             chiorb1_loop: do j=1,Nchi
                chiorb2_loop:do i=1,Nchi
                   unitmp=uni(ol(j,1),l,qshift(k))*conjg(uni(ol(i,1),l,qshift(k)))&
                        *uni(ol(i,2),m,k)*conjg(uni(ol(j,2),m,k))
                   phi(i,j)=phi(i,j)-unitmp*(1.0d0-ffermi(l,qshift(k))-ffermi(m,k))&
                        /cmplx(w-eig(m,k)-eig(l,qshift(k)+2.0d0*mu),idelta)
                end do chiorb2_loop
             end do chiorb1_loop
          end do band2_loop
       end do band1_loop
    end do kloop
    calc_phi=phi(:,:)/Nk
  end function calc_phi
end module calc_irr_phi

subroutine get_phi_irr(phi,uni,eig,ffermi,qshift,ol,wl,Nchi,Norb,Nk,Nw,idelta,eps,mu,temp) bind(C)
  use calc_irr_phi
  implicit none
  integer(int64),intent(in):: Nk,Norb,Nw,Nchi
  integer(int64),intent(in),dimension(Nk):: qshift
  integer(int64),intent(in),dimension(Nchi,2):: ol
  real(real64),intent(in):: temp,mu,eps,idelta
  real(real64),intent(in),dimension(Norb,Nk):: eig,ffermi
  real(real64),intent(in),dimension(Nw):: wl
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(real64),intent(out),dimension(Nchi,Nchi,Nw):: phi

  integer(int64) i
    
  !$omp parallel do private(i)
  wloop: do i=1,Nw
     phi(:,:,i)=calc_phi(Nk,Norb,Nchi,uni,eig,ffermi,ol,mu,temp,qshift,wl(i),idelta,eps)
  end do wloop
  !$omp end parallel do
end subroutine get_phi_irr

subroutine get_tr_phi(trphi,phi_orb,phi,olist,Nw,Nchi,Norb) bind(C)
  use,intrinsic:: iso_fortran_env, only:int32,int64,real64
  implicit none
  integer(int64),intent(in):: Nchi,Nw,Norb
  integer(int64),intent(in),dimension(Nchi,2):: olist
  complex(real64),intent(in),dimension(Nchi,Nchi,Nw):: phi
  complex(real64),intent(out),dimension(Nw):: trphi
  complex(real64),intent(out),dimension(Norb+2,Nw):: phi_orb
  
  integer(int64) i,j,k
  
  !$omp parallel do private(j,k)
  wloop:do i=1,Nw
     orb_lop1:do j=1,Nchi
        trphi(i)=trphi(i)+phi(j,j,i)
        if(olist(j,1)==olist(j,2))then
           orb_loop2:do k=1,Nchi
              if(olist(k,1)==olist(k,2))then
                 if(olist(j,1)==olist(k,1))then
                    phi_orb(olist(j,1),i)=phi(k,j,i)
                 end if
                 if(olist(k,1)==2 .and. olist(j,1)==3)phi_orb(Norb+1,i)=phi(k,j,i)
                 if(olist(k,1)==2 .and. olist(j,1)==4)phi_orb(Norb+2,i)=phi(k,j,i)
              end if
           end do orb_loop2
        end if
     end do orb_lop1
  end do wloop
  !$omp end parallel do
end subroutine get_tr_phi

subroutine phiq_map(trphi,uni,eig,ffermi,klist,ol,mu,temp,ecut,idelta,eps,Nx,Ny,Nk,Norb,Nchi) bind(C)
  use calc_irr_phi
  implicit none
  integer(int64),intent(in):: Nx,Ny,Nk,Norb,Nchi
  integer(int64),intent(in),dimension(Nchi,2):: ol
  real(real64),intent(in):: ecut,idelta,eps,temp,mu
  real(real64),intent(in),dimension(3,Nk):: klist
  real(real64),intent(in),dimension(Norb,Nk):: eig,ffermi
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(real64),intent(out),dimension(Ny,Nx):: trphi
  
  integer(int32) info
  integer(int64) i,j,l,m,n
  integer(int64),dimension(Nk):: qshift
  integer(int32),dimension(Nchi):: ipiv
  real(real64),dimension(3):: qpoint
  complex(real64),dimension(Nchi,Nchi):: phi
  complex(real64),dimension(2*Nchi):: work
  
  !$omp parallel
  !$omp workshare
  trphi(:,:)=0.0d0
  !$omp end workshare
  !$omp do private(i,j,l,m,n,phi,qpoint,qshift,ipiv,work,info)
  do i=1,Nx
     do j=1,Ny
        qpoint(1)=dble(i-1)/Nx
        qpoint(2)=dble(j-1)/Ny
        qpoint(3)=0.0d0
        call set_iqshift(qpoint,klist,qshift,Nk)
        phi(:,:)=calc_phi(Nk,Norb,Nchi,uni,eig,ffermi,ol,mu,temp,qshift,ecut,idelta,eps)
        do l=1,Nchi
           trphi(j,i)=trphi(j,i)+phi(l,l)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine phiq_map
