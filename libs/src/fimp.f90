subroutine gen_imp_ham(ham_imp,ham_r,rvec,ham_i,imp_list,rlist,eps,Nimp,Nsite,Nr,Norb) bind(C)
  use,intrinsic:: iso_fortran_env, only: int64,real64
  implicit none
  integer(int64),intent(in):: Nimp,Nsite,Norb,Nr
  integer(int64),intent(in),dimension(Nimp):: imp_list
  real(real64),intent(in):: eps
  real(real64),intent(in),dimension(3,Nsite):: rlist
  real(real64),intent(in),dimension(3,Nr):: rvec
  complex(real64),intent(in),dimension(Norb,Norb,Nr):: ham_r,ham_i
  complex(real64),intent(out),dimension(Norb*Nsite,Norb*Nsite):: ham_imp

  integer(int64) i,j,k,l,m,n,k_onsite,Nrx,Nry,Nrz
  logical sw_imp1,sw_imp2
  real(real64) tmpr(3)

  Nrx=maxval(rlist(1,:))+1
  Nry=maxval(rlist(2,:))+1
  Nrz=maxval(rlist(3,:))+1
  !detect onsite_rvec
  do k=1,Nr
     if(sum(abs(rvec(:,k)))<eps)then
        k_onsite=k
        exit
     end if
  end do

  !$omp parallel do private(j,k,l,m,n,sw_imp1,sw_imp2,tmpr)
  site_loop1:do i=1,Nsite
     sw_imp1=.false.
     do n=1,Nimp
        if(i==imp_list(n)+1)then
           sw_imp1=.true.
           exit
        end if
     end do
     if(sw_imp1)then
        do l=1,Norb
           do m=l,Norb
              ham_imp(m+(i-1)*Norb,l+(i-1)*Norb)=ham_i(m,l,k_onsite)
              ham_imp(l+(i-1)*Norb,m+(i-1)*Norb)=conjg(ham_i(m,l,k_onsite))
           end do
           ham_imp(l+(i-1)*Norb,l+(i-1)*Norb)=ham_imp(l+(i-1)*Norb,l+(i-1)*Norb)
        end do
     else
        do l=1,Norb
           do m=l,Norb
              ham_imp(m+(i-1)*Norb,l+(i-1)*Norb)=ham_r(m,l,k_onsite)
              ham_imp(l+(i-1)*Norb,m+(i-1)*Norb)=conjg(ham_r(m,l,k_onsite))
           end do
           ham_imp(l+(i-1)*Norb,l+(i-1)*Norb)=ham_imp(l+(i-1)*Norb,l+(i-1)*Norb)
        end do
     end if
     site_loop2: do j=i+1,Nsite
        sw_imp2=.false.
        do n=1,Nimp
           if(j==imp_list(n))then
              sw_imp2=.true.
              exit
           end if
        end do
        tmpr(:)=rlist(:,j)-rlist(:,i) !get_distance i to j
        if(tmpr(1)>Nrx/2)tmpr(1)=tmpr(1)-Nrx
        if(tmpr(2)>Nry/2)tmpr(2)=tmpr(2)-Nry
        if(tmpr(3)>Nrz/2)tmpr(3)=tmpr(3)-Nrz
        if(tmpr(1)< -Nrx/2)tmpr(1)=tmpr(1)+Nrx
        if(tmpr(2)< -Nry/2)tmpr(2)=tmpr(2)+Nry
        if(tmpr(3)< -Nrz/2)tmpr(3)=tmpr(3)+Nrz
        !print'(3F4.1,2I2)',tmpr(:),i,j
        hop_loop:do k=1,Nr
           if(sum(abs(tmpr(:)-rvec(:,k)))<eps)then
              do l=1,Norb
                 !$omp simd
                 do m=1,Norb
                    ham_imp(m+(j-1)*Norb,l+(i-1)*Norb)=ham_r(m,l,k)
                    ham_imp(l+(i-1)*Norb,m+(j-1)*Norb)=conjg(ham_r(m,l,k))
                 end do
                 !$omp end simd
              end do
              exit
           end if
        end do hop_loop
     end do site_loop2
  end do site_loop1
  !$omp end parallel do
end subroutine gen_imp_ham

subroutine get_dft_imp_ham(ham_k,ham_imp,klist,rlist,Nk,Nsite,Norb) bind(C)
  use constant  
  implicit none
  integer(int64),intent(in):: Nk,Nsite,Norb
  real(real64),intent(in),dimension(3,Nsite):: rlist
  real(real64),intent(in),dimension(3,Nk):: klist
  complex(real64),intent(in),dimension(Norb*Nsite,Norb*Nsite):: ham_imp
  complex(real64),intent(out),dimension(Norb*Nk,Norb*Nk):: ham_k

  integer i,j,k,l,m,n
  real(real64) phase

  !$omp parallel do private(j,k,l,m,n,phase)
  k_loop1: do i=1,Nk
     k_loop2: do j=i,Nk
        site_loop1: do k=1,Nsite
           site_loop2: do l=1,Nsite
              phase=2*pi*(sum(klist(:,j)*rlist(:,l))-sum(klist(:,i)*rlist(:,k)))
              orb_loop1: do m=1,Norb
                 orb_loop2: do n=1,Norb
                    ham_k(n+Nk*(j-1),m+Nk*(i-1))=ham_k(n+Nk*(j-1),m+Nk*(i-1))&
                         +ham_imp(n+Nsite*(l-1),m+Nsite*(k-1))*cmplx(cos(phase),-sin(phase))
                 end do orb_loop2
              end do orb_loop1
           end do site_loop2
        end do site_loop1
        orbk_loop: do m=1,Norb !orb_loop for hermite setting for Hamiltonian
           orbk_loop2: do n=1,Norb
              ham_k(m+Nk*(i-1),n+Nk*(j-1))=conjg(ham_k(n+Nk*(j-1),m+Nk*(i-1)))
           end do orbk_loop2
        end do orbk_loop
     end do k_loop2
  end do k_loop1
  !$omp end parallel do
end subroutine get_dft_imp_ham

subroutine get_spectrum_spagehtti(spa,uni,eigs,klist,rlist,wlist,Nw,Nk,Nsite,Norb,mu,eta) bind(C)
  use constant
  implicit none
  integer(int64),intent(in):: Nw,Nk,Nsite,Norb
  real(real64),intent(in):: eta,mu
  real(real64),intent(in),dimension(Nw):: wlist
  real(real64),intent(in),dimension(3,Nsite):: rlist
  real(real64),intent(in),dimension(3,Nk):: klist
  real(real64),intent(in),dimension(Norb*Nsite):: eigs
  complex(real64),intent(in),dimension(Norb*Nsite,Norb*Nsite):: uni
  complex(real64),intent(out),dimension(Nw,Nk):: spa

  integer i,j,k,l,m,n
  real(real64) phase

  !$omp parallel
  !$omp workshare
  spa(:,:)=0.0d0
  !$omp end workshare
  !$omp do private(i,j,k,l,n,m,phase)
  k_loop: do i=1,Nk
     site_loop1: do j=1,Nsite
        site_loop2: do k=1,Nsite
           phase=2*pi*sum(klist(:,i)*(rlist(:,j)-rlist(:,k)))
           orb_loop: do l=1,Norb
              eig_loop: do n=1,Nsite*Norb
                 w_loop:do m=1,Nw
                    spa(m,i)=spa(m,i)+conjg(uni(l+(j-1)*Nsite,n))*uni(l+(k-1)*Nsite,n)&
                         *cmplx(cos(phase),-sin(phase))/cmplx(wlist(m)-eigs(n)+mu,eta)
                 end do w_loop
              end do eig_loop
           end do orb_loop
        end do site_loop2
     end do site_loop1
  end do k_loop
  !$omp end do
  !$omp end parallel
end subroutine get_spectrum_spagehtti
