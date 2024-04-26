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

  integer(int64) i,j,k,l,m,n,k_onsite
  logical sw_imp1,sw_imp2
  real(real64) tmpr(3)

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
        if(i==imp_list(n))then
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
        end do
     else
        do l=1,Norb
           do m=l,Norb
              ham_imp(m+(i-1)*Norb,l+(i-1)*Norb)=ham_r(m,l,k_onsite)
              ham_imp(l+(i-1)*Norb,m+(i-1)*Norb)=conjg(ham_r(m,l,k_onsite))
           end do
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
  do i=1,Nk
     do j=i,Nk
        do k=1,Nsite
           do l=1,Nsite
              phase=2*pi*(sum(klist(:,j)*rlist(:,l))-sum(klist(:,i)*rlist(:,k)))
              do m=1,Norb
                 do n=1,Norb
                    ham_k(n+Nk*(j-1),m+Nk*(i-1))=ham_k(n+Nk*(j-1),m+Nk*(i-1))&
                         +ham_imp(n+Nsite*(l-1),m+Nsite*(k-1))*cmplx(cos(phase),-sin(phase))
                 end do
              end do
           end do
        end do
        do m=1,Norb
           do n=1,Norb
              ham_k(m+Nk*(i-1),n+Nk*(j-1))=conjg(ham_k(n+Nk*(j-1),m+Nk*(i-1)))
           end do
        end do
     end do
  end do
  !$omp end parallel do
end subroutine get_dft_imp_ham
