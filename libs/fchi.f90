subroutine get_qshift(qpoint,klist,qshift,Nk) bind(C)
  use,intrinsic:: iso_fortran_env, only:int64,real64
  implicit none
  integer(int64),intent(in):: Nk
  real(real64),intent(in),dimension(3,Nk):: klist
  real(real64),intent(in),dimension(3):: qpoint
  integer(int64),intent(out),dimension(Nk):: qshift

  call set_qshift(qpoint,klist,qshift,Nk)
end subroutine get_qshift

subroutine set_qshift(qpoint,klist,qshift,Nk)
  !shift k to k+q
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nk
  real(real64),intent(in),dimension(3,Nk):: klist
  real(real64),intent(in),dimension(3):: qpoint
  integer(int64),intent(out),dimension(Nk):: qshift

  integer(int32) i,j,k
  real(real64) tmp
  real(real64),dimension(3,Nk):: kqlist

  !$omp parallel
  !$omp workshare
  kqlist(:,:)=0.0d0
  !$omp end workshare
  !$omp do private(j)
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
  !$omp do private(j,tmp)
  kq_loop: do i=1,Nk
     k_loop: do j=1,Nk
        tmp=sum(abs(klist(:,j)-kqlist(:,i)))
        if(tmp<1.0d-9)then !may be Nk=O(10^9) calc is too difficult
           qshift(i)=j
           exit
        end if
     end do k_loop
  end do kq_loop
  !$omp end do
  !$omp end parallel
end subroutine set_qshift

module calc_irr_chi
  use,intrinsic:: iso_fortran_env, only:int32,int64,real64
  implicit none
contains
  function calc_chi(Nk,Norb,Nchi,uni,eig,ffermi,ol,temp,qshift,w,idelta,eps)
    integer(int64),intent(in):: Nk,Norb,Nchi
    integer(int64),intent(in),dimension(Nk):: qshift
    integer(int64),intent(in),dimension(Nchi,2):: ol
    real(real64),intent(in):: temp,eps,idelta,w
    real(real64),intent(in),dimension(Norb,Nk):: eig,ffermi
    complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni

    integer(int32) i,j,k,l,m
    complex(real64) unitmp
    complex(real64),dimension(Nchi,Nchi):: chi,calc_chi

    chi(:,:)=0.0d0
    kloop: do k=1,Nk
       band1_loop: do l=1,Norb
          band2_loop: do m=1,Norb
             chiorb1_loop: do j=1,Nchi
                chiorb2_loop:do i=1,Nchi
                   unitmp=uni(ol(j,1),l,qshift(k))*conjg(uni(ol(i,1),l,qshift(k)))&
                        *uni(ol(i,2),m,k)*conjg(uni(ol(j,2),m,k))
                   if(abs(w)==0.0d0 .and. abs(eig(m,k)-eig(l,qshift(k)))<1.0d-9)then
                      chi(i,j)=chi(i,j)+unitmp*ffermi(m,k)*(1.0d0-ffermi(m,k))/temp
                   else if(abs(ffermi(l,qshift(k))-ffermi(m,k))>eps)then
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

subroutine get_tr_chi(trchis,trchi0,chis_orb,chis,chi0,olist,Nw,Nchi,Norb) bind(C)
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nchi,Nw,Norb
  integer(int64),intent(in),dimension(Nchi,2):: olist
  complex(real64),intent(in),dimension(Nchi,Nchi,Nw):: chis,chi0
  complex(real64),intent(out),dimension(Nw):: trchis,trchi0
  complex(real64),intent(out),dimension(Norb+2,Nw):: chis_orb

  integer(int32) i,j,k
 
  !$omp parallel do private(j,k)
  wloop:do i=1,Nw
     orb_lop1:do j=1,Nchi
        if(olist(j,1)==olist(j,2))then
           orb_loop2:do k=1,Nchi
              if(olist(k,1)==olist(k,2))then
                 trchis(i)=trchis(i)+chis(k,j,i)
                 trchi0(i)=trchi0(i)+chi0(k,j,i)
                 if(olist(j,1)==olist(k,1))then
                    chis_orb(olist(j,1),i)=chis(k,j,i)
                 end if
                 if(olist(k,1)==2 .and. olist(j,1)==3)chis_orb(Norb+1,i)=chis(k,j,i)
                 if(olist(k,1)==2 .and. olist(j,1)==4)chis_orb(Norb+2,i)=chis(k,j,i)
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
  integer(int64),intent(in):: Nk,Norb,Nw,Nchi
  integer(int64),intent(in),dimension(Nk):: qshift
  integer(int64),intent(in),dimension(Nchi,2):: ol
  real(real64),intent(in):: temp,eps,idelta
  real(real64),intent(in),dimension(Norb,Nk):: eig,ffermi
  real(real64),intent(in),dimension(Nw):: wl
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(real64),intent(out),dimension(Nchi,Nchi,Nw):: chi

  integer(int32) i
 
  !$omp parallel do private(i)
  wloop: do i=1,Nw
     chi(:,:,i)=calc_chi(Nk,Norb,Nchi,uni,eig,ffermi,ol,temp,qshift,wl(i),idelta,eps)
  end do wloop
  !$omp end parallel do
end subroutine get_chi_irr

subroutine chiq_map(trchis,trchi,uni,eig,ffermi,klist,Smat,ol,temp,ecut,idelta,eps,Nx,Ny,Nk,Norb,Nchi) bind(C)
  use calc_irr_chi
  implicit none
  integer(int64),intent(in):: Nx,Ny,Nk,Norb,Nchi
  integer(int64),intent(in),dimension(Nchi,2):: ol
  real(real64),intent(in):: ecut,idelta,eps,temp
  real(real64),intent(in),dimension(3,Nk):: klist
  real(real64),intent(in),dimension(Norb,Nk):: eig,ffermi
  real(real64),intent(in),dimension(Nchi,Nchi):: Smat
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(real64),intent(out),dimension(Ny,Nx):: trchis,trchi

  integer(int32) i,j,l,m,n,info
  integer(int64),dimension(Nk):: qshift
  integer(int32),dimension(Nchi):: ipiv
  real(real64),dimension(3):: qpoint
  complex(real64),dimension(Nchi,Nchi):: chi,tmp,tmp2
  complex(real64),dimension(2*Nchi):: work

  !$omp parallel
  !$omp workshare
  trchi(:,:)=0.0d0
  trchis(:,:)=0.0d0
  !$omp end workshare
  !$omp do private(j,l,m,n,chi,tmp,tmp2,qpoint,qshift,ipiv,work,info)
  do i=1,Nx
     do j=1,Ny
        qpoint(1)=dble(i-1)/Nx
        qpoint(2)=dble(j-1)/Ny
        qpoint(3)=0.0d0
        call set_qshift(qpoint,klist,qshift,Nk)
        chi(:,:)=calc_chi(Nk,Norb,Nchi,uni,eig,ffermi,ol,temp,qshift,ecut,idelta,eps)
        tmp(:,:)=0.0d0
        do l=1,Nchi
           do m=1,Nchi
              !$omp simd
              do n=1,Nchi
                 tmp(m,l)=tmp(m,l)-chi(m,n)*Smat(n,l)
              end do
              !$omp end simd
           end do
           tmp(l,l)=tmp(l,l)+1.0d0
        end do
        call zgetrf(Nchi,Nchi,tmp,Nchi,ipiv,info)
        call zgetri(Nchi,tmp,Nchi,ipiv,work,2*Nchi,info)
        tmp2(:,:)=0.0d0
        do l=1,Nchi
           do m=1,Nchi
              !$omp simd
              do n=1,Nchi
                 tmp2(m,l)=tmp2(m,l)+tmp(m,n)*chi(n,l)
              end do
              !$omp end simd
           end do
        end do
        !take chis_llmm
        do l=1,Nchi
           if(ol(l,1)==ol(l,2))then
              do m=1,Nchi
                 if(ol(m,2)==ol(m,2))then
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
  use,intrinsic:: iso_fortran_env, only:int32,int64,real64
  implicit none
  integer(int64),intent(in):: Nchi,Nw
  real(real64),dimension(Nchi,Nchi):: Smat
  complex(real64),intent(in),dimension(Nchi,Nchi,Nw):: chi0
  complex(real64),intent(out),dimension(Nchi,Nchi,Nw):: chis

  integer(int32) i,l,m,n
  integer(int32) info
  integer(int32),dimension(Nchi):: ipiv
  complex(real64),dimension(Nchi):: work
  complex(real64),dimension(Nchi,Nchi):: tmp

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
              chis(m,l,i)=chis(m,l,i)+tmp(m,n)*chi0(n,l,i)
           end do
           !$omp end simd
        end do
     end do
  end do
  !$omp end parallel do
end subroutine get_chis

subroutine get_smat(Smat,ol,Uval,Jval,Nchi,Norb) bind(C)
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nchi,Norb
  integer(int64),intent(in),dimension(Nchi,2):: ol
  real(real64),intent(in):: Uval,Jval
  real(real64),intent(out),dimension(Nchi,Nchi):: Smat
 
  integer(int32) i,j
 
  !$omp parallel
  !$omp workshare
  Smat(:,:)=0.0d0
  !$omp end workshare
  !$omp do private(j)
  do i=1,Nchi
     do j=1,Nchi
        if((ol(i,1)==ol(i,2)).and.(ol(j,1)==ol(j,2)))then
           if(ol(i,1)==ol(j,1))then
              Smat(j,i)=Uval
           else
              Smat(j,i)=Jval
           end if
        else if((ol(i,1)==ol(j,1)).and.(ol(i,2)==ol(j,2)))then
           Smat(j,i)=Uval-2*Jval
        else if((ol(i,1)==ol(j,2)).and.(ol(i,2)==ol(j,1)))then
           Smat(j,i)=Jval
        end if
     end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine get_smat
