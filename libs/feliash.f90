subroutine get_V_delta_nsoc_flex(chi,Smat,Cmat,Nk,Nw,Nchi,sw_pair)
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nk,Nw,Nchi
  logical(1),intent(in):: sw_pair
  real(real64),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  complex(real64),intent(inout),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(int32) i,j,l,m,n,info
  integer(int32),dimension(Nchi):: ipiv
  complex(real64),dimension(Nchi,Nchi):: cmat1,cmat2,cmat3,cmat4,cmat5
  complex(real64),dimension(2*Nchi):: work

  do j=1,Nw
     do i=1,Nk
        !$omp parallel
        !$omp workshare
        cmat1(:,:)=0.0d0
        cmat2(:,:)=0.0d0
        !$omp end workshare
        !$omp do private(m,n)
        do l=1,Nchi
           do m=1,Nchi
              do n=1,Nchi
                 cmat1(m,l)=cmat1(m,l)-chi(i,j,m,n)*Smat(n,l) !-chi0S
                 cmat2(m,l)=cmat2(m,l)+chi(i,j,m,n)*Cmat(n,l) !chi0C
              end do
           end do
        end do
        !$omp end do
        !$omp workshare
        cmat3(:,:)=-cmat1(:,:) !chi0S
        cmat4(:,:)=cmat2(:,:)  !chi0C
        !$omp end workshare
        !$omp do 
        do l=1,Nchi
           cmat1(l,l)=cmat1(l,l)+1.0d0 !I-chi0S
           cmat2(l,l)=cmat2(l,l)+1.0d0 !I+chi0C
        end do
        !$omp end do
        !$omp end parallel
        call zgetrf(Nchi,Nchi,cmat1,Nchi,ipiv,info)
        call zgetri(Nchi,cmat1,Nchi,ipiv,work,2*Nchi,info)
        call zgetrf(Nchi,Nchi,cmat2,Nchi,ipiv,info)
        call zgetri(Nchi,cmat2,Nchi,ipiv,work,2*Nchi,info)
        !$omp parallel
        !$omp workshare
        cmat5(:,:)=0.0d0
        !$omp end workshare
        !$omp do private(m,n)
        do l=1,Nchi
           do m=1,Nchi
              do n=1,Nchi
                 cmat5(m,l)=cmat5(m,l)+cmat1(m,n)*cmat3(n,l) !(1-chi0S)^-1chi0S
              end do
           end do
        end do
        !$omp end do
        !$omp workshare
        cmat1(:,:)=cmat5(:,:)
        cmat5(:,:)=0.0d0
        !$omp end workshare
        !$omp do private(m,n)
        do l=1,Nchi
           do m=1,Nchi
              do n=1,Nchi
                 cmat5(m,l)=cmat5(m,l)+cmat2(m,n)*cmat4(n,l) !(1-chi0S)^-1chi0S
              end do
           end do
        end do
        !$omp end do
        !$omp workshare
        cmat2(:,:)=cmat5(:,:)
        !$omp end workshare
        if(sw_pair)then !singlet
           !$omp workshare
           cmat4(:,:)=0.5d0*(Smat(:,:)+Cmat(:,:))
           !$omp end workshare
           !$omp do private(m,n)
           do l=1,Nchi
              do m=1,Nchi
                 do n=1,Nchi
                    cmat4(m,l)=cmat4(m,l)+1.5d0*Smat(m,n)*cmat1(n,l)-0.5d0*Cmat(m,n)*cmat2(n,l)
                 end do
              end do
           end do
           !$omp end do
        else !triplet
           !$omp workshare
           cmat4(:,:)=0.5d0*(Smat(:,:)-Cmat(:,:))
           !$omp end workshare
           !$omp do private(m,n)
           do l=1,Nchi
              do m=1,Nchi
                 do n=1,Nchi
                    cmat4(m,l)=cmat4(m,l)-0.5d0*(Smat(m,n)*cmat1(n,l)+Cmat(m,n)*cmat2(n,l))
                 end do
              end do
           end do
           !$omp end do
        end if
        !$omp workshare
        chi(i,j,:,:)=cmat4(:,:)
        !$omp end workshare
        !$omp end parallel
     end do
  end do
end subroutine get_V_delta_nsoc_flex
