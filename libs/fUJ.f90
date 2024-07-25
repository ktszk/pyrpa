subroutine get_scmat(Smat,Cmat,ol,Uval,Jval,Nchi) bind(C)
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nchi
  integer(int64),intent(in),dimension(Nchi,2):: ol
  real(real64),intent(in):: Uval,Jval
  real(real64),intent(out),dimension(Nchi,Nchi):: Smat,Cmat
   
  integer(int32) i,j
  real(real64) Upval

  Upval=Uval-2*Jval
  !$omp parallel
  !$omp workshare
  Smat(:,:)=0.0d0
  Cmat(:,:)=0.0d0
  !$omp end workshare
  !$omp do private(j)
  do i=1,Nchi
     do j=1,Nchi
        if((ol(i,1)==ol(i,2)).and.(ol(j,1)==ol(j,2)))then
           if(ol(i,1)==ol(j,1))then
              Smat(j,i)=Uval !iiii>U
              Cmat(j,i)=Uval !iiii>U
           else
              Smat(j,i)=Jval         !iijj>J
              Cmat(j,i)=2*Upval-Jval !iijj>2U'-J
           end if
        else if((ol(i,1)==ol(j,1)).and.(ol(i,2)==ol(j,2)))then
           Smat(j,i)=Upval         !ijij>U'
           Cmat(j,i)=-Upval+2*Jval !ijij>2J-U'
        else if((ol(i,1)==ol(j,2)).and.(ol(i,2)==ol(j,1)))then
           Smat(j,i)=Jval    !ijji>J'
           Cmat(j,i)=Jval    !ijji>J'
        end if
     end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine get_scmat
  
