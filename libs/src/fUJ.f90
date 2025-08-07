subroutine get_scmat(Smat,Cmat,ol,Uval,Jval,Nchi) bind(C)
  !
  !> make S and C that is onsite interaction verteces without SOC
  !
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

subroutine get_Vmat_soc(Vmat,osl,Uval,Jval,Nchi,Norb)
  !
  !> make onsite interaction vertex V with SOC
  !
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nchi,Norb
  integer(int64),intent(in),dimension(Nchi,4):: osl !1:orb1, 2:orb2, 3:spin pair, 4:site
  real(real64),intent(in):: Uval,Jval
  real(real64),intent(out),dimension(Nchi,Nchi):: Vmat

  integer(int32) i,j
  real(real64) Upval

  Upval=Uval-2*Jval
  !$omp parallel
  !$omp workshare
  Vmat(:,:)=0.0d0
  !$omp end workshare
  !$omp do private(j)
  do i=1,Nchi
     do j=1,Nchi
        if(osl(i,4)==osl(j,4))then !interaction consider only same site
           if((osl(i,3)==osl(j,3)) .and. (osl(i,3)==1 .or. osl(i,3)==4))then !Vuuuu or Vdddd ((C-S)/2)
              if(osl(i,1)==osl(j,2) .and. osl(j,1)==osl(i,2))then
                 continue               !Viiii and Vijji is zero
              else if(osl(i,1)==osl(i,2) .and. osl(j,1)==osl(j,2))then
                 Vmat(j,i)=Upval-Jval   !Viijj=U'-J
              else if(osl(i,1)==osl(j,1) .and. osl(i,2)==osl(j,2))then
                 Vmat(j,i)=Jval-Upval   !Vijij=-Viijj=J-U'
              end if
           else if((osl(i,3)==2 .and. osl(j,3)==3) .or. (osl(i,3)==3 .and. osl(j,3)==2))then !Vuddn or Vduud (-S)
              if(osl(i,1)==osl(i,2) .and. osl(j,1)==osl(j,2))then
                 if(osl(i,1)==osl(j,1))then
                    Vmat(j,i)=-Uval !Viiii=-U
                 else
                    Vmat(j,i)=-Jval !Viijj=-U'
                 end if
              else if(osl(i,1)==osl(j,1) .and. osl(i,2)==osl(j,2))then
                 Vmat(j,i)=-Upval   !Vijij=-J
              else if(osl(i,1)==osl(j,2) .and. osl(j,1)==osl(i,2))then
                 Vmat(j,i)=-Jval    !Vijji=-J
              end if
           else if((osl(i,3)==1 .and. osl(j,3)==4) .or. (osl(i,3)==4 .and. osl(j,3)==1))then !Vuudd or Vdduu ((C+S)/2)
              if(osl(i,1)==osl(i,2) .and. osl(j,1)==osl(j,2))then
                 if(osl(i,1)==osl(j,1))then
                    Vmat(j,i)=Uval  !Viiii=U
                 else
                    Vmat(j,i)=Upval !Viijj=U'
                 end if
              else if(osl(i,1)==osl(j,1) .and. osl(i,2)==osl(j,2))then
                 Vmat(j,i)=Jval     !Vijij=J
              else if(osl(i,1)==osl(j,2) .and. osl(j,1)==osl(i,2))then
                 Vmat(j,i)=Jval     !Vijji=J'
              end if
           end if
        end if
     end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine get_Vmat_soc
