subroutine get_scmat(Smat,Cmat,ol,site,Uval,Jval,Nchi) bind(C)
  !> make S and C that is onsite interaction verteces without SOC
  !!@param Smat,out: S-matrix (vertex matrix for spin), size(Nchi,Nchi)
  !!@param Cmat,out: C-matrix (vertex matrix for charge), size(Nchi,Nchi)
  !!@param    ol,in: list of orbitals, size(Nchi,2)
  !!@param  Uval,in: U value
  !!@param  Jval,in: J value
  !!@param  Nchi,in: The number of orbitals for chi
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nchi
  integer(int64),intent(in),dimension(Nchi,2):: ol
  integer(int64),intent(in),dimension(Nchi):: site
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
        if(site(i)==site(j))then
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
        end if
     end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine get_scmat

subroutine get_scmat_orb(Smat,Cmat,ol,site,Umat,Jmat,Nchi,Norb) bind(C)
  !> make S and C that is orbital dependent onsite interaction verteces without SOC
  !!@param Smat,out: S-matrix (vertex matrix for spin), size(Nchi,Nchi)
  !!@param Cmat,out: C-matrix (vertex matrix for charge), size(Nchi,Nchi)
  !!@param    ol,in: list of orbitals, size(Nchi,2)
  !!@param  Umat,in: orbital dependent U
  !!@param  Jmat,in: orbital dependent J
  !!@param  Nchi,in: The number of footnote of chi
  !!@param  Norb,in: The number of orbitals
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nchi,Norb
  integer(int64),intent(in),dimension(Nchi,2):: ol
  integer(int64),intent(in),dimension(Nchi):: site
  real(real64),intent(in),dimension(Norb,Norb):: Umat,Jmat
  real(real64),intent(out),dimension(Nchi,Nchi):: Smat,Cmat

  integer(int32) i,j

  !$omp parallel
  !$omp workshare
  Smat(:,:)=0.0d0
  Cmat(:,:)=0.0d0
  !$omp end workshare
  !$omp do private(j)
  do i=1,Nchi
     do j=1,Nchi
        if(site(i)==site(j))then
           if((ol(i,1)==ol(i,2)).and.(ol(j,1)==ol(j,2)))then
              if(ol(i,1)==ol(j,1))then
                 Smat(j,i)=Umat(ol(i,1),ol(i,1)) !iiii>U
                 Cmat(j,i)=Umat(ol(i,1),ol(i,1)) !iiii>U
              else
                 Smat(j,i)=Jmat(ol(i,1),ol(j,1)) !iijj>J
                 Cmat(j,i)=2*Umat(ol(i,1),ol(j,1))-Jmat(ol(i,1),ol(j,1)) !iijj>2U'-J
              end if
           else if((ol(i,1)==ol(j,1)).and.(ol(i,2)==ol(j,2)))then
              Smat(j,i)=Umat(ol(i,1),ol(i,2)) !ijij>U'
              Cmat(j,i)=-Umat(ol(i,1),ol(i,2))+2*Jmat(ol(i,1),ol(i,2)) !ijij>2J-U'
           else if((ol(i,1)==ol(j,2)).and.(ol(i,2)==ol(j,1)))then
              Smat(j,i)=Jmat(ol(i,1),ol(i,2))  !ijji>J'
              Cmat(j,i)=Jmat(ol(i,1),ol(i,2))  !ijji>J'
           end if
        end if
     end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine get_scmat_orb

subroutine get_Vmat_soc(Vmat,ol,sl,site,invs,Uval,Jval,Nchi,Norb) bind(C)
  !> make onsite interaction vertex V with SOC
  !!@param Vmat,out: V-matrix (vertex matrix for soc calculation), size(Nchi,Nchi)
  !!@param    ol,in: list of orbitals, size(Nchi,2)
  !!@param    sl,in: list of spin, size(Norb)
  !!@param  site,in: list of atomic site, size(Nchi)
  !!@param  Uval,in: U value
  !!@param  Jval,in: J value
  !!@param  Nchi,in: The number of orbitals for chi
  !!@param  Norb,in: The number of orbitals
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nchi,Norb
  integer(int64),intent(in),dimension(Nchi,2):: ol !1:orb1, 2:orb2
  integer(int64),intent(in),dimension(Nchi):: site
  integer(int64),intent(in),dimension(Norb):: sl,invs
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
        if(site(i)==site(j))then !interaction consider only same site
           if(sl(ol(i,1))==sl(ol(i,2)) .and. sl(ol(j,1))==sl(ol(j,2)))then
              if(sl(ol(i,1))==sl(ol(j,1)))then !Vuuuu or Vdddd ((C-S)/2)
                 if(ol(i,1)==ol(j,2) .and. ol(j,1)==ol(i,2))then
                    continue             !Viiii and Vijji is zero (intra,pair hoppings)
                 else if(ol(i,1)==ol(i,2) .and. ol(j,1)==ol(j,2))then
                    Vmat(j,i)=Upval-Jval !Viijj=U'-J (inter-hund)
                 else if(ol(i,1)==ol(j,1) .and. ol(i,2)==ol(j,2))then
                    Vmat(j,i)=Jval-Upval !Vijij=-Viijj=J-U' (inter-hund)
                 end if
              else !Vuudd or Vdduu ((C+S)/2)
                 if(ol(i,1)==ol(i,2) .and. ol(j,1)==ol(j,2))then
                    if(ol(i,1)==invs(ol(j,1)))then
                       Vmat(j,i)=Uval  !Viiii=U (intra)
                    else
                       Vmat(j,i)=Upval !Viijj=U' (inter)
                    end if
                 else if(ol(i,1)==invs(ol(j,1)) .and. ol(i,2)==invs(ol(j,2)))then
                    Vmat(j,i)=Jval     !Vijij=J (hund)
                 else if(ol(i,1)==invs(ol(j,2)) .and. ol(j,1)==invs(ol(i,2)))then
                    Vmat(j,i)=Jval     !Vijji=J' (pair hoppings)
                 end if
              end if
           else if(sl(ol(i,1))==sl(ol(j,1)) .and. sl(ol(i,2))==sl(ol(j,2)))then !Vudud or Vdudu (-S)
              if(ol(i,1)==invs(ol(i,2)) .and. ol(j,1)==invs(ol(j,2)))then
                 if(ol(i,1)==ol(j,1))then
                    Vmat(j,i)=-Uval !Viiii=-U (intra)
                 else
                    Vmat(j,i)=-Jval !Viijj=-J (hund)
                 end if
              else if(ol(i,1)==ol(j,1) .and. ol(i,2)==ol(j,2))then
                 Vmat(j,i)=-Upval   !Vijij=-U' (inter)
              else if(ol(i,1)==invs(ol(j,2)) .and. ol(j,1)==invs(ol(i,2)))then
                 Vmat(j,i)=-Jval    !Vijji=-J' (pair hoppings)
              end if
           end if
        end if
     end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine get_Vmat_soc

subroutine get_Vmat_soc_orb(Vmat,ol,sl,site,invs,Umat,Jmat,Nchi,Norb) bind(C)
  !> make orbital dependent onsite interaction vertex V with SOC
  !!@param Vmat,out: V-matrix (vertex matrix for soc calculation), size(Nchi,Nchi)
  !!@param    ol,in: list of orbitals, size(Nchi,2)
  !!@param    sl,in: list of spin, size(Norb)
  !!@param  site,in: list of atomic site, size(Nchi)
  !!@param  Umat,in: orbital dependent U
  !!@param  Jmat,in: orbital dependent J
  !!@param  Nchi,in: The number of orbitals for chi
  !!@param  Norb,in: The number of orbitals
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nchi,Norb
  integer(int64),intent(in),dimension(Nchi,2):: ol !1:orb1, 2:orb2
  integer(int64),intent(in),dimension(Nchi):: site
  integer(int64),intent(in),dimension(Norb):: sl,invs
  real(real64),intent(in),dimension(Norb,Norb):: Umat,Jmat
  real(real64),intent(out),dimension(Nchi,Nchi):: Vmat

  integer(int32) i,j

  !$omp parallel
  !$omp workshare
  Vmat(:,:)=0.0d0
  !$omp end workshare
  !$omp do private(j)
  do i=1,Nchi
     do j=1,Nchi
        if(site(i)==site(j))then !interaction consider only same site
           if(sl(ol(i,1))==sl(ol(i,2)) .and. sl(ol(j,1))==sl(ol(j,2)))then
              if(sl(ol(i,1))==sl(ol(j,1)))then !Vuuuu or Vdddd ((C-S)/2)
                 if(ol(i,1)==ol(j,2) .and. ol(j,1)==ol(i,2))then
                    continue                   !Viiii and Vijji is zero (intra & pair hoppings)
                 else if(ol(i,1)==ol(i,2) .and. ol(j,1)==ol(j,2))then
                    Vmat(j,i)=Umat(ol(i,1),ol(j,1))-Jmat(ol(i,1),ol(j,1)) !Viijj=U'-J (inter-hund)
                 else if(ol(i,1)==ol(j,1) .and. ol(i,2)==ol(j,2))then
                    Vmat(j,i)=Jmat(ol(i,1),ol(i,2))-Umat(ol(i,1),ol(i,2)) !Vijij=-Viijj=J-U' (inter-hund)
                 end if
              else !Vuudd or Vdduu ((C+S)/2)
                 if(ol(i,1)==ol(i,2) .and. ol(j,1)==ol(j,2))then
                    if(ol(i,1)==invs(ol(j,1)))then
                       Vmat(j,i)=Umat(ol(i,1),ol(i,1)) !Viiii=U (intra)
                    else
                       Vmat(j,i)=Umat(ol(i,1),ol(j,1)) !Viijj=U' (inter)
                    end if
                 else if(ol(i,1)==invs(ol(j,1)) .and. ol(i,2)==invs(ol(j,2)))then
                    Vmat(j,i)=Jmat(ol(i,1),ol(i,2))    !Vijij=J (hund)
                 else if(ol(i,1)==invs(ol(j,2)) .and. ol(j,1)==invs(ol(i,2)))then
                    Vmat(j,i)=Jmat(ol(i,1),ol(i,2))    !Vijji=J' (pair hoppings)
                 end if
              end if
           else if(sl(ol(i,1))==sl(ol(j,1)) .and. sl(ol(i,2))==sl(ol(j,2)))then !Vudud or Vdudu (-S)
              if(ol(i,1)==invs(ol(i,2)) .and. ol(j,1)==invs(ol(j,2)))then
                 if(ol(i,1)==ol(j,1))then
                    Vmat(j,i)=-Umat(ol(i,1),ol(i,1)) !Viiii=-U (intra)
                 else
                    Vmat(j,i)=-Jmat(ol(i,1),ol(j,1)) !Viijj=-U' (inter)
                 end if
              else if(ol(i,1)==ol(j,1) .and. ol(i,2)==ol(j,2))then
                 Vmat(j,i)=-Umat(ol(i,1),ol(i,2))    !Vijij=-J (hund)
              else if(ol(i,1)==invs(ol(j,2)) .and. ol(j,1)==invs(ol(i,2)))then
                 Vmat(j,i)=-Jmat(ol(i,1),ol(i,2))    !Vijji=-J (pair hoppings)
              end if
           end if
        end if
     end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine get_Vmat_soc_orb
