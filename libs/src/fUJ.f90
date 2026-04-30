subroutine get_scmat(Smat,Cmat,ol,site,Uval,Jval,Nchi) bind(C)
  !> make S and C that is onsite interaction verteces without SOC
  !!@param Smat,out: S-matrix (vertex matrix for spin), size(Nchi,Nchi)
  !!@param Cmat,out: C-matrix (vertex matrix for charge), size(Nchi,Nchi)
  !!@param    ol,in: list of orbitals, size(Nchi,2)
  !!@param  Uval,in: U value
  !!@param  Jval,in: J value
  !!@param  Nchi,in: The number of orbitals for chi
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nchi
  integer(c_int64_t),intent(in),dimension(Nchi,2):: ol
  integer(c_int64_t),intent(in),dimension(Nchi):: site
  real(c_double),intent(in):: Uval,Jval
  real(c_double),intent(out),dimension(Nchi,Nchi):: Smat,Cmat

  integer(c_int32_t) i,j
  real(c_double) Upval

  ! Kanamori Hamiltonian: U'=U-2J (standard relation for t2g orbitals)
  ! Smat = spin vertex (enters spin susceptibility chi_s)
  ! Cmat = charge vertex (enters charge susceptibility chi_c)
  ! Matrix elements follow from the multi-orbital Hubbard model:
  !   S[iijj] = J,  C[iijj] = 2U'-J   (i/=j, density-density)
  !   S[ijij] = U', C[ijij] = 2J-U'   (Hund's exchange)
  !   S[ijji] = J,  C[ijji] = J        (pair hopping)
  Upval=Uval-2*Jval
  !$omp parallel
  !$omp workshare
  Smat(:,:)=0.0d0
  Cmat(:,:)=0.0d0
  !$omp end workshare
  !$omp do private(j)
  do i=1,Nchi
     do j=1,Nchi
        if(site(i)==site(j))then   ! on-site interactions only
           if((ol(i,1)==ol(i,2)).and.(ol(j,1)==ol(j,2)))then
              if(ol(i,1)==ol(j,1))then
                 Smat(j,i)=Uval !iiii: intra-orbital U
                 Cmat(j,i)=Uval !iiii: intra-orbital U
              else
                 Smat(j,i)=Jval         !iijj: Hund J
                 Cmat(j,i)=2*Upval-Jval !iijj: 2U'-J
              end if
           else if((ol(i,1)==ol(j,1)).and.(ol(i,2)==ol(j,2)))then
              Smat(j,i)=Upval         !ijij: inter-orbital U'
              Cmat(j,i)=-Upval+2*Jval !ijij: 2J-U'
           else if((ol(i,1)==ol(j,2)).and.(ol(i,2)==ol(j,1)))then
              Smat(j,i)=Jval    !ijji: pair hopping J'=J
              Cmat(j,i)=Jval    !ijji: pair hopping J'=J
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
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nchi,Norb
  integer(c_int64_t),intent(in),dimension(Nchi,2):: ol
  integer(c_int64_t),intent(in),dimension(Nchi):: site
  real(c_double),intent(in),dimension(Norb,Norb):: Umat,Jmat
  real(c_double),intent(out),dimension(Nchi,Nchi):: Smat,Cmat

  integer(c_int32_t) i,j

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
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nchi,Norb
  integer(c_int64_t),intent(in),dimension(Nchi,2):: ol !1:orb1, 2:orb2
  integer(c_int64_t),intent(in),dimension(Nchi):: site
  integer(c_int64_t),intent(in),dimension(Norb):: sl,invs
  real(c_double),intent(in):: Uval,Jval
  real(c_double),intent(out),dimension(Nchi,Nchi):: Vmat

  integer(c_int32_t) i,j
  real(c_double) Upval

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
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nchi,Norb
  integer(c_int64_t),intent(in),dimension(Nchi,2):: ol !1:orb1, 2:orb2
  integer(c_int64_t),intent(in),dimension(Nchi):: site
  integer(c_int64_t),intent(in),dimension(Norb):: sl,invs
  real(c_double),intent(in),dimension(Norb,Norb):: Umat,Jmat
  real(c_double),intent(out),dimension(Nchi,Nchi):: Vmat

  integer(c_int32_t) i,j

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

subroutine get_V_delta_nsoc_flex(chi,Smat,Cmat,Nk,Nw,Nchi,sw_pair)
  !> This function obtains pairing interaction V_delta without soc
  !!@param  chi,inout: irreducible susceptibility and pairing interaction
  !!@param    Smat,in: S-matrix
  !!@param    Cmat,in: C-matrix
  !!@param      Nk,in: Number of k-points
  !!@param      Nw,in: Number of Matsubara frequencies
  !!@param    Nchi,in: Number of footnote of chi
  !!@param sw_pair,in: switch of singlet or triplet paring interacton
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,Nchi
  logical(1),intent(in):: sw_pair
  real(c_double),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  complex(c_double),intent(inout),dimension(Nk,Nw,Nchi,Nchi):: chi

   integer(c_int32_t) i,j,l,info
  integer(c_int32_t),dimension(Nchi):: ipiv
   complex(c_double),dimension(Nchi,Nchi):: cmat1,cmat2,cmat3,cmat4,cmat5,Smat_c,Cmat_c,V0_c

   Smat_c(:,:)=cmplx(Smat(:,:),0.0d0,kind=c_double)
   Cmat_c(:,:)=cmplx(Cmat(:,:),0.0d0,kind=c_double)
   ! Bare (static) pairing vertex from Kanamori model:
   !   singlet: V0 = (S+C)/2 = Vud  (spin-singlet ↑↓ channel)
   !   triplet: V0 = (S-C)/2 = Vuu  (spin-triplet ↑↑ channel)
   if(sw_pair)then
      V0_c(:,:)=cmplx(0.5d0*(Smat(:,:)+Cmat(:,:)),0.0d0,kind=c_double) !bare Vud=(C+S)/2
   else
      V0_c(:,:)=cmplx(0.5d0*(Smat(:,:)-Cmat(:,:)),0.0d0,kind=c_double) !bare Vuu=(S-C)/2
   end if

  !$omp parallel do collapse(2) private(i,cmat1,cmat2,cmat3,cmat4,cmat5,ipiv,info,l)
  wloop:do j=1,Nw
     qloop:do i=1,Nk
        call zgemm('N','N',Nchi,Nchi,Nchi,(-1.0d0,0.0d0),chi(i,j,:,:),Nchi,Smat_c,Nchi,(0.0d0,0.0d0),cmat1,Nchi) !-chi0S
        call zgemm('N','N',Nchi,Nchi,Nchi,(1.0d0,0.0d0),chi(i,j,:,:),Nchi,Cmat_c,Nchi,(0.0d0,0.0d0),cmat2,Nchi)  !chi0C
        cmat3(:,:)=-cmat1(:,:) !chi0S (RHS for first solve)
        cmat4(:,:)=cmat2(:,:)  !chi0C (RHS for second solve)
        do l=1,Nchi
           cmat1(l,l)=cmat1(l,l)+1.0d0 !I-chi0S
           cmat2(l,l)=cmat2(l,l)+1.0d0 !I+chi0C
        end do
        call zgesv(Nchi,Nchi,cmat1,Nchi,ipiv,cmat3,Nchi,info) !(I-chi0S)X=chi0S -> cmat3=chiS
        if(info/=0)then; print*,'zgesv failed: info=',info; stop; end if
        call zgesv(Nchi,Nchi,cmat2,Nchi,ipiv,cmat4,Nchi,info) !(I+chi0C)X=chi0C -> cmat4=chiC
        if(info/=0)then; print*,'zgesv failed: info=',info; stop; end if
        cmat1(:,:)=cmat3(:,:) !cmat1=chiS
        cmat2(:,:)=cmat4(:,:) !cmat2=chiC
        cmat4(:,:)=V0_c(:,:)   !start from bare static vertex
        ! FLEX pairing vertex (RPA-dressed):
        !   singlet: V_Δ = (S+C)/2 + 3/2·S·χ_s - 1/2·C·χ_c
        !   triplet: V_Δ = (S-C)/2 - 1/2·S·χ_s - 1/2·C·χ_c
        if(sw_pair)then !singlet
           call zgemm('N','N',Nchi,Nchi,Nchi,(1.5d0,0.0d0),Smat_c,Nchi,cmat1,Nchi,(1.0d0,0.0d0),cmat4,Nchi)  !+3/2·S·χ_s
           call zgemm('N','N',Nchi,Nchi,Nchi,(-0.5d0,0.0d0),Cmat_c,Nchi,cmat2,Nchi,(1.0d0,0.0d0),cmat4,Nchi) !-1/2·C·χ_c
        else !triplet
           call zgemm('N','N',Nchi,Nchi,Nchi,(-0.5d0,0.0d0),Smat_c,Nchi,cmat1,Nchi,(1.0d0,0.0d0),cmat4,Nchi) !-1/2·S·χ_s
           call zgemm('N','N',Nchi,Nchi,Nchi,(-0.5d0,0.0d0),Cmat_c,Nchi,cmat2,Nchi,(1.0d0,0.0d0),cmat4,Nchi) !-1/2·C·χ_c
        end if
        chi(i,j,:,:)=cmat4(:,:)
     end do qloop
  end do wloop
  !$omp end parallel do
end subroutine get_V_delta_nsoc_flex

subroutine get_vsigma_flex_nosoc(chi,Smat,Cmat,Nk,Nw,Nchi) bind(C,name='get_vsigma_flex_nosoc_')
  !> This function obtains interaction V_sigma without soc
  !!@param  chi,inout: irreducible susceptibility and pairing interaction
  !!@param    Smat,in: S-matrix
  !!@param    Cmat,in: C-matrix
  !!@param      Nk,in: Number of k-points
  !!@param      Nw,in: Number of Matsubara frequencies
  !!@param    Nchi,in: Number of footnote of chi
   use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,Nchi
  real(c_double),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  complex(c_double),intent(inout),dimension(Nk,Nw,Nchi,Nchi):: chi

   integer(c_int32_t) i,j,l,info
  integer(c_int32_t),dimension(Nchi):: ipiv
   complex(c_double),dimension(Nchi,Nchi):: cmat1,cmat2,cmat3,cmat4,cmat5,Smat_c,Cmat_c,SC_c

   Smat_c(:,:)=cmplx(Smat(:,:),0.0d0,kind=c_double)
   Cmat_c(:,:)=cmplx(Cmat(:,:),0.0d0,kind=c_double)
   SC_c(:,:)=Smat_c(:,:)+Cmat_c(:,:)

  !$omp parallel do collapse(2) private(i,cmat1,cmat2,cmat3,cmat4,cmat5,ipiv,info,l)
  do j=1,Nw
     do i=1,Nk
        call zgemm('N','N',Nchi,Nchi,Nchi,(-1.0d0,0.0d0),chi(i,j,:,:),Nchi,Smat_c,Nchi,(0.0d0,0.0d0),cmat1,Nchi) !-chi0S
        call zgemm('N','N',Nchi,Nchi,Nchi,(1.0d0,0.0d0),chi(i,j,:,:),Nchi,Cmat_c,Nchi,(0.0d0,0.0d0),cmat2,Nchi)  !chi0C
        cmat3(:,:)=-cmat1(:,:) !chi0S (RHS for first solve)
        cmat4(:,:)=cmat2(:,:)  !chi0C (RHS for second solve)
        cmat5(:,:)=cmat3(:,:)+cmat4(:,:) !chi0S+chi0C (save before zgesv overwrites)
        do l=1,Nchi
           cmat1(l,l)=cmat1(l,l)+1.0d0 !I-chi0S
           cmat2(l,l)=cmat2(l,l)+1.0d0 !I+chi0C
        end do
        call zgesv(Nchi,Nchi,cmat1,Nchi,ipiv,cmat3,Nchi,info) !(I-chi0S)X=chi0S -> cmat3=chiS
        if(info/=0)then; print*,'zgesv failed: info=',info; stop; end if
        call zgesv(Nchi,Nchi,cmat2,Nchi,ipiv,cmat4,Nchi,info) !(I+chi0C)X=chi0C -> cmat4=chiC
        if(info/=0)then; print*,'zgesv failed: info=',info; stop; end if
        ! FLEX self-energy kernel: V_σ = 3/2·S·χ_s + 1/2·C·χ_c - 1/4·(S+C)·(χ_s+χ_c) + static (3/2·S - 1/2·C)
        cmat1(:,:)=1.5d0*Smat(:,:)-0.5d0*Cmat(:,:)   !static (HF-like) bare vertex
        call zgemm('N','N',Nchi,Nchi,Nchi,(1.5d0,0.0d0),Smat_c,Nchi,cmat3,Nchi,(1.0d0,0.0d0),cmat1,Nchi)  !+3/2·S·χ_s
        call zgemm('N','N',Nchi,Nchi,Nchi,(0.5d0,0.0d0),Cmat_c,Nchi,cmat4,Nchi,(1.0d0,0.0d0),cmat1,Nchi)  !+1/2·C·χ_c
        call zgemm('N','N',Nchi,Nchi,Nchi,(-0.25d0,0.0d0),SC_c,Nchi,cmat5,Nchi,(1.0d0,0.0d0),cmat1,Nchi) !-1/4·(S+C)·(χ_s+χ_c), subtract double count
        chi(i,j,:,:)=cmat1(:,:)
     end do
  end do
  !$omp end parallel do
end subroutine get_vsigma_flex_nosoc

subroutine get_V_soc_flex(chi,Vmat,sgnsig2,Nk,Nw,Nchi)
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,Nchi
  real(c_double),intent(in),dimension(Nchi,Nchi)::sgnsig2
  real(c_double),intent(in),dimension(Nchi,Nchi):: Vmat
  complex(c_double),intent(inout),dimension(Nk,Nw,Nchi,Nchi):: chi
  
   integer(c_int32_t) i,j,l,info
  integer(c_int32_t),dimension(Nchi):: ipiv
   complex(c_double),dimension(Nchi,Nchi):: cmat1,cmat2,cmat3,cmatv,neg_cmatv

   cmatv(:,:)=cmplx(Vmat(:,:),0.0d0,kind=c_double)
   neg_cmatv(:,:)=-cmatv(:,:)

  ! SOC-RPA pairing vertex: V_Δ = -V + V·χ,  χ = (I+χ₀·V)^{-1}·χ₀·V
  !$omp parallel do collapse(2) private(i,cmat1,cmat2,cmat3,ipiv,info,l)
  wloop:do j=1,Nw
     qloop:do i=1,Nk
        call zgemm('N','N',Nchi,Nchi,Nchi,(1.0d0,0.0d0),chi(i,j,:,:),Nchi,cmatv,Nchi,(0.0d0,0.0d0),cmat1,Nchi) !chi0V
        cmat2(:,:)=cmat1(:,:) !chi0V (RHS for zgesv)
        do l=1,Nchi
           cmat1(l,l)=cmat1(l,l)+1.0d0 !I+chi0V
        end do
        call zgesv(Nchi,Nchi,cmat1,Nchi,ipiv,cmat2,Nchi,info) !(I+chi0V)X=chi0V -> cmat2=chi (RPA)
        if(info/=0)then; print*,'zgesv failed: info=',info; stop; end if
        cmat3(:,:)=neg_cmatv(:,:)   !start from -V
        call zgemm('N','N',Nchi,Nchi,Nchi,(1.0d0,0.0d0),cmatv,Nchi,cmat2,Nchi,(1.0d0,0.0d0),cmat3,Nchi) !-V + V*chi
        chi(i,j,:,:)=cmat3(:,:)
     end do qloop
  end do wloop
  !$omp end parallel do
end subroutine get_V_soc_flex