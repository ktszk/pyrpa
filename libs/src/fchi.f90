subroutine get_qshift(qpoint,klist,qshift,Nk) bind(C,name="get_qshift_")
  !> shift k to k+q (hash-based O(Nk) mapping)
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nk
  real(real64),intent(in),dimension(3,Nk):: klist
  real(real64),intent(in),dimension(3):: qpoint
  integer(int64),intent(out),dimension(Nk):: qshift

  integer(int32) :: i,j
  integer(int64) :: jj
  integer(int64) :: scale_int, m, h, step
  integer(int64) :: i1, i2, i3, j1, j2, j3
  integer(int64), allocatable :: ht_key(:), ht_val(:)
  integer(int64) :: key, key_t
  real(real64) :: tmpf(3)

  ! parameters for integerized grid
  scale_int = 1048576_int64  ! 2^20 resolution
  m = max(3_int64, 2_int64*Nk + 3_int64)

  allocate(ht_key(m))
  allocate(ht_val(m))
  ht_key(:) = -1_int64
  ht_val(:) = -1_int64

  ! Build hash table from klist (single-threaded)
  do jj = 1, int(Nk)
     i1 = int(klist(1,jj)*dble(scale_int)+0.5d0, int64)
     i2 = int(klist(2,jj)*dble(scale_int)+0.5d0, int64)
     i3 = int(klist(3,jj)*dble(scale_int)+0.5d0, int64)
     i1 = mod(i1, scale_int)
     i2 = mod(i2, scale_int)
     i3 = mod(i3, scale_int)
     key = i1 + i2*scale_int + i3*scale_int*scale_int
     h = mod(abs(key), m) + 1_int64
     step = 1_int64
     do while(ht_key(h) /= -1_int64 .and. ht_key(h) /= key)
        h = h + step
        if(h > m) h = h - m
     end do
     ht_key(h) = key
     ht_val(h) = jj
  end do

  ! Lookup for shifted k = k + qpoint (periodic mod 1)
  !$omp parallel do private(i1,i2,i3,key_t,h,step,tmpf,j1,j2,j3)
  do i = 1, Nk
     tmpf(:) = klist(:,i) + qpoint(:)
     do j = 1, 3
        if(tmpf(j) >= 1.0d0) then
           tmpf(j) = tmpf(j) - 1.0d0
        else if(tmpf(j) < 0.0d0) then
           tmpf(j) = tmpf(j) + 1.0d0
        end if
     end do
     j1 = int(tmpf(1)*dble(scale_int)+0.5d0, int64)
     j2 = int(tmpf(2)*dble(scale_int)+0.5d0, int64)
     j3 = int(tmpf(3)*dble(scale_int)+0.5d0, int64)
     j1 = mod(j1, scale_int)
     j2 = mod(j2, scale_int)
     j3 = mod(j3, scale_int)
     key_t = j1 + j2*scale_int + j3*scale_int*scale_int
     h = mod(abs(key_t), m) + 1_int64
     step = 1_int64
     qshift(i) = 1
     do while(ht_key(h) /= -1_int64)
        if(ht_key(h) == key_t) then
           qshift(i) = ht_val(h)
           exit
        end if
        h = h + step
        if(h > m) h = h - m
     end do
  end do
  !$omp end parallel do

  deallocate(ht_key)
  deallocate(ht_val)
end subroutine get_qshift

module calc_irr_chi
  use,intrinsic:: iso_fortran_env, only:int32,int64,real64
  implicit none
contains
  function calc_chi(Nk,Norb,Nchi,uni,eig,ffermi,ol,temp,qshift,w,idelta,eps)
    !> This function obtain irreducible susceptibility chi_0
    !!@param        Nk: The number of k-points
    !!@param      Norb: The number of orbitals
    !!@param      Nchi: The footnote of chi
    !!@param       uni: unitary matrix

      !!@param       eig: energies of bands
    !!@param    ffermi: fermi distribute function
    !!@param        ol: the list of the properties of orbitals at footnote of chi
    !!@param      temp: temperature
    !!@param    qshift: q-shifted klist
    !!@param         w: frequency
    !!@param    idelta: dumping factor
    !!@param       eps: threshold of calculation value
    !!@return calc_chi: irreducible susceptibility matrix
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
  !> This function take trace of chi
  !!@param   trchis,out: trace value of chi_s
  !!@param   trchi0,out: trace value of chi_0
  !!@param chis_orb,out: The array of chi_s values at orbital diagonal elements
  !!@param      chis,in: spin susceptibility
  !!@param      chi0,in: irreducible susceptibility
  !!@param     olist,in: the list of  properties of chis,chi0 footnote
  !!@param        Nw,in: The number of frequency mesh
  !!@param      Nchi,in: The number of footnote of chis,chi0
  !!@param      Norb,in: The number of orbitals
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nchi,Nw,Norb
  integer(int64),intent(in),dimension(Nchi,2):: olist
  complex(real64),intent(in),dimension(Nchi,Nchi,Nw):: chis,chi0
  complex(real64),intent(out),dimension(Nw):: trchis,trchi0
  complex(real64),intent(out),dimension(Norb+2,Nw):: chis_orb

  integer(int32) i,j,k

  trchis(:)    = (0.0d0, 0.0d0)
  trchi0(:)    = (0.0d0, 0.0d0)
  chis_orb(:,:)= (0.0d0, 0.0d0)

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
  !> This function obtains irreducible susceptibility at q-point
  !!@param   chi,out: irreducible susceptibility
  !!@param    uni,in: unitary matrix
  !!@param    eig,in: energies of bands
  !!@param ffermi,in: fermi distribute function
  !!@param qshift,in: q-shifted klist
  !!@param     ol,in: the list of the properties of orbitals at footnote of chi
  !!@param     wl,in: the list of frequency
  !!@param   Nchi,in: The footnote of chi
  !!@param   Norb,in: The number of orbitals
  !!@param     Nk,in: The number of k-points
  !!@param     Nw,in: The number of frequency mesh
  !!@param idelta,in: dumping factor
  !!@param    eps,in: threshold of calculation value
  !!@param   temp,in: temperature
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
  !> This function obtains irreducible susceptibility at w=ecut
  !!@param trchis,out: trace value of chi_s
  !!@param  trchi,out: trace value of chi_0
  !!@param     uni,in: unitary matrix
  !!@param     eig,in: energies of bands
  !!@param  ffermi,in: fermi distribute function
  !!@param   klist,in: klist
  !!@param    Smat,in: S-matrix
  !!@param      ol,in: the list of the properties of orbitals at footnote of chi
  !!@param    temp,in: temperature
  !!@param    ecut,in: energy value of cut energy plane
  !!@param  idelta,in: dumping factor
  !!@param     eps,in: threshold of calculation value
  !!@param      Nx,in: The number of kx mesh
  !!@param      Ny,in: The number of ky mesh
  !!@param      Nk,in: The number of k-points
  !!@param    Norb,in: The number of orbitals
  !!@param    Nchi,in: The footnote of chi
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

  integer(int32) i,j,l,m,info
  integer(int64),dimension(Nk):: qshift
  integer(int32),dimension(Nchi):: ipiv
  real(real64),dimension(3):: qpoint
  complex(real64),dimension(Nchi,Nchi):: chi,tmp,tmp2,Smat_c
  complex(real64),dimension(2*Nchi):: work

  Smat_c = Smat

  !$omp parallel
  !$omp workshare
  trchi(:,:)=0.0d0
  trchis(:,:)=0.0d0
  !$omp end workshare
  !$omp do private(j,l,chi,tmp,tmp2,qpoint,qshift,ipiv,work,info)
  do i=1,Nx
     do j=1,Ny
        qpoint(1)=dble(i-1)/Nx
        qpoint(2)=dble(j-1)/Ny
        qpoint(3)=0.0d0
        call get_qshift(qpoint,klist,qshift,Nk)
        chi(:,:)=calc_chi(Nk,Norb,Nchi,uni,eig,ffermi,ol,temp,qshift,ecut,idelta,eps)
        ! tmp = I - chi * Smat
        call zgemm('N','N',Nchi,Nchi,Nchi,(-1.0d0,0.0d0),chi,Nchi,Smat_c,Nchi,(0.0d0,0.0d0),tmp,Nchi)
        do l=1,Nchi
           tmp(l,l)=tmp(l,l)+(1.0d0,0.0d0)
        end do
        call zgetrf(Nchi,Nchi,tmp,Nchi,ipiv,info)
        if(info/=0)then; print*,'zgetrf failed: info=',info; stop; end if
        call zgetri(Nchi,tmp,Nchi,ipiv,work,2*Nchi,info)
        if(info/=0)then; print*,'zgetri failed: info=',info; stop; end if
        ! tmp2 = tmp * chi
        call zgemm('N','N',Nchi,Nchi,Nchi,(1.0d0,0.0d0),tmp,Nchi,chi,Nchi,(0.0d0,0.0d0),tmp2,Nchi)
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
  !> This function obtain spin susceptibility
  !!@param chis,out: spin susceptibility
  !!@param  chi0,in: irreducible susceptibility
  !!@param  Smat,in: S-matrix
  !!@param  Nchi,in: The footnote of chi
  !!@param    Nw,in: The number of w mesh
  use,intrinsic:: iso_fortran_env, only:int32,int64,real64
  implicit none
  integer(int64),intent(in):: Nchi,Nw
  real(real64),dimension(Nchi,Nchi):: Smat
  complex(real64),intent(in),dimension(Nchi,Nchi,Nw):: chi0
  complex(real64),intent(out),dimension(Nchi,Nchi,Nw):: chis

  integer(int32) i,l,info
  integer(int32),dimension(Nchi):: ipiv
  complex(real64),dimension(2*Nchi):: work
  complex(real64),dimension(Nchi,Nchi):: tmp,Smat_c

  Smat_c = Smat

  !$omp parallel do private(tmp,l,work,ipiv,info)
  do i=1,Nw
     ! tmp = I - chi0(:,:,i) * Smat
     call zgemm('N','N',Nchi,Nchi,Nchi,(-1.0d0,0.0d0),chi0(:,:,i),Nchi,Smat_c,Nchi,(0.0d0,0.0d0),tmp,Nchi)
     do l=1,Nchi
        tmp(l,l)=tmp(l,l)+(1.0d0,0.0d0)
     end do
     call zgetrf(Nchi,Nchi,tmp,Nchi,ipiv,info)
     if(info/=0)then; print*,'zgetrf failed: info=',info; stop; end if
     call zgetri(Nchi,tmp,Nchi,ipiv,work,2*Nchi,info)
     if(info/=0)then; print*,'zgetri failed: info=',info; stop; end if
     ! chis(:,:,i) = tmp * chi0(:,:,i)
     call zgemm('N','N',Nchi,Nchi,Nchi,(1.0d0,0.0d0),tmp,Nchi,chi0(:,:,i),Nchi,(0.0d0,0.0d0),chis(:,:,i),Nchi)
  end do
  !$omp end parallel do
end subroutine get_chis
