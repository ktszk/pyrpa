subroutine get_qshift(qpoint,klist,qshift,Nk) bind(C,name="get_qshift_")
  !> shift k to k+q (hash-based O(Nk) mapping)
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nk
  real(c_double),intent(in),dimension(3,Nk):: klist
  real(c_double),intent(in),dimension(3):: qpoint
  integer(c_int64_t),intent(out),dimension(Nk):: qshift

  integer(c_int32_t) :: i,j
  integer(c_int64_t) :: jj
  integer(c_int64_t) :: scale_int, m, h, step
  integer(c_int64_t) :: i1, i2, i3, j1, j2, j3
  integer(c_int64_t), allocatable :: ht_key(:), ht_val(:)
  integer(c_int64_t) :: key, key_t
  integer(c_int64_t) :: nmiss
  logical :: found
  real(c_double) :: tmpf(3)

  ! Hash-table approach for O(Nk) k+q lookup (vs O(Nk^2) brute force)
  ! k-coordinates are quantized to integer keys at 2^20 resolution
  scale_int = 1048576_c_int64_t  ! 2^20 resolution (~10^-6 fractional precision)
  m = max(3_c_int64_t, 2_c_int64_t*Nk + 3_c_int64_t)  ! table size > 2*Nk to keep load factor < 0.5

  allocate(ht_key(m))
  allocate(ht_val(m))
  ht_key(:) = -1_c_int64_t
  ht_val(:) = -1_c_int64_t

  ! Build hash table from klist (single-threaded)
  do jj = 1, int(Nk)
     i1 = int(klist(1,jj)*dble(scale_int)+0.5d0, c_int64_t)
     i2 = int(klist(2,jj)*dble(scale_int)+0.5d0, c_int64_t)
     i3 = int(klist(3,jj)*dble(scale_int)+0.5d0, c_int64_t)
     i1 = mod(i1, scale_int)
     i2 = mod(i2, scale_int)
     i3 = mod(i3, scale_int)
     key = i1 + i2*scale_int + i3*scale_int*scale_int
     h = mod(abs(key), m) + 1_c_int64_t
     step = 1_c_int64_t
     do while(ht_key(h) /= -1_c_int64_t .and. ht_key(h) /= key)
        h = h + step
        if(h > m) h = h - m
     end do
     ht_key(h) = key
     ht_val(h) = jj
  end do

  ! Lookup for shifted k = k + qpoint (periodic mod 1)
  nmiss = 0
  !$omp parallel do private(i1,i2,i3,key_t,h,step,tmpf,j1,j2,j3,found) reduction(+:nmiss)
  do i = 1, Nk
     tmpf(:) = klist(:,i) + qpoint(:)
     do j = 1, 3
        if(tmpf(j) >= 1.0d0) then
           tmpf(j) = tmpf(j) - 1.0d0
        else if(tmpf(j) < 0.0d0) then
           tmpf(j) = tmpf(j) + 1.0d0
        end if
     end do
     j1 = int(tmpf(1)*dble(scale_int)+0.5d0, c_int64_t)
     j2 = int(tmpf(2)*dble(scale_int)+0.5d0, c_int64_t)
     j3 = int(tmpf(3)*dble(scale_int)+0.5d0, c_int64_t)
     j1 = mod(j1, scale_int)
     j2 = mod(j2, scale_int)
     j3 = mod(j3, scale_int)
     key_t = j1 + j2*scale_int + j3*scale_int*scale_int
     h = mod(abs(key_t), m) + 1_c_int64_t
     step = 1_c_int64_t
     qshift(i) = 1
     found = .false.
     do while(ht_key(h) /= -1_c_int64_t)
        if(ht_key(h) == key_t) then
           qshift(i) = ht_val(h)
           found = .true.
           exit
        end if
        h = h + step
        if(h > m) h = h - m
     end do
     if(.not. found) nmiss = nmiss + 1
  end do
  !$omp end parallel do
  if(nmiss > 0)then
     print'(A,I0,A,I0,A)','WARNING(get_qshift): ',nmiss,' of ',Nk,' k+q points not found on the k-mesh.'
     print'(A)','  klist must be a periodic [0,1) mesh without endpoints (e.g. gen_klist with sw_pp=False).'
     print'(A)','  Missing points fall back to index 1; the resulting chi/phi is unreliable.'
  end if

  deallocate(ht_key)
  deallocate(ht_val)
end subroutine get_qshift

module calc_irr_chi
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
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
    integer(c_int64_t),intent(in):: Nk,Norb,Nchi
    integer(c_int64_t),intent(in),dimension(Nk):: qshift
    integer(c_int64_t),intent(in),dimension(Nchi,2):: ol
    real(c_double),intent(in):: temp,eps,idelta,w
    real(c_double),intent(in),dimension(Norb,Nk):: eig,ffermi
    complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni

    integer(c_int32_t) i,j,k,l,m,nchi32
    real(c_double) temp_safe,w_eps
    complex(c_double) weight
    complex(c_double),dimension(Nchi):: A_vec,B_vec
    complex(c_double),dimension(Nchi,Nchi):: chi,calc_chi

    temp_safe=max(temp,1.0d-12)
    w_eps=1.0d-12
    nchi32=int(Nchi,c_int32_t)
    chi(:,:)=0.0d0
    kloop: do k=1,Nk
       band1_loop: do l=1,Norb
          band2_loop: do m=1,Norb
             ! compute scalar weight once per (k,l,m)
             if(abs(w)<w_eps .and. abs(eig(m,k)-eig(l,qshift(k)))<1.0d-9)then
                weight=cmplx(ffermi(m,k)*(1.0d0-ffermi(m,k))/temp_safe,0.0d0,kind=c_double)
             else if(abs(ffermi(l,qshift(k))-ffermi(m,k))>eps)then
                weight=(ffermi(l,qshift(k))-ffermi(m,k))&
                     /cmplx(w+eig(m,k)-eig(l,qshift(k)),idelta,kind=c_double)
             else
                cycle band2_loop
             end if
             ! A_vec(j) = uni(ol(j,1),l,qshift(k)) * conjg(uni(ol(j,2),m,k))
             do j=1,Nchi
                A_vec(j)=uni(ol(j,1),l,qshift(k))*conjg(uni(ol(j,2),m,k))
             end do
             ! B_vec(i) = conjg(uni(ol(i,1),l,qshift(k))) * uni(ol(i,2),m,k)
             do i=1,Nchi
                B_vec(i)=conjg(uni(ol(i,1),l,qshift(k)))*uni(ol(i,2),m,k)
             end do
             ! chi(i,j) += weight * B_vec(i) * A_vec(j)
             call zgeru(nchi32,nchi32,weight,B_vec,1,A_vec,1,chi,nchi32)
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
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nchi,Nw,Norb
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  complex(c_double),intent(in),dimension(Nchi,Nchi,Nw):: chis,chi0
  complex(c_double),intent(out),dimension(Nw):: trchis,trchi0
  complex(c_double),intent(out),dimension(Norb+2,Nw):: chis_orb

  integer(c_int32_t) i,j,k

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
  integer(c_int64_t),intent(in):: Nk,Norb,Nw,Nchi
  integer(c_int64_t),intent(in),dimension(Nk):: qshift
  integer(c_int64_t),intent(in),dimension(Nchi,2):: ol
  real(c_double),intent(in):: temp,eps,idelta
  real(c_double),intent(in),dimension(Norb,Nk):: eig,ffermi
  real(c_double),intent(in),dimension(Nw):: wl
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(c_double),intent(out),dimension(Nchi,Nchi,Nw):: chi

  integer(c_int32_t) i
 
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
  integer(c_int64_t),intent(in):: Nx,Ny,Nk,Norb,Nchi
  integer(c_int64_t),intent(in),dimension(Nchi,2):: ol
  real(c_double),intent(in):: ecut,idelta,eps,temp
  real(c_double),intent(in),dimension(3,Nk):: klist
  real(c_double),intent(in),dimension(Norb,Nk):: eig,ffermi
  real(c_double),intent(in),dimension(Nchi,Nchi):: Smat
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(c_double),intent(out),dimension(Ny,Nx):: trchis,trchi

  integer(c_int32_t) i,j,l,m,info
  integer(c_int64_t),dimension(Nk):: qshift
  integer(c_int32_t),dimension(Nchi):: ipiv
  real(c_double),dimension(3):: qpoint
  complex(c_double),dimension(Nchi,Nchi):: chi,tmp,tmp2,Smat_c

  Smat_c = Smat

  !$omp parallel
  !$omp workshare
  trchi(:,:)=0.0d0
  trchis(:,:)=0.0d0
  !$omp end workshare
  !$omp do private(j,l,chi,tmp,tmp2,qpoint,qshift,ipiv,info)
  do i=1,Nx
     do j=1,Ny
        qpoint(1)=dble(i-1)/Nx
        qpoint(2)=dble(j-1)/Ny
        qpoint(3)=0.0d0
        call get_qshift(qpoint,klist,qshift,Nk)
        chi(:,:)=calc_chi(Nk,Norb,Nchi,uni,eig,ffermi,ol,temp,qshift,ecut,idelta,eps)
        ! RPA: chi_s = (I - chi0*S)^{-1} * chi0
        tmp2(:,:)=chi(:,:)
        call zgemm('N','N',Nchi,Nchi,Nchi,(-1.0d0,0.0d0),chi,Nchi,Smat_c,Nchi,(0.0d0,0.0d0),tmp,Nchi)
        do l=1,Nchi
           tmp(l,l)=tmp(l,l)+(1.0d0,0.0d0)
        end do
        call zgesv(Nchi,Nchi,tmp,Nchi,ipiv,tmp2,Nchi,info)
        if(info/=0)then; print*,'zgesv failed: info=',info; stop; end if
        !take chis_llmm
        do l=1,Nchi
           if(ol(l,1)==ol(l,2))then
              do m=1,Nchi
                 if(ol(m,1)==ol(m,2))then
                    ! Sum chi_{ll,mm} components (same convention as get_tr_chi)
                    trchis(j,i)=trchis(j,i)+tmp2(m,l)
                    trchi(j,i)=trchi(j,i)+chi(m,l)
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
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nchi,Nw
  real(c_double),dimension(Nchi,Nchi):: Smat
  complex(c_double),intent(in),dimension(Nchi,Nchi,Nw):: chi0
  complex(c_double),intent(out),dimension(Nchi,Nchi,Nw):: chis

  integer(c_int32_t) i,l,info
  integer(c_int32_t),dimension(Nchi):: ipiv
  complex(c_double),dimension(Nchi,Nchi):: tmp,Smat_c

  Smat_c = Smat

  !$omp parallel do private(tmp,l,ipiv,info)
  do i=1,Nw
     chis(:,:,i)=chi0(:,:,i)
     call zgemm('N','N',Nchi,Nchi,Nchi,(-1.0d0,0.0d0),chi0(:,:,i),Nchi,Smat_c,Nchi,(0.0d0,0.0d0),tmp,Nchi)
     do l=1,Nchi
        tmp(l,l)=tmp(l,l)+(1.0d0,0.0d0)
     end do
     call zgesv(Nchi,Nchi,tmp,Nchi,ipiv,chis(:,:,i),Nchi,info)
     if(info/=0)then; print*,'zgesv failed: info=',info; stop; end if
  end do
  !$omp end parallel do
end subroutine get_chis

subroutine ckchi_impl(chi,Smat,Cmat,kmap,invk,Nk,Nkall,Nchi,Nw,maxchi0s_out)
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nkall,Nchi,Nw
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
  real(c_double),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  real(c_double),intent(out):: maxchi0s_out
  complex(c_double),intent(in),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(c_int32_t) i,info,chisk,chick,chiskall,chickall
  real(c_double) maxchi0s,maxchi0c,maxchi0s2,maxchi0c2
  real(c_double),dimension(2*Nchi):: rwork
  complex(c_double),dimension(Nchi*Nchi*4+1):: work
  complex(c_double),dimension(Nchi):: eigs,eigc
  complex(c_double),dimension(Nchi,Nchi):: chi0s,chi0c,tmp1,tmp2,Smat_c,Cmat_c,chi_tmp

  Smat_c = cmplx(Smat, 0.0d0, kind=c_double)
  Cmat_c = cmplx(Cmat, 0.0d0, kind=c_double)
  maxchi0s2=-1.0d5
  maxchi0c2=-1.0d5
  chisk=1; chick=1; chiskall=1; chickall=1  ! defaults in case no eigenvalue exceeds the initial max
  do i=1,Nk
     chi_tmp = chi(i,1,:,:)
     call zgemm('N','N',Nchi,Nchi,Nchi,(1.0d0,0.0d0),chi_tmp,Nchi,Smat_c,Nchi,(0.0d0,0.0d0),chi0s,Nchi)
     call zgemm('N','N',Nchi,Nchi,Nchi,(-1.0d0,0.0d0),chi_tmp,Nchi,Cmat_c,Nchi,(0.0d0,0.0d0),chi0c,Nchi)
     call zgeev('N','N',Nchi,chi0s,Nchi,eigs,tmp1,Nchi,tmp2,Nchi,work,Nchi*Nchi*4+1,rwork,info)
     call zgeev('N','N',Nchi,chi0c,Nchi,eigc,tmp1,Nchi,tmp2,Nchi,work,Nchi*Nchi*4+1,rwork,info)
     maxchi0s=maxval(dble(eigs))
     maxchi0c=maxval(dble(eigc))
     if(maxchi0s>maxchi0s2)then
        chisk=i
        maxchi0s2=maxchi0s
     end if
     if(maxchi0c>maxchi0c2)then
        chick=i
        maxchi0c2=maxchi0c
     end if
  end do
  do i=1,Nkall !get kmap footnote
     if(invk(2,i)==0)then
        if(invk(1,i)==chisk)chiskall=i
        if(invk(1,i)==chick)chickall=i
     end if
  end do
  print'(A3,3I4,F12.8)','SDW',kmap(:,chiskall),maxchi0s2
  print'(A3,3I4,F12.8)','CDW',kmap(:,chickall),maxchi0c2
  maxchi0s_out=maxchi0s2
end subroutine ckchi_impl

subroutine get_chi0(chi,Smat,Cmat,Gk,kmap,invk,olist,temp,Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi,maxchi0s_out) bind(C)
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
  real(c_double),intent(in):: temp
  real(c_double),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(c_double),intent(out),dimension(Nk,Nw,Nchi,Nchi):: chi
  real(c_double),intent(out):: maxchi0s_out

  integer(c_int32_t),dimension(Nchi,Nchi,2)::chi_map
  integer(c_int32_t),dimension(Nchi*(Nchi+1)/2,2)::irr_chi

  call get_chi_map(chi_map,irr_chi,olist,Nchi)
  call get_chi0_conv(chi,Gk,kmap,invk,irr_chi,chi_map,olist,temp,Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi)
  call ckchi_impl(chi,Smat,Cmat,kmap,invk,Nk,Nkall,Nchi,Nw,maxchi0s_out)
end subroutine get_chi0

subroutine get_chi_map(chi_map,irr_chi,olist,Nchi)
  !> This function generate index of exchange symmetry chi1234(q,iw)=chi*4321(q,iw).
  !> This symmetry can use system has no spin dependence and TRS.
  !!@param chi_map,out: mapping list of chi index
  !!@param irr_chi,out: irreducible index of chii
  !!@param    olist,in: orbital index list of chi index
  !!@param     Nchi,in: Number of chi index
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nchi
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int32_t),intent(out),dimension(Nchi,Nchi,2):: chi_map
  integer(c_int32_t),intent(out),dimension(Nchi*(Nchi+1)/2,2):: irr_chi

  integer(c_int32_t) l1,m1,l2,m2,iter
  integer(c_int32_t),dimension(Nchi,Nchi):: chi_irr
  integer(c_int32_t),dimension(4):: tmp1,tmp2
  logical ck

  chi_irr(:,:)=0
  chi_map(:,:,:)=0
  do l1=1,Nchi
     tmp1(1)=olist(l1,1) !1 of 1234
     tmp1(2)=olist(l1,2) !2 of 1234
     do m1=1,Nchi
        tmp1(3)=olist(m1,1) !3 of 1234
        tmp1(4)=olist(m1,2) !4 of 1234
        ck=.false.
        do l2=1,Nchi
           tmp2(3)=olist(l2,2) !2 of 4321
           tmp2(4)=olist(l2,1) !1 of 4321
           do m2=1,Nchi
              tmp2(1)=olist(m2,2) !4 of 4321
              tmp2(2)=olist(m2,1) !3 of 4321
              if(sum(abs(tmp1(:)-tmp2(:)))==0)then !4321 correspond to 1234
                 chi_map(m1,l1,1)=m2
                 chi_map(m1,l1,2)=l2
                 if(chi_irr(m2,l2)==0)then !4321 is not irreducible
                    chi_irr(m1,l1)=1 !1 is irreducible index
                 end if
                 ck=.true.
                 exit
              end if
           end do
           if(ck)exit
        end do
     end do
  end do

  !get list of irreducible chi index
  iter=1
  do l1=1,Nchi
     do m1=1,Nchi
        if(chi_irr(m1,l1)==1)then
           irr_chi(iter,1)=l1
           irr_chi(iter,2)=m1
           iter=iter+1
        end if
     end do
  end do
end subroutine get_chi_map

subroutine chi0_threshold(chi,Nk,Nw,Nchi)
  !> Zero out numerically tiny real/imaginary parts of chi (cosmetic cleanup,
  !> same eps as the original get_chi0_conv).
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,Nchi
  complex(c_double),intent(inout),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(c_int32_t) i,j,l,m
  real(c_double),parameter:: eps=1.0d-9

  !$omp parallel do collapse(2) private(i,j)
  do l=1,Nchi
     do m=1,Nchi
        do j=1,Nw
           do i=1,Nk
              if(abs(dble(chi(i,j,m,l)))<eps) chi(i,j,m,l)=cmplx(0.0d0,imag(chi(i,j,m,l)),kind=c_double)
              if(abs(imag(chi(i,j,m,l)))<eps) chi(i,j,m,l)=cmplx(dble(chi(i,j,m,l)),0.0d0,kind=c_double)
           end do
        end do
     end do
  end do
  !$omp end parallel do
end subroutine chi0_threshold

subroutine chi0_conv_acc(chi,Gk,kmap,invk,irr_chi,chi_map,olist,temp,coef,&
     Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi)
  !> Accumulate coef * (FFT-convolution bubble of Gk) into chi.
  !> Body of the original get_chi0_conv, made additive (chi += coef*conv) so the
  !> tail-corrected driver can combine conv[G] - conv[G0] + analytic reference
  !> without a second chi-sized work array.  The small-value thresholding is
  !> applied by the callers (chi0_threshold) after the final combination.
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nw,Norb,Nchi,Nkall,Nk,Nx,Ny,Nz
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
  integer(c_int32_t),intent(in),dimension(Nchi,Nchi,2):: chi_map
  integer(c_int32_t),intent(in),dimension(Nchi*(Nchi+1)/2,2):: irr_chi
  real(c_double),intent(in):: temp,coef
  complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(c_double),intent(inout),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(c_int32_t) i,j,k,l,m,n,iorb
  integer(c_int32_t) ii(0:Nx-1),ij(0:Ny-1),ik(0:Nz-1),iw(2*Nw)
  real(c_double) weight
  complex(c_double),dimension(0:Nx-1,0:Ny-1,0:Nz-1,2*Nw):: tmp,tmpgk13,tmpgk42

  weight=coef*temp/dble(Nkall)
  ii(0)=0
  ij(0)=0
  ik(0)=0
  iw(1)=1
  !$omp parallel
  !$omp do
  do i=1,Nx-1
     ii(i)=Nx-i
  end do
  !$omp end do
  !$omp do
  do i=1,Ny-1
     ij(i)=Ny-i
  end do
  !$omp end do
  !$omp do
  do i=1,Nz-1
     ik(i)=Nz-i
  end do
  !$omp end do
  !$omp do
  do i=2,2*Nw
     iw(i)=2*Nw-i+2
  end do
  !$omp end do
  !$omp end parallel
  call init_fft_plans(tmpgk13,tmp,Nx,Ny,Nz,2*Nw)
  orb_loop:do iorb=1,Nchi*(Nchi+1)/2
     l=irr_chi(iorb,1)
     m=irr_chi(iorb,2)
     !use symmetry G^lm(k,iw)=G^ml(k,-iw) from Hermitian symmetry of Hamiltonian
     !$omp parallel do private(i)
     w_loop_Gk_to_tmp:do j=1,Nw
        k_loop_Gk_to_tmp:do i=1,Nkall
           if(invk(2,i)==0)then
              !iw
              tmpgk13(kmap(1,i),kmap(2,i),kmap(3,i),j)=Gk(invk(1,i),j,olist(l,1),olist(m,1)) !G13(k,iw)
              tmpgk42(kmap(1,i),kmap(2,i),kmap(3,i),j)=Gk(invk(1,i),j,olist(m,2),olist(l,2)) !G42(k,iw)
              !-iw
              tmpgk13(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=conjg(Gk(invk(1,i),j,olist(m,1),olist(l,1))) !G13(k,-iw)=G^*31(k,iw)
              tmpgk42(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=conjg(Gk(invk(1,i),j,olist(l,2),olist(m,2))) !G42(k,-iw)=G^*24(k,iw)
           else if(invk(2,i)==1)then
              !iw
              tmpgk13(kmap(1,i),kmap(2,i),kmap(3,i),j)=Gk(invk(1,i),j,olist(m,1),olist(l,1)) !G13(-k,iw)=G^31(k,iw)
              tmpgk42(kmap(1,i),kmap(2,i),kmap(3,i),j)=Gk(invk(1,i),j,olist(l,2),olist(m,2)) !G42(-k,iw)=G^24(k,iw)
              !-iw
              tmpgk13(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=conjg(Gk(invk(1,i),j,olist(l,1),olist(m,1))) !G13(-k,-iw)=G^*13(k,iw)
              tmpgk42(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=conjg(Gk(invk(1,i),j,olist(m,2),olist(l,2))) !G42(-k,-iw)=G^*42(k,iw)
           end if
        end do k_loop_Gk_to_tmp
     end do w_loop_Gk_to_tmp
     !$omp end parallel do
     call FFT(tmpgk13,tmp,Nx,Ny,Nz,2*Nw,.true.)
     call FFT(tmpgk42,tmp,Nx,Ny,Nz,2*Nw,.true.)
     !calculate G(r)G(-r)
     !$omp parallel do private(i,j,k)
     w_loop_conv:do n=1,2*Nw
        z_loop:do k=0,Nz-1
           y_loop:do j=0,Ny-1
              !$omp simd
              x_loop:do i=0,Nx-1
                 tmp(i,j,k,n)=-tmpgk13(i,j,k,n)*tmpgk42(ii(i),ij(j),ik(k),iw(n))
              end do x_loop
              !$omp end simd
           end do y_loop
        end do z_loop
     end do w_loop_conv
     !$omp end parallel do
     call FFT(tmp,tmpgk13,Nx,Ny,Nz,2*Nw,.false.)
     !$omp parallel do private(i,j)
     w_loop_tmp_to_chi:do j=1,Nw
        k_loop_tmp_to_chi:do i=1,Nkall
           if(invk(2,i)==0)then
              chi(invk(1,i),j,m,l)=chi(invk(1,i),j,m,l)+tmp(kmap(1,i),kmap(2,i),kmap(3,i),j)*weight
           end if
        end do k_loop_tmp_to_chi
     end do w_loop_tmp_to_chi
     !$omp end parallel do
     chi(:,:,chi_map(m,l,1),chi_map(m,l,2))=conjg(chi(:,:,m,l))
  end do orb_loop
  call destroy_fft_plans()
end subroutine chi0_conv_acc

subroutine get_chi0_conv(chi,Gk,kmap,invk,irr_chi,chi_map,olist,temp,&
     Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi) bind(C,name='get_chi0_conv_')
  !> This function obtains chi_0 using convolution (sharp Matsubara cutoff).
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nw,Norb,Nchi,Nkall,Nk,Nx,Ny,Nz
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
  integer(c_int32_t),intent(in),dimension(Nchi,Nchi,2):: chi_map
  integer(c_int32_t),intent(in),dimension(Nchi*(Nchi+1)/2,2):: irr_chi
  real(c_double),intent(in):: temp
  complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(c_double),intent(out),dimension(Nk,Nw,Nchi,Nchi):: chi

  !$omp parallel workshare
  chi(:,:,:,:)=(0.0d0,0.0d0)
  !$omp end parallel workshare
  call chi0_conv_acc(chi,Gk,kmap,invk,irr_chi,chi_map,olist,temp,1.0d0,&
       Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi)
  call chi0_threshold(chi,Nk,Nw,Nchi)
end subroutine get_chi0_conv

subroutine chi0_ref_tau(chi,eig,uni,mu,kmap,invk,irr_chi,chi_map,olist,temp,&
     Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi)
  !> Accumulate into chi the ANALYTIC reference bubble chi0[G0,G0](q,i nu_m).
  !>
  !> The non-interacting particle-hole bubble is evaluated through the exact
  !> imaginary-time product
  !>     chi0(q,tau) = -(1/N) sum_k G0(k+q,tau) G0(k,-tau)
  !> with G0(k,tau) built analytically per band (no Matsubara truncation), the
  !> k-convolution done by batched spatial FFT and the bosonic tau -> i nu_m
  !> transform by DFT on the uniform tau grid (M = 2*Nw points).
  !>
  !> Companion of chi0_conv_acc inside chi0_tail_impl:
  !>     chi0 = conv[G] - conv[G0] + (this reference)
  !> By bilinearity the FFT then only ever sees the fast-decaying residual
  !> (G-G0 ~ 1/w^2), removing the O(1/Nw) sharp-cutoff tail error of the plain
  !> convolution (residual O(1/Nw^2)).
  !>
  !> tau=0 jump handling: for orbital pairs with a nonzero equal-time
  !> commutator the bubble chi0(tau) has a step at tau=0 (the 1/(i nu)
  !> high-frequency component of chi0).  The jump J is computed analytically
  !> from the free density matrix and handled by the sawtooth trick: the DFT
  !> of J*(1/2 - tau/beta) is replaced by its exact transform -J/(i nu_m).
  !>
  !> Assumes TRS without SOC (H(-k) = H(k)^*), the same restriction as
  !> get_chi0_conv / gen_irr_k_TRS.
  use constant
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nw,Norb,Nchi,Nkall,Nk,Nx,Ny,Nz
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
  integer(c_int32_t),intent(in),dimension(Nchi,Nchi,2):: chi_map
  integer(c_int32_t),intent(in),dimension(Nchi*(Nchi+1)/2,2):: irr_chi
  real(c_double),intent(in):: temp,mu
  real(c_double),intent(in),dimension(Norb,Nk):: eig
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(c_double),intent(inout),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(c_int32_t) i,j,k,l,m,n,s,iorb,kk,M32,a13,b13,a42,b42
  integer(c_int32_t) ii(0:Nx-1),ij(0:Ny-1),ik(0:Nz-1)
  real(c_double) beta,dtau,xi,tau,fn,weight,nu
  complex(c_double) Jlm,zsumA,zsumB,ph
  real(c_double),allocatable:: gp(:,:,:),gm(:,:,:),ff(:,:)
  complex(c_double),allocatable:: U13(:,:),U42(:,:),rho(:,:,:)
  complex(c_double),dimension(Norb,Norb):: rhobar
  complex(c_double),dimension(Nw):: sawc
  complex(c_double),dimension(0:Nx-1,0:Ny-1,0:Nz-1,2*Nw):: tmpA,tmpB,tmp

  M32=int(2*Nw,c_int32_t)
  beta=1.0d0/temp
  dtau=beta/dble(M32)
  weight=1.0d0/(temp*dble(Nkall))

  ii(0)=0
  ij(0)=0
  ik(0)=0
  do i=1,Nx-1
     ii(i)=Nx-i
  end do
  do i=1,Ny-1
     ij(i)=Ny-i
  end do
  do i=1,Nz-1
     ik(i)=Nz-i
  end do

  ! --- stable analytic band kernels on the tau grid ------------------------
  ! gp(s,n,k) = g_n(k, tau_{s-1}) = -e^{-xi tau}(1-f)   (0 <= tau < beta)
  ! gm(s,n,k) = g_n(k,-tau_{s-1}) = +e^{ xi tau} f
  allocate(gp(2*Nw,Norb,Nk),gm(2*Nw,Norb,Nk),ff(Norb,Nk))
  !$omp parallel do private(n,s,xi,tau,fn)
  do i=1,Nk
     do n=1,Norb
        xi=eig(n,i)-mu
        if(xi>=0.0d0)then
           fn=exp(-beta*xi)/(1.0d0+exp(-beta*xi))
        else
           fn=1.0d0/(1.0d0+exp(beta*xi))
        end if
        ff(n,i)=fn
        do s=1,M32
           tau=dble(s-1)*dtau
           if(xi>=0.0d0)then
              gp(s,n,i)=-exp(-xi*tau)/(1.0d0+exp(-beta*xi))
              gm(s,n,i)= exp(-xi*(beta-tau))/(1.0d0+exp(-beta*xi))
           else
              gp(s,n,i)=-exp(xi*(beta-tau))/(1.0d0+exp(beta*xi))
              gm(s,n,i)= exp(xi*tau)/(1.0d0+exp(beta*xi))
           end if
        end do
     end do
  end do
  !$omp end parallel do

  ! --- BZ-averaged free density matrix (for the tau=0 jump) ----------------
  allocate(rho(Norb,Norb,Nk))
  !$omp parallel do private(l,m,n)
  do i=1,Nk
     do l=1,Norb
        do m=1,Norb
           rho(m,l,i)=(0.0d0,0.0d0)
           do n=1,Norb
              rho(m,l,i)=rho(m,l,i)+uni(m,n,i)*conjg(uni(l,n,i))*ff(n,i)
           end do
        end do
     end do
  end do
  !$omp end parallel do
  rhobar(:,:)=(0.0d0,0.0d0)
  do i=1,Nkall
     kk=int(invk(1,i),c_int32_t)
     if(invk(2,i)==0)then
        rhobar(:,:)=rhobar(:,:)+rho(:,:,kk)
     else if(invk(2,i)==1)then
        rhobar(:,:)=rhobar(:,:)+transpose(rho(:,:,kk))
     end if
  end do
  rhobar(:,:)=rhobar(:,:)/dble(Nkall)
  deallocate(rho)

  ! --- sawtooth: exact transform minus its rectangle-rule DFT --------------
  ! saw(tau) = 1/2 - tau/beta carries the unit tau=0 jump; exact coefficients
  ! are -1/(i nu_m) (0 at nu=0).  sawc = exact - numeric replaces the poorly
  ! converging jump part of the sampled bubble by its analytic value.
  do j=1,Nw
     if(j==1)then
        sawc(j)=(0.0d0,0.0d0)
     else
        nu=2.0d0*pi*dble(j-1)*temp
        sawc(j)=cmplx(0.0d0,1.0d0/nu,kind=c_double)
     end if
     do s=0,M32-1
        ph=exp(cmplx(0.0d0,2.0d0*pi*dble(j-1)*dble(s)/dble(M32),kind=c_double))
        sawc(j)=sawc(j)-dtau*ph*(0.5d0-dble(s)/dble(M32))
     end do
  end do

  allocate(U13(Norb,Nk),U42(Norb,Nk))
  call init_fft_plans(tmpA,tmp,Nx,Ny,Nz,2*Nw)
  orb_loop:do iorb=1,Nchi*(Nchi+1)/2
     l=irr_chi(iorb,1)
     m=irr_chi(iorb,2)
     a13=int(olist(l,1),c_int32_t)  ! G13 = G0_{a13,b13}(k+q, tau)
     b13=int(olist(m,1),c_int32_t)
     a42=int(olist(m,2),c_int32_t)  ! G42 = G0_{a42,b42}(k, -tau)
     b42=int(olist(l,2),c_int32_t)
     !$omp parallel do private(n)
     do i=1,Nk
        do n=1,Norb
           U13(n,i)=uni(a13,n,i)*conjg(uni(b13,n,i))
           U42(n,i)=uni(a42,n,i)*conjg(uni(b42,n,i))
        end do
     end do
     !$omp end parallel do
     ! full-grid tau kernels; G0(-k,tau) = G0(k,tau)^T (TRS) -> conj(U) coefficients
     !$omp parallel do private(i,kk,n,zsumA,zsumB)
     do s=1,M32
        do i=1,Nkall
           kk=int(invk(1,i),c_int32_t)
           zsumA=(0.0d0,0.0d0)
           zsumB=(0.0d0,0.0d0)
           if(invk(2,i)==0)then
              do n=1,Norb
                 zsumA=zsumA+U13(n,kk)*gp(s,n,kk)
                 zsumB=zsumB+U42(n,kk)*gm(s,n,kk)
              end do
           else if(invk(2,i)==1)then
              do n=1,Norb
                 zsumA=zsumA+conjg(U13(n,kk))*gp(s,n,kk)
                 zsumB=zsumB+conjg(U42(n,kk))*gm(s,n,kk)
              end do
           end if
           tmpA(kmap(1,i),kmap(2,i),kmap(3,i),s)=zsumA
           tmpB(kmap(1,i),kmap(2,i),kmap(3,i),s)=zsumB
        end do
     end do
     !$omp end parallel do
     ! spatial-only FFT (tau axis untouched), product with spatial reversal,
     ! then one 4D backward = spatial q-transform + bosonic tau -> nu DFT
     call fft3d_batch(tmpA,tmp,Nx,Ny,Nz,2*Nw,.true.)
     call fft3d_batch(tmpB,tmp,Nx,Ny,Nz,2*Nw,.true.)
     !$omp parallel do private(i,j,k)
     do s=1,M32
        do k=0,Nz-1
           do j=0,Ny-1
              !$omp simd
              do i=0,Nx-1
                 tmp(i,j,k,s)=-tmpA(i,j,k,s)*tmpB(ii(i),ij(j),ik(k),s)
              end do
              !$omp end simd
           end do
        end do
     end do
     !$omp end parallel do
     call FFT(tmp,tmpB,Nx,Ny,Nz,2*Nw,.false.)
     ! tau=0 jump of this pair: J = delta_{13} rhobar_{42} - delta_{42} rhobar_{13}
     Jlm=(0.0d0,0.0d0)
     if(a13==b13) Jlm=Jlm+rhobar(a42,b42)
     if(a42==b42) Jlm=Jlm-rhobar(a13,b13)
     !$omp parallel do private(i,j)
     do j=1,Nw
        do i=1,Nkall
           if(invk(2,i)==0)then
              chi(invk(1,i),j,m,l)=chi(invk(1,i),j,m,l)&
                   +tmp(kmap(1,i),kmap(2,i),kmap(3,i),j)*weight+Jlm*sawc(j)
           end if
        end do
     end do
     !$omp end parallel do
     chi(:,:,chi_map(m,l,1),chi_map(m,l,2))=conjg(chi(:,:,m,l))
  end do orb_loop
  call destroy_fft_plans()
  deallocate(U13,U42,gp,gm,ff)
end subroutine chi0_ref_tau

subroutine chi0_tail_impl(chi,Gk,eig,uni,mu,kmap,invk,irr_chi,chi_map,olist,temp,&
     Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi)
  !> Tail-corrected chi0: chi0 = conv[G] - conv[G0] + chi0_ref[G0].
  !> By bilinearity of the bubble the FFT convolution only acts on terms
  !> containing dG = G - G0 (decaying like 1/w^2 or faster), while the slowly
  !> decaying 1/(iw) reference part is evaluated analytically in imaginary time
  !> (chi0_ref_tau).  Sharp-cutoff error: O(1/Nw) -> O(1/Nw^2).
  !> G0 is built from (eig, uni) at the SAME chemical potential mu as Gk so the
  !> 1/(iw) coefficients cancel exactly.
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nw,Norb,Nchi,Nkall,Nk,Nx,Ny,Nz
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
  integer(c_int32_t),intent(in),dimension(Nchi,Nchi,2):: chi_map
  integer(c_int32_t),intent(in),dimension(Nchi*(Nchi+1)/2,2):: irr_chi
  real(c_double),intent(in):: temp,mu
  real(c_double),intent(in),dimension(Norb,Nk):: eig
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(c_double),intent(out),dimension(Nk,Nw,Nchi,Nchi):: chi

  complex(c_double),allocatable:: G0(:,:,:,:)

  !$omp parallel workshare
  chi(:,:,:,:)=(0.0d0,0.0d0)
  !$omp end parallel workshare
  call chi0_conv_acc(chi,Gk,kmap,invk,irr_chi,chi_map,olist,temp,1.0d0,&
       Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi)
  allocate(G0(Nk,Nw,Norb,Norb))
  !$omp parallel workshare
  G0(:,:,:,:)=(0.0d0,0.0d0)
  !$omp end parallel workshare
  call gen_green0(G0,eig,uni,mu,temp,Nk,Nw,Norb)
  call chi0_conv_acc(chi,G0,kmap,invk,irr_chi,chi_map,olist,temp,-1.0d0,&
       Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi)
  deallocate(G0)
  call chi0_ref_tau(chi,eig,uni,mu,kmap,invk,irr_chi,chi_map,olist,temp,&
       Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi)
  call chi0_threshold(chi,Nk,Nw,Nchi)
end subroutine chi0_tail_impl

subroutine get_chi0_tail(chi,Smat,Cmat,Gk,eig,uni,kmap,invk,olist,mu,temp,&
     Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi,maxchi0s_out) bind(C)
  !> Tail-corrected variant of get_chi0 (see chi0_tail_impl); needs the band
  !> data (eig, uni) and the chemical potential of Gk in addition.
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
  real(c_double),intent(in):: temp,mu
  real(c_double),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  real(c_double),intent(in),dimension(Norb,Nk):: eig
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(c_double),intent(out),dimension(Nk,Nw,Nchi,Nchi):: chi
  real(c_double),intent(out):: maxchi0s_out

  integer(c_int32_t),dimension(Nchi,Nchi,2)::chi_map
  integer(c_int32_t),dimension(Nchi*(Nchi+1)/2,2)::irr_chi

  call get_chi_map(chi_map,irr_chi,olist,Nchi)
  call chi0_tail_impl(chi,Gk,eig,uni,mu,kmap,invk,irr_chi,chi_map,olist,temp,&
       Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi)
  call ckchi_impl(chi,Smat,Cmat,kmap,invk,Nk,Nkall,Nchi,Nw,maxchi0s_out)
end subroutine get_chi0_tail

subroutine get_chi0_sum(chi,Gk,klist,invk,irr_chi,chi_map,olist,temp,Nw,Nk,Nkall,Norb,Nchi) bind(C)
  !> It obtains chi_0 using summation. Its cost is O(Nk^2), so it is heavy. You should use get_chi0_conv.
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nw,Norb,Nchi,Nk,Nkall
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),dimension(3,Nkall):: invk
  integer(c_int32_t),intent(in),dimension(Nchi,Nchi,2):: chi_map
  integer(c_int32_t),intent(in),dimension(Nchi*(Nchi+1)/2,2):: irr_chi
  real(c_double),intent(in),dimension(3,Nkall):: klist
  real(c_double),intent(in):: temp
  complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(c_double),intent(out),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(c_int32_t) i,j,l,m,iq,iw,iorb
  integer(c_int32_t),dimension(2*Nw):: wshift
  integer(c_int64_t),dimension(Nkall):: qshift
  real(c_double) weight
  complex(c_double),dimension(Nkall,2*Nw):: tmpgk13,tmpgk42

  !$omp parallel workshare
  chi(:,:,:,:)=0.0d0
  !$omp end parallel workshare

  weight=temp/dble(Nkall)
  orb_loop:do iorb=1,Nchi*(Nchi+1)/2
     l=irr_chi(iorb,1)
     m=irr_chi(iorb,2)
     !$omp parallel do private(i)
     do j=1,Nw
        do i=1,Nkall
           if(invk(2,i)==0)then
              tmpgk13(i,j)=Gk(i,j,olist(m,1),olist(l,1))
              tmpgk42(i,j)=Gk(i,j,olist(l,2),olist(m,2))
              tmpgk13(i,2*Nw-j+1)=conjg(Gk(i,j,olist(l,1),olist(m,1)))
              tmpgk42(i,2*Nw-j+1)=conjg(Gk(i,j,olist(m,2),olist(l,2)))
           else if(invk(2,i)==1)then
              tmpgk13(i,j)=Gk(i,j,olist(l,1),olist(m,1))
              tmpgk42(i,j)=Gk(i,j,olist(m,2),olist(l,2))
              tmpgk13(i,2*Nw-j+1)=conjg(Gk(i,j,olist(m,1),olist(l,1)))
              tmpgk42(i,2*Nw-j+1)=conjg(Gk(i,j,olist(l,2),olist(m,2)))
           end if
        end do
     end do
     !$omp end parallel do
     wloop: do iw=1,Nw
        wl_loop: do j=1,2*Nw
           wshift(j)=mod(j+iw-1,2*Nw)
           if(wshift(j)==0)then
              wshift(j)=2*Nw
           end if
        end do wl_loop
        qloop: do iq=1,Nkall
           if(invk(2,iq)==0)then
              call get_qshift(klist(:,iq),klist,qshift,Nk)
              !$omp parallel do private(i) reduction(+: chi)
              wloop2: do j=1,2*Nw
                 !$omp simd
                 kloop: do i=1,Nkall
                    chi(invk(1,iq),iw,m,l)=chi(invk(1,iq),iw,m,l)-tmpgk13(i,j)*tmpgk42(qshift(i),wshift(j))
                 end do kloop
                 !$omp end simd
              end do wloop2
              !$omp end parallel do
           end if
        end do qloop
     end do wloop
     chi(:,:,chi_map(m,l,1),chi_map(m,l,2))=conjg(chi(:,:,m,l))
  end do orb_loop
  chi(:,:,:,:)=chi(:,:,:,:)*weight
end subroutine get_chi0_sum