subroutine get_iqshift(qpoint,klist,qshift,Nk) bind(C,name="get_iqshift_")
  !>shift k to -k+q
  !!@param  qpoint,in: q-vector
  !!@param   klist,in: list of all k-point
  !!@param qshift,out: list of footnote of klist that correspond to k+q shift
  !!@param      Nk,in: number of k-points
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk
  real(c_double),intent(in),dimension(3,Nk):: klist
  real(c_double),intent(in),dimension(3):: qpoint
  integer(c_int64_t),intent(out),dimension(Nk):: qshift
  
  integer(c_int64_t) i,j,k,ck
  real(c_double) tmp
  real(c_double),dimension(3,Nk):: kqlist
  
  integer(c_int64_t) :: scale_int, i1, i2, i3, j1, j2, j3
  integer(c_int64_t), allocatable :: ht_key(:), ht_val(:)
  integer(c_int64_t) :: m, h, step
  integer(c_int64_t) :: jj
  integer(c_int64_t) :: key, key_t
  real(c_double) :: tmpf(3)

  scale_int=1048576_c_int64_t
  m = max(3_c_int64_t, 2_c_int64_t*Nk + 3_c_int64_t)
  allocate(ht_key(m))
  allocate(ht_val(m))
  ht_key(:) = -1_c_int64_t
  ht_val(:) = -1_c_int64_t

  ! build hash table from klist (single-threaded)
  do jj=1,int(Nk)
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
        if(h>m) h = h - m
     end do
     ht_key(h) = key
     ht_val(h) = jj
  end do

  ! Lookup for -k+q: compute 1-k+q = (-k+q) mod 1 (folds into [0,1))
  !$omp parallel do private(i1,i2,i3,key_t,h,step,tmpf,j1,j2,j3)
  do i=1,Nk
     tmpf(:)=1.0d0-klist(:,i)+qpoint(:)
     do j=1,3
        if(tmpf(j)>=1.0d0) then
           tmpf(j)=tmpf(j)-1.0d0
        else if(tmpf(j)<0.0d0) then
           tmpf(j)=tmpf(j)+1.0d0
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
     do while(ht_key(h) /= -1_c_int64_t)
        if(ht_key(h) == key_t) then
           qshift(i) = ht_val(h)
           exit
        end if
        h = h + step
        if(h>m) h = h - m
     end do
  end do

  deallocate(ht_key)
  deallocate(ht_val)
end subroutine get_iqshift

module calc_irr_phi
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
contains
  function calc_phi(Nk,Norb,Nchi,uni,eig,ffermi,ol,mu,temp,qshift,w,idelta,eps)
    !> This function obtain irreducible sc susceptibility phi_0
    !!@param        Nk: The number of k-points
    !!@param      Norb: The number of orbitals
    !!@param      Nchi: The footnote of chi
    !!@param       uni: unitary matrix
    !!@param       eig: energies of bands
    !!@param    ffermi: fermi distribute function
    !!@param        ol: the list of the properties of orbitals at footnote of chi
    !!@param        mu: chemical potential
    !!@param      temp: temperature
    !!@param    qshift: q-shifted klist
    !!@param         w: frequency
    !!@param    idelta: dumping factor
    !!@param       eps: threshold of calculation value
    !!@return calc_phi: irreducible sc susceptibility matrix
    integer(c_int64_t),intent(in):: Nk,Norb,Nchi
    integer(c_int64_t),intent(in),dimension(Nk):: qshift
    integer(c_int64_t),intent(in),dimension(Nchi,2):: ol
    real(c_double),intent(in):: mu,temp,eps,idelta,w
    real(c_double),intent(in),dimension(Norb,Nk):: eig,ffermi
    complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni
  
    integer(c_int32_t) i,j,k,l,m,nchi32
    complex(c_double) weight
    complex(c_double),dimension(Nchi):: A_vec,B_vec
    complex(c_double),dimension(Nchi,Nchi):: phi,calc_phi

    nchi32=int(Nchi,c_int32_t)
    phi(:,:)=0.0d0
    kloop: do k=1,Nk
       band1_loop: do l=1,Norb
          band2_loop: do m=1,Norb
             ! skip singular denominator
             if(abs(w-eig(m,k)-eig(l,qshift(k))+2.0d0*mu)<eps .and. idelta<eps)then
                cycle band2_loop
             end if
             ! Cooper pair propagator (particle-particle channel):
             ! weight = -(1-f_l(-k+q) - f_m(k)) / (ω - E_m(k) - E_l(-k+q) + 2μ)
             weight=-(1.0d0-ffermi(l,qshift(k))-ffermi(m,k))&
                  /cmplx(w-eig(m,k)-eig(l,qshift(k))+2.0d0*mu,idelta,kind=c_double)
             ! A_vec(j) = uni(ol(j,1),l,qshift(k)) * conjg(uni(ol(j,2),m,k))
             do j=1,Nchi
                A_vec(j)=uni(ol(j,1),l,qshift(k))*conjg(uni(ol(j,2),m,k))
             end do
             ! B_vec(i) = conjg(uni(ol(i,1),l,qshift(k))) * uni(ol(i,2),m,k)
             do i=1,Nchi
                B_vec(i)=conjg(uni(ol(i,1),l,qshift(k)))*uni(ol(i,2),m,k)
             end do
             ! phi(i,j) += weight * B_vec(i) * A_vec(j)
             call zgeru(nchi32,nchi32,weight,B_vec,1,A_vec,1,phi,nchi32)
          end do band2_loop
       end do band1_loop
    end do kloop
    calc_phi=phi(:,:)/Nk
  end function calc_phi
end module calc_irr_phi

subroutine get_phi_irr(phi,uni,eig,ffermi,qshift,ol,wl,Nchi,Norb,Nk,Nw,idelta,eps,mu,temp) bind(C)
  !> This function obtains irreducible sc susceptibility at q-point
  !!@param   phi,out: irreducible sc susceptibility
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
  !!@param     mu,in: chemical potential
  !!@param   temp,in: temperature
  use calc_irr_phi
  implicit none
  integer(c_int64_t),intent(in):: Nk,Norb,Nw,Nchi
  integer(c_int64_t),intent(in),dimension(Nk):: qshift
  integer(c_int64_t),intent(in),dimension(Nchi,2):: ol
  real(c_double),intent(in):: temp,mu,eps,idelta
  real(c_double),intent(in),dimension(Norb,Nk):: eig,ffermi
  real(c_double),intent(in),dimension(Nw):: wl
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(c_double),intent(out),dimension(Nchi,Nchi,Nw):: phi

  integer(c_int64_t) i

  !$omp parallel do private(i)
  wloop: do i=1,Nw
     phi(:,:,i)=calc_phi(Nk,Norb,Nchi,uni,eig,ffermi,ol,mu,temp,qshift,wl(i),idelta,eps)
  end do wloop
  !$omp end parallel do
end subroutine get_phi_irr

subroutine get_tr_phi(trphi,phi_orb,phi,olist,Nw,Nchi,Norb) bind(C)
  !> This function take trace of phi
  !!@param   trphi,out: trace value of phi
  !!@param phi_orb,out: The array of phi values at orbital diagonal elements
  !!@param      phi,in: sc susceptibility
  !!@param    olist,in: the list of  properties of chis,chi0 footnote
  !!@param       Nw,in: The number of frequency mesh
  !!@param     Nchi,in: The number of footnote of chis,chi0
  !!@param     Norb,in: The number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nchi,Nw,Norb
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  complex(c_double),intent(in),dimension(Nchi,Nchi,Nw):: phi
  complex(c_double),intent(out),dimension(Nw):: trphi
  complex(c_double),intent(out),dimension(Norb+2,Nw):: phi_orb
  
  integer(c_int64_t) i,j,k

  trphi(:)     = (0.0d0, 0.0d0)
  phi_orb(:,:) = (0.0d0, 0.0d0)

  !$omp parallel do private(j,k)
  wloop:do i=1,Nw
     orb_lop1:do j=1,Nchi
        trphi(i)=trphi(i)+phi(j,j,i)
        if(olist(j,1)==olist(j,2))then
           orb_loop2:do k=1,Nchi
              if(olist(k,1)==olist(k,2))then
                 if(olist(j,1)==olist(k,1))then
                    phi_orb(olist(j,1),i)=phi(k,j,i)
                 end if
                 if(olist(k,1)==2 .and. olist(j,1)==3)phi_orb(Norb+1,i)=phi(k,j,i)
                 if(olist(k,1)==2 .and. olist(j,1)==4)phi_orb(Norb+2,i)=phi(k,j,i)
              end if
           end do orb_loop2
        end if
     end do orb_lop1
  end do wloop
  !$omp end parallel do
end subroutine get_tr_phi

subroutine phiq_map(trphi,uni,eig,ffermi,klist,ol,mu,temp,ecut,idelta,eps,Nx,Ny,Nk,Norb,Nchi,sw_omega) bind(C)
  !> This function obtains irreducible susceptibility at w=ecut
  !!@param   trphi,out: trace value of phi
  !!@param      uni,in: unitary matrix
  !!@param      eig,in: energies of bands
  !!@param   ffermi,in: fermi distribute function
  !!@param    klist,in: klist
  !!@param       ol,in: the list of the properties of orbitals at footnote of chi
  !!@param       mu,in: chemical potential
  !!@param     temp,in: temperature
  !!@param     ecut,in: energy value of cut energy plane
  !!@param   idelta,in: dumping factor
  !!@param      eps,in: threshold of calculation value
  !!@param       Nx,in: The number of kx mesh
  !!@param       Ny,in: The number of ky mesh
  !!@param       Nk,in: The number of k-points
  !!@param     Norb,in: The number of orbitals
  !!@param     Nchi,in: The footnote of chi
  !!@param sw_omega,in: switch matsubara or real frequency
  use calc_irr_phi
  implicit none
  integer(c_int64_t),intent(in):: Nx,Ny,Nk,Norb,Nchi
  integer(c_int64_t),intent(in),dimension(Nchi,2):: ol
  real(c_double),intent(in):: ecut,idelta,eps,temp,mu
  real(c_double),intent(in),dimension(3,Nk):: klist
  real(c_double),intent(in),dimension(Norb,Nk):: eig,ffermi
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(c_double),intent(out),dimension(Ny,Nx):: trphi
  logical(1),intent(in):: sw_omega

  integer(c_int32_t) info
  integer(c_int64_t) i,j,l,m,n
  integer(c_int64_t),dimension(Nk):: qshift
  integer(c_int32_t),dimension(Nchi):: ipiv
  real(c_double) wre,wim
  real(c_double),dimension(3):: qpoint
  complex(c_double),dimension(Nchi,Nchi):: phi
  complex(c_double),dimension(2*Nchi):: work

  ! sw_omega=.true.: real frequency ω+iδ;  .false.: Matsubara iω_n (w=0, δ=0)
  if(sw_omega)then !set omega=cmplex(w,idelta)
     wre=ecut
     wim=idelta
  else             !set omega=cmplex(0,omega_n)
     wre=0.0d0
     wim=0.0d0
  end if
  !$omp parallel
  !$omp workshare
  trphi(:,:)=0.0d0
  !$omp end workshare
  !$omp do private(i,j,l,m,n,phi,qpoint,qshift,ipiv,work,info)
  do i=1,Nx
     do j=1,Ny
        qpoint(1)=dble(i-1)/Nx
        qpoint(2)=dble(j-1)/Ny
        qpoint(3)=0.0d0
        call get_iqshift(qpoint,klist,qshift,Nk)
        phi(:,:)=calc_phi(Nk,Norb,Nchi,uni,eig,ffermi,ol,mu,temp,qshift,wre,wim,eps)
        do l=1,Nchi
           trphi(j,i)=trphi(j,i)+phi(l,l)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine phiq_map
