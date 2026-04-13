subroutine get_iqshift(qpoint,klist,qshift,Nk) bind(C,name="get_iqshift_")
  !>shift k to -k+q
  !!@param  qpoint,in: q-vector
  !!@param   klist,in: list of all k-point
  !!@param qshift,out: list of footnote of klist that correspond to k+q shift
  !!@param      Nk,in: number of k-points
  use,intrinsic:: iso_fortran_env, only:int64,real64
  implicit none
  integer(int64),intent(in):: Nk
  real(real64),intent(in),dimension(3,Nk):: klist
  real(real64),intent(in),dimension(3):: qpoint
  integer(int64),intent(out),dimension(Nk):: qshift
  
  integer(int64) i,j,k,ck
  real(real64) tmp
  real(real64),dimension(3,Nk):: kqlist
  
  integer(int64) :: scale_int, i1, i2, i3, j1, j2, j3
  integer(int64), allocatable :: ht_key(:), ht_val(:)
  integer(int64) :: m, h, step
  integer(int64) :: jj
  integer(int64) :: key, key_t
  real(real64) :: tmpf(3)

  scale_int=1048576_int64
  m = max(3_int64, 2_int64*Nk + 3_int64)
  allocate(ht_key(m))
  allocate(ht_val(m))
  ht_key(:) = -1_int64
  ht_val(:) = -1_int64

  ! build hash table from klist (single-threaded)
  do jj=1,int(Nk)
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
        if(h>m) h = h - m
     end do
     ht_key(h) = key
     ht_val(h) = jj
  end do

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
        if(h>m) h = h - m
     end do
  end do

  deallocate(ht_key)
  deallocate(ht_val)
end subroutine get_iqshift

module calc_irr_phi
  use,intrinsic:: iso_fortran_env, only:int32,int64,real64
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
    integer(int64),intent(in):: Nk,Norb,Nchi
    integer(int64),intent(in),dimension(Nk):: qshift
    integer(int64),intent(in),dimension(Nchi,2):: ol
    real(real64),intent(in):: mu,temp,eps,idelta,w
    real(real64),intent(in),dimension(Norb,Nk):: eig,ffermi
    complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  
    integer(int32) i,j,k,l,m
    complex(real64) unitmp
    complex(real64),dimension(Nchi,Nchi):: phi,calc_phi
  
    phi(:,:)=0.0d0
    kloop: do k=1,Nk
       band1_loop: do l=1,Norb
          band2_loop: do m=1,Norb
             chiorb1_loop: do j=1,Nchi
                chiorb2_loop:do i=1,Nchi
                   if(abs(w-eig(m,k)-eig(l,qshift(k))+2.0d0*mu)<eps .and. idelta<eps)then
                      continue
                   else
                      unitmp=uni(ol(j,1),l,qshift(k))*conjg(uni(ol(i,1),l,qshift(k)))&
                           *uni(ol(i,2),m,k)*conjg(uni(ol(j,2),m,k))
                      phi(i,j)=phi(i,j)-unitmp*(1.0d0-ffermi(l,qshift(k))-ffermi(m,k))&
                           /cmplx(w-eig(m,k)-eig(l,qshift(k))+2.0d0*mu,idelta)
                   end if
                end do chiorb2_loop
             end do chiorb1_loop
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
  integer(int64),intent(in):: Nk,Norb,Nw,Nchi
  integer(int64),intent(in),dimension(Nk):: qshift
  integer(int64),intent(in),dimension(Nchi,2):: ol
  real(real64),intent(in):: temp,mu,eps,idelta
  real(real64),intent(in),dimension(Norb,Nk):: eig,ffermi
  real(real64),intent(in),dimension(Nw):: wl
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(real64),intent(out),dimension(Nchi,Nchi,Nw):: phi

  integer(int64) i

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
  use,intrinsic:: iso_fortran_env, only:int32,int64,real64
  implicit none
  integer(int64),intent(in):: Nchi,Nw,Norb
  integer(int64),intent(in),dimension(Nchi,2):: olist
  complex(real64),intent(in),dimension(Nchi,Nchi,Nw):: phi
  complex(real64),intent(out),dimension(Nw):: trphi
  complex(real64),intent(out),dimension(Norb+2,Nw):: phi_orb
  
  integer(int64) i,j,k

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
  integer(int64),intent(in):: Nx,Ny,Nk,Norb,Nchi
  integer(int64),intent(in),dimension(Nchi,2):: ol
  real(real64),intent(in):: ecut,idelta,eps,temp,mu
  real(real64),intent(in),dimension(3,Nk):: klist
  real(real64),intent(in),dimension(Norb,Nk):: eig,ffermi
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(real64),intent(out),dimension(Ny,Nx):: trphi
  logical(1),intent(in):: sw_omega

  integer(int32) info
  integer(int64) i,j,l,m,n
  integer(int64),dimension(Nk):: qshift
  integer(int32),dimension(Nchi):: ipiv
  real(real64) wre,wim
  real(real64),dimension(3):: qpoint
  complex(real64),dimension(Nchi,Nchi):: phi
  complex(real64),dimension(2*Nchi):: work

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
