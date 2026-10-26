module mkself_rn
  !> Electron number n(mu) for mkself's chemical-potential search (flex_scf's solve_mu).
  !> It is passed to solve_mu as an argument, so it is a module procedure reading state that
  !> mkself points at, not an internal procedure of mkself: passing an internal procedure needs
  !> a trampoline on the stack, i.e. an executable stack, which ld.lld refuses to link (LLVM
  !> flang) and glibc >= 2.41 refuses to dlopen.
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double
  implicit none
  integer(c_int64_t):: rn_Nk,rn_Nw,rn_Norb,rn_Nkall
  real(c_double):: rn_temp
  real(c_double),pointer,contiguous:: rn_eig(:,:)=>null(),rn_kmult(:)=>null()
  complex(c_double),pointer,contiguous:: rn_eig_gkinv(:,:,:)=>null()
contains
  subroutine get_rn(rn,rmu)
    real(c_double),intent(in):: rmu
    real(c_double),intent(out):: rn
    call get_rn_mkself(rn,rmu,rn_eig,rn_kmult,rn_eig_gkinv,rn_temp,rn_Nk,rn_Nw,rn_Norb,rn_Nkall)
  end subroutine get_rn
end module mkself_rn

subroutine gen_green0(Gk,eig,uni,mu,temp,Nk,Nw,Norb)
  !> This function obtains green function G0
  !!@param  Gk,out: green function
  !!@param  eig,in: energies at k-points
  !!@param  uni,in: unitary matrix
  !!@param   mu,in: chemical potential
  !!@param temp,in: Temperature
  !!@param   Nk,in: The number of k-points
  !!@param   Nw,in: The number of energies mesh
  !!@param Norb,in: The number of orbitals
  use freq_grid, only:fg_omega,FG_FERMION
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  use constant
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,norb
  real(c_double),intent(in):: mu,temp
  real(c_double),intent(in),dimension(Norb,Nk):: eig
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(c_double),intent(out),dimension(Nk,Nw,Norb,Norb):: Gk
  
  integer(c_int32_t) i,j,l,m,n
  complex(c_double) iw

  Gk(:,:,:,:)=0.0d0

  !$omp parallel private(l,m)
  do l=1,Norb
     do m=1,Norb
        !$omp do private(iw,i,j,n)
        wloop: do j=1,Nw !ien=pi(2l+1)/beta l=0,1,... beta=(kBT)^-1
           iw=cmplx(mu,fg_omega(FG_FERMION,j,temp),kind=c_double) !iω_j = μ + i(2j-1)πT (j=1: ω_0=πT)
           band_loop: do n=1,Norb
              kloop: do i=1,Nk
                 Gk(i,j,m,l)=Gk(i,j,m,l)+uni(m,n,i)*conjg(uni(l,n,i))/(iw-eig(n,i))
              end do kloop
           end do band_loop
        end do wloop
        !$omp end do
     end do
  end do
  !$omp end parallel
end subroutine gen_green0

subroutine gen_green_inv(Gk,self,hamk,mu,temp,Nk,Nw,Norb)
  !> This function obtains inverse of green function G^-1
  !!@param  Gk,out: inverse of green function
  !!@param self,in: self energies
  !!@param hamk,in: hamiltonian at k-point
  !!@param   mu,in: chemical potential
  !!@param temp,in: Temperature
  !!@param   Nk,in: The number of k-points
  !!@param   Nw,in: The number of energies mesh
  !!@param Norb,in: The number of orbitals
  use freq_grid, only:fg_omega,FG_FERMION
  use,intrinsic:: iso_c_binding, only: c_int64_t,c_double,c_int32_t
  use constant
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,Norb
  real(c_double),intent(in):: mu,temp
  complex(c_double),intent(in),dimension(norb,norb,Nk):: hamk
  complex(c_double),intent(in),dimension(Nk,Nw,norb,norb):: self
  complex(c_double),intent(out),dimension(Nk,Nw,norb,norb):: Gk

  integer(c_int32_t)i,j,l,m
  complex(c_double) iw

  !G^-1=G^-1_0-sigma (=iwI-Hk-sigma)
  !$omp parallel private(l,m)
  do l=1,Norb
     do m=1,Norb
        !$omp do private(iw,i)
        do j=1,Nw
           iw=cmplx(mu,fg_omega(FG_FERMION,j,temp),kind=c_double)
           do i=1,Nk
              Gk(i,j,m,l)=-(hamk(m,l,i)+self(i,j,m,l))
              if(l==m)then
                 Gk(i,j,l,l)=Gk(i,j,l,l)+iw
              end if
           end do
        end do
        !$omp end do
     end do
  end do
  !$omp end parallel
end subroutine gen_green_inv

subroutine gen_green_inv_from_eig(Gk,self,uni,eig,mu,temp,Nk,Nw,Norb) bind(C)
  !> This function obtains inverse of green function G^-1 form energy and unitary matrix
  !!@param  Gk,out: inverse of green function
  !!@param self,in: self energies
  !!@param  uni,in: unitary matrix
  !!@param  eig,in: energies at k-points
  !!@param   mu,in: chemical potential
  !!@param temp,in: Temperature
  !!@param   Nk,in: The number of k-points
  !!@param   Nw,in: The number of energies mesh
  !!@param Norb,in: The number of orbitals
  use freq_grid, only:fg_omega,FG_FERMION
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  use constant
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,norb
  real(c_double),intent(in):: mu,temp
  real(c_double),intent(in),dimension(norb,Nk):: eig
  complex(c_double),intent(in),dimension(norb,norb,Nk):: uni
  complex(c_double),intent(in),dimension(Nk,Nw,norb,norb):: self
  complex(c_double),intent(out),dimension(Nk,Nw,norb,norb):: Gk

  integer(c_int32_t)i,j,l,m,n
  complex(c_double) iw

  !$omp parallel private(l,m)
  !$omp workshare
  Gk(:,:,:,:)=0.0d0
  !$omp end workshare
  do l=1,norb
     do m=1,norb
        !$omp do private(iw,i,n)
        wloop: do j=1,Nw
           iw=cmplx(mu,fg_omega(FG_FERMION,j,temp),kind=c_double) !(2(j-1)+1)pi/beta
           kloop: do i=1,Nk
              band_loop: do n=1,norb
                 Gk(i,j,m,l)=Gk(i,j,m,l)+uni(m,n,i)*conjg(uni(l,n,i))*(iw-eig(n,i))
              end do band_loop
              Gk(i,j,m,l)=Gk(i,j,m,l)-self(i,j,m,l)
           end do kloop
        end do wloop
        !$omp end do
     end do
  end do
  !$omp end parallel
end subroutine gen_green_inv_from_eig

subroutine getinv(Gk,Nk,Nw,Norb)
  !> This function obtain green function form G^-1
  !!@param Gk,inout: green function
  !!@param    Nk,in: The number of k-points
  !!@param    Nw,in: The number of energies mesh
  !!@param  Norb,in: The number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,Norb
  complex(c_double),intent(inout),dimension(Nk,Nw,Norb,Norb):: Gk

  integer(c_int32_t) i,j,ipiv(Norb),info
  complex(c_double) tmp(Norb,Norb),work(2*Norb)

  !$omp parallel do private(i,j,tmp,work,ipiv,info)
  do i=1,Nw
     do j=1,Nk
        tmp(:,:)=Gk(j,i,:,:)
        call zgetrf(Norb,Norb,tmp,Norb,ipiv,info)
        if(info/=0)then; print*,'zgetrf failed in getinv: info=',info; stop; end if
        call zgetri(Norb,tmp,Norb,ipiv,work,2*Norb,info)
        if(info/=0)then; print*,'zgetri failed in getinv: info=',info; stop; end if
        Gk(j,i,:,:)=tmp(:,:)
     end do
  end do
  !$omp end parallel do
end subroutine getinv

subroutine calc_sigma(sigmak,Gk,Vsigma,Smat,Cmat,kmap,invk,olist,temp,Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz)
  !> This function obtain new gap function Delta=sum VF
  !!@param sigmak,out: self energy
  !!@param      Gk,in: green function
  !!@param  Vsigma,in: irreducible susceptibility and interaction
  !!@param    Smat,in: S-matrix
  !!@param    Cmat,in: C-matrix
  !!@param    kmap,in: property of k-point
  !!@param    invk,in: list of reverse k-points
  !!@param   olist,in: property of chi footnote
  !!@param    temp,in: Temperature
  !!@param   Nkall,in: Number of all k-points
  !!@param      Nk,in: Number of k-points
  !!@param      Nw,in: Number of Matsubara frequencies
  !!@param    Nchi,in: Number of footnote of chi
  !!@param    Norb,in: Number of orbitals
  !!@param      Nx,in: Number of kx mesh
  !!@param      Ny,in: Number of ky mesh
  !!@param      Nz,in: Number of kz mesh
  use numerics, only:drop_roundoff
  use freq_grid
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  real(c_double),intent(in):: temp
  real(c_double),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  complex(c_double),intent(in),dimension(Nk,Nw,Nchi,Nchi):: Vsigma
  complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(c_double),intent(out),dimension(Nk,Nw,Norb,Norb):: sigmak

  integer(c_int32_t) i,j,k,n,l,m
  integer(c_int64_t) jp,jm,j0
  complex(c_double) vinf,st
  real(c_double) weight
  complex(c_double),allocatable,dimension(:,:,:,:):: tmpVsigma,tmp,tmpgk

  allocate(tmpVsigma(0:Nx-1,0:Ny-1,0:Nz-1,fg_nt(Nw)),tmp(0:Nx-1,0:Ny-1,0:Nz-1,fg_nt(Nw)),&
       tmpgk(0:Nx-1,0:Ny-1,0:Nz-1,fg_nt(Nw)))
  weight=temp/dble(Nkall)
  !$omp parallel workshare
  sigmak(:,:,:,:)=0.0d0
  !$omp end parallel workshare
  call fg_plans_init(tmpVsigma,tmp,Nx,Ny,Nz,Nw)
  j0=fg_pos(FG_BOSON,1_c_int64_t)   ! v=0
  do l=1,Nchi
     do m=1,Nchi
        !$omp parallel
        !$omp do
        do i=1,Nkall !j=1>iw=0, j=Nw>iw=inf
           if(invk(2,i)==0)then
              tmpVsigma(kmap(1,i),kmap(2,i),kmap(3,i),j0)=Vsigma(invk(1,i),1,m,l) !Vsigma(k,iw)
           else if(invk(2,i)==1)then
              tmpVsigma(kmap(1,i),kmap(2,i),kmap(3,i),j0)=Vsigma(invk(1,i),1,l,m) !Vsigma(-k,iw)
           end if
        end do
        !$omp end do
        !$omp do private(i,jp,jm)
        do j=2,Nw
           jp=fg_pos(FG_BOSON,j)
           jm=fg_neg(FG_BOSON,j,Nw)
           do i=1,Nkall
              if(invk(2,i)==0)then
                 tmpVsigma(kmap(1,i),kmap(2,i),kmap(3,i),jp)=Vsigma(invk(1,i),j,m,l) !Vsigma(k,iw)
                 tmpVsigma(kmap(1,i),kmap(2,i),kmap(3,i),jm)=conjg(Vsigma(invk(1,i),j,l,m)) !Vsigma(k,-iw)ml=Vsigma^*lm(k,iw)
              else if(invk(2,i)==1)then
                 tmpVsigma(kmap(1,i),kmap(2,i),kmap(3,i),jp)=Vsigma(invk(1,i),j,l,m) !Vsigma(-k,iw)=Vsigma(-k,iw)^T
                 tmpVsigma(kmap(1,i),kmap(2,i),kmap(3,i),jm)=conjg(Vsigma(invk(1,i),j,m,l)) !Vsigma(-k,-iw)=Vsigma^*(k,iw)
              end if
           end do
        end do
        !$omp end do
        !$omp do private(i,jp,jm)
        do j=1,Nw
           jp=fg_pos(FG_FERMION,j)
           jm=fg_neg(FG_FERMION,j,Nw)
           do i=1,Nkall
              if(invk(2,i)==0)then
                 tmpgk(kmap(1,i),kmap(2,i),kmap(3,i),jp)=Gk(invk(1,i),j,olist(m,2),olist(l,2)) !G42(iw)
                 tmpgk(kmap(1,i),kmap(2,i),kmap(3,i),jm)=conjg(Gk(invk(1,i),j,olist(l,2),olist(m,2))) !G42(k,-iw)=G^*24(k,iw)
              else if(invk(2,i)==1)then
                 tmpgk(kmap(1,i),kmap(2,i),kmap(3,i),jp)=Gk(invk(1,i),j,olist(l,2),olist(m,2)) !G42(-k,iw)=G24(k,iw)
                 tmpgk(kmap(1,i),kmap(2,i),kmap(3,i),jm)=conjg(Gk(invk(1,i),j,olist(m,2),olist(l,2))) !G42(-k,-iw)=G^*42(k,iw)
              end if
           end do
        end do
        !$omp end do
        !$omp end parallel
        vinf=cmplx(1.5d0*Smat(m,l)-0.5d0*Cmat(m,l),0.0d0,kind=c_double)   ! bare vertex (v -> infinity)
        call fg_set_static(tmpVsigma,vinf,Nx,Ny,Nz,Nw)
        call fg_forward(FG_BOSON,tmpVsigma,tmp,Nx,Ny,Nz,Nw)
        st=fg_static_term(tmpgk,Nx,Ny,Nz,Nw)   ! 0 on the uniform grid
        call fg_forward(FG_FERMION,tmpgk,tmp,Nx,Ny,Nz,Nw)
        !$omp parallel
        !$omp do private(i,j,k,n)
        do n=1,int(fg_nt(Nw),c_int32_t)
           do k=0,Nz-1
              do j=0,Ny-1
                 do i=0,Nx-1
                    tmp(i,j,k,n)=tmpVsigma(i,j,k,n)*tmpgk(i,j,k,n)
                 end do
              end do
           end do
        end do
        !$omp end do
        !$omp workshare
      tmpgk=cmplx(0.0d0,0.0d0,kind=c_double)
        !$omp end workshare
        !$omp end parallel
        call fg_backward(FG_FERMION,tmp,tmpgk,Nx,Ny,Nz,Nw)
        !$omp parallel do private(i,j,jp)
        do j=1,Nw
           jp=fg_pos(FG_FERMION,j)
           do i=1,Nkall
              if(invk(2,i)==0)then
                 sigmak(invk(1,i),j,olist(m,1),olist(l,1))=sigmak(invk(1,i),j,olist(m,1),olist(l,1))&
                      +tmp(kmap(1,i),kmap(2,i),kmap(3,i),jp)
              end if
           end do
        end do
        !$omp end parallel do
        if(.not.fg_uniform())then   ! IR: the bare vertex was kept out of the transform
           sigmak(:,:,olist(m,1),olist(l,1))=sigmak(:,:,olist(m,1),olist(l,1))+vinf*st
        end if
     end do
  end do

  call fg_plans_free()
  !$omp parallel workshare
  sigmak(:,:,:,:)=sigmak(:,:,:,:)*weight
  !$omp end parallel workshare
  call drop_roundoff(sigmak,Nk*Nw*Norb*Norb)
end subroutine calc_sigma

subroutine mkself(sigmak,mu,Smat,Cmat,kmap,invk,olist,hamk,eig,uni,mu_init,rfill,temp,&
     scf_loop,pp,eps,Nkall,Nk,Nw,Norb,Nchi,Nx,Ny,Nz,sub_sigma,sw_out,sw_in,m_diis,sw_rescale,&
     sw_tail,sigma_scale) bind(C)
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  use flex_scf, only:solve_mu,diis_state,diis_init,diis_free,mix_sigma
  use mkself_rn
  use freq_grid, only:fg_uniform
  implicit none
  integer(c_int64_t),intent(in):: Nw,Norb,Nchi,Nkall,Nk,Nx,Ny,Nz,scf_loop,m_diis
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
  logical(1),intent(in):: sw_in,sw_out,sw_rescale
  logical(1),intent(in):: sw_tail  !< tail-corrected chi0 (chi0_tail_impl) instead of the sharp-cutoff convolution
  integer(c_int64_t),intent(in):: sub_sigma
  real(c_double),intent(in):: temp,eps,pp,rfill,mu_init
  real(c_double),intent(in):: sigma_scale  !< multiply the sigma seed loaded via sw_in (U-annealing: ~(U_new/U_old)^2)
  real(c_double),intent(in),target,dimension(Norb,Nk):: eig
  real(c_double),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  real(c_double),intent(out):: mu
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni,hamk
  complex(c_double),intent(out),dimension(Nk,Nw,Norb,Norb):: sigmak

  integer(c_int32_t) scf_i
  integer(c_int32_t),dimension(Nchi,Nchi,2)::chi_map
  integer(c_int32_t),dimension(Nchi*(Nchi+1)/2,2)::irr_chi
  real(c_double)esterr,mu_OLD,eps_sgm,maxchi0s_global,rescale_t
  complex(c_double),allocatable,dimension(:,:,:,:):: Gk,sigmak0
  complex(c_double),allocatable,dimension(:,:,:,:):: chi
  ! DIIS
  type(diis_state):: diis   !! Pulay history for the sigma mixing (module flex_scf)
  ! bracket cache for renew_mu
  real(c_double):: muS_cache,muL_cache
  logical:: bracket_valid
  ! eigenvalues of H(k)+Sigma(k,iw) precomputed before mu search
  complex(c_double),allocatable,target:: eig_gkinv(:,:,:)
  ! multiplicity of each irreducible k in the full BZ. The TRS reduction pairs k with -k,
  ! but the time-reversal invariant momenta (Gamma, X, M, Z, ...) are their own partner and
  ! therefore have multiplicity 1 while every other k has 2. An unweighted average over the
  ! wedge would over-weight the TRIMs; chi0 and Sigma already expand to the full BZ, so the
  ! filling has to be weighted the same way to stay consistent.
  real(c_double),allocatable,target,dimension(:):: kmult
  integer(c_int32_t) ik_mult

  if(.not.fg_uniform())then
     ! sub_sigma=2 averages over a uniform grid, and sigma.bin is read back as one (plibs.regrid_sigma_bin)
     if(sub_sigma==2 .or. sw_in .or. sw_out)then
        print'(A)','ERROR(mkself): sub_sigma=2 and sw_in/sw_out are for the uniform Matsubara grid only'
        stop
     end if
  end if
  allocate(Gk(Nk,Nw,Norb,Norb),sigmak0(Nk,Nw,Norb,Norb),chi(Nk,Nw,Nchi,Nchi),kmult(Nk))
  kmult(:)=0.0d0
  do ik_mult=1,int(Nkall,c_int32_t)
     kmult(invk(1,ik_mult))=kmult(invk(1,ik_mult))+1.0d0
  end do

  eps_sgm=1.0d-10
  mu=mu_init
  if(sw_in)then
     print*,"load self"
     call io_sigma(.false.)
     if(sigma_scale/=1.0d0)then
        print'(A,F10.6)','[FLEX] scale loaded sigma seed by',sigma_scale
        !$omp parallel workshare
        sigmak(:,:,:,:)=sigmak(:,:,:,:)*sigma_scale
        !$omp end parallel workshare
     end if
     !$omp parallel
     !$omp workshare
     sigmak0(:,:,:,:)=sigmak(:,:,:,:)
     !$omp end workshare
     !$omp end parallel
     call gen_green_inv(Gk,sigmak,hamk,mu,temp,Nk,Nw,Norb)
     call getinv(Gk,Nk,Nw,Norb)
  else
     mu_old=mu*1.2
     !$omp parallel workshare
     sigmak0(:,:,:,:)=0.0d0
     Gk(:,:,:,:)=0.0d0 !gen_green0 need to initialization of Gk
     !$omp end parallel workshare
     call gen_green0(Gk,eig,uni,mu,temp,Nk,Nw,Norb)
  end if
  call get_chi_map(chi_map,irr_chi,olist,Nchi)
  call diis_init(diis,m_diis,Nk,Nw,Norb)
  allocate(eig_gkinv(Norb,Nk,Nw))
  ! state of get_rn (module mkself_rn)
  rn_Nk=Nk; rn_Nw=Nw; rn_Norb=Norb; rn_Nkall=Nkall; rn_temp=temp
  rn_eig=>eig; rn_kmult=>kmult; rn_eig_gkinv=>eig_gkinv
  bracket_valid=.false.
  muS_cache=0.0d0
  muL_cache=0.0d0
  iter_loop: do scf_i=1,scf_loop
     print'(A5,I5)','iter=',scf_i
     if(sw_tail)then
        ! tail-corrected bubble: conv[G]-conv[G0]+analytic reference at the CURRENT mu
        call chi0_tail_impl(chi,Gk,eig,uni,mu,kmap,invk,irr_chi,chi_map,olist,temp,&
             Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi)
     else
        call get_chi0_conv(chi,Gk,kmap,invk,irr_chi,chi_map,olist,temp,Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi)
     end if
     call ckchi_impl(chi,Smat,Cmat,kmap,invk,Nk,Nkall,Nchi,Nw,maxchi0s_global)
     if(sw_rescale .and. maxchi0s_global>=1.0d0)then
        ! the clamp target must stay just below 1: the resulting overshoot of the
        ! next Sigma is what kicks the loop off the rescaled manifold.  Any softer
        ! clamp (fixed 0.98, adaptive 1-(S-1), ...) admits a spurious self-consistent
        ! rescaled state and converges to it even when a genuine unrescaled fixed
        ! point with S < 1 exists (verified on a near-critical two-orbital model)
        rescale_t=1.0d0-1.0d-4
        print'(A,F10.6,A,F7.4)','[FLEX] Stoner factor=',maxchi0s_global,'>= 1: rescaling chi0 to',rescale_t
        chi(:,:,:,:)=chi(:,:,:,:)*rescale_t/maxchi0s_global
     end if
     call get_Vsigma_flex_nosoc(chi,Smat,Cmat,Nk,Nw,Nchi)
     print'(A16,E12.4,A5,E12.4)','Re V_sigma: max:',maxval(dble(chi)),' min:',minval(dble(chi))
     call calc_sigma(sigmak,Gk,chi,Smat,Cmat,kmap,invk,olist,temp,Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz)
     ! sub_sigma: subtract a static, k-DEPENDENT part of Sigma from every
     ! frequency.  This is DELIBERATE and goes beyond removing the Hartree-Fock
     ! double counting: pinning Sigma(k,iw_0) (mode 1) or its frequency average
     ! (mode 2) to zero keeps the Fermi-surface SHAPE at the non-interacting one,
     ! which is the intended behaviour here - do not "fix" it back to a pure
     ! constant subtraction.
     !
     ! It also disposes of the only O(1/Nw) Matsubara-truncation error in Sigma.
     ! V_sigma = V_inf + dV with dV ~ 1/nu^2, so the dV*G convolution truncates at
     ! O(1/Nw^2) like the tail-corrected chi0; the O(1/Nw) piece sits entirely in
     ! the V_inf (HF) term, which is k- and w-independent and therefore removed
     ! here (or absorbed by mu).  Hence no tail correction is needed for Sigma
     ! even when sw_tail is on for chi0.
     if(sub_sigma>0)then
        sub_self:block
          integer(c_int32_t) iw,l,m
          complex(c_double),allocatable,dimension(:,:,:):: sub_sigmak
          allocate(sub_sigmak(Nk,Norb,Norb))
          if(sub_sigma==1)then
             !$omp parallel
             do l=1,Norb
                do m=l,Norb
                   !$omp workshare
                   sub_sigmak(:,l,m)=(sigmak(:,1,l,m)+conjg(sigmak(:,1,m,l)))*0.5d0
                   !$omp end workshare
                   if(l.ne.m)sub_sigmak(:,m,l)=conjg(sub_sigmak(:,l,m))
                end do
             end do
             !$omp end parallel
          else  !sub_sigma==2: frequency-average HF subtraction
             !$omp parallel
             do l=1,Norb
                do m=l,Norb
                   !$omp workshare
                   sub_sigmak(:,l,m)=sum(sigmak(:,1:Nw,l,m)+conjg(sigmak(:,1:Nw,m,l)),dim=2)*(0.5d0/dble(Nw))
                   !$omp end workshare
                   if(l.ne.m)sub_sigmak(:,m,l)=conjg(sub_sigmak(:,l,m))
                end do
             end do
             !$omp end parallel
          end if
          !$omp parallel do private(iw)
          do iw=1,Nw
             sigmak(:,iw,:,:)=sigmak(:,iw,:,:)-sub_sigmak(:,:,:)
          end do
          !$omp end parallel do
        end block sub_self
     end if
     call compare_sigma()
     if(esterr<eps)then
        exit
     end if
     call renew_mu()
     call gen_green_inv(Gk,sigmak,hamk,mu,temp,Nk,Nw,Norb)
     call getinv(Gk,Nk,Nw,Norb)
     !$omp parallel
     !$omp workshare
     sigmak0(:,:,:,:)=sigmak(:,:,:,:)
     !$omp end workshare
     !$omp end parallel
  end do iter_loop
  if(maxchi0s_global>=1.0d0)then
     print'(A)','[FLEX] WARNING: Stoner factor >= 1 at the final iteration:'
     print'(A,F10.6,A)','[FLEX]   max eig[chi0*S] =',maxchi0s_global,&
          ' -> the normal state is magnetically unstable at this T/U.'
     if(sw_rescale)then
        print'(A)','[FLEX]   chi0 rescaling was still ACTIVE: the returned sigma/chi correspond to'
        print'(A)','[FLEX]   an artificially weakened interaction, NOT a physical FLEX solution.'
     else
        print'(A)','[FLEX]   chi_s is singular/negative there; the returned results are unreliable.'
     end if
  end if
  call diis_free(diis)
  if(sw_out)then
     call io_sigma(.true.)
  end if
  call renew_mu()
  nullify(rn_eig,rn_kmult,rn_eig_gkinv)
  deallocate(eig_gkinv)
contains
  subroutine compare_sigma()
    !> DIIS update of sigmak plus the convergence measure esterr (shared, see module flex_scf)
    call mix_sigma(diis,sigmak,sigmak0,esterr,pp,eps_sgm,invk,kmap,Nk,Nkall,Nw,Norb)
  end subroutine compare_sigma

  subroutine renew_mu()
    !> update mu so that get_rn(mu) = rfill (shared Brent solver, see module flex_scf)
    call precompute_eig_gkinv()   ! eigenvalues of H+Sigma, so get_rn needs no matrix inverse
    call solve_mu(mu,mu_OLD,bracket_valid,muS_cache,muL_cache,rfill,esterr,get_rn)
  end subroutine renew_mu

  subroutine precompute_eig_gkinv()
    use constant
    use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
    integer(c_int32_t) i,j,info_eig
    complex(c_double),dimension(Norb,Norb):: A_tmp
    complex(c_double),dimension(Norb):: eig_tmp
    complex(c_double),dimension(Norb,Norb):: vl_dum,vr_dum
    real(c_double),dimension(2*Norb):: rwork_eig
    complex(c_double),dimension(2*Norb):: work_eig

    !$omp parallel do private(A_tmp,eig_tmp,vl_dum,vr_dum,rwork_eig,work_eig,info_eig)
    do j=1,Nw
       do i=1,Nk
          A_tmp(:,:)=hamk(:,:,i)+sigmak(i,j,:,:)
          call zgeev('N','N',Norb,A_tmp,Norb,eig_tmp,vl_dum,Norb,vr_dum,Norb, &
                     work_eig,2*Norb,rwork_eig,info_eig)
          if(info_eig/=0)then
             print*,'zgeev failed in precompute_eig_gkinv: info=',info_eig
             stop
          end if
          eig_gkinv(:,i,j)=eig_tmp(:)
       end do
    end do
    !$omp end parallel do
  end subroutine precompute_eig_gkinv
  
  subroutine io_sigma(sw)
    logical,intent(in):: sw !True: out, False: in
    open(55,file='sigma.bin',form='unformatted')
    if(sw)then
       write(55)mu
       write(55)mu_OLD
       write(55)sigmak
       ! metadata footer (temp, Nw, Norb, Nk) for plibs.regrid_sigma_bin;
       ! trailing record, so readers of only the first 3 records are unaffected
       write(55)temp,dble(Nw),dble(Norb),dble(Nk)
    else
       read(55)mu
       read(55)mu_OLD
       read(55)sigmak
       call ck_sigma_meta()
    end if
    close(55)
  end subroutine io_sigma

  subroutine ck_sigma_meta()
    !> Verify the metadata footer (if present) against the current mesh.
    !> Without this check a sigma.bin with a LARGER Nw is read without any
    !> Fortran I/O error but with a silently scrambled layout.
    real(c_double) meta(4)
    integer(c_int32_t) ios
    read(55,iostat=ios)meta
    if(ios/=0)return  !no footer (old file): nothing to check
    if(nint(meta(2))/=Nw .or. nint(meta(3))/=Norb .or. nint(meta(4))/=Nk)then
       print'(A)','io_sigma ERROR: sigma.bin was written on a different mesh:'
       print'(A,I6,A,I4,A,I8)','   file: Nw=',nint(meta(2)),'  Norb=',nint(meta(3)),'  Nk=',nint(meta(4))
       print'(A,I6,A,I4,A,I8)','   run:  Nw=',int(Nw,c_int32_t),'  Norb=',int(Norb,c_int32_t),'  Nk=',int(Nk,c_int32_t)
       print'(A)','   convert the seed with plibs.regrid_sigma_bin first'
       stop
    end if
    if(abs(meta(1)-temp)>1.0d-10)then
       print'(A,F10.6,A,F10.6,A)','io_sigma: seed was written at T=',meta(1),&
            ' eV, reused raw at T=',temp,' eV (T-annealing)'
    end if
  end subroutine ck_sigma_meta
end subroutine mkself

subroutine get_rn_mkself(rn,rmu,eig,kmult,eig_gkinv,temp,Nk,Nw,Norb,Nkall)
  !> n(mu) of mkself; called through get_rn (module mkself_rn), which supplies mkself's state
  use constant
  use freq_grid, only:fg_uniform,fg_dens_coef,fg_omega,FG_FERMION
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,Norb,Nkall
  real(c_double),intent(in):: rmu,temp
  real(c_double),intent(in),dimension(Norb,Nk):: eig
  real(c_double),intent(in),dimension(Nk):: kmult
  complex(c_double),intent(in),dimension(Norb,Nk,Nw):: eig_gkinv
  real(c_double),intent(out):: rn

  integer(c_int32_t) i,j,n
  real(c_double) tmp,deltagk
  complex(c_double) iw_im,cp,cm,x
  integer(c_int64_t) j64

  ! n(μ) = (1/Nkall){Σ_{k,n} w_k f(ε_{kn}) + 2T·Σ_{k,j} w_k [G(k,iω_j) - G_0(k,iω_j)]}
  ! w_k = kmult(k) is the multiplicity of the irreducible k in the full BZ (1 on the TRIMs,
  ! 2 elsewhere); sum(w_k) = Nkall, so this is the exact full-BZ average.
  ! non-interacting Fermi sum
  tmp=0.0d0
  do i=1,Nk
     tmp=tmp+kmult(i)*sum(0.5d0*(1.0d0-tanh(0.5d0*(eig(:,i)-rmu)/temp)))
  end do
  ! correction from self-energy: use precomputed eigenvalues of H(k)+Sigma(k,iw): no getinv needed
  deltagk=0.0d0
  if(.not.fg_uniform())then
     ! IR: T sum_n [G - G_0](i w_n) is [G - G_0](tau = 0), continuous there (no 1/iw tail), read off
     ! the sampled values by the basis (fg_dens_coef); G - G_0 at -i w is the conjugate
     !$omp parallel do reduction(+:deltagk) private(n,i,iw_im,cp,cm,x,j64)
     do j=1,int(Nw,c_int32_t)
        j64=j
        iw_im=cmplx(0.0d0,fg_omega(FG_FERMION,j64,temp),kind=c_double)
        call fg_dens_coef(j64,cp,cm)
        do i=1,Nk
           do n=1,Norb
              x=1.0d0/(rmu+iw_im-eig_gkinv(n,i,j))-1.0d0/(rmu+iw_im-eig(n,i))
              deltagk=deltagk+kmult(i)*dble(cp*x+cm*conjg(x))
           end do
        end do
     end do
     !$omp end parallel do
     rn=(tmp+deltagk)/Nkall
     return
  end if
  !$omp parallel do reduction(+:deltagk) private(n,i,iw_im)
  do j=1,Nw
     iw_im=cmplx(0.0d0,dble(2*j-1)*pi*temp,kind=c_double)
     do i=1,Nk
        do n=1,Norb
           deltagk=deltagk+kmult(i)*dble(1.0d0/(rmu+iw_im-eig_gkinv(n,i,j)) &
                                         -1.0d0/(rmu+iw_im-eig(n,i)))
        end do
     end do
  end do
  !$omp end parallel do
  rn=(tmp+2*temp*deltagk)/Nkall
end subroutine get_rn_mkself

subroutine get_chis_chic(chis,chic,chi,Smat,Cmat,Nk,Nw,Nchi) bind(C)
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,Nchi
  real(c_double),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  complex(c_double),intent(in),dimension(Nk,Nw,Nchi,Nchi):: chi
  complex(c_double),intent(out),dimension(Nk,Nchi,Nchi):: chis,chic

  integer(c_int32_t) i,l,info
  integer(c_int32_t),dimension(Nchi):: ipiv
  complex(c_double),dimension(Nchi,Nchi):: cmat1,cmat2,cmat3,cmat4,Smat_c,Cmat_c

  Smat_c=cmplx(Smat,0.0d0,kind=c_double)
  Cmat_c=cmplx(Cmat,0.0d0,kind=c_double)

  qloop:do i=1,Nk
     call zgemm('N','N',Nchi,Nchi,Nchi,(-1.0d0,0.0d0),chi(i,1,:,:),Nchi,Smat_c,Nchi,(0.0d0,0.0d0),cmat1,Nchi) !-chi0S
     call zgemm('N','N',Nchi,Nchi,Nchi,(1.0d0,0.0d0),chi(i,1,:,:),Nchi,Cmat_c,Nchi,(0.0d0,0.0d0),cmat2,Nchi)  !chi0C
     do l=1,Nchi
        cmat1(l,l)=cmat1(l,l)+1.0d0 !I-chi0S
        cmat2(l,l)=cmat2(l,l)+1.0d0 !I+chi0C
     end do
     cmat3(:,:)=chi(i,1,:,:) !RHS copy for spin susceptibility
     cmat4(:,:)=chi(i,1,:,:) !RHS copy for charge susceptibility
     ! RPA spin susceptibility:   chi_s = (I - chi0·S)^{-1} chi0
     call zgesv(Nchi,Nchi,cmat1,Nchi,ipiv,cmat3,Nchi,info) !(I-chi0S)X=chi0 -> cmat3=chiS
     if(info/=0)then; print*,'zgesv failed: info=',info; stop; end if
     ! RPA charge susceptibility: chi_c = (I + chi0·C)^{-1} chi0
     call zgesv(Nchi,Nchi,cmat2,Nchi,ipiv,cmat4,Nchi,info) !(I+chi0C)X=chi0 -> cmat4=chiC
     if(info/=0)then; print*,'zgesv failed: info=',info; stop; end if
     chis(i,:,:)=cmat3(:,:)
     chic(i,:,:)=cmat4(:,:)
  end do qloop
end subroutine get_chis_chic

subroutine get_eig_or_tr_chi(chiqout,chi,invk,Nkall,Nk,Nchi,sw_eig) bind(C)
   use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
   implicit none
   integer(c_int64_t),intent(in):: Nkall,Nk,Nchi
   integer(c_int64_t),intent(in),dimension(3,Nkall):: invk
   logical,intent(in):: sw_eig
   complex(c_double),intent(in),dimension(Nk,Nchi,Nchi):: chi
   complex(c_double),intent(out),dimension(Nkall):: chiqout

   integer(c_int32_t) i,l,info
   real(c_double),dimension(2*Nchi):: rwork
   complex(c_double),dimension(Nchi*Nchi*4+1):: work
   complex(c_double),dimension(Nchi):: eig
   complex(c_double),dimension(Nchi,Nchi):: tmp,tmp1,tmp2

   !$omp parallel
   !$omp workshare
   chiqout(:)=0.0d0
   !$omp end workshare
   !$omp do private(tmp,tmp1,tmp2,rwork,work,eig,i,l,info)
   do i=1,Nkall
      if(sw_eig)then
         tmp(:,:)=chi(invk(1,i),:,:)
         call zgeev('N','N',Nchi,tmp,Nchi,eig,tmp1,Nchi,tmp2,Nchi,work,Nchi*Nchi*4+1,rwork,info)
         chiqout(i)=maxval(dble(eig))
      else
         do l=1,Nchi
            chiqout(i)=chiqout(i)+chi(invk(1,i),l,l)
         end do
      end if
   end do
   !$omp end do
   !$omp end parallel
end subroutine