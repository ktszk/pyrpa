module mkself_soc_rn
  !> Electron number n(mu) for mkself_soc's chemical-potential search (flex_scf's solve_mu).
  !> A module procedure reading state that mkself_soc points at, not an internal procedure of
  !> mkself_soc, for the reason given in module mkself_rn (fself.f90): no executable stack.
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double
  implicit none
  integer(c_int64_t):: rn_Nk,rn_Nw,rn_Norb,rn_Nkall
  real(c_double):: rn_temp
  real(c_double),pointer,contiguous:: rn_eig(:,:)=>null(),rn_kmult(:)=>null()
  complex(c_double),pointer,contiguous:: rn_Gk(:,:,:,:)=>null(),rn_sigmak(:,:,:,:)=>null()
  complex(c_double),pointer,contiguous:: rn_hamk(:,:,:)=>null(),rn_uni(:,:,:)=>null()
contains
  subroutine get_rn(rn,rmu)
    real(c_double),intent(in):: rmu
    real(c_double),intent(out):: rn
    call get_rn_mkself_soc(rn,rmu,rn_Gk,rn_sigmak,rn_hamk,rn_uni,rn_eig,rn_kmult,rn_temp,&
         rn_Nk,rn_Nw,rn_Norb,rn_Nkall)
  end subroutine get_rn
end module mkself_soc_rn

subroutine mkself_soc(sigmak,mu,Vmat,kmap,invk,invs,olist,slist,hamk,eig,uni,mu_init,&
     rfill,temp,scf_loop,pp,eps,Nkall,Nk,Nw,Norb,Nchi,Nx,Ny,Nz,sub_sigma,sw_out,&
     sw_in,m_diis,sw_rescale,sigma_scale) bind(C)
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t,c_bool
  use flex_scf, only:solve_mu,diis_state,diis_init,diis_free,mix_sigma
  use mkself_soc_rn
  use freq_grid, only:fg_uniform
  implicit none
  integer(c_int64_t),intent(in):: Nw,Norb,Nchi,Nkall,Nk,Nx,Ny,Nz,scf_loop,m_diis
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),dimension(Norb):: slist,invs
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
  logical(c_bool),intent(in):: sw_in,sw_out,sw_rescale
  !> The self-energy channel always subtracts the FLEX second-order double counting: the RPA
  !> matrix series reaches that diagram through two spin channels, and keeping both overestimates
  !> Sigma (by 2x for one orbital, 1.5x for two orbitals at J=0, measured against the exact second
  !> order).  See get_V_soc_flex; the pairing channel passes .false. instead, having no such term.
  logical(1),parameter:: dcsub=.true.
  integer(c_int64_t),intent(in):: sub_sigma
  real(c_double),intent(in):: temp,eps,pp,rfill,mu_init
  real(c_double),intent(in):: sigma_scale  !< multiply the sigma seed loaded via sw_in (U-annealing: ~(U_new/U_old)^2)
  real(c_double),intent(in),target,dimension(Norb,Nk):: eig
  real(c_double),intent(in),dimension(Nchi,Nchi):: Vmat
  real(c_double),intent(out):: mu
  complex(c_double),intent(in),target,dimension(Norb,Norb,Nk):: uni,hamk
  complex(c_double),intent(out),target,dimension(Nk,Nw,Norb,Norb):: sigmak

  integer(c_int32_t) scf_i
  integer(c_int32_t),dimension(Nchi,Nchi,2)::chi_map
  integer(c_int32_t),dimension(Nchi*(Nchi+1)/2,2)::irr_chi
  integer(c_int64_t),dimension(Nchi):: invschi
  real(c_double)esterr,mu_OLD,eps_sgm,maxchi0_global,rescale_t
  real(c_double),dimension(Norb,Norb):: sgnsig
  real(c_double),dimension(Nchi,Nchi):: sgnsig2
  complex(c_double),allocatable,target,dimension(:,:,:,:):: Gk
  complex(c_double),allocatable,dimension(:,:,:,:):: sigmak0
  complex(c_double),allocatable,dimension(:,:,:,:):: chi
  ! DIIS
  type(diis_state):: diis   !! Pulay history for the sigma mixing (module flex_scf)
  ! bracket cache for renew_mu
  real(c_double):: muS_cache,muL_cache
  logical:: bracket_valid
  ! multiplicity of each irreducible k in the full BZ (1 on the time-reversal invariant
  ! momenta, 2 elsewhere) -- see the note in fself.f90's mkself
  real(c_double),allocatable,target,dimension(:):: kmult
  integer(c_int32_t) ik_mult

  if(.not.fg_uniform())then
     ! sub_sigma=2 averages over a uniform grid, and sigma.bin is read back as one (plibs.regrid_sigma_bin)
     if(sub_sigma==2 .or. sw_in .or. sw_out)then
        print'(A)','ERROR(mkself_soc): sub_sigma=2 and sw_in/sw_out are for the uniform Matsubara grid only'
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
  call get_invschi_soc(invschi,olist,invs,Nchi)
  call get_sgnsig()
  call get_chi_map_soc(chi_map,irr_chi,olist,invs,Nchi,Norb)
  call diis_init(diis,m_diis,Nk,Nw,Norb)
  ! state of get_rn (module mkself_soc_rn)
  rn_Nk=Nk; rn_Nw=Nw; rn_Norb=Norb; rn_Nkall=Nkall; rn_temp=temp
  rn_Gk=>Gk; rn_sigmak=>sigmak; rn_hamk=>hamk; rn_uni=>uni; rn_eig=>eig; rn_kmult=>kmult
  bracket_valid=.false.
  muS_cache=0.0d0
  muL_cache=0.0d0
  iter_loop: do scf_i=1,scf_loop
     print'(A5,I5)','iter=',scf_i
     call get_chi0_conv_soc(chi,Gk,kmap,invk,invs,irr_chi,chi_map,olist,sgnsig,&
          sgnsig2,temp,Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi)
     call ckchi_soc_impl(chi,Vmat,kmap,invk,Nk,Nkall,Nchi,Nw,maxchi0_global)
     if(sw_rescale .and. maxchi0_global>=1.0d0)then
        ! the clamp target must stay just below 1: see the note in mkself (fself.f90)
        rescale_t=1.0d0-1.0d-4
        print'(A,F10.6,A,F7.4)','[FLEX] Stoner factor=',maxchi0_global,'>= 1: rescaling chi0 to',rescale_t
        chi(:,:,:,:)=chi(:,:,:,:)*rescale_t/maxchi0_global
     end if
     call get_V_soc_flex(chi,Vmat,olist,slist,Nk,Nw,Nchi,dcsub)
     print'(A16,E12.4,A5,E12.4)','Re V_sigma: max:',maxval(dble(chi)),' min:',minval(dble(chi))
     call calc_sigma_soc(sigmak,Gk,chi,Vmat,kmap,invk,invs,olist,slist,sgnsig,sgnsig2,temp,Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz)
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
  if(maxchi0_global>=1.0d0)then
     print'(A)','[FLEX] WARNING: Stoner factor >= 1 at the final iteration:'
     print'(A,F10.6,A)','[FLEX]   max eig[chi0*V] =',maxchi0_global,&
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
  nullify(rn_Gk,rn_sigmak,rn_hamk,rn_uni,rn_eig,rn_kmult)
contains
  subroutine get_sgnsig()
    integer(c_int32_t) i,j
    do j=1,Norb
       do i=j,Norb
          sgnsig(i,j)=slist(i)*slist(j)
          sgnsig(j,i)=sgnsig(i,j)
       end do
    end do

    do j=1,Nchi
       do i=j,Nchi
          sgnsig2(i,j)=slist(olist(i,1))*slist(olist(i,2))&
               *slist(olist(j,1))*slist(olist(j,2))
          sgnsig2(j,i)=sgnsig2(i,j)
       end do
    end do
  end subroutine get_sgnsig



  subroutine compare_sigma()
    !> DIIS update of sigmak plus the convergence measure esterr (shared, see module flex_scf)
    call mix_sigma(diis,sigmak,sigmak0,esterr,pp,eps_sgm,invk,kmap,Nk,Nkall,Nw,Norb)
  end subroutine compare_sigma

  subroutine renew_mu()
    !> update mu so that get_rn(mu) = rfill (shared Brent solver, see module flex_scf)
    call solve_mu(mu,mu_OLD,bracket_valid,muS_cache,muL_cache,rfill,esterr,get_rn)
  end subroutine renew_mu
  
  subroutine io_sigma(sw)
    !> Same single-record layout as the no-SOC io_sigma (fself.f90):
    !> mu, mu_OLD, sigmak, metadata footer (temp, Nw, Norb, Nk).
    !> The old element-wise SOC format is no longer read.
    logical,intent(in):: sw !True: out, False: in
    open(55,file='sigma.bin',form='unformatted')
    if(sw)then
       write(55)mu
       write(55)mu_OLD
       write(55)sigmak
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
    !> Verify the metadata footer (if present) against the current mesh;
    !> see the note in fself.f90
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
end subroutine mkself_soc

subroutine get_rn_mkself_soc(rn,rmu,Gk,sigmak,hamk,uni,eig,kmult,temp,Nk,Nw,Norb,Nkall)
  !> n(mu) of mkself_soc; called through get_rn (module mkself_soc_rn), which supplies the
  !> state.  Gk is mkself_soc's work array and is overwritten.
  use constant
  use freq_grid, only:fg_uniform,fg_dens_coef,fg_omega,FG_FERMION
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,Norb,Nkall
  real(c_double),intent(in):: rmu,temp
  real(c_double),intent(in),dimension(Norb,Nk):: eig
  real(c_double),intent(in),dimension(Nk):: kmult
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni,hamk
  complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: sigmak
  complex(c_double),intent(inout),dimension(Nk,Nw,Norb,Norb):: Gk
  real(c_double),intent(out):: rn

  integer(c_int32_t) l,i,j,n
  real(c_double) tmp,deltagk
  complex(c_double):: Gk0,iw,cp,cm,x
  integer(c_int64_t) j64

  ! weight each irreducible k by its full-BZ multiplicity; sum(kmult) = Nkall
  tmp=0.0d0
  do i=1,Nk
     tmp=tmp+kmult(i)*sum(0.5d0*(1.0d0-tanh(0.5d0*(eig(:,i)-rmu)/temp)))
  end do
  call gen_green_inv(Gk,sigmak,hamk,rmu,temp,Nk,Nw,Norb)
  call getinv(Gk,Nk,Nw,Norb)
  deltagk=0.0d0
  ! IR: T sum_n [G - G_0](i w_n) is [G - G_0](tau = 0) from the sampled values (see get_rn_mkself);
  ! uniform: 2T Re sum_n.  G_ll(-i w) = G_ll(i w)^*.
  !$omp parallel do reduction(+:deltagk) private(l,i,n,iw,Gk0,cp,cm,x,j64)
  do j=1,int(Nw,c_int32_t)
     j64=j
     if(.not.fg_uniform())call fg_dens_coef(j64,cp,cm)
     iw=cmplx(rmu,fg_omega(FG_FERMION,j64,temp),kind=c_double)
     do l=1,Norb
        do i=1,Nk
           Gk0=0.0d0
           do n=1,Norb
              Gk0=Gk0+uni(l,n,i)*conjg(uni(l,n,i))/(iw-eig(n,i))
           end do
           x=Gk(i,j,l,l)-Gk0
           if(fg_uniform())then
              deltagk=deltagk+kmult(i)*dble(x)
           else
              deltagk=deltagk+kmult(i)*dble(cp*x+cm*conjg(x))
           end if
        end do
     end do
  end do
  !$omp end parallel do
  if(fg_uniform())deltagk=2*temp*deltagk
  rn=(tmp+deltagk)/Nkall
end subroutine get_rn_mkself_soc

subroutine calc_sigma_soc(sigmak,Gk,Vsigma,Vmat,kmap,invk,invs,olist,slist,sgnsig,sgnsig2,temp,Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz)
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
  integer(c_int64_t),intent(in),dimension(Norb):: slist,invs
  real(c_double),intent(in):: temp
  real(c_double),intent(in),dimension(Nchi,Nchi):: Vmat
  real(c_double),dimension(Norb,Norb):: sgnsig
  real(c_double),dimension(Nchi,Nchi):: sgnsig2
  complex(c_double),intent(in),dimension(Nk,Nw,Nchi,Nchi):: Vsigma
  complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(c_double),intent(out),dimension(Nk,Nw,Norb,Norb):: sigmak

  integer(c_int32_t) i,j,k,n,l,m
  integer(c_int64_t) jp,jm,j0
  real(c_double) weight
  complex(c_double) vinf,st
  complex(c_double),allocatable,dimension(:,:,:,:):: tmpVsigma,tmp,tmpgk

  allocate(tmpVsigma(0:Nx-1,0:Ny-1,0:Nz-1,fg_nt(Nw)),tmp(0:Nx-1,0:Ny-1,0:Nz-1,fg_nt(Nw)),&
       tmpgk(0:Nx-1,0:Ny-1,0:Nz-1,fg_nt(Nw)))
  j0=fg_pos(FG_BOSON,1_c_int64_t)   ! v=0
  weight=temp/dble(Nkall)
  !$omp parallel workshare
  sigmak(:,:,:,:)=0.0d0
  !$omp end parallel workshare
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
                 tmpgk(kmap(1,i),kmap(2,i),kmap(3,i),jp)=sgnsig(olist(l,2),olist(m,2))*Gk(invk(1,i),j,invs(olist(l,2)),invs(olist(m,2))) !G42(-k,iw)=G24(k,iw)
                 tmpgk(kmap(1,i),kmap(2,i),kmap(3,i),jm)=sgnsig(olist(l,2),olist(m,2))*conjg(Gk(invk(1,i),j,invs(olist(m,2)),invs(olist(l,2)))) !G42(-k,-iw)=G^*42(k,iw)
              end if
           end do
        end do
        !$omp end do
        !$omp end parallel
        vinf=cmplx(-Vmat(m,l),0.0d0,kind=c_double)
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

  !$omp parallel workshare
  sigmak(:,:,:,:)=sigmak(:,:,:,:)*weight
  !$omp end parallel workshare
  call drop_roundoff(sigmak,Nk*Nw*Norb*Norb)
end subroutine calc_sigma_soc
