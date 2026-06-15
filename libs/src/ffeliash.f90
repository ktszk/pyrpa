!non-linear FLEX-Eliashberg equations solver (without SOC)
subroutine eliashberg(delta,sigmak,Gk,hamk,Smat,Cmat,olist,prt,kmap,invk,mu,temp,eps,&
     Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,itemax,gap_sym,sw_sigma,sw_Vconst,m_diis,gap_min,&
     sw_amp_newton) bind(C)
  !> Non-linear self-consistent FLEX-Eliashberg loop (no SOC).
  !>
  !> One iteration:
  !>   1. Δ_new = (T/Nkall) · K[F]                          (gap equation, mkdelta_nsoc)
  !>   2. Δ ← Pulay/DIIS extrapolation of past (Δ_in, Δ_out) pairs    (or linear mixing fallback)
  !>   3. (sw_sigma=true) Σ ← FLEX self-energy from V_σ·G    (calc_sigma)
  !>      G_0_dressed = (G_0^-1 - Σ)^-1                       (gen_green_inv + getinv)
  !>   4. G = G_0_dressed (I + Δ·F†)                          (mkgreliah_trs_nsoc, Dyson)
  !>      F = -G_0_dressed · Δ · G^†                          (mkfkeliash_trs_nsoc)
  !>   5. χ_0_S = χ^GG + χ^FF†,  χ_0_C = χ^GG - χ^FF†         (get_chi0sc)
  !>   6. V_σ, V_Δ ← FLEX vertex from χ_0_S, χ_0_C            (mkV_flex_nosoc)
  !>
  !> The caller is responsible for supplying a non-zero initial delta (e.g. via
  !> get_initial_delta or a previously-converged linearized-Eliashberg gap function).
  !>
  !> DIIS (Pulay) update for delta:
  !>   keep circular buffer of recent (Δ_in, Δ_out) pairs of length m_diis,
  !>   residual r_i = Δ_out_i - Δ_in_i,
  !>   solve [B 1; 1ᵀ 0][c; λ] = [0; 1] with B_ij = Re<r_i, r_j>,
  !>   set Δ ← Σ_i c_i · Δ_out_i  (true Pulay extrapolation when n_cur ≥ 2;
  !>   linear mixing pp·Δ_out + (1-pp)·Δ_in otherwise).
  !>
  !!@param      delta,inout: gap function Δ(k,iω); user-supplied initial value, overwritten with converged Δ
  !!@param     sigmak,out: self-energy Σ(k,iω); zeroed when sw_sigma=false
  !!@param       hamk,in: normal-state Hamiltonian H(k)
  !!@param       Smat,in: S-matrix (Stoner)
  !!@param       Cmat,in: C-matrix (charge)
  !!@param      olist,in: chi-index orbital pairs
  !!@param        prt,in: parity factors per orbital
  !!@param       kmap,in: k-point property
  !!@param       invk,in: inverse-k mapping
  !!@param         mu,in: chemical potential
  !!@param       temp,in: temperature
  !!@param        eps,in: convergence threshold on |Δ_new-Δ|/|Δ|
  !!@param     itemax,in: maximum SCF iterations
  !!@param    gap_sym,in: gap symmetry (>=0 singlet, <0 triplet)
  !!@param   sw_sigma,in: include FLEX self-energy correction
  !!@param   sw_Vconst,in: use constant pairing interaction
  !!@param     m_diis,in: DIIS history depth (m_diis<=1 → linear mixing only)
  !!@param    gap_min,in: gap-collapse threshold as a fraction of the initial |Δ|max.
  !!                      When |Δ_new|max drops below gap_min·|Δ_init|max the gap is
  !!                      decaying toward the trivial solution (T >= Tc): stop early.
  !!                      Set gap_min<=0 to disable the check.
  !!@param sw_amp_newton,in: enable amplitude-direction Newton (secant) acceleration.
  !!                      Near the bifurcated SC fixed point the slow mode is the gap
  !!                      amplitude a; the shape is relaxed by plain fixed-point
  !!                      (u=M[Δ]/‖M[Δ]‖) while a is driven by a secant root-find on the
  !!                      Rayleigh gain g(a)=<Δ,M[Δ]>/<Δ,Δ>=1. Bypasses DIIS when true.
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t,c_bool
  implicit none
  integer(c_int64_t),intent(in):: Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,itemax,gap_sym,m_diis
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
  logical(c_bool),intent(in):: sw_sigma,sw_Vconst,sw_amp_newton
  real(c_double),intent(in):: temp,eps,mu,gap_min
  real(c_double),intent(in),dimension(Norb):: prt
  real(c_double),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: hamk
  complex(c_double),intent(inout),dimension(Nk,Nw,Norb,Norb):: Gk, sigmak, delta

  integer(c_int32_t) i_iter
  integer(c_int32_t),dimension(Nchi,Nchi,2)::chi_map
  integer(c_int32_t),dimension(Nchi*(Nchi+1)/2,2)::irr_chi
  logical(1) sw_pair
  logical do_pulay
  real(c_double) weight,maxchi0s_global,delta_diff,delta_norm,delta0_max,new_max
  ! amplitude accelerator state (v0-anchored cubic-pitchfork solve, Phase 1)
  integer(c_int32_t) amp_phase,npts,n_polish,n_anchor
  real(c_double) a_cur,a_nxt,g_cur,a1,a2,g1,g2,proj,cfit,lam0
  real(c_double),parameter:: amp_tol=2.0d-2   ! |g-1| below which the amplitude is "solved"
  integer(c_int32_t),parameter:: n_polish_max=10   ! DIIS shape-polish steps between re-anchors
  real(c_double),parameter:: pp=0.3d0,eps_reg=1.0d-12   ! linear mixing rate / DIIS regularization
  ! DIIS state
  integer(c_int32_t) m_diis_eff,n_hist,i_hist,n_cur,ih,jh,idx_i,idx_j,info_diis
  complex(c_double),allocatable,dimension(:,:,:,:,:):: xout_hist,res_hist
  real(c_double),allocatable,dimension(:,:):: B_diis
  real(c_double),allocatable,dimension(:):: rhs_diis
  integer(c_int32_t),allocatable,dimension(:):: ipiv_diis
  complex(c_double),dimension(Nk,Nw,Norb,Norb):: newdelta,fk,Gk0,v0
  complex(c_double),allocatable,dimension(:,:,:,:):: Vdelta,Vsigma

  ! All FLEX vertex / chi0 routines expect (Nk,Nw,Nchi,Nchi) ordering.
  allocate(Vdelta(Nk,Nw,Nchi,Nchi),Vsigma(Nk,Nw,Nchi,Nchi))
  ! DIIS: history length is at least 1 (slot for current iterate, used for linear-mixing fallback).
  m_diis_eff=max(int(m_diis,c_int32_t),1)
  allocate(xout_hist(Nk,Nw,Norb,Norb,m_diis_eff))
  allocate(res_hist (Nk,Nw,Norb,Norb,m_diis_eff))
  allocate(B_diis(m_diis_eff+1,m_diis_eff+1))
  allocate(rhs_diis(m_diis_eff+1))
  allocate(ipiv_diis(m_diis_eff+1))
  n_hist=0
  i_hist=0
  ! amplitude accelerator: anchor the gap shape to the (normalized) seed eigenvector v0.
  amp_phase=1
  npts=0
  n_polish=0
  n_anchor=0
  a1=0.0d0; a2=0.0d0; g1=0.0d0; g2=0.0d0
  a_cur=sqrt(real(sum(conjg(delta)*delta),kind=c_double))+1.0d-300  ! seed amplitude ‖Δ_init‖
  v0(:,:,:,:)=delta(:,:,:,:)/a_cur                                  ! unit-norm gap shape
  if(.not.sw_amp_newton) amp_phase=2   ! skip Phase 1; use DIIS throughout

  if(gap_sym>=0)then
     sw_pair=.true.
     print'(A7)','singlet'
  else
     sw_pair=.false.
     print'(A7)','triplet'
  end if
  weight=temp/dble(Nkall)
  if(m_diis_eff<2)then
     print'(A,I0,A,F4.2,A)','DIIS disabled (m_diis=',m_diis,'); pure linear mixing pp=',pp,'.'
  else
     print'(A,I0,A,F4.2,A)','DIIS enabled, m_diis=',m_diis_eff,'; linear-mixing fallback pp=',pp,'.'
  end if

  ! ===== initial setup ====================================================
  ! NOTE: caller must supply a non-zero initial delta. We DO NOT call
  ! get_initial_delta here — using e.g. a converged linear-Eliashberg gap as
  ! initial value is the recommended high-accuracy fast-start.
  Gk0(:,:,:,:)=Gk(:,:,:,:) ! Gk0 starts as bare G_0, may be updated with Σ if sw_sigma
  delta0_max=maxval(abs(delta))+1.0d-30 ! seed |Δ|max, reference for gap-collapse (T>=Tc) check
  call get_chi_map(chi_map,irr_chi,olist,Nchi)                          ! chi-index symmetry table for FFT loops
  if(sw_Vconst)then
     ! With constant pairing interaction, we only need to compute χ_0_S, χ_0_C once at the very beginning.
     fk(:,:,:,:)=0.0d0
     call get_chi0sc(Vsigma,Vdelta,Gk,fk,kmap,invk,irr_chi,chi_map,olist,prt,temp,&
                     Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi,sw_pair)
     call mkfk_trs_nsoc(fk,Gk,delta,Nk,Nw,Norb) ! linearized seed F = -G·Δ·G^†
  else
     ! With FLEX pairing interaction, we need χ_0_S, χ_0_C from the start to build V_σ, V_Δ for the first Σ and Δ update.
     ! This is because the FLEX vertices depend on χ_0, and the self-energy depends on the vertices.
     call mkfk_trs_nsoc(fk,Gk,delta,Nk,Nw,Norb) ! linearized seed F = -G·Δ·G^†
     call get_chi0sc(Vsigma,Vdelta,Gk,fk,kmap,invk,irr_chi,chi_map,olist,prt,temp,&
                     Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi,sw_pair)                            ! Vsigma=χ_0_S, Vdelta=χ_0_C
  end if
  call ckchi_impl(Vsigma,Smat,Cmat,kmap,invk,Nk,Nkall,Nchi,Nw,maxchi0s_global)   ! SDW Stoner check on χ_0_S
  call mkV_flex_nosoc(Vdelta,Vsigma,Smat,Cmat,Nk,Nw,Nchi,sw_pair)        ! χ_0 → V_σ (Vsigma), V_Δ (Vdelta)
  if(.not.sw_sigma)then
     !$omp parallel workshare
     sigmak(:,:,:,:)=(0.0d0,0.0d0)
     !$omp end parallel workshare
  end if
  ! ===== self-consistent loop ============================================
  iter_loop: do i_iter=1,itemax
     print'(A5,I5)','iter=',i_iter
     ! 1. Δ_new = K[F] (gap equation)
     call mkdelta_nsoc(newdelta,fk,Vdelta,Smat,Cmat,kmap,invk,prt,olist,&
          Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,sw_pair)
     !$omp parallel workshare
     newdelta(:,:,:,:)=newdelta(:,:,:,:)*weight                          ! T/Nkall normalization
     !$omp end parallel workshare
     ! 2. convergence check (raw residual, before mixing)
     delta_diff=maxval(abs(newdelta-delta))
     delta_norm=maxval(abs(delta))+1.0d-30
     new_max=maxval(abs(newdelta))
     print'(A,E12.4,A,E12.4,A,E12.4)','  |Δ_new-Δ|max=',delta_diff,'  rel=',delta_diff/delta_norm,&
          '  |Δ_new|max=',new_max
     if(delta_diff/delta_norm<eps)then
        !$omp parallel workshare
        delta(:,:,:,:)=newdelta(:,:,:,:)
        !$omp end parallel workshare
        print'(A,I5)','converged at iter=',i_iter
        exit iter_loop
     end if
     ! 2b. gap-collapse check: Δ decaying toward the trivial solution (T >= Tc)
     if(gap_min>0.0d0 .and. new_max<gap_min*delta0_max)then
        !$omp parallel workshare
        delta(:,:,:,:)=newdelta(:,:,:,:)
        !$omp end parallel workshare
        print'(A,E12.4,A,E12.4,A)','gap collapsed: |Δ_new|max=',new_max,' < ',gap_min*delta0_max,&
             ' (likely T >= Tc; no SC solution), stopping'
        exit iter_loop
     end if

     ! 3. adaptive re-anchor: after n_polish_max DIIS shape-polish steps, adopt the
     !    relaxed shape as the new anchor v0 and re-solve the amplitude. This re-locks
     !    the amplitude (which otherwise wobbles during the DIIS polish) and projects out
     !    any negative-mode contamination, while v0 keeps improving toward the true shape.
     if(sw_amp_newton .and. amp_phase==2 .and. n_polish>=n_polish_max)then
        a_cur=sqrt(real(sum(conjg(delta)*delta),kind=c_double))+1.0d-300
        !$omp parallel workshare
        v0(:,:,:,:)=delta(:,:,:,:)/a_cur                   ! re-anchor to current relaxed shape
        !$omp end parallel workshare
        npts=0; a1=0.0d0; a2=0.0d0; g1=0.0d0; g2=0.0d0
        amp_phase=1
        n_anchor=n_anchor+1
        print'(A,I0,A)','  re-anchor #',n_anchor,': v0 <- relaxed shape, re-solving amplitude'
     end if
     ! 3'. update of delta
     if(amp_phase==1)then
        ! ---- Phase 1: v0-anchored cubic-pitchfork amplitude solve ---------------
        ! The gap shape is held at the linearized eigenvector v0 (Δ=a·v0), which kills
        ! the drift toward the unphysical negative eigenmode. Only the scalar amplitude
        ! a is solved. The v0-projected gain g(a)=<v0,M[a v0]>/a follows the pitchfork
        ! form g(a)=λ0-c·a^2 near threshold; the SC amplitude is the root of g(a)=1,
        ! a*=sqrt((λ0-1)/c). λ0 and c are fit from two probed amplitudes, then a jumps
        ! to a*. Once |g-1|<amp_tol the amplitude is solved and we switch to DIIS
        ! (Phase 2) to relax the shape (the v0-perpendicular component) to convergence.
        ! newdelta = M[a_cur·v0]; current iterate Δ = a_cur·v0, so <v0,Δ>=a_cur.
        proj =real(sum(conjg(v0)*newdelta),kind=c_double)
        g_cur=proj/a_cur                                  ! gain g(a_cur) along v0
        a2=a1; g2=g1; a1=a_cur; g1=g_cur; npts=min(npts+1,2)
        if(npts>=2 .and. abs(a1*a1-a2*a2)>1.0d-30)then
           cfit=(g1-g2)/(a2*a2-a1*a1)                      ! g=λ0-c·a^2  (slope in a^2)
           lam0=g1+cfit*a1*a1
           if(cfit>1.0d-12 .and. lam0>1.0d0)then
              a_nxt=sqrt((lam0-1.0d0)/cfit)                ! cubic-pitchfork root g(a*)=1
              a_nxt=min(max(a_nxt,0.3d0*a_cur),3.0d0*a_cur)
              print'(A,F9.4,A,F9.4,A,F9.4,A,E11.4,A,E11.4)','  amp[P1 fit]: g=',g_cur,&
                   ' lam0=',lam0,' c=',cfit,'  a=',a_cur,' -> ',a_nxt
           else
              a_nxt=0.5d0*a_cur                            ! no positive root (T>=Tc): shrink
              print'(A,F9.4,A,E11.4,A,E11.4)','  amp[P1 shrink]: g=',g_cur,'  a=',a_cur,' -> ',a_nxt
           end if
        else
           if(abs(g_cur-1.0d0)<0.2d0)then
              a_nxt=0.9d0*a_cur                            ! gentle probe near solution (re-anchor)
           else
              a_nxt=0.6d0*a_cur                            ! probe a second amplitude
           end if
           print'(A,F9.4,A,E11.4,A,E11.4)','  amp[P1 probe]: g=',g_cur,'  a=',a_cur,' -> ',a_nxt
        end if
        if(abs(g_cur-1.0d0)<amp_tol)then
           amp_phase=2                                     ! amplitude solved -> DIIS shape polish
           a_nxt=a_cur
           n_polish=0
           n_hist=0; i_hist=0                              ! fresh DIIS history on the (re)anchored shape
           print'(A,F9.5,A)','  amplitude solved (g=',g_cur,'); DIIS shape polish'
        end if
        !$omp parallel workshare
        delta(:,:,:,:)=a_nxt*v0(:,:,:,:)                   ! force shape=v0, amplitude=a_nxt
        !$omp end parallel workshare
        a_cur=a_nxt
        go to 100   ! skip DIIS block this iteration
     end if
     ! 3'. DIIS / linear mixing update of delta (Phase 2, or sw_amp_newton=False)
     i_hist=mod(i_hist,m_diis_eff)+1
     n_hist=min(n_hist+1,m_diis_eff)
     n_cur=n_hist
     !$omp parallel workshare
     xout_hist(:,:,:,:,i_hist)=newdelta(:,:,:,:)
     res_hist (:,:,:,:,i_hist)=newdelta(:,:,:,:)-delta(:,:,:,:)
     !$omp end parallel workshare
     do_pulay=(n_cur>=2 .and. m_diis_eff>=2)
     if(do_pulay)then
        ! ---- Pulay matrix: B_ij = Re<r_i, r_j>, augmented with Lagrange row/col ----
        B_diis(1:n_cur+1,1:n_cur+1)=0.0d0
        do ih=1,n_cur
           idx_i=modulo(i_hist-n_cur+ih-1,m_diis_eff)+1
           do jh=1,n_cur
              idx_j=modulo(i_hist-n_cur+jh-1,m_diis_eff)+1
              B_diis(ih,jh)=real(sum(conjg(res_hist(:,:,:,:,idx_i))*res_hist(:,:,:,:,idx_j)),kind=c_double)
           end do
           B_diis(ih,ih)=B_diis(ih,ih)+eps_reg                            ! regularize diagonal
           B_diis(ih,n_cur+1)=1.0d0                                       ! Lagrange constraint
           B_diis(n_cur+1,ih)=1.0d0
        end do
        rhs_diis(1:n_cur)=0.0d0
        rhs_diis(n_cur+1)=1.0d0                                           ! Σ c_i = 1
        call dgesv(n_cur+1,1,B_diis(1:n_cur+1,1:n_cur+1),n_cur+1, &
             ipiv_diis(1:n_cur+1),rhs_diis(1:n_cur+1),n_cur+1,info_diis)
        if(info_diis/=0)then
           print'(A,I0,A,F4.2,A)','DIIS: dgesv failed (info=',info_diis,'), falling back to linear mixing pp=',pp,'.'
           do_pulay=.false.
        end if
     end if
     if(do_pulay)then
        ! ---- Pulay extrapolation: Δ = Σ_i c_i · Δ_out_i ----
        !$omp parallel workshare
        delta(:,:,:,:)=(0.0d0,0.0d0)
        !$omp end parallel workshare
        do ih=1,n_cur
           idx_i=modulo(i_hist-n_cur+ih-1,m_diis_eff)+1
           !$omp parallel workshare
           delta(:,:,:,:)=delta(:,:,:,:)+rhs_diis(ih)*xout_hist(:,:,:,:,idx_i)
           !$omp end parallel workshare
        end do
     else
        ! linear mixing (n_cur=1 first iteration, DIIS disabled, or dgesv-failed fallback)
        !$omp parallel workshare
        delta(:,:,:,:)=pp*newdelta(:,:,:,:)+(1.0d0-pp)*delta(:,:,:,:)
        !$omp end parallel workshare
     end if
     n_polish=n_polish+1   ! count DIIS shape-polish steps since last (re)anchor

100  continue   ! resume point for amplitude-Newton path (skips DIIS block)
     ! 4. update Σ and dressed G_0 (FLEX self-energy)
     if(sw_sigma)then
        call calc_sigma(sigmak,Gk,Vsigma,Smat,Cmat,kmap,invk,olist,temp,&
             Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz)
        call gen_green_inv(Gk0,sigmak,hamk,mu,temp,Nk,Nw,Norb)             ! Gk0 := G_0^-1 - Σ
        call getinv(Gk0,Nk,Nw,Norb)                                        ! invert in place
     end if
     ! 5. SC Dyson: G = Gk0(I+Δ·F†), F = -Gk0·Δ·G^†
     call mkgreliah_trs_nsoc(Gk,Gk0,fk,delta,Nk,Nw,Norb)
     call mkfkeliash_trs_nsoc(fk,Gk0,Gk,delta,Nk,Nw,Norb)
     if(.not. sw_Vconst)then
        ! 6. χ_0_S, χ_0_C from new G, F
        call get_chi0sc(Vsigma,Vdelta,Gk,fk,kmap,invk,irr_chi,chi_map,olist,prt,temp,&
                        Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi,sw_pair)
        call ckchi_impl(Vsigma,Smat,Cmat,kmap,invk,Nk,Nkall,Nchi,Nw,maxchi0s_global)
        ! 7. update FLEX vertices V_σ (Vsigma), V_Δ (Vdelta)
        call mkV_flex_nosoc(Vdelta,Vsigma,Smat,Cmat,Nk,Nw,Nchi,sw_pair)
     end if
  end do iter_loop

  deallocate(xout_hist,res_hist,B_diis,rhs_diis,ipiv_diis)
  deallocate(Vdelta,Vsigma)
end subroutine eliashberg

subroutine mkgreliah_trs_nsoc(Gk,Gk0,fk,delta,Nk,Nw,Norb)
  !> Update normal Green's function in the SC state (TRS, no SOC):
  !>   G(k,iω) = Gk0(k,iω) · ( I + Δ(k,iω) · F†(k,iω) )
  !> derived from the Nambu Dyson equation closed under TRS.
  !> Gk0 is the (possibly Σ-dressed) normal Green's function.
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,Norb
  complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk0,fk,delta
  complex(c_double),intent(out),dimension(Nk,Nw,Norb,Norb):: Gk

  integer(c_int32_t) i,j,l,m
  complex(c_double),dimension(Norb,Norb):: cmat1

  !$omp parallel do collapse(2) private(i,j,l,m,cmat1)
  do j=1,Nw
     do i=1,Nk
        ! cmat1 = Δ(k,iω) · F†(k,iω)   (zgemm 'N','C' applies Hermitian conjugate to fk)
        call zgemm('N','C',Norb,Norb,Norb,(1.0d0,0.0d0),delta(i,j,:,:),Norb,fk(i,j,:,:),Norb,(0.0d0,0.0d0),cmat1,Norb)
        do l=1,Norb
           cmat1(l,l)=cmat1(l,l)+1.0d0   ! cmat1 := I + Δ·F†
        end do
        ! G = Gk0 · (I + Δ·F†)
        call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),Gk0(i,j,:,:),Norb,cmat1,Norb,(0.0d0,0.0d0),Gk(i,j,:,:),Norb)
     end do
  end do
  !$omp end parallel do
end subroutine mkgreliah_trs_nsoc

subroutine mkfkeliash_trs_nsoc(fk,Gk0,Gk,delta,Nk,Nw,Norb)
  !> Update anomalous Green's function in the SC state (TRS, no SOC):
  !>   F(k,iω) = -Gk0(k,iω) · Δ(k,iω) · G^†(k,iω)
  !> Gk should be the up-to-date normal G(k,iω) returned by mkgreliah_trs_nsoc.
  !> The diagonal real / off-diagonal Hermitian projection at the end enforces the
  !> Hermitian symmetry of F that holds under TRS without SOC.
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,Norb
  complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk,Gk0,delta
  complex(c_double),intent(out),dimension(Nk,Nw,Norb,Norb):: fk

  integer(c_int32_t) i,j,l,m
  complex(c_double),dimension(Norb,Norb):: cmat1

  !$omp parallel do collapse(2) private(i,j,l,m,cmat1)
  do j=1,Nw
     do i=1,Nk
        ! cmat1 = Δ(k,iω) · G†(k,iω)
        call zgemm('N','C',Norb,Norb,Norb,(1.0d0,0.0d0),delta(i,j,:,:),Norb,Gk(i,j,:,:),Norb,(0.0d0,0.0d0),cmat1,Norb)
        ! fk = -Gk0 · Δ · G†
        call zgemm('N','N',Norb,Norb,Norb,(-1.0d0,0.0d0),Gk0(i,j,:,:),Norb,cmat1,Norb,(0.0d0,0.0d0),fk(i,j,:,:),Norb)
        ! enforce F^† = F (TRS+no-SOC symmetry)
        do l=1,Norb
           fk(i,j,l,l)=dble(fk(i,j,l,l))
           do m=l+1,Norb
              fk(i,j,l,m)=conjg(fk(i,j,m,l))
           end do
        end do
     end do
  end do
  !$omp end parallel do
end subroutine mkfkeliash_trs_nsoc

subroutine mkV_flex_nosoc(Vdelta,Vsigma,Smat,Cmat,Nk,Nw,Nchi,sw_pair)
  !> Combined FLEX vertex builder for SC iteration.
  !> Solves the RPA series for χ_s, χ_c per (k,iω), then assembles V_σ and V_Δ.
  !>   χ_s = (I - χ_0_S·S)^-1 · χ_0_S·S
  !>   χ_c = (I + χ_0_C·C)^-1 · χ_0_C·C
  !>   V_σ = 3/2·S - 1/2·C + 3/2·S·χ_s + 1/2·C·χ_c - 1/4·(S+C)·(χ_s+χ_c)
  !>         (last term subtracts the double-counted second-order contribution)
  !>   V_Δ_singlet = (S+C)/2 + 3/2·S·χ_s - 1/2·C·χ_c
  !>   V_Δ_triplet = (S-C)/2 - 1/2·S·χ_s - 1/2·C·χ_c
  !>
  !> I/O semantics (in-place):
  !>   on entry  Vsigma = χ_0_S, Vdelta = χ_0_C
  !>   on exit   Vsigma = V_σ,   Vdelta = V_Δ
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,Nchi
  logical(1),intent(in):: sw_pair
  real(c_double),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  ! Storage convention shared with all FLEX routines: (Nk, Nw, Nchi, Nchi).
  complex(c_double),intent(inout),dimension(Nk,Nw,Nchi,Nchi):: Vdelta,Vsigma

  integer(c_int32_t) i,j,l,info
  integer(c_int32_t),dimension(Nchi):: ipiv
  complex(c_double),dimension(Nchi,Nchi):: cmat1,cmat2,cmat3,cmat4,cmat5,Smat_c,Cmat_c,V0_c,SC_c

  Smat_c(:,:)=cmplx(Smat(:,:),0.0d0,kind=c_double)
  Cmat_c(:,:)=cmplx(Cmat(:,:),0.0d0,kind=c_double)
  SC_c(:,:)=Smat_c(:,:)+Cmat_c(:,:)
  ! Bare (static) pairing vertex from Kanamori model
  if(sw_pair)then
     V0_c(:,:)=cmplx(0.5d0*(Smat(:,:)+Cmat(:,:)),0.0d0,kind=c_double) !bare Vud=(C+S)/2
  else
     V0_c(:,:)=cmplx(0.5d0*(Smat(:,:)-Cmat(:,:)),0.0d0,kind=c_double) !bare Vuu=(S-C)/2
  end if
  !$omp parallel do collapse(2) private(i,cmat1,cmat2,cmat3,cmat4,cmat5,ipiv,info,l)
  wloop:do j=1,Nw
     qloop:do i=1,Nk
        ! cmat1 = -χ_0_S·S, cmat2 = χ_0_C·C
        call zgemm('N','N',Nchi,Nchi,Nchi,(-1.0d0,0.0d0),Vsigma(i,j,:,:),Nchi,Smat_c,Nchi,(0.0d0,0.0d0),cmat1,Nchi)
        call zgemm('N','N',Nchi,Nchi,Nchi,(1.0d0,0.0d0), Vdelta(i,j,:,:),Nchi,Cmat_c,Nchi,(0.0d0,0.0d0),cmat2,Nchi)
        cmat3(:,:)=-cmat1(:,:) ! χ_0_S·S (RHS for first zgesv)
        cmat4(:,:)= cmat2(:,:) ! χ_0_C·C (RHS for second zgesv)
        cmat5(:,:)=cmat3(:,:)+cmat4(:,:) ! χ_0_S·S + χ_0_C·C  (cache before zgesv overwrites cmat3,cmat4)
        do l=1,Nchi
           cmat1(l,l)=cmat1(l,l)+1.0d0 ! I - χ_0_S·S
           cmat2(l,l)=cmat2(l,l)+1.0d0 ! I + χ_0_C·C
        end do
        ! χ_s, χ_c via RPA inversion
        call zgesv(Nchi,Nchi,cmat1,Nchi,ipiv,cmat3,Nchi,info) ! (I-χ_0_S·S) X = χ_0_S·S  -> cmat3 = χ_s
        if(info/=0)then; print*,'zgesv failed: info=',info; stop; end if
        call zgesv(Nchi,Nchi,cmat2,Nchi,ipiv,cmat4,Nchi,info) ! (I+χ_0_C·C) X = χ_0_C·C  -> cmat4 = χ_c
        if(info/=0)then; print*,'zgesv failed: info=',info; stop; end if
        cmat1(:,:)=cmat3(:,:) ! χ_s
        cmat2(:,:)=cmat4(:,:) ! χ_c
        ! ----- V_σ (self-energy vertex) -------------------------------------
        cmat3(:,:)=1.5d0*Smat(:,:)-0.5d0*Cmat(:,:)   ! static HF-like bare vertex
        call zgemm('N','N',Nchi,Nchi,Nchi,(1.5d0,0.0d0),Smat_c,Nchi,cmat1,Nchi,(1.0d0,0.0d0),cmat3,Nchi)  ! +3/2·S·χ_s
        call zgemm('N','N',Nchi,Nchi,Nchi,(0.5d0,0.0d0),Cmat_c,Nchi,cmat2,Nchi,(1.0d0,0.0d0),cmat3,Nchi)  ! +1/2·C·χ_c
        call zgemm('N','N',Nchi,Nchi,Nchi,(-0.25d0,0.0d0),SC_c,Nchi,cmat5,Nchi,(1.0d0,0.0d0),cmat3,Nchi)  ! -1/4·(S+C)·(χ_s+χ_c)  double-count fix
        ! ----- V_Δ (pairing vertex) -----------------------------------------
        cmat4(:,:)=V0_c(:,:)   ! start from bare static vertex
        if(sw_pair)then ! singlet
           call zgemm('N','N',Nchi,Nchi,Nchi,(1.5d0,0.0d0),Smat_c,Nchi,cmat1,Nchi,(1.0d0,0.0d0),cmat4,Nchi)  ! +3/2·S·χ_s
           call zgemm('N','N',Nchi,Nchi,Nchi,(-0.5d0,0.0d0),Cmat_c,Nchi,cmat2,Nchi,(1.0d0,0.0d0),cmat4,Nchi) ! -1/2·C·χ_c
        else            ! triplet
           call zgemm('N','N',Nchi,Nchi,Nchi,(-0.5d0,0.0d0),Smat_c,Nchi,cmat1,Nchi,(1.0d0,0.0d0),cmat4,Nchi) ! -1/2·S·χ_s
           call zgemm('N','N',Nchi,Nchi,Nchi,(-0.5d0,0.0d0),Cmat_c,Nchi,cmat2,Nchi,(1.0d0,0.0d0),cmat4,Nchi) ! -1/2·C·χ_c
        end if
        Vsigma(i,j,:,:)=cmat3(:,:)   ! V_σ
        Vdelta(i,j,:,:)=cmat4(:,:)   ! V_Δ
     end do qloop
  end do wloop
  !$omp end parallel do
end subroutine mkV_flex_nosoc

subroutine get_chi0sc(Vsigma,Vdelta,Gk,fk,kmap,invk,irr_chi,chi_map,olist,prt,temp,Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi,sw_pair)
  !> Build SC irreducible susceptibility in the FLEX spin/charge channels.
  !>   χ^GG  = G(k+q)·G(k)              (normal bubble)
  !>   χ^FF  = -F(k+q)·F^*(k)           (anomalous bubble; sign/structure handled by sw_pair)
  !>   χ_0_S = χ^GG + χ^FF              -> Vsigma  (spin-channel χ_0)
  !>   χ_0_C = χ^GG - χ^FF              -> Vdelta  (charge-channel χ_0)
  !> Output ordering (Nk,Nw,Nchi,Nchi) matches the rest of the FLEX pipeline.
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nx,Ny,Nz,Nk,Nw,Norb,Nkall,Nchi
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
  integer(c_int32_t),dimension(Nchi,Nchi,2)::chi_map
  integer(c_int32_t),dimension(Nchi*(Nchi+1)/2,2)::irr_chi
  logical(1),intent(in):: sw_pair
  real(c_double),intent(in):: temp
  real(c_double),intent(in),dimension(Norb):: prt
  complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk,fk
  complex(c_double),intent(out),dimension(Nk,Nw,Nchi,Nchi):: Vsigma,Vdelta

  integer(c_int32_t) i,j
  complex(c_double),dimension(Nchi,Nchi):: cmat1,cmat2

  call get_chi0_conv(Vsigma,Gk,kmap,invk,irr_chi,chi_map,olist,temp,Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi)            ! Vsigma <- χ^GG
  call get_chi0_conv_ff(Vdelta,fk,kmap,invk,irr_chi,chi_map,olist,prt,temp,Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi,sw_pair) ! Vdelta <- χ^FF
  !$omp parallel do collapse(2) private(i,j,cmat1,cmat2)
  do j=1,Nw
     do i=1,Nk
        cmat1(:,:)=Vsigma(i,j,:,:)+Vdelta(i,j,:,:) ! χ_0_S = χ^GG + χ^FF
        cmat2(:,:)=Vsigma(i,j,:,:)-Vdelta(i,j,:,:) ! χ_0_C = χ^GG - χ^FF
        Vsigma(i,j,:,:)=cmat1(:,:)
        Vdelta(i,j,:,:)=cmat2(:,:)
     end do
  end do
  !$omp end parallel do
end subroutine get_chi0sc

subroutine get_chi0_conv_ff(chi,Fk,kmap,invk,irr_chi,chi_map,olist,prt,temp,Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi,sw_pair)
  !> Anomalous bubble χ^FF(q,iω) = -Σ_k F(k+q, iω+iω') · F^*(k, iω') via real-space FFT.
  !>
  !> Algorithm:
  !>   - Map F(k,iω) to a 2*Nw frequency grid using TRS / particle-hole symmetry, with a
  !>     spin-symmetry sign sgn = 1 (singlet, even-ω F) or -1 (triplet d_z, odd-ω F).
  !>   - tmpfk14 stores F†_{lm}(k,iω) in real space (complex-conjugate component for one leg).
  !>     tmpfk23 stores F_{lm}(k,iω) directly (other leg).
  !>   - FFT both to real space, multiply tmpfk14(r,τ)·tmpfk23(-r,-τ), inverse FFT.
  !>   - Result is the χ^FF bubble at chi-orbital pair (l,m) and its mirror (m,l) via Hermitian conjugate.
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nw,Norb,Nchi,Nkall,Nk,Nx,Ny,Nz
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
  integer(c_int32_t),intent(in),dimension(Nchi,Nchi,2):: chi_map
  integer(c_int32_t),intent(in),dimension(Nchi*(Nchi+1)/2,2):: irr_chi
  logical(1),intent(in):: sw_pair
  real(c_double),intent(in):: temp
  real(c_double),intent(in),dimension(Norb):: prt
  complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: Fk
  complex(c_double),intent(out),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(c_int32_t) i,j,k,l,m,n,iorb
  integer(c_int32_t) ii(0:Nx-1),ij(0:Ny-1),ik(0:Nz-1),iw(2*Nw)
  real(c_double) weight,sgn,dprt
  real(c_double),parameter:: eps=1.0d-9
  complex(c_double),dimension(0:Nx-1,0:Ny-1,0:Nz-1,2*Nw):: tmp,tmpfk14,tmpfk23
  
  ! sgn: F(-k,ω) = sgn·F^T(k,ω);  dprt: Δ(k,-ω) = dprt·Δ*(k,ω)
  ! singlet: F(-k)=+F^T(k), even-ω gap; triplet: F(-k)=-F^T(k), odd-ω gap
  if(sw_pair)then
     sgn=1.0d0 !singlet Fk=F^T-k
     dprt=1.0d0 !gap parity even freq singlet is even
  else
     sgn=-1.0d0 !triplet Fk=-F^T-k
     dprt=-1.0d0 !gap parity even freq triplet is odd
  end if
  weight=temp/dble(Nkall)
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
  call init_fft_plans(tmpfk14,tmp,Nx,Ny,Nz,2*Nw)
  orb_loop:do iorb=1,Nchi*(Nchi+1)/2
     l=irr_chi(iorb,1)
     m=irr_chi(iorb,2)
     !use symmetry G^lm(k,iw)=G^ml(k,-iw) from Hermitian symmetry of Hamiltonian
     !$omp parallel do private(i)
     w_loop_Gk_to_tmp:do j=1,Nw
        k_loop_Gk_to_tmp:do i=1,Nkall
           if(invk(2,i)==0)then
              !iw
              tmpfk14(kmap(1,i),kmap(2,i),kmap(3,i),j)=conjg(Fk(invk(1,i),j,olist(m,2),olist(l,1))) !F^+14(k,iw)
              tmpfk23(kmap(1,i),kmap(2,i),kmap(3,i),j)=Fk(invk(1,i),j,olist(l,2),olist(m,1)) !F23(k,iw)
              !-iw
              tmpfk14(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=prt(olist(l,1))*prt(olist(m,2))*sgn*dprt*conjg(Fk(invk(1,i),j,olist(l,1),olist(m,2))) !F^+14(k,-iw)=sgnF^*14(k,iw)
              tmpfk23(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=prt(olist(m,1))*prt(olist(l,2))*sgn*dprt*Fk(invk(1,i),j,olist(m,1),olist(l,2)) !F23(k,-iw)=sgnF32(k,iw)
           else if(invk(2,i)==1)then
              !iw
              tmpfk14(kmap(1,i),kmap(2,i),kmap(3,i),j)=conjg(Fk(invk(1,i),j,olist(m,1),olist(l,1))) !F^+14(-k,iw)=F41(k,iw)
              tmpfk23(kmap(1,i),kmap(2,i),kmap(3,i),j)=Fk(invk(1,i),j,olist(l,2),olist(m,2)) !F23(-k,iw)=F^23(k,iw)
              !-iw
              tmpfk14(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=prt(olist(l,1))*prt(olist(m,2))*sgn*dprt*conjg(Fk(invk(1,i),j,olist(l,1),olist(m,2))) !F^+14(-k,-iw)
              tmpfk23(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=prt(olist(m,1))*prt(olist(l,2))*sgn*dprt*Fk(invk(1,i),j,olist(l,2),olist(m,1)) !F23(-k,-iw)
           end if
        end do k_loop_Gk_to_tmp
     end do w_loop_Gk_to_tmp
     !$omp end parallel do
     call FFT(tmpfk14,tmp,Nx,Ny,Nz,2*Nw,.true.)
     call FFT(tmpfk23,tmp,Nx,Ny,Nz,2*Nw,.true.)
     !calculate G(r)G(-r)
     !$omp parallel do private(i,j,k)
     w_loop_conv:do n=1,2*Nw
        z_loop:do k=0,Nz-1
           y_loop:do j=0,Ny-1
              !$omp simd
              x_loop:do i=0,Nx-1
                 tmp(i,j,k,n)=-tmpfk14(i,j,k,n)*tmpfk23(ii(i),ij(j),ik(k),iw(n))
              end do x_loop
              !$omp end simd
           end do y_loop
        end do z_loop
     end do w_loop_conv
     !$omp end parallel do
     call FFT(tmp,tmpfk14,Nx,Ny,Nz,2*Nw,.false.)
     !$omp parallel do private(i,j)
     w_loop_tmp_to_chi:do j=1,Nw
        k_loop_tmp_to_chi:do i=1,Nkall
           if(invk(2,i)==0)then
              chi(invk(1,i),j,m,l)=tmp(kmap(1,i),kmap(2,i),kmap(3,i),j)*weight
              if(abs(dble(chi(invk(1,i),j,m,l)))<eps) chi(invk(1,i),j,m,l)=cmplx(0.0d0,imag(chi(invk(1,i),j,m,l)),kind=c_double)
              if(abs(imag(chi(invk(1,i),j,m,l)))<eps) chi(invk(1,i),j,m,l)=cmplx(dble(chi(invk(1,i),j,m,l)),0.0d0,kind=c_double)
           end if
        end do k_loop_tmp_to_chi
     end do w_loop_tmp_to_chi
     !$omp end parallel do
     chi(:,:,chi_map(m,l,1),chi_map(m,l,2))=conjg(chi(:,:,m,l))
  end do orb_loop
  call destroy_fft_plans()
end subroutine get_chi0_conv_ff