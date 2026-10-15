! Linearized Eliashberg equation WITH spin-orbit coupling, using TRS and (optionally) inversion.
!
! Basis and storage
!   Norb counts SPIN-orbitals: slist(a) = +-1 is the spin of a, invs(a) its opposite-spin partner.
!   Gk and the pairing vertex chi live on the Nk irreducible TRS k-points; the gap lives on the
!   Ndk points chosen by get_delta_kmap_soc, and idk/ipk say where each full-BZ point i went:
!     inversion is a symmetry (gap_sym /= gap_unit, par_mix)
!         parity is good, so Delta(-k) = pp P Delta(k) P is REBUILT instead of stored:
!         Ndk = Nk, idk = invk(1,:), ipk = invk(2,:) (=1 on the rebuilt -k half)
!     parity mixing (gap_sym = par_mix)
!         nothing relates k and -k: Ndk = Nkall, idk = identity, ipk = 0
!   lin_eliash_soc expands the result back onto the full BZ in delta_full, so callers see no
!   difference.  Only w_n > 0 is stored; w_n < 0 comes from Pauli, F(k,-w) = -F^T(-k,w).
!   pp = +-1 is the parity EIGENVALUE of the gap (from gap_sym), prt(a) = +-1 the parity of
!   orbital a: Delta_ab(-k) = pp prt(a) prt(b) Delta_ab(k).  (Unlike lin_eliash there is no iprm,
!   so inversion must not permute the orbitals; use par_mix for cells where it does.)
!
! What gap_sym fixes
!   Its SIGN is the inversion parity of the gap (>0 even, <0 odd) and its magnitude picks the
!   band-basis seed (get_initial_delta_soc).  It is NOT a singlet/triplet switch: with SOC a parity
!   sector holds BOTH the (pseudospin singlet, even w) and the (pseudospin triplet, odd w) pairing,
!   and no exact symmetry separates them.  Two mechanisms keep the result meaningful:
!     project_spin   only when G is SU(2) invariant (no SOC, no Zeeman; detected at run time by
!                    is_spin_rotation_invariant).  There singlet and triplet ARE separate sectors,
!                    so the seed's one is kept and round-off into the other removed, and the
!                    singlet/triplet of lin_eliash are reproduced exactly.
!     select_branch  with SOC, Arnoldi returns the largest eigenvalue whose Ritz vector carries the
!                    spin character of gap_sym rather than the largest of the parity sector, so
!                    lambda stays continuous as the SOC strength goes to zero.
!   The power method (arnoldi_m = 0) has no branch tracking and converges to the dominant state of
!   the sector, so use Arnoldi when SOC is present.
!
! Real-linear kernel
!   mkfk_trs_soc and mkdelta_soc compute part of F/Delta and fill the rest by conjugation
!   (TRS: Delta_uu = -Delta_dd^+, Delta_ud Hermitian), so apply_op is linear over the REALS only.
!   Arnoldi therefore uses a real inner product (get_inner) and real Ritz coefficients (real_ritz);
!   a complex Krylov space mixes the two and silently returns wrong eigenvalues.

subroutine lin_eliash_soc(delta_full,chi,Gk,uni,init_delta,Vmat,sgnsig,sgnsig2,prt,olist,slist,kmap,invk,invs,invschi,&
     temp,eps,Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,itemax,gap_sym,arnoldi_m,lambda_out) bind(C)
  !> calculate linearized eliashberg equations with soc and TRS
  !!@param   lambda_out,out: physical (largest) Eliashberg eigenvalue
  !!@param delta_full,out: gap function on the full BZ.  Internally Δ is stored on Ndk k-points
  !!                       only (see get_delta_kmap_soc): with inversion the -k half is rebuilt
  !!                       from parity, so Ndk=Nk; for parity mixing (gap_sym=-10) Ndk=Nkall
  !!@param         Gk,in: normal green function
  !!@param        uni,in: unitary matrix
  !!@param init_delta,in: band basis initial gap function
  !!@param       Vmat,in: V-matrix
  !!@param      olist,in: property of chi index
  !!@param      slist,in: spin property of orbital index
  !!@param       kmap,in: property of k-points
  !!@param       invk,in: list of reverse k-points
  !!@param       invs,in: list of spin flipped orbital index
  !!@param       temp,in: Temperature
  !!@param        eps,in: Threshold of calculation
  !!@param      Nkall,in: Number of all k-points
  !!@param         Nk,in: Number of k-points
  !!@param         Nw,in: Number of Matsubara frequencies
  !!@param       Nchi,in: Number of footnote of chi
  !!@param       Norb,in: Number of orbitals
  !!@param         Nx,in: Number of kx mesh
  !!@param         Ny,in: Number of ky mesh
  !!@param         Nz,in: Number of kz mesh
  !!@param     itemax,in: maximum iteration of power method
  !!@param    gap_sym,in: gap symmetry number.  Its SIGN is the inversion parity of the gap
  !!                      (>0 even, <0 odd) and its magnitude picks the seed; 0 = unit-matrix seed,
  !!                      -10 = parity mixing.  See the notes at the top of this file.
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  use eliash_krylov
  implicit none
  integer(c_int64_t),intent(in):: Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,itemax,gap_sym
  integer(c_int64_t),intent(in):: arnoldi_m  !! Krylov subspace dimension (0=power method, >0=Arnoldi)
  real(c_double),intent(out):: lambda_out
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),dimension(Norb):: slist,invs
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
  integer(c_int64_t),intent(in),dimension(Nchi):: invschi
  real(c_double),intent(in):: temp,eps
  real(c_double),intent(in),dimension(Nchi,Nchi):: Vmat
  real(c_double),intent(in),dimension(Norb,Norb):: sgnsig
  real(c_double),intent(in),dimension(Nchi,Nchi):: sgnsig2
  real(c_double),intent(in),dimension(Norb):: prt
  real(c_double),intent(in),dimension(Nk,Norb):: init_delta
  complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(c_double),intent(out),dimension(Nkall,Nw,Norb,Norb):: delta_full
  complex(c_double),intent(inout),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(c_int32_t) i_iter,i_eig,count,i,l,m
  integer(c_int32_t),parameter:: eig_max=2
  !! special values of gap_sym (its sign is the parity of the gap, see the notes at the top)
  logical(1),parameter:: no_dcsub=.false.     !! pairing channel: see get_V_soc_flex
  integer(c_int64_t),parameter:: gap_unit=0   !! unit-matrix seed: no parity and no spin projection
  integer(c_int64_t),parameter:: par_mix=-10  !! parity mixing: no parity, gap kept on the full BZ
  integer(c_int64_t) Ndk
  integer(c_int64_t),dimension(Nkall):: idk,ipk
  real(c_double) norm,inorm,norm2,weight,lambda_rq,vec_err,pp
  logical su2
  complex(c_double),allocatable,dimension(:,:,:,:):: delta,newdelta,fk
  !> vector geometry for the shared Krylov routines (idk says which stored point carries each
  !> full-BZ point, so norms and inner products still run over the whole BZ)
  type(kspace) ks

  lambda_out=0.0d0
  ! Δ is stored on Ndk points: idk(i) = stored index of full k-point i, ipk(i)=1 if Δ(k_i) is
  ! rebuilt by parity, Δ(-k) = pp·P Δ(k) P  (pp=+1 singlet, -1 triplet)
  call get_delta_kmap_soc(idk,ipk,Ndk,invk,Nkall,Nk,gap_sym)
  pp=merge(1.0d0,-1.0d0,gap_sym>0)
  allocate(delta(Ndk,Nw,Norb,Norb),newdelta(Ndk,Nw,Norb,Norb),fk(Ndk,Nw,Norb,Norb))
  call kspace_init(ks,idk,Nkall,Ndk,Nw,Norb)
  ! Without SOC (G spin-diagonal and spin-independent) singlet and triplet are separate sectors:
  ! the kernel keeps the seed's spin structure, but round-off would grow the other sector, so it
  ! is projected out (project_spin).  With SOC they mix and only the parity is imposed.
  su2=is_spin_rotation_invariant()
  if(gap_sym/=gap_unit .and. gap_sym/=par_mix)then
     if(su2)then
        print'(A)','no SOC: gap restricted to the seed spin sector (singlet / S_z=0 triplet)'
     else
        print'(A)','SOC: singlet and triplet mix, only the parity of the gap is fixed'
     end if
  end if
  weight=temp/dble(Nkall)
  ! norm2 is the shift used in the power iteration: newdelta = -weight*K*delta + norm2*delta
  ! converges to (λ + norm2)*delta, so true eigenvalue λ = Rayleigh quotient - norm2.
  ! After the 1st eigenvalue is found, norm2 is set to its norm for deflation (shift-deflation).
  norm2=0.0d0
  call get_V_soc_flex(chi,Vmat,olist,slist,Nk,Nw,Nchi,no_dcsub)   !pairing: no double-counting subtraction
  ! pairing vertex V_Delta = -V/2 + V chi V: get_V_soc_flex returns -V + V chi V (the form the
  ! self-energy uses), but in the gap equation the antisymmetrized bare vertex is both the direct
  ! and the exchange term and both orderings F_{l2 m1} of the pair are summed, so -V would count
  ! the bare term twice (the RPA ladder V chi V is direct-channel only).  With 1/2 the singlet
  ! channel of a spin-conserving system reduces to (S+C)/2 + 3/2 S chi_s S - 1/2 C chi_c C,
  ! i.e. lin_eliash.  The w_n->inf slot in mkdelta_soc is -Vmat/2 accordingly.
  !$omp parallel do private(l,m)
  do l=1,int(Nchi,c_int32_t)
     do m=1,int(Nchi,c_int32_t)
        chi(:,:,m,l)=chi(:,:,m,l)+0.5d0*Vmat(m,l)
     end do
  end do
  !$omp end parallel do
  print'(A15,2E16.8)','V_delta max is ',maxval(dble(chi)),maxval(aimag(chi))
  print'(A15,2E16.8)','V_delta min is ',minval(dble(chi)),minval(aimag(chi))
  if(arnoldi_m>0)then !--- Arnoldi solver (shared with lin_eliash, see module eliash_krylov)
     call get_initial_delta_soc(newdelta,init_delta,uni,kmap,slist,invk,invs,idk,ipk,&
          Nkall,Nk,Ndk,Nw,Norb,gap_sym)
     if(gap_sym==gap_unit .or. gap_sym==par_mix)then
        call arnoldi_two_pass(lambda_out,delta,newdelta,ks,arnoldi_m,apply_op)
     else
        ! with SOC the parity sector holds two branches, so the solver is told which one to take
        call arnoldi_two_pass(lambda_out,delta,newdelta,ks,arnoldi_m,apply_op,classify_branch)
     end if
  else !--- Power method
     eigenval_loop:do i_eig=1,eig_max !solve eig_val using power method, 1st eig is usually large negative value
     call get_initial_delta_soc(delta,init_delta,uni,kmap,slist,invk,invs,idk,ipk,Nkall,Nk,Ndk,Nw,Norb,gap_sym)
     count=0 !count too small eigenvalue
     iter_loop:do i_iter=1, itemax !iteration
        call mkfk_trs_soc(fk,Gk,delta,sgnsig,slist,invk,invs,idk,ipk,Nkall,Nk,Ndk,Nw,Norb,gap_sym)
        call mkdelta_soc(newdelta,fk,chi,Vmat,sgnsig,sgnsig2,prt,kmap,invk,invs,invschi,olist,slist,idk,ipk,&
             Nkall,Nk,Ndk,Nw,Nchi,Norb,Nx,Ny,Nz,gap_sym)
        call project_parity(newdelta)
        call project_spin(newdelta)
        !$omp parallel workshare
        newdelta(:,:,:,:)=-newdelta(:,:,:,:)*weight+delta(:,:,:,:)*norm2
        !$omp end parallel workshare
        norm=kv_norm(ks,newdelta)
        if(norm <= 0.0d0)then
           print*,'Error: norm <= 0 in iter_loop'
           stop
        end if
        inorm=1.0d0/norm
        ! Rayleigh quotient: lambda_rq = Re(<delta, newdelta>) assuming ||delta||=1
        ! shift is +norm2, so true eigenvalue = lambda_rq - norm2
        lambda_rq=kv_dot(ks,delta,newdelta)
        ! Vector convergence: ||newdelta/norm - delta||
        vec_err=kv_err(ks,newdelta,delta,inorm)
        if(abs(lambda_rq-norm2)>=1.0d2 .or. abs(lambda_rq-norm2)<1.0d-6)then
           print'(I3,A13,2E16.8)',i_iter,' lambda_elsh=',lambda_rq-norm2
        else
           print'(I3,A13,2F12.8)',i_iter,' lambda_elsh=',lambda_rq-norm2
        end if
        if(vec_err<eps)then
           if((lambda_rq-norm2)>1.0d-1)then !do not finish until eig>0.1
              exit
           else if((lambda_rq-norm2)<0.0d0)then !negative eigenvalue: exit immediately
              exit
           else !consider small positive eig
              if(vec_err<eps*1.0d-2)then
                 count=count+1
              end if
              if(count>30)exit !if eigenvalue <0.1  until 30 count exit
           end if
        end if
        !$omp parallel workshare
        delta(:,:,:,:)=newdelta(:,:,:,:)*inorm
        !$omp end parallel workshare
     end do iter_loop
     print*,'eliash=',lambda_rq-norm2
     lambda_out=lambda_rq-norm2
     if(i_eig==1 .and. (lambda_rq-norm2)>0.0d0)then
        print*,'1st eigenvalue is positive: skipping 2nd loop'
        exit
     end if
     norm2=norm
  end do eigenval_loop
  norm=kv_norm(ks,newdelta)
  inorm=1.0d0/norm
  !$omp parallel workshare
  delta(:,:,:,:)=newdelta(:,:,:,:)*inorm
  !$omp end parallel workshare
  end if
  ! expand the stored Δ to the full BZ (-k half from parity)
  !$omp parallel do private(l,m)
  do i=1,int(Nkall,c_int32_t)
     do m=1,int(Norb,c_int32_t)
        do l=1,int(Norb,c_int32_t)
           delta_full(i,:,l,m)=merge(pp*prt(l)*prt(m),1.0d0,ipk(i)==1)*delta(idk(i),:,l,m)
        end do
     end do
  end do
  !$omp end parallel do
  deallocate(delta,newdelta,fk)
contains
  subroutine project_parity(d)
    !> Impose the parity of the seed where it is not built in: Δ(-k) of the other k-points is
    !> rebuilt from parity (Δ(-k) = pp·P Δ(k) P), but a k=-k point (TRIM) is stored once, so there
    !> Δ_ab(K) <- [1 + pp·prt_a·prt_b]/2 · Δ_ab(K).  Without it an odd-parity triplet keeps a
    !> spurious TRIM component.  Skipped for parity mixing, which is stored on the full BZ.
    complex(c_double),intent(inout),dimension(Ndk,Nw,Norb,Norb):: d
    integer(c_int32_t) ii,a,b
    if(gap_sym==gap_unit .or. gap_sym==par_mix)return
    !$omp parallel do private(a,b)
    do ii=1,int(Nkall,c_int32_t)
       if(ipk(ii)==0 .and. invk(3,ii)==ii)then
          do b=1,int(Norb,c_int32_t)
             do a=1,int(Norb,c_int32_t)
                d(idk(ii),:,a,b)=0.5d0*(1.0d0+pp*prt(a)*prt(b))*d(idk(ii),:,a,b)
             end do
          end do
       end if
    end do
    !$omp end parallel do
  end subroutine project_parity

  logical function is_spin_rotation_invariant()
    !> true if G has no spin-flip part and equal up/down blocks (no SOC, no Zeeman field)
    integer(c_int32_t) a,b
    real(c_double) tol,dev
    tol=1.0d-12*maxval(abs(Gk))
    dev=0.0d0
    do b=1,int(Norb,c_int32_t)
       do a=1,int(Norb,c_int32_t)
          if(slist(a)/=slist(b))then
             dev=max(dev,maxval(abs(Gk(:,:,a,b))))
          else if(slist(a)==1)then
             dev=max(dev,maxval(abs(Gk(:,:,a,b)-Gk(:,:,invs(a),invs(b)))))
          end if
       end do
    end do
    is_spin_rotation_invariant=(dev<=tol)
  end function is_spin_rotation_invariant

  subroutine project_spin(d)
    !> Spin-rotation-invariant case only: keep the spin structure of the seed, Δ = Δ_orb ⊗ χ with
    !> χ = iσ_y (singlet, gap_sym>0) or σ_z iσ_y (S_z=0 triplet d_z, gap_sym<0):
    !> Δ_{a↓,b↑} = ∓Δ_{a↑,b↓} (- singlet, + triplet), same-spin parts 0.  Δ↓↑ is copied from Δ↑↓
    !> rather than averaged with it: Δ↑↓ carries the TRS gauge (its lower triangle is filled by
    !> Hermiticity in mkdelta_soc), and mixing in the independently computed Δ↓↑ breaks that.
    complex(c_double),intent(inout),dimension(Ndk,Nw,Norb,Norb):: d
    integer(c_int32_t) a,b
    real(c_double) sgs
    if(.not.su2 .or. gap_sym==gap_unit .or. gap_sym==par_mix)return
    sgs=merge(-1.0d0,1.0d0,gap_sym>0)
    !$omp parallel do private(a)
    do b=1,int(Norb,c_int32_t)
       do a=1,int(Norb,c_int32_t)
          if(slist(a)==slist(b))then
             d(:,:,a,b)=0.0d0
          else if(slist(a)==1)then
             d(:,:,invs(a),invs(b))=sgs*d(:,:,a,b)
          end if
       end do
    end do
    !$omp end parallel do
  end subroutine project_spin

  subroutine apply_op(w_out,v_in)
    !> one application of the Eliashberg kernel, w_out = (T/Nkall)*K*v_in; uses host fk as scratch.
    !> Flat dummies: the shared driver stores Krylov vectors as plain vectors, and the rank-4
    !> routines below pick up the (Ndk,Nw,Norb,Norb) shape by sequence association.
    complex(c_double),intent(in):: v_in(*)
    complex(c_double),intent(out):: w_out(*)
    integer(c_int32_t) q
    call mkfk_trs_soc(fk,Gk,v_in,sgnsig,slist,invk,invs,idk,ipk,Nkall,Nk,Ndk,Nw,Norb,gap_sym)
    call mkdelta_soc(w_out,fk,chi,Vmat,sgnsig,sgnsig2,prt,kmap,invk,invs,invschi,olist,slist,idk,ipk,&
         Nkall,Nk,Ndk,Nw,Nchi,Norb,Nx,Ny,Nz,gap_sym)
    call project_parity(w_out)
    call project_spin(w_out)
    !$omp parallel do private(q)
    do q=1,ks%nvec
       w_out(q)=-w_out(q)*weight
    end do
    !$omp end parallel do
  end subroutine apply_op

  subroutine classify_branch(ok,x)
    !> Branch test handed to the shared Arnoldi driver: does x have the spin character gap_sym
    !> asks for?  With SOC a parity sector holds both the (pseudospin singlet, even w) and the
    !> (pseudospin triplet, odd w) pairing and no exact symmetry separates them, so without this
    !> the solver would return whichever is larger and lambda would jump as SOC goes to zero.
    complex(c_double),intent(in):: x(*)
    logical,intent(out):: ok
    real(c_double) ws,wt
    call get_spin_weight(ws,wt,x)
    ok=merge(ws>wt,wt>ws,gap_sym>0)
  end subroutine classify_branch

  subroutine get_spin_weight(ws,wt,d)
    !> Split |Delta|^2 into its spin-singlet weight ws and spin-triplet weight wt, by projecting
    !> the opposite-spin block on Delta_du = -Delta_ud (singlet) and Delta_du = +Delta_ud (S_z=0
    !> triplet); the equal-spin blocks are the S_z=+-1 triplet and go to wt.  Only ws vs wt is
    !> used (which one dominates), so the overall normalisation does not matter.
    complex(c_double),intent(in),dimension(Ndk,Nw,Norb,Norb):: d
    real(c_double),intent(out):: ws,wt
    integer(c_int32_t) a,b
    ws=0.0d0
    wt=0.0d0
    do b=1,int(Norb,c_int32_t)
       do a=1,int(Norb,c_int32_t)
          if(slist(a)/=slist(b))then
             if(slist(a)==1)then
                ws=ws+sum(abs(d(:,:,a,b)-d(:,:,invs(a),invs(b)))**2)
                wt=wt+sum(abs(d(:,:,a,b)+d(:,:,invs(a),invs(b)))**2)
             end if
          else
             wt=wt+2.0d0*sum(abs(d(:,:,a,b))**2)   ! equal spin: S_z=+-1 triplet
          end if
       end do
    end do
  end subroutine get_spin_weight

end subroutine lin_eliash_soc

subroutine mkfk_trs_soc(fk,Gk,delta,sgnsig,slist,invk,invs,idk,ipk,Nkall,Nk,Ndk,Nw,Norb,gap_sym)
  !
  !>calculate linearized anomalous green function Fk with TRS
  !>if we considder TRS,F_ab(k)=G_ac(k)Delta_cd(k)Gbd(-k)
  !> =G_ac(k)Delta_cd(k)G^*_dbss'(k)ss'
  !> delta and fk live on the Ndk stored k-points (get_delta_kmap_soc); F is computed only there
  !
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nkall,Nk,Ndk,Nw,Norb,gap_sym
  integer(c_int64_t),intent(in),dimension(3,Nkall):: invk
  integer(c_int64_t),intent(in),dimension(Nkall):: idk,ipk
  integer(c_int64_t),intent(in),dimension(Norb):: invs,slist
  real(c_double),intent(in),dimension(Norb,Norb):: sgnsig
  complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(c_double),intent(in),dimension(Ndk,Nw,Norb,Norb):: delta
  complex(c_double),intent(out),dimension(Ndk,Nw,Norb,Norb):: fk
   
  integer(c_int32_t) i,j,l,m,n
  complex(c_double),dimension(Norb,Norb):: gmat1,gmat2,cmat1

  !$omp parallel do collapse(2) private(i,j,l,m,n,gmat1,gmat2,cmat1)
  do j=1,Nw
     do i=1,Nkall
        if(ipk(i)==1)cycle !Δ(k_i) is rebuilt from parity, not stored
        ! F(k) = G(k)·Δ(k)·G(-k)^T, using TRS: G(-k)_{ln} = s_n s_l G(k)_{ī n̄}
        ! invk(2,i)==0: k is direct; invk(2,i)==1: k is TRS-mapped from irr k
        if(invk(2,i)==0)then
           gmat1(:,:)=Gk(invk(1,i),j,:,:)   !G(k)
           do l=1,Norb
              do n=1,Norb
                 gmat2(l,n)=sgnsig(n,l)*Gk(invk(1,i),j,invs(l),invs(n))  !G(-k)_{ln}=s_n*s_l*G(k)_{ī n̄}
              end do
           end do
        else
           do l=1,Norb
              do n=1,Norb
                 gmat1(l,n)=sgnsig(n,l)*Gk(invk(1,i),j,invs(n),invs(l))  !G(k)_{ln} from G(-irr_k)
                 gmat2(l,n)=Gk(invk(1,i),j,n,l)                           !G(-k)=G(irr_k)^T
              end do
           end do
        end if

        ! F = -G(k) · Δ(k) · G(-k)^†  (linearized anomalous Green function)
        call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),gmat1,Norb,delta(idk(i),j,:,:),Norb,(0.0d0,0.0d0),cmat1,Norb)
        call zgemm('N','C',Norb,Norb,Norb,(-1.0d0,0.0d0),cmat1,Norb,gmat2,Norb,(0.0d0,0.0d0),fk(idk(i),j,:,:),Norb)
     end do
  end do
  !$omp end parallel do
  
  ! TRS (Kramers) relations only; they hold with SOC.  No singlet/triplet spin structure is
  ! imposed (with SOC the two mix within a parity sector); gap_sym only fixes the parity.
  do l=1,Norb
     if(slist(l)==1)then
        do m=1,Norb
           !$omp do private(i,j)
           do j=1,Nw
              do i=1,Ndk
                 if(slist(m)==1)then
                    fk(i,j,invs(l),invs(m))=-conjg(fk(i,j,m,l))
                 else
                    if(invs(m)==l)then
                       fk(i,j,invs(l),invs(m))=dble(fk(i,j,m,l))
                    end if
                 end if
              end do
           end do
           !$omp end do
        end do
     end if
  end do
end subroutine mkfk_trs_soc
  
subroutine mkdelta_soc(newdelta,delta,Vdelta,Vmat,sgnsig,sgnsig2,prt,kmap,invk,invs,invschi,olist,slist,idk,ipk,&
     Nkall,Nk,Ndk,Nw,Nchi,Norb,Nx,Ny,Nz,gap_sym)
  !> new gap function Delta = sum V F.  delta (here F) and newdelta live on the Ndk stored k-points
  !> (get_delta_kmap_soc); F at the other points is rebuilt by parity, F(-k) = pp·P F(k) P.
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nkall,Nk,Ndk,Nw,Nchi,Norb,Nx,Ny,Nz,gap_sym
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
  integer(c_int64_t),intent(in),dimension(Nkall):: idk,ipk
  real(c_double),intent(in),dimension(Norb):: prt
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),dimension(Norb):: slist,invs
  integer(c_int64_t),intent(in),dimension(Nchi):: invschi
  real(c_double),intent(in),dimension(Norb,Norb):: sgnsig
  real(c_double),intent(in),dimension(Nchi,Nchi):: Vmat,sgnsig2
  complex(c_double),intent(in),dimension(Nk,Nw,Nchi,Nchi):: Vdelta
  complex(c_double),intent(in),dimension(Ndk,Nw,Norb,Norb):: delta
  complex(c_double),intent(out),dimension(Ndk,Nw,Norb,Norb):: newdelta

  integer(c_int32_t) i,j,k,n,l,m
  real(c_double) pp,sp,sm
  complex(c_double),dimension(0:Nx-1,0:Ny-1,0:Nz-1,2*Nw):: tmpVdelta,tmpfk,tmp
  
  pp=merge(1.0d0,-1.0d0,gap_sym>0) !parity used for k-points with ipk=1
  !$omp parallel workshare
  newdelta(:,:,:,:)=0.0d0
  !$omp end parallel workshare
  do l=1,Nchi
     do m=1,Nchi
        ! all spin components are kept (with SOC singlet and triplet mix within a parity sector):
        ! Duu, Ddu computed, Dud upper triangle; Ddd and the rest of Dud from TRS below
        if(slist(olist(l,1))*slist(olist(m,2))<0 .or. slist(olist(l,1))==1)then !calc Duu, Dud, Ddu
              if(slist(olist(m,2))==1 .or. olist(l,1)>=invs(olist(m,2)))then !Duu all, Dud upper triangle
                 !$omp parallel
                 !$omp do private(i)
                 do i=1,Nkall
                    if(invk(2,i)==0)then !k,0
                       tmpVdelta(kmap(1,i),kmap(2,i),kmap(3,i),1)=Vdelta(invk(1,i),1,m,l) !j=1 corresponds to w_n=0
                    else if(invk(2,i)==1)then !-k,0
                       ! V(-q)_{ml} = s_{m1}s_{m2}s_{l1}s_{l2} * V(q)_{invschi(l),invschi(m)}  (TRS+SOC)
                       tmpVdelta(kmap(1,i),kmap(2,i),kmap(3,i),1)=sgnsig2(m,l)*Vdelta(invk(1,i),1,invschi(l),invschi(m)) !V^ml(-q,0)=V^lm(q,0)
                    end if
                    tmpVdelta(kmap(1,i),kmap(2,i),kmap(3,i),Nw+1)=-0.5d0*Vmat(m,l) !Nw+1 consider w_n=>inf limit (bare -V/2)
                 end do
                 !$omp end do
                 !$omp do private(i,j)
                 do j=2,Nw
                    do i=1,Nkall
                       if(invk(2,i)==0)then !k,iw,k,-iw
                          tmpVdelta(kmap(1,i),kmap(2,i),kmap(3,i),j)=Vdelta(invk(1,i),j,m,l)
                          tmpVdelta(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+2)=conjg(Vdelta(invk(1,i),j,l,m)) !V^ml(q,-iw)=V^*lm(q,iw)
                       else if(invk(2,i)==1)then !-k,iw, -k,-iw
                          tmpVdelta(kmap(1,i),kmap(2,i),kmap(3,i),j)=&
                               sgnsig2(m,l)*Vdelta(invk(1,i),j,invschi(l),invschi(m)) !V^ml(-q,iw)=V^lm(q,iw)
                          tmpVdelta(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+2)=&
                               sgnsig2(m,l)*conjg(Vdelta(invk(1,i),j,invschi(m),invschi(l))) !V^ml(-q,-iw)=V^*ml(q,iw)
                       end if
                    end do
                 end do
                 !$omp end do
                 !$omp do private(i,j,sp,sm)
                 do j=1,Nw
                    do i=1,Nkall
                       ! ipk=1 marks a k-point whose F is not stored but rebuilt from its partner
                       ! by parity, F_ab(k) = pp*prt(a)*prt(b)*F_ab(-k); sp does that for F(k) and
                       ! sm for the F(-k) that supplies the negative frequency.  ipk is 0 for every
                       ! point in the parity-mixing mode, where both factors are 1.
                       sp=merge(pp*prt(olist(l,2))*prt(olist(m,1)),1.0d0,ipk(i)==1)
                       sm=merge(pp*prt(olist(l,2))*prt(olist(m,1)),1.0d0,ipk(invk(3,i))==1)
                       tmpfk(kmap(1,i),kmap(2,i),kmap(3,i),j)=sp*delta(idk(i),j,olist(l,2),olist(m,1))
                       tmpfk(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=-sm*delta(idk(invk(3,i)),j,olist(m,1),olist(l,2)) !F(k,-ω)=-F^T(-k,ω): anomalous GF Matsubara symmetry
                    end do
                 end do
                 !$omp end do
                 !$omp end parallel
                 call FFT(tmpVdelta,tmp,Nx,Ny,Nz,2*Nw,.true.)
                 call FFT(tmpfk,tmp,Nx,Ny,Nz,2*Nw,.true.)
                 !$omp parallel
                 !$omp workshare
                 tmp(:,:,:,:)=0.0d0
                 !$omp end workshare
                 !$omp do private(i,j,k)
                 do n=1,2*Nw
                    do k=0,Nz-1
                       do j=0,Ny-1
                          do i=0,Nx-1
                             tmp(i,j,k,n)=tmpVdelta(i,j,k,n)*tmpfk(i,j,k,n)
                          end do
                       end do
                    end do
                 end do
                 !$omp end do
                 !$omp workshare
                 tmpfk(:,:,:,:)=0.0d0
                 !$omp end workshare
                 !$omp end parallel
                 call FFT(tmp,tmpfk,Nx,Ny,Nz,2*Nw,.false.)
                 !$omp parallel do private(i)
                 do j=1,Nw
                    do i=1,Nkall
                       if(ipk(i)==1)cycle !only stored points
                       newdelta(idk(i),j,olist(l,1),olist(m,2))=newdelta(idk(i),j,olist(l,1),olist(m,2))&
                            +tmp(kmap(1,i),kmap(2,i),kmap(3,i),j)
                    end do
                 end do
                 !$omp end parallel do
              end if
        end if
     end do
  end do
  !$omp parallel
  do l=1,Norb
     if(slist(l)==1)then
        do m=1,Norb
           !$omp do private(i,j)
           do j=1,Nw
              do i=1,Ndk
                 if(slist(m)==1)then !delta_uu(k)=-delta^+_dd(k)
                    newdelta(i,j,invs(l),invs(m))=-conjg(newdelta(i,j,m,l))
                 else
                    if(invs(m)==l)then !Dud(k)=D^+ud(k)
                       newdelta(i,j,l,m)=dble(newdelta(i,j,l,m))
                    else if(l>invs(m))then !Dud(k)=D^+ud(k)
                       newdelta(i,j,invs(m),invs(l))=conjg(newdelta(i,j,l,m))
                    end if
                 end if
              end do
           end do
           !$omp end do
        end do
     end if
  end do
  !$omp end parallel
end subroutine mkdelta_soc

subroutine get_initial_delta_soc(delta,init_delta,uni,kmap,slist,invk,invs,idk,ipk,Nkall,Nk,Ndk,Nw,Norb,gap_sym)
  !> make orbital basis initial gap function
  !!@param     delta,out: gap function
  !!@param init_delta,in: band basis initial gap function
  !!@param        uni,in: unitary matrix
  !!@param       kmap,in: property of k-points
  !!@param       invk,in: footnote of reverse k
  !!@param      Nkall,in: Number of all k-points
  !!@param         Nk,in: Number of irreducible k-points
  !!@param       Norb,in: Number of orbitals
  !!@param    gap_sym,in: gap symmetry number
  !!@param    idk,ipk,in: stored index / parity-rebuilt flag of each k-point (get_delta_kmap_soc);
  !!                      delta is set on the Ndk stored points only
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  use constant
  implicit none
  integer(c_int64_t),intent(in):: Nw,Norb,Nkall,Nk,Ndk,gap_sym
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
  integer(c_int64_t),intent(in),dimension(Nkall):: idk,ipk
  integer(c_int64_t),intent(in),dimension(Norb):: slist,invs
  real(c_double),intent(in),dimension(Nk,Norb):: init_delta
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(c_double),intent(out),dimension(Ndk,Nw,Norb,Norb):: delta

  integer(c_int32_t) i,j,k,l,m,n
  integer(c_int64_t),parameter:: par_mix=-10
  real(c_double) norm
  complex(c_double),dimension(Norb,Norb):: matL,matR,tmpm

  !$omp parallel do private(i,k,l,m,n,matL,matR,tmpm)
  do i=1,Nkall
     if(ipk(i)==1)cycle !rebuilt from parity
     k=invk(1,i)
     if(invk(2,i)==0)then
        matL=uni(:,:,k)
        do n=1,Norb
           matL(:,n)=matL(:,n)*init_delta(k,n)
        end do
        do l=1,Norb
           matR(:,l)=slist(l)*conjg(uni(invs(l),:,k))
        end do
        call zgemm('N','N',Norb,Norb,Norb,(-1.0d0,0.0d0),matL,Norb,matR,Norb,(0.0d0,0.0d0),delta(idk(i),1,:,:),Norb)
     else if(invk(2,i)==1)then
        do m=1,Norb
           matL(m,:)=slist(m)*conjg(uni(invs(m),:,k))
        end do
        do n=1,Norb
           matL(:,n)=matL(:,n)*init_delta(k,n)
        end do
        matR=transpose(uni(:,:,k))
        call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),matL,Norb,matR,Norb,(0.0d0,0.0d0),delta(idk(i),1,:,:),Norb)
     else
        delta(idk(i),1,:,:)=0.0d0
     end if

     ! The Kramers-partner pairing above is a (pseudo)spin SINGLET, U f (i sigma_y) U^T, and its -k half
     ! is filled with the even-parity singlet rule, whatever the parity of init_delta.  For a triplet
     ! seed turn it into d_z: sigma_z (i sigma_y) = real-spin sign on the row index, and give the -k half
     ! odd parity.  Otherwise a p-wave seed is an even-parity singlet, and the triplet symmetrization in
     ! mkdelta_soc turns it into an even-parity (odd-frequency) triplet.
     if(gap_sym<0 .and. gap_sym/=par_mix)then
        do m=1,Norb
           delta(idk(i),1,m,:)=slist(m)*delta(idk(i),1,m,:)
        end do
        if(invk(2,i)==1)delta(idk(i),1,:,:)=-delta(idk(i),1,:,:)
     end if

     if(gap_sym>0)then
        do l=1,Norb
           do m=1,Norb
              if(slist(m)*slist(l)>=0) delta(idk(i),1,m,l)=0.0d0
           end do
        end do
     end if
  end do
  !$omp end parallel do

  !$omp parallel do private(j)
  do j=2,Nw
     delta(:,j,:,:)=delta(:,1,:,:)
  end do
  !$omp end parallel do

  norm=0.0d0
  !$omp parallel
  do l=1,Norb
     do m=1,Norb
        !$omp do private(j,i) reduction(+:norm)
        do j=1,Nw
           do i=1,Nkall
              norm=norm+abs(delta(idk(i),j,m,l))*abs(delta(idk(i),j,m,l)) !full-BZ norm
           end do
        end do
        !$omp end do
     end do
  end do
  !$omp workshare
  delta(:,:,:,:)=delta(:,:,:,:)/sqrt(2.0d0*norm)
  !$omp end workshare
  !$omp end parallel
end subroutine get_initial_delta_soc

subroutine get_delta_kmap_soc(idk,ipk,Ndk,invk,Nkall,Nk,gap_sym)
  !> k-point storage of the SOC gap function.  With inversion (any gap_sym but 0 and the parity
  !> mixing -10) parity is a good quantum number, Δ(-k) = pp·P Δ(k) P, so Δ is stored on the
  !> irreducible TRS list (Ndk=Nk) and the -k half rebuilt; for parity mixing all Nkall points.
  !!@param   idk,out: stored index of each full k-point
  !!@param   ipk,out: 1 if Δ at that k-point is rebuilt from parity, 0 if stored
  !!@param   Ndk,out: number of stored k-points
  !!@param   invk,in: list of reverse k-points
  !!@param  Nkall,in: Number of all k-points
  !!@param     Nk,in: Number of irreducible k-points
  !!@param gap_sym,in: gap symmetry number
  use,intrinsic:: iso_c_binding, only:c_int64_t
  implicit none
  integer(c_int64_t),intent(in):: Nkall,Nk,gap_sym
  integer(c_int64_t),intent(in),dimension(3,Nkall):: invk
  integer(c_int64_t),intent(out),dimension(Nkall):: idk,ipk
  integer(c_int64_t),intent(out):: Ndk

  integer(c_int64_t),parameter:: gap_unit=0,par_mix=-10
  integer(c_int64_t) i

  if(gap_sym/=gap_unit .and. gap_sym/=par_mix)then
     Ndk=Nk
     idk(:)=invk(1,:)
     ipk(:)=invk(2,:)
  else
     Ndk=Nkall
     idk(:)=[(i,i=1,Nkall)]
     ipk(:)=0
  end if
end subroutine get_delta_kmap_soc

subroutine conv_delta_orb_to_band_soc(deltab,delta,uni,invk,invs,slist,Norb,Nkall,Nk,Nw,gap_sym) bind(C)
  !> This function transform orbital basis gap function into band basis gap function
  !> Band n at k is paired with its Kramers partner (inverse of get_initial_delta_soc), so the result is
  !> gauge covariant for a (pseudo)spin singlet.  A triplet (gap_sym<0, except parity mixing -10) is
  !> first brought to singlet form by sigma_z on the row index (sigma_z sigma_z i sigma_y = i sigma_y);
  !> the output is then its d_z (Sz=0) component in the band basis.  No extra sign on the -k half
  !> (unlike the seed): the -k half then shows d_z(-k) = -d_z(k), i.e. the physical odd parity.
  !!@param deltab,out: band basis gap function
  !!@param   delta,in: orbital basis gap function
  !!@param     uni,in: unitary matrix
  !!@param    invk,in: footnote of reverse k
  !!@param    Norb,in: Number of orbital
  !!@param   Nkall,in: Number of all k-points
  !!@param      Nk,in: Number of irreducible k-points
  !!oparam      Nw,in: Number of Matsubara frequencies
  !!@param gap_sym,in: gap symmetry number (sign selects singlet/triplet)
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nw,Norb,Nk,Nkall,gap_sym
  integer(c_int64_t),intent(in),dimension(Norb):: slist,invs
  integer(c_int64_t),intent(in),dimension(3,Nkall):: invk
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(c_double),intent(in),dimension(Nkall,Nw,Norb,Norb):: delta
  complex(c_double),intent(out),dimension(Nkall,Norb,Norb):: deltab

  integer(c_int32_t) i,k,l,n
  integer(c_int64_t),parameter:: par_mix=-10
  complex(c_double),dimension(Norb,Norb):: matL,matR,tmpm,dtmp

  !$omp parallel do private(i,k,l,n,matL,matR,tmpm,dtmp)
  do i=1,Nkall
     k=invk(1,i)
     if(invk(2,i)==0)then
        matL=conjg(transpose(uni(:,:,k)))
        do n=1,Norb
           matR(n,:)=slist(n)*uni(invs(n),:,k)
        end do
     else if(invk(2,i)==1)then
        do n=1,Norb
           matL(:,n)=-slist(n)*uni(invs(n),:,k)
        end do
        matR=conjg(uni(:,:,k))
     else
        deltab(i,:,:)=0.0d0
        cycle
     end if

     dtmp(:,:)=delta(i,1,:,:)
     if(gap_sym<0 .and. gap_sym/=par_mix)then !triplet d_z -> singlet form: sigma_z on the row index
        do n=1,Norb
           dtmp(n,:)=slist(n)*dtmp(n,:)
        end do
     end if
     call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),matL,Norb,dtmp,Norb,(0.0d0,0.0d0),tmpm,Norb)
     call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),tmpm,Norb,matR,Norb,(0.0d0,0.0d0),deltab(i,:,:),Norb)
  end do
  !$omp end parallel do
end subroutine conv_delta_orb_to_band_soc
