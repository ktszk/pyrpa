module constant
  use,intrinsic:: iso_c_binding, only:c_double
  implicit none
  real(c_double),parameter:: pi=3.141592653589793238462643383279d0
end module constant

module numerics
  use,intrinsic:: iso_c_binding, only:c_double,c_int64_t
  implicit none
contains
  subroutine drop_roundoff(z,n)
    !> Zero the real/imaginary parts of z that are pure FFT round-off.
    !>
    !> The cut is RELATIVE to the largest element: the convolution carries a
    !> round-off floor of ~1e-16*max|z|, so 1e-13*max|z| removes it with three
    !> decades of margin and touches nothing else.  Every call site used to use
    !> a fixed ABSOLUTE 1e-9 instead, which is silently scale dependent: for a
    !> small U, a low T or a weakly polarizable system - where chi or sigma is
    !> itself O(1e-8..1e-6) - it deleted physical structure rather than noise.
    !>
    !> n is passed explicitly so a contiguous slice (e.g. chi(:,:,m,l)) can be
    !> handed over by sequence association.
    integer(c_int64_t),intent(in):: n
    complex(c_double),intent(inout):: z(n)
    real(c_double) eps,amax
    integer(c_int64_t) i
    amax=maxval(abs(z))
    if(amax<=0.0d0)return
    eps=1.0d-13*amax
    !$omp parallel do private(i)
    do i=1,n
       if(abs(dble(z(i)))<eps) z(i)=cmplx(0.0d0,imag(z(i)),kind=c_double)
       if(abs(imag(z(i)))<eps) z(i)=cmplx(dble(z(i)),0.0d0,kind=c_double)
    end do
    !$omp end parallel do
  end subroutine drop_roundoff

  subroutine reverse_index(ii,ij,ik,iw,Nx,Ny,Nz,Nw)
    !> Index tables for the reflections used by every FFT convolution here:
    !> ii/ij/ik send a real-space mesh point to -r (0 -> 0, i -> N-i) and iw sends a Matsubara
    !> index on the 2*Nw circular grid to -tau / -v_n (1 -> 1, i -> 2*Nw-i+2).
    !> Same table in fchi, fchi0withsoc, ffeliash and f_eliash_nosym, hence shared.
    use,intrinsic:: iso_c_binding, only:c_int64_t,c_int32_t
    integer(c_int64_t),intent(in):: Nx,Ny,Nz,Nw
    integer(c_int32_t),intent(out):: ii(0:Nx-1),ij(0:Ny-1),ik(0:Nz-1),iw(2*Nw)
    integer(c_int32_t) i
    ii(0)=0
    ij(0)=0
    ik(0)=0
    iw(1)=1
    do i=1,int(Nx,c_int32_t)-1
       ii(i)=int(Nx,c_int32_t)-i
    end do
    do i=1,int(Ny,c_int32_t)-1
       ij(i)=int(Ny,c_int32_t)-i
    end do
    do i=1,int(Nz,c_int32_t)-1
       ik(i)=int(Nz,c_int32_t)-i
    end do
    do i=2,2*int(Nw,c_int32_t)
       iw(i)=2*int(Nw,c_int32_t)-i+2
    end do
  end subroutine reverse_index
end module numerics

module eliash_krylov
  !> Krylov machinery shared by the two TIME-REVERSAL-SYMMETRIC linearized Eliashberg solvers,
  !> lin_eliash (feliash.f90, no SOC) and lin_eliash_soc (feliashwithsoc.f90).
  !>
  !> Both build only part of the gap function and fill the rest by conjugation (Delta is Hermitian
  !> in the TRS gauge; with SOC also Delta_uu = -Delta_dd^+), so their kernels are linear over the
  !> REALS only: K(i*Delta) /= i*K(Delta).  Everything here is therefore built on a REAL inner
  !> product (kv_dot) with REAL Ritz coefficients (real_ritz).  A complex Krylov space silently
  !> mixes the two halves and returns wrong eigenvalues, which is why this is shared rather than
  !> written out twice.  f_eliash_nosym.f90 has NO time-reversal symmetry, its kernel is genuinely
  !> complex-linear, and it keeps its own complex Arnoldi -- do not route it through this module.
  !>
  !> A gap-function vector is an array (nd,Nw,Norb,Norb) of stored k-points; which k-points are
  !> stored differs between the solvers (irreducible list, or the full BZ), so the geometry is
  !> carried by a kspace and the arrays are taken as flat assumed-size dummies.
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  private
  public:: kspace,kspace_init,kv_norm,kv_dot,kv_err,real_ritz,arnoldi_two_pass
  public:: eliash_op_if,eliash_class_if

  type:: kspace
     !> Shape of one gap-function vector, plus how its stored k-points cover the full BZ.
     integer(c_int32_t):: nd=0     !< number of stored k-points (first array dimension)
     integer(c_int32_t):: nrest=0  !< Nw*Norb*Norb, i.e. everything after the k index
     integer(c_int32_t):: nvec=0   !< nd*nrest, the length of one vector
     !> wk(p) = how many full-BZ points are represented by stored point p.  Sums and norms must
     !> run over the whole BZ even though only nd points are kept, so every term is weighted
     !> with wk; a solver that stores the full BZ simply has wk = 1 everywhere.
     real(c_double),allocatable:: wk(:)
  end type kspace

  abstract interface
     subroutine eliash_op_if(w_out,v_in)
       !> one application of the linearized Eliashberg kernel, w_out = (T/Nkall) K v_in
       import:: c_double
       complex(c_double),intent(in):: v_in(*)
       complex(c_double),intent(out):: w_out(*)
     end subroutine eliash_op_if

     subroutine eliash_class_if(ok,x)
       !> does the eigenvector x belong to the solution branch the caller asked for?
       import:: c_double
       complex(c_double),intent(in):: x(*)
       logical,intent(out):: ok
     end subroutine eliash_class_if
  end interface
contains

  subroutine kspace_init(ks,idk,nkall,nd,nw,norb)
    !> Set up the vector geometry.  idk(i) is the stored index of full-BZ point i, so counting how
    !> often each stored index appears gives the multiplicity wk.
    type(kspace),intent(out):: ks
    integer(c_int64_t),intent(in):: nkall,nd,nw,norb
    integer(c_int64_t),intent(in):: idk(nkall)
    integer(c_int64_t) i
    ks%nd=int(nd,c_int32_t)
    ks%nrest=int(nw*norb*norb,c_int32_t)
    ks%nvec=ks%nd*ks%nrest
    allocate(ks%wk(ks%nd))
    ks%wk(:)=0.0d0
    do i=1,nkall
       ks%wk(idk(i))=ks%wk(idk(i))+1.0d0
    end do
  end subroutine kspace_init

  real(c_double) function kv_dot(ks,u,v)
    !> Real inner product 2*Re<u,v>, summed over the full BZ (weight wk) and over the stored
    !> w_n>0 half; the factor 2 accounts for the w_n<0 half that is not stored.
    !> Real, not complex: see the note at the top of the module.
    type(kspace),intent(in):: ks
    complex(c_double),intent(in):: u(*),v(*)
    integer(c_int32_t) r,p,o
    real(c_double) s
    s=0.0d0
    !$omp parallel do private(r,p,o) reduction(+:s)
    do r=1,ks%nrest
       o=(r-1)*ks%nd            ! start of the (w_n,orbital,orbital) block r
       do p=1,ks%nd
          s=s+ks%wk(p)*dble(conjg(u(o+p))*v(o+p))
       end do
    end do
    !$omp end parallel do
    kv_dot=2.0d0*s
  end function kv_dot

  real(c_double) function kv_norm(ks,x)
    !> norm induced by kv_dot
    type(kspace),intent(in):: ks
    complex(c_double),intent(in):: x(*)
    kv_norm=sqrt(max(kv_dot(ks,x,x),0.0d0))
  end function kv_norm

  real(c_double) function kv_err(ks,xnew,xold,inrm)
    !> Vector convergence of the power iteration: min(||xnew*inrm - xold||, ||xnew*inrm + xold||).
    !> The minimum of the two signs handles an eigenvector that flips sign every iteration, which
    !> is what a negative eigenvalue does.
    type(kspace),intent(in):: ks
    complex(c_double),intent(in):: xnew(*),xold(*)
    real(c_double),intent(in):: inrm
    integer(c_int32_t) r,p,o
    real(c_double) sp,sm
    sp=0.0d0
    sm=0.0d0
    !$omp parallel do private(r,p,o) reduction(+:sp,sm)
    do r=1,ks%nrest
       o=(r-1)*ks%nd
       do p=1,ks%nd
          sp=sp+ks%wk(p)*abs(xnew(o+p)*inrm-xold(o+p))**2
          sm=sm+ks%wk(p)*abs(xnew(o+p)*inrm+xold(o+p))**2
       end do
    end do
    !$omp end parallel do
    kv_err=sqrt(2.0d0*min(sp,sm))
  end function kv_err

  subroutine real_ritz(y,n)
    !> Rotate the Ritz coefficients y of a real Hessenberg matrix to be real, which for a real
    !> Ritz value only removes the arbitrary phase ZGEEV leaves behind.  Necessary because the
    !> operator is real-linear: a complex y would combine the Krylov vectors in a way the operator
    !> does not reproduce, so the "eigenvector" would not be one.
    integer(c_int32_t),intent(in):: n
    complex(c_double),intent(inout):: y(n)
    complex(c_double) ph
    ph=y(maxloc(abs(y),1))
    y(:)=cmplx(dble(y(:)*conjg(ph)/abs(ph)),0.0d0,c_double)
  end subroutine real_ritz

  subroutine ritz_vector(x,V,y,ks,m_act)
    !> x = sum_ii y(ii) * v_{ii-1}, the Ritz vector in the original space
    type(kspace),intent(in):: ks
    integer(c_int32_t),intent(in):: m_act
    complex(c_double),intent(in):: V(ks%nvec,0:*),y(*)
    complex(c_double),intent(out):: x(*)
    integer(c_int32_t) ii,p
    !$omp parallel do private(p)
    do p=1,ks%nvec
       x(p)=(0.0d0,0.0d0)
    end do
    !$omp end parallel do
    do ii=1,m_act
       !$omp parallel do private(p)
       do p=1,ks%nvec
          x(p)=x(p)+y(ii)*V(p,ii-1)
       end do
       !$omp end parallel do
    end do
  end subroutine ritz_vector

  subroutine pick_branch(idx,x,eigenvals,VR,ldvr,V,ks,m_act,shift,classify)
    !> Choose which Ritz pair to return and leave its (normalized) Ritz vector in x.
    !>
    !> Without a classify callback this is simply the largest physical eigenvalue.  With one, it
    !> is the largest eigenvalue whose eigenvector the caller accepts: with SOC a parity sector
    !> holds both the (pseudospin singlet, even w) and the (pseudospin triplet, odd w) pairing and
    !> no exact symmetry separates them, so the largest eigenvalue of the sector need not be the
    !> branch that was asked for, and lambda would jump as the SOC strength goes to zero.
    !>
    !> idx=0 means no Ritz value is positive at all; the caller decides what to do.
    type(kspace),intent(in):: ks
    integer(c_int32_t),intent(out):: idx
    integer(c_int32_t),intent(in):: m_act
    integer(c_int32_t),intent(in):: ldvr   !< leading dimension of VR as the caller allocated it;
                                           !< m_act can be smaller (early stop), so it must be
                                           !< passed rather than assumed equal to m_act
    complex(c_double),intent(in):: eigenvals(*),VR(ldvr,*),V(ks%nvec,0:*)
    complex(c_double),intent(out):: x(*)
    real(c_double),intent(in):: shift
    procedure(eliash_class_if),optional:: classify
    integer(c_int32_t),parameter:: ncand=50   ! Ritz vectors built at most; each is cheap next to
                                              ! one apply_op, but not free
    integer(c_int32_t) ii,jj,n,itry,itmp
    integer(c_int32_t),dimension(m_act):: iord
    complex(c_double),dimension(m_act):: y    ! Ritz coefficients of the candidate, phase fixed
    logical ok

    ! candidates: a physical eigenvalue must be positive (lambda<0 is a repulsive solution)
    n=0
    do ii=1,m_act
       if(dble(eigenvals(ii))-shift>0.0d0)then; n=n+1; iord(n)=ii; end if
    end do
    idx=0
    if(n==0)return
    ! order them by descending eigenvalue (selection sort: m_act is at most a few hundred)
    do ii=1,n-1
       do jj=ii+1,n
          if(dble(eigenvals(iord(jj)))>dble(eigenvals(iord(ii))))then
             itmp=iord(ii); iord(ii)=iord(jj); iord(jj)=itmp
          end if
       end do
    end do

    ! the Ritz coefficients must be made real before the Ritz vector is formed: with a
    ! real-linear operator a complex combination of the Krylov vectors is not an eigenvector
    do itry=1,merge(min(n,ncand),1,present(classify))
       idx=iord(itry)
       y(1:m_act)=VR(1:m_act,idx)
       call real_ritz(y,m_act)
       call ritz_vector(x,V,y,ks,m_act)
       if(.not.present(classify))exit
       call classify(ok,x)
       if(ok)then
          if(itry>1)print'(A,I0,A,F12.6)','  branch: took Ritz ',itry,' of the sector, below ',&
               dble(eigenvals(iord(1)))-shift
          exit
       end if
       if(itry==min(n,ncand))then
          print'(A)','  WARNING(pick_branch): no Ritz vector of the requested branch; taking the largest'
          idx=iord(1)
          y(1:m_act)=VR(1:m_act,idx)
          call real_ritz(y,m_act)
          call ritz_vector(x,V,y,ks,m_act)
       end if
    end do
  end subroutine pick_branch

  subroutine arnoldi_two_pass(lambda,x,seed,ks,m_arnoldi,tol,apply_op,classify)
    !> Two-pass Arnoldi for the physical (largest positive) Eliashberg eigenvalue.
    !>
    !> The kernel's dominant eigenvalue is a large NEGATIVE one, so a plain Arnoldi converges to
    !> the wrong end of the spectrum.  Pass 1 therefore finds lambda_- (most negative) and its
    !> eigenvector; pass 2 runs (A + |lambda_-|) deflated against that eigenvector, which makes
    !> lambda_+ the dominant eigenvalue of the shifted operator.
    !>
    !> Arnoldi factorization used in both passes: with an orthonormal basis V_m,
    !> A V_m = V_m H_m + v_m h_{m+1}^T, and the eigenpairs (theta, y) of H_m give the Ritz pairs
    !> (theta, V_m y) with residual |h_{m+1}.y|.  All inner products are real (kv_dot) because
    !> the operator is only real-linear -- see the note at the top of the module.
    !>
    !> Restarts (Krylov-Schur): one Krylov space of fixed size m is not enough in general.  For a
    !> small lambda_+ next to a large |lambda_-| an m=10 space left lambda_+ wrong by ~1%, and
    !> two physical eigenvalues 1% apart cannot be told apart by restarting from a single Ritz
    !> vector (it converges to whichever dominates that vector, not to the largest).  So each
    !> pass is restarted THICK: H_m is brought to real Schur form with the wanted end of its
    !> spectrum (most negative for pass 1, most positive for pass 2) in the leading k ~ m/2 block,
    !> those k Schur vectors are kept, and the space is expanded again to m, until the residual
    !> of the wanted Ritz pair is below tol*|lambda|.  The memory is m+1+k vectors.
    !!@param  lambda,out: physical eigenvalue
    !!@param       x,out: its eigenvector, normalized
    !!@param      seed,in: starting vector (v_0), need not be normalized
    !!@param        ks,in: vector geometry
    !!@param m_arnoldi,in: Krylov dimension
    !!@param       tol,in: relative residual at which a Ritz pair counts as converged
    !!@param  apply_op,in: one application of the kernel
    !!@param  classify,in: optional branch test, see pick_branch
    type(kspace),intent(in):: ks
    integer(c_int64_t),intent(in):: m_arnoldi
    real(c_double),intent(in):: tol
    real(c_double),intent(out):: lambda
    complex(c_double),intent(out):: x(*)
    complex(c_double),intent(in):: seed(*)
    procedure(eliash_op_if):: apply_op
    procedure(eliash_class_if),optional:: classify

    integer(c_int32_t) m_dim,p,m_act,idx,lwork,info,irs,j0
    integer(c_int32_t),parameter:: max_restart=200
    real(c_double) beta,shift,lambda1,hij,res,thr,keep_ar,keep_tol
    logical want_max,sw_keep
    complex(c_double),allocatable:: V(:,:),w(:),x1(:),yd(:),vres(:),Vk(:,:)
    complex(c_double),allocatable:: H_mat(:,:),A_mat(:,:),eigenvals(:)
    complex(c_double),allocatable:: VL_dum(:,:),VR(:,:),zwork(:)
    real(c_double),allocatable:: rwork(:)

    m_dim=max(int(m_arnoldi,c_int32_t),2)
    lwork=max(64,2*m_dim)
    allocate(V(ks%nvec,0:m_dim-1),w(ks%nvec),x1(ks%nvec),yd(m_dim),vres(ks%nvec),Vk(ks%nvec,m_dim))
    allocate(H_mat(m_dim+1,m_dim),A_mat(m_dim,m_dim),eigenvals(m_dim))
    allocate(VL_dum(1,1),VR(m_dim,m_dim),zwork(lwork),rwork(2*m_dim))

    ! ===== Pass 1: no shift, find lambda_- ==================================
    print'(A,I4,A)','Arnoldi pass 1 (no shift): Krylov dim =',m_dim,' ...'
    call set_v0(seed)
    j0=0
    want_max=.false.
    sw_keep=.false.
    pass1: do irs=1,max_restart
       call arnoldi_loop(0.0d0,.false.,j0)
       call diagonalize('pass 1')
       if(irs==1)then
          ! early exit: if the first space already resolved the wanted branch at a clearly
          ! positive, CONVERGED Ritz value, that is lambda_+ and the shifted pass is unnecessary
          call pick_branch(idx,x,eigenvals,VR,m_dim,V,ks,m_act,0.0d0,classify)
          if(idx>0)then
             if(dble(eigenvals(idx))>=0.1d0 .and. &
                  ritz_res(idx)<=tol*abs(dble(eigenvals(idx))))then
                lambda=dble(eigenvals(idx))
                print'(A,F12.6)','  eliash   =',lambda
                call normalize(x)
                return
             end if
          end if
       end if
       idx=minloc(dble(eigenvals(1:m_act)),1)     ! most negative Ritz value
       lambda1=dble(eigenvals(idx))
       res=ritz_res(idx)
       if(res<=tol*abs(lambda1))exit pass1
       if(irs==max_restart)then
          print'(A,ES10.3)','  WARNING(arnoldi): lambda_- not converged, relative residual',res/abs(lambda1)
          exit pass1
       end if
       call thick_restart(j0)
    end do pass1
    ! deflation direction x1 = Ritz vector of lambda_-
    yd(1:m_act)=VR(1:m_act,idx)
    call real_ritz(yd,m_act)
    call ritz_vector(x1,V,yd,ks,m_act)
    call normalize(x1)
    print'(A,F12.6,A,I0,A,ES9.2,A)','  lambda_- =',lambda1,'  (',irs,' Krylov spaces, residual ',&
         res/abs(lambda1),')'

    ! ===== Pass 2: shift + deflation, find lambda_+ =========================
    ! B = (A + shift) restricted to the complement of x1; eigenvalue lambda_i -> lambda_i + shift,
    ! so with shift = -lambda_- the wanted lambda_+ becomes the dominant eigenvalue of B.
    shift=-lambda1
    print'(A,F10.6,A,I4,A)','Arnoldi pass 2 (shift=',shift,'): Krylov dim =',m_dim,' ...'
    call set_v0(seed)
    call deflate(V(1,0))
    call normalize(V(1,0))
    j0=0
    want_max=.true.
    pass2: do irs=1,max_restart
       call arnoldi_loop(shift,.true.,j0)
       call diagonalize('pass 2')
       call pick_branch(idx,x,eigenvals,VR,m_dim,V,ks,m_act,shift,classify)
       if(idx==0)then
          idx=maxloc(dble(eigenvals(1:m_act)),1)  ! nothing positive: largest Ritz value
          yd(1:m_act)=VR(1:m_act,idx)
          call real_ritz(yd,m_act)
          call ritz_vector(x,V,yd,ks,m_act)
       end if
       lambda=dble(eigenvals(idx))-shift
       call normalize(x)
       ! the residual is that of the shifted operator, but it bounds the error of lambda_+ itself,
       ! so it is measured against |lambda_+| (floored, for a lambda_+ near zero)
       res=ritz_res(idx)
       if(res<=tol*max(abs(lambda),1.0d-3))exit pass2
       if(irs==max_restart)then
          print'(A,ES10.3)','  WARNING(arnoldi): lambda_+ not converged, relative residual',&
               res/max(abs(lambda),1.0d-3)
          exit pass2
       end if
       ! the restart must keep the Ritz pair that was picked: with classify it can lie below the
       ! top half of the spectrum, which the threshold alone would throw away at every restart
       sw_keep=.true.
       keep_ar=dble(eigenvals(idx))
       call thick_restart(j0)
    end do pass2
    print'(A,F12.6,A,I0,A,ES9.2,A)','  eliash   =',lambda,'  (',irs,' Krylov spaces, residual ',&
         res/max(abs(lambda),1.0d-3),')'
    ! x solves the DEFLATED operator P A P (P = 1 - x1 x1^T), not A itself:
    ! A x = lambda x + c x1 with c = <x1, A x>, so the eigenvector of A is x + c/(lambda_+ - lambda_-) x1
    if(abs(lambda-lambda1)>1.0d-12)then    ! lambda_+ = lambda_- (e.g. V=0): nothing to undo
       call apply_op(w,x)
       hij=kv_dot(ks,x1,w)/(lambda-lambda1)
       !$omp parallel do private(p)
       do p=1,ks%nvec
          x(p)=x(p)+hij*x1(p)
       end do
       !$omp end parallel do
       call normalize(x)
    end if
    deallocate(V,w,x1,yd,vres,Vk,H_mat,A_mat,eigenvals,VL_dum,VR,zwork,rwork)
  contains
    real(c_double) function ritz_res(i)
      !> residual norm of Ritz pair i of the current space, |h_{m+1,m}| |y_m| / ||y|| (the last
      !> row of H has only that element after an expansion step); 0 after a lucky breakdown,
      !> where the space is invariant and the Ritz pair exact
      integer(c_int32_t),intent(in):: i
      if(m_act<m_dim)then
         ritz_res=0.0d0
      else
         ritz_res=abs(H_mat(m_act+1,m_act))*abs(VR(m_act,i))/max(sqrt(sum(abs(VR(1:m_act,i))**2)),1.0d-300)
      end if
    end function ritz_res

    logical function sel_wanted(ar,ai)
      !> dgees selection: the wanted end of the spectrum (threshold thr set by thick_restart);
      !> ai is not needed, the two members of a complex pair share ar
      real(c_double),intent(in):: ar,ai
      if(want_max)then
         sel_wanted=(ar>=thr)
      else
         sel_wanted=(ar<=thr)
      end if
      if(sw_keep)sel_wanted=sel_wanted .or. abs(ar-keep_ar)<=keep_tol
    end function sel_wanted

    subroutine thick_restart(jnext)
      !> Krylov-Schur restart: real Schur form H = Q T Q^T with the wanted eigenvalues in the
      !> leading k x k block; keep V Q(:,1:k), put the old residual direction next, and continue
      !> the expansion from column k.  A V_k = V_k T_k + v_m b^T with b = h_{m+1,m} Q(m,1:k).
      integer(c_int32_t),intent(out):: jnext
      integer(c_int32_t) k,sdim,linfo,lw,ii
      real(c_double),allocatable:: Hr(:,:),Q(:,:),wr(:),wi(:),work(:),srt(:)
      logical,allocatable:: bwork(:)
      complex(c_double),allocatable:: Qc(:,:)
      real(c_double) hm1
      allocate(Hr(m_dim,m_dim),Q(m_dim,m_dim),wr(m_dim),wi(m_dim),srt(m_dim),bwork(m_dim),Qc(m_dim,m_dim))
      lw=max(64,6*m_dim)
      allocate(work(lw))
      hm1=dble(H_mat(m_dim+1,m_dim))
      Hr(:,:)=dble(H_mat(1:m_dim,1:m_dim))
      ! threshold: the k-th value from the wanted end, k ~ m/2 (at least 1, at most m-2)
      k=min(max(m_dim/2,1),m_dim-2)
      srt(:)=dble(eigenvals(1:m_dim))
      call sort_real(srt,m_dim)                      ! ascending
      if(want_max)then
         thr=srt(m_dim-k+1)-1.0d-12*max(abs(srt(1)),abs(srt(m_dim)))
      else
         thr=srt(k)+1.0d-12*max(abs(srt(1)),abs(srt(m_dim)))
      end if
      ! eigenvalues of the real H from dgees vs those of zgeev: equal up to rounding
      keep_tol=1.0d-8*max(abs(srt(1)),abs(srt(m_dim)),1.0d-300)
      call dgees('V','S',sel_wanted,m_dim,Hr,m_dim,sdim,wr,wi,Q,m_dim,work,lw,bwork,linfo)
      if(linfo/=0)then; print*,'DGEES failed in thick_restart: info=',linfo; stop; end if
      k=max(sdim,1)
      if(k<m_dim .and. abs(Hr(min(k+1,m_dim),k))>0.0d0)k=k+1   ! do not split a 2x2 block
      if(k>=m_dim)then
         ! at most m-1 vectors are kept (one new direction is needed); step back over a whole
         ! 2x2 block rather than splitting it, which would break A V_k = V_k T_k + v b^T
         k=m_dim-1
         if(k>1 .and. abs(Hr(k+1,k))>0.0d0)k=k-1
      end if
      Qc(:,1:k)=cmplx(Q(:,1:k),0.0d0,c_double)
      call zgemm('N','N',ks%nvec,k,m_dim,(1.0d0,0.0d0),V,ks%nvec,Qc,m_dim,(0.0d0,0.0d0),Vk,ks%nvec)
      do ii=1,k
         !$omp parallel do private(p)
         do p=1,ks%nvec
            V(p,ii-1)=Vk(p,ii)
         end do
         !$omp end parallel do
      end do
      !$omp parallel do private(p)
      do p=1,ks%nvec
         V(p,k)=vres(p)
      end do
      !$omp end parallel do
      H_mat(:,:)=(0.0d0,0.0d0)
      H_mat(1:k,1:k)=cmplx(Hr(1:k,1:k),0.0d0,c_double)
      H_mat(k+1,1:k)=cmplx(hm1*Q(m_dim,1:k),0.0d0,c_double)
      jnext=k
      deallocate(Hr,Q,wr,wi,srt,bwork,Qc,work)
    end subroutine thick_restart

    subroutine sort_real(a,n)
      !> ascending insertion sort (n is a Krylov dimension, tens at most)
      integer(c_int32_t),intent(in):: n
      real(c_double),intent(inout):: a(n)
      integer(c_int32_t) ii,jj
      real(c_double) t
      do ii=2,n
         t=a(ii)
         jj=ii-1
         do while(jj>=1)
            if(a(jj)<=t)exit
            a(jj+1)=a(jj)
            jj=jj-1
         end do
         a(jj+1)=t
      end do
    end subroutine sort_real

    subroutine set_v0(s)
      !> v_0 = normalized seed
      complex(c_double),intent(in):: s(*)
      integer(c_int32_t) q
      !$omp parallel do private(q)
      do q=1,ks%nvec
         V(q,0)=s(q)
      end do
      !$omp end parallel do
      call normalize(V(1,0))
    end subroutine set_v0

    subroutine normalize(y)
      complex(c_double),intent(inout):: y(*)
      integer(c_int32_t) q
      real(c_double) nrm
      nrm=kv_norm(ks,y)
      if(nrm<=0.0d0)then; print*,'Error(arnoldi_two_pass): zero vector'; stop; end if
      !$omp parallel do private(q)
      do q=1,ks%nvec
         y(q)=y(q)*(1.0d0/nrm)
      end do
      !$omp end parallel do
    end subroutine normalize

    subroutine deflate(y)
      !> remove the x1 component, so the iteration stays in the complement of lambda_-
      complex(c_double),intent(inout):: y(*)
      integer(c_int32_t) q
      real(c_double) c
      c=kv_dot(ks,x1,y)
      !$omp parallel do private(q)
      do q=1,ks%nvec
         y(q)=y(q)-c*x1(q)
      end do
      !$omp end parallel do
    end subroutine deflate

    subroutine arnoldi_loop(shft,sw_deflate,jstart)
      !> expand the basis from column jstart to m: w = (A + shft) v_j, optionally deflated, then
      !> Gram-Schmidt against v_0..v_j applied twice (a thick restart keeps old vectors, whose
      !> orthogonality a single sweep would slowly lose).  The last new direction is kept in
      !> vres for the next restart.  jstart=0 starts a fresh factorization.
      real(c_double),intent(in):: shft
      logical,intent(in):: sw_deflate
      integer(c_int32_t),intent(in):: jstart
      integer(c_int32_t) jj,kk,q,isw
      real(c_double) hjj
      if(jstart==0)H_mat=(0.0d0,0.0d0)
      m_act=m_dim
      do jj=jstart,m_dim-1
         call apply_op(w,V(1,jj))
         if(shft/=0.0d0)then
            !$omp parallel do private(q)
            do q=1,ks%nvec
               w(q)=w(q)+shft*V(q,jj)
            end do
            !$omp end parallel do
         end if
         if(sw_deflate)call deflate(w)
         do isw=1,2
            do kk=0,jj
               hjj=kv_dot(ks,V(1,kk),w)
               H_mat(kk+1,jj+1)=H_mat(kk+1,jj+1)+cmplx(hjj,0.0d0,c_double)
               !$omp parallel do private(q)
               do q=1,ks%nvec
                  w(q)=w(q)-hjj*V(q,kk)
               end do
               !$omp end parallel do
            end do
         end do
         beta=kv_norm(ks,w)
         H_mat(jj+2,jj+1)=cmplx(beta,0.0d0,c_double)
         if(beta<1.0d-14)then                   ! lucky breakdown: the Krylov space is invariant
            m_act=jj+1
            print'(A,I5)','  Krylov space closed early at dim =',m_act
            exit
         end if
         if(jj<m_dim-1)then
            !$omp parallel do private(q)
            do q=1,ks%nvec
               V(q,jj+1)=w(q)*(1.0d0/beta)
            end do
            !$omp end parallel do
         else
            !$omp parallel do private(q)
            do q=1,ks%nvec
               vres(q)=w(q)*(1.0d0/beta)
            end do
            !$omp end parallel do
         end if
      end do
    end subroutine arnoldi_loop

    subroutine diagonalize(tag)
      !> Ritz values and coefficients of the m_act x m_act projected matrix
      character(*),intent(in):: tag
      integer(c_int32_t) kk
      A_mat(1:m_act,1:m_act)=H_mat(1:m_act,1:m_act)   ! ZGEEV destroys its input
      call zgeev('N','V',m_act,A_mat,m_dim,eigenvals,VL_dum,1,VR,m_dim,zwork,lwork,rwork,info)
      if(info/=0)then; print*,'ZGEEV failed ',tag,': info=',info; stop; end if
      print'(A,A,A)','Ritz values (',tag,', real):'
      do kk=1,m_act
         print'(I4,F14.6)',kk,dble(eigenvals(kk))
      end do
    end subroutine diagonalize
  end subroutine arnoldi_two_pass
end module eliash_krylov

module flex_scf
  !> Pieces of the FLEX self-consistency loop shared by mkself (fself.f90, no SOC) and
  !> mkself_soc (fselfwithsoc.f90).  The two drivers differ in how G, chi and Sigma are built,
  !> but the chemical-potential search around them is the same algorithm, so it lives here
  !> instead of being copied (it was, and the copies had already started to drift).
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  private
  public:: flex_rn_if,solve_mu
  public:: diis_state,diis_init,diis_free,mix_sigma

  abstract interface
     subroutine flex_rn_if(rn,rmu)
       !> electron number per k-point at chemical potential rmu, for the current Sigma
       import:: c_double
       real(c_double),intent(out):: rn
       real(c_double),intent(in):: rmu
     end subroutine flex_rn_if
  end interface

  type:: diis_state
     !> Pulay (DIIS) history for a self-consistency loop on a (Nk,Nw,Norb,Norb) field:
     !> the last m outputs and their residuals, in a circular buffer.
     integer(c_int32_t):: m=0        !< history length
     integer(c_int32_t):: i_hist=0   !< slot written last
     integer(c_int32_t):: n_hist=0   !< slots filled so far (<= m)
     complex(c_double),allocatable:: xout(:,:,:,:,:),res(:,:,:,:,:)
     real(c_double),allocatable:: B(:,:),rhs(:)
     integer(c_int32_t),allocatable:: ipiv(:)
  end type diis_state

contains

  subroutine diis_init(st,m,Nk,Nw,Norb)
    type(diis_state),intent(out):: st
    integer(c_int64_t),intent(in):: m,Nk,Nw,Norb
    st%m=int(m,c_int32_t)
    st%i_hist=0
    st%n_hist=0
    allocate(st%xout(Nk,Nw,Norb,Norb,m),st%res(Nk,Nw,Norb,Norb,m))
    allocate(st%B(m+1,m+1),st%rhs(m+1),st%ipiv(m+1))
  end subroutine diis_init

  subroutine diis_free(st)
    type(diis_state),intent(inout):: st
    if(allocated(st%xout))deallocate(st%xout,st%res,st%B,st%rhs,st%ipiv)
    st%m=0
    st%i_hist=0
    st%n_hist=0
  end subroutine diis_free

  subroutine mix_sigma(st,sigmak,sigmak0,esterr,pp,eps_sgm,invk,kmap,Nk,Nkall,Nw,Norb)
    !> One DIIS (Pulay) update of the FLEX self-energy, plus the convergence estimate.
    !>
    !> Pulay: with residuals r_i = Sigma_out_i - Sigma_in_i, solve [B 1; 1^T 0][c; lambda] = [0; 1]
    !> for B_ij = Re<r_i,r_j> and extrapolate Sigma = sum_i c_i Sigma_out_i.  The first step (or
    !> m_diis=1) has nothing to extrapolate from and falls back to linear mixing with rate pp,
    !> and so does a singular B (dgesv failure), which weights the most recent entry only.
    !>
    !> esterr is the largest relative change of any Sigma element above eps_sgm; elements below
    !> that are skipped so that near-zero components do not dominate the convergence test.
    !!@param          st,inout: DIIS history
    !!@param      sigmak,inout: in: Sigma_out of this iteration; out: the mixed Sigma
    !!@param       sigmak0,in: Sigma_in of this iteration
    !!@param        esterr,out: largest relative change, i.e. the convergence measure
    !!@param            pp,in: linear-mixing rate used when DIIS cannot extrapolate yet
    !!@param       eps_sgm,in: ignore elements smaller than this in esterr
    !!@param invk,kmap,in: k-point tables, used only to label the worst element in the printout
    type(diis_state),intent(inout):: st
    integer(c_int64_t),intent(in):: Nk,Nkall,Nw,Norb
    integer(c_int64_t),intent(in),dimension(3,Nkall):: invk,kmap
    complex(c_double),intent(inout),dimension(Nk,Nw,Norb,Norb):: sigmak
    complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: sigmak0
    real(c_double),intent(out):: esterr
    real(c_double),intent(in):: pp,eps_sgm

    integer(c_int32_t) i,j,l,m,kerr,iwerr,lerr,merr,ih,jh,idx_i,idx_j,info,n_cur
    real(c_double) est,eps_reg
    complex(c_double) cksigm
    complex(c_double),dimension(Nk,Nw,Norb,Norb):: sigma_diis

    ! --- store this iteration's output and residual in the circular buffer ---
    st%i_hist=mod(st%i_hist,st%m)+1
    st%n_hist=min(st%n_hist+1,st%m)
    st%xout(:,:,:,:,st%i_hist)=sigmak(:,:,:,:)
    st%res(:,:,:,:,st%i_hist)=sigmak(:,:,:,:)-sigmak0(:,:,:,:)
    n_cur=st%n_hist

    ! --- Pulay matrix B of size (n_cur+1 x n_cur+1), last row/column the Lagrange constraint ---
    st%B(1:n_cur+1,1:n_cur+1)=0.0d0
    do ih=1,n_cur
       idx_i=modulo(st%i_hist-n_cur+ih-1,st%m)+1
       do jh=1,n_cur
          idx_j=modulo(st%i_hist-n_cur+jh-1,st%m)+1
          st%B(ih,jh)=dble(sum(conjg(st%res(:,:,:,:,idx_i))*st%res(:,:,:,:,idx_j)))
       end do
       st%B(ih,n_cur+1)=1.0d0
       st%B(n_cur+1,ih)=1.0d0
    end do
    eps_reg=1.0d-12                      ! regularize the diagonal against exact singularity
    if(n_cur>0)then
       do ih=1,n_cur
          st%B(ih,ih)=st%B(ih,ih)+eps_reg
       end do
    end if
    st%rhs(1:n_cur)=0.0d0
    st%rhs(n_cur+1)=1.0d0                ! sum c_i = 1
    call dgesv(n_cur+1,1,st%B(1:n_cur+1,1:n_cur+1),n_cur+1, &
               st%ipiv(1:n_cur+1),st%rhs(1:n_cur+1),n_cur+1,info)
    if(info/=0)then
       print*,'DIIS: dgesv failed (info=',info,'), fallback to most-recent entry'
       st%rhs(1:n_cur)=0.0d0
       if(n_cur>=1) st%rhs(n_cur)=1.0d0
    end if

    ! --- extrapolation sigma_diis = sum_i c_i * xout_i ---
    sigma_diis(:,:,:,:)=0.0d0
    do ih=1,n_cur
       idx_i=modulo(st%i_hist-n_cur+ih-1,st%m)+1
       sigma_diis=sigma_diis+st%rhs(ih)*st%xout(:,:,:,:,idx_i)
    end do

    ! --- convergence estimate and mixing ---
    esterr=0.0d0
    est=100.0d0
    kerr=1; iwerr=1; lerr=1; merr=1      ! defaults in case every |sigma| < eps_sgm
    do l=1,Norb
       do m=1,Norb
          do j=1,Nw
             do i=1,Nk
                cksigm=sigma_diis(i,j,m,l)
                if(abs(cksigm)>eps_sgm)then
                   est=abs((sigmak0(i,j,m,l)-cksigm)/cksigm)
                   if(est>esterr)then
                      esterr=est
                      kerr=i
                      iwerr=j
                      lerr=l
                      merr=m
                   end if
                end if
                if(n_cur>=2)then         ! DIIS active: keep the extrapolation as it is
                   sigmak(i,j,m,l)=cksigm
                else                     ! nothing to extrapolate from yet: linear mixing
                   sigmak(i,j,m,l)=pp*cksigm+(1.0d0-pp)*sigmak0(i,j,m,l)
                end if
             end do
          end do
       end do
    end do
    do i=1,Nkall
       if(invk(1,i)==kerr)then
          print '(A7,E9.2,A14,3(1X,I3),I5,2I3,I5)','esterr=',esterr,' at k,iw,m,l=', &
               & kmap(1,i),kmap(2,i),kmap(3,i),iwerr,merr,lerr,kerr
          exit
       end if
    end do
  end subroutine mix_sigma

  subroutine solve_mu(mu,mu_OLD,bracket_valid,muS_cache,muL_cache,rfill,esterr,get_rn)
    !> Move the chemical potential so that get_rn(mu) = rfill.
    !>
    !> Two stages: widen [muS,muL] until it brackets the target (reusing the bracket cached from
    !> the previous call when it is still valid, which is the common case once the loop settles),
    !> then Brent's method (inverse quadratic interpolation with a bisection fallback) on
    !> get_rn(mu) - rfill.  The tolerance follows esterr, the current Sigma convergence error:
    !> there is no point solving mu far more accurately than Sigma itself is converged.
    !!@param          mu,inout: chemical potential, updated in place
    !!@param      mu_OLD,inout: previous mu; sets the initial bracket width
    !!@param bracket_valid,inout: whether muS_cache/muL_cache hold a usable bracket
    !!@param muS_cache,muL_cache,inout: cached bracket
    !!@param           rfill,in: target electron number
    !!@param          esterr,in: current Sigma convergence error, sets the mu tolerance
    !!@param          get_rn,in: electron number at a given mu
    real(c_double),intent(inout):: mu,mu_OLD,muS_cache,muL_cache
    logical,intent(inout):: bracket_valid
    real(c_double),intent(in):: rfill,esterr
    procedure(flex_rn_if):: get_rn

    integer(c_int32_t) i_iter
    integer(c_int32_t),parameter:: itemax=100
    logical flag
    logical bracket_found
    real(c_double) rnS,rnL,rnc,rnM,muc,mud,muL,muS,muM,eps,dmu

      if(esterr>1.0d-2)then
         eps= 1.0d-8
      else
         eps= esterr*1.0d-1
      end if
      dmu= abs(mu-mu_OLD)*2.0d0
      if (dmu<eps*4.0d0) dmu= eps*4.0d0
      mu_OLD=mu

      bracket_found=.false.
      if(bracket_valid)then
         call get_rn(rnS,muS_cache)
         call get_rn(rnL,muL_cache)
         if(rnS<rfill .and. rnL>rfill)then
            muS=muS_cache
            muL=muL_cache
            bracket_found=.true.
            print'(A)','[mu] cached bracket valid'
         else
            print'(A)','[mu] cached bracket invalid, re-searching'
         end if
      end if

      if(.not.bracket_found)then
         muL= mu+dmu
         muS= mu-dmu
         upper_lim: do i_iter=1,itemax
            mu=muL
            rnL=0.0d0
            call get_rn(rnL,mu)
            if(rnL<rfill)then
               if(abs(rfill-rnL)>0.5d0)then
                  muL=muL+1.0d0
               else if(abs(rfill-rnL)>dmu)then
                  muL=muL+0.5d0
               else
                  muL= muL +dmu
               end if
            else
               exit
            end if
            if(i_iter==itemax)then
               print*,'Too many'
               stop
            end if
         end do upper_lim

         lower_lim: do i_iter=1,itemax
            mu=muS
            rnS=0.0d0
            call get_rn(rnS,mu)
            if(rnS>rfill)then
               if(abs(rnS-rfill)>0.5d0)then
                  muS=muS-1.0d0
               else if(abs(rnS-rfill)>dmu)then
                  muS=muS-0.5d0
               else
                  muS=muS-dmu
               end if
            else
               exit
            end if
            if(i_iter==itemax)then
               print*,'Too many'
               stop
            end if
         end do lower_lim
      end if

      ! save wide bracket before Brent narrows it
      muS_cache=muS
      muL_cache=muL
      bracket_valid=.true.

      rnL=rnL-rfill
      rnS=rnS-rfill
      rnc=rnS
      muc=muS
      mud=0.0d0
      flag=.false.
      brent_loop: do i_iter=1,itemax
         if(rnc/=rnS .and. rnc/=rnL)then
            muM=(muL*rnS*rnc*(rnS-rnc)+muS*rnL*rnc*(rnc-rnL)+muc*rnL*rnS*(rnL-rnS))&
                 /((rnL-rnS)*(rnc-rnL)*(rnc-rnS))
         else
            muM=muL-rnL*(muL-muS)/(rnL-rnS)
         end if
         if((0.25d0*(3.0d0*muS+muL) > muM .or. muM > muL) .or. &
              (flag .and. abs(muM-muL) >= abs(muL-muc)*0.5d0) .or. &
              (flag .eqv. .false. .and. abs(muM-muL) >= abs(muc-mud)*0.5d0) .or. &
              (flag .and. abs(muL-muc)<1.0d-8) .or. &
              (flag.eqv. .false. .and. abs(muc-mud)<1.0d-8))then
            muM= (muL+muS)*0.5d0
            flag=.true.
         else
            flag=.false.
         end if
         if(abs(muL-muM)<eps)exit
         mu=muM
         rnM=0
         call get_rn(rnM,muM)
         !print '(1x,a,2f22.16,l2)','muM,rnM=   ',muM,rnM,flag
         mud=muc
         muc=muL
         rnc=rnL
         if(rnS*(rnM-rfill)<0)then
            muL=muM
            rnL=rnM-rfill
         else
            muS=muM
            rnS=rnM-rfill
         end if
         if(abs(rnS)<abs(rnL))then
            muM=muL
            muL=muS
            muS=muM
            rnM=rnL
            rnL=rnS
            rnS=rnM
         end if
         if(i_iter==itemax)then
            print *,'Too many loop!'
            stop
         end if
      end do brent_loop
      if(rnL==rnS)then
         mu=(muS+muL)*0.5d0
      else
         mu= (muS*rnL-muL*rnS)/(rnL-rnS)
      end if
      call get_rn(rnM,mu)
      mu_OLD=mu
      print'(A4,F8.4,A10,F8.4)','mu  =',mu,' filling =',rnM
  end subroutine solve_mu
end module flex_scf

module occupation
  !> Smooth Lindhard occupation factor shared by the susceptibility kernels
  !> (fchi.f90 / fchisc.f90) and the Kubo transport kernel (fcond.f90).
  use,intrinsic:: iso_c_binding, only:c_double
  implicit none
contains
  pure function occ_factor(el,em,fl,fm,temp) result(p)
    !> p = (f_l - f_m)/(e_m - e_l), evaluated so that it stays CONTINUOUS through
    !> the degenerate point e_m = e_l, where it tends to -df/de = f(1-f)/T.
    !>
    !> Writing the Lindhard weight as p/(w + e_m - e_l + i*delta) makes the
    !> intra-band (Pauli / Drude) term and the inter-band term one and the same
    !> expression.  A hard "degenerate if |de| < 1e-9, else divide by
    !> de + i*delta" split instead loses the whole -df/de weight for every pair
    !> whose splitting lies between that threshold and delta - which for a
    !> weakly split band pair is a factor-of-two error in chi0(q,0).
    !>
    !> The finite difference is replaced by the analytic derivative well before
    !> catastrophic cancellation matters: at |de| = 1e-6 T the neglected term is
    !> O((de/T)^2) ~ 1e-12 relative, the cancellation error only ~1e-10.
    !!@param el,em: band energies (em is the one the derivative is taken at)
    !!@param fl,fm: Fermi factors at el, em
    !!@param  temp: temperature, ALREADY regularized by the caller (temp > 0)
    real(c_double),intent(in):: el,em,fl,fm,temp
    real(c_double):: p,de,de_min
    de=em-el
    de_min=max(1.0d-6*temp,1.0d-12)
    if(abs(de)<de_min)then
       p=fm*(1.0d0-fm)/temp
    else
       p=(fl-fm)/de
    end if
  end function occ_factor

  pure function pair_factor(xl,xm,fl,fm,temp) result(p)
    !> Particle-particle (Cooper) analogue of occ_factor:
    !>     p = (1 - f_l - f_m)/(xi_m + xi_l)
    !> with xi measured from mu.  Both numerator and denominator vanish at
    !> xi_l = -xi_m (the pair-degenerate line that carries the BCS log), and
    !> since 1 - f(xi_l) - f(xi_m) = [tanh(b xi_l/2) + tanh(b xi_m/2)]/2 the
    !> ratio tends to -df/dxi = f(1-f)/T there - finite, and in fact the LARGEST
    !> weight in the whole sum.  Dropping those terms (as the old
    !> |denominator| < eps guard in calc_phi did) removes exactly the states
    !> that drive the pairing instability.
    !!@param xl,xm: band energies measured from mu
    !!@param fl,fm: Fermi factors at xl, xm
    !!@param  temp: temperature, ALREADY regularized by the caller (temp > 0)
    real(c_double),intent(in):: xl,xm,fl,fm,temp
    real(c_double):: p,s,s_min
    s=xm+xl
    s_min=max(1.0d-6*temp,1.0d-12)
    if(abs(s)<s_min)then
       p=fm*(1.0d0-fm)/temp
    else
       p=(1.0d0-fl-fm)/s
    end if
  end function pair_factor
end module occupation

subroutine openmp_params(omp_num,omp_check) bind(C)
  !$ use omp_lib
  use,intrinsic:: iso_c_binding, only:c_int64_t
  implicit none
  integer(c_int64_t),intent(out):: omp_num
  logical(1),intent(out):: omp_check
  !$ if(.true.)then
  !$   omp_check=.true.
  !$   omp_num=omp_get_max_threads()
  !$ else
       omp_check=.false.
       omp_num=0
  !$ end if
end subroutine openmp_params

! Caches forward/backward FFTW plans across calls within a loop.
! Call init_fft_plans before the orbital loop and destroy_fft_plans after.
!
! Thread-safety contract:
!  - The cache is a single global pair of plans: init/destroy must bracket one
!    serial region of use and must NOT be nested or used by two loops at once.
!  - FFTW plan creation/destruction is NOT thread-safe; all plan create/destroy
!    calls below are serialized with a named critical section.
!  - Executing a cached plan via dfftw_execute_dft (new-array execute) IS
!    thread-safe, but the arrays must have the same size/alignment as the ones
!    used at planning time.
module fftw_plan_cache
  use,intrinsic::iso_c_binding, only:c_int64_t
  implicit none
  integer(c_int64_t),save :: plan_fwd=0, plan_bwd=0
end module fftw_plan_cache

subroutine init_fft_plans(cmat,tmp,Nx,Ny,Nz,Nw)
  use fftw_plan_cache
  use,intrinsic::iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nx,Ny,Nz,Nw
  complex(c_double),intent(inout),dimension(Nx,Ny,Nz,Nw):: cmat,tmp
  integer(c_int32_t),dimension(4):: Nlist
  Nlist=(/int(Nx,c_int32_t),int(Ny,c_int32_t),int(Nz,c_int32_t),int(Nw,c_int32_t)/)
  !$omp critical(fftw_plan_ops)
  if(plan_fwd/=0 .or. plan_bwd/=0)then
     print*,'WARNING(init_fft_plans): plan cache already in use; overwriting (plans leak).'
  end if
  call dfftw_plan_dft(plan_fwd,4,Nlist,cmat,tmp,-1,64)
  call dfftw_plan_dft(plan_bwd,4,Nlist,cmat,tmp, 1,64)
  !$omp end critical(fftw_plan_ops)
end subroutine init_fft_plans

subroutine destroy_fft_plans()
  use fftw_plan_cache
  implicit none
  !$omp critical(fftw_plan_ops)
  call dfftw_destroy_plan(plan_fwd)
  call dfftw_destroy_plan(plan_bwd)
  plan_fwd=0
  plan_bwd=0
  !$omp end critical(fftw_plan_ops)
end subroutine destroy_fft_plans

subroutine FFT(cmat,tmp,Nx,Ny,Nz,Nw,SW)
  !> 4D Fast Fourier transform function
  !!@param cmat,inout: FFT source and results
  !!@param  tmp,inout: temporary
  !!@param      Nx,in: x mesh
  !!@param      Ny,in: y mesh
  !!@param      Nz,in: z mesh
  !!@param      Nw,in: w mesh
  !!@param      SW,in: .true. forward (sign -1, unnormalized);
  !!                   .false. backward (sign +1, normalized by 1/(Nx*Ny*Nz*Nw))
  use fftw_plan_cache
  use,intrinsic::iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nx,Ny,Nz,Nw
  logical,intent(in):: SW
  complex(c_double),intent(inout),dimension(Nx,Ny,Nz,Nw):: cmat,tmp

  integer(c_int64_t) plan
  integer(c_int32_t) Inv
  integer(c_int32_t),dimension(4):: Nlist

  Nlist=(/int(Nx,c_int32_t),int(Ny,c_int32_t),int(Nz,c_int32_t),int(Nw,c_int32_t)/)
  if(plan_fwd /= 0)then
     ! Reuse cached plans: avoids heap alloc/free on every call
     if(SW)then
        call dfftw_execute_dft(plan_fwd,cmat,tmp)
     else
        call dfftw_execute_dft(plan_bwd,cmat,tmp)
     end if
  else
     ! No cached plan: plan/execute/destroy per call.
     ! Plan creation/destruction is serialized (FFTW requires it); execution is thread-safe.
     if(SW)then
        Inv=-1
     else
        Inv=1
     end if
     !$omp critical(fftw_plan_ops)
     call dfftw_plan_dft(plan,4,Nlist,cmat,tmp,Inv,64)
     !$omp end critical(fftw_plan_ops)
     call dfftw_execute(plan)
     !$omp critical(fftw_plan_ops)
     call dfftw_destroy_plan(plan)
     !$omp end critical(fftw_plan_ops)
  end if
  if(.not. SW)then
     cmat(:,:,:,:)=tmp(:,:,:,:)/product(Nlist)
  else
     cmat(:,:,:,:)=tmp(:,:,:,:)
  end if
end subroutine FFT

subroutine fft3d_batch(cmat,tmp,Nx,Ny,Nz,M,SW)
  !> Batched 3D spatial FFT: transforms the first three (k-space) axes of
  !> cmat(Nx,Ny,Nz,M) independently for each of the M trailing slices, leaving
  !> the 4th axis untouched.  Used by the chi0 tail-corrected reference bubble,
  !> where the 4th axis is imaginary time and must NOT be transformed.
  !!@param cmat,inout: FFT source; overwritten with the result
  !!@param  tmp,inout: work array of the same shape
  !!@param Nx,Ny,Nz,in: spatial mesh
  !!@param       M,in: number of trailing slices (batch size)
  !!@param      SW,in: .true. forward (sign -1, unnormalized);
  !!                   .false. backward (sign +1, normalized by 1/(Nx*Ny*Nz))
  use,intrinsic::iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nx,Ny,Nz,M
  logical,intent(in):: SW
  complex(c_double),intent(inout),dimension(Nx,Ny,Nz,M):: cmat,tmp

  integer(c_int64_t) plan
  integer(c_int32_t) Inv,dist,howmany
  integer(c_int32_t),dimension(3):: Nlist

  Nlist=(/int(Nx,c_int32_t),int(Ny,c_int32_t),int(Nz,c_int32_t)/)
  dist=int(Nx*Ny*Nz,c_int32_t)
  howmany=int(M,c_int32_t)
  if(SW)then
     Inv=-1
  else
     Inv=1
  end if
  !$omp critical(fftw_plan_ops)
  call dfftw_plan_many_dft(plan,3,Nlist,howmany,cmat,Nlist,1,dist,&
       tmp,Nlist,1,dist,Inv,64)
  !$omp end critical(fftw_plan_ops)
  call dfftw_execute(plan)
  !$omp critical(fftw_plan_ops)
  call dfftw_destroy_plan(plan)
  !$omp end critical(fftw_plan_ops)
  if(SW)then
     cmat(:,:,:,:)=tmp(:,:,:,:)
  else
     cmat(:,:,:,:)=tmp(:,:,:,:)/dble(dist)
  end if
end subroutine fft3d_batch

subroutine gen_ham(ham_k,klist,ham_r,rvec,Nk,Nr,Norb) bind(C)
  !> This function generate model hamiltonian from hoppings
  !!@param ham_k,out: k-space Hamiltonian
  !!@param  klist,in: k-points list
  !!@param  ham_r,in: hopping integrals
  !!@param   rvec,in: The list of lattice site
  !!@param     Nk,in: The number of k-points
  !!@param     Nr,in: The number of r-vector
  !!@param   Norb,in: The number of orbitals
  use constant
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nr,Norb
  real(c_double),intent(in),dimension(3,Nk):: klist
  real(c_double),intent(in),dimension(3,Nr):: rvec
  complex(c_double),intent(in),dimension(Norb,Norb,Nr):: ham_r
  complex(c_double),intent(out),dimension(Norb,Norb,Nk):: ham_k

  integer(c_int32_t) i,j,l,m
  real(c_double) phase
  complex(c_double) cphase

  ham_k(:,:,:)=(0.0d0,0.0d0)

  ! H(k) = sum_R  H(R) * exp(-2pi*i * k.R)
  !$omp parallel do private(l,m,j,phase,cphase)
  klop: do i=1,Nk
     rloop: do j=1,Nr
        phase=2*pi*sum(klist(:,i)*rvec(:,j))
        cphase=cmplx(cos(phase),-sin(phase),kind=c_double)   ! exp(-i*phi)
        do l=1,Norb
           do m=l,Norb
              ham_k(m,l,i)=ham_k(m,l,i)+ham_r(m,l,j)*cphase
           end do
        end do
     end do rloop
     ! Enforce Hermitian symmetry: diagonal must be real; upper = conj(lower)
     do l=1,Norb
        ham_k(l,l,i)=dble(ham_k(l,l,i))
        do m=l+1,Norb
           ham_k(l,m,i)=conjg(ham_k(m,l,i))
        end do
     end do
  end do klop
  !$omp end parallel do
end subroutine gen_ham

subroutine get_eig(eig,uni,ham_k,Nk,Norb) bind(C)
  !> This function obtains energies
  !!@param  eig,out: energies at k-points
  !!@param  uni,out: unitary matrix at k-points
  !!@param ham_k,in: k-space Hamiltonian
  !!@param    Nk,in: The number of k-points
  !!@param  Norb,in: The number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,Norb
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: ham_k
  real(c_double),intent(out),dimension(Norb,Nk):: eig
  complex(c_double),intent(out),dimension(Norb,Norb,Nk):: uni

  integer(c_int32_t) i,info
  integer(c_int64_t) lwork
  real(c_double) rwork(3*Norb-2),eq(Norb)
  complex(c_double) work_query(1),en(Norb,Norb)
  complex(c_double),allocatable :: work(:)

  call zheev('V','U',Norb,en,Norb,eq,work_query,-1_c_int64_t,rwork,info)
  lwork = int(dble(work_query(1)), c_int64_t)

  !$omp parallel private(en,eq,work,rwork,info)
  allocate(work(lwork))
  !$omp do
  kloop: do i=1,Nk
     en(:,:)=ham_k(:,:,i)
     call zheev('V','U',Norb,en,Norb,eq,work,lwork,rwork,info)
     if(info/=0)then
        print*,'zheev failed in get_eig: info=',info
        stop
     end if
     uni(:,:,i)=en(:,:)
     eig(:,i)=eq(:)
  end do kloop
  !$omp end do
  deallocate(work)
  !$omp end parallel
end subroutine get_eig

subroutine get_eig_mlo(eig,uni,ham_k,Ovlk,Nk,Norb) bind(C)
  !> This function obtains energies
  !!@param  eig,out: energies at k-points
  !!@param  uni,out: unitary matrix at k-points
  !!@param ham_k,in: k-space Hamiltonian
  !!@param  Ovlk,in: k-space Overlap integrals
  !!@param    Nk,in: The number of k-points
  !!@param  Norb,in: The number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,norb
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: ham_k,Ovlk
  real(c_double),intent(out),dimension(Norb,Nk):: eig
  complex(c_double),intent(out),dimension(Norb,Norb,Nk):: uni

  integer(c_int32_t) i,j,k,l,m,info
  integer(c_int64_t) lwork
  real(c_double) rwork(3*Norb-2),eq(Norb),norm
  complex(c_double) work_query(1)
  complex(c_double),allocatable :: work(:)
  complex(c_double),dimension(Norb,Norb):: tmp,tmp2,tmp3
  real(c_double),parameter:: ovl_thresh=1.0d-8

  call zheev('V','U',norb,tmp,norb,eq,work_query,-1_c_int64_t,rwork,info)
  lwork = int(dble(work_query(1)), c_int64_t)

  !$omp parallel private(tmp,tmp2,tmp3,norm,eq,work,rwork,info)
  allocate(work(lwork))
  !$omp do
  kloop: do i=1,Nk
     tmp(:,:)=Ovlk(:,:,i)
     call zheev('V','U',norb,tmp,norb,eq,work,lwork,rwork,info)
     if(info/=0)then
        print*,'zheev failed in get_eig_mlo (overlap): info=',info
        stop
     end if
     ! Canonical orthogonalization: rescale eigenvectors by 1/sqrt(eigenvalue) of S
     do j=1,Norb
        if(eq(j)>ovl_thresh)then
           tmp2(:,j)=tmp(:,j)/sqrt(cmplx(eq(j),kind=c_double))
        else
           tmp2(:,j)=(0.0d0,0.0d0) !discard near-null basis vector (near-linear dependence)
        end if
     end do
     ! Transform H to orthonormal basis: H_orth = tmp2^H * H * tmp2
     call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),ham_k(:,:,i),Norb,tmp2,Norb,(0.0d0,0.0d0),tmp3,Norb)
     call zgemm('C','N',Norb,Norb,Norb,(1.0d0,0.0d0),tmp2,Norb,tmp3,Norb,(0.0d0,0.0d0),tmp,Norb)
     call zheev('V','U',norb,tmp,norb,eq,work,lwork,rwork,info)
     if(info/=0)then
        print*,'zheev failed in get_eig_mlo (ham): info=',info
        stop
     end if
     eig(:,i)=eq(:)
     ! uni = tmp2 * tmp  (compose transformations)
     call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),tmp2,Norb,tmp,Norb,(0.0d0,0.0d0),tmp3,Norb)
     uni(:,:,i)=tmp3(:,:)
  end do kloop
  !$omp end do
  deallocate(work)
  !$omp end parallel
end subroutine get_eig_mlo

subroutine get_ffermi(ffermi,eig,mu,temp,Nk,Norb) bind(C)
  !> This function obtains fermi functions
  !!@param ffermi,out: fermi distribute functions
  !!@param     eig,in: energies at k-points
  !!@param      mu,in: chemical potential
  !!@param    temp,in: Temperature
  !!@param      Nk,in: The number of k-points
  !!@param    Norb,in: The number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,Norb
  real(c_double),intent(in):: mu,temp
  real(c_double),intent(in),dimension(Norb,Nk):: eig
  real(c_double),intent(out),dimension(Norb,Nk):: ffermi

  integer(c_int32_t) i,j
  real(c_double) itemp,temp_safe

  temp_safe=max(temp,1.0d-12)
  itemp=0.5d0/temp_safe
  ! f(e) = 1/(exp((e-mu)/T)+1) = 0.5 - 0.5*tanh((e-mu)/(2T))  [numerically stable form]
  !$omp parallel do private(j)
  do i=1,Nk
     do j=1,Norb
        ffermi(j,i)=0.5d0-0.5d0*tanh((eig(j,i)-mu)*itemp)
     end do
  end do
  !$omp end parallel do
end subroutine get_ffermi

subroutine get_imass0(imk,klist,ham_r,rvec,Nk,Nr,Norb) bind(C)
  use constant
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nr,norb
  real(c_double),intent(in),dimension(3,Nk):: klist
  real(c_double),intent(in),dimension(3,Nr):: rvec
  complex(c_double),intent(in),dimension(Norb,Norb,Nr):: ham_r
  complex(c_double),intent(out),dimension(3,3,Norb,Norb,Nk):: imk

  integer(c_int32_t) i,j,k,l,m,n
  real(c_double) phase
  complex(c_double) cphase

  imk(:,:,:,:,:)=0.0d0

  !$omp parallel do private(l,m,j,k,n,phase,cphase)
  kloop: do i=1,Nk
     rloop: do j=1,Nr
        phase=2*pi*sum(klist(:,i)*rvec(:,j))
        cphase=cmplx(cos(phase),-sin(phase),kind=c_double)
        do l=1,Norb
           do m=l,Norb
              !$omp simd
              axis1: do k=1,3
                 axis2: do n=k,3
                    imk(n,k,m,l,i)=imk(n,k,m,l,i)-rvec(n,j)*rvec(k,j)*ham_r(m,l,j)*cphase
                 end do axis2
              end do axis1
              !$omp end simd
           end do
        end do
     end do rloop
     ! The r-loop only accumulates the axis lower triangle (n>=k) of the orbital
     ! lower triangle (m>=l).  Fill the rest from the two independent symmetries:
     !   orbital: M_{lm} = conjg(M_{ml})   (H(R) Hermitian after the R-sum)
     !   axis   : M_{kn} = M_{nk}          (d2/dk_n dk_k is symmetric, NOT conjugated)
     ! Both the (m,l) and the (l,m) orbital blocks need the axis copy; conjugating
     ! the axis transpose instead would corrupt every m/=l element and leave
     ! imk(k,n,l,m) (k<n) at zero.
     do l=1,Norb
        do m=l,Norb
           axis12: do k=1,3
              axis22: do n=k,3
                 imk(n,k,l,m,i)=conjg(imk(n,k,m,l,i))
                 if(n/=k)then
                    imk(k,n,m,l,i)=imk(n,k,m,l,i)
                    imk(k,n,l,m,i)=imk(n,k,l,m,i)
                 end if
              end do axis22
           end do axis12
        end do
     end do
  end do kloop
  !$omp end parallel do
end subroutine get_imass0

subroutine get_imassk(imk,imk0,vk0,eig,mrot,uni,Nk,Norb) bind(C)
  !> This function obtains the band-basis inverse effective mass tensor
  !!   d2 eps_n/dk_a dk_b = (U^H d2H/dk_a dk_b U)_nn
  !!                      + 2 sum_{m/=n} Re[H~^a_nm H~^b_mn]/(eps_n - eps_m),  H~^a = U^H dH/dk_a U
  !! The second (interband / level-repulsion) term is 2nd-order perturbation theory and is
  !! the SAME order as the first one in a multiband model -- it used to be missing here.
  !! It diverges at a band crossing, where a band-resolved mass simply does not exist, so
  !! pairs closer than deg_thresh are dropped rather than blown up.
  !!@param  imk,out: band-basis inverse effective mass tensor
  !!@param imk0,in: orbital-basis d2H/dk_a dk_b (get_imass0)
  !!@param  vk0,in: orbital-basis dH/dk_a (get_vlm0)
  !!@param  eig,in: band energies at k-points
  !!@param mrot,in: rotation matrix for real space
  !!@param  uni,in: unitary matrix orbital to band
  !!@param   Nk,in: The number of k-points
  !!@param Norb,in: The number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,Norb
  real(c_double),intent(in),dimension(3,3):: mrot
  real(c_double),intent(in),dimension(Norb,Nk):: eig
  complex(c_double),intent(in),dimension(3,3,Norb,Norb,Nk):: imk0
  complex(c_double),intent(in),dimension(3,Norb,Norb,Nk):: vk0
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni
  real(c_double),intent(out),dimension(3,3,Norb,Nk):: imk

  integer(c_int32_t) i,j_ax,k_ax,l,m,n
  real(c_double) de,acc
  complex(c_double) tmp(3,3,Norb)
  complex(c_double) Mtmp(Norb,Norb),Wtmp(Norb,Norb),res_mat(Norb,Norb)
  complex(c_double) vb(Norb,Norb,3)
  real(c_double),parameter:: deg_thresh=1.0d-6

  !$omp parallel do private(tmp,Mtmp,Wtmp,res_mat,vb,de,acc,j_ax,k_ax,l,m,n)
  kloop: do i=1,Nk
     !band-basis velocity matrices vb(:,:,a) = U^H (dH/dk_a) U  (needed by the interband term)
     do k_ax=1,3
        Mtmp(:,:)=vk0(k_ax,:,:,i)
        call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),Mtmp,Norb,uni(:,:,i),Norb,(0.0d0,0.0d0),Wtmp,Norb)
        call zgemm('C','N',Norb,Norb,Norb,(1.0d0,0.0d0),uni(:,:,i),Norb,Wtmp,Norb,(0.0d0,0.0d0),res_mat,Norb)
        vb(:,:,k_ax)=res_mat(:,:)
     end do
     !rotate orb to band via zgemm: tmp(k_ax,j_ax,n)=diag(U^H * imk0(k_ax,j_ax,:,:,i) * U)_n
     do k_ax=1,3
        do j_ax=1,3
           Mtmp(:,:)=imk0(k_ax,j_ax,:,:,i)
           call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),Mtmp,Norb,uni(:,:,i),Norb,(0.0d0,0.0d0),Wtmp,Norb)
           call zgemm('C','N',Norb,Norb,Norb,(1.0d0,0.0d0),uni(:,:,i),Norb,Wtmp,Norb,(0.0d0,0.0d0),res_mat,Norb)
           do n=1,norb
              ! 2 sum_{m/=n} Re[H~^a_nm H~^b_mn]/(eps_n-eps_m).  Re of the PRODUCT: taking the
              ! real part of each factor first would drop the interference between the axes.
              acc=0.0d0
              do m=1,norb
                 de=eig(n,i)-eig(m,i)
                 if(abs(de)<deg_thresh) cycle       ! also skips m==n
                 acc=acc+2.0d0*dble(vb(n,m,k_ax)*vb(m,n,j_ax))/de
              end do
              tmp(k_ax,j_ax,n)=res_mat(n,n)+cmplx(acc,0.0d0,kind=c_double)
           end do
        end do
     end do
     !rotate axis
     imk(:,:,:,i)=0.0d0
     band_loop2: do n=1,norb
        do l=1,3
           do m=1,3
              do j_ax=1,3
                 do k_ax=1,3
                    imk(m,l,n,i)=imk(m,l,n,i)+mrot(m,k_ax)*mrot(l,j_ax)*dble(tmp(k_ax,j_ax,n))
                 end do
              end do
           end do
        end do
     end do band_loop2
  end do kloop
  !$omp end parallel do
end subroutine get_imassk

subroutine get_imassk_mlo(imk,imk0,smk0,vk0,sk0,eig,mrot,uni,Nk,Norb) bind(C)
  !> This function obtains the band-basis inverse effective mass tensor for a
  !! NON-ORTHOGONAL (MLO) basis, H(k)C = eps S(k)C with C^H S C = 1.  Writing
  !!   H~^a = C^H dH/dk_a C,  S~^a = C^H dS/dk_a C,  G^a_ln = H~^a_ln - eps_n S~^a_ln
  !! (so G^a_nn = d eps_n/dk_a is the get_veloc_mlo velocity),
  !!   d2 eps_n/dk_a dk_b = [C^H (d2H/dk_a dk_b - eps_n d2S/dk_a dk_b) C]_nn
  !!                      - (d eps_n/dk_a) S~^b_nn - (d eps_n/dk_b) S~^a_nn
  !!                      + 2 sum_{l/=n} Re[conjg(G^b_ln) G^a_ln]/(eps_n - eps_l)
  !! The two extra terms beyond the orthogonal case come from the k-dependence of the
  !! normalization C^H S C = 1 (which fixes Re D^b_nn = -S~^b_nn/2 for dC = C D) and from
  !! the eigenvector response D^b_ln = G^b_ln/(eps_n - eps_l).  Note eps_n (the band whose
  !! mass is wanted), NOT eps_l, sits inside G -- that asymmetry is real.
  !! Setting S=1 reduces every line to get_imassk.
  !!@param  imk,out: band-basis inverse effective mass tensor
  !!@param imk0,in: orbital-basis d2H/dk_a dk_b (get_imass0 on ham_r)
  !!@param smk0,in: orbital-basis d2S/dk_a dk_b (get_imass0 on S_r)
  !!@param  vk0,in: orbital-basis dH/dk_a (get_vlm0 on ham_r)
  !!@param  sk0,in: orbital-basis dS/dk_a (get_vlm0 on S_r)
  !!@param  eig,in: generalized eigenvalues (get_eig_mlo)
  !!@param mrot,in: rotation matrix for real space
  !!@param  uni,in: generalized eigenvectors, S-orthonormal (get_eig_mlo)
  !!@param   Nk,in: The number of k-points
  !!@param Norb,in: The number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,Norb
  real(c_double),intent(in),dimension(3,3):: mrot
  real(c_double),intent(in),dimension(Norb,Nk):: eig
  complex(c_double),intent(in),dimension(3,3,Norb,Norb,Nk):: imk0,smk0
  complex(c_double),intent(in),dimension(3,Norb,Norb,Nk):: vk0,sk0
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni
  real(c_double),intent(out),dimension(3,3,Norb,Nk):: imk

  integer(c_int32_t) i,j_ax,k_ax,l,m,n
  real(c_double) de,acc,en
  complex(c_double) gk,gj
  complex(c_double) tmp(3,3,Norb)
  complex(c_double) Mtmp(Norb,Norb),Wtmp(Norb,Norb),res_h(Norb,Norb),res_s(Norb,Norb)
  complex(c_double) vb(Norb,Norb,3),sb(Norb,Norb,3)
  real(c_double) vd(3,Norb)
  real(c_double),parameter:: deg_thresh=1.0d-6

  !$omp parallel do private(tmp,Mtmp,Wtmp,res_h,res_s,vb,sb,vd,gk,gj,de,acc,en,j_ax,k_ax,l,m,n)
  kloop: do i=1,Nk
     !band-basis first derivatives: vb(:,:,a)=C^H (dH/dk_a) C, sb(:,:,a)=C^H (dS/dk_a) C
     do k_ax=1,3
        Mtmp(:,:)=vk0(k_ax,:,:,i)
        call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),Mtmp,Norb,uni(:,:,i),Norb,(0.0d0,0.0d0),Wtmp,Norb)
        call zgemm('C','N',Norb,Norb,Norb,(1.0d0,0.0d0),uni(:,:,i),Norb,Wtmp,Norb,(0.0d0,0.0d0),res_h,Norb)
        vb(:,:,k_ax)=res_h(:,:)
        Mtmp(:,:)=sk0(k_ax,:,:,i)
        call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),Mtmp,Norb,uni(:,:,i),Norb,(0.0d0,0.0d0),Wtmp,Norb)
        call zgemm('C','N',Norb,Norb,Norb,(1.0d0,0.0d0),uni(:,:,i),Norb,Wtmp,Norb,(0.0d0,0.0d0),res_s,Norb)
        sb(:,:,k_ax)=res_s(:,:)
        ! band velocity d eps_n/dk_a = G^a_nn (S~^a_nn is real, H~ and S~ are Hermitian)
        do n=1,norb
           vd(k_ax,n)=dble(vb(n,n,k_ax))-eig(n,i)*dble(sb(n,n,k_ax))
        end do
     end do
     do k_ax=1,3
        do j_ax=1,3
           Mtmp(:,:)=imk0(k_ax,j_ax,:,:,i)
           call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),Mtmp,Norb,uni(:,:,i),Norb,(0.0d0,0.0d0),Wtmp,Norb)
           call zgemm('C','N',Norb,Norb,Norb,(1.0d0,0.0d0),uni(:,:,i),Norb,Wtmp,Norb,(0.0d0,0.0d0),res_h,Norb)
           Mtmp(:,:)=smk0(k_ax,j_ax,:,:,i)
           call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),Mtmp,Norb,uni(:,:,i),Norb,(0.0d0,0.0d0),Wtmp,Norb)
           call zgemm('C','N',Norb,Norb,Norb,(1.0d0,0.0d0),uni(:,:,i),Norb,Wtmp,Norb,(0.0d0,0.0d0),res_s,Norb)
           do n=1,norb
              en=eig(n,i)
              acc=0.0d0
              do l=1,norb
                 de=en-eig(l,i)
                 if(abs(de)<deg_thresh) cycle       ! also skips l==n
                 gk=vb(l,n,k_ax)-en*sb(l,n,k_ax)
                 gj=vb(l,n,j_ax)-en*sb(l,n,j_ax)
                 acc=acc+2.0d0*dble(conjg(gj)*gk)/de
              end do
              acc=acc-vd(k_ax,n)*dble(sb(n,n,j_ax))-vd(j_ax,n)*dble(sb(n,n,k_ax))
              tmp(k_ax,j_ax,n)=res_h(n,n)-en*res_s(n,n)+cmplx(acc,0.0d0,kind=c_double)
           end do
        end do
     end do
     !rotate axis
     imk(:,:,:,i)=0.0d0
     band_loop6: do n=1,norb
        do l=1,3
           do m=1,3
              do j_ax=1,3
                 do k_ax=1,3
                    imk(m,l,n,i)=imk(m,l,n,i)+mrot(m,k_ax)*mrot(l,j_ax)*dble(tmp(k_ax,j_ax,n))
                 end do
              end do
           end do
        end do
     end do band_loop6
  end do kloop
  !$omp end parallel do
end subroutine get_imassk_mlo

subroutine get_vlm0(vk,klist,ham_r,rvec,Nk,Nr,Norb) bind(C)
  !> This function obtain orbital basis velocity
  !!@param   vk,out: orbital-basis group velocities
  !!@param klist,in: k-points list
  !!@param ham_r,in: hopping integrals
  !!@param  rvec,in: The list of lattice site
  !!@param    Nk,in: The number of k-points
  !!@param    Nr,in: The number of r-vector
  !!@param  Norb,in: The number of orbitals
  use constant
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nr,Norb
  real(c_double),intent(in),dimension(3,Nk):: klist
  real(c_double),intent(in),dimension(3,Nr):: rvec
  complex(c_double),intent(in),dimension(Norb,Norb,Nr):: ham_r
  complex(c_double),intent(out),dimension(3,Norb,Norb,Nk):: vk

  integer(c_int32_t) i,j,k,l,m
  real(c_double) phase
  complex(c_double) cphase

  vk(:,:,:,:)=0.0d0

  !$omp parallel do private(l,m,j,k,phase,cphase)
  kloop: do i=1,Nk
     rloop: do j=1,Nr
        phase=2*pi*sum(klist(:,i)*rvec(:,j))
        ! d/dk_alpha exp(-i*phi) = -i*R_alpha * exp(-i*phi)
        ! Re/Im part: -(sin+i*cos) = -i*(cos-i*sin) = -i*exp(-i*phi)
        cphase=cmplx(sin(phase),cos(phase),kind=c_double)   ! = -i*exp(-i*phi) in real/imag components
        do l=1,Norb
           do m=l,Norb
              !$omp simd
              vaxis: do k=1,3
                 vk(k,m,l,i)=vk(k,m,l,i)-rvec(k,j)*ham_r(m,l,j)*cphase
              end do vaxis
              !$omp end simd
           end do
        end do
     end do rloop
     do l=1,Norb
        do m=l,Norb
           vaxis2: do k=1,3
              vk(k,l,m,i)=conjg(vk(k,m,l,i))
           end do vaxis2
        end do
     end do
  end do kloop
  !$omp end parallel do
end subroutine get_vlm0

subroutine get_veloc(vk,vk0,mrot,uni,Nk,Norb) bind(C)
  !> This function obtain group velocity at bands
  !!@param  vk,out: band-basis group velocities
  !!@param  vk0,in: orbital-basis group velocities
  !!@param mrot,in: rotation matrix for real space
  !!@param  uni,in: unitary matrix orbital to band
  !!@param   Nk,in: The number of k-points
  !!@param Norb,in: The number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,Norb
  real(c_double),intent(in),dimension(3,3):: mrot
  complex(c_double),intent(in),dimension(3,Norb,Norb,Nk):: vk0
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni
  real(c_double),intent(out),dimension(3,Norb,Nk):: vk

  integer(c_int32_t) i,j_dir,l,m,n
  complex(c_double) tmp(3,Norb)
  complex(c_double) vtmp(Norb,Norb),Wtmp(Norb,Norb),tmp_j(Norb,Norb)

  !$omp parallel do private(tmp,vtmp,Wtmp,tmp_j,j_dir,l,m,n)
  kloop: do i=1,Nk
     !rotate orb to band via zgemm: tmp_j = U^H * vk0(j,:,:,i) * U, take diagonal
     do j_dir=1,3
        vtmp(:,:)=vk0(j_dir,:,:,i)
        call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),vtmp,Norb,uni(:,:,i),Norb,(0.0d0,0.0d0),Wtmp,Norb)
        call zgemm('C','N',Norb,Norb,Norb,(1.0d0,0.0d0),uni(:,:,i),Norb,Wtmp,Norb,(0.0d0,0.0d0),tmp_j,Norb)
        do n=1,Norb
           tmp(j_dir,n)=tmp_j(n,n)
        end do
     end do
     !rotate axis
     vk(:,:,i)=0.0d0
     band_loop2: do n=1,Norb
        do l=1,3
           do m=1,3
              vk(m,n,i)=vk(m,n,i)+mrot(m,l)*dble(tmp(l,n))
           end do
        end do
     end do band_loop2
  end do kloop
  !$omp end parallel do
end subroutine get_veloc

subroutine get_vnm(vk,vk0,mrot,uni,Nk,Norb) bind(C)
  !> This function obtain band-basis group velocity with inter-band
  !!@param  vk,out: band-basis group velocities
  !!@param  vk0,in: orbital-basis group velocities
  !!@param mrot,in: rotation matrix for real space
  !!@param  uni,in: unitary matrix orbital to band
  !!@param   Nk,in: The number of k-points
  !!@param Norb,in: The number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,Norb
  real(c_double),intent(in),dimension(3,3):: mrot
  complex(c_double),intent(in),dimension(3,Norb,Norb,Nk):: vk0
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(c_double),intent(out),dimension(3,Norb,Norb,Nk):: vk

  integer(c_int32_t) i,j_dir,l,m,n,k
  complex(c_double) tmp(3,Norb,Norb)
  complex(c_double),dimension(Norb,Norb):: vtmp,Wtmp,tmp_j

  !$omp parallel
  !$omp workshare
  vk(:,:,:,:)=0.0d0
  !$omp end workshare
  !$omp do private(tmp,vtmp,Wtmp,tmp_j,j_dir,l,m,n,k)
  kloop: do i=1,Nk
     ! tmp(j,:,:) = uni^H * vk0(j,:,:,i) * uni  for j=1..3
     do j_dir=1,3
        vtmp(:,:)=vk0(j_dir,:,:,i)
        call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),vtmp,Norb,uni(:,:,i),Norb,(0.0d0,0.0d0),Wtmp,Norb)
        call zgemm('C','N',Norb,Norb,Norb,(1.0d0,0.0d0),uni(:,:,i),Norb,Wtmp,Norb,(0.0d0,0.0d0),tmp_j,Norb)
        tmp(j_dir,:,:)=tmp_j(:,:)
     end do
     !rotate axis
     band_loop3: do n=1,Norb
        do k=1,Norb
           do l=1,3
              do m=1,3
                 vk(m,k,n,i)=vk(m,k,n,i)+mrot(m,l)*dble(tmp(l,k,n))
              end do
           end do
        end do
     end do band_loop3
  end do kloop
  !$omp end do
  !$omp end parallel
end subroutine get_vnm

subroutine get_veloc_mlo(vk,vk0,sk0,eig,mrot,uni,Nk,Norb) bind(C)
  !> This function obtain group velocity at bands for a NON-ORTHOGONAL (MLO) basis
  !! The generalized problem H(k)C = eps S(k)C with C^H S C = 1 has the extra
  !! overlap-derivative term in the Hellmann-Feynman velocity:
  !!   v_n = (C^H dH/dk C)_nn - eps_n (C^H dS/dk C)_nn
  !! Dropping the second term (i.e. calling get_veloc) makes v_n wrong wherever
  !! dS/dk /= 0, which is everywhere for MLO/ecalj (ftype 3,4) hoppings.
  !!@param  vk,out: band-basis group velocities
  !!@param  vk0,in: orbital-basis dH/dk (get_vlm0 on ham_r)
  !!@param  sk0,in: orbital-basis dS/dk (get_vlm0 on S_r)
  !!@param  eig,in: generalized eigenvalues (get_eig_mlo)
  !!@param mrot,in: rotation matrix for real space
  !!@param  uni,in: generalized eigenvectors, S-orthonormal (get_eig_mlo)
  !!@param   Nk,in: The number of k-points
  !!@param Norb,in: The number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,Norb
  real(c_double),intent(in),dimension(3,3):: mrot
  real(c_double),intent(in),dimension(Norb,Nk):: eig
  complex(c_double),intent(in),dimension(3,Norb,Norb,Nk):: vk0,sk0
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni
  real(c_double),intent(out),dimension(3,Norb,Nk):: vk

  integer(c_int32_t) i,j_dir,l,m,n
  complex(c_double) tmp(3,Norb)
  complex(c_double),dimension(Norb,Norb):: vtmp,Wtmp,tmp_j,tmp_s

  !$omp parallel do private(tmp,vtmp,Wtmp,tmp_j,tmp_s,j_dir,l,m,n)
  kloop: do i=1,Nk
     !rotate orb to band via zgemm: tmp_j = U^H * vk0(j,:,:,i) * U, likewise tmp_s for sk0
     do j_dir=1,3
        vtmp(:,:)=vk0(j_dir,:,:,i)
        call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),vtmp,Norb,uni(:,:,i),Norb,(0.0d0,0.0d0),Wtmp,Norb)
        call zgemm('C','N',Norb,Norb,Norb,(1.0d0,0.0d0),uni(:,:,i),Norb,Wtmp,Norb,(0.0d0,0.0d0),tmp_j,Norb)
        vtmp(:,:)=sk0(j_dir,:,:,i)
        call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),vtmp,Norb,uni(:,:,i),Norb,(0.0d0,0.0d0),Wtmp,Norb)
        call zgemm('C','N',Norb,Norb,Norb,(1.0d0,0.0d0),uni(:,:,i),Norb,Wtmp,Norb,(0.0d0,0.0d0),tmp_s,Norb)
        do n=1,Norb
           tmp(j_dir,n)=tmp_j(n,n)-eig(n,i)*tmp_s(n,n)
        end do
     end do
     !rotate axis
     vk(:,:,i)=0.0d0
     band_loop4: do n=1,Norb
        do l=1,3
           do m=1,3
              vk(m,n,i)=vk(m,n,i)+mrot(m,l)*dble(tmp(l,n))
           end do
        end do
     end do band_loop4
  end do kloop
  !$omp end parallel do
end subroutine get_veloc_mlo

subroutine get_vnm_mlo(vk,vk0,sk0,Ovlk,eig,mrot,uni,Nk,Norb) bind(C)
  !> This function obtain band-basis group velocity with inter-band elements for a
  !! NON-ORTHOGONAL (MLO) basis.  The model is read as an orthogonal tight-binding
  !! model in the Loewdin basis, Htilde = S^-1/2 H S^-1/2, whose velocity operator is
  !! dHtilde/dk.  With X = S^1/2 and E = diag(eig) that is exactly
  !!   M = C^H (dH/dk) C - G E - E G^H ,      G = C^H (dX/dk) X C
  !! and dX/dk follows from X dX + dX X = dS/dk, i.e. in the eigenbasis of S
  !!   (dX)_nm = (U^H (dS/dk) U)_nm / (sqrt(s_n) + sqrt(s_m)) .
  !! On the diagonal dX X + X dX = dS collapses the two G terms to eps_n (dS/dk)_nn,
  !! reproducing the get_veloc_mlo group velocity, so the two routines stay consistent.
  !!
  !! The cheaper ansatz M = H~ - (eps_n+eps_m)/2 S~ is exact only when [S, dS/dk] = 0.
  !! For the ecalj Cu sdp model it overshoots sum|v_nm|^2 by 1.5x at 10-15 eV transition
  !! energy and 2.7x at 20-30 eV (against 0.5x undershoot when dS/dk is dropped), which
  !! shows up directly as excess eps_2 above 10 eV.  Setting S=1 reduces this to get_vnm.
  !!@param   vk,out: band-basis group velocities
  !!@param   vk0,in: orbital-basis dH/dk (get_vlm0 on ham_r)
  !!@param   sk0,in: orbital-basis dS/dk (get_vlm0 on S_r)
  !!@param  Ovlk,in: k-space overlap S(k) (gen_ham on S_r)
  !!@param   eig,in: generalized eigenvalues (get_eig_mlo)
  !!@param  mrot,in: rotation matrix for real space
  !!@param   uni,in: generalized eigenvectors, S-orthonormal (get_eig_mlo)
  !!@param    Nk,in: The number of k-points
  !!@param  Norb,in: The number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,Norb
  real(c_double),intent(in),dimension(3,3):: mrot
  real(c_double),intent(in),dimension(Norb,Nk):: eig
  complex(c_double),intent(in),dimension(3,Norb,Norb,Nk):: vk0,sk0
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: Ovlk,uni
  complex(c_double),intent(out),dimension(3,Norb,Norb,Nk):: vk

  integer(c_int32_t) i,j_dir,l,m,n,k,info
  integer(c_int64_t) lwork
  real(c_double) rwork(3*Norb-2),eq(Norb),sq(Norb)
  complex(c_double) work_query(1)
  complex(c_double),allocatable :: work(:)
  complex(c_double) tmp(3,Norb,Norb)
  complex(c_double),dimension(Norb,Norb):: Umat,Xmat,Gmat,vtmp,Wtmp,tmp_j,tmp_s
  real(c_double),parameter:: ovl_thresh=1.0d-8

  call zheev('V','U',Norb,Umat,Norb,eq,work_query,-1_c_int64_t,rwork,info)
  lwork = int(dble(work_query(1)), c_int64_t)

  !$omp parallel private(tmp,Umat,Xmat,Gmat,vtmp,Wtmp,tmp_j,tmp_s, &
  !$omp                  eq,sq,rwork,work,info,j_dir,l,m,n,k)
  allocate(work(lwork))
  !$omp do
  kloop: do i=1,Nk
     ! --- S(k) = U diag(eq) U^H, then X = S^1/2 = U diag(sqrt(eq)) U^H ---
     Umat(:,:)=Ovlk(:,:,i)
     call zheev('V','U',Norb,Umat,Norb,eq,work,lwork,rwork,info)
     if(info/=0)then
        print*,'zheev failed in get_vnm_mlo (overlap): info=',info
        stop
     end if
     do n=1,Norb
        sq(n)=sqrt(max(eq(n),ovl_thresh))
     end do
     do n=1,Norb
        Wtmp(:,n)=Umat(:,n)*sq(n)
     end do
     call zgemm('N','C',Norb,Norb,Norb,(1.0d0,0.0d0),Wtmp,Norb,Umat,Norb,(0.0d0,0.0d0),Xmat,Norb)

     do j_dir=1,3
        ! --- dX/dk : (U^H dS U)_nm / (sq_n + sq_m), rotated back ---
        vtmp(:,:)=sk0(j_dir,:,:,i)
        call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),vtmp,Norb,Umat,Norb,(0.0d0,0.0d0),Wtmp,Norb)
        call zgemm('C','N',Norb,Norb,Norb,(1.0d0,0.0d0),Umat,Norb,Wtmp,Norb,(0.0d0,0.0d0),tmp_s,Norb)
        do m=1,Norb
           do n=1,Norb
              tmp_s(n,m)=tmp_s(n,m)/(sq(n)+sq(m))
           end do
        end do
        call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),Umat,Norb,tmp_s,Norb,(0.0d0,0.0d0),Wtmp,Norb)
        call zgemm('N','C',Norb,Norb,Norb,(1.0d0,0.0d0),Wtmp,Norb,Umat,Norb,(0.0d0,0.0d0),vtmp,Norb)

        ! --- G = C^H (dX) X C ---
        call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),vtmp,Norb,Xmat,Norb,(0.0d0,0.0d0),tmp_s,Norb)
        call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),tmp_s,Norb,uni(:,:,i),Norb,(0.0d0,0.0d0),Wtmp,Norb)
        call zgemm('C','N',Norb,Norb,Norb,(1.0d0,0.0d0),uni(:,:,i),Norb,Wtmp,Norb,(0.0d0,0.0d0),Gmat,Norb)

        ! --- M = C^H (dH/dk) C - G E - E G^H ---
        vtmp(:,:)=vk0(j_dir,:,:,i)
        call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),vtmp,Norb,uni(:,:,i),Norb,(0.0d0,0.0d0),Wtmp,Norb)
        call zgemm('C','N',Norb,Norb,Norb,(1.0d0,0.0d0),uni(:,:,i),Norb,Wtmp,Norb,(0.0d0,0.0d0),tmp_j,Norb)
        do n=1,Norb
           do k=1,Norb
              tmp(j_dir,k,n)=tmp_j(k,n)-Gmat(k,n)*eig(n,i)-eig(k,i)*conjg(Gmat(n,k))
           end do
        end do
     end do

     !rotate axis
     vk(:,:,:,i)=(0.0d0,0.0d0)
     band_loop5: do n=1,Norb
        do k=1,Norb
           do l=1,3
              do m=1,3
                 vk(m,k,n,i)=vk(m,k,n,i)+mrot(m,l)*dble(tmp(l,k,n))
              end do
           end do
        end do
     end do band_loop5
  end do kloop
  !$omp end do
  deallocate(work)
  !$omp end parallel
end subroutine get_vnm_mlo

subroutine gen_tr_greenw_0(trGk,wl,eig,mu,delta,Nk,Nw,Norb) bind(C)
  !> This function obtain orbital(band) trace of green function
  !!@param trGk,out: trace of green function
  !!@param    wl,in: list of energies
  !!@param   eig,in: energies at k-points
  !!@param    mu,in: chemical potential
  !!@param delta,in: dumping factor
  !!@param    Nk,in: The number of k-points
  !!@param    Nw,in: The number of energies mesh
  !!@param  Norb,in: The number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,Norb
  real(c_double),intent(in):: mu,delta
  real(c_double),intent(in),dimension(Norb,Nk):: eig
  real(c_double),intent(in),dimension(Nw):: wl
  complex(c_double),intent(out),dimension(Nw,Nk):: trGk

  integer(c_int32_t) i,j,n

  trGk(:,:) = (0.0d0, 0.0d0)

  !$omp parallel do private(i,n)
  kloop: do j=1,Nk
     !$omp simd
     wloop: do i=1,Nw
        band_loop: do n=1,Norb
           trGk(i,j)=trGk(i,j)+1./cmplx(wl(i)-eig(n,j)+mu,delta,kind=c_double)
        end do band_loop
     end do wloop
     !$omp end simd
  end do kloop
  !$omp end parallel do
end subroutine gen_tr_greenw_0

subroutine gen_dos(Dos,wl,eig,uni,mu,delta,Nk,Nw,Norb) bind(C)
  !> This function obtain partial Dos (and sum of them is total)
  !!@param  Dos,out: partial Density of states
  !!@param    wl,in: list of energies
  !!@param   eig,in: energies at k-points
  !!@param   uni,in: unitary matrix at k-points  
  !!@param    mu,in: chemical potential
  !!@param delta,in: dumping factor
  !!@param    Nk,in: The number of k-points
  !!@param    Nw,in: The number of energies mesh
  !!@param  Norb,in: The number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,Norb
  real(c_double),intent(in):: mu,delta
  real(c_double),intent(in),dimension(Norb,Nk):: eig
  real(c_double),intent(in),dimension(Nw):: wl
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(c_double),intent(out),dimension(Nw,Norb):: Dos

  integer(c_int32_t) i,j,k,n

  !$omp parallel
  !$omp workshare
  Dos(:,:)=0.0d0
  !$omp end workshare
  orb_loop: do j=1,Norb
     !$omp do private(i,k,n)
     wloop: do i=1,Nw
        kloop: do k=1,Nk
           bandloop: do n=1,Norb
              Dos(i,j)=Dos(i,j)+uni(j,n,k)*conjg(uni(j,n,k))/cmplx(wl(i)-eig(n,k)+mu,delta,kind=c_double)
           end do bandloop
        end do kloop
     end do wloop
     !$omp end do
  end do orb_loop
  !$omp workshare
  Dos(:,:)=Dos(:,:)/Nk
  !$omp end workshare
  !$omp end parallel
end subroutine gen_dos

subroutine get_parity_prop(Pmn,rvec,ham_r,Norb,Nr) bind(C)
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Norb,Nr
  real(c_double),intent(in),dimension(3,Nr):: rvec
  complex(c_double),intent(in),dimension(Norb,Norb,Nr):: ham_r
  real(c_double),intent(inout),dimension(Norb,Norb):: Pmn

  integer(c_int32_t) i,j,l,m
  real(c_double) diff_r
  real(c_double),parameter:: eps=1e-6
  real(c_double),dimension(Norb,Norb):: parity_den,parity_num
  
  !$omp parallel
  !$omp workshare
  parity_num(:,:)=0.0d0
  parity_den(:,:)=0.0d0
  !$omp end workshare
  !$omp do private(i,j,diff_r) reduction(+: parity_num,parity_den)
  ! Find inversion-partner pairs: R_j = -R_i (i.e. R_i + R_j = 0)
  ! Compare H(R) and H(-R): if Re(<H(R),H(-R)>) > 0 → even parity (+1), else odd (-1)
  do i=1,Nr
     do j=i+1,Nr
        diff_r=sum(abs(rvec(:,i)+rvec(:,j)))
        if(diff_r<eps)then      ! found partner: rvec(:,j) = -rvec(:,i)
           do l=1,Norb
              do m=1,Norb
                 parity_num(m,l)=parity_num(m,l)+dble(ham_r(m,l,j)*conjg(ham_r(m,l,i)))
                 parity_den(m,l)=parity_den(m,l)+abs(ham_r(m,l,i))**2
              end do
           end do
           exit
        end if
     end do
  end do
  !$omp end do
  !$omp do private(l,m)
  do l=1,Norb
     do m=1,Norb
        if(parity_den(m,l)>1e-10)then
           Pmn(m,l)=parity_num(m,l)/parity_den(m,l)
        else
           Pmn(m,l)=0.0d0
        end if
     end do
  end do
  !$omp end do
  !$omp do private(l,m)
  do l=1,Norb
     do m=1,Norb
        if(Pmn(m,l)>=0.0d0)then
           Pmn(m,l)=1.0d0
        else
           Pmn(m,l)=-1.0d0
        end if
     end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine get_parity_prop
