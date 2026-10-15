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
  !!@param    gap_sym,in: gap symmetry number
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
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
  integer(c_int64_t) Ndk
  integer(c_int64_t),dimension(Nkall):: idk,ipk
  real(c_double) norm,inorm,norm2,weight,lambda_rq,vec_err,pp
  logical su2
  complex(c_double),allocatable,dimension(:,:,:,:):: delta,newdelta,fk

  lambda_out=0.0d0
  ! Δ is stored on Ndk points: idk(i) = stored index of full k-point i, ipk(i)=1 if Δ(k_i) is
  ! rebuilt by parity, Δ(-k) = pp·P Δ(k) P  (pp=+1 singlet, -1 triplet)
  call get_delta_kmap_soc(idk,ipk,Ndk,invk,Nkall,Nk,gap_sym)
  pp=merge(1.0d0,-1.0d0,gap_sym>0)
  allocate(delta(Ndk,Nw,Norb,Norb),newdelta(Ndk,Nw,Norb,Norb),fk(Ndk,Nw,Norb,Norb))
  ! Without SOC (G spin-diagonal and spin-independent) singlet and triplet are separate sectors:
  ! the kernel keeps the seed's spin structure, but round-off would grow the other sector, so it
  ! is projected out (project_spin).  With SOC they mix and only the parity is imposed.
  su2=is_spin_rotation_invariant()
  if(gap_sym/=0 .and. gap_sym/=-10)then
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
  call get_V_soc_flex(chi,Vmat,sgnsig2,Nk,Nw,Nchi)
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
  if(arnoldi_m>0)then !--- Arnoldi solver
     call solve_arnoldi(arnoldi_m)
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
        call get_norm(norm,newdelta)
        if(norm <= 0.0d0)then
           print*,'Error: norm <= 0 in iter_loop'
           stop
        end if
        inorm=1.0d0/norm
        ! Rayleigh quotient: lambda_rq = Re(<delta, newdelta>) assuming ||delta||=1
        ! shift is +norm2, so true eigenvalue = lambda_rq - norm2
        call get_rayleigh(lambda_rq,delta,newdelta)
        ! Vector convergence: ||newdelta/norm - delta||
        call get_vec_err(vec_err,newdelta,delta,inorm)
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
  call get_norm(norm,newdelta)
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
    !> spurious TRIM component.  Skipped for parity mixing (gap_sym=-10), stored on the full BZ.
    complex(c_double),intent(inout),dimension(Ndk,Nw,Norb,Norb):: d
    integer(c_int32_t) ii,a,b
    if(gap_sym==0 .or. gap_sym==-10)return
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
    if(.not.su2 .or. gap_sym==0 .or. gap_sym==-10)return
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

  subroutine get_norm(norm,func)
    complex(c_double),intent(in),dimension(Ndk,Nw,Norb,Norb):: func
    real(c_double),intent(out):: norm

    integer(c_int32_t) i,j,l,m
    real(c_double) tmp

    tmp=0.0d0
    !$omp parallel
    do l=1,Norb
       do m=1,Norb
          !$omp do private(i,j),reduction(+:tmp)
          do j=1,Nw
             do i=1,Nkall
                tmp=tmp+abs(func(idk(i),j,m,l))*abs(func(idk(i),j,m,l))
             end do
          end do
          !$omp end do
       end do
    end do
    !$omp end parallel
    norm=sqrt(2.0d0*tmp)
  end subroutine get_norm

  subroutine get_rayleigh(rq,del,newdel)
    !> Rayleigh quotient: rq = Re(<del, newdel>) assuming ||del||=1
    complex(c_double),intent(in),dimension(Ndk,Nw,Norb,Norb):: del,newdel
    real(c_double),intent(out):: rq

    integer(c_int32_t) i,j,l,m
    real(c_double) tmp

    tmp=0.0d0
    !$omp parallel
    do l=1,Norb
       do m=1,Norb
          !$omp do private(i,j),reduction(+:tmp)
          do j=1,Nw
             do i=1,Nkall
                tmp=tmp+dble(conjg(del(idk(i),j,m,l))*newdel(idk(i),j,m,l))
             end do
          end do
          !$omp end do
       end do
    end do
    !$omp end parallel
    rq=2.0d0*tmp
  end subroutine get_rayleigh

  subroutine get_vec_err(verr,newdel,del,inrm)
    !> Vector convergence: min(||newdel*inrm - del||, ||newdel*inrm + del||)
    !> Taking the minimum handles sign-flipping convergence for negative eigenvalues
    complex(c_double),intent(in),dimension(Ndk,Nw,Norb,Norb):: newdel,del
    real(c_double),intent(in):: inrm
    real(c_double),intent(out):: verr

    integer(c_int32_t) i,j,l,m
    real(c_double) tmp_pos,tmp_neg

    tmp_pos=0.0d0
    tmp_neg=0.0d0
    !$omp parallel
    do l=1,Norb
       do m=1,Norb
          !$omp do private(i,j),reduction(+:tmp_pos,tmp_neg)
          do j=1,Nw
             do i=1,Nkall
                tmp_pos=tmp_pos+abs(newdel(idk(i),j,m,l)*inrm-del(idk(i),j,m,l))**2
                tmp_neg=tmp_neg+abs(newdel(idk(i),j,m,l)*inrm+del(idk(i),j,m,l))**2
             end do
          end do
          !$omp end do
       end do
    end do
    !$omp end parallel
    verr=sqrt(2.0d0*min(tmp_pos,tmp_neg))
  end subroutine get_vec_err

  subroutine solve_arnoldi(m_arnoldi)
    !> Two-pass Arnoldi: pass 1 (no shift) finds lambda_-, pass 2 (shift+deflation) finds lambda_+.
    !>
    !> Algorithm sketch (Arnoldi factorization of operator A := apply_op):
    !>   - Build orthonormal Krylov basis V_m = [v_0, v_1, ..., v_{m-1}] with v_0 = initial Δ.
    !>   - Each step: w = A v_j, modified Gram-Schmidt against V(:,0..j) -> H_mat(:,j+1),
    !>                beta = ||w_residual||, v_{j+1} = w/beta.
    !>   - Result satisfies A V_m = V_m H_m + beta v_m e_m^T  (H_m: m x m upper Hessenberg).
    !>   - Diagonalize H_m via ZGEEV; eigenvalues = Ritz values (approx eigenvalues of A),
    !>     Ritz vector = V_m * y where H_m y = theta y.
    !>   - Pass 2 replaces A with (A + shift*I) and projects out delta1 each step,
    !>     making lambda_+ the dominant Ritz value.
    integer(c_int64_t),intent(in):: m_arnoldi
    integer(c_int32_t) :: m_dim,ip,j,ii,m_act,idx,lwork,info
    real(c_double) :: beta,shift,lambda1
    complex(c_double) :: hij
    ! V         : Krylov basis (last index = basis vector number 0..m-1)
    ! delta1    : Pass-1 Ritz vector for lambda_-, used as deflation direction in Pass 2
    ! H_mat     : (m+1) x m upper Hessenberg matrix from Arnoldi factorization
    ! A_mat     : working copy of H_mat(1:m,1:m) for ZGEEV (destroyed on output)
    ! eigenvals : Ritz values (eigenvalues of A_mat)
    ! VR        : right eigenvectors of A_mat (coefficients of Ritz vectors in V basis)
    ! m_act     : actual Krylov dimension reached (may be < m_dim on lucky breakdown)
    complex(c_double),allocatable :: V(:,:,:,:,:),delta1(:,:,:,:)
    complex(c_double),allocatable :: H_mat(:,:),A_mat(:,:),eigenvals(:)
    complex(c_double),allocatable :: VL_dum(:,:),VR(:,:),zwork(:)
    real(c_double),allocatable :: rwork(:)

    m_dim=int(m_arnoldi,c_int32_t)
    lwork=max(64,2*m_dim)
    allocate(V(Ndk,Nw,Norb,Norb,0:m_dim-1),delta1(Ndk,Nw,Norb,Norb))
    allocate(H_mat(m_dim+1,m_dim))
    allocate(A_mat(m_dim,m_dim),eigenvals(m_dim))
    allocate(VL_dum(1,1),VR(m_dim,m_dim),zwork(lwork),rwork(2*m_dim))

    ! ===== Pass 1: no shift, find lambda_- =====================================
    H_mat=(0.0d0,0.0d0)
    ! v_0: normalized initial gap function (already unit-norm from get_initial_delta_soc)
    call get_initial_delta_soc(V(:,:,:,:,0),init_delta,uni,kmap,slist,invk,invs,idk,ipk,Nkall,Nk,Ndk,Nw,Norb,gap_sym)
    print'(A,I4,A)','Arnoldi pass 1 (no shift): Krylov dim =',m_dim,' ...'
    m_act=m_dim
    ! ----- Arnoldi iteration: extend the Krylov basis one vector per step -----
    do j=0,m_dim-1
       call apply_op(newdelta,V(:,:,:,:,j))                ! w = A v_j
       ! modified Gram-Schmidt: orthogonalize w against v_0..v_j
       do ii=0,j
          call get_inner(hij,V(:,:,:,:,ii),newdelta)        ! H_{ii,j} = <v_ii, w>
          H_mat(ii+1,j+1)=hij
          !$omp parallel workshare
          newdelta(:,:,:,:)=newdelta(:,:,:,:)-hij*V(:,:,:,:,ii)  ! w := w - H_{ii,j} v_ii
          !$omp end parallel workshare
       end do
       call get_norm(beta,newdelta)                         ! beta = ||w_residual||
       H_mat(j+2,j+1)=cmplx(beta,0.0d0,c_double)            ! sub-diagonal entry H_{j+1,j}
       ! lucky breakdown: residual is essentially zero -> K_{j+1} is invariant
       if(beta<1.0d-14)then; m_act=j+1; exit; end if
       if(j<m_dim-1)then
          !$omp parallel do
          do ip=1,int(Ndk,c_int32_t)
             V(ip,:,:,:,j+1)=newdelta(ip,:,:,:)*(1.0d0/beta)   ! v_{j+1} = w / beta
          end do
          !$omp end parallel do
       end if
    end do
    ! ----- diagonalize the m_act x m_act Hessenberg block -----
    A_mat(1:m_act,1:m_act)=H_mat(1:m_act,1:m_act)           ! ZGEEV destroys input, copy first
    call zgeev('N','V',m_act,A_mat,m_dim,eigenvals,VL_dum,1,VR,m_dim,zwork,lwork,rwork,info)
    if(info/=0)then; print*,'ZGEEV failed pass 1: info=',info; stop; end if
    print'(A)','Pass 1 Ritz values (real):'
    do ii=1,m_act
       print'(I4,F14.6)',ii,dble(eigenvals(ii))
    end do
    idx=minloc(dble(eigenvals(1:m_act)),1)                  ! pick most negative Ritz value
    lambda1=dble(eigenvals(idx))
    call real_ritz(VR(1:m_act,idx),m_act)
    print'(A,F12.6)','  lambda_- =',lambda1
    ! early exit: if Pass 1 already resolved a sufficiently positive Ritz value,
    ! it is the physical eliashberg eigenvalue lambda_+ — skip Pass 2.
    if(maxval(dble(eigenvals(1:m_act)))>=0.1d0)then
       idx=maxloc(dble(eigenvals(1:m_act)),1)
       print'(A,F12.6)','  eliash   =',dble(eigenvals(idx))
       lambda_out=dble(eigenvals(idx))
       call real_ritz(VR(1:m_act,idx),m_act)
       ! reconstruct Ritz vector x = sum_ii VR(ii,idx) * v_{ii-1}
       !$omp parallel workshare
       delta(:,:,:,:)=(0.0d0,0.0d0)
       !$omp end parallel workshare
       do ii=1,m_act
          !$omp parallel workshare
          delta(:,:,:,:)=delta(:,:,:,:)+VR(ii,idx)*V(:,:,:,:,ii-1)
          !$omp end parallel workshare
       end do
       call get_norm(beta,delta)
       !$omp parallel workshare
       delta(:,:,:,:)=delta(:,:,:,:)*(1.0d0/beta)            ! normalize result
       !$omp end parallel workshare
       deallocate(V,delta1,H_mat,A_mat,eigenvals,VL_dum,VR,zwork,rwork)
       return
    end if
    ! build delta1 = Ritz vector for lambda_- (= V_m * y_-), store for deflation in Pass 2
    !$omp parallel do
    do ip=1,int(Ndk,c_int32_t)
       delta1(ip,:,:,:)=(0.0d0,0.0d0)
    end do
    !$omp end parallel do
    do ii=1,m_act
       !$omp parallel do
       do ip=1,int(Ndk,c_int32_t)
          delta1(ip,:,:,:)=delta1(ip,:,:,:)+VR(ii,idx)*V(ip,:,:,:,ii-1)
       end do
       !$omp end parallel do
    end do
    call get_norm(beta,delta1)
    !$omp parallel do
    do ip=1,int(Ndk,c_int32_t)
       delta1(ip,:,:,:)=delta1(ip,:,:,:)*(1.0d0/beta)         ! normalize deflation vector
    end do
    !$omp end parallel do

    ! ===== Pass 2: shift + deflation, find lambda_+ ============================
    ! Operator: B := (A + shift*I) restricted to span{delta1}^perp.
    ! Spectrum mapping: eigenvalue lambda_i of A -> lambda_i + shift for B.
    ! With shift = -lambda_-, lambda_+ + shift becomes the largest-magnitude eigenvalue,
    ! so Arnoldi converges to it preferentially.
    shift=-lambda1   ! shift = |lambda_-| > 0; all eigenvalues become lambda_i + shift
    H_mat=(0.0d0,0.0d0)
    call get_initial_delta_soc(V(:,:,:,:,0),init_delta,uni,kmap,slist,invk,invs,idk,ipk,Nkall,Nk,Ndk,Nw,Norb,gap_sym)
    ! orthogonalize initial vector against delta1, then renormalize -> v_0 in delta1^perp
    call get_inner(hij,delta1,V(:,:,:,:,0))
    !$omp parallel do
    do ip=1,int(Ndk,c_int32_t)
       V(ip,:,:,:,0)=V(ip,:,:,:,0)-hij*delta1(ip,:,:,:)
    end do
    !$omp end parallel do
    call get_norm(beta,V(:,:,:,:,0))
    !$omp parallel do
    do ip=1,int(Ndk,c_int32_t)
       V(ip,:,:,:,0)=V(ip,:,:,:,0)*(1.0d0/beta)
    end do
    !$omp end parallel do
    print'(A,F10.6,A,I4,A)','Arnoldi pass 2 (shift=',shift,'): Krylov dim =',m_dim,' ...'
    m_act=m_dim
    ! ----- Arnoldi iteration with operator B = P_{delta1^perp} (A + shift*I) -----
    do j=0,m_dim-1
       ! Step (1): w = A v_j
       call apply_op(newdelta,V(:,:,:,:,j))
       ! Step (2): w := w + shift * v_j  (apply spectral shift)
       !$omp parallel workshare
       newdelta(:,:,:,:)=newdelta(:,:,:,:)+shift*V(:,:,:,:,j)   !shift
       !$omp end parallel workshare
       ! Step (3): project out delta1 (deflate lambda_- direction) -> w in delta1^perp
       call get_inner(hij,delta1,newdelta)                       !deflation
       !$omp parallel workshare
       newdelta(:,:,:,:)=newdelta(:,:,:,:)-hij*delta1(:,:,:,:)
       !$omp end parallel workshare
       ! Step (4): modified Gram-Schmidt against v_0..v_j (same as Pass 1)
       do ii=0,j  !Gram-Schmidt
          call get_inner(hij,V(:,:,:,:,ii),newdelta)
          H_mat(ii+1,j+1)=hij
          !$omp parallel workshare
          newdelta(:,:,:,:)=newdelta(:,:,:,:)-hij*V(:,:,:,:,ii)
          !$omp end parallel workshare
       end do
       call get_norm(beta,newdelta)
       H_mat(j+2,j+1)=cmplx(beta,0.0d0,c_double)
       if(beta<1.0d-14)then; m_act=j+1; exit; end if   ! lucky breakdown
       if(j<m_dim-1)then
          !$omp parallel do
          do ip=1,int(Ndk,c_int32_t)
             V(ip,:,:,:,j+1)=newdelta(ip,:,:,:)*(1.0d0/beta)        ! v_{j+1} = w / beta
          end do
          !$omp end parallel do
       end if
    end do
    ! ----- diagonalize the deflated/shifted Hessenberg matrix -----
    A_mat(1:m_act,1:m_act)=H_mat(1:m_act,1:m_act)
    call zgeev('N','V',m_act,A_mat,m_dim,eigenvals,VL_dum,1,VR,m_dim,zwork,lwork,rwork,info)
    if(info/=0)then; print*,'ZGEEV failed pass 2: info=',info; stop; end if
    ! physical eigenvalue = Ritz value - shift (undo the shift to get original spectrum of A)
    print'(A)','Pass 2 Ritz values - physical (real):'
    do ii=1,m_act
       print'(I4,F14.6)',ii,dble(eigenvals(ii))-shift
    end do
    idx=maxloc(dble(eigenvals(1:m_act)),1)                ! largest Ritz value -> lambda_+ + shift
    print'(A,F12.6)','  eliash   =',dble(eigenvals(idx))-shift
    lambda_out=dble(eigenvals(idx))-shift
    call real_ritz(VR(1:m_act,idx),m_act)
    ! reconstruct Ritz vector for lambda_+: delta = V_m * y_+
    !$omp parallel workshare
    delta(:,:,:,:)=(0.0d0,0.0d0)
    !$omp end parallel workshare
    do ii=1,m_act
       !$omp parallel workshare
       delta(:,:,:,:)=delta(:,:,:,:)+VR(ii,idx)*V(:,:,:,:,ii-1)
       !$omp end parallel workshare
    end do
    call get_norm(beta,delta)
    !$omp parallel workshare
    delta(:,:,:,:)=delta(:,:,:,:)*(1.0d0/beta)               ! normalize gap function
    !$omp end parallel workshare
    ! y above is an eigenvector of the deflated P A P (P = 1 - delta1 delta1^T), not of A:
    ! A y = lambda_+ y + c delta1 with c = <delta1, A y>, so x = y + c/(lambda_+ - lambda_-) delta1
    call apply_op(newdelta,delta)
    call get_inner(hij,delta1,newdelta)
    if(abs(lambda_out-lambda1)>1.0d-12)then     ! lambda_+ = lambda_- (e.g. V=0): nothing to undo
       hij=hij/(lambda_out-lambda1)
       !$omp parallel workshare
       delta(:,:,:,:)=delta(:,:,:,:)+hij*delta1(:,:,:,:)
       !$omp end parallel workshare
       call get_norm(beta,delta)
       !$omp parallel workshare
       delta(:,:,:,:)=delta(:,:,:,:)*(1.0d0/beta)
       !$omp end parallel workshare
    end if

    deallocate(V,delta1,H_mat,A_mat,eigenvals,VL_dum,VR,zwork,rwork)
  end subroutine solve_arnoldi

  subroutine apply_op(w_out,v_in)
    !> apply Eliashberg operator: w_out = (T/Nkall)*K*v_in; uses host fk as scratch
    complex(c_double),intent(in),dimension(Ndk,Nw,Norb,Norb):: v_in
    complex(c_double),intent(out),dimension(Ndk,Nw,Norb,Norb):: w_out
    call mkfk_trs_soc(fk,Gk,v_in,sgnsig,slist,invk,invs,idk,ipk,Nkall,Nk,Ndk,Nw,Norb,gap_sym)
    call mkdelta_soc(w_out,fk,chi,Vmat,sgnsig,sgnsig2,prt,kmap,invk,invs,invschi,olist,slist,idk,ipk,&
         Nkall,Nk,Ndk,Nw,Nchi,Norb,Nx,Ny,Nz,gap_sym)
    call project_parity(w_out)
    call project_spin(w_out)
    !$omp parallel workshare
      w_out(:,:,:,:)=-w_out(:,:,:,:)*weight
    !$omp end parallel workshare
  end subroutine apply_op

  subroutine get_inner(h,u,v)
    !> real inner product h = 2*Re<u,v> summed over full BZ.  mkfk_trs_soc/mkdelta_soc fill
    !> components by conjugation (TRS: D_uu=-D_dd^+, D_ud Hermitian), so apply_op is only
    !> real-linear and Arnoldi must run over the real vector space it acts on.
    complex(c_double),intent(in),dimension(Ndk,Nw,Norb,Norb):: u,v
    complex(c_double),intent(out):: h
    integer(c_int32_t) i,j,l,m
    real(c_double) tmp_r,tmp_i
    tmp_r=0.0d0
    tmp_i=0.0d0
    !$omp parallel
    do l=1,Norb
       do m=1,Norb
          !$omp do private(i,j),reduction(+:tmp_r,tmp_i)
          do j=1,Nw
             do i=1,Nkall
                tmp_r=tmp_r+dble(conjg(u(idk(i),j,m,l))*v(idk(i),j,m,l))
                tmp_i=tmp_i+aimag(conjg(u(idk(i),j,m,l))*v(idk(i),j,m,l))
             end do
          end do
          !$omp end do
       end do
    end do
    !$omp end parallel
    h=cmplx(2.0d0*tmp_r,0.0d0,c_double)
  end subroutine get_inner

  subroutine real_ritz(y,n)
    !> Ritz coefficients y of the real Hessenberg matrix: fix the arbitrary ZGEEV phase so that
    !> y is real (exact for a real Ritz value), keeping the Ritz vector in the TRS-real space
    integer(c_int32_t),intent(in):: n
    complex(c_double),intent(inout):: y(n)
    complex(c_double) ph
    ph=y(maxloc(abs(y),1))
    y(:)=cmplx(dble(y(:)*conjg(ph)/abs(ph)),0.0d0,c_double)
  end subroutine real_ritz
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
                       sp=merge(pp*prt(olist(l,2))*prt(olist(m,1)),1.0d0,ipk(i)==1)        !F(k) from stored point
                       sm=merge(pp*prt(olist(l,2))*prt(olist(m,1)),1.0d0,ipk(invk(3,i))==1) !F(-k) from stored point
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

  integer(c_int64_t) i

  if(gap_sym/=0 .and. gap_sym/=-10)then
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
