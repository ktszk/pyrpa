subroutine lin_eliash(delta,chi,Gk,uni,init_delta,Smat,Cmat,olist,prt,kmap,invk,temp,eps,&
     Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,itemax,gap_sym,arnoldi_m) bind(C)
  !> calculate linearized eliashberg equations with TRS without soc
  !!@param     delta,out: gap function
  !!@param         Gk,in: normal green function
  !!@note: invk(:,i) index validity check:
  !!        - invk(1,i) must be in range [1, Nk]
  !!        - Enable bounds checking: ifx -check=bounds or gfortran -fbounds-check
  !!@param        uni,in: unitary matrix
  !!@param init_delta,in: band basis initial gap function
  !!@param       Smat,in: S-matrix
  !!@param       Cmat,in: C-matrix
  !!@param      olist,in: property of chi footnote
  !!@param       kmap,in: property of k-points
  !!@param       invk,in: list of reverse k-points
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
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,itemax,gap_sym
  integer(int64),intent(in):: arnoldi_m  !! Krylov subspace dimension (0=power method, >0=Arnoldi)
  integer(int64),intent(in),dimension(Nchi,2):: olist
  integer(int64),intent(in),dimension(3,Nkall):: kmap,invk
  real(real64),intent(in):: temp,eps
  real(real64),intent(in),dimension(Norb):: prt
  real(real64),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  real(real64),intent(in),dimension(Nk,Norb):: init_delta
  complex(real64),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(real64),intent(out),dimension(Nk,Nw,Norb,Norb):: delta
  complex(real64),intent(inout),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(int32) i_iter,i_eig,i
  integer(int32),parameter:: eig_max=2
  logical(1) sw_pair
  real(real64) norm,inorm,weight,lambda_rq,lambda_prev,lambda_phys,lambda1,vec_err,proj
  complex(real64),dimension(Nk,Nw,Norb,Norb):: newdelta,fk,delta1

  if(gap_sym>=0)then
     sw_pair=.true.
     print'(A7)','singlet'
  else
     sw_pair=.false.
     print'(A7)','triplet'
  end if
  weight=temp/dble(Nkall)
  lambda1=0.0d0
  call get_V_delta_nsoc_flex(chi,Smat,Cmat,Nk,Nw,Nchi,sw_pair)
  print'(A15,2E16.8)','V_delta max is ',maxval(dble(chi)),maxval(aimag(chi))
  print'(A15,2E16.8)','V_delta min is ',minval(dble(chi)),minval(aimag(chi))
  if(arnoldi_m>0)then !--- Arnoldi solver
     call solve_arnoldi(arnoldi_m)
  else !--- Power method (shift + deflation)
     eigenval_loop:do i_eig=1,eig_max !solve eig_val using power method, 1st eig is usually large negative value
        call get_initial_delta(delta,init_delta,uni,kmap,invk,Nkall,Nk,Nw,Norb,gap_sym)
        if(i_eig==2)then !orthogonalize initial vector against 1st eigenvector
           call get_rayleigh(proj,delta1,delta)
           !$omp parallel workshare
           delta(:,:,:,:)=delta(:,:,:,:)-proj*delta1(:,:,:,:)
           !$omp end parallel workshare
           call get_norm(norm,delta)
           !$omp parallel workshare
           delta(:,:,:,:)=delta(:,:,:,:)*(1.0d0/norm)
           !$omp end parallel workshare
        end if
        lambda_prev=0.0d0
        iter_loop:do i_iter=1,itemax !iteration
           call mkfk_trs_nsoc(fk,Gk,delta,Nk,Nw,Norb)
           call mkdelta_nsoc(newdelta,fk,chi,Smat,Cmat,kmap,invk,prt,olist,Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,sw_pair)
           !$omp parallel workshare
           newdelta(:,:,:,:)=newdelta(:,:,:,:)*weight
           !$omp end parallel workshare
           ! explicit orthogonal deflation + shift (i_eig==2):
           ! shift by -lambda1 (=|lambda1|>0) makes most-positive eigenvalue dominant
           if(i_eig==2)then
              !$omp parallel workshare
              newdelta(:,:,:,:)=newdelta(:,:,:,:)-lambda1*delta(:,:,:,:)
              !$omp end parallel workshare
              call get_rayleigh(proj,delta1,newdelta)
              !$omp parallel workshare
              newdelta(:,:,:,:)=newdelta(:,:,:,:)-proj*delta1(:,:,:,:)
              !$omp end parallel workshare
           end if
           call get_norm(norm,newdelta)
           if(norm <= 0.0d0)then
              print*,'Error: norm <= 0 in iter_loop'
              stop
           end if
           inorm=1.0d0/norm
           call get_rayleigh(lambda_rq,delta,newdelta)
           lambda_phys=lambda_rq+lambda1
           call get_vec_err(vec_err,newdelta,delta,inorm)
           if(abs(lambda_phys)>=1.0d2 .or. abs(lambda_phys)<1.0d-6)then
              print'(I4,A13,2E16.8,A9,E12.4)',i_iter,' lambda_elsh=',lambda_phys,abs(lambda_phys-lambda_prev),'  vec_err=',vec_err
           else
              print'(I4,A13,2F12.8,A9,E12.4)',i_iter,' lambda_elsh=',lambda_phys,abs(lambda_phys-lambda_prev),'  vec_err=',vec_err
           end if
           if(lambda_phys<0.0d0 .and. vec_err<eps)exit
           if(vec_err<eps .and. abs(lambda_phys-lambda_prev)<eps*(abs(lambda_phys)+1.0d-10))exit
           lambda_prev=lambda_phys
           !$omp parallel workshare
           delta(:,:,:,:)=newdelta(:,:,:,:)*inorm
           !$omp end parallel workshare
        end do iter_loop
        if(i_eig==1)then
           if(lambda_phys>0.0d0)then
              print*,'1st eigenvalue is positive: skipping 2nd loop'
              print*,'eliash=',lambda_phys
              exit
           end if
           lambda1=lambda_rq !save first eigenvalue for shift in i_eig=2
           delta1(:,:,:,:)=newdelta(:,:,:,:)*inorm !store 1st eigenvector (normalized) for deflation
        else
           print*,'eliash=',lambda_phys
        end if
     end do eigenval_loop
     call get_norm(norm,newdelta)
     if(norm <= 0.0d0)then
        print*,'Error: final norm <= 0'
        stop
     end if
     inorm=1.0d0/norm
     !$omp parallel workshare
     delta(:,:,:,:)=newdelta(:,:,:,:)*inorm
     !$omp end parallel workshare
  end if
contains
  subroutine get_norm(norm,func)
    complex(real64),intent(in),dimension(Nk,Nw,Norb,Norb):: func
    real(real64),intent(out):: norm

    integer(int32) i,j,l,m
    real(real64) tmp

    tmp=0.0d0
    !$omp parallel
    do l=1,Norb
       do m=1,Norb
          !$omp do private(i,j),reduction(+:tmp)
          do j=1,Nw
             do i=1,Nkall
                tmp=tmp+abs(func(invk(1,i),j,m,l))*abs(func(invk(1,i),j,m,l))
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
    complex(real64),intent(in),dimension(Nk,Nw,Norb,Norb):: del,newdel
    real(real64),intent(out):: rq

    integer(int32) i,j,l,m
    real(real64) tmp

    tmp=0.0d0
    !$omp parallel
    do l=1,Norb
       do m=1,Norb
          !$omp do private(i,j),reduction(+:tmp)
          do j=1,Nw
             do i=1,Nkall
                tmp=tmp+dble(conjg(del(invk(1,i),j,m,l))*newdel(invk(1,i),j,m,l))
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
    !> Note: if both tmp_pos and tmp_neg are 0, verr becomes 0, indicating convergence
    complex(real64),intent(in),dimension(Nk,Nw,Norb,Norb):: newdel,del
    real(real64),intent(in):: inrm
    real(real64),intent(out):: verr

    integer(int32) i,j,l,m
    real(real64) tmp_pos,tmp_neg

    tmp_pos=0.0d0
    tmp_neg=0.0d0
    !$omp parallel
    do l=1,Norb
       do m=1,Norb
          !$omp do private(i,j),reduction(+:tmp_pos,tmp_neg)
          do j=1,Nw
             do i=1,Nkall
                tmp_pos=tmp_pos+abs(newdel(invk(1,i),j,m,l)*inrm-del(invk(1,i),j,m,l))**2
                tmp_neg=tmp_neg+abs(newdel(invk(1,i),j,m,l)*inrm+del(invk(1,i),j,m,l))**2
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
    !> Uses host delta1 as deflation vector and lambda1 to store the pass-1 eigenvalue.
    integer(int64),intent(in):: m_arnoldi
    integer(int32) :: m_dim,j,ii,m_act,idx,lwork,info
    real(real64) :: beta,shift
    complex(real64) :: hij
    complex(real64),allocatable :: V(:,:,:,:,:)
    complex(real64),allocatable :: H_mat(:,:),A_mat(:,:),eigenvals(:)
    complex(real64),allocatable :: VL_dum(:,:),VR(:,:),zwork(:)
    real(real64),allocatable :: rwork(:)

    m_dim=int(m_arnoldi,int32)
    lwork=max(64,2*m_dim)
    allocate(V(Nk,Nw,Norb,Norb,0:m_dim-1),H_mat(m_dim+1,m_dim))
    allocate(A_mat(m_dim,m_dim),eigenvals(m_dim))
    allocate(VL_dum(1,1),VR(m_dim,m_dim),zwork(lwork),rwork(2*m_dim))

    ! ===== Pass 1: no shift, find lambda_- =====================================
    H_mat=(0.0d0,0.0d0)
    call get_initial_delta(V(:,:,:,:,0),init_delta,uni,kmap,invk,Nkall,Nk,Nw,Norb,gap_sym)
    print'(A,I4,A)','Arnoldi pass 1 (no shift): Krylov dim =',m_dim,' ...'
    m_act=m_dim
    do j=0,m_dim-1
       call apply_op(newdelta,V(:,:,:,:,j))
       do ii=0,j
          call get_inner(hij,V(:,:,:,:,ii),newdelta)
          H_mat(ii+1,j+1)=hij
          !$omp parallel workshare
          newdelta(:,:,:,:)=newdelta(:,:,:,:)-hij*V(:,:,:,:,ii)
          !$omp end parallel workshare
       end do
       call get_norm(beta,newdelta)
       H_mat(j+2,j+1)=cmplx(beta,0.0d0,real64)
       if(beta<1.0d-14)then; m_act=j+1; exit; end if
       if(j<m_dim-1)then
          !$omp parallel workshare
          V(:,:,:,:,j+1)=newdelta(:,:,:,:)*(1.0d0/beta)
          !$omp end parallel workshare
       end if
    end do
    A_mat(1:m_act,1:m_act)=H_mat(1:m_act,1:m_act)
    call zgeev('N','V',m_act,A_mat,m_dim,eigenvals,VL_dum,1,VR,m_dim,zwork,lwork,rwork,info)
    if(info/=0)then; print*,'ZGEEV failed pass 1: info=',info; stop; end if
    print'(A)','Pass 1 Ritz values (real):'
    do ii=1,m_act
       print'(I4,F14.6)',ii,dble(eigenvals(ii))
    end do
    idx=minloc(dble(eigenvals(1:m_act)),1)
    lambda1=dble(eigenvals(idx))   !save lambda_- in host variable
    print'(A,F12.6)','  lambda_- =',lambda1
    ! early exit: if largest Ritz value >= 0.1, adopt it as lambda_+
    if(maxval(dble(eigenvals(1:m_act)))>=0.1d0)then
       idx=maxloc(dble(eigenvals(1:m_act)),1)
       print'(A,F12.6)','  eliash   =',dble(eigenvals(idx))
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
       delta(:,:,:,:)=delta(:,:,:,:)*(1.0d0/beta)
       !$omp end parallel workshare
       deallocate(V,H_mat,A_mat,eigenvals,VL_dum,VR,zwork,rwork)
       return
    end if
    ! build delta1 = Ritz vector for lambda_-, store in host delta1 for deflation
    !$omp parallel workshare
    delta1(:,:,:,:)=(0.0d0,0.0d0)
    !$omp end parallel workshare
    do ii=1,m_act
       !$omp parallel workshare
       delta1(:,:,:,:)=delta1(:,:,:,:)+VR(ii,idx)*V(:,:,:,:,ii-1)
       !$omp end parallel workshare
    end do
    call get_norm(beta,delta1)
    !$omp parallel workshare
    delta1(:,:,:,:)=delta1(:,:,:,:)*(1.0d0/beta)
    !$omp end parallel workshare

    ! ===== Pass 2: shift + deflation, find lambda_+ ============================
    shift=-lambda1   ! shift = |lambda_-| > 0; all eigenvalues become lambda_i + shift
    H_mat=(0.0d0,0.0d0)
    call get_initial_delta(V(:,:,:,:,0),init_delta,uni,kmap,invk,Nkall,Nk,Nw,Norb,gap_sym)
    ! orthogonalize initial vector against delta1
    call get_inner(hij,delta1,V(:,:,:,:,0))
    !$omp parallel workshare
    V(:,:,:,:,0)=V(:,:,:,:,0)-hij*delta1(:,:,:,:)
    !$omp end parallel workshare
    call get_norm(beta,V(:,:,:,:,0))
    !$omp parallel workshare
    V(:,:,:,:,0)=V(:,:,:,:,0)*(1.0d0/beta)
    !$omp end parallel workshare
    print'(A,F10.6,A,I4,A)','Arnoldi pass 2 (shift=',shift,'): Krylov dim =',m_dim,' ...'
    m_act=m_dim
    do j=0,m_dim-1
       ! w = (A + shift*I)*v_j - <delta1, (A+shift*I)*v_j>*delta1
       call apply_op(newdelta,V(:,:,:,:,j))
       !$omp parallel workshare
       newdelta(:,:,:,:)=newdelta(:,:,:,:)+shift*V(:,:,:,:,j)   !shift
       !$omp end parallel workshare
       call get_inner(hij,delta1,newdelta)                       !deflation
       !$omp parallel workshare
       newdelta(:,:,:,:)=newdelta(:,:,:,:)-hij*delta1(:,:,:,:)
       !$omp end parallel workshare
       do ii=0,j  !Gram-Schmidt
          call get_inner(hij,V(:,:,:,:,ii),newdelta)
          H_mat(ii+1,j+1)=hij
          !$omp parallel workshare
          newdelta(:,:,:,:)=newdelta(:,:,:,:)-hij*V(:,:,:,:,ii)
          !$omp end parallel workshare
       end do
       call get_norm(beta,newdelta)
       H_mat(j+2,j+1)=cmplx(beta,0.0d0,real64)
       if(beta<1.0d-14)then; m_act=j+1; exit; end if
       if(j<m_dim-1)then
          !$omp parallel workshare
          V(:,:,:,:,j+1)=newdelta(:,:,:,:)*(1.0d0/beta)
          !$omp end parallel workshare
       end if
    end do
    A_mat(1:m_act,1:m_act)=H_mat(1:m_act,1:m_act)
    call zgeev('N','V',m_act,A_mat,m_dim,eigenvals,VL_dum,1,VR,m_dim,zwork,lwork,rwork,info)
    if(info/=0)then; print*,'ZGEEV failed pass 2: info=',info; stop; end if
    ! physical eigenvalue = Ritz value - shift (undo the shift)
    print'(A)','Pass 2 Ritz values - physical (real):'
    do ii=1,m_act
       print'(I4,F14.6)',ii,dble(eigenvals(ii))-shift
    end do
    idx=maxloc(dble(eigenvals(1:m_act)),1)
    print'(A,F12.6)','  eliash   =',dble(eigenvals(idx))-shift
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
    delta(:,:,:,:)=delta(:,:,:,:)*(1.0d0/beta)
    !$omp end parallel workshare

    deallocate(V,H_mat,A_mat,eigenvals,VL_dum,VR,zwork,rwork)
  end subroutine solve_arnoldi

  subroutine apply_op(w_out,v_in)
    !> apply Eliashberg operator: w_out = (T/Nkall)*K*v_in; uses host fk as scratch
    complex(real64),intent(in),dimension(Nk,Nw,Norb,Norb):: v_in
    complex(real64),intent(out),dimension(Nk,Nw,Norb,Norb):: w_out
    call mkfk_trs_nsoc(fk,Gk,v_in,Nk,Nw,Norb)
    call mkdelta_nsoc(w_out,fk,chi,Smat,Cmat,kmap,invk,prt,olist,Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,sw_pair)
    !$omp parallel workshare
    w_out(:,:,:,:)=w_out(:,:,:,:)*weight
    !$omp end parallel workshare
  end subroutine apply_op

  subroutine get_inner(h,u,v)
    !> complex inner product h = 2*<u,v> summed over IBZ (factor 2 for k/-k symmetry)
    complex(real64),intent(in),dimension(Nk,Nw,Norb,Norb):: u,v
    complex(real64),intent(out):: h
    integer(int32) i,j,l,m
    real(real64) tmp_r,tmp_i
    tmp_r=0.0d0
    tmp_i=0.0d0
    !$omp parallel
    do l=1,Norb
       do m=1,Norb
          !$omp do private(i,j),reduction(+:tmp_r,tmp_i)
          do j=1,Nw
             do i=1,Nkall
                tmp_r=tmp_r+dble(conjg(u(invk(1,i),j,m,l))*v(invk(1,i),j,m,l))
                tmp_i=tmp_i+aimag(conjg(u(invk(1,i),j,m,l))*v(invk(1,i),j,m,l))
             end do
          end do
          !$omp end do
       end do
    end do
    !$omp end parallel
    h=cmplx(2.0d0*tmp_r,2.0d0*tmp_i,real64)
  end subroutine get_inner
end subroutine lin_eliash

subroutine get_V_delta_nsoc_flex(chi,Smat,Cmat,Nk,Nw,Nchi,sw_pair)
  !> This function obtains pairing interaction V_delta without soc
  !!@param  chi,inout: irreducible susceptibility and pairing interaction
  !!@param    Smat,in: S-matrix
  !!@param    Cmat,in: C-matrix
  !!@param      Nk,in: Number of k-points
  !!@param      Nw,in: Number of Matsubara frequencies
  !!@param    Nchi,in: Number of footnote of chi
  !!@param sw_pair,in: switch of singlet or triplet paring interacton
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nk,Nw,Nchi
  logical(1),intent(in):: sw_pair
  real(real64),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  complex(real64),intent(inout),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(int32) i,j,l,m,n,info
  integer(int32),dimension(Nchi):: ipiv
  complex(real64),dimension(Nchi,Nchi):: cmat1,cmat2,cmat3,cmat4,cmat5
  complex(real64),dimension(2*Nchi):: work

  wloop:do j=1,Nw
     qloop:do i=1,Nk
        !$omp parallel
        !$omp workshare
        cmat1(:,:)=0.0d0
        cmat2(:,:)=0.0d0
        !$omp end workshare
        !$omp do private(m,n)
        do l=1,Nchi
           do m=1,Nchi
              do n=1,Nchi
                 cmat1(m,l)=cmat1(m,l)-chi(i,j,m,n)*Smat(n,l) !-chi0S
                 cmat2(m,l)=cmat2(m,l)+chi(i,j,m,n)*Cmat(n,l) !chi0C
              end do
           end do
        end do
        !$omp end do
        !$omp workshare
        cmat3(:,:)=-cmat1(:,:) !chi0S
        cmat4(:,:)=cmat2(:,:)  !chi0C
        !$omp end workshare
        !$omp do 
        do l=1,Nchi
           cmat1(l,l)=cmat1(l,l)+1.0d0 !I-chi0S
           cmat2(l,l)=cmat2(l,l)+1.0d0 !I+chi0C
        end do
        !$omp end do
        !$omp end parallel
        call zgetrf(Nchi,Nchi,cmat1,Nchi,ipiv,info)
        if(info/=0)then; print*,'zgetrf failed: info=',info; stop; end if
        call zgetri(Nchi,cmat1,Nchi,ipiv,work,2*Nchi,info)
        if(info/=0)then; print*,'zgetri failed: info=',info; stop; end if
        call zgetrf(Nchi,Nchi,cmat2,Nchi,ipiv,info)
        if(info/=0)then; print*,'zgetrf failed: info=',info; stop; end if
        call zgetri(Nchi,cmat2,Nchi,ipiv,work,2*Nchi,info)
        if(info/=0)then; print*,'zgetri failed: info=',info; stop; end if
        !$omp parallel
        !$omp workshare
        cmat5(:,:)=0.0d0
        !$omp end workshare
        !$omp do private(m,n)
        do l=1,Nchi
           do m=1,Nchi
              do n=1,Nchi
                 cmat5(m,l)=cmat5(m,l)+cmat1(m,n)*cmat3(n,l) !(1-chi0S)^-1chi0S
              end do
           end do
        end do
        !$omp end do
        !$omp workshare
        cmat1(:,:)=cmat5(:,:) !cmat1=(I-chi_0S)^-1chi_0S
        cmat5(:,:)=0.0d0
        !$omp end workshare
        !$omp do private(m,n)
        do l=1,Nchi
           do m=1,Nchi
              do n=1,Nchi
                 cmat5(m,l)=cmat5(m,l)+cmat2(m,n)*cmat4(n,l) !(1+chi0C)^-1chi0C
              end do
           end do
        end do
        !$omp end do
        !$omp workshare
        cmat2(:,:)=cmat5(:,:) !cmat2=(1+chi_0C)^--1chi_0C
        !$omp end workshare
        if(sw_pair)then !singlet
           !$omp workshare
           cmat4(:,:)=0.5d0*(Smat(:,:)+Cmat(:,:)) !bare Vud=(C+S)/2
           !$omp end workshare
           !$omp do private(m,n)
           do l=1,Nchi
              do m=1,Nchi
                 do n=1,Nchi
                    cmat4(m,l)=cmat4(m,l)+1.5d0*Smat(m,n)*cmat1(n,l)-0.5d0*Cmat(m,n)*cmat2(n,l) !(C+S)/2+3/2SchisS-1/2CchicC
                 end do
              end do
           end do
           !$omp end do
        else !triplet
           !$omp workshare
           cmat4(:,:)=0.5d0*(Smat(:,:)-Cmat(:,:)) !bare Vuu=(C-S)/2
           !$omp end workshare
           !$omp do private(m,n)
           do l=1,Nchi
              do m=1,Nchi
                 do n=1,Nchi
                    cmat4(m,l)=cmat4(m,l)-0.5d0*(Smat(m,n)*cmat1(n,l)+Cmat(m,n)*cmat2(n,l)) !(C-S)/2-1/2SchisS-1/2CchicC
                 end do
              end do
           end do
           !$omp end do
        end if
        !$omp single
        chi(i,j,:,:)=cmat4(:,:)
        !$omp end single
        !$omp end parallel
     end do qloop
  end do wloop
end subroutine get_V_delta_nsoc_flex

subroutine mkfk_trs_nsoc(fk,Gk,delta,Nk,Nw,Norb) bind(C,name="mkfk_trs_nsoc_")
  !>calculate linearized anomalous green function Fk with TRS
  !>if we considder TRS,F_ab(k)=G_ac(k)Delta_cd(k)Gbd(-k)
  !> =G_ac(k)Delta_cd(k)G^*_dbss'(k)ss'
  !!@param   fk,out: anomalous green function
  !!@param    Gk,in: normal green function
  !!@param delta,in: gap function
  !!@param  invk,in: list of reverse k-points
  !!@param Nkall,in: Number of all k-points
  !!@param    Nk,in: Number of irreducible k-points
  !!@param    Nw,in: Number of Matsubara frequencies
  !!@param  Norb,in: Number of orbitals
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nk,Nw,Norb
  complex(real64),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(real64),intent(in),dimension(Nk,Nw,Norb,Norb):: delta
  complex(real64),intent(out),dimension(Nk,Nw,Norb,Norb):: fk
 
  integer(int32) i,j,l,m,n
  complex(real64),dimension(Nk,Nw,Norb,Norb):: cmat1
  !$omp parallel
  !$omp workshare
  cmat1(:,:,:,:)=0.0d0
  fk(:,:,:,:)=0.0d0
  !$omp end workshare
  do l=1,Norb
     do m=1,Norb
        do n=1,Norb
           !$omp do private(i,j)
           do j=1,Nw
              do i=1,Nk
                 cmat1(i,j,m,l)=cmat1(i,j,m,l)+Gk(i,j,m,n)*delta(i,j,n,l)
              end do
           end do
           !$omp end do
        end do
     end do
  end do

  do l=1,Norb
     do m=l,Norb
        do n=1,Norb
           !$omp do private(i,j)
           do j=1,Nw
              do i=1,Nk
                 fk(i,j,m,l)=fk(i,j,m,l)-cmat1(i,j,m,n)*conjg(Gk(i,j,l,n))
              end do
           end do
           !$omp end do
        end do
        if(l/=m)then !if TRS and no SOC delta^+(k,iw)=delta(k,iw)
           !$omp do private(i,j)
           do j=1,Nw
              !$omp simd
              do i=1,Nk
                 fk(i,j,l,m)=conjg(fk(i,j,m,l))
              end do
              !$omp end simd
           end do
           !$omp end do
        else
           !$omp do private(i,j)
           do j=1,Nw
              !$omp simd
              do i=1,Nk
                 fk(i,j,l,l)=dble(fk(i,j,l,l))
              end do
              !$omp end simd
           end do
           !$omp end do
        end if
     end do
  end do
  !$omp end parallel
end subroutine mkfk_trs_nsoc

subroutine mkdelta_nsoc(newdelta,delta,Vdelta,Smat,Cmat,kmap,invk,prt,olist,Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,sw_pair)
  !> This function obtain new gap function Delta=sum VF
  !!@param newdelta,out: new gap function
  !!@param     delta,in: anomalous green function
  !!@param    Vdelta,in: irreducible susceptibility and pairing interaction
  !!@param      Smat,in: S-matrix
  !!@param      Cmat,in: C-matrix
  !!@param      kmap,in: property of k-point
  !!@param      invk,in: list of reverse k-points
  !!@param     olist,in: property of chi footnote
  !!@param     Nkall,in: Number of all k-points
  !!@param        Nk,in: Number of k-points
  !!@param        Nw,in: Number of Matsubara frequencies
  !!@param      Nchi,in: Number of footnote of chi
  !!@param      Norb,in: Number of orbitals
  !!@param        Nx,in: Number of kx mesh
  !!@param        Ny,in: Number of ky mesh
  !!@param        Nz,in: Number of kz mesh
  !!@param sw_pair,in: switch of singlet or triplet paring interacton
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz
  integer(int64),intent(in),dimension(3,Nkall):: kmap,invk
  integer(int64),intent(in),dimension(Nchi,2):: olist
  logical(1),intent(in):: sw_pair
  real(real64),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  real(real64),intent(in),dimension(Norb):: prt
  complex(real64),intent(in),dimension(Nk,Nw,Nchi,Nchi):: Vdelta
  complex(real64),intent(in),dimension(Nk,Nw,Norb,Norb):: delta
  complex(real64),intent(out),dimension(Nk,Nw,Norb,Norb):: newdelta
  
  integer(int32) i,j,k,n,l,m
  real(real64) sgn,dprt
  complex(real64),dimension(0:Nx-1,0:Ny-1,0:Nz-1,2*Nw):: tmpVdelta,tmpfk,tmp
  logical(1),parameter:: sw_odd_freq=.false.

  if(sw_pair)then
     sgn=1.0d0 !singlet Fk=F^T-k
     dprt=1.0d0 !gap parity even freq singlet is even
  else
     sgn=-1.0d0 !triplet Fk=-F^T-k
     dprt=-1.0d0 !gap parity even freq triplet is odd
  end if

  if(sw_odd_freq)then !consider odd frequency
     dprt=-dprt
  end if
  !$omp parallel workshare
  newdelta(:,:,:,:)=0.0d0
  !$omp end parallel workshare
  do l=1,Nchi
     do m=1,Nchi
        if(olist(l,1)>=olist(m,2))then
           !$omp parallel
           !$omp do private(i)
           do i=1,Nkall
              if(invk(2,i)==0)then !k,0
                 tmpVdelta(kmap(1,i),kmap(2,i),kmap(3,i),1)=Vdelta(invk(1,i),1,m,l) !j=1 corresponds to w_n=0
              else if(invk(2,i)==1)then !-k,0
                 tmpVdelta(kmap(1,i),kmap(2,i),kmap(3,i),1)=Vdelta(invk(1,i),1,l,m) !V^ml(-q,0)=V^lm(q,0)
              end if
              tmpVdelta(kmap(1,i),kmap(2,i),kmap(3,i),Nw+1)=0.5d0*(Cmat(m,l)+sgn*Smat(m,l)) !Nw+1 consider w_n=>inf limit
           end do
           !$omp end do
           !$omp do private(i,j)
           do j=2,Nw
              do i=1,Nkall
                 if(invk(2,i)==0)then !k,iw,k,-iw
                    tmpVdelta(kmap(1,i),kmap(2,i),kmap(3,i),j)=Vdelta(invk(1,i),j,m,l)
                    tmpVdelta(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+2)=conjg(Vdelta(invk(1,i),j,l,m)) !V^ml(q,-iw)=V^*lm(q,iw)
                 else if(invk(2,i)==1)then !-k,iw, -k,-iw
                    tmpVdelta(kmap(1,i),kmap(2,i),kmap(3,i),j)=Vdelta(invk(1,i),j,l,m) !V^ml(-q,iw)=V^lm(q,iw)
                    tmpVdelta(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+2)=conjg(Vdelta(invk(1,i),j,m,l)) !V^ml(-q,-iw)=V^*ml(q,iw)
                 end if
              end do
           end do
           !$omp end do
           !$omp do private(i,j)
           do j=1,Nw
              do i=1,Nkall
                 if(invk(2,i)==0)then !k
                    tmpfk(kmap(1,i),kmap(2,i),kmap(3,i),j)=delta(invk(1,i),j,olist(l,2),olist(m,1))
                    tmpfk(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=sgn*dprt*prt(olist(l,2))*prt(olist(m,1))*&
                         conjg(delta(invk(1,i),j,olist(l,2),olist(m,1))) !F(k,-w)=sgn*dprt*p*p*F^*(k,w) (no soc only)
                 else if(invk(2,i)==1)then !-k
                    tmpfk(kmap(1,i),kmap(2,i),kmap(3,i),j)=prt(olist(l,2))*prt(olist(m,1))*dprt*delta(invk(1,i),j,olist(l,2),olist(m,1)) !F(-k,w)=p*p*dprt*F(k,w)
                    tmpfk(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=sgn*delta(invk(1,i),j,olist(m,1),olist(l,2)) !F(-k,-w)=sgn*F(k,w)
                 end if
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
                    !$omp simd
                    do i=0,Nx-1
                       tmp(i,j,k,n)=tmpVdelta(i,j,k,n)*tmpfk(i,j,k,n)
                    end do
                    !$omp end simd
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
                 if(invk(2,i)==0)then
                    newdelta(invk(1,i),j,olist(l,1),olist(m,2))=newdelta(invk(1,i),j,olist(l,1),olist(m,2))&
                         +tmp(kmap(1,i),kmap(2,i),kmap(3,i),j)
                 end if
              end do
           end do
           !$omp end parallel do
        end if
     end do
  end do

  !$omp parallel
  do l=1,Norb
     do m=l,Norb
        !$omp do private(i,j)
        do j=1,Nw
           do i=1,Nk
              if(l/=m)then
                 newdelta(i,j,l,m)=conjg(newdelta(i,j,m,l))
              else
                 newdelta(i,j,l,l)=dble(newdelta(i,j,l,l))
              end if
           end do
        end do
        !$omp end do
     end do
  end do
  !$omp end parallel
end subroutine mkdelta_nsoc

subroutine get_initial_delta(delta,init_delta,uni,kmap,invk,Nkall,Nk,Nw,Norb,gap_sym)
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
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  use constant
  implicit none
  integer(int64),intent(in):: Nw,Norb,Nkall,Nk,gap_sym
  integer(int64),intent(in),dimension(3,Nkall):: kmap,invk
  real(real64),intent(in),dimension(Nk,Norb):: init_delta
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(real64),intent(out),dimension(Nk,Nw,Norb,Norb):: delta

  integer(int32) i,j,l,m,n
  real(real64) norm,sgn
  complex(real64),dimension(Nk,Norb,Norb):: delta0

  if(gap_sym==0)then
     !$omp parallel
     !$omp workshare
     delta(:,:,:,:)=0.0d0
     !$omp end workshare
     do l=1,Norb
        !$omp workshare
        delta(:,:,l,l)=1.0d0
        !$omp end workshare
     end do
     !$omp end parallel
  else
     if(gap_sym>0)then
        sgn=1.0d0
     else
        sgn=-1.0d0
     end if
     !$omp parallel
     !$omp workshare
     delta0(:,:,:)=0.0d0
     !$omp end workshare
     do l=1,Norb
        do m=l,Norb
           !$omp do private(i,n)
           do i=1,Nk
              do n=1,Norb
                    delta0(i,m,l)=delta0(i,m,l)+uni(m,n,i)*init_delta(i,n)*conjg(uni(l,n,i))
              end do
           end do
           !$omp end do
           !$omp do private(i)
           do i=1,Nk
              if(l/=m)then !TRS gap function is Hermite
                 delta0(i,l,m)=conjg(delta0(i,m,l))
              else !TRS diagonal gap is real
                 delta0(i,l,l)=dble(delta0(i,l,l))
              end if
           end do
           !$omp end do
        end do
     end do

     do l=1,Norb
        do m=1,Norb
           !$omp do private(i,j)
           do j=1,Nw
              !$omp simd
              do i=1,Nk
                 delta(i,j,m,l)=delta0(i,m,l)
              end do
              !$omp end simd
           end do
           !$omp end do
        end do
     end do
     !$omp end parallel
  end if

  norm=0.0d0
  !$omp parallel
  do l=1,Norb
     do m=1,Norb
        !$omp do private(j,i) reduction(+:norm)
        do j=1,Nw
           do i=1,Nkall
              norm=norm+abs(delta(invk(1,i),j,m,l))*abs(delta(invk(1,i),j,m,l))
           end do
        end do
        !$omp end do
     end do
  end do
  !$omp end parallel
  ! norm=0 check
  if(norm <= 0.0d0)then
     print*,'Error: get_initial_delta norm <= 0'
     stop
  end if
  !$omp parallel
  !$omp workshare
  delta(:,:,:,:)=delta(:,:,:,:)/sqrt(2.0d0*norm)
  !$omp end workshare
  !$omp end parallel
end subroutine get_initial_delta

subroutine conv_delta_orb_to_band(deltab,delta,uni,prt,invk,Norb,Nkall,Nk,Nw,gap_sym) bind(C)
  !> This function transform orbital basis gap function into band basis gap function
  !!@param deltab,out: band basis gap function
  !!@param   delta,in: orbital basis gap function
  !!@param     uni,in: unitary matrix
  !!@param    invk,in: footnote of reverse k
  !!@param    Norb,in: Number of orbital
  !!@param   Nkall,in: Number of all k-points
  !!@param      Nk,in: Number of irreducible k-points
  !!oparam      Nw,in: Number of Matsubara frequencies
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nw,Norb,Nk,Nkall,gap_sym
  integer(int64),intent(in),dimension(3,Nkall):: invk
  real(real64),intent(in),dimension(Norb):: prt
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(real64),intent(in),dimension(Nk,Nw,Norb,Norb):: delta
  complex(real64),intent(out),dimension(Nkall,Norb,Norb):: deltab

  integer(int32) i,l,m,n
  real(real64) dprt
  complex(real64),dimension(Nkall,Norb,Norb):: tmp

  if(gap_sym<0)then
     dprt=-1.0d0 !odd
  else
     dprt=1.0d0 !even
  end if

  !$omp parallel
  !$omp workshare
  deltab(:,:,:)=0.0d0
  tmp(:,:,:)=0.0d0
  !$omp end workshare
  do l=1,Norb
     do m=1,Norb
        do n=1,Norb
           !$omp do private(i)
           do i=1,Nkall
              if(invk(2,i)==0)then
                 tmp(i,m,l)=tmp(i,m,l)+conjg(uni(n,m,invk(1,i)))*delta(invk(1,i),1,n,l)
              else if(invk(2,i)==1)then
                 tmp(i,m,l)=tmp(i,m,l)+uni(n,m,invk(1,i))*dprt*prt(n)*prt(l)*delta(invk(1,i),1,n,l)
              end if
           end do
           !$omp end do
        end do
     end do
  end do

  do l=1,Norb
     do m=1,Norb
        do n=1,Norb
           !$omp do private(i)
           do i=1,Nkall
              if(invk(2,i)==0)then
                 deltab(i,m,l)=deltab(i,m,l)+tmp(i,m,n)*uni(n,l,invk(1,i))
              else if(invk(2,i)==1)then
                 deltab(i,m,l)=deltab(i,m,l)+tmp(i,m,n)*conjg(uni(n,l,invk(1,i)))
              end if
           end do
           !$omp end do
        end do
     end do
  end do
  !$omp end parallel
end subroutine conv_delta_orb_to_band

subroutine remap_delta(delta,delta0,prt,invk,Nkall,Nk,Nw,Norb,gap_sym) bind(C,name='remap_delta_')
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nkall,Nk,Nw,Norb,gap_sym
  integer(int64),intent(in),dimension(3,Nkall):: invk
  real(real64),intent(in),dimension(Norb):: prt
  complex(real64),intent(in),dimension(Nk,Nw,Norb,Norb):: delta0
  complex(real64),intent(out),dimension(Nkall,Nw,Norb,Norb):: delta

  integer(int32) i,j,l,m
  real(real64) dprt
  if(gap_sym<0)then
     dprt=-1.0d0
  else
     dprt=1.0d0
  end if
  
  !$omp parallel
  do l=1,Norb
     do m=1,Norb
        !$omp do private(i,j)
        do j=1,Nw
           do i=1,Nkall
              if(invk(2,i)==0)then
                 delta(i,j,m,l)=delta0(invk(1,i),j,m,l)
              else if(invk(2,i)==1)then
                 delta(i,j,m,l)=dprt*prt(m)*prt(l)*delta0(invk(1,i),j,m,l)
              end if
           end do
        end do
        !$omp end do
     end do
  end do
  !$omp end parallel
end subroutine remap_delta
