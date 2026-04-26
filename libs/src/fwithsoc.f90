subroutine lin_eliash_soc(delta,chi,Gk,uni,init_delta,Vmat,sgnsig,sgnsig2,prt,olist,slist,kmap,invk,invs,invschi,&
     temp,eps,Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,itemax,gap_sym,arnoldi_m) bind(C)
  !> calculate linearized eliashberg equations with soc and TRS
  !!@param     delta,out: gap function
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
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,itemax,gap_sym
  integer(int64),intent(in):: arnoldi_m  !! Krylov subspace dimension (0=power method, >0=Arnoldi)
  integer(int64),intent(in),dimension(Nchi,2):: olist
  integer(int64),intent(in),dimension(Norb):: slist,invs
  integer(int64),intent(in),dimension(3,Nkall):: kmap,invk
  integer(int64),intent(in),dimension(Nchi):: invschi
  real(real64),intent(in):: temp,eps
  real(real64),intent(in),dimension(Nchi,Nchi):: Vmat
  real(real64),intent(in),dimension(Norb,Norb):: sgnsig
  real(real64),intent(in),dimension(Nchi,Nchi):: sgnsig2
  real(real64),intent(in),dimension(Norb):: prt
  real(real64),intent(in),dimension(Nk,Norb):: init_delta
  complex(real64),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(real64),intent(out),dimension(Nkall,Nw,Norb,Norb):: delta
  complex(real64),intent(inout),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(int32) i_iter,i_eig,count,i
  integer(int32),parameter:: eig_max=2
  real(real64) norm,inorm,norm2,weight,lambda_rq,vec_err
  complex(real64),dimension(Nkall,Nw,Norb,Norb):: newdelta,fk

  weight=temp/dble(Nkall)
  ! norm2 is the shift used in the power iteration: newdelta = -weight*K*delta + norm2*delta
  ! converges to (λ + norm2)*delta, so true eigenvalue λ = Rayleigh quotient - norm2.
  ! After the 1st eigenvalue is found, norm2 is set to its norm for deflation (shift-deflation).
  norm2=0.0d0
  call get_V_soc_flex(chi,Vmat,sgnsig2,Nk,Nw,Nchi)
  print'(A15,2E16.8)','V_delta max is ',maxval(dble(chi)),maxval(aimag(chi))
  print'(A15,2E16.8)','V_delta min is ',minval(dble(chi)),minval(aimag(chi))
  if(arnoldi_m>0)then !--- Arnoldi solver
     call solve_arnoldi(arnoldi_m)
  else !--- Power method
  eigenval_loop:do i_eig=1,eig_max !solve eig_val using power method, 1st eig is usually large negative value
     call get_initial_delta_soc(delta,init_delta,uni,kmap,slist,invk,invs,Nkall,Nk,Nw,Norb,gap_sym)
     count=0 !count too small eigenvalue
     iter_loop:do i_iter=1, itemax !iteration
        call mkfk_trs_soc(fk,Gk,delta,sgnsig,slist,invk,invs,Nkall,Nk,Nw,Norb,gap_sym)
        call mkdelta_soc(newdelta,fk,chi,Vmat,sgnsig,sgnsig2,kmap,invk,invs,invschi,olist,slist,Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,gap_sym)
        !$omp parallel workshare
        newdelta(:,:,:,:)=-newdelta(:,:,:,:)*weight+delta(:,:,:,:)*norm2
        !$omp end parallel workshare
        call get_norm(norm,newdelta)
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
contains
  subroutine get_norm(norm,func)
    complex(real64),intent(in),dimension(Nkall,Nw,Norb,Norb):: func
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
                tmp=tmp+abs(func(i,j,m,l))*abs(func(i,j,m,l))
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
    complex(real64),intent(in),dimension(Nkall,Nw,Norb,Norb):: del,newdel
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
                tmp=tmp+dble(conjg(del(i,j,m,l))*newdel(i,j,m,l))
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
    complex(real64),intent(in),dimension(Nkall,Nw,Norb,Norb):: newdel,del
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
                tmp_pos=tmp_pos+abs(newdel(i,j,m,l)*inrm-del(i,j,m,l))**2
                tmp_neg=tmp_neg+abs(newdel(i,j,m,l)*inrm+del(i,j,m,l))**2
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
    integer(int64),intent(in):: m_arnoldi
    integer(int32) :: m_dim,ip,j,ii,m_act,idx,lwork,info
    real(real64) :: beta,shift,lambda1
    complex(real64) :: hij
    complex(real64),allocatable :: V(:,:,:,:,:),delta1(:,:,:,:)
    complex(real64),allocatable :: H_mat(:,:),A_mat(:,:),eigenvals(:)
    complex(real64),allocatable :: VL_dum(:,:),VR(:,:),zwork(:)
    real(real64),allocatable :: rwork(:)

    m_dim=int(m_arnoldi,int32)
    lwork=max(64,2*m_dim)
    allocate(V(Nkall,Nw,Norb,Norb,0:m_dim-1),delta1(Nkall,Nw,Norb,Norb))
    allocate(H_mat(m_dim+1,m_dim))
    allocate(A_mat(m_dim,m_dim),eigenvals(m_dim))
    allocate(VL_dum(1,1),VR(m_dim,m_dim),zwork(lwork),rwork(2*m_dim))

    ! ===== Pass 1: no shift, find lambda_- =====================================
    H_mat=(0.0d0,0.0d0)
    call get_initial_delta_soc(V(:,:,:,:,0),init_delta,uni,kmap,slist,invk,invs,Nkall,Nk,Nw,Norb,gap_sym)
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
          !$omp parallel do
          do ip=1,int(Nkall,int32)
             V(ip,:,:,:,j+1)=newdelta(ip,:,:,:)*(1.0d0/beta)
          end do
          !$omp end parallel do
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
    lambda1=dble(eigenvals(idx))
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
       deallocate(V,delta1,H_mat,A_mat,eigenvals,VL_dum,VR,zwork,rwork)
       return
    end if
    ! build delta1 = Ritz vector for lambda_-, store for deflation
    !$omp parallel do
    do ip=1,int(Nkall,int32)
       delta1(ip,:,:,:)=(0.0d0,0.0d0)
    end do
    !$omp end parallel do
    do ii=1,m_act
       !$omp parallel do
       do ip=1,int(Nkall,int32)
          delta1(ip,:,:,:)=delta1(ip,:,:,:)+VR(ii,idx)*V(ip,:,:,:,ii-1)
       end do
       !$omp end parallel do
    end do
    call get_norm(beta,delta1)
    !$omp parallel do
    do ip=1,int(Nkall,int32)
       delta1(ip,:,:,:)=delta1(ip,:,:,:)*(1.0d0/beta)
    end do
    !$omp end parallel do

    ! ===== Pass 2: shift + deflation, find lambda_+ ============================
    shift=-lambda1   ! shift = |lambda_-| > 0; all eigenvalues become lambda_i + shift
    H_mat=(0.0d0,0.0d0)
    call get_initial_delta_soc(V(:,:,:,:,0),init_delta,uni,kmap,slist,invk,invs,Nkall,Nk,Nw,Norb,gap_sym)
    ! orthogonalize initial vector against delta1
    call get_inner(hij,delta1,V(:,:,:,:,0))
    !$omp parallel do
    do ip=1,int(Nkall,int32)
       V(ip,:,:,:,0)=V(ip,:,:,:,0)-hij*delta1(ip,:,:,:)
    end do
    !$omp end parallel do
    call get_norm(beta,V(:,:,:,:,0))
    !$omp parallel do
    do ip=1,int(Nkall,int32)
       V(ip,:,:,:,0)=V(ip,:,:,:,0)*(1.0d0/beta)
    end do
    !$omp end parallel do
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
          !$omp parallel do
          do ip=1,int(Nkall,int32)
             V(ip,:,:,:,j+1)=newdelta(ip,:,:,:)*(1.0d0/beta)
          end do
          !$omp end parallel do
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

    deallocate(V,delta1,H_mat,A_mat,eigenvals,VL_dum,VR,zwork,rwork)
  end subroutine solve_arnoldi

  subroutine apply_op(w_out,v_in)
    !> apply Eliashberg operator: w_out = (T/Nkall)*K*v_in; uses host fk as scratch
    complex(real64),intent(in),dimension(Nkall,Nw,Norb,Norb):: v_in
    complex(real64),intent(out),dimension(Nkall,Nw,Norb,Norb):: w_out
    call mkfk_trs_soc(fk,Gk,v_in,sgnsig,slist,invk,invs,Nkall,Nk,Nw,Norb,gap_sym)
    call mkdelta_soc(w_out,fk,chi,Vmat,sgnsig,sgnsig2,kmap,invk,invs,invschi,olist,slist,Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,gap_sym)
    !$omp parallel workshare
      w_out(:,:,:,:)=-w_out(:,:,:,:)*weight
    !$omp end parallel workshare
  end subroutine apply_op

  subroutine get_inner(h,u,v)
    !> complex inner product h = 2*<u,v> summed over full BZ
    complex(real64),intent(in),dimension(Nkall,Nw,Norb,Norb):: u,v
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
                tmp_r=tmp_r+dble(conjg(u(i,j,m,l))*v(i,j,m,l))
                tmp_i=tmp_i+aimag(conjg(u(i,j,m,l))*v(i,j,m,l))
             end do
          end do
          !$omp end do
       end do
    end do
    !$omp end parallel
    h=cmplx(2.0d0*tmp_r,2.0d0*tmp_i,real64)
  end subroutine get_inner
end subroutine lin_eliash_soc

subroutine get_chi0_soc(chi,sgnsig,sgnsig2,invschi,Vmat,Gk,kmap,invk,invs,olist,slist,temp,&
     Nx,Ny,Nz,Nw,Nk,Nkall,Nchi,Norb) bind(C)
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz
  integer(int64),intent(in),dimension(Nchi,2):: olist
  integer(int64),intent(in),dimension(Norb):: slist,invs
  integer(int64),intent(in),dimension(3,Nkall):: kmap,invk
  integer(int64),intent(out),dimension(Nchi):: invschi
  real(real64),intent(in):: temp
  real(real64),intent(in),dimension(Nchi,Nchi):: Vmat
  real(real64),intent(out),dimension(Norb,Norb):: sgnsig
  real(real64),intent(out),dimension(Nchi,Nchi):: sgnsig2
  complex(real64),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(real64),intent(out),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(int32),dimension(Nchi,Nchi,2)::chi_map
  integer(int32),dimension(Nchi*(Nchi+1)/2,2)::irr_chi

  call get_invschi()
  call get_sgnsig()
  call get_chi_map_soc(chi_map,irr_chi,olist,invs,Nchi,Norb)
  call get_chi0_conv_soc(chi,Gk,kmap,invk,invs,irr_chi,chi_map,olist,sgnsig,sgnsig2,temp,Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi)
  call ckchi()
contains
  subroutine get_sgnsig()
    ! sgnsig(i,j)  = s_i * s_j  (orbital spin signature for TRS G(-k) relation)
    ! sgnsig2(i,j) = s_{i1}*s_{i2}*s_{j1}*s_{j2}  (chi index spin signature)
    integer(int32) i,j
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

    subroutine ckchi()
    integer(int32) i,l,m,n,info,chik,chikall,iorb
    real(real64) maxchi0,maxchi02
    real(real64),dimension(2*Nchi):: rwork
    complex(real64),dimension(Nchi*Nchi*4+1):: work
    complex(real64),dimension(Nchi):: eigc
    complex(real64),dimension(Nchi,Nchi):: chi0,tmp,tmp1

    maxchi02=-1.0d5
    do i=1,Nk
       !$omp parallel
       !$omp workshare
       chi0(:,:)=0.0d0
       !$omp end workshare
       !$omp do private(l,m,n)
       do l=1,Nchi
          do m=1,Nchi
             do n=1,Nchi
                chi0(m,l)=chi0(m,l)-chi(i,1,m,n)*Vmat(n,l)
             end do
          end do
       end do
       !$omp end do
       !$omp end parallel
       call zgeev('N','N',Nchi,chi0,Nchi,eigc,tmp1,Nchi,tmp,Nchi,work,Nchi*Nchi*4+1,rwork,info)
       if(info/=0)then; print*,'zgeev failed: info=',info; stop; end if
       maxchi0=maxval(dble(eigc))
       if(maxchi0>maxchi02)then
          chik=i
          maxchi02=maxchi0
       end if
    end do
    do i=1,Nkall
       if(invk(2,i)==0 .and. invk(1,i)==chik)then
          chikall=i
          exit
       end if
    end do
    print'(A7,3I4,F12.8)','SDW/CDW',kmap(:,chikall),maxchi02
  end subroutine ckchi

  subroutine get_invschi()
    integer(int32) l,m,i,j
    invschi(:)=0
    do l=1,Nchi
       i=invs(olist(l,1))
       j=invs(olist(l,2))
       do m=1,Nchi
          if(i==olist(m,1) .and. j==olist(m,2))then
             invschi(l)=m
             exit
          end if
       end do
    end do
  end subroutine get_invschi
end subroutine get_chi0_soc

subroutine get_V_soc_flex(chi,Vmat,sgnsig2,Nk,Nw,Nchi)
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nk,Nw,Nchi
  real(real64),intent(in),dimension(Nchi,Nchi)::sgnsig2
  real(real64),intent(in),dimension(Nchi,Nchi):: Vmat
  complex(real64),intent(inout),dimension(Nk,Nw,Nchi,Nchi):: chi
  
   integer(int32) i,j,l,info
  integer(int32),dimension(Nchi):: ipiv
   complex(real64),dimension(Nchi,Nchi):: cmat1,cmat2,cmat3,cmatv,neg_cmatv

   cmatv(:,:)=cmplx(Vmat(:,:),0.0d0,kind=real64)
   neg_cmatv(:,:)=-cmatv(:,:)

  ! SOC-RPA pairing vertex: V_Δ = -V + V·χ,  χ = (I+χ₀·V)^{-1}·χ₀·V
  !$omp parallel do collapse(2) private(i,cmat1,cmat2,cmat3,ipiv,info,l)
  wloop:do j=1,Nw
     qloop:do i=1,Nk
        call zgemm('N','N',Nchi,Nchi,Nchi,(1.0d0,0.0d0),chi(i,j,:,:),Nchi,cmatv,Nchi,(0.0d0,0.0d0),cmat1,Nchi) !chi0V
        cmat2(:,:)=cmat1(:,:) !chi0V (RHS for zgesv)
        do l=1,Nchi
           cmat1(l,l)=cmat1(l,l)+1.0d0 !I+chi0V
        end do
        call zgesv(Nchi,Nchi,cmat1,Nchi,ipiv,cmat2,Nchi,info) !(I+chi0V)X=chi0V -> cmat2=chi (RPA)
        if(info/=0)then; print*,'zgesv failed: info=',info; stop; end if
        cmat3(:,:)=neg_cmatv(:,:)   !start from -V
        call zgemm('N','N',Nchi,Nchi,Nchi,(1.0d0,0.0d0),cmatv,Nchi,cmat2,Nchi,(1.0d0,0.0d0),cmat3,Nchi) !-V + V*chi
        chi(i,j,:,:)=cmat3(:,:)
     end do qloop
  end do wloop
  !$omp end parallel do
end subroutine get_V_soc_flex

subroutine mkfk_trs_soc(fk,Gk,delta,sgnsig,slist,invk,invs,Nkall,Nk,Nw,Norb,gap_sym) bind(C,name="mkfk_trs_soc_")
  !
  !>calculate linearized anomalous green function Fk with TRS
  !>if we considder TRS,F_ab(k)=G_ac(k)Delta_cd(k)Gbd(-k)
  !> =G_ac(k)Delta_cd(k)G^*_dbss'(k)ss'
  !
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nkall,Nk,Nw,Norb,gap_sym
  integer(int64),intent(in),dimension(3,Nkall):: invk
  integer(int64),intent(in),dimension(Norb):: invs,slist
  real(real64),intent(in),dimension(Norb,Norb):: sgnsig
  complex(real64),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(real64),intent(in),dimension(Nkall,Nw,Norb,Norb):: delta
  complex(real64),intent(out),dimension(Nkall,Nw,Norb,Norb):: fk
   
  integer(int32) i,j,l,m,n
  complex(real64),dimension(Norb,Norb):: gmat1,gmat2,cmat1

  !$omp parallel do collapse(2) private(i,j,l,m,n,gmat1,gmat2,cmat1)
  do j=1,Nw
     do i=1,Nkall
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
        call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),gmat1,Norb,delta(i,j,:,:),Norb,(0.0d0,0.0d0),cmat1,Norb)
        call zgemm('N','C',Norb,Norb,Norb,(-1.0d0,0.0d0),cmat1,Norb,gmat2,Norb,(0.0d0,0.0d0),fk(i,j,:,:),Norb)
     end do
  end do
  !$omp end parallel do
  
  do l=1,Norb
     if(slist(l)==1)then
        do m=1,Norb
           !$omp do private(i,j)
           do j=1,Nw
              do i=1,Nkall
                 if(slist(m)==1)then
                    if(gap_sym<0)then
                       fk(i,j,invs(l),invs(m))=-conjg(fk(i,j,m,l))
                    else
                       fk(i,j,m,l)=0.0d0
                       fk(i,j,invs(l),invs(m))=0.0d0
                    end if
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
  do l=1,Norb
     if(slist(l)==1)then
        do m=1,Norb
           if(slist(m)==-1)then
              !$omp do private(i,j)
              do j=1,Nw
                 do i=1,Nkall
                    if(gap_sym==-10)then
                       continue
                    else if(gap_sym<0)then
                       fk(i,j,invs(m),invs(l))=fk(i,j,m,l)
                    else
                       fk(i,j,invs(m),invs(l))=-fk(i,j,m,l)
                    end if
                 end do
              end do
              !$omp end do
           end if
        end do
     end if
  end do
end subroutine mkfk_trs_soc
  
subroutine mkdelta_soc(newdelta,delta,Vdelta,Vmat,sgnsig,sgnsig2,kmap,invk,invs,invschi,olist,slist,Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,gap_sym)
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,gap_sym
  integer(int64),intent(in),dimension(3,Nkall):: kmap,invk
  integer(int64),intent(in),dimension(Nchi,2):: olist
  integer(int64),intent(in),dimension(Norb):: slist,invs
  integer(int64),intent(in),dimension(Nchi):: invschi
  real(real64),intent(in),dimension(Norb,Norb):: sgnsig
  real(real64),intent(in),dimension(Nchi,Nchi):: Vmat,sgnsig2
  complex(real64),intent(in),dimension(Nk,Nw,Nchi,Nchi):: Vdelta
  complex(real64),intent(in),dimension(Nkall,Nw,Norb,Norb):: delta
  complex(real64),intent(out),dimension(Nkall,Nw,Norb,Norb):: newdelta

  integer(int32) i,j,k,n,l,m
  integer(int64),parameter::  par_mix=-10
  complex(real64),dimension(0:Nx-1,0:Ny-1,0:Nz-1,2*Nw):: tmpVdelta,tmpfk,tmp
  
  !$omp parallel workshare
  newdelta(:,:,:,:)=0.0d0
  !$omp end parallel workshare
  do l=1,Nchi
     do m=1,Nchi
        ! triplet (gap_sym<0): all spin components needed;
        ! singlet (gap_sym>0): only ↑↓/↓↑ components (same-spin pairs vanish for singlet)
        if((gap_sym<0) .or. (gap_sym>0 .and. slist(olist(l,1))*slist(olist(m,2))<0))then !triplet gap need all D, singlet gap only Dud(du)
           if((gap_sym==par_mix .and. slist(olist(l,1))*slist(olist(m,2))<0) .or. slist(olist(l,1))==1)then !calc only Duu, Dud
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
                    tmpVdelta(kmap(1,i),kmap(2,i),kmap(3,i),Nw+1)=-Vmat(m,l) !Nw+1 consider w_n=>inf limit
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
                 !$omp do private(i,j)
                 do j=1,Nw
                    do i=1,Nkall
                       tmpfk(kmap(1,i),kmap(2,i),kmap(3,i),j)=delta(i,j,olist(l,2),olist(m,1))
                       tmpfk(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=-delta(invk(3,i),j,olist(m,1),olist(l,2)) !F(k,-ω)=-F(-k,ω): anomalous GF Matsubara symmetry
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
                       newdelta(i,j,olist(l,1),olist(m,2))=newdelta(i,j,olist(l,1),olist(m,2))&
                            +tmp(kmap(1,i),kmap(2,i),kmap(3,i),j)
                    end do
                 end do
                 !$omp end parallel do
              end if
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
              do i=1,Nkall
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
  do l=1,Norb
     if (slist(l)==1)then
        do m=1,Norb
           if(slist(m)==-1)then
              !$omp do private(j,i)
              do j=1,Nw
                 do i=1,Nkall
                    if(gap_sym>0)then !if system has parity, singlet Delta_ud(k)=-Delta_du(k)
                       newdelta(i,j,invs(l),invs(m))=-newdelta(i,j,l,m)
                    else if(gap_sym==par_mix)then !parity mixing
                       continue
                    else !if system has parity, triplet Delta_ud(k)=Delta_du(k)
                       newdelta(i,j,invs(l),invs(m))=newdelta(i,j,l,m)
                    end if
                 end do
              end do
              !$omp end do
           end if
        end do
     end if
  end do
  !$omp end parallel
end subroutine mkdelta_soc

subroutine get_chi0_conv_soc(chi,Gk,kmap,invk,invs,irr_chi,chi_map,olist,&
     sgnsig,sgnsig2,temp,Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi)
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nw,Norb,Nchi,Nkall,Nk,Nx,Ny,Nz
  integer(int64),intent(in),dimension(Nchi,2):: olist
  integer(int64),intent(in),dimension(Norb):: invs
  integer(int64),intent(in),dimension(3,Nkall):: kmap,invk
  integer(int32),intent(in),dimension(Nchi,Nchi,2):: chi_map
  integer(int32),intent(in),dimension(Nchi*(Nchi+1)/2,2):: irr_chi
  real(real64),intent(in):: temp
  real(real64),intent(in),dimension(Norb,Norb):: sgnsig
  real(real64),intent(in),dimension(Nchi,Nchi):: sgnsig2
  complex(real64),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(real64),intent(out),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(int32) i,j,k,l,m,n,iorb
  integer(int32) ii(0:Nx-1),ij(0:Ny-1),ik(0:Nz-1),iw(2*Nw)
  real(real64) weight
  real(real64),parameter:: eps=1.0d-9
  complex(real64),dimension(0:Nx-1,0:Ny-1,0:Nz-1,2*Nw):: tmp,tmpgk13,tmpgk42
   
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
  orb_loop:do iorb=1,Nchi*(Nchi+1)/2
     l=irr_chi(iorb,1)
     m=irr_chi(iorb,2)
     !use symmetry G^lm(k,iw)=G^ml(k,-iw) from Hermitian symmetry of Hamiltonian
     !$omp parallel do private(i,j)
     w_loop_Gk_to_tmp:do j=1,Nw
        k_loop_Gk_to_tmp:do i=1,Nkall
           if(invk(2,i)==0)then !irredicible k
              !iw
              tmpgk13(kmap(1,i),kmap(2,i),kmap(3,i),j)=Gk(invk(1,i),j,olist(l,1),olist(m,1)) !G13(k,iw)
              tmpgk42(kmap(1,i),kmap(2,i),kmap(3,i),j)=Gk(invk(1,i),j,olist(m,2),olist(l,2)) !G42(k,iw)
              !-iw
              tmpgk13(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=conjg(Gk(invk(1,i),j,olist(m,1),olist(l,1))) !G13(k,-iw)=G^*31(k,iw)
              tmpgk42(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=conjg(Gk(invk(1,i),j,olist(l,2),olist(m,2))) !G42(k,-iw)=G^*24(k,iw)
           else if(invk(2,i)==1)then !-k
              !iw
              tmpgk13(kmap(1,i),kmap(2,i),kmap(3,i),j)=sgnsig(olist(m,1),olist(l,1))&
                   *Gk(invk(1,i),j,invs(olist(m,1)),invs(olist(l,1))) !G13ss'(-k,iw)=ss'G^31-s'-s(k,iw)
              tmpgk42(kmap(1,i),kmap(2,i),kmap(3,i),j)=sgnsig(olist(l,2),olist(m,2))&
                   *Gk(invk(1,i),j,invs(olist(l,2)),invs(olist(m,2))) !G42ss'(-k,iw)=ss'G^24-s's(k,iw)
              !-iw
              tmpgk13(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=sgnsig(olist(l,1),olist(m,1))&
                   *conjg(Gk(invk(1,i),j,invs(olist(l,1)),invs(olist(m,1)))) !G13ss'(-k,-iw)=ss'G^*13-s-s'(k,iw)
              tmpgk42(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=sgnsig(olist(m,2),olist(l,2))&
                   *conjg(Gk(invk(1,i),j,invs(olist(m,2)),invs(olist(l,2)))) !G42(-k,-iw)=ss'G^*42-s-s'(k,iw)
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
              chi(invk(1,i),j,m,l)=tmp(kmap(1,i),kmap(2,i),kmap(3,i),j)*weight
              if(abs(dble(chi(invk(1,i),j,m,l)))<eps) chi(invk(1,i),j,m,l)=cmplx(0.0d0,imag(chi(invk(1,i),j,m,l)),kind=real64)
              if(abs(imag(chi(invk(1,i),j,m,l)))<eps) chi(invk(1,i),j,m,l)=cmplx(dble(chi(invk(1,i),j,m,l)),0.0d0,kind=real64)
           end if
        end do k_loop_tmp_to_chi
     end do w_loop_tmp_to_chi
     !$omp end parallel do
     chi(:,:,chi_map(m,l,1),chi_map(m,l,2))=sgnsig2(m,l)*conjg(chi(:,:,m,l))
  end do orb_loop
end subroutine get_chi0_conv_soc

subroutine get_initial_delta_soc(delta,init_delta,uni,kmap,slist,invk,invs,Nkall,Nk,Nw,Norb,gap_sym)
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
  integer(int64),intent(in),dimension(Norb):: slist,invs
  real(real64),intent(in),dimension(Nk,Norb):: init_delta
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(real64),intent(out),dimension(Nkall,Nw,Norb,Norb):: delta

   integer(int32) i,j,k,l,m,n
  real(real64) norm
   complex(real64),dimension(Norb,Norb):: matL,matR,tmpm

  !$omp parallel do private(i,k,l,m,n,matL,matR,tmpm)
  do i=1,Nkall
     k=invk(1,i)
     if(invk(2,i)==0)then
        matL=uni(:,:,k)
        do n=1,Norb
           matL(:,n)=matL(:,n)*init_delta(k,n)
        end do
        do l=1,Norb
           matR(:,l)=slist(l)*conjg(uni(invs(l),:,k))
        end do
        call zgemm('N','N',Norb,Norb,Norb,(-1.0d0,0.0d0),matL,Norb,matR,Norb,(0.0d0,0.0d0),delta(i,1,:,:),Norb)
     else if(invk(2,i)==1)then
        do m=1,Norb
           matL(m,:)=slist(m)*conjg(uni(invs(m),:,k))
        end do
        do n=1,Norb
           matL(:,n)=matL(:,n)*init_delta(k,n)
        end do
        matR=transpose(uni(:,:,k))
        call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),matL,Norb,matR,Norb,(0.0d0,0.0d0),delta(i,1,:,:),Norb)
     else
        delta(i,1,:,:)=0.0d0
     end if

     if(gap_sym>0)then
        do l=1,Norb
           do m=1,Norb
              if(slist(m)*slist(l)>=0) delta(i,1,m,l)=0.0d0
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
              norm=norm+abs(delta(i,j,m,l))*abs(delta(i,j,m,l))
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

subroutine get_chi_map_soc(chi_map,irr_chi,olist,invs,Nchi,Norb)
  !> This function generate index of exchange symmetry chi1234(q,iw)=sgn*chi*4321(q,iw).
  !> This symmetry can use system has TRS iwth soc.
  !!@param chi_map,out: mapping list of chi index
  !!@param irr_chi,out: irreducible index of chii
  !!@param    olist,in: orbital index list of chi index
  !!@param     invs,in: spin flipped index of Hamiltonian
  !!@param     Nchi,in: Number of chi index
  !!@param     Norb,in: Number of orbitals
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nchi,Norb
  integer(int64),intent(in),dimension(Nchi,2):: olist
  integer(int64),intent(in),dimension(Norb):: invs
  integer(int32),intent(out),dimension(Nchi,Nchi,2):: chi_map
  integer(int32),intent(out),dimension(Nchi*(Nchi+1)/2,2):: irr_chi

  integer(int32) l1,m1,l2,m2,iter
  integer(int32),dimension(Nchi,Nchi):: chi_irr
  integer(int32),dimension(4):: tmp1,tmp2
  logical(int32) ck

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
           tmp2(3)=invs(olist(l2,2)) !2 of 4321 with spin flip
           tmp2(4)=invs(olist(l2,1)) !1 of 4321 with spin flip
           do m2=1,Nchi
              tmp2(1)=invs(olist(m2,2)) !4 of 4321 with spin flip
              tmp2(2)=invs(olist(m2,1)) !3 of 4321 with spin flip
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
end subroutine get_chi_map_soc

subroutine get_chis_chic_soc(chic,chiszz,chispm,chi,Vmat,orb_list,olist,slist,invs,Nk,Nw,Nchi,Norb) bind(C)
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nk,Nw,Nchi,Norb
  integer(int64),intent(in),dimension(Nchi,2):: olist
  integer(int64),intent(in),dimension(Norb):: slist,invs
  integer(int64),dimension(Nchi):: orb_list
  real(real64),intent(in),dimension(Nchi,Nchi):: Vmat
  complex(real64),intent(in),dimension(Nk,Nw,Nchi,Nchi):: chi
  complex(real64),intent(out),dimension(Nk,Nchi/4,Nchi/4):: chiszz,chic,chispm

   integer(int32) i,l,m,info
  integer(int32),dimension(Nchi):: ipiv
   complex(real64),dimension(Nchi,Nchi):: cmat1,cmat2,cmatv

   cmatv(:,:)=cmplx(Vmat(:,:),0.0d0,kind=real64)

  qloop:do i=1,Nk
     call zgemm('N','N',Nchi,Nchi,Nchi,(1.0d0,0.0d0),chi(i,1,:,:),Nchi,cmatv,Nchi,(0.0d0,0.0d0),cmat1,Nchi) !chi0V
     cmat2(:,:)=chi(i,1,:,:) !RHS for zgesv = chi0
     do l=1,Nchi
        cmat1(l,l)=cmat1(l,l)+1.0d0 !I-chi0V
     end do
     call zgesv(Nchi,Nchi,cmat1,Nchi,ipiv,cmat2,Nchi,info) !(I-chi0V)X=chi0 -> cmat2=(I-chi0V)^-1chi0
     if(info/=0)then; print*,'zgesv failed: info=',info; stop; end if
     do l=1,Nchi
        do m=1,Nchi
           if(slist(olist(l,1))==slist(olist(m,1)) .and. slist(olist(l,2))==slist(olist(m,2)))then
              if(slist(olist(l,1))==-slist(olist(l,2)) .and. slist(olist(l,1))==1)then
                 chispm(i,orb_list(l),orb_list(m))=chispm(i,orb_list(l),orb_list(m))+cmat2(m,l)
              end if
           end if
           if(slist(olist(l,1))==slist(olist(l,2)) .and. slist(olist(m,1))==slist(olist(m,2)))then
              chic(i,orb_list(l),orb_list(m))=chic(i,orb_list(l),orb_list(m))+cmat2(m,l)
              if(slist(olist(l,1))==slist(olist(m,1)))then
                 chiszz(i,orb_list(l),orb_list(m))=chiszz(i,orb_list(l),orb_list(m))+cmat2(m,l)
              else
                 chiszz(i,orb_list(l),orb_list(m))=chiszz(i,orb_list(l),orb_list(m))-cmat2(m,l)
              end if
           end if
        end do
     end do
  end do qloop
  chic(:,:,:)=chic(:,:,:)*0.5d0
  chiszz(:,:,:)=chiszz(:,:,:)*0.25
end subroutine get_chis_chic_soc

subroutine conv_delta_orb_to_band_soc(deltab,delta,uni,invk,invs,slist,Norb,Nkall,Nk,Nw) bind(C)
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
  integer(int64),intent(in):: Nw,Norb,Nk,Nkall
  integer(int64),intent(in),dimension(Norb):: slist,invs
  integer(int64),intent(in),dimension(3,Nkall):: invk
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(real64),intent(in),dimension(Nkall,Nw,Norb,Norb):: delta
  complex(real64),intent(out),dimension(Nkall,Norb,Norb):: deltab

   integer(int32) i,k,l,n
   complex(real64),dimension(Norb,Norb):: matL,matR,tmpm

  !$omp parallel do private(i,k,l,n,matL,matR,tmpm)
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

     call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),matL,Norb,delta(i,1,:,:),Norb,(0.0d0,0.0d0),tmpm,Norb)
     call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),tmpm,Norb,matR,Norb,(0.0d0,0.0d0),deltab(i,:,:),Norb)
  end do
  !$omp end parallel do
end subroutine conv_delta_orb_to_band_soc

subroutine mkself_soc(sigmak,mu,Vmat,kmap,invk,invs,olist,slist,hamk,eig,uni,mu_init,&
     rfill,temp,scf_loop,pp,eps,Nkall,Nk,Nw,Norb,Nchi,Nx,Ny,Nz,sub_sigma,sw_out,&
     sw_in,m_diis,sw_rescale) bind(C)
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nw,Norb,Nchi,Nkall,Nk,Nx,Ny,Nz,scf_loop,m_diis
  integer(int64),intent(in),dimension(Nchi,2):: olist
  integer(int64),intent(in),dimension(Norb):: slist,invs
  integer(int64),intent(in),dimension(3,Nkall):: kmap,invk
  logical(1),intent(in):: sw_in,sw_out,sw_rescale
  integer(int64),intent(in):: sub_sigma
  real(real64),intent(in):: temp,eps,pp,rfill,mu_init
  real(real64),intent(in),dimension(Norb,Nk):: eig
  real(real64),intent(in),dimension(Nchi,Nchi):: Vmat
  real(real64),intent(out):: mu
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni,hamk
  complex(real64),intent(out),dimension(Nk,Nw,Norb,Norb):: sigmak

  integer(int32) scf_i
  integer(int32),dimension(Nchi,Nchi,2)::chi_map
  integer(int32),dimension(Nchi*(Nchi+1)/2,2)::irr_chi
  integer(int64),dimension(Nchi):: invschi
  real(real64)esterr,mu_OLD,eps_sgm,maxchi0_global
  real(real64),dimension(Norb,Norb):: sgnsig
  real(real64),dimension(Nchi,Nchi):: sgnsig2
  complex(real64),dimension(Nk,Nw,Norb,Norb):: Gk,sigmak0
  complex(real64),dimension(Nk,Nw,Nchi,Nchi):: chi
  ! DIIS
  integer(int32):: n_hist,i_hist
  complex(real64),allocatable:: xout_hist(:,:,:,:,:),res_hist(:,:,:,:,:)
  real(real64),allocatable:: B_diis(:,:),rhs_diis(:)
  integer(int32),allocatable:: ipiv_diis(:)
  ! bracket cache for renew_mu
  real(real64):: muS_cache,muL_cache
  logical:: bracket_valid

  eps_sgm=1.0d-10
  mu=mu_init
  if(sw_in)then
     print*,"load self"
     call io_sigma(.false.)
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
  call get_invschi()
  call get_sgnsig()
  call get_chi_map_soc(chi_map,irr_chi,olist,invs,Nchi,Norb)
  allocate(xout_hist(Nk,Nw,Norb,Norb,m_diis))
  allocate(res_hist (Nk,Nw,Norb,Norb,m_diis))
  allocate(B_diis(m_diis+1,m_diis+1))
  allocate(rhs_diis(m_diis+1))
  allocate(ipiv_diis(m_diis+1))
  n_hist=0
  i_hist=0
  bracket_valid=.false.
  muS_cache=0.0d0
  muL_cache=0.0d0
  iter_loop: do scf_i=1,scf_loop
     print'(A5,I5)','iter=',scf_i
     call get_chi0_conv_soc(chi,Gk,kmap,invk,invs,irr_chi,chi_map,olist,sgnsig,&
          sgnsig2,temp,Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi)
     call ckchi()
     if(sw_rescale .and. maxchi0_global>=1.0d0)then
        print'(A,F10.6,A)','[FLEX] Stoner factor=',maxchi0_global,'>= 1: rescaling chi0'
        chi(:,:,:,:)=chi(:,:,:,:)*(1.0d0-1.0d-4)/maxchi0_global
     end if
     call get_V_soc_flex(chi,Vmat,sgnsig2,Nk,Nw,Nchi)
     print'(A16,E12.4,A5,E12.4)','Re V_sigma: max:',maxval(dble(chi)),' min:',minval(dble(chi))
     call calc_sigma_soc(sigmak,Gk,chi,Vmat,kmap,invk,invs,olist,slist,sgnsig,sgnsig2,temp,Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz)
     if(sub_sigma>0)then
        sub_self:block
          integer(int32) iw,l,m
          complex(real64),dimension(Nk,Norb,Norb):: sub_sigmak
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
  deallocate(xout_hist,res_hist,B_diis,rhs_diis,ipiv_diis)
  if(sw_out)then
     call io_sigma(.true.)
  end if
  call renew_mu()
contains
  subroutine get_sgnsig()
    integer(int32) i,j
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

  subroutine get_invschi()
    integer(int32) l,m,i,j
    invschi(:)=0
    do l=1,Nchi
       i=invs(olist(l,1))
       j=invs(olist(l,2))
       do m=1,Nchi
          if(i==olist(m,1) .and. j==olist(m,2))then
             invschi(l)=m
             exit
          end if
       end do
    end do
  end subroutine get_invschi

  subroutine ckchi()
    integer(int32) i,l,m,n,info,chik,chikall,iorb
    real(real64) maxchi0,maxchi02
    real(real64),dimension(2*Nchi):: rwork
    complex(real64),dimension(Nchi*Nchi*4+1):: work
    complex(real64),dimension(Nchi):: eigc
    complex(real64),dimension(Nchi,Nchi):: chi0,tmp,tmp1

    maxchi02=-1.0d5
    do i=1,Nk
       !$omp parallel
       !$omp workshare
       chi0(:,:)=0.0d0
       !$omp end workshare
       !$omp do private(l,m,n)
       do l=1,Nchi
          do m=1,Nchi
             do n=1,Nchi
                chi0(m,l)=chi0(m,l)-chi(i,1,m,n)*Vmat(n,l)
             end do
          end do
       end do
       !$omp end do
       !$omp end parallel
       call zgeev('N','N',Nchi,chi0,Nchi,eigc,tmp1,Nchi,tmp,Nchi,work,Nchi*Nchi*4+1,rwork,info)
       if(info/=0)then; print*,'zgeev failed: info=',info; stop; end if
       maxchi0=maxval(dble(eigc))
       if(maxchi0>maxchi02)then
          chik=i
          maxchi02=maxchi0
       end if
    end do
    do i=1,Nkall
       if(invk(2,i)==0 .and. invk(1,i)==chik)then
          chikall=i
          exit
       end if
    end do
    print'(A7,3I4,F12.8)','SDW/CDW',kmap(:,chikall),maxchi02
    maxchi0_global=maxchi02
  end subroutine ckchi

  subroutine compare_sigma()
      integer(int32) i,j,l,m,kerr,iwerr,lerr,merr,ih,jh,idx_i,idx_j,info,n_cur
      real(real64) est
      real(real64) eps_reg
      complex(real64) cksigm
      complex(real64),dimension(Nk,Nw,Norb,Norb):: sigma_diis

    ! --- Store current output and residual in circular buffer ---
    i_hist=mod(i_hist,int(m_diis,int32))+1
    n_hist=min(n_hist+1,int(m_diis,int32))
    xout_hist(:,:,:,:,i_hist)=sigmak(:,:,:,:)
    res_hist(:,:,:,:,i_hist)=sigmak(:,:,:,:)-sigmak0(:,:,:,:)

    ! determine current DIIS dimension
    n_cur = n_hist

    B_diis(1:n_cur+1,1:n_cur+1)=0.0d0
    do ih=1,n_cur
       idx_i=modulo(i_hist-n_cur+ih-1,int(m_diis,int32))+1
       do jh=1,n_cur
          idx_j=modulo(i_hist-n_cur+jh-1,int(m_diis,int32))+1
          B_diis(ih,jh)=dble(sum(conjg(res_hist(:,:,:,:,idx_i))*res_hist(:,:,:,:,idx_j)))
       end do
       B_diis(ih,n_cur+1)=1.0d0
       B_diis(n_cur+1,ih)=1.0d0
    end do
    ! regularize diagonal slightly to avoid exact singularity
    eps_reg=1.0d-12
    if(n_cur>0) then
       do ih=1,n_cur
          B_diis(ih,ih)=B_diis(ih,ih)+eps_reg
       end do
    end if

    ! --- Right-hand side vector (Lagrange constraint: sum c_i = 1) ---
    rhs_diis(1:n_cur)=0.0d0
    rhs_diis(n_cur+1)=1.0d0

    ! --- Solve B*c = rhs via LAPACK dgesv ---
    call dgesv(n_cur+1,1,B_diis(1:n_cur+1,1:n_cur+1),n_cur+1, &
               ipiv_diis(1:n_cur+1),rhs_diis(1:n_cur+1),n_cur+1,info)
    ! If dgesv fails (singular matrix), fall back to linear mixing.
    ! Set weight 1 on the most recent entry (index n_cur in the compacted rhs).
    if(info/=0)then
       print*,'DIIS: dgesv failed (info=',info,'), fallback to most-recent entry'
       rhs_diis(1:n_cur)=0.0d0
       if(n_cur>=1) then
          rhs_diis(n_cur)=1.0d0
       end if
    end if

    ! --- DIIS extrapolation: sigma_diis = sum_i c_i * xout_hist_i ---
    sigma_diis(:,:,:,:)=0.0d0
    do ih=1,n_cur
       idx_i=modulo(i_hist-n_cur+ih-1,int(m_diis,int32))+1
       sigma_diis=sigma_diis+rhs_diis(ih)*xout_hist(:,:,:,:,idx_i)
    end do

    ! --- Convergence check and mixing ---
    ! n_cur=1 (first step or m_diis=1): fall back to linear mixing
    ! n_cur>=2 (DIIS active): use sigma_diis directly to preserve optimal extrapolation
    esterr=0.0d0
    est=100.0d0
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
                if(n_cur>=2)then
                   sigmak(i,j,m,l)=cksigm
                else
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
  end subroutine compare_sigma

  subroutine renew_mu()
    integer(int32) i_iter
    integer(int32),parameter:: itemax=100
    logical(int32) flag
    logical bracket_found
    real(real64) rnS,rnL,rnc,rnM,muc,mud,muL,muS,muM,eps,dmu

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
  end subroutine renew_mu
  
  subroutine get_rn(rn,rmu)
    use constant
    real(real64),intent(in):: rmu
    real(real64),intent(out):: rn

    integer(int32) l,i,j,n
    real(real64) tmp,deltagk
    complex(real64):: Gk0,iw

    tmp=sum(0.5d0*(1.0d0-tanh(0.5d0*(eig(:,:)-rmu)/temp)))
    call gen_green_inv(Gk,sigmak,hamk,rmu,temp,Nk,Nw,Norb)
    call getinv(Gk,Nk,Nw,Norb)
    deltagk=0.0d0
    !$omp parallel do reduction(+:deltagk) private(l,i,n,iw,Gk0)
    do j=1,Nw
       do l=1,Norb
          iw=cmplx(rmu,dble(2*(j-1)+1)*pi*temp,kind=real64)
          do i=1,Nk
             Gk0=0.0d0
             do n=1,Norb
                Gk0=Gk0+uni(l,n,i)*conjg(uni(l,n,i))/(iw-eig(n,i))
             end do
             deltagk=deltagk+dble(Gk(i,j,l,l)-Gk0)
          end do
       end do
    end do
    !$omp end parallel do
    rn=(tmp+2*temp*deltagk)/Nk
  end subroutine get_rn
  
  subroutine io_sigma(sw)
    logical(int32),intent(in):: sw !True: out, False: in
    integer(int32)i,j,l,m
    open(55,file='sigma.bin',form='unformatted')
    if(sw)then
       write(55)mu
       write(55)mu_OLD
    else
       read(55)mu
       read(55)mu_OLD
    end if
    do l=1,Norb
       do m=1,Norb
          do j=1,Nw
             do i=1,Nk
                if(sw)then
                   write(55)sigmak(i,j,m,l)
                else
                   read(55)sigmak(i,j,m,l)
                end if
             end do
          end do
       end do
    end do
    close(55)
  end subroutine io_sigma
end subroutine mkself_soc

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
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz
  integer(int64),intent(in),dimension(3,Nkall):: kmap,invk
  integer(int64),intent(in),dimension(Nchi,2):: olist
  integer(int64),intent(in),dimension(Norb):: slist,invs
  real(real64),intent(in):: temp
  real(real64),intent(in),dimension(Nchi,Nchi):: Vmat
  real(real64),dimension(Norb,Norb):: sgnsig
  real(real64),dimension(Nchi,Nchi):: sgnsig2
  complex(real64),intent(in),dimension(Nk,Nw,Nchi,Nchi):: Vsigma
  complex(real64),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(real64),intent(out),dimension(Nk,Nw,Norb,Norb):: sigmak

  integer(int32) i,j,k,n,l,m
  real(real64) weight
  real(real64),parameter:: eps=1.0d-9
  complex(real64),dimension(0:Nx-1,0:Ny-1,0:Nz-1,2*Nw):: tmpVsigma,tmp,tmpgk

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
              tmpVsigma(kmap(1,i),kmap(2,i),kmap(3,i),1)=Vsigma(invk(1,i),1,m,l) !Vsigma(k,iw)
           else if(invk(2,i)==1)then
              tmpVsigma(kmap(1,i),kmap(2,i),kmap(3,i),1)=Vsigma(invk(1,i),1,l,m) !Vsigma(-k,iw)
           end if
           tmpVsigma(kmap(1,i),kmap(2,i),kmap(3,i),Nw+1)=-Vmat(m,l)
        end do
        !$omp end do
        !$omp do private(i)
        do j=2,Nw
           do i=1,Nkall
              if(invk(2,i)==0)then
                 tmpVsigma(kmap(1,i),kmap(2,i),kmap(3,i),j)=Vsigma(invk(1,i),j,m,l) !Vsigma(k,iw)
                 tmpVsigma(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+2)=conjg(Vsigma(invk(1,i),j,l,m)) !Vsigma(k,-iw)ml=Vsigma^*lm(k,iw)
              else if(invk(2,i)==1)then
                 tmpVsigma(kmap(1,i),kmap(2,i),kmap(3,i),j)=Vsigma(invk(1,i),j,l,m) !Vsigma(-k,iw)=Vsigma(-k,iw)^T
                 tmpVsigma(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+2)=conjg(Vsigma(invk(1,i),j,m,l)) !Vsigma(-k,-iw)=Vsigma^*(k,iw)
              end if
           end do
        end do
        !$omp end do
        !$omp do private(i)
        do j=1,Nw
           do i=1,Nkall
              if(invk(2,i)==0)then
                 tmpgk(kmap(1,i),kmap(2,i),kmap(3,i),j)=Gk(invk(1,i),j,olist(m,2),olist(l,2)) !G42(iw)
                 tmpgk(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=conjg(Gk(invk(1,i),j,olist(l,2),olist(m,2))) !G42(k,-iw)=G^*24(k,iw)
              else if(invk(2,i)==1)then
                 tmpgk(kmap(1,i),kmap(2,i),kmap(3,i),j)=sgnsig(olist(l,2),olist(m,2))*Gk(invk(1,i),j,invs(olist(l,2)),invs(olist(m,2))) !G42(-k,iw)=G24(k,iw)
                 tmpgk(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=sgnsig(olist(l,2),olist(m,2))*conjg(Gk(invk(1,i),j,invs(olist(m,2)),invs(olist(l,2)))) !G42(-k,-iw)=G^*42(k,iw)
              end if
           end do
        end do
        !$omp end do
        !$omp end parallel
        call FFT(tmpVsigma,tmp,Nx,Ny,Nz,2*Nw,.true.)
        call FFT(tmpgk,tmp,Nx,Ny,Nz,2*Nw,.true.)
        !$omp parallel
        !$omp do private(i,j,k,n)
        do n=1,2*Nw
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
      tmpgk=cmplx(0.0d0,0.0d0,kind=real64)
        !$omp end workshare
        !$omp end parallel
        call FFT(tmp,tmpgk,Nx,Ny,Nz,2*Nw,.false.)
        !$omp parallel do private(i,j)
        do j=1,Nw
           do i=1,Nkall
              if(invk(2,i)==0)then
                 sigmak(invk(1,i),j,olist(m,1),olist(l,1))=sigmak(invk(1,i),j,olist(m,1),olist(l,1))&
                      +tmp(kmap(1,i),kmap(2,i),kmap(3,i),j)
              end if
           end do
        end do
        !$omp end parallel do
     end do
  end do

  do l=1,Norb
     do m=1,Norb
        !$omp parallel do private(i,j)
        do j=1,Nw
           do i=1,Nk
              sigmak(i,j,m,l)=sigmak(i,j,m,l)*weight
              if(abs(dble(sigmak(i,j,m,l)))<eps) sigmak(i,j,m,l)=cmplx(0.0d0,imag(sigmak(i,j,m,l)),kind=real64)
              if(abs(imag(sigmak(i,j,m,l)))<eps) sigmak(i,j,m,l)=cmplx(dble(sigmak(i,j,m,l)),0.0d0,kind=real64)
           end do
        end do
        !$omp end parallel do
     end do
  end do
end subroutine calc_sigma_soc
