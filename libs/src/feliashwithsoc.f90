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
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,itemax,gap_sym
  integer(c_int64_t),intent(in):: arnoldi_m  !! Krylov subspace dimension (0=power method, >0=Arnoldi)
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
  complex(c_double),intent(out),dimension(Nkall,Nw,Norb,Norb):: delta
  complex(c_double),intent(inout),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(c_int32_t) i_iter,i_eig,count,i
  integer(c_int32_t),parameter:: eig_max=2
  real(c_double) norm,inorm,norm2,weight,lambda_rq,vec_err
  complex(c_double),dimension(Nkall,Nw,Norb,Norb):: newdelta,fk

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
    complex(c_double),intent(in),dimension(Nkall,Nw,Norb,Norb):: func
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
    complex(c_double),intent(in),dimension(Nkall,Nw,Norb,Norb):: del,newdel
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
    complex(c_double),intent(in),dimension(Nkall,Nw,Norb,Norb):: newdel,del
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
    integer(c_int64_t),intent(in):: m_arnoldi
    integer(c_int32_t) :: m_dim,ip,j,ii,m_act,idx,lwork,info
    real(c_double) :: beta,shift,lambda1
    complex(c_double) :: hij
    complex(c_double),allocatable :: V(:,:,:,:,:),delta1(:,:,:,:)
    complex(c_double),allocatable :: H_mat(:,:),A_mat(:,:),eigenvals(:)
    complex(c_double),allocatable :: VL_dum(:,:),VR(:,:),zwork(:)
    real(c_double),allocatable :: rwork(:)

    m_dim=int(m_arnoldi,c_int32_t)
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
       H_mat(j+2,j+1)=cmplx(beta,0.0d0,c_double)
       if(beta<1.0d-14)then; m_act=j+1; exit; end if
       if(j<m_dim-1)then
          !$omp parallel do
          do ip=1,int(Nkall,c_int32_t)
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
    do ip=1,int(Nkall,c_int32_t)
       delta1(ip,:,:,:)=(0.0d0,0.0d0)
    end do
    !$omp end parallel do
    do ii=1,m_act
       !$omp parallel do
       do ip=1,int(Nkall,c_int32_t)
          delta1(ip,:,:,:)=delta1(ip,:,:,:)+VR(ii,idx)*V(ip,:,:,:,ii-1)
       end do
       !$omp end parallel do
    end do
    call get_norm(beta,delta1)
    !$omp parallel do
    do ip=1,int(Nkall,c_int32_t)
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
    do ip=1,int(Nkall,c_int32_t)
       V(ip,:,:,:,0)=V(ip,:,:,:,0)-hij*delta1(ip,:,:,:)
    end do
    !$omp end parallel do
    call get_norm(beta,V(:,:,:,:,0))
    !$omp parallel do
    do ip=1,int(Nkall,c_int32_t)
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
       H_mat(j+2,j+1)=cmplx(beta,0.0d0,c_double)
       if(beta<1.0d-14)then; m_act=j+1; exit; end if
       if(j<m_dim-1)then
          !$omp parallel do
          do ip=1,int(Nkall,c_int32_t)
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
    complex(c_double),intent(in),dimension(Nkall,Nw,Norb,Norb):: v_in
    complex(c_double),intent(out),dimension(Nkall,Nw,Norb,Norb):: w_out
    call mkfk_trs_soc(fk,Gk,v_in,sgnsig,slist,invk,invs,Nkall,Nk,Nw,Norb,gap_sym)
    call mkdelta_soc(w_out,fk,chi,Vmat,sgnsig,sgnsig2,kmap,invk,invs,invschi,olist,slist,Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,gap_sym)
    !$omp parallel workshare
      w_out(:,:,:,:)=-w_out(:,:,:,:)*weight
    !$omp end parallel workshare
  end subroutine apply_op

  subroutine get_inner(h,u,v)
    !> complex inner product h = 2*<u,v> summed over full BZ
    complex(c_double),intent(in),dimension(Nkall,Nw,Norb,Norb):: u,v
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
                tmp_r=tmp_r+dble(conjg(u(i,j,m,l))*v(i,j,m,l))
                tmp_i=tmp_i+aimag(conjg(u(i,j,m,l))*v(i,j,m,l))
             end do
          end do
          !$omp end do
       end do
    end do
    !$omp end parallel
    h=cmplx(2.0d0*tmp_r,2.0d0*tmp_i,c_double)
  end subroutine get_inner
end subroutine lin_eliash_soc

subroutine mkfk_trs_soc(fk,Gk,delta,sgnsig,slist,invk,invs,Nkall,Nk,Nw,Norb,gap_sym) bind(C,name="mkfk_trs_soc_")
  !
  !>calculate linearized anomalous green function Fk with TRS
  !>if we considder TRS,F_ab(k)=G_ac(k)Delta_cd(k)Gbd(-k)
  !> =G_ac(k)Delta_cd(k)G^*_dbss'(k)ss'
  !
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nkall,Nk,Nw,Norb,gap_sym
  integer(c_int64_t),intent(in),dimension(3,Nkall):: invk
  integer(c_int64_t),intent(in),dimension(Norb):: invs,slist
  real(c_double),intent(in),dimension(Norb,Norb):: sgnsig
  complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(c_double),intent(in),dimension(Nkall,Nw,Norb,Norb):: delta
  complex(c_double),intent(out),dimension(Nkall,Nw,Norb,Norb):: fk
   
  integer(c_int32_t) i,j,l,m,n
  complex(c_double),dimension(Norb,Norb):: gmat1,gmat2,cmat1

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
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,gap_sym
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),dimension(Norb):: slist,invs
  integer(c_int64_t),intent(in),dimension(Nchi):: invschi
  real(c_double),intent(in),dimension(Norb,Norb):: sgnsig
  real(c_double),intent(in),dimension(Nchi,Nchi):: Vmat,sgnsig2
  complex(c_double),intent(in),dimension(Nk,Nw,Nchi,Nchi):: Vdelta
  complex(c_double),intent(in),dimension(Nkall,Nw,Norb,Norb):: delta
  complex(c_double),intent(out),dimension(Nkall,Nw,Norb,Norb):: newdelta

  integer(c_int32_t) i,j,k,n,l,m
  integer(c_int64_t),parameter::  par_mix=-10
  complex(c_double),dimension(0:Nx-1,0:Ny-1,0:Nz-1,2*Nw):: tmpVdelta,tmpfk,tmp
  
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
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  use constant
  implicit none
  integer(c_int64_t),intent(in):: Nw,Norb,Nkall,Nk,gap_sym
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
  integer(c_int64_t),intent(in),dimension(Norb):: slist,invs
  real(c_double),intent(in),dimension(Nk,Norb):: init_delta
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(c_double),intent(out),dimension(Nkall,Nw,Norb,Norb):: delta

   integer(c_int32_t) i,j,k,l,m,n
  real(c_double) norm
   complex(c_double),dimension(Norb,Norb):: matL,matR,tmpm

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
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nw,Norb,Nk,Nkall
  integer(c_int64_t),intent(in),dimension(Norb):: slist,invs
  integer(c_int64_t),intent(in),dimension(3,Nkall):: invk
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(c_double),intent(in),dimension(Nkall,Nw,Norb,Norb):: delta
  complex(c_double),intent(out),dimension(Nkall,Norb,Norb):: deltab

   integer(c_int32_t) i,k,l,n
   complex(c_double),dimension(Norb,Norb):: matL,matR,tmpm

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
