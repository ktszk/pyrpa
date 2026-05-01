subroutine mkself_soc(sigmak,mu,Vmat,kmap,invk,invs,olist,slist,hamk,eig,uni,mu_init,&
     rfill,temp,scf_loop,pp,eps,Nkall,Nk,Nw,Norb,Nchi,Nx,Ny,Nz,sub_sigma,sw_out,&
     sw_in,m_diis,sw_rescale) bind(C)
   use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t,c_bool
  implicit none
  integer(c_int64_t),intent(in):: Nw,Norb,Nchi,Nkall,Nk,Nx,Ny,Nz,scf_loop,m_diis
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),dimension(Norb):: slist,invs
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
   logical(c_bool),intent(in):: sw_in,sw_out,sw_rescale
  integer(c_int64_t),intent(in):: sub_sigma
  real(c_double),intent(in):: temp,eps,pp,rfill,mu_init
  real(c_double),intent(in),dimension(Norb,Nk):: eig
  real(c_double),intent(in),dimension(Nchi,Nchi):: Vmat
  real(c_double),intent(out):: mu
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni,hamk
  complex(c_double),intent(out),dimension(Nk,Nw,Norb,Norb):: sigmak

  integer(c_int32_t) scf_i
  integer(c_int32_t),dimension(Nchi,Nchi,2)::chi_map
  integer(c_int32_t),dimension(Nchi*(Nchi+1)/2,2)::irr_chi
  integer(c_int64_t),dimension(Nchi):: invschi
  real(c_double)esterr,mu_OLD,eps_sgm,maxchi0_global
  real(c_double),dimension(Norb,Norb):: sgnsig
  real(c_double),dimension(Nchi,Nchi):: sgnsig2
  complex(c_double),dimension(Nk,Nw,Norb,Norb):: Gk,sigmak0
  complex(c_double),dimension(Nk,Nw,Nchi,Nchi):: chi
  ! DIIS
  integer(c_int32_t):: n_hist,i_hist
  complex(c_double),allocatable:: xout_hist(:,:,:,:,:),res_hist(:,:,:,:,:)
  real(c_double),allocatable:: B_diis(:,:),rhs_diis(:)
  integer(c_int32_t),allocatable:: ipiv_diis(:)
  ! bracket cache for renew_mu
  real(c_double):: muS_cache,muL_cache
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
          integer(c_int32_t) iw,l,m
          complex(c_double),dimension(Nk,Norb,Norb):: sub_sigmak
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

  subroutine get_invschi()
    integer(c_int32_t) l,m,i,j
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
    integer(c_int32_t) i,l,m,n,info,chik,chikall,iorb
    real(c_double) maxchi0,maxchi02
    real(c_double),dimension(2*Nchi):: rwork
    complex(c_double),dimension(Nchi*Nchi*4+1):: work
    complex(c_double),dimension(Nchi):: eigc
    complex(c_double),dimension(Nchi,Nchi):: chi0,tmp,tmp1

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
      integer(c_int32_t) i,j,l,m,kerr,iwerr,lerr,merr,ih,jh,idx_i,idx_j,info,n_cur
      real(c_double) est
      real(c_double) eps_reg
      complex(c_double) cksigm
      complex(c_double),dimension(Nk,Nw,Norb,Norb):: sigma_diis

    ! --- Store current output and residual in circular buffer ---
    i_hist=mod(i_hist,int(m_diis,c_int32_t))+1
    n_hist=min(n_hist+1,int(m_diis,c_int32_t))
    xout_hist(:,:,:,:,i_hist)=sigmak(:,:,:,:)
    res_hist(:,:,:,:,i_hist)=sigmak(:,:,:,:)-sigmak0(:,:,:,:)

    ! determine current DIIS dimension
    n_cur = n_hist

    B_diis(1:n_cur+1,1:n_cur+1)=0.0d0
    do ih=1,n_cur
       idx_i=modulo(i_hist-n_cur+ih-1,int(m_diis,c_int32_t))+1
       do jh=1,n_cur
          idx_j=modulo(i_hist-n_cur+jh-1,int(m_diis,c_int32_t))+1
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
       idx_i=modulo(i_hist-n_cur+ih-1,int(m_diis,c_int32_t))+1
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
  end subroutine renew_mu
  
  subroutine get_rn(rn,rmu)
    use constant
    use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
    real(c_double),intent(in):: rmu
    real(c_double),intent(out):: rn

    integer(c_int32_t) l,i,j,n
    real(c_double) tmp,deltagk
    complex(c_double):: Gk0,iw

    tmp=sum(0.5d0*(1.0d0-tanh(0.5d0*(eig(:,:)-rmu)/temp)))
    call gen_green_inv(Gk,sigmak,hamk,rmu,temp,Nk,Nw,Norb)
    call getinv(Gk,Nk,Nw,Norb)
    deltagk=0.0d0
    !$omp parallel do reduction(+:deltagk) private(l,i,n,iw,Gk0)
    do j=1,Nw
       do l=1,Norb
          iw=cmplx(rmu,dble(2*(j-1)+1)*pi*temp,kind=c_double)
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
    logical,intent(in):: sw !True: out, False: in
    integer(c_int32_t)i,j,l,m
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
  real(c_double) weight
  real(c_double),parameter:: eps=1.0d-9
  complex(c_double),dimension(0:Nx-1,0:Ny-1,0:Nz-1,2*Nw):: tmpVsigma,tmp,tmpgk

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
      tmpgk=cmplx(0.0d0,0.0d0,kind=c_double)
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
              if(abs(dble(sigmak(i,j,m,l)))<eps) sigmak(i,j,m,l)=cmplx(0.0d0,imag(sigmak(i,j,m,l)),kind=c_double)
              if(abs(imag(sigmak(i,j,m,l)))<eps) sigmak(i,j,m,l)=cmplx(dble(sigmak(i,j,m,l)),0.0d0,kind=c_double)
           end do
        end do
        !$omp end parallel do
     end do
  end do
end subroutine calc_sigma_soc
