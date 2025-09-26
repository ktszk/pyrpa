subroutine lin_eliash_soc(delta,chi,Gk,uni,init_delta,Vmat,sgnsig,sgnsig2,olist,slist,kmap,invk,invs,invschi,&
     temp,eps,Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,itemax,gap_sym) bind(C)
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
  integer(int64),intent(in),dimension(Nchi,2):: olist
  integer(int64),intent(in),dimension(Norb):: slist,invs
  integer(int64),intent(in),dimension(3,Nkall):: kmap,invk
  integer(int64),intent(in),dimension(Nchi):: invschi
  real(real64),intent(in):: temp,eps
  real(real64),intent(in),dimension(Nchi,Nchi):: Vmat
  real(real64),intent(in),dimension(Norb,Norb):: sgnsig
  real(real64),intent(in),dimension(Nchi,Nchi):: sgnsig2
  real(real64),intent(in),dimension(Nk,Norb):: init_delta
  complex(real64),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(real64),intent(out),dimension(Nkall,Nw,Norb,Norb):: delta
  complex(real64),intent(inout),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(int32) i_iter,i_eig,count,i
  integer(int32),parameter:: eig_max=2
  real(real64) norm,normb,inorm,norm2,weight
  complex(real64),dimension(Nkall,Nw,Norb,Norb):: newdelta,fk

  weight=temp/dble(Nkall)
  norm2=0.0d0
  normb=0.0d0
  call get_V_soc_flex(chi,Vmat,sgnsig2,Nk,Nw,Nchi)
  print'(A15,2E16.8)','V_delta max is ',maxval(dble(chi)),maxval(aimag(chi))
  print'(A15,2E16.8)','V_delta min is ',minval(dble(chi)),minval(aimag(chi))
  eigenval_loop:do i_eig=1,eig_max !solve eig_val using power method, 1st eig is usually large negative value
     call get_initial_delta_soc(delta,init_delta,uni,kmap,slist,invk,invs,Nkall,Nk,Nw,Norb,gap_sym)
     count=0 !count too small eigenvalue
     iter_loop:do i_iter=1,itemax !iteration
        call mkfk_trs_soc(fk,Gk,delta,sgnsig,invk,invs,Nkall,Nk,Nw,Norb)
        !something worng mkdelta_soc
        call mkdelta_soc(newdelta,fk,chi,Vmat,sgnsig,sgnsig2,kmap,invk,invs,invschi,olist,slist,Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,gap_sym)
        !$omp parallel workshare
        newdelta(:,:,:,:)=newdelta(:,:,:,:)*weight-delta(:,:,:,:)*norm2
        !$omp end parallel workshare
        call get_norm(norm,newdelta)
        inorm=1.0d0/norm
        if(abs(norm-norm2)>=1.0d2 .or. abs(norm-norm2)<1.0d-6)then
           print'(I3,A13,2E16.8)',i_iter,' lambda_elsh=',norm-norm2
        else
           print'(I3,A13,2F12.8)',i_iter,' lambda_elsh=',norm-norm2
        end if
        if(abs(norm-normb)*inorm<eps)then
           if(norm-norm2>1.0d-1)then !do not finish until eig>0.1
              exit
           else if(.true.)then !consider small eig
              if((norm-norm2)>0.0 .and. abs((norm-normb)*inorm)<eps*1.0d-2)then
                 count=count+1
              end if
              if(count>30)exit !if eigenvalue <0.1  until 30 count exit
           end if
        end if
        normb=norm
        !$omp parallel workshare
        delta(:,:,:,:)=newdelta(:,:,:,:)*inorm
        !$omp end parallel workshare
     end do iter_loop
     if(i_eig==2)then
        print*,'eliash=',norm-norm2
     end if
     norm2=norm
  end do eigenval_loop
  call get_norm(norm,newdelta)
  !$omp parallel workshare
  delta(:,:,:,:)=newdelta(:,:,:,:)*inorm
  !$omp end parallel workshare
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
  
  integer(int32) i,j,l,m,n,info,l2,m2,iorb
  integer(int32),dimension(Nchi):: ipiv
  complex(real64),dimension(Nchi,Nchi):: cmat1,cmat2,cmat3
  complex(real64),dimension(2*Nchi):: work

  wloop:do j=1,Nw
     qloop:do i=1,Nk
        !$omp parallel
        !$omp workshare
        cmat1(:,:)=0.0d0
        !$omp end workshare
        !$omp do private(m,n)
        do l=1,Nchi
           do m=1,Nchi
              do n=1,Nchi
                 cmat1(m,l)=cmat1(m,l)+chi(i,j,m,n)*Vmat(n,l) !-chi0V
              end do
           end do
        end do
        !$omp end do
        !$omp workshare
        cmat2(:,:)=cmat1(:,:) !chi0V
        !$omp end workshare
        !$omp do 
        do l=1,Nchi
           cmat1(l,l)=cmat1(l,l)+1.0d0 !I-chi0V
        end do
        !$omp end do
        !$omp end parallel
        call zgetrf(Nchi,Nchi,cmat1,Nchi,ipiv,info)
        call zgetri(Nchi,cmat1,Nchi,ipiv,work,2*Nchi,info)
        !$omp parallel
        !$omp workshare
        cmat3(:,:)=0.0d0
        !$omp end workshare
        !$omp do private(l,m,n)
        do l=1,Nchi
           do m=1,Nchi
              do n=1,Nchi
                 cmat3(m,l)=cmat3(m,l)+cmat1(m,n)*cmat2(n,l) !(1-chi0V)^-1chi0V
              end do
           end do
        end do
        !$omp end do
        !$omp workshare
        cmat2(:,:)=-Vmat(:,:)
        !$omp end workshare
        !$omp do private(l,m,n)
        do l=1,Nchi
           do m=1,Nchi
              do n=1,Nchi
                 cmat2(m,l)=cmat2(m,l)+Vmat(m,n)*cmat3(n,l)
              end do
           end do
        end do
        !$omp end do
        !$omp end parallel
        chi(i,j,:,:)=cmat2(:,:)
     end do qloop
  end do wloop
end subroutine get_V_soc_flex

subroutine mkfk_trs_soc(fk,Gk,delta,sgnsig,invk,invs,Nkall,Nk,Nw,Norb)
  !
  !>calculate linearized anomalous green function Fk with TRS
  !>if we considder TRS,F_ab(k)=G_ac(k)Delta_cd(k)Gbd(-k)
  !> =G_ac(k)Delta_cd(k)G^*_dbss'(k)ss'
  !
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nkall,Nk,Nw,Norb
  integer(int64),intent(in),dimension(3,Nkall):: invk
  integer(int64),intent(in),dimension(Norb):: invs
  real(real64),intent(in),dimension(Norb,Norb):: sgnsig
  complex(real64),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(real64),intent(in),dimension(Nkall,Nw,Norb,Norb):: delta
  complex(real64),intent(out),dimension(Nkall,Nw,Norb,Norb):: fk
   
  integer(int32) i,j,l,m,n
  complex(real64),dimension(Nkall,Nw,Norb,Norb):: cmat1
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
              do i=1,Nkall
                 if(invk(2,i)==0)then !cmat1(k,iw)=G(k,iw)Delta(k,iw)
                    cmat1(i,j,m,l)=cmat1(i,j,m,l)+Gk(invk(1,i),j,m,n)*delta(i,j,n,l)
                 else if(invk(2,i)==1)then !G(k,iw)=G^T(-k,iw)
                    cmat1(i,j,m,l)=cmat1(i,j,m,l)+sgnsig(n,m)*Gk(invk(1,i),j,invs(n),invs(m))*delta(i,j,n,l)
                 end if
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
              do i=1,Nkall
                 if(invk(2,i)==0)then
                    fk(i,j,m,l)=fk(i,j,m,l)-cmat1(i,j,m,n)*conjg(Gk(invk(1,invk(3,i)),j,l,n))
                 else if(invk(2,i)==1)then
                    fk(i,j,m,l)=fk(i,j,m,l)-cmat1(i,j,m,n)*sgnsig(n,l)*conjg(Gk(invk(1,invk(3,i)),j,invs(n),invs(l)))
                 end if
              end do
           end do
           !$omp end do
        end do
        if(l/=m)then !if TRS and no SOC delta^+(k,iw)=delta(k,iw)
           !$omp do private(i,j)
           do j=1,Nw
              do i=1,Nkall
                 fk(i,j,l,m)=conjg(fk(i,j,m,l))
              end do
           end do
           !$omp end do
        else
           !$omp do private(i,j)
           do j=1,Nw
              do i=1,Nkall
                 fk(i,j,l,l)=dble(fk(i,j,l,l))
              end do
           end do
           !$omp end do
        end if
     end do
  end do
  !$omp end parallel
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
        if((gap_sym<0) .or. (gap_sym>0 .and. slist(olist(l,1))*slist(olist(m,2))<0))then !triplet gap need all D, singlet gap only Dud(du)
           if((gap_sym==par_mix .and. slist(olist(l,1))*slist(olist(m,2))<0) .or. slist(olist(l,1))==1)then !calc only Duu, Dud
              if(slist(olist(m,2))==1 .or. olist(l,1)>=invs(olist(m,2)))then !Duu all, Dud upper triangle
                 !$omp parallel
                 !$omp do private(i)
                 do i=1,Nkall
                    if(invk(2,i)==0)then !k,0
                       tmpVdelta(kmap(1,i),kmap(2,i),kmap(3,i),1)=Vdelta(invk(1,i),1,m,l) !j=1 corresponds to w_n=0
                    else if(invk(2,i)==1)then !-k,0
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
                       tmpfk(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=-delta(invk(3,i),j,olist(m,1),olist(l,2)) !F(k,-w)=F(-k,w)
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
              if(abs(dble(chi(invk(1,i),j,m,l)))<eps) chi(invk(1,i),j,m,l)=cmplx(0.0d0,imag(chi(invk(1,i),j,m,l)))
              if(abs(imag(chi(invk(1,i),j,m,l)))<eps) chi(invk(1,i),j,m,l)=cmplx(dble(chi(invk(1,i),j,m,l)),0.0d0)
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
  real(real64),intent(in),dimension(Nkall,Norb):: init_delta
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(real64),intent(out),dimension(Nkall,Nw,Norb,Norb):: delta

  integer(int32) i,j,l,m,n
  real(real64) norm

  !$omp parallel
  !$omp workshare
  delta(:,:,:,:)=0.0d0
  !$omp end workshare
  do l=1,Norb !Delta_b(k)=U(k)Delta_o(k)U^T(-k)=-sU(k)Delta_o(k)U*^(nl)_(-s)(k)
     do m=1,Norb
        if(gap_sym<0 .or. (gap_sym>0 .and. slist(m)*slist(l)<0))then !singlet gap only ud,du spin
           !$omp do private(i,n)
           do i=1,Nkall
              if(invk(2,i)==0)then
                 do n=1,Norb
                    delta(i,1,m,l)=delta(i,1,m,l)-uni(m,n,invk(1,i))*init_delta(i,n)*slist(l)*conjg(uni(invs(l),n,invk(1,i)))
                 end do
              else if(invk(2,i)==1)then
                 do n=1,Norb
                    delta(i,1,m,l)=delta(i,1,m,l)+slist(m)*conjg(uni(invs(m),n,invk(1,i)))*init_delta(i,n)*uni(l,n,invk(1,i))
                 end do
              end if
           end do
           !$omp end do
           !$omp do private(i,j)
           do j=2,Nw
              !$omp simd
              do i=1,Nkall
                 delta(i,j,m,l)=delta(i,1,m,l)
              end do
              !$omp end simd
           end do
           !$omp end do
        end if
     end do
  end do
  !$omp end parallel

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

  integer(int32) i,l,m,n,info
  integer(int32),dimension(Nchi):: ipiv
  complex(real64),dimension(Nchi,Nchi):: cmat1,cmat2
  complex(real64),dimension(2*Nchi):: work

  qloop:do i=1,Nk
     !$omp parallel
     !$omp workshare
     cmat1(:,:)=0.0d0
     !$omp end workshare
     !$omp do private(m,n)
     do l=1,Nchi
        do m=1,Nchi
           do n=1,Nchi
              cmat1(m,l)=cmat1(m,l)+chi(i,1,m,n)*Vmat(n,l) !-chi0V
           end do
        end do
     end do
     !$omp end do
     !$omp do
     do l=1,Nchi
        cmat1(l,l)=cmat1(l,l)+1.0d0 !I-chi0V
     end do
     !$omp end do
     !$omp end parallel
     call zgetrf(Nchi,Nchi,cmat1,Nchi,ipiv,info)
     call zgetri(Nchi,cmat1,Nchi,ipiv,work,2*Nchi,info)
     !$omp parallel
     !$omp workshare
     cmat2(:,:)=0.0d0
     !$omp end workshare
     !$omp do private(l,m,n)
     do l=1,Nchi
        do m=1,Nchi
           do n=1,Nchi
              cmat2(m,l)=cmat2(m,l)+cmat1(m,n)*chi(i,1,n,l) !(1-chi0V)^-1chi0V
           end do
        end do
     end do
     !$omp end do
     !$omp end parallel
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

  integer(int32) i,l,m,n
  complex(real64),dimension(Nkall,Norb,Norb):: tmp

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
                 tmp(i,m,l)=tmp(i,m,l)+conjg(uni(n,m,invk(1,i)))*delta(i,1,n,l)
              else if(invk(2,i)==1)then
                 tmp(i,m,l)=tmp(i,m,l)+slist(m)*uni(invs(n),m,invk(1,i))*delta(i,1,n,l)
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
                 deltab(i,m,l)=deltab(i,m,l)+tmp(i,m,n)*slist(n)*uni(invs(n),l,invk(1,i))
              else if(invk(2,i)==1)then
                 deltab(i,m,l)=deltab(i,m,l)+tmp(i,m,n)*conjg(uni(n,l,invk(1,i)))
              end if
           end do
           !$omp end do
        end do
     end do
  end do
  !$omp end parallel
end subroutine conv_delta_orb_to_band_soc

subroutine mkself_soc(sigmak,mu,Vmat,kmap,invk,invs,olist,slist,hamk,eig,uni,mu_init,&
     rfill,temp,scf_loop,pp,eps,Nkall,Nk,Nw,Norb,Nchi,Nx,Ny,Nz,sw_sub_sigma,sw_out,&
     sw_in) bind(C)
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nw,Norb,Nchi,Nkall,Nk,Nx,Ny,Nz,scf_loop
  integer(int64),intent(in),dimension(Nchi,2):: olist
  integer(int64),intent(in),dimension(Norb):: slist,invs
  integer(int64),intent(in),dimension(3,Nkall):: kmap,invk
  logical(1),intent(in):: sw_in,sw_out,sw_sub_sigma
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
  real(real64)esterr,mu_OLD,eps_sgm
  real(real64),dimension(Norb,Norb):: sgnsig
  real(real64),dimension(Nchi,Nchi):: sgnsig2
  complex(real64),dimension(Nk,Nw,Norb,Norb):: Gk,sigmak0
  complex(real64),dimension(Nk,Nw,Nchi,Nchi):: chi

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
  iter_loop: do scf_i=1,scf_loop
     print'(A5,I5)','iter=',scf_i
     call get_chi0_conv_soc(chi,Gk,kmap,invk,invs,irr_chi,chi_map,olist,sgnsig,&
          sgnsig2,temp,Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi)
     call ckchi()
     call get_V_soc_flex(chi,Vmat,sgnsig2,Nk,Nw,Nchi)
     print'(A16,E12.4,A5,E12.4)','Re V_sigma: max:',maxval(dble(chi)),' min:',minval(dble(chi))
     call calc_sigma_soc(sigmak,Gk,chi,Vmat,kmap,invk,olist,temp,Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz)
     if(sw_sub_sigma)then
        sub_self:block
          integer(int32) iw,l,m
          complex(real64),dimension(Nk,Norb,Norb):: sub_sigmak
          !$omp parallel
          do l=1,Norb
             do m=l,Norb
                !$omp workshare
                sub_sigmak(:,l,m)=(sigmak(:,1,l,m)+conjg(sigmak(:,1,m,l)))*0.5d0
                !$omp end workshare
                if(l.ne.m)sub_sigmak(:,m,l)=conjg(sub_sigmak(:,l,m))
             end do
          end do
          !$omp do
          do iw=1,Nw
             sigmak(:,iw,:,:)=sigmak(:,iw,:,:)-sub_sigmak(:,:,:)
          end do
          !$omp end do
          !$omp end parallel
        end block sub_self
     end if
     call compair_sigma()
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

  subroutine compair_sigma()
    integer(int32) i,j,l,m, kerr,iwerr,lerr,merr
    real(real64) est
    complex(real64) cksigm

    esterr=0.0d0
    est=100
    do l=1,Norb
       do m=1,Norb
          do j=1,Nw
             do i=1,Nk
                cksigm=sigmak(i,j,m,l)
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
                sigmak(i,j,m,l)=pp*cksigm+(1.0d0-pp)*sigmak0(i,j,m,l)
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
  end subroutine compair_sigma

  subroutine renew_mu()
    integer(int32) i_iter
    integer(int32),parameter:: itemax=100
    logical(int32) flag
    real(real64) rnS,rnL,rnc,rnM,muc,mud,muL,muS,muM,eps,dmu

    if(esterr>1.0d-2)then
       eps= 1.0d-8
    else
       eps= esterr*1.0d-1
    end if
    dmu= abs(mu-mu_OLD)*2.0d0
    if (dmu<eps*4.0d0) dmu= eps*4.0d0
    mu_OLD=mu
    muL= mu+dmu
    muS= mu-dmu
    upper_lim: do i_iter=1,itemax
       mu=muL
       rnL=00d0
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
    do l=1,Norb
       do j=1,Nw
          iw=cmplx(mu,dble(2*(j-1)+1)*pi*temp)          
          do i=1,Nk
             Gk0=0.0d0
             do n=1,Norb
                Gk0=Gk0+uni(l,n,i)*conjg(uni(l,n,i))/(iw-eig(n,i))
             end do
             deltagk=deltagk+dble(Gk(i,j,l,l)-Gk0)
          end do
       end do
    end do
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

subroutine calc_sigma_soc(sigmak,Gk,Vsigma,Vmat,kmap,invk,olist,temp,Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz)
  !> This function obtain new gap function Delta=sum VF
  !!@param sigmak,out: self energy
  !!@param      Gk,in: green function
  !!@param  Vsigma,in: irreducible susceptibility and interaction
  !!@param    Smat,in: S-matrix
  !!@param    Cmat,in: C-matrix
  !!@param    kmap,in: property of k-point
  !!@param    invk,in: list of reverse k-points
  !!@param   olist,in: property of chi footnote
  !!@param    temp,in: Tempearature
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
  real(real64),intent(in):: temp
  real(real64),intent(in),dimension(Nchi,Nchi):: Vmat
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
           tmpVsigma(kmap(1,i),kmap(2,i),kmap(3,i),Nw+1)=Vmat(m,l)
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
                 tmpgk(kmap(1,i),kmap(2,i),kmap(3,i),j)=Gk(invk(1,i),j,olist(l,2),olist(m,2)) !G42(-k,iw)=G24(k,iw)
                 tmpgk(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=conjg(Gk(invk(1,i),j,olist(m,2),olist(l,2))) !G42(-k,-iw)=G^*42(k,iw)
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
        tmpgk=cmplx(0.0d0,0.0d0)
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
              if(abs(dble(sigmak(i,j,m,l)))<eps) sigmak(i,j,m,l)=cmplx(0.0d0,imag(sigmak(i,j,m,l)))
              if(abs(imag(sigmak(i,j,m,l)))<eps) sigmak(i,j,m,l)=cmplx(dble(sigmak(i,j,m,l)),0.0d0)
           end do
        end do
        !$omp end parallel do
     end do
  end do
end subroutine calc_sigma_soc
