subroutine lin_eliash_soc(delta,Gk,uni,init_delta,Vmat,olist,slist,kmap,invk,invs,temp,eps,&
     Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,itemax,gap_sym) bind(C)
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
  real(real64),intent(in):: temp,eps
  real(real64),intent(in),dimension(Nchi,Nchi):: Vmat
  real(real64),intent(in),dimension(Nk,Norb):: init_delta
  complex(real64),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(real64),intent(out),dimension(Nkall,Nw,Norb,Norb):: delta

  integer(int32) i_iter,i_eig,count,i
  integer(int32),parameter:: eig_max=2
  integer(int32),dimension(Nchi):: invschi
  integer(int32),dimension(Nchi,Nchi,2)::chi_map
  integer(int32),dimension(Nchi*(Nchi+1)/2,2)::irr_chi
  real(real64) norm,normb,inorm,norm2,weight
  real(real64),dimension(Norb,Norb):: sgnsig
  real(real64),dimension(Nchi,Nchi):: sgnsig2
  complex(real64),dimension(Nk,Nw,Nchi,Nchi):: chi
  complex(real64),dimension(Nkall,Nw,Norb,Norb):: newdelta,fk

  call get_invschi()
  call get_sgnsig()
  call get_chi_map_soc(chi_map,irr_chi,olist,invs,Nchi,Norb)
  weight=temp/dble(Nkall)
  norm2=0.0d0
  normb=0.0d0
  call get_chi0_conv_soc(chi,Gk,kmap,invk,invs,irr_chi,chi_map,olist,sgnsig,sgnsig2,temp,Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi)
  call ckchi()
  call get_V_delta_soc_flex(chi,Vmat,sgnsig2,Nk,Nw,Nchi)
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
        newdelta(:,:,:,:)=newdelta(:,:,:,:)*weight+delta(:,:,:,:)*norm2
        !$omp end parallel workshare
        call get_norm(norm,newdelta)
        inorm=1.0d0/norm
        if(abs(norm-norm2)>=1.0d2 .or. abs(norm-norm2)<1.0d-6)then
           print'(I3,A13,2E16.8)',i_iter,' lambda_elsh=',norm-norm2
        else
           print'(I3,A13,2F12.8)',i_iter,' lambda_elsh=',norm-norm2
        end if
        if(abs(norm-normb)*inorm<eps)then
           if((norm-norm2)>1.0d-1)then !do not finish until eig>0.1
              exit
           else if(.true.)then !consider small eig
              if(abs((norm-normb)*inorm)<eps*1.0d-2)then
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
 
  subroutine ckchi()
    integer(int32) i,l,m,n,info,chik,chikall,iorb
    real(real64) maxchi0,maxchi02
    real(real64),dimension(2*Nchi):: rwork
    complex(real64),dimension(Nchi*Nchi*4+1):: work
    complex(real64),dimension(Nchi):: eigc
    complex(real64),dimension(Nchi,Nchi):: chi0,tmp,tmp1

    maxchi02=-1.0d5
    do i=1,Nk
       chi0(:,:)=0.0d0
       do l=1,Nchi
          do m=1,Nchi
             do n=1,Nchi
                chi0(m,l)=chi0(m,l)+chi(i,1,m,n)*Vmat(n,l)
             end do
          end do
       end do
       call zgeev('N','N',Nchi,chi0,Nchi,eigc,tmp1,Nchi,tmp,Nchi,work,Nchi*Nchi*4+1,rwork,info)
       maxchi0=maxval(abs(dble(eigc)))
       if(maxchi0>maxchi02)then
          chik=i
          maxchi02=maxchi0
       end if
    end do
    do i=1,Nkall
       if(invk(2,i)==0)then
          if(invk(1,i)==chik)chikall=i
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
end subroutine lin_eliash_soc

subroutine get_V_delta_soc_flex(chi,Vmat,sgnsig2,Nk,Nw,Nchi)
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
end subroutine get_V_delta_soc_flex

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
  integer(int32),intent(in),dimension(Nchi):: invschi
  real(real64),intent(in),dimension(Norb,Norb):: sgnsig
  real(real64),intent(in),dimension(Nchi,Nchi):: Vmat,sgnsig2
  complex(real64),intent(in),dimension(Nk,Nw,Nchi,Nchi):: Vdelta
  complex(real64),intent(in),dimension(Nkall,Nw,Norb,Norb):: delta
  complex(real64),intent(out),dimension(Nkall,Nw,Norb,Norb):: newdelta

  integer(int32) i,j,k,n,l,m
  complex(real64),dimension(0:Nx-1,0:Ny-1,0:Nz-1,2*Nw):: tmpVdelta,tmpfk,tmp
  
  !$omp parallel workshare
  newdelta(:,:,:,:)=0.0d0
  !$omp end parallel workshare
  do l=1,Nchi
     do m=1,Nchi
        if((gap_sym<0) .or. (gap_sym>0 .and. slist(olist(l,1))*slist(olist(m,2))<0))then
           if(slist(olist(l,1))==1)then
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
                    newdelta(i,j,invs(olist(m,2)),invs(olist(l,1)))=newdelta(i,j,invs(olist(m,2)),invs(olist(l,1)))&
                         -sgnsig(olist(l,1),olist(m,2))*conjg(tmp(kmap(1,i),kmap(2,i),kmap(3,i),j))
                 end do
              end do
              !$omp end parallel do
           end if
        end if
     end do
  end do
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

subroutine get_vsigma_flex_soc(chi,Vmat,Nk,Nw,Nchi)
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nk,Nw,Nchi
  real(real64),intent(in),dimension(Nchi,Nchi):: Vmat
  complex(real64),intent(inout),dimension(Nk,Nw,Nchi,Nchi):: chi
 
  integer(int32) i,j,l,m,n,info
  integer(int32),dimension(Nchi):: ipiv
  complex(real64),dimension(Nchi,Nchi):: cmat1,cmat2,cmat3
  complex(real64),dimension(2*Nchi):: work
 
  do j=1,Nw
     do i=1,Nk
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
        cmat2(:,:)=cmat1(:,:)  !chi0V
        !$omp end workshare
        !$omp do 
        do l=1,Nchi
           cmat1(l,l)=cmat1(l,l)+1.0d0 !I-chi0S
        end do
        !$omp end do
        !$omp end parallel
        call zgetrf(Nchi,Nchi,cmat1,Nchi,ipiv,info)
        call zgetri(Nchi,cmat1,Nchi,ipiv,work,2*Nchi,info)
        !$omp parallel
        !$omp workshare
        cmat3(:,:)=0.0d0
        !$omp end workshare
        !$omp do private(m,n)
        do l=1,Nchi
           do m=1,Nchi
              do n=1,Nchi
                 cmat3(m,l)=cmat3(m,l)+cmat1(m,n)*cmat2(n,l) !(1-chi0S)^-1chi0S
              end do
           end do
        end do
        !$omp end do
        !$omp workshare
        cmat2(:,:)=-Vmat(:,:)
        !$omp end workshare
        !$omp do private(m,n)
        do l=1,Nchi
           do m=1,Nchi
              do n=1,Nchi
                 cmat2(m,l)=cmat2(m,l)+Vmat(m,n)*cmat3(n,l) !subtract double count 2nd order buble(ladder)
              end do
           end do
        end do
        !$omp end do
        !!$omp workshare
        chi(i,j,:,:)=cmat2(:,:)
        !!$omp end workshare
        !$omp end parallel
     end do
  end do
end subroutine get_vsigma_flex_soc

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
  do l=1,Norb
     do m=1,Norb
        if(gap_sym<0 .or. (gap_sym>0 .and. slist(m)*slist(l)<0))then !singlet gap only ud,du spin
           !$omp do private(i,n)
           do i=1,Nkall
              if(invk(2,i)==0)then
                 do n=1,Norb
                    delta(i,1,m,l)=delta(i,1,m,l)+uni(m,n,invk(1,i))*init_delta(i,n)*conjg(uni(l,n,invk(1,i)))
                 end do
              else if(invk(2,i)==1)then
                 do n=1,Norb
                    delta(i,1,m,l)=delta(i,1,m,l)+uni(l,n,invk(1,i))*init_delta(i,n)*slist(m)*conjg(uni(invs(m),n,invk(1,i)))
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

  chi_irr(:,:)=0
  chi_map(:,:,:)=0
  do l1=1,Nchi
     tmp1(1)=olist(l1,1) !1 of 1234
     tmp1(2)=olist(l1,2) !2 of 1234
     do m1=1,Nchi
        tmp1(3)=olist(m1,1) !3 of 1234
        tmp1(4)=olist(m1,2) !4 of 1234
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
                 exit
              end if
           end do
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
