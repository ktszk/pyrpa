subroutine lin_eliash(delta,Gk,uni,init_delta,Smat,Cmat,olist,kmap,invk,temp,eps,&
     Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,itemax,gap_sym) bind(C)
  !> calculate linearized eliashberg equations
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,itemax,gap_sym
  integer(int64),intent(in),dimension(Nchi,2):: olist
  integer(int64),intent(in),dimension(3,Nkall):: kmap,invk
  real(real64),intent(in):: temp,eps
  real(real64),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  real(real64),intent(in),dimension(Nkall,Norb):: init_delta
  complex(real64),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(real64),intent(out),dimension(Nkall,Nw,Norb,Norb):: delta

  integer(int32) i_iter,i_eig,count,i
  integer(int32),parameter:: eig_max=2
  integer(int32),dimension(Nchi,Nchi,2)::chi_map
  integer(int32),dimension(Nchi*(Nchi+1)/2,2)::irr_chi
  logical(1) sw_pair
  real(real64) norm,normb,inorm,norm2,weight
  complex(real64),dimension(Nk,Nw,Nchi,Nchi):: chi
  complex(real64),dimension(Nkall,Nw,Norb,Norb):: newdelta,fk

  if(gap_sym>=0)then
     sw_pair=.true.
     print'(A7)','singlet'
  else
     sw_pair=.false.
     print'(A7)','triplet'
  end if
  call get_chi_map(chi_map,irr_chi,olist,Nchi)
  weight=temp/dble(Nkall)
  norm2=0.0d0
  normb=0.0d0
  call get_chi0_conv(chi,Gk,kmap,invk,irr_chi,chi_map,olist,temp,Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi)
  call ckchi()
  call get_V_delta_nsoc_flex(chi,Smat,Cmat,Nk,Nw,Nchi,sw_pair)
  print'(A15,2E16.8)','V_delta max is ',maxval(dble(chi)),maxval(aimag(chi))
  print'(A15,2E16.8)','V_delta min is ',minval(dble(chi)),minval(aimag(chi))
  eigenval_loop:do i_eig=1,eig_max !solve eig_val using power method, 1st eig is usually large negative value
     call get_initial_delta(delta,init_delta,uni,kmap,invk,Nkall,Nk,Nw,Norb,gap_sym)
     count=0 !count too small eigenvalue
     iter_loop:do i_iter=1,itemax !iteration
        call mkfk_trs_nsoc(fk,Gk,delta,invk,Nkall,Nk,Nw,Norb)
        call mkdelta_nsoc(newdelta,fk,chi,Smat,Cmat,kmap,invk,olist,Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,sw_pair)
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
    integer(int32) i,l,m,n,info,chisk,chick,chiskall,chickall
    real(real64) maxchi0s,maxchi0c,maxchi0s2,maxchi0c2
    real(real64),dimension(2*Nchi):: rwork
    complex(real64),dimension(Nchi*Nchi*4+1):: work
    complex(real64),dimension(Nchi):: eigs,eigc
    complex(real64),dimension(Nchi,Nchi):: chi0s,chi0c,tmp1,tmp2

    maxchi0s2=-1.0d5
    maxchi0c2=-1.0d5
    do i=1,Nk
       chi0s(:,:)=0.0d0
       chi0c(:,:)=0.0d0
       do l=1,Nchi
          do m=1,Nchi
             do n=1,Nchi
                chi0s(m,l)=chi0s(m,l)+chi(i,1,m,n)*Smat(n,l)
                chi0c(m,l)=chi0c(m,l)-chi(i,1,m,n)*Cmat(n,l)
             end do
          end do
       end do
       call zgeev('N','N',Nchi,chi0s,Nchi,eigs,tmp1,Nchi,tmp2,Nchi,work,Nchi*Nchi*4+1,rwork,info)
       call zgeev('N','N',Nchi,chi0c,Nchi,eigc,tmp1,Nchi,tmp2,Nchi,work,Nchi*Nchi*4+1,rwork,info)
       maxchi0s=maxval(dble(eigs))
       maxchi0c=maxval(dble(eigc))
       if(maxchi0s>maxchi0s2)then
          chisk=i
          maxchi0s2=maxchi0s
       end if
       if(maxchi0c>maxchi0c2)then
          chick=i
          maxchi0c2=maxchi0c
       end if
    end do
    do i=1,Nkall
       if(invk(2,i)==0)then
          if(invk(1,i)==chisk)chiskall=i
          if(invk(1,i)==chick)chickall=i
       end if
    end do
    print'(A3,3I4,F12.8)','SDW',kmap(:,chiskall),maxchi0s2
    print'(A3,3I4,F12.8)','CDW',kmap(:,chickall),maxchi0c2
  end subroutine ckchi
end subroutine lin_eliash

subroutine get_V_delta_nsoc_flex(chi,Smat,Cmat,Nk,Nw,Nchi,sw_pair)
  !
  !> obtain pairing interaction V_delta without soc
  !
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
        call zgetri(Nchi,cmat1,Nchi,ipiv,work,2*Nchi,info)
        call zgetrf(Nchi,Nchi,cmat2,Nchi,ipiv,info)
        call zgetri(Nchi,cmat2,Nchi,ipiv,work,2*Nchi,info)
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
                 cmat5(m,l)=cmat5(m,l)+cmat2(m,n)*cmat4(n,l) !(1+chi0C)^-1chi0S
              end do
           end do
        end do
        !$omp end do
        !$omp workshare
        cmat2(:,:)=cmat5(:,:) !cmat2=(1+chi_0C)^--1chi_0C
        !$omp end workshare
        if(sw_pair)then !singlet
           !$omp workshare
           cmat4(:,:)=0.5d0*(Smat(:,:)+Cmat(:,:)) !Vud=(C+S)/2
           !$omp end workshare
           !$omp do private(m,n)
           do l=1,Nchi
              do m=1,Nchi
                 do n=1,Nchi
                    cmat4(m,l)=cmat4(m,l)+1.5d0*Smat(m,n)*cmat1(n,l)-0.5d0*Cmat(m,n)*cmat2(n,l)
                 end do
              end do
           end do
           !$omp end do
        else !triplet
           !$omp workshare
           cmat4(:,:)=0.5d0*(Smat(:,:)-Cmat(:,:)) !Vuu=(C-S)/2
           !$omp end workshare
           !$omp do private(m,n)
           do l=1,Nchi
              do m=1,Nchi
                 do n=1,Nchi
                    cmat4(m,l)=cmat4(m,l)-0.5d0*(Smat(m,n)*cmat1(n,l)+Cmat(m,n)*cmat2(n,l))
                 end do
              end do
           end do
           !$omp end do
        end if
        !!$omp workshare
        chi(i,j,:,:)=cmat4(:,:)
        !!$omp end workshare
        !$omp end parallel
     end do qloop
  end do wloop
end subroutine get_V_delta_nsoc_flex

subroutine mkfk_trs_nsoc(fk,Gk,delta,invk,Nkall,Nk,Nw,Norb) bind(C,name="mkfk_trs_nsoc_")
  !
  !>calculate linearized anomalous green function Fk with TRS
  !>if we considder TRS,F_ab(k)=G_ac(k)Delta_cd(k)Gbd(-k)
  !> =G_ac(k)Delta_cd(k)G^*_dbss'(k)ss'
  !
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nkall,Nk,Nw,Norb
  integer(int64),intent(in),dimension(3,Nkall):: invk
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
                    cmat1(i,j,m,l)=cmat1(i,j,m,l)+Gk(invk(1,i),j,n,m)*delta(i,j,n,l)
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
                    fk(i,j,m,l)=fk(i,j,m,l)-cmat1(i,j,m,n)*conjg(Gk(invk(1,invk(3,i)),j,n,l))
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
end subroutine mkfk_trs_nsoc

subroutine mkdelta_nsoc(newdelta,delta,Vdelta,Smat,Cmat,kmap,invk,olist,Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,sw_pair)
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz
  integer(int64),intent(in),dimension(3,Nkall):: kmap,invk
  integer(int64),intent(in),dimension(Nchi,2):: olist
  logical(1),intent(in):: sw_pair
  real(real64),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  complex(real64),intent(in),dimension(Nk,Nw,Nchi,Nchi):: Vdelta
  complex(real64),intent(in),dimension(Nkall,Nw,Norb,Norb):: delta
  complex(real64),intent(out),dimension(Nkall,Nw,Norb,Norb):: newdelta
  
  integer(int32) i,j,k,n,l,m
  real(real64) sgn
  complex(real64),dimension(0:Nx-1,0:Ny-1,0:Nz-1,2*Nw):: tmpVdelta,tmpfk,tmp

  if(sw_pair)then
     sgn=1.0d0 !singlet Fk=F^T-k
  else
     sgn=-1.0d0 !triplet Fk=-F^T-k
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
                 tmpfk(kmap(1,i),kmap(2,i),kmap(3,i),j)=delta(i,j,olist(l,2),olist(m,1))
                 tmpfk(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=sgn*delta(invk(3,i),j,olist(m,1),olist(l,2)) !F(k,-w)=sgn*F(-k,w) (no soc only)
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
     end do
  end do

  !$omp parallel
  do l=1,Norb
     do m=l,Norb
        !$omp do private(i,j)
        do j=1,Nw
           do i=1,Nkall
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
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  use constant
  implicit none
  integer(int64),intent(in):: Nw,Norb,Nkall,Nk,gap_sym
  integer(int64),intent(in),dimension(3,Nkall):: kmap,invk
  real(real64),intent(in),dimension(Nkall,Norb):: init_delta
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(real64),intent(out),dimension(Nkall,Nw,Norb,Norb):: delta

  integer(int32) i,j,l,m,n
  real(real64) norm

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
     !$omp parallel
     !$omp workshare
     delta(:,:,:,:)=0.0d0
     !$omp end workshare
     do l=1,Norb
        do m=l,Norb
           !$omp do private(i,n)
           do i=1,Nkall
              if(invk(2,i)==0)then
                 do n=1,Norb
                    delta(i,1,m,l)=delta(i,1,m,l)+uni(m,n,invk(1,i))*init_delta(i,n)*conjg(uni(l,n,invk(1,i)))
                 end do
              else if(invk(2,i)==1)then
                 do n=1,Norb
                    delta(i,1,m,l)=delta(i,1,m,l)+uni(l,n,invk(1,i))*init_delta(i,n)*conjg(uni(m,n,invk(1,i)))
                 end do
              end if
           end do
           !$omp end do
           !$omp do private(i,j)
           do j=2,Nw
              do i=1,Nkall
                 delta(i,j,m,l)=delta(i,1,m,l)
              end do
           end do
           !$omp end do
           !$omp do private(i,j)
           do j=1,Nw
              do i=1,Nkall
                 if(l/=m)then !TRS gap function is Hermite
                    delta(i,j,l,m)=conjg(delta(i,j,m,l))
                 else !TRS diagonal gap is real
                    delta(i,j,l,l)=dble(delta(i,j,l,l))
                 end if
              end do
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
end subroutine get_initial_delta

subroutine conv_delta_orb_to_band(deltab,delta,uni,invk,Norb,Nkall,Nk,Nw) bind(C)
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nw,Norb,Nk,Nkall
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
                 tmp(i,m,l)=tmp(i,m,l)+uni(n,m,invk(1,i))*delta(i,1,n,l)
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
