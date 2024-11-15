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
        !$omp workshare
        chi(i,j,:,:)=cmat4(:,:)
        !$omp end workshare
        !$omp end parallel
     end do qloop
  end do wloop
end subroutine get_V_delta_nsoc_flex

subroutine get_V_delta_soc_flex(chi,Vmat,Nk,Nw,Nchi)
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nk,Nw,Nchi
  real(real64),intent(in),dimension(Nchi,Nchi):: Vmat
  complex(real64),intent(inout),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(int32) i,j,l,m,n,info
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
        !$omp do private(m,n)
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
        !$omp do private(m,n)
        do l=1,Nchi
           do m=1,Nchi
              do n=1,Nchi
                 cmat2(m,l)=cmat2(m,l)+Vmat(m,n)*cmat3(n,l)
              end do
           end do
        end do
        !$omp end do
        !$omp workshare
        chi(i,j,:,:)=cmat2(:,:)
        !$omp end workshare
        !$omp end parallel
     end do qloop
  end do wloop
end subroutine get_V_delta_soc_flex

subroutine lin_eliash(delta,Gk,uni,Smat,Cmat,olist,kmap,invk,temp,eps,&
     Nk,Nw,Nchi,Norb,Nx,Ny,Nz,itemax,gap_sym) bind(C)
  !
  !> calculate linearized eliashberg equations
  !
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nk,Nw,Nchi,Norb,Nx,Ny,Nz,itemax,gap_sym
  integer(int64),intent(in),dimension(Nchi,2):: olist
  integer(int64),intent(in),dimension(3,Nk):: kmap
  integer(int64),intent(in),dimension(Nk):: invk
  real(real64),intent(in):: temp,eps
  real(real64),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  complex(real64),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(real64),intent(out),dimension(Nk,Nw,Norb,Norb):: delta

  integer(int32) i_iter,i_eig,count,i
  integer(int32),parameter:: eig_max=2
  logical(1) sw_pair
  real(real64) norm,normb,inorm,norm2,weight
  complex(real64),dimension(Nk,Nw,Nchi,Nchi):: chi
  complex(real64),dimension(Nk,Nw,Norb,Norb):: newdelta,fk
 
  if(gap_sym>=0)then
     sw_pair=.true.
     print'(A7)','singlet'
  else
     sw_pair=.false.
     print'(A7)','triplet'
  end if
  weight=temp/dble(Nk)
  norm2=0.0d0
  normb=0.0d0
  call get_chi0_conv(chi,Gk,kmap,olist,temp,Nx,Ny,Nz,Nw,Nk,Norb,Nchi)
  call ckchi()
  call get_V_delta_nsoc_flex(chi,Smat,Cmat,Nk,Nw,Nchi,sw_pair)
  print'(A15,2E16.8)','V_delta max is ',maxval(dble(chi)),maxval(aimag(chi))
  print'(A15,2E16.8)','V_delta min is ',minval(dble(chi)),minval(aimag(chi))
  eigenval_loop:do i_eig=1,eig_max !solve eig_val using power method, 1st eig is usually large negative value
     call get_initial_delta(delta,uni,kmap,Nk,Nw,Norb,Nx,Ny,Nz,gap_sym)
     count=0 !count too small eigenvalue
     iter_loop:do i_iter=1,itemax !iteration
        call mkfk_trs_nsoc(fk,Gk,delta,Nk,Nw,Norb)
        call mkdelta_nsoc(newdelta,fk,chi,Smat,Cmat,kmap,invk,olist,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,sw_pair)
        !$omp parallel workshare
        newdelta(:,:,:,:)=newdelta(:,:,:,:)*weight+delta(:,:,:,:)*norm2
        !$omp end parallel workshare
        call get_norm()
        inorm=1.0d0/norm
        if(abs(norm-norm2)>=1.0d2 .or. abs(norm-norm2)<1.0d-6)then
           print'(I3,A13,E16.8)',i_iter,' lambda_elsh=',norm-norm2
        else
           print'(I3,A13,F12.8)',i_iter,' lambda_elsh=',norm-norm2
        end if
        if(abs((norm-normb)*inorm)<eps)then
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
contains
  subroutine get_norm()
    integer(int32) i,j,l,m
    real(real64) tmp

    tmp=0.0d0
    !$omp parallel
    do l=1,Norb
       do m=1,Norb
          !$omp do private(i)
          do j=1,Nw
             do i=1,Nk
                tmp=tmp+2.0d0*abs(newdelta(i,j,m,l))*abs(newdelta(i,j,m,l))
             end do
          end do
          !$omp end do
       end do
    end do
    !$omp end parallel
    norm=sqrt(tmp)
  end subroutine get_norm
  
  subroutine ckchi()
    integer(int32) i,l,m,n,info,chisk,chick
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
    print'(A3,3I4,F12.8)','SDW',kmap(:,chisk),maxchi0s2
    print'(A3,3I4,F12.8)','CDW',kmap(:,chick),maxchi0c2
  end subroutine ckchi
end subroutine lin_eliash

subroutine mkfk_trs_nsoc(fk,Gk,delta,Nk,Nw,Norb)
  !
  !>calculate linearized anomalous green function Fk with TRS
  !>if we considder TRS,F_ab(k)=G_ac(k)Delta_cd(k)Gbd(-k)
  !> =G_ac(k)Delta_cd(k)G^*_dbss'(k)ss'
  !
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
  !$omp end workshare
  do l=1,Norb
     do m=1,Norb
        do n=1,Norb
           !$omp do private(i)
           do j=1,Nw
              do i=1,Nk
                 cmat1(i,j,m,l)=cmat1(i,j,m,l)+Gk(i,j,m,n)*delta(i,j,n,l)
              end do
           end do
           !$omp end do
        end do
     end do
  end do
  !$omp workshare
  fk(:,:,:,:)=0.0d0
  !$omp end workshare
  do l=1,Norb
     do m=1,Norb
        do n=1,Norb
           !$omp do private(i)
           do j=1,Nw
              do i=1,Nk
                 fk(i,j,m,l)=fk(i,j,m,l)-cmat1(i,j,m,n)*conjg(Gk(i,j,l,n))
              end do
           end do
           !$omp end do
        end do
     end do
  end do
  !$omp end parallel
end subroutine mkfk_trs_nsoc

subroutine mkdelta_nsoc(newdelta,delta,Vdelta,Smat,Cmat,kmap,invk,olist,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,sw_pair)
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nk,Nw,Nchi,Norb,Nx,Ny,Nz
  integer(int64),intent(in),dimension(3,Nk):: kmap
  integer(int64),intent(in),dimension(Nk):: invk
  integer(int64),intent(in),dimension(Nchi,2):: olist
  logical(1),intent(in):: sw_pair
  real(real64),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  complex(real64),intent(in),dimension(Nk,Nw,Nchi,Nchi):: Vdelta
  complex(real64),intent(in),dimension(Nk,Nw,Norb,Norb):: delta
  complex(real64),intent(out),dimension(Nk,Nw,Norb,Norb):: newdelta
  
  integer(int32) i,j,k,n,l,m
  real(real64) sgn
  complex(real64),dimension(0:Nx-1,0:Ny-1,0:Nz-1,2*Nw):: tmpVdelta,tmpfk,tmp

  if(sw_pair)then
     sgn=1.0d0 !singlet Fk=F-k
  else
     sgn=-1.0d0 !triplet Fk=-F-k
  end if
  newdelta(:,:,:,:)=0.0d0
  do l=1,Nchi
     do m=1,Nchi
        !$omp parallel
        !$omp do
        do i=1,Nk
           tmpVdelta(kmap(1,i),kmap(2,i),kmap(3,i),1)=Vdelta(i,1,m,l) !j=1 corresponds to w_n=0
           tmpVdelta(kmap(1,i),kmap(2,i),kmap(3,i),Nw+1)=0.5d0*(Cmat(m,l)+sgn*Smat(m,l)) !Nw+1 consider w_n=>inf limit
        end do
        !$omp end do
        !$omp do private(i)
        do j=2,Nw
           do i=1,Nk
              tmpVdelta(kmap(1,i),kmap(2,i),kmap(3,i),j)=Vdelta(i,j,m,l)
              tmpVdelta(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+2)=conjg(Vdelta(i,j,l,m))
           end do
        end do
        !$omp end do
        !$omp do private(i)
        do j=1,Nw
           do i=1,Nk
              tmpfk(kmap(1,i),kmap(2,i),kmap(3,i),j)=delta(i,j,olist(m,1),olist(l,2))
              tmpfk(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=sgn*delta(invk(i),j,olist(l,2),olist(m,1)) !F(k,-w)=sgn*F(-k,w) (no soc only)
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
           do i=1,Nk
              newdelta(i,j,olist(m,2),olist(l,1))=newdelta(i,j,olist(m,2),olist(l,1))&
                   +tmp(kmap(1,i),kmap(2,i),kmap(3,i),j)
           end do
        end do
        !$omp end parallel do
     end do
  end do
end subroutine mkdelta_nsoc

subroutine get_initial_delta(delta,uni,kmap,Nk,Nw,Norb,Nx,Ny,Nz,gap_sym)
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  use constant
  implicit none
  integer(int64),intent(in):: Nw,Norb,Nk,Nx,Ny,Nz,gap_sym
  integer(int64),intent(in),dimension(3,Nk):: kmap
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(real64),intent(out),dimension(Nk,Nw,Norb,Norb):: delta

  integer(int32) i,j,l,m,n
  real(real64) norm
  complex(real64),dimension(Nk,Norb):: deltab
  if(gap_sym==0)then
     deltab(:,:)=1.0d0
  else
     do l=1,Norb
        do i=1,Nk
           select case(gap_sym)
           case(1)
              deltab(i,l)=cos(2*pi*dble(kmap(1,i))/Nx)-cos(2*pi*dble(kmap(2,i))/Ny) !dx2-y2
           case(2)
              deltab(i,l)=2.0d0*cos(2*pi*dble(kmap(1,i))/Nx)*cos(2*pi*dble(kmap(2,i))/Ny) !spm
           case default !same as 0
              deltab(i,l)=1.0d0
           end select
        end do
     end do
  end if

  do l=1,Norb
     do m=1,Norb
        do n=1,Norb
           do j=1,Nw
              do i=1,Nk
                 delta(i,j,m,l)=delta(i,j,m,l)+conjg(uni(m,n,i))*uni(l,n,i)*deltab(i,n)
              end do
           end do
        end do
     end do
  end do
  norm=2.0d0*sum(abs(delta(:,:,:,:))**2)
  delta(:,:,:,:)=delta(:,:,:,:)/norm
end subroutine get_initial_delta

subroutine conv_delta_orb_to_band(deltab,delta,uni,Norb,Nk,Nw)
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nw,Norb,Nk
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(real64),intent(in),dimension(Nk,Nw,Norb,Norb):: delta
  complex(real64),intent(out),dimension(Nk,Norb,Norb):: deltab

  integer(int32) i,l,m,n
  complex(real64),dimension(Nk,Norb,Norb):: tmp

  !$omp parallel
  !$omp workshare
  deltab(:,:,:)=0.0d0
  tmp(:,:,:)=0.0d0
  !$omp end workshare
  !$omp do private(m,n,i)
  do l=1,Norb
     do m=1,Norb
        do n=1,Norb
           do i=1,Nk
              tmp(i,m,l)=tmp(i,m,l)+conjg(uni(m,n,i))*delta(i,1,n,l)
           end do
        end do
     end do
  end do
  !$omp end do
  !$omp do private(m,n,i)
  do l=1,Norb
     do m=1,Norb
        do n=1,Norb
           do i=1,Nk
              deltab(i,m,l)=deltab(i,m,l)+tmp(i,m,n)*uni(n,l,i)
           end do
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine conv_delta_orb_to_band
