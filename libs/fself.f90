subroutine FFT(cmat,tmp,Nx,Ny,Nz,Nw,SW)
  use,intrinsic::iso_fortran_env, only:int32,int64,real64
  implicit none
  integer(int64),intent(in):: Nx,Ny,Nz,Nw
  logical,intent(in):: SW
  complex(real64),intent(inout),dimension(Nx,Ny,Nz,Nw):: cmat,tmp
  
  integer(int64) plan
  integer(int32) Inv
  integer(int32),dimension(4):: Nlist

  Nlist=(/Nx,Ny,Nx,Nw/)
  if(SW)then
     Inv=-1
  else
     Inv=1
  end If
  print*,'init fft'
  call dfftw_plan_dft(plan,4,Nlist,cmat,tmp,Inv,64)
  print*,'execute'
  call dfftw_execute(plan)
  print*,'destroy'
  call dfftw_destroy_plan(plan)
  if(.not. SW) cmat=tmp/product(Nlist)
end subroutine FFT

subroutine gen_green0(Gk,eig,uni,mu,temp,Nk,Nw,Norb) bind(C,name="gen_green0_")
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  use constant
  implicit none
  integer(int64),intent(in):: Nk,Nw,norb
  real(real64),intent(in):: mu,temp
  real(real64),intent(in),dimension(Norb,Nk):: eig
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(real64),intent(out),dimension(Nk,Nw,Norb,Norb):: Gk
  
  integer(int32) i,j,l,m,n
  complex(real64) iw

  !$omp parallel private(l,m,n)
  do l=1,Norb
     do m=1,Norb
        band_loop: do n=1,Norb
           !$omp do private(iw,i,j)
           wloop: do j=1,Nw !ien=pi(2l+1)/beta l=0,1,...
              iw=cmplx(mu,dble(2*(j-1)+1)*pi*temp)
              kloop: do i=1,Nk
                 Gk(i,j,m,l)=Gk(i,j,m,l)+uni(l,n,i)*conjg(uni(m,n,i))/(iw-eig(n,i))
              end do kloop
           end do wloop
           !$omp end do
        end do band_loop
     end do
  end do
  !$omp end parallel
end subroutine gen_green0

subroutine gen_green_inv(Gk,self,hamk,mu,temp,Nk,Nw,Norb) bind(C,name="gen_green_inv_")
  use,intrinsic:: iso_fortran_env, only: int64,real64,int32
  use constant
  implicit none
  integer(int64),intent(in):: Nk,Nw,Norb
  real(real64),intent(in):: mu,temp
  complex(real64),intent(in),dimension(norb,norb,Nk):: hamk
  complex(real64),intent(in),dimension(Nk,Nw,norb,norb):: self
  complex(real64),intent(out),dimension(Nk,Nw,norb,norb):: Gk

  integer(int32)i,j,l,m
  complex(real64) iw

  !G^-1=G^-1_0-sigma (=iwI-Hk-sigma)
  !$omp parallel private(l,m)
  do l=1,Norb
     do m=1,Norb
        !$omp do private(iw,i)
        do j=1,Nw
           iw=cmplx(mu,dble(2*(j-1)+1)*pi*temp)
           do i=1,Nk
              Gk(i,j,m,l)=-(hamk(m,l,i)+self(i,j,m,l))
              if(l==m)then
                 Gk(i,j,l,l)=Gk(i,j,l,l)+iw
              end if
           end do
        end do
        !$omp end do
     end do
  end do
  !$omp end parallel
end subroutine gen_green_inv

subroutine gen_green_inv_from_eig(Gk,self,uni,eig,mu,temp,Nk,Nw,Norb) bind(C)
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  use constant
  implicit none
  integer(int64),intent(in):: Nk,Nw,norb
  real(real64),intent(in):: mu,temp
  real(real64),intent(in),dimension(norb,Nk):: eig
  complex(real64),intent(in),dimension(norb,norb,Nk):: uni
  complex(real64),intent(in),dimension(Nk,Nw,norb,norb):: self
  complex(real64),intent(out),dimension(Nk,Nw,norb,norb):: Gk

  integer(int32)i,j,l,m,n
  complex(real64) iw

  !$omp parallel private(l,m)
  do l=1,norb
     do m=1,norb
        !$omp do private(iw,i,n)
        wloop: do j=1,Nw
           iw=cmplx(mu,dble(2*(j-1)+1)*pi*temp)
           kloop: do i=1,Nk
              band_loop: do n=1,norb
                 Gk(i,j,m,l)=uni(m,n,i)*conjg(uni(l,n,i))*(iw-eig(n,i))
              end do band_loop
              Gk(i,j,m,l)=Gk(i,j,m,l)-self(i,j,m,l)
           end do kloop
        end do wloop
        !$omp end do
     end do
  end do
  !$omp end parallel
end subroutine gen_green_inv_from_eig

subroutine getinv(Gk,Nk,Nw,Norb) bind(C,name="getinv_")
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nk,Nw,Norb
  complex(real64),intent(inout),dimension(Nk,Nw,Norb,Norb):: Gk

  integer(int32) i,j,ipiv(Norb),info
  complex(real64) tmp(Norb,Norb),work(2*Norb)

  !$omp parallel do private(i,j,tmp,work,ipiv,info)
  do i=1,Nw
     do j=1,Nk
        tmp(:,:)=Gk(j,i,:,:)
        call zgetrf(Norb,Norb,tmp,Norb,ipiv,info)
        call zgetri(Norb,tmp,Norb,ipiv,work,2*Norb,info)
     end do
  end do
  !$omp end parallel do
end subroutine getinv

subroutine get_chi0_conv(chi,Gk,kmap,olist,temp,Nx,Ny,Nz,Nw,Nk,Norb,Nchi) bind(C)
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nw,Norb,Nchi,Nk,Nx,Ny,Nz
  integer(int64),intent(in),dimension(Nchi,2):: olist
  integer(int64),intent(in),dimension(3,Nk):: kmap
  real(real64),intent(in):: temp
  complex(real64),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(real64),intent(out),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(int32) i,j,k,l,m,n
  integer(int32) ii(0:Nx-1),ij(0:Ny-1),ik(0:Nz-1)
  real(real64) weight
  complex(real64),dimension(0:Nx-1,0:Ny-1,0:Nz-1,2*Nw):: tmpchi,tmp,tmp1,tmp2

  weight=temp/dble(Nk)
  ii(0)=0
  ij(0)=0
  ik(0)=0
  !$omp parallel
  !$omp do
  do i=1,Nx-1
     ii(i)=Nx-i
  end do
  !$omp end do
  !$omp do
  do i=1,Ny-1
     ij(i)=Ny-j
  end do
  !$omp end do
  !$omp do
  do i=1,Nz-1
     ik(i)=Nz-k
  end do
  !$omp end do
  !$omp end parallel
  chi_orb_loop1:do l=1,Nchi
     chi_orb_loo2:do m=1,Nchi
        !use symmetry G^lm(k,iw)=G^ml(k,-iw)
        !$omp parallel do private(i)
        w_loop_Gk_to_tmp:do j=1,Nw
           k_loop_Gk_to_tmp:do i=1,Nk
              tmp1(kmap(1,i),kmap(2,i),kmap(3,i),j)=Gk(i,j,olist(l,1),olist(m,1)) !G13(iw)
              tmp2(kmap(1,i),kmap(2,i),kmap(3,i),j)=Gk(i,j,olist(m,2),olist(l,2)) !G42(iw)
              tmp1(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=conjg(Gk(i,j,olist(m,1),olist(l,1))) !G13(-iw)=G^*31(iw)
              tmp2(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=conjg(Gk(i,j,olist(l,2),olist(m,2))) !G42(-iw)=G^*24(iw)
           end do k_loop_Gk_to_tmp
        end do w_loop_Gk_to_tmp
        !$omp end parallel do
        call FFT(tmp1,tmp,Nx,Ny,Nz,2*Nw,.true.)
        call FFT(tmp2,tmp,Nx,Ny,Nz,2*Nw,.true.)
        !calculate G(r)G(-r)
        !$omp parallel do private(i,j,k)
        w_loop_conv:do n=1,2*Nw
           z_loop:do k=0,Nz-1
              y_loop:do j=0,Ny-1
                 !$omp simd
                 x_loop:do i=0,Nx-1
                    tmpchi(i,j,k,n)=tmp1(i,j,k,n)*tmp2(ii(i),ij(j),ik(k),2*nw-n+1)
                 end do x_loop
                 !$omp end simd
              end do y_loop
           end do z_loop
        end do w_loop_conv
        !$omp end parallel do
        !call FFT(tmpchi,tmp,Nx,Ny,Nz,2*Nw,.false.)
        !$omp parallel do private(i,j)
        w_loop_tmp_to_chi:do j=1,Nw
           k_loop_tmp_to_chi:do i=1,Nk
              chi(i,j,m,l)=tmpchi(kmap(1,i),kmap(2,i),kmap(3,i),j)
           end do k_loop_tmp_to_chi
        end do w_loop_tmp_to_chi
        !$omp end parallel do
     end do chi_orb_loo2
  end do chi_orb_loop1
  chi(:,:,:,:)=chi(:,:,:,:)*weight
end subroutine get_chi0_conv

subroutine get_chi0_sum(chi,Gk,klist,olist,temp,Nw,Nk,Norb,Nchi) bind(C)
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nw,Norb,Nchi,Nk
  integer(int64),intent(in),dimension(Nchi,2):: olist
  real(real64),intent(in),dimension(3,Nk):: klist
  real(real64),intent(in):: temp
  complex(real64),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(real64),intent(out),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(int32) i,j,l,m,iq,iw
  integer(int32),dimension(2*Nw):: wshift
  integer(int64),dimension(Nk):: qshift
  real(real64) weight
  complex(real64),dimension(Nk,2*Nw):: tmpgk13,tmpgk42

  weight=temp/dble(Nk)
  chi(:,:,:,:)=0.0d0
  do l=1,Nchi
     do m=1,Nchi
        !$omp parallel
        !$omp do private(i,j)
        do j=1,Nw
           do i=1,Nk
              tmpgk13(i,j)=Gk(i,j,olist(m,1),olist(l,1))
              tmpgk42(i,j)=Gk(i,j,olist(l,2),olist(m,2))
              tmpgk13(i,2*Nw-j+1)=conjg(Gk(i,j,olist(l,1),olist(m,1)))
              tmpgk42(i,2*Nw-j+1)=conjg(Gk(i,j,olist(m,2),olist(l,2)))
           end do
        end do
        !$omp end do
        !$omp do private(iw,j,iq,i,wshift,qshift)
        wloop: do iw=1,Nw
           do j=1,2*Nw
              wshift(j)=mod(j+iw-1,2*Nw)
              if(wshift(j)==2*Nw)then
                 wshift(j)=2*Nw
              end if
           end do
           qloop: do iq=1,Nk
              call get_qshift(klist(:,iq),klist,qshift,Nk)
              wloop2: do j=1,2*Nw
                 kloop: do i=1,Nk
                    chi(iq,iw,m,l)=chi(iq,iw,m,l)-tmpgk13(i,j)*tmpgk42(qshift(i),wshift(j))
                 end do kloop
              end do wloop2
           end do qloop
        end do wloop
        !$omp end do
        !$omp end parallel
     end do
  end do
  chi(:,:,:,:)=chi(:,:,:,:)*weight
end subroutine get_chi0_sum

subroutine get_trace_chi(trchi,chi,Nk,Nw,Nchi) bind(C)
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nk,Nw,Nchi
  complex(real64),intent(in),dimension(Nk,Nw,Nchi,Nchi):: chi
  complex(real64),intent(out),dimension(Nk,Nw):: trchi

  integer(int32) i,j,k

  !trchi(:,:)=0.0d0
  do k=1,Nchi
     do j=1,Nw
        do i=1,Nk
           trchi(i,j)=trchi(i,j)+chi(i,j,k,k)
        end do
     end do
  end do
end subroutine get_trace_chi

subroutine mkself(selfen,hamk,eig,uni,mu,temp,scf_loop,eps,Nk,Nw,Norb) bind(C)
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nk,Nw,Norb,scf_loop
  real(real64),intent(in):: mu,temp,eps
  real(real64),intent(in),dimension(Norb,Nk):: eig
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni,hamk
  complex(real64),intent(out),dimension(Nk,Nw,Norb,Norb):: selfen

  integer(int32) scf_i
  complex(real64),dimension(Nk,Nw,Norb,Norb):: Gk,selfen0
  
  selfen0(:,:,:,:)=0.0d0
  call gen_green0(Gk,eig,uni,mu,temp,Nk,Nw,Norb)
  do scf_i=1,scf_loop
     print*,'iter=',scf_i
     selfen(:,:,:,:)=-1.5*Gk(:,:,:,:)
     print*,'gen self energy'
     if(sum(abs(selfen(:,:,:,:)-selfen0(:,:,:,:)))<eps)then
        exit
     end if
     call gen_green_inv(Gk,selfen,hamk,mu,temp,Nk,Nw,Norb)
     print*,'calc gk inv'
     call getinv(Gk,Nk,Nw,Norb)
     print*,'calc gk'
     selfen0(:,:,:,:)=selfen(:,:,:,:)
  end do
end subroutine mkself
