subroutine FFT(cmat,tmp,Nx,Ny,Nz,Nw,SW)
  use,intrinsic::iso_fortran_env, only:int32,int64,real64
  implicit none
  integer(int64),intent(in):: Nx,Ny,Nz,Nw
  logical,intent(in):: SW
  complex(real64),intent(inout),dimension(Nx,Ny,Nz,Nw):: cmat,tmp
  
  integer(int64) plan
  integer(int32) Inv
  integer(int32),dimension(4):: Nlist

  Nlist=(/Nx,Ny,Nz,Nw/)
  if(SW)then
     Inv=-1
  else
     Inv=1
  end If
  call dfftw_plan_dft(plan,4,Nlist,cmat,tmp,Inv,64)
  call dfftw_execute(plan)
  call dfftw_destroy_plan(plan)
  if(.not. SW)then
     cmat(:,:,:,:)=tmp(:,:,:,:)/product(Nlist)
  else
     cmat(:,:,:,:)=tmp(:,:,:,:)
  end if
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
           wloop: do j=1,Nw !ien=pi(2l+1)/beta l=0,1,... beta=(kBT)^-1
              iw=cmplx(mu,dble(2*(j-1)+1)*pi*temp) !j=1=>l=0
              kloop: do i=1,Nk
                 Gk(i,j,m,l)=Gk(i,j,m,l)+uni(m,n,i)*conjg(uni(l,n,i))/(iw-eig(n,i))
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
  !$omp workshare
  Gk(:,:,:,:)=0.0d0
  !$omp end workshare
  do l=1,norb
     do m=1,norb
        !$omp do private(iw,i,n)
        wloop: do j=1,Nw
           iw=cmplx(mu,dble(2*(j-1)+1)*pi*temp)
           kloop: do i=1,Nk
              band_loop: do n=1,norb
                 Gk(i,j,m,l)=Gk(i,j,m,l)+uni(m,n,i)*conjg(uni(l,n,i))*(iw-eig(n,i))
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
        Gk(j,i,:,:)=tmp(:,:)
     end do
  end do
  !$omp end parallel do
end subroutine getinv

subroutine get_chi0_conv(chi,Gk,kmap,invk,olist,temp,Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi) bind(C,name='get_chi0_conv_')
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nw,Norb,Nchi,Nkall,Nk,Nx,Ny,Nz
  integer(int64),intent(in),dimension(Nchi,2):: olist
  integer(int64),intent(in),dimension(3,Nkall):: kmap
  integer(int64),intent(in),dimension(2,Nkall):: invk
  real(real64),intent(in):: temp
  complex(real64),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(real64),intent(out),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(int32) i,j,k,l,m,n
  integer(int32) ii(0:Nx-1),ij(0:Ny-1),ik(0:Nz-1)
  real(real64) weight
  complex(real64),dimension(0:Nx-1,0:Ny-1,0:Nz-1,2*Nw):: tmp,tmpgk13,tmpgk42
  
  weight=temp/dble(Nkall)
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
     ij(i)=Ny-i
  end do
  !$omp end do
  !$omp do
  do i=1,Nz-1
     ik(i)=Nz-i
  end do
  !$omp end do
  !$omp end parallel
  l1_l2_loop1:do l=1,Nchi !olist(l)=(l1,l2)
     l3_l4_loo2:do m=1,Nchi !olist(m)=(l3,l4)
        !use symmetry G^lm(k,iw)=G^ml(k,-iw) from Hermitian symmetry of Hamiltonian
        !$omp parallel do private(i)
        w_loop_Gk_to_tmp:do j=1,Nw
           k_loop_Gk_to_tmp:do i=1,Nkall
              if(invk(2,i)==0)then
                 tmpgk13(kmap(1,i),kmap(2,i),kmap(3,i),j)=Gk(invk(1,i),j,olist(l,1),olist(m,1)) !G13(k,iw)
                 tmpgk42(kmap(1,i),kmap(2,i),kmap(3,i),j)=Gk(invk(1,i),j,olist(m,2),olist(l,2)) !G42(k,iw)
                 tmpgk13(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=conjg(Gk(invk(1,i),j,olist(m,1),olist(l,1))) !G13(k,-iw)=G^*31(k,iw)
                 tmpgk42(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=conjg(Gk(invk(1,i),j,olist(l,2),olist(m,2))) !G42(k,-iw)=G^*24(k,iw)
              else if(invk(2,i)==1)then
                 tmpgk13(kmap(1,i),kmap(2,i),kmap(3,i),j)=Gk(invk(1,i),j,olist(m,1),olist(l,1)) !G13(-k,iw)=G^31(k,iw)
                 tmpgk42(kmap(1,i),kmap(2,i),kmap(3,i),j)=Gk(invk(1,i),j,olist(l,2),olist(m,2)) !G42(-k,iw)=G^24(k,iw)
                 tmpgk13(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=conjg(Gk(invk(1,i),j,olist(l,1),olist(m,1))) !G13(-k,-iw)=G^13(k,iw)
                 tmpgk42(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=conjg(Gk(invk(1,i),j,olist(m,2),olist(l,2))) !G42(-k,-iw)=G^42(k,iw)
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
                    tmp(i,j,k,n)=-tmpgk13(i,j,k,n)*tmpgk42(ii(i),ij(j),ik(k),2*Nw-n+1)
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
              end if
           end do k_loop_tmp_to_chi
        end do w_loop_tmp_to_chi
        !$omp end parallel do
     end do l3_l4_loo2
  end do l1_l2_loop1
end subroutine get_chi0_conv

subroutine get_chi0_sum(chi,Gk,klist,olist,temp,Nw,Nk,Nkall,Norb,Nchi) bind(C)
  !
  !> It obtains chi_0 using summation. Its cost is O(Nk^2), so it is heavy. You should use get_chi0_conv.
  !
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nw,Norb,Nchi,Nk,Nkall
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

  weight=temp/dble(Nkall)
  do l=1,Nchi
     do m=1,Nchi
        !$omp parallel do private(i)
        do j=1,Nw
           do i=1,Nk
              tmpgk13(i,j)=Gk(i,j,olist(m,1),olist(l,1))
              tmpgk42(i,j)=Gk(i,j,olist(l,2),olist(m,2))
              tmpgk13(i,2*Nw-j+1)=conjg(Gk(i,j,olist(l,1),olist(m,1)))
              tmpgk42(i,2*Nw-j+1)=conjg(Gk(i,j,olist(m,2),olist(l,2)))
           end do
        end do
        !$omp end parallel do
        wloop: do iw=1,Nw
           wl_loop: do j=1,2*Nw
              wshift(j)=mod(j+iw-1,2*Nw)
              if(wshift(j)==0)then
                 wshift(j)=2*Nw
              end if
           end do wl_loop
           qloop: do iq=1,Nk
              call get_qshift(klist(:,iq),klist,qshift,Nk)
              !$omp parallel do private(i) reduction(+: chi)
              wloop2: do j=1,2*Nw
                 !$omp simd
                 kloop: do i=1,Nk
                    chi(iq,iw,m,l)=chi(iq,iw,m,l)-tmpgk13(i,j)*tmpgk42(qshift(i),wshift(j))
                 end do kloop
                 !$omp end simd
              end do wloop2
              !$omp end parallel do
           end do qloop
        end do wloop
     end do
  end do
  chi(:,:,:,:)=chi(:,:,:,:)*weight
end subroutine get_chi0_sum

subroutine get_vsigma_flex_nosoc(chi,Smat,Cmat,Nk,Nw,Nchi) bind(C,name='get_vsigma_flex_nosoc_')
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nk,Nw,Nchi
  real(real64),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  complex(real64),intent(inout),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(int32) i,j,l,m,n,info
  integer(int32),dimension(Nchi):: ipiv
  complex(real64),dimension(Nchi,Nchi):: cmat1,cmat2,cmat3,cmat4,cmat5
  complex(real64),dimension(2*Nchi):: work

  do j=1,Nw
     do i=1,Nk
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
        cmat1(:,:)=cmat5(:,:)
        cmat5(:,:)=0.0d0
        !$omp end workshare
        !$omp do private(m,n)
        do l=1,Nchi
           do m=1,Nchi
              do n=1,Nchi
                 cmat5(m,l)=cmat5(m,l)+cmat2(m,n)*cmat4(n,l) !(1-chi0S)^-1chi0S
              end do
           end do
        end do
        !$omp end do
        !$omp workshare
        cmat2(:,:)=cmat5(:,:)
        cmat5(:,:)=cmat3(:,:)+cmat4(:,:)
        cmat4(:,:)=1.5d0*Smat(:,:)-0.5d0*Cmat(:,:)
        !$omp end workshare
        !$omp do private(m,n)
        do l=1,Nchi
           do m=1,Nchi
              do n=1,Nchi
                 cmat4(m,l)=cmat4(m,l)+1.5d0*Smat(m,n)*cmat1(n,l)+0.5d0*Cmat(m,n)*cmat2(n,l)&
                      -0.25*(Cmat(m,n)+Smat(m,n))*cmat5(n,l) !subtract double count 2nd order buble(ladder)
              end do
           end do
        end do
        !$omp end do
        !$omp workshare
        chi(i,j,:,:)=cmat4(:,:)
        !$omp end workshare
        !$omp end parallel
     end do
  end do
end subroutine get_vsigma_flex_nosoc

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
        !$omp workshare
        chi(i,j,:,:)=cmat2(:,:)
        !$omp end workshare
        !$omp end parallel
     end do
  end do
end subroutine get_vsigma_flex_soc

subroutine calc_sigma(sigmak,Gk,Vsigma,Smat,Cmat,kmap,invk,olist,temp,Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz)
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz
  integer(int64),intent(in),dimension(3,Nkall):: kmap
  integer(int64),intent(in),dimension(2,Nkall):: invk
  integer(int64),intent(in),dimension(Nchi,2):: olist
  real(real64),intent(in):: temp
  real(real64),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  complex(real64),intent(in),dimension(Nk,Nw,Nchi,Nchi):: Vsigma
  complex(real64),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(real64),intent(out),dimension(Nk,Nw,Norb,Norb):: sigmak

  integer(int32) i,j,k,n,l,m
  real(real64) weight
  complex(real64),dimension(0:Nx-1,0:Ny-1,0:Nz-1,2*Nw):: tmpVsigma,tmp,tmpgk

  weight=temp/dble(Nkall)
  sigmak(:,:,:,:)=0.0d0
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
           tmpVsigma(kmap(1,i),kmap(2,i),kmap(3,i),Nw+1)=1.5d0*Smat(m,l)-0.5d0*Cmat(m,l)
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
        tmpgk=0.0d0
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
  !$Omp parallel workshare
  sigmak(:,:,:,:)=sigmak(:,:,:,:)*weight
  !$omp end parallel workshare
end subroutine calc_sigma

subroutine mkself(sigmak,Smat,Cmat,kmap,invk,olist,hamk,eig,uni,mu,rfill,temp,&
     scf_loop,pp,eps,Nkall,Nk,Nw,Norb,Nchi,Nx,Ny,Nz,sw_out,sw_in) bind(C)
  use,intrinsic:: iso_fortran_env, only:int64,real64,int32
  implicit none
  integer(int64),intent(in):: Nw,Norb,Nchi,Nkall,Nk,Nx,Ny,Nz,scf_loop
  integer(int64),intent(in),dimension(Nchi,2):: olist
  integer(int64),intent(in),dimension(3,Nkall):: kmap
  integer(int64),intent(in),dimension(2,Nkall):: invk
  logical(1),intent(in):: sw_in,sw_out
  real(real64),intent(in):: temp,eps,pp,rfill
  real(real64),intent(in),dimension(Norb,Nk):: eig
  real(real64),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  real(real64),intent(inout):: mu
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: uni,hamk
  complex(real64),intent(out),dimension(Nk,Nw,Norb,Norb):: sigmak

  integer(int32) scf_i,i
  real(real64)esterr,mu_old
  complex(real64),dimension(Nk,Nw,Norb,Norb):: Gk,sigmak0
  complex(real64),dimension(Nk,Nw,Nchi,Nchi):: chi

  if(sw_in)then
     call io_sigma(.false.)
     call gen_green_inv(Gk,sigmak,hamk,mu,temp,Nk,Nw,Norb)
     call getinv(Gk,Nk,Nw,Norb)
  else
     sigmak0(:,:,:,:)=0.0d0
     mu_old=mu*1.2
     Gk(:,:,:,:)=0.0d0 !gen_green0 need to initialization of Gk
     call gen_green0(Gk,eig,uni,mu,temp,Nk,Nw,Norb)
  end if
  do scf_i=1,scf_loop
     print*,'iter=',scf_i
     call get_chi0_conv(chi,Gk,kmap,invk,olist,temp,Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi)
     call get_Vsigma_flex_nosoc(chi,Smat,Cmat,Nk,Nw,Nchi)
     print'(A16,E12.4,A5,E12.4)','Re V_sigma: max:',maxval(dble(chi)),' min:',minval(dble(chi))
     call calc_sigma(sigmak,Gk,chi,Smat,Cmat,kmap,invk,olist,temp,Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz)
     call compair_sigma()
     if(esterr<eps)then
        exit
     end if
     call renew_mu()
     call gen_green_inv(Gk,sigmak,hamk,mu,temp,Nk,Nw,Norb)
     call getinv(Gk,Nk,Nw,Norb)
     sigmak0(:,:,:,:)=sigmak(:,:,:,:)
  end do
  if(sw_out)then
     call io_sigma(.true.)
  end if
  call renew_mu()
contains
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
                if(abs(cksigm)>1.0d-10)then
                   est=abs((sigmak0(i,j,m,l)-cksigm)/cksigm)
                   if(est>esterr)then
                      esterr=est
                      kerr=i
                      iwerr=j
                      lerr=l
                      merr=m
                   end if
                end if
                sigmak(i,j,m,l)=pp*cksigm+(1-pp)*sigmak0(i,j,m,l)
             end do
          end do
       end do
    end do
    print '(A7,E12.4,A14,2I5,2I3)','esterr=',esterr,' at ik,iw,m,l=',kerr,iwerr,merr,lerr
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
    mu_old=mu
    print'(A4,F8.4,A5,F8.4)','mu  =',mu,' rn =',rnM
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
    logical(int32),intent(in):: sw
    integer(int32)i,j,l,m
    open(55,file='sigma.bin',form='unformatted')
    if(sw)then
       write(55)mu
       write(55)mu_old
    else
       read(55)mu
       read(55)mu_old
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
end subroutine mkself
