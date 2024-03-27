subroutine FFT(cmat,tmp,Nx,Ny,Nz,Nw,SW)
    implicit none
    integer(8),intent(in):: Nx,Ny,Nz,Nw
    logical(4),intent(in):: SW
    complex(8),intent(inout),dimension(Nx,Ny,Nz,Nw):: cmat,tmp
  
    integer(8) plan
    integer(4) Inv
    integer(4),dimension(4):: Nlist

    Nlist=(/Nx,Ny,Nx,Nw/)
    if(SW)then
       Inv=-1
    else
       Inv=1
    end If
    call dfftw_plan_dft(plan,4,Nlist,cmat,tmp,Inv,64)
    !call dfftw_execute(plan)
    !call dfftw_destroy_plan(plan)
    !if(.not. SW) cmat=tmp/product(Nlist)
end subroutine FFT

subroutine gen_green0(Gk,eig,uni,mu,temp,Nk,Nw,norb) bind(C)
   use constant
   implicit none
   integer(8),intent(in):: Nk,Nw,norb
   real(8),intent(in):: mu,temp
   real(8),intent(in),dimension(norb,Nk):: eig
   complex(8),intent(in),dimension(norb,norb,Nk):: uni
   complex(8),intent(out),dimension(Nk,Nw,norb,norb):: Gk
   
   integer(8) i,j,l,m,n
   complex(8) iw
 
   do l=1,norb
      do m=1,norb
         band_loop: do n=1,norb
            !$omp parallel do private(iw,i,j)
            wloop: do j=1,Nw
               iw=cmplx(mu,dble(2*(j-1)+1)*pi*temp)
               !!$omp simd
               kloop: do i=1,Nk
                  Gk(i,j,m,l)=Gk(i,j,m,l)+uni(l,n,i)*conjg(uni(m,n,i))/(iw-eig(n,i))
               end do kloop
               !!$omp end simd
            end do wloop
            !$omp end parallel do
         end do band_loop
      end do
   end do
 end subroutine gen_green0

subroutine get_chi0_comb(chi,Gk,kmap,olist,Nx,Ny,Nz,Nw,Nk,Norb,Nchi) bind(C)
  implicit none
  integer(8),intent(in):: Nw,Norb,Nchi,Nk,Nx,Ny,Nz
  integer(8),intent(in),dimension(Nchi,2):: olist
  integer(8),intent(in),dimension(3,Nk):: kmap
  complex(8),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(8),intent(out),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(8) i,j,k,l,m,n
  integer(8) ii(0:Nx-1),ij(0:Ny-1),ik(0:Nz-1)
  complex(8),dimension(0:Nx-1,0:Ny-1,0:Nz-1,2*Nw):: tmpchi,tmp,tmp1,tmp2

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
  do l=1,Nchi
     do m=1,Nchi
        !use symmetry G^lm(k,iw)=G^ml(k,-iw)
        !$omp parallel do private(i)
        do j=1,Nw
           do i=1,Nk
              tmp1(kmap(1,i),kmap(2,i),kmap(3,i),j)=Gk(i,j,olist(l,1),olist(m,2)) !G13(iw)
              tmp2(kmap(1,i),kmap(2,i),kmap(3,i),j)=Gk(i,j,olist(l,2),olist(m,1)) !G42(iw)
              tmp1(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=conjg(Gk(i,j,olist(m,2),olist(l,1))) !G13(-iw)
              tmp2(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=conjg(Gk(i,j,olist(m,1),olist(l,2))) !G42(-iw)
           end do
        end do
        !$omp end parallel do
        call FFT(tmp1,tmp,Nk,1,1,Nw,.true.)
        call FFT(tmp2,tmp,Nk,1,1,Nw,.true.)
        !calculate G(r)G(-r)
        !$omp parallel do private(i,j,k)
        do n=1,2*Nw
           do k=0,Nz-1
              do j=0,Ny-1
                 !$omp simd
                 do i=0,Nx-1
                    tmpchi(i,j,k,n)=tmp1(i,j,k,n)*tmp2(ii(i),ij(j),ik(k),2*nw-n+1)
                 end do
                 !$omp end simd
              end do
           end do
        end do
        !$omp end parallel do
        call FFT(tmpchi,tmp,Nx,Ny,Nz,2*Nw,.false.)
        !$omp parallel do private(i,j)
        do j=1,Nw
           do i=1,Nk
              chi(i,j,m,l)=tmpchi(kmap(1,i),kmap(2,i),kmap(3,i),j)
           end do
        end do
        !$omp end parallel do
     end do
  end do
end subroutine get_chi0_comb
