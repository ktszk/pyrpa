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
    call dfftw_execute(plan)
    call dfftw_destroy_plan(plan)
    if(.not. SW) cmat=tmp/product(Nlist)
end subroutine FFT

subroutine get_chi0_comb(chi,Gk,kmap,olist,Nx,Ny,Nz,Nw,Nk,Norb,Nchi) bind(C)
  implicit none
  integer(8),intent(in):: Nx,Ny,Nz,Nw,Norb,Nchi,Nk
  integer(8),intent(in),dimension(Nchi,2):: olist
  integer(8),intent(in),dimension(3,Nk):: kmap
  complex(8),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(8),intent(out),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(8) i,j,k,l,m,n
  integer(8) ii(0:Nx-1),ij(0:Ny-1),ik(0:Nz-1)
  complex(8),dimension(0:Nx-1,0:Ny-1,0:Nz-1,2*Nw):: tmp1,tmp2,tmpchi,tmp

  do i=0,Nx-1
     ii(i)=mod(Nx-i,Nx)
  end do
  
  do j=0,Ny-1
     ij(j)=mod(Ny-j,Ny)
  end do
  
  do k=0,Nz-1
     ik(k)=mod(Nz-k,Nz)
  end do

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
        call FFT(tmp1,tmp,Nx,Ny,Nz,2*Nw,.true.)
        call FFT(tmp2,tmp,Nx,Ny,Nz,2*Nw,.true.)
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