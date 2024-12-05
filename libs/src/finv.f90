subroutine generate_irr_kpoint_inv(klist,kmap,invk_ft_list,Nk,Nx,Ny,Nz) bind(C)
  use constant
  implicit none
  integer(int64),intent(in):: Nx,Ny,Nz,Nk
  integer(int64),intent(out),dimension(2,Nx*Ny*Nz):: invk_ft_list
  integer(int64),intent(out),dimension(3,Nx*Ny*Nz):: kmap
  real(real64),intent(out),dimension(3,Nk)::klist

  real(real64),dimension(3,Nx*Ny*Nz):: all_k

  call gen_allk(all_k,kmap)
  call gen_invk(klist)
  gen_inv_ft:block
    integer(int32) Nkall,i,j,k
    real(real64) tmp(3),iktmp(3),eps
    eps=1.0d0/max(Nx,Ny,Nz)
    Nkall=Nx*Ny*Nz
    !$omp parallel do private(i,j,tmp,iktmp)
    do i=1,Nkall
       do j=1,Nk
          tmp(:)=all_k(:,i)-klist(:,j)
          if(sum(abs(tmp))<eps)then
             invk_ft_list(1,i)=j
             invk_ft_list(2,i)=0
             exit
          end if
          do k=1,3
             if(klist(k,j)==0.0d0)then
                iktmp(k)=0.0d0
             else
                iktmp(k)=1.0d0-klist(k,j)
             end if
          end do
          tmp(:)=all_k(:,i)-iktmp(:)
          if(sum(abs(tmp))<eps)then
             invk_ft_list(1,i)=j
             invk_ft_list(2,i)=1
             exit
          end if
       end do
    end do
    !$omp end parallel do
  end block gen_inv_ft
contains
  subroutine gen_allk(klist,kmap)
    integer(int64),intent(out),dimension(3,Nx*Ny*Nz):: kmap
    real(real64),intent(out),dimension(3,Nx*Ny*Nz):: klist
    integer(int32) i,j,k,iter_k
    real(real64) dx,dy,dz
    dx=1.0d0/dble(Nx)
    dy=1.0d0/dble(Ny)
    dz=1.0d0/dble(Nz)
    !$omp parallel
    !$omp workshare
    klist(:,:)=0.0d0
    !$omp end workshare
    !$omp do private(i,j,k,iter_k)
    do k=0,Nz-1
       do j=0,Ny-1
          do i=0,Nx-1
             iter_k=i+1+j*Nx+k*Nx*Ny
             kmap(1,iter_k)=i
             kmap(2,iter_k)=j
             kmap(3,iter_k)=k
             klist(1,iter_k)=i*dx
             klist(2,iter_k)=j*dy
             klist(3,iter_k)=k*dz
          end do
       end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine gen_allk

  subroutine gen_invk(klist)
    real(real64),intent(out),dimension(3,Nk)::klist

    integer(int32) i,j,k,iter_k
    klist(:,:)=0.0d0
    iter_k=1
    do j=0,int(Ny/2) !kz=0 plane
       do i=0,int(Nx/2)
          klist(1,iter_k)=dble(i)/dble(Nx)
          klist(2,iter_k)=dble(j)/dble(Ny)
          iter_k=iter_k+1
       end do
    end do
  
    do j=int(Ny/2)+1,Ny-1
       do i=1,int(Nx/2)-1
          klist(1,iter_k)=dble(i)/dble(Nx)
          klist(2,iter_k)=dble(j)/dble(Ny)
          iter_k=iter_k+1
       end do
    end do

    if(Nz>1)then
       do k=1,int(Nz/2)-1 !kz\=0,pi plane
          do j=0,Ny-1
             do i=0,Nx-1
                klist(1,iter_k)=dble(i)/dble(Nx)
                klist(2,iter_k)=dble(j)/dble(Ny)
                klist(3,iter_k)=dble(k)/dble(Nz)
                iter_k=iter_k+1
             end do
          end do
       end do

       do j=0,int(Ny/2) !kz=pi/2 plane
          do i=0,int(Nx/2)
             klist(1,iter_k)=dble(i)/dble(Nx)
             klist(2,iter_k)=dble(j)/dble(Ny)
             klist(3,iter_k)=0.5d0
             iter_k=iter_k+1
          end do
       end do

       do j=int(Ny/2)+1,Ny-1
          do i=1,int(Nx/2)-1
             klist(1,iter_k)=dble(i)/dble(Nx)
             klist(2,iter_k)=dble(j)/dble(Ny)
             klist(3,iter_k)=0.5d0
             iter_k=iter_k+1
          end do
       end do
    end if
  end subroutine gen_invk
end subroutine generate_irr_kpoint_inv
