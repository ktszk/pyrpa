subroutine generate_irr_kpoint_inv(klist,invk_ft_list,Nx,Ny,Nz) bind(C)
  use constant
  implicit none
  integer(int64),intent(in):: Nx,Ny,Nz
  integer(int64),intent(out),dimension(Nx*Ny*Nz):: invk_ft_list
  real(real64),intent(out),dimension(3,int(Nx*Ny*Nz/2)+4)::klist

  real(real64),dimension(3,Nx*Ny*Nz):: all_k

  call gen_allk(all_k)
  call gen_invk(klist)
  gen_inv_ft:block
    integer(int32) Nkall,Nkinv,i,j,k
    real(real64) tmp(3),iktmp(3),eps
    eps=1.0d0/max(Nx,Ny,Nz)
    Nkall=Nx*Ny*Nz
    Nkinv=int(Nkall*0.5d0)+4
    !$omp parallel do private(i,j,tmp,iktmp)
    do i=1,Nkall
       do j=1,Nkinv
          tmp(:)=all_k(:,i)-klist(:,j)
          if(sum(abs(tmp))<eps)then
             invk_ft_list(i)=j
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
             invk_ft_list(i)=j
             exit
          end if
       end do
    end do
    !$omp end parallel do
  end block gen_inv_ft
contains
  subroutine gen_allk(klist)
    real(real64),intent(out),dimension(3,Nx*Ny*Nz):: klist
    integer(int32) i,j,k,iter_k

    !$omp parallel
    !$omp workshare
    klist(:,:)=0.0d0
    !$omp end workshare
    !$omp do private(i,j,k,iter_k)
    do k=0,Nz-1
       do j=0,Ny-1
          do i=0,Nx-1
             iter_k=i+1+j*Nx+k*Nx*Ny
             klist(1,iter_k)=dble(i)/dble(Nx)
             klist(2,iter_k)=dble(j)/dble(Ny)
             klist(3,iter_k)=dble(k)/dble(Nz)
          end do
       end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine gen_allk

  subroutine gen_invk(klist)
    real(real64),intent(out),dimension(3,int(Nx*Ny*Nz/2)+4)::klist

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
  end subroutine gen_invk
end subroutine generate_irr_kpoint_inv
