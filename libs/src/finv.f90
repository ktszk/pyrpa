subroutine generate_irr_kpoint_inv(klist,kmap,invk_ft_list,Nk,Nx,Ny,Nz) bind(C)
  !> This function obtain irreducible klist considering reverse symmetry of k-point
  !!@param        klist,out: irreducible klist
  !!@param         kmap,out: all k properties
  !!@param invk_ft_list,out: footnote of reverse k
  !!@param            Nk,in: The number of irreducible k-point
  !!@param            Nx,in: kx mesh
  !!@param            Ny,in: ky mesh
  !!@param            Nz,in: kz mesh
  use constant
  implicit none
  integer(int64),intent(in):: Nx,Ny,Nz,Nk
  integer(int64),intent(out),dimension(3,Nx*Ny*Nz):: invk_ft_list,kmap
  real(real64),intent(out),dimension(3,Nk)::klist

  real(real64),dimension(3,Nx*Ny*Nz):: all_k

  call gen_allk(all_k,kmap)
  call gen_invk(klist)

  gen_inv_ft:block
    integer(int32) Nkall,i,j,k
    real(real64) tmp(3),iktmp(3),eps
    eps=1.0d0/max(Nx,Ny,Nz)
    Nkall=Nx*Ny*Nz
    !$omp parallel do private(i,j,k,tmp,iktmp)
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

       do j=1,Nkall
          do k=1,3
             if(all_k(k,j)==0.0d0)then
                iktmp(k)=0.0d0
             else
                iktmp(k)=1.0d0-all_k(k,j)
             end if
          end do
          tmp(:)=all_k(:,i)-iktmp(:)
          if(sum(abs(tmp))<eps)then
             invk_ft_list(3,i)=j
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

    integer(int32) i,j,k,iter_k,iter_k_ini
    !$omp parallel workshare
    klist(:,:)=0.0d0
    !$omp end parallel workshare
    if(mod(Nx,2)==0 .and. mod(Ny,2)==0)then !Nx,Ny are even
       !$omp parallel
       !$omp do private(j,i,iter_k)
       do j=0,int(Ny/2) !kz=0 plane
          do i=0,int(Nx/2)
             iter_k=i+1+int(Nx/2+1)*j
             klist(1,iter_k)=dble(i)/dble(Nx)
             klist(2,iter_k)=dble(j)/dble(Ny)
          end do
       end do
       !$omp end do
       !$omp single
       iter_k_ini=int((Nx/2+1)*(Ny/2+1))
       !$omp end single
       !$omp do private(i,j,iter_k)
       do j=int(Ny/2)+1,Ny-1
          do i=1,int(Nx/2)-1
             iter_k=iter_k_ini+i+int(Nx/2-1)*(j-int(Ny/2)-1)
             klist(1,iter_k)=dble(i)/dble(Nx)
             klist(2,iter_k)=dble(j)/dble(Ny)
          end do
       end do
       !$omp end do
       !$omp end parallel
       iter_k_ini=int(Nx*Ny/2)+3
    else
       iter_k=1
       if(mod(Nx*Ny,2)==0)then
          if(mod(Nx,2)==0)then
             do j=0,int((Ny-1)/2)
                do i=0,int(Nx/2)
                   klist(1,iter_k)=dble(i)/dble(Nx)
                   klist(2,iter_k)=dble(j)/dble(Ny)
                   iter_k=iter_k+1
                end do
             end do

             do j=int((Ny+1)/2),Ny-1
                do i=1,int(Nx/2)-1
                   klist(1,iter_k)=dble(i)/dble(Nx)
                   klist(2,iter_k)=dble(j)/dble(Ny)
                   iter_k=iter_k+1
                end do
             end do
             iter_k_ini=iter_k
          else
             continue
          end if
       else !Nx,Ny are odd
          !$omp parallel
          !$omp do private(j,i,iter_k)
          do j=0,int((Ny-1)/2) !kz=0 plane
             do i=0,int((Nx-1)/2)
                iter_k=i+1+int((Nx+1)/2)*j
                klist(1,iter_k)=dble(i)/dble(Nx)
                klist(2,iter_k)=dble(j)/dble(Ny)
             end do
          end do
          !$omp end do
          !$omp single
          iter_k_ini=int((Nx+1)*(Ny+1)/4)
          !$omp end single
          !$omp do private(i,j,iter_k)
          do j=int((Ny+1)/2),Ny-1
             do i=1,int((Nx-1)/2)
                iter_k=iter_k_ini+i+int((Nx-1)/2)*(j-int((Ny+1)/2))
                klist(1,iter_k)=dble(i)/dble(Nx)
                klist(2,iter_k)=dble(j)/dble(Ny)
             end do
          end do
          !$omp end do
          !$omp end parallel
          iter_k_ini=int((Nx*Ny+1)/2)+1
       end if
    end if

    if(Nz>1)then
       if(mod(Nz,2)==0)then !kz=pi plane (consider only Nz is even)
          !$omp parallel private(k)
          do k=1,int(Nz/2)-1 !kz\=0,pi plane
             !$omp do private(j,i,iter_k)
             do j=0,Ny-1
                do i=0,Nx-1
                   iter_k=iter_k_ini+i+Nx*j+Nx*Ny*(k-1)
                   klist(1,iter_k)=dble(i)/dble(Nx)
                   klist(2,iter_k)=dble(j)/dble(Ny)
                   klist(3,iter_k)=dble(k)/dble(Nz)
                end do
             end do
             !$omp end do
          end do
          !$omp end parallel
          iter_k_ini=iter_k_ini+Nx*Ny*int(Nz/2-1)
          if(mod(Nx,2)==0 .and. mod(Ny,2)==0)then !Nx,Ny are even
             !$omp parallel
             !$omp do private(i,j,iter_k)
             do j=0,int(Ny/2) !kz=pi/2 plane
                do i=0,int(Nx/2)
                   iter_k=iter_k_ini+i+int(Nx/2+1)*j
                   klist(1,iter_k)=dble(i)/dble(Nx)
                   klist(2,iter_k)=dble(j)/dble(Ny)
                   klist(3,iter_k)=0.5d0
                end do
             end do
             !$omp end do
             !$omp single
             iter_k_ini=iter_k_ini+int((Nx/2+1)*(Ny/2+1))-1
             !$omp end single
             !$omp do private(i,j,iter_k)
             do j=int(Ny/2)+1,Ny-1
                do i=1,int(Nx/2)-1
                   iter_k=iter_k_ini+i+int(Nx/2-1)*(j-int(Ny/2)-1)
                   klist(1,iter_k)=dble(i)/dble(Nx)
                   klist(2,iter_k)=dble(j)/dble(Ny)
                   klist(3,iter_k)=0.5d0
                end do
             end do
             !$omp end do
             !$omp end parallel
          else
             iter_k=iter_k_ini
             if(mod(Nx*Ny,2)==0)then
                if(mod(Nx,2)==0)then
                   do j=0,int((Ny-1)/2)
                      do i=0,int(Nx/2)
                         klist(1,iter_k)=dble(i)/dble(Nx)
                         klist(2,iter_k)=dble(j)/dble(Ny)
                         klist(3,iter_k)=0.5d0
                         iter_k=iter_k+1
                      end do
                   end do

                   do j=int((Ny+1)/2),Ny-1
                      do i=1,int(Nx/2)-1
                         klist(1,iter_k)=dble(i)/dble(Nx)
                         klist(2,iter_k)=dble(j)/dble(Ny)
                         klist(3,iter_k)=0.5d0
                         iter_k=iter_k+1
                      end do
                   end do
                else
                   continue
                end if
             else !Nx,Ny are odd
                do j=0,int((Ny-1)/2) !kz=0 plane
                   do i=0,int((Nx-1)/2)
                      klist(1,iter_k)=dble(i)/dble(Nx)
                      klist(2,iter_k)=dble(j)/dble(Ny)
                      klist(3,iter_k)=0.5d0
                      iter_k=iter_k+1
                   end do
                end do

                do j=int((Ny+1)/2),Ny-1
                   do i=1,int((Nx-1)/2)
                      klist(1,iter_k)=dble(i)/dble(Nx)
                      klist(2,iter_k)=dble(j)/dble(Ny)
                      klist(3,iter_k)=0.5d0
                      iter_k=iter_k+1
                   end do
                end do
             end if
          end if
       else !kz is odd
          !$omp parallel private(k)
          do k=1,int((Nz-1)/2) !kz\=0,pi plane
             !$omp do private(i,j,iter_k)
             do j=0,Ny-1
                do i=0,Nx-1
                   iter_k=iter_k_ini+i+Nx*j+Nx*Ny*(k-1)
                   klist(1,iter_k)=dble(i)/dble(Nx)
                   klist(2,iter_k)=dble(j)/dble(Ny)
                   klist(3,iter_k)=dble(k)/dble(Nz)
                end do
             end do
             !$omp end do
          end do
          !$omp end parallel
       end if !kz even or not
    end if !Nz>1
  end subroutine gen_invk
end subroutine generate_irr_kpoint_inv
