subroutine generate_irr_kpoint_inv(klist,kmap,invk_ft_list,Nk,Nx,Ny,Nz) bind(C)
  !> This function obtain irreducible klist considering reverse symmetry of k-point
  !!@param        klist,out: irreducible klist
  !!@param         kmap,out: all k properties
  !!@param invk_ft_list,out: footnote of reverse k. 1:index of irreducible k-points, 2: normal or reverse, 3: index of -k at all-ipoints
  !!@param            Nk,in: The number of irreducible k-point
  !!@param            Nx,in: kx mesh
  !!@param            Ny,in: ky mesh
  !!@param            Nz,in: kz mesh
  use constant
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nx,Ny,Nz,Nk
  integer(c_int64_t),intent(out),dimension(3,Nx*Ny*Nz):: invk_ft_list,kmap
  real(c_double),intent(out),dimension(3,Nk)::klist

  real(c_double),dimension(3,Nx*Ny*Nz):: all_k

  call gen_allk(all_k,kmap) !get info of all k-points
  call gen_irr_k(klist) !get irreducible k-points' info

  gen_inv_ft:block
    integer(c_int64_t) Nkall,i,ix,iy,iz,jx,jy,jz,jd,jt,nbad,ifl
    integer(c_int64_t),allocatable:: irr_of(:),irr_neg_of(:)

    Nkall=Nx*Ny*Nz
    allocate(irr_of(Nkall),irr_neg_of(Nkall))
    irr_of(:)=0
    irr_neg_of(:)=0

    ! Reverse lookup, built once in O(Nk).  gen_allk lays the full BZ out as the
    ! regular [0,1)^3 mesh flat = ix+1 + iy*Nx + iz*Nx*Ny, so "which full-grid
    ! point is this irreducible point" and "... is its TRS partner -k" are both
    ! plain integer arithmetic.  This replaces two linear searches -- one over
    ! the Nk irreducible points, one over all Nkall points -- that ran once per
    ! full-grid point and made the routine O(Nkall^2): measured 0.17 s at 32^3
    ! but ~10 s at 64^3 and minutes at 100^3.
    !
    ! Working in integers also retires the old coordinate comparison
    ! sum|k_i - k_j| < eps with eps = 1/max(Nx,Ny,Nz): that tolerance is exactly
    ! the smallest mesh spacing, so it was correct only by virtue of the strict
    ! < and of dble(i)/dble(N) round-tripping exactly through the multiply.
    !
    ! Serial on purpose: the loop is O(Nk) with a ten-flop body, and keeping it
    ! ordered makes "first index wins" exact if gen_irr_k ever emits the same
    ! point twice (it does for a few meshes, e.g. 1x2x1 lists Gamma twice).
    do i=1,Nk
       ix=modulo(nint(klist(1,i)*Nx,c_int64_t),Nx)
       iy=modulo(nint(klist(2,i)*Ny,c_int64_t),Ny)
       iz=modulo(nint(klist(3,i)*Nz,c_int64_t),Nz)
       ifl=ix+1+iy*Nx+iz*Nx*Ny
       if(irr_of(ifl)==0) irr_of(ifl)=i
       jx=modulo(-ix,Nx); jy=modulo(-iy,Ny); jz=modulo(-iz,Nz)
       ifl=jx+1+jy*Nx+jz*Nx*Ny
       if(irr_neg_of(ifl)==0) irr_neg_of(ifl)=i
    end do

    nbad=0
    !$omp parallel do private(i,ix,iy,iz,jx,jy,jz,jd,jt) reduction(+:nbad)
    do i=1,Nkall
       ix=kmap(1,i); iy=kmap(2,i); iz=kmap(3,i)
       jx=modulo(-ix,Nx); jy=modulo(-iy,Ny); jz=modulo(-iz,Nz)
       invk_ft_list(3,i)=jx+1+jy*Nx+jz*Nx*Ny      ! index of -k in the full list
       jd=irr_of(i)
       jt=irr_neg_of(i)
       ! The old loop scanned j upward testing k == k_j first and k == -k_j
       ! second, so the smaller index wins and a self-paired point (k == -k,
       ! i.e. Gamma and the zone-boundary points) is reported as direct.
       if(jd>0 .and. (jt==0 .or. jd<=jt))then
          invk_ft_list(1,i)=jd
          invk_ft_list(2,i)=0 !0 = direct mapping (k == irr_k)
       else if(jt>0)then
          invk_ft_list(1,i)=jt
          invk_ft_list(2,i)=1 !1 = TRS-related mapping (k == -irr_k mod 1)
       else
          invk_ft_list(1,i)=0
          invk_ft_list(2,i)=0
          nbad=nbad+1
       end if
    end do
    !$omp end parallel do

    deallocate(irr_of,irr_neg_of)
    ! Never silently hand back a zero index: it would be used as a 1-based
    ! subscript downstream and read out of bounds.  Reaching here means
    ! gen_irr_k built an incomplete wedge for this mesh (known to happen when
    ! Nx is odd and Ny is even).
    if(nbad>0)then
       write(*,'(a,i0,a)') 'generate_irr_kpoint_inv: ',nbad, &
            ' full-grid k-points have no irreducible representative'
       write(*,'(a,i0,a,i0,a,i0)') '  mesh Nx=',Nx,' Ny=',Ny,' Nz=',Nz
       write(*,'(a)') '  gen_irr_k did not cover this mesh; use a different k-mesh'
       error stop 1
    end if
  end block gen_inv_ft
contains
  subroutine gen_allk(klist,kmap)
    integer(c_int64_t),intent(out),dimension(3,Nx*Ny*Nz):: kmap
    real(c_double),intent(out),dimension(3,Nx*Ny*Nz):: klist
    integer(c_int32_t) i,j,k,iter_k
    real(c_double) dx,dy,dz
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

  subroutine gen_irr_k(klist)
    real(c_double),intent(out),dimension(3,Nk)::klist

    integer(c_int32_t) i,j,k,iter_k,iter_k_ini
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
          else !Nx is odd, Ny is even
             ! Same convention as the branch above -- fold kx into [0,Nx/2] and
             ! fold ky as well on the self-paired kx lines -- but with Nx odd
             ! there is no kx=pi, so kx=0 is the ONLY self-paired line and only
             ! that one loses its upper half of ky.  The count is the same
             ! Nx*Ny/2+1 as the Nx-even/Ny-odd case (TRS_irr.typ derives only
             ! that one), yet the construction is not its mirror image: the
             ! removal is one line's worth, not two.
             do j=0,int(Ny/2) !ky<=Ny/2: every kx of the half zone
                do i=0,int((Nx-1)/2)
                   klist(1,iter_k)=dble(i)/dble(Nx)
                   klist(2,iter_k)=dble(j)/dble(Ny)
                   iter_k=iter_k+1
                end do
             end do

             do j=int(Ny/2)+1,Ny-1 !ky>Ny/2: drop kx=0, the self-paired line
                do i=1,int((Nx-1)/2)
                   klist(1,iter_k)=dble(i)/dble(Nx)
                   klist(2,iter_k)=dble(j)/dble(Ny)
                   iter_k=iter_k+1
                end do
             end do
             iter_k_ini=iter_k
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
                else !Nx is odd, Ny is even (see the kz=0 plane above)
                   do j=0,int(Ny/2)
                      do i=0,int((Nx-1)/2)
                         klist(1,iter_k)=dble(i)/dble(Nx)
                         klist(2,iter_k)=dble(j)/dble(Ny)
                         klist(3,iter_k)=0.5d0
                         iter_k=iter_k+1
                      end do
                   end do

                   do j=int(Ny/2)+1,Ny-1
                      do i=1,int((Nx-1)/2)
                         klist(1,iter_k)=dble(i)/dble(Nx)
                         klist(2,iter_k)=dble(j)/dble(Ny)
                         klist(3,iter_k)=0.5d0
                         iter_k=iter_k+1
                      end do
                   end do
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
  end subroutine gen_irr_k
end subroutine generate_irr_kpoint_inv

subroutine get_kweight(weight,invk,Nk,Nkall) bind(C)
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nkall,Nk
  integer(c_int64_t),intent(in),dimension(3,Nkall):: invk
  real(c_double),intent(out),dimension(Nk):: weight

  integer(c_int32_t) i

  weight(:)=0.0d0
  do i=1,Nkall
     weight(invk(1,i))=weight(invk(1,i))+1.0d0
  end do
end subroutine get_kweight
