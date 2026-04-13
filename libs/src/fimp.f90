subroutine gen_imp_ham(ham_imp,ham_r,rvec,ham_i,imp_list,rlist,eps,Nimp,Nsite,Nr,Norb) bind(C)
  use,intrinsic:: iso_fortran_env, only: int64,real64
  implicit none
  integer(int64),intent(in):: Nimp,Nsite,Norb,Nr
  integer(int64),intent(in),dimension(Nimp):: imp_list
  real(real64),intent(in):: eps
  real(real64),intent(in),dimension(3,Nsite):: rlist
  real(real64),intent(in),dimension(3,Nr):: rvec
  complex(real64),intent(in),dimension(Norb,Norb,Nr):: ham_r,ham_i
  complex(real64),intent(out),dimension(Norb*Nsite,Norb*Nsite):: ham_imp

  integer(int64) i,j,k,l,m,n,k_onsite,Nrx,Nry,Nrz
  logical sw_imp1,sw_imp2
  real(real64) tmpr(3)

  Nrx=maxval(rlist(1,:))+1
  Nry=maxval(rlist(2,:))+1
  Nrz=maxval(rlist(3,:))+1
  !detect onsite_rvec
  do k=1,Nr
     if(sum(abs(rvec(:,k)))<eps)then
        k_onsite=k
        exit
     end if
  end do

  !$omp parallel do private(j,k,l,m,n,sw_imp1,sw_imp2,tmpr)
  site_loop1:do i=1,Nsite
     sw_imp1=.false.
     do n=1,Nimp
        if(i==imp_list(n)+1)then
           sw_imp1=.true.
           exit
        end if
     end do
     if(sw_imp1)then
        do l=1,Norb
           do m=l,Norb
              ham_imp(m+(i-1)*Norb,l+(i-1)*Norb)=ham_i(m,l,k_onsite)
              ham_imp(l+(i-1)*Norb,m+(i-1)*Norb)=conjg(ham_i(m,l,k_onsite))
           end do
           ham_imp(l+(i-1)*Norb,l+(i-1)*Norb)=ham_imp(l+(i-1)*Norb,l+(i-1)*Norb)
        end do
     else
        do l=1,Norb
           do m=l,Norb
              ham_imp(m+(i-1)*Norb,l+(i-1)*Norb)=ham_r(m,l,k_onsite)
              ham_imp(l+(i-1)*Norb,m+(i-1)*Norb)=conjg(ham_r(m,l,k_onsite))
           end do
           ham_imp(l+(i-1)*Norb,l+(i-1)*Norb)=ham_imp(l+(i-1)*Norb,l+(i-1)*Norb)
        end do
     end if
     site_loop2: do j=i+1,Nsite
        sw_imp2=.false.
        do n=1,Nimp
           if(j==imp_list(n))then
              sw_imp2=.true.
              exit
           end if
        end do
        tmpr(:)=rlist(:,j)-rlist(:,i) !get_distance i to j
        if(tmpr(1)>Nrx/2)tmpr(1)=tmpr(1)-Nrx
        if(tmpr(2)>Nry/2)tmpr(2)=tmpr(2)-Nry
        if(tmpr(3)>Nrz/2)tmpr(3)=tmpr(3)-Nrz
        if(tmpr(1)< -Nrx/2)tmpr(1)=tmpr(1)+Nrx
        if(tmpr(2)< -Nry/2)tmpr(2)=tmpr(2)+Nry
        if(tmpr(3)< -Nrz/2)tmpr(3)=tmpr(3)+Nrz
        !print'(3F4.1,2I2)',tmpr(:),i,j
        hop_loop:do k=1,Nr
           if(sum(abs(tmpr(:)-rvec(:,k)))<eps)then
              do l=1,Norb
                 !$omp simd
                 do m=1,Norb
                    ham_imp(m+(j-1)*Norb,l+(i-1)*Norb)=ham_r(m,l,k)
                    ham_imp(l+(i-1)*Norb,m+(j-1)*Norb)=conjg(ham_r(m,l,k))
                 end do
                 !$omp end simd
              end do
              exit
           end if
        end do hop_loop
     end do site_loop2
  end do site_loop1
  !$omp end parallel do
end subroutine gen_imp_ham

subroutine get_dft_imp_ham(ham_k,ham_imp,klist,rlist,Nk,Nsite,Norb) bind(C)
  use constant  
  implicit none
  integer(int64),intent(in):: Nk,Nsite,Norb
  real(real64),intent(in),dimension(3,Nsite):: rlist
  real(real64),intent(in),dimension(3,Nk):: klist
  complex(real64),intent(in),dimension(Norb*Nsite,Norb*Nsite):: ham_imp
  complex(real64),intent(out),dimension(Norb*Nk,Norb*Nk):: ham_k

  integer i,j,k,l,m,n
  real(real64) phase

  !$omp parallel do private(j,k,l,m,n,phase)
  k_loop1: do i=1,Nk
     k_loop2: do j=i,Nk
        site_loop1: do k=1,Nsite
           site_loop2: do l=1,Nsite
              phase=2*pi*(sum(klist(:,j)*rlist(:,l))-sum(klist(:,i)*rlist(:,k)))
              orb_loop1: do m=1,Norb
                 orb_loop2: do n=1,Norb
                    ham_k(n+Nk*(j-1),m+Nk*(i-1))=ham_k(n+Nk*(j-1),m+Nk*(i-1))&
                         +ham_imp(n+Nsite*(l-1),m+Nsite*(k-1))*cmplx(cos(phase),-sin(phase))
                 end do orb_loop2
              end do orb_loop1
           end do site_loop2
        end do site_loop1
        orbk_loop: do m=1,Norb !orb_loop for hermite setting for Hamiltonian
           orbk_loop2: do n=1,Norb
              ham_k(m+Nk*(i-1),n+Nk*(j-1))=conjg(ham_k(n+Nk*(j-1),m+Nk*(i-1)))
           end do orbk_loop2
        end do orbk_loop
     end do k_loop2
  end do k_loop1
  !$omp end parallel do
end subroutine get_dft_imp_ham

subroutine get_spectrum_spagehtti(spa,uni,eigs,klist,rlist,wlist,Nw,Nk,Nsite,Norb,mu,eta) bind(C)
  use constant
  implicit none
  integer(int64),intent(in):: Nw,Nk,Nsite,Norb
  real(real64),intent(in):: eta,mu
  real(real64),intent(in),dimension(Nw):: wlist
  real(real64),intent(in),dimension(3,Nsite):: rlist
  real(real64),intent(in),dimension(3,Nk):: klist
  real(real64),intent(in),dimension(Norb*Nsite):: eigs
  complex(real64),intent(in),dimension(Norb*Nsite,Norb*Nsite):: uni
  complex(real64),intent(out),dimension(Nw,Nk):: spa

  integer i,j,k,l,m,n
  real(real64) phase

  !$omp parallel
  !$omp workshare
  spa(:,:)=0.0d0
  !$omp end workshare
  !$omp do private(i,j,k,l,n,m,phase)
  k_loop: do i=1,Nk
     site_loop1: do j=1,Nsite
        site_loop2: do k=1,Nsite
           phase=2*pi*sum(klist(:,i)*(rlist(:,j)-rlist(:,k)))
           orb_loop: do l=1,Norb
              eig_loop: do n=1,Nsite*Norb
                 w_loop:do m=1,Nw
                    spa(m,i)=spa(m,i)+conjg(uni(l+(j-1)*Nsite,n))*uni(l+(k-1)*Nsite,n)&
                         *cmplx(cos(phase),-sin(phase))/cmplx(wlist(m)-eigs(n)+mu,eta)
                 end do w_loop
              end do eig_loop
           end do orb_loop
        end do site_loop2
     end do site_loop1
  end do k_loop
  !$omp end do
  !$omp end parallel
end subroutine get_spectrum_spagehtti

subroutine calc_gloc(Gloc,hamk,sigma,z,Nk,Norb) bind(C)
  !> calc_gloc
  !> Compute local Green's function: Gloc = (1/Nk) sum_k [z*I - H(k) - Sigma]^(-1)
  !!@param   Gloc,out: local Green's function [Norb,Norb] complex128
  !!@param   hamk, in: k-space Hamiltonian [Norb,Norb,Nk] complex128
  !!@param  sigma, in: self-energy [Norb,Norb] complex128
  !!@param      z, in: complex frequency (iω_n or ω+iδ)
  !!@param     Nk, in: number of k-points
  !!@param   Norb, in: number of orbitals
  use,intrinsic:: iso_fortran_env, only:int32,int64,real64
  use omp_lib, only: omp_get_active_level
  implicit none
  integer(int64),intent(in):: Nk,Norb
  complex(real64),intent(in):: z
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: hamk
  complex(real64),intent(in),dimension(Norb,Norb):: sigma
  complex(real64),intent(out),dimension(Norb,Norb):: Gloc

  complex(real64),dimension(Norb,Norb):: gk
  complex(real64),dimension(Norb*Norb):: work
  integer(int32),dimension(Norb):: ipiv
  integer(int64) :: ik
  integer(int32) :: i,j,info

  Gloc(:,:) = (0.0d0, 0.0d0)

  !$omp parallel if(omp_get_active_level()==0) private(gk,ipiv,work,info,i,j)
  !$omp do reduction(+:Gloc)
  do ik = 1, Nk
     ! gk = z*I - H(k) - Sigma
     do j = 1, Norb
        do i = 1, Norb
           gk(i,j) = -hamk(i,j,ik) - sigma(i,j)
        end do
        gk(j,j) = gk(j,j) + z
     end do
     ! invert gk
     call zgetrf(Norb, Norb, gk, Norb, ipiv, info)
     if (info /= 0) cycle
     call zgetri(Norb, gk, Norb, ipiv, work, Norb*Norb, info)
     if (info /= 0) cycle
     ! accumulate
     do j = 1, Norb
        do i = 1, Norb
           Gloc(i,j) = Gloc(i,j) + gk(i,j)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

  Gloc(:,:) = Gloc(:,:) / dble(Nk)
end subroutine calc_gloc

subroutine calc_tmat_cpa(sigma_new,Gloc,sigma,VA,VB,x,Norb) bind(C)
  !> calc_tmat_cpa
  !> Compute CPA t-matrices and update self-energy:
  !>   delta_A = V_A - Sigma, delta_B = V_B - Sigma
  !>   t_A = delta_A [I - Gloc*delta_A]^(-1)
  !>   t_B = delta_B [I - Gloc*delta_B]^(-1)
  !>   Sigma_new = Sigma + x*t_A + (1-x)*t_B
  !!@param sigma_new,out: updated self-energy [Norb,Norb]
  !!@param     Gloc,  in: local Green's function [Norb,Norb]
  !!@param    sigma,  in: current self-energy [Norb,Norb]
  !!@param       VA,  in: onsite potential of species A [Norb,Norb]
  !!@param       VB,  in: onsite potential of species B [Norb,Norb]
  !!@param        x,  in: concentration of species A
  !!@param     Norb,  in: number of orbitals
  use,intrinsic:: iso_fortran_env, only:int32,int64,real64
  implicit none
  integer(int64),intent(in):: Norb
  real(real64),intent(in):: x
  complex(real64),intent(in),dimension(Norb,Norb):: Gloc,sigma,VA,VB
  complex(real64),intent(out),dimension(Norb,Norb):: sigma_new

  complex(real64),dimension(Norb,Norb):: dA,dB,matA,matB,tA,tB
  complex(real64),dimension(Norb*Norb):: work
  integer(int32),dimension(Norb):: ipiv
  integer(int32) :: i,j,info

  ! delta_alpha = V_alpha - Sigma
  dA(:,:) = VA(:,:) - sigma(:,:)
  dB(:,:) = VB(:,:) - sigma(:,:)

  ! matA = I - Gloc * delta_A
  call zgemm('N','N',Norb,Norb,Norb,(-1.0d0,0.0d0),Gloc,Norb,dA,Norb,(0.0d0,0.0d0),matA,Norb)
  do i = 1, Norb
     matA(i,i) = matA(i,i) + (1.0d0, 0.0d0)
  end do
  ! invert matA
  call zgetrf(Norb, Norb, matA, Norb, ipiv, info)
  if (info /= 0) then; print*, 'CPA: zgetrf(A) failed, info=', info; stop; end if
  call zgetri(Norb, matA, Norb, ipiv, work, Norb*Norb, info)
  if (info /= 0) then; print*, 'CPA: zgetri(A) failed, info=', info; stop; end if
  ! tA = delta_A * matA^(-1)
  call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),dA,Norb,matA,Norb,(0.0d0,0.0d0),tA,Norb)

  ! matB = I - Gloc * delta_B
  call zgemm('N','N',Norb,Norb,Norb,(-1.0d0,0.0d0),Gloc,Norb,dB,Norb,(0.0d0,0.0d0),matB,Norb)
  do i = 1, Norb
     matB(i,i) = matB(i,i) + (1.0d0, 0.0d0)
  end do
  ! invert matB
  call zgetrf(Norb, Norb, matB, Norb, ipiv, info)
  if (info /= 0) then; print*, 'CPA: zgetrf(B) failed, info=', info; stop; end if
  call zgetri(Norb, matB, Norb, ipiv, work, Norb*Norb, info)
  if (info /= 0) then; print*, 'CPA: zgetri(B) failed, info=', info; stop; end if
  ! tB = delta_B * matB^(-1)
  call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),dB,Norb,matB,Norb,(0.0d0,0.0d0),tB,Norb)

  ! sigma_new = sigma + x*tA + (1-x)*tB
  do j = 1, Norb
     do i = 1, Norb
        sigma_new(i,j) = sigma(i,j) + x*tA(i,j) + (1.0d0-x)*tB(i,j)
     end do
  end do
end subroutine calc_tmat_cpa

subroutine solve_cpa(sigma_cpa,hamk,VA,VB,x,z,pp,Nk,Norb,maxiter,tol) bind(C)
  !> solve_cpa
  !> Self-consistent CPA loop for a single frequency z.
  !>   1. G_loc = (1/Nk) sum_k [z - H(k) - Sigma_CPA]^(-1)
  !>   2. t_A, t_B from V_A, V_B and G_loc
  !>   3. Sigma_CPA <- Sigma_CPA + x*t_A + (1-x)*t_B  (with linear mixing)
  !>   4. Repeat until convergence.
  !!@param sigma_cpa,inout: CPA self-energy [Norb,Norb] — supply initial guess, returns converged value
  !!@param      hamk,   in: k-space Hamiltonian [Norb,Norb,Nk]
  !!@param        VA,   in: onsite potential of species A [Norb,Norb]
  !!@param        VB,   in: onsite potential of species B [Norb,Norb]
  !!@param         x,   in: concentration of species A
  !!@param         z,   in: complex frequency
  !!@param        pp,   in: linear mixing parameter (0<pp<=1)
  !!@param        Nk,   in: number of k-points
  !!@param      Norb,   in: number of orbitals
  !!@param   maxiter,   in: maximum iterations
  !!@param       tol,   in: convergence tolerance
  use,intrinsic:: iso_fortran_env, only:int64,real64
  implicit none
  integer(int64),intent(in):: Nk,Norb,maxiter
  real(real64),intent(in):: x,pp,tol
  complex(real64),intent(in):: z
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: hamk
  complex(real64),intent(in),dimension(Norb,Norb):: VA,VB
  complex(real64),intent(inout),dimension(Norb,Norb):: sigma_cpa

  complex(real64),dimension(Norb,Norb):: Gloc,sigma_new
  real(real64) diff
  integer(int64) iter,i,j

  interface
     subroutine calc_gloc(Gloc,hamk,sigma,z,Nk,Norb) bind(C)
       import:: int64,real64
       integer(int64),intent(in):: Nk,Norb
       complex(real64),intent(in):: z
       complex(real64),intent(in),dimension(Norb,Norb,Nk):: hamk
       complex(real64),intent(in),dimension(Norb,Norb):: sigma
       complex(real64),intent(out),dimension(Norb,Norb):: Gloc
     end subroutine calc_gloc
     subroutine calc_tmat_cpa(sigma_new,Gloc,sigma,VA,VB,x,Norb) bind(C)
       import:: int64,real64
       integer(int64),intent(in):: Norb
       real(real64),intent(in):: x
       complex(real64),intent(in),dimension(Norb,Norb):: Gloc,sigma,VA,VB
       complex(real64),intent(out),dimension(Norb,Norb):: sigma_new
     end subroutine calc_tmat_cpa
  end interface

  do iter = 1, maxiter
     ! step 1: compute local Green's function
     call calc_gloc(Gloc, hamk, sigma_cpa, z, Nk, Norb)

     ! step 2-3: compute t-matrices and updated self-energy
     call calc_tmat_cpa(sigma_new, Gloc, sigma_cpa, VA, VB, x, Norb)

     ! convergence check
     diff = 0.0d0
     do j = 1, Norb
        do i = 1, Norb
           diff = diff + abs(sigma_new(i,j) - sigma_cpa(i,j))**2
        end do
     end do
     diff = sqrt(diff)

     ! linear mixing
     do j = 1, Norb
        do i = 1, Norb
           sigma_cpa(i,j) = (1.0d0-pp)*sigma_cpa(i,j) + pp*sigma_new(i,j)
        end do
     end do

     if (diff < tol) exit
  end do
end subroutine solve_cpa

subroutine solve_cpa_array(sigma_cpa_w,hamk,VA,VB,x,zlist,pp,Nk,Norb,Nw,maxiter,tol) bind(C)
  !> solve_cpa_array
  !> Run CPA self-consistent loop for an array of frequencies (松原 or 実軸).
  !> Each frequency is solved independently — parallelized over ω.
  !!@param sigma_cpa_w,inout: CPA self-energy [Norb,Norb,Nw] — initial guess in, converged out
  !!@param        hamk,   in: k-space Hamiltonian [Norb,Norb,Nk]
  !!@param          VA,   in: onsite of A [Norb,Norb]
  !!@param          VB,   in: onsite of B [Norb,Norb]
  !!@param           x,   in: concentration of A
  !!@param       zlist,   in: complex frequency array [Nw]
  !!@param          pp,   in: mixing parameter
  !!@param          Nk,   in: number of k-points
  !!@param        Norb,   in: number of orbitals
  !!@param          Nw,   in: number of frequencies
  !!@param     maxiter,   in: max iterations per frequency
  !!@param         tol,   in: convergence tolerance
  use,intrinsic:: iso_fortran_env, only:int64,real64
  implicit none
  integer(int64),intent(in):: Nk,Norb,Nw,maxiter
  real(real64),intent(in):: x,pp,tol
  complex(real64),intent(in),dimension(Nw):: zlist
  complex(real64),intent(in),dimension(Norb,Norb,Nk):: hamk
  complex(real64),intent(in),dimension(Norb,Norb):: VA,VB
  complex(real64),intent(inout),dimension(Norb,Norb,Nw):: sigma_cpa_w

  integer(int64) iw

  interface
     subroutine solve_cpa(sigma_cpa,hamk,VA,VB,x,z,pp,Nk,Norb,maxiter,tol) bind(C)
       import:: int64,real64
       integer(int64),intent(in):: Nk,Norb,maxiter
       real(real64),intent(in):: x,pp,tol
       complex(real64),intent(in):: z
       complex(real64),intent(in),dimension(Norb,Norb,Nk):: hamk
       complex(real64),intent(in),dimension(Norb,Norb):: VA,VB
       complex(real64),intent(inout),dimension(Norb,Norb):: sigma_cpa
     end subroutine solve_cpa
  end interface

  !$omp parallel do private(iw)
  do iw = 1, Nw
     call solve_cpa(sigma_cpa_w(:,:,iw), hamk, VA, VB, x, zlist(iw), pp, Nk, Norb, maxiter, tol)
  end do
  !$omp end parallel do
end subroutine solve_cpa_array
