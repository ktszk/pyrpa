! =============================================================================
! fimp.f90  —  Impurity and CPA (Coherent Potential Approximation) routines
! =============================================================================
!
! INDEX CONVENTION (used throughout this file)
! ---------------------------------------------
! All packed matrices follow orbital-fast / site(or k-point)-slow ordering:
!
!   ham_imp(Norb*Nsite, Norb*Nsite)
!     row/col for orbital m at site i  →  m + (i-1)*Norb   (m=1..Norb, i=1..Nsite)
!
!   ham_k(Norb*Nk, Norb*Nk)
!     row/col for orbital m at k-point i  →  m + (i-1)*Norb  (m=1..Norb, i=1..Nk)
!
!   uni(Norb*Nsite, Norb*Nsite)   — eigenvectors of ham_imp
!     same packing as ham_imp: orbital m at site i  →  m + (i-1)*Norb
!
!   sigma_cpa_w(Norb, Norb, Nw)  — CPA self-energy per frequency point
!
! IMP_LIST CONVENTION
! -------------------
!   imp_list uses 0-based indices (Python/C convention).
!   Inside the Fortran loops the conversion is:  Fortran_index = imp_list(n) + 1
!
! WANNIER HAMILTONIAN CONVENTION (used in gen_imp_ham)
! -----------------------------------------------------
! Following the Wannier90 convention, all H(R) arrays are defined as:
!
!   H(m, l, k)  =  <m, rvec(:,k) | H | l, 0>
!
! where m and l are orbital indices and rvec(:,k) is the lattice displacement.
! This means for a supercell pair (site i, site j) with displacement
! tmpr = R_j - R_i, the Hamiltonian block (row=m at j, col=l at i) is:
!
!   <m, R_j | H | l, R_i>  =  <m, tmpr | H | l, 0>  =  H(m, l, k) at rvec(:,k)=tmpr
!
! Three H(R) arrays cover the four possible species combinations:
!
!   ham_r(m,l,k):  host-host    H_rr(R),  <m_host, R | H | l_host, 0>
!   ham_i(m,l,k):  imp-imp      H_ii(R),  <m_imp,  R | H | l_imp,  0>  (R=0 → on-site)
!   ham_ri(m,l,k): cross-type   H_ri(R),  <m_imp,  R | H | l_host, 0>
!
! Hopping for pair (i<j), tmpr = R_j - R_i:
!
!   host_i, host_j  →  ham_r(m,l,k)          at rvec = +tmpr
!   imp_i,  imp_j   →  ham_i(m,l,k)          at rvec = +tmpr
!   host_i, imp_j   →  ham_ri(m,l,k)         at rvec = +tmpr
!   imp_i,  host_j  →  conjg(ham_ri(l,m,k))  at rvec = -tmpr
!                       (derived: <m_host,tmpr|H|l_imp,0>
!                               = conjg(<l_imp,-tmpr|H|m_host,0>)
!                               = conjg(ham_ri(l,m,k)) at rvec=-tmpr)
!
! The lower triangle is always set as the Hermitian conjugate of the upper.
!
! =============================================================================

! -----------------------------------------------------------------------------
subroutine gen_imp_ham(ham_imp,ham_r,rvec,ham_i,ham_ri,imp_list,rlist,eps,Nimp,Nsite,Nr,Norb) bind(C)
!> gen_imp_ham
!> Build the real-space supercell Hamiltonian for a host crystal with embedded
!> impurities.  All three Wannier Hamiltonians (host-host, imp-imp, cross) are
!> used to fill both on-site and hopping blocks.
!!@param  ham_imp,out: supercell Hamiltonian [Norb*Nsite,Norb*Nsite] complex128
!!@param    ham_r, in: host-host Wannier H_rr(R) [Norb,Norb,Nr] complex128
!!@param     rvec, in: Wannier R-vectors in fractional coords [3,Nr]
!!@param    ham_i, in: imp-imp Wannier H_ii(R) [Norb,Norb,Nr] complex128; R=0 → on-site
!!@param   ham_ri, in: cross-species Wannier H_ri(R) [Norb,Norb,Nr] complex128
!!@param imp_list, in: 0-based indices of impurity sites [Nimp] c_int64_t
!!@param    rlist, in: site positions in fractional coords [3,Nsite]
!!@param      eps, in: tolerance for R-vector matching
!!@param     Nimp, in: number of impurity sites
!!@param    Nsite, in: total number of sites in supercell
!!@param       Nr, in: number of Wannier R-vectors
!!@param     Norb, in: number of orbitals per site
!
! Wannier convention: H(m,l,k) = <m, rvec(:,k) | H | l, 0>
!   ham_r(m,l,k)  = <m_host, R | H | l_host, 0>
!   ham_i(m,l,k)  = <m_imp,  R | H | l_imp,  0>  (R=0 → on-site)
!   ham_ri(m,l,k) = <m_imp,  R | H | l_host, 0>  (bra=imp, ket=host)
!
! Hopping selection for pair (i<j), tmpr = R_j - R_i:
!
!   host_i, host_j  →  ham_r(m,l,k)          at rvec = +tmpr
!   imp_i,  imp_j   →  ham_i(m,l,k)          at rvec = +tmpr
!   host_i, imp_j   →  ham_ri(m,l,k)         at rvec = +tmpr
!   imp_i,  host_j  →  conjg(ham_ri(l,m,k))  at rvec = -tmpr
!                       (by Hermitian symmetry:
!                        <m_host,tmpr|H|l_imp,0>
!                          = conjg(<l_imp,-tmpr|H|m_host,0>)
!                          = conjg(ham_ri(l,m,k)) at rvec=-tmpr)
!
! The lower triangle is always set as the Hermitian conjugate of the upper.
! If a required R-vector is not found in the table, that pair is skipped.
! -----------------------------------------------------------------------------
  use,intrinsic:: iso_c_binding, only: c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nimp,Nsite,Norb,Nr
  integer(c_int64_t),intent(in),dimension(Nimp):: imp_list  ! 0-based impurity site indices
  real(c_double),intent(in):: eps                          ! tolerance for R-vector matching
  real(c_double),intent(in),dimension(3,Nsite):: rlist     ! fractional coords of each site
  real(c_double),intent(in),dimension(3,Nr):: rvec         ! Wannier R-vectors (fractional)
  complex(c_double),intent(in),dimension(Norb,Norb,Nr):: ham_r   ! host-host H(R)
  complex(c_double),intent(in),dimension(Norb,Norb,Nr):: ham_i   ! imp-imp   H(R)
  complex(c_double),intent(in),dimension(Norb,Norb,Nr):: ham_ri  ! cross     H(R), bra=imp
  complex(c_double),intent(out),dimension(Norb*Nsite,Norb*Nsite):: ham_imp

  integer(c_int64_t) i,j,k,l,m,n,k_onsite,k_hop,k_hop_neg,Nrx,Nry,Nrz
  logical sw_imp1,sw_imp2
  real(c_double) tmpr(3)

  ! Supercell dimensions (used for minimum-image wrapping below)
  Nrx=maxval(rlist(1,:))+1
  Nry=maxval(rlist(2,:))+1
  Nrz=maxval(rlist(3,:))+1

  ! Locate the on-site R-vector (R=0) in the rvec table
  k_onsite=0
  do k=1,Nr
     if(sum(abs(rvec(:,k)))<eps)then
        k_onsite=k
        exit
     end if
  end do
  if(k_onsite==0)then
     print*,'Error: no onsite R-vector (R=0) found in gen_imp_ham'
     stop
  end if

  ! Parallelise over sites.  All scalar temporaries are thread-private.
  !$omp parallel do private(j,k,l,m,n,sw_imp1,sw_imp2,tmpr,k_hop,k_hop_neg)
  site_loop1:do i=1,Nsite

     ! --- Determine whether site i is an impurity (imp_list is 0-based) ---
     sw_imp1=.false.
     do n=1,Nimp
        if(i==imp_list(n)+1)then
           sw_imp1=.true.
           exit
        end if
     end do

     ! --- Fill the on-site (diagonal) block for site i ---
     ! Only the upper triangle m>=l is written explicitly; the lower triangle
     ! is set as its Hermitian conjugate.
     ! Index: orbital m at site i  →  m + (i-1)*Norb
     if(sw_imp1)then
        ! Impurity on-site: ham_i(R=0)
        do l=1,Norb
           do m=l,Norb
              ham_imp(m+(i-1)*Norb,l+(i-1)*Norb)=ham_i(m,l,k_onsite)
              ham_imp(l+(i-1)*Norb,m+(i-1)*Norb)=conjg(ham_i(m,l,k_onsite))
           end do
        end do
     else
        ! Host on-site: ham_r(R=0)
        do l=1,Norb
           do m=l,Norb
              ham_imp(m+(i-1)*Norb,l+(i-1)*Norb)=ham_r(m,l,k_onsite)
              ham_imp(l+(i-1)*Norb,m+(i-1)*Norb)=conjg(ham_r(m,l,k_onsite))
           end do
        end do
     end if

     ! --- Fill off-diagonal (hopping) blocks between site i and site j>i ---
     site_loop2: do j=i+1,Nsite

        ! Determine whether site j is an impurity
        sw_imp2=.false.
        do n=1,Nimp
           if(j==imp_list(n)+1)then
              sw_imp2=.true.
              exit
           end if
        end do

        ! Displacement tmpr = R_j - R_i in fractional coordinates.
        ! Apply minimum-image convention so that tmpr falls within
        ! [-Nrα/2, +Nrα/2] along each axis.
        tmpr(:)=rlist(:,j)-rlist(:,i)
        if(tmpr(1)>Nrx/2)tmpr(1)=tmpr(1)-Nrx
        if(tmpr(2)>Nry/2)tmpr(2)=tmpr(2)-Nry
        if(tmpr(3)>Nrz/2)tmpr(3)=tmpr(3)-Nrz
        if(tmpr(1)< -Nrx/2)tmpr(1)=tmpr(1)+Nrx
        if(tmpr(2)< -Nry/2)tmpr(2)=tmpr(2)+Nry
        if(tmpr(3)< -Nrz/2)tmpr(3)=tmpr(3)+Nrz

        ! Look up R-vector index for +tmpr (used by all cases except imp_i→host_j)
        k_hop=0
        do k=1,Nr
           if(sum(abs(tmpr(:)-rvec(:,k)))<eps)then
              k_hop=k
              exit
           end if
        end do

        ! For (imp_i, host_j): additionally look up index for -tmpr.
        ! Required because the Hermitian identity pulls out ham_ri at -tmpr;
        ! see the hopping table in the file header for the derivation.
        k_hop_neg=0
        if(sw_imp1 .and. .not.sw_imp2)then
           do k=1,Nr
              if(sum(abs(tmpr(:)+rvec(:,k)))<eps)then  ! rvec(:,k) == -tmpr
                 k_hop_neg=k
                 exit
              end if
           end do
        end if

        ! --- Fill hopping block and its Hermitian conjugate ---
        ! Upper triangle element: ham_imp(m+(j-1)*Norb, l+(i-1)*Norb)
        !   = <m_j, tmpr | H | l_i, 0>   (Wannier convention, see file header)
        ! Lower triangle: Hermitian conjugate, always set explicitly here.
        if(.not.sw_imp1 .and. .not.sw_imp2)then
           ! host_i, host_j  →  ham_r at +tmpr
           if(k_hop==0) cycle site_loop2
           do l=1,Norb
              !$omp simd
              do m=1,Norb
                 ham_imp(m+(j-1)*Norb,l+(i-1)*Norb)=ham_r(m,l,k_hop)
                 ham_imp(l+(i-1)*Norb,m+(j-1)*Norb)=conjg(ham_r(m,l,k_hop))
              end do
              !$omp end simd
           end do

        else if(sw_imp1 .and. sw_imp2)then
           ! imp_i, imp_j  →  ham_i at +tmpr
           if(k_hop==0) cycle site_loop2
           do l=1,Norb
              !$omp simd
              do m=1,Norb
                 ham_imp(m+(j-1)*Norb,l+(i-1)*Norb)=ham_i(m,l,k_hop)
                 ham_imp(l+(i-1)*Norb,m+(j-1)*Norb)=conjg(ham_i(m,l,k_hop))
              end do
              !$omp end simd
           end do

        else if(.not.sw_imp1 .and. sw_imp2)then
           ! host_i, imp_j  →  ham_ri(m,l,k) at +tmpr
           ! ham_ri(m,l,k) = <m_imp, rvec(:,k)|H|l_host,0>, so direct use at rvec=+tmpr
           if(k_hop==0) cycle site_loop2
           do l=1,Norb
              !$omp simd
              do m=1,Norb
                 ham_imp(m+(j-1)*Norb,l+(i-1)*Norb)=ham_ri(m,l,k_hop)
                 ham_imp(l+(i-1)*Norb,m+(j-1)*Norb)=conjg(ham_ri(m,l,k_hop))
              end do
              !$omp end simd
           end do

        else
           ! imp_i, host_j  →  conjg(ham_ri(l,m,k)) at rvec=-tmpr
           ! Derivation: <m_host,tmpr|H|l_imp,0>
           !           = conjg(<l_imp,-tmpr|H|m_host,0>)  [H Hermitian + translation]
           !           = conjg(ham_ri(l,m,k_hop_neg))     [note: l,m index swap]
           if(k_hop_neg==0) cycle site_loop2
           do l=1,Norb
              !$omp simd
              do m=1,Norb
                 ham_imp(m+(j-1)*Norb,l+(i-1)*Norb)=conjg(ham_ri(l,m,k_hop_neg))
                 ham_imp(l+(i-1)*Norb,m+(j-1)*Norb)=ham_ri(l,m,k_hop_neg)
              end do
              !$omp end simd
           end do
        end if

     end do site_loop2
  end do site_loop1
  !$omp end parallel do
end subroutine gen_imp_ham

! -----------------------------------------------------------------------------
subroutine gen_bdg_ham(ham_bdg,ham_imp,gap_imp,Nsite,Norb) bind(C)
!> gen_bdg_ham
!> Build the real-space Bogoliubov-de Gennes (BdG) Hamiltonian from a
!> normal-state supercell Hamiltonian and a gap (pairing) function, both
!> produced by gen_imp_ham (or the same packing convention).
!>
!> Nambu basis ordering (orbital-fast, site-slow):
!>   particle sector: orbital m at site i  →  m + (i-1)*Norb         (row/col 1..Norb*Nsite)
!>   hole    sector: orbital m at site i  →  m + (i-1)*Norb + Norb*Nsite (row/col Norb*Nsite+1..2*Norb*Nsite)
!>
!> Block structure of the 2*Norb*Nsite × 2*Norb*Nsite matrix:
!>
!>   ham_bdg = |  ham_imp          gap_imp        |
!>             |  gap_imp†       - conjg(ham_imp)  |
!>
!!@param  ham_bdg, out: BdG Hamiltonian [2*Norb*Nsite, 2*Norb*Nsite] complex128
!!@param  ham_imp,  in: normal-state supercell H [Norb*Nsite, Norb*Nsite] complex128
!!@param  gap_imp,  in: real-space gap function  [Norb*Nsite, Norb*Nsite] complex128
!!@param    Nsite,  in: number of sites in supercell
!!@param     Norb,  in: number of orbitals per site
! -----------------------------------------------------------------------------
  use,intrinsic:: iso_c_binding, only: c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nsite,Norb
  complex(c_double),intent(in), dimension(Norb*Nsite,Norb*Nsite):: ham_imp
  complex(c_double),intent(in), dimension(Norb*Nsite,Norb*Nsite):: gap_imp
  complex(c_double),intent(out),dimension(2*Norb*Nsite,2*Norb*Nsite):: ham_bdg

  integer(c_int64_t) N

  N = Norb*Nsite
  ham_bdg(:,:) = (0.0d0,0.0d0)
  ! Top-left block:  H
  ham_bdg(1:N, 1:N) = ham_imp
  ! Top-right block:  Δ
  ham_bdg(1:N, N+1:2*N) = gap_imp
  ! Bottom-left block:  Δ†  =  conjg(transpose(gap_imp))
  ham_bdg(N+1:2*N, 1:N) = conjg(transpose(gap_imp))
  ! Bottom-right block:  -H*  =  -conjg(ham_imp)
  ham_bdg(N+1:2*N, N+1:2*N) = -conjg(ham_imp)
end subroutine gen_bdg_ham

! -----------------------------------------------------------------------------
subroutine get_dft_imp_ham(ham_k,ham_imp,klist,rlist,Nk,Nsite,Norb) bind(C)
! Fourier-transform the real-space impurity supercell Hamiltonian ham_imp to
! a (k,k')-space representation ham_k.
!
! The transformation is:
!
!   ham_k(n+(j-1)*Norb, m+(i-1)*Norb)
!     = Σ_{k,l=1}^{Nsite}  ham_imp(n+(l-1)*Norb, m+(k-1)*Norb)
!       × exp(-i 2π (k_j·R_l - k_i·R_k))
!
! where k_i, k_j are k-points and R_k, R_l are site positions (both in
! fractional coordinates).
!
! Only the upper triangle (j≥i in k-point index) is accumulated; the lower
! triangle is set from Hermitian symmetry at the end of each (i,j) pair.
!
! Parallel accumulation over the outer k-point loop uses an OMP reduction
! on ham_k to avoid write races.
!
! Index convention (orbital-fast, k-point-slow):
!   ham_k row/col for orbital n at k-point j  →  n + (j-1)*Norb
! -----------------------------------------------------------------------------
  use constant
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nsite,Norb
  real(c_double),intent(in),dimension(3,Nsite):: rlist   ! site positions (fractional)
  real(c_double),intent(in),dimension(3,Nk):: klist      ! k-point coordinates (fractional)
  complex(c_double),intent(in),dimension(Norb*Nsite,Norb*Nsite):: ham_imp
  complex(c_double),intent(out),dimension(Norb*Nk,Norb*Nk):: ham_k

  integer i,j,k,l,m,n
  real(c_double) phase

  ! Must initialise to zero; the loops below accumulate into ham_k
  ham_k(:,:) = (0.0d0, 0.0d0)
  ! reduction(+:ham_k) is required because multiple iterations of the outer
  ! loop (over i) write to overlapping columns of ham_k
  !$omp parallel do private(j,k,l,m,n,phase) reduction(+:ham_k)
  k_loop1: do i=1,Nk
     k_loop2: do j=i,Nk   ! upper triangle only; lower filled by Hermitian symmetry below
        site_loop1: do k=1,Nsite
           site_loop2: do l=1,Nsite
              ! Phase = 2π (k_j·R_l - k_i·R_k), all in fractional coordinates
              ! The factor e^{-i·phase} implements the Bloch sum
              phase=2*pi*(sum(klist(:,j)*rlist(:,l))-sum(klist(:,i)*rlist(:,k)))
              orb_loop1: do m=1,Norb
                 orb_loop2: do n=1,Norb
                    ! ham_imp index: orbital-fast, site-slow
                    !   row: orbital n at site l  →  n + (l-1)*Norb
                    !   col: orbital m at site k  →  m + (k-1)*Norb
                  ham_k(n+Norb*(j-1),m+Norb*(i-1))=ham_k(n+Norb*(j-1),m+Norb*(i-1))&
                     +ham_imp(n+Norb*(l-1),m+Norb*(k-1))*cmplx(cos(phase),-sin(phase),kind=c_double)
                 end do orb_loop2
              end do orb_loop1
           end do site_loop2
        end do site_loop1
        ! Fill lower triangle: ham_k(i,j) = conjg( ham_k(j,i) )
        orbk_loop: do m=1,Norb
           orbk_loop2: do n=1,Norb
              ham_k(m+Norb*(i-1),n+Norb*(j-1))=conjg(ham_k(n+Norb*(j-1),m+Norb*(i-1)))
           end do orbk_loop2
        end do orbk_loop
     end do k_loop2
  end do k_loop1
  !$omp end parallel do
end subroutine get_dft_imp_ham

! -----------------------------------------------------------------------------
subroutine get_spectrum_spaghetti(spa,uni,eigs,klist,rlist,wlist,Nw,Nk,Nsite,Norb,mu,eta) bind(C)
! Compute the k-resolved spectral function (spaghetti plot) from the
! eigenstates of the impurity supercell Hamiltonian.
!
! The spectral function is obtained by Bloch-unfolding the supercell
! Green's function:
!
!   A(k, ω) = -(1/π) Im G(k, ω)
!
!   G(k, ω) = Σ_{R,R'} Σ_l Σ_n
!              e^{i 2π k·(R-R')}
!              × conjg(uni(l+(R-1)*Norb, n)) × uni(l+(R'-1)*Norb, n)
!              / (ω - ε_n + μ + i·η)
!
! where:
!   uni(:,n)   = n-th eigenvector of ham_imp (columns of the unitary matrix)
!   eigs(n)    = n-th eigenvalue
!   l          = orbital index (1..Norb)
!   R, R'      = site indices used as real-space labels (j, k in the code)
!   spa(m,i)   stores the full complex Green's function; the caller takes -Im/π
!
! Index convention for uni (same as ham_imp, orbital-fast / site-slow):
!   element for orbital l at site j  →  l + (j-1)*Norb
! -----------------------------------------------------------------------------
  use constant
  implicit none
  integer(c_int64_t),intent(in):: Nw,Nk,Nsite,Norb
  real(c_double),intent(in):: eta,mu
  real(c_double),intent(in),dimension(Nw):: wlist
  real(c_double),intent(in),dimension(3,Nsite):: rlist   ! site positions (fractional)
  real(c_double),intent(in),dimension(3,Nk):: klist      ! k-points (fractional)
  real(c_double),intent(in),dimension(Norb*Nsite):: eigs ! supercell eigenvalues
  complex(c_double),intent(in),dimension(Norb*Nsite,Norb*Nsite):: uni  ! eigenvectors (columns)
  complex(c_double),intent(out),dimension(Nw,Nk):: spa   ! spectral function G(k,ω)

  integer i,j,k,l,m,n
  real(c_double) phase

  !$omp parallel
  !$omp workshare
  spa(:,:)=0.0d0
  !$omp end workshare
  !$omp do private(i,j,k,l,n,m,phase)
  k_loop: do i=1,Nk
     site_loop1: do j=1,Nsite       ! site R  (ket side)
        site_loop2: do k=1,Nsite    ! site R' (bra side)
           ! Bloch phase: e^{i 2π k·(R_j - R_k)}
           phase=2*pi*sum(klist(:,i)*(rlist(:,j)-rlist(:,k)))
           orb_loop: do l=1,Norb
              eig_loop: do n=1,Nsite*Norb
                 w_loop:do m=1,Nw
                    ! Contribution from eigenstate n:
                    !   <k,l|n> = Σ_j e^{ik·R_j} conjg(uni(l+(j-1)*Norb, n))
                    !   <n|k,l> = Σ_k' e^{-ik·R_k'} uni(l+(k'-1)*Norb, n)
                    ! uni index: orbital l at site j  →  l + (j-1)*Norb
                  spa(m,i)=spa(m,i)+conjg(uni(l+(j-1)*Norb,n))*uni(l+(k-1)*Norb,n)&
                     *cmplx(cos(phase),-sin(phase),kind=c_double)/cmplx(wlist(m)-eigs(n)+mu,eta,kind=c_double)
                 end do w_loop
              end do eig_loop
           end do orb_loop
        end do site_loop2
     end do site_loop1
  end do k_loop
  !$omp end do
  !$omp end parallel
end subroutine get_spectrum_spaghetti

! -----------------------------------------------------------------------------
subroutine calc_gloc(Gloc,hamk,sigma,z,Nk,Norb) bind(C)
!> calc_gloc
!> Compute local Green's function: Gloc = (1/Nk) sum_k [z*I - H(k) - Sigma]^(-1)
!!@param   Gloc,out: local Green's function [Norb,Norb] complex128
!!@param   hamk, in: k-space Hamiltonian [Norb,Norb,Nk] complex128
!!@param  sigma, in: self-energy [Norb,Norb] complex128
!!@param      z, in: complex frequency (iω_n or ω+iδ)
!!@param     Nk, in: number of k-points
!!@param   Norb, in: number of orbitals
!
! Algorithm:
!   For each k, form gk = z*I - H(k) - Sigma and invert it via LAPACK
!   (zgetrf + zgetri).  The inverses are averaged over k.
!
! OMP strategy:
!   Parallelised over k with a reduction on Gloc and the failure counter
!   nfail.  The 'if(omp_get_active_level()==0)' guard prevents nested
!   parallelism when this routine is called from inside solve_cpa_array,
!   which already runs an OMP loop over frequencies.
! -----------------------------------------------------------------------------
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  use omp_lib, only: omp_get_active_level
  implicit none
  integer(c_int64_t),intent(in):: Nk,Norb
  complex(c_double),intent(in):: z
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: hamk
  complex(c_double),intent(in),dimension(Norb,Norb):: sigma
  complex(c_double),intent(out),dimension(Norb,Norb):: Gloc

  complex(c_double),dimension(Norb,Norb):: gk
  complex(c_double),dimension(Norb*Norb):: work  ! LAPACK workspace (Norb^2 is sufficient)
  integer(c_int32_t),dimension(Norb):: ipiv         ! LAPACK pivot indices
  integer(c_int64_t) :: ik, nfail
  integer(c_int32_t) :: i,j,info

  Gloc(:,:) = (0.0d0, 0.0d0)
  nfail = 0

  !$omp parallel if(omp_get_active_level()==0) private(gk,ipiv,work,info,i,j)
  !$omp do reduction(+:Gloc,nfail)
  do ik = 1, Nk
     ! Build gk = z*I - H(k) - Sigma
     do j = 1, Norb
        do i = 1, Norb
           gk(i,j) = -hamk(i,j,ik) - sigma(i,j)
        end do
        gk(j,j) = gk(j,j) + z
     end do
     ! LU factorisation, then inversion
     call zgetrf(Norb, Norb, gk, Norb, ipiv, info)
     if (info /= 0) then; nfail = nfail + 1; cycle; end if
     call zgetri(Norb, gk, Norb, ipiv, work, Norb*Norb, info)
     if (info /= 0) then; nfail = nfail + 1; cycle; end if
     ! Accumulate gk^{-1} into Gloc
     do j = 1, Norb
        do i = 1, Norb
           Gloc(i,j) = Gloc(i,j) + gk(i,j)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

  if (nfail > 0) print*, 'Warning: calc_gloc: matrix inversion failed for', nfail, '/', Nk, 'k-points'

  Gloc(:,:) = Gloc(:,:) / dble(Nk)
end subroutine calc_gloc

! -----------------------------------------------------------------------------
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
!
! Physical background — single-site CPA:
!   delta_alpha = V_alpha - Sigma  is the scattering potential of species alpha
!   relative to the effective CPA medium described by Sigma.
!
!   The single-site t-matrix in the CPA medium is:
!     t_alpha = delta_alpha (I - G_loc delta_alpha)^{-1}
!
!   The CPA self-consistency condition is:
!     x * t_A + (1-x) * t_B = 0    (average t-matrix vanishes)
!
!   sigma_new = sigma + x*t_A + (1-x)*t_B  is one step of the fixed-point
!   iteration toward this condition.  At convergence sigma_new → sigma.
!
! Matrix ordering note:
!   t_alpha = delta_alpha * M_alpha^{-1}   (delta on the LEFT of the inverse)
!   This is the standard BEB / Velicky-KE single-site form.
! -----------------------------------------------------------------------------
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Norb
  real(c_double),intent(in):: x
  complex(c_double),intent(in),dimension(Norb,Norb):: Gloc,sigma,VA,VB
  complex(c_double),intent(out),dimension(Norb,Norb):: sigma_new

  complex(c_double),dimension(Norb,Norb):: dA,dB,matA,matB,tA,tB
  complex(c_double),dimension(Norb*Norb):: work
  integer(c_int32_t),dimension(Norb):: ipiv
  integer(c_int32_t) :: i,j,info

  ! delta_alpha = V_alpha - Sigma
  dA(:,:) = VA(:,:) - sigma(:,:)
  dB(:,:) = VB(:,:) - sigma(:,:)

  ! --- Species A ---
  ! matA = I - Gloc * delta_A
  call zgemm('N','N',Norb,Norb,Norb,(-1.0d0,0.0d0),Gloc,Norb,dA,Norb,(0.0d0,0.0d0),matA,Norb)
  do i = 1, Norb
     matA(i,i) = matA(i,i) + (1.0d0, 0.0d0)
  end do
  ! Invert matA in-place → matA = (I - Gloc*delta_A)^{-1}
  call zgetrf(Norb, Norb, matA, Norb, ipiv, info)
  if (info /= 0) then; print*, 'CPA: zgetrf(A) failed, info=', info; stop; end if
  call zgetri(Norb, matA, Norb, ipiv, work, Norb*Norb, info)
  if (info /= 0) then; print*, 'CPA: zgetri(A) failed, info=', info; stop; end if
  ! tA = delta_A * (I - Gloc*delta_A)^{-1}
  call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),dA,Norb,matA,Norb,(0.0d0,0.0d0),tA,Norb)

  ! --- Species B ---
  ! matB = I - Gloc * delta_B
  call zgemm('N','N',Norb,Norb,Norb,(-1.0d0,0.0d0),Gloc,Norb,dB,Norb,(0.0d0,0.0d0),matB,Norb)
  do i = 1, Norb
     matB(i,i) = matB(i,i) + (1.0d0, 0.0d0)
  end do
  ! Invert matB in-place → matB = (I - Gloc*delta_B)^{-1}
  call zgetrf(Norb, Norb, matB, Norb, ipiv, info)
  if (info /= 0) then; print*, 'CPA: zgetrf(B) failed, info=', info; stop; end if
  call zgetri(Norb, matB, Norb, ipiv, work, Norb*Norb, info)
  if (info /= 0) then; print*, 'CPA: zgetri(B) failed, info=', info; stop; end if
  ! tB = delta_B * (I - Gloc*delta_B)^{-1}
  call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),dB,Norb,matB,Norb,(0.0d0,0.0d0),tB,Norb)

  ! sigma_new = sigma + x*tA + (1-x)*tB
  ! (One fixed-point step; converges when x*tA + (1-x)*tB → 0)
  do j = 1, Norb
     do i = 1, Norb
        sigma_new(i,j) = sigma(i,j) + x*tA(i,j) + (1.0d0-x)*tB(i,j)
     end do
  end do
end subroutine calc_tmat_cpa

! -----------------------------------------------------------------------------
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
!!@param        pp,   in: linear mixing parameter (0<pp<=1; pp=1 is no mixing)
!!@param        Nk,   in: number of k-points
!!@param      Norb,   in: number of orbitals
!!@param   maxiter,   in: maximum iterations
!!@param       tol,   in: convergence tolerance
!
! Convergence metric:
!   diff = ||sigma_new - sigma_cpa||_F  =  ||x*tA + (1-x)*tB||_F
!   This is the residual of the CPA self-consistency condition directly.
!
! Linear mixing:
!   sigma_cpa <- (1-pp)*sigma_cpa + pp*sigma_new
!   Applied after convergence check; pp=1 gives the unmodified fixed-point
!   update, smaller pp stabilises iteration near singular points.
! -----------------------------------------------------------------------------
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,Norb,maxiter
  real(c_double),intent(in):: x,pp,tol
  complex(c_double),intent(in):: z
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: hamk
  complex(c_double),intent(in),dimension(Norb,Norb):: VA,VB
  complex(c_double),intent(inout),dimension(Norb,Norb):: sigma_cpa

  complex(c_double),dimension(Norb,Norb):: Gloc,sigma_new
  real(c_double) diff
  integer(c_int64_t) iter,i,j

  interface
     subroutine calc_gloc(Gloc,hamk,sigma,z,Nk,Norb) bind(C)
       import:: c_int64_t,c_double
       integer(c_int64_t),intent(in):: Nk,Norb
       complex(c_double),intent(in):: z
       complex(c_double),intent(in),dimension(Norb,Norb,Nk):: hamk
       complex(c_double),intent(in),dimension(Norb,Norb):: sigma
       complex(c_double),intent(out),dimension(Norb,Norb):: Gloc
     end subroutine calc_gloc
     subroutine calc_tmat_cpa(sigma_new,Gloc,sigma,VA,VB,x,Norb) bind(C)
       import:: c_int64_t,c_double
       integer(c_int64_t),intent(in):: Norb
       real(c_double),intent(in):: x
       complex(c_double),intent(in),dimension(Norb,Norb):: Gloc,sigma,VA,VB
       complex(c_double),intent(out),dimension(Norb,Norb):: sigma_new
     end subroutine calc_tmat_cpa
  end interface

  diff = huge(diff)  ! ensure the non-convergence warning fires if maxiter=0
  do iter = 1, maxiter
     ! Step 1: local Green's function at current Sigma_CPA
     call calc_gloc(Gloc, hamk, sigma_cpa, z, Nk, Norb)

     ! Steps 2–3: compute t-matrices and candidate new self-energy
     call calc_tmat_cpa(sigma_new, Gloc, sigma_cpa, VA, VB, x, Norb)

     ! Convergence check: Frobenius norm of the CPA residual x*tA + (1-x)*tB
     diff = 0.0d0
     do j = 1, Norb
        do i = 1, Norb
           diff = diff + abs(sigma_new(i,j) - sigma_cpa(i,j))**2
        end do
     end do
     diff = sqrt(diff)

     ! Linear mixing before the next iteration
     do j = 1, Norb
        do i = 1, Norb
           sigma_cpa(i,j) = (1.0d0-pp)*sigma_cpa(i,j) + pp*sigma_new(i,j)
        end do
     end do

     if (diff < tol) exit
  end do
  if (diff >= tol) print*, 'Warning: solve_cpa: not converged after', maxiter, 'iterations, diff=', diff
end subroutine solve_cpa

! -----------------------------------------------------------------------------
subroutine solve_cpa_array(sigma_cpa_w,hamk,VA,VB,x,zlist,pp,Nk,Norb,Nw,maxiter,tol) bind(C)
!> solve_cpa_array
!> Run CPA self-consistent loop for an array of frequencies (Matsubara or real frequency).
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
!
! Parallelism:
!   The OMP loop is over the Nw frequencies.  Each frequency is independent,
!   so no reduction or synchronisation is needed.  Each thread works on a
!   contiguous slice sigma_cpa_w(:,:,iw) (Fortran column-major: the first two
!   indices are contiguous in memory), so there are no write races.
!
!   Inside solve_cpa → calc_gloc, a second OMP region is guarded by
!   omp_get_active_level()==0, so it is disabled here, keeping the nesting
!   flat (frequency-level parallelism only).
! -----------------------------------------------------------------------------
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,Norb,Nw,maxiter
  real(c_double),intent(in):: x,pp,tol
  complex(c_double),intent(in),dimension(Nw):: zlist
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: hamk
  complex(c_double),intent(in),dimension(Norb,Norb):: VA,VB
  complex(c_double),intent(inout),dimension(Norb,Norb,Nw):: sigma_cpa_w

  integer(c_int64_t) iw

  interface
     subroutine solve_cpa(sigma_cpa,hamk,VA,VB,x,z,pp,Nk,Norb,maxiter,tol) bind(C)
       import:: c_int64_t,c_double
       integer(c_int64_t),intent(in):: Nk,Norb,maxiter
       real(c_double),intent(in):: x,pp,tol
       complex(c_double),intent(in):: z
       complex(c_double),intent(in),dimension(Norb,Norb,Nk):: hamk
       complex(c_double),intent(in),dimension(Norb,Norb):: VA,VB
       complex(c_double),intent(inout),dimension(Norb,Norb):: sigma_cpa
     end subroutine solve_cpa
  end interface

  !$omp parallel do private(iw)
  do iw = 1, Nw
     ! sigma_cpa_w(:,:,iw) is contiguous in memory (Fortran column-major);
     ! each thread writes to a distinct frequency slice — no race condition.
     call solve_cpa(sigma_cpa_w(:,:,iw), hamk, VA, VB, x, zlist(iw), pp, Nk, Norb, maxiter, tol)
  end do
  !$omp end parallel do
end subroutine solve_cpa_array
