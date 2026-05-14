subroutine mkBdGhamk(hamBdGk,hamk,delta,Nk,Norb) bind(C)
  !> Build BdG Hamiltonian in Nambu basis Ψ_k = (c_{k,↑,1..N}, c†_{-k,↓,1..N})^T
  !>
  !>   H_BdG(k) = [  H(k)        Δ(k)    ]
  !>              [  Δ†(k)     -H^T(-k)  ]
  !>
  !> The lower-right block uses -conjg(H(k)) which equals -H^T(-k) ONLY when
  !> H is inversion-symmetric (H(-k)=H(k)) combined with Hermiticity. This
  !> implementation is therefore restricted to centrosymmetric systems
  !> (no Rashba-type SOC with broken inversion).
  !>
  !> The Nambu basis (c_↑, c†_↓) covers S_z=0 pairing only:
  !>   - singlet (Δ_↑↓ = -Δ_↓↑)
  !>   - triplet d-vector ‖ ẑ  (Δ_↑↓ = +Δ_↓↑)
  !> S_z=±1 triplets (d_x, d_y, chiral p, etc.) require a different basis.
  !>
  !!@param hamBdGk,out: 2Norb x 2Norb BdG Hamiltonian on k-mesh
  !!@param    hamk,in: normal-state Hamiltonian (typically with -mu*I subtracted)
  !!@param   delta,in: gap function Δ(k) in the (Δ_↑↓) component
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,Norb
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: hamk,delta
  complex(c_double),intent(out),dimension(2*Norb,2*Norb,Nk):: hamBdGk

  integer(c_int32_t) k,l,m

  do k=1,Nk
     do l=1,Norb
        do m=1,Norb
           hamBdGk(l,m,k)=hamk(l,m,k)                          !  H(k)            (particle block)
           hamBdGk(l+Norb,m+Norb,k)=-conjg(hamk(l,m,k))         ! -H*(k) = -H^T(-k) under inversion sym
           hamBdGk(l,m+Norb,k)=delta(l,m,k)                    !  Δ_{lm}(k)        (anomalous, upper-right)
           hamBdGk(l+Norb,m,k)=conjg(delta(m,l,k))              !  (Δ†)_{lm}(k)     (anomalous, lower-left)
        end do
     end do
  end do
end subroutine mkBdGhamk

module calc_irr_chi_sc
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
contains
  function irr_chi_sc(Nk,Norb,Nchi,uni,eig,ffermi,ol,temp,qshift,w,idelta,eps,sw_spsym)
    !> Irreducible spin susceptibility in the SC state, evaluated in BdG eigenbasis.
    !>
    !> Splits into a normal (G·G) and anomalous (F·F†) contribution:
    !>   chi_{ij}(q,w) = (1/Nk) Σ_k Σ_{l,m=1..2Norb} weight(l,m,k,q,w)
    !>                * [ B_i A_j + sgn * R_i P_j ]
    !>   l : BdG band at k+q,  m : BdG band at k
    !>   weight = [f_l(k+q) - f_m(k)] / [w + E_m(k) - E_l(k+q) + iδ]   (= -Π_0 convention)
    !>   degenerate static limit: lim → f(1-f)/T = -∂f/∂E   (intra-band)
    !>
    !> Sign convention from spin-space symmetry:
    !>   sw_spsym = .false. → singlet     → sgn = +1  (Yosida: χ_zz → 0 at T=0)
    !>   sw_spsym = .true.  → triplet d∥ẑ → sgn = -1  (χ_zz preserved at T=0)
    !> At Δ→0, BdG eigenvectors become block-diagonal in Nambu space, so P/R
    !> cannot both be non-zero for the same (l,m); FF vanishes and GG reduces to
    !> the normal-state expression.
    !!@param        Nk: The number of k-points
    !!@param      Norb: The number of orbitals
    !!@param      Nchi: The footnote of chi
    !!@param       uni: unitary that diagonalizes H_BdG; stored as uni(orbital,band,k)
    !!@param       eig: energies of bogoliubov quasi-particle bands [2Norb,Nk]
    !!@param    ffermi: fermi distribution evaluated at eig [2Norb,Nk]
    !!@param        ol: the list of the properties of orbitals at footnote of chi
    !!@param      temp: temperature (regularized by temp_safe to avoid /0)
    !!@param    qshift: q-shifted klist; qshift(k) returns the index of k+q
    !!@param         w: frequency
    !!@param    idelta: dumping factor (Lorentzian broadening)
    !!@param       eps: skip threshold for tiny |Δf| contributions (performance)
    !!@param   sw_spsym: symmetry of spin space (true: triplet_dz, false: singlet)
    !!@return irr_chi_sc: irreducible susceptibility matrix [Nchi,Nchi]
    integer(c_int64_t),intent(in):: Nk,Norb,Nchi
    integer(c_int64_t),intent(in),dimension(Nk):: qshift
    integer(c_int64_t),intent(in),dimension(Nchi,2):: ol
    real(c_double),intent(in):: temp,eps,idelta,w
    logical,intent(in):: sw_spsym
    real(c_double),intent(in),dimension(2*Norb,Nk):: eig,ffermi
    complex(c_double),intent(in),dimension(2*Norb,2*Norb,Nk):: uni

    integer(c_int32_t) i,k,l,m,nchi32
    real(c_double) temp_safe,w_eps,sgn
    complex(c_double) weight
    ! A,B : matrix-element vectors for the normal G·G contribution.
    !       chi += weight * B(i)*A(j) accumulates u*(o1_i,l,k+q) u(o1_j,l,k+q)
    !                                           * u(o2_i,m,k)   u*(o2_j,m,k)
    ! P,R : matrix-element vectors for the anomalous F·F† contribution
    !       (mix upper Norb / lower Norb across k+q and k).
    complex(c_double),dimension(Nchi):: A_vec,B_vec,P_vec,R_vec
    complex(c_double),dimension(Nchi,Nchi):: chi,irr_chi_sc

    if(sw_spsym)then
       sgn=-1.0d0 !triplet_dz: FF doubles GG -> χ_zz unsuppressed for d∥ẑ
    else
       sgn=+1.0d0 !singlet:    FF cancels GG at T=0 -> Yosida function
    end if
    temp_safe=max(temp,1.0d-12)   ! guard against T=0 in degenerate static limit
    w_eps=1.0d-12                 ! threshold to detect static (w≈0) branch
    nchi32=int(Nchi,c_int32_t)
    chi(:,:)=0.0d0
    !$omp parallel do reduction(+:chi) private(l,m,i,weight,A_vec,B_vec,P_vec,R_vec)
    kloop: do k=1,Nk
       band1_loop: do l=1,2*Norb       ! BdG band index at k+q
          band2_loop: do m=1,2*Norb    ! BdG band index at k
             ! ----- weight: Lindhard-like factor ----------------------------------
             if(abs(w)<w_eps .and. abs(eig(m,k)-eig(l,qshift(k)))<1.0d-9)then
                ! degenerate intra-band static limit: 0/0 -> -∂f/∂E = f(1-f)/T
                weight=cmplx(ffermi(m,k)*(1.0d0-ffermi(m,k))/temp_safe,0.0d0,kind=c_double)
             else if(abs(ffermi(l,qshift(k))-ffermi(m,k))>eps)then
                ! generic case; idelta provides the retarded prescription
                weight=(ffermi(l,qshift(k))-ffermi(m,k))&
                     /cmplx(w+eig(m,k)-eig(l,qshift(k)),idelta,kind=c_double)
             else
                cycle band2_loop       ! both Fermi factors equal -> negligible
             end if
             ! ----- BdG matrix elements -------------------------------------------
             ! uni(a, n, k):  a∈[1..Norb]   = particle component  (c_↑)
             !                a∈[Norb+1..2Norb] = hole component  (c†_↓)
             ! ol(i,1)=o1_i, ol(i,2)=o2_i are the orbital indices for chi index i.
             do i=1,nchi32
                A_vec(i)=uni(ol(i,1),l,qshift(k))*conjg(uni(ol(i,2),m,k))         ! u(o1,l,k+q) u*(o2,m,k)
                B_vec(i)=conjg(uni(ol(i,1),l,qshift(k)))*uni(ol(i,2),m,k)         ! conjugate of A_vec
                P_vec(i)=uni(ol(i,1),l,qshift(k))*conjg(uni(ol(i,2)+Norb,m,k))    ! particle@k+q × hole@k
                R_vec(i)=uni(ol(i,2)+Norb,l,qshift(k))*conjg(uni(ol(i,1),m,k))    ! hole@k+q × particle@k
             end do
             ! ----- accumulate rank-1 outer products ------------------------------
             ! zgeru: chi := chi + alpha * x * y^T  (no conjugation on y)
             call zgeru(nchi32,nchi32,weight,B_vec,1,A_vec,1,chi,nchi32)          ! G·G part
             call zgeru(nchi32,nchi32,sgn*weight,R_vec,1,P_vec,1,chi,nchi32)      ! F·F† part
          end do band2_loop
       end do band1_loop
    end do kloop
    !$omp end parallel do
    irr_chi_sc=chi(:,:)/Nk            ! k-mesh average (matches normal-state convention)
  end function irr_chi_sc
end module calc_irr_chi_sc

subroutine get_chi_irr_sc(chi,uni,eig,ffermi,qshift,ol,wl,Nchi,Norb,Nk,Nw,idelta,eps,temp,sw_spsym) bind(C)
  !> Driver: evaluates irr_chi_sc on a frequency mesh wl(1..Nw) at a fixed q.
  !> The eigenvalues / eigenvectors / Fermi factors are computed once outside
  !> this routine and reused for every ω; only the energy denominator depends on ω.
  !!@param   chi,out: irreducible susceptibility
  !!@param    uni,in: unitary matrix
  !!@param    eig,in: energies of bogoliubov quasi-particle bands
  !!@param ffermi,in: fermi distribute function
  !!@param qshift,in: q-shifted klist
  !!@param     ol,in: the list of the properties of orbitals at footnote of chi
  !!@param     wl,in: the list of frequency
  !!@param   Nchi,in: The footnote of chi
  !!@param   Norb,in: The number of orbitals
  !!@param     Nk,in: The number of k-points
  !!@param     Nw,in: The number of frequency mesh
  !!@param idelta,in: dumping factor
  !!@param    eps,in: threshold of calculation value
  !!@param   temp,in: temperature
  !!@param sw_spsym,in: symmetry of spin space (true: triplet_dz, false: singlet)
  use calc_irr_chi_sc
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,Norb,Nw,Nchi
  integer(c_int64_t),intent(in),dimension(Nk):: qshift
  integer(c_int64_t),intent(in),dimension(Nchi,2):: ol
  logical,intent(in):: sw_spsym
  real(c_double),intent(in):: temp,eps,idelta
  real(c_double),intent(in),dimension(2*Norb,Nk):: eig,ffermi
  real(c_double),intent(in),dimension(Nw):: wl
  complex(c_double),intent(in),dimension(2*Norb,2*Norb,Nk):: uni
  complex(c_double),intent(out),dimension(Nchi,Nchi,Nw):: chi

  integer(c_int32_t) i
 
  !$omp parallel do private(i)
  wloop: do i=1,Nw
     chi(:,:,i)=irr_chi_sc(Nk,Norb,Nchi,uni,eig,ffermi,ol,temp,qshift,wl(i),idelta,eps,sw_spsym)
  end do wloop
  !$omp end parallel do
end subroutine get_chi_irr_sc

subroutine get_band_to_orb_delta(delta,init_delta,uni,Nk,Norb) bind(C)
  !> Transform the gap function from band-diagonal representation to orbital basis.
  !>
   !> For the row-eigenvector convention used in this codebase,
   !>
   !>   U(k) H(k) = E(k) U(k)
   !>
   !> and the orbital-basis pairing matrix is
  !>
   !>   delta_orb(k)_{ab} = sum_n U_{n,a}(k) delta_n(k) U_{n,b}(k)
   !>                     = U(k)^T * diag(delta_n(k)) * U(k)
  !>
   !> where delta_n(k) are complex band-diagonal gap amplitudes and U(k) are
   !> row-eigenvectors. This is the inverse of the band-projection step in
   !> get_initial_delta / conv_delta_orb_to_band: given a band-diagonal gap,
   !> it returns the full Norb x Norb matrix in orbital basis.
  !>
  !!@param     delta,out: gap function in orbital basis [Norb, Norb, Nk]
  !!@param init_delta,in: band-diagonal gap amplitudes [Nk, Norb] complex128
   !!@param       uni,in: row-eigenvectors; uni(band, orbital, k) [Norb, Norb, Nk] complex128
  !!@param        Nk,in: number of k-points
  !!@param      Norb,in: number of orbitals (= number of bands)
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Norb,Nk
  complex(c_double),intent(in),dimension(Nk,Norb):: init_delta
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(c_double),intent(out),dimension(Norb,Norb,Nk):: delta

   integer(c_int32_t) i,a,b,n
   complex(c_double) zsum
  !$omp parallel
  !$omp workshare
  delta(:,:,:)=0.0d0
  !$omp end workshare
  !$omp end parallel

  !$omp parallel do private(i,a,b,n,zsum)
  do i=1,Nk
     do a=1,Norb
        do b=1,Norb
           zsum=(0.0d0,0.0d0)
           do n=1,Norb
              zsum=zsum+uni(n,a,i)*init_delta(i,n)*uni(n,b,i)
           end do
           delta(a,b,i)=zsum
        end do
     end do
  end do
  !$omp end parallel do
end subroutine get_band_to_orb_delta
