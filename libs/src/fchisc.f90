subroutine mkBdGhamk(hamBdGk,hamk,delta,Nk,Norb) bind(C)
  !> Build BdG Hamiltonian in Nambu basis Ψ_k = (c_{k,↑,1..N}, c†_{-k,↓,1..N})^T
  !>
  !>   H_BdG(k) = [  H(k)        Δ(k)    ]
  !>              [  Δ†(k)      -H(k)    ]
  !>
  !> The hole block is -H^T(-k) by construction (the Nambu doubling reorders
  !> c_{-k} H(-k) c†_{-k}), and for a spin-round time-reversal-symmetric model that
  !> IS -H(k), with no further assumption:
  !>       [H^T(-k)]_lm = [H(-k)]_ml = [H_ml(k)]*  (TRS: H(-k) = H*(k))
  !>                                  = H_lm(k)    (hermiticity: H_ml = H_lm*)
  !> so H^T(-k) = H(k) exactly.  In particular NO inversion symmetry is needed and
  !> -k is never referenced.  An earlier version used -conjg(H(k)) = -H*(k), which by
  !> the same TRS relation is -H(-k): that differs from -H^T(-k) by a transpose and is
  !> correct only when H(k) is real, i.e. H(-k)=H(k).  That excluded exactly the
  !> multi-site cells this code targets -- with several sites per cell, inversion
  !> EXCHANGES sites, so H(-k) = P^-1 H(k) P with a non-trivial signed permutation P
  !> (see plibs.inversion_op) and H(-k) != H(k) even though the crystal is
  !> centrosymmetric.  Measured on an SSH chain (t1=1.0, t2=0.6, |Δ|=0.25) the old
  !> form moved the BdG eigenvalues by 0.064, i.e. 26% of the gap.
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
           hamBdGk(l+Norb,m+Norb,k)=-hamk(l,m,k)                ! -H(k) = -H^T(-k) by TRS + hermiticity
           hamBdGk(l,m+Norb,k)=delta(l,m,k)                    !  Δ_{lm}(k)        (anomalous, upper-right)
           hamBdGk(l+Norb,m,k)=conjg(delta(m,l,k))              !  (Δ†)_{lm}(k)     (anomalous, lower-left)
        end do
     end do
  end do
end subroutine mkBdGhamk

module calc_irr_chi_sc
  use occupation, only:occ_factor
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double,c_bool
  implicit none
contains
  function irr_chi_sc(Nk,Norb,Nchi,uni,eig,ffermi,ol,temp,qshift,w,idelta,eps,sw_spsym)
    !> Irreducible spin susceptibility in the SC state, evaluated in BdG eigenbasis.
    !>
    !> In the Nambu basis ψ_k = (c_{k↑}, c†_{-k↓}) the spin density S_z(q) is
    !> Σ_k ψ†_{k+q} ψ_k: the ↓ part, -c†_{k+q,a↓} c_{k,b↓} = +c_{k,b↓} c†_{k+q,a↓}, is
    !> a hole-hole bilinear (orbital pair swapped, k relabelled -k-q).  Its matrix
    !> element between BdG state l at k+q and m at k, for chi index i = (a,b), is
    !>   X_i = p*(a,l,k+q) p(b,m,k)  +  h*(b,l,k+q) h(a,m,k)
    !> (p = particle rows 1..Norb, h = hole rows Norb+1..2Norb of uni), and
    !>   chi_{ij}(q,w) = (1/Nk) Σ_k Σ_{l,m=1..2Norb} weight(l,m,k,q,w) * X_i X_j* / 2
    !>   l : BdG band at k+q,  m : BdG band at k
    !>   weight = [f_l(k+q) - f_m(k)] / [w + E_m(k) - E_l(k+q) + iδ]   (= -Π_0 convention)
    !>   degenerate static limit: lim → f(1-f)/T = -∂f/∂E   (intra-band)
    !> The 1/2 keeps the old normalisation (one spin species, = chi0 at Δ=0).
    !> chi(1) holds the pp·pp + hh·hh ("G·G") terms, chi(2) the pp·hh cross ("F·F†")
    !> terms, which vanish at Δ=0 because each BdG state is then pure p or pure h.
    !> Every state enters each product once plain and once conjugated, so the result
    !> is invariant under the arbitrary phase LAPACK gives each eigenvector.  The
    !> earlier F·F† term, h(l,k+q)p(l,k+q) * p*(m,k)h*(m,k), was not: it picked up
    !> e^{2iθ_l - 2iθ_m}, harmless only for real eigenvectors (real H(k), as in the
    !> one-orbital tests) and gauge noise of either sign for a complex multi-orbital
    !> H(k), where it drove 1/T1T negative at low T.
    !>
    !> Sign convention from spin-space symmetry.  What sgn selects is unambiguous:
    !>   sw_spsym = .false. → sgn = +1  FF cancels GG at T=0 → Yosida suppression
    !>   sw_spsym = .true.  → sgn = -1  FF adds to GG        → response survives at T=0
    !> Mapping onto the field direction, via χ_ij = χ_N[δ_ij - d̂_i d̂_j (1-Y(T))]:
    !>   singlet   : χ is isotropic and Yosida-suppressed for every field direction,
    !>               so .false. is the only meaningful setting.
    !>   triplet   : the component PARALLEL to d is Yosida-suppressed (.false.), the
    !>               components PERPENDICULAR to d are unaffected (.true.).
    !>               For d∥ẑ that means H∥z → .false., H⊥z → .true.  (This is the
    !>               Sr2RuO4 reading: d∥ẑ leaves the IN-PLANE Knight shift flat
    !>               through Tc while χ_zz still drops.)
    !> NOTE: an earlier version of this comment labelled sgn=-1 as "χ_zz preserved
    !> for d∥ẑ".  That is the wrong component -- for d∥ẑ it is χ_xx=χ_yy that is
    !> preserved and χ_zz that is suppressed, exactly like the singlet.  Only the
    !> label was wrong; the arithmetic below is unchanged.
    !> sgn multiplies the cross terms; the arithmetic choice it encodes is unchanged.
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
    !!@param   sw_spsym: spin channel: .false. = Yosida-suppressed (singlet, or H ∥ d),
    !!                     .true.  = preserved at T=0 (triplet, H ⊥ d)
    !!@return irr_chi_sc: irreducible susceptibility matrix [Nchi,Nchi]
    integer(c_int64_t),intent(in):: Nk,Norb,Nchi
    integer(c_int64_t),intent(in),dimension(Nk):: qshift
    integer(c_int64_t),intent(in),dimension(Nchi,2):: ol
    real(c_double),intent(in):: temp,eps,idelta,w
    logical(c_bool),intent(in):: sw_spsym
    real(c_double),intent(in),dimension(2*Norb,Nk):: eig,ffermi
    complex(c_double),intent(in),dimension(2*Norb,2*Norb,Nk):: uni

    integer(c_int32_t) i,k,l,m,nchi32
    real(c_double) temp_safe,w_eps,sgn,el,em,fl,fm,pocc
    complex(c_double) weight
    ! B,P : the particle and hole halves of X_i (see above); A,R their conjugates,
    !       so chi += w/2 * x * conj(y)^T over (x,y) in {B,P}x{B,P}.
    complex(c_double),dimension(Nchi):: A_vec,B_vec,P_vec,R_vec
    complex(c_double) whalf
    complex(c_double),dimension(2,Nchi,Nchi):: chi,irr_chi_sc

    if(sw_spsym)then
       sgn=-1.0d0 !FF adds to GG    -> χ survives at T=0 (triplet, H ⊥ d)
    else
       sgn=+1.0d0 !FF cancels GG    -> Yosida function   (singlet, or H ∥ d)
    end if
    temp_safe=max(temp,1.0d-12)   ! guard against T=0 in degenerate static limit
    w_eps=1.0d-12                 ! threshold to detect static (w≈0) branch
    nchi32=int(Nchi,c_int32_t)
    chi(:,:,:)=0.0d0
    !$omp parallel do reduction(+:chi) private(l,m,i,weight,whalf,A_vec,B_vec,P_vec,R_vec,el,em,fl,fm,pocc)
    kloop: do k=1,Nk
       band1_loop: do l=1,2*Norb       ! BdG band index at k+q
          band2_loop: do m=1,2*Norb    ! BdG band index at k
             ! ----- weight: Lindhard-like factor ----------------------------------
             el=eig(l,qshift(k)); em=eig(m,k)
             fl=ffermi(l,qshift(k)); fm=ffermi(m,k)
             if(abs(w)<w_eps)then
                ! Static limit: real and idelta-free (see calc_chi in fchi.f90).
                ! occ_factor -> -∂f/∂E continuously, so a BdG pair split by less
                ! than idelta keeps its weight.  This matters more here than in
                ! the normal state: the +-E Bogoliubov branches become degenerate
                ! wherever the gap has a node.
                pocc=occ_factor(el,em,fl,fm,temp_safe)
                if(abs(pocc)<=eps)cycle band2_loop
                weight=cmplx(pocc,0.0d0,kind=c_double)
             else if(abs(fl-fm)>eps)then
                ! generic case; idelta provides the retarded prescription
                weight=(fl-fm)/cmplx(w+em-el,idelta,kind=c_double)
             else
                cycle band2_loop       ! both Fermi factors equal -> negligible
             end if
             ! ----- BdG matrix elements -------------------------------------------
             ! uni(a, n, k):  a∈[1..Norb]   = particle component  (c_↑)
             !                a∈[Norb+1..2Norb] = hole component  (c†_↓)
             ! ol(i,1)=o1_i, ol(i,2)=o2_i are the orbital indices for chi index i.
             do i=1,nchi32
                B_vec(i)=conjg(uni(ol(i,1),l,qshift(k)))*uni(ol(i,2),m,k)              ! p*(o1,l,k+q) p(o2,m,k)
                P_vec(i)=conjg(uni(ol(i,2)+Norb,l,qshift(k)))*uni(ol(i,1)+Norb,m,k)    ! h*(o2,l,k+q) h(o1,m,k)
             end do
             A_vec(:)=conjg(B_vec(:))
             R_vec(:)=conjg(P_vec(:))
             ! ----- accumulate rank-1 outer products ------------------------------
             ! zgeru: chi := chi + alpha * x * y^T  (no conjugation on y)
             whalf=0.5d0*weight
             call zgeru(nchi32,nchi32,whalf,B_vec,1,A_vec,1,chi(1,:,:),nchi32)          ! pp·pp
             call zgeru(nchi32,nchi32,whalf,P_vec,1,R_vec,1,chi(1,:,:),nchi32)          ! hh·hh
             call zgeru(nchi32,nchi32,sgn*whalf,B_vec,1,R_vec,1,chi(2,:,:),nchi32)      ! pp·hh cross
             call zgeru(nchi32,nchi32,sgn*whalf,P_vec,1,A_vec,1,chi(2,:,:),nchi32)      ! hh·pp cross
          end do band2_loop
       end do band1_loop
    end do kloop
    !$omp end parallel do
    irr_chi_sc(:,:,:)=chi(:,:,:)/Nk            ! k-mesh average (matches normal-state convention)
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
  !!@param sw_spsym,in: spin channel: .false. = Yosida-suppressed (singlet, or H ∥ d),
  !!                   .true.  = preserved at T=0 (triplet, H ⊥ d)
  use calc_irr_chi_sc
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double,c_bool
  implicit none
  integer(c_int64_t),intent(in):: Nk,Norb,Nw,Nchi
  integer(c_int64_t),intent(in),dimension(Nk):: qshift
  integer(c_int64_t),intent(in),dimension(Nchi,2):: ol
  logical(c_bool),intent(in):: sw_spsym
  real(c_double),intent(in):: temp,eps,idelta
  real(c_double),intent(in),dimension(2*Norb,Nk):: eig,ffermi
  real(c_double),intent(in),dimension(Nw):: wl
  complex(c_double),intent(in),dimension(2*Norb,2*Norb,Nk):: uni
  complex(c_double),intent(out),dimension(2,Nchi,Nchi,Nw):: chi

  integer(c_int32_t) i
 
  !$omp parallel do private(i)
  wloop: do i=1,Nw
     chi(:,:,:,i)=irr_chi_sc(Nk,Norb,Nchi,uni,eig,ffermi,ol,temp,qshift,wl(i),idelta,eps,sw_spsym)
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
   !>   delta_orb(k)_{ab} = sum_n conjg(U_{n,a}(k)) delta_n(k) U_{n,b}(k)
   !>                     = U(k)^dagger * diag(delta_n(k)) * U(k)
  !>
   !> where delta_n(k) are complex band-diagonal gap amplitudes and U(k) are
   !> row-eigenvectors.  This is the form for which the BdG matrix of mkBdGhamk,
   !> [[H, delta^dagger], [delta, -H]] in Python index order, splits into the bands
   !> with E = +-sqrt(xi_n^2 + |delta_n|^2).  (U^T diag U, used before, differs from it
   !> whenever the eigenvectors are complex.)
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
              zsum=zsum+conjg(uni(n,a,i))*init_delta(i,n)*uni(n,b,i)
           end do
           delta(a,b,i)=zsum
        end do
     end do
  end do
  !$omp end parallel do
end subroutine get_band_to_orb_delta
