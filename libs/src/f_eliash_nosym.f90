! Linearized Eliashberg equation WITHOUT any lattice/time-reversal symmetry.
!
! Every array lives on the full Brillouin zone (Nkall points, laid out by kmap),
! in the spin-orbital basis of the SOC code (Norb includes spin, Vmat from
! get_Vmat_soc).  Nothing here relates k to -k through TRS or inversion, so
! Zeeman fields, magnetic order, chiral states, etc. are all allowed.
!
! Only two relations are used to recover negative Matsubara frequencies from
! the stored positive ones; neither needs TRS:
!   (1) Hermiticity of H and of the bare vertex:
!         G(k,-iw) = G(k,iw)^+ ,  chi(q,-iv) = chi(q,iv)^+  (pair-index matrix)
!   (2) Fermion antisymmetry (Pauli), for pairs (k+Q,-k):
!         F_ab(k,-iw) = -F_ba(-k-Q,iw) ,  Delta_ab(k,-iw) = -Delta_ba(-k-Q,iw)
! so Gk/delta/fk are stored as (Nkall,Nw,Norb,Norb) with w_n>0 only, and
! chi/V_Delta as (Nkall,Nw,Nchi,Nchi) with v_n>=0 only.
!
! Finite pair momentum (FFLO / FF): the pair is labelled (k+Q,-k) rather than
! (k+Q/2,-k+Q/2), so Q only has to be a mesh vector (iQ, integer mesh units),
! not Q/2.  The kernel stays a convolution in k-k' because both legs carry the
! same transfer.  chi0/V_Delta do not depend on Q and are reused for every Q.
! iQ=(0,0,0) is the uniform BCS pairing.
!
! Conventions (pair-index storage, F = -G Delta G^T, overall sign of the kernel)
! follow feliashwithsoc.f90, including the pairing vertex V_Delta = -V/2 + V chi V
! (Nyquist slot -V/2, see get_vdelta_nosym), so in a TRS system the result
! reproduces lin_eliash_soc, and its singlet channel reproduces lin_eliash.

subroutine get_kshift_nosym(iks,kmap,ksgn,iQ,Nkall,Nx,Ny,Nz)
  !> index of ksgn*k+Q on the full k-list, built from the integer mesh coordinates
  !!@param  iks,out: full-list index of ksgn*k+Q for each k
  !!@param  kmap,in: integer mesh coordinates of each k-point
  !!@param  ksgn,in: sign of k (+1 or -1)
  !!@param    iQ,in: shift Q in integer mesh units
  !!@param Nkall,in: Number of all k-points
  !!@param    Nx,in: Number of kx mesh
  !!@param    Ny,in: Number of ky mesh
  !!@param    Nz,in: Number of kz mesh
  use,intrinsic:: iso_c_binding, only:c_int64_t
  implicit none
  integer(c_int64_t),intent(in):: Nkall,Nx,Ny,Nz,ksgn
  integer(c_int64_t),intent(in),dimension(3):: iQ
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap
  integer(c_int64_t),intent(out),dimension(Nkall):: iks

  integer(c_int64_t) i
  integer(c_int64_t),allocatable,dimension(:,:,:):: iflat

  allocate(iflat(0:Nx-1,0:Ny-1,0:Nz-1))
  if(Nkall/=Nx*Ny*Nz)then
     print*,'Error(get_kshift_nosym): Nkall /= Nx*Ny*Nz; nosym routines need the full BZ'
     stop
  end if
  iflat(:,:,:)=0
  do i=1,Nkall
     iflat(kmap(1,i),kmap(2,i),kmap(3,i))=i
  end do
  !$omp parallel do private(i)
  do i=1,Nkall
     iks(i)=iflat(modulo(ksgn*kmap(1,i)+iQ(1),Nx),modulo(ksgn*kmap(2,i)+iQ(2),Ny),modulo(ksgn*kmap(3,i)+iQ(3),Nz))
  end do
  !$omp end parallel do
  if(any(iks(:)==0))then
     print*,'Error(get_kshift_nosym): kmap does not cover the full mesh'
     stop
  end if
end subroutine get_kshift_nosym

subroutine get_chi0_nosym(chi,Gk,kmap,olist,temp,Nx,Ny,Nz,Nw,Nkall,Norb,Nchi) bind(C)
  !> irreducible susceptibility chi0 on the full BZ without TRS
  !> chi(q,iv,m,l) = chi0_{(l1 l2),(m1 m2)} = -(T/N) sum_k G_{l1 m1}(k+q,iw+iv) G_{m2 l2}(k,iw)
  !> Only Hermiticity is used: one FFT of pair (l,m) gives chi(.,v>=0,m,l) and,
  !> from its v<0 half, chi(.,v>=0,l,m)=conjg(chi(.,-v,m,l)).
  !!@param   chi,out: irreducible susceptibility (full BZ, v_n>=0)
  !!@param     Gk,in: Green function on the full BZ (w_n>0)
  !!@param   kmap,in: integer mesh coordinates of each k-point
  !!@param  olist,in: orbital pair of each chi index
  !!@param   temp,in: Temperature
  !!@param     Nx,in: Number of kx mesh
  !!@param     Ny,in: Number of ky mesh
  !!@param     Nz,in: Number of kz mesh
  !!@param     Nw,in: Number of Matsubara frequencies
  !!@param  Nkall,in: Number of all k-points
  !!@param   Norb,in: Number of orbitals (spin included)
  !!@param   Nchi,in: Number of chi index
  use numerics, only:drop_roundoff
  use freq_grid
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nw,Norb,Nchi,Nkall,Nx,Ny,Nz
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap
  real(c_double),intent(in):: temp
  complex(c_double),intent(in),dimension(Nkall,Nw,Norb,Norb):: Gk
  complex(c_double),intent(out),dimension(Nkall,Nw,Nchi,Nchi):: chi

  integer(c_int32_t) i,j,k,l,m,n
  integer(c_int32_t) ii(0:Nx-1),ij(0:Ny-1),ik(0:Nz-1),iw(fg_nt(Nw))
  real(c_double) sg(fg_nt(Nw))   ! reflection sign of a fermionic factor
  integer(c_int64_t) jp,jm
  real(c_double) weight
  complex(c_double),allocatable,dimension(:,:,:,:):: tmp,tmpgk13,tmpgk42

  weight=temp/dble(Nkall)
  allocate(tmp(0:Nx-1,0:Ny-1,0:Nz-1,fg_nt(Nw)),tmpgk13(0:Nx-1,0:Ny-1,0:Nz-1,fg_nt(Nw)),&
       tmpgk42(0:Nx-1,0:Ny-1,0:Nz-1,fg_nt(Nw)))
  ! ii,ij,ik,iw: index of -r and -tau
  call fg_reverse(ii,ij,ik,iw,Nx,Ny,Nz,Nw,sg)
  call fg_plans_init(tmpgk13,tmp,Nx,Ny,Nz,Nw)
  l_loop:do l=1,Nchi
     m_loop:do m=1,l
        !$omp parallel do private(i,j,jp,jm)
        do j=1,Nw
           jp=fg_pos(FG_FERMION,j)
           jm=fg_neg(FG_FERMION,j,Nw)
           do i=1,Nkall
              tmpgk13(kmap(1,i),kmap(2,i),kmap(3,i),jp)=Gk(i,j,olist(l,1),olist(m,1)) !G13(k,iw)
              tmpgk42(kmap(1,i),kmap(2,i),kmap(3,i),jp)=Gk(i,j,olist(m,2),olist(l,2)) !G42(k,iw)
              tmpgk13(kmap(1,i),kmap(2,i),kmap(3,i),jm)=conjg(Gk(i,j,olist(m,1),olist(l,1))) !G13(k,-iw)=G^*31(k,iw)
              tmpgk42(kmap(1,i),kmap(2,i),kmap(3,i),jm)=conjg(Gk(i,j,olist(l,2),olist(m,2))) !G42(k,-iw)=G^*24(k,iw)
           end do
        end do
        !$omp end parallel do
        call fg_forward(FG_FERMION,tmpgk13,tmp,Nx,Ny,Nz,Nw)
        call fg_forward(FG_FERMION,tmpgk42,tmp,Nx,Ny,Nz,Nw)
        !calculate G(r)G(-r)
        !$omp parallel do private(i,j,k)
        do n=1,int(fg_nt(Nw),c_int32_t)
           do k=0,Nz-1
              do j=0,Ny-1
                 !$omp simd
                 do i=0,Nx-1
                    tmp(i,j,k,n)=-sg(n)*tmpgk13(i,j,k,n)*tmpgk42(ii(i),ij(j),ik(k),iw(n))
                 end do
                 !$omp end simd
              end do
           end do
        end do
        !$omp end parallel do
        call fg_backward(FG_BOSON,tmp,tmpgk13,Nx,Ny,Nz,Nw)
        !$omp parallel do private(i,j,jp,jm)
        do j=1,Nw
           jp=fg_pos(FG_BOSON,j)
           jm=fg_neg(FG_BOSON,j,Nw)   ! slot of -v_j (v_0 = 0 is its own partner: jm = jp)
           do i=1,Nkall
              chi(i,j,m,l)=tmp(kmap(1,i),kmap(2,i),kmap(3,i),jp)*weight
              chi(i,j,l,m)=conjg(tmp(kmap(1,i),kmap(2,i),kmap(3,i),jm))*weight !chi^lm(q,iv)=chi^*ml(q,-iv)
           end do
        end do
        !$omp end parallel do
        ! relative round-off cut on this (m,l) block (see drop_roundoff)
        call drop_roundoff(chi(1,1,m,l),Nkall*Nw)
        if(l/=m)call drop_roundoff(chi(1,1,l,m),Nkall*Nw)
     end do m_loop
  end do l_loop
  call fg_plans_free()
  deallocate(tmp,tmpgk13,tmpgk42)
end subroutine get_chi0_nosym

subroutine mkfk_nosym(fk,Gk,delta,ikq,ikm,Nkall,Nw,Norb)
  !> linearized anomalous green function without TRS for pairs (k+Q,-k)
  !> F(k,iw) = -G(k+Q,iw) Delta(k,iw) G^T(-k,-iw) = -G(k+Q,iw) Delta(k,iw) G^*(-k,iw)
  !!@param    fk,out: anomalous green function (w_n>0)
  !!@param     Gk,in: normal green function on the full BZ
  !!@param  delta,in: gap function
  !!@param    ikq,in: full-list index of k+Q
  !!@param    ikm,in: full-list index of -k
  !!@param  Nkall,in: Number of all k-points
  !!@param     Nw,in: Number of Matsubara frequencies
  !!@param   Norb,in: Number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nkall,Nw,Norb
  integer(c_int64_t),intent(in),dimension(Nkall):: ikq,ikm
  complex(c_double),intent(in),dimension(Nkall,Nw,Norb,Norb):: Gk,delta
  complex(c_double),intent(out),dimension(Nkall,Nw,Norb,Norb):: fk

  integer(c_int32_t) i,j
  complex(c_double),dimension(Norb,Norb):: cmat1,gmat2

  !$omp parallel do collapse(2) private(i,j,cmat1,gmat2)
  do j=1,Nw
     do i=1,Nkall
        gmat2(:,:)=conjg(Gk(ikm(i),j,:,:)) !G^T(-k,-iw)=G^*(-k,iw)
        ! cmat1 = G(k+Q,iw) * Delta(k,iw)
        call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),Gk(ikq(i),j,:,:),Norb,delta(i,j,:,:),Norb,(0.0d0,0.0d0),cmat1,Norb)
        ! fk = -cmat1 * G^*(-k,iw)
        call zgemm('N','N',Norb,Norb,Norb,(-1.0d0,0.0d0),cmat1,Norb,gmat2,Norb,(0.0d0,0.0d0),fk(i,j,:,:),Norb)
     end do
  end do
  !$omp end parallel do
end subroutine mkfk_nosym

subroutine mkdelta_nosym(newdelta,fk,Vdelta,Vmat,kmap,ikmq,olist,Nkall,Nw,Nchi,Norb,Nx,Ny,Nz)
  !> new gap function Delta_{l1 m2}(k) = sum_q V^{ml}(q) F_{l2 m1}(k-q)  (weight T/N applied by caller)
  !> F_{ab}(r,tau) is transformed once per orbital pair and reused for every chi pair,
  !> so one call costs Nchi^2 + 2*Norb^2 FFTs instead of 3*Nchi^2.  Pairs (l,m) whose V_Delta
  !> block and bare vertex vanish (e.g. spin-flip blocks when Sz is conserved) are skipped.
  !!@param newdelta,out: new gap function
  !!@param       fk,in: anomalous green function
  !!@param   Vdelta,in: pairing interaction V_Delta (full BZ, v_n>=0)
  !!@param     Vmat,in: bare vertex (w_n=>inf limit of V_Delta is -Vmat/2)
  !!@param     kmap,in: integer mesh coordinates of each k-point
  !!@param     ikmq,in: full-list index of -k-Q (Pauli partner of k)
  !!@param    olist,in: orbital pair of each chi index
  !!@param    Nkall,in: Number of all k-points
  !!@param       Nw,in: Number of Matsubara frequencies
  !!@param     Nchi,in: Number of chi index
  !!@param     Norb,in: Number of orbitals
  !!@param       Nx,in: Number of kx mesh
  !!@param       Ny,in: Number of ky mesh
  !!@param       Nz,in: Number of kz mesh
  use freq_grid
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nkall,Nw,Nchi,Norb,Nx,Ny,Nz
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap
  integer(c_int64_t),intent(in),dimension(Nkall):: ikmq
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  real(c_double),intent(in),dimension(Nchi,Nchi):: Vmat
  complex(c_double),intent(in),dimension(Nkall,Nw,Nchi,Nchi):: Vdelta
  complex(c_double),intent(in),dimension(Nkall,Nw,Norb,Norb):: fk
  complex(c_double),intent(out),dimension(Nkall,Nw,Norb,Norb):: newdelta

  integer(c_int32_t) i,j,k,n,l,m,a,b,c,d
  integer(c_int64_t) jp,jm
  integer(c_int32_t),dimension(Norb,Norb):: ichi
  real(c_double) vmax
  logical,dimension(Nchi,Nchi):: nzv
  complex(c_double),allocatable,dimension(:,:,:,:):: tmpVdelta,tmpdel,tmp
  complex(c_double),allocatable,dimension(:,:,:,:,:,:):: fr
  complex(c_double),dimension(Norb,Norb):: static_sum   ! IR: raw sum of F_ab over (k, w), for the bare vertex
  complex(c_double) static_corr

  ! ichi(a,b): chi index of orbital pair (a,b), 0 if the pair is not in olist
  ! nzv(m,l): block (m,l) of V_Delta contributes (relative cut far below FFT round-off)
  vmax=max(maxval(abs(Vdelta)),maxval(abs(Vmat)))
  !$omp parallel do private(l,m)
  do l=1,Nchi
     do m=1,Nchi
        nzv(m,l)=(abs(Vmat(m,l))>1.0d-14*vmax .or. maxval(abs(Vdelta(:,:,m,l)))>1.0d-14*vmax &
             .or. maxval(abs(Vdelta(:,:,l,m)))>1.0d-14*vmax)
     end do
  end do
  !$omp end parallel do
  ichi(:,:)=0
  do l=1,Nchi
     ichi(olist(l,1),olist(l,2))=l
  end do
  allocate(tmpVdelta(0:Nx-1,0:Ny-1,0:Nz-1,fg_nt(Nw)),tmpdel(0:Nx-1,0:Ny-1,0:Nz-1,fg_nt(Nw)),&
       tmp(0:Nx-1,0:Ny-1,0:Nz-1,fg_nt(Nw)))
  allocate(fr(0:Nx-1,0:Ny-1,0:Nz-1,fg_nt(Nw),Norb,Norb))
  call fg_plans_init(tmpVdelta,tmp,Nx,Ny,Nz,Nw)

  ! F_ab(r,tau) for every orbital pair
  do b=1,Norb
     do a=1,Norb
        !$omp parallel do private(i,j,jp,jm)
        do j=1,Nw
           jp=fg_pos(FG_FERMION,j)
           jm=fg_neg(FG_FERMION,j,Nw)
           do i=1,Nkall
              tmpdel(kmap(1,i),kmap(2,i),kmap(3,i),jp)=fk(i,j,a,b)
              tmpdel(kmap(1,i),kmap(2,i),kmap(3,i),jm)=-fk(ikmq(i),j,b,a) !F_ab(k,-iw)=-F_ba(-k-Q,iw)
           end do
        end do
        !$omp end parallel do
        static_sum(a,b)=fg_static_term(tmpdel,Nx,Ny,Nz,Nw)   ! 0 on the uniform grid
        call fg_forward(FG_FERMION,tmpdel,tmp,Nx,Ny,Nz,Nw)
        !$omp parallel workshare
        fr(:,:,:,:,a,b)=tmpdel(:,:,:,:)
        !$omp end parallel workshare
     end do
  end do

  ! Delta_ad(r,tau) = sum_{b,c} V^{ml}(r,tau) F_bc(r,tau),  l=(a,b), m=(c,d)
  d_loop:do d=1,Norb
     a_loop:do a=1,Norb
        !$omp parallel workshare
        tmpdel(:,:,:,:)=0.0d0
        !$omp end parallel workshare
        static_corr=0.0d0
        do b=1,Norb
           l=ichi(a,b)
           if(l==0)cycle
           do c=1,Norb
              m=ichi(c,d)
              if(m==0)cycle
              if(.not.nzv(m,l))cycle
              !$omp parallel
              !$omp do private(i,j,jp)
              do j=1,Nw
                 jp=fg_pos(FG_BOSON,j)
                 do i=1,Nkall
                    tmpVdelta(kmap(1,i),kmap(2,i),kmap(3,i),jp)=Vdelta(i,j,m,l)
                 end do
              end do
              !$omp end do
              !$omp do private(i,j,jm)
              do j=2,Nw
                 jm=fg_neg(FG_BOSON,j,Nw)
                 do i=1,Nkall
                    tmpVdelta(kmap(1,i),kmap(2,i),kmap(3,i),jm)=conjg(Vdelta(i,j,l,m)) !V^ml(q,-iv)=V^*lm(q,iv)
                 end do
              end do
              !$omp end do
              !$omp end parallel
              call fg_set_static(tmpVdelta,cmplx(-0.5d0*Vmat(m,l),0.0d0,kind=c_double),Nx,Ny,Nz,Nw) !bare -V/2
              static_corr=static_corr-0.5d0*Vmat(m,l)*static_sum(b,c)
              call fg_forward(FG_BOSON,tmpVdelta,tmp,Nx,Ny,Nz,Nw)
              !$omp parallel do private(i,j,k)
              do n=1,int(fg_nt(Nw),c_int32_t)
                 do k=0,Nz-1
                    do j=0,Ny-1
                       !$omp simd
                       do i=0,Nx-1
                          tmpdel(i,j,k,n)=tmpdel(i,j,k,n)+tmpVdelta(i,j,k,n)*fr(i,j,k,n,b,c)
                       end do
                       !$omp end simd
                    end do
                 end do
              end do
              !$omp end parallel do
           end do
        end do
        call fg_backward(FG_FERMION,tmpdel,tmp,Nx,Ny,Nz,Nw)
        !$omp parallel do private(i,j,jp)
        do j=1,Nw
           jp=fg_pos(FG_FERMION,j)
           do i=1,Nkall
              newdelta(i,j,a,d)=tmpdel(kmap(1,i),kmap(2,i),kmap(3,i),jp)
           end do
        end do
        !$omp end parallel do
        if(.not.fg_uniform())newdelta(:,:,a,d)=newdelta(:,:,a,d)+static_corr   ! IR: bare vertex kept out of the transform
     end do a_loop
  end do d_loop
  call fg_plans_free()
  deallocate(tmpVdelta,tmpdel,tmp,fr)
end subroutine mkdelta_nosym

subroutine get_initial_delta_nosym(delta,init_delta,Nkall,Nw,Norb)
  !> copy a frequency-independent orbital-basis seed to all w_n and normalize it
  !!@param      delta,out: gap function
  !!@param init_delta,in: orbital basis initial gap function on the full BZ
  !!@param      Nkall,in: Number of all k-points
  !!@param         Nw,in: Number of Matsubara frequencies
  !!@param       Norb,in: Number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nkall,Nw,Norb
  complex(c_double),intent(in),dimension(Nkall,Norb,Norb):: init_delta
  complex(c_double),intent(out),dimension(Nkall,Nw,Norb,Norb):: delta

  integer(c_int32_t) j
  real(c_double) norm

  !$omp parallel do private(j)
  do j=1,Nw
     delta(:,j,:,:)=init_delta(:,:,:)
  end do
  !$omp end parallel do
  norm=sqrt(2.0d0*sum(abs(delta(:,:,:,:))**2)) !factor 2: w_n<0 half
  if(norm <= 0.0d0)then
     print*,'Error: get_initial_delta_nosym norm <= 0'
     stop
  end if
  !$omp parallel workshare
  delta(:,:,:,:)=delta(:,:,:,:)/norm
  !$omp end parallel workshare
end subroutine get_initial_delta_nosym

subroutine get_vdelta_nosym(chi,Vmat,olist,slist,Nkall,Nw,Nchi,sw_dcsub) bind(C)
  !> interaction from chi0, in place (Q independent: compute once and reuse for every Q).
  !> Same two channels as the SOC code, since both use the same spin-orbital vertex Vmat:
  !>   sw_dcsub=.false. pairing     V_Delta = -V/2 + V chi_RPA          (the usual use, hence the name)
  !>   sw_dcsub=.true.  self-energy V_sigma = -V + V chi_RPA - V_od chi0 V_od
  !> The -V/2 vs -V and the double-counting subtraction are exactly the difference between
  !> lin_eliash_soc and mkself_soc; see get_V_soc_flex.  This solver itself only does pairing.
  !> The bare term carries 1/2: the antisymmetrized bare vertex is both the direct and the
  !> exchange term, and the gap equation sums over both orderings F_{l2 m1} of the pair, so
  !> -V alone would count it twice, whereas the RPA ladder V chi V is direct-channel only.
  !> With 1/2 the singlet channel of a spin-conserving system reduces to the standard
  !> V_s = (S+C)/2 + 3/2 S chi_s S - 1/2 C chi_c C (lin_eliash).  get_V_soc_flex itself keeps -V
  !> (the self-energy uses it); lin_eliash_soc applies the same 1/2 as here.
  !!@param chi,inout: in chi0 from get_chi0_nosym, out V_Delta
  !!@param   Vmat,in: V-matrix
  !!@param  Nkall,in: Number of all k-points
  !!@param     Nw,in: Number of Matsubara frequencies
  !!@param   Nchi,in: Number of chi index
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_bool
  implicit none
  integer(c_int64_t),intent(in):: Nkall,Nw,Nchi
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),dimension(*):: slist
  logical(c_bool),intent(in):: sw_dcsub
  real(c_double),intent(in),dimension(Nchi,Nchi):: Vmat
  complex(c_double),intent(inout),dimension(Nkall,Nw,Nchi,Nchi):: chi

  integer(c_int64_t) l,m
  logical(1) dcsub

  dcsub=logical(sw_dcsub,1)
  call get_V_soc_flex(chi,Vmat,olist,slist,Nkall,Nw,Nchi,dcsub)
  if(.not.dcsub)then   ! pairing only: -V -> -V/2, the bare term counted once
     !$omp parallel do private(l,m)
     do l=1,Nchi
        do m=1,Nchi
           chi(:,:,m,l)=chi(:,:,m,l)+0.5d0*Vmat(m,l)
        end do
     end do
     !$omp end parallel do
  end if
  print'(A15,2E16.8)','V max is ',maxval(dble(chi)),maxval(aimag(chi))
  print'(A15,2E16.8)','V min is ',minval(dble(chi)),minval(aimag(chi))
end subroutine get_vdelta_nosym

subroutine lin_eliash_nosym(delta,chi,Gk,init_delta,Vmat,olist,kmap,iQ,temp,eps,&
     Nkall,Nw,Nchi,Norb,Nx,Ny,Nz,itemax,arnoldi_m,lambda_out) bind(C)
  !> calculate linearized eliashberg equations without TRS (full BZ, spin-orbital basis)
  !> for pairs (k+Q,-k) with a fixed pair momentum Q on the k-mesh (iQ=0: uniform pairing)
  !> No singlet/triplet or parity projection is applied: the leading eigenvector of the
  !> full kernel is returned, and the channel follows from the seed and the physics.
  !!@param      delta,out: gap function (full BZ, w_n>0)
  !!@param  lambda_out,out: physical (largest real part) Eliashberg eigenvalue
  !!@param        chi,in: pairing interaction V_Delta from get_vdelta_nosym (not modified)
  !!@param         Gk,in: normal green function on the full BZ
  !!@param init_delta,in: orbital basis initial gap function on the full BZ
  !!@param       Vmat,in: V-matrix
  !!@param      olist,in: orbital pair of each chi index
  !!@param       kmap,in: integer mesh coordinates of each k-point
  !!@param         iQ,in: pair momentum Q in integer mesh units (Q=(iQ(1)/Nx,iQ(2)/Ny,iQ(3)/Nz))
  !!@param       temp,in: Temperature
  !!@param        eps,in: Threshold of calculation
  !!@param      Nkall,in: Number of all k-points
  !!@param         Nw,in: Number of Matsubara frequencies
  !!@param       Nchi,in: Number of chi index
  !!@param       Norb,in: Number of orbitals (spin included)
  !!@param         Nx,in: Number of kx mesh
  !!@param         Ny,in: Number of ky mesh
  !!@param         Nz,in: Number of kz mesh
  !!@param     itemax,in: maximum iteration of power method
  !!@param  arnoldi_m,in: Krylov subspace dimension (0=power method, >0=Arnoldi)
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nkall,Nw,Nchi,Norb,Nx,Ny,Nz,itemax,arnoldi_m
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap
  integer(c_int64_t),intent(in),dimension(3):: iQ
  real(c_double),intent(in):: temp,eps
  real(c_double),intent(in),dimension(Nchi,Nchi):: Vmat
  complex(c_double),intent(in),dimension(Nkall,Norb,Norb):: init_delta
  complex(c_double),intent(in),dimension(Nkall,Nw,Norb,Norb):: Gk
  complex(c_double),intent(out),dimension(Nkall,Nw,Norb,Norb):: delta
  complex(c_double),intent(in),dimension(Nkall,Nw,Nchi,Nchi):: chi
  real(c_double),intent(out):: lambda_out

  integer(c_int32_t) i_iter,i_eig
  integer(c_int32_t),parameter:: eig_max=2
  integer(c_int64_t),allocatable,dimension(:):: ikq,ikm,ikmq
  integer(c_int64_t),dimension(3),parameter:: iQ0=(/0_c_int64_t,0_c_int64_t,0_c_int64_t/)
  real(c_double) norm,inorm,weight,lambda_rq,lambda_prev,lambda_phys,lambda1,vec_err,proj
  complex(c_double),allocatable,dimension(:,:,:,:):: newdelta,fk,delta1

  allocate(ikq(Nkall),ikm(Nkall),ikmq(Nkall))
  allocate(newdelta(Nkall,Nw,Norb,Norb),fk(Nkall,Nw,Norb,Norb),delta1(Nkall,Nw,Norb,Norb))
  weight=temp/dble(Nkall)
  lambda1=0.0d0
  lambda_out=0.0d0
  call get_kshift_nosym(ikq,kmap,1_c_int64_t,iQ,Nkall,Nx,Ny,Nz)     !k+Q
  call get_kshift_nosym(ikm,kmap,-1_c_int64_t,iQ0,Nkall,Nx,Ny,Nz)   !-k
  call get_kshift_nosym(ikmq,kmap,-1_c_int64_t,-iQ,Nkall,Nx,Ny,Nz)  !-k-Q
  print'(A,3I5)','pair momentum iQ =',iQ
  if(arnoldi_m>0)then !--- Arnoldi solver
     call solve_arnoldi(arnoldi_m)
  else !--- Power method (shift + deflation)
     eigenval_loop:do i_eig=1,eig_max !solve eig_val using power method, 1st eig is usually large negative value
        call get_initial_delta_nosym(delta,init_delta,Nkall,Nw,Norb)
        if(i_eig==2)then !orthogonalize initial vector against 1st eigenvector
           call get_rayleigh(proj,delta1,delta)
           !$omp parallel workshare
           delta(:,:,:,:)=delta(:,:,:,:)-proj*delta1(:,:,:,:)
           !$omp end parallel workshare
           call get_norm(norm,delta)
           !$omp parallel workshare
           delta(:,:,:,:)=delta(:,:,:,:)*(1.0d0/norm)
           !$omp end parallel workshare
        end if
        lambda_prev=0.0d0
        iter_loop:do i_iter=1,itemax !iteration
           call apply_op(newdelta,delta)
           ! explicit orthogonal deflation + shift (i_eig==2):
           ! shift by -lambda1 (=|lambda1|>0) makes most-positive eigenvalue dominant
           if(i_eig==2)then
              !$omp parallel workshare
              newdelta(:,:,:,:)=newdelta(:,:,:,:)-lambda1*delta(:,:,:,:)
              !$omp end parallel workshare
              call get_rayleigh(proj,delta1,newdelta)
              !$omp parallel workshare
              newdelta(:,:,:,:)=newdelta(:,:,:,:)-proj*delta1(:,:,:,:)
              !$omp end parallel workshare
           end if
           call get_norm(norm,newdelta)
           if(norm <= 0.0d0)then
              print*,'Error: norm <= 0 in iter_loop'
              stop
           end if
           inorm=1.0d0/norm
           call get_rayleigh(lambda_rq,delta,newdelta)
           lambda_phys=lambda_rq+lambda1
           call get_vec_err(vec_err,newdelta,delta,inorm)
           if(abs(lambda_phys)>=1.0d2 .or. abs(lambda_phys)<1.0d-6)then
              print'(I4,A13,2E16.8,A9,E12.4)',i_iter,' lambda_elsh=',lambda_phys,abs(lambda_phys-lambda_prev),'  vec_err=',vec_err
           else
              print'(I4,A13,2F12.8,A9,E12.4)',i_iter,' lambda_elsh=',lambda_phys,abs(lambda_phys-lambda_prev),'  vec_err=',vec_err
           end if
           if(lambda_phys<0.0d0 .and. vec_err<eps)exit
           if(vec_err<eps .and. abs(lambda_phys-lambda_prev)<eps*(abs(lambda_phys)+1.0d-10))exit
           lambda_prev=lambda_phys
           !$omp parallel workshare
           delta(:,:,:,:)=newdelta(:,:,:,:)*inorm
           !$omp end parallel workshare
        end do iter_loop
        if(i_eig==1)then
           if(lambda_phys>0.0d0)then
              print*,'1st eigenvalue is positive: skipping 2nd loop'
              print*,'eliash=',lambda_phys
              lambda_out=lambda_phys
              exit
           end if
           lambda1=lambda_rq !save first eigenvalue for shift in i_eig=2
           delta1(:,:,:,:)=newdelta(:,:,:,:)*inorm !store 1st eigenvector (normalized) for deflation
        else
           print*,'eliash=',lambda_phys
           lambda_out=lambda_phys
        end if
     end do eigenval_loop
     call get_norm(norm,newdelta)
     if(norm <= 0.0d0)then
        print*,'Error: final norm <= 0'
        stop
     end if
     inorm=1.0d0/norm
     !$omp parallel workshare
     delta(:,:,:,:)=newdelta(:,:,:,:)*inorm
     !$omp end parallel workshare
  end if
  deallocate(newdelta,fk,delta1)
contains
  subroutine get_norm(norm,func)
    !> norm over the full BZ; factor 2 for the w_n<0 half
    complex(c_double),intent(in),dimension(Nkall,Nw,Norb,Norb):: func
    real(c_double),intent(out):: norm

    integer(c_int32_t) i,j,l,m
    real(c_double) tmp

    tmp=0.0d0
    !$omp parallel
    do l=1,Norb
       do m=1,Norb
          !$omp do private(i,j),reduction(+:tmp)
          do j=1,Nw
             do i=1,Nkall
                tmp=tmp+abs(func(i,j,m,l))*abs(func(i,j,m,l))
             end do
          end do
          !$omp end do
       end do
    end do
    !$omp end parallel
    norm=sqrt(2.0d0*tmp)
  end subroutine get_norm

  subroutine get_rayleigh(rq,del,newdel)
    !> Rayleigh quotient: rq = Re(<del, newdel>) assuming ||del||=1
    complex(c_double),intent(in),dimension(Nkall,Nw,Norb,Norb):: del,newdel
    real(c_double),intent(out):: rq

    complex(c_double) h
    call get_inner(h,del,newdel)
    rq=dble(h)
  end subroutine get_rayleigh

  subroutine get_vec_err(verr,newdel,del,inrm)
    !> Vector convergence: min(||newdel*inrm - del||, ||newdel*inrm + del||)
    !> Taking the minimum handles sign-flipping convergence for negative eigenvalues
    complex(c_double),intent(in),dimension(Nkall,Nw,Norb,Norb):: newdel,del
    real(c_double),intent(in):: inrm
    real(c_double),intent(out):: verr

    integer(c_int32_t) i,j,l,m
    real(c_double) tmp_pos,tmp_neg

    tmp_pos=0.0d0
    tmp_neg=0.0d0
    !$omp parallel
    do l=1,Norb
       do m=1,Norb
          !$omp do private(i,j),reduction(+:tmp_pos,tmp_neg)
          do j=1,Nw
             do i=1,Nkall
                tmp_pos=tmp_pos+abs(newdel(i,j,m,l)*inrm-del(i,j,m,l))**2
                tmp_neg=tmp_neg+abs(newdel(i,j,m,l)*inrm+del(i,j,m,l))**2
             end do
          end do
          !$omp end do
       end do
    end do
    !$omp end parallel
    verr=sqrt(2.0d0*min(tmp_pos,tmp_neg))
  end subroutine get_vec_err

  subroutine solve_arnoldi(m_arnoldi)
    !> Two-pass Arnoldi: pass 1 (no shift) finds lambda_-, pass 2 (shift+deflation) finds lambda_+.
    !> Same algorithm as lin_eliash; see the comments there.  Without TRS the kernel is in
    !> general not similar to a real matrix, so Ritz values are printed with their imaginary part.
    integer(c_int64_t),intent(in):: m_arnoldi
    integer(c_int32_t) :: m_dim,j,ii,m_act,idx,lwork,info
    real(c_double) :: beta,shift
    complex(c_double) :: hij
    complex(c_double),allocatable :: V(:,:,:,:,:)
    complex(c_double),allocatable :: H_mat(:,:),A_mat(:,:),eigenvals(:)
    complex(c_double),allocatable :: VL_dum(:,:),VR(:,:),zwork(:)
    real(c_double),allocatable :: rwork(:)

    m_dim=int(m_arnoldi,c_int32_t)
    lwork=max(64,2*m_dim)
    allocate(V(Nkall,Nw,Norb,Norb,0:m_dim-1),H_mat(m_dim+1,m_dim))
    allocate(A_mat(m_dim,m_dim),eigenvals(m_dim))
    allocate(VL_dum(1,1),VR(m_dim,m_dim),zwork(lwork),rwork(2*m_dim))

    ! ===== Pass 1: no shift, find lambda_- =====================================
    H_mat=(0.0d0,0.0d0)
    call get_initial_delta_nosym(V(:,:,:,:,0),init_delta,Nkall,Nw,Norb)
    print'(A,I4,A)','Arnoldi pass 1 (no shift): Krylov dim =',m_dim,' ...'
    m_act=m_dim
    do j=0,m_dim-1
       call apply_op(newdelta,V(:,:,:,:,j))                ! w = A v_j
       do ii=0,j                                            ! modified Gram-Schmidt
          call get_inner(hij,V(:,:,:,:,ii),newdelta)
          H_mat(ii+1,j+1)=hij
          !$omp parallel workshare
          newdelta(:,:,:,:)=newdelta(:,:,:,:)-hij*V(:,:,:,:,ii)
          !$omp end parallel workshare
       end do
       call get_norm(beta,newdelta)
       H_mat(j+2,j+1)=cmplx(beta,0.0d0,c_double)
       if(beta<1.0d-14)then; m_act=j+1; exit; end if      ! lucky breakdown
       if(j<m_dim-1)then
          !$omp parallel workshare
          V(:,:,:,:,j+1)=newdelta(:,:,:,:)*(1.0d0/beta)
          !$omp end parallel workshare
       end if
    end do
    A_mat(1:m_act,1:m_act)=H_mat(1:m_act,1:m_act)
    call zgeev('N','V',m_act,A_mat,m_dim,eigenvals,VL_dum,1,VR,m_dim,zwork,lwork,rwork,info)
    if(info/=0)then; print*,'ZGEEV failed pass 1: info=',info; stop; end if
    print'(A)','Pass 1 Ritz values (real, imag):'
    do ii=1,m_act
       print'(I4,2F14.6)',ii,dble(eigenvals(ii)),aimag(eigenvals(ii))
    end do
    idx=minloc(dble(eigenvals(1:m_act)),1)                  ! most negative Ritz value
    lambda1=dble(eigenvals(idx))
    print'(A,F12.6)','  lambda_- =',lambda1
    if(maxval(dble(eigenvals(1:m_act)))>=0.1d0)then         ! lambda_+ already resolved
       idx=maxloc(dble(eigenvals(1:m_act)),1)
       print'(A,2F12.6)','  eliash   =',dble(eigenvals(idx)),aimag(eigenvals(idx))
       lambda_out=dble(eigenvals(idx))
       call ritz_vector(delta,V,VR(:,idx),m_act,m_dim)
       deallocate(V,H_mat,A_mat,eigenvals,VL_dum,VR,zwork,rwork)
       return
    end if
    call ritz_vector(delta1,V,VR(:,idx),m_act,m_dim) ! deflation vector for pass 2

    ! ===== Pass 2: shift + deflation, find lambda_+ ============================
    shift=-lambda1
    H_mat=(0.0d0,0.0d0)
    call get_initial_delta_nosym(V(:,:,:,:,0),init_delta,Nkall,Nw,Norb)
    call get_inner(hij,delta1,V(:,:,:,:,0))
    !$omp parallel workshare
    V(:,:,:,:,0)=V(:,:,:,:,0)-hij*delta1(:,:,:,:)
    !$omp end parallel workshare
    call get_norm(beta,V(:,:,:,:,0))
    !$omp parallel workshare
    V(:,:,:,:,0)=V(:,:,:,:,0)*(1.0d0/beta)
    !$omp end parallel workshare
    print'(A,F10.6,A,I4,A)','Arnoldi pass 2 (shift=',shift,'): Krylov dim =',m_dim,' ...'
    m_act=m_dim
    do j=0,m_dim-1
       call apply_op(newdelta,V(:,:,:,:,j))
       !$omp parallel workshare
       newdelta(:,:,:,:)=newdelta(:,:,:,:)+shift*V(:,:,:,:,j)   !shift
       !$omp end parallel workshare
       call get_inner(hij,delta1,newdelta)                       !deflation
       !$omp parallel workshare
       newdelta(:,:,:,:)=newdelta(:,:,:,:)-hij*delta1(:,:,:,:)
       !$omp end parallel workshare
       do ii=0,j  !Gram-Schmidt
          call get_inner(hij,V(:,:,:,:,ii),newdelta)
          H_mat(ii+1,j+1)=hij
          !$omp parallel workshare
          newdelta(:,:,:,:)=newdelta(:,:,:,:)-hij*V(:,:,:,:,ii)
          !$omp end parallel workshare
       end do
       call get_norm(beta,newdelta)
       H_mat(j+2,j+1)=cmplx(beta,0.0d0,c_double)
       if(beta<1.0d-14)then; m_act=j+1; exit; end if
       if(j<m_dim-1)then
          !$omp parallel workshare
          V(:,:,:,:,j+1)=newdelta(:,:,:,:)*(1.0d0/beta)
          !$omp end parallel workshare
       end if
    end do
    A_mat(1:m_act,1:m_act)=H_mat(1:m_act,1:m_act)
    call zgeev('N','V',m_act,A_mat,m_dim,eigenvals,VL_dum,1,VR,m_dim,zwork,lwork,rwork,info)
    if(info/=0)then; print*,'ZGEEV failed pass 2: info=',info; stop; end if
    print'(A)','Pass 2 Ritz values - physical (real, imag):'
    do ii=1,m_act
       print'(I4,2F14.6)',ii,dble(eigenvals(ii))-shift,aimag(eigenvals(ii))
    end do
    idx=maxloc(dble(eigenvals(1:m_act)),1)
    print'(A,2F12.6)','  eliash   =',dble(eigenvals(idx))-shift,aimag(eigenvals(idx))
    lambda_out=dble(eigenvals(idx))-shift
    call ritz_vector(delta,V,VR(:,idx),m_act,m_dim)
    deallocate(V,H_mat,A_mat,eigenvals,VL_dum,VR,zwork,rwork)
  end subroutine solve_arnoldi

  subroutine ritz_vector(x,V,y,m_act,m_dim)
    !> x = V_m * y, normalized (Ritz vector from Arnoldi basis V and eigenvector y of H_m)
    integer(c_int32_t),intent(in):: m_act,m_dim
    complex(c_double),intent(in),dimension(Nkall,Nw,Norb,Norb,0:m_dim-1):: V
    complex(c_double),intent(in),dimension(m_dim):: y
    complex(c_double),intent(out),dimension(Nkall,Nw,Norb,Norb):: x
    integer(c_int32_t) jj,l
    real(c_double) nrm
    !$omp parallel workshare
    x(:,:,:,:)=(0.0d0,0.0d0)
    !$omp end parallel workshare
    do jj=1,m_act
       !$omp parallel do private(l)
       do l=1,Norb
          x(:,:,:,l)=x(:,:,:,l)+y(jj)*V(:,:,:,l,jj-1)
       end do
       !$omp end parallel do
    end do
    call get_norm(nrm,x)
    !$omp parallel workshare
    x(:,:,:,:)=x(:,:,:,:)*(1.0d0/nrm)
    !$omp end parallel workshare
  end subroutine ritz_vector

  subroutine apply_op(w_out,v_in)
    !> apply Eliashberg operator: w_out = -(T/Nkall)*K*v_in; uses host fk as scratch
    !> (overall minus as in lin_eliash_soc, matching the SOC vertex V_Delta=-V+V*chi)
    complex(c_double),intent(in),dimension(Nkall,Nw,Norb,Norb):: v_in
    complex(c_double),intent(out),dimension(Nkall,Nw,Norb,Norb):: w_out
    call mkfk_nosym(fk,Gk,v_in,ikq,ikm,Nkall,Nw,Norb)
    call mkdelta_nosym(w_out,fk,chi,Vmat,kmap,ikmq,olist,Nkall,Nw,Nchi,Norb,Nx,Ny,Nz)
    !$omp parallel workshare
    w_out(:,:,:,:)=-w_out(:,:,:,:)*weight
    !$omp end parallel workshare
  end subroutine apply_op

  subroutine get_inner(h,u,v)
    !> complex inner product h = 2*<u,v> over the full BZ (factor 2 for the w_n<0 half)
    complex(c_double),intent(in),dimension(Nkall,Nw,Norb,Norb):: u,v
    complex(c_double),intent(out):: h
    integer(c_int32_t) i,j,l,m
    real(c_double) tmp_r,tmp_i
    tmp_r=0.0d0
    tmp_i=0.0d0
    !$omp parallel
    do l=1,Norb
       do m=1,Norb
          !$omp do private(i,j),reduction(+:tmp_r,tmp_i)
          do j=1,Nw
             do i=1,Nkall
                tmp_r=tmp_r+dble(conjg(u(i,j,m,l))*v(i,j,m,l))
                tmp_i=tmp_i+aimag(conjg(u(i,j,m,l))*v(i,j,m,l))
             end do
          end do
          !$omp end do
       end do
    end do
    !$omp end parallel
    h=cmplx(2.0d0*tmp_r,2.0d0*tmp_i,c_double)
  end subroutine get_inner
end subroutine lin_eliash_nosym

subroutine conv_delta_orb_to_band_nosym(deltab,delta,uni,kmap,iQ,Norb,Nkall,Nw,Nx,Ny,Nz) bind(C)
  !> transform the lowest-w_n orbital basis gap into the band basis without TRS
  !> Delta_band(k) = U^+(k+Q) Delta_orb(k) U^*(-k): row n = band n at k+Q, column n' = band n' at -k
  !!@param deltab,out: band basis gap function
  !!@param   delta,in: orbital basis gap function (full BZ)
  !!@param     uni,in: unitary matrix on the full BZ
  !!@param    kmap,in: integer mesh coordinates of each k-point
  !!@param      iQ,in: pair momentum Q in integer mesh units
  !!@param    Norb,in: Number of orbitals
  !!@param   Nkall,in: Number of all k-points
  !!@param      Nw,in: Number of Matsubara frequencies
  !!@param      Nx,in: Number of kx mesh
  !!@param      Ny,in: Number of ky mesh
  !!@param      Nz,in: Number of kz mesh
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Norb,Nkall,Nw,Nx,Ny,Nz
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap
  complex(c_double),intent(in),dimension(Norb,Norb,Nkall):: uni
  complex(c_double),intent(in),dimension(Nkall,Nw,Norb,Norb):: delta
  complex(c_double),intent(out),dimension(Nkall,Norb,Norb):: deltab

  integer(c_int32_t) i
  integer(c_int64_t),intent(in),dimension(3):: iQ
  integer(c_int64_t),allocatable,dimension(:):: ikq,ikm
  integer(c_int64_t),dimension(3),parameter:: iQ0=(/0_c_int64_t,0_c_int64_t,0_c_int64_t/)
  complex(c_double),dimension(Norb,Norb):: matL,matR,tmpm

  allocate(ikq(Nkall),ikm(Nkall))
  call get_kshift_nosym(ikq,kmap,1_c_int64_t,iQ,Nkall,Nx,Ny,Nz)
  call get_kshift_nosym(ikm,kmap,-1_c_int64_t,iQ0,Nkall,Nx,Ny,Nz)
  !$omp parallel do private(i,matL,matR,tmpm)
  do i=1,Nkall
     matL=conjg(transpose(uni(:,:,ikq(i))))
     matR=conjg(uni(:,:,ikm(i)))
     call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),matL,Norb,delta(i,1,:,:),Norb,(0.0d0,0.0d0),tmpm,Norb)
     call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),tmpm,Norb,matR,Norb,(0.0d0,0.0d0),deltab(i,:,:),Norb)
  end do
  !$omp end parallel do
end subroutine conv_delta_orb_to_band_nosym
