! Linearized Eliashberg equation WITHOUT spin-orbit coupling, using TRS and inversion.
!
! Spin is not an index here: the singlet and triplet channels are separate kernels, picked by the
! sign of gap_sym (sw_pair), and differ in the pairing vertex and in the sign sgn of the Pauli
! relation F(-k,-w) = sgn F^T(k,w)  (sgn = +1 singlet, -1 triplet).
!
! Storage: Gk, chi and Delta live on the Nk irreducible TRS k-points, with w_n > 0 only.
!   -k   comes from inversion parity,  Delta(-k,w) = dprt P Delta(k,w) P  (dprt = +1 even, -1 odd),
!   -w_n from combining that with Pauli, Delta(k,-w) = sgn dprt P Delta^T(k,w) P.
!   P is a SIGNED PERMUTATION of the orbitals: prt(a) = +-1 is the parity of orbital a and iprm(a)
!   is the orbital inversion sends it to (identity for one site, plibs.inversion_op otherwise).
!   k = -k points (TRIM) are stored once, so parity cannot rebuild them and is imposed explicitly
!   at the end of mkdelta_nsoc; without that an odd-parity (triplet) gap keeps a parity-violating
!   component there, which gives a wrong eigenvalue in multi-orbital models.
!
! Real-linear kernel: mkdelta_nsoc computes the lower triangle of Delta and fills the upper one by
! conjugation (Delta is Hermitian in the TRS gauge), so apply_op is linear over the REALS only.
! Arnoldi therefore runs with a real inner product (get_inner) and real Ritz coefficients
! (real_ritz); a complex Krylov space mixes the two and silently returns wrong eigenvalues.

module lin_eliash_op
  !> The kernel of lin_eliash for the shared Arnoldi driver (eliash_krylov's arnoldi_two_pass).
  !> It is passed to the driver as an argument, so it is a module procedure reading state that
  !> lin_eliash points at, not an internal procedure of lin_eliash: passing an internal procedure
  !> needs a trampoline on the stack, i.e. an executable stack, which ld.lld refuses to link (LLVM
  !> flang) and glibc >= 2.41 refuses to dlopen.
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double
  implicit none
  integer(c_int64_t):: op_Nkall,op_Nk,op_Nw,op_Nchi,op_Norb,op_Nx,op_Ny,op_Nz,op_nvec
  logical(1):: op_sw_pair
  real(c_double):: op_weight
  integer(c_int64_t),pointer,contiguous:: op_olist(:,:)=>null(),op_kmap(:,:)=>null(),&
       op_invk(:,:)=>null(),op_iprm(:)=>null()
  real(c_double),pointer,contiguous:: op_prt(:)=>null(),op_Smat(:,:)=>null(),op_Cmat(:,:)=>null()
  complex(c_double),pointer,contiguous:: op_Gk(:,:,:,:)=>null(),op_chi(:,:,:,:)=>null(),&
       op_fk(:,:,:,:)=>null()
contains
  subroutine apply_op(w_out,v_in)
    !> one application of the Eliashberg kernel, w_out = (T/Nkall)*K*v_in; uses op_fk as scratch.
    !> Flat dummies: the shared driver stores Krylov vectors as plain vectors, and the rank-4
    !> kernels below pick up the (Nk,Nw,Norb,Norb) shape by sequence association.
    complex(c_double),intent(in):: v_in(*)
    complex(c_double),intent(out):: w_out(*)
    integer(c_int64_t) q
    call mkfk_trs_nsoc(op_fk,op_Gk,v_in,op_Nk,op_Nw,op_Norb)
    call mkdelta_nsoc(w_out,op_fk,op_chi,op_Smat,op_Cmat,op_kmap,op_invk,op_prt,op_iprm,op_olist,&
         op_Nkall,op_Nk,op_Nw,op_Nchi,op_Norb,op_Nx,op_Ny,op_Nz,op_sw_pair)
    !$omp parallel do private(q)
    do q=1,op_nvec
       w_out(q)=w_out(q)*op_weight
    end do
    !$omp end parallel do
  end subroutine apply_op
end module lin_eliash_op

subroutine lin_eliash(delta,chi,Gk,uni,init_delta,Smat,Cmat,vph,olist,prt,iprm,kmap,invk,temp,eps,&
     Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,itemax,gap_sym,arnoldi_m,lambda_out) bind(C)
  !> calculate linearized eliashberg equations with TRS without soc
  !!@param     delta,out: gap function
  !!@param   lambda_out,out: physical (largest) Eliashberg eigenvalue
  !!@param         Gk,in: normal green function
  !!@note: invk(:,i) index validity check:
  !!        - invk(1,i) must be in range [1, Nk]
  !!        - Enable bounds checking: ifx -check=bounds or gfortran -fbounds-check
  !!@param        uni,in: unitary matrix
  !!@param init_delta,in: band basis initial gap function
  !!@param       Smat,in: S-matrix
  !!@param       Cmat,in: C-matrix
  !!@param        vph,in: k-independent pairing term added to V_Delta(q,iv_j) [Nw,Nchi,Nchi] in
  !!                       the singlet/triplet channel gap_sym selects (j=1: v=0), e.g. the
  !!                       phonon-mediated attraction of plibs.epa_pairing_vertex; 0 = none
  !!@param      olist,in: property of chi footnote
  !!@param       kmap,in: property of k-points
  !!@param       invk,in: list of reverse k-points
  !!@param       temp,in: Temperature
  !!@param        eps,in: Threshold of calculation
  !!@param      Nkall,in: Number of all k-points
  !!@param         Nk,in: Number of k-points
  !!@param         Nw,in: Number of Matsubara frequencies
  !!@param       Nchi,in: Number of footnote of chi
  !!@param       Norb,in: Number of orbitals
  !!@param         Nx,in: Number of kx mesh
  !!@param         Ny,in: Number of ky mesh
  !!@param         Nz,in: Number of kz mesh
  !!@param     itemax,in: maximum iteration of power method
  !!@param    gap_sym,in: gap symmetry number.  Its SIGN picks the channel and the parity that
  !!                      goes with it: >0 singlet (even parity), <0 triplet (odd parity); its
  !!                      magnitude picks the band-basis seed (get_initial_delta), 0 = unit matrix.
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  use eliash_krylov
  use lin_eliash_op
  use grid_work, only:grid_work_free
  implicit none
  integer(c_int64_t),intent(in):: Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,itemax,gap_sym
  integer(c_int64_t),intent(in):: arnoldi_m  !! Krylov subspace dimension (0=power method, >0=Arnoldi)
  integer(c_int64_t),intent(in),target,dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),target,dimension(3,Nkall):: kmap,invk
  !! Inversion in the orbital basis is the signed permutation P_{m,iprm(m)}=prt(m):
  !! Delta_ml(-k) = dprt*prt(m)*prt(l)*Delta_{iprm(m),iprm(l)}(k).
  !! iprm=identity recovers the old diagonal-parity behaviour bit for bit.
  integer(c_int64_t),intent(in),target,dimension(Norb):: iprm
  real(c_double),intent(in):: temp,eps
  real(c_double),intent(in),target,dimension(Norb):: prt
  real(c_double),intent(in),target,dimension(Nchi,Nchi):: Smat,Cmat
  real(c_double),intent(in),dimension(Nw,Nchi,Nchi):: vph
  real(c_double),intent(in),dimension(Nk,Norb):: init_delta
  complex(c_double),intent(in),target,dimension(Nk,Nw,Norb,Norb):: Gk
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(c_double),intent(out),dimension(Nk,Nw,Norb,Norb):: delta
  complex(c_double),intent(inout),target,dimension(Nk,Nw,Nchi,Nchi):: chi
  real(c_double),intent(out):: lambda_out

  integer(c_int32_t) i_iter,i_eig,i
  integer(c_int64_t) ivj,ivm,ivl
  integer(c_int32_t),parameter:: eig_max=2
  logical(1) sw_pair
  real(c_double) norm,inorm,weight,lambda_rq,lambda_prev,lambda_phys,lambda1,vec_err,proj
  complex(c_double),allocatable,dimension(:,:,:,:):: newdelta,delta1
  complex(c_double),allocatable,target,dimension(:,:,:,:):: fk
  !> vector geometry for the shared Krylov routines: Delta is kept on the Nk irreducible
  !> k-points, and invk(1,i) says which of them carries full-BZ point i
  type(kspace) ks

  allocate(newdelta(Nk,Nw,Norb,Norb),fk(Nk,Nw,Norb,Norb),delta1(Nk,Nw,Norb,Norb))
  if(gap_sym>=0)then
     sw_pair=.true.
     print'(A7)','singlet'
  else
     sw_pair=.false.
     print'(A7)','triplet'
  end if
  weight=temp/dble(Nkall)
  call kspace_init(ks,invk(1,:),Nkall,Nk,Nw,Norb)
  lambda1=0.0d0
  lambda_out=0.0d0
  call get_V_delta_nsoc_flex(chi,Smat,Cmat,Nk,Nw,Nchi,sw_pair)
  if(any(vph/=0.0d0))then !extra k-independent pairing (e.g. phonons), vanishing as v -> inf
     !$omp parallel do collapse(2) private(ivm,ivl,ivj)
     do ivl=1,Nchi
        do ivm=1,Nchi
           do ivj=1,Nw
              chi(:,ivj,ivm,ivl)=chi(:,ivj,ivm,ivl)+vph(ivj,ivm,ivl)
           end do
        end do
     end do
     !$omp end parallel do
  end if
  print'(A15,2E16.8)','V_delta max is ',maxval(dble(chi)),maxval(aimag(chi))
  print'(A15,2E16.8)','V_delta min is ',minval(dble(chi)),minval(aimag(chi))
  if(arnoldi_m>0)then !--- Arnoldi solver (shared with lin_eliash_soc, see module eliash_krylov)
     ! no branch test is passed: singlet and triplet are separate kernels here (sw_pair), so the
     ! largest eigenvalue of this kernel is already the one gap_sym asks for
     call get_initial_delta(newdelta,init_delta,uni,kmap,invk,Nkall,Nk,Nw,Norb,gap_sym)
     ! state of apply_op (module lin_eliash_op)
     op_Nkall=Nkall; op_Nk=Nk; op_Nw=Nw; op_Nchi=Nchi; op_Norb=Norb; op_Nx=Nx; op_Ny=Ny; op_Nz=Nz
     op_nvec=ks%nvec; op_sw_pair=sw_pair; op_weight=weight
     op_olist=>olist; op_kmap=>kmap; op_invk=>invk; op_iprm=>iprm; op_prt=>prt
     op_Smat=>Smat; op_Cmat=>Cmat; op_Gk=>Gk; op_chi=>chi; op_fk=>fk
     call arnoldi_two_pass(lambda_out,delta,newdelta,ks,arnoldi_m,eps,apply_op)
     nullify(op_olist,op_kmap,op_invk,op_iprm,op_prt,op_Smat,op_Cmat,op_Gk,op_chi,op_fk)
  else !--- Power method (shift + deflation)
     eigenval_loop:do i_eig=1,eig_max !solve eig_val using power method, 1st eig is usually large negative value
        call get_initial_delta(delta,init_delta,uni,kmap,invk,Nkall,Nk,Nw,Norb,gap_sym)
        if(i_eig==2)then !orthogonalize initial vector against 1st eigenvector
           proj=kv_dot(ks,delta1,delta)
           !$omp parallel workshare
           delta(:,:,:,:)=delta(:,:,:,:)-proj*delta1(:,:,:,:)
           !$omp end parallel workshare
           norm=kv_norm(ks,delta)
           !$omp parallel workshare
           delta(:,:,:,:)=delta(:,:,:,:)*(1.0d0/norm)
           !$omp end parallel workshare
        end if
        lambda_prev=0.0d0
        iter_loop:do i_iter=1,itemax !iteration
           call mkfk_trs_nsoc(fk,Gk,delta,Nk,Nw,Norb)
           call mkdelta_nsoc(newdelta,fk,chi,Smat,Cmat,kmap,invk,prt,iprm,olist,Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,sw_pair)
           !$omp parallel workshare
           newdelta(:,:,:,:)=newdelta(:,:,:,:)*weight
           !$omp end parallel workshare
           ! explicit orthogonal deflation + shift (i_eig==2):
           ! shift by -lambda1 (=|lambda1|>0) makes most-positive eigenvalue dominant
           if(i_eig==2)then
              !$omp parallel workshare
              newdelta(:,:,:,:)=newdelta(:,:,:,:)-lambda1*delta(:,:,:,:)
              !$omp end parallel workshare
              proj=kv_dot(ks,delta1,newdelta)
              !$omp parallel workshare
              newdelta(:,:,:,:)=newdelta(:,:,:,:)-proj*delta1(:,:,:,:)
              !$omp end parallel workshare
           end if
           norm=kv_norm(ks,newdelta)
           if(norm <= 0.0d0)then
              print*,'Error: norm <= 0 in iter_loop'
              stop
           end if
           inorm=1.0d0/norm
           lambda_rq=kv_dot(ks,delta,newdelta)
           lambda_phys=lambda_rq+lambda1
           vec_err=kv_err(ks,newdelta,delta,inorm)
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
     norm=kv_norm(ks,newdelta)
     if(norm <= 0.0d0)then
        print*,'Error: final norm <= 0'
        stop
     end if
     inorm=1.0d0/norm
     !$omp parallel workshare
     delta(:,:,:,:)=newdelta(:,:,:,:)*inorm
     !$omp end parallel workshare
  end if
  call grid_work_free()   ! mkdelta_nsoc's work arrays, kept across the iterations
end subroutine lin_eliash

subroutine mkfk_trs_nsoc(fk,Gk,delta,Nk,Nw,Norb)
  !>calculate linearized anomalous green function Fk with TRS
  !>if we considder TRS,F_ab(k)=G_ac(k)Delta_cd(k)Gbd(-k)
  !> =G_ac(k)Delta_cd(k)G^*_dbss'(k)ss'
  !!@param   fk,out: anomalous green function
  !!@param    Gk,in: normal green function
  !!@param delta,in: gap function
  !!@param  invk,in: list of reverse k-points
  !!@param Nkall,in: Number of all k-points
  !!@param    Nk,in: Number of irreducible k-points
  !!@param    Nw,in: Number of Matsubara frequencies
  !!@param  Norb,in: Number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,Norb
  complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk,delta
  complex(c_double),intent(out),dimension(Nk,Nw,Norb,Norb):: fk
 
  integer(c_int32_t) i,j,l,m
  complex(c_double),dimension(Norb,Norb):: cmat1

  !$omp parallel do collapse(2) private(i,j,l,m,cmat1)
  do j=1,Nw
     do i=1,Nk
        ! cmat1 = G(k,iw) * Delta(k,iw)
        call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),Gk(i,j,:,:),Norb,delta(i,j,:,:),Norb,(0.0d0,0.0d0),cmat1,Norb)
        ! fk = -cmat1 * G(k,iw)^H
        call zgemm('N','C',Norb,Norb,Norb,(-1.0d0,0.0d0),cmat1,Norb,Gk(i,j,:,:),Norb,(0.0d0,0.0d0),fk(i,j,:,:),Norb)

        ! enforce Hermitian form used in the original implementation
        do l=1,Norb
           fk(i,j,l,l)=dble(fk(i,j,l,l))
           do m=l+1,Norb
              fk(i,j,l,m)=conjg(fk(i,j,m,l))
           end do
        end do
     end do
  end do
  !$omp end parallel do
end subroutine mkfk_trs_nsoc

subroutine mkdelta_nsoc(newdelta,delta,Vdelta,Smat,Cmat,kmap,invk,prt,iprm,olist,Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz,sw_pair)
  !> This function obtain new gap function Delta=sum VF
  !!@param newdelta,out: new gap function
  !!@param     delta,in: anomalous green function
  !!@param    Vdelta,in: irreducible susceptibility and pairing interaction
  !!@param      Smat,in: S-matrix
  !!@param      Cmat,in: C-matrix
  !!@param      kmap,in: property of k-point
  !!@param      invk,in: list of reverse k-points
  !!@param     olist,in: property of chi footnote
  !!@param     Nkall,in: Number of all k-points
  !!@param        Nk,in: Number of k-points
  !!@param        Nw,in: Number of Matsubara frequencies
  !!@param      Nchi,in: Number of footnote of chi
  !!@param      Norb,in: Number of orbitals
  !!@param        Nx,in: Number of kx mesh
  !!@param        Ny,in: Number of ky mesh
  !!@param        Nz,in: Number of kz mesh
  !!@param sw_pair,in: switch of singlet or triplet paring interacton
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  use grid_work, only:grid_work_get,tmpVdelta=>gw1,tmpfk=>gw2,tmp=>gw3  ! kept between calls
  use freq_grid
  implicit none
  integer(c_int64_t),intent(in):: Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  !! orbital permutation of the inversion operator (identity = old behaviour)
  integer(c_int64_t),intent(in),dimension(Norb):: iprm
  logical(1),intent(in):: sw_pair
  real(c_double),intent(in),dimension(Nchi,Nchi):: Smat,Cmat
  real(c_double),intent(in),dimension(Norb):: prt
  complex(c_double),intent(in),dimension(Nk,Nw,Nchi,Nchi):: Vdelta
  complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: delta
  complex(c_double),intent(out),dimension(Nk,Nw,Norb,Norb):: newdelta
  
  integer(c_int32_t) i,j,k,n,l,m
  integer(c_int64_t) jp,jm,j0
  real(c_double) sgn,dprt
  complex(c_double),dimension(Norb,Norb):: dtrim
  logical(1),parameter:: sw_odd_freq=.false.
  complex(c_double) vinf,static_sum

  ! sgn: Pauli, F(-k,-ω) = sgn·F^T(k,ω);  dprt: inversion parity, F(-k,ω) = dprt·P F(k,ω) P
  ! -ω is rebuilt as F(k,-ω) = sgn·dprt·P F^*(k,ω) P, using F^T=F^* (Δ Hermitian: TRS lets the
  ! gauge be chosen so, exactly, for singlet d0 and triplet dz).  No frequency parity is assumed,
  ! so the odd-ω (orbital-antisymmetric) admixture of multi-orbital gaps is kept.  Since only the
  ! lower triangle is computed and the rest filled by conjugation, the map is real-linear only:
  ! lin_eliash therefore runs Arnoldi over the real space of Hermitian Δ (real inner product).
  call grid_work_get(Nx,Ny,Nz,fg_nt(Nw))
  j0=fg_pos(FG_BOSON,1_c_int64_t)   ! v=0
  if(sw_pair)then
     sgn=1.0d0 !singlet Fk=F^T-k
     dprt=1.0d0 !gap parity even freq singlet is even
  else
     sgn=-1.0d0 !triplet Fk=-F^T-k
     dprt=-1.0d0 !gap parity even freq triplet is odd
  end if

  if(sw_odd_freq)then !consider odd frequency
     dprt=-dprt
  end if
  !$omp parallel workshare
  newdelta(:,:,:,:)=0.0d0
  !$omp end parallel workshare
  do l=1,Nchi
     do m=1,Nchi
        if(olist(l,1)>=olist(m,2))then
           !$omp parallel
           !$omp do private(i)
           do i=1,Nkall
              if(invk(2,i)==0)then !k,0
                 tmpVdelta(kmap(1,i),kmap(2,i),kmap(3,i),j0)=Vdelta(invk(1,i),1,m,l) !stored j=1 is v=0 (slot j0)
              else if(invk(2,i)==1)then !-k,0
                 tmpVdelta(kmap(1,i),kmap(2,i),kmap(3,i),j0)=Vdelta(invk(1,i),1,l,m) !V^ml(-q,0)=V^lm(q,0)
              end if
           end do
           !$omp end do
           !$omp do private(i,j,jp,jm)
           do j=2,Nw
              jp=fg_pos(FG_BOSON,j)
              jm=fg_neg(FG_BOSON,j,Nw)
              do i=1,Nkall
                 if(invk(2,i)==0)then !k,iw,k,-iw
                    tmpVdelta(kmap(1,i),kmap(2,i),kmap(3,i),jp)=Vdelta(invk(1,i),j,m,l)
                    tmpVdelta(kmap(1,i),kmap(2,i),kmap(3,i),jm)=conjg(Vdelta(invk(1,i),j,l,m)) !V^ml(q,-iw)=V^*lm(q,iw)
                 else if(invk(2,i)==1)then !-k,iw, -k,-iw
                    tmpVdelta(kmap(1,i),kmap(2,i),kmap(3,i),jp)=Vdelta(invk(1,i),j,l,m) !V^ml(-q,iw)=V^lm(q,iw)
                    tmpVdelta(kmap(1,i),kmap(2,i),kmap(3,i),jm)=conjg(Vdelta(invk(1,i),j,m,l)) !V^ml(-q,-iw)=V^*ml(q,iw)
                 end if
              end do
           end do
           !$omp end do
           !$omp do private(i,j,jp,jm)
           do j=1,Nw
              jp=fg_pos(FG_FERMION,j)
              jm=fg_neg(FG_FERMION,j,Nw)
              do i=1,Nkall
                 ! The single-inversion terms below carry P = signed permutation, i.e.
                 ! p_m p_l on the amplitude AND iprm on the orbital indices. The
                 ! F(-k,-w) term applies P twice (P^2=1), so it needs neither.
                 if(invk(2,i)==0)then !k
                    tmpfk(kmap(1,i),kmap(2,i),kmap(3,i),jp)=delta(invk(1,i),j,olist(l,2),olist(m,1))
                    tmpfk(kmap(1,i),kmap(2,i),kmap(3,i),jm)=sgn*dprt*prt(olist(l,2))*prt(olist(m,1))*&
                         conjg(delta(invk(1,i),j,iprm(olist(l,2)),iprm(olist(m,1)))) !F(k,-w)=sgn*dprt*p*p*F^*(k,w) (no soc only)
                 else if(invk(2,i)==1)then !-k
                    tmpfk(kmap(1,i),kmap(2,i),kmap(3,i),jp)=prt(olist(l,2))*prt(olist(m,1))*dprt*delta(invk(1,i),j,iprm(olist(l,2)),iprm(olist(m,1))) !F(-k,w)=p*p*dprt*F(k,w)
                    tmpfk(kmap(1,i),kmap(2,i),kmap(3,i),jm)=sgn*delta(invk(1,i),j,olist(m,1),olist(l,2)) !F(-k,-w)=sgn*F(k,w)
                 end if
              end do
           end do
           !$omp end do
           !$omp end parallel
           vinf=cmplx(0.5d0*(Smat(m,l)+sgn*Cmat(m,l)),0.0d0,kind=c_double) !bare vertex
           call fg_set_static(tmpVdelta,vinf,Nx,Ny,Nz,Nw)
           call fg_forward(FG_BOSON,tmpVdelta,tmp,Nx,Ny,Nz,Nw)
           static_sum=fg_static_term(tmpfk,Nx,Ny,Nz,Nw)   ! 0 on the uniform grid
           call fg_forward(FG_FERMION,tmpfk,tmp,Nx,Ny,Nz,Nw)
           !$omp parallel
           !$omp workshare
           tmp(:,:,:,:)=0.0d0
           !$omp end workshare
           !$omp do private(i,j,k)
           do n=1,int(fg_nt(Nw),c_int32_t)
              do k=0,Nz-1
                 do j=0,Ny-1
                    !$omp simd
                    do i=0,Nx-1
                       tmp(i,j,k,n)=tmpVdelta(i,j,k,n)*tmpfk(i,j,k,n)
                    end do
                    !$omp end simd
                 end do
              end do
           end do
           !$omp end do
           !$omp workshare
           tmpfk(:,:,:,:)=0.0d0
           !$omp end workshare
           !$omp end parallel
           call fg_backward(FG_FERMION,tmp,tmpfk,Nx,Ny,Nz,Nw)
           !$omp parallel do private(i,jp)
           do j=1,Nw
              jp=fg_pos(FG_FERMION,j)
              do i=1,Nkall
                 if(invk(2,i)==0)then
                    newdelta(invk(1,i),j,olist(l,1),olist(m,2))=newdelta(invk(1,i),j,olist(l,1),olist(m,2))&
                         +tmp(kmap(1,i),kmap(2,i),kmap(3,i),jp)
                 end if
              end do
           end do
           !$omp end parallel do
           if(.not.fg_uniform())then   ! IR: the bare vertex was kept out of the transform
              newdelta(:,:,olist(l,1),olist(m,2))=newdelta(:,:,olist(l,1),olist(m,2))+vinf*static_sum
           end if
        end if
     end do
  end do

  !$omp parallel
  do l=1,Norb
     do m=l,Norb
        !$omp do private(i,j)
        do j=1,Nw
           do i=1,Nk
              if(l/=m)then
                 newdelta(i,j,l,m)=conjg(newdelta(i,j,m,l))
              else
                 newdelta(i,j,l,l)=dble(newdelta(i,j,l,l))
              end if
           end do
        end do
        !$omp end do
     end do
  end do
  !$omp end parallel

  ! k=-k points (TRIM) are stored once and never rebuilt from parity, so impose it there:
  ! Δ(K) <- [Δ(K) + dprt·P Δ(K) P]/2.  Without this an odd-parity (triplet) gap keeps a spurious,
  ! parity-violating component at TRIM, which in multi-orbital models gives a wrong eigenvalue.
  ! The test is 2K = 0 modulo a reciprocal lattice vector, i.e. K and -K are the same mesh point.
  !$omp parallel do private(i,j,l,m,dtrim)
  do i=1,Nkall
     if(invk(2,i)==0 .and. modulo(2*kmap(1,i),Nx)==0 .and. modulo(2*kmap(2,i),Ny)==0&
          .and. modulo(2*kmap(3,i),Nz)==0)then
        do j=1,Nw
           do m=1,Norb
              do l=1,Norb
                 dtrim(l,m)=dprt*prt(l)*prt(m)*newdelta(invk(1,i),j,iprm(l),iprm(m))
              end do
           end do
           newdelta(invk(1,i),j,:,:)=0.5d0*(newdelta(invk(1,i),j,:,:)+dtrim(:,:))
        end do
     end if
  end do
  !$omp end parallel do
end subroutine mkdelta_nsoc

subroutine get_initial_delta(delta,init_delta,uni,kmap,invk,Nkall,Nk,Nw,Norb,gap_sym)
  !> make orbital basis initial gap function
  !!@param     delta,out: gap function
  !!@param init_delta,in: band basis initial gap function
  !!@param        uni,in: unitary matrix
  !!@param       kmap,in: property of k-points
  !!@param       invk,in: footnote of reverse k
  !!@param      Nkall,in: Number of all k-points
  !!@param         Nk,in: Number of irreducible k-points
  !!@param       Norb,in: Number of orbitals
  !!@param    gap_sym,in: gap symmetry number; its SIGN is the channel (>0 singlet, <0 triplet)
  !!                      and 0 means the unit-matrix seed
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  use constant
  implicit none
  integer(c_int64_t),intent(in):: Nw,Norb,Nkall,Nk,gap_sym
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
  real(c_double),intent(in),dimension(Nk,Norb):: init_delta
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(c_double),intent(out),dimension(Nk,Nw,Norb,Norb):: delta

  integer(c_int32_t) i,j,l,m,n
  real(c_double) norm
  complex(c_double),allocatable,dimension(:,:,:):: delta0
  complex(c_double),dimension(Norb,Norb):: tmpu

  allocate(delta0(Nk,Norb,Norb))
  if(gap_sym==0)then
     !$omp parallel
     !$omp workshare
     delta(:,:,:,:)=0.0d0
     !$omp end workshare
     do l=1,Norb
        !$omp workshare
        delta(:,:,l,l)=1.0d0
        !$omp end workshare
     end do
     !$omp end parallel
  else
     !$omp parallel
     !$omp workshare
     delta0(:,:,:)=0.0d0
     !$omp end workshare
     !$omp end parallel

     !$omp parallel do private(i,l,m,n,tmpu)
     do i=1,Nk
        do n=1,Norb
           tmpu(:,n)=uni(:,n,i)*init_delta(i,n)
        end do
        call zgemm('N','C',Norb,Norb,Norb,(1.0d0,0.0d0),tmpu,Norb,uni(:,:,i),Norb,(0.0d0,0.0d0),delta0(i,:,:),Norb)
        do l=1,Norb
           delta0(i,l,l)=dble(delta0(i,l,l))
           do m=l+1,Norb
              delta0(i,l,m)=conjg(delta0(i,m,l))
           end do
        end do
     end do
     !$omp end parallel do

     !$omp parallel do private(j)
     do j=1,Nw
        delta(:,j,:,:)=delta0(:,:,:)
     end do
     !$omp end parallel do
  end if

  norm=0.0d0
  !$omp parallel
  do l=1,Norb
     do m=1,Norb
        !$omp do private(j,i) reduction(+:norm)
        do j=1,Nw
           do i=1,Nkall
              norm=norm+abs(delta(invk(1,i),j,m,l))*abs(delta(invk(1,i),j,m,l))
           end do
        end do
        !$omp end do
     end do
  end do
  !$omp end parallel
  ! norm=0 check
  if(norm <= 0.0d0)then
     print*,'Error: get_initial_delta norm <= 0'
     stop
  end if
  !$omp parallel
  !$omp workshare
  delta(:,:,:,:)=delta(:,:,:,:)/sqrt(2.0d0*norm)
  !$omp end workshare
  !$omp end parallel
end subroutine get_initial_delta

subroutine conv_delta_orb_to_band(deltab,delta,uni,prt,iprm,invk,Norb,Nkall,Nk,Nw,gap_sym) bind(C)
  !> This function transform orbital basis gap function into band basis gap function
  !!@param deltab,out: band basis gap function
  !!@param   delta,in: orbital basis gap function
  !!@param     uni,in: unitary matrix
  !!@param    invk,in: footnote of reverse k
  !!@param    Norb,in: Number of orbital
  !!@param   Nkall,in: Number of all k-points
  !!@param      Nk,in: Number of irreducible k-points
  !!oparam      Nw,in: Number of Matsubara frequencies
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nw,Norb,Nk,Nkall,gap_sym
  integer(c_int64_t),intent(in),dimension(3,Nkall):: invk
  !! orbital permutation of the inversion operator (identity = old behaviour)
  integer(c_int64_t),intent(in),dimension(Norb):: iprm
  real(c_double),intent(in),dimension(Norb):: prt
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: uni
  complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: delta
  complex(c_double),intent(out),dimension(Nkall,Norb,Norb):: deltab

  integer(c_int32_t) i,k,l,n
  real(c_double) dprt
  complex(c_double),dimension(Norb,Norb):: matL,matR,tmpm,dtmp

  if(gap_sym<0)then
     dprt=-1.0d0 !odd
  else
     dprt=1.0d0 !even
  end if

  !$omp parallel do private(i,k,l,n,matL,matR,tmpm,dtmp)
  do i=1,Nkall
     k=invk(1,i)
     if(invk(2,i)==0)then
        matL=conjg(transpose(uni(:,:,k)))
        call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),matL,Norb,delta(k,1,:,:),Norb,(0.0d0,0.0d0),tmpm,Norb)
        call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),tmpm,Norb,uni(:,:,k),Norb,(0.0d0,0.0d0),deltab(i,:,:),Norb)
     else if(invk(2,i)==1)then
        ! Delta_orb(-k) = dprt * P Delta_orb(k) P^T with P the signed permutation
        ! (iprm,prt). Forming it explicitly is identical to the old column scaling
        ! when iprm is the identity, and stays correct when it is not.
        do l=1,Norb
           do n=1,Norb
              dtmp(n,l)=dprt*prt(n)*prt(l)*delta(k,1,iprm(n),iprm(l))
           end do
        end do
        matL=transpose(uni(:,:,k))
        call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),matL,Norb,dtmp,Norb,(0.0d0,0.0d0),tmpm,Norb)
        matR=conjg(uni(:,:,k))
        call zgemm('N','N',Norb,Norb,Norb,(1.0d0,0.0d0),tmpm,Norb,matR,Norb,(0.0d0,0.0d0),deltab(i,:,:),Norb)
     else
        deltab(i,:,:)=0.0d0
     end if
  end do
  !$omp end parallel do
end subroutine conv_delta_orb_to_band

subroutine remap_delta(delta,delta0,prt,iprm,invk,Nkall,Nk,Nw,Norb,gap_sym) bind(C)
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nkall,Nk,Nw,Norb,gap_sym
  integer(c_int64_t),intent(in),dimension(3,Nkall):: invk
  !! orbital permutation of the inversion operator (identity = old behaviour)
  integer(c_int64_t),intent(in),dimension(Norb):: iprm
  real(c_double),intent(in),dimension(Norb):: prt
  complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: delta0
  complex(c_double),intent(out),dimension(Nkall,Nw,Norb,Norb):: delta

  integer(c_int32_t) i,j,l,m
  real(c_double) dprt
  if(gap_sym<0)then
     dprt=-1.0d0
  else
     dprt=1.0d0
  end if
  
  !$omp parallel
  do l=1,Norb
     do m=1,Norb
        !$omp do private(i,j)
        do j=1,Nw
           do i=1,Nkall
              if(invk(2,i)==0)then
                 delta(i,j,m,l)=delta0(invk(1,i),j,m,l)
              else if(invk(2,i)==1)then
                 delta(i,j,m,l)=dprt*prt(m)*prt(l)*delta0(invk(1,i),j,iprm(m),iprm(l))
              end if
           end do
        end do
        !$omp end do
     end do
  end do
  !$omp end parallel
end subroutine remap_delta
