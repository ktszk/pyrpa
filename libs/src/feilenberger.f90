! -----------------------------------------------------------------------------
! Quasiclassical Eilenberger / Riccati hot loops (scalar + 2x2 spin matrix).
!
! These are the innermost trajectory integrations shared by the surface, vortex,
! lattice and d-vector solvers in libs/plibs/_eilenberger*.py.  Each is an
! unconditionally stable analytic step (scalar tanh / matrix Mobius) along a
! quasiclassical chord, vectorized/parallelized over the independent chords and
! Matsubara (or real) frequencies.  ctypes-callable (bind(C)); arrays follow the
! numpy C-order <-> Fortran reversed-dimension convention used throughout pyrpa:
! a numpy array [Ns,Nchord,Nw] is declared here as dimension(Nw,Nchord,Ns).
! -----------------------------------------------------------------------------

subroutine riccati_chords(g,f,om,dd,hvf,ds,Ns,Nchord,Nw) bind(C)
!> riccati_chords
!> Scalar quasiclassical g, f along many straight chords, via the unconditionally
!> stable tanh step.  gamma is integrated forward from the upstream bulk root,
!> gamma-tilde backward from the downstream bulk root, then
!>   g = (1 - gamma*gammat)/(1 + gamma*gammat),  f = 2 gamma/(1 + gamma*gammat).
!> Mirrors _chord_gf / _chord_gf_pos / _integrate_vec (midpoint piecewise-constant
!> step) exactly; omega may be position dependent (Doppler / impurity).
!!@param   g,out: normal propagator [Ns,Nchord,Nw] complex128
!!@param   f,out: anomalous propagator [Ns,Nchord,Nw] complex128
!!@param   om, in: (renormalized) frequency along each chord [Ns,Nchord,Nw] complex128
!!@param   dd, in: order parameter along each chord [Ns,Nchord,Nw] complex128
!!@param  hvf, in: hbar |v_F|
!!@param   ds, in: arc-length step per chord [Nchord] (e.g. dx/|cos beta|)
!!@param   Ns, in: points per chord
!!@param Nchord,in: number of chords
!!@param   Nw, in: number of frequencies
  use,intrinsic:: iso_c_binding, only: c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Ns,Nchord,Nw
  real(c_double),intent(in):: hvf
  real(c_double),intent(in),dimension(Nchord):: ds
  complex(c_double),intent(in),dimension(Nw,Nchord,Ns):: om,dd
  complex(c_double),intent(out),dimension(Nw,Nchord,Ns):: g,f
  integer(c_int64_t) c,i,w
  real(c_double) t
  complex(c_double) D1,om1,R,Tt,gm,gt,Dmid,ommid,den
  complex(c_double),allocatable:: gamf(:,:),gamt(:,:)
  !$omp parallel do default(none) shared(g,f,om,dd,hvf,ds,Ns,Nchord,Nw) &
  !$omp   private(i,w,t,D1,om1,R,Tt,gm,gt,Dmid,ommid,den,gamf,gamt)
  do c=1,Nchord
     t=ds(c)/hvf
     allocate(gamf(Nw,Ns),gamt(Nw,Ns))
     ! forward gamma: upstream bulk root at i=1
     do w=1,Nw
        D1=dd(w,c,1); om1=om(w,c,1)
        R=sqrt(om1*om1+D1*conjg(D1))
        if(abs(D1)>0d0)then
           gamf(w,1)=(R-om1)/conjg(D1)
        else
           gamf(w,1)=(0d0,0d0)
        end if
     end do
     do i=1,Ns-1
        do w=1,Nw
           Dmid=0.5d0*(dd(w,c,i)+dd(w,c,i+1))
           ommid=0.5d0*(om(w,c,i)+om(w,c,i+1))
           R=sqrt(ommid*ommid+Dmid*conjg(Dmid))
           Tt=tanh(R*t)/R
           gm=gamf(w,i)
           gamf(w,i+1)=(gm+Tt*(Dmid-ommid*gm))/(1d0+Tt*(ommid+conjg(Dmid)*gm))
        end do
     end do
     ! backward gamma-tilde: downstream bulk root at i=Ns
     do w=1,Nw
        D1=dd(w,c,Ns); om1=om(w,c,Ns)
        R=sqrt(om1*om1+D1*conjg(D1))
        if(abs(D1)>0d0)then
           gamt(w,Ns)=(R-om1)/D1
        else
           gamt(w,Ns)=(0d0,0d0)
        end if
     end do
     do i=Ns,2,-1
        do w=1,Nw
           ! integrate with gap conj(dd) along the reversed path (== _integrate_vec)
           Dmid=0.5d0*(conjg(dd(w,c,i))+conjg(dd(w,c,i-1)))
           ommid=0.5d0*(om(w,c,i)+om(w,c,i-1))
           R=sqrt(ommid*ommid+Dmid*conjg(Dmid))
           Tt=tanh(R*t)/R
           gt=gamt(w,i)
           gamt(w,i-1)=(gt+Tt*(Dmid-ommid*gt))/(1d0+Tt*(ommid+conjg(Dmid)*gt))
        end do
     end do
     ! combine
     do i=1,Ns
        do w=1,Nw
           gm=gamf(w,i); gt=gamt(w,i)
           den=1d0+gm*gt
           g(w,c,i)=(1d0-gm*gt)/den
           f(w,c,i)=2d0*gm/den
        end do
     end do
     deallocate(gamf,gamt)
  end do
end subroutine riccati_chords


subroutine expm4(N,M)
!> 4x4 complex matrix exponential via eigendecomposition (zgeev): exp(N) =
!> V diag(exp(lambda)) V^{-1}.  Adequate for the (diagonalizable) Riccati Mobius
!> generators here; matches scipy.linalg.expm to working precision.
  use,intrinsic:: iso_c_binding, only: c_double
  implicit none
  complex(c_double),intent(in):: N(4,4)
  complex(c_double),intent(out):: M(4,4)
  complex(c_double) A(4,4),VR(4,4),VL(1,1),w(4),Vi(4,4),tmp(4,4)
  complex(c_double) work(64)
  real(c_double) rwork(8)
  integer ipiv(4),info,k
  A=N
  call zgeev('N','V',4,A,4,w,VL,1,VR,4,work,64,rwork,info)
  ! Vi = inv(VR)
  Vi=VR
  call zgetrf(4,4,Vi,4,ipiv,info)
  call zgetri(4,Vi,4,ipiv,work,64,info)
  ! tmp = diag(exp(w)) * Vi
  do k=1,4
     tmp(k,:)=exp(w(k))*Vi(k,:)
  end do
  M=matmul(VR,tmp)
end subroutine expm4


subroutine matrix_riccati_batch(g,f,om,Dpath,hvf,ds,h,Ns,Nw) bind(C)
!> matrix_riccati_batch
!> 2x2 spin-matrix quasiclassical g, f along ONE inhomogeneous trajectory, batched
!> over frequencies.  a is integrated forward from the upstream bulk root, b backward
!> (reversed path) from the downstream bulk root, each via the exact fractional-linear
!> (Mobius) step M=exp(N ds/hvf); then g=(I+ab)^{-1}(I-ab), f=(I+ab)^{-1} 2a.  Mirrors
!> matrix_trajectory_gf exactly.  The gap matrix path Dpath is frequency-independent;
!> OpenMP parallelizes over the frequencies.
!!@param   g,out: normal 2x2 propagator [Ns,Nw,2,2] complex128
!!@param   f,out: anomalous 2x2 propagator [Ns,Nw,2,2] complex128
!!@param  om, in: (renormalized) frequencies [Nw] complex128
!!@param Dpath,in: 2x2 gap matrix along the path [Ns,2,2] complex128
!!@param hvf,ds,h: hbar|v_F|, arc-length step, Zeeman energy
!!@param  Ns,Nw, in: points per path, number of frequencies
  use,intrinsic:: iso_c_binding, only: c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Ns,Nw
  real(c_double),intent(in):: hvf,ds,h
  complex(c_double),intent(in),dimension(Nw):: om
  complex(c_double),intent(in),dimension(2,2,Ns):: Dpath
  complex(c_double),intent(out),dimension(2,2,Nw,Ns):: g,f
  integer(c_int64_t) iw,i
  real(c_double) t
  complex(c_double) Dmat(2,2,Ns),Dmid(2,2),Nmat(4,4),Mexp(4,4)
  complex(c_double) a(2,2,Ns),b(2,2,Ns),aa(2,2),bb(2,2),ab(2,2),P(2,2),res(2,2)
  complex(c_double) sz(2,2),I2(2,2)
  t=ds/hvf
  sz=reshape([(1d0,0d0),(0d0,0d0),(0d0,0d0),(-1d0,0d0)],[2,2])
  I2=reshape([(1d0,0d0),(0d0,0d0),(0d0,0d0),(1d0,0d0)],[2,2])
  ! natural (row,col) gap matrices: numpy[i,p,q] = Dpath(q,p,i)
  do i=1,Ns
     Dmat(:,:,i)=transpose(Dpath(:,:,i))
  end do
  !$omp parallel do default(none) shared(g,f,om,Dmat,t,sz,I2,Ns,Nw,h) &
  !$omp   private(i,Dmid,Nmat,Mexp,a,b,aa,bb,ab,P,res)
  do iw=1,Nw
     ! forward a from upstream bulk root
     call bulk_root(om(iw),Dmat(:,:,1),.true.,a(:,:,1))
     do i=1,Ns-1
        Dmid=0.5d0*(Dmat(:,:,i)+Dmat(:,:,i+1))
        call riccati_gen(om(iw),Dmid,h,1d0,sz,Nmat)
        call expm4(Nmat*t,Mexp)
        call mobius(a(:,:,i),Mexp,a(:,:,i+1))
     end do
     ! backward b from downstream bulk root
     call bulk_root(om(iw),Dmat(:,:,Ns),.false.,b(:,:,Ns))
     do i=Ns,2,-1
        Dmid=0.5d0*(Dmat(:,:,i)+Dmat(:,:,i-1))
        call riccati_gen(om(iw),Dmid,h,-1d0,sz,Nmat)
        call expm4(Nmat*t,Mexp)
        call mobius(b(:,:,i),Mexp,b(:,:,i-1))
     end do
     ! combine: g=(I+ab)^{-1}(I-ab), f=(I+ab)^{-1} 2a
     do i=1,Ns
        aa=a(:,:,i); bb=b(:,:,i)
        ab=matmul(aa,bb)
        call inv2(I2+ab,P)
        res=matmul(P,I2-ab);  g(:,:,iw,i)=transpose(res)
        res=matmul(P,2d0*aa); f(:,:,iw,i)=transpose(res)
     end do
  end do
end subroutine matrix_riccati_batch


subroutine riccati_gen(om,D,h,sgn,sz,N)
!> 4x4 Mobius generator (matches _riccati_N): N=[[-2 om I + i sgn h sz, D],
!> [conj(D)^T, i sgn h sz]].
  use,intrinsic:: iso_c_binding, only: c_double
  implicit none
  complex(c_double),intent(in):: om,D(2,2),sz(2,2)
  real(c_double),intent(in):: h,sgn
  complex(c_double),intent(out):: N(4,4)
  complex(c_double) zee(2,2),I2(2,2)
  I2=reshape([(1d0,0d0),(0d0,0d0),(0d0,0d0),(1d0,0d0)],[2,2])
  zee=(0d0,1d0)*sgn*h*sz
  N(1:2,1:2)=-2d0*om*I2+zee
  N(1:2,3:4)=D
  N(3:4,1:2)=transpose(conjg(D))
  N(3:4,3:4)=zee
end subroutine riccati_gen


subroutine bulk_root(om,D,is_a,a)
!> Homogeneous (unitary) matrix Riccati root (matches riccati_matrix_bulk):
!> a = src/(om+E), E=sqrt(om^2+0.5 tr(D D^dag)), src = D (is_a) else conj(D)^T.
  use,intrinsic:: iso_c_binding, only: c_double
  implicit none
  complex(c_double),intent(in):: om,D(2,2)
  logical,intent(in):: is_a
  complex(c_double),intent(out):: a(2,2)
  complex(c_double) DDd(2,2),E,d2
  DDd=matmul(D,transpose(conjg(D)))
  d2=0.5d0*(DDd(1,1)+DDd(2,2))
  E=sqrt(om*om+d2)
  if(is_a)then
     a=D/(om+E)
  else
     a=transpose(conjg(D))/(om+E)
  end if
end subroutine bulk_root


subroutine mobius(a,M,anew)
!> One fractional-linear step: anew=(M11 a + M12)(M21 a + M22)^{-1}, M=exp(N t).
  use,intrinsic:: iso_c_binding, only: c_double
  implicit none
  complex(c_double),intent(in):: a(2,2),M(4,4)
  complex(c_double),intent(out):: anew(2,2)
  complex(c_double) u(2,2),w(2,2),wi(2,2)
  u=matmul(M(1:2,1:2),a)+M(1:2,3:4)
  w=matmul(M(3:4,1:2),a)+M(3:4,3:4)
  call inv2(w,wi)
  anew=matmul(u,wi)
end subroutine mobius


subroutine inv2(A,Ai)
!> 2x2 complex matrix inverse.
  use,intrinsic:: iso_c_binding, only: c_double
  implicit none
  complex(c_double),intent(in):: A(2,2)
  complex(c_double),intent(out):: Ai(2,2)
  complex(c_double) det
  det=A(1,1)*A(2,2)-A(1,2)*A(2,1)
  Ai(1,1)= A(2,2)/det
  Ai(2,2)= A(1,1)/det
  Ai(1,2)=-A(1,2)/det
  Ai(2,1)=-A(2,1)/det
end subroutine inv2


subroutine matrix_riccati_chords(g,f,om,Dpath,hvf,ds,h,Ns,Nchord,Nw) bind(C)
!> matrix_riccati_chords
!> Like matrix_riccati_batch but for MANY trajectories (chords): the 2x2 spin-matrix
!> Riccati g, f along each chord, batched over frequencies, OpenMP over chords.  Used
!> by the 2D vortex d-vector solver (one call per Fermi-surface direction).
!!@param   g,out: normal 2x2 propagator [Ns,Nchord,Nw,2,2] complex128
!!@param   f,out: anomalous 2x2 propagator [Ns,Nchord,Nw,2,2] complex128
!!@param  om, in: (renormalized) frequency along each chord [Ns,Nchord,Nw] complex128
!!                (position dependent, e.g. om + i v_F.Q Doppler shift)
!!@param Dpath,in: 2x2 gap matrix along each chord [Ns,Nchord,2,2] complex128
!!@param hvf,ds,h: hbar|v_F|, arc-length step, Zeeman energy
!!@param Ns,Nchord,Nw,in: points per chord, chords, frequencies
  use,intrinsic:: iso_c_binding, only: c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Ns,Nchord,Nw
  real(c_double),intent(in):: hvf,ds,h
  complex(c_double),intent(in),dimension(Nw,Nchord,Ns):: om
  complex(c_double),intent(in),dimension(2,2,Nchord,Ns):: Dpath
  complex(c_double),intent(out),dimension(2,2,Nw,Nchord,Ns):: g,f
  integer(c_int64_t) c,iw,i
  real(c_double) t
  complex(c_double) Dmid(2,2),Nmat(4,4),Mexp(4,4)
  complex(c_double) a(2,2,Ns),b(2,2,Ns),aa(2,2),bb(2,2),ab(2,2),P(2,2),res(2,2)
  complex(c_double) sz(2,2),I2(2,2)
  complex(c_double),allocatable:: Dmat(:,:,:,:)         ! (2,2,Ns,Nchord), natural (row,col)
  t=ds/hvf
  sz=reshape([(1d0,0d0),(0d0,0d0),(0d0,0d0),(-1d0,0d0)],[2,2])
  I2=reshape([(1d0,0d0),(0d0,0d0),(0d0,0d0),(1d0,0d0)],[2,2])
  allocate(Dmat(2,2,Ns,Nchord))                          ! precompute so (chord,freq) collapse
  do c=1,Nchord
     do i=1,Ns
        Dmat(:,:,i,c)=transpose(Dpath(:,:,c,i))
     end do
  end do
  ! parallelize over the whole (chord, frequency) space for full load balance
  !$omp parallel do collapse(2) default(none) shared(g,f,om,Dmat,t,sz,I2,Ns,Nchord,Nw,h) &
  !$omp   private(i,Dmid,Nmat,Mexp,a,b,aa,bb,ab,P,res)
  do c=1,Nchord
     do iw=1,Nw
        call bulk_root(om(iw,c,1),Dmat(:,:,1,c),.true.,a(:,:,1))
        do i=1,Ns-1
           Dmid=0.5d0*(Dmat(:,:,i,c)+Dmat(:,:,i+1,c))
           call riccati_gen(0.5d0*(om(iw,c,i)+om(iw,c,i+1)),Dmid,h,1d0,sz,Nmat)
           call expm4(Nmat*t,Mexp)
           call mobius(a(:,:,i),Mexp,a(:,:,i+1))
        end do
        call bulk_root(om(iw,c,Ns),Dmat(:,:,Ns,c),.false.,b(:,:,Ns))
        do i=Ns,2,-1
           Dmid=0.5d0*(Dmat(:,:,i,c)+Dmat(:,:,i-1,c))
           call riccati_gen(0.5d0*(om(iw,c,i)+om(iw,c,i-1)),Dmid,h,-1d0,sz,Nmat)
           call expm4(Nmat*t,Mexp)
           call mobius(b(:,:,i),Mexp,b(:,:,i-1))
        end do
        do i=1,Ns
           aa=a(:,:,i); bb=b(:,:,i)
           ab=matmul(aa,bb)
           call inv2(I2+ab,P)
           res=matmul(P,I2-ab);  g(:,:,iw,c,i)=transpose(res)
           res=matmul(P,2d0*aa); f(:,:,iw,c,i)=transpose(res)
        end do
     end do
  end do
  deallocate(Dmat)
end subroutine matrix_riccati_chords
