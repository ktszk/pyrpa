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
!!@param   ds, in: arc-length step
!!@param   Ns, in: points per chord
!!@param Nchord,in: number of chords
!!@param   Nw, in: number of frequencies
  use,intrinsic:: iso_c_binding, only: c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Ns,Nchord,Nw
  real(c_double),intent(in):: hvf,ds
  complex(c_double),intent(in),dimension(Nw,Nchord,Ns):: om,dd
  complex(c_double),intent(out),dimension(Nw,Nchord,Ns):: g,f
  integer(c_int64_t) c,i,w
  real(c_double) t
  complex(c_double) D1,om1,R,Tt,gm,gt,Dmid,ommid,den
  complex(c_double),allocatable:: gamf(:,:),gamt(:,:)
  t=ds/hvf
  !$omp parallel do private(i,w,D1,om1,R,Tt,gm,gt,Dmid,ommid,den,gamf,gamt)
  do c=1,Nchord
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
