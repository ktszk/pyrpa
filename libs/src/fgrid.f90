! -----------------------------------------------------------------------------
! Grid kernels for the vortex / vortex-lattice quasiclassical solvers.
!
! These are the innermost real-space operations of libs/plibs/_eilenberger_vortex.py
! and _eilenberger_spin.py: the order parameter, the impurity self-energy, the
! vector potential and the propagators all live on small 2D grids but are sampled
! at O(10^4-10^6) scattered chord points, for every Fermi-velocity direction of
! every self-consistency iteration.  Done with scipy's RegularGridInterpolator
! (two calls per complex field, plus temporaries) that sampling dominated the
! solvers; here it is one cache-friendly pass with the frequency axis innermost.
!
! ctypes-callable (bind(C)); arrays follow the numpy C-order <-> Fortran
! reversed-dimension convention used throughout pyrpa: a numpy array [n,m] is
! declared here as dimension(m,n).
! -----------------------------------------------------------------------------

subroutine bilinear_eval(out,fld,x0,dx,xhi,ng,px,py,npt,nw,fill) bind(C)
!> bilinear_eval
!> Bilinear interpolation of a complex field on a UNIFORM square grid
!> x(i)=y(i)=x0+(i-1)*dx, evaluated at scattered points (px,py) and vectorized
!> over a trailing (frequency) axis.  Points outside the grid get `fill`, which
!> reproduces scipy RegularGridInterpolator(method='linear', bounds_error=False,
!> fill_value=fill).
!!@param  out,out: interpolated values [npt,nw] complex128
!!@param  fld, in: field on the grid [ng,ng,nw] complex128 (first axis = x)
!!@param   x0, in: first grid coordinate (both axes)
!!@param   dx, in: grid spacing (both axes, uniform)
!!@param  xhi, in: LAST grid coordinate, passed explicitly rather than rebuilt as
!!                 x0+(ng-1)*dx: on a linspace axis the two differ by an ulp, and
!!                 whole rows of chord points sit exactly on the boundary, which
!!                 would then fall outside and be filled instead of interpolated
!!@param   ng, in: grid points per axis
!!@param px,py,in: sample coordinates [npt]
!!@param  npt, in: number of sample points
!!@param   nw, in: length of the trailing axis (1 for a plain 2D field)
!!@param fill, in: value returned outside the grid
  use,intrinsic:: iso_c_binding, only: c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: ng,npt,nw
  real(c_double),intent(in):: x0,dx,xhi
  complex(c_double),intent(in):: fill
  real(c_double),intent(in),dimension(npt):: px,py
  complex(c_double),intent(in),dimension(nw,ng,ng):: fld
  complex(c_double),intent(out),dimension(nw,npt):: out
  integer(c_int64_t) p,w,i0,j0
  real(c_double) u,v,wu,wv,c00,c10,c01,c11
  !$omp parallel do default(none) &
  !$omp   shared(out,fld,x0,dx,xhi,ng,px,py,npt,nw,fill) &
  !$omp   private(p,w,i0,j0,u,v,wu,wv,c00,c10,c01,c11)
  do p=1,npt
     if(px(p)<x0 .or. px(p)>xhi .or. py(p)<x0 .or. py(p)>xhi)then
        do w=1,nw
           out(w,p)=fill
        end do
        cycle
     end if
     u=(px(p)-x0)/dx
     v=(py(p)-x0)/dx
     i0=int(floor(u),c_int64_t)+1                 ! 1-based lower cell corner
     j0=int(floor(v),c_int64_t)+1
     if(i0<1) i0=1
     if(j0<1) j0=1
     if(i0>ng-1) i0=ng-1                          ! upper edge falls in the last cell
     if(j0>ng-1) j0=ng-1
     wu=u-real(i0-1,c_double)
     wv=v-real(j0-1,c_double)
     c00=(1d0-wu)*(1d0-wv); c10=wu*(1d0-wv)
     c01=(1d0-wu)*wv;       c11=wu*wv
     do w=1,nw
        out(w,p)=c00*fld(w,j0,i0)+c10*fld(w,j0,i0+1) &
                +c01*fld(w,j0+1,i0)+c11*fld(w,j0+1,i0+1)
     end do
  end do
end subroutine bilinear_eval


subroutine bilinear_cell_eval(out,fld,ng,px,py,npt,Dlx,Dly,zeta) bind(C)
!> bilinear_cell_eval
!> Bilinear interpolation of a complex CELL-CENTRED field on the (oblique)
!> magnetic unit cell, with periodic wrap: grid point k sits at fractional
!> coordinate (k+0.5)/ng-0.5, and Cartesian points are mapped to fractional ones
!> by the inverse of x=(r1(1-zeta)+r2 zeta)Dlx, y=(-r1+r2)Dly.  Mirrors
!> _periodic_eval / _to_frac in _eilenberger_vortex.py.
!!@param out,out: interpolated values [npt] complex128
!!@param fld, in: cell field [ng,ng] complex128
!!@param  ng, in: grid points per axis
!!@param px,py,in: Cartesian sample coordinates [npt]
!!@param npt, in: number of sample points
!!@param Dlx,Dly,zeta,in: cell geometry (widths and obliqueness)
  use,intrinsic:: iso_c_binding, only: c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: ng,npt
  real(c_double),intent(in):: Dlx,Dly,zeta
  real(c_double),intent(in),dimension(npt):: px,py
  complex(c_double),intent(in),dimension(ng,ng):: fld
  complex(c_double),intent(out),dimension(npt):: out
  integer(c_int64_t) p,i0,j0,i1,j1
  real(c_double) f1,f2,u,v,wu,wv,gn
  gn=real(ng,c_double)
  !$omp parallel do default(none) &
  !$omp   shared(out,fld,ng,gn,px,py,npt,Dlx,Dly,zeta) &
  !$omp   private(p,i0,j0,i1,j1,f1,f2,u,v,wu,wv)
  do p=1,npt
     f1=px(p)/Dlx-zeta*py(p)/Dly
     f2=px(p)/Dlx+(1d0-zeta)*py(p)/Dly
     u=modulo(f1+0.5d0,1d0)*gn-0.5d0
     v=modulo(f2+0.5d0,1d0)*gn-0.5d0
     i0=int(floor(u),c_int64_t)
     j0=int(floor(v),c_int64_t)
     wu=u-real(i0,c_double)
     wv=v-real(j0,c_double)
     i1=modulo(i0+1,ng); j1=modulo(j0+1,ng)
     i0=modulo(i0,ng);   j0=modulo(j0,ng)
     out(p)=(1d0-wu)*(1d0-wv)*fld(j0+1,i0+1)+wu*(1d0-wv)*fld(j0+1,i1+1) &
           +(1d0-wu)*wv*fld(j1+1,i0+1)      +wu*wv*fld(j1+1,i1+1)
  end do
end subroutine bilinear_cell_eval


subroutine abrikosov_z(out,px,py,npt,Dlx,Dly,zeta,nsum) bind(C)
!> abrikosov_z
!> Abrikosov (lowest-Landau-level) quasi-periodic order parameter of the magnetic
!> cell, evaluated at scattered points:
!>   Z = (2 Dr)^(1/4) e^{-i pi xl yl}
!>       sum_p e^{-pi Dr (yl+Dy0+p)^2} e^{-2i pi (p(Dx0+p zeta/2)+(Dy0+p) xl)}
!> with xl=x/Dlx, yl=y/Dly, Dr=Dly/Dlx, Dx0=-(1+zeta)/2, Dy0=-1/2.  Mirrors
!> _abrikosov_z; it seeds the lattice gap and, through its unit phase, supplies
!> the winding and the London vector potential, so it is evaluated on every
!> chord of every direction.
!!@param out,out: Z at the sample points [npt] complex128
!!@param px,py,in: Cartesian sample coordinates [npt]
!!@param npt, in: number of sample points
!!@param Dlx,Dly,zeta,in: cell geometry (widths and obliqueness)
!!@param nsum, in: LLL sum truncation (p = -nsum..nsum)
  use,intrinsic:: iso_c_binding, only: c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: npt,nsum
  real(c_double),intent(in):: Dlx,Dly,zeta
  real(c_double),intent(in),dimension(npt):: px,py
  complex(c_double),intent(out),dimension(npt):: out
  integer(c_int64_t) p,k
  real(c_double) Dr,Dx0,Dy0,xl,yl,pref,ph,pi
  complex(c_double) z
  pi=4d0*atan(1d0)
  Dr=Dly/Dlx
  Dx0=-0.5d0*(1d0+zeta)
  Dy0=-0.5d0
  pref=sqrt(sqrt(2d0*Dr))
  !$omp parallel do default(none) &
  !$omp   shared(out,px,py,npt,Dlx,Dly,zeta,nsum,Dr,Dx0,Dy0,pref,pi) &
  !$omp   private(p,k,xl,yl,ph,z)
  do p=1,npt
     xl=px(p)/Dlx
     yl=py(p)/Dly
     z=(0d0,0d0)
     do k=-nsum,nsum
        ph=-2d0*pi*(real(k,c_double)*(Dx0+real(k,c_double)*zeta*0.5d0) &
                    +(Dy0+real(k,c_double))*xl)
        z=z+exp(-pi*Dr*(yl+Dy0+real(k,c_double))**2)*cmplx(cos(ph),sin(ph),c_double)
     end do
     ph=-pi*xl*yl
     out(p)=pref*z*cmplx(cos(ph),sin(ph),c_double)
  end do
end subroutine abrikosov_z
