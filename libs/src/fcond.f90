subroutine calc_lij(L11,L22,L12,vk,eig,ffermi,Norb,Nk,mu,w,idelta,eps,temp) bind(C)
  !> calc_Lij
  !!@param    L11,out: L11@linear response theory
  !!@param    L22,out: L22@linear response theory
  !!@param    L12,out: L12(L21^T)@linear response theory
  !!@param      vk,in: group velocity:
  !!@param     eig,in: energy of bands
  !!@param  ffermi,in: fermi distribute function
  !!@param    Norb,in: The number of orbitals
  !!@param      Nk,in: The number of k-points
  !!@param      mu,in: chemical potential
  !!@param       w,in: frequency (energy)
  !!@param  idelta,in: dumping factor
  !!@param     eps,in: threshold of energy
  !!@param    temp,in: Temperature
  use occupation, only:occ_factor
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,Norb
  real(c_double),intent(in):: temp,eps,w,idelta,mu
  real(c_double),intent(in),dimension(Norb,Nk):: eig,ffermi
  complex(c_double),intent(in),dimension(3,Norb,Norb,Nk):: vk
  complex(c_double),intent(out),dimension(3,3):: L11,L12,L22

  integer(c_int64_t) i,j,k,l,m
  complex(c_double) tmp,denom
  real(c_double) ebar,pocc,temp_safe
  complex(c_double),parameter::ii=(0.0d0,1.0d0)

  temp_safe=max(temp,1.0d-12)
  !$omp parallel
  !$omp workshare
  L11(:,:)=0.0d0
  L12(:,:)=0.0d0
  L22(:,:)=0.0d0
  !$omp end workshare
  !$omp do reduction(+: L11,L12,L22) private(i,l,m,j,k,tmp,ebar,pocc,denom)
  ! Kubo formula, ONE expression for every band pair (l,m):
  !   L_ij = (i/Nk) sum_{k,l,m} v_ml v_lm * p(l,m) / (w + e_m - e_l + i*delta)
  ! with p(l,m) = (f_l-f_m)/(e_m-e_l) -> -df/de at e_m = e_l (occ_factor).
  ! The l=m term is then exactly the Drude weight f(1-f)/(T(w+i*delta)), so the
  ! intra- and inter-band branches are the same function and the crossover is
  ! smooth.  The old split (hard 1e-9 degeneracy test plus a |f_l-f_m| > eps
  ! skip) silently dropped every pair whose splitting was below ~eps*4T, i.e.
  ! exactly the weakly split pairs that carry a full intra-band weight.
  ! The thermoelectric/thermal moments use the symmetrized heat-current vertex
  ! ebar = (e_l+e_m)/2 - mu (reduces to (e-mu) intraband; keeps L12 = L21^T / Onsager).
  k_loop: do i=1,Nk
     band_loop1: do l=1,Norb
        band_loop2: do m=1,Norb
           pocc=occ_factor(eig(l,i),eig(m,i),ffermi(l,i),ffermi(m,i),temp_safe)
           ebar=0.5d0*(eig(l,i)+eig(m,i))-mu
           ! Negligibility is judged on the LARGEST moment (L22 carries an extra
           ! ebar^2, so a pair far from mu can have a tiny pocc and still matter
           ! there) and NEVER for l==m: that is the Drude weight, which the
           ! Boltzmann kernel calc_kn sums without any threshold, so skipping it
           ! would break the exact dc Kubo <-> Boltzmann correspondence.
           ! (The old |f_l-f_m| > eps test did both wrong: it biased L12 by 13%
           ! and L22 by 5% on a two-orbital model.)
           if(l/=m .and. abs(pocc)*max(1.0d0,ebar*ebar)<=eps)cycle band_loop2
           denom=cmplx(w+eig(m,i)-eig(l,i),idelta,kind=c_double)
           do j=1,3
              do k=1,3
                 tmp=vk(k,m,l,i)*vk(j,l,m,i)*pocc/denom
                 L11(k,j)=L11(k,j)+tmp
                 L12(k,j)=L12(k,j)+tmp*ebar
                 L22(k,j)=L22(k,j)+tmp*ebar*ebar
              end do
           end do
        end do band_loop2
     end do band_loop1
  end do k_loop
  !$omp end do
  !$omp workshare
  L11(:,:)=ii*L11(:,:)/Nk
  L12(:,:)=ii*L12(:,:)/Nk
  L22(:,:)=ii*L22(:,:)/Nk
  !$omp end workshare
  !$omp end parallel
end subroutine calc_lij

subroutine calc_lij_wl(L11,L22,L12,vk,eig,ffermi,Norb,Nk,Nw,mu,wl,idelta,eps,temp) bind(C)
  !> calc_Lij for all Nw frequencies in a single k-loop pass.
  !> Precomputes frequency-independent prefactors per (k,l,m) pair,
  !> then applies the w-dependent denominator in the inner loop.
  !!@param  L11,out: L11 [Nw,3,3]
  !!@param  L22,out: L22 [Nw,3,3]
  !!@param  L12,out: L12 [Nw,3,3]
  !!@param   vk,in: group velocity [3,Norb,Norb,Nk]
  !!@param  eig,in: energy of bands [Norb,Nk]
  !!@param   ff,in: fermi distribute function [Norb,Nk]
  !!@param Norb,in: The number of orbitals
  !!@param   Nk,in: The number of k-points
  !!@param   Nw,in: The number of frequency points
  !!@param   mu,in: chemical potential
  !!@param   wl,in: frequency mesh [Nw]
  !!@param idelta,in: dumping factor
  !!@param  eps,in: threshold of energy
  !!@param temp,in: Temperature
  use occupation, only:occ_factor
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,Norb,Nw
  real(c_double),intent(in):: temp,eps,idelta,mu
  real(c_double),intent(in),dimension(Nw):: wl
  real(c_double),intent(in),dimension(Norb,Nk):: eig,ffermi
  complex(c_double),intent(in),dimension(3,Norb,Norb,Nk):: vk
  complex(c_double),intent(out),dimension(Nw,3,3):: L11,L12,L22

  integer(c_int64_t) i,j,k,l,m,iw
  real(c_double) de,fl,fm,el,em,ebar,pocc,temp_safe
  complex(c_double) p11,p12,p22,denom
  complex(c_double),parameter::ii=(0.0d0,1.0d0)

  L11(:,:,:)=(0.0d0,0.0d0)
  L12(:,:,:)=(0.0d0,0.0d0)
  L22(:,:,:)=(0.0d0,0.0d0)
  temp_safe=max(temp,1.0d-12)

  ! Same single smooth kernel as calc_lij (see the comment there); the
  ! frequency-independent prefactor is hoisted out of the w-loop.
  !$omp parallel do reduction(+: L11,L12,L22) private(i,l,m,j,k,iw,de,fl,fm,el,em,ebar,pocc,p11,p12,p22,denom)
  k_loop: do i=1,Nk
     band_loop1: do l=1,Norb
        fl=ffermi(l,i); el=eig(l,i)
        band_loop2: do m=1,Norb
           fm=ffermi(m,i); em=eig(m,i)
           de=em-el
           pocc=occ_factor(el,em,fl,fm,temp_safe)
           ! symmetrized heat-current vertex ebar = (e_l+e_m)/2 - mu
           ebar=0.5d0*(el+em)-mu
           if(l/=m .and. abs(pocc)*max(1.0d0,ebar*ebar)<=eps)cycle band_loop2   ! see calc_lij
           do j=1,3
              do k=1,3
                 p11=vk(k,m,l,i)*vk(j,l,m,i)*pocc
                 p12=p11*ebar
                 p22=p12*ebar
                 do iw=1,Nw
                    denom=cmplx(wl(iw)+de,idelta,kind=c_double)
                    L11(iw,j,k)=L11(iw,j,k)+p11/denom
                    L12(iw,j,k)=L12(iw,j,k)+p12/denom
                    L22(iw,j,k)=L22(iw,j,k)+p22/denom
                 end do
              end do
           end do
        end do band_loop2
     end do band_loop1
  end do k_loop
  !$omp end parallel do

  L11(:,:,:)=ii*L11(:,:,:)/Nk
  L12(:,:,:)=ii*L12(:,:,:)/Nk
  L22(:,:,:)=ii*L22(:,:,:)/Nk
end subroutine calc_lij_wl

subroutine calc_kn(K0,K1,K2,eig,veloc,kweight,tau,temp,mu,Nk,Norb) bind(C)
  !> calc_Kn
  !> Kn_ij=sum_k(v_ki*v_kj*(e_k-mu)^n*(-df(e_k)/de))
  !!@param     K0,out: corresponds to sigma@Bolzmann theory
  !!@param     K1,out: corresponds to sigmaS@Bolzmann theory
  !!@param     K2,out: corresponds to kappa@Bolzmann theory
  !!@param     eig,in: energy of bands
  !!@param   veloc,in: group velocity:
  !!@param kweight,in: weight of k-points
  !!@param     tau,in: relaxation time
  !!@param    temp,in: Temperature
  !!@param      mu,in: chemical potential
  !!@param      Nk,in: The number of k-points
  !!@param    Norb,in: The number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nk,Norb
  real(c_double),intent(in):: temp,mu
  real(c_double),intent(in),dimension(Norb,Nk):: eig,tau
  real(c_double),intent(in),dimension(Nk):: kweight
  real(c_double),intent(in),dimension(3,Norb,Nk):: veloc
  real(c_double),intent(out),dimension(3,3):: K0,K1,K2

  real(c_double),dimension(Norb,Nk):: dfermi
  integer(c_int32_t) i,j,l,m
  real(c_double) tmp,temp_safe

  temp_safe=max(temp,1.0d-12)

  ! -df/de = 0.25*(1-tanh^2((e-mu)/2T))/T  [numerically stable derivative of Fermi-Dirac]
  !$omp parallel
  !$omp do private(j)
  get_dfermi: do i=1,Nk
     do j=1,Norb
        dfermi(j,i)=0.25d0*(1.0d0-tanh(0.5d0*(eig(j,i)-mu)/temp_safe)**2)/temp_safe
     end do
  end do get_dfermi
  !$omp end do
  !$omp workshare
  K0(:,:)=0.0d0
  K1(:,:)=0.0d0
  K2(:,:)=0.0d0
  !$omp end workshare
  
  !$omp do private(j,l,m,tmp) reduction(+:K0,K1,K2)
  ! K_n[ij] = sum_{k,band} v_i*v_j * tau * (-df/de) * (e-mu)^n * weight
  ! K0 ~ sigma (charge conductivity kernel)
  ! K1 ~ sigmaS (thermoelectric kernel)
  ! K2 ~ kappa (thermal conductivity kernel)
  get_Kn: do i=1,Nk
     band_loop: do j=1,Norb
        axis1: do l=1,3
           !$omp simd
           axis2: do m=1,3
              tmp=veloc(m,j,i)*veloc(l,j,i)*dfermi(j,i)*tau(j,i)*kweight(i)
              K0(m,l)=K0(m,l)+tmp
              K1(m,l)=K1(m,l)+tmp*(eig(j,i)-mu)
              K2(m,l)=K2(m,l)+tmp*(eig(j,i)-mu)*(eig(j,i)-mu)
           end do axis2
           !$omp end simd
        end do axis1
     end do band_loop
  end do get_Kn
  !$omp end do
  !$omp end parallel
end subroutine calc_kn

subroutine calc_sigma_hall(shall,sabs,eig,veloc,imass,kweight,tau,temp,mu,Nk,Norb) bind(C)
  !> calc_sigma_hall
  !> Weak-field (B-linear) Boltzmann Hall kernel for ALL THREE field directions,
  !> band resolved, antisymmetrized in the two current indices:
  !>
  !>   A(a,b,g) = sum_{k,n} w_k tau^2 (-df/de) eps_{g,u,v} v_a M^-1_{b,u} v_v
  !>   S(a,b,g) = 0.5*(A(a,b,g) - A(b,a,g))
  !>
  !> sigma^(1)_{ab}(B||g)/B is S(a,b,g) times the prefactor that turns calc_kn's
  !> K0 into sigma, so the normalization matches calc_kn: a kweight-weighted SUM
  !> over k with no 1/Nk.
  !>
  !> Returning the whole tensor rather than the single B||z scalar is what lets
  !> the caller build rho^(1) = -sigma^-1 sigma^(1) sigma^-1 properly.  Reading
  !> R_H off sigma_xy/(sigma_xx sigma_yy) instead assumes sigma is diagonal in
  !> xyz, which is false for any cell whose axes are not mutually orthogonal --
  !> a monoclinic cell has sigma_xy /= 0 already at B=0.  The band resolution
  !> feeds the two-fluid decomposition needed in a compensated metal, and the
  !> three field directions are what experiment actually varies.
  !>
  !> Contracting at (a,b,g)=(x,y,z) gives, term by term,
  !>   -0.5 w (vx^2 M^-1_yy + vy^2 M^-1_xx - 2 vx vy M^-1_xy),
  !> i.e. minus the textbook one-sided scalar, PROVIDED M^-1 is symmetric (it is,
  !> being d2e/dk_a dk_b).  An asymmetric M^-1 means a non-hermitian H(k) or a
  !> broken get_imassk, and breaks that equivalence.
  !>
  !> sabs accumulates sum |contribution| with the same indexing.  The Hall
  !> integrand changes sign around the Fermi surface, and |S|/sabs is the
  !> cancellation ratio: where it is small the total is a tiny residue of large
  !> opposing parts and R_H is dominated by k-mesh noise, which is what makes its
  !> SIGN move with the mesh.
  !!@param    shall,out: Hall kernel [3,3,3,Norb] (fortran (g,b,a,n))
  !!@param     sabs,out: sum of |contribution|, same indexing
  !!@param       eig,in: energy of bands
  !!@param     veloc,in: group velocity
  !!@param     imass,in: inverse of effective mass
  !!@param   kweight,in: weight of k-points
  !!@param       tau,in: relaxation time
  !!@param      temp,in: Temperature
  !!@param        mu,in: chemical potential
  !!@param        Nk,in: The number of k-points
  !!@param      Norb,in: The number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nk,Norb
  real(c_double),intent(in):: temp,mu
  real(c_double),intent(in),dimension(Norb,Nk):: eig,tau
  real(c_double),intent(in),dimension(Nk):: kweight
  real(c_double),intent(in),dimension(3,Norb,Nk):: veloc
  real(c_double),intent(in),dimension(3,3,Norb,Nk):: imass
  real(c_double),intent(out),dimension(3,3,3,Norb):: shall,sabs

  integer(c_int32_t) i,j,ia,ib,ig,iu,iv
  real(c_double) w,temp_safe,t,cv(3),d(3)
  real(c_double) eps(3,3,3)

  eps(:,:,:)=0.0d0
  eps(1,2,3)=1.0d0; eps(2,3,1)=1.0d0; eps(3,1,2)=1.0d0
  eps(1,3,2)=-1.0d0; eps(3,2,1)=-1.0d0; eps(2,1,3)=-1.0d0
  temp_safe=max(temp,1.0d-12)
  shall(:,:,:,:)=0.0d0
  sabs(:,:,:,:)=0.0d0

  !$omp parallel do private(i,j,ia,ib,ig,iu,iv,w,t,cv,d) reduction(+:shall,sabs)
  k_loop: do i=1,Nk
     band_loop: do j=1,Norb
        ! -df/de = 0.25*(1-tanh^2((e-mu)/2T))/T, times tau^2 and the k-weight
        w=0.25d0*(1.0d0-tanh(0.5d0*(eig(j,i)-mu)/temp_safe)**2)/temp_safe &
             *tau(j,i)*tau(j,i)*kweight(i)
        if(abs(w)<1.0d-300) cycle band_loop
        field_loop: do ig=1,3
           do iu=1,3                       ! cv(u) = eps_{g,u,v} v_v
              cv(iu)=0.0d0
              do iv=1,3
                 cv(iu)=cv(iu)+eps(ig,iu,iv)*veloc(iv,j,i)
              end do
           end do
           do ib=1,3                       ! d(b) = M^-1_{b,u} cv(u)
              d(ib)=0.0d0
              do iu=1,3
                 d(ib)=d(ib)+imass(ib,iu,j,i)*cv(iu)
              end do
           end do
           do ia=1,3
              do ib=1,3
                 t=0.5d0*w*(veloc(ia,j,i)*d(ib)-veloc(ib,j,i)*d(ia))
                 shall(ig,ib,ia,j)=shall(ig,ib,ia,j)+t
                 sabs(ig,ib,ia,j)=sabs(ig,ib,ia,j)+abs(t)
              end do
           end do
        end do field_loop
     end do band_loop
  end do k_loop
  !$omp end parallel do
end subroutine calc_sigma_hall

subroutine calc_k0_band(K0n,eig,veloc,kweight,tau,temp,mu,Nk,Norb) bind(C)
  !> calc_k0_band
  !> Per-band charge transport kernel K0n(i,j,n) = sum_k v_i v_j tau (-df/de) w_k.
  !> calc_kn returns only the band SUM; keeping the bands apart is what allows
  !> the two-fluid decomposition sigma = sum_n sigma_n that explains a Hall
  !> coefficient in a compensated metal.  Sums over n to calc_kn's K0.
  !!@param      K0n,out: per-band charge kernel [3,3,Norb]
  !!@param       eig,in: energy of bands
  !!@param     veloc,in: group velocity
  !!@param   kweight,in: weight of k-points
  !!@param       tau,in: relaxation time
  !!@param      temp,in: Temperature
  !!@param        mu,in: chemical potential
  !!@param        Nk,in: The number of k-points
  !!@param      Norb,in: The number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nk,Norb
  real(c_double),intent(in):: temp,mu
  real(c_double),intent(in),dimension(Norb,Nk):: eig,tau
  real(c_double),intent(in),dimension(Nk):: kweight
  real(c_double),intent(in),dimension(3,Norb,Nk):: veloc
  real(c_double),intent(out),dimension(3,3,Norb):: K0n

  integer(c_int32_t) i,j,l,m
  real(c_double) w,temp_safe

  temp_safe=max(temp,1.0d-12)
  K0n(:,:,:)=0.0d0
  !$omp parallel do private(i,j,l,m,w) reduction(+:K0n)
  k_loop: do i=1,Nk
     band_loop: do j=1,Norb
        w=0.25d0*(1.0d0-tanh(0.5d0*(eig(j,i)-mu)/temp_safe)**2)/temp_safe &
             *tau(j,i)*kweight(i)
        do l=1,3
           do m=1,3
              K0n(m,l,j)=K0n(m,l,j)+veloc(m,j,i)*veloc(l,j,i)*w
           end do
        end do
     end do band_loop
  end do k_loop
  !$omp end parallel do
end subroutine calc_k0_band

subroutine calc_tdf(tdf,eig,veloc,kweight,tau,Nw,Nk,Norb) bind(C)
  !> calc tdf function
  !> sum_k(v_ki*v_kj*tau)
  !!@param    tdf,out: transport distribution function
  !!@param     eig,in: energy of bands
  !!@param   veloc,in: group velocity:
  !!@param kweight,in: weight of k-points
  !!@param     tau,in: relaxation time
  !!@param      Nw,in: The number of energy mesh
  !!@param      Nk,in: The number of k-points
  !!@param    Norb,in: The number of orbitals
  use constant
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nk,Norb,Nw
  real(c_double),intent(in),dimension(Norb,Nk):: eig,tau
  real(c_double),intent(in),dimension(Nk):: kweight
  real(c_double),intent(in),dimension(3,Norb,Nk):: veloc
  real(c_double),intent(out),dimension(3,3,Nw):: tdf

  integer(c_int32_t) i,j,l,m,iw
  real(c_double) tmp,emax,emin,id,dw
  emax=maxval(eig)
  emin=minval(eig)
  dw=(emax-emin)/dble(Nw)
  ! Lorentzian width: must stay wider than the energy grid, otherwise the delta
  ! function falls between the bins and tdf becomes a grid-dependent comb.
  id=max(1.0d-3,2.0d0*dw)
  ! Parallelize over the energy grid: each thread owns distinct iw slices,
  ! so no reduction array is needed (a per-iteration reduction copy of tdf
  ! overflows the stack for large Nw).
  !$omp parallel do private(l,m,i,j,tmp)
  omega_loop: do iw=1,Nw
     axis1: do l=1,3
        axis2: do m=l,3
           tmp=0.0d0
           k_loop: do i=1,Nk
              band_loop: do j=1,Norb
                 ! bin CENTRE E_iw = emin + (iw-1/2)*dw (same binning as get_tau)
                 tmp=tmp+veloc(m,j,i)*veloc(l,j,i)*tau(j,i)*kweight(i) &
                      /(((dble(iw)-0.5d0)*dw+emin-eig(j,i))**2+id*id)
              end do band_loop
           end do k_loop
           ! Lorentzian delta: (1/pi)*id/((E-e)^2+id^2).
           ! NO 1/Nk here: like calc_kn this returns the kweight-weighted SUM, so
           ! that int dE (-df/dE) tdf(E) reproduces calc_kn's K0 and the caller
           ! applies 1/(Nk*Vuc) exactly once.
           tdf(m,l,iw)=tmp*id/pi
           tdf(l,m,iw)=tdf(m,l,iw)
        end do axis2
     end do axis1
  end do omega_loop
  !$omp end parallel do
end subroutine calc_tdf

subroutine get_tau(tau,tauw,eig,tau_max,eps,tau_mode,Nk,Nw,Norb) bind(C)
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,Norb,tau_mode
  real(c_double),intent(in):: eps,tau_max
  real(c_double),intent(in),dimension(Nw):: tauw
  real(c_double),intent(in),dimension(Norb,Nk):: eig
  real(c_double),intent(out),dimension(Norb,Nk):: tau

  integer(c_int32_t) i,j,iter_w
  real(c_double) Emax,Emin,Elength
  Emin=minval(eig(:,:))
  Emax=maxval(eig(:,:))
  Elength=Emax-Emin
  !$omp parallel do private(i,j,iter_w)
  do i=1,Nk
     do j=1,Norb
        if (Elength<1.0d-30) then
           iter_w=1
        else
           iter_w=int(Nw*(eig(j,i)-Emin)/Elength)+1
           iter_w=max(1,min(iter_w,int(Nw)))
        end if
        if(tau_mode==1)then
           if(tauw(iter_w)<eps)then
              tau(j,i)=tau_max
           else
              tau(j,i)=tau_max*eps/tauw(iter_w)
           end if
        else if(tau_mode==2)then
           tau(j,i)=tauw(iter_w)
        end if
     end do
  end do
  !$omp end parallel do
end subroutine get_tau

subroutine calc_tau_epa(tau,gavg,wavg,eig,edge,step,mu,temp,&
     Nk,Norb,nmodes,nbin,ngrid,nbin_max) bind(C)
  !> calc_tau_epa
  !> Compute EPA relaxation time from epa.x output (job='egrid'), following the
  !> energy-conserving EPA rate (Samsonidze & Kozinsky, Adv. Energy Mater. 8,
  !> 1800246 (2018)):
  !>
  !>   1/tau(e) = 2*pi * sum_nu { <|g_nu(e, e+w_nu)|^2> * [nB(w_nu) +   f(e+w_nu)] * rho(e+w_nu)
  !>                            + <|g_nu(e, e-w_nu)|^2> * [nB(w_nu) + 1-f(e-w_nu)] * rho(e-w_nu) }
  !>
  !> The final state is fixed by energy conservation to e +- w_nu (phonon
  !> absorption / emission) and weighted by the electronic DOS rho there.
  !> rho is per spin per unit cell, computed here from ``eig`` by linear
  !> (cloud-in-cell) binning onto the EPA energy bins.  With hbar = 1 and
  !> energies in eV, tau is returned in units of hbar/eV.
  !!@param        tau,out: relaxation time [Norb,Nk] (hbar/eV)
  !!@param       gavg, in: EPA averaged |g|^2 [nmodes,nbin_max,nbin_max,ngrid] (eV^2)
  !!@param       wavg, in: averaged phonon freq per mode [nmodes] (eV)
  !!@param        eig, in: electronic eigenvalues [Norb,Nk] (eV)
  !!@param       edge, in: grid edges [ngrid] (eV)
  !!@param       step, in: grid steps [ngrid] (eV)
  !!@param         mu, in: chemical potential (eV)
  !!@param       temp, in: temperature kB*T (eV)
  !!@param         Nk, in: number of k-points
  !!@param       Norb, in: number of orbitals
  !!@param     nmodes, in: number of phonon modes
  !!@param       nbin, in: number of bins per grid [ngrid]
  !!@param      ngrid, in: number of energy grids (typically 2)
  !!@param   nbin_max, in: max(nbin) — leading dimension for gavg
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nk,Norb,nmodes,ngrid,nbin_max
  integer(c_int64_t),intent(in),dimension(ngrid):: nbin
  real(c_double),intent(in):: mu,temp
  real(c_double),intent(in),dimension(ngrid):: edge,step
  real(c_double),intent(in),dimension(nmodes):: wavg
  real(c_double),intent(in),dimension(nmodes,nbin_max,nbin_max,ngrid):: gavg
  real(c_double),intent(in),dimension(Norb,Nk):: eig
  real(c_double),intent(out),dimension(Norb,Nk):: tau

  integer(c_int32_t) ik,ib,ig,ig_found,jbin,kk,nu,isgn
  real(c_double) eps,gamma,w,nB,ff,xb,xf,g2,ef,dosv,xc,wlo
  real(c_double),parameter:: pi=acos(-1.0d0)
  real(c_double),parameter:: tau_max=1.0d+15
  real(c_double),parameter:: xcut=500.0d0
  real(c_double),dimension(nbin_max,ngrid):: dos

  ! --- electronic DOS per spin per unit cell on the EPA bins -----------------
  ! rho(bin) = (1/Nk) sum_{k,band} delta(e_bin - e_kn), delta -> linear
  ! (cloud-in-cell) weighting between the two adjacent bin centers, so the
  ! energy integral of rho over the grid is preserved.
  dos(:,:)=0.0d0
  do ig=1,int(ngrid)
     do ik=1,int(Nk)
        do ib=1,int(Norb)
           xc=(eig(ib,ik)-edge(ig))/step(ig)-0.5d0   ! fractional bin-center coordinate
           kk=int(floor(xc))+1                       ! lower bin
           wlo=1.0d0-(xc-dble(kk-1))                 ! weight of lower bin
           if(kk>=1 .and. kk<=int(nbin(ig)))   dos(kk,ig)=dos(kk,ig)+wlo
           if(kk+1>=1 .and. kk+1<=int(nbin(ig))) dos(kk+1,ig)=dos(kk+1,ig)+(1.0d0-wlo)
        end do
     end do
  end do
  dos(:,:)=dos(:,:)/dble(Nk)
  do ig=1,int(ngrid)
     dos(:,ig)=dos(:,ig)/step(ig)                    ! states / (eV spin cell)
  end do

  !$omp parallel do private(ik,ib,ig,ig_found,jbin,kk,nu,isgn,eps,gamma,w,nB,ff,xb,xf,g2,ef,dosv)
  do ik=1,Nk
     do ib=1,Norb
        eps=eig(ib,ik)
        ! find which grid & bin this state belongs to
        jbin=-1
        ig_found=-1
        do ig=1,ngrid
           jbin=int((eps-edge(ig))/step(ig))+1
           if(jbin>=1 .and. jbin<=int(nbin(ig)))then
              ig_found=ig
              exit
           end if
           jbin=-1
        end do
        if(jbin<1)then
           tau(ib,ik)=tau_max
           cycle
        end if

        gamma=0.0d0
        do nu=1,nmodes
           w=wavg(nu)
           if(w<1.0d-8) cycle

           ! Bose distribution nB(w)
           xb=w/temp
           if(xb>xcut)then; nB=0.0d0
           else;            nB=1.0d0/(exp(xb)-1.0d0)
           end if

           ! isgn=+1: phonon absorption (final e+w, factor nB+f)
           ! isgn=-1: phonon emission   (final e-w, factor nB+1-f)
           do isgn=1,-1,-2
              ef=eps+dble(isgn)*w
              kk=int((ef-edge(ig_found))/step(ig_found))+1
              if(kk<1 .or. kk>int(nbin(ig_found))) cycle  ! final state outside grid
              g2=gavg(nu,kk,jbin,ig_found)
              if(abs(g2)<1.0d-30) cycle
              dosv=dos(kk,ig_found)
              if(dosv<1.0d-30) cycle

              ! Fermi distribution f(ef)
              xf=(ef-mu)/temp
              if(xf>xcut)then;      ff=0.0d0
              else if(xf<-xcut)then; ff=1.0d0
              else;                  ff=1.0d0/(exp(xf)+1.0d0)
              end if

              if(isgn>0)then
                 gamma=gamma+g2*(nB+ff)*dosv        ! absorption
              else
                 gamma=gamma+g2*(nB+1.0d0-ff)*dosv  ! emission
              end if
           end do
        end do

        gamma=2.0d0*pi*gamma
        if(gamma>1.0d-30)then
           tau(ib,ik)=1.0d0/gamma
        else
           tau(ib,ik)=tau_max
        end if
     end do
  end do
  !$omp end parallel do
end subroutine calc_tau_epa
