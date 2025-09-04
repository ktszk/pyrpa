subroutine calc_lij(L11,L22,L12,vk,eig,ffermi,Norb,Nk,mu,w,idelta,eps,temp) bind(C)
  !> calc_Lij
  !!@param    L11,out: L11@linear response theory
  !!@param    L22,out: L22@linear response theory
  !!@param    L12,out: L12(L21^T)@linear response theory
  !!@param      vk,in: group velocity:
  !!@param     eig,in: energy of bands
  !!@param  ffermi,in: fermi distribute function
  !!@param    Norb,in: The numbrer of orbitals
  !!@param      Nk,in: The number of k-points
  !!@param      mu,in: chemical potential
  !!@param       w,in: frequency (energy)
  !!@param  idelta,in: dumping factor
  !!@param     eps,in: threshold of energy
  !!@param    temp,in: Temperature
  use,intrinsic:: iso_fortran_env, only:int64,real64
  implicit none
  integer(int64),intent(in):: Nk,Norb
  real(real64),intent(in):: temp,eps,w,idelta,mu
  real(real64),intent(in),dimension(Norb,Nk):: eig,ffermi
  complex(real64),intent(in),dimension(3,Norb,Norb,Nk):: vk
  complex(real64),intent(out),dimension(3,3):: L11,L12,L22
  
  integer(int64) i,j,k,l,m
  complex(real64) tmp
  complex(real64),parameter::ii=(0.0d0,1.0d0)

  !$omp parallel
  !$omp workshare
  L11(:,:)=0.0d0
  L12(:,:)=0.0d0
  L22(:,:)=0.0d0
  !$omp end workshare
  !$omp do reduction(+: L11,L12,L22) private(i,l,m,j,k,tmp)
  k_loop: do i=1,Nk
     band_loop1: do l=1,Norb
        band_loop2: do m=1,Norb
           do j=1,3
              do k=1,3
                 if(abs(eig(m,i)-eig(l,i))<1.0d-9)then
                    tmp=vk(k,m,m,i)*vk(j,m,m,i)*ffermi(m,i)*(1.0d0-ffermi(m,i))/(temp*cmplx(w,idelta))
                    L11(k,j)=L11(k,j)+tmp
                    L12(k,j)=L12(k,j)+tmp*(eig(m,i)-mu)
                    L22(k,j)=L22(k,j)+tmp*(eig(m,i)-mu)*(eig(m,i)-mu)
                 else if(abs(ffermi(l,i)-ffermi(m,i))>eps)then
                    tmp=vk(k,m,l,i)*vk(j,l,m,i)*(ffermi(l,i)-ffermi(m,i))/((eig(m,i)-eig(l,i))&
                         *cmplx(w+eig(m,i)-eig(l,i),idelta))
                    L11(k,j)=L11(k,j)+tmp
                    L12(k,j)=L12(k,j)+tmp*(eig(l,i)-mu)
                    L22(k,j)=L22(k,j)+tmp*(eig(m,i)-mu)*(eig(l,i)-mu)
                 end if
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
  !!@param    Norb,in: The numbrer of orbitals
  use,intrinsic:: iso_fortran_env, only:int64,real64
  implicit none
  integer(int64),intent(in):: Nk,Norb
  real(real64),intent(in):: temp,mu
  real(real64),intent(in),dimension(Norb,Nk):: eig,tau
  real(real64),intent(in),dimension(Nk):: kweight
  real(real64),intent(in),dimension(3,Norb,Nk):: veloc
  real(real64),intent(out),dimension(3,3):: K0,K1,K2

  real(real64),dimension(Norb,Nk):: dfermi
  integer(int64) i,j,l,m
  real(real64) tmp

  !$omp parallel
  !$omp do private(j)
  get_dfermi: do i=1,Nk
     do j=1,Norb
        dfermi(j,i)=0.25d0*(1.0d0-tanh(0.5d0*(eig(j,i)-mu)/temp)**2)/temp
     end do
  end do get_dfermi
  !$omp end do
  !$omp workshare
  K0(:,:)=0.0d0
  K1(:,:)=0.0d0
  K2(:,:)=0.0d0
  !$omp end workshare
  
  !$omp do private(j,l,m,tmp) reduction(+:K0,K1,K2)
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
  !!@param    Norb,in: The numbrer of orbitals
  use,intrinsic:: iso_fortran_env, only:int64,real64
  implicit none
  integer(int64),intent(in):: Nk,Norb,Nw
  real(real64),intent(in),dimension(Norb,Nk):: eig,tau
  real(real64),intent(in),dimension(Nk):: kweight
  real(real64),intent(in),dimension(3,Norb,Nk):: veloc
  real(real64),intent(out),dimension(3,3,Nw):: tdf

  integer(int64) i,j,l,m,iw
  real(real64) tmp,emax,emin,id,dw
  id=1.0d-3
  emax=maxval(eig)
  emin=minval(eig)
  dw=(emax-emin)/Nw
  !$omp parallel
  !$omp workshare
  tdf(:,:,:)=0.0d0
  !$omp end workshare
  omega_loop: do iw=1,Nw
     axis1: do l=1,3
        axis2: do m=l,3
           !$omp do private(i,j) reduction(+:tdf)
           k_loop: do i=1,Nk
              band_loop: do j=1,Norb
                 tdf(m,l,iw)=tdf(m,l,iw)+veloc(m,j,i)*veloc(l,j,i)*tau(j,i)*kweight(i)/((iw*dw+emin-eig(j,i))**2+id*id)
              end do band_loop
           end do k_loop
           !$omp end do
           tdf(l,m,iw)=tdf(m,l,iw)
        end do axis2
     end do axis1
  end do omega_loop
  !$omp workshare
  tdf(:,:,:)=tdf(:,:,:)*id/Nk
  !$omp end workshare
  !$omp end parallel
end subroutine calc_tdf
