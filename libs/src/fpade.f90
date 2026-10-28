subroutine get_a(a,xn,inp_data,Np)
  !> calculate pade's constants a
  !!@param       a,out:
  !!@param       xn,in: variables of origin data
  !!@param inp_data,in: The origin data
  !!@param       Np,in: The number of origin data
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Np !length of data point
  complex(c_double),intent(in),dimension(Np):: xn !variables of data
  complex(c_double),intent(in),dimension(Np):: inp_data !data
  complex(c_double),intent(out),dimension(Np):: a

  integer(c_int32_t) i,j
  complex(c_double),dimension(Np):: g1,g0

  ! Thiele's reciprocal difference algorithm for Padé continued fraction coefficients
  ! g_{i}(z) = (g_{i-1}(z_{i-1}) - g_{i-1}(z)) / ((z - z_{i-1}) * g_{i-1}(z))
  g0(:)=inp_data(:) !g_1=inp_data
  !$omp parallel workshare
  a(:)=0.0d0
  !$omp end parallel workshare
  a(1)=g0(1)  !a_1=g_1(z1)
  do i=2,Np !g0=g_(i-1), g1=g_i
     !$omp parallel
     !$omp do simd private(j)
     do j=1,Np !evaluate g_i at all data points
        ! Guard against g0(j)=0 to avoid division by zero; merge selects zero branch
        g1(j)=merge((0.0d0,0.0d0),(g0(i-1)-g0(j))/ &
             merge((1.0d0,0.0d0),(xn(j)-xn(i-1))*g0(j),g0(j)==(0.0d0,0.0d0)), &
             g0(j)==(0.0d0,0.0d0))
     end do
     !$omp end do simd
     !$omp workshare
     g0(:)=g1(:) !advance: g_{i-1} <- g_i
     !$omp end workshare
     !$omp end parallel
     a(i)=g1(i) !a_i = g_i(z_i), the i-th continued fraction coefficient
  end do

  return
end subroutine get_a

subroutine get_QP(P,Q,a,xn,wlist,Nw,Np)
  !>  calculate rational function P,Q for pade
  !!@param    P,out:
  !!@param    Q,out:
  !!@param     a,in:
  !!@param    xn,in: variables of origin data
  !!@param wlist,in: new variables of pade data
  !!@param    Nw,in: The number of pade data
  !!@param    Np,in: The number of origin data
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nw,Np
  complex(c_double),intent(in),dimension(Np):: xn,a
  complex(c_double),intent(in),dimension(Nw):: wlist
  complex(c_double),intent(out),dimension(Nw):: P,Q !output P_N,Q_N

  integer(c_int32_t) i,j
  complex(c_double),dimension(Nw):: P0,P1,Q0,Q1

  Q0(:)=1.0d0
  Q1(:)=1.0d0
  P0(:)=0.0d0
  P1(:)=a(1)
  do i=2,Np
     !$omp parallel
     !$omp do simd private(j)
     do j=1,Nw !calc Pi(w)=P_i-1(w)+a_i(w-z_i-1)P_i-2(w) (Q is same as P)
        P(j)=P1(j)+a(i)*(wlist(j)-xn(i-1))*P0(j)
        Q(j)=Q1(j)+a(i)*(wlist(j)-xn(i-1))*Q0(j)
     end do
     !$omp end do simd
     !$omp workshare
     P0(:)=P1(:)
     P1(:)=P(:)
     Q0(:)=Q1(:)
     Q1(:)=Q(:)
     !$omp end workshare
     !$omp end parallel
  end do

  return
end subroutine get_QP

subroutine pade_analytic_continuation_arrays(arrayin,arrayout,&
     iwlist,wlist,Nk,Niw,Nw)
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,Niw
  complex(c_double),intent(in),dimension(Nk,Niw)::arrayin
  complex(c_double),intent(in),dimension(Nw):: wlist
  complex(c_double),intent(in),dimension(Niw):: iwlist
  complex(c_double),intent(out),dimension(Nk,Nw):: arrayout

  integer(c_int32_t) i
  complex(c_double),dimension(Niw):: a
  complex(c_double),dimension(Nw):: P,Q

  do i=1,Nk
     call get_a(a,iwlist,arrayin(i,:),Niw)
     call get_QP(P,Q,a,iwlist,wlist,Nw,Niw)
     arrayout(i,:)=P(:)/Q(:)
  end do
end subroutine pade_analytic_continuation_arrays

subroutine pade_with_trace(arrayin,arrayout,iwlist,wlist,Nk,Niw,Nw,Norb) bind(C)
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,Niw,Norb
  complex(c_double),intent(in),dimension(Nk,Niw,Norb,Norb)::arrayin
  complex(c_double),intent(in),dimension(Nw):: wlist
  complex(c_double),intent(in),dimension(Niw):: iwlist
  complex(c_double),intent(out),dimension(Nw,Nk):: arrayout

  integer(c_int32_t) i,j,k
  complex(c_double),allocatable,dimension(:,:):: tmp

  allocate(tmp(Nk,Nw))
  arrayout(:,:)=(0.0d0,0.0d0)  ! accumulated below; do not rely on the caller zeroing it
  do i=1,Norb
     call pade_analytic_continuation_arrays(arrayin(:,:,i,i),tmp,iwlist,wlist,Nk,Niw,Nw)
     do j=1,Nk
        do k=1,Nw
           arrayout(k,j)=arrayout(k,j)+tmp(j,k)
        end do
     end do
  end do
end subroutine pade_with_trace

subroutine pade_batch(arrayout,arrayin,iwlist,zlist,Nd,Niw,Nz) bind(C)
  !> Pade continuation of Nd independent functions, each evaluated at its OWN
  !> Nz points (zlist(:,d)).  Same Thiele continued fraction as get_a/get_QP,
  !> but the recursion runs serially per function and OpenMP is over the Nd
  !> functions: get_a/get_QP open a parallel region per recursion step, which
  !> costs far more than the O(Niw^2) work of one small fit.
  !!@param arrayout,out: continued values f_d(zlist(:,d)) [Nz,Nd]
  !!@param   arrayin,in: data f_d(iwlist) [Niw,Nd]
  !!@param    iwlist,in: Matsubara points [Niw]
  !!@param     zlist,in: evaluation points per function [Nz,Nd]
  !!@param        Nd,in: The number of functions
  !!@param       Niw,in: The number of Matsubara points
  !!@param        Nz,in: The number of evaluation points per function
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nd,Niw,Nz
  complex(c_double),intent(in),dimension(Niw,Nd):: arrayin
  complex(c_double),intent(in),dimension(Niw):: iwlist
  complex(c_double),intent(in),dimension(Nz,Nd):: zlist
  complex(c_double),intent(out),dimension(Nz,Nd):: arrayout

  integer(c_int64_t) d,i,j,nterm
  complex(c_double),dimension(Niw):: a,g0,g1
  complex(c_double) p0,p1,q0,q1,p,q

  !$omp parallel do private(d,i,j,nterm,a,g0,g1,p0,p1,q0,q1,p,q)
  do d=1,Nd
     ! continued-fraction coefficients (Thiele reciprocal differences, as get_a).
     ! Unlike get_a the fraction is TRUNCATED at the first vanishing coefficient
     ! (|a_i z_i| <= 1e-12, a_i*z is dimensionless for i>=2): data that is exactly
     ! a low-order rational function (a constant or linear Sigma, the CPA/impurity
     ! limits) otherwise drives the later reciprocal differences through 0/0 and
     ! the round-off garbage behind it turns the continuation into NaN.
     g0(:)=arrayin(:,d)
     a(:)=0.0d0
     a(1)=g0(1)
     nterm=Niw
     do i=2,Niw
        do j=1,Niw
           if(g0(j)==(0.0d0,0.0d0))then
              g1(j)=0.0d0
           else
              g1(j)=(g0(i-1)-g0(j))/((iwlist(j)-iwlist(i-1))*g0(j))
           end if
        end do
        g0(:)=g1(:)
        a(i)=g1(i)
        if(abs(a(i)*iwlist(i))<=1.0d-12)then
           nterm=i-1
           exit
        end if
     end do
     ! P_N/Q_N three-term recursion (as get_QP)
     do j=1,Nz
        q0=1.0d0; q1=1.0d0
        p0=0.0d0; p1=a(1)
        p=p1; q=q1
        do i=2,nterm
           p=p1+a(i)*(zlist(j,d)-iwlist(i-1))*p0
           q=q1+a(i)*(zlist(j,d)-iwlist(i-1))*q0
           p0=p1; p1=p
           q0=q1; q1=q
        end do
        arrayout(j,d)=p/q
     end do
  end do
  !$omp end parallel do
end subroutine pade_batch
