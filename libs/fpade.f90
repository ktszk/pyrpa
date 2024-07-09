subroutine get_a(a,xn,inp_data,Np) bind(C)
  !
  !> calculate pade's constants a
  !
  use,intrinsic:: iso_fortran_env, only:int32,int64,real64
  implicit none
  integer(int64),intent(in):: Np !length of data point
  complex(real64),intent(in),dimension(Np):: xn !variables of data
  complex(real64),intent(in),dimension(Np):: inp_data !data
  complex(real64),intent(out),dimension(Np):: a

  integer(int32) i,j
  complex(real64),dimension(Np):: g1,g0

  !$omp parallel workshare
  g0(:)=inp_data(:) !g_1=inp_data
  a(:)=0.0d0
  !$omp end parallel workshare
  a(1)=g0(1)  !a_1=g_1(z1)
  do i=2,Np !g0=g_(i-1),g1=g_i
     !$omp parallel
     !$omp workshare
     g1(:)=0.0d0
     !$omp end workshare
     !$omp do private(j)
     do j=1,Np !calc all g1 at data point
        if(g0(j)==0)then
           continue
        else
           g1(j)=(g0(i-1)-g0(j))/((xn(j)-xn(i-1))*g0(j))
        end if
     end do
     !$omp end do
     !$omp workshare
     g0(:)=g1(:) !renew g0->g_i
     !$omp end workshare
     !$omp end parallel
     a(i)=g1(i) !a_i=g_i(z_i)
  end do

  return
end subroutine get_a

subroutine get_QP(P,Q,a,xn,wlist,Nw,Np) bind(C)
  !
  !>  calculate rational function P,Q for pade
  !
  use,intrinsic:: iso_fortran_env, only:int32,int64,real64
  implicit none
  integer(int64),intent(in):: Nw,Np
  complex(real64),intent(in),dimension(Np):: xn,a
  complex(real64),intent(in),dimension(Nw):: wlist
  complex(real64),intent(out),dimension(Nw):: P,Q !output P_N,Q_N

  integer(int32) i,j,k
  complex(real64),dimension(Nw):: P0,P1,Q0,Q1

  Q0(:)=1.0d0
  Q1(:)=1.0d0
  P0(:)=0.0d0
  P1(:)=a(1)
  do i=2,Np
     !$omp parallel
     !$omp workshare
     Q(:)=0.0d0
     P(:)=0.0d0
     !$omp end workshare
     !$omp do private(j) reduction(+: P,Q)
     do j=1,Nw !calc Pi(w)=P_i-1(w)+a_i(w-z_i-1)P_i-2(w) (Q is same as P)
        P(j)=P1(j)+a(i)*(wlist(j)-xn(i-1))*P0(j)
        Q(j)=Q1(j)+a(i)*(wlist(j)-xn(i-1))*Q0(j)
     end do
     !$omp end do
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
