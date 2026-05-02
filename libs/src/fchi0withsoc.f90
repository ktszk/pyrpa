subroutine get_chi0_conv_soc(chi,Gk,kmap,invk,invs,irr_chi,chi_map,olist,&
     sgnsig,sgnsig2,temp,Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi)
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nw,Norb,Nchi,Nkall,Nk,Nx,Ny,Nz
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),dimension(Norb):: invs
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
  integer(c_int32_t),intent(in),dimension(Nchi,Nchi,2):: chi_map
  integer(c_int32_t),intent(in),dimension(Nchi*(Nchi+1)/2,2):: irr_chi
  real(c_double),intent(in):: temp
  real(c_double),intent(in),dimension(Norb,Norb):: sgnsig
  real(c_double),intent(in),dimension(Nchi,Nchi):: sgnsig2
  complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(c_double),intent(out),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(c_int32_t) i,j,k,l,m,n,iorb
  integer(c_int32_t) ii(0:Nx-1),ij(0:Ny-1),ik(0:Nz-1),iw(2*Nw)
  real(c_double) weight
  real(c_double),parameter:: eps=1.0d-9
  complex(c_double),dimension(0:Nx-1,0:Ny-1,0:Nz-1,2*Nw):: tmp,tmpgk13,tmpgk42
   
  weight=temp/dble(Nkall)
  ii(0)=0
  ij(0)=0
  ik(0)=0
  iw(1)=1
  !$omp parallel
  !$omp do
  do i=1,Nx-1
     ii(i)=Nx-i
  end do
  !$omp end do
  !$omp do
  do i=1,Ny-1
     ij(i)=Ny-i
  end do
  !$omp end do
  !$omp do
  do i=1,Nz-1
     ik(i)=Nz-i
  end do
  !$omp end do
  !$omp do
  do i=2,2*Nw
     iw(i)=2*Nw-i+2
  end do
  !$omp end do
  !$omp end parallel
  orb_loop:do iorb=1,Nchi*(Nchi+1)/2
     l=irr_chi(iorb,1)
     m=irr_chi(iorb,2)
     !use symmetry G^lm(k,iw)=G^ml(k,-iw) from Hermitian symmetry of Hamiltonian
     !$omp parallel do private(i,j)
     w_loop_Gk_to_tmp:do j=1,Nw
        k_loop_Gk_to_tmp:do i=1,Nkall
           if(invk(2,i)==0)then !irredicible k
              !iw
              tmpgk13(kmap(1,i),kmap(2,i),kmap(3,i),j)=Gk(invk(1,i),j,olist(l,1),olist(m,1)) !G13(k,iw)
              tmpgk42(kmap(1,i),kmap(2,i),kmap(3,i),j)=Gk(invk(1,i),j,olist(m,2),olist(l,2)) !G42(k,iw)
              !-iw
              tmpgk13(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=conjg(Gk(invk(1,i),j,olist(m,1),olist(l,1))) !G13(k,-iw)=G^*31(k,iw)
              tmpgk42(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=conjg(Gk(invk(1,i),j,olist(l,2),olist(m,2))) !G42(k,-iw)=G^*24(k,iw)
           else if(invk(2,i)==1)then !-k
              !iw
              tmpgk13(kmap(1,i),kmap(2,i),kmap(3,i),j)=sgnsig(olist(m,1),olist(l,1))&
                   *Gk(invk(1,i),j,invs(olist(m,1)),invs(olist(l,1))) !G13ss'(-k,iw)=ss'G^31-s'-s(k,iw)
              tmpgk42(kmap(1,i),kmap(2,i),kmap(3,i),j)=sgnsig(olist(l,2),olist(m,2))&
                   *Gk(invk(1,i),j,invs(olist(l,2)),invs(olist(m,2))) !G42ss'(-k,iw)=ss'G^24-s's(k,iw)
              !-iw
              tmpgk13(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=sgnsig(olist(l,1),olist(m,1))&
                   *conjg(Gk(invk(1,i),j,invs(olist(l,1)),invs(olist(m,1)))) !G13ss'(-k,-iw)=ss'G^*13-s-s'(k,iw)
              tmpgk42(kmap(1,i),kmap(2,i),kmap(3,i),2*Nw-j+1)=sgnsig(olist(m,2),olist(l,2))&
                   *conjg(Gk(invk(1,i),j,invs(olist(m,2)),invs(olist(l,2)))) !G42(-k,-iw)=ss'G^*42-s-s'(k,iw)
           end if
        end do k_loop_Gk_to_tmp
     end do w_loop_Gk_to_tmp
     !$omp end parallel do
     call FFT(tmpgk13,tmp,Nx,Ny,Nz,2*Nw,.true.)
     call FFT(tmpgk42,tmp,Nx,Ny,Nz,2*Nw,.true.)
     !calculate G(r)G(-r)
     !$omp parallel do private(i,j,k)
     w_loop_conv:do n=1,2*Nw
        z_loop:do k=0,Nz-1
           y_loop:do j=0,Ny-1
              !$omp simd
              x_loop:do i=0,Nx-1
                 tmp(i,j,k,n)=-tmpgk13(i,j,k,n)*tmpgk42(ii(i),ij(j),ik(k),iw(n))
              end do x_loop
              !$omp end simd
           end do y_loop
        end do z_loop
     end do w_loop_conv
     !$omp end parallel do
     call FFT(tmp,tmpgk13,Nx,Ny,Nz,2*Nw,.false.)
     !$omp parallel do private(i,j)
     w_loop_tmp_to_chi:do j=1,Nw
        k_loop_tmp_to_chi:do i=1,Nkall
           if(invk(2,i)==0)then
              chi(invk(1,i),j,m,l)=tmp(kmap(1,i),kmap(2,i),kmap(3,i),j)*weight
              if(abs(dble(chi(invk(1,i),j,m,l)))<eps) chi(invk(1,i),j,m,l)=cmplx(0.0d0,imag(chi(invk(1,i),j,m,l)),kind=c_double)
              if(abs(imag(chi(invk(1,i),j,m,l)))<eps) chi(invk(1,i),j,m,l)=cmplx(dble(chi(invk(1,i),j,m,l)),0.0d0,kind=c_double)
           end if
        end do k_loop_tmp_to_chi
     end do w_loop_tmp_to_chi
     !$omp end parallel do
     chi(:,:,chi_map(m,l,1),chi_map(m,l,2))=sgnsig2(m,l)*conjg(chi(:,:,m,l))
  end do orb_loop
end subroutine get_chi0_conv_soc

subroutine get_chi0_soc(chi,sgnsig,sgnsig2,invschi,Vmat,Gk,kmap,invk,invs,olist,slist,temp,&
     Nx,Ny,Nz,Nw,Nk,Nkall,Nchi,Norb) bind(C)
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nkall,Nk,Nw,Nchi,Norb,Nx,Ny,Nz
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),dimension(Norb):: slist,invs
  integer(c_int64_t),intent(in),dimension(3,Nkall):: kmap,invk
  integer(c_int64_t),intent(out),dimension(Nchi):: invschi
  real(c_double),intent(in):: temp
  real(c_double),intent(in),dimension(Nchi,Nchi):: Vmat
  real(c_double),intent(out),dimension(Norb,Norb):: sgnsig
  real(c_double),intent(out),dimension(Nchi,Nchi):: sgnsig2
  complex(c_double),intent(in),dimension(Nk,Nw,Norb,Norb):: Gk
  complex(c_double),intent(out),dimension(Nk,Nw,Nchi,Nchi):: chi

  integer(c_int32_t),dimension(Nchi,Nchi,2)::chi_map
  integer(c_int32_t),dimension(Nchi*(Nchi+1)/2,2)::irr_chi

  call get_invschi()
  call get_sgnsig()
  call get_chi_map_soc(chi_map,irr_chi,olist,invs,Nchi,Norb)
  call get_chi0_conv_soc(chi,Gk,kmap,invk,invs,irr_chi,chi_map,olist,sgnsig,sgnsig2,temp,Nx,Ny,Nz,Nw,Nk,Nkall,Norb,Nchi)
  call ckchi()
contains
  subroutine get_sgnsig()
    ! sgnsig(i,j)  = s_i * s_j  (orbital spin signature for TRS G(-k) relation)
    ! sgnsig2(i,j) = s_{i1}*s_{i2}*s_{j1}*s_{j2}  (chi index spin signature)
    integer(c_int32_t) i,j
    do j=1,Norb
       do i=j,Norb
          sgnsig(i,j)=slist(i)*slist(j)
          sgnsig(j,i)=sgnsig(i,j)
       end do
    end do

    do j=1,Nchi
       do i=j,Nchi
          sgnsig2(i,j)=slist(olist(i,1))*slist(olist(i,2))&
               *slist(olist(j,1))*slist(olist(j,2))
          sgnsig2(j,i)=sgnsig2(i,j)
       end do
    end do
  end subroutine get_sgnsig

    subroutine ckchi()
    integer(c_int32_t) i,l,m,n,info,chik,chikall,iorb
    integer(c_int32_t) Nchi_i,lwork
    real(c_double) maxchi0,maxchi02
    real(c_double),dimension(2*Nchi):: rwork
    complex(c_double),dimension(Nchi*Nchi*4+1):: work
    complex(c_double),dimension(Nchi):: eigc
    complex(c_double),dimension(Nchi,Nchi):: chi0,tmp,tmp1

    Nchi_i=int(Nchi,c_int32_t)
    lwork=int(Nchi*Nchi*4+1,c_int32_t)
    maxchi02=-1.0d5
    do i=1,Nk
       !$omp parallel
       !$omp workshare
       chi0(:,:)=0.0d0
       !$omp end workshare
       !$omp do private(l,m,n)
       do l=1,Nchi
          do m=1,Nchi
             do n=1,Nchi
                chi0(m,l)=chi0(m,l)-chi(i,1,m,n)*Vmat(n,l)
             end do
          end do
       end do
       !$omp end do
       !$omp end parallel
       call zgeev('N','N',Nchi_i,chi0,Nchi_i,eigc,tmp1,Nchi_i,tmp,Nchi_i,work,lwork,rwork,info)
       if(info/=0)then; print*,'zgeev failed: info=',info; stop; end if
       maxchi0=maxval(dble(eigc))
       if(maxchi0>maxchi02)then
          chik=i
          maxchi02=maxchi0
       end if
    end do
    do i=1,Nkall
       if(invk(2,i)==0 .and. invk(1,i)==chik)then
          chikall=i
          exit
       end if
    end do
    print'(A7,3I4,F12.8)','SDW/CDW',kmap(:,chikall),maxchi02
  end subroutine ckchi

  subroutine get_invschi()
    integer(c_int32_t) l,m,i,j
    invschi(:)=0
    do l=1,Nchi
       i=invs(olist(l,1))
       j=invs(olist(l,2))
       do m=1,Nchi
          if(i==olist(m,1) .and. j==olist(m,2))then
             invschi(l)=m
             exit
          end if
       end do
    end do
  end subroutine get_invschi
end subroutine get_chi0_soc

subroutine get_chi_map_soc(chi_map,irr_chi,olist,invs,Nchi,Norb)
  !> This function generate index of exchange symmetry chi1234(q,iw)=sgn*chi*4321(q,iw).
  !> This symmetry can use system has TRS iwth soc.
  !!@param chi_map,out: mapping list of chi index
  !!@param irr_chi,out: irreducible index of chii
  !!@param    olist,in: orbital index list of chi index
  !!@param     invs,in: spin flipped index of Hamiltonian
  !!@param     Nchi,in: Number of chi index
  !!@param     Norb,in: Number of orbitals
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nchi,Norb
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),dimension(Norb):: invs
  integer(c_int32_t),intent(out),dimension(Nchi,Nchi,2):: chi_map
  integer(c_int32_t),intent(out),dimension(Nchi*(Nchi+1)/2,2):: irr_chi

  integer(c_int32_t) l1,m1,l2,m2,iter
  integer(c_int32_t),dimension(Nchi,Nchi):: chi_irr
  integer(c_int32_t),dimension(4):: tmp1,tmp2
  logical ck

  chi_irr(:,:)=0
  chi_map(:,:,:)=0
  do l1=1,Nchi
     tmp1(1)=olist(l1,1) !1 of 1234
     tmp1(2)=olist(l1,2) !2 of 1234
     do m1=1,Nchi
        tmp1(3)=olist(m1,1) !3 of 1234
        tmp1(4)=olist(m1,2) !4 of 1234
        ck=.false.
        do l2=1,Nchi
           tmp2(3)=invs(olist(l2,2)) !2 of 4321 with spin flip
           tmp2(4)=invs(olist(l2,1)) !1 of 4321 with spin flip
           do m2=1,Nchi
              tmp2(1)=invs(olist(m2,2)) !4 of 4321 with spin flip
              tmp2(2)=invs(olist(m2,1)) !3 of 4321 with spin flip
              if(sum(abs(tmp1(:)-tmp2(:)))==0)then !4321 correspond to 1234
                 chi_map(m1,l1,1)=m2
                 chi_map(m1,l1,2)=l2
                 if(chi_irr(m2,l2)==0)then !4321 is not irreducible
                    chi_irr(m1,l1)=1 !1 is irreducible index
                 end if
                 ck=.true.
                 exit
              end if
           end do
           if(ck)exit
        end do
     end do
  end do

  !get list of irreducible chi index
  iter=1
  do l1=1,Nchi
     do m1=1,Nchi
        if(chi_irr(m1,l1)==1)then
           irr_chi(iter,1)=l1
           irr_chi(iter,2)=m1
           iter=iter+1
        end if
     end do
  end do
end subroutine get_chi_map_soc

subroutine get_chis_chic_soc(chic,chiszz,chispm,chi,Vmat,orb_list,olist,slist,invs,Nk,Nw,Nchi,Norb) bind(C)
  use,intrinsic:: iso_c_binding, only:c_int64_t,c_double,c_int32_t
  implicit none
  integer(c_int64_t),intent(in):: Nk,Nw,Nchi,Norb
  integer(c_int64_t),intent(in),dimension(Nchi,2):: olist
  integer(c_int64_t),intent(in),dimension(Norb):: slist,invs
  integer(c_int64_t),dimension(Nchi):: orb_list
  real(c_double),intent(in),dimension(Nchi,Nchi):: Vmat
  complex(c_double),intent(in),dimension(Nk,Nw,Nchi,Nchi):: chi
  complex(c_double),intent(out),dimension(Nk,Nchi/4,Nchi/4):: chiszz,chic,chispm

  integer(c_int32_t) i,l,m,info
  integer(c_int32_t),dimension(Nchi):: ipiv
  complex(c_double),dimension(Nchi,Nchi):: cmat1,cmat2,cmatv

  cmatv(:,:)=cmplx(Vmat(:,:),0.0d0,kind=c_double)

  qloop:do i=1,Nk
     call zgemm('N','N',Nchi,Nchi,Nchi,(1.0d0,0.0d0),chi(i,1,:,:),Nchi,cmatv,Nchi,(0.0d0,0.0d0),cmat1,Nchi) !chi0V
     cmat2(:,:)=chi(i,1,:,:) !RHS for zgesv = chi0
     do l=1,Nchi
        cmat1(l,l)=cmat1(l,l)+1.0d0 !I-chi0V
     end do
     call zgesv(Nchi,Nchi,cmat1,Nchi,ipiv,cmat2,Nchi,info) !(I-chi0V)X=chi0 -> cmat2=(I-chi0V)^-1chi0
     if(info/=0)then; print*,'zgesv failed: info=',info; stop; end if
     do l=1,Nchi
        do m=1,Nchi
           if(slist(olist(l,1))==slist(olist(m,1)) .and. slist(olist(l,2))==slist(olist(m,2)))then
              if(slist(olist(l,1))==-slist(olist(l,2)) .and. slist(olist(l,1))==1)then
                 chispm(i,orb_list(l),orb_list(m))=chispm(i,orb_list(l),orb_list(m))+cmat2(m,l)
              end if
           end if
           if(slist(olist(l,1))==slist(olist(l,2)) .and. slist(olist(m,1))==slist(olist(m,2)))then
              chic(i,orb_list(l),orb_list(m))=chic(i,orb_list(l),orb_list(m))+cmat2(m,l)
              if(slist(olist(l,1))==slist(olist(m,1)))then
                 chiszz(i,orb_list(l),orb_list(m))=chiszz(i,orb_list(l),orb_list(m))+cmat2(m,l)
              else
                 chiszz(i,orb_list(l),orb_list(m))=chiszz(i,orb_list(l),orb_list(m))-cmat2(m,l)
              end if
           end if
        end do
     end do
  end do qloop
  chic(:,:,:)=chic(:,:,:)*0.5d0
  chiszz(:,:,:)=chiszz(:,:,:)*0.25
end subroutine get_chis_chic_soc
