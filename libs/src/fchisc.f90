subroutine mkBdGhamk(hamBdGk,hamk,delta,Nk,Norb)
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
  integer(c_int64_t),intent(in):: Nk,Norb
  complex(c_double),intent(in),dimension(Norb,Norb,Nk):: hamk,delta
  complex(c_double),intent(out),dimension(2*Norb,2*Norb,Nk):: hamBdGk

  integer(c_int32_t) k,l,m

  do k=1,Nk
     do l=1,Norb
        do m=1,Norb
           hamBdGk(l,m,k)=hamk(l,m,k)
           hamBdGk(l+Norb,m+Norb,k)=-conjg(hamk(l,m,k))
           hamBdGk(l,m+Norb,k)=delta(l,m,k)
           hamBdGk(l+Norb,m,k)=conjg(delta(m,l,k))
        end do
     end do
  end do
end subroutine mkBdGhamk

module calc_irr_chi_sc
  use,intrinsic:: iso_c_binding, only:c_int32_t,c_int64_t,c_double
  implicit none
contains
  function irr_chi_sc(Nk,Norb,Nchi,uni,eig,ffermi,ol,temp,qshift,w,idelta,eps,sw_spsym)
    !!@param        Nk: The number of k-points
    !!@param      Norb: The number of orbitals
    !!@param      Nchi: The footnote of chi
    !!@param       uni: unitary matrix
    !!@param       eig: energies of bands
    !!@param    ffermi: fermi distribute function
    !!@param        ol: the list of the properties of orbitals at footnote of chi
    !!@param      temp: temperature
    !!@param    qshift: q-shifted klist
    !!@param         w: frequency
    !!@param    idelta: dumping factor
    !!@param       eps: threshold of calculation value
    !!@param   sw_spsym: symmetry flag
    !!@return irr_chi_sc: irreducible susceptibility matrix
    integer(c_int64_t),intent(in):: Nk,Norb,Nchi
    integer(c_int64_t),intent(in),dimension(Nk):: qshift
    integer(c_int64_t),intent(in),dimension(Nchi,2):: ol
    real(c_double),intent(in):: temp,eps,idelta,w
    logical,intent(in):: sw_spsym
    real(c_double),intent(in),dimension(2*Norb,Nk):: eig,ffermi
    complex(c_double),intent(in),dimension(2*Norb,2*Norb,Nk):: uni

   integer(c_int32_t) i,k,l,m,nchi32
   real(c_double) temp_safe,w_eps,sgn
   complex(c_double) weight
   complex(c_double),dimension(Nchi):: A_vec,B_vec,P_vec,R_vec
   complex(c_double),dimension(Nchi,Nchi):: chi,irr_chi_sc

   if(sw_spsym)then
     sgn=-1.0d0 !triplet_dz
   else
     sgn=+1.0d0 !singlet
   end if
   temp_safe=max(temp,1.0d-12)
   w_eps=1.0d-12
   nchi32=int(Nchi,c_int32_t)
   chi(:,:)=0.0d0
   !$omp parallel do reduction(+:chi) private(l,m,i,weight,A_vec,B_vec,P_vec,R_vec)
   kloop: do k=1,Nk
      band1_loop: do l=1,2*Norb
         band2_loop: do m=1,2*Norb
            if(abs(w)<w_eps .and. abs(eig(m,k)-eig(l,qshift(k)))<1.0d-9)then
               weight=cmplx(ffermi(m,k)*(1.0d0-ffermi(m,k))/temp_safe,0.0d0,kind=c_double)
            else if(abs(ffermi(l,qshift(k))-ffermi(m,k))>eps)then
               weight=(ffermi(l,qshift(k))-ffermi(m,k))&
                    /cmplx(w+eig(m,k)-eig(l,qshift(k)),idelta,kind=c_double)
            else
               cycle band2_loop
            end if
            do i=1,nchi32
               ! uni is stored as uni(orbital,band,k)
               A_vec(i)=uni(ol(i,1),l,qshift(k))*conjg(uni(ol(i,2),m,k))
               B_vec(i)=conjg(uni(ol(i,1),l,qshift(k)))*uni(ol(i,2),m,k)
               P_vec(i)=uni(ol(i,1),l,qshift(k))*conjg(uni(ol(i,2)+Norb,m,k))
               R_vec(i)=uni(ol(i,2)+Norb,l,qshift(k))*conjg(uni(ol(i,1),m,k))
            end do
            call zgeru(nchi32,nchi32,weight,B_vec,1,A_vec,1,chi,nchi32)
            call zgeru(nchi32,nchi32,sgn*weight,R_vec,1,P_vec,1,chi,nchi32)
         end do band2_loop
      end do band1_loop
   end do kloop
   !$omp end parallel do
   irr_chi_sc=chi(:,:)/Nk
  end function irr_chi_sc
end module calc_irr_chi_sc