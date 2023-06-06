      subroutine swap_RS(i1,i2)
	use params
      !use openacc
      use chainm
      use lengths
      implicit integer (i-z)
!      parameter(nrep=100)
!      parameter(ndim=1500)	!maximum length of chain-length
!      parameter(nvec=416)
!      common/lengths/Lch,Lch1,Lch2
      common/nswap/bNSa(100),bNSt(100)

      common/sw1/aT_rep(nrep),E_rep(nrep)
      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      common/sw3/icarep(ndim,nrep)
      common/sw4/exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep)
      common/sw5/egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep)
      common/sw7/ecxrep(ndim,nrep),ecyrep(ndim,nrep),eczrep(ndim,nrep)
      common/sw8/ebxrep(ndim,nrep),ebyrep(ndim,nrep),ebzrep(ndim,nrep)
      common/sw9/etxrep(ndim,nrep),etyrep(ndim,nrep),etzrep(ndim,nrep)
      common/chain0/ras(ndim),nfl
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
!      common/chainm/mv(ndim)

!$acc data copy(icap_rep,xrep,yrep,zrep,ezrep,icarep,
!$acc&           etyrep,etxrep,egzrep,ebyrep,ebxrep,ecxrep,
!$acc&           ecyrep, eczrep,exrep,eyrep,exrep,egxrep,
!$acc&           etzrep,egyrep,abzrep) present (mv)
      if(metro_swap(i1,i2).eq.1)then
         bNSa(i1)=bNSa(i1)+1    !number of accepted swaps
! !$acc parallel loop num_gangs(100) vector_length(128)
!$acc parallel loop
         do i=1,Lch
            t_ica=icarep(i,i1)             !ica(i)
            icarep(i,i1)=icarep(i,i2)
            icarep(i,i2)=t_ica

            if(mv(i).gt.0)then
               t_x=xrep(i,i1)              !x(i)
               t_y=yrep(i,i1)
               t_z=zrep(i,i1)
               xrep(i,i1)=xrep(i,i2)
               yrep(i,i1)=yrep(i,i2)
               zrep(i,i1)=zrep(i,i2)
               xrep(i,i2)=t_x
               yrep(i,i2)=t_y
               zrep(i,i2)=t_z
  
          else
               e_x=exrep(i,i1)           !ex(i)
               e_y=eyrep(i,i1)
               e_z=ezrep(i,i1)
               exrep(i,i1)=exrep(i,i2)
               eyrep(i,i1)=eyrep(i,i2)
               ezrep(i,i1)=ezrep(i,i2)
               exrep(i,i2)=e_x
               eyrep(i,i2)=e_y
               ezrep(i,i2)=e_z
               e_x=egxrep(i,i1)          !egx(i)
               e_y=egyrep(i,i1)
               e_z=egzrep(i,i1)
               egxrep(i,i1)=egxrep(i,i2)
               egyrep(i,i1)=egyrep(i,i2)
               egzrep(i,i1)=egzrep(i,i2)
               egxrep(i,i2)=e_x
               egyrep(i,i2)=e_y
               egzrep(i,i2)=e_z
               e_x=ecxrep(i,i1)          !ecx(i)
               e_y=ecyrep(i,i1)
               e_z=eczrep(i,i1)
               ecxrep(i,i1)=ecxrep(i,i2)
               ecyrep(i,i1)=ecyrep(i,i2)
               eczrep(i,i1)=eczrep(i,i2)
               ecxrep(i,i2)=e_x
               ecyrep(i,i2)=e_y
               eczrep(i,i2)=e_z
               e_x=ebxrep(i,i1)          !ebx(i)
               e_y=ebyrep(i,i1)
               e_z=ebzrep(i,i1)
               ebxrep(i,i1)=ebxrep(i,i2)
               ebyrep(i,i1)=ebyrep(i,i2)
               ebzrep(i,i1)=ebzrep(i,i2)
               ebxrep(i,i2)=e_x
               ebyrep(i,i2)=e_y
               ebzrep(i,i2)=e_z
               e_x=etxrep(i,i1)          !etx(i)
               e_y=etyrep(i,i1)
               e_z=etzrep(i,i1)
               etxrep(i,i1)=etxrep(i,i2)
               etyrep(i,i1)=etyrep(i,i2)
               etzrep(i,i1)=etzrep(i,i2)
               etxrep(i,i2)=e_x
               etyrep(i,i2)=e_y
               etzrep(i,i2)=e_z
            endif
         enddo
! !$acc end kernels
***   swap E ----------------------------->
         attt=E_rep(i1)         !exchange E_rep
         E_rep(i1)=E_rep(i2)
         E_rep(i2)=attt
***   swap moveable points---------------->
      endif
      bNSt(i1)=bNSt(i1)+1       !number of total swaps.

c^^^^^^^^^^^^^^^^^ swap_RS finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!$acc end data     
      return
      end
