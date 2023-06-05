      subroutine swap(i1,i2)
      use params
      !use openacc
      use lengths
      use trackn
      implicit integer (i-z)
!      parameter(nrep=100)
!      parameter(ndim=1500)	!maximum length of chain-length
!      parameter(nvec=416)
!      common/lengths/Lch,Lch1,Lch2
      common/nswap/bNSa(100),bNSt(100)

      common/sw1/aT_rep(nrep),E_rep(nrep)
      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      common/sw4/exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep)
      common/sw5/egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep)
      common/chain0/ras(ndim),nfl
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      
!      common/trackn/n_tem(100)

      if(metro_swap(i1,i2).eq.1)then
         bNSa(i1)=bNSa(i1)+1    !number of accepted swaps
! !$acc data copy(tempx,tempy,tempz,xrep,yrep,zrep)
!$acc data copyin(xrep(:lch,:),zrep(:lch,:),yrep(:lch,:)) !copyout()
!$acc kernels
         do i=1,Lch
*     x:
            tempx=xrep(i,i1)
            tempy=yrep(i,i1)
            tempz=zrep(i,i1)
            xrep(i,i1)=xrep(i,i2)
            yrep(i,i1)=yrep(i,i2)
            zrep(i,i1)=zrep(i,i2)
            xrep(i,i2)=tempx
            yrep(i,i2)=tempy
            zrep(i,i2)=tempz
         enddo
!$acc end kernels
***   swap E ----------------------------->
         attt=E_rep(i1)         !exchange E_rep
         E_rep(i1)=E_rep(i2)
         E_rep(i2)=attt
***   swap n_tem ------------------------->
         temp=n_tem(i1)
         n_tem(i1)=n_tem(i2)
         n_tem(i2)=temp
***   swap moveable points---------------->
!$acc end data
      endif
      bNSt(i1)=bNSt(i1)+1       !number of total swaps.

      return
      end
