      subroutine set_current
	use params
      use chain1
      use echain1
      use chainm 
      use lengths
      use echain2
      use echain4
      use echain5
      implicit integer(i-z)
!      parameter(ndim=1500)
!      parameter(nrep=100)
!      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
!      common/lengths/Lch,Lch1,Lch2
      common/sw1/aT_rep(nrep),E_rep(nrep)
      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      common/sw3/icarep(ndim,nrep)
      common/sw4/exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep)
      common/sw5/egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep)
      common/sw7/ecxrep(ndim,nrep),ecyrep(ndim,nrep),eczrep(ndim,nrep)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)

      common/chain0/ras(ndim),nfl
!      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
!      common/echain1/ex(ndim),ey(ndim),ez(ndim)
!      common/echain2/egx(ndim),egy(ndim),egz(ndim)
!      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim)
!      common/echain5/ebx(ndim),eby(ndim),ebz(ndim)

*     put the itemp'th replica conformation as current conformation------>
*     moveable points:
      do i=1,Lch
         x(i)=xrep(i,itemp)     !initial coordinate of itemp's replica
         y(i)=yrep(i,itemp)
         z(i)=zrep(i,itemp)
      enddo

      do i=1,Lch1
         j=i+1
         wx=x(j)-x(i)
         wy=y(j)-y(i)
         wz=z(j)-z(i)
         ica(i)=vector(wx,wy,wz) !identify order number of each bond-vector
      enddo
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)

c^^^^^^^^^^^ set current (x,y,z,ica) finished ^^^^^^^^^^^^^^^^^^^
      return
      end
