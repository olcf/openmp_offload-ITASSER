      subroutine set_current_RS
	use params
      use chainm
      use chain1
      use echain1
      use lengths
      use echain2
      use echain4
      use echain5
      use echain6
      implicit integer(i-z)
!      parameter(ndim=1500)
!      parameter(nrep=100)
!      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
!      common/lengths/Lch,Lch1,Lch2
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)

      common/sw1/aT_rep(nrep),E_rep(nrep)
      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      common/sw3/icarep(ndim,nrep)
      common/sw4/exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep) !Ca
      common/sw5/egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep) !SG
      common/sw7/ecxrep(ndim,nrep),ecyrep(ndim,nrep),eczrep(ndim,nrep) !cc
      common/sw8/ebxrep(ndim,nrep),ebyrep(ndim,nrep),ebzrep(ndim,nrep) !hb
      common/sw9/etxrep(ndim,nrep),etyrep(ndim,nrep),etzrep(ndim,nrep) !CB

      common/chain0/ras(ndim),nfl
!      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
!      common/echain1/ex(ndim),ey(ndim),ez(ndim) !Ca
!      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
!      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
!      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
!      common/echain6/etx(ndim),ety(ndim),etz(ndim) !CB
!      common/chainm/mv(ndim)

      do i=1,Lch
         x(i)=xrep(i,itemp)
         y(i)=yrep(i,itemp)
         z(i)=zrep(i,itemp)
         ica(i)=icarep(i,itemp)
         if(mv(i).lt.0)then
            ex(i)=exrep(i,itemp) !Ca
            ey(i)=eyrep(i,itemp)
            ez(i)=ezrep(i,itemp)
            egx(i)=egxrep(i,itemp) !SG
            egy(i)=egyrep(i,itemp)
            egz(i)=egzrep(i,itemp)
            ecx(i)=ecxrep(i,itemp) !cc
            ecy(i)=ecyrep(i,itemp)
            ecz(i)=eczrep(i,itemp)
            ebx(i)=ebxrep(i,itemp) !hb
            eby(i)=ebyrep(i,itemp)
            ebz(i)=ebzrep(i,itemp)
            etx(i)=etxrep(i,itemp) !CB
            ety(i)=etyrep(i,itemp)
            etz(i)=etzrep(i,itemp)
         endif
      enddo
      ica(0)=ica(2)             !only useful for moveable pieces
      ica(Lch)=ica(Lch2)        !only useful for moveable pieces
!$acc update device(ex,ey,ez,egx,egy,egz,ecx,ecy,ecz,
!$acc&       ebx,eby,ebz,etx,ety,etz,x,y,z,ica)
      return
      end
