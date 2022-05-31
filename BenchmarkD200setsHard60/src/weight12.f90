      function weight12(dE)
	use params
!      parameter(nrep=100)
      common/temperature/itemp,atemp
      common/aTs/aTs1,aTs2,aTs_rep(nrep),aTs
      common/aTTs/aTTs1,aTTs2,aTTs_rep(nrep),aTTs
      common/ichos/ichos

      if(ichos.eq.1)then
         weight12=exp(-dE/atemp)
         ichos=2
      else
         weight12=exp(-square2(dE)/aTs)
         ichos=1
      endif
c     write(*,*)ichos,aTs,aTTs,dE,weight12

      return
      end
