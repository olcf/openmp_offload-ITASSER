      subroutine set_temperature
	use params
      implicit integer(i-z)
!      parameter(nrep=100)
      common/commonuse1/random,ncycle
      common/commonuse2/atemp1,atemp2,N_rep,phot
      common/sw1/aT_rep(nrep),E_rep(nrep)
      common/initialinput/switch,k_cycle,k_phot,N_ann
      common/thrtem/aT_ann(100)
      common/aTs/aTs1,aTs2,aTs_rep(nrep),aTs
      common/aTTs/aTTs1,aTTs2,aTTs_rep(nrep),aTTs

***********for normal run ***********************************
      do i=1,N_REP
         aT_rep(i)=atemp1*(atemp2/atemp1)**(float(i-1)/(N_rep-1))
         aTs_rep(i)=aTs1*(aTs2/aTs1)**(float(i-1)/(N_rep-1))
         aTTs_rep(i)=aTTs1*(aTTs2/aTTs1)**(float(i-1)/(N_rep-1))
c     write(*,*)i,aT_rep(i)
      enddo

c      stop
c^^^^^^^^^^^^^^^^^ set aT_rep(i) finished ^^^^^^^^^^^^^^^
      return
      end
