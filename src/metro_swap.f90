      integer function metro_swap(i,j)
	use params
!      parameter(nrep=100)
      common/sw1/aT_rep(nrep),E_rep(nrep)
      common/ranzy/nozy

! !$acc data copy(aT_rep,E_rep,aweight)

c     w_12=wi(j)wj(i)/wi(i)wj(j) ------------->
      aaa=(1/aT_rep(i)-1/aT_rep(j))*(E_rep(i)-E_rep(j))
      aweight=exp(aaa)

      if(aranzy(nozy).le.aweight)then
! !$acc kernels
         metro_swap=1		!swap accepted
! !$acc end  kernels
      else
! !$acc kernels
         metro_swap=3		!swap rejected
! !$acc end kernels
      endif
! !$acc end data
      return
      end
