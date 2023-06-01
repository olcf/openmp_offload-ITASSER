      subroutine test_overlap
	use params
      use lengths
      use aax_mod
      use aay_mod
      use aaz_mod
      implicit integer(i-z)
!      common/lengths/Lch,Lch1,Lch2
      common/looks/exc,exc1,exc2

      do i=1,Lch
         do j=i+3,Lch
            dis2=(aax(i)-aax(j))**2+(aay(i)-aay(j))**2+
     $           (aaz(i)-aaz(j))**2
            if(dis2.lt.exc)then
               write(*,*)i,j,dis2,'  Ca-Ca overlap'
            endif
         enddo
      enddo

      return
      end
