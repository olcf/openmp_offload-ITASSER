	module aay_mod

      contains

      function aay(i)
	use params
!$acc routine seq
! !$OMP declare target 
      use chainm
      use chain1
      use echain1
!$OMP declare target      
      implicit none
      integer, value :: i
      real :: aay

      if(mv(i).gt.0)then
         aay=y(i)
      else
         aay=ey(i)
      endif
      return
      end function aay
      end module aay_mod
