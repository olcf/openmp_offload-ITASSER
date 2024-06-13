      module aax_mod

      contains

      function aax(i)
      use params
      use chainm
      use chain1
      use echain1
! !$acc routine seq
      implicit none
!$OMP declare target      
      integer, value :: i
      real :: aax

      if(mv(i).gt.0)then
         aax=x(i)
      else
         aax=ex(i)
      endif
      return
! !$OMP end declare target 
      end function aax

      end module aax_mod
