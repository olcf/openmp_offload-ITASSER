      module aaz_mod

      contains

      function aaz(i)
	use params
      use chainm
      use echain1
      use chain1
!$acc routine seq
!$OMP declare target 
      implicit none
      integer, value :: i
      real :: aaz
      if(mv(i).gt.0)then
         aaz=z(i)
      else
         aaz=ez(i)
      endif
      return
      end function aaz
      end module aaz_mod
