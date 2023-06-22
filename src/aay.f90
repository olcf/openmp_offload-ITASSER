        module aay_mod

        contains

        function aay(i)
        use params
        use chainm
        use chain1
        use echain1
!$acc routine seq
        implicit none
        integer, value :: i
        real :: aay
!$OMP declare target 

      if(mv(i).gt.0)then
         aay=y(i)
      else
         aay=ey(i)
      endif
      return
      end function aay
      end module aay_mod
