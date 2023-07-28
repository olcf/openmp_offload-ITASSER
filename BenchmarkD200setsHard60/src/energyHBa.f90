      module energyHBa_mod

      contains
      function energyHBa(i,k,a,axki,ayki,azki,dxi,dyi,dzi,dxk,dyk,dzk)
      use params
      use ENERGY
      use hba
      use hbb
      implicit integer(i-z)
      integer, value :: i,k
      real, value :: a,axki,ayki,azki,dxi,dyi,dzi,dxk,dyk,dzk
!// !$acc routine seq
!$OMP declare target

! !$OMP declare target
!c     because helix is right-hand, bi=v(i-1)(x)v(i) is always
!c     same direction as v(i---k), if (k>i);
!c     reverse direction as v(i---k), if (k<i); 
      energyHBa=0

      bxk=5.7*dxk
      byk=5.7*dyk
      bzk=5.7*dzk
      bxi=5.7*dxi
      byi=5.7*dyi
      bzi=5.7*dzi
      if(k.gt.i)then            !bki,bxk,axki are in same directory
         br=(bxk-axki)**2+(byk-ayki)**2+(bzk-azki)**2
         if(br.lt.acut_hh)then
            b=1/(1+abs(br-3.2)) !br->3.4
            energyHBa=energyHBa-EHBIJ(i,k)*a*b
         endif
         br=(bxi-axki)**2+(byi-ayki)**2+(bzi-azki)**2
         if(br.lt.acut_hh)then
            b=1/(1+abs(br-3.2))
            energyHBa=energyHBa-EHBIJ(i,k)*a*b
         endif
      else                      !bki,bxk, and axki are in reverse directory
         br=(bxk+axki)**2+(byk+ayki)**2+(bzk+azki)**2
         if(br.lt.acut_hh)then
            b=1/(1+abs(br-3.2))
            energyHBa=energyHBa-EHBIJ(i,k)*a*b
         endif
         br=(bxi+axki)**2+(byi+ayki)**2+(bzi+azki)**2
         if(br.lt.acut_hh)then
            b=1/(1+abs(br-3.2))
            energyHBa=energyHBa-EHBIJ(i,k)*a*b
         endif
      endif 
! c^^^^^^^^^^^^ H-bond energy finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end function energyHBa
      end module energyHBa_mod
