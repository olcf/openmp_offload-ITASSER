      module energyHBb_mod
  
      contains

      function energyHBb(i,k,a,axki,ayki,azki,dxi,dyi,dzi,dxk,dyk,dzk)
	use params
      use ENERGY
      use hba
      use hbb
      implicit none
!$acc routine seq
!$OMP declare target 
!      implicit none
      integer, value :: i, k
      real, value :: a,axki,ayki,azki,dxi,dyi,dzi,dxk,dyk,dzk
      real :: bxk,byk,bzk,bxi,byi,bzi,br,b
      real :: energyHBb


      energyHBb=0

      bxk=5.3*dxk
      byk=5.3*dyk
      bzk=5.3*dzk
      bxi=5.3*dxi
      byi=5.3*dyi
      bzi=5.3*dzi
      br=(bxk-axki)**2+(byk-ayki)**2+(bzk-azki)**2
      if(br.lt.bcut_hh) then
         b=2./(2.+br)           !br->0
         energyHBb=energyHBb-EHBIJ(i,k)*a*b
      endif
      br=(bxk+axki)**2+(byk+ayki)**2+(bzk+azki)**2
      if(br.lt.bcut_hh) then
         b=2./(2.+br)
         energyHBb=energyHBb-EHBIJ(i,k)*a*b
      endif
      br=(bxi-axki)**2+(byi-ayki)**2+(bzi-azki)**2
      if(br.lt.bcut_hh) then
         b=2./(2.+br)
         energyHBb=energyHBb-EHBIJ(i,k)*a*b
      endif
      br=(bxi+axki)**2+(byi+ayki)**2+(bzi+azki)**2
      if(br.lt.bcut_hh) then
         b=2./(2.+br)
         energyHBb=energyHBb-EHBIJ(i,k)*a*b
      endif

c^^^^^^^^^^^^ H-bond energy finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end function energyHBb
      end module energyHBb_mod
