      module ei5_mod
   
      contains

      function ei5(i,idist)
      use params
      use lengths
      use short1
      use chainm
      use chain1
      use echain1
      implicit none
      integer, value :: i, idist
      integer :: im2, ip2
      real :: axm, aym, azm
      real :: axp, ayp, azp
      real :: ei5
      integer :: xxx
!$acc routine seq
!$OMP declare target
      ei5=0
      if(idist.gt.5)then
         im2=i-2
         ip2=i+2
         if(im2.ge.1.and.ip2.le.Lch)then
            if(mv(im2).gt.0)then
               axm=x(im2)
               aym=y(im2)
               azm=z(im2)
            else
               axm=ex(im2)
               aym=ey(im2)
               azm=ez(im2)
            endif
            if(mv(ip2).gt.0)then
               axp=x(ip2)
               ayp=y(ip2)
               azp=z(ip2)
            else
               axp=ex(ip2)
               ayp=ey(ip2)
               azp=ez(ip2)
            endif
            xxx=nint((axm-axp)**2+(aym-ayp)**2+(azm-azp)**2)
            if(xxx.gt.500)xxx=500
            ei5=acops(im2,jbin(xxx)) !i2-residue, jbin. es<0.
         endif
      endif
      return
      end function ei5
      end module ei5_mod
