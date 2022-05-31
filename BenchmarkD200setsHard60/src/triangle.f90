      subroutine triangle(bx,by,bz,br,bthita,bphi)
	use params
      parameter(bpi=3.1415926)

      br=sqrt(bx*bx+by*by+bz*bz)
      bthita=acos(bz/br)
      if(abs(bx).gt.0.00001)then
         bphi=atan(abs(by/bx))
      else
         bphi=0
      endif
      if(bx.gt.0)then
         if(by.lt.0)bphi=bphi+bpi*1.5
      else
         if(by.gt.0)then
            bphi=bphi+bpi*0.5
         else
            bphi=bphi+bpi
         endif
      endif

c^^^^^^^^^^^^^^^^^^ triangle done ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
