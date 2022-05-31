      function mcheck(i0,bx0,by0,bz0)
	use params
      use lengths
!      parameter(ndim=1500)
      implicit integer(i-z)
!      parameter(nvec=416)
!      common/lengths/Lch,Lch1,Lch2
      common/initialstr/cx(ndim),cy(ndim),cz(ndim)
      common/mcheck_dis/amcheck_dis

      mcheck=1
      do i=1,Lch
         if(i0.ne.i)then
            dis=di(bx0,by0,bz0,cx(i),cy(i),cz(i)) !distance
            if(dis.le.amcheck_dis)then
               mcheck=3
               return
            endif
         endif
      enddo

      return
      end
