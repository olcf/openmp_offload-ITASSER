      subroutine beta_sheet(x,y,z,n)
	use params
!      parameter(ndim=1500)
      dimension x(ndim),y(ndim),z(ndim)

      dist2=6.70098             !distance between i and i+2
      high=1.792820             !distance between i+1 and middle of (i,i+2)
      do i=1,n
         x(i)=dist2/2*(i-1)
         if(int(i/2)*2.eq.i)then
            y(i)=0
         else
            y(i)=high
         endif
         z(i)=0
      enddo

      return
      end
