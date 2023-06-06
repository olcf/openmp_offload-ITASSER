      subroutine alpha_helix(x,y,z,n)
	use params
!      parameter(ndim=1500)
      dimension x(ndim),y(ndim),z(ndim)

      rad=2.3                   !redius of helix
      do i=1,n
         angle=100*3.1415926/180*(i-1) !100 degree per residues
         x(i)=rad*cos(angle)
         y(i)=rad*sin(angle)
         z(i)=1.5*(i-1)            !increase 1.5 per residue
      enddo

      return
      end
