      subroutine set_move_retio !set movement percentage
	use params
      implicit integer(i-z)
      common/commonuse2/atemp1,atemp2,N_rep,phot

      common/moveretio1/h2,h3s,h3d,h4s,h4d,h5s,h5d,h6,h7,hend
      common/moveretio2/bh2,bh3s,bh3d,bh4s,bh4d,bh5s,bh5d,bh6
      common/moveretio3/bh7a,bh7b,bhendn,bhendc

      hsum=h2+h3s+h3d+h4s+h4d+h5s+h5d+h6+hend
      bh2=h2/hsum
      bh3s=bh2+h3s/hsum
      bh3d=bh3s+h3d/hsum
      bh4s=bh3d+h4s/hsum
      bh4d=bh4s+h4d/hsum
      bh5s=bh4d+h5s/hsum
      bh5d=bh5s+h5d/hsum
      bh6=bh5d+h6/hsum
      bhendn=bh6+hend/2./hsum
      bhendc=bhendn+hend/2./hsum

c      write(*,1)h2,h3s,h3d,h4s,h4d,h5s,h5d,h6,hend
c      write(*,1)bh2,bh3s,bh3d,bh4s,bh4d,bh5s,bh5d,bh6,
c     $     bhendn,bhendc
c      write(*,1)bh2,bh3s-bh2,bh3d-bh3s,bh4s-bh3d,bh4d-bh4s,bh5s-bh4d,
c     $     bh5d-bh5s,bh6-bh5d,bhendn-bh6,bhendc-bhendn
c 1    format(13f6.3)
c      stop

c ^^^^^^^^^^^^^^^^^^^^^^^ set retio finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
