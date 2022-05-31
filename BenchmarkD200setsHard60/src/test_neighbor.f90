      subroutine test_neighbor
	use params
      use chain1
      use lengths
      implicit integer(i-z)
!      parameter(ndim=1500)
!      parameter(nrep=100)
!      parameter(nvec=416)
!      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
!      common/lengths/Lch,Lch1,Lch2
      common/logica/goodc
      logical goodc(nvec,nvec)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)

      do j=1,Lch
         i=j-1
         if(j.ne.1.and.j.ne.Lch) then
            ii=ica(i)
            jj=ica(j)
            if(.not.goodc(ii,jj)) then
           write(20,8112)i,j,vx(ii),vy(ii),vz(ii),vx(jj),vy(jj),vz(jj)
 8112          format(5x,'warning2 -wrong input chain - vectors ',8i4)
               stop
            endif
         endif
      enddo 

      return
      end
