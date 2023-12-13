      function zyrmsd(arms)
	use params
      use chain1
      implicit integer (i-z)
!      parameter(ndim=1500)	!maximum length of chain-length
      common/rmsdrange/nca1,nca2
      common/CA/dx(ndim),dy(ndim),dz(ndim)
      double precision r_m(3,ndim),r_d(3,ndim),w(ndim)
      double precision u(3,3),t(3),rms !armsd is real
      data w /ndim*1.0/

      k=0
      do i=nca1,nca2
         k=k+1
         r_m(1,k)=x(i+1)*0.87   !model conformation
         r_m(2,k)=y(i+1)*0.87
         r_m(3,k)=z(i+1)*0.87
         r_d(1,k)=dx(i)         !native conformation
         r_d(2,k)=dy(i)
         r_d(3,k)=dz(i)
      enddo

      nn=k                      !number of data points
      call u3b(w,r_m,r_d,nn,0,rms,u,t,ier)
      arms=dsqrt(rms/nn)	!RMSD is real, rms is double precision
      zyrmsd=arms

      return
      end
