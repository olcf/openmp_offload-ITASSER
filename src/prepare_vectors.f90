      subroutine prepare_vectors
	use params
      implicit integer(i-z)
!      parameter(nvec=416)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      dimension n(100)
      common/lattice/m_latt,latt1,latt2

      do i=1,100
         n(i)=0
      enddo

      nwmax=0
      aaa=0
      nn=5
      do x=-nn,nn
         do y=-nn,nn
            do z=-nn,nn
               vector(x,y,z)=0
               r=x*x+y*y+z*z
               if(r.ge.latt1.and.r.le.latt2) then
                  nwmax=nwmax+1
                  vx(nwmax)=x
                  vy(nwmax)=y
                  vz(nwmax)=z
                  vector(x,y,z)=nwmax
c                  write(*,*)nwmax,vx(nwmax),vy(nwmax),vz(nwmax),r,
c     $                 sqrt(float(r)*0.87*0.87)
                  n(r)=n(r)+1
                  aaa=aaa+sqrt(float(r)*0.87*0.87)
               endif
            enddo
         enddo
      enddo
      aaa=aaa/float(nwmax)
      write(20,*)'n1_all=',nwmax,'  <vr>=',aaa

c      do i=1,nwmax
c         write(*,*)i,vx(i),vy(i),vz(i)
c      enddo
c      stop

c      do i=10,30
c         write(*,*)i,n(i),sqrt(float(i)*0.87*0.87)
c      enddo
c      stop
      
c     i=1,5
c           10          24   2.751182    
c           11          24   2.885463    
c           12           8   3.013768    +
c           13          24   3.136830    +
c           14          48   3.255242    x
c           15           0   3.369496    
c           16           6   3.480000    x
c           17          48   3.587102    x
c           18          36   3.691097    x
c           19          24   3.792242    x
c           20          24   3.890758    x
c           21          48   3.986841    x
c           22          24   4.080662    x
c           23           0   4.172373    
c           24          24   4.262112    x
c           25          30   4.350000    x
c           26          72   4.436147    +
c           27          32   4.520653    
c           28           0   4.603607    
c           29          72   4.685093    
c           30          48   4.765186    
c nwmax=         616  <vr>=   3.982909    
c nwmax=         312  <vr>=   3.809868    

c      stop
      return
      end
