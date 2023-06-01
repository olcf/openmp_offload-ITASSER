      subroutine prepare_neighbors
	use params
      use lengths
      implicit integer(i-z)
!                parameter(nvec=416)
      common/logica/goodc
      logical goodc(nvec,nvec)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/three/angle(nvec,nvec)
!      common/lengths/Lch,Lch1,Lch2

      common/shift/u21(nvec,nvec),m12(nvec),u1(nvec,100),u2(nvec,100)
      common/nwmax1/nvecr

      max_m12=0
      do i=1,nvecr
         m12(i)=0
      enddo

      mmm=0
      nnn=0
      kkk=0
      do i=1,nvecr
      do j=1,nvecr
         u21(i,j)=0
         a2=vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)
         b2=vx(j)*vx(j)+vy(j)*vy(j)+vz(j)*vz(j)
         c2=(vx(i)+vx(j))**2+(vy(i)+vy(j))**2+(vz(i)+vz(j))**2
         cosangle=(a2+b2-c2)/(2*sqrt(a2*b2))
         angle(i,j)=acos(cosangle)*180/3.1415926
c     in database, angle is in [65,165];
         if(angle(i,j).gt.65.and.angle(i,j).lt.165)then
            goodc(i,j)=.true.
            mmm=mmm+1
            ijx=vx(i)+vx(j)
            ijy=vy(i)+vy(j)
            ijz=vz(i)+vz(j)
            do k=1,nvecr
               if(vx(k).eq.ijx.and.vy(k).eq.ijy.and.vz(k).eq.ijz)then
                  kkk=kkk+1
                  u21(i,j)=k    !vi+vj=vk
                  m12(k)=m12(k)+1
                  u1(k,m12(k))=i
                  u2(k,m12(k))=j
                  if(max_m12.lt.m12(k))max_m12=m12(k)
                  goto 10
               endif
            enddo
 10         continue
         else
            goodc(i,j)=.false.
         endif
         nnn=nnn+1
c         write(*,*)i,j,mmm,nnn,angle(i,j),goodc(i,j)
      enddo
      enddo

      n=0
      do i=1,nvecr
         r=vx(i)**2+vy(i)**2+vz(i)**2
         if(m12(i).gt.0)n=n+1
c     if(r.gt.17)write(*,*)i,r,m12(i)
      enddo
c      write(*,*)'n2_all=',nnn
c      write(*,*)'n2good=',mmm,'  n21=',kkk
c      write(*,*)'n1_all=',nvecr,'  n12=',n

c      stop
      return
      end
