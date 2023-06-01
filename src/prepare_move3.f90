      subroutine prepare_move3
	use params
      implicit integer(i-z)
!      parameter(nvec=416)
      common/logica/goodc
      logical goodc(nvec,nvec)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/move31/v31(-15:15,-15:15,-15:15,6)
      common/move32/v32(-15:15,-15:15,-15:15,6)
      common/move33/v33(-15:15,-15:15,-15:15,6)
      common/move34/Np3(-15:15,-15:15,-15:15)
      common/nwmax1/nvecr

c      dimension v31(-15:15,-15:15,-15:15,100)
c      dimension v32(-15:15,-15:15,-15:15,100)
c      dimension v33(-15:15,-15:15,-15:15,100)
c      dimension Np3(-15:15,-15:15,-15:15)

      do i=-15,15
         do j=-15,15
            do k=-15,15
               Np3(i,j,k)=0
            enddo
         enddo
      enddo


      num3=0
      max=0
      do 101 i=1,nvecr
         do 102 j=1,nvecr
            if(goodc(i,j))then
               do 103 k=1,nvecr
                  if(goodc(j,k))then
                     num3=num3+1
                     rx=vx(i)+vx(j)+vx(k)
                     ry=vy(i)+vy(j)+vy(k)
                     rz=vz(i)+vz(j)+vz(k)
                     Np3(rx,ry,rz)=Np3(rx,ry,rz)+1
                     v31(rx,ry,rz,Np3(rx,ry,rz))=i
                     v32(rx,ry,rz,Np3(rx,ry,rz))=j
                     v33(rx,ry,rz,Np3(rx,ry,rz))=k
                     if(max.le.Np3(rx,ry,rz)) max=Np3(rx,ry,rz)
                     write(*,*)i,j,k,max
                  endif
 103           continue
            endif
 102     continue
 101  continue

      write(*,*)'number of move3=',num3,max

      stop
      return
      end
