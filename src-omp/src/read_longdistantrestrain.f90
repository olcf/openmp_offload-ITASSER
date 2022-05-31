      subroutine read_longdistantrestrain
	use params
      use lengths
      use longdist
      implicit integer(i-z)
!!      parameter(ndim=1500)
!!      parameter(nvec=416)
!      common/lengths/Lch,Lch1,Lch2
!      COMMON/longdist/MdisL(ndim),kdisL(ndim,500),distL(ndim,500)
c     maximum N=4*L/10, if L=100, maximum N=400
      dimension r1(400000),r2(400000),dis(400000)

      read(22,*)NdisL
      if(NdisL.gt.400000)NdisL=400000
      do i=1,NdisL
         read(22,*)r1(i),r2(i),dis(i)
         dis(i)=dis(i)/0.87
      enddo
      
      do i=1,Lch
         MdisL(i)=0             !number of prediction for 'i'
         do j=1,NdisL
            if(r1(j).eq.i)then
               MdisL(i)=MdisL(i)+1
               kdisL(i,MdisL(i))=r2(j) !r2(j) with 'i'
               distL(i,MdisL(i))=dis(j) !predicted distance for i<->r2(j)
               if(MdisL(i).ge.500)goto 101
            endif
            if(r2(j).eq.i)then
               MdisL(i)=MdisL(i)+1
               kdisL(i,MdisL(i))=r1(j)
               distL(i,MdisL(i))=dis(j) !predicted distance for i<->r2(j)
               if(MdisL(i).ge.500)goto 101
            endif
         enddo
 101     continue
         if(MdisL(i).gt.500)then !size>5000/4
            write(*,*)i,MdisL(i),'N_res>500 for one residue, exit'
            stop
         endif
      enddo
!$acc update device(MdisL,kdisL,distL)

c      write(*,*)
c      write(*,*)'----------- distant map ---------------'
c      nnc=0
c      do i=1,Lch
c        nnc=nnc+MdisL(i)
c        write(*,12)i,MdisL(i),(KdisL(i,j),distL(i,j)*0.87,j=1,MdisL(i))
c 12     format(i4,':',i3,20(i4,'-'f5.2))
c      enddo
c      write(*,*)'Number of distmap=',nnc,' Lch=',Lch
c      stop
      
c^^^^^^^^^^^^^^^^^^ read contact restrains finished ^^^^^^^^^^^^^^^^^
      return
      end
