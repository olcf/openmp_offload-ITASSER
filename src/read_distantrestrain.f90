      subroutine read_distantrestrain
	use params
      use lengths
      use distres
      use RCN
      implicit integer(i-z)
!      parameter(ndim=1500)
!      parameter(nvec=416)
!      COMMON/RCN/Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
      COMMON/RCN1/er1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
!      common/lengths/Lch,Lch1,Lch2
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
!      common/distres/er4,es3c
      common/resnumber/Ncom,Ndis,accur
      common/cutoff/cut1a,cut2a,cut3a,cut4a,cut1b,cut2b,cut3b,cut4b

      dimension r1(10000),r2(10000),dis(10000),deviation(10000)
      
      read(17,*)Ndis
      if(Ndis.gt.10000)then
         write(*,*)'warning: too many short range distance restraints!'
         Ndis=10000
      endif
      do i=1,Ndis
         read(17,*)r1(i),r2(i),nothing,dis(i),deviation(i)
         dis(i)=dis(i)/0.87
         deviation(i)=deviation(i)/0.87
         if(deviation(i).lt.0.5)deviation(i)=0.5
      enddo

      do i=1,Lch
         Mdis(i)=0              !number of prediction for 'i'
         do j=1,Ndis
            if(r1(j).eq.i)then
               Mdis(i)=Mdis(i)+1
               kdis(i,Mdis(i))=r2(j) !r2(j) with 'i'
               dist(i,Mdis(i))=dis(j) !predicted distance for i<->r2(j)
               dev(i,Mdis(i))=deviation(j)
            endif
            if(r2(j).eq.i)then
               Mdis(i)=Mdis(i)+1
               kdis(i,Mdis(i))=r1(j)
               dist(i,Mdis(i))=dis(j) !predicted distance for i<->r2(j)
               dev(i,Mdis(i))=deviation(j)
            endif
         enddo
      enddo
c      dilim=1.5*Ndis          !background number for derviation
      dilim=0.5*Ndis          !background number for derviation

!$acc update device(Mdis,kdis,dist,dev)

c      write(*,*)
c      write(*,*)'----------- distant map ---------------'
c      nnc=0
c      do i=1,Lch
c         nnc=nnc+Mdis(i)
c         write(*,12)i,Mdis(i),(Kdis(i,j),dist(i,j)*0.87,j=1,Mdis(i))
c 12      format(i4,':',i3,20(i4,'-'f5.2))
c      enddo
c      write(*,*)'Number of distmap=',nnc,' Lch=',Lch

c^^^^^^^^^^^^^^^^^^ read contact restrains finished ^^^^^^^^^^^^^^^^^
      return
      end
