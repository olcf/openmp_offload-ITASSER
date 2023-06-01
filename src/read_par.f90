      subroutine read_par
	use params
      use lengths
      use pair1
      use par
      implicit integer(i-z)
!      parameter(ndim=1500)
!      common/lengths/Lch,Lch1,Lch2
!      common/par/apar(ndim,ndim)
      common/fractpair/fract_pair1,fract_pair3
      common/pair33/i_pair3,mk_pair3
      common/weight/chuan
!      common/pair1/eh2,eh1b,eh1c

      dimension apar1(ndim,ndim) !from par.dat
      dimension apar2(ndim,ndim) !from pair1.dat

c     read from 'par.dat'-------------->
      Nline=1000
      do i=1,Lch
         read(15,*)
         do i_line=1,Nline
            line_end=min(10*i_line,Lch) !25,50,75,100,...., ending point
            read(15,*)(apar1(i,j),j=(i_line-1)*10+1,line_end)
            if(line_end.ge.Lch) go to 1
         enddo
 1       continue
      enddo

c     read from 'pair1.dat'-------------->
      if(i_pair3.gt.0)then      ! 'pair1.dat' exist
         Nline=100              !maximum Lch is 2500 !!!
         do i=1,Lch
            read(26,*)
            do i_line=1,Nline
               line_end=min(25*i_line,Lch) !25,50,75,100,...., ending point
               read(26,*)(apar2(i,j),j=(i_line-1)*25+1,line_end)
               if(line_end.eq.Lch) go to 2
            enddo
 2          continue
            read(26,*)
         enddo
      endif
      
***   combine and rescale the pair interaction ---------------->
      apar1_max=-100
      apar2_max=-100
      do i=1,Lch
         do j=1,Lch
            if(abs(apar1(i,j)).gt.apar1_max)apar1_max=abs(apar1(i,j))
            if(abs(apar2(i,j)).gt.apar2_max)apar2_max=abs(apar2(i,j))
         enddo
      enddo
      do i=1,Lch
         do j=1,Lch
            apar(i,j)=0
            if(i_pair3.gt.0)then !'pair3.dat' and 'pair1.dat' exist
               if(apar1_max.lt.0.00001)then !par.dat is wrong
                  apar(i,j)=apar2(i,j) !'pair1.dat' only
               else
                  apar(i,j)=apar1(i,j)*apar2_max/apar1_max*
     &                 (1-fract_pair1)+apar2(i,j)*fract_pair1
               endif
            else
               apar(i,j)=apar1(i,j) !'pair1.dat' not exist, 'par.dat' only
            endif
         enddo
      enddo
      
c      write(*,*)apar1_max,apar2_max
c      do i=1,Lch
c         do j=1,Lch
c            write(*,*)i,j,apar1(i,j),apar2(i,j),apar(i,j)
c         enddo
c      enddo
c      stop

c^^^^^^^^^^^^^^^^^ pair-wise potential finished ^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
