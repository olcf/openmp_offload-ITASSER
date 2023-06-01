      subroutine read_exp
	use params
      use lengths
      use expose
      implicit integer(i-z)
!      parameter(ndim=1500)
!      common/lengths/Lch,Lch1,Lch2
!      common/expose/mp(20,ndim),area(ndim)
      common/nana1/nana

      na=nana
      read(27,*)
      do i=1,Lch
         read(27,*)itmp,(mp(j,i),j=1,na)
         area(i)=0
         do j=1,na
            if(mp(j,i).eq.0)mp(j,i)=-1 !-1, bury; 1, expose
            area(i)=area(i)+mp(j,i)
         enddo
c         write(*,*)i,mp(1,i),mp(2,i),mp(12,i)
      enddo
      close(27)
c      stop
!$acc update device(area)

c^^^^^^^^^^^^^^^^^^ read contact restrains finished ^^^^^^^^^^^^^^^^^
      return
      end
