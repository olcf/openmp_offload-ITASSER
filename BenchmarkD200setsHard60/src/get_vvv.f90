      subroutine get_vvv
	use params
      use chainm
      use lengths
      use aax_mod
      use aay_mod
      use aaz_mod
      implicit integer (i-z)
!      parameter(ndim=1500)      !maximum length of chain-length
      common/excluded/vvv(ndim,ndim)
!      common/chainm/mv(ndim)
!      common/lengths/Lch,Lch1,Lch2
      common/looks/exc,exc1,exc2

      do i=1,Lch
         do j=1,Lch
            vvv(i,j)=1          !every pair should be checked
            if(mv(i).lt.0.and.mv(j).lt.0)then
               if(mv(i).eq.mv(j))then !belong to the same fragment
                  vvv(i,j)=-1   !Not be checked
               else
                  dist2=(aax(i)-aax(j))**2+(aay(i)-aay(j))**2
     $                 +(aaz(i)-aaz(j))**2
                  if(dist2.lt.exc)then !initially violated
                     vvv(i,j)=-2
                  endif
               endif
            endif
c            write(*,*)i,j,vvv(i,j)
         enddo
      enddo

c^^^^^^^^^^^^^^^^^^^ vvv(i,j) finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      return
      end
