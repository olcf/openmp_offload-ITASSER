      subroutine move_point
	use params
      use chainm
      use lengths
      implicit integer(i-z)
!      parameter(ndim=1500)
!      parameter(nrep=10)
!      parameter(nvec=416)
!      common/lengths/Lch,Lch1,Lch2
      common/chain0/ras(ndim),nfl
      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
!      common/chainm/mv(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      
      dimension mv1(5000)
      
      do i=1,Lch+10
         mv1(i)=-1              !frozen point
      enddo
      
      do i=1,nfl
         mv1(ras(i))=1          !moveable point
      enddo
      
c     re-decide the movement range of tremendicy:
      Mend_N=0                  !Mend_N points from 1 are moveable
      k=0
      do i=1,Lch
         k=k+1
         if(k.gt.Mend.or.mv1(i).lt.0) goto 111
         Mend_N=Mend_N+1
      enddo
 111  continue
      Mend_C=0                  !Mend_C points from Lch are moveable
      k=0
      do i=Lch,1,-1
         k=k+1
         if(k.gt.Mend.or.mv1(i).lt.0) goto 222
         Mend_C=Mend_C+1
      enddo
 222  continue

c     find 6-bond-move position -------------->
      nfl6=0                    !total number of 6-bond-move points
      do i=1,Lch                !i: fixed border of moving block
         if(mv1(i).gt.0)then
            if(mv1(i+1).gt.0)then
               if(mv1(i+2).gt.0)then
                  if(mv1(i+3).gt.0)then
                     if(mv1(i+4).gt.0)then
                        if(mv1(i+5).gt.0)then
                           if(mv1(i+6).gt.0)then !i+6: fixed border
                              nfl6=nfl6+1
                              ras6(nfl6)=i
                           endif
                        endif
                     endif
                  endif
               endif
            endif
         endif
      enddo

c     find 5-bond-move position -------------->
      nfl5=0                    !total number of 5-bond-move points
      do i=1,Lch                !i: fixed border of moving block
         if(mv1(i).gt.0)then
            if(mv1(i+1).gt.0)then
               if(mv1(i+2).gt.0)then
                  if(mv1(i+3).gt.0)then
                     if(mv1(i+4).gt.0)then
                        if(mv1(i+5).gt.0)then !i+5: fixed border
                           nfl5=nfl5+1
                           ras5(nfl5)=i
                        endif
                     endif
                  endif
               endif
            endif
         endif
      enddo

c     find 4-bond-move position -------------->
      nfl4=0                    !total number of 4-bond-move points
      do i=1,Lch                !i: fixed border of moving block
         if(mv1(i).gt.0)then
            if(mv1(i+1).gt.0)then
               if(mv1(i+2).gt.0)then
                  if(mv1(i+3).gt.0)then
                     if(mv1(i+4).gt.0)then !i+4: fixed border
                        nfl4=nfl4+1
                        ras4(nfl4)=i
                     endif
                  endif
               endif
            endif
         endif
      enddo

c     find 3-bond-move position -------------->
      nfl3=0                    !total number of 3-bond-move points
      do i=1,Lch                !i: fixed border of moving block
         if(mv1(i).gt.0)then
            if(mv1(i+1).gt.0)then
               if(mv1(i+2).gt.0)then
                  if(mv1(i+3).gt.0)then !i+3 fixed border
                     nfl3=nfl3+1
                     ras3(nfl3)=i
                  endif
               endif
            endif
         endif
      enddo

c     find 2-bond-move position -------------->
      nfl2=0                    !total number of 2-bond-move points
      do i=1,Lch                !i: fixed border of moving block
         if(mv1(i).gt.0)then
            if(mv1(i+1).gt.0)then
               if(mv1(i+2).gt.0)then !i+2 fixed border
                  nfl2=nfl2+1
                  ras2(nfl2)=i
               endif
            endif
         endif
      enddo
      
      do i=1,Lch
         mv(i)=mv1(i)
      enddo
      
c^^^^^^^^^^^^^^^^^^^^ move_point finished ^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
