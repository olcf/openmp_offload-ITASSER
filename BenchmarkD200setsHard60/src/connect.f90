      subroutine connect(i1,i2,pass)
      use params
      use lengths
      implicit integer(i-z)
!      parameter(ndim=1500)
!      parameter(nvec=416)
!      common/lengths/Lch,Lch1,Lch2
      common/initialstr/cx(ndim),cy(ndim),cz(ndim)
      common/mcheck_dis/amcheck_dis

      bdis0=3.5
      amcheck_dis=3.1            !3.1*0.87=2.7
      pass=1
!!$OMP target update to( i1 )      
! !$OMP target update to(i1, i2, pass) 
! !$OMP target update to(i2, i1, pass) !, pass)
      if(i1.eq.1)then           !!!!!N-terminal or whole structure, random walk
! !$OMP target update to(i2, pass)
         do i=i2,1,-1
            n_check=0
 10         call get_bond(axx,ayy,azz,3.8)
            cx(i)=cx(i+1)+axx
            cy(i)=cy(i+1)+ayy
            cz(i)=cz(i+1)+azz
            n_check=n_check+1
            if(n_check.gt.100000)then
               pass=3
               return
            endif
            if(mcheck(i,cx(i),cy(i),cz(i)).eq.3) goto 10 !excluded volumn
         enddo
      elseif(i2.eq.Lch)then     !!!!!!!!!!!!!!!!!!C_terminal,
         do i=i1,Lch
            n_check=0
 11         call get_bond(axx,ayy,azz,3.8)
            cx(i)=cx(i-1)+axx
            cy(i)=cy(i-1)+ayy
            cz(i)=cz(i-1)+azz
            n_check=n_check+1
            if(n_check.gt.100000)then
               pass=4
               return
            endif
            if(mcheck(i,cx(i),cy(i),cz(i)).eq.3) goto 11
         enddo
      else                      !!!!!!!!!!!!!interval
         diss=di(cx(i2+1),cy(i2+1),cz(i2+1),cx(i1-1),cy(i1-1),cz(i1-1))
         adis=diss/float((i2+1)-(i1-1))
***   linear connect for big gap---->
         adis0=3.5
         x_check=0
 15      x_check=x_check+1
         if(adis.ge.adis0)then
            dex=3.5*(cx(i1-1)-cx(i2+1))/diss !3.5*cos(thita)
            dey=3.5*(cy(i1-1)-cy(i2+1))/diss
            dez=3.5*(cz(i1-1)-cz(i2+1))/diss
            do j=i2,i1,-1
               cx(j)=cx(j+1)+dex
               cy(j)=cy(j+1)+dey
               cz(j)=cz(j+1)+dez
            enddo
            return              !end of connection
         endif
***   random walk from i1 to i2--------->
         m_check=0              !try 2000 times of whole walk
 13      m_check=m_check+1
         do 14 i=i1,i2
            n_check=0           !each point try 2000 times
 12         call get_bond(axx,ayy,azz,3.8)
            cx(i)=cx(i-1)+axx
            cy(i)=cy(i-1)+ayy
            cz(i)=cz(i-1)+azz
            n_check=n_check+1
            if(int(n_check/1000)*1000.eq.n_check)then 
               bdis0=bdis0*1.03
               if(bdis0.ge.3.8)bdis0=3.8
               amcheck_dis=2.8  !2.5*0.87=2.18
               if(i1.eq.i2)then !artificial gap, skip when initial model has error
                  amcheck_dis=2.0
               endif
            endif
            if(n_check.eq.4000)then
               goto 13
            endif
            if(m_check.ge.2000)then !can not pass the connection
               if(x_check.le.6)then
                  adis0=adis0*0.995
                  goto 15
               endif
               pass=5
               return
            endif
            if(mcheck(i,cx(i),cy(i),cz(i)).eq.3) goto 12 !check excluded V
            aaa=float(i2+1-i)
            bdis=di(cx(i),cy(i),cz(i),cx(i2+1),cy(i2+1),cz(i2+1))/aaa
            if(i.lt.i2.and.bdis.ge.bdis0) goto 12 !check remain steps
            if(i.eq.i2)then
               if(bdis.gt.4.2.or.bdis.lt.3.4) goto 12 !last step
            endif
            bdis0=3.5
            amcheck_dis=3.1
 14      continue
      endif

c^^^^^^^^^^^^^^^^ connect of gap [i1,i2] finished ^^^^^^^^^^^^^^^^^^
      return
      end
