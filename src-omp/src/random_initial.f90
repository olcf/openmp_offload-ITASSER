      subroutine random_initial
	use params
      use chain1
      use chainm
      use echain1
      use lengths
      implicit integer(i-z)
!      parameter(ndim=1500)
!      parameter(nrep=10)
!                parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
!      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
!      common/lengths/Lch,Lch1,Lch2
      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      common/commonuse2/atemp1,atemp2,N_rep,phot
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      dimension ip(ndim)
      common/sw3/icarep(ndim,nrep)
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      common/chain0/ras(ndim),nfl
      common/looks/exc,exc1,exc2
      common/ranzy/nozy
      common/nwmax1/nvecr

      do 102 k=1,N_REP
 88      x(1)=0
         y(1)=0
         z(1)=0
         m=0
         do 101 i=2,Lch
 99         ip(i)=int(aranzy(nozy)*nvecr)+1
            m=m+1
            if(m.gt.1000000)then
               write(*,*) 'UNSOLVABLE STERIC PROBLEM > EXIT_2'
               goto 88          !unsolvable steric problem
            endif
            if(i.gt.2)then      !check neighbor
               if(.not.goodc(ip(i-1),ip(i)))goto 99
            endif
            x(i)=x(i-1)+vx(ip(i))
            y(i)=y(i-1)+vy(ip(i))
            z(i)=z(i-1)+vz(ip(i))
            do j=1,i-1          !check excluded volumn for Ca
               ir=(x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2
               if(ir.lt.exc) goto 99
            enddo
 101     continue

         do i=1,Lch
            xrep(i,k)=x(i)
            yrep(i,k)=y(i)
            zrep(i,k)=z(i)
         enddo
 102  continue

******prepare movement for normal movement *****************
      nfl=Lch
      do i=1,nfl
         ras(i)=i
      enddo
      call move_point           !decide movement point for notmal run
      nfr=0                     !number of frozen fragments.

c^^^^^^^^^^^^ initial chains finished ^^^^^^^^^^^^^^^^^^^^^
      return
      end
