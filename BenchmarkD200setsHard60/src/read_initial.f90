      subroutine read_initial
	use params
      use chainm
      use chain1
      use echain1
      use lengths
      use echain2
      use seqe
      use stick
      implicit integer(i-z)
!      parameter(ndim=1500)
!      parameter(nrep=10)
!      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
!      common/lengths/Lch,Lch1,Lch2
      common/commonuse2/atemp1,atemp2,N_rep,phot
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/initialinput/switch,k_cycle,k_phot,N_ann
      character c1*4,c2*2,c3*3,aaaaa*10,text,text1*22
      character*80 head
      character sign(ndim)
      common/initialstr/cx(ndim),cy(ndim),cz(ndim)
!      common/seqe/seq(ndim),sec(ndim)
      common/looks/exc,exc1,exc2

      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      dimension ax(ndim),ay(ndim),az(ndim)
      dimension x0(ndim),y0(ndim),z0(ndim)
      dimension M_i(ndim),M_f(ndim)

      common/chain0/ras(ndim),nfl
      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6

      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
!      common/echain2/egx(ndim),egy(ndim),egz(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C

c      common/chainm/mv(ndim)
      dimension q_bk(ndim)

      character*3 sequ
      common/aminoacid/sequ(ndim)
      common/bigbond/i_bigbond,teco
      common/rmsdder/i_chunk(ndim),ex0(ndim),ey0(ndim),ez0(ndim)

      common/consensus1/cx0(0:ndim),cy0(0:ndim),cz0(0:ndim),q(0:ndim)
      dimension n_i(ndim),n_f(ndim)
      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)

c      common/stick1/nstick,astick,nrmsd,ermsd
c      common/stick2/iq(ndim,nrep)
c      common/stick3/ax00(ndim,nrep),ay00(ndim,nrep),az00(ndim,nrep)
      common/nwmax1/nvecr

      i_rep=0                   !replicas

 104  rewind(24)
      read(24,*)n_thr           !Number of real templates in 'init.dat'
      do 102 i_init=1,n_thr

****************************************************
c     read q(i), cx0(i) from 'init.dat'
****************************************************
         do i=1,Lch
            q(i)=0
         enddo
         read(24,*)L_ali
         do i=1,L_ali
            read(24,1237)text,ii,text,a1,a2,a3
            cx0(ii)=a1
            cy0(ii)=a2
            cz0(ii)=a3
            q(ii)=1
         enddo
         read(24,*)text         !TER
 1237    format(A22,I4,A4,3F8.3)
c^^^^^^^^^^q(i) and cx0(i) obtained ^^^^^^^^^^^^^^^^^^^^

********************************************************
c     remove small segments, so that we have a better random walk
********************************************************
         M_a=0                  !number of aligned segments
         q(0)=0
         do i=1,Lch
            if(q(i).ne.q(i-1).and.q(i).eq.1)then
               M_a=M_a+1
               M_i(M_a)=i
            endif
            if(q(i).eq.1)M_f(M_a)=i
         enddo
         do i=1,M_a
            L_a=M_f(i)-M_i(i)+1
            if(L_a.le.L_cut)then
               do j=M_i(i),M_f(i)
                  q(j)=0
               enddo
            endif
         enddo
c^^^^^^^^^^^remove small segment finished ^^^^^^^^^^^^^

********************************************************
c     check and find GAP, only according to q(i)
********************************************************
         N_g=0                  !number of gaps
         q(0)=1
         do i=1,Lch
            if(q(i).ne.q(i-1).and.q(i).eq.0)then
               N_g=N_g+1
               n_i(N_g)=i       !initial point of the gap
            endif
            if(q(i).eq.0)n_f(N_g)=i !final point of the gap
         enddo
c^^^^^^^^^^^^^^^check gap finished ^^^^^^^^^^^^^^^^^

********************************************************
c     fill GAP: cx0(i) -> cx(i)
********************************************************
         n_walk=0               !for number of reject by excluded volumn
         exc_eff=exc            !for excluded volumn
 70      n_walk=n_walk+1
         if(n_walk.gt.100)then
            exc_eff=exc_eff*0.99
         endif
         if(n_walk.gt.10000)then  !10000
            write(20,*)'unsolvable structure',i
            write(*,*)'unsolvable structure',i
            stop
         endif
         do i=1,Lch
            if(q(i).eq.1)then
               cx(i)=cx0(i)
               cy(i)=cy0(i)
               cz(i)=cz0(i)
            else
               cx(i)=1000000.   !for checking excluded volumn
               cy(i)=1000000.
               cz(i)=1000000.
            endif
         enddo
         do i=2,n_g
            i1=n_i(i)
            i2=n_f(i)
            call connect(i1,i2,pass) !fill missed cooridinates, no move others
            if(pass.ge.3)goto 70 !re-walk
         enddo
         if(n_g.ge.1)then
            i1=n_i(1)
            i2=n_f(1)
            call connect(i1,i2,pass)
            if(pass.ge.3)goto 70 !re-walk
         endif
*^^^^^^^^^^^^^^Fill gap done, cx(i) is continuous ^^^^^^^^^^^^^^^^^^^^^^^

*************************************************************************
c     project chain onto lattices, decide (x,y,z):
*************************************************************************
***   ax(i)=cx(i)/0.87, for the decision of (x,y,z):
         do i=1,Lch
            ax(i)=cx(i)/0.87    !C_alpha scaled by 0.87
            ay(i)=cy(i)/0.87
            az(i)=cz(i)/0.87
         enddo
c     (ax,ay,az) into (x,y,z)----------------------->
         x(1)=nint(ax(1))
         y(1)=nint(ay(1))
         z(1)=nint(az(1))
         xm=x(1)
         ym=y(1)
         zm=z(1)
         do 101 i=2,Lch
            dis_min=100000000000. !minimum distance between c_alpha and lattice
            j_ch=0              !chosen vector
            do 100 j=1,nvecr
               if(i.gt.2)then   !check good neighbor
                  if(.not.goodc(jm,j))goto 100
               endif
               x_tmp=xm+vx(j)
               y_tmp=ym+vy(j)
               z_tmp=zm+vz(j)
c     check excluded volumn---->
               do m=1,i-3       !!!!! on-lattice part.
                  disaa=(x_tmp-x(m))**2+(y_tmp-y(m))**2+(z_tmp-z(m))**2
                  if(disaa.lt.exc_eff) goto 100
               enddo
c     ^^^^^^^^^^^^^^^^^^^^^^^^^^
               dis=(ax(i)-x_tmp)**2+(ay(i)-y_tmp)**2+(az(i)-z_tmp)**2
               if(dis.lt.dis_min)then
                  j_ch=j
                  x_ch=x_tmp
                  y_ch=y_tmp
                  z_ch=z_tmp
                  dis_min=dis
               endif
 100        continue
            if(j_ch.lt.1)goto 70 !refill the gaps
            x(i)=x_ch           !get (x(i),y(i),z(i)) here
            y(i)=y_ch
            z(i)=z_ch
            jm=j_ch
            xm=x(i)
            ym=y(i)
            zm=z(i)
 101     continue
c^^^^^^^^^^^^^^^^^^^project lattice done ^^^^^^^^^^^^^^^^^^^^^^^

***   record the initial conformation of k'th replica--->
         i_rep=i_rep+1
         do i=1,Lch
            xrep(i,i_rep)=x(i)  !ica will be calculated in set_current
            yrep(i,i_rep)=y(i)
            zrep(i,i_rep)=z(i)
            iq(i,i_rep)=q(i)
            if(q(i).eq.1)then
               ax00(i,i_rep)=cx0(i)/0.87
               ay00(i,i_rep)=cy0(i)/0.87
               az00(i,i_rep)=cz0(i)/0.87
            endif
         enddo
         if(i_rep.ge.N_rep)goto 105
 102  continue
      if(i_rep.lt.N_rep)goto 104

******prepare movement for normal movement *****************
 105  nfl=Lch
      do i=1,nfl
         ras(i)=i
      enddo
      call move_point           !decide movement point for notmal run
      nfr=0                     !number of frozen fragments.

!$acc update device(ax00,ay00,az00,iq,x,y,z)

c      do i=1,nfl
c         write(*,*)i,ras(i),ras2(i),ras3(i),ras4(i),ras5(i)
c      enddo
c      open(unit=70,file='initial.pdb',status='unknown')
c      write(70,*)N_rep
c      do j=1,N_rep
c         write(70,*)Lch
c         do i=1,Lch
c            write(70,1037)i,sequ(i),i,xrep(i,j)*0.87,yrep(i,j)*0.87,zrep(i,j)*0.87
c         enddo
c         write(70,*)'TER'
c      enddo
c 1037 format('ATOM  ',i5,'  CA  ',A3,I6,4X,3F8.3)
c      close(70)
c      stop

c^^^^^^^^^^^^ read initial chain finished ^^^^^^^^^^^^^^^^^^^^^
      return
      end
