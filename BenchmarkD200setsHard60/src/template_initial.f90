      subroutine template_initial
	use params
      use chainm
      use chain1
      use echain1
      use lengths
      use distres
      use echain2
      use seqe
      use RES
      use aax_mod
      use aay_mod
      use aaz_mod
      implicit integer(i-z)
!      parameter(ndim=1500)
!      parameter(nrep=10)
!      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
!      use echain2
!      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
!      common/lengths/Lch,Lch1,Lch2
      common/xyzrs/xrs(ndim,40,40),yrs(ndim,40,40),zrs(ndim,40,40)
      common/commonuse2/atemp1,atemp2,N_rep,phot
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/rs/i_thr0
      common/initialinput/switch,k_cycle,k_phot,N_ann
      character c1*4,c2*2,c3*3,aaaaa*10,text,text1*22
      character*80 head
      character sign(ndim)
      common/initialstr/cx(ndim),cy(ndim),cz(ndim)
!      common/seqe/seq(ndim),sec(ndim)
      common/looks/exc,exc1,exc2

      COMMON/RCN1/er1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
!      common/distres/er4,es3c
!      COMMON/RES/ER3,er5,er6,er7,Mcom(ndim),Kcom(ndim,100)

      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      common/sw3/icarep(ndim,nrep)
      common/sw4/exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep)
      common/sw5/egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep)
      common/sw7/ecxrep(ndim,nrep),ecyrep(ndim,nrep),eczrep(ndim,nrep)
      common/sw8/ebxrep(ndim,nrep),ebyrep(ndim,nrep),ebzrep(ndim,nrep)
      common/sw9/etxrep(ndim,nrep),etyrep(ndim,nrep),etzrep(ndim,nrep)
      common/excluded/vvv(ndim,ndim)

      common/consensus1/cx0(0:ndim),cy0(0:ndim),cz0(0:ndim),q(0:ndim)

      common/mng/m_g(100)
      dimension n_i(ndim),n_f(ndim)
      dimension cx0b(ndim),cy0b(ndim),cz0b(ndim)
      dimension cx0g(ndim),cy0g(ndim),cz0g(ndim)
      dimension ax(ndim),ay(ndim),az(ndim)
      dimension x0(ndim),y0(ndim),z0(ndim)

      common/beta1/axalf(0:19),ayalf(0:19),azalf(0:19)
      common/beta2/axbet(0:19),aybet(0:19),azbet(0:19)
      common/beta3/bxalf(0:19),byalf(0:19),bzalf(0:19)
      common/beta4/bxbet(0:19),bybet(0:19),bzbet(0:19)

      dimension M_i(ndim),M_f(ndim)

      dimension mf(0:ndim)

      common/chain0/ras(ndim),nfl
      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)
      common/rotat/ar0(ndim),athita0(ndim),aphi0(ndim)

      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
!      common/echain1/ex(ndim),ey(ndim),ez(ndim)
!      common/echain2/egx(ndim),egy(ndim),egz(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/ssp/ssp

!      common/chainm/mv(ndim)
      dimension q_bk(ndim)

      character*3 sequ
      common/aminoacid/sequ(ndim)
      common/bigbond/i_bigbond,teco
      common/rmsdder/i_chunk(ndim),ex0(ndim),ey0(ndim),ez0(ndim)

      common/nwmax1/nvecr

      dimension loop_i(ndim),loop_f(ndim),loop_p(ndim)
      common/mloopf/mloop
      dimension q_old(ndim)

******************************************************
c     set-up the template to use
******************************************************
      exc_eff=exc               !for excluded volumn
      M_consensus=0             !not use consensus
      if(i_thr0.eq.0)then       !means use consensus of top-2 template
         i_thr0=1               !if failed, using the first template
         M_consensus=1          !use consensus in the following
      endif

****************************************************
c     read q(i), cx0(i) from 'init.dat'
****************************************************
      rewind(24)
      read(24,*)n_thr           !Number of real templates in 'init.dat'
      if(i_thr0.gt.n_thr)then
         write(20,*)'without threading structure at this i_thr0!'
         write(*,*)'without threading structure at this i_thr0!'
         stop
      endif
      do k=1,i_thr0
         do i=1,Lch
            q(i)=0
         enddo
         read(24,*)Nal
         do i=1,Nal
            read(24,1237)text,ii,text,a1,a2,a3
            cx0(ii)=a1
            cy0(ii)=a2
            cz0(ii)=a3
            q(ii)=1
         enddo
         read(24,*)text
      enddo
 1237 format(A22,I4,A4,3F8.3)
      if(M_consensus.eq.1)call get_consensus !get q(i),cx0(i) from consensus
c^^^^^^^^^^q(i) and cx0(i) obtained ^^^^^^^^^^^^^^^^^^^^

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     cut 5 loops if the template is full-length
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      mloop=5
c     read(19,*) mloop          !how many broken points if full-length
      write(20,*)'mloop=',mloop
      do i=1,Lch
         q_old(i)=q(i)
      enddo
***** check whether the template is full ------------->
      nq0=0                     !number of un-aligned regions
      nbigb=0                   !number of big-bonds
      do i=1,Lch
         if(q(i).eq.0)then
            nq0=nq0+1
         endif
         if(i.ge.2)then
            if(q(i-1).eq.1.and.q(i).eq.1)then
               dis=di(cx0(i-1),cy0(i-1),cz0(i-1),cx0(i),cy0(i),cz0(i))
               if(dis.gt.5.0)then
                  nbigb=nbigb+1
               endif
            endif
         endif
      enddo
      write(20,*)'number of un-aligned regions, nq0=',nq0
      write(20,*)'number of big-bonds, nbigb=',nbigb
      if(nq0.eq.0.and.nbigb.eq.0)then !it is a full length model
ccc   check loop position & loop number
         nloop=0
         if(sec(1).eq.1)then
            nloop=nloop+1
            loop_i(nloop)=1
            loop_f(nloop)=1
         endif
         do i=2,Lch
            if(sec(i).eq.1.and.sec(i-1).ne.sec(i))then
               nloop=nloop+1
               loop_i(nloop)=i
            endif
            if(sec(i).eq.1)then
               loop_f(nloop)=i
            endif
         enddo
         do i=1,nloop
            loop_p(i)=int(float(loop_f(i)+loop_i(i))/2.0)
         enddo
ccc   cut mloop by setting q(i)=0
         if(mloop.gt.0)then
            delt=float(nloop)/float(mloop+1)
            if(delt.lt.1)delt=1
            kk=1
            do i=1,nloop
               if(i.ge.delt*kk)then
                  m=loop_p(i)
                  q(m)=0
                  if(kk.ge.mloop)goto 50
                  kk=kk+1
               endif
            enddo
         endif
      endif
 50   continue
c^^^^^^^^^^^^^^^^ cut loop is finished ^^^^^^^^^^^^^^^^^^^

c      do i=1,Lch
c         write(*,*)i,sec(i),q_old(i),q(i)
c      enddo
c      stop
**********************************************************
c     using predicted secondary structures
**********************************************************
      if(ssp.eq.1)then
         call secondary         !redifine q(i) and cx0(i)
      endif

********************************************************
c     remove small segments, so that we have a better random walk
********************************************************
      M_a=0                     !number of aligned segments
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
c^^^^^^^^^^^ remove small segment finished ^^^^^^^^^^^^^

********************************************************
c     check and find GAP, only according to q(i)
********************************************************
      N_g=0                     !number of gaps
      q(0)=1
      do i=1,Lch
         if(q(i).ne.q(i-1).and.q(i).eq.0)then
            N_g=N_g+1
            n_i(N_g)=i          !initial point of the gap
         endif
         if(q(i).eq.0)n_f(N_g)=i !final point of the gap
      enddo
c^^^^^^^^^^^^^^^ check gap finished ^^^^^^^^^^^^^^^^^

      n_check_bond=0
 103  continue                  !re-walk and try to make dis(i,i+1)<7A
********************************************************
c     fill GAP: cx0(i) -> cx(i)
********************************************************
      n_walk=0                  !for number of reject by excluded volumn
 70   n_walk=n_walk+1
      if(n_walk.gt.100)then
         exc_eff=exc_eff*0.99
      endif
      if(n_walk.gt.10000)then !10000
         write(20,*)'unsolvable structure',pass,i,exc_eff
         write(*,*)'unsolvable structure',pass,i,exc_eff
         stop
      endif
      do i=1,Lch
         if(q(i).eq.1)then
            cx(i)=cx0(i)
            cy(i)=cy0(i)
            cz(i)=cz0(i)
         else
            cx(i)=1000000.      !for checking excluded volumn
            cy(i)=1000000.
            cz(i)=1000000.
         endif
      enddo
      do i=2,n_g
         i1=n_i(i)
         i2=n_f(i)
         call connect(i1,i2,pass) !fill missed cooridinates, no move others
         if(pass.ge.3)goto 70             !re-walk
      enddo
      if(n_g.ge.1)then
         i1=n_i(1)
         i2=n_f(1)
         call connect(i1,i2,pass)
         if(pass.ge.3)goto 70   !re-walk
      endif
*^^^^^^^^^^^^^^Fill gap done, cx(i) is continuous ^^^^^^^^^^^^^^^^^^^^^^^

**********************************************************************
c     decide moveable points nfl, ras(i), according to 3 factors:
**********************************************************************
***   1: ras => gap +- 2, smallest moveable-length is 3
      nfl=0                     !number of moveable points
      do i=1,Lch
         ras(i)=0               !residue name of i'th moveable point.
      enddo
      do i=1,n_g
         i1=n_i(i)-2
         i2=n_f(i)+2
         if(i1.lt.1)i1=1
         if(i2.gt.Lch)i2=Lch
         do ii=i1,i2
            nfl=nfl+1           !number of moveable residues
            ras(nfl)=ii         !locations
         enddo
      enddo
      if(nfl.gt.0)call sort_ras
      if(i_bigbond.eq.1)then    !only for medium
***   2a: ras => bigbond on  +-2 for templates:
         do i=2,Lch
            if(q(i-1).eq.1.and.q(i).eq.1)then
               bondlen=di(cx(i-1),cy(i-1),cz(i-1),cx(i),cy(i),cz(i))
               if(bondlen.gt.4.6)then
                  ik_f=i+1      !ending point
                  ik_i=i-2      !starting point of big-bond
                  if(ik_i.lt.1)ik_i=1
                  if(ik_f.gt.Lch)ik_f=Lch
                  do ik=ik_i,ik_f
                     nfl=nfl+1
                     ras(nfl)=ik
                  enddo
               endif
            endif
         enddo
         if(nfl.gt.0)call sort_ras
      else
***   2b: ras => bigbond +-1
         bond_max=0
         do i=2,Lch
            j=i-1               !starting point of big-bond
            bondlen=di(cx(i),cy(i),cz(i),cx(j),cy(j),cz(j))
            if(bondlen.gt.bond_max) bond_max=bondlen
            if(bondlen.gt.5.0)then
               do t=i,Lch
                  dh=di(cx(j),cy(j),cz(j),cx(t),cy(t),cz(t))
                  adh=dh/float(t-j)
                  if(adh.lt.3.498)goto 23 !can be walked to
               enddo
 23            ik_f=t+1         !ending point + 1
               ik_i=j-1         !starting point of big-bond - 1
               if(ik_i.lt.1)ik_i=1
               if(ik_f.gt.Lch)ik_f=Lch
               do ik=ik_i,ik_f
                  nfl=nfl+1
                  ras(nfl)=ik
               enddo
            endif
         enddo
         if(nfl.gt.0)call sort_ras
      endif
***   3: ras0 => smallpiece:
      do i=1,Lch
         mf(i)=0                !frozen
      enddo
      do i=1,nfl
         mf(ras(i))=1           !moveable points
      enddo
      nfr=0                     !number of frozen fragments
      mf(0)=1
      do i=1,Lch
         if(mf(i).ne.mf(i-1).and.mf(i).eq.0)then
            nfr=nfr+1
            nfr_i(nfr)=i        !starting point of nfr'th fragment
         endif
         if(mf(i).eq.0)nfr_f(nfr)=i !ending point of nfr'th fragment
      enddo
c     'nfr' frozen fragments in [nfr_i(i),nfr_f(i)], convert small piece ---->
      do i=1,nfr
         l_fr=nfr_f(i)-nfr_i(i)+1 !length of the frozen fragment
         if(l_fr.le.L_cut)then
            do j=nfr_i(i),nfr_f(i)
               nfl=nfl+1
               ras(nfl)=j
            enddo
         endif
      enddo
      if(nfl.gt.0)call sort_ras
*^^^^^^^^^^^^^^^^^^^^ nfl, ras(i) finished ^^^^^^^^^^^^^^^^^^^^^^

*****************************************************************
c     decide nfl2, ras2(i), mv(i):
*****************************************************************
      if(nfl.lt.1)then
         write(20,*)'without moveable points in this template!'
         write(*,*)'without moveable points in this template!'
      endif
      call move_point           !decide movement set from ras(i)
*^^^^^^^^^^ decision of movement set finished ^^^^^^^^^^^^^^^^^^^

****************************************************************************
c     decide nfr_i(i), nfr_f(i), i=1,...,nfr, according to nfl,ras(nfl):
****************************************************************************
      do i=1,Lch
         mf(i)=0
         sign(i)="f"            !frozen
      enddo
      do i=1,nfl
         mf(ras(i))=1
         sign(ras(i))="m"       !moveable points
      enddo
      mf(0)=1
      nfr=0                     !number of frozen fragments
      do i=1,100
         nfr_i(i)=0
         nfr_f(i)=0
      enddo
      do i=1,Lch
         if(mf(i).ne.mf(i-1).and.mf(i).eq.0)then
            nfr=nfr+1
            nfr_i(nfr)=i         !starting point of nfr'th fragment
         endif
         if(mf(i).eq.0)nfr_f(nfr)=i !ending point of nfr'th fragment
      enddo
c^^^^^^^^^^^^^^^^^^^ nfr, nfr_i(i), nfr_f(i) finished ^^^^^^^^^^^^^^^^^^

****************************************************************
c     redefine mv(i) for frozen fragments (for excluded volumn)
****************************************************************
      do i=1,nfr
         do j=nfr_i(i),nfr_f(i)
            mv(j)=-i
         enddo
      enddo
c^^^^^^^^^^^^^^^^^ mv(i) finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

****************************************************************
c     translation distance d_xyz00(i), rotation angle00(i)
      do i=1,nfr
         siz=nfr_f(i)-nfr_i(i)+1
         angle00(i)=angle0*10/float(siz)
         d_xyz00(i)=d_xyz0*10/float(siz)
         if(siz.gt.40)d_xyz00(i)=d_xyz00(i)/2.0
      enddo
c^^^ d_xyz00(i), angle00(i) done ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

****************************************************************
c     Record frozen parts:
c     there are 4 arraies needed to record and transfer:
c     common/sw4/exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep)    CA
c     common/sw5/egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep) SG
c     common/sw7/ecxrep(ndim,nrep),ecyrep(ndim,nrep),eczrep(ndim,nrep) cc
c     common/sw8/ebxrep(ndim,nrep),ebyrep(ndim,nrep),ebzrep(ndim,nrep) bb
c     common/sw9/etxrep(ndim,nrep),etyrep(ndim,nrep),etzrep(ndim,nrep) CB
****************************************************************
      cx0(0)=cx(1)+(cx(2)-cx(3))
      cy0(0)=cy(1)+(cy(2)-cy(3))
      cz0(0)=cz(1)+(cz(2)-cz(3))
      cx0(Lch+1)=cx(Lch)+(cx(Lch1)-cx(Lch2))
      cy0(Lch+1)=cy(Lch)+(cy(Lch1)-cy(Lch2))
      cz0(Lch+1)=cz(Lch)+(cz(Lch1)-cz(Lch2))
      do itemp=1,N_rep
         do i=1,Lch
            if(mv(i).lt.0)then  !frozen points
c     Ca:
               exrep(i,itemp)=cx0(i)/0.87
               eyrep(i,itemp)=cy0(i)/0.87
               ezrep(i,itemp)=cz0(i)/0.87
c     uniform vector:
               amx=cx0(i)-cx0(i-1)
               amy=cy0(i)-cy0(i-1)
               amz=cz0(i)-cz0(i-1)
               aaa=sqrt(amx**2+amy**2+amz**2)
               amx=amx/aaa      !ax(i)-ax(i-1)
               amy=amy/aaa
               amz=amz/aaa
               apx=cx0(i+1)-cx0(i)
               apy=cy0(i+1)-cy0(i)
               apz=cz0(i+1)-cz0(i)
               aaa=sqrt(apx**2+apy**2+apz**2)
               apx=apx/aaa      !ax(i+1)-ax(i)
               apy=apy/aaa
               apz=apz/aaa
               ang=acos(-(amx*apx+amy*apy+amz*apz))*180/3.1415926 !
               aaax=amx+apx
               aaay=amy+apy
               aaaz=amz+apz
               aaa=sqrt(aaax**2+aaay**2+aaaz**2)
               aaax=aaax/aaa
               aaay=aaay/aaa
               aaaz=aaaz/aaa
c     cc:
               ccx=amx-apx
               ccy=amy-apy
               ccz=amz-apz
               aaa=sqrt(ccx**2+ccy**2+ccz**2)
               ccx=ccx/aaa
               ccy=ccy/aaa
               ccz=ccz/aaa
               ecxrep(i,itemp)=ccx !cc
               ecyrep(i,itemp)=ccy
               eczrep(i,itemp)=ccz
c     Hb:
               bbx=amy*apz-amz*apy
               bby=amz*apx-amx*apz
               bbz=amx*apy-amy*apx
               aaa=sqrt(bbx**2+bby**2+bbz**2)
               bbx=bbx/aaa
               bby=bby/aaa
               bbz=bbz/aaa
               ebxrep(i,itemp)=bbx !Hb
               ebyrep(i,itemp)=bby
               ebzrep(i,itemp)=bbz
c     CB:
               k=seq(i)
               if(ang.lt.105)then !alpha
                  dx=(bxalf(k)*aaax+byalf(k)*bbx+bzalf(k)*ccx)/0.87
                  dy=(bxalf(k)*aaay+byalf(k)*bby+bzalf(k)*ccy)/0.87
                  dz=(bxalf(k)*aaaz+byalf(k)*bbz+bzalf(k)*ccz)/0.87
               else
                  dx=(bxbet(k)*aaax+bybet(k)*bbx+bzbet(k)*ccx)/0.87
                  dy=(bxbet(k)*aaay+bybet(k)*bby+bzbet(k)*ccy)/0.87
                  dz=(bxbet(k)*aaaz+bybet(k)*bbz+bzbet(k)*ccz)/0.87
               endif
               etxrep(i,itemp)=exrep(i,itemp)+dx
               etyrep(i,itemp)=eyrep(i,itemp)+dy
               etzrep(i,itemp)=ezrep(i,itemp)+dz
c     SG:
               k=seq(i)
               if(ang.lt.105)then !alpha
                  dx=(axalf(k)*aaax+ayalf(k)*bbx+azalf(k)*ccx)/0.87
                  dy=(axalf(k)*aaay+ayalf(k)*bby+azalf(k)*ccy)/0.87
                  dz=(axalf(k)*aaaz+ayalf(k)*bbz+azalf(k)*ccz)/0.87
               else
                  dx=(axbet(k)*aaax+aybet(k)*bbx+azbet(k)*ccx)/0.87
                  dy=(axbet(k)*aaay+aybet(k)*bby+azbet(k)*ccy)/0.87
                  dz=(axbet(k)*aaaz+aybet(k)*bbz+azbet(k)*ccz)/0.87
               endif
               egxrep(i,itemp)=exrep(i,itemp)+dx
               egyrep(i,itemp)=eyrep(i,itemp)+dy
               egzrep(i,itemp)=ezrep(i,itemp)+dz
            endif
         enddo
      enddo
c     current ex(i) for check of excluded volumn:
      do i=1,Lch
         if(mv(i).lt.0)then     !frozen points
            ex(i)=cx0(i)/0.87
            ey(i)=cy0(i)/0.87
            ez(i)=cz0(i)/0.87
         endif
      enddo
c^^^^^^^^^^^^^^^ record of frozen-fragment are set^^^^^^^^^^^^^^^^^^^^^^

*************************************************************************
c     project chain onto lattices, decide (x,y,z):
*************************************************************************
***   ax(i)=cx(i)/0.87, for the decision of (x,y,z):
      do i=1,Lch
         ax(i)=cx(i)/0.87       !C_alpha scaled by 0.87
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
         dis_min=100000000000.  !minimum distance between c_alpha and lattice
         j_ch=0                 !chosen vector
         do 100 j=1,nvecr
            if(i.gt.2)then      !check good neighbor
               if(.not.goodc(jm,j))goto 100
            endif
            x_tmp=xm+vx(j)
            y_tmp=ym+vy(j)
            z_tmp=zm+vz(j)
c     check excluded volumn---->
            if(mv(i).gt.0)then  !on-lattice
               do m=1,Lch       !!!!! off-lattice part.
                  if(abs(i-m).ge.3.and.mv(m).lt.0)then
                     disaa=(x_tmp-ex(m))**2+(y_tmp-ey(m))**2+
     $                    (z_tmp-ez(m))**2
                     if(disaa.lt.exc_eff) goto 100
                  endif
               enddo
               do m=1,i-3       !!!!! on-lattice part.
                  if(mv(m).gt.0)then
                     disaa=(x_tmp-x(m))**2+(y_tmp-y(m))**2+
     $                    (z_tmp-z(m))**2
                     if(disaa.lt.exc_eff) goto 100
                  endif
               enddo
            endif
c     ^^^^^^^^^^^^^^^^^^^^^^^^^^
            dis=(ax(i)-x_tmp)**2+(ay(i)-y_tmp)**2+(az(i)-z_tmp)**2
            if(dis.lt.dis_min)then
               j_ch=j
               x_ch=x_tmp
               y_ch=y_tmp
               z_ch=z_tmp
               dis_min=dis
            endif
 100     continue
         if(j_ch.lt.1)goto 70   !refill the gaps
         x(i)=x_ch              !get (x(i),y(i),z(i)) here
         y(i)=y_ch
         z(i)=z_ch
         jm=j_ch
         xm=x(i)
         ym=y(i)
         zm=z(i)
 101  continue
c^^^^^^^^^^^^^^^^^^^project lattice done ^^^^^^^^^^^^^^^^^^^^^^^

***********************************************************
*     check the bond-length
***********************************************************
      do i=1,Lch1
         r2=
     $        (nint(aax(i))-nint(aax(i+1)))**2+
     $        (nint(aay(i))-nint(aay(i+1)))**2+
     $        (nint(aaz(i))-nint(aaz(i+1)))**2
         if(r2.gt.65)then       !dis>7A
            n_check_bond=n_check_bond+1
            if(n_check_bond.lt.100)goto 103
         endif
c         write(*,*)i,r2,sqrt(float(r2))*.87
      enddo
cccc  

*************************************************************************
c     record replcas of (x,y,z), icarep:
*************************************************************************
      do itemp=1,N_rep
         do i=1,Lch
            xrep(i,itemp)=x(i)  !x(i)
            yrep(i,itemp)=y(i)
            zrep(i,itemp)=z(i)
            if(i.lt.Lch)then
               wx=x(i+1)-x(i)   !ica(i)
               wy=y(i+1)-y(i)
               wz=z(i+1)-z(i)
               icarep(i,itemp)=vector(wx,wy,wz)
            endif
         enddo
         icarep(Lch,itemp)=icarep(Lch2,itemp) !just for unity
      enddo

      call get_vvv !get vvv(i,j). vvv(i,j)>0 check; vvv(i,j)<0 not check

***********************************************************************
c      decide the frequence of bulk movements: f=1/n_fra
***********************************************************************
      n_fra=nint(8*sqrt(float(Lch)/100.0)) !move bulk once in n_fra local move
cccccc n_fra decrease with number of frozen residues
      nnn=0
      do i=1,Lch
        if(mv(i).lt.0)nnn=nnn+1
      enddo
      n_fra=nint(n_fra*(float(Lch)/(2.0*float(nnn)+0.001)))
cccccc n_fra decrease with number of frozen fragments:
      np=1
      if(nfr.eq.1)np=2
      n_fra=nint(n_fra*(2.8/(float(nfr)+0.001))**np)
      if(n_fra.lt.1) n_fra=1 !nfr big and nnn big
      if(n_fra.gt.20) n_fra=20 !nfr small and nnn small
c^^^^^^^^^^^ n_fra decided ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

********************************************************************
*     output the threading results:
********************************************************************
      write(20,*)'Lch=',Lch
      write(20,*)"#align=",Nal,'   #unalign=',Lch-Nal,"n_g=",n_g
      write(20,*)'#frozen=',Lch-nfl,'    #moveable=',nfl
      write(20,*)
      write(20,*)'number of frozen pieces: nfr=',nfr
      do i=1,nfr
         siz=nfr_f(i)-nfr_i(i)+1
         write(20,42)i,siz,nfr_i(i),nfr_f(i),d_xyz00(i),angle00(i)*57.3
      enddo
 42   format(i5,i5,' [',i3,',',i3,']',f8.3,f8.3)
      write(20,*)'---number of local movement for each bulk move:',n_fra
      write(20,*)
      write(20,*)'#chain fix/move  SEC   SEQ'
      do i=1,Lch
         write(20,41)i,sign(i),sec(i),sequ(i)
      enddo
 41   format(i7,a5,i8,a8)
      do i=1,Lch
         do j=1,Lch
            if(vvv(i,j).eq.-2)then
               write(20,*)i,j,'  not be checked for excluded volumn'
            endif
         enddo
      enddo
      write(20,*)'bond_max=',bond_max

*******for calculation of RMSD of fragment------------------------>
      do i=1,Lch
         i_chunk(i)=-1
      enddo
      do i=1,nfr
         do j=nfr_i(i),nfr_f(i)
            i_chunk(j)=i
            ex0(j)=ex(j)
            ey0(j)=ey(j)
            ez0(j)=ez(j)
         enddo
      enddo
*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

c      open(23,file='tmp',status='unknown')
c      do i=1,Lch
c         write(*,*)i,mv(i),q(i)
c         if(mv(i).le.0)then
c            write(23,1037)i,sequ(i),i,aax(i)*0.87,aay(i)*0.87,aaz(i)*0.87
c         endif
c      enddo
c 1037 format('ATOM  ',i5,'  CA  ',A3,I6,4X,3F8.3)
c      close(23)

c      stop
c^^^^^^^^^^^^^^^^^^^^^^^^^ template_initial finished ^^^^^^^^^^^^^^^^^^^^^
      return
      end
