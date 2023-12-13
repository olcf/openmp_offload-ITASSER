      subroutine secondary
	use params
      use lengths
      implicit integer(i-z)
!      parameter(ndim=1500)
!      parameter(nrep=10)
!      parameter(nvec=416)
      character*3 sequ
!      common/lengths/Lch,Lch1,Lch2
      common/consensus1/cx0(0:ndim),cy0(0:ndim),cz0(0:ndim),q(0:ndim)
      common/aminoacid/sequ(ndim)
      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)

      dimension sec_high(0:ndim)
      dimension nse_i(200),nse_f(200),nse_type(200)

      dimension cx00(0:ndim),cy00(0:ndim),cz00(0:ndim)
      dimension ax(ndim),ay(ndim),az(ndim)

ccc   RMSD:
      double precision r_1(3,ndim),r_2(3,ndim),r_3(3,ndim),w(ndim)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /ndim*1.0/
ccc   

ccc   read SSP --------------------------------->
      open(unit=23,file='seq.dat',status='old') !ssp from native
      do i=1,Lch
         read(23,707)k,name,sec_high(i)
c     write(*,*)i,k,name,sec_high(i)
      enddo
      close(23)
 707  format(i5,3x,a3,2i5)

ccc   find continuous SSP pieces -------------------->
      nse=0                     !number of fragments
      do i=1,200
         nse_i(i)=0
         nse_f(i)=0
      enddo
      sec_high(0)=1
      do i=1,Lch
         if(sec_high(i).ge.2.and.sec_high(i).ne.sec_high(i-1))then
            nse=nse+1
            nse_i(nse)=i        !initial location of nse'th ssp
            nse_type(nse)=sec_high(i) !type of ssp
         endif
         if(sec_high(i).ne.1)nse_f(nse)=i !final location of nse'th ssp
      enddo

ccc   remove small secondary structures -------------->
      nse_old=nse
      nse=0
      do i=1,nse_old
         len=nse_f(i)-nse_i(i)+1
         if(len.ge.L_cut)then
            nse=nse+1
            nse_i(nse)=nse_i(i)
            nse_f(nse)=nse_f(i)
            nse_type(nse)=nse_type(i)
         endif
      enddo

c      do i=1,nse
c         write(*,*)i,nse_i(i),nse_f(i),nse_f(i)-nse_i(i)+1,nse_type(i)
c      enddo

ccc   generate new (cx0,cy0,cz0) --------------------------->
      do i=1,Lch
         if(q(i).eq.1)then
            cx00(i)=cx0(i)
            cy00(i)=cy0(i)
            cz00(i)=cz0(i)
         endif
      enddo
      do 10 i=1,nse             !nse of SSP
         len=nse_f(i)-nse_i(i)+1 !length of SSP

ccc   generate standard fragments --------------------->
         if(nse_type(i).eq.2)then
            call alpha_helix(ax,ay,az,len) !alpha-helix
         else
            call beta_sheet(ax,ay,az,len) !beta-sheet
         endif
         
ccc   rotation matrix based on initial template ------->
*1: aligned
         n=0
         do j=nse_i(i),nse_f(i)
            if(q(j).eq.1)then
               n=n+1
               k=j-nse_i(i)+1
               r_1(1,n)=ax(k)
               r_1(2,n)=ay(k)
               r_1(3,n)=az(k)
               r_2(1,n)=cx00(j)
               r_2(2,n)=cy00(j)
               r_2(3,n)=cz00(j)
            endif
         enddo
*2: gapped at middle
         if(n.lt.2)then         !it is a gap at ssp region
            r_1(1,1)=ax(1)
            r_1(2,1)=ay(1)
            r_1(3,1)=az(1)
            r_1(1,2)=ax(nse_f(i)-nse_i(i)+1)
            r_1(2,2)=ay(nse_f(i)-nse_i(i)+1)
            r_1(3,2)=az(nse_f(i)-nse_i(i)+1)
            n=0
            j=nse_i(i)
 11         j=j-1
            if(q(j).ne.1.and.j.gt.1)goto 11
            if(q(j).eq.1)then
               n=n+1
               ax1=cx00(j)
               ay1=cy00(j)
               az1=cz00(j)
               j0=j
               mm=1
            endif
            j=nse_f(i)
 12         j=j+1
            if(q(j).ne.1.and.j.lt.Lch)goto 12
            if(q(j).eq.1)then
               n=n+1
               ax2=cx00(j)
               ay2=cy00(j)
               az2=cz00(j)
               mm=2
             endif
             if(n.eq.2)then
               aaa=sqrt((ax2-ax1)**2+(ay2-ay1)**2+(az2-az1)**2)
               bbb=sqrt((ax(len)-ax(1))**2+(ay(len)-ay(1))**2+
     &              (az(len)-az(1))**2)
               ccc=abs(nse_i(i)-j0)*1
               r_2(1,1)=ax1+(ax2-ax1)/aaa*ccc
               r_2(2,1)=ay1+(ay2-ay1)/aaa*ccc
               r_2(3,1)=az1+(az2-az1)/aaa*ccc
               r_2(1,2)=ax1+(ax2-ax1)/aaa*(ccc+bbb)
               r_2(2,2)=ay1+(ay2-ay1)/aaa*(ccc+bbb)
               r_2(3,2)=az1+(az2-az1)/aaa*(ccc+bbb)
             endif
         endif
*3: gap at terminal
         if(n.lt.2)then         !n=1, the ssp gap is at terminal
            n=2
            if(mm.eq.1)then     !C-terminal
*3a: C-terminal
               m=0
               j=nse_i(i)
 13            j=j-1
               if(q(j).eq.1)m=m+1
               if(m.eq.1.and.q(j).eq.1)then
                  j0=j
                  ax2=cx00(j)
                  ay2=cy00(j)
                  az2=cz00(j)
               endif
               if(m.lt.5)goto 13 !the minimum alignment length is 5
               if(m.eq.5.and.q(j).eq.1)then
                  ax1=cx00(j)
                  ay1=cy00(j)
                  az1=cz00(j)
               endif
               aaa=sqrt((ax2-ax1)**2+(ay2-ay1)**2+(az2-az1)**2)
               bbb=sqrt((ax(len)-ax(1))**2+(ay(len)-ay(1))**2+
     &              (az(len)-az(1))**2)
               ccc=abs(nse_i(i)-j0)*2
               r_2(1,1)=ax2+(ax2-ax1)/aaa*ccc
               r_2(2,1)=ay2+(ay2-ay1)/aaa*ccc
               r_2(3,1)=az2+(az2-az1)/aaa*ccc
               r_2(1,2)=ax2+(ax2-ax1)/aaa*(ccc+bbb)
               r_2(2,2)=ay2+(ay2-ay1)/aaa*(ccc+bbb)
               r_2(3,2)=az2+(az2-az1)/aaa*(ccc+bbb)
            else
*3a: N-terminal
               m=0
               j=nse_f(i)
 14            j=j+1
               if(q(j).eq.1)m=m+1
               if(m.eq.1.and.q(j).eq.1)then
                  j0=j
                  ax1=cx00(j)
                  ay1=cy00(j)
                  az1=cz00(j)
               endif
               if(m.lt.5)goto 14 !the minimum alignment length is 5
               if(m.eq.5.and.q(j).eq.1)then
                  ax2=cx00(j)
                  ay2=cy00(j)
                  az2=cz00(j)
               endif
               aaa=sqrt((ax2-ax1)**2+(ay2-ay1)**2+(az2-az1)**2)
               bbb=sqrt((ax(len)-ax(1))**2+(ay(len)-ay(1))**2+
     &              (az(len)-az(1))**2)
               ccc=abs(j0-nse_f(i))*2
               r_2(1,2)=ax1+(ax1-ax2)/aaa*ccc
               r_2(2,2)=ay1+(ay1-ay2)/aaa*ccc
               r_2(3,2)=az1+(az1-az2)/aaa*ccc
               r_2(1,1)=ax1+(ax1-ax2)/aaa*(ccc+bbb)
               r_2(2,1)=ay1+(ay1-ay2)/aaa*(ccc+bbb)
               r_2(3,1)=az1+(az1-az2)/aaa*(ccc+bbb)
            endif
         endif
c         write(*,*)'n=',n,'  i=',i
c         do j=1,n
c            write(*,*)j,r_2(1,j),r_2(2,j),r_2(3,j)
c         enddo
         call u3b(w,r_1,r_2,n,1,rms,u,t,ier) !u rotate r_1 to r_2
ccc   rotate ssp onto initial template ------------>
         do j=nse_i(i),nse_f(i)
            k=j-nse_i(i)+1
            cx0(j)=t(1)+u(1,1)*ax(k)+u(1,2)*ay(k)+u(1,3)*az(k)
            cy0(j)=t(2)+u(2,1)*ax(k)+u(2,2)*ay(k)+u(2,3)*az(k)
            cz0(j)=t(3)+u(3,1)*ax(k)+u(3,2)*ay(k)+u(3,3)*az(k)
         enddo
 10   continue

ccc   re-define q(i) ------------------------------>
      do i=1,Lch
         q(i)=0
      enddo
      do i=1,nse
         do j=nse_i(i),nse_f(i)
            q(j)=1
         enddo
      enddo

c      open(23,file='tmp',status='unknown')
c      do i=1,Lch
c         if(q(i).eq.1)then
c            write(23,1237)i,sequ(i),i,cx0(i),cy0(i),cz0(i)
c         endif
c      enddo
c 1237 format('ATOM  ',i5,'  CA  ',A3,I6,4X,3F8.3)
c      close(23)

c      stop
      return
      end
