      subroutine read_quarsi3
      use params
      use lengths
      use pair
      use seqe
      use size
      use sizea
      use pair1
      implicit integer(i-z)
!      parameter(ndim=1500)
!                parameter(nvec=416)
      character*3 NAME
!      COMMON/size/ arla(0:19,0:19),arlm(0:19,0:19),arlp(0:19,0:19)
!      COMMON/sizea/ ala(0:19,0:19),alm(0:19,0:19),alp(0:19,0:19)
!      COMMON/pair/ apa(ndim,ndim),app(ndim,ndim),apm(ndim,ndim)
!      common/lengths/Lch,Lch1,Lch2
!      common/seqe/seq(ndim),sec(ndim)
!      common/pair1/eh2,eh1b,eh1c
      common/pairmap2/apabla(0:19,0:19),apablp(0:19,0:19)
      common/pairmap3/apablm(0:19,0:19),amap(ndim,ndim),nmap
      common/paircut/ash
      common/zscore/izscore
      common/fractpair/fract_pair1,fract_pair3
      common/pair33/i_pair3,mk_pair3
      common/weight/chuan

      dimension apba(ndim,ndim),apbp(ndim,ndim),apbm(ndim,ndim)

c     Pairwise interactions apablp ... and cut-off parmeters
c     arlp, orientation dependent, pairwise specific, sequence
c     independent

c     read contact-pair potential from 'quarsi3.comm' ------->
      read(3,*)
      read(3,*)
      do i=0,19
         read(3,725) NAME, (apablp(i,j),j=0,19) !for app
      enddo
      read(3,*)
      read(3,*)
      do i=0,19
         read(3,725) NAME, (apablm(i,j),j=0,19) !for apm
      enddo
      read(3,*)
      read(3,*)
      do i=0,19
         read(3,725) NAME, (apabla(i,j),j=0,19) !for apa
      enddo
c     read distance-range from 'quarsi3.comm' ------->
      read(3,*)
      read(3,*)
      do i=0,19
         read(3,725) NAME, (arlp(i,j),j=0,19) !max distance for parallel
      enddo
      read(3,*)
      read(3,*)
      do i=0,19
         read(3,725) NAME, (arlm(i,j),j=0,19) !for perpendicular contact
      enddo
      read(3,*)
      read(3,*)
      do i=0,19
         read(3,725) NAME, (arla(i,j),j=0,19) !for antiparellel pair
      enddo
 725  format(a3,1x,20f5.1)

c     width of energy-well: [arla,ala]
      ash=0.17
      ash_min=1-ash
      ash_max=1+ash
      do i=0,19
         do j=0,19
            ala(i,j)=(arla(i,j)*ash_max/0.87)**2
            alm(i,j)=(arlm(i,j)*ash_max/0.87)**2
            alp(i,j)=(arlp(i,j)*ash_max/0.87)**2
            arla(i,j)=(arla(i,j)*ash_min/0.87)**2
            arlm(i,j)=(arlm(i,j)*ash_min/0.87)**2
            arlp(i,j)=(arlp(i,j)*ash_min/0.87)**2
         enddo
      enddo
c     E=EH1/2, for r in [0,arlp];
c     E=app-es*fs,  for [0,alp];
c     E=0,     for r in [alp,00].
c^^^^^^^^^^^^^contact interaction range finished ^^^^^^^^^^^^^^^^^

*>>>>>>>>>>>>>
c     read from 'pair3.dat'-------------->
      i_pair3=-1                !without pair3.dat exist
      read(25,*,end=1000)
      rewind(25)
      i_pair3=1
      Nline=100                 !maximum Lch is 2500 !!!
      do i=1,Lch
         read(25,*)
         do i_line=1,Nline
            line_end=min(25*i_line,Lch) !25,50,75,100,...., ending point
            read(25,*)(apba(i,j),j=(i_line-1)*25+1,line_end)
            if(line_end.eq.Lch) go to 4074
         enddo
 4074    continue
         read(25,*)
      enddo
      do i=1,Lch
         read(25,*)
         do i_line=1,Nline
            line_end=min(25*i_line,Lch) !25,50,75,100,...., ending point
            read(25,*)(apbm(i,j),j=(i_line-1)*25+1,line_end)
            if(line_end.eq.Lch) go to 4075
         enddo
 4075    continue
         read(25,*)
      enddo
      do i=1,Lch
         read(25,*)
         do i_line=1,Nline
            line_end=min(25*i_line,Lch) !25,50,75,100,...., ending point
            read(25,*)(apbp(i,j),j=(i_line-1)*25+1,line_end)
            if(line_end.eq.Lch) go to 4076
         enddo
 4076    continue
         read(25,*)
      enddo
      close(25)
 1000 continue
*<<<<<<<<<<<<<<<<

c     combine the data in 'quarsi3.comm' and 'pair3.dat' to get
c     contact potential-------------------->
      do i=1,Lch
         ii=SEQ(i)
         do j=1,Lch
            jj=SEQ(j)
            if(iabs(i-j).lt.5) then
               dd=0.0
            else
               dd=0.25          !encourage contact of distant residues.
            endif

            apa(i,j)=apabla(ii,jj)-dd !quasi3.comm only
            apm(i,j)=apablm(ii,jj)-dd !quasi3.comm only
            app(i,j)=apablp(ii,jj)-dd !quasi3.comm only
            if(i_pair3.gt.0)then
               apa(i,j)=fract_pair3*apba(i,j)+(1-fract_pair3)*
     &              apabla(ii,jj)-dd
               apm(i,j)=fract_pair3*apbm(i,j)+(1-fract_pair3)*
     &              apablm(ii,jj)-dd
               app(i,j)=fract_pair3*apbp(i,j)+(1-fract_pair3)*
     &              apablp(ii,jj)-dd
            endif
         enddo
      enddo

c^^^^^^^^^^^^^^^ pair-potential is obtained ^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
