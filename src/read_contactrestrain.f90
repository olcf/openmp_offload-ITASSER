      subroutine read_contactrestrain
	use params
      use lengths
      use distres
      use RES
      use freg
      implicit integer(i-z)
!      parameter(ndim=1500)
!      parameter(nvec=416)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev    
!      common/lengths/Lch,Lch1,Lch2
!      common/distres/er4,es3c
!      COMMON/RES/ER3,er5,er6,er7,Mcom(ndim),Kcom(ndim,100)
      common/resnumber/Ncom,Ndis,accur
      common/cutoff/cut1a,cut2a,cut3a,cut4a,cut1b,cut2b,cut3b,cut4b
!      common/freg/aweig(ndim,ndim)
      common/zscore/izscore
      
      DIMENSION r1(8000),r2(8000)

c     READS Side group - side group contacts 
c     (from NMR or therading predictions or clusters)

      if(izscore.eq.1)then      !easy target
         cut_min=0.2
         cut0=0.6
      elseif(izscore.eq.2)then  !medium target
         cut_min=0.1
         cut0=0.5
      else                      !hard target
         cut_min=0.1
         cut0=0.4
      endif
ccc   pool all the restraints into r1(i),r2(i)--------------->
 11   continue
      read(16,*)ntmp
      i_c=0
      do i=1,ntmp
         read(16,*)i1,i2,conf
         if(conf.gt.1)conf=1
         if(conf.gt.cut_min)then
            i_c=i_c+1
            if(i_c.gt.8000)then
               write(*,*)'Too many COMB restraints!!!!!!!!'
               cut_min=cut_min*1.1
               rewind(16)
               goto 11
            endif
            r1(i_c)=i1
            r2(i_c)=i2
            if(conf.gt.cut0)then
               aweig(i1,i2)=1+abs(conf-cut0)*4
            else
               aweig(i1,i2)=1-abs(conf-cut0)*2
            endif
            aweig(i2,i1)=aweig(i1,i2)
c     write(*,*)i1,i2,conf,aweig(i1,i2)
         endif
      enddo
      Ncom=i_c

ccc   map r1,2(i) into Mcom(i),Kcom(i,Mcom(i))------------>
      do i=1,Lch
         Mcom(i)=0              !number of contacts with 'i'
         do j=1,Ncom
            if(r1(j).eq.i)then
               Mcom(i)=Mcom(i)+1
               Kcom(i,Mcom(i))=r2(j)
            endif
            if(r2(j).eq.i)then
               Mcom(i)=Mcom(i)+1
               Kcom(i,Mcom(i))=r1(j)
            endif
         enddo
      enddo

      colim=1.5*Ncom           !background number for derviation
c     the larger 'colim' is, the weaker the contact restrain is.

ccc   output restraints------------->
      write(20,*)'Number of restraints:',Ncom
      write(20,*)'----------- contact restraints ---------------'
      nnc=0
      do i=1,Lch
         nnc=nnc+Mcom(i)
c FAILS reason unknown
cmec!!         write(20,12)i,Mcom(i),(Kcom(i,j),j=1,Mcom(i))
 12      format(i4,'(',i2,'):',20i4)
      enddo
      write(20,*)'Number of contact=',nnc,' Lch=',Lch
      write(20,*)'fc=',float(nnc)/Lch
      write(20,*)

c^^^^^^^^^^^^^^^^^^ read contact restrains finished ^^^^^^^^^^^^^^^^^
!$acc update device(Mcom,Kcom(:,:),aweig) 
      return
      end
