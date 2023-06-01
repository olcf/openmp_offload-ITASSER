      subroutine read_CAcontact8
	use params
      use lengths
      use CAcontact
      implicit integer(i-z)
!      parameter(ndim=1500)
!      common/lengths/Lch,Lch1,Lch2
      common/zscore/izscore
!      common/CA8/McomCA8(ndim),KcomCA8(ndim,100),aweigCA8(ndim,ndim)
      common/CAcontact2/dist_CA_cut8
      DIMENSION r1(8000),r2(8000)
      
      dist_CA_cut8_tmp=8.0
      if(izscore.eq.1)then
         dist_CA_cut8=(dist_CA_cut8_tmp/0.87)**2
      elseif(izscore.eq.2)then
         dist_CA_cut8=(dist_CA_cut8_tmp/0.87)**2
      else
         dist_CA_cut8=(dist_CA_cut8_tmp/0.87)**2
      endif
      if(izscore.eq.1)then      !easy target
         cut_min=0.2
         cut0=0.5
      elseif(izscore.eq.2)then  !medium target
         cut_min=0.1
         cut0=0.4
      else                      !hard target
         cut_min=0.1
         cut0=0.3
      endif
ccc   pool all the restraints into r1(i),r2(i)--------------->
 11   continue
      read(28,*)ntmp
      i_c=0
      do i=1,ntmp
         read(28,*)i1,i2,conf
         if(conf.gt.1)conf=1
         if(conf.gt.cut_min)then
            i_c=i_c+1
            if(i_c.gt.8000)then
               write(*,*)'Too many COMB8CA restraints!!!!!'
               cut_min=cut_min*1.1
               rewind(28)
               goto 11
            endif
            r1(i_c)=i1
            r2(i_c)=i2
            if(conf.gt.cut0)then
               aweigCA8(i1,i2)=1+abs(conf-cut0)*4
            else
               aweigCA8(i1,i2)=1-abs(conf-cut0)*2
            endif
            aweigCA8(i2,i1)=aweigCA8(i1,i2)
         endif
      enddo
      NcomCA=i_c

ccc   map r1,2(i) into McomCA(i),KcomCA(i,McomCA(i))------------>
      do i=1,Lch
         McomCA8(i)=0            !number of contacts with 'i'
         do j=1,NcomCA
            if(r1(j).eq.i)then
               McomCA8(i)=McomCA8(i)+1
               KcomCA8(i,McomCA8(i))=r2(j)
            endif
            if(r2(j).eq.i)then
               McomCA8(i)=McomCA8(i)+1
               KcomCA8(i,McomCA8(i))=r1(j)
            endif
         enddo
      enddo
!$acc update device(McomCA8,KcomCA8,aweigCA8)

ccc   output restraints------------->
c      write(*,*)'Number of restraints:',NcomCA
c      write(*,*)'----------- CAcontact restraints ---------------'
c      nnc=0
c      do i=1,Lch
c         nnc=nnc+McomCA(i)
c         write(*,12)i,McomCA(i),(KcomCA(i,j),j=1,McomCA(i))
c 12      format(i4,'(',i2,'):',20i4)
c         write(*,13)i,McomCA(i),(aweigCA(i,KcomCA(i,j)),j=1,McomCA(i))
c 13      format(i4,'(',i2,'):',20f8.5)
c      enddo
c      write(*,*)'Number of CAcontact=',nnc,' Lch=',Lch
c      write(*,*)
c      stop

c^^^^^^^^^^^^^^^^^^ read CAcontact restrains finished ^^^^^^^^^^^^^^^^^
      return
      end
