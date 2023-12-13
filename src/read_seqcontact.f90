      subroutine read_seqcontact
      use params
      use lengths
      use svm1
      use svm2
      implicit integer(i-z)
!      parameter(ndim=1500)
!      common/lengths/Lch,Lch1,Lch2
      
!      common/svm1/Mcon(15,ndim),Kcon(15,ndim,100),awei(15,ndim,ndim)
!      common/svm2/acc_cut,dist_svm2(15),nk_svm(3),ik_svm(3,15)
      common/svm4/acc_cutB,acc_cutG
      DIMENSION r1(8000),r2(8000),conf0(4,15),idist(15)
      dimension mmm(3,15),nnn(3,15)
      
      acc_cut=acc_cut+0.01      !shift one bar
ccc   read conf.comm ------->
 10   read(50,*)acc0
      if(acc0+0.001.lt.acc_cut)then
         do i=1,4
            read(50,*)
            do j=1,13
               read(50,*)
            enddo
         enddo
         goto 10
      endif
      nk_svm(1)=0               !number of contact files: C-alpha (5)
      nk_svm(2)=0               !number of contact files: C-beta (3)
      nk_svm(3)=0               !number of contact files: SG (3)
      do k1=1,3                 !contact order is only needed in this subroutine
         read(50,*)
         do k2=1,13
            read(50,*)conf0(k1,k2),dist,it
c     adjust conf0:
            if(it.eq.2)then     !CB
               conf0(k1,k2)=conf0(k1,k2)*acc_cutB
            endif
            if(it.eq.3)then     !SG
               conf0(k1,k2)=conf0(k1,k2)*acc_cutG
            endif
            
c     record number of contact type:
            if(k1.eq.1)then
               idist(k2)=dist
               dist_svm2(k2)=(dist/0.87)**2
               nk_svm(it)=nk_svm(it)+1 !number of CA contacts
               ik_svm(it,nk_svm(it))=k2 !type of contact (13 contacts)
            endif
         enddo
      enddo
ccc   ^^^^^^^^^ read conf.com completed ^^^^^^^^^^

c      do i=1,3
c         do j=1,nk_svm(i)
c            write(*,*)i,j,ik_svm(i,j)
c         enddo
c      enddo
c      do i=1,13
c         write(*,*)i,idist(i),dist_svm2(i)
c      enddo
      
ccc   set-up minimum contact taken from each program
ccc   decided based on table at /home/yzhang/pdbinput/contact/readme
      nc0_5=Lch*0.05            !total is Lch*0.15
      nc1=Lch*0.1               !total is Lch*0.3
      nc1_7=Lch*0.17            !total is Lch/2=Lch*0.5
      nc2=Lch*0.2               !total is Lch*0.6
      do k1=1,3                 !short/medm/long
         do k2=1,13             !betacon, SPcon etc
            mmm(k1,k2)=0
            if(k1==1)then       !short-range, |i-j| in [6,11]
               if(idist(k2).eq.8.and.k2.ne.12)then !all except for psicov
                  mmm(k1,k2)=nc1_7
               endif
            elseif(k1==2)then   !medm-range, |i-j| in [12-24]
               if(k2.eq.3.or.k2.eq.6.or.k2.eq.13)then !svmseqca8,svmseqcb8,spcon
                  mmm(k1,k2)=nc1_7
               elseif(k2.eq.11)then !betacon with good acc but bad cov
                  mmm(k1,k2)=nc1
               endif
            else                !long-range, |i-j|>24
               if(k2.eq.12.or.k2.eq.13)then !psicov,spcon
                  mmm(k1,k2)=nc1_7
               elseif(k2.eq.3.or.k2.eq.6.or.k2.eq.9.or.k2.eq.10.or.
     &                 k2.eq.11)then !svmseqCA8,CB8,SG8,svmcon,betacon
                  mmm(k1,k2)=nc0_5
               endif
            endif
            mmm(k1,k2)=nint(mmm(k1,k2)*0.8)
         enddo
      enddo
      
ccc   read contact predictions ------->
      do k=1,13
         i_c=0
         do i=1,3
            nnn(i,k)=0
         enddo
         
         read(50+k,*,end=14)ntmp
         if(ntmp.gt.2000)ntmp=2000
         do i=1,ntmp
            
c     read contact:
            read(50+k,*)i1,i2,conf
            if(conf.gt.1)conf=1
            if(i1.gt.i2)then
               i3=i1
               i1=i2
               i2=i3
            endif
            
c     contact order:
            ico=i2-i1
            if(ico.gt.24)then   ![25-inf]
               jco=3
            elseif(ico.gt.11)then ![12,24]
               jco=2
            elseif(ico.gt.5)then ![6,11]
               jco=1
            else
               goto 120
            endif
            
c     decide if take the contact prediction
            mk=0                !donot take
            if(nnn(jco,k).lt.mmm(jco,k))then !specific cut for programs
               mk=1
            elseif(i_c.lt.nc2)then !total i_c > Lch*0.2
               mk=1
            elseif(conf.gt.conf0(jco,k))then !based on confidence score
               mk=1
            endif
            
c     record contact and weight
            if(mk.eq.1)then
               nnn(jco,k)=nnn(jco,k)+1
               i_c=i_c+1
               r1(i_c)=i1
               r2(i_c)=i2
               awei(k,i1,i2)=2.5*(1+(conf-conf0(jco,k)))
               awei(k,i2,i1)=awei(k,i1,i2)
            endif
 120        continue
         enddo
 14      NcomCA=i_c
         
ccc   map r1,2(i) into Mcon(k,i),Kcon(k,i,Mcon(i))------------>
         do i=1,Lch
            Mcon(k,i)=0         !number of contacts with 'i'
            do j=1,NcomCA
               if(r1(j).eq.i)then
                  Mcon(k,i)=Mcon(k,i)+1
                  Kcon(k,i,Mcon(k,i))=r2(j)
               endif
               if(r2(j).eq.i)then
                  Mcon(k,i)=Mcon(k,i)+1
                  Kcon(k,i,Mcon(k,i))=r1(j)
               endif
            enddo
         enddo
         
         goto 130
c     output restraints------------->
         write(*,*)'acc0,cutoff.conf_short,medm,long=',acc0,
     &        k,conf0(1,k),conf0(2,k),conf0(3,k)
         write(*,*)'Number of restraints:',NcomCA
         nnc=0
         do i=1,Lch
            nnc=nnc+Mcon(k,i)
            write(*,12)i,Mcon(k,i),(Kcon(k,i,j),j=1,Mcon(k,i))
            write(*,13)i,Mcon(k,i),(awei(k,i,Kcon(k,i,j)),j=1,Mcon(k,i))
 12         format(i4,'(',i2,'):',20i4)
 13         format(i4,'(',i2,'):',20f8.5)
         enddo
         write(*,*)'Number of CAcontact=',nnc,NcomCA*2,' Lch=',Lch
         write(*,*)'N_short=',nnn(1,k),mmm(1,k),k
         write(*,*)'N_medm =',nnn(2,k),mmm(2,k),k
         write(*,*)'N_long =',nnn(3,k),mmm(3,k),k
         write(*,*)
 130     continue
      enddo
      
c      stop
!$acc update device(Mcon,Kcon,dist_svm2,awei,nk_svm,ik_svm,acc_cut)
c^^^^^^^^^^^^^^^^^^ read CAcontact restrains finished ^^^^^^^^^^^^^^^^^
      return
      end
