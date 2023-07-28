      subroutine set_common
	use params
      use lengths
      use distres
      use pair1
      use svm2
      use one
      use RES
      implicit integer(i-z)
!      parameter(ndim=1500)
!      parameter(nrep=100)       !number of replicas
!      parameter(nvec=416)
      character protein*10
!      common/lengths/Lch,Lch1,Lch2
!      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1
      common/maxdis2/maxdis2(ndim)
      common/arandom/  aarand,abrand,acrand,adrand
      common/nswap/bNSa(100),bNSt(100)
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      common/commonuse1/random,ncycle
      common/commonuse2/atemp1,atemp2,N_rep,phot
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
      common/sw1/aT_rep(nrep),E_rep(nrep)
      character*6 mname
      character*80 line
      common/movename/mname(100)
      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)
      COMMON/RCN1/er1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
!      common/distres/er4,es3c
!      COMMON/RES/ER3,er5,er6,er7,Mcom(ndim),Kcom(ndim,100)
      common/zscore/izscore
      common/excluded/vvv(ndim,ndim)
      common/countres/t_comb,t_dist,t_combCA,s_comb,s_dist,s_combCA
      common/countres1/t_combCA8,t_distL,s_combCA8,s_distL,N_resc
!      common/pair1/eh2,eh1b,eh1c
      common/weight/chuan
      character type*10
      common/bigbond/i_bigbond,teco
      common/ssp/ssp
      common/fractpair/fract_pair1,fract_pair3
      common/pair33/i_pair3,mk_pair3
      common/ichos/ichos
      common/ranzy/nozy
      common/lattice/m_latt,latt1,latt2
      common/nwmax1/nvecr
      common/stick1/nstick,astick,nrmsd,ermsd
      common/stick2/iq(ndim,nrep)
      common/stick3/ax00(ndim,nrep),ay00(ndim,nrep),az00(ndim,nrep)
      common/stick6/bx00(ndim),by00(ndim),bz00(ndim),armsd_min0
      common/res2/er14,er15,er16,er17
!      common/svm2/acc_cut,dist_svm2(15),nk_svm(3),ik_svm(3,15)
      common/svm3/er21,er22,er23
      common/svm4/acc_cutB,acc_cutG
      
cccccccccccccc set the random generator cccccccccccccccccccccccc
      nozy=random
      if(nozy.gt.0)nozy=-nozy
      firstrandom=aranzy(nozy)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      read(21,*)
      read(21,*)Lch             !length of chain
      read(21,*)protein
      Lch1=Lch-1
      Lch2=Lch-2
      anvec=nvec-0.00001
      contt=1.5*float(Lch)      !TARGET NUMBER OF CONTACTS   1.5*N
      do i=1,Lch
         maxdis2(i)=20*i*i      !maximum distance of walk in i steps
      enddo
      write(20,*)'Target:',protein
      write(20,*)'Length:',Lch
      ichos=1

      if(m_latt.eq.1)then
         latt1=14
         latt2=25
         nvecr=312
      else
         latt1=12
         latt2=26
         nvecr=416
      endif
      anvec=nvecr-0.00001
      
      do i=1,Lch
         do j=1,Lch
            vvv(i,j)=1          !every pair should be checked
         enddo
      enddo

      armsd_min0=1000

c      if(Lch.ge.300)ncycle=int(ncycle*0.8)
c      if(Lch.ge.400)ncycle=int(ncycle*0.8)
c      if(Lch.ge.500)ncycle=int(ncycle*0.8)
c      if(Lch.ge.600)ncycle=int(ncycle*0.8)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      N_resc=0
      t_comb=0
      t_dist=0
      t_distL=0
      t_combCA=0
      t_combCA8=0
      s_comb=0
      s_dist=0
      s_distL=0
      s_combCA=0
      s_combCA8=0

ccccccccccccccccccccccc Temperature ccccccccccccccccccccccccccc
c     [80,130] is the standard:
      if(Lch.lt.80)then
         atemp1=atemp1*0.97
         atemp2=atemp2*0.91
         if(Lch.lt.55)then
            atemp1=atemp1*0.97
            atemp2=atemp2*0.91
         endif
      endif
      if(Lch.gt.130)then
         atemp1=atemp1*1.05
         atemp2=atemp2*1.1
         if(Lch.gt.165)then
            atemp1=atemp1*1.05
            atemp2=atemp2*1.1
            if(Lch.gt.200)then
               atemp1=atemp1*1.05
               atemp2=atemp2*1.1
            endif
         endif
      endif

ccccccccccccc Number of replicas #################
      if(Lch.gt.165)then        !50
         N_rep=N_rep+10
      endif
      if(Lch.gt.240)then        !60
         N_rep=N_rep+10
      endif
      if(Lch.gt.300)then        !70
         N_rep=N_rep+10
      endif
      if(Lch.gt.400)then      !80
         N_rep=N_rep+10
      endif
      if(N_rep.gt.80)N_rep=80

ccccccccccccccc movement name cccccccccccccccccccccccc
      mname(2)='move2a'
      mname(3)='move3s'
      mname(4)='move3d'
      mname(5)='move4s'
      mname(6)='move4d'
      mname(7)='move8'
      mname(8)='move5s'
      mname(9)='move5d'
      mname(10)='move6'
      mname(11)='move_n'
      mname(12)='move_c'
      mname(13)='move7a'
      mname(14)='move7b'
      mname(15)='move9'
      mname(16)='tran_N'
      mname(17)='tran_M'
      mname(18)='tran_C'
      mname(19)='rot_N' !no
      mname(20)='rot_M'
      mname(21)='rot_C' !no
      mname(22)='trot_N'
      mname(23)='trot_M'
      mname(24)='trot_C'
      mname(25)='defo_N'
      mname(26)='defo_M'
      mname(27)='defo_C'

ccccccccccreset restraints weights according to zscore cccccccccccc
      rewind(24)
      read(24,*)n_thr,type
      if(type.eq.'easy')then    !--------------->easy
         izscore=1
         er1=chuan*3.6          !for dist.dat
         er3=chuan*0.765        !for comb.dat
         er4=chuan*0.45         !for comb.dat of deviation
         er5=chuan*2.7          !for combCA.dat
         eh1c=chuan*1.8         !for par.dat
         er6=chuan*0.45         !for distL.dat
         er7=chuan*500          !for RMSD
         if(ssp.eq.1)er7=chuan*100 !for RMSD
         i_bigbond=3            !decide fragment base on distance
         fract_pair1=0.4
         fract_pair3=0.3

         acc_cut=0.4
         acc_cutB=1.1
         acc_cutG=1.1
         er21=2.5
         er22=2.5
         er23=2.5
      elseif(type.eq.'medm')then !--------------->medium
         izscore=2
         er1=chuan*4.05         !for dist.dat
         er3=chuan*0.81         !for comb.dat
         er4=chuan*0.405        !for comb.dat of deviation
         er5=chuan*1.08         !for combCA.dat
         eh1c=chuan*1.0         !for par.dat
         er6=chuan*1.0          !for distL.dat
         er7=chuan*100          !for RMSD
         if(ssp.eq.1)er7=chuan*10 !for RMSD
         i_bigbond=1            !decide fragment base on +-2
         fract_pair1=0.7
         fract_pair3=0.3

         acc_cut=0.4
         acc_cutB=1.1
         acc_cutG=1.1
         er21=2.5
         er22=2.5
         er23=2.5
      else                      !--------------->hard
         izscore=3
         er1=chuan*2.7          !for dist.dat
         er3=chuan*0.4          !for comb.dat
         er4=chuan*0.27         !for comb.dat of deviation
         er5=chuan*0.4          !for combCA.dat
         eh1c=chuan*1.5         !for par.dat
         er6=chuan*0.5          !for distL.dat
         er7=chuan*5            !for RMSD
         i_bigbond=3            !decide fragment base on distance
         fract_pair1=0.3
         fract_pair3=0.7
         
         acc_cut=0.375
         acc_cutB=1.1
         acc_cutG=1.1
         er21=2.0
         er22=2.0
         er23=2.0
      endif
      if(teco.eq.1) i_bigbond=1 !teco=1, template from teco
      if(chuan.le.0.0001)then   !without using restraints
         eh1c=1                 !for par.dat
         fract_pair1=1          !do not use par.dat
      endif
c      er14=er5
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c^^^^^^^^^^^^^^^^ common parameters finished ^^^^^^^^^^^^^^^^^^^^^^^      
      return
      end
