              program TASSER
              use params
              !use openacc
              use backup2
              use hb
              use sizea
              use size
              use bisec
              use chainm
              use chain1
              use echain1
              use short1
              use lengths
              use ENERGY
              use hba
              use hbb
              use hb
              use pair
              use icgg
              use sg
              use distres
              use ehbenergy
              use concutt
              use pair1
              use shortcom
              use echain2
              use echain4
              use order
              use echain5
              use seqe
              use one
              use RCN
              use short
              use RES
              use svm1
              use svm2
              use fr
              use stick
              use short
              use shape
              use backup1
              use ehbc
              use trackn 
              implicit integer(i-z)
              common/logica/goodc
              logical look, goodc(nvec,nvec)
              common/three/angle(nvec,nvec)
        !      common/lengths/Lch,Lch1,Lch2
        !      common/seqe/seq(ndim),sec(ndim)
        common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
        !      common/bisec/cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
        !      common/hb/hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
              common/cutoff/cut1a,cut2a,cut3a,cut4a,cut1b,cut2b,cut3b,cut4b
        !      common/hopp/eonehw(0:19)
              common/looks/exc,exc1,exc2
        !      common/hba/eh5a,Cr2a,acut_bb,acut_cc,acut_vv,acut_hh
        !      common/hbb/eh5b,Cr2b,bcut_bb,bcut_cc,bcut_vv,bcut_hh
        !      common/concutt/concut(0:19,0:19),concut2(0:19,0:19),concut_sc

              common/arandom/  aarand,abrand,acrand,adrand
        !      COMMON/ENERGY/EH5,ES3,ES3a,ES3b,EH1,EH1a,EH4,EHBIJ(ndim,ndim)
        !      COMMON/pair/ apa(ndim,ndim),app(ndim,ndim),apm(ndim,ndim)
        !      COMMON/size/ arla(0:19,0:19),arlm(0:19,0:19),arlp(0:19,0:19)
        !      COMMON/short/ IBIN(-300:300),asr(ndim,-12:12) !safe when vr^2<30
        !    COMMON/short1/ JBIN(0:500),bsr(ndim,16),acops(ndim,16)
        !      COMMON/sizea/ ala(0:19,0:19),alm(0:19,0:19),alp(0:19,0:19)
        !      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1
        !      common/shape/amx,amy,amz,afs(ndim),afsn(ndim) !  common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1 !  common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
        !      common/fr/frga(ndim),frgb(ndim)
        !      COMMON/RES/ER3,er5,er6,er7,Mcom(ndim),Kcom(ndim,100)
        !      COMMON/RCN/Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
              COMMON/RCN1/er1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
        !      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
        !      common/shortcom/eh3,es4,es5,es6,es7,es7a,es7b,es7c
        !      common/sg/gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
              common/maxi/maxin,vect1,vect2
        !      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
              common/forpreparemove4/ asrr(0:19,0:19,-12:12)
              common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev    
        !      common/distres/er4,es3c
              common/rmsdrange/nca1,nca2
              common/CA/dx(ndim),dy(ndim),dz(ndim)
              common/msichores/msicho
        !      common/ehbenergy/EHB1,EHB1a,EHB1b,EHB1c,EHB2,EHB3,EHB4,EHB5,EHB6
        !      common/ehbenergy1/EHB5a,EHB5b
      common/eshortenergy1/ESHORT1,ESHORT2,ESHORT3,ESHORT4,ESHORT11
      common/eshortenergy2/ESHORT4a,ESHORT5,ESHORT5a,ESHORT5b,ESHORT5c
      common/eshortenergy3/ESHORT6,ESHORT7,ESHORT8,ESHORT9,ESHORT10
      common/eshortenergy4/ESHORT12
      common/otherenergy/E_cord,E_cnum
      common/resnumber/Ncom,Ndis,accur
      common/initialinput/switch,k_cycle,k_phot,N_ann
      common/thrtem/aT_ann(100)
      common/outputxyz/fxyz(3,ndim)

      common/temperature/itemp,atemp
      common/beta1/axalf(0:19),ayalf(0:19),azalf(0:19)
      common/beta2/axbet(0:19),aybet(0:19),azbet(0:19)
!      common/pair1/eh2,eh1b,eh1c
      dimension E_s(nrep),E_ss(nrep)
      common/paircut/ash
      common/countres/t_comb,t_dist,t_combCA,s_comb,s_dist,s_combCA
      common/countres1/t_combCA8,t_distL,s_combCA8,s_distL,N_resc

      common/nswap/bNSa(100),bNSt(100)
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)

!      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
!      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
!      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
!      common/backup2/eprofo,eprofn,energ
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/pairmap2/apabla(0:19,0:19),apablp(0:19,0:19)
      common/pairmap3/apablm(0:19,0:19),amap(ndim,ndim),nmap

      common/moveretio1/h2,h3s,h3d,h4s,h4d,h5s,h5d,h6,h7,hend
      common/moveretio2/bh2,bh3s,bh3d,bh4s,bh4d,bh5s,bh5d,bh6
      common/moveretio3/bh7a,bh7b,bhendn,bhendc
!      common/icgg/ icg(ndim), EH6  
      common/rs/i_thr0
      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)
      common/nrepfile/n_repf

      common/commonuse1/random,ncycle
      common/commonuse2/atemp1,atemp2,N_rep,phot
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/mng/m_g(100)
      common/acct/accept0
      character*6 mname
      character fn
      common/movename/mname(100)
      common/readinitial/m_initial
      common/eall/N_sum(100),energ_sum(100),energ_sum2(100),E_min

      common/chain0/ras(ndim),nfl
!      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
!      common/echain1/ex(ndim),ey(ndim),ez(ndim) !CA
!      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
!      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
!      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
!      common/chainm/mv(ndim)

      common/sw1/aT_rep(nrep),E_rep(nrep)
      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      common/sw3/icarep(ndim,nrep)
      common/sw4/exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep) !Ca
      common/sw5/egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep) !SG
      common/sw7/ecxrep(ndim,nrep),ecyrep(ndim,nrep),eczrep(ndim,nrep) !cc
      common/sw8/ebxrep(ndim,nrep),ebyrep(ndim,nrep),ebzrep(ndim,nrep) !hb
      common/weight/chuan
      common/bigbond/i_bigbond,teco
      common/ssp/ssp
      common/defoangle/defo_angle
      common/fractpair/fract_pair1,fract_pair3
      common/zscore/izscore
      common/aTs/aTs1,aTs2,aTs_rep(nrep),aTs
      common/aTTs/aTTs1,aTTs2,aTTs_rep(nrep),aTTs
      common/ichos/ichos
      common/nana1/nana
      common/ranzy/nozy
      common/hours/hour_max
c      common/stick1/nstick,astick,nrmsd,ermsd
c      common/stick2/iq(ndim,nrep)
c      common/stick3/ax00(ndim,nrep),ay00(ndim,nrep),az00(ndim,nrep)
c      common/stick6/bx00(ndim),by00(ndim),bz00(ndim),armsd_min0
!      common/trackn/n_tem(100)
      dimension eshort12_a(nrep)
      common/lattice/m_latt,latt1,latt2
      common/nwmax1/nvecr
      common/iter/n_run
      common/aminoacid/sequ(ndim)
      character*3 sequ
c      common/stick7/itemp0,icycle,icycle0
      common/mloopf/mloop
      common/res2/er14,er15,er16,er17
!      common/svm1/Mcon(15,ndim),Kcon(15,ndim,100),awei(15,ndim,ndim)
!      common/svm2/acc_cut,dist_svm2(15),nk_svm(3),ik_svm(3,15)
      common/svm3/er21,er22,er23
      common/svm4/acc_cutB,acc_cutG
      
cccc  RMSD:
      double precision r_1(3,ndim),r_2(3,ndim),r_3(3,ndim),w(ndim)
      double precision u(3,3),t(3),rms,drms !armsd is real
      character(len=40) :: dateBuf
      data w /ndim*1.0/
ccc   
      real :: timeA,TimeB,start_datE,end_datR
ccccccccccccccccccccccccc common input files cccccccccccccccccccccccccccccc
      open(unit=1,file='contact.comm', status='old') !cutoff of contact predi.
      open(unit=2,file='profile3.comm',status='old') !envir
      open(unit=3,file='quasi3.comm',  status='old') !SG-potential, format
      open(unit=4,file='sidechain.comm',status='old') !for Sc position.
      open(unit=5,file='r13.comm',    status='old') !E_short of (i,i+2)
      open(unit=12,file='r14.comm',   status='old') !E_short of (i,i+3)
      open(unit=7,file='r14h.comm',   status='old') !stress helical-stru.
      open(unit=8,file='r14e.comm',   status='old') !stress extended-stru.
      open(unit=9,file='r15.comm',   status='old') !E_short(i,i+4)
      open(unit=10,file='r15h.comm',   status='old')
      open(unit=11,file='r15e.comm',   status='old')

ccccccccccccccc sequence specified input files cccccccccccccccc///////////
      open(unit=14,file='seq.dat',     status='old') !note format
      open(unit=15,file='par.dat',   status='old') !SG-pair potential
      open(unit=16,file='comb.dat',   status='old') !SG-contact restraints
      open(unit=17,file='dist.dat',    status='old') !CA-distant restraints
      open(unit=18,file='combCA.dat',status='old') !CA-contact restraints
      open(unit=28,file='comb8CA.dat',status='old') !CA-contact at 8A
      open(unit=21,file='rmsinp',status='old') !chain length for Lch
      open(unit=22,file='distL.dat',status='old') !long range CA dist restraint
      open(unit=25,file='pair3.dat',status='unknown') !orientation specific
      open(unit=26,file='pair1.dat',status='unknown') !no orientation
      open(unit=27,file='exp.dat',status='unknown')
      open(unit=24,file='init.dat',status='old')
      
      open(unit=19, file='in.dd',     status='old')
      open(unit=20, file='out.d',     status='unknown')
ccccccccccccccccfor E-t cccccccccccccccccccccccccccccccccccccccccccccccccc
c     open(unit=91, file='swepa.d',    status='unknown') !$$
c     open(unit=92, file='swepb.d',    status='unknown') !$$
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     30-40 -->rep
cccccccccccccccc sequence-based contact predictions --------------->
      open(unit=50,file='conf3.comm',status='unknown')
      open(unit=51,file='svmseqca6.dat',status='unknown') !1
      open(unit=52,file='svmseqca7.dat',status='unknown') !2
      open(unit=53,file='svmseqca8.dat',status='unknown') !3

      open(unit=54,file='svmseqcb6.dat',status='unknown') !4
      open(unit=55,file='svmseqcb7.dat',status='unknown') !5
      open(unit=56,file='svmseqcb8.dat',status='unknown') !6

      open(unit=57,file='svmseqsg6.dat',status='unknown') !7
      open(unit=58,file='svmseqsg7.dat',status='unknown') !8
      open(unit=59,file='svmseqsg8.dat',status='unknown') !9
      
      open(unit=60,file='svmcon.dat',status='unknown') !10, beta
      open(unit=61,file='betacon.dat',status='unknown') !11, CA
      open(unit=62,file='psicov.dat',status='unknown') !12, beta
      open(unit=63,file='spcon.dat',status='unknown') !13, beta
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      read(19,*) hour_max
      read(19,*) random,ncycle,phot,N_rep,n_run
      read(19,*) h2,h3s,h3d,h4s,h4d,h5s,h5d,h6,hend
      read(19,*) atemp2,atemp1,exc,exc1,exc2,Mend,defo_angle
      read(19,*) d_xyz0,angle0,L_cut,teco
      
      read(19,*) switch,i_thr0,ssp !i_thr0: 0->consensus; 1->top-1; 2->2th
      
      read(19,*) eh5a,Cr2a,acut_bb,acut_cc,acut_vv,acut_hh !alpha-type HB
      read(19,*) eh5b,Cr2b,bcut_bb,bcut_cc,bcut_vv,bcut_hh !alpha-type HB
      
      read(19,*) eh1,eh2,eh3,eh4 !for EHBs (6)
      read(19,*) eh1a,eh1b
      read(19,*) es2,es3,es4,es5,es6     !for ESHORTs (6)
      read(19,*) es3a,es3b,es3c
      read(19,*) en1,en2,en3    !for ensemble (3)
      read(19,*) chuan
      read(19,*) aTs2,aTs1    !SQ2
      read(19,*) aTTs2,aTTs1  !ARCSH
      read(19,*) nana           !number of columns in exp.dat
      read(19,*) nstick,astick,nrmsd,ermsd !whether stick to templates
      read(19,*) m_latt         !how many vectors
      
***   
       call  CPU_TIME(start_datE )
!       write(20,*)'starting time: ',fdate() !pgf77 has problem on fdate()
       write(20,*)'starting time: ', start_datE   !date(dateBuf) 
!      write(*,*)'starting time: ',fdate()!
!       start_datE = cpu_time( )
       write(20,*)'starting time:' ,start_ datE
***   
      
      call set_common           !set common parameters
      call read_seq             !read seq(i),sec(i) from 'seq.dat'
      call read_centro          !read eoinp(i,dis) from 'centro.comm'
      call read_profile         !read envir(ia,im,ip,i,j) from 'profile3.comm'
      call read_E13             !read 1-3 short-range E from 'r13.comm'
      call read_E14             !read 1-4 potential from 'r14*.dat'
      call read_E15             !read 1-5 potential from 'r15*.dat'
      call read_quarsi3         !read 'quarsi3.comm'
      call read_par             !read 'par.dat', pair-wise potential
      call read_concut          !read cut-off for contact prediction
      call read_contactrestrain !read contact restrains from 'comb.dat'
      call read_distantrestrain !read distant restrains from 'dist.dat'
      call read_longdistantrestrain !read long distant restrain from 'distL.dat'
      call read_exp             !read slovent expose prediction
      call read_CAcontact       !read CAcontact restrains from 'combCA.dat'
      call read_CAcontact8      !read CAcontact restrains from 'comb8CA.dat'
      call read_seqcontact      !read seq_based_contact 'svmseqca6.dat'
      call reset_temperature    !reset temperature according N_rest
      call set_temperature      !set temperature for different replic
      call set_EHB              !set structure-specitic H-bond, EHBIJ(i,j)
      
      call prepare_vectors      !prepare all possible bond-vectors
      call prepare_neighbors    !define goodc(i,j), angle(i,j), prod(i,j)
      call prepare_beta         !define C_beta, C_group, and hydrogen-bond
! !$acc update device(hbx,hby,hbz,gx,gy,
! !acc& gz,frga,frgb,McomCA8,KcomCA8,
! !acc& aweigCA8,area,nopp,noaa,nomm)
      call prepare_frg          !compute the secondary fragment biases

      call get_acorder          !calculate contact order
      call write_parameter      !print out initial parameters

      call set_move_retio       !set movement percentage
      call prepare_move2        !2-bond move, num2=26784
      
      n_repf=5                  !number of output replicas
      if(izscore.eq.1)then      !easy target
         n_repf=8               !number of output replicas
      endif
      if(switch.gt.1)then
         n_repf=3               !number of output replicas
         call template_initial  !initial model from templates
      endif
      
ccccccccccccccc trajectory files cccccccccccccccccccccccccccccccccccccc
      if(n_repf.gt.N_rep)n_repf=N_rep
      do i=1,n_repf
         if(i.lt.10)then
            fn=char(48+i)
            open(unit=30+i,file='rep'//fn//'.tra',status='unknown')
         else
            fn=char(48+(i-10))
            open(unit=30+i,file='rep1'//fn//'.tra',status='unknown')
         endif
      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if(switch.gt.1)then
         call template_simulation !simulate the templates
         stop
      endif

      do i=1,100
         bNSa(i)=0              !aceptance for swep
         bNSt(i)=0

         bNa(i)=0               !acceptance for move2,3,4,5,6,7
         bNt(i)=0

         bNNa(i)=0              !acceptance for different temperature.
         bNNt(i)=0

         N_sum(i)=0
         energ_sum(i)=0         !<E_tot>
         energ_sum2(i)=0        !<E_tot^2>
         eshort12_a(i)=0
      enddo
      i_tr=0                    !order number of output trajectory
      E_min=10000
      mcycle=0
      
      i_run=0
 15   i_run=i_run+1
c     call random_initial       !produce initial structure randomly.
      call read_initial         !read initial (x,y,z) from 'init.dat'

      do i=1,n_rep
         n_tem(i)=i             !i'th replica from n_tem(i)'th template
      enddo
!$acc update device(n_tem)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc       The main cycle start from here !                         ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! !$OMP parallel 
      do 1111 icycle=1,ncycle
         do 2222 itemp=1,N_rep  !iterate for all the replicas
            atemp=aT_rep(itemp)	!current temperature
            aTs=aTs_rep(itemp)
            aTTs=aTTs_rep(itemp)
            call set_current    !get current (x,y,z,ica)
            call initial_move   !update center, axis, energy
ccc
!!$OMP parallel do 
            do 3333 iphot=1,phot !N_swap, iterate at fixed temperature
               do 4444 i_lch=1,Lch
                  fff=aranzy(nozy)
                  if(fff.le.bh2)then
                     call move2
                  elseif(fff.le.bh3s)then
                     call move3s
                  elseif(fff.le.bh3d)then
                     call move3d
                  elseif(fff.le.bh4s)then
                     call move4s
                  elseif(fff.le.bh4d)then
                     call move4d
                  elseif(fff.le.bh5s)then
                     call move5s
                  elseif(fff.le.bh5d)then
                     call move5d
                  elseif(fff.le.bh6)then
                     call move6
c                     call move7a
c                     call move7b
c                     call move9
                  elseif(fff.le.bhendn)then
                     call move_n_end
                  else
                     call move_c_end
                  endif
                  call CPU_TIME(timeA)
                  atime=timeA/3600.0
                  !atime=second()/3600.0 !pgf77:etime(tarray),real tarray(2)
                  if(atime.gt.hour_max)goto 901
 4444          continue
 3333       continue
!!$OMP end parallel do                   
ccc   
ccccccccccrecord energy and (x,y,z) cccccccccccccccccccc
            E_rep(itemp)=energy_tot() !whole energy
c            energy_total=
c     $            eh1*EHB1       !general soft-core energy for Ca-SC
c     $           +eh1a*EHB1a    !general soft-core energy for CA-CA 
c     $           +eh1b*EHB1b    !general soft-core energy for SC-SC 
c     $           +eh1c*EHB1c    !pair-wise potential for SC-SC 
c     $           +eh2*EHB2      !soft-core energy for SC-SC (quarsi3)
c     $           +eh3*EHB3      !coupling of secondary structure and pairwise
c     $           +eh4*EHB4      !enhanced quarsi3
c     $           +eh5a*EHB5a    !H-bond (alpha)
c     $           +eh5b*EHB5b    !H-bond (beta)
c     $           +es2*ESHORT2   !bury potential for SC
c     $           +er1*ESHORT3   !distmap
c     $           +er3*ESHORT4   !contact restrain
c     $           +er4*ESHORT4a  !deviation of contact restrain
c     $           +er5*ESHORT9   !CAcontact restrain
c     $           +er6*ESHORT10   !longrange CA-dist restrain
c     $           +er7*ESHORT11   !derivation to template
c     $           +astick*ESHORT12   !derivation to ax00
c     $           +es3*ESHORT5   !general bias to protein-like structure
c     $           +es3a*ESHORT5a  !panality on crumpling
c     $           +es3b*ESHORT5b  !bias to predicted alpha/beta structures
c     $           +es3c*ESHORT5c  !bias to possible alpha/beta structures
c     $           +es4*ESHORT6   !R13
c     $           +es5*ESHORT7   !R14
c     $           +es6*ESHORT8   !R15
c     $           +en1*eprofo    !E_environment
c     $           +en2*E_cord    !Contact order
c     $           +en3*E_cnum    !Contact number
c     parameters with input: eh1,eh1a,eh1b,eh2,eh3,eh4,eh5a,eh5b,
c                            es2,es3,es3a,es3b,es3c,es4,es5,
c                            en1,en2,en3
c     parameters mandtoried: eh1c,er1,er2,er3,er4,er5,er6,er7
c            write(*,*)E_rep(itemp),energy_total
c            write(*,*)icycle,itemp,E_rep(itemp),eshort12,'=='
            if(E_rep(itemp).lt.E_min) E_min=E_rep(itemp)
            do i=1,Lch
               xrep(i,itemp)=x(i)
               yrep(i,itemp)=y(i)
               zrep(i,itemp)=z(i)
            enddo
            eshort12_a(itemp)=eshort12_a(itemp)+eshort12
 2222    continue
         
cccccccccccccccccc print out 'swep.d' cccccccccccccccccccccccccccc
c         write(91,91)icycle,(E_rep(i),i=1,20)       !$$
c         write(92,91)icycle,(E_rep(i),i=21,N_rep)   !$$
c 91      format(i10,21f9.1)

ccccccccccccccccc<RMSD>, <E> cccccccccccccccccccccccccc
         do i=1,N_rep
            energ_sum(i)=energ_sum(i)+E_rep(i)
            energ_sum2(i)=energ_sum2(i)+E_rep(i)*E_rep(i)
            N_sum(i)=N_sum(i)+1
         enddo
         
ccccccccccccccccccccc snapshots of E(1), E(N_rep) ccccccccccccc
         if(icycle.eq.icycle/1*1)then
            i_tr=i_tr+1
            do k=1,n_repf
               write(30+k,401)Lch,E_rep(k),i_tr,icycle
               do i=1,Lch
                  abx=xrep(i,k)*0.87
                  aby=yrep(i,k)*0.87
                  abz=zrep(i,k)*0.87
                  write(30+k,402)abx,aby,abz
               enddo
            enddo
         endif
 401     format(i8,1x,f10.1,2i8)
 402     format(f10.3,1x,f10.3,1x,f10.3)

         call count_restrains   !count number of satisfied restraints
         
ccccccccccccccccc swap replicas cccccccccccccccccccccccccccccccc
         if(icycle.eq.icycle/2*2)then
            do i=1,N_rep-1,2    !swap odd replicas
               call swap(i,i+1)
            enddo
         else
            do i=2,N_rep-1,2
               call swap(i,i+1) !swap even replicas
            enddo
          endif
          mcycle=mcycle+1
 1111 continue
! !$OMP end parallel           
      if(i_run.lt.n_run)goto 15
 901  continue
c--------------------------Main cycle ended here!!---------------

      call test_neighbor        !check all the neighboring residues
      call test_overlap         !test the overlap of C_a and C_b

      energy_tot_tmp=energy_tot()
      write(20,*)'E_final=',energy_tot_tmp

      write(20,*)
      write(20,*)'<s_comb>=',s_comb/float(N_resc)
      write(20,*)'<t_comb>=',t_comb/float(N_resc)
      write(20,*)'s_comb/t_comb=',float(s_comb)/(t_comb+0.001)
      write(20,*)
      write(20,*)'<s_dist>=',s_dist/float(N_resc)
      write(20,*)'<t_dist>=',t_dist/float(N_resc)
      write(20,*)'s_dist/t_dist=',float(s_dist)/(t_dist+0.001)
      write(20,*)
      write(20,*)'<s_distL>=',s_distL/float(N_resc)
      write(20,*)'<t_distL>=',t_distL/float(N_resc)
      write(20,*)'s_distL/t_distL=',float(s_distL)/(t_distL+0.001)
      write(20,*)
      write(20,*)'<s_combCA>=',s_combCA/float(N_resc)
      write(20,*)'<t_combCA>=',t_combCA/float(N_resc)
      write(20,*)'s_combCA/t_combCA=',float(s_combCA)/(t_combCA+0.001)
      write(20,*)
      write(20,*)'<s_combCA8>=',s_combCA8/float(N_resc)
      write(20,*)'<t_combCA8>=',t_combCA8/float(N_resc)
      write(20,*)'s_combCA8/t_combCA8=',
     &     float(s_combCA8)/(t_combCA8+0.001)
      write(20,*)
      
cccccc output 'stick.pdb' ccccccccccccccccccccccccccccccc
c     sticki.pdb is the conformation closest to i'th structure in init.pdb
      if(nrmsd.eq.1)then
        write(20,*)' ------- RMSD to templates ------------'
        write(20,*)'rmsd_min0=',armsd_min0*0.87
        write(20,*)'icycle0=',icycle0
        write(20,*)'itemp0=',itemp0
        
        open(unit=70,file='stick.pdb',status='unknown')
        do i=1,Lch
          write(70,1037)i,sequ(i),i,bx00(i)*0.87,by00(i)*0.87,
     &      bz00(i)*0.87
        enddo
        write(70,*)'TER'
        close(70)
      endif
 1037 format('ATOM  ',i5,'  CA  ',A3,I6,4X,3F8.3)
      
cccccccccccccccccccccccc Na/Nt cccccccccccccccccccccccccc
      WRITE(20,*)
      WRITE(20,*) 'i_move    move   Na(i)  Nt(i)   Na(i)/Nt(i)'
      do i=2,15
         if(bNt(i).gt.1)then
            write(20,5004) i,mname(i),bNa(i),bNt(i),bNa(i)/bNt(i)
         else
            write(20,5004) i,mname(i),bNa(i),bNt(i)
         endif
      enddo
 5004 format(I4,A9,2f15.1,f11.6)
      
ccccccccccccccccccccccccccE_final, NSa/NSt ccccccccccccccccccccccc
      WRITE(20,*)
      WRITE(20,*)'-------------- E_final, Na_swap/Nt_swap ---------'
      WRITE(20,*) 'i  T(i) final_E(i)  Nsa(i)  Nst(i)  Nsa(i)/Nst(i)'
      do i=1, n_rep
         if(bNSt(i).gt.1)then
            write(20,5005) i,aT_rep(i),E_rep(i),
     $           bNSa(i),bNSt(i),bNSa(i)/bNSt(i)
         else
            write(20,5005) i,aT_rep(i),E_rep(i),bNSa(i),bNSt(i)
         endif
      enddo
 5005 format(I4,f7.2,f8.1,2f15.1,f11.6)

ccccccccccccccccccccccc <E>, NNa/NNt ccccccccccccccccccccccccccc
      WRITE(20,*)
      WRITE(20,*)'------------ <energy>, Na(i)/Nt(i) ----------------'
      write(20,*)'i_rep  T   <E>   NNa(i_temp)  NNt(i_temp)  Na/Nt
     &     <ESHORT12>'
      do i=1,N_rep
         energ_sum(i)=energ_sum(i)/N_sum(i)
         energ_sum2(i)=energ_sum2(i)/N_sum(i)
         if(bNNt(i).gt.1)then
            cheat=(energ_sum2(i)-energ_sum(i)**2)/(aT_rep(i)**2)
            write(20,5006)i,aT_rep(i),energ_sum(i),cheat,bNNa(i)
     $           ,bNNt(i),bNNa(i)/bNNt(i),eshort12_a(i)/N_sum(i)
         else
            write(20,5006)i,aT_rep(i),energ_sum(i),cheat,bNNa(i)
     $           ,bNNt(i)
         endif
      enddo
 5006 format(I4,f7.2,f8.1,f15.3,2f12.1,f11.6,6f8.1)
      write(20,*)'E_min=',E_min

      write(20,*)
      write(20,*)'ncycle_max=',ncycle*n_run
      write(20,*)'ncycle_real=',mcycle
      write(20,*)
      write(20,*)'hour_max=',hour_max
      write(20,*)'hour_real=',atime
      write(20,*)
!      write(20,*)'ending time: ',fdate()
      call  CPU_TIME(end_datR )
      write(20,*)'ending time: ',end_datR !date(dateBuf) 

      write(*,*)
      write(*,*)'ncycle_max=',ncycle*n_run
      write(*,*)'ncycle_real=',mcycle
      write(*,*)
      write(*,*)'hour_max=',hour_max
      write(*,*)'hour_real=',atime
      write(*,*)
!      write(*,*)'ending time: ',fdate()
      call  CPU_TIME(end_datR )
      write(20,*)'ending time: ',end_datR  !fdate() !date(dateBuf) 
     
      STOP
      END
ccccccc=======================================================cccccccccc
cc          The main program ended!
ccccccc=======================================================cccccccccc
