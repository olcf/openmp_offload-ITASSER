
*****************************************************************************
*     This program is to generate protein structural decoys by Monte Carlo  *
*     simulations under an on-and-off-lattice CAS model. It is illigal to   *
*     distribute any part of this code without writing permission from the  *
*     author. Please address comments/bug-reports to: zhng@umich.edu        *
*****************************************************************************
*
*     This program should be compiled by 
*     'gfortran -static -O3 -ffast-math -lm -o cas cas.f'
*
*     Last update by yzhang on March 18 2013
*
!!      common/shape/amx,amy,amz,afs(ndim),afsn(ndim) ! !common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1 !!common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
c        1         2         3         4         5         6         7 !
c 3456789012345678901234567890123456789012345678901234567890123456789012345678
       module params
       parameter ndim = 1500
       parameter nrep = 100
       parameter nvec = 416
       end module params


      module trackn
          use params
          integer n_tem(100)
!$acc declare create(n_tem)
      end module trackn
 

        
      module chainm
         use params
         real mv(ndim)
!$acc declare create(mv)
      end module chainm
             
      module eigen
         real AA(3,3), EE(3), HH(3,3)
!$acc declare create(AA,EE,HH)         
      end module


      module ehbc
          use params
          real envir(0:15,0:15,0:15,0:19,4),en1
!$acc declare create(envir,en1)   
      end module ehbc

      module backup1
         use params
         integer nopp(ndim),nomm(ndim),noaa(ndim),nhbnn(ndim)
!$acc declare create(nopp,nomm,noaa,nhbnn)
      end module backup1


      module shape 
         use params
         real amx,amy,amz,afs(ndim),afsn(ndim)
!$acc declare create(amx,amy,amz,afs,afsn)
      end module shape

      module stick
	  use params
          integer nstick,astick,nrmsd,ermsd
          integer iq(ndim,nrep)
          real ax00(ndim,nrep),ay00(ndim,nrep),az00(ndim,nrep)
          real bx00(ndim),by00(ndim),bz00(ndim),armsd_min0
          integer itemp0,icycle,icycle0
!$acc declare create (ax00,ay00,az00,bx00,by00,bz00)
      end module

      module short
         use params 
          integer  IBIN(-300:300)
          real asr(ndim,-12:12),csr(ndim,2)
          real codevsum, didevsum 
!$acc declare create(IBIN,asr,csr)
      end module short 

      module ENERGY
         use params
         real EH5,ES3,ES3a,ES3b,EH1,EH1a,EH4,EHBIJ(ndim,ndim)
!$acc declare create(EH5,ES3,ES3a,ES3b,EH1,EH1a,EH4,EHBIJ)
      end module ENERGY
      
      module svm1
         use params
            integer Mcon(15,ndim),Kcon(15,ndim,100)
            real awei(15,ndim,ndim) 
!$acc declare create(Mcon,Kcon,awei)
      end module svm1
     
      module svm2
        use params
            integer acc_cut, nk_svm(3)
            integer ik_svm(3,15)   
            real dist_svm2(15)   
!$acc declare create(acc_cut,nk_svm,ik_svm,dist_svm2)
      end module svm2

      module RES
        use params 
        integer ER3,er5,er6,er7,Mcom(ndim)
        real Kcom(ndim,100)
!$acc declare create(ER3,er5,er6,er7,Mcom(ndim),Kcom(ndim,100))
      end module RES 

      module hba
             integer eh5a,Cr2a,acut_bb,acut_cc,acut_vv,acut_hh
!$acc declare create(eh5a,Cr2a,acut_bb,acut_cc,acut_vv,acut_hh)
      end module hba

      module lengths
             integer Lch, Lch1, Lch2
!$acc declare create(Lch,Lch1,Lch2)
      end module lengths

      module short1
             integer  JBIN(0:500) 
             real  acops(1500,16)
             real  bsr(1500,16)
!$acc declare create (JBIN,bsr,acops)
      end module short1
  
      module echain1
             real ex(1500),ey(1500),ez(1500)
!$acc declare create(ex,ey,ez)
      end module echain1
      
      module shortcom
              integer eh3,es4,es5,es6,es7,es7a,es7b,es7c
!$acc declare create(eh3,es4,es5,es6,es7,es7a,es7b,es7c)
      end module shortcom

      module hbb
             integer eh5b,Cr2b,bcut_bb,bcut_cc,bcut_vv,bcut_hh
!$acc declare create (eh5b,Cr2b,bcut_bb,bcut_cc,bcut_vv,bcut_hh)
      end module hbb

      module chain1
             integer  x(1500),y(1500),z(1500)
             integer ica(0:1500)
!$acc declare create(x,y,z,ica)
      end module chain1
      
      module seqe
             integer  seq(1500), sec(1500)
!$acc declare create(seq,sec)
      end module seqe
      
      module hb
              real hbx(416,416), hby(416,416),hbz(416,416)
!$acc declare create(hbx,hby,hbz)
      end module hb
      
      module bisec
             real cax(416,416), cay(416,416), caz(416,416)
!$acc declare create(cax,cay,caz)
      end module bisec
      
      module sizea
             real ala(0:19,0:19),alm(0:19,0:19),alp(0:19,0:19)
!$acc declare create(ala,alm,alp)
      end module sizea
      
      module pair 
             real apa(1500,1500),app(1500,1500),apm(1500,1500)
!$acc declare create(apa,app,apm)
      end module pair
      
      module size
             real arla(0:19,0:19), arlm(0:19,0:19),arlp(0:19,0:19)
!$acc declare create(arla,arlm,arlp)
      end module size
      
      module ehbenergy
             integer EHB1,EHB1b,EHB1c,EHB2,EHB3,EHB5,EHB6  !EHB4,EHB1a
!$acc declare create(EHB1,EHB1b,EHB1c,EHB2,EHB3,EHB5,EHB6) !EHB4,EHB1a
      end module ehbenergy
      
      module ehbenergy1
             integer EHB5a,EHB5b
!$acc declare create(EHB5a,EHB5b)
      end module ehbenergy1
         
      module tempArrays
              use params
              dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
              dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

              dimension ex_o(ndim),ey_o(ndim),ez_o(ndim) !CA
              dimension egx_o(ndim),egy_o(ndim),egz_o(ndim) !SG
              dimension ecx_o(ndim),ecy_o(ndim),ecz_o(ndim) !cc
              dimension ebx_o(ndim),eby_o(ndim),ebz_o(ndim) !Hb
              dimension etx_o(ndim),ety_o(ndim),etz_o(ndim) !CB
              dimension ex_n(ndim),ey_n(ndim),ez_n(ndim) !CA
              dimension egx_n(ndim),egy_n(ndim),egz_n(ndim) !SG
              dimension ecx_n(ndim),ecy_n(ndim),ecz_n(ndim) !cc
              dimension ebx_n(ndim),eby_n(ndim),ebz_n(ndim) !Hb
              dimension etx_n(ndim),ety_n(ndim),etz_n(ndim) !CB     
!$acc declare create(ox,nx,ex_o,egx_o,ecx_o,ebx_o,etx_o,
!$acc&                ex_n,egx_n,ecx_n,ebx_n,etx_n,oy,oz,oo,
!$acc&                ny,nz,nn,ey_o,ez_o,egy_o,egz_o,
!$acc&                ecy_o,ecz_o,eby_o,ebz_o,ety_o,etz_o,
!$acc&                ety_o,etz_o,ey_n,ez_n,egy_n,egz_n,
!$acc&                ecy_n,ecz_n,eby_n,ebz_n,ety_n,etz_n)
      end module tempArrays
 
      module par
             real apar(1500,1500)
!$acc declare create(apar)
      end module par
     
      module expose
         use params
	 integer :: mp(20,ndim),area(ndim)
!$acc declare create(mp,area)
      end module expose

      module CAcontact
c Merge CAcontact, CA8
	use params
        integer McomCA(ndim),KcomCA(ndim,100)
	integer McomCA8(ndim),KcomCA8(ndim,100)
	real aweigCA(ndim,ndim),aweigCA8(ndim,ndim)
!$acc declare create(McomCA,KcomCA,McomCA8,KcomCA8,aweigCA,aweigCA8)
      end module CAcontact
 
      module sg
             real hx(416,416,0:19),hy(416,416,0:19),hz(416,416,0:19)
             real gx(416,416,0:19),gy(416,416,0:19),gz(416,416,0:19)
!$acc declare create(gx,gy,gz,hx,hy,hz)
      end module sg
      
      module icgg
             integer EH6,icg(1500)
!$acc declare create(icg,EH6)
      end module icgg
      
      module envir1
             real nop(1500),nom(1500),noa(1500),nhbn(1500)
!$acc declare create(nop,nom,noa,nhbn)
      end module envir1
      
      module distres
             integer er4,es3c
!$acc declare create(er4,es3c)
      end module distres
      
      module  order
              integer acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
!$acc declare create(acorder,en2,sumct,sumcto,icnt,icnto,dord,en3) 
      end module order 
      
      module echain2
              real egx(1500),egy(1500),egz(1500)
!$acc declare create (egx,egy,egz)
      end module echain2
      
      module echain4
             real ecx(1500),ecy(1500),ecz(1500)
!$acc declare create (ecx,ecy,ecz)
      end module echain4
      
      module echain5
             real ebx(1500),eby(1500),ebz(1500)
!$acc declare create (ebx,eby,ebz)
      end module echain5
     
      module fr
         use params
         real frga(ndim),frgb(ndim)
!$acc declare create (frga,frgb)
      end module fr

      module freg
          use params
          real aweig(ndim,ndim)
!$acc declare create (aweig)
      end module freg

      module echain6
             real etx(1500),ety(1500),etz(1500)
!$acc declare create (etx,ety,etz)
      end module echain6
      
      module concutt
          real concut(0:19,0:19),concut2(0:19,0:19)
          integer concut_sc
!$acc declare create(concut,concut2,concut_sc)
      end module concutt
      
      module pair1
          integer eh2,eh1b,eh1c
!$acc declare create(eh2,eh1b,eh1c)
      end module pair1

      module longdist
           use params
           integer MdisL(ndim),kdisL(ndim,500),distL(ndim,500)
!$acc declare create(MdisL(:ndim),kdisL(:ndim,:500),distL(:ndim,:500))
      end module longdist

!      COMMON/RCN/Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
      module RCN
        use params
        integer Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
!$acc declare create(Mdis(:ndim),kdis(:ndim,:100),dist(:ndim,:100),
!$acc&   dev(:ndim,:100))
      end module RCN

      module one
        use lengths
             integer acrit,es2,es1
             real eoinp(0:19,0:100),eonekd(0:19),contt
             ! real, PRIVATE :: eonekd(0:19)
!$acc declare create(acrit,contt,eoinp,es2,es1,eonekd)
cccccccccccccccccccccc read centrosymmetric potential
cccccccccccccccccccc
c     eonekd(A) controls centrosymmetric potential of C_g.
      contains 
      subroutine read_centro
!      use lengths
!      use one
      implicit integer(i-z)
      parameter(nvec=416)
      character*3 NAME
!      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1
!      common/lengths/Lch,Lch1,Lch2
      common/hopp/eonehw(0:19)

c     read hydrophobic potential for Sg, positive for hydrophobic
c     residue--->
      data eonekd /-0.4, 1.8, -0.8, 2.5, 4.2, -0.7, 4.5,
     &     -1.6, 1.9,  -3.5, -3.5, 3.8,
     &     -3.9, -3.5, -3.5, -4.5,
     &     -3.2, 2.8, -1.3, -0.9/
c            ^          ^     ^     !contradict with 'centro.comm'
c     read hydrophilic potential for Sg, positive for hydrophilic
c     residue--->
      data eonehw /0.0, -0.5, 0.3, -1.0, -1.5, -0.4, -1.8,
     &     0.0, -1.3, 3.0, 0.2, -1.8,
     &     3.0, 3.0, 0.2, 3.0,
     &     -0.5, -2.5, -2.3, -3.4/

c     expected gyration radius:
      acrit=2.2*exp(0.38*alog(float(Lch)))/0.87 !gyrat-radius~2.2*l^0.38
*     Defination of gyration-radius: acrit=sqrt(<(r-r0)^2>)

c^^^^^^^^^^^^^^^^^ read centersymmetric potential finished ^^^^^^^^^^^
      end subroutine read_centro
! !$acc declare create(acrit,contt,eoinp,es2,es1,eonekd)   
      end module one 
      
      program TASSER
      use params
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
      common/hopp/eonehw(0:19)
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
      common/backup2/eprofo,eprofn,energ
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
       real :: timeA,TimeB
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
!      write(20,*)'starting time: ',fdate() !pgf77 has problem on fdate()
       write(20,*)'starting time: ',date(dateBuf) 
!      write(*,*)'starting time: ',fdate()
       write(20,*)'starting time: ',date(dateBuf) 
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
      do 1111 icycle=1,ncycle
         do 2222 itemp=1,N_rep  !iterate for all the replicas
            atemp=aT_rep(itemp)	!current temperature
            aTs=aTs_rep(itemp)
            aTTs=aTTs_rep(itemp)
            call set_current    !get current (x,y,z,ica)
            call initial_move   !update center, axis, energy
ccc
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
      write(20,*)'ending time: ',date(dateBuf) 

      write(*,*)
      write(*,*)'ncycle_max=',ncycle*n_run
      write(*,*)'ncycle_real=',mcycle
      write(*,*)
      write(*,*)'hour_max=',hour_max
      write(*,*)'hour_real=',atime
      write(*,*)
!      write(*,*)'ending time: ',fdate()
      write(20,*)'ending time: ',date(dateBuf) 
     
      STOP
      END
ccccccc=======================================================cccccccccc
cc          The main program ended!
ccccccc=======================================================cccccccccc
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
      subroutine read_seq
	use params
      use lengths
      use seqe
      use icgg
      implicit integer(i-z)
!      parameter(ndim=1500)
!                parameter(nvec=416)
      character*3 aa(-1:20), NAME,sequ
!      common/seqe/seq(ndim),sec(ndim)
!      common/lengths/Lch,Lch1,Lch2
!      common/icgg/ icg(ndim), EH6  
      common/aminoacid/sequ(ndim)

      data aa/ 'BCK','GLY','ALA','SER','CYS','VAL','THR','ILE',
c                -1    0     1     2     3     4     5     6
     &     'PRO','MET','ASP','ASN','LEU',
c            7     8     9    10    11
     &     'LYS','GLU','GLN','ARG',
c           12    13    14    15
     &     'HIS','PHE','TYR','TRP','CYX'/
c           16    17    18    19    20

      do 121 i=1,Lch
         read(14,707) k,NAME,SEC(I),tmp
         do j=0,19
            if(NAME.eq.aa(j)) then
               SEQ(i)=j
               icg(i)=0
               sequ(i)=name
               if(NAME.eq.'ASP'.or.NAME.eq.'GLU') icg(i)=-1
               if(NAME.eq.'LYS'.or.NAME.eq.'ARG') icg(i)= 1	
               go to 121
            endif
         enddo
         SEQ(i)=0
         icg(i)=0
         sequ(i)='GLY'
 121  continue
 707  format(i5,3x,a3,2i5)
      close(14)

c^^^^^^^^^^^^^^^^^^^^^ read sequence finished ^^^^^^^^^^^^^^^^^^
      return
      end
      subroutine read_profile
	use params
      use one
      use ehbc
      implicit integer(i-z)
!      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
!      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1

      do kkk=1,4
      do i=0,19
      do im=0,15
      do ip=0,15
      do ia=0,15
         envir(im,ip,ia,i,kkk)=2.0
      end do
      end do
      end do
      end do
      enddo

c     PROFILE3 potential =envir(#of antiparalel contacts,
c     #of orthogonal, # of parallel, aminoacid's type)
c     ia,im, ip - taken modulo 2, i.e. 0-1 contact, 2-3,...
c     profile3.comm is a score table, i.e. envir(envir_class,A)
c     here environment class is number of contacts on residue A.

      do i=0,19                 !from column
         read(2,*)
         do im=0,8
         do ia=0,8
            read(2,*)(envir(ia,im,ip,i,3),ip=0,8) !this is used
         end do
         read(2,*)
         end do
         read(2,*)
      enddo

c      do i=0,19                 !from column
c         read(2,*)
c         do im=0,8
c         do ia=0,8
c            read(2,*)(envir(ia,im,ip,i,1),ip=0,8)
c         end do
c         read(2,*)
c         end do
c         read(2,*)
c      end do

c      do i=0,19                 !from column
c         read(2,*)
c         do im=0,8
c         do ia=0,8
c            read(2,*)(envir(ia,im,ip,i,2),ip=0,8)
c         end do
c         read(2,*)
c         end do
c         read(2,*)
c      end do

c      do i=0,19                 !from column
c         read(2,*)
c         do im=0,8
c         do ia=0,8
c            read(2,*)(envir(ia,im,ip,i,4),ip=0,8)
c         end do
c         read(2,*)
c         end do
c         read(2,*)
c      end do

c^^^^^^^^^^^^^^^^^^^^^^^^ read profile finished ^^^^^^^^^^^^^^^^^^^^
      return
      end
      subroutine read_E13
	use params
      use lengths
      use shortcom
      use seqe
      use short
      implicit integer(i-z)
!      parameter(ndim=1500)
!                parameter(nvec=416)
!      common/short2/codevsum,didevsum,csr(ndim,2)
!      common/lengths/Lch,Lch1,Lch2
!      common/seqe/seq(ndim),sec(ndim)
!      common/shortcom/eh3,es4,es5,es6,es7,es7a,es7b,es7c
      dimension csre(0:19,0:19,2)

c     R13 potential - two bins only (helical and expanded)
c     r2<48, E=csre(i,j,1); r2>48, E=csre(i,j,2)
      do i=0,19
         do j=0,19
            read(5,*)
            read(5,*) (csre(i,j,k),k=1,2)
         enddo
      enddo

      do i=1,Lch2
         do k=1,2
            csr(i,k)=2.0*csre(seq(i),seq(i+2),k)
         enddo
      enddo

!$acc update device(csr)

c^^^^^^^^^^^^^^^^^ read E13 finished ^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
      subroutine read_E14
	use params
      use lengths
      use seqe
      use shortcom
      use short
      implicit integer(i-z)
!      parameter(ndim=1500)
!                parameter(nvec=416)
!      COMMON/short/ IBIN(-300:300),asr(ndim,-12:12)
!      common/lengths/Lch,Lch1,Lch2
!      common/seqe/seq(ndim),sec(ndim)
      common/forpreparemove4/asrr(0:19,0:19,-12:12)
      DIMENSION asrh(0:19,0:19,-12:12),asre(0:19,0:19,-12:12)
!      common/shortcom/eh3,es4,es5,es6,es7,es7a,es7b,es7c

ccccccccc read asrr,asrh,asre------------------>
      do i=0,19                 !asrr(Ai,Bi,dis) from 'r14.comm'
         do j=0,19
            read(12,*)
            read(12,*) (asrr(i,j,k),k=-12,-5)
            read(12,*) (asrr(i,j,k),k=-4,3) !without k=4
            read(12,*) (asrr(i,j,k),k=5,12)
            do k=4,1,-1
               asrr(i,j,k)=asrr(i,j,k-1) !without k=0
            enddo
         enddo
      enddo
      do i=0,19                 !asrh(Ai,Bi,dis) from 'r14h.comm'
         do j=0,19
            read(7,*)
            read(7,*) (asrh(i,j,k),k=-12,-5)
            read(7,*) (asrh(i,j,k),k=-4,3)
            read(7,*) (asrh(i,j,k),k=5,12)
            do k=4,1,-1
               asrh(i,j,k)=asrh(i,j,k-1)
            enddo
         enddo
      enddo
      do i=0,19                 !asre(Ai,Bi,dis) from 'r14e.comm'
         do j=0,19
            read(8,*)
            read(8,*) (asre(i,j,k),k=-12,-5)
            read(8,*) (asre(i,j,k),k=-4,3)
            read(8,*) (asre(i,j,k),k=5,12)
            do k=4,1,-1
               asre(i,j,k)=asre(i,j,k-1)
            enddo
         enddo
      enddo
c^^^^^^^^^ read asrr,asrh,asre finished ^^^^^^^^^^^^^^^^^

      do i=1,Lch-3
         do k=-12,12
            asr(i,k)=asrr(seq(i+1),seq(i+2),k) !general
            if(sec(i+1).eq.2.AND.sec(i+2).eq.2.AND.sec(i+3).eq.2)then
               if(sec(i).eq.2) then !helix
c                  asr(i,k)=(asr(i,k)+asrh(seq(i+1),seq(i+2),k))/2.0
                  asr(i,k)=asrh(seq(i+1),seq(i+2),k)
               endif
            endif
            if(sec(i+1).eq.4.AND.sec(i+2).eq.4.AND.sec(i+3).eq.4)then
               if(sec(i).eq.4) then !sheet
c                  asr(i,k)=(asr(i,k)+asre(seq(i+1),seq(i+2),k))/2.0
                  asr(i,k)=asre(seq(i+1),seq(i+2),k)*1.5
               endif
            endif
         enddo
      enddo
c^^^^^^^^^^^^ asr(i,ibin(r14)) finished ^^^^^^^^^^^^^^^^^^^^^
c     r(i,i+3)=k, E=asr(i,k), 12 bins (24 bins when considering chiral)
      do i=1,300
         kk=int((sqrt(float(i))*0.87))+1
         if(kk.gt.12) kk=12
         IBIN(I) = kk           !convert lattice r^2 into real r
         IBIN(-I)=-kk
      ENDDO
      IBIN(0)=IBIN(1)
!$acc update device(IBIN,asr(:,:))
c^^^^^^^^^^^^^^^^^ read E14 finished ^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
      subroutine read_E15
	use params
      use short1
      use lengths
      use seqe
      use shortcom
      implicit integer(i-z)
!      parameter(ndim=1500)
!                parameter(nvec=416)
!      COMMON/short1/ JBIN(0:500),bsr(ndim,16),acops(ndim,16)
!      common/lengths/Lch,Lch1,Lch2
!      common/seqe/seq(ndim),sec(ndim)
      DIMENSION bsrh(0:19,0:19,16)
      dimension bsre(0:19,0:19,16)
      dimension bsrr(0:19,0:19,16)
!      common/shortcom/eh3,es4,es5,es6,es7,es7a,es7b,es7c

cccccccc read bsrr,bsrh,bsre ----------------->
      do i=0,19                 !read bsrr from 'r15.dat'
         do j=0,19
            read(9,*)
            read(9,*) (bsrr(i,j,k),k=1,8)
            read(9,*) (bsrr(i,j,k),k=9,16)
         enddo
      enddo
      do i=0,19                 !read bsrh from 'r15h.dat'
         do j=0,19
            read(10,*)
            read(10,*) (bsrh(i,j,k),k=1,8)
            read(10,*) (bsrh(i,j,k),k=9,16)
         enddo
      enddo	
      do i=0,19                 !read bsre from 'r15e.dat'
         do j=0,19
            read(11,*)
            read(11,*) (bsre(i,j,k),k=1,8)
            read(11,*) (bsre(i,j,k),k=9,16)
         enddo
      enddo	

      do i=1,Lch-4
         do k=1,16
            bsr(i,k)=bsrr(seq(i+1),seq(i+3),k)
            if(sec(i+1).eq.2.AND.sec(i+2).eq.2.AND.sec(i+3).eq.2)then
c               bsr(i,k)=(bsr(i,k)+bsrh(seq(i+1),seq(i+3),k))/2.0 !helix
               bsr(i,k)=bsrh(seq(i+1),seq(i+3),k)
            endif
            if(sec(i+1).eq.4.AND.sec(i+2).eq.4.AND.sec(i+3).eq.4)then
c               bsr(i,k)=(bsr(i,k)+bsre(seq(i+1),seq(i+3),k))/2.0 !sheet
               bsr(i,k)=bsre(seq(i+1),seq(i+3),k)*1.5
            endif
         enddo
      enddo
c^^^^^^^^^^^^^^^^^^^^ E_15(Ai,Aj,dis) prepared ^^^^^^^^^^^^^^^^^^^

c     prepare distance bin-------------------->
      do i=0,500
         kk=int((sqrt(float(i))*0.87))+1 !i, lattice-dist; kk, real distance
         if(kk.gt.16) kk=16
         JBIN(I) = kk           !jbin: real distance
      ENDDO

ccccc acops(i,jbin) to enhance the contacts between gragments cccc
      do i=1,Lch-4
         acops(i,1)=(min(bsr(i,1),0.0))/2.0 !acpos<0
         do k=2,15
            acops(i,k)=min(0.0,bsr(i,k-1)+2.0*bsr(i,k)+bsr(i,k+1)) !<0
         enddo
         acops(i,16)=(min(bsr(i,16),0.0))/2.0
      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!$acc update device(bsr,acops,JBIN)

c^^^^^^^^^^^^^^^^^ read E15 finished ^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
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
      subroutine read_par
	use params
      use lengths
      use pair1
      use par
      implicit integer(i-z)
!      parameter(ndim=1500)
!      common/lengths/Lch,Lch1,Lch2
!      common/par/apar(ndim,ndim)
      common/fractpair/fract_pair1,fract_pair3
      common/pair33/i_pair3,mk_pair3
      common/weight/chuan
!      common/pair1/eh2,eh1b,eh1c

      dimension apar1(ndim,ndim) !from par.dat
      dimension apar2(ndim,ndim) !from pair1.dat

c     read from 'par.dat'-------------->
      Nline=1000
      do i=1,Lch
         read(15,*)
         do i_line=1,Nline
            line_end=min(10*i_line,Lch) !25,50,75,100,...., ending point
            read(15,*)(apar1(i,j),j=(i_line-1)*10+1,line_end)
            if(line_end.ge.Lch) go to 1
         enddo
 1       continue
      enddo

c     read from 'pair1.dat'-------------->
      if(i_pair3.gt.0)then      ! 'pair1.dat' exist
         Nline=100              !maximum Lch is 2500 !!!
         do i=1,Lch
            read(26,*)
            do i_line=1,Nline
               line_end=min(25*i_line,Lch) !25,50,75,100,...., ending point
               read(26,*)(apar2(i,j),j=(i_line-1)*25+1,line_end)
               if(line_end.eq.Lch) go to 2
            enddo
 2          continue
            read(26,*)
         enddo
      endif
      
***   combine and rescale the pair interaction ---------------->
      apar1_max=-100
      apar2_max=-100
      do i=1,Lch
         do j=1,Lch
            if(abs(apar1(i,j)).gt.apar1_max)apar1_max=abs(apar1(i,j))
            if(abs(apar2(i,j)).gt.apar2_max)apar2_max=abs(apar2(i,j))
         enddo
      enddo
      do i=1,Lch
         do j=1,Lch
            apar(i,j)=0
            if(i_pair3.gt.0)then !'pair3.dat' and 'pair1.dat' exist
               if(apar1_max.lt.0.00001)then !par.dat is wrong
                  apar(i,j)=apar2(i,j) !'pair1.dat' only
               else
                  apar(i,j)=apar1(i,j)*apar2_max/apar1_max*
     &                 (1-fract_pair1)+apar2(i,j)*fract_pair1
               endif
            else
               apar(i,j)=apar1(i,j) !'pair1.dat' not exist, 'par.dat' only
            endif
         enddo
      enddo
      
c      write(*,*)apar1_max,apar2_max
c      do i=1,Lch
c         do j=1,Lch
c            write(*,*)i,j,apar1(i,j),apar2(i,j),apar(i,j)
c         enddo
c      enddo
c      stop

c^^^^^^^^^^^^^^^^^ pair-wise potential finished ^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
      subroutine read_concut
	use params
      use concutt
      implicit integer(i-z)
      character*3 NAME
      dimension cut(0:19,0:19),cut_dev(0:19,0:19)
!      common/concutt/concut(0:19,0:19),concut2(0:19,0:19),concut_sc
      common/zscore/izscore

      rewind(1)
      read(1,*)
      do i=0,19
         read(1,*)NAME,(cut(i,j),j=0,19)
      enddo
      read(1,*)
      do i=0,19
         read(1,*)NAME,(cut_dev(i,j),j=0,19)
      enddo

      do i=0,19
         do j=0,19
            if(izscore.eq.1)then
               concut(i,j)=7
            elseif(izscore.eq.2)then
               concut(i,j)=7.5
            else
               concut(i,j)=cut(i,j)+cut_dev(i,j)*2.5 !real cut-off
            endif
            concut(i,j)=concut(i,j)/0.87 !cut-off on lattice
            concut2(i,j)=concut(i,j)**2 !cut-off squared on lattice
         enddo
      enddo
!$acc update device (concut,concut2)
      return
      end
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
      subroutine read_distantrestrain
	use params
      use lengths
      use distres
      use RCN
      implicit integer(i-z)
!      parameter(ndim=1500)
!      parameter(nvec=416)
!      COMMON/RCN/Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
      COMMON/RCN1/er1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
!      common/lengths/Lch,Lch1,Lch2
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
!      common/distres/er4,es3c
      common/resnumber/Ncom,Ndis,accur
      common/cutoff/cut1a,cut2a,cut3a,cut4a,cut1b,cut2b,cut3b,cut4b

      dimension r1(10000),r2(10000),dis(10000),deviation(10000)
      
      read(17,*)Ndis
      if(Ndis.gt.10000)then
         write(*,*)'warning: too many short range distance restraints!'
         Ndis=10000
      endif
      do i=1,Ndis
         read(17,*)r1(i),r2(i),nothing,dis(i),deviation(i)
         dis(i)=dis(i)/0.87
         deviation(i)=deviation(i)/0.87
         if(deviation(i).lt.0.5)deviation(i)=0.5
      enddo

      do i=1,Lch
         Mdis(i)=0              !number of prediction for 'i'
         do j=1,Ndis
            if(r1(j).eq.i)then
               Mdis(i)=Mdis(i)+1
               kdis(i,Mdis(i))=r2(j) !r2(j) with 'i'
               dist(i,Mdis(i))=dis(j) !predicted distance for i<->r2(j)
               dev(i,Mdis(i))=deviation(j)
            endif
            if(r2(j).eq.i)then
               Mdis(i)=Mdis(i)+1
               kdis(i,Mdis(i))=r1(j)
               dist(i,Mdis(i))=dis(j) !predicted distance for i<->r2(j)
               dev(i,Mdis(i))=deviation(j)
            endif
         enddo
      enddo
c      dilim=1.5*Ndis          !background number for derviation
      dilim=0.5*Ndis          !background number for derviation

!$acc update device(Mdis,kdis,dist,dev)

c      write(*,*)
c      write(*,*)'----------- distant map ---------------'
c      nnc=0
c      do i=1,Lch
c         nnc=nnc+Mdis(i)
c         write(*,12)i,Mdis(i),(Kdis(i,j),dist(i,j)*0.87,j=1,Mdis(i))
c 12      format(i4,':',i3,20(i4,'-'f5.2))
c      enddo
c      write(*,*)'Number of distmap=',nnc,' Lch=',Lch

c^^^^^^^^^^^^^^^^^^ read contact restrains finished ^^^^^^^^^^^^^^^^^
      return
      end
      subroutine read_longdistantrestrain
	use params
      use lengths
      use longdist
      implicit integer(i-z)
!!      parameter(ndim=1500)
!!      parameter(nvec=416)
!      common/lengths/Lch,Lch1,Lch2
!      COMMON/longdist/MdisL(ndim),kdisL(ndim,500),distL(ndim,500)
c     maximum N=4*L/10, if L=100, maximum N=400
      dimension r1(400000),r2(400000),dis(400000)

      read(22,*)NdisL
      if(NdisL.gt.400000)NdisL=400000
      do i=1,NdisL
         read(22,*)r1(i),r2(i),dis(i)
         dis(i)=dis(i)/0.87
      enddo
      
      do i=1,Lch
         MdisL(i)=0             !number of prediction for 'i'
         do j=1,NdisL
            if(r1(j).eq.i)then
               MdisL(i)=MdisL(i)+1
               kdisL(i,MdisL(i))=r2(j) !r2(j) with 'i'
               distL(i,MdisL(i))=dis(j) !predicted distance for i<->r2(j)
               if(MdisL(i).ge.500)goto 101
            endif
            if(r2(j).eq.i)then
               MdisL(i)=MdisL(i)+1
               kdisL(i,MdisL(i))=r1(j)
               distL(i,MdisL(i))=dis(j) !predicted distance for i<->r2(j)
               if(MdisL(i).ge.500)goto 101
            endif
         enddo
 101     continue
         if(MdisL(i).gt.500)then !size>5000/4
            write(*,*)i,MdisL(i),'N_res>500 for one residue, exit'
            stop
         endif
      enddo
!$acc update device(MdisL,kdisL,distL)

c      write(*,*)
c      write(*,*)'----------- distant map ---------------'
c      nnc=0
c      do i=1,Lch
c        nnc=nnc+MdisL(i)
c        write(*,12)i,MdisL(i),(KdisL(i,j),distL(i,j)*0.87,j=1,MdisL(i))
c 12     format(i4,':',i3,20(i4,'-'f5.2))
c      enddo
c      write(*,*)'Number of distmap=',nnc,' Lch=',Lch
c      stop
      
c^^^^^^^^^^^^^^^^^^ read contact restrains finished ^^^^^^^^^^^^^^^^^
      return
      end
      subroutine read_exp
	use params
      use lengths
      use expose
      implicit integer(i-z)
!      parameter(ndim=1500)
!      common/lengths/Lch,Lch1,Lch2
!      common/expose/mp(20,ndim),area(ndim)
      common/nana1/nana

      na=nana
      read(27,*)
      do i=1,Lch
         read(27,*)itmp,(mp(j,i),j=1,na)
         area(i)=0
         do j=1,na
            if(mp(j,i).eq.0)mp(j,i)=-1 !-1, bury; 1, expose
            area(i)=area(i)+mp(j,i)
         enddo
c         write(*,*)i,mp(1,i),mp(2,i),mp(12,i)
      enddo
      close(27)
c      stop
!$acc update device(area)

c^^^^^^^^^^^^^^^^^^ read contact restrains finished ^^^^^^^^^^^^^^^^^
      return
      end
      subroutine read_CAcontact
	use params
      use lengths
      use CAcontact
      implicit integer(i-z)
!      parameter(ndim=1500)
!      common/lengths/Lch,Lch1,Lch2
      common/zscore/izscore
!      common/CAcontact/McomCA(ndim),KcomCA(ndim,100),aweigCA(ndim,ndim)
      common/CAcontact1/dist_CA_cut
      DIMENSION r1(8000),r2(8000)

      dist_CA_cut_tmp=6.5
      if(izscore.eq.1)then
         dist_CA_cut=(dist_CA_cut_tmp/0.87)**2
      elseif(izscore.eq.2)then
         dist_CA_cut=(dist_CA_cut_tmp/0.87)**2
      else
         dist_CA_cut=(dist_CA_cut_tmp/0.87)**2
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
      read(18,*)ntmp
      i_c=0
      do i=1,ntmp
         read(18,*)i1,i2,conf
c         conf=float(i3)/float(i4)
         if(conf.gt.1)conf=1
         if(conf.gt.cut_min)then
            i_c=i_c+1
            if(i_c.gt.8000)then
               write(*,*)'Too many COMBCA restraints!!!!!'
               cut_min=cut_min*1.1
               rewind(18)
               goto 11
            endif
            r1(i_c)=i1
            r2(i_c)=i2
            if(conf.gt.cut0)then
               aweigCA(i1,i2)=1+abs(conf-cut0)*4
            else
               aweigCA(i1,i2)=1-abs(conf-cut0)*2
            endif
            aweigCA(i2,i1)=aweigCA(i1,i2)
c            write(*,*)i1,i2,conf,aweigCA(i1,i2)
         endif
      enddo
      NcomCA=i_c

ccc   map r1,2(i) into McomCA(i),KcomCA(i,McomCA(i))------------>
      do i=1,Lch
         McomCA(i)=0              !number of contacts with 'i'
         do j=1,NcomCA
            if(r1(j).eq.i)then
               McomCA(i)=McomCA(i)+1
               KcomCA(i,McomCA(i))=r2(j)
            endif
            if(r2(j).eq.i)then
               McomCA(i)=McomCA(i)+1
               KcomCA(i,McomCA(i))=r1(j)
            endif
         enddo
      enddo
!$acc update device(McomCA,KcomCA,aweigCA)

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
      subroutine reset_temperature
	use params
      use lengths
      use RES
      use distres
      implicit integer(i-z)
!      parameter(ndim=1500)      !number of residues
!      parameter(nrep=100)       !number of replicas
!      parameter(nvec=416)       !number of vectors
!      common/lengths/Lch,Lch1,Lch2
      common/resnumber/Ncom,Ndis,accur
      common/commonuse2/atemp1,atemp2,N_rep,phot
      COMMON/RCN1/er1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
!      common/distres/er4,es3c
!      COMMON/RES/ER3,er5,er6,er7,Mcom(ndim),Kcom(ndim,100)
      real r_dis,r_con,r_dev,T10,T20,T1a,T2a

      if(er1+er3+er4.lt.0.1)return !without restrains
      
      a_rest=float(Ncom)/(1.3*Lch)
      a_rest=sqrt(a_rest)
      if(a_rest.lt.0.875)a_rest=0.875 !because of 80/70 -> rest/without_rest
      atemp2=atemp2*a_rest
      if(atemp2.gt.115) atemp2=115.0

      return
      end
      subroutine set_temperature
	use params
      implicit integer(i-z)
!      parameter(nrep=100)
      common/commonuse1/random,ncycle
      common/commonuse2/atemp1,atemp2,N_rep,phot
      common/sw1/aT_rep(nrep),E_rep(nrep)
      common/initialinput/switch,k_cycle,k_phot,N_ann
      common/thrtem/aT_ann(100)
      common/aTs/aTs1,aTs2,aTs_rep(nrep),aTs
      common/aTTs/aTTs1,aTTs2,aTTs_rep(nrep),aTTs

***********for normal run ***********************************
      do i=1,N_REP
         aT_rep(i)=atemp1*(atemp2/atemp1)**(float(i-1)/(N_rep-1))
         aTs_rep(i)=aTs1*(aTs2/aTs1)**(float(i-1)/(N_rep-1))
         aTTs_rep(i)=aTTs1*(aTTs2/aTTs1)**(float(i-1)/(N_rep-1))
c     write(*,*)i,aT_rep(i)
      enddo

c      stop
c^^^^^^^^^^^^^^^^^ set aT_rep(i) finished ^^^^^^^^^^^^^^^
      return
      end
      subroutine set_EHB
	use params
      use lengths
      use ENERGY
      use seqe
      implicit integer(i-z)
!      parameter(ndim=1500)
!      parameter(nvec=416)
!      COMMON/ENERGY/EH5,ES3,ES3a,ES3b,EH1,EH1a,EH4,EHBIJ(ndim,ndim)
!      common/lengths/Lch,Lch1,Lch2
!      common/seqe/seq(ndim),sec(ndim)

c
c     EHBIJ - set-up secondary structure dependent
c     strength of the hyrogen bond network - stronger for helices
c     and beta-sheets
c

      do i=1,Lch
         is=sec(i)
         do j=1,Lch
            js=sec(j)
            EHBIJ(i,j)=1
            if(iabs(i-j).eq.3.and.is.eq.2.and.js.eq.2)then
               EHBIJ(i,j)=EHBIJ(i,j)+0.5 !helix, enhanced
            endif
            if(is.eq.4.or.js.eq.4) then
               if(is*js.ne.8.and.iabs(i-j).gt.4)then
                  EHBIJ(i,j)=EHBIJ(i,j)+0.5 !beta-beta, enhanced
               endif
            endif
         enddo
      enddo

c^^^^^^^^^^^ set H_bond finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
      subroutine prepare_vectors
	use params
      implicit integer(i-z)
!      parameter(nvec=416)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      dimension n(100)
      common/lattice/m_latt,latt1,latt2

      do i=1,100
         n(i)=0
      enddo

      nwmax=0
      aaa=0
      nn=5
      do x=-nn,nn
         do y=-nn,nn
            do z=-nn,nn
               vector(x,y,z)=0
               r=x*x+y*y+z*z
               if(r.ge.latt1.and.r.le.latt2) then
                  nwmax=nwmax+1
                  vx(nwmax)=x
                  vy(nwmax)=y
                  vz(nwmax)=z
                  vector(x,y,z)=nwmax
c                  write(*,*)nwmax,vx(nwmax),vy(nwmax),vz(nwmax),r,
c     $                 sqrt(float(r)*0.87*0.87)
                  n(r)=n(r)+1
                  aaa=aaa+sqrt(float(r)*0.87*0.87)
               endif
            enddo
         enddo
      enddo
      aaa=aaa/float(nwmax)
      write(20,*)'n1_all=',nwmax,'  <vr>=',aaa

c      do i=1,nwmax
c         write(*,*)i,vx(i),vy(i),vz(i)
c      enddo
c      stop

c      do i=10,30
c         write(*,*)i,n(i),sqrt(float(i)*0.87*0.87)
c      enddo
c      stop
      
c     i=1,5
c           10          24   2.751182    
c           11          24   2.885463    
c           12           8   3.013768    +
c           13          24   3.136830    +
c           14          48   3.255242    x
c           15           0   3.369496    
c           16           6   3.480000    x
c           17          48   3.587102    x
c           18          36   3.691097    x
c           19          24   3.792242    x
c           20          24   3.890758    x
c           21          48   3.986841    x
c           22          24   4.080662    x
c           23           0   4.172373    
c           24          24   4.262112    x
c           25          30   4.350000    x
c           26          72   4.436147    +
c           27          32   4.520653    
c           28           0   4.603607    
c           29          72   4.685093    
c           30          48   4.765186    
c nwmax=         616  <vr>=   3.982909    
c nwmax=         312  <vr>=   3.809868    

c      stop
      return
      end
      subroutine prepare_neighbors
	use params
      use lengths
      implicit integer(i-z)
!                parameter(nvec=416)
      common/logica/goodc
      logical goodc(nvec,nvec)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/three/angle(nvec,nvec)
!      common/lengths/Lch,Lch1,Lch2

      common/shift/u21(nvec,nvec),m12(nvec),u1(nvec,100),u2(nvec,100)
      common/nwmax1/nvecr

      max_m12=0
      do i=1,nvecr
         m12(i)=0
      enddo

      mmm=0
      nnn=0
      kkk=0
      do i=1,nvecr
      do j=1,nvecr
         u21(i,j)=0
         a2=vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)
         b2=vx(j)*vx(j)+vy(j)*vy(j)+vz(j)*vz(j)
         c2=(vx(i)+vx(j))**2+(vy(i)+vy(j))**2+(vz(i)+vz(j))**2
         cosangle=(a2+b2-c2)/(2*sqrt(a2*b2))
         angle(i,j)=acos(cosangle)*180/3.1415926
c     in database, angle is in [65,165];
         if(angle(i,j).gt.65.and.angle(i,j).lt.165)then
            goodc(i,j)=.true.
            mmm=mmm+1
            ijx=vx(i)+vx(j)
            ijy=vy(i)+vy(j)
            ijz=vz(i)+vz(j)
            do k=1,nvecr
               if(vx(k).eq.ijx.and.vy(k).eq.ijy.and.vz(k).eq.ijz)then
                  kkk=kkk+1
                  u21(i,j)=k    !vi+vj=vk
                  m12(k)=m12(k)+1
                  u1(k,m12(k))=i
                  u2(k,m12(k))=j
                  if(max_m12.lt.m12(k))max_m12=m12(k)
                  goto 10
               endif
            enddo
 10         continue
         else
            goodc(i,j)=.false.
         endif
         nnn=nnn+1
c         write(*,*)i,j,mmm,nnn,angle(i,j),goodc(i,j)
      enddo
      enddo

      n=0
      do i=1,nvecr
         r=vx(i)**2+vy(i)**2+vz(i)**2
         if(m12(i).gt.0)n=n+1
c     if(r.gt.17)write(*,*)i,r,m12(i)
      enddo
c      write(*,*)'n2_all=',nnn
c      write(*,*)'n2good=',mmm,'  n21=',kkk
c      write(*,*)'n1_all=',nvecr,'  n12=',n

c      stop
      return
      end
      subroutine prepare_beta
	use params
      use hb
      use bisec
      use sg
      implicit integer(i-z)
!                parameter(nvec=416)
      common/logica/goodc
      logical goodc(nvec,nvec)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
!      common/bisec/cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
!      common/hb/hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
      common/three/angle(nvec,nvec)
!      common/sg/gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
!      common/cb/hx(nvec,nvec,0:19),hy(nvec,nvec,0:19),hz(nvec,nvec,0:19)
      common/beta1/axalf(0:19),ayalf(0:19),azalf(0:19)
      common/beta2/axbet(0:19),aybet(0:19),azbet(0:19)
      common/beta3/bxalf(0:19),byalf(0:19),bzalf(0:19)
      common/beta4/bxbet(0:19),bybet(0:19),bzbet(0:19)
      common/nwmax1/nvecr


      do k=0,19
         read(4,*)
         read(4,*)axalf(k),ayalf(k),azalf(k),axbet(k),aybet(k),azbet(k)
      enddo
      do k=0,19
         read(4,*)
         read(4,*)bxalf(k),byalf(k),bzalf(k),bxbet(k),bybet(k),bzbet(k)
      enddo
      CLOSE(4)                  !sidecent.comm

ccccccccccc define hydrogen-bond, C_beta, C_group for good (i,j)----->
      esp=0.000001
      do 101 i=1,nvecr
      do 102 j=1,nvecr
         avi=sqrt(float(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)))
         avxi=vx(i)/avi
         avyi=vy(i)/avi
         avzi=vz(i)/avi
         avj=sqrt(float(vx(j)*vx(j)+vy(j)*vy(j)+vz(j)*vz(j)))
         avxj=vx(j)/avj
         avyj=vy(j)/avj
         avzj=vz(j)/avj
         
ccc   if vi and vj is parallel, a is ok but b=c=0 ------->
         if(abs(avxi-avxj).lt.esp)then
            if(abs(avyi-avyj).lt.esp)then
               if(abs(avzi-avzj).lt.esp)then
                  ax=avxi+avxj
                  ay=avyi+avyj
                  az=avzi+avzj
                  aaa=sqrt(ax*ax+ay*ay+az*az)
                  ax=ax/aaa     !(vi+vj)/|vi+vj|
                  ay=ay/aaa     !(vi+vj)/|vi+vj|
                  az=az/aaa     !(vi+vj)/|vi+vj|
                  
c     calculate b=a(x)u, u=(1,1,1):
                  bx=ay*1-az*1
                  by=az*1-ax*1
                  bz=ax*1-ay*1
                  bbb=sqrt(bx*bx+by*by+bz*bz)
                  bx=bx/bbb     ! a(x)1/|a(x)1|
                  by=by/bbb     ! a(x)1/|a(x)1|
                  bz=bz/bbb     ! a(x)1/|a(x)1|
                  
c     calculate c=a(x)b:
                  cx=ay*bz-az*by
                  cy=az*bx-ax*bz
                  cz=ax*by-ay*bx
                  ccc=sqrt(cx*cx+cy*cy+cz*cz)
                  cx=cx/ccc     ! a(x)b/|a(x)b|
                  cy=cy/ccc     ! a(x)b/|a(x)b|
                  cz=cz/ccc     ! a(x)b/|a(x)b|
                  cax(i,j)=cx
                  cay(i,j)=cy
                  caz(i,j)=cz
                  goto 29
               endif
            endif
         endif
ccc   check if vi and vj is anti-parallel, c is ok, b=a=0 ------->
         if(abs(avxi+avxj).lt.esp)then
            if(abs(avyi+avyj).lt.esp)then
               if(abs(avzi+avzj).lt.esp)then
                  cx=avxi-avxj
                  cy=avyi-avyj
                  cz=avzi-avzj
                  ccc=sqrt(cx*cx+cy*cy+cz*cz)
                  cx=cx/ccc     !(vi-vj)/|vi-vj|
                  cy=cy/ccc     !(vi-vj)/|vi-vj|
                  cz=cz/ccc     !(vi-vj)/|vi-vj|
                  cax(i,j)=cx   !(vi-vj)/|vi-vj|
                  cay(i,j)=cy   !(vi-vj)/|vi-vj|
                  caz(i,j)=cz   !(vi-vj)/|vi-vj|

c     calculate a=c(x)u, u=(1,1,1):
                  ax=cy*1-cz*1
                  ay=cz*1-cx*1
                  az=cx*1-cy*1
                  aaa=sqrt(ax*ax+ay*ay+az*az)
                  ax=ax/aaa     ! c(x)1/|c(x)1|
                  ay=ay/aaa     ! c(x)1/|c(x)1|
                  az=az/aaa     ! c(x)1/|c(x)1|
                  
c     calculate b=c(x)a:
                  bx=cy*az-cz*ay
                  by=cz*ax-cx*az
                  bz=cx*ay-cy*ax
                  bbb=sqrt(bx*bx+by*by+bz*bz)
                  bx=bx/bbb     ! c(x)a/|c(x)a|
                  by=by/bbb     ! c(x)a/|c(x)a|
                  bz=bz/bbb     ! c(x)a/|c(x)a|
                  goto 29
               endif
            endif
         endif
         
         ax=avxi+avxj
         ay=avyi+avyj
         az=avzi+avzj
         aaa=sqrt(ax*ax+ay*ay+az*az)
         ax=ax/aaa              !(vi+vj)/|vi+vj|
         ay=ay/aaa              !(vi+vj)/|vi+vj|
         az=az/aaa              !(vi+vj)/|vi+vj|
         
         bx=avyi*avzj-avzi*avyj
         by=avzi*avxj-avxi*avzj
         bz=avxi*avyj-avyi*avxj
         bbb=sqrt(bx*bx+by*by+bz*bz)
         bx=bx/bbb             ! vi(x)vj/|vi(x)vj|
         by=by/bbb             ! vi(x)vj/|vi(x)vj|
         bz=bz/bbb             ! vi(x)vj/|vi(x)vj|

         cx=avxi-avxj
         cy=avyi-avyj
         cz=avzi-avzj
         ccc=sqrt(cx*cx+cy*cy+cz*cz)
         cx=cx/ccc              !(vi-vj)/|vi-vj|
         cy=cy/ccc              !(vi-vj)/|vi-vj|
         cz=cz/ccc              !(vi-vj)/|vi-vj|
         cax(i,j)=cx            !(vi-vj)/|vi-vj|
         cay(i,j)=cy            !(vi-vj)/|vi-vj|
         caz(i,j)=cz            !(vi-vj)/|vi-vj|
         
 29      continue

         goto 39
c     check if aaa/bbb/ccc=0 ------->
         if(aaa.lt.esp.or.bbb.lt.esp.or.ccc.lt.esp)then
            write(*,19)i,j,avxi,avyi,avzi,aaa,bbb,ccc
         endif
c     check if a.b=b.c=a.c=0 ------->
         ab=ax*bx+ay*by+az*bz
         bc=bx*cx+by*cy+bz*cz
         ac=ax*cx+ay*cy+az*cz
         if(ab.gt.esp.or.ab.gt.esp.or.ac.gt.esp)then
            write(*,19)i,j,avxi,avyi,avzi,ab,bc,ac
         endif
c         write(*,19)i,j,avxi,avyi,avzi,ab,bc,ac
 19      format(2i5,10f8.3)
 39       continue

c     H-bond (unit vector):
         hbx(i,j)=bx
         hby(i,j)=by
         hbz(i,j)=bz

c     side-chain coordinate from C_a to side-chain ---------------->
         do k=0,19
            if(angle(i,j).lt.105) then ! alpha-helix or turn like
               gx(i,j,k)=(axalf(k)*ax+ayalf(k)*bx+azalf(k)*cx)/0.87
               gy(i,j,k)=(axalf(k)*ay+ayalf(k)*by+azalf(k)*cy)/0.87
               gz(i,j,k)=(axalf(k)*az+ayalf(k)*bz+azalf(k)*cz)/0.87

               hx(i,j,k)=(bxalf(k)*ax+byalf(k)*bx+bzalf(k)*cx)/0.87
               hy(i,j,k)=(bxalf(k)*ay+byalf(k)*by+bzalf(k)*cy)/0.87
               hz(i,j,k)=(bxalf(k)*az+byalf(k)*bz+bzalf(k)*cz)/0.87
            else                ! beta-sheet
               gx(i,j,k)=(axbet(k)*ax+aybet(k)*bx+azbet(k)*cx)/0.87
               gy(i,j,k)=(axbet(k)*ay+aybet(k)*by+azbet(k)*cy)/0.87
               gz(i,j,k)=(axbet(k)*az+aybet(k)*bz+azbet(k)*cz)/0.87

               hx(i,j,k)=(bxbet(k)*ax+bybet(k)*bx+bzbet(k)*cx)/0.87
               hy(i,j,k)=(bxbet(k)*ay+bybet(k)*by+bzbet(k)*cy)/0.87
               hz(i,j,k)=(bxbet(k)*az+bybet(k)*bz+bzbet(k)*cz)/0.87
            endif
         enddo
 102  continue
 101  continue

!$acc update device(hbx,hby,hbz,gx,gy,gz)

c      stop
      return
      end
      subroutine prepare_frg
	use params
      use lengths
      use fr
      use seqe
      implicit integer(i-z)
!      parameter(ndim=1500)
!                parameter(nvec=416)
!      common/fr/frga(ndim),frgb(ndim)
!      common/seqe/seq(ndim),sec(ndim)
!      common/lengths/Lch,Lch1,Lch2

      do i=1,Lch
         frga(i)=0.0
         frgb(i)=0.0
      enddo

      do i=1,Lch-7
         q=0
         do j=i,i+7
            if(sec(j).eq.2) q=q+2 !helix structure.
         enddo
         if(q.eq.16)then        !8 continue alpha-residues
            frga(i)=10.5/0.87   !distance for 7 alpha-bonds
         endif
      enddo

      do i=1,Lch-6
         q=0
         do j=i+1,i+5
            if(sec(j).eq.4) q=q+4 !beta structure
         enddo
         if(q.eq.20)then        !5 continue beta-residues
            if(sec(i).ne.2.and.sec(i+6).ne.2)then
               frgb(i)=19.1/0.87 !distance for 6 beta-bonds
            endif
         endif
      enddo

c      do i=1,Lch
c         write(*,*)i,seq(i),sec(i),frga(i),frgb(i)
c      enddo
!$acc update device(frga,frgb)
      return
      end
      subroutine count_restrains
	use params
      use lengths
      use RES
      use concutt
      use sg
      use seqe
      use longdist
      use RCN
      use CAcontact
      implicit integer (i-z)
!!      parameter(nrep=100)
!!      parameter(ndim=1500)	!maximum length of chain-length
!!      parameter(nvec=416)

!      common/lengths/Lch,Lch1,Lch2
      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
!      common/sg/gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
!      common/seqe/seq(ndim),sec(ndim)
!      common/concutt/concut(0:19,0:19),concut2(0:19,0:19),concut_sc
!      COMMON/RCN/Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
!      COMMON/RES/ER3,er5,er6,er7,Mcom(ndim),Kcom(ndim,100)
!      common/CAcontact/McomCA(ndim),KcomCA(ndim,100),aweigCA(ndim,ndim)
      common/CAcontact1/dist_CA_cut
!      COMMON/longdist/MdisL(ndim),kdisL(ndim,500),distL(ndim,500)
!      common/CA8/McomCA8(ndim),KcomCA8(ndim,100),aweigCA8(ndim,ndim)
      common/CAcontact2/dist_CA_cut8
 
      dimension ax(ndim),ay(ndim),az(ndim) !SG
      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim) !CA
      
      common/countres/t_comb,t_dist,t_combCA,s_comb,s_dist,s_combCA
      common/countres1/t_combCA8,t_distL,s_combCA8,s_distL,N_resc
      
      do i=1,3
         N_resc=N_resc+1
         do j=1,Lch
            x(j)=xrep(j,i)
            y(j)=yrep(j,i)
            z(j)=zrep(j,i)
         enddo
         do j=1,Lch1
            wx=x(j+1)-x(j)
            wy=y(j+1)-y(j)
            wz=z(j+1)-z(j)
            ica(j)=vector(wx,wy,wz) !identify order number of each bond-vector
         enddo
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         do j=1,Lch
            ax(j)=x(j)+gx(ica(j-1),ica(j),seq(j))
            ay(j)=y(j)+gy(ica(j-1),ica(j),seq(j))
            az(j)=z(j)+gz(ica(j-1),ica(j),seq(j))
         enddo

         do j=1,Lch
c     comb.dat--------->
            do k=1,Mcom(j)
               m=Kcom(j,k)      !k'th contact with j
               if(m.gt.j)then
                  t_comb=t_comb+1
                  dis=(ax(j)-ax(m))**2+(ay(j)-ay(m))**2
     $                 +(az(j)-az(m))**2
                  if(dis.le.concut2(seq(j),seq(m)))then
                     s_comb=s_comb+1
                  endif
               endif
            enddo
c     dist.dat--------->
            do k=1,Mdis(j)
               m=kdis(j,k)      !j-m restraints
               if(m.gt.j) then  !to avoid repeat
                  t_dist=t_dist+1
                  dis2=(x(j)-x(m))**2+(y(j)-y(m))**2+(z(j)-z(m))**2
                  dis=sqrt(dis2)
                  err=abs(dis-dist(j,k)) !dist: predicted dis
                  if(err.lt.dev(j,k)) then !dev: deviation for arca
                     s_dist=s_dist+1
                  endif
               endif
            enddo
c     distL.dat--------->
            do k=1,MdisL(j)
               m=kdisL(j,k)     !j-m restraints
               if(m.gt.j) then  !to avoid repeat
                  t_distL=t_distL+1
                  dis2=(x(j)-x(m))**2+(y(j)-y(m))**2+(z(j)-z(m))**2
                  dis=sqrt(dis2)
                  err=abs(dis-distL(j,k)) !dist: predicted dis
                  if(err.lt.1.5) then !dev: deviation for arca
                     s_distL=s_distL+1
                  endif
               endif
            enddo
c     combCA.dat--------->
            do k=1,McomCA(j)
               m=KcomCA(j,k)      !k'th contact with j
               if(m.gt.j)then
                  t_combCA=t_combCA+1
                  dis=(x(j)-x(m))**2+(y(j)-y(m))**2+(z(j)-z(m))**2
                  if(dis.le.dist_CA_cut)then
                     s_combCA=s_combCA+1
                  endif
               endif
            enddo
c     combCA8.dat--------->
            do k=1,McomCA8(j)
               m=KcomCA8(j,k)   !k'th contact with j
               if(m.gt.j)then
                  t_combCA8=t_combCA8+1
                  dis=(x(j)-x(m))**2+(y(j)-y(m))**2+(z(j)-z(m))**2
                  if(dis.le.dist_CA_cut8)then
                     s_combCA8=s_combCA8+1
                  endif
               endif
            enddo
ccc   
         enddo
      enddo

c^^^^^^^^^^^^^^^^^ restraints count finished ^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
      subroutine prepare_move2
	use params
      implicit integer(i-z)
!      parameter(nvec=416)
      common/logica/goodc
      logical goodc(nvec,nvec)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)
      common/w2a/w21(-10:10,-10:10,-10:10,80)
      common/w2b/w22(-10:10,-10:10,-10:10,80)
      common/w2c/Nw(-10:10,-10:10,-10:10)
      common/nwmax1/nvecr
      common/lattice/m_latt,latt1,latt2

      do i=-10,10
         do j=-10,10
            do k=-10,10
               Nw(i,j,k)=0
            enddo
         enddo
      enddo

      max=0
      maxa=0
      nnn=0
      mmm=0

      do 101 i=1,nvecr
         do 102 j=1,nvecr
            if(goodc(i,j))then
               ijx=vx(i)+vx(j)  !vx in (-5,5)
               ijy=vy(i)+vy(j)
               ijz=vz(i)+vz(j)
***   based on ijx:
               Nw(ijx,ijy,ijz)=Nw(ijx,ijy,ijz)+1 !based on r
               w21(ijx,ijy,ijz,Nw(ijx,ijy,ijz))=i
               w22(ijx,ijy,ijz,Nw(ijx,ijy,ijz))=j
               if(maxa.lt.Nw(ijx,ijy,ijz))maxa=Nw(ijx,ijy,ijz)
***   based on i,j:
               Np2(i,j)=0
               do 103 ii=1,nvecr
                  jx=ijx-vx(ii)
                  jy=ijy-vy(ii)
                  jz=ijz-vz(ii)
                  jr=jx*jx+jy*jy+jz*jz
                  if(jr.ge.latt1.and.jr.le.latt2)then
                     jj=vector(jx,jy,jz)
                     if(goodc(ii,jj))then
                        Np2(i,j)=Np2(i,j)+1 !based on i,j
                        v21(i,j,Np2(i,j))=ii
                        v22(i,j,Np2(i,j))=jj
                        mmm=mmm+1 !number of total memory occupied by vx21.
                     endif
                  endif
 103           continue
               if(max.lt.Np2(i,j))max=Np2(i,j)
c               write(*,*)i,j,Np2(i,j),max,maxa
               nnn=nnn+1        !number of possible pairs
            endif
 102     continue
 101  continue

ccc among all 312*312=97344 pairs, 67272 pairs are legal (~70%).
ccc all Np2(i,j)>=2, i.e. there are at least one other path for any pair.
ccc <Np2(i,j)>=27.
      
      write(20,*)'maximum of Np2(i,j)=',max,maxa
      write(20,*)'number of possible pair=',nnn
      write(20,*)'sum of Np2(i,j), total memory=',mmm
cccc  the following is the cases without goodc on prefabracated pairs:
c              N_v  N_pare N_pp   mmm
c     [14,25], 312, 97344, 67272, 1830000 (25% memory, within memory)

cccc the following is cases without goodc limitation:
c     [14,25], 306, 93636, 3872214 (90% memory, beyond memory)
c     [14,24], 282, 79524, 2923758 (80% memory, beyond memory)
c      stop

      return
      end
      subroutine prepare_move3
	use params
      implicit integer(i-z)
!      parameter(nvec=416)
      common/logica/goodc
      logical goodc(nvec,nvec)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/move31/v31(-15:15,-15:15,-15:15,6)
      common/move32/v32(-15:15,-15:15,-15:15,6)
      common/move33/v33(-15:15,-15:15,-15:15,6)
      common/move34/Np3(-15:15,-15:15,-15:15)
      common/nwmax1/nvecr

c      dimension v31(-15:15,-15:15,-15:15,100)
c      dimension v32(-15:15,-15:15,-15:15,100)
c      dimension v33(-15:15,-15:15,-15:15,100)
c      dimension Np3(-15:15,-15:15,-15:15)

      do i=-15,15
         do j=-15,15
            do k=-15,15
               Np3(i,j,k)=0
            enddo
         enddo
      enddo


      num3=0
      max=0
      do 101 i=1,nvecr
         do 102 j=1,nvecr
            if(goodc(i,j))then
               do 103 k=1,nvecr
                  if(goodc(j,k))then
                     num3=num3+1
                     rx=vx(i)+vx(j)+vx(k)
                     ry=vy(i)+vy(j)+vy(k)
                     rz=vz(i)+vz(j)+vz(k)
                     Np3(rx,ry,rz)=Np3(rx,ry,rz)+1
                     v31(rx,ry,rz,Np3(rx,ry,rz))=i
                     v32(rx,ry,rz,Np3(rx,ry,rz))=j
                     v33(rx,ry,rz,Np3(rx,ry,rz))=k
                     if(max.le.Np3(rx,ry,rz)) max=Np3(rx,ry,rz)
                     write(*,*)i,j,k,max
                  endif
 103           continue
            endif
 102     continue
 101  continue

      write(*,*)'number of move3=',num3,max

      stop
      return
      end
      subroutine get_acorder
	use params
      use lengths
      use seqe
      use order
      implicit integer(i-z)
!      parameter(ndim=1500)
!      parameter(nvec=416)
!      common/seqe/seq(ndim),sec(ndim)
!      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
!      common/lengths/Lch,Lch1,Lch2
      
c     data struct /' coil','helix',' turn',' beta','    ?'/

************* number of predicted structures *************
      n_H=0                     !number of Helix
      n_E=0                     !number of Extension
      do i=1,Lch
         if(sec(i).eq.2) n_H=n_H+1
         if(sec(i).eq.4) n_E=n_E+1
      enddo

      if(n_H+n_E.lt.2)then      !use alpha=H/E
         if(Lch.lt.50)then
            alph=0.252
         else if(Lch.lt.70)then
            alph=0.230
         else if(Lch.lt.90)then
            alph=0.212
         else if(Lch.lt.110)then
            alph=0.198
         else if(Lch.lt.130)then
            alph=0.168
         else if(Lch.lt.150)then
            alph=0.176
         else if(Lch.lt.170)then
            alph=0.183
         else if(Lch.lt.190)then
            alph=0.160
         else
            alph=0.134
         endif
      else                      !use alpha=aH+bE
         a1=float(n_H)/float(n_H+n_E)
         a2=float(n_E)/float(n_H+n_E)
         if(Lch.lt.50)then
            alph=0.116*a1+0.324*a2
         else if(Lch.lt.70)then
            alph=0.119*a1+0.357*a2
         else if(Lch.lt.90)then
            alph=0.115*a1+0.280*a2
         else if(Lch.lt.110)then
            alph=0.105*a1+0.259*a2
         else if(Lch.lt.130)then
            alph=0.132*a1+0.269*a2
         else if(Lch.lt.150)then
            alph=0.105*a1+0.272*a2
         else if(Lch.lt.170)then
            alph=0.114*a1+0.186*a2
         else if(Lch.lt.190)then
            alph=0.116*a1+0.197*a2
         else
            alph=0.107*a1+0.184*a2
         endif
      endif

      acorder=alph*Lch
c^^^^^^^^^^^^^^^^^^ contact order done ^^^^^^^^^^^^^^^^^^
      return
      end
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
         if(n_walk.gt.10000)then
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
