
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
          integer,  parameter :: ndim = 1500
          integer,  parameter ::  nrep = 100
          integer,  parameter :: nvec = 416
      end module params


      module trackn
          use params
            integer n_tem(100)
! !$OMP declare target(n_tem)
!$OMP declare target(n_tem)
      end module trackn
 

        
      module chainm
         use params
          real :: mv(ndim)
! !$OMP declare target(mv)
!$OMP declare target(mv)
! !$OMP end declare target  
      end module chainm
             
      module eigen
         real :: AA(3,3), EE(3), HH(3,3)
! !$OMP declare target(AA,EE,HH)
!$OMP declare target(AA,EE,HH)         
      end module


      module ehbc
          use params
          real :: envir(0:15,0:15,0:15,0:19,4),en1
!$OMP declare target(envir,en1)   
      end module ehbc

      module backup1
         use params
         integer nopp(ndim),nomm(ndim),noaa(ndim),nhbnn(ndim)
!$OMP declare target(nopp,nomm,noaa,nhbnn)
      end module backup1


      module shape 
         use params
         real :: amx,amy,amz,afs(ndim),afsn(ndim)
!$OMP declare target(amx,amy,amz,afs,afsn)
      end module shape

      module stick
	  use params
          integer nstick,nrmsd
          integer iq(ndim,nrep)
          real :: ax00(ndim,nrep),ay00(ndim,nrep),az00(ndim,nrep)
          real :: bx00(ndim),ermsd,by00(ndim),bz00(ndim),armsd_min0,astick
          integer itemp0,icycle,icycle0
!$OMP declare target (ax00,iq,ay00,az00,bx00,by00,bz00)
      end module

      module short
         use params 
          integer  IBIN(-300:300)
          real :: asr(ndim,-12:12),csr(ndim,2)
          real :: codevsum, didevsum 
!$OMP declare target(IBIN,asr,csr)
      end module short 

      module ENERGY
         use params
         real :: EH5,ES3,ES3a,ES3b,EH1,EH1a,EH4,EHBIJ(ndim,ndim)
!$OMP declare target(EH5,ES3,ES3a,ES3b,EH1,EH1a,EH4,EHBIJ)
      end module ENERGY
      
      module svm1
         use params
            integer Mcon(15,ndim),Kcon(15,ndim,100)
            real :: awei(15,ndim,ndim) 
!$OMP declare target(Mcon,Kcon,awei)
      end module svm1
   
      module backup2 
         use params
             integer eprofo,eprofn,energ
!$OMP declare target(eprofo,eprofn,energ)

      end module backup2


      module svm2
        use params
            integer acc_cut, nk_svm(3)
            integer ik_svm(3,15)   
            real :: dist_svm2(15)   
!$OMP declare target(acc_cut,nk_svm,ik_svm,dist_svm2)
      end module svm2

      module RES
        use params 
        real ::  ER3,er5,er6,er7
        integer :: Mcom(ndim)
        real :: Kcom(ndim,100)
!$OMP declare target(ER3,er5,er6,er7,Mcom,Kcom)
      end module RES 

      module hba
             real :: eh5a,Cr2a,acut_bb,acut_cc,acut_vv,acut_hh
!$OMP declare target(eh5a,Cr2a,acut_bb,acut_cc,acut_vv,acut_hh)
      end module hba

      module lengths
           integer :: Lch, Lch1, Lch2
!$OMP declare target(Lch,Lch1,Lch2)
      end module lengths

      module short1
             integer  JBIN(0:500) 
             real ::  acops(1500,16)
             real ::  bsr(1500,16)
!$OMP declare target (JBIN,bsr,acops)
      end module short1
  
      module echain1
             real :: ex(1500),ey(1500),ez(1500)
! !$OMP declare target(ex,ey,ez)
!$OMP declare target(ex,ey,ez)
      end module echain1
      
      module shortcom
              real ::  eh3,es4,es5,es6,es7,es7a,es7b,es7c
!$OMP declare target(eh3,es4,es5,es6,es7,es7a,es7b,es7c)
      end module shortcom

      module hbb
             real :: eh5b,Cr2b,bcut_bb,bcut_cc,bcut_vv,bcut_hh
!$OMP declare target (eh5b,Cr2b,bcut_bb,bcut_cc,bcut_vv,bcut_hh)
      end module hbb

      module chain1
             integer  x(1500),y(1500),z(1500)
             integer ica(0:1500)
! !$OMP declare target(x,y,z,ica)
!$OMP declare target(x,y,z,ica)
      end module chain1
      
      module seqe
             integer  seq(1500), sec(1500)
!$OMP declare target(seq,sec)
      end module seqe
      
      module hb
              real :: hbx(416,416), hby(416,416),hbz(416,416)
!$OMP declare target(hbx,hby,hbz)
      end module hb
      
      module bisec
             real :: cax(416,416), cay(416,416), caz(416,416)
!$OMP declare target(cax,cay,caz)
      end module bisec
      
      module sizea
             real :: ala(0:19,0:19),alm(0:19,0:19),alp(0:19,0:19)
!$OMP declare target(ala,alm,alp)
      end module sizea
      
      module pair 
             real :: apa(1500,1500),app(1500,1500),apm(1500,1500)
!$OMP declare target(apa,app,apm)
      end module pair
      
      module size
             real :: arla(0:19,0:19), arlm(0:19,0:19),arlp(0:19,0:19)
!$OMP declare target(arla,arlm,arlp)
      end module size
      
      module ehbenergy
             integer EHB1,EHB1b,EHB1c,EHB2,EHB3,EHB5,EHB6  !EHB4,EHB1a
!$OMP declare target(EHB1,EHB1b,EHB1c,EHB2,EHB3,EHB5,EHB6) !EHB4,EHB1a
      end module ehbenergy
      
      module ehbenergy1
             integer EHB5a,EHB5b
!$OMP declare target(EHB5a,EHB5b)
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
!$OMP declare target(ox,nx,ex_o,egx_o,ecx_o,ebx_o,etx_o,
!$OMP&                ex_n,egx_n,ecx_n,ebx_n,etx_n,oy,oz,oo,
!$OMP&                ny,nz,nn,ey_o,ez_o,egy_o,egz_o,
!$OMP&                ecy_o,ecz_o,eby_o,ebz_o,ety_o,etz_o,
!$OMP&                ey_n,ez_n,egy_n,egz_n,
!$OMP&                ecy_n,ecz_n,eby_n,ebz_n,ety_n,etz_n)
      end module tempArrays
 
      module par
             real :: apar(1500,1500)
!$OMP declare target(apar)
      end module par
     
      module expose
         use params
	 integer :: mp(20,ndim),area(ndim)
!$OMP declare target(mp,area)
      end module expose

      module CAcontact
c Merge CAcontact, CA8
	use params
        integer McomCA(ndim),KcomCA(ndim,100)
	integer McomCA8(ndim),KcomCA8(ndim,100)
	real :: aweigCA(ndim,ndim),aweigCA8(ndim,ndim)
!$OMP declare target(McomCA,KcomCA,McomCA8,KcomCA8,aweigCA,aweigCA8)
      end module CAcontact
 
      module sg
             real :: hx(416,416,0:19),hy(416,416,0:19),hz(416,416,0:19)
             real :: gx(416,416,0:19),gy(416,416,0:19),gz(416,416,0:19)
!$OMP declare target(gx,gy,gz,hx,hy,hz)
      end module sg
      
      module icgg
             integer EH6,icg(1500)
!$OMP declare target(icg,EH6)
      end module icgg
      
      module envir1
             real :: nop(1500),nom(1500),noa(1500),nhbn(1500)
!$OMP declare target(nop,nom,noa,nhbn)
      end module envir1
      
      module distres
             real :: er4,es3c
!$OMP declare target(er4,es3c)
      end module distres
      
      module  order
              integer sumct,sumcto,icnt,icnto
              real :: en2,en3,dord,acorder
!$OMP declare target(acorder,en2,sumct,sumcto,icnt,icnto,dord,en3) 
      end module order 
      
      module echain2
              real :: egx(1500),egy(1500),egz(1500)
!$OMP declare target (egx,egy,egz)
      end module echain2
      
      module echain4
             real :: ecx(1500),ecy(1500),ecz(1500)
!$OMP declare target (ecx,ecy,ecz)
      end module echain4
      
      module echain5
             real :: ebx(1500),eby(1500),ebz(1500)
!$OMP declare target (ebx,eby,ebz)
      end module echain5
     
      module fr
         use params
         real :: frga(ndim),frgb(ndim)
!$OMP declare target (frga,frgb)
      end module fr

      module freg
          use params
          real :: aweig(ndim,ndim)
!$OMP declare target (aweig)
      end module freg

      module echain6
             real :: etx(1500),ety(1500),etz(1500)
!$OMP declare target (etx,ety,etz)
      end module echain6
      
      module concutt
          real :: concut(0:19,0:19),concut2(0:19,0:19)
          integer concut_sc
!$OMP declare target(concut,concut2,concut_sc)
      end module concutt
      
      module pair1
          real :: eh2,eh1b,eh1c
!$OMP declare target(eh2,eh1b,eh1c)
      end module pair1

      module longdist
           use params
           integer MdisL(ndim),kdisL(ndim,500),distL(ndim,500)
!$OMP declare target(MdisL,kdisL,distL)
      end module longdist

!      COMMON/RCN/Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
      module RCN
        use params
        integer Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
!$OMP declare target(Mdis,kdis,dist,dev)
      end module RCN

      module one
        use lengths
             integer acrit
             real :: eoinp(0:19,0:100),es1,es2,contt,eonekd(0:19)
             real :: eonehw(0:19)
         
       DATA eonekd /-0.4, 1.8, -0.8, 2.5, 4.2, -0.7, 4.5, -1.6, 1.9
     $             , -3.5, -3.5, 3.8, -3.9, -3.5, -3.5, -4.5, -3.2, 2.8
     $             , -1.3, -0.9/
c            ^          ^     ^     !contradict with 'centro.comm'
c     read hydrophilic potential for Sg, positive for hydrophilic
c     residue--->
       DATA eonehw /0.0, -0.5, 0.3, -1.0, -1.5, -0.4, -1.8, 0.0, -1.3
     $             , 3.0, 0.2, -1.8, 3.0, 3.0, 0.2, 3.0, -0.5, -2.5
     $             , -2.3, -3.4 /
             ! real, PRIVATE :: eonekd(0:19)
!$OMP declare target(acrit,contt,eoinp,es2,es1,eonekd,eonehw)
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
!      common/hopp/eonehw(0:19)

c     read hydrophobic potential for Sg, positive for hydrophobic
c     residue--->
!      data eonekd /-0.4, 1.8, -0.8, 2.5, 4.2, -0.7, 4.5,
!     &     -1.6, 1.9,  -3.5, -3.5, 3.8,
!     &     -3.9, -3.5, -3.5, -4.5,
!     &     -3.2, 2.8, -1.3, -0.9/

!       DATA eonekd /-0.4, 1.8, -0.8, 2.5, 4.2, -0.7, 4.5, -1.6, 1.9 
!     $             ,-3.5, -3.5, 3.8, -3.9, -3.5, -3.5, -4.5, -3.2, 2.8
!     $             ,-1.3, -0.9/
c            ^          ^     ^     !contradict with 'centro.comm'
c     read hydrophilic potential for Sg, positive for hydrophilic
c     residue--->
!      DATA eonehw /0.0, -0.5, 0.3, -1.0, -1.5, -0.4, -1.8, 0.0, -1.3 
!     $            ,3.0, 0.2, -1.8, 3.0, 3.0, 0.2, 3.0, -0.5, -2.5
!     $            ,-2.3, -3.4 /

c     expected gyration radius:
      acrit=2.2*exp(0.38*alog(float(Lch)))/0.87 !gyrat-radius~2.2*l^0.38
*     Defination of gyration-radius: acrit=sqrt(<(r-r0)^2>)

c^^^^^^^^^^^^^^^^^ read centersymmetric potential finished ^^^^^^^^^^^
      end subroutine read_centro
! !$OMP declare target(acrit,contt,eoinp,es2,es1,eonekd)   
      end module one 
      
