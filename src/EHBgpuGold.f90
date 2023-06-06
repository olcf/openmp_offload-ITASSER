      FUNCTION EHB(jjjj,kkkk,ISTAT)
	use params
      use chainm
      !use openacc
      use chain1
      use echain1 
      use short1
      use lengths
      use ENERGY
      use hba
      use hbb
      use hb
      use bisec
      use sg
      use pair
      use size
      use sizea
      use icgg
      use distres
!      use ehbenergy
!      use ehbenergy1
      use par
      use seqe
      use concutt
      use pair1
      use shortcom
      use one
      use echain2
      use echain4
      use echain5 
      use envir1
      use order
      use aax_mod
      use aay_mod
      use aaz_mod
      use avv_mod
      use ei5_mod
      use energyHBa_mod
      use energyHBb_mod
     
      IMPLICIT INTEGER(I-Z)
      EHB1=0                    !+1/r of Ca-SC
      EHB1a=0                   !+1/r for non-parallel contact of Ca-Ca
      EHB1b=0                   !excluded volume of SC-SC 
      EHB1c=0                   !pair potential of SC-SC, paor1.dat+ par.dat

      EHB2=0                    !quarsi3 for SC-SC, pair3.dat+quarsi3.comm
      EHB3=0                    !enhance good-piece contacts
      EHB4=0                    !-1/r for parallel contact of Ca-Ca

      EHB5a=0                    !H-bond energy for alpha-type
      EHB5b=0                    !H-bond energy for beta-type

      if(ISTAT.gt.0)THEN        !Enew
         ICNT=ICNTO
         SUMCT=SUMCTO
      ENDIF

c     coupling of secondary structure with pairwise interactions
c     included - thus expanded range
      jj=jjjj-1
      if(jj.lt.1)jj=1
      kk=kkkk+1
      if(kk.gt.Lch)kk=Lch
      ICNT1=ICNT
      SUMCNT1=SUMCNT
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*****************************************************************
********** circle for movement involved window ******************
*****************************************************************
! !$acc data copy (ehb4,ehb5a,ehb5b,ehb3,ehb2,
! !$acc&              ehb1c,ehb1b,ehb1a,ehb1,SUMCT1,ICNT1)
!$acc parallel loop gang default (present) ! reduction(+:EHB1,EHB1a,EHB1b,EHB1c,EHB2,EHB3,EHB4,
! !$acc&                        EHB5,EHB5a,EHB5b,ICNT1,SUMCT1)
! !$acc&    default(present) 
      DO 1002 k=jj,kk
         kseq=seq(k)
         if(mv(k).gt.0)then     !moveable point
            axk=x(k)            !Ca
            ayk=y(k)
            azk=z(k)
            nv1=ica(k-1)
            nv2=ica(k)
            agxk=axk+GX(nv1,nv2,kseq) !Cg_k
            agyk=ayk+GY(nv1,nv2,kseq)
            agzk=azk+GZ(nv1,nv2,kseq)
            bxk=HBX(nv1,nv2)    !H-bond direction
            byk=HBY(nv1,nv2)
            bzk=HBZ(nv1,nv2)
            cxk=CAX(nv1,nv2)    !vertical vector
            cyk=CAY(nv1,nv2)
            czk=CAZ(nv1,nv2)	
         else
            axk=ex(k)           !Ca
            ayk=ey(k)
            azk=ez(k)
            agxk=egx(k)         !SG
            agyk=egy(k)
            agzk=egz(k)
            bxk=ebx(k)          !Hb
            byk=eby(k)
            bzk=ebz(k)
            cxk=ecx(k)          !cc
            cyk=ecy(k)
            czk=ecz(k)
         endif
         km2=k-2
         kp2=k+2
         if(km2.ge.1.and.kp2.le.Lch)then
            xxx=nint((aax(km2)-aax(kp2))**2+
     $           (aay(km2)-aay(kp2))**2+(aaz(km2)-aaz(kp2))**2) !real dist
            if(xxx.gt.500)xxx=500
            ek5=acops(km2,jbin(xxx)) !k2-residue, jbin. es<0.
         else
            ek5=0
         endif
*****************************************************************
********************** start i-cicle ****************************
*****************************************************************
!$acc loop gang reduction(+:EHB1,EHB1a,EHB1b,EHB1c,EHB2,EHB3,EHB4,
!$acc&                        EHB5,EHB5a,EHB5b,ICNT1,SUMCT1)
!acc&  private(i,fact,axki,ayki,azki,dist,bxi,byi,bzi,
!acc&          axi,ayi,azi,idist)
         do 1001 i=1,Lch
            iend=max(k+1,kk)
            if(i.ge.k-1.and.i.le.iend)goto 1001 !to avoid repeat
            if(mv(i).gt.0)then
               axi=x(i)
               ayi=y(i)
               azi=z(i)
            else
               axi=ex(i)
               ayi=ey(i)
               azi=ez(i)
            endif
            axki=(axk-axi)      !r_k-r_i
            ayki=(ayk-ayi)
            azki=(azk-azi)
            Cr2=axki**2+ayki**2+azki**2+0.000001 !(r_k-r_i)**2
            if(Cr2.lt.120)then  !120-->9.53, there is pairwise interaction
               idist=iabs(i-k)
               iseq=seq(i)
               if(mv(i).gt.0)then
                  nv1=ica(i-1)
                  nv2=ica(i)
                  agxi=axi+GX(nv1,nv2,iseq) !Cg_k
                  agyi=ayi+GY(nv1,nv2,iseq)
                  agzi=azi+GZ(nv1,nv2,iseq)
                  bxi=HBX(nv1,nv2) !H-bond direction
                  byi=HBY(nv1,nv2)
                  bzi=HBZ(nv1,nv2)
                  cxi=CAX(nv1,nv2) !outside 
                  cyi=CAY(nv1,nv2)
                  czi=CAZ(nv1,nv2)	
               else
                  agxi=egx(i)   !SG
                  agyi=egy(i)
                  agzi=egz(i)
                  bxi=ebx(i)    !H-bond
                  byi=eby(i)
                  bzi=ebz(i)
                  cxi=ecx(i)    !bisector vector
                  cyi=ecy(i)
                  czi=ecz(i)
               endif

c     1/r excluded for Ca-SC pair ----------------->
               if(kseq.gt.0) then !not GLY, we have SG
                  aks=(agxk-axi)**2+(agyk-ayi)**2+(agzk-azi)**2+.0001 !SG_kCa_i
                  if(aks.lt.36) then !36-->5.22A
                     if(aks.lt.13.0) aks=13.0 !13-->3.17A
                     EHB1=EHB1+((13.0/aks)-0.5)
                  endif
               endif
               if(iseq.gt.0) then
                  aks=(axk-agxi)**2+(ayk-agyi)**2+(azk-agzi)**2+.0001 !Ca_kSG_i
                  if(aks.lt.36.0) then
                     if(aks.lt.13.0) aks=13.0
                     EHB1=EHB1+((13.0/aks)-0.5)
                  endif
               endif

c     pair-wise potential of SC-SC---------->
               Gr2=(agxi-agxk)**2+(agyi-agyk)**2+(agzi-agzk)**2 !SG_i,SG_k
               if(Gr2.lt.concut2(iseq,kseq))then
                  EHB1c=EHB1c+apar(i,k)
               endif

c     quarsi3 for SC-SC pair, 1/r for Ca-Ca ----------------->
               cc=cxi*cxk+cyi*cyk+czi*czk !c_i*c_k
               IF(cc.gt.0.5)THEN !c_i//c_k
                  if(Cr2.lt.60)EHB4=EHB4-(30/(max(30.,Cr2))-0.5) !6.74A
                  if(Gr2.lt.alp(iseq,kseq))then
                     EHB3=EHB3-ek5*ei5(i,idist)
!$acc atomic update
                     NOP(k)=NOP(k)+ISTAT
!$acc atomic update
                     NOP(i)=NOP(i)+ISTAT
                     ICNT1=ICNT1+istat*idist
                     SUMCT1=SUMCT1+ISTAT
                     if(Gr2.gt.arlp(iseq,kseq))THEN
                        EHB2=EHB2+app(i,k) !quarsi3
                     else
                        EHB1b=EHB1b+1 !excluded volume
                     endif
                  endif
               ELSE
                  if(Cr2.lt.33)EHB1a=EHB1a+(16.0/Cr2-0.5) !5A
                  IF(cc.lt.-0.5) THEN !antiparallel pair-wise of (i,k)
                     if(Gr2.lt.ala(iseq,kseq))then
                        EHB3=EHB3-ek5*ei5(i,idist)
!$acc atomic update
                        NOA(k)=NOA(k)+ISTAT
!$acc atomic update
                        NOA(i)=NOA(i)+ISTAT
                        ICNT1=ICNT1+istat*idist
                        SUMCT1=SUMCT1+ISTAT
                        if(Gr2.gt.arla(iseq,kseq))THEN
                           EHB2=EHB2+apa(i,k)
                        else
                           EHB1b=EHB1b+1 !excluded volume
                        endif
                     endif
                  ELSE          !neither parallel nor antiparallel
                     if(Gr2.lt.alm(iseq,kseq)) then
                        EHB3=EHB3-ek5*ei5(i,idist)
!$acc atomic update
                        NOM(k)=NOM(k)+ISTAT
!$acc atomic update
                        NOM(i)=NOM(i)+ISTAT
                        ICNT1=ICNT1+istat*idist !distance of pairs
                        SUMCT1=SUMCT1+ISTAT !number of contact pairs
                        if(Gr2.gt.arlm(iseq,kseq))THEN
                           EHB2=EHB2+apm(i,k)
                        else
                           EHB1b=EHB1b+1 !excluded volume
                        endif
                     endif
                  ENDIF
               ENDIF
! !$acc end data
c^^^^^^^^^^^^^ interaction of pair-wise finished ^^^^^^^^^^^^^^^^^

ccccccccccccccccccccccccc Hydrogen-bond cccccccccccccccccccccccccc
***** 1, cc; 2, bb; 3, distance; 4, vij; 5, sec(i).
c               if(idist.eq.3)then
c                  alpha-H-bond
c               elseif(idist.lt.20)then
c                  anti-parallel-sheet-H-bond
c               else
c                  if(bb<0)then
c                     antiparallel sheet-H-bond
c                  else
c                     parallel sheet-H-bond
c                  endif
c               endif
*HHHHH***************** alpha ******************************************
***   bb   bb_an  cc    cc_an vv   vv_an hbl r2i  r2k  Cr2  Cr2m  Hdis
***   0.82  34.9  0.42  65.3  0.52  58.8 5.7 3.29 3.53 35.9 31.62 4.19 (alpha)
***   0.81  35.3  0.38  67.3  0.48  61.0 5.7 3.55 3.39 35.5 30.45 4.10 (alpha1)
***   0.85  31.9  0.50  59.5  0.58  54.0 5.7 3.03 2.92 34.1 30.29 4.13 (alpha2)
************************************************************************
*********************** beta *******************************************
***   bb    bb_an   cc  cc_an vv12 vv12_a vv12 vv12_a hbl r2i r2k Cr2 Cr2m
***   -0.81 144.8 0.79  36.7 -0.94 160.6 -0.95 162.4 5.3 8.38 1.37 30 30 (r)
***   -0.84 149.1 0.83  31.8 -0.92 161.1 -0.94 162.7 5.3 3.91 4.43 33 31 (r1)
***   -0.64 135.0 0.46  61.8 -0.75 139.6 -0.74 137.6 5.1 7.81 6.55 33 31 (r2)
***    0.98   8.6 0.92  22.2  0.92  22.5  0.93  18.8 5.3 2.49 0.99 31 30 (p)
************************************************************************
               if(i.eq.1.or.i.eq.Lch.or.k.eq.1.or.k.eq.Lch)goto 1003
               if(idist.eq.3)then !!!!!!!!!!!!!!!!!!!!!!!!!!!
***   alpha-helix:
                  if(sec(k).ne.4.and.sec(i).ne.4)then
                  if(Cr2.lt.Cr2a.and.cc.gt.acut_cc)then
                  bb=bxi*bxk+byi*byk+bzi*bzk !cos(Hbi.Hbk)
                  if(bb.gt.acut_bb)then
                  av11=avv(aax(i-1),aay(i-1),aaz(i-1),axi,ayi,azi,
     $                    aax(k-1),aay(k-1),aaz(k-1),axk,ayk,azk)
                  if(av11.gt.acut_vv)then
                  av22=avv(axi,ayi,azi,aax(i+1),aay(i+1),aaz(i+1),
     $                    axk,ayk,azk,aax(k+1),aay(k+1),aaz(k+1))
                  if(av22.gt.acut_vv)then
***   ->
                     fact=(1-abs(cc-0.4))*(1-abs(bb-0.815))
                     EHB5a=EHB5a+energyHBa(i,k,fact,
     $                    axki,ayki,azki,bxi,byi,bzi,bxk,byk,bzk)
                  endif
                  endif
                  endif
                  endif
                  endif
               elseif(idist.gt.4.and.idist.lt.20)then !!!!!!!!!!!!!!!!!!!!!!!
***   antiparallel-sheet
                  if(sec(k).ne.2.and.sec(i).ne.2)then
                  if(Cr2.lt.Cr2b.and.cc.gt.bcut_cc)then
                  bb=bxi*bxk+byi*byk+bzi*bzk !cos(Hbi.Hbk)
                  if(bb.lt.-bcut_bb)then !antiparallel
                  av12=avv(aax(i-1),aay(i-1),aaz(i-1),axi,ayi,azi,
     $                    axk,ayk,azk,aax(k+1),aay(k+1),aaz(k+1))
                  if(av12.lt.-bcut_vv)then
                  av21=avv(axi,ayi,azi,aax(i+1),aay(i+1),aaz(i+1),
     $                    aax(k-1),aay(k-1),aaz(k-1),axk,ayk,azk)
                  if(av21.lt.-bcut_vv)then
***   ->
                     fact=abs(bb)*cc !bb->-1,cc->1
                     EHB5b=EHB5b+energyHBb(i,k,fact,
     $                    axki,ayki,azki,bxi,byi,bzi,bxk,byk,bzk)
                  endif
                  endif
                  endif
                  endif
                  endif
               elseif(idist.ge.20)then !!!!!!!!!!!!!!!!!!!!!!!!
                  if(sec(k).ne.2.and.sec(i).ne.2)then
                  if(Cr2.lt.Cr2b.and.cc.gt.bcut_cc)then
                  bb=bxi*bxk+byi*byk+bzi*bzk !cos(Hbi.Hbk) 
                  if(bb.lt.-bcut_bb)then !antiparallel
***   antiparallel-sheet:
                     av12=avv(aax(i-1),aay(i-1),aaz(i-1),axi,ayi,azi,
     $                    axk,ayk,azk,aax(k+1),aay(k+1),aaz(k+1))
                     if(av12.lt.-bcut_vv)then
                     av21=avv(axi,ayi,azi,aax(i+1),aay(i+1),aaz(i+1),
     $                       aax(k-1),aay(k-1),aaz(k-1),axk,ayk,azk)
                     if(av21.lt.-bcut_vv)then
***   ->
                        fact=abs(bb)*cc !bb->1,cc->1
                        EHB5b=EHB5b+energyHBb(i,k,fact,
     $                       axki,ayki,azki,bxi,byi,bzi,bxk,byk,bzk)
                     endif
                     endif
                  elseif(bb.gt.bcut_bb)then !bb>0, parallel H-bond
***   parallel-sheet:
                     av11=avv(aax(i-1),aay(i-1),aaz(i-1),axi,ayi,azi,
     $                    aax(k-1),aay(k-1),aaz(k-1),axk,ayk,azk)
                     if(av11.gt.bcut_vv)then
                     av22=avv(axi,ayi,azi,aax(i+1),aay(i+1),aaz(i+1),
     $                       axk,ayk,azk,aax(k+1),aay(k+1),aaz(k+1))
                     if(av22.gt.bcut_vv)then
***   ->
                        fact=bb*cc !bb->1,cc->1
                        EHB5b=EHB5b+energyHBb(i,k,fact,
     $                       axki,ayki,azki,bxi,byi,bzi,bxk,byk,bzk)
                     endif
                     endif
                  endif
                  endif
                  endif
               endif            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 1003          continue
c^^^^^^^^^^^^^H-bond energy finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            endif               !Ca-Ca<120
! !$acc update host(EHB1,EHB1a,EHB1b,EHB1c,EHB2,EHB3,EHB4,
! !$acc&                        EHB5,EHB5a,EHB5b,ICNT1,SUMCT1)
 1001    continue               !i -> [1,Lch]
 1002 continue                  !k -> [jj,kk]
! this 1!$acc end data
c     ICNT/SUMCT: average distance of residues in each contact-pairs.
c     ICNT/SUMCT are calculated when call Enew
c     ICNTO/SUMCTO are not changed when call Eold
! !$acc update host(EHB1,EHB1a,EHB1b,EHB1c,EHB2,EHB3,EHB4,
! !$acc&                        EHB5,EHB5a,EHB5b,ICNT1,SUMCT1)
      if(istat.lt.0) then       !Eold
         b=abs(ICNT1/(0.000001+float(SUMCT1))-acorder)
         a=abs(ICNTO/(0.000001+float(SUMCTO))-acorder)
         d=abs(float(SUMCT1)-contt) !deviation of contact number on new conform
         c=abs(float(SUMCTO)-contt) !deviation of contact number on new confor
         dord=en2*(b-a)+en3*(d-c) !not included in EHB
      endif
      ICNT = INCT1     
      SUMCT = SUMCT1
      
      EHB
     $     =eh1*EHB1            !+1/r of Ca-SC
     $     +eh1a*EHB1a          !+1/r for non-parallel of Ca-Ca
     $     +eh1b*EHB1b          !excluded volumn of SC-SC
     $     +eh1c*EHB1c          !pair-wise potential of SC-SC
     $     +eh2*EHB2            !quarsi3 for SC-SC
     $     +eh3*EHB3            !enhance good piece
     $     +eh4*EHB4            !-1/r for parallel contact of Ca-Ca
     $     +eh5a*EHB5a          !H-bond energy (alpha)
     $     +eh5b*EHB5b          !H-bond energy (beta)

c      write(*,*)EHB,eh1c,EHB1c,jjjj,kkkk
c ^^^^^^^^^^ EHB finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      RETURN
      END
