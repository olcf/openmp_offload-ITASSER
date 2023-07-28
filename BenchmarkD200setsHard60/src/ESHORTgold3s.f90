      FUNCTION ESHORT(iiii,jjjj,ISTAT)
      use params
      use chainm
      use chain1
      use echain1
      use short1
      use lengths
      use ENERGY
      use hb
      use echain5 
      use bisec
      use distres
      use concutt
      use sg
      use shortcom
      use echain2
      use echain4
      use seqe
      use one
      use longdist
      use RCN
      use echain6
      use short
      use RES
      use svm1 
      use svm2
      use fr
      use freg
      use aax_mod
      use aay_mod
      use aaz_mod
      use expose
      use stick
      use eigen
      use CAcontact 
      use shape
      use trackn 
      IMPLICIT INTEGER(I-Z)
c      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
C FIX afs
!      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
c      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/three/angle(nvec,nvec)
c      COMMON/RCN1/er1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/msichores/msicho
      common/eshortenergy1/ESHORT1,ESHORT2,ESHORT3,ESHORT4,ESHORT11
      common/eshortenergy2/ESHORT4a,ESHORT5,ESHORT5a,ESHORT5b,ESHORT5c
      common/eshortenergy3/ESHORT6,ESHORT7,ESHORT8,ESHORT9,ESHORT10
      common/eshortenergy4/ESHORT12
c      common/hopp/eonehw(0:19)
c      common/hom1/asr2(ndim,4),asr3(ndim,4),asr4(ndim,14),asr5(ndim,8)
c      common/hom2/ibb2(0:999),ibb3(0:999),ibb4(-999:999),ibb5(0:999)
c      common/CAcontact/McomCA(ndim),KcomCA(ndim,100),aweigCA(ndim,ndim)
c      common/CA8/McomCA8(ndim),KcomCA8(ndim,100),aweigCA8(ndim,ndim)
      common/CAcontact2/dist_CA_cut8
      common/zscore/izscore
      common/CAcontact1/dist_CA_cut
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      common/rmsdder/i_chunk(ndim),ex0(ndim),ey0(ndim),ez0(ndim)
      common/E_defo/i_E_defo
      common/initialinput/switch,k_cycle,k_phot,N_ann
      common/rmsd_defo/armsd,armsd_sum(100),N_rmsd(100)
!      common/expose/mp(20,ndim),area(ndim)
      common/res2/er14,er15,er16,er17
      common/sidechain/bgx(ndim),bgy(ndim),bgz(ndim)
      dimension cgx(ndim),cgy(ndim),cgz(ndim)
      common/center/cex,cey,cez
c      common/eigen/AA(3,3),EE(3),HH(3,3)
c      common/stick1/nstick,astick,nrmsd,ermsd
c      common/stick2/iq(ndim,nrep)
c      common/stick3/ax00(ndim,nrep),ay00(ndim,nrep),az00(ndim,nrep)
c      common/stick6/bx00(ndim),by00(ndim),bz00(ndim),armsd_min0
c      common/stick7/itemp0,icycle,icycle0
!      common/trackn/n_tem(100)
      common/temperature/itemp,atemp
      common/svm3/er21,er22,er23
      
cccc   RMSD:
      double precision r_1(3,ndim),r_2(3,ndim),r_3(3,ndim),w(ndim)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /ndim*1.0/
ccc   
! !$acc update device(ica)

!$acc data present(distl(:,:), kdisl(:,:),mdisl(:),n_tem(:) 
!$acc&  ,dist(:,:),mdis(:),kdis(:,:), dev(:,:)
!$acc&  ,hbz(:,:),hby(:,:),hbx(:,:)
!$acc&  ,bsr(:,:),jbin(:),Kcon(:,:,:),Mcon(:,:)
!$acc$  ,Mcom(:),Kcom(:,:),awei(:,:,:) 
!$acc&  ,hx(:,:,:), hy(:,:,:), hz(:,:,:)
!$acc&  ,gx(:,:,:), gy(:,:,:), gz(:,:,:)
!$acc&  ,x(:), y(:), z(:), acc_cut
!$acc&  ,dist_svm2(:),nk_svm(:),ik_svm(:,:)
!$acc&  ,ex(:), ey(:), ez(:),IBIN,asr
!$acc&  ,egx(:), egy(:), egz(:)
!$acc&  ,ebx(:), eby(:), ebz(:)
!$acc&  ,ecx(:), ecy(:), ecz(:)
!$acc&  ,etx(:), ety(:), etz(:),mcomca(:)
!$acc&  ,concut(:,:),concut2(:,:),afs(:),afsn(:)
!$acc&  ,area(:),ee(:),hh(:,:),ax00,ay00,az00) 

c      ESHORT=0.0
      ESHORT2=0                 !bury potential for SC

      ESHORT3=0  !distance restrain from both threading (for C_a)
      ESHORT4=0  !contact restrain from threading (for SC)
      ESHORT4a=0 !deviation of contact restrain

      ESHORT5=0  !bias2,3: v(i)-v(i+4) anti/parallel; c(i)-c(i+2) anit/paralel
      ESHORT5a=0 !crumpling
      ESHORT5b=0 !bias4 to predicted alpha/beta structure.
      ESHORT5c=0 !bias1 to possible alpha/beta structure. 

      ESHORT6=0  !correlation of E13 of Ca
      ESHORT7=0  !correlation of E14, from both common and 2th specific data
      ESHORT8=0  !correlation of E15, from both common and 2th specific data

      ESHORT9=0  !contact restraints of CA
      ESHORT10=0 !Long-range distance restraints of CA

      ESHORT11=0 !RMSD deviation
      ESHORT12=0 !distance from template

      ESHORT13=0 !rmsd to template
c      ESHORT14=0 !contact restraints of CA at 8A
      
c      ESHORT17=0                !deviation to comb8CA.dat
      ESHORT21=0                !CA_SVM
      ESHORT22=0                !CB_SVM
      ESHORT23=0                !SG_SVM
      
      if(istat.lt.0)THEN       !Eold 
         diold=0.0
         coold=0.0
      else                      !Enew
         conew=0.0   !penality for incorrect contact (from contact restrain)
         dinew=0.0   !penality for incorrect distance (from distant restain)
      endif

!$acc parallel loop gang async(1) reduction(+:ESHORT2) default(present) 
      do 1 i=iiii,jjjj
         iseq=seq(i)
         if(mv(i).gt.0)then
            axi=x(i)
            ayi=y(i)
            azi=z(i)
            nv1=ica(i-1)
            nv2=ica(i)
            agxi=axi+GX(nv1,nv2,iseq) !Cg_k
            agyi=ayi+GY(nv1,nv2,iseq)
            agzi=azi+GZ(nv1,nv2,iseq)
         else
            axi=ex(i)           !Ca
            ayi=ey(i)
            azi=ez(i)
            agxi=egx(i)         !SG
            agyi=egy(i)
            agzi=egz(i)
         endif
***** Bury energy ocf SG --------------------->
         gxn=HH(1,1)*(agxi-cex)+HH(1,2)*(agyi-cey)+HH(1,3)*(agzi-cez)
         gyn=HH(2,1)*(agxi-cex)+HH(2,2)*(agyi-cey)+HH(2,3)*(agzi-cez)
         gzn=HH(3,1)*(agxi-cex)+HH(3,2)*(agyi-cey)+HH(3,3)*(agzi-cez)
         fff=gxn**2/EE(1)+gyn**2/EE(2)+gzn**2/EE(3) !=5 for ellipsoid surphase
         if(fff.lt.2.5)then
           aaa=fff-2.5
            if(aaa.lt.-1)aaa=-1
            ESHORT2=ESHORT2-aaa*area(i) !area(i)=-1, buried; +1, exposed
         endif
c^^^^^^^^^^^^^^^^^^ bury energy finished ^^^^^^^^^^^^^^^^^^

c     deepth factor: afs=1, r<r0; [0.5,1], r0<r<2r0; 0.5, r>2r0
         ar=sqrt((axi-cex)**2+(ayi-cey)**2+(azi-cez)**2+0.01) !r
         afs(i)=acrit/ar
         if(afs(i).gt.1)afs(i)=1
         if(afs(i).lt.0.5)afs(i)=0.5
         
****  restrain from 'dist.dat'---------------------->
         N_dis=Mdis(i)          !number of distance restraints on 'i'
!$acc loop vector  reduction(+:ESHORT3) 
         do k=1,N_dis
            j=kdis(i,k)         !i-j restraints
            if(j.lt.i.or.j.gt.jjjj) then !to avoid repeat
               dij2=(axi-aax(j))**2+(ayi-aay(j))**2+(azi-aaz(j))**2
               dij=sqrt(dij2)
               err=abs(dij-dist(i,k)) !dist: predicted dis
               if(err.gt.dev(i,k)) then !dev: deviation for arca
                  ESHORT3=ESHORT3+1
               endif
            endif
         enddo
c^^^^^^^^^^^^^^^^^ distant restrain finished ^^^^^^^^^^^^^^^^^^^^^^^^^

****  long-range dist-restrain from 'distL.dat'---------------------->
         N_disL=MdisL(i)        !number of long-range dist restraints on 'i'
!$acc loop vector reduction(+:ESHORT10) 
         do k=1,N_disL
            j=kdisL(i,k)        !i-j restraints
            if(j.lt.i.or.j.gt.jjjj) then !to avoid repeat
               dij2=(axi-aax(j))**2+(ayi-aay(j))**2+(azi-aaz(j))**2
               dij=sqrt(dij2)
               err=abs(dij-distL(i,k))
               if(err.lt.1)err=1
               ESHORT10=ESHORT10-1/err
            endif
         enddo
c^^^^^^^^^^^^^^^^^distant restrain finished ^^^^^^^^^^^^^^^^^^^^^^^^^
 1    continue
!$acc parallel loop gang default(present) async(1)  !reduction(+:ESHORT2) 
ccc   Contact restrain from 'comb.dat' -------------------------------->
       do 11 i=iiii,jjjj
          N_com=Mcom(i)         !number of contacts on i
          mk=0
          do kk=1,nk_svm(3)
             mk=mk+Mcon(ik_svm(3,kk),i)
          enddo
          if(N_com.ge.1.or.mk.ge.1)then
             iseq=seq(i)
             if(mv(i).gt.0)then
                agxi=x(i)+GX(ica(i-1),ica(i),iseq) !SG
                agyi=y(i)+GY(ica(i-1),ica(i),iseq) !SG
                agzi=z(i)+GZ(ica(i-1),ica(i),iseq) !SG
             else
                agxi=egx(i)     !SG
                agyi=egy(i)
                agzi=egz(i)
             endif
             
ccc   'comb.dat':
!$acc loop vector reduction(+:ESHORT4,coold,conew) 
             do k=1,N_com
                j=Kcom(i,k)     !k'th contact with i
                if(j.lt.i.OR.j.gt.jjjj)then
                   jseq=seq(j)
                   if(mv(j).gt.0)then
                      agxj=x(j)+gx(ica(j-1),ica(j),jseq)
                      agyj=y(j)+gy(ica(j-1),ica(j),jseq)
                      agzj=z(j)+gz(ica(j-1),ica(j),jseq)
                   else
                      agxj=egx(j)
                      agyj=egy(j)
                      agzj=egz(j)
                   endif
                   aij2=(agxi-agxj)**2+(agyi-agyj)**2+(agzi-agzj)**2
                   if(aij2.gt.concut2(iseq,jseq))then
                      ESHORT4=ESHORT4+2*aweig(i,j) !no contact
                      if(istat.lt.0)then
                         coold=coold+sqrt(aij2)-concut(iseq,jseq) !penality
                      else
                         conew=conew+sqrt(aij2)-concut(iseq,jseq) !panelity
                      endif
                   endif
                endif
             enddo
             
ccc   'svmseqSG678.dat':
!$acc loop worker ! reduction(+:ESHORT23)  !it was vector/ change 2
!worker and back to vector
             do kk=1,nk_svm(3)  !number contact files
                ik=ik_svm(3,kk) !original file order
!$acc loop vector  reduction(+:ESHORT23)  
                do k=1,Mcon(ik,i) !number of contacts in ik'th file
                   j=Kcon(ik,i,k) !k'th contact with i in ik'th file
                   if(j.lt.i.OR.j.gt.jjjj)then
                      jseq=seq(j)
                      if(mv(j).gt.0)then
                         agxj=x(j)+gx(ica(j-1),ica(j),jseq)
                         agyj=y(j)+gy(ica(j-1),ica(j),jseq)
                         agzj=z(j)+gz(ica(j-1),ica(j),jseq)
                      else
                         agxj=egx(j)
                         agyj=egy(j)
                         agzj=egz(j)
                      endif
                      aij2=(agxi-agxj)**2+(agyi-agyj)**2+(agzi-agzj)**2
                      if(aij2.lt.dist_svm2(ik))then
                         ESHORT23=ESHORT23-awei(ik,i,j)
                      endif
                   endif
                enddo
             enddo
! !$acc end kernels
          endif
 11    continue

!$acc parallel loop gang default(present) async(1)
ccc   Contact restrain from 'svmseqCB678.dat' ---------------------->
       do 12 i=iiii,jjjj
          mk=0
!$acc loop vector  reduction(+:mk) ! it was vector
          do kk=1,nk_svm(2)
             mk=mk+Mcon(ik_svm(2,kk),i)
          enddo
          if(mk.ge.1)then
             iseq=seq(i)
             if(mv(i).gt.0)then
                abxi=x(i)+HX(ica(i-1),ica(i),iseq) !CB
                abyi=y(i)+HY(ica(i-1),ica(i),iseq) !CB
                abzi=z(i)+HZ(ica(i-1),ica(i),iseq) !CB
             else
                abxi=etx(i)     !CB
                abyi=ety(i)
                abzi=etz(i)
             endif
             
ccc   'svmseqCB678.dat':
!$acc loop worker ! reduction(+:ESHORT22)
             do kk=1,nk_svm(2)  !number contact files
                ik=ik_svm(2,kk) !original file order
!$acc loop vector reduction(+:ESHORT22)
                do k=1,Mcon(ik,i) !number of contacts in ik'th file
                   j=Kcon(ik,i,k) !k'th contact with i in ik'th file
                   if(j.lt.i.OR.j.gt.jjjj)then
                      jseq=seq(j)
                      if(mv(j).gt.0)then
                         abxj=x(j)+hx(ica(j-1),ica(j),jseq)
                         abyj=y(j)+hy(ica(j-1),ica(j),jseq)
                         abzj=z(j)+hz(ica(j-1),ica(j),jseq)
                      else
                         abxj=etx(j)
                         abyj=ety(j)
                         abzj=etz(j)
                      endif
                      aij2=(abxi-abxj)**2+(abyi-abyj)**2+(abzi-abzj)**2
                      if(aij2.lt.dist_svm2(ik))then
                         ESHORT22=ESHORT22-awei(ik,i,j)
                      endif
                   endif
                enddo
             enddo
          endif
 12    continue

!$acc parallel loop gang async(1) default(present)
ccc   Contact restrain from 'combCA.dat' ------------------------->
c      dist_CA_cut=(6.0/0.87)**2
      do 33 i=iiii,jjjj
         N_com=McomCA(i)        !number of contacts on i
         mk=0
!$acc loop vector  reduction(+:mk) 
         do kk=1,nk_svm(1)
            mk=mk+Mcon(ik_svm(1,kk),i)
         enddo
         if(N_com.ge.1.or.mk.ge.1)then
            if(mv(i).gt.0)then
               axi=x(i)
               ayi=y(i)
               azi=z(i)
            else
               axi=ex(i)
               ayi=ey(i)
               azi=ez(i)
            endif

ccc   'combCA.dat':
!$acc loop vector reduction(+:ESHORT9) 
            do k=1,N_com
               j=KcomCA(i,k)    !k'th contact with i
               if(j.lt.i.OR.j.gt.jjjj)then
                  aij2=(axi-aax(j))**2+(ayi-aay(j))**2+(azi-aaz(j))**2
                  if(aij2.lt.dist_CA_cut)then
                     ESHORT9=ESHORT9-aweigCA(i,j)
                  endif
               endif
            enddo
            
!$acc loop worker ! reduction(+:ESHORT21)
ccc   'svmseqCA678.dat':
            do kk=1,nk_svm(1)
               ik=ik_svm(1,kk)
!$acc loop vector  reduction(+:ESHORT21)
               do k=1,Mcon(ik,i)
                  j=Kcon(ik,i,k) !k'th contact with i
                  if(j.lt.i.OR.j.gt.jjjj)then
                     aij2=(axi-aax(j))**2+(ayi-aay(j))**2+
     &                    (azi-aaz(j))**2
                     if(aij2.lt.dist_svm2(ik))then
                        ESHORT21=ESHORT21-awei(ik,i,j)
                     endif
                  endif
               enddo
            enddo
            
c            goto 123
c            do k=1,N_com8
c               j=KcomCA8(i,k)   !k'th contact with i
c               if(j.lt.i.OR.j.gt.jjjj)then
c                  aij2=(axi-aax(j))**2+(ayi-aay(j))**2+(azi-aaz(j))**2
c                  err=aij2
c                  if(aij2.lt.dist_CA_cut8)then
c                     ESHORT14=ESHORT14-aweigCA8(i,j)
c                     err=dist_CA_cut8
c                  endif
c                  err=sqrt(err)
c                  if(err.le.dist_CA_cut8*1.5)then
c                     ESHORT17=ESHORT17-aweigCA8(i,j)/err
c                  endif
c               endif
c            enddo
c 123        continue
         endif
 33   continue
! !!!!!$acc end data
c^^^^^^^^^^^^^^^^^^ CAcontact restrains finished ^^^^^^^^^^^^^^^^^^^^^
      
*********E13 --------------------------------------->
      i1=max(iiii-1,1)
      i2=min(jjjj-1,Lch-2)
!$acc parallel loop gang reduction(+:ESHORT6) default(present) async(1)
!I put a gang  here
      do i=i1,i2
         ar13=(aax(i)-aax(i+2))**2+(aay(i)-aay(i+2))**2+
     $        (aaz(i)-aaz(i+2))**2
         if(ar13.lt.48) then    !6.03A
            ESHORT6=ESHORT6+csr(i,1)
         else
            ESHORT6=ESHORT6+csr(i,2)
         endif
      enddo
c^^^^^^^^^^^^^^ E_13 finished ^^^^^^^^^^^^^^^^^^^^^^^^^
      
*********E14, E15, proteinlike bias1 to possible alpha/beta --------------->
      i1=max(iiii-4,1)
      i2=min(jjjj,Lch-4)
! !$acc  data copyin(asr(i1:i2,:),ibin(:),hbz(:,:),sec(i1+1:i2+3)
! !$acc& ,jbin(:),hby(:,:),hbx(:,:),ebz(i1+1:i2+3)
! !$acc& ,eby(i1+1:i2+3),ebx(i1+1:i2+3),bsr(i1:i2,:))
!$acc parallel loop gang async(1) default(present)
!$acc&   reduction(+:ESHORT7,ESHORT8,ESHORT5c)
      do 2 i=i1,i2
         if(mv(i).gt.0)then
            ax1=x(i)
            ay1=y(i)
            az1=z(i)
            agx1=x(i)+gx(ica(i-1),ica(i),seq(i))
            agy1=y(i)+gy(ica(i-1),ica(i),seq(i))
            agz1=z(i)+gz(ica(i-1),ica(i),seq(i))
         else
            ax1=ex(i)
            ay1=ey(i)
            az1=ez(i)
            agx1=egx(i)
            agy1=egy(i)
            agz1=egz(i)
         endif
         if(mv(i+1).gt.0)then
            ax2=x(i+1)
            ay2=y(i+1)
            az2=z(i+1)
            agx2=x(i+1)+gx(ica(i),ica(i+1),seq(i+1))
            agy2=y(i+1)+gy(ica(i),ica(i+1),seq(i+1))
            agz2=z(i+1)+gz(ica(i),ica(i+1),seq(i+1))
         else
            ax2=ex(i+1)
            ay2=ey(i+1)
            az2=ez(i+1)
            agx2=egx(i+1)
            agy2=egy(i+1)
            agz2=egz(i+1)
         endif
         if(mv(i+2).gt.0)then
            ax3=x(i+2)
            ay3=y(i+2)
            az3=z(i+2)
            agx3=x(i+2)+gx(ica(i+1),ica(i+2),seq(i+2))
            agy3=y(i+2)+gy(ica(i+1),ica(i+2),seq(i+2))
            agz3=z(i+2)+gz(ica(i+1),ica(i+2),seq(i+2))
         else
            ax3=ex(i+2)
            ay3=ey(i+2)
            az3=ez(i+2)
            agx3=egx(i+2)
            agy3=egy(i+2)
            agz3=egz(i+2)
         endif
         if(mv(i+3).gt.0)then
            ax4=x(i+3)
            ay4=y(i+3)
            az4=z(i+3)
            agx4=x(i+3)+gx(ica(i+2),ica(i+3),seq(i+3))
            agy4=y(i+3)+gy(ica(i+2),ica(i+3),seq(i+3))
            agz4=z(i+3)+gz(ica(i+2),ica(i+3),seq(i+3))
         else
            ax4=ex(i+3)
            ay4=ey(i+3)
            az4=ez(i+3)
            agx4=egx(i+3)
            agy4=egy(i+3)
            agz4=egz(i+3)
         endif
         if(mv(i+4).gt.0)then
            ax5=x(i+4)
            ay5=y(i+4)
            az5=z(i+4)
            agx5=x(i+4)+gx(ica(i+3),ica(i+4),seq(i+4))
            agy5=y(i+4)+gy(ica(i+3),ica(i+4),seq(i+4))
            agz5=z(i+4)+gz(ica(i+3),ica(i+4),seq(i+4))
         else
            ax5=ex(i+4)
            ay5=ey(i+4)
            az5=ez(i+4)
            agx5=egx(i+4)
            agy5=egy(i+4)
            agz5=egz(i+4)
         endif
ccccccE14:
         ax=ax2-ax1
         ay=ay2-ay1
         az=az2-az1
         bx=ax3-ax2
         by=ay3-ay2
         bz=az3-az2
         cx=ax4-ax3
         cy=ay4-ay3
         cz=az4-az3
         abx=ay*bz-az*by
         aby=az*bx-ax*bz
         abz=ax*by-ay*bx
         hand=abx*cx+aby*cy+abz*cz !a(x)b.c, chirality, >0, right-hand
         ar14=(ax1-ax4)**2+(ay1-ay4)**2+(az1-az4)**2
         r14=nint(ar14)
         if(r14.gt.300) r14=300
         if(hand.lt.0) r14=-r14 !<0, left-hand three-bond.
         ESHORT7=ESHORT7+asr(i,ibin(r14)) !asr(i,dis(i,i+4)) from r14*.dat
cccccccccE15:
         ar15=(ax1-ax5)**2+(ay1-ay5)**2+(az1-az5)**2
         r15=nint(ar15)
         if(r15.gt.500) r15=500
         ESHORT8=ESHORT8+bsr(i,jbin(r15))

cccccccccbias1: encourage helix/sheet-like structure
         if(ar15.lt.75)THEN     !7.53A: helix
            if(sec(i+1).ne.4.AND.sec(i+2).ne.4.AND.sec(i+3).ne.4)then
               if(ibin(r14).gt.4.AND.ibin(r14).lt.8)then !right-hand,3A->8A
                  dot13=(ax4-ax3)*(ax2-ax1)
     $                 +(ay4-ay3)*(ay2-ay1)+(az4-az3)*(az2-az1)
                  if(dot13.lt.0)then
                     dot24=(ax5-ax4)*(ax3-ax2)
     $                    +(ay5-ay4)*(ay3-ay2)+(az5-az4)*(az3-az2)
                     if(dot24.lt.0)then
                        dot14=(ax5-ax4)*(ax2-ax1)
     $                       +(ay5-ay4)*(ay2-ay1)+(az5-az4)*(az2-az1)
                        if(dot14.gt.0)then
                           ff=(afs(i)*afs(i+1)+afs(i+3)*afs(i+4))/2
                           ESHORT5c=ESHORT5c-2-ff
                        endif
                     endif
                  endif
               endif
            endif
         elseif(ar15.gt.160)then !11A, beta
            if(sec(i+1).ne.2.AND.sec(i+2).ne.2.AND.sec(i+3).ne.2)then
               if(mv(i+1).gt.0)then
                  bx2=hbx(ica(i),ica(i+1))
                  by2=hby(ica(i),ica(i+1))
                  bz2=hbz(ica(i),ica(i+1))
               else
                  bx2=ebx(i+1)  !Hb
                  by2=eby(i+1)
                  bz2=ebz(i+1)
               endif
               if(mv(i+3).gt.0)then
                  bx4=hbx(ica(i+2),ica(i+3))
                  by4=hby(ica(i+2),ica(i+3))
                  bz4=hbz(ica(i+2),ica(i+3))
               else
                  bx4=ebx(i+3)
                  by4=eby(i+3)
                  bz4=ebz(i+3)
               endif
               dot24=bx2*bx4+by2*by4+bz2*bz4 !cos(H1.H3)
               if(dot24.gt.0.71)then !angle<45
                  if(mv(i+2).gt.0)then
                     bx3=hbx(ica(i+1),ica(i+2))
                     by3=hby(ica(i+1),ica(i+2))
                     bz3=hbz(ica(i+1),ica(i+2))
                  else
                     bx3=ebx(i+2)
                     by3=eby(i+2)
                     bz3=ebz(i+2)
                  endif
                  dot23=bx2*bx3+by2*by3+bz2*bz3 !cos(H1.H2)
                  if(dot23.lt.-0.71)then !angle >135
                     ff=(afs(i)*afs(i+1)+afs(i+3)*afs(i+4))/2
c                     ESHORT5c=ESHORT5c-2-ff
                     ESHORT5c=ESHORT5c-2-ff*2
                  endif
               endif
            endif
         endif
 2    continue
! !$acc end data
c^^^^^^^^^^E14, E_15, bias1 finished ^^^^^^^^^^^^^^^^^^^^^^^^^^

****  Protein-like bias2: b(i), b(i+4) parallel------->
      i1=max(iiii-4,1)
      i2=min(jjjj-1,Lch-5)
!$acc parallel loop gang  reduction(+:ESHORT5) async(1) default(present)
      do i=i1,i2
         ff=(afs(i)*afs(i+1)+afs(i+3)*afs(i+4))/2 !ff=1, r<r0; =0.25, r>2r0.
         if(mv(i).gt.0)then
            bx1=hbx(ica(i-1),ica(i))
            by1=hby(ica(i-1),ica(i))
            bz1=hbz(ica(i-1),ica(i))
         else
            bx1=ebx(i)
            by1=eby(i)
            bz1=ebz(i)
         endif
         if(mv(i+4).gt.0)then
            bx5=hbx(ica(i+3),ica(i+4))
            by5=hby(ica(i+3),ica(i+4))
            bz5=hbz(ica(i+3),ica(i+4))
         else
            bx5=ebx(i+4)
            by5=eby(i+4)
            bz5=ebz(i+4)
         endif
         b15=bx1*bx5+by1*by5+bz1*bz5
         if(seq(i).eq.2.and.seq(i+2).eq.2.and.seq(i+4).eq.2)then
            if(b15.gt.0.9)then
               ESHORT5=ESHORT5-ff !alpha
            endif
         else
            if(b15.gt.0.5.or.b15.lt.-0.3)then !beta or turn
c               ESHORT5=ESHORT5-ff
               ESHORT5=ESHORT5-ff*2
            endif
         endif
      enddo
c^^^^^^^^^^^^^^^^^^ bias2 finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

****  Protein-like bias3: c(i),c(i+2) anti/parallel:
      i1=max(iiii-2,1)
      i2=min(jjjj,Lch-2)
! !$acc data  copyin(ebx(i1:i2+4),eby(i1:i2+4),ebz(i1:i2+4),hbx(:,:)
! !$acc& ,hbz(:,:),hby(:,:))
!$acc parallel loop gang worker vector reduction(+:ESHORT5) async(1) 
!$acc& default(present) 
      do i=i1,i2
         ff=(afs(i)+afs(i+1)+afs(i+2))/3
         if(mv(i).gt.0)then
            cx1=cax(ica(i-1),ica(i))
            cy1=cay(ica(i-1),ica(i))
            cz1=caz(ica(i-1),ica(i))
         else
            cx1=ecx(i)
            cy1=ecy(i)
            cz1=ecz(i)
         endif
         if(mv(i+2).gt.0)then
            cx3=cax(ica(i+1),ica(i+2))
            cy3=cay(ica(i+1),ica(i+2))
            cz3=caz(ica(i+1),ica(i+2))
         else
            cx3=ecx(i+2)
            cy3=ecy(i+2)
            cz3=ecz(i+2)
         endif
         c13=abs(cx1*cx3+cy1*cy3+cz1*cz3)
         c13=min(0.71,c13)/0.71 !c13 is the same in [0,45]
         ESHORT5=ESHORT5-ff*c13
      enddo

c^^^^^^^^^^^^^^^^^^ bias3 finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ccccccbias4a: to predicted alpha fragment --------------->
      i1=max(iiii-6,1)
      i2=min(jjjj-1,Lch-7)
!$acc parallel loop reduction(+:ESHORT5b) async(1) default(present) gang(1024)
      do i=i1,i2
         if(frga(i).gt.1)then
            dis=sqrt((aax(i)-aax(i+7))**2+
     $           (aay(i)-aay(i+7))**2+(aaz(i)-aaz(i+7))**2)
            ESHORT5b=ESHORT5b+abs(dis-frga(i))
         endif
      enddo
ccccccbias4b: to predicted beta fragment --------------->
      i1=max(iiii-5,1)
      i2=min(jjjj-1,Lch-6)
!$acc parallel loop reduction(+:ESHORT5b) async(1) default(present)
      do i=i1,i2
         if(frgb(i).gt.1)then
            dis=sqrt((aax(i)-aax(i+6))**2+
     $           (aay(i)-aay(i+6))**2+(aaz(i)-aaz(i+6))**2)
c            ESHORT5b=ESHORT5b+abs(dis-frgb(i))
            ESHORT5b=ESHORT5b+abs(dis-frgb(i))*2
         endif
      enddo
c^^^^^^^^^^^^^^^ predicted alpha/beta bias finished ^^^^^^^^^^^^^^^

ccc   penality for crumpling structures ---------------------------->
      i1=max(iiii-11,1)
      i2=min(jjjj-1,Lch-12)
!$acc parallel loop reduction(+:ESHORT5a) async(1) default(present)
      do i=i1,i2
         if(mv(i).gt.0)then
            ax0=x(i)
            ay0=y(i)
            az0=z(i)
         else
            ax0=ex(i)
            ay0=ey(i)
            az0=ez(i)
         endif
         if(mv(i+4).gt.0)then
            ax4=x(i+4)
            ay4=y(i+4)
            az4=z(i+4)
         else
            ax4=ex(i+4)
            ay4=ey(i+4)
            az4=ez(i+4)
         endif
         if(mv(i+8).gt.0)then
            ax8=x(i+8)
            ay8=y(i+8)
            az8=z(i+8)
         else
            ax8=ex(i+8)
            ay8=ey(i+8)
            az8=ez(i+8)
         endif
         avx1=ax4-ax0
         avy1=ay4-ay0
         avz1=az4-az0
         avx2=ax8-ax4
         avy2=ay8-ay4
         avz2=az8-az4
         aaa=avx1*avx2+avy1*avy2+avz1*avz2
         if(aaa.lt.0)then
            if(mv(i+12).gt.0)then
               ax12=x(i+12)
               ay12=y(i+12)
               az12=z(i+12)
            else
               ax12=ex(i+12)
               ay12=ey(i+12)
               az12=ez(i+12)
            endif
            avx3=ax12-ax8
            avy3=ay12-ay8
            avz3=az12-az8
            bbb=avx1*avx3+avy1*avy3+avz1*avz3
            if(bbb.gt.0)then
               ccc=avx2*avx3+avy2*avy3+avz2*avz3
               if(ccc.lt.0)then
                  ESHORT5a=ESHORT5a+1 !crumpling
               endif
            endif
         endif
      enddo
!$acc wait
! !$acc end data
c^^^^^^^^^^^^^ penality of bizard structure finished (ESC1) ^^^^^^^^^^^^

c     Further penalize deriviation of restrain, if larger than (colim, dilim):
      IF(istat.eq.10) then      !calculate E from the beginning
cc         if(dinew.gt.dilim) ESHORT=ESHORT+er1*(dinew-dilim)
cc         if(conew.gt.colim) ESHORT=ESHORT+er3*(conew-colim)
         if(conew.gt.colim) ESHORT4a=conew-colim !colim=number of restrains.
      endif
c     codevsum and didevsum are backuped as total deviation for old 
c     conformation after the first istat=10 are finished

c     the panalty was not counted in Enew, the following is the difference
c     of panlity on new and old conformation, i.e. dE=Enew-Eold. Since
c     this part energy appear only at Eold, the sign is oppsite:
      if(istat.lt.0) then       !return to old energy
         codev=codevsum+conew-coold !total deviation for new conformation
         didev=didevsum+dinew-diold !total deviation for new conformation
c     total penalty-energy beacuse of restrain-deviation on old conformation
cc         if(codevsum.gt.colim) ESHORT=ESHORT+er3*(codevsum-colim)
cc         if(didevsum.gt.dilim) ESHORT=ESHORT+er1*(didevsum-dilim)
         if(codevsum.gt.colim) ESHORT4a=codevsum-colim
c     total penalty-energy beacuse of restrain-deviation on new conformation:
cc         if(codev.gt.colim) ESHORT=ESHORT-er3*(codev-colim)
cc         if(didev.gt.dilim) ESHORT=ESHORT-er1*(didev-dilim)
         if(codev.gt.colim) ESHORT4a=ESHORT4a-(codev-colim)
      endif

***   panelty for the deviation from template---------------->
      if(switch.eq.5)then
         if(ISTAT.eq.10)then    !total energy
            do i=1,nfr
               nn=0
               do j=nfr_i(i),nfr_f(i)
                  nn=nn+1
                  r_1(1,nn)=ex0(j)
                  r_1(2,nn)=ey0(j)
                  r_1(3,nn)=ez0(j)
                  r_2(1,nn)=ex(j)
                  r_2(2,nn)=ey(j)
                  r_2(3,nn)=ez(j)
               enddo
               call u3b(w,r_1,r_2,nn,0,rms,u,t,ier)
               armsd=dsqrt(rms/nn) !RMSD is real, rms is double precision
               ESHORT11=ESHORT11+armsd**4
            enddo
         else                   !fragment energy
            if(i_E_defo.eq.1)then !for defo_M but not trot_M
               if(iiii.eq.1)then
                  i=i_chunk(1)
               elseif(iiii.eq.Lch)then
                  i=i_chunk(Lch)
               else
                  i=i_chunk(iiii+2) !iiii+2 is first residue of the fragment
               endif
               if(i.gt.0)then
                  nn=0
                  do j=nfr_i(i),nfr_f(i)
                     nn=nn+1
                     r_1(1,nn)=ex0(j)
                     r_1(2,nn)=ey0(j)
                     r_1(3,nn)=ez0(j)
                     r_2(1,nn)=ex(j)
                     r_2(2,nn)=ey(j)
                     r_2(3,nn)=ez(j)
                  enddo
                  call u3b(w,r_1,r_2,nn,0,rms,u,t,ier)
                  armsd=dsqrt(rms/nn) !RMSD is real, rms is double precision
                  ESHORT11=ESHORT11+armsd**4
               endif
            endif
         endif
      endif
***^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

********distance restraints------------------------->
      if(nstick.eq.1)then
!$acc parallel loop reduction(+:eshort12) async(1) default(present)
         do i=iiii,jjjj
            if(iq(i,n_tem(itemp)).eq.1)then
               adis=(aax(i)-ax00(i,n_tem(itemp)))**2+
     &              (aay(i)-ay00(i,n_tem(itemp)))**2+
     &              (aaz(i)-az00(i,n_tem(itemp)))**2
               adis=sqrt(adis)
               ESHORT12=ESHORT12+adis 
            endif
         enddo
      endif
!$acc wait
!$acc end data
*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

********distance restraints for RMSD ------------------------->
      if(nrmsd.eq.1)then
         nn=0
         do i=1,Lch
            if(iq(i,n_tem(itemp)).eq.1)then
               nn=nn+1
               r_1(1,nn)=ax00(i,n_tem(itemp))
               r_1(2,nn)=ay00(i,n_tem(itemp))
               r_1(3,nn)=az00(i,n_tem(itemp))
               r_2(1,nn)=aax(i)
               r_2(2,nn)=aay(i)
               r_2(3,nn)=aaz(i)
            endif
         enddo
         call u3b(w,r_1,r_2,nn,0,rms,u,t,ier)
         armsd=dsqrt(rms/nn)    !RMSD is real, rms is double precision
         ESHORT13=armsd**2
         if(armsd.lt.armsd_min0)then
           armsd_min0=armsd
           itemp0=itemp
           icycle0=icycle
           do i=1,Lch
             bx00(i)=aax(i)
             by00(i)=aay(i)
             bz00(i)=aaz(i)
c             write(*,*)i,bx00(i),by00(i),bz00(i)
           enddo
         endif
      else
         ESHORT13=0
      endif
*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      ESHORT=
     $     +es2*ESHORT2
     $     +er1*ESHORT3
     $     +er3*ESHORT4
     $     +er4*ESHORT4a
     $     +es3*ESHORT5
     $     +es3a*ESHORT5a
     $     +es3b*ESHORT5b
     $     +es3c*ESHORT5c
     $     +es4*ESHORT6
     $     +es5*ESHORT7
     $     +es6*ESHORT8
     $     +er5*ESHORT9
     $     +er6*ESHORT10
     $     +er7*ESHORT11
     $     +astick*ESHORT12
     $     +ermsd*ESHORT13
     $     +er21*ESHORT21
     $     +er22*ESHORT22
     $     +er23*ESHORT23

c     $     +er17*ESHORT17
      

c      write(*,*)iiii,jjjj,ESHORT,ESHORT21,er21
c     write(*,*)eshort,eshort12,astick,itemp

c ^^^^^^^^^^ E_short finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! !$acc end data
      RETURN
      END
