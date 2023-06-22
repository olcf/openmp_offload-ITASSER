      subroutine tran_C
      use params
      !use openacc
      use backup2
      use chainm
      use chain1
      use echain1
      use lengths
      use echain2
      use echain4
      use seqe
      use envir1
      use order
      use echain5
      use echain6
      use short
      use shape
      use backup1
      use ehbc
      use tempArrays
      implicit integer(i-z)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
!      common/lengths/Lch,Lch1,Lch2
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
!      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
!      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
!      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
!      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
!      common/backup2/eprofo,eprofn,energ
!      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
!      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
!      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)

!      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
!      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

!      dimension ex_o(ndim),ey_o(ndim),ez_o(ndim) !CA
!      dimension egx_o(ndim),egy_o(ndim),egz_o(ndim) !SG
!      dimension ecx_o(ndim),ecy_o(ndim),ecz_o(ndim) !cc
!      dimension ebx_o(ndim),eby_o(ndim),ebz_o(ndim) !Hb
!      dimension etx_o(ndim),ety_o(ndim),etz_o(ndim) !CB
!      dimension ex_n(ndim),ey_n(ndim),ez_n(ndim) !CA
!      dimension egx_n(ndim),egy_n(ndim),egz_n(ndim) !SG
!      dimension ecx_n(ndim),ecy_n(ndim),ecz_n(ndim) !cc
!      dimension ebx_n(ndim),eby_n(ndim),ebz_n(ndim) !Hb
!      dimension etx_n(ndim),ety_n(ndim),etz_n(ndim) !CB

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/rotat/ar0(ndim),athita0(ndim),aphi0(ndim)

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)
      common/w2a/w21(-10:10,-10:10,-10:10,80)
      common/w2b/w22(-10:10,-10:10,-10:10,80)
      common/w2c/Nw(-10:10,-10:10,-10:10)

!      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
!      common/echain1/ex(ndim),ey(ndim),ez(ndim) !CA
!      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
!      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
!      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
!      common/echain6/etx(ndim),ety(ndim),etz(ndim) !CB
!      common/chainm/mv(ndim)
      common/ranzy/nozy
      common/lattice/m_latt,latt1,latt2

      m=nfr_i(ifr)-2            !starting point of movement
      m1=m+1
      m2=m+2                    !starting point of fragment
      m3=m+3
************** rotate the fragment **************************
***   
c     d_xyz0=0.5               !range of translation movement
      ax0=(aranzy(nozy)*2-1)*d_xyz00(ifr) !translation of coordinates, [-det,det]
      ay0=(aranzy(nozy)*2-1)*d_xyz00(ifr)
      az0=(aranzy(nozy)*2-1)*d_xyz00(ifr)
***   check m-terminal of fragment--------------->
      ax1=ax0+ex(m2)
      ay1=ay0+ey(m2)
      az1=az0+ez(m2)
***   draw distant fragment closer ----->
      r2_bond=(x(m1)-nint(ex(m2)))**2+(y(m1)-nint(ey(m2)))**2
     $     +(z(m1)-nint(ez(m2)))**2
      mbig=0
      if(r2_bond.gt.latt2)then
         dis2_old=(x(m1)-ex(m2))**2+(y(m1)-ey(m2))**2+(z(m1)-ez(m2))**2
 31      t1=int(aranzy(nozy)*anvec)+1 !random vector
         if(.not.goodc(ica(m-1),t1))goto 31
         xm1=x(m)+vx(t1)
         ym1=y(m)+vy(t1)
         zm1=z(m)+vz(t1)
         dis2_new=(xm1-ax1)**2+(ym1-ay1)**2+(zm1-az1)**2
         if(dis2_new.lt.dis2_old)then
 32         t2=int(aranzy(nozy)*anvec)+1 !random vector
            if(.not.goodc(t1,t2))goto 32
            mbig=1
            goto 203            !take t1, t2
         endif
      endif
***
      x1=nint(ax1)
      y1=nint(ay1)
      z1=nint(az1)
      ijx=x1-x(m)
      if(abs(ijx).gt.10)goto 202 !without correct movement, quit
      ijy=y1-y(m)
      if(abs(ijy).gt.10)goto 202 !without correct movement, quit
      ijz=z1-z(m)
      if(abs(ijz).gt.10)goto 202 !without correct movement, quit
      nc=Nw(ijx,ijy,ijz)
      if(nc.lt.1)goto 202       !without correct movement, quit
      p=int(aranzy(nozy)*nc)+1    ![1,nc]
      t1=W21(ijx,ijy,ijz,p)
      t2=W22(ijx,ijy,ijz,p)
      if(.not.goodc(ica(m-1),t1))goto 202

c     backup old path ------------>
 203  oo(m)=ica(m)
      oo(m1)=ica(m1)
      ox(m1)=x(m1)
      oy(m1)=y(m1)
      oz(m1)=z(m1)
!$acc kernels 
!$OMP target teams distribute parallel do simd !nowait
!$OMP& map(to:ex,ey,ez,egx,egy,egz,ecx,ecy,ecz,ebx,eby,ebz,etx,ety,etz) 
!$OMP& map(from:ex_o,ey_o,ez_o,egx_o,egy_o,egz_o,ecx_o,ecy_o,
!$OMP& ecz_o,ebx_o,eby_o,ebz_o,etx_o,ety_o,etz_o)
      do i=m2,Lch
         ex_o(i)=ex(i)          !CA
         ey_o(i)=ey(i)
         ez_o(i)=ez(i)
         egx_o(i)=egx(i)        !SG
         egy_o(i)=egy(i)
         egz_o(i)=egz(i)
         ecx_o(i)=ecx(i)        !cc
         ecy_o(i)=ecy(i)
         ecz_o(i)=ecz(i)
         ebx_o(i)=ebx(i)        !Hb
         eby_o(i)=eby(i)
         ebz_o(i)=ebz(i)
         etx_o(i)=etx(i)        !CB
         ety_o(i)=ety(i)
         etz_o(i)=etz(i)
      enddo
!$OMP end target teams distribute parallel do simd 
!$acc end kernels
c     prepare the new path------------>
      nn(m)=t1
      nn(m1)=t2
      nx(m1)=x(m)+vx(nn(m))
      ny(m1)=y(m)+vy(nn(m))
      nz(m1)=z(m)+vz(nn(m))
!$acc kernels
      do i=m2,Lch
         ex_n(i)=ax0+ex(i)      !CA
         ey_n(i)=ay0+ey(i)
         ez_n(i)=az0+ez(i)
         egx_n(i)=ax0+egx(i)    !SG
         egy_n(i)=ay0+egy(i)
         egz_n(i)=az0+egz(i)
         ecx_n(i)=ecx(i)        !cc
         ecy_n(i)=ecy(i)
         ecz_n(i)=ecz(i)
         ebx_n(i)=ebx(i)        !Hb
         eby_n(i)=eby(i)
         ebz_n(i)=ebz(i)
         etx_n(i)=ax0+etx(i)    !CB
         ety_n(i)=ay0+ety(i)
         etz_n(i)=az0+etz(i)
      enddo
!$acc end kernels
      d2=(nx(m1)-ex_n(m3))**2+(ny(m1)-ey_n(m3))**2+(nz(m1)-ez_n(m3))**2
      if((d2.lt.23.or.d2.gt.78).and.mbig.eq.0)goto 202

c     change conformation to new path -------------->
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
!$acc kernels
      do i=m2,Lch
         ex(i)=ex_n(i)          !CA
         ey(i)=ey_n(i)
         ez(i)=ez_n(i)
         egx(i)=egx_n(i)        !SG
         egy(i)=egy_n(i)
         egz(i)=egz_n(i)
         ecx(i)=ecx_n(i)        !cc
         ecy(i)=ecy_n(i)
         ecz(i)=ecz_n(i)
         ebx(i)=ebx_n(i)        !Hb
         eby(i)=eby_n(i)
         ebz(i)=ebz_n(i)
         etx(i)=etx_n(i)        !CB
         ety(i)=ety_n(i)
         etz(i)=etz_n(i)
      enddo
!$acc end kernels
      if(look(m,Lch))then         ! check excluded volumn for passage of [m,n]
c     calculate E_new--------------->
!$acc kernels
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m,Lch,1)+ESHORT(m,Lch,1) !use ifs
         do kkk=m,Lch
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo
!$acc end kernels
c     return back the conformation and calculate E_old --------->
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
!$acc kernels
         do i=m2,Lch
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
            etx(i)=etx_o(i)     !CB
            ety(i)=ety_o(i)
            etz(i)=etz_o(i)
         enddo
         Eold=EHB(m,Lch,-1)+ESHORT(m,Lch,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo
!$acc end kernels
c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c     E=energ
c     Et=E+de

         call metro(dE,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(18)=bNa(18)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
!$acc kernels
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m,Lch
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
!$acc end kernels
            codevsum=codev
            didevsum=didev	
c     energ=energ+de	
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
!$acc kernels
            do i=m2,Lch
               ex(i)=ex_n(i)    !CA
               ey(i)=ey_n(i)
               ez(i)=ez_n(i)
               egx(i)=egx_n(i)  !SG
               egy(i)=egy_n(i)
               egz(i)=egz_n(i)
               ecx(i)=ecx_n(i)  !cc
               ecy(i)=ecy_n(i)
               ecz(i)=ecz_n(i)
               ebx(i)=ebx_n(i)  !Hb
               eby(i)=eby_n(i)
               ebz(i)=ebz_n(i)
               etx(i)=etx_n(i)  !CB
               ety(i)=ety_n(i)
               etz(i)=etz_n(i)
            enddo
!$acc end kernels
            call get_center     !update (amx,amy,amz)
c            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
!$acc kernels
         do i=m2,Lch
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
            etx(i)=etx_o(i)     !CB
            ety(i)=ety_o(i)
            etz(i)=etz_o(i)
         enddo
!$acc end kernels
      endif                     !for look(i,m)

 202  continue
      bNt(18)=bNt(18)+1
      bNNt(itemp)=bNNt(itemp)+1

! !$acc update device(ex,ey,ez,egx,egy,egz,ecx,ecy,ecz,afs
! !$acc&    ebx,eby,ebz,etx,ety,etz,x,y,z,afsn)

c ^^^^^^^^^^^^^^^^^ tran_C finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
