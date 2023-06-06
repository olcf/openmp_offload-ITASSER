      subroutine trot_N
      use params
      use backup2
      !use openacc
      use chainm
      use chain1
      use echain1
      use lengths
      use echain4
      use echain2
      use envir1
      use order
      use seqe
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

!      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
!      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

!      dimension ex_o(ndim),ey_o(ndim),ez_o(ndim) !CA
!      dimension egx_o(ndim),egy_o(ndim),egz_o(ndim) !SG
!      dimension ecx_o(ndim),ecy_o(ndim),ecz_o(ndim) !cc
!      dimension ebx_o(ndim),eby_o(ndim),ebz_o(ndim) !Hb
!      dimension etx_o(ndim),ety_o(ndim),etz_o(ndim) !CB
!     dimension ex_n(ndim),ey_n(ndim),ez_n(ndim) !CA
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
      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)
      common/ranzy/nozy
      common/lattice/m_latt,latt1,latt2

      n=nfr_f(ifr)+2            !ending point of movement
      n1=n-1
      n2=n-2                    !ending point of fragment
      n3=n-3
************** rotate the fragment **************************
***   
c     d_xyz0=0.5               !range of translation movement
      ax0=(aranzy(nozy)*2-1)*d_xyz00(ifr) !translation of coordinates, [-det,det]
      ay0=(aranzy(nozy)*2-1)*d_xyz00(ifr)
      az0=(aranzy(nozy)*2-1)*d_xyz00(ifr)
c     rotation axis (awx,awy,awz):
      acos_theta=1-2.*aranzy(nozy)
      asin_theta=sqrt(1.0-acos_theta*acos_theta)
      aphi=2.*3.1415926*aranzy(nozy) ![0,2*pai]
      awx=asin_theta*cos(aphi)
      awy=asin_theta*sin(aphi)
      awz=acos_theta
c     rotation angle           !in arc-angle, 1=57o, rotation angle.
      ang=(aranzy(nozy)*2-1)*angle00(ifr) ![-angle00,angle00]
c     rotation points:
      ax=ex(n2)
      ay=ey(n2)
      az=ez(n2)
***   rotation matrix:
      acos=cos(ang)
      asin=sin(ang)
      a11=awx*awx*(1-acos)+acos
      a12=awx*awy*(1-acos)-awz*asin
      a13=awx*awz*(1-acos)+awy*asin
      a21=awx*awy*(1-acos)+awz*asin
      a22=awy*awy*(1-acos)+acos
      a23=awy*awz*(1-acos)-awx*asin
      a31=awx*awz*(1-acos)-awy*asin
      a32=awy*awz*(1-acos)+awx*asin
      a33=awz*awz*(1-acos)+acos
***   check n-terminal of fragment--------------->
      ax1=ax0+ax+(ex(n2)-ax)*a11+(ey(n2)-ay)*a12+(ez(n2)-az)*a13 !v1=v0+(a)v
      ay1=ay0+ay+(ex(n2)-ax)*a21+(ey(n2)-ay)*a22+(ez(n2)-az)*a23
      az1=az0+az+(ex(n2)-ax)*a31+(ey(n2)-ay)*a32+(ez(n2)-az)*a33
***   draw distant fragment closer ----->
      r2_bond=(x(n1)-nint(ex(n2)))**2+(y(n1)-nint(ey(n2)))**2
     $     +(z(n1)-nint(ez(n2)))**2
      nbig=0
      if(r2_bond.gt.latt2)then
         dis2_old=(x(n1)-ex(n2))**2+(y(n1)-ey(n2))**2+(z(n1)-ez(n2))**2
 33      t4=int(aranzy(nozy)*anvec)+1 !random vector
         if(.not.goodc(t4,ica(n)))goto 33
         xn1=x(n)-vx(t4)
         yn1=y(n)-vy(t4)
         zn1=z(n)-vz(t4)
         dis2_new=(xn1-ax1)**2+(yn1-ay1)**2+(zn1-az1)**2
         if(dis2_new.lt.dis2_old)then
 34         t3=int(aranzy(nozy)*anvec)+1 !random vector
            if(.not.goodc(t3,t4))goto 34
            nbig=1
            goto 204            !take t1, t2
         endif
      endif
***
      x1=nint(ax1)
      y1=nint(ay1)
      z1=nint(az1)
      ijx=x(n)-x1
      if(abs(ijx).gt.10)goto 202 !without correct movement, quit
      ijy=y(n)-y1
      if(abs(ijy).gt.10)goto 202 !without correct movement, quit
      ijz=z(n)-z1
      if(abs(ijz).gt.10)goto 202 !without correct movement, quit
      nc=Nw(ijx,ijy,ijz)
      if(nc.lt.1)goto 202       !without correct movement, quit
      p=int(aranzy(nozy)*nc)+1    ![1,nc]
      t3=W21(ijx,ijy,ijz,p)
      t4=W22(ijx,ijy,ijz,p)
      if(.not.goodc(t4,ica(n)))goto 202

c     backup old path ------------>
 204  oo(n2)=ica(n2)
      oo(n1)=ica(n1)
      ox(n1)=x(n1)
      oy(n1)=y(n1)
      oz(n1)=z(n1)
!$acc kernels
      do i=1,n2
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
!$acc end kernels
c     prepare the new path------------>
      nn(n2)=t3
      nn(n1)=t4
      nx(n1)=x(n)-vx(nn(n1))
      ny(n1)=y(n)-vy(nn(n1))
      nz(n1)=z(n)-vz(nn(n1))
!$acc kernels
      do i=1,n2
        ex_n(i)=ax0+ax+(ex(i)-ax)*a11+(ey(i)-ay)*a12+(ez(i)-az)*a13 !CA
        ey_n(i)=ay0+ay+(ex(i)-ax)*a21+(ey(i)-ay)*a22+(ez(i)-az)*a23
        ez_n(i)=az0+az+(ex(i)-ax)*a31+(ey(i)-ay)*a32+(ez(i)-az)*a33
        egx_n(i)=ax0+ax+(egx(i)-ax)*a11+(egy(i)-ay)*a12+(egz(i)-az)*a13 !SG
        egy_n(i)=ay0+ay+(egx(i)-ax)*a21+(egy(i)-ay)*a22+(egz(i)-az)*a23
        egz_n(i)=az0+az+(egx(i)-ax)*a31+(egy(i)-ay)*a32+(egz(i)-az)*a33
        ecx_n(i)=ecx(i)*a11+ecy(i)*a12+ecz(i)*a13 !cc
        ecy_n(i)=ecx(i)*a21+ecy(i)*a22+ecz(i)*a23
        ecz_n(i)=ecx(i)*a31+ecy(i)*a32+ecz(i)*a33
        ebx_n(i)=ebx(i)*a11+eby(i)*a12+ebz(i)*a13 !bb
        eby_n(i)=ebx(i)*a21+eby(i)*a22+ebz(i)*a23
        ebz_n(i)=ebx(i)*a31+eby(i)*a32+ebz(i)*a33
        etx_n(i)=ax0+ax+(etx(i)-ax)*a11+(ety(i)-ay)*a12+(etz(i)-az)*a13 !CB
        ety_n(i)=ay0+ay+(etx(i)-ax)*a21+(ety(i)-ay)*a22+(etz(i)-az)*a23
        etz_n(i)=az0+az+(etx(i)-ax)*a31+(ety(i)-ay)*a32+(etz(i)-az)*a33
      enddo
!$acc end kernels
      d2=(nx(n1)-ex_n(n3))**2+(ny(n1)-ey_n(n3))**2+(nz(n1)-ez_n(n3))**2
      if((d2.lt.23.or.d2.gt.78).and.nbig.eq.0)goto 202

c     change conformation to new path -------------->
      ica(n2)=nn(n2)
      ica(n1)=nn(n1)
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)
      x(n1)=nx(n1)
      y(n1)=ny(n1)
      z(n1)=nz(n1)
!$acc kernels
      do i=1,n2
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
      if(look(1,n))then         ! check excluded volumn for passage of [m,n]
c     calculate E_new--------------->
!$acc kernels
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
! !$acc update device(ex,ey,ez,egx,egy,egz,ecx,ecy,ecz,
! !$acc&   ebx,eby,ebz,etx,ety,etz,x,y,z)
         Enew=EHB(1,n,1)+ESHORT(1,n,1) !use ifs
         do kkk=1,n
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo
!$acc end kernels
c     return back the conformation and calculate E_old --------->
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
!$acc kernels
         do i=1,n2
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
            etx(i)=ebx_o(i)     !CB
            ety(i)=eby_o(i)
            etz(i)=ebz_o(i)
         enddo
! !$acc update device(ex,ey,ez,egx,egy,egz,ecx,ecy,ecz,
! !$acc&   ebx,eby,ebz,etx,ety,etz,x,y,z)
         Eold=EHB(1,n,-1)+ESHORT(1,n,-1)

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
            bNa(22)=bNa(22)+1
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
            do kkk=1,n
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
!$acc end kernels
            codevsum=codev
            didevsum=didev	
c     energ=energ+de	
            ica(n2)=nn(n2)
            ica(n1)=nn(n1)
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            x(n1)=nx(n1)
            y(n1)=ny(n1)
            z(n1)=nz(n1)
!$acc kernels
            do i=1,n2
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
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
!$acc kernels
         do i=1,n2
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
      bNt(22)=bNt(22)+1
      bNNt(itemp)=bNNt(itemp)+1

! !$acc update device(ex,ey,ez,egx,egy,egz,ecx,ecy,ecz,noaa,
! !$acc&       ebx,eby,ebz,etx,ety,etz,x,y,z,nopp,nomm)

c ^^^^^^^^^^^^^^^^^ trot_N finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
