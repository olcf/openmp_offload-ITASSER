      subroutine defo_M
	use params
      use backup2
      use chainm
      !use openacc
      use chain1
      use echain1
      use lengths
      use echain2
      use echain4
      use order
      use envir1
      use seqe
      use echain5
      use echain6
      use short
      use shape
      use backup1
      use ehbc
      use tempArrays
      implicit integer(i-z)
!      parameter(ndim=1500)
!      parameter(nrep=100)       !number of replicas
!      parameter(nvec=416)
      parameter(apai=3.1415926)
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
      common/defoangle/defo_angle
      common/rmsd_defo/armsd,armsd_sum(100),N_rmsd(100)
      common/ranzy/nozy
      common/lattice/m_latt,latt1,latt2
!$acc wait
      m=nfr_i(ifr)-2            !starting point of movement
      n=nfr_f(ifr)+2            !ending point of movement
      m1=m+1
      m2=m+2                    !starting point of fragment
      m3=m+3
      n1=n-1
      n2=n-2                    !ending point of fragment
      n3=n-3
************** rotate the fragment **************************
***   
ccc   d_xyz0=0.5                  !maximum translation coordinates
      ax0=(aranzy(nozy)*2-1)*d_xyz00(ifr) !translation of coordinates, [-det,det]
      ay0=(aranzy(nozy)*2-1)*d_xyz00(ifr)
      az0=(aranzy(nozy)*2-1)*d_xyz00(ifr)
ccc   rotation axis (awx,awy,awz)->N-terminal; (bwx,bwy,bwz)->C-terminal:
      acos_theta=1-2.*aranzy(nozy)
      asin_theta=sqrt(1.0-acos_theta*acos_theta)
      aphi=2.*3.1415926*aranzy(nozy) ![0,2*pai]
      awx=asin_theta*cos(aphi)
      awy=asin_theta*sin(aphi)
      awz=acos_theta
      acos_theta=1-2.*aranzy(nozy)
      asin_theta=sqrt(1.0-acos_theta*acos_theta)
      aphi=2.*3.1415926*aranzy(nozy) ![0,2*pai]
      bwx=asin_theta*cos(aphi)
      bwy=asin_theta*sin(aphi)
      bwz=acos_theta
ccc   rotation angle           !in arc-angle, 1=57o, rotation angle.
      ang=(aranzy(nozy)*2-1)*angle00(ifr)/defo_angle ![-angle00/d,angle00/d]
      bng=(aranzy(nozy)*2-1)*angle00(ifr)/defo_angle ![-angle00/d,angle00/d]
ccc   rotation points:
      n_rot=m2+int((n2-m2+1)*aranzy(nozy)) ![m2,n2]
      ax=ex(n_rot)
      ay=ey(n_rot)
      az=ez(n_rot)
ccc   rotation matrix:
      acos=cos(ang)             !for N-terminal of the fragments--->
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
      acos=cos(bng)             !for N-terminal of the fragments--->
      asin=sin(bng)
      b11=bwx*bwx*(1-acos)+acos
      b12=bwx*bwy*(1-acos)-bwz*asin
      b13=bwx*bwz*(1-acos)+bwy*asin
      b21=bwx*bwy*(1-acos)+bwz*asin
      b22=bwy*bwy*(1-acos)+acos
      b23=bwy*bwz*(1-acos)-bwx*asin
      b31=bwx*bwz*(1-acos)-bwy*asin
      b32=bwy*bwz*(1-acos)+bwx*asin
      b33=bwz*bwz*(1-acos)+acos
***   check m-terminal of fragment--------------->
      ax1=ax0+ax+(ex(m2)-ax)*a11+(ey(m2)-ay)*a12+(ez(m2)-az)*a13 !v1=v0+(a)v
      ay1=ay0+ay+(ex(m2)-ax)*a21+(ey(m2)-ay)*a22+(ez(m2)-az)*a23
      az1=az0+az+(ex(m2)-ax)*a31+(ey(m2)-ay)*a32+(ez(m2)-az)*a33
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
***   check n-terminal of fragment--------------->
 203  ax1=ax0+ax+(ex(n2)-ax)*b11+(ey(n2)-ay)*b12+(ez(n2)-az)*b13
      ay1=ay0+ay+(ex(n2)-ax)*b21+(ey(n2)-ay)*b22+(ez(n2)-az)*b23
      az1=az0+az+(ex(n2)-ax)*b31+(ey(n2)-ay)*b32+(ez(n2)-az)*b33
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
 204  oo(m)=ica(m)
      oo(m1)=ica(m1)
      ox(m1)=x(m1)
      oy(m1)=y(m1)
      oz(m1)=z(m1)
      oo(n2)=ica(n2)
      oo(n1)=ica(n1)
      ox(n1)=x(n1)
      oy(n1)=y(n1)
      oz(n1)=z(n1)
!$acc kernels async(1)
      do i=m2,n2
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
      nn(m)=t1
      nn(m1)=t2
      nx(m1)=x(m)+vx(nn(m))
      ny(m1)=y(m)+vy(nn(m))
      nz(m1)=z(m)+vz(nn(m))
      nn(n2)=t3
      nn(n1)=t4
      nx(n1)=x(n)-vx(nn(n1))
      ny(n1)=y(n)-vy(nn(n1))
      nz(n1)=z(n)-vz(nn(n1))
!$acc kernels  async(1)
      do i=m2,n_rot
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
      do i=n_rot+1,n2
         ex_n(i)=ax0+ax+(ex(i)-ax)*b11+(ey(i)-ay)*b12+(ez(i)-az)*b13 !CA
         ey_n(i)=ay0+ay+(ex(i)-ax)*b21+(ey(i)-ay)*b22+(ez(i)-az)*b23
         ez_n(i)=az0+az+(ex(i)-ax)*b31+(ey(i)-ay)*b32+(ez(i)-az)*b33
         egx_n(i)=ax0+ax+(egx(i)-ax)*b11+(egy(i)-ay)*b12+(egz(i)-az)*b13 !SG
         egy_n(i)=ay0+ay+(egx(i)-ax)*b21+(egy(i)-ay)*b22+(egz(i)-az)*b23
         egz_n(i)=az0+az+(egx(i)-ax)*b31+(egy(i)-ay)*b32+(egz(i)-az)*b33
         ecx_n(i)=ecx(i)*b11+ecy(i)*b12+ecz(i)*b13 !cc
         ecy_n(i)=ecx(i)*b21+ecy(i)*b22+ecz(i)*b23
         ecz_n(i)=ecx(i)*b31+ecy(i)*b32+ecz(i)*b33
         ebx_n(i)=ebx(i)*b11+eby(i)*b12+ebz(i)*b13 !bb
         eby_n(i)=ebx(i)*b21+eby(i)*b22+ebz(i)*b23
         ebz_n(i)=ebx(i)*b31+eby(i)*b32+ebz(i)*b33
         etx_n(i)=ax0+ax+(etx(i)-ax)*b11+(ety(i)-ay)*b12+(etz(i)-az)*b13 !CB
         ety_n(i)=ay0+ay+(etx(i)-ax)*b21+(ety(i)-ay)*b22+(etz(i)-az)*b23
         etz_n(i)=az0+az+(etx(i)-ax)*b31+(ety(i)-ay)*b32+(etz(i)-az)*b33
      enddo
!$acc end kernels
      d2=(nx(m1)-ex_n(m3))**2+(ny(m1)-ey_n(m3))**2+(nz(m1)-ez_n(m3))**2
      if((d2.lt.23.or.d2.gt.78).and.mbig.eq.0)goto 202
      d2=(nx(n1)-ex_n(n3))**2+(ny(n1)-ey_n(n3))**2+(nz(n1)-ez_n(n3))**2
      if((d2.lt.23.or.d2.gt.78).and.nbig.eq.0)goto 202

c     change conformation to new path -------------->
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      ica(n2)=nn(n2)
      ica(n1)=nn(n1)
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)
      x(n1)=nx(n1)
      y(n1)=ny(n1)
      z(n1)=nz(n1)
!$acc kernels  async(1)
      do i=m2,n2
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
      if(look(m,n))then         ! check excluded volumn for passage of [m,n]
c     calculate E_new--------------->
!$acc kernels  async(1)
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
!$acc end kernels ! Just added this
         Enew=EHB(m,n,1)+ESHORT(m,n,1) !use ifs

!$acc kernels async(1) ! Just added this
         do kkk=m,n
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo
!$acc end kernels
c     return back the conformation and calculate E_old --------->
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
!$acc kernels
         do i=m2,n2
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
!$acc end kernels ! Just added this
         Eold=EHB(m,n,-1)+ESHORT(m,n,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
!$acc kernels  async(1)
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
            bNa(26)=bNa(26)+1
            bNNa(itemp)=bNNa(itemp)+1
            armsd_sum(itemp)=armsd_sum(itemp)+armsd
            N_rmsd(itemp)=N_rmsd(itemp)+1
            icnto=icnt
            sumcto=sumct 
!$acc kernels async(1)
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m,n
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
!$acc end kernels
            codevsum=codev
            didevsum=didev	
c     energ=energ+de	
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            ica(n2)=nn(n2)
            ica(n1)=nn(n1)
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            x(n1)=nx(n1)
            y(n1)=ny(n1)
            z(n1)=nz(n1)
!$acc kernels  async(1)
            do i=m2,n2
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
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
!$acc kernels  async(1)
         do i=m2,n2
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
!$acc wait(2)
 202  continue
      bNt(26)=bNt(26)+1
      bNNt(itemp)=bNNt(itemp)+1
! !$acc update device(egx,egy,egz)

c ^^^^^^^^^^^^^^^^^ defo_M finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
