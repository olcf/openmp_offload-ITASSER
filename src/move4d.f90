      subroutine move4d
	use params
      use backup2
      !use openacc
      use chainm
      use chain1
      use echain1
      use lengths
      use envir1
      use seqe
      use order
      use short
      use ehbc
      use shape
      use backup1
      implicit integer(i-z)
!      parameter(nrep=100)       !number of replicas
!      parameter(ndim=1500)
!      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
!      common/lengths/Lch,Lch1,Lch2
!      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
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

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/ranzy/nozy

c     choose the position of (m) ----------->
      i=int(aranzy(nozy)*nfl4)+1
      m=ras4(i)
      m1=m+1
      m2=m+2
      m3=m+3
      m4=m+4

ccccccccccccccc 1th temperor 2-bond movement cccccccccccccccc
      nc=Np2(ica(m),ica(m1))
 201  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      t=v21(ica(m),ica(m1),p)
      if(.not.goodc(ica(m-1),t))goto 201
      t1=v22(ica(m),ica(m1),p)
      if(.not.goodc(t1,ica(m2)))goto 201
ccccccccccccccc 2th temperor 2-bond movement cccccccccccccccc
      nc=Np2(t1,ica(m2))
 202  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      tt1=v21(t1,ica(m2),p)
      if(.not.goodc(t,tt1))goto 202
      t2=v22(t1,ica(m2),p)
      if(.not.goodc(t2,ica(m3)))goto 202
ccccccccccccccc 3th temperor 2-bond movement cccccccccccccccc
      nc=Np2(t2,ica(m3))
 203  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      tt2=v21(t2,ica(m3),p)
      if(.not.goodc(tt1,tt2))goto 203
      t3=v22(t2,ica(m3),p)
      if(.not.goodc(t3,ica(m4)))goto 203

ccccccccccccc first 2-bond movement cccccccccccccccccccccccccc
      nc=Np2(t,tt1)
 204  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      nn(m)=v21(t,tt1,p)
      if(.not.goodc(ica(m-1),nn(m)))goto 204
      ttt1=v22(t,tt1,p)
      if(.not.goodc(ttt1,tt2))goto 204
ccccccccccccc second 2-bond movement cccccccccccccccccccccccccc
      nc=Np2(ttt1,tt2)
 205  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      nn(m1)=v21(ttt1,tt2,p)
      if(.not.goodc(nn(m),nn(m1)))goto 205
      ttt2=v22(ttt1,tt2,p)
      if(.not.goodc(ttt2,t3))goto 205
ccccccccccccc third 2-bond movement cccccccccccccccccccccccccc
      nc=Np2(ttt2,t3)
 206  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      nn(m2)=v21(ttt2,t3,p)
      if(.not.goodc(nn(m1),nn(m2)))goto 206
      nn(m3)=v22(ttt2,t3,p)
      if(.not.goodc(nn(m3),ica(m4)))goto 206
c^^^^^^^^ four-bonds movement conformation finished ^^^^^^^^^^
      if(nn(m).eq.ica(m).and.nn(m1).eq.ica(m1).
     $     and.nn(m2).eq.ica(m2))goto 113

c     prepare new path and backup old path ------------>
      ox(m1)=x(m1)              !memory of old path
      oy(m1)=y(m1)              !memory of old path
      oz(m1)=z(m1)              !memory of old path
      ox(m2)=x(m2)              !memory of old path
      oy(m2)=y(m2)              !memory of old path
      oz(m2)=z(m2)              !memory of old path
      ox(m3)=x(m3)              !memory of old path
      oy(m3)=y(m3)              !memory of old path
      oz(m3)=z(m3)              !memory of old path
      oo(m)=ica(m)              !memory of old path
      oo(m1)=ica(m1)            !memory of old path
      oo(m2)=ica(m2)            !memory of old path
      oo(m3)=ica(m3)            !memory of old path
      nx(m1)=x(m)+vx(nn(m))     !memory of new path
      ny(m1)=y(m)+vy(nn(m))     !memory of new path
      nz(m1)=z(m)+vz(nn(m))     !memory of new path
      nx(m2)=nx(m1)+vx(nn(m1))  !memory of new path
      ny(m2)=ny(m1)+vy(nn(m1))  !memory of new path
      nz(m2)=nz(m1)+vz(nn(m1))  !memory of new path
      nx(m3)=nx(m2)+vx(nn(m2))  !memory of new path
      ny(m3)=ny(m2)+vy(nn(m2))  !memory of new path
      nz(m3)=nz(m2)+vz(nn(m2))  !memory of new path

c     change conformation to new path -------------->
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      x(m2)=nx(m2)
      y(m2)=ny(m2)
      z(m2)=nz(m2)
      x(m3)=nx(m3)
      y(m3)=ny(m3)
      z(m3)=nz(m3)
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      ica(m2)=nn(m2)
      ica(m3)=nn(m3)
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)

      if(look(m,m4))then       ! check excluded volumn for passage of [i,m]
c     calculate E_new--------------->
!$acc kernels
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m,m4,1)+ESHORT(m,m4,1) !use ifs
         do kkk=m,m4
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo
!$acc end kernels
c     return back the conformation and calculate E_old --------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         x(m2)=ox(m2)
         y(m2)=oy(m2)
         z(m2)=oz(m2)
         x(m3)=ox(m3)
         y(m3)=oy(m3)
         z(m3)=oz(m3)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(m2)=oo(m2)
         ica(m3)=oo(m3)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         Eold=EHB(m,m4,-1)+ESHORT(m,m4,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
!$acc parallel loop
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c         E=energ
c         Et=E+de

         call metro(dE,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(6)=bNa(6)+1
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
            do kkk=m,m4
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
!$acc end kernels
            codevsum=codev
            didevsum=didev	
c            energ=energ+de	
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            x(m2)=nx(m2)
            y(m2)=ny(m2)
            z(m2)=nz(m2)
            x(m3)=nx(m3)
            y(m3)=ny(m3)
            z(m3)=nz(m3)
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            ica(m2)=nn(m2)
            ica(m3)=nn(m3)
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            call get_center       !update (amx,amy,amz)
c            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         x(m2)=ox(m2)
         y(m2)=oy(m2)
         z(m2)=oz(m2)
         x(m3)=ox(m3)
         y(m3)=oy(m3)
         z(m3)=oz(m3)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(m2)=oo(m2)
         ica(m3)=oo(m3)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
      endif                     !for look(i,m)
113   continue
      bNt(6)=bNt(6)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move4d finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
