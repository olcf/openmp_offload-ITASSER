      subroutine move3d
	use params
      use backup2
!      !use openacc
      use chainm
      use chain1
      use echain1
      use lengths
      use seqe
      use order
      use envir1
      use short
      use shape
      use ehbc
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
      i=int(aranzy(nozy)*nfl3)+1
      m=ras3(i)
      m1=m+1
      m2=m+2
      m3=m+3

ccccccccccccc temporal 1th 2-bond movement cccccccccccccccccc
      nc=Np2(ica(m),ica(m1))
 201  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      t=v21(ica(m),ica(m1),p)
      if(.not.goodc(ica(m-1),t))goto 201
      t1=v22(ica(m),ica(m1),p)
      if(.not.goodc(t1,ica(m2)))goto 201
ccccccccccccc temporal 2th 2-bond movement cccccccccccccccccc
      nc=Np2(t1,ica(m2))
 202  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      tt1=v21(t1,ica(m2),p)
      if(.not.goodc(t,tt1))goto 202
      t2=v22(t1,ica(m2),p)
      if(.not.goodc(t2,ica(m3)))goto 202

ccccccccccccc 1th 2-bond movement cccccccccccccccccccccccccc
      nc=Np2(t,tt1)
 203  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      nn(m)=v21(t,tt1,p)
      if(.not.goodc(ica(m-1),nn(m)))goto 203
      ttt1=v22(t,tt1,p)
      if(.not.goodc(ttt1,t2))goto 203
ccccccccccccc 2th 2-bond movement cccccccccccccccccccccccccc
      nc=Np2(ttt1,t2)
 204  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      nn(m1)=v21(ttt1,t2,p)
      if(.not.goodc(nn(m),nn(m1)))goto 204
      nn(m2)=v22(ttt1,t2,p)
      if(.not.goodc(nn(m2),ica(m3)))goto 204
c^^^^^^^^ three-bonds movement conformation finished ^^^^^^^^^^
      if(nn(m).eq.ica(m).and.nn(m1).eq.ica(m1))goto 113


c     prepare new path and backup old path ------------>
      ox(m1)=x(m1)              !memory of old path
      oy(m1)=y(m1)              !memory of old path
      oz(m1)=z(m1)              !memory of old path
      ox(m2)=x(m2)              !memory of old path
      oy(m2)=y(m2)              !memory of old path
      oz(m2)=z(m2)              !memory of old path
      oo(m)=ica(m)              !memory of old path
      oo(m1)=ica(m1)            !memory of old path
      oo(m2)=ica(m2)            !memory of old path
      nx(m1)=x(m)+vx(nn(m))     !memory of new path
      ny(m1)=y(m)+vy(nn(m))     !memory of new path
      nz(m1)=z(m)+vz(nn(m))     !memory of new path
      nx(m2)=nx(m1)+vx(nn(m1))  !memory of new path
      ny(m2)=ny(m1)+vy(nn(m1))  !memory of new path
      nz(m2)=nz(m1)+vz(nn(m1))  !memory of new path

c     change conformation to new path -------------->
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      x(m2)=nx(m2)
      y(m2)=ny(m2)
      z(m2)=nz(m2)
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      ica(m2)=nn(m2)
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)

      if(look(m,m3))then       ! check excluded volumn for passage of [i,m]
c     calculate E_new--------------->
!$acc kernels
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
! !$acc end kernels
         Enew=EHB(m,m3,1)+ESHORT(m,m3,1) !use ifs
         do kkk=m,m3
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
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(m2)=oo(m2)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         Eold=EHB(m,m3,-1)+ESHORT(m,m3,-1)

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
            bNa(4)=bNa(4)+1
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
            do kkk=m,m3
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
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            ica(m2)=nn(m2)
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
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(m2)=oo(m2)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
      endif                     !for look(i,m)
113   continue
      bNt(4)=bNt(4)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move3d finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
