      subroutine move2
      use params
      use backup2
      !!use openacc
      use chainm
      use chain1
      use echain1
      use lengths
      use envir1
      use order
      use seqe
      use short
      use shape
      use backup1
      use ehbc
      use tempArrays
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
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      double precision bNa,bNt,bNNa,bNNt
!      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev    
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)

!      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
!      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/ranzy/nozy

c     choose the position of (m) ----------->
      i=int(aranzy(nozy)*nfl2)+1  ![1,nfl2]
      m=ras2(i)
      m1=m+1
      m2=m+2

cccc  all the pairs have at least one another pair, so nc>=2
      nc=Np2(ica(m),ica(m1))    !number of pathes from m to m2

c     choose p'th new path from (m) to (m2) -------------->
 201  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      nn(m)=v21(ica(m),ica(m1),p)
      if(.not.goodc(ica(m-1),nn(m))) goto 201
      nn(m1)=v22(ica(m),ica(m1),p)
      if(.not.goodc(nn(m1),ica(m2))) goto 201
c^^^^^^^^^^ new conformation chosen finished ^^^^^^^^^^^^^^^^^^^^^^^^
      if(nn(m).eq.ica(m))goto 113

c     prepare new path and backup old path ------------>
      ox(m1)=x(m1)              !memory of old path
      oy(m1)=y(m1)              !memory of old path
      oz(m1)=z(m1)              !memory of old path
      oo(m)=ica(m)              !memory of old path
      oo(m1)=ica(m1)            !memory of old path
      nx(m1)=x(m)+vx(nn(m))     !memory of new path
      ny(m1)=y(m)+vy(nn(m))     !memory of new path
      nz(m1)=z(m)+vz(nn(m))     !memory of new path

c     change conformation to new path -------------->
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)

      if(look(m,m2))then        ! check excluded volumn for passage of [m,m2]
c     calculate E_new--------------->
!$acc data present(seq(:),nop(:),nom(:),noa(:)
!$acc& ,afs(:),noaa(:),nomm(:),envir(:,:,:,:,:)) !copy(eprofn)
!$acc kernels async(2)  
! !$acc loop gang(1024)
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
!$acc end kernels !Just added this 
         Enew=EHB(m,m2,1)+ESHORT(m,m2,1) !icnt,nop are repeated.
! !$acc loop vector(256)
!$acc kernels  async(2) ! Just added this one
         do kkk=m,m2
            afsn(kkk)=afs(kkk)  !afs(i) is useful in ESHORT
         enddo
!$acc end kernels
! !$acc end data
! !$acc update host(x,y,z)
c     return back the conformation and calculate E_old --------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         Eold=EHB(m,m2,-1)+ESHORT(m,m2,-1) !repeated part of icnt, nop removed

c     calculate eprofn while dord was calculated when call EHB(m,m2,1)--->
         eprofn=0.0
!$acc kernels loop gang reduction(+:eprofn) async(2) !gang I just took out
         do pp=1,Lch
            is=seq(pp)          
            ia=noa(pp)          
            ip=nop(pp)          !nopp(old)+1(new)-1(old)=new
            im=nom(pp)
c            if(ip.lt.0.or.ip.gt.15.or
c     &           .ia.lt.0.or.ia.gt.15.or
c     &           .im.lt.0.or.im.gt.15)then
c               write(*,*)'ia,ip,im=',ia,ip,im
c            endif
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo
!$acc end kernels
!$acc wait(2)
!$acc end data
c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c         E=energ
c         Et=E+de

         call metro(dE,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(2)=bNa(2)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt          !backup1
            sumcto=sumct        !backup2
! !$acc kernels ! loop gang(12) vector(1024) 
            do pp=1,Lch
               nopp(pp)=nop(pp) !backup3
               nomm(pp)=nom(pp) !backup4
               noaa(pp)=noa(pp) !backup5
            enddo
            eprofo=eprofn       !backup6
            codevsum=codev      !backup7
            didevsum=didev      !backup8
            do kkk=m,m2
               afs(kkk)=afsn(kkk) !backup9
            enddo
! !$acc end kernels
! !$acc end data
c            energ=energ+de
c     change the conformation to the new position--------->
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            call get_center     !update (amx,amy,amz)
c            call get_axis
         endif                  !for id
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
      endif                     !for look(i,m)
 113  continue
      bNt(2)=bNt(2)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move2 finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
