      subroutine move6
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
      use shape
      use ehbc
      use backup1
      use tempArrays

      implicit integer(i-z)
!      parameter(ndim=1500)
!      parameter(nrep=100)       !number of replicas
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

!      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
!      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/ranzy/nozy
      common/lattice/m_latt,latt1,latt2

c      nob=6+int(aranzy(nozy)*6.999) !number of moved bonds, [6,12]
      nob=6

c     choose the position (m) to be moved
      i=int(aranzy(nozy)*nfl6)+1 ![1,nfl6]
      m=ras6(i)
      m1=m+1
      mnob=m+nob
      mnob1=mnob-1
      mnob2=mnob-2

 201  mx=int(aranzy(nozy)*2.999999)-1 ![-1,0,1]
      my=int(aranzy(nozy)*2.999999)-1 ![-1,0,1]
      mz=int(aranzy(nozy)*2.999999)-1 ![-1,0,1]
      if((mx*mx+my*my+mz*mz).eq.0) goto 201

      wx=vx(ica(m))+mx
      wy=vy(ica(m))+my
      wz=vz(ica(m))+mz
      ir=wx*wx+wy*wy+wz*wz
      if(ir.lt.latt1)goto 201
      if(ir.gt.latt2)goto 201
      nn(m)=vector(wx,wy,wz)    !new vector of vertix 'm'
      if(.not.goodc(ica(m-1),nn(m))) goto 202
      if(.not.goodc(nn(m),ica(m1))) goto 202

      wx=vx(ica(mnob1))-mx
      wy=vy(ica(mnob1))-my
      wz=vz(ica(mnob1))-mz
      ir=wx*wx+wy*wy+wz*wz
      if(ir.lt.latt1) goto 202
      if(ir.gt.latt2) goto 202
      nn(mnob1)=vector(wx,wy,wz) !new vector of vertix 'm+5'
      if(.not.goodc(ica(mnob2),nn(mnob1))) goto 202
      if(.not.goodc(nn(mnob1),ica(mnob))) goto 202

c     prepare new path and backup old path ------------>
!$acc data present(afs(:),afsn(:),noa(:),nop(:),nom(:),
!$acc&              noaa(:),nomm(:),nopp(:),envir(:,:,:,:,:),oy,oz,
!$acc&              nn,nx,ny,nz,oo,ox)
!$acc parallel loop
      do i=m1,mnob1
         ox(i)=x(i)             !memory of old path
         oy(i)=y(i)             !memory of old path
         oz(i)=z(i)             !memory of old path
         nx(i)=x(i)+mx          !memory of new path
         ny(i)=y(i)+my          !memory of new path
         nz(i)=z(i)+mz          !memory of new path
      enddo
      oo(m)=ica(m)              !memory of old path
      oo(mnob1)=ica(mnob1)      !memory of old path

c     change conformation to new path -------------->
      do i=m1,mnob1
         x(i)=nx(i)
         y(i)=ny(i)
         z(i)=nz(i)
      enddo
      ica(m)=nn(m)
      ica(mnob1)=nn(mnob1)
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)

      if(look(m,mnob))then      ! check excluded volumn for passage of [i,m]
c     calculate E_new--------------->
!$acc kernels
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m,mnob,1)+ESHORT(m,mnob,1) !use ifs
         do kkk=m,mnob
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo
c     return back the conformation and calculate E_old --------->
         do i=m1,mnob1
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
         enddo
!$acc end kernels
         ica(m)=oo(m)
         ica(mnob1)=oo(mnob1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         Eold=EHB(m,mnob,-1)+ESHORT(m,mnob,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
!$acc kernels
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
            bNa(10)=bNa(10)+1
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
            do kkk=m,mnob
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
!$acc end kernels
            codevsum=codev
            didevsum=didev	
c            energ=energ+de	
!$acc parallel loop
            do i=m1,mnob1
               x(i)=nx(i)
               y(i)=ny(i)
               z(i)=nz(i)
            enddo
            ica(m)=nn(m)
            ica(mnob1)=nn(mnob1)
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            call get_center     !update (amx,amy,amz)
c            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
!$acc parallel loop
         do i=m1,mnob1
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
         enddo
         ica(m)=oo(m)
         ica(mnob1)=oo(mnob1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
      endif                     !for look(i,m)
!$acc end data
 202  continue
      bNt(10)=bNt(10)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move6 finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
