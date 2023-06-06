      subroutine move9
	use params
      !use openacc
      use backup2
      use chainm
      use chain1
      use echain1
      use lengths
      use order
      use seqe
      use one
      use envir1
      use short
      use shape
      use ehbc
      use backup1
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
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
!      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
!      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)

      dimension fax(ndim),fay(ndim),faz(ndim)
!      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
!      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/ranzy/nozy
      common/nwmax1/nvecr

******backup old conformation --------------->
! !$acc kernels
      do i=1,lenf
         ox(i)=x(i)
         oy(i)=y(i)
         oz(i)=z(i)
         oo(i)=ica(i)
      enddo
! !$acc end kernels
******prepare new confromation ---------------->
c     (x,y,z)->(fax,fay,faz) based on [x(k),y[k],z[k]]:
      k=int(lenf*aranzy(nozy))+1   ![1,lenf]
      ddxyz=1.5
 99   fdx=(aranzy(nozy)*ddxyz*2.0)-ddxyz ![-1.5, 1.5]
      fdy=(aranzy(nozy)*ddxyz*2.0)-ddxyz
      fdz=(aranzy(nozy)*ddxyz*2.0)-ddxyz
      ar=fdx*fdx+fdy*fdy+fdz*fdz
      if(ar.lt.0.75) go to 99
!$acc data present(nop(:),nom(:),noa(:),noaa(:),nomm(:),
!$acc&              nopp(:),envir(:,:,:,:,:))
!$acc kernels
      do i=1,lenf
         ar=sqrt(float((x(k)-x(i))**2+(y(k)-y(i))**2+(z(k)-z(i))**2))
         fax(i)=x(i)+fdx*(1.0-ar/(1.5*acrit))
         fay(i)=y(i)+fdy*(1.0-ar/(1.5*acrit))
         faz(i)=z(i)+fdz*(1.0-ar/(1.5*acrit))
      enddo
!$acc end kernels
******project the conformation onto lattices, i.e. (fax,fay,faz)->(nx,ny,nz):
      px=nint(fax(1))
      py=nint(fay(1))
      pz=nint(faz(1))
      nx(1)=px
      ny(1)=py
      nz(1)=pz
! !$acc kernels loop
      DO 101 i=2,lenf
         armin=10000.
         kk=0
! !$acc kernels
         do 1009 k=1,nvecr
            if(i.ne.2) then
               if(.not.goodc(kkkk,k))goto 1009
            endif
            kx=px+vx(k)
            ky=py+vy(k)
            kz=pz+vz(k)
            bx=fax(i)-float(kx)
            by=fay(i)-float(ky)
            bz=faz(i)-float(kz)
            ar=bx*bx+by*by+bz*bz
            if(ar.lt.armin) then
               kk=k
               mx=kx
               my=ky
               mz=kz
               armin=ar
            endif
! !$acc end kernels
 1009    continue
         kkkk=kk
         if(kk.EQ.0) then
            write(20,*)' 	ERROR in the CHAIN global move' 
            GO TO  113          !do not move
         endif
         nx(i)=mx
         ny(i)=my
         nz(i)=mz
         px=mx
         py=my
         pz=mz
 101  continue
***************** new ica(i):
      ic=0
      do i=1,lenf1
         j=i+1
         wx=nx(j)-nx(i)
         wy=ny(j)-ny(i)
         wz=nz(j)-nz(i)
         nn(i)=vector(wx,wy,wz)
         if(nn(i).ne.oo(i)) ic=ic+1
      enddo
      if(ic.eq.0) go to 113     !do not move
c^^^^^^^^^^^^^^^^^^ new conformation obtained ^^^^^^^^^^^^^^^^^^^^^^^^^

c     change conformation to new path -------------->
      do i=1,lenf
         x(i)=nx(i)
         y(i)=ny(i)
         z(i)=nz(i)
         ica(i)=nn(i)
      enddo
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)
      
      if(look(2,lenf1))then
c     calculate E_new--------------->
!$acc kernels
         do pp=2,lenf1
            nop(pp)=nopp(pp)
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(2,lenf1,1)+ESHORT(2,lenf1,1)
         do kkk=2,lenf1
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo
!$acc end kernels
c     return back the conformation and calculate E_old --------->
         do i=1,lenf
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
            ica(i)=oo(i)
         enddo
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         Eold=EHB(2,lenf1,-1)+ESHORT(2,lenf1,-1)
      
c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
!$acc parallel loop
         do pp=2,lenf1
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
            bNa(15)=bNa(15)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 	
!$acc parallel loop
            do pp=2,lenf1
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
!$acc parallel loop
            do kkk=2,lenf1
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c            energ=energ+de	
            do i=1,lenf
               x(i)=nx(i)
               y(i)=ny(i)
               z(i)=nz(i)
               ica(i)=nn(i)
            enddo
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            call get_center     !update (amx,amy,amz)
c            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         do i=1,lenf
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
            ica(i)=oo(i)
         enddo
! !$acc end data
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
      endif
 113  continue
!$acc end data
      bNt(15)=bNt(15)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move9 finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
