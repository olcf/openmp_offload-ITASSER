      subroutine move8
	use params
      !use openacc
      use backup2
      use chainm
      use chain1
      use echain1
      use lengths
      use envir1
      use order
      use seqe
      use short
      use shape
      use ehbc
      use backup1
      use tempArrays
      implicit integer(i-z)
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
      common/maxdis2/maxdis2(ndim)
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
!      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)
      common/shift/u21(nvec,nvec),m12(nvec),u1(nvec,100),u2(nvec,100)

      common/move31/v31(-15:15,-15:15,-15:15,6)
      common/move32/v32(-15:15,-15:15,-15:15,6)
      common/move33/v33(-15:15,-15:15,-15:15,6)
      common/move34/Np3(-15:15,-15:15,-15:15)

!      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
!      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)
      common/ranzy/nozy
      common/nwmax1/nvecr

c     q bonds will be moved, q is in [Mran1,Mran2]
      Mran1=4
      M_ran2=6
      q=int(aranzy(nozy)*(Mran2-Mran1+1))+Mran1
      q3=q-3

c     choose the position (m1) to be moved
c     m1>1, m2<lenf, otherwise can not do goodc()!
c     [m1+1,m2-1] will be really moved ---------------->
 111  m1=int(aranzy(nozy)*(lenf-q-2.00001))+2 !m1=[2,lenf-q-1]
      m2=m1+q                   !m2<=lenf-1

      mx=x(m2)-x(m1)
      my=y(m2)-y(m1)
      mz=z(m2)-z(m1)
cccccccccccccccFirst q-3 bonds cccccccccccccccccccccccccc
      nn(m1-1)=ica(m1-1)
      do i=1,q3
         uu1=0
 112     uu=1+int(nvecr*aranzy(nozy)) ! uu in [1,nvec]
         uu1=uu1+1
         if(uu1.gt.nvecr)goto 111 ! path maybe wrong!
         xx=mx-vx(uu)
         yy=my-vy(uu)
         zz=mz-vz(uu)

         xx2=xx*xx
         if(xx2.ge.maxdis2(q-i))goto 112
         yy2=yy*yy
         if(yy2.ge.maxdis2(q-i))goto 112
         rr=xx2+yy2+zz*zz   ! rr=xx*xx+yy*yy+zz*zz
         if(rr.ge.maxdis2(q-i))goto 112 !too distant for remaining walks.
         if(rr.lt.12)goto 112 !too close for overlap
         if(.not.goodc(uu,nn(m1+i-2)))goto 112 !check neighbor

         nn(m1+i-1)=uu
         mx=xx
         my=yy
         mz=zz
      enddo
      if(mx.gt.12.or.mx.lt.-12) goto 111 !no defined vector in 3-bond move
      if(my.gt.12.or.my.lt.-12) goto 111 !no defined vector in 3-bond move
      if(mz.gt.12.or.mz.lt.-12) goto 111 !no defined vector in 3-bond move
      
ccccccccccccccccccc Last 3 bonds ccccccccccccccccccccccccc
      nc=Np3(mx,my,mz)          !number of pathes from i to i+m ???
      if(nc.eq.0)then
         write(*,*)'absent q-movement in move5',q,mx,my,mz,m1,m2
         goto 114
      endif
      uu1=0
 113  p=int(aranzy(nozy)*(nc-0.00001))+1    ![1,nc]
      uu1=uu1+1
      if(uu1.gt.nc)goto 111
      nn(m2-3)=v31(mx,my,mz,p)
      if(.not.goodc(nn(m2-4),nn(m2-3))) goto 113
      nn(m2-1)=v33(mx,my,mz,p)  !???
      if(.not.goodc(nn(m2-1),ica(m2))) goto 113
      nn(m2-2)=v32(mx,my,mz,p)  !???

c     backup old path ------------>
!$acc data present(nop(:),nom(:),noa(:),noaa(:),nomm(:),
!$acc&              nopp(:),envir(:,:,:,:,:)) copyin(vx,vy,vz)
!$acc kernels
      do i=1,q
         ox(m1+i)=x(m1+i)
         oy(m1+i)=y(m1+i)
         oz(m1+i)=z(m1+i)
         oo(m1+i-1)=ica(m1+i-1)
      enddo
!$acc end kernels
c     prepare new path ------------>
      nx(m1+1)=x(m1)+vx(nn(m1))
      ny(m1+1)=y(m1)+vy(nn(m1))
      nz(m1+1)=z(m1)+vz(nn(m1))
! !$acc data copyin(vx,vy,vz)
!$acc kernels  
      do i=2,q
         nx(m1+i)=nx(m1+i-1)+vx(nn(m1+i-1))
         ny(m1+i)=ny(m1+i-1)+vy(nn(m1+i-1))
         nz(m1+i)=nz(m1+i-1)+vz(nn(m1+i-1))
      enddo
! !$acc end data
! !$acc end kernel
c     change conformation to new path -------------->
! !$acc kernels
      do i=1,q
         x(m1+i)=nx(m1+i)
         y(m1+i)=ny(m1+i)
         z(m1+i)=nz(m1+i)
         ica(m1+i-1)=nn(m1+i-1)
      enddo
!$acc end kernels
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)
      
      if(look(m1,m2))then       ! check excluded volumn for passage of [i,m]
c     calculate E_new--------------->
!$acc kernels
         do pp=2,lenf1
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
c            nhbn(pp)=nhbnn(pp)
         enddo
         Enew=EHB(m1,m2,1)+ESHORT(m1,m2,1) !use ifs
         do kkk=m1,m2
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo
!$acc end kernels
c     return back the conformation and calculate E_old --------->
!$acc parallel loop
         do i=1,q
            x(m1+i)=ox(m1+i)
            y(m1+i)=oy(m1+i)
            z(m1+i)=oz(m1+i)
            ica(m1+i-1)=oo(m1+i-1)
         enddo
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         Eold=EHB(m1,m2,-1)+ESHORT(m1,m2,-1)

c     calculate eprofn while dord was calculated when call EHB(i,m,1)--->
         eprofn=0.0
!$acc parallel loop
         do pp=2,lenf1
            is=seq(pp)
            ia=noa(pp)        !noa is now 
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+eprofn-eprofo
c         E=energ
c         Et=E+de

         call metro(dE,id) !id=1, accepted; id=3, rejected

         if(id.eq.3)then        !rejected
            icnt=icnto
            sumct=sumcto 		
         else                   !accepted
            bNa(7)=bNa(7)+1
            bN5a(q)=bN5a(q)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
!$acc kernels
            do pp=2,lenf1
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
c               nhbnn(pp)=nhbn(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m1,m2
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
!$acc end kernels
            codevsum=codev
            didevsum=didev	
c            energ=energ+de	
!$acc parallel loop
            do i=1,q
               x(m1+i)=nx(m1+i)
               y(m1+i)=ny(m1+i)
               z(m1+i)=nz(m1+i)
               ica(m1+i-1)=nn(m1+i-1)
            enddo
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            call get_center       !update (amx,amy,amz)
c            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
!$acc parallel loop
         do i=1,q
            x(m1+i)=ox(m1+i)
            y(m1+i)=oy(m1+i)
            z(m1+i)=oz(m1+i)
            ica(m1+i-1)=oo(m1+i-1)
         enddo
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
      endif                     !for look(i,m)
!$acc end data
 114  continue
      bNt(7)=bNt(7)+1
      bN5t(q)=bN5t(q)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move8 finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
