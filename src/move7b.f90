      subroutine move7b
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
      use backup1
      use ehbc
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
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
!      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)
      common/shift/u21(nvec,nvec),m12(nvec),u1(nvec,100),u2(nvec,100)

!      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
!      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)
      common/ranzy/nozy

      mm=12

      do i=1,lenf
         i21=2+int(aranzy(nozy)*(Lch-6)) ![2,lenf-5]
         if(u21(ica(i21),ica(i21+1)).ne.0)goto 31
      enddo
      goto 33                   !can not find a two-bond can merg into one-bond

 31   i0=min(lenf2,i21+mm)
      do i=1,5
         i12=i21+3+int(aranzy(nozy)*(i0-i21-2.0001)) ![i21+3,i0]
         nc=m12(ica(i12))       !number of ways of extending 1-->2
         if(nc.gt.0)goto 32
      enddo
      goto 33                !can not find a one-bond can extend into two-bond

c     ------------- 2 -> 1 ----------------------------
 32   nn(i21)=u21(ica(i21),ica(i21+1))
      if(.not.goodc(ica(i21-1),nn(i21)))goto33
      if(.not.goodc(nn(i21),ica(i21+2)))goto33
c     ------------- 1 -> 2 ----------------------------
      p=int(aranzy(nozy)*nc)+1
      nn(i12-1)=u1(ica(i12),p)
      nn(i12)=u2(ica(i12),p)
      if(.not.goodc(ica(i12-1),nn(i12-1)))goto33
      if(.not.goodc(nn(i12),ica(i12+1)))goto33

c     backup old path ------------>
      oo(i21)=ica(i21)
!$acc data present(nop(:),nom(:),noa(:),noaa(:),nomm(:),
!$acc&              nopp(:),envir(:,:,:,:,:),oo,ox,oy,oz,
!$acc&              nn,nx,ny,nz)
!$acc parallel loop
      do i=i21+1,i12
         ox(i)=x(i)
         oy(i)=y(i)
         oz(i)=z(i)
         oo(i)=ica(i)
      enddo
c     prepare new path ------------>
      do i=i21+1,i12-2
         nn(i)=ica(i+1)
      enddo
      nx(i21+1)=x(i21)+vx(nn(i21))
      ny(i21+1)=y(i21)+vy(nn(i21))
      nz(i21+1)=z(i21)+vz(nn(i21))
!$acc kernels 
      do i=i21+2,i12
         nx(i)=nx(i-1)+vx(nn(i-1))
         ny(i)=ny(i-1)+vy(nn(i-1))
         nz(i)=nz(i-1)+vz(nn(i-1))
      enddo
!$acc end kernels
c     change conformation to new path -------------->
      ica(i21)=nn(i21)
!$acc parallel loop
      do i=i21+1,i12
         x(i)=nx(i)
         y(i)=ny(i)
         z(i)=nz(i)
         ica(i)=nn(i)
      enddo
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)

      m1=i21       !fixed
      m2=i12+1     !fixed
      if(look(m1,m2))then       ! check excluded volumn for passage of [m1,m2]
c     calculate E_new--------------->
!$acc kernels
         do pp=2,lenf1
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m1,m2,1)+ESHORT(m1,m2,1) !use ifs

         do kkk=m1,m2
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo
!$acc end kernels
c     return back the conformation and calculate E_old --------->
         ica(i21)=oo(i21)
!$acc parallel loop
         do i=i21+1,i12
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
            ica(i)=oo(i)
         enddo
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         Eold=EHB(m1,m2,-1)+ESHORT(m1,m2,-1)

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
c     E=energ
c     Et=E+de

         call metro(dE,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(14)=bNa(14)+1
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
            do kkk=m1,m2
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c            energ=energ+de	
            ica(i21)=nn(i21)
!$acc parallel loop
            do i=i21+1,i12
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
         ica(i21)=oo(i21)
!$acc parallel loop
         do i=i21+1,i12
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
            ica(i)=oo(i)
         enddo
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
      endif                     !for look(i,m)
!$acc end data
 33   continue
      bNt(14)=bNt(14)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move7b finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
