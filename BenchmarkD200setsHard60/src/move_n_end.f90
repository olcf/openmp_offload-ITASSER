      subroutine move_n_end
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
      use ehbc
      use backup1
      use shape
      implicit integer(i-z)
!      parameter(ndim=1500)
!                parameter(nvec=416)
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

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/ranzy/nozy

      m=int(aranzy(nozy)*(Mend_N-1))+1 !m=[1,Mend_N-1], [1,m] will be moved.
      
c     back-up old path---------->
      do i=1,m
         ox(i)=x(i)
         oy(i)=y(i)
         oz(i)=z(i)
         oo(i)=ica(i)
      enddo
c     prepare new path and move to new path -------------->
!$acc data present(nop(:),nom(:),noa(:),noaa(:),nomm(:),
!$acc&              nopp(:),envir(:,:,:,:,:))
! !$acc kernels
      do i=m,1,-1
 111     iv=int(aranzy(nozy)*anvec)+1
         if(.not.goodc(iv,ica(i+1))) goto 111
         x(i)=x(i+1)-vx(iv)
         y(i)=y(i+1)-vy(iv)
         z(i)=z(i+1)-vz(iv)
         ica(i)=iv
         nx(i)=x(i)             !memory of new conformation
         ny(i)=y(i)             !memory of new conformation
         nz(i)=z(i)             !memory of new conformation
         nn(i)=ica(i)           !memory of new conformation
      enddo
! !$acc end kernels
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)
      
      if(look(1,m+1))then       ! check excluded volumn for passage of [2,m+1]
c     calculate E_new--------------->
!$acc parallel loop
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(1,m+1,1)+ESHORT(1,m+1,1) !use ifs
!$acc parallel loop
         do kkk=1,m+1
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo
         
c     return back the conformation and calculate E_old --------->
         do i=1,m
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
            ica(i)=oo(i)
         enddo
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         
         Eold=EHB(1,m+1,-1)+ESHORT(1,m+1,-1) !cut off old nop by istat=-1

c     calculate eprofn while dord was calculated when call EHB(i,m,1)--->
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
            bNa(11)=bNa(11)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
!$acc parallel loop
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
c            energ=energ+de
            do i=1,m
               x(i)=nx(i)
               y(i)=ny(i)
               z(i)=nz(i)
               ica(i)=nn(i)
            enddo
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
!$acc parallel loop
            do kkk=1,m+1
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev
            call get_center     !update (amx,amy,amz)
c            call get_axis
         endif
      else                      !for look(i,m) loop
         do i=1,m
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
            ica(i)=oo(i)
         enddo
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
      endif                     !for look(i,m)
!$acc end data
113   continue
      bNt(11)=bNt(11)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move_n_end finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
