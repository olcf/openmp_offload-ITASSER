      subroutine initial_move
	use params
      use chainm
      use backup2
      use chain1
      use lengths
      use seqe
      use echain1
      use one
      use sg
      use echain2
      use echain4
      use order
      use echain5
      use envir1
      use short
      use shape
      use backup1
      use ehbc
      use eigen
      implicit integer(i-z)
!      parameter(ndim=1500)
!      parameter(nrep=100)
!      parameter(nvec=416)
!      common/lengths/Lch,Lch1,Lch2
!      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
!      common/seqe/seq(ndim),sec(ndim)
!      common/sg/gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
!      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
!      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1
!      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
!      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
!      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
!      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
!      common/backup2/eprofo,eprofn,energ
!      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev

!      common/echain1/ex(ndim),ey(ndim),ez(ndim) !CA
!      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
!      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
!      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
!      common/chainm/mv(ndim)

      dimension axx(ndim),ayy(ndim),azz(ndim)

*     check excluded volumn
      call get_vvv
*     move coordinates of the chain to its center of mass ---->
      call get_center

******* get template ready *********************
      

*     calculate (amx,amy,amz) -------------->
c      call get_axis

cccccccc prepare all the initial parameters of this replica ccccccccc
      icnto=0                   !need when ISTAT>0
      sumcto=0                  !need when ISTAT>0
      do i=1,Lch
         nop(i)=0               !number of parallel contacts
         noa(i)=0
         nom(i)=0
      enddo

      energ=EHB(1,Lch,1)+ESHORT(1,Lch,10) !initialize all below parameters
ccc   eprofo,eprofn: calculated each time;
ccc   icnt: always start from icnto; (need initial)
ccc   sumct: always start from sumcto; (need initial)
ccc   nop(): start from 0, or nopp; (need initial)
ccc   noa(): start from 0, or noaa; (need initial)
ccc   nom(): start from 0, or nomm; (need initial)
ccc   codevsum: only used when istat=-1, i.e. in Enew.
ccc   didevsum: only used when istat=-1, i.e. in Enew.
ccc   afs(): keep in store.

cccc  backup all the parameters, calculated in EHB() and ESHORT() cccccccccc
      eprofo=0.0
      do k=1,Lch
         is=seq(k)
         ia=noa(k)              !number of antiparallel contact-apirs
         ip=nop(k)              !number of parallel contact-apirs ^^
         im=nom(k)              !number of orgonal contact-apirs
         eprofo=eprofo+envir(ia,im,ip,is,3)
      enddo

      icnto=icnt                !backup of contact-order
      sumcto=sumct              !backup of contact-order
      do i=1,Lch
         nopp(i)=nop(i)         !backup of number of contact-pair
         noaa(i)=noa(i)         !backup of number of contact-pair
         nomm(i)=nom(i)         !backup of number of contact-pair
      enddo
      codevsum=conew            !backup  panelity for total deviation of comb
      didevsum=dinew            !backup  panelity for total deviation of dist
!$acc update device(nopp,noaa,nomm)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c^^^^^^^^^^^^^^^^^ initialization finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
