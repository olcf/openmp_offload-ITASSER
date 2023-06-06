      subroutine get_center
	use params
      use chainm
      use chain1
      use lengths
      use echain1
      use sg
      use echain2
      use echain4
      use seqe
      use echain5
      use eigen
      implicit integer(i-z)
!      parameter(ndim=1500)
!      parameter(nrep=100)
!      parameter(nvec=416)
!      common/lengths/Lch,Lch1,Lch2
      common/initialinput/switch,k_cycle,k_phot,N_ann

      common/chain0/ras(ndim),nfl
!      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
!      common/echain1/ex(ndim),ey(ndim),ez(ndim) !CA
!      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
!      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
!      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
!      common/chainm/mv(ndim)
!      common/seqe/seq(ndim),sec(ndim)
!      common/sg/gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      
      common/center/cex,cey,cez
!      common/eigen/AA(3,3),EE(3),HH(3,3)

*********************************************************
*     record centriod------------------------->
      cex=0
      cey=0
      cez=0
      if(switch.eq.1)then
         do i=1,Lch
            cex=cex+x(i)
            cey=cey+y(i)
            cez=cez+z(i)
         enddo
      else                      !!!!!!!!template based run
         do i=1,Lch
            if(mv(i).gt.0)then
               cex=cex+x(i)
               cey=cey+y(i)
               cez=cez+z(i)
            else
               cex=cex+ex(i)
               cey=cey+ey(i)
               cez=cez+ez(i)
            endif
         enddo
      endif
      cex=cex/float(Lch)
      cey=cey/float(Lch)
      cez=cez/float(Lch)
      
*********************************************************************
*     calculate axis ----------------------->
      AA(1,1)=0
      AA(1,2)=0
      AA(1,3)=0
      AA(2,2)=0
      AA(2,3)=0
      AA(3,3)=0
      do i=1,Lch
         if(mv(i).gt.0)then
            agxi=x(i)+GX(ica(i-1),ica(i),seq(i))-cex
            agyi=y(i)+GY(ica(i-1),ica(i),seq(i))-cey
            agzi=z(i)+GZ(ica(i-1),ica(i),seq(i))-cez
         else
            agxi=egx(i)-cex
            agyi=egy(i)-cey
            agzi=egz(i)-cez
         endif
         AA(1,1)=AA(1,1)+agxi*agxi
         AA(1,2)=AA(1,2)+agxi*agyi
         AA(1,3)=AA(1,3)+agxi*agzi
         AA(2,2)=AA(2,2)+agyi*agyi
         AA(2,3)=AA(2,3)+agyi*agzi
         AA(3,3)=AA(3,3)+agzi*agzi
      enddo
      AA(1,1)=AA(1,1)/float(Lch)
      AA(1,2)=AA(1,2)/float(Lch)
      AA(1,3)=AA(1,3)/float(Lch)
      AA(2,2)=AA(2,2)/float(Lch)
      AA(2,3)=AA(2,3)/float(Lch)
      AA(3,3)=AA(3,3)/float(Lch)
      AA(2,1)=AA(1,2)
      AA(3,1)=AA(1,3)
      AA(3,2)=AA(2,3)

c     E(1)=<x^2>=Ax^2/5; E(2)=<y^2>=Ay^2/5; E(3)=<z^2>=Az^2/5
      call eigenvalue           !get rotation matrix TT

! !$acc update device(AA)

c^^^^^^^^^^^^^^^^^ (amx,amy,amz) finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
