      subroutine read_E15
	use params
      use short1
      use lengths
      use seqe
      use shortcom
      implicit integer(i-z)
!      parameter(ndim=1500)
!                parameter(nvec=416)
!      COMMON/short1/ JBIN(0:500),bsr(ndim,16),acops(ndim,16)
!      common/lengths/Lch,Lch1,Lch2
!      common/seqe/seq(ndim),sec(ndim)
      DIMENSION bsrh(0:19,0:19,16)
      dimension bsre(0:19,0:19,16)
      dimension bsrr(0:19,0:19,16)
!      common/shortcom/eh3,es4,es5,es6,es7,es7a,es7b,es7c

cccccccc read bsrr,bsrh,bsre ----------------->
      do i=0,19                 !read bsrr from 'r15.dat'
         do j=0,19
            read(9,*)
            read(9,*) (bsrr(i,j,k),k=1,8)
            read(9,*) (bsrr(i,j,k),k=9,16)
         enddo
      enddo
      do i=0,19                 !read bsrh from 'r15h.dat'
         do j=0,19
            read(10,*)
            read(10,*) (bsrh(i,j,k),k=1,8)
            read(10,*) (bsrh(i,j,k),k=9,16)
         enddo
      enddo	
      do i=0,19                 !read bsre from 'r15e.dat'
         do j=0,19
            read(11,*)
            read(11,*) (bsre(i,j,k),k=1,8)
            read(11,*) (bsre(i,j,k),k=9,16)
         enddo
      enddo	

      do i=1,Lch-4
         do k=1,16
            bsr(i,k)=bsrr(seq(i+1),seq(i+3),k)
            if(sec(i+1).eq.2.AND.sec(i+2).eq.2.AND.sec(i+3).eq.2)then
c               bsr(i,k)=(bsr(i,k)+bsrh(seq(i+1),seq(i+3),k))/2.0 !helix
               bsr(i,k)=bsrh(seq(i+1),seq(i+3),k)
            endif
            if(sec(i+1).eq.4.AND.sec(i+2).eq.4.AND.sec(i+3).eq.4)then
c               bsr(i,k)=(bsr(i,k)+bsre(seq(i+1),seq(i+3),k))/2.0 !sheet
               bsr(i,k)=bsre(seq(i+1),seq(i+3),k)*1.5
            endif
         enddo
      enddo
c^^^^^^^^^^^^^^^^^^^^ E_15(Ai,Aj,dis) prepared ^^^^^^^^^^^^^^^^^^^

c     prepare distance bin-------------------->
      do i=0,500
         kk=int((sqrt(float(i))*0.87))+1 !i, lattice-dist; kk, real distance
         if(kk.gt.16) kk=16
         JBIN(I) = kk           !jbin: real distance
      ENDDO

ccccc acops(i,jbin) to enhance the contacts between gragments cccc
      do i=1,Lch-4
         acops(i,1)=(min(bsr(i,1),0.0))/2.0 !acpos<0
         do k=2,15
            acops(i,k)=min(0.0,bsr(i,k-1)+2.0*bsr(i,k)+bsr(i,k+1)) !<0
         enddo
         acops(i,16)=(min(bsr(i,16),0.0))/2.0
      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!$acc update device(bsr,acops,JBIN)

c^^^^^^^^^^^^^^^^^ read E15 finished ^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
