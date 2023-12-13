      subroutine read_E14
	use params
      use lengths
      use seqe
      use shortcom
      use short
      implicit integer(i-z)
!      parameter(ndim=1500)
!                parameter(nvec=416)
!      COMMON/short/ IBIN(-300:300),asr(ndim,-12:12)
!      common/lengths/Lch,Lch1,Lch2
!      common/seqe/seq(ndim),sec(ndim)
      common/forpreparemove4/asrr(0:19,0:19,-12:12)
      DIMENSION asrh(0:19,0:19,-12:12),asre(0:19,0:19,-12:12)
!      common/shortcom/eh3,es4,es5,es6,es7,es7a,es7b,es7c

ccccccccc read asrr,asrh,asre------------------>
      do i=0,19                 !asrr(Ai,Bi,dis) from 'r14.comm'
         do j=0,19
            read(12,*)
            read(12,*) (asrr(i,j,k),k=-12,-5)
            read(12,*) (asrr(i,j,k),k=-4,3) !without k=4
            read(12,*) (asrr(i,j,k),k=5,12)
            do k=4,1,-1
               asrr(i,j,k)=asrr(i,j,k-1) !without k=0
            enddo
         enddo
      enddo
      do i=0,19                 !asrh(Ai,Bi,dis) from 'r14h.comm'
         do j=0,19
            read(7,*)
            read(7,*) (asrh(i,j,k),k=-12,-5)
            read(7,*) (asrh(i,j,k),k=-4,3)
            read(7,*) (asrh(i,j,k),k=5,12)
            do k=4,1,-1
               asrh(i,j,k)=asrh(i,j,k-1)
            enddo
         enddo
      enddo
      do i=0,19                 !asre(Ai,Bi,dis) from 'r14e.comm'
         do j=0,19
            read(8,*)
            read(8,*) (asre(i,j,k),k=-12,-5)
            read(8,*) (asre(i,j,k),k=-4,3)
            read(8,*) (asre(i,j,k),k=5,12)
            do k=4,1,-1
               asre(i,j,k)=asre(i,j,k-1)
            enddo
         enddo
      enddo
c^^^^^^^^^ read asrr,asrh,asre finished ^^^^^^^^^^^^^^^^^

      do i=1,Lch-3
         do k=-12,12
            asr(i,k)=asrr(seq(i+1),seq(i+2),k) !general
            if(sec(i+1).eq.2.AND.sec(i+2).eq.2.AND.sec(i+3).eq.2)then
               if(sec(i).eq.2) then !helix
c                  asr(i,k)=(asr(i,k)+asrh(seq(i+1),seq(i+2),k))/2.0
                  asr(i,k)=asrh(seq(i+1),seq(i+2),k)
               endif
            endif
            if(sec(i+1).eq.4.AND.sec(i+2).eq.4.AND.sec(i+3).eq.4)then
               if(sec(i).eq.4) then !sheet
c                  asr(i,k)=(asr(i,k)+asre(seq(i+1),seq(i+2),k))/2.0
                  asr(i,k)=asre(seq(i+1),seq(i+2),k)*1.5
               endif
            endif
         enddo
      enddo
c^^^^^^^^^^^^ asr(i,ibin(r14)) finished ^^^^^^^^^^^^^^^^^^^^^
c     r(i,i+3)=k, E=asr(i,k), 12 bins (24 bins when considering chiral)
      do i=1,300
         kk=int((sqrt(float(i))*0.87))+1
         if(kk.gt.12) kk=12
         IBIN(I) = kk           !convert lattice r^2 into real r
         IBIN(-I)=-kk
      ENDDO
      IBIN(0)=IBIN(1)
!$acc update device(IBIN,asr(:,:))
c^^^^^^^^^^^^^^^^^ read E14 finished ^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
