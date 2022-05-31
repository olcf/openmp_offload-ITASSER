      subroutine prepare_frg
	use params
      use lengths
      use fr
      use seqe
      implicit integer(i-z)
!      parameter(ndim=1500)
!                parameter(nvec=416)
!      common/fr/frga(ndim),frgb(ndim)
!      common/seqe/seq(ndim),sec(ndim)
!      common/lengths/Lch,Lch1,Lch2

      do i=1,Lch
         frga(i)=0.0
         frgb(i)=0.0
      enddo

      do i=1,Lch-7
         q=0
         do j=i,i+7
            if(sec(j).eq.2) q=q+2 !helix structure.
         enddo
         if(q.eq.16)then        !8 continue alpha-residues
            frga(i)=10.5/0.87   !distance for 7 alpha-bonds
         endif
      enddo

      do i=1,Lch-6
         q=0
         do j=i+1,i+5
            if(sec(j).eq.4) q=q+4 !beta structure
         enddo
         if(q.eq.20)then        !5 continue beta-residues
            if(sec(i).ne.2.and.sec(i+6).ne.2)then
               frgb(i)=19.1/0.87 !distance for 6 beta-bonds
            endif
         endif
      enddo

c      do i=1,Lch
c         write(*,*)i,seq(i),sec(i),frga(i),frgb(i)
c      enddo
!$acc update device(frga,frgb)
      return
      end
