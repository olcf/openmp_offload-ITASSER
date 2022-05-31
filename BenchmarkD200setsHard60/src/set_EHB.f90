      subroutine set_EHB
	use params
      use lengths
      use ENERGY
      use seqe
      implicit integer(i-z)
!      parameter(ndim=1500)
!      parameter(nvec=416)
!      COMMON/ENERGY/EH5,ES3,ES3a,ES3b,EH1,EH1a,EH4,EHBIJ(ndim,ndim)
!      common/lengths/Lch,Lch1,Lch2
!      common/seqe/seq(ndim),sec(ndim)

c
c     EHBIJ - set-up secondary structure dependent
c     strength of the hyrogen bond network - stronger for helices
c     and beta-sheets
c

      do i=1,Lch
         is=sec(i)
         do j=1,Lch
            js=sec(j)
            EHBIJ(i,j)=1
            if(iabs(i-j).eq.3.and.is.eq.2.and.js.eq.2)then
               EHBIJ(i,j)=EHBIJ(i,j)+0.5 !helix, enhanced
            endif
            if(is.eq.4.or.js.eq.4) then
               if(is*js.ne.8.and.iabs(i-j).gt.4)then
                  EHBIJ(i,j)=EHBIJ(i,j)+0.5 !beta-beta, enhanced
               endif
            endif
         enddo
      enddo

c^^^^^^^^^^^ set H_bond finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
