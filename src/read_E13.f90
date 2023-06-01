      subroutine read_E13
	use params
      use lengths
      use shortcom
      use seqe
      use short
      implicit integer(i-z)
!      parameter(ndim=1500)
!                parameter(nvec=416)
!      common/short2/codevsum,didevsum,csr(ndim,2)
!      common/lengths/Lch,Lch1,Lch2
!      common/seqe/seq(ndim),sec(ndim)
!      common/shortcom/eh3,es4,es5,es6,es7,es7a,es7b,es7c
      dimension csre(0:19,0:19,2)

c     R13 potential - two bins only (helical and expanded)
c     r2<48, E=csre(i,j,1); r2>48, E=csre(i,j,2)
      do i=0,19
         do j=0,19
            read(5,*)
            read(5,*) (csre(i,j,k),k=1,2)
         enddo
      enddo

      do i=1,Lch2
         do k=1,2
            csr(i,k)=2.0*csre(seq(i),seq(i+2),k)
         enddo
      enddo

!$acc update device(csr)

c^^^^^^^^^^^^^^^^^ read E13 finished ^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
