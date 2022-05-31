      subroutine sort_ras
	use params
      use lengths
      implicit integer(i-z)
!      parameter(ndim=1500)
!      parameter(nrep=10)
!      parameter(nvec=416)
!      common/lengths/Lch,Lch1,Lch2
      common/chain0/ras(ndim),nfl

      dimension ras0(ndim)

      ras_min=10000
      ras_max=-10000
      do i=1,nfl
         if(ras(i).gt.ras_max)ras_max=ras(i)
         if(ras(i).lt.ras_min)ras_min=ras(i)
      enddo

      nfl0=1
      ras0(nfl0)=ras_min
      do 1 while(ras0(nfl0).lt.ras_max)
         ras_min=10000
         do i=1,nfl
            if(ras(i).gt.ras0(nfl0))then
               if(ras(i).lt.ras_min)then
                  ras_min=ras(i)
               endif
            endif
         enddo
         nfl0=nfl0+1
         ras0(nfl0)=ras_min
 1    continue

      nfl=nfl0
      do i=1,nfl
         ras(i)=ras0(i)
      enddo

*^^^^^^^^^^^^^^ sort ras done ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
