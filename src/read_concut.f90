      subroutine read_concut
	use params
      use concutt
      implicit integer(i-z)
      character*3 NAME
      dimension cut(0:19,0:19),cut_dev(0:19,0:19)
!      common/concutt/concut(0:19,0:19),concut2(0:19,0:19),concut_sc
      common/zscore/izscore

      rewind(1)
      read(1,*)
      do i=0,19
         read(1,*)NAME,(cut(i,j),j=0,19)
      enddo
      read(1,*)
      do i=0,19
         read(1,*)NAME,(cut_dev(i,j),j=0,19)
      enddo

      do i=0,19
         do j=0,19
            if(izscore.eq.1)then
               concut(i,j)=7
            elseif(izscore.eq.2)then
               concut(i,j)=7.5
            else
               concut(i,j)=cut(i,j)+cut_dev(i,j)*2.5 !real cut-off
            endif
            concut(i,j)=concut(i,j)/0.87 !cut-off on lattice
            concut2(i,j)=concut(i,j)**2 !cut-off squared on lattice
         enddo
      enddo
!$acc update device (concut,concut2)
      return
      end
