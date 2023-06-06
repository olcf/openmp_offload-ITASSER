      subroutine read_seq
	use params
      use lengths
      use seqe
      use icgg
      implicit integer(i-z)
!      parameter(ndim=1500)
!                parameter(nvec=416)
      character*3 aa(-1:20), NAME,sequ
!      common/seqe/seq(ndim),sec(ndim)
!      common/lengths/Lch,Lch1,Lch2
!      common/icgg/ icg(ndim), EH6  
      common/aminoacid/sequ(ndim)

      data aa/ 'BCK','GLY','ALA','SER','CYS','VAL','THR','ILE',
c                -1    0     1     2     3     4     5     6
     &     'PRO','MET','ASP','ASN','LEU',
c            7     8     9    10    11
     &     'LYS','GLU','GLN','ARG',
c           12    13    14    15
     &     'HIS','PHE','TYR','TRP','CYX'/
c           16    17    18    19    20

      do 121 i=1,Lch
         read(14,707) k,NAME,SEC(I),tmp
         do j=0,19
            if(NAME.eq.aa(j)) then
               SEQ(i)=j
               icg(i)=0
               sequ(i)=name
               if(NAME.eq.'ASP'.or.NAME.eq.'GLU') icg(i)=-1
               if(NAME.eq.'LYS'.or.NAME.eq.'ARG') icg(i)= 1	
               go to 121
            endif
         enddo
         SEQ(i)=0
         icg(i)=0
         sequ(i)='GLY'
 121  continue
 707  format(i5,3x,a3,2i5)
      close(14)

c^^^^^^^^^^^^^^^^^^^^^ read sequence finished ^^^^^^^^^^^^^^^^^^
      return
      end
