      subroutine prepare_move2
	use params
      implicit integer(i-z)
!      parameter(nvec=416)
      common/logica/goodc
      logical goodc(nvec,nvec)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)
      common/w2a/w21(-10:10,-10:10,-10:10,80)
      common/w2b/w22(-10:10,-10:10,-10:10,80)
      common/w2c/Nw(-10:10,-10:10,-10:10)
      common/nwmax1/nvecr
      common/lattice/m_latt,latt1,latt2

      do i=-10,10
         do j=-10,10
            do k=-10,10
               Nw(i,j,k)=0
            enddo
         enddo
      enddo

      max=0
      maxa=0
      nnn=0
      mmm=0

      do 101 i=1,nvecr
         do 102 j=1,nvecr
            if(goodc(i,j))then
               ijx=vx(i)+vx(j)  !vx in (-5,5)
               ijy=vy(i)+vy(j)
               ijz=vz(i)+vz(j)
***   based on ijx:
               Nw(ijx,ijy,ijz)=Nw(ijx,ijy,ijz)+1 !based on r
               w21(ijx,ijy,ijz,Nw(ijx,ijy,ijz))=i
               w22(ijx,ijy,ijz,Nw(ijx,ijy,ijz))=j
               if(maxa.lt.Nw(ijx,ijy,ijz))maxa=Nw(ijx,ijy,ijz)
***   based on i,j:
               Np2(i,j)=0
               do 103 ii=1,nvecr
                  jx=ijx-vx(ii)
                  jy=ijy-vy(ii)
                  jz=ijz-vz(ii)
                  jr=jx*jx+jy*jy+jz*jz
                  if(jr.ge.latt1.and.jr.le.latt2)then
                     jj=vector(jx,jy,jz)
                     if(goodc(ii,jj))then
                        Np2(i,j)=Np2(i,j)+1 !based on i,j
                        v21(i,j,Np2(i,j))=ii
                        v22(i,j,Np2(i,j))=jj
                        mmm=mmm+1 !number of total memory occupied by vx21.
                     endif
                  endif
 103           continue
               if(max.lt.Np2(i,j))max=Np2(i,j)
c               write(*,*)i,j,Np2(i,j),max,maxa
               nnn=nnn+1        !number of possible pairs
            endif
 102     continue
 101  continue

ccc among all 312*312=97344 pairs, 67272 pairs are legal (~70%).
ccc all Np2(i,j)>=2, i.e. there are at least one other path for any pair.
ccc <Np2(i,j)>=27.
      
      write(20,*)'maximum of Np2(i,j)=',max,maxa
      write(20,*)'number of possible pair=',nnn
      write(20,*)'sum of Np2(i,j), total memory=',mmm
cccc  the following is the cases without goodc on prefabracated pairs:
c              N_v  N_pare N_pp   mmm
c     [14,25], 312, 97344, 67272, 1830000 (25% memory, within memory)

cccc the following is cases without goodc limitation:
c     [14,25], 306, 93636, 3872214 (90% memory, beyond memory)
c     [14,24], 282, 79524, 2923758 (80% memory, beyond memory)
c      stop

      return
      end
