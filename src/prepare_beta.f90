      subroutine prepare_beta
	use params
      use hb
      use bisec
      use sg
      implicit integer(i-z)
!                parameter(nvec=416)
      common/logica/goodc
      logical goodc(nvec,nvec)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
!      common/bisec/cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
!      common/hb/hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
      common/three/angle(nvec,nvec)
!      common/sg/gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
!      common/cb/hx(nvec,nvec,0:19),hy(nvec,nvec,0:19),hz(nvec,nvec,0:19)
      common/beta1/axalf(0:19),ayalf(0:19),azalf(0:19)
      common/beta2/axbet(0:19),aybet(0:19),azbet(0:19)
      common/beta3/bxalf(0:19),byalf(0:19),bzalf(0:19)
      common/beta4/bxbet(0:19),bybet(0:19),bzbet(0:19)
      common/nwmax1/nvecr


      do k=0,19
         read(4,*)
         read(4,*)axalf(k),ayalf(k),azalf(k),axbet(k),aybet(k),azbet(k)
      enddo
      do k=0,19
         read(4,*)
         read(4,*)bxalf(k),byalf(k),bzalf(k),bxbet(k),bybet(k),bzbet(k)
      enddo
      CLOSE(4)                  !sidecent.comm

ccccccccccc define hydrogen-bond, C_beta, C_group for good (i,j)----->
      esp=0.000001
      do 101 i=1,nvecr
      do 102 j=1,nvecr
         avi=sqrt(float(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)))
         avxi=vx(i)/avi
         avyi=vy(i)/avi
         avzi=vz(i)/avi
         avj=sqrt(float(vx(j)*vx(j)+vy(j)*vy(j)+vz(j)*vz(j)))
         avxj=vx(j)/avj
         avyj=vy(j)/avj
         avzj=vz(j)/avj
         
ccc   if vi and vj is parallel, a is ok but b=c=0 ------->
         if(abs(avxi-avxj).lt.esp)then
            if(abs(avyi-avyj).lt.esp)then
               if(abs(avzi-avzj).lt.esp)then
                  ax=avxi+avxj
                  ay=avyi+avyj
                  az=avzi+avzj
                  aaa=sqrt(ax*ax+ay*ay+az*az)
                  ax=ax/aaa     !(vi+vj)/|vi+vj|
                  ay=ay/aaa     !(vi+vj)/|vi+vj|
                  az=az/aaa     !(vi+vj)/|vi+vj|
                  
c     calculate b=a(x)u, u=(1,1,1):
                  bx=ay*1-az*1
                  by=az*1-ax*1
                  bz=ax*1-ay*1
                  bbb=sqrt(bx*bx+by*by+bz*bz)
                  bx=bx/bbb     ! a(x)1/|a(x)1|
                  by=by/bbb     ! a(x)1/|a(x)1|
                  bz=bz/bbb     ! a(x)1/|a(x)1|
                  
c     calculate c=a(x)b:
                  cx=ay*bz-az*by
                  cy=az*bx-ax*bz
                  cz=ax*by-ay*bx
                  ccc=sqrt(cx*cx+cy*cy+cz*cz)
                  cx=cx/ccc     ! a(x)b/|a(x)b|
                  cy=cy/ccc     ! a(x)b/|a(x)b|
                  cz=cz/ccc     ! a(x)b/|a(x)b|
                  cax(i,j)=cx
                  cay(i,j)=cy
                  caz(i,j)=cz
                  goto 29
               endif
            endif
         endif
ccc   check if vi and vj is anti-parallel, c is ok, b=a=0 ------->
         if(abs(avxi+avxj).lt.esp)then
            if(abs(avyi+avyj).lt.esp)then
               if(abs(avzi+avzj).lt.esp)then
                  cx=avxi-avxj
                  cy=avyi-avyj
                  cz=avzi-avzj
                  ccc=sqrt(cx*cx+cy*cy+cz*cz)
                  cx=cx/ccc     !(vi-vj)/|vi-vj|
                  cy=cy/ccc     !(vi-vj)/|vi-vj|
                  cz=cz/ccc     !(vi-vj)/|vi-vj|
                  cax(i,j)=cx   !(vi-vj)/|vi-vj|
                  cay(i,j)=cy   !(vi-vj)/|vi-vj|
                  caz(i,j)=cz   !(vi-vj)/|vi-vj|

c     calculate a=c(x)u, u=(1,1,1):
                  ax=cy*1-cz*1
                  ay=cz*1-cx*1
                  az=cx*1-cy*1
                  aaa=sqrt(ax*ax+ay*ay+az*az)
                  ax=ax/aaa     ! c(x)1/|c(x)1|
                  ay=ay/aaa     ! c(x)1/|c(x)1|
                  az=az/aaa     ! c(x)1/|c(x)1|
                  
c     calculate b=c(x)a:
                  bx=cy*az-cz*ay
                  by=cz*ax-cx*az
                  bz=cx*ay-cy*ax
                  bbb=sqrt(bx*bx+by*by+bz*bz)
                  bx=bx/bbb     ! c(x)a/|c(x)a|
                  by=by/bbb     ! c(x)a/|c(x)a|
                  bz=bz/bbb     ! c(x)a/|c(x)a|
                  goto 29
               endif
            endif
         endif
         
         ax=avxi+avxj
         ay=avyi+avyj
         az=avzi+avzj
         aaa=sqrt(ax*ax+ay*ay+az*az)
         ax=ax/aaa              !(vi+vj)/|vi+vj|
         ay=ay/aaa              !(vi+vj)/|vi+vj|
         az=az/aaa              !(vi+vj)/|vi+vj|
         
         bx=avyi*avzj-avzi*avyj
         by=avzi*avxj-avxi*avzj
         bz=avxi*avyj-avyi*avxj
         bbb=sqrt(bx*bx+by*by+bz*bz)
         bx=bx/bbb             ! vi(x)vj/|vi(x)vj|
         by=by/bbb             ! vi(x)vj/|vi(x)vj|
         bz=bz/bbb             ! vi(x)vj/|vi(x)vj|

         cx=avxi-avxj
         cy=avyi-avyj
         cz=avzi-avzj
         ccc=sqrt(cx*cx+cy*cy+cz*cz)
         cx=cx/ccc              !(vi-vj)/|vi-vj|
         cy=cy/ccc              !(vi-vj)/|vi-vj|
         cz=cz/ccc              !(vi-vj)/|vi-vj|
         cax(i,j)=cx            !(vi-vj)/|vi-vj|
         cay(i,j)=cy            !(vi-vj)/|vi-vj|
         caz(i,j)=cz            !(vi-vj)/|vi-vj|
         
 29      continue

         goto 39
c     check if aaa/bbb/ccc=0 ------->
         if(aaa.lt.esp.or.bbb.lt.esp.or.ccc.lt.esp)then
            write(*,19)i,j,avxi,avyi,avzi,aaa,bbb,ccc
         endif
c     check if a.b=b.c=a.c=0 ------->
         ab=ax*bx+ay*by+az*bz
         bc=bx*cx+by*cy+bz*cz
         ac=ax*cx+ay*cy+az*cz
         if(ab.gt.esp.or.ab.gt.esp.or.ac.gt.esp)then
            write(*,19)i,j,avxi,avyi,avzi,ab,bc,ac
         endif
c         write(*,19)i,j,avxi,avyi,avzi,ab,bc,ac
 19      format(2i5,10f8.3)
 39       continue

c     H-bond (unit vector):
         hbx(i,j)=bx
         hby(i,j)=by
         hbz(i,j)=bz

c     side-chain coordinate from C_a to side-chain ---------------->
         do k=0,19
            if(angle(i,j).lt.105) then ! alpha-helix or turn like
               gx(i,j,k)=(axalf(k)*ax+ayalf(k)*bx+azalf(k)*cx)/0.87
               gy(i,j,k)=(axalf(k)*ay+ayalf(k)*by+azalf(k)*cy)/0.87
               gz(i,j,k)=(axalf(k)*az+ayalf(k)*bz+azalf(k)*cz)/0.87

               hx(i,j,k)=(bxalf(k)*ax+byalf(k)*bx+bzalf(k)*cx)/0.87
               hy(i,j,k)=(bxalf(k)*ay+byalf(k)*by+bzalf(k)*cy)/0.87
               hz(i,j,k)=(bxalf(k)*az+byalf(k)*bz+bzalf(k)*cz)/0.87
            else                ! beta-sheet
               gx(i,j,k)=(axbet(k)*ax+aybet(k)*bx+azbet(k)*cx)/0.87
               gy(i,j,k)=(axbet(k)*ay+aybet(k)*by+azbet(k)*cy)/0.87
               gz(i,j,k)=(axbet(k)*az+aybet(k)*bz+azbet(k)*cz)/0.87

               hx(i,j,k)=(bxbet(k)*ax+bybet(k)*bx+bzbet(k)*cx)/0.87
               hy(i,j,k)=(bxbet(k)*ay+bybet(k)*by+bzbet(k)*cy)/0.87
               hz(i,j,k)=(bxbet(k)*az+bybet(k)*bz+bzbet(k)*cz)/0.87
            endif
         enddo
 102  continue
 101  continue

!$acc update device(hbx,hby,hbz,gx,gy,gz)

c      stop
      return
      end
