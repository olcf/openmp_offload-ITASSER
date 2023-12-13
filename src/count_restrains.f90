      subroutine count_restrains
	use params
      use lengths
      use RES
      use concutt
      use sg
      use seqe
      use longdist
      use RCN
      use CAcontact
      implicit integer (i-z)
!!      parameter(nrep=100)
!!      parameter(ndim=1500)	!maximum length of chain-length
!!      parameter(nvec=416)

!      common/lengths/Lch,Lch1,Lch2
      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
!      common/sg/gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
!      common/seqe/seq(ndim),sec(ndim)
!      common/concutt/concut(0:19,0:19),concut2(0:19,0:19),concut_sc
!      COMMON/RCN/Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
!      COMMON/RES/ER3,er5,er6,er7,Mcom(ndim),Kcom(ndim,100)
!      common/CAcontact/McomCA(ndim),KcomCA(ndim,100),aweigCA(ndim,ndim)
      common/CAcontact1/dist_CA_cut
!      COMMON/longdist/MdisL(ndim),kdisL(ndim,500),distL(ndim,500)
!      common/CA8/McomCA8(ndim),KcomCA8(ndim,100),aweigCA8(ndim,ndim)
      common/CAcontact2/dist_CA_cut8
 
      dimension ax(ndim),ay(ndim),az(ndim) !SG
      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim) !CA
      
      common/countres/t_comb,t_dist,t_combCA,s_comb,s_dist,s_combCA
      common/countres1/t_combCA8,t_distL,s_combCA8,s_distL,N_resc
      
      do i=1,3
         N_resc=N_resc+1
         do j=1,Lch
            x(j)=xrep(j,i)
            y(j)=yrep(j,i)
            z(j)=zrep(j,i)
         enddo
         do j=1,Lch1
            wx=x(j+1)-x(j)
            wy=y(j+1)-y(j)
            wz=z(j+1)-z(j)
            ica(j)=vector(wx,wy,wz) !identify order number of each bond-vector
         enddo
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         do j=1,Lch
            ax(j)=x(j)+gx(ica(j-1),ica(j),seq(j))
            ay(j)=y(j)+gy(ica(j-1),ica(j),seq(j))
            az(j)=z(j)+gz(ica(j-1),ica(j),seq(j))
         enddo

         do j=1,Lch
c     comb.dat--------->
            do k=1,Mcom(j)
               m=Kcom(j,k)      !k'th contact with j
               if(m.gt.j)then
                  t_comb=t_comb+1
                  dis=(ax(j)-ax(m))**2+(ay(j)-ay(m))**2
     $                 +(az(j)-az(m))**2
                  if(dis.le.concut2(seq(j),seq(m)))then
                     s_comb=s_comb+1
                  endif
               endif
            enddo
c     dist.dat--------->
            do k=1,Mdis(j)
               m=kdis(j,k)      !j-m restraints
               if(m.gt.j) then  !to avoid repeat
                  t_dist=t_dist+1
                  dis2=(x(j)-x(m))**2+(y(j)-y(m))**2+(z(j)-z(m))**2
                  dis=sqrt(dis2)
                  err=abs(dis-dist(j,k)) !dist: predicted dis
                  if(err.lt.dev(j,k)) then !dev: deviation for arca
                     s_dist=s_dist+1
                  endif
               endif
            enddo
c     distL.dat--------->
            do k=1,MdisL(j)
               m=kdisL(j,k)     !j-m restraints
               if(m.gt.j) then  !to avoid repeat
                  t_distL=t_distL+1
                  dis2=(x(j)-x(m))**2+(y(j)-y(m))**2+(z(j)-z(m))**2
                  dis=sqrt(dis2)
                  err=abs(dis-distL(j,k)) !dist: predicted dis
                  if(err.lt.1.5) then !dev: deviation for arca
                     s_distL=s_distL+1
                  endif
               endif
            enddo
c     combCA.dat--------->
            do k=1,McomCA(j)
               m=KcomCA(j,k)      !k'th contact with j
               if(m.gt.j)then
                  t_combCA=t_combCA+1
                  dis=(x(j)-x(m))**2+(y(j)-y(m))**2+(z(j)-z(m))**2
                  if(dis.le.dist_CA_cut)then
                     s_combCA=s_combCA+1
                  endif
               endif
            enddo
c     combCA8.dat--------->
            do k=1,McomCA8(j)
               m=KcomCA8(j,k)   !k'th contact with j
               if(m.gt.j)then
                  t_combCA8=t_combCA8+1
                  dis=(x(j)-x(m))**2+(y(j)-y(m))**2+(z(j)-z(m))**2
                  if(dis.le.dist_CA_cut8)then
                     s_combCA8=s_combCA8+1
                  endif
               endif
            enddo
ccc   
         enddo
      enddo

c^^^^^^^^^^^^^^^^^ restraints count finished ^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
