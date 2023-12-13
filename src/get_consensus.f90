      subroutine get_consensus
	use params
      use lengths
      implicit integer(i-z)
!      parameter(ndim=1500)
!      parameter(nrep=10)
!      parameter(nvec=416)
!      common/lengths/Lch,Lch1,Lch2
      common/consensus1/cx0(0:ndim),cy0(0:ndim),cz0(0:ndim),q(0:ndim)
      dimension cx1(ndim),cy1(ndim),cz1(ndim),q1(ndim),ip1(ndim)
      dimension cx2(ndim),cy2(ndim),cz2(ndim),q2(ndim),ip2(ndim)
      dimension qq(ndim)
      real xx,yy,zz
ccc   RMSD:
      double precision r_1(3,ndim),r_2(3,ndim),r_3(3,ndim),w(ndim)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /ndim*1.0/
ccc   

******************************************************************
***   read templates------------->
      rewind(24)
      read(24,*)N_tmp
      if(N_tmp.lt.2)then
         write(20,*)'There is only one template'
         write(20,*)'Top-1 template is used'
         return                 !use the first template
      else
***   read first template:
         read(24,*)N_al
         do i=1,N_al
            read(24,1237)text,ii,text,a1,a2,a3
            cx1(ii)=a1
            cy1(ii)=a2
            cz1(ii)=a3
            q1(ii)=1
         enddo
         read(24,*)text
***   read second template:
         read(24,*)N_al
         do i=1,N_al
            read(24,1237)text,ii,text,a1,a2,a3
            cx2(ii)=a1
            cy2(ii)=a2
            cz2(ii)=a3
            q2(ii)=1
         enddo
         read(24,*)text
***
      endif
 1237 format(A22,I4,A4,3F8.3)

******************************************************************
***   decided qq(i)------------->
      do i=1,Lch
         qq(i)=0
      enddo
      k=0
      do i=1,Lch
         if(q1(i).eq.1.and.q2(i).eq.1)then
            k=k+1
            ip1(k)=i
            r_1(1,k)=cx1(i)
            r_1(2,k)=cy1(i)
            r_1(3,k)=cz1(i)
            r_2(1,k)=cx2(i)
            r_2(2,k)=cy2(i)
            r_2(3,k)=cz2(i)
         endif
      enddo
      write(*,*)"#common aligned=",k
      if(k.lt.10)then
         write(20,*)'There is less than 10 common aligned residues'
         write(20,*)'Top-1 template is used'
         return                 !no common aligned points, using template1
      endif
      call u3b(w,r_1,r_2,k,1,rms,u,t,ier) !u rotate r_1 to r_2
      armsd=dsqrt(rms/k)        !RMSD12
      write(20,*)'RMSD1=',armsd,k !RMSD between template1,2 for common align
      kk=0
      do j=1,k
         xx=t(1)+u(1,1)*r_1(1,j)+u(1,2)*r_1(2,j)+u(1,3)*r_1(3,j)
         yy=t(2)+u(2,1)*r_1(1,j)+u(2,2)*r_1(2,j)+u(2,3)*r_1(3,j)
         zz=t(3)+u(3,1)*r_1(1,j)+u(3,2)*r_1(2,j)+u(3,3)*r_1(3,j)
         dis=sqrt((xx-r_2(1,j))**2+(yy-r_2(2,j))**2+(zz-r_2(3,j))**2)
c         write(*,*)j,ip1(j),dis
         if(dis.lt.5)then
            kk=kk+1
            qq(ip1(j))=1
         endif
      enddo
      if(kk.lt.10)then
         write(20,*)'There is less than 10 close common residues'
         write(20,*)'Top-1 template is used'
         return                 !no consensus points, using template1
      endif

******************************************************************
***   q(i)=qq(i), cx0(i)=(cx1+cx2)/2--------------->
      do i=1,Lch
         q(i)=qq(i)
         cx0(i)=1000000.        !for checking excluded volumn
         cy0(i)=1000000.
         cz0(i)=1000000.
      enddo
      k=0
      do i=1,Lch
         if(q(i).eq.1)then
            k=k+1
            ip2(k)=i
            r_1(1,k)=cx1(i)
            r_1(2,k)=cy1(i)
            r_1(3,k)=cz1(i)
            r_2(1,k)=cx2(i)
            r_2(2,k)=cy2(i)
            r_2(3,k)=cz2(i)
         endif
      enddo
      call u3b(w,r_1,r_2,k,1,rms,u,t,ier) !u rotate r_1 to r_2
      do j=1,k
         xx=t(1)+u(1,1)*r_1(1,j)+u(1,2)*r_1(2,j)+u(1,3)*r_1(3,j)
         yy=t(2)+u(2,1)*r_1(1,j)+u(2,2)*r_1(2,j)+u(2,3)*r_1(3,j)
         zz=t(3)+u(3,1)*r_1(1,j)+u(3,2)*r_1(2,j)+u(3,3)*r_1(3,j)
         cx0(ip2(j))=(xx+r_2(1,j))/2
         cy0(ip2(j))=(yy+r_2(2,j))/2
         cz0(ip2(j))=(zz+r_2(3,j))/2
      enddo
********************************************************************
      L_cut=2                   !small segment to be shrowed away.

c^^^^^^^^^^q(i) and cx0(i) obtained ^^^^^^^^^^^^^^^^^^^^

      return
      end
