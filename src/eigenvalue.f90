      subroutine eigenvalue
	use params
        use eigen
c      common/eigen/AA(3,3),EE(3),HH(3,3)

c     #### EEigenvalues ############
      pi=3.1415926
      p1=-1
      p2=AA(1,1)+AA(2,2)+AA(3,3) 
      p3=-AA(1,1)*AA(2,2)-AA(1,1)*AA(3,3)-AA(2,2)*AA(3,3)+
     &     AA(2,3)*AA(2,3)+AA(1,3)*AA(1,3)+AA(1,2)*AA(1,2)
      p4=AA(1,1)*AA(2,2)*AA(3,3)+AA(1,2)*AA(2,3)*AA(3,1)+
     &     AA(1,3)*AA(2,1)*AA(3,2)-AA(1,3)*AA(2,2)*AA(3,1)-
     &     AA(1,1)*AA(2,3)*AA(3,2)-AA(1,2)*AA(2,1)*AA(3,3)
      p5=(-(1.0/3)*(p2/p1)**2+p3/p1)/3 
      ap5=sqrt(-p5)
      p6=((2.0/27)*(p2/p1)**3-(1.0/3)*(p2*p3/p1**2)+p4/p1)/2
      p7=acos(-p6/sqrt(-p5**3))
      p8=2*ap5*cos(p7/3.0)
      p9=-2*ap5*cos((p7+pi)/3.0)
      p10=-2*ap5*cos((p7-pi)/3.0) 
      p11=p2/(3*p1)
      EE(1)=p8-p11               !eigenvalue
      EE(2)=p9-p11               !eigenvalue
      EE(3)=p10-p11              !eigenvalue
      
c     ##### normalized eigenvectors #########
      do i=1,3
         fnorm1=AA(2,1)*AA(1,2)-(AA(1,1)-EE(i))*(AA(2,2)-EE(i))
         x=((AA(2,2)-EE(i))*AA(1,3)-AA(1,2)*AA(2,3))/fnorm1
         y=((AA(1,1)-EE(i))*AA(2,3)-AA(2,1)*AA(1,3))/fnorm1
         HH(i,3)=1/sqrt(x*x+y*y+1)
         HH(i,1)=x*HH(i,3)
         HH(i,2)=y*HH(i,3)
c     write(*,*)i,EE(i),HH(i,1),HH(i,2),HH(i,3)
      enddo

      return
      end
