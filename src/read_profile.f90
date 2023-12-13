      subroutine read_profile
	use params
      use one
      use ehbc
      implicit integer(i-z)
!      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
!      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1

      do kkk=1,4
      do i=0,19
      do im=0,15
      do ip=0,15
      do ia=0,15
         envir(im,ip,ia,i,kkk)=2.0
      end do
      end do
      end do
      end do
      enddo

c     PROFILE3 potential =envir(#of antiparalel contacts,
c     #of orthogonal, # of parallel, aminoacid's type)
c     ia,im, ip - taken modulo 2, i.e. 0-1 contact, 2-3,...
c     profile3.comm is a score table, i.e. envir(envir_class,A)
c     here environment class is number of contacts on residue A.

      do i=0,19                 !from column
         read(2,*)
         do im=0,8
         do ia=0,8
            read(2,*)(envir(ia,im,ip,i,3),ip=0,8) !this is used
         end do
         read(2,*)
         end do
         read(2,*)
      enddo

c      do i=0,19                 !from column
c         read(2,*)
c         do im=0,8
c         do ia=0,8
c            read(2,*)(envir(ia,im,ip,i,1),ip=0,8)
c         end do
c         read(2,*)
c         end do
c         read(2,*)
c      end do

c      do i=0,19                 !from column
c         read(2,*)
c         do im=0,8
c         do ia=0,8
c            read(2,*)(envir(ia,im,ip,i,2),ip=0,8)
c         end do
c         read(2,*)
c         end do
c         read(2,*)
c      end do

c      do i=0,19                 !from column
c         read(2,*)
c         do im=0,8
c         do ia=0,8
c            read(2,*)(envir(ia,im,ip,i,4),ip=0,8)
c         end do
c         read(2,*)
c         end do
c         read(2,*)
c      end do

c^^^^^^^^^^^^^^^^^^^^^^^^ read profile finished ^^^^^^^^^^^^^^^^^^^^
      return
      end
