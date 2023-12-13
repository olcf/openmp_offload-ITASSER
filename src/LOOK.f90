      FUNCTION LOOK(jj,kk)
	use params
      use chainm
      use chain1
      use seqe
      use echain1
      use lengths
      IMPLICIT INTEGER(I-Z)
      LOGICAL LOOK
!      parameter(ndim=1500)
!      parameter(nvec=416)
!      common/seqe/seq(ndim),sec(ndim) 	 	
!      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
!      common/lengths/Lch,Lch1,Lch2
      common/looks/exc,exc1,exc2
!      common/chainm/mv(ndim)
!      common/echain1/ex(ndim),ey(ndim),ez(ndim)
      common/initialinput/switch,k_cycle,k_phot,N_ann
      common/excluded/vvv(ndim,ndim)

c     Ca-Ca r2>13;   r>sqrt(13)*0.87=3.14A
c     Ca-Ca r2>14;   r>sqrt(14)*0.87=3.26A
c     Ca-Ca r2>15;   r>sqrt(15)*0.87=3.37A
c     Ca-Ca r2>16;   r>sqrt(16)*0.87=3.48A
c     Ca-Ca r2>17;   r>sqrt(17)*0.87=3.59A
c     Ca-Ca r2>18;   r>sqrt(18)*0.87=3.69A
c     Ca-Ca r2>19;   r>sqrt(19)*0.87=3.79A
c     Ca-Ca r2>20;   r>sqrt(20)*0.87=3.89A
c     Ca-Ca r2>21;   r>sqrt(21)*0.87=3.99A
c     Ca-Ca r2>22;   r>sqrt(22)*0.87=4.08A
c     Ca-Ca r2>23;   r>sqrt(23)*0.87=4.17A *
c     Ca-Ca r2>24;   r>sqrt(24)*0.87=4.26A
c     Ca-Ca r2>25;   r>sqrt(25)*0.87=4.35A
c     Ca-Ca r2>26;   r>sqrt(26)*0.87=4.44A
c     Ca-Ca r2>27;   r>sqrt(27)*0.87=4.52A

      if(switch.eq.1)then       !normal lattice running
         do k=jj,kk
            i1=k-3
            i2=max(k+3,kk+1)
            do i=1,Lch
               if(i.le.i1.or.i.ge.i2)then
                  ir2=(x(k)-x(i))**2+(y(k)-y(i))**2+(z(k)-z(i))**2
                  if(ir2.lt.exc) then
                     LOOK=.FALSE.
                     RETURN
                  endif
               endif
            enddo
         enddo
      else                      !off-lattice
         do k=jj,kk
            if(mv(k).gt.0)then  !off-lattice
               axk=x(k)
               ayk=y(k)
               azk=z(k)
            else                !on-lattice
               axk=ex(k)
               ayk=ey(k)
               azk=ez(k)
            endif
            i1=k-3
            i2=max(k+3,kk+1)
            do i=1,Lch
               if(i.le.i1.or.i.ge.i2)then
                  if(vvv(i,k).gt.0)then
                     if(mv(i).gt.0)then
                        axi=x(i)
                        ayi=y(i)
                        azi=z(i)
                     else
                        axi=ex(i)
                        ayi=ey(i)
                        azi=ez(i)
                     endif
                     ar2=(axk-axi)**2+(ayk-ayi)**2+(azk-azi)**2
                     if(ar2.lt.exc) then
                        LOOK=.FALSE.
                        RETURN
                     endif
                  endif
               endif
            enddo
         enddo
      endif
      LOOK=.TRUE.

c ^^^^^^^^^^ Look finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      RETURN
      END
