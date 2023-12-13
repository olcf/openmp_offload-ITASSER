      subroutine get_acorder
	use params
      use lengths
      use seqe
      use order
      implicit integer(i-z)
!      parameter(ndim=1500)
!      parameter(nvec=416)
!      common/seqe/seq(ndim),sec(ndim)
!      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
!      common/lengths/Lch,Lch1,Lch2
      
c     data struct /' coil','helix',' turn',' beta','    ?'/

************* number of predicted structures *************
      n_H=0                     !number of Helix
      n_E=0                     !number of Extension
      do i=1,Lch
         if(sec(i).eq.2) n_H=n_H+1
         if(sec(i).eq.4) n_E=n_E+1
      enddo

      if(n_H+n_E.lt.2)then      !use alpha=H/E
         if(Lch.lt.50)then
            alph=0.252
         else if(Lch.lt.70)then
            alph=0.230
         else if(Lch.lt.90)then
            alph=0.212
         else if(Lch.lt.110)then
            alph=0.198
         else if(Lch.lt.130)then
            alph=0.168
         else if(Lch.lt.150)then
            alph=0.176
         else if(Lch.lt.170)then
            alph=0.183
         else if(Lch.lt.190)then
            alph=0.160
         else
            alph=0.134
         endif
      else                      !use alpha=aH+bE
         a1=float(n_H)/float(n_H+n_E)
         a2=float(n_E)/float(n_H+n_E)
         if(Lch.lt.50)then
            alph=0.116*a1+0.324*a2
         else if(Lch.lt.70)then
            alph=0.119*a1+0.357*a2
         else if(Lch.lt.90)then
            alph=0.115*a1+0.280*a2
         else if(Lch.lt.110)then
            alph=0.105*a1+0.259*a2
         else if(Lch.lt.130)then
            alph=0.132*a1+0.269*a2
         else if(Lch.lt.150)then
            alph=0.105*a1+0.272*a2
         else if(Lch.lt.170)then
            alph=0.114*a1+0.186*a2
         else if(Lch.lt.190)then
            alph=0.116*a1+0.197*a2
         else
            alph=0.107*a1+0.184*a2
         endif
      endif

      acorder=alph*Lch
c^^^^^^^^^^^^^^^^^^ contact order done ^^^^^^^^^^^^^^^^^^
      return
      end
