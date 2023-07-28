      function energy_tot()
	use params
      use lengths
      use backup2
      use ENERGY
      use envir1
      use ehbenergy
      use ehbenergy1
      use RES
      use pair1
      use shortcom
      use seqe
      use order
      use one
      use short
      use ehbc
      use backup1
      use shape
      
      implicit integer(i-z)
!      parameter(ndim=1500)
!      parameter(nvec=416)
!      common/lengths/Lch,Lch1,Lch2
!      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev    
!      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
!      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
!      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
!      common/backup2/eprofo,eprofn,energ
!      common/seqe/seq(ndim),sec(ndim)
!      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
!      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1
!      COMMON/ENERGY/EH5,ES3,ES3a,ES3b,EH1,EH1a,EH4,EHBIJ(ndim,ndim)
!      common/pair1/eh2,eh1b,eh1c
!      common/shortcom/eh3,es4,es5,es6,es7,es7a,es7b,es7c
      COMMON/RCN1/er1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
!      COMMON/RES/ER3,er5,er6,er7,Mcom(ndim),Kcom(ndim,100)
      common/otherenergy/E_cord,E_cnum

!      common/ehbenergy/EHB1,EHB1a,EHB1b,EHB1c,EHB2,EHB3,EHB4,EHB5,EHB6
!      common/ehbenergy1/EHB5a,EHB5b

      ICNTO=0
      SUMCTO=0
      DO i=1,Lch
         NOP(i)=0
         NOA(i)=0
         NOM(i)=0
      enddo

      energy_tot=EHB(1,Lch,1)+ESHORT(1,Lch,10)

      eprofo=0.0
!$acc data copy(ia,im,ip,is,eprofo,en2,E_cord,en3,E_cnum)
!$acc kernels
      DO k=1,Lch
         is=seq(k)
         ia=NOA(k)
         ip=NOP(k)
         im=NOM(k)
         eprofo=eprofo+envir(ia,im,ip,is,3)
      enddo
!$acc end kernels 

      E_cord=abs(ICNT/(0.00001+float(SUMCT))-acorder)
      E_cnum=abs(float(SUMCT)-contt)

      energy_tot=energy_tot+en1*eprofo+en2*E_cord+en3*E_cnum

c     following is for next movement ---------->
      icnto=icnt                !backup of contact-order
      sumcto=sumct              !backup of contact-order
      do i=1,Lch
         nopp(i)=nop(i)         !backup of number of contact-pair
         noaa(i)=noa(i)         !backup of number of contact-pair
         nomm(i)=nom(i)         !backup of number of contact-pair
      enddo
      codevsum=conew            !backup of panelity for total deviation of comb
      didevsum=dinew            !backup of panelity for total deviation of dist

c^^^^^^^^^^^^^^^^^^^^ Energy_tot finished ^^^^^^^^^^^^^^^^^^^^^^
!$acc end data       
      return
      end
