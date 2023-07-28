      subroutine reset_temperature
	use params
      use lengths
      use RES
      use distres
      implicit integer(i-z)
!      parameter(ndim=1500)      !number of residues
!      parameter(nrep=100)       !number of replicas
!      parameter(nvec=416)       !number of vectors
!      common/lengths/Lch,Lch1,Lch2
      common/resnumber/Ncom,Ndis,accur
      common/commonuse2/atemp1,atemp2,N_rep,phot
      COMMON/RCN1/er1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
!      common/distres/er4,es3c
!      COMMON/RES/ER3,er5,er6,er7,Mcom(ndim),Kcom(ndim,100)
      real r_dis,r_con,r_dev,T10,T20,T1a,T2a

      if(er1+er3+er4.lt.0.1)return !without restrains
      
      a_rest=float(Ncom)/(1.3*Lch)
      a_rest=sqrt(a_rest)
      if(a_rest.lt.0.875)a_rest=0.875 !because of 80/70 -> rest/without_rest
      atemp2=atemp2*a_rest
      if(atemp2.gt.115) atemp2=115.0

      return
      end
