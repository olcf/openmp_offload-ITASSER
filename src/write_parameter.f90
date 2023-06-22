      subroutine write_parameter
	use params
      use lengths
      use distres
      use pair1
      use order
      use RCN
      use RES
      implicit integer(i-z)
!      parameter(ndim=1500)
!      parameter(nvec=416)
!      parameter(nrep=100)       !number of replicas
      common/commonuse1/random,ncycle
      common/commonuse2/atemp1,atemp2,N_rep,phot
!      common/lengths/Lch,Lch1,Lch2
!      COMMON/RES/ER3,er5,er6,er7,Mcom(ndim),Kcom(ndim,100)
!      COMMON/RCN/Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
      COMMON/RCN1/er1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
!      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
!      common/distres/er4,es3c
      common/resnumber/Ncom,Ndis,accur
      common/initialinput/switch,k_cycle,k_phot,N_ann
      common/looks/exc,exc1,exc2
!      common/pair1/eh2,eh1b,eh1c
      common/ranzy/nozy
      common/iter/n_run

      write(20,*)'Protein real length..................',Lch
      write(20,*)
      write(20,*)'Ncycle...............................',ncycle
      write(20,*)
      write(20,*)'Number of runs.......................',n_run
      write(20,*)
      write(20,*)'number of N_rep......................',N_REP
      write(20,*)'N_swap...............................',phot
      write(20,*)'ncycle*N_swap........................',ncycle*phot
      write(20,*)'local moves each replica:',ncycle*phot*float(Lch)
      write(20,*)'total moves in MC..',ncycle*phot*float(Lch)*N_rep
      write(20,*) 
      write(20,*)'.......N_dist........................',Ndis
      write(20,*)'.......N_comb........................',Ncom
      write(20,*)'.......N_dist/length........',Ndis/float(Lch)
      write(20,*)'.......N_comb/length........',Ncom/float(Lch)
      write(20,*) 
      write(20,*)'maximum temperture...................',atemp2
      write(20,*)'minimum temperture...................',atemp1
      write(20,*)
      write(20,*)'Excluded volumn parameter=',exc,sqrt(exc)*.87
      write(20,*)
      write(20,*)'................er1........',er1
      write(20,*)'................er3........',er3
      write(20,*)'................er4........',er4
      write(20,*)'................er5........',er5
      write(20,*)'................er6........',er6
      write(20,*)'................er7........',er7
      write(20,*)'................eh1c.......',eh1c
      write(20,*)
      write(20,*)'initial random number seed...........',random
      write(20,*)'contact order........................',acorder
      write(20,*)
      write(20,*)'the first arandom number.............',aranzy(nozy)
      write(20,*) 
      write(20,*)'*******************************************'
      if(switch.gt.1)then
         write(20,*)'This is RS running of threading structure!'
      else
         write(20,*)'This is a normal simulation.'
      endif
      write(20,*)'*******************************************'
      write(20,*) 

c ^^^^^^^^^^ write parameter finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
