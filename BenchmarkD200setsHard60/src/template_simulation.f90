      subroutine template_simulation
	use params
      use chainm
      use echain1
      use lengths
      use chain1
      use echain2
      use echain4
      use order
      use echain5
      use echain6
      use eigen
      use short
      use backup1
      use omp_lib
      use iso_fortran_env

      implicit integer(i-z)
!      use chain1
!      parameter(ndim=1500)
!      parameter(nrep=10)
!      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
!      common/lengths/Lch,Lch1,Lch2
!      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/xyzrs/xrs(ndim,40,40),yrs(ndim,40,40),zrs(ndim,40,40)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/commonuse1/random,ncycle
      common/commonuse2/atemp1,atemp2,N_rep,phot
!      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev    
!      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
!      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      dimension x0(ndim),y0(ndim),z0(ndim)  !for E_min
      character*6 mname
      common/movename/mname(100)
      common/nrepfile/n_repf

      common/temperature/itemp,atemp

      common/chain0/ras(ndim),nfl
      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/mng/m_g(100)

      common/initialinput/switch,k_cycle,k_phot,N_ann
      common/thrtem/aT_ann(100)
      common/eall/N_sum(100),energ_sum(100),energ_sum2(100),E_min
      common/acct/accept0
      common/moveretio2/bh2,bh3s,bh3d,bh4s,bh4d,bh5s,bh5d,bh6
      common/moveretio3/bh7a,bh7b,bhendn,bhendc

      common/nswap/bNSa(100),bNSt(100)
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t

      dimension aNNa(100),aNNt(100),N_out(100) !for output annealing
       real :: TimeB 
     
      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)

      common/sw1/aT_rep(nrep),E_rep(nrep)
      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      common/sw3/icarep(ndim,nrep)
      common/sw4/exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep)
      common/sw5/egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep)
      common/sw7/ecxrep(ndim,nrep),ecyrep(ndim,nrep),eczrep(ndim,nrep)
      common/sw8/ebxrep(ndim,nrep),ebyrep(ndim,nrep),ebzrep(ndim,nrep)
      common/sw9/etxrep(ndim,nrep),etyrep(ndim,nrep),etzrep(ndim,nrep)

      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
!      common/echain1/ex(ndim),ey(ndim),ez(ndim)
!      common/echain2/egx(ndim),egy(ndim),egz(ndim)
!      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim)
!      common/echain5/ebx(ndim),eby(ndim),ebz(ndim)
!      common/echain6/etx(ndim),ety(ndim),etz(ndim)
!      common/chainm/mv(ndim)
      common/outputxyz/fxyz(3,ndim)
      character*3 sequ
      common/aminoacid/sequ(ndim)
      character text
      common/E_defo/i_E_defo
      common/rmsd_defo/armsd,armsd_sum(100),N_rmsd(100)
      common/aTs/aTs1,aTs2,aTs_rep(nrep),aTs
      common/aTTs/aTTs1,aTTs2,aTTs_rep(nrep),aTTs
      common/ranzy/nozy
      common/hours/hour_max
      real (real64):: mvTs, mvTe, mvT,ehbTs,ehbTe,ehbT
! !$acc data copyout (energ_sum2(:),energ_sum(:),bnt(:)
! !$acc&              ,n_rmsd(:),bnna(:),bna(:),bnst(:)
! !$acc&               ,armsd_sum(:),bnsa(:),bnnt(:),n_sum(:))  
! !$acc&  copy(ezrep(:,:),
! !$acc&    etyrep(:,:),etxrep(:,:),egzrep(:,:),
! !$acc&    ebyrep(:,:),ebxrep(:,:),ecyrep(:,:),
! !$acc&    ecxrep(:,:),ebzrep(:,:),egyrep(:,:),egxrep(:,:),
! !$acc&    eczrep(:,:),eyrep(:,:),exrep(:,:),etzrep(:,:))
! !$acc kernels loop gang private(i)
      do i=1,100
         bNSa(i)=0              !aceptance for swep
         bNSt(i)=0
         bNa(i)=0               !acceptance for move2,3,4,5,6,7
         bNt(i)=0
         bNNa(i)=0              !acceptance for different temperature.
         bNNt(i)=0
         N_sum(i)=0
         energ_sum(i)=0         !<E_tot>
         energ_sum2(i)=0        !<E_tot^2>
         armsd_sum(i)=0
         N_rmsd(i)=0
      enddo
! !$acc end kernels
      E_min=10000

      i_tr=0                    !order number of output trajectory
      i_fra=0
      mcycle=0
!$acc update device(ica,AA)
      do 1111 icycle=1,ncycle
         do 2222 itemp=1,N_rep  !iterate for all the replicas
            atemp=aT_rep(itemp) !current temperature
            aTs=aTs_rep(itemp)
            aTTs=aTTs_rep(itemp)
            call set_current_RS !get current (x,y,z,ica,T)
            call initial_move   !update center, axis, energy
ccc   
            do 3333 iphot=1,phot !N_swap, iterate at fixed temperature
               do 4444 i_nfl=1,nfl
                  fff=aranzy(nozy)
                  if(fff.le.bh2)then
                     call move2
                  elseif(fff.le.bh3s)then
                     if(nfl3.ge.1)then
                        call move3s
                     else
                        call move2
                     endif
                  elseif(fff.le.bh3d)then
                     if(nfl3.ge.1)then
                        call move3d
                     else
                        call move2
                     endif
                  elseif(fff.le.bh4s)then
                     if(nfl4.ge.1)then
                        call move4s
                     else
                        call move2
                     endif
                  elseif(fff.le.bh4d)then
                     if(nfl4.ge.1)then
                        call move4d
                     else
                        call move2
                     endif
                  elseif(fff.le.bh5s)then
                     if(nfl5.ge.1)then
                        call move5s
                     else
                        call move2
                     endif
                  elseif(fff.le.bh5d)then
                     if(nfl5.ge.1)then
                        call move5d
                     else
                        call move2
                     endif
                  elseif(fff.le.bh6)then
                     if(nfl6.ge.1)then
                        call move6
                     else
                        call move2
                     endif
                  elseif(fff.le.bhendn)then
                     if(ras(1).eq.1)then
                        call move_n_end
                     else
                        call move2
                     endif
                  else
                     if(ras(nfl).eq.Lch)then
                        call move_c_end
                     else
                        call move2
                     endif
                  endif
ccccccccccccccccccccfragment movement ccccccccccccccccccccc
                  if(switch.eq.3)goto 4441 !freeze the template
                  i_fra=i_fra+1
                  if(i_fra.ge.n_fra.and.nfr.gt.0)then !1 chunk move in nfr move
                     i_fra=0
                     ifr=int(aranzy(nozy)*nfr)+1 ![1,nfr]
                     if(switch.eq.2)then !rotation+translation
                        if(nfr_i(ifr).eq.1)then
                           call trot_N !translate+rotate N-terminal
                        elseif(nfr_f(ifr).eq.Lch)then
                           call trot_C !translate+rotate C-terminal
                        else
                         mvTs=omp_get_wtime()     
                           call trot_M !translate+rotate Middle-fragment
                         mvTe=omp_get_wtime() 
                        endif
                     elseif(switch.eq.4)then !translation
                        if(nfr_i(ifr).eq.1)then
                           call tran_N !translate N-terminal
                        elseif(nfr_f(ifr).eq.Lch)then
                           call tran_C !translate C-terminal
                        else
                           call tran_M !translate Middle-fragment
                        endif
                     elseif(switch.eq.5)then !rotation+translation+deformation
                        aran_num=aranzy(nozy)
                        if(nfr_i(ifr).eq.1)then
                           if(aran_num.gt.0.5)then
                              i_E_defo=3
                              call trot_N !translate+rotate N-terminal
                           else
                              i_E_defo=1
                              call defo_N !translate+deform N-terminal
                           endif
                        elseif(nfr_f(ifr).eq.Lch)then
                           if(aran_num.gt.0.5)then
                              i_E_defo=3
                              call trot_C !translate+rotate C-terminal
                           else
                              i_E_defo=1
                              call defo_C !translate+deform C-terminal
                           endif
                        else
                           if(aran_num.gt.0.5)then
                              i_E_defo=3
                              call trot_M !translate+rotate Middle-fragment
                           else
                              i_E_defo=1
                              call defo_M !translate+deform Middle-fragment
                           endif
                        endif
                     endif
                  endif
 4441             continue

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           call CPU_TIME(TimeB)
                  atime=TimeB/3600.0               
                   !atime=second()/3600.0
                    if(atime.gt.hour_max)goto 901
 4444           continue
 3333         continue
ccc
ccccccccccrecord energy and (x,y,z) cccccccccccccccccccccccc
            E_rep(itemp)=energy_tot() !whole energy
            if(E_rep(itemp).lt.E_min) E_min=E_rep(itemp)
! !$acc kernels loop  private(i)\
! !$OMP target teams distribute parallel do
! !$OMP parallel do simd private(i) shared(ica,x,y,z,ex,ey,ez)            
            do i=1,Lch
               icarep(i,itemp)=ica(i)
               if(mv(i).gt.0)then
                  xrep(i,itemp)=x(i)
                  yrep(i,itemp)=y(i)
                  zrep(i,itemp)=z(i)
               else
                  exrep(i,itemp)=ex(i) !Ca
                  eyrep(i,itemp)=ey(i)
                  ezrep(i,itemp)=ez(i)
                  egxrep(i,itemp)=egx(i) !SG
                  egyrep(i,itemp)=egy(i)
                  egzrep(i,itemp)=egz(i)
                  ecxrep(i,itemp)=ecx(i) !cc
                  ecyrep(i,itemp)=ecy(i)
                  eczrep(i,itemp)=ecz(i)
                  ebxrep(i,itemp)=ebx(i) !Hb
                  ebyrep(i,itemp)=eby(i)
                  ebzrep(i,itemp)=ebz(i)
                  etxrep(i,itemp)=etx(i) !CB
                  etyrep(i,itemp)=ety(i)
                  etzrep(i,itemp)=etz(i)
               endif
            enddo
! !$OMP end parallel do             
!!$OMP end target teams distribute parallel do
! !$acc end kernels
c^^^^^^^^^^^^record done ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 2222    continue

ccccccccccccccccc<RMSD>, <E> cccccccccccccccccccccccccc
! !$acc kernels loop gang private(i)
!!$OMP target teams distribute parallel do       
!$OMP parallel do simd           
         do i=1,N_rep
            energ_sum(i)=energ_sum(i)+E_rep(i)
            energ_sum2(i)=energ_sum2(i)+E_rep(i)*E_rep(i)
            N_sum(i)=N_sum(i)+1
         enddo
!$OMP end parallel do simd         
!!$OMP end target teams distribute parallel do
! !$acc end kernels
!
ccccccccccccccccccccc snapshots of E(1), E(N_rep) ccccccccccccc
         if(icycle.eq.icycle/1*1)then
            i_tr=i_tr+1
            do k=1,n_repf
               write(30+k,401)Lch,E_rep(k),i_tr,icycle
               do i=1,Lch
                  if(mv(i).gt.0)then
                     abx=xrep(i,k)*0.87
                     aby=yrep(i,k)*0.87
                     abz=zrep(i,k)*0.87
                  else
                     abx=exrep(i,k)*0.87
                     aby=eyrep(i,k)*0.87
                     abz=ezrep(i,k)*0.87
                  endif
                  write(30+k,402)abx,aby,abz
               enddo
            enddo
         endif
 401     format(i8,1x,f10.1,2i8)
 402     format(f10.3,1x,f10.3,1x,f10.3)

ccccccccccccccccc swap replicas cccccccccccccccccccccccccccccccc
         if(icycle.eq.icycle/2*2)then
            do i=1,N_rep-1,2 !swap odd replicas
               call swap_RS(i,i+1)
            enddo
         else
            do i=2,N_rep-1,2
               call swap_RS(i,i+1) !swap even replicas
            enddo
          endif
          mcycle=mcycle+1
 1111 continue
 901  continue
! !$acc end data
c--------------------------Main cycle ended here!!---------------
      
cccccccccccccccccccccccc Na/Nt cccccccccccccccccccccccccc
      WRITE(20,*)
      WRITE(20,*) 'RS-->i_move   move   Na(i)  Nt(i)   Na(i)/Nt(i)'
      do i=2,27
         if(bNt(i).gt.1)then
            write(20,5004) i,mname(i),bNa(i),bNt(i),bNa(i)/bNt(i)
         else
            write(20,5004) i,mname(i),bNa(i),bNt(i)
         endif
      enddo
 5004 format(I4,A9,2f15.1,f11.6)
      
ccccccccccccccccccccccccccE_final, NSa/NSt ccccccccccccccccccccccc
      WRITE(20,*)
      WRITE(20,*)'RS-->------------ E_final, Na_swap/Nt_swap ---------'
      WRITE(20,*) 'i  T(i) final_E(i)  Nsa(i)  Nst(i)  Nsa(i)/Nst(i)'
      do i=1, N_rep
         if(bNSt(i).gt.1)then
            write(20,5005) i,aT_rep(i),E_rep(i),
     $           bNSa(i),bNSt(i),bNSa(i)/bNSt(i)
         else
            write(20,5005) i,aT_rep(i),E_rep(i),bNSa(i),bNSt(i)
         endif
      enddo
 5005 format(I4,f7.2,f10.1,2f15.1,f11.6)
      energy_tot_tmp=energy_tot()
      write(20,*)'E_final=',energy_tot_tmp

ccccccccccccccccccccccc<E>, NNa/NNt ccccccccccccccccccccccccccc
      WRITE(20,*)
      WRITE(20,*)'RS---------- <energy>, Na(i)/Nt(i) ----------------'
      write(20,*)'i_rep   T    <E>   <RMSD>   c   Na   Nt   Na/Nt'
      do i=1,N_rep
        energ_sum(i)=energ_sum(i)/N_sum(i)
        energ_sum2(i)=energ_sum2(i)/N_sum(i)
        cheat=(energ_sum2(i)-energ_sum(i)**2)/(aT_rep(i)**2)
        write(20,5006)i,aT_rep(i),energ_sum(i),armsd_sum(i)
     &       /(N_rmsd(i)+0.00001),
     &       cheat,bNNa(i),bNNt(i),bNNa(i)/(bNNt(i)+0.00001)
      enddo
 5006 format(I4,f7.2,f10.1,f8.3,f8.1,2f12.1,f11.6)
      write(20,*)'E_min=',E_min

      write(20,*)
      write(20,*)'ncycle_max=',ncycle
      write(20,*)'ncycle_real=',mcycle
      write(20,*)
      write(20,*)'hour_max=',hour_max
      write(20,*)'hour_real=',atime
      write(20,*)
!!FIX      write(20,*)'ending time: ',fdate()
      
      write(*,*)
      write(*,*)'ncycle_max=',ncycle
      write(*,*)'ncycle_real=',mcycle
      write(*,*)
      write(*,*)'hour_max=',hour_max
      write(*,*)'hour_real=',atime
      write(*,*)
      write(*,* )'Start-Time',mvTs,mvTe,ehbTs,ehbTe
       print '( a15, es10.3e2 )', 'Speedups      : ',ehbTe
!!FIX      write(*,*)'ending time: ',fdate()
      
      stop
      return
      end
