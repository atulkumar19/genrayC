
c        ********************** outpt***********************
c        *                      -----                       *
c        *   output is used by drkgs2 as an	                *
c        *   output subroutine. it has not to change the    *
c        *   values of its formal input  parameters.  if    *
c        *   prmt(5) is not equal to zero,  the  control    *
c        *   will be transferred to the main program.       *
c        *       its formal parameters are:                 *
c        *         x, y, dery, ihlf, ndim, prmt             *
c        *                                                  *
c        ****************************************************
c        !   this code prints: values of u(ndim);	    !
c        !   creates the data for 3D code; 		    !
c        !   it writes array u(6) in file( the name of      !
c        !   file is given as  the first parameter in	    !
c        !   genray.in file),                   	    !
c        !   eps- value of Hamiltonian. 		    !
c        !   It controls the conditions: if the ray point   !
c        !   is inside the plasma or not .		    !
c        !   It controls the conditions: if the refractive  !
c        !   index less the cnmax or not .		    !
c        !   It writes the output date in mnemonic.txt file	    !
c        !   It corrects the trajectory for Hamiltonian     !
c        !   conservation 
c        !   It scatters the perpendicular refractive index
c----------------------------------------------------------------------
c         input parameters:,t,u,deru,ihlf,ndim,prmt,iflagh,ht
c         output parameters:
c            if iraystop=1  then end of the j_ray calculation
c            if iraystop=0  then continuation of the j_ray calculations
c            iflagh=(2 rays outside  the plasma after the correction;
c                    1 ray is near the plasma boundary after Runge-Kutta
c                      procedure (it is after reflection)
c                    3 ordinary situation ray is inside the plasma
c                      after correction procedure)
c            u(i) after correction and reflection
c            it prepares and writes  parameters for 3d and onetwo
c            codes
c            it changes the u():
c            1)due to correction procedure that conserves Hamiltonian D=0
c            2)due to boundary reflection in reflection points
c            3)due to n_perp scattering after reflection
c
c
c           For i_ox.eq.1 it calculates antenna vertex coordinates
c           z_st_ox,r_st_ox,phi_st_ox,alpha_st_ox,beta_st_ox
c           and puts these coordinates into  cone_ec
c           For i_ox.eq.2 it creates OX mode conversion jump
c           in the small radial direction and calculates the 
c           OX transmission coefficient: transm_ox
c-----------------------------------------------------------------*
c           this program uses the following functions and	  *
c           subroutins b,gamma1,hamilt1,  bound,prep3d,write3d    *
c                      scatperp                                   *
c-----------------------------------------------------------------*
      subroutine outpt(t,u,deru,ihlf,ndim,prmt,iflagh,iraystop)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'write.i'
      include 'cefield.i'
      include 'cone.i'
c-----the following variable is for N_perp scatterin procedure   
      !!! include 'scatnper.i'
c-----the following to get rma,zma for ox_conversion
      include 'three.i'

      dimension u(*),deru(*),prmt(*),pu(6)
c-----the following include is for hot dispersion with electron+ion, 
c     it has the arrays for the correction subroutines 
      include 'nperpcom.i'
      include 'ions.i'
      
c for test only
      double complex chamilt
cend test

      include 'output.i'
      include 'oxb.i'
    
      real*8 r_old, !major radius at previous point
     &t_old,        !t=poloidal length at previus point
     &ham_old       ! Hamiltonian at previous step
      save r_old,t_old      

      if (first) then
        first=.false.
        r_old=u(2)
      endif

      iraystop=0
      iflagh=3
c-----------------------------------------------------------

      z1=u(1)
      r1=u(2)
      phi1=u(3)
      cnz1=u(4)
      cnr1=u(5)
      cm1=u(6)
      iter=0

      cnphi1=cm1/r1
      cnpar=(bz*cnz1+br*cnr1+bphi*cnphi1)/bmod
     
      cnpar2=cnpar*cnpar
      cn2=cnz1*cnz1+cnr1*cnr1+cnphi1*cnphi1
      cnper=dsqrt(dabs(cn2-cnpar2))
c--------------------------------------------------------------     

c------------------------------------------------------------
c     if the number of the time step nstep_rk is bigger that maxsteps_rk
c     then stop ray calculations
      if(nstep_rk.gt.maxsteps_rk) then
        write(*,*)'**********nstep_rk.gt.maxsteps_rk****************'
	iraystop=1
	return
      end if
      nstep_rk=nstep_rk+1
c--------------------------------------------------------------------
c     if ray is close to resonance point then stop ray calculations
      cnmode=cn(r1,cnz1,cnr1,cm1)
      cnmax=10000.d0
      cnmax=1.d+6

      if(cnmode.gt.cnmax) then
	write(*,*)'***************nmode.gt.cnmax***********'
	iraystop=1
	return
      end if

c------------------------------------------------------------------
c      check that cold plasma X mode is close to UH resonance
c      ann can be EBW conversion 
c-------------------------------------------------------------
      xe=x(u(1),u(2),u(3),1)
      ye=y(u(1),u(2),u(3),1)
      uh=dsqrt(xe+ye*ye)
      uh=xe+ye*ye
      del_uh=1.d-2
      temp_e=temperho(rho,1) 
      rme=9.1094d-28 
      vt_e=dsqrt(2.d0*temp_e*1.6022d-9/rme) !(cm/sec)
                                             ! thermal velocity
      clight=2.99792458d10
      cnper_max_ebw=ye*(clight/vt_e)
      
      if((((id.eq.1).or.(id.eq.2)).and.(uh.gt.1.d0).and.
     & ((uh-1.d0).lt.del_uh))
     &.and.(cnper.gt.cnper_max_ebw)) then
        write(*,*)'********************************************'
        write(*,*)'For cold plasma dispersion relation id=',id
        write(*,*)'the ray is close to upper hybrid resonance'
        write(*,*)' (uh-1.d0).lt.del_uh, uh=',uh,'del_uh=',del_uh
        write(*,*)'N>1, N=',dsqrt(cn2),'cnper=',cnper
        write(*,*)'EBW condition: cnper*(vt_e/clight)/ye=',
     &  cnper*(vt_e/clight)/ye
        write(*,*)'It can be Xmode to close to the UH resonance'
        write(*,*)'Xmode can be transformed to EBW'
        write(*,*)'For Xmode- EBW convertion hot plasma dispersion'
        write(*,*)'can be used id=6,9'
        write(*,*)'The ray calculation stopped' 
        write(*,*)'********************************************'
        iraystop=1
	return
      endif

      if (t.gt.poldist_mx) then
	 write(*,*)'***************t.gt.poldist_mx***********'
      	 iraystop=1
	 return
      end if

c     if Group Velocity, normalized to c, is greater than 1.1
c     (give 10 percent grace) then stop the ray.  [RWH: 030427].

      vgrmods=deru(1)**2+deru(2)**2+(u(2)*deru(3))**2

      if (vgrmods.gt.1.1) then
         write(*,*)
         write(*,*) '*************************************************'
         write(*,*) 'outpt: Stopping ray.'
         write(*,*)' vgroup>1.1,   abs(vgroup) = ',dsqrt(vgrmods)  
         write(*,*) '*************************************************'
         write(*,*)

CENM 1Sep05 -- best to really stop the ray when using the relativistic
C    dispersion relation, otherwise it goes on and on without much progress.
C    vgrmods.gt.1.1 usually when it has nearly reached full depletion of power
C    in the ray anyway.
      	 if (id.eq.14) iraystop=1
      end if

c---------------------------------------------------------------------
c     change of the dispersion function near the cyclotron resonance 
      bmod=b(z1,r1,phi1)
      
      yj=y(z1,r1,phi1,jy_d)

      if(iswitch.eq.1) then
c-------  change of the dispersion function and the absorption subroutines
c         near the cyclotron resonance points       
          call switch_da(yj,del_y,id,iabsorp,idswitch,iabswitch)
      endif
      
c---------------------------------------------------------
c     correction
c---------------------------------------------------------
c     the  switch off the Hamiltonian correction procedure
      if(icorrect.eq.0) goto 11
c---------------------------------------------------------
      epscor=prmt(4)
      bmod=b(z1,r1,phi1)
      gam=gamma1(z1,r1,phi1,cnz1,cnr1,cm1)
      eps=hamilt1(u(1),u(2),u(3),u(4),u(5),u(6))
c      write(*,*)'output after hamilt1 eps',eps
      eps1=eps
  
      cnphi1=cm1/r1
      cnpar=(bz*cnz1+br*cnr1+bphi*cnphi1)/bmod
     
      cnpar2=cnpar*cnpar
      cn2=cnz1*cnz1+cnr1*cnr1+cnphi1*cnphi1
      cnper=dsqrt(dabs(cn2-cnpar2))

      if(dabs(eps).lt.epscor) goto 53

      if ((cnper.ge.0.d0).and.(id.eq.6)) then
c--------correction from the solution n_perp=n_perp(n_parallel)
c        Now it is for id=4,5,6,7
         
        do j=1,nbulk
          massc(j)=dmas(j)
          xc(j)=x(z1,r1,phi1,j)
          yc(j)=y(z1,r1,phi1,j)
          if(j.eq.1) yc(1)=-yc(1)
          tec(j)=tempe(z1,r1,phi1,j)*1.d+3 !(eV) averaged temperature
          tpopc(j)=tpoprho(rho,j)
          vflowc(j)=vflowrho(rho,j)         
        enddo
 
        cnpar=(bz*cnz1+br*cnr1+bphi*cm1/r1)/bmod
        accurcy=epscor 
        naccurc=5
        cnper=dsqrt(dabs((cnz1**2+cnr1**2+(cm1/r1)**2)-cnpar**2))
c        write(*,*)'output before solvnperp cnpar,cnper',cnpar,cnper
        ihermloc=iherm
        
        call solvnperp(nbulk,massc,xc,yc,tec,tpopc,vflowc,cnpar,id,
     *  ihermloc,accurcy,naccurc,cnprim,cnper) !-> get cnper
        write(*,*)'output after solvnperp cnpar,cnper',cnpar,cnper
        write(*,*)'output before correct2 cnz1,cnr1',cnz1,cnr1
        
        call correct2(cnpar,cnper,
     *  cnz1,cnr1,cm1,r1,bz,br,bphi,bmod,  cnznew,cnrnew) !-> get cnznew,cnrnew
        goto 21 ! YuP: skipping correct3()?

cYuP?        eps=1.d-8
cYuP?        itermax=20
c        write(*,*)'output before correct3 z1,r1,phi1',z1,r1,phi1
cYuP?        call correct3(cnpar,cnper,cnz1,cnr1,cm1,z1,r1,phi1,
cYuP?     .  eps,itermax,cnznew,cnrnew,rnew)
c        write(*,*)'output after correct3 rnew,cnznew,cnrnew',
c     .  rnew,cnznew,cnrnew        
cYuP?         u(2)=rnew

 21     continue
        write(*,*)'output after correct2 cnznew,cnrnew',cnznew,cnrnew 
        !pause
        u(4)=cnznew
        u(5)=cnrnew

        goto 11
      endif 
cSAP080829
 60   continue

c--------correction from the problem
c        sum{delt(x_i)**2}=min
c        disp_func(x+delt(x))=0 
	 iter=1
	 icor=0 
c----------------------------------------------------------------------  
ctest
         bmod=b(z1,r1,phi1)
         gam=gamma1(z1,r1,phi1,cnz1,cnr1,cm1)
         eps=hamilt1(u(1),u(2),u(3),u(4),u(5),u(6))

 52	 continue
         dip=1.d0
         call dddrz1(t,u,deru)
         sum=deru(1)**2+deru(2)**2
 12      dlambd=2.d0*eps/sum*dip
c-----------------------------------------------------------
	 do 50 i=1,2
 50	   pu(i)=u(i)

	 do 51 i=4,5
 51      pu(i)=u(i)-0.5*dlambd*deru(i-3)

         pu(3)=u(3)
         pu(6)=u(6)
c-----------------------------------------------------------
         z1=pu(1)
         r1=pu(2)
         phi1=pu(3)
         cnz1=pu(4)
         cnr1=pu(5)
         cm1=pu(6)
c         write(*,*)'output  before boundc r1,z1',r1,z1
c         write(*,*)'output  cnz1,cnr1,cm1', cnz1,cnr1,cm1
         call boundc(z1,r1,iboundc)

c         write(*,*)'output after boundc iboundc',iboundc

	 if (iboundc.eq.1) then
c           write(*,*)'in output iboundc=1 lflagh=2'
c	   write(*,*)'outside the plasma after the correction'
c---------------------------------------------------------------------
c          corrected ray point is outside the plasma
c          it is necessary to reduce the time step in the Runge -Kutta
c          procedure  and to recalculate the array u(i)
c----------------------------------------------------------------------
           iflagh=2
c           goto 120

c          after the correction step ray is out of plasma in up()
c          stop the correction procedure at the last correction step u()
           iflagh=3
           goto 11
	 endif

         bmod=b(z1,r1,phi1)
         gam=gamma1(z1,r1,phi1,cnz1,cnr1,cm1)
         eps=hamilt1(pu(1),pu(2),pu(3),pu(4),pu(5),pu(6))

       if(eps.gt.1.d120) then
           write(*,*)'in output eps.gt.1.d20 iter=',iter
           write(*,*)'Hamiltonian correction procedure stoped'
           write(*,*)'Hamiltonian=',eps
           goto 53
	 end if
  
	 if (dabs(eps).gt.dabs(eps1)) then
	     write(*,*)'in output eps.gt.eps1,dip',eps,eps1,dip             
	     dip=0.5d0*dip
	     goto 12
	 end if
         
c-------------------------------------------------------------
	 eps1=eps
c-------------------------------------------------------------
	 do 14 i=1,ndim
 14      u(i)=pu(i)

	 iter=iter+1
	 if(dabs(eps).lt.epscor) goto 53
         
	 if(iter.gt.20)then
	   write(*,*)'Hamiltonian correction procedure
     1	   made 20 iterations and stopped , Hamiltonian=',eps
	   goto 53
	 end if

	 goto 52
 53      continue

c--------------------------------------------------------------
c     end of correction
c--------------------------------------------------------------
 11   continue
c--------------------------------------------------------------
c     measure error in the dispersion realation
c     If D/(N|gradD|) > toll_hamilt stop ray calculation
c--------------------------------------------------------------

c      write(*,*)'output.f before refractive_index_relative_error'

      call refractive_index_relative_error(u(1),u(2),u(3),u(4),u(5),u(6)
     &,iraystop)

c      write(*,*)'output.f after refractive_index_relative_error iraystop
c     &',iraystop

      if(iraystop.eq.1) then
        write(*,*)
         write(*,*) '*************************************************'
         write(*,*) 'outpt: Stopping ray.'
         write(*,*) ' D/(N|gradD|) > toll_hamilt'  
         write(*,*) '*************************************************'
         write(*,*)
         return
      endif
c------------------------------------------------------------------
c     Creates the jump of the ray point throw the OX mode conversion
c     area where X_e=1 (V_perp=0)
c-------------------------------------------------------------------
      if (i_ox.eq.2) then
        write(*,*)'was_not_ox_conversion',was_not_ox_conversion
        i_call_prep3d_in_output_at_i_ox_conversion_eq_1=0
        if (was_not_ox_conversion) then
          rma_loc=rma
          zma_loc=zma
          eps_xe_loc=eps_xe
          
c---------calculate the coordinates of X mode cutoff point
c---------r_x,z_x,phi_x,cnr_x,cnz_x,cm_x 
cyup          write(*,*)'output'
cyup          write(*,*)'before ox_conversion u(1),u(2),u(3),u(4),u(5),u(6)'
cyup     &    ,u(1),u(2),u(3),u(4),u(5),u(6)

          call ox_conversion(u(2),u(1),u(3),u(5),u(4),u(6),
     &    rma_loc,zma_loc, !temporally
     &    eps_xe_loc,
     &    r_x,z_x,phi_x,cnr_x,cnz_x,cm_x,i_ox_conversion_loc)
          i_ox_conversion=i_ox_conversion_loc

          write(*,*)'i_ox_conversion=',i_ox_conversion
          write(*,*)'after ox_conversion u(1),u(2),u(3),u(4),u(5),u(6)'
     &    ,u(1),u(2),u(3),u(4),u(5),u(6)
          
          if (i_ox_conversion.eq.1) then
c           calculate transmission coefficient for OX mode
c           conversion:transm_ox

c--------------------------------------------------------
c           The following data are for write.i to prepare 
c           the output data for mnemonic.nc file            
              
            Y_abs=dabs(y(u(1),u(2),u(3),1))
            cn_par_optimal=dsqrt(Y_abs/(Y_abs+1)) !optimal N parallel
                                                  !for OX conversion
            bmod=b(u(1),u(2),u(3))
            cnpar_ox=(bz*u(4)+br*u(5)+bphi*u(6)/u(2))/bmod ! N_parallel
                                                           ! before OX
                                                           ! conversion
            cn_b_gradpsi=(u(5)*bphi*dpdzd-
     &      (u(6)/u(2))*(br*dpdzd-bz*dpdrd)-u(4)*bphi*dpdrd)/
     &   dsqrt((bphi*dpdzd)**2+(br*dpdzd-bz*dpdrd)**2+(bphi*dpdrd)**2)


c-------------------------------------------------------
   
c-----------transmission coefficient in the point before the jump 
            call transmit_coef_ox(u(1),u(2),u(3),u(4),u(5),u(6),
     &                           transm_ox)
            write(*,*)'output OX mode conversion coefficient transm_ox'
     &      ,transm_ox

c--------------------------------------------------------------------
c           write the data fo O mode  before jump
            i_call_prep3d_in_output_at_i_ox_conversion_eq_1=1
c            write(*,*)'output before ox jump before prep3d'
c            write(*,*)'u(1),u(2)',u(1),u(2)
            call prep3d(t,u,deru,iraystop)
c--------------------------------------------------------------------

            u(1)=z_x
            u(2)=r_x
            u(3)=phi_x
            u(4)=cnz_x
            u(5)=cnr_x
            u(6)=cm_x
            write(*,*)'output after jump before prep3d u= ',
     &     u(1),u(2),u(3),u(4),u(5),u(6)

c------------------------------------------------------------------------
c           write the data for x mode  after jump
            i_call_prep3d_in_output_at_i_ox_conversion_eq_1=2
            call prep3d(t,u,deru,iraystop)
c--------------------------------------------------------------------
            was_not_ox_conversion= .false.
          endif ! i_ox_conversion.eq.1
        endif   ! was_not_ox_conversion 
      endif     ! i_ox.eq.2
c-------------------------------------------------------------------
    
 40   continue
     
      
 41   continue

      r_old=u(2)
      t_old=t
      i_output_data=0 !this time step is not the output step

     
      if(t.lt.prmt(7))goto 20
      i_output_data=1 !this time step is the output step
c-------------------------------------------------------
      write(*,*)'tau=',t,' 1/y(1)=',1.d0/y(z1,r1,phi1,1),
     &'1/y(2)=',1.d0/y(z1,r1,phi1,2)

      
c     The work on the data preparation for the output
c     at given time steps
c---------------------------------------------------------
c     denormalization of ray point coordinates
c---------------------------------------------------------
      uout1=r0x*u(1)
      uout2=r0x*u(2)
      uout3=u(3)
      uout4=u(4)
      uout5=u(5)
      uout6=u(6)*r0x
c-----------------------------------------------------------
c      write(*,30)uout2,uout1,uout3,uout4,uout5,uout6
      write(*,*)uout2,uout1,uout3,uout4,uout5,uout6

c      write(*,130)uout2,uout1,uout3,uout4,uout5,uout6
130   format(3x,6(' ',e16.9))

c      write(i1_,110)uout2,uout1,uout3,uout4,uout5,uout6
      
c      prmt(7)=prmt(7)+prmt(6)
      prmt(7)=(dint(t/prmt(6))+1)*prmt(6)

cSAP090515
      write(*,*)'t,prmt(6),prmt(7)',t,prmt(6),prmt(7)

 110  format(3x,6(' ',e13.6))
c***************************************************
      bmod=b(z1,r1,phi1)

c      write(*,*)'output bef eps z1,r1,phi1',z1,r1,phi1
c      write(*,*)'output bef eps bz,br,bphi,bmod',bz,br,bphi,bmod
c      write(*,*)'output bef eps cnz1,cnr1,cm1',cnz1,crn1,cm1
c       write(*,*)'output rho',rho
c     gam=gamma1(z1,r1,phi1,cnz1,cnr1,cm1)
      gam=gamma1(u(1),u(2),u(3),u(4),u(5),u(6))
      ds2=dsin(gam)**2
      dc2=dcos(gam)**2
c      write(*,*)'output bef eps gam,ds,dc',gam,dsin(gam),dcos(gam)    
      call abc(u(1),u(2),u(3),ds2,dc2,ad,bd,cd)

c      write(*,*)'output bef eps u(1),u(2),u(3)',u(1),u(2),u(3)
c      write(*,*)'output bef eps u(4),u(5),u(6)',u(4),u(5),u(6)
c      write(*,*)'output bef eps bz,br,bphi,bmod',bz,br,bphi,bmod
      
      eps=hamilt1(u(1),u(2),u(3),u(4),u(5),u(6))
      write(*,*)'in output epshamilt1=',eps

      cn2=cnz1*cnz1+cnr1*cnr1+(cm1/r1)**2

c testc      write(*,*)'output bef eps z1,r1,phi1',z1,r1,phi1
c      write(*,*)'output bef eps bz,br,bphi,bmod',bz,br,bphi,bmod
c      write(*,*)'output bef eps cnz1,cnr1,cm1',cnz1,crn1,cm1
c       write(*,*)'output rho',rho
      cnpar=(cnz1*bz+cnr1*br+(cm1/r1)*bphi)/bmod
      cnper2=cn2*ds2
      cnper=dsqrt(cnper2)
      xe=x(u(1),u(2),u(3),1)
      ye=y(u(1),u(2),u(3),1)
      uh=dsqrt(xe+ye*ye)
      uh=xe+ye*ye
      write(*,*)'output xe,ye,uh',xe,ye,uh

      if (i_uh_switch.eq.1) then
         !change the output step prmt(6) for prmt6_uh_switch 
         if(uh.lt.uh_switch) then
           prmt(6)=prmt6_uh_switch
         else
         ! regular output step
           prmt(6)=prmt6
         endif
      endif


      if (i_power_switch_resonance.eq.1) then
         !change the output step prmt(6) for prmt6_power_switch_resonance
         do k=1,n_power_switch_resonance 
           if(dabs(ye-y_power_switch_resonance(k)).lt.
     &       del_y_power_switch_resonance) then
             prmt(6)=prmt6_power_switch_resonance
           else
           ! regular output step
             prmt(6)=prmt6
           endif
         enddo
      endif  


c      write(*,*)'output cnper,cnpar',cnper,cnpar
      do i=2,nbulk
         xi=x(u(1),u(2),u(3),i)
         yi=y(u(1),u(2),u(3),i)
c         write(*,*)'output i,xi,yi',i,xi,yi
      enddo      
c end test
 568  format(3x,'eps=',d13.6)
     
      call prep3d(t,u,deru,iraystop)

      if (iraystop.eq.1) then
	 return
      end if

      if (nrayelt.eq.nrelt) then
	 write(*,*)'***************nrayelt=nrelt***********'
      	 iraystop=1
	 return
      end if

  20  continue
c---------------------------------------------------
c     end of output data preparation at the given time steps
c--------------------------------------------------
  30  format(1x,6(1x,e11.4))
c-----------------------------------------------------------
c     control of the reflection moment
c-----------------------------------------------------------
c     this call is for the calculation of
c     dzdt=deru(1) and drdt=deru(2).
c     subroutine bound uses deru(1) and deru(2)	to detemine
c     if the ray point goes into or out the plasma.
      call rside1(u,deru)
c-----------------------------------------------------------
      write(*,*)' 20 output before bound u(2),u(1)',u(2),u(1)
      call bound(u(1),u(2),u(3),u(4),u(5),u(6),iflref,
     &      z_ref,r_ref,phi_ref,cnzref,cnrref,cmref,
     &      ibound,deru(1),deru(2))

      write(*,*)'20 output after bound ibound,iflref,u(2),u(1)',
     &           ibound,iflref,u(2),u(1)

      if (iflref.eq.1) then

        write(*,*)'the data in reflection point before reflection'
        write(*,*)'z=',u(1),'r=',u(2),'phi=',u(3)
        write(*,*)'cnz=',u(4),'cnr=',u(5),'cm=',u(6)
        bmod=b(u(1),u(2),u(3))    
        write(*,*)'cnpar=',(u(4)*bz+u(5)*br+u(6)*bphi/u(2))/bmod
        write(*,*)'cn2=',u(4)**2+u(5)**2+(u(6)/u(2))**2
        write(*,*)'cnper=',dsqrt((u(4)**2+u(5)**2+(u(6)/u(2))**2)-
     &                 ((u(4)*bz+u(5)*br+u(6)*bphi/u(2))/bmod)**2)
        write(*,*)'after reflection'

        u(1)= z_ref
        u(2)= r_ref
        u(3)= phi_ref

        write(*,*)'z_ref=',z_ref,'r_ref=',z_ref,'pfi_ref=',phi_ref
        write(*,*)'cnzref=',cnzref,'cnrref=',cnrref,'cmref=',cmref
        write(*,*)'cnpar=',(cnzref*bz+cnrref*br+u(6)*bphi/u(2))/bmod
        write(*,*)'cn2=',cnzref**2+cnrref**2+(u(6)/u(2))**2
        write(*,*)'cn2=',cnzref**2+cnrref**2+(cmref/u(2))**2
        write(*,*)'cnper=',dsqrt((cnzref**2+cnrref**2+(u(6)/u(2))**2)-
     &                 ((cnzref*bz+cnrref*br+u(6)*bphi/u(2))/bmod)**2)
        write(*,*)'cnper=',dsqrt((cnzref**2+cnrref**2+(cmref/u(2))**2)-
     &                 ((cnzref*bz+cnrref*br+cmref*bphi/u(2))/bmod)**2)

      endif

      if ((i_ox.eq.1).and.(iflref.eq.1)) then
c--------calculate the EC cone vertex coordinates
c        for the optimal OX mode conversion
           bmod=b(u(1),u(2),u(3))     
           b_z0=bz
           b_r0=br
           b_phi0=bphi
           dpsi_dz=dpdzd
           dpsi_dr=dpdrd
           z0=u(1)
           r0=u(2)
           phi0=u(3)
           cnz0=u(4)
           cnr0=u(5)
           cm0=u(6)
           cnphi0=cm0/r0              
           write(*,*)'output i_ox=1'
           write(*,*)'z0,r0,phi0',z0,r0,phi0
           write(*,*)'cnz0,cnr0,cm0,cnphi0',cnz0,cnr0,cm0,cnphi0
           
c------test
c           write(*,*)'bz,br,bphi',bz,br,bphi
            write(*,*)'output.f before call antenna_vertex'
           cnpar_test=(cnz0*bz+cnr0*br+cnphi0*bphi)/bmod
           write(*,*)'output cnpar_test',cnpar_test
           cn0_s=cnz0**2+cnr0**2+cnphi0**2
           cnper_test=dsqrt(cn0_s-cnpar_test**2)
c           write(*,*)'dsqrt(cn0_s)',dsqrt(cn0_s)
           write(*,*)'output cnper_test',cnper_test
           sin_phi=dsin(phi0)
           cos_phi=dcos(phi0)
           cnx=cnr0*cos_phi-cnphi0*sin_phi
           cny=cnr0*sin_phi+cnphi0*cos_phi
           cnz=cnz0  
           write(*,*)'cnx,cny,cnz',cnx,cny,cnz
           call ninit_ec(z0,r0,phi0,cnx,cny,cnz,cnteta,cnphi0)
           write(*,*)'cnteta,cnphi0',cnteta,cnphi0
           write(*,*)'cnrho',dsqrt(cn0_s-cnteta**2-cnphi0**2) 
           ioxm_loc=ioxm
           ioxm=1  
           call cninit12(z0,r0,phi0,cnpar_test,cnteta,cnphi0,
     1                  cnz_t,cnr_t,cm_t,iraystop)
           ioxm=ioxm_loc 
           write(*,*)'output.f after cninit12 cnz_t,cnr_t,cm_t',
     &     cnz_t,cnr_t,cm_t
c----------endtest
           icone=icone_loc
           call antenna_vertex(u(2),u(3),u(1),-u(5),-u(6)/u(2),-u(4),
     &     b_r0,b_phi0,b_z0,dpsi_dz,dpsi_dr,
     &     r_st,phi_st,z_st,alpha_st,beta_st,icone)

c----------put the data into cone.i
           r_st_ox=r_st
           phi_st_ox=phi_st
           z_st_ox=z_st
           alpha_st_ox=alpha_st
           beta_st_ox=beta_st

           write(*,*)'oxb r_st,phi_st,z_st,alpha_st,beta_st',
     &     r_st,phi_st,z_st,alpha_st,beta_st
 
      endif ! ((i_ox.eq.1).and.(iflref.eq.1)  

      u(4)=cnzref
      u(5)=cnrref
      u(6)=cmref

cSAP081010
      if (iflref.eq.1) call prep3d(t,u,deru,iraystop)
c----------------------------------------------------------
c     call b() to calculate the small radius rho (inside b())
c     the resulting will be in common block  one.i
c-----------------------------------------------------------
      bmod=b(u(1),u(2),u(3)) 
      if (iflref.eq.1) then         
         iflagh=1
      else
c---------------------------------------------------------------
      end if

      rhooldsc=rho

      if (irefl.ge.ireflm) then
	 write(*,*)'irefl.ge.ireflm u(2)',u(2)
         write(*,*)'i_output_data ',i_output_data 
         if (i_output_data.eq.0) then
             call prep3d(t,u,deru,iraystop)
         endif
 	 iraystop=1            
      else
 	 iraystop=0
      end if

 120  continue
c----------------------------------------------------
c     for control phi_deviation for one time step
c     write(*,*)'in output phiold,u(3)',phiold,u(3)
      phiold=u(3)  ! phiold is in common write.i
c----------------------------------------------------
c      write(*,*)'end outpt'
      return
      end




      subroutine correct2(cnpar,cnper,
     *cnz,cnr,cm,r,bz,br,bphi,bmod,cnz1,cnr1)

c      subroutine correct2(xe,ye,te_kev,tpope,vflowe,cnpar,id,
c     *ihermloc,accurcy,naccurc,cnprim,cnper,
c     *cnz,cnr,cm,r,bz,br,bphi,bmod,cnz1,cnr1)
c-----------------------------------------------------------
c     it calculates new values of N_perp coordinate
c     to get the Hamiltonian conservation
c     from the solution Nperp(N_parallel)
      implicit none
c-----input
      double precision cnpar,cnper,
     *cnz,cnr,cm,r,bz,br,bphi,bmod

c-----output
      double precision cnr1,cnz1
c-----local
      double precision cnper2,cnpar2,d1,d2,a0,a1,a2,det,
     *delnm,delnp,cnz1p,cnz1m,cnr1p,cnr1m,cnperold
      integer inewn
c     inewn=0 ! There were used the old values(cnz,cnr) inside 
c               this subbroutine
c          =1 ! There were used the new corrected
c               values(cnz1,cnr1) inside this subroutine

cfor test
      double precision oldnpol,newnpol,cn2new,bpol,bplnpl,
     *cosnb_pl
      
c------------------------------------------------------------  
c     vector{cnper}={nz-(bz/bmod)N_par,
c                    nr-(br/bmod)N_par
c                    nphi-(bphi/bmod)N_par}
c     we will determine the new (cnz1,cnr1) from the following system
c     cnr1*br+cnz1*bz=cnpar*bmod-(cm/r)*bphi=d1
c     cnr1**2+cnz1**2=cnpar**2+cnper**2-(cm/r)**2=d2
c-----------------------------------------------------------------
      cnper2=cnper*cnper
      cnpar2=cnpar*cnpar    
      d1=cnpar*bmod-(cm/r)*bphi
      d2=cnpar2+cnper2-(cm/r)**2

      oldnpol=dsqrt(cnz*cnz+cnr*cnr)
      cn2new=cnper2+cnpar2
      newnpol=dsqrt(cn2new-(cm/r)**2)
      bpol=dsqrt(bz*bz+br*br)/bmod
      bplnpl=cnpar-(cm/r)*(bphi/bmod)
c      write(*,*)'correct2 oldnpol,newnpol',oldnpol,newnpol
      cosnb_pl=bplnpl/newnpol
      inewn=0
      if (d2.lt.0.d0) then
         write(*,*)'correct2 d2=cnpar2+cnper2-(cm/r)**2<0
     *   correction was not made'       
         goto 10
      else
         if ((br*br+bz*bz).eq.0.d0) then
c           write(*,*)'correct2 (br*br+bz*bz).eq.0.d0'
            cnperold=dsqrt(cnz*cnz+cnr*cnr)
c           write(*,*)'coorect2 cnper,cnperold',cnper,cnperold
            if(cnperold.eq.0.d0) then
               cnz1=cnper/dsqrt(2.d0)
               cnr1=cnper/dsqrt(2.d0)
            else 
               cnz1=cnz*cnper/cnperold
               cnr1=cnr*cnper/cnperold
            endif
            inewn=1
c            write(*,*)'correct2 cnz1,cnr1',cnz1,cnr1
            goto 10
         endif

         if(br.ne.0.d0) then
c----------equation a2*cnz1**2-2*a1*cnz1+a0=0

           a2=1.d0+(bz/br)**2
           a1=(bz*d1)/(br*br)
           a0=(d1/br)**2-d2 
           det=a1*a1-a2*a0

           if(det.lt.0.d0) then
             write(*,*)'correct2 br.ne.0 det<0
     *       correction in correct2 was not made'
             goto 11
           else
             cnz1m=(a1-dsqrt(det))/a2
             cnz1p=(a1+dsqrt(det))/a2
             cnr1m=(d1-bz*cnz1m)/br
             cnr1p=(d1-bz*cnz1p)/br
             inewn=1        
           endif
         endif
 11      continue

         if(bz.ne.0.d0) then
c----------equation a2*cnr1**2-2*a1*cnr1+a0=0
           a2=1.d0+(br/bz)**2
           a1=(br*d1)/(bz*bz)
           a0=(d1/bz)**2-d2 
           det=a1*a1-a2*a0

           if(det.lt.0.d0) then
             write(*,*)'correct2 bz.ne.0 det<0
     *       correction in correct2 was not made'
             goto 10
           else
             cnr1m=(a1-dsqrt(det))/a2
             cnr1p=(a1+dsqrt(det))/a2
             cnz1m=(d1-br*cnr1m)/bz
             cnz1p=(d1-br*cnr1p)/bz
             inewn=1	     
           endif
         endif

      endif

      delnm=(cnz-cnz1m)**2+(cnr-cnr1m)**2
      delnp=(cnz-cnz1p)**2+(cnr-cnr1p)**2
c      write(*,*)'output_correct2 delnm,delnp',delnm,delnp
      if (delnm.le.delnp) then
        cnz1=cnz1m
        cnr1=cnr1m
      else
        cnz1=cnz1p
        cnr1=cnr1p
      endif

 10   continue

      if (inewn.eq.0) then
c-------the correction was not made
        cnz1=cnz
        cnr1=cnr
      endif

      return
      end


c-------------------------------------------------------------
c     calculation n_perp(n_par) for Mazzucato and Forest codes
c     It uses the input cnper=Re(n_perp) for the first iteration
c     It calculates cnprim=Im(N_perp) for forest case (id=6)
      subroutine solvnperp(nbulk,dmas,x_ar,y_ar,te_ev_ar,tpop_ar,
     .vflow_ar,
     .cnpar,id,
     .ihermloc,accurcy,naccurc,cnprim,cnper)
      implicit none
c-----input
      integer nbulk
      double precision dmas(*),x_ar(*),y_ar(*),te_ev_ar(*),tpop_ar(*),
     .vflow_ar(*),
     .cnpar,cnper,cnprim,accurcy
      integer id,ihermloc,naccurc
c-----output 
c     cnper,cnprim   
c-----local
      double precision step,hamr,hamrp,hamrm,dhamrdnr
      double precision cnperp,cnperm,dcnper,hamold,cnpernew,cnpernw1
      double complex chamilt
      integer iter
c     calculation cnper from the solution of the equation
c     dhamr/d_nper*delnper=-hamr

      cnprim=0.d0
      step=1.d-7
      iter=0
c      write(*,*)'solvnper xe,ye,te_kev',xe,ye,te_kev
c      write(*,*)'solvnper cnpar,cnper,id,naccurc',cnpar,cnper,id,naccurc
      cnpernew=cnper !initial value
      cnpernw1=cnpernew

 10   continue

c      write(*,*)'1 solvnper cnpar,cnper,id',cnpar,cnper,id
c      write(*,*)'solvnperp bef dispfun1 cnpar,iter,cnpernew',
c     *cnpar,iter,cnpernew

      call dispfun_output(nbulk,dmas,x_ar,y_ar,te_ev_ar,tpop_ar,
     &vflow_ar,
     *cnpar,cnpernew,cnprim,id,ihermloc,
     *chamilt)
      hamr=dreal(chamilt)

c      write(*,*)'solvnperp after dispfun1 iter,hamr,cnpernw1,cnpernew'
c     *,iter,hamr,cnpernw1,cnpernew

      if (iter.eq.0) then 
        hamold=hamr
      else        
        if (dabs(hamold).le.dabs(hamr)) then
          write(*,*)'solvnperp dabs(hamold)<dabs(hamr) iter=',iter
          write(*,*)'solution Nperp not found hamold,hamr',hamold,hamr
          write(*,*)'cnper=cnpernw1',cnpernw1
          cnper=cnpernw1     
          goto 20
        else
          cnper=cnpernew
          cnpernw1=cnpernew
        endif
      endif

      if (dabs(hamr).le.accurcy) goto 20

      cnperp=cnper+step*cnper
      call dispfun_output(nbulk,dmas,x_ar,y_ar,te_ev_ar,tpop_ar,
     &vflow_ar,
     *cnpar,cnperp,cnprim,id,ihermloc,chamilt)
      hamrp=dreal(chamilt)
      
      cnperm=cnper-step*cnper
      call dispfun_output(nbulk,dmas,x_ar,y_ar,te_ev_ar,tpop_ar,
     &vflow_ar,
     *cnpar,cnperm,cnprim,id,ihermloc,chamilt)
      hamrm=dreal(chamilt)
         
      dhamrdnr=(hamrp-hamrm)/(2.d0*step*cnper)
       
c      write(*,*)'dhamrdnr,hamr',dhamrdnr,hamr
      if (dhamrdnr.eq.0.d0) then
        write(*,*)'solvnperp dhamdnnr=0'
        goto 20
      else    
        dcnper=-hamr/dhamrdnr
        cnpernew=cnper+dcnper
        iter=iter+1

        if(cnpernew.lt.0.d0) then
          write(*,*)'output.f in solvnper cnpernew<0'
          cnpernew=cnper
          goto 20
        endif

        if (iter.gt.naccurc) then
           write(*,*)'solvnper iter>naccurc'
           cnper=cnpernew
           goto 20
        endif

      endif
      hamold=hamr
      goto 10

 20   continue
c      write(*,*)'solvnper bef end cnper,iter',cnper,iter
      return
      end



      subroutine dispfun_output(nbulk,dmas,x_ar,y_ar,te_ev_ar,tpop_ar,
     .vflow_ar,
     .cnpar,cnper,cnprim,id,
     .ihermloc,chamilt)
c     It calculates the complex dispersion function chamilt
c     for different tensors (NOw only for id=4,5 and 6)
c     id=1
c     id=2
c     id=3 Appleton-Hartry
c     id=6 Forest code
c-----it calculates Imaginary part cnprim=Im(N_perp) for Forest code
      implicit none
c-----input
      integer nbulk
      double precision dmas(*),x_ar(*),y_ar(*),te_ev_ar(*),
     .tpop_ar(*),vflow_ar(*)      
      double precision xe,ye,te_kev,tpope,vflowe
      double precision cnpar,cnper,cnprim
      integer id,ihermloc
c-----output
      double complex chamilt
c     cnprim : is a solution for forest code
c-----external dhot,complx1
      double complex dhot_sum       
c-----local
      double precision hamr,hami,mode
      double precision mox2de,mu
      double precision cnper2,te_ev,cnparp,cnperp
      double complex cpp,sol
      double complex ceps(3,3)
   
      if (((id.eq.1).or.(id.eq.2)).or.(id.eq.3)) then
         write(*,*)'dispfun id=',id,'now id can be id=4,=5 or 6'
         stop         
      endif
      
      if (id.eq.6) then
c--------Hot dispersion function
c        electron + ions non-relativistic plasma
         cnparp=cnpar
         cnperp=cnper
         hamr=dreal(dhot_sum(nbulk,dmas,x_ar,y_ar,te_ev_ar,tpop_ar,
     .   vflow_ar,cnpar,cnper,1,ceps))
         hami=0.d0
         chamilt=dcmplx(hamr,hami)
         goto 100
      endif !id=6 Hot non-relativistic dispersion

 100  continue
      return
      end




      subroutine switch_da(y,del_y,id,iabsorp,idswitch,iabswitch)
c-----switches on(off) the dispersion function and the absorption near the cyclotron
c     resonance points 
c     input:
c     y    (omega_c/omega) algebraic
c     del_y the distance in Y from the resonace to swith on dispersion and absorption
c     id        the dispersion far from reconace
c     iabsorp   the absorpsion far from the resonance
c     idswitch  the dispersion near resonance
c     iabswitch  the absorpsion near resonance
c     output: as a local save variable
c     iswitch !=0 before the first change of  dispersion
c             !=1 after first changing of the dispersion
      implicit none
c-----input
      double precision y,del_y !
      integer id,iabsorp,idswitch,iabswitch     
c-----local
      integer iy
      integer iswitch,id_old,iabs_old 
      save iswitch,id_old,iabs_old

      iy=nint(1.d0/y)
c      write(*,*)'1 iswitch',iswitch,'y',y,'iy',iy,'del_y',del_y
      if((dabs(1.d0/y-dfloat(iy)).le.del_y).and.(iswitch.ne.1)) then
c------ the change of the dispersion function and absorption near the cyclotron resonance
        iswitch=1  
        id_old=id
        id=idswitch       ! new
        iabs_old=iabsorp
        iabsorp=iabswitch ! new 
c        write(*,*)'in output:  switch_da'
c        write(*,*)'id=',id,'iabsorp=',iabsorp
      endif
c      write(*,*)'2 iswitch',iswitch,'y',y
      if ((dabs(1.d0/y-dfloat(iy)).ge.del_y).and.(iswitch.eq.1)) then
c------ the back switch of the dispersion function after ec resonance 
        iswitch=0
        id=id_old
        iabsorp=iabs_old
c        write(*,*)'in output: switch_da back switch of dispersion'
c        write(*,*)'id=',id,'iabsorp=',iabsorp
      endif         
c      write(*,*)'output.f 3 iswitch,id',iswitch,id
      return
      end




      subroutine correct3(cnpar,cnper,cnz,cnr,cm,z,r,phi,
     .eps,itermax,cnz1,cnr1,r1)
c     calculates the corrected values 
c     cnz1=cnz+delnz,cnr1=cnr+delnr,r1=r+delr from the problem
c     J=delnz**2+delnr**2_delr**2
c     min(J)
c     N_par(Nz+delnz,Nr+delnz,r+delr)=cnpar
c     N_perp(Nz+delnz,Nr+delnr,r+delnr)=cnperp

c-----input
      double precision cnpar,cnper,cnz,cnr,cm,z,r,phi
      double precision eps      !the accuracy of iterations
      integer itermax           !max number of steps in iterations
c-----output
      double precision cnz1,cnr1,r1
c-----locals
      double precision dr,dnz,dnr
      double precision f1,df1ddnz,df1ddnr,df1ddr,
     .f2,df2ddnz,df2ddnr,df2ddr       
      integer iter
      double precision derg(3,3),g123ar(3),det,det1,det2,det3,norma
      double precision bzl,brl,bphil,bmodl,
     .dbzdzl,dbzdrl,dbrdzl,dbrdrl,
     .dbphdzl,dbphdrl,dbmdzl,dbmdrl

      iter=0
      cnz1=cnz
      cnr1=cnr
      r1=r
      dnz=0.d0
      dnz=0.d0
      dr=0.d0
      
 10   continue
      write(*,*)'output correct3 z,r,phi',z,r,phi
      write(*,*)'output correct3 cnz,cnr,cm',cnz,cnr,cm

      call bcomp(z,r,phi,bzl,brl,bphil,bmodl,
     .dbzdzl,dbzdrl,dbrdzl,dbrdrl,
     .dbphdzl,dbphdrl,dbmdzl,dbmdrl)
  
      call g123(cnz1,cnr1,cm,z,r1,phi,dnr,dnz,dr,cnpar,cnper,
     .bzl,brl,bphil,bmodl,dbzdrl,dbrdrl,dbphdrl,dbmdrl,
     .g123ar)
      norma=g123ar(1)**2+g123ar(1)**2+g123ar(1)**2

      if(norma.lt.eps) goto 20

      call drvg123(cnz1,cnr1,cm,z,r1,phi,dnr,dnz,dr,cnpar,cnper,
     .bzl,brl,bphil,bmodl,dbzdrl,dbrdrl,dbphdrl,dbmdrl,
     .derg)
      
      det=derg(1,1)*(derg(2,2)*derg(3,3)-derg(2,3)*derg(3,2))-
     .derg(1,2)*(derg(2,1)*derg(3,3)-derg(2,3)*derg(3,1))+
     .derg(1,3)*(derg(2,1)*derg(3,2)-derg(2,2)*derg(3,1))
      
      det1=-g123ar(1)*(derg(2,2)*derg(3,3)-derg(2,3)*derg(3,2))+
     .g123ar(2)*(derg(1,2)*derg(3,3)-derg(1,3)*derg(3,2))-
     .g123ar(3)*(derg(1,2)*derg(2,3)-derg(1,3)*derg(2,2))

      det2=g123ar(1)*(derg(2,1)*derg(3,3)-derg(2,3)*derg(3,1))-
     .g123ar(2)*(derg(1,1)*derg(3,3)-derg(1,3)*derg(3,1))+
     .g123ar(3)*(derg(1,1)*derg(2,3)-derg(2,1)*derg(1,3))

      det3=-g123ar(1)*(derg(2,1)*derg(3,2)-derg(2,2)*derg(3,1))+
     .g123ar(2)*(derg(1,1)*derg(3,2)-derg(1,2)*derg(3,1))-
     .g123ar(3)*(derg(1,1)*derg(2,2)-derg(1,2)*derg(2,1))

      if (det.eq.0.d0) then 
        write(*,*)' in output.d correct3 det=0'
        stop
      endif

      dnz=det1/det
      dnr=det2/det
      dr=det3/det
      iter=iter+1

      if(iter.gt.itermax) then     
        write(*,*)'in correct3 iter>itermax'
        goto 20
      endif

      goto 10
      
 20   continue

      cnz1=cnz+dnz
      cnr1=cnr+dnr
      r1=r+dr

      return
      end




      double precision function f1(cnz,cnr,cm,r,bz,br,bphi,
     .bmod,cnpar)
      implicit none
      double precision cnz,cnr,cm,r,bz,br,bphi,bmod,cnpar
      f1=cnz*bz+cnr*br+cm*bphi/r-cnpar*bmod
      return
      end

      double precision function f2(cnz,cnr,cm,r,cn)
      implicit none
c-----cn**2=Npar**2+Nper**2  
      double precision cnz,cnr,cm,r,cn
      f2=cnz**2+cnr**2+(cm/r)**2-cn**2
      return
      end

      double precision function df1dnz(bz)
      implicit none
      double precision bz
      df1dnz=bz
      return
      end

      double precision function df1dnr(br)
      implicit none
      double precision br
      df1dnr=br
      return
      end

      double precision function df1dr(cnz,cnr,cm,r,bphi,bmod,cnpar,
     .dbzdr,dbrdr,dbphdr,dbmdr)
      implicit none
      double precision cnz,cnr,cm,r,bphi,bmod,cnpar,
     .dbzdr,dbrdr,dbphdr,dbmdr 
      df1dr=cnz*dbzdr+cnr*dbrdr-cm*bphi/r**2+
     .cm*dbphdr/r-cnpar*dbmdr
      return
      end

      double precision function df2dnz(cnz)
      implicit none
      double precision cnz,cnpar,bz,bmod
      df2dnz=2*cnz
      return
      end

      double precision function df2dnr(cnr)
      implicit none
      double precision cnr,cnpar,br,bmod
      df2dnr=2*cnr
      return
      end

      double precision function df2dr(cm,r)
      implicit none
      double precision cm,r
      df2dr=-cm*cm/r**3
      return
      end
    
      subroutine al12(dnz,dnr,df1dnz,df1dnr,df2dnz,df2dnr,al1,al2)
c     calculates al1 and al2
      implicit none
c-----input
      double precision dnz,dnr,df1dnz,df1dnr,df2dnz,df2dnr
c-----output
      double precision al1,al2
c-----local 
      double precision det,det1,det2

      det=df1dnz*df2dnr-df1dnr*df2dnz
      det1=2.d0*(dnz*df2dnr-dnr*df2dnz)
      det2=2.d0*(-dnz*df1dnr+dnr*df1dnz)

      if (det.eq.0.d0) then
         write(*,*)'output.f in al12 det=0 stop'
         stop
      else
         al1=det1/det
         al2=det2/det
      endif

      return
      end




      double precision function dlagrdr(dr,al1,al2,df1dr,df2dr)
c     calculates derivative d(Lagrange function)/d(deltar) 
      implicit none
c-----input
      double precision dr,al1,al2,df1dr,df2dr 
      dlagrdr=2.d0*dr+al1*df1dr+al2*df2dr
      return
      end
  
  
  
  
      subroutine g123(cnz,cnr,cm,z,r,phi,dnr,dnz,dr,cnpar,cnper,
     .bz,br,bphi,bmod,dbzdr,dbrdr,dbphdr,dbmdr,
     .g123ar)
c     input magnetic field should be done in
c     (r+delr) point
c     calculates
c     g123ar(1)=f1(Nz+delNz,Nr+delNr,r+delr)
c     g123ar(2)=f2(Nz+delNz,Nr+delNr,r+delr)
c     g12ar(3)=2*delr+al1(delNz,delNr,delr)*df1(delNz,delNr,delr)/d(delr)+
c               al2(delNz,delNr,delr)*df2(delNz,delNr,delr)/d(delr)

      implicit none
c-----input 
      double precision cnz,cnr,cm,z,r,phi,dnz,dnr,dr,cnpar,cnper,
     .bz,br,bphi,bmod,dbzdr,dbrdr,dbphdr,dbmdr
c-----external
      double precision dlagrdr,df1dr,df2dr,f1,f2,
     .df1dnz,df1dnr,df2dnz,df2dnr
c-----output
      double precision g123ar(3)
c-----local
      double precision al1,al2,df1drl,df2drl,cn,cnzl,cnrl,rl,
     .df1dnzl,df1dnrl,df2dnzl,df2dnrl
      cnzl=cnz+dnz
      cnrl=cnr+dnr
      rl=r+dr
      df1drl=df1dr(cnzl,cnrl,cm,rl,bphi,bmod,cnpar,
     .dbzdr,dbrdr,dbphdr,dbmdr)

      df2drl=df2dr(cm,rl)

      cn=dsqrt(cnper**2+cnpar**2)

      g123ar(1)=f1(cnzl,cnrl,cm,rl,bz,br,bphi,bmod,cnpar)
      g123ar(2)=f2(cnzl,cnrl,cm,rl,cn)
     
      df1dnzl = df1dnz(bz)
      df1dnrl = df1dnr(br)
      df2dnzl = df2dnz(cnz)
      df2dnrl = df2dnr(cnr)

      call al12(dnz,dnr,df1dnzl,df1dnrl,df2dnzl,df2dnrl,al1,al2)
c      write(*,*)'output al1,al2',al1,al2

      g123ar(3)=dlagrdr(dr,al1,al2,df1drl,df2drl)
c       write(*,*)'output g123 g123ar(3)',g123ar(3)
      return
      end




      subroutine drvg123(cnz,cnr,cm,z,r,phi,dnr,dnz,dr,cnpar,cnper,
     .bz,br,bphi,bmod,dbzdr,dbrdr,dbphdr,dbmdr,
     .derg)
c     calculates derivatives derg from g123
      implicit none
c-----input 
      double precision cnz,cnr,cm,z,r,phi,dnz,dnr,dr,cnpar,cnper,
     .bz,br,bphi,bmod,dbzdr,dbrdr,dbphdr,dbmdr
c-----external
      double precision dlagrdr,df1dr,df2dr,f1,f2,
     .df1dnz,df1dnr,df2dnz,df2dnr
c-----output
      double precision derg(3,3)
c-----local
      double precision al1,al2,df1drl,df2drl, 
     .cnzl,cnrl,rl,dnzl,dnrl,drl,
     .df1dnzl,df1dnrl,df2dnzl,df2dnrl,
     .df1dnzp,df1dnrp,df2dnzp,df2dnrp,
     .df1dnzm,df1dnrm,df2dnzm,df2dnrm,
     .df1drp,df1drm,df2drp,df2drm,
     .step,cnzp,cnzm,cnrp,cnrm,rp,rm,g3p,g3m,
     .dnzp,dnzm,dnrp,dnrm,drp,drm,
     .bzl,brl,bphil,bmodl,
     .dbzdzl,dbzdrl,dbrdzl,dbrdrl,
     .dbphdzl,dbphdrl,dbmdzl,dbmdrl

      cnzl=cnz+dnz
      cnrl=cnr+dnr
      rl=r+dr

      df1drl=df1dr(cnzl,cnrl,cm,rl,bphi,bmod,cnpar,
     .dbzdr,dbrdr,dbphdr,dbmdr)

      df2drl=df2dr(cm,r)

      derg(1,1)=df1dnz(bz)                                  !dg1/ddnz
      derg(1,2)=df1dnr(br)                                  !dg1/ddnr
      derg(1,3)=df1dr(cnzl,cnrl,cm,rl,bphi,bmod,cnpar,      !dg1/ddr
     .dbzdr,dbrdr,dbphdr,dbmdr)

      derg(2,1)=df2dnz(cnzl)                                !dg2/ddnz
      derg(2,2)=df2dnr(cnrl)                                !dg2/ddnz
      derg(2,3)=df2dr(cm,rl)                                !dg2/ddr

      step=1.d-7   
c-----d(g3)/d(dnz)
      dnzp=dnz+step
      cnzp=cnz+dnzp

      df1drp=df1dr(cnzp,cnrl,cm,rl,bphi,bmod,cnpar,
     .dbzdr,dbrdr,dbphdr,dbmdr)

      df2dnzp = df2dnz(cnzp)

      call al12(dnzp,dnr,df1dnzl,df1dnrl,df2dnzp,df2dnrl,al1,al2)

      g3p=dlagrdr(dr,al1,al2,df1drp,df2drl)
      
      dnzm=dnz-step
      cnzm=cnz+dnzm

      df1drm=df1dr(cnzm,cnrl,cm,rl,bphi,bmod,cnpar,
     .dbzdr,dbrdr,dbphdr,dbmdr)

      df2dnzm = df2dnz(cnzm)

      call al12(dnzm,dnr,df1dnzl,df1dnrl,df2dnzm,df2dnrl,al1,al2)

      g3m=dlagrdr(dr,al1,al2,df1drm,df2drl)
      
      derg(3,1)=(g3p-g3m)/(2.d0*step)      !dg3ddnz

c-----d(g3)/d(dnr)
      dnrp=dnr+step
      cnrp=cnr+dnr

      df1drp=df1dr(cnzl,cnrp,cm,rl,bphi,bmod,cnpar,
     .dbzdr,dbrdr,dbphdr,dbmdr)

      df2dnrp = df2dnr(cnrp)

      call al12(dnz,dnrp,df1dnzl,df1dnrl,df2dnzl,df2dnrp,al1,al2)

      g3p=dlagrdr(dr,al1,al2,df1drp,df2drl)
      
      dnrm=dnr-step
      cnrm=cnr+dnrm

      df1drm=df1dr(cnzl,cnrm,cm,rl,bphi,bmod,cnpar,
     .dbzdr,dbrdr,dbphdr,dbmdr)

      df2dnrm = df2dnr(cnrm)

      call al12(dnz,dnrm,df1dnzl,df1dnrl,df2dnzl,df2dnrm,al1,al2)

      g3m=dlagrdr(dr,al1,al2,df1drm,df2drl)
      
      derg(3,2)=(g3p-g3m)/(2.d0*step)         !dg3ddnr

c-----dgr/d(dr)
      drp=dr+step
      rp=r+drp
  
      call bcomp(z,rp,phi,bzl,brl,bphil,bmodl,
     .dbzdzl,dbzdrl,dbrdzl,dbrdrl,
     .dbphdzl,dbphdrl,dbmdzl,dbmdrl)

      df1drp=df1dr(cnzl,cnrl,cm,rp,bphil,bmodl,cnpar,
     .dbzdrl,dbrdrl,dbphdrl,dbmdrl)

      df2drp=df2dr(cm,rp)

      df1dnzp = df1dnz(bzl)
      df1dnrp = df1dnr(brl)

      call al12(dnz,dnr,df1dnzp,df1dnrp,df2dnzl,df2dnrl,al1,al2)

      g3p=dlagrdr(drp,al1,al2,df1drp,df2drp)

      drm=dr-step
      rm=r+drm
  
      call bcomp(z,rm,phi,bzl,brl,bphil,bmodl,
     .dbzdzl,dbzdrl,dbrdzl,dbrdrl,
     .dbphdzl,dbphdrl,dbmdzl,dbmdrl)

      df1drm=df1dr(cnzl,cnrl,cm,rm,bphil,bmodl,cnpar,
     .dbzdrl,dbrdrl,dbphdrl,dbmdrl)

      df2drm=df2dr(cm,rm)

      df1dnzm = df1dnz(bzl)
      df1dnrm = df1dnr(brl)

      call al12(dnz,dnr,df1dnzm,df1dnrm,df2dnzp,df2dnrl,al1,al2)

      g3m=dlagrdr(drm,al1,al2,df1drp,df2drp)

      derg(3,3)=(g3p-g3m)/(2.d0*step)         !dg3ddr

      return
      end



    
      subroutine bcomp(z,r,phi,bzl,brl,bphil,bmodl,
     .dbzdzl,dbzdrl,dbrdzl,dbrdrl,
     .dbphdzl,dbphdrl,dbmdzl,dbmdrl)
c     calculates the magnetic field and its derivatives 
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      double precision b      
      bmodl=b(z,r,phi)
      bzl=bz
      brl=br
      bphil=bphi
      dbzdzl=dbzdz
      dbzdrl=dbzdr
      dbrdzl=dbrdz
      dbrdrl=dbrdr
      dbphdzl=dbphdz
      dbphdrl=dbphdr
      dbmdzl=dbmdz
      dbmdrl=dbmdr

      return
      end




      subroutine set_output ! called by dinit_1ray_xyz
      include 'param.i'
      include 'output.i'
        
      first=.true.
      was_not_ox_conversion=.true. !initialize by set_output/dinit_1ray_xyz

      return
      end
       


 



      subroutine refractive_index_relative_error(z,r,phi,cnz,cnr,cm,
     &iraystop)
c-------------------------------------------------------------------
c     Calculate the relative error of delta_N/N=delta_n_devide_n
c     of the dispersion relation D(cnz,cnr,cm)=0 at the given point
c     (z,z,phi,cnz,cnr,cm)
c     delta_N/N=(D/N)/|gradD|
c     Here
c     D=D(z,z,phi,cnz,cnr,cm)
c     N=|N|=sqrt(Nz**2+Nr**2+(cm/r)**2)
c          N_phi=cm/r
c     |gradD|=sqrt((dD/dN_z)**2+(dD/dN_r)**2+(dD/dN_phi)**2)=
c            =sqrt((dD/dN_z)**2+(dD/dN_r)**2+(dD/dcm)**2*r**2)
c
c     If delta_N/N=(D/N)/|gradD| > toll_hamilt it will put iraystop=1
c     else         iraystop=0
c     Variable toll_hamilt is set in one.i
c------------------------------------------------------------
      implicit none
      include 'param.i'
      include 'one.i'
c-----input
      real*8 z,r,phi, !space coordinates
     *cnz,cnr,cm      !refractive index coordinates
c-----output
      integer iraystop 
        
c-----externals dddrz1
      real*8 b,gamma1,hamilt1

c-----locals
      real*8 u(6),deru(6),
     &grad_d,  ! |gradD in N space|
     &cn,d,     
     &delta_n_devide_n ! delta_N/N=(D/N)/|gradD|

      u(1) = z
      u(2) = r                                                    
      u(3) = phi                                                  
      u(4) = cnz                                                  
      u(5) = cnr                                                  
      u(6) = cm
c      write(*,*)'in  refractive_index_relative_error'
      bmod=b(z,r,phi)
      gam=gamma1(z,r,phi,cnz,cnr,cm)

c      write(*,*)'in  refractive_index_relative_error gam',gam

      call dddrz1(0.d0,u,deru)

c      write(*,*)'in refractive_index_relative_error after dddrz1'

      d=hamilt1(z,r,phi,cnz,cnr,cm)
 
c      write(*,*)'in  refractive_index_relative_error d=',d


      grad_d=dsqrt(deru(1)**2+deru(2)**2+deru(3)**2*r**2)
      cn=dsqrt(cnz**2+cnr**2+(cm/r)**2)
       
      delta_n_devide_n=d/(cn*grad_d)

c      write(*,*)'delta_n_devide_n=',delta_n_devide_n
c      write(*,*)'toll_hamilt=',toll_hamilt

      if( delta_n_devide_n.gt.toll_hamilt) then 
         iraystop=1
      else
         iraystop=0
      endif

      if(iraystop.eq.1) then
         write(*,*)
         write(*,*) '*************************************************'
         write(*,*) 'outpt: in  refractive_index_relative_error'
         write(*,*) ' D/(N|gradD|) > toll_hamilt'
         write(*,*) 'D=',d,'N=',cn,'grad_d=',grad_d
         write(*,*) 'D/(N|gradD|)=delta_n_devide_n',delta_n_devide_n
         write(*,*) 'toll_hamilt=',toll_hamilt
         write(*,*) '*************************************************'
         write(*,*)
      endif
   
      return
      end


