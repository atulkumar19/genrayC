


c        **********************_dinit_mr************************
c        *                        -                           *
c        * this subroutine reads the data from genray.in     *
c        ******************************************************
c
c-------------------------dinit_mr---------------------------------
c        It creates data for multiple ray case. 
c        for spline approximations                       	   
c        density,temperature and Z_effective profiles.             
c        Calculates initial data for all rays.                     
c        It uses input data from genray.in or genray.dat 
c        file for multiple ray case.                               
c        Data from genray.in or genray.dat file were read          
c        previously in genray.f file   
c        It reads or creates the non-maxwellian electron distribution
c
c        it uses the following functions and subroutines           
c        bmod,spldens1                                             

c      	 This program reads data from genray.in file for	   !
c        multiple ray case.                                        !
c        It creates the files for spline approximations 	   !
c        density,temperature and Z_effective profiles.             !
c        it uses the following functions and subroutines           !
c        bmod,spldens1                                             !
c------------------------------------------------------------------
c        output parameters:					   !
c                          ndim-number of the ray-tracing equations!
c                          nray-number of the rays at antenna      !
c------------------------------------------------------------------

      subroutine dinit_mr(ndim,nray) ! converted to xyz

      implicit none

      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'three.i'
      include 'five.i'
      include 'cone.i'
      include 'grill.i'
      include 'rkutta.i'
      include 'six.i'
      include 'write.i'
      include 'writencdf.i'
      include 'onetwo.i'
      include 'output.i'
      !!! include 'scatnper.i'
      !!! include 'emissa.i'
      double precision thetax,phix,xe0
      common /convers/ thetax,phix,xe0

c-----output
      integer
     &ndim,  !number of the ray-tracing equations
     &nray   !number of the rays at antenna
c..............................................................
c     these two arrays are for namelist work only


c      dimension prof2(nbulka,ndensa)
c      dimension prof(nbulka*ndensa)
c     Following are working arrays for the namelist input of density 
c     and temperature,... at nonuniform radii grids
      real*8 prof2_nonuniform(nbulka,ndensa),prof_radii(nbulka,ndensa)
c      integer nj_tab(nbulka) !the number of profile points for each species
c..............................................................
      include 'dinit_nml.i'
c..............................................................
      integer nbulk1,i,j,k,iraystop,i1,j1,ifreq,ii,initial,nray1,icone,
     &imax

      real*8 trnspi,h,v0,w0,hfreq,delt,zefftest,zion,psi,denstot,
     &szini,szi2ni,stini,pressure,den_test,prestest,tem_test,
     &energy,pitch,fdist,dfdx,dfdpitch,dfdp,psi_mag,xi,yi,tetan,wpw2in,
     &zconv,rconv,theta,rhoconv,phiconv,xconv,yconv,cnparopt,
     &h_rho,dens,temp,psi_loc,
     &tpop,
cSAP090311 for test
     &zeff_loc,vflow_loc

c-----externals
      real*8 b,psi_rho,prespsi,densrho,temperho,psif,x,y,tpoprho,
cSAP090311
     & zeffrho,vflowrho
cSAP080711BH080714
c-----genray input file: genray.in or genray.dat
      character*10 genray_in_dat

      iraystop=0
      pi=4*datan(1.d0)
      trnspi=pi/180.d0

c---------------------------------------------------------------
c     read all namelists from input genray.in or genray.dat file
c---------------------------------------------------------------
cSAP080731
c      call read_all_namelists(genray_in_dat,ndim,nray)
c      write(*,*)'in dinit_mr, genray_in_dat=', genray_in_dat
c-----If input file is genray.in, then it will change
c     input data to genray.dat file format 
c      if (genray_in_dat.eq.'genray.in')  then
c         write(*,*)'genray.f before transform_genray_in_to_dat'
c         call transform_genray_in_to_dat
c      endif
c---------------------------------------------------------------------
      if (n_wall.gt.0) then    
c-------calculate poloidal angles [at radian] and small radius
c       of wall and limiter points
c       thetapol_wall(i=1,..,n_wall)
c       thetapol_limiter(i=1,...,n_limiter(j),j=1,max_limiters)
c       rho_wall(i=1,..,n_wall)
c     rho_limiter(i=1,..,n_limiter(j),j=1,max_limiters))

c  
c       These arrays will be in common /fourb/ in file fourb.i
        call wall_limiter_theta_pol_rho

      endif

c------------------------------------------------------------------------
      write(*,*)'Absorption'
c------------------------------------------------------------------------
      write(*,*)'iabsorp=',iabsorp
      if(iabsorp.le.4  .or.
     +   iabsorp.eq.6  .or.  iabsorp.eq.7  .or.  iabsorp.eq.12) then
      ! continue
      else
      stop 'Not setup for this value of iabsorp. Choose from 1-4,6,7,12'
      endif     
c     iabsorp=2 for LH waves
c     iabsorp=3 for FW waves
c-----------------------------------------------------

      write(*,*)'in dinit.f n_relt_harm1',n_relt_harm1
      write(*,*)'in dinit.f n_relt_harm',n_relt_harm
      write(*,*)'in dinit.f n_relt_harma',n_relt_harma
      if (n_relt_harm1.eq.9999)then
        n_relt_harm1=-n_relt_harm
        n_relt_harm2= n_relt_harm
      else
        n_relt_harm2=n_relt_harm1+n_relt_harm-1
      endif
       
      if(n_relt_harm1.lt.n_relt_harm1a) then
         write(*,*)'n_relt_harm1<n_relt_harm1a'
         write(*,*)'it should be n_relt_harm1=>n_relt_harm1a'
         write(*,*)'please change n_relt_harm1a'        
         write(*,*)'in param.i and recompile the code'
         write(*,*)'n_relt_harm1,n_relt_harm1a',
     &               n_relt_harm1,n_relt_harm1a
         stop 'in dinit.f' 
      endif
 
      if(n_relt_harm2.gt.n_relt_harma) then
         write(*,*)'n_relt_harm2>n_relt_harma'
         write(*,*)'it should be n_relt_harm2=<n_relt_harma'
         write(*,*)'please change n_relt_harma'
         write(*,*)'in param.i and recompile the code'
         write(*,*)'n_relt_harm2,n_relt_harm2a',
     &               n_relt_harm2,n_relt_harm2a
         stop 'in dinit.f' 
      endif
      write(*,*)'dinit.f 1 nbulk=',nbulk 
      do i=1,nbulk 
        if ((i_salphal(i).ne.1).and.(i_salphal(i).ne.0)) then
          write(*,*)'(i_salphal(i).ne.1).or.(i_salphal(i).ne.0)'
          write(*,*)'It should i_salphal(i) =0 or =1'
          write(*,*)'i=',i,'i_salphal(i)',i_salphal(i)
          write(*,*)'Please chagne i_saplhal(i)'
          write(*,*)'in input genray.in or genray.dat file'
          stop 'in dinit.f i_salphal problem'
        endif
      enddo 
           
      write(*,*)'Plasma parameters'
      write(*,*)'izeff=',izeff
      do i=1,nbulk
         write(*,*)'i,temp_scale(i),den_scale(i)',
     .   i,temp_scale(i),den_scale(i)
      enddo 

      do i=1,nbulk
         te0(i)=ate0(i)
         teb(i)=ateb(i)
      enddo
c-----------------------------------------------------
c      /species/
c   plasma component charges charge(i)=mod(charge(i)/charge_electron)
c----------------------------------------------------

c-----------------------------------------------------
c  plasma components mass dmas(i)=Mass(i)/Mass_electron
c-----------------------------------------------------
      do i=1,nbulk
         write(*,*)'i, dmas(i)',i,dmas(i)
      enddo
c--------------------------------------------------------------

c     calculation of nbulk1
      if(((izeff.eq.0).or.(izeff.eq.2)).or.(izeff.eq.3)) then
c        izeff=0, zeff will be calculated using the given ions;
c                 the electron density will be calculated using ion's densities;
c             =1  ion's densities(i) i=nbulk and i=nbulk-1 will be calculated  using
c                 Zeff, the electon density and ion's densities(i), i=2,nbulk-1;
c        izeff=2, zeff will not coincide with the plasma components
c             =3  it uses eqdsk pres (pressure) and ions densities_i
c                 for i=2,... nbulk
c                 Let temperature T_E=T_i
c                 pres=dens1(k,1)*temp1(k,1)+
c                      Sum(i=2,nbulk)(dens1(k,i)*temp1(k,i)
c                 In this case we will calculate Zeff(rho),
c                 dens_electron(rho) and T_e(rho)=T_i(rho)
c             =4  it uses eqdsk pres (pressure), zeff,ions densities
c                 for i=2,... nbulk-2 (nbulk>2) and T_e(rho),T_i(rho)
c                 pres=dens1(k,1)*temp1(k,1)+
c                      Sum(i=2,nbulk)(dens1(k,i)*temp1(k,i)
c                 In this case we will calculate dens_electron(rho) and
c                 ion densities for i=nbulk and i=nbulk-1
         nbulk1=nbulk         
      else
c        (izeff=1 or izeff=4), zeff is given, the ions component will be calculated
         if (nbulk.le.2) nbulk1=nbulk
         if (nbulk.eq.2) then
	    write(*,*)'nbulk=2, Zeff must be equal charge(2)'
	    write(*,*)'Please check it or use the option izeff=0'
cSAP090801
c	    stop
	 endif
         if (nbulk.gt.2) nbulk1=nbulk-2
      endif !izeff
      write(*,*)'nbulk1=',nbulk1
c------------------------------------------------------------------
      h=1.d0/(ndens-1)
      do i=1,ndens
        rhom(i)=h*(i-1)
      enddo
c------------------------------------------------------------------
c------------------------------------------------------------------
      if(idens.eq.0) then        
         write(*,*)'Analytical radial profiles'
c-dense_(i)=(dense0(i)-denseb(i))*(1-rho**rn1de(i))**rn2de(i)+denseb(i)
c------------------------------------------------------------------\     
	 do i=1,nbulk1
	   write(*,*)'i, dense0(i)',i,dense0(i)
	 enddo

	 do i=1,nbulk1
	   write(*,*)'i, denseb(i)',i,denseb(i)
	 enddo

	 do i=1,nbulk1
	   write(*,*)'i, rn1de(i)',i,rn1de(i)
	 enddo

	 do i=1,nbulk1
	   write(*,*)'i, rn2de(i)',i,rn2de(i)
	 enddo

	 do i=1,nbulk1
	   write(*,*)'i, dense0(i)',i,dense0(i)
	 enddo

c--------creation the array dens1(ndensa,nbulka)
cSm080118
         do i=1,nbulk1
            do k=1,ndens
               rho=rhom(k)
	       dens1(k,i)=(dense0(i)-denseb(i))*
     1	                  (1-rho**rn1de(i))**rn2de(i)+denseb(i)
            enddo
         enddo

	 do i=nbulk1,1,-1
	    do k=1,ndens
	       rho=rhom(k)

	       if (((izeff.eq.0).or.(izeff.eq.3)).and.(i.eq.1)) then
	          dens1(k,1)=0.d0
		  do j=2,nbulk
		    dens1(k,1)=dens1(k,1)+charge(j)*dens1(k,j)
		  enddo
	       else
	          dens1(k,i)=(dense0(i)-denseb(i))*
     1	                     (1-rho**rn1de(i))**rn2de(i)+denseb(i)
cSm070426
c-----------------multiply density profiles by den_scale
                  dens1(k,i)=dens1(k,i)*den_scale(i)
                  write(*,*)'dinit i,k,dens1(k,i)',i,k,dens1(k,i)
                write(*,*)'dense0(i),denseb(i),rho,rn1de(i),rn2de(i)',
     *          dense0(i),denseb(i),rho,rn1de(i),rn2de(i)
      
               endif
            enddo
         enddo
	 write(*,*)'end of analytical density profiles input'

c------------------------------------------------------------------         
c         /tpopprof/
c         /vflprof/
c--------------------------------------------------------------------       
c        creation the arrays for analytical profoliles
c        tpop1(ndensa,nbulka), vflow1(ndensa,nbulka)
	 do i=1,nbulk
	    do k=1,ndens
	       rho=rhom(k)	
	       tpop1(k,i)=(tp0(i)-tpb(i))*
     1                    (1-rho**rn1tp(i))**rn2tp(i)+tpb(i)

               vflow1(k,i)=(vfl0(i)-vflb(i))*
     1                    (1-rho**rn1vfl(i))**rn2vfl(i)+vflb(i)
	    enddo
	 enddo

c---------------------------------------------------------------------
c        /zprof/          
	 if (izeff.eq.3) goto 10
c        /tprof/
c----------------------------------------------------------
cTemperature
ctempe_(i)=(te0(i)-teb(i))*(1-rho**rn1te(i))**rn2te(i)+teb(i)
c-----------------------------------------------------------
	 do i=1,nbulk
           write(*,*)'i, te0(i)',i,te0(i)	
	   write(*,*)'i, teb(i)',i,teb(i)
	 enddo

	 do i=1,nbulk
	   write(*,*)'i, rn1te(i)',i,rn1te(i)
	 enddo

	 do i=1,nbulk
	   write(*,*)'i, rn2te(i)',i,rn2te(i)
	 enddo

         do i=1,nbulk
          write(*,*)'i,tp0(i),tpb(i)',i,tp0(i),tpb(i)
          write(*,*)'rn1tp(i),rn2tp(i)',rn1tp(i),rn2tp(i)
         enddo

c------- creation of array temp1(ndensa,nbulka)
	 do i=1,nbulk
	    do k=1,ndens
	       rho=rhom(k)
	       temp1(k,i)=(te0(i)-teb(i))*
     1                    (1-rho**rn1te(i))**rn2te(i)+teb(i)
cSm070426
c--------------multiply temperature profiles by temp_scale
               temp1(k,i)=temp1(k,i)*temp_scale(i)
	    enddo
	 enddo
10	 continue         


         write(*,*)'zeff0,zeffb,rn1zeff,rn2zeff'
         write(*,*)zeff0,zeffb,rn1zeff,rn2zeff
	 if(((izeff.eq.1).or.(izeff.eq.2)).or.(izeff.eq.4)) then
c           the given analytical Zeff profile
c-------------------------------------------
c           zeff=(zeff0-zeffb)*(1-rho**rn1zeff)**rn2zeff+zeffb
c-------------------------------------------
c           the creation of array zeff1(ndens)
	    do k=1,ndens
	       rho=rhom(k)
               zeff1(k)=(zeff0-zeffb)*(1-rho**rn1zeff)**rn2zeff+zeffb
	    enddo
	 endif

	 if((izeff.eq.0).or.(izeff.eq.3)) then
	    write(*,*)'izeff=1 zeff will be calculated using the given
     1	    ions densities'
	 endif
      endif ! idens analytical
c--------------------------------------------------------------------

c--------------------------------------------------------------
      if (idens.eq.1) then
c        spline approximation of the density, temperature, zeff,
c        tpop=T_perp/T_parallel, vflow
c        radial profiles
c        input of the arrays on the radial mesh	from dtzprof.dat

        write(*,*)'dinit_mr: Spline radial profiles '
        write(*,*)'nbulka,ndensa,nbulk,ndens',nbulka,ndensa,nbulk,ndens
        
        write(*,*)'dinit_mr: before mult by den_scale'
        do i=1,nbulk
          write(*,*)'i',i,'dens1(k,i)'
          write(*,*)(dens1(k,i),k=1,ndens)
        enddo
c---------multiply density profiles by den_scale(i)
cBh080102          do i=i1,nbulk          
        do i=1,nbulk
        do k=1,ndens
           dens1(k,i)=den_scale(i)*dens1(k,i)
        enddo
        enddo
        
        write(*,*)'dinit_mr: nbulk1',nbulk1,'ndens',ndens
        if ((izeff.eq.0).or.(izeff.eq.3)) then
           i1=2
        else
           i1=1
        endif

        write(*,*)'dinit_mr: after mult by den_scale' 
        do i=i1,nbulk1
           write(*,*)'i',i,'dens1(k,i)'
           write(*,*)(dens1(k,i),k=1,ndens)
        enddo
         
        !c---------------------------------------
        ! calculation of the electron density from
        ! the charge neutrality 
        !c---------------------------------------
        if ((izeff.eq.0).or.(izeff.eq.3)) then
           do k=1,ndens
              dens1(k,1)=0.d0
           do j=2,nbulk
              dens1(k,1)=dens1(k,1)+charge(j)*dens1(k,j)
           enddo
           enddo
         endif

         write(*,*)'dinit dens1(k,1)'
         write(*,*)(dens1(k,1),k=1,ndens)

c--------multiply temperature profiles by temp_scale
         do i=1,nbulk
         do k=1,ndens
c               write(*,*)'i,k,temp1(k,i),temp_scale(i)',
c     &                    i,k,temp1(k,i),temp_scale(i)
               temp1(k,i)=temp1(k,i)*temp_scale(i)
c               write(*,*)'temp1(k,i)',temp1(k,i)
         enddo
         enddo

         do i=1,nbulk
            write(*,*)'i',i,'temp1(k,i)'
            write(*,*) (temp1(k,i),k=1,ndens)
         enddo

c---------tpoptab
         do i=1,nbulk
            write(*,*)'i',i,'tpop1(k,i)'
            write(*,*) (tpop1(k,i),k=1,ndens)
         enddo

c--------vflowtab
         do i=1,nbulk
            write(*,*)'i',i,'vflow1(k,i)'
            write(*,*) (vflow1(k,i),k=1,ndens)
         enddo

         write(*,*)'dinit_mr: izeff=',izeff
         if(((izeff.eq.1).or.(izeff.eq.2)).or.(izeff.eq.4)) then
c           the given Zeff profile is in the table form
            write(*,*)'dinit_mr: zeff1(k)'
         else
            write(*,*)'dinit_mr: uniform zeff1',zeff1
         endif !izeff              
      endif ! idens=1
c-----------------------------------------------------------


 20   continue 
      write(*,*)'dinit_mr after    20   continue '
c------------------------------------------------------------

c---------------------------------------------------------
c     read the data for for EC cone vertex coordinates calculations
c     This case is used for the optimal OX mode conversion.
      
c      /ox/


      if(i_ox.eq.1) then
        istart=3
        prmt(3)=-prmt3 !to create the negative time
        i_vgr_ini=+1
        ireflm=1   
        write(*,*)'dinit.f i_ox.eq.1 prmt(3)',prmt(3)
      endif

      if(((i_ox.ne.0).and.(i_ox.ne.1)).and.(i_ox.ne.2)) then
         write(*,*)'i_ox can  =0 or =1 or =2'
         write(*,*)'in namelist /ox/ i_ox=',i_ox
         write(*,*)'please change i_ox in input file'
         stop 'in dinit /ox/'
      endif  
         
    

c------------------------------------------------------------
c       no emission calculations, set one frequency
        nfreq=1
c-----------------------------------------------
c     for the emission multi frequency case
c     calculate the electron gyro-frequency freqncy0 at the plasma center
c---------------------------------------------
      bmod=b(zma,rma,0.d0) !TL
      freqncy0=28.0*b0*bmod !GHZ

      if (nfreq.eq.1) then
c--------no emission or only one frequency in the emission calculations   
         v0=806.2/(frqncy*frqncy)
         w0=28.0*b0/frqncy 
c-----------------------------------------------------
c        set the arrays for v and w
         v(1)=v0
         w(1)=w0      
         do i=2,nbulk
            v(i)=v0*charge(i)**2/dmas(i)
            w(i)=w0*charge(i)/dmas(i)
            write(*,*)'dinit: i charge(i),dmas(i),v(i),w(i)'
     +      ,i,charge(i),dmas(i),v(i),w(i)
         enddo
      else
      endif
     
c-------------------------------------------------------------
      if ((izeff.eq.0).or.(izeff.eq.3))then
c---------------------------------------------------------
c        calculation of the table for the radial profile zeff1(ndens)
         call zeffcalc
c        zeff1(ndens) is in common six.i
	 do i=1,ndens
	   write(*,*)'i',i,'zeff1(i)',zeff1(i)
         enddo
      endif !izeff=0
c---------------------------------------------------------
      if(izeff.eq.1) then
c---------------------------------------------------------
c        calculation of the table for the ion densities profiles
c---------------------------------------------------------
         if (nbulk.lt.3) then
            write(*,*)'nbulk.lt.3, Zeff must be equal charge(2),
     1	    control it and use the option izeff=0'
            stop
         else
c           nbulk.ge.3
c           calculation of the tables for the radial profile
c           dens1(ndens,nbulk) and dens1(ndens,nbulk-1)
            if( charge(nbulk).eq.charge(nbulk-1)) then
              write(*,*)'Warning in dinit: nbulk(.ge.3)=',nbulk
              write(*,*)'in dinit: charge(nbulk)=charge(nbulk-1)'
              write(*,*)'it is impossible to find the ions densities'
              write(*,*)'change charge(nulk) or charge(nbulk-1)'
              write(*,*)'it should be charge(nulk)>charge(nbulk-1)'
              write(*,*)'or use the option izeff=0'
              stop
            endif

            call denscalc
c------------------------------------------------
c for test
            do i1=1,nbulk
               write(*,*)'in dinit after call denscalc i1=',i1
               do j1=1,ndens
	          write(*,*)'j1=',j1,'dens1(j1,i1)',dens1(j1,i1)
               enddo
	    enddo
	    do j1=1,ndens
	       zefftest=0.d0
	       zion=0.d0
	       do i1=2,nbulk
	          if(dens1(j1,1).ne.0.d0) then
	             zefftest=zefftest+(dens1(j1,i1)/dens1(j1,1))*
     1                        charge(i1)*charge(i1)
		     zion=zion+charge(i1)*dens1(j1,i1)/dens1(j1,1)
		  else
		     write(*,*)'dens1(j1,1)=0'
		     zefftest=zefftest+1.d0*charge(i1)*charge(i1)
		  endif
	       enddo
	       write(*,*)'j1',j1,'zefftest',zefftest,'zion',zion
	    enddo
c end test
c------------------------------------------------

	 endif! nbulk
      endif ! izeff=1
      if (izeff.eq.3) then
	 do i=1,nbulk
	    do k=1,ndens
	       rho=rhom(k)
	       psi=psi_rho(rho)
	       denstot=dens1(k,1)
	       do ii=2,nbulk
		  denstot=denstot+dens1(k,ii)
	       enddo
	       temp1(k,i)=prespsi(psi)/denstot/(1.6d3)
	       write(*,*)'in dinit izeff=3,k,rho',k,rho
	       write(*,*)'prespsi(psi),denstot',prespsi(psi),denstot
	       write(*,*)'in dinit i,temp1(k,i)',i,temp1(k,i)
	       write(*,*)'in dinit i,dens1(k,i)',i,dens1(k,i)
	       write(*,*)'in dinit i,dens1(k,1)',i,dens1(k,1)

	    enddo
	 enddo
      endif !izeff=3

      if(izeff.eq.4) then

cSAP090801 
c        if(nbulk.lt.3) then
c	  write(*,*)'in dinit: izeff=4 nbulk=',nbulk
c	  write(*,*)'in dinit: it should be nbulk>2, use another '
c	  write(*,*)'in dinit: izeff option '
c	  stop
c	else

c         nbulk.ge.3
c         calculation of the tables for the radial profiles
c         dens1(ndens,nbulk),dens1(ndens,nbulk-1) and
c         dens1(ndens,1)
	  call denscalp
c test izeff=4	beg
cSAP090801
        if (nbulk.ge.3) then
	  do j=1,ndens
	    szini=0.d0	   !sum(i=2,nbulk){charge(i)*dens1(j,i)}
	    szi2ni=0.d0	   !sum(i=2,nbulk){charge(i)**2*dens1(j,i)}
	    stini=0.d0	   !sum(i=2,nbulk){temp1(j,i)*dens1(j,i})
	    rho=rhom(j)
	    psi=psi_rho(rho)
	    pressure=prespsi(psi)/1.6d3
	    do i=2,nbulk
	       szini=szini+charge(i)*dens1(j,i)
	       szi2ni=szi2ni+charge(i)*charge(i)*dens1(j,i)
	       stini=stini+temp1(j,i)*dens1(j,i)
	    enddo
	    prestest=temp1(j,1)*dens1(j,1)+stini
	    zefftest=szi2ni/dens1(j,1)
	    write(*,*)'in dinit izeff=4,j,rho',j,rho
	    write(*,*)'pressure,prestest',pressure,prestest
	    write(*,*)'zeff,zefftest',zeff1(j),zefftest
	    write(*,*)'dens1(j,1),szini',dens1(j,1),szini
	  enddo !j
c test izeff=4	end
	endif !nbulk.ge.3
      endif !izeff=4

 30   continue ! if (partner.eq. ....) goto 20

c     creation of the density,temperature,zeff and
c     tpop, vflow
c     spline coefficients
      call spldens1
      write(*,*)'dinit.f after spldens1'
c-----test printing density,temperature,zeff
      do j=1,ndens
         write(*,*)'j,rhom(j)',j,rhom(j)
         do i=1,nbulk
           den_test=densrho(rhom(j),i)
           tem_test=temperho(rhom(j),i)
           write(*,*)'i,dens(i,rhom(j)),temp(i,rhom(j))',
     &     i,den_test,tem_test
         enddo
      enddo

      do i=1,ncone ! from genray.in
         phist(i)=phist(i)*trnspi ! degree to radians
         betast(i)=betast(i)*trnspi
         alfast(i)=alfast(i)*trnspi
      enddo
c--------------------------------------------------------------
c  normalization of the start parameters
      do i=1,ncone
         zst(i)=zst(i)/r0x ! from genray.in
         rst(i)=rst(i)/r0x ! from genray.in
         xst(i)=rst(i)*dcos(phist(i))
         yst(i)=rst(i)*dsin(phist(i))
      enddo
c  normaliszation 'diskdisk' parameters  
      d_disk=d_disk/r0x
      rho_launching_disk=rho_launching_disk/r0x
      rho_focus_disk=rho_focus_disk/r0x
      sigma_launching_disk=sigma_launching_disk/r0x
c-----------------------------------------------------------
      psi_mag=psif(zma,rma)
      write(*,*)'psi on the magnetic axis psi_mag=',psi_mag
ctest
      if (nfreq.eq.1) then
c--------no emission or only one frequency in the emission calculations   
         xi=x(zma,rma,0.d0,1)
         yi=y(zma,rma,0.d0,1)
         write(*,*)'dinit.f at magnetic axis Xe,Ye ',xi,yi
c---------------------------------------------------------------
c     the creation the data for the contours X_e=const,Y_e=const
c     B_tot,B_tor, B_pol on the plate (rho,theta) 
c     These data will have the form of the tables in 
c     xybrhoth.dat: rho(i),theta(j),xe,ye,(xe+ye*ye),bmod,bphi,
c     *              dsqrt(bz**2+br**2)
c      nrhomap=30
c      nthetmap=100
c      call mapxyb(nrhomap,nthetmap)
c      write(*,*)'mapxyb was created'
c      stop
      endif
c---------------------------------------------------------------
c     if ray start point is outside the plasma (ECR -case)
c     then:
c     determination of arrays :1)for ray coordinates on the ECR cone
c     alphaj(nray),betaj(nray)-spherical and cylindrical angles
c                              2)for wave power	angle distribution
c     powj(nray) -power flowing in the ray chanel at antenna
c     with normalized  Sum(i=1,nray)delpw0(i)=1
c---------------------------------------------------------------
      if (istart.eq.1) then
c---------EC waves---------------------------------
         write(*,*)'in dinit istart=',istart
         write(*,*)'raypatt=',raypatt
         write(*,*)'ncone=',ncone

         nray1=1                ! Counter for position in ray arrays
         do icone=1,ncone
            if (raypatt.ne.'toray') then
            else
               tetan=betast(icone)
            endif
            write(*,*)'alpha1,na1,na2,alpha2,phist,alfast,tetan',
     1           alpha1(icone),na1,na2,alpha2(icone),phist(icone),
     1           alfast(icone),tetan
c            if (raypatt.ne.'toray') then
            if (raypatt.eq.'genray') then
               alpha1(icone)=alpha1(icone)*trnspi
               alpha2(icone)=alpha2(icone)*trnspi
               tetan=0.5d0*pi-betast(icone) !Polar angle (radians)
               write(*,*)'in dinit_mr before cone_ec'
               !-----------
               call cone_ec(alpha1(icone),na1,na2,alpha2(icone),
     1              phist(icone),alfast(icone),tetan,powtot(icone),nray,
     1              alphaj(nray1),betaj(nray1),powj(nray1)) !->out
               !-----------
               do i=nray1,(nray1-1)+nray
                  zstj(i)=zst(icone)
                  rstj(i)=rst(icone)
                  phistj(i)=phist(icone)
                  xstj(i)=rstj(i)*dcos(phistj(i))
                  ystj(i)=rstj(i)*dsin(phistj(i))
               enddo
               write(*,*)'in dinit_mr after cone_ec:'
               do i=nray1,(nray1-1)+nray
                  write(*,*)' i,powj(i)',i,powj(i)
                  write(*,*)' i,betaj(i),alphaj(i)',i,betaj(i),alphaj(i)
               enddo
            endif
            if(raypatt.eq.'toray') then
               alfast(icone)=alfast(icone)/trnspi
               tetan=betast(icone)/trnspi
               write(*,*)'dinit bef raypat: tetan,alfast,gzone,nray_in',
     1              tetan,alfast(icone),gzone,nray_in
               call raypat(tetan,alfast(icone),alpha1(icone),cr,nray_in,
     1              gzone,mray,betaj(nray1),alphaj(nray1))
               
               write(*,*)'dinit after raypat: gzone,nray_in',
     1              gzone,nray_in
               nray=nray_in
               do i=nray1,(nray1-1)+nray
                  zstj(i)=zst(icone)
                  rstj(i)=rst(icone)
                  phistj(i)=phist(icone)
                  xstj(i)=rstj(i)*dcos(phistj(i))
                  ystj(i)=rstj(i)*dsin(phistj(i))
               enddo
               
               do i=nray1,(nray1-1)+nray
                  write(*,*)'Raypat ray starting angles (degrees):'
                  write(*,*)' i,betaj(i),alphaj(i)',i,betaj(i),alphaj(i)
               enddo
               do i=nray1,(nray1-1)+nray
                  powj(i)=(powtot(icone)/nray)*1.e13
                  alphaj(i)=alphaj(i)*trnspi
                  betaj(i)=(90.-betaj(i))*trnspi
               enddo
               write(*,*)'in dinit_mr after raypat powj(i)'
               do i=nray1,(nray1-1)+nray
                  write(*,*)' i,powj(i)',i,powj(i)
                  write(*,*)' i,betaj(i),alphaj(i)',i,betaj(i),alphaj(i)
               enddo
            endif               !(On raypatt)
             
            if (raypatt.eq.'diskdisk') then
               call disk_to_disk_rays_initial_launching_data(nray)
            endif !diskdisk
c           stop 'dinit.f after disk_to_disk_rays_initial_launching_dat'

            if (raypatt.eq.'diskbeam') then

               write(*,*)'dinit.f before'
               write(*,*)'disk_beam_rays_initial_launching_data'

               call disk_beam_rays_initial_launching_data(nray)

               write(*,*)'dinit.f after'
               write(*,*)'disk_beam_rays_initial_launching_data'

            endif !diskdisk
            powtott=powtott+powtot(icone)*1.e13   !(erg/sec)
            nray1=nray1+nray
         enddo                  !(On icone)
         nray=nray1-1
        write(*,*)'dinit, after EC starting conditions, total nray=',
     1        nray
      endif                     !(On istart.eq.1, EC)

      if (istart.eq.2) then
c---------LH or WF----------
	  write(*,*)'before grill_lh '
	  write(*,*)'ngrilla,ngrill,rma,zma',ngrilla,ngrill,rma,zma
c--------------------------------------------------------------
          call grill_lh(rhopsi0,ngrilla,ngrill,thgrill,phigrill,
     1    height,nthin,nthinmax,
     1    anmin,anmax,nnkpar,powers,powtott,
     1    antormin,antormax,nnktor,anpolmin,anpolmax,nnkpol,
     1    n_theta_pol,
     1    rma,zma,psimag,
     1    fpwth,wtheta,wtheta_edge,anzin,anzin_edge,pwcpl,nnkprmax,
     1    anztorin,anzpolin,pwcpl_tp,nnktormax,nnkpolmax,
     1    anztorin_edge,anzpolin_edge,pwcpl_t_ar,pwcpl_p_ar,
     1    nray,arzu0,arru0,arphiu0,arntheta,arnphi,powinilh,
     1    nraymaxl,wdnpar0,igrillpw,igrilltw,i_n_poloidal,
     1    ilaunch,r0launch,phi0launch,z0launch,
     1    i_grill_pol_mesh,i_grill_npar_ntor_npol_mesh)

	   write(*,*)'in dinitmr after grill nray=',nray

          do iray=1,nray         
             write(*,*)'iray,arzu0(iray),arru0(iray),arphiu0(iray)',
     &                  iray,arzu0(iray),arru0(iray),arphiu0(iray)
          enddo

      endif !LH or FW

      if (istart.eq.3) then
c---------ECR O_X_EBW mode conversion case----------
c         It uses the wave input data from the grill form
c         It sets i_n_poloidal=1
c	  It calculates: 
c	  1) the space point M_conv=(rhoconv,theta)=(zconv,rconv)
c            at which x_e=wpw2in=1. The value of the poloidal
c            angle theta (degree) is given bellow
c            if theta=0  point is on the outer part of the magnetic surface
c            if theta=90 point is on the top of the magnetic surface
c         2) the optimal value N_parallel_optimal=cnparopt
c            for O_X mode conversion
c         3) gives
c            anmin(1)=cnparopt-0.01d0 
c            anmax(1)=cnparopt+0.01d0 
c         4) call grill_lh to calculate the initial data
c         5) gives the coordinates of the initial O_X mode conversion poin 
c            arzu0(1)=zconv
c            arru0(1)=rconv
c-------------------------------------------------------------
          i_n_poloidal=1
c         these lines are to create the start point inside the 
c         plasma for the given point xe=wpw2in
cSAP050510	  wpw2in=1.d0
          wpw2in=1.002d0
c          theta=0.d0   !poloidal angle  (degree)
c          theta=-30.d0
          theta=thgrill(1)
          write(*,*)'dinit before owconvr theta=',theta
          phix=0.d0 ! YuP Added. 
          call owconvr(theta,wpw2in,rhoconv,zconv,rconv)
    	    write(*,*)'dinit rhoconv,zconv,rconv',rhoconv,zconv,rconv
          rhopsi0(1)=rhoconv
          phiconv=phix ! YuP was 0.d0
          bmod=b(zconv,rconv,phiconv)
          xconv=x(zconv,rconv,0.d0,1)
          yconv=y(zconv,rconv,0.d0,1)
c---------calculation of the optimal value N_parallel_optimal
c         for O_X mode conversion
          cnparopt=dsqrt(yconv/(1.d0+yconv))
    	    write(*,*)'dinit xconv,yconv,cnparopt',xconv,yconv,cnparopt
c         write(*,*)'dinit old value of rhopsi0(1)',rhopsi0(1) 
          rhopsi0(1)=rhoconv
          write(*,*)'dinit new rhopsi0(1)',rhopsi0(1) 
          write(*,*)'dinit old anmin(1),anmax(1)',anmin(1),anmax(1)
          anmin(1)=cnparopt-0.01d0 
          anmax(1)=cnparopt+0.01d0 
          write(*,*)'dinit new anmin(1),anmax(1)',anmin(1),anmax(1)
c---------------------------------------------------------------
          write(*,*)'dinit before grill_lh ngrilla,ngrill',
     1    ngrilla,ngrill
          call grill_lh(rhopsi0,ngrilla,ngrill,thgrill,phigrill,
     1    height,nthin,nthinmax,
     1    anmin,anmax,nnkpar,powers,powtott,
     1    antormin,antormax,nnktor,anpolmin,anpolmax,nnkpol,
     1    n_theta_pol,
     1    rma,zma,psimag,
     1    fpwth,wtheta,wtheta_edge,anzin,anzin_edge,pwcpl,nnkprmax,
     1    anztorin,anzpolin,pwcpl_tp,nnktormax,nnkpolmax,
     1    anztorin_edge,anzpolin_edge,pwcpl_t_ar,pwcpl_p_ar,
     1    nray,arzu0,arru0,arphiu0,arntheta,arnphi,powinilh,
     1    nraymaxl,wdnpar0,igrillpw,igrilltw,i_n_poloidal,
     1    ilaunch,r0launch,phi0launch,z0launch,
     1    i_grill_pol_mesh,i_grill_npar_ntor_npol_mesh)


c for ecr internal case 
          arzu0(1)=zconv
          arru0(1)=rconv
c--end ecr internal case
     
	   write(*,*)'dinit_mr: after grill'
      endif ! O_X_EBW  mode conversion case

ctest density profile
      imax=101
      h_rho=1.d0/dfloat(imax-1)
      open(88,file='dens.bin',form='unformatted')
 200  format(2(1pe10.3))

      do i=0,imax-1
         rho=h_rho*i
         dens=densrho(rho,1)
         temp=temperho(rho,1)
         tpop=tpoprho(rho,1)
         zeff_loc=zeffrho(rho)
         vflow_loc=vflowrho(rho,1) 
         psi_loc=psi_rho(rho)
         write(88) real(rho),real(dens),real(temp)
      enddo
      write(88) 
      close(88)
c-----allocate pointers at writencdf.i and write_i
      write(*,*)'dinit_mr:  before  ainalloc_writencdf nray',nray
      call ainalloc_writencdf_i(nray)
      call ainalloc_write_i(nray)

      return ! dinit_mr
      end











c        **********************dinit_1ray**********************
c        * this subroutine                                    *
c        * 1)if ray was launched outside the plasma           *
c        *   (istart=1, ECR wave)			      *
c        *   it detemines the point where the ray inter-    *
c        *   sects the plasma boundary(zu0,ru0,phiu0),and the *
c        *   tangent components of the refractive index	      *
c        *   cnteta,cnphi} .				      *
c        *   \hat{theta}=grad(psi) x \hat{phi}/abs(grad(psi)).*
c        *   Here \hat{theta} is the poloidal unit vector.    *
c        *   Here x is the vector product.		      *
c        *   If istart.ne.1 it launch the ray                 *
c        *   in the given point inside the plasma .	      *
c        * 2)	Then it calculates the parallel to magnetic   *
c        *   field component of the refractive index  cnpar.  *
c        *      Then it determinates the normal   to   magne- *
c        *   tic surface component of refractive index cnrho. *
c        *   It is directed inside the plasma		      *
c        *      Then it determinates the initial components   *
c        *   of the refractive index:			      *
c        *              cnz,cnr,cm=cnphi*r		      *
c	 *   Then it calculate and prints the value of        *
c        *   dispersion function d=eps                        *
c        *   The result is  the initial condition  for the    *
c        *   ray equations				      *
c        ******************************************************
c
c------------------------------------------------------------------
c								   !
c        input parameters					   !
c        zst,rst,phist                                            !
c	 EC wave:						   !
c        alfast- toroidal angle(radian) of the ray at antenna      !
c        betast- angle(radian) between the  horizontal             !
c                          plane and the ray at antenna            !
c        LH and FW waves:					   !
c        cnteta,cnphi						   !
c        output parameters					   !
c        u(1)=zu0,u(2)=ru0,u(3)=phiu0,u(4)=cnz,u(5)=cnr,u(6)=cm	   !
c        iraystop- index to stop the ray determination(=1)	   !
c       	   or make the ray determination (=0)		   !
c------------------------------------------------------------------
c        it uses the following functions and subroutines           !
c        ias1r,bmod,y,x,gamma1,s,abc,hamilt1,plasmaray,ninit_ec     !
c        cinit                                                     !
c------------------------------------------------------------------
      subroutine dinit_1ray(zst,rst,phist,alfast,betast,
     1                      cnteta,cnphi,u,iraystop)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'three.i'
      include 'five.i'
      include 'grill.i'
      include 'write.i'
      include 'rkutta.i'

      dimension u(6),deru(6)

cfor test_to plot D(ReN_perp,ImN_perp)
      integer n_param
      parameter (n_param=3)
      character*(15) name_param(n_param) !names of the input parameters
      real param(n_param)               !values of the input parameters     
cendtest

c-----for hot plasma roots
      double precision 
     &N_perp_root_ar(n_hot_roots_a)     !hot plasma roots
c---------------------------------------------------------------
c------for  rho_ini_LHFW
      integer  i_rho_ini_LHWH_found


      iraystop=0
      pi=4*datan(1.d0)
      trnspi=pi/180.d0
c---------------------------------------------------------------
c     if ray start point is outside the plasma (ECR -case)
c     then:
c          determine:    1) where the vacuum ray
c                           intersects the plasma boundary
c                        2) cnteta,cnphi-toroidal and poloidal
c                           components of the vacuum
c                           refractive index
c---------------------------------------------------------------
      call set_output ! initialize output.i
                      ! for each new ray
      
c      if (i_ox.eq.1) then
          call set_oxb          ! initialize oxb.i
                                ! for each new ray
c      endif
      if (istart.eq.1) then
c--------EC wave
         call plasmray(zst,rst,phist,alfast,betast,
     1                 zu0,ru0,phiu0,iraystop, raypatt)
         if (iraystop.eq.1) then
	      return
         end if
c	 -------------------------------
c        the shift of the initial point inside the plasma from the boundary
	   z=zu0
	   r=ru0
         call edgcor(z,r,zu0,ru0)
c        end of the shift
c	 -------------------------------
c        nx,ny,nz in start point
         cnzst=dsin(betast)
         cnxst=dcos(betast)*dcos(alfast+phist)
         cnyst=dcos(betast)*dsin(alfast+phist)
c----------------------------------------------------------------
	   bmod=b(zu0,ru0,phiu0)
c----------------------------------------------------------------
c        ninit_ec creates the tangent to magnetic surface
c        components  of the refractive index cnteta,cnp                 
c        in the initial point (zu0,ru0,phiu0) for ECR wave
c-----------------------------------------------------------------
         call ninit_ec(zu0,ru0,phiu0,cnxst,cnyst,cnzst,cnteta,cnphi)
      endif

      if ((istart.eq.2).or.(istart.eq.3)) then
c--------LH and FW wave, OX-conversion pt.
         zu0=zst
         ru0=rst
         phiu0=phist
         bmod=b(zu0,ru0,phiu0)
         write(*,*)' in dinit_1ray istart=',istart
         write(*,*)'magnetic field in initial point'
         write(*,*)'bmod =',bmod,'zu0,ru0,phiu0',zu0,ru0,phiu0
         write(*,*)'in dinit_1ray ksi_nperp', ksi_nperp
ctest for initial conditions
	   xe=x(zu0,ru0,phiu0,1)
	   ye=y(zu0,ru0,phiu0,1)
cendtest
      end if
c--------------------------------------------------------------
      z=zu0
      r=ru0
      phi=phiu0
      cosphi=dcos(phiu0)
      sinphi=dsin(phiu0)

c---------------------------------------------------------------
c     calculations of the parallel (to the magnetic field)
c     refractive index component cnpar
c     N_par=(B_phi*N_phi+N_theta*e_theta*{e_z*B_z+e_r*B_r})/bmod
c     e_theta=(e_z*dpsi/dr-e_r*dpsi/dz)/abs(grad(psi))
c--------------------------------------------------------------
      ppp=(dpdzd*dpdzd+dpdrd*dpdrd)
      gradpsi=dsqrt(dpdzd*dpdzd+dpdrd*dpdrd)
      cnpar1=(cnphi*bphi+cnteta*(bz*dpdrd-br*dpdzd)/gradpsi)/bmod
      cnpar2=cnpar1**2
c Smirnov 961210 beg
c     test of the initial conditions
c      btheta0=(bz*dpdrd-br*dpdzd)/gradpsi
c      write(*,*)'in dinit1_ray btheta0',btheta0,'bphi',bphi,'bmod',bmod
c      write(*,*)'in dinit1_ray cnpat1,cnpar2',cnpar1,cnpar2
c Smirnov 961210 end
c-------------------------------------------------------------------
       write(*,*)'in 1ray cnpar1',cnpar1,'cnteta',cnteta,'cnphi',cnphi
       write(*,*)'i_rho_find_hot_nperp_roots',i_rho_find_hot_nperp_roots
c---------------------------------------------------------------
      if (i_rho_find_hot_nperp_roots.eq.1) then
c-----------------------------------------------------------------
c       finds the small radius rho_ini > rho_min_find_hot_nperp_roots
c       at the vector rho^ where
c       hot plasma dispersdion function D_hot(npar) has three roots.
c       The vector rho^ is starting at the edge point (r_edge,z_edge,phi_edge),
c       and directed to the magnetic axis O(rma,zma,phi_edge)
c-------------------------------------------------------------------
       write(*,*)'dinit.f r,z,phi,cnpar1', r,z,phi,cnpar1
 
        call rho_ini_hot_nperp_roots(r,z,phi,cnpar1)     
c     &  rho_ini,z_ini,r_ini)
        stop 'dinit.f after call rho_ini_hot_nperp_roots'
      endif
c------------------------------------------------------------------    

c-------------------------------------------------------------
c     fit the initial value of rho_ini fort LH or FW cutoff
c--------------------------------------------------------------
      write(*,*)'dinit.f i_rho_cutoff ',i_rho_cutoff
      write(*,*)'cnteta,cnphi',cnteta,cnphi  
      write(*,*)'dinit.f before i_rho_cutoff=1 z,r',z,r
      if (i_rho_cutoff.eq.1) then
         costheta=(r-rma)/dsqrt((z-zma)**2+(r-rma)**2)
         sintheta=(z-zma)/dsqrt((z-zma)**2+(r-rma)**2)
         if(sintheta.gt.0d0) then
            theta=dacos(costheta)
         else
            theta=2*pi-dacos(costheta)
         endif  
         write(*,*)'dinit_1ray before rho_ini i_n_poloidal,n_theta_pol'
     &   ,i_n_poloidal,n_theta_pol
         write(*,*)'cnteta,cnphi',cnteta,cnphi

         call rho_ini_LHFW(theta,phi,cnpar1,
     &   i_n_poloidal,n_theta_pol,cnphi,
     &   rho_ini,z_ini,r_ini,cntheta_ini,cnphi_ini,
     &   i_rho_ini_LHFW_found)

         write(*,*)'dinit.f after rho_ini_LHFW i_rho_ini_LHFW_found=',
     &              i_rho_ini_LHFW_found

         if (i_rho_ini_LHFW_found.eq.1) then
c-----------cutoff point with new z,r coordinates was found
            z=z_ini 
            r=r_ini
            write(*,*)'dinit.f after  rho_ini_LHFW z,r',z,r
         else
c-----------cutoff point was not found
            write(*,*)'cutoff point was not found'
            iraystop=1
            return
         endif
      endif

c      write(*,*)'dinit.f after rho_ini_LHFW,cntheta_ini,cnphi_ini',
c     &cntheta_ini,cnphi_ini
               
      if (i_rho_cutoff.eq.0) then
         cntheta_ini=cnteta
         cnphi_ini=cnphi
      endif

      write(*,*)'dinit.f after rho_ini_LHFW,cntheta_ini,cnphi_ini',
     &cntheta_ini,cnphi_ini

c--------------------------------------------------------------------
c     cninit solves the dispersion relation N=N(n_par)
c     Then subroutine calculates the initial components
c     of the refractive index  cnz,cnr,cm
c---------------------------------------------------------
      write(*,*)'dinit z,r,phi,cnpar1,cntheta_ini,cnphi_ini',
     &z,r,phi,cnpar1,cntheta_ini,cnphi_ini

cBH070123 start
      psi=psif(z,r)
      rho=rhopsi(psi)
      bmod=b(z,r,phi)
!      Rho= 9.535660582E-01  Saveliev starting condition (temporary)
      write(*,*)'dinit r,z,rho,dens,temp ',r,z,rho,densrho(rho,1),
     +          temperho(rho,1)

cSAP090518
      do i=1,nbulk
        write(*,*)'dinit i',i
        dens_i=dense(z,r,phi,i)
        write(*,*)'dinit i,dens_i',i,dens_i
      enddo
      x_e=x(z,r,phi,1)
      y_e=y(z,r,phi,1)
      write(*,*)'x_e,y_e',x_e,y_e

      if(nbulk.ge.2) then
        x_2=x(z,r,phi,2)
        y_2=y(z,r,phi,2)
        write(*,*)'x_2,y_2',x_2,y_2
      endif

      if (nbulk.ge.3)then
        x_3=x(z,r,phi,3) 
        y_3=y(z,r,phi,3)
        write(*,*)'x_3,y_3',x_3,y_3
      endif

      if(nbulk.eq.2) then
         w_lh_d_w=dsqrt(x_2)
         write(*,*)'w_lh_d_w',w_lh_d_w
       endif
       if(nbulk.eq.3) then
         w_lh_d_w=dsqrt(x_2+x_3)
         write(*,*)'w_lh_d_w',w_lh_d_w
       endif
cBH070123 end

      ioxm_n_npar_loc= ioxm_n_npar
      ioxm_n_npar=1
       call nper_npar_ioxm_n_npar(2,z,r,phi,cnpar1,
     & cnper,iraystop) ! ioxm_n_npar will be set in one.i
      write(*,*)'iraystop,ioxm_n_npar,cnper',iraystop,ioxm_n_npar,cnper

      write(*,*)'dinit.f subroutine dinit_1ray ifreq_write',ifreq_write

      if(iraystop.eq.0) then
         wn_perp_ioxm_p(ifreq_write)=cnper
      else
         wn_perp_ioxm_p(ifreq_write)=-1.d0
      endif

       ioxm_n_npar=-1
       call nper_npar_ioxm_n_npar(2,z,r,phi,cnpar1,
     & cnper,iraystop) ! ioxm_n_npar will be set in one.i
      write(*,*)'iraystop,ioxm_n_npar,cnper',iraystop,ioxm_n_npar,cnper

      write(*,*)'dinit.f subroutine dinit_1ray nfreq',nfreq

      write(*,*)'dinit.f subroutine dinit_1ray ifreq_write',ifreq_write

      if(iraystop.eq.0) then
         wn_perp_ioxm_m(ifreq_write)=cnper
      else
         wn_perp_ioxm_m(ifreq_write)=-1.d0
      endif
      ioxm_n_npar= ioxm_n_npar_loc

      write(*,*)'dinit.f i_look_roots',i_look_roots
      if (i_look_roots.eq.1)then   
c-----------------------------------------------------------------
c       plot ReD_hot(nperp) at given npar
c-----------------------------------------------------------------
        name_param(1)='xe'
        name_param(2)='ye'
        name_param(3)='npar'
        param(1)=xe
        param(2)=ye
        param(3)=cnpar1
        write(*,*)'dinit.f before map_d_hot n_param', n_param
        write(*,*)'name_param ',name_param 
        write(*,*)'param ', param
        call map_d_cold(z,r,phi,cnpar1,n_nperp_plot,cnperp_plot_min,
     &  cnperp_plot_max,
     &  name_param,param,n_param)   
        call map_d_hot(z,r,phi,cnpar1,n_nperp_plot,cnperp_plot_min,
     &  cnperp_plot_max,
     &  name_param,param,n_param)
        call map_d_hot(z,r,phi,cnpar1,n_nperp_plot,cnperp_plot_min,
     &  cnperp_plot_max,
     &  name_param,param,n_param)
        write(*,*)'dinit.f after map_d_hot'
        call rho_ini_hot_nperp_roots(r,z,phi,cnpar1)  
c       stop 'after map_d_hot'
c--------------------------------------------------------------
c       calculate all roots N_perpendicular of the hot plasma
c       dispersion function D_hot(N_perp=0) at the interval
c       0 < N_perpendicular < cN_perp_root_max
c----------------------------------------------------------
        call calculate_hot_nperp_roots(z,r,phi,cnpar1,
     &  n_hot_roots,N_perp_root_ar)

c        stop 'dinit.f after calculate_hot_nperp_roots'
        iraystop=1
        return
      endif !   i_look_roots=1    
c-----------------------------------------------------------------
      if ((id.eq.1).or.(id.eq.2)) then
c-------calculate two cold plasma roots for different ioxm
        ioxm_loc=ioxm
        ioxm=+1
        call cninit(z,r,phi,cnpar1,cntheta_ini,cnphi_ini,
     1            cnz,cnr,cm,iraystop)

        if (iraystop.eq.0) then
           cnper=dsqrt(cnz**2+cnr**2+(cm/r)**2-cnpar1**2)   
        else
           cnper =-1.d0
        endif
        wn_perp_ioxm_p(ifreq_write)=cnper

        ioxm=-1
        call cninit(z,r,phi,cnpar1,cntheta_ini,cnphi_ini,
     1            cnz,cnr,cm,iraystop)

        if (iraystop.eq.0) then
           cnper=dsqrt(cnz**2+cnr**2+(cm/r)**2-cnpar1**2)   
        else
           cnper =-1.d0
        endif
        wn_perp_ioxm_m(ifreq_write)=cnper

        ioxm= ioxm_loc
        wye_0(ifreq_write)=y(z,r,phi,1)
        wxe_0(ifreq_write)=x(z,r,phi,1)
cend_test
      endif !id=1,2

      write(*,*)'dinit_1ray before cninit ksi_nperp',ksi_nperp
      call plot_cold_n_perp_omega_npar(z,r,phi,cnpar1,
     &0.1d0,3.d0,1000)

      call cninit(z,r,phi,cnpar1,cntheta_ini,cnphi_ini,
     1            cnz,cnr,cm,iraystop)



      if (iraystop.eq.1) then
         return
      end if

      write(*,*)'dinit after cninit cn**2',cnz**2+cnr**2+(cm/r)**2,
     &'cn=',dsqrt(cnz**2+cnr**2+(cm/r)**2)

      

      cn2p=cnz**2+cnr**2+cnphi**2
      cnphi=cm/r
      cn2=cnz**2+cnr**2+cnphi**2
      cnper=dsqrt(cn2-cnpar1**2)
      write(*,*)'in dinit cn2',cn2,'cnper,cnpar1',cnper,cnpar1

      bmod=b(z,r,phi)
      gam=gamma1(z,r,phi,cnz,cnr,cm)
      write(*,*)'1ray before d=hamilt1'
      write(*,*)'z,r,phi before d=hamilt1 z,r,phi',z,r,phi
      write(*,*)'z,r,phi before d=hamilt1 cnz,cnr,cm',cnz,cnr,cm

      dh=hamilt1(z,r,phi,cnz,cnr,cm)

      write(*,*)'dinit_1ray after d=hamilt1 dh=',dh
c The check of the Hamiltonian value for the found initial conditions.
      epshamin=1.d-6
      epshamin=1.d-2
      epshamin=1.d-1
     
c-------------------------------------------------------------------

      write(*,*)'before outinit'
      call outinit(u)
      write(*,*)'in dinit_1ray after call outini nrayelt= ',nrayelt
      irefl=0
      u(1)=z
      u(2)=r
      u(3)=phi
      u(4)=cnz
      u(5)=cnr
      u(6)=cm
      write(*,*)'dinit_1ray before prep3d powini',powini
      write(*,*)'dinit_1ray before prep3d u',u

      write(*,*)'dinit_1ray before prep3d rma,zma',rma,zma
      call prep3d(0.d0,u,deru,iraystop) ! YuP: 0.0->0.d0

    
      write(*,*)'initial data'
      write(*,*)'z=',z,'r=',r,'phi=',phi
      !write(*,*)'cnz=',cnz,'cnr=',cnr,'cm=',cm
      write(*,*)'rho=',rho

c-----safety factor calculations
      psi_initial=psi_rho(rho)
      q_initial=qsafety_psi(psi_initial)
      write(*,*)'psi_initial,q_initial',psi_initial,q_initial
      write(*,*)'b(z,r,phi)',b(z,r,phi)
cyup      b_av=b_average(psi_initial)
cyup      write(*,*)'b_average(psi_initial)',b_av
c     stop 'dinit_1ray'
      write(*,*)'end dinit_1ray'
      return
      end




c        **********************zeffcalc************************
c        *                        -                           *
c        * this subroutine calculates zeff(ndens)             *
c        * using the ions densities 			      *
c        ******************************************************
c
c------------------------------------------------------------------
c        output parameters: array zeff(ndens)
c------------------------------------------------------------------

      subroutine zeffcalc
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'six.i'
c---------------------------------------------------------
c     calculation of the table for the radial profile zeff1(ndens)
c---------------------------------------------------------
      if (nbulk.eq.1) then
         write(*,*)'Warning dinit: nbulk=1, it will be created zeff=1'
      endif

      do j=1,ndens
         if (nbulk.eq.1) then
            zeff1(j)=1.d0
            goto 10
         endif
c---------------------------------------------------------
c        electron density : dens_el
         dens_el=dens1(j,1)
         zeff1(j)=0.d0
	 do i=2,nbulk
           dens_ion=dens1(j,i)
           if (dens_el.ne.0.d0) then
              denside=dens_ion/dens_el
           else
              if (dens_ion.eq.0.d0) then
                 denside=1.d0
              else
                 write(*,*)'in zeffcalc j',j
	         write(*,*)'densel=0 but dens_ion(i).ne.0'
	         write(*,*)'i number of plasma component=',i
	     write(*,*)'change the parameters of the density profiles'
	         stop
	      endif
	   endif

	   zeff1(j)=zeff1(j)+charge(i)*charge(i)*denside
         enddo !i=2,nbulk
 10      continue
c---------------------------------------------------------
      enddo !j=1,ndens
c---------------------------------------------------------
c     end of calculation of the table for the radial profile
      return
      end
c        **********************zeffcal1************************
c        *                        -                           *
c        * this subroutine calculates zeff(ndens)             *
c        * using the ions densities. 			      *
c        *It is the old version FOR the charge neutrality control
c        ******************************************************
c
c------------------------------------------------------------------
c        output parameters: array zeff(ndens)
c------------------------------------------------------------------

      subroutine zeffcal1
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'six.i'
c---------------------------------------------------------
c     calculation of the table for the radial profile zeff1(ndens)
c---------------------------------------------------------
      h=1.d0/(ndens-1.d0)
      if (nbulk.eq.1) then
        write(*,*)'zeffcal1: WARNING nbulk=1, it will be set zeff=1'
      endif

      do j=1,ndens
         rho=h*(j-1)
         if (nbulk.eq.1) then
            zeff1(j)=1.d0
            goto 10
         endif
c---------------------------------------------------------
c        electron density : dens_el
         if (idens.eq.0) then
c           analytical representation of the electron density profile
            dens_el=(dense0(1)-denseb(1))*
     1              (1-rho**rn1de(1))**rn2de(1)+denseb(1)
         else
c           table representation of the electron density profile
            dens_el=dens1(j,1)
         endif
c---------------------------------------------------------
         zeff1(j)=0.d0
         qusneutr=0.d0 ! for the control of the charge-neutrality
         do i=2,nbulk
            if (idens.eq.0) then
c              analytical representation of the ion density profiles
               dens_ion=(dense0(i)-denseb(i))*
     1                  (1-rho**rn1de(i))**rn2de(i)+denseb(i)
            else
c              table representation of the ion density profiles
               dens_ion=dens1(j,i)
            endif

            if (dens_el.ne.0.d0) then
               denside=dens_ion/dens_el
            else
               if (dens_ion.eq.0.d0) then
                 denside=1.d0
               else
               write(*,*)'in zeffcalc rho=',rho
               write(*,*)'densel=0 but dens_ion(i).ne.0'
               write(*,*)'i number of plasma component=',i
               write(*,*)'change the parameters of the density profiles'
               stop
               endif
            endif
c---------------------------------------------------------
c   control of the charge-neutrality condition in the point rho(i) 
c---------------------------------------------------------
            qusneutr=qusneutr+charge(i)*denside
            zeff1(j)=zeff1(j)+charge(i)*charge(i)*denside
         enddo !i=2,nbulk
c---------------------------------------------------------
         if (qusneutr.ne.1.d0) then
            write(*,*)'in zeffcalc the bad quasi-neutrality'
            write(*,*)'in rho(j)=',rho,'j=',j
            write(*,*)'sum({i=2,nbulk}(charge(i)*dens(i(i)/dens_e)',
     1	    qusneutr
            write(*,*)'change the charge(j) or the densities profiles'
    	    stop
         endif
 10      continue
      enddo !j=1,ndens
c---------------------------------------------------------
c     end of calculation of the table for the radial profile
      return
      end


c     **********************denscalc************************
c     *                        -                           *
c     * this subroutine calculates ions densities          *
c     * dens1(j,nbulk) and dens1(j,nbulk-1)                *
c     * (here j=1,ndens),				   *
c     * using the zeff and dens1(ndens,i) i=2,nbulk-2      *
c     ******************************************************
c
c------------------------------------------------------------------
c     output parameters: array dens1(ndens,nbulk)
c------------------------------------------------------------------

      subroutine denscalc
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'six.i'
c---------------------------------------------------------
c     calculation of the table for the radial profile
c     dens1(ndens,nbulk), dens1(ndens,nbulk-1)
c---------------------------------------------------------
cc      write(*,*)'in denscalc ndens,nbulk,idens,izeff'
cc      write(*,*)ndens,nbulk,idens,izeff
      if (nbulk.le.2) then
         write(*,*)'in denscalc nbulk.le.2'
         write(*,*)'zeff must be equal charge(2)'
         write(*,*)'and dense(2)=dense(1)/charge(2)'
         write(*,*)'use the option izeff=0 and these parameres'
	 stop
      else
c        nbulk.ge.3
         if( charge(nbulk).eq.charge(nbulk-1)) then
	    write(*,*)'Warning in denscalc: nbulk(.ge.3)=',nbulk
	    write(*,*)'in denscalc: charge(nbulk)=charge(nbulk-1)'
	    write(*,*)'it is impossible to find the ions denscalc'
	    write(*,*)'change charge(nulk) or charge(nbulk-1)'
	    write(*,*)'it should be charge(nulk)>charge(nbulk-1)'
	    write(*,*)'or use the option izeff=0'
	    stop
	 endif
      endif

      h=1.d0/(ndens-1.d0)
      do j=1,ndens
         rho=h*(j-1)
cc         write(*,*)'in denscalc j,rho',j,rho
c---------------------------------------------------------
c        electron density : dens_el
	 if (idens.eq.0) then
c           analytical representation of the electron density profile
	    dens_el=(dense0(1)-denseb(1))*
     1              (1-rho**rn1de(1))**rn2de(1)+denseb(1)
	 else
c           table representation of the electron density profile
	    dens_el=dens1(j,1)
	 endif
c---------------------------------------------------------
         sum1=0.d0
         sum2=0.d0
cc	 write(*,*)'in denscalc j=',j,'rho',rho,'dens_el',dens_el
	 if(nbulk.ge.4) then
c           sum1=sum{i=2,nbulk-2}(density(i)/electron_density*
c                charge(i)*(charge(nbulk)-charge(i))
c           sum2=sum{i=2,nbulk-2}(density(i)/electron_density*
c                charge(i)*(charge(nbulk-1)-charge(i))
            do i=2,nbulk-2
	       if ((charge(nbulk)-charge(i)).lt.0) then
	         write(*,*)'in denscalc charge(nbulk).lt.charge(i)'
		 write(*,*)'change array charge in genray.in'
		 write(*,*)'in array charge(i+1) must be >charge(i)'
		 write(*,*)'change the array charge(i) in genray.in'
		 stop
	       endif
	       if ((charge(nbulk-1)-charge(i)).lt.0) then
	         write(*,*)'in denscalc charge(nbulk-1).lt.charge(i)'
		 write(*,*)'change array charge in genray.in'
		 write(*,*)'in array charge(i+1) must be >charge(i)'
		 write(*,*)'change the array charge(i) in genray.in'
		 stop
	       endif
	       if (idens.eq.0) then
c                 analytical representation of the ion density profiles
	          dens_ion=(dense0(i)-denseb(i))*
     1                     (1-rho**rn1de(i))**rn2de(i)+denseb(i)
	       else
c                 table representation of the ion density profiles
	          dens_ion=dens1(j,i)
	       endif

               if (dens_el.ne.0.d0) then
	          denside=dens_ion/dens_el
	       else
	          if(dens_ion.eq.0.d0) then
	            denside=1.d0
		  else
		    write(*,*)'in denscalc dens_el=0,
     1		    but dens_ion.ne.0, j=',j,'i=',i
                    write(*,*)'change density profiles parameters'
		    stop
	          endif
	       endif
	       sum1=sum1+denside*charge(i)*(charge(nbulk)-charge(i))
	       sum2=sum2+denside*charge(i)*(charge(nbulk-1)-charge(i))
	    enddo !nbulk
	 endif ! nbulk.ge.4
	 p3=charge(nbulk)-charge(nbulk-1)
	 if (p3.lt.0.d0) then
	    write(*,*)'in denscalc charge(nbulk).lt.charge(nbulk-1)'
	    write(*,*)'change array charge in genray.in'
	    write(*,*)'in array charge(nbulk) must be >charge(nbulk-1)'
	    write(*,*)'change the array charge(i) in genray.in'
	    stop
	 endif
	 p1=charge(nbulk)-zeff1(j)-sum1
	 if (p1.lt.0.d0) then
	    write(*,*)'in denscalc charge(nbulk)-zeff1(j)-sum1.lt.0'
	    write(*,*)'j=',j,'it is impossible to calculate
     1	    the density(j,nbulk-1) using zeff profile'
	    write(*,*)'change the array charge(i) in genray.in'
	    stop
	 endif
	 p2=zeff1(j)-charge(nbulk-1)+sum2
	 if (p2.lt.0.d0) then
	    write(*,*)'in denscalc charge(nbulk)-zeff1(j)-sum1.lt.0'
	    write(*,*)'j=',j,'it is impossible to calculate
     1	    the density(j,nbulk) using zeff profile'
	    write(*,*)'change the array charge(i) in genray.in'
	    stop
	 endif
	 dens1(j,nbulk-1)=p1/(charge(nbulk-1)*p3)*dens1(j,1)
	 dens1(j,nbulk)=p2/(charge(nbulk)*p3)*dens1(j,1)
cc         write(*,*)'in denscalc dens1(j,nbulk-1),dens1(j,nbulk)'
cc         write(*,*)dens1(j,nbulk-1),dens1(j,nbulk)
      enddo ! j
 10   continue
      return
      end
      
      
      
c     **********************denscalp************************
c     * this subroutine calculates electron and		   *
c     * ions densities profiles: dense(j,1),               *
c     * dens1(j,nbulk) and dens1(j,nbulk-1)                *
c     * (here j=1,ndens),				   *
c     * using the zeff and tempe1(ndens,i) i=2,nbulk       *
c     * and eqdsk pres (pressure)
c     ******************************************************
c
c------------------------------------------------------------------
c     output parameters: array dens1(ndens,nbulk)
c------------------------------------------------------------------

      subroutine denscalp
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'six.i'
c---------------------------------------------------------
c     calculation of the table for the radial profile
c     dens1(ndens,nbulk), dens1(ndens,nbulk-1)
c---------------------------------------------------------
cc      write(*,*)'in denscalp ndens,nbulk,idens,izeff'
cc      write(*,*)ndens,nbulk,idens,izeff

      if (nbulk.eq.2) then
         write(*,*)'in denscalp nbulk.eq.2'
         write(*,*)'use another option izeff or another parameters'
         stop
      else
        if (nbulk.ge.3) then
c        nbulk.ge.3
         if( charge(nbulk).eq.charge(nbulk-1)) then
	    write(*,*)'Warning in denscalc: nbulk(.ge.3)=',nbulk
	    write(*,*)'in denscalp: charge(nbulk)=charge(nbulk-1)'
	    write(*,*)'it is impossible to find the ions denscalp'
	    write(*,*)'change charge(nulk) or charge(nbulk-1)'
	    write(*,*)'it should be charge(nulk)>charge(nbulk-1)'
	    write(*,*)'or use another option izeff'
	    stop
         endif
        endif
      endif

      h=1.d0/(ndens-1.d0)
      do j=1,ndens
         rho=h*(j-1)
         psi=psi_rho(rho)
         pressure=prespsi(psi)/1.6d3
         if(nbulk.eq.1) then
c----------electron density n=p/2T)
           dens1(j,1)=pressure/(2.d0*temp1(j,1))
           goto 10
         endif
c---------------------------------------------------------
         sum1=0.d0
         sum2=0.d0
         sum4=0.d0
cc	 write(*,*)'in denscalc j=',j,'rho',rho,'dens_el',dens_el
	 if(nbulk.ge.4) then
c           sum1=sum{i=2,nbulk-2}(density(i)*charge(i)*
c                                 (charge(nbulk)-charge(i))
c           sum2=sum{i=2,nbulk-2}(density(i)*charge(i)*
c                                 (charge(nbulk-1)-charge(i))
c           sum4=sum{i=2,nbulk-2}(density(i)*temp(i))
            do i=2,nbulk-2
	       if ((charge(nbulk)-charge(i)).lt.0) then
	         write(*,*)'in denscalp charge(nbulk).lt.charge(i)'
		 write(*,*)'change array charge in genray.in'
		 write(*,*)'in array charge(i+1) must be >charge(i)'
		 write(*,*)'change the array charge(i) in genray.in'
		 stop
	       endif
	       if ((charge(nbulk-1)-charge(i)).lt.0) then
	         write(*,*)'in denscalp charge(nbulk-1).lt.charge(i)'
		 write(*,*)'change array charge in genray.in'
		 write(*,*)'in array charge(i+1) must be >charge(i)'
		 write(*,*)'change the array charge(i) in genray.in'
		 stop
	       endif
	       if (idens.eq.0) then
c                 analytical representation of the ion density profiles
	          dens_ion=(dense0(i)-denseb(i))*
     1                     (1-rho**rn1de(i))**rn2de(i)+denseb(i)
	       else
c                 table representation of the ion density profiles
	          dens_ion=dens1(j,i)
	       endif

	       sum1=sum1+dens_ion*charge(i)*(charge(nbulk)-charge(i))
	       sum2=sum2+dens_ion*charge(i)*(charge(nbulk-1)-charge(i))
	       sum4=sum4+dens_ion*temp1(j,i)
	    enddo !nbulk
	 endif ! nbulk.ge.4
	 p3=charge(nbulk)-charge(nbulk-1)
	 if (p3.lt.0.d0) then
	    write(*,*)'in denscalp charge(nbulk).lt.charge(nbulk-1)'
	    write(*,*)'change array charge in genray.in'
	    write(*,*)'in array charge(nbulk) must be >charge(nbulk-1)'
	    write(*,*)'change the array charge(i) in genray.in'
	    stop
	 endif
c----------------------------------------------------------------
c the formula for dense(nbulk-1) and dense(nbulk)
c dense(nbulk-1)=((charge(nbulk)-zeff)*dense(1)-sum1))/
c                (charge(nbulk-1)*(charge(nbulk)-charge(nbulk-1)))
c dense(nbulk)  =((zeff-charge(nbulk-1))*dense(1)+sum2))/
c                (charge(nbulk)*(charge(nbulk)-charge(nbulk-1)))
c----------------------------------------------------------------
c calculation of the electron density dens1(j,1)
c from the the plasma pressure:
c pressure=dens1(j,1)*temp1(j,1)+sum(i=2,nbulk1)(dens1(j,i)*temp1(j,i))
c    +dens1(j,nbulk-1)*temp1(j,nbulk-1)+dens1(j,nbulk)*temp1(j,nbulk)
         p4=temp1(j,1)+
     1   temp1(j,nbulk-1)*(charge(nbulk)-zeff1(j))/(charge(nbulk-1)*p3)
     1   +temp1(j,nbulk)*(zeff1(j)-charge(nbulk-1))/(charge(nbulk)*p3)
         p5=pressure-sum4+
     1   temp1(j,nbulk-1)*sum1/(charge(nbulk-1)*p3)-
     1   temp1(j,nbulk)*sum2/(charge(nbulk)*p3)
         if(p4.eq.0.d0) then
	   write(*,*)'in denscalp:izeff=4 p4=0 bad conditions j=',j
	   stop
	 else
	   dens1(j,1)=p5/p4
	   if(dens1(j,1).lt.0.d0) then
	     write(*,*)'in denscalp:izeff=4, des1(j,1)<0,j=',j
	     write(*,*)'bad conditions'
	     stop
	   endif
	 endif
         p1=dens1(j,1)*(charge(nbulk)-zeff1(j))-sum1
	 if (p1.lt.0.d0) then
	    write(*,*)'in denscalp j=',j
	    write(*,*)'dens1(j,1)*(charge(nbulk)-zeff1(j))-sum1.lt.0'
	    write(*,*)'j=',j,'it is impossible to calculate
     1	    the density(j,nbulk-1) using zeff profile'
	    write(*,*)'change the array charge(i) in genray.in'
	    stop
	 endif
	 p2=dens1(j,1)*(zeff1(j)-charge(nbulk-1))+sum2
	 if (p2.lt.0.d0) then
	    write(*,*)'in denscalp j=',j
	    write(*,*)'dens1(j,1)*(charge(nbulk)-zeff1(j))-sum1.lt.0'
	    write(*,*)'j=',j,'it is impossible to calculate
     1	    the density(j,nbulk) using zeff profile'
	    write(*,*)'change the array charge(i) in genray.in'
	    stop
	 endif
	 dens1(j,nbulk-1)=p1/(charge(nbulk-1)*p3)
	 dens1(j,nbulk)=p2/(charge(nbulk)*p3)
         write(*,*)'in denscalp: j=',j
	 write(*,*)'dens1(j,1),dens1(j,nbulk-1),dens1(j,nbulk)'
	 write(*,*)dens1(j,1),dens1(j,nbulk-1),dens1(j,nbulk)

 10      continue

      enddo ! j
      return
      end




c--------------mapxyb---------------------------------
      subroutine mapxyb(nrhomap,nthetmap)
   
c     (1) It creates the data for plots
c     X_e, Y_e B_tot, B_ptor, B_pol (rho,theta)
c     the results are in the output files
c     xybrhoth.dat: rho(i),theta(j),xe,ye,(xe+ye*ye),bmod,bphi,
c     *              dsqrt(bz**2+br**2)
c     Here (0<rho<1,0<theta<pi) are the points of mesh

c      implicit none
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'

c-----input
      integer nrhomap,nthetmap  ! The number of mesh points (rho,theta) 
     
      double precision psi_rho,b,x,y
      external psi_rho  ! calculates psi(rho)
      external zr_psith ! calculates (z,r) at given (psi,theta_poloidal)
      external b,x,y
      
c-----local
      double precision hthetmap,hrhomap,
     * thetmap,rhomap,psix,
     * rmap,zmap,phi,xe,ye
      integer i,j

      write(*,*)'dinit begin mapxyb nrhomap,nthetmap',nrhomap,nthetmap

c-----steps of  mesh      
      hrhomap=(1.d0-1.d-4)/(nrhomap-1)
      write(*,*)'mapxyb hrhomap',hrhomap
        
      pi=4.d0*datan(1.d0) 
      hthetmap=pi/(nthetmap-1)
      write(*,*)'mapxyb thetmapn',hthetmap

      open(1,file='xybrhoth.dat')
     
1     format (8(' ',1pe11.4))      
      write(1,*)' rho theta xe ye uh bmod bphi bpol '
      phi=0.d0
      do i=1,nrhomap
        rhomap=hrhomap*(i-1)
        write(*,*)'i,rhomap',i,rhomap
        psix=psi_rho(rhomap)

        do j=1,nthetmap
           thetmap=hthetmap*(j-1)
           write(*,*)'j,thetmap',j,thetmap
           call zr_psith(psix,thetmap,zmap,rmap)
           bmod=b(zmap,rmap,phi)
           write(*,*)'bmod',bmod
           xe=x(zmap,rmap,phi,1)
           ye=y(zmap,rmap,phi,1)
           write(1,1)rhomap,thetmap,xe,ye,(xe+ye*ye),bmod,bphi,
     *              dsqrt(bz**2+br**2)
        enddo
      enddo
                    
      close(1)

      return
      end





