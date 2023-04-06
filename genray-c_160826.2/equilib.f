c        *********************EQUILIB************************
c        *  EQUILIB is the code for reading data for       *
c        *  toroidally symmetric plasma equilibrium given   *
c        *  the EQDSK (ref. Lang Lao) format.              *
c        *  It transforms these data for input data for    *
c        *  GENRAY code.				   *
c        *  RAYINP normalizes the eqdsk parameters by	   *
c        *  r0x(in m) and b0(in Tl) and			   *
c        *  creates the arrays for spline codes which      *
c        *  calculates: the (psi) and (feqd) magnetic field   *
c        *  functions; zlimit upper(r) and zlimit under(r) *
c        ***************************************************
c
c-------------------------------------------------------------------!
c        input data are read from:				    !
c        equilib.dat=       eqdsk  file    			    !
c        genray.in=        a namelist file with the normalization  !
c                           parameters r0x and b0                   !
c-------------------------------------------------------------------!
c======> Yu.P. 2011 Added: model_b=1 option.
c        Magnetic field is set as  Bz=uniform, Br=0, Bphi~1/r
c        Skip reading eqdsk in this case.
c======> Yu.P. 2011 Added: model_b=2 option.
c        Magnetic field is calculated from a set of magnetic coils
c        (current loops). 
c        Toroidal Bphi~1/r can be added, similar to model_b=1.
c        Skip reading eqdsk in this case.
c======> Yu.P. 2012   Added: model_b=4 option.
c        FRC-type magnetic field. Components are set as
c         b_x= 0.d0
c         b_y= 0.d0
c         b_z= bz0*tanh( akappa*[2*(r/rs_frc)**2  - 1.] )
c        where
c         rs_frc  - separatrix radius (specify in genray.in) [m]
c         bz0     - magnetic field at r=+INF, approximately at wall [T]
c         akappa  - such "Kappa" that  beta_frc*Kappa-tanh(Kappa)=0
c         beta_frc= 1.d0 - 0.5*(rs_frc/wall_rmax)**2
c        The value of akappa is found in subroutine input :
c        akappa= rtbis(func_kappa, 1.d-5, 10.d0, 1.d-14)
c        Last three arguments: initial guesses that bracket the root
c        and the accuracy for searching the root. Adjust if needed.
c        Also, the Toroidal Bphi~1/r can be added, similar to model_b=1.
c        Such FRC-magnetic field is associated with the density profile
c        (model_rho_dens=4)  (FRC-like plasma)
c        !-1-> Rigid Rotor profile:
c         dens_rr= dens0rr*( sech(akappa*rho) )**2
c        where rho=abs(2*(r/rs_frc)**2 -1.) 
c        ! With such definition, rho=0 at r=rs_frs/sqrt(2), 
c        ! and rho=1 at r=0 or r=rs_frc .
c        ! Note: At r=0 and r=rs_frc, dens_rr= dens0rr*[sech(akappa)]^2
c        !-2-> Uniform background density profile is also added:
c         dens_ub= dens0ub 
c        if (r .gt. rs_frc) then 
c           ! Linear drop in region  rs_frc < r < wall_rmax
c           dens_ub= dens0ub*(1.d0 -(r-rs_frc)/(wall_rmax-rs_frc))
c        endif
c        dense_xyz= dens_rr + dens_ub !electrons only, for now
c-------------------------------------------------------------------!
c********************************************************************
c  This program uses following external files:
c  input , dinitr , limitr 
c*********************************************************************
      subroutine equilib
      implicit none 
        call input
        call dinitr      
      return
      end




c        ********************INPUT**************************
c        *   subroutine INPUT reads the eqdsk data , 	   *
c        *   normalizes the eqdsk parameters by 	   *
c        *   r0x(in m) and b0(in Tl)                       *
c        ***************************************************
c
c-------------------------------------------------------------------!
c        input data are read from:				    !
c        eqdsk_cq           file                 		    !
c        genray.in         it is a  file to read  the normalization!
c                           parameters r0x and b0                   !
c        output data  for subroutine limitr() are writen in:	    !
c        common 'three' and 'fourb'                                 !
c-------------------------------------------------------------------!
c======> Yu.P. 2011 Added: model_b=1 option.
c        Magnetic field is set as  Bz=uniform, Br=0, Bphi~1/r
c        Skip reading eqdsk in this case.
c======> Yu.P. 2011 Added: model_b=2 option.
c        Magnetic field is calculated from a set of magnetic coils
c        (current loops). 
c        Toroidal Bphi~1/r can be added, similar to model_b=1.
c        Skip reading eqdsk in this case.
c======> Yu.P. 2016 Added: model_b=3 option. See GENRAY-c_help.txt.
c======> Yu.P. 2012   Added: model_b=4 option.
c        FRC-type magnetic field. Components are set as
c         b_x= 0.d0
c         b_y= 0.d0
c         b_z= bz0*tanh( akappa*[2*(r/rs_frc)**2  - 1.] )
c        where
c         rs_frc  - separatrix radius (specify in genray.in) [m]
c         bz0     - magnetic field at r=+INF, approximately at wall [T]
c         akappa  - such "Kappa" that  beta_frc*Kappa-tanh(Kappa)=0
c         beta_frc= 1.d0 - 0.5*(rs_frc/wall_rmax)**2
c        The value of akappa is found in subroutine input :
c        akappa= rtbis(func_kappa, 1.d-5, 10.d0, 1.d-14)
c        Last three arguments: initial guesses that bracket the root
c        and the accuracy for searching the root. Adjust if needed.
c        Also, the Toroidal Bphi~1/r can be added, similar to model_b=1.
c        Such FRC-magnetic field is associated with the density profile
c        (model_rho_dens=4)  (FRC-like plasma)
c        !-1-> Rigid Rotor profile:
c         dens_rr= dens0rr*( sech(akappa*rho) )**2
c        where rho=abs(2*(r/rs_frc)**2 -1.) 
c        ! With such definition, rho=0 at r=rs_frs/sqrt(2), 
c        ! and rho=1 at r=0 or r=rs_frc .
c        ! Note: At r=0 and r=rs_frc, dens_rr= dens0rr*[sech(akappa)]^2
c        !-2-> Uniform background density profile is also added:
c         dens_ub= dens0ub 
c        if (r .gt. rs_frc) then 
c           ! Linear drop in region  rs_frc < r < wall_rmax
c           dens_ub= dens0ub*(1.d0 -(r-rs_frc)/(wall_rmax-rs_frc))
c        endif
c        dense_xyz= dens_rr + dens_ub !electrons only, for now
c-------------------------------------------------------------------!
      subroutine input

c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none 

      include 'param.i'
      include 'three.i'
      include 'fourb.i'
      include 'one.i'
      include 'onetwo.i'

* nreqd,nzeqd      are the numbers of points in the horisontal (r)
*                  and vertical (z) directions.
* rdimeqd,zdimeqd  are the horisontal and vertical full-widths
*                  of the  rectangle [meters]
* reqd             is the nominal major radius of the torus.
* redeqd           is the major radius of the inner edge
*                  of rectangular grid.
* zmideqd          is the vertical shift of the rectangular box
*                  up-down simmetry plane.
* rma,zma          are the major and vertical height of magnetic axis.
* psimag,psilim    are the poloidal flux function values at the
*                  magnetic axis and the last closed flux surface
*                   (touching the limiter or the separatrix).
* beqd             is the toroidal magnetic field at reqd.
* toteqd           is the toroidal current.
* psimx1,psixm2,xax1,xax2,zax1,zax2,psisep,xsep,ysep - OBSOLETE.
*  ATTENTION:
* psimx1=psimag
* xax1=rma
* zax1=zma
* psisep=psilim
c       feqd(nreqda),pres(nreqda)
* feqd = r * B_phi  at nreqd equispaced points in psi
*                   from psimag to psilim.
* p                 are the pressure values at the same points.
c       ffpeqd(nreqda),ppeqd(nreqda)
* ffpeqd           (=  feqd_prime) at the same points.
* ppeqd            (=  p_prime) at the same points.
c       peqd(nreqda,nzeqda)
* peqd(nreqd,nzeqd)  are the psi values on the nreqd * nzeqd
*                     equispaced grid.
* qpsi(nreqda)        are q characteristic.
* nnlim,nnves        are the numbers of point at limiters and
*                    vacuum vessel wall.
c       rlimit(nnlim),zlimit(nnlim)
* rlimit,zlimit      is the r,z location of the limiters wall.
c       rves(nnves),zves(nnves)
* rves,zves          is the r,z location of the limiters
*                    vacuum vessel wall.


      character*12 namfil
      
      include 'name_genr.i'
      include 'name_tokamak.i'
      integer kode,i,j
      real*8  dpsiar,dflux,pr,pb,dpsi,dstep, x,y,r,z,dra,
     +        b_x,b_y,b_z,sumpsi, rtbis, func_kappa, bxyz
      real*8 R_b0, PSI, PSI_b0, dreq
      integer ir_b0

      external length_char, func_kappa
      integer length_char

c     Data from /genr/ and /tokamak/ namelists were read
c     in genray.f using read_all_namelists
      goto 10
c-----------------------------------------------------------
c     input /genr/ and /tokamak/ namelists
c-----------------------------------------------------------
      open(1,file='genray.in',status='old',iostat=kode)
      if (kode.ne.0) then
         open(1,file='genray.dat',status='old',iostat=kode)
         if (kode.ne.0) then
            write(*,*)'dinit:Neither genray.in or genray.dat r present'
            stop
         endif
      endif

      write(*,*)'equilib.f before read genr'
      rewind(unit=1)
      read(1,genr,iostat=kode)
      write(*,*)'equilib.f after read genr kode=',kode
      write(*,genr)
      call check_read(kode,'genr')
 
      rewind(unit=1)
      read(1,tokamak,iostat=kode)
      write(*,*)'equilib.f tokamak kode=',kode
      write(*,tokamak)
      call check_read(kode,'tokamak')

c-----check the input data in namelist /tokamak/
      if ((indexrho.lt.1).or.(indexrho.gt.6)) then
         write(*,*)'in equilib.f in reading namelist /tokamak/'
         write(*,*)'It should be 0<indexrho<7, but indexrho =',indexrho
         write(*,*)'Change indexrho in genray.in file'
         stop
      endif
      
      if(model_b.eq.0)then
       if((indexrho.le.3).or.(indexrho.eq.6)) then
       write(*,*) 'For model_b=0, only indexrho=4 or 5 are recommended.'
       write(*,*) 'Values of rho>1 are only defined for indexrho=4,5.'
       endif
      else
       write(*,*) 'indexrho set to 4 (the only option for model_b.ne.0)'
       indexrho=4 ! rho=sqrt((psi-psimag)/psilim-psimag))
       ! Technically, the value of indexrho 
       ! is not used in this case. 
       ! But it is set to 4 anyway, for a print-out or storage.
      endif

      if ((ipsi.lt.0).or.(ipsi.gt.1)) then
         write(*,*)'in equilib.f in reading namelist /tokamak/'
         write(*,*)'It should be -1<ipsi<2, but ipsi =',ipsi
         write(*,*)'Change ipsi in genray.in file'
         stop
      endif

      if ((ionetwo.lt.0).or.(ionetwo.gt.1)) then
         write(*,*)'in equilib.f in reading namelist /tokamak/'
         write(*,*)'ionetwo can be 0 or 1, but ionetwo=',ionetwo
         write(*,*)'Resetting ionetwo to 1'
         ionetwo=1
      endif

      if ((ieffic.lt.1).or.(ieffic.gt.6)) then
         write(*,*)'in equilib.f in reading namelist /tokamak/'
         write(*,*)'It should be 0<ieffic<7, but ieffic =',ieffic
         write(*,*)'Change ieffic in genray.in file'
         stop
      endif

      if ((psifactr.le.0).or.(psifactr.gt.1)) then
         write(*,*)'in equilib.f in reading namelitsst /tokamak/'
         write(*,*)'It should be 0<psifactr=<1, but psifactr=',psifactr
         write(*,*)'Change psifactr in genray.in file'
         stop
      endif

      if (length_char(eqdskin).gt.512) then
         write(*,1001)
 1001    format('STOP: eqdskin spec too long (max:512)')
         STOP
      endif
    
      if (length_char(eqdsktype).gt.16) then
         write(*,*) 'STOP: eqdsktype spec too long (max:16)'
         STOP
      endif

      if (NR.gt.NRA) then
        write(*,*)'NR > NRA'
        write(*,*)'Should be NR.le.NRA'
        write(*,*)'Change NR in genray.dat or NRA in param.i' 
        STOP 
      endif

      close(1)

c-----end check of namelist /tokamak/

      write(*,*)'in equilib r0x,b0',r0x,b0
      write(*,*)'ipsi',ipsi

      write(*,*) ' End of input'

 10   continue
 
      if(model_b.ne.0) goto 12 ! skip reading eqdskin
      
c-----------------------------------------------------------
c     Read the EQDSK file
c-----------------------------------------------------------
      write(*,*)'eqdskin= ',eqdskin

      open(30,file=eqdskin)

cBH020822  Adding nveqd (.le. nreqd) for different number of
cBH020822  flux surfaces on which p,feqd,p',feqd',q are tabulated.
cBH020822  The standard EQDSK does not incorporate this feature,
cBH020822  although is is available in cql3d.
cBH020822  Bonoli uses it for ACCOME eqdsk output.

c      read(30,2) nreqd,nzeqd
      nveqd=0
      write(*,*)'equilib.f before read (30,2)'
      read(30,2) nreqd,nzeqd,nveqd
      write(*,*)'nreqd,nzeqd,nveqd',nreqd,nzeqd,nveqd
2     format(52x,3i4)

      nxeqd= 2*nreqd  ! YuP: x-grid size in cartesian coords
      nyeqd= 2*nreqd  ! YuP: y-grid size in cartesian coords

      if (nveqd.gt.nreqd) stop 'nveqd.gt.nreqd NOT ENABLED'
      if (nveqd.eq.0) nveqd=nreqd

      call check_param(1)
c      if ((nreqd.gt.nreqda).or.(nzeqd.gt.nzeqda)) then
c        write(*,1000) nreqd,nreqda,nzeqd,nzeqda
c 1000   format('in equilib.dat in input',/,
c     .  'the dimensions of eqdsk (in eqilib.dat) nreqd or nzeqd',/,
c     .  'are bigger than the parameters nreqda or nzeqda in param.i'
c     .  ,/'nreqd=',I5,'nreqda=',I5
c     .  ,/'nzeqd=',I5,'nzeqda=',I5
c     .  ,/'Please change nreqda or nzeqda in param.i')
c        stop
c      endif

      read(30,3) rdimeqd,zdimeqd,reqd,redeqd,zmideqd
      write(*,*)'rdimeqd,zdimeqd,reqd,redeqd,zmideqd'
      write(*,*) rdimeqd,zdimeqd,reqd,redeqd,zmideqd
3     format(5e16.9)

      xeqmax= redeqd+rdimeqd ! YuP max for x-grid of cartesian grid
      xeqmin=-xeqmax
      yeqmax= redeqd+rdimeqd ! YuP max for y-grid of cartesian grid
      yeqmin=-yeqmax
      zeqmax= zmideqd + 0.5d0*zdimeqd
      zeqmin= zmideqd - 0.5d0*zdimeqd
      write(*,*)'xeqmin,xeqmax=',xeqmin,xeqmax
      write(*,*)'yeqmin,yeqmax=',yeqmin,yeqmax
      write(*,*)'zeqmin,zeqmax=',zeqmin,zeqmax
      
      peqd=0.d0 ! initialize, just in case
      beq=0.d0  ! initialize, just in case
      
      read(30,3) rma,zma,psimag,psilim,beqd
      read(30,3) toteqd,psimx1,psimx2,xax1,xax2
      read(30,3) zax1,zax2,psisep,xsep,ysep
      read(30,3) (feqd(i),i=1,nveqd)
      read(30,3) (pres(i),i=1,nveqd)
      read(30,3) (ffpeqd(i),i=1,nveqd)
      read(30,3) (ppeqd(i),i=1,nveqd)
      read(30,3) ((peqd(i,j),i=1,nreqd),j=1,nzeqd)
c------------------------------------------------------------
c     creation of the poloidal flux peqd(i,j) with the minimum value
c     on the magnetic axis (psimag<psilim)
      write(*,*)' psimag, psilim in eqdsk :', psimag,psilim
      if(eqdsktype.eq.'TAE') then
        ! For some reason, psilim in TAE eqdsk is the value at grid edge, 
        ! far outside of LCFS.
        write(*,*)' WARNING: resetting psilim to 0.'
        write(*,*)' psilim should be the pol.flux at LCFS!'
        psilim=0.d0 !overwriting datafile to have psilim at LCFS, 
        !instead of grid edge! ESSENTIAL for TAE runs!
      endif
      
      write(*,*)'rma,zma,psimag,psilim,beqd'
      write(*,*) rma,zma,psimag,psilim,beqd
      if(psimag.gt.psilim) then
        dpsimax=-1
      else
        dpsimax=1
      endif
      psimag=dpsimax*psimag
      psilim=dpsimax*psilim
      write(*,*)'in equilib psimag,psilim',psimag,psilim
      do i=1,nreqd
        do j=1,nzeqd
	   peqd(i,j)=dpsimax*peqd(i,j) ! in equilib/input: reading_eqdsk
	enddo
      enddo
      write(*,*)'MIN/MAX of peqd(i,j) after mult. by dpsimax:',
     +  minval(peqd), maxval(peqd) 
c------------------------------------------------------------
      write(*,*)' eqdsktype=', eqdsktype
      !pause
      write(*,*)' If there are problems with reading eqdsk file'
      write(*,*)
     + ' try to adjust eqdsktype. Options: "TAE","tokamak","mirror" '
      if(eqdsktype.eq.'TAE') then
         ! For TAE (Tri Alpha Energy) FRC case, the eqdsk file 
         ! does not contain qpsi array, but instead it contains 
         ! necut and tmcut arrays (electron density [1/m^3] and
         ! el.temperature [eV]). So, one extra read block:
         read(30,3) ((necut(i,j),i=1,nreqd),j=1,nzeqd) !e.density [1/m^3]
         read(30,3) ((tmcut(i,j),i=1,nreqd),j=1,nzeqd) !Te [eV]
         write(*,*)'MIN/MAX of necut[1/m^3]',minval(necut),maxval(necut)
         write(*,*)'MIN/MAX of tmcut[eV]   ',minval(tmcut),maxval(tmcut)
         necut=necut*1.d-19 !n(R,Z) in units [1e19 m^-3] 
      else ! usual 'tokamak' eqdsk or a 'mirror' machine 
         read(30,3) (qpsi(i),i=1,nveqd)
      endif
      
      read(30,4) nnlim,nnves
      write(*,*)'in equilib after nnlim,nnves',nnlim,nnves
 4    format(2i5)
c-----------------------------
      nnlim=0 ! uncomment to ignore limiter points from eqdsk
      if (nnlim.ne.0) then
         allocate(rlimit(nnlim))
         allocate(zlimit(nnlim))
         read(30,3) (rlimit(i),zlimit(i),i=1,nnlim)
      else ! nnlim=0
c        ------------------------------------
c        The eqdsk data without limiter points.
c        ------------------------------------
         allocate(rlimit(nlimit)) ! these arrays can be in use:
         allocate(zlimit(nlimit)) ! set to LCFS. nlimit is in param.i
         rlimit=0.d0
         zlimit=0.d0  ! initialize
         nnlim=nlimit ! Reset for usage of these two arrays.
         write(*,*)'nnlim=0 eqdsk data without limiter points'
         psisep=psilim
         write(*,*)'psisep,psilim',psisep,psilim
      endif
      
      if (nnves.ne.0) then
         allocate(rves(nnves))
         allocate(zves(nnves))
         read(30,3) (rves(i),zves(i),i=1,nnves)
         write(*,*)' rves, zves:'
         do i=1,nnves
            write(*,*) rves(i),zves(i)
         enddo
      else ! nnves=0
         allocate(rves(nves)) ! nves is in param.i
         allocate(zves(nves)) ! Allocate just in case of table-data def.
         rves=0.d0
         zves=0.d0  ! initialize
         nnves=nves ! Reset for usage of these two arrays.
      endif
      
      close(30)
      write(*,*)' nreqd,nzeqd,nveqd'
      write(*,2) nreqd,nzeqd,nveqd

cc      write(*,*)'rdimeqd,zdimeqd,reqd,redeqd,zmideqd'
cc      write(*,3) rdimeqd,zdimeqd,reqd,redeqd,zmideqd
cc      write(*,*)'rma,zma,psimag,psilim,beqd'
cc      write(*,3) rma,zma,psimag,psilim,beqd
cc      write(*,*)'in equilib psimag,psilim',psimag,psilim
cc      write(*,*)'toteqd,psimx1,psixm2,xax1,xax2'
cc      write(*,3) toteqd,psimx1,psixm2,xax1,xax2
cc      write(*,*)'zax1,zax2,psisep,xsep,ysep'
cc      write(*,3) zax1,zax2,psisep,xsep,ysep
c      write(*,*) 'feqd(i),i=1,nveqd'
c      write(*,3) (feqd(i),i=1,nveqd)
c      write(*,*) 'pres(i),i=1,nveqd'
c      write(*,3) (pres(i),i=1,nveqd)
c      write(*,*) 'ffpeqd(i),i=1,nveqd'
c      write(*,3) (ffpeqd(i),i=1,nveqd)
c      write(*,*) 'ppeqd(i),i=1,nveqd'
c      write(*,3) (ppeqd(i),i=1,nveqd)

c      do i=1,nreqd
c         do j=1,nzeqd
c            write(*,*)'i,j,peqd(i,j)',i,j,peqd(i,j)
c         enddo
c      enddo
c      write(*,*) '(peqd(i,j),i=1,nreqd),j=1,nzeqd)'
cc      write(*,3) ((peqd(i,j),i=1,nreqd),j=1,nzeqd)
c      write(*,*) 'qpsi(i),i=1,nveqd'
c      write(*,3) (qpsi(i),i=1,nveqd)
c      write(*,*)'nnlim,nnves'
c      write(*,4) nnlim,nnves
cc      write(*,*) 'rlimit(i),zlimit(i),i=1,nnlim'
cc      write(*,3) (rlimit(i),zlimit(i),i=1,nnlim)
cc      write(*,*) 'rves(i),zves(i),i=1,nnves'
cc      write(*,3) (rves(i),zves(i),i=1,nnves)



      if (nveqd.ne.nreqd) then
c-----------------------------------------------------------
c     Interpolate feqd,pres,ffpeqd,ppeqd,qpsi to 
c     nreqd equispaced points.
c     The nveqd option has been enabled in CQL3D, and we
c     enable it also for GENRAY.  It is not a part of the
c     standard EQDSK  (BobH, 020722).
c
c     We use the standard, well documented spline package,
c     in zcunix.f, rather than use the uncommented subroutines
c     used in the majority of GENRAY cases (BobH, 020722).
c-----------------------------------------------------------

cAdding new arrays: psiar,d2feqd,d2pres,d2ffpeqd,d2ppeqd,d2qpsi,
c                   workk,r8temp,tabl,i1p,itabl
c     Create cubic spline arrays
c
         dpsiar=(psilim-psimag)/(nveqd-1)
         do i=1,nveqd
            psiar(i)=psimag+(i-1)*dpsiar
         enddo
         write(*,*) 'psiar(i),i=1,nveqd'
         write(*,3) (psiar(i),i=1,nveqd)
c     Use cubic splines at end points to obtain coefficients:
         i1p(1)=4 ! fitting cubic over first four points in psiar array
         i1p(2)=4 ! fitting cubic over last four points in psiar array
         call coeff1(nveqd,psiar,feqd,  d2feqd,  i1p,1,workk)
         call coeff1(nveqd,psiar,pres,  d2pres,  i1p,1,workk)
         call coeff1(nveqd,psiar,ffpeqd,d2ffpeqd,i1p,1,workk)!d2ffpeqd not used?
         call coeff1(nveqd,psiar,ppeqd, d2ppeqd, i1p,1,workk)!d2ppeqd not used?
         call coeff1(nveqd,psiar,qpsi,  d2qpsi,  i1p,1,workk)!d2qpsi not used?
         
c     Create nreqd long array and spline onto it.
         dflux=(psilim-psimag)/(nreqd-1)
         do i=1,nreqd
            flux(i)=psimag+(i-1)*dflux
         enddo
         
         itabl(1)=1
         itabl(2)=0
         itabl(3)=0
         do i=1,nreqd
            call terp1(nveqd,psiar,feqd,d2feqd,flux(i),1,tabl,itabl)
            r8temp(i)=tabl(1)
         enddo
         do i=1,nreqd
            feqd(i)=r8temp(i)
         enddo
         
         do i=1,nreqd
            call terp1(nveqd,psiar,pres,d2pres,flux(i),1,tabl,itabl)
            r8temp(i)=tabl(1)
         enddo
         do i=1,nreqd
            pres(i)=r8temp(i)
         enddo
         
         do i=1,nreqd
            call terp1(nveqd,psiar,ffpeqd,d2ffpeqd,flux(i),1,tabl,itabl)
            r8temp(i)=tabl(1)
         enddo
         do i=1,nreqd
            ffpeqd(i)=r8temp(i)
         enddo
         
         do i=1,nreqd
            call terp1(nveqd,psiar,ppeqd,d2ppeqd,flux(i),1,tabl,itabl)
            r8temp(i)=tabl(1)
         enddo
         do i=1,nreqd
            ppeqd(i)=r8temp(i)
         enddo
         
         do i=1,nreqd
            call terp1(nveqd,psiar,qpsi,d2qpsi,flux(i),1,tabl,itabl)
            r8temp(i)=tabl(1)
         enddo
         do i=1,nreqd
            qpsi(i)=r8temp(i)
         enddo
         
         write(*,*) 'flux(i),i=1,nreqd'
         write(*,3) (flux(i),i=1,nreqd)
         write(*,*) 'feqd(i),i=1,nreqd'
         write(*,3) (feqd(i),i=1,nreqd)
         write(*,*) 'pres(i),i=1,nreqd'
         write(*,3) (pres(i),i=1,nreqd)
         write(*,*) 'ffpeqd(i),i=1,nreqd'
         write(*,3) (ffpeqd(i),i=1,nreqd)
         write(*,*) 'ppeqd(i),i=1,nreqd'
         write(*,3) (ppeqd(i),i=1,nreqd)
         write(*,*) 'qpsi(i),i=1,nreqd'
         write(*,3) (qpsi(i),i=1,nreqd)
      endif ! nveqd.ne.nreqd
      
  12  continue ! handle for model_b.ne.0

      ! For any model_b:
      ! Overwrite wall_rmax, wall_zmax, etc, 
      ! to match r_wall(),z_wall() 
      if(n_wall.gt.0) then
         wall_rmin= minval(abs(r_wall))  ! r_wall(i) can be a negative
         wall_rmax= maxval(abs(r_wall))
         wall_zmin= minval(z_wall)
         wall_zmax= maxval(z_wall)
      endif

      if(model_b.ne.0) then
         redeqd= wall_rmin
         rdimeqd= (wall_rmax-wall_rmin) !R extension of cyl.chamber
         zdimeqd= (wall_zmax-wall_zmin) !Z-axial extension
         ! zmideqd is set to 0 in default_in
         zeqmax= wall_zmax
         zeqmin= wall_zmin
         xeqmax= wall_rmax ! max for x-grid of cartesian grid [meters]
         xeqmin=-xeqmax
         yeqmax= wall_rmax ! max for y-grid of cartesian grid [meters]
         yeqmin=-yeqmax  
!         if(model_b.eq.2) then ! B-field from external coils
!           redeqd= maxval(radc) ! set to largest r of all coils
!           wall_rmax=redeqd
!           rdimeqd= 2*redeqd !R extension 
!           xeqmax= redeqd ! max for x-grid of cartesian grid [meters]
!           xeqmin=-xeqmax
!           yeqmax= redeqd ! max for y-grid of cartesian grid [meters]
!           yeqmin=-yeqmax  
!         endif         
         if(model_rho_dens.eq.1  .or. model_rho_dens.eq.2) then 
           ! ensure that rho=1.1 surface is within the grid box
            xeqmax= max(xeqmax, elx0+elax*1.1)
            xeqmin= min(xeqmin, elx0-elax*1.1)
            yeqmax= max(yeqmax, ely0+elay*1.1)
            yeqmin= min(yeqmin, ely0-elay*1.1)
            zeqmax= max(zeqmax, elz0+elaz*1.1)
            zeqmin= min(zeqmin, elz0-elaz*1.1)
            zdimeqd= (zeqmax-zeqmin) !z-axial extension
            zmideqd=elz0 ! midplane
            zma=elz0 ! should not matter
            rma=0.d0 ! should not matter 
         endif
         if(model_rho_dens.eq.3) then 
           ! (x,y)-spline based on dengrid(i,j) data from file.
           ! Grids' limits: xdenmin, xdenmax, ydenmin, ydenmax
           ! Note: in this model, density is uniform in z-direction.
           ! Set boundaries within dengrid grid limits:
            xeqmax= min(xeqmax, xdenmax) 
            xeqmin= max(xeqmin, xdenmin)
            yeqmax= min(yeqmax, ydenmax)
            yeqmin= max(yeqmin, ydenmin)
         endif
         if(model_rho_dens.eq.4) then 
            zmideqd=0.d0 ! midplane
            zma=0.d0 
            rma=rs_frc/sqrt(2.) ! Magnetic Axis
         endif
         dpsimax=1 ! +1 means that psilim>psimag
         psimag=0.d0
         
         if(model_b.eq.1) then
            ! Uniform Bz, Br=0, Bphi~1/r
            ! wall_*** values are from namelist /tokamak/
            rma=0.d0 ! model_b=1  no magnetic axis 
            zma=0.d0 ! model_b=1  
            psimag=0.d0
            psilim=0.d0 ! will be determined below.
            rlim= rlim_wall_fraction*
     +            min(xeqmax,yeqmax,abs(xeqmin),abs(yeqmin))
            psilim=0.5*bz0*rlim**2
         endif
         if(model_b.eq.2) then ! B-field from external coils
            y=0.d0 ! phi=0 (axial symmetry is assumed)
            z=0.d0 ! midplane
            rma=0.d0 ! model_b=2  no magnetic axis 
            zma=0.d0 ! model_b=2  
            psimag=0.d0
            psilim=0.d0 ! will be determined below.
            rlim= rlim_wall_fraction*wall_rmax
            write(*,*)'equilib/input/model_b=2:  rlim=',rlim
            ! rlim_wall_fraction is in namelist 
            ! (set it in genray.in)
            ! It is used to define the value of rlim as a fraction
            ! of chamber wall radius:
            !    rlim= rlim_wall_fraction*wall_rmax
            ! rlim corresponds to psilim and to "rho"=1.
            ! Although for model_b=1 or 2 there are no flux surfaces,
            ! we still need to define the edge of plasma, 
            ! for definition of plasma profiles n and T.  
            ! This is done by setting the value of rlim; 
            ! the value of psilim is found by integration 
            ! of poloidal (Z) flux over radius R (at fixed Z=0 coord.)
            ! psilim corresponds to "rho" such that rho(psilim)=1.
            ! So, the value of rlim sets the edge of hot plasma.
            ! The profiles of n and T are set to drop exponentially 
            ! or linearly at psi>psilim ("rho">1).
            dra= rlim/(nreqd-1)
            ! Find psilim by integration (no need for accuracy)
            do i=2,nreqd
               x= (i-1)*dra ! [0+dra; rlim]
               call bfield_coils(x,y,z, b_x, b_y, b_z) !-> get b_z
               psilim= psilim+ x*b_z*dra ! integral
            enddo
         endif ! model_b=2
         if(model_b.eq.3) then !mirror machine ("mirror1" model)
            rma=0.d0 ! model_b=3  no magnetic axis 
            zma=0.d0 ! model_b=3  or, formally, at (R,Z)=(0,0)
            psimag=0.d0 ! pol flux at (R,Z)=(0,0)
            !Next: define glb_mirror and psilim. 
            ! 1. Define variable glb_mirror (the scale length, [m])
            ! used in the model equations:
            !     Bz(R,Z)=  B00*J0(R/glb)*cosh(Z/glb)
            !     Br(R,Z)= -B00*J1(R/glb)*sinh(Z/glb)
            glb_mirror=0.5*zbox_mirror/acosh(rmirror) 
            ! Note: zbox_mirror=2*Zmax.
            ! set rmirror in genray.in; default is 2.
            ! set zbox_mirror in genray.in: 
            ! effectively zbox_mirror defines the location of 
            ! magnetic throats at Z= +zbox_mirror/2 and -zbox_mirror/2.
            ! 2. Define psilim.
            ! Select proper coordinate (x,y,z) for psilim:
            y=0.d0 ! phi=0 (axial symmetry is assumed: select any y)
            z=zma  ! the midplane
            rlim= rbox_mirror ! this can be used for x coord (when y=0)
            ! rlim corresponds to psilim and to "rho"=1.
            ! Although for model_b=1,2,3 there is no "LCFS" ,
            ! we still need to define the edge of plasma, 
            ! for definition of plasma profiles n and T.  
            ! This is done by setting the value of rlim; 
            ! So, the value of rlim sets the edge of hot plasma.
            ! The profiles of n and T are set to drop exponentially 
            ! or linearly at psi>psilim ("rho">1).
            ! But check that rbox_mirror is smaller than R_b0
            ! where R_b0 is such that J0(R_b0/glb)=0.
            ! From the model equations, Bz(R,Z=0) becomes 0 
            ! at R/glb= R_b0/glb = 2.4048, the 1st null of J0(R/glb),
            R_b0= glb_mirror*2.4048
            ! This is also the point where poloidal flux PSI(R) 
            ! reaches max value (when going along Z=0 line). 
            ! The point (R_b0,Z=0) is a saddle point
            ! for PSI in (R,Z) plane. 
            ! The magnetic field line in (R,Z) plane, defined as
            ! PSI(R,Z)=const, becomes "broken" at (R_b0,Z=0),
            ! i.e. it cannot follow from Z=-0.5*zbox to Z=+0.5*zbox region.
            ! Check and reset the value of rbox_mirror:
            if(rbox_mirror.gt.0.95*R_b0)then
               rbox_mirror= 0.95*R_b0 ! "0.95", to leave some space 
               ! for definition of "last surface not hitting chamber"
               ! that could be at R=0.98*R_b0, for example.
               WRITE(*,*)
            WRITE(*,*)'equilib: RESETTING rbox_mirror to ', rbox_mirror
            WRITE(*,*)'equilib: Limitation of mirror1(model_b=3) model:'
               WRITE(*,*)'equilib: rbox_mirror cannot exceed R_b0 .'
               WRITE(*,*)'equilib: See description of model_b=3 .'
               WRITE(*,*)
               rlim= rbox_mirror ! reset, too
            endif
            x=rlim
            call eq_mirror1(x,y,z, PSI, b_x, b_y, b_z) !-> PSI and B
            psilim= PSI ! actual psilim, not yet adjusted for the code.
            write(*,*)'equilib/model_b=3: rlim,psilim,glb_mirror=',
     +                rlim,psilim,glb_mirror
         endif ! model_b=3
         if(model_b.eq.4) then  ! model_b=4 (FRC-like plasma)
            ! Bz= bz0*tanh(Kappa*[2*(r/rs)^2 -1]), Br=0, Bphi~1/r
            ! bz0= magnetic field at r=+INF, approximately at wall [T]
            beta_frc= 1.d0 - 0.5*(rs_frc/wall_rmax)**2  !-> to one.i
            ! Find "kappa" such that  beta_frc*kappa-tanh(kappa)=0
            akappa= rtbis(func_kappa, 1.d-5, 10.d0, 1.d-14) !-> to one.i
            ! Last three arguments: initial guesses that bracket the root
            ! and the accuracy for searching the root.
            ! Now scan along x, calculate pol.flux psilim at r=rlim
            !psilim can be used to define "rho" such that rho(psilim)=1.
            rlim= rs_frc ! == Separatrix radius.
            rma=  rs_frc/sqrt(2.) ! model_b=4  magn.axis is at rs/sqrt(2)
            zma=0.d0 ! model_b=4  
            psilim=0.d0 ! will be determined below, at r=rlim=rs_frc
            dra= rlim/(nreqd-1)
            y=0.d0 ! phi=0 (axial symmetry is assumed)
            z=zma  ! midplane
            ! Find psilim by integration (no need for accuracy)
            do i=2,nreqd
               x= (i-1)*dra ! [0+dra; rlim]
               call bfield_frc(x,y,z, b_x, b_y, b_z) !-> get b_z
               psilim= psilim+ x*b_z*dra ! integral from r=0 to r=rlim
            enddo
            psimag=0.d0 ! will be determined below, at r=rma
            dra= rma/(nreqd-1)
            y=0.d0 ! phi=0 (axial symmetry is assumed)
            z=zma  ! midplane
            ! Find psimag by integration (no need for accuracy)
            do i=2,nreqd
               x= (i-1)*dra ! [0+dra; rma]
               call bfield_frc(x,y,z, b_x, b_y, b_z) !-> get b_z
               psimag= psimag+ x*b_z*dra ! integral from r=0 to r=rma
            enddo
         endif ! model_b=4 (FRC-like plasma)

         if(psimag.gt.psilim) then ! In the code, peqd() array should be
            dpsimax=-1             ! an ascending function of rho.
         else                      ! The physical sign of PSI is tracked
            dpsimax=1              ! by sign of dpsimax .
         endif
         ! Negative dpsimax corresponds to Bz<0 (R*Bz=dPSI/dR)
         psimag=dpsimax*psimag !Adjusted for the code: Now psilim>psimag
         psilim=dpsimax*psilim !Adjusted: Now psilim>psimag (model_b.ne.0)
      endif ! model_b .ne. 0
  
c---------------------------------------------------------
c  normalization of the eqdsk data
c--------------------------------------------------------
      pr=1.d0/r0x
      pb=1.d0/b0
      rdimeqd=rdimeqd*pr
      zdimeqd=zdimeqd*pr
      reqd=reqd*pr
      redeqd=redeqd*pr
      zmideqd=zmideqd*pr
      xeqmin=xeqmin*pr
      xeqmax=xeqmax*pr
      yeqmin=yeqmin*pr
      yeqmax=yeqmax*pr
      zeqmin=zeqmin*pr
      zeqmax=zeqmax*pr
      rma=rma*pr
      zma=zma*pr
      psimag=psimag*pr*pr*pb
      psilim=psilim*pr*pr*pb
      beqd=beqd*pb
      psimx1=psimx1*pr*pr*pb
      psimx2=psimx2*pr*pr*pb
      xax1=xax1*pb
      xax2=xax2*pb
      zax1=zax1*pr
      zax2=zax2*pr
      psisep=psisep*pr*pr*pb
      xsep=xsep*pr
      ysep=ysep*pr
      write(*,*)'in equilib psimag,psilim',psimag,psilim
      write(*,*)'nveqd,nreqda',nveqd,nreqda
      
c----- GRIDS ---------------------------------------------
      dreq= rdimeqd/(nreqd-1)
      write(*,*) 'rdimeqd, redeqd=', rdimeqd, redeqd
      do 100 i=1,nreqd
        req(i)= redeqd+dreq*(i-1)
 100  continue
 
c      dstep=zdimeqd/(nzeqd-1)
c      do 200 i=1,nzeqd
c        zeq(i)= zmideqd+dstep*(i-1) -zdimeqd*0.5d0
c 200  continue

      dstep= (zeqmax-zeqmin)/(nzeqd-1)
      do i=1,nzeqd
        zeq(i)= zeqmin+dstep*(i-1) ! YuP [zeqmin; zeqmax]
      enddo

      dstep= (xeqmax-xeqmin)/(nxeqd-1)
      do i=1,nxeqd
        xeq(i)= xeqmin+dstep*(i-1) ! YuP [xeqmin; xeqmax]
      enddo
 
      dstep= (yeqmax-yeqmin)/(nyeqd-1)
      do i=1,nyeqd
        yeq(i)= yeqmin+dstep*(i-1) ! YuP [yeqmin; yeqmax]
      enddo

      if(model_b.eq.1) then ! Uniform Bz=bz0, Br=0.
        ! Define peqd(i,j)  == psi/2pi
        do  j=1,nzeqd
        do  i=1,nreqd
          peqd(i,j)= 0.5d0*bz0*req(i)**2 ! no dep. on z 
        enddo
        enddo
      endif

      if(model_b.eq.2) then ! Field from coils (current loops)
         ! Define peqd(i,j)  == psi/2pi
         y=0.d0
         do j=1,nzeqd
            z=zeq(j)
            ! Find psi by integration (no need for accuracy)
            sumpsi=0.d0
            dra= req(2)-req(1)
            do i=1,nreqd
               x= req(i) 
               call bfield_coils(x,y,z, b_x, b_y, b_z) !-> get b_z
               sumpsi= sumpsi+ x*b_z*dra ! integral up to req(i)
               peqd(i,j)= sumpsi  !for model_b=2 (coils)
            enddo
         enddo
         do j=1,nzeqd
         do i=1,nreqd
            peqd(i,j)=dpsimax*peqd(i,j) !for model_b=2 (coils)
         enddo
         enddo
      endif ! model_b=2

      if(model_b.eq.3) then !mirror machine ("mirror1" model)
         ! Define peqd(i,j)==psi/2pi (should be asc.func. of rho; see dpsimax)
         y=0.d0 ! x will be scanned, as R.
         do j=1,nzeqd
            z=zeq(j)
            do i=1,nreqd
               x= req(i) 
               call eq_mirror1(x,y,z, PSI, b_x, b_y, b_z) !-> PSI and B
               peqd(i,j)=dpsimax*PSI  !for model_b=3 (mirror1)
               !Mult. by dpsimax to be consist. with psilim>psimag for the code
            enddo
         enddo
         !Note: req() grid goes beyond R_b0
         ! where R_b0 is such that J0(R_b0/glb)=0,
         ! and PSI(R,Z) has a saddle point.
         !When going along Z=0, PSI(R) grows up to the R=R_b0 line, 
         ! then goes down; we want to avoid such values of PSI. 
         !To avoid problems with T(rho) and n(rho) profiles,
         ! consider modifying peqd() outside of R_b0.
         R_b0= glb_mirror*2.4048
         ir_b0= INT( (R_b0-redeqd)/dreq ) ! lower-nearest ir index
         if(ir_b0.lt.nreqd)then
            do j=1,nzeqd
               z=zeq(j)
               !do i= 1,ir_b0 ! At R<R_b0
               !   peqd(i,j)=peqd(i,j)    No changes at R<R_b0
               !enddo
               call eq_mirror1(R_b0,y,z, PSI_b0, b_x, b_y, b_z) !-> PSI
               do i= ir_b0+1, nreqd ! At R>R_b0 (and any Z)
                  peqd(i,j)=dpsimax*PSI_b0  !or simple: peqd(ir_b0,j)
               enddo
               ! Now peqd is a growing function from R=0 up to the line 
               ! R=R_b0 (at any Z), and flat outside.
            enddo
         endif
      endif ! model_b=3

      if(model_b.eq.4) then  ! model_b=4 (FRC-like plasma)
         ! Bz= bz0*tanh(Kappa*[2*(r/rs)^2 -1]), Br=0, Bphi~1/r
         ! bz0= magnetic field at r=+INF, approximately at wall [T]
         ! Define peqd(i,j)  == psi/2pi
         y=0.d0 ! phi=0 (axial symmetry is assumed)
         do j=1,nzeqd
            z=zeq(j)
            ! Find psi by integration (no need for accuracy)
            sumpsi=0.d0
            dra= req(2)-req(1)
            do i=1,nreqd
               x= req(i) ! r=[0;wall_rmax]
               call bfield_frc(x,y,z, b_x, b_y, b_z) !-> get b_z
               sumpsi= sumpsi+ x*b_z*dra ! integral up to req(i)
               peqd(i,j)= sumpsi !for model_b=4 (FRC model)
               !if(j.eq.1)write(*,'(3e13.4)')x,b_z,sumpsi
            enddo
         enddo
         do j=1,nzeqd
         do i=1,nreqd
            peqd(i,j)=dpsimax*peqd(i,j) !for model_b=4 (FRC model)
         enddo
         enddo
      endif ! model_b=4 (FRC-like plasma)

      if(model_b.ne.0) then
         dreq=req(2)-req(1)
         ! Define feqd(i)  == R*Bphi
         do i=1,nreqd
            r= req(i) 
            if(r.gt.dreq) then
               feqd(i)= rbphi0 ! Bphi~1/r
            else ! limit Bphi near r=0; assume Bphi~r
               feqd(i)= rbphi0*(r/dreq)**2  ! Note: req(1) is zero
            endif
         enddo
      endif ! model_b>0 (analytical models)
      

      do 30 i=1,nveqd 
        feqd(i)=feqd(i)*pr*pb
        pres(i)=pres(i)*pr*pr*pb
 30   continue
      do 31 i=1,nreqd
      do 32 j=1, nzeqd
 32      peqd(i,j)=peqd(i,j)*pr*pr*pb
 31   continue
      if (nnlim.gt.0) then
         do 33 i=1,nnlim
         rlimit(i)=rlimit(i)*pr
 33	   zlimit(i)=zlimit(i)*pr
      endif
c------------------------------------------------------------------
      dpsi=(psilim-psimag)/(nreqd-1) ! positive
      do 15 i=1,nreqd
        flux(i)=psimag+dpsi*(i-1) ! ascending (growth)
15    continue

c-------------------------------------------------------------
c     test map psi(i,j) on (R,Z) plane
c------------------------------------
c      open(40,file='psirz.dat')
c41     format(3(' ',e16.9))
c      write(40,*)' r z psi'         
c      do i=1,nreqd
c         do j=1,nzeqd
c         write(40,41)req(i),zeq(j) ,peqd(i,j)
c         enddo
c      enddo
c      close(40)
c      stop
c     end test map
c---------------------------------------------------

      return ! input
      end




c*************************DINITR***************************************
c  creation of arrays for eqdsk spline coefficients		      *
c  for psieqd:  tx(nr4a),ty(nz4a),cxy(nr4a,nz4a)			      *
c  for feqd  :  txf(nr4a),cx(nr4a)				      *
c  for pres  :  tpres(nr4a),cpres(nr4a)	        		      *
c  creation of arrays for limiter spline coefficients		      *
c  for zlimmiter plus(r)  : trlimp(nlim4),cxlimp(nlim4)		      *
c  for zlimmiter minus(r) : trlimm(nlim4),cxlimm(nlim4)		      *
c								      *
c  this program uses the following subroutines:  		      *
c								      *
c  limitr(zzp,rrp,zzm,rrm,ip,im,				      *
c         rmax,zrmax,rmin,zrmin,zmax,rzmax,zmin,rzmin)		      *
c  for creation arrays:{zzp(ip),rrp(ip)} for zlimmiter upper (r),     *
c                      {zzm(im),rrm(im)} for zlimmiter under (r)      *
c  from zlimit,rlimit						      *
c								      *
c  iac1r(rlimr,ip,ip4,zlimr,lx,mxa,fxa,mxb,fxb,trlimp,cxlimp)	      *
c  -calculates the spline coefficients for 1d function		      *
c								      *
c  IAC2R(rrr,nx,zzr,ny,peqdr,nfx,nfy,lx,ly,mxa,arfxa,mxb,arfxb,	      *
c  mya,arfya,myb,arfyb,tx,ty,cxy,ncx,ncy,nry,cy)		      *
c  - calculates the spline coefficients for 2d function		      *
c**********************************************************************

      subroutine dinitr
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'fourb.i'
      include 'five.i' ! rmax, etc.
      include 'gr.i' ! zpsi, rpsi (npsi,nteta1)
      real*8,allocatable :: zzp(:),rrp(:),zzm(:),rrm(:) ! (nnlim) (was nlimit)
      real*8,allocatable :: zlimr(:),rlimr(:) ! (nnlim) (was nlimit)
      double precision
     1     arfxa(nzeqda),arfxb(nzeqda),arfya(nreqda),arfyb(nreqda),
     1	   fxa,fxb,fr(nreqda),fluxr(nreqda),rrr(nreqda),zzr(nzeqda),
     2	   peqdr(nreqda,nzeqda)
      double precision	ias1r,ias2r,ias2r_Sm
      character*12 namfil
      
      IF (.NOT. ALLOCATED(zzp)) then
      allocate(zzp(nnlim))
      allocate(rrp(nnlim))
      allocate(zzm(nnlim))
      allocate(rrm(nnlim))
      allocate(zlimr(nnlim))
      allocate(rlimr(nnlim))
      endif
      
      if(nteta1.lt.nnlim) then
       write(*,*)'dinitr: nnlim,nteta+1=', nnlim,nteta1
       write(*,*)'dinitr: increase nteta in param.i to make nteta>nnlim'
       stop
      endif
      if(nlim4.lt.nnlim) then
       write(*,*)'dinitr: nnlim,nlimit+4=', nnlim,nlim4
       write(*,*)'increase nlimit in param.i to make nlimit+4 > nnlim'
       stop
      endif
c--------------------------------------------------------------------
c     determination of rmax,rmin,zmax,zmin
c--------------------------------------------------------------------      
      if(model_b.ne.0) then
         ! wall_*** values are from namelist /tokamak/
         rmax= wall_rmax
         rmin= abs(wall_rmin)
         zmax= wall_zmax
         zmin= wall_zmin
         if(model_rho_dens.eq.1  .or. model_rho_dens.eq.2) then 
           ! ensure that rho_density=1.1 surface is within the grid box
            rmax= max(rmax, elx0+elax*1.1)
            rmax= max(rmax, ely0+elay*1.1)
            zmax= max(zmax, elz0+elaz*1.1)
            zmin= min(zmin, elz0-elaz*1.1)
         endif
         ! Check that it is within equilib.grid
         rmax= min(rmax,xeqmax,yeqmax,abs(xeqmin),abs(yeqmin))
         rzmax=rmax
         rzmin=rmin
         zrmax=zmax
         zrmin=zmin
      endif
      
      if(model_b.eq.0) then ! when using eqdsk data
         if (maxval(rlimit).gt.0.d0)  then
        write(*,*)'before limitr'
            call limitr(nnlim,zzp,rrp,zzm,rrm,ip,im,
     1           rmax,zrmax,rmin,zrmin,zmax,rzmax,zmin,rzmin)
        write(*,*)'after limitr'
         else
            call makelacn(zmax1,rmax1,zmin1,rmin1,
     1           rzmax1,zrmax1,rzmin1,zrmin1)
            rmax=rmax1
            zmax=zmax1
            rmin=rmin1
            zmin=zmin1
            zrmax=zrmax1
            rzmax=rzmax1
            zrmin=zrmin1
            rzmin=rzmin1
         endif
      endif
c--------------------------------------------------------------------
c     spline coefficients for psi, feqd and pres creation 
c--------------------------------------------------------------------
      nx=nreqd
      ny=nzeqd
      nr4=nx+4
      nz4=ny+4
      lx=1
      mxa=0
      fxa=0.d0
      mxb=1
      fxb=0.d0
      do 1111 i=1,nx
        fluxr(i)=flux(i)
        fr(i)=feqd(i) ! r*B_phi  at nreqd equispaced points in psi
 1111	rrr(i)=req(i)
      do 1112 j=1,ny
 	zzr(j)=zeq(j)
 1112 continue

      do 1113 i=1,nx
      do 1113 j=1,ny
 1113 peqdr(i,j)=peqd(i,j)
      			    
      call iac1r(fluxr,nx,nr4,fr,lx,mxa,fxa,mxb,fxb,txf,cx) !fr=feqd
c-------------------------------------------------------------------	
      call iac1r(fluxr,nx,nr4,pres,lx,mxa,fxa,mxb,fxb,tpres,cpres) !pres
c-------------------------------------------------------------------
      ny=nzeqd
      lx=1
      ly=1
      mxa=0
      mxb=0
      if(eqdsktype.eq.'mirror')then !YuP[12-2016] enforce dPSI/dR=0 
         !at R=0 and all Z.
         mxa=1 ! mxa=1 means: First derivative is given at rrr(1),
         !which is R=0 in a mirror machine.
         arfxa(1:nzeqda)=0.d0 ! Derivative=0 at all Z
      endif
      mya=0  ! yup: presumably 3 means periodic conditions (?)
      myb=0
      do 92 i=1,ny !==nzeqd
        arfxa(i)=0.
 92     arfxb(i)=0.
      do 94 i=1,nx !==nreqd
        arfya(i)=0.
 94     arfyb(i)=0.
      nfx=nx
      nfy=ny
      ncx=nx+4
      ncy=ny+4
      nry=max(nx,ny)+4
      call iac2r_Sm(rrr,nx,zzr,ny,peqdr,nfx,nfy,lx,ly,mxa,arfxa,mxb,
     &                 arfxb,
     &                 mya,arfya,myb,arfyb,tx,ty,cxy,ncx,ncy,nry,cy,
     &                 nreqda,nr4a) ! peqd  (==psi, poloidal_flux/2pi)
      write(*,*)
     1     'spline coefficients for psi, feqd, and pres were created'

      if(model_b.ne.0) then
         return
      endif
ctest     
      cnorm=0.d0     
      do i=1,nx
         do j=1,ny
           psi_t=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,0,0,rrr(i),zzr(j),
     &     nr4a)
           cnorm=cnorm+(peqdr(i,j)-psi_t)**2
         enddo
      enddo
      cnorm=dsqrt(cnorm)/dfloat(nx*ny)
      write(*,*)'equilib cnorm of psi',cnorm
cendtest

c--------------------------------------------------------------------
c  The calculations of the coordinates for flux functions
c          r(psi,teta)   z(psi,teta)
c          arrays rpsi(j,i) zpsi(j,i)
c          j=1,npsi(number of counturs poloidal flux=constant)
c          i=1,nteta+1(number of point in poloidal angle)
c  and creates arrays ar(nl,nteta+1),az(nl,nteta+1) for nl surfaces
c  Used for spline (see rhospl) and plotting
      if(model_b.eq.0) then ! when using eqdsk data
         call gr2new
      endif
c--------------------------------------------------------------------
c     creation of spline coefficients for limiter
c--------------------------------------------------------------------
      if (maxval(rlimit).eq.0.d0) then
c        ------------------------------------------------
c        creation of arrays zlimit(), rlimit()
c        They set equal to rpsi(npsi,i) zpsi(npsi,i) (the
c        coordinates of the magnetic surface psi(r,s)=psilim.
c        In this case it must be nlimitr=nteta1
c      	 -------------------------------------------------
c        write(*,*)'nnlim,npsi',nnlim,npsi
	 do i=1,nnlim !nlimit
	   zlimit(i)=zpsi(npsi,i)
	   rlimit(i)=rpsi(npsi,i)
c	   write(*,*)'i,zlimit(i),rlimit(i)',i,zlimit(i),rlimit(i)
	 enddo
         call limitr(nnlim,zzp,rrp,zzm,rrm,ip,im,
     1   rmax,zrmax,rmin,zrmin,zmax,rzmax,zmin,rzmin)
c         write(*,*)'in equilib/dinitr  after call limitr 
c     1   ip,im,rmax,zrmax,rmin,zrmin,zmax,rzmax,zmin,rzmin',
c     1   ip,im,rmax,zrmax,rmin,zrmin,zmax,rzmax,zmin,rzmin
      endif
      ip4=ip+4
      lx=1
      mxa=0
      fxa=0.0d0
      mxb=0
      fxb=0.0d0
      if (rrp(1).lt.rrp(ip)) then
        do 111 i=1,ip
          zlimr(i)=zzp(i)
 111      rlimr(i)=rrp(i)
      else
        do 2112 i=1,ip
          zlimr(i)=zzp(ip-i+1)
 2112     rlimr(i)=rrp(ip-i+1)
      end if
 17   format(2e16.9)
      call iac1r(rlimr,ip,ip4,zlimr,lx,mxa,fxa,mxb,fxb,trlimp,cxlimp)
      im4=im+4
      if (rrm(1).lt.rrm(im)) then
        do 112 i=1,im
          zlimr(i)=zzm(i)
 112      rlimr(i)=rrm(i)
      else
        do 1122 i=1,im
          zlimr(i)=zzm(im-i+1)
 1122     rlimr(i)=rrm(im-i+1)
      end if
c      write(*,*)'in equilib.for rlimr,zlimr,i=1,im'
c      read(*,*)
c      write(*,17)(rlimr(i),zlimr(i),i=1,im)
      call iac1r(rlimr,im,im4,zlimr,lx,mxa,fxa,mxb,fxb,trlimm,cxlimm)

        return
        end





c*************LIMITR**************************************************
c   this subroutine creates the arrays {zzp(ip),rrp(ip)},
c   {zzm(im),rrp(im)} for function zlimit upper(r) and zlimit under(r)
c   and determinates the coordinates of:
c      top limitter point  (zmax,rzmax),
c      botom limitter point(zmin,rzmin),
c      inner limitter point(zrmin,rmin),
c      outer limitter point(zrmax,rmax)
c**********************************************************************
c   input data are in common 'three' and 'fourb'
c   output data are all parameters of the limitr()
c----------------------------------------------------------------------
      subroutine limitr(nlim,zzp,rrp,zzm,rrm,ip,im,
     1 rmax,zrmax,rmin,zrmin,zmax,rzmax,zmin,rzmin)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'three.i'
      include 'fourb.i'
      include 'limit.i'
      real*8 zzp(nlim),rrp(nlim),zzm(nlim),rrm(nlim) ! (nnlim) (was nlimit)
      real*8 zlimr(nlim),rlimr(nlim) ! (nnlim) (was nlimit)

c------------------------------------------------------------------
c   determination rmin and rmax  of limitter and
c   coordinates of these points: (zrmax,rmax),(zrmin,rmin)
c------------------------------------------------------------------
      irmax=1
      irmin=1
      rmax=rlimit(1)
      rmin=rlimit(1)
      
      do 10 i=2,nnlim  ! was nlimit
        if (rmax.lt.rlimit(i)) then
          rmax=rlimit(i)
          irmax=i
          goto 10  !YuP  Why?
        end if

        if (rmin.gt.rlimit(i)) then
          rmin=rlimit(i)
          irmin=i
        end if
 10   continue
      zrmin=zlimit(irmin)
      zrmax=zlimit(irmax)
c        write(*,*)'irmin=',irmin,'rmin=',rmin,'zrmin=',zrmin
c        write(*,*)'irmax=',irmax,'rmax=',rmax,'zrmax=',zrmax
c---------------------------------------------------------------------
c  begin of arrays rrp,zzp  and  rrm,zzm  creation
c---------------------------------------------------------------------
       if (irmin.gt.irmax) then
         do 20 i=1,nnlim-irmin
	     rrp(i)=rlimit(irmin+i-1)
 	     zzp(i)=zlimit(irmin+i-1)
 20      continue

         do 21 i=nnlim-irmin+1,nnlim-1
	     rrp(i)=rlimit(i-nnlim+irmin)
 	     zzp(i)=zlimit(i-nnlim+irmin)
 21      continue
       else
c--------- if irmin.le.irmax then--------------------------
         do 22 i=1,nnlim-irmin
	     rrp(i)=rlimit(i+irmin-1)
 	     zzp(i)=zlimit(i+irmin-1)
 22      continue
         do 23 i=nnlim-irmin+1,nnlim-1
	     rrp(i)=rlimit(i-nnlim+irmin)
 	     zzp(i)=zlimit(i-nnlim+irmin)
 23      continue
       end if
c------------------------------------------------------------------
       do 24 i=1, nnlim-1
         rlimit(i)=rrp(i)
         zlimit(i)=zzp(i)
 24    continue
c------------------------------------------------------------------
c   determination rmin and rmax ,zmin and zmax of limitter and
c   coordinates of these points:(zrmax,rmax),(zrmin,rmin),
c   (zmax,rzmax),(zmin,rzmin)
c------------------------------------------------------------------
	irmax=1
	irmin=1
	izmax=1
	izmin=1
	rmax=rlimit(1)
	rmin=rlimit(1)
	zmax=zlimit(1)
	zmin=zlimit(1)

	do 30 i=1,nnlim-1

	  if (rmax.lt.rlimit(i)) then
	    rmax=rlimit(i)
	    irmax=i
	    zrmax=zlimit(irmax)
	  end if

	  if (rmin.gt.rlimit(i)) then
	    rmin=rlimit(i)
	    irmin=i
	    zrmin=zlimit(irmin)
	  end if

	  if (zmax.lt.zlimit(i)) then
	    zmax=zlimit(i)
	    izmax=i
	    rzmax=rlimit(izmax)
	  end if

	  if (zmin.gt.zlimit(i)) then
	    zmin=zlimit(i)
	    izmin=i
	    rzmin=rlimit(izmin)
	  end if
 30    continue

c        write(*,*)'irmin=',irmin,'rmin=',rmin,'zrmin=',zrmin
c        write(*,*)'irmax=',irmax,'rmax=',rmax,'zrmax=',zrmax
c        write(*,*)'izmin=',izmin,'zmin=',zmin,'rzmin=',rzmin
c        write(*,*)'izmax=',izmax,'zmax=',zmax,'rzmax=',rzmax
c------------------------------------------------------------------
c       rrp(1) must be less then rrp(ip)
c       rrm(1) must be less then rrm(im)
c       following part of program creates arrays {rrp,zzp} and
c       {rrm,zzm}  in  which rrp(1).lt.rrp(ip),
c                            rrm(1).lt.rrm(im)
c------------------------------------------------------------------
	 if ((izmax.ge.1).and.(izmax.le.irmax)) then
c-----------------------------------------------------------------
c rrp(irmax),zzp(irmax),rrm(nnlim-irmax+1),zzm(nnlim-irmax+1)
c -----------------------------------------------------------------

	   ip=irmax
	   im=nnlim-irmax+1
	   do 40 i=1,ip
	     rrp(i)=rlimit(i)
  40	     zzp(i)=zlimit(i)
           do 41 i=1,im-1
	     rrm(i)=rlimit(i+irmax-1)
  41	     zzm(i)=zlimit(i+irmax-1)
             rrm(im)=rrp(1)
	     zzm(im)=zzp(1)
	 else
c
c rrm(irmax),zzm(irmax),rrp(nnlim-irmax+1),zzp(nnlim-irmax+1)
c
	   im=irmax
	   ip=nnlim-irmax+1
	   do 42 i=1,im
	     rrm(i)=rlimit(i)
  42	     zzm(i)=zlimit(i)
           do 43 i=1,ip-1
	     rrp(i)=rlimit(i+irmax-1)
  43	     zzp(i)=zlimit(i+irmax-1)
             rrp(ip)=rrm(1)
	     zzp(ip)=zzm(1)
	 end if
c-----------------------------------------------------------------
      if (rrp(1).gt.rrp(ip)) then
        do i=1,ip
          zlimr(i)=zzp(i)
          rlimr(i)=rrp(i)
	enddo
        do i=1,ip
          zzp(i)=zlimr(ip-i+1)
          rrp(i)=rlimr(ip-i+1)
	enddo
      end if
      if (rrm(1).gt.rrm(im)) then
        do i=1,im
          zlimr(i)=zzm(i)
          rlimr(i)=rrm(i)
	enddo
        do i=1,im
          zzm(i)=zlimr(im-i+1)
          rrm(i)=rlimr(im-i+1)
	enddo
      end if
c-----------------------------------------------------------------
c  creation of the monotonic limiter array zzp  near rmin and rmax
      izmax=1
      zmax=zzp(1)
      do i=2,ip
        if (zmax.lt.zzp(i)) then
          zmax=zzp(i)
          izmax=i
	endif
      enddo
c      write(*,*)'ip,izmax,zzp(izmax),rrp(izmax)',
c     1 ip,izmax,zzp(izmax),rrp(izmax)

c--------near rmin --------------------------------
c     the correction for the second point (begin)
c      write(*,*)'before correction rrp(2),rrp(1)',rrp(2),rrp(1)
c      write(*,*)'before correction zzp(2),zzp(1)',zzp(2),zzp(1)
      if (rrp(2).lt.rrp(1)+1.d-5) then
         do j=3,izmax
            if (rrp(j).gt.rrp(2)) then
	       j0=j
	       goto 45
	    endif
	 enddo
 45      continue
         if(zzp(j0).eq.zzp(1)) then
	    rrp(2)=rrp(1)+1.d-5
	 else
cSmirnov970104 beg
	    if(dabs(rrp(1)-rrp(2)).lt.1.d-6) then
	       rrp(2)=0.5d0*(rrp(j0)+rrp(1))
	    else
	    rrp(2)=rrp(1)+(zzp(2)-zzp(1))*(rrp(j0)-rrp(1))/
     1	       	      (zzp(j0)-zzp(1))
	    endif
cSmirnov961226 end
	 endif
      endif
c      write(*,*)'after correction rrp(2),rrp(1)',rrp(2),rrp(1)
c      write(*,*)'after correction zzp(2),zzp(1)',zzp(2),zzp(1)
c     the correction for the second point (end)
c     --------------------------------------
      do i=3,izmax-1
c         if(rrp(i).le.rrp(i-1)) then
         if(rrp(i).le.(rrp(i-1)+5.d-3)) then
	   j0=izmax
	   do j=i+1,izmax
	     if (rrp(j).gt.rrp(i-1)) then
	       j0=j
	       goto 50
	     endif
	   enddo
 50        continue
c        -----------------------------------------------------
c        calculation of the new value rrp(i) using the equation:
c        (zzp(i)-zzp(i-1))/(rrp(i)-rrp(i-1))=
c        (zzp(j0)-zzp(i-1))/(rrp(i0)-rrp(i-1))
c        -----------------------------------------------------
            if(zzp(j0).eq.zzp(i-1)) then
	       rrp(i)=rrp(i)+1.d-5
	    else
	       rrp(i)=rrp(i-1)+(zzp(i)-zzp(i-1))*(rrp(j0)-rrp(i-1))/
     1	       	      (zzp(j0)-zzp(i-1))
	    endif
	 endif
      enddo
c--------near rmax
c     ------------------------------------
c     the correction for the ip-1 point (begin)
c      write(*,*)'before correction ip',ip
c      write(*,*)'before correction rrp(ip-1),rrp(ip)',rrp(ip-1),rrp(ip)
c      write(*,*)'before correction zzp(ip-1),zzp(ip)',zzp(ip-1),zzp(ip)
      if (rrp(ip-1).gt.rrp(ip)-1.d-5) then
         do j=ip-2,izmax,-1
            if (rrp(j).lt.rrp(ip-1)) then
	       j0=j
	       goto 55
	    endif
	 enddo
 55      continue
c         write(*,*)'j0,zzp(j0),zzp(ip)',j0,zzp(j0),zzp(ip)
         if(zzp(j0).eq.zzp(ip)) then
	    rrp(ip-1)=rrp(ip)-1.d-5
c	    write(*,*)'1 cor rrp(ip-1)',rrp(ip-1)
	 else
cSmirnov961226 beg
	    if(dabs(rrp(ip-1)-rrp(ip)).lt.1.d-6) then
	       rrp(ip-1)=0.5d0*(rrp(j0)+rrp(ip))
	    else
	       rrp(ip-1)=rrp(ip)+(zzp(ip-1)-zzp(ip))*(rrp(j0)-rrp(ip))/
     1	       	      (zzp(j0)-zzp(ip))
	    endif
cSmirnov961226 end
c	    write(*,*)'2 cor rrp(ip-1)',rrp(ip-1)
	 endif
      endif
c      write(*,*)'after correction rrp(ip-1),rrp(ip)',rrp(ip-1),rrp(ip)
c      write(*,*)'after correction zzp(ip-1),zzp(ip)',zzp(ip-1),zzp(ip)
c     the correction for the ip-1 point (end)
c     ------------------------------------

      do i=ip-2,izmax+1,-1
c         if(rrp(i).ge.rrp(i+1)) then
         if(rrp(i).ge.(rrp(i+1)-5.d-3)) then
	   j0=izmax
	   do j=i-1,izmax,-1
	     if (rrp(j).lt.rrp(i+1)) then
	       j0=j
	       goto 60
	     endif
	   enddo
 60        continue
c        -----------------------------------------------------
c        calculation of the new value rrp(i) using the equation:
c        (zzp(i)-zzp(i+1))/(rrp(i)-rrp(i+1))=
c        (zzp(j0)-zzp(i+1))/(rrp(j0)-rrp(i+1))
c        -----------------------------------------------------
            if(zzp(j0).eq.zzp(i+1)) then
	       rrp(i)=rrp(i)-1.d-5
	    else
	       rrp(i)=rrp(i+1)+(zzp(i)-zzp(i+1))*(rrp(j0)-rrp(i+1))/
     1	       	      (zzp(j0)-zzp(i+1))
	    endif
	 endif
      enddo

c-----------------------------------------------------------------
c  creation of monotonic limiter array zzm  near rmin and rmax

      izmin=1
      zmin=zzm(1)
      do i=2,im
        if (zmin.gt.zzm(i)) then
          zmin=zzm(i)
          izmin=i
	endif
      enddo
c      write(*,*)'im,izmin,zzm(izmin),rrm(izmin)',
c     1 im,izmin,zzm(izmin),rrm(izmin)

c--------near rmin
c     the correction for the second point (begin)
c      write(*,*)'before correction rrm(2),rrm(1)',rrm(2),rrm(1)
c      write(*,*)'before correction zzm(2),zzm(1)',zzm(2),zzm(1)
      if (rrm(2).lt.rrm(1)+1.d-5) then
         do j=3,izmax
            if (rrm(j).gt.rrm(2)) then
	       j0=j
	       goto 65
	    endif
	 enddo
 65      continue
         if(zzm(j0).eq.zzm(1)) then
	    rrm(2)=rrm(1)+1.d-5
	 else
cSmirnov961226 beg
	    if(dabs(rrm(1)-rrm(2)).lt.1.d-6) then
	       rrm(2)=0.5d0*(rrm(j0)+rrp(1))
	    else
	       rrm(2)=rrm(1)+(zzm(2)-zzm(1))*(rrm(j0)-rrm(1))/
     1	       	      (zzm(j0)-zzm(1))
	    endif
cSmirnov961226 end
	 endif
      endif
c      write(*,*)'after correction rrm(2),rrm(1)',rrm(2),rrm(1)
c      write(*,*)'after correction zzm(2),zzm(1)',zzm(2),zzm(1)
c     the correction for the second point (end)

      do i=3,izmin-1
c         if(rrm(i).le.rrm(i-1)) then
         if(rrm(i).le.(rrm(i-1)+5.d-3)) then
	   j0=izmin
	   do j=i+1,izmin
	     if (rrm(j).gt.rrm(i-1)) then
	       j0=j
	       goto 70
	     endif
	   enddo
 70        continue
c        -----------------------------------------------------
c        calculation of the new value rrp(i) using the equation:
c        (zzm(i)-zzm(i-1))/(rrm(i)-rrm(i-1))=
c        (zzm(j0)-zzm(i-1))/(rrm(i0)-rrm(i-1))
c        -----------------------------------------------------
            if(zzm(j0).eq.zzm(i-1)) then
	       rrm(i)=rrm(i)+1.d-5
	    else
	       rrm(i)=rrm(i-1)+(zzm(i)-zzm(i-1))*(rrm(j0)-rrm(i-1))/
     1	       	      (zzm(j0)-zzm(i-1))
	    endif
	 endif
      enddo
c--------near rmax
c     ------------------------------------
c     the correction for the im-1 point (begin)
c      write(*,*)'before correction im',im
c      write(*,*)'before correction rrm(im-1),rrm(im)',rrm(im-1),rrm(im)
c      write(*,*)'before correction zzm(im-1),zzm(im)',zzm(im-1),zzm(im)
      if (rrm(ip-1).gt.rrm(im)-1.d-5) then
         do j=im-2,izmax,-1
            if (rrm(j).lt.rrm(im-1)) then
	       j0=j
	       goto 75
	    endif
	 enddo
 75      continue
         if(zzm(j0).eq.zzm(im)) then
	    rrm(im-1)=rrm(im)-1.d-5
	 else
cSmirnov961226 beg
	    if(dabs(rrm(im-1)-rrm(im)).lt.1.d-6) then
	       rrm(im-1)=0.5d0*(rrm(j0)+rrm(im))
	    else
	       rrm(im-1)=rrm(im)+(zzm(im-1)-zzm(im))*(rrm(j0)-rrm(im))/
     1	       	      (zzm(j0)-zzm(im))
	    endif
cSmirnov961226 end
	 endif
      endif
c      write(*,*)'after correction rrm(im-1),rrm(im)',rrm(im-1),rrm(im)
c      write(*,*)'after correction zzm(im-1),zzm(im)',zzm(im-1),zzm(im)
c     the correction for the im-1 point (end)
c     ------------------------------------
      do i=im-2,izmin+1,-1
c         if(rrm(i).ge.rrm(i+1)) then
         if(rrm(i).ge.(rrm(i+1)-5.d-3)) then
	   j0=izmin
	   do j=i-1,izmin,-1
	     if (rrm(j).lt.rrm(i+1)) then
	       j0=j
	       goto 80
	     endif
	   enddo
 80        continue
c        -----------------------------------------------------
c        calculation of the new value rrm(i) using the equation:
c        (zzm(i)-zzm(i+1))/(rrm(i)-rrm(i+1))=
c        (zzm(j0)-zzm(i+1))/(rrm(j0)-rrm(i+1))
c        -----------------------------------------------------
            if(zzm(j0).eq.zzm(i+1)) then
	       rrm(i)=rrm(i)-1.d-5
	    else
	       rrm(i)=rrm(i+1)+(zzm(i)-zzm(i+1))*(rrm(j0)-rrm(i+1))/
     1	       	      (zzm(j0)-zzm(i+1))
	    endif
	 endif
      enddo

c-----------------------------------------------------------------
c  end of arrays rrp,zzp  and  rrm,zzm  creation
c------------------------------------------------------------------
      write(*,*)'in limitr after monotonization'
      do i=1,ip
        rpl(i)=rrp(i)
        zpl(i)=zzp(i)
c        write(*,*)'i,rrp(i),zzp(i)',i,rrp(i),zzp(i)
      enddo
      do i=1,im
        rml(i)=rrm(i)
        zml(i)=zzm(i)
c        write(*,*)'in i,rrm(i),zzm(i)',i,rrm(i),zzm(i)
      enddo
      return
      end







c***********************MAKELACN****************************************
c      it calculates the parameters
c      zmax1,rmaz1,zmin1,rmin1 for the eqdsk without
c      limiter data (nnlim=0)
c      zmax1,rmaz1,zmin1,rmin1 are closed but inside the
c      boundaries of the Lackner rectangle
c             ----------------------
c            |         . zmax1      |
c            |	                    |
c            |	   the points       |
c            |	     of the	        |
c     rmin1 .|      last            |.rmax1
c            |	     closed	        |
c            |	     flux           |
c            |	    surface         | Lackner Rectangle
c            |                      |
c             ----------------------
c                       . zmin1
c
c     (these points are inside the  separatrix line)
c----------------------------------------------------------------------
c      input data:
c        peqd(nreqda,nzeqda) is the psi function ,in common/fourb/
c        psisep,psilim     in common/three/
c        rma,zma           coordinates of magnetic axis in common/three/
c        nreqd,nzeqd       the number of points in the horizontal(r)
c                          and vertical(z) directions,  in common/three/
c        req(nreqda),zeq(nzeqda) eqdsk mesh points, in common/fourb/
c-----------------------------------------------------------------------
c      output data:
c         zmax1,rmax1,zmin1,rmin1 -boundaries of Lackner rectangle
c         rzmax1,zrmax1,rzmin1,zrmin1
c-----------------------------------------------------------------------
      subroutine makelacn(zmax1,rmax1,zmin1,rmin1,
     1 rzmax1,zrmax1,rzmin1,zrmin1)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'three.i'
      include 'fourb.i'
      dimension ix(nreqda)
      dimension rminj(nzeqda),rmaxj(nzeqda)
      dimension irminj(nzeqda),irmaxj(nzeqda)
      dimension psiminj(nzeqda),rpsiminj(nzeqda)
c      write(*,*)'in equilib in makelacn'
c------------------------------------------------------------
c     determination of eqdsk mesh point close to magnetic axis
c     zeq(imag) close to zma, req(jmag) close to rma
c------------------------------------------------------------
c      write(*,*)'in makelacn nreqd,nzeqd',nreqd,nzeqd
      dif=10.d15
      imag=1
      do i=1,nreqd
	difnew=dabs(req(i)-rma)
	if (difnew.lt.dif) then
	   dif=difnew
	   imag=i
	endif
      enddo !i

      dif=10.d15
      jmag=1
      do j=1,nzeqd
	difnew=dabs(zeq(j)-zma)
	if (difnew.lt.dif) then
	   dif=difnew
	   jmag=j
	endif
      enddo !i

      psiminj(jmag)=peqd(imag,jmag)
      rpsiminj(jmag)=req(imag)
      write(*,*)'imag,req(imag),rma',imag,req(imag),rma
      write(*,*)'jmag,zeq(jmag),zma',jmag,zeq(jmag),zma
      write(*,*)'peqd(imag,jmag),psimag',peqd(imag,jmag),psimag
c------------------------------------------------------------
c     Determination of rminj(nzeqd),rmaxj(nzeqd) -min and max major
c     radius(at given j) of the eqdsk points which are outside the plasma
c     and have the minimum distance from the plasma.
c     Determination of irminj(nzeqd),irmaxj(nzeqd) -the numbers in array r(i)
c     rminj(j)=r(irminj(j)), rmaxj(j)=r(irmax(j))
c     Determination of psiminj(nzeqd) and rpsimin(nzeqd)
c     psiminj(j) is the minimum value of psi along the major radius
c     from r=rminj(j) to r =rmaxj(j).
c     psiminj(j)=min{i=irminj(j)+1,irmaxj(j)-1} peqd(i,j)
c     rpsimin(j) is the value of major radius in which peqd(i,j)=psiminj(j)
c
c           {z(j),r(irminj(j)}             	{z(j),r(irmaxj(j)}
c     z(j) : .  .  .  *|  .  .   .   .   .  .  .  | *  .  .
c	        rminj(j)         plasma	            rmaxj(j)
c	       irminj(j)         	            irmaxj(j)
c------------------------------------------------------------
      jzmin=nzeqd
      jzmax=1
      do j=jmag,nzeqd
c     -------------------------------------------------------
c     for the eqdsk points with z.ge.zma
c     -------------------------------------------------------
	 ki=0
c        -------------------------------------------------------
c        ix(ki) is the array of the numbers i of eqdsk mesh req(i) for which
c        peqd(ix,j) .lt. psilim	, ki is the quantity of the such points
c        -------------------------------------------------------
         do i=1,nreqd
            if (peqd(i,j).le.psilim) then
               ki=ki+1
               ix(ki)=i
            endif
         enddo   !i

         if (ki.ge.1) then
c           ---------------------------------------------------
c           in this case there are exist the points in which
c           peqd(i,j).le.psilim for the given j
c           determination of the point of
c           min{i=irmin,irmax}peqd(i,j)  	for the given j
c           ---------------------------------------------------
	    irmin=ix(1)-1
	    irmax=ix(ki)+1
	    if (irmin.lt.1) irmin=1
	    if (irmax.gt.nreqd) irmax=nreqd
	    jzmax=j
	    psiminj(j)=peqd(irmin,j)
	    rpsiminj(j)=req(irmin)
	    do i=1,ki
	       ip=ix(i)
	       if(peqd(ip,j).lt.psiminj(j)) then
	          psiminj(j)=peqd(ip,j)
	          rpsiminj(j)=req(ip)
	       endif
	    enddo !i
c	    write(*,*)'in equilib irmin',irmin
	    rminj(j)=req(irmin)
	    irminj(j)=irmin
c	    write(*,*)'in equilib irmax',irmax
	    rmaxj(j)=req(irmax)
	    irmaxj(j)=irmax
c         write(*,*)'j,rminj(j),rmaxj(j),psiminj(j),rpsiminj(j)',
c     1              j,rminj(j),rmaxj(j),psiminj(j),rpsiminj(j)
	 else
c           ---------------------------------------------------
c           in this case there are no points in which
c           peqd(i,j).le.psilim for the given j
c           ---------------------------------------------------
	    rminj(j)=req(nreqd)
	    irminj(j)=nreqd
	    rmaxj(j)=req(1)
	    irmaxj(j)=1

	    jzmax=j
	    rzmax1=rpsiminj(jzmax)
            goto 10
         endif

      enddo !j
 10   continue
      do j=jmag-1,1,-1
c     -------------------------------------------------------
c     for the eqdsk points with z.lt.zma
c     -------------------------------------------------------
         ki=0
c        -------------------------------------------------------
c        ix(ki) is the array of the numbers i of eqdsk mesh req(i) for which
c        peqd(ix,j) .lt. psilim	, ki is the quantity of the such points
c        -------------------------------------------------------
         do i=1,nreqd
	    if (peqd(i,j).le.psilim) then
	       ki=ki+1
	       ix(ki)=i
	    endif
	 enddo   !i

	 if (ki.ge.1) then
c           ---------------------------------------------------
c           in this case there are exist points where
c           peqd(i,j).le.psilim for the given j	,
c           determination of the point of
c           min{i=irmin,irmax}peqd(i,j)  	for the given j
c           ---------------------------------------------------
	    irmin=ix(1)-1
	    irmax=ix(ki)+1
	    if (irmin.lt.1) irmin=1
	    if (irmax.gt.nreqd) irmax=nreqd
	    jzmin=j
	    psiminj(j)=peqd(irmin,j)
	    rpsiminj(j)=req(irmin)
	    do i=1,ki
	       ip=ix(i)
	       if(peqd(ip,j).lt.psiminj(j)) then
	          psiminj(j)=peqd(ip,j)
	          rpsiminj(j)=req(ip)
	       endif
	    enddo !i

	    rminj(j)=req(irmin)
	    irminj(j)=irmin
	    rmaxj(j)=req(irmax)
	    irmaxj(j)=irmax
c         write(*,*)'j,rminj(j),rmaxj(j),psiminj(j),rpsiminj(j)',
c     1              j,rminj(j),rmaxj(j),psiminj(j),rpsiminj(j)
	 else
c           ---------------------------------------------------
c           in this case there are no points in which
c           peqd(i,j).le.psilim for the given j
c           ---------------------------------------------------
cc	      write(*,*)'out of plasmaj,zeq(j)',j,zeq(j)
	    rminj(j)=req(nreqd)
	    irminj(j)=nreqd
	    rmaxj(j)=req(1)
	    irmaxj(j)=1

	    jzmin=j
	    rzmin1=rpsiminj(jzmin)
	    goto 20
         endif

      enddo !j
 20   continue
c-----------------------------------------------------------------------
c     determination of zmax1 and zmin1
c-----------------------------------------------------------------------
cc      write(*,*)'jzmin,jzmax,nzeqd',jzmin,jzmax,nzeqd
      if(jzmax.lt.nzeqd) then
         zmax1=zeq(jzmax)
	 rzmax1=rpsiminj(jzmax-1)
      else
c        ---------------------------------------------------
c        it is necessary to calculate the jzmax1
c        (in this case the separatrix area exist close to plasma)
c        ---------------------------------------------------
         do j=jmag+1,nzeqd
cc	    write(*,*)'j,psiminj(j),psiminj(j-1)'
cc	    write(*,*)j,psiminj(j),psiminj(j-1)
            if(psiminj(j).gt.psiminj(j-1))then
	      zmax1=zeq(j)
	      jzmax=j
cc	      write(*,*)'j,zmax1,jzmax',j,zmax1,jzmax
	    else
	       goto 30
            endif
         enddo
 30      continue
cc	 jzmax=jzmax-1
	 jzmax=jzmax
	 zmax1=zeq(jzmax)
	 rzmax1=rpsiminj(jzmax)
cc         write(*,*)'jzmax,zmax1',jzmax,zmax1
      endif

      if(jzmin.gt.1) then
         zmin1=zeq(jzmin)
         rzmin1=rpsiminj(jzmin+1)
      else
c        ---------------------------------------------------
c        it is necessary to calculate the jzmin1
c        (in this case the separatrix area exist close to plasma)
c        ---------------------------------------------------
         do j=jmag-1,1,-1
cc	    write(*,*)'j,psiminj(j),psiminj(j+1)'
cc	    write(*,*)j,psiminj(j),psiminj(j+1)
            if(psiminj(j).gt.psiminj(j+1))then
	       zmin1=zeq(j)
	       jzmin=j
cc	       write(*,*)'j,zmin1,jzmin',j,zmin1,jzmin
	    else
	       goto 40
            endif
         enddo
 40      continue
cc	 jzmin=jzmin+1
	 jzmin=jzmin
	 zmin1=zeq(jzmin)
	 rzmin1=rpsiminj(jzmin)
      endif
cc      write(*,*)'jzmin,zmin1',jzmin,zmin1
c---------------------------------------------------------------
      rmax1=req(1)
      rmin1=req(nreqd)
c-------------------------------------------
c      determination of rmax1 and rmin1
c-------------------------------------------
      do j=jzmin,jzmax
	 if (rminj(j).lt.rmin1) then
	    rmin1=rminj(j)
	    zrmin1=zeq(j)
	 endif
	 if (rmaxj(j).gt.rmax1) then
	    rmax1=rmaxj(j)
	    zrmax1=zeq(j)
	 endif
      enddo !j

c      write(*,*)'rmax1,rmin1,zmax1,zmin1',rmax1,rmin1,zmax1,zmin1
      return
      end






c     ****************PRESPSI********************************
c     * this function calculates plasma pressure on psi	     
c     *******************************************************
c     ! input parameter: psi				    !
      double precision
     1function prespsi(psi)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'five.i'
      double precision ias1r
c     idx derivativs order 0.ge.idx.le.3
      idx=0
cSm030224
      nr4=nx+4
      prespsi=ias1r(tpres,nx,nr4,cpres,idx,psi)
      return
      end




           
      subroutine wrtnetcdf_wall_limiter_data(netcdfnml)
 
c-------------------------------------------------------------
c     writes wall and limiter coordinates to the existing
c     netcdf file   netcdfnml
c
c     The wall and limiters data are in /one.nml/ and /fourb/ 

c     n_wall is a number of wall points
c
c     max_limiters is a number of limiters
c
c     wall and limiter (r,z) coordinates in [m]
c     r_wall(1:n_wall)   
c     z_wall(1:n_wall)  
c
c     n_limiter(1:max_limiters) number of limiter points
c
c     r_limiter(1:n_limiter_a,1:max_limiters_a)   
c     z_limiter(1:n_limiter_a,1:max_limiters_a)
c--------------------------------------------------------------
      implicit none
      include 'param.i'
      include 'one.i'
      include 'fourb.i'

      include 'netcdf.inc'
       

c-----input
      character*(*) netcdfnml      !name of output nc file
    
c-----locals----------------------------------------  
      integer i,j,n_limiter_max

c-----Storage tem1 is used in netcdf writes
      integer n1n2
      parameter (n1n2=n_limiter_a*max_limiters_a)
      real*8 tem1(n1n2)
      real*8 tem12(2*max_limiters_a)

c --- some stuff for netCDF file ---
      integer ncid,vid,istatus,n_wall_id,max_limiters_id,ncvdef2,ncdid2
     + ,ncddef2
      integer n_limiter_max_id,i_two_id,two

      integer limiter_dims(2),start(2),limiter_count(2),   
     &phi_limiter_dims(2),phi_limiter_count(2)
    
      data start/1,1/

c-----calculate maximal value of limiter points n_limiter_max
c     at all limiters (from 1 to max_limiters)
      n_limiter_max=0 
      n_limiter_max_id=0
      max_limiters_id=0

      write(*,*)'in sub wrtnetcdf_wall_limiter_data'
 
      do i=1,max_limiters
          n_limiter_max=max(n_limiter_max,n_limiter(i))
         write(*,*)'wrtnetcdf_wall_limiter_data'
         write(*,*)'i,n_limiter(i),n_limiter_max',
     &              i,n_limiter(i),n_limiter_max
      enddo
      write(*,*)'n_limiter_max',n_limiter_max


      limiter_count(1)=n_limiter_max
      limiter_count(2)=max_limiters
      write(*,*)'limiter_count',limiter_count
      two=2
      phi_limiter_count(1)=two
      phi_limiter_count(2)=max_limiters
c.......................................................................
cl    1. Initialize part
c
c-----------------------------------------------------------------------
c
c     1.1 open existing netCDF file: netcdfnml, 
c         and define additional;
c         dimensions,variables and attributes
c-------------------------------------------------------------
c      create net CDF file  define dimensions, variables
c      and attributes
c------------------------------------------------------------

c.......................................................................
c     1.1.1 open existing netCDF file: netcdfnml,
c     create netCDF filename (Entering define mode.)
c     integer function ncopn(filename,write,error_code)
c     Ref to page 40 of NetCDF-2 manual.
c     CLOBber old file, if it exists.
c     istatus is 0, if no errors.

      write(*,*)'in netcdf_wall_limiter_data netcdfnml=',netcdfnml

      call ncopn2(netcdfnml,NCWRITE,ncid,istatus) ! Open existing netCDF file

      write(*,*)'in netcdf_wall_limiter_data after ncopn(n netcdfnml'

      call check_err(istatus)
c,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,      
c     1.1.2 
c     Put Open NetCDF file into Define Mode
c     subroutine ncredf(integer ncid,integer error_code)
c     p. 41 of netcdf manual
c

      call ncredf3(ncid,istatus)

      call check_err(istatus)

c-------------------------------------------------------------------
c     1.1.3 determine dimension ID of 
c
c.......................................................................
cl    1.1.2 define dimensions
c.......................................................................
c      write(*,*)'before ncddef(ncid,n_wall nwall=',n_wall

      if (n_wall.gt.0)
     &n_wall_id=ncddef2(ncid,'n_wall',n_wall,istatus)
      call check_err(istatus)

      if (n_limiter_max.gt.0)
     &n_limiter_max_id=ncddef2(ncid,'n_limiter_max',n_limiter_max,
     &                        istatus)

      write(*,*)'netcdf_wall_limiter_data'
      write(*,*)'n_wall',n_wall
      write(*,*)'max_limiters',max_limiters

      call check_err(istatus)

      i_two_id=ncddef2(ncid,'i_two_id',two,istatus)
      call check_err(istatus)

      limiter_dims(1)=n_limiter_max_id
      limiter_dims(2)=max_limiters_id
      write(*,*)'limiter_dims',limiter_dims
      phi_limiter_dims(1)=i_two_id
      phi_limiter_dims(2)=max_limiters_id
c.......................................................................
cl    1.1.3 define variables
c.......................................................................
      vid=ncvdef2(ncid,'n_wall',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,21,
     +            'Number of wall points',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'max_limiters',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,18,
     +            'Number of limiters',istatus)
      call check_err(istatus)

     
      if (n_wall.gt.0) then

         vid=ncvdef2(ncid,'r_wall',NCDOUBLE,1,n_wall_id,istatus)
         call check_err(istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,17,
     &               'wall r coordinate',istatus)
         call check_err(istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,1,
     &              'm',istatus)
         call check_err(istatus)

         vid=ncvdef2(ncid,'z_wall',NCDOUBLE,1,n_wall_id,istatus)
         call check_err(istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,17,
     &               'wall z coordinate',istatus)
         call check_err(istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,1,
     &              'm',istatus)
         call check_err(istatus)
      endif

      if(max_limiters.gt.0)then

         vid=ncvdef2(ncid,'n_limiter',NCLONG,1,max_limiters_id,istatus)
         call check_err(istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +            'Numbers of limiter points',istatus)
         call check_err(istatus)

        vid=ncvdef2(ncid,'r_limiter',NCDOUBLE,2,limiter_dims,istatus)
        call check_err(istatus)
        call ncaptc2(ncid,vid,'long_name',NCCHAR,20,
     +           'limiter r coordinate',istatus)
        call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'm',istatus)
        call check_err(istatus)

        vid=ncvdef2(ncid,'z_limiter',NCDOUBLE,2,limiter_dims,istatus)
        call check_err(istatus)
        call ncaptc2(ncid,vid,'long_name',NCCHAR,20,
     +           'limiter z coordinate',istatus)
        call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'm',istatus)
        call check_err(istatus)

       write(*,*)'before vid=ncvdef2(phi_limiter'
        vid=ncvdef2(ncid,'phi_limiter',NCDOUBLE,2,phi_limiter_dims,
     &            istatus)
        call check_err(istatus)
        call ncaptc2(ncid,vid,'long_name',NCCHAR,31,
     +           'toroidal bounadries of limiters',istatus)
        call ncaptc2(ncid,vid,'units',NCCHAR,6,
     +           'degree',istatus)
        call check_err(istatus)
        write(*,*)'after vid=ncvdef2(phi_limiter'
      endif

cl    1.1.5 end the define-mode and start the data-mode
c     p. 51-2 of manual

      call ncendf3(ncid,istatus)
      call check_err(istatus)
          
c----------------------------------------------------------------
c     write data
c-----------------------------------------------------------------
      call ncvid2(vid,ncid,'n_wall',istatus)
      call ncvpt_int2(ncid,vid,0,0,n_wall,istatus)

      call ncvid2(vid,ncid,'max_limiters',istatus)
      call ncvpt_int2(ncid,vid,0,0,max_limiters,istatus)

      write(*,*)'n_wall',n_wall

      if (n_wall.gt.0) then
         do i=1,n_wall
            write(*,*)'i,r_wall(i)',i,r_wall(i)
         enddo

         write(*,*)'ncid',ncid 
         write(*,*)'before vid=ncvid(ncid,r_wall'
         call ncvid2(vid,ncid,'r_wall',istatus)
         call ncvpt_doubl2(ncid,vid,1,n_wall,r_wall,istatus)
         call check_err(istatus)

          write(*,*)'before vid=ncvid(ncid,z_wall'
         call ncvid2(vid,ncid,'z_wall',istatus)
         call ncvpt_doubl2(ncid,vid,1,n_wall,z_wall,istatus) 
         call check_err(istatus)
      endif

c---------------------------------------------------
c         4. Close netCDF file
c        Page 53-54 netCDF-2 Manual
c
      call ncclos3(ncid,istatus)
      call check_err(istatus)

      return
      end       

      

     
      subroutine netcdf_eqdsk_data(netcdfnml)
c-------------------------------------------------------------
c     writes eqdsk poloidal flux array 
c     (peqd(i,j),i=1,nreqd),j=1,nzeqd)
c     and xeq(nxeqda),yeq(nyeqda),req(nreqda),zeq(nzeqda) mesh 
c     to the existing netcdf file   netcdfnml
c
c--------------------------------------------------------------
      implicit none
      include 'param.i'
      include 'fourb.i' ! contains peqd()
      include 'netcdf.inc'
      include 'three.i' ! contains psimag,psilim
c-----input
      character*(*) netcdfnml      !name of output nc file
    
c-----locals----------------------------------------  
      integer i
c --- some stuff for netCDF file ---
      integer ncid,vid,istatus,nxeqd_id,nyeqd_id,nreqd_id,nzeqd_id,
     + ncvdef2,ncdid2,ncddef2
      integer psi_dims(2),start(2),psi_count(2)

c-----Storage tem1 is used in netcdf writes.
      integer nxny
      parameter (nxny=nreqda*nzeqda)
      real*8 tem1(nxny)
 
      data start/1,1/
c.......................................................................
cl    1. Initialize part
c
c-----------------------------------------------------------------------
c
c     1.1 open existing netCDF file: netcdfnml, 
c         and define additional;
c         dimensions,variables and attributes
c-------------------------------------------------------------
c      create net CDF file  define dimensions, variables
c      and attributes
c------------------------------------------------------------

c.......................................................................
c     1.1.1 open existing netCDF file: netcdfnml,
c     create netCDF filename (Entering define mode.)
c     integer function ncopn(filename,write,error_code)
c     Ref to page 40 of NetCDF-2 manual.
c     CLOBber old file, if it exists.
c     istatus is 0, if no errors.

      write(*,*)'in  netcdf_eqdsk_data=',netcdfnml

      call ncopn2(netcdfnml,NCWRITE,ncid,istatus) ! Open existing netCDF file
      call check_err(istatus)
c,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,      
c     1.1.2 
c     Put Open NetCDF file into Define Mode
c     subroutine ncredf(integer ncid,integer error_code)
c     p. 41 of netcdf manual
c
      call ncredf3(ncid,istatus)
      call check_err(istatus)
c-------------------------------------------------------------------
c     1.1.3 determine dimension ID of 
c
c.......................................................................
cl    1.1.2 define dimensions
c.......................................................................
      nxeqd_id=ncddef2(ncid,'nxeqd',nxeqd,istatus) ! YuP added xyz-grid
      nyeqd_id=ncddef2(ncid,'nyeqd',nyeqd,istatus) ! YuP added xyz-grid
      nreqd_id=ncddef2(ncid,'nreqd',nreqd,istatus)
      nzeqd_id=ncddef2(ncid,'nzeqd',nzeqd,istatus)
      
      write(*,*)'nxeqd',nxeqd,nxeqd_id
      write(*,*)'nyeqd',nyeqd,nyeqd_id
      write(*,*)'nreqd',nreqd,nreqd_id
      write(*,*)'nzeqd',nzeqd,nzeqd_id

c.......................................................................
cl    1.1.3 define variables
c.......................................................................
      vid=ncvdef2(ncid,'nxeqd',NCLONG,0,0,istatus) ! YuP added xyz-grid
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,18,
     +            'Number of x points',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'nyeqd',NCLONG,0,0,istatus) ! YuP added xyz-grid
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,18,
     +            'Number of y points',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'nreqd',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,18,
     +            'Number of r points',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'nzeqd',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,18,
     +            'Number of z points',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'psimag',NCDOUBLE,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,34,
     +            'poloidal flux/2pi on magnetic axis',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,2,
     &           'Wb',istatus)

      vid=ncvdef2(ncid,'psilim',NCDOUBLE,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,26,
     +            'poloidal flux/2pi at rho=1',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,2,
     &           'Wb',istatus)


      vid=ncvdef2(ncid,'eqdsk_x',NCDOUBLE,1,nxeqd_id,istatus) !added xyz-grid
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,19,
     &            'equilib. B  x-array',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,1,
     &           'm',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'eqdsk_y',NCDOUBLE,1,nyeqd_id,istatus) !added xyz-grid
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,19,
     &            'equilib. B  y-array',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,1,
     &           'm',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'eqdsk_r',NCDOUBLE,1,nreqd_id,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,19,
     &            'equilib. B  r-array',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,1,
     &           'm',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'eqdsk_z',NCDOUBLE,1,nzeqd_id,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,19,
     &            'equilib. B  z-array',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,1,
     &           'm',istatus)
      call check_err(istatus)

      psi_dims(1)=nreqd_id
      psi_dims(2)=nzeqd_id

      psi_count(1)=nreqd
      psi_count(2)=nzeqd

      write(*,*)'psi_dims',psi_dims
      write(*,*)'psi_count',psi_count

      vid=ncvdef2(ncid,'eqdsk_psi',NCDOUBLE,2,psi_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,19,
     +         'poloidal flux /2pi ',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,2,
     &           'Wb',istatus)
      call check_err(istatus)
      
      vid=ncvdef2(ncid,'eqdsk_beq',NCDOUBLE,2,psi_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,27,
     +         'Total equil. Magnetic field',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     &           'Tesla',istatus)
      call check_err(istatus)

cl    1.1.5 end the define-mode and start the data-mode
c     p. 51-2 of manual

      call ncendf3(ncid,istatus)
      call check_err(istatus)

c----------------------------------------------------------------
c     write data
c-----------------------------------------------------------------
      write(*,*)'ncid=',ncid
      
      call ncvid2(vid,ncid,'nreqd',istatus)
      call ncvpt_int2(ncid,vid,0,0,nreqd,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'nxeqd',istatus) ! YuP added xyz-grid 
      call ncvpt_int2(ncid,vid,0,0,nxeqd,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'nyeqd',istatus) ! YuP added xyz-grid
      call ncvpt_int2(ncid,vid,0,0,nyeqd,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'nzeqd',istatus)
      call ncvpt_int2(ncid,vid,0,0,nzeqd,istatus)
      call check_err(istatus)
      
      !do i=1,nreqd
      !write(*,*)'equilib/netcdf: req=',i,req(i)
      !pause
      !enddo
      call ncvid2(vid,ncid,'psimag',istatus)
      call ncvpt_doubl2(ncid,vid,0,0,psimag,istatus)
      call check_err(istatus)
      call ncvid2(vid,ncid,'psilim',istatus)
      call ncvpt_doubl2(ncid,vid,0,0,psilim,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'eqdsk_r',istatus)
      call ncvpt_doubl2(ncid,vid,1,nreqd,req,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'eqdsk_x',istatus) ! YuP added xyz-grid
      call ncvpt_doubl2(ncid,vid,1,nxeqd,xeq,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'eqdsk_y',istatus) ! YuP added xyz-grid
      call ncvpt_doubl2(ncid,vid,1,nyeqd,yeq,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'eqdsk_z',istatus)
      call ncvpt_doubl2(ncid,vid,1,nzeqd,zeq,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'eqdsk_psi',istatus)
      call check_err(istatus)
        tem1=0.d0
        call pack21(peqd,1,nreqda,1,nzeqda,tem1,nreqd,nzeqd)
        call ncvpt_doubl2(ncid,vid,start,psi_count,tem1,istatus)
      !call ncvpt_doubl2(ncid,vid,start,psi_count,peqd(1,1),istatus)
      call check_err(istatus)
      write(*,*)'MIN/MAX of peqd(i,j) in saving to genray.nc:',
     +  minval(peqd), maxval(peqd) 

c---------------------------------------------------
c         4. Close netCDF file
c        Page 53-54 netCDF-2 Manual
c
      call ncclos3(ncid,istatus)
      call check_err(istatus)
      

      return
      end       

      

     
