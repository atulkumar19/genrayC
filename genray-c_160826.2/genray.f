
 
c        ********************* GENRAY-c *********************
c        *                     ------                       *
c        * GENRAY-c is the computer code for obtaining  the *
c        * trajectories of different  wave  modes  by  ray- *
c        * tracing techniques; using cartesian coords.      *
c        *                                                  *
c        * Authors: Alexander P. Smirnov (primary)          *
c        *          Moscow State University                 *
c        *          (sap@ns.cnt.ru, sap@cs.msu.su)          *
c        *          R.W. Harvey          (secondary)        *
c        *          CompX                                   *
c        *          (bobh@compxco.com)                      *
c        *          Yuri Petrov                             *
c        *          CompX                                   *
c        *          (petrov@compxco.com)                    *
c        *                                                  *
c        * Manual: GENRAY, Report CompX-01-2000 (2000)      *
c        *                 CompX, PO Box 2672, Del Mar, CA  * 
c        *                                                  *
c        ****************************************************
c
c-----------------------------------------------------------------!
c								  !
c           method of solution and coordinate system.		  !
c								  !
c        the wave trajectories are obtained from the solution     !
c        of  geometrical  optics  equations  where   independent  !
c        space variables are: 					  
c								  
c           x,y,z - Cartesian coordinates			  
c								  
c        and canonically conjugate momenta are:		          
c								  
c           k_x, k_y, k_z.					  
c								  
c        Geometrical optics equation in these  variables  are     !
c        hamiltonian in form and they are  solved  by  4th-order  !
c        runge-kutta method.					  
c        The code maintains conservation of the hamiltonian       ! 
c        function with given accuracy eps by using two different  !
c        numerical methods :                                      !
c        for isolv=1 - the solution of six ray equations          !
c         with correction of coordinates at each time step.       ! 
c        for isolv=2                                              !
c         the solution of the dispersion equation for one ray	  !
c         variable and the determination the other 5 variables	  !
c         from 5 geometric optics equations. The code determines  !
c         automatically which ray variable which must be          !
c         determined from dispertion relation.			  !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
c                                                                 !
c         the code is written in fortran-77 and does not require  !
c       any  special  subroutines  and  functions  from  fortran  !
c       libraries.						  !
c                                                                 !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
c                                                                 !
c                   version   19/09/00                            !
c                                                                 !
c-----------------------------------------------------------------!
c        it reads the input files: equilib.dat,genray.in  	  !
c-----------------------------------------------------------------!


CMPIINSERTPOSITION PROGRAMSTART


      PROGRAM GENRAY

c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none

      include 'param.i' ! specifies code parameters 
      include 'commons.i'

CMPIINSERTPOSITION DECLARATION

      real*8  u(6),deru(6),aux(8,6),
     &energy,pitch,fdist,dfdx,dfdpitch,dfdp, !for call dskin
     &r,z,phi,fdens,fdens0,fmaxw,
     &xst1,yst1,zst1,rst1,phist1,alfast1,betast1,
     &cnteta,cnphi,      
     &win_0dsmax,ds,r0_em,z0_em,
     &psi_loc,rho_loc,clight,sum_emission,sum_emission_wall


      integer iraystop,nray,ndim,
     &initial,ihlf,n_r,n_z,ifreq,
     &i_bad_initial_conditions,n,
     &nrayelt_o_cutoff_emis,n0,igenray,
     &i_geom_optic_loc

c-----externals
      external outpt, outpt_xyz
      external rside1, rside_xyz, dense_xyz,
     &fdens_fdist,fdens0_fdist,fdens_fkin2,
     &length_char,tempe,temperho,rhopsi,r_2nd_harm,
     & ddnsrho,densrho

      real*8 thetapol,psif,drhodz,dthetadz,dthetadr,drhodr,ddnsrho,b1,
     + den
      real*8 b,bxyz,dense,density_r_z_i,d_density_r_z_i_d_r,x,dxdr,dxdz,
     &d_density_r_z_i_d_z,drhopsi,dense_xyz,
     &fdens_fdist,fdens0_fdist,fdens_fkin2,tempe,temperho,rhopsi,
     &r_2nd_harm,zeffrho,
     &dense_no_RZ_spline
      real*8 xstart, ystart, zstart, rstart, phistart, 
     + cnzstart,cnrstart,cmstart, cnxstart,cnystart

      integer length_char

      real*4 time_loop_ifreq_1,time_loop_ifreq_2,time_loop_ifreq,
     &time_before_rk,time_after_rk,time_rk,
     &time_genray_1,time_genray_2,
     &time_emission_2,
     &time_emission_1
c------------------------------------------------------------
cfor test CD_adj_LH_efficiency
      real*8 lh_cd_efficiency,cnpar,cnper,thetapol_l,
     &u_ph_karney,u1,efficien,temp_kev,unorm

      integer n_radial
cfor interpolation_chi
      real*8 u0,theta0,
     &chi,d_chi_d_u_out,d_chi_d_theta_out
      integer n_radial0
cSAP080711
c-----genray input file: genray.in or genray.dat
      character*10 genray_in_dat

cSAP081122
      real*8 
     &length_b

cSAP090205 to check subroutine sigma_edge_n_theta_pol
      integer i,i0r,j0r,i0z,j0z

      real*8  theta_pol_radian,
     &sigma_edge_n,d_sigma_edge_n_d_theta_pol,
     &d2_sigma_edge_n_d2_theta_pol

ctest for spline density art rz mesh
      real*8 dens_rho_l,dens_rz,d_norm,diff_dens,
     &dif_d_dens_dr,dif_d_dens_dz,
     &d_norm_r,d_norm_z,d_dens_rho_r,d_dens_rho_z,
     &d_dens_spl_r,d_dens_spl_z,
     & dens_rho_theta,d_dens_rho_theta_d_rho,
     & d_dens_rho_theta_d_theta,
     & dn_drho,dro_dpsi,
     & step_rz,d_p,d_m,deriv_l,r_p,r_m,z_p,z_m,
     & x_l,x_p,x_m,dxdr_l,dxdz_l,dxdr_d,dxdz_d,
     & rho_z_p,rho_z_m,rho_r_p,rho_r_m,     
     & dens_z_p,dens_z_m,dens_r_p,dens_r_m
      integer j,k,m
      
CMPIINSERTPOSITION INITIALIZATION

        call cpu_time(time_genray_1)

c---------------------------------------------------
c     Write out the version number
c--------------------------------------------------

      write(t_,1000) version
 1000 format(15x,"GENRAY VERSION: ",a)
      write(6,'(//16x,
     +          "=================================================")')
      if (length_char(t_) .gt. 512) stop 'GENRAY:Adjust length of t_'
      write(6,*) t_(1:length_char(t_))
      write(6,'(16x,
     +          "=================================================",/)')


CMPIINSERTPOSITION STARTBARRIER

c---------------------------------------------------
c     check the parameters in param.i
c--------------------------------------------------
      if (nrya.ne.(max(nreqda,nzeqda)+4)) then
         write(*,*)'it should be nrya.eq.(max(nreqda,nzeqda)+4)'
         write(*,*)'but nrya,nreqda,nzeqda', nrya,nreqda,nzeqd
         write(*,*)'change the parameter nrya in param.i'
         stop
      endif

      if (nzy.ne.(max(npsi,nteta1)+4)) then
         write(*,*)'it should be nzy.eq.(max(npsi,nteta1)+4)'
         write(*,*)'but nzy,npsi,nteta1', nzy,npsi,nteta1
         write(*,*)'change the parameter nzy in param.i'
         stop
      endif
            
c-----the creation of default input data like in genray.dat file
      call default_in
      write(*,*)'genray.f after default_in'
      write(*,*)'i_resonance_curve_integration_method=',
     &           i_resonance_curve_integration_method

c-----check if input file is genray.in, then it will change
c     default data to genray.in file (MKSA) format  
      call transform_input_data_to_MKSA

c---- read all namelists from input genray.in or genray.dat file
      call read_all_namelists(genray_in_dat,ndim,nray)
      write(*,*)'in genray.f genray_in_dat=', genray_in_dat
      write(*,*)'genray.f after read_all_ prmt ',prmt

      write(*,*)'genray.f after read_all_ ioxm',ioxm   
      write(*,*)'genray.f after read_all_ nray',nray  

      write(*,*)'genray.f after read_all_namelists data in /genr/'
      write(*,*)'r0x=',r0x
      write(*,*)'b0=',b0
      write(*,*)'outdat=',outdat
      write(*,*)'stat=',stat
      write(*,*)'mnemonic=',mnemonic
      write(*,*)'rayop=',rayop
      write(*,*)'dielectric_op=',dielectric_op
      write(*,*)' =======>   ixyz=',ixyz
      
c-----allocate data 
      call ainalloc

c-----test sigma_edge_n
      pi=4*datan(1.d0)
      !write(*,*)' i_edge_dens_anal',i_edge_dens_anal
      if(i_edge_dens_anal.eq.2) then
c--------spline data for  sigmedgn_ar
         do i=1,n_pol_edge_dens
            theta_pol_radian=2.d0*pi/( n_pol_edge_dens-1)*(i-1)
            call sigma_edge_n_theta_pol(theta_pol_radian,
     &         sigma_edge_n,d_sigma_edge_n_d_theta_pol,
     &         d2_sigma_edge_n_d2_theta_pol)
            write(*,*)'i,theta_pol_radian,sigma_edge_n',
     &                 i,theta_pol_radian,sigma_edge_n
            write(*,*)'d_sigma_edge_n_d_theta_pol',
     &                 d_sigma_edge_n_d_theta_pol
            write(*,*)'d2_sigma_edge_n_d2_theta_pol',
     &                 d2_sigma_edge_n_d2_theta_pol
         enddo
       endif

c---------------------------------------------------------------
c     read the density profile data from file
c     and setup coefficients for 2D spline approximation
      if(model_rho_dens.eq.3)then
         call density_profile_read
         write(*,*)' genray.f after call density_profile_read'
      endif
c---------------------------------------------------------------

c-----If input file is genray.in, then it will change
c     input data to genray.dat file format 
      if (genray_in_dat.eq.'genray.in')  then
         write(*,*)'genray.f before transform_genray_in_to_dat'
         call transform_genray_in_to_dat 
      endif
      write(*,*)'genray.f after transform_genray  prmt=',prmt
c---------------------------------------------------------------
c     reading the eqdsk data
c     and creation of the coefficients for spline approximation
c     for: feqd, psi -magnetic field functions
c          zpllim(r),zminlim (r) - limmiter boundary
      !write(*,*)' genray.f before call equilib'
      call equilib
      !write(*,*)' genray.f after call equilib'
      !write(*,*)'NR',NR
c---------------------------------------------------------------
      if(model_b.eq.0) then
c         write(*,*) 'genray.f before call rhospl'
         call rhospl     
c         write(*,*) 'genray.f after call rhospl'
      else
         write(*,*) 'genray.f:  rhospl is not called for model_b>0'
      endif
c      write(*,*)'arpsi',arpsi
c----------------------------------------------------
c--------------------------------------------------------------------
c     reading the name of the output file
c--------------------------------------------------------------------
      i1_=91
      open(i1_,file=outdat)
      iraystop=0
c-------------------------------------------------------------------
c     reading input data from genray.in file
      write(*,*) 'genray.f before call dinit_mr'
      write(*,*) 'NR=',NR
      
      if(ixyz.eq.0) then
         call dinit_mr(ndim,nray)
      else
         call dinit_mr_xyz(ndim,nray) ! ixyz=1 Cartesian
         ! will give arxu0(iray), aryu0(iray), arzu0(iray),
         ! arnpar(iray), arnper_tang(iray)
      endif

c-----allocate pointers at writencdf.i
c      write(*,*)'genray.f before  ainalloc_writencdf nray',nray
c      call ainalloc_writencdf_i(nray)
c      call ainalloc_write_i(nray)

      write(*,*)'genray.f after dinit_mr prmt=',prmt
      write(*,*)'genray !!!!!ndim',ndim
      write(*,*)'genray after dinit_mr nray=',nray
      write(*,*)'genray after dinit_mr nbulk',nbulk
      write(*,*)'v',v
      write(*,*)'w',w
      write(*,*)'genray.f after dinit_mr ioxm ',ioxm
      write(*,*)'genray.f after dinit_mr freqncy0',freqncy0
      do iray=1,nray         
         write(*,*)'1 iray,arzu0(iray),arru0(iray),arphiu0(iray)',
     &              iray,arzu0(iray),arru0(iray),arphiu0(iray)
      enddo
      write(*,*)'genray.f before n_wall.eq.1  n_wall',n_wall 

      if (n_wall.gt.1)then
c------------------------------------------------------------------
c        create additional points at the chamber wall
c-----------------------------------------------------------------
c        interpolate wall points at the wall mesh with additinal ponts 
c        r_wall_add(n_wall_add) z_wall_add(n_wall_add)
c        and calculate the number of points of this mesh: n_wall_add)
c        Result will be in fourbe.i
         call create_fine_mesh_for_chamber_wall_limiter_coordinates
cyup         call add_horisontal_limiter_walls
         call MK_graph_chamber_wall
c         stop 'genray,f after  MK_graph_chamber_wall'
c--------------------------------------------------------------------
c        create RZ meshes: rr_add(nreqd_add),zz_add(nzeqd_add)
c-------------------------------------------------------------------
         if(ixyz.eq.0) call creat_fine_2D_poloidal_mesh
c--------------------------------------------------------------------
c        calculate distance from RZ mesh points to chamber wall:
c        distance_to_wall(nreqd_add,nzeqd_add)
c-------------------------------------------------------------------
         write(*,*)'genray.f before distance_for_wall'
         if(ixyz.eq.0) call distance_for_wall  
         write(*,*)'genray.f after distance_for_wall'
c--------test fast subroutine: it gave wrong result 
c         call distance_for_wall_1
c         stop 'genray.f after distance_for_wall_1'
c-------------------------------------------------------------------
c        Calculate density array:  density_r_z(nreqd_add,nzeqd_add,nbulk)
c        at (rr_add,zz_add) mesh
c        For( n_wall.gt.1) case create density fall near the chamber wall.
c-------------------------------------------------------------------
          if(ixyz.eq.0) call density_at_zr_plane 
c-------------------------------------------------------------------
c         calculate density spline coefficients 
c         at RZ mesh using array density_r_z
c-------------------------------------------------------------------
         write(*,*)'genray.f before splcoef_density_r_z'
         if(ixyz.eq.0) call splcoef_density_r_z
         write(*,*)'genray.f after splcoef_density_r_z'
c         stop 'genray.f after splcoef_density_r_z'
       endif
c-------------------------------------------------------------------
c     set zero to arrays power and current  for subroutine p_c_prof
c------------------------------------------
      do iray=1,nray         
         write(*,*)'3 iray,arzu0(iray),arru0(iray),arphiu0(iray)',
     &              iray,arzu0(iray),arru0(iray),arphiu0(iray)
      enddo

      write(*,*)'in genray ionetwo',ionetwo
      !if(ionetwo.eq.1) then
         call onetwoini ! set spower*(), scurr*() arrays to 0.
      !endif
c-------------------------------------------------------------------
      call output_con1
c---------------------------------------------------------------
c     the loop on all rays
c---------------------------------------------------------------
c     Initialize arrays for Runge-Kutta subroutine
      call arrays(ndim,deru,prmt,ihlf)
c-----------------------------------------------------------------
      if (isolv.eq.1) then
        write(*,*)'the runge-kutta solution of 6 hamiltonian equations
     1  with correction which gives hamiltonian conservation'
        write(*,*)'accuracy of the hamiltonian conservation
     1	in correction procedure epscor=',prmt(4)
      endif

      if (isolv.eq.2) then
        write(*,*)'the runge-kutta solution of 5 hamiltonian equations
     1  and the solution of the dispersion relation for the sixth
     2  ray variable'
      endif


CMPIINSERTPOSITION ENDBARRIER
c
c-----Construct names of .txt and .nc ray data files
      if( length_char(mnemonic).gt.124)
     1                         stop 'Adjust mnemonic in genray.f'
      write(filetxt,1001) mnemonic(1:length_char(mnemonic))
 1001 format(a,".txt")
      write(filenc,1002) mnemonic(1:length_char(mnemonic))
 1002 format(a,".nc")
c
c     open text file for 3d FP code
      if(rayop.eq."text" .or. rayop.eq."both") then
        i_=92
        open(i_,file=filetxt)
      endif
c-----preparing data for output mnemonic.txt and mnemonic.nc file
      write(*,*)'in genray.f before write3d1 nray',nray
      call write3d1(nray)      
      write(*,*)'in genray.f sub after write3d1 nrayl=',nrayl
      
      if(rayop.eq."netcdf" .or. rayop.eq."both") then
         !--------------------!
         call netcdf_create   ! NetCDF
         !--------------------!
         ! YuP 120505 Create file here, 
         ! then write data that is not dependent on rays
         write(*,*)'genray.f netcdf_create: netcdf file is created'  
         ! Now write peqdsk and x,y,z eqdsk mesh :
         call netcdf_eqdsk_data(filenc)
         if (n_wall.ne.0) then ! write wall and limiter coordinates 
            call wrtnetcdf_wall_limiter_data(filenc)
         endif
         call wrtnetcdf_plasma_prof(filenc) ! plasma profiles, dispersion roots
         write(*,*)'after wrtnetcdf_plasma_prof'
      endif
      
c----------------------------------------------------------------
c     calculate the OX optimal directions of the EC cone central ray
      write(*,*)'genray.f before  gr_OX_optimal_direction'
      if (i_ox.eq.1) then
         ifreq_write=1 !it is used in dinit_1ray
         if(ixyz.eq.0) then
           call gr_OX_optimal_direction(ndim) ! calls drkgs2
         else
           call gr_OX_optimal_direction_xyz(ndim) ! calls drkgs2_xyz
         endif
      endif
c----------------------------------------------------------------

CMPIINSERTPOSITION JUSTBEFORETHELOOP
      write(*,*)'genray.f before do 20 nray=',nray

      do iray=1,nray         
         write(*,*)'iray,arxu0(iray),aryu0(iray),arzu0(iray)',
     &              iray,arxu0(iray),aryu0(iray),arzu0(iray)
      enddo

      i_total_bad_initial_conditions=0
      power_launched=0.d0

       do 20 iray=1,nray !-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o

CMPIINSERTPOSITION STARTTHELOOP

c--------the loop over all rays        
         write(i1_,11)iray
 11      format('####  iray=',i3)
         write(*,'(///)')
         write(*,*)'### Starting iray ========================== ',iray
         write(*,*)'istart=',istart,'   ioxm=',ioxm

CMPIINSERTPOSITION PRINTHEADER

         nfreq=1
         i_geom_optic_loc=i_geom_optic
         do 25 ifreq=1,nfreq
            i_geom_optic=i_geom_optic_loc
            ifreq_write=ifreq ! to set in write.i
            call cpu_time(time_loop_ifreq_1)

            if(istart.eq.1) then
c ------------EC wave
              xst1=xstj(iray) 
              yst1=ystj(iray) 
              zst1=zstj(iray) 
              rst1=rstj(iray) 
              phist1=phistj(iray) 
              alfast1=alphaj(iray)
              betast1=betaj(iray)
              write(*,*)'genray before dinit_1ray ,zst1,rst1,phist1',
     1	             zst1,rst1,phist1,'alfast1,betast1',alfast1,betast1
              write(*,*)'in genray iray,powj(iray)',iray,powj(iray)
              powini=powj(iray)
               
              if(ixyz.eq.0)then
                 call dinit_1ray(zst1,rst1,phist1,alfast1,betast1,
     1                       cnteta,cnphi,u,iraystop)
              else
                 call dinit_1ray_xyz(xst1,yst1,zst1, 
     +              alfast1,betast1, cnteta,cnphi,
     +              u,iraystop) !->out
              endif
       
              i_bad_initial_conditions=0
              write(*,*)'genray.f after dinit_1ray',
     &        'i_bad_initial_conditions',i_bad_initial_conditions

              if (iraystop.eq.1) then
	        write(*,*)'iray=',iray,'iraystop=1
     1                  bad initial conditions'
                nrayelt=0
                i_bad_initial_conditions=1
	        goto 24
	      end if
	    endif !istart.eq.1

          if((istart.eq.2).or.(istart.eq.3))then
c ------------LH and FW waves >istart=2
c ------------ECR O_X mode conversion special case: istart=31 
              powini=powinilh(iray)
              write(*,*)'genray before dinit_1ray arxu0,aryu0,arzu0=',
     1	            arxu0(iray),aryu0(iray),arzu0(iray)
              write(*,*)'alfast1,betast1,arntheta(iray),arnphi(iray)',
     1                 alfast1,betast1,arntheta(iray),arnphi(iray)
              write(*,*)'genray.f before dinit_1ray ioxm ',ioxm
              if(ixyz.eq.0)then
                call dinit_1ray(arzu0(iray),arru0(iray),arphiu0(iray),
     1	              alfast1,betast1,arntheta(iray),arnphi(iray),
     +                u,iraystop)  !->out
              else
                call dinit_1ray_xyz(arxu0(iray),aryu0(iray),arzu0(iray),
     +                alfast1,betast1, arntheta(iray),arnphi(iray),
     +                u,iraystop) !->out
              endif
              write(*,*)'genray.f after dinit_1ray iraystop=',iraystop
              !pause
	    endif !istart.eq.2 .or. istart.eq.3

          i_bad_initial_conditions=0

	    if (iraystop.eq.1) then
               write(*,*)'iray=',iray,'iraystop=1
     1                  bad initial conditions'
               i_bad_initial_conditions=1
               nrayelt=0
               !pause
	       goto 24
	    end if
  
	    prmt(7)=prmt(6) ! Genray: prmt(7) initialized
     
c----------------------------------------------------------
c           call b() to calculate the small radius rho (inside b())
c           the result will be in common block  one.i
c-----------------------------------------------------------
            if(ixyz.eq.0)then
               bmod=b(u(1),u(2),u(3))
            else
cyup               bmod= bxyz(u(1),u(2),u(3)) ! what for? rho?
               den=dense_xyz(u(1),u(2),u(3),1) !-> get rho
            endif
CMPIINSERTPOSITION STARTRUNGEKUTTA
           write(*,*)' x,y,z,br,bz,bphi',u(1),u(2),u(3),br,bz,bphi
           write(*,*)' toteqd,psimag,psilim,dpsimax',
     +                 toteqd,psimag,psilim,dpsimax

            call cpu_time(time_before_rk)

            nstep_rk=1 ! initialize the number of Runge-Kutta time step 
            
            if (isolv.eq.1) then
c--------------------------------------------------------------
c              The Runge-Kutta solution of 6 hamiltonian equations
c              with correction which gives 
c              hamiltonian conservation with
c              accuracy epscor=prmt(4)
c--------------------------------------------------------------
               if(irkmeth.eq.0) then
c                4_th order Runge-Kutta method with constant time step
ccc                 call drkgs(prmt,u,deru,ndim,ihlf,rside1,outpt,aux)
                 stop 'set irkmeth=2 or 3; other values not supported'
               endif
              
               if(irkmeth.eq.1) then
c                 5-th order Runge-Kutta method with variable time step,
ccc                  call drkgs1(prmt,u,deru,ndim,ihlf,rside1,outpt,aux)
                 stop 'set irkmeth=2 or 3; other values not supported'
               endif
               
               if(irkmeth.eq.2) then !one of main options: drkgs2_xyz
c                 4th order Runge-Kutta with variable time step,
c                          time step can be reduced or enlarged
                  write(*,*)'genray.f before  drkgs2 ioxm=',ioxm
                  if(ixyz.eq.0)then
                     call drkgs2(prmt,  u,  deru,  ndim, ihlf,
     &                        rside1,  outpt,  aux,  i_output)
                  else
                     call drkgs2_xyz(prmt, u, deru, ndim, ihlf,
     &                        rside_xyz, outpt_xyz, aux, i_output)
                  endif
               endif
               
               if(irkmeth.eq.3) then ! main suggested option for now.
c                 4th order Runge-Kutta with variable time step,
c                          automatically selected. Need to set (example)
                 ! dL_step=1.d-3 ! [m]  max allowed change in r.
                 ! dN_step=1.d-2 ! max allowed change in refraction index.
                 ! The code will set the time step h=dt_code for integration
                 ! in such a way that the change in configuration space
                 ! is not larger than dL_step, and also
                 ! the change in refr. index |N| is not larger than dN_step
                     call drkgs_auto(prmt, u, deru, ndim, ihlf,
     &                      rside_xyz, outpt_xyz, aux, dL_step,dN_step)
               endif
            endif ! isolv.eq.1
          
            if (isolv.eq.2) then
c-------------------------------------------------------------------
c              The Runge-Kutta solution of 5 hamiltonian equations
c              and the solution of the dispersion relation for the sixth
c              ray variable'
c-------------------------------------------------------------------------
ccc              call rkb1(prmt,u,deru,ndim,ihlf,rsideb1,outptb1,aux)
c------------------------------------------------------------------------
                 stop 'set isolv=1; other values not supported'
            endif ! isov.eq.2
                          
 10         format(a/a)

            call cpu_time(time_after_rk)
            time_rk=time_after_rk-time_before_rk
            write(*,*)'time_initial conditions',
     &                 time_before_rk-time_loop_ifreq_1
            write(*,*)'time_rk',time_rk, '  nstep_rk=',nstep_rk 

CMPIINSERTPOSITION ENDRUNGEKUTTA
           
c------------------------------------------------------------------
c           creation of the file: genray.bin for xdraw
c           input data from common blocks gr.cb and write
               write(*,*)'genray.f before mk_graph'
               write(*,*)'genray. f (i_emission.eq.0) before mk_graph'
               call mk_graph(iray,nray,ifreq)
               call mk_gr3d(iray,nray)
c              creation of the file: npar.bin for xdraw
               call GRAPHnpr
               if (iwcntr.eq.1) then
c                 calculation of contours wb_c = consts
                  call mk_grapc(iray,iwopen,iwj)
               endif
c--------------------------------------------------------------------
c           calculation of power(array spower(NR)) and
c           current(array scurrent(NR)) radial profiles
c           as sum of profiles Sum(i=1,iray)
c---------------------------------------------------------------------
          !write(*,*)'ionetwo=',ionetwo
	    !if(ionetwo.eq.1) then
               call sonetwo ! add power*(i) to spower*(i), 
                            ! current*(i) to scurrent*(i) 
	    !endif
c------------------------------------------------------------
            
 24         continue !it was iraystop=1 bad initial condition for the ray
        
CMPIINSERTPOSITION ENDEMISSION
            if(i_bad_initial_conditions.eq.0) then
c--------------------------------------------------------------------
c             calculates total power launched along rays with 
c             good initial conditions
c--------------------------------------------------------------------
              power_launched=power_launched+powini
c---------------------------------------------------------------------
c             writing ray data to mnemonic.txt and/or 
c             saving data for mnemonic.nc
c-------------------------------------------------------------------  
	        call write3d
c------------------------------------------------------------
            else
c-------------calculate the total number of rays having 
c             i_bad initial conditions=1
               i_total_bad_initial_conditions= 
     &         i_total_bad_initial_conditions+1 
            endif !i_bad_initial_conditions 
c--------------------------------------------------------------

 25         continue ! ifreq       
          
 20         continue ! iray=1,nray !-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o
   
CMPIINSERTPOSITION JUSTAFTERTHELOOP

      write(*,*)'genray after 20 continue'
c     end iray loop

      if(rayop.eq."text" .or. rayop.eq."both") then
c--------close file for 3d FP code
         close(i_)
      endif

      write(*,*)'genray: powtot_e,powtot_i',powtot_e,powtot_i
      
      if(rayop.eq."netcdf" .or. rayop.eq."both") then
         call wrtnetcdf(1)
         write(*,*)'genray.f wrtnetcdf(1): define dimensions'  
         call wrtnetcdf(0)              
         write(*,*)'genray.f wrtnetcdf(0): write data' 
      endif
  
      if ((istart.eq.2).or.(istart.eq.3)) then
c--------write ray starting coordinates in filenc.nc file
         write(*,*)'genray.f before wrtnetcdf_grill_launch filenc',
     &   filenc
         call wrtnetcdf_grill_launch(filenc)
      endif

      if (istart.eq.1) then
c--------write ray starting coordinates in filenc.nc file
         write(*,*)'genray.f before wrtnetcdf_EC_launch filenc',
     &   filenc
         call wrtnetcdf_EC_launch(filenc)
      endif

      close(i1_)

c--------------------------------------------------------------------
c     calculation power and current density profiles at radius rho
c---------------------------------------------------------------------
      write(*,*)'genray.f before dnonetwo ionetwo',ionetwo
      if(ixyz.eq.0)then
        call dnonetwo         
      else
        !if(model_b.eq.0) 
        call dnonetwo_xyz
      endif
      
 1010 format('total power absorbed at reflections(watt)=',1pe14.6)
      write(*,1010) w_tot_pow_absorb_at_refl*1.d-7

c--------------------------------------------------------------------
c     Output power and current density profiles at radius rho
c--------------------------------------------------------------------
      if(ionetwo.eq.1) then
         write(*,*)'genray.f before mk_gronetwo ionetwo',ionetwo
         call mk_gronetwo
         call mk_gronetwo_1
      endif

      !write(*,*)'genray.f before wrtnetcdf_prof 1'
      call wrtnetcdf_prof(filenc,1)
      !write(*,*)'genray.f before wrtnetcdf_prof 0'
      call wrtnetcdf_prof(filenc,0)
      !write(*,*)'genray.f after wrtnetcdf_prof 0'

c--------------------------------------------------------------------
c     creates some data for drawing
c--------------------------------------------------------------------
      !write(*,*)'genray before mkgrtool' 
      if(itools.eq.1) call mkgrtool
      !write(*,*)'genray after mkgrtool' 
c--------------------------------------------------------------------
c     write dielectric tensor to mnemonic.nc file
c--------------------------------------------------------------------
      if (dielectric_op.eq.'enabled') then
         call wrtnetcdf_eps(filenc)
      endif

      call cpu_time(time_genray_2)
      write(*,1003) time_genray_2-time_genray_1
 1003 format('genray.f runtime [sec] ',1pd14.5)
      write(*,1004)
 1004 format('genray.f: Normal end of program')

CMPIINSERTPOSITION ENDMPI

      stop 
      end
c======================================================================
c======================================================================

c======================================================================
c======================================================================


CMPIINSERTPOSITION SUBS



c *****************************************************************
c ********* Tokamak data output to con1 ***************************
c *****************************************************************
      subroutine output_con1
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      open(30,file='con1')
      write(30,90) nreqd,nzeqd
90    format(8x,' the data of input file',/,
     . 2x,'the numbers of points in the r direction:',i3,/,
     . 2x,'the numbers of points in the z direction:',i3)
      write(30,11) rdimeqd*r0x,zdimeqd*r0x,reqd*r0x,redeqd*r0x,
     1             zmideqd*r0x
11    format(2x,'the full-width of rectangle: dx=',f12.6,' m',/,
     .  30x,' dy=',f12.6,' m',/,
     .  2x,'the major radius of the torus:',f12.6,' m',/,
     .  2x,'the major radius of the inner edge',/,
     .  2x,'of rectangular grid:',f12.6,' m',/,
     .  2x,'the vertical shift up-down symmetry plane:',f12.6,' m')
      write(30,12) rma*r0x,zma*r0x,psimag*b0*r0x**2,
     1             psilim*b0*r0x**2,beqd*b0
12    format(2x,'the major radius of magnetic axis:',f12.6,' m',/,
     .  2x,'the vertical height of magnetic axis:',f12.6,' m',/,
     .  2x,'the poloidal flux function values',/,
     .  2x,'at the magnetic axis:',f12.6,/,
     .  2x,'and the last closed flux surface:',f12.6,/,
     .  2x,'the toroidal magnetic field',/,
     .  2x,'at major radius of the torus:',f12.6,' t')
      write(30,13) toteqd
13    format(2x,'the toroidal curent:',e17.8,' a')
      close(30)

      return
      end




      subroutine onetwoini
c--------------------------------------------
c     set zero to arrays power and current
c     for subroutine p_c_prof
c YuP[10-1-2014] Added: 
c     definition of rho_bin(), rho_bin_center(), binvol, etc.
c------------------------------------------
      implicit none
c      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'onetwo.i'
      include 'five.i' ! rmax,...
      include 'three.i' ! eq.grid
      include 'rho.i' ! voltot, areatot,...
      include 'fourb.i' ! contains req,zeq,xeq,yeq grids
c-----locals
      integer i,j,l, kk
      real*8 den, rmx, hrho, rholeft, rhoright, x, y, z
      real*8 dense_xyz, rhov, rhos, rho_lrho
      real*8 PSI_in, rhox, dz_line, ravg, volum(NR),areac(NR)
      real*8 psi_rho
      real*8 dRgrid, dZgrid, Rgrid_mn,Rgrid_mx, Zgrid_mn,Zgrid_mx
      
      real*8 Rline(nzeqda),Zline(nzeqda) ! for tracing field line
      ! and calculating binvol volume; could be other size, but seems ok.
      integer iline,nline
      
      real*8 psi0,tini,tm,htini,costet,sintet,zt0,rt0,t0
      integer ierr

      write(*,*)'onetwoini NR,NRA', NR,NRA
      do i=1,NR
         spower(i)=0.0d0
         spower_e(i)=0.0d0
         spower_i(i)=0.0d0
         spower_cl(i)=0.0d0
         scurrent(i)=0.0d0       
      enddo
      do i=1,NR-1
         s_cur_den_parallel(i)=0.d0
      enddo
      do kk=1,nbulk
         do i=1,NR-1
            spower_s(i,kk)=0.0d0
         enddo
      enddo

c YuP[10-1-2014] Added: definition of rho_bin(), rho_bin_center(), binvol, etc.
c These 1D arrays are stored in onetwo_no_nml.i, and reused when needed.
c Before, many of them where calculated locally in different subroutines.     
c-----create (for profiles):
c     small radius array rho_bin_center(i) for bin(i) centers
      !if(model_b.ne.0)then
      ![09-24-2014] Now for any model_b:
      ! Find the value of rho at wall (can be >1; open flux surfaces)
      ! But be sure not to exceed the grid:
      rmx= min(wall_rmax, xeqmax,yeqmax,abs(xeqmin),abs(yeqmin))
      
      ! Find rhowall
      if(model_b.eq.2 .or. model_b.eq.3)then ! mirror machine: 
      !pol.flux goes to INF at coils;
      !flat outside of R>R_b0 for mirror1 model
         den=dense_xyz(rmx,0.d0,0.d0,1)
         rhowall=rho
         write(*,*) 'onetwoini: rmx, rho(rmx)=', rmx, rhowall
      else ! model_b.eq.0,1,4
         ! Scan equilibrium B grid, find largest rho - set rho_bin(NR)
         ! WARNING: if the grids are too refined, it may take for ever...
         ! To save cpu time, switched to scan over (Y,Z) grid only
         ! assuming toroidal symmetry (the problem could be if ray
         ! is started at a corner of (X,Y) grid, but this is not likely).
         write(*,*) nxeqd,xeq(1),xeq(nxeqd)
         write(*,*) nyeqd,yeq(1),yeq(nyeqd)
         write(*,*) nzeqd,zeq(1),zeq(nzeqd)
         rhowall=0.d0
c         do i=1,nxeqd
c           x=xeq(i) != xeqmin+dstep*(i-1) ! YuP [xeqmin; xeqmax]
         x=0.d0 !To reduce cpu time, scan only within (Y,Z) grid, at x=0
         do j=1,nyeqd
           y=yeq(j) != yeqmin+dstep*(i-1) ! YuP [yeqmin; yeqmax]
         do l=1,nzeqd
           z=zeq(l) != zeqmin+dstep*(i-1) ! YuP [zeqmin; zeqmax]
           den=dense_xyz(x,y,z,1) !-> get rho 
           rhowall=max(rhowall,rho) ! rhowall can be >1
         enddo
         enddo
c         enddo
      endif
      
      hrho=rhowall/(NR-1)  ! rhowall can be >1
      write(*,*) 'onetwoini: rhowall=', rhowall
      !pause !!!
      !else
      !  hrho= 1.0/(NR-1) ! gives rho=[0;1]
      !endif
      do i=1,NR
         rho_bin(i)=hrho*(i-1) ! [0;rhowall]    can be >1
      enddo
      do i=1,NR-1
         rho_bin_center(i)=0.5d0*(rho_bin(i)+rho_bin(i+1))
      enddo

      volum(:)=0.d0
      areac(:)=0.d0
      if(model_b.eq.3)then ! trace field lines, define volume
         do i=1,NR
             ! Use the volume of a solid of revolution around Z-axis:
             ! Volum= pi*INTEGRAL(R^2*dZ) where R==R(Z) field line
             ! The area is Area= INTEGRAL(R*dZ)
             nline= size(Rline)
             rhox= rho_bin(i)
             PSI_in= dpsimax*psi_rho(rhox)
             !Note: psi_rho=psimag+(psilim-psimag)*rhox*rhox
             !In Genray, psilim>psimag, so that psi_rho is ascending
             !from R=0 to R=edge.
             !But in the model, PSI can be ascending or descending,
             !depending on b00_mirror sign (=sign of Bz).
             !To have the proper input value of PSI_in (proper sign),
             !psi_rho is multiplied by dpsimax.
             call eq_mirror1_field_line(PSI_in,Rline,Zline,nline)
             !line is traced from -0.5*zbox_mirror to +0.5*zbox_mirror
             !The formula in the above subr. fails at R=R_b0 or beyond,
             !but that region is not important (outside of plasma)
             do iline=2,nline
                dz_line= Zline(iline)-Zline(iline-1)
                ravg= 0.5*(Rline(iline)+Rline(iline-1)) ! R(Z) at field line
                volum(i)= volum(i) +dz_line*ravg*ravg ! INTEGRAL(R^2*dZ)
                areac(i)= areac(i) +dz_line*ravg
             enddo
             volum(i)=abs(volum(i)*pi)*1.e6 !1e6 to convert to cm^3
             areac(i)=abs(areac(i))*2.e4 !2 here is for left|right symmetry
             !(line is traced over |)) region; area is over ((|)) region).
             !In other words, the field line was traced only over right 
             !region r>rmag, not the whole closed surface as in tokamaks. 
         enddo
      endif    

      write(*,*)'rmx,zma,zmax,rmax,psilim=',
     +    rmx,zma,zmax,rmax,psilim

      if(model_b.eq.0 .and. eqdsktype.eq.'mirror')then !YuP[11-2016]
         ! trace field lines, define volume
         ! Use the volume of a solid of revolution around Z-axis:
         ! Volum= pi*INTEGRAL(R^2*dZ) where R==R(Z) field line
         ! The area is Area= INTEGRAL(R*dZ)
         ! The design is similar to YuP[11-2016] in gr2new
         ! but here the set of psi0 "surfaces" is different.
         ! Here, they correspond to the set of rho_bin(j).
         do j=2,NR  ! start with j=2;   j=1 is R=0 line, so volume(1)=0
             nline= size(Rline)
             rhox= rho_bin(j)
             psi0= psi_rho(rhox)
             do iline=1,nline
                Zline(iline)= zeqmin+(iline-1)*(zeqmax-zeqmin)/(nline-1) ![zeqmin; zeqmax]
                tini=0.d0 ! [m] initial t along R(t)=rma+t
                tm=rmx   ! [m] largest t (see zrcntrbin)
                htini=tm*0.002d0 ! initial t-step
                ! For a given psi0, and given Z=zpsi, find R=rpsi
                ! such that psif_xyz(x=Rline,y=0,z=Zline) = psi0.
                ! For searching of rpsi, we use same subr.zrcntrbin  
                ! as for tokamaks, but with a trick:
                ! The scanning direction now is along Z=const lines,
                ! (rather than constant pol.angle)
                ! so that only R is scanned, for each given Z.
                ! The values of Z are set by Zline(iline), 
                ! and they are the same for any psi0.
                ! This is based on assumption that all field lines
                ! cover the whole range of Z=[zmin;zmax],
                ! which is usually the case in a mirror machine.
                ! Index j labels different field lines, 
                ! and index i labels points along field lines.
                !binary method for solution of equation psi(r(t),z(t))=psi0
                rma=0.d0
                zma=Zline(iline) !to make proper input for zrcntrbin
                costet=1.d0 ! which means scanning along rt1= rma+1*t1
                sintet=0.d0 ! which means scanning along zt1= zma+0*t1
                call zrcntrbin(tini,psi0,costet,sintet,tm,htini,
     1          ierr,zt0,rt0,t0)
                if (ierr.eq.1) then
                  Rline(iline)=rt0
c       write(*,'(a,2i4,2e12.4)')'j,iline, Rline,Zline=', 
c     +  j,iline,Rline(iline),Zline(iline)
                else
                  ! Failed to find the solution of psi(r(t),z(t))=psi0.
                  ! It may happen in the far corners of (R,Z) equil grid:
                  ! the field lines that start at the near-corner
                  !  at Z=zmin (or zmax) would go outside the grid
                  ! at smaller Z. So the volume cannot be defined 
                  ! in this case. But it's not really important.
                  !write(*,*)'Rline(iline),Rline(iline-1),rt0=',
     +            ! Rline(iline),Rline(iline-1),rt0
c                  write(*,*)' onetwoini: zrcntrbin gave ierr ',ierr
c                  write(*,*)' it is impossible to find the flux surface'
c                  write(*,*)' with rho_bin=',rhox, '  psi0=', psi0
                   Rline(iline)=rmax ! just set it to the largest R.
c       write(*,'(a,2i4,2e12.4)')'j,iline, Rline,Zline=', 
c     +  j,iline,Rline(iline),Zline(iline)
c                  write(*,*)'----------------------------------------'
                endif
                if(iline.ge.2)then
                dz_line= Zline(iline)-Zline(iline-1)
                ravg= 0.5*(Rline(iline)+Rline(iline-1)) ! R(Z) at field line
                volum(j)= volum(j) +dz_line*ravg*ravg ! INTEGRAL(R^2*dZ)
                areac(j)= areac(j) +dz_line*ravg
                endif
             enddo ! iline == step along field line
             volum(j)=abs(volum(j)*pi)*1.e6 !1e6 to convert to cm^3
             areac(j)=abs(areac(j))*2.e4 !2 here is for left|right symmetry
             !(line is traced over |)) region; area is over ((|)) region).
             !In other words, the field line was traced only over right 
             !region r>rmag, not the whole closed surface as in tokamaks. 
         enddo ! j=1:NR  step along rho_bin
         ! restore:
         rma=0.d0 ! mirror machine
         zma=0.d0 ! mirror machine
      endif    



      write(*,*)'voltot*1.d6, areatot*1.d4=',voltot*1.d6,areatot*1.d4
      do i=1,NR-1
         rholeft= rho_bin(i)   !YuP: was: (i-1)/dble(NR-1)
         rhoright=rho_bin(i+1) !YuP: was: rholeft+ 1.d0/dble(NR-1)
         if(model_b.eq.0)then
         
            if((eqdsktype.eq.'tokamak').or.(eqdsktype.eq.'TAE'))then
              !write(*,*)'rhov,rhos=',rhov(rhoright),rhos(rhoright)
              if(rhoright.le.1.d0) then
               binvol(i)=voltot*(rhov(rhoright)**2-rhov(rholeft)**2)
     1                   *1.d6
               binarea(i)=areatot*(rhos(rhoright)**2-rhos(rholeft)**2)
     1                   *1.d4
               pollen(i)=
     1           totlength*50.d0*(rho_lrho(rhoright)+rho_lrho(rholeft))
              else ! rho>1
               ! "Open" flux surfaces. 
               ! Difficult to define volume in general case.
               ! Set it equal to the volume of the previous bin:
               binvol(i)= binvol(i-1)
               binarea(i)=binarea(i-1)
               pollen(i)= pollen(i-1)
               ! Not accurate, but at least  
               ! the collisional power density p_cl
               ! can be calculated in this region (rho>1).
              endif
            endif
            
            if(eqdsktype.eq.'mirror')then ! and model_b=0 here
               ! Open field lines
               binvol(i)=  (volum(i+1)-volum(i)) ! cm^3 
               binarea(i)= (areac(i+1)-areac(i)) ! cm^2    
            endif
            
         elseif(model_b.eq.3)then
            binvol(i)=  (volum(i+1)-volum(i)) ! cm^3 
            binarea(i)= (areac(i+1)-areac(i)) ! cm^2    
         else
            binvol(i)=1.
            binarea(i)=1.
            ! Generally, no flux surfaces and (not even "open" surfaces)
            ! (example: magnetic mirror machine).
            ! The powden_* arrays  will store 
            ! the locally absorbed power [Watt] for a given 
            ! rho-bin (which is defined based on model_rho_dens)
cyup           write(*,*)'binvol is not defined for model_b>0'
         endif
         write(*,'(a,i4,5e12.4)')
     +   ' onetwoini: i, rho_bin_center, binvol,volum, binarea,areac=',
     +        i,rho_bin_center(i),binvol(i),volum(i),binarea(i),areac(i)
      enddo   
      write(*,*)'sum(binvol),sum(binarea)=',sum(binvol),sum(binarea)
      write(*,*)'volum(NR),areac(NR)=     ',volum(NR),areac(NR)
      !pause
      
      ! YuP[Nov-2014] 
      ! For power deposition profiles over (R,Z) rectangular grid.
      ! (R,Z) grid for power deposition profiles:
      Rgrid_mn= 0.d0
      Rgrid_mx= rmx
      dRgrid= (Rgrid_mx-Rgrid_mn)/(NRgrid-1)
      ! Note: for Z-direction, the values of zeqmax and |zeqmin|,
      ! which are based on eqdsk grid, could be quite large,
      ! example: TAE-C2 FRC. 
      ! The waves are not likely to be propagating
      ! at Z far beyond separatrix.
      ! So, for now, we choose zmin/zmax corresponding to the 
      ! lower and upper X-points,  plus a bit (say, 10%):   
      Zgrid_mn= zmin-0.10*(zmax-zmin)   
      Zgrid_mx= zmax+0.10*(zmax-zmin)   
      !For full range, use: 
      Zgrid_mn=min(zmin,zeqmin)
      Zgrid_mx=max(zmax,zeqmax)
      dZgrid= (Zgrid_mx-Zgrid_mn)/(NZgrid-1)
      do i=1,NRgrid
         Rgrid(i)= Rgrid_mn +dRgrid*(i-1)
      enddo
      do i=1,NZgrid
         Zgrid(i)= Zgrid_mn +dZgrid*(i-1)
      enddo
      write(*,*) 'Rgrid min/max:', Rgrid(1), Rgrid(NRgrid)
      write(*,*) 'Zgrid min/max:', Zgrid(1), Zgrid(NZgrid)
      
      Write(*,*)'onetwoini: NRgrid, dRgrid=', NRgrid, dRgrid
      Write(*,*)'onetwoini: NZgrid, dZgrid=', NZgrid, dZgrid
      
      if(dZgrid .gt. 1.5*dRgrid)then
      write(*,*)'onetwoini: dZgrid>>dRgrid. Consider increasing NZgrid'
      endif
      
      if(dRgrid .gt. 1.5*dZgrid)then
      write(*,*)'onetwoini: dRgrid>>dZgrid. Consider increasing NRgrid'
      endif

      !pause
      
      return
      end
      
      
      
***********charnumb*******************************************
*     The transformation of the integer j to character chj
*     It is assumed that the number of the decimal numbers in j<=8
******************************************************
c     input : j   is integer 
c     ATTENTION: the number of decimal positions in j should be<=8
c     output: chj is character*8 
c-------------------------------------------------------------
      subroutine charnumb(j,chj)
      integer nj
      parameter (nj=8) ! the max number of the decimal numbers in j 
      character chj*8
      character ikch
      dimension ikch(nj)
      ich0=ichar('0')

      k=1
      jk=j/10

      do while(jk.ne.0) 
         k=k+1
	 jk=jk/10
      enddo	 
c      write(*,*)'k is the number of the decimal nimbers in j',k
      if(k.gt.nj)then
        write(*,*)'in charnumb k>nj'
	stop
      endif
      do i=1,nj
         ikch(i)='0'
      enddo

      jk=j
      do i=1,k
         jkn=jk/10
	 ik=jk-jkn*10
	 jk=jkn
c	 ikch(k) is the decimal number in the j in the k position  
c        The position numeration is from the right side
	 ikch(i)=char(ik+ich0)
      enddo
      chj=ikch(8)//ikch(7)//ikch(6)//ikch(5)//ikch(4)//ikch(3)//
     +ikch(2)//ikch(1)
c      write(*,*)'in charnumb chj=',chj
      return
      end




c*************************contrb1 *************************************
c  It calculates contours coordinates for the open contours R=R(z):   *
c          modb(r,z)=const and y_(r,z,phi,2)=1./n, n=1,2,...          *      			      *
c          arrays ry(j,i) zpsi(i)  				      *
c          j=1,100 npsi(number of countours n=1,100)        	      *
c          i=0,20 ( a number of points in z direction)                *
c---------------------------------------------------------------------
c  input data are in common one.i,five.i                	      *
c---------------------------------------------------------------------*
c* Output: arrays (r,z,) for the given values 1/Y_i=n to               *
c  files: outputy=n.dat         				      *
c**********************************************************************

      subroutine contrb1
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
      include 'gr.i'
      include 'write.i'
      dimension	 zy(20),ry(100,20)
      double precision ias1r
      character chj*8
      character outdaty*20
      iw_j=iwj
c      write(*,*)'in contrb1 iw_j',iw_j
      phi=0.d0

      nnr=40
      nz=40

      hr=rdimeqd/(nnr-1.d0)
      hz=zdimeqd/(nz-1.d0)
      bmin=b(zrmin,rmin,phi)
c      write(*,*)'0 bmin',bmin
     
      open(5,file='wdwci.doc')
      do i=1,nnr+1
	r=4.d0+0.05d0*(dfloat(i)-1.d0)
        do j=1,nz
	  z=-2.d0+0.1d0*(dfloat(j)-1.d0)
          bmod=b(z,r,phi)
	  dyci=1.d0/y(z,r,phi,iw_j)
	  write(5,10)r,z,dyci
	enddo
      enddo
      close(5)
     
      nnr=40
      nz=30
      nz=71
      hr=rdimeqd/(nnr-1.d0)
      hz=zdimeqd/(nz-1.d0)
      hz=14.d0/(nz-1.d0)
c      hr=0.05d0
c      nz=0.1d0
      open(4,file='btot.doc')
      do i=1,nnr
	r=redeqd +hr*(i-1.d0)
c	r=4.d0+hr*(i-1.d0)
        do j=1,nz
          z=zmideqd+hz*(j-1)-zdimeqd*0.5
          z=-7.0d0+hz*(j-1)
c	  write(*,*)'in contrb1 i,z,j,r',i,z,j,r
          btot=b(z,r,phi)
c	  write(*,*)'in contrb1 i,z,j,r, btot',i,z,j,r,btot
	  write(4,10)r,z,btot
        enddo
      enddo
      close(4)

      btotl=b(zrmin,rmin,phi)
      bmod=btotl
      yl_i=y(zrmin,rmin,phi,iw_j)
      btotr=b(zrmax,rmax,phi)
      bmod=btotr
      yr_i=y(zrman,rmax,phi,iw_j)
      write(*,*)'in contourb btotl btotr',btotl,btotr
      write(*,*)'in contourb yl_i yr_i',yl_i,yr_i
      nll=int(1./yl_i)
      nlr=int(1./yl_i)+1
      nrl=int(1./yr_i)
      nrr=int(1./yr_i)+1
      write(*,*)'in contrb1 nll,nlr,nrl,nrr',nll,nlr,nrl,nrr
      nz=20
      hz=(zmax-zmin)/dfloat(nz)
      epsy=1d-5 ! accuracy 
      phi=0.d0
      do i=1,nz
        z=zmin+i*hz
	zy(i)=z
	do j=nlr,nrl
c          determination of ry(j,i) where
c          y(z,ry,phi,iw_j)=1/j (using the binary method)
	   tl=rmin
	   tr=rmax
	   c=1.d0/float(j)
           do while ((tr-tl).gt.epsy)
             t=tl+(tr-tl)*0.5d0
             r=t
	     bmod=b(z,r,phi)
             y1=y(z,r,phi,iw_j)-c
             rtr=tr
	     bmod=b(z,rtr,phi)
             y2=y(z,rtr,phi,iw_j)-c
             if ((y1*y2).gt.0) then
               tr=t
             else
               tl=t
             end if
           end do
c          -----------------------------------------------
c          end of the binary method
c          -----------------------------------------------

           ry(j,i)=r
         enddo !j	
      enddo !i	
c         CALL ASSIGN("assign -F f77 -N ieee u:84",ier)
      write(*,*)'contours 1/Y_i=n will be plotted for n=',nlr,'...',nrl 
 10   format(3(1pe11.3))
 11   format(2(1pe11.3))
 12   format(34(1pe11.3))
      do j=nlr,nrl
        call charnumb(j,chj)
        outdaty=chj//'.dat'
	write(*,*)outdaty
	open(1,file=outdaty)
        do i=1,20
	  bmod=b(zy(i),ry(j,i),phi)
c	  write(*,*)'j,y=1/j,y(zy(i),ry(j,i),phi,iw_j),zy(i),ry(j,i)',
c     6	  j,1.d0/dfloat(j),y(zy(i),ry(j,i),phi,iw_j),zy(i),ry(j,i)
     
c---------------- idx derivativs order 0.ge.idx.le.3---------------
          idx=0
          ipx=ip
          ipx4=ip+4
          rrr=ry(j,i)
          zzrp=ias1r(trlimp,ipx,ipx4,cxlimp,idx,rrr)
          zp=zzrp
          ipx=im
          ipx4=im+4
          zzrm=ias1r(trlimm,ipx,ipx4,cxlimm,idx,rrr)
          zm=zzrm
c         write(*,*)'in bound r,zm,z,zp',rrr,zm,zy(i),zp
          if ((zy(i).le.zp).and.(zy(i).ge.zm)) then
c           write(*,*)'in write r,zm,z,zp',rrr,zm,zy(i),zp
c           write(*,*)'ry(j,i)*100,zy(i)*100)'
c           write(*,*)ry(j,i)*100,zy(i)*100

c            write(1,10)ry(j,i),zy(i),dfloat(j)
c            write(2,11)ry(j,i),zy(i)
            write(1,11)ry(j,i)*100.d0,zy(i)*100.d0 ! *100 to get cm
         WRITE(84) REAL(ry(j,i)*100),REAL(zy(i)*100),
     1   REAL(XT(3,NP+1)),REAL(YT(3,NP+1)),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt))
         WRITE(83,5)ry(j,i)*100,zy(i)*100,
     1   XT(3,NP+1),YT(3,NP+1),
     2   ws(nrayelt),delpwr(nrayelt),rez(nrayelt),
     3   spsi(nrayelt),wnpar(nrayelt),wnper(nrayelt),
     4   salphal(nrayelt)
          endif
	enddo
        WRITE(84)
        WRITE(83,2)
 2    format(/)
	close(1)
      enddo
 5    format(11(1pe10.3))
      return
      end





c*************************contrb2 *************************************
c  It calculates contours coordinastes for contours:		      *
c          modb(r,z)=const and y_(r,z,phi,2)=1./n, n=1,2,...          *
c          for the case then the closed contours exist.                *
c          arrays ry(j,i) thetac(i)  				      *
c          j=1,100 npsi(number of countours n=1,100)        	      *
c          i=0,nthetac (number of points in thetac direction)      	      *
c---------------------------------------------------------------------
c  input data are in common one.i,five.i                	      *
c---------------------------------------------------------------------*
c* Output: arrays (r,z) for the give values 1/Y_i=n to               *
c  files: outputy=n.dat         				      *
c**********************************************************************

      subroutine contrb2
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
      include 'gr.i'
      include 'write.i'
      parameter (nthetac=20, nthetac1=nthetac+1)
      dimension	 thetac(nthetac),zy(100,nthetac1),ry(100,nthetac1)
      double precision ias1r
      character outdaty*20,chj*8
      iw_j=iwj
      phi=0.d0
c------------------------------------------------------
c     creation of the btot.doc file with the coordinates (r,z) and
c     the values of mod(b(r,z))
c     calculation the coordinates (rbmin zbmin) of the point
c     inside the plasma where btot(rbmin,zbmin)=bmin
c------------------------------------------------------
c      write(*,*)'in contrb2 before nz=40'
      nnr=40
      nz=40
      hr=rdimeqd/(nnr-1.d0)
      hz=zdimeqd/(nz-1.d0)
      bmin=b(zrmin,rmin,phi)
c      write(*,*)'0 bmin',bmin
      open(4,file='btot.doc')
      open(5,file='wdwci.doc')
      do i=1,41
	r=4.d0+0.05d0*(dfloat(i)-1.d0)
        do j=1,40
	  z=-2.d0+0.1d0*(dfloat(j)-1.d0)
          bmod=b(z,r,phi)
	  dyci=1.d0/y(z,r,phi,iw_j)
	  write(5,10)r,z,dyci
	enddo
      enddo
      close(5)
      close(4)
c      stop
      do i=1,nnr
	r=redeqd +hr*(i-1.d0)
        do j=1,nz
          z=zmideqd+hz*(j-1)-zdimeqd*0.5
          btot=b(z,r,phi)
cc          bmod=b(z,r,phi)
cc	  dyci=1.d0/y(z,r,phi,iw_j)
c	  write(*,*)'in contrb2 i,z,j,r,btot,bmin',i,z,j,r,btot	,bmin
c---------------- idx derivativs order 0.ge.idx.le.3---------------
          idx=0
          ipx=ip
          ipx4=ip+4
          rrr=r
          zzrp=ias1r(trlimp,ipx,ipx4,cxlimp,idx,rrr)
          zp=zzrp
          ipx=im
          ipx4=im+4
          zzrm=ias1r(trlimm,ipx,ipx4,cxlimm,idx,rrr)
          zm=zzrm
          if ((z.le.zp).and.(z.ge.zm)) then
	    if(bmin.gt.btot) then
	      bmin=btot
	      rbmin=r
	      zbmin=z
	    endif
	  endif
	  write(4,10)r,z,btot
cc	  write(5,10)r,z,dyci
        enddo
      enddo
      close(4)
c      close(5)
      write(*,*)'in contrb2 bmin, zbmin,rbmin',bmin, zbmin,rbmin
c      stop
      btotmin=b(zbmin,rbmin,phi)
      bmod=btotmin
      ybmin_i=y(zbmin,rbin,phi,iw_j)
c      write(*,*)'in contrb2 ybmin_i',ybmin_i
c------------------------------------------------------
c     determination of yl_i(on the left side of plasma)
c     determination of yr_i(on the right side of plasma)
      btotl=b(zrmin,rmin,phi)
      bmod=btotl
      yl_i=y(zrmin,rmin,phi,iw_j)
      btotr=b(zrmax,rmax,phi)
      bmod=btotr
      yr_i=y(zrman,rmax,phi,iw_j)
c      write(*,*)'in contrb2 btotl btotr',btotl,btotr
c      write(*,*)'in contrb2 yl_i yr_i',yl_i,yr_i
c---------------------------------------------------------
c     determination of the ysep=1/Y_ci on the the last closed contour  
      if(yl_i.gt.yr_i) then
         ysep=yl_i
      else
         ysep=yr_i
      endif
c      write(*,*)'in contrb2 ysep',ysep
      ny_min=int(1.d0/ysep+2.d0)
      ny_max=int(1.d0/ybmin_i)
c      write(*,*)'in contrb2 ny_min,ny_max',ny_min,ny_max
c---------------------------------------------------------
c  Calculations coordinates of contours 1/Yci=n				      *
c          r(n,thetac)   z(n,thetac)     			      
c          arrays ry(j,i) zy(j,i)  				      
c          j=ny_min,ny_max(number of contours )	      
c          i=1,nthetac+1(number of points in the poloidal angle)  
      pi=4.d0*datan(1.d0)
      hteta=2.d0*pi/dble(nthetac)
      epsy=1d-3 ! accuracy 
c------------------------------------------------------
c     
c     2)we will create the limiter points using the close flux
c      surface psi(r,z)=psilim*psifactr, here
c      psifactr is a parameter (it must be .le.1) to avoide the
c      problems with the nonmonotonic psi function near the separatrix.
c      psifactr is given in genray.in file (It is in common/one/)
c     ------------------------------------
c      psilim=psimag+(psilim-psimag)*psifactr
c     
c-------------------------------------------------------
c        ipsi=1 !  to calculate contours
c        ipsi=0 !  to read contours data from file:psi.bin
c       if (ipsi.eq.0) then
c         open(1,file='psi.bin',form='unformatted',status='old')
c         do i=1,npsi
c	  do j=1,nteta1
c             read(1)zpsi(i,j),rpsi(i,j)
c	  end do
c         end do
c	 close(1)
c	 go to 200
c       endif
c------------------------------------------------------------------
      do 100 i=1,nthetac
         theta=hteta*(dble(i)-0.5d0)
         sintet=dsin(theta)
         costet=dcos(theta)
	 tini=0.d0
	 tm=3.5d0
	 htini=tm*0.02d0
         maxiter=10
	 htmin=1.d-3
         do 20 j=ny_max,ny_min,-1
           write(*,*)'genray.f contrb2 number of contour Y=1/n j=',j      
           n0=j 
c------------------------------------------------
c          binary method for solution of equation 1/Y(r(t),z(t))=n0
           call zrcntr2(tini,n0,costet,sintet,tm,htini,
     1     ierr,zt0,rt0,t0,zbmin,rbmin,epsy)
c------------------------------------------------
	   tini=t0
	   if (ierr.eq.1) then
              zy(j,i)=zt0
              ry(j,i)=rt0
              write(*,*)'j,i,zy(j,i),ry(j,i)',j,i,zy(j,i),ry(j,i)
           else
              write(*,*)'zrcontor gave ierr ',ierr
              write(*,*)'it is impossible to find the 1/y surface '
	      write(*,*)' with n',n0,'j=',j
	      write(*,*)' i=',i,'theta',theta
	      write(*,*)' with n0-1',n0-1
	      write(*,*)' try to change ny_min=arpsi(j-1) 
     1        or reduce the factor psifactr in subr. equilib'
       	      stop
	   endif
20       continue
100   continue

      do 40 j=ny_max,ny_min,-1
           zy(j,nthetac1)=zy(j,1)
           ry(j,nthetac1)=ry(j,1)
40    continue
c----------------------------------------------------------
c     if ipsi=1 then continue,write file y.bin
      open(3,file='y.bin',form='unformatted')
      do i=1,npsi
        do j=1,nthetac1
           write(3)zy(i,j),ry(i,j)
        end do
      end do
      close(3)
c------------------------------------------------------------
c     if ipsi=1 then continue,ipsi=0 then read file psi.bin
200    continue
c-------------------------------------------------------------
 10   format(3(1pe11.3))
 11   format(2(1pe11.3))
 12   format(34(1pe11.3))
      do j=ny_min,ny_max
        call charnumb(j,chj)
        outdaty=chj//'.dat'
c	write(*,*)outdaty
	open(1,file=outdaty)
c	write(j,*)'r ',j
c	goto 17 !!!!
        do i=1,nthetac1
	  bmod=b(zy(j,i),ry(j,i),phi)
	  write(*,*)'j,y=1/j,y(zy(i),ry(j,i),phi,iw_j),zy(j,i),ry(j,i)',
     6	  j,1.d0/dfloat(j),y(zy(j,i),ry(j,i),phi,iw_j),zy(j,i),ry(j,i)
     
          WRITE(84) REAL(ry(j,i)*100.),REAL(zy(j,i)*100.),
     1    REAL(XT(3,NP+1)),REAL(YT(3,NP+1)),
     2    REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3    REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4    REAL(salphal(nrayelt))
          WRITE(83,5)ry(j,i)*100,zy(j,i)*100,
     1    XT(3,NP+1),YT(3,NP+1),
     2    ws(nrayelt),delpwr(nrayelt),rez(nrayelt),
     3    spsi(nrayelt),wnpar(nrayelt),wnper(nrayelt),
     4    salphal(nrayelt)
	enddo !i
        WRITE(84)
        WRITE(83,2)
 2    format(/)
 17   continue
      close(1)
      enddo !j
 5    format(11(1pe10.3))
c      write(*,*)'end of contrb2'
c      stop
      return
      end
      
      
      
      
c*************************zrcntr2*********************************** *
c  It calculates  contours coordinates of the countor point 	       *
c   r(n,teta)   z(n,teta) for the given:			       *
c   1/y_ci(z,r)=n0 and poloidal angle teta0(in radians)                *
c   It using the binary methode                                        *
c---------------------------------------------------------------------
c  input data are in common blocks /three/,/five/,/gr.cb/	       *
c  tini -initial (left) value of the t( ray parameter )		       *
c  n0   -given value of the 1/y_ci       			       *
c  costet0,sintet0 -for the given poloidal angle		       *
c  rbmin major radius of point mod(b)=min (inside the plasma)	       *
c  zbmin Z coordinate of point mod(b)=min (inside the plasma)	       *
c  tm   -maximal value of t (right)				       *
c  htini   ininial step of the t 	 			       *
c  epsy-accuracy of equation solution (=max(abs( 1/y_n-1/y_n+1))
c----------------------------------------------------------------------*
c  Output: ,zt0(n0,teta0),rt0(n0,teta0) and parameter t0               *
c           rt0=rma+t0*costet0 ,  zt0=zma+t0*sintet0		       *
c  ierr -index if the solution was obtained =1	else=0		       *
c**********************************************************************
      subroutine zrcntr2(tini,n0,costet0,sintet0,tm,htini,
     1 ierr,zt0,rt0,t0,zbmin,rbmin,epsy)
      implicit real*8 (a-h,o-z)
      include 'param.i'
      include 'three.i'
      include 'one.i'
      iw_j=iwj
      ierr=1
      ht=htini
      phi=0.d0
      n0=-n0
c      write(*,*)'in zrcntr2 n0,costet0,sintet0',
c     1 n0,costet0,sintet0
c      write(*,*)'tini,tm,htini,epsy'
c     1 ,tini,tm,htini,epsy
      t1=tini
      rt1=rbmin+t1*costet0
      zt1=zbmin+t1*sintet0
      bmod=b(zt1,rt1,phi)
      rn1=-1.d0/y(zt1,rt1,phi,iw_j)
 10   t2=t1+ht
c      write(*,*)'10 t2,t1,rt1,zt1,rn1',t2,t1,rt1,zt1,rn1
      t2=dmin1(t2,tm)
c      write(*,*)'t2',t2
      if((t2.eq.tm.and.t1.eq.tm).or.
     1   (t2.eq.0.d0.and.t1.eq.0.d0)) then
	 ierr=0
c	 write(*,*)'in zrcntr2 error exit ierr',ierr
c         write(*,*)'t2,t1,rt1,zt1,rt2,zt2',t2,t1,rt1,zt1,rt2,zt2
c	 write(*,*)'in zrcntr2 error rn1,rn2,n0',rn1,rn2,n0
c        error exit
         goto 30
      endif

      rt2=rbmin+t2*costet0
      zt2=zbmin+t2*sintet0
      bmod=b(zt2,rt2,phi)
      rn2=-1.d0/y(zt2,rt2,phi,iw_j)
      dn=(dfloat(n0)-rn1)*(dfloat(n0)-rn2)
c      write(*,*)'10 t2,rt2,zt2,rn2,dn',t2,rt2,zt2,rn2,dn
      if(dn.le.0.d0)go to 20
      t1=t2
      rn1=rn2
      goto 10
c-----------------------------------------------------------
c     n0 is between rn1 and rn2
c     binary iteration methode
c-----------------------------------------------------------
 20   continue
      tr=t2
      tl=t1
      do while ((tr-tl).gt.epsy)
         t=tl+(tr-tl)*0.5d0
         r=rbmin+t*costet0
         z=zbmin+t*sintet0
         bmod=b(z,r,phi)
         rn1=-1.d0/y(z,r,phi,iw_j)-n0
         rtr=rbmin+tr*costet0
         ztr=zbmin+tr*sintet0
         bmod=b(ztr,rtr,phi)
         rn2=-1.d0/y(ztr,rtr,phi,iw_j)-n0
         if ((rn1*rn2).gt.0) then
            tr=t
         else
            tl=t
         end if
      end do
c     -----------------------------------------------
c          end of the binary methode
c     -----------------------------------------------
      t0=t
      zt0=z
      rt0=r
 30   continue
c      write(*,*)'the end of zrctr2 t0',t0
      return
      end


      subroutine check_param(i_op)
c-----check parameters in param.i
c     i_op=1 check the parameters for eqdsk
c     i_op=2 check the parameters for grill in genray.in

      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i' ! gives the input parameters 
      include 'commons.i'
   
      if ((i_op.lt.1).or.(i_op.gt.2)) then
         write(*,*)'in genray.f,in check_param: wrong value i_op=',i_op
         write(*,*)'It shoud be (i_op.eq.1).or.(i_op.eq.2)'
         write(*,*)'Change i_op in call check_param(i_op)'
         stop
      endif

      if (nraymax.lt.0) then
         write(*,*)'in param.i nraymax <0,it should be >0'
         write(*,*)'Change nraymax in param.i'
         stop
      endif 

      if (i_op.eq.1) then
c--------check the parameters for eqdsk 
              
         if (nveqd.gt.nreqd) then
            write(*,10110)
            stop
         endif
10110    format("subroutine equilib-nveqd > nreqd")
cSm030224
         if ( (nreqd.gt.nreqda).or.(nzeqd.gt.nzeqda)
     +    .or.(nxeqd.gt.nxeqda).or.(nyeqd.gt.nyeqda) ) then
           write(*,1000) nreqd,nreqda,nzeqd,nzeqda
 1000 format('in equilib.dat in input',/,
     .'the dimensions of eqdsk (in eqilib.dat) nreqd or nzeqd',/,
     .'are bigger then the parameters nreqda or nzeqda in param.i'
     .,/,'nreqd=',I5,'nreqda=',I5
     .,/,'nzeqd=',I5,'nzeqda=',I5
     .,/,'Change nreqda or nzeqda or nxeqda or nyeqda in param.i')
           stop
         endif

         if (nrya.ne.(max0(nreqda,nzeqda)+4)) then
            write(*,*)'in param.i nry.ne.(max0(nreqda,nzeqda)+4)'
            write(*,*)'nrya=',nrya
            write(*,*)'max0(nreqda,nzeqda)',max0(nreqda,nzeqda)
            write(*,*)'Change nrya in param.i'
            stop
         endif 
      
         if (nzy.ne.(max0(npsi4,nteta1)+4)) then
            write(*,*)'in param.i nzy.ne.max0(npsi4,nteta1)+4'
            write(*,*)'nzy= ',nzy
            write(*,*)'max0(npsi4,nteta1)',max0(npsi4,nteta1)
            write(*,*)'Change nzy in param.i'
            stop
         endif
          
         goto 10

      endif !i_op=1

      if (i_op.eq.2) then
c--------check the parameters for grill in genray.in
         if (ngrill.gt.ngrilla) then
 20   format('equilib in check_param ngrill>ngrilla',/,
     .'it should be ngrilla.ge,ngrill',/,
     .'ngrilla= ',I4,'ngrill= ',I4,/,
     .'Change ngrilla in param.i or grilld in genray.in')
              write(*,20)ngrilla,ngrill
            stop
         endif
 
         nmax=0
         do i=1,ngrill
           if (nnkpar(i).gt.nmax) nmax=nnkpar(i)
         enddo

         if (nnkprmax.ne.nmax)then
            write(*,30)nnkprmax,nmax
 30          format('genray.f in check_param: it should be',/,
     .      'nnkprmax.ge.max{i=1,ngrill}nnkpar(i)',/,
     .      'nnkprmax= ',I4,'max{i=1,ngrill}nnkpar(i)= ',I4)
            goto 10             
         endif 
          
         nmax=0
         do i=1,ngrill
           if (nthin(i).gt.nmax) nmax=nthin(i)
         enddo

         if (nthinmax.ne.nmax)then
            write(*,40)nthinmax,nmax
 40          format('genray.f in check_param: it should be',/,
     .      'nthinmax.ge.max{i=1,ngrill}nthin(i)',/,
     .      'nthinmax= ',I4,'max{i=1,ngrill}nthin(i)= ',I4)
            goto 10             
         endif 

      endif !i_op=2

 10   continue

      return
      end
              



