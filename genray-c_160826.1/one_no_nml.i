
c This include file contains the variables which are used 
c in most of the subroutines and functions:  
c magnetic field; small radius; wave frequency;
c sin and cos of the angle gam between the refractive index and
c the magnetic field;
c the constants w and v  for the calculations Y and X;
c the characteristics of the profiles of the density, temperature, tpop,
c vflow, zeff;
c the characteristics for the toroidal density variations; 
c the values of the different switches
c No  namelist data.

   
      character t_*512        !Character variable for temporary use.
      character filetxt*128,filenc*128
        
      real*8  
     &            a,q0,pq,rho,rhowall, 
     &            drhodzb, drho_dx, drho_dy, drho_dz,
     +            vgr_x, vgr_y, vgr_z,
     &            ax,pi,gam, beta_frc, akappa, bx,by,bz,
     &            br,bphi,bmod,o_bmod, psid,
     +            dpdxd,dpdyd,dpdzd,dpdrd,dpdzzd,dpdrrd,dpdzrd,
     &            dbzdz,dbzdr,dbzdph,dbrdz,dbrdr,dbrdph,
     &            dbphdz,dbphdr,dbphdph,
     +            dbxdx,dbxdy,dbxdz, dbydx,dbydy,dbydz, dbzdx,dbzdy,
     &            dqdrho,dbmdx,dbmdy,dbmdz,dbmdr,dbmdph,
     &         	  dc2dz,dc2dx,dc2dy, dc2dr,dc2dph,
     &            dc2dnz,dc2dnx,dc2dny, dc2dnr,dc2dm,
     &            w,v,
     &            te0,teb,
     &            powtott,powini,
     &	          dpsimax,
     &            power_launched,
     +            glb_mirror
  
      integer
     &            icount,
     &            irefl,
     &            i_,i1_,iboundb,       
     &            nstep_rk,
     &            n_relt_harm2,
     &            k_root,               !the number of eigenvalues
     &            i_total_bad_initial_conditions,
     &            lrzmax

      common/one_no_nml/ a,q0,pq,rho,rhowall,
     &            drhodzb,drho_dx,drho_dy,drho_dz,
     +            vgr_x, vgr_y, vgr_z,
     &            ax,pi,gam, beta_frc, akappa, bx,by,bz,
     &            br,bphi,bmod,o_bmod, psid,
     +            dpdxd,dpdyd,dpdzd,dpdrd,dpdzzd,dpdrrd,dpdzrd,
     &            dbzdz,dbzdr,dbzdph,dbrdz,dbrdr,dbrdph,
     &            dbphdz,dbphdr,dbphdph,
     +            dbxdx,dbxdy,dbxdz, dbydx,dbydy,dbydz, dbzdx,dbzdy,
     &            dqdrho,dbmdx,dbmdy,dbmdz,dbmdr,dbmdph,
     &         	  dc2dz,dc2dx,dc2dy, dc2dr,dc2dph,
     &            dc2dnz,dc2dnx,dc2dny, dc2dnr,dc2dm,
     &            w(nbulka),v(nbulka),         
     &            te0(nbulka),teb(nbulka),
     &            powtott,powini,
     &	          dpsimax,power_launched,
     +            glb_mirror,
     &            icount,
     &            irefl,
     &            i_,i1_,iboundb,   
     &            nstep_rk,
     &            n_relt_harm2,  
     &            k_root,               !the number of eigenvalues
     &            i_total_bad_initial_conditions,
     &            filetxt,filenc,t_,
     &            lrzmax
  
c--------------------------------------------------------------------
c     a           - small radius(normalized)
c     b0          - characteristic magnetic field (tl)
c     q0,pq	  -
c     rho	  -normalized small radius
c     psifactr    =<1, is used in gr2new to find the last flux surface
c                  psi(r,z)=psilim*psifactr
c     ax	  -
c     r0x         - characteristic radius(m)
c     frqncy      - wave friquency (GHz)
c     pi	  -	  4*datan(1.d0)
c     gam         - angle between magnetic field and refractive index(rad)
c     bz,br,bphi  - components of the magnetic field
c     bmod        - module of the magnetic field
c     psid        - poloidal flux
c     dpdxd,dpdyd,dpdzd,dpdrd - derivatives of psi 
c     dbzdz,dbzdr,dbzdph -derivatives from z-component of magnetic field
c     dbrdz,dbrdr,dbrdph -derivatives from r-component of magnetic field
c     dbphdz,dbphdr,dbphdph- derivatives from toroidal magnetic field
c     dbxdx,dbxdy,dbxdz, dbydx,dbydy,dbydz, dbzdx,dbzdy - derivatives of b-components
c     dqdrho
c     dbmdx,dbmdy,dbmdz,dbmdr,dbmdph -derivatives from the module of magnetic field
c     dc2dz,dc2dr,dc2dph -derivatives from cos(gam)**2 by z,r and phi
c     dc2dnz,dc2dnr,dc2dm -derivatives from cos(gam)**2 coordinates of the
c                         refractive index n_z,n_r and m
c     dc2dz,dc2dx,dc2dy = d[cos^2(gamma)]/dz, d[cos^2(gamma)]/dx, d[cos^2(gamma)]/dy
c     dc2dnz,dc2dnx,dc2dny = d[cos^2(gamma)]/dnz, d[cos^2(gamma)]/dnx, d[cos^2(gamma)]/dny
c            v0=806.2/(frqncy*frqncy)     ! set in dinit.f
c            v(i)=v0*charge(i)**2/dmas(i) ! set in dinit.f
c            w0=28.0*b0/frqncy            ! set in dinit.f;  b0=1.d0 (normalization)
c            w(i)=w0*charge(i)/dmas(i)    ! set in dinit.f
c     dense0,denseb - central and bound electron densities (in 10**13/cm**3)
c     den(rho)=denseb+(dense0-denseb)*(1-rho**rn1de)**rn2de
c     YuP[Added,2014]: Now, the density profile is set as den(rho)*gn(z),
c                      where gn(z) is the dropoff function,
c                      such that at |z|<zbegin_den_dropoff , 
c                      the density profile is just den(rho) (with gn(z)==1),
c                      but at |z|>zbegin_den_dropoff ,
c                      the density drops exponentially over the 
c                      effective length = zlength_den_dropoff .
c                      By default, zlength_den_dropoff is set to 0.d0, 
c                      which will set gn(z) to 1.0 everywhere.
c                      This option is only available for model_rho_dens=0, for now.
c     te0,teb - central and bound electron temperatures (keV)
c     tempe(rho)=teb+(te0-teb)*(1-rho**rn1te)**rn2te !average temperatura
c     tpop(rho)=tpb+(tp0-tpb)*(1-rho**rn1tp)**rn2tp  !T_perp/T_parallel
c     vfow(rho)=vflb+(vfl0-vflb)*(1-rho**rn1vfl)**rn2vfl  !drift velocity || B
c     zeff0,zeffb - central and bound effective charge
c     zeff0(rho)=zeffb+(zeff0-zeffb)*(1-rho**rn1zeff)**rn2zeff
c     powtott  - total power summed over antennas
c     powini   - initial power in ray_iray chanal
c     delpwrmn - Minimum power in each ray, as a fraction of
c                starting power in the ray, after which ray is
c                stopped.  
c     nbulk -number of the plasma components
c     ioxm - sign before square root in dispersion relation
c           N=N(theta)
c           root=(-b+ioxm*sqrt(b**2-4a*c))/2a
c     ioxm_n_npar - sign before square root in dispersion relation
c           N=N(N_parallel)
c           root=(-g+ioxm_n_npar*sqrt(g**2-4f*w))/2g
c     icount -
c     ib - index of plasma component for which delib=1-yib may be
c                equal zero inside the plasma,
c                dispersion relation is multiplied by delib
c     id - index of the form of the dispersion relation
c     isolv -index of the solution method of ray trasing equations
c     idif  -index of analytical(=1) or numerical (=2) differentiation
c                  of Hamiltonian
c     irkmeth- index of the modification of Runge-Kutta  procedure
c     ireflm - maximum number of ray reflections
c     irefl  - number of ray reflections  
c                with reflection absorption
c     istart - index to shoose wave start is outside the plasma (EC)
c                    or inside plasma (LH,FW)
c     idens  - index of spline or quasiparobolic density,temperature
c                     and z_eff  approximation
c     indexrho -index to choose radial coordinate
c     ionetwo - index to calculate power and current drive profiles
c     jwave   -	number of the wave harmonic
c     i_      - number of the file for wr! i_look_roots=0    !do not plot D(N_perp) and do not calculate all 
!                   !hot roots
!             =1    !plot D(N_perp) and do calculate all 
!                   !hot roots and do not calculate ray   
!             =2    !calculate hot root, use the root with number k_hot_root
!                   !as the initial ray condition and calculate rayiting 3D data
c     i1_     - number of the file for writing ray data
c     izeff   =0,then zeff will be calculated from the given plasma
c              components; =1 then zeff are given, but the ions
c              components(density) will be calculated
c		iboundb  this flag controls (in subroutine b.for)
c              is the point inside the plasma or not:
c              =-1 the ray point is inside th   
c	       =2 r<rmin+epsbnd or r>rmax-epsbnd  (outside the plasma)
c	       =1 z<zlimiter_min(r) or z>zlimiter_plus(r)(outside the plasma)
c		drhodzb  is a derivative ! i_look_roots=0    !do not plot D(N_perp) and do not calculate all 
!                   !hot roots
!             =1    !plot D(N_perp) and do calculate all 
!                   !hot roots and do not calculate ray   
!             =2    !calculate hot root, use the root with number k_hot_root
!                   !as the initial ray condition and calculate ray of the small radius from z.
c              (For the case when the point is outside the plasma).
c		It was calculated in b.for.
c     iabsorp -choice of Imag(N_perp) 
c             (N_perp is the perpendicular refractive index)
c              iabsorp=1 for EC waves from Mazzucato solver
c               =2 for LH waves, =3 for FW waves
c For use with Abhay Ram's dispersion relation (id=14) parameters to
c control the integration routipne for the Trubnikov integral.
c
c relres offers 4 choices for the resolution to use for the relativistic
c dispersion relation using the Trubnikov integral:
c relres=1 low resolution, errabs0=1.d-4, errrel0=1.d-4, navg=3, diff_err=0.1
c       =2 medium res., errabs0=1.d-5, errrel0=1.d-5, navg=12, diff_err=1.d-3
c       =3 high res., errabs0=1.d-6, errrel0=1.d-6, navg=25, diff_err=1.d-6
c       =4 user-defined res., set the following parameters manually:
c The Trubnikov one-dimensional (complex) integral is performed by splitting
c up a region from 0 to 1.d8 into 10^6 pieces, and each piece is integrated
c using the SLATEC adaptive quadrature routine dqag. errabs0 and errrel0 are
c the absolute and relative error tolerances passed directly to dqaq.
c Then the adjacent pieces are compared (it is an oscillatory integrand)
c and using navg number of pieces, when the average difference between them
c are less then diff_err, the integration is presumed finished (Thus it may
c finish long before the upper limit of 1.d8).
c
c errabs0 - absolute error for dqag integration routine
c errrel0 - relative error for dqag integration routine
c navg - number of adjacent integration intervals to use in comparison
c diff_err - error tolerance using navg pieces, when the average difference
c         is less than diff_err, then the integration is done.
c--------------------------------------------------------
c       flux=B~.B+E~.d(omega*eps_herm)/(domega).E
c    iflux=1 the flux will be calculculated using the the group velocity from
c           the choosed disperson relation (with given id) and the electric
c           field calculated for choosed iabsorp                  
c    iflux=2 the flux will be calculated using V_gr for the electron cold
c            plasma dispersion and polarization (using subroutine  grpde2)     
c-------------------------------------------------------
c     ipsi    -=1 calculation of conturs psi(z,r)=const,
c              =0 reading these contours from psi.bin file 
c                 (see genray.in file)
c     dpsimax -if in eqdsk (psimag.gt.psilim) then dpsimax=-1. else dpsimax=1.
c               It was determined in equilib.for: subroutine input
c     nrelt     Maximum number of ray elements per ray.
c               Must be .le. nrelta (a parameter)
c     iwcntr   =1 genray will call mk_grapc(iray) for the  
c                 calculations of contours
c               wb_c=i,
c              =0 call mk_grapc will be commented
c iwcntr =1 genray.f will calculate the contours wb_c=n
c        =0 genray.f will not do it
c iwopen =1 mk_grapc will calculte open contours wb_c=n	(using contrb1)
c         2 mk_grapc will calculte close contours wb_c=n (using contrb2)
c iwj    =mk_grapc will calculate contours wb_cj=n, here j is a kind of plasma
c         component must be.le.nbulk, j=1 for the electron hirofrequency
c         j=2,  for the ion (j kind) hirofrequency
c ieffic  chooses the formula for the current drive efficiency
c        =1 asymptotic symple formula (homogenius, nonrelativictic)
c        =2 asymptotic symple formula (East-Karney )
c        =3 asymptotic symple formula (curba subroutine)
c ibw    =0 the launch is not for Bernstein waves
c        =1 the launch of electron Bernstein wave from dhot tensor
c           The last case works only for istart=2 and grill_lh conditions 
c------------------------------------------------------------------------
c The change of the dispersion relation and absorption
c near the gyro-frequency points
c-------------------------------------------------
c iswitch=1   To use the change of the dispersion relation and
c             absorption
c        =0   Do not use the change of the dispersion relation and
c             absorption 
c     del_y   If the difference |1-nY(jy)|<del_y 
c             (jy=1-nbulk ,n=...-2,-1,0, 1,2,3,...)
c             then switch on the new 
c             given type of the disperion and absorption.
c   jy_d      is the type of plasma species 1<=jy_d<=nbulk
c   idswitch  is the type of the dispersion function near the 
c             gyro-frequency points
c             It can be equal 1,2,3,4,5,6
c   iabswitch is the type of the absorption near the gyro-frequency point  
c-------------------------------------------------------------------
c   The parameters for the density variation     
c   var0,denn,denm,an,sigman (see x.f: vardens)
c   var0 is an amplitude of the density variation (del_n_0) (see 3.37...)
c   denm is the number of the poloidal mode in the density variation(l_theta)
c   denn is the number of the toroidal mode in the density variation(l_phi)
c   an   is the radial localization of the variation (rho_0)
c   sigman is the parameter that characterizes the radial thickness
c          of the density fluctuation    
c----------------------------------------------------------------
c   the switch of the correction in outpt
c   icorrect=0 switch off the correction
c           =1 switch on the correction
c----------------------------------------------------------------
c  itools =0 do not use mkgrtool in genray.f
c         =1 use mkgrtool in genray.f
c---------------------------------------------------------------
c  nstep_rk    is the number of the Runge-Kutta time step 
c
c  maxsteps_r  is the maximal number of the time steps of Runge-Kutta
c              solver
c----------------------------------------------------------------
c iout3d    ='enable'  to write output 3d.dat wile   [DEFUNCT]
c           ='disable' not write 3d.dat file         [DEFUNCT]
c--------------------------------------------------
c i_vgr_ini =+1 the wave is directed into the plasma (in the initial point)
c           =-1 the wave is directed out the plasma (in the initial point) 
c--------------------------------------------------
c outdat is the name of the output file that contains the data along the ray
c        r,z,phi,n_r,n_z,m (usually it outdat='zrn.dat')
c stat   is the status of outdat file (usually stat='new')        
c netcdfnm !name of the input *.nc file with 3d non-maxwellian distribution
c------------------------------------------------------------------------
c The data for the creation of the analytical
c non-Maxwellian electron distribution:
c     rtem0,
c     rtail,r1t,r2t,ttail,rhot,r1h,r2h,thotpar,thotper,
c     hotexp,hotmnpar,hotmxpar,hotmnper,hotmxper,
c     rbeam,r1b,r2b,tbeam,ebeam,thbeam
c     jx,iym,lrz,ngen,lrzmax 
c The data fornalytic calculation of continuous non-Maxwellian distribution
c with  three  temperaures in three energy  ranges
c     rvtovte1,rvtovte2,rtemp1,rtemp2,rtemp3
c--------------------------------------------------
c  n_relt_harm1 is the lowest, i.e., minimum harmonic used in the
c               anti-hermitian dielectric tensor calculations.
c               It man be positive of negative.
c               Default value is +9999, in which case this input
c               is ignored.
c   n_relt_harm (.ge.1) gives the number of EC harmonics used 
c               in anti-hermitian dielectric tensor calculations
c               If n_relt_harm1=9999, then harmonics from 
c                 -n_relt_harm to +n_relt_harm are used.
c               If n_relt_harm1.ne.9999, then harmonics from
c                  n_relt_harm1 to n_relt_harm1+n_relt_harm are used.
c     It is necessary that the harmonics used in this calculation
c     be within the range of parameters [n_relt_harm1a,n_relt_harm2a]
c-----------------------------------------------------------
c     the data for the emission calculations
c     i_emission=0 no emission
c               =1 emission calculations
c---------------------------------------------------------
c   0<tol_emis=<1 tolerance parameter to add the new mesh point s_n
c                if in_0(n)>tol_emis*i_0 
c--------------------------------------------------------
c   nfreq=<nfreqa (see param.i) is the number of frequences
c         nfreq=1 gives detailed plots of emission from a single ray
c         ngraq.gt.1 gives spectra covering the specified frequency
c          range
c   freq00=<freq01 are ratios of the minimal and maximal emission
c         frequencies to the central electron gyro-frequency
c         f_ce(rho=0) 
c------------------------------------------------------------- 
c   nharm1=<nharm=<nharm2 give the number of used EC emission haremonics 
c
c   wallr is the wall reflection coefficient
c
c     i_rrind chooses the subroutine to calculate N_ray 
c     i_rrind=0  N_ray=N
c     i_rrind =1 from cold electron plasma using rrind   
c     i-rrind =2 from hot non-relativistic dispersion relation 
c---------------------------------------------------------------
c     i_r_2nd_harm=1 to calculate the major radius of the EC 2nd harmonic
c                 =0 do not calculate (default =0) 
c---------------------------------------------------------------
c poldist_mx is the maximal polidal distance (cm) along the ray
c            default=1.d+5
c---------------------------------------------------
c   temp_scale(ncomp),den_scale(ncomp) are the parameters to multiply
c   the given temperature and density profiles
c--------------------------------------------------------------------
c i_plot_b =1 create figures for the magnetic field,density and temperature 
c             profiles in plot.ps file using subroutine map_b based on PGplot
c i_plot_b =0 do not write the b,n,T figures to plot.ps file 
c
c i_plot_d =1 create the dispersion function contours d(ReNperp,ImN_perp)
c             in plot.ps file using PGplot
c i_plot_d =0 do not create the dispersion function contours
c-----------------------------------------------------------
c i_im_nperp chois of the metode to find Im_N_perp for hot plasma (iabsorp=4)
c-----------------------------------------------------------------------
c i_im_nperp=1 Im_N_perp=abs(ImD_full/(dD_hermitian/dReN_perp)) 
c i_im_nperp=2 (Re_N_perp,Im_N_perp) is the complex root 
c              (of the complex dispersion relation)
c              calculated by Newton iterations with the numerical
c              derivatives (the chord method)
c---------------------------------------------------------------
c  ndens is the number of point for the radial profiles 
c        density and temperature
c---------------------------------------------------------------
c k_root,       the number of the dispersion tensor
c               eigenvalue, used for id=10, determined in cninit.f
c----------------------------------------------------------- 
! i_geom_optic sets  the form of the ray equations
!              =1  integration in time (default):
!                  ray-tracing equations right hand side=
!                  dr^/dt=-(dD/dN^)/(dD/domega)
!                  dN^/dt=+(dD/dr^)/(dD/domega)
!                  In this case rside1 gives v_group.
!              =2  integration is space,
!                  ray-tracing equations right hand side=
!                  dr^/dl=- ray_direction * (dD/dN^)p
!                  dN^/dl=  ray_direction * (dD/dr^)p
!                  p=1.d0/dsqrt(deru(1)**2+deru(2)**2+deru(3)**2)=
!                  deru(1)=dD/dNx,deru(2)=dD/dNy,deru(3)=dD/dNz,
c----------------------------------------------------------------------
c ray_direction =+1 as default it
c               -1
c It is a multiplier in right hand side of ray-tracing equations
c It is used for i_geom_optic=2 case
c----------------------------------------------------------------------
c i_ox =0 /default/ do not use the vertex coordinates calculations 
c      =1 use EC cone vertex cordinates calculations
c      =-2 creates ray jump at OX conversion  area
c-----------------------------------------------------------------------------
c theta_bot(icone)<theta_top(icone) 
c are the poloidal angle boundaries (degree) at
c Xe=1 surface.They are used in genray.f to find the optimal ray
c direction at the given EC cone vertex icone=1,...,ncone
!----------------------------------------------------------
! eps_antenna 
! is the vertical (in Z) extension of antenna;  needed for i_ox=1 case.
! In i_ox=1 case, several rays are traced from O-cutoff layer 
! back to plasma edge, to the surface that contains antenna
! (the shape of this surface is determined by i_ant, 
! see the beginning of subroutine outpt_xyz). 
! If the end point (r_st_ox,z_st_ox)
! is close enough to the antenna position, the ray is selected as
! the "optimal" ray to be launched from antenna in the run with i_ox=2.
! More specifically, for the antenna positioned  
! at rst(icone),zst(icone)  (specified in genray.in),  
! the ray is selected if(abs(z_st_ox-zst(icone)) .le. eps_antenna*0.5)
!---------------------------------------------------------
c i_ox_poloidal_max is the maximal number of the used poloidal angles
c                   in the bisection method. This method calculates the
c                   the poloidal angle theta_pol of the point M at the
c                   flux surface Xe(rho_=1 : M(poloidal_angle,rho).  
c                   The ray lauched from M to the plasma edge will
c                   go to the EC cone vertex.
c---------------------------------------------------------
c eps_xe   Is the parameter which sets the vicinitity
c          of the O-mode cutoff surface
c          If xe > (1-eps_xe) then the subroutine ox_conversion  
c          makes the ray jump in small radius direction
c          and finds X mode.
c-------------------------------------------------------------------            
c partner  = 'disabled' to use input profiles from genray.dat/genray.in
c          = 'onetwo'   to use input plasma profile data from the file: 
c                       genray_profs_in created by transport code
c                       and to write the output file of power absorption
c                       and current drive genray_profs_out
c                       for the transport code.
c-------------------------------------------------------------------
c i_output =1 output at equal poloidal distance
c          =2 output at equal total distance
c------------------------------------------------------------------
c prmt6 distance step [m] for output of ray data
c------------------------------------------------------------------
c i_uh_switch=1    if uh=dsrt(xe+ye**2) < uh_switch then change 
c                  the output step prmt(6) for prmt6_uh_switch 
c            =0    do not change the output step prmt(6)
c prmt6_uh_switch  is the output step for i_uh_switch=1 case
c uh_switch        if uh<uh_switch then change the output step for    
c                  i_uh_switch=1 case
c---------------------------------------------------------------------
c      iabsorp_collisional =0 no additional collisional absorption
c                          =1 collisional absorption  using formula
c                             Im(N)=dabs(nu_ei/(gr_perp))*clight/omega)
c      coll_mult =1.d0(default), multiplies above coll absorp expression
c---------------------------------------------------------------------------
c i_salphal (nbulka)  sets which damping will be in salphal_nc
c	     Default:i_salphal(1)=1,i_salphal(2:nbulk)=0 electron damping only  
c            A particular species contribution to salphal_nc is added if
c            i_salphal(species_number)=1. That is, damping coefficients
c            for all species with i_salphal(k).ne.0 are summed into saplhal_nc 
c--------------------------------------------------------------------------- 
c i_plot_disp_cold    ! =1 to plot D_cold(N_perp) to plot.ps file 
c                     ! =0 do not plot  
c                     ! It is used only  in grill_lh to create
c                     ! the plot in the initial point
c-----------------------------------------------------------------------
c i_plot_wave_normal_cold! for plotting cold wavenormal to plot.ps file
c                        !  =1 to plot wavenormal
c                        !  =0 do not plot
c------------------------------------------------------------------------
c nonuniform_profile_mesh= 'enabled' use nonuniform small radius mesh for input
c                                spline profiles (works for idens=1 only)
c                = 'disabled'    do not use nonuniform mesh
c--------------------------------------------------------------------------------  
c  Calculation of the small radius value near the plasma edge
c  where LH or FW have cutoff:  
c  i_rho_cutoff=0 (default) no calculations
c              =1 use these calculations
c--------------------------------------------------------------------
c  rho_step_find_LHFW_cutoff  is the non dimensional small radius step
c                            used in subroutine  rho_ini_LHFW                            
c                            It is used at i_rho_cutoff=1
c---------------------------------------------------------------------
c  rho_initial_find_LHFW_cutoff  is the initial small radius point. 
c                            As default rho=1- rho_step_find_LHFW_cutoff
c                            It is used at i_rho_cutoff=1
c--------------------------------------------------------------------------------
c           i_emission_spectrum      !to calculate emisssion spectrum   
c-----------------------------------------------------------------------------
c           cnperp_plot_min,cnperp_plot_max !max and min Nperp to plt D(Nperp)      
c           n_nperp_plot,                   !number of Nperp points to plot D(Nperp)
c-----------------------------------------------------------------------------------
c            cN_perp_root_max               !max value of n_perp to
c                                           !seek hot roots 
c            n_points_root                  !number of  N_perp mesh points
c                                           !to find hot plasma roots
c-------------------------------------------------------------------------------------
c i_look_roots=0    !do not plot D(N_perp) and do not calculate all 
c                   !hot roots
c             =1    !plot D(N_perp) and do calculate all 
c                   !hot roots and do not calculate ray   
c             =2    !calculate hot root, use the root with number k_hot_root
c                   !as the initial ray condition and calculate ray
c k_hot_root        is the number of the hot plasma root
c                   N_perp_root_ar(k_hot_root)
c                   which will be used for ray initial condition
c                   It works for i_look_roots=2 case only
c-------------------------------------------------------------------------
c For subroutine  rho_ini_hot_nperp_roots to find hot plasma roots
c i_rho_find_hot_nperp_roots,
c rho_step_find_hot_nperp_roots,
c rho_min_find_hot_nperp_roots,
c--------------------------------------------------------------
c  measure error in the dispersion realation
c  If D/(N|gradD|) > toll_hamilt stop ray calculation
c
c  toll_hamilt 
c--------------------------------------------------------------
c total number of rays of rays with bad initial conditions, having 
c i_bad initial conditions=1:
c
c i_total_bad_initial_conditions
c---------------------------------------------------------------
c power_launched ! total power launched along rays with good initial
c                ! conditions
c----------------------------------------------------------------
c thetapol_wall(n_wall_a),      ! poloidal angles [radians] of
c thetapol_limiter(n_limiter_a) ! wall and limiter points     
c------------------------------------------------------------------
