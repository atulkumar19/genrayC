
c     a_change.i
c
c
c***********************************************************************
c
c     This file records changes in the code, starting from 2014. 
c
c     Records start from the bottom of the file.
c***********************************************************************


c[35] YuP[11-2016] Added capability to read eqdsk in case of a mirror
machine. Use model_b=0, specify eqdskin=(your eqdsk file name),
and specify eqdsktype='mirror'.
The bin volume is properly defined now, based on pol.flux levels,
open field lines. See subr.gr2new, subr.rhospl and subr.onetwoini.
Look for (model_b.eq.0 .and. eqdsktype.eq.'mirror').



c[34] YuP[11-2016] Adjustment in subr. absorpfd_xyz().
The limits L1,L2 of summations over Bessel functions 
can be larger than nb (set to 500 in the subroutine),
when ray approaches B~0 (it can happen in FRC).
Although the value of rho_larm is limited (see rho_larm_max),
the value of Nperp (Kperp) can be very large in case of HHFW
with Nperp~c/Valfven  ~ c/B.
Now the ray is allowed to continue, with abs(L1) and L2 reset to be
smaller than nb.



c[33] YuP[11-2016] Correction in subr. dwpw_2() which calculates
(numerically) the derivatives of (omega_pl_i/omega)**2.
Now the step of differentiation is taken from namelist:
hh=der_r 
Before [11-2016] it was "hard-wired": hh=rmax*1.d-4


c[32] YuP[11-2016] Modified subr. ZFUN_vkw() - now it has two arguments
related to value of Npar: vkw_in and vkw_adj 
(vkw is Npar*beta== Vthermal*Kpar/omega ).
One value (vkw_adj) is based on adjusted Npar, 
with lower-limit value set for |nll|==|Npar|. 
This is done to avoid a jump in plasma dispersion function 
(Real part of CZ0) when Npar~0.
However, for calculation of damping (Im part of CZ0), we use 
the original nll_in, i.e. using vkw_in.
The adjustment of nll_in is done in subr. DHOT_s().


c[31] YuP[11-2016] Added a new namelist input, to be set in genray.in,  
  in &ox section :
  ox_step_dir== type of stepping in O-X mode coupling runs,
       when looking for X-mode across evanescent layer.
       'gradne' - step along grad(n) direction, or
       'vgroup' - along Vgroup velocity direction which is
                  taken just before O-X jump is initiated.



c [30] YuP[11-2016]  Added a new namelist input, to be set in genray.in,  
  in &dispers section :
  rho_larm_max== Upper limit for Larmor radius [cm]
  It is used for absorption calculation 
  (for iabsorp=3, i.e. when using subr. absorpfd_xyz, for now).
  This is done to address the inflation of argument of bessel function
  when ray goes over B~0 point (in FRC).
  Recommended: In FRC run (when magn. field B can go through zero)
  set it to the distance between null point r0 and separatrix rs. 
  Default value is rho_larm_max=10.d10, which effectively means 
  no upper limit.
  Physically, even at B=0 the Larmor radius cannot be zero because 
  a particle travels away from the localized B=0 layer
  and spends most of time in regions where B>0.
  An accurate value of rho_larm_max can be found from particle tracing
  (using some other codes), although the trajectory may not look like
  a circular orbit, especially for fast particles.
  If the value of rho_larm_max is set to a very small value, the damping
  of waves on fast ions (finite-Larmor radius effects) will be reduced 
  because the value of Kperp*rho_larm will be too small.
  Setting rho_larm_max to the recommended distance (rs-r0) improves 
  the stability of runs.
  
             


c----------------------------------------
c version="Genray-c_160826"
c----------------------------------------

c [29] Many adjustments in plotting script genray_c_plot.py.

c [28] Fixed a bug in grill_lh_xyz, related to power setting 
  for igrillpw.eq.1 - the powers for each ray should be same,
  but they were different.   YuP[08-2016]

c [27] YuP[07-2016] Added values of psimag and psilim into *.nc file
c and adjusted genray_c_plot.py accordingly:
c now the levels of PSI(R,Z)=eqdsk_psi are only plotted up to psilim.

c----------------------------------------
c version="Genray-c_160718"
c----------------------------------------

c [26] YuP[07-2016] Fixed a couple of bugs in netcdfr3d
c related to defining 'powtot_s' netcdf name 
c (rdims was not properly defined),
c and definition of 'eqdskin' name for the nc file
c (char512id was not defined).
c Also added check_err(istatus) after each 
c vid=ncvdef2() call. It will help to detect an error
c if there any problem with definition of netcdf names.


c----------------------------------------
c version="Genray-c_160405"
c----------------------------------------

c [25]  YuP[03-2016,04-2016] Added new option model_b=3. 
 In genray.in (namelist &tokamak), need to specify:
 the value of rmirror, rbox_mirror, zbox_mirror, b00_mirror. 
 See GENRAY-c_help for a detailed description.
 The model is for a mirror machine; the same model is present in CQL3D
 (named "mirror1" model) and can be used for a coupled run GENRAY-CQL3D 
 for a mirror machine: set same values for the model, 
 same T(rho) and n(rho) profiles; run Genray-C first, copy the produced 
 genray.nc file with rays' data, run CQL3D using urfmod='enabled' option
 that allows reading the genray.nc file.
 


c----------------------------------------
c [2016-03-21] version="Genray-c_160321"
c----------------------------------------

c [24] YuP [03-2016] Added new namelist variable: rho_reflect.
It is used for checking that rho>rho_reflect. 
If such condition is met, then ray is reflected back towards plasma.
The default value is rho_reflect=1.d10 .
Any big number means that effectively this type of reflection 
will NOT be triggered, as the condition rho>rho_reflect is never met.
Why needed: in a magnetic mirror machine, 
the definition of rho=1 is quite arbitrary;
there is no "last closed surface".
So, it is useful to define rho_reflect in such a way
that it corresponds to the plasma edge;    
it could be quite different from 1.0.
It is recommended to set rho_reflect to be in a region where
the gradient of density is not zero, in fact it should be a "normal"
negative value (normal in a sense that density goes down 
with higher values of rho). This is done because not only 
the condition rho>rho_reflect is checked, but also the direction 
of rays' propagation (the group velocity) relative to 
the grad(ne) direction. Thus, when the rays are started outside 
of rho_reflect, the reflection is not triggered as the rays cross 
the boundary rho=rho_reflect from outside,
but it will be triggered as they cross such boundary 
from inside the plasma.
The desired value of rho_reflect can be deduced from looking at plot
like "genray_wcw_vs_R.png", which shows the graph of rho coordinate
vs x coordinate (R coordinate), usually at (y=0,z=0) plane.
This variable can also be used for a general magnetic equilibrium,
not only magnetic mirror.


c [23] YuP[03-2016] NEW option for ray integration: irkmeth=3 (used with isolv=1)
!     For irkmeth=3 option (usage of drkgs_auto), set (Example):
!      dL_step=1.d-3 ! [m]  max allowed change in cartesian coords.
!      dN_step=1.d-2 ! max allowed change in refraction index.
c     The code will set the time step h = dt_code for integration
c     in such a way that the change |dr| in configuration space
c     is not larger than dL_step, and also
c     the change in refr. index |N| is not larger than dN_step.
!     Also needed (example):
!      prmt6=1.d-3 [m] ! distance along ray [m] for saving data.
!     With the irkmeth=3 option the value of prmt6 
!     does NOT affect the step of Ruge-Kutta integration,
!     while with irkmeth=2 it does.
!     Other prmt* values are not needed.
!     Value of toll_hamilt is still operational, as usually.

c [22] YuP[03-2016] Added more options for numerical methods control,
!     through new variables in namelist: 
!     Set the steps for numerical derivatives of dispersion equation
!     over cartesian coordinates (der_r) and over refraction index (der_n).
!     This is only needed for idif=2 option.
!     Derivatives are needed for the right-hand side of ODE,
!     dD/dx, etc., dD/dn_x, etc.
!     where D(x,y,z,n_x,n_y,n_z,omega)=0 is the dispersion function,
!     or the "Hamiltonian" defined through function hamilt_xyz().
!     The derivatives are calculated in rside_xyz() and dddrz1_xyz().
!     For example, the derivative of D over x will be found as
!     [ D(x+der_r,...) - D(x-der_r,...) ] / (2*der_r)
!     Similarly - for y and z directions (using same step der_r).
!     And for the refraction index, the derivative is 
!     [ D(...,n_x+der_n,...) - D(...,n_x-der_n,...) ] / (2*der_n)
!     Similarly - for n_y and n_z components (using same step der_n).
!      der_r= 1.d-4     ! [m]
!      der_n= 1.d-4     ! [no units for refraction index] 
!     Note: If used together with irkmeth=3 option (new in 2016),
!     it is recommended to "coordinate" the values of der_r and der_n
!     with values of dL_step and dN_step;  For example, 
!     set der_r= 0.1*dL_step and der_n= 0.1*dN_step.
!     Generally, der_r should be smaller than dL_step
!     and der_n should be smaller than dN_step.
!     Smaller values are supposed to yield better accuracy,
!     however, if der_r or der_n are too small, 
!     the derivative may become zero because of computer accuracy,
!     i.e. it may happen that, within rounding error,
!     D(x+der_r,...) - D(x-der_r,...) = 0
!     although physically there should be a difference.
!     So, the steps for derivatives should not be too large
!     but also not too small.
!     The optimal values can only be found by trials -
!     they depend on scale of gradients (B and density),
!     also on wave type, proximity of resonance, etc.
!     Another derivative step:
!      der_f= 1.d-4 ! Units: fraction of frqncy f.
!     This is needed for calc. of dD/domega derivative,
!     [ D(...,f*(1+der_f)) - D(...,f*(1-der_f)) ] / (2*f*der_f)
!     The result is not very sensitive to the value of der_f.
!     It is recommended to keep it within 1.d-6...1.d-4 range. 



c----------------------------------------
c [2016-02-18] version="Genray-c_160218"
c----------------------------------------


c [21] YuP [Feb-2016] New plots are added in genray_c_plot.py:
Profile of B(z) (along the axis x=0,y=0);
Profiles of Alfven speed V_Alfven, and associated Nperp=c/V_Alfven,
and Kperp*rho_Ti based on Nperp=c/V_Alfven and thermal ion speed
 (this plot is only made if ions are present as species);
Plots of Eplus, Eminus, Ki, Nperp as a func. of distance along ray,
 one plot under another; 
Plots of resonance velocity V||res/V_Ti = (1-nres*wci_w)/(NII*Vti/clight)
as a func. of distance along ray. This is meant for ion cyclotron range,
and the resonance number nres can be specified in genray_c_plot.py,
for two possible resonances (two plots are made one under another).


c [20] YuP [Feb-2016] Some bugs are fixed in starting rays when more than 
one species is used (see subroutine cninit12_n_gam_xyz).

c [19] YuP [Jan-2016] These values for a density drop in the Z-direction
are set for each species now 
(in earlier version - for electrons only).
zbegin_den_dropoff(i)
zlength_den_dropoff(i)   
Default values: 0 => no drop-off is applied.
c     YuP[2014]: Now, the density profile is set as den(rho)*gn(z),
c     where gn(z) is the dropoff function,
c     such that at |z|<zbegin_den_dropoff , 
c     the density profile is just den(rho) (with gn(z)==1),
c     but at |z|>zbegin_den_dropoff ,
c     the density drops exponentially over the 
c     effective length = zlength_den_dropoff .
c     By default, zlength_den_dropoff is set to 0.d0, 
c     which will set gn(z) to 1.0 everywhere.
c     This option is only available for model_rho_dens=0, for now.


c [18] YuP [Feb-2016] Adjustment in subroutine cnxyz(..),
at the very end where a "Reversing Nrho (N_normal_to_surface)"
is done.
Sometimes Vgroup is nearly perp to grad(PSI),
i.e. ray travels along magnetic field lines 
or along poloidal flux = const surface.
In this case do not make a reversal of Nrho.
                
         
c [17] YuP [01-2016] Bug fix in function rhopsi(psi) and drhopsi(psi).
Sometimes, at r=0, spline gives value of
psi not exactly equal to psimag.
It can be somewhat smaller than psimag.
Example: psimag=0, and psi=-0.2d-22.
It results in negative value under dsqrt.
Added abs(), as in
 rhopsi=dsqrt(abs((psi-psimag)/(psilim-psimag))) 



c [16] Fixed a bug related to a false power absorption at OX conversion.
c It is caused by a jump of ray power at O-X conversion layer
c (if the transmission coefficient is less than 1.0)
c and by application of a "usual" adjustment of deposited power.
c This adjustment procedure should be only 
c valid at points of ray that had no O-X jump. 
c For details, see subr. p_c_prof_xyz, look for these comment lines:
c  ! The previous time step (is-1) was the last point of O-mode
c  ! and this step is the first point of X-mode.
c  ! There is generally a big change in delpwr from is-1 to is
c  ! because of transm_ox coefficient.
c  ! Do not include this jump into calculation of power to e,i,cl
c 
c It is recommended to rerun all cases where a power absorption
c on e,i,cl was observed in the vicinity of O-X transition layer.


c [15] Added model_rho_dens=5 option.
c It only works together with model_b=0 (reading eqdsk file for B field)
c and only for eqdsktype='TAE' (special form of eqdsk file
c that contains profiles of electron density).
c The profile of temperature is still setup through analytic shapes
c (using idens=0). 

c [14] Fixed bugs related to NaNs in power numbers, 
c when a ray gets to the edge of grids.
c (See subr. p_c_prof_xyz, 
c case of rhobegin and rhoend 
c being at the last point of rho_bin()
c or outside of rho_bin grid ).

c [13] Changed definition of rho_bin(). 
c Previously, the largest rho_bin(NR) was defined through rmax - 
c major radius related to chamber wall.
c With such definition, the value of rho_bin was not defined 
c at points like (R=rmax, any_large_Z).
c The rho_bin (and rho_bin_center) array is used 
c to calculate the power absorption profiles.
c Now the scan is done over (Y,Z) grid (keeping X=0) 
c and the largest rho is determined over the whole (Y,Z) grid,
c which is normally at the corners of (Y,Z) grid.
c This is done for a far-off-midplane ray launch (large Z)
c where rho_bin was not defined in the previous setup.
c See subroutine onetwoini for details.
c The exclusion from this new procedure is the case of model_b=2:
c rho_bin is setup using the old procedure, because the pol.flux 
c can go to nearly INF near the position of magnetic coils
c and so the value of rho (based on pol.flux) can be very large.
c Any absorption that happens at rho>rho_bin(NR) is attributed
c to rho_bin(NR) (lagest rho in rho_bin grid).
c The value of NRA (which is used as a default value for NR)
c is increased from 201 to 401 in param.i. 
c This is done to keep the grid size hrho=rho_bin(2)-rho_bin(1)
c about the same as before, when the range of rho_bin grid was smaller.
c It is recommended to explicitly define the value of NR in genray.in
c in section &tokamak.

c [12] Modified the procedure for "jumping" across OX cutoff,
c to find the X-mode.
c Now it is done in direction of group velocity
c just before the jump. 
c The transition from O to X-mode looks smooth now.
c The original version was - along gradient of density.
c Sometimes, when O-mode ray was approaching the cutoff
c at nearly tangential angle to magnetic line (surface),
c and the logic for the O-X jump was triggered,
c the ray was forced to step-in across the cutoff layer
c in nearly perpendicular direction 
c (that is, along gradient of density), 
c which does not seem to be natural.
c Now the ray steps-in along the Vgroup direction,
c and so it continues smoothly along same path.
c It might happen, though, 
c that the ray will glide tangentially away
c from the surface, i.e. going away from Xe=1 layer;
c in such case it will continue as an O-mode.
c User can still revert the procedure to
c stepping along grad(ne) across the cutoff layer:
c see step_dir=... option in subroutine find_rho_X_xyz .

c [11] Added the plots of power deposition in (R,Z) plane.
c The local power deposition is saved over uniform (R,Z) grid.
c The grid size could be explicitly defined in &tokamak namelist,
c through NRgrid and NZgrid names. The default values are 
c NRgrid=60 and NZgrid=190.
c User might want to use zoom function in Python figure (find a
c pictogram of magnifying glass in the plot menu).
c When the power depositions are strongly localized,
c it is difficult to notice the profile of deposition,
c which is plotted over ray trajectories, 
c and so it could be masked by them.
c Also, it should be noted that in these plots
c (genray_profiles_power_RZ.png)
c the rays are plotted using (R(t),Z(t)) arrays ,
c where R(t) = sqrt(X(t)^2 + Y(t)^2) 
c is the major radius coordinate along ray,
c which is always positive.
c This is done because the power of ray deposition 
c is mapped onto the (R,Z) plane, 
c no matter what the toroidal angle is.
c (Otherwise, the power deposition 
c had to be saved over (X,Y,Z) 3D grid.)
c Because of always-positive value of R(t),
c sometimes the rays may look "reflected" from 
c R=0 line, when in fact they cross it and continue
c on the other side.

c [10] Added plots related to O-X transmission coefficient.
c See plots genray_rays_OX_transmission.png
c produced by genray_plot.py.
c Subplot "Indicator of OX conversion" 
c shows value of 0.0 for rays that were not converted to X-mode,
c and value of 1.0 if ray was converted 
c (even if the transmission coefficient is low). 
c Each ray is identified by value of Npar at ray starting point
c (which is the horizontal axis),  and also by color.
c The value of Npar just before the O-X "jump"
c is shown at the bottom subplot. 
c The value of transmission coefficient is shown at the right-top
c subplot; value of 1.0 corresponds to complete conversion
c of O-mode into X-mode (no reflected power).

c----------------------------------------
c [2014-10-22] version="Genray-c_141022"
c----------------------------------------

c [9] Plotting script genray_plot.py is reworked - 
c now it is compatible with Python on Mac
c (there were some issues with contour plots).
c Also, few more plots are added. 

c [8] Subroutine DHOT_s (used for id=6)
c is reworked to eliminate 1/npar divergence.
c New SUBROUTINE ZFUN_vkw(resn, vkw, CZ0_vkw, CZ0_resn_vkw)
c is added to calculate Zn/(npar*beta) and resn*Zn/(npar*beta)
c to avoid the divergence ( resn=(omega-n*omega_c)/omega, 
c and vkw=(npar*beta)== Vthermal*Kpar/omega ).
c Also, alternative definition of summation range over harmonics
c is added into DHOT_s, which reduces the number of required harmonics.
c Now the rays with very small npar~0.01 can be traced (O-X-EBW).
c There is still an accumulated loss of accuracy for smaller npar:
c although the rays with a very small npar (0.001) can be launched,
c they are stopped at a resonance (such as w=wc) because of degraded
c accuracy. Usually such rays are of no significance because
c most of their power is lost at O-X conversion layer 
c (low transmission coefficient).

c [7] Subroutine find_rho_X_xyz for searching the X-root  
c for O-X conversion is improved.
c The procedure recalculates grad(n) at each iteration during
c stepping towards X=1.0 (until it finds X-mode).
c The control is done by checking the value of local rho,
c not r as before.

c [6] New input in the namelist (section &tokamak ):
c "eqdsktype=" (character*16) Type of eqdsk data file 
c (for model_b=0, i.e., when reading the eqdsk data file).
c For TAE (Tri Alpha Energy) FRC case, the eqdsk file 
c does not contain qpsi array, but instead it contains 
c necut and tmcut arrays (electron density [1/m^3] and
c el.temperature [eV]). So, one extra read block is implemented.
c Options so far: "TAE" (default) or "tokamak" .

c [5] New input in the namelist (section &tokamak ):
c "rlim_wall_fraction=" !==rlim/wall_rmax (for model_b= 1 or 2)
c It is used to define the value of rlim as a fraction
c of chamber wall radius:
c rlim= rlim_wall_fraction*wall_rmax
c rlim corresponds to psilim and to "rho"=1.
c Although for model_b=1 or 2 there are no flux surfaces,
c we still need to define the edge of plasma, 
c for definition of plasma profiles n and T.  
c This is done by setting the value of rlim; 
c the value of psilim is found by integration 
c of poloidal (Z) flux over radius R (at fixed Z=0 coord.)
c psilim corresponds to "rho" such that rho(psilim)=1.
c So, the value of rlim sets the edge of hot plasma.
c The profiles of n and T are set to drop exponentially 
c or linearly at psi>psilim ("rho">1).
c Default value: 0.82 .

c [4] Definition of rho_bin(), rho_bin_center(), binvol, etc. is moved
c to subroutine onetwoini.
c These 1D arrays are stored in onetwo_no_nml.i and reused when needed.
c Before, many of them where calculated locally in different subroutines.    
c These arrays are defined even when there are no flux surfaces 
c (mirror machine, FRC, etc),
c with rho_bin extended to the chamber wall (rho can be larger than 1.0)
c or to the edge of equilibrium-B grid, whichever is smaller.

c [3] New option is added: density profile can be set as den(rho)*gn(z)
c where gn(z) is the dropoff function in Z-direction.
c See GENRAY-c_help.txt, about definition of
c zbegin_den_dropoff  and  zlength_den_dropoff .

c [2] Now (ws(is)-ws(is-1))) is dl along ray.
c Before, it was the "poloidal" projection of the ray element,
c which may not be always relevant, if Bpol=0.
c The plots of "delpwr" (change of power == amplitude^2 of a ray)
c are made versus distance along ray (integral of dl);
c before - they were made versus poloidal distance.

c [1] Procedure for calculating the plasma dispersion roots 
c is revised. Now all values can be specified in the namelist,
c in section &output. 
c See GENRAY-c_help.txt, look for "i_save_disp=' and other 
c names following this option (they all have suffix "_disp").

c [0] Starting with version="Genray-c_120405"
