
c      the data for the profiles of the absorbed power and current drive.
c      PARAMETER (NRA) is set in param.i
c      No namelist data
   
c-----real*8
      real*8 powtot_e,powtot_i,powtot_cl,allpower,
     &       allpw_e,allpw_i,allpw_cl,
     &       allcur,currtot,
     &       parallel_cur_total,toroidal_cur_total,poloidal_cur_total
      real*8,  pointer ::
     &         powtot_s(:),       !(nbulka),
     &         power(:),          !(NRA)
     &         current(:),        !(NRA)
     &         temparr(:),        !(NRA)
     &         zeffarr(:),        !(NRA)
     &	       spower(:),         !(NRA)
     &         scurrent(:),       !(NRA)
     &         powden(:),         !(NRA)
     &         currden(:),        !(NRA) 
     &         power_e(:),        !(NRA) 
     &         power_i(:),        !(NRA) 
     &         power_s(:,:),      !(NRA,nbulka),
     &         power_cl(:),       !(NRA) 
     &         spower_e(:),       !(NRA) 
     &         spower_i(:),       !(NRA) 
     &         spower_s(:,:),     !(NRA,nbulka),
     5	       spower_cl(:),      !(NRA) 
     &         powden_e(:),       !(NRA) 
     &         powden_i(:),       !(NRA) 
     &         powden_s(:,:),     !(NRA,nbulka)
     &         powden_cl(:),      !(NRA) 
     &         currden_s(:),      !(NRA) 
     & 	       powden_e_s(:),     !(NRA) 
     &         powden_i_s(:),     !(NRA) 
     &         powden_cl_s(:),    !(NRA)         
     & 	       densprof(:,:),     !(NRA,nbulka)
     &         temprof(:,:),      !(NRA,nbulka)
     &         zefprof(:),        !(NRA)           
     +         bmodprofxz(:,:), !(nxeqd,nzeqd)
     +         densprofxz(:,:), !(nxeqd,nzeqd)
     +         densprofyz(:,:), !(nyeqd,nzeqd)
     +         densprofxy(:,:), !(nyeqd,nzeqd)
     &         rho_bin(:),        !(NRA)     
     &         rho_bin_center(:), !(NRA-1)     
     &         binvol(:),         !(NRA-1) 
     &         binarea(:),        !(NRA-1) 
     &         binarea_pol(:),    !(NRA-1) 
     &         pollen(:),         !(NRA-1) 
     &         cur_den_parallel(:),          !(NRA)   
     &         s_cur_den_parallel(:),        !(NRA-1)   
     &         s_cur_den_onetwo(:),          !(NRA-1)   
     &         s_cur_den_toroidal(:),        !(NRA)   
     &         s_cur_den_poloidal(:),        !(NRA)  
     &         Rgrid(:),                          !(NRgrid)
     &         Zgrid(:),                          !(NZgrid)
     &         pwr_rz_e(:,:),  spwr_rz_e(:,:),  !(NRgrid,NZgrid)
     &         pwr_rz_i(:,:),  spwr_rz_i(:,:),  !(NRgrid,NZgrid)
     &         pwr_rz_cl(:,:),  spwr_rz_cl(:,:),  !(NRgrid,NZgrid)
     &         pwr_rz_s(:,:,:), spwr_rz_s(:,:,:)  !(NRgrid,NZgrid,nbulk)
  


      common/onetwo_no_nml/ powtot_e,powtot_i,powtot_cl,allpower,
     &         allpw_e,allpw_i,allpw_cl,
     &         allcur,currtot,  
     &         parallel_cur_total,toroidal_cur_total,poloidal_cur_total,
c-----pointers
     &         powtot_s,
     &         power,
     &         current,
     &         temparr,
     &         zeffarr,
     &	       spower,
     &         scurrent,
     &         powden,
     &         currden,
     &         power_e,
     &         power_i,
     &         power_s,
     &         power_cl,
     &         spower_e,
     &         spower_i,
     &         spower_s,
     &	       spower_cl,
     &         powden_e,
     &         powden_i,
     &         powden_s,
     &         powden_cl,
     &         currden_s,
     & 	       powden_e_s,
     &         powden_i_s,
     &         powden_cl_s,
     8	       densprof,
     &         temprof,
     &         zefprof,
     +         bmodprofxz,densprofxz,densprofyz,densprofxy,
     &         rho_bin,
     &         rho_bin_center,
     &         binvol,
     &         binarea,
     &         binarea_pol,
     &         pollen,
     &         cur_den_parallel,
     &         s_cur_den_parallel,
     &         s_cur_den_onetwo,
     &         s_cur_den_toroidal,
     &         s_cur_den_poloidal,
     &         Rgrid,     
     &         Zgrid,     
     &         pwr_rz_e, spwr_rz_e, 
     &         pwr_rz_i, spwr_rz_i, 
     &         pwr_rz_cl, spwr_rz_cl, 
     &         pwr_rz_s,  spwr_rz_s
    
c-----------------------------------------------------------------------------
c    NR is the maximal number of points in the small radius direction 
c      for the
c      calculations of the power and current drive radial profiles 
c      in each radial bin.
c    NRA= maximal value of NR. It is set in param.i
c    power_e is the absorbed power due to collisionless electron damping
c    power_i is the absorbed power due to collisionless ion damping
c    power_s is the absorbed power due to coll'less ion damping on each ion
c    power_cl is the absorbed power due to collisional damping
c    spower_e summed absorbed power due to coll'less electron damping
c    spower_i summed absorbed power due to collisionless ion damping
c    spower_s summed absorbed power due to coll'less ion damping on each ion
c    spower_cl summed absorbed power due to collisional damping
c    allcur     total current
c    densprof,temprof,zefprof   for output of radial profiles to .nc file
c
c    For output of electron density profiles in xz and yz-plane to .nc file
c    densprofxz(:,:) !(nxeqd,nzeqd)
c    densprofyz(:,:) !(nyeqd,nzeqd)
c    densprofxy(:,:) !(nxeqd,nyeqd)
c    For plots of total B (equilibrium) :
c    bmodprofxz(:,:), !(nxeqd,nzeqd)
c
c
c050406
c    currtot is a total current along B field
c
c    cur_den_parallel is the flux surface averaged parallel current density
c                     from one ray
c                    
c   s_cur_den_parallel is the flux surface averaged parallel current density
c                     from all rays
c
c   binarea_pol is bin area of the poloidal cross-section at theta poloidal=0 

c   YuP[Nov-2014] 
c    For power deposition profiles over (R,Z) rectangular grid.
c    For each ray (re-used for each ray): 
c   pwr_rz_s(NRgrid,NZgrid,nbulk) is coll-less damping for each species
c   pwr_rz_cl(NRgrid,NZgrid) is collisional damping.
c    Sum over all rays: 
c   spwr_rz_e(NRgrid,NZgrid)        !coll-less damping for e
c   spwr_rz_i(NRgrid,NZgrid)        !coll-less damping for i
c   spwr_rz_cl(NRgrid,NZgrid)       !collisional damping
c   spwr_rz_s(NRgrid,NZgrid,nbulk)  !coll-less damping for each species
c    (R,Z) grid for power deposition profiles:
c   Rgrid(NRgrid)
c   Zgrid(NZgrid)
