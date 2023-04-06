c--------------------------------------------------------
c    namelist /tokamak/
c-------------------------------------------------------------
      namelist /tokamak/ model_b,indexrho,ipsi,ionetwo,ieffic,
     + psifactr, bz0, rbphi0, wall_rmin,wall_rmax,wall_zmin,wall_zmax,
     + rlim_wall_fraction,
     + rs_frc, ! separatrix radius for FRC-like plasmas
     + curc,radc,zcoil,
     + rmirror,zbox_mirror,rbox_mirror,b00_mirror,
     + eqdskin,eqdsktype,NR, NRgrid, NZgrid,
     & ncoils, n_wall,max_limiters,n_limiter,r_wall,z_wall,
     & r_limiter,z_limiter,phi_limiter,
     & h_add_wall
