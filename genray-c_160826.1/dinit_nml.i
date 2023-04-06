c----------------------------------------
c     declaration of namelist variables which
c     are used in subroutine dinit_mr and
c     are not used in any common block
c-----------------------------------
      integer ndim1,                ! /numercl/
     &nj_tab(nbulka)                ! number of radial points of nonuniform
                                    ! profiles 

      real*8
     &prof(nbulka*ndensa),          !Working array for uniform and non-uniform
                                    !profiles
     &prof2(nbulka,ndensa),         !Working array for uniform and
                                    !and non-uniform profile          
     &prof_2d(ndensa,nbulka),       !Nonuniform table profile
     &radii_2d(ndensa,nbulka)       !given by lines at radii_2d mesh
c-----------------------------------------------------------------------------
c     ndim1=6 (number of the ray tracing equations)
c     prmt1=prmt(1)= initial time for a ray ! Not needed. 0 by default.
c     prmt2=prmt(2)= largest allowed time for ray advancing  (Not used)
c     prmt3=prmt(3)= initial time step for integration (normalized units)
c     prmt4=prmt(4)= required accuracy
c     prmt6=prmt(6)  [m] distance step for saving ray data; 
c                        MAY AFFECT the step of integration!
c-----------------------------------------------------------------
c     namelists for all table plasma profiles at uniform
c     radial mesh written by columns:
c-----------------------------------------------------------------
c     namelist /dentab/ prof 
c     namelist /temtab/ prof
c     namelist /tpoptab/ prof
c     namelist /vflowtab/ prof
c     namelist /zeftab/ zeff1
c-----------------------------------------------------------------
c     namelists for all table plasma profiles at non uniform
c     radial mesh written by columns:
c-----------------------------------------------------------------
c     namelist /dentab_nonuniform/   nj_tab,prof,prof_radii
c     namelist /temtab_nonuniform/   nj_tab,prof,prof_radii
c     namelist /tpoptab_nonuniform/  nj_tab,prof,prof_radii
c     namelist /vflowtab_nonuniform / nj_tab,prof,prof_radii
c     namelist /zeftab_nonuniform /   nj_tab,prof,prof_radii
c-----------------------------------------------------------------
c     namelists for all table plasma profiles at non uniform
c     radial mesh written by lines:
c-----------------------------------------------------------------
c     namelist /dentab_nonuniform_line/   nj_tab,prof_2d,radii_2d
c     namelist /temtab_nonuniform_line/   nj_tab,prof_2d,radii_2d
c     namelist /tpoptab_nonuniform_line/  nj_tab,prof_2d,radii_2d
c     namelist /vflowtab_nonuniform_line / nj_tab,prof_2d,radii_2d
c     namelist /zeftab_nonuniform_line /   nj_tab,prof_2d,radii_2d
c------------------------------------------------------------------
