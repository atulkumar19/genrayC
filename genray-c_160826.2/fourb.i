
c the arrays from the input eqdskin (default="equilib.dat") file
      real*8 necut,tmcut, peqd,feqd,pres,ffpeqd,
     1 ppeqd,qpsi, req,xeq,yeq,zeq, flux, beq
     
      real*8, pointer :: rlimit(:), zlimit(:) ! (nlimit) or (nnlim)
      real*8, pointer :: rves(:),   zves(:)   ! (nves) or (nnves)

c	Set up some variables for re-gridding of input f(1:nveqd), etc.,
c       arrays onto a equispaced mesh of dimension nreqd.
c       BobH020822.
      integer iworkka,itabl,i1p,n_wall_add,n_limiter_add

      parameter(iworkka=3*nreqda+1)
      real*8 psiar,d2feqd,d2pres,d2ffpeqd,d2ppeqd,d2qpsi,workk,r8temp,
     &tabl,
     & thetapol_wall,thetapol_limiter,
     & rho_wall,rho_limiter,
     &r_wall_add,z_wall_add,
     &rr_add,zz_add,distance_to_wall,
     &density_r_z,                   !density at rz mesh
     &density_r_z_rr,density_r_z_zz, !second derivatives
     &density_r_z_rrzz               !4th derivatives
    
     

      common/fourb/ feqd(nreqda),pres(nreqda),ffpeqd(nreqda),
     1             ppeqd(nreqda),
     +             necut(nreqda,nzeqda),tmcut(nreqda,nzeqda),
     +             peqd(nreqda,nzeqda), beq(nreqda,nzeqda),
     2             qpsi(nreqda),
     5             req(nreqda),xeq(nxeqda),yeq(nyeqda),zeq(nzeqda),
     6             flux(nreqda),
c------BobH02022
     &             psiar(nreqda),d2feqd(nreqda),d2pres(nreqda),
     1             d2ffpeqd(nreqda),d2ppeqd(nreqda),d2qpsi(nreqda),
     2            r8temp(nreqda),i1p(2),workk(iworkka),tabl(3),itabl(3),
c-----for wall and limiter
     &             thetapol_wall(n_wall_a),
     &             thetapol_limiter(n_limiter_a,max_limiters_a),
     &             rho_wall(n_wall_a),
     &             rho_limiter(n_limiter_a,max_limiters_a),
     &             r_wall_add( n_wall_add_a,0:max_limiters_a),
     &             z_wall_add( n_wall_add_a,0:max_limiters_a),
c-----for density fall near the wall
     &             rr_add(nreqd_add_a),zz_add(nzeqd_add_a),
     &   distance_to_wall(nreqd_add_a,nzeqd_add_a,0:max_limiters_a),
     &   density_r_z(nreqd_add_a,nzeqd_add_a,nbulka,0:max_limiters_a),
     &  density_r_z_rr(nreqd_add_a,nzeqd_add_a,nbulka,0:max_limiters_a),
     &  density_r_z_zz(nreqd_add_a,nzeqd_add_a,nbulka,0:max_limiters_a),
     &density_r_z_rrzz(nreqd_add_a,nzeqd_add_a,nbulka,0:max_limiters_a),
     &   n_wall_add(0:max_limiters_a)

      common/fourb_pointer/ rlimit,zlimit,rves,zves
      
c
c  nreqda must be >= nreqd
c  nzeqda must be >= nzeqd
c
c----------------------------------------------------------------
c thetapol_wall,       ! poloidal angles [radians] of
c thetapol_limiter     ! wall and limiter points     
c rho_wall,rho_limiter !small radius
c


c Input data from dendsk file:
      integer        nxden,nyden
      real*8         denmin,denmax
      real*8         xdenmin, xdenmax, ydenmin, ydenmax
      real*8         xden(nxdena),yden(nydena),dengrid(nxdena,nydena)
      common/dendsk_xy/
     1               nxden,nyden,
     1               denmin,denmax, !To be found: dengrid min/max [m^-3]
     2               xdenmin, xdenmax, ydenmin, ydenmax,
     3               xden, ! x-uniform-grid
     3               yden, ! y-uniform-grid
     4               dengrid