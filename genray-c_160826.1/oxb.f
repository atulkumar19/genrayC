

c     The subroutines for the optimaL OXB launch.
      subroutine antenna_surface(theta,z_ant,r_ant)
c     determine the possible antenna position 

c     z_ant(theta), r_ant(theta) [m]
c     theta is a poloidal angle [rad]
  
c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none

      include 'three.i'
c----------------------------------------
      include 'param.i'                !for i_ant=2 case
      include 'five.i' !zmax,rzmax     !
      include 'cone.i' !zst(),rst()    !
c------------------------------------------
      include 'antenna.i'
  
      logical first !used in i_ant=1 case
      data first / .true./
      save first

c-----input
      real*8 theta    ! major radius [m]
    
c-----output
      real*8 z_ant,r_ant

c-----local
      real*8 p,p1,p2,det, det_a,det_b,a_sq,b_sq,
     & accuracy,p_intersection,r_intersection,z_intersection,
     &distance,theta_l,pi

      
c-----externals
      real*8 rtbis,psilim_psif_zp_rp
      external psilim_psif_zp_rp
 
      pi=4.d0*datan(1.d0)
      write(*,*)'antenna_surface theta[degree]',theta*180.d0/pi
      
      if((i_ant.lt.1).or.(i_ant.gt.2)) then
        write(*,*)'oxb.f in function z_antenna(r)'
        write(*,*)'(i_ant.lt.1).or.(i_ant.gt.2)' 
        write(*,*)'i_ant=',i_ant
        write(*,*)'Please check i_ant'
        stop 'in oxb.f'
      endif


      if(i_ant.eq.2) then
c------------------------------------------------------------------------
c       EC cone vertex is at the ellipse curve curve: z_antenna(r)
c
c       (z_antenna-zma)**2  (r-rma)**2 
c       --------------    + ---------- =1
c            A**2               B**2
c
c       Let the following points are at the ellipse curve:
c       (zst,rst) the given EC cone vertex
c       (zmax+delta_z,rzmax) This point is vertically shifted over the point 
c                           with with maximal Z at the limiter line
c                           (zmax,rzmax)
c                           
c       (zst-zma)**2       (rst-rma)**2 
c       --------------  + -------------- =1
c          A**2               B**2
c
c
c       (zmax+delta_z-zma)**2       (rzmax-rma)**2 
c       ---------------------    + --------------- =1
c          A**2                        B**2
c
c
c       It is the system for the coefficients A**2 and B**2
c
c       Det=(zst-zma)**2 * (rzmax-rma)**2 -
c           (rst-rma)**2 * (zmax+delta_z-zma)**2
c
c       Det_A=  (rzmax-rma)**2 - (rst-rma)**2
c
c       Det_B=  (zst-zma)**2 - (zmax+delta_z-zma)**2  
c
c       (1/A**2)= Det_A / Det 
c       (1/B**2)= Det_B / Det
c
c        Z(r)= zma +/- sqrt{A**2 -(r-rma)**2*(A**2/B**2)}  
c
c        z(theta)=zma+A*sin(theta)
c        r(theta)=rma+B*cos(theta)

c-------------------------------------------------------------------------
c        i_cone is the number of the EC cone i_cone=1,...,nconea

        det=(zst(i_cone)-zma)**2*(rzmax-rma)**2-
     &      (rst(i_cone)-rma)**2*(zmax+delta_z-zma)**2
 
        det_a=(rzmax-rma)**2-(rst(i_cone)-rma)**2
        det_b=(zst(i_cone)-zma)**2-(zmax+delta_z-zma)**2  
        write(*,*)'z_antenna det_a,det_b',det_a,det_b

        a_sq= det/det_a
        b_sq= det/det_b
        write(*,*)'z_antenna a_sq,b_sq',a_sq,b_sq

        z_ant=zma+dsqrt(a_sq)*dsin(theta)
        r_ant=rma+dsqrt(b_sq)*dcos(theta)
      endif  !i_ant.eq.2

      if (i_ant.eq.1) then
c---------------------------------------------------------------------
c        The antenna surface will be at the same given distance from
c        the plasma edge. This distances is along the straight line 
c        from the magnetic axis. The distance is equal to the
c        distance of the given EC cone vertex from the plasma edge.
c-------------------------------------------------------------------
c         r-rma      z-zma
c        ------- = --------=p is the line through the magnetic axis (rma,zma)
c        rst-rma    zst-zma   and the EC cone vertex (rst,zst)
c
c        r=rma + (rst-rma)p 
c        z=zma + (zst-zma)p
c
c        or   z=zma + (r-rma)*(zst-zma)/(rst-rma) 
c
c        Plasma boundary: rmin < r < rmax
c
c        The point of this line intersection with the plasma boundary
c
c        psilim-psif(z_intersecsion,r_intersection) =0 
c
c        Here z_intersection=rma+(rst-rma)p_intersection
c             r_intersection=zma+(zst-zma)p_intersection
c---------------------------------------------------------------------
c         write(*,*)'antenna_surface first=',first
        
         if (first) then
c-----------create the antenna surface          
            do iray=1,ncone
c--------------the loop by all EC cone vertexs
                write(*,*)'antenna_surface iray=',iray
c-----------------------------------------------------
c              calculate intersection [r_intersection,z_intersection]
c              of plasma boundary with the stright line
c              from the magnetic axis to the given EC cone vertex
c------------------------------------------------------
               i_cone=iray !set this parameter to antenna.i
               accuracy=1.d-7
               p1=0.d0
               p2=1.d0
               write(*,*)'antenna surface before p_intersection=rtbis'
               p_intersection= rtbis(psilim_psif_zp_rp,p1,p2,accuracy)
c               write(*,*)'antenna surface after p_intersection=',
c     &         p_intersection
c               write(*,*)'rma,rst(i_cone),p_intersection',
c     &         rma,rst(i_cone),p_intersection
               r_intersection=rma+(rst(i_cone)-rma)*p_intersection
               z_intersection=zma+(zst(i_cone)-zma)*p_intersection
c               write(*,*)'antenna surface z_intersection,
c     &         r_intersection',z_intersection,r_intersection
c            write(*,*)'antenna surface i_cone,zst(i_cone),rst(i_cone)',
c     &      i_cone,zst(i_cone),rst(i_cone)
               distance=dsqrt((zst(i_cone)-z_intersection)**2+
     &                        (rst(i_cone)-r_intersection)**2)          
               write(*,*)'antenna surface z_intersection,
     &         r_intersection,distance',z_intersection,
     &         r_intersection,distance
c------------------------------------------------------------------
c              calculate spline coefficients for the antenna points
c              at the surface which is at the equal distance from
c              the plasma {z_an(theta),r_ant(theta)}
c------------------------------------------------------------------  
c               write(*,*)'antenna surface before antenna_equal_dist'
c               write(*,*)'first,i_cone,distance,theta,z_ant,r_ant',
c     &         first,i_cone,distance,theta,z_ant,r_ant
               theta_l=theta
               call antenna_equal_dist(first,i_cone,distance,
     &                                 theta_l,z_ant,r_ant)
               write(*,*)'********************************'
               write(*,*)'antenna surface was created antenna surface'
               write(*,*)'********************************'
            enddo !iray   
             
            first=.false.
         endif !first
c--------calculate the antenna vertex coordinates z_ant,r_ant
c        at the antenna surface for the given poloidal angle theta
c        using spline

c         write(*,*)'antenna surface before antenna_equal_dist'
c         write(*,*)'******************************************'

c         write(*,*)'first,i_cone,distance,theta',
c     &         first,i_cone,distance,theta 
         theta_l=theta
         call antenna_equal_dist(first,i_cone,distance,
     &                           theta_l,z_ant,r_ant)
c         write(*,*)'antenna surface after antenna_equal_dist'
c         write(*,*)'theta,z_ant,r_ant',theta,z_ant,r_ant
      endif  !i_ant.eq.3

      return      
      end


c======================================================================
c======================================================================


      subroutine transform_N_rphiz_to_xyz(n_r,n_phi,phi,n_x,n_y)
c-----transform the vector coordinates N_r,N_phi to 
c     N_x,N_y

      implicit none

c-----input     
      real*8 n_r,n_phi,     !cylindrical coordinates
     &phi                   !toroidal angle [radians]

c-----output
      real*8 n_x,n_y        !x,y coordinates

c-----local
      real*8 sin_phi,cos_phi

      sin_phi=dsin(phi)
      cos_phi=dcos(phi)

      n_x=n_r*cos_phi-n_phi*sin_phi
      n_y=n_r*sin_phi+n_phi*cos_phi

      return
      end

c======================================================================
c======================================================================


      subroutine edg_vac_refr_index_rphiz(n_r,n_phi,n_z,dpsi_dz,dpsi_dr,
     &n_vac_r,n_vac_phi,n_vac_z)
c-----calculate vacuum refractive index at the plasma edge
c     n_vac_r,n_vac_phi,n_vac_z
c     using the refractive index in plasma at the plasma edge
c     n_r,n_phi,n_z 
c     and the poloidal flux gradient
c     dpsi_dz,dpsi_dr,

      implicit none
c-----input
      real*8 n_r,n_phi,n_z,  !refractive index at the plasma boundary 
                             !in plasma side
     &dpsi_dz,dpsi_dr        !gradiend of the poloidal flux        
                             !at the plasma boundary
c-----output
      real*8 n_vac_r,n_vac_phi,n_vac_z  ! refractive index at the plasma
                                        ! boundary in vacuum side

c-----local
      real*8 p,

     &gradpsi_mod,
     &n_gradpsi_pl_r,n_gradpsi_pl_z,      ! refractive index along the vector
     &n_gradpsi_pl_mod,                   ! grad_psi at plasma side

     &n_psi_r,n_psi_z,                    ! refractive index along 
     &n_psi_mod_s,                        ! the magnetic surface

     &n_gradpsi_vac_r,n_gradpsi_vac_z,    ! refractive index along the vector
     &n_gradpsi_vac_mod                   ! grad_psi at vacuum side

c-----refractive index along grad_psi at plasma side
      gradpsi_mod=dsqrt(dpsi_dr**2+dpsi_dz**2)
      n_gradpsi_pl_mod=(n_r*dpsi_dr+n_z*dpsi_dz)/gradpsi_mod
      n_gradpsi_pl_r=n_gradpsi_pl_mod*dpsi_dr/gradpsi_mod
      n_gradpsi_pl_z=n_gradpsi_pl_mod*dpsi_dz/gradpsi_mod

c-----refractive index along the magnetic surface 
c     at plasma side and at vacuum side
      n_psi_r=n_r-n_gradpsi_pl_r
      n_psi_z=n_z-n_gradpsi_pl_z
      n_psi_mod_s=n_psi_r**2+n_psi_z**2+n_phi**2

c-----refractive index along grad_psi at vacuum side
      p=1.d0-n_psi_mod_s
      if (p.gt.0.d0)then
          n_gradpsi_vac_mod=dsqrt(p)
      else
          write(*,*)'*****WARNING****'
          write(*,*)'oxb  edg_vac_refr_index np_psi_mod_s>1'
          write(*,*)'n_psi_mod_s',n_psi_mod_s
          n_gradpsi_vac_mod=0.d0
c          stop 'oxb'
      endif

c      write(*,*)'edg vac_refr_index_rphiz'
c      write(*,*)'n_gradpsi_vac_mod**2+n_psi_mod_s',
c     &n_gradpsi_vac_mod**2+n_psi_mod_s

      n_gradpsi_vac_r=n_gradpsi_vac_mod*dpsi_dr/gradpsi_mod
      n_gradpsi_vac_z=n_gradpsi_vac_mod*dpsi_dz/gradpsi_mod

c-----vacuum refractive index
      n_vac_r= n_psi_r + n_gradpsi_vac_r
      n_vac_z= n_psi_z + n_gradpsi_vac_z
      n_vac_phi=n_phi

c     write(*,*)'edg_vac_refr_index_rphiz'
c     write(*,*)'dsqrt(n_vac_r**2+n_vac_z**2+n_vac_phi**2)',
c    &dsqrt(n_vac_r**2+n_vac_z**2+n_vac_phi**2)
      return
      end

c======================================================================
c======================================================================


      subroutine edg_vac_refr_index(n_r,n_phi,n_z,phi,
     &dpsi_dz,dpsi_dr,
     &n_vac_x,n_vac_y,n_vac_z)
c-----calculate (x,y,z) components of the vacuum refractive index
c     at the plasma boundary
     
      implicit none

c-----input
      real*8 n_r,n_phi,n_z,  !refractive index at the plasma boundary 
                             !in plasma side
     &phi,                   !toroidal angle [radian] 
     &dpsi_dz,dpsi_dr        !gradiend of the poloidal flux        
                             !at the plasma boundary
c-----output
      real*8 n_vac_x,n_vac_y,n_vac_z  !refractive index at the plasma boundary 
                                      !in vacuum side
c-----local
      real*8 n_vac_r,n_vac_phi        ! refractive index at the plasma
                                      ! boundary in vacuum side

c-----calculate vacuum refractive index at the plasma edge
c     n_vac_r,n_vac_phi,n_vac_z
c     using
c     the refractive index in plasma side at the plasma edge
c     n_r,n_phi,n_z  
c     and the poloidal flux gradient components at the plasma edge
c     dpsi_dz,dpsi_dr    

      call edg_vac_refr_index_rphiz(n_r,n_phi,n_z,dpsi_dz,dpsi_dr,
     &n_vac_r,n_vac_phi,n_vac_z)

c-----calculate vacuum refractive index at the plasma edge
c     n_vac_x,n_vac_y,n_vac_z
c     using 
c     the cylindrical coordinates: n_vac_r,n_vac_phi
c     and the toroidal angle: phi
     
      call transform_N_rphiz_to_xyz(n_vac_r,n_vac_phi,phi,
     &n_vac_x,n_vac_y)
       
      return
      end


c======================================================================
c======================================================================


      subroutine antenna_vertex(r_0,phi_0,z_0,n0_r,n0_phi,n0_z,
     &b_r0,b_phi0,b_z0,dpsi_dz,dpsi_dr,
     & r_st,phi_st,z_st,alpha_st,beta_st,icone)
c-----calculate O mode cone vertex
c     r_st,phi_st,z_st
c     and the central ray direction
c     alpha_st,beta_st
c
c     As the intersection of the vacume ray (strigth line)
c     starting in the point (r_0,phi_0,z_0) in the direction    
c     determined by vacuum refractive index
c
      implicit none

c-----input
      real*8 r_0,phi_0,z_0, !space location of the O-mode 
                            !at the plasma boundary [m]
     &n0_r,n0_phi,n0_z,     !o-mode refractive index at the plasma boundary
                            ! in plasma side
     &b_r0,b_phi0,b_z0,     !the magnetic field in  M_0(r_0,phi_0,z_0)
     &dpsi_dz,dpsi_dr       !gradiend of the poloidal flux        
                            !at the plasma boundary
      integer icone         !the number of EC cone vertex =1,...,ncone
                            !ncone =< nconea
c-----------------------------------------------------------------
      real*8 r_min,r_max
c-------------------------------------------------------------------
      include 'antenna.i'      
c--------------------------------------------------------------------   

c-----output
      real*8 r_st,phi_st,z_st, !EC cone vertex [m, degree]
     &alpha_st,!central ray toroidal angle [degree]
     &beta_st !central ray poloidal angle
     
c-----externals
      real*8 x_bin_min,f_antenna_min
      external f_antenna_min

c-----locals
      real*8 n_vac_x,n_vac_y,n_vac_z,
     &x_st,y_st,pi,gamma,delta,trans,
     &x_0,y_0
      double precision
     &pacc,p_antenna,p1,p2
     
c-----for test
      real*8 r_p,z_p,z_antenna_p,z_antenna 
      real*8 cn0_s, cnteta,cnphi,cnrho

      write(*,*)'in antenna_vertex'
c------------------------------------
      i_ant=1 ! equal distance from antenna until plasma
c------------------------------------
c      i_ant=2       ! ellipse curve for EC cone vertex positions
c      i_cone=1      ! the number od EC cone =1,...,nconea
      i_cone=icone  
      delta_z=0.2d0 ! vertical shift over (zmax,rzmax) point
c--------------------------------------
c      write(*,*)'r_0,phi_0,z_0,n0_r,n0_phi,n0_z',
c     &r_0,phi_0,z_0,n0_r,n0_phi,n0_z
c      write(*,*)'b_r0,b_phi0,b_z0,dpsi_dz,dpsi_dr',
c     &b_r0,b_phi0,b_z0,dpsi_dz,dpsi_dr
c-----------------------------------------------------------------
c     calculate (x,y,z) components of the vacuum refractive index
c     n_vac_x,n_vac_y,n_vac_z
c     at the plasma boundary: M_0(r_0,phi_0,z_0)
c-----------------------------------------------------------------
      call edg_vac_refr_index(n0_r,n0_phi,n0_z,phi_0,
     &dpsi_dz,dpsi_dr,  
     &n_vac_x,n_vac_y,n_vac_z)
c      write(*,*)'after edg_vac_refr_index'
c      write(*,*)'n0_r,n0_phi,n0_z,phi_0',n0_r,n0_phi,n0_z,phi_0
c      write(*,*)'dpsi_dz,dpsi_dr',dpsi_dz,dpsi_dr
c      write(*,*)'n_vac_x,n_vac_y,n_vac_z',n_vac_x,n_vac_y,n_vac_z
c      write(*,*)'total vac n',dsqrt(n_vac_x**2+n_vac_y**2+n_vac_z**2)
c---------------------------------------------------------------------
c     test calculations cnteta,cnphi at vacuum side
c---------------------------------------------------------------------
cyup      call ninit_ec(z_0,r_0,phi_0,n_vac_x,n_vac_y,n_vac_z,cnteta,cnphi)
c      write(*,*)'cnteta,cnphi',cnteta,cnphi
cyup      cn0_s=n_vac_x**2+n_vac_y**2+n_vac_z**2
c      write(*,*)'cn0_s',cn0_s
cyup      cnrho=dsqrt(cn0_s-cnteta**2-cnphi**2)
c      write(*,*)'cnrho', cnrho
c----------------------------------------------------------------
c     calculate the vertex coordinates
c----------------------------------------------------------------
      call set_common_ox(r_0,phi_0,z_0,
     &n_vac_x,n_vac_y,n_vac_z)
     
      x_0=r_0*dcos(phi_0)
      y_0=r_0*dsin(phi_0) 
     
      pacc=1.d-10
       
      p1=0.d0
      p2=10.d0 !to check how to find the left boundary
c      write(*,*)'p1,p2',p1,p2

cfor test
c      z_p=z_0+n_vac_z*p1
c      r_p=dsqrt((x_0+p1*n_vac_x)**2+(y_0+p1*n_vac_y)**2)
c      write(*,*)'p1,z_p,r_p',p1,z_p,r_p

c      z_p=z_0+n_vac_z*p2
c      r_p=dsqrt((x_0+p2*n_vac_x)**2+(y_0+p2*n_vac_y)**2)
c      write(*,*)'p2,z_p,r_p',p2,z_p,r_p

      write(*,*)'antenna_vertex before x_bin_min p1,p2,pacc',
     &p1,p2,pacc

      p_antenna= x_bin_min(f_antenna_min,p1,p2,pacc)

c      write(*,*)'after x_bin_min p_antenna',p_antenna
      
      x_st=r_0*dcos(phi_0)+p_antenna*n_vac_x
      y_st=r_0*dsin(phi_0)+p_antenna*n_vac_y 
      z_st=z_0+p_antenna*n_vac_z
      r_st=dsqrt(x_st**2+y_st**2)

      phi_st=dacos(x_st/r_st)
      if (y_st.lt.0.d0) phi_st=-phi_st

c-----calculation the central ray direction

c--------poloidal angle beta_st calculation
      pi=4.d0*datan(1.d0)
      if(z_0.le.z_st) then
         if(r_0.le.r_st) then
            !  0=<beta_st<pi/2
            beta_st=dasin(n_vac_z)
         else  
            ! pi/2<beta_st<pi
            beta_st=pi-dasin(n_vac_z)
         endif
      else
         if(r_0.le.r_st) then
            !  -pi/2<beta_st=<0
            beta_st=dasin(n_vac_z)
         else  
            ! -pi<beta_st<-pi/2
            beta_st=-pi-dasin(n_vac_z)
         endif
      endif
      beta_st=-beta_st

c-------toroidal angle calculation alph-a_st
      delta=dsqrt(r_st**2+r_0**2-2.d0*r_st*r_0*dcos(phi_st-phi_0))
      gamma=dacos((r_st**2+delta**2-r_0**2)/(2.d0*r_st*delta))

      if (phi_st.ge.phi_0) then
        alpha_st=pi+gamma
      else
        alpha_st=pi-gamma
      endif 

c-----transformation from radians to degrees
c      write(*,*)'alpha_st,phi_st',alpha_st,phi_st
c      write(*,*)'alpha_st+phi_st',alpha_st+phi_st
      trans=180.d0/pi
      beta_st=  beta_st*trans
      alpha_st= alpha_st*trans
      phi_st=   phi_st*trans
c      write(*,*)'alpha_st,phi_st degree',alpha_st,phi_st
c      write(*,*)'alpha_st+phi_st degree',alpha_st+phi_st
      return
      end


c======================================================================
c======================================================================



      subroutine set_common_ox(r_0_l,phi_0_l,z_0_l,
     &n_vac_x_l,n_vac_y_l,n_vac_z_l)
c-----put the variables to common/ox/
      implicit none
c-----input
      real*8 r_0_l,phi_0_l,z_0_l,
     &n_vac_x_l,n_vac_y_l,n_vac_z_l

      real*8 r_0,phi_0,z_0,
     &n_vac_x,n_vac_y,n_vac_z
      common /ox/r_0,phi_0,z_0,
     &n_vac_x,n_vac_y,n_vac_z

      r_0=r_0_l
      phi_0=phi_0_l
      z_0=z_0_l
      n_vac_x= n_vac_x_l
      n_vac_y= n_vac_y_l
      n_vac_z= n_vac_z_l

c      write(*,*)'in set_common_ox'
c      write(*,*)'r_0,phi_0,z_0,n_vac_x,n_vac_y,n_vac_z',
c     &r_0,phi_0,z_0,n_vac_x,n_vac_y,n_vac_z

      return
      end


c======================================================================
c======================================================================


      subroutine ox_conversion(r_in,z_in,phi_in,nr_in,nz_in,m_in,
     &rma,zma, !temporally
     &eps_xe,  ! The parameter which sets the vicinitity
               ! of the O-mode cutoff surface
               ! If xe >(1-eps_xe) then this subroutine  
               ! makes the ray jump in small radius direction
               ! and finds X mode.
c-------the point is near the OX mode conversion area
     &r_x,z_x,phi_x,nr_x,nz_x,m_x,i_ox_conversion)
c     Creates the jump of the ray point throw the OX mode conversion
c     area where X_e=1 (V_perp=0)

      implicit none               

c-----input
      real*8 r_in,z_in,phi_in,nr_in,nz_in,m_in !the cordinates before
                                               !O mode cutoff
      real*8 eps_xe

c-----output
      real*8 r_x,z_x,phi_x,nr_x,nz_x,m_x !the coordinates of
                                         !X mode cutoff point
      integer i_ox_conversion !=0 was not OX conversion
                              !=1 was the jump in the radial direction
c-----locals
      real*8 xe_0,xe_in,bmod,
     &nphi_in,nx_in,ny_in,ntheta,nphi,
     &costheta,sintheta,theta,pi,
     &rma,zma,
     &cnpar,
     &cnpar_x,cnper_x,rho_x,m_in_loc,
     &transm_ox_loc !for test

      integer iraystop

c-----externals
      real*8 x,b

      pi=4.d0*datan(1.d0)
 
      i_ox_conversion=0
      xe_0=1.d0

      bmod=b(z_in,r_in,phi_in)
      xe_in=x(r_in,z_in,phi_in,1)

      write(*,*)'before find_maximal_OX transmission'
      write(*,*)'r_in,z_in,phi_in',r_in,z_in,phi_in
c      write(*,*)'oxb r_in,z_i,bmod',r_in,z_in,bmod
c      write(*,*)'oxb xe_0,xe_in',xe_0,xe_in 
c      write(*,*)'before find_maximal_OX_transmission'
       write(*,*)'ox_conversion 0 ph_in',phi_in

      call find_maximal_OX_transmission(eps_xe,xe_0,
     & r_in,z_in,phi_in,nr_in,nz_in,m_in,
     & i_ox_conversion)
      
cyup      write(*,*)'after find_maximal_OX_transmission
cyup     & i_ox_conversion',i_ox_conversion
      write(*,*)'r_in,z_in,phi_in',r_in,z_in,phi_in
      write(*,*)'ox_conversion 1 ph_in',phi_in
      call transmit_coef_ox(z_in,r_in,phi_in,nz_in,nr_in,m_in,
     &                           transm_ox_loc)
      write(*,*)'transm_ox_loc',transm_ox_loc
      write(*,*)'ox_conversion 2 ph_in',phi_in
cendtest

      if(i_ox_conversion.eq.1) then
c-------the point is near the OX mode conversion area
        nphi_in=m_in/r_in

c-------trasnform the vector coordinates Nr_in,Nphi_in to Nx_in,Ny_in
        write(*,*)'oxb.f before transform_N_rphiz_to_xyz nr_in,nphi_in'
     & ,nr_in,nphi_in
        call transform_N_rphiz_to_xyz(nr_in,nphi_in,phi_in,nx_in,ny_in)
        write(*,*)'ox_conversion 3 ph_in',phi_in
c----------------------------------------------------------------
c       ninit_ec creates the tangent to magnetic surface
c       components  of the refractive index cnteta,cnphi                 
c       in the initial point (z_in,r_in,phi_in) for ECR wave
c-----------------------------------------------------------------
        write(*,*)'ox_conversion 4 ph_in',phi_in
        call ninit_ec(z_in,r_in,phi_in,nx_in,ny_in,nz_in,ntheta,nphi)
        write(*,*)'ox_conversion 5 ph_in',phi_in

        costheta=(r_in-rma)/dsqrt((z_in-zma)**2+(r_in-rma)**2)
        sintheta=(z_in-zma)/dsqrt((z_in-zma)**2+(r_in-rma)**2)

        write(*,*)'in ox_conversion sintheta,costheta',sintheta,costheta

        if (costheta.gt.1.d0)  costheta= 1.d0
        if (costheta.lt.-1.d0) costheta=-1.d0
        if (sintheta.ge.0.d0) then
           theta=+dacos(costheta)
        else  
           theta=-dacos(costheta)
        endif
        write(*,*)'**************************************'
        write(*,*)'r_in,z_in',r_in,z_in
        write(*,*)'in ox_conversion theta[degree]',theta*180.d0/pi

c--------------------------------------------------------------------
c       cninit solves the dispersion relation N=N(n_par)
c       Then subroutine calculates the initial components
c       of the refractive index  cnz,cnr,cm
c---------------------------------------------------------
c
c        call cninit(z_in,r_in,phi_in,cnpar,ntheta,nphi,
c     &              z_out,r_out,m_out,iraystop)
c--------------------------------------------------------
c     find the small radius rho_x at the vector rho^ where
c     X mode with ntheta, nphi has the cutoff.
c     The vector rho^ starting at the magnetic axis O(rma,zma).
c     The poloidal angle of this vector is theta (radians).
c-----------------------------------------------------------

       write(*,*)'oxb. in  ox_conversion before call find_rho_x'
       write(*,*)'ox_converion before find_rho z_in,r_in,m_in',
     & z_in,r_in,m_in
 
      m_in_loc= m_in
      call find_rho_x(theta,phi_in,
     &  z_in,r_in,ntheta,nphi,m_in_loc,
     &  rho_x,z_x,r_x,nz_x,nr_x,m_x,cnpar_x,cnper_x)
       bmod=b(z_x,r_x,phi_in)
      endif
      phi_x=phi_in
      return
      end
      

c======================================================================
c======================================================================

         
      subroutine find_rho_x(theta,phi,
     &z_ini,r_ini,cntheta_ini,cnphi_ini,cm_ini,
     &rho_x,z_x,r_x,cnz_x,cnr_x,cm_x,cnpar_x,cnper_x)

c-----finds the small radius rho_x at the vector rho^ where
c     X mode with ntheta, nphi has the cutoff.
c     The vector rho^ starting at the magnetic axis O(rma,zma).
c     The poloidal angle of this vector is theta (radians).

      implicit none
      include 'param.i'
      include 'one.i'
c-----input
      real*8 theta              ! poloidal angle(radians) measured from the
                                ! midplane in clockwise direction 
      real*8 phi                ! toroidal angle (radians) 
      real*8 z_ini,r_ini        ! z,r coordinates before O cutoff
c      real*8 rho_ini            ! the initial radius before O cutoff
      real*8 cntheta_ini        ! N poloidal in rho_x point
      real*8 cnphi_ini          ! N_toroidal in rho_x point
      real*8 cm_ini             ! N_phi*r_ini 
c     real*8 rho is the small radius in (z_ini,r_ini) point .
c     It is in common/one/ in one.i
c-----output
      real*8 rho_x,             ! normalized small radius in X cutoff point
     *r_x,z_x
      real*8 cnpar_x,cnper_x    ! N_parallel N_perp in X_cutoff point
      real*8 cnz_x,cnr_x,cm_x   ! refractive index coordinates of
                                ! X cutoff point      
c-----external zr_psith
      real*8 psi_rho_loc,b
      real*8 psi_rho,rhopsi,x
c
c-----local
      real*8 psi,hstep,ppp,gradpsi
      integer iraystop,id_loc,ioxm_loc
      real*8 z,r,cnz,cnr,cnpar
      real*8 rho_loc,cnper2p,cnper2m,cnperp_loc,cnperm_loc,cnper2,cn2,
     &cntang2,cnrho2,cnrho,cnperp,cm,cnpar_loc,cn_loc,cnper_loc,
     &xe

c      save was_not_o_cutoff
      logical was_not_o_cutoff
c      data was_not_o_cutoff /.true./

      hstep=1.d-3
      hstep=1.d-5
      pi=4.d0*datan(1.d0)
      
      stop 'find_rho_x: not used?'
       
      rho_loc=rho !initialization
      write(*,*)'in find_rho_x rho_loc initial',rho_loc
      write(*,*)'in find_rho_x phi',phi
      was_not_o_cutoff= .true.
      id_loc=id
      ioxm_loc=ioxm 
      
 10   continue
      write(*,*)'oxb.f in  find_rho_x rho_loc',rho_loc
      psi=psi_rho(rho_loc)
           
      call zr_psith(psi,theta,z,r) !-> z(psi,theta),r(psi,theta)

      write(*,*)'oxb.f in  find_rho_x after zr_psith'     
      write(*,*)'psi,theta,z,r',psi,theta,z,r

c---------------------------------------------------------------
c     calculations of the parallel (to the magnetic field)
c     refractive index component cnpar
c     N_par=(B_phi*N_phi+N_theta*e_theta*{e_z*B_z+e_r*B_r})/bmod
c     e_theta=(e_z*dpsi/dr-e_r*dpsi/dz)/abs(grad(psi))
c--------------------------------------------------------------
      bmod=b(z,r,phi)
      ppp=(dpdzd*dpdzd+dpdrd*dpdrd)
      gradpsi=dsqrt(dpdzd*dpdzd+dpdrd*dpdrd)
      cnpar=(cnphi_ini*bphi+cntheta_ini*(bz*dpdrd-br*dpdzd)/gradpsi)/
     &bmod

      id=2
      
c      write(*,*)'find_rho_x was_not_o_cutoff= ',was_not_o_cutoff 
c      write(*,*)'id_loc,id',id_loc,id

      if (was_not_o_cutoff) then
         ioxm=1
      else
         ioxm=-1
      endif

      call npernpar(z,r,phi,cnpar,cnper2p,cnper2m)
c      write(*,*)'cnper2p,cnper2m',cnper2p,cnper2m
      if(cnper2p.ge.0.d0) then
          cnperp_loc=dsqrt(cnper2p)
c          write(*,*)'cnperp_loc=',cnperp_loc
      endif
      if(cnper2m.ge.0.d0) then
          cnperm_loc=dsqrt(cnper2m)
c          write(*,*)'cnperm_loc=',cnperm_loc
      endif

      if (.not.was_not_o_cutoff) then
        if ((cnper2p.le.0.d0).or.(cnper2m.lt.0.d0)) then
           rho_loc=rho_loc-hstep
           goto 10
        else           
           if (cnper2p.lt.cnper2m) then
             cnper2=cnper2p
           else
             cnper2=cnper2m
           endif 

           cn2=cnper2+cnpar**2
           cntang2=cntheta_ini**2+cnphi_ini**2

           if(cn2.lt.cntang2)then
             rho_loc=rho_loc-hstep
             goto 10
           endif

           cn2=cnpar**2+cnper2
           cnrho2=cn2-cntheta_ini**2-cnphi_ini**2
           cnrho=dsqrt(cnrho2)
           cnperp=cnrho
           write(*,*)'cn2',cn2
           call cnzcnr(z,r,phi,cntheta_ini,cnphi_ini,cnrho,cnz,cnr,cm)

           iraystop=0 
           write(*,*)'oxb.f find_rho_x bef goto30 id_loc,id',id_loc,id
           goto 30
        endif
      endif

      call cninit12(z,r,phi,cnpar,cntheta_ini,cnphi_ini,
     &cnz,cnr,cm_ini,iraystop)

 30   continue

      id=id_loc
      ioxm=ioxm_loc

      if (was_not_o_cutoff) then
c        write(*,*)'iraystop',iraystop
        if(iraystop.eq.0) then
c---------ioxm=1 O-mode exists in this point
          rho_loc=rho_loc-hstep
          if (rho_loc.lt.0d0) then 
             write(*,*)'rho_x was not found the rho_loc point'
             stop
          endif
          go to 10                   
        else 
         was_not_o_cutoff=.false.
         write(*,*)'was_not_o_cutoff',was_not_o_cutoff

 20      continue
         xe=x(z,r,phi,1)
         if (xe.le.1) then
           rho_loc=rho_loc-hstep
           psi=psi_rho(rho_loc)      
           call zr_psith(psi,theta,z,r) !-> z(psi,theta),r(psi,theta)   
           bmod=b(z,r,phi)
           goto 20
        endif
  
         go to 10 
        endif
      else
        if(iraystop.eq.1) then
c---------ioxm=1 O-mode does not exist in this point
          rho_loc=rho_loc-hstep
          if (rho_loc.lt.0d0) then 
             write(*,*)'rho_x was not found the rho_loc point'
             stop
          endif
          go to 10   
        else
          rho_x=rho_loc
          r_x=r
          z_x=z
c          write(*,*)'r_x,z_x',r_x,z_x
          cnr_x=cnr
          cnz_x=cnz
          cm_x=cm_ini
          cnpar_x=cnpar
          cnper_x=dsqrt(cnz_x**2+cnr_x**2+(cm_x/r_x)**2-cnpar_x**2)
c          write(*,*)'in rho_x rho_x,cnper_x',rho_x,cnper_x
          id=id_loc
c          write(*,*)'id=',id
        endif
      endif

      write(*,*)'phi',phi
      return
      end


c======================================================================
c======================================================================



      subroutine ox_conversion_grill_in_poloidal_point(theta_pol,
     &i_n_optimal)

c    It calculates: 
c	  1) the space point M_conv=(rhoconv,theta)=(zconv,rconv)
c            at which x_e=x0=1. The value of the poloidal
c            angle theta (degree)=theta_pol is given bellow
c            if theta=0  point is on the outer part of the magnetic surface
c            if theta=90 point is on the top of the magnetic surface
c         2) the optimal value N_parallel_optimal=cnparopt
c            for O_X mode conversion
c            if i_n_optimal=1 cnparopt >0
c            if i_n_optimal=2 cnparopt <0
c            
c         3) gives
c            anmin(1)=cnparopt-0.01d0 
c            anmax(1)=cnparopt+0.01d0 
c         4) call grill_lh to calculate the initial data
c         5) gives the coordinates of the initial O_X mode conversion poin 
c            arzu0(1)=zconv
c            arru0(1)=rconv
c 
c     output data are in 'grill.i'
c-------------------------------------
c     that gives the optimal N_paralle at X_e=1 surface 
     
      implicit integer (i-n), real*8 (a-h,o-z)

c-----input 
      real*8 theta_pol
      integer i_n_optimal
     
      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'three.i'
      include 'five.i'
      include 'cone.i'
      include 'grill.i'
      include 'rkutta.i'
      include 'six.i'
      !!! include 'scatnper.i'
      include 'write.i'

      istart=3
      
c-----ECR O_X_EBW mode conversion case----------
c         It uses the wave input data from the grill form
c         It sets i_n_poloidal=1
c	  It calculates: 
c	  1) the space point M_conv=(rhoconv,theta)=(zconv,rconv)
c            at which x_e=x0=1. The value of the poloidal
c            angle theta (degree) is given bellow
c            if theta=0  point is on the outer part of the magnetic surface
c            if theta=90 point is on the top of the magnetic surface
c         2) the optimal value N_parallel_optimal=cnparopt
c            for O_X mode conversion
c            if i_n_optimal=1 cnparopt >0
c            if i_n_optimal=2 cnparopt <0
c            
c         3) gives
c            anmin(1)=cnparopt-0.01d0 
c            anmax(1)=cnparopt+0.01d0 
c         4) call grill_lh to calculate the initial data
c         5) gives the coordinates of the initial O_X mode conversion poin 
c            arzu0(1)=zconv
c            arru0(1)=rconv
c-------------------------------------------------------------
      i_n_poloidal=1
c     the following lines are to create the start point inside the 
c     plasma for the given point xe=x0
      x0=1.d0
c      x0=0.998d0
      x0=0.999d0
c      x0=0.98d0

c      theta=0.d0   !poloidal angle  (degree)
c      theta=-30.d0
c      theta=thgrill(1)

      theta=theta_pol
      thgrill(1)=theta

     
      write(*,*)'ox_central_ray_direction before owconvr theta='
     &,theta

      call owconvr (theta,x0,rhoconv,zconv,rconv)

      write(*,*)'ox_central_ray_direction rhoconv,zconv,rconv',
     &           rhoconv,zconv,rconv

      rhopsi0(1)=rhoconv
      phiconv=0.d0
      bmod=b(zconv,rconv,phiconv)
      xconv=x(zconv,rconv,0.d0,1)
      yconv=y(zconv,rconv,0.d0,1)
c-----calculation of the optimal value N_parallel_optimal
c     for O_X mode conversion
      cnparopt=dsqrt(yconv/(1.d0+yconv))
cSm050529
      if  (i_n_optimal.eq.2) cnparopt=-cnparopt

      write(*,*)'ox_central_ray_direction xconv,yconv,cnparopt',
     &           xconv,yconv,cnparopt

c      write(*,*)'ox_central_ray_direction old value of rhopsi0(1)',
c     &rhopsi0(1) 

      rhopsi0(1)=rhoconv

      write(*,*)'ox_central_ray_direction new rhopsi0(1)',rhopsi0(1) 
      write(*,*)'ox_central_ray_direction anmin(1),anmax(1)',
     &anmin(1),anmax(1)

      anmin(1)=cnparopt-0.01d0 
      anmax(1)=cnparopt+0.01d0
 
      write(*,*)' ox_central_ray_direction new anmin(1),anmax(1)',
     &anmin(1),anmax(1)

c---------------------------------------------------------------
      write(*,*)' ox_central_ray_direction before grilk_lh
     1ngrilla,ngrill',
     1ngrilla,ngrill

      call grill_lh(rhopsi0,ngrilla,ngrill,thgrill,phigrill,
     1height,nthin,nthinmax,
     1anmin,anmax,nnkpar,powers,powtott,
     1antormin,antormax,nnktor,anpolmin,anpolmax,nnkpol,
     1n_theta_pol,
     1rma,zma,psimag,
cSm050309
     1fpwth,wtheta,wtheta_edge,anzin,anzin_edge,pwcpl,nnkprmax,
     1anztorin,anzpolin,pwcpl_tp,nnktormax,nnkpolmax,
     1anztorin_edge,anzpolin_edge,pwcpl_t_ar,pwcpl_p_ar,
     1nray,arzu0,arru0,arphiu0,arntheta,arnphi,powinilh,
     1nraymaxl,wdnpar0,igrillpw,igrilltw,i_n_poloidal,
     1ilaunch,r0launch,phi0launch,z0launch,
cSm050309
     1i_grill_pol_mesh,i_grill_npar_ntor_npol_mesh)

c for ecr internal case 
      arzu0(1)=zconv
      arru0(1)=rconv
c--end ecs internal case

      return
      end
      

c======================================================================
c======================================================================


      real*8 function transmission_ox(n_par,n_pol,Y,L_n,freqncy)
c-----calculate the transmission function for OX mode conversion
c     using Preinhaelter and Kopecky formula (1973)
         
      implicit none

c-----input
      real*8 n_par,   !parallel to magnetic field refractive index
     &n_pol,          !poloidal refactive index 
     &Y,              !omega_ce/omega
     &L_n,            !electron density scale length [cm]
     &freqncy         !wave frequency [GHZ]
c-----local
      real*8 pi,n_par_optimal,omega,clight,Y_abs
      real*8 transmission_ox_p,transmission_ox_m
      pi=4.d0*datan(1.d0)

      omega=2.d0*pi*freqncy*1.d9

      clight=2.9979d10 !cm/sec
      
      Y_abs=dabs(Y)
      n_par_optimal=dsqrt(Y_abs/(Y_abs+1))   
c      write(*,*)'oxb.f in transmission_ox L_n,Y,freqncy',L_n,Y,freqncy
c      write(*,*)'oxb.f in transmission_ox n_par_optimal,n_par,n_pol',
c     & n_par_optimal,n_par,n_pol
      transmission_ox_p=dexp(-pi*omega*L_n/clight*dsqrt(0.5d0*Y_abs)*
     &        (2.d0*(1.d0+Y_abs)*(n_par_optimal-n_par)**2+n_pol**2))
     
      transmission_ox_m=dexp(-pi*omega*L_n/clight*dsqrt(0.5d0*Y_abs)*
     &        (2.d0*(1.d0+Y_abs)*(-n_par_optimal-n_par)**2+n_pol**2))

cSm050529
c      transmission_ox=transmission_ox_p+transmission_ox_m
cyup      write(*,*)'oxb.f in transmission_ox L_n,n_par_optimal,n_par',
cyup     & L_n,n_par_optimal,n_par

      transmission_ox=dmax1(transmission_ox_p,transmission_ox_m)
      
cyup      write(*,*)'transmission_ox_p,transmission_ox_m,transmission_ox',
cyup     &transmission_ox_p,transmission_ox_m,transmission_ox

      return
      end


c======================================================================
c======================================================================


      subroutine transmit_coef_ox(z,r,phi,cnz,cnr,cm,transm_ox)
c     calculate transmission coefficient for OX mode conversion 

c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none
      include 'param.i' 
      include 'one.i'
c-----input
      real*8 z,r,phi, !the coordinates [m,m,radians]
     &cnr,cnz,cm      !refractive index coordinates

c-----output
      real*8 transm_ox ! transmission coefficient of OX mode conversion

c-----locals
      real*8 y_e,L_n,grad_x,cnpar,cnpol,b_pol,freqncy
   
c-----externals
      real*8 b,y,x,dxdr,dxdz,dxdphi,transmission_ox 
     
c      write(*,*)'transmit_coef_ox z,r,phi,cnz,cnr,cm',
c     &z,r,phi,cnz,cnr,cm

      bmod=b(z,r,phi)
c      write(*,*)'transmit_coef_ox bz,br,bphi,bmod',bz,br,bphi,bmod
      y_e=y(z,r,phi,1)
      
c-----L_n=X/grad_X
     
      grad_x=dsqrt(dxdr(z,r,phi,1)**2+dxdz(z,r,phi,1)**2+
     &             (dxdphi(z,r,phi,1)/r)**2)

c      write(*,*)'oxb.f dxdr(z,r,phi,1)',dxdr(z,r,phi,1)
c      write(*,*)'oxb.f dxdz(z,r,phi,1)',dxdz(z,r,phi,1)      
c      write(*,*)'oxb.f dxdphi(z,r,phi,1)',dxdphi(z,r,phi,1)

      L_n=x(z,r,phi,1)/grad_x*100.d0        !cm
c      write(*,*)'oxb.f z,r,phi,grad_x',z,r,phi,grad_x
c      write(*,*)'oxb.f x(z,r,phi,1)',x(z,r,phi,1)
c      write(*,*)'oxb.f L_n',L_n

      cnpar=(cnz*bz+cnr*br+cm*bphi/r)/bmod

c      write(*,*)'transmit_coef_ox cnpar',cnpar
c      write(*,*)'cnz,bz,cnr,br,cm,bphi,r,bmod',
c     &cnz,bz,cnr,br,cm,bphi,r,bmod

      freqncy=frqncy                         ! GHZ
  
      b_pol=dsqrt(bz**2+br**2)
      cnpol=(cnz*bz+cnr*br)/b_pol     !poloidal refractive index

c-----Instead of N_poloidal we will use the refractive index
c     N^*[b*gradPsi]/|[b*gradPsi]|
c     [b*gradPsi] is a vector production
c
c     b^=e^r*br+e^phi*bphi+e^z*bz
c     gradPsi^=e^r*dpdrd+e^phi*0+e^z*dpdzd
c
c                 | e^r    e^phi   e^z  |
c     [b*gradPsi]=| br     bphi    bz   |=
c                 | dpdrd  0       dpdzd|
c     =e^r*bphi*dpdzd - e^phi(br*dpdzd-bz*dpdrd) - e^z*bphi*dpdrd
c
c     N_[b*gradPsi]=N^*[b*gradPsi]/|[b*gradPsi]|

      cnpol=(cnr*bphi*dpdzd-(cm/r)*(br*dpdzd-bz*dpdrd)-cnz*bphi*dpdrd)/
     &dsqrt((bphi*dpdzd)**2+(br*dpdzd-bz*dpdrd)**2+(bphi*dpdrd)**2)

      transm_ox=transmission_ox(cnpar,cnpol,y_e,L_n,freqncy)
cyup      write(*,*)'oxb.f in transmit_coef_ox transm_ox=',transm_ox
      return
      end
      

c======================================================================
c======================================================================


      subroutine OX_power_transmission(is,i_ox,i_ox_conversion,
     & delpwr_o,transm_ox,delpwr_x)
c----------------------------------------------------------
c     If i_ox.ne.2 then it does nothing
c     delpwr_X=delpwr_X
c
c     It works for i_ox.eq.2 case.
c     It calculates the power in the ray chanel 
c     after OX conversion jump 
c     delpwr_X=delpwr_O*transm_ox
c     Here:
c     delpwr_O is the power in the ray chanel 
c              in the first point after after OX conversion
c
c     transm_ox is OX transmission coefficient
c
c     delpwr_X is the power of X-mode in the ray
c              chanel after OX conversion
c
c----------------------------------------------------------
      implicit none

c-----input:
      real*8 delpwr_o,! power before trasmission
     &transm_ox 

      integer is,i_ox,
     &i_ox_conversion !=0 was not OX conversion
                      !=1 was OX conversion 
                      !calculated in output.f 
                      !in call ox_conversion()
c-----output
      real*8 delpwr_x   !after transmission

c-----local
      logical first
c      data first /.true./
      save first
     
      if(is.eq.1) first=.true.

cyup      write(*,*)'OX_power_transmission is,i_ox,first,i_ox_conversion',
cyup     &is,i_ox,first,i_ox_conversion
     
      if(i_ox.ne.2) then
        delpwr_x=delpwr_o
        return
      else 
        if(i_ox_conversion.eq.1) then
           if (first) then
             delpwr_x=delpwr_o*transm_ox
             first=.false.
           else
             delpwr_x=delpwr_o
           endif
        else
          delpwr_x=delpwr_o
        endif
      endif

      return
      end


c======================================================================
c======================================================================



      subroutine write_optimal_O_cone_central_ray_multivertex(icone,
     & i_n_optimal)
c-----write the coordinates of EC central ray
c     with the optimal direction for OX mode conversion

c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none

      include 'param.i'
      include 'cone.i'

c-----input
      integer icone ! EC cone vetrex number 
      integer i_n_optimal !1 or 2
                          ! if i_n_optimal=1 N_par_optimal=N_par_optimal
                          ! if i_n_optimal=2 N_par_optimal=-N_par_optimal
c-----local
      character*40 ch_zst,ch_rst,ch_phist,ch_betast,ch_alfast
      character*40 ch_i_n_optimal   
      character*3  text
     
      character*6 format,format_i
      logical first
      data first /.true./
      save first

      write(*,*)'in write_optimal_O_cone_central_ray_multivertex'
            

      if (first) then 
         open(50,file='ECcone_optimal.dat')
         first=.false.
      endif

      write(*,*)'icone',icone,'i_n_optimal',i_n_optimal
      
cSAP0890423 It was a problem under gfortran here
c      write(text,*)icone
c      write(*,*)'text ',text
      format='d21.15'
       
      format_i='i2'
cSAP080423
c      ch_i_n_optimal='(1X,"i_n_optimal('//text//')=",'//format_i//')'
c      ch_zst='(1X,"zst('//text//')=",'//format//')'
c      ch_rst='(1X,"rst('//text//')=",'//format//')'
c      ch_phist='(1X,"phist('//text//')=",'//format//')'
c      ch_betast='(1X,"betast('//text//')=",'//format//')'
c      ch_alfast='(1X,"alfast('//text//')=",'//format//')'

      ch_i_n_optimal=
     &'(1X,"i_n_optimal("'//format_i//'")=",'//format_i//')'
      ch_zst=
     &'(1X,"zst("'//format_i//'")=",'//format//')'
      ch_rst=
     &'(1X,"rst("'//format_i//'")=",'//format//')'
      ch_phist=
     &'(1X,"phist("'//format_i//'")=",'//format//')'
      ch_betast=
     &'(1X,"betast("'//format_i//'")=",'//format//')'
      ch_alfast=
     &'(1X,"alfast("'//format_i//'")=",'//format//')'
      
      write(*,*)z_st_ox
      write(*,*)r_st_ox
      write(*,*)phi_st_ox
      write(*,*)beta_st_ox
      write(*,*)alpha_st_ox
cSAP080423
c     write(50,fmt=ch_i_n_optimal)i_n_optimal
c      write(50,fmt=ch_zst)z_st_ox
c      write(50,fmt=ch_rst)r_st_ox
c      write(50,fmt=ch_phist)phi_st_ox

      write(50,fmt=ch_i_n_optimal)icone,i_n_optimal
      write(50,fmt=ch_zst)icone,z_st_ox
      write(50,fmt=ch_rst)icone,r_st_ox
      write(50,fmt=ch_phist)icone,phi_st_ox

c      if (raypatt.eq.'genray') then
cSAP080423
c         write(50,fmt=ch_betast)beta_st_ox
c         write(50,fmt=ch_alfast)alpha_st_ox
         write(50,fmt=ch_betast)icone,beta_st_ox
         write(50,fmt=ch_alfast)icone,alpha_st_ox
c      endif

c      if (raypatt.eq.'toray') then
c         write(50,fmt=ch_betast)(90.d0-beta_st_ox)
c         write(50,fmt=ch_alfast)(90.d0-alpha_st_ox)
c      endif
           
      if ((i_n_optimal.eq.2).and.(icone.eq.ncone)) close(50)

      return
      end


c======================================================================
c======================================================================


      subroutine write_optimal_O_cone_central_ray
c-----write the coordinates of EC central ray
c     with the optimal direction for OX mode conversion

c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none

      include 'param.i'
      include 'cone.i'

      write(*,*)'in write_optimal_O_cone_central_ray'
      write(*,1)z_st_ox
      write(*,2)r_st_ox
      write(*,3)phi_st_ox
      write(*,4)beta_st_ox
      write(*,5)alpha_st_ox

      open(50,file='ECcone_optimal.dat')
      write(50,1)z_st_ox
      write(50,2)r_st_ox
      write(50,3)phi_st_ox
      write(50,4)beta_st_ox
      write(50,5)alpha_st_ox

      close(50)

 1    format('zst(1)= ',d21.15)
 2    format('rst(1)= ',d21.15)
 3    format('phist(1)= ',d21.15)
 4    format('betast(1)= ',d21.15)
 5    format('alfast(1)= ',d21.15)

      return
      end



c======================================================================
c======================================================================


      subroutine gr_OX_optimal_direction(ndim)
c-----calculate the OX optimal direction of the EC cone central ray
      

      implicit integer (i-n), real*8 (a-h,o-z)
  
      include 'param.i'
      include 'one.i'
      include 'cone.i'
      include 'rkutta.i'
      include 'grill.i'
      include 'three.i' 
      dimension u(6),deru(6),aux(8,6)
      external rside1,b

c-----input integer ndim ! the number of ODE equations     

 
      write(*,*)'in gr_OX_optimal_direction ncone',ncone

      do icone=1,ncone

         icone_loc=icone ! for cone.i file
         
         do i_n_optimal=1,2 ! the loop for pos. and neg. N_parallel_optimal

           i_ox_poloidal=0 ! the counter of the poloidal rays for OX
                           ! optimal direction calculations
              
40         continue ! begin of the loop for the calculations 
                    ! of the optimal direction in the given EC cone vertex
                    ! for OX mode conversion

           if(i_ox.eq.1) then
             i_ox_poloidal=i_ox_poloidal+1    
             write(*,*)'genray.f oxb  i_ox_poloidal ',i_ox_poloidal

             if(i_ox_poloidal.eq.1) then
               theta_pol_0=theta_bot(icone) !theta_bot is given 
                                            !in genray.dat file 
                                           !It is the boundary of the poloidal 
                                           !angle at Xe=1 surface
             else
               if(i_ox_poloidal.eq.2)then
                 theta_pol_0=theta_top(icone) !theta_top is given 
                                              !in genray.dat file
                                            !It is the boundary of the poloidal
                                            !angle at Xe=1 surface
               else
                 theta_pol_0=0.5d0*(theta_pol_left+theta_pol_right)
               endif
             endif

             theta_pol=theta_pol_0              
c--------------------------------------------------------
c            for the given poloidal angle theta_pol
c            it will calculate
c            the coordinates of the initial O_X mode conversion point 
c            the optimal value N_parallel_optimal=(+,-) cnparopt
c            If i_n_optimal.eq.1 then  N_parallel_optimal=+cnparopt
c            If i_n_optimal.eq.2 then  N_parallel_optimal=-cnparopt
c            set the n_parallel grill spectrum
c            anmin(1)=cnparopt-0.01d0 
c            anmax(1)=cnparopt+0.01d0 
c--------------------------------------------------------
             write(*,*)'genray.f oxbv'
             write(*,*)'before ox_conversion_grill_in_poloidal_point'
             write(*,*)'theta_pol ',theta_pol
             call ox_conversion_grill_in_poloidal_point(theta_pol,
     &       i_n_optimal)
           endif !i_ox.eq.1
c-----------------------------------------------------------------1end

           iray=1

           write(*,*)'oxb.f sub gr_OX_optimal_direction bef dinit_1ray'

           call dinit_1ray(arzu0(iray),arru0(iray),arphiu0(iray),
     1	   alfast1,betast1,arntheta(iray),arnphi(iray),u,iraystop)
           write(*,*)'!!!oxb  dinit_1ray'	  
          
	   if(iraystop.eq.1) then
	      write(*,*)'icone=',icone,'iray=',iray,'iraystop=1
     &                bad initial conditions'
              nrayelt=0
	      goto 20
	   endif 
  
	   prmt(7)=prmt(1)+prmt(6)
c----------------------------------------------------------
c          call b() to calculate the small radius rho (inside b())
c          the result will be in common block  one.i
c-----------------------------------------------------------
           bmod=b(u(1),u(2),u(3))

CMPIINSERTPOSITION STARTRUNGEKUTTA

           nstep_rk=1 ! initialize the number of Runge-Kutta time step 
            
c--------------------------------------------------------------beg2
           if(isolv.eq.1) then
c--------------------------------------------------------------
c            The Runge-Kutta solution of 6 hamiltonian equations
c            with correction which gives 
c            hamiltonian conservation with
c            accuracy epscor=prmt(4)
c--------------------------------------------------------------
             if(irkmeth.eq.0) then
c              4_th order Runge-Kutta method with constant time step
ccc               call drkgs(prmt,u,deru,ndim,ihlf,rside1,outpt,aux)
                 stop 'set irkmeth=2; other values not supported'
             endif
             if(irkmeth.eq.1) then
c              5-th order Runge-Kutta method with variable time step,
ccc               call drkgs1(prmt,u,deru,ndim,ihlf,rside1,outpt,aux)
                 stop 'set irkmeth=2; other values not supported'
             endif
             if(irkmeth.eq.2) then
c              4th order Runge-Kutta with variable time step,
c              time step can be reduce or enlarge
               write(*,*)'prmt',prmt
               call drkgs2(prmt,u,deru,ndim,ihlf,rside1,outpt,aux,
     &                     i_output)
             endif
           endif ! isolv.eq.1
c-------------------------------------------------------------------end2


           if(isolv.eq.2) then
c-------------------------------------------------------------------
c            The Runge-Kutta solution of 5 hamiltonian equations
c            and the solution of the dispersion relation for the sixth
c            ray variable'
c-------------------------------------------------------------------------
ccc              call rkb1(prmt,u,deru,ndim,ihlf,rsideb1,outptb1,aux)
c------------------------------------------------------------------------
                 stop 'set isolv=1; other values not supported'
           endif ! isov.eq.2

c---------------------------------------------------------------------beg5
           if(i_ox.eq.1) then
c------------set the boundaries for the calculations of
c            OX optimal direction in EC cone vertex           

c---------------------------------------------------------------------beg3
             if(i_ox_poloidal.eq.1) then 
c---------------calculate the poloidal angle theta_st of the given
c               cone vertex
cSm050525
                rho_st=dsqrt((rst(icone)-rma)**2+(zst(icone)-zma)**2)
                cos_theta_st=(rst(icone)-rma)/rho_st
                sin_theta_st=(zst(icone)-zma)/rho_st
                if (cos_theta_st .gt.  1.d0)  cos_theta_st= 1.d0
                if (cos_theta_st .lt. -1.d0)  cos_theta_st=-1.d0
                if (sin_theta_st.ge.0.d0) then
                   theta_st=+dacos(cos_theta_st)
                else  
                   theta_st=-dacos(cos_theta_st)
                endif
               
                write(*,*)'***************************************'
                write(*,*)'rst(icone),zst(icone)',rst(icone),zst(icone)
                write(*,*)'rma,zma',rma,zma
                write(*,*)'rho_st',rho_st
                write(*,*)'rst(icone)-rma,cos_theta_st',
     &                    rst(icone)-rma,cos_theta_st
                write(*,*)'zst(icone)-zma,sin_theta_st',
     &                     zst(icone)-zma,sin_theta_st
                write(*,*)'theta_st[degree]',theta_st*180.d0/pi
                write(*,*)'***************************************'

                rho_st_ox=dsqrt((r_st_ox-rma)**2+(z_st_ox-zma)**2)
                cos_theta=(r_st_ox-rma)/rho_st_ox
                sin_theta=(z_st_ox-zma)/rho_st_ox
                if (cos_theta .gt.  1.d0)  cos_theta= 1.d0
                if (cos_theta .lt. -1.d0)  cos_theta=-1.d0
                if (sin_theta.ge.0.d0) then
                  theta_st_ox=+dacos(cos_theta)
                else  
                  theta_st_ox=-dacos(cos_theta)
                endif
                
                theta_st_left=theta_st_ox
                f_left=theta_st_left-theta_st
                write(*,*)'theta_st_left[degree],theta_st[degree]',
     &          theta_st_left*180.d0/pi,theta_st*180.d0/pi
                write(*,*)'f_left[degree],rho_st',
     &                     f_left*180.d0/pi,rho_st

                r_st_left=r_st_ox
                theta_pol_left=theta_pol_0
c               f_left=r_st_left-rst(icone)
             
                if(dabs(f_left*rho_st*1.d+2).le.eps_antenna) then
                   theta_ox=theta_pol_0
c                  write(*,*)'r_st_left,rst(icone),eps_antenna,theta_ox'
c     &                      ,r_st_left,rst(icone),eps_antenna,theta_ox
                 write(*,*)'theta_st_left,theta_st,eps_antenna,theta_ox'
     &             ,theta_st_left*180.d0/pi,theta_st*180.d0/pi,
     &             eps_antenna,theta_ox*180.d0/pi

c                  call write_optimal_O_cone_central_ray
                call write_optimal_O_cone_central_ray_multivertex(icone,
     &             i_n_optimal)
                   goto 20
c                  stop 'genray,f'
                else
                   goto 40
                endif
             endif !i_ox_poloidal.eq.1 
c---------------------------------------------------------------------end3                 

c---------------------------------------------------------------------beg4
             if(i_ox_poloidal.eq.2) then 
cSm050525
                rho_st_ox=dsqrt((r_st_ox-rma)**2+(z_st_ox-zma)**2)
                cos_theta=(r_st_ox-rma)/rho_st_ox
                sin_theta=(z_st_ox-zma)/rho_st_ox
                if (cos_theta .gt.  1.d0)  cos_theta= 1.d0
                if (cos_theta .lt. -1.d0)  cos_theta=-1.d0

                if (sin_theta.ge.0.d0) then
                   theta_st_ox=+dacos(cos_theta)
                else  
                   theta_st_ox=-dacos(cos_theta)
                endif

                theta_st_right=theta_st_ox
                f_right=theta_st_right-theta_st
                write(*,*)'theta_st_right[degree],theta_st[degree]',
     &                     theta_st_right*180.d0/pi,theta_st*180.d0/pi
                write(*,*)'f_right[degree],rho_st',
     &          f_right*180.d0/pi,rho_st
                r_st_right=r_st_ox
                theta_pol_right=theta_pol_0
c               f_right=r_st_right-rst(icone)

                if(dabs(f_right*rho_st*1.d+2).le.eps_antenna) then
                  theta_ox=theta_pol_0
                  write(*,*)'r_st_right,rst(icone),eps_antenna,theta_ox'
     &                      ,r_st_right,rst(icone),eps_antenna,theta_ox
                write(*,*)'theta_st_right,theta_st,eps_antenna,theta_ox'
     &                    ,theta_st_right,theta_st,eps_antenna,theta_ox 
               call write_optimal_O_cone_central_ray_multivertex(icone,
     &            i_n_optimal)
                  goto 20
                  stop 'genray,f'
                else
                  if(f_right*f_left.lt.0.d0) then
                    goto 40
                  else
                    write(*,*)'i_ox_poloidal=2'
                    write(*,*)'f_right*f_left < 0.d0 no roots'
                    write(*,*)'Please change theta_bop or theta_top'
                    goto 20
c                   stop 'genray oxb.f'
                  endif
                endif

             else ! i_ox_poloidal.gt.2

                write(*,*)'i_ox_poloidal,r_st_left,r_st_right',
     &                i_ox_poloidal,r_st_left,r_st_right
                write(*,*)'i_ox_poloidal,theta_st_left,theta_st_right',
     &                i_ox_poloidal,theta_st_left,theta_st_right

                write(*,*)'rst(icone),r_st_ox',rst(icone),r_st_ox
                write(*,*)'zst(icone),z_st_ox',zst(icone),z_st_ox

                rho_st_ox=dsqrt((r_st_ox-rma)**2+(z_st_ox-zma)**2)
                cos_theta=(r_st_ox-rma)/rho_st_ox
                sin_theta=(z_st_ox-zma)/rho_st_ox
                if (cos_theta .gt.  1.d0)  cos_theta= 1.d0
                if (cos_theta .lt. -1.d0)  cos_theta=-1.d0

                if (sin_theta.ge.0.d0) then
                  theta_st_ox=+dacos(cos_theta)
                else  
                  theta_st_ox=-dacos(cos_theta)
                endif
                write(*,*)'theta_st,theta_st_ox',theta_st,theta_st_ox

                f_center=theta_st_ox-theta_st
                write(*,*)'f_center,rho_st',f_center,rho_st

                if(dabs(f_center*rho_st*1.d+2).lt.eps_antenna) then
                  theta_ox=theta_pol_0
                  write(*,*)'The solution is found'
                  write(*,*)'r_st_ox,rst(icone),eps_antenna,theta_ox'
     &                      ,r_st_ox,rst(icone),eps_antenna,theta_ox
                  write(*,*)'theta_st_ox,theta_st,eps_antenna,theta_ox'
     &                      ,theta_st_ox,theta_st,eps_antenna,theta_ox
               call write_optimal_O_cone_central_ray_multivertex(icone,
     &            i_n_optimal)
                  goto 20
                endif

                if(i_ox_poloidal.gt.i_ox_poloidal_max) then
                  theta_ox=theta_pol_0
                  write(*,*)'genray.f'
                  write(*,*)'i_ox_poloidal.gt.i_ox_poloidal_max'
                  write(*,*)'r_st_ox,rst(icone),eps_antenna,theta_ox',
     &                       r_st_ox,rst(icone),eps_antenna,theta_ox
                  write(*,*)'theta_st_ox,theta_st,eps_antenna,theta_ox',
     &                       theta_st_ox,theta_st,eps_antenna,theta_ox
                call write_optimal_O_cone_central_ray_multivertex(icone,
     &            i_n_optimal)
                  goto 20
                endif

                if(f_left*f_center.lt.0.d0)then
                  theta_pol_right=theta_pol_0
                  f_right=f_center
                  write(*,*)'f_left*f_center.lt.0.d0'
                  write(*,*)'theta_pol_left,theta_pol_right',
     &                       theta_pol_left,theta_pol_right
                  goto 40
                else
                  !f_right*f_center.lt.0.d0
                  theta_pol_left=theta_pol
                  f_left=f_center
                  write(*,*)'f_right*f_center.lt.0.d0'
                  write(*,*)'theta_pol_left,theta_pol_right',
     &                       theta_pol_left,theta_pol_right
                  goto 40
                endif

             endif  !i_ox_poloidal.eq.2
           endif !i_ox.eq.1

 20      continue  
         enddo ! i_n_optimal=1,2 
      enddo   ! icone 
      
      stop 'oxb'                    
      return
      end


c======================================================================
c======================================================================


      real*8 function psilim_psif_zp_rp(p)
c---------------------------------------------------------------
c     It is used to find the itersection of the sraigth line
c     with the plasma boundary. 
c     It calculates the difference: psilim-psi(z,r)
c     Here z=z(p), r=r(p) are the coordinates of the point along 
c     the line:
c     rp=rma + (rst(i_cone)-rma)*p 
c     zp=zma + (zst(i_cone)-zma)*p
c
c     i_cone is the number of the EC cone
c           It is in antenna.i common block
c---------------------------------------------------------------       
c      implicit real*8 (a-h,o-z)
      implicit none
c-----input 
      real*8 p ! the parameter which determines the point
             ! on the straight line            

      include 'three.i'   ! rma,zma,psilim
      include 'param.i'     
      include 'cone.i'    !zst(ncone),rst(ncone) 
      include 'antenna.i' ! gives i_cone the number

c-----locals 
      real*8 rp,zp
c-----external
      real*8  psif  

      rp=rma + (rst(i_cone)-rma)*p 
      zp=zma + (zst(i_cone)-zma)*p
c      write(*,*)'psilim_psif_zp_rp p,i_cone', p,i_cone
c      write(*,*)'psilim_psif_zp_rp rp,zp', rp,zp
      psilim_psif_zp_rp= psilim-psif(zp,rp)
c      write(*,*)'psilim_psif_zp_rp end'
      return
      end


c======================================================================
c======================================================================


      subroutine antenna_equal_dist(first,i_cone,distance,
     &theta,z_ant,r_ant)
c-----------------------------------------------------------
      implicit none 
      include 'param.i' 
      include 'three.i'
      include 'cone.i'
      
c-----input      
      logical first !
      integer i_cone !the index of a given EC cone 
      real*8 theta, ! poloidal angle (radians)
     &distance      ! the distance from the plasma to antenna
                    ! along the sraight line 
                    ! rp=rma + (rst(icone)-rma)*p 
                    ! zp=zma + (zst(icone)-zma)*p      
c-----output
      real*8 z_ant,r_ant ! coordinates of the antenna vertex at
                         ! the given poloidal angle theta
                         ! for the given i_cone 
 
c-----local
      real*8 theta_ar(nteta),
     &z_ant_ar(nconea,nteta),r_ant_ar(nconea,nteta),
     &z_bound(nteta),r_bound(nteta),tabl(3),
     &temp_function(nteta),temp_deriv(nteta),pi,step,rho_bound,p,
     &workk(3*nteta+1),two_pi,theta_loc
      integer iop(2),itabl(3),i

      real*8 work(3*nteta+1),d2z_ant(nconea,nteta),d2r_ant(nconea,nteta)
c-----externals period_argument
      real*8 period_argument
      
      save d2z_ant,d2r_ant,theta_ar,r_ant_ar,z_ant_ar
    
      pi=4.d0*datan(1.d0) 

c      write(*,*)'antenna_equal_dist first,i_cone,distance,theta',
c     &first,i_cone,distance,theta

      step=2.d0*pi/(nteta-1.d0)
c      write(*,*)'antenna_equal_dist nteta,i_cone',nteta,i_cone
      if (first) then
c-------------------------------------------------------------------
c       creation of spline coefficients for the antenna surface
c-------------------------------------------------------------------
        do i=1,nteta-1
          theta_ar(i)=(i-1)*step
c          write(*,*)'antenna_equal_dist i,theta_ar(i)',i,theta_ar(i)
c---------creation mesh [z_bound,r_bound] along plasma boundary 
          !-> Get z(psilim,theta_ar),r(psilim,theta_ar)
          call zr_psith(psilim,theta_ar(i),z_bound(i),r_bound(i))
          rho_bound=dsqrt((z_bound(i)-zma)**2+(r_bound(i)-rma)**2)
          p=1.d0+distance/rho_bound
c---------creation mesh [z_ant_ar,r_ant_ar] along antenna surface
          z_ant_ar(i_cone,i)=z_bound(i)+(p-1.d0)*(z_bound(i)-zma)
          r_ant_ar(i_cone,i)=r_bound(i)+(p-1.d0)*(r_bound(i)-rma)
c          write(*,*)'antenna_equal_dist z_ant_ar(i_cone,i),
c     &    r_ant_ar(i_cone,i)',z_ant_ar(i_cone,i),r_ant_ar(i_cone,i)
        enddo !i

        theta_ar(nteta)= 2.d0*pi
        z_ant_ar(i_cone,nteta)= z_ant_ar(i_cone,1)
        r_ant_ar(i_cone,nteta)= r_ant_ar(i_cone,1)

c        do i=1,nteta
c        write(*,*)'i,theta_ar(i),z_ant_ar(i_cone,i),r_ant_ar(i_cone,i)',
c     &  i,theta_ar(i),z_ant_ar(i_cone,i),r_ant_ar(i_cone,i)
c        enddo
c----------------------------------------------------------------
c     
        iop(1)=3 ! periodic spline boundary conditions
        iop(2)=3

c-------creation spline coefficients [d2z_ant,d2r_ant] for the antenna surface
        call take_1d_array_from_2d_array(temp_function,z_ant_ar,
     &  nconea,ncone,nteta,nteta,i_cone)
        call coeff1(nteta,theta_ar,temp_function,temp_deriv,
     &              iop,1,workk)
        call put_1d_array_to_2d_array(temp_deriv,d2z_ant,
     &  nconea,ncone,nteta,nteta,i_cone)

        call take_1d_array_from_2d_array(temp_function,r_ant_ar,
     &  nconea,ncone,nteta,nteta,i_cone)
        call coeff1(nteta,theta_ar,r_ant_ar(i_cone,1),d2r_ant(i_cone,1),
     &              iop,1,workk)
        call put_1d_array_to_2d_array(temp_deriv,d2r_ant,
     &  nconea,ncone,nteta,nteta,i_cone)
     
        write(*,*)'antenna_equal_dist after spline coeff1' 
 
ctest : check splines
c         itabl(1)=1
c         itabl(2)=0
c         itabl(3)=0

c       call take_1d_array_from_2d_array(temp_function,r_ant_ar,
c     & nconea,ncone,nteta,nteta,i_cone)
c       call take_1d_array_from_2d_array(temp_deriv,d2r_ant,
c     & nconea,ncone,nteta,nteta,i_cone)
c       do i=1,nteta
c        call terp1(nteta,theta_ar,temp_function,temp_deriv,
c     &              theta_ar(i),1,tabl,itabl)
c         write(*,*)'i,theta_ar(i),tabl(1), r_ant_ar(i_cone,i)',
c     &   i,theta_ar(i),tabl(1), r_ant_ar(i_cone,i)
c       enddo

c       call take_1d_array_from_2d_array(temp_function,z_ant_ar,
c     & nconea,ncone,nteta,nteta,i_cone)
c       call take_1d_array_from_2d_array(temp_deriv,d2z_ant,
c     & nconea,ncone,nteta,nteta,i_cone)
c       do i=1,nteta
c        call terp1(nteta,theta_ar,temp_function,temp_deriv,
c     &              theta_ar(i),1,tabl,itabl)
c         write(*,*)'i,theta_ar(i),tabl(1), z_ant_ar(i_cone,i)',
c     &   i,theta_ar(i),tabl(1), z_ant_ar(i_cone,i)
c       enddo

cendtest      
      
      endif !first
     
      itabl(1)=1
      itabl(2)=0
      itabl(3)=0

      two_pi=8.d0*atan(1.d0)         
      theta_loc=period_argument(theta,two_pi)
     
      call take_1d_array_from_2d_array(temp_function,z_ant_ar,
     &nconea,ncone,nteta,nteta,i_cone)
      call take_1d_array_from_2d_array(temp_deriv,d2z_ant,
     &nconea,ncone,nteta,nteta,i_cone)
      call terp1(nteta,theta_ar,temp_function,temp_deriv,
     &           theta_loc,1,tabl,itabl)
      z_ant=tabl(1)

      call take_1d_array_from_2d_array(temp_function,r_ant_ar,
     &nconea,ncone,nteta,nteta,i_cone)
      call take_1d_array_from_2d_array(temp_deriv,d2r_ant,
     &nconea,ncone,nteta,nteta,i_cone)
      call terp1(nteta,theta_ar,temp_function,temp_deriv,
     &           theta_loc,1,tabl,itabl)
      r_ant=tabl(1)

c      write(*,*)'antenna_equal_dist theta,r_ant,z_ant',
c     &theta,r_ant,z_ant

      return
      end


c======================================================================
c======================================================================


      real*8 function x_bin_min(f,x1,x2,xacc)
c---------------------------------------------------------------
c     Using the binary metod founds the point x_min=x_bin_min
c     which gives the minimal value for f(x)
c     at (x1,x2). 
c     xacc is the accuracy of x_min determination
c------------------------------------------------------------- 
      implicit none
c-----input
      real*8 x1,x2,xacc
c-----externals
      real*8 f
      external f   
c-----locals
      real*8 delta,f1,f2,u1,u2
      integer i,i_max
     
      i=0
      i_max=30

      write(*,*)'in x_bin_min'

      delta=1.d-1*xacc

 10   continue

      u1=0.5d0*(x1+x2-delta)
      u2=0.5d0*(x1+x2+delta)
      
      f1=f(u1)
c      write(*,*)'x_bin_min after f1,u1,f1',u1,f1

      f2=f(u2)

      if (f1.le.f2) then
        x1=x1
        x2=u2
      else
        x1=u2
        x2=x2
      endif

      i=i+1
      if(i.gt.i_max) then
        write(*,*)'***************************************'
        write(*,*)'in x_bin_min' 
        write(*,*)'i>i_max i_max=',i_max
        write(*,*)'dabs(x2-x1)',dabs(x2-x1)
        write(*,*)'***************************************'
        goto 20
      endif 

      x_bin_min=0.5d0*(x1+x2)
      if (dabs(f(x_bin_min)).lt.xacc) goto 20

      if (dabs(x2-x1).gt.xacc) goto 10
          
 20   continue
      x_bin_min=0.5d0*(x1+x2)
      write(*,*)'x_bin_min=',x_bin_min
cSAP090814 may be it will help
c      write(*,*)'f(x_bin_min)=',f(x_bin_min)
      
      return
      end


c======================================================================
c======================================================================

    
      real*8 function f_antenna_min(p)
c----------------------------------------------------------------
c     Calulate the function:
c     f_antenna_min(theta)=(r_ray(p)-r_ant(theta(p)))**2+
c                         +(z_ray(p)-z_ant(theta(p)))**2)
c
c     Here
c     r_ant(theta),z_ant(theta) are the coordinates of the antenna
c     surface
c
c     r_ray(p),z_ray(p) are the coordinates at the vacume ray
c     - straight line starting at the plasma edge z_0,r_0, with the
c     direction n_0_r,n_0_z
c     
c     The ray equations:
c     r_ray(p)=rma+rho_ray(p)cos(theta_ray(p))
c     z_ray(p)=zma+rho_ray(p)sin(theta_ray(p))
c     rho_ray(p)sqrt[(r_ray(p)-rma)**2+(z_ray(p)-zma)**2]
c
c     Special cases.
c     1) n_0_r=0
c        In this case the ray equation is: r_ray=r_0
c     2) n_0_z=0
c        In this case the ray equation is: z_ray=z_0

c      implicit real*8 (a-h,o-z)
      implicit none

c-----input
      real*8 p !parameter along the ray

      real*8  r_0,phi_0,z_0,
     &n_vac_x,n_vac_y,n_vac_z
      common /ox/r_0,phi_0,z_0,
     &n_vac_x,n_vac_y,n_vac_z
      include 'three.i' !rma,zma

c-----locals
      real*8 pi, sin_phi,cos_phi,n_vac_r_s,r_ray,z_ray,rho_ray,
     &cos_theta,sin_theta, theta_ray,z_ant,r_ant,
     &f_antenna_min_l !for test

c      write(*,*)'in f_antenna_min(p)'

      pi=4.d0*datan(1.d0)
c------------------------------------------------------
      sin_phi=dsin(phi_0)
      cos_phi=dcos(phi_0)

      n_vac_r_s=n_vac_x**2+n_vac_y**2

      r_ray=dsqrt(r_0**2+p**2*n_vac_r_s+2.d0*p*r_0*
     &            (n_vac_x*cos_phi+n_vac_y*sin_phi))

      z_ray=z_0 + p * n_vac_z
      rho_ray=dsqrt((r_ray-rma)**2+(z_ray-zma)**2)

      cos_theta=(r_ray-rma)/rho_ray
      sin_theta=(z_ray-zma)/rho_ray
      if (cos_theta .gt.  1.d0)  cos_theta= 1.d0
      if (cos_theta .lt. -1.d0)  cos_theta=-1.d0

      if (sin_theta.ge.0.d0) then
         theta_ray=+dacos(cos_theta)
      else  
         theta_ray=-dacos(cos_theta)
      endif
      
      call antenna_surface(theta_ray,z_ant,r_ant)
      f_antenna_min=(r_ray-r_ant)**2+(z_ray-z_ant)**2
      return
      end

c======================================================================
c======================================================================


      subroutine take_1d_array_from_2d_array(ar_1d,ar_2d,
     &n1_max,n1,n2_max,n2,j)
c-----take 1D array ar_1d(i=1,n2)= ar_2d(j,i=1,n2), j=const
c     from 2D array ar_2d(n1,n2)
c-----input
      integer n1,n2,j,n1_max,n2_max     
      real*8 ar_2d(n1_max,n2_max)

c-----output
      real*8 ar_1d(n2_max)        

c-----locals`
      integer i 
     
      do i=1,n2
         ar_1d(i)=ar_2d(j,i)
      enddo
     
      return
      end


c======================================================================
c======================================================================


      subroutine put_1d_array_to_2d_array(ar_1d,ar_2d,
     &n1_max,n1,n2_max,n2,j)
c-----put 1D array ar_1d(i=1,n2)=> ar_2d(j,i=1,n2), j=const
c     to 2D array ar_2d(n1,n2)
c-----input
      integer n1,n2,j,n1_max,n2_max     
      real*8 ar_1d(n2_max)

c-----output
      real*8 ar_2d(n1_max,n2_max)        

c-----locals
      integer i 
     
      do i=1,n2
         ar_2d(j,i)=ar_1d(i)
      enddo
     
      return
      end


c======================================================================
c======================================================================


      real*8 function period_argument(x,period)
c-----calculate the argment in the period area
c     i= i=x/period
c     if (x.gt.0.d0) then      
c       period_argument=x-i*period
c     else
c       period_argument=x+(-i+1)*period
c     endif

      implicit none
c-----input
      real*8 x
      real*8 period ! >0

c-----local
      integer i     
 
      if (period.le.0.d0) then
         write(*,*)'in period_argument period.le.0'
         write(*,*)'it should period>0'
         stop 'in period_argument p'
      endif
  
      i=x/period
      if (x.gt.0.d0) then      
        period_argument=x-i*period
      else
        period_argument=x+(-i+1)*period
      endif

      return
      end


c======================================================================
c======================================================================


      subroutine find_maximal_OX_transmission(eps_xe,x0,
     & r_in,z_in,phi_in,nr_in,nz_in,m_in,
     & i_ox_conversion)
c-----if((xe.gt.(x0-eps_xe)).and.(xe.le.x0)) then it
c     calculates OX transmission coefficient transm_ox in u() point
c     and compares this coefficient with the transmission coeifficient
c     transm_ox_old at the previous point u_old()
c     Finds the point where the transmission coefficient has the maximal
c     value and put i_ox_conversion=1
      implicit none
c-----input
      real*8 eps_xe,x0,r_in,z_in,phi_in,nr_in,nz_in,m_in
     
      include 'param.i'  
      include 'output.i'
      include 'oxb.i'
c-----output 
c     r_in,z_in,phi_in,nr_in,nz_in,m_in the ray coordinates at
c             the point with the maxima transmission coefficient
      integer i_ox_conversion !=0 the ray is far from 
                              !   the OX conversion point
                              !=1 the ray is near the OX
                              !   conversion point
c-----externals
      real*8 x,b

c-----local
      real*8 xe,transm_ox,u(6),u_old(6),bmod
      integer i

      save  u_old

      write(*,*)'  find_maximal_OX_transmission'
      write(*,*)'eps_xe,r_in,z_in',eps_xe,r_in,z_in


      u(1)=z_in
      u(2)=r_in
      u(3)=phi_in
      u(4)=nz_in
      u(5)=nr_in 
      u(6)=m_in

      bmod=b(z_in,r_in,phi_in)
      xe=x(u(1),u(2),u(3),1)
      i_ox_conversion=0

      write(*,*)'xe,x0-eps_xe', xe,x0-eps_xe
      write(*,*)'was_in_ox_vicinity',was_in_ox_vicinity
      if((xe.gt.(x0-eps_xe)).and.(xe.le.x0)) then
         write(*,*)'before transmit_coef_ox'
         was_in_ox_vicinity=.true.

         call transmit_coef_ox(u(1),u(2),u(3),u(4),u(5),u(6),
     &                           transm_ox)

         write(*,*)'after transmit_coef_ox transm_ox,transm_ox_old',
     &    transm_ox,transm_ox_old

         if (transm_ox.lt.transm_ox_old) then
            i_ox_conversion=1
            z_in=u_old(1)
            r_in=u_old(2)
            phi_in=u_old(3)
            nz_in=u_old(4)
            nr_in=u_old(5) 
            m_in=u_old(6)
            goto 10
         endif

         do i=1,6
            u_old(i)=u(i)
         enddo

         transm_ox_old=transm_ox 
      else
         if((xe.gt.x0).and.was_in_ox_vicinity) then
           write(*,*)'xe.gt.x0).and.was_in_ox_vicinity'
           transm_ox=transm_ox_old 
           i_ox_conversion=1
            z_in=u_old(1)
            r_in=u_old(2)
            phi_in=u_old(3)
            nz_in=u_old(4)
            nr_in=u_old(5) 
            m_in=u_old(6)
            write(*,*)'transm_ox',transm_ox
            goto 10
          endif
      endif

 10   continue
      write(*,*)'after 10 i_ox_conversion',i_ox_conversion
      return
      end     


c======================================================================
c======================================================================


      subroutine set_oxb
c-----initialise oxb.i at each new ray
      implicit none
      include 'param.i'
      include 'oxb.i'
      
      was_in_ox_vicinity=.false.
      transm_ox_old=0.d0
      i_call_prep3d_in_output_at_i_ox_conversion_eq_1=0
      return
      end
