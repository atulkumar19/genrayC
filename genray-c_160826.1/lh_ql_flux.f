






     
      real*8 function averaged_density(nmax,imax,ua,th0a,fm,lam1a)
c-------------------------------------------------------------
c     calculates flux averaged density from maxwellian distriburtion fm
c     dens0=2pi*integral{0,pi}d_theta_0*sin(theta_0)
c               integral{0,infinity}u_0**2d_u_0
c               lambda(rho,u_0,theta_0)*f_maxwell(T,u_0)
c-------------------------------------------------------------
      implicit none
c-----input
      integer
     & nmax,                       !number of points in momentum mesh
     & imax                        !number of points in pitsch angle mesh
      real*8
     &ua(0:nmax - 1),              !momentum mesh 
     &th0a(0:imax - 1),            !pitch angle mesh         
     &fm(0:nmax-1),                ! relativistic maxwellian distribution
     &lam1a(0:imax)                !cos(theta_0)*integral[dl/cos(theta)]
                                   !cos(theta)=dsqrt(1-b(l)sin(theta_0)**2)
c-----locals
      integer i,n
      real*8 pi,du,dth0,sum_theta,sum_u

      pi=4.d0*datan(1.d0)
      
      du=ua(2)-ua(1)
      dth0=th0a(2)-th0a(1)
 
      write(*,*)'in fuction averaged_density nmax,imax',nmax,imax
c      write(*,*)'ua',ua
c      write(*,*)'th0a',th0a
      write(*,*)'du,dth0,pi,imax',du,dth0,pi,imax  
      sum_theta=0.d0
      do i= 1, imax-1
         write(*,*)'i,lam1a(i),th0a(i),dsin(th0a(i))',
     &               i,lam1a(i),th0a(i),dsin(th0a(i))
         sum_theta=sum_theta+lam1a(i)*dsin(th0a(i))
      enddo
      sum_theta=sum_theta*2.d0 
      write(*,*)'sum_theta',sum_theta
      write(*,*)'nmax',nmax
      write(*,*)'ua',ua
      write(*,*)'fm',fm 

      sum_u = dot_product(ua(:nmax-1)**2,fm(:nmax-1))

      write(*,*)'sum_u',sum_u     

      averaged_density=2.d0*pi*sum_theta*sum_u*du*dth0
      write(*,*)'averaged_density', averaged_density

      return
      end

      function current1 (du, ua, c, nmax, dth0, th0a, imax, imax1,
     1   fm, chi, cur)
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  09:55:56   6/14/01
C...Switches: -p4 -yb
      use kind_spec

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(kind=dp) current
      integer nmax, imax, imax1
      real(kind=dp) du, c, dth0
      real(kind=dp), dimension(0:nmax - 1) :: ua
      real(kind=dp), dimension(0:imax - 1) :: th0a
      real(kind=dp), dimension(0:nmax - 1) :: fm
      real(kind=dp), dimension(0:imax1 - 1,0:nmax - 1) :: chi
      real(kind=dp), dimension(0:imax - 1) :: cur
      real(kind=dp)  current1 !SAP080101
      real(kind=dp), dimension(0:imax - 1) :: cur1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: n, i
      real(kind=dp) :: pi
ctest
      integer :: j
      real(kind=dp) :: chi_norm
C-----------------------------------------------

      pi = 4*datan(1.0d0)
      cur(:imax-1) = 0.d0
      cur1(:imax-1) = 0.d0

c      write(*,*)'current1 imax,imax1,nmax',imax,imax1,nmax

c      write(*,*)'in current1 chi(0,130)',chi(0,130) 

      do n = 0, nmax - 1
         cur(:imax-1) = cur(:imax-1) + chi(:imax-1,n)*fm(n)*ua(n)**3/
     1      dsqrt(1 + (ua(n)/c)**2)           
      end do

      do i=0,imax-1
        do n=0,nmax-1
          cur1(i) = cur1(i) + chi(i,n)*fm(n)*ua(n)**3/
     1      dsqrt(1 + (ua(n)/c)**2)
c          write(*,*)'i,n,ua(n),fm(n),chi(i,n)',i,n,ua(n),fm(n),chi(i,n)
c          write(*,*)'cur1(i)',cur1(i)  
        enddo
      enddo

c      do i=0,imax
c        write(*,*)'i,cur(i),cur1(i)',i,cur(i),cur1(i)
c      enddo   

      current1=0.d0  
      do i=0,imax-1
         current1=current1+cur1(i)*dsin(th0a(i))*dcos(th0a(i))
      enddo

      current=dot_product(dsin(th0a(:imax-1))*dcos(th0a(:imax-1)),cur(:
     1   imax-1)) 
c      write(*,*)'current,current1',current,current1
       
      current = 4*pi*du*dth0*current

      current1 = 4*pi*du*dth0*current1
      write(*,*)'current,current1',current,current1

      return
      end function current1

     
   





      subroutine absorbed_power_using_ql_flux(cnpar,cnper,z_m,r_m,phi,
     &fluxn,power,del_s_poloidal,
     &absorbed_power)

      implicit none
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'cefield.i'    

c-----input
      real*8
     &cnpar,cnper,       !parallel and perpendicular refractive index
     &z_m,r_m,           !space ray coordinates [meter]
     &phi,               !toroidal angle or the ray point [radian]
     &power,             !power in ray channel [erg/sec]
     &del_s_poloidal     !poloidal length of the ray element []

      real*8 fluxn       !power flux at unit electric field |E|=1 
                         !flux=B~(i)*B(i)+E~(i)*d(omega*eps(i,j))/domega*E(j)
                         !fluxn=0.5*dreal(flux*dconjg(flux))*vgrpol 
                         !vgrpol is a poloidal group velocity
                         !normalized to clight   
c-----output          
      real*8
     &absorbed_power     !erg/sec the power abdsorbed at the ray element
   
c-----external
      real*8
     &psi_rho,
     &bmin_psi

c-----locals
      real*8
     &rho_loc,cos_theta_pol,sin_theta_pol,theta_pol,
     &psi_loc,unorm,
     &cd_efficiency,pow_dens_0_tilda,cd_dens,
     &b_pol,bmin,
     &clight,s_poloidal_tilda,
     &absorbed_power_dev_s_pol               !erg/[sec*cm] the power abdsorbed at the ray element 
                                             !of poloidal length delta_s_pol [cm]
                                             
c     &vgr_poloidal
c     &length_b

      integer n_theta_pol,
     &n_radial !can be arbitrary < n_radial_a

      pi=4.d0*atan(1.d0)
c      write(*,*)'pi',pi
      rho_loc=dsqrt((r_m-rma)**2+(z_m-zma)**2) 

      write(*,*)'in absorbed_power_using_ql_flux rho_loc',rho_loc

      cos_theta_pol=(r_m-rma)/rho_loc
      sin_theta_pol=(z_m-zma)/rho_loc

      if (cos_theta_pol .gt. 1.d0)  cos_theta_pol= 1.d0
      if (cos_theta_pol .lt. -1.d0)  cos_theta_pol=-1.d0
      
      if (sin_theta_pol.ge.0.d0) then
         theta_pol=+dacos(cos_theta_pol)
      else  
c---------it should be 0 =< theta_pol< 2pi
         theta_pol=2.d0*pi-dacos(cos_theta_pol)
      endif 

c      write(*,*)'theta_pol',theta_pol

      if (theta_pol.lt.0.d0) then
         theta_pol=theta_pol+2.d0*pi
      endif

      if (theta_pol.gt.2.d0*pi) then
         n_theta_pol=theta_pol/(2.d0*pi)
         theta_pol=theta_pol-2.d0*pi*n_theta_pol
      endif

      n_radial=1
      psi_loc=psi_rho(rho_loc)  ! YuP: error? was rho
      !pause

      call QL_power_1(
     &psi_loc,
     &theta_pol,cnpar,cnper,cex,cey,cez,
     &unorm,
     &pow_dens_0_tilda)

      bmin= bmin_psi(psi_loc)
c-----lengt along B field: length_b [meters]
c      call calc_length_b(psi_loc,length_b)
c-----total group velocity
      clight=2.99792458d10              !light speed [cm/sec]

c      vgr_poloidal=dsqrt(vgr(1)**2+vgr(2)**2)! poloidal group velocity /clight
c      vgr_poloidal=vgr_poloidal*clight !sm/sec]

      s_poloidal_tilda= fluxn*clight/(8.d0*pi) 

      absorbed_power_dev_s_pol=power*(bmod/bmin)*
     &         (pow_dens_0_tilda/s_poloidal_tilda)
           
      absorbed_power=absorbed_power_dev_s_pol*del_s_poloidal ![erg/sec]
      return
      end





      subroutine QL_power_1(
     &psi_in,
     &theta_pol,n_parallel,n_perp,e_x,e_y,e_z,
     &unorm, 
     &pow_dens_0_tilda)
c------------------------------------------------------
c     calculates bounce averaged power:  pow_dens 
c     at unit wave elecric field vector |E|=1
c------------------------------------------------------
      implicit none
      include 'param.i'
      include 'one_nml.i'
      !!! include 'adj.i'
      include 'write.i'
      
c-----input
 
      real*8  
     &theta_pol,           !poloidal angle [radian]    
     &n_parallel,          !parallel refractive index
     &n_perp,              !perpendicular refractive index 
     &psi_in               !poloidal flux
      complex*16 
     &e_x,
     &e_y,                 !electric field polarization in Stix frame
     &e_z                  !at |E|=1

c-----output
      real*8
     &unorm,               !sqrt(T/m) cm/sec
     &cd_efficiency,
     &pow_dens_0_tilda,     ![1/sec], see manual: absorbed power calculation
c                          !Power density normalized to m*u_norm**2*n/tau_n 
     &cd_dens              !current drive density normalized to q*u_norm*n
                           !Here n is the density  
c-----external
      real*8 temperho,densrho,
     &bmax_psi,bmin_psi,psi_rho,rhopsi,y,b

c-----local
      integer n_radial0,n_harm_adj,
c     &i_resonance_curve_integration_method,
     &kmax,
     &i_calculate_CD     ! =0, calculate power only (no CD) 
                         ! using input arguments computed at
                         ! given space point
                         ! In this case n_radial0 can be arbitary
      real*8 u0,     
     &rho_small,           !small radius
     &temp_kev,            !electron temperatura [KeV]
     &dense,               !electron dencity 10**13/cm**3     
     &b_ratio,             !B/B_min
     &pi,bmin,bmax,deltb,psi,th0max,sin_trap,bmod,
     &d_chi_du0,d_chi_dtheta0,u0max,chi,
     &clight,             !light speed [cm/sec]
     &arg1,cln,
     &z,r,phi,            !bellow phi is set arbitratry phi=0
     &y_loc,              !electron y at point (z,r,phi)     
     &power_nharm_adj,CD_nharm_adj, !under integral functions
     

     &p_par_rl(2),         !for check max and min p_par
     &k_1_kev,             !egrs in 1 KeV      (erg)
     &charge_electron,     !electron charge (statcoulomb)
     &mass_e,              !electron mass [g]
     &theta_temperature,   !m_e*c^2/T_e
     &pow_dens,            
                           !using QL flux
     &omega                !2*pi*frqncy*1.d9 [1/sec]

      k_1_kev=1.6022d-9          !egrs in 1 KeV      (erg)
      charge_electron=4.8032d-10 !electron charge (statcoulomb)
      pi=4.d0*datan(1.d0)
      clight=2.99792458d10              !light speed [cm/sec]    
      mass_e=9.1094d-28         !electron rest mass (g)

      rho_small=rhopsi(psi_in)

      psi=psi_in

      temp_kev=temperho(rho_small,1)
      dense=densrho(rho_small,1) 
ccc      unorm=dsqrt(temp_kev/t)*dsqrt(1.6022d0/0.91095d0)*1.d9 ![sm/sec] sqrt(T/m)
      bmax= bmax_psi(psi)
      bmin= bmin_psi(psi)
	     
      deltb = (bmax-bmin)/bmin
      th0max = atan2(1.0d0,sqrt(deltb))
      sin_trap=dsqrt(bmin/bmax)


      pow_dens=0.d0
      cd_dens=0.d0

      phi=0.d0
      call zr_psith(psi_in,theta_pol,z,r) !get z(psi,theta),r(psi,theta)

      bmod=b(z,r,phi)  
      b_ratio=dmax1(bmod/bmin,1.d0)

      y_loc=y(z,r,phi,1)

      pow_dens=2.d0*pi*pow_dens
      cd_dens=-2.d0*pi*cd_dens

c      write(*,*)'pow_dens',pow_dens
c      write(*,*)'cd_dens',cd_dens
     
      if(pow_dens.ne.0.d0) then
        cd_efficiency=cd_dens/pow_dens
      else
        cd_efficiency=0.d0
      endif

      write(*,*)'in CD_adj_efficiency_1 cd_efficiency=',
     &cd_efficiency
         
c-----cln -Coulomb Ln
      arg1 = 1.d3/temp_kev*dsqrt(10.d0*dense)	!temp KeV, dense 10**13/cm**3
      cln=24.d0-dlog(arg1)

c-----adj efficiency is normalized to eta_0_adj=e*tau_n/(m*u_t) 
c     tau_n=T**(3/2)*sqrt(m)/(4*pi*e**4*density*cln)
c     u_t=cqrt(T/m)
c     eta_0_adj=T/[4*pi*e**3*n*cln]|cgs=T_kev/(n_13*cln)*
c           [k_1_kev/(4*pi*charge_electron**3*10**13)]           !cgs
c      coef=1.6022d+8/(4*pi*4.8032**3)          !cgs  (statampere/cm**2)/(erg/(sec*cm**3))
c      coef=1.6022d-1/(4*pi*4.8032**3*3)=3.8352d-5!     (A/cm**2)/(erg/(sec*cm**3))
c      eta_0_adj=(temp_kev/(cln*dense))*coef

      cd_efficiency=cd_efficiency*
     &3.84d0*temp_kev/(cln*dense)	       !  (A/m**2)/(joule/(sec*m**3))
c      write(*,*)'temp_kev,dense,cln,cd_efficiency',
c     &temp_kev,dense,cln,cd_efficiency
      cd_efficiency=cd_efficiency*1.d-5  !  (A/cm**2)/(erg/(sec*cm**3))
      write(*,*)'temp_kev,dense,cln,cd_efficiency',
     &temp_kev,dense,cln,cd_efficiency

c---------------------------------------------------------------------
cSAP081201
c     
c     The used relativistic Maxwellian distribution was normalized
c     to  unit density electron density:
c     4*pi*integral{f_m*(p/mc)^2*d(p/mc)}=1
c---------------------------------------------------------------------
  
      omega=2.d0*pi*frqncy*1.d9 ![1/sec]
c      write(*,*)'frqncy,omega',frqncy,omega
c---- power_dens_0_tilda has dimension [1/sec]
c      write(*,*)'pow_dens',pow_dens
c      write(*,*)'dense*1.d13',dense*1.d13
      pow_dens_0_tilda=pow_dens*
     &dense*1.d13*      !the power was calculated for unit density
                        !This term transforms power to local density
     &(unorm/clight)*   !it was the integration by (p_perp/mc)*d(p_perp/mc)
     &charge_electron**2/mass_e/omega !in subroutine intgr_relt_power_cd_adj
                        !This transforms this integration to
                        !(p_perp/(m*unorm)*)d(p_perp/(m*unorm))
                        !and uses all other normalizations
                        !in under the integral term
c      write(*,*)'charge_electron**2/mass_e/omega [cm^3/sec] ',
c     &           charge_electron**2/mass_e/omega
c       write(*,*)'pow_dens_0_tilda',pow_dens_0_tilda,'[1/sec]'
      return
      end

     

      subroutine calc_length_b(psi,length_b)
c-----calculate length at flux surface along B field line
c     when the poloidal angle is changed at 2*PI
      implicit none 
      include 'param.i'
      include 'gr.i'
      include 'one.i'
      include 'three.i'
      include 'rho.i'
c      real*8 length_b(npsi) !length along B
c-----input
      real*8 psi !polidal flyx
c-----output
      real*8  length_b ! length of B field line
c-----externals
      real*8 b
c-----locals
      integer j,i
      real*8 theta,z,r,phi,b_pol,dl_pol,dl_b,psi_loc
      logical first

c----------------------------------------------------------
c     to use spline functions from zcunix.f:
c     coeff1 and terp1
      integer i1p(2)
      real*8, dimension(1:3*npsi+1) :: work_l
      integer itabl(3)
      real*8  tabl(3)

      data first /.true./
      save first

          
      if (first) then
     
         length_b_ar(1)=0.d0  !length along B field
    
         phi=0.d0   
         do j=2,npsi           
            psi_loc=arpsi(j)            
            length_b_ar(j)=0.d0  !initialization of the length 
                           !along B field
	    do i=1,nteta
	       theta=arteta(i)
               call zr_psith(psi_loc,theta,z,r)
               bmod=b(z,r,phi)
               b_pol=dsqrt(bz**2+br**2)

c               write(*,*)'j,i,bmod,b_pol,bmod/b_pol',
c     &                    j,i,bmod,b_pol,bmod/b_pol

               dl_pol=dsqrt((rpsi(j,i+1)-rpsi(j,i))**2+
     &             (zpsi(j,i+1)-zpsi(j,i))**2)
               dl_b=dl_pol*bmod/b_pol        
               length_b_ar(j)=length_b_ar(j)+dl_b  ![meter]
            enddo !i 
           
         enddo !j
 
c------------------------------------------------------------------
c        spline coefficient calculation for length_b(psi) using spline
c        from zcunix.f
c---------------------------------------------------------------
         i1p(1)=4 !first derivatives at the left boundary will be fitted

         i1p(2)=4 !first derivatives at the right boundary will be fitted

         call coeff1(npsi,arpsi,length_b_ar,d2_length_b_psi,i1p,1,
     &               work_l)
 
         first = .false.
      else
         if (psi.lt.psimag) then
            length_b=0.d0
            goto 10
         endif
c---------------------------------------------------
c        b field line calculation using spline from zcunix.f
c---------------------------------------------------
         itabl(1)=1
         itabl(2)=0
         itabl(3)=0
          
         call terp1(npsi,arpsi,length_b_ar,d2_length_b_psi,
     &             psi,1,tabl,itabl)
         
         length_b=tabl(1)
      endif

 10   continue
      return
      end
     



