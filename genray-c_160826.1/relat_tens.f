
      subroutine anth_rlt(X,Y,T_kev,nll_in,np_in,
     +n_relt_harm1,n_relt_harm2,n_relt_intgr,
     +i_resonance_curve_integration_method,epsi,
     +i_fkin,r,z,phi,
     +aK)
c     calculates anti-hermition relativistic dielectric tensor aK
c     for electron plasma
c
c     INPUTS:
c
c      X = (fpe/f)**2
c      Y = fce/f
c      T_kev  - electron temperature
c      nll_in - parallel index of refraction N.
c      np_in  - perpendicular refractive index N
c     n_relt_harm1 min number of EC harmonics used in anti-hermitian
c               dielectric tensor calculations
c     n_relt_harm2 max number of EC harmonics used in anti-hermitian
c               dielectric tensor calculations
c               n_relt_harm1 =< n= <n_relt_harm2
c      n_relt_intgr - the number of points for the integration over p_perp
c      i_fkin =0 the usage of the analytical relativistic distributin
c             =1 the usage of the numerical 3D distribution from diskf or netcdfnm.nc 
c                written be CQL3D code or created analytically at mesh points
c      r      - the major radius (it is used for i_fkin=1)
c      z      - the vertical coordinate  (it is used for i_fkin=1)
c      phi    - toroidal angle
c-------------------------------------------------------------------
!       i_resonance_curve_integration_method=1 !rectangle integration
!                                              !over angle,
!                                              !for ellipse case only
!       i_resonance_curve_integration_method=2 !rectangle formula,
!                                              !p_perp integration
!                                              !elliptic or hyperbolic case
!       i_resonance_curve_integration_method=3 !trapezoidal formula,
!                                              !p_perp integration
!                                              !elliptic or hyperbolic case
!       i_resonance_curve_integration_method=4 !adaptive Simpson integration 
!                                              !p_perp integration
!                                              !elliptic or hyperbolic case
!
!       i_resonance_curve_integration_method is used in subroutine intgr_rl
!       to choose the numerical integration method for 
!       anti-hermitian relativistic dielectric tensor calculation.
!       This applies for iabsorp=6,7 and for emission calculations.
c--------------------------------------------------------------------
c     OUTPUT:
c      aK(3,3):  the nine components of the anti-hermition
c                part dielectric relativistic tensor
c                evaluated at (X,Y,Te,nll,np,n)

      implicit none
c     input 
      double precision X,Y,T_kev,nll_in,np_in,epsi
      integer n_relt_harm1,n_relt_harm2,n_relt_intgr,i_fkin,
     &i_resonance_curve_integration_method

      double precision r,z,phi

c     output
      double complex aK(3,3)
c     local
      double precision c,mass_e,k_1_kev,nll,np,nlls,p_perp0,theta,pi,
     .dens,xpidens,rho,psi
      double complex integral(3,3)
      integer jn,ires

c     external zeroK, npnllmin, ec_cond, intgr_rl,fdens_fdist
      double precision fdens_fdist,densrho,rhopsi,psif
 
      c =2.99792458d10          !speed of light     (cm/sec)
      mass_e=9.1094d-28         !electron rest mass (g)
      k_1_kev=1.6022d-9         !egrs in 1 KeV      (erg)

      theta=mass_e*c**2/(k_1_kev*T_kev)

c      write(*,*)'relt_tens.f sub anth_rlt: r,z,T_kev,theta',
c     &r,z,T_kev,theta
cyup      write(*,*)'relt_tens.f sub anth_rlt:'
cyup      write(*,*)'i_resonance_curve_integration_method',
cyup     & i_resonance_curve_integration_method

      pi=4.d0*datan(1.d0)

      call npnllmin(nll_in,np_in,nll,np)
      nlls = nll**2
 
c-----initialization aK=0
      call zeroK(aK)

cyup      write(*,*)'n_relt_harm1,n_relt_harm2',n_relt_harm1,n_relt_harm2
c-----the loop over the cyclotron harmonics 
      do jn=n_relt_harm1,n_relt_harm2               
c-------control the EC resonance condition
cSm060315 -jn         
        call ec_cond(-jn,Y,nll,ires,p_perp0)

cyup        write(*,*)'anth_rl jn,ires,p_perp0',jn,ires,p_perp0

        if(ires.eq.0) then
c---------no resonace
          goto 10
        endif
cSm060315 -jn      
c        write(*,*)'relt_tens.f sub anth_rlt jn=',jn
cyup        write(*,*)'bef intgr_rl i_resonance_curve_integration_method',
cyup     &  i_resonance_curve_integration_method
   
        call intgr_rl(-jn,nll,np,Y,theta,ires,p_perp0,n_relt_intgr,
     &  i_resonance_curve_integration_method,epsi,
     +  i_fkin,r,z,phi,
     +  integral)

c        write(*,*)'after intgr_rl i_resonance_curve_integration_method',
c     &  i_resonance_curve_integration_method
   
c        write(*,*)'anth_rl,integral',integral
       
        aK(1,1)=aK(1,1)+integral(1,1)
        aK(1,2)=aK(1,2)+integral(1,2)
        aK(1,3)=aK(1,3)+integral(1,3)
        aK(2,2)=aK(2,2)+integral(2,2)
        aK(2,3)=aK(2,3)+integral(2,3)
        aK(3,3)=aK(3,3)+integral(3,3)

10      continue      
      enddo   !jn

      if (i_fkin.eq.0) then       
        dens=1.d0
      else      
        psi=psif(z,r)
        rho=rhopsi(psi)
        dens=densrho(rho,1)
      endif
      
c-----normalization for the unit electron density
      xpidens=X*pi/dens
      
      aK(1,1)=-xpidens*aK(1,1)
      aK(1,2)=-xpidens*aK(1,2)
      aK(1,3)=-xpidens*aK(1,3)
      aK(2,2)=-xpidens*aK(2,2)
      aK(2,3)=-xpidens*aK(2,3)
      aK(3,3)=-xpidens*aK(3,3)
      
      aK(2,1)=-aK(1,2)
      aK(3,1)= aK(1,3)
      aK(3,2)=-aK(2,3)

      return
      end

      double complex function dhot_rlt(K,aK,nll,np,nprim)
c     calculates complex dispersion function dhot_rl
c     using the non-relativistic hermition tesor K and 
c     relativistic anti-hermition dielectric tensor aK
c     for electron plasma for complex Im_N_perp=(np,nprim)
c     input
c       K(3,3)   complex hermition tensor
c       aK(3,3)  complex anti-hermition relativistic tensor
c       nll      Re(N_parallel)
c       np       Re(N_perpendicular)
c       nprim    Im(N_perpendicular)

      implicit none
c     input
      double complex K(3,3), aK(3,3)
      double precision nll,np,nprim
c     local
      double complex Kxx,Kxy,Kxz,Kyx,Kyy,Kyz,Kzx,Kzy,Kzz
      double complex cnp,nps
      double precision nlls
      double complex K_herm(3,3)
      integer i1,j1	
      double complex i

c      write(*,*)'dhot_rlt nll,np,nprim',nll,np,nprim		
c      write(*,*)'dhot_rlt K',K
c      write(*,*)'dhot_rlt aK',aK

      i=dcmplx(0.d0,1.d0)

c       compute k_herm the hermitian part of the dielectric tensor k_sum
        call herm(K,K_herm)

c	Kxx=K(1,1)+aK(1,1)*i
c	Kxy=K(1,2)+aK(1,2)*i
c	Kxz=K(1,3)+aK(1,3)*i
c	Kyx=K(2,1)+aK(2,1)*i
c	Kyy=K(2,2)+aK(2,2)*i
c	Kyz=K(2,3)+aK(2,3)*i
c	Kzx=K(3,1)+aK(3,1)*i
c	Kzy=K(3,2)+aK(3,2)*i
c	Kzz=K(3,3)+aK(3,3)*i
cSm060707
	Kxx=K_herm(1,1)+aK(1,1)*i
	Kxy=K_herm(1,2)+aK(1,2)*i
	Kxz=K_herm(1,3)+aK(1,3)*i
	Kyx=K_herm(2,1)+aK(2,1)*i
	Kyy=K_herm(2,2)+aK(2,2)*i
	Kyz=K_herm(2,3)+aK(2,3)*i
	Kzx=K_herm(3,1)+aK(3,1)*i
	Kzy=K_herm(3,2)+aK(3,2)*i
	Kzz=K_herm(3,3)+aK(3,3)*i

	nlls=nll*nll
        cnp=dcmplx(np,nprim)
	nps=cnp*cnp

      !write(*,*)'nlls,nps',nlls,nps 
      !write(*,*)'Kxx,Kxy,Kxz',Kxx,Kxy,Kxz 
      !write(*,*)'Kyx,Kyy,Kyz',Kyx,Kyy,Kyz 
      !write(*,*)'Kzx,Kzy,KXz',Kzx,Kzy,Kzz    
         	
	dhot_rlt =(Kxx-nlls) * (Kyy-nlls-nps) * (Kzz-nps)
     .+  Kxy * Kyz * (Kzx+cnp*nll) 
     .+ (Kxz+cnp*nll) * Kyx * Kzy
     .- (Kzx+cnp*nll) * (Kyy-nlls-nps) * (Kxz+cnp*nll)
     .-  Kzy * Kyz * (Kxx-nlls)
     .- (Kzz - nps) * Kyx * Kxy
      !write(*,*)'dhot_rlt',dhot_rlt
      return
	end

      subroutine ec_cond(n,Y,nll,ires,p_perp0)
c     controls EC resonance condition gamma=N_par*p_par+nY
c     and calcultes the max value p_perp0 of the perpendicular moment
c     for the ellipse case
c     input 
c      n    - the number of EC harmonic
c      Y    - omega_ce/omega for the electron rest mass
c      nll  - N_parallel to the magnetic field
c     output
c      ires =0 no EC resonance
c           =1 resonance ellipse
c           =2 resonance parabola
c           =3 resonance hyperbole
c      p_perp0 - is the max perpendicular to the magnetic field momentum 
c                divided by (mc) for the ellipse case;
c                m is the electron rest mass; c is the light speed. 
 
      implicit none
c     input
      integer n
      double precision Y,nll
c     output
      integer ires
      double precision p_perp0
c     local
      double precision nlls,nY,nYs

      nY=n*Y 
      nYs=nY*nY     
      nlls=nll*nll

      if(nlls.lt.1.d0)then
c       resonace ellipse
c-------------------------------------------------------
        if(((nlls+nYs).lt.1.d0).or.(nY.le.0.d0)) then
c         no resonace
          ires=0
        else    
          ires=1   
          p_perp0=dsqrt((nlls+nYs-1.d0)/(1.d0-nlls))
        endif
c-------------------------------------------------------
      else  
        if(nlls.gt.1.d0)then  
c         resonance hyperbole
c-------------------------------------------------------
          ires=3
c-------------------------------------------------------
        else 
c         nlls.eq.1.d0
c         resonace parabola 
c-------------------------------------------------------
          if(nY.le.0.d0) then
c           no resonace 
            ires=0
          else
            ires=2 
          endif
c--------------------------------------------------------
        endif 
      endif

      return
      end

      subroutine g_n_calc(p_perp,y,nll,np,theta,n,i_fkin,r,z,phi,g_n)
c     under integral complex marix function G_n(3,3)
c     input
c       p_perp  - momentum divided by (mc)
c       y       = omega_ce/omega for the electron rest mass    
c       nll     - N_parallel
c       np      - N_perpendicular 
c       theta   - mc**2/T
c       n       - n number of the given resonance harmonic
c       i_fkin =0 the usage of the analytical relativistic distributin
c              =1 the usage of the numerical 3D distribution from diskf 
c                 written be CQL3D code
c       r       - the major radius (it is used for i_fkin=1)
c       z       - the vertical coordinate  (it is used for i_fkin=1)
c       phi     - toroidal angle(it is used for i_fkin=1)
c     output
c       g_n(3,3) under integral matrix double complex function
c-------------------------------------------------------
      implicit none
c-----input        
      double precision p_perp,y,nll,np,theta,r,z,phi
      integer n,i_fkin
c-----external root_res, s_calcn, zeroK, u_n
      double precision u_n    

c-----output
      double complex g_n(3,3) 

c-----local
      double complex sn_k(3,3),g_nk(3,3)
      double precision gamma, p_par_rl(2),coeff,eps,ps,ps_max,pi,
     +dgdp_par,u_nk
      integer k,kmax,i,j
c     kmax - total number of the resonace condition root p_perp_k(p_perp)
      
    

      pi=4.d0*datan(1.d0)
c-----calculations ot the roots p_par_rl of the resonance condition
c-----gamma=N_par*p_par_rl+nY

      call root_res(p_perp,nll,n,Y,kmax,p_par_rl)

c      write(*,*)'in g_n_calc: p_perp, nll,n,Y,kmax,p_par_rl'
c     &,p_perp, nll,n,Y,kmax,p_par_rl
c-----initialize g_n
      call zeroK(g_n)

      if (kmax.eq.0) goto 20
  
c      eps=1.d-9 ! the min value of Maxwell exponent
c      p s_max=(1.d0-dlog(eps)/theta)**2-1.d0
c      write(*,*)'g_n_calc kmax',kmax      
      do k=1,kmax ! the sum over all P_parallel_k(perp) roots 
   
         ps=p_perp*p_perp+p_par_rl(k)*p_par_rl(k)
         gamma=dsqrt(1.d0+ps)

c         if(ps.gt.ps_max) then
c           write(*,*)'relat_tens g_n k,ps,ps_max',k,ps,ps_max
c           goto 10
c         endif
 
c        write(*,*)'g_n_calc before s_calcn: k',k  
         call s_calcn(n,p_perp,p_par_rl(k),np,y,sn_k)
c         write(*,*)'g_n_calc after s_calcn: k,sn_k',k,sn_k  

c--------resonance condition uses the delta function with argument
c        g_delta=gamma-nll*p_par-n*y
c        the derivative from this argument d(g_delta)/dp_par=dgdp_par

c         dgdp_par=(p_par_rl(k)-nll*y)/gamma
         dgdp_par=dabs((p_par_rl(k)-nll*gamma)/gamma)

c         write(*,*)'in g_n_calc before u_nk theta',theta

         u_nk=u_n(p_perp,p_par_rl(k),y,nll,theta,n,i_fkin)

cSm0960306
c        write(*,*)'g_n   dgdp_par',  dgdp_par

         coeff=2.d0*pi*u_nk/dgdp_par

c         write(*,*)'g_n k,p_perp,p_par_rl(k)',k,p_perp,p_par_rl(k)
c         write(*,*)'g_n dsqrt(ps),u_nk',dsqrt(ps),u_nk
c         write(*,*)'g_n u_nk,dgdp_par,coeff',u_nk,dgdp_par,coeff 
c         write(*,*)'coeff',coeff           

         g_nk(1,1)=coeff*sn_k(1,1)  
         g_nk(1,2)=coeff*sn_k(1,2)
         g_nk(1,3)=coeff*sn_k(1,3)
         g_nk(2,2)=coeff*sn_k(2,2)
         g_nk(2,3)=coeff*sn_k(2,3)
         g_nk(3,3)=coeff*sn_k(3,3)

         g_nk(2,1)=-g_nk(1,2)
         g_nk(3,1)= g_nk(1,3)
         g_nk(3,2)=-g_nk(2,3)        

         do i=1,3
            do j=1,3
               g_n(i,j)=g_n(i,j)+g_nk(i,j)
            enddo
         enddo

c         write(*,*)'in g_n: g_nk',g_nk

 10      continue
       
c         write(*,*)'in g_n_calc g_nk',g_nk
      enddo !kmax

 20   continue

      return
      end


      subroutine intgr_rl(n,nll,np,Y,theta,ires,p_perp0,
     +n_relt_intgr,i_resonance_curve_integration_method,epsi,
     +i_fkin,r,z,phi,
     +integral)
c-----------------------------------------------------------------
c     calculates the matrix: double complex integral
c     for the relativistic electron plasma 
c     I(n)^=integral(0<=p_perp<=p_perp0)G_n(p_perp)
c     G_n={sum_k=1,2,G_nk(p_perp)
c     for the EC harmonic with number 'n'
c----------------------------------------------------------------- 
c     input
c       n        - EC harmonic number
c       nll      - N_parallel
c       np       - N_perpendicular
c       Y        = omega_ce/omega for the electron rest mass
c       theta    = mc**2/T
c       ires     =0 no resonance,=1 ellipse,=2 parabola,=3 hyperbole 
c       n_relt_intgr - the number of points for the integration over p_perp
c       p_perp0  - max value of the perpendicular momentum divided by mc
c                  on the resonanse ellipse
c       i_fkin   =0 the usage of the analytical relativistic distributin
c                =1 the usage of the numerical 3D distribution from diskf 
c                   written be CQL3D code
c       r        - the major radius (it is used for i_fkin=1)
c       z        - the vertical coordinate  (it is used for i_fkin=1)
c       phi        - toroidal angle (it is used for i_fkin=1)
c----------------------------------------------------------
c     i_resonance_curve_interation_method
c
c     i_resonance_curve_integration_method=1 angle integration
c     now it works for ellipse case only                             c
c
c     i_resonance_curve_integration_method=2 rectangle formula
c     for p_perp integration
c     
c     i_resonance_curve_integration_method=3 trapezoidal formula
c     for p_perp integration
c
c     i_resonance_curve_integration_method=4 !adaptive Simpson 
c     for p_perp integration
c
c       epsi is the accuracy used in integration by adaptive Simpson
c------------------------------------------------------------ 
c     
c     output
c       integral(3,3) double complex integral from G_n over 0<=p_perp=<p_perpmax      
c----------------------------------------------------------------    
c      The integration method is specified by the variable
c      i_resonance_curve_integration_method.
c      Now this variable is set inside this subroutine.
c----------------------------------------------------------------

      implicit none
c     input
      double precision nll,np,Y,theta,p_perp0,epsi
      integer n,ires,n_relt_intgr,i_fkin,
     &i_resonance_curve_integration_method

      double precision r,z,phi
      double precision vnormloc,massloc ! cm/sec, g
      COMMON /dskin1/vnormloc,massloc


c     output
      double complex integral(3,3)
c     local
cSm060725
c      double precision eps,p_permax,h,p_perp,p,p_t,jmax,clight
      double precision eps,p_permax,h,p_perp,p,p_t,clight
      double complex i,g_n(3,3)
cSm060725
      integer j,jmax
      double precision vmax_d_vt 
ctest begin`
      double precision t_intgr,p_int
ctest end
c-----for integration along ellipse by angle
      double precision rme,rtem0,vper,xint,cs,sn,thet1,thet2,p_par,
     & p_par_min,p_par_pl,p_perp_min,p_perp_pl,v0dc,vmax1dc,vpar0dc,
     +vpar1dc,vpar2dc,dth,vt0,cdvt,tem0,dvper,vmaxdt
     
c     external
c     ec_cond, zeroK, g_n 
      double precision temperho
 
      write(*,*)'relt_tens.f sub intgr_rl:'
c      write(*,*)'i_resonance_curve_integration_method',
c     & i_resonance_curve_integration_method

      i = ( 0.0d0,1.0d0)        !imaginary number     
      jmax=n_relt_intgr
      

      vmax_d_vt=10.d0

      call p_perp_max_calc(i_fkin,theta,n,Y,nll,vmax_d_vt,
     &vnormloc,p_permax,ires)

      write(*,*)'in intgr_rl after p_perp_max_calc p_permax,ires',
     &p_permax,ires

      if(ires.eq.4) then
c-------the resonance curve is outside the grid
        call zeroK(integral)
        goto 10
      else
         h=p_permax/(n_relt_intgr)   ! step of integration over p_perp
      endif

      write(*,*)'in intgr_rl h',h

c-----calculations of the integrals over p_perp
      call zeroK(integral)
      write(*,*)'in intgr_rl set integral=zero integral',integral

c----------------------------------------------------------
c     i_resonance_curve_interation_method
c
c     i_resonance_curve_integration_method=1 angle integration
c     now it works for ellipse case only                             c
c     i_resonance_curve_integration_method=2 rectangle formula
c     for p_perp integration
c     
c     i_resonance_curve_integration_method=3 trapezoidal formula
c     for p_perp integration
c     
c------------------------------------------------------------ 
c      i_resonance_curve_integration_method=1
c      i_resonance_curve_integration_method=2
c      i_resonance_curve_integration_method=3 !trapezoidal 
c      i_resonance_curve_integration_method=4 !adaptive Simpson 
     
      goto (1,2,3,4) i_resonance_curve_integration_method

 1    continue
cSm060327       
      if(dabs(nll).ge.1.d0) goto 3 !to trapezoidal formula 
c--------------------------------------------------------------------
c     integration along the resonance curve using angle (thet)
c     thet1<thet<thet2  
c     begin
c     This case now works for ellipce case only
c     |N_parallel| <1
c--------------------------------------------------------------------
      rme=9.1094d-28 
ccc      call get_rtem0_from_one(rtem0)
      tem0=temperho(0.d0,1)*rtem0
      vt0=dsqrt(2.d0*tem0*1.6022d-9/rme) !(cm/sec)the  central
                                         ! thermal velocity
      clight=2.99792458d10
      cdvt=clight/vt0
      vmaxdt=dsqrt(rtem0)
ccc      call ec_condh(n,Y,nll,vmaxdt/cdvt,ires,v0dc,vmax1dc,vpar0dc,
ccc     +vpar1dc,vpar2dc,thet1,thet2)
     
      call zeroK(integral)     
c-----Subdivide theta range of integration
      dth=(thet2-thet1)/(jmax-1)
      do j=1,jmax
         xint=thet1+(j-1)*dth
         cs=dcos(xint)
         sn=dsin(xint)
         p_par=(vpar0dc-v0dc*cs)       !vper/c
         p_perp=vmax1dc*sn             !vpar/c
         if(p_perp.lt.1.d-12) p_perp=1.d-3
        
ccc         call g_n_calc_theta(p_perp,p_par,y,nll,np,theta,n,i_fkin,
ccc     &   r,z,phi,g_n)
         p_int=1.d0
        
         if((j.ne.jmax).and.(j.ne.1)) then         
           sn=dsin(xint+dth)
           p_perp_pl=vmax1dc*sn             !vper/c
           sn=dsin(xint-dth)
           p_perp_min=vmax1dc*sn            !vper/c  
           h=0.5d0*(p_perp_pl-p_perp_min)
         else
           if(j.eq.1) then
             sn=dsin(xint+dth)
             p_perp_pl=vmax1dc*sn             !vper/c
             h=0.5d0*(p_perp_pl-p_perp)
           else
             !j=jmax
             sn=dsin(xint-dth)
             p_perp_min=vmax1dc*sn             !vper/c
             h=0.5d0*(p_perp-p_perp_min)
           endif
         endif
                 
         integral(1,1)=integral(1,1)+h*g_n(1,1)
         integral(1,2)=integral(1,2)+h*g_n(1,2)
         integral(1,3)=integral(1,3)+h*g_n(1,3)
         integral(2,2)=integral(2,2)+h*g_n(2,2)
         integral(2,3)=integral(2,3)+h*g_n(2,3)
         integral(3,3)=integral(3,3)+h*g_n(3,3)
      enddo 
      goto 20
c--------------------------------------------------------------------
c     integration along the resonance curve using angle (thet)
c     thet1<thet<thet2  
c     end
c--------------------------------------------------------------------

 2    continue
c-------------------------------------------------------------------
c     integration along the resonance curve using rectangle
c     formula
c     begin
c------------------------------------------------------------------- 
      do j=1,jmax
        p_perp=h*(j-0.5)

        call g_n_calc(p_perp,y,nll,np,theta,n,i_fkin,r,z,phi,g_n)
        
        integral(1,1)=integral(1,1)+g_n(1,1)
        integral(1,2)=integral(1,2)+g_n(1,2)
        integral(1,3)=integral(1,3)+g_n(1,3)
        integral(2,2)=integral(2,2)+g_n(2,2)
        integral(2,3)=integral(2,3)+g_n(2,3)
        integral(3,3)=integral(3,3)+g_n(3,3)
      enddo 

      integral(1,1)=integral(1,1)*h
      integral(1,2)=integral(1,2)*h
      integral(1,3)=integral(1,3)*h
      integral(2,2)=integral(2,2)*h
      integral(2,3)=integral(2,3)*h
      integral(3,3)=integral(3,3)*h  
      goto 20
c ------------------------------------------------------------------
c     integration along the resonance curve using rectangle
c     formula
c     end
c-------------------------------------------------------------------

 3    continue
cSm060306
c -------------------------------------------------------------------
c     integration along the resonance curve using trapezoidal
c     formula
c     begin
c--------------------------------------------------------------------
      do j=0,jmax
         p_int=1.d0
         if((j.eq.1).or.(j.eq.jmax)) p_int=0.5d0
         p_perp=h*j
         if(p_perp.lt.1.d-3) p_perp=1.d-3
         if(j.eq.jmax) p_perp=p_perp-1.d-3
         
         call g_n_calc(p_perp,y,nll,np,theta,n,i_fkin,r,z,phi,g_n)
                  
         integral(1,1)=integral(1,1)+p_int*g_n(1,1)
         integral(1,2)=integral(1,2)+p_int*g_n(1,2)
         integral(1,3)=integral(1,3)+p_int*g_n(1,3)
         integral(2,2)=integral(2,2)+p_int*g_n(2,2)
         integral(2,3)=integral(2,3)+p_int*g_n(2,3)
         integral(3,3)=integral(3,3)+p_int*g_n(3,3)
      enddo 

      integral(1,1)=integral(1,1)*h
      integral(1,2)=integral(1,2)*h
      integral(1,3)=integral(1,3)*h
      integral(2,2)=integral(2,2)*h
      integral(2,3)=integral(2,3)*h
      integral(3,3)=integral(3,3)*h

      goto 20
c-------------------------------------------------------------------
c     end integration along the resonance curve using trapezoidal
c     formula
c     end
c-------------------------------------------------------------------

 4    continue
cSm070417
c -------------------------------------------------------------------
c     integration along the resonance curve using 
c     the adaptive Simpson function
c     begin
c--------------------------------------------------------------------
      write(*,*)'before calc_integral_array_by_adaptive_simpson'

      call calc_absorption_integral_array_by_adaptive_simpson
     &(y,nll,np,theta,n,i_fkin,r,z,phi,p_permax,epsi,integral)
      write(*,*)'in intgr_rl after'
      write(*,*)'calc_absorption_integral_array_by_adaptive_simpson'
      write(*,*)'integral',integral

c----------------------------------------------------------------------                  
 20   continue

      integral(2,1)=-integral(2,1)
      integral(3,1)= integral(3,1)
      integral(3,2)=-integral(3,2)
           
 10   continue

      return
      end   
  

      subroutine root_res(p_perp,nll,n,Y,kmax,p_par_rl)
c     calculates the roots p_par_k(p_perp) of the relativistic
c     resonance condition gamma=N_par*p_par+nY
c     and kmax- the number of roots  k=0,1,..kmax  
c     input
c       p_perp      - perpendicular momentum divided by (mc)
c       nll         - N_parallel
c       Y           - omega_ce/omega for the electron rest mass
c       n           - EC harmonic number    
c     output
c       kmax        - total number of roots
c       p_par_rl(2) - roots {parallel momentum divided by (mc)}
c-----------------------------------------------------------
      implicit none      
c     input
      double precision p_perp,nll,Y
      integer n
c     output
      double precision p_par_rl(2)
      integer kmax
c     local 
      double precision p_perp0,nlls,p_perps,ny,nys,det,sqrdet
      integer ires
c     external
c     ec_cond

      call ec_cond(n,Y,nll,ires,p_perp0)
c      write(*,*)'root_res ires,p_perp0',ires,p_perp0
      if(ires.eq.0) then
c-------no roots
        kmax=0
        goto 10
      endif

      nlls=nll*nll
      p_perps=p_perp*p_perp
      ny=n*Y
      nys=ny*ny

      if(ires.eq.1) then
c-------resonance ellipse
c        write(*,*)'root_res in ires=1 p_perp,p_perp0',p_perp,p_perp0
        if(p_perp.lt.p_perp0) then
c---------two roots
          kmax=2
          det=nlls+nys-1.d0-(1.d0-nlls)*p_perps

          if(det.lt.0.d0)then
             write(*,*)'in root_res det<0'
             ires=0
             kmax=0
             goto 10
          endif

          sqrdet=dsqrt(det)                     
          p_par_rl(1)=(nll*ny+sqrdet)/(1.d0-nlls)
          p_par_rl(2)=(nll*ny-sqrdet)/(1.d0-nlls)
c          write(*,*)'root_res in ires=1 kmax',kmax
        else
          if (p_perp.eq.p_perp0) then
c------------one root
             kmax=1
ctest begin
c             det=nlls+nys-1.d0-(1.d0-nlls)*p_perps
c             write(*,*)'on root_res p_perp,det',p_perp,det  
c             det=nlls+nys-1.d0-(1.d0-nlls)*p_perp0**2
c             write(*,*)'on root_res p_perp0,det(p_perp0)',p_perp0,det         
ctest end           
             p_par_rl(1)=nll*ny/(1.d0-nlls)
          else
c------------no roots
             kmax=0
          endif             
        endif
c        write(*,*)'root_res 2 in ires=1 kmax',kmax
        goto 10
      endif !end ellipse

      if(ires.eq.2) then
c-------resonance parabola, one root
        kmax=1               
        p_par_rl(1)=(p_perps+1.d0-nys)/(2.d0*nll*ny)
        goto 10
      endif !end parabola

      if(ires.eq.3) then
c-------resonance hyperbole, one root   
        kmax=1
        det=nlls+nys-1.d0-(1.d0-nlls)*p_perps
        sqrdet=dsqrt(det)
        if(nll.lt.-1.d0)then 
c---------p||_1=>-infinity                     
          p_par_rl(1)=(nll*ny+sqrdet)/(1.d0-nlls)
        else
c---------nll.gt.1 and p||_2=>+infinity
          p_par_rl(1)=(nll*ny-sqrdet)/(1.d0-nlls)
        endif
      endif
     
 10   continue
c      write(*,*)'root_res after 10 kmax',kmax
      return
      end
      


      subroutine s_calcn(n,p_perp,p_par,np,y,sn)
c     calculates tensor sn=p_perp*S(n)
c     input
c      n      - the number of the given EC harmonic
c      p_perp,p_par are the components of momentum divided by mc
c      np     =N_perpendicular
c      y      =algebraic omega_ce/omega for the electron rest mass
c     output
c      sn(3,3)- double complex  sn=p_perp*S(n)

      implicit none
c     input
      integer n
      double precision p_perp,p_par,np,y
c     output  
      double complex sn(3,3)
c     external
c     besj,DBESJ
     
c     locals ATTENTION!!!: the local variables are real and complex

c      complex i    
c      real  bj,bj1,bj_prim,b,d,b_abs
c      real p_perps,p_perpar,bjs
c      integer n_abs,ier,k
c      complex snc(3,3)

      double  complex i    
      double precision  bj,bj1,bj_prim,b,b_abs,bj_prim1
      double precision  p_perps,p_perpar,bjs
      integer n_abs,ier,k
      double complex snc(3,3)
      integer nz ! should be =0, If nz=1 the bessel function(dBESJ) will be J=0

      i=dcmplx(0.0d0,1.0d0)

c-----Calculation of the Bessel function bj= J_n(b) and 
c     its derivative bj_prim=dJ_n(b)/db.
c-----d is the relative error (input),ier -error switch,ier=0 is OK

c      write(*,*)'s_calncl n,p_perp,p_par,np,y',n,p_perp,p_par,np,y
        
      n_abs=iabs(n)            ! the order of the Bessel function

      k=1                      ! coefficient =1 for J_n with not-negative n      
      if(n.lt.0) k=(-1)**n_abs ! coefficient for Bessel function J_n with n<0

      b=np*p_perp/y            ! Bessel function argument

      b_abs=dabs(b)

c      write(*,*)'s_calncl b_abs',b_abs
cSAP080104
      goto 10 

cSm060318       
      call DBESJ(b_abs,dble(n_abs),1,bj,nz)
c      write(*,*)'relat_tens old n_abs,bj',n_abs,bj
c      call bes_calc(b_abs,n_abs,bj,bj_prim)
c      write(*,*)'relat_tens new n_abs,bj',n_abs,bj
      if(b.lt.0.d0)  bj=bj*(-1)**n_abs   !negative argument b for J_(n_abs)(b)
 
cSm060318
      call DBESJ(b_abs,dble(n_abs+1),1,bj1,nz)
c      write(*,*)'relat_tens,old n_abs+1 bj1',bj1
      if(b.lt.0.d0) bj1=bj1*(-1)**(n_abs+1) !negative argument b for J_(n_abs+1)(b)      
     
      if(n.eq.0) then 
        bj_prim=-bj1
      else
c-------n_abs.ne.0 
c        write(*,*)'relat_tens before b=0 b',b
        if(b.eq.0.d0)then 
          if(n_abs.eq.1)then
            bj_prim=0.5d0
          else
c-----------n_abs.ge.2
            bj_prim=0.0  d0   
          endif
        else 
c---------argument b.ne.0
          bj_prim=-bj1+n_abs*bj/b
        endif
c         write(*,*)'relat_tens bj_prim=',bj_prim
      endif
    
      bj=k*bj               ! for negative n
      bj_prim=k*bj_prim     ! for negative n

c      write(*,*)'in s_n  old bj,bj_prim',bj,bj_prim
 
10    continue
cSAP080104
      call bes_calc(b,n,bj,bj_prim)
c      write(*,*)'relat_tens new n_abs,bj,bj_prim',n_abs,bj,bj_prim

c-----calculaion sn=p_perp*s
      p_perps=p_perp*p_perp
      p_perpar=p_perp*p_par
      bjs=bj*bj

c-----complex 
c      snc(1,1)=p_perps*(n*bj/b)**2
cc      snc(1,2)=-i*p_perps*(n*bj*bj_prim/b)
c      snc(1,2)=i*p_perps*(n*bj*bj_prim/b)
c      snc(1,3)=p_perpar*(n*bjs/b)

cSAP080105 using b=np*p_perp/y  
      snc(1,1)=(n*bj*y/np)**2
      snc(1,2)=i*p_perp*n*bj*bj_prim*y/np
      snc(1,3)=p_par*n*bjs*y/np


      snc(2,2)=p_perps*(bj_prim*bj_prim)
c      snc(2,3)=i*p_perpar*(bj*bj_prim)
      snc(2,3)=-i*p_perpar*(bj*bj_prim)
      snc(3,3)=p_par*p_par*(bjs)
      
c-----double complex sn(3,3)   
      sn(1,1)=snc(1,1)
      sn(1,2)=snc(1,2)
      sn(1,3)=snc(1,3)
      sn(2,2)=snc(2,2)
      sn(2,3)=snc(2,3)
      sn(3,3)=snc(3,3)

      sn(2,1)=-sn(1,2)
      sn(3,1)= sn(1,3)
      sn(3,2)=-sn(2,3)

c      write(*,*)'relat_tens.f in s_calcn() sn',sn
            
      return
      end
 
      subroutine zeroK(K)
c     puts zero to K(i,j)
      implicit none
c     output
      double complex K(3,3)
c     local
      integer i,j
      double complex zero 
      zero=dcmplx(0.d0,0.d0)
      do i=1,3
        do j=1,3
          K(i,j)=zero
        enddo
      enddo
      return
      end

      double precision function u_n(p_perp, p_par, wcew, nll,
     +theta, n, i_fkin)
c     calculates the function u_n=U_n(np=p_per,p_perp=p_parallel)=(1/gamma)*
c     (n*Y*df/dp^_perpendicular+N_parallel*p^_perpendicular*df/dp^_parallel)
c     u_n is used in under integral complex marix function G_nk(3,3)
c     Here p^=p/mc
c     input
c       p_perp  - perpendicular momentum divided by (mc)
c       p_par     parallel momentum/mc
c       wcew    = omega_ce/omega for the electron rest mass    
c       nll     - N_parallel 
c       theta   - mc**2/T
c       n       - n number of the given resonance harmonic
c       i_fkin =0 the usage of the analytical relativistic Maxwellian distributin
c-------------------------------------------------------
      implicit none
c     input        
      double precision p_perp,p_par,wcew,nll,theta
      integer n,i_fkin
c     external besk2as
c     local
      double precision gamma, bk2,f_maxw,ps,p,ps_max,eps,pi,d_maxw_dp,
     + m_e,clight,energy,psi,rho,pitch,pitch0,bmin,btotal,
     + c_pitch,s_pitch,s_pitch0,c_pitch0,
     + dptc0dpt,                          !d(pitch0)/d(pitch)
     + dens 
c     distributin function and its derivatives
      double precision fdist0,dfdx,dfdpitch,dfdpitc0,dfdp
      integer initial     
      pi=4*datan(1.d0)
      ps= p_perp*p_perp+p_par*p_par
      p=dsqrt(ps)
      gamma= dsqrt(1.d0+ps)
      u_n=0.d0 ! to initialize

      if (i_fkin.eq.0) then ! it is always 0 now
c--------usage of the analytical relativistic Maxwellian distribution

c        calculation the Mackdonalds function bk2=K_2(theta)*EXP(theta)	
         call besk2as(theta,bk2)
  
cyup         eps=1.d-9 ! the min value of Maxwell exponent
cyup         ps_max=(1.d0-dlog(eps)/theta)**2-1.d0

         f_maxw=theta*dexp(theta*(1.d0-gamma))/(4.d0*pi*bk2) 
cyup ?         u_n=-p_perp*theta*f_maxw*(n*wcew+nll*p_par)/(gamma*gamma)

         f_maxw=theta*dexp(theta*(1.d0-gamma))/(4.d0*pi*bk2)
         d_maxw_dp=p*theta*f_maxw/gamma

         u_n=p_perp*d_maxw_dp*(n*wcew+nll*p_par)/gamma/p  
         
      else
         stop 'u_n is setup for i_fkin=0 only (i_diskf=0)' 
      endif
      
      return
      end    

c-----------------




      subroutine root_res_test(p_perp,nll,n,Y,kmax,p_par_rl,
     &p_0,del_p)
c     calculates the roots p_par_k(p_perp) of the relativistic
c     resonance condition gamma=N_par*p_par+nY
c     and kmax- the number of roots  k=0,1,..kmax  
c     input
c       p_perp      - perpendicular momentum divided by (mc)
c       nll         - N_parallel
c       Y           - omega_ce/omega for the electron rest mass
c       n           - EC harmonic number    
c     output
c       kmax        - total number of roots
c       p_par_rl(2) - roots {parallel momentum divided by (mc)}
c       p_0         -ellipse center for ires=1 case
c       del_p       -shift from the ellipse center for  ires=1 case
c-----------------------------------------------------------
      implicit none      
c     input
      double precision p_perp,nll,Y
      integer n
c     output
      double precision p_par_rl(2),p_0,del_p 
       
      integer kmax
      
c     local 
      double precision p_perp0,nlls,p_perps,ny,nys,det,sqrdet
      integer ires
c     external
c     ec_cond

      call ec_cond(n,Y,nll,ires,p_perp0)
c      write(*,*)'root_res ires,p_perp0',ires,p_perp0
      if(ires.eq.0) then
c-------no roots
        kmax=0
        goto 10
      endif

      nlls=nll*nll
      p_perps=p_perp*p_perp
      ny=n*Y
      nys=ny*ny

      p_0=nll*ny/(1.d0-nlls)
      del_p=dsqrt(nlls+nys-1.d0)/(nlls-1.d0)

      if(ires.eq.1) then
c-------resonance ellipse
        p_0=nll*ny/(1.d0-nlls)
        del_p=dsqrt(nlls+nys-1.d0)/(1.d0-nlls)
        if(p_perp.lt.p_perp0) then
c---------two roots
          kmax=2
          det=nlls+nys-1.d0-(1.d0-nlls)*p_perps

          if(det.lt.0.d0)then
             write(*,*)'in root_res det<0'
             ires=0
             kmax=0
             goto 10
          endif

          sqrdet=dsqrt(det)                     
          p_par_rl(1)=(nll*ny+sqrdet)/(1.d0-nlls)
          p_par_rl(2)=(nll*ny-sqrdet)/(1.d0-nlls)

        else
          if (p_perp.eq.p_perp0) then
c------------one root
             kmax=1
ctest begin
c             det=nlls+nys-1.d0-(1.d0-nlls)*p_perps
c             write(*,*)'on root_res p_perp,det',p_perp,det  
c             det=nlls+nys-1.d0-(1.d0-nlls)*p_perp0**2
c             write(*,*)'on root_res p_perp0,det(p_perp0)',p_perp0,det         
ctest end           
             p_par_rl(1)=nll*ny/(1.d0-nlls)
           
          else
c------------no roots
             kmax=0
          endif             
        endif
        goto 10
      endif !end ellipse

      if(ires.eq.2) then
c-------resonance parabola, one root
        kmax=1               
        p_par_rl(1)=(p_perps+1.d0-nys)/(2.d0*nll*ny)
        goto 10
      endif !end parabola

      if(ires.eq.3) then
c-------resonance hyperbole, one root 
        kmax=1
        det=nlls+nys-1.d0-(1.d0-nlls)*p_perps
        sqrdet=dsqrt(det)
        if(nll.lt.-1.d0)then 
c---------p||_1=>-infinity                     
          p_par_rl(1)=(nll*ny+sqrdet)/(1.d0-nlls)
        else
c---------nll.gt.1 and p||_2=>+infinity
          p_par_rl(1)=(nll*ny-sqrdet)/(1.d0-nlls)
        endif
      endif

 10   continue
      return
      end
      
      subroutine p_perp_max_calc(i_fkin,theta,n,Y,npar,vmax_d_vt,
     &vnorm,p_perp_max,ires)
c-------------------------------------------------------------
c     calculate maximal value of p_perp =p_perp_max
c     and ires.If ires=4 then the resonance curve is outside the
c     grid for i_fkin=1 case
c-------------------------------------------------------------
      implicit none
c-----input
      double precision npar, ! parallel refractive index
     &theta,                 ! = mc**2/T
     &Y,                     ! omega_ce/omega for the electron rest mass
     &vmax_d_vt,             ! maximal [momentum/pc] devided by the thermal
                             ! momentum  used in the integration for analytical
                             ! Maxwellian distribution (i_fkin=1 case)
     &vnorm                  ! maximal momentum/m at the grid (momentum-per-mass) 
		             ! [cms/sec]
      integer n,  ! harmonic number
     &i_fkin      !  =0 the usage of the analytic relativistic distribution
                  !  =1 the usage of the numerical 3D distribution from diskf 
                  !  written be CQL3D code 
                  !  or created inside genray in forrer.f
 
c-----output
      double precision p_perp_max ! max p_perp in integration along
                                  ! the resonace curve in [p/mc]
      integer ires !ires =4 outside the greed 
c-----locals            
      double precision p_perp0,gammax,clight,ny,vnorm_d_c,
     &p_par_bound,p_perp_bound,
     & p_par_0,p_par_dl,p_par_1,p_par_2,p_par_perp0
     
      clight=2.99792458d10          !cm/sec           
      call ec_cond(n,Y,npar,ires,p_perp0)
cyup      write(*,*)'in sub p_perp_max_calc after ec_cond'
cyup      write(*,*)'i_fkin,ires,p_perp0',i_fkin,ires,p_perp0

      if (i_fkin.eq.0) then
c--------analytical Maxwellian distribution
         p_perp_max=vmax_d_vt*dsqrt((1.d0/theta+1.d0)**2-1.d0) ! vmax_d_vt*
                                                               ! thermal momentum
                                                               ! [p/mc]
cyup         write(*,*)'in sub p_perp_max_calc 1 p_perp_max',p_perp_max
         if (ires.eq.1) p_perp_max=dmin1(p_perp_max,p_perp0)
cyup         write(*,*)'p_perp_max=dmin1(p_perp_max,p_perp0)',p_perp_max
      elseif (i_fkin.eq.1) then
c--------distribution was set on the grid
         vnorm_d_c=vnorm/clight      
         ny=n*y    
         gammax=dsqrt(1.d0+vnorm_d_c**2) 
         p_par_bound=(gammax-ny)/npar    ![p_parallel/mc]                    

c         write(*,*)'vnorm_d_c,p_par_bound,n,y',
c     &              vnorm_d_c,p_par_bound,n,y

         if(ires.eq.1) then !ellipse case, abs(npar)<1
            call ec_cond(n,Y,npar,ires,p_perp0)

            p_par_0=npar*ny/(1.d0-npar**2)      ! ellipse center p_parallel/mc coordinate
            p_par_dl=dsqrt(npar**2+(ny)**2-1.d0)/(1.d0-npar**2)                                  
            p_par_1=p_par_0-dabs(p_par_dl)      ! ellipse points p_parallel/mc at p_perp=0
            p_par_2=p_par_0+dabs(p_par_dl)
c            write(*,*)'p_par_0,p_perp0,p_par_dl,p_par_1,p_par_2',
c     &                 p_par_0,p_perp0,p_par_dl,p_par_1,p_par_2

            if (vnorm_d_c.ge.dabs(p_par_bound)) then
c--------------resonance ellipse has the intersection with the circle vnorm
               p_perp_bound=dsqrt(vnorm_d_c**2-p_par_bound**2)
c               write(*,*)'p_par_bound,p_perp_bound',
c     &                    p_par_bound,p_perp_bound             
               if (p_par_0.ge.0.d0)then 
c-----------------the ellipse is shifted to the right side
                  if (p_par_bound.le.p_par_0) then
c--------------------max P_perp in the intersection point
c                    (p_par_bound,p_perp_bound) 
                     p_perp_max=p_perp_bound
                  else
c--------------------max P_perp is in the top ellipse point 
c                    (p_par_0,p_perp0)  
                     p_perp_max=p_perp0
                  endif
               else
                  !v_par_0.lt.0.d0 
c-----------------the ellipse is shifted to the left side
                  if (p_par_bound.gt.p_par_0) then
c--------------------max P_perp in the intersection point
c                    (p_par_bound,p_perp_bound) 
                     p_perp_max=p_perp_bound
                  else
c--------------------max P_perp is in the top ellipse point 
c                    (p_par_0,p_perp0)  
                     p_perp_max=p_perp0
                  endif
               endif
            else
c--------------resonance ellipse has not intersection with the circle vnorm
               if((dabs(p_par_1).lt.vnorm_d_c).and.
     &            (dabs(p_par_2).lt.vnorm_d_c)) then
c-----------------the resonance ellipse is inside the circle vnorm 
c                 max P_perp is in the pot ellipse point   
c                            (p_par_0,p_perp0)
                  p_perp_max=p_perp0
                else
c-----------------the resonance ellipse is outside the circle vnorm
c                 It means no resonanse inside the circle 
                  ires=4
                endif  
            endif
         endif !ires=1 

         if(ires.eq.2) then !resonance parabola, abs npar=1
             p_par_perp0=(1.d0-ny**2)/(2.d0*ny*npar) !p_par on the parabla at p_perp=0
             if(dabs(p_par_perp0).lt.vnorm_d_c) then 
               p_perp_max=dsqrt((vnorm_d_c)**2-p_par_bound**2) ![p_perp/mc]
             else
               ires=4 ! the resonance curve is outside the grid 
               goto 10
             endif
         endif !ires=2
    
         if(ires.eq.3) then !resonance hyperbole, abs(npar)>1
             p_par_0=npar*ny/(1.d0-npar**2)      ! the center point p_parallel/mc coordinate
                                                 ! betwee0n the hyperbole branches
             p_par_dl=dsqrt(npar**2+(ny)**2-1.d0)/(1.d0-npar**2)
             if(npar.gt.1.d0) then
               !the resonance curve is
               !the positive hyperbole branche: p_par_2=p_par_0-p_par_dl < p_par < +infinity
               p_par_2=p_par_0-p_par_dl
               if (dabs(p_par_2).lt.vnorm_d_c) then
                  p_perp_max=dsqrt(vnorm_d_c**2-p_par_bound**2) ![p_perp/mc]
               else
                  ires=4 ! the resonanse curve is outside the grid 
                  goto 10
               endif
             else
               ! npar <-1 and the resonance curve is
               ! the negatibe hyperbole branche: -infinity < p_par < p_par_1=p_par_0+p_par_dl
               p_par_1=p_par_0+p_par_dl
               if (dabs(p_par_1).lt.vnorm_d_c) then
                  p_perp_max=dsqrt(vnorm_d_c**2-p_par_bound**2) ![p_perp/mc]
               else
                  ires=4 ! the resonanse curve is outside the grid 
                  goto 10
               endif
             endif 
         endif !ires=3 
      endif !i_fkin
 
 10   continue

      return
      end  

c============================================================================
c     Adaptive Simpson functions.
C
c     Modificated method for vector function.
c     
c     Original description and functions are from:
c
C PAGE 187-190: NUMERICAL MATHEMATICS AND COMPUTING, CHENEY/KINCAID, 1985
C
C FILE: SIMP.FOR
C
C ADAPTIVE SCHEME FOR SIMPSON'S RULE (SIMP,ASMP,PUSH,POP,FCN)
c============================================================================

      subroutine fcn_absorption_vector(p_perp,fcn)
c---------------------------------------------------------------------------
c     culculates under integral function for anti-hermitian dielectric tesor
c---------------------------------------------------------------------------

      implicit none
c      include 'param_kmax.i'
      integer k_max !max number of elements
      parameter (k_max=6) 
c----------------------------------------------------------------
c     input
c----------------------------------------------------------------- 
      real*8 p_perp !under integral function argument
c-----------------------------------------------------------------
c     output
c-----------------------------------------------------------------
      real*8 fcn(k_max) !array of under integral fuctions values
c-----------------------------------------------------------------
      integer n_l,i_fkin_l
      real*8 y_l,nll_l,np_l,theta_l,z_l,r_l,phi_l 
      
      common /fcn_input/y_l,nll_l,np_l,theta_l,n_l,i_fkin_l,
     &z_l,r_l,phi_l 
c-----------------------------------------------------------------
c     locals
c----------------------------------------------------------------
      double complex g_n(3,3) !comples under integral functions

      call g_n_calc(p_perp,y_l,nll_l,np_l,theta_l,n_l,i_fkin_l,
     &r_l,z_l,phi_l,g_n)

      FCN(1) =dreal(g_n(1,1)) 
      FCN(2) =dimag(g_n(1,2))
      FCN(3) =dreal(g_n(1,3))
      FCN(4) =dreal(g_n(2,2))
      FCN(5) =dimag(g_n(2,3))
      FCN(6) =dreal(g_n(3,3))

      return
      end


      subroutine PUSH(A,B,LEVEL,STACK,IST,I)
 
      implicit none
      integer level,ist,i
      real*8 a,b    
      REAL*8 STACK(IST,3)     
  
      IF(I .LT. IST) THEN   
        I = I+1 
        STACK(I,1) = A      
        STACK(I,2) = B
c        STACK(I,3) = dREAL(LEVEL)
         STACK(I,3) = dfloat(LEVEL)
      
        RETURN
      ELSE
        STOP 'STACK OVERFLOW IN PUSH' 
      END IF      
      END 
  
      subroutine POP(A,B,LEVEL,STACK,IST,I)
      implicit none
      real*8 a,b
      integeri,ist,level     
      REAL*8 STACK(IST,3)
     
      IF(I .GT. 0) THEN 
        A = STACK(I,1)
        B = STACK(I,2)
        LEVEL = INT(STACK(I,3))       
        I = I-1 
        RETURN
      ELSE
        STOP 'STACK UNDERFLOW IN POP' 
      END IF      
      END 
  
      subroutine ASMP(FCN_vector,k_max,fcn,
     &A,B,EPSI,LVMAX,STACK,IST,FSTACK,IFS,
     &SUM,IFLAG)
      implicit none

c-----input
      integer k_max !max number of element

      real*8 a,b,epsi,sum(k_max) 
      integer lvmax,level,itop,iflag,loop,ist,ifs,lmax,k,
     &i     
      REAL*8 STACK(IST,3),FSTACK(IFS,3) 
      real*8 fcn(k_max)
      EXTERNAL FCN_vector
      COMMON /BLK/ ITOP,LEVEL,LMAX    
    
      ITOP = 0
      LEVEL = 0   
      IFLAG = 0 
      do k=1,k_max
        SUM(k) = 0.0d0
      enddo 
      LMAX = 2**LVMAX 
      IF(IST .GE. LVMAX+1 .AND. IFS .GE. LMAX) THEN 
        CALL PUSH(A,B,LEVEL,STACK,IST,ITOP)     
        DO 2 LOOP=1,2*LMAX-1

c          write(*,*)'ASMP LOOP',LOOP

          IF(ITOP .EQ. 0) RETURN      
          CALL POP(A,B,LEVEL,STACK,IST,ITOP)    
          CALL SIMP(FCN_vector,k_max,fcn,
     &             A,B,EPSI,LVMAX,STACK,IST,FSTACK,IFS,
     &             SUM,IFLAG)
c          write(*,*)'loop,level,iflag,lvmax,lmax',
c     &    loop,level,iflag,lvmax,lmax
c          write(*,*)'ist,ifs',ist,ifs
c          write(*,*)'sum(k)',(sum(k),k=1,6)
c          do i=1,ist
c          write(*,*)'i',i,'stack(i,k)',(stack(i,k),k=1,3)
c          enddo
c          do i=1,ifs
c          write(*,*)'i',i,'fstack(i,k)',(fstack(i,k),k=1,3)
c          enddo

          

    2   CONTINUE  
      ELSE
        PRINT 3   
      END IF      
      RETURN
   3  FORMAT(//5X,'NOT ENOUGH WORKSPACE IN STACK OR FSTACK')
      END 
  
      subroutine SIMP(FCN_vector,k_max,fcn,
     &A,B,EPSI,LVMAX,STACK,IST,FSTACK,IFS,
     &SUM,IFLAG)
      implicit none
      integer k_max !max number of elements 

      real*8 a,b,epsi,sum(k_max),h,sum2(k_max),sum4(k_max),h4      
      integer itop,level,lmax,iflag,ifs,i,ist,lvmax,icondition,k
      real*8 fcn(k_max) 
      REAL*8 STACK(IST,3),FSTACK(IFS,3),X(0:4),F(0:4,k_max) 
      COMMON /BLK/ ITOP,LEVEL,LMAX    
      H = B - A 
      H4 = H/4.0d0  
           
      DO 2 I = 0,4
c        X(I) = A + dREAL(I)*H4
        X(I) = A + dfloat(I)*H4
        call fcn_vector(X(I),fcn)
c        write(*,*)'i,X(i)',i,X(i)
        do k=1,k_max 
           F(I,k) = FCN(k)
c           write(*,*)'k,F(I,k)',k,F(I,k)
        enddo  
   2  CONTINUE

      do k=1,k_max    
        SUM2(k) = H*(F(0,k) + 4.0d0*F(2,k) + F(4,k))/6.0d0 
        SUM4(k) = H*(F(0,k) + 4.0d0*F(1,k) + 2.0d0*F(2,k)+
     &          4.0d0*F(3,k)+F(4,k))/12.d0
      enddo

c      icondition=0 
c      do k=1,k_max
c        IF(dABS(SUM4(k) - SUM2(k)).LE.15.0d0*EPSI/dREAL(2**LEVEL))THEN
c        write(*,*)'k',k
c        write(*,*)'dABS(SUM4(k) - SUM2(k))',dABS(SUM4(k) - SUM2(k))
c        write(*,*)'15.0d0*EPSI/dfloat(2**LEVEL)',
cc     &             15.0d0*EPSI/dfloat(2**LEVEL)
c        IF(dABS(SUM4(k) - SUM2(k)).LE.15.0d0*EPSI/dfloat(2**LEVEL))THEN
c           icondition=icondition+0
c        else
c           icondition=icondition+1
c        endif
c        write(*,*)'0 k,icondition',k,icondition
c      enddo

      icondition=0 
      do k=1,k_max
c        IF(dABS(SUM4(k) - SUM2(k)).LE.15.0d0*EPSI/dREAL(2**LEVEL))THEN
c        write(*,*)'k',k
c        write(*,*)'dABS(SUM4(k) - SUM2(k))',dABS(SUM4(k) - SUM2(k))
c        write(*,*)'15.0d0*EPSI/dfloat(2**LEVEL)',
c     &             15.0d0*EPSI/dfloat(2**LEVEL)
        IF(dABS(SUM4(k) - SUM2(k)).GT.15.0d0*EPSI/dfloat(2**LEVEL))THEN
           icondition=1
        endif
c        write(*,*)'1 k,icondition',k,icondition
      enddo
      IF (LEVEL.LT.3) icondition=1
      

c      IF(dABS(SUM4(k) - SUM2(k)) .LE. 15.0d0*EPSI/dREAL(2**LEVEL)) THEN
c      write(*,*)'icondition',icondition  
      if (icondition.eq.0) then       
        do k=1,k_max
          SUM(k) = SUM(k) + SUM4(k)
        enddo
      ELSE
        IF(LEVEL .LT. LVMAX) THEN     
          LEVEL = LEVEL + 1 
          CALL PUSH(X(2),X(4),LEVEL,STACK,IST,ITOP)       
          CALL PUSH(X(0),X(2),LEVEL,STACK,IST,ITOP)       
        ELSE      
          CALL PUSH(X(0),X(4),LEVEL,FSTACK,IFS,IFLAG)     
        END IF    
      END IF      
      RETURN
      END 

      subroutine calc_absorption_integral_array_by_adaptive_simpson
     &(y,nll,np,theta,n,i_fkin,r,z,phi,p_permax,epsi,integral)
c------------------------------------------------------------------------
c     calculate intergals integral(3,3) for anti-hermitian diectric tensor
c     for the EC harmonic with number 'n'
c-------------------------------------------------------------------------
       implicit none
c------------------------------------------------------------------------
c      input
c------------------------------------------------------------------------
      integer n,i_fkin
      real*8 y,nll,np,theta,z,r,phi,p_permax,epsi 
c     n        is EC harmonic number
c     i_fkin   =0 the usage of the analytical relativistic distributin
c              =1 the usage of the numerical 3D distribution from diskf 
c                   written be CQL3D code
c     nll      is N_parallel
c     theta    = mc**2/T
c     r        is the major radius (it is used for i_fkin=1)
c     z        is  the vertical coordinate  (it is used for i_fkin=1)
c     phi      is  the toroidal angle (it is used for i_fkin=1)
c     p_permax is  the maximal boundary of integration
c     epsi     is the accuracy used in integration by adaptive Simpson
c------------------------------------------------------------------------
c     output
c------------------------------------------------------------------------
      double complex integral(3,3)
c-------------------------------------------------------------------------
c     locals
c--------------------------------------------------------------------------
      integer k_max !max number of elements
      parameter (k_max=6)
      integer ist,ifs,k   
c      PARAMETER (IST=5,IFS=16)
c      PARAMETER (IST=9,IFS=256)  
c       PARAMETER (IST=13,IFS=10048)
c       PARAMETER (IST=16,IFS=100048)
c       PARAMETER (IST=17,IFS=200048)!070723-old
      PARAMETER (IST=20,IFS=200048)

      REAL*8 STACK(IST,3),FSTACK(IFS,3)
      real*8 a,b,sum(k_max)
      integer lvmax,iflag,i
      real*8 fcn(k_max)
    
      integer n_l,i_fkin_l
      real*8 y_l,nll_l,np_l,theta_l,z_l,r_l,phi_l 
      common /fcn_input/y_l,nll_l,np_l,theta_l,n_l,i_fkin_l,
     &z_l,r_l,phi_l
 
      EXTERNAL FCN_absorption_vector
c      DATA EPSI/5.0d-5/, LVMAX/4/
c      DATA EPSI/5.0d-5/, LVMAX/7/
c      DATA EPSI/5.0d-5/, LVMAX/15/
c      DATA EPSI/5.0d-5/, LVMAX/12/
c      DATA EPSI/5.0d-5/, LVMAX/16/
c       DATA EPSI/5.0d-6/, LVMAX/16/  !070722-new
c       DATA EPSI/1.0d-6/, LVMAX/16/  !070722-new

c       DATA EPSI/1.0d-3/, LVMAX/16/ !070722-old
c      DATA EPSI/1.0d-3/, LVMAX/19/
      DATA LVMAX/16/  !070722-new



c-----------------------------------------------------------------
c     set common /fcn_input/
c-----------------------------------------------------------------
      y_l=y
      nll_l=nll
      np_l=np
      theta_l=theta
      n_l=n
      i_fkin_l=i_fkin
      z_l=z
      r_l=r
      phi_l=phi
c-----------------------------------------------------------------
c     integration boundaries
c------------------------------------------------------------------
c      A = 0.0d0
      A = 1.d-5
c      A = 1.d-3

      B = p_permax-1.d-5
     
      CALL ASMP(FCN_absorption_vector,k_max,fcn,
     &A,B,EPSI,LVMAX,STACK,IST,FSTACK,IFS,
     &SUM,IFLAG) 

c      do k=1,k_max 
c        write(*,*)'k=',k 
c        PRINT 3,SUM(k)
c      enddo  
     

      IF(IFLAG .EQ. 0) THEN 
c$$$        PRINT 4   
      ELSE
        write(*,*)'          WITH BAD SUBINTERVALS:'   
        DO 2 I=1,IFLAG      
          write(*,*) FSTACK(I,1),FSTACK(I,2),INT(FSTACK(I,3))
    2   CONTINUE  
      END IF      
   3  FORMAT(//5X,'APPROXIMATE INTEGRAL =',E22.14/)   
   4  FORMAT(10X,'WITH NO BAD SUBINTERVALS')    
   5  FORMAT(10X,'WITH BAD SUBINTERVALS:')      
   6  FORMAT(10X,'[',F10.5,',',F10.5,']',2X,'LEVEL =',I5) 

c-------------------------------------------------------
      integral(1,1)=dcmplx(Sum(1),0.d0) 
      integral(1,2)=dcmplx(0.d0,Sum(2))
      integral(1,3)=dcmplx(Sum(3),0.d0) 
      integral(2,2)=dcmplx(Sum(4),0.d0)
      integral(2,3)=dcmplx(0.d0,Sum(5))
      integral(3,3)=dcmplx(Sum(6),0.d0) 

     
      return
      END 

 
