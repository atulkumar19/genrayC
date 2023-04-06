c        **********************   x   ***********************
c        *                        -                         *
c        * this function calculates the ratio of squares of *
c        * plasma frequency of species i and wave frequency *
c        * x=(omega_pl_i/omega)**2    
c        * small radius rho is in common one
c
c        * The minimal X value is set inside this function
c        * Xmin=1.d-6. If X<Xmin then X=Xmin                        
c        ****************************************************
c
c-------------------------------------------------------------------
c								    !
c        input parameters					    !
c      								    !
c      z, r phi - coordinates of the  point  where  the  ratio  is  !
c                 calculated.      		         	    !
c		  rho is in common one						    !
c-------------------------------------------------------------------
      double precision
     1function x(z,r,phi,i)
c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none
     
      include 'param.i'
      include 'one.i'
c-----input
      real*8 z,r,phi !space coordinates
c     rho is in common  /one/
      integer i !number of plasma species

c-----local 
      real*8 den

c-----external
      real*8 dense  
      
      if(ixyz.eq.1) then
        stop 'x(z,r,phi,i): should not be called for ixyz=1'
      endif
          
      den=dense(z,r,phi,i)
      if(den.lt.0.d0)then
         den=0.d0
      endif
      x=v(i)*den
      x=dmax1(x,1.d-6)      
      return
      end

c     **********************vardens***********************
c     *                        -                         *
c     * this function calculates density variations      *
c     ****************************************************
c
c-------------------------------------------------------------------
c								    !
c        input parameters					    !
c      								    !
c      z, r phi - coordinates of the  point  where  the  ratio  is  !
c                 calculated.      		         	    !
c		  rho is in common one						    !
c-------------------------------------------------------------------
      double precision
     1function vardens(z,r,phi)
cSAP090201
c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none
      include 'param.i'
      include 'one.i'
c-----------------------------------------------------
c     parameters of the density variations
c     are in one.i (var0,denn,denm,an,sigman)
c-----------------------------------------------------
c-----input
      real*8 z,r,phi !space coordinates
c-----locals  
      real*8 theta, !poloioidal angle
     &varrho,varphtht,varphi,vartheta
c-----externals
      real*8 thetapol 

      pi=4.d0*datan(1.d0)
ccc      varrho=var0*0.5d0*(1.d0+datan((rho-an)/sigman**2)/(0.5*pi))
ccc     1 *(rho-an)/(1.d0-an)
       varrho=0.d0
cc      write(*,*)'varrho',varrho
c-----------------------------------------------------
      theta=thetapol(z,r)
c      write(*,*)'in vardens ater thetapol theta',theta
      varphi=dcos(denn*phi)
      vartheta=1.d0+dcos(denm*theta)
      varphtht=varphi*vartheta
c      varphtht=dcos(denm*theta+denn*phi)
c-----------------------------------------------------
      vardens=varrho*varphtht
cc      write(*,*)'in vardens varrho',varrho
cc      write(*,*)'varphi,vartheta,vardens',varphi,vartheta,vardens
      return
      end

      double precision
     1function dvarddz(z,r,phi)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
c-----------------------------------------------------
c      varrho=var0
c      dvarrdr=0.d0
c-----------------------------------------------------
      psi_xr=psif(z,r)
      dro_dpsi=drhopsi(psi_xr)
c-----------------------------------------------------
cc      varrho=var0*dexp(-((1.d0-an/rho)/sigman)**2)
cc      dvarrdz=-2.d0*varrho*an*(1.d0-an/rho)/(rho*sigman)**2*
cc     1dro_dpsi*dpdzd
c-----------------------------------------------------
      pi=4.d0*datan(1.d0)
      var0=0.d0
      varrho=var0*0.5d0*(1.d0+datan((rho-an)/sigman**2)/(0.5*pi))
     1 *(rho-an)/(1.d0-an)
      dvarrdz=var0*(0.5d0*(1.d0+datan((rho-an)/sigman**2)/(0.5d0*pi))/
     1 (1.d0-an)+0.5d0*(rho-an)/(1.d0+((rho-an)/sigman**2)**2)/
     1 (0.5d0*pi*sigman**2)/(1.d0-an))*
     1 dro_dpsi*dpdzd
c-----------------------------------------------------

      theta=thetapol(z,r)
      varphi=dcos(denn*phi)
      vartheta=1.d0+dcos(denm*theta)
      varphtht=varphi*vartheta
      dvarthdz=-denm*dsin(denm*theta)*dthetadz(z,r)
      dvarptdz=varphi*dvarthdz
      dvarddz=varrho*dvarptdz+varphtht*dvarrdz
cc      write(*,*)'in dvarddz varrho',varrho,'dvarptdz',dvarptdz
cc      write(*,*)'in dvarddz varphtht',varphtht,'dvarrdz',dvarrdz
cc      write(*,*)'in dvarddz=',dvarddz
      return
      end

      double precision
     1function dvarddr(z,r,phi)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
c-----------------------------------------------------
c      varrho=var0
c      dvarrdr=0.d0
c-----------------------------------------------------
      psi_xr=psif(z,r)
      dro_dpsi=drhopsi(psi_xr)
c-----------------------------------------------------
cc      varrho=var0*dexp(-((1.d0-an/rho)/sigman)**2)
cc      dvarrdr=-2.d0*varrho*an*(1.d0-an/rho)/(rho*sigman)**2*
cc     1dro_dpsi*dpdrd
c-----------------------------------------------------
      pi=4.d0*datan(1.d0)
      var0=0.d0
      varrho=var0*0.5d0*(1.d0+datan((rho-an)/sigman**2)/(0.5*pi))
     1 *(rho-an)/(1.d0-an)
      dvarrdr=var0*(0.5d0*(1.d0+datan((rho-an)/sigman**2)/(0.5d0*pi))/
     1 (1.d0-an)+0.5d0*(rho-an)/(1.d0+((rho-an)/sigman**2)**2)/
     1 (0.5d0*pi*sigman**2)/(1.d0-an))*
     1 dro_dpsi*dpdrd
c-----------------------------------------------------
      theta=thetapol(z,r)
      varphi=dcos(denn*phi)
      vartheta=1.d0+dcos(denm*theta)
      varphtht=varphi*vartheta
      dvarthdr=-denm*dsin(denm*theta)*dthetadr(z,r)
      dvarptdr=varphi*dvarthdr
c      dvarddr=varrho*dvarptdr
      dvarddr=varrho*dvarptdr+varphtht*dvarrdr
cc      write(*,*)'in dvarddr varrho',varrho,'dvarptdr',dvarptdr
cc      write(*,*)'in dvarddr varphtht',varphtht,'dvarrdr',dvarrdr
cc      write(*,*)'in dvarddr=',dvarddr
      return
      end

      double precision
     1function dvarddph(z,r,phi)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
c------------------------------------------------
c      varrho=var0
c------------------------------------------------
c      varrho=var0*dexp(-((1.d0-an/rho)/sigman)**2)
c------------------------------------------------
      pi=4.d0*datan(1.d0)
      var0=0.d0
      varrho=var0*0.5d0*(1.d0+datan((rho-an)/sigman**2)/(0.5*pi))
     1 *(rho-an)/(1.d0-an)
c------------------------------------------------
      theta=thetapol(z,r)
      vartheta=1.d0+dcos(denm*theta)
      dvarphdp=-denn*dsin(denn*phi)
      dvarptdp=dvarphdp*vartheta
      dvarddph=varrho*dvarptdp
cc      write(*,*)'in dvarddph vartheta',vartheta,'dvarphdp',dvarphdp
cc      write(*,*)'in dvarddph dvarptdp',dvarptdp,'varrho',varrho
cc      write(*,*)'in dvarddph',dvarddph
      return
      end

c************  thetapol********
      double precision
     1function thetapol(z,r)
c     calculation of the poloidal angle theta (rad)
c     -pi.gt.theta.le.pi
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      double precision modr,cos_p
      zcomp=z-zma
      rcomp=r-rma
      modr=dsqrt(zcomp*zcomp+rcomp*rcomp)  !=|r|
      if(modr.lt.1.d-9) then
        theta=0.d0
      else
        cos_p=rcomp/modr
         
        if(zcomp.ge.0.d0) then
c         theta=dacos(rcomp/modr)
          if(dabs(cos_p).le.1.d0) then
             theta=dacos(cos_p)
          else
             if(cos_p.ge.0.d0) then
               theta=0.d0
             else
               theta=pi
             endif
          endif 
        else
c         theta=-dacos(rcomp/modr)
          if(dabs(cos_p).le.1.d0) then
             theta=-dacos(cos_p)
          else
             if(cos_p.ge.0.d0) then
               theta=0.d0
             else
               theta=-pi
             endif
          endif 
        endif
      endif
      thetapol=theta
      return
      end


c************ dthetadz********
      double precision
     1function dthetadz(z,r)
c     calculation of the derivative D(theta)/Dz
c     -pi.gt.theta.le.pi
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      double precision modr,modr2,modr3
      zcomp=z-zma
      rcomp=r-rma
      modr=dsqrt(zcomp*zcomp+rcomp*rcomp)  !=|r|
      modr2=modr*modr
      dthetadz=rcomp/modr2
      return
      end


c************ dthetadr********
      double precision
     1function dthetadr(z,r)
c     calculation of the derivative D(theta)/Dr
c     -pi.gt.theta.le.pi
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      double precision modr,modr2,modr3
      zcomp=z-zma
      rcomp=r-rma
      modr=dsqrt(zcomp*zcomp+rcomp*rcomp)  !=|r|
      modr2=modr*modr
      dthetadr=-zcomp/modr2
      return
      end


