c        ********************** dddrz1_xyz ******************
c        *                      -----                       *
c        * this subroutine calculates the derivatives
c        * from Hamiltonian  				    *
c        * for 6 ray equations                              *
c        * The Hamiltonian derivatives  are calculated	    *
c        * analytically(idif=1) or numerically (idif=2)	    *
c        * This subroutine is used for the Hamiltonian
c        * correction precedure
c        ****************************************************
c-------------------------------------------------------------------
c								    !
c        input parameters					    !
c      								    !
c      u - solution of geometrical optics equations at ray point    !
c      u(1) = x                                                     !
c      u(2) = y                                                     !
c      u(3) = z                                                     !
c      u(4) = n_x                                                   !
c      u(5) = n_y                                                   !
c      u(6) = n_z                                                   !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
c      							   	 
c        output :					    
c                                      
c     Here, deru(i) are  (different from rside_xyz !!!)
c     deru(1)= d(hamiltonian)/d(N_x)    
c     deru(2)= d(hamiltonian)/d(N_y)   
c     deru(3)= d(hamiltonian)/d(N_z)    
c     deru(4)=-d(hamiltonian)/d(x)      
c     deru(5)=-d(hamiltonian)/d(y)       
c     deru(6)=-d(hamiltonian)/d(z)      
c-------------------------------------------------------------------
c     this program uses following functions and subroutines	    !
c     bxyz, gamma1_xyz,
c     dwpw_2,dwcw,s_xyz,	
c     hamilt_xyz    						
c-------------------------------------------------------------------
      subroutine dddrz1_xyz(u,deru)
c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none
      include 'param.i'
      include 'one.i'
c this line only for call control
      include 'eps.i'
c-----input
      double precision u(6), x,y,z, cn_x,cn_y,cn_z, cnx,cny,cnz, cn2
c-----output
      double precision deru(6)
c-----externals
      double precision bxyz, gamma1_xyz, hamilt_xyz,
     &cn, wpw_2, wcw
c-----locals
      double precision vp(nbulka),wp(nbulka)
      double precision wf, r,
     &hw,hfrqnc,
     &hp,hm,h00,
     &dddz,
     &dddx,
     &dddy,
     &frqncpl,df,frqncmn,dddw,p,
     &ds,dc,ds2,dc2, 
     &dwpw2dz,dwpw2dx,dwpw2dy,dwcwdz,dwcwdx,dwcwdy,
     &xi,yi,py2,py3,py4,px,ds4,det,zer,sqrdet,p1,
     &ddetdx,ddetdy,pz,pz2,dpzdx,dpzdy,dddxi,dddyi,dddc2,
     &dxdw,dydw,
     &pc,s1,s2,s3,s4,s6,s7,xib,yob,delib,
     &peym,peyp,a0e,a1e,b0e,b1e,c0e,c1e,dele,ad,bd,cd,
     &px2,yib,peym2,
     &da1edc,da0edc,db1edc,db0edc,
     &da1ed,da0ed,db1ed,db0ed,dc1edc,dc0edc,dadc2,dbdc2,dcdc2,
     &ppe,ppb,pbym,pbyp,pbym2,a0b,a1b,b0b,b1b,c0b,c1b,
     &da1bdc,da0bdc,db1bdc,db0bdc,dc1bdc,dc0bdc,
     &dadz,dadx,dady,dbdz,dbdx,dbdy,dcdz,dcdx,dcdy,
     &dadw,dbdw,dcdw,pyim,pyip,pyim2,
     &ds1dxi,ds2dxi,ds1dyi,ds2dyi,ds3dxi,ds4dxi,ds3dyi,ds4dyi,
     &ds6dxi,ds7dxi,ds6dyi,ds7dyi,
     &da1exi,da0exi,db1exi,db0exi,dc1exi,dc0exi,
     &da1eyi,da0eyi,db1eyi,db0eyi,dc1eyi,dc0eyi,
     &dadxi,dbdxi,dcdxi,dadyi,dbdyi,dcdyi,
     &da1bxi,da0bxi,db1bxi,db0bxi,dc1bxi,dc0bxi,
     &da1byi,da0byi,db1byi,db0byi,dc1byi,dc0byi,
     &ds1dxb,ds2dxb,ds1dyb,ds2dyb,ds3dxb,ds4dxb,ds3dyb,ds4dyb,
     &ds2dxe,ds2dye,
     &ds6dxb,ds7dxb,ds6dyb,ds7dyb,
     &da1bxb,da0bxb,db1bxb,db0bxb,dc1bxb,dc0bxb,
     &da1byb,da0byb,db1byb,db0byb,dc1byb,dc0byb,
     &xe,ye,
     &ds1dxe,ds3dxe,ds4dxe,ds6dxe,ds7dxe,
     &ds1dye,ds3dye,ds4dye,ds6dye,ds7dye,
     &da1bxe,da0bxe,db1bxe,db0bxe,dc1bxe,dc0bxe,
     &da1bye,da0bye,db1bye,db0bye,dc1bye,dc0bye,
     &da1exe,da0exe,db1exe,db0exe,dc1exe,dc0exe,
     &da1eye,da0eye,db1eye,db0eye,dc1eye,dc0eye,
     &dcn,dcn2,dcn4,dadnz,dbdnz,dcdnz,
     &dadnx,dbdnx,dcdnx,dadny,dbdny,dcdny,
     &p4,p5,p6,p7,
     &sum,dddn,
     +dddcnx,dddcny,dddcnz,
     +xp,xm,yp,ym,zp,zm,cn2plus,cn2minus

      ! dwpw2dx,dwpw2dy,dwpw2dz ! derivs of (wp/w)^2
      ! dwcwdx,dwcwdy,dwcwdz ! derivs of wc/w
      
      integer i

      x= u(1) 
      y= u(2) 
      z= u(3) 
      cn_x= u(4) ! Nx
      cn_y= u(5) ! Ny
      cn_z= u(6) ! Nz
      cn2= cn_x*cn_x + cn_y*cn_y + cn_z*cn_z
      r=sqrt(x*x+y*y)
c-------------------------------------------------------
      wf=frqncy
      do 11 i=1,nbulk
         vp(i)=v(i) !  just to save these values
         wp(i)=w(i)
 11   continue
c-------------------for analytical derivatives-----------
c	idif=1
c-------------------for numerical derivatives------------
c	idif=2
c------------------------------------------------------------------      
      if (idif.eq.1) goto 1955 ! 
c----------------------------------------------------------------
c idif=2 numerical hamiltonian derivative calculation
c----------------------------------------------------------------
      hw=der_f !! for frequency
      pi=3.1415926d0
      hfrqnc=hw
      hfrqnc=hfrqnc*frqncy ! =der_f*frqncy
c----------------------------------------------------------------
      bmod=bxyz(x,y,z) ! to find bx,by,bz,bmod -> /one.i/
      cn2= cn_x*cn_x + cn_y*cn_y + cn_z*cn_z
      gam= gamma1_xyz(x,y,z, cn_x,cn_y,cn_z) ! uses bx,by,bz,bmod /one.i/
      h00=  hamilt_xyz(x,y,z, cn2) ! uses gam; calls b(), wpw_2(), wcw()

      cnx= cn_x+der_n
      cn2= cnx*cnx + cn_y*cn_y + cn_z*cn_z
      gam= gamma1_xyz(x,y,z, cnx,cn_y,cn_z) ! uses bx,by,bz,bmod /one.i/
      hp=  hamilt_xyz(x,y,z, cn2) ! uses gam; calls b(), wpw_2(), wcw()
      cnx= cn_x-der_n
      cn2= cnx*cnx + cn_y*cn_y + cn_z*cn_z
      gam= gamma1_xyz(x,y,z, cnx,cn_y,cn_z) 
      hm=  hamilt_xyz(x,y,z, cn2) ! uses gam
      dddcnx=(hp-hm)/(2.d0*der_n)
      deru(1)=dddcnx
      !--- Another version: one-side derivative:
      !dddcnx=(hp-h00)/(der_n)
      !dddcnx=(h00-hm)/(der_n)
      !deru(1)=dddcnx
      !----------------------------------------

      cny= cn_y+der_n
      cn2= cn_x*cn_x + cny*cny + cn_z*cn_z
      gam= gamma1_xyz(x,y,z, cn_x,cny,cn_z) 
      hp=  hamilt_xyz(x,y,z, cn2) ! uses gam
      cny= cn_y-der_n
      cn2= cn_x*cn_x + cny*cny + cn_z*cn_z
      gam= gamma1_xyz(x,y,z, cn_x,cny,cn_z) 
      hm=  hamilt_xyz(x,y,z, cn2) ! uses gam
      dddcny=(hp-hm)/(2.d0*der_n)
      deru(2)=dddcny

      cnz= cn_z+der_n
      cn2= cn_x*cn_x + cn_y*cn_y + cnz*cnz
      gam= gamma1_xyz(x,y,z, cn_x,cn_y,cnz) 
      hp=  hamilt_xyz(x,y,z, cn2) ! uses gam
      cnz= cn_z-der_n
      cn2= cn_x*cn_x + cn_y*cn_y + cnz*cnz
      gam= gamma1_xyz(x,y,z, cn_x,cn_y,cnz) 
      hm=  hamilt_xyz(x,y,z, cn2) ! uses gam
      dddcnz=(hp-hm)/(2.d0*der_n)
      deru(3)=dddcnz

      cn2= cn_x*cn_x + cn_y*cn_y + cn_z*cn_z

      xp=  x+der_r
      bmod=bxyz(xp,y,z) !-> get b and derivs of b
      gam= gamma1_xyz(xp,y,z, cn_x,cn_y,cn_z) ! uses bx,by,bz,bmod /one.i/
      hp=  hamilt_xyz(xp,y,z, cn2) ! uses gam; calls b(), wpw_2(), wcw()
      xm=  x-der_r
      bmod=bxyz(xm,y,z) !-> get b and derivs of b
      gam= gamma1_xyz(xm,y,z, cn_x,cn_y,cn_z)
      hm=  hamilt_xyz(xm,y,z, cn2)
      dddx=(hp-hm)/(2.d0*der_r)
      deru(4)=-dddx
      !--- Another version: one-side derivative:
      !dddx=(hp-h00)/(der_r) ! (H(x+dr)-H00)/dr
      !dddx=(h00-hm)/(der_r) ! (H00-H(x-dr))/dr
      !deru(4)=-dddx
      !----------------------------------------
           
      yp=  y+der_r
      bmod=bxyz(x,yp,z) !-> get b and derivs of b
      gam= gamma1_xyz(x,yp,z, cn_x,cn_y,cn_z)
      hp=  hamilt_xyz(x,yp,z, cn2)
      ym=  y-der_r
      bmod=bxyz(x,ym,z) !-> get b and derivs of b
      gam= gamma1_xyz(x,ym,z, cn_x,cn_y,cn_z)
      hm=  hamilt_xyz(x,ym,z, cn2)
      dddy=(hp-hm)/(2.d0*der_r)
      deru(5)=-dddy

      zp=  z+der_r
      bmod=bxyz(x,y,zp) !-> get b and derivs of b
      gam= gamma1_xyz(x,y,zp, cn_x,cn_y,cn_z)
      hp=  hamilt_xyz(x,y,zp, cn2)
      zm=  z-der_r
      bmod=bxyz(x,y,zm) !-> get b and derivs of b
      gam= gamma1_xyz(x,y,zm, cn_x,cn_y,cn_z)
      hm=  hamilt_xyz(x,y,zm, cn2)
      dddz=(hp-hm)/(2.d0*der_r)
      deru(6)=-dddz

c-----------------------------------------------------------
      bmod=bxyz(x,y,z) ! to find bx,by,bz,bmod -> /one.i/
c*************************************************
!      frqncpl=frqncy+hfrqnc ! =frqncy*(1.d0+der_f)
!      df=frqncy/frqncpl
!      do 12 i=1,nbulk
!         v(i)=vp(i)*df*df ! used through common one.i
!         w(i)=wp(i)*df
! 12   continue
!      cnx=cn_x*df
!      cny=cn_y*df
!      cnz=cn_z*df
!      cn2plus=cn2*df*df
!      gam= gamma1_xyz(x,y,z, cnx,cny,cnz)  
!      hp=  hamilt_xyz(x,y,z, cn2plus) ! uses gam
c*************************************************
!      frqncmn=frqncy-hfrqnc ! =frqncy*(1.d0-der_f)
!      df=frqncy/frqncmn
!      do 15 i=1,nbulk
!         v(i)=vp(i)*df*df ! used through common one.i
!         w(i)=wp(i)*df
! 15   continue
!      cnx=cn_x*df
!      cny=cn_y*df
!      cnz=cn_z*df
!      cn2minus=cn2*df*df
!      gam= gamma1_xyz(x,y,z, cnx,cny,cnz)  
!      hm=  hamilt_xyz(x,y,z, cn2minus) ! uses gam
!c*************************************************
!      !write(*,'(6e12.3)')deru
!      dddw=(hp-hm)/(2.0d0*hw)  
!      if(dddw.eq.0.d0)then
!        write(*,*) 'dddrz1_xyz:  dddw=0'
!        stop
!      else
!        p=-1.d0/dddw
!      endif
      
      do 14 i=1,nbulk
          v(i)=vp(i) ! restore these arrays to original values
          w(i)=wp(i)
 14   continue

ccc YuP Lines with ccc are present in rside1
ccc           if (i_geom_optic.eq.2) then
c------------the right hand side of ray-tracing equations 
c            (dD/dN)/dl and (dD/dr)/dl
c            It gives the ray equations in the form dr/dl and dN/dl
c            l is the lenth along the ray            
ccc              p=1.d0/dsqrt(deru(1)**2+deru(2)**2+deru(3)**2)
ccc              p=ray_direction*p
ccc          endif
          
ccc        do 16 i=1,6
ccc 16	  deru(i)=deru(i)*p
          
c-----------------------------------------------------------
      go to 1953 !-> return/end
c----------------------------------------------------------------
c  end of numerical calculation of hamiltonian derivatives
c------------------------------------------------------------------
 1955	  continue

      x= u(1) 
      y= u(2) 
      z= u(3) 
      cn_x= u(4) ! Nx
      cn_y= u(5) ! Ny
      cn_z= u(6) ! Nz
      cnx=cn_x
      cny=cn_y
      cnz=cn_z
      bmod=bxyz(x,y,z) !-> get b and derivs of b
      gam=gamma1_xyz(x,y,z, cnx,cny,cnz)
      ds=dsin(gam)
      dc=dcos(gam)
      ds2=ds*ds
      dc2=dc*dc

      if (id.eq.3) then
c
c         Appleton-Hartree dispersion relation
c
         if (idif.eq.1) then ! analytical derivs
            stop 'dddrz1_xyz: id=3/idif=1  Not ready yet. Choose idif=2'
            ! For now, only numerical derivs.
         endif
         call dwpw_2(x,y,z,1, dwpw2dx,dwpw2dy,dwpw2dz) ! derivs of (wp/w)^2
         call dwcw(x,y,z,1, dwcwdx,dwcwdy,dwcwdz) ! derivs of wc/w
         ! need to define dxidr,dxidz,etc.
         xi=wpw_2(x,y,z,1)
         yi=wcw(x,y,z,1)

         py2=yi*yi
         py3=py2*yi
         py4=py2*py2
         px=1.d0-xi
         px2=px*px
         ds4=ds2*ds2

         det=py4*ds4+4.d0*py2*px2*dc2
         zer=0.d0
         if (det.lt.zer) then
            write(*,*)'det in rside less then 0 det=',det
            stop
         end if
         sqrdet=dsqrt(det)
         p1=0.5d0/sqrdet
         ddetdx=-8.d0*py2*px*dc2
         ddetdy=4.d0*py3*ds4+8.d0*yi*px2*dc2
         pz=2.d0*px-py2*ds2+ioxm*sqrdet
         pz2=1.d0/(pz*pz)
         dpzdx=-2.d0+ioxm*p1*ddetdx
         dpzdy=-2.d0*yi*ds2+ioxm*p1*ddetdy
         dddxi=-2.d0*((2.d0*xi-1.d0)*pz+xi*px*dpzdx)*pz2
         dddyi=-2.d0*xi*px*dpzdy*pz2
         dddc2=-2.d0*xi*px*pz2*(py2+ioxm*(-2.d0*py4*(1.d0-dc2)+
     1                    4.d0*py2*px2)*p1)
c-----------------------------------------------------------------
         dddz=dddxi*dwpw2dz+dddyi*dwcwdz+dddc2*dc2dz
         dddx=dddxi*dwpw2dx+dddyi*dwcwdx+dddc2*dc2dx
         dddy=dddxi*dwpw2dy+dddyi*dwcwdy+dddc2*dc2dy
         dddcnz=2.d0*cnz+dddc2*dc2dnz
         dddcnx=2.d0*cnx+dddc2*dc2dnx
         dddcny=2.d0*cny+dddc2*dc2dny
c********************************************
ccc YuP Lines with ccc are present in rside1()
ccc	     dxdw=-2.d0*xi
ccc	     dydw=-yi
ccc         dddw=dddxi*dxdw+dddyi*dydw
ccc     1		  -2.d0*(cnz*cnz+cnx*cnx+cny*cny)
ccc             dddw=dddw/wf 
c--------------------------------------------------------------------
         goto 50
      end if ! id=3
c   end if Appleton - Hartree
c-----------------------------------------------------------------
           if ((id.eq.1).or.(id.eq.2)) then
c  ---------------------------------------------------------------
c  cold plasma dispersion relation with electrons and ions
c  wpw_2(i=1),wcw(i=1) electrons component
c  wpw_2(i may be=2,nbulk),wcw(i may be=2,nbulk) ions components
c  ib number of component for which delib=1-yib may be equal
c  zero inside the plasma,ib  may be=from 1 to nbulk
c  dispersion relation is multiplied by delib
c
c         if (idif.eq.1) then ! analytical derivs
c            stop 'dddrz1_xyz: id=1,2/idif=1  Not ready. Choose idif=2' 
c            ! For now, only numerical derivs.
c         endif

	   pc=1.d0+dc2

         call s_xyz(x,y,z,s1,s2,s3,s4,s6,s7)

         xib=wpw_2(x,y,z,ib)
         yib=wcw(x,y,z,ib)
         delib=1.d0-yib

c  ib =1 (cyclotron resonance conditions dele=0 may be
c         realised in plasma, dispersion relation is
c         multiplied by dele )
c         ad,bd,cd calculations in dispersion relation
c       d=a*n**4+b*n**2+c
c------------------------------------------------------------------
c  if 2 begin
      if (ib.eq.1) then
       xe=xib
       ye=yib
         peym=xe/(1.d0-ye)
	 peyp=xe/(1.d0+ye)
	   a0e=-peyp*ds2
	   a1e=s7*ds2+s4*dc2
	   b0e=s4*peyp*pc+xe*(s6-peyp)*ds2
	   b1e=-s4*s7*pc-s3*(s6-peyp)*ds2
	   c0e=-xe*s4*(s6-peyp)
	   c1e=s4*s3*(s6-peyp)
	   dele=1.d0-ye
	     ad=dele*a1e+a0e
	     bd=dele*b1e+b0e
	     cd=dele*c1e+c0e
c----------------------------------------------------------------------
	     da1edc=-s7+s4
	     da0edc=peyp
	     db1edc=-s7*s4+s3*(s6-peyp)
	     db0edc=s4*peyp-xe*(s6-peyp)
	     dc1edc=0.d0
	     dc0edc=0.d0

	     dadc2=dele*da1edc+da0edc
	     dbdc2=dele*db1edc+db0edc
	     dcdc2=dele*dc1edc+dc0edc
c-------------------------------------------------------------------
      end if
c  if 2 end
c-----------------------------------------------------------------
c  ib .gt.1 (cyclotron resonance conditions delib=0 may be
c            realised in plasma, dispersion relation is
c            multiplied by delib )
c            ad,bd,cd calculations
c	     d=a*n**4+b*n**2+c
c
c  if 3 begin
      if (ib.gt.1) then
       xe=wpw_2(x,y,z,1)
       ye=wcw(x,y,z,1)
         ppe=1.d0/(1.d0-ye)
         peym=xe/(1.d0-ye)
	 peyp=xe/(1.d0+ye)
	 peym2=peyp*ppe

         ppb=1.d0/(1.d0-yib)
         pbym=xib/(1.d0-yib)
	 pbyp=xib/(1.d0+yib)
	 pbym2=pbyp*ppb

	   a0b=-pbyp*ds2
	   a1b=(s1-peym2)*ds2+s4*dc2
	   b0b=s4*pbyp*pc+xib*(s3-peym)*ds2
	   b1b=-s4*(s1-peym2)*pc-(s2-peyp)*(s3-peym)*ds2
	   c0b=-xib*s4*(s3-peym)
	   c1b=s4*(s2-peyp)*(s3-peym)
	   delib=1.d0-yib
	     ad=delib*a1b+a0b
	     bd=delib*b1b+b0b
	     cd=delib*c1b+c0b
c--------------------------------------------------------------------
	     da1bdc=-(s1-peym2)+s4
	     da0bdc=pbyp
	     db1bdc=-s4*(s1-peym2)+(s2-peyp)*(s3-peym)
	     db0bdc=s4*pbyp-xib*(s3-peym)
	     dc1bdc=0.d0
	     dc0bdc=0.d0

	     dadc2=delib*da1bdc+da0bdc
	     dbdc2=delib*db1bdc+db0bdc
 	     dcdc2=delib*dc1bdc+dc0bdc
c-------------------------------------------- -----------0-----------
      end if
c  if 3 end
c-----------------------------------------------------------------

c  derivatives calculations
c   dadz,dadx,dady
c   dbdz,dbdx,dbdy
c   dcdz,dcdx,dcdy
c   dadw,dbdw,dcdw
c----------------------------------------------------------------
       dadz=0.d0
       dadx=0.d0
       dady=0.d0
       dbdz=0.d0
       dbdx=0.d0
       dbdy=0.d0
       dcdz=0.d0
       dcdx=0.d0
       dcdy=0.d0
       dadw=0.0d0
       dbdw=0.0d0
       dcdw=0.0d0

c  if nbulk.gt.1  (ions components are in plasma )
c
c  if 4 begin
      if (nbulk.gt.1) then
       do 10 i=2,nbulk
c
c  derivatives calculations
c   dadx,dady
c   dbdx,dbdy
c   dcdx,dcdj

         call dwpw_2(x,y,z,i, dwpw2dx,dwpw2dy,dwpw2dz) ! derivs of (wp/w)^2
         call dwcw(x,y,z,i, dwcwdx,dwcwdy,dwcwdz) ! derivs of wc/w

         if ( i.eq.ib) goto 20
c
c  i.ne.ib
c
         xi=wpw_2(x,y,z,i)
         yi=wcw(x,y,z,i)

	   pyim=1.d0/(1.d0-yi)
	   pyip=1.d0/(1.d0+yi)
	   pyim2=pyim*pyim

c if 5 begin
	     if (ib.gt.1) then
	       ds1dxi=pyim*pyip
	       ds2dxi=pyim
	       ds1dyi=2.d0*xi*yi*ds1dxi*ds1dxi
	       ds2dyi=xi*pyim2
	     end if
c if 5 end
             ds3dxi=pyip
	     ds4dxi=1.d0
             ds3dyi=-xi*pyip*pyip
             ds4dyi=0.d0
c if 6 begin
	     if (ib.eq.1) then
	       ds6dxi=pyim
	       ds7dxi=pyim*pyip
	       ds6dyi=xi*pyim2
	       ds7dyi=2.d0*xi*yi*ds7dxi*ds7dxi
	     end if
c if 6 end

c if 7 begin
	     if (ib.eq.1) then

               da1exi=-ds7dxi*ds2-ds4dxi*dc2
	       da0exi=0.d0
	       db1exi=(ds4dxi*s7+ds7dxi*s4)*pc+
     1		      (ds3dxi*(s6-xe/(1.d0+ye))+ds6dxi*s3)*ds2
               db0exi=-ds4dxi*xe/(1.d0+ye)*pc-
     1                 xe*ds6dxi*ds2
               dc1exi=-ds4dxi*s3*(s6-xe/(1.d0+ye))-
     1                 s4*ds3dxi*(s6-xe/(1.d0+ye))-s4*s3*ds6dxi
               dc0exi=xe*ds4dxi*(s6-xe/(1.d0+ye))+xe*s4*ds6dxi

	       da1eyi=-ds7dyi*ds2-ds4dyi*dc2
	       da0eyi=0.d0
	       db1eyi=(ds4dyi*s7+ds7dyi*s4)*pc+
     1                (ds3dyi*(s6-xe/(1.d0+ye))+ds6dyi*s3)*ds2
               db0eyi=-ds4dyi*xe/(1.d0+ye)*pc-xe*ds6dyi*ds2
	       dc1eyi=-ds4dyi*s3*(s6-xe/(1.d0+ye))-
     1    	       ds3dyi*s4*(s6-xe/(1.d0+ye))-ds6dyi*s4*s3
               dc0eyi=xe*ds4dyi*(s6-xe/(1.d0+ye))+ds6dyi*xe*s4

	         dadxi=dele*da1exi+da0exi
	         dbdxi=dele*db1exi+db0exi
	         dcdxi=dele*dc1exi+dc0exi
                 dadyi=dele*da1eyi+da0eyi
                 dbdyi=dele*db1eyi+db0eyi
                 dcdyi=dele*dc1eyi+dc0eyi

	       goto 30
	     end if
c if 7 end

c if 8 begin
	     if (ib.gt.1) then

               da1bxi=-ds1dxi*ds2-ds4dxi*dc2
	       da0bxi=0.d0
	       db1bxi=(ds4dxi*(s1-xe/(1.d0-ye*ye))+ds1dxi*s4)*pc+
     1		  (ds2dxi*(s3-xe/(1.d0-ye))+ds3dxi*(s2-xe/(1.d0+ye)))*
     2                 ds2
               db0bxi=-ds4dxi*xib/(1.d0+yib)*pc-
     1                 xib*ds3dxi*ds2
               dc1bxi=-ds4dxi*(s2-xe/(1.d0+ye))*(s3-xe/(1.d0-ye))-
     1                 ds2dxi*s4*(s3-xe/(1.d0-ye))-
     2                 ds3dxi*s4*(s2-xe/(1.d0+ye))
               dc0bxi=xib*ds4dxi*(s3-xe/(1.d0-ye))+xib*s4*ds3dxi

	       da1byi=-ds1dyi*ds2-ds4dyi*dc2
	       da0byi=0.d0
	       db1byi=(ds4dyi*(s1-xe/(1.d0-ye*ye))+s4*ds1dyi)*pc+
     1          (ds2dyi*(s3-xe/(1.d0-ye))+ds3dyi*(s2-xe/(1.d0+ye)))*
     2                 ds2
               db0byi=-ds4dyi*xib/(1.d0+yib)*pc-xib*ds3dyi*ds2
	       dc1byi=s4*(-ds2dyi*(s3-xe/(1.d0-ye))-
     1    	       ds3dyi*(s2-xe/(1.d0+ye)))
               dc0byi=xib*s4*ds3dyi

                 dadxi=delib*da1bxi+da0bxi
                 dbdxi=delib*db1bxi+db0bxi
                 dcdxi=delib*dc1bxi+dc0bxi
                 dadyi=delib*da1byi+da0byi
                 dbdyi=delib*db1byi+db0byi
                 dcdyi=delib*dc1byi+dc0byi

	       goto 30
	     end if
c if 8 end
 20     continue
c
c         i=ib, i.ne.1
c

c if 15 begin
	     if (ib.gt.1) then
	       ds1dxb=0.d0
	       ds2dxb=0.d0
	       ds1dyb=0.d0
	       ds2dyb=0.d0
	     end if
c if 15 end
             ds3dxb=1.d0/(1.d0+yib)
	     ds4dxb=1.d0
             ds3dyb=-xib/(1.d0+yib)**2
             ds4dyb=0.d0
c if 16 begin
	     if (ib.eq.1) then
	       ds6dxb=0.d0
	       ds7dxb=0.d0
	       ds6dyb=0.d0
	       ds7dyb=0.d0
	     end if
c if 16 end

c if 9 begin
	     if (ib.gt.1) then

               da1bxb=-ds1dxb*ds2-ds4dxb*dc2
               da0bxb=-1.d0/(1.d0+yib)*ds2
               db1bxb=(ds4dxb*(s1-xe/(1.d0-ye*ye))+ds1dxb*s4)*pc+
     1            (ds2dxb*(s3-xe/(1.d0-ye))+ds3dxb*(s2-xe/(1.d0+ye)))*
     2                 ds2
               db0bxb=(-ds4dxb*xib/(1.d0+yib)+s4/(1.d0+yib))*pc+
     1                ((s3-xe/(1.d0-ye))-xib*ds3dxb)*ds2
               dc1bxb=-ds4dxb*(s2-xe/(1.d0+ye))*(s3-xe/(1.d0-ye))-
     1                 ds2dxb*s4*(s3-xe/(1.d0-ye))-
     2                 ds3dxb*s4*(s2-xe/(1.d0+ye))
            dc0bxb=-s4*(s3-xe/(1.d0-ye))+xib*ds4dxb*(s3-xe/(1.d0-ye))+
     +                xib*s4*ds3dxb

               da1byb=0.d0
               da0byb=xib/(1.d0+yib)**2*ds2
               db1byb=(ds4dyb*(s1-xe/(1.d0-ye*ye))+s4*ds1dyb)*pc+
     1          (ds2dyb*(s3-xe/(1.d0-ye))+ds3dyb*(s2-xe/(1.d0+ye)))*
     2                 ds2
               db0byb=(-ds4dyb*xib/(1.d0+yib)-s4*xib/(1.d0+yib)**2)*pc-
     -                xib*ds3dyb*ds2
               dc1byb=s4*(-ds2dyb*(s3-xe/(1.d0-ye))-
     1                 ds3dyb*(s2-xe/(1.d0+ye)))
               dc0byb=xib*s4*ds3dyb

                 dadxi=delib*da1bxb+da0bxb
                 dbdxi=delib*db1bxb+db0bxb
                 dcdxi=delib*dc1bxb+dc0bxb
                 dadyi=delib*da1byb+da0byb-a1b
                 dbdyi=delib*db1byb+db0byb-b1b
                 dcdyi=delib*dc1byb+dc0byb-c1b

	     end if
c if 9 end

 30       continue

          dadz=dadz+dadxi*dwpw2dz+dadyi*dwcwdz
          dadx=dadx+dadxi*dwpw2dx+dadyi*dwcwdx
          dady=dady+dadxi*dwpw2dy+dadyi*dwcwdy

          dbdz=dbdz+dbdxi*dwpw2dz+dbdyi*dwcwdz
          dbdx=dbdx+dbdxi*dwpw2dx+dbdyi*dwcwdx
          dbdy=dbdy+dbdxi*dwpw2dy+dbdyi*dwcwdy

          dcdz=dcdz+dcdxi*dwpw2dz+dcdyi*dwcwdz
          dcdx=dcdx+dcdxi*dwpw2dx+dcdyi*dwcwdx
          dcdy=dcdy+dcdxi*dwpw2dy+dcdyi*dwcwdy

         xi=wpw_2(x,y,z,i)
         yi=wcw(x,y,z,i)
c*************************************c
	  dadw=dadw-2*dadxi*xi-dadyi*yi
	  dbdw=dbdw-2*dbdxi*xi-dbdyi*yi
	  dcdw=dcdw-2*dcdxi*xi-dcdyi*yi
c new--------------------------------------

 10       continue
         end if
c if 4 nbulk > 1 end


c
c     i = 1
c
             ds1dxe=0.d0
             ds2dxe=0.d0
             ds3dxe=0.d0
             ds4dxe=1.d0
             ds6dxe=0.d0
             ds7dxe=0.d0

             ds1dye=0.d0
             ds2dye=0.d0
             ds3dye=0.d0
             ds4dye=0.d0
             ds6dye=0.d0
             ds7dye=0.d0

c if 10 begin
c
c   ib > 1, i = 1
c
	     if (ib.gt.1) then
               da1bxe=-1.d0/(1.d0-ye*ye)*ds2-ds4dxe*dc2
               da0bxe=0.d0
               db1bxe=(ds4dxe*(s1-xe/(1.d0-ye*ye))+s4/(1.d0-ye*ye))*pc+
     1                (1.d0/(1.d0+ye)*(s3-xe/(1.d0-ye))+
     +                 1.d0/(1.d0-ye)*(s2-xe/(1.d0+ye)))*ds2
               db0bxe=-ds4dxe*xib/(1.d0+yib)*pc-xib/(1.d0-ye)*ds2
               dc1bxe=s4*(-1.d0/(1.d0+ye)*(s3-xe/(1.d0-ye))-
     -                     1.d0/(1.d0-ye)*(s2-xe/(1.d0+ye)))-
     -                 ds4dxe*((s2-xe/(1.d0+ye))*(s3-xe/(1.d0-ye)))
               dc0bxe=xib*s4/(1.d0-ye)+xib*ds4dxe*(s3-xe/(1.d0-ye))

               da1bye=-2.d0*xe*ye/(1.d0-ye*ye)**2*ds2-ds4dye*dc2
               da0bye=0.d0
               db1bye=s4*2.*xe*ye/(1.d0-ye*ye)**2*pc+
     +                (-xe/(1.d0+ye)**2*(s3-xe/(1.d0-ye))+
     +                  xe/(1.d0-ye)**2*(s2-xe/(1.d0+ye)))*ds2
               db0bye=-xib*xe/(1.d0-ye)**2*ds2
               dc1bye=s4*(xe/(1.d0+ye)**2*(s3-xe/(1.d0-ye))-
     -                    xe/(1.d0-ye)**2*(s2-xe/(1.d0+ye)))
               dc0bye=xib*s4*xe/(1.d0-ye)**2

                 dadxi=delib*da1bxe+da0bxe
                 dbdxi=delib*db1bxe+db0bxe
                 dcdxi=delib*dc1bxe+dc0bxe
                 dadyi=delib*da1bye+da0bye
                 dbdyi=delib*db1bye+db0bye
                 dcdyi=delib*dc1bye+dc0bye

               goto 40
	     end if
c if 10 end

c if 11 begin
c
c   ib = 1, i = 1
c
             if (ib.eq.1) then

               da1exe=-ds4dxe*dc2
               da0exe=-1.d0/(1.d0+ye)*ds2
               db1exe=ds4dxe*s7*pc+s3/(1.d0+ye)*ds2
               db0exe=(-ds4dxe*xe/(1.d0+ye)+s4/(1.d0+ye))*pc+
     +                 (s6-2.d0*xe/(1.d0+ye))*ds2
               dc1exe=-ds4dxe*s3*(s6-xe/(1.d0+ye))-
     -                s4*s3/(1.d0+ye)
               dc0exe=-s4*(s6-2.d0*xe/(1.d0+ye))+xe*(s6-xe/(1.d0+ye))

               da1eye=0.d0
               da0eye=xe/(1.d0+ye)**2*ds2
               db1eye=-s3*xe/(1.d0+ye)**2*ds2
               db0eye=-s4*xe/(1.d0+ye)**2*pc+xe*xe/(1.d0+ye)**2*ds2
               dc1eye=-ds4dye*s3*(s6-xe/(1.d0+ye))+s4*s3*xe/(1.d0+ye)**2
               dc0eye=-xe*xe*s4/(1.d0+ye)**2

                 dadxi=dele*da1exe+da0exe
                 dbdxi=dele*db1exe+db0exe
                 dcdxi=dele*dc1exe+dc0exe
                 dadyi=dele*da1eye+da0eye-a1e
                 dbdyi=dele*db1eye+db0eye-b1e
                 dcdyi=dele*dc1eye+dc0eye-c1e


	     end if
c if 11 end

 40     continue
         call dwpw_2(x,y,z,1, dwpw2dx,dwpw2dy,dwpw2dz) ! derivs of (wp/w)^2
         call dwcw(x,y,z,1, dwcwdx,dwcwdy,dwcwdz) ! derivs of wc/w
c------------------------------------------------------------------
        dadz=dadz+dadxi*dwpw2dz+dadyi*dwcwdz+dadc2*dc2dz
        dadx=dadx+dadxi*dwpw2dx+dadyi*dwcwdx+dadc2*dc2dx
        dady=dady+dadxi*dwpw2dy+dadyi*dwcwdy+dadc2*dc2dy

        dbdz=dbdz+dbdxi*dwpw2dz+dbdyi*dwcwdz+dbdc2*dc2dz
        dbdx=dbdx+dbdxi*dwpw2dx+dbdyi*dwcwdx+dbdc2*dc2dx
        dbdy=dbdy+dbdxi*dwpw2dy+dbdyi*dwcwdy+dbdc2*dc2dy

        dcdz=dcdz+dcdxi*dwpw2dz+dcdyi*dwcwdz+dcdc2*dc2dz
        dcdx=dcdx+dcdxi*dwpw2dx+dcdyi*dwcwdx+dcdc2*dc2dx
        dcdy=dcdy+dcdxi*dwpw2dy+dcdyi*dwcwdy+dcdc2*dc2dy
c---------------------------------------------------------------------
c new-------------------------------------
       xi=wpw_2(x,y,z,1)
       yi=wcw(x,y,z,1)
	  dadw=dadw-2*dadxi*xi-dadyi*yi
	  dbdw=dbdw-2*dbdxi*xi-dbdyi*yi
	  dcdw=dcdw-2*dcdxi*xi-dcdyi*yi

c new--------------------------------------

        dcn=dsqrt(cnz*cnz+cnx*cnx+cny*cny)
        dcn2=dcn*dcn
        dcn4=dcn2*dcn2

	dadnz=dadc2*dc2dnz
	dbdnz=dbdc2*dc2dnz
	dcdnz=dcdc2*dc2dnz

	dadnx=dadc2*dc2dnx
	dbdnx=dbdc2*dc2dnx
	dcdnx=dcdc2*dc2dnx

	dadny=dadc2*dc2dny
	dbdny=dbdc2*dc2dny
	dcdny=dcdc2*dc2dny
c------------------------------------------------------------------

        if (id.eq.1) then

c
c  dispersion relation a*n**4+b*n**2+c=0
c

          dddz=dcn4*dadz+dcn2*dbdz+dcdz
          dddx=dcn4*dadx+dcn2*dbdx+dcdx
          dddy=dcn4*dady+dcn2*dbdy+dcdy

          dddn=4.d0*ad*dcn2+2.d0*bd
          dddcnz=dddn*cnz+dcn4*dadnz+dcn2*dbdnz+dcdnz
          dddcnx=dddn*cnx+dcn4*dadnx+dcn2*dbdnx+dcdnx
          dddcny=dddn*cny+dcn4*dadny+dcn2*dbdny+dcdny
ccc YuP Lines with ccc are present in rside1
ccc-new-----------------------------------------
ccc          dddw=dcn4*dadw+dcn2*dbdw+dcdw
ccc     1         -4.d0*dcn4*ad-2.d0*dcn2*bd
ccc          dddw=dddw/wf
ccc	  write(*,*)'dddw=',dddw,'dddcnr=',dddcnr,'cnr=',cnr
ccc-new-----------------------------------------

        end if

        if (id.eq.2) then

c
c  dispersion relation n**2-(-b+iom*sqrt(b*2-4*a*c))/(2*a)=0
c

         det=dsqrt(bd*bd-4.d0*ad*cd)
         p4=0.5d0/(ad*ad)
         p5=1.d0/det
         p6=-bd+ioxm*det
c------------------------------------------------------------------
         p7=ioxm*(bd*dbdz-2.d0*dadz*cd-2.d0*ad*dcdz)*p5
         dddz=-((-dbdz+p7)*ad-dadz*p6)*p4
         p7=ioxm*(bd*dbdx-2.d0*dadx*cd-2.d0*ad*dcdx)*p5
         dddx=-((-dbdx+p7)*ad-dadx*p6)*p4
         p7=ioxm*(bd*dbdy-2.d0*dady*cd-2.d0*ad*dcdy)*p5
         dddy=-((-dbdy+p7)*ad-dady*p6)*p4

         p7=ioxm*(bd*dbdnz-2.d0*dadnz*cd-2.d0*ad*dcdnz)*p5
          dddcnz=2.d0*cnz-
     1      	 ((-dbdnz+p7)*ad-dadnz*p6)*p4
         p7=ioxm*(bd*dbdnx-2.d0*dadnx*cd-2.d0*ad*dcdnx)*p5
          dddcnx=2.d0*cnx-
     1      	 ((-dbdnx+p7)*ad-dadnx*p6)*p4
         p7=ioxm*(bd*dbdny-2.d0*dadny*cd-2.d0*ad*dcdny)*p5
          dddcny=2.d0*cny-
     1      	 ((-dbdny+p7)*ad-dadny*p6)*p4
c--------------------------------------------------------------------
ccc YuP Lines with ccc are present in rside1()
ccc new--------------------------------------------------------------
ccc         p7=ioxm*(bd*dbdw-2.d0*dadw*cd-2.d0*ad*dcdw)*p5
ccc         dddw=-((-dbdw+p7)*ad-dadw*p6)*p4
ccc     1	      -2.d0*dcn2
ccc         dddw=dddw/wf
ccc--------------------------------------------------------------------
c--------------------------------------------------------------------

        end if
c--------------------------------------------------------------------
        goto 50
           end if
c   end of cold plasma  dispersion
c-----------------------------------------------------------------
        if (id.eq.6)then
c-----------------------------------------------------------------
c       hot non-relativistic plasma dispersion from Forest code
c-----------------------------------------------------------------
cyup          call hotdervs(nbulk,u,wf,dddcnz,dddcnr,dddcm,
cyup     .                  dddz,dddr,dddphi,dddw)
           call hotdervs_xyz(u,wf,dddcnx,dddcny,dddcnz,
     .                       dddx,dddy,dddz,dddw)
           goto 50
        end if
c       end of hot non-relativistic plasma dispersion from Forest code

c-----------------------------------------------------------------
c-----------------------------------------------------------------
 50     continue
           deru(1)=dddcnx
           deru(2)=dddcny
           deru(3)=dddcnz
           deru(4)=-dddx
           deru(5)=-dddy
           deru(6)=-dddz

           
ccc           dddw=-dddw*wf
            
ccc           if (i_geom_optic.eq.1) then
ccc             deru(1)=deru(1)/dddw
ccc             deru(2)=deru(2)/dddw
ccc             deru(3)=deru(3)/dddw
ccc             deru(4)=deru(4)/dddw
ccc             deru(5)=deru(5)/dddw
ccc             deru(6)=deru(6)/dddw
ccc           else
ccc              if(i_geom_optic.eq.2) then
ccc                p=1.d0/dsqrt(deru(1)**2+deru(2)**2+(r*deru(3))**2)
ccc                p=ray_direction*p
ccc                deru(1)=p*deru(1)
ccc                deru(2)=p*deru(2)
ccc                deru(3)=p*deru(3)
ccc                deru(4)=p*deru(4)
ccc                deru(5)=p*deru(5)
ccc                deru(6)=p*deru(6)
ccc               endif
ccc           endif
ccc YuP: lines with ccc are present in rside1()
c      write(*,*)'dddrz1_xyz: id,ib,ioxm=',id,ib,ioxm

 1953   continue
 
c           write(*,'(a,3e14.5)')' dddrz1_xyz: R, rho, dddw=', 
c     +                            R,rho,(-dddw*wf)
c           write(*,'(a,3e12.3)')' dddrz1_xyz: dddcnx,dddcny,dddcnz',
c     +                            dddcnx,dddcny,dddcnz
c           write(*,'(a,3e12.3)')' dddrz1_xyz      : dddx,dddy,dddz',
c     +                            dddx,dddy,dddz
           !pause !!!

        !write(*,*)' END OF dddrz1_xyz->hamilt_xyz:'
        !hp=  hamilt_xyz(x,y,z, cn2) ! test
        !write(*,'(a,8e12.3)')'r,hp,deru(1:6)=',r,hp,deru(1:6)
        !pause
      return
      end


