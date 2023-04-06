c
c        *********************dense_no_RZ_spline***************
c        *                        -                           *
c        * this function calculates the density profile       *
c        * as a function  from psi(z,r,phi) 		      *
c        *  rho=psi-psimag 				      *
c        *  a=psilim-psimag inside the function b()           *
c        *  these variables are in common 'one'
c        * This function has not chamber wall effects.
c        ******************************************************
c
c------------------------------------------------------------------
c								   !
c        input parameters					   !
c      								   !
c      z, r, phi - point which the components of magnetic field	   !
c                  are calculated in.				   !
c       rho is from common/one/                                    !
c------------------------------------------------------------------
c      uses
c         constants dense0,denseb,rn1de,rn2de,idens are in common 'one'
c         they were read in subroutine dinit from the file genray.in
c         rho is  in common 'one'
c      functions
c         _densrho_(rho)    finds dense from spline
c------------------------------------------------------------------
      real*8 FUNCTION dense_no_RZ_spline(z,r,phi,i)
c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none
      include 'param.i'
      include 'one.i'
c-----input
      real*8 z,r,phi!space coordinates
c     rho is in common  /one/
      integer i !plasma species number
c-----local
      real*8  theta_pol, !0=< theta_pol <2pi
     & rho_small,dens_rho_theta,
     & d_dens_rho_theta_d_rho,d_dens_rho_theta_d_theta
c-----externals
      real*8 thetapol,   ! -pi <thetapol =<pi
     & vardens,densrho
c      write(*,*)'in dense z,r,phi,i,nbulk',z,r,phi,i,nbulk

      if(ixyz.eq.1) then
        stop 'dense_no_RZ_spline() should not be used for ixyz=1'
      endif

      if(i.gt.nbulk)then
        write(*,*)'in dense i.gt.nbulk: i,nbulk',i,nbulk
	stop
      endif
c      write(*,*)'in dense rho',rho
      if (rho.gt.1.d0-1.d-10) then
cSAP090206
c         dense=densrho(rho,i)

         theta_pol=thetapol(z,r) ! -pi <thetapol =<pi
         if (theta_pol.lt.0d0) then
            theta_pol=theta_pol+2*pi !pi< theta_pol<2pi
         endif

c--------calculate density outside LCFS
c         write(*,*)'dense.f in function dense before'
c         write(*,*)'dens_rho_theta_LCFS i',i
         rho_small=rho

c         write(*,*)'dense.f z,r,phi,i,rho_small',z,r,phi,i,rho_small

         call dens_rho_theta_LCFS(rho_small,theta_pol,i,
     &        dens_rho_theta,d_dens_rho_theta_d_rho,
     &        d_dens_rho_theta_d_theta)
          dense_no_RZ_spline = dens_rho_theta 
      else
          dense_no_RZ_spline=densrho(rho,i)
      endif
      return
      END


c        *********************tempe ***************************
c        *                        -                           *
c        * this function calculates the electron temperature  *
c        * profile as a function  from rho       	      *
c        * rho=rhopsi(psi) 				      *
c        * it calculates in function b()		      *
c        * rho is  in common 'one'			      *
c        ******************************************************
c
c------------------------------------------------------------------
c								   !
c        input parameters					   !
c      								   !
c      z, r, phi - point which the components of magnetic field	   !
c                  are calculated in.				   !
c       rho is from common/one/
c------------------------------------------------------------------
	double precision
     1function tempe(z,r,phi,i)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
c------------------------------------------------------------------
c     te0,teb,rn1te,rn2te  are in common 'one'
c     they were read in subroutine dinit from the file genray.in
c           rho is  in common 'one'
c------------------------------------------------------------------
cc     analytic form of the i component temperature profile
c     if(idens.eq.0) then
c	  x1=(rho)**rn1te(i)
c          tempe=(te0(i)-teb(i))*(1.d0-x1)**rn2te(i)+teb(i)
cc                         spline form of the electron temperature profile
c      else
      if(ixyz.eq.1) then
        stop 'tempe(z,r,phi,i) should not be used for ixyz=1'
      endif
      
      tempe=temperho(rho,i)
      return
      end



c        *********************zeff  ***************************
c        *                        -                           *
c        * this function calculates the Z_eff profile         *
c        * as a function  from psi(z,r,phi) 		      *
c        *  rho=psi-psimag 				      *
c        *  a=psilim-psimag inside the function b()           *
c        *  these variables are in common 'one'
c        ******************************************************
c--------------------------------- ---------------------------------
c								   !
c        input parameters					   !
c      								   !
c      z, r, phi - point which the components of magnetic field	   !
c                  are calculated in.				   !
c       rho is from common/one/
c------------------------------------------------------------------
      double precision FUNCTION zeffi(z,r,phi)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'

c------------------------------------------------------------------
c      zeff0,zeffb,rn1zeff,rn2zeff  are in common 'one'
c      they were read in subroutine dinit_mr from the file genray.in
c------------------------------------------------------------------
	  zeffi=zeffrho(rho)
      return
      END


      subroutine sigma_edge_n_theta_pol(theta_pol_radian,
     &sigma_edge_n,d_sigma_edge_n_d_theta_pol,
     &d2_sigma_edge_n_d2_theta_pol)
c-----calculates sigma_edge_n and its derivative by poloidal angle
c     theta_pol [radians]
c     
c     It will use the analytical formula or spline approximation
c     depends on the variable i_edge_dens_anal given in edge_prof_nml.i
c     i_edge_dens_anal, ! to choose analytic formula or spline:
c                       ! =0 sigmedgn=const which was set in 
c                       !    the input namelist/write/
c                       ! =1 analytic formula for sigmedgn(theta_pol)
c                       ! =2 table data for  sigmedgn(theta_pol)
      implicit none
      include 'param.i'
      include 'edge_prof.i'  
      include 'one.i'
c-----input
      real*8 
     &theta_pol_radian       !poloidal angle [radians]
c-----output
      real*8
     &sigma_edge_n,
     &d_sigma_edge_n_d_theta_pol,  !first derivate
     &d2_sigma_edge_n_d2_theta_pol !second derivative
c-----locals 
     
      real*8, dimension(3*n_pol_edge_dens+1) :: workk
      real*8 tabl(3),exp1,exp2,del1,del2,
     & d_arg_exp1_d_theta, d_arg_exp2_d_theta,
     & d2_arg_exp1_d2_theta,d2_arg_exp2_d2_theta

      integer i,iop(2),itabl(3)

      logical first 
      data  first /.true./
      save first
  
      pi=4.d0*datan(1.d0)
     
c      write(*,*)'sigma_edge_n_theta_pol theta_pol_radian',
c     &           theta_pol_radian

c      write(*,*)'i_edge_dens_anal',i_edge_dens_anal
c      write(*,*)'first',first

      if ( n_pol_edge_dens.eq.1) then 
         i_edge_dens_anal=1
          sigma_edgen_0= sigmedgn_ar(1)
      endif

cSAP090222
      if (i_edge_dens_anal.eq.0) then
c        sigma_edge_n=sigmedgn=const set in the input namelist /edge_prof_nml/
         sigma_edge_n=sigmedgn
         d_sigma_edge_n_d_theta_pol=0.d0
         d2_sigma_edge_n_d2_theta_pol=0.d0
         goto 10
      endif

      if (i_edge_dens_anal.eq.1) then
c--------analytical formula---------------
c        
c         sigma_edge_n=sigma_edgen_0
c         d_sigma_edge_n_d_theta_pol=0.d0
c         d2_sigma_edge_n_d2_theta_pol=0.d0

         if (first) then
c-----------transform from degree to radian

           sigma_theta_pol_edge_1_radian=
     &           sigma_theta_pol_edge_1_degree*pi/180.d0    
 
c           write(*,*)'dense.f sigma_theta_pol_edge_1_radian',
c     &                        sigma_theta_pol_edge_1_radian

           sigma_theta_pol_edge_2_radian=
     &           sigma_theta_pol_edge_2_degree*pi/180.d0

           theta_pol_edge_1_radian=theta_pol_edge_1_degree*pi/180.d0
           theta_pol_edge_2_radian=theta_pol_edge_2_degree*pi/180.d0

c          write(*,*)'dense.f sigma_theta_pol_edge_1_radian',
c     &                        sigma_theta_pol_edge_1_radian

           first=.false.  
         endif

c         write(*,*)'dense.f sigma_edge_n_theta_pol'
c         write(*,*)'theta_pol_radian',theta_pol_radian
c         write(*,*)'theta_pol_edge_1_radian',theta_pol_edge_1_radian
c         write(*,*)'theta_pol_edge_2_radian',theta_pol_edge_2_radian
c         write(*,*)'sigma_theta_pol_edge_1_radian',
c     &              sigma_theta_pol_edge_1_radian
c         write(*,*)'sigma_theta_pol_edge_2_radian',
c     &              sigma_theta_pol_edge_2_radian

         exp1=dexp(-((theta_pol_radian-theta_pol_edge_1_radian)/
     &                sigma_theta_pol_edge_1_radian)**2)

         exp2=dexp(-((theta_pol_radian-theta_pol_edge_2_radian)/
     &                sigma_theta_pol_edge_2_radian)**2)
 
c        write(*,*)'exp1,exp2',exp1,exp2

c         write(*,*)'sigma_edgen_max_1',sigma_edgen_max_1
c         write(*,*)'sigma_edgen_max_2',sigma_edgen_max_2
c         write(*,*)'sigma_edgen_min',sigma_edgen_min

         del1=sigma_edgen_1-sigma_edgen_0

         del2=sigma_edgen_2-sigma_edgen_0

c         write(*,*)'del1,del2',del1,del2

         d_arg_exp1_d_theta=
     &         -2.d0*(theta_pol_radian-theta_pol_edge_1_radian)/
     &         sigma_theta_pol_edge_1_radian**2

         d_arg_exp2_d_theta=
     &         -2.d0*(theta_pol_radian-theta_pol_edge_2_radian)/
     &         sigma_theta_pol_edge_2_radian**2

         d2_arg_exp1_d2_theta=-2.d0/sigma_theta_pol_edge_1_radian**2
         d2_arg_exp2_d2_theta=-2.d0/sigma_theta_pol_edge_2_radian**2
       
c         write(*,*)'sigmedgn',sigmedgn

         sigma_edge_n=sigma_edgen_0+del1*exp1+del2*exp2

c         write(*,*)'sigma_edge_n',sigma_edge_n

         d_sigma_edge_n_d_theta_pol=
     &        +del1*exp1*d_arg_exp1_d_theta
     &        +del1*exp2*d_arg_exp2_d_theta

c         write(*,*)'d_sigma_edge_n_d_theta_pol',
c     &              d_sigma_edge_n_d_theta_pol

         d2_sigma_edge_n_d2_theta_pol=
     &        del1*exp1*(d_arg_exp1_d_theta**2+d2_arg_exp1_d2_theta)
     &       +del2*exp2*(d_arg_exp2_d_theta**2+d2_arg_exp2_d2_theta)
     
c        write(*,*)'d2_sigma_edge_n_d2_theta_pol',
c     &             d2_sigma_edge_n_d2_theta_pol

      endif

      if(i_edge_dens_anal.eq.2) then
c-------using spline ---------------
        if (first) then
c----------create spline coefficients
c
c----------check input data

c           write(*,*)'theta_pol_edge_dens_ar_degree(1)',
c     &     theta_pol_edge_dens_ar_degree(1)
          

           if (dabs(theta_pol_edge_dens_ar_degree(1)).gt.1.d-13) 
     &         then
               write(*,*)'in the input file genray.dat or genray.in'
               write(*,*)'in namelist edge_prof_nml'
               write(*,*)'theta_pol_edge_dens_ar_degree(1).ne.0'
               write(*,*)'theta_pol_edge_dens_ar_degree(1)=',
     &                    theta_pol_edge_dens_ar_degree(1)
               write(*,*)'Please set theta_pol_edge_dens_ar_degree(1)=0'
               write(*,*)'in the input file' 
               stop 'check input edge_prof_nml'
           endif

c           write(*,*)'n_pol_edge_dens',n_pol_edge_dens
c           write(*,*)'theta_pol_edge_dens_ar_degree(n_pol_edge_dens)',
c     &                theta_pol_edge_dens_ar_degree(n_pol_edge_dens)

           if(dabs(theta_pol_edge_dens_ar_degree(n_pol_edge_dens)-
     &              360.d0).gt.1.d-13) then
               write(*,*)'in the input file genray.dat or genray.in'
               write(*,*)'in namelist edge_prof_nml'
               write(*,*)'theta_pol_edge_dens_ar_degree(n_pol_edge_dens)
     &      .ne.360'
            write(*,*)'theta_pol_edge_dens_ar_degree(n_pol_edge_dens)=',
     &      theta_pol_edge_dens_ar_degree(n_pol_edge_dens)
               write(*,*)'Please set 
     &       theta_pol_edge_dens_ar_degree(n_pol_edge_dens)=360'
               write(*,*)'in the input file' 
               stop 'check input edge_prof_nml'
           endif

           do i=1,n_pol_edge_dens
             
             if(theta_pol_edge_dens_ar_degree(i).lt.0.d0) then
               write(*,*)'in the input file genray.dat or genray.in'
               write(*,*)'in namelist edge_prof_nml'
               write(*,*)'theta_pol_edge_dens_ar_degree(i).le.0'
               write(*,*)'at i=',i
               write(*,*)'theta_pol_edge_dens_ar_degree(i)=',
     &                    theta_pol_edge_dens_ar_degree(i)
              write(*,*)'Please change theta_pol_edge_dens_ar_degree(i)'
               write(*,*) 'in input file'
               stop 'check input edge_prof_nml'
             endif

             if(theta_pol_edge_dens_ar_degree(i).gt.360.d0) then
               write(*,*)'in the input file genray.dat or genray.in'
               write(*,*)'in namelist edge_prof_nml'
               write(*,*)'theta_pol_edge_dens_ar_degree(i).gt.360'
               write(*,*)'at i=',i
               write(*,*)'theta_pol_edge_dens_ar_degree(i)=',
     &                    theta_pol_edge_dens_ar_degree(i)
              write(*,*)'Please change theta_pol_edge_dens_ar_degree(i)'
               write(*,*) 'in the input file'
               stop 'check input edge_prof_nml'
             endif

             if (i.gt.1) then
c----------------check that the given poloidal angle array is increasing
                if(theta_pol_edge_dens_ar_degree(i-1).ge.
     &             theta_pol_edge_dens_ar_degree(i)) then
                   write(*,*)'in the input file genray.dat or genray.in'
                   write(*,*)'in namelist edge_prof_nml'
                   write(*,*)'it was given non-momotonic array'
                   write(*,*)'theta_pol_edge_dens_ar_degree(i))'
                   write(*,*)'at i=',i
              write(*,*)'Please change theta_pol_edge_dens_ar_degree(i)'
                   write(*,*) 'in the input file'
                   stop 'check input edge_prof_nml'
               endif
             endif !i.gt.1

             if (sigmedgn_ar(i).le.0.d0) then
                write(*,*)'in the input file genray.dat or genray.in'
                write(*,*)'in namelist edge_prof_nml'
                write(*,*)'sigmedgn_ar(i).le.0 at i=',i
                write(*,*)'sigmedgn_ar(i)=',sigmedgn_ar(i)
                write(*,*)'it shoud be positive'
                write(*,*)'Please change sigmedgn_ar(i)'
                write(*,*) 'in the input file'
                stop 'check input edge_prof_nml'
             endif

           enddo 

           do i=1,n_pol_edge_dens
               theta_pol_edge_dens_ar_radian(i)=
     &         theta_pol_edge_dens_ar_degree(i)*pi/180.d0
           enddo
c-------------------------------------------------------------------
c          creation of spline coefficients for sigmedgn
c-------------------------------------------------------------------     

           iop(1)=3 ! periodic spline boundary conditions
           iop(2)=3  

           call coeff1(n_pol_edge_dens,
     &                 theta_pol_edge_dens_ar_radian,
     &                 sigmedgn_ar,sigmedgn_deriv,
     &                 iop,1,workk)

c----test
c          itabl(1)=1
c          itabl(2)=1
c          itabl(3)=1
c          write(*,*)'theta_pol_edge_dens_ar_radian',
c     &               theta_pol_edge_dens_ar_radian
c          write(*,*)'sigmedgn_ar',sigmedgn_ar
c          write(*,*)'sigmedgn_deriv',sigmedgn_deriv

c          do i=1,n_pol_edge_dens
             
c             call terp1(n_pol_edge_dens,
c     &       theta_pol_edge_dens_ar_radian,
c     &       sigmedgn_ar,sigmedgn_deriv,
c     &       theta_pol_edge_dens_ar_radian(i),1,tabl,itabl)

c            sigma_edge_n=tabl(1)
c            d_sigma_edge_n_d_theta_pol=tabl(2)
c            d2_sigma_edge_n_d2_theta_pol=tabl(3)
c            write(*,*)'i,theta_pol_edge_dens_ar_radian(i)',
c     &                 i,theta_pol_edge_dens_ar_radian(i)
c            write(*,*)'sigma_edge_n',sigma_edge_n
c            write(*,*)'d_sigma_edge_n_d_theta_pol',
c     &                 d_sigma_edge_n_d_theta_pol
c            write(*,*)'d2_sigma_edge_n_d2_theta_pol',
c     &                 d2_sigma_edge_n_d2_theta_pol
c          enddo
c----end test

          first=.false. 
        endif ! if first: create spline coefficients
c
c-------calculate sigmedgen and its derivative by poloidal angle
c  

        itabl(1)=1
        itabl(2)=1
        itabl(3)=1

c        write(*,*)'theta_pol_edge_dens_ar_radian',
c     &             theta_pol_edge_dens_ar_radian
c        write(*,*)'sigmedgn_ar',sigmedgn_ar
c        write(*,*)'sigmedgn_deriv',sigmedgn_deriv
c        write(*,*)'theta_pol_radian',theta_pol_radian
 
        call terp1(n_pol_edge_dens,
     &  theta_pol_edge_dens_ar_radian,
     &  sigmedgn_ar,sigmedgn_deriv,
     &  theta_pol_radian,1,tabl,itabl)

c        write(*,*)'itabl',itabl
c        write(*,*)'tabl',tabl

        sigma_edge_n=tabl(1)
        d_sigma_edge_n_d_theta_pol=tabl(2)
        d2_sigma_edge_n_d2_theta_pol=tabl(3)  

      endif !i_edge_dens_anal.eq.2

 10   continue

      return
      end

      subroutine dens_rho_theta_LCFS(rho_small,theta_pol,i,
     &dens_rho_theta,d_dens_rho_theta_d_rho,
     &d_dens_rho_theta_d_theta)
c-----------------------------------------------------------------------
c     Calculate density and derivatives:d_dentsity/d_rho
c     d_dentsity/d_theta_poloidal
c     versus radius and poloidal angle outside LCFS 
c     1 < rho_small, 0=< theta_pol <2pi
c
c     for i plasma specie
c-----------------------------------------------------------------------
      implicit none

      include 'param.i'
      include 'one.i'
      include 'six.i'
      include 'edge_prof_nml.i'

c-----input
      real*8
     &rho_small, !small radius [non-dimensional]
     &theta_pol  !poloidal angle [radians] 0=< theta <2pi
                 ! =zero at outer side at the equatorial plane
                 ! =pi   at inner side at the equatorial plane
                 ! =pi/2 at the top of the polidal cross-section
!     &dens_min_edge   !minimal density it is set in 'edge_prof.nml'
!                 ! if (dens_rho_theta < dens_min_edge) then
!                 !    dens_rho_theta=dens_min_edge
!                 !    d_dens_rho_theta_d_theta=0

      integer i  ! number of plasma species i=1,..,nbulk
c-----output 
      real*8
     &dens_rho_theta,                !density
     &d_dens_rho_theta_d_rho,        !derivative d_density/d_rho
     &d_dens_rho_theta_d_theta       !derivative d_density/d_theta_poloidal

c-----locals
      real*8 rhoedge,densedge,d2_dens_drho(ndensa),
     & sigma_edge_n,d_sigma_edge_n_d_theta_pol,
     & d2_sigma_edge_n_d2_theta_pol,tabl(3),
     & densedge_e,dens_ratio

      integer itabl(3),k

      integer iwarn
      data iwarn/0/
      save iwarn

      rhoedge=1.d0 !small radius at LCFS

      itabl(1)=1
      itabl(2)=0
      itabl(3)=0 

c      write(*,*)'dense.f dens_rho_theta_LCFS ndens,i,nbulk',
c     & ndens,i,nbulk

      do k=1,ndens
         densm(k)=dens1(k,i)
         d2_dens_drho(k)=d2_dens_drho1(k,i)
      enddo

      call terp1(ndens,rhom,densm,d2_dens_drho,rhoedge,1,tabl,itabl)

      densedge=tabl(1) 

      call sigma_edge_n_theta_pol(theta_pol,
     &sigma_edge_n,d_sigma_edge_n_d_theta_pol,
     &d2_sigma_edge_n_d2_theta_pol)

      dens_rho_theta=densedge*dexp(-(rho_small-1.d0)/sigma_edge_n)

cSAP090222
c      write(*,*)'densedge,i,rho_small,sigma_edge_n,dens_rho_theta',
c     &           i,densedge,rho_small,sigma_edge_n,dens_rho_theta


      d_dens_rho_theta_d_theta = dens_rho_theta*
     &   ((rho_small-1.d0)/sigma_edge_n**2)*d_sigma_edge_n_d_theta_pol

      d_dens_rho_theta_d_rho =- dens_rho_theta/sigma_edge_n
     
     
cSAP090513
c      if (densedge.lt.dens_min_edge) then
c         write(*,*)'WARNINNG:in dens_rho_theta_LCFS'
c         write(*,*)'densedge.lt.dens_min_edge'
c         write(*,*)'It wil be dens_rho_theta=densedge'
c         dens_rho_theta=densedge 
c         d_dens_rho_theta_d_rho=0.d0
c         d_dens_rho_theta_d_theta=0.d0
c      else
c         if (dens_rho_theta.lt.dens_min_edge) then
c         write(*,*)'WARNINNG:in dens_rho_theta_LCFS'
c         write(*,*)'dens_rho_theta.lt.dens_min_edge'
c         write(*,*)'It wil be dens_rho_theta=dens_min_edge '
c            dens_rho_theta=dens_min_edge 
c            d_dens_rho_theta_d_theta=0.d0
c            d_dens_rho_theta_d_rho=0.d0
c         endif
c      endif

c------------------------------------------------------------
cSAP090526
      dens_ratio=1.d0

      if (i.gt.1) then
c----------------------------------------------------------
c        calculates the ratio density(i)/density_e at rho=1
c-----------------------------------------------------------
c        calculate the electron density at rho=1
         do k=1,ndens
            densm(k)=dens1(k,1)
            d2_dens_drho(k)=d2_dens_drho1(k,1)
         enddo
         call terp1(ndens,rhom,densm,d2_dens_drho,rhoedge,1,tabl,itabl)
         densedge_e=tabl(1) 
         dens_ratio=densedge/densedge_e
      endif

cSAP090526     
c      if (dens_rho_theta.lt.dens_min_edge) then
      if (dens_rho_theta.lt.dens_min_edge*dens_ratio) then
         if (iwarn.lt.50) then
         iwarn=iwarn+1   
         write(*,*)'WARNINNG:in dens_rho_theta_LCFS'
         write(*,*)'densedge.lt.dens_min_edge*dens_ratio'
         write(*,*)'It wil be dens_rho_theta=dens_min_edge*dens_ratio'
         endif
cSAP090526
c         dens_rho_theta=dens_min_edge
         dens_rho_theta=dens_min_edge*dens_ratio
         d_dens_rho_theta_d_rho=0.d0
         d_dens_rho_theta_d_theta=0.d0
      endif      

      return
      end
!
c
c        *********************dense ***************************
c        *                        -                           
c        * this function calculates the density profile       
c        * as a function  from psi(z,r,phi) 		      
c        *  rho=psi-psimag 				      
c        *  a=psilim-psimag inside the function b()           
c        *  these variables are in common 'one'
c        *                                                    
c        * This function uses 2D spline at RZ mesh outside LCFS
c        * which has density fall near the chamber wall.      
c        *
c        * For n_wall=0 case this function should give density
c        * close to the density created by
*        * FUNCTION dense_no_RZ_spline(z,r,phi,i)
c        ******************************************************
c
c------------------------------------------------------------------
c								   !
c        input parameters					   !
c      								   !
c      z, r, phi - point which the components of magnetic field	   !
c                  are calculated in.				   !
c       rho is from common/one/                                    !
c------------------------------------------------------------------
c      uses
c         constants dense0,denseb,rn1de,rn2de,idens are in common 'one'
c         they were read in subroutine dinit from the file genray.in
c         rho is  in common 'one'
c      functions
c         _densrho_(rho)    finds dense from spline
c------------------------------------------------------------------
      real*8 FUNCTION dense(z,r,phi,i)
c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none
      include 'param.i'
      include 'one.i'
c-----input
      real*8 z,r,phi!space coordinates
c     rho is in common  /one/
      integer i !plasma species number
c-----local
      real*8  theta_pol, !0=< theta_pol <2pi
     &rho_small,dens_rho_theta,
     &d_dens_rho_theta_d_rho,d_dens_rho_theta_d_theta
c-----externals 
      real*8 density_r_z_i
      real*8 thetapol,   ! -pi <thetapol =<pi
     & vardens,densrho


      if(ixyz.eq.1) then
        stop 'dense(z,r,phi,i) should not be used for ixyz=1'
      endif
      
      if(i.gt.nbulk)then
	stop
      endif
      if (rho.gt.1.d0-1.d-10) then
c----------------------------------------------------------
c         density outside LCFS
c----------------------------------------------------------
          if(n_wall.gt.1) then
c-----------------------------------------------------------
c           calculate density using 2D spline at RZ mesh
c-----------------------------------------------------------
            dense=density_r_z_i(z,r,phi,i)
          else
c----------------------------------------------------------
c----------------------------------------------------------------
c           calculate density using formula versus small radius
c           and poloidal angle
c----------------------------------------------------------------
            theta_pol=thetapol(z,r) ! -pi <thetapol =<pi
            if (theta_pol.lt.0d0) then
              theta_pol=theta_pol+2*pi !pi< theta_pol<2pi
            endif

            rho_small=rho
   
            call dens_rho_theta_LCFS(rho_small,theta_pol,i,
     &        dens_rho_theta,d_dens_rho_theta_d_rho,
     &        d_dens_rho_theta_d_theta)
            dense = dens_rho_theta 

         endif
      else
         dense=densrho(rho,i)
      endif
      return
      END


