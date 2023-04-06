c        **********************cninit**************************
c        *                        -                           *
c        * It solves the dispersion relation N=N(n_par)       *
c        * Then subroutine calculates the initial components  *
c        * of the refractive index N: cnz=N_z,cnr=N_r,cm=M    *
c        ******************************************************
c
c------------------------------------------------------------------
c					                          !
c        input parameters				          !
c        z,r,phi,cnpar,cnteta,cnphi 		       	          !
c        if i_n_poloidal=4 cnpar=N_parallel will be calculated
c        inside this subroutine using cnteta,cnphi
c        istart from common block/ one/                           !
c        output parameters				          !
c        u(4)=cnz,u(5)=cnr,u(6)=cm 			          !
c        iraystop=1 is switch to stop the ray calculation         !
c------------------------------------------------------------------
c        it uses the following functions and subroutines          !
c        ias1r,  b ,y,x,gamma1,s,abc,hamilt1,                     !
c------------------------------------------------------------------
      subroutine cninit(z,r,phi, cnpar,cnteta,cnphi,
     1                  cnz,cnr,cm, iraystop)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'      
      include 'grill.i'
      
      complex*16 hamilt_c !complex hamiltonian used for id=10
      save
      double complex hotnperp,cmplnper,relativistic_nperp
      double precision b
      external hotnperp,relativistic_nperp,b      

c-----for hot plasma roots
      double precision 
     &N_perp_root_ar(n_hot_roots_a)     !hot plasma roots

      double complex K(3,3) ! YuP: added

      iraystop=0
      pi=4*datan(1.d0)
      cnpar2=cnpar*cnpar
      cntang2=cnteta*cnteta+cnphi*cnphi
      bmod=b(z,r,phi)

      if(istart.eq.2) then !grill conditions
      
        if (i_n_poloidal.eq.3) then
c----------refractive index is specified by N_parallel 
c          and ksi_nperp the angle between grad(psi) and_ N_perpendicular
c
c          set zero values for following variables
c          to get the solution of the dispersion function N_perp(N_parallel)
c          In this case N_perp=N
           cnteta=0.d0
           cnphi=0.d0
           cntang2=0.d0
        endif

        if (i_n_poloidal.eq.4) then
c---------refractive index is specified by N_toroidal and N_poloidal
c         calculation of the parallel refractive index component
c         N_parallel=cnpar
          gradpsi=dsqrt(dpdrd*dpdrd+dpdzd*dpdzd)
          b_teta=(bz*dpdrd-br*dpdzd)/gradpsi !poloidal magnetic field
          cnpar=(cnphi*bphi+cnteta*b_teta)/bmod !N_parallel
          write(*,*)'cninit.f cnteta,cnphi,cnpar',cnteta,cnphi,cnpar
        endif

      endif !istart.eq.2
c-----------------------------------------------------------------
      iroot=0
c------------------------------------------------------------------
c     epsmode is the accuracy  for mode selection (O or X) on
c             the plasma boundary
c------------------------------------------------------------------
      epsmode=0.0000001d0
c--------------------------------------------------------------------
c     the initial condition in plasma boundary
c     for cn2 from cnpar by using dispersion relation
c     f*n**4 +g*n**2+w=0
c     in the following form n**2=(-g+ioxm*sqrt(g**2-4*f*w))/(2*f)
c------------------------------------------------------------------
c------------------------------------------------------------------
c     Appleton - Hartree dispersion relation
c
c if 1
      if (id.eq.3) then

         call cninit3(z,r,phi,cnpar,cnteta,cnphi,
     1                  cnz,cnr,cm,iraystop)
         cnper=dsqrt(cnz**2+cnr**2+(cm/r)**2-cnpar**2)
         goto 111
      end if !id=3
c end if 1
c-----------------------------------------------------------------
c     initial condition in plasma boundary
c     for cn2 from cnpar by using dispersin relation
c     f*n**4 +g*n**2+w=0
c     in the following form n**2=(-g+ioxm*sqrt(g**2-4*f*w))/(2*f)
c------------------------------------------------------------------
c     cold plasma dispresion relation
c------------------------------------------------------------------
c if 0
      if ((id.eq.1).or.(id.eq.2)) then

c         write(*,*)'cninit before cninit12 z,r,phi,cnpar,cnteta,cnphi',
c     &   z,r,phi,cnpar,cnteta,cnphi

         call cninit12(z,r,phi,cnpar,cnteta,cnphi,
     1                  cnz,cnr,cm,iraystop)
         if (iraystop.eq.1)then
           write(*,*)'the given conditions did not give root' 
           return
         endif
         cnper=dsqrt(cnz**2+cnr**2+(cm/r)**2-cnpar**2)
         goto 111
      end if !id=1 or id=2
c end if 0

c if 4    

      if(id.eq.6) then
c-----------------------------------------------------------------
c       initial condition on plasma boundary
c       for cnper2 from cnpar by using the hot dispersion relation
c-----------------------------------------------------------------
        if (i_look_roots.eq.2) then     
          call calculate_hot_nperp_roots(z,r,phi,cnpar,
     &    n_hot_roots,N_perp_root_ar)
          cnper=N_perp_root_ar(k_hot_root)
          cnprim=0.d0
        else
          ihermloc=iherm
          write(*,*)'cninit.f before hotnperp'
          cmplnper=hotnperp(z,r,phi,cnpar,cnteta,cnphi,K,iraystop) ! YuP: K was not decl.
          write(*,*)'cninit.f after hotnperp cmplnper',cmplnper
          if (iraystop.eq.1)then
           write(*,*)'the given conditions did not give hot root' 
          return
          endif
          cnper=dreal(cmplnper)
          cnprim=dimag(cmplnper)            
          write(*,*)'in cninit aft hotnperp cnper,cnprim',cnper,cnprim
        endif

      endif !id=6,8,9

c-----------------------------------------------------------------------    

      if (id.eq.14) then
c--------------------------------------------------------------------
c         id =11 Eric Nelson-Melby relativistic tensor
c         Dispersion function = Re(Det) 
c         id=14 Abhay Ram relativistic electron dispersion function
c----------------------------------------------------------------------
c       initial condition by using dispersion relation
c       which was used in "Mazzucato's" code 
c       and double complex function relativistic_nperp
c----------------------------------------------------------------------
        cmplnper=relativistic_nperp(z,r,phi,cnpar,cnteta,cnphi,K,
     &                              iraystop)
        cnper=dreal(cmplnper)
        cnprim=dimag(cmplnper)

        write(*,*)'cninit id14 z,r,phi,cnpar,cnteta,cnphi,cnper,cnprim',
     &  z,r,phi,cnpar,cnteta,cnphi,cnper,cnprim
      endif ! (id.eq.11).or.(id.eq.14)
c----------------------------------------------------------------        

c--------------------------------------------------------------------
c-----calculations of initial values cnz,cnr,cm from cnper for id=4,5,6,7    
      cn2=cnper**2+cnpar**2   
      cnrho2=cn2-cnteta**2-cnphi**2
      if (cnrho2.lt.0.d0) cnrho2=1.d-16
      cnrho=dsqrt(cnrho2)

      write(*,*)'in cninit before cnzcnr cnper,cnpar,cnrho',
     & cnper,cnpar,cnrho
      write(*,*) 'cn=',dsqrt(cn2)
c---------------------------------------------------------------------
      call cnzcnr(z,r,phi,cnteta,cnphi,cnrho,cnz,cnr,cm)
c---------------------------------------------------------------------
      write(*,*)'in cninit after cnzcnr cnz,cnr,cm',cnz,cnr,cm
      cm=cnphi*r
      gam=gamma1(z,r,phi,cnz,cnr,cm)
      write(*,*)'cninit.f before ddd=hamilt1'
      ddd=hamilt1(z,r,phi,cnz,cnr,cm)
      write(*,*)'in cninit ddd=',ddd


  111 continue

c----------------------------------------------------------------
c     for the grill conditions it will use the different grill types
c---------------------------------------------------------------
      write(*,*)'cninit istart,i_n_poloidal',istart,i_n_poloidal
      if(istart.eq.2) then !grill conditions
         if (i_n_poloidal.eq.2) then !input N_parallel, N_poloidal
c            calculate N_phi,N_theta,N_rho
c-----------------------------------------------------------------------
c            dpdrd=dpsidr dpdzr=d psidr were calculated by b(z,r,phi)
c            They are in  common/one/
c---------------------------------------------------------------------
             gradpsi=dsqrt(dpdrd*dpdrd+dpdzd*dpdzd)
             b_teta=(bz*dpdrd-br*dpdzd)/gradpsi !poloidal magnetic field
             cnteta_full=n_theta_pol
             write(*,*)'cninit i_n_poloidal,n_theta_pol',
     &       i_n_poloidal,n_theta_pol
             write(*,*)'cnteta_full',cnteta_full

             alpha_teta=(cnteta_full-cnpar*b_teta/bmod)/cnper
             write(*,*)'cnpar,b_teta,bmod,cnper,alpha_teta',
     &       cnpar,b_teta,bmod,cnper,alpha_teta
     
             cnphi_full=cnpar*bphi/bmod-alpha_teta*cnper*b_teta/bphi
             write(*,*)'cnphi_full',cnphi_full
c_test
             cnteta_full=cnpar*b_teta/bmod+cnper*alpha_teta
             write(*,*)'2 cnteta_full=',cnteta_full

             arg= 1.d0-(alpha_teta*bmod/bphi)**2
             if (arg.lt.0.d0) then
               write(*,*)'cninit.f 1.d0-(alpha_teta*bmod/bphi)**2<0'
               write(*,*)'cninit.f change n_theta_pol'
               stop
             endif

             cnrho_full=cnper*dsqrt(arg)
             write(*,*)'cnteta_full,cnphi_full,cnrho_full,cn**2',
     &       cnteta_full,cnphi_full,cnrho_full,
     &       cnteta_full**2+cnphi_full**2+cnrho_full**2

         endif !i_n_poloidal=2

         if (i_n_poloidal.eq.4) then !input N_toroidal and N_poloidal
c-------------calculate N_phi,N_theta,N_rho
              cnteta_full=cnteta
              cnphi_full=cnphi
              write(*,*)'cninit i_n_poloidal=4 cnteta,cnphi,cnper,cnpar'
     &                   ,cnteta,cnphi,cnper,cnpar
     
              if (cnper.ne.0.d0) then
                 sin_ksi=(cnphi*b_teta-cnteta*bphi)/(bmod*cnper)
                 if (dabs(sin_ksi).le.1.d0) then
                    cos_ksi=dsqrt(1.d0-sin_ksi**2)
                    cnrho_full=cnper*cos_ksi*i_vgr_ini
                 else
                    write(*,1010)
 1010               format('in cninit.f i_n_poloidal=4 case',/,
     &              'dabs(sin_ksi)>1,it is imposible to find N_rho')
                    iraystop=1
                    return                
                 endif
              else
                 if ((cnphi/bphi).ne.(cnteta/b_teta)) then
                    write(*,1000)
 1000               format('in cninit.f i_n_poloidal=4 case',/,
     &             'cnper=0 and it is impossible to find n_parallel')
                   iraystop=1
                   return
                 else
                   cnrho_full=0.d0
                 endif
              endif
              
              write(*,*)'cninit.f i_n_poloidal=4 cnrho_full',cnrho_full
         endif ! i_n_pol=4

         if((i_n_poloidal.ne.1).and.(i_n_poloidal.ne.3)) then
             call cnzcnr(z,r,phi,cnteta_full,cnphi_full,cnrho_full,
     &       cnz,cnr,cm)

             write(*,*)'cninit end grill condition'
             write(*,*)'cnteta_full,cnphi_full,cnrho_full,cn**2',
     &       cnteta_full,cnphi_full,cnrho_full,
     &       cnteta_full**2,cnphi_full**2+cnrho_full**2
             write(*,*)'cnz**2+cnr**2+(cm/r)**2',
     &       cnz**2+cnr**2+(cm/r)**2
             write(*,*)'cnz,cnr,cm',cnz,cnr,cm
             write(*,*)'(cnz*bz+cnr*br+cm*bphi/r)/bmod',
     &        (cnz*bz+cnr*br+cm*bphi/r)/bmod
         endif

      endif !grill conditions


      return
      end
      
      
      
      
      
c_____________________________________________________________
c        **********************cnzcnr**************************
c        This subroutine calculates the initial value		 *
c        cnz,cnr,cm					         *
c        It directs the wave into or out the plasma              *
c        the input parameter i_vgr_ini is in common /one/        *
c        i_vgr_ini =+1 the wave is directed into the plasma      *
c                      (in the initial point)                    *
c                  =-1 the wave is directed out the plasma       *
c-----------------------------------------------------------------
c      * input parameters:z,r,phi ,cnteta,cnphi,cnrho	      	 *
c      *                  cnrho is directed inside the plasma	 *
c________________________________________________________________
c     * output parameters: cnz,cnr,cm -refractive index components*
c-----------------------------------------------------------------
      subroutine cnzcnr(z,r,phi,cnteta,cnphi,cnrho, cnz,cnr,cm)
      implicit integer (i-n), real*8 (a-h,o-z)
c-----------------------------------------------------------------
c	           z       r      phi
c       e_phi= {   0   ,   0    ,   1    }
c       e_psi= { dpsidz, dpsidr ,   0    }/mod(grad(psi))
c       e_teta={ dpsidr,-dpsidz ,   0    }/mod(grad(psi))
c
c-----------------------------------------------------------------
      include 'param.i'
      include 'one.i'
      include 'grill.i'
      dimension u(6),deru(6)
c---------------------------------------------------------------------
c     dpdrd=dpsidr dpdzr=dpsidr were calculated by b(z,r,phi)
c     They are in  common/one/
c---------------------------------------------------------------------
      gradpsi=dsqrt(dpdrd*dpdrd+dpdzd*dpdzd)
c-------------------------------------------------------------------
c     the initialization of cnrho 
      cirho=1.d0
c------------------------------------------------------------------
10    continue
      write(*,*)'cninit.f in cnzcnr cnteta,cirho,cnrho',
     &cnteta,cirho,cnrho
c     write(*,*)'cninit.f in cnzcnr dpdzd,dpdrd,gradpsi',
c    &dpdzd,dpdrd,gradpsi
      cnz=(cnteta*dpdrd-cirho*cnrho*dpdzd)/gradpsi
      cnr=(-cnteta*dpdzd-cirho*cnrho*dpdrd)/gradpsi
      cm=cnphi*r
      write(*,*)'cninit.f in cnzcnr cnz,cnr,cm',cnz,cnr,cm
cSm030515
c      if((i_n_poloidal.eq.1).or.(i_n_poloidal.eq.2)) then
cSm040415
c     I tried it for id=4 ECR case when rside1 used the derivatives
c     by the trajectory d/dl (not by time d/dt).
c     In some case d/dl runs the ray to the plasma edge.
c     As I understand this changing of the cnrho sign is not correct
c     fo EC launch. It changes the initial NR and NZ signes. 

      if(((i_n_poloidal.eq.1).or.(i_n_poloidal.eq.2)).
     &    or.(istart.eq.1)) then 

  
c*****************************************************************
c        determination cirho to obey the situation
c        when groop velocity is directed inside or outside  the plasma
c        in the initial point
c-----------------------------------------------------------------
         u(1)=z
         u(2)=r
         u(3)=phi
         u(4)=cnz
         u(5)=cnr
         u(6)=cm
         call rside1(u,deru) 
c----------------------------------------------------------------
c        cmultpl is the scalar multiplication V_groop*grad(psi)
c        gradient(psi) is directed outside the plasma
c----------------------------------------------------------------
         cmultpl=dpdzd*deru(1)+dpdrd*deru(2)
        write(*,*)'cninit in cnzcnr cmultpl,i_vgr_ini',cmultpl,i_vgr_ini

         if ((cmultpl*i_vgr_ini).gt.0.d0) then
c-----------the poloidal direction of the group velocity is opposite
c           to the direction determined by the parameter i_vgr_ini     
c           We change the sign of the poloidal group velocity 
            cirho=-1.d0
            go to 10
         end if
c****************************************************************
c         end irho determination
      endif !i_n_poloidal =1 or =2
c----------------------------------------------------------------
      write(*,*)'cninit.f at end, cnzcnr cnz,cnr,cm:',cnz,cnr,cm
      return
      end





c        **********************npernpar************************
c        *                        -                           *
c        * It solves the dispersion relation Nper=N(n_par)    *
c        * Then subroutine calculates the initial components  *
c        * of the refractive index  cnz,cnr,cm                *
c        ******************************************************
c
c---------------------------------------------------
c     							   !
c     input parameters					   !
c     z,r,phi,cnpar 	                                   !
c     output parameters:cnper2p,cnper2m	                   !
c---------------------------------------------------
c     it uses the following functions and subroutines      !
c     ias1r,b   ,y,x,gamma1,s,abc,hamilt1,                 !
c-----------------------------------------------------------
      subroutine npernpar(z,r,phi,cnpar,cnper2p,cnper2m)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'

      pi=4*datan(1.d0)
      cnpar2=cnpar*cnpar
      cnpar4=cnpar2*cnpar2
      bmod=b(z,r,phi)
c-------------------------------------------------------------
c     calculations of  cnper2p and cnper2m
c     from cnpar2 by using the dispersin relation
c     f*cnper**4 +(2f*cnpar**2+g)*cnper**2+(f*cnpar**4+g*cnpar**2+w)=0
c     f*cnper**4 +gnew*cnper**2+wnew=0
c     gnew=2fcnpar**2+g, wnew=f*cnpar**4+g*cnpar**2+w
c     in the following form cnper**2=(-gnew+,-*sqrt(gnew**2-4*f*wnew))/(2*f)
c------------------------------------------------------------------
c     Appleton - Hartree dispersion relation
c
c if 1
      if (id.eq.3) then
         xi=x(z,r,phi,1)
         yi=y(z,r,phi,1)
         py2=yi*yi
         py4=py2*py2
         px=1.-xi
         px2=px*px
c------------------------------------------------------------------
         pyp=xi/(1.+yi)
         f1e=1.d0
         f0e=-pyp
         g1e=cnpar2*(-xi)+pyp+xi-2.
         g0e=cnpar2*pyp+xi*(1.-pyp)+pyp*(1.-xi)
         w1e=cnpar2*(xi-pyp)+(1.-pyp)*(1.-xi)
         w0e=cnpar2*(pyp*(1-xi)-xi*(1.-pyp))-(1.-xi)*(1-pyp)*xi
         dele=1.-yi
         fd=f1e*dele+f0e
         gd=g1e*dele+g0e
         wd=w1e*dele+w0e
c new coefficients
         gnew=gd+2.d0*fd*cnpar2
         wnew=wd+gd*cnpar2+fd*cnpar4
         gd=gnew
         wd=wnew

         detin=gd**2-4.*fd*wd
         if (detin.lt.0d0) then
            write(*,*)' 1 in npernpar detin  less then zero '
            return
         endif
         cnper2p=(-gd+dsqrt(detin))/(2.*fd)
         cnper2m=(-gd-dsqrt(detin))/(2.*fd)
c         WRITE(*,*)'Aplt cnpernpar cnper2p,cnper2m',cnper2p,cnper2m
      end if
c end if 1
c------------------------------------------------------------------
c     cold plasma dispresion relation
c------------------------------------------------------------------
c if 0
      if ((id.eq.1).or.(id.eq.2)) then
        call s(z,r,phi,s1,s2,s3,s4,s6,s7)

c       ib=1 electron resonance condition may be
c  if 2
        if (ib.eq.1) then
c         write(*,*)'cnint cold plasma ib=1 '
          xe=x(z,r,phi,1)
          ye=y(z,r,phi,1)
          pyp=xe/(1.d0+ye)
          dele=1.d0-ye
          f1e=s7
          f0e=-pyp
          g1e=cnpar2*(s4-s7)-(s6-pyp)*s3-s4*s7
          g0e=cnpar2*pyp+xe*(s6-pyp)+pyp*s4
          w1e=cnpar2*(-s7*s4+(s6-pyp)*s3)+(s6-pyp)*s3*s4
          w0e=cnpar2*(pyp*s4-xe*(s6-pyp))-s4*(s6-pyp)*xe
          fd=f1e*dele+f0e
          gd=g1e*dele+g0e
          wd=w1e*dele+w0e
c new coefficients
          gnew=gd+2.d0*fd*cnpar2
          wnew=wd+gd*cnpar2+fd*cnpar4
          gd=gnew
          wd=wnew

          detin=gd**2-4.d0*fd*wd
          if (detin.lt.0d0) then
             write(*,*)' 2 in npenpar detin  less then zero'
             cnper2p=-1.d0
             cnper2m=-1.d0
             return
	  end if
          cnper2p=(-gd+dsqrt(detin))/(2.d0*fd)
          cnper2m=(-gd-dsqrt(detin))/(2.d0*fd)
c          write(*,*)'in cninit.f  1 npernpar cnper2p,cnper2m',
c     .    cnper2p,cnper2m
          goto 111
        end if
c end if 2
c
c     ib.gt.1 iones resonance condition may be
c  if 3
        if (ib.gt.1) then
c          write(*,*)'cold plasma ib .gt.1  '
           xb=x(z,r,phi,ib)
           yb=y(z,r,phi,ib)
           xe=x(z,r,phi,1)
           ye=y(z,r,phi,1)
           pype=xe/(1.d0+ye)
           pyme=xe/(1.d0-ye)
           pyme2=pype/(1.d0-ye)
           pypb=xb/(1.d0+yb)
           delib=1.d0-yb
	   f1b=(s1-pyme2)
	   f0b=-pypb

	   g1b=cnpar2*(s4-(s1-pyme2))-(s2-pype)*(s3-pyme)-s4*(s1-pyme2)
           g0b=cnpar2*pypb+xb*(s3-pyme)+pypb*s4

	   w1b=cnpar2*(-s4*(s1-pyme2)+(s2-pype)*(s3-pyme))+
     1	     s4*(s2-pype)*(s3-pyme)
           w0b=cnpar2*(pypb*s4-xb*(s3-pyme))-s4*(s3-pyme)*xb
           fd=f1b*delib+f0b
           gd=g1b*delib+g0b
           wd=w1b*delib+w0b
c new coefficients
           gnew=gd+2.d0*fd*cnpar2
           wnew=wd+gd*cnpar2+fd*cnpar4
           gd=gnew
           wd=wnew

           detin=gd**2-4.d0*fd*wd
           if (detin.lt.0d0) then
              write(*,*)' 3 in dinit detin  less then zero '
              return
           end if
           cnper2p=(-gd+dsqrt(detin))/(2.d0*fd)
           cnper2m=(-gd-dsqrt(detin))/(2.d0*fd)
c        write(*,*)'in npernpar ib.qt.1 cnper2p,cnper2m',cnper2p,cnper2m
           goto 111
        end if
c end if 3
      end if
c end if 0
  111 continue

      return
      end
      
      
      
      
      
c        **********************cninit3*************************
c        *                        -                           *
c        * It solves the dispersion relation N=N(n_par)       *
c        * for Appleton-Hartree disperstion relation
c        * Then subroutine calculates the initial components  *
c        * of the refractive index  cnz,cnr,cm        	      *
c        ******************************************************
c
c------------------------------------------------------------------
c					                          !
c        input parameters				          !
c        z,r,phi,cnpar,cnpar,cnphi 		       	          !
c        istart from common block/ one/                           !
c        output parameters				          !
c        u(4)=cnz,u(5)=cnr,u(6)=cm 			          !
c        iraystop=1 end ray calculation                           !
c------------------------------------------------------------------
c        it uses the following functions and subroutines          !
c        ias1r,b   ,y,x,gamma1,s,abc,hamilt1,                     !
c------------------------------------------------------------------
      subroutine cninit3(z,r,phi,cnpar,cnteta,cnphi,
     1                  cnz,cnr,cm,iraystop)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'

      save
      double complex cmplnper      

      iraystop=0
      pi=4*datan(1.d0)
      cnpar2=cnpar*cnpar
      cntang2=cnteta*cnteta+cnphi*cnphi
      bmod=b(z,r,phi)

c-----------------------------------------------------------------
      iroot=0
c------------------------------------------------------------------
c     epsmode is the accuracy  for selection mode (O or X) on
c             the plasma boundary
c------------------------------------------------------------------
      epsmode=0.0000001d0
c--------------------------------------------------------------------
c     initial condition in plasma boundary
c     for cn2 from cnpar by using dispersion relation
c     f*n**4 +g*n**2+w=0
c     in the following form n**2=(-g+ioxm*sqrt(g**2-4*f*w))/(2*f)
c------------------------------------------------------------------
c------------------------------------------------------------------
c     Appleton - Hartree dispersion relation
c
      xi=x(z,r,phi,1)
      yi=y(z,r,phi,1)
      py2=yi*yi
      py4=py2*py2
      px=1.d0-xi
      px2=px*px
c------------------------------------------------------------------
      pyp=xi/(1.d0+yi)
      f1e=1.d0
      f0e=-pyp
      g1e=cnpar2*(-xi)+pyp+xi-2.d0
      g0e=cnpar2*pyp+xi*(1.-pyp)+pyp*(1.d0-xi)
      w1e=cnpar2*(xi-pyp)+(1.d0-pyp)*(1.d0-xi)
      w0e=cnpar2*(pyp*(1.d0-xi)-xi*(1.d0-pyp))-(1.d0-xi)*(1-pyp)*xi
      dele=1.d0-yi
      fd=f1e*dele+f0e
      gd=g1e*dele+g0e
      wd=w1e*dele+w0e
      detin=gd**2-4.d0*fd*wd
      if (detin.lt.0d0) then
	 write(*,*)' 3 in cninit detin  less then zero '
	 iraystop=1
	 return
      end if

c     cn2=(-gd+ioxm*dsqrt(detin))/(2.d0*fd)
      cn2p=(-gd+dsqrt(detin))/(2.d0*fd)
      cn2m=(-gd-dsqrt(detin))/(2.d0*fd)
      WRITE(*,*)'apl cninit cn2p,cn2m',cn2p,cn2m
      if((cn2m.lt.0.d0).and.(cn2p.lt.0d0)) then
            write(*,*)'in cninit2 two roots of the dispersion < 0'
            write(*,*)'cn2m,cn2p',cn2m,cn2p
            write(*,*)'the given wave can not exist in plasma'
	    iraystop=1
	    return
      endif
10    iroot=iroot+1
      if(iroot.eq.1) then
        cn2=cn2p
	if(cn2.lt.cntang2)then
           write(*,*)'in cninit2 cn2p.lt.cntang2'
           go to 10
        end if
      else
        cn2=cn2m
        if(cn2.lt.cntang2)then
           write(*,*)'in cninit2 cn2m.lt.cntang2'
           write(*,*)'the given mode can not exist in plasma'
           iraystop=1
           return
        end if
      endif
c     write(*,*)'in cninit cn2=',cn2
      cnrho2=cn2-cnteta**2-cnphi**2
      cnrho=dsqrt(cnrho2)
      cnperp=cnrho
cSm050826
      cnperp=dsqrt(cn2-cnpar2)
c------------------------------------------------------------------
      call cnzcnr(z,r,phi,cnteta,cnphi,cnrho,cnz,cnr,cm)
c-------------------------------------------------------------------
      cm=cnphi*r
      gam=gamma1(z,r,phi,cnz,cnr,cm)
      ds=dsin(gam)
      dc=dcos(gam)
      ds2=ds*ds
      dc2=dc*dc
      ds4=ds2*ds2
c--------------------------------------------------------------------
c     control that cn2 and gam are the solution of the dispersion
c     relation  n**2=(-b+ioxm*sqrt(b**2-4*a*c))/(2*a)
c--------------------------------------------------------------------
      sqrdet=dsqrt(py4*ds4+4.d0*py2*px2*dc2)
      pz=2.d0*px-py2*ds2+ioxm*sqrdet
      cn2new=1.d0-2.d0*xi*px/pz
c     write(*,*)'cn2new=',cn2new
      if (iroot.eq.1) then
	dnp=dabs(cn2-cn2new)
	if (dnp.gt.epsmode)then
	   goto 10
	end if
      else
	dnm=dabs(cn2-cn2new)
	if (dnm.gt.epsmode)then
           write(*,*)'the given mode can not exist in plasma'
	   iraystop=1
	   return
	end if
      end if
      goto 111
      
  111 continue
      return
      end
      
      
      

c        **********************cninit12_n_gam******************
c        *                        -                           *
c        * It solves the dispersion relation 
c        & N=N(n_par)=N(gam,ioxm)                             *
c        * for cold plasma for given ioxm                     *
c        * Then subroutine calculates the initial components  *
c        * of the refractive index  cnz,cnr,cm         	      *
c        ******************************************************
c
c------------------------------------------------------------------
c					                          !
c        input parameters				          !
c        z,r,phi,cnpar,cnteta,cnphi 		       	          !
c        istart from common block/ one/                           !
c        output parameters				          !
c        u(4)=cnz,u(5)=cnr,u(6)=cm 			          !
c        iraystop=1 end ray calculation                           !
c------------------------------------------------------------------
c        it uses the following functions and subroutines          !
c        ias1r,b   ,y,x,gamma1,s,abc,hamilt1,                     !
c------------------------------------------------------------------
      subroutine cninit12_n_gam(z,r,phi,cnpar,cnteta,cnphi,
     1                  cnz,cnr,cm,iraystop)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'

      save
      double complex cmplnper      
      
      iraystop=0
      pi=4*datan(1.d0)
      cnpar2=cnpar*cnpar
      cntang2=cnteta*cnteta+cnphi*cnphi
      bmod=b(z,r,phi)
c      write(*,*)'cninit12_n_gam ioxm=',ioxm
c      write(*,*)'z,r,phi,cnpar,cnteta,cnphi',z,r,phi,cnpar,cnteta,cnphi
c-----------------------------------------------------------------
      iroot=0
c------------------------------------------------------------------
c     epsmode is the accuracy  for selection mode (O or X) on
c             the plasma boundary
c------------------------------------------------------------------
      epsmode=0.0000001d0
      epsmode=1.d-7
c--------------------------------------------------------------------
c     initial condition in plasma boundary
c     for cn2 from cnpar by using dispersion relation
c     f*n**4 +g*n**2+w=0
c     in the following form n**2=(-g+ioxm*sqrt(g**2-4*f*w))/(2*f)
c------------------------------------------------------------------
c     cold plasma dispresion relation
c------------------------------------------------------------------
c if 0
c      write(*,*)'cninit12 z,r,phi',z,r,phi
      if ((id.eq.1).or.(id.eq.2)) then
        call s(z,r,phi,s1,s2,s3,s4,s6,s7)

c       ib=1 electron resonance condition may be
c  if 2
        if (ib.eq.1) then
          xe=x(z,r,phi,1)
          ye=y(z,r,phi,1)
c          write(*,*)'cninit12 xe,ye',xe,ye
	    pyp=xe/(1.d0+ye)
	    dele=1.d0-ye
	    f1e=s7
	    f0e=-pyp
	    g1e=cnpar2*(s4-s7)-(s6-pyp)*s3-s4*s7
	    g0e=cnpar2*pyp+xe*(s6-pyp)+pyp*s4
	    w1e=cnpar2*(-s7*s4+(s6-pyp)*s3)+(s6-pyp)*s3*s4
            w0e=cnpar2*(pyp*s4-xe*(s6-pyp))-s4*(s6-pyp)*xe
	    fd=f1e*dele+f0e
	    gd=g1e*dele+g0e
	    wd=w1e*dele+w0e
	    detin=gd**2-4.d0*fd*wd
c            write(*,*)'cninit12 fd,gd,wd,detin',fd,gd,wd,detin
	    if (detin.lt.0d0) then
	       write(*,*)' 2 in dinit detin  less then zero '
	       iraystop=1
	       return
	    end if
c	    cn2=(-gd+ioxm*dsqrt(detin))/(2.d0*fd)
	    cn2p=(-gd+dsqrt(detin))/(2.d0*fd)
	    cn2m=(-gd-dsqrt(detin))/(2.d0*fd)
c	    write(*,*)'in cninit fd,gd,wd,detin'
	    write(*,*)fd,gd,wd,detin
	    write(*,*)'in cninit12 cn2p,cn2m',cn2p,cn2m
            write(*,*)'cninit12 cn2p-cnpar2,cn2m-cnpar2',
     .      (cn2p-cnpar2),(cn2m-cnpar2)

ctest angle
            if (cn2p.gt.0.d0)then
               if(cnpar2.le.cn2p) then
                 dc_l=cnpar/dsqrt(cn2p)
                 gam_l=dacos(dc_l)
                 write(*,*)'in cninit12_n_gam cn2p,gam_l',cn2p,gam_l
                else
                 write(*,*)'in cninit12_n_gam cnpar2>cn2p',cnpar2,cn2p
                endif
            endif  

            if (cn2m.gt.0.d0)then
              if(cnpar2.le.cn2m) then
                dc_l=cnpar/dsqrt(cn2m)
                gam_l=dacos(dc_l)
                write(*,*)'in cninit12_n_gam cn2m,gam_l',cn2m,gam_l
              else
                write(*,*)'in cninit12_n_gam cnpar2>cn2m',cnpar2,cn2m
              endif
            endif


            if((cn2m.lt.0.d0).and.(cn2p.lt.0d0)) then
             write(*,*)'in cninit2 two roots of the dispersion < 0'
             write(*,*)'cn2m,cn2p',cn2m,cn2p
             write(*,*)'the given wave can not exist in plasma'
	     iraystop=1
	     return
	    endif

20	    iroot=iroot+1
            write(*,*)'in cninit12 iroot,cntang2',iroot,cntang2

	    if(iroot.eq.1) then
              cn2=cn2p
	      if(cn2.lt.cntang2)then
	        write(*,*)'in cninit2 cn2p.lt.cntang2'
	        go to 20
	      end if
	    else
              cn2=cn2m
	      if(cn2.lt.cntang2)then
	        write(*,*)'in cninit2 cn2m.lt.cntang2'
	        write(*,*)'the given wave can not exist in plasma'
	        iraystop=1
	        return
	      end if
	    end if
            write(*,*)'in cninit cn2=',cn2
	    cnrho2=cn2-cnteta**2-cnphi**2
	    cnrho=dsqrt(cnrho2)
            cnperp=cnrho
cSm060906
            cnperp=dsqrt(cn2-cnpar2)
            write(*,*)'cnpar,dsqrt(cn2-cnpar2)',cnpar,dsqrt(cn2-cnpar2)
c-------------------------------------------------------------------
          write(*,*)'cninit12 cn2',cn2
          call cnzcnr(z,r,phi,cnteta,cnphi,cnrho,cnz,cnr,cm)

          write(*,*)'cninit12  after cnzcnr z,r,phi', z,r,phi
          write(*,*)'cnteta,cnphi,cnrho',cnteta,cnphi,cnrho
          write(*,*)'cnz,cnr,cm',cnz,cnr,cm
c-------------------------------------------------------------------
          cm=cnphi*r
	  gam=gamma1(z,r,phi,cnz,cnr,cm)
          ds=dsin(gam)
          dc=dcos(gam)
          ds2=ds*ds
          dc2=dc*dc
          ds4=ds2*ds2
c---------------------------------------------------------------------
c         controle that cn2 and gam are the solution of the dispersion
c         relation  n**2=(-b+ioxm*sqrt(b**2-4*a*c))/(2*a)
c---------------------------------------------------------------------
          call abc(z,r,phi,ds2,dc2,ad,bd,cd)

          d4=ad
	  d2=bd
	  d0=cd
          det=d2*d2-4.d0*d4*d0
          write(*,*)'gam,d4,d2,d0,det',gam,d4,d2,d0,det
c-------------------------------------------------------------
          cn1=dsqrt(d2*d2-4.d0*d4*d0)
          write(*,*)'(-d2+cn1)/(2.d0*d4)',(-d2+cn1)/(2.d0*d4)
          write(*,*)'(-d2-cn1)/(2.d0*d4)',(-d2-cn1)/(2.d0*d4)
          cn2new=(-d2+ioxm*cn1)/(2.d0*d4)
	  write(*,*)'cninit12 cn2new=',cn2new,'iroot=',iroot
	  if(iroot.eq.1) then
	      dnp=dabs(cn2-cn2new)
c              write(*,*)'cninit12 dnp,epsmode',dnp,epsmode
	      if (dnp.gt.epsmode)then
	         goto 20
	      end if
	  else
	      dnm=dabs(cn2-cn2new)
	      if (dnm.gt.epsmode)then
                write(*,*)'the given mode can not exist in plasma'
	        iraystop=1
	        return
	      end if
	  end if
          goto 111
        end if
c end if 2
c
c       ib.gt.1 iones resonance condition may be
c  if 3
        if (ib.gt.1) then
c         write(*,*)'in cninit: cold plasma ib .gt.1 x,r.phi,ib',
c     +   z,r,phi,ib
          xb=x(z,r,phi,ib)
          yb=y(z,r,phi,ib)
          xe=x(z,r,phi,1)
          ye=y(z,r,phi,1)

c          write(*,*)'in cninit12_n_gam z,r,ib,xb,yb',z,r,ib,xb,yb
c          write(*,*)'in cninit12_n_gam cnpar',cnpar
c          write(*,*)'in cninit12_n_gam s1,s2,s3,s4',s1,s2,s3,s4

	  pype=xe/(1.d0+ye)
	  pyme=xe/(1.d0-ye)
	  pyme2=pype/(1.d0-ye)
	  pypb=xb/(1.d0+yb)
	  delib=1.d0-yb
	  f1b=(s1-pyme2)
	  f0b=-pypb

	  g1b=cnpar2*(s4-(s1-pyme2))-(s2-pype)*(s3-pyme)-s4*(s1-pyme2)
	  g0b=cnpar2*pypb+xb*(s3-pyme)+pypb*s4

	  w1b=cnpar2*(-s4*(s1-pyme2)+(s2-pype)*(s3-pyme))+
     1            s4*(s2-pype)*(s3-pyme)
          w0b=cnpar2*(pypb*s4-xb*(s3-pyme))-s4*(s3-pyme)*xb
	  fd=f1b*delib+f0b
	  gd=g1b*delib+g0b
	  wd=w1b*delib+w0b
	  detin=gd**2-4.d0*fd*wd

c       write(*,*)'in cninit12_n_gam f1b,delib,f0b,fd',f1b,delib,f0b,fd
c       write(*,*)'in cninit12_n_gam g1b,delib,g0b,gd',g1b,delib,g0b,gd
c       write(*,*)'in cninit12_n_gam w1b,delib,w0b,wd',w1b,delib,w0b,wd

c          write(*,*)'in cninit12_n_gam gd,fd,wd,detin',gd,fd,wd,detin

	  if (detin.lt.0d0) then
	      write(*,*)' 3 in cninit detin  less then zero '
	      iraystop=1
	      return
	  end if
c	  cn2=(-gd+ioxm*dsqrt(detin))/(2.*fd)
c          write(*,*)'cninit: gd,fd,wd,detin',gd,fd,wd,detin
	  cn2p=(-gd+dsqrt(detin))/(2.d0*fd)
	  cn2m=(-gd-dsqrt(detin))/(2.d0*fd)
         
          if((cn2m.lt.0.d0).and.(cn2p.lt.0.d0)) then
            write(*,*)'in cninit2 two roots of the dispersion < 0'
            write(*,*)'cn2m,cn2p',cn2m,cn2p
            write(*,*)'the given wave can not exist in plasma'
	    iraystop=1
	    return
	   endif

	  write(*,*)'in cninit ib.qt.1 cn2p,cn2m',cn2p,cn2m
c	  write(*,*)'in cninit fd,gd,wd',fd,gd,wd
30	  iroot=iroot+1

          write(*,*)'in cninit12 iroot,cntang2',iroot,cntang2

	  if(iroot.eq.1) then
             cn2=cn2p
	     if(cn2.lt.cntang2)then
	        write(*,*)'in cninit2 cn2p.lt.cntang2'
	        go to 30
	     end if
	  else
             cn2=cn2m
	     if(cn2.lt.cntang2)then
	        write(*,*)'in cninit2 cn2m.lt.cntang2'
	        write(*,*)'the given wave can not exist in plasma'
	        iraystop=1
	        return
	     end if
	  end if

          write(*,*)'in cninit cn2=',cn2
	  cnrho2=cn2-cnteta**2-cnphi**2
	  cnrho=dsqrt(cnrho2)   
          cnperp=cnrho
          write(*,*)'cnperp=cnrho',cnperp
cSm050826
          cnperp=dsqrt(cn2-cnpar2)
          write(*,*)'cnpar,dsqrt(cn2-cnpar2)',cnpar,dsqrt(cn2-cnpar2)
c------------------------------------------------------------------
c          write(*,*)'cnini12 ib>1 before call cnzcnr'
c          write(*,*)'z,r,phi,cnteta,cnphi,cnrho',
c     &    z,r,phi,cnteta,cnphi,cnrho 
          call cnzcnr(z,r,phi,cnteta,cnphi,cnrho,cnz,cnr,cm)
c          write(*,*)'cninit12 after cnzcnr cnz,cnr,cm',cnz,cnr,cm

c---------test begin
cSm050906
          cnpar_test=(bz*cnz+br*cnr+bphi*cm/r)/bmod
          write(*,*)'cninit12 cnpar,cnpar_test', cnpar,cnpar_test
          write(*,*)'cn2,cnz**2+cnr**2+(cm/r)**2',
     &    cn2,cnz**2+cnr**2+(cm/r)**2
c---------test_end

c---------------------------------------------------------------------
          cm=cnphi*r
	  gam=gamma1(z,r,phi,cnz,cnr,cm)
          write(*,*)'cninit12 from gamma1 gam',gam 

cSm050825
c---------test begin
c          gam=dacos(cnpar/dsqrt(cn2))
c---------test_end

          ds=dsin(gam)
          dc=dcos(gam)
          ds2=ds*ds
          dc2=dc*dc
          ds4=ds2*ds2
c--------------------------------------------------------------------
c test hamilt
          write(*,*)'in cninit cnpar2,cn2*dc2',cnpar2,cn2*dc2
c--------------------------------------------------------------------
c         controle that cn2 and gam are the solution of the dispersion
c         relation  n**2=(-b+ioxm*sqrt(b**2-4*a*c))/(2*a)
c---------------------------------------------------------------------
          call abc(z,r,phi,ds2,dc2,ad,bd,cd)
          d4=ad
	  d2=bd
	  d0=cd
          det=d2*d2-4.d0*d4*d0
          write(*,*)'gam,d4,d2,d0,det',gam,d4,d2,d0,det
c---------------------------------------------------
c test hamilt=0?
          hamtest=fd*cn2*cn2+gd*cn2+wd

         write(*,*)'in cninit fd,gd,wd,cn2*cn2,cn2,hamtest'
         write(*,*) fd,gd,wd,cn2*cn2,cn2,hamtest

          hamtest=d4*cn2*cn2+d2*cn2+d0
          pt4=d4*cn2*cn2
          pt2=d2*cn2
          pt=pt4+pt2+d0
          ptm=pt4-pt2+d0
c          write(*,*)'cninit12 pt4,pt2,d0',pt4,pt2,d0
c          write(*,*)'pt,ptm',pt,ptm
          write(*,*)'in cninit d4,d2,d0,cn4,cn2,hamtest'
          write(*,*)d4,d2,d0,cn2*cn2,cn2,hamtest
c         write(*,*)'(fd-d4),(gd-d2),(wd-d0)'
c         write(*,*)(fd-d4),(gd-d2),(wd-d0)
          gdt=dc2*cn2*((1.d0-xe-xb)
     1       -(1.d0-xb/(1.d0-yb*yb)-xe/(1.d0-ye*ye)))
     1       +(-(1.d0-xb/(1.d0-yb)-xe/(1.d0+ye))*
     1       (1.d0-xb/(1.d0+yb)-xe/(1.d0-ye))-
     1       (1.d0-xb/(1.d0-yb*yb)-xe/(1.d0-xe*xe))*(1.d0-xe-xb))
          gdt=gdt*delib
c         write(*,*)'gd,gdt',gd,gdt
          btest=-(1.d0-xb/(1.d0-yb*yb)-
     1         xe/(1.d0-ye*ye))*(1.d0-xe-xb)*
     1         (1.d0+dc2)-(1.d0-xb/(1.d0-yb)-xe/(1.d0+ye))*
     1         (1.d0-xb/(1.d0+yb)-xe/(1.d0-ye))*ds2
c         write(*,*)'bd,btest',bd,btest
          ptt=(xe*ye*ye*delib/(1.d0-ye*ye)+xb*yb*yb/(1.d0+yb))
          fmina=-dc2*ptt
c         write(*,*)'xe,ye,xb,yb',xe,ye,xb,yb
          ptt1=(-(1.d0-xb/(1.d0-yb)-xe/(1.d0+ye))*
     1         (1.d0-xb/(1.d0+yb)-xe/(1.d0-ye))+
     1    (1.d0-xe-xb)*(1.d0-xe/(1.d0-ye*ye)-xb/(1.d0-yb*yb)))*delib
          gminb=dc2*(-cn2*ptt-ptt1)
          wminc=-dc2*cn2*ptt1
c         write(*,*)'fmina,gminb,wminc'
c         write(*,*)fmina,gminb,wminc
c-------------------------------------------------------------
          cn1=dsqrt(d2*d2-4.d0*d4*d0)
          write(*,*)'(-d2+cn1)/(2.d0*d4)',(-d2+cn1)/(2.d0*d4)
          write(*,*)'(-d2-cn1)/(2.d0*d4)',(-d2-cn1)/(2.d0*d4)
          cn2new=(-d2+ioxm*cn1)/(2.d0*d4)
	  write(*,*)'cn2new=',cn2new
	  if(iroot.eq.1) then
	     dnp=dabs(cn2-cn2new)
	     if (dnp.gt.epsmode)then
	        goto 30
	     end if
	  else
	     dnm=dabs(cn2-cn2new)
	     if (dnm.gt.epsmode)then
	        write(*,*)' the given mode cannot exist in plasma '
	        iraystop=1
	        return
	     end if
	  end if
          goto 111
        end if
c end if 3
      end if
c end if 0
 
    
  111 continue
                                                              
      return
      end



      subroutine n_cold_gam(z,r,phi,gam,cn_p,cn_m,iraystop_p,iraystop_m)
c----------------------------------------------------------------
c     Calculates cold plasma dispersion relalation roots
c     cn=N(gam,ioxm) in the given space point (z,r,phi).
c     as a root of the equation a*N**4+b**2+c=0,
c     cn=(-b+ioxm*sqrt(b**2-4ac))/2a
c     cn_p=(-b+sqrt(b**2-4ac))/2a for iocm=+1
c     cn_m=(-b-sqrt(b**2-4ac))/2a for ioxm=-1
c
c     Here: 
c     coefficients a=A*delta**3,b=B*delta**3,c=C*delta**3,
c
c     See book  Kroll,Travelspils
c     A=eps_1*sin(gam)**2+eps_3*cos(gam)**2
c     B=-eps_1*eps_2(1+cos(gam)**2)-(eps_1**2-eps_2**2)sin(gam)**2
c     C=eps_3*(eps_1**2-eps_2**2)
c
c     delta=1-y_e for electrons ib=1 (ib is used in subroutine
c                                     abc) 
c          =1-y_i for ions      ib=i>1
c----------------------------------------------------------------     
      implicit none
c-----input
      real*8 z,r,phi, ! space coordinates of the given point
     *gam             ! the angle [radians] between the wave vector
                      ! and the magnetic field                      
c-----output
      real*8 cn_p,cn_m    !roots cn_p=N(gam,ioxm=+1), cn_m=N(gam,ioxm=-1),  
      integer iraystop_p, !=0 the root with ioxm=+1 found
                          !=1 no root  with ioxm=+1 did not find
     &        iraystop_m  !=0 the root with ioxm=-1 found
                          !=1 no root  with ioxm=-1 did not find
 
c-----locals
      real*8 ds,dc,ds2,dc2,
     &a,b,c,det,cn2_p,cn2_m

      iraystop_p=0
      iraystop_m=0
      
      ds=dsin(gam)
      dc=dcos(gam)
      ds2=ds*ds
      dc2=dc*dc

      call abc(z,r,phi,ds2,dc2,a,b,c)

      det=b**2-4.d0*a*c
      if(det.lt.0.d0)then
        write(*,*)'in subroutine n_cold_gam det<0'
        write(*,*)'the cold plasma dispersion eq.'
        write(*,*)'has not root N_perp(gam)'
        write(*,*)'for the given angle gam=',gam
        iraystop_p=1
        iraystop_m=1
        return
      endif


      cn2_p=(-b+dsqrt(det))/(2.d0*a) !N**2
      if(cn2_p.gt.0.d0)then
        cn_p=dsqrt(cn2_p) 
      else
        write(*,*)'in subroutine n_gam negative N**2=cn2_p<0'
        write(*,*)'The cold plasma dispersion eq.'
        write(*,*)'has not positive root N_perp(gam,ioxm)'
        write(*,*)'for the given angle gam=',gam,' and ioxm=+1'
        iraystop_p=1
      endif

      cn2_m=(-b-dsqrt(det))/(2.d0*a) !N**2
      if(cn2_m.gt.0.d0)then
        cn_m=dsqrt(cn2_m) 
      else
        write(*,*)'in subroutine n_cold_gam negative N**2=cn2_m<0'
        write(*,*)'The cold plasma dispersion eq.'
        write(*,*)'has not positive root N_perp(gam,ioxm)'
        write(*,*)'for the given angle gam=',gam,' and ioxm=-1'
        iraystop_p=1
      endif

      return
      end

      

c        **********************nper_npar_ioxm_n_npar**********
c                                -                            
c         It solves the cold plasma id=1,2 or
c         Appleton-Harty id=3 dispersion relation N_per=N_per(n_par)
c         for given ioxm_n_nper                              
c         Then subroutine calculates the initial components  
c         of the refractive index  cnz,cnr,cm                
c        ******************************************************
c
c---------------------------------------------------
c     							   !
c     input parameters					   !
c     z,r,phi,cnpar 	                                   !
c     ioxm_n_npar is inside one.i common block             !
c     It works for cold plasma dispesion functions         !
c     id_loc=1,2,3                                         !
c     output parameters:N_per(N_par)=cnper                 !
c     iraystop=o the root was found                        !
c     iraystop=1 root was not found                        !
c----------------------------------------------------------
c     it uses the following functions and subroutines      !
c     ias1r,b,y,x,gamma1,s,abc                             !
c-----------------------------------------------------------

      subroutine nper_npar_ioxm_n_npar(id_loc,z,r,phi,cnpar,
     &cnper,iraystop)

      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'

      iraystop=0

      pi=4*datan(1.d0)
      cnpar2=cnpar*cnpar
      cnpar4=cnpar2*cnpar2
      bmod=b(z,r,phi)
c-------------------------------------------------------------
c     calculations of  cnrer2p and cnper2m
c     from cnpar2 by using the dispersin relation
c     f*cnper**4 +(2f*cnpar**2+g)*cnper**2+(f*cnpar**4+g*cnpar**2+w)=0
c     f*cnper**4 +gnew*cnper**2+wnew=0
c     gnew=2fcnpar**2+g, wnew=f*cnpar**4+g*cnpar**2+w
c     in the following form cnper**2=(-gnew+,-*sqrt(gnew**2-4*f*wnew))/(2*f)
c------------------------------------------------------------------
c     Appleton - Hartry dispersion relation
      if (id_loc.eq.3) then
         xi=x(z,r,phi,1)
         yi=y(z,r,phi,1)
         py2=yi*yi
         py4=py2*py2
         px=1.d0-xi
         px2=px*px
c------------------------------------------------------------------
         pyp=xi/(1.d0+yi)
         f1e=1.d0
         f0e=-pyp
         g1e=cnpar2*(-xi)+pyp+xi-2.d0
         g0e=cnpar2*pyp+xi*(1.d0-pyp)+pyp*(1.d0-xi)
         w1e=cnpar2*(xi-pyp)+(1.d0-pyp)*(1.d0-xi)
         w0e=cnpar2*(pyp*(1.d0-xi)-xi*(1.d0-pyp))-
     &       (1.d0-xi)*(1.d0-pyp)*xi
         dele=1.-yi
         fd=f1e*dele+f0e
         gd=g1e*dele+g0e
         wd=w1e*dele+w0e
c new coefficients
         gnew=gd+2.d0*fd*cnpar2
         wnew=wd+gd*cnpar2+fd*cnpar4
         gd=gnew
         wd=wnew

         detin=gd**2-4.d0*fd*wd
         if(detin.lt.0.d0) then
           write(*,*)'1 nper_npar_ioxm_n_nparn detin less then zero'
           write(*,*)'nper_npar_ioxm_n_npar no roots'
           iraystop=1
           return
         endif

         cnper2=(-gd+ioxm_n_npar*dsqrt(detin))/(2.d0*fd)
         if(cnper2.lt.0d0) then
            write(*,*)'2 nper_npar_ioxm_n_nparn cnper2<0'
            write(*,*)'the root N_perp**2 is negative'
            write(*,*)'nper_npar_ioxm_n_npar no positive root'    
            iraystop=1
            return
         else
            cnper=dsqrt(cnper2)
         endif           
      end if !id_loc.eq.3
c------------------------------------------------------------------
c     cold plasma dispresion relation
c------------------------------------------------------------------
      if ((id_loc.eq.1).or.(id_loc.eq.2)) then
        call s(z,r,phi,s1,s2,s3,s4,s6,s7)

        if (ib.eq.1) then
c---------ib=1 electron resonance condition may be
c         write(*,*)'nper_npar_ioxm_n_par cold plasma ib=1 '
          xe=x(z,r,phi,1)
          ye=y(z,r,phi,1)
          pyp=xe/(1.d0+ye)
          dele=1.d0-ye
          f1e=s7
          f0e=-pyp
          g1e=cnpar2*(s4-s7)-(s6-pyp)*s3-s4*s7
          g0e=cnpar2*pyp+xe*(s6-pyp)+pyp*s4
          w1e=cnpar2*(-s7*s4+(s6-pyp)*s3)+(s6-pyp)*s3*s4
          w0e=cnpar2*(pyp*s4-xe*(s6-pyp))-s4*(s6-pyp)*xe
          fd=f1e*dele+f0e
          gd=g1e*dele+g0e
          wd=w1e*dele+w0e
c new coefficients
          gnew=gd+2.d0*fd*cnpar2
          wnew=wd+gd*cnpar2+fd*cnpar4
          gd=gnew
          wd=wnew

          detin=gd**2-4.d0*fd*wd
          if(detin.lt.0.d0) then           
            write(*,*)'2 nper_npar_ioxm_n_npar detin less then zero'
            write(*,*)'nper_npar_ioxm_n_npar no roots'
            iraystop=1
            return
	  end if

          
          cnper2=(-gd+ioxm_n_npar*dsqrt(detin))/(2.d0*fd)
          if(cnper2.lt.0.d0)then
            write(*,*)'2 nper_npar_ioxm_n_nparn cnper2<0'
            write(*,*)'the root N_perp**2 is negative'
            write(*,*)'nper_npar_ioxm_n_npar no positive root'    
            iraystop=1
            return
          else
            cnper=dsqrt(cnper2)
          endif           
          goto 111
        end if !ib.eq.1


        if(ib.gt.1) then
c---------ib.gt.1 iones resonance condition may be
c         write(*,*)'cold plasma ib .gt.1  '
          xb=x(z,r,phi,ib)
          yb=y(z,r,phi,ib)
          xe=x(z,r,phi,1)
          ye=y(z,r,phi,1)
          pype=xe/(1.d0+ye)
          pyme=xe/(1.d0-ye)
          pyme2=pype/(1.d0-ye)
          pypb=xb/(1.d0+yb)
          delib=1.d0-yb
	  f1b=(s1-pyme2)
	  f0b=-pypb

	  g1b=cnpar2*(s4-(s1-pyme2))-(s2-pype)*(s3-pyme)-s4*(s1-pyme2)
          g0b=cnpar2*pypb+xb*(s3-pyme)+pypb*s4

	  w1b=cnpar2*(-s4*(s1-pyme2)+(s2-pype)*(s3-pyme))+
     1	     s4*(s2-pype)*(s3-pyme)
          w0b=cnpar2*(pypb*s4-xb*(s3-pyme))-s4*(s3-pyme)*xb
          fd=f1b*delib+f0b
          gd=g1b*delib+g0b
          wd=w1b*delib+w0b
c new coefficients
          gnew=gd+2.d0*fd*cnpar2
          wnew=wd+gd*cnpar2+fd*cnpar4
          gd=gnew
          wd=wnew

          detin=gd**2-4.d0*fd*wd
          if(detin.lt.0.d0) then           
            write(*,*)'2 nper_npar_ioxm_n_npar detin less then zero'
            write(*,*)'nper_npar_ioxm_n_npar no roots'
            iraystop=1
            return
	  end if

          cnper2=(-gd+ioxm_n_npar*dsqrt(detin))/(2.d0*fd)
          if(cnper2.lt.0.d0) then
            write(*,*)'3 nper_npar_ioxm_n_nparn cnper2<0'
            write(*,*)'the root N_perp**2 is negative'
            write(*,*)'nper_npar_ioxm_n_npar no positive root'    
            iraystop=1
            return
          else
            cnper=dsqrt(cnper2)
          endif           
          goto 111
        end if !ib>1
      end if !(id_loc.eq.1).or.(id_loc.eq.2)
 
  111 continue

      return
      end


c        **********************cninit12***************************
c        *                        -                               *
c        * It solves the dispersion relation N=N(n_par)           *
c        * for cold plasma                                        *
c        * Then subroutine calculates the initial components      *
c        * of the refractive index  cnz,cnr,cm         	          *
c        * If ioxm_n_npar=0 it uses root N(gam,ioxm)              *
c        * If ioxm_n_npar=1 or -1 it uses root N(npar ioxm_n_npar)*
c        *                then calculates angle gam(N_par,N_per)  *
c        *                then finds ioxm which gives the root    *
c        *                N(gam,ioxm)=N(N_par,ioxm_n_npar)        *
c        
c        *********************************************************
c
c------------------------------------------------------------------
c					                          !
c     input parameters		         		          !
c     z,r,phi,cnpar,cnteta,cnphi 		       	          !
c                                                                 !
c     output parameters				                  !
c     cnz,cnr,cm are the components of the refractive index       !        
c     iraystop=1 end ray calculation                              !
c------------------------------------------------------------------
c     it uses the following functions                             !
c     b,gamma1                                                    !
c------------------------------------------------------------------
      subroutine cninit12(z,r,phi,cnpar,cnteta,cnphi,
     1                  cnz,cnr,cm,iraystop)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
cSAP090504
      include 'grill.i'

      save
      double complex cmplnper      
      
      iraystop=0
      pi=4*datan(1.d0)
      cnpar2=cnpar*cnpar
      cntang2=cnteta*cnteta+cnphi*cnphi
      bmod=b(z,r,phi)
      write(*,*)'cninit12 ioxm_n_npar=',ioxm_n_npar
      if(ioxm_n_npar.eq.0) then
c-------------------------------------------------------------
c       calculates the root N(Npar)=N(gam,ioxm) for givem ioxm
c-------------------------------------------------------------
        call cninit12_n_gam(z,r,phi,cnpar,cnteta,cnphi,
     &                      cnz,cnr,cm,iraystop)
      else
c-------------------------------------------------------------------
c       calculates the root N(N_par,ioxm_n_npar) for given ioxm_n_npar
c       finds ioxm to get  N(N_par,ioxm_n_npar)=N(gam,ioxm)
c-------------------------------------------------------------------
        write(*,*)'z,r,phi,cnpar,cnteta,cnphi',
     &             z,r,phi,cnpar,cnteta,cnphi
c-----------------------------------------------------------------
        iroot=0
c------------------------------------------------------------------
c       epsmode is the accuracy  for selection mode (O or X) on
c             the plasma boundary
c------------------------------------------------------------------
        epsmode=1.d-8

        id_loc=2
c------------------------------------------------------------------
c       solve the cold plasma id=1,2  
c       dispersion relation cnpar=N_per(n_par)
c------------------------------------------------------------------        
        call nper_npar_ioxm_n_npar(id_loc,z,r,phi,cnpar,
     &  cnper,iraystop) !ioxm_n_npar was set in one.i
        write(*,*)'cninit.f in cninit12 after nper_npar_ioxm_n_npar'
        write(*,*)'ioxm_n_npar,cnper,iraystop',
     &             ioxm_n_npar,cnper,iraystop
c-------------------------------------------------------------------
        if(iraystop.eq.1) then
          write(*,*)'in cninit12 no root N_(npar,ioxm_n_npar)'
          return
        else 
cSAP090504  
           if (i_n_poloidal.eq.3) then !input N_parallel, ksi_nperp
c-------------calculate N_phi,N_theta,N_rho
              gradpsi=dsqrt(dpdrd*dpdrd+dpdzd*dpdzd)
              b_teta=(bz*dpdrd-br*dpdzd)/gradpsi !poloidal magnetic field
              write(*,*)'cninit.f i_n_poloidal=3 ksi_nperp',ksi_nperp
              rad_ksi_nperp=ksi_nperp*pi/180.d0 !transfrm degrees to radians
              write(*,*)'rad_ksi_nperp',rad_ksi_nperp
              cnteta=(cnpar*b_teta+cnper*bphi*dsin(rad_ksi_nperp))
     &                    /bmod
              cnphi=(cnpar*bphi-cnper*b_teta*dsin(rad_ksi_nperp))
     &                    /bmod   
              cnrho=cnper*dcos(rad_ksi_nperp)
              write(*,*)'cninit.f i_n_poloidal=3 cnteta,cnphi',
     &                                            cnteta,cnphi
              write(*,*)'cninit.f i_n_poloidal=3 cnrho_full',cnrho_full
              cntang2=cnteta*cnteta+cnphi*cnphi
           endif ! i_n_pol=3

          
          cn2_npar=cnper**2+cnpar**2
          if(cn2_npar.lt.cntang2)then
            write(*,*)'in cninit12 (cn2_npar.lt.cntang2)'
            write(*,*)'no root N(npar,ioxm_n_npar) > N_tang'
            iraystop=1
            return
          else
            cnrho2=cn2_npar-cnteta**2-cnphi**2
	    cnrho=dsqrt(cnrho2)
c--------------------------------------------------------------
c           calculate refractive index components cnz,cnr,cm
c           for given z,r,phi,cntheta,cnphi,cnrho
c-------------------------------------------------------------- 
            write(*,*)'cninit.f in cninit12 before cnzcnr'
            write(*,*)'r,phi,cnteta,cnphi,cnrho',
     &                 r,phi,cnteta,cnphi,cnrho

cSAP090601
c           In this case ioxm_n_npar.ne.0 and
c           ioxm is not determined
c           Subroutine cnzcnr uses subroutine rside1, 
c           which for id=2 case needs the determined value ioxm 
c           To avoid this problem we will use id=1 in cnzcnr
c           May be signs of the radial group velocity
c           for id=1 and id=2 can be different?
c           If yes, then it can creat the new problem. 

            if (id.eq.2) then
              id_loc_2=id
              id=1
            endif

            call cnzcnr(z,r,phi,cnteta,cnphi,cnrho,cnz,cnr,cm)
cSAP090601          
            id=id_loc_2 

            write(*,*)'cninit.f in cninit12 after cnzcnr'
            write(*,*)'cnz,cnr,cm',cnz,cnr,cm

c--------------------------------------------------------------
c           calculate the angle gam between the refractive index 
c           and the magnetic field
c-------------------------------------------------------------
            gam=gamma1(z,r,phi,cnz,cnr,cm)
c-------------------------------------------------------------
c           calculate cold plasma roots roots 
c           cn_p=N(gam,ioxm=+1), cn_m=N(gam,ioxm=-1)
c-------------------------------------------------------------          
            call n_cold_gam(z,r,phi,gam,cn_p,cn_m,
     &                      iraystop_p,iraystop_m)
c-------------------------------------------------------------
c           choose ioxm for which
c           N(N_par,ioxm_n_npar)=N(gam,ioxm)
c-------------------------------------------------------------
c           check ioxm=+1 root
c-------------------------------------------------------------
            if(iraystop_p.eq.0) then
              delta=dabs(cn2_npar-cn_p**2)
              if(delta.gt.epsmode)then
c---------------root N(N_par,ioxm_n_npar).ne.N(gam,ioxm=+1)
              else
c---------------root N(N_par,ioxm_n_npar)=N(gam,ioxm=+1)
                ioxm=+1
                iroot=iroot+1
                write(*,*)'in cninit12 found ioxm=',ioxm 
c               goto 10 
              endif                
	    else
c-------------no root for N(gam,ioxm=+1)
            endif !iraystop_p.eq.0

c-------------------------------------------------------------
c           check ioxm=-1 root
c-------------------------------------------------------------
            if(iraystop_m.eq.0) then
              delta=dabs(cn2_npar-cn_m**2)
              if(delta.gt.epsmode)then
c---------------root N(N_par,ioxm_n_npar).ne.N(gam,ioxm=-1)
              else
c---------------root N(N_par,ioxm_n_npar)=N(gam,ioxm=-1)
                ioxm=-1                  
                write(*,*)'in cninit12 found ioxm=',ioxm 
                iroot=iroot+1
c               goto 10 
              endif                
	    else
c-------------no root for N(gam,ioxm=-1)
            endif !iraystop_m.eq.0

c------------------------------------------------------------
            if(iroot.eq.0) then
              write(*,*)'in cninit12 no ioxm was found'
              write(*,*)'to get N(N_par,ioxm_n_par)=N(gam,ioxm)'
              iraystop=1
              return 
            else
              if(iroot.eq.2) then
                write(*,*)'*******WARNING******************'
                write(*,*)'in cninit12 two ioxm was found'
                write(*,*)'to get N(N_par,ioxm_n_par)=N(gam,ioxm)'
              endif !iroot.eq.2
            endif !iroot.eq.0
          endif !cn2_npar.lt.cntang2
        endif !iraystop.eq.1
      endif !ioxm_n_npar=0

c  10 continue
                                                              
      return
      end

