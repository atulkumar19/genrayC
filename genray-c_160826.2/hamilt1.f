
c        ********************** hamilt1 **********************
c        *                      ------                      *
c        * this function calculates the hamiltonian of the  *
c        * the system of geometrical optics equations       *
c        * It will calculate the dielectric tensor
c        * reps(3,3) and will put this tensor to eps.i file
c        ****************************************************
c
c-------------------------------------------------------------------
c								    !
c        input parameters					    !
c      								    !
c      z, r phi - coordinates of the point where the hamiltonian is !
c                 calculated.      		         	    !
c								    !
c      cnz, cnr, cm - n_z, n_r, r*n_phi - components of  wave  ref- !
c                     ractive index at this point.                  !
c      the angle gam between refractive index n and magnetic field  !
c      from common 'one'					    !
c-------------------------------------------------------------------
      double precision
     1function hamilt1(z,r,phi,cnz,cnr,cm)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'eps.i'
      include 'ions.i'
      double precision dshkarof
      double complex dhot,dhot_sum
      double complex ceps(3,3),hamiltc
      external x,y,tempe,dshkarof
      external dhot_sum

      double precision x_ar(nbulka),y_ar(nbulka),t_av_ar(nbulka),
     .tpop_ar(nbulka),vflow_ar(nbulka)
      
      double complex K(3,3),dK_dx_ar(3,3,nbulka),dK_dy_ar(3,3,nbulka),
     &dK_dt_ar(3,3,nbulka),dK_dnper(3,3),dK_dnpar(3,3)


      double complex compl_nper !for Eric tensor
      double complex eps_weiss(3,3)
c-----external
      double complex det       

      bmod=b(z,r,phi)
      r2=r*r
      cn2=(cnz*cnz+cnr*cnr+cm*cm/r2)
      cnt=dsqrt(cn2)
      cn4=cn2*cn2
           
      ds=dsin(gam)
      dc=dcos(gam)
         cnpar=cnt*dc
         cnper=cnt*ds
      ds2=ds*ds
      dc2=dc*dc
      if (((id.eq.1).or.(id.eq.2)).or.(id.eq.3)) then
        call tensrcld(z,r,phi)
      end if
c---------------------------------------------------------
c
c     Appleton-Hartree dispersion relation
c---------------------------------------------------------
      if (id.eq.3) then
         ds4=ds2*ds2
         xi=x(z,r,phi,1)
         yi=y(z,r,phi,1)
         py2=yi*yi
         py4=py2*py2
         px=1.d0-xi
         px2=px*px
         sqrdet=dsqrt(py4*ds4+4.*py2*px2*dc2)
         pz=2.d0*px-py2*ds2+ioxm*sqrdet
         hamilt1=cn2-(1.d0-2.d0*xi*px/pz)
         goto 10
      end if

c--------------------------------------------------------
      if ((id.eq.1).or.(id.eq.2)) then
c--------------------------------------------------------
c       id=1 or id=2 cold plasma dispersion relation
c--------------------------------------------------------
        call abc(z,r,phi,ds2,dc2,ad,bd,cd)
        d4=ad
       	d2=bd
       	d0=cd
c--------------------------------------------------------
        if (id.eq.1) then
           hamilt1=d4*cn4+d2*cn2+d0
        end if

        if (id.eq.2) then
           hamilt1=cn2+(d2-ioxm*dsqrt(d2*d2-4.d0*d4*d0))/
     *                (2.d00*d4)
        end if
        go to 10
      end if
c     end cold plasma dispersion relation
c-----------------------------------------------------------
c     det=dsqrt(d2*d2-4.d00*d4*d0)
c     cn2od=(-d2+det)/(2.d00*d4)
c     cn2ex=(-d2-det)/(2.d00*d4)
c     write(*,*)'cn2od=',cn2od,'cn2ex=',cn2ex
c     cn2c=cnz*cnz+cnr*cnr+cm*cm/(r*r)
c     write(*,*)'cn2c=',cn2c
c-----------------------------------------------------------
c-------------------------------------------------------------
c     Hot non-relativistic plasma
      if (id.eq.6) then

         do i=1,nbulk
           x_ar(i)=x(z,r,phi,i)
           y_ar(i)=y(z,r,phi,i)
           if(i.eq.1) y_ar(1)=-y_ar(1)
           te=tempe(z,r,phi,i) ! kev
           t_av_ar(i)=te*1000.d0      ! ev 
           tpop_ar(i)=tpoprho(rho,i)
           vflow_ar(i)=vflowrho(rho,i)
         enddo

         hamiltc=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     .   vflow_ar,cnpar,cnper,1,reps)

         hamilt1=dreal(hamiltc)
         go to 10
      end if
c **  end if forest
c-------------------------------------------------------------
c-----------------------------------------------------------
c     id=14 dispersion function
c     D=real part(D relativistic dispersion function)
c     from Eric Nelson_Melby dielectric tensor ,if npar.le. 0.38D0
c     from Ram Abhay dielectric tensor ,if npar.ge. 0.38D0
C 1Sep2005 -- Now Disp_combined is entirely Abhay's tensor, which 
C works just as well or faster than Nelson-Melby's, when the resolution
C is lower than the maximum like it used to be.
C
C ENM 15Mar2006 -- After finding that the combined version which jumps from
C     one to another dispersion relation depending on n_parallel didn't work
C     very well, and finding that the Ram version works well, even for fairly
C     small n_parallel, as long as you adjust the resolution parameters (see
C     genray.dat template), now id=14 is just the Ram tensor (id=11 is just the
C     Nelson-Melby tensor)
C
      if (id.eq.14) then
         x_ar(1)=x(z,r,phi,1)
         y_ar(1)=y(z,r,phi,1)    !question for electron -y? 
         te=tempe(z,r,phi,1) !kev

         compl_nper=dcmplx(cnper,0.d0)

         call Disp_Ram(te,cnpar,x_ar(1),y_ar(1),
     +   compl_nper,K,hamiltc) !K is in z-y plane
                               !in Stix coordinates

         if (iherm.eq.1) then
c-----------calcualte hermitian part reps of the complex tensor K
            call herm(K,reps)
            hamiltc=det(reps,cnpar,compl_nper)   
         endif
         hamilt1=dreal(hamiltc)
         go to 10
      end if
c **  end if id=14

 10   continue
      return
      end
