c        ********************** abc   ***********************
c        * this subroutine calculates the coefficients      *
c        * a,b,c, for the cold dispersion relation          *
c        * d=a*n**4+b*n**2+c                                *
c        ****************************************************
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
c      	input parameters					    !	
c	point coordinates: z,r,phi    				    !
c       ds2=sin(gam)**2,dc2=cos(gam)**2                             !
c       gam is the angle between magnetic field and refractiv index !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
      subroutine abc(z,r,phi,ds2,dc2,ad,bd,cd)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
c  ---------------------------------------------------------------
c  cold plasma dispersion relation with electrons and ions
c  x(i=1),y(i=1) for electron component
c  x(i may be=2,nbulk),y(i may be=2,...,nbulk) ions components.
c  ib is a number of component for which delib=1-yib may be equal
c  zero inside the plasma,ib may be=from 1 to nbulk
c  The dispersion relation is multiplied by delib
c  --------------------------------------------------------------
      pc=1.d0+dc2
      call s(z,r,phi,s1,s2,s3,s4,s6,s7)
      xib=x(z,r,phi,ib)
      yib=y(z,r,phi,ib)
      delib=1.d0-yib

c----------------------------------------------------------------
c  ib =1 (the cyclotron resonance conditions dele=0 may be
c         realised in plasma, dispersion relation is
c         multiplied by dele )
c         ad,bd,cd calculations
c------------------------------------------------------------------
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
      end if
c-----------------------------------------------------------------
c  ib .gt.1 (the cyclotron resonance conditions delib=0 may be
c            realised in plasma, dispersion relation is
c            multiplied by delib )
c            ad,bd,cd calculations
c------------------------------------------------------------------
      if (ib.gt.1) then
         xe=x(z,r,phi,1)
         ye=y(z,r,phi,1)
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
      end if

      return
      end



