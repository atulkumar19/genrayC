c        **********************   y   ***********************
c        *                        -                         *
c        * this function calculates the ratio of  cyclotron *
c        * frequency of species i and wave frequency        *
c        * y=omega_cyclotron_i/omega                        *
c        ****************************************************
c
c-------------------------------------------------------------------
c								    !
c        input parameters					    !
c      								    !
c      z, r phi - coordinates of the  point  where  the  ratio  is  !
c                 calculated.      		         	    !
c      small radius rho is in common one
c-------------------------------------------------------------------
      double precision
     1function y(z,r,phi,i)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      if(ixyz.eq.1) then
        stop 'y(z,r,phi,i): should not be called for ixyz=1'
      endif
      y=w(i)*bmod
      if(dabs(y).lt.1.d-8) y=1.d-8
      return
      end

