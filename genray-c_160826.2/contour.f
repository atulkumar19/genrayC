
      real function rbound(r8)
c
c     Converts a real*8 argument to a real number,
c     equal to 0. (if r8=0.) or,
c     bounded in absolute value by 1.e-33 and 1.e+33.
c     This can be used to convert real*8 numbers to
c     real numbers, and to keep the resulting numbers
c     within the specified bounds.  This is necessary
c     for the PGPLOT library running on a 32-bit machine.
c     (1.e-35 was found to be too small in some cases,
c      on the DEC Alpha).
c     For a 64-bit machine, one might consider appropriate
c     adjustment of em33/ep33.
c
      real*8 r8,r8sign,r8abs
      real*8 em33,ep33,zero,one
      data em33/1.d-33/, ep33/1.d+33/, zero/0.d0/, one/1.d0/

      r8abs=abs(r8)
      if (r8abs.ne.zero) then
         r8sign=sign(one,r8)
         r8abs=min(r8abs,ep33)
         rbound=r8sign*max(r8abs,em33)
      else
         rbound=0.
      endif
      return
      end


      subroutine aminmx(x,n1,n2,n,xmin,xmax)
      double precision x,xmin,xmax
      dimension x(n2)
      xmin=+1.d100
      xmax=-1.d100
      do 1  i=n1,n2,n
      xmin=dmin1(xmin,x(i))
 1    xmax=dmax1(xmax,x(i))
      return
      end


      subroutine rminmx(x,n1,n2,n,xmin,xmax)
      real x,xmin,xmax
      dimension x(n2)
cSm030515 these lines work for g77 under linux at PC
      xmin=+1.e30
      xmax=-1.e30
cSm030515 these lines work for digital fortran under windows at PC
c      xmin=+1.e38
c      xmax=-1.e38
      do 1  i=n1,n2,n
      xmin=amin1(xmin,x(i))
 1    xmax=amax1(xmax,x(i))
      return
      end






