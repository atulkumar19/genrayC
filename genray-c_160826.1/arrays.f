c        ********************** arrays ***********************
c        *                      ------                       *
c        * this subroutine initializes the arrays  deru  and *
c        * prmt that drkgs  (runge-kutta  subroutine)  later *
c        * uses.					     *
c        *****************************************************
c
c------------------------------------------------------------------
c								   !
c        input parameters					   !
c      								   !
c      ndim - number of differential equations			   !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c								   !
c        output parameters					   !
c								   !
c      deru - vector of weight coefficients of solution error.	   !
c             later it becomes the vector of right hand side of    !
c	      the geometrical optics equations.			   !
c      prmt - vector of parameters for drkgs (initial step of	   !
c             integration, accuracy of solution etc. see also	   !
c             drkgs.for).					   !
c------------------------------------------------------------------
      subroutine arrays(ndim,deru,prmt,ihlf)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      dimension deru(*),prmt(*)

      t=1.d0/ndim
      do 100 i=1,ndim
         deru(i)=t
 100  continue
      prmt(7)=prmt(6) !arrays/Initialize: for internal control 
      !!!prmt(8)=prmt(1)+prmt(3) ! NOT USED
c      write(*,*)'in arrays prmt(1)=',prmt(1)
c      write(*,*)'in arrays prmt(2)=',prmt(2)
c      write(*,*)'in arrays tau=prmt(3)=',prmt(3)
c      write(*,*)'in arrays acccuracy=prmt(4)=',prmt(4)
c      write(*,*)'in arrays hprint=prmt(6)=',prmt(6)
c      write(*,*)'in arrays hamiltonian accuracy=prmt(9)=',prmt(9)
 20   format(f6.3/f8.3/d7.1)
 25   format(d7.1/d7.1)
 30   format(i3)

      return
      end

