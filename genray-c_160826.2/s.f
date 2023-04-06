c        **********************   s   ************************
c        *                        -                          *
c        * this subroutine calculates  six  different  sums  *
c        * that are used by rside for cold plasma dispersion *
c        *****************************************************
c
c-------------------------------------------------------------------
c								    !
c        input parameters					    !
c      								    !
c      z, r phi - coordinates of  the  point  where  the  angle  is !
c                 calculated.      		         	    !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
c								    !
c        output parameters                                          !
c								    !
c     s1 = 1 - sum_(i.ne.ib) { x_i/[1-(y_i)**2] }, at    ib.ne.1    !
c								    !
c     s2 = 1 - sum_(i.ne.ib) { x_i/[1-y_i] },      at    ib.ne.1    !
c  								    !
c     s3 = 1 - sum_(i) { x_i/[1+y_i] },      		            !
c								    !
c     s4 = 1 - sum_(i) { x_i } - x_1 ,      		            !
c								    !
c     s6 = 1 - sum_(i) {x_i/[1-y_i] ,	           at    ib.eq.1    !
c								    !
c     s7 = 1 - sum_(i) {x_i/[1-(y_i)**2] ,	   at    ib.eq.1    !
c								    !
c     i=2,nbulk; ib=1,nbulk; nbulk - number of bulk particles       !
c-------------------------------------------------------------------
      subroutine s(z,r,phi,s1,s2,s3,s4,s6,s7)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
        ds1=0.d0
        ds2=0.d0
        ds3=0.d0
        ds4=0.d0
        ds6=0.d0
        ds7=0.d0

      if (nbulk.gt.1) then
        do 10 i=2,nbulk

          xi=x(z,r,phi,i)
          yi=y(z,r,phi,i)
          ds3=ds3+xi/(1.d0+yi)
          ds4=ds4+xi
         
          if (ib.eq.1) then
            ds6=ds6+xi/(1.d0-yi)
            ds7=ds7+xi/(1.d0-yi*yi)             
          end if
          
          if ((i.eq.ib).or.(ib.eq.1)) goto 10
          
          ds1=ds1+xi/(1.d0-yi*yi)
          ds2=ds2+xi/(1.d0-yi)  
                       
 10     continue
      end if ! finish if ef nbulk.gt.1
      
        s1=1.d0-ds1
        s2=1.d0-ds2
        s3=1.d0-ds3      
        pps=x(z,r,phi,1)      
        s4=1.d0-ds4-pps
        s6=1.d0-ds6
        s7=1.d0-ds7
       
      return
      end
