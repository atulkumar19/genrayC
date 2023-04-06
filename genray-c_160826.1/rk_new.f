
c======================================================================
c======================================================================


c        4_th order  Runge-Kutta method with automatic 
c        time step selection
c
c        **********************     drkgs2   ****************
c        *                      -----                       
c        * this Runge-Kutta subroutine finds the solution of the 
c        *      system of ordinary differential equations          
c        ****************************************************
c
c
c-----------------------------------------------------------------
c                                                                !
c       1. method of solution                                    !
c                                                                !
c             runge-kutta method  of  4th-order                  !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
c                                                                !
c       2. input parameters                                      !
c                                                                !
c          prmt  - the vector of input and output data:          !
c              prmt(1) - the lower boundary of interval of       !
c                        integration;                            !
c              prmt(2) - the upper boundary ...;                 !
c              prmt(3) - the initial step of integration;        !
c              prmt(4) - the accuracy of solution;               !
c              prmt(5) - if it is not equal to zero the control  !
c                        is transferred to the main program;     !
c                                                                !
c          y     - the input vector of initial values. later     !
c                  it becomes the vector of results obtained     !
c                  at the intermediate points x;                 !
c                                                                !
c          dery  - the  vector of derivatives of y at point x	 !
c                                                                !
c          ndim  - the number of system equations;               !
c                                                                !
c          aux   - the auxiliary array.                          !
c                                                                !
c          i_output =1 the output at the poloidal distance steps !
c                   =2 the output at the total distance steps    !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
c                                                                !
c       3. subroutines and functions used                        !
c                                                                !
c          fct - calculates the right hand sides of  the  ode    !
c                system. it has not to change the  values  of    !
c                x and y. its formal parameters are:             !
c                x, y, dery.                                     !
c                                                                !
c          outp - output subroutine. it has not to change the    !
c                 values of all  its  formal  parameters.  if    !
c                 prmt(5) is not equal to zero,  the  control    !
c                 will be transferred to the main program.       !
c                 its formal parameters are:			 !
c                 x, y, dery, ihlf, ndim, prmt                   !
c          output date						 !
c    u(1)=z,u(2)=r,u(3)=phi,u(4)=nz,u(5)=nr,u(6)=m               !                                         !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
      subroutine drkgs2(prmt,u,deru,ndim,ihlf,fct,outp,aux,i_output)
      implicit integer (i-n), real*8 (a-h,o-z)
      dimension u(6),deru(6),aux(8,6),prmt(9),up(6),uu(6),startu(6)
      integer i_output
      dimension uplus(6)
      double precision prmt,u,deru,aux,t,tend,h,up,tt,hh,uu,startu,
     1 dd,uplus,eps
      double precision us,uz,ur,dels1,dels2,uphi
      integer iflag
      
      external fct !=rside1 (here: in original drkgs2)

      t=prmt(1)
      tt=t
      tend=prmt(2)
      h=prmt(3)
      hh=h/2.0
      prmt(5)=0.d0
      eps=prmt(4)
      iflag=0  ! control of the RK method accuracy
      iflagh=3 ! control of the plasma edg intersecsion
      us=0.d0  ! poloidal length (initial value)
c---------------------------------------------------------------
 10   continue

      if (dabs(h).lt.1.d-11) then
         write(*,*)'***** In Runge-Kutta subroutine drkgs2 *********'
         write(*,*)'***** it was too small time step h.lt.1.d-11 ***' 
         write(*,*)'***** can not to get the given accuracy prmt4 **' 
         write(*,*)'h',h
         iraystop=1
         goto 100
      endif    

      uz=u(1)
      ur=u(2)
      uphi=u(3)
      do i=1,ndim
	uu(i)=u(i)
	startu(i)=u(i)
      enddo
      tt=t
      hh=h/2.0d0
      call fct(u,deru) !fct(t,u,deru)
      do i=1,ndim	       
        aux(1,i)=h*deru(i)
	  aux(5,i)=u(i)+0.5*aux(1,i)
 	  up(i)=aux(5,i)
      enddo
c------------------------------------------------------------------
      call fct(up,deru) !fct(t+0.5*h,up,deru)
      do i=1,ndim
         aux(2,i)=h*deru(i)
         aux(5,i)=u(i)+0.5*aux(2,i)
         up(i)=aux(5,i)
      enddo

      call fct(up,deru) !fct(t+0.5*h,up,deru)
      do i=1,ndim
         aux(3,i)=h*deru(i)
         aux(5,i)=u(i)+aux(3,i)
         up(i)=aux(5,i)
      enddo
      call fct(up,deru) !fct(t+h,up,deru)
      do i=1,ndim
        aux(4,i)=h*deru(i)
      enddo
c--------------------------------------------------------------------
      do i=1,ndim
        up(i)=u(i)+1.d0/6.d0*
     1        (aux(1,i)+2.d0*aux(2,i)+2.d0*aux(3,i)+aux(4,i))
      enddo
c--------------------------------------------------------------------
      do i=1,ndim
         u(i)=up(i)
      enddo
c=====================================================================
      do 320 j=1,2
      call fct(uu,deru) ! fct(tt,uu,deru)

      do i=1,ndim
         aux(1,i)=hh*deru(i)
	 aux(5,i)=uu(i)+0.5*aux(1,i)
 	 up(i)=aux(5,i)
      enddo
c--------------------------------------------------------------------
      call fct(up,deru) !fct(tt+0.5*hh,up,deru)
      do i=1,ndim
         aux(2,i)=hh*deru(i)
         aux(5,i)=uu(i)+0.5*aux(2,i)
 	 up(i)=aux(5,i)
      enddo
c--------------------------------------------------------------------
      call fct(up,deru) !fct(tt+0.5*hh,up,deru)
      do i=1,ndim
         aux(3,i)=hh*deru(i)
         aux(5,i)=uu(i)+aux(3,i)
 	 up(i)=aux(5,i)
      enddo
c--------------------------------------------------------------------
      call fct(up,deru) !fct(tt+hh,up,deru)
      do i=1,ndim
        aux(4,i)=hh*deru(i)
      enddo
c--------------------------------------------------------------------
      do i=1,ndim
        uu(i)=uu(i)+1./6.*(aux(1,i)+2.*aux(2,i)+2.*aux(3,i)+aux(4,i))
      enddo
c--------------------------------------------------------------------
      tt=tt+hh
c--------------------------------------------------------------------
 320  continue

c=====================================================================
      dd=0.0d0
      
      do 80 i=1,ndim

c        if(i.eq.3) goto 80
cc--------relative error
         dp=dabs(u(i))
         if (dp.gt.1.d0) then
            dp=1.d0/dp
         else
            dp=1.d0
         endif
         dd=dd+(u(i)-uu(i))*(u(i)-uu(i))*dp
cc--------absolute error
c         dd=dd+(u(i)-uu(i))*(u(i)-uu(i))

 80   continue  
      
      dd=dsqrt(dd)

      eps_ham =prmt(9)
      call callc_eps_ham(uu,ham)

      if ((dd .lt. eps).and.(dabs(ham).lt.eps_ham)) goto 189
c      if (dd .lt. eps) goto 189
         
      if (iflag .le. 0) then
        h=hh
        do i=1,ndim
           u(i)=startu(i)
        enddo
c         write(*,*)'!!!!!!!!!!!!!!!! back !!!!!!!!!!!!!!!'
        iflag=-1
        goto 10
      else
        h=hh
        do i=1,ndim
           u(i)=uplus(i)
        enddo
c         write(*,*)'!!!!!!!!!!!!!!!! done !!!!!!!!!!!!!!!'
        iflag=2
      endif

 189  continue
      if (iflagh.eq. 2) goto 191
      if (iflag .lt. 0) goto 191
cSAP080807
       if ((dd .lt. (eps/16)).and. (dabs(ham).lt.prmt(9))) then
       h=h*2.0d0
        do i=1,ndim
           uplus(i)=u(i)
        enddo
        do i=1,ndim
           u(i)=startu(i)
        enddo
c         write(*,*)'!!!!!!!!!!!!!!!! back++ !!!!!!!!!!!!!!!'
        iflag=1
        goto 10
      endif

c=====================================================================

191   continue

c      t=t+h
      iflag=0
      uout1=u(1)
      uout2=u(2)
      uout3=u(3)
      uout4=u(4)
      uout5=u(5)
      uout6=u(6)
c-----------------------------------------------------------
c      write(*,*)'in rk before outpt deru(1),deru(2)',deru(1),deru(2)

      if (i_output.eq.1) then
c        poloidal distance
         dels2=(uz-u(1))*(uz-u(1))+(ur-u(2))*(ur-u(2))
      endif

      if (i_output.eq.2) then
c       total distance
         dels2=(uz-u(1))**2+(ur-u(2))**2+(ur*uphi-u(2)*u(3))**2
      endif

      dels1=dsqrt(dels2)

      usold=us
 
      us=us+dels1 ! the distance after the Runge-Kutta with the variable time step  

c-----check that the ray point is close to the output point
      delh=1.d-1     
 
      i_delh_1=0
 20   continue
      i_delh=0

      if(us.gt.(prmt(7)+delh*prmt(6))) then

c-------trajectory jumped past the output point us > prmt(7)+delh*prmt(6)
c-------We will reduse the time step to get the point close to the given output point 
        i_delh=1
        hold=h 
 
        if (i_delh_1.ne.1) hnew=h*(prmt(7)+delh*prmt(6)-usold)/(dels1)
        i_delh_1=1

c-------one step using the Runge-Kutta with the constant time step h=hnew

c       go back (to the previous time step) to start of the Runge-Kutta procedure 
     
        do i=1,ndim
	   u(i)=startu(i)
        enddo
c------ one step Runge-Kutta procedure    
        write(*,*)'rk_new before drkgs0 hnew,u(1),u(2)',hnew,u(1),u(2)
        call drkgs0(hnew,u,deru,ndim,fct,aux)
        write(*,*)'rk_new after drkgs0 u(1),u(2)',u(1),u(2)
        call boundc(u(1),u(2),iboundc )
      endif
c-----------------------------------------------------------
      if (i_output.eq.1) then
c        poloidal distance after one step using drkgs0
         dels2=(uz-u(1))*(uz-u(1))+(ur-u(2))*(ur-u(2))
      endif

      if (i_output.eq.2) then
c        total distance
         dels2=(uz-u(1))**2+(ur-u(2))**2+(ur*uphi-u(2)*u(3))**2
      endif

      dels1=dsqrt(dels2)

      us=usold
      us=us+dels1 ! the distance after one step with h=hnew
      if (i_delh.eq.0) goto 60
      goto 60  
      if(us.lt.(prmt(7)+delh*prmt(6))) then 
c--------hnew is too shot to jump past prmt(7)+delhprmt(6).
c        Now we will increase hnew.
          write(*,*)' hold,hnew',hold,hnew
          hnew=0.5d0*(hold+hnew)
          write(*,*)'increase hnew hnew,hold',hnew,hold
          goto 20      
       endif
 60   continue
      call outp(us,u,deru,ihlf,ndim,prmt,iflagh,iraystop)
c-----------------------------------------------------------
      if (iflagh.eq.1) then !if 1
c----------------------------------------------------------------
c        ray is near the plasma boundary after Runge-Kutta
c        procedure (it is after reflection)
c-----------------------------------------------------------------
      else
        if (iflagh.eq.2) then !if 2
c-------------------------------------------------------------------
c         Ray is outside  the plasma after the correction.
c         It is necessary to reduce the time step in the Runge-Kutta
c         procedure and to recalculate the array u(i)
c----------------------------------------------------------------------
          t=t-h
          h=h*0.5d0
	  dtstep=prmt(3)/h
	  write(*,*)'ray is outside the plasma after correction'
c------------------------------------------------------------------
	  if (dtstep.lt.100.d0)then !if 2.1
	    do i=1,ndim
	      u(i)=startu(i)
	    enddo
	    goto 10
          else
            write(*,*)'dtstep.gt.100.0'
	    do i=1,ndim
	      u(i)=startu(i)
	    enddo
	    goto 10
	  endif !end if 2.1
c-------------------------------------------------------------------
        else
c-------------------------------------------------------------------
c     ordinary ray point	  iflagh=3
c     value of u(i) after correction procedure
c-------------------------------------------------------------------
        endif !end if 2
      endif !end if 1
c--------------------------------------------------------------------
      if (iraystop.eq.1) then
        goto 100
      end if     
      if (prmt(5).gt.0) goto 100
  30      format(1x,6(1x,e11.4))
c------------------------------------------------------------------
      goto 10
  100 continue
      return
      end ! end of drkgs2 ============================================
      
      
      

c======================================================================
c======================================================================


c        ********************** drkgs0    *******************
c         4_th order Runge-Kutta subroutine with constant time 
c         step finds the solution of
c         the system of the ordinary differential equations
c         only for one time step         
c        ****************************************************
c
c-----------------------------------------------------------------
c                                                                !
c       1. method of solution                                    !
c                                                                !
c          Runge-Kutta method  of  4th-order                     !
c          with the constant step integration
c          It makes only one step!
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
c                                                                !
c       2. input parameters                                      !
c                                                                !
c          h     - time step                                     ! 
c                                                                !
c          u     - the input vector of initial values. later     !
c                  it becomes the vector of results obtained     !
c                  after one time stept                         !
c                                                                !
c          deru  - the  vector of derivatives of u at point t	 !
c                                                                !
c          ndim  - the number of system equations;               !
c                                                                !
c          aux   - the auxiliary array.                          !
c                                                                !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
c                                                                !
c       3. subroutines and functions used                        !
c                                                                !
c          fct - calculates the right hand sides of  the  ode    !
c                system.   It does not change the  values  of    !
c                u. its formal parameters are:             
c                u, deru.    NO t dependence anymore YuP[03-2016]                       !
c                                                                !
c       			  		 
c    u(1)=z,u(2)=r,u(3)=phi,u(4)=nz,u(5)=nr,u(6)=m   (for fct=rside1)
c    or
c    u(1)=x,u(2)=y,u(3)=z,u(4)=nx,u(5)=ny,u(6)=nz   (for fct=rside_xyz)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
      subroutine drkgs0(h,u,deru,ndim,fct,aux) ! fct=rside_xyz or rside1
      implicit integer (i-n), real*8 (a-h,o-z)

      dimension u(*),deru(*),aux(8,*),up(6)
      double precision u,deru,aux,h,up

      call fct(u,deru)  !fct(t,u,deru)
      do i=1,ndim
         aux(1,i)=h*deru(i)
         aux(5,i)=u(i)+0.5d0*aux(1,i)
         up(i)=aux(5,i)
      enddo
c------------------------------------------------------------------
      call fct(up,deru) !fct(t+0.5d0*h,up,deru)
      do i=1,ndim
         aux(2,i)=h*deru(i)
         aux(5,i)=u(i)+0.5d0*aux(2,i)
         up(i)=aux(5,i)
      enddo
c--------------------------------------------------------------------
      call fct(up,deru) !fct(t+0.5d0*h,up,deru)
      do i=1,ndim
         aux(3,i)=h*deru(i)
         aux(5,i)=u(i)+aux(3,i)
         up(i)=aux(5,i)
      enddo
c--------------------------------------------------------------------
      call fct(up,deru) !fct(t+h,up,deru)
      do i=1,ndim
         aux(4,i)=h*deru(i)
      enddo
c--------------------------------------------------------------------
      do i=1,ndim
         up(i)=u(i)+1.d0/6.d0*(aux(1,i)+2.d0*aux(2,i)+2.d0*aux(3,i)
     *          +aux(4,i))
      enddo

      do i=1,ndim
         u(i)=up(i)
      enddo

      return
      end !  drkgs0

c======================================================================
c======================================================================

      subroutine callc_eps_ham(u,ham)
c-------------------------------------
c     calculates hamiltonian
c     to check rk accuracy
c-----------------------------
      implicit none
      include 'param.i'
      include 'one.i'
c-----input
      real*8 u(6)    !ray coordinates
c-----output
      real*8 ham
c-----externals
      real*8 b,gamma1,hamilt1

c-----local

      bmod=b(u(1),u(2),u(3))
      gam=gamma1(u(1),u(2),u(3),u(4),u(5),u(6))
      ham=hamilt1(u(1),u(2),u(3),u(4),u(5),u(6))

      return
      end

c======================================================================
c======================================================================
