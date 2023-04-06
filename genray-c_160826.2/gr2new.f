c--------------------------------------------------------------------
c***********************zrlacner*****************************************
c      it calculates the parameter t for which
c        (z(t),r(t)) are the coordinates of the point in which
c      line: z=zmag+sin(teta)*t
c            z=rmag+cos(teta)*t
c      intersects the Lackner rectangle
c----------------------------------------------------------------------
c      input data:sintet=sin(teta),costet=cos(teta)
c                       teta-poloidal angle
c                       It must be not equal 0;0.5pi;pi;1.5pi;2pi
c                  zmag,rmag -magnetic axis coordinates
c                  zmax,rmax,zmin,rmin -boundaries of Lackner rectangle
c-----------------------------------------------------------------------
c      output data:t -parameter (length along the line)
c      these data are used in     to find the contours
c-----------------------------------------------------------------------
      subroutine zrlacner(sintet,costet,zmag,rmag,zmax,rmax,rmin,zmin,
     1                    t)
      implicit integer (i-n), real*8 (a-h,o-z)
c-------------------------------------------
      zmaxd=zmax
      zmind=zmin
      rmaxd=rmax
      rmind=rmin
      zmax=zmaxd+0.01d0
      rmax=rmaxd+0.01d0
      zmin=zmind-0.01d0
      rmin=rmind-0.01d0
c----------------------------------------------------------
      if ((costet.gt.0.d0).and.(sintet.gt.0.d0)) then
c     1 quarter of the (r,z) plane
        tz=(zmax-zmag)/sintet
        tr=(rmax-rmag)/costet
	go to 1
      end if
c----------------------------------------------------------
      if ((costet.lt.0.d0).and.(sintet.gt.0.d0)) then
c     2 quarter of the (r,z) plane
        tz=(zmax-zmag)/sintet
        tr=(rmin-rmag)/costet
	go to 1
      end if
c----------------------------------------------------------
      if ((costet.lt.0.d0).and.(sintet.lt.0.d0)) then
c     3 quarter of the (r,z) plane
        tz=(zmin-zmag)/sintet
        tr=(rmin-rmag)/costet
	go to 1
      end if
c----------------------------------------------------------
      if ((costet.gt.0.d0).and.(sintet.lt.0.d0)) then
c     4 quarter of the (r,z) plane
        tz=(zmin-zmag)/sintet
        tr=(rmax-rmag)/costet
	go to 1
      end if
c----------------------------------------------------------
1       continue
        t=dmin1(tz,tr)
c----------------------------------------------------------
      rmax=rmaxd
      rmin=rmind
      zmax=zmaxd
      zmin=zmind

      return
      end





c*************************gr2new***************************************
c  It calculates coordinates for the flux functions:		
c          r(psi,teta)   z(psi,teta)     			      *
c          arrays rpsi(j,i) zpsi(j,i)  				      *
c          j=1,npsi(the number of the countours poloidal_flux=constant)	      
c          i=1,nteta+1(number of the point in the poloidal angle)    	      
c  and creates arrays ar(nl,nteta+1),az(nl,nteta+1) for nl surfaces   *
c  Used for spline (see rhospl) and plotting
c---------------------------------------------------------------------
c  input data are in common blocks /three.i/,/five.i/,/gr.i/	      *
c  nl-the number of flux contours for plotting
c---------------------------------------------------------------------*
c  Output: arrays AR,AZ, zpsi,rpsi into common block /gr/             *
c          (file gr.i)                                               *
c**********************************************************************
c
      subroutine gr2new  ! only for model_b.eq.0
      implicit real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
      include 'fourb.i' ! nnlim, rlimit
      include 'gr.i'
      integer*2 it1,it2,it3,it4,itf1,itf2,itf3,itf4
      
       
      if (maxval(rlimit).eq.0.d0) then
c        In this case:
c        2)we will create the limiter points using the close flux
c        surface psi(r,z)=psilim*psifactr, here
c        psifactr is a parameter (it must be .le.1) to avoide the
c        problems with the nonmonotonic psi function near the separatrix.
c        psifactr is given in genray.in file (it is in common/one/
c        3) the number of limiter points nlimit let equal to the number
c        of points on the magnetic surface nlimit=nteta	(see common
c        blocks: gr.cb, fourb.cb  and five
c        ------------------------------------
         write(*,*)'nnlim=0 eqdsk data without limiter points'
         psilim=psimag+(psilim-psimag)*psifactr
         write(*,*)'psisep,psilim,psifactr',psisep,psilim,psifactr
      endif
      
c-------------------------------------------------------
      hpsi=(psilim-psimag)/dble(npsi-1)
      write(*,*)'gr2new: psilim,psimag,rma,zma',psilim,psimag,rma,zma
      !pause
      do j=1,npsi
        arpsi(j)=psimag+hpsi*(j-1)
      enddo
       
      pi=4.d0*datan(1.d0)
      if( (eqdsktype.eq.'tokamak') .or. (eqdsktype.eq.'TAE') )then
        hteta=2.d0*pi/dble(nteta) 
        ! poloidal angle: [0;2pi] in a closed-surface geometry
        do i=1,nteta1
           arteta(i)=hteta*(dble(i)-0.5d0)
        enddo
      elseif(eqdsktype.eq.'mirror')then
        hteta=pi/nteta 
        ! poloidal angle: [-pi/2; +pi/2] in a mirror machine
        do i=1,nteta1 ! = 1:(nteta+1)
           arteta(i)= -0.5*pi + hteta*(i-1) ! [-pi/2; +pi/2]
        enddo
      else
        write(*,*)
     +   'gr2new: only setup for eqdsktype="tokamak","TAE","mirror" '
        stop
      endif
       
c----------------------------------------------------------
c      if ipsi=1 then continue,ipsi=0 then read file psi.bin
cbegin
c        ipsi=1 !  to calculate contours (default value)
c        ipsi=0 !  to read contours data from file:psi.bin
       if (ipsi.eq.0) then
         open(1,file='psi.bin',form='unformatted',status='old')
         do i=1,npsi
           do j=1,nteta1
             read(1)zpsi(i,j),rpsi(i,j)            
             write(*,*)'read i,j',i,j
             write(*,*)'zpsi(i,j),rpsi(i,j)',zpsi(i,j),rpsi(i,j)
           end do
         end do
         close(1)
         go to 200
       endif
       
c------------------------------------------------------------------
       if((eqdsktype.eq.'tokamak') .or. (eqdsktype.eq.'TAE'))then
         ! applicable to closed field-lines only
       
         j=1 ! magnetic axis; set (rpsi,zpsi) to (rmag,zmag) at all teta
         do  i=1,nteta1
           zpsi(j,i)=zma
           rpsi(j,i)=rma 
         end do
 
         do 10 i=1,nteta1
           teta=arteta(i)
           sintet=dsin(teta)
           costet=dcos(teta)
           pp=1.d0
           pm=1.d0
           rmaxp=rmax*pp
           rminp=rmin*pm
           zmaxp=zmax*pp
           zminp=zmin*pm

           call zrlacner(sintet,costet,zma,rma,
     1	   zmaxp,rmaxp,rminp,zminp,tm)
           !write(*,*)'gr2new-zrlacner: i,teta,tm=',i,teta,tm
           tini=0.d0
           htini=tm*0.02d0
           maxiter=10
           epspsit=epspsi
           htmin=1.d-5
         
           do 20 j=2,npsi
             psi0=arpsi(j)
c	   write(*,*)'gr2new bef zrcontor i,teta,j,psi0',i,teta,j,psi0
c------------------------------------------------
c          Newton method for solution of equation psi(r(t),z(t))=psi0
cc           call zrcontor(tini,psi0,costet,sintet,tm,htini,
cc     1     epspsit,maxiter,htmin,ierr,zt0,rt0,t0)
c------------------------------------------------
c          binary method for solution of equation psi(r(t),z(t))=psi0
c           write(*,*)'in gr2new before zrcntrbin tm,tini',tm,tini
           call zrcntrbin(tini,psi0,costet,sintet,tm,htini,
     1     ierr,zt0,rt0,t0)
c           write(*,*)'in gr2new after zrcntrbin zt0,rt0,t0',zt0,rt0,t0
c------------------------------------------------
           tini=t0
           if (ierr.eq.1) then
              zpsi(j,i)=zt0
              rpsi(j,i)=rt0
           else
              write(*,*)' gr2new: zrcntrbin gave ierr ',ierr
              write(*,*)' it is impossible to find the flux surface '
              write(*,*)' with arpsi(j)',arpsi(j),'j=',j
              write(*,*)' i=',i,'arteta(i)=',arteta(i)
              write(*,*)' with arpsi(j-1)',arpsi(j-1)
              write(*,*)' Probably got to an open field line.'
              write(*,*)' Suggestions: '
              write(*,*)'  Change psilim to arpsi(j-1) (in eqdsk file)'
              write(*,*)'  or reduce psifactr (section &tokamak)'
              !YuP[03-2016] Most likely - an open field line
              !going outside of the box [rmin,rmax;zmin,zmax].
              !Instead of stopping the code, set (rpsi,zpsi) 
              !to the edge point (of the equilibrium box) :
              zpsi(j,i)= zma+sintet*tm
              rpsi(j,i)= rma+costet*tm
              write(*,*)' Setting (rpsi,zpsi) to the edge point.'
              write(*,*)' (rpsi,zpsi)=', rpsi(j,i),zpsi(j,i)
              write(*,*)'----------------------------------------'
              !YuP: A typical situation for a mirror machine, 
              !at teta near -pi/2 or +pi/2 .
           endif
20       continue ! j (psi)
10       continue   ! i (pol.angle)

         ! closed surfaces: set the last point equal to 1st
         ! but not in a mirror machine.
         ! YuP: Actually, this part maybe not needed anymore
         ! since [03-2016] I set nteta to nteta1 in 'do 10' loop above.
         do 40 j=1,npsi
           zpsi(j,nteta+1)=zpsi(j,1)
           rpsi(j,nteta+1)=rpsi(j,1)
40       continue

       endif ! (eqdsktype.eq.'tokamak') .or. (eqdsktype.eq.'TAE')
c----------------------------------------------------------


       if(eqdsktype.eq.'mirror')then !open field lines only YuP[11-2016]
         write(*,*)'rhospl/mirror: zmin,zmax=',zmin,zmax
         write(*,*)'rhospl/mirror: rmin,rmax=',rmin,rmax
         write(*,*)'rhospl/mirror: rma,zma=',rma,zma
         dzpsi= (zmax-zmin)/(nteta)
         j=1 ! R=0 axis
         do i=1,nteta1
            zpsi(j,i)= zmin+dzpsi*(i-1) ![zmin;zmax]
            rpsi(j,i)= 0.d0
         end do
         tini=rmin ! [m] initial t along R(t)=rma+t
         tm=rmax   ! [m] largest t (see zrcntrbin)
         htini=tm*0.02d0 ! initial t-step
         do j=2,npsi
           psi0=arpsi(j) ! =psimag+hpsi*(j-1)  uniform grid
           ! For a given psi0, and given Z=zpsi, find R=rpsi
           ! such that psif_xyz(x=rpsi,y=0,z=zpsi) = psi0.
           ! For searching of rpsi, we use same subr.zrcntrbin  
           ! as for tokamaks, but with a trick:
           ! The scanning direction now is along Z=const lines,
           ! (rather than constant pol.angle)
           ! so that only R is scanned, for each given Z.
           ! The values of Z are set by zpsi(j=1,i), 
           ! and they are the same for any j (any psi0).
           ! This is based on assumption that all field lines
           ! cover the whole range of Z=[zmin;zmax],
           ! which is usually the case in a mirror machine.
           ! Index j labels different field lines, 
           ! and index i labels points along field lines.
           do i=1,nteta1
             !binary method for solution of equation psi(r(t),z(t))=psi0
             rma=0.d0
             zma=zpsi(1,i) !to make proper input for zrcntrbin
             costet=1.d0 ! which means scanning along rt1= rma+1*t1
             sintet=0.d0 ! which means scanning along zt1= zma+0*t1
             call zrcntrbin(tini,psi0,costet,sintet,tm,htini,
     1       ierr,zt0,rt0,t0)
             zpsi(j,i)=zpsi(1,i) ! same as along the axis ![zmin;zmax]
             if (ierr.eq.1) then
               !zpsi(j,i)=zt0
               rpsi(j,i)=rt0
             else
              write(*,*)' gr2new: zrcntrbin gave ierr ',ierr
              write(*,*)' it is impossible to find the flux surface '
              write(*,*)' with arpsi(j)',arpsi(j),'j=',j
              write(*,*)' i=',i,'arteta(i)=',arteta(i)
              write(*,*)' with arpsi(j-1)',arpsi(j-1)
              write(*,*)' Can be non-monotonic PSI array in R-direction'
              write(*,*)' Suggestions: '
              write(*,*)'  Change psilim to arpsi(j-1) (in eqdsk file)'
              write(*,*)'----------------------------------------'
              pause
             endif
             !write(*,'(a,2i4,2e12.4)')'j,i, rt0,zt0=', j,i,rt0,zt0
           end do ! i=1:nteta1
         end do ! j =2:npsi
         ! restore:
         rma=0.d0 ! mirror machine
         zma=0.d0 ! mirror machine
       endif ! eqdsktype.eq.'mirror'
       

c      if ipsi=1 then continue,write file psi.bin
        open(1,file='psi.bin',form='unformatted')
        do i=1,npsi
           do j=1,nteta1
             write(1)zpsi(i,j),rpsi(i,j)
           end do
        end do
        close(1)
c------------------------------------------------------------
c      if ipsi=1 then continue,ipsi=0 then read file psi.bin
200    continue
c-------------------------------------------------------------
c      do 100 j=1,npsi
c	j=npsi
c	write(*,*)'in gr2new test'
c        do 100 i=1,nteta1
c 	  write(*,*)'j,arpsi(j),i,arteta(i)',j,arpsi(j),i,arteta(i)
c 	  write(*,*)'zpsi,rpsi',zpsi(j,i),rpsi(j,i)
c	  read(*,*)
 100   continue
c----------------------------------------------------
c  the arrays for contours ploting
       jstep=npsi/NL
       do 50 j=1,NL
         do 60 i=1,nteta+1
c   multiplication to get length in cm and denormalization
	   AZ(j,i)=zpsi(jstep*j,i)*100.d0*r0x
	   AR(j,i)=rpsi(jstep*j,i)*100.d0*r0x
60	 continue
50     continue
c------------------------------------------------------
       return
       end





c*************************zrcontor*********************************** *
c  It calculates  contours coordinates of the countor point 	       *
c   r(psi,teta)   z(psi,teta) for the given:			       *
c   flux function psi(z,r)=psi0 and poloidal angle teta0(in radians)   *     			      *
c---------------------------------------------------------------------
c  input data are in common blocks /three/,/five/,/gr.cb/	       *
c  tini -initial (left) value of the t( ray parameter )		       *
c  psi0 -given value of the poloidal flux psi			       *
c  costet0,sintet0 -for the given poloidal angle		       *
c  tm   -maximal value of t (right)					       *
c  htini   ininial step of the t 	 			       *
c  maxiter - maximal number of iteratios for the equation solution     *
c  epspsit-accuracy of equation solution (=max(abs( t_n-t_n+1))        *
c  epspsi-accuracy of equation solution (=max(abs( psi_n-psi_n+1) in param.i)     *
c  htmin - the minimal parameter t step
c----------------------------------------------------------------------*
c  Output: ,zt0(psi0,teta0),rt0(psi0,teta0) and parameter t0           *
c           rt0=rma+t0*costet0 ,  zt0=zma+t0*sintet0		       *
c  ierr -index if the solution was obtained =1	else=0		       *
c**********************************************************************
      subroutine zrcontor(tini,psi0,costet0,sintet0,tm,htini,
     1 epspsit,maxiter,htmin,ierr,zt0,rt0,t0)
      implicit real*8 (a-h,o-z)
      include 'param.i'
      include 'three.i'
      ierr=1
      ht=htini
c      write(*,*)'in zrcontor psi0,costet0,sintet0',
c     1 psi0,costet0,sintet0
c      write(*,*)'tini,tm,htini,epspsit,epspsi'
c     1 ,tini,tm,htini,epspsit,epspsi
c      write(*,*)'maxiter,htmin',maxiter,htmin
 10   continue
      t1=tini
      rt1=rma+t1*costet0
      zt1=zma+t1*sintet0
      psi1=psif(zt1,rt1)
c      write(*,*)'10 t1,rt1,zt1,psi1',t1,rt1,zt1,psi1
 20   t2=t1+ht
c      write(*,*)'20 t2,t1,rt1,zt1,psi1',t2,t1,rt1,zt1,psi1
      t2=dmin1(t2,tm)
c     write(*,*)'t2',t2
      if((t2.eq.tm.and.t1.eq.tm).or.
     1   (t2.eq.0.d0.and.t1.eq.0.d0)) then
	 ierr=0
c        error exit
         goto 50
      endif

      rt2=rma+t2*costet0
      zt2=zma+t2*sintet0
      psi2=psif(zt2,rt2)
      dpsi=(psi0-psi1)*(psi0-psi2)
c      write(*,*)'20 t2,rt2,zt2,psi2,dpsi',t2,rt2,zt2,psi2,dpsi
      if(dpsi.le.0.d0)go to 30
      t1=t2
      psi1=psi2
      goto 20
c-----------------------------------------------------------
c     psi0 is between psi1 and psi2
c     iteration methode
c-----------------------------------------------------------
 30   numbiter=0
      tn=t1+ht*0.5d0
      rtn=rma+tn*costet0
      ztn=zma+tn*sintet0
c      write(*,*)'30 tn,rtn,ztn',tn,rtn,ztn
 40   continue
      psin=psif(ztn,rtn)
      call dpsidzdr(z,r,phi,dpsidz,dpsidr)
c      write(*,*)'40 tn,rtn,ztn,psin',tn,rtn,ztn,psin
      dzdt=sintet0
      drdt=costet0
      dfpsidt=dpsidz*dzdt+dpsidr*drdt
c      write(*,*)'40 dpsidz,dpsidr,dfdpidt',dpsidz,dpsidr,dfpsidt
      if (dfpsidt.eq.0.d0) then
c         write(*,*)'in zrcontur dfdpsidt=0 psi0=',psi0,
c     1    'costet0=',costet0
         stop
      endif
      dpsi=psin-psi0
      delt=-dpsi/dfpsidt
      tn=tn+delt
      rtn=rma+tn*costet0
      ztn=zma+tn*sintet0
      if(dabs(delt).lt.epspsit) then
c        write(*,*)'1 delt,epspsit',delt,epspsit
c	   write(*,*)'delt,tn',delt,tn
c	   write(*,*)'psin-psi0,psin,psi0',psin-psi0,psin,psi0
c	   write(*,*)'rtn,ztn',rtn,ztn
        go to 50
      endif
      if(psi0.eq.0.and.dabs(dpsi).lt.epspsi) then
c      	write(*,*)'2 psi0,dpsi,epspsi',psi0,dpsi,epspsi
        go to 50
      endif
c      if(dabs(psin-psi0).lt.epspsi)then
c        write(*,*)'4 psin-psi0',psin-psi0
c	   write(*,*)'delt,tn',delt,tn
c	   write(*,*)'psin-psi0,psin,psi0',psin-psi0,psin,psi0
c	   write(*,*)'rtn,ztn',rtn,ztn
c        go to 50
c      endif
      if(psi0.ne.0.)  then
        if(dabs(dpsi/psi0).lt.epspsi)then
c	  write(*,*)'3 dpsi/psi0,epspsi',dpsi/psi0,epspsi
c          write(*,*)'psin-psi0',psin-psi0
	  go to 50
	endif
      endif
c------------------------------------------------
      numbiter=numbiter+1
      if(numbiter.ge.maxiter)then
         if(ht.le.htmin) stop 'zrcontur'
	 ht=0.5d0*ht
	 goto 10
      endif
      goto 40
c------------------------------------------------
c     end of iterations
c------------------------------------------------
 50   continue
      t0=tn
      zt0=ztn
      rt0=rtn
      return
      end

c*************************zrcntrbin*********************************** *
c  It calculates  the contours coordinates of the countour points      *
c   r(psi,teta)   z(psi,teta) for the given:			       *
c   flux function psi(z,r)=psi0 and poloidal angle teta0(in radians)   *
c   It using the binary methode                                        *
c---------------------------------------------------------------------
c  input data are in common blocks /three/,/five/,/gr.cb/	       *
c  tini -initial (left) value of the t( ray parameter )		       *
c  psi0 -given value of the poloidal flux psi			       *
c  costet0,sintet0 -for the given poloidal angle		       *
c  tm   -maximal value of t (right)				       *
c  htini   initial step of the t 	 			       *
c  epspsi-accuracy of equation solution (=max(abs( psi_n-psi_n+1))(in param.i)
c----------------------------------------------------------------------*
c  Output: ,zt0(psi0,teta0),rt0(psi0,teta0) and parameter t0           *
c           rt0=rma+t0*costet0 ,  zt0=zma+t0*sintet0		       *
c  ierr -index if the solution was obtained =1	else=0		       *
c**********************************************************************
      subroutine zrcntrbin(tini,psi0,costet0,sintet0,tm,htini,
     1 ierr,zt0,rt0,t0)
      implicit real*8 (a-h,o-z)
      include 'param.i'
      include 'three.i'
      ierr=1
      ht=htini
c      write(*,*)'in zrcontor !1 psi0,costet0,sintet0',
c     1 psi0,costet0,sintet0
c      write(*,*)'tini,tm,htini,epspsi'
c     1 ,tini,tm,htini,epspsi
      t1=tini
      rt1=rma+t1*costet0
      zt1=zma+t1*sintet0
      psi1=psif_xyz(rt1,0.d0,zt1) ! y=0  ok for now?
 10   t2=t1+ht
c      write(*,*)'10 t2,t1,rt1,zt1,psi1',t2,t1,rt1,zt1,psi1
      t2=dmin1(t2,tm)
c      write(*,*)'t2',t2
      if((t2.eq.tm.and.t1.eq.tm).or.
     1   (t2.eq.0.d0.and.t1.eq.0.d0)) then
	 ierr=0
c	 write(*,*)'in zrcntrbin error t2,tm,t1',t2,tm,t1
c	 write(*,*)'in zrcntrbin error exit ierr',ierr
c        error exit
         goto 30
      endif

      rt2=rma+t2*costet0
      zt2=zma+t2*sintet0
      psi2=psif_xyz(rt2,0.d0,zt2)  ! y=0  ok for now?
      dpsi=(psi0-psi1)*(psi0-psi2)
c      write(*,*)'10 t2,rt2,zt2,psi2,dpsi',t2,rt2,zt2,psi2,dpsi
      if(dpsi.le.0.d0)go to 20
      t1=t2
      psi1=psi2
      goto 10
c-----------------------------------------------------------
c     psi0 is between psi1 and psi2
c     binary iteration methode
c-----------------------------------------------------------
 20   continue
      tr=t2
      tl=t1
      do while ((tr-tl).gt.epspsi)
         t=tl+(tr-tl)*0.5d0
         r=rma+t*costet0
         z=zma+t*sintet0
         psi1=psif_xyz(r,0.d0,z)-psi0 ! y=0  ok for now?
         rtr=rma+tr*costet0
         ztr=zma+tr*sintet0
         psi2=psif_xyz(rtr,0.d0,ztr)-psi0 ! y=0  ok for now?
         if ((psi1*psi2).gt.0) then
            tr=t
         else
            tl=t
         end if
      end do
c     -----------------------------------------------
c          end of the binary methode
c     -----------------------------------------------
      t0=t
      zt0=z
      rt0=r
 30   continue
c      write(*,*)'the end of zrcotrbin t0',t0
      return
      end





c        **********************dpsidzdr***********************
c        *                        -                           *
c        * this subroutine calculates the derivatives         *
c        * output:dpsidz and dpsidr           		      *
c        ******************************************************
c
c------------------------------------------------------------------
c								   !
c        input parameters					   !
c      								   !
c      z, r, phi - point which the components of magnetic field	   !
c                  are calculated in.				   !
c------------------------------------------------------------------
      subroutine dpsidzdr(z,r,phi,dpsidz,dpsidr)
      implicit integer (i-n), real*8 (a-h,o-z)

      include 'param.i'
      include 'five.i'
      double precision ias1r,ias2r,ias2r_Sm
c
c nr4 , nz4 and nry  are given in common five as parameters
c
cSm030224
      nr4=nx+4
      nz4=ny+4
      ncx=nr4
      ncy=nz4
c      psi=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,0,0,r,z)
c      dpsidr=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,1,0,r,z)
c      dpsidz=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,0,1,r,z)
      dpsidr=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,1,0,r,z,nr4a)
      dpsidz=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,0,1,r,z,nr4a)
      return
      end

