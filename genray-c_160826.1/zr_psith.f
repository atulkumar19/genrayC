c*************************zr_psith************************************
c  this subroutine calculates z(psi,theta) and r(psi,theta) -they are *      
c            bicubic spline approximation			      *
c								      *
c  these program uses the following subroutines:  		      *
c      z=ias2r(txz,npsi,tyz,nteta1,czxy,npsi4,nteta14,0,0,psi,teta)   *
c**********************************************************************
c  input parameters:
c  psi   poloidal flux,
c  teta  poloidal angle(radian)   
c  input data from common blocks:gr.i, rho.i     
c**********************************************************************
c  output parameters: z(psi,theta),r(psi,theta)			      *
c----------------------------------------------------------------------
      subroutine zr_psith(psi,teta,z,r)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'gr.i'
      include 'rho.i'

      double precision ias2r
      double precision ias1r

c      write(*,*)'in zr_psith psi,teta',psi,teta

      if(model_b.ne.0) then
        stop 'zr_psith: Not setup for model_b>0'
      endif

      pi=4.d0*datan(1.d0)

c-----transform the angle teta to the region 0=< teta1 <2*pi
      ratio=teta/(2.d0*pi)
      teta1=teta

      if (ratio.gt.1) then
          k=int(ratio)
          teta1=teta-2.d0*pi*k
      endif

      if (ratio.lt.0) then
          k=int(ratio)-1
          teta1=teta-2.d0*pi*k
      endif
       
c      write(*,*)'in zr_psith psi,teta,teta1',psi,teta,teta1
c-----------------------------------------------------------------
c     z(psi,teta) ,dz/dpsi(psi,teta), dz/dteta(pis,teta)
c-----------------------------------------------------------------
              z=ias2r(txz,npsi,tyz,nteta1,czxy,npsi4,nteta14,
     1   	          0,0,psi,teta1)
c-----------------------------------------------------------------
c     r(psi,teta) ,dr/dpsi(psi,teta), dr/dteta(pis,teta)
c-----------------------------------------------------------------
              r=ias2r(txr,npsi,tyr,nteta1,crxy,npsi4,nteta14,
     1	              0,0,psi,teta1)
c------------------------------------------------------------------
c    Jacobian=ABS(dr/dpsi*dz/dteta-dr/dteta*dz/dpsi)
c------------------------------------------------------------------
c      	      cjacobn=dabs(drdpsi*dzdteta-drdteta*dzdpsi)
c	      write(*,*)'cjacobn',cjacobn
c-------------------------------------------------------------------
c     toroidal magnetic field b_tor=feqd(psi)/r
c-------------------------------------------------------------------
c         nr4 is in common five
c         nx=npsi is in common five
c-------------------------------------------------------------------
c         idx=0
c         feqd=ias1r(txf,nx,nr4,cx,idx,psi)
c         btor=feqd/r
c         write(*,*)'b_tor',btor
c------------------------------------------------------------------
     
        return
        end




c*************************rbmax**********************************
c  The calculations of the arrays
c  ar_min(npsi)=min{teta}(r(psi,teta))-min major radius at the given
c                                     flux surface psi(r,z)=const 				      
c  ar_max(npsi)=max{teta}(r(psi,teta))-max major radius at the given
c                                     flux surface psi(r,z)=const 				      
c  ab_max(npsi)=max{teta}(b(psi,teta))-max b_total(magnetic field)
c                        at the given flux surface psi(r,z)=const 
c  ab_min(npsi)=max{teta}(b(psi,teta))-min b_total(magnetic field)
c                        at the given flux surface psi(r,z)=const 
c******************************************************************
c     input :  from gr.i					      *
c     It uses contours coordinates for flux functions:				      *
c          r(psi,teta)   z(psi,teta)     			      *
c          arrays rpsi(j,i) zpsi(j,i)  				      *
c          j=1,npsi(number of counturs poloidal flux=constant)	      *
c          i=1,nteta+1(number of point in poloidal angle mesh)        *
c          this arays were calculated in gr2new			      *
c          rpsi(j,i)=r(arpsi(j),arteta(i))                 	      *
c          zpsi(j,i)=z(arpsi(j),arteta(i))                 	      *
c          arpsi(npsi),arteta(nteta1)				      *
c     It uses  function b(z,r,phi), ATTENTION  !!!!!!		      *
c     All data for b() must be calculated before CALL fluxrb	      *			      *
c------------------------------------------------------------------
c     output: arrays ar_min(npsi),ar_max(npsi),ab_max(npsi) in gr.i
c------------------------------------------------------------------
      subroutine rbmax
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'gr.i'

      phi=0.d0
c      write(*,*)'in rbmax'
      do j=1,npsi
        z=zpsi(j,1)
        r=rpsi(j,1)
c        write(*,*)'in rrbmax j,z,r',j,z,r
        bmod=b(z,r,phi)
c        write(*,*)'in rbmax 1 bmod=',bmod
        ar_max(j)=r
        ar_min(j)=r
        ab_max(j)=bmod
        ab_min(j)=bmod
        do i=1,nteta
c          teta=arteta(i)
	   z=zpsi(j,i)
	   r=rpsi(j,i)
c           write(*,*)'i z,r',i,z,r
c           write(*,*)'j,i,arpsi(j),arteta(i)',j,i,arpsi(j),arteta(i)
	   bmod=b(z,r,phi)
	   if (ar_min(j).gt.r) ar_min(j)=r
	   if (ar_max(j).lt.r) ar_max(j)=r
	   if (ab_max(j).lt.bmod) ab_max(j)=bmod
	   if (ab_min(j).gt.bmod) ab_min(j)=bmod
c	   write(*,*)'j,psi',j,arpsi(j)
c	   write(*,*)'i,teta',i,arteta(i)
c	   write(*,*)'ar_min(j),r,ar_max(j)',ar_min(j),r,ar_max(j)
c	   write(*,*)'ab_max(j),bmod',ab_max(j),bmod
c	   write(*,*)'ab_min(j),bmod',ab_min(j),bmod
        enddo ! i=1,nteta
      enddo ! j=1,npsi

c     arrays ar_min(npsi),ar_min(npsi),ab_max(npsi),ab_min(npsi)
c     have been created.They are in gr.i
      return
      end




       double precision
     1 function rmax_psi(psi)	!(in (m))
c------------------------------------------------------------------
c     this function calculates max major radius on the given 
c     flux surface 'psi' if psi>psilim it gives rmax_psi(psilim)=rmax 
c---------------------------------------------------------------------
c     input parameter: poloidal flux 'psi'
c     trmaxpsi(npsi4),crmaxpsi(npsi4) in rho.i
c---------------------------------------------------------------------
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'gr.i'
      include 'rho.i'
      include 'three.i'
      include 'five.i'
      double precision ias1r
      
      if(model_b.ne.0) then
        stop 'rmax_psi: Not setup for model_b>0'
      endif

      if (psi.gt.psilim) then
c         the point is outside the plasma
	  rmax_psi=rmax
      else
	  idx=0
	  rmax_psi=ias1r(trmaxpsi,npsi,npsi4,crmaxpsi,idx,psi)
      endif
      return
      end




       double precision
     1 function rmin_psi(psi)  ! in (m)
c------------------------------------------------------------------
c     this function calculates min major radius on the given 
c     flux surface 'psi' if psi>psilim it gives rmin_psi(psilim)=rmax 
c---------------------------------------------------------------------
c     input parameter: poloidal flux 'psi'
c     trminpsi(npsi4),crminxpsi(npsi4) in rho.i
c---------------------------------------------------------------------
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'gr.i'
      include 'rho.i'
      include 'three.i'
      include 'five.i'
      double precision ias1r

      if(model_b.ne.0) then
        stop 'rmax_psi: Not setup for model_b>0'
      endif

      if (psi.gt.psilim) then
c        the point is outside the plasma
	 rmin_psi=rmin
      else
	 idx=0
	 rmin_psi=ias1r(trminpsi,npsi,npsi4,crminpsi,idx,psi)
      endif
      return
      end




       double precision
     1 function bmax_psi(psi)	 ! (in Tl)
c------------------------------------------------------------------
c     this function calculates max b_tot given flux surface 'psi' 
c---------------------------------------------------------------------
c     input parameter: poloidal flux 'psi'
c     tbmaxpsi(npsi4),cbmaxpsi(npsi4) in rho.i
c---------------------------------------------------------------------
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'gr.i'
      include 'rho.i'
      double precision ias1r

      if(model_b.ne.0) then
        stop 'rmax_psi: Not setup for model_b>0'
      endif

      idx=0
      bmax_psi=ias1r(tbmaxpsi,npsi,npsi4,cbmaxpsi,idx,psi)
      return
      end

      double precision
     1 function bmin_psi(psi)	 ! (in Tl)
c------------------------------------------------------------------
c     this function calculates min b_tot given flux surface 'psi' 
c---------------------------------------------------------------------
c     input parameter: poloidal flux 'psi'
c     tbminpsi(npsi4),cbminpsi(npsi4) in rho.i
c---------------------------------------------------------------------
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'gr.i'
      include 'rho.i'
      double precision ias1r

      if(model_b.ne.0) then
        stop 'rmax_psi: Not setup for model_b>0'
      endif

      idx=0
      bmin_psi=ias1r(tbminpsi,npsi,npsi4,cbminpsi,idx,psi)
      return
      end




      double precision function psi0_minus_psi_zr(r)
c--------------------------------------------------
c     calculates function: psi0_minus_psi_ar=psi_0-psi(z_0,r)
c
c     Input data:
c     psi_0 is the value of the poloidal flux at the limiter
c           It is given common common_psi_minus_psi_zr
c     z_0   is the given vertical coordinate
c           It is given common common_psi_minus_psi_zr
c--------------------------------------------------
      implicit none
c-----input
       double precision
     & r          ! major radius

       double precision
     & psi_0,     !poloidal flux from /common_psi_m_psi_zr/ 
     & z_0        !vertical coordinate from /common_psi_m_psi_zr/ 
      common /common_psi_minus_psi_zr/ psi_0,z_0 ! input data
c-----externals
      double precision psif
 
      psi0_minus_psi_zr=psi_0-psif(z_0,r)

      return
      end




cSm070131
      subroutine r_min_r_max_from_psi_z(psi0,z0,accuracy,
     & rmin_psi_z,rmax_psi_z)
c----------------------------------------------------------
c     calculates minimal and maximal major radiuses: 
c     r=rmin_psi_z < r=rmax_psi_z
c     at the flux surface with the given poloidal flux: psi0
c     at the given vertical coordinate: z0
c----------------------------------------------------------
      implicit none
c-----input
      double precision
     & psi0,     ! poloidal flux
     & z0,       ! vertical coordinate
     & accuracy  ! accuracy of roots
                 ! {r=rmin_psi_z, r=rmax_psi_z}
                 ! calculation 
                 ! for the equation: psi_0-psi(z_0,r)=0
     

      include 'param.i'
      include 'three.i'       
      include 'five.i'
c-----output
      double precision rmin_psi_z,rmax_psi_z !roots

c-----externals
      double precision  rtbis,psi0_minus_psi_zr
      external psi0_minus_psi_zr
c-----local
      double precision X1,X2, 
c---------------------------------------------------------
     & psi_0,   ! poloidal flux in /common_psi_m_psi_zr/ 
     & z_0      ! vertical coordinate in 
                !  common /common_psi_minus_psi_zr/ 
c-----------------------------------------------------------
c-----external
      double precision psif         

      common /common_psi_minus_psi_zr/ psi_0,z_0 ! input data

c----------------------------------------------------------
c     put parameters to common / common_psi_minus_psi_zr/
c----------------------------------------------------------
      psi_0=psi0
      z_0=z_0
     
c----------------------------------------------------------
c     calculate maximal root: rmax_psi_z
c----------------------------------------------------------
      if(psi0.gt.psif(z0,rmax*1.001d0)) then
        rmax_psi_z=rmax
      else
        X1=rmax*1.001d0
        X2=rma 
        rmax_psi_z=rtbis(psi0_minus_psi_zr,x1,x2,accuracy)   
      endif

c----------------------------------------------------------
c     calculate minimal root: rmin_psi_z
c----------------------------------------------------------
      if(psi0.gt.psif(z0,rmin*0.999d0)) then
        rmin_psi_z=rmin
      else
        X1=rmin*0.999d0
        X2=rma 
         rmin_psi_z=rtbis(psi0_minus_psi_zr,x1,x2,accuracy)   
      endif

      return
      end


