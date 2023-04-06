c*************************spldens**************************************
c  creation of arrays for density spline coefficients		      *
c  trdens(ndens4a),cxdens(ndens4a)                  		      *
c  creation of arrays for temperature, zeff,tpop,vflow                *
c  spline coefficients	                                              *
c                                                                     *
c  this program uses the following subroutines:  		      *
c								      *
c  iac1r(rhom,ip,ip4,denm,lx,mxa,fxa,mxb,fxb,trdens,cxdens)	      *
c  -calculates the spline coefficients for 1d function		      *
c								      *
c**********************************************************************
c  input data from common 'six.i' and 'one.i'			      *
c  output data  for GENRAY are in files: spldens.dat     	      *
c**********************************************************************

      subroutine spldens1
cSAP090311
      implicit none
c      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'six.i'
      include 'one.i'
 
c-----locals
      real*8 tempearr(ndensa),zeffarr(ndensa)
      real*8 tpoparr(ndensa),vflowarr(ndensa)

cSm061205 to use spline functions coeff1
      real*8 d2_f_drho(ndensa),work(3*ndensa+1),tabl(3),
     &h,fxb,fxa,rhoj,sdens, stemp,svflow,stpop
      integer i1p(2),itabl(3),
     &ip,i,ip4,j,k,lx,mxa,mxb,nbulk1,ndens4

c-----externals
      real*8 ias1r,
     & temperho,tpoprho,vflowrho,densrho

      if(nbulka.lt.nbulk) then
        write(*,*)'nbulka is given in param.i=',nbulka
        write(*,*)'nbulk is given in in genray.in=',nbulk
	write(*,*)'nbulka must be .ge. nbulk,1
     1  please change these parameters'
	stop
      endif
      write(*,*)'in spldens1 izeff',izeff
c--------------------------------------------------------------------
c     creation of the spline coefficients for density
c     temperature and Zeff
c--------------------------------------------------------------------
      ! h=1.d0/(ndens-1)
c for test
c      do i=1,ndens-1
c          x1=(h*(i-1))
c	  densm(i)=(dense0-denseb)*(1.d0-x1**rn1de)**rn2de+denseb
c	  tempearr(i)=(te0-teb)*(1.d0-x1**rn1te)**rn2te+teb
c	  zeffarr(i)=(zeff0-zeffb)*(1.d0-x1**rn1zeff)+zeffb
c      enddo
c      densm(ndensa)=denseb
c      tempearr(ndensa)=teb
c      zeffarr(ndensa)=zeffb
c-----------------------------------------------------------------
c     creation of the data file: dentemz.dat .This file has
c     the density,temperature and Zeff radial profiles
c----------------------------------------------------------------
cc      open(10,file='dentemz.dat')
cc      write(10,3) (densm(i),i=1,ndens)
cc      write(10,3) (tempearr(i),i=1,ndens)
cc      write(10,3) (zeffarr(i),i=1,ndens)
cc      close(10)
c--------------------------------------------------------------------
c     reading of the file dentemz.dat with the density,temperature
c     and Zeff radial profiles
c---------------------------------------------------------------------
cc      open(10,file='dentemz.dat')
cc3     format(5e16.9)
cc      read(10,3) (densm(i),i=1,ndens)
cc      read(10,3) (tempearr(i),i=1,ndens)
cc      read(10,3) (zeffarr(i),i=1,ndens)
cc      close(10)
c---------------- densm in 10**13/ cm**3----------------------------
c     calculations of the spline coefisients for the spline
c     approximation of density, temperature ,Zeff
c     tpop,vflow
c-------------------------------------------------------------------
      ip=ndens
      do i=1,ndens
        !rhom(i)=h*(i-1) ! should be defined by now
        write(*,*) 'spldens1: i, rhom(i)=', i, rhom(i)
      enddo
      !pause !!!
c-------------------------------------------------------------------
c     calculation of nbulk1
      if(izeff.eq.0) then
c        zeff will be calculated using the given ions;
         nbulk1=nbulk
      else
c        =1 zeff is given, the ions component will be calculated
         if (nbulk.eq.1) nbulk1=1
         if (nbulk.eq.2) nbulk1=2
         if (nbulk.gt.2) nbulk1=nbulk-2
      endif
cc      write(*,*)'in spldens nbulk1=',nbulk1
c-------------------------------------------------------------------
     
c     calculation of the spline coefficients for density and temperature
c     tpop and vflow  using old spline
      ip4=ndens+4
cSm030226
      ndens4=ndens+4
      lx=1
      mxa=0
      fxa=0.0d0
c      mxb=0
cSm050202 changed mbx=1 it removed the negative temperature at plasma edge
      mxb=1
      fxb=0.0d0
c------------------------------------------------
c     test dens1(i,k)
      do k=1,nbulk
         write(*,*)'spldens1 k'
	 do i=1,ndens
            write(*,*)'k,i,dens1(i,k)',k,i,dens1(i,k)	    
	 enddo
         do i=1,ndens
            write(*,*)'k,i,temp1(i,k)',k,i,temp1(i,k)	    
	 enddo
         do i=1,ndens
            write(*,*)'k,i,tpop1(i,k)',k,i,tpop1(i,k)	    
	 enddo

         do i=1,ndens
            write(*,*)'k,i,vflow1(i,k)',k,i,vflow1(i,k)	    
	 enddo
      enddo
c     end test dens1(i,k)
c-----------------------------------------------
      
      do k=1,nbulk
	 do i=1,ndens
	    densm(i)=dens1(i,k)
	    tempearr(i)=temp1(i,k)
            tpoparr(i)=tpop1(i,k)
            vflowarr(i)=vflow1(i,k) !m/sec
	 enddo
cSm070103
         goto 10

c-----------------------------------------------------------------
c        old spline coefficients
         call iac1r(rhom,ip,ip4,densm,lx,mxa,fxa,mxb,fxb,trdens,cxdens)
         call iac1r(rhom,ip,ip4,tempearr,lx,mxa,fxa,mxb,fxb,trtempe,
     1           cxtempe)
         call iac1r(rhom,ip,ip4,tpoparr,lx,mxa,fxa,mxb,fxb,trtpop,
     1           cxtpop)
         call iac1r(rhom,ip,ip4,vflowarr,lx,mxa,fxa,mxb,fxb,trvflow,
     1           cxvflow)

c_test old spline
c         write(*,*)'in subroutine spldens1 form ias1r'
c         write(*,*)'k',k
c         idx=0
c         do i=1,ndens
c           dens_test=ias1r(trdens,ndens,ndens4,cxdens,idx,rhom(i))
c           write(*,*)'i,rhom(i),densm(i),dens_test',
c     &             i,rhom(i),densm(i),dens_test
c         enddo
c         imax=101
c         h_rho=1.d0/dfloat(imax-1)
c         do i=1,imax-1
c           rho=h_rho*i
c           dens_test=ias1r(trdens,ndens,ndens4,cxdens,idx,rho)
c           write(*,*)'i,rho,dens_test',i,rho,dens_test
c         enddo
c_endtest old spline

         do i=1,ndens4
	   trdens1(i,k)=trdens(i)
	   trtemp1(i,k)=trtempe(i)
	   cxdens1(i,k)=cxdens(i)
	   cxtemp1(i,k)=cxtempe(i)
           trtpop1(i,k)=trtpop(i)
	   trvflow1(i,k)=trvflow(i)
	   cxtpop1(i,k)=cxtpop(i)
	   cxvflow1(i,k)=cxvflow(i)
         enddo

c-----------------------------------------------------------------
cSm061205 spline coefficients using call coeff1 for
c         density,temperature,tpop,vflow
c----------------------------------------------------------------
 10     continue

c        i1p(1)=4
c        i1p(2)=1 !second derivative given at x(ndens).  place
c                 !value of second derivative in w(1)
c                 !before call to coeff1.
c        d2_f_drho(ndens)=0.d0
c---------------------------------------------------------------
c         i1p(2)=2 !first derivative given at x(ndens).  place
c                 !value of second derivative in w(1)
c                 !before call to coeff1.
c         d2_f_drho(ndens)=(densm(ndens)-densm(ndens-1))/
c     &                    (rhom(ndens)-rhom(ndens-1))
c--------------------------------------------------------------
        i1p(2)=4 ! boundary condition at rho=1
                 ! the first derivative at rho(ndens) is
                 !  calculated by fitting a cubic to 4 points
                 !  rho(nneds) through rho(ndens-3).
cSm070102
c------- boundary condition at rho=0: df/drho=0
        i1p(1)=2 
        d2_f_drho(1)=0.d0
c--------------------------------------------------------
        call coeff1(ndens,rhom,densm,d2_f_drho,i1p,1,work) !density
        do i=1,ndens
	  d2_dens_drho1(i,k)=d2_f_drho(i)
        enddo

        d2_f_drho(1)=0.d0
        call coeff1(ndens,rhom,tempearr,d2_f_drho,i1p,1,work) !temperature
        do i=1,ndens
	  d2_temp_drho1(i,k)=d2_f_drho(i)
        enddo    

        d2_f_drho(1)=0.d0
        call coeff1(ndens,rhom,tpoparr,d2_f_drho,i1p,1,work) !tpop
        do i=1,ndens
	  d2_tpop_drho1(i,k)=d2_f_drho(i)
        enddo

        d2_f_drho(1)=0.d0
        call coeff1(ndens,rhom,vflowarr,d2_f_drho,i1p,1,work) !vflow
        do i=1,ndens
	  d2_vflow_drho1(i,k)=d2_f_drho(i)
        enddo
c------------------------------------------------------------------
      enddo ! k

c----------------------------------------------------------------------
c-------------------------------------------------------------------
c     calculation of the spline coefficients for zeff
      do i=1,ndens
         zeffarr(i)=zeff1(1)
      enddo

cSm061206
      goto 20
c-----old spline coefficients calculations for zeff  
      call iac1r(rhom,ip,ip4,zeffarr,lx,mxa,fxa,mxb,fxb,trzeff,cxzeff)

cSm061206
 20   continue
c-----new spline coefficient calculations using coeff1 for zeff
      d2_f_drho(1)=0.d0
      call coeff1(ndens,rhom,zeffarr,d2_zeff_drho1,i1p,1,work) 

      return
      end


      real*8 FUNCTION temperho(rho_small,i)
c-------------------------------------------------------
c     temperature on radius (from spline)
c     i=species (1,electrons)
c-------------------------------------------------------
      implicit none
c      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'

cSm030226
      include 'one.i'
      include 'six.i'
cSAP090204
      include 'edge_prof_nml.i'

c-----input
      real*8 rho_small ! small radius normalized
      integer i        ! number of plasma species
c-----locals  
cSm061206
      real*8 tabl(3),d2_temp_drho(ndensa),temp(ndensa), 
     &rhoedg,temperho1,tempedg, rho095,temp095, dln
      integer itabl(3),idx,j,k,ndens4

c-----externals
cSAP090304
c      double precision ias1r
      real*8 ias1r


      idx=0
cSm030226
      ndens4=ndens+4
      do j=1,ndens4
	 trtempe(j)=trtemp1(j,i)
	 cxtempe(j)=cxtemp1(j,i)
      enddo
      if (rho_small.gt.(1.d0-1.d-5)) then
c        the point is outside the plasma
         rhoedg=1.d0
         goto 10
cold         tempedg=ias1r(trtempe,ndens,ndens4,cxtempe,idx,rhoedg) !old spline
cSm061206
c-------------------------------------------------------------
 10      continue               !new spline using coef1, terp1       
         itabl(1)=1
         itabl(2)=0
         itabl(3)=0 
         do k=1,ndens
            temp(k)=temp1(k,i)
            d2_temp_drho(k)=d2_temp_drho1(k,i)
         enddo
         call terp1(ndens,rhom,temp,d2_temp_drho,rhoedg,1,tabl,
     &              itabl)
         tempedg=tabl(1)
         ! YuP-120412: redefine sigmedgt to make a smooth derivative 
         ! in profile at rho=1.
         rho095=0.99
         call terp1(ndens,rhom,temp,d2_temp_drho,rho095,1,tabl,
     &              itabl)
         temp095=tabl(1)
         dln=log(temp095/tempedg) !should be positive for most profiles
         if(dln.gt.0.d0)then
            sigmedgt=(1.0-rho095)/dln
         else ! dln=0., which means a flat profile near rho~1.
            sigmedgt=2.0 ! a very slow drop after rho>1
         endif
c--------------------------------------------------------------
         ! temp: Exponential drop at rho>rho095
         temperho1=tempedg*dexp(-(rho_small-1.d0)/sigmedgt)
c         temperho=dmax1(temperho1,tempedg*1.d-1)
         temperho=dmax1(temperho1,temp_min_edge)
c         write(*,*)'temperho  rho_small, tempedg, temperho',
c     +                        rho_small,tempedg,temperho
      else 
         goto 20       
c old         temperho=ias1r(trtempe,ndens,ndens4,cxtempe,idx,rho_small)!old spline
cSm061206 new spline: coeff1 terp1      
c-------------------------------------------------------------
 20      continue                
         itabl(1)=1
         itabl(2)=0
         itabl(3)=0 
         do k=1,ndens
            temp(k)=temp1(k,i)
            d2_temp_drho(k)=d2_temp_drho1(k,i)
         enddo
         call terp1(ndens,rhom,temp,d2_temp_drho,rho_small,1,tabl,
     &              itabl)
         temperho=tabl(1)
c         write(*,*)'  temperho',temperho
c--------------------------------------------------------------
      endif

      return

      END


      real*8 FUNCTION zeffrho(rho_small)
c-------------------------------------------------------
c     Zeff on radius (from spline)
c-------------------------------------------------------
      implicit none
c      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
cSm030226
      include 'one.i'      
      include 'six.i'
c-----input
      real*8 rho_small ! small radius normalized
c-----locals
cSm061206
      real*8 tabl(3),rhoedg
      integer itabl(3),ndens4,idx
c-----externsals
      real*8 ias1r

cSm030226
      ndens4=ndens+4
      idx=0
      if (rho_small.gt.(1.d0-1.d-5)) then
c        the point is outside the plasma
         rhoedg=1.d0
cSm061206
         goto 10
         zeffrho=ias1r(trzeff,ndens,ndens4,cxzeff,idx,rhoedg) !old spline

cSm061206-------------------------------------------------------------
 10      continue !new spline : coeff1 terp1
         itabl(1)=1
         itabl(2)=0
         itabl(3)=0          
         call terp1(ndens,rhom,zeff1,d2_zeff_drho1,rhoedg,1,tabl,itabl)
         zeffrho=tabl(1)
c----------------------------------------------------------------------

c	 write(*,*)'in edg zeffrho zeffrho',zeffrho
      else
cSm061206
         goto 20
         zeffrho=ias1r(trzeff,ndens,ndens4,cxzeff,idx,rho_small) !old spline

cSm061206-------------------------------------------------------------------
 20      continue !new spline : coeff1 terp1
         itabl(1)=1
         itabl(2)=0
         itabl(3)=0          
         call terp1(ndens,rhom,zeff1,d2_zeff_drho1,rho_small,1,tabl,
     &   itabl)
         zeffrho=tabl(1)
c-------------------------------------------------------------------------

      endif
      return
      END


      real*8 FUNCTION densrho(rho_small,i)
c-------------------------------------------------------
c     density on radius (from spline)
c-------------------------------------------------------
      implicit none
c      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'     
      include 'six.i'
      include 'edge_prof_nml.i'
c-----input
      real*8 rho_small ! small radius normalized
      integer i        ! number of plasma species
      real*8 tabl(3),d2_dens_drho(ndensa),
     & rhoedg,densedg,
     & densedge_e,dens_ratio, dens095, rho095, dln

      integer itabl(3),
     & ndens4,idx,j,k

c-----externals
      real*8 ias1r
      
      if(model_rho_dens.ne.0) then
        stop 'densrho() is only valid for model_rho_dens=0'
      endif

      ndens4=ndens+4
      idx=0
c      write(*,*)'in densrho(rho_small,i), rho_small,i',rho_small,i
c      write(*,*)'in densrho ndens4',ndens4
      do j=1,ndens4
	 trdens(j)=trdens1(j,i)
	 cxdens(j)=cxdens1(j,i)
      enddo
      if (rho_small.gt.(1.d0-1.d-5)) then
c        the point is outside the plasma
         rhoedg=1.d0
         goto 10
cyup         densedg=ias1r(trdens,ndens,ndens4,cxdens,idx,rhoedg) !old spline
c         write(*,*)'ias1 rhoedg,densedg', rhoedg,densedg
cSm061204 new spline using coeff1 terp1-----------------------
 10      continue
         itabl(1)=1
         itabl(2)=0
         itabl(3)=0 
         do k=1,ndens
            densm(k)=dens1(k,i)
            d2_dens_drho(k)=d2_dens_drho1(k,i)
         enddo
         call terp1(ndens,rhom,densm,d2_dens_drho,rhoedg,1,tabl,
     &              itabl)
         densedg=tabl(1)  
c--------------------------------------------------------------------------
         ! YuP-120412: redefine sigmedgn to make a smooth derivative 
         ! in profile at rho=1.
         rho095=0.99
         call terp1(ndens,rhom,densm,d2_dens_drho,rho095,1,tabl,
     &              itabl)
         dens095=tabl(1)
         dln=log(dens095/densedg) !should be positive for most profiles
         if(dln.gt.0.d0)then
            sigmedgn=(1.0-rho095)/dln
         else ! dln=0., which means a flat profile near rho~1.
            sigmedgn=2.0 ! a very slow drop after rho>1
         endif
         ! dens: Exponential drop at rho>rho095
         densrho=densedg*dexp(-(rho_small-1.d0)/sigmedgn)
         !write(*,*)'densrho: sigmedgn,sigmedgt=',sigmedgn,sigmedgt
         !pause
         if (i.gt.1) then
c----------------------------------------------------------
c           calculates the ratio density(i)/density_e at rho=1
c-----------------------------------------------------------
c           calculate the electron density at rho=1
            do k=1,ndens
               densm(k)=dens1(k,1)
               d2_dens_drho(k)=d2_dens_drho1(k,1)
            enddo
            call terp1(ndens,rhom,densm,d2_dens_drho,rhoedg,1,tabl,
     &               itabl)
            densedge_e=tabl(1) 
            dens_ratio=densedg/densedge_e
         endif
         densrho=dmax1(densrho,dens_min_edge*dens_ratio)

      else ! in plasma (rho<1)

         goto 20
cyup         densrho=ias1r(trdens,ndens,ndens4,cxdens,idx,rho_small) !old spline
cyup         write(*,*)'in densrho ias1r: rho_small densrho',rho_small,
cyup     &   densrho
cSm061204 new spline: coeff1 terp1
c----------------------------------------------------------------
 20      continue
         itabl(1)=1
         itabl(2)=0
         itabl(3)=0 
         do k=1,ndens
            densm(k)=dens1(k,i)
            d2_dens_drho(k)=d2_dens_drho1(k,i)
         enddo
         call terp1(ndens,rhom,densm,d2_dens_drho,rho_small,1,tabl,
     &              itabl)
         densrho=tabl(1)
c------------------------------------------------------------------
c	 write(*,*)'in densrho terp1: rho_small,densrho',rho_small,
c     &   densrho
      endif
      return
      END





      real*8 FUNCTION tpoprho(rho_small,i)
c-------------------------------------------------------
c     tpop=T_perp/T_parallel on radius (from spline),i is a tipe of the plasma component
c-------------------------------------------------------
cSAP090311
      implicit none
c      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
cSm030226
      include 'one.i'
      include 'six.i'

c-----input
      real*8 rho_small !normalized small radius
      integer i        !number of plasma specie
c-----locals
cSm061205
      real*8 tabl(3),d2_tpop_drho(ndensa),tpop(ndensa),rhoedg
      integer itabl(3),ndens4,idx,j,k

c-----externals
      real*8 ias1r
cSm030226
      ndens4=ndens+4
      idx=0

cSAP090311
c      write(*,*)'spldens.f tpoprho i, rho_small',i,rho_small 

      do j=1,ndens4
	 trtpop(j)=trtpop1(j,i)
	 cxtpop(j)=cxtpop1(j,i)
      enddo
      if (rho_small.gt.(1.d0-1.d-5)) then
c        the point is outside the plasma
cSAP080831
	 rhoedg=1.d0-1.d-10
	 rhoedg=1.d0
         goto 10
         tpoprho=ias1r(trtpop,ndens,ndens4,cxtpop,idx,rhoedg) !old spline

cSm061204 
c--------new spline using coeff1 terp1-----------------------
 10      continue
         itabl(1)=1
         itabl(2)=0
         itabl(3)=0

cSAP090311
c         write(*,*)'spldens.f tpoprho i, rho_small',i,rho_small 
cSAP090311
         do k=1,ndensa
            tpop(k)= 1.d0
            d2_tpop_drho(k)=0.d0
         enddo

         do k=1,ndens
            tpop(k)=tpop1(k,i)
            d2_tpop_drho(k)=d2_tpop_drho1(k,i)
         enddo

c        write(*,*)'tpop',tpop
c        write(*,*)'d2_tpop_drho',d2_tpop_drho
c        write(*,*)'rhom',rhom
  
         call terp1(ndens,rhom,tpop,d2_tpop_drho,rhoedg,1,tabl,
     &              itabl)
         tpoprho=tabl(1)

c         write(*,*)'spldens.f in tpoprho rhoedg,tpoprho',rhoedg,tpoprho

c--------------------------------------------------------------------------

      else
         goto 20
         tpoprho=ias1r(trtpop,ndens,ndens4,cxtpop,idx,rho_small) !old spline

cSm061206 new spline: coeff1 terp1      
c-------------------------------------------------------------
 20      continue                
         itabl(1)=1
         itabl(2)=0
         itabl(3)=0 
         do k=1,ndens
            tpop(k)=tpop1(k,i)
            d2_tpop_drho(k)=d2_tpop_drho1(k,i)
         enddo
         call terp1(ndens,rhom,tpop,d2_tpop_drho,rho_small,1,tabl,
     &              itabl)
         tpoprho=tabl(1)
c--------------------------------------------------------------
      endif
      return
      END

      real*8 FUNCTION vflowrho(rho_small,i) 
c-------------------------------------------------------
c     vflow(cm/sec) on radius (from spline),
c     i is type of the plasma component
c-------------------------------------------------------
      implicit none
c      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
cSm030226      
      include 'one.i'
      include 'six.i'
c-----input
      real*8 rho_small ! small radius normalized
      integer i        ! number of plasma species
c-----locals
cSm061205
      real*8 tabl(3),d2_vflow_drho(ndensa),vflow(ndensa),
     &rhoedg
      integer itabl(3),
     &ndens4,idx,j,k

c-----externals
      real*8 ias1r

cSm030226
      ndens4=ndens+4

      idx=0

      do j=1,ndens4
	 trvflow(j)=trvflow1(j,i)
	 cxvflow(j)=cxvflow1(j,i)
      enddo
      if (rho_small.gt.(1.d0-1.d-5)) then
c        the point is outside the plasma
	 rhoedg=1.d0
         goto 10
         vflowrho=ias1r(trvflow,ndens,ndens4,cxvflow,idx,rhoedg)	!m/sec
cSm061204 
c--------new spline using coeff1 terp1-----------------------
 10      continue
         itabl(1)=1
         itabl(2)=0
         itabl(3)=0 
         do k=1,ndens
            vflow(k)=vflow1(k,i)                 !m/sec
            d2_vflow_drho(k)=d2_vflow_drho1(k,i)
         enddo
         call terp1(ndens,rhom,vflow,d2_vflow_drho,rhoedg,1,tabl,
     &              itabl)
         vflowrho=tabl(1)
c--------------------------------------------------------------------------
      else
         goto 20
         vflowrho=ias1r(trvflow,ndens,ndens4,cxvflow,idx,rho_small) !   m/sec old soline
cSm061204 
c--------new spline using coeff1 terp1-----------------------
 20      continue
         itabl(1)=1
         itabl(2)=0
         itabl(3)=0 
         do k=1,ndens
            vflow(k)=vflow1(k,i)               !m/sec
            d2_vflow_drho(k)=d2_vflow_drho1(k,i)
         enddo
         call terp1(ndens,rhom,vflow,d2_vflow_drho,rho_small,1,tabl,!  m/sec
     &              itabl)
         vflowrho=tabl(1)
c--------------------------------------------------------------------------
      endif

      vflowrho=1.d2*vflowrho !sm/sec

      return
      END




      real*8 FUNCTION ddnsrho(rho_small,i)
c-------------------------------------------------------
c     derivative density by radius (from spline)
c-------------------------------------------------------
      implicit none
c      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'six.i'
      include 'edge_prof_nml.i'
c-----input
      real*8 rho_small !small radius, normalized
      integer i        !number of plasma specie
c-----locals
      real*8 tabl(3),d2_dens_drho(ndensa),
     & densedg,densrho,rhoedg,
     & densedge_e ,dens_ratio
      integer itabl(3),ndens4,j,idx,k
c-----externals
      real*8 ias1r

      if(ixyz.eq.1) then
       ! stop 'ddnsrho:  Not setup for ixyz=1'
      endif
      
      if(model_rho_dens.ne.0) then
        stop 'ddnsrho: Only valid for model_rho_dens=0'
      endif
      
      ndens4=ndens+4
      do j=1,ndens4
	 trdens(j)=trdens1(j,i)
	 cxdens(j)=cxdens1(j,i)
      enddo

      if (rho_small.gt.(1.d0-1.d-5)) then
         idx=0
c        the point is outside the plasma
         rhoedg=1.d0
         goto 10
         densedg=ias1r(trdens,ndens,ndens4,cxdens,idx,rhoedg) !old spline

cSm061204-------new soline coeff1, terp1
 10       continue
         itabl(1)=1
         itabl(2)=0
         itabl(3)=0 
         do k=1,ndens
            densm(k)=dens1(k,i)
            d2_dens_drho(k)=d2_dens_drho1(k,i)
         enddo
         call terp1(ndens,rhom,densm,d2_dens_drho,rhoedg,1,tabl,
     &              itabl)
         densedg=tabl(1)
c---------------------------------------------------------------------
         densrho=densedg*dexp(-(rho_small-1.d0)/sigmedgn)
         dens_ratio=1.d0
         if (i.gt.1) then
c----------------------------------------------------------
c           calculates the ratio density(i)/density_e at rho=1
c-----------------------------------------------------------
c           calculate the electron density at rho=1
            do k=1,ndens
               densm(k)=dens1(k,1)
               d2_dens_drho(k)=d2_dens_drho1(k,1)
            enddo
            call terp1(ndens,rhom,densm,d2_dens_drho,rhoedg,1,tabl,
     &                 itabl)
            densedge_e=tabl(1) 
            dens_ratio=densedg/densedge_e
         endif

         if(densrho.ge.dens_min_edge*dens_ratio) then
           ddnsrho=-densrho/sigmedgn
         else
           densrho=dens_min_edge*dens_ratio
           ddnsrho=0.d0
         endif

      else ! rho_small<1
         idx=1

         goto 20
         ddnsrho=ias1r(trdens,ndens,ndens4,cxdens,idx,rho_small) !old spline

cSm061204-------new soline coeff1, terp1
 20       continue
         itabl(1)=0
         itabl(2)=1
         itabl(3)=0 
         do k=1,ndens
            densm(k)=dens1(k,i)
            d2_dens_drho(k)=d2_dens_drho1(k,i)
         enddo
 
         call terp1(ndens,rhom,densm,d2_dens_drho,rho_small,1,tabl,
     &              itabl)
         ddnsrho=tabl(2)
c----------------------------------------------------------------
      endif
      return
      END




      real*8 FUNCTION dtemdrho(rho_small,i)
c-------------------------------------------------------
c     derivative: temperature by radius (from spline) T-Kev
c-------------------------------------------------------
      implicit none
c      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
cSm030226      
      include 'one.i'
      include 'six.i'
cSAP090304
      include 'edge_prof_nml.i'
c-----input
      real*8 rho_small !small radius, normalized
      integer i        !number of plasma specie
c-----locals
cSm061206
      real*8 tabl(3),d2_temp_drho(ndensa),temp(ndensa),
     &rhoedg,temprho1,temprho,temperho,tempedg
      integer itabl(3),ndens4,idx,j,k
c-----externals
      real*8 ias1r
      
cSm030226
      ndens4=ndens+4

      do j=1,ndens4
	 trtempe(j)=trtemp1(j,i)
	 cxtempe(j)=cxtemp1(j,i)
      enddo
      if (rho_small.gt.(1.d0-1.d-5)) then
         idx=0
c        the point is outside the plasma
	 rhoedg=1.d0
cBH080731	 sigmedg=0.02d0  !Using sigmedgt
         goto 10
         tempedg=ias1r(trtempe,ndens,ndens4,cxtempe,idx,rhoedg) !old spline

cSm061206
c-------------------------------------------------------------
 10      continue               !new spline using coef1, terp1       
         itabl(1)=1
         itabl(2)=0
         itabl(3)=0 
         do k=1,ndens
            temp(k)=temp1(k,i)
            d2_temp_drho(k)=d2_temp_drho1(k,i)
         enddo
         call terp1(ndens,rhom,temp,d2_temp_drho,rhoedg,1,tabl,
     &              itabl)
         tempedg=tabl(1)
c--------------------------------------------------------------
	 temprho1=tempedg*dexp(-(rho_small-1.d0)/sigmedgt)

cSAP090304
c         if (temprho1.ge.tempedg*1.d-1) then
         if (temprho1.ge.temp_min_edge) then
            temprho=temprho1
            dtemdrho=-temprho/sigmedgt
         else
c            temperho=tempedg*1.d-1
            temperho=temp_min_edge 
            dtemdrho=0.d0
         endif         
      else
         idx=1

         goto 20
         dtemdrho=ias1r(trtempe,ndens,ndens4,cxtempe,idx,rho_small)!old spline

cSm061206 new spline: coeff1 terp1      
c-------------------------------------------------------------
 20      continue                
         itabl(1)=0
         itabl(2)=1
         itabl(3)=0 
         do k=1,ndens
            temp(k)=temp1(k,i)
            d2_temp_drho(k)=d2_temp_drho1(k,i)
         enddo
         call terp1(ndens,rhom,temp,d2_temp_drho,rho_small,1,tabl,
     &              itabl)
         dtemdrho=tabl(2)
c--------------------------------------------------------------
      endif
      return
      END

      real*8 FUNCTION dtpoprho(rho_small,i)
c-------------------------------------------------------
c     derivative tpop=T_perp/T_parallel by radius (from spline)
c-------------------------------------------------------
      implicit none
c      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
cSm030226
      include 'one.i'
      include 'six.i'

c-----input
      real*8 rho_small !normalized small radius
      integer i        !number of plasma specie
c-----locals
cSm061205
      real*8 tabl(3),d2_tpop_drho(ndensa),tpop(ndensa),rhoedg,tpoprho
      integer itabl(3),ndens4,idx,j,k

c-----externals
      real*8 ias1r

cSm030226
      ndens4=ndens+4
      
      do j=1,ndens4
	 trtpop(j)=trtpop1(j,i)
	 cxtpop(j)=cxtpop1(j,i)
      enddo
      idx=1
      if (rho_small.gt.(1.d0-1.d-5)) then
c        the point is outside the plasma
         dtpoprho=0.d0	  
      else

         goto 20
         dtpoprho=ias1r(trtpop,ndens,ndens4,cxtpop,idx,rho_small)

cSm061206 new spline: coeff1 terp1      
c-------------------------------------------------------------
 20      continue                
         itabl(1)=0
         itabl(2)=1
         itabl(3)=0 
         do k=1,ndens
            tpop(k)=tpop1(k,i)
            d2_tpop_drho(k)=d2_tpop_drho1(k,i)
         enddo
         call terp1(ndens,rhom,tpop,d2_tpop_drho,rho_small,1,tabl,
     &              itabl)
         tpoprho=tabl(2)
c--------------------------------------------------------------
      endif
      return
      END

      real*8 FUNCTION dvflwrho(rho_small,i)
c-------------------------------------------------------
c     derivative vflow by radius (from spline) d(cm/sec)/d_rho
c-------------------------------------------------------
      implicit none
c      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
cSm030226
      include 'one.i'
      include 'six.i'
c-----input
      real*8 rho_small ! small radius normalized
      integer i        ! number of plasma species
c-----locals
cSm061205
      real*8 tabl(3),d2_vflow_drho(ndensa),vflow(ndensa)
      integer itabl(3),k
c-----externals
      real*8 ias1r
cSm030226
c      ndens4=ndens+4
c      do j=1,ndens4
c         trvflow(j)=trvflow1(j,i)
c	 cxvflow(j)=cxvflow1(j,i)
c      enddo
c      idx=1
      
      if (rho_small.gt.(1.d0-1.d-5)) then
c        the point is outside the plasma
         dvflwrho=0.d0	  
      else
         goto 20
c         dvflwrho=ias1r(trvflow,ndens,ndens4,cxvflow,idx,rho_small)

cSm061204 
c--------new spline using coeff1 terp1-----------------------
 20      continue
         itabl(1)=0
         itabl(2)=1
         itabl(3)=0 
         do k=1,ndens
            vflow(k)=vflow1(k,i)
            d2_vflow_drho(k)=d2_vflow_drho1(k,i)
         enddo
         
         call terp1(ndens,rhom,vflow,d2_vflow_drho,rho_small,1,tabl,!  m/sec
     &              itabl)
         dvflwrho=tabl(2)
c--------------------------------------------------------------------------\       
      endif

      dvflwrho=1.d2*dvflwrho !d(sm/sec)/drho
c      write(*,*)'dvflwrho', dvflwrho
      return
      END



c        ************************ rhof(z,r,phi)  **************
c        *                        -                           *
c        * this function calculates small radius rho          *
c        * it controls if the ray point is inside the plasma  *
c        * if the ray point is outside the plasma then the    *
c        * it chooses the specific value for rho	      *
c        ******************************************************
c
c-------------------------------------------------------------
c        input parameters				     !
c      							     !
c      z, r, phi    
c------------------------------------------------------------
      real*8
     1function rhof(z,r,phi) ! for GENRAY-C, check rhof_xyz
      implicit none
c      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
      include 'six.i'
c-----input
      real*8 z,r,phi

c-----externals
      real*8
     1 ias1r,ias2r_Sm,b,rhopsi,
     & thetapol
c-----locals
      real*8 epsbnd,rrr,zzrp,zp,zzrm,zm,res,delr,delz,
     &theta_l,r_b,z_b,
     &x_l,y_l,l_zr,l_z_b_r_b
      integer idx,ipx,ipx4,nr4,nz4

      bmod=b(z,r,phi)     
c--------------------------------------------------------------
c     control that the point is inside the plasma
      epsbnd=1.d-10
      iboundb=-1

c      write(*,*)'spldens.g in rhof r,z',r,z
c-----------------------
cSAP080817
      x_l=r-rma
      y_l=z-zma

      if ((x_l**2+y_l**2).lt.1.d-10) then
         goto 10
      endif

c-----calculate poloidal angle  0=<theta_l<2pi

c      write(*,*)'spldens.f rhof(z,r,phi) r,z,rma,zma', r,z,rma,zma

      call theta_rz(x_l,y_l,theta_l) 

c      write(*,*)'spldens.f rhof(z,r,phi)x_l,y_l,theta_l ',
c     &x_l,y_l,theta_ l

c      write(*,*)'spldens.g in theta_l', theta_l

c-----calculate coordinates x_b,z_b at the limmiter surface psi=psilim
c     and the poloidal angle thetapol
      call zr_psith(psilim,theta_l,z_b,r_b)

c      write(*,*)'spldens.g in psilim,z_b,r_b',psilim,z_b,r_b

      l_zr=dsqrt(x_l**2+y_l**2)
      l_z_b_r_b=dsqrt((r_b-rma)**2+(z_b-zma)**2)
      if (l_zr.ge.l_z_b_r_b) then
         rho = l_zr / l_z_b_r_b
         rhof=rho
         iboundb=3
      endif
      goto 10
c-----------------------
      if ((r.le.rmin+epsbnd).or.(r.ge.rmax-epsbnd)) then
c         write(*,*)'in rhof rmin+epsbnd,r,rmax-epsbnd',
c     1   rmin+epsbnd,r,rmax-epsbnd
         iboundb=2
c	 delr=0.5d0*(rmax-rmin)
c         if (r.le.rmin+epsbnd) rhof=1.d0+(rmin-r+epsbnd)/delr
c         if (r.ge.rmax-epsbnd) rhof=1.d0+(r-rmax+epsbnd)/delr
c         write(*,*)'rhof delr,rhof',delr,rhof
c	 respsi=psimag+(psilim-psimag)*rhof*rhof
cSm070325 
c         write(*,*)'b.f rho',rho
         theta_l=thetapol(z,r) 
         call zr_psith(psilim,theta_l,z_b,r_b)
         rhof=dsqrt(((r-rma)**2+(z-zma)**2)/((r_b-rma)**2+(z_b-zma)**2))
         goto 10
      end if
c---------------- idx derivativs order 0.ge.idx.le.3---------------
      rrr=r
      idx=0
      ipx=ip
      ipx4=ip+4
      zzrp=ias1r(trlimp,ipx,ipx4,cxlimp,idx,rrr)
      zp=zzrp
      ipx=im
      ipx4=im+4
      zzrm=ias1r(trlimm,ipx,ipx4,cxlimm,idx,rrr)
      zm=zzrm
      if ((z.ge.zp-epsbnd).or.(z.le.zm+epsbnd)) then
c	 delz=0.5d0*(zmax-zmin)
c         if (z.ge.zp-epsbnd) then
c	    rhof=1.d0+(z-zp+epsbnd)/delz
c	 else
c            rhof=1.d0+(zm-z+epsbnd)/delz
c	 endif
         iboundb=1
cSm070324 
c         write(*,*)'b.f rho',rho
         theta_l=thetapol(z,r) 
         call zr_psith(psilim,theta_l,z_b,r_b)
         rhof=dsqrt(((r-rma)**2+(z-zma)**2)/((r_b-rma)**2+(z_b-zma)**2))
       
         goto 10
      end if
c---------------------------------------------------------------------
 10     continue
c end of the control that the point is inside the plasma
c--------------------------------------------------------------
c
c nr4a , nz4a and nrya  are given in param.i 
c nreqda,nzeqda,nlimit are given in param.i
c
cSm030224 
      nr4=nx+4
      nz4=ny+4
 
      ncx=nr4
      ncy=nz4

cSm030224
c      res=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,0,0,r,z)
      res=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,0,0,r,z,nr4a)
c------------------------------------------------------------------
c this  part of the programm calculates rhof (if the point is inside the plasma).


      if (iboundb.eq.-1) rhof=rhopsi(res)

      rho=rhof ! the data for common block one.i
                          
      return
      end




c     ********************** drhodr **********************
c     derivative d(rho)/dr
c------------------------------------------------------------------
c								   !
c     input parameters   					   !
c      								   !
c     z, r phi -  the coordinates of the point where               !
c                 the derivative will be calculated.               !
c------------------------------------------------------------------
c     the input value rho was calculated in b, rho is in one.i 
c     uses
c             rho,iboundb,dpdrd   from common 'one.i'
c     dpdrd was calucated in b ( in one.i)
c     b(z,r,phi) should be run before this function
c------------------------------------------------------------------
c     if the point is outside the plasma 
c     then use a special calculations of drhodr
c------------------------------------------------------------------            
c     functions psif,drhopsi
c------------------------------------------------------------------
      real*8 FUNCTION drhodr(z,r,phi)
      implicit none
c      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'

c-----input
      real*8 z,r,phi !space coordinates

c-----locals
cBH090803      Some compilers don't like declaration of drhodr above and
cMH090803      also here.
cBH090803      real*8 theta_l,rho_l,drhodr,psi_xr,dro_dpsi,z_b,r_b,
      real*8 theta_l,rho_l,psi_xr,dro_dpsi,z_b,r_b,
     &rho_LCFS,rho_non_normalized,
     &x_l,y_l

      real*8 step,rho_l_p,rho_l_m,numer_drhodr,
     &z_b_p,r_b_p,z_b_m,r_b_m,
     &d_z_b_dr,d_r_b_dr,d_rho_LCFS_dr


c-----externals
      real*8 thetapol,psif,drhopsi
      
c      if(ixyz.eq.1) then
c        write(*,*)'drhodr: should not be called if ixyz=1'
c        stop
c      endif

      if(iboundb.ge.1)then
c--------the point is outside LCFS 
         
        theta_l=thetapol(z,r) ! -pi<theta_l<pi 
ctest
c        write(*,*)'in drhodr 1 theta_l',theta_l
c        x_l=r-rma
c        y_l=z-zma
c        call theta_rz(x_l,y_l,theta_l) !0 =< theta_l < 2*pi
c        write(*,*)'in drhodr 2 theta_l',theta_l
c-end_test

        call zr_psith(psilim,theta_l,z_b,r_b)
        
cSAP090328       
c        rho_l=dsqrt(((r-rma)**2+(z-zma)**2)/((r_b-rma)**2+(z_b-zma)**2))
c        drhodr=(r-rma)/(((r_b-rma)**2+(z_b-zma)**2)*rho_l)
        rho_LCFS=dsqrt((r_b-rma)**2+(z_b-zma)**2)
        rho_non_normalized=dsqrt((r-rma)**2+(z-zma)**2)
        rho_l=rho_non_normalized/rho_LCFS

c       drhodr=(r-rma)/(((r_b-rma)**2+(z_b-zma)**2)*rho_l) 
c       Terms with d_r_b/dr and d_r_b/dr should be added
c       
c-------Calculate numerical derivatives d/dr from r_b, z_b
        x_l=r-rma
        y_l=z-zma
        step=1.d-4
        call theta_rz(x_l+step,y_l,theta_l)
        call zr_psith(psilim,theta_l,z_b_p,r_b_p)
        call theta_rz(x_l-step,y_l,theta_l)
        call zr_psith(psilim,theta_l,z_b_m,r_b_m)
        d_z_b_dr=(z_b_p-z_b_m)/(2.d0*step)     
        d_r_b_dr=(r_b_p-r_b_m)/(2.d0*step)
        d_rho_LCFS_dr=((r_b-rma)*d_r_b_dr+(z_b-zma)*d_z_b_dr)/rho_LCFS
ctest----------
c        write(*,*)'d_z_b_dr,d_r_b_dr,d_rho_LCFS_dr',
c     &             d_z_b_dr,d_r_b_dr,d_rho_LCFS_dr
c        drhodr=(r-rma)/(rho_non_normalized*rho_LCFS)
c        write(*,*)'in drhodr old drhodr',drhodr
c-end----- test 
        drhodr=(r-rma)/(rho_non_normalized*rho_LCFS)-
     &  (rho_l/rho_LCFS)*d_rho_LCFS_dr

c        write(*,*)'in drhodr new drhodr',drhodr
c-------test numerical derivative drhodr
c        step=1.d-4
c-------r+step
c        call theta_rz(x_l+step,y_l,theta_l)
c        call zr_psith(psilim,theta_l,z_b,r_b)
c        rho_LCFS=dsqrt((r_b-rma)**2+(z_b-zma)**2)        
c        rho_non_normalized=dsqrt((r+step-rma)**2+(z-zma)**2) 
c        rho_l_p=rho_non_normalized/rho_LCFS

c-------r-step 
c        call theta_rz(x_l-step,y_l,theta_l)
c        call zr_psith(psilim,theta_l,z_b,r_b)
c        rho_LCFS=dsqrt((r_b-rma)**2+(z_b-zma)**2)
c        rho_non_normalized=dsqrt((r-step-rma)**2+(z-zma)**2) 
c        rho_l_m=rho_non_normalized/rho_LCFS
        
c        numer_drhodr=(rho_l_p-rho_l_m)/(2.d0*step) !numerical derivativ drho_dr

c        write(*,*)'in drhodr drhodr,numer_drhodr',drhodr,numer_drhodr
      else
        psi_xr=psif(z,r)
        dro_dpsi=drhopsi(psi_xr)
        drhodr=dro_dpsi*dpdrd    !dpdrd was calucated in b ( in one.i)
      endif

      return
      END


c     ********************** drhodz **********************
c     derivative d(rho)/dz

c------------------------------------------------------------------
c								   !
c        input parameters					   !
c      								   !
c      z, r phi - coordinates of the point where the derivative is !
c                 calculated.        		         	   !
c------------------------------------------------------------------
c     the input value rho was calculated in b, rho is in one.i 
c     uses
c             rho,iboundb,dpdzd   from common 'one.i'
c     dpdrd was calculated in b ( in one.i)
c     b(z,r,phi) should be run before this function
c------------------------------------------------------------------
c         uses
c           rho,iboundb,dpdzd   from common 'one'
c           functions psif,drhopsi
c----------------------------------------------------------------------

      real*8 FUNCTION drhodz(z,r,phi)
      implicit none
c      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'five.i'
      include 'three.i'

c-----input
      real*8 z,r,phi !space coordinates

c-----locals
cBH090803      Some compilers don't like declaration of drhodz above and
cMH090803      also here.
cBH090803      real*8 theta_l,rho_l,drhodz,psi_xr,dro_dpsi,z_b,r_b,
      real*8 theta_l,rho_l,psi_xr,dro_dpsi,z_b,r_b,
     &rho_LCFS,rho_non_normalized,
     &x_l,y_l

      real*8 step,rho_l_p,rho_l_m,numer_drhodz,
     &z_b_p,r_b_p,z_b_m,r_b_m,
     &d_z_b_dz,d_r_b_dz,d_rho_LCFS_dz


c-----externals
      real*8 thetapol,psif,drhopsi

      if(iboundb.ge.1)then
c--------the point is outside LCFS 

        theta_l=thetapol(z,r) 
        call zr_psith(psilim,theta_l,z_b,r_b)

cSAP090328    
c        rho_l=dsqrt(((r-rma)**2+(z-zma)**2)/((r_b-rma)**2+(z_b-zma)**2))
c        drhodz=(z-zma)/(((r_b-rma)**2+(z_b-zma)**2)*rho_l)
        rho_LCFS=dsqrt((r_b-rma)**2+(z_b-zma)**2)
        rho_non_normalized=dsqrt((r-rma)**2+(z-zma)**2)
        rho_l=rho_non_normalized/rho_LCFS

c       drhodz=(z-zma)/(rho_non_normalized*rho_LCFS)
c       Terms with d_r_b/dz and d_r_b/dz should be added
c
c-------Calculate numerical derivatives d/dz from r_b, z_b
        x_l=r-rma
        y_l=z-zma
        step=1.d-4
        call theta_rz(x_l,y_l+step,theta_l)
        call zr_psith(psilim,theta_l,z_b_p,r_b_p)
        call theta_rz(x_l,y_l-step,theta_l)
        call zr_psith(psilim,theta_l,z_b_m,r_b_m)
        d_z_b_dz=(z_b_p-z_b_m)/(2.d0*step)     
        d_r_b_dz=(r_b_p-r_b_m)/(2.d0*step)
        d_rho_LCFS_dz=((r_b-rma)*d_r_b_dz+(z_b-zma)*d_z_b_dz)/rho_LCFS
ctest----------
c        write(*,*)'d_z_b_dz,d_r_b_dz,d_rho_LCFS_dz',
c     &             d_z_b_dz,d_r_b_dz,d_rho_LCFS_dz
c        drhodz=(z-zma)/(rho_non_normalized*rho_LCFS)
c        write(*,*)'in drhodr old drhodz',drhodz
c-end----- test  
        drhodz=(z-zma)/(rho_non_normalized*rho_LCFS)-
     &  (rho_l/rho_LCFS)*d_rho_LCFS_dz

c        write(*,*)'in drhodz new drhodz',drhodz

c-------test numerical derivative drhodz
c        step=1.d-4
c-------z+step
c        call theta_rz(x_l,y_l+step,theta_l)
c        call zr_psith(psilim,theta_l,z_b,r_b)
c        rho_LCFS=dsqrt((r_b-rma)**2+(z_b-zma)**2)        
c        rho_non_normalized=dsqrt((r-rma)**2+(z+step-zma)**2) 
c        rho_l_p=rho_non_normalized/rho_LCFS

c-------z-step 
c        call theta_rz(x_l,y_l-step,theta_l)
c        call zr_psith(psilim,theta_l,z_b,r_b)
c        rho_LCFS=dsqrt((r_b-rma)**2+(z_b-zma)**2)
c        rho_non_normalized=dsqrt((r-rma)**2+(z-step-zma)**2) 
c        rho_l_m=rho_non_normalized/rho_LCFS
        
c        numer_drhodz=(rho_l_p-rho_l_m)/(2.d0*step) !numerical derivativ drho_dz

c        write(*,*)'in drhodz drhodz,numer_drhodz',drhodz,numer_drhodz
      else
        psi_xr=psif(z,r)
        dro_dpsi=drhopsi(psi_xr)
        drhodz=dro_dpsi*dpdzd ! dpdzd was calucated in b ( in one.i)      
      endif
      return
      END

      real*8 FUNCTION dtempdr(z,r,phi,i)
c-------------------------------------------------------
c     derivative temperature by r (from spline) T-Kev
c-------------------------------------------------------
c     b(z,r,phi) should be run before this function
      implicit none
c     input
      double precision z,r,phi
      integer i
c     external rhof,dtemdrho,drgodr
      double precision rhof,dtemdrho,drhodr
c     local
      double precision rholoc

      rholoc=rhof(z,r,phi)
      dtempdr=dtemdrho(rholoc,i)*drhodr(z,r,phi)
      
      return
      end

      real*8 FUNCTION dtempdz(z,r,phi,i)
c-------------------------------------------------------
c     derivative temperature by z (from spline) T-Kev
c-------------------------------------------------------
c     b(z,r,phi) should be run before this function
      implicit none
c     input
      double precision z,r,phi
      integer i
c     external rhof,dtemdrho,drhodz
      double precision rhof,dtemdrho,drhodz
c     local
      double precision rholoc

      rholoc=rhof(z,r,phi)
      dtempdz=dtemdrho(rholoc,i)*drhodz(z,r,phi)
      
      return
      end


      real*8 FUNCTION dvflowdz(z,r,phi,i)
c-------------------------------------------------------
c     derivative vflow by z (from spline) cm/sec
c-------------------------------------------------------
c     b(z,r,phi) should be run before this function
      implicit none
c     input
      double precision z,r,phi
      integer i
c     external rhof,dvflwrho,drhodz
      double precision rhof,dvflwrho,drhodz
c     local
      double precision rholoc

      rholoc=rhof(z,r,phi)
      dvflowdz=dvflwrho(rholoc,i)*drhodz(z,r,phi)
      
      return
      end

      real*8 FUNCTION dvflowdr(z,r,phi,i)
c-------------------------------------------------------
c     derivative vflow by r (from spline) cm/sec
c-------------------------------------------------------
c     b(z,r,phi) should be run before this function
      implicit none
c     input
      double precision z,r,phi
      integer i
c     external rhof,dvflwrho,drhodr
      double precision rhof,dvflwrho,drhodr
c     local
      double precision rholoc

      rholoc=rhof(z,r,phi)
      dvflowdr=dvflwrho(rholoc,i)*drhodr(z,r,phi)
      
      return
      end

      real*8 FUNCTION dtpopdz(z,r,phi,i)
c-------------------------------------------------------
c     derivative tpop by z (from spline) 
c     i is the number of the plasma component
c-------------------------------------------------------
c     b(z,r,phi) should be run before this function
      implicit none
c     input
      double precision z,r,phi
      integer i
c     external rhof,dtpoprho,drhodz
      double precision rhof,dtpoprho,drhodz
c     local
      double precision rholoc

      rholoc=rhof(z,r,phi)
      dtpopdz=dtpoprho(rholoc,i)*drhodz(z,r,phi)
      
      return
      end

      real*8 FUNCTION dtpopdr(z,r,phi,i)
c-------------------------------------------------------
c     derivative tpop by r (from spline) 
c     i is the number of the plasma component
c-------------------------------------------------------
c     b(z,r,phi) should be run before this function
      implicit none
c     input
      double precision z,r,phi
      integer i
c     external rhof,dtpoprho,drhodr
      double precision rhof,dtpoprho,drhodr
c     local
      double precision rholoc

      rholoc=rhof(z,r,phi)
      dtpopdr=dtpoprho(rholoc,i)*drhodr(z,r,phi)
      
      return
      end


