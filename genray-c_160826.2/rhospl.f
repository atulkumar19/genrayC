c*************************rhospl**************************************
c  creation of arrays for  spline coefficients	         	      *
c  for rho_s:    trhos(npsi4),crhspsi(npsi4)			      *
c  for rho_bt:   trhobt(npsi4),crhbtpsi(npsi4)			      *
c  for rho_v:    trhov(npsi4),crhvpsi(npsi4)			      *
c  for rho^2_s:  trhos2(npsi4),crhsps2(npsi4)			      *
c  for rho^2_bt: trhobt2(npsi4),crhbtps2(npsi4)			      *
c  for rho^2_v:  trhov2(npsi4),crhvps2(npsi4)			      *
c  Sm050302
c  for rho_r:    trhor(npsi4),crhrpsi(npsi4)
c  for rho^2_r:  trhor2(npsi4),crhrps2(npsi4)
c  ---------------------------------------------		      *
c  for psi(rho): tpsi(npsi4),cpsirho(npsi4)			      *
c  ---------------------------------------------                      *
c  for z(psi,theta) : txz,tyz,czxy ,(czy-working array)
c  for r(psi,theta) : txr,tyr,crxy ,(czy-working array)
c  ---------------------------------------------                      *
c  rmax_psi(psi) : trmaxpsi(npsi4),crmaxpsi(npsi4)
c  rmin_psi(psi) : trminpsi(npsi4),crminpsi(npsi4)
c  bmax_psi(psi) : tbmaxpsi(npsi4),cbmaxpsi(npsi4)
c  bmin_psi(psi) : tbminpsi(npsi4),cbminpsi(npsi4)
c  arrays ar_min(npsi), ar_max(npsi),ab_max(npsi),ab_min(npsi)
c  will be calculated here inside subroutine rbmax
c  ---------------------------------------------                      *
c  this program uses the following subroutines:  		      *
c								      *
c  iac1r(rlimr,ip,ip4,zlimr,lx,mxa,fxa,mxb,fxb,trlimp,cxlimp)	      *
c  -calculates the spline coefficients for 1d function		      *
c
c								      *
c  IAC2R(rrr,nx,zzr,ny,peqdr,nfx,nfy,lx,ly,mxa,arfxa,mxb,arfxb,	      *
c  mya,arfya,myb,arfyb,tx,ty,cxy,ncx,ncy,nry,cy)		      *
c  - calculates the spline coefficients for 2d function		     
c
c new version of this subroutine (080815) uses spline subroutines
c coeff1 and terp1 from zcunix.f
c creation of arrys for spline coefficients
c for rho_s:  d2_rhos_psi(npsi)
c for rho_bt: d2_rhotf_psi(npsi)
c for rho_v:  d2_rhovl_psi(npsi)
c for rho_l:  d2_rhol_psi(npsi)
c for rho_r:  d2_rhor_psi(npsi)

c---------------------------------------------------------------------
 
c**********************************************************************
c  input data from common 'gr.i'        			      *
c  output data to common 'rho.i'        			      *
c**********************************************************************

      subroutine rhospl
c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none

      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
      include 'rho.i'
      include 'gr.i'

      double precision arfxa(nteta1),arfxb(nteta1),arfya(npsi),
     1 arfyb(npsi)
      real*8 area(npsi),torflux(npsi),volume(npsi)
c      real*8 arrho_s(npsi),arrho_tf(npsi),arrho_vl(npsi)
      real*8 arrho_psi(npsi) !for rho=sqrt(|(psi-psi_mag)/(psi_lim-psimag)|)
cSm010212
      real*8 arrho_psi2(npsi)!for rho=|(psi-psi_mag)/(psi_lim-psimag)|

      real*8 pollength(npsi)
cSAP081122
c      real*8 arrho_l(npsi)
         
      real*8 arho2_s(npsi),arho2_tf(npsi),arho2_vl(npsi),
     &arho2_r(npsi)                  !arrho_r()**2
cSm050304
c      double precision arrho_r(npsi),!for rho=(r_max(psi)-rmin(psi))/
c                                     !        (r_max(psi_lim)-rmin(psi_lim))
c     &arho2_r(npsi)                  !arrho_r()**2

      real*8
     &hteta,hpsi,sintd2,fxa,fxb,fold,psi,f,teta,btor,btorold,rho2ji,
     &rho2joi,p,pv,pb,pl,
     &rhopsi_,psi_l, dz_line, ravg

      integer
     &lx,ly,mxa,mxb,mya,myb,nfx,nfy,idx,nr4,j,i
c------------------------------------------------------
      external  ias1r,rmax_psi,rmin_psi,rhopsi,b
      real*8 ias1r,rmax_psi,rmin_psi,rhopsi,b
c----------------------------------------------------------
c     to use spline functions from zcunix.f:
c     coeff1 and terp1
c-----------------------------------------------------------
      integer i1p(2)
      real*8, dimension(1:3*npsi+1) :: work_l
   
      write(*,*)'rhospl: indexrho,npsi=',indexrho,npsi

      hteta=arteta(2)-arteta(1)
      hpsi=arpsi(2)-arpsi(1)
ctest
      do i=1,npsi
       write(*,*)'rhospl: i,arpsi(i)',i,arpsi(i)
      enddo
cendtest
       
      pi=4.d0*datan(1.d0)
      sintd2=dsin(0.5d0*hteta)
c--------------------------------------------------------------------
c     spline coefficients for r(psi,teta) and z(psi,teta) creation
c--------------------------------------------------------------------
      lx=1
      ly=1
      mxa=0
      mxb=0
c       mya=3 ! yup: presumably 3 means periodic conditions (?)
       mya=0
cc       myb=3 ! yup: presumably 3 means periodic conditions (?)
      myb=0
      nfx=npsi
      nfy=nteta1
c-----------------------------------------
c     ncx=npsi+4
c     ncy=nteta1+4
c     nzy=max0(npsi,nteta1)+4
c-----------------------------------------

      call iac2r(arpsi,npsi,arteta,nteta1,rpsi,nfx,nfy,lx,ly,mxa,arfxa,
     1 mxb,arfxb,mya,arfya,myb,arfyb,txr,tyr,crxy,npsi4,nteta14,nzy,czy)
     
      call iac2r(arpsi,npsi,arteta,nteta1,zpsi,nfx,nfy,lx,ly,mxa,arfxa,
     1 mxb,arfxb,mya,arfya,myb,arfyb,txz,tyz,czxy,npsi4,nteta14,nzy,czy)
c--------------------------------------------------------------------
c     spline coefficients for z(psi,teta) r(psi,teta) were created
c--------------------------------------------------------------------
CSm050302
c------------------------------------------------------------------
c  calculations of arrays
c  ar_min(npsi), ar_max(npsi), ab_max(npsi) ab_min(npei)
c------------------------------------------------------------------
c      write(*,*)'rhospl before call rbmax' 
      if(ixyz.eq.0) then ! r-phi coords
         call rbmax  
      else
         call rbmax_xyz
      endif
ctest
      do i=1,npsi
         write(*,*)'rhospl: i,ar_min,ar_max',i,ar_min(i),ar_max(i)
      enddo
cendtest
c------------------------------------------------------------------
c     spline coefficient calculation for rmax_psi(psi)
c------------------------------------------------------------------
      call iac1r(arpsi,npsi,npsi4,ar_max,lx,mxa,fxa,mxb,fxb,trmaxpsi
     1           ,crmaxpsi)
c------------------------------------------------------------------
c     spline coefficient calculation for rmin_psi(psi)
c------------------------------------------------------------------
      call iac1r(arpsi,npsi,npsi4,ar_min,lx,mxa,fxa,mxb,fxb,trminpsi
     1           ,crminpsi)

ctest
      write(*,*)'in rhospl.f'
c      write(*,*)'rmin_psi(psimag)',rmin_psi(psimag)
c      write(*,*)'rmax_psi(psimag)',rmax_psi(psimag)
c      write(*,*)'rma',rma
      write(*,*)'psilim,psimag',psilim,psimag
c      write(*,*)'arpsi',arpsi
c      stop 'rhospl.f'
cendtest
c--------------------------------------------------------------------
c     spline coefficients for r_max_psi(psi) rmin_psi) were created
c--------------------------------------------------------------------
c     spline coefficient calculation for bmax_psi(psi)
c------------------------------------------------------------------
      call iac1r(arpsi,npsi,npsi4,ab_max,lx,mxa,fxa,mxb,fxb,tbmaxpsi
     1           ,cbmaxpsi)

c     spline coefficient calculation for bmin_psi(psi)
c-----------------------------------------------------------------
      call iac1r(arpsi,npsi,npsi4,ab_min,lx,mxa,fxa,mxb,fxb,tbminpsi
     1           ,cbminpsi)
c--------------------------------------------------------------------
c     spline coefficients for bmax_psi(psi,) bmin_psi(psi) were created
c--------------------------------------------------------------------
c             nx=npsi
c	      ny=nteta1
c	      ncx=nx+4=npsi4
c	      ncy=ny+4=nteta14
c-------------------------------------------------------------------
       area(1)=0.d0
       torflux(1)=0.d0
       volume(1)=0.d0
       pollength(1)=0.d0 ! poloidal length of the magnetic surface
       arrho_r(1)=0.d0 !R_max(psi)-R(min(psi)/
                       ! (rmax_psi(arpsi(npsi))-rmin_psi(arpsi(npsi)))
                       
                       
      if(eqdsktype.eq.'mirror')then ! open field lines
        j=1 ! rho=0, i.e. the axis R=0
        pollength(j)=0.d0 ! initialization of the poloidal length
        psi=arpsi(j)
        arrho_r(j)=(rmax_psi(psi)-rmin_psi(psi))/
     &             (rmax_psi(arpsi(npsi))-rmin_psi(arpsi(npsi)))
        do i=1,nteta
            btor=0.d0 !f/rpsi(j,i) 
            btorold=0.d0 !fold/rpsi(j-1,i)
            ! For a given field line j: 
            ! consider a section (cut) of the line [i;i+1]
            dz_line= zpsi(j,i+1)-zpsi(j,i) ! Z-size of one section
            ravg= 0.5*(rpsi(j,i+1)+rpsi(j,i)) ! avr R at field line
            !for j=1 line, ravg is zero at any i-section.
            !Area to the left of the field-line section: dz_line*ravg
            p= dz_line*ravg*2 ![m^2] ! 2 here is for left|right symmetry
             !(line is traced over |)) region; area is over ((|)) region).
             !In other words, the field line was traced only over right 
             !region r>rmag, not the whole closed surface as in tokamaks. 
            !Volume of the plasma:
            pv= pi*dz_line*ravg*ravg ! INTEGRAL(pi*R^2*dZ) [m^3]
            area(j)=area(j)+p 
            volume(j)=volume(j)+pv
            torflux(j)=torflux(j)+p*0.5d0*(btor+btorold)
            pl=dsqrt((rpsi(j,i+1)-rpsi(j,i))**2+
     1               (zpsi(j,i+1)-zpsi(j,i))**2)
            pollength(j)=pollength(j)+pl
            write(*,'(a,4e12.4)')'rhospl: rpsi,zpsi, pollength,volume',
     +       rpsi(j,i),zpsi(j,i), pollength(j),volume(j) 
        enddo
        !write(*,*)'done with j=',j, '  psi=',psi
        !pause
        
        do j=2,npsi
        pollength(j)=0.d0 ! initialization of the poloidal length
        psi=arpsi(j)
        arrho_r(j)=(rmax_psi(psi)-rmin_psi(psi))/
     &             (rmax_psi(arpsi(npsi))-rmin_psi(arpsi(npsi)))
        do i=1,nteta
            btor=0.d0 !f/rpsi(j,i) 
            btorold=0.d0 !fold/rpsi(j-1,i)
            ! For a given field line j: 
            ! consider a section (cut) of the line [i;i+1]
            dz_line= zpsi(j,i+1)-zpsi(j,i) ! Z-size of one section
            ravg= 0.5*(rpsi(j,i+1)+rpsi(j,i)) ! avr R at field line
            !for j=1 line, ravg is zero at any i-section.
            !Area to the left of the field-line section: dz_line*ravg
            p= dz_line*ravg*2 ![m^2] ! 2 here is for left|right symmetry
             !(line is traced over |)) region; area is over ((|)) region).
             !In other words, the field line was traced only over right 
             !region r>rmag, not the whole closed surface as in tokamaks. 
            !Volume of the plasma:
            pv= pi*dz_line*ravg*ravg ! INTEGRAL(pi*R^2*dZ) [m^3]
            area(j)=area(j)+p 
            volume(j)=volume(j)+pv
            torflux(j)=torflux(j)+p*0.5d0*(btor+btorold)
            pl=dsqrt((rpsi(j,i+1)-rpsi(j,i))**2+
     1               (zpsi(j,i+1)-zpsi(j,i))**2)
            pollength(j)=pollength(j)+pl
            write(*,'(a,4e12.4)')'rhospl: rpsi,zpsi, pollength,volume',
     +       rpsi(j,i),zpsi(j,i), pollength(j),volume(j) 
        enddo
        !write(*,*)'done with j=',j, '  psi=',psi
        !pause
        enddo ! j=2:npsi
        write(*,*)'volume(npsi),area(npsi)= ',volume(npsi),area(npsi)
        write(*,*)'---------------------------------------------------'
      endif ! eqdsktype.eq.'mirror'
         
         
      if((eqdsktype.eq.'tokamak') .or. (eqdsktype.eq.'TAE'))then
         ! applicable to closed field-lines only
      idx=0
      fold=ias1r(txf,nx,nr4,cx,idx,psimag)
      do 180 j=2,npsi
         area(j)=area(j-1)
         torflux(j)=torflux(j-1)
         volume(j)=volume(j-1)
         pollength(j)=0.d0 ! initialization of the poloidal length
         psi=arpsi(j)
         arrho_r(j)=(rmax_psi(psi)-rmin_psi(psi))/
     &              (rmax_psi(arpsi(npsi))-rmin_psi(arpsi(npsi)))
c                write(*,*)'rhospl.f j,psi,arrho_r(j)',j,psi,arrho_r(j)
c                write(*,*)'rmax_psi(psi)',rmax_psi(psi)
         idx=0
         f=ias1r(txf,nx,nr4,cx,idx,psi)
         
         do 190 i=1,nteta
            teta=arteta(i)
c------------------------------------------------------------------
            btor=f/rpsi(j,i) 
            btorold=fold/rpsi(j-1,i)
c------------------------------------------------------------------
            rho2ji=(rpsi(j,i)-rma)**2+(zpsi(j,i)-zma)**2
            rho2joi=(rpsi(j-1,i)-rma)**2+(zpsi(j-1,i)-zma)**2
            p=(rho2ji-rho2joi)*sintd2
            area(j)=area(j)+p
            pv=0.5d0*(rpsi(j-1,i)+rpsi(j,i))*2.d0*pi
            pb=0.5d0*(btor+btorold)
            volume(j)=volume(j)+p*pv
            torflux(j)=torflux(j)+p*pb
            pl=dsqrt((rpsi(j,i+1)-rpsi(j,i))**2+
     1                 (zpsi(j,i+1)-zpsi(j,i))**2)
            pollength(j)=pollength(j)+pl
            write(*,'(a,4e12.4)')'rhospl: rpsi,zpsi, pollength,volume',
     +       rpsi(j,i),zpsi(j,i), pollength(j),volume(j) 
190     continue ! poloidal iteta
      write(*,*)'done with j=',j, '  psi=',psi
      fold=f
180   continue ! 2:npsi
      endif !((eqdsktype.eq.'tokamak') .or. (eqdsktype.eq.'TAE'))
      

             
c              write(*,*)'in rhospl.f'
c              write(*,*)'arrho_r',arrho_r
              
c------------------------------------------------------------------------
c             areatot in m**2
	      areatot=area(npsi)
              areatot=areatot*r0x**2
c------------------------------------------------------------------------
c             voltot in m**3
	      voltot=volume(npsi)
              voltot=voltot*r0x**3
              write(*,*)'rhospl voltot m**3',voltot
c------------------------------------------------------------------------
c             torflux in tesla*m**2
              torftot=torflux(npsi)
              torftot=torftot*r0x**2*b0
              write(*,*)'rhospl torftot tesla*m**2',torftot
c------------------------------------------------------------------------
c             poloidal length of the last magnetic surface in m
              totlength=pollength(npsi)
              totlength=totlength*r0x
              write(*,*)'rhospl totlength m',totlength
c------------------------------------------------------------------------
              do 200 j=1,npsi

                arrho_s(j)=dsqrt(area(j)/area(npsi))

                arrho_tf(j)=dsqrt(torflux(j)/torflux(npsi))
                arrho_vl(j)=dsqrt(volume(j)/volume(npsi))
                arrho_psi(j)=dsqrt((arpsi(j)-arpsi(1))/
     1                             (arpsi(npsi)-arpsi(1)))

                arrho_psi2(j)=(arpsi(j)-arpsi(1))/
     1                        (arpsi(npsi)-arpsi(1))

                arrho_l(j)=pollength(j)/pollength(npsi)

                arho2_s(j)=area(j)/area(npsi)
                arho2_tf(j)=torflux(j)/torflux(npsi)
                arho2_vl(j)=volume(j)/volume(npsi)
                arho2_r(j)=arrho_r(j)**2
200           continue
            
c-----------------------
c      ip4=npsi4
c----------------------
      lx=1
      mxa=0
      fxa=0.0d0
      mxb=0
      fxb=0.0d0
c------------------------------------------------------------------
c     spline coefficient calculation for rho(psi) using old spline
c------------------------------------------------------------------
      call iac1r(arpsi,npsi,npsi4,arrho_s,lx,mxa,fxa,mxb,fxb,trhos,
     1            crhspsi)
      call iac1r(arpsi,npsi,npsi4,arrho_tf,lx,mxa,fxa,mxb,fxb,trhobt,
     1            crhbtpsi)
      call iac1r(arpsi,npsi,npsi4,arrho_vl,lx,mxa,fxa,mxb,fxb,trhov,
     1            crhvpsi)
      call iac1r(arpsi,npsi,npsi4,arrho_l,lx,mxa,fxa,mxb,fxb,trholpol,
     1            crhlpsi)
      call iac1r(arpsi,npsi,npsi4,arrho_r,lx,mxa,fxa,mxb,fxb,trhor,
     1            crhrpsi)

c------------------------------------------------------------------
c     spline coefficient calculation for rho(psi) using spline
c     from zcunix.f
c------------------------------------------------------------------
      i1p(1)=4
      i1p(2)=4
      call coeff1(npsi,arpsi,arrho_s,d2_rhos_psi,i1p,1,work_l)
      call coeff1(npsi,arpsi,arrho_tf,d2_rhotf_psi,i1p,1,work_l)
      call coeff1(npsi,arpsi,arrho_vl,d2_rhovl_psi,i1p,1,work_l) 
      call coeff1(npsi,arpsi,arrho_l,d2_rhol_psi,i1p,1,work_l)
      call coeff1(npsi,arpsi,arrho_r,d2_rhor_psi,i1p,1,work_l)

c------------------------------------------------------------------
c     spline coefficient calculation for rho**2(psi)
c------------------------------------------------------------------
      call iac1r(arpsi,npsi,npsi4,arho2_s,lx,mxa,fxa,mxb,fxb,trhos2,
     1            crhsps2)
      call iac1r(arpsi,npsi,npsi4,arho2_tf,lx,mxa,fxa,mxb,fxb,trhobt2,
     1            crhbtps2)
      call iac1r(arpsi,npsi,npsi4,arho2_vl,lx,mxa,fxa,mxb,fxb,trhov2,
     1            crhvps2)
      call iac1r(arpsi,npsi,npsi4,arho2_r,lx,mxa,fxa,mxb,fxb,trhor2,
     1            crhrps2)
   
c------------------------------------------------------------------
c     spline coefficient calculation for rho_s(rho) and rho_v(rho)
cSm990109 
c     and rho_lpol(rho)  
cSm050304
c     and rho_r(rho)
c------------------------------------------------------------------
      if(indexrho.eq.1) then
      call iac1r(arrho_s,npsi,npsi4,arrho_vl,lx,mxa,fxa,mxb,fxb,trhov_rh
     1           ,crhvrho)
      call iac1r(arrho_s,npsi,npsi4,arrho_l,lx,mxa,fxa,mxb,fxb,trhol_rh
     1           ,crhlrho)
      call iac1r(arrho_s,npsi,npsi4,arrho_r,lx,mxa,fxa,mxb,fxb,trhor_rh
     1           ,crhrrho)
      endif

      if(indexrho.eq.2) then
         call iac1r(arrho_tf,npsi,npsi4,arrho_s,lx,mxa,fxa,mxb,fxb,
     1              trhos_rh,crhsrho)
         call iac1r(arrho_tf,npsi,npsi4,arrho_vl,lx,mxa,fxa,mxb,fxb,
     1              trhov_rh,crhvrho)
         call iac1r(arrho_tf,npsi,npsi4,arrho_l,lx,mxa,fxa,mxb,fxb,
     1              trhol_rh,crhlrho)
         call iac1r(arrho_tf,npsi,npsi4,arrho_r,lx,mxa,fxa,mxb,fxb,
     1              trhor_rh,crhrrho)
      endif

      if(indexrho.eq.3) then
         call iac1r(arrho_vl,npsi,npsi4,arrho_s,lx,mxa,fxa,mxb,fxb,
     1              trhos_rh,crhsrho)
         call iac1r(arrho_vl,npsi,npsi4,arrho_l,lx,mxa,fxa,mxb,fxb,
     1              trhol_rh,crhlrho)
         call iac1r(arrho_vl,npsi,npsi4,arrho_r,lx,mxa,fxa,mxb,fxb,
     1              trhor_rh,crhrrho)
      endif

      if(indexrho.eq.4) then
         call iac1r(arrho_psi,npsi,npsi4,arrho_s,lx,mxa,fxa,mxb,fxb,
     1              trhos_rh,crhsrho)
         call iac1r(arrho_psi,npsi,npsi4,arrho_vl,lx,mxa,fxa,mxb,fxb,
     1              trhov_rh,crhvrho)
         call iac1r(arrho_psi,npsi,npsi4,arrho_l,lx,mxa,fxa,mxb,fxb,
     1              trhol_rh,crhlrho)
         call iac1r(arrho_psi,npsi,npsi4,arrho_r,lx,mxa,fxa,mxb,fxb,
     1              trhor_rh,crhrrho)
      endif

      if(indexrho.eq.5) then
         call iac1r(arrho_psi2,npsi,npsi4,arrho_s,lx,mxa,fxa,mxb,fxb,
     1              trhos_rh,crhsrho)
         call iac1r(arrho_psi2,npsi,npsi4,arrho_vl,lx,mxa,fxa,mxb,fxb,
     1              trhov_rh,crhvrho)
         call iac1r(arrho_psi2,npsi,npsi4,arrho_l,lx,mxa,fxa,mxb,fxb,
     1              trhol_rh,crhlrho)
         call iac1r(arrho_psi2,npsi,npsi4,arrho_r,lx,mxa,fxa,mxb,fxb,
     1              trhor_rh,crhrrho)
      endif

      if(indexrho.eq.6) then
         call iac1r(arrho_r,npsi,npsi4,arrho_s,lx,mxa,fxa,mxb,fxb,
     1              trhos_rh,crhsrho)
         call iac1r(arrho_r,npsi,npsi4,arrho_vl,lx,mxa,fxa,mxb,fxb,
     1              trhov_rh,crhvrho)
         call iac1r(arrho_r,npsi,npsi4,arrho_l,lx,mxa,fxa,mxb,fxb,
     1              trhol_rh,crhlrho)
      endif
c------------------------------------------------------------------
c     spline coefficient calculation for psi_rho(rho)
c------------------------------------------------------------------
      if(indexrho.eq.1) then
         call iac1r(arrho_s,npsi,npsi4,arpsi,lx,mxa,fxa,mxb,fxb,
     1              tpsi,cpsirho)
      endif

      if(indexrho.eq.2) then
         call iac1r(arrho_tf,npsi,npsi4,arpsi,lx,mxa,fxa,mxb,fxb,
     1              tpsi,cpsirho)
      endif
      
      if(indexrho.eq.3) then
         call iac1r(arrho_vl,npsi,npsi4,arpsi,lx,mxa,fxa,mxb,fxb,
     1              tpsi,cpsirho)
      endif

      if(indexrho.eq.6) then
         call iac1r(arrho_r,npsi,npsi4,arpsi,lx,mxa,fxa,mxb,fxb,
     1              tpsi,cpsirho)
      endif
    
      return
      end

c====================================================================
c====================================================================

       real*8
     1 function rhopsi(psi)
c------------------------------------------------------------------
c     this function calculates normalized radial coordinate
c     from poloidal flux 'psi'=< psilim
c     if indexrho=1 rho=sqrt(area inside the flux surface)
c     if indexrho=2 rho=sqrt(toroidal flux)
c     if indexrho=3 rho=sqrt(volume inside the flux surface)
c     if indexrho=4 rho=sqrt((psi-psimag)/psilim-psimag))
c     if indexrho=5 rho=((psi-psimag)/psilim-psimag))
c     if indexrho=6 rho=(r_max(psi)-rmin(psi))/
c                       (r_max(psi_lim)-rmin(psi_lim))
c---------------------------------------------------------------------
c     input parameter: poloidal flux 'psi'
c---------------------------------------------------------------------
      implicit none

      include 'param.i'
      include 'gr.i'
      include 'one.i'
      include 'six.i'
      include 'rho.i'
      include 'three.i'

c-----input 
      real*8 psi !poloidal flux
c-----external
      real*8 ias1r
c-----locals
      integer idx
c-----for zcunix spline-------------------
      integer itabl(3)
      real*8  tabl(3)


      if(model_rho_dens.eq.0  .or.  model_rho_dens.eq.5) then
         ! proceed: only works for 0 and 5
      else
         stop 'rhopsi: Only valid for model_rho_dens=0 or 5'
      endif

        ! YuP[2016] Spline may give a value of psi 
        ! a little smaller than psimag (at r=0) => 
        ! may cause error in rhopsi
        if(psilim.gt.psimag)then
          if(psi.lt.psimag)then
            psi=psimag
          endif
        else ! psilim<psimag
          if(psi.gt.psimag)then
            psi=psimag
          endif
        endif
      
      if(model_b.ne.0) then
c         write(*,'(a,3e12.3)')
c     +    'rhopsi: psi,psimag,psilim=',psi,psimag,psilim
         ! YuP[01-2016] Sometimes, at r=0, spline gives value of
         ! psi not exactly equal to psimag.
         ! It can be somewhat smaller than psimag.
         ! Example: psimag=0, and psi=-0.2d-22.
         ! It results in negative value under dsqrt.
         ! Added abs()
         rhopsi=dsqrt(abs((psi-psimag)/(psilim-psimag))) ! can be >1.0
         return
      endif
      
      ! model_b=0:
      goto 20 !to rho calculations using spline from zcinix.f
      
      if (psi.lt.psimag) then
         rhopsi=0.d0
         goto 10
      endif

      if (psi.gt.psilim) then !outside the (hot) plasma: set rho>1 values
          if (indexrho.ne.5) then
             rhopsi=dsqrt(abs((psi-psimag)/(psilim-psimag)))
             !YuP: technically, the above is only valid for indexrho=4
          else ! indexrho=5 :
             rhopsi=abs((psi-psimag)/(psilim-psimag))
          endif
      else
	  idx=0
	  if(indexrho.eq.1) then
	       rhopsi=ias1r(trhos,npsi,npsi4,crhspsi,idx,psi)
	  endif
	  if(indexrho.eq.2) then
	       rhopsi=ias1r(trhobt,npsi,npsi4,crhbtpsi,idx,psi)
	  endif
	  if(indexrho.eq.3) then
	       rhopsi=ias1r(trhov,npsi,npsi4,crhvpsi,idx,psi)
	  endif
	  if(indexrho.eq.4) then
	       rhopsi=dsqrt(abs((psi-psimag)/(psilim-psimag)))
	  endif
	  if(indexrho.eq.5) then
	       rhopsi=abs((psi-psimag)/(psilim-psimag))
	  endif
	  if(indexrho.eq.6) then
	       rhopsi=ias1r(trhor,npsi,npsi4,crhrpsi,idx,psi)
	  endif
      endif

 20   continue
c---------------------------------------------------
c     rho calculation using spline from zcunix.f
c---------------------------------------------------
      itabl(1)=1
      itabl(2)=0
      itabl(3)=0

      if (psi.lt.psimag) then
         rhopsi=0.d0
         goto 10
      endif

      if (psi.gt.psilim) then
c         the point is outside the plasma
          if (indexrho.ne.5) then
             rhopsi=dsqrt(abs((psi-psimag)/(psilim-psimag)))
             !YuP: technically, the above is only valid for indexrho=4
          else ! indexrho=5 :
             rhopsi=abs((psi-psimag)/(psilim-psimag))
          endif
      else
	  idx=0
	  if(indexrho.eq.1) then
c	       rhopsi=ias1r(trhos,npsi,npsi4,crhspsi,idx,psi)
               call terp1(npsi,arpsi,arrho_s,d2_rhos_psi,
     &                    psi,1,tabl,itabl)
               rhopsi=tabl(1)
c               write(*,*)'rhospl.f rhopsi indexrho=1,rhopsi',rhopsi
	  endif
	  if(indexrho.eq.2) then
c	       rhopsi=ias1r(trhobt,npsi,npsi4,crhbtpsi,idx,psi)
               call terp1(npsi,arpsi,arrho_tf,d2_rhotf_psi,
     &                    psi,1,tabl,itabl)
               rhopsi=tabl(1)
c               write(*,*)'in rhospl psi,rhopsi', psi,rhopsi
	  endif
	  if(indexrho.eq.3) then
c	       rhopsi=ias1r(trhov,npsi,npsi4,crhvpsi,idx,psi)
               call terp1(npsi,arpsi,arrho_vl,d2_rhovl_psi,
     &                    psi,1,tabl,itabl)
               rhopsi=tabl(1)
	  endif
	  if(indexrho.eq.4) then
	       rhopsi=dsqrt(abs((psi-psimag)/(psilim-psimag)))
	  endif
	  if(indexrho.eq.5) then
	       rhopsi=abs((psi-psimag)/(psilim-psimag))
	  endif

	  if(indexrho.eq.6) then
c	       rhopsi=ias1r(trhor,npsi,npsi4,crhrpsi,idx,psi)
               call terp1(npsi,arpsi,arrho_r,d2_rhor_psi,
     &                    psi,1,tabl,itabl)
               rhopsi=tabl(1)
	  endif
       endif

 10    continue
     
       return
       end



       double precision
     1 function rhos(rhox)
c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none

      include 'param.i'
      include 'gr.i'
      include 'one.i'
      include 'six.i'
      include 'rho.i'
c-----input
      double precision rhox !small radius
c-----external
      double precision ias1r
c-----local
      integer idx

      if(model_b.ne.0) then
        stop 'rhos: Not setup for model_b>0'
      endif

c          write(*,*)'rhospl.f in rhos rhox=',rhox
	  idx=0
	  if(indexrho.eq.1) then
	       rhos=rhox
	  else
	       rhos=ias1r(trhos_rh,npsi,npsi4,crhsrho,idx,rhox)
	  endif
c          write(*,*)'rhos',rhos
       return
       end



	double precision
     1 function rhov(rhox)

c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none


      include 'param.i'
      include 'gr.i'
      include 'one.i'
      include 'six.i'
      include 'rho.i'

c-----input
      double precision rhox !small radius
c-----external
      double precision ias1r
c-----local
      integer idx

      idx=0
      if(model_b.ne.0) then
        stop 'rhov: Not setup for model_b>0'
      endif

      if(indexrho.eq.3) then
	 rhov=rhox
      else
	 rhov=ias1r(trhov_rh,npsi,npsi4,crhvrho,idx,rhox)
      endif

      return
      end



      double precision
     1 function rho_lrho(rhox)

c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none

      include 'param.i'
      include 'gr.i'
      include 'one.i'
      include 'six.i'
      include 'rho.i'
      
c-----input
      double precision rhox !small radius
c-----external
      double precision ias1r
c-----local
      integer idx

      if(model_b.ne.0) then
        stop 'rho_lrho: Not setup for model_b>0'
      endif

      idx=0	  
      rho_lrho=ias1r(trhol_rh,npsi,npsi4,crhlrho,idx,rhox)
	  
      return
      end

cSm050304
      double precision
     1 function rho_rrho(rhox)

c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none


      include 'param.i'
      include 'gr.i'
      include 'one.i'
      include 'six.i'
      include 'rho.i'
      
c-----input
      double precision rhox !small radius
c-----external
      double precision ias1r
c-----local
      integer idx

      if(model_b.ne.0) then
        stop 'rho_rrho: Not setup for model_b>0'
      endif

      idx=0

      if(indexrho.eq.6) then
	 rho_rrho=rhox
      else  
         rho_rrho=ias1r(trhor_rh,npsi,npsi4,crhrrho,idx,rhox)
      endif
	  
      return
      end


       double precision
     1 function rho_lpsi(psi)
c------------------------------------------------------------------
c     this function calculates normalized radial coordinate
c     rho_length_poloidal from poloidal flux 'psi'
c---------------------------------------------------------------------
c     input parameter: poloidal flux 'psi'
c---------------------------------------------------------------------
c       implicit integer (i-n), real*8 (a-h,o-z)
      implicit none

      include 'param.i'
      include 'gr.i'
      include 'one.i'
      include 'six.i'
      include 'rho.i'
      include 'three.i'

c-----input
      double precision psi !poloidal flux
c-----external
      double precision ias1r
c-----local
      integer idx

       stop 'rho_lpsi  Not setup for xyz?'
       
      if (psi.gt.psilim) then
c         the point is outside the plasma
cSm010212
          if (indexrho.ne.5) then
             rho_lpsi=dsqrt(abs((psi-psimag)/(psilim-psimag)))
             !YuP: technically, the above is only valid for indexrho=4   
          else
             rho_lpsi=abs((psi-psimag)/(psilim-psimag))
          endif

      else
	  idx=0
          rho_lpsi=ias1r(trholpol,npsi,npsi4,crhlpsi,idx,psi)	  
      endif

      return
      end




      double precision FUNCTION drhopsi(psi)

c------------------------------------------------------------------
c     this function calculates the derivatives from normalized radial
c     coordinate by the poloidal flux 'psi'
c     if indexrho=1 rho=sqrt(area inside the flux surface)
c     if indexrho=2 rho=sqrt(toroidal flux)
c     if indexrho=3 rho=sqrt(volume inside the flux surface)
c     if indexrho=4 rho=sqrt((psi-psimag)/psilim-psimag))
c     if indexrho=5 rho=((psi-psimag)/psilim-psimag))
c     if indexrho=6 rho=(r_max(psi)-rmin(psi))/
c                       (r_max(psi_lim)-rmin(psi_lim))
c---------------------------------------------------------------------
c     input parameter: poloidal flux 'psi'
c---------------------------------------------------------------------
c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none

      include 'param.i'
      include 'gr.i'
      include 'one.i'
      include 'six.i'
      include 'rho.i'
      include 'three.i'
      
c-----input
      double precision psi !poloidal flux
c-----external
      double precision ias1r
c-----local
      integer idx
c-----for zcunix spline-------------------
      integer itabl(3)
      real*8  tabl(3)

        ! YuP[2016] Spline may give a value of psi 
        ! a little smaller than psimag (at r=0) => 
        ! may cause error in rhopsi
        if(psilim.gt.psimag)then
          if(psi.lt.psimag)then
            psi=psimag
          endif
        else ! psilim<psimag
          if(psi.gt.psimag)then
            psi=psimag
          endif
        endif

      if(abs(psi-psimag).eq.0.d0)then
         drhopsi=0.d0 ! Technically, it can be INF at psi=psimag
         ! but other derivatives are 0 anyway.
         return
      endif
      
      if(model_b.ne.0) then
         drhopsi=0.5d0/dsqrt((psi-psimag)*(psilim-psimag))
         return
      endif

      goto 20 !to rho calculations using spline from zcinix.f

      if (psi.gt.psilim) then
c         the point is outside the plasma
          if (indexrho.ne.5) then
	     drhopsi=0.5d0/dsqrt((psi-psimag)*(psilim-psimag))
             !YuP: technically, the above is only valid for indexrho=4
          else
	     drhopsi=1.d0/(psilim-psimag)
          endif
      else
          idx=1
	  if(indexrho.eq.1) then
	       drhopsi=ias1r(trhos,npsi,npsi4,crhspsi,idx,psi)      
	  endif
	  if(indexrho.eq.2) then
	       drhopsi=ias1r(trhobt,npsi,npsi4,crhbtpsi,idx,psi)
          endif
	  if(indexrho.eq.3) then
	       drhopsi=ias1r(trhov,npsi,npsi4,crhvpsi,idx,psi)       
	  endif
	  if(indexrho.eq.4) then
	       drhopsi=0.5d0/dsqrt((psi-psimag)*(psilim-psimag))
	  endif
	  if(indexrho.eq.5) then
	       drhopsi=1.d0/(psilim-psimag)
	  endif
	  if(indexrho.eq.6) then
	       drhopsi=ias1r(trhor,npsi,npsi4,crhrpsi,idx,psi) 
          endif
      endif

 20   continue
c---------------------------------------------------
c    d_rho/d_psi calculation using spline from zcunix.f
c---------------------------------------------------
      itabl(1)=0
      itabl(2)=1
      itabl(3)=0

      if (psi.gt.psilim) then
c         the point is outside the plasma
          if (indexrho.ne.5) then
             drhopsi=0.5d0/dsqrt((psi-psimag)*(psilim-psimag))
             !YuP: technically, the above is only valid for indexrho=4
          else
             drhopsi=1.d0/(psilim-psimag)
          endif
      else
          idx=1
	  if(indexrho.eq.1) then
c	       drhopsi=ias1r(trhos,npsi,npsi4,crhspsi,idx,psi)
               call terp1(npsi,arpsi,arrho_s,d2_rhos_psi,
     &                    psi,1,tabl,itabl)
               drhopsi=tabl(2)
	  endif
	  if(indexrho.eq.2) then
c	       drhopsi=ias1r(trhobt,npsi,npsi4,crhbtpsi,idx,psi)
               call terp1(npsi,arpsi,arrho_tf,d2_rhotf_psi,
     &                    psi,1,tabl,itabl)
               drhopsi=tabl(2)
	  endif
	  if(indexrho.eq.3) then
c	       drhopsi=ias1r(trhov,npsi,npsi4,crhvpsi,idx,psi)
               call terp1(npsi,arpsi,arrho_vl,d2_rhovl_psi,
     &                    psi,1,tabl,itabl)
               drhopsi=tabl(2)
	  endif
	  if(indexrho.eq.4) then
	       drhopsi=0.5d0/dsqrt((psi-psimag)*(psilim-psimag))
	  endif
	  if(indexrho.eq.5) then
	       drhopsi=1.d0/(psilim-psimag)
	  endif

	  if(indexrho.eq.6) then
c	       drhopsi=ias1r(trhor,npsi,npsi4,crhrpsi,idx,psi)
               call terp1(npsi,arpsi,arrho_r,d2_rhor_psi,
     &                    psi,1,tabl,itabl)
               drhopsi=tabl(2)
	  endif
      endif
      

      return
      END





      double precision
     1 function drho2psi(psi)

c------------------------------------------------------------------
c     this function calculates the derivatives from normalized radial
c     coordinate rho**2 by the poloidal flux 'psi'
c     if indexrho=1 rho=sqrt(area inside the flux surface)
c     if indexrho=2 rho=sqrt(toroidal flux)
c     if indexrho=3 rho=sqrt(volume inside the flux surface)
c     if indexrho=4 rho=sqrt((psi-psimag)/(psilim-psimag))
c     if indexrho=5 rho=(psi-psimag)/(psilim-psimag)
c     if indexrho=6 rho=(r_max(psi)-rmin(psi))/
c                       (r_max(psi_lim)-rmin(psi_lim)
c---------------------------------------------------------------------
c     input parameter: poloidal flux 'psi'
c---------------------------------------------------------------------
c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none

      include 'param.i'
      include 'gr.i'
      include 'one.i'
      include 'six.i'
      include 'rho.i'
      include 'three.i'
      
c-----input
      double precision psi !poloidal flux
c-----external
      double precision ias1r
c-----local
      integer idx

        !stop 'drho2psi  Not setup for xyz?'
        
      idx=1
      
      if(model_b.ne.0) then
         drho2psi=1.d0/(psilim-psimag)
         return
      endif

	  if(indexrho.eq.1) then
	       drho2psi=ias1r(trhos2,npsi,npsi4,crhsps2,idx,psi)
	  endif
	  if(indexrho.eq.2) then
	       drho2psi=ias1r(trhobt2,npsi,npsi4,crhbtps2,idx,psi)
	  endif
	  if(indexrho.eq.3) then
	       drho2psi=ias1r(trhov2,npsi,npsi4,crhvps2,idx,psi)
	  endif
	  if(indexrho.eq.4) then
	       drho2psi=1.d0/(psilim-psimag)
	  endif
          if(indexrho.eq.5) then
	       drho2psi=2.d0*(psilim-psimag)/(psilim-psimag)**2
	  endif
cSm050304
          if(indexrho.eq.6) then
	       drho2psi=ias1r(trhor2,npsi,npsi4,crhrps2,idx,psi)
	  endif
      return
      end


      double precision
     1 function psi_rho(rhox)
c------------------------------------------------------------------
c     this function calculates poloidal flux 'psi'
c     on  normalized radial coordinate 'rho'
c------------------------------------------------------------------
c     if indexrho=1 rho=sqrt(area inside the flux surface)
c     if indexrho=2 rho=sqrt(toroidal flux)
c     if indexrho=3 rho=sqrt(volume inside the flux surface)
c     if indexrho=4 rho=sqrt((psi-psimag)/(psilim-psimag))
c     if indexrho=5 rho=(psi-psimag)/(psilim-psimag)
c     if indexrho=6 rho=(r_max(psi)-rmin(psi))/
c                       (r_max(psi_lim)-rmin(psi_lim)
c---------------------------------------------------------------------
c     input parameter: rho -radius
c---------------------------------------------------------------------
c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none

      include 'param.i'
      include 'gr.i'
      include 'one.i'
      include 'six.i'
      include 'rho.i'
      include 'three.i'
     
c-----input
      double precision rhox !small radius
c-----external
      double precision ias1r
c-----local
      integer idx
      
      if(model_b.ne.0) then
         psi_rho=psimag+(psilim-psimag)*rhox*rhox
         return
      endif

      ! For model_b=0 :
      if (rhox.gt.1.d0) then
c         the point is outside the plasma
          if (indexrho.ne.5) then
             psi_rho=psimag+(psilim-psimag)*rhox*rhox
             !YuP: technically, the above is only valid for indexrho=4
          else
             psi_rho=psimag+(psilim-psimag)*rhox
          endif
      else ! rhox<=1
          if(indexrho.eq.4) then  
             psi_rho=psimag+(psilim-psimag)*rhox*rhox
             return
          endif
          if(indexrho.eq.5) then
             psi_rho=psimag+(psilim-psimag)*rhox
             return
          endif
          idx=0
          psi_rho=ias1r(tpsi,npsi,npsi4,cpsirho,idx,rhox)
          return
      endif ! rhox
      
      return
      end




      double precision function qsafety_psi(psi)
c-----calculates safety factor q=d(toroidal flux)/d(poloidal flux)
c     tesla*m**2/[psi]

c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none

      include 'param.i'
      include 'one.i'
      include 'rho.i'
c-----input
      double precision  psi !poloidal flux
c-----externals
      double precision drho2psi
c     this function drho2psi(psi) calculates the derivatives from normalized radial
c     coordinate rho**2 by the poloidal flux 'psi'
c     if indexrho=1 rho=sqrt(area inside the flux surface)
c     if indexrho=2 rho=sqrt(toroidal flux)
c     if indexrho=3 rho=sqrt(volume inside the flux surface)
c     if indexrho=4 rho=sqrt((psi-psimag)/(psilim-psimag))
c     if indexrho=5 rho=(psi-psimag)/(psilim-psimag)
c     if indexrho=6 rho=(r_max(psi)-rmin(psi))/
c                       (r_max(psi_lim)-rmin(psi_lim)
c-----local
      integer indexrho_loc,indexrho_original
      double precision d_torflux_d_psi

      pi=4.d0*datan(1.d0)
      indexrho_original=indexrho

      indexrho=2   !toroidal flux
      d_torflux_d_psi=drho2psi(psi)*torftot !tesla*m**2/[psi?]
      indexrho=indexrho_original

      qsafety_psi=d_torflux_d_psi/(2.d0*pi)
      return
      end




      double precision function dvol_dpsi(psi)
c-----calculates the derivative d_volume/d_psi
c     m**3/[psi]

c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none
      include 'param.i'
      include 'one.i'
      include 'rho.i'
c-----input
      double precision  psi !poloidal flux
c-----externals
      double precision drho2psi
c     this function drho2psi(psi) calculates the derivatives from normalized radial
c     coordinate rho**2 by the poloidal flux 'psi'
c     if indexrho=1 rho=sqrt(area inside the flux surface)
c     if indexrho=2 rho=sqrt(toroidal flux)
c     if indexrho=3 rho=sqrt(volume inside the flux surface)
c     if indexrho=4 rho=sqrt((psi-psimag)/(psilim-psimag))
c     if indexrho=5 rho=(psi-psimag)/(psilim-psimag)
c     if indexrho=6 rho=(r_max(psi)-rmin(psi))/
c                       (r_max(psi_lim)-rmin(psi_lim)
c-----local
      integer indexrho_loc,indexrho_original

       !stop 'dvol_dpsi   Not setup for xyz?'
       
      indexrho_original=indexrho

      indexrho=3   !volume
      dvol_dpsi=drho2psi(psi)*voltot !m**3/[psi?]
      indexrho=indexrho_original

      return
      end


      double precision function b_average(psi_in)
c---------------------------------------------------------------------
c     calculates the afveraged magnetic field 
c     at the flux surface psi_in
c     It uses the function b th calcultes the mugnetic field.
c     It should be used after first call of subroutine rhospl
c     that calculates the cubic spline coefficients for small radius
c--------------------------------------------------------------------

c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none


      include 'param.i'
      include 'gr.i'
      include 'rho.i'
      include 'one.i'

      integer n_work
      parameter (n_work=3*npsi+1)
c-----input
      double precision  psi_in !poloidal flux
c-----external
       double precision b
c-----local
      double precision pollength(npsi),bmod_av(npsi),
     &work(3*npsi+1),d2bmod_psi(npsi),tabl(3),z,r,pl,psi_l

      integer i1p(2),itabl(3),i_first,i,j

      data i_first /0/
      save i_first,d2bmod_psi,bmod_av

c      write(*,*)'b_average i_first',i_first

      if (i_first.eq.0) then
c-------------------------------------------------------------------
c        calculate the array bmod_av(npsi) of the averaged magnetic field
c        along the magnetic surfaces
c        and calulate cubic spline coefficients
c--------------------------------------------------------------------
         pollength(1)=0.d0 ! poloidal length of the mugnetic surface
c         write(*,*)'b_average zpsi(1,1),rpsi(1,1)',zpsi(1,1),rpsi(1,1)
         z=zpsi(1,1)
         r=rpsi(1,1)
         bmod_av(1)=b(z,r,0.d0)
c         write(*,*)'bmod_av(1)',bmod_av(1)
         do 10 j=2,npsi             
            pollength(j)=0.d0 ! initialization of the poloidal length  
	    psi_l=arpsi(j)
            bmod_av(j)=0
c            write(*,*)'j,psi',j,psi
	    do 20 i=1,nteta
               z=0.5d0*(zpsi(j,i)+zpsi(j,i+1))
               r=0.5d0*(rpsi(j,i)+rpsi(j,i+1))
c               write(*,*)'i,r,z',i,r,z
               bmod=b(z,r,0.d0)
               pl=dsqrt((rpsi(j,i+1)-rpsi(j,i))**2+
     1                 (zpsi(j,i+1)-zpsi(j,i))**2)
               pollength(j)=pollength(j)+pl
               bmod_av(j)=bmod_av(j)+bmod*pl
c               write(*,*)'i,pl,bmod,bmod_av(j)',i,pl,bmod,bmod_av(j)
 20         continue
            bmod_av(j)=bmod_av(j)/pollength(j) !normalize at the poloidal lenght
c            write(*,*)'bmod_av(j),pollength(j)',bmod_av(j),pollength(j)
 10      continue

ctest
c        do j=1,npsi
c          write(*,*)'b_average j,arpsi(j),bmod_av(j)',
c     &    j,arpsi(j),bmod_av(j)
c          write(*,*)'bmin_psi(arpsi(j)),bmax_psi(arpsi(j))',
c     &    bmin_psi(arpsi(j)),bmax_psi(arpsi(j))
c        enddo !test
      

c----------------------------------------------------------
c        calculate cubic spline coefficients for the function
c        that gives the averaged magnetic field:  b_av_psi(psi)              
c---------------------------------------------------------       
         i1p(1)=4
         i1p(2)=4 
         call coeff1(npsi,arpsi,bmod_av,d2bmod_psi,i1p,1,work)
c-----------------------------------------------------------
        i_first=1
      endif !spline coefficients are calculated 

c-----------------------------------------------------------
c     calculation b_averaged 
c------------------------------------------------------------
      itabl(1)=1
      itabl(2)=0
      itabl(3)=0
      call terp1(npsi,arpsi,bmod_av,d2bmod_psi,psi_in,1,tabl,itabl)
      b_average=tabl(1)

      return
      end




      subroutine average_variables(psi_in,b_av,bs_av,r_av,dr_av,drs_av)
c---------------------------------------------------------------------
c     calculates the flux surface averaged variables
c     <B>, <B**2>, <r>, <1/r>, <1/r**2>
c     at the flux surface psi
c     <A>(psi_in) =Int{(A/|B_pol|)dl_pol}
c     
c     It uses the function b to calculates the magnetic field.
c     It should be used after first call of subroutine rhospl
c     that calculates the cubic spline coefficients for small radius
c--------------------------------------------------------------------

c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none
      include 'param.i'
      include 'gr.i'
      include 'rho.i'
      include 'one.i'
      include 'three.i'

      integer n_work
      parameter (n_work=3*npsi+1)

c-----input
      real*8 psi_in !poloidal flux

c-----output
      real*8
     &b_av,  !<B>
     &bs_av, !<B**2>
     &r_av,  !<R>
     &dr_av, !<1/R>
     &drs_av !<1/r**2>

c-----externals
      real*8 b
c-----local
      real*8 
     &pollength(npsi),
     &work(3*npsi+1),tabl(3),
     &b_av_ar(npsi),  d2_b_psi(npsi),
     &bs_av_ar(npsi), d2_bs_psi(npsi),
     &r_av_ar(npsi),  d2_r_psi(npsi),
     &dr_av_ar(npsi), d2_dr_psi(npsi),
     &drs_av_ar(npsi),d2_drs_psi(npsi),
     &dbpol_av_ar(npsi),
     &z,r,pl,b_pol,
     &psi_l

        
      integer i1p(2),itabl(3)
      integer i_first,i,j

      data i_first /0/
      save i_first,
     &b_av_ar,d2_b_psi,
     &bs_av_ar,d2_bs_psi,
     &r_av_ar,d2_r_psi,
     &dr_av_ar,d2_dr_psi,
     &drs_av_ar,d2_drs_psi

c      write(*,*)'average_variables i_first',i_first

      if (i_first.eq.0) then
c-------------------------------------------------------------------
c        calculate the array b_av(npsi) of the averaged magnetic field
c        along the magnetic surfaces
c        and calulate cubic spline coefficients
c--------------------------------------------------------------------
         pollength(1)=0.d0 ! poloidal length of the mugnetic surface
c         write(*,*)'b_average zpsi(1,1),rpsi(1,1)',zpsi(1,1),rpsi(1,1)
         z=zpsi(1,1)
         r=rpsi(1,1)
c         write(*,*)'rma,zma',rma,zma
c         write(*,*)'rpsi(1,1),zpsi(1,1)',rpsi(1,1),zpsi(1,1)
         b_av_ar(1)=b(z,r,0.d0)
         bs_av_ar(1)=b_av_ar(1)**2
         r_av_ar(1)= r
         dr_av_ar(1)= 1.d0/r
         drs_av_ar(1)= 1.d0/r**2
c         write(*,*)'b_av(1)',b_av(1)
         do 10 j=2,npsi             
            pollength(j)=0.d0 ! initialization of the poloidal length  
	    psi_l=arpsi(j)
            b_av_ar(j)=0.d0
            dbpol_av_ar(j)=0.d0
            bs_av_ar(j)=0.d0
            r_av_ar(j)=0.d0
            dr_av_ar(j)=0.d0
            drs_av_ar(j)=0.d0
            

c            write(*,*)'j,psi',j,psi
	    do 20 i=1,nteta
               z=0.5d0*(zpsi(j,i)+zpsi(j,i+1))
               r=0.5d0*(rpsi(j,i)+rpsi(j,i+1))
c               write(*,*)'i,r,z',i,r,z
               bmod=b(z,r,0.d0)
               b_pol=dsqrt(bz**2+br**2)
               pl=dsqrt((rpsi(j,i+1)-rpsi(j,i))**2+
     1                 (zpsi(j,i+1)-zpsi(j,i))**2)
               pollength(j)=pollength(j)+pl

               dbpol_av_ar(j)= dbpol_av_ar(j)+pl/b_pol
               b_av_ar(j)=b_av_ar(j)+bmod*pl/b_pol
               bs_av_ar(j)=bs_av_ar(j)+bmod**2*pl/b_pol
               r_av_ar(j)= r_av_ar(j)+r*pl/b_pol
               dr_av_ar(j)=dr_av_ar(j)+pl/(r*b_pol)
               drs_av_ar(j)=drs_av_ar(j)+pl/(r**2*b_pol)

c               write(*,*)'i,pl,bmod,b_av(j)',i,pl,bmod,b_av(j)
 20         continue
            b_av_ar(j)=b_av_ar(j)/dbpol_av_ar(j) !normalize
            bs_av_ar(j)=bs_av_ar(j)/dbpol_av_ar(j)
            r_av_ar(j)=r_av_ar(j)/dbpol_av_ar(j)
            dr_av_ar(j)=dr_av_ar(j)/dbpol_av_ar(j)
            drs_av_ar(j)=drs_av_ar(j)/dbpol_av_ar(j)

c           write(*,*)'j,b_av_ar(j),bs_av_ar(j)',j,b_av_ar(j),bs_av_ar(j) 
c           write(*,*)'r_av_ar(j),dr_av_ar(j),drs_av_ar(j)',
c     &                r_av_ar(j),dr_av_ar(j),drs_av_ar(j)
 10      continue

ctest
c        do j=1,npsi
c          write(*,*)'b_average j,arpsi(j),b_av(j)',
c     &    j,arpsi(j),b_av(j)
c          write(*,*)'bmin_psi(arpsi(j)),bmax_psi(arpsi(j))',
c     &    bmin_psi(arpsi(j)),bmax_psi(arpsi(j))
c        enddo !test
      

c----------------------------------------------------------
c        calculate cubic spline coefficients for the function
c        that gives the averaged magnetic field:  b_av_psi(psi)              
c---------------------------------------------------------       
         i1p(1)=4
         i1p(2)=4 
         call coeff1(npsi,arpsi, b_av_ar,  d2_b_psi,  i1p,1,work)
         call coeff1(npsi,arpsi, bs_av_ar, d2_bs_psi, i1p,1,work)
         call coeff1(npsi,arpsi, r_av_ar,  d2_r_psi,  i1p,1,work)
         call coeff1(npsi,arpsi, dr_av_ar, d2_dr_psi, i1p,1,work)
         call coeff1(npsi,arpsi, drs_av_ar,d2_drs_psi,i1p,1,work)
c-----------------------------------------------------------
        i_first=1

ctest_spline
c      itabl(1)=1
c      itabl(2)=0
c      itabl(3)=0
c      do j=1,npsi
c         psi_l=arpsi(j)
c         call terp1(npsi,arpsi,b_av_ar,d2_b_psi,psi_l,1,tabl,itabl)
c         b_av=tabl(1)
c         write(*,*)'j,psi_l,b_av_ar(j),b_av',j,psi_l,b_av_ar(j),b_av
c      enddo 
c      psi_l=-0.3988948224555455d0
c      call terp1(npsi,arpsi,b_av_ar,d2_b_psi,psi_l,1,tabl,itabl)
c      b_av=tabl(1)
c      write(*,*)'itabl',itabl
c      write(*,*)'rhospl.f test psi_l,b_av',psi_l,b_av
c      write(*,*)'arpsi',arpsi
c      write(*,*)'b_av_ar',b_av_ar
c      write(*,*)'d2_b_psi',d2_b_psi
c      write(*,*)'tabl',tabl
cendtest 

      endif !spline coefficients are calculated 

c-----------------------------------------------------------
c     calculation b_averaged 
c------------------------------------------------------------
      itabl(1)=1
      itabl(2)=0
      itabl(3)=0
c      write(*,*)'subroutine average_variables npsi,psi_in',npsi,psi_in
c      write(*,*)'before terp1_test'
c      call terp1_test(npsi,arpsi,b_av_ar,d2_b_psi,psi_in,1,tabl,itabl)
c      b_av=tabl(1)
c      write(*,*)'npsi',npsi
c      write(*,*)'itabl',itabl
c      write(*,*)'arpsi',arpsi
c      write(*,*)'b_av_ar',b_av_ar
c      write(*,*)'d2_b_psi',d2_b_psi
c      write(*,*)'tabl',tabl

c      call terp1_Sm(npsi,arpsi,b_av_ar,d2_b_psi,psi_in,1,tabl,itabl)
      call terp1(npsi,arpsi,b_av_ar,d2_b_psi,psi_in,1,tabl,itabl)
      b_av=tabl(1)
      call terp1(npsi,arpsi,bs_av_ar,d2_bs_psi,psi_in,1,tabl,itabl)
      bs_av=tabl(1) 
      call terp1(npsi,arpsi,r_av_ar,d2_r_psi,psi_in,1,tabl,itabl)
      r_av=tabl(1)
      call terp1(npsi,arpsi,dr_av_ar,d2_dr_psi,psi_in,1,tabl,itabl)
      dr_av=tabl(1)
      call terp1(npsi,arpsi,drs_av_ar,d2_drs_psi,psi_in,1,tabl,itabl)
      drs_av=tabl(1)
c      write(*,*)'b_av,bs_av,r_av,dr_av,drs_av',
c     &b_av,bs_av,r_av,dr_av,drs_av
      return
      end




      subroutine theta_rz(r,z,theta)
c-------------------------------------------------
c     calculate cilindrical angle theta (radian) at given
c     input coordinate r,z
c     0 =< theta < 2*pi
c--------------------------------------------------
      implicit none
c-----input 
      real*8 r,z
c-----output
      real*8 theta
c-----local 
      real*8 pi

      pi=4.d0*datan(1.d0)

      if (z.ge.0.d0) then
         theta=dacos(r/dsqrt(r**2+z**2)) ! [0; pi]
      else
         theta=2*pi-dacos(r/dsqrt(r**2+z**2))
      endif

      return
      end      


       
