
c======================================================================
c======================================================================
c Yuri Petrov   CompX  May-2011       All subroutines and functions
c in this file are used with ixyz=1 option (cartesian coords Genray-c)
c======================================================================
c======================================================================

      subroutine unit_vectors_perp_b(x,y,z, 
     +                     uvtx, uvty, uvtz,  uvnx, uvny, uvnz) !->out
c===> Yu.P. 2011
c--------------------------------------------------------------------
c     1. Calculate the components of unit vector  which is
c     perpendicular to b-field and TANGENTIAL to psi-flux-surface  
c--------------------------------------------------------------------
c     (uvtx) =                                   
c     (uvty) =     [grad(psi) x b] / |grad(psi)||b| = 
c     (uvtz) =                                   
c
c               = (dpsi/dy bz - dpsi/dz by) / |grad(psi)||b| 
c               = (dpsi/dz bx - dpsi/dx bz) / |grad(psi)||b| 
c               = (dpsi/dx by - dpsi/dy bx) / |grad(psi)||b| 
c--------------------------------------------------------------------
c     2. Calculate the components of unit vector  which is
c     perpendicular to b-field and NORMAL to psi-flux-surface  
c--------------------------------------------------------------------
c     (uvnx) =                       = (uvtz*by - uvty*bz) / |b|         
c     (uvny) =     [b x uvt] / |b|   = (uvtx*bz - uvtz*bx) / |b|
c     (uvnz) =                       = (uvty*bx - uvtx*by) / |b|               
c
c      Note: It can be shown that 
c            [b x uvt] / |b| = grad(psi)/|grad(psi)|
c            (based on (b . grad(psi)) =0)
c--------------------------------------------------------------------
      implicit none  
      include 'param.i'
      include 'one.i'
c-----input
      real*8 x,y,z
c-----output
      real*8 uvtx, uvty, uvtz ! tangential to psi
      real*8 uvnx, uvny, uvnz ! normal to psi
c-----externals
      real*8 bxyz
c-----locals
      real*8 grpsi, o_bgrpsi
      real*8 uvbx,uvby,uvbz      
      !---------------------
      bmod= bxyz(x,y,z) !-> components of b and grad(psi) -> /one.i/
      uvbx= bx*o_bmod ! unit vector component bx/b
      uvby= by*o_bmod ! unit vector component by/b
      uvbz= bz*o_bmod ! unit vector component bz/b
      grpsi= dsqrt(dpdxd**2 + dpdyd**2 + dpdzd**2) ! |grad(psi)|
      if (grpsi .ne. 0.d0) then
        o_bgrpsi= 1.d0/(grpsi*bmod) !   1 / |grad(psi)||b|
        ! Tangential to flux surface :
        uvtx= (dpdyd*bz - dpdzd*by)*o_bgrpsi
        uvty= (dpdzd*bx - dpdxd*bz)*o_bgrpsi
        uvtz= (dpdxd*by - dpdyd*bx)*o_bgrpsi
        ! Normal to flux surface :
        uvnx= (uvtz*by - uvty*bz)*o_bmod ! equivalent to dpdxd/grpsi
        uvny= (uvtx*bz - uvtz*bx)*o_bmod ! equivalent to dpdyd/grpsi
        uvnz= (uvty*bx - uvtx*by)*o_bmod ! equivalent to dpdzd/grpsi
c        write(*,'(a,3f10.3)') '   x,    y,    z=', x, y, z
c        write(*,'(a,3f10.3)') 'dpdxd,dpdyd,dpdzd',
c     +                         dpdxd/grpsi, dpdyd/grpsi, dpdzd/grpsi
c        write(*,'(a,3f10.3)') 'uvbx, uvby, uvbz=',uvbx, uvby, uvbz
c        write(*,'(a,3f10.3)') 'uvtx, uvty, uvtz=',uvtx, uvty, uvtz
c        write(*,'(a,3f10.3)') 'uvnx, uvny, uvnz=',uvnx, uvny, uvnz      
      else
        write(*,*) 'unit_vectors_perp_b:  x,y,z=',x,y,z
        stop 'unit_vectors_perp_b: grpsi==|grad(psi)|=0.0'  
        ! needs work for this case?
      endif
      return
      end


c======================================================================
c======================================================================

      subroutine nparper_to_nxyz(cnpar, cnper_tang, cnper_norm,
     +                           uvtx, uvty, uvtz,  uvnx, uvny, uvnz,
     +                           cnx,cny,cnz) !->out
c===> Yu.P. 2011
c--------------------------------------------------------------------
c     Calculate components of refractive index (Nx,Ny,Nz) 
c     in cartesian coordinates. 
c     Known: Npar, Nperp_tangential_to_psi, Nperp_normal_to_psi,
c     unit vector (uvtx,uvty,uvtz) tangential to psi flux surface,
c     unit vector (uvnx,uvny,uvnz) normal to psi flux surface.
c     (Nx,Ny,Nz) are found by solving a linear system
c     A*(Nx,Ny,Nz)' = (Npar, Nperp_tang, Nperp_norm)'
c     where matrix A is 
c     (bx/b   by/b   bz/b)
c     (uvtx   uvty   uvtz)
c     (uvnx   uvny   uvnz)
c--------------------------------------------------------------------
      implicit none  
      include 'param.i'
      include 'one.i' ! Contains bx,by,bz,bmod,o_bmod
c-----input
      real*8 cnpar, cnper_tang, cnper_norm,
     +       uvtx, uvty, uvtz,  uvnx, uvny, uvnz
c-----output
      real*8 cnx,cny,cnz ! (Nx,Ny,Nz) 
c-----locals
      real*8 uvbx,uvby,uvbz, deta, detx, dety, detz
      !---------------------
      uvbx= bx*o_bmod ! unit vector component bx/b
      uvby= by*o_bmod ! unit vector component by/b
      uvbz= bz*o_bmod ! unit vector component bz/b
      deta= uvbx*(uvty*uvnz-uvtz*uvny)
     +     +uvby*(uvtz*uvnx-uvtx*uvnz)
     +     +uvbz*(uvtx*uvny-uvty*uvnx)  ! determinant of A
      if(deta .ne. 0.d0) then
        ! Use Cramer's rule to find (Nx,Ny,Nz)
        detx= cnpar*(uvty*uvnz - uvtz*uvny)
     +       +uvby*(uvtz*cnper_norm - cnper_tang*uvnz)
     +       +uvbz*(cnper_tang*uvny - uvty*cnper_norm)  
     
        dety= uvbx*(cnper_tang*uvnz - uvtz*cnper_norm)
     +       +cnpar*(uvtz*uvnx - uvtx*uvnz)
     +       +uvbz*(uvtx*cnper_norm - cnper_tang*uvnx)  

        detz= uvbx*(uvty*cnper_norm - cnper_tang*uvny)
     +       +uvby*(cnper_tang*uvnx - uvtx*cnper_norm)
     +       +cnpar*(uvtx*uvny - uvty*uvnx)
     
        cnx= detx/deta
        cny= dety/deta
        cnz= detz/deta
      else ! deta=0.
        write(*,*) 'bx,by,bz=',bx,by,bz
        write(*,*) 'cnpar, cnper_tang, cnper_norm=',
     +              cnpar, cnper_tang, cnper_norm
        stop 'nparper_to_nxyz:  detA==0.0' ! needs work for this case?
      endif
      ! test
      if( abs( cnpar**2+cnper_tang**2+cnper_norm**2 -
     +   (cnx**2+cny**2+cnz**2) ) .gt. 1.d-5) then
        write(*,'(a,4f10.3)') 'cnpar, cnper_tang, cnper_norm, cn2=',
     +              cnpar, cnper_tang, cnper_norm,
     +              cnpar**2+cnper_tang**2+cnper_norm**2
        write(*,'(a,4f10.3)') 'cnx, cny, cnz, cn2=',
     +              cnx, cny, cnz, cnx**2+cny**2+cnz**2
        stop 'problem in nparper_to_nxyz'
      endif
      ! end test
      return
      end

c======================================================================
c======================================================================
      subroutine zrphi_to_xyz(r,phi, x,y)
c-----transform the cylindrical coords to cartesian coords.
      implicit integer (i-n), real*8 (a-h,o-z)
      !--------------
      x= r*dcos(phi)
      y= r*dsin(phi)
      ! z=z, so no work 
      return
      end

c======================================================================
c======================================================================
      subroutine xyz_to_zrphi(x,y, r,phi)
c-----transform the cartesian coords to cylindrical coords.
      implicit integer (i-n), real*8 (a-h,o-z)
c-----local
      data twopi/6.283185307179586476925287d0/
      !--------------
      phi= atan2(y,x)  ! Returns value in range [-pi,pi]
c      if(x.lt.0.d0) phi=phi+twopi ! [0,2pi]
      r= dsqrt(x**2+y**2)
      ! z=z, so no work 
      return
      end

c======================================================================
c======================================================================
      subroutine nzrphi_to_nxyz(r,phi,cnr,cm, cn_x,cn_y)
c-----transform the cylindrical coords to cartesian coords.
      implicit integer (i-n), real*8 (a-h,o-z)
      !--------------
      cnphi=cm/r
      sin_phi=dsin(phi)
      cos_phi=dcos(phi)
      cn_x= cnr*cos_phi - cnphi*sin_phi ! YuP added for now
      cn_y= cnr*sin_phi + cnphi*cos_phi ! YuP added for now
      return
      end

c======================================================================
c======================================================================
      subroutine nxyz_to_nzrphi(cn_x,cn_y,r,phi, cnr,cm)
c-----transform the cartesian coords to cylindrical coords.
      implicit integer (i-n), real*8 (a-h,o-z)
c-----local
      real*8 sin_phi,cos_phi
      !--------------
      sin_phi=dsin(phi)
      cos_phi=dcos(phi)
      cnr=   cn_x*cos_phi + cn_y*sin_phi
      cnphi=-cn_x*sin_phi + cn_y*cos_phi
      cm= r*cnphi
      return
      end
      
c======================================================================
c======================================================================      
      
c        ************************ bxyz ************************     
c        this function calculates the components  of  the   
c        magnetic field bx,by,bz, 
c        the absolute value of the total magnetic field bmod
c        and the derivatives of bx,by,bz,bmod by x,y,z. 
c	 						      
c        It controls if the ray point is inside the plasma  
c        if the ray point is outside the plasma then the    
c        calculations will be finished 		      
c        ******************************************************
c
c------------------------------------------------------------------
c								   
c        input parameters					   
c      								   
c      x,y,z - cartesian coords. where the components of 	   
c                 magnetic field are calculated.				   
c------------------------------------------------------------------
      double precision
     1function bxyz(x,y,z)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'fourb.i' ! contains req,zeq,xeq,yeq grids
      include 'five.i' ! contains rmax
      include 'six.i'
      real*8
     1 ias1r,ias2r_Sm
      real*8
     & thetapol, l_zr, l_z_b_r_b
c-----local
      real*8 spol, sin_phi,cos_phi

      data ifirstc/1/
      save ifirstc

c----------------------------------------------------------------
      r2= x*x+y*y
      r4= r2*r2
      r=  dsqrt(r2)
      dreq=req(2)-req(1) ! to use as a criteria of being close to r=0
                          ! if(r.lt.dreq) --> avoid 1/r division
c----------------------------------------------------------------
      ! Define spol == the sign of Bpoloidal
      spol=-1.d0   ! if(toteqd > 0.d0, along phi-direction) 
      if(toteqd.lt.0.d0) spol=1.d0 ! toteqd is plasma current from eqdsk
      spol=-spol  ! reversed again?
      if(toteqd.eq.0.d0) then ! This is the case for model_b>0
           if (ifirstc.eq.1) then
              write(*,*)'current is not given (toteqd=0)' 
              ifirstc=0
           endif
c          The total current is not given. The sign of the B_pol will be
c          determined from psi, using the original poloidal flux from eqdsk.
c          In this case the poloidal flux has max on the magnetic axis
c          (and dpsimax=-1), or it has min (and dpsimax=+1).
           spol=-dpsimax  ! YuP error? Changed to -dpsimax 05-05-2011
      endif
c----------------------------------------------------------------o    
      epsbnd=2.d-3 ![m]
      iboundb=-1
c--------------------------------------------------------------
      goto 10 ! YuP Added: skip boundary control and def. of rho
c--------------------------------------------------------------
c     control that the point is inside the plasma
!      r_l=r-rma
!      z_l=z-zma
!      if ((r_l**2+z_l**2).lt.1.d-10) then ! close to mag.axis
!         goto 10
!      endif
!
!c-----calculate poloidal angle  thetapol
!      call theta_rz(r_l,z_l,theta_l) !0 =< theta_l < 2*pi
!
!c-----calculate coordinates x_b,z_b at the limiter surface (or LCFS)
!c      psi=psilim  and the poloidal angle thetapol
!      call zr_psith(psilim,theta_l,z_b,r_b)
!      
!      l_zr=dsqrt(r_l**2+z_l**2)
!      l_z_b_r_b=dsqrt((r_b-rma)**2+(z_b-zma)**2)
!
!      if (l_zr.ge.l_z_b_r_b) then
!         rho = l_zr / l_z_b_r_b         
!         iboundb=3
!      endif
!      goto 10
!c----------------------------------------------------
!      if ((r.lt.rmin).or.(r.ge.rmax-epsbnd)) then
!         iboundb=2
!         theta_l=thetapol(z,r) 
!         call zr_psith(psilim,theta_l,z_b,r_b)
!         rho=dsqrt(((r-rma)**2+(z-zma)**2)/((r_b-rma)**2+(z_b-zma)**2))
!         goto 10
!      end if
!c---------------- idx derivativs order 0.ge.idx.le.3---------------
!      rrr=r
!      idx=0
!      ipx=ip
!      ipx4=ip+4
!      zzrp=ias1r(trlimp,ipx,ipx4,cxlimp,idx,rrr)
!      zp=zzrp
!      ipx=im
!      ipx4=im+4
!      zzrm=ias1r(trlimm,ipx,ipx4,cxlimm,idx,rrr)
!      zm=zzrm
!      if ((z.ge.zp-epsbnd).or.(z.le.zm+epsbnd)) then 
!         iboundb=1
!         theta_l=thetapol(z,r) 
!         call zr_psith(psilim,theta_l,z_b,r_b)
!         rho=dsqrt(((r-rma)**2+(z-zma)**2)/((r_b-rma)**2+(z_b-zma)**2))
!         goto 10
!      end if
c---------------------------------------------------------------------
c end of the control that the point is inside the plasma
 10   continue
c--------------------------------------------------------------
c
c nr4a , nz4a and nrya  are given in common five as parameters
c nreqda,nzeqda,nlimit are given in common four as parameters
c
      nr4=nx+4
      nz4=ny+4
      ncx=nr4
      ncy=nz4
      res=   ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,0,0,r,z,nr4a)
      dpdrd=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,1,0,r,z,nr4a) ! dpsi/dr
      dpdzd=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,0,1,r,z,nr4a) ! dpsi/dz
      psid=res
      dpdph= 0.d0 ! dpsi/dphi ! Zero for now, generalize later
      dpdzrd=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,1,1,r,z,nr4a) ! d2psi/dzdr
      dpdrrd=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,2,0,r,z,nr4a) ! d2psi/drdr
      dpdzzd=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,0,2,r,z,nr4a) ! d2psi/dzdz
c--------------------------------------------------------------------
      if(model_b.eq.0) then ! using eqdsk
         nr4=nx+4  
         if (iboundb.ge.1) then
            !the point is outside the psi>psilim 
            idx=0 !    idx derivativs order 0.ge.idx.le.3
            ffr=ias1r(txf,nx,nr4,cx,idx,psilim)
            ffd=ffr ! to define bphi
            dffr=0.d0
            dres=dffr ! to define derivs of bphi
         else
            !the point is inside the psi<psilim
            idx=0
            ffr=ias1r(txf,nx,nr4,cx,idx,res)
            ffd=ffr ! to define bphi
            idx=1
            dffr=ias1r(txf,nx,nr4,cx,idx,res)
            dres=dffr ! to define derivs of bphi
         endif
         if(r.ge.dreq) then
           pp=1.d0/r ! presumably never gets to r=0
           cos_phi=x*pp
           sin_phi=y*pp      
           dpdxd= cos_phi*dpdrd - (sin_phi*pp)*dpdph ! dpsi/dx ! Added for now
           dpdyd= sin_phi*dpdrd + (cos_phi*pp)*dpdph ! dpsi/dy ! Added for now
           bz= -dpdrd*pp*spol
           br=  dpdzd*pp*spol
           bphi=ffd*pp
           bx= br*cos_phi - bphi*sin_phi ! YuP added for now
           by= br*sin_phi + bphi*cos_phi ! YuP added for now
         else ! r=0
           dpdxd= 0.d0 ! dpsi/dx 
           dpdyd= 0.d0 ! dpsi/dy 
           !bz= -dpdrrd*spol !YuP[12-2016] not always good because of inaccuracy 
           !in second derivative of PSI. Find Bz from a small R close to 0:
           rtiny=dreq
           dpdrd0=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,1,0,rtiny,z,nr4a) ! dpsi/dr
           bz=-(dpdrd0/rtiny)*spol !YuP[12-2016] new vers. of Bz at R~0
           br=  0.d0
           bphi=0.d0
           bx=  0.d0
           by=  0.d0
         endif
      endif ! model_b=0 (eqdsk)
c----------------------------------------------------------------o
cYu.P. Analytical version for model_b=1. (No need for spline)
      if(model_b.eq.1) then ! uniform bz, no br, and bphi~1/r (if any)
        br= 0.d0
        bz= bz0      ! bz   Uniform in r 
        dbzdx= 0.d0  ! dbz/dx
        dbzdy= 0.d0  ! dbz/dy
        dbzdz= 0.d0  ! dbz/dz
        if(r.gt.dreq) then
          bphi= rbphi0/r
          bx= -rbphi0*y/r2             ! bx
          by=  rbphi0*x/r2             ! by
          dbxdx=-rbphi0*(-2.*x*y/r4)   ! dbx/dx
          dbxdy=-rbphi0*((x*x-y*y)/r4) ! dbx/dy
          dbxdz= 0.d0                  ! dbx/dz
          dbydx= rbphi0*((y*y-x*x)/r4) ! dby/dx == -dbx/dy
          dbydy= rbphi0*(-2.*x*y/r4)   ! dby/dy == -dbx/dx
          dbydz= 0.d0                  ! dby/dz
        else ! limit bphi near r=0; assume bphi~r
          bphi= rbphi0*r/dreq**2  ! Note: req(1) can be zero
          bx=  -rbphi0*y/dreq**2  ! bx 
          by=   rbphi0*x/dreq**2  ! by
          dbxdx= 0.d0               ! dbx/dx
          dbxdy=-rbphi0/dreq**2   ! dbx/dy
          dbxdz= 0.d0               ! dbx/dz
          dbydx= rbphi0/dreq**2   ! dby/dx
          dbydy= 0.d0               ! dby/dy
          dbydz= 0.d0               ! dby/dz
        endif
        dpdxd= x*bz  ! dpsi/dx
        dpdyd= y*bz  ! dpsi/dy
        dpdzd= 0.d0  ! dpsi/dz
        psid=  0.5d0*bz0*r2 ! psi  poloidal flux ! store in one.i
        goto 20 !-> define |b|, 1/|b|, derivs of |b|;  finish
      endif ! model_b.eq.1

cYu.P. Analytical version for model_b=2. (No need for spline)
      if(model_b.eq.2) then 
        ! field from coils (loop currents), and bphi~1/r (if any)
        call bfield_coils(x,y,z, bx, by, bz) !-> get bx,by,bz 
        if(R.gt.dreq) then
          pp=1.d0/R ! never gets to R=0
          cos_phi=x*pp
          sin_phi=y*pp      
          br= bx*cos_phi+by*sin_phi !was: sqrt(bx*bx+by*by) which is |br|
          ! Add toroidal field:
          bphi= rbphi0/r 
          bx= bx -rbphi0*y/r2          ! bx
          by= by +rbphi0*x/r2          ! by
        else ! limit bphi near R=0; assume bphi~R
          br=0.d0
          ! Add toroidal field:
          bphi= rbphi0*r/dreq**2  ! Note: req(1) can be zero
          bx=  bx -rbphi0*y/dreq**2  ! bx 
          by=  by +rbphi0*x/dreq**2  ! by
        endif
        ! psid  found above from spline
        dpdxd= -x*bz*spol  ! dpsi/dx 
        dpdyd= -y*bz*spol  ! dpsi/dy 
        dpdzd=  r*br*spol  ! dpsi/dz
        dpdrd= -r*bz*spol  ! dpsi/dr
        ! For now, find the derivatives of |B| field numerically,
        ! because sometimes it's difficult to define them analytically.
        ! But usually B is a smooth function of (x,y,z)
        hh=rmax*1.d-4
        odhh=0.5/hh
        call bfield_coils(x+hh,y,z, bxp, byp, bzp) !-> get bx,by,bz 
        call bfield_coils(x-hh,y,z, bxm, bym, bzm) !-> get bx,by,bz 
        dbxdx= (bxp-bxm)*odhh
        dbydx= (byp-bym)*odhh
        dbzdx= (bzp-bzm)*odhh
        call bfield_coils(x,y+hh,z, bxp, byp, bzp) !-> get bx,by,bz 
        call bfield_coils(x,y-hh,z, bxm, bym, bzm) !-> get bx,by,bz 
        dbxdy= (bxp-bxm)*odhh
        dbydy= (byp-bym)*odhh
        dbzdy= (bzp-bzm)*odhh
        call bfield_coils(x,y,z+hh, bxp, byp, bzp) !-> get bx,by,bz 
        call bfield_coils(x,y,z-hh, bxm, bym, bzm) !-> get bx,by,bz 
        dbxdz= (bxp-bxm)*odhh
        dbydz= (byp-bym)*odhh
        dbzdz= (bzp-bzm)*odhh
        goto 20 !-> define |b|, 1/|b|, derivs of |b|;  finish
      endif ! model_b.eq.2
      
cYu.P. Analytical version for model_b=3. (No need for spline)
      if(model_b.eq.3) then 
        ! model "mirror1", and bphi~1/r (if any)
        call eq_mirror1(x,y,z, PSI, bx, by, bz) !-> PSI and B
        if(r.gt.dreq) then
          pp=1.d0/r ! never gets to r=0
          cos_phi=x*pp
          sin_phi=y*pp      
          br= bx*cos_phi+by*sin_phi !was: sqrt(bx*bx+by*by) which is |br|
          ! Add toroidal field:
          bphi= rbphi0/r 
          bx= bx -rbphi0*y/r2          ! bx
          by= by +rbphi0*x/r2          ! by
        else ! limit bphi near r=0; assume bphi~r
          br=0.d0
          ! Add toroidal field:
          bphi= rbphi0*r/dreq**2  ! Note: req(1) can be zero
          bx=  bx -rbphi0*y/dreq**2  ! bx 
          by=  by +rbphi0*x/dreq**2  ! by
        endif
        ! psid  found above from spline
        dpdxd= -x*bz*spol  ! dpsi/dx 
        dpdyd= -y*bz*spol  ! dpsi/dy 
        dpdzd=  r*br*spol  ! dpsi/dz
        dpdrd= -r*bz*spol  ! dpsi/dr
        ! For now, find the derivatives of |B| field numerically,
        ! because sometimes it's difficult to define them analytically.
        ! But usually B is a smooth function of (x,y,z)
        hh=rmax*1.d-4
        odhh=0.5/hh
        call eq_mirror1(x+hh,y,z,PSI, bxp, byp, bzp) !-> get bx,by,bz 
        call eq_mirror1(x-hh,y,z,PSI, bxm, bym, bzm) !-> get bx,by,bz 
        dbxdx= (bxp-bxm)*odhh
        dbydx= (byp-bym)*odhh
        dbzdx= (bzp-bzm)*odhh
        call eq_mirror1(x,y+hh,z,PSI, bxp, byp, bzp) !-> get bx,by,bz 
        call eq_mirror1(x,y-hh,z,PSI, bxm, bym, bzm) !-> get bx,by,bz 
        dbxdy= (bxp-bxm)*odhh
        dbydy= (byp-bym)*odhh
        dbzdy= (bzp-bzm)*odhh
        call eq_mirror1(x,y,z+hh,PSI, bxp, byp, bzp) !-> get bx,by,bz 
        call eq_mirror1(x,y,z-hh,PSI, bxm, bym, bzm) !-> get bx,by,bz 
        dbxdz= (bxp-bxm)*odhh
        dbydz= (byp-bym)*odhh
        dbzdz= (bzp-bzm)*odhh
        goto 20 !-> define |b|, 1/|b|, derivs of |b|;  finish
      endif ! model_b.eq.3

cYu.P. Analytical version for model_b=4. (No need for spline)
      if(model_b.eq.4) then 
        ! FRC-like plasma with Bz-field only, plus bphi~1/r (if any)
        call bfield_frc(x,y,z, bx, by, bz) !-> get bx,by,bz 
        if(r.gt.dreq) then
          pp=1.d0/r ! never gets to r=0
          cos_phi=x*pp
          sin_phi=y*pp      
          br= bx*cos_phi+by*sin_phi ! Note: bx=0, by=0 for model_b=4
          dbzdz= 0.d0  ! dbz/dz
          dbzdr= bz0*akappa*4*r /
     /        ( rs_frc*cosh( akappa*(2.*r2/rs_frc**2 -1.d0) ) )**2
          dbzdx= dbzdr*cos_phi
          dbzdy= dbzdr*sin_phi
          ! Add toroidal field (~ 1/r):
          bphi= rbphi0/r 
          bx= bx -rbphi0*y/r2          ! bx
          by= by +rbphi0*x/r2          ! by
          dbxdx=-rbphi0*(-2.*x*y/r4)   ! dbx/dx
          dbxdy=-rbphi0*((x*x-y*y)/r4) ! dbx/dy
          dbxdz= 0.d0                  ! dbx/dz
          dbydx= rbphi0*((y*y-x*x)/r4) ! dby/dx == -dbx/dy
          dbydy= rbphi0*(-2.*x*y/r4)   ! dby/dy == -dbx/dx
          dbydz= 0.d0                  ! dby/dz
        else ! limit near r=0
          br=0.d0 ! For FRC-field, this model
          dbzdz= 0.d0  ! dbz/dz
          dbzdr= 0.d0
          dbzdx= 0.d0
          dbzdy= 0.d0
          ! Add toroidal field: ( limit bphi near r=0: assume bphi~r )
          bphi= rbphi0*r/dreq**2  ! Note: req(1) can be zero
          bx=  bx -rbphi0*y/dreq**2  ! bx 
          by=  by +rbphi0*x/dreq**2  ! by
          dbxdx= 0.d0               ! dbx/dx
          dbxdy=-rbphi0/dreq**2   ! dbx/dy
          dbxdz= 0.d0               ! dbx/dz
          dbydx= rbphi0/dreq**2   ! dby/dx
          dbydy= 0.d0               ! dby/dy
          dbydz= 0.d0               ! dby/dz
        endif
        dpdxd= -x*bz*spol  ! dpsi/dx 
        dpdyd= -y*bz*spol  ! dpsi/dy 
        dpdzd=  r*br*spol  ! dpsi/dz
        dpdrd= -r*bz*spol  ! dpsi/dr
        goto 20 !-> define |b|, 1/|b|, derivs of |b|;  finish
      endif ! model_b.eq.4 (FRC-like plasma)

c------------------------------------------------------------------
c this  part of the program calculates rho_magnetic (based on pol.flux)
c for the functions dense,dxdr,dxdz
cyup       if (iboundb.eq.-1) rho=rhopsi(res)  ! in bxyz
c------------------------------------------------------------------
c     \vector{B}=B_phi+\vector{B_pol}
c     B_phi=(Feqd(psi)/R)\hat{phi}
c     \vector{B_pol}=(1/R)*grad(psi) x \hat{phi}
c     here x is the vector product
c     then  B_z=(1/R)*dpsi/dr and B_r=-(1/R)*dpsi/dz  (*).
c     If the total current is not given, then the sign of the B_pol
c     determined from psi, using the previous formula (*).
c     If the current 'toteqd'is given, then the sign of it will indicate 
c     the direction of B_pol, regardless of the convention for psi(R,Z)
c     toteqd is positive if it is directed along \hat{phi}.
c-------------------------------------------------------------------
c     The poloidal flux	array peqd(i,j) was transformed in equilib.f
c     in subroutine input. It was multiplied by dpsimax.
c     If (psimag.gt.psilim) then dpsimax=-1.d0 else  dpsimax=1.d0
c     So after this multiplication poloidal flux peqd(i,j)
c     always has the minimum on the magnetic axis.
c-------------------------------------------------------------------
c     Derivatives of magnetic field - needed for analytical dD/dr...
c     (Not ready yet)
      goto 20 ! skip for now
c------------------------------------------------------------------
      dbzdz=-pp*dpdzrd
      dbzdr=pp*pp*dpdrd-pp*dpdrrd
      dbzdph=0.d0
      
      dbrdz=pp*dpdzzd
      dbrdr=-pp*pp*dpdzd+pp*dpdzrd
      dbrdph=0.d0
      
      dbphdz=pp*dres*dpdzd
      dbphdr=-pp*bphi+pp*dres*dpdrd
      dbphdph=0.d0
      
      dbzdz=dbzdz*spol
      dbzdr=dbzdr*spol
      dbrdz=dbrdz*spol
      dbrdr=dbrdr*spol      
      
      dbxdx= cos_phi*(cos_phi*dbrdr -sin_phi*dbphdr) ! YuP added for now
     +      -sin_phi*(cos_phi*dbrdph-sin_phi*dbphdph
     +               -sin_phi*br    -cos_phi*bphi    )*pp
     
      dbxdy= sin_phi*(cos_phi*dbrdr -sin_phi*dbphdr) ! YuP added for now
     +      +cos_phi*(cos_phi*dbrdph-sin_phi*dbphdph
     +               -sin_phi*br    -cos_phi*bphi    )*pp
      
      dbxdz= dbrdz*cos_phi - dbphdz*sin_phi ! YuP added for now
      
      dbydx= cos_phi*(sin_phi*dbrdr +cos_phi*dbphdr) ! YuP added for now
     +      -sin_phi*(sin_phi*dbrdph+cos_phi*dbphdph
     +               +cos_phi*br    -sin_phi*bphi    )*pp
     
      dbydy= sin_phi*(sin_phi*dbrdr +cos_phi*dbphdr) ! YuP added for now
     +      +cos_phi*(sin_phi*dbrdph+cos_phi*dbphdph
     +               +cos_phi*br    -sin_phi*bphi    )*pp

      dbydz= dbrdz*sin_phi + dbphdz*cos_phi ! YuP added for now
      dbzdx= cos_phi*dbzdr - sin_phi*pp*dbzdph ! YuP added for now
      dbzdy= sin_phi*dbzdr + cos_phi*pp*dbzdph ! YuP added for now
c------------------------------------------------------------------
 20   continue
c------------------------------------------------------------------
      bxyz= dsqrt(bx*bx + by*by + bz*bz)
      ob=   1.d0/bxyz
      o_bmod=ob !-> store in one.i
      dbmdx= ob*(bx*dbxdx + by*dbydx + bz*dbzdx) 
      dbmdy= ob*(bx*dbxdy + by*dbydy + bz*dbzdy) 
      dbmdz= ob*(bx*dbxdz + by*dbydz + bz*dbzdz) 
c---------------------------------------------------------------------
      return
      end


c======================================================================
c======================================================================


c        ************************ rhof_xyz(x,y,z)  **************
c        *                        -                           *
c        * this function calculates small radius rho          *
c        * it controls if the ray point is inside the plasma;  *
c        * if the ray point is outside the plasma then     *
c        * it chooses the specific value for rho	      *
c        ******************************************************
c
c-------------------------------------------------------------
c      input:				     !
c      x,y,z  cartesian coords.    
c------------------------------------------------------------
      real*8
     1function rhof_xyz(x,y,z)  ! needs work
      implicit none
c      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
      include 'six.i'
c-----input
      real*8 x,y,z
c-----externals
      real*8
     1 ias1r,ias2r_Sm,bxyz,rhopsi,
     & thetapol
c-----locals
      real*8 epsbnd,rrr,zzrp,zp,zzrm,zm,res,delr,delz,
     & r, theta_l, r_b, z_b,
     & r_l, z_l,  l_zr, l_z_b_r_b
      integer idx,ipx,ipx4,nr4,nz4
      
      stop 'rhof_xyz  is not used yet'

      bmod=bxyz(x,y,z)     
c--------------------------------------------------------------
c     for checking that the point is inside the plasma
      epsbnd=2.d-3 ![m]
      iboundb=-1
c-----------------------
      r= dsqrt(x*x+y*y)
      r_l=r-rma
      z_l=z-zma

      if ((r_l**2+z_l**2).lt.1.d-10) then
         goto 10
      endif

      call theta_rz(r_l,z_l,theta_l) ! get pol.angle 0=< theta <2*pi

c-----calculate coordinates r_b,z_b at the limiter surface psi=psilim
c     and the poloidal angle thetapol
      call zr_psith(psilim,theta_l, z_b,r_b)

c      write(*,*)'spldens.g in psilim,z_b,r_b',psilim,z_b,r_b

      l_zr=dsqrt(r_l**2+z_l**2)
      l_z_b_r_b=dsqrt((r_b-rma)**2+(z_b-zma)**2)
      if (l_zr.ge.l_z_b_r_b) then
         rhof_xyz= l_zr / l_z_b_r_b
         iboundb=3
      endif
      goto 10
c-----------------------
      if ((r.lt.rmin).or.(r.ge.rmax-epsbnd)) then
         iboundb=2
         theta_l=thetapol(z,r) 
         call zr_psith(psilim,theta_l,z_b,r_b)
         rhof_xyz
     +    =dsqrt(((r-rma)**2+(z-zma)**2)/((r_b-rma)**2+(z_b-zma)**2))
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
         theta_l=thetapol(z,r) 
         call zr_psith(psilim,theta_l,z_b,r_b)
         rhof_xyz=
     +     dsqrt(((r-rma)**2+(z-zma)**2)/((r_b-rma)**2+(z_b-zma)**2))
       
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


      if (iboundb.eq.-1) rhof_xyz=rhopsi(res)

      rho=rhof_xyz ! the data for common block one.i
                          
      return
      end


c======================================================================
c======================================================================
c        **********************   wcw   *********************
c        *                        -                         *
c        * this function calculates the ratio of  cyclotron *
c        * frequency of species i and wave frequency        *
c        * wcw= omega_cyclotron_i/omega                        *
c        ****************************************************
c
c-------------------------------------------------------------------
c								    
c        input parameters					    
c      								    
c      x,y,z - cartesian coordinates of the  point    
c              where  the  ratio  is calculated.      		         	    
c-------------------------------------------------------------------
      double precision
     1function wcw(x,y,z,i)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i' ! contains w(i), bmod
      !--------------
c     w0=28.0*b0/frqncy            ! set in dinit.f;  b0=1.d0 (normalization)
c     w(i)= w0*charge(i)/dmas(i)   ! set in dinit.f
      isp=min(i,nbulk) ! not to exceed nbulk
      wcw=w(isp)*bmod
      if(dabs(wcw).lt.1.d-8) wcw=1.d-8
      return
      end
      
      
c======================================================================
c======================================================================
c        **********************   wpw_2   *******************
c        *                        -                         *
c        * this function calculates the ratio of squares of *
c        * plasma frequency of species i and wave frequency *
c        * wpw_2= (omega_pl_i/omega)**2    
c        * The minimal wpw_2 value is set inside this function
c        * wpw_2_min=1.d-6. If wpw_2<wpw_2_min then wpw_2=wpw_2_min                        
c        ****************************************************
c
c-------------------------------------------------------------------
c								    
c        input parameters					    
c      								    
c      x,y,z - cartesian coordinates of the  point    
c              where  the  ratio  is calculated.      		         	    
c      small radius rho is in common one
c-------------------------------------------------------------------
      double precision
     1function wpw_2(x,y,z,i)
      implicit integer (i-n), real*8 (a-h,o-z)
c      implicit none
      include 'param.i'
      include 'one.i' ! contains v(i), rho, etc.
c-----input
      real*8 x,y,z ! cartesian space coordinates
c     rho is in common  /one/
      integer i,isp !number of plasma species
c-----local 
      real*8 den
c-----external
      real*8 dense_xyz      
      !dense_xyz finds rho, stores in one.i; rho is set by model_rho_dens
      !write(*,*)'wpw_2: i,v(i)=',i,v(i)
      isp=min(i,nbulk) ! not to exceed nbulk
      den= dense_xyz(x,y,z,isp) 
c     v0=806.2/(frqncy*frqncy)     ! set in dinit.f
c     v(i)=v0*charge(i)**2/dmas(i) ! set in dinit.f
      wpw_2= v(isp)*den
      wpw_2= dmax1(wpw_2,1.d-6)      
      return
      end
c======================================================================
c======================================================================


c        *********************dense_xyz ***********************
c        *                        -                           
c        * this function calculates the density profile       
c        * as a function of x,y,z (mostly as a function of rho)
c        *  rho= function of poloidal flux (model_rho_dens=0 or 5)
c        *     or model-based surfaces (model_rho_dens=1-4)
c        *  rho is saved into common 'one.i'
c        *                                                    
c        ******************************************************
c
c------------------------------------------------------------------
c								   
c        input parameters					   
c      								   
c       x,y,z - cartesian coordinates of the  point    
c               where  the  density  is calculated.      		         	    
c       rho is saved to common one.i                                    
c------------------------------------------------------------------
c      uses
c         constants dense0,denseb,rn1de,rn2de,idens are in common 'one'
c         For model_rho_dens=1,2; and 4 (partially):
c         elx0 ely0 elz0  elax elay elaz 
c         For model_rho_dens=2; and 4 (partially):
c         dens0rr dens0es dens0ub eltheta Rm0rr Rm0es rtau r_ub_edge
c         The values are from the file genray.in (reading: dinit)
c      functions
c         _densrho_(rho)    finds dense_xyz from spline
c------------------------------------------------------------------
      real*8 function dense_xyz(x,y,z,i)
c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none
      include 'param.i'
      include 'one.i' ! stores rho
      include 'three.i'
      include 'fourb.i'
      include 'five.i'
c-----input
      real*8 x,y,z ! cartesian space coordinates
      integer i !plasma species number
c-----externals 
      real*8 densrho,rho_dens_xyz,ias2r_Sm
c-----local
      real*8 arg,sech,rho2,Rs0rr,Rs,r,r2,rrk,dens_ub,dens_rr,dens_es
      real*8 gn, zabs, zj
      real*8 rn,zn,  da11,da12,da21,da22, drn,dzn
      real*8 dreq, dzeq, req_mn,req_mx, zeq_mn,zeq_mx
      integer ir,ir1,ir2, iz,iz1,iz2
      
      sech(arg)= 1.d0/cosh(arg) !2.d0/(exp(arg)+exp(-arg)) ! define sech
      
      if(i.gt.nbulk) stop 'dense_xyz: i>nbulk violation'
      
      if(model_rho_dens.eq.0) then
      
         rho=rho_dens_xyz(x,y,z) 
         rho2=rho*rho
         
         dense_xyz= densrho(rho,i) ! uses spline(rho)
         ! Exponential drop at rho>rho095 is already applied in densrho
         
         if( zlength_den_dropoff(i).ne.0.d0) then
         ! Apply additional drop in z-axis: gn(z) factor
         zabs=abs(z)
         if( zabs.gt.zbegin_den_dropoff(i) ) then
          gn=
     +    exp(-((zabs-zbegin_den_dropoff(i))/zlength_den_dropoff(i))**2)
          dense_xyz= dense_xyz*gn !set a density drop-off in z-coord
         endif
         endif

      elseif(model_rho_dens.eq.1) then
         rho=rho_dens_xyz(x,y,z) 
         rho2=rho*rho
         dense_xyz= (dense0(i)-denseb(i))*den_scale(i)*
     1	                  (1-rho**rn1de(i))**rn2de(i)+denseb(i)

         if ((izeff.eq.0).or.(izeff.eq.3)) then
           write(*,*)'dense_xyz: model_rho_dens=',model_rho_dens
           stop 'dense_xyz: Not setup for izeff=0 or 3'
         endif

      elseif(model_rho_dens.eq.2) then
         rho=rho_dens_xyz(x,y,z) 
         rho2=rho*rho
         !-1-> Rigid Rotor profile:
         Rs0rr= max(elax,elay)
         rrk=(Rm0rr/Rs0rr)**2  ! K in Eq.(38)
         Rs=Rs0rr ! Assume no dependence in z, for now
         dens_rr= dens0rr*( sech(rho2-rrk)/sech(rrk) )**2
         ! Note: At rho=0 -> nrr= n0rr
         ! Note: At rho=1 -> nrr= n0rr*[sech(1-K)/sech(K)]^2
         !-2-> Ellipsoidal Spindle profile:
         dens_es= dens0es* exp( -(Rs*rho-Rm0es)**2 / (2*rtau**2) )
         !-3-> Uniform background density profile:
         dens_ub= dens0ub 
         r= sqrt(x*x+y*y)
         if (r .gt. r_ub_edge) then 
            ! Linear drop in region  r_ub_edge < r < wall_rmax
            dens_ub= dens0ub*(1.d0 -(r-r_ub_edge)/(wall_rmax-r_ub_edge))
         endif
         ! Set density to zero outside of chamber wall radius wall_rmax:
c         if (r .gt. wall_rmax) then
c            dens_rr=0.d0    ! Skip this option 
c            dens_es=0.d0    !  because it creates a jump at wall_rmax
c            dens_ub=0.d0    !  and affects accuracy during reflection
c         endif
         dense_xyz= dens_rr + dens_es + dens_ub !electrons only, for now

      elseif(model_rho_dens.eq.3) then 
         ! (x,y)-spline based on dengrid(i,j) data from file.
         ! Note: in this model, dn/dz=0; uniform density in z-direction.
         dense_xyz=ias2r_Sm(tx_den,nx_den,ty_den,ny_den,
     +                  cxy_den,ncx_den,ncy_den,0,0,x,y,nx4a)
         ! Define rho here; needed for boundary control, etc.
         rho=1.d0-(dense_xyz-denmin)/(denmax-denmin) !-> to one.i
         ! With such definition, rho=0 at n(x,y)=denmax, 
         ! and rho=1 at plasma edge 
         ! (which is the set of points with n(x,y)=denmin,
         ! where denmin is the smallest positive density on the grid)

      elseif(model_rho_dens.eq.4) then !FRC-like plasma, uniform Bz in Z
         r2=x*x+y*y
         r= sqrt(r2)
         rho= rho_dens_xyz(x,y,z)  ! rho=abs(2.*(r/rs_frc)^2 -1.d0) 
         ! Note: rho=0 at r=rs_frs/sqrt(2);  rho=1 at r=0 or r=rs_frc
         !-1-> Rigid Rotor profile:
         dens_rr= dens0rr*( sech(akappa*rho) )**2
         ! Note: At r=0 and r=rs_frc,   nrr= n0rr*[sech(akappa)]^2
         !-2-> Uniform background density profile:
         dens_ub= dens0ub 
         if (r .gt. rs_frc) then 
            ! Linear drop in region  rs_frc < r < wall_rmax
            dens_ub= dens0ub*(1.d0 -(r-rs_frc)/(wall_rmax-rs_frc))
         endif
         dense_xyz= dens_rr + dens_ub !electrons only, for now
         
      elseif(model_rho_dens.eq.5) then !TAE-FRC, model_b=0 (eqdsk)
         if(eqdsktype.ne.'TAE') stop
         ! Only for eqdsktype='TAE'.
         ! Profile of ne is in eqdsk file, 
         ! over same grid as B or PSI==peqd(nreqd,nzeqd).
         ! necut(nreqd,nzeqd) [10^19 /m^3] is stored in fourb.i
         rho=rho_dens_xyz(x,y,z) ! rho is same as for model_rho_dens=0
         rho2=rho*rho
         r= sqrt(x*x+y*y) ! major R radius
         zj=z
         ! Bilinear interpolation from four nearest nodes:
         req_mn=req(1)
         req_mx=req(nreqd)
         dreq=req(2)-req(1)
         zeq_mn=zeq(1)
         zeq_mx=zeq(nzeqd)
         dzeq=zeq(2)-zeq(1)
         ! Make sure (r,z) point is within (req,zeq):
         r= max(r, req_mn ) 
         r= min(r, req_mx-dreq) !stay away from upper edge
         zj= max(zj, zeq_mn ) 
         zj= min(zj, zeq_mx-dzeq) !stay away from upper edge
         ! Find the cell in (req,zeq) that contains (r,z):
         rn= (r-req_mn)/dreq  ! normalized R coordinate
         zn= (zj-zeq_mn)/dzeq  ! normalized Z coordinate
         ir=  rn    ! lower integer: can be from 0 to NRgrid-2
         ir1= ir+1
         ir2= ir+2  ! rn is between ir1 and ir2 nodes of req
         iz=  zn    ! lower integer: can be from 0 to NZgrid-2
         iz1= iz+1
         iz2= iz+2  ! zn is between iz1 and iz2 nodes of zeq
         ! In normalized units:
         ! distance between (r,z) point and left-bottom-corner:
         drn= rn-ir  
         dzn= zn-iz 
         da11 = necut(IR1,IZ1)
         da21 = necut(IR2,IZ1) - da11
         da12 = necut(IR1,IZ2) - da11
         da22 = necut(IR2,IZ2) - necut(IR1,IZ2) - da21
         dense_xyz= da11 + drn*da21 + dzn*(da12 + drn*da22) 
         !Correction for edge (important for rho>1). Units 10^19 m^-3
         dense_xyz= dense_xyz +0.05*exp(-(rho-1.d0)**2) 
         ! exp() has a peak at rho=1, and exp drop to both sides
         ! from rho=1. (For rho<1, this adjustment is not essential
         ! because the density is large enough at rho<1.)
         ! YuP[Nov-2014] Done

      else
        PRINT*,'Set model_rho_dens to 0 or 1 or 2 or 3 or 4 or 5'
        PRINT*,'model_rho_dens=1,2,3,4 works with ixyz=1 only'
        PRINT*,'model_rho_dens=0 works with ixyz=0 or 1'
        PRINT*,'model_rho_dens=5 only for model_b=0, eqdsktype=TAE'
        stop 'dense_xyz'      
      endif
      
      dense_xyz=max(1.d-5,dense_xyz) !Set Lower limit. Units 10^19 m^-3
      
      return
      end ! dense_xyz
c======================================================================
c======================================================================


c        *********************tempe_xyz ***********************
c        *                        -                           *
c        * this function calculates the electron temperature  *
c        * profile as a function  of rho       	              *
c        *  rho= function of poloidal flux (model_rho_dens=0 or 5)
c        *     or model-based surfaces (model_rho_dens=1-4)
c        * rho is in common 'one'                             *
c        ******************************************************
c
c------------------------------------------------------------------
c								   
c        input parameters					   
c      								   
c      x,y,z - cartesian coordinates of the  point    
c              where  the  temperature  is calculated.      		         	    
c      rho is from common one.i
c------------------------------------------------------------------
	double precision
     1function tempe_xyz(x,y,z,i)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i' ! stores rho
      !--------------
      ! get rho based on model for density (model_rho_dens=1,2,4), 
      ! 2D spline of data (model_rho_dens=3),
      ! or based on magnetic flux (model_rho_dens=0,5):
      den=dense_xyz(x,y,z,1) !-> get rho (stored in one.i) 
      !--------------
	tempe_xyz=temperho(rho,i)
      return
      end
c======================================================================
c======================================================================


c        ********************** gamma1_xyz ***********************
c        *                      -----                        *
c        * this function calculates the angle  between  the  *
c        * wave vector and the magnetic field and the        *
c        * derivatives from cos^2(gamma) by x,y,z,Nx,Ny,Nz   *
c        *****************************************************
c
c-------------------------------------------------------------------
c								    
c        input parameters					    
c      								    
c      x,y,z - cartesian coordinates of the  point    
c              where  the  angle  is calculated.      		         	    
c								    
c      cn_x,cn_y,cn_z are N_x, N_y, N_z components of a wave 
c                     refractive index.                           
c-------------------------------------------------------------------
      double precision
     1function gamma1_xyz(x,y,z, cn_x,cn_y,cn_z)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      !--------------
      cnt2= cn_x*cn_x + cn_y*cn_y + cn_z*cn_z
      cnt= dsqrt(cnt2)    ! == |N|
      gg=  cn_x*bx + cn_y*by + cn_z*bz ! == N.b
      cnbm= cnt*bmod
      arg= gg/cnbm         ! == N.b/(|N||b|)
      if(dabs(arg).ge.1.d0) then
          if  (abs(arg).gt.1.d+100) arg=sign(1.d+100,arg)
          if(arg.ge.1.d0)then 
            arg=1.d0
            dc=1.d0
            gamma1_xyz=0.d0
          else
            arg=-1.d0
            dc=-1.d0
            pi=4.d0*datan(1.d0)
            gamma1_xyz=pi
          endif
      else ! |arg|<1.0 
          gamma1_xyz=dacos(arg)
          dc=dcos(gamma1_xyz)
      endif
       

      dc2dx= 2.d0*dc*((cn_z*dbzdx+cn_x*dbxdx+cn_y*dbydx)/cnbm-
     1               dc*dbmdx*o_bmod) !== d[cos^2(gamma)]/dx

      dc2dy= 2.d0*dc*((cn_z*dbzdy+cn_x*dbxdy+cn_y*dbydy)/cnbm-
     1               dc*dbmdy*o_bmod) !== d[cos^2(gamma)]/dy

      dc2dz= 2.d0*dc*((cn_z*dbzdz+cn_x*dbxdz+cn_y*dbydz)/cnbm-
     1               dc*dbmdz*o_bmod) !== d[cos^2(gamma)]/dz
     
      dc2dnx=2.d0*dc*(bx/cnbm-dc*cn_x/cnt2) !== d[cos^2(gamma)]/dNx
      dc2dny=2.d0*dc*(by/cnbm-dc*cn_y/cnt2) !== d[cos^2(gamma)]/dNy
      dc2dnz=2.d0*dc*(bz/cnbm-dc*cn_z/cnt2) !== d[cos^2(gamma)]/dNz

      return
      end
      
c======================================================================
c======================================================================

c-----------------------------------------------------------------------
c     this subroutine DNONETWO_xyz calculates power and current density
c     profiles on rho. Arrays: powden(NR) (erg/(cm**3*c)
c                           and currden(NR).
c-----------------------------------------------------------------------
      subroutine dnonetwo_xyz
!      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none
      include 'param.i'
      include 'one.i'
      include 'onetwo.i'
      include 'gr.i'
      include 'rho.i'
      include 'three.i'
      include 'five.i'
c-----externals
      real*8 rhov,rhos,rho_lrho,psi_rho,ias1r, bxyz
c-----locals
      real*8 pwtot_s
      dimension pwtot_s(nbulka)  !Added for indiv ion contrib[BH041009]
      real*8 rholeft,rhoright,binplarea,psi,f_eqd,
     &b_av,bs_av,r_av,dr_av,drs_av,theta,r_bmin,x_bmin,y_bmin,z_bmin,
     + b_pol_bmin,
     & powertot,rho0,rho0_pol,poloidlen,
     &powrtot1,currtot1,pwtot_e,pwtot_i,pwtot_cl,
     &rho_l,psi_l  
      integer i,kk,k,idx,nr4,j
     
      pi=4.0d0*datan(1.0d0)
c-----------------------------------------------------------------------
c----------------------------------------------------------------------
c     spower (erg/sec), scurrent(A), binvol(cm**3),binarea(cm**2)
c     powden(erg/(sec*cm**3)),currden(A/cm**2)
c     voltot (m**3), areatot (m**2)
c----------------------------------------------------------------------
      do i=1,NR-1
         powden(i)=spower(i)/binvol(i)
         powden_e(i)=spower_e(i)/binvol(i)
         powden_i(i)=spower_i(i)/binvol(i)
         !YuP if(iabsorp.eq.3) then ! powden_s
            do kk=1,nbulk !YuP was: 2,nbulk
               powden_s(i,kk)=spower_s(i,kk)/binvol(i)
            enddo
         !YuP endif
         powden_cl(i)=spower_cl(i)/binvol(i)
cyup         currden(i)=scurrent(i)/binarea(i) !toroidal current density
      enddo  ! i=1,NR-1

c      do k=2,nbulk
c        write(*,*)'k=',k,'powden_s'
c        write(*,*)(powden_s(i,k), i=1,NR)
c      enddo

c-----------------------------------------------------------
c     CD calculation using GA memo
c-----------------------------------------------------------
c     cur_den_onetwo=<j_parallel.B>/B_0=<j_parallel><B**2>/<B>/B_0
c----------------------------------------------------------
      idx=0 
      nr4=nx+4
 
cyup      write(*,1010)'cur_den_tor=p_tor*cur_den_par, '
cyup     &   // 'p_tor = drs_av*f_eqd/(b_av*dr_av)'
cyup      write(*,1010)'cur_den_pol=p_pol*cur_den_par, '
cyup     &   // 'p_pol = b_pol_bmin/b_av'

cyup      write(*,1010)'i rho_bin_center powden_e   cur_den_par  p_tor'
cyup     & // '      cur_den_tor'
cyup     & // '   p_pol    cur_den_pol'
    
 1010    format(/,1x,a)
 1011    format(i3,7(1pe12.4))

      if(model_b.eq.0) then ! eqdsk file, tokamak geom.
      
      do i=1,NR-1
        rho_l= rho_bin_center(i) ! YuP: was: (i-0.5d0)/(NR-1)   
        psi_l=psi_rho(rho_l)
        !write(*,*)'i,rho_l,psi_l',i,rho_l,psi_l
cSAP080321 argument psi in f_eqd was changed from psi to psi_l
        f_eqd=ias1r(txf,nx,nr4,cx,idx,psi_l)

        call average_variables_xyz(psi_l,b_av,bs_av,r_av,dr_av,drs_av)
        
        if(beqd.ne.0.d0) then ! Bphi(R0) is not zero
          s_cur_den_onetwo(i)= s_cur_den_parallel(i)*bs_av/(b_av*beqd)
        endif
        
        s_cur_den_toroidal(i)= s_cur_den_parallel(i)*drs_av*f_eqd/
     &                          (b_av*dr_av)

        theta=0.d0
        call zr_psith(psi_l,theta,z_bmin,r_bmin)
        x_bmin=r_bmin  ! needs work ? or not important?
        y_bmin=0.d0
        bmod= bxyz(x_bmin,y_bmin,z_bmin) !-> get bz
        b_pol_bmin=dabs(bz) ! poloidal B at the point
                                      ! with minimal B
                                      ! with theta poloidal =0


        s_cur_den_poloidal(i) = s_cur_den_parallel(i)*b_pol_bmin/
     &                          b_av
c----------------------------------------------------------------
c         write(*,1011)i,rho_l,powden_e(i),s_cur_den_parallel(i),
c     &   drs_av*f_eqd/(b_av*dr_av),s_cur_den_toroidal(i),
c     &   b_pol_bmin/b_av,s_cur_den_poloidal(i)
      enddo 
      
      
c-----------------------------------------------------------
c     powertot (erg/sec), currtot(A)
c-----------------------------------------------------------     
      powertot=0.0d0
      powtot_e=0.0
      powtot_i=0.0d0
      !YuP if (iabsorp.eq.3) then ! powtot_s
         do kk=1,nbulk !YuP was: 2,nbulk
            powtot_s(kk)=0.0d0
         enddo
      !YuP endif
      powtot_cl=0.0d0
      currtot=0.0d0
   
      do i=1,NR-1
c     integration formulas ****INT1**
         powertot=powertot+powden(i)*binvol(i)
         powtot_e=powtot_e+powden_e(i)*binvol(i)
         powtot_i=powtot_i+powden_i(i)*binvol(i)
         !YuP if (iabsorp.eq.3) then ! powden_s
            do kk=1,nbulk !YuP was: 2,nbulk
               powtot_s(kk)=powtot_s(kk)+powden_s(i,kk)*binvol(i)
            enddo
         !YuP endif
         powtot_cl=powtot_cl+powden_cl(i)*binvol(i)
         currtot=currtot+currden(i)*binarea(i)
      enddo
     
      parallel_cur_total=0.d0
      toroidal_cur_total=0.d0
      poloidal_cur_total=0.d0
      do j=1,NR-1
cyup        parallel_cur_total=parallel_cur_total+
cyup     &                     s_cur_den_parallel(j)*binarea(j)
cyup        toroidal_cur_total=toroidal_cur_total+
cyup     &                     s_cur_den_toroidal(j)*binarea(j)
cyup        poloidal_cur_total=poloidal_cur_total+
cyup     &                     s_cur_den_poloidal(j)*binarea_pol(j)
      enddo


cyup      write(*,*)'parallel_cur_total, toroidal_cur_total ',
cyup     &parallel_cur_total, toroidal_cur_total 
cyup      write(*,*)'poloidal_cur_total ',poloidal_cur_total

cyup      write(*,*)'INT1 testing DNONETWO. powertot=erg/sec',powertot,
cyup     1' currtot=A',currtot
 
cyup      write(*,*)'testing 1 DNONETWO powtot_e,powtot_i,powtot_cl',
cyup     &     powtot_e,powtot_i,powtot_cl
cyup      write(*,*)'powtot_s(1:nbulk)=',(powtot_s(kk),kk=1,nbulk)

      endif ! (model_b.eq.0) then ! eqdsk file, tokamak geom.

      powrtot1=0.d0
      currtot1=0.d0      
      pwtot_e=0.d0
      pwtot_i=0.d0
      !YuP if (iabsorp.eq.3) then ! pwtot_s
         do kk=1,nbulk !YuP was:2,nbulk
            pwtot_s(kk)=0.d0
         enddo
      !YuP endif
      pwtot_cl=0.d0
      do i=1,NR
         powrtot1=powrtot1+spower(i)
         pwtot_e=pwtot_e+spower_e(i)
         pwtot_i=pwtot_i+spower_i(i)
         !YuP if (iabsorp.eq.3) then ! spower_s pwtot_s
            do kk=1,nbulk !YuP was: 2,nbulk
               pwtot_s(kk)=pwtot_s(kk)+spower_s(i,kk)
            enddo
         !YuP endif
         pwtot_cl=pwtot_cl+spower_cl(i)
         currtot1=currtot1+scurrent(i) !total toroidal current        
      enddo
cyup      write(*,*)'INT2 testing DNONETWO. powrtot1=erg/sec',powrtot1,
cyup     1' currtot1=A',currtot1
cyup      write(*,*)'testing 2 DNONETWO pwtot_e,pwtot_i,pwtot_cl',
cyup     &     pwtot_e,pwtot_i,pwtot_cl
cyup      write(*,*)'pwtot_s(1:nbulk)=',(pwtot_s(kk),kk=1,nbulk)
cyup      write(*,*)'powden: ',powden
cyup      write(*,*)'powden_e: ',powden_e
cyup      write(*,*)'powden_i: ',powden_i
cyup      do kk=2,nbulk
cyup         write(*,*)'kk,powden_s(,kk): ',kk,(powden_s(i,kk),i=1,NR-1)
cyup      enddo

      if(model_b.ne.0) then ! model b-field, open surfaces.         
         write(*,*)'dnonetwo_xyz: binvol is not defined for model_b>0'
         do i=1,NR
            powertot=powertot+spower(i)
            powtot_e=powtot_e+spower_e(i)
            powtot_i=powtot_i+spower_i(i)
            !YuP if (iabsorp.eq.3) then ! spower_s
               do kk=1,nbulk !YuP was: 2,nbulk
                  powtot_s(kk)=powtot_s(kk)+spower_s(i,kk)
               enddo
            !YuP endif
            powtot_cl=powtot_cl+spower_cl(i)
         enddo
      endif

c -----------------------------------------------------------------
c     normalization of density profiles
c     it dependents on numerical integration formulas
c     In the case  ***INT**** we have the following normalization:
      do i=1,NR
         if (powertot.ne.0.d0) then
             powden(i)=powden(i)*powrtot1/powertot
         endif
         if (powtot_e.ne.0.d0) then
             powden_e(i)=powden_e(i)*pwtot_e/powtot_e
         endif
         if (powtot_i.ne.0.d0) then
             powden_i(i)=powden_i(i)*pwtot_i/powtot_i
         endif
         !YuP if (iabsorp.eq.3) then ! powden_s
         do kk=1,nbulk  ! YuP was: 2,nbulk
            if (powtot_s(kk).ne.0.d0) then
               powden_s(i,kk)=powden_s(i,kk)*pwtot_s(kk)/powtot_s(kk)
            endif
         enddo
         !YuP endif
         if (powtot_cl.ne.0.d0) then
             powden_cl(i)=powden_cl(i)*pwtot_cl/powtot_cl
         endif
         if (currtot.ne.0.d0) then
cyup             currden(i)=currden(i)*currtot1/currtot
         endif        
      enddo
c------------------------------------------------------------------

      write(*,998)
 998  format(//)
cSAP090306
c      write(*,999)powtott
c 999  format(' total injected power(erg/sec)   =',1pe14.6)
      write(*,999)powtott*1.d-7
 999  format(' total injected power(watt)   =',1pe14.6)
      write(*,1005)i_total_bad_initial_conditions
cSAP090306  
c      write(*,1006)power_launched
c      write(*,1000)powrtot1
c 1000 format(' total absorbed power(erg/sec)   =',1pe14.6)
c      write(*,1001)pwtot_e
c 1001 format(' total absorbed power_e(erg/sec) =',1pe14.6)
c      write(*,1002)pwtot_i
c 1002 format(' total absorbed power_i(erg/sec) =',1pe14.6)
      write(*,1006)power_launched*1.d-7
      write(*,1000)powrtot1*1.d-7
 1000 format(' total absorbed power(watt)   =',1pe14.6)
      write(*,1001)pwtot_e*1.d-7
 1001 format(' total absorbed power_e(watt) =',1pe14.6)
      write(*,1002)pwtot_i*1.d-7
 1002 format(' total absorbed power_i(watt) =',1pe14.6)
      !YuP if (iabsorp.eq.3) then
         write(*,10021)(kk-1,pwtot_s(kk)*1.d-7,kk=2,nbulk)
10021    format(' total absrbd power_s(power) for ion species',
     1       i2,' =',1pe14.6)
      !endif
      write(*,1003)pwtot_cl*1.d-7
 1003 format(' total absorbed power_cl(power)=',1pe14.6)
 1004 format(' total tor curr drive per LL memo[Cohen/(1+eps)] (A)=',
     1     1pe14.6)
 1005 format(' number of rays with bad initial conditions ='1i6)
 1006 format(' total launched power with good initial conditions',/,
     1 '                        (watt)=',1pe14.6)
      write(*,1007)parallel_cur_total
      write(*,1008)toroidal_cur_total
      write(*,1009)poloidal_cur_total
 1007 format('total parallel current parallel_cur_total (A)=',
     &    1pe14.6)
 1008 format('total toroidal current toroidal_cur_total (A)=',
     &    1pe14.6)
 1009 format('total poloidal current poloidal_cur_total (A)=',
     &    1pe14.6)

      return
      end ! dnonetwo_xyz

c======================================================================
c======================================================================

      subroutine efKarney_xyz(z,r,r0x,z_eff,temp,den,jwave,cnpar,
     + effKarn) !->out
c------------------------------------------------------------------
c  RF current drive efficiency(asymptotic formula, nonrelativistic case)  
c  D.A. Ehst and C.F.F.Karney Nucl.Fus. Vol.31, No.10 (1991), p 1933-1938
c  A subroutine to calculate the local current density driven by RF waves
c--Uses CD efficiency empirical formula based on numerical Fokker-
c  Planck bounce-averaged calculations
c  This formula is for
c  1)Landau damping of lower hybrid (LH) slow waves resonant  at parallel
c    velocities	above the electron thermal velocity
c  2)Slow frequincy fast (compressional Alfven) wave (AW) may resonant with
c    low phase velocity electrons via combined Landay damping and transit
c    time magnetic damping  
c------------------------------------------------------------------
c     input parameters: z,r -coordinates of the ray point(cm!!!!!)ATTENTION
c                       z_eff-effective charge
c                       temp-temperature (keV)
c                       den-density(10**13 cm**-3)
c                       jwave (=islofa)=-1 ALfen wave , 0- Landau damp.
c                       cnpar -paralell to magnetic field refractive
c                              index
c                       r0x character length
c     output parameter:  J/P=efficient  in (A/m**2)/(Wt/m**3)
c------------------------------------------------------------------
c     output parameter: efficien  in (A/cm**2)/(erg/(c*cm**3))
c------------------------------------------------------------------
c     It uses:
c     double precision function: psif_xyz(x,y,z) and
c     subroutine: zr_psith(psi,theta,z,r)
c     double precision functions from zr_psith.f: rmax_psi(psi),
c     rmin_psi(psi), bmax_psi(psi)
c------------------------------------------------------------------
      implicit integer (i-n), real*8 (a-h,o-z)
      islofa=jwave
      zeff=z_eff

c      write(*,*)'in efKarney jwave,z_eff',jwave,z_eff

      if(islofa.eq.-1) then	!Alfven damping
        akcd = 11.91d0/(0.678d0+zeff)
        c0cd = 4.13d0/zeff**0.707d0
        amcd = 2.48d0
        ccd = 0.0987d0
        acd = 12.3d0
      else			!Landau damping
        akcd = 3.d0/zeff
        c0cd = 3.83d0/zeff**0.707d0
        amcd = 1.38d0
        ccd = 0.389d0
        acd = 0.d0
      endif
c      write(*,*)'in efKarney (K,D,m,c,a) kcd,c0cd,amcd,ccd,acd',
c     &                       akcd,c0cd,amcd,ccd,acd
c--------------------------------------------------------------------
c     aspct is inverse aspct ratio
c--------------------------------------------------------------------
c     bmax,bmin  in toray
c     if(dabs(bmax-bmin).lt.1.0d-8) bmax=bmin+1.0d-8
c     aspct=(bmax-bmin)/(bmax+bmin)
c--------------------------------------------------------------------
c     Here :aspct is calculating using rmax_psi and rmin_psi
c     conversion from cm to m
      zd=z*0.01d0/r0x
      rd=r*0.01d0/r0x
      phi=0.d0 !   ????
      xd=rd
      yd=0.d0 ! ok for now ?
c----------------------------------------
      psid=psif_xyz(xd,yd,zd)
c      write(*,*)'efKarney psid',psid
c----- rmaxpsi and rminpsi are the largest and and the least
c     values of the radius (rmaxpsi and rminpsi in (m))
c     on the given flux surface psi=psi(zd,rd)
      rmaxpsi=rmax_psi(psid)
      rminpsi=rmin_psi(psid)
      if (rmaxpsi.lt.rd)rmaxpsi=rd
      if (rminpsi.gt.rd)rminpsi=rd
      if(dabs(rmaxpsi-rminpsi).lt.1.0d-8) rminpsi=rmaxpsi-1.0d-8
      aspct=(rmaxpsi-rminpsi)/(rmaxpsi+rminpsi)
      epsil=aspct
 
c      write(*,*)'efKarney eplis',epsil
c----------------------------------------
      bmod=bxyz(xd,yd,zd) !-> get bmod
      bmax=bmax_psi(psid)
      if (bmod.gt.bmax) bmod=bmax
c----------------------------------------
c     cvac - light velocity in (cm/sec)
      cvac=2.997930d+10
c     ve=(sqrt(Te/me)),    ve in (cm/sec),Te(keV)     !article p.1935
      ve=1.32d9*dsqrt(temp)	                      !temp keV
c     normalized resonance parallel velocity:
c     wte=u_res=v_par/v_e=cvac/(ve*cnpar)
      wte=cvac/(ve*cnpar)			      !article p.1934

c      write(*,*)'cvac,temp,ve,cnpar,wte',cvac,temp,ve,cnpar,wte
      write(*,*)'Ehst Karney wte',wte
      u1=wte/dsqrt(2.d0) 
cSAP080905
      wte=dabs(wte)
c----------------------------------------
      alt2=1.d0-bmod/bmax			      !article p.1935 
      alt2=dabs(alt2)
      alt1 =dsqrt(alt2)
c     write(*,*)'efKarney (labmda_t) alt1',alt1
c----------------------------------------
      if (alt2.ne.0.d0) then
         ytt = (1.d0-alt2)*wte**2/alt2		      !article (11)
         rprof = 1.0d0-(epsil**0.77d0*dsqrt(12.25d0+wte**2))/
     1  (3.5d0*epsil**0.77d0+wte)					      !article (7)

c         write(*,*)'efKarney epsil,wte,epsil**0.77d0',
c     &                       epsil,wte,epsil**0.77d0
c         write(*,*)'efKarney epsil**0.77d0*dsqrt(12.25d0+wte**2))',
c     &                       epsil**0.77d0*dsqrt(12.25d0+wte**2)
c         write(*,*)'efKarney (3.5d0*epsil**0.77d0+wte)',
c     &                       (3.5d0*epsil**0.77d0+wte)
c         write(*,*)'efKarney in R second term',
c     &  (epsil**0.77d0*dsqrt(12.25d0+wte**2))/
c     &  (3.5d0*epsil**0.77d0+wte)

         arg = (ccd*ytt)**amcd
         cprof = 1.d0-dexp(-arg)	      	       !article (9)
         amprof = 1.d0+acd*(alt1/wte)**3	       !article (10)
      else
        rprof=1.d0
        cprof=1.d0
        amprof=1.d0
      endif
c      write(*,*)'efKarney (R) rprof',rprof
      eta0 = akcd/(dabs(wte))+c0cd+4.d0*wte**2/(5.d0+zeff) !article (14)
      eta=cprof*amprof*eta0*rprof 		      !article (13)
c      write(*,*)'4.d0*wte**2/(5.d0+zeff)',4.d0*wte**2/(5.d0+zeff) 
c      write(*,*)'eta,cprof,amprof,eta0,rprof ',
c     &           eta,cprof,amprof,eta0,rprof 
       write(*,*)'eta0,eta',eta0,eta
c-----------------------------------------------------------
c     efficiency  in (A/m**2)/(joule/(c*m**3))
c     temperature temp in kev
c     density     dens in 10**19 /m**3
c-----------------------------------------------------------
c     cln -Coulomb Ln
      arg1 = 1.d3/temp*dsqrt(10.d0*den)	!temp KeV, den 10**13/cm**3
      cln=24.d0-dlog(arg1)
      effKarn=eta*3.84d0*temp/(cln*den)	!(a/m**2)/joule/(c*m**3)
c      write(*,*)'temp,den,cln,eta,effKarn',temp,den,cln,eta,effKarn
      write(*,*)'effKarn   (a/m**2)/joule/(c*m**3)', effKarn
c-----------------------------------------------------------
c     efficiency  in (A/cm**2)/(erg/(sec*cm**3))
c     temperature temp in kev
c     density     den in 10**19/m**3
c-----------------------------------------------------------
ctest beg
c     LH wave
      efficien=2.d0/(5.d0+z_eff)*
     1         (u1*u1+(7.0d0/4.0d0+9.0d0/4.0d0/(3.d0+z_eff)))
      efficien=efficien*temp/den*(17.d0/cln)*9.d0*1.d-6
      write(*,*)'u1',u1  
c     write(*,*)'asimptotic efficiency A/cm/(egr(c*cm**3))',efficien
      efficien=8.d0/(5.d0+z_eff)*u1*u1 !test
      write(*,*)'test efficien ',efficien

c      efficien=efficien*temp/den*(17.d0/cln)*9.d0*1.d-6 
c     write(*,*)'efficiency A/cm/(egr(c*cm**3))',efficien
      efficien=efficien*3.84d0*temp/den/cln
      write(*,*)'asymptotic efficiency(A/m**2)/(joule/cm**3))',efficien
ctest end    
c-------------------------------------------------------------
c     determination the of the current drive sign
      if (u1.ge.0.d0) then
         s=1.d0
      else
         s=-1.d0
      endif
      efficien=efficien*s
c-------------------------------------------------------------
cSm050923
c      effKarn=effKarn*s*1.d-5  !  A/cm**2/(erg/(sec*cm**3)
      effKarn=-effKarn*s*1.d-5  !  A/cm**2/(erg/(sec*cm**3)
c      write(*,*)'effcient A/cm**2/(egr/(sec*cm**3))',efficien
c      write(*,*)'effKarn  A/cm**2/(erg/(sec*cm**3))',effKarn
c     stop
      return
      END


c======================================================================
c======================================================================



      subroutine average_variables_xyz(psi_in,
     + b_av,bs_av,r_av,dr_av,drs_av) !->out
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
      real*8 bxyz
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
     &x,y,z,r,pl,b_pol,
     &psi_l, ppp,gradpsi,cs,sn

        
      integer i1p(2),itabl(3)
      integer i_first,i,j

      data i_first /0/
      save i_first,
     &b_av_ar,d2_b_psi,
     &bs_av_ar,d2_bs_psi,
     &r_av_ar,d2_r_psi,
     &dr_av_ar,d2_dr_psi,
     &drs_av_ar,d2_drs_psi

      if (i_first.eq.0) then
c-------------------------------------------------------------------
c        calculate the array b_av(npsi) of the averaged magnetic field
c        along the magnetic surfaces
c        and calulate cubic spline coefficients
c--------------------------------------------------------------------
         pollength(1)=0.d0 ! poloidal length of the magnetic surface
c         write(*,*)'b_average zpsi(1,1),rpsi(1,1)',zpsi(1,1),rpsi(1,1)
         z=zpsi(1,1)
         r=rpsi(1,1) ! can be 0 
c         write(*,*)'rma,zma',rma,zma
c         write(*,*)'rpsi(1,1),zpsi(1,1)',rpsi(1,1),zpsi(1,1)
         x=r
         y=0.d0 ! ok for now ? needs work to do average over all flux surf.
         b_av_ar(1)= bxyz(x,y,z) !->get b
         bs_av_ar(1)=b_av_ar(1)**2
         r_av_ar(1)= r
         if(r.gt.1.d-33)then
          dr_av_ar(1)= 1.d0/r
          drs_av_ar(1)= 1.d0/r**2
         else
          dr_av_ar(1)= 0.d0
          drs_av_ar(1)= 0.d0
         endif
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
            
	    do 20 i=1,nteta
              z=0.5d0*(zpsi(j,i)+zpsi(j,i+1))
              r=0.5d0*(rpsi(j,i)+rpsi(j,i+1))
              x=r
              y=0.d0 ! ok for now ? needs work to do average over all flux surf.
              bmod= bxyz(x,y,z) !-> get b and grad(psi)
              ! Bpoloid = (b.[grad(psi) X e_phi]) / |gradpsi| 
              ! e_phi= {-y/r      ,   x/r     ,     0    }
              ! Find poloidal magnetic field
              ppp=(dpdxd*dpdxd+dpdyd*dpdyd+dpdzd*dpdzd)
              gradpsi=dsqrt(ppp)
              cs=x/r
              sn=y/r
              b_pol=(bz*(dpdxd*cs+dpdyd*sn)-(bx*cs+by*sn)*dpdzd)/gradpsi
               pl=dsqrt((rpsi(j,i+1)-rpsi(j,i))**2+
     1                 (zpsi(j,i+1)-zpsi(j,i))**2)
               pollength(j)=pollength(j)+pl

               dbpol_av_ar(j)= dbpol_av_ar(j)+pl/b_pol
               b_av_ar(j)=b_av_ar(j)+bmod*pl/b_pol
               bs_av_ar(j)=bs_av_ar(j)+bmod**2*pl/b_pol
               r_av_ar(j)= r_av_ar(j)+r*pl/b_pol
               dr_av_ar(j)=dr_av_ar(j)+pl/(r*b_pol)
               drs_av_ar(j)=drs_av_ar(j)+pl/(r**2*b_pol)
 20         continue
            b_av_ar(j)=b_av_ar(j)/dbpol_av_ar(j) !normalize
            bs_av_ar(j)=bs_av_ar(j)/dbpol_av_ar(j)
            r_av_ar(j)=r_av_ar(j)/dbpol_av_ar(j)
            dr_av_ar(j)=dr_av_ar(j)/dbpol_av_ar(j)
            drs_av_ar(j)=drs_av_ar(j)/dbpol_av_ar(j)
 10      continue

c----------------------------------------------------------
c        calculate cubic spline coefficients for the function
c        that gives the averaged magnetic field:  b_av_psi(psi)              
c---------------------------------------------------------       
         i1p(1)=4
         i1p(2)=4 
         call coeff1(npsi,arpsi,b_av_ar,d2_b_psi,i1p,1,work)
         call coeff1(npsi,arpsi,bs_av_ar,d2_bs_psi,i1p,1,work)
         call coeff1(npsi,arpsi,r_av_ar,d2_r_psi,i1p,1,work)
         call coeff1(npsi,arpsi,dr_av_ar,d2_dr_psi,i1p,1,work)
         call coeff1(npsi,arpsi,drs_av_ar,d2_drs_psi,i1p,1,work)
c-----------------------------------------------------------
        i_first=1

      endif !spline coefficients are calculated 

c-----------------------------------------------------------
c     calculation b_averaged 
c------------------------------------------------------------
      itabl(1)=1
      itabl(2)=0
      itabl(3)=0
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
      return
      end



c======================================================================
c======================================================================



c        ********************** hamilt_xyz ******************
c        *                      ------                      *
c        * this function calculates the hamiltonian of the  *
c        * the system of geometrical optics equations       *
c        * It will calculate the dielectric tensor
c        * reps(3,3) and will put this tensor to eps.i file
c        ****************************************************
c
c-------------------------------------------------------------------
c								    
c        input parameters					    
c      								    
c      x,y,z - cartesian coordinates of the  point    
c              where  the  hamiltonian  is calculated.      		         	    
c								    !
c      cnt2= cnper**2+cnpar**2 == cn_x*cn_x + cn_y*cn_y + cn_z*cn_z
c      (Square of total refractive index at this point)
c      The angle gam between refractive index n and magnetic field
c      is taken from common 'one'					    !
c-------------------------------------------------------------------
      double precision
     1function hamilt_xyz(x,y,z,cnt2)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'eps.i'
      include 'ions.i'
      double complex dhot,dhot_sum
      double complex ceps(3,3),hamiltc
      external wcw,wpw_2,tempe_xyz
      external dhot_sum
      !--------------
      double precision x_ar(nbulka),y_ar(nbulka),t_av_ar(nbulka),
     .tpop_ar(nbulka),vflow_ar(nbulka)
      double complex K(3,3),dK_dx_ar(3,3,nbulka),dK_dy_ar(3,3,nbulka),
     &dK_dt_ar(3,3,nbulka),dK_dnper(3,3),dK_dnpar(3,3)
      double complex compl_nper !for Eric tensor
      double complex eps_weiss(3,3)
c-----external
      double complex det       
      !--------------

      hamilt_xyz=0.d0 ! to initialize for all id except those below.
      
      bmod=bxyz(x,y,z) !-> get bmod for wcw
      
      cn=dsqrt(cnt2) ! |N|
      cn2=cnt2       ! N^2 (input)
      cn4=cn2*cn2    ! N^4
           
      ds=dsin(gam)
      dc=dcos(gam)
      cnpar=cn*dc ! |Npar|
      cnper=cn*ds ! Nperp
      ds2=ds*ds
      dc2=dc*dc
      if (id.eq.1 .or. id.eq.2 .or. id.eq.3) then
        call tensrcld_xyz(x,y,z)
      end if
c---------------------------------------------------------
c     Appleton-Hartree dispersion relation
c---------------------------------------------------------
      if (id.eq.3) then
         ds4=ds2*ds2
         xi=wpw_2(x,y,z,1)
         yi=wcw(x,y,z,1)
         py2=yi*yi
         py4=py2*py2
         px=1.d0-xi
         px2=px*px
         sqrdet=dsqrt(py4*ds4+4.*py2*px2*dc2)
         pz=2.d0*px-py2*ds2+ioxm*sqrdet
         hamilt_xyz=cn2-(1.d0-2.d0*xi*px/pz)
         goto 10
      end if
c-----------------------------------------------------------

c--------------------------------------------------------
c       id=1 or id=2 cold plasma dispersion relation
c--------------------------------------------------------
      if ((id.eq.1).or.(id.eq.2)) then
        call abc_xyz(x,y,z,ds2,dc2,ad,bd,cd)
        d4=ad
       	d2=bd
       	d0=cd
        if (id.eq.1) then
           hamilt_xyz= d4*cn4+d2*cn2+d0
        end if

        if (id.eq.2) then
           hamilt_xyz= cn2+(d2-ioxm*dsqrt(d2*d2-4.d0*d4*d0))/
     *                        (2.d00*d4)
        end if
        go to 10
      end if
c--------------------------------------------------------

c-----------------------------------------------------------
c     Hot non-relativistic plasma
c-----------------------------------------------------------
      if (id.eq.6) then
         do i=1,nbulk
           x_ar(i)=wpw_2(x,y,z,i)
           y_ar(i)=wcw(x,y,z,i)
           if(i.eq.1) y_ar(1)=-y_ar(1)
           te= tempe_xyz(x,y,z,i) ! kev ! as a function of rho, for now
           t_av_ar(i)=te*1000.d0      ! ev 
           tpop_ar(i)= tpoprho(rho,i) ! rho is in one.i
           vflow_ar(i)= vflowrho(rho,i)
         enddo
         !write(*,*)'before hamiltc: x,y,z, Y=',x,y,z, y_ar(1)
         hamiltc=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     .                     vflow_ar,cnpar,cnper,1,reps)

         hamilt_xyz=dreal(hamiltc)
         go to 10
      end if
c-----------------------------------------------------------

c-------------------------------------------------------------
c     id=14 dispersion function
c-----------------------------------------------------------
c     D=real part(D relativistic dispersion function)
c     from Eric Nelson_Melby dielectric tensor ,if npar.le. 0.38D0
c     from Ram Abhay dielectric tensor ,if npar.ge. 0.38D0
C 1Sep2005 -- Now Disp_combined is entirely Abhay's tensor, which 
C works just as well or faster than Nelson-Melby's, when the resolution
C is lower than the maximum like it used to be.
C
C ENM 15Mar2006 -- After finding that the combined version which jumps from
C     one to another dispersion relation depending on n_parallel didnot work
C     very well, and finding that the Ram version works well, even for fairly
C     small n_parallel, as long as you adjust the resolution parameters (see
C     genray.dat template), now id=14 is just the Ram tensor (id=11 is just the
C     Nelson-Melby tensor)
C
      if (id.eq.14) then
         x_ar(1)=wpw_2(x,y,z,1)
         y_ar(1)=wcw(x,y,z,1)    !question for electron -y? 
         te= tempe_xyz(x,y,z,1) !kev
         compl_nper=dcmplx(cnper,0.d0)

         call Disp_Ram(te,cnpar,x_ar(1),y_ar(1),
     +      compl_nper,K,hamiltc) !K is in z-y plane in Stix coordinates

         if (iherm.eq.1) then
c-----------calcualte hermitian part reps of the complex tensor K
            call herm(K,reps)
            hamiltc= det(reps,cnpar,compl_nper)   
         endif
         hamilt_xyz=dreal(hamiltc)
         go to 10
      end if ! end if id=14

 10   continue
c      write(*,*)'hamilt_xyz: id,ioxm=',id,ioxm,hamilt_xyz
      return
      end
c======================================================================
c======================================================================



c        ********************tensrcld_xyz**********************
c        * 						      *
c        * this subroutine  calculates the components of      *
c        * the complex dielectric tensor reps(3,3)            *
c        * for cold plasma                                    *
c        *****************************************************
c-------------------------------------------------------------
c       input parameters				    
c     	x,y,z coordinates  of ray point		    
c       .....						     
c       rho-small radius is calculated in subroutine: b   
c            rho is in common one    			     
c       output parameter				     
c       reps(3,3)-complex dielectric tensor in common /eps/  
c------------------------------------------------------------
      subroutine tensrcld_xyz(x,y,z)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'eps.i'
      !--------------
       xe= wpw_2(x,y,z,1) ! electrons
       bmod=bxyz(x,y,z) !-> get b (needed for wcw)
       ye= wcw(x,y,z,1)
       ye1= (1.d0-ye*ye)
       if(ye1.eq.0.d0) then
         write(*,*)'tensrcld_xyz: WARNING w=wce resonance. ye=',ye
         pause
       endif
       epsper= 1.d0-xe/(1.d0-ye*ye)
       g= xe*ye/(1.d0-ye*ye)
       epspar= 1.d0-xe
	 do i=2,nbulk ! ions
          xi= wpw_2(x,y,z,i)
          yi= wcw(x,y,z,i)
          epsper= epsper-xi/(1.d0-yi*yi)
          epspar= epspar-xi
          g= g-xi*yi/(1.d0-yi*yi)
	 enddo
c-------------------------------------------------
	 reps(1,1)=dcmplx(epsper,0.d0)
	 reps(1,2)=dcmplx(0.d0,g)
	 reps(2,1)=-reps(1,2)
	 reps(2,2)= reps(1,1)
	 reps(3,3)=dcmplx(epspar,0.d0)
	 reps(1,3)=dcmplx(0.d0,0.d0)
	 reps(3,1)=dcmplx(0.d0,0.d0)
	 reps(2,3)=dcmplx(0.d0,0.d0)
	 reps(3,2)=dcmplx(0.d0,0.d0)
c-------------------------------------------------
	 return
	 end
c======================================================================
c======================================================================




c        ********************** abc_xyz   *******************
c        * this subroutine calculates the coefficients      *
c        * a,b,c, for the cold dispersion relation          *
c        * d=a*n**4+b*n**2+c                                *
c        ****************************************************
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
c      	input parameters					    	
c       point coordinates: x,y,z    				    
c       ds2=sin(gam)**2,  dc2=cos(gam)**2                 
c       gam is the angle between magnetic field and refractiv index !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
      subroutine abc_xyz(x,y,z,ds2,dc2, ad,bd,cd)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i' ! contains gam
c  ---------------------------------------------------------------
c  cold plasma dispersion relation with electrons and ions
c  wpw_2(i=1),wcw(i=1) for electron component
c  wpw_2(i may be=2,nbulk),wcw(i may be=2,...,nbulk) ions components.
c  ib is a number of component for which delib=1-yib may be equal
c  zero inside the plasma,ib may be=from 1 to nbulk
c  The dispersion relation is multiplied by delib
c  --------------------------------------------------------------
      pc=1.d0+dc2
      call s_xyz(x,y,z,s1,s2,s3,s4,s6,s7)
      xib=wpw_2(x,y,z,ib)
      bmod=bxyz(x,y,z) !-> get b (needed for wcw)
      yib=wcw(x,y,z,ib)
      delib=1.d0-yib
      if(delib.eq.0.d0)then
        write(*,*)'abc_xyz: WARNING: w=wc resonance i=',ib
        pause
      endif
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
         xe=wpw_2(x,y,z,1)
         ye=wcw(x,y,z,1)
         ppe=1.d0/(1.d0-ye)
         peym=xe/(1.d0-ye)
         peyp=xe/(1.d0+ye)
         peym2=peyp*ppe
         !---------------
         ppb=1.d0/(1.d0-yib)
         pbym=xib/(1.d0-yib)
         pbyp=xib/(1.d0+yib)
         pbym2=pbyp*ppb
         !---------------
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
c------------------------------------------------------------------
      return
      end
c======================================================================
c======================================================================

c        **********************   s_xyz   ********************
c        *                        -                          *
c        * this subroutine calculates  six  different  sums  *
c        * that are used by rside for cold plasma dispersion *
c        *****************************************************
c
c-------------------------------------------------------------------
c								    
c        input parameters					    
c      								    
c     x,y,z - cartesian coordinates of the point 	    
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
c								    !
c        output parameters                                          !
c								    !
c     s1 = 1 - sum_(i.ne.ib) { wpw_2_i/[1-(wcw_i)**2] }, at ib.ne.1 
c								    !
c     s2 = 1 - sum_(i.ne.ib) { wpw_2_i/[1-wcw_i] },      at ib.ne.1 
c  								    !
c     s3 = 1 - sum_(i) { wpw_2_i/[1+wcw_i] },    
c								    !
c     s4 = 1 - sum_(i) { wpw_2_i } - wpw_2_1 , 
c								    !
c     s6 = 1 - sum_(i) {wpw_2_i/[1-wcw_i] ,	             at ib.eq.1
c								    !
c     s7 = 1 - sum_(i) {wpw_2_i/[1-(wcw_i)**2] ,	     at ib.eq.1
c								    !
c     i=2,nbulk; ib=1,nbulk; nbulk - number of bulk particles       !
c-------------------------------------------------------------------
      subroutine s_xyz(x,y,z, s1,s2,s3,s4,s6,s7)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      ds1=0.d0
      ds2=0.d0
      ds3=0.d0
      ds4=0.d0
      ds6=0.d0
      ds7=0.d0
      bmod=bxyz(x,y,z) !-> get b (needed for wcw)
      !------------------------
      if (nbulk.gt.1) then
        do 10 i=2,nbulk
          xi= wpw_2(x,y,z,i)
          yi= wcw(x,y,z,i)
          ds3=ds3+xi/(1.d0+yi)
          ds4=ds4+xi
          !------------------------
          if (ib.eq.1) then
            ds6=ds6+xi/(1.d0-yi)
            ds7=ds7+xi/(1.d0-yi*yi)             
          end if
          !------------------------
          if ((i.eq.ib).or.(ib.eq.1)) goto 10
          ds1=ds1+xi/(1.d0-yi*yi)
          ds2=ds2+xi/(1.d0-yi)  
 10     continue
      end if ! finish if nbulk.gt.1
      !------------------------
      pps= wpw_2(x,y,z,1) ! electrons
      !write(*,*)'s_xyz: x,y,z,nbulk,wpw_2',x,y,z,nbulk,pps 
      s1=1.d0-ds1
      s2=1.d0-ds2
      s3=1.d0-ds3      
      s4=1.d0-ds4-pps
      s6=1.d0-ds6
      s7=1.d0-ds7
      !write(*,*)'s_xyz: s1-s7',s1,s2,s3,s4,s6,s7 
      return
      end ! s_xyz()/END
      
c======================================================================
c======================================================================

c        ********************** rside_xyz *********************
c        *                      -----                         *
c        * this subroutine calculates the right hand sides    *
c        * of the system of geometrical optics equations      *
c        * for cold, Appelton-Hartree thermal plasma          *
c        * for 6 ray equations                                *
c        * The hamiltonian derivatives  are calculated	      *
c        * analytically(idif=1) or numerically (idif=2)	      *
c        ******************************************************
c
c-------------------------------------------------------------------
c								    
c        input parameters					    
c      								    
c      t - parameter of trajectory 
c          Not an actual physical time, but 
c          normalized as t= t_physical[sec] *c/omega                               
c                                  		         	    
c      u - solution of geometrical optics equations at the point t  
c      u(1) = x                                                     !
c      u(2) = y                                                     !
c      u(3) = z                                                     !
c      u(4) = n_x                                                   !
c      u(5) = n_y                                                   !
c      u(6) = n_z                                                   !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
c      							   	 
c        output parameters					  
c                                                                   !
c     deru(i) are right  hand  sides  of  geometrical  optics       !
c     equations							    
c     deru(1)=-{d(hamiltonian)/d(N_x)}/{d(hamiltonian)/d(omega)}    !
c     deru(2)=-{d(hamiltonian)/d(N_y)}/{d(hamiltonian)/d(omega)}    !
c     deru(3)=-{d(hamiltonian)/d(N_z)}/{d(hamiltonian)/d(omega)}    !
c     deru(4)={d(hamiltonian)/d(x)}/{d(hamiltonian)/d(omega)}       !
c     deru(5)={d(hamiltonian)/d(y)}/{d(hamiltonian)/d(omega)}       !
c     deru(6)={d(hamiltonian)/d(z)}/{d(hamiltonian)/d(omega)}       !
c     
c-------------------------------------------------------------------
c     this program uses following functions and subroutines	 
c     bxyz, gamma1_xyz,
c     dwpw_2,dwcw,s_xyz,	
c     hamilt_xyz    						
c-------------------------------------------------------------------
      subroutine rside_xyz(u,deru)
c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none
      include 'param.i'
      include 'one.i'
c this line only for call control
      include 'eps.i'
      include 'ions.i'
      
c-----input
      double precision t,u(6), x,y,z, cn_x,cn_y,cn_z, cnx,cny,cnz, cn2
c-----output
      double precision deru(6)
c-----externals
      double precision bxyz, gamma1_xyz, hamilt_xyz,
     &cn, wpw_2, wcw
c-----locals
      double precision vp(nbulka),wp(nbulka)
      double precision wf, r,
     &hw,hfrqnc,
     &hp,hm,h00,
     &dddz,
     &dddx,
     &dddy,
     &frqncpl,df,frqncmn,dddw,p,
     &ds,dc,ds2,dc2, 
     &dwpw2dz,dwpw2dx,dwpw2dy,dwcwdz,dwcwdx,dwcwdy,
     &xi,yi,py2,py3,py4,px,ds4,det,zer,sqrdet,p1,
     &ddetdx,ddetdy,pz,pz2,dpzdx,dpzdy,dddxi,dddyi,dddc2,
     &dxdw,dydw,
     &pc,s1,s2,s3,s4,s6,s7,xib,yob,delib,
     &peym,peyp,a0e,a1e,b0e,b1e,c0e,c1e,dele,ad,bd,cd,
     &px2,yib,peym2,
     &da1edc,da0edc,db1edc,db0edc,
     &da1ed,da0ed,db1ed,db0ed,dc1edc,dc0edc,dadc2,dbdc2,dcdc2,
     &ppe,ppb,pbym,pbyp,pbym2,a0b,a1b,b0b,b1b,c0b,c1b,
     &da1bdc,da0bdc,db1bdc,db0bdc,dc1bdc,dc0bdc,
     &dadz,dadx,dady,dbdz,dbdx,dbdy,dcdz,dcdx,dcdy,
     &dadw,dbdw,dcdw,pyim,pyip,pyim2,
     &ds1dxi,ds2dxi,ds1dyi,ds2dyi,ds3dxi,ds4dxi,ds3dyi,ds4dyi,
     &ds6dxi,ds7dxi,ds6dyi,ds7dyi,
     &da1exi,da0exi,db1exi,db0exi,dc1exi,dc0exi,
     &da1eyi,da0eyi,db1eyi,db0eyi,dc1eyi,dc0eyi,
     &dadxi,dbdxi,dcdxi,dadyi,dbdyi,dcdyi,
     &da1bxi,da0bxi,db1bxi,db0bxi,dc1bxi,dc0bxi,
     &da1byi,da0byi,db1byi,db0byi,dc1byi,dc0byi,
     &ds1dxb,ds2dxb,ds1dyb,ds2dyb,ds3dxb,ds4dxb,ds3dyb,ds4dyb,
     &ds2dxe,ds2dye,
     &ds6dxb,ds7dxb,ds6dyb,ds7dyb,
     &da1bxb,da0bxb,db1bxb,db0bxb,dc1bxb,dc0bxb,
     &da1byb,da0byb,db1byb,db0byb,dc1byb,dc0byb,
     &xe,ye,
     &ds1dxe,ds3dxe,ds4dxe,ds6dxe,ds7dxe,
     &ds1dye,ds3dye,ds4dye,ds6dye,ds7dye,
     &da1bxe,da0bxe,db1bxe,db0bxe,dc1bxe,dc0bxe,
     &da1bye,da0bye,db1bye,db0bye,dc1bye,dc0bye,
     &da1exe,da0exe,db1exe,db0exe,dc1exe,dc0exe,
     &da1eye,da0eye,db1eye,db0eye,dc1eye,dc0eye,
     &dcn,dcn2,dcn4,dadnz,dbdnz,dcdnz,
     &dadnx,dbdnx,dcdnx,dadny,dbdny,dcdny,
     &p4,p5,p6,p7,
     &sum,dddn,
     +dddcnx,dddcny,dddcnz,
     +xp,xm,yp,ym,zp,zm,cn2plus,cn2minus

      ! dwpw2dx,dwpw2dy,dwpw2dz ! derivs of (wp/w)^2
      ! dwcwdx,dwcwdy,dwcwdz ! derivs of wc/w
      real*8 Vgr_c_2,Vgr_c,scale_vgr 
      
      integer i

      x= u(1) 
      y= u(2) 
      z= u(3) 
      cn_x= u(4) ! Nx
      cn_y= u(5) ! Ny
      cn_z= u(6) ! Nz
      cn2= cn_x*cn_x + cn_y*cn_y + cn_z*cn_z
      r=sqrt(x*x+y*y)
c-------------------------------------------------------
      wf=frqncy
      do 11 i=1,nbulk
         vp(i)=v(i) !  just to save these values
         wp(i)=w(i)
 11   continue
c-------------------for analytical derivatives-----------
c	idif=1
c-------------------for numerical derivatives------------
c	idif=2
c------------------------------------------------------------------      
      if (idif.eq.1) goto 1955 ! skip idif.eq.2 below
c----------------------------------------------------------------
c idif=2 numerical hamiltonian derivative calculation
c----------------------------------------------------------------
      hw=der_f !! for frequency
      pi=3.1415926d0
      hfrqnc= der_f*frqncy ! hfrqnc*frqncy ! 
c----------------------------------------------------------------
      bmod=bxyz(x,y,z) ! to find bx,by,bz,bmod -> /one.i/
      cn2= cn_x*cn_x + cn_y*cn_y + cn_z*cn_z
      gam= gamma1_xyz(x,y,z, cn_x,cn_y,cn_z) ! uses bx,by,bz,bmod /one.i/
      h00=  hamilt_xyz(x,y,z, cn2) ! uses gam; calls b(), wpw_2(), wcw()

      cnx= cn_x+der_n
      cn2= cnx*cnx + cn_y*cn_y + cn_z*cn_z
      gam= gamma1_xyz(x,y,z, cnx,cn_y,cn_z) ! uses bx,by,bz,bmod /one.i/
      hp=  hamilt_xyz(x,y,z, cn2) ! uses gam; calls b(), wpw_2(), wcw()
      cnx= cn_x-der_n
      cn2= cnx*cnx + cn_y*cn_y + cn_z*cn_z
      gam= gamma1_xyz(x,y,z, cnx,cn_y,cn_z) 
      hm=  hamilt_xyz(x,y,z, cn2) ! uses gam
      dddcnx=(hp-hm)/(2.d0*der_n)
      deru(1)=dddcnx
      !--- Another version: one-side derivative:
      !dddcnx=(hp-h00)/(der_n)
      !dddcnx=(h00-hm)/(der_n)
      !deru(1)=dddcnx
      !----------------------------------------

      cny= cn_y+der_n
      cn2= cn_x*cn_x + cny*cny + cn_z*cn_z
      gam= gamma1_xyz(x,y,z, cn_x,cny,cn_z) 
      hp=  hamilt_xyz(x,y,z, cn2) ! uses gam
      cny= cn_y-der_n
      cn2= cn_x*cn_x + cny*cny + cn_z*cn_z
      gam= gamma1_xyz(x,y,z, cn_x,cny,cn_z) 
      hm=  hamilt_xyz(x,y,z, cn2) ! uses gam
      dddcny=(hp-hm)/(2.d0*der_n)
      deru(2)=dddcny

      cnz= cn_z+der_n
      cn2= cn_x*cn_x + cn_y*cn_y + cnz*cnz
      gam= gamma1_xyz(x,y,z, cn_x,cn_y,cnz) 
      hp=  hamilt_xyz(x,y,z, cn2) ! uses gam
      cnz= cn_z-der_n
      cn2= cn_x*cn_x + cn_y*cn_y + cnz*cnz
      gam= gamma1_xyz(x,y,z, cn_x,cn_y,cnz) 
      hm=  hamilt_xyz(x,y,z, cn2) ! uses gam
      dddcnz=(hp-hm)/(2.d0*der_n)
      deru(3)=dddcnz

      cn2= cn_x*cn_x + cn_y*cn_y + cn_z*cn_z

      xp=  x+der_r
      bmod=bxyz(xp,y,z) !-> get b and derivs of b
      gam= gamma1_xyz(xp,y,z, cn_x,cn_y,cn_z) ! uses bx,by,bz,bmod /one.i/
      hp=  hamilt_xyz(xp,y,z, cn2) ! uses gam; calls b(), wpw_2(), wcw()
      xm=  x-der_r
      bmod=bxyz(xm,y,z) !-> get b and derivs of b
      gam= gamma1_xyz(xm,y,z, cn_x,cn_y,cn_z)
      hm=  hamilt_xyz(xm,y,z, cn2)
      dddx=(hp-hm)/(2.d0*der_r)
      deru(4)=-dddx
      !--- Another version: one-side derivative:
      !dddx=(hp-h00)/(der_r) ! (H(x+dr)-H00)/dr
      !dddx=(h00-hm)/(der_r) ! (H00-H(x-dr))/dr
      !deru(4)=-dddx
      !----------------------------------------
           
      yp=  y+der_r
      bmod=bxyz(x,yp,z) !-> get b and derivs of b
      gam= gamma1_xyz(x,yp,z, cn_x,cn_y,cn_z)
      hp=  hamilt_xyz(x,yp,z, cn2)
      ym=  y-der_r
      bmod=bxyz(x,ym,z) !-> get b and derivs of b
      gam= gamma1_xyz(x,ym,z, cn_x,cn_y,cn_z)
      hm=  hamilt_xyz(x,ym,z, cn2)
      dddy=(hp-hm)/(2.d0*der_r)
      deru(5)=-dddy

      zp=  z+der_r
      bmod=bxyz(x,y,zp) !-> get b and derivs of b
      gam= gamma1_xyz(x,y,zp, cn_x,cn_y,cn_z)
      hp=  hamilt_xyz(x,y,zp, cn2)
      zm=  z-der_r
      bmod=bxyz(x,y,zm) !-> get b and derivs of b
      gam= gamma1_xyz(x,y,zm, cn_x,cn_y,cn_z)
      hm=  hamilt_xyz(x,y,zm, cn2)
      dddz=(hp-hm)/(2.d0*der_r)
      deru(6)=-dddz

c-----------------------------------------------------------
      bmod=bxyz(x,y,z) ! to find bx,by,bz,bmod -> /one.i/
c*************************************************
      frqncpl=frqncy+hfrqnc ! ! = frqncy*(1.d0+der_f)
      df=frqncy/frqncpl
      do 12 i=1,nbulk
         v(i)=vp(i)*df*df ! used through common one.i
         w(i)=wp(i)*df
 12   continue
      cnx=cn_x*df
      cny=cn_y*df
      cnz=cn_z*df
      cn2plus=cn2*df*df
      gam= gamma1_xyz(x,y,z, cnx,cny,cnz)  
      hp=  hamilt_xyz(x,y,z, cn2plus) ! uses gam
c*************************************************
      frqncmn=frqncy-hfrqnc ! = frqncy*(1.d0-der_f)
      df=frqncy/frqncmn
      do 15 i=1,nbulk
         v(i)=vp(i)*df*df ! used through common one.i
         w(i)=wp(i)*df
 15   continue
      cnx=cn_x*df
      cny=cn_y*df
      cnz=cn_z*df
      cn2minus=cn2*df*df
      gam= gamma1_xyz(x,y,z, cnx,cny,cnz)  
      hm=  hamilt_xyz(x,y,z, cn2minus) ! uses gam
c*************************************************
      !write(*,'(6e12.3)')deru
      dddw=(hp-hm)/(2.0d0*hw)
      if(dddw.eq.0.d0)then
        write(*,*) 'rside_xyz:  dddw=0'
        stop
      else
        p=-1.d0/dddw
      endif
      
      do 14 i=1,nbulk
          v(i)=vp(i) ! restore these arrays to original values
          w(i)=wp(i)
 14   continue

      if (i_geom_optic.eq.2) then
c------------the right hand side of ray-tracing equations 
c            (dD/dN)/dl and (dD/dr)/dl
c            It gives the ray equations in the form dr/dl and dN/dl
c            l is the lenth along the ray            
              p=1.d0/dsqrt(deru(1)**2+deru(2)**2+deru(3)**2)
              p=ray_direction*p
      endif
          
      do 16 i=1,6
 16   deru(i)=deru(i)*p
          
c-----------------------------------------------------------
      go to 1953 !-> idif.eq.2 done; skip idif.eq.1 part, return/end
c----------------------------------------------------------------
c  end of numerical calculation of hamiltonian derivatives
c------------------------------------------------------------------
 1955	  continue ! idif.eq.1 starts here
c      write(*,*)'rside_xyz: id,ib,ioxm=',id,ib,ioxm

      x= u(1) 
      y= u(2) 
      z= u(3) 
      cn_x= u(4) ! Nx
      cn_y= u(5) ! Ny
      cn_z= u(6) ! Nz
      cnx=cn_x
      cny=cn_y
      cnz=cn_z
      bmod=bxyz(x,y,z) !-> get b and derivs of b
      gam=gamma1_xyz(x,y,z, cnx,cny,cnz)
      ds=dsin(gam)
      dc=dcos(gam)
      ds2=ds*ds
      dc2=dc*dc

c
c         Appleton-Hartree dispersion relation
c
      if (id.eq.3) then
         if (idif.eq.1) then ! analytical derivs
            stop 'rside_xyz: id=3/idif=1  Not ready yet. Choose idif=2' 
            ! For now, only numerical derivs.
         endif
         call dwpw_2(x,y,z,1, dwpw2dx,dwpw2dy,dwpw2dz) ! derivs of (wp/w)^2
         call dwcw(x,y,z,1, dwcwdx,dwcwdy,dwcwdz) ! derivs of wc/w
         ! need to define dxidr,dxidz,etc.
         xi=wpw_2(x,y,z,1)
         yi=wcw(x,y,z,1)
         py2=yi*yi
         py3=py2*yi
         py4=py2*py2
         px=1.d0-xi
         px2=px*px
         ds4=ds2*ds2
         det=py4*ds4+4.d0*py2*px2*dc2
         zer=0.d0
         if (det.lt.zer) then
            write(*,*)'det in rside less then 0 det=',det
            stop
         end if
         sqrdet=dsqrt(det)
         p1=0.5d0/sqrdet
         ddetdx=-8.d0*py2*px*dc2
         ddetdy=4.d0*py3*ds4+8.d0*yi*px2*dc2
         pz=2.d0*px-py2*ds2+ioxm*sqrdet
         pz2=1.d0/(pz*pz)
         dpzdx=-2.d0+ioxm*p1*ddetdx
         dpzdy=-2.d0*yi*ds2+ioxm*p1*ddetdy
         dddxi=-2.d0*((2.d0*xi-1.d0)*pz+xi*px*dpzdx)*pz2
         dddyi=-2.d0*xi*px*dpzdy*pz2
         dddc2=-2.d0*xi*px*pz2*(py2+ioxm*(-2.d0*py4*(1.d0-dc2)+
     1                    4.d0*py2*px2)*p1)
         dddz=dddxi*dwpw2dz+dddyi*dwcwdz+dddc2*dc2dz
         dddx=dddxi*dwpw2dx+dddyi*dwcwdx+dddc2*dc2dx
         dddy=dddxi*dwpw2dy+dddyi*dwcwdy+dddc2*dc2dy
         dddcnz=2.d0*cnz+dddc2*dc2dnz
         dddcnx=2.d0*cnx+dddc2*dc2dnx
         dddcny=2.d0*cny+dddc2*dc2dny
         dxdw=-2.d0*xi
         dydw=-yi
         dddw=dddxi*dxdw+dddyi*dydw
     1		  -2.d0*(cnz*cnz+cnx*cnx+cny*cny)
         dddw=dddw/wf 
         goto 50
      end if ! id=3
c   end if Appleton - Hartree

c-----------------------------------------------------------------
      if ((id.eq.1).or.(id.eq.2)) then
c  ---------------------------------------------------------------
c  cold plasma dispersion relation with electrons and ions
c  wpw_2(i=1),wcw(i=1) electrons component
c  wpw_2(i may be=2,nbulk),wcw(i may be=2,nbulk) ions components
c  ib number of component for which delib=1-yib may be equal
c  zero inside the plasma, ib may be=from 1 to nbulk
c  dispersion relation is multiplied by delib
c
c         if (idif.eq.1) then ! analytical derivs
c            stop 'rside_xyz: id=1,2/idif=1  Not ready. Choose idif=2' 
c            ! For now, only numerical derivs.
c         endif

	   pc=1.d0+dc2

         call s_xyz(x,y,z,s1,s2,s3,s4,s6,s7)

         xib=wpw_2(x,y,z,ib)
         yib=wcw(x,y,z,ib)
         delib=1.d0-yib

c  ib =1 (cyclotron resonance conditions dele=0 may be
c         realised in plasma, dispersion relation is
c         multiplied by dele )
c         ad,bd,cd calculations in dispersion relation
c       d=a*n**4+b*n**2+c
c------------------------------------------------------------------
c  if 2 begin
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
         da1edc=-s7+s4
         da0edc=peyp
         db1edc=-s7*s4+s3*(s6-peyp)
         db0edc=s4*peyp-xe*(s6-peyp)
         dc1edc=0.d0
         dc0edc=0.d0
         dadc2=dele*da1edc+da0edc
         dbdc2=dele*db1edc+db0edc
         dcdc2=dele*dc1edc+dc0edc
      endif ! ib=1
c  if 2 end
c-----------------------------------------------------------------
c  ib .gt.1 (cyclotron resonance conditions delib=0 may be
c            realised in plasma, dispersion relation is
c            multiplied by delib )
c            ad,bd,cd calculations
c	     d=a*n**4+b*n**2+c
c
c  if 3 begin
      if (ib.gt.1) then
         xe=wpw_2(x,y,z,1)
         ye=wcw(x,y,z,1)
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
         da1bdc=-(s1-peym2)+s4
         da0bdc=pbyp
         db1bdc=-s4*(s1-peym2)+(s2-peyp)*(s3-peym)
         db0bdc=s4*pbyp-xib*(s3-peym)
         dc1bdc=0.d0
         dc0bdc=0.d0
         dadc2=delib*da1bdc+da0bdc
         dbdc2=delib*db1bdc+db0bdc
         dcdc2=delib*dc1bdc+dc0bdc
      endif ! ib>1
c  if 3 end
c-----------------------------------------------------------------

c  derivatives calculations
c   dadz,dadx,dady
c   dbdz,dbdx,dbdy
c   dcdz,dcdx,dcdy
c   dadw,dbdw,dcdw
c----------------------------------------------------------------
      dadz=0.d0
      dadx=0.d0
      dady=0.d0
      dbdz=0.d0
      dbdx=0.d0
      dbdy=0.d0
      dcdz=0.d0
      dcdx=0.d0
      dcdy=0.d0
      dadw=0.0d0
      dbdw=0.0d0
      dcdw=0.0d0

c  if nbulk.gt.1  (ions components are in plasma )
c
c  if 4 begin
      if (nbulk.gt.1) then
      
        do 10 i=2,nbulk
c
c  derivatives calculations
c   dadx,dady
c   dbdx,dbdy
c   dcdx,dcdj

          call dwpw_2(x,y,z,i, dwpw2dx,dwpw2dy,dwpw2dz) ! derivs of (wp/w)^2
          call dwcw(x,y,z,i, dwcwdx,dwcwdy,dwcwdz) ! derivs of wc/w

          if ( i.eq.ib) goto 20
c
c  i.ne.ib
c
          xi=wpw_2(x,y,z,i)
          yi=wcw(x,y,z,i)

          pyim=1.d0/(1.d0-yi)
          pyip=1.d0/(1.d0+yi)
          pyim2=pyim*pyim

c if 5 begin
          if (ib.gt.1) then
             ds1dxi=pyim*pyip
             ds2dxi=pyim
             ds1dyi=2.d0*xi*yi*ds1dxi*ds1dxi
             ds2dyi=xi*pyim2
          end if
c if 5 end
          ds3dxi=pyip
          ds4dxi=1.d0
          ds3dyi=-xi*pyip*pyip
          ds4dyi=0.d0
c if 6 begin
          if (ib.eq.1) then
             ds6dxi=pyim
             ds7dxi=pyim*pyip
             ds6dyi=xi*pyim2
             ds7dyi=2.d0*xi*yi*ds7dxi*ds7dxi
          end if
c if 6 end

c if 7 begin
          if (ib.eq.1) then
             da1exi=-ds7dxi*ds2-ds4dxi*dc2
             da0exi=0.d0
             db1exi=(ds4dxi*s7+ds7dxi*s4)*pc+
     1		      (ds3dxi*(s6-xe/(1.d0+ye))+ds6dxi*s3)*ds2
             db0exi=-ds4dxi*xe/(1.d0+ye)*pc-
     1                 xe*ds6dxi*ds2
             dc1exi=-ds4dxi*s3*(s6-xe/(1.d0+ye))-
     1                 s4*ds3dxi*(s6-xe/(1.d0+ye))-s4*s3*ds6dxi
             dc0exi=xe*ds4dxi*(s6-xe/(1.d0+ye))+xe*s4*ds6dxi

             da1eyi=-ds7dyi*ds2-ds4dyi*dc2
             da0eyi=0.d0
             db1eyi=(ds4dyi*s7+ds7dyi*s4)*pc+
     1                (ds3dyi*(s6-xe/(1.d0+ye))+ds6dyi*s3)*ds2
             db0eyi=-ds4dyi*xe/(1.d0+ye)*pc-xe*ds6dyi*ds2
             dc1eyi=-ds4dyi*s3*(s6-xe/(1.d0+ye))-
     1    	       ds3dyi*s4*(s6-xe/(1.d0+ye))-ds6dyi*s4*s3
             dc0eyi=xe*ds4dyi*(s6-xe/(1.d0+ye))+ds6dyi*xe*s4

             dadxi=dele*da1exi+da0exi
             dbdxi=dele*db1exi+db0exi
             dcdxi=dele*dc1exi+dc0exi
             dadyi=dele*da1eyi+da0eyi
             dbdyi=dele*db1eyi+db0eyi
             dcdyi=dele*dc1eyi+dc0eyi
             goto 30
          end if
c if 7 end

c if 8 begin
          if (ib.gt.1) then
             da1bxi=-ds1dxi*ds2-ds4dxi*dc2
             da0bxi=0.d0
             db1bxi=(ds4dxi*(s1-xe/(1.d0-ye*ye))+ds1dxi*s4)*pc+
     1		  (ds2dxi*(s3-xe/(1.d0-ye))+ds3dxi*(s2-xe/(1.d0+ye)))*
     2                 ds2
             db0bxi=-ds4dxi*xib/(1.d0+yib)*pc-
     1                 xib*ds3dxi*ds2
             dc1bxi=-ds4dxi*(s2-xe/(1.d0+ye))*(s3-xe/(1.d0-ye))-
     1                 ds2dxi*s4*(s3-xe/(1.d0-ye))-
     2                 ds3dxi*s4*(s2-xe/(1.d0+ye))
             dc0bxi=xib*ds4dxi*(s3-xe/(1.d0-ye))+xib*s4*ds3dxi

             da1byi=-ds1dyi*ds2-ds4dyi*dc2
             da0byi=0.d0
             db1byi=(ds4dyi*(s1-xe/(1.d0-ye*ye))+s4*ds1dyi)*pc+
     1          (ds2dyi*(s3-xe/(1.d0-ye))+ds3dyi*(s2-xe/(1.d0+ye)))*
     2                 ds2
             db0byi=-ds4dyi*xib/(1.d0+yib)*pc-xib*ds3dyi*ds2
             dc1byi=s4*(-ds2dyi*(s3-xe/(1.d0-ye))-
     1    	       ds3dyi*(s2-xe/(1.d0+ye)))
             dc0byi=xib*s4*ds3dyi

             dadxi=delib*da1bxi+da0bxi
             dbdxi=delib*db1bxi+db0bxi
             dcdxi=delib*dc1bxi+dc0bxi
             dadyi=delib*da1byi+da0byi
             dbdyi=delib*db1byi+db0byi
             dcdyi=delib*dc1byi+dc0byi
             goto 30
          end if
c if 8 end
 20       continue
c
c         i=ib, i.ne.1
c

c if 15 begin
          if (ib.gt.1) then
             ds1dxb=0.d0
             ds2dxb=0.d0
             ds1dyb=0.d0
             ds2dyb=0.d0
          end if
c if 15 end
          ds3dxb=1.d0/(1.d0+yib)
          ds4dxb=1.d0
          ds3dyb=-xib/(1.d0+yib)**2
          ds4dyb=0.d0
c if 16 begin
          if (ib.eq.1) then
             ds6dxb=0.d0
             ds7dxb=0.d0
             ds6dyb=0.d0
             ds7dyb=0.d0
          end if
c if 16 end

c if 9 begin
          if (ib.gt.1) then
            da1bxb=-ds1dxb*ds2-ds4dxb*dc2
            da0bxb=-1.d0/(1.d0+yib)*ds2
            db1bxb=(ds4dxb*(s1-xe/(1.d0-ye*ye))+ds1dxb*s4)*pc+
     1            (ds2dxb*(s3-xe/(1.d0-ye))+ds3dxb*(s2-xe/(1.d0+ye)))*
     2                 ds2
            db0bxb=(-ds4dxb*xib/(1.d0+yib)+s4/(1.d0+yib))*pc+
     1                ((s3-xe/(1.d0-ye))-xib*ds3dxb)*ds2
            dc1bxb=-ds4dxb*(s2-xe/(1.d0+ye))*(s3-xe/(1.d0-ye))-
     1                 ds2dxb*s4*(s3-xe/(1.d0-ye))-
     2                 ds3dxb*s4*(s2-xe/(1.d0+ye))
            dc0bxb=-s4*(s3-xe/(1.d0-ye))+xib*ds4dxb*(s3-xe/(1.d0-ye))+
     +                xib*s4*ds3dxb
            da1byb=0.d0
            da0byb=xib/(1.d0+yib)**2*ds2
            db1byb=(ds4dyb*(s1-xe/(1.d0-ye*ye))+s4*ds1dyb)*pc+
     1          (ds2dyb*(s3-xe/(1.d0-ye))+ds3dyb*(s2-xe/(1.d0+ye)))*
     2                 ds2
            db0byb=(-ds4dyb*xib/(1.d0+yib)-s4*xib/(1.d0+yib)**2)*pc-
     -                xib*ds3dyb*ds2
            dc1byb=s4*(-ds2dyb*(s3-xe/(1.d0-ye))-
     1                 ds3dyb*(s2-xe/(1.d0+ye)))
            dc0byb=xib*s4*ds3dyb
            dadxi=delib*da1bxb+da0bxb
            dbdxi=delib*db1bxb+db0bxb
            dcdxi=delib*dc1bxb+dc0bxb
            dadyi=delib*da1byb+da0byb-a1b
            dbdyi=delib*db1byb+db0byb-b1b
            dcdyi=delib*dc1byb+dc0byb-c1b
          end if
c if 9 end

 30       continue

          dadz=dadz+dadxi*dwpw2dz+dadyi*dwcwdz
          dadx=dadx+dadxi*dwpw2dx+dadyi*dwcwdx
          dady=dady+dadxi*dwpw2dy+dadyi*dwcwdy

          dbdz=dbdz+dbdxi*dwpw2dz+dbdyi*dwcwdz
          dbdx=dbdx+dbdxi*dwpw2dx+dbdyi*dwcwdx
          dbdy=dbdy+dbdxi*dwpw2dy+dbdyi*dwcwdy

          dcdz=dcdz+dcdxi*dwpw2dz+dcdyi*dwcwdz
          dcdx=dcdx+dcdxi*dwpw2dx+dcdyi*dwcwdx
          dcdy=dcdy+dcdxi*dwpw2dy+dcdyi*dwcwdy

          xi=wpw_2(x,y,z,i)
          yi=wcw(x,y,z,i)
c*************************************c
          dadw=dadw-2*dadxi*xi-dadyi*yi
          dbdw=dbdw-2*dbdxi*xi-dbdyi*yi
          dcdw=dcdw-2*dcdxi*xi-dcdyi*yi
c new--------------------------------------

 10     continue
 
      end if
c if 4 nbulk > 1 end


c
c     i = 1
c
        ds1dxe=0.d0
        ds2dxe=0.d0
        ds3dxe=0.d0
        ds4dxe=1.d0
        ds6dxe=0.d0
        ds7dxe=0.d0
        ds1dye=0.d0
        ds2dye=0.d0
        ds3dye=0.d0
        ds4dye=0.d0
        ds6dye=0.d0
        ds7dye=0.d0

c if 10 begin
c
c   ib > 1, i = 1
c
        if (ib.gt.1) then
           da1bxe=-1.d0/(1.d0-ye*ye)*ds2-ds4dxe*dc2
           da0bxe=0.d0
           db1bxe=(ds4dxe*(s1-xe/(1.d0-ye*ye))+s4/(1.d0-ye*ye))*pc+
     1                (1.d0/(1.d0+ye)*(s3-xe/(1.d0-ye))+
     +                 1.d0/(1.d0-ye)*(s2-xe/(1.d0+ye)))*ds2
           db0bxe=-ds4dxe*xib/(1.d0+yib)*pc-xib/(1.d0-ye)*ds2
           dc1bxe=s4*(-1.d0/(1.d0+ye)*(s3-xe/(1.d0-ye))-
     -                     1.d0/(1.d0-ye)*(s2-xe/(1.d0+ye)))-
     -                 ds4dxe*((s2-xe/(1.d0+ye))*(s3-xe/(1.d0-ye)))
           dc0bxe=xib*s4/(1.d0-ye)+xib*ds4dxe*(s3-xe/(1.d0-ye))

           da1bye=-2.d0*xe*ye/(1.d0-ye*ye)**2*ds2-ds4dye*dc2
           da0bye=0.d0
           db1bye=s4*2.*xe*ye/(1.d0-ye*ye)**2*pc+
     +                (-xe/(1.d0+ye)**2*(s3-xe/(1.d0-ye))+
     +                  xe/(1.d0-ye)**2*(s2-xe/(1.d0+ye)))*ds2
           db0bye=-xib*xe/(1.d0-ye)**2*ds2
           dc1bye=s4*(xe/(1.d0+ye)**2*(s3-xe/(1.d0-ye))-
     -                    xe/(1.d0-ye)**2*(s2-xe/(1.d0+ye)))
           dc0bye=xib*s4*xe/(1.d0-ye)**2
           dadxi=delib*da1bxe+da0bxe
           dbdxi=delib*db1bxe+db0bxe
           dcdxi=delib*dc1bxe+dc0bxe
           dadyi=delib*da1bye+da0bye
           dbdyi=delib*db1bye+db0bye
           dcdyi=delib*dc1bye+dc0bye
           goto 40
        end if ! ib>1
c if 10 end

c if 11 begin
c
c   ib = 1, i = 1
c
        if (ib.eq.1) then
           da1exe=-ds4dxe*dc2
           da0exe=-1.d0/(1.d0+ye)*ds2
           db1exe=ds4dxe*s7*pc+s3/(1.d0+ye)*ds2
           db0exe=(-ds4dxe*xe/(1.d0+ye)+s4/(1.d0+ye))*pc+
     +                 (s6-2.d0*xe/(1.d0+ye))*ds2
           dc1exe=-ds4dxe*s3*(s6-xe/(1.d0+ye))-
     -                s4*s3/(1.d0+ye)
           dc0exe=-s4*(s6-2.d0*xe/(1.d0+ye))+xe*(s6-xe/(1.d0+ye))
           da1eye=0.d0
           da0eye=xe/(1.d0+ye)**2*ds2
           db1eye=-s3*xe/(1.d0+ye)**2*ds2
           db0eye=-s4*xe/(1.d0+ye)**2*pc+xe*xe/(1.d0+ye)**2*ds2
           dc1eye=-ds4dye*s3*(s6-xe/(1.d0+ye))+s4*s3*xe/(1.d0+ye)**2
           dc0eye=-xe*xe*s4/(1.d0+ye)**2
           dadxi=dele*da1exe+da0exe
           dbdxi=dele*db1exe+db0exe
           dcdxi=dele*dc1exe+dc0exe
           dadyi=dele*da1eye+da0eye-a1e
           dbdyi=dele*db1eye+db0eye-b1e
           dcdyi=dele*dc1eye+dc0eye-c1e
        end if ! ib=1
c if 11 end

 40     continue
        call dwpw_2(x,y,z,1, dwpw2dx,dwpw2dy,dwpw2dz) ! derivs of (wp/w)^2
        call dwcw(x,y,z,1, dwcwdx,dwcwdy,dwcwdz) ! derivs of wc/w
        dadz=dadz+dadxi*dwpw2dz+dadyi*dwcwdz+dadc2*dc2dz
        dadx=dadx+dadxi*dwpw2dx+dadyi*dwcwdx+dadc2*dc2dx
        dady=dady+dadxi*dwpw2dy+dadyi*dwcwdy+dadc2*dc2dy

        dbdz=dbdz+dbdxi*dwpw2dz+dbdyi*dwcwdz+dbdc2*dc2dz
        dbdx=dbdx+dbdxi*dwpw2dx+dbdyi*dwcwdx+dbdc2*dc2dx
        dbdy=dbdy+dbdxi*dwpw2dy+dbdyi*dwcwdy+dbdc2*dc2dy

        dcdz=dcdz+dcdxi*dwpw2dz+dcdyi*dwcwdz+dcdc2*dc2dz
        dcdx=dcdx+dcdxi*dwpw2dx+dcdyi*dwcwdx+dcdc2*dc2dx
        dcdy=dcdy+dcdxi*dwpw2dy+dcdyi*dwcwdy+dcdc2*dc2dy

        xi=wpw_2(x,y,z,1)
        yi=wcw(x,y,z,1)
        dadw=dadw-2*dadxi*xi-dadyi*yi
        dbdw=dbdw-2*dbdxi*xi-dbdyi*yi
        dcdw=dcdw-2*dcdxi*xi-dcdyi*yi

        dcn=dsqrt(cnz*cnz+cnx*cnx+cny*cny)
        dcn2=dcn*dcn
        dcn4=dcn2*dcn2

        dadnz=dadc2*dc2dnz
        dbdnz=dbdc2*dc2dnz
        dcdnz=dcdc2*dc2dnz
        dadnx=dadc2*dc2dnx
        dbdnx=dbdc2*dc2dnx
        dcdnx=dcdc2*dc2dnx

        dadny=dadc2*dc2dny
        dbdny=dbdc2*dc2dny
        dcdny=dcdc2*dc2dny

c  dispersion relation a*n**4+b*n**2+c=0
        if (id.eq.1) then
          dddz=dcn4*dadz+dcn2*dbdz+dcdz
          dddx=dcn4*dadx+dcn2*dbdx+dcdx
          dddy=dcn4*dady+dcn2*dbdy+dcdy
          dddn=4.d0*ad*dcn2+2.d0*bd
          dddcnz=dddn*cnz+dcn4*dadnz+dcn2*dbdnz+dcdnz
          dddcnx=dddn*cnx+dcn4*dadnx+dcn2*dbdnx+dcdnx
          dddcny=dddn*cny+dcn4*dadny+dcn2*dbdny+dcdny
          dddw=dcn4*dadw+dcn2*dbdw+dcdw
     1         -4.d0*dcn4*ad-2.d0*dcn2*bd
          dddw=dddw/wf
        end if ! id=1

c  dispersion relation n**2-(-b+iom*sqrt(b*2-4*a*c))/(2*a)=0
        if (id.eq.2) then
          det=dsqrt(bd*bd-4.d0*ad*cd)
          p4=0.5d0/(ad*ad)
          p5=1.d0/det
          p6=-bd+ioxm*det
          p7=ioxm*(bd*dbdz-2.d0*dadz*cd-2.d0*ad*dcdz)*p5
          dddz=-((-dbdz+p7)*ad-dadz*p6)*p4
          p7=ioxm*(bd*dbdx-2.d0*dadx*cd-2.d0*ad*dcdx)*p5
          dddx=-((-dbdx+p7)*ad-dadx*p6)*p4
          p7=ioxm*(bd*dbdy-2.d0*dady*cd-2.d0*ad*dcdy)*p5
          dddy=-((-dbdy+p7)*ad-dady*p6)*p4

          p7=ioxm*(bd*dbdnz-2.d0*dadnz*cd-2.d0*ad*dcdnz)*p5
          dddcnz=2.d0*cnz-
     1      	 ((-dbdnz+p7)*ad-dadnz*p6)*p4
          p7=ioxm*(bd*dbdnx-2.d0*dadnx*cd-2.d0*ad*dcdnx)*p5
          dddcnx=2.d0*cnx-
     1      	 ((-dbdnx+p7)*ad-dadnx*p6)*p4
          p7=ioxm*(bd*dbdny-2.d0*dadny*cd-2.d0*ad*dcdny)*p5
          dddcny=2.d0*cny-
     1      	 ((-dbdny+p7)*ad-dadny*p6)*p4
          p7=ioxm*(bd*dbdw-2.d0*dadw*cd-2.d0*ad*dcdw)*p5
          dddw=-((-dbdw+p7)*ad-dadw*p6)*p4
     1	      -2.d0*dcn2
          dddw=dddw/wf
        end if ! id=2
        goto 50
        
      end if
c   end of cold plasma  dispersion

c-----------------------------------------------------------------
c     hot non-relativistic plasma dispersion from Forest code
      if (id.eq.6)then
          call hotdervs_xyz(u,wf,dddcnx,dddcny,dddcnz,
     .                           dddx,dddy,dddz,dddw)
          goto 50
      end if
c-----------------------------------------------------------------

 50   continue
 
      deru(1)=dddcnx
      deru(2)=dddcny
      deru(3)=dddcnz
      deru(4)=-dddx
      deru(5)=-dddy
      deru(6)=-dddz
      dddw=-dddw*wf ! wf=frqncy
      if (i_geom_optic.eq.1) then
         deru(1)=deru(1)/dddw
         deru(2)=deru(2)/dddw
         deru(3)=deru(3)/dddw
         deru(4)=deru(4)/dddw
         deru(5)=deru(5)/dddw
         deru(6)=deru(6)/dddw
      else
         if(i_geom_optic.eq.2) then
            p=1.d0/dsqrt(deru(1)**2+deru(2)**2+deru(3)**2)
            p=ray_direction*p
            deru(1)=p*deru(1)
            deru(2)=p*deru(2)
            deru(3)=p*deru(3)
            deru(4)=p*deru(4)
            deru(5)=p*deru(5)
            deru(6)=p*deru(6)
         endif
      endif

1953  continue ! to skip idif.eq.1 in the above, when idif.eq.2 was used

c           write(*,'(a,3e14.5)')'rside_xyz: R, rho, dddw=', R,rho,dddw
c           write(*,'(a,3e12.3)')'rside_xyz: dddcnx,dddcny,dddcnz',
c     +                           dddcnx,dddcny,dddcnz
c           write(*,'(a,3e12.3)')'rside_xyz      : dddx,dddy,dddz',
c     +                           dddx,dddy,dddz
c           write(*,'(a,3e12.3)')'rside_xyz Vgr/c: Vx/c,Vy/c,Vz/c',
c     +                           deru(1),deru(2),deru(3)
           !pause !!!

      !YuP[04-2016] Set limitations for Vgroup/c components
      !deru(1)=max(deru(1),-1.0)
      !deru(1)=min(deru(1),+1.0)
      !deru(2)=max(deru(2),-1.0)
      !deru(2)=min(deru(2),+1.0)
      !deru(3)=max(deru(3),-1.0)
      !deru(3)=min(deru(3),+1.0)
c      Vgr_c_2= deru(1)**2+deru(2)**2+deru(3)**2 ! |Vgr/c|^2
c      Vgr_c= dsqrt(Vgr_c_2)
c      if(Vgr_c.ge.1.001)then 
         ! Vgr/c >1    - usually it happens at close proximity 
         ! of a cyclotron resonance (for EBW). 
         ! As ray approaches the res.layer, Vgr drops to nearly zero,
         ! and then it can jump (numerical instability?)
c        hp=  hamilt_xyz(x,y,z, cn2) ! test: is it large?
c        write(*,'(a,3e12.3)')' END OF rside_xyz:  x,hp, (Vgroup/c)^2=',
c     +   x,hp,Vgr_c_2
        ! Rescale components of Vgr/c to keep it below 1.0:
        !scale_vgr= 1.0/Vgr_c
        !deru(1)=deru(1)*scale_vgr
        !deru(2)=deru(2)*scale_vgr
        !deru(3)=deru(3)*scale_vgr
        ! Now |Vgr/c| is 1.0
        ! This is a crude method, and may result in unstable 
        ! behaviour of ray. But it forces it to be absorbed more fully.
c      endif
      
      return ! rside_xyz
      end
      
c======================================================================
c======================================================================

c        4_th order  Runge-Kutta method with automatic 
c        time step selection
c
c        **********************  drkgs2_xyz  ****************
c        *                      -----                       
c        * this Runge-Kutta subroutine finds the solution of the 
c        *      system of ordinary differential equations          
c        ****************************************************
c
c        drkgs2_xyz: called if(isolv=1 & irkmeth=2)
!  irkmeth=2: Most developed method of ray equation integration.
!             Time or space step in the equations
!             (according to setting of i_geom_optic) is controlled
!             so that output is at intervals prmt6 (meters).
!             As ray approaches the plasma edge, it is reflected
!             at the last closed flux surface.
!     *** NOTE: irkmeth=2 works best for fully relativistic dispersion
!        relations (id=11,12,14,15).
c
c The set of equations to be solved (depends on i_geom_optic):
!  i_geom_optic=1  Integration in time (default):
!                  ray-tracing equations right hand side=
!                  dr^/dt=-(dD/dN^)/(dD/domega)
!                  dN^/dt=+(dD/dr^)/(dD/domega)
!                  In this case deru(1:3) gives v_group (normalized).
c                  Note: dt here is NOT an actual physical time, but
c                  normalized as dt_code= dt_physical[sec]*c/omega  
c                                               
!  i_geom_optic=2  Integration is space,
!                  ray-tracing equations right hand side=
!                  dr^/dl=- ray_direction * (dD/dN^)p
!                  dN^/dl=  ray_direction * (dD/dr^)p
!                  p=1.d0/dsqrt(deru(1)**2+deru(2)**2+deru(3)**2)=
!                  deru(1)=dD/dNx,deru(2)=dD/dNy,deru(3)=dD/dNz,
c-----------------------------------------------------------------
c                                                                !
c     1. method of solution                                    !
c                                                                !
c             runge-kutta method  of  4th-order                  !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
c                                                                !
c     2. input parameters                                      !
c                                                                !
c     prmt  - the vector of input (and output) data:          !
c     prmt1=prmt(1)= initial time for a ray ! Not needed. 0 by default.
c     prmt2=prmt(2)= largest allowed time for ray advancing  (Not used)
c     prmt3=prmt(3)= initial time step for integration (normalized units)
c     prmt4=prmt(4)=required accuracy of solution;
c     !!!prmt(5)=0. Do not change: internal control:
!                   if it is not equal to zero the control  
!                   is transferred to the main program after one step.     
!     prmt(6)=prmt6 [m] distance step for saving ray data
!                   MAY AFFECT the h step of integration !
c     !!!prmt(7)=prmt(6) ! do not change: control logic
c     !!!prmt(8)=prmt(1)+prmt(3) ! NOT USED
c     prmt(9)=prmt9  ! accuracy of Hamiltonian
!                it will reduce Runge-Kutta time step 
!                if dabs(ham).ge.eps_ham)
!                 =1.d15 by default 
!                (for such big value the comparison
!                of the hamiltonian (hamilt<prmt9) in Runge-Kutta
!                practically does no change time step
!                and Runge-Kutta works as if no check is done)
!     An example of what prmt parameters work
!        well for relativistic EBW and irkmeth=2:
!                prmt3=prmt(3)=1.0d-04
!                prmt4=prmt(4)=5.0d-04
!                prmt6=prmt(6)=1.0d-03     [m]
c
c     2.1 Other input parameters                                      !
c          u()   - the input vector of initial values. later     !
c                  it becomes the vector of results obtained     !
c                  at the intermediate points u;                 !
c        
c                                                                !
c          deru() - the  vector of derivatives at point x,y,z	 !
c                                                                !
c          ndim  - the number of system equations;               !
c                                                                !
c          aux   - the auxiliary array.                          !
c                                                                !
c          i_output =1 the output at the poloidal distance steps !
c                   =2 the output at the total distance steps    !
c              (This selection sets the meaning of prmt6 control.)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
c                                                                !
c       3. subroutines and functions used                        !
c                                                                !
c          fct - calculates the right hand sides of  the  ODE    !
c                system. It has not to change the  values  of    !
c                u().     Its formal parameters are: u, deru.                                     !
c                                                                !
c          outp_xyz - output subroutine. it should not change the!
c                 values of all  its  formal  parameters.  
c              output data                                  
c              u(1)=x, u(2)=y, u(3)=z, u(4)=Nx, u(5)=Ny, u(6)=Nz                                         !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
      subroutine drkgs2_xyz(prmt,u,deru,ndim,ihlf,fct,outp_xyz,
     +  aux,i_output)
      implicit integer (i-n), real*8 (a-h,o-z)
      dimension u(6),deru(6),aux(8,6),prmt(9),up(6),uu(6),startu(6)
      integer i_output
      dimension uplus(6)
      double precision prmt,u,deru,aux,t,tend,h,up,hh,uu,startu,
     1 dd,uplus,eps4
      double precision us, ur,uphi, dels1,dels2, ux,uy,uz
      integer iflag
      
      external fct !=rside_xyz

      !!!t=prmt(1)    !drkgs2_xyz: (initial/start) time; not used.
      !!!tt= 0.d0 !t!   starting time not used
      !!!tend=prmt(2) !drkgs2_xyz:  final time; not used.
      h=prmt(3) ! initial value for h; will be adjusted below
      hh=h/2.0
      prmt(5)=0.d0 !if =1.d0, it will stop the ray after 1 step
      eps4=prmt(4)
      iflag=0  ! control of the RK method accuracy
      iflagh=3 ! control of the plasma edge intersection
      us=0.d0  ! distance or length along ray (initial value)
      dels2= 0.d0 !  to initialize
c---------------------------------------------------------------
 10   continue ! handle for time advances (in a loop)

      if (dabs(h).lt.1.d-11) then
         write(*,*)'***** In Runge-Kutta subroutine drkgs2_xyz *****'
         write(*,*)'***** time step is too small prmt3<1.e-11  *****' 
         write(*,*)'***** cannot achieve the given accuracy prmt4 **' 
         write(*,*)'h=prmt3=',h
         write(*,*)'  iraystop->1'
         iraystop=1
         goto 100 ! exit
      endif    

      ux=u(1) !=x
      uy=u(2) !=y
      uz=u(3) !=z
      do i=1,ndim
        uu(i)=u(i)
        startu(i)=u(i)
      enddo
      !!!tt=t
      hh=h/2.0d0
      call fct(u,deru)  !fct(t,u,deru) ! rhs of ODE
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
         u(i)=up(i) ! new u() is obtained
      enddo
c=====================================================================
      !for accuracy control: two steps starting with old uu 
      !but using h/2
      do 320 j=1,2 
      ! u->uu was saved before step 1
      call fct(uu,deru) !fct(tt,uu,deru)
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
 320  continue
      ! Done: uu() after 2steps using h/2 is obtained
!So, (ux,uy,uz) is before step, u(i) - after 1 step with h, 
!                               uu(i) - after 2 steps with h/2

c=====================================================================
      dd=0.0d0
      
      ru=  sqrt(u(1)**2  + u(2)**2  + u(3)**2 ) ! |r| after 1step with h
      ruu= sqrt(uu(1)**2 + uu(2)**2 + uu(3)**2) ! |r| after 2steps with h/2
      rave=0.5*(ru+ruu)
      dru=sqrt( (u(1)-uu(1))**2 +(u(2)-uu(2))**2 +(u(3)-uu(3))**2 )
      !write(*,'(a,2e12.3)')'drkgs2_xyz: rave,dru', rave,dru
      dd= dru/rave ! relative error
      !YuP: this is not a good idea to use rave (=r radius)
      !as a reference value for dd, because r could be close to 0.
      !Better use a characteristic scale, like diameter of device
      !or radius r to ray starting point. Suggestion: rmax
      
      eps_ham =prmt(9)
      xu=  uu(1)
      yu=  uu(2)
      zu=  uu(3)
      cnxu=uu(4) !=Nx refr.index
      cnyu=uu(5) !=Ny
      cnzu=uu(6) !=Nz
      cn2u= cnxu*cnxu + cnyu*cnyu + cnzu*cnzu
      bmod=bxyz(xu,yu,zu) !-> get b and derivs of b
      gam= gamma1_xyz(xu,yu,zu, cnxu,cnyu,cnzu) ! uses bx,by,bz,bmod
      ham= hamilt_xyz(xu,yu,zu, cn2u) ! uses gam

      !write(*,*)'drkgs2_xyz: dd,eps4=', dd,eps4
      ! dd= dru/rave ! relative error;    eps_ham=prmt(9)
      if ((dd .lt. eps4).and.(dabs(ham).lt.eps_ham)) then 
         goto 189 ! accuracy achieved
         !(But still need to check that dd is not exactly zero)
      endif
         
      if (iflag .le. 0) then
        h=hh
        do i=1,ndim
           u(i)=startu(i)
        enddo
c         write(*,*)'!!!!!!!!!!!!!!!! back !!!!!!!!!!!!!!!'
        iflag=-1
        !write(*,*)'dd>eps4  goto 10',dd,eps4
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

      if ((dd .lt. eps4/16.d0) .and. (dabs(ham).lt.prmt(9))) then
        ! The time step might be too small (dd can be 0)
        ! Double the time step and try again:
        h=h*2.0d0
        do i=1,ndim
           uplus(i)=u(i)
        enddo
        do i=1,ndim
           u(i)=startu(i)
        enddo
c         write(*,*)'!!!!!!!!!!!!!!!! back++ !!!!!!!!!!!!!!!'
        iflag=1
        !write(*,*)'dd<eps4/16  h->2h, goto 10'
        goto 10
      endif

c=====================================================================

191   continue

c      t=t+h
      iflag=0
c-----------------------------------------------------------
      if (i_output.eq.1) then
c        poloidal distance squared
         ur= sqrt(ux**2+uy**2) ! r original
         ura=sqrt(u(1)**2+u(2)**2) ! r new
         dels2= (uz-u(3))*(uz-u(3)) + (ur-ura)*(ur-ura)
      endif

      if (i_output.eq.2) then
c        total distance squared
         dels2= (ux-u(1))**2 + (uy-u(2))**2 + (uz-u(3))**2
      endif

      dels1=dsqrt(dels2) ! distance (one step)
      usold=us ! accumulated  dist., from previous step
      us=us+dels1 ! accumulated distance   

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
        !write(*,*)'before drkgs0 hnew,r',hnew,sqrt(u(1)**2+u(2)**2)
        call drkgs0(hnew,u,deru,ndim,fct,aux) ! uses fct()=rside_xyz()
        ! check that the ray is within xyz-grid and within r-grid, 
        dxdt=deru(1)
        dydt=deru(2)
        dzdt=deru(3)
        call boundc_xyz(u(1),u(2),u(3),dxdt,dydt,dzdt, iboundc) 
      endif
c-----------------------------------------------------------
      if (i_output.eq.1) then
c        poloidal distance after one step using drkgs0
         ur= sqrt(ux**2+uy**2) ! r original
         ura=sqrt(u(1)**2+u(2)**2) ! r new
         dels2= (uz-u(3))*(uz-u(3)) + (ur-ura)*(ur-ura)
      endif

      if (i_output.eq.2) then
c        total distance
         dels2= (ux-u(1))**2 + (uy-u(2))**2 + (uz-u(3))**2
      endif

      dels1=dsqrt(dels2)

      us=usold
      us=us+dels1 ! the distance after one step with h=hnew
      if (i_delh.eq.0) goto 60
      goto 60  
      if(us.lt.(prmt(7)+delh*prmt(6))) then 
c--------hnew is too shot to jump past prmt(7)+delh*prmt(6).
c        Now we will increase hnew.
          hnew=0.5d0*(hold+hnew)
          goto 20      
       endif
 60   continue
 
      !========================================================
      call outp_xyz(us,u,deru,ihlf,ndim,prmt,iflagh,iraystop)
      !outpt_xyz, among other control, increments the time step
      !nstep_rk ; if nstep_rk is bigger than maxsteps_rk
      !it will stop the ray.
      !Among other staff, the outpt_xyz checks: if(us>prmt(7)),
      !which is initially prmt(7)=prmt(6), 
      !then save ray data and also increase prmt(7): add prmt(6), 
      !so the new value of prmt(7) will be used 
      !during next call of outpt_xyz.
      !This way, the data is saved each prmt(6) steps.
      !========================================================
      
c-----------------------------------------------------------
      if (iflagh.eq.1) then !if 1
c----------------------------------------------------------------
c        ray is near the plasma boundary after Runge-Kutta
c        procedure (it is after reflection)
c-----------------------------------------------------------------
      else
        if (iflagh.eq.2) then !if 2 ! NOT USED ?
c-------------------------------------------------------------------
c         Ray is outside  the plasma after the correction.
c         It is necessary to reduce the time step in the Runge-Kutta
c         procedure and to recalculate the array u(i)
c----------------------------------------------------------------------
          !!!t=t-h NOT USED
          h=h*0.5d0
          write(*,*)'ray is outside the plasma after correction'
          pause
c------------------------------------------------------------------
          if ((prmt(3)/h).lt.100.d0)then !if 2.1
             do i=1,ndim
                u(i)=startu(i)
             enddo
             goto 10
          else
             write(*,*)' Initial_step/adjusted_step: prmt(3)/h > 100'
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
         goto 100 !-> finish and stop the ray
      end if     
      if (prmt(5).gt.0) then ! handle to exit integration along ray
         !write(*,*)'prmt(5)>0  goto 100 finish',prmt(5)
         goto 100 !-> finish here
      endif
  30  format(1x,6(1x,e11.4))
c------------------------------------------------------------------
      goto 10 ! go back for another step along ray.
  100 continue
      return
      end ! end of drkgs2 ============================================
c======================================================================
c======================================================================


c======================================================================
c======================================================================

c        4_th order  Runge-Kutta method with automatic 
c        time step selection
c
c        **********************  drkgs_auto  ****************
c        *                      -----                       
c        * this Runge-Kutta subroutine finds the solution of the 
c        *      system of ordinary differential equations          
c        ****************************************************
c
c        drkgs_auto: called if(isolv=1 & irkmeth=3), for ixyz=1 only
c
c The set of equations to be solved :
!        Integration in time (default):
!        ray-tracing equations right hand side=
!        dr^/dt_code = -(dD/dN^)/(dD/domega)
!        dN^/dt_code = +(dD/dr^)/(dD/domega)
!        deru(1:3) gives v_group (normalized).
c        Note: dt_code here is NOT an actual physical time, but
c        normalized as dt_code= dt_physical[sec]*c/omega  
c                                               
c-----------------------------------------------------------------
c                                                                !
c       1. method of solution                                    !
c                                                                !
c          runge-kutta method  of  4th-order                  
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
c                                                                !
c       2. input parameters                                      !
c                         
c          dL_step
c          dN_step  
c     dL_step=1.d-3 ! [m]  max allowed change in cartesian coords.
c     dN_step=1.d-2 ! max allowed change in refraction index.
c     The code will set the time step h = dt_code for integration
c     in such a way that the change |dr| in configuration space
c     is not larger than dL_step, and also
c     the change in refr. index |N| is not larger than dN_step.
c          prmt(3)=initial step of integration (not really needed)  
c              Step h is set below, based on dL_step and dN_step.                                  
c          prmt(6)= distance step along ray [m] for saving data.
c                  Does NOT affect the step of integration 
!                  in this subroutine!
c
c          u()   - the input vector of initial values. later     !
c                  it becomes the vector of results obtained     !
c                  at the intermediate points u;                 !
c        
c          deru() - the  vector of derivatives at point x,y,z	 !
c                                                                !
c          ndim  - the number of system equations;               !
c                                                                !
c          aux   - the auxiliary array (working array).  
c                                                                !
c      
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
c                                                                !
c       3. subroutines and functions used                        !
c                                                                !
c          fct - calculates the right hand sides of  the  ode    !
c                system. It has not to change the  values  of    !
c                u().     Its formal parameters are: u, deru.                                     !
c                                                                !
c          outp_xyz - output subroutine. it should not change the!
c                 values of all its  formal  parameters.  
c                 output data:                                   !
c          u(1)=x, u(2)=y, u(3)=z, u(4)=Nx, u(5)=Ny, u(6)=Nz     !                                         !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
      subroutine drkgs_auto(prmt,u,deru,ndim,ihlf,fct,outp_xyz,
     +  aux,dL_step,dN_step)
      implicit integer (i-n), real*8 (a-h,o-z)
      dimension u(6),deru(6),aux(8,6),prmt(9),up(6),uu(6),startu(6)
      integer i_output
      dimension uplus(6)
      double precision prmt,u,deru,aux,t,tend,h,up,hh,uu,startu,
     1 dd,uplus,eps4
      double precision us, ur,uphi, dels1,dels2, ux,uy,uz
      integer iflag
      
      external fct !=rside_xyz

      h=prmt(3) !step of integration; not really needed here:
      !h is set below, based on dL_step and dN_step.
      !(But reserve for future, maybe will be used?)
      prmt(5)=0.d0 ! YuP: if =1.d0, it will stop the ray after 1 step
      iflag=0  ! control of the RK method accuracy
      iflagh=3 !drkgs_auto: initialize (control of the plasma edge intersection)
      us=0.d0  !cumulative distance along ray (initial value)
      dels2= 0.d0 !  to initialize
      
c---------------------------------------------------------------
 10   continue ! handle for time-advancing (loop)
      ! Save values from previous time step:
      ux=u(1) !=x
      uy=u(2) !=y
      uz=u(3) !=z
      do i=1,ndim
        uu(i)=u(i)
        startu(i)=u(i)
      enddo
      call fct(u,deru) ! rhs of ODE: get dr/dt and dN/dt
      !Set the time step.
      ! This is the speed of ray propagation (normalized Vgroup):
      drdt= sqrt(deru(1)*deru(1)+deru(2)*deru(2)+deru(3)*deru(3))
      ! This is the "speed" of change of refr. index:
      dNdt= sqrt(deru(4)*deru(4)+deru(5)*deru(5)+deru(6)*deru(6))
      
      ! consider the inverse time steps 
      ! 1/dt ~ (dr/dt)/dL_step
      ! and 
      ! 1/dt ~ (dN/dt)/dN_step
      ! Select the larger of the two - it will correspond to the smaller 
      ! value of dt
      step_inv= max( drdt/dL_step , dNdt/dN_step )
      ! Why we need to consider the inverse time step:
      ! because one of these speeds, dr/dt or dN/dt, can be exactly 0
      ! or at least very close to 0.

      if(step_inv.gt.0.d0)then
        dt_code= 1.d0/step_inv !Normalized units: dt_code= dt[sec]*c/omega
        h= dt_code
      else ! both drdt and dNdt are 0 
           !(could be very small, below rounding error)
        write(*,*)'***** In Runge-Kutta subroutine drkgs_auto *****'
        write(*,*)'***** time step is too small.'
        write(*,*)'***** dr/dt, dN/dt=',drdt,dNdt
        write(*,*)'***** dL_step=',dL_step
        write(*,*)'***** dN_step=',dN_step
        write(*,*)'  iraystop->1'
        iraystop=1
        goto 100 ! exit
      endif
      ! The advantage of using both drdt/dL_step and dNdt/dN_step is:
      ! At plasma edge, for example, when starting as O-mode,
      ! the value of N is ~1, change of N is also small, dN ~1 
      ! (going down to 0 at O-X cutoff), while the group velocity 
      ! is large Vgr~c. In this case the time step is determined by 
      ! step_inv= drdt/dL_step, so that dt_code= dL_step/drdt.
      ! In opposite case of slow electrostatic waves, 
      ! the group velocity can be very low, so that drdt~0.
      ! On the other hand, the change in N can be very steep,
      ! dN>10 over a short travel distance.
      ! In this case the time step is determined by dNdt/dN_step term
      ! so that dt_code= dN_step/dNdt .
      
      do i=1,ndim	       
        aux(1,i)=h*deru(i)
        aux(5,i)=u(i)+0.5*aux(1,i)
        up(i)=aux(5,i)
      enddo

      call fct(up,deru)
      do i=1,ndim
         aux(2,i)=h*deru(i)
         aux(5,i)=u(i)+0.5*aux(2,i)
         up(i)=aux(5,i)
      enddo

      call fct(up,deru)
      do i=1,ndim
         aux(3,i)=h*deru(i)
         aux(5,i)=u(i)+aux(3,i)
         up(i)=aux(5,i)
      enddo
      
      call fct(up,deru)
      do i=1,ndim
        aux(4,i)=h*deru(i)
      enddo
      
      do i=1,ndim
        up(i)=u(i)+1.d0/6.d0*
     1        (aux(1,i)+2.d0*aux(2,i)+2.d0*aux(3,i)+aux(4,i))
      enddo
      ! Done: new u(i) after one step is obtained.
      do i=1,ndim
        u(i)=up(i)
      enddo
c--------------------------------------------------------------------
c     This section can be used for hamiltonian (dispersion D=0) control,
c     but we skip it for now.
c      xu=  u(1)
c      yu=  u(2)
c      zu=  u(3)
c      cnxu=u(4) !=Nx refr.index
c      cnyu=u(5) !=Ny
c      cnzu=u(6) !=Nz
c      cn2u= cnxu*cnxu + cnyu*cnyu + cnzu*cnzu
c      bmod=bxyz(xu,yu,zu) !-> get b and derivs of b
c      gam= gamma1_xyz(xu,yu,zu, cnxu,cnyu,cnzu) ! uses bx,by,bz,bmod
c      ham= hamilt_xyz(xu,yu,zu, cn2u) ! uses gam
c--------------------------------------------------------------------      

c-----------------------------------------------------------
c     total distance traveled in one step
      dels2= (ux-u(1))**2 + (uy-u(2))**2 + (uz-u(3))**2
      dels1= dsqrt(dels2)
      us= us+dels1 ! accumulate the distance with each step. 
      iflagh=3 !drkgs_auto/before outp_xyz : initialize
      !========================================================
      call outp_xyz(us,u,deru,ihlf,ndim,prmt,iflagh,iraystop)
      !outpt_xyz, among other control, increments the time step
      !nstep_rk ; if nstep_rk is bigger than maxsteps_rk
      !it will stop the ray.
      !Among other staff, the outpt_xyz checks: if(us>prmt(7)),
      !which is initially prmt(7)=prmt(6), 
      !then save ray data and also increase prmt(7): add prmt(6), 
      !so the new value of prmt(7) will be used 
      !during next call of outpt_xyz.
      !This way, the data is saved each prmt(6) steps.
      !========================================================
      !Ray can also be stopped if Vgroup/c >2, or 
      ! D/(N|gradD|) > toll_hamilt
      ! Also, the value of iflagh can be changed 3->1
      ! if(iflref.eq.1)
      if(iraystop.eq.1) then
         goto 100 ! finish and stop the ray
      endif
      
      goto 10  !drkgs_auto: another time step
      
  100 continue !drkgs_auto: exit handle
      return
      end ! end of drkgs_auto =========================================
c======================================================================
c======================================================================



      subroutine transmit_coef_ox_xyz(x,y,z, cn_x,cn_y,cn_z, transm_ox)
c     calculate transmission coefficient for OX mode conversion 
c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none
      include 'param.i' 
      include 'one.i'
c-----input
      real*8 r, x,y,z,         !the coordinates [m,m,radians]
     &cnr,cnz, cn_x,cn_y,cn_z   !refractive index coordinates

c-----output
      real*8 transm_ox ! transmission coefficient of OX mode conversion

c-----locals
      real*8 y_e,L_n,grad_x,cnpar,cnpol,freqncy
      real*8 dwpw2dx,dwpw2dy,dwpw2dz
   
c-----externals
      real*8 bxyz,bgrp2,wpw_2,wcw,transmission_ox 
     
      bmod=bxyz(x,y,z) !-> get bx,by,bz, and some derivatives
      y_e= wcw(x,y,z,1)
      
c-----L_n=wpw2/grad_X
      call dwpw_2(x,y,z,1, dwpw2dx,dwpw2dy,dwpw2dz) ! derivs of (wp/w)^2
      grad_x= dsqrt(dwpw2dx**2 + dwpw2dy**2 + dwpw2dz**2)

      L_n= wpw_2(x,y,z,1)/grad_x*100.d0        !cm

      cnpar= (cn_x*bx + cn_y*by + cn_z*bz)*o_bmod

      freqncy=frqncy                         ! GHZ
  
c-----Instead of N_poloidal we will use the refractive index
c     N^*[b*gradPsi]/|[b*gradPsi]|
c     [b*gradPsi] is a vector production
c
c     b^=e^r*br+e^phi*bphi+e^z*bz
c     gradPsi^=e^r*dpdrd+e^phi*0+e^z*dpdzd
c
c                 | e^r    e^phi   e^z  |
c     [b*gradPsi]=| br     bphi    bz   |=
c                 | dpdrd  0       dpdzd|
c     =e^r*bphi*dpdzd - e^phi(br*dpdzd-bz*dpdrd) - e^z*bphi*dpdrd
c
c     N_[b*gradPsi]=N^*[b*gradPsi]/|[b*gradPsi]|
     
      bgrp2= (by*dpdzd-bz*dpdyd)**2 + 
     +       (bz*dpdxd-bx*dpdzd)**2 + 
     +       (bx*dpdyd-by*dpdxd)**2
     
      cnpol= ( cn_x*(by*dpdzd-bz*dpdyd) + 
     +         cn_y*(bz*dpdxd-bx*dpdzd) + 
     +         cn_z*(bx*dpdyd-by*dpdxd)  ) / dsqrt(bgrp2)
     
      transm_ox=transmission_ox(cnpar,cnpol,y_e,L_n,freqncy)
c      write(*,*)'oxb.f in transmit_coef_ox transm_ox=',transm_ox
      return
      end

c======================================================================
c======================================================================

c/*
c        ********************** dwpw_2 **********************
c        *                      ----                        *
c        * this subroutine calculates the derivatives  of   *
c        * (omega_pl_i/omega)**2  with  respect to          *
c        * x,y,z cartesian coords.  (i - type of particles) *
c        ****************************************************
c
c------------------------------------------------------------------
c								   
c        input:					   
c      x,y,z - cartesian coordinates of the point  
c              where the derivatives are calculated.
c      i = 1 - electrons,                                          !
c        > 1 - ions.                                               !
c      rho is from common/one/  
c
c        output:
c      dwpw2dx,dwpw2dy,dwpw2dz   (derivatives)
c------------------------------------------------------------------
c         uses
c          v,dense0,denseb,rn1de,rn2de,rho,idens   from common 'one'
c           psilim,psimag            from common 'three'
c           functions drhopsi, ddnsdrho
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      subroutine dwpw_2(x,y,z,i, dwpw2dx,dwpw2dy,dwpw2dz)
c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none
      include 'param.i'
      include 'one.i'
      include 'fourb.i'
      include 'five.i' ! contains rmax
c-----input
      real*8 x,y,z, r !space coordinates
      integer i       !the plasma species number
c-----locals
      real*8  
     &  den,dn_drho,dro_dpsi,psi_xr, drhodx,drhody,drhodz, wpw2,
     + xth,yth, rho2, Rs0rr, rrk,Rs, dnrr_drho, dnes_drho, dnub_dr,
     + sech, arg, dndx,dndy, grad_n_est, r2, dnrr_dr, hh,wpw2p,wpw2m
c-----output:
      real*8 dwpw2dx,dwpw2dy,dwpw2dz
c-----externals
      real*8 drhopsi,ddnsrho,densrho,wpw_2,psif_xyz, dense_xyz,ias2r_Sm
      
      sech(arg)= 1.d0/cosh(arg) !2.d0/(exp(arg)+exp(-arg)) ! define sech

      if (model_rho_dens.eq.0 .or. model_rho_dens.eq.5) then
      ! For now, find the derivatives of (omegap/omega)^2 numerically,
      ! because sometimes it's difficult to define them analytically.
      ! But usually omegap is a smooth function of (x,y,z)
        hh=der_r ! before [11-2016]: hh=rmax*1.d-4
        wpw2p= wpw_2(x+hh,y,z,i)
        wpw2m= wpw_2(x-hh,y,z,i)
        dwpw2dx= (wpw2p-wpw2m)/(2.*hh)
        wpw2p= wpw_2(x,y+hh,z,i)
        wpw2m= wpw_2(x,y-hh,z,i)
        dwpw2dy= (wpw2p-wpw2m)/(2.*hh)
        wpw2p= wpw_2(x,y,z+hh,i)
        wpw2m= wpw_2(x,y,z-hh,i)
        dwpw2dz= (wpw2p-wpw2m)/(2.*hh)
        wpw2p= wpw_2(x,y,z,i) !call with original (x,y,z) to restore rho, etc.
        return
      endif

      !------------------------------------------------------------
      ! get rho based on model for density (model_rho_dens=1,2,4), 
      ! 2D spline of data (model_rho_dens=3),
      ! or based on magnetic flux (model_rho_dens=0,5):
      den=dense_xyz(x,y,z,i) !-> get rho (stored in one.i) 
      !------------------------------------------------------------
      
      if (model_rho_dens.eq.0) then ! rho is based on magnetic flux
      
        psi_xr=psif_xyz(x,y,z) ! spline at (r,z)
        dro_dpsi=drhopsi(psi_xr) ! not acurate !
        dn_drho=ddnsrho(rho,i) ! derivative density by rho; spline
         if(abs(dn_drho) .lt. 0.d0) then
           ! The |grad(n)| from spline was too small - 
           ! Probably because of loss of accuracy.
           dn_drho= -dense_xyz(0.d0,0.d0,0.d0,i) ! (peak density)/rho=1
         endif
        dwpw2dx= v(i)*dn_drho*dro_dpsi*dpdxd
        dwpw2dy= v(i)*dn_drho*dro_dpsi*dpdyd
        dwpw2dz= v(i)*dn_drho*dro_dpsi*dpdzd
        
      elseif(model_rho_dens.eq.1) then
      
        ! rho is defined on ellipsoid surfaces
        ! rho= sqrt(((x-elx0)/elax)^2 +((y-ely0)/elay)^2 +((z-elz0)/elaz)^2)
        if(rho.ne. 0.d0) then
          xth=  (x-elx0)*costt + (y-ely0)*sintt
          yth= -(x-elx0)*sintt + (y-ely0)*costt
          drhodx= ( xth*costt/elax**2 -yth*sintt/elay**2 )/rho
          drhody= ( xth*sintt/elax**2 +yth*costt/elay**2 )/rho
          drhodz= (z-elz0)/(rho*elaz*elaz)
          ! derivative of density by rho:
          ! dense_xyz= (dense0(i)-denseb(i))*den_scale(i)*
          !                  (1-rho**rn1de(i))**rn2de(i)+denseb(i)
          dn_drho= (dense0(i)-denseb(i))*den_scale(i)*
     *             (-rn1de(i)*rn2de(i)*rho**(rn1de(i)-1.))*
     *             (1.-rho**rn1de(i))**(rn2de(i)-1.)
          dwpw2dx= v(i)*dn_drho*drhodx
          dwpw2dy= v(i)*dn_drho*drhody
          dwpw2dz= v(i)*dn_drho*drhodz
        else  ! rho=0
          dwpw2dx= 0.d0
          dwpw2dy= 0.d0
          dwpw2dz= 0.d0
        endif ! rho
        
      elseif(model_rho_dens.eq.2) then
      
        ! rho is defined on ellipse surfaces in x-y; uniform in z
        ! rho= sqrt( ((x-elx0)/elax)^2 +((y-ely0)/elay)^2 )
        if(rho.ne. 0.d0) then
          xth=  (x-elx0)*costt + (y-ely0)*sintt
          yth= -(x-elx0)*sintt + (y-ely0)*costt
          drhodx= ( xth*costt/elax**2 -yth*sintt/elay**2 )/rho
          drhody= ( xth*sintt/elax**2 +yth*costt/elay**2 )/rho
          drhodz= 0.d0 
          rho2=rho*rho
          !===> Derivative of density by rho <===
          !-1-> Rigid Rotor profile:
          Rs0rr= max(elax,elay)
          rrk=(Rm0rr/Rs0rr)**2  ! K in Eq.(38)
          Rs=Rs0rr ! Assume no dependence in z, for now
          ! For dens_rr= dens0rr*( sech(rho2-rrk)/sech(rrk) )**2
          dnrr_drho= dens0rr*( -4*rho*tanh(rho2-rrk))*
     *                       ( sech(rho2-rrk)/sech(rrk) )**2
          !-2-> Ellipsoidal Spindle profile:
          ! For dens_es= dens0es* exp( -(Rs*rho-Rm0es)**2 / (2*rtau**2) )
          dnes_drho= dens0es* ( -(Rs*rho-Rm0es)*Rs/rtau**2 )*
     *               exp( -(Rs*rho-Rm0es)**2 / (2*rtau**2) )
          !-3-> Uniform background density profile:
          dnub_dr= 0.d0 ! within r < r_ub_edge
          r= sqrt(x*x+y*y)
          if (r .gt. r_ub_edge) then 
             ! Linear drop in region  r_ub_edge < r < wall_rmax
             ! For dens_ub= dens0ub*(1.d0 -(r-r_ub_edge)/(wall_rmax-r_ub_edge))
             dnub_dr= -dens0ub/(wall_rmax-r_ub_edge)
          endif
          dn_drho= dnrr_drho + dnes_drho ! Sum of RR and ES
c          if (r .gt. wall_rmax) then
c             ! Outside of chamber wall radius wall_rmax:
c             dn_drho=0.d0  !  Skip this option,
c             dnub_dr=0.d0  !  because of a jump of derivative
c          endif
          dwpw2dx= v(i)*( dn_drho*drhodx + dnub_dr*x/r )
          dwpw2dy= v(i)*( dn_drho*drhody + dnub_dr*y/r )
          dwpw2dz= v(i)*( dn_drho*drhodz ) 
        else  ! rho=0
          dwpw2dx= 0.d0
          dwpw2dy= 0.d0
          dwpw2dz= 0.d0
        endif ! rho

      elseif(model_rho_dens.eq.3) then
         ! (x,y)-spline based on dengrid(i,j) data from file.
         ! Note: in this model, dn/dz=0; uniform density in z-direction.
         ! Note: If the initial data is not accurate (not smooth),  
         ! the grad(n) of such profile can be quite bad !!!
         grad_n_est= (denmax-denmin)/xdenmax ! estimate of |grad(n)|
         r= sqrt(x*x+y*y)
         if(x.le.xdenmax .and. x.ge.xdenmin .and.
     +      y.le.ydenmax .and. y.ge.ydenmin ) then
           ! Find dn/dx:
           dndx=ias2r_Sm(tx_den,nx_den,ty_den,ny_den,
     +                  cxy_den,ncx_den,ncy_den,1,0,x,y,nx4a)
           ! Find dn/dy:
           dndy=ias2r_Sm(tx_den,nx_den,ty_den,ny_den,
     +                  cxy_den,ncx_den,ncy_den,0,1,x,y,nx4a)
         else ! out of grid
           dndx= -(x/r)*grad_n_est
           dndy= -(y/r)*grad_n_est
         endif
         if(abs(dndx)+abs(dndy) .lt. grad_n_est*0.01) then
           ! The |grad(n)| from spline was too small - 
           ! Probably because of loss of accuracy.
           dndx= -(x/r)*grad_n_est
           dndy= -(y/r)*grad_n_est
         endif
cc         write(*,'(a,4e12.3)')'x,y, dn/dx,dn/dy==',x,y,dndx,dndy
         dwpw2dx= v(i)*dndx
         dwpw2dy= v(i)*dndy
         dwpw2dz= 0.d0

      elseif(model_rho_dens.eq.4) then
          r2=x*x+y*y
          r= sqrt(r2)
          !===> Derivative of density by r <===
          !-1-> Rigid Rotor profile:
          arg= akappa*(2.*r2/rs_frc**2 -1.d0)
          ! For dens_rr= dens0rr*( sech(akappa*(2*r2/rs_frc**2 -1)) )**2
          dnrr_dr= -dens0rr*tanh(arg)*8*akappa*r*(sech(arg)/rs_frc)**2
          !-2-> Uniform background density profile:
          dnub_dr= 0.d0 ! within r < rs_frc
          if (r .gt. rs_frc) then 
             ! Linear drop in region  rs_frc < r < wall_rmax
             ! For dens_ub= dens0ub*(1.d0 -(r-rs_frc)/(wall_rmax-rs_frc))
             dnub_dr= -dens0ub/(wall_rmax-rs_frc)
          endif
          if(r.gt.0.d0)then 
            dwpw2dx= v(i)*(dnrr_dr + dnub_dr)*x/r 
            dwpw2dy= v(i)*(dnrr_dr + dnub_dr)*y/r 
          else ! r=0
            dwpw2dx= 0.d0 
            dwpw2dy= 0.d0 
          endif
          dwpw2dz= 0.d0 ! uniform in z

      else
      
        PRINT*,'Set model_rho_dens to 0 or 1 or 2 or 3 or 4 or 5'
        PRINT*,'model_rho_dens=1-5 works with ixyz=1 only (cartesian)'
        PRINT*,'model_rho_dens=0 works with ixyz=0 or 1'
        stop 'dwpw_2'
        
      endif ! model_rho_dens
      
      return
      END
c======================================================================
c======================================================================


c        ********************** dwcw ****************************
c        *                      ------                          *
c        * this subroutine calculates  the  derivatives  of     *
c        * (omega_cyclotron_i/omega) with  respect to  x,y,z    *
c        * (i - type of particles)                              *
c        ********************************************************
c
c------------------------------------------------------------------
c								   
c        input:					   
c      x,y,z - cartesian coordinates of the point  
c              where derivatives are calculated.      
c      i = 1 - electrons,                                         
c        > 1 - ions.    
c                                           
c        output:
c      dwcwdx,dwcwdy,dwcwdz   (derivatives)
c------------------------------------------------------------------
      subroutine dwcw(x,y,z,i, dwcwdx,dwcwdy,dwcwdz)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i' ! contains w(i),dbmdx,dbmdy,dbmdz
      ! Define derivatives of (omegac(i)/omega) over (x,y,z):
	dwcwdx= w(i)*dbmdx
	dwcwdy= w(i)*dbmdy
	dwcwdz= w(i)*dbmdz
	!write(*,*)'dwcw:',dwcwdx,dwcwdy,dwcwdz,rmax
      return
      end
c======================================================================
c======================================================================


      subroutine refractive_index_relative_error_xyz(x,y,z,
     +                            cn_x,cn_y,cn_z,  iraystop)
c-------------------------------------------------------------------
c     Calculate the relative error of delta_N/N=delta_n_devide_n
c     of the dispersion relation D(Nx,Ny,Nz)=0 at the given point
c     (x,y,z, Nx,Ny,Nz)
c     delta_N/N=(D/N)/|gradD|
c     Here
c     N=|N|=sqrt(Nx**2+Ny**2+Nz**2)
c     |gradD|=sqrt((dD/dNx)**2+(dD/dNy)**2+(dD/dNz)**2)
c
c     If delta_N/N=(D/N)/|gradD| > toll_hamilt it will put iraystop=1
c     else         iraystop=0
c     Variable toll_hamilt is set in one.i
c------------------------------------------------------------
      implicit none
      include 'param.i'
      include 'one.i'
c-----input
      real*8 x,y,z,   r,    !space coordinates
     *cn_x,cn_y,cn_z      !refractive index coordinates
c-----output
      integer iraystop 
        
c-----externals 
      real*8 bxyz,gamma1_xyz,hamilt_xyz, wpw_2

c-----locals
      real*8 u(6),deru(6),
     &grad_d,  ! |gradD in N space|
     &cn2,cn,d,     
     &delta_n_devide_n ! delta_N/N=(D/N)/|gradD|

      u(1) = x
      u(2) = y                                                    
      u(3) = z                                                  
      u(4) = cn_x                                                  
      u(5) = cn_y                                                 
      u(6) = cn_z
      
      cn2= cn_x*cn_x + cn_y*cn_y + cn_z*cn_z
      cn=  dsqrt(cn2)

      bmod=bxyz(x,y,z) !-> get b and derivs of b
      gam= gamma1_xyz(x,y,z, cn_x,cn_y,cn_z) ! uses bx,by,bz,bmod
      d= hamilt_xyz(x,y,z, cn2) ! uses gam

      call dddrz1_xyz(u,deru) ! Here, deru(1:3) is dHamiltonian/dN
      grad_d=dsqrt(deru(1)**2 +deru(2)**2 +deru(3)**2) !dHamiltonian/dN!

c      write(*,*)'refractive_index_relative_error_xyz: id,ib,ioxm=',
c     + id,ib,ioxm
       
      delta_n_devide_n= d/(cn*grad_d)

c      write(*,'(a,4e13.3)')'D/(N|gradD|)=',delta_n_devide_n, d,cn,grad_d
c      write(*,'(a,6e13.3)')'deru(1:6)=',deru(1:6)
      if( (delta_n_devide_n.gt.toll_hamilt) .and. (cn.gt.1.0)) then 
         !YuP[11-2016] Added: (cn.gt.1.0) condition
         ! Apply this check in case of N>1 only (away from cut-off)
         ! i.e., do not check the toll_hamilt condition
         ! if ray is close to the cut-off layer (where N can be ~0)
         write(*,*)
         write(*,*) '*************************************************'
         write(*,*)'refractive_index_relative_error_xyz:  iraystop->1'
         write(*,*)'  rho=',rho, '  Xe=',wpw_2(x,y,z,1)
         iraystop=1
      else
         iraystop=0
      endif

      if(iraystop.eq.1) then
         write(*,*) 'STOPPED because of D/(N|gradD|) > toll_hamilt'
         write(*,'(3(a,e12.4))') '  D=',d,'  N=',cn,'  grad_d=',grad_d
         write(*,*) 'D/(N|gradD|)=delta_n_devide_n',delta_n_devide_n
         write(*,*) 'toll_hamilt=',toll_hamilt
         write(*,*)'TRY TO INCREASE toll_hamilt in genray.in'
         write(*,*) '*************************************************'
         !write(*,*)
         !pause
      endif
   
      return
      end
      
      

c======================================================================
c======================================================================

c        **********************  cnxyz  **************************
c        This subroutine calculates the initial value		     *
c        cnx,cny,cnz                                             *
c        It directs the wave into or out the plasma              *
c        the input parameter i_vgr_ini is in common /one/        *
c        i_vgr_ini =+1 the wave is directed into the plasma      *
c                      (in the initial point)                    *
c                  =-1 the wave is directed out the plasma       *
c-----------------------------------------------------------------
c      * input: x,y,z,  cnpar, cnper_tang, cnper_norm,	      	
c________________________________________________________________
c      * output: cnx,cny,cnz -refractive index components
c-----------------------------------------------------------------
      subroutine cnxyz(x,y,z, cnpar, cnper_tang, cnper_norm,
     +                 cnx,cny,cnz) !->out
      implicit integer (i-n), real*8 (a-h,o-z)
c-----------------------------------------------------------------
      include 'param.i'
      include 'one.i'
      include 'grill.i'
      dimension u(6),deru(6)
      
c-----------------------------------------------------------------YuP
c     Calculate the components of unit vector  which is
c     perpendicular to b-field and TANGENTIAL to psi-flux-surface;
c     and the components of unit vector  which is
c     perpendicular to b-field and NORMAL to psi-flux-surface;
      call unit_vectors_perp_b(x,y,z, 
     +                     uvtx, uvty, uvtz,  uvnx, uvny, uvnz) !->out
c     Now find Nx,Ny,Nz 
      call nparper_to_nxyz(cnpar, cnper_tang, cnper_norm,
     +                     uvtx, uvty, uvtz,  uvnx, uvny, uvnz,
     +                     cnx,cny,cnz) !->out    Nx,Ny,Nz found
      write(*,'(a,3e12.3)')' cnxyz(1)after nparper_to_nxyz:',cnx,cny,cnz
c------------------------------------------------------------------YuP
c-------------------------------------------------------------------
c     initialization of cirho - the sign of Nperp_normal_to_surface
      cirho=1.d0 ! try +1 first, reverse if needed.
c------------------------------------------------------------------
10    continue

      call unit_vectors_perp_b(x,y,z, 
     +                     uvtx, uvty, uvtz,  uvnx, uvny, uvnz) !->out
c     Now find Nx,Ny,Nz 
      call nparper_to_nxyz(cnpar, cnper_tang, cnper_norm*cirho,
     +                     uvtx, uvty, uvtz,  uvnx, uvny, uvnz,
     +                     cnx,cny,cnz) !->out    Nx,Ny,Nz found
      write(*,'(a,3e12.3)')' cnxyz(2)after nparper_to_nxyz:',cnx,cny,cnz
      
c     I tried it for id=4 ECR case when rside1 used the derivatives
c     by the trajectory d/dl (not by time d/dt).
c     In some case d/dl runs the ray to the plasma edge.
c     As I understand this changing of the cnrho sign is not correct
c     fo EC launch.  

      if(((i_n_poloidal.eq.1).or.(i_n_poloidal.eq.2)).
     &    or.(istart.eq.1)) then  
c*****************************************************************
c        determination of cirho to obey the situation
c        when groop velocity is directed inside or outside the plasma
c        at the initial point
c-----------------------------------------------------------------
         u(1)=x
         u(2)=y
         u(3)=z
         u(4)=cnx
         u(5)=cny
         u(6)=cnz
         ! rside_xyz calls hamilt_xyz, uses ioxm !
         call rside_xyz(u,deru) !-> get deru() 
c----------------------------------------------------------------
c        cmultpl is the scalar multiplication V_group*grad(psi)
c        gradient(psi_code) is directed outside the plasma.
c        The grad(psi_code) is positive (pointing outside of plasma),
c        so V_group*grad(psi_code) should be negative.
c----------------------------------------------------------------
         cmultpl= dpdxd*deru(1) + dpdyd*deru(2) + dpdzd*deru(3)
         !write(*,'(6e11.2)')bx,bz,cnx,cnz,deru(1),deru(3)
         !YuP [02-2016] sometimes Vgroup is nearly perp to grad(psi_code).
         ! In this case do not make a reversal of Nrho.
         grpsi= dsqrt(dpdxd**2 + dpdyd**2 + dpdzd**2)   ! |grad(PSI)|
         vgroup=dsqrt(deru(1)**2+deru(2)**2+deru(3)**2) ! |Vgroup/c|
         if ( cmultpl.gt. 0.01*grpsi*vgroup ) then
            ! This condition is added on [Feb-2016] 
            ! If it is not satisfied, 
            ! i.e. (Vgroup.gradPSI) < 0.01*|Vgroup|*|gradPSI|
            ! which is an indication that 
            ! Vgroup is nearly perp to grad(PSI)
            ! then do not make a reversal of Nrho.
         if ( cmultpl*i_vgr_ini.gt.0.d0 ) then
            write(*,*)' cnxyz: cnpar, grpsi*(du/dt)==cmultpl=',
     +        cnpar,cmultpl
            write(*,*)' Vgroup= ',deru(1),deru(2),deru(3)
            write(*,*)' gradPSI=',dpdxd,  dpdyd,  dpdzd
            write(*,*)' Reversing Nrho (N_normal_to_surface)...'
c-----------the poloidal direction of the group velocity is opposite
c           to the direction determined by the parameter i_vgr_ini     
c           We change the sign of the poloidal group velocity 
            cirho=-1.d0
            go to 10
         end if
         endif
         
      endif !i_n_poloidal =1 or =2
c*****************************************************************
      !pause
      return
      end

c======================================================================
c======================================================================

c        ********************* cninit12_n_gam_xyz *************
c        *                        -                           *
c        * It solves the dispersion relation 
c        * N=N(n_par)=N(gam,ioxm)                             *
c        * for cold plasma for given ioxm                     *
c        * Then subroutine calculates the initial components  *
c        * of the refractive index  cnx,cny,cnz               *
c        ******************************************************
c
c------------------------------------------------------------------
c					                          !
c        input parameters				          !
c        x,y,z,cnpar,
c        cnper_tang - refractive index tangential to psi flux surface
c                     and perp to b-field (contributes to Nperp) 		       	  
c        istart from common block/ one/                           !
c        output parameters				          !
c        u(4)=cnx,u(5)=cny,u(6)=cnz 			          !
c        iraystop=1 end ray calculation                           !
c------------------------------------------------------------------
c        it uses the following functions and subroutines          !
c        ias1r,bxyz,wpw_2,wcw,gamma1_xyz,s_xyz,abc_xyz,hamilt_xyz !
c------------------------------------------------------------------
      subroutine cninit12_n_gam_xyz(x,y,z,cnpar,cnper_tang,
     1                  cnx,cny,cnz,iraystop)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
      save
      double complex cmplnper      
      !--------------
      write(*,*)'cninit12_n_gam_xyz: ioxm===',ioxm
      !pause
      
      iraystop=0
      pi=4*datan(1.d0)
      cnpar2=cnpar*cnpar
      cntang2= cnper_tang*cnper_tang + cnpar2 !cnteta*cnteta+cnphi*cnphi
      bmod=bxyz(x,y,z) ! get bmod for wcw
c-----------------------------------------------------------------
c------------------------------------------------------------------
c     epsmode is the accuracy  for selection mode (O or X) on
c             the plasma boundary
c------------------------------------------------------------------
      epsmode=1.d-3 !1.d-4 !1.d-7
      !cnx=0.d0
      !cny=0.d0
      !cnz=0.d0 ! Initialize [2016]
      ! WRONG: Do not initialize! 
      ! If no roots are found, at least cnx,xny,cnz will remain intact!
c if 0
      if ((id.eq.1).or.(id.eq.2)) then
c--------------------------------------------------------------------
c     initial condition at plasma boundary:
c     Find n^2 as a func. of npar by using dispersion relation
c     f*n^4 + g*n^2 + w = 0
c     Two roots are  n^2= [-g +/-sqrt(g**2-4*f*w)]/(2*f)
c------------------------------------------------------------------
c     cold plasma dispersion relation
c------------------------------------------------------------------
        call s_xyz(x,y,z,s1,s2,s3,s4,s6,s7)

        iroot=1  ! initialize; will be incremented

c       ib=1 electron resonance condition may be in plasma
c  if 2
        if (ib.eq.1) then
          xe=wpw_2(x,y,z,1)
          ye=wcw(x,y,z,1)
          pyp=xe/(1.d0+ye)
          dele=1.d0-ye
          f1e=s7
          f0e=-pyp
          g1e=cnpar2*(s4-s7)-(s6-pyp)*s3-s4*s7 ! Function of Npar
          g0e=cnpar2*pyp+xe*(s6-pyp)+pyp*s4    ! Function of Npar
          w1e=cnpar2*(-s7*s4+(s6-pyp)*s3)+(s6-pyp)*s3*s4 ! Func. of Npar
          w0e=cnpar2*(pyp*s4-xe*(s6-pyp))-s4*(s6-pyp)*xe ! Func. of Npar
          fd=f1e*dele+f0e
          gd=g1e*dele+g0e
          wd=w1e*dele+w0e
          detin=gd**2-4.d0*fd*wd
          if (detin.lt.0d0) then ! negative discriminant
             !YuP[11-2016] adjustment for a very small neg.discriminant
             if(abs(detin).lt. 1.d-8*(gd**2+4.d0*fd*wd) ) then
               ! Could be a very small negative value, from rounding.
               ! Then, set it to zero:
               detin=0.d0
             else ! large negative discriminant: no wave here.
               write(*,*)'cninit12_n_gam: detin<0. detin,gd**2,Xe=',
     +                                           detin,gd**2,xe
               write(*,*)'cninit12_n_gam: detin<0.    iraystop->1'
               !pause
               iraystop=1
               return
             endif
          end if
          cn2p=(-gd+dsqrt(detin))/(2.d0*fd)
          cn2m=(-gd-dsqrt(detin))/(2.d0*fd)

c--------- Check that N^2 is not smaller than tangential-to-surf Ntang^2
          irootp=0 ! initialize
          irootm=0 ! initialize
          if (cn2p.ge.cntang2)then
             dc_l=cnpar/dsqrt(cn2p)
             gam_l=dacos(dc_l)
             write(*,*)'cninit12_n_gam cn2p',cn2p
             irootp=1 ! means: cn2p is ok
          else
             write(*,*)'cninit12_n_gam cn2p<cntang2', cn2p,cntang2
          endif  

          if (cn2m.ge.cntang2)then
             dc_l=cnpar/dsqrt(cn2m)
             gam_l=dacos(dc_l)
             write(*,*)'cninit12_n_gam cn2m',cn2m
             irootm=1 ! means: cn2m is ok
          else
             write(*,*)'cninit12_n_gam cn2m<cntang2', cn2m,cntang2
          endif

          if(irootp+irootm .eq. 0) then ! none of them is good
             write(*,*)'cninit12_n_gam two roots of dispersion <cntang2'
             write(*,*)'cn2m,cn2p',cn2m,cn2p
             write(*,*)'cninit12_n_gam: WAVE CANNOT EXIST.  iraystop->1'
             iraystop=1
             return
          endif  

20        continue ! handle for two trials of roots
          write(*,*)'cninit12_n_gam: iroot=',iroot
          
	    if(irootp+irootm .eq.2)then !both cn2p and cn2m are good (>cntang2)
	      if(iroot.eq.1) then ! first trial
              if(ioxm.eq. 1) cn2=cn2p !first, try this root
              if(ioxm.eq.-1) cn2=cn2m !or this, depending on ioxm
	      else     !   iroot=2: second trial (max two trials)
              if(ioxm.eq. 1) cn2=cn2m !second, try this root
              if(ioxm.eq.-1) cn2=cn2p !or this, depending on ioxm
            endif
	    end if ! irootp+irootm .eq.2
	    
	    if(irootp+irootm .eq.1)then 
	      !only one of cn2p,cn2m was ok (>cntang2)
            if(irootp.eq. 1) cn2=cn2p ! This root was ok
            if(irootm.eq. 1) cn2=cn2m ! This root was ok
            iroot=2 ! to indicate that this is the last trial
	    end if ! irootp+irootm .eq.1

          cnrho2=cn2-cntang2
          cnrho=dsqrt(cnrho2)
          write(*,*)'cninit12_n_gam before cnxyz: cn2,cnrho=',cn2,cnrho
c-------------------------------------------------------------------
          cnper_norm=cnrho ! found above
          !Careful! cnxyz calls rside_xyz -> hamilt_xyz, which uses ioxm
          ! Here - call for ib=1 (electrons only):
          call cnxyz(x,y,z, cnpar, cnper_tang, cnper_norm,
     +                 cnx,cny,cnz) !->out
          cn2=cnx*cnx+cny*cny+cnz*cnz
          write(*,*)'cninit12_n_gam after cnxyz: 
     +               cn2,cnper_norm=',cn2,cnper_norm
c-----------------------------------------------------------------YuP          
          gam=gamma1_xyz(x,y,z,cnx,cny,cnz)
          ds=dsin(gam)
          dc=dcos(gam)
          ds2=ds*ds
          dc2=dc*dc
          ds4=ds2*ds2
c---------------------------------------------------------------------
c         control that cn2 and gam are the solution of the dispersion
c         relation  n**2=(-b+ioxm*sqrt(b**2-4*a*c))/(2*a)
c         This part below follows hamilt definition for id=2.
c---------------------------------------------------------------------
          call abc_xyz(x,y,z,ds2,dc2, ad,bd,cd)
          d4=ad
          d2=bd
          d0=cd
          det=d2*d2-4.d0*d4*d0
          cn1=dsqrt(d2*d2-4.d0*d4*d0)
          cn2new=  (-d2+ioxm*cn1)/(2.d0*d4) ! YuP: 
          write(*,*)'cninit12_n_gam after abc: cn2new  =',cn2new,ioxm
          dnp=dabs(cn2-cn2new)/dabs(cn2+cn2new)
          if (dnp.gt.epsmode) then ! Too far from ioxm-target-root
             if(iroot.eq.1) then
                write(*,*)'cninit12_n_gam: dnp>epsmode. Try another.'
                iroot=iroot+1 !->2
                goto 20 ! Try another
             else ! iroot=2: already all trials were made
                write(*,*)'cninit12_n_gam: cn2new far from cn2p or cn2m'
                write(*,*)'No roots N^2(gam) for this ioxm. iraystop->1'
                iraystop=1
                return                
             endif
          end if
          goto 111 ! finish/return
        end if ! ib=1
c-------------------------------------------------------------
c end if 2
c
c       ib.gt.1 ions resonance condition may be
c  if 3
        if (ib.gt.1) then
          xb=wpw_2(x,y,z,ib)
          yb=wcw(x,y,z,ib)
          xe=wpw_2(x,y,z,1)
          ye=wcw(x,y,z,1)
          pype=xe/(1.d0+ye)
          pyme=xe/(1.d0-ye)
          pyme2=pype/(1.d0-ye)
          pypb=xb/(1.d0+yb)
          delib=1.d0-yb
          f1b=(s1-pyme2)
          f0b=-pypb

          g1b=cnpar2*(s4-(s1-pyme2))-(s2-pype)*(s3-pyme)-s4*(s1-pyme2)
          g0b=cnpar2*pypb+xb*(s3-pyme)+pypb*s4

          w1b=cnpar2*(-s4*(s1-pyme2)+(s2-pype)*(s3-pyme))+
     1            s4*(s2-pype)*(s3-pyme)
          w0b=cnpar2*(pypb*s4-xb*(s3-pyme))-s4*(s3-pyme)*xb
          fd=f1b*delib+f0b
          gd=g1b*delib+g0b
          wd=w1b*delib+w0b
          detin=gd**2-4.d0*fd*wd

        if (detin.lt.0d0) then ! negative discriminant
             !YuP[11-2016] adjustment for a very small neg.discriminant
             if(abs(detin).lt. 1.d-8*(gd**2+4.d0*fd*wd) ) then
               ! Could be a very small negative value, from rounding.
               ! Then, set it to zero:
               detin=0.d0
             else ! large negative discriminant: no wave here.
               write(*,*)'cninit12_n_gam: detin<0. detin,gd**2,Xe=',
     +                                             detin,gd**2,xe
               write(*,*)'iraystop->1'
               iraystop=1
               return
             endif
        end if
	  cn2p=(-gd+dsqrt(detin))/(2.d0*fd)
	  cn2m=(-gd-dsqrt(detin))/(2.d0*fd)
         
!        if((cn2m.lt.0.d0).and.(cn2p.lt.0.d0)) then
!            write(*,*)'in cninit2 two roots of the dispersion < 0'
!            write(*,*)'cn2m,cn2p',cn2m,cn2p
!            write(*,*)'the given wave can not exist in plasma'
!            write(*,*)'  iraystop->1'
!            iraystop=1
!            return
!        endif

c--------- Check that N^2 is not smaller than tangential-to-surf Ntang^2
          irootp=0 ! initialize
          irootm=0 ! initialize
          if (cn2p.ge.cntang2)then
             dc_l=cnpar/dsqrt(cn2p)
             gam_l=dacos(dc_l)
             write(*,*)'cninit12_n_gam cn2p',cn2p
             irootp=1 ! means: cn2p is ok
          else
             write(*,*)'cninit12_n_gam cn2p<cntang2', cn2p,cntang2
          endif  

          if (cn2m.ge.cntang2)then
             dc_l=cnpar/dsqrt(cn2m)
             gam_l=dacos(dc_l)
             write(*,*)'cninit12_n_gam cn2m',cn2m
             irootm=1 ! means: cn2m is ok
          else
             write(*,*)'cninit12_n_gam cn2m<cntang2', cn2m,cntang2
          endif

          if(irootp+irootm .eq. 0) then ! none of them is good
             write(*,*)'cninit12_n_gam two roots of dispersion <cntang2'
             write(*,*)'cn2m,cn2p',cn2m,cn2p
             write(*,*)'cninit12_n_gam: WAVE CANNOT EXIST.  iraystop->1'
             iraystop=1
             return
          endif  

30        continue ! handle for two trials of roots
          write(*,*)'cninit12_n_gam: iroot=',iroot
          
	    if(irootp+irootm .eq.2)then !both cn2p and cn2m are good (>cntang2)
	      if(iroot.eq.1) then ! first trial
              if(ioxm.eq. 1) cn2=cn2p !first, try this root
              if(ioxm.eq.-1) cn2=cn2m !or this, depending on ioxm
	      else     !   iroot=2: second trial (max two trials)
              if(ioxm.eq. 1) cn2=cn2m !second, try this root
              if(ioxm.eq.-1) cn2=cn2p !or this, depending on ioxm
            endif
	    end if ! irootp+irootm .eq.2
	    
	    if(irootp+irootm .eq.1)then 
	      !only one of cn2p,cn2m was ok (>cntang2)
            if(irootp.eq. 1) cn2=cn2p ! This root was ok
            if(irootm.eq. 1) cn2=cn2m ! This root was ok
            iroot=2 ! to indicate that this is the last trial
	    end if ! irootp+irootm .eq.1


!	  if(iroot.eq.1) then
!           cn2=cn2p
!	     if(cn2.lt.cntang2)then
!	        write(*,*)'in cninit2 cn2p.lt.cntang2'
!	        go to 30
!	     end if
!	  else ! iroot=2
!           cn2=cn2m
!	     if(cn2.lt.cntang2)then
!	        write(*,*)'in cninit2 cn2m.lt.cntang2'
!	        write(*,*)'the given wave can not exist in plasma'
!              write(*,*)'  iraystop->1'
!	        iraystop=1
!	        return
!	     end if
!	  end if

          cnrho2=cn2-cntang2
          cnrho=dsqrt(cnrho2)   
c------------------------------------------------------------------
          cnper_norm=cnrho ! found above
      write(*,'(a,3e11.3)')' before cnxyz:cnpar,cnper_tang,cnper_norm=',
     +cnpar, cnper_tang, cnper_norm
          ! Here - call for ib>1 (electrons and ions):
          !Careful! cnxyz calls rside_xyz -> hamilt_xyz, which uses ioxm
          call cnxyz(x,y,z, cnpar, cnper_tang, cnper_norm,
     +                 cnx,cny,cnz) !->out
      write(*,'(a,3e11.3)')' after  cnxyz: cnx,cny,cnz=',cnx,cny,cnz
c------------------------------------------------------------------YuP
          gam=gamma1_xyz(x,y,z,cnx,cny,cnz)
          ds=dsin(gam)
          dc=dcos(gam)
          ds2=ds*ds
          dc2=dc*dc
          ds4=ds2*ds2
c--------------------------------------------------------------------
c test hamilt
c--------------------------------------------------------------------
c         control that cn2 and gam are the solution of the dispersion
c         relation  n**2=(-b+ioxm*sqrt(b**2-4*a*c))/(2*a)
c---------------------------------------------------------------------
          call abc_xyz(x,y,z,ds2,dc2,ad,bd,cd)
          d4=ad
          d2=bd
          d0=cd
          det=d2*d2-4.d0*d4*d0
          write(*,*)'cninit12_n_gam_xyz: after abc_xyz: det=',det
c---------------------------------------------------
c test hamilt=0?
!          hamtest=fd*cn2*cn2+gd*cn2+wd
!          hamtest=d4*cn2*cn2+d2*cn2+d0
!          pt4=d4*cn2*cn2
!          pt2=d2*cn2
!          pt=pt4+pt2+d0
!          ptm=pt4-pt2+d0
!          gdt=dc2*cn2*((1.d0-xe-xb)
!     1       -(1.d0-xb/(1.d0-yb*yb)-xe/(1.d0-ye*ye)))
!     1       +(-(1.d0-xb/(1.d0-yb)-xe/(1.d0+ye))*
!     1       (1.d0-xb/(1.d0+yb)-xe/(1.d0-ye))-
!     1       (1.d0-xb/(1.d0-yb*yb)-xe/(1.d0-xe*xe))*(1.d0-xe-xb))
!          gdt=gdt*delib
!c         write(*,*)'gd,gdt',gd,gdt
!          btest=-(1.d0-xb/(1.d0-yb*yb)-
!     1         xe/(1.d0-ye*ye))*(1.d0-xe-xb)*
!     1         (1.d0+dc2)-(1.d0-xb/(1.d0-yb)-xe/(1.d0+ye))*
!     1         (1.d0-xb/(1.d0+yb)-xe/(1.d0-ye))*ds2
!c         write(*,*)'bd,btest',bd,btest
!          ptt=(xe*ye*ye*delib/(1.d0-ye*ye)+xb*yb*yb/(1.d0+yb))
!          fmina=-dc2*ptt
!c         write(*,*)'xe,ye,xb,yb',xe,ye,xb,yb
!          ptt1=(-(1.d0-xb/(1.d0-yb)-xe/(1.d0+ye))*
!     1         (1.d0-xb/(1.d0+yb)-xe/(1.d0-ye))+
!     1    (1.d0-xe-xb)*(1.d0-xe/(1.d0-ye*ye)-xb/(1.d0-yb*yb)))*delib
!          gminb=dc2*(-cn2*ptt-ptt1)
!          wminc=-dc2*cn2*ptt1
c-------------------------------------------------------------
          cn1=dsqrt(d2*d2-4.d0*d4*d0)
          write(*,*)'cninit12_n_gam_xyz:(-d2+cn1)/(2*d4)',
     +     (-d2+cn1)/(2.d0*d4)
          write(*,*)'cninit12_n_gam_xyz:(-d2-cn1)/(2*d4)',
     +     (-d2-cn1)/(2.d0*d4)
          cn2new=(-d2+ioxm*cn1)/(2.d0*d4)
          dnp=dabs(cn2-cn2new)/dabs(cn2+cn2new)
          if (dnp.gt.epsmode) then ! Too far from ioxm-target-root
             if(iroot.eq.1) then
                write(*,*)'cninit12_n_gam: dnp>epsmode. Try another.'
                iroot=iroot+1 !->2
                goto 30 ! Try another
             else ! iroot=2: already all trials were made
                write(*,*)'cninit12_n_gam: cn2new far from cn2p or cn2m'
                write(*,*)'No roots N^2(gam) for this ioxm. iraystop->1'
                iraystop=1
                return                
             endif
          end if
          goto 111 ! finish
        end if ! ib>1
c end if 3
      end if ! (id.eq.1).or.(id.eq.2)
c end if 0
 
  111 continue
      write(*,*)'cninit12_n_gam_xyz: Selected cn2new==',cn2new
      return
      end

c======================================================================
c======================================================================


c        **********************cninit12_xyz ***********************
c        *                        -                               *
c        * It solves the dispersion relation N=N(n_par)           *
c        * for cold plasma                                        *
c        * Then subroutine calculates the initial components      *
c        * of the refractive index  cnx,cny,cnz         	          *
c        * If ioxm_n_npar=0 it uses root N(gam,ioxm)              *
c        * If ioxm_n_npar=1 or -1 it uses root N(npar ioxm_n_npar)*
c        *                then calculates angle gam(N_par,N_per)  *
c        *                then finds ioxm which gives the root    *
c        *                N(gam,ioxm)=N(N_par,ioxm_n_npar)        *
c        
c        *********************************************************
c
c------------------------------------------------------------------
c					                          !
c     input parameters		         		          !
c     x,y,z,cnpar,cnper_tang	       	          !
c                                                                 !
c     output parameters				                  !
c     cnx,cny,cnz are the components of the refractive index       !        
c     iraystop=1 end ray calculation                              !
c------------------------------------------------------------------
c     it uses the following functions                             !
c     bxyz,gamma1_xyz                                                    !
c------------------------------------------------------------------
      subroutine cninit12_xyz(x,y,z,cnpar,cnper_tang,
     1                  cnx,cny,cnz,iraystop) !->out
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
      include 'grill.i'
      save
      double complex cmplnper      
      
      write(*,*)'cninit12_xyz -----------> ioxm_n_npar=',ioxm_n_npar
      !pause

      iraystop=0
      pi=4*datan(1.d0)
      cnpar2=cnpar*cnpar
      cntang2=cnper_tang*cnper_tang + cnpar2 !cnteta*cnteta+cnphi*cnphi
      bmod=bxyz(x,y,z) !-> for gamma_xyz?
               
      if(ioxm_n_npar.eq.0) then
c-------------------------------------------------------------
c       calculates the root N(Npar)=N(gam,ioxm) for givem ioxm
c-------------------------------------------------------------
        call cninit12_n_gam_xyz(x,y,z,cnpar,cnper_tang, 
     &                      cnx,cny,cnz,iraystop)
      else
c-------------------------------------------------------------------
c       calculates the root N(N_par,ioxm_n_npar) for given ioxm_n_npar
c       finds ioxm to get  N(N_par,ioxm_n_npar)=N(gam,ioxm)
c-------------------------------------------------------------------
        iroot=0
c------------------------------------------------------------------
c       epsmode is the accuracy  for selection mode (O or X) on
c             the plasma boundary
c------------------------------------------------------------------
        epsmode=1.d-3 !1.d-4 !1.d-8
        id_loc=2
c------------------------------------------------------------------
c       solve the cold plasma id=1,2  
c       dispersion relation cnper=N_per(n_par)
c------------------------------------------------------------------        
        call nper_npar_ioxm_n_npar_xyz(id_loc,x,y,z,cnpar,
     &  cnper,iraystop) !ioxm_n_npar was set in one.i
c-------------------------------------------------------------------
        if(iraystop.eq.1) then
        
          write(*,*)'cninit12_xyz No root Nperp(npar,ioxm_n_npar)'
          return
          
        else ! iraystop=0    Nperp(npar) was found 
        
          if (i_n_poloidal.eq.3) then !input N_parallel, ksi_nperp
            rad_ksi_nperp=ksi_nperp*pi/180.d0 ! degrees to radians
            cnrho=      cnper*dcos(rad_ksi_nperp)
            cnper_tang= cnper*dsin(rad_ksi_nperp)
            cntang2= cnper_tang*cnper_tang + cnpar2
          endif ! i_n_pol=3

          cn2_npar=cnper**2+cnpar**2
          if(cn2_npar.lt.cntang2)then
            write(*,*)'cninit12_xyz (cn2_npar<cntang2)'
            write(*,*)'no root N(npar,ioxm_n_npar) > N_tang'
            write(*,*)'  iraystop->1'
            iraystop=1
            return
          else
            cnrho2=cn2_npar-cntang2
            cnrho=dsqrt(cnrho2)
c--------------------------------------------------------------
c           calculate refractive index components cnx,cny,cnz
c           for given x,y,z,cnpar,cnper_tang
c-------------------------------------------------------------- 
cSAP090601
c           In this case ioxm_n_npar.ne.0 and
c           ioxm is not determined
c           Subroutine cnxyz uses subroutine rside1, 
c           which for id=2 case needs the determined value ioxm 
c           To avoid this problem we will use id=1 in cnxyz
c           May be signs of the radial group velocity
c           for id=1 and id=2 can be different?
c           If yes, then it can create the new problem. 

            if (id.eq.2) then
              id_loc_2=id
              id=1
            endif

            cnper_norm=cnrho ! found above
            call cnxyz(x,y,z, cnpar, cnper_tang, cnper_norm,
     +                 cnx,cny,cnz) !->out
c-----------------------------------------------------------------YuP
            id=id_loc_2 
c--------------------------------------------------------------
c           calculate the angle gam between the refractive index 
c           and the magnetic field
c-------------------------------------------------------------
            gam=gamma1_xyz(x,y,z,cnx,cny,cnz)
c-------------------------------------------------------------
c           calculate cold plasma roots roots 
c           cn_p=N(gam,ioxm=+1), cn_m=N(gam,ioxm=-1)
c-------------------------------------------------------------          
            call n_cold_gam_xyz(x,y,z,gam,cn_p,cn_m,
     &                      iraystop_p,iraystop_m)
c-------------------------------------------------------------
c           choose ioxm for which
c           N(N_par,ioxm_n_npar)=N(gam,ioxm)
c-------------------------------------------------------------
c           check ioxm=+1 root
c-------------------------------------------------------------
            if(iraystop_p.eq.0) then
              delta=dabs(cn2_npar-cn_p**2)/dabs(cn2_npar+cn_p**2)
              if(delta.gt.epsmode)then
c---------------root N(N_par,ioxm_n_npar).ne.N(gam,ioxm=+1)
              else
c---------------root N(N_par,ioxm_n_npar)=N(gam,ioxm=+1)
                ioxm=+1
                iroot=iroot+1
                write(*,*)'cninit12_xyz found ioxm=',ioxm 
c               goto 10 
              endif                
            else
c-------------no root for N(gam,ioxm=+1)
            endif !iraystop_p.eq.0

c-------------------------------------------------------------
c           check ioxm=-1 root
c-------------------------------------------------------------
            if(iraystop_m.eq.0) then
              delta=dabs(cn2_npar-cn_m**2)/dabs(cn2_npar+cn_m**2)
              if(delta.gt.epsmode)then
c---------------root N(N_par,ioxm_n_npar).ne.N(gam,ioxm=-1)
              else
c---------------root N(N_par,ioxm_n_npar)=N(gam,ioxm=-1)
                ioxm=-1                  
                write(*,*)'cninit12_xyz found ioxm=',ioxm 
                iroot=iroot+1
c               goto 10 
              endif                
            else
c-------------no root for N(gam,ioxm=-1)
            endif !iraystop_m.eq.0

c------------------------------------------------------------
            if(iroot.eq.0) then
              write(*,*)'cninit12_xyz no ioxm was found'
              write(*,*)'to get N(N_par,ioxm_n_npar)=N(gam,ioxm)'
              write(*,*)'  iraystop->1'
              iraystop=1
              return 
            else
              if(iroot.eq.2) then
                write(*,*)'*******WARNING******************'
                write(*,*)'cninit12_xyz two ioxm was found'
                write(*,*)'to get N(N_par,ioxm_n_npar)=N(gam,ioxm)'
              endif !iroot.eq.2
            endif !iroot.eq.0
          endif !cn2_npar.lt.cntang2
        endif !iraystop.eq.1
      endif !ioxm_n_npar=0
                                                              
      return
      end

c======================================================================
c======================================================================


      subroutine n_cold_gam_xyz(x,y,z,gam,
     +                          cn_p,cn_m, iraystop_p,iraystop_m)
c----------------------------------------------------------------
c     Calculates cold plasma dispersion relalation roots
c     cn=N(gam,ioxm) in the given space point (x,y,z).
c     as a root of the equation a*N**4+b**2+c=0,
c     cn=(-b+ioxm*sqrt(b**2-4ac))/2a
c     cn_p=(-b+sqrt(b**2-4ac))/2a for ioxm=+1
c     cn_m=(-b-sqrt(b**2-4ac))/2a for ioxm=-1
c
c     Here: 
c     coefficients a=A*delta**3,b=B*delta**3,c=C*delta**3,
c
c     See book  Kroll,Travelspils
c     A=eps_1*sin(gam)**2+eps_3*cos(gam)**2
c     B=-eps_1*eps_2(1+cos(gam)**2)-(eps_1**2-eps_2**2)sin(gam)**2
c     C=eps_3*(eps_1**2-eps_2**2)
c
c     delta=1-y_e for electrons ib=1 (ib is used in subroutine
c                                     abc) 
c          =1-y_i for ions      ib=i>1
c----------------------------------------------------------------     
      implicit none
c-----input
      real*8 x,y,z,   ! space coordinates of the given point
     *gam             ! the angle [radians] between the wave vector
                      ! and the magnetic field                      
c-----output
      real*8 cn_p,cn_m    !roots cn_p=N(gam,ioxm=+1), cn_m=N(gam,ioxm=-1),  
      integer iraystop_p, !=0 the root with ioxm=+1 found
                          !=1 no root  with ioxm=+1 did not find
     &        iraystop_m  !=0 the root with ioxm=-1 found
                          !=1 no root  with ioxm=-1 did not find
 
c-----locals
      real*8 ds,dc,ds2,dc2,
     &a,b,c,det,cn2_p,cn2_m

      iraystop_p=0
      iraystop_m=0
      
      ds=dsin(gam)
      dc=dcos(gam)
      ds2=ds*ds
      dc2=dc*dc

      call abc_xyz(x,y,z,ds2,dc2,a,b,c)

      det=b**2-4.d0*a*c
      if(det.lt.0.d0)then
        write(*,*)'in subr. n_cold_gam_xyz det<0'
        write(*,*)'the cold plasma dispersion eq.'
        write(*,*)'has not root N_perp(gam)'
        write(*,*)'for the given angle gam=',gam
        iraystop_p=1
        iraystop_m=1
        return
      endif

      cn2_p=(-b+dsqrt(det))/(2.d0*a) !N**2
      if(cn2_p.gt.0.d0)then
        cn_p=dsqrt(cn2_p) 
      else
        write(*,*)'in subr. n_gam negative N**2=cn2_p<0'
        write(*,*)'The cold plasma dispersion eq.'
        write(*,*)'has not positive root N_perp(gam,ioxm)'
        write(*,*)'for the given angle gam=',gam,' and ioxm=+1'
        iraystop_p=1
      endif

      cn2_m=(-b-dsqrt(det))/(2.d0*a) !N**2
      if(cn2_m.gt.0.d0)then
        cn_m=dsqrt(cn2_m) 
      else
        write(*,*)'in subr. n_cold_gam_xyz negative N**2=cn2_m<0'
        write(*,*)'The cold plasma dispersion eq.'
        write(*,*)'has not positive root N_perp(gam,ioxm)'
        write(*,*)'for the given angle gam=',gam,' and ioxm=-1'
        iraystop_p=1
      endif

      return
      end

      
c======================================================================
c======================================================================


c        **********************nper_npar_ioxm_n_npar_xyz ******
c                                -                            
c         It solves the cold plasma id=1,2 or
c         Appleton-Hartree id=3 dispersion relation N_per=N_per(n_par)
c         for given ioxm_n_npar                              
c        ******************************************************
c
c---------------------------------------------------
c     							   !
c     input parameters					   !
c     x,y,z,cnpar 	                                   !
c     ioxm_n_npar is inside one.i common block             !
c     It works for cold plasma dispersion functions         !
c     id_loc=1,2,3                                         !
c     output parameters: N_per(N_par)=cnper                 !
c     iraystop=0 the root was found                        !
c     iraystop=1 the root was not found                        !
c----------------------------------------------------------
c     it uses the following functions and subroutines  
c     ias1r,bxyz,wpw_2,wcw,gamma1_xyz,s_xyz,abc_xyz         
c-----------------------------------------------------------

      subroutine nper_npar_ioxm_n_npar_xyz(id_loc,x,y,z,cnpar,
     &cnper,iraystop)

      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'

      iraystop=0 ! remains 0 if the root (cnper) is found
      
      write(*,*)
      write(*,*)'nper_npar_ioxm_n_npar_xyz: Nper(Npar) for cold plasma'

      pi=4*datan(1.d0)
      cnpar2=cnpar*cnpar
      cnpar4=cnpar2*cnpar2
      cnper=0.d0 ! to initialize ! Should be found below.
      bmod=bxyz(x,y,z) !-> get bmod for wcw
c-------------------------------------------------------------
c     calculations of  cnrer2p and cnper2m
c     from cnpar2 by using the dispersin relation
c     f*cnper**4 +(2f*cnpar**2+g)*cnper**2+(f*cnpar**4+g*cnpar**2+w)=0
c     f*cnper**4 +gnew*cnper**2+wnew=0
c     gnew=2fcnpar**2+g, wnew=f*cnpar**4+g*cnpar**2+w
c     in the following form cnper**2=(-gnew+,-*sqrt(gnew**2-4*f*wnew))/(2*f)
c------------------------------------------------------------------
c     Appleton - Hartree dispersion relation
      if (id_loc.eq.3) then
         xi=wpw_2(x,y,z,1)
         yi=wcw(x,y,z,1)
         py2=yi*yi
         py4=py2*py2
         px=1.d0-xi
         px2=px*px
c------------------------------------------------------------------
         pyp=xi/(1.d0+yi)
         f1e=1.d0
         f0e=-pyp
         g1e=cnpar2*(-xi)+pyp+xi-2.d0
         g0e=cnpar2*pyp+xi*(1.d0-pyp)+pyp*(1.d0-xi)
         w1e=cnpar2*(xi-pyp)+(1.d0-pyp)*(1.d0-xi)
         w0e=cnpar2*(pyp*(1.d0-xi)-xi*(1.d0-pyp))-
     &       (1.d0-xi)*(1.d0-pyp)*xi
         dele=1.-yi
         fd=f1e*dele+f0e
         gd=g1e*dele+g0e
         wd=w1e*dele+w0e
c new coefficients
         gnew=gd+2.d0*fd*cnpar2
         wnew=wd+gd*cnpar2+fd*cnpar4
         gd=gnew
         wd=wnew

         detin=gd**2-4.d0*fd*wd
         if(detin.lt.0.d0) then ! negative discriminant
             !YuP[11-2016] adjustment for a very small neg.discriminant
             if(abs(detin).lt. 1.d-8*(gd**2+4.d0*fd*wd) ) then
               ! Could be a very small negative value, from rounding.
               ! Then, set it to zero:
               detin=0.d0
             else ! large negative discriminant: no wave here.
               write(*,*)'nper_npar_ioxm_n_npar_xyz: detin,gd**2,Xe=',
     +                                           detin,gd**2,xe
               write(*,*)'nper_npar_ioxm_n_npar_xyz detin<0'
               write(*,*)'nper_npar_ioxm_n_npar_xyz no roots'
               write(*,*)'  iraystop->1'
               iraystop=1
               return
           endif
         endif

         cnper2=(-gd+ioxm_n_npar*dsqrt(detin))/(2.d0*fd)
         if(cnper2.lt.0d0) then
            write(*,*)'nper_npar_ioxm_n_npar_xyz cnper2<0'
            write(*,*)'the root N_perp**2 is negative'
            write(*,*)'nper_npar_ioxm_n_npar_xyz no positive root'    
            write(*,*)'  iraystop->1'
            iraystop=1
            return
         else
            cnper=dsqrt(cnper2)
         endif           
      end if !id_loc.eq.3
c------------------------------------------------------------------
c     cold plasma dispresion relation
c------------------------------------------------------------------
      if ((id_loc.eq.1).or.(id_loc.eq.2)) then
      
        call s_xyz(x,y,z,s1,s2,s3,s4,s6,s7)

        if (ib.eq.1) then
c---------ib=1 electron resonance condition may be present in plasma
          xe=wpw_2(x,y,z,1)
          ye=wcw(x,y,z,1)
          pyp=xe/(1.d0+ye)
          dele=1.d0-ye
          f1e=s7
          f0e=-pyp
          g1e=cnpar2*(s4-s7)-(s6-pyp)*s3-s4*s7
          g0e=cnpar2*pyp+xe*(s6-pyp)+pyp*s4
          w1e=cnpar2*(-s7*s4+(s6-pyp)*s3)+(s6-pyp)*s3*s4
          w0e=cnpar2*(pyp*s4-xe*(s6-pyp))-s4*(s6-pyp)*xe
          fd=f1e*dele+f0e
          gd=g1e*dele+g0e
          wd=w1e*dele+w0e
          ! new coefficients
          gnew=gd+2.d0*fd*cnpar2
          wnew=wd+gd*cnpar2+fd*cnpar4
          gd=gnew
          wd=wnew
          write(*,*)'nper_npar_ioxm_n_npar_xyz ioxm_n_npar=',ioxm_n_npar
          detin=gd**2-4.d0*fd*wd
          if(detin.lt.0.d0) then ! negative discriminant
             !YuP[11-2016] adjustment for a very small neg.discriminant	  
             if(abs(detin).lt. 1.d-8*(gd**2+4.d0*fd*wd) ) then
               ! Could be a very small negative value, from rounding.
               ! Then, set it to zero:
               detin=0.d0
             else ! large negative discriminant: no wave here.
               write(*,*)'nper_npar_ioxm_n_npar_xyz: detin,gd**2,Xe=',
     +                                           detin,gd**2,xe           
               write(*,*)'nper_npar_ioxm_n_npar_xyz detin<0'
               write(*,*)'nper_npar_ioxm_n_npar_xyz No roots'
               write(*,*)'  iraystop->1'
               iraystop=1
               return
            endif
          end if

          cnper2=(-gd+ioxm_n_npar*dsqrt(detin))/(2.d0*fd)
          if(cnper2.lt.0.d0)then
            write(*,*)'  nper_npar_ioxm_n_npar_xyz:  N_perp**2 <0'
            write(*,*)'  iraystop->1'
            iraystop=1
            return
          else
            cnper=dsqrt(cnper2)
            write(*,*)'  nper_npar_ioxm_n_npar_xyz:  N_perp=',cnper    
          endif  
                   
          goto 111 !-> finish/exit
        end if !ib.eq.1
  

        if(ib.gt.1) then
c---------ib.gt.1 both electrons ions resonance condition may be present
          xb=wpw_2(x,y,z,ib)
          yb=wcw(x,y,z,ib)
          xe=wpw_2(x,y,z,1)
          ye=wcw(x,y,z,1)
          pype=xe/(1.d0+ye)
          pyme=xe/(1.d0-ye)
          pyme2=pype/(1.d0-ye)
          pypb=xb/(1.d0+yb)
          delib=1.d0-yb
          f1b=(s1-pyme2)
          f0b=-pypb

          g1b=cnpar2*(s4-(s1-pyme2))-(s2-pype)*(s3-pyme)-s4*(s1-pyme2)
          g0b=cnpar2*pypb+xb*(s3-pyme)+pypb*s4

          w1b=cnpar2*(-s4*(s1-pyme2)+(s2-pype)*(s3-pyme))+
     1        s4*(s2-pype)*(s3-pyme)
          w0b=cnpar2*(pypb*s4-xb*(s3-pyme))-s4*(s3-pyme)*xb
          fd=f1b*delib+f0b
          gd=g1b*delib+g0b
          wd=w1b*delib+w0b
c new coefficients
          gnew=gd+2.d0*fd*cnpar2
          wnew=wd+gd*cnpar2+fd*cnpar4
          gd=gnew
          wd=wnew

          detin=gd**2-4.d0*fd*wd
          write(*,*)'ib=',ib, '   Trying ioxm_n_npar=',ioxm_n_npar
          if(detin.lt.0.d0) then ! negative discriminant
             !YuP[11-2016] adjustment for a very small neg.discriminant	  
             if(abs(detin).lt. 1.d-8*(gd**2+4.d0*fd*wd) ) then
               ! Could be a very small negative value, from rounding.
               ! Then, set it to zero:
               detin=0.d0
             else ! large negative discriminant: no wave here.
               write(*,*)'nper_npar_ioxm_n_npar_xyz: detin,gd**2,Xe=',
     +                                           detin,gd**2,xe        
               write(*,*)'nper_npar_ioxm_n_npar_xyz detin<0'
               write(*,*)'nper_npar_ioxm_n_npar_xyz no roots'
               write(*,*)'  iraystop->1'
               iraystop=1
               return
            endif
          end if

          cnper2=(-gd+ioxm_n_npar*dsqrt(detin))/(2.d0*fd)
          if(cnper2.lt.0.d0) then
            write(*,*)'  cnper2<0.  The root N_perp**2 is negative'
                   write(*,*)'  iraystop->1'
            iraystop=1
            return
          else
            cnper=dsqrt(cnper2)
            write(*,*)'  nper_npar_ioxm_n_npar_xyz:  N_perp=',cnper    
          endif   
                  
          goto 111
        end if !ib>1
      end if !(id_loc.eq.1).or.(id_loc.eq.2)
 
  111 continue

      return
      end

c======================================================================
c======================================================================


c        **********************npernpar_xyz *******************
c        *                        -                           *
c        * It solves the dispersion relation Nper=N(n_par)    *
c        ******************************************************
c
c---------------------------------------------------
c     input parameters					   !
c     x,y,z,cnpar 	                                   !
c     output parameters: cnper2p,cnper2m (two roots)              
c---------------------------------------------------
c     it uses the following functions and subroutines      !
c     bxyz, wpw_2, wcw, s_xyz
c-----------------------------------------------------------
      subroutine npernpar_xyz(x,y,z,cnpar, cnper2p,cnper2m)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'

      pi=4*datan(1.d0)
      cnpar2= cnpar*cnpar
      cnpar4= cnpar2*cnpar2
      cnper2p=0.d0
      cnper2m=0.d0
      bmod= bxyz(x,y,z) ! get bmod; needed for wcw
c      write(*,*)'npernpar id,ib= ',id,ib

c-------------------------------------------------------------
c     calculations of  cnper2p and cnper2m
c     from cnpar2 by using the dispersion relation
c     f*cnper**4 +(2f*cnpar**2+g)*cnper**2+(f*cnpar**4+g*cnpar**2+w)=0
c     f*cnper**4 +gnew*cnper**2+wnew=0
c     gnew=2fcnpar**2+g, wnew=f*cnpar**4+g*cnpar**2+w
c     in the following form cnper**2=(-gnew+,-*sqrt(gnew**2-4*f*wnew))/(2*f)
c------------------------------------------------------------------
c     Appleton - Hartree dispersion relation
c
c if 1
      if (id.eq.3) then
         xi=wpw_2(x,y,z,1)
         yi=wcw(x,y,z,1)
         py2=yi*yi
         py4=py2*py2
         px=1.-xi
         px2=px*px
c------------------------------------------------------------------
         pyp=xi/(1.+yi)
         f1e=1.d0
         f0e=-pyp
         g1e=cnpar2*(-xi)+pyp+xi-2.
         g0e=cnpar2*pyp+xi*(1.-pyp)+pyp*(1.-xi)
         w1e=cnpar2*(xi-pyp)+(1.-pyp)*(1.-xi)
         w0e=cnpar2*(pyp*(1-xi)-xi*(1.-pyp))-(1.-xi)*(1-pyp)*xi
         dele=1.-yi
         fd=f1e*dele+f0e
         gd=g1e*dele+g0e
         wd=w1e*dele+w0e
c new coefficients
         gnew=gd+2.d0*fd*cnpar2
         wnew=wd+gd*cnpar2+fd*cnpar4
         gd=gnew
         wd=wnew

         detin=gd**2-4.*fd*wd
         if (detin.lt.0d0) then ! negative discriminant
             !YuP[11-2016] adjustment for a very small neg.discriminant	 
             if(abs(detin).lt. 1.d-8*(gd**2+4.d0*fd*wd) ) then
               ! Could be a very small negative value, from rounding.
               ! Then, set it to zero:
               detin=0.d0
             else ! large negative discriminant: no wave here.
               write(*,*)'npernpar_xyz: detin,gd**2,Xe=',
     +                                  detin,gd**2,xe
               write(*,*)' 1 in npernpar_xyz detin<0'
               return
            endif
         endif
         cnper2p=(-gd+dsqrt(detin))/(2.*fd)
         cnper2m=(-gd-dsqrt(detin))/(2.*fd)
      end if
c end if 1
c------------------------------------------------------------------
c     cold plasma dispresion relation
c------------------------------------------------------------------
c if 0
      if ((id.eq.1).or.(id.eq.2)) then
        call s_xyz(x,y,z,s1,s2,s3,s4,s6,s7)

c       ib=1 electron resonance condition may be present
c  if 2
        if (ib.eq.1) then
          xe=wpw_2(x,y,z,1)
          ye=wcw(x,y,z,1)
          pyp=xe/(1.d0+ye)
          dele=1.d0-ye
          f1e=s7
          f0e=-pyp
          g1e=cnpar2*(s4-s7)-(s6-pyp)*s3-s4*s7
          g0e=cnpar2*pyp+xe*(s6-pyp)+pyp*s4
          w1e=cnpar2*(-s7*s4+(s6-pyp)*s3)+(s6-pyp)*s3*s4
          w0e=cnpar2*(pyp*s4-xe*(s6-pyp))-s4*(s6-pyp)*xe
          fd=f1e*dele+f0e
          gd=g1e*dele+g0e
          wd=w1e*dele+w0e
c new coefficients
          gnew=gd+2.d0*fd*cnpar2
          wnew=wd+gd*cnpar2+fd*cnpar4
          gd=gnew
          wd=wnew

          detin=gd**2-4.d0*fd*wd
          if (detin.lt.0d0) then ! negative discriminant
             !YuP[11-2016] adjustment for a very small neg.discriminant
             if(abs(detin).lt. 1.d-8*(gd**2+4.d0*fd*wd) ) then
               ! Could be a very small negative value, from rounding.
               ! Then, set it to zero:
               detin=0.d0
             else ! large negative discriminant: no wave here.
               write(*,*)'npernpar_xyz: detin<0. detin,gd**2,Xe=',
     +                                           detin,gd**2,xe
               !pause
               cnper2p=-1.d0 ! just to indicate: non-propagating
               cnper2m=-1.d0 ! just to indicate: non-propagating
               return
             endif
          end if
          cnper2p=(-gd+dsqrt(detin))/(2.d0*fd)
          cnper2m=(-gd-dsqrt(detin))/(2.d0*fd)
          goto 111
        end if
c end if 2
c
c     ib.gt.1 -> ion resonance condition may be present
c  if 3
        if (ib.gt.1) then
c          write(*,*)'cold plasma ib .gt.1  '
           xb=wpw_2(x,y,z,ib)
           yb=wcw(x,y,z,ib)
           xe=wpw_2(x,y,z,1)
           ye=wcw(x,y,z,1)
           pype=xe/(1.d0+ye)
           pyme=xe/(1.d0-ye)
           pyme2=pype/(1.d0-ye)
           pypb=xb/(1.d0+yb)
           delib=1.d0-yb
           f1b=(s1-pyme2)
           f0b=-pypb

           g1b=cnpar2*(s4-(s1-pyme2))-(s2-pype)*(s3-pyme)-s4*(s1-pyme2)
           g0b=cnpar2*pypb+xb*(s3-pyme)+pypb*s4

           w1b=cnpar2*(-s4*(s1-pyme2)+(s2-pype)*(s3-pyme))+
     1	     s4*(s2-pype)*(s3-pyme)
           w0b=cnpar2*(pypb*s4-xb*(s3-pyme))-s4*(s3-pyme)*xb
           fd=f1b*delib+f0b
           gd=g1b*delib+g0b
           wd=w1b*delib+w0b
c new coefficients
           gnew=gd+2.d0*fd*cnpar2
           wnew=wd+gd*cnpar2+fd*cnpar4
           gd=gnew
           wd=wnew

           detin=gd**2-4.d0*fd*wd
           if (detin.lt.0d0) then ! negative discriminant
             !YuP[11-2016] adjustment for a very small neg.discriminant
             if(abs(detin).lt. 1.d-8*(gd**2+4.d0*fd*wd) ) then
               ! Could be a very small negative value, from rounding.
               ! Then, set it to zero:
               detin=0.d0
             else ! large negative discriminant: no wave here.
               write(*,*)'npernpar_xyz: detin,gd**2,Xe=',
     +                                  detin,gd**2,xe
               write(*,*)'npernpar_xyz detin<0'
               return
             endif
           end if
           cnper2p=(-gd+dsqrt(detin))/(2.d0*fd)
           cnper2m=(-gd-dsqrt(detin))/(2.d0*fd)
c        write(*,*)'in npernpar ib.qt.1 cnper2p,cnper2m',cnper2p,cnper2m
           goto 111
        end if ! ib>1
c end if 3
      end if ! npernpar_xyz/END: (id.eq.1).or.(id.eq.2)
c end if 0
  111 continue

      return
      end ! npernpar_xyz()/END

c======================================================================
c======================================================================

      subroutine correct2_xyz(cnpar,cnpero,cnpern,
     * cnxo,cnyo,cnzo, bx,by,bz, 
     + cnxn,cnyn,cnzn) ! out
c-----------------------------------------------------------
c     it calculates new values of cnxn,cnyn,cnzn,
c     based on new N_perp(N_parallel) found by solvnperp() 
c     such that the Hamiltonian (disperion determinant) is still zero.
c     Note: Nper was adjusted, but Npar is the same.
c-----------------------------------------------------------
      implicit none
c-----input
      double precision 
     *cnpar,cnpero,cnpern, ! Npar, old Nper and new Nper(Npar)
     *cnxo,cnyo,cnzo, ! old values of N-vector components
     *bx,by,bz ! magnetic field components
c-----output
      double precision cnxn,cnyn,cnzn !New values of N-vector components

c-----local
      double precision cnpar2,cn2, bmod, a,b,c,d,det,
     *cnxm,cnym,cnzm, cnxp,cnyp,cnzp, delnm,delnp
      integer inewn
c------------------------------------------------------------  
c     inewn=0 ! Could not find new (cnxn,cnyn,cnzn);
c               continue using old (cnxo,cnyo,cnzo)
c          =1 ! New values (cnxn,cnyn,cnzn) are found.
c------------------------------------------------------------  
c     Find components (cnxn,cnyn,cnzn) of Nnew from three Eqs:
c     1. (Nnew.B) = Npar|B|   [(.) is the scalar product of vectors]
c        which is written as
c        cnxn*bx + cnyn*by + cnzn*bz = cnpar*bmod
c     2. cnxn^2 + cnyn^2 + cnzn^2 = cnper^2 + cnpar^2
c     3. Assume that the new Npern is in the same direction as old Npero,
c        (Npern.Npero) = |Npern|*|Npero|  
c        (angle between Npern and Npero is 0)
c        which can be written as
c        cnxn*cnxo + cnyn*cnyo + cnzn*cnzo = cnpero*cnpern + cnpar^2
c-----------------------------------------------------------------
      cnpar2= cnpar*cnpar   
      cn2 = cnpar2 + cnpern*cnpern  ! N^2 new 
      bmod= dsqrt(bx*bx + by*by + bz*bz)
      
c      a= (cnxo*bz-cnzo*bx)/(cnzo*by-cnyo*bz)
c      c= (cnyo*bx-cnxo*by)/(cnzo*by-cnyo*bz)
c      b= ( cnzo*cnpar*bmod-(cnpern*cnpero+cnpar2)*bz)/(cnzo*by-cnyo*bz)
c      d= (-cnyo*cnpar*bmod+(cnpern*cnpero+cnpar2)*by)/(cnzo*by-cnyo*bz)
      
      inewn=0
      det= -1.d0 ! For now; finish later
      
      if (det.lt.0.d0) then
         write(*,*)'correct2_xyz det<0; 
     +    new (cnxn,cnyn,cnzn) cannot be found;
     *    correction was not made'       
         goto 10
      else
        ! needs work
      endif

      delnm= (cnxm-cnxo)**2 + (cnym-cnyo)**2 + (cnzm-cnzo)**2
      delnp= (cnxp-cnxo)**2 + (cnyp-cnyo)**2 + (cnzp-cnzo)**2
      if (delnm.le.delnp) then
        cnxn=cnxm
        cnyn=cnym
        cnzn=cnzm
      else
        cnxn=cnxp
        cnyn=cnyp
        cnzn=cnzp
      endif

 10   continue

      if (inewn.eq.0) then ! could not find new values.
c-------the correction was not made; use old values:
        cnxn=cnxo
        cnyn=cnyo
        cnzn=cnzo
      endif

      return
      end

c======================================================================
c======================================================================

c     *****  absorp_relativist_disp_combined_xyz  ***********
c 
c     *  this subroutine calculates imaginary part       
c     *  of refractive index using formula from:   	 
c     *  2k^_i_s*(P^+T^)=P_abs_s {Stix p.74 (17)}        
c     *  and relativistic tensor from Disp_combined         
c     *   05/03/18                                       
c     ****************************************************
c         input parameters: x,y,z  (m) 
c                           cnpar -N_parallel (to magnetic field) 
c                           cnper -N_perpendicular (real part)       
c
c         output parameter: cnprim_e -imaginary part of N_perp    
c                                     (electron absorption)	  
c
c                           It calculates the electric field polarization
c                           using the hot plasma tensor.The polarization
c                           (cex,cey,cez) will be in common /efield.i/ 
c------------------------------------------------------------------
      subroutine absorp_relativist_disp_combined_xyz(x,y,z,cnpar,cnper,
     &cnprim_e)
c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none
      include 'param.i'
      include 'eps.i'
      include 'cefield.i'
      
c-----input 
      real*8  x,y,z,cnper,cnpar
c-----output
      real*8 cnprim_e
    
c-----local
      real*8 X_e,Y_e,T_e,
     & cnprim,power_abs_e,omega,bmod,pi,
     & T_perp_r,d_reps_herm_d_n_perp,
     & ex,ey,ez,eplus,eminus,epar,pp,step,
     & Poynting_vector_perp_r

      complex*16 D_rel,! relativistic dispersion function Determinant
     & K_s(3,3),      ! susceptibilties of the given specie  
     & cnprp,           ! complex n_perp
     & ci,
     & kappa_e(3,3), !anti-Hermitian part of eps for electron specie
     & T_perp_c,power_abs_e_c,
     & reps_p_c(3,3),reps_m_c(3,3),d_reps_herm_d_n_perp_c,
     & cbx,cby,cbz,ce(3),cec(3),cb(3),cbc(3),
     & Poynting_vector_x_c, c_cnper,cnper_p,cnper_m
C ENM 26jul06 -- added c_cnper to be a complex version of cnper (with 0 imaginary part)
C  because disp_combined expects a complex n_perp as input (Also moved cnper_p and cnper_m to the complex*16 list
      integer iherm_loc,i,j,j1,j2
 
c-----external
      real*8 bxyz,wpw_2,wcw,tempe_xyz,tpoprho

      pi=4.d0*dtan(1.d0)
      c_cnper=dcmplx(cnper,0.0d0)
c----------------------------------------------
c     calculate relativistic plasma tensor reps(). 
c     reps() will be in common /eps/ in eps.i file.
c---------------------------------------------
      bmod=bxyz(x,y,z) ! get bmod; needed for wcw
     
      X_e=wpw_2(x,y,z,1) ! rho is calc.here
      Y_e=wcw(x,y,z,1)
c     if(i.eq.1) y_ar(1)=-y_ar(1) ? question
      T_e=tempe_xyz(x,y,z,1)        ! kev
             
      call Disp_combined(T_e,cnpar,X_e,Y_e,c_cnper,reps,D_rel)
    
c--------------------------------------------------------------
c     Calculate electric field polarization cex,cey,cez
c     using relativistic dielectric tensor reps()
c     and real N_perp=cnper.
c     The polarization will be in common /cefield/ in cefield.i file.
c-------------------------------------------------------------------   
      cnprp=dcmplx(cnper,0.d0)
      call efield1(cnpar,cnprp,ex,ey,ez,eplus,eminus,epar)
      ce(1)=cex
      ce(2)=cey   
      ce(3)=cez
c-------------------------------------------------------------------
c     complex components of the wave magnetic field
c
c     ce - wave electric field
c     cec- comlex conjugate wave electric field
c     cb - wave magnetic field
c     cbc- comlex conjugate wave magnetic field
c-------------------------------------------------------------------
      cbx=-cnpar*cey
      cby=-cnper*cez+cnpar*cex
      cbz=cnper*cey
 
      cb(1)=cbx
      cb(2)=cby
      cb(3)=cbz
      
      do j=1,3
         cec(j)=dconjg(ce(j))
         cbc(j)=dconjg(cb(j))
      enddo      
c------------------------------------------------------------------
c     Calculate apsorped power for all species
c    
c     Power=omega/(8Pi)[E~(i) . eps_a_herm(i,j) 'E~(j)] 
c------------------------------------------------------------------
c     The code calculates absorped powes for electron specie
c     using anti-hermitian part of the susceptibilty
c     The code uses normalized Power=Power/omega
c------------------------------------------------------------------
c     1. Calculate susceptibilities K_s() for all species.      
c     2. Calculate absorped power  power_abs_s(i) for
c        all "i=1,...,nbulk" species.      
c-----------------------------------------------------------------
      ci=dcmplx(0.d0,1.d0)

      cnprim=0.d0
c     omega=2.d0*pi*frqncy*1.d9       !frqncy in GHZ

      power_abs_e_c=0.d0

      do j1=1,3
         do j2=1,3
            kappa_e(j1,j2)=0.5d0*(reps(j1,j2)-dconjg(reps(j2,j1)))/ci
            power_abs_e_c=power_abs_e_c+cec(j1)*kappa_e(j1,j2)*ce(j2)
         enddo
      enddo

      power_abs_e_c=power_abs_e_c/(8.d0*pi)
      power_abs_e=dreal(power_abs_e_c)

      write(*,*)'absorp__relativist_disp_combined'
      write(*,*)'power_abs_e_c,power_abs_e',
     &          power_abs_e_c,power_abs_e
     
c------------------------------------------------------------------
c     calculations of the Poynting flux 
c
c                      c_light
c     Poynting_flux^=  ------- [ E^~ (vectr*) B^+ E^ (vectr*) B^~ ]
c                       16 Pi 
c-------------------------------------------------------------------
c     x direction of the Stix sytem is parallel to N_perp^.
c     So, x component of the Poynting flux is parallel to N_perp^.
c
c     Poynting_vector_x_c is complex x component of 
c     Poynting flux in the Stix coordinates. 
c------------------------------------------------------------
c     In the code it was used the normalized variable 
c     Poynting_vector_perp_r= Poynting /c_light
c-------------------------------------------------------------
      pp=1.d0/(16.d0*pi)
      Poynting_vector_x_c = pp*((cec(2)*cb(3)-cec(3)*cb(2))+
     &                          (ce(2)*cbc(3)-ce(3)*cbc(2)))     

      Poynting_vector_perp_r=dreal(Poynting_vector_x_c)
c-------------------------------------------------------------------
c     calculations T=T_vector the flux of nonelectromagnetic energy
c     Stix book p.74 (19)
c           omega        d                c_light       d
c     T^= - ----- E^~.----(eps_h^^).E^= - -------- E^~. ---(eps_h^^) .E^
c           16Pi         dk^              16Pi          dN^
c
c     For calculations of Im (K_perp) we need T^_perp.
c     T^perp it is parallel to Re( k^_perp)
c
c                c_light          d
c     T^_perp= - -------- E^~(i). -----------(eps_h(i,j) . E(j)
c                 16Pi            dRe(N_perp)              
c   
c-------------------------------------------------------------------
c     In the code it was used the normalized variable
c     T_perp_r = (T/c_light) 
c------------------------------------------------------------------
      step=1.d-5 !1.d-3 absorp_relativist_disp_combined_xyz

      cnper_p=c_cnper*(1.d0+step)
      cnper_m=c_cnper*(1.d0-step)
         
      call Disp_combined(T_e,cnpar,X_e,Y_e,cnper_p,reps_p_c,D_rel)

      call Disp_combined(T_e,cnpar,X_e,Y_e,cnper_m,reps_m_c,D_rel)
      
      T_perp_c=dcmplx(0.d0,0.d0)

      do j1=1,3
        do j2=1,3
          d_reps_herm_d_n_perp_c=(reps_p_c(j1,j2)-
     &                            reps_m_c(j1,j2))/(2.d0*step*cnper)  

          T_perp_c=T_perp_c+cec(j1)*d_reps_herm_d_n_perp_c*ce(j2)
          enddo
      enddo

      T_perp_r=-dreal(T_perp_c)/(16.d0*pi)
c      write(*,*)'T_perp_c,T_perp_r',T_perp_c,T_perp_r 
      
      cnprim_e=0.5d0*power_abs_e/(T_perp_r+Poynting_vector_perp_r)

      cnprim_e=dabs(cnprim_e)

      if(cnprim_e.lt.0.d0) then
        write(*,*)'!!!cnprim_e<0 cnprim_e', cnprim_e 
      endif
     
      write(*,*)'T_perp_r,Poynting_vector_perp_r,T/P',
     &T_perp_r,Poynting_vector_perp_r,T_perp_r/Poynting_vector_perp_r

      write(*,*)'cnprim_e',cnprim_e
         
c--------------------------------------------------------
c     calculated comples relativistic dielectric tensor reps
c     reps will be in common /eps.i/
c--------------------------------------------------------     
      call Disp_combined(T_e,cnpar,X_e,Y_e,c_cnper,reps,D_rel)
    

      return
      end

c======================================================================
c======================================================================


c     ********************** absorpfd_xyz ****************
c     *                      -----                       *
c     *  this subroutine calculates imaginary part       *
c     *  of refractive index using formula from:   	 *
c     *  S.C.Chiu,V.S.Chan,R.W.Harvey,M.Porkolab	 *
c     *  Theory of fast wave current drive for tokamak   *
c     *  plasma,					 *
c     *  Nuclear fusion,1989,vol.29,No.12, p.2175-2186   *
c     !!the differens with absorpwf in the modified bessel function 
c     !!the real*8 bessel function it gives the good ions absorption
c     *                                                  *
c     *  Smirnov021018: Two corrections are introduced:  *
c     *    Eq. (6) for eps23 has "-" sign                *
c     *    Substitute Z_0(zeta)|Stix for all Z(zeta_e)   *
c     ****************************************************
c         input parameters: x,y,z  (m) 
c                           cnpar -N_parallel (to magnetic field) 
c                           cnper -N_perpendicular (real part)       
c                           tempe- electron temperature  (keV)
c                           dense- electron density (10*13/cm**3
c                           nbulk is the number of the plasma species
c                           bmod is the total magnetic field (Tl)
c                           tempiar(2:nbulka) ions temperature(kev)
c                           nbulk is the number of the plasma species
c                           bmod is the total magnetic field (Tl)

c         The input parameter nb (the max number of the modified Bessel
c         functions) is set inside absorpfw using operator 'parameter'. 
      
c         output parameter: cnprim_e -imaginary part of N_perp    
c                                     (electron absorption)	  
c                           cnprim_i -imaginary part of N_perp    
c                                     (ion absorption) 
c
cBH041009:  Added breakdown of ion damping into components, cnprim_s(2:nbulk)
c 	  
c-----------------------------------------------------------------*
c    it uses beslci1 modified double precision Bessel functions   *
c------------------------------------------------------------------
      subroutine absorpfd_xyz(x,y,z,f000,rho_larm_max0,cnpar,cnper,
     1 tempe,dense,tempiar,
     1 nbulk,bmod,cnprim_e,cnprim_i,cnprim_s)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'ions.i'
      include 'eps.i'

      dimension cnprim_s(*)   !BH041009, for separate ion contributions
                              !cnprim_s(1:nbulk), cnprim_s(1)=cnprim_e
      dimension tempiar(nbulka)
      double precision mu,mu2
      parameter (nb=500) ! the max order of the modified Bessel function
      dimension br(nb),bi(nb)

      double complex argz,zf,zfp
      double precision rpr,rpim

      double complex ci,ckappae2,czero
         
      czero=dcmplx(0.d0,0.d0)
c********************************************************************
c electron absorption:
c     Imag(N_perp)=abs(N_perp)*sqrt(pi)/4*beta_e*ksi_e*exp(-ksi_e**2)*G
c--------------------------------------------------------------------
      pi=4.d0*datan(1.d0)
      spi=dsqrt(pi)
c----------------------------------------
c     cvac - light velocity in (cm/sec)
      cvac=2.997930d+10
c     ve=(sqrt(2Te/me)),    ve in (cm/sec),Te(keV)
      ve=1.87d9*dsqrt(tempe)	    !tempe keV
      omega=2.d0*pi*f000*1.d9 ![1/sec]      
c----------------------------------------
c     cksi_e=omega/(k_par*ve)=cvac/(N_par*ve)
      if(cnpar.ne.0.d0) then
         cksi_e=cvac/(cnpar*ve)
      else
         cksi_e=0.d0
      endif
cSm021022
c      cksi_e=dabs(cksi_e)
c----------------------------------------
c     beta_e=p_e/(btot**2/(8pi)),p_e electron pressure
c     p_e=dense*tempe
c     bmod      ! tl
      btot=bmod*1.d4       ! gauss
      p_e=tempe*dense*1.6d+4	!(dyne/cm)
      p_b=btot**2/(8.d0*pi)	!magnetic pressure  (dyne/cm)
      beta_e=p_e/p_b

      mu=510.94d0/tempe
      mu2=mu*mu
c----------------------------------------
c     G=(N_par**2-S)/A2+
c       {1+omeg_pe**2*Z_r/(omega_e**2*abs(Z)**2*(N_par**2-S)}*
c       (m_e*c**2/T_e)**2*(omega_e*D/omega*abs(eps33))**2/(A1*A2)
c----------------------------------------
c     p_=P,d_=D,s_=S,  P,D,S are like in the article
      xe=wpw_2(x,y,z,1)
      ye=wcw(x,y,z,1)
      s_=1.d0+xe/(ye*ye)
      d_=0.d0
      do i=2,nbulk
        xi=wpw_2(x,y,z,i)
c        write(*,*)'absorpfd i,xi',i,xi
        yi=wcw(x,y,z,i)
        p=1.d0/(1.d0-yi*yi)
        s_=s_-xi*p
        d_=d_+xi*p/yi
      enddo
      p_=-xe
c---------------------------------------
c      rpr=real(cksi_e)
c      argz=cmplx(rpr,0.)
      argz=dcmplx(cksi_e,0.d0) ! cksi_e=cvac/(cnpar*ve) Can be 1/0 ?
c      call zfunc(argz,zf,zfp)
cSm021022
c      call czeta(argz,zf,zfp,ierror)
      call CZETA0(cnpar,argz,zf,zfp,ierror)
c      rpr=real(zf)
c      rpr=dble(zf)
      rpr=dreal(zf)
      
      zf_r=dble(rpr)      !real(Z(cksi_e))
c      rpim=aimag(zf)      
      rpim=dimag(zf)      
      zf_im=dble(rpim)	  !imag(Z(cksi_e))

      cksi_e2=cksi_e*cksi_e
      z_r=2.d0*cksi_e2*(1.d0+cksi_e*zf_r) ! article formula (18)
      z_im=2.d0*cksi_e2*cksi_e*zf_im	  ! article formula (18)
      cmodez_2=z_r*z_r+z_im*z_im
      eps33_r=xe*z_r			  ! article formula (5)
      eps33_im=xe*z_im			  ! article formula (5)
      eps33_2=eps33_r*eps33_r+eps33_im*eps33_im
      cnpar2=cnpar*cnpar
      p1=cnpar2-s_			  ! for (15),(16),(17)
      dp1=1.d0/p1			  ! for	(15),(16),(17)
      z_rdz2=z_r/(cmodez_2)
      p2=xe*z_rdz2*dp1/(ye*ye)
      a1_=p1+p2*cnpar2			  ! (16)
      a2_=p1+p2*(cnpar2+s_)		  ! (17)
      g_=p1/a2_+(1.d0+p2)*mu2*ye*ye*d_*d_/(eps33_2*a1_*a2_) !(15)
      cnprim_e=0.25d0*spi*cnper*beta_e*cksi_e*dexp(-cksi_e2)*g_	!(13)
cSm021021
      cnprim_e=dabs(cnprim_e)     

      cn1=p1/a2_ !(15)
      cn2=(1.d0+p2)*mu2*ye*ye*d_*d_/(eps33_2*a1_*a2_) !(15)

c end of the electron absorption calculations
      ci=dcmplx(0.d0,1.d0)
      cmue=0.5d0*(cnper*ve/(ye*cvac))**2
      ckappae2=xe*2.d0*cmue*cksi_e*(zf_r+ci*zf_im)  ! formula (7)
   
      xe=wpw_2(x,y,z,1)
      ye=wcw(x,y,z,1)

cSm021018
c-----tensor for the electrons with the hot correction
c     The hot correction is in eps(2,2),eps(3,3),eps(2,3)
c-----------------------------------------------------
      do i=1,3
         do j=1,3
         reps(i,j)=czero
         enddo
      enddo

      reps(1,1)=1.d0+xe/ye**2 !part of formula(2)
cSm021018 
c      reps(1,2)=-ci*xe/ye      !part of formula(4) 
cSm050207
c      reps(1,2)=xe/ye      !part of formula(4)
cSm080916
      reps(1,2)=-ci*xe/ye      !part of formula(4)
      reps(2,2)=reps(1,1)+ckappae2                  !part of formula (3)
      reps(3,3)=xe*(z_r+ci*z_im)                            !formula (5)
c      reps(2,3)=-ci*0.5*dsqrt(2.d0*cmue)*reps(3,3)/cksi_e  !formula (6)
cSm021018
      reps(2,3)=-ci*0.5*dsqrt(2.d0*cmue)*reps(3,3)/cksi_e
cSm021022
c      if(cnpar.lt.0.d0) reps(2,3)=-reps(2,3)
         
      cnprim_s(1)=cnprim_e

c      do i=2,nbulk
c        xi=wpw_2(x,y,z,i)
c        yi=wcw(x,y,z,i)
c        reps(1,1)=reps(1,1)+xi/(yi**2-1.d0)
c        reps(2,2)=reps(2,2)+xi/(yi**2-1.d0)
c        reps(1,2)=reps(1,2)-ci*yi*xi/(yi**2-1.d0)
c      enddo
c***********************************************************************
c ion absorption:
c-----------------------------------------------------------------
c     cksi_il=(1-l*yi)/(N_par*Vi/cvac)

c     eps11_i, eps22_i, eps12_r, - only with ion contributions
      cnprim_i=0.d0
      do i=2,nbulk

         eps11_i=0.d0
         eps22_i=0.d0
         eps12_r=0.d0

cSm050209 The elements eps11_r, eps22_r, eps12_i are used only
c        in the dielectric tensor reps()
c        They are not used in absorption (cnprim) calculations
         eps11_r=0.d0
         eps22_r=0.d0
         eps12_i=0.d0

c	 vi=sqrt(2Ti/mi),    vi in (cm/sec),Ti(keV)
         vi=1.87d9*dsqrt(tempiar(i)/dmas(i))
         xi=wpw_2(x,y,z,i)
         yi=wcw(x,y,z,i)
         rho_larm= vi/(yi*omega) ! [cm] Larmor radius of ion species
         !YuP[11-2016] limit the value of rho_larm by max value 
         !(set in namelist &dispers in genray.in):
         !rho_larm_max [cm] is the Upper limit for Larmor radius, 
         ! used for absorption calculation.
         ! A very large number means no limit.
         ! Recommended: In FRC run (when B goes through zero)
         ! set it to the distance between null point r0 and separatrix rs. 
         rho_larm= min(rho_larm, rho_larm_max0)
         !cmui=0.5d0*(cnper*vi/(yi*cvac))**2 ! =0.5(Kperp*rho_i)^2
         !Alternatively (equivalently), expressed through rho_larm:
         cmui=0.5*(cnper*omega/cvac *rho_larm)**2
         expmui=dexp(-cmui)
         if(cnpar.ne.0.d0) then
            cksi_0=cvac/(cnpar*vi) ! can be 1/0 ?
         else
            cksi_e=0.d0
         endif
            
         if( abs(cnpar*vi).lt.(cvac*1.e-8) )then ! w/(|Kpar|*Vti) >1e8
          write(*,*)'absorpfd_xyz/WARNING: tiny npar, npar*vi/c=',
     +     cnpar,cnpar*vi/cvac
          !pause
         endif
cSm050209
c         cksi_0=dabs(cksi_0)

         p22=xi*cksi_0*expmui	       ! for (2),(3),(4)
         p11=p22/cmui		       ! for (2)
c------------------------------------------
c        izi=1 to calculate modified Bessel functions
c        BESSEL FUNCTIONS AND ARGUMENTS ARE REAL
         ize=1
         call beslci1(cmui,0.d0,nb,ize,br,bi,ncalc)
         p=10.d0*dabs(cnpar)*vi/cvac
         l0=1/yi ! only for printout: l0=[omega/omega_ci]
         l1=(1-p)/yi 
         l2=(1+p)/yi+1 ! approximately 10*|Kpar|*rho_Larmor_ion
         !YuP[04-2016] The above determination for the range of harmonics
         !seems to be not good sometimes. In particular, when npar~0,
         !then l1=omega/omega_c, and l2=l1+1, so only two integers.
         !Better to extend the range by one harmonic  number:
         l1=l1-1 ! YuP
         l2=l2+1 ! YuP
         !YuP: no detectable change from this l1,l2 extension?
         !write(*,*)'absorpfd_xyz l1,l2=',l1,l2
         if(abs(L2).gt.nb-2) then
           write(*,*)'absorpfd_xyz: npar, npar*vi/c=',
     +     cnpar,cnpar*vi/cvac
           write(*,*)'absorpfd_xyz: L2>nb-2; Increase
     1	nb in absorpfd_xyz; L0=1/yi, L2, rho_larm,Kper*rho_larm=',
     +  L0,L2, rho_larm, cnper*omega/cvac *rho_larm
           !stop  ! YuP[10-2016] commented: let it continue.
           l2=min(l2,nb-2)
           l2=max(l2,-nb+2) ! lower limit for a negative l2 is -(nb-2)
         endif 

	 if(abs(L1).gt.nb-2) then
     	   write(*,*)'in absorpfd_xyz L1>nb-2 Increase
     1	   nb in the subr. absorpfd_xyz; L1,l0=',L1,l0 
           !stop  ! YuP[10-2016] commented: let it continue.
           l1=min(l1,nb-2)
           l1=max(l1,-nb+2) ! lower limit for a negative l1 is -(nb-2)
	 endif 

      do l=l1,l2
c          cksi_l=(1-l*yi)/(N_par*Vi/cvac)
         if(l.ge.0) then
            lmode=l
         else
            lmode=-l
         endif
         besm_l=br(lmode+1)	! |l| order modified Bessel function
         besm_lp=br(lmode+2)
         if (l.eq.0)then
             cksi_l=cksi_0 ! cksi_0=cvac/(cnpar*vi)
             cksi_l2=cksi_l*cksi_l
             zfi_l=spi*dexp(-cksi_l2)	 !Im(Z_func(cksi_l))
c-------------------------------------------------------------------------
cSm050209    The calculatios for eps11_r,eps22_r,eps12_i
c            They are not used in absorption calculations
             argz=dcmplx(cksi_l,0.d0)             
             call CZETA0(cnpar,argz,zf,zfp,ierror)
c             write(*,*)'absorpfd i,l,zfi_l,imag(zf)',i,l,zfi_l,dimag(zf)
             zfr_l=dreal(zf)
c-------------------------------------------------------------------------
c	     besmp_l is a derivative from modified Bessel function I_0
c            by argument cmui :I_0 prime=I_1
             besmp_l=besm_lp
             eps22_i=eps22_i+2.d0*p22*cmui*(besm_l-besmp_l)*zfi_l  !(3)

cSm050209    The calculations for eps22_r
c            They are not used in absorption calculations
             eps22_r=eps22_r+2.d0*p22*cmui*(besm_l-besmp_l)*zfr_l   !(3)
	   else ! l<0 or l>0
             if(cnpar.ne.0.d0) then
               cksi_l=(1.d0-l*yi)*cvac/(cnpar*vi) ! Can be 1/0 ?
             else
               cksi_l=0.d0
             endif               
             if( abs(cnpar*vi).lt.abs(cvac*(1.0-l*yi)*1.e-8) )then 
               ! abs(w-l*wc)/(|Kpar|*Vti) >1e8
               write(*,*)'absorpfd_xyz/WARNING: tiny npar',cnpar,' l=',l
               !pause
             endif
             cksi_l2=cksi_l*cksi_l
             zfi_l=spi*dexp(-cksi_l2)	 !Im(Z_func(cksi_l))
c-------------------------------------------------------------------------
cSm050209    The calculatios for eps11_r,eps22_r,eps12_i
c            They are not used in absorption calculations
             argz=dcmplx(cksi_l,0.d0)             
             call CZETA0(cnpar,argz,zf,zfp,ierror)
c             write(*,*)'absorpfd i,l,zfi_l,dimag(zf),dimag(zf)',
c     &                           i,l,zfi_l,dimag(zf),dimag(zf)
             zfr_l=dreal(zf)
             zfi_l=dimag(zf)
c-------------------------------------------------------------------------
             eps11_i=eps11_i+p11*l*l*besm_l*zfi_l		   !(2)
c	     besmp_l is a derivative from modified Bessel function I_l
c            by argument cmui :I_l prime=0.5*(I_(l-1)+I_(l+1))
             besm_lm=br(lmode)
             besmp_l=0.5d0*(besm_lp+besm_lm)
             eps22_i=eps22_i+2.d0*p22*cmui*(besm_l-besmp_l)*zfi_l  !(3)
             eps12_r=eps12_r+p22*l*(besm_l-besmp_l)*zfi_l	   !(4)

cSm050209    The calculatios for epss11_r,eps22_r,eps12_i
c            They are not used in absorption calculations
             eps11_r=eps11_r+p11*l*l*besm_l*zfr_l                  !(2)
             eps22_r=eps22_r+2.d0*p22*cmui*(besm_l-besmp_l)*zfr_l  !(3)
             eps12_i=eps12_i-p22*l*(besm_l-besmp_l)*zfr_l	   !(4)

         endif ! l=0 or not
        enddo !l  (bessel function order)

        eps22_i=eps22_i+eps11_i					   !(3)

        cnprim_s(i)=0.5d0*cnper*(
     1  (eps11_i-2.d0*d_*eps33_r*eps12_r/eps33_2)/a2_+
     1  (-p1*(eps22_i+eps11_i)-2.d0*d_*eps12_r)*a1_/((p1*p1-d_*d_)*a2_))
     
        !if(cnprim_s(i).lt.0.d0)then
        !write(*,*)'absorpfd_xyz: i,cnprim_s(i)=',i,cnprim_s(i)
        !YuP[04-2016] Sometimes cnprim_s are negative,
        !especially for ion damping (case of HHFW, for example)
        !If they are negative, they are reset to abs() in prep3d_xyz
        !just after calling absorpfd_xyz.
        !Before 04-11-2016: the negativity of total(e+i) cnprim was checked.
        !After  04-11-2016: each of cnprim_e, cnprim_i is checked. 
        !The total can be positive, at the same time cnprim_i can be neg. 
        !endif

        cnprim_i=cnprim_i+cnprim_s(i)

c end of the ion absorption calculations

c-----add the hot ions terms to the tensor
cSm050209
c      reps(1,1)=reps(1,1)+eps11_i
c      reps(2,2)=reps(2,2)+eps22_i
c      reps(1,2)=reps(1,2)+ci*eps12_r

        reps(1,1)=reps(1,1) + ci*eps11_i + eps11_r
        reps(2,2)=reps(2,2) + ci*eps22_i + eps22_r
        reps(1,2)=reps(1,2) + eps12_r    + ci*eps12_i


      enddo !i  (ion species)


c      write(*,*)'after cnprim_i',cnprim_i

      reps(2,1)=-reps(1,2)
      reps(3,1)=reps(1,3)
      reps(3,2)=-reps(2,3)

c***********************************************************************
      return
      end


c======================================================================
c======================================================================

c        ********************** absorplh_xyz ****************
c        *                      -----                       *
c        *  this subroutine calculates the imaginary part   *
c        *  of refractive index using formula from:   	    *
c        *  P.Bonoli, Linear theory of lower hybrid heating *
c        *  IEEE Transaction on Plasma Science,Vol. PS-12,  *
c        *  No.2,(1984), p 95-107                           *
c        *  See also Paul T. Bonoli and Ronald C. Englade,  *
c        *  Phys. Fluids 29, 2937 (1986).                   *
c        ****************************************************
c         input parameters: u(6)  -x,y,z,Nx,Ny,Nz  
c                           cnpar -N_par			  
c                           cnper -N_per (real part)   
c                           temper -electron temperature  (keV)  
c                           densty -electron density   (10*13/cm**3
c                           tempiar -ions temperature  (keV)      *
c                           b_x,b_y,b_z,bmod-magnetic field	  *
c                           nbulk - number of plasma components	  *
c                           frqncy- wave friquency [GHz]	  *
c                           zeff  - plasma charge      		  *
c         output parameter: cnprim_e -imaginary part of N_perp    *
c                             (uncollisional electron absorption) *
c                           cnprim_i -imaginary part of N_perp    *
c                                (uncollisional ion absorption)   *
c                           cnprim_cl-imagionary part of N_perp	  *
c                                collisional(ions and electron)	  *
c-----------------------------------------------------------------*
      subroutine absorplh_xyz(u,cnpar,cnper,temper,densty,tempiar,
     1 b_x,b_y,b_z,nbulk,bmod,frqncy,zeff,
     1 cnprim_e,cnprim_i,cnprim_cl)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'ions.i'
      dimension u(6),deruu(6)
      dimension tempiar(nbulka)
********************************************************************
c electron absorption (Landau damping, formula (21a) in the Bonoli article)
c Di_e=2*sqrt(pi)*(omegape/omega)**2*(N_par*N_perp)**2*abs(x_oe**3)
c      *exp(-x_oe**2)
c--------------------------------------------------------------------
      x= u(1) 
      y= u(2) 
      z= u(3) 
      cn_x= u(4) ! Nx
      cn_y= u(5) ! Ny
      cn_z= u(6) ! Nz
c--------------------------------------------------------------------
      pi=4.d0*datan(1.d0)
      spi=dsqrt(pi)
      cnpar2=cnpar*cnpar
      cnpar4=cnpar2*cnpar2
      cnper2=cnper*cnper
      cnper4=cnper2*cnper2
c----------------------------------------
c     cvac - light velocity in (cm/sec)
      cvac=2.997930d+10
c----------------------------------------
c     ve=(sqrt(Te/me)),    ve in (cm/sec),Te(keV)
      ve=1.87d9*dsqrt(0.5d0*temper)
      sqrt2=dsqrt(2.d0)
c----------------------------------------
c     x_oe=omega/((sqrt(2))*k_par*ve)=cvac/(sqrt(2)*N_par*ve)
      x_oe=cvac/(sqrt2*cnpar*ve)
      x_oe2=x_oe*x_oe
      x_oe3=x_oe2*dabs(x_oe)
      ye=wcw(x,y,z,1)
      xe=wpw_2(x,y,z,1)
c--------------------------------------
c     di_e is imaginary part of the dispersion function(uncollisional)
      di_e=2.d0*spi*xe*cnper2*cnpar2*x_oe3*dexp(-x_oe2)
c*********************************************************************
c     di_i is imaginary part of the dispersion function(uncollisional),
c     formula (21b) in the Bonoli article.
c----------------------------------------------------------
      di_i=0.d0

c      write(*,*)' absorplh.f nbulk',nbulk
      do i=2,nbulk
        tempi=tempiar(i)
c        write(*,*)'i,tempi,dmas(i)',i,tempi,dmas(i)
c-------vi=(sqrt(Ti/mi)),    vi in (cm/sec),Ti(keV)
        vi=1.87d9*dsqrt(0.5d0*tempi/dmas(i))
       
c-------x_oi=omega/((sqrt(2))*k_par*vi)=cvac/(sqrt(2)*N_par*vi)
        x_oi=cvac/(sqrt2*cnpar*vi)
        x_oi2=x_oi*x_oi
        x_oi3=x_oi2*dabs(x_oi)
        x_oi4=x_oi2*x_oi2
        yi=wcw(x,y,z,i)
        xi=wpw_2(x,y,z,i)

        sum=0.d0
        
        pdi_i=2.d0*spi*xi*cnper4*x_oi3*dexp(-x_oi2)
        psum=yi*dabs(x_oi)/spi
	if (pdi_i.lt.1.d-10) goto 20
        n0=1/yi
        pn=10.d0*sqrt2*vi*dabs(cnpar)/cvac
        x_ni0=(1.d0-n0*yi)*cvac/(sqrt2*cnpar*vi)
        n1=(1-pn)/yi-1
        x_ni1=(1.d0-n1*yi)*cvac/(sqrt2*cnpar*vi)
        n2=(1+pn)/yi+1
        x_ni2=(1.d0-n2*yi)*cvac/(sqrt2*cnpar*vi)
        
c        write(*,*)'absorplh.f n1,n2',n1,n2 
        do n=n1,n2
c--------------------------------------
c         x_ni=(omega-n*omegas_ci)/sqrt(2)*k_par*vi)
          x_ni=(1.d0-n*yi)*cvac/(sqrt2*cnpar*vi)
          x_ni2=x_ni*x_ni
          if (x_ni2.gt.100)then
            goto 10
          endif
          sum=sum+dexp(-x_ni2)
 10       continue
        enddo !n
c        write(*,*)'absorplh.f,pdi_i,sum',pdi_i,sum
 20     sum=sum*psum

        di_i=di_i+pdi_i*sum

      enddo !i
c*********************************************************************
c     di_ic is the imaginary part of the dispersion function(collisional)
c     Collisional, ions and electron damping coefficients: (cgs units)
c     formula (26) from the Bonoli article
c     1.47e-9 = (4/3)*sqrt(2*pi)*16*q**4/(sqrt(me)*(1.6e-9 erg/keV)**1.5)
      rnuei = 1.47d-9*1.d+13*densty*zeff/temper**1.5 ! [1/sec]
      omega=2.d0*pi*frqncy*1.d9	  ! [1/sec], frqncy in GHz
      di_ic=rnuei*(xe*cnper2/(ye*ye)+xe*cnpar2)*cnper2/omega
c--------------------------------------------------------------
c     calculation of (dD_real/dN_par)
c     D is from the P.Bonoli article  , formula (7),(9) and (17) 
      xi=0.d0
      if(nbulk.gt.1) xi=wpw_2(x,y,z,2)
      epsper=1.d0+xe/(ye*ye)-xi
      epspar=1.d0-xe-xi
      epsxy=xe/ye
c------------------------------------------
c p6 calculations
      ye2=ye*ye
      ye4=ye2*ye2
      pe=ve/cvac
      pe2=pe*pe
c     vi=(sqrt(Ti/mi)),    vi in (cm/sec),Ti(keV)
      tempi=tempiar(2)
      vi=1.87d9*dsqrt(0.5d0*tempi/dmas(2))
      pi=vi/cvac
      pi2=pi*pi
      p6=-(3.d0*xi*pi2+0.75d0*xe/ye4*pe2)
c------------------------------------------
      p4=epsper
      p2=(epsper+epspar)*(cnpar2-epsper)+epsxy*epsxy
      p0=epspar*((cnpar2-epsper)*(cnpar2-epsper)-epsxy*epsxy)
c------------------------------------------
      dddnper=(6.d0*p6*cnper4+4.d0*p4*cnper2+2.d0*p2)*cnper
      cnprim_e=-(di_e/dddnper)
      cnprim_i=-(di_i/dddnper)
      cnprim_cl=-(di_ic/dddnper)
      cnprim_e=dabs(cnprim_e)
      cnprim_i=dabs(cnprim_i)
      cnprim_cl=dabs(cnprim_cl)
      return
      end

c======================================================================
c======================================================================

      subroutine absorbed_power_using_ql_flux_xyz(cnpar,cnper,
     + x,y,z, fluxn, power, del_s_poloidal,
     & absorbed_power) !-> out

      implicit none
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'cefield.i'    

c-----input
      real*8
     &cnpar,cnper,       !parallel and perpendicular refractive index
     &x,y,z,             !ray cartesian coordinates [meter]
     &power,             !power in ray channel [erg/sec]
     &del_s_poloidal     !poloidal length of the ray element []

      real*8 fluxn       !power flux at unit electric field |E|=1 
                         !flux=B~(i)*B(i)+E~(i)*d(omega*eps(i,j))/domega*E(j)
                         !fluxn=0.5*dreal(flux*dconjg(flux))*vgrpol 
                         !vgrpool is a poloidal group velocity
                         !normalized to clight   
c-----output          
      real*8
     &absorbed_power     !erg/sec the power abdsorbed at the ray element
   
c-----external
      real*8
     &psi_rho, psif_xyz,
     &bmin_psi

c-----locals
      real*8
     &psi_loc,unorm,
     &cd_efficiency,pow_dens_0_tilda,cd_dens,
     &b_pol,bmin, r,
     &clight,s_poloidal_tilda,
     &absorbed_power_dev_s_pol               !erg/[sec*cm] the power abdsorbed at the ray element 
                                             !of poloidal length delta_s_pol [cm]
                                             
      pi=4.d0*atan(1.d0)
      clight=2.99792458d10              !light speed [cm/sec]

      psi_loc=psif_xyz(x,y,z)

      call QL_power_xyz(
     &psi_loc, x,y,z,
     &cnpar,cnper,cex,cey,cez,
     &unorm,
     &pow_dens_0_tilda)

      bmin= bmin_psi(psi_loc)
c-----total group velocity
c      vgr_poloidal=dsqrt(vgr(1)**2+vgr(2)**2)! poloidal group velocity /clight
c      vgr_poloidal=vgr_poloidal*clight !sm/sec]

      s_poloidal_tilda= fluxn*clight/(8.d0*pi) 

      absorbed_power_dev_s_pol=power*(bmod/bmin)*
     &         (pow_dens_0_tilda/s_poloidal_tilda)
           
      absorbed_power=absorbed_power_dev_s_pol*del_s_poloidal ![erg/sec]
      return
      end


c======================================================================
c======================================================================

      subroutine QL_power_xyz(
     &psi_in, x,y,z,
     &n_parallel,n_perp,e_x,e_y,e_z,
     &unorm, 
     &pow_dens_0_tilda)
c------------------------------------------------------
c     calculates bounce averaged power:  pow_dens 
c     at unit wave elecric field vector |E|=1
c------------------------------------------------------
      implicit none
      include 'param.i'
      include 'one_nml.i'
      !!! include 'adj.i'
      include 'write.i'      
c-----input
      real*8  
     &x,y,z,               !cartesian coords. [m]    
     &n_parallel,          !parallel refractive index
     &n_perp,              !perpendicular refractive index 
     &psi_in               !poloidal flux
      complex*16 
     &e_x,
     &e_y,                 !electric field polarization in Stix frame
     &e_z                  !at |E|=1
c-----output
      real*8
     &unorm,               !sqrt(T/m) cm/sec
     &cd_efficiency,
     &pow_dens_0_tilda,     ![1/sec], see manual: absorbed power calculation
c                          !Power density normalized to m*u_norm**2*n/tau_n 
     &cd_dens              !current drive density normalized to q*u_norm*n
                           !Here n is the density  
c-----external
      real*8 temperho,densrho,
     &bmax_psi,bmin_psi,psi_rho,rhopsi,
     + bxyz,wcw,wpw_2, tempe_xyz,dense_xyz

c-----local
      integer n_radial0,n_harm_adj,
c     &i_resonance_curve_integration_method,
     &kmax,
     &i_calculate_CD     ! =0, calculate power only (no CD) 
                         ! using input arguments computed at
                         ! given space point
                         ! In this case n_radial0 can be arbitary
      real*8 u0,     
     &rho_small,           !small radius
     &temp_kev,            !electron temperatura [KeV]
     &dense,               !electron dencity 10**13/cm**3     
     &b_ratio,             !B/B_min
     &pi,bmin,bmax,deltb,psi,th0max,sin_trap,bmod,
     &d_chi_du0,d_chi_dtheta0,u0max,chi,
     &clight,             !light speed [cm/sec]
     &arg1,cln,
     &y_loc,              !electron wc/w at point (x,y,z)     
     &power_nharm_adj,CD_nharm_adj, !under integral functions
     

     &p_par_rl(2),         !for check max and min p_par
     &k_1_kev,             !egrs in 1 KeV      (erg)
     &charge_electron,     !electron charge (statcoulomb)
     &mass_e,              !electron mass [g]
     &theta_temperature,   !m_e*c^2/T_e
     &pow_dens,            
                           !using QL flux
     &omega                !2*pi*frqncy*1.d9 [1/sec]

      k_1_kev=1.6022d-9          !egrs in 1 KeV      (erg)
      charge_electron=4.8032d-10 !electron charge (statcoulomb)
      pi=4.d0*datan(1.d0)
      clight=2.99792458d10              !light speed [cm/sec]    
      mass_e=9.1094d-28         !electron rest mass (g)

cyup      rho_small=rhopsi(psi_in)
cyup      temp_kev=temperho(rho_small,1)
cyup      dense=densrho(rho_small,1) 
      temp_kev= tempe_xyz(x,y,z,1)  ! YuP added
      dense=    dense_xyz(x,y,z,1)  ! YuP added
ccc      unorm=dsqrt(temp_kev/t)*dsqrt(1.6022d0/0.91095d0)*1.d9 ![sm/sec] sqrt(T/m)
      psi=psi_in
      bmax= bmax_psi(psi)
      bmin= bmin_psi(psi)
	     
      deltb = (bmax-bmin)/bmin
      th0max = atan2(1.0d0,sqrt(deltb))
      sin_trap=dsqrt(bmin/bmax)

      pow_dens=0.d0
      cd_dens=0.d0

      bmod= bxyz(x,y,z)  !-> get bmod
      b_ratio=dmax1(bmod/bmin,1.d0)

      y_loc= wcw(x,y,z,1)

      pow_dens=2.d0*pi*pow_dens
      cd_dens=-2.d0*pi*cd_dens

c      write(*,*)'pow_dens',pow_dens
c      write(*,*)'cd_dens',cd_dens
     
      if(pow_dens.ne.0.d0) then
        cd_efficiency=cd_dens/pow_dens
      else
        cd_efficiency=0.d0
      endif

      write(*,*)'in CD_adj_efficiency_1 cd_efficiency=',
     &cd_efficiency
         
c-----cln -Coulomb Ln
      arg1 = 1.d3/temp_kev*dsqrt(10.d0*dense)	!temp KeV, dense 10**13/cm**3
      cln=24.d0-dlog(arg1)

c-----adj efficiency is normalized to eta_0_adj=e*tau_n/(m*u_t) 
c     tau_n=T**(3/2)*sqrt(m)/(4*pi*e**4*density*cln)
c     u_t=cqrt(T/m)
c     eta_0_adj=T/[4*pi*e**3*n*cln]|cgs=T_kev/(n_13*cln)*
c           [k_1_kev/(4*pi*charge_electron**3*10**13)]           !cgs
c      coef=1.6022d+8/(4*pi*4.8032**3)          !cgs  (statampere/cm**2)/(erg/(sec*cm**3))
c      coef=1.6022d-1/(4*pi*4.8032**3*3)=3.8352d-5!     (A/cm**2)/(erg/(sec*cm**3))
c      eta_0_adj=(temp_kev/(cln*dense))*coef

      cd_efficiency=cd_efficiency*
     &3.84d0*temp_kev/(cln*dense)	       !  (A/m**2)/(joule/(sec*m**3))
c      write(*,*)'temp_kev,dense,cln,cd_efficiency',
c     &temp_kev,dense,cln,cd_efficiency
      cd_efficiency=cd_efficiency*1.d-5  !  (A/cm**2)/(erg/(sec*cm**3))
      write(*,*)'temp_kev,dense,cln,cd_efficiency',
     &temp_kev,dense,cln,cd_efficiency

c---------------------------------------------------------------------
cSAP081201
c     
c     The used relativistic Maxwellian distribution was normalized
c     to  unit density electron density:
c     4*pi*integral{f_m*(p/mc)^2*d(p/mc)}=1
c---------------------------------------------------------------------
  
      omega=2.d0*pi*frqncy*1.d9 ![1/sec]
c      write(*,*)'frqncy,omega',frqncy,omega
c---- power_dens_0_tilda has dimension [1/sec]
c      write(*,*)'pow_dens',pow_dens
c      write(*,*)'dense*1.d13',dense*1.d13
c      write(*,*)'unorm,clight,unorm/clight',unorm,clight,unorm/clight
      pow_dens_0_tilda=pow_dens*
     &dense*1.d13*      !the power was calculated for unit density
                        !This term transforms power to local density
     &(unorm/clight)*   !it was the integration by (p_perp/mc)*d(p_perp/mc)
     &charge_electron**2/mass_e/omega !in subroutine intgr_relt_power_cd_adj
                        !This transforms this integration to
                        !(p_perp/(m*unorm)*)d(p_perp/(m*unorm))
                        !and uses all other normalizations
                        !in under the integral term
c      write(*,*)'charge_electron**2/mass_e/omega [cm^3/sec] ',
c     &           charge_electron**2/mass_e/omega
c       write(*,*)'pow_dens_0_tilda',pow_dens_0_tilda,'[1/sec]'
      return
      end

     

c======================================================================
c======================================================================

c        ********************** flown_xyz ********************
c        *                      ------                       *
c        * this subroutine calculates the wave energy flow   *
c        * for normalized electric field modeE=1             *
c        *****************************************************
c
c--------------------------------------------------------------------
c								    
c        input parameters					    
c      								    
c      x,y,z - coordinates of the point where the electric field 
c                 is calculated.      		         	    
c								    
c      cnx,cny,cnz - components of  wave  refractive index. 
c       							    
c      complex components of dielectric tensor reps(i,j) are 	    !
c      in common block 'eps'.These components are created in        !
c      subroutine hamilt_xyz					    !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
c								    !
c        output parameters					    !
c								    !
c        cflown                               			    !
c        cflown=B~(i)*B(i)+E~(i)*d(omega*eps(i,j))/domega*E(j)	    !
c        B~ =conjg(B), E~=conjg(E)				    !
c--------------------------------------------------------------------
      subroutine flown_xyz(x,y,z, cnx,cny,cnz, cflown)
c     cex,cey,cez - complex electr.field polarization E_x/E,E_y/E,E_z/E
c     from common bloc /cefield/
c     cbx,cby,cbz - complex magnet.field polarization B=N*E modeE=1
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'eps.i'
      include 'cefield.i'
      dimension vp(nbulka),wp(nbulka)
      double complex cbx,cby,cbz
      double complex cepsp,cepsm,dwepsdw,ce,cec,cb,cbc,cflown
      dimension cepsp(3,3),cepsm(3,3),dwepsdw(3,3)
      dimension ce(3),cec(3),cb(3),cbc(3)
      complex*16 ceps_herm(3,3) !Hermitian part of reps
c--------------------------------------------------------------
c     ce - wave electric field
c     cec- comlex conjugate wave electric field
c     cb - wave magnetic field
c     cbc- comlex conjugate wave magnetic field
c---------------------------------------------------------------
      bmod=bxyz(x,y,z) ! get b and derivs of b
      gam1=gamma1_xyz(x,y,z,cnx,cny,cnz)
      ds=dsin(gam1)
      dc=dcos(gam1)
      dcn= dsqrt(cnx*cnx+cny*cny+cnz*cnz)
      cnpar=dcn*dc
      cnper=dcn*ds
c      write(*,*)'in flown dcn,cnpar,cnper',dcn,cnpar,cnper
c      write(*,*)'cex,cey,cez',cex,cey,cez
c----------------------------------------------------------
c    complex components of the wave magnetic field
      cbx=-cnpar*cey
      cby=-cnper*cez+cnpar*cex
      cbz=cnper*cey
c      write(*,*)'cbx,cby,cbz',cbx,cby,cbz
c----------------------------------------------------------
      cb(1)=cbx
      cb(2)=cby
      cb(3)=cbz
      ce(1)=cex
      ce(2)=cey
      ce(3)=cez
c      write(*,*)'cb(i)',cb(1),cb(2),cb(3)
c      write(*,*)'ce(i)',ce(1),ce(2),ce(3)
c-----------------------------------------------------------
      !step=1.d-7 ! prep3d_xyz->flown_xyz()
      hw=der_f ! ! prep3d_xyz->flown_xyz() was hw=step
      hw=hw*frqncy
      pi=3.1415926d0
      hfrqnc=hw  ! =der_f*frqncy
      w0=frqncy
c-----------------------------------------------------------
      gam=gamma1_xyz(x,y,z,cnx,cny,cnz)
      do 1 i=1,nbulk
         vp(i)=v(i)
         wp(i)=w(i)
 1    continue
      frqncpl=frqncy+hfrqnc
      df=frqncy/frqncpl
      do 2 i=1,nbulk
         v(i)=vp(i)*df*df
         w(i)=wp(i)*df
 2    continue

c************************************************
      cnxplus=cnx*df
      cnyplus=cny*df
      cnzplus=cnz*df
      cnt2= cnxplus**2 + cnyplus**2 + cnzplus**2
      hp=hamilt_xyz(x,y,z,cnt2)
c*************************************************
      w0p=w0+hw
cSAP081122 calculates hermitian part of the dielectric tensor reps
      call hermitian_part(reps,ceps_herm)
      do 3 i=1,3
        do 3 j=1,3
           cepsp(i,j)=w0p*ceps_herm(i,j)
 3    continue
c----------------------------------------------------------
      frqncmn=frqncy-hfrqnc  ! ! = frqncy*(1.d0-der_f)
      df=frqncy/frqncmn
      do 4 i=1,nbulk
        v(i)=vp(i)*df*df
        w(i)=wp(i)*df
 4    continue
c************************************************
      cnxminus=cnx*df
      cnyminus=cny*df
      cnzminus=cnz*df
      cnt2= cnxminus**2 + cnyminus**2 + cnzminus**2
      hm=hamilt_xyz(x,y,z,cnt2)
c*************************************************
      w0m=w0-hw
       call hermitian_part(reps,ceps_herm)

cSAP081122 calculates hermitian part of the dielectric tensor reps
      do 5 i=1,3
         do 5 j=1,3
           cepsm(i,j)=w0m*ceps_herm(i,j)
c          write(*,*)'after hm i,j,cepsm(i,j)',i,j,cepsm(i,j)
 5    continue
c-----------------------------------------------------------
      do 6 i=1,nbulk
        v(i)=vp(i)
        w(i)=wp(i)
 6    continue
c-----------------------------------------------------------
      do 7 i=1,3
         do 7 j=1,3
           dwepsdw(i,j)=(cepsp(i,j)-cepsm(i,j))/(2.d0*hw)
 7    continue
      do 8 i=1,3
           cec(i)=dconjg(ce(i))
           cbc(i)=dconjg(cb(i))
 8    continue

      cflown=dcmplx(0.0d00,0.0d00)
      do 9 i=1,3
         cflown=cflown+cb(i)*cbc(i)
         do 9 j=1,3
 	    cflown=cflown+cec(i)*dwepsdw(i,j)*ce(j)
 9    continue

      return
      end


c======================================================================
c======================================================================

c        ********************** prepebw_xyz ***********************
c        *                      -----                             *
c        *  prepebw -subroutine to prepare the  data for studying  *
C        *  ebw                                                   *
c        **********************************************************
      subroutine prepebw_xyz(u,is)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'write.i'
      include 'three.i'
      include 'n_parb.i'
      double precision modr,modrold
      dimension u(6)
c---------------------------------------------
      x= u(1) 
      y= u(2) 
      z= u(3) 
      cnx= u(4) ! Nx
      cny= u(5) ! Ny
      cnz= u(6) ! Nz
c---------------------------------------
c     the creation of the array wthet(i): poloidal angle
c     along the ray trajectory (ir radians)
      r= sqrt(x*x+y*y)
      zcomp=z-zma
      rcomp=r-rma
      modr=dsqrt(zcomp*zcomp+rcomp*rcomp)  

      btot(is)=bmod
      gradpsi(is)=dsqrt(dpdxd**2+dpdyd**2+dpdzd**2)
      gradpdr(is)=gradpsi(is)/r
c--------------------------------------------------------------------
c     e_theta=e_psi x e_phi		   (x is a vector production)
c     e_psi=(e_z*dpsidz+e_r*dpsidr)/sqrt(dpsidz**2+dpsidr**2)
c     e_theta=(e_z*dpsidr-e_r*dpsidz)/sqrt(dpsidz**2+dpsidr**2)
c     Npol=N * e_theta
c--------------------------------------------------------------------
      !cosphi=x/r
      !sinphi=y/r
      ! bpoloid = (B.[grad(psi) X e_phi]) / |gradpsi| 
      ! For now ok; needs work?
      bpoloid(is)=(bz*(dpdxd*x+dpdyd*y)-(bx*x+by*y)*dpdzd)
     /               /(r*gradpsi(is))
      wnpol(is)=(cnz*(dpdxd*x+dpdyd*y)-(cnx*x+cny*y)*dpdzd)
     /               /(r*gradpsi(is))
c--------------------------------------------------------------------
      if (is.eq.1) then
         if(zcomp.ge.0.d0) then
            wthet(1)=dacos(rcomp/modr)
         else
            wthet(1)=-dacos(rcomp/modr)
         endif
      else

      rold= sqrt(xold*xold+yold*yold)
      zoldcomp=zold-zma
      roldcomp=rold-rma
      modrold=dsqrt(zoldcomp*zoldcomp+roldcomp*roldcomp) !|r_old|

c--------------------------------------------------------------------
c            [r*r_old]=|r|*|r_old|*sin(deltheta)
c            |e_z        e_r       e_phi |
c            |zcomp	 rcomp	   0	 |=e_phi*|r|*||r_old|*sin(deltheta)
c            |zoldcomp   roldcomp  0	 |
	 sindelth=(zcomp*roldcomp-rcomp*zoldcomp)/(modr*modrold)
c--------------------------------------------------------------------
c            (r*r_old)=|r|*|r_old|*cos(delttheta)
	 cosdelth=(zcomp*zoldcomp+rcomp*roldcomp)/(modr*modrold)
	 if(cosdelth.gt.1.d0-1.d-15) cosdelth=1.d0-1.d-15
	 if(cosdelth.lt.-1.d0+1.d-15) cosdelth=-1.d0+1.d-15

c--------------------------------------------------------------------
         if (sindelth.gt.0.d0) then
            deltheta=dacos(cosdelth)
         else
            deltheta=-dacos(cosdelth)
         endif
         wthet(is)=wthet(is-1)+deltheta
      endif ! is

      wmdevr(is)= (by*x-bx*y)*(cny*x-cnx*y)/(r*r*bmod)  ! what for?

c--------------------------------------------------------------------
      xe=wpw_2(x,y,z,1)
      wxe(is)=dsqrt(xe)
c     wnrho is the radial component of the refructive index.
c     It is positive if it is directed along the gradient 
c     of the flux surface. 
      wnrho(is)=(cnx*dpdxd+cny*dpdyd+cnz*dpdzd)/gradpsi(is)
      gnpar(is)=bpoloid(is)/btot(is)
      wp(is)=gnpar(is)*wxe(is)
          
 10   continue
      return
      end ! prepebw_xyz

c======================================================================
c======================================================================

      subroutine set_nperpcom_xyz(nll,x,y,z,dmas)
c-----set the common block /nperpcom.i for the function hotnp
      implicit none
c-----input 
      double precision nll !the parallel refractive index N_par
      double precision x,y,z,r
      double precision dmas(*) 
c-----externals
      double precision tempe_xyz,tpoprho,vflowrho,rhopsi,wpw_2,wcw
c-----output the data for common block 
      include 'param.i' 
      include 'one.i' ! contains rho
      include 'nperpcom.i'
c-----local 
      integer j
      double precision psi

      nbulkc=nbulk
      nllc=nll

      if(nbulka.lt.nbulk) then
        write(*,*)'in forest.f in set_nperpcom nbulka.lt.nbulk'
        write(*,*)'nbulka=',nbulka,'nbulk=',nbulk
	write(*,*)'change parameter nbulka in file param.i'
        stop  
      endif
      
      do j=1,nbulk
         massc(j)=dmas(j)
         xc(j)=wpw_2(x,y,z,j) !-> finds rho (through dense_xyz)
         yc(j)=wcw(x,y,z,j)
         if(j.eq.1) yc(1)=-yc(1) ! negative wccw=(omega_ce/omega) for electrons
         tec(j)=tempe_xyz(x,y,z,j)*1.d+3 !(eV) averaged temperature
         tpopc(j)=tpoprho(rho,j)
         vflowc(j)=vflowrho(rho,j)         
      enddo
      return
      end

c======================================================================
c======================================================================

      subroutine p_c_prof_xyz(is,rhobegin,rhoend,cnpar,cnper,
     &cefldx,cefldy,cefldz)
c----------------------------------------------------------------
c     this subroutine calculates absorpt power profiles
c     calculates RF current  profiles	 (from deposited power)
c-----------------------------------------------------------------
c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none
      include 'param.i'
      include 'one.i'
      INCLUDE 'three.i'
      INCLUDE 'onetwo.i'
      INCLUDE 'write.i'
      INCLUDE 'gr.i'
      INCLUDE 'rho.i'
      INCLUDE 'five.i' ! contains rmax
      include 'oxb.i'  ! contains nrayelt_o_cutoff,
                       ! i_call_prep3d_in_output_at_i_ox_conversion_eq_1
     
c-----input
      real*8 rhobegin,rhoend,cnpar,cnper
      integer is
      complex*16 cefldx,cefldy,cefldz !polarization

c-----locals
      real*8 delpow_s,ppow_s, bavr
      dimension delpow_s(nbulka) !Added for indiv ion contrib[BH041009]
      dimension ppow_s(nbulka)  !Added for indiv ion contrib[BH041009]

      real*8 r_r,zmar,rmar,z_effr,tempr,denr,cnparr,yer,effic_r

      real*8 r_m,x_m,y_m,z_m, x_r,y_r,z_r, wss,
     +temp,ye,u1,den,z_eff,rholeft,rhoright,
     &rhomax,rhomin,eff_rho_max,eff_rho_min,hrho,delpower,delpow_i,
     &delpow_e,delpow_cl,del, delpower_del,
     +r0,x0,y0,z0,r0m,x0m,y0m,z0m,psiloc,zfacgeom,aspct,
     &delcurr,rho0,poloidlen,rho0_pol,delrho,
     &ppow,pcur,ppow_e,ppow_i,ppow_cl,
     &r_bmin_right,r_bmin_left,z_bmin_right,z_bmin_left,theta,psi

      integer i,kk,jbinmin,jbinmax,j
cSAP070831
      real*8 psi_loc,cos_theta_pol,sin_theta_pol,rho_loc,theta_pol
      integer n_theta_pol
      
      real*8 rmid, zmid, rn,zn,  da11,da12,da21,da22, dpn, drn,dzn
      real*8 dRgrid, dZgrid, Rgrid_mn,Rgrid_mx, Zgrid_mn,Zgrid_mx
      integer ir,ir1,ir2, iz,iz1,iz2
      
c-----external
      real*8  tempe_xyz,dense_xyz,zeff,psi_rho,bxyz,efficien,rhov,rhos,
     &qsafety_psi,bmin_psi,bmax_psi,dvol_dpsi,rho_lrho,u_res,
     &b_average_xyz, zeffi, zeffrho, wcw, wpw_2, psif_xyz

      pi=4.d0*datan(1.d0)

c     wx,wy and wz are in (cm) 
      x_m=wx(is)*0.01d0/r0x ! normalization
      y_m=wy(is)*0.01d0/r0x ! normalization
      z_m=wz(is)*0.01d0/r0x !     ! Local Z of ray element 'is'
      r_m= dsqrt(x_m**2 + y_m**2) ! Local R of ray element 'is'
      bmod=bxyz(x_m,y_m,z_m) !-> get b (needed for wcw)
      
      ! Make sure rhobegin,rhoend are not exceeding rho_bin(NR).
      ! In such a case, the ray is supposed to be stopped, 
      ! but just in case:
      rhobegin= min(rhobegin,rho_bin(NR))
      rhoend=   min(rhoend,rho_bin(NR))

c  radius rho is in common/one/, calculated in subroutine: dense_xyz
c  rho is used in functions:z_eff
c  The magnetic field bmode is in common/one/ .bmode is used
c  in function wcw(x,y,z)

c begin if is=1      

      if(is.eq.1) then

c        calculation of efficiency on the first step
         temp=tempe_xyz(x_m,y_m,z_m,1)
         ye=wcw(x_m,y_m,z_m,1)
         u1=u_res(jwave,cnpar,temp,ye)
         den=dense_xyz(x_m,y_m,z_m,1) ! rho is calc. inside
         z_eff=zeffrho(rho) !
         
         if (rho.ge.1.d0) then ! no CD
c---------------------------------------------------------------
c          zero CD efficiency outside the plasma
c----------------------------------------------------------------
           eff(is)=0.d0
           goto 100
         endif
 
         
         if (ieffic.eq.1 .and. model_b.eq.0) then
c -------------------------------------------------------------------
c     calculation of current drive efficiency using asymptotic formula
c--------------------------------------------------------------------
            eff(is)=efficien(z_eff,u1,jwave,temp,den) 
         endif

         if (ieffic.eq.2 .and. model_b.eq.0) then
c -------------------------------------------------------------------
c     calculation of current drive efficiency using Ehst-Karney formula
c--------------------------------------------------------------------
            call efKarney_xyz(z_m*100.d0*r0x,r_m*100.0d0*r0x,r0x,
     +                 z_eff,temp,den,jwave,
     1                 cnpar,eff(is))
c----------------------------------------------------
c          ef_Karney using Bonoli subroutine:
c          For the FW NSTX case the Bonoli subroutine gave efficiency
c          which coincided exactly with the efficiency obtained 
c          calculated bu subroutine call efKarney
c-------------------------------------------------------
c          call efKarney_Bonoli(z_m*100.d0*r0x,r_m*100.0d0*r0x,r0x,
c     +                 z_eff,temp,den,jwave,
c     1                 cnpar,eff(is))
c------------------------------------------------------------
         endif   ! ieffic.eq.2

         if (ieffic.eq.3 .and. model_b.eq.0) then
c -------------------------------------------------------------------
c          calculation of current drive efficiency using curba
c--------------------------------------------------------------------
           x_r=(wx(is))
           y_r=(wy(is))
           z_r=(wz(is))
           r_r= dsqrt(x_r**2 + y_r**2)
           zmar=(zma*100.d0*r0x) !cm
           rmar=(rma*100.d0*r0x) !cm
           z_effr=(z_eff)
           tempr=(temp)
           denr=(den)
           cnparr=(cnpar)
           yer=(ye)
           call effcurb_xyz(x_r,y_r,z_r,zmar,rmar,r0x,z_effr,tempr,denr,
     +            jwave, cnparr, yer, effic_r)
           write(*,*)'in p_c_prof after effcurb efffic_r',effic_r
           eff(is)=(effic_r)
         endif     ! ieffic.eq.3

         if (ieffic.eq.4 .and. model_b.eq.0) then
c -------------------------------------------------------------------
c          calculation of current drive efficiency using Lin_liu
c          TorGA_curgap
c--------------------------------------------------------------------
           x_r=(wx(is))
           y_r=(wy(is))
           z_r=(wz(is))
           r_r= dsqrt(x_r**2 + y_r**2)
           zmar=(zma*100.d0*r0x) !cm
           rmar=(rma*100.d0*r0x) !cm
           z_effr=(z_eff)
           tempr=(temp)
           denr=(den)
           cnparr=(cnpar)
           yer=(ye)
           call eff_Lin_Liu_xyz(x_r,y_r,z_r,zmar,rmar,r0x,z_effr,
     +                 tempr,denr,
     &                 jwave,
     &                 cnparr,
     &                 cnper,ioxm,ye,
     &                 cefldx,cefldy,cefldz,
     &                 effic_r)
           write(*,*)'in p_c_prof after eff_Lin_Liu efffic_r',effic_r
           eff(is)=(effic_r)
         endif    ! ieffic.eq.4

c--------------------------------------------------------------------
c        toroidal and poloidal current drives efficiencies for is=1
 100     continue
         bmod=bxyz(x_m,y_m,z_m) ! for what?
c--------------------------------------------------------------------
         allpower=0.0d0
         allpw_e=0.0d0
         allpw_i=0.0d0
         allpw_cl=0.0d0
         allcur=0.0d0
c------- initialization arrays
c        for power and current
         do i=1,NR
            power(i)=0.0d0
            current(i)=0.0d0           
            power_e(i)=0.0d0
            power_i(i)=0.0d0
            power_cl(i)=0.0d0
            cur_den_parallel(i)=0.d0 
         enddo

         !YuP if (iabsorp.eq.3) then ! power_s set to 0.
            do kk=1,nbulk ! YuP was: 2,nbulk
               do i=1,NR
                  power_s(i,kk)=0.0d0
               enddo
            enddo
         !YuP endif
         
         pwr_rz_e=0.d0  ! initialize at is=1
         pwr_rz_i=0.d0
         pwr_rz_cl=0.d0

	 return !goto 30
      endif
c end if is=1

c-------------------------------------------------------------------------
c-------------------------------------------------------------------------

c-----------------------------------------------------------
c     here delpower and allpower are in (erg/sec)
c------------------------------------------------------------
      wss= ws(is)-ws(is-1) ! dl length along ray
      
            
      if(i_call_prep3d_in_output_at_i_ox_conversion_eq_1.eq.0)then
        ! Normal change of power along ray (power depos. to e,i,cl)
        delpower= delpwr(is-1)-delpwr(is)
        delpow_i= .5*wss*(delpwr(is-1)*sdpwr(is-1)+delpwr(is)*sdpwr(is))
        !YuP if (iabsorp.eq.3) then ! delpow_s
          do kk=2,nbulk
            delpow_s(kk)= 0.5d0*wss*
     1         (delpwr(is-1)*salphas(is-1,kk)+delpwr(is)*salphas(is,kk))
          enddo
        !endif
        delpow_e= .5*wss*
     1            (delpwr(is-1)*salphal(is-1)+delpwr(is)*salphal(is))
        delpow_cl=.5*wss*
     1            (delpwr(is-1)*salphac(is-1)+delpwr(is)*salphac(is))

        del= delpow_i +delpow_e +delpow_cl
        !if (del.gt.1.d-100) then ! YuP [Nov-2014] Changed to:
        if (del.gt.1.d-100*delpwr(1)) then ! compare dP with P(t=0)
           delpower_del= delpower/del ! For correction of depos.powers
        else ! del~0 almost no absorption
           goto 30 !-> exit: nothing to do
        endif
      elseif(i_call_prep3d_in_output_at_i_ox_conversion_eq_1.eq.2)then
        ! The previous time step (is-1) was the last point of O-mode
        ! and this step is the first point of X-mode.
        ! There is generally a big change in delpwr from is-1 to is
        ! because of transm_ox coefficient.
        ! Do not include this jump into calculation of power to e,i,cl 
        delpow_i= wss*delpwr(is)*sdpwr(is)
        !YuP if (iabsorp.eq.3) then ! delpow_s
          do kk=2,nbulk
            delpow_s(kk)= wss*delpwr(is)*salphas(is,kk)
          enddo
        !endif
        delpow_e=  wss*delpwr(is)*salphal(is)
        delpow_cl= wss*delpwr(is)*salphac(is)
        del= delpow_i+delpow_e+delpow_cl !Not equal to delpower in this case
        delpower= (delpwr(is-1)-delpwr(is)) ! Includes the O->X loss
        delpower_del=1.d0 ! No correction in powers to e,i,cl
      else
        ! Do nothing: this case is not possible. But just in case:
        del=0.d0
        delpower= del
        delpower_del=1.d0
        goto 30 !-> exit: nothing to do
      endif
                 
      ! Adjust powers 
      !(Normally, delpower=del, but sometimes they can differ 
      ! in 2nd digit if absorption is strong, as in case of EBW)
      delpow_e=delpow_e*delpower_del
      delpow_i=delpow_i*delpower_del
      
      if(delpow_e.lt.0.d0 .or. 
     +   delpow_i.lt.0.d0 .or. delpwr(is).lt.0.d0)then
      write(*,*)
     +'p_c_prof_xyz: is, delpow_e, delpow_i, delpower_del,
     + delpwr(is), sdpwr(is)=',
     +               is, delpow_e, delpow_i, delpower_del,
     + delpwr(is), sdpwr(is)
      !pause
      endif
      
      !YuP if (iabsorp.eq.3) then
         do kk=1,nbulk ! YuP was: 2,nbulk ! delpow_s
            delpow_s(kk)=delpow_s(kk)*delpower_del
         enddo
      !YuP endif
      delpow_cl=delpow_cl*delpower_del
      
      allpower= allpower + (delpwr(is-1)-delpwr(is))
      !allpower sums-up any change of power along ray, including OX-loss
      !(power "lost" at O->X mode transmission, which is 
      ! technically not a lost power but a reflected O-mode power)
      allpw_e=allpw_e+delpow_e
      allpw_i=allpw_i+delpow_i
      allpw_cl=allpw_cl+delpow_cl

c      if(abs(del-delpower).gt.1.d-2*del .and. delpower.gt.1.d-4)then
c        write(*,'(a,i5,2e13.4)')'is,delpwr(is-1),delpwr(is)=',
c     +                           is,delpwr(is-1),delpwr(is)
c        write(*,'(a,3e13.4)')'r_m, del, delpower=', r_m, del, delpower
c        pause !!!
c      endif
    
c-----------------------------------------------------------


      ! YuP[Nov-2014] Added: Saving Local absorbed power over (R,Z) grid.
      ! The power absorbed in a ray element from (is-1) to (is) point
      ! is attributed to four corners of the grid cell
      ! which contains the point (rmid,zmid) of the ray element.
      ! This is done by bilinear weighting, similar to PIC methods.
      !wx,wy,wz are in (cm) 
      x0=wx(is-1)*0.01d0/r0x ! normalization
      y0=wy(is-1)*0.01d0/r0x ! normalization
      z0=wz(is-1)*0.01d0/r0x ! ! Local Z of ray point 'is-1'
      r0= dsqrt(x0**2 + y0**2) ! Local R of ray point 'is-1'
      rmid=0.5*(r_m+r0) ! R mid-point between (is-1) and (is) of the ray
      zmid=0.5*(z_m+z0) ! Z mid-point between (is-1) and (is) of the ray
      ! This is from the grid definition:
      Rgrid_mn= Rgrid(1)
      Rgrid_mx= Rgrid(NRgrid)
      Zgrid_mn= Zgrid(1)
      Zgrid_mx= Zgrid(NZgrid)
      dRgrid= Rgrid(2)-Rgrid(1)   ! because of uniform grid in R
      dZgrid= Zgrid(2)-Zgrid(1)   ! because of uniform grid in Z
      ! Make sure (rmid,zmid) point is within (Rgrid,Zgrid):
      rmid= max(rmid, Rgrid(1) ) 
      rmid= min(rmid, Rgrid_mx-dRgrid*1.d-3) !stay away from upper edge
      zmid= max(zmid, Zgrid(1) ) 
      zmid= min(zmid, Zgrid_mx-dZgrid*1.d-3) !stay away from upper edge
      ! Find the cell in (Rgrid,Zgrid) that contains (rmid,zmid):
      rn= (rmid-Rgrid_mn)/dRgrid  ! normalized R coordinate
      zn= (zmid-Zgrid_mn)/dZgrid  ! normalized Z coordinate
      ir=  rn    ! lower integer: can be from 0 to NRgrid-2
      ir1= ir+1
      ir2= ir+2  ! rn is between ir1 and ir2 nodes of Rgrid
      iz=  zn    ! lower integer: can be from 0 to NZgrid-2
      iz1= iz+1
      iz2= iz+2  ! zn is between iz1 and iz2 nodes of Zgrid
      ! In normalized units:
      ! distance between (rmid,zmid) point and left-bottom-corner 
      drn= rn-ir  
      dzn= zn-iz 
      ! Four areas within the cell marked by a given point (rn,zn)
      !    _____
      !  Z| |   |     Left-top area: da12 ,   Right-top:    da22
      !   |-o---|
      !   |_|___|     Left-bottom:   da11 ,   Right-bottom: da21
      !        R
      ! Note that in normalized units, the area of the grid cell is 1.0
      da11= drn*dzn
      da12= drn -da11 ! == drn*(1-dzn)
      da21= dzn -da11 ! == dzn*(1-drn)
      da22= 1.d0 -dzn -da12 ! == 1-dzn -drn*(1-dzn) = (1-drn)*(1-dzn)
      ! Note that da11+da12+da21+da22 = 1.0
      ! The above form uses only one multiplication !
      ! Attribute the weighting factor proportional to the area 
      ! to the opposite corner:
      !- e -
      dpn= delpow_e !absorbed power in the grid cell
      pwr_rz_e(ir1,iz1)= pwr_rz_e(ir1,iz1) +da22*dpn
      pwr_rz_e(ir1,iz2)= pwr_rz_e(ir1,iz2) +da21*dpn
      pwr_rz_e(ir2,iz2)= pwr_rz_e(ir2,iz2) +da11*dpn
      pwr_rz_e(ir2,iz1)= pwr_rz_e(ir2,iz1) +da12*dpn
      !- i -
      dpn= delpow_i !absorbed power in the grid cell
      pwr_rz_i(ir1,iz1)= pwr_rz_i(ir1,iz1) +da22*dpn
      pwr_rz_i(ir1,iz2)= pwr_rz_i(ir1,iz2) +da21*dpn
      pwr_rz_i(ir2,iz2)= pwr_rz_i(ir2,iz2) +da11*dpn
      pwr_rz_i(ir2,iz1)= pwr_rz_i(ir2,iz1) +da12*dpn
      !- cl -
      dpn= delpow_cl !absorbed power in the grid cell
      pwr_rz_cl(ir1,iz1)= pwr_rz_cl(ir1,iz1) +da22*dpn
      pwr_rz_cl(ir1,iz2)= pwr_rz_cl(ir1,iz2) +da21*dpn
      pwr_rz_cl(ir2,iz2)= pwr_rz_cl(ir2,iz2) +da11*dpn
      pwr_rz_cl(ir2,iz1)= pwr_rz_cl(ir2,iz1) +da12*dpn
      !( pwr_rz_s() could be added here. Not needed, for now)
      ! YuP[Nov-2014] Done: Saving Local absorbed power over (R,Z) grid.


      if(rhoend.gt.rhobegin) then
         rhomax=rhoend
         rhomin=rhobegin
      else
         rhomin=rhoend
         rhomax=rhobegin
      endif
          
cyup YuP[10-03-2014] Allow for rho>1, for any model_b
      hrho= rho_bin(2)-rho_bin(1) ! YuP: was: hrho=1.d0/(NR-1)
 
      if (rhobegin.ge.rho_bin(NR) .and. rhoend.ge.rho_bin(NR) ) then 
         !Special case: 
         !  rhobegin and rhoend are at the last point of rho_bin()
         !  (or outside of rho_bin grid)
         power(NR)= power(NR)+delpower 
cyup         cur_den_parallel(NR)=cur_den_parallel(NR)+
cyup     &     delpow_e*0.5d0*(eff_rho_min+eff_rho_max)/binvol(NR)
cyup         current(NR)=current(NR)+delcurr   !toroidal current from bin
         power_e(NR)= power_e(NR)+delpow_e
         power_i(NR)= power_i(NR)+delpow_i
            do kk=1,nbulk ! YuP was: 2,nbulk
               power_s(NR,kk)= power_s(NR,kk)+delpow_s(kk)
            enddo
         power_cl(NR)= power_cl(NR)+delpow_cl     
         goto 30
      endif !  
    
      ! Find jbinmin and jbinmax such that rhomin is in the bin #jbinmin
      ! and rhomax is in the bin #jbinmax
      jbinmin=1 ! Initial guess
      jbinmax=NR-1
      do j=1,NR-1
         if(rhomin.lt.rho_bin(j+1)) then
           jbinmin=j
           goto 10
         endif
      enddo
10    continue
      do j=jbinmin,NR-1
         if(rhomax.lt.rho_bin(j+1)) then
            jbinmax=j
            goto 20
         endif
      enddo
20    continue

      if (jbinmin.eq.NR) jbinmin=NR-1
      
c     r0 (cm)
      r0=0.5d0*(wr(is)+wr(is-1)) ! for CD calculations

c---- calculation of the efficiency 
      temp=tempe_xyz(x_m,y_m,z_m,1)
      ye=wcw(x_m,y_m,z_m,1)
      u1=u_res(jwave,cnpar,temp,ye)
      den=dense_xyz(x_m,y_m,z_m,1) ! rho is calc. inside
      z_eff=zeffrho(rho) 
      if (rho.ge.1.d0) then ! no CD
c---------------------------------------------------------------
c        zero CD efficiency outside the plasma
c----------------------------------------------------------------
         eff(is)=0.d0
         goto 110
      endif
   
      if (ieffic.eq.1 .and. model_b.eq.0) then
c -------------------------------------------------------------------
c       calculation of current drive efficiency using asymptotic formula
c--------------------------------------------------------------------
        eff(is)=efficien(z_eff,u1,jwave,temp,den)
      endif

      if (ieffic.eq.2 .and. model_b.eq.0) then
c -------------------------------------------------------------------
c       calculation of current drive efficiency using Ehst-Karney formula
c--------------------------------------------------------------------
      call efKarney_xyz(z_m*100.0d0*r0x,r_m*100.0d0*r0x,r0x,
     & z_eff,temp,den,jwave,cnpar,
     &              eff(is))
c----------------------------------------------------
c          ef_Karney using Bonoli subroutine:
c          For the FW NSTX case the Bonoli subroutine gave efficiency
c          which coinsided exectly wiht the efficiency obtained 
c          calculated bu subroutine call efKarney
c-------------------------------------------------------
c      call efKarney_Bonoli(z_m*100.0d0*r0x,r_m*100.0d0*r0x,r0x,
c     & z_eff,temp,den,jwave,cnpar,
c     &              eff(is))
      endif  !ieffic.eq.2

      if (ieffic.eq.3 .and. model_b.eq.0) then
c -------------------------------------------------------------------
c       calculation of current drive efficiency using curba
c-------------------------------------------------------------------
           x_r=(wx(is))
           y_r=(wy(is))
           z_r=(wz(is))
        r_r= dsqrt(x_r**2 + y_r**2)
        zmar=(zma*100.d0*r0x)
        rmar=(rma*100.d0*r0x)
        z_effr=(z_eff)
        tempr=(temp)
        denr=(den)
        cnparr=(cnpar)
        yer=real(ye)

ctest
cyup        u1=u_res(jwave,cnparr,tempr,yer)
cyup        write(*,*)'jwave,cnpar,tempr,yer,u1',jwave,cnparr,tempr,yer,u1
cyup        effic_r=efficien(z_effr,u1,jwave,tempr,denr)
cyup        write(*,*)'asimptotic: z_effr,denr,effic_r',z_effr,denr,effic_r
cendtest
        call effcurb_xyz(x_r,y_r,z_r,zmar,rmar,r0x,z_effr,tempr,denr,
     +  jwave,cnparr,yer,effic_r)
c        write(*,*)'in p_c_prof after effcurb is,effic_r',is,effic_r
        eff(is)=effic_r
      endif !ieffic.eq.3

      if (ieffic.eq.4 .and. model_b.eq.0) then
c -------------------------------------------------------------------
c       calculation of current drive efficiency using eff_Lin_Liu
c-------------------------------------------------------------------
           x_r=(wx(is))
           y_r=(wy(is))
           z_r=(wz(is))
        r_r= dsqrt(x_r**2 + y_r**2)
        zmar=(zma*100.d0*r0x)
        rmar=(rma*100.d0*r0x)
        z_effr=(z_eff)
        tempr=(temp)
        denr=(den)
        cnparr=(cnpar)
        yer=real(ye)

        call effcurb_xyz(x_r,y_r,z_r,zmar,rmar,r0x,z_effr,tempr,denr,
     +  jwave,cnparr,yer,effic_r)

        write(*,*)'in p_c_prof after effcurb is,effic_r',is,effic_r

           call eff_Lin_Liu_xyz(x_r,y_r,z_r,zmar,rmar,r0x,z_effr,
     +                 tempr,denr,
     &                 jwave,
     &                 cnparr,
     &                 cnper,ioxm,ye,
     &                 cefldx,cefldy,cefldz,
     &                 effic_r)
        write(*,*)'in p_c_prof after eff_Lin_Liu efffic_r',effic_r
   
        eff(is)=effic_r
      endif !4


 110  continue

c--------------------------------------------------------------------
c     delpower(erg/sec),delcurr(Ampere),r0(cm)
c-----------------------------------------------------------

c-----calculate parallel CD using the geometric factor 1/(2pi*r0) 
c      delcurr=delpow_e*0.5d0*(eff(is-1)+eff(is))/(2*pi*r0)
c      write(*,*)'1/(2*pi*r0)',1.d0/(2.d0*pi*r0)

c------geometric factor 1/(r*pi*R_mag_axis)
c      zfacgeom = 1.0d0 / (2.0d0 * pi * rma)/100.d0
c      write(*,*)'1/(2*pi*rma)/100.d0',1.d0/(2.d0*pi*rma)/100.d0

c-----calculate toroidal CD
      x0=0.5d0*(wx(is)+wx(is-1))
      y0=0.5d0*(wy(is)+wy(is-1))
      z0=0.5d0*(wz(is)+wz(is-1))
      r0m=r0*1.d-2
      x0m=x0*1.d-2   
      y0m=y0*1.d-2   
      z0m=z0*1.d-2   
c-----geometric factor ~ 1/b_averaged 
      if(model_b.ge.1) then
         bavr=bxyz(x0m,y0m,z0m) ! needs work
      else ! (model_b.eq.0)
         psiloc=psif_xyz(x0m,y0m,z0m) ! YuP Added instead of psi_rho
         bavr=b_average_xyz(psiloc)
         zfacgeom = -2.d0 * pi * qsafety_psi(psiloc)/
     &     (dvol_dpsi(psiloc)*dpsimax*bavr)/100.d0
         delcurr=delpow_e*0.5d0*(eff(is-1)+eff(is))*zfacgeom !toroidal current
      endif

                                                          !created by delpow_e
c
cSAP080902 to plot delta power and delta current along the ray
      delpow_e_ar(is)=delpow_e
      delcur_par_ar(is)=delpow_e*0.5d0*(eff(is-1)+eff(is))
      
c-----toroidal and poloidal CD from old genray version      
cyup      rho0=0.5*(spsi(is)+spsi(is-1)) ! the small radius
cyup      rho0_pol=rho_lrho(rho0)
cyup      poloidlen=rho0_pol*totlength*100.d0     ! poloidal length cm
            
      allcur=allcur+delcurr                     !total toroidal current
      
      if(rhoend.gt.rhobegin) then
          eff_rho_max=eff(is)
          eff_rho_min=eff(is-1)
      else
          eff_rho_max=eff(is-1)
          eff_rho_min=eff(is)
      endif

      if (jbinmin.eq.jbinmax) then ! rhobegin and rhoend are within same bin
         power(jbinmin)=power(jbinmin)+delpower 
cyup         cur_den_parallel(jbinmin)=cur_den_parallel(jbinmin)+
cyup     &     delpow_e*0.5d0*(eff_rho_min+eff_rho_max)/binvol(jbinmin)
cyup         current(jbinmin)=current(jbinmin)+delcurr   !toroidal current from bin
         power_e(jbinmin)=power_e(jbinmin)+delpow_e
         power_i(jbinmin)=power_i(jbinmin)+delpow_i
            do kk=1,nbulk ! YuP was: 2,nbulk
               power_s(jbinmin,kk)=power_s(jbinmin,kk)+delpow_s(kk)
            enddo
         power_cl(jbinmin)=power_cl(jbinmin)+delpow_cl     
         goto 30
      endif ! if (jbinmin.eq.jbinmax) 
      delrho=rhomax-rhomin
      ppow=delpower/delrho 
      pcur=delcurr/delrho
      ppow_e=delpow_e/delrho
      ppow_i=delpow_i/delrho
      !YuP if (iabsorp.eq.3) then
         do kk=1,nbulk ! YuP was: 2,nbulk ! ppow_s
            ppow_s(kk)=delpow_s(kk)/delrho
         enddo
      !YuP endif
      ppow_cl=delpow_cl/delrho

c------------------------------------------------------------
c     power (erg/sec), current(A)
c------------------------------------------------------------

      power(jbinmin)=power(jbinmin)+ppow*(rho_bin(jbinmin+1)-rhomin)
      power_e(jbinmin)=  power_e(jbinmin)
     +                  +ppow_e*(rho_bin(jbinmin+1)-rhomin)
      power_i(jbinmin)=  power_i(jbinmin)
     +                  +ppow_i*(rho_bin(jbinmin+1)-rhomin)
      power_cl(jbinmin)= power_cl(jbinmin)
     +                  +ppow_cl*(rho_bin(jbinmin+1)-rhomin)
      !YuP if (iabsorp.eq.3) then ! power_s
         do kk=1,nbulk ! YuP was: 2,nbulk
            power_s(jbinmin,kk)= power_s(jbinmin,kk)
     &                          +ppow_s(kk)*(rho_bin(jbinmin+1)-rhomin)
         enddo
      !YuP endif
      
cyup      cur_den_parallel(jbinmin)=cur_den_parallel(jbinmin)+
cyup     &(delpow_e/delrho)*(rho_bin(jbinmin+1)-rhomin)/binvol(jbinmin)*
cyup     &0.5d0*(eff_rho_min+
cyup     &      (eff_rho_min+(eff_rho_max-eff_rho_min)*
cyup     &                   (rho_bin(jbinmin+1)-rhomin)/(delrho)))
cyup      current(jbinmin)=
cyup     1 current(jbinmin)+pcur*(rho_bin(jbinmin+1)-rhomin)  !toroidal current from bin
 
      power(jbinmax)=power(jbinmax)+ppow*(rhomax-rho_bin(jbinmax))
      power_e(jbinmax)=power_e(jbinmax)+ppow_e*(rhomax-rho_bin(jbinmax))
      power_i(jbinmax)=power_i(jbinmax)+ppow_i*(rhomax-rho_bin(jbinmax))
      power_cl(jbinmax)=power_cl(jbinmax)+
     1                  ppow_cl*(rhomax-rho_bin(jbinmax))
      !YuP if (iabsorp.eq.3) then ! power_s
         do kk=1,nbulk ! YuP was: 2,nbulk
            power_s(jbinmax,kk)=
     &          power_s(jbinmax,kk)+ppow_s(kk)*(rhomax-rho_bin(jbinmax))
         enddo
      !YuP endif

cyup      cur_den_parallel(jbinmax)=cur_den_parallel(jbinmax)+
cyup     &(delpow_e/delrho)*(rhomax-rho_bin(jbinmax))/binvol(jbinmax)*
cyup     &0.5d0*(eff_rho_max+
cyup     &   (eff_rho_min+(eff_rho_max-eff_rho_min)*
cyup     &   (rhomax-rho_bin(jbinmax))/(delrho)))  
cyup      current(jbinmax)=
cyup     1	  current(jbinmax)+pcur*(rhomax-rho_bin(jbinmax)) !toroidal current from bin

      if(jbinmax.gt.(jbinmin+1)) then
         do j=(jbinmin+1),(jbinmax-1)
            power(j)=power(j)+ppow*hrho
            power_e(j)=power_e(j)+ppow_e*hrho
            power_i(j)=power_i(j)+ppow_i*hrho
            power_cl(j)=power_cl(j)+ppow_cl*hrho
            !YuP if (iabsorp.eq.3) then ! power_s
               do kk=1,nbulk ! YuP was: 2,nbulk
                  power_s(j,kk)=power_s(j,kk)+ppow_s(kk)*hrho
               enddo
            !YuP endif
cyup            cur_den_parallel(j)=cur_den_parallel(j)+
cyup     &      (delpow_e/delrho)*hrho/binvol(j)*
cyup     &      (eff_rho_min+(eff_rho_max-eff_rho_min)*
cyup     &      (rho_bin_center(j)-rhomin)/delrho)
cyup            current(j)=current(j)+pcur*hrho !toroidal current from bins
         enddo
      endif

30    continue


      
      end ! p_c_prof_xyz
      
      

c======================================================================
c======================================================================


      subroutine anth_rlt_xyz(wpw2,wccw,T_kev,nll_in,np_in,
     +n_relt_harm1,n_relt_harm2,n_relt_intgr,
     +i_resonance_curve_integration_method,epsi,
     +i_fkin,x,y,z,
     +aK)
c     calculates anti-hermitian relativistic dielectric tensor aK
c     for electron plasma
c
c     INPUTS:
c
c      wpw2 = (fpe/f)**2
c      wccw = fce/f
c      T_kev  - electron temperature
c      nll_in - parallel index of refraction N.
c      np_in  - perpendicular refractive index N
c     n_relt_harm1 min number of EC harmonics used in anti-hermitian
c               dielectric tensor calculations
c     n_relt_harm2 max number of EC harmonics used in anti-hermitian
c               dielectric tensor calculations
c               n_relt_harm1 =< n= <n_relt_harm2
c      n_relt_intgr - the number of points for the integration over p_perp
c      i_fkin =0 the usage of the analytical relativistic distributin
c             =1 the usage of the numerical 3D distribution from diskf or netcdfnm.nc 
c                written be CQL3D code or created analytically at mesh points
c      x,y,z  - cartesian coords.
c-------------------------------------------------------------------
!       i_resonance_curve_integration_method=1 !rectangle integration
!                                              !over angle,
!                                              !for ellipse case only
!       i_resonance_curve_integration_method=2 !rectangle formula,
!                                              !p_perp integration
!                                              !elliptic or hyperbolic case
!       i_resonance_curve_integration_method=3 !trapezoidal formula,
!                                              !p_perp integration
!                                              !elliptic or hyperbolic case
!       i_resonance_curve_integration_method=4 !adaptive Simpson integration 
!                                              !p_perp integration
!                                              !elliptic or hyperbolic case
!
!       i_resonance_curve_integration_method is used in subroutine intgr_rl
!       to choose the numerical integration method for 
!       anti-hermitian relativistic dielectric tensor calculation.
!       This applies for iabsorp=6,7 and for emission calculations.
c--------------------------------------------------------------------
c     OUTPUT:
c      aK(3,3):  the nine components of the anti-hermitian
c                part dielectric relativistic tensor
c                evaluated at (wpw2,wccw,Te,nll,np,n)

      implicit none
c     input 
      double precision wpw2,wccw,T_kev,nll_in,np_in,epsi
      integer n_relt_harm1,n_relt_harm2,n_relt_intgr,i_fkin,
     &i_resonance_curve_integration_method

      double precision r,x,y,z

c     output
      double complex aK(3,3)
c     local
      double precision c,mass_e,k_1_kev,nll,np,nlls,p_perp0,theta,pi,
     .dens,xpidens,rho,psi
      double complex integral(3,3)
      integer jn,ires

c     external zeroK, npnllmin, ec_cond, intgr_rl,fdens_fdist
      double precision fdens_fdist,densrho,rhopsi, dense_xyz
 
      c =2.99792458d10          !speed of light     (cm/sec)
      mass_e=9.1094d-28         !electron rest mass (g)
      k_1_kev=1.6022d-9         !ergs in 1 KeV      (erg)

      theta=mass_e*c**2/(k_1_kev*T_kev)

c      write(*,*)'relt_tens.f sub anth_rlt: r,z,T_kev,theta',
c     &r,z,T_kev,theta
cyup      write(*,*)'relt_tens.f sub anth_rlt:'
cyup      write(*,*)'i_resonance_curve_integration_method',
cyup     & i_resonance_curve_integration_method

      pi=4.d0*datan(1.d0)

c YuP 120430      call npnllmin(nll_in,np_in,nll,np)
      nll=nll_in ! YuP 120430
      np=np_in   ! YuP 120430
      
      nlls = nll**2
 
c-----initialization aK=0
      call zeroK(aK)

cyup      write(*,*)'n_relt_harm1,n_relt_harm2',n_relt_harm1,n_relt_harm2
c-----the loop over the cyclotron harmonics 
      do jn=n_relt_harm1,n_relt_harm2               
c-------control the EC resonance condition
cSm060315 -jn         
        call ec_cond(-jn,wccw,nll,ires,p_perp0)

        if(ires.eq.0) then
c---------no resonace
          goto 10
        endif
cSm060315 -jn      
c        write(*,*)'relt_tens.f sub anth_rlt jn=',jn
cyup        write(*,*)'bef intgr_rl i_resonance_curve_integration_method',
cyup     &  i_resonance_curve_integration_method
   
        call intgr_rl_xyz(-jn,nll,np,wccw,theta,ires,p_perp0,
     &  n_relt_intgr,i_resonance_curve_integration_method,epsi,
     +  i_fkin, x,y,z,
     +  integral)

c        write(*,*)'after intgr_rl i_resonance_curve_integration_method',
c     &  i_resonance_curve_integration_method
   
c        write(*,*)'anth_rl,integral',integral
       
        aK(1,1)=aK(1,1)+integral(1,1)
        aK(1,2)=aK(1,2)+integral(1,2)
        aK(1,3)=aK(1,3)+integral(1,3)
        aK(2,2)=aK(2,2)+integral(2,2)
        aK(2,3)=aK(2,3)+integral(2,3)
        aK(3,3)=aK(3,3)+integral(3,3)

10      continue      
      enddo   !jn

      if (i_fkin.eq.0) then       
        dens=1.d0
      else      
        dens=dense_xyz(x,y,z,1) ! YuP added
      endif
      
c-----normalization for the unit electron density
      xpidens=wpw2*pi/dens
      
      aK(1,1)=-xpidens*aK(1,1)
      aK(1,2)=-xpidens*aK(1,2)
      aK(1,3)=-xpidens*aK(1,3)
      aK(2,2)=-xpidens*aK(2,2)
      aK(2,3)=-xpidens*aK(2,3)
      aK(3,3)=-xpidens*aK(3,3)
      
      aK(2,1)=-aK(1,2)
      aK(3,1)= aK(1,3)
      aK(3,2)=-aK(2,3)

      return
      end
      
c======================================================================
c======================================================================



      subroutine intgr_rl_xyz(n,nll,np,wccw,theta,ires,p_perp0,
     +n_relt_intgr,i_resonance_curve_integration_method,epsi,
     +i_fkin, x,y,z,
     +integral)
c-----------------------------------------------------------------
c     calculates the matrix: double complex integral
c     for the relativistic electron plasma 
c     I(n)^=integral(0<=p_perp<=p_perp0)G_n(p_perp)
c     G_n={sum_k=1,2,G_nk(p_perp)
c     for the EC harmonic with number 'n'
c----------------------------------------------------------------- 
c     input
c       n        - EC harmonic number
c       nll      - N_parallel
c       np       - N_perpendicular
c       wccw     = omega_ce/omega for the electron rest mass
c       theta    = mc**2/T
c       ires     =0 no resonance,=1 ellipse,=2 parabola,=3 hyperbole 
c       n_relt_intgr - the number of points for the integration over p_perp
c       p_perp0  - max value of the perpendicular momentum divided by mc
c                  on the resonanse ellipse
c       i_fkin   =0 the usage of the analytical relativistic distributin
c                =1 the usage of the numerical 3D distribution from diskf 
c                   written be CQL3D code
c       x,y,z   - cartesian coords.
c----------------------------------------------------------
c     i_resonance_curve_interation_method
c
c     i_resonance_curve_integration_method=1 angle integration
c     now it works for ellipse case only                             c
c
c     i_resonance_curve_integration_method=2 rectangle formula
c     for p_perp integration
c     
c     i_resonance_curve_integration_method=3 trapezoidal formula
c     for p_perp integration
c
c     i_resonance_curve_integration_method=4 !adaptive Simpson 
c     for p_perp integration
c
c       epsi is the accuracy used in integration by adaptive Simpson
c------------------------------------------------------------ 
c     
c     output
c       integral(3,3) double complex integral from G_n over 0<=p_perp=<p_perpmax      
c----------------------------------------------------------------    
c      The integration method is specified by the variable
c      i_resonance_curve_integration_method.
c      Now this variable is set inside this subroutine.
c----------------------------------------------------------------

      implicit none
c     input
      double precision nll,np,wccw,theta,p_perp0,epsi
      integer n,ires,n_relt_intgr,i_fkin,
     &i_resonance_curve_integration_method

      double precision r,x,y,z
      double precision vnormloc,massloc ! cm/sec, g
      COMMON /dskin1/vnormloc,massloc

c     output
      double complex integral(3,3)
c     local
      double precision eps,p_permax,h,p_perp,p,p_t,clight
      double complex i,g_n(3,3)
cSm060725
      integer j,jmax
      double precision vmax_d_vt 
ctest begin`
      double precision t_intgr,p_int
ctest end
c-----for integration along ellipse by angle
      double precision rme,rtem0,vper,xint,cs,sn,thet1,thet2,p_par,
     & p_par_min,p_par_pl,p_perp_min,p_perp_pl,v0dc,vmax1dc,vpar0dc,
     +vpar1dc,vpar2dc,dth,vt0,cdvt,tem0,dvper,vmaxdt
     
c     external
c     ec_cond, zeroK, g_n 
      double precision temperho
 
c      write(*,*)'relt_tens.f sub intgr_rl:'
c      write(*,*)'i_resonance_curve_integration_method',
c     & i_resonance_curve_integration_method

      i = ( 0.0d0,1.0d0)        !imaginary number     
      jmax=n_relt_intgr
      vmax_d_vt=10.d0

      call p_perp_max_calc(i_fkin,theta,n,wccw,nll,vmax_d_vt,
     &vnormloc,p_permax,ires)

      if(ires.eq.4) then
c-------the resonance curve is outside the grid
        call zeroK(integral)
        goto 10
      else
         h=p_permax/(n_relt_intgr)   ! step of integration over p_perp
      endif

cyup      write(*,*)'in intgr_rl h',h

c-----calculations of the integrals over p_perp
      call zeroK(integral)
cyup      write(*,*)'in intgr_rl set integral=zero integral',integral

c----------------------------------------------------------
c     i_resonance_curve_interation_method
c
c     i_resonance_curve_integration_method=1 angle integration
c     now it works for ellipse case only                             c
c     i_resonance_curve_integration_method=2 rectangle formula
c     for p_perp integration
c     
c     i_resonance_curve_integration_method=3 trapezoidal formula
c     for p_perp integration
c     
c------------------------------------------------------------ 
c      i_resonance_curve_integration_method=1
c      i_resonance_curve_integration_method=2
c      i_resonance_curve_integration_method=3 !trapezoidal 
c      i_resonance_curve_integration_method=4 !adaptive Simpson 
     
      goto (1,2,3,4) i_resonance_curve_integration_method

 1    continue
cSm060327       
      if(dabs(nll).ge.1.d0) goto 3 !to trapezoidal formula 
c--------------------------------------------------------------------
c     integration along the resonance curve using angle (thet)
c     thet1<thet<thet2  
c     begin
c     This case now works for ellipce case only
c     |N_parallel| <1
c--------------------------------------------------------------------
      rme=9.1094d-28 
ccc      call get_rtem0_from_one(rtem0)
      tem0=temperho(0.d0,1)*rtem0
      vt0=dsqrt(2.d0*tem0*1.6022d-9/rme) !(cm/sec)the  central
                                         ! thermal velocity
      clight=2.99792458d10
      cdvt=clight/vt0
      vmaxdt=dsqrt(rtem0)
ccc      call ec_condh(n,wccw,nll,vmaxdt/cdvt,ires,v0dc,vmax1dc,vpar0dc,
ccc     +vpar1dc,vpar2dc,thet1,thet2)
     
      call zeroK(integral)     
c-----Subdivide theta range of integration
      dth=(thet2-thet1)/(jmax-1)
      do j=1,jmax
         xint=thet1+(j-1)*dth
         cs=dcos(xint)
         sn=dsin(xint)
         p_par=(vpar0dc-v0dc*cs)       !vper/c
         p_perp=vmax1dc*sn             !vpar/c
         if(p_perp.lt.1.d-12) p_perp=1.d-3
        
ccc         call g_n_calc_theta(p_perp,p_par,y,nll,np,theta,n,i_fkin,
ccc     &   r,z,phi,g_n)
         p_int=1.d0
        
         if((j.ne.jmax).and.(j.ne.1)) then         
           sn=dsin(xint+dth)
           p_perp_pl=vmax1dc*sn             !vper/c
           sn=dsin(xint-dth)
           p_perp_min=vmax1dc*sn            !vper/c  
           h=0.5d0*(p_perp_pl-p_perp_min)
         else
           if(j.eq.1) then
             sn=dsin(xint+dth)
             p_perp_pl=vmax1dc*sn             !vper/c
             h=0.5d0*(p_perp_pl-p_perp)
           else
             !j=jmax
             sn=dsin(xint-dth)
             p_perp_min=vmax1dc*sn             !vper/c
             h=0.5d0*(p_perp-p_perp_min)
           endif
         endif
                 
         integral(1,1)=integral(1,1)+h*g_n(1,1)
         integral(1,2)=integral(1,2)+h*g_n(1,2)
         integral(1,3)=integral(1,3)+h*g_n(1,3)
         integral(2,2)=integral(2,2)+h*g_n(2,2)
         integral(2,3)=integral(2,3)+h*g_n(2,3)
         integral(3,3)=integral(3,3)+h*g_n(3,3)
      enddo 
      goto 20
c--------------------------------------------------------------------
c     integration along the resonance curve using angle (thet)
c     thet1<thet<thet2  
c     end
c--------------------------------------------------------------------

 2    continue
c-------------------------------------------------------------------
c     integration along the resonance curve using rectangle
c     formula
c     begin
c------------------------------------------------------------------- 
      do j=1,jmax
        p_perp=h*(j-0.5)

        call g_n_calc_xyz(p_perp,wccw,nll,np,theta,n,i_fkin,x,y,z,g_n)
        
        integral(1,1)=integral(1,1)+g_n(1,1)
        integral(1,2)=integral(1,2)+g_n(1,2)
        integral(1,3)=integral(1,3)+g_n(1,3)
        integral(2,2)=integral(2,2)+g_n(2,2)
        integral(2,3)=integral(2,3)+g_n(2,3)
        integral(3,3)=integral(3,3)+g_n(3,3)
      enddo 

      integral(1,1)=integral(1,1)*h
      integral(1,2)=integral(1,2)*h
      integral(1,3)=integral(1,3)*h
      integral(2,2)=integral(2,2)*h
      integral(2,3)=integral(2,3)*h
      integral(3,3)=integral(3,3)*h  
      goto 20
c ------------------------------------------------------------------
c     integration along the resonance curve using rectangle
c     formula
c     end
c-------------------------------------------------------------------

 3    continue
cSm060306
c -------------------------------------------------------------------
c     integration along the resonance curve using trapezoidal
c     formula
c     begin
c--------------------------------------------------------------------
      do j=0,jmax
         p_int=1.d0
         if((j.eq.1).or.(j.eq.jmax)) p_int=0.5d0
         p_perp=h*j
         if(p_perp.lt.1.d-3) p_perp=1.d-3
         if(j.eq.jmax) p_perp=p_perp-1.d-3
         
         call g_n_calc_xyz(p_perp,wccw,nll,np,theta,n,i_fkin,x,y,z,g_n)
                  
         integral(1,1)=integral(1,1)+p_int*g_n(1,1)
         integral(1,2)=integral(1,2)+p_int*g_n(1,2)
         integral(1,3)=integral(1,3)+p_int*g_n(1,3)
         integral(2,2)=integral(2,2)+p_int*g_n(2,2)
         integral(2,3)=integral(2,3)+p_int*g_n(2,3)
         integral(3,3)=integral(3,3)+p_int*g_n(3,3)
      enddo 

      integral(1,1)=integral(1,1)*h
      integral(1,2)=integral(1,2)*h
      integral(1,3)=integral(1,3)*h
      integral(2,2)=integral(2,2)*h
      integral(2,3)=integral(2,3)*h
      integral(3,3)=integral(3,3)*h

      goto 20
c-------------------------------------------------------------------
c     end integration along the resonance curve using trapezoidal
c     formula
c     end
c-------------------------------------------------------------------

 4    continue
cSm070417
c -------------------------------------------------------------------
c     integration along the resonance curve using 
c     the adaptive Simpson function
c     begin
c--------------------------------------------------------------------
cyup      write(*,*)'before calc_integral_array_by_adaptive_simpson'

      call calc_absorption_integral_array_by_adaptive_simpson_xyz
     &(wccw,nll,np,theta,n,i_fkin,x,y,z,p_permax,epsi,integral)
cyup      write(*,*)'in intgr_rl after'
cyup      write(*,*)'calc_absorption_integral_array_by_adaptive_simpson'
cyup      write(*,*)'integral',integral

c----------------------------------------------------------------------                  
 20   continue

      integral(2,1)=-integral(2,1)
      integral(3,1)= integral(3,1)
      integral(3,2)=-integral(3,2)
           
 10   continue

      return
      end   

c======================================================================
c======================================================================

      subroutine calc_absorption_integral_array_by_adaptive_simpson_xyz
     &(wccw,nll,np,theta,n,i_fkin,x,y,z,p_permax,epsi,integral)
c------------------------------------------------------------------------
c     calculate intergals integral(3,3) for anti-hermitian diectric tensor
c     for the EC harmonic with number 'n'
c-------------------------------------------------------------------------
       implicit none
c------------------------------------------------------------------------
c      input
c------------------------------------------------------------------------
      integer n,i_fkin
      real*8 wccw,nll,np,theta,x,y,z,p_permax,epsi 
c     n        is EC harmonic number
c     i_fkin   =0 the usage of the analytical relativistic distributin
c              =1 the usage of the numerical 3D distribution from diskf 
c                   written be CQL3D code
c     nll      is N_parallel
c     theta    = mc**2/T
c     x,y,z    cartesian coords
c     p_permax is  the maximal boundary of integration
c     epsi     is the accuracy used in integration by adaptive Simpson
c------------------------------------------------------------------------
c     output
c------------------------------------------------------------------------
      double complex integral(3,3)
c-------------------------------------------------------------------------
c     locals
c--------------------------------------------------------------------------
      integer k_max !max number of elements
      parameter (k_max=6)
      integer ist,ifs,k   
c      PARAMETER (IST=5,IFS=16)
c      PARAMETER (IST=9,IFS=256)  
c       PARAMETER (IST=13,IFS=10048)
c       PARAMETER (IST=16,IFS=100048)
c       PARAMETER (IST=17,IFS=200048)!070723-old
      PARAMETER (IST=20,IFS=200048)

      REAL*8 STACK(IST,3),FSTACK(IFS,3)
      real*8 a,b,sum(k_max)
      integer lvmax,iflag,i
      real*8 fcn(k_max)
    
      integer n_l,i_fkin_l
      real*8 wccw_l, nll_l, np_l, theta_l, x_l, y_l, z_l 
      common /fcn_input_xyz/wccw_l,nll_l,np_l,theta_l,n_l,i_fkin_l,
     &x_l,y_l,z_l
 
      EXTERNAL FCN_absorption_vector_xyz
c      DATA EPSI/5.0d-5/, LVMAX/4/
c      DATA EPSI/5.0d-5/, LVMAX/7/
c      DATA EPSI/5.0d-5/, LVMAX/15/
c      DATA EPSI/5.0d-5/, LVMAX/12/
c      DATA EPSI/5.0d-5/, LVMAX/16/
c       DATA EPSI/5.0d-6/, LVMAX/16/  !070722-new
c       DATA EPSI/1.0d-6/, LVMAX/16/  !070722-new
c       DATA EPSI/1.0d-3/, LVMAX/16/ !070722-old
c      DATA EPSI/1.0d-3/, LVMAX/19/
      DATA LVMAX/16/  !070722-new
c-----------------------------------------------------------------
c     set common /fcn_input_xyz/
c-----------------------------------------------------------------
      wccw_l=wccw
      nll_l=nll
      np_l=np
      theta_l=theta
      n_l=n
      i_fkin_l=i_fkin
      x_l=x
      y_l=y
      z_l=z
c-----------------------------------------------------------------
c     integration boundaries
c------------------------------------------------------------------
c      A = 0.0d0
      A = 1.d-5
c      A = 1.d-3

      B = p_permax-1.d-5
     
      CALL ASMP(FCN_absorption_vector_xyz,k_max,fcn,
     &A,B,EPSI,LVMAX,STACK,IST,FSTACK,IFS,
     &SUM,IFLAG) 

      IF(IFLAG .EQ. 0) THEN 
c$$$        PRINT 4   
      ELSE
        write(*,*)'          WITH BAD SUBINTERVALS:'   
        DO 2 I=1,IFLAG      
          write(*,*) FSTACK(I,1),FSTACK(I,2),INT(FSTACK(I,3))
    2   CONTINUE  
      END IF      
   3  FORMAT(//5X,'APPROXIMATE INTEGRAL =',E22.14/)   
   4  FORMAT(10X,'WITH NO BAD SUBINTERVALS')    
   5  FORMAT(10X,'WITH BAD SUBINTERVALS:')      
   6  FORMAT(10X,'[',F10.5,',',F10.5,']',2X,'LEVEL =',I5) 
c-------------------------------------------------------
      integral(1,1)=dcmplx(Sum(1),0.d0) 
      integral(1,2)=dcmplx(0.d0,Sum(2))
      integral(1,3)=dcmplx(Sum(3),0.d0) 
      integral(2,2)=dcmplx(Sum(4),0.d0)
      integral(2,3)=dcmplx(0.d0,Sum(5))
      integral(3,3)=dcmplx(Sum(6),0.d0) 
      return
      END 

c======================================================================
c======================================================================

c======================================================================
c     Adaptive Simpson functions.
C
c     Modificated method for vector function.
c     
c     Original description and functions are from:
c
C PAGE 187-190: NUMERICAL MATHEMATICS AND COMPUTING, CHENEY/KINCAID,1985
C
C FILE: SIMP.FOR
C
C ADAPTIVE SCHEME FOR SIMPSON'S RULE (SIMP,ASMP,PUSH,POP,FCN)
c======================================================================

      subroutine fcn_absorption_vector_xyz(p_perp,fcn)
c-----------------------------------------------------------------------
c     calc. under-integral function for anti-hermitian dielectric tensor
c-----------------------------------------------------------------------
      implicit none
c      include 'param_kmax.i'
      integer k_max !max number of elements
      parameter (k_max=6) 
c-----------------------------------------------------------------
c     input
c----------------------------------------------------------------- 
      real*8 p_perp !under integral function argument
c-----------------------------------------------------------------
c     output
c-----------------------------------------------------------------
      real*8 fcn(k_max) !array of under integral fuctions values
c-----------------------------------------------------------------
      integer n_l,i_fkin_l
      
      real*8 wccw_l, nll_l, np_l, theta_l, x_l, y_l, z_l 
      common /fcn_input_xyz/wccw_l,nll_l,np_l,theta_l,n_l,i_fkin_l,
     &x_l,y_l,z_l
c-----------------------------------------------------------------
c     locals
c-----------------------------------------------------------------
      double complex g_n(3,3) !complex under-integral functions

      call g_n_calc_xyz(p_perp,wccw_l,nll_l,np_l,theta_l,n_l,i_fkin_l,
     & x_l, y_l, z_l, 
     + g_n) !->out

      FCN(1) =dreal(g_n(1,1)) 
      FCN(2) =dimag(g_n(1,2))
      FCN(3) =dreal(g_n(1,3))
      FCN(4) =dreal(g_n(2,2))
      FCN(5) =dimag(g_n(2,3))
      FCN(6) =dreal(g_n(3,3))

      return
      end



c======================================================================
c======================================================================



      subroutine g_n_calc_xyz(p_perp,wccw,nll,np,theta,n,i_fkin,
     +  x,y,z,g_n)
c     under integral complex marix function G_n(3,3)
c     input
c       p_perp  - momentum divided by (mc)
c       wccw    = omega_ce/omega for the electron rest mass    
c       nll     - N_parallel
c       np      - N_perpendicular 
c       theta   - mc**2/T
c       n       - n number of the given resonance harmonic
c       i_fkin =0 the usage of the analytical relativistic distributin
c              =1 the usage of the numerical 3D distribution from diskf 
c                 written be CQL3D code
c       x,y,z   - cartesian coords
c     output
c       g_n(3,3) under integral matrix double complex function
c-------------------------------------------------------
      implicit none
c-----input        
      double precision p_perp,wccw,nll,np,theta,x,y,z
      integer n,i_fkin
c-----external root_res, s_calcn, zeroK, u_n
      double precision u_n    

c-----output
      double complex g_n(3,3) 

c-----local
      double complex sn_k(3,3),g_nk(3,3)
      double precision gamma, p_par_rl(2),coeff,eps,ps,ps_max,pi,
     +dgdp_par,u_nk
      integer k,kmax,i,j
c     kmax - total number of the resonace condition root p_perp_k(p_perp)

      pi=4.d0*datan(1.d0)
c-----calculations ot the roots p_par_rl of the resonance condition
c-----gamma=N_par*p_par_rl+nY

      kmax=0 ! to initialize
      call root_res(p_perp,nll,n,wccw,kmax,p_par_rl)

c-----initialize g_n
      call zeroK(g_n)

      if (kmax.eq.0) goto 20
  
c      eps=1.d-9 ! the min value of Maxwell exponent
c      p s_max=(1.d0-dlog(eps)/theta)**2-1.d0
c      write(*,*)'g_n_calc kmax',kmax      
      do k=1,kmax ! the sum over all P_parallel_k(perp) roots 
   
         ps=p_perp*p_perp+p_par_rl(k)*p_par_rl(k)
         gamma=dsqrt(1.d0+ps)

c         if(ps.gt.ps_max) then
c           write(*,*)'relat_tens g_n k,ps,ps_max',k,ps,ps_max
c           goto 10
c         endif
 
c        write(*,*)'g_n_calc before s_calcn: k',k  
         call s_calcn(n,p_perp,p_par_rl(k),np,wccw,sn_k)
c         write(*,*)'g_n_calc after s_calcn: k,sn_k',k,sn_k  

c--------resonance condition uses the delta function with argument
c        g_delta=gamma-nll*p_par-n*wccw
c        the derivative from this argument d(g_delta)/dp_par=dgdp_par

c         dgdp_par=(p_par_rl(k)-nll*wccw)/gamma
         dgdp_par=dabs((p_par_rl(k)-nll*gamma)/gamma)

         u_nk=u_n(p_perp,p_par_rl(k),wccw,nll,theta,n,i_fkin)

         coeff=2.d0*pi*u_nk/dgdp_par

         g_nk(1,1)=coeff*sn_k(1,1)  
         g_nk(1,2)=coeff*sn_k(1,2)
         g_nk(1,3)=coeff*sn_k(1,3)
         g_nk(2,2)=coeff*sn_k(2,2)
         g_nk(2,3)=coeff*sn_k(2,3)
         g_nk(3,3)=coeff*sn_k(3,3)

         g_nk(2,1)=-g_nk(1,2)
         g_nk(3,1)= g_nk(1,3)
         g_nk(3,2)=-g_nk(2,3)        

         do i=1,3
            do j=1,3
               g_n(i,j)=g_n(i,j)+g_nk(i,j)
            enddo
         enddo
 10      continue
       
      enddo !kmax

 20   continue

      return
      end


c======================================================================
c======================================================================



c        ********************** prep3d_xyz ******************
c        *                      -----                       
c        *  prep3d_xyz - to prepare the data 
c        *  for	output files: netcdf, etc.  
c        ****************************************************
c         input parameters: iray -number of the ray from antenna  
c                                 is in common/cone/ 
c output parameter: iraystop (if power in the ray channel 
c   delpwr(is).lt. delpwrmn*delpwr(1), then iraystop=1). 
c Also, stop ray if fluxn(is).lt.0.0, iraystop=1, and set fluxn 
c  to previous value fluxn(is-1).   [RWH:030428]
c-----------------------------------------------------------------------
      subroutine prep3d_xyz(u,deru,iraystop)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'write.i'
      include 'cefield.i'
      include 'cone.i'
      include 'grill.i'
      include 'eps.i'
      include 'fourb.i'
      include 'oxb.i'
      include 'output.i'
      include 'three.i'
cSm070128 for test
      include 'six.i'

      dimension u(*),deru(*),vgr(3),bf(3)
c      dimension u(6),deru(6),vgr(3),bf(3)
      
      dimension cnprim_s(nbulka)  !Added for indiv ion contrib[BH041009]
      dimension ckvipl_s(nbulka)  !Added for indiv ion contrib[BH041009]
      
      dimension tempiar(nbulka)
C------ 2 is maximum number of roots to look for in the muller root finding
      double complex cflown,cnprp,fr_func_noeps,dfrfunc,roots(2)
      integer info(2),ier,nsig
      double complex dhot,dhot_rlt,dhot_sum

      double precision x_ar(nbulka),y_ar(nbulka),t_av_ar(nbulka),
     .tpop_ar(nbulka),vflow_ar(nbulka),image_d
      
      double complex hotnp
      double complex dhot_sumt,dhot_sum_e,dhot_sum_i
cfor_test
      double complex eps1(3,3),cez1,eps2(3,3),cez2
      double complex cd 
CENM For the iabsorp=12 damping option to work (find full complex solution of
C    D(nperp)=0), the muller root-finding subroutine is called, and either
C    Ram's relativistic function is used (for id=14) or Nelson-Melby's (for id=11)
      external fr_func_noeps ! for dispersion function root finder (id=14)

cend_for_test

      double complex integral(3,3),fluctcur_n(3,3)
cfor emission test

cfor test ono tensor
      double complex K_dx_ar(3,3,nbulka),dK_dy_ar(3,3,nbulka),
     &dK_dt_ar(3,3,nbulka),dK_dnper(3,3),dK_dnpar(3,3)
      double precision mass_ar(nbulka)
cendtest
c     to calculate data for ploting the resonance ellipse boundary
      double precision p_par_rl(2)

!      external dhot

cSm030226      
c      integer nbulkaa
c      parameter (nbulkaa=5)
      double complex K(3,3),dK(3,3,7),dd(5),d,ddn,aK(3,3)
c      double complex dd5(5,nbulkaa),ddnp_h,ddnll_h,ddnp
      double complex dd5(5,nbulka),ddnp_h,ddnll_h,ddnp
 
c-----------------------------------------------------------------------
c     for test of relativistic absorption
      double complex eps_test(3,3)
c-----------------------------------------------------------------------
c     for test grpde2
      complex exde,eyde,ezde
      real rnpar,rnper,romega,romegpe,rvgrpdc,redenfac
c     for cold plasma+relativistic tensor
      double complex disp_func,dcold_rlt

c-----for cpu_time
      real*4 time_prep3d_1,time_prep3d_2,time_prep3d
     &time_prep3d_emis_1,time_prep3d_emis_2,time_prep3d_emis

c-----to check eigenvalues
      complex*16   eigenvalue(3)

c-----to check hermitian part of K
      complex*16 K_herm(3,3)

      double precision optical_depth
c-----real*8 p_flux !for flux calculations


c-----for reflection lost 
      integer irefl_old
      real*8  tot_pow_absorb_at_refl
 
c-----for test N perpendicular coinside at the given ray pint
c     with the dispersion relation solution
      real*8
     &cnper_p,cnper_m      
     
      data irefl_old /0/
      data tot_pow_absorb_at_refl/0.d0/
c-------------------------------------------
      data optical_depth/0.d0/

      data time_prep3d /0./
      data time_prep3d_em /0./

      save cnprim_old
      save time_prep3d,time_prep3d_em
      save  optical_depth
      save irefl_old
      save tot_pow_absorb_at_refl

      save ckvipold
      call cpu_time(time_prep3d_1)

      nrayelt=nrayelt+1
      if (nrayelt.gt.nrelta) then
         write(*,*)'in prep3d nrayelt.gt.nrelta nrayelt,nrelta',
     .   nrayelt,nrelta
         write(*,*)'it should be nrayelt=<nrelta'
         write(*,*)'increase nrelta in param.i'
         stop
      endif  
      is=nrayelt

      iraystop=0
      pi=4.d0*datan(1.d0)
      cnprim=0.d0 ! to initialize
      cnper2p=0.d0
      cnper2m=0.d0 ! to initialize
      cnper=0.d0 ! to initialize
      cflown=(0.d0,0.d0) ! to initialize
      f000=frqncy
c----------------------------------------
c     cvac (cm/sec)
      cvac=2.997930d+10
c----------------------------------------
c     cld (cm),frgncy(GHz)
      cld=cvac/(2.d0*pi*frqncy*1.0d+09)
c----------------------------------------
c     now proposed that r0x=1 m
      r00=100.d0*r0x
      t00=cvac/(2.d0*pi*frqncy*r00)
c-----------------------------------------
      x= u(1) 
      y= u(2) 
      z= u(3)   ! [m]
      cnx= u(4) ! Nx
      cny= u(5) ! Ny
      cnz= u(6) ! Nz
c----------------------------------------------------
c     poloidal angle 0< theta_pol [degree] <360
c-----------------------------------------------------
c      write(*,*)'prerp3d r,z,rma,zma', r,z,rma,zma
      r=dsqrt(x*x+y*y)
      call theta_rz((r-rma),(z-zma),wtheta_pol(is))
      wtheta_pol(is)=(wtheta_pol(is))*180d0/pi ! Yup: only for write/xdraw?
      
      if (is.gt.1) then ! Yup: only for write/xdraw?
         if ((wtheta_pol(is-1).le.90d0).and.
     &       (wtheta_pol(is).ge.180.d0)) then
            wtheta_pol(is)=wtheta_pol(is)-2*180.d0
         endif
      endif
c---------------------------------------
c     bmod,bx,by,bz (Tl)
      bmod=bxyz(x,y,z)  !-> get b and derivs of b
      bf(1)=bx
      bf(2)=by
      bf(3)=bz
      gam=gamma1_xyz(x,y,z,cnx,cny,cnz)
      ds=dsin(gam)
      dc=dcos(gam)
      cnt= dsqrt(cnx*cnx+cny*cny+cnz*cnz)

      cnpar=cnt*dc
      cnper=cnt*ds
      
      dens_e= dense_xyz(x,y,z,1) !-> get rho and dens_e

c-----waves absorption and electric field calculations----
      if ((iabsorp.eq.3).or.(iabsorp.eq.2)) then
c	absorption for lh and fw
c------------------------------------------------------------
c       electric field using the cold plasma dielectric tensor
        call tensrcld_xyz(u(1),u(2),u(3))
        cnprp=dcmplx(cnper,0.d0)
        call efield1(cnpar,cnprp,ex,ey,ez,eplus,eminus,epar)
c       electric field parallel to wave
        enp=(ex*cnper+ez*cnpar)/dsqrt(cnper**2+cnpar**2)
        
c-------put cold plasma tensor to w_ceps array
        do i=1,3
           do j=1,3
              w_ceps(i,j,is)=reps(i,j) !from  eps.i
           enddo
        enddo
           
c-------------------------------------------------------------
        temp_e=tempe_xyz(x,y,z,1)
        do i=2,nbulk
          tempiar(i)=tempe_xyz(x,y,z,i)
        enddo

        z_eff=zeffrho(rho) 

        if (iabsorp.eq.3) then
c----------FW absorption
           cnprim_cl=0.d0
c----------absorpfd uses complex function modified bessel zfunc(argz,zf,zfp) 
c           call absorpf1(z,r,phi,cnpar,cnper,temp_e,dens_e,tempiar,
c----------absorpfd uses double complex function modified bessel 
c          czeta(argz,zf,zfp,ierror) and calculates the dielectric tensor
c          reps )(in common eps.i)
c          for the electron plasma with the hot correction using
c          Chiu et al, Nucl.Fus Vol. 29, No.12(1989) p.2175
c          formula (2),(3),(4),(5),(6) and (7)
           rho_larm_max0=rho_larm_max !YuP[11-2016]
           call absorpfd_xyz(x,y,z,f000,rho_larm_max0,cnpar,cnper,
     1      temp_e,dens_e,tempiar,nbulk,bmod,cnprim_e,cnprim_i,cnprim_s)

           if (ion_absorption.eq.'enabled') then 
              cnprim=cnprim_e+cnprim_i
           else
              cnprim=cnprim_e !only electron absorption
           endif 

           if (cnprim_e.lt.0.d0 .or. cnprim_i.lt.0.d0) then
             !YuP[04-2016] Sometimes cnprim_s are negative,
             !especially for ion damping (case of HHFW, for example).
             !If they are negative, they are reset to abs().
             !Before 04-11-2016: 
             !         the negativity of total(e+i) cnprim was checked.
             !After  04-11-2016: each of cnprim_e, cnprim_i is checked.
             write(*,'(a,4e12.3,a)')
     +        'prep3d_xyz: cnprim<0. cnprim_e,cnprim_i,cnpar,cnper='
     +        ,cnprim_e,cnprim_i, cnpar,cnper,
     +        '  Resetting to abs(cnprim)'
             !if(cnprim_i.lt.0.d0) cnprim_i=0.d0 ! another version
             !if(cnprim_e.lt.0.d0) cnprim_e=0.d0 ! another version
             cnprim_e=dabs(cnprim_e)
             cnprim_i=dabs(cnprim_i)
             do i=1,nbulk
                cnprim_s(i)=dabs(cnprim_s(i))
             enddo
             !YuP: But is it a valid change? Maybe set to 0, if cnprim<0?
           endif

           if (ion_absorption.eq.'enabled') then 
              ! repeated, after possible resetting of
              ! cnprim_i to abs(cnprim_i)
              cnprim=cnprim_e+cnprim_i
           else
              cnprim=cnprim_e !only electron absorption
           endif 

c----------electric field calculations using the dielectric tensor
c          from absorpfd (the electron plasma with the thermal correction)
c
c          electric field using the cold plasma dielectric tensor
           cnprp=dcmplx(cnper,cnprim)

           do i=1,3
              do j=1,3
                 w_ceps(i,j,is)=reps(i,j) !cold plasma with thermal correction
              enddo
           enddo
         
           call efield1(cnpar,cnprp,ex,ey,ez,eplus,eminus,epar)
c          electric field parallel to wave
           enp=(ex*cnper+ez*cnpar)/dsqrt(cnper**2+cnpar**2)

!           if (is.eq.1000000) then
!c------------ test N_perp and polarization in one given ray point WF case
!              write(*,*)'is=81'
!              write(*,*)'cnpar,cnper',cnpar,cnper
!              write(*,*)'cex,cey,cez',cex,cey,cez
!
!              call npernpar_xyz(x,y,z,cnpar, cnper2p,cnper2m)
!
!              if (cnper2p.ge.0.d0) then
!                 cnper_p=dsqrt(cnper2p)
!                 write(*,*)'cnper_p,cnper',cnper_p,cnper
!              endif
!
!              if (cnper2m.ge.0.d0) then
!                 cnper_m=dsqrt(cnper2m)
!                 write(*,*)'cnper_m,cnper',cnper_m,cnper
!              endif
!
!cyup              psi_loc=psi_rho(rho)
!              psi_loc=psif_xyz(x,y,z) ! YuP Added instead of above
!              x_m=x  
!              y_m=y  
!              z_m=z  
!              r_m=dsqrt(x*x+y*y)
!              rho_loc=dsqrt((r_m-rma)**2+(z_m-zma)**2) ! needs work
!              cos_theta_pol=(r_m-rma)/rho_loc
!              sin_theta_pol=(z_m-zma)/rho_loc
!              if (cos_theta_pol .gt. 1.d0)  cos_theta_pol= 1.d0
!              if (cos_theta_pol .lt. -1.d0)  cos_theta_pol=-1.d0
!              if (sin_theta_pol.ge.0.d0) then
!                 theta_pol=+dacos(cos_theta_pol)
!              else  
!                 theta_pol=-dacos(cos_theta_pol)
!              endif
!
!              if (theta_pol.lt.0.d0) then
!                 theta_pol=theta_pol+2.d0*pi
!              endif
!
!              if (theta_pol.gt.2.d0*pi) then
!                 n_theta_pol=theta_pol/(2.d0*pi)
!                 theta_pol=theta_pol-2.d0*pi*n_theta_pol
!              endif
!
!c-------------calculate the hot plasma full dielectric tensor and the 
!c             electric field
!                 x_ar(i)=wpw_2(x,y,z,i)
!            endif !is.eq.81)
!          
!c-----------end test

        endif ! iabsorp=3

        if (iabsorp.eq.2) then
c----------LH absorption
           !write(*,*)'prep3_xyz: before absorplh_xyz'
           call absorplh_xyz(u,cnpar,cnper,temp_e,dens_e,tempiar
     1                  ,bx,by,bz,nbulk,bmod,frqncy,z_eff,
     1                   cnprim_e,cnprim_i,cnprim_cl)
           !write(*,*)'prep3_xyz: after absorplh_xyz'
        endif !iabsorp=2

         if (ion_absorption.eq.'enabled') then
            cnprim=cnprim_e+cnprim_i
         else
            cnprim=cnprim_e
         endif
         
         cnprim=cnprim+cnprim_cl

      endif !iabsorp=2 or =3



      if (iabsorp.eq.1) then
c-------EC wave case.The complex electric field calculations using
c       hermitian or full mazzucato tensor (with antihermitian part)
c       reps(3,3), and N_per=N_per(N_par) from mazzucato solver
c       ihermloc=iherm
c       ihermloc=1
        ihermloc=2
c-------EC absorption (from Mazzucato code)
        ioptmaz=1 ! estimation of cnper1 from coldm
        cnper1=cnper
        ihermloc=2
        call cnpermuz_xyz(cnpar,ihermloc,x,y,z,cnper1,cnprim,ioptmaz)
        ioptmaz=2 ! estimation of cnper1 from input parameter cnper
        cnprim_cl=0.d0
        cnprim_e=cnprim
        cnprim_i=0.d0
c-------electric field for mazzucato tensor
c        ihermloc=iherm
        ihermloc=2 !full hermitian + antihermitian Mazzucato tens. 
        call hamiltmuz_ful_xyz(cnpar,ihermloc,x,y,z,cnper1,
     .  cnprim,hamiltmz)
c-------put Mazzucato tensor from complex N perpendicular to w_ceps array
        do i=1,3
           do j=1,3
              w_ceps(i,j,is)=reps(i,j) !from  eps.i
           enddo
        enddo
        cnprp=dcmplx(cnper1,cnprim)
        call efield1(cnpar,cnprp,ex,ey,ez,eplus,eminus,epar)
cyup        enp=(ex*cnper+ez*cnpar)/dsqrt(cnper**2+cnpar**2) !print-out?
      endif ! (iabsorp.eq.1) 


      if(iabsorp.eq.4) then
c       EC wave case.The complex electric field calculations
c       using Forest tensor (with antihermition part)
c       reps(3,3), and N_per=N_per(N_par) from forest solver
c-------EC absorption (from Forest code)

        cnparp=cnpar
        cnperp=cnper
     
        if (((id.eq.1).or.(id.eq.2)).or.(id.eq.3)) then
cc----------calculate N_perp(N_par) from hot plasma disp.relation   
cc          calculates two roots from the cold plasma as the initial
cc          approximation 
c           call npernpar(z,r,phi,cnpar,cnper2p,cnper2m)

cc----------set the data to common npercom.i for hotnp function
            call set_nperpcom_xyz(cnpar,x,y,z,dmas)
cc          calculates Nper(Npar) from hot plasma with the initial 
cc          iteration cnper from cold plasma
c           hotnperp_xyz=hotnperp_xyz(nbulk,ibw,cnperp,cnper2p,cnper2m,K,iraystop)
c           cnperp=hotnperp_xyz
        endif      
 
        do i=1,nbulk
           x_ar(i)=wpw_2(x,y,z,i)
           y_ar(i)=wcw(x,y,z,i)
           if(i.eq.1) y_ar(1)=-y_ar(1)
           te=tempe_xyz(x,y,z,i)        ! kev
           t_av_ar(i)=te*1000.d0      ! ev 
           tpop_ar(i)=tpoprho(rho,i)
           vflow_ar(i)=vflowrho(rho,i)
        enddo

        if(i_im_nperp.eq.1) then  !......................
c---------calculate ImN_perp using the formula
c         ImN_perp=abs(ImD_full/dD_hermitian/dReN_perp))
c         write(*,*)'in prep3d_xyz-1: x,y,z=',x,y,z
          !write(*,*)'in prep3d_xyz-1: x,y,z, Y=',x,y,z, y_ar(1)
          d=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     .    vflow_ar,cnparp,cnperp,2,reps)
          dham=dreal(d)
          call Ddhot(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,vflow_ar
     .    ,cnparp,cnperp,reps,dd5,ddnp_h,ddnll_h,ddnp)
          cnprim = dabs(DIMAG(D) / DREAL(ddnp))
          write(*,*)'cnprim=(ImD/dD/dn_perp)= ',cnprim
          call hot_nperp_muller(nbulk,dmas,x_ar,y_ar,t_av_ar,
     &    tpop_ar,vflow_ar,cnparp,cnperp,cnprim)
          write(*,*)'muller cnprim',cnprim 
          do i=1,3
             do j=1,3
                w_ceps(i,j,is)=reps(i,j) !hot plasma tensor
             enddo
          enddo
          cnper_new=cnperp !it will be used in the electric field calculations
        endif  ! i_im_nperp.eq.1   !.....................

        if(i_im_nperp.eq.2) then   !:::::::::::::::::::::
c---------find (Im_N_perp,ReN_perp) the root of the complex dispersion relation
c         using the Newton method with numerical derivatives (the chord method)
          iter_max=100 !max number of the iterations
          iter_max=3 !max number of the iterations
c---------initial values of Im_N_perp=cnprim Re_N_perp=cnperp
          cnprim=0.d0   
          write(*,*)'prep3d before call solv_nperp_hot cnprim=',cnprim
          write(*,*)'nbulk,dmas,x_ar,y_ar',
     &               nbulk,dmas,x_ar,y_ar
          write(*,*)'t_av_ar',t_av_ar
          write(*,*)'tpop_ar',tpop_ar      
          write(*,*)'vflow_ar',vflow_ar
          write(*,*)'cnparp,cnperp,iter_max,cnper_new,cnprim',
     &               cnparp,cnperp,iter_max,cnper_new,cnprim
          call solv_nperp_hot(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     .    vflow_ar,cnparp,cnperp,iter_max,cnper_new,cnprim)
          write(*,*)'Newton cnperp,cnper_new,cnprim',
     .    cnperp,cnper_new,cnprim
          cnprim_old=cnprim
          cnprim=dabs(cnprim)
        endif !  i_im_nperp.eq.2 ::::::::::::::::::::::::

        cnprim_cl=0.d0
        cnprim_e=dabs(cnprim)
        cnprim_i=0.d0

 30     continue
        cnprim=dabs(cnprim)

c-------electric field for Forest tensor
c       cnprp=dcmplx(cnper,cnprim)
        cnprp=dcmplx(cnper,0.d0)

c-------full tensor with new n_perp calculated by solver
          !write(*,*)'in prep3d_xyz-2: x,y,z, Y=',x,y,z, y_ar(1)
        d=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     .            vflow_ar,cnparp,cnper,2,reps)      
        cnprp=dcmplx(cnper_new,cnprim)
        cnprp=dcmplx(cnper,0.d0)
 
        call efield1(cnpar,cnprp,ex,ey,ez,eplus,eminus,epar)
        enp=(ex*cnper+ez*cnpar)/dsqrt(cnper**2+cnpar**2)         
      endif ! iabsorp.eq.4

c------------------------------------------------------------------
      if(iabsorp.eq.6) then
c       EC wave case.The complex electric field calculations
c       using Forest tensor (with antihermition part)
c       reps(3,3), and N_per=N_per(N_par) from forest solver
c-------EC relativistic absorption 
c       dielectric tensor=hermitian part(from Forest code)+
c                         anti-hermitian part(full relativistic) 
	
        cnparp=cnpar
        cnperp=cnper
        
c-------Hermitian non-relativistic tensor reps        

        do i=1,nbulk
           x_ar(i)=wpw_2(x,y,z,i)
           y_ar(i)=wcw(x,y,z,i)
           if(i.eq.1) y_ar(1)=-y_ar(1)
           te=tempe_xyz(x,y,z,i) ! kev
           t_av_ar(i)=te*1000.d0      ! ev 
           tpop_ar(i)=tpoprho(rho,i)
           vflow_ar(i)=vflowrho(rho,i)
        enddo
 
cyup        d=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
cyup     .   vflow_ar,cnparp,cnperp,1,reps) ! YuP: not used?
        if ((i_diskf.eq.0).or.(i_diskf.eq.5)) then
c----------usage of the analytical relativistic function and its derivatives 
           i_fkin=0
        else 
c----------usage of the mech relativistic function and its derivatives
           i_fkin=1
        endif

        call anth_rlt_xyz(x_ar(1),y_ar(1), t_av_ar(1)*1.d-3,
     +   cnparp,cnperp,
     +  n_relt_harm1,n_relt_harm2,n_relt_intgr,
     +  i_resonance_curve_integration_method,epsi,
     +  i_fkin,x,y,z,
     +  aK)
         
        cnprimp=0.d0
        d=dhot_rlt(reps,aK,cnparp,cnperp,cnprimp)
        dham=dreal(d)
        !write(*,*)'reps',reps
        !write(*,*)'aK',aK
        !write(*,*)'d',d
        call Ddhot(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,vflow_ar
     . ,cnparp,cnperp,reps,dd5,ddnp_h,ddnll_h,ddnp)
        
        cnprim = dabs(DIMAG(D) / DREAL(ddnp))

        !write(*,*)'dimag(d)',dimag(d)
        !write(*,*)'dreal(ddnp),cnprim',dreal(ddnp),cnprim

        cnprim_cl=0.d0
        cnprim_e=cnprim
        cnprim_i=0.d0

c-------electric field for Forest tensor
        cnprp=dcmplx(cnper,cnprim)
        cnprp=dcmplx(cnper,0.d0)

          !write(*,*)'in prep3d_xyz-3: x,y,z, Y=',x,y,z, y_ar(1)
        d=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     .  vflow_ar,cnparp,cnperp,2,reps)

        call efield1(cnpar,cnprp,ex,ey,ez,eplus,eminus,epar)
        enp=(ex*cnper+ez*cnpar)/dsqrt(cnper**2+cnpar**2)
      endif ! iabsorp.eq.6
c-------------------------------------------------------------
      if(iabsorp.eq.7) then
c       EC wave case.The complex electric field calculations
c       using Cold plasma tensor +antihermition relativistic tensor
c-------EC relativistic absorption 
c       dielectric tensor=hermitian part(cold plasma)+
c                         anti-hermitian part(full relativistic) 
	
        cnparp=cnpar
        cnperp=cnper
        
c-------calculate Hermitian cold plasma complex tensor reps. 
c       It will be in eps.i        
        call tensrcld_xyz(x,y,z)

        do i=1,nbulk
           x_ar(i)=wpw_2(x,y,z,i)
           y_ar(i)=wcw(x,y,z,i)
           if(i.eq.1) y_ar(1)=-y_ar(1)
           te=tempe_xyz(x,y,z,i) ! kev
           t_av_ar(i)=te*1000.d0      ! ev 
           tpop_ar(i)=tpoprho(rho,i)
           vflow_ar(i)=vflowrho(rho,i)
        enddo
        
        if ((i_diskf.eq.0).or.(i_diskf.eq.5)) then
c----------usage of the analytical relativistic function and its derivatives 
           i_fkin=0
        else 
c----------usage of the mech relativistic function and its derivatives
           i_fkin=1
        endif
         
        call anth_rlt_xyz(x_ar(1),y_ar(1), t_av_ar(1)*1.d-3,
     +    cnparp,cnperp,
     +    n_relt_harm1,n_relt_harm2,n_relt_intgr,
     +    i_resonance_curve_integration_method,epsi,
     +    i_fkin,x,y,z,
     +    aK)
        
c------ complex dispersion function calculated from the sum of
c       of the cold electron plasma dielectric tensor eps_h
c       and the relativistic electron anti-hermition dielectric tensor eps_a

        disp_func=dcold_rlt(reps,aK,cnparp,cnperp)

c-------calculate the derivative d(D_hermitian)/d(ReN_perp)
c       from the electron cold plasma dispersion function D
        ddnp=dDcold(reps,cnpar,cnper)
        
        cnprim = dabs(DIMAG(disp_func) / DREAL(ddnp))
	
        cnprim_cl=0.d0
        cnprim_e=cnprim
        cnprim_i=0.d0
     
c-------electric field for cold plasma
        cnprp=dcmplx(cnper,cnprim)
        cnprp=dcmplx(cnper,0.d0)
        
        call efield1(cnpar,cnprp,ex,ey,ez,eplus,eminus,epar)
        enp=(ex*cnper+ez*cnpar)/dsqrt(cnper**2+cnpar**2)

      endif ! iabsorp.eq.7


c----------------------------------------------------------------
      if(iabsorp.eq.10) then
c--------------------------------------------------------------
c        The absorption is calculated for relativistic dispersion
c        (combined E. Nelson-Melby  and  A.Ram)
c        using the formula from Stix book p.74 (17,18,21)
c        Im(k_perp)= 0.5*Power_abs/(P^+T^) 
c
c        Here 
c    
c        Power_abs=omega/(8pi)[ E~(i) . (eps_a_herm(i,j) . E(j)]
c
c        P^ = (c/16pi)[E~*B+E*B~] is Poining vector,calculated
c             using hot plasma complex dieletric tensor.
c
c        T^ = -omega/(16pi)[ E~(i) . d/dk^(eps_herm(i,j) . E(j)]
c             Is a flux of nonelectromagnetic energy
c----------------------------------------------------------------
         call absorp_relativist_disp_combined_xyz(x,y,z,cnpar,cnper,
     &   cnprim_e)
 
         cnprim_i=0.d0     
         cnprim=cnprim_e+cnprim_i
         cnprim_cl=0.d0
        
         write(*,*)'after absorp_relativist_disp_combined'
         write(*,*)'cnprim_e',cnprim_e

c----------------------------------------------------------------
c        electric field polarization (cex,cey,cez)
c        calculated using the relativistic dielectric tensor
c        will be common /efield.i/
c---------------------------------------------------------------
c        relativistic (full, non-hermitian) dielectric tensor reps(3,3)
c        will be in common /eps.i/
c---------------------------------------------------------------        
        cnprp=dcmplx(cnper,cnprim)
        call efield1(cnpar,cnprp,ex,ey,ez,eplus,eminus,epar)
         
      endif !iabsorp=10	

C---------------------------------------------
      if (iabsorp.eq.12) then
C------------------------------------------
C       Test basic elements needed for iabsorp

C----- use Muller algorithm to find dispersion relation root given
C----- n_parallel and frequency and plasma parameters.
         
C************ Here are some hard-coded numbers for the accuracy with
C************ which to search for the root for the iabsorp.eq.12 option
C************ Since it should be sticking very close, usually the iterations
C************ would not be much of an issue.
         errabs=1.d-6
         nsig=6
         nknown=0
         nrts=1                 ! just look for one root
         nguess=nrts 
         nnew=nrts
         itmax=50
C******* For initial guess, just search with the real part from before and just
C******* 0 imaginary part. Usually, if the ray is propagating, the damping
C******* is low anyway, so it should not be too large imaginary part.
         roots(1)=dcmplx(cnper,0.0d0)
      write(*,*)'()()(()()(()()initial data'
      write(*,*)'x=',x,'y=',y,'z=',z
      !write(*,*)'cnx=',cnx,'cny=',cny,'cnz=',cnz
      write(*,*)'rho=',rho
      print *,'====== is=',is
      print *,'cnprim: ',cnprim,' cnpar: ',cnpar,' cnper: ',cnper
      
         cn2=cnx**2+cny**2+cnz**2

         print *,'cn2=',cn2,' cnper(calculated) ',dsqrt(cn2-cnpar**2)
         if (is.lt.1) then
            print *,'****** is=',is
         else
            print *,'nper before: ',wnper(is-1)
         endif

C id.eq.14, call fr_func_noeps, using Ram's dielectric function
c-----------print relativistic tensor for testing
            X_e=wpw_2(x,y,z,1)
            Y_e=wcw(x,y,z,1)
            T_e=tempe_xyz(x,y,z,1)  
            cnprp=dcmplx(cnper,0.d0)
            write(*,*)'X_e,Y_e ',X_e,Y_e
            write(*,*)'cnpar',cnpar
            write(*,*)'cnprp',cnprp
            call Disp_Ram(T_e,cnpar,X_e,Y_e,cnprp,K,d) !K is in z-y plane
                               !in Stix coordinates
            write(*,*)'prep3d_xyz K',K 
            call herm(K,K_herm)
            write(*,*)'prep3d_xyz K_herm',K_herm
            write(*,*)'D ',D 
c           end print relativistic tensor for testing
c-------------------------------------------------------------

         call muller(fr_func_noeps,errabs,nsig,nknown,nguess,nnew,
     +     roots,itmax,info,ier)

         cnprim=abs(imag(roots(1)))
         cnper=abs(dble(roots(1)))

         cnprim_cl=0.d0
         cnprim_e=cnprim
         cnprim_i=0.d0

C******* To be consistent with all other methods of calculating
C******* force cnprim and cnper to be positive.
         print *,'&&&&&&&&&&&& cnprim: ',cnprim,' cnper: ',cnper

         cnprp=dcmplx(cnper,cnprim)

c-------------------------------------------------------------
cSm060719
c--------calculate dielectric tensor reps for electric field calculations
         call Disp_combined(T_e,cnpar,X_e,Y_e,cnprp,reps,d)

         call efield1(cnpar,cnprp,ex,ey,ez,eplus,eminus,epar)
      endif   !iabsorp=12


C----------------------------- END OF IABSORP MODULES -----------------

c--------------------------------------------------------------------
c        put dielectric tensor reps into w_ceps 
c--------------------------------------------------------------------
         do i=1,3
           do j=1,3
               w_ceps(i,j,is)=reps(i,j) !Put the tensor reps for writing
           enddo
         enddo
    
      is=nrayelt
      wye(is)=wcw(x,y,z,1)
      if(nbulk.gt.1) wyi(is)=wcw(x,y,z,2)

      wxe(is)=wpw_2(x,y,z,1)
      if(nbulk.gt.1) wxi(is)=wpw_2(x,y,z,2)

      if (nbulk.gt.3) then
         wyi2(is)=wcw(x,y,z,4)
         wxi2(is)=wpw_2(x,y,z,4)
      endif
      r=dsqrt(x*x+y*y)

c---------------------------------------------------------
c     ws (cm)
      if (is.eq.1) then
         ws(1)=0.d0  ! ray distance 
         rhoold=rho ! rho is found at beginning of prep3d_xyz
         xold=0.d0
         yold=0.d0
         zold=0.d0
         i_ox_conversion=0
      else
         rold= sqrt(xold*xold+yold*yold)
         !delws=dsqrt((z-zold)**2+(r-rold)**2)*r0x*100. ! dl_pol
         delws=dsqrt((z-zold)**2+(x-xold)**2+(y-yold)**2)*r0x*100. !dl
         !YuP[09-26-2014]: Now (ws(is)-ws(is-1))) is dl along ray.
         ws(is)= ws(is-1) +delws
      end if

cyup      psi_s=psi_rho(rho)
cyup      q_s=qsafety_psi(psi_s) ! yup: just to print?
cyup      wphi(is)=u(3)
              call xyz_to_zrphi(x,y, r,phi) ! remove later
              wphi(is)=phi ! remove later

      call prepebw_xyz(u,is)
           
      xold=x
      yold=y
      zold=z
c--------------------------------------------------------
c     cflown - dimensionless for E_x/E,E_y/E,E_z/E
c  !!!!now flown is calculated using Mazzucato dielectric tensor
c  !!!!and electric field was calculated using Mazzucato tensor
c  !!!!it is only for EC wave .For LH and FW it is necessery
c !!!!!to create new subroutine flown ?what tensor?

      if (iflux.eq.1) then
          call flown_xyz(u(1),u(2),u(3),u(4),u(5),u(6),cflown)
      endif

      if (iflux.eq.2) then
c------- flux from the cold plasma    
         xe=wpw_2(x,y,z,1)
         ye=wcw(x,y,z,1)
         rnpar=cnpar
         rnper=cnper
         romega=1/ye
         romegpe=dsqrt(xe)/ye
         nsigma=ioxm
c        nsigma=1
c        nsigma=-1
         call grpde2(rnpar, rnper, romega, romegpe, nsigma,
     .                   rvgrpdc, redenfac, exde, eyde, ezde)
c        write(*,*)'+ redenfac',redenfac
         cflown=2.d0*redenfac
      endif !cold electron plasma flux
    
c-------------------------------------------------------

c-----calculate the group velocity
      id_old=id

      i_geom_optic_loc=i_geom_optic
      i_geom_optic=1 !to get group velocity in deru 
      call rside_xyz(u,deru)
      i_geom_optic=i_geom_optic_loc

      vgr(1)=deru(1) ! Vgroup/c
      vgr(2)=deru(2)
      vgr(3)=deru(3)
      
      vgrmods= vgr(1)**2 + vgr(2)**2 + vgr(3)**2
      if(vgrmods.gt.1.1d0) then
         write(*,*)
         write(*,*) '*************************************************'
         write(*,*)'WARNING vgroup/c>1,   |vgroup/c|=',dsqrt(vgrmods)  
         write(*,*) '*************************************************'
         write(*,*)
         if(vgrmods.gt. 2.) then
           write(*,*)'  iraystop->1'
           iraystop=1
           nrayelt=nrayelt-1 ! YuP[12-2016] Remove the last element - can be a "jump"
         endif
      endif

c-----the data for mnemonic.nc output file
c     the group velocity normalized to c
      wvgr_x(is)=vgr(1)
      wvgr_y(is)=vgr(2)
      wvgr_z(is)=vgr(3)
c-----refractive index
      wn_x(is)=cnx
      wn_y(is)=cny
      wn_z(is)=cnz
      
        call xyz_to_zrphi(x,y, r,phi) ! remove later
        call nxyz_to_nzrphi(cnx,cny,r,phi, cnr,cm) ! remove later
        wn_r(is)=cnr ! remove later
        wn_phi(is)=cm/r ! remove later
        
      ! Magnitude of group velocity, squared:
      vgrs=vgr(1)**2+vgr(2)**2+vgr(3)**2 !(Vgr,x)^2+(Vgr,y)^2+(Vgr,z)^2
      
c-----------------------------------------------------------
c     vdotb -projection of the group velocity  on the
c            magnetic field multiplited by bmod
c-----------------------------------------------------------
      vdotb= vgr(1)*bx + vgr(2)*by + vgr(3)*bz ! (Vgr.B) as vectors
c-----------------------------------------------------------
c     vgperps -perpendicular (to magnetic field)
c              component of the group velocity, squared
c-----------------------------------------------------------
c      vgperps=0.0d0
c      do i=1,3
c        vgperps=vgperps+(vgr(i)-vdotb*bf(i)*o_bmod**2)**2
c      enddo
      ! The above is same as 
      vgperps= vgrs - (vdotb*o_bmod)**2
c----------------------------------------------------------
c     collisional damping 
c----------------------------------------------------------
      if( iabsorp_collisional.eq.1) then
         temp_e=tempe_xyz(x,y,z,1)
         !dens_e=dense_xyz(x,y,z,1) ! rho and dens_e were found above
         z_eff=zeffrho(rho) 
         frqncy_l=frqncy
         v_gr_perp= dsqrt(vgperps)
         call absorp_collisional(temp_e,dens_e,frqncy_l,z_eff,
     &   v_gr_perp,coll_mult,
     &   cnprim_cl)
 2       cnprim=cnprim+cnprim_cl
      endif !iabsorp_collisional=1 
c----------------------------------------------------------
      !ckvi=dsqrt(vgperps/vgrmods)*cnprim/cld ! not used?
      
      ! needs work for a non-axisymmetric case:
      vgrrad= (vgr(1)*x + vgr(2)*y)/r ! radial (R) group vel.
      vgrpls= vgrrad**2 + vgr(3)**2   ! (Vgr,r)^2 + (Vgr,z)^2
      vgrpol= dsqrt(vgrpls) ! group vel. in poloidal cross-section
      
      wf=frqncy
c----------------------------------------------------------
c     ckvipol  (1/cm)
cyup      vratio=dsqrt(vgperps/vgrpls) !(Vgr_perp/Vgr_pol)
cyup      ckvipol=vratio*cnprim/cld !(Vgr_perp/Vgr_pol)*Im(n_perp)*omega/c [1/cm]
cyup      ckvipl_e=vratio*cnprim_e/cld
cyup      ckvipl_i=vratio*cnprim_i/cld
cyup      ckvipl_cl=vratio*cnprim_cl/cld
      ! With the above definitions, damping exponent will be defined as
      ! pwexp= exp[ -2*(Vgr_perp/Vgr_pol)*(Im(n_perp)*omega/c) *dl_pol ]
      
      ! But here is an alternative version:
      vratio=dsqrt(vgperps/vgrs) !(Vgr_perp/Vgr)
      ckvipol=  vratio*cnprim/cld !(Vgr_perp/Vgr)*Im(n_perp)*omega/c [1/cm]
      ckvipl_e= vratio*cnprim_e/cld
      ckvipl_i= vratio*cnprim_i/cld
      ckvipl_cl=vratio*cnprim_cl/cld
      ! With these definitions, damping exponent is defined as
      ! pwexp= exp[ -2*(Vgr_perp/Vgr)*(Im(n_perp)*omega/c) *dl ]
      ! where dl is along ray now (rather than a poloidal projection)  

      !Checking that none of cnprim is negative:
      if(ckvipol.lt.0.d0    .or.
     +   ckvipl_e.lt.0.d0   .or.
     +   ckvipl_i.lt.0.d0   .or.
     +   ckvipl_cl.lt.0.d0       )then
         write(*,*)'prep3d_xyz: cnprim,cnprim_e,cnprim_i,cnprim_cl=',
     +                          cnprim,cnprim_e,cnprim_i,cnprim_cl
         pause
      endif
      
 
cBH041009  Only germaine if iabsorp.eq.3:
      do kk=2,nbulk
         ckvipl_s(kk)=vratio*cnprim_s(kk)/cld
      enddo
c--------------------------------------------------------
c----------------------------------------------------------
      seikon(is)=0.d0
      spsi(is)=rho/a
c---------------------------------------------------------
c     wx,wy,wz (cm)  Coordinates along ray trajectory
      wx(is)=x*r00
      wy(is)=y*r00
      wz(is)=z*r00
      wr(is)=dsqrt(wx(is)**2 +wy(is)**2) 
      
      call xyz_to_zrphi(x,y, r,phi) ! remove later
      wphi(is)=phi ! remove later

      wnpar(is)=cnpar
      wnper(is)=cnper

      wmtor(is)=(-cnx*y+cny*x)*r00 ! cm= r*cnphi ! needs work? only for xdraw

c-------------------------------------------------------------------
c     Here delpwr(is) is the power(erg/c) in the ray
c     channel.It is equal powini_iray at antenna.
      if (is.eq.1) then
c        powj(iray) (erg/c) was calculated in cone_ec
c        powinilh(iray) (erg/c) was calculated in grill_lh
c        powini=powj or powinilh
         p=powini
cSAP081202
c        delpwr(is)=dexp(-2.d0*ckvipol*(ws(is)))*p
         delpwr(is)=p
         write(*,*)'is=1 delpwr(1)',delpwr(1)
         optical_depth=0.d0
cSm021101
         iflref_old=0

cSAP081111
         ckvipold=ckvipol
      else
        if (iabsorp_ql.eq.0) then
c----------do not use QL flux for absorption 
cSmirnov970105
c          argexp=(-2.d0*ckvipol*(ws(is)-ws(is-1)))
cSAP081111
c	   ckvipold=0.5d0*(sdpwr(is-1)+salphal(is-1)+salphac(is-1))
c           write(*,*)'ckvipol,ckvipold',ckvipol,ckvipold
c           write(*,*)'0.5d0*(sdpwr(is-1)+salphal(is-1)+salphac(is-1))',
c     &               0.5d0*(sdpwr(is-1)+salphal(is-1)+salphac(is-1))

           argexp=(-(ckvipol+ckvipold)*(ws(is)-ws(is-1))) 
           !YuP[09-26-2014]: Now (ws(is)-ws(is-1))) is dl along ray.
           ! Originally: it was dl_pol
cSAPO81111
           ckvipold=ckvipol
cSmirnov970105
           if (dabs(argexp).gt.90.d0) then ! argexp < -90.
             ! strong damping: all power is damped in one step
             delpwr(is)=0.d0
             pwexp=0.d0
           else 
             pwexp=dexp(argexp)
             if (pwexp.lt.1.d-50) then
                ! strong damping: all power is damped in one step
                 delpwr(is)=0.d0
             endif
           endif
           optical_depth=optical_depth+dabs(argexp)
c          write(*,*)'optical_depth',optical_depth
                  
cSm050225
c-YuP-130605: changed ray-stopping criterion from argexp>0.d0 to this: 
           if(argexp.gt.1.d-30)then
              write(*,*)'******************************************'
              write(*,*)'WARNING in prep3d_xyz argexp>0' 
              write(*,*)'It will give the growing ray power'
              write(*,*)'Stopping the ray:  iraystop->1'
              write(*,*)'******************************************'
              argexp=0.d0
              pwexp=1.d0
              iraystop=1
              nrayelt=nrayelt-1 ! YuP[12-2016] Remove the last element - can be a "jump"
              !pause !!!
           endif
c-YuP-130605: From print-out: 
c Even though sometimes ckvipol and ckvipold both are zero,
c yet the value of  argexp= -(ckvipol+ckvipold)*(ws(is)-ws(is-1))
c is not zero (but rather a small number ~ 1e-321).
c Because of this seemingly insignificant rounding error,
c the rays were stopped prematurely.
c It only happens on Hopper/PGI, not on IntelVisualFortran. 
cSm021101
c          write(*,*)'===delpwr(is-1),argexp',delpwr(is-1),argexp
           delpwr(is)=delpwr(is-1)*pwexp
c          delpwr(is)=delpwr(is-1)*dexp(-2.d0*dabs(wi)*
c     &               (ws(is)-ws(is-1))/(dsqrt(vgrs)*cvac)*    
c     &              (2.d0*pi*frqncy*1.d+9))
         else
c----------to use QL flux for absorption calculations
c          iabsorp_ql=1

           call absorbed_power_using_ql_flux_xyz
     +      ( wnpar(is-1),wnper(is-1),
     &     wx(is-1),wy(is-1),wz(is-1),
     &     fluxn(is-1),delpwr(is-1),(ws(is)-ws(is-1)),
     &     absorbed_power_ql )
           !YuP[09-26-2014] Need to check (ws(is)-ws(is-1)) input here:
           ! it was meant to be dl_pol originally.
     
           delpwr(is)=delpwr(is-1)-absorbed_power_ql

c           write(*,*)'QL absorption'
c           write(*,*)'delpwr(is-1),absorbed_power_ql,delpwr(is)',
c     &                delpwr(is-1),absorbed_power_ql,delpwr(is)

         endif !iabsorp_ql  
c-------------------------------------------------------------------
c        reflection lost at the plasma edge
c------------------------------------------------------------------
c         write(*,*)'prep3d_xyz refl_loss,irefl,irefl_old',
c     &             refl_loss,irefl,irefl_old
         tot_pow_absorb_at_refl=tot_pow_absorb_at_refl+
     &           delpwr(is)*refl_loss*(irefl-irefl_old)

c         write(*,*)'prep3d_xyz before refl_looss delpwr(i)',delpwr(is)

         delpwr(is)=delpwr(is)*(1.d0-refl_loss*(irefl-irefl_old))
         irefl_old=irefl 
         w_tot_pow_absorb_at_refl=tot_pow_absorb_at_refl
c-------------------------------------------------------------------    
          
cSAP090603
c         write(*,*)'prep3d_xyz delpwr(is)',delpwr(is)
c         write(*,*)'delpwr(1)',delpwr(1)

cSAP090603
cyup         if((delpwr(is).gt.1.d-200).and.(delpwr(1).gt.1.d-200))then
cyup           if(i_ox.ne.1) write(*,*)'-dlog(delpwr(is)/delpwr(1))',
cyup     &             -dlog(delpwr(is)/delpwr(1))
cyup         endif
 
c         if(argexp.gt.0.d0)then
c           write(*,*)'******************************************'
c           write(*,*)'WARNING in prep3d_xyz argexp>0' 
c           write(*,*)'It will give the growing ray power'
c          write(*,*)'******************************************'
c         endif
           
         if(delpwr(is).lt.delpwrmn*delpwr(1))then
            write(*,*)'in prep3d_xyz delpwr(is).lt.delpwrmn*delpwr(1)**'
c           stop ray_iray calculations
            write(*,*)'  iraystop->1'
            iraystop=1
         endif
      end if


      if(i_ox.eq.2) then
c---------------------------------------------------------------
c       It works for i_ox=2 case after OX mode conversion point,
c       where i_ox_conversion=1
c       It will reduce the power from O mode to X mode using 
c       transmission coefficient transm_ox

        if (i_call_prep3d_in_output_at_i_ox_conversion_eq_1.eq.0)goto 20
        ! No OX jump yet (not reached yet)

        if (i_call_prep3d_in_output_at_i_ox_conversion_eq_1.eq.1) then
        !The FIRST call of prep3d_xyz in outpt_xyz after OX conversion 
          nrayelt=nrayelt-1
          is=is-1 ! last point with O-mode 
          nrayelt_o_cutoff=is  !ray element where O cutoff was found
          write(*,'(a,i5)')'prep3d_xyz: nrayelt_o_cutoff=',
     +                                  nrayelt_o_cutoff
          goto 20 
        endif

        if (i_call_prep3d_in_output_at_i_ox_conversion_eq_1.eq.2) then
        !The SECOND call of prep3d_xyz in outpt_xyz after OX conversion
          transm_ox_loc=transm_ox
          ! Set the ray power at first point of X-mode:
          delpwr_o=delpwr(is-1) !last point of O-mode (just before jump)
          delpwr_x=delpwr_o*transm_ox
          delpwr(is)=delpwr_x   !first point of X-mode (just after jump)
          write(*,'(a,i5,2e13.4)')'prep3d_xyz is,transm_ox,delpwr(is)=',
     +                                        is,transm_ox,delpwr(is)
        endif
 20   continue
c---------------------------------------------------------------
      endif ! i_ox.eq.2


c     sdpwr(is)=0.d0
c----------------------------------------------------------------------
cSmirnov961122
      if(istart.eq.1) then
c        electron cyclotron waves
         wdnpar(is)=dabs(0.05d0*wnpar(1))
      else
cSmirnov961205
c        grill conditions for the waves (LH or FW)
         wdnpar(is)=wdnpar0(iray)
      endif

      cwexde(is)=cex
      cweyde(is)=cey
      cwezde(is)=cez
c---------------------------------------------------------
c     vgrpol(cm/sec)=vgrpol*r00/t00
c     vrgpol/c=vgrpol/wf
c---------------------------------------

c---------------------------------------
      wf=frqncy
      p_flux=cflown*dconjg(cflown)    
      fluxn(is)=dsqrt(p_flux)*vgrpol*0.5d0

c----------------------------------------------------------------------
c   Stop ray if flux.lt.0.0 (set fluxn previous value)  [RWH:030428]
c----------------------------------------------------------------------
         if(fluxn(is).lt.0.0)then
            write(*,*)
            write(*,*) '***********************************************'
            write(*,*) 'in prep3d_xyz fluxn(is).lt.0.0. Set iraystop=1'
            write(*,*) 'in prep3d_xyz Set fluxn(is)=fluxn(is-1)'
            write(*,*) '***********************************************'
            write(*,*)
            if(is.gt.1) then
               fluxn(is)=fluxn(is-1)
            else
               fluxn(is)=1.
            endif
c           stop ray_iray calculations
            iraystop=1
           nrayelt=nrayelt-1 ! YuP[12-2016] Remove the last element - can be a "jump"
         endif
     
c-----------------------
      one=1d0
      !feqd can be 0 in model_b.ne.0
      if( abs(feqd(1)).lt.1.d-33 ) then ! ~ zero
        sbtot(is)=bmod*10000.d0*b0
      else
        sbtot(is)=bmod*10000.d0*b0*sign(one,feqd(1))
      endif
cBH040915:  Magnetic field components
      sb_x(is)=1.e4*bx
      sb_y(is)=1.e4*by
      sb_z(is)=1.e4*bz
       sb_r(is)=1.e4*br 
       sb_phi(is)=1.e4*bphi ! remove later
c-----------------------------------------------------------
c     dense - dimensionless
      sene(is)=dense_xyz(x,y,z,1)*1.0d+13
      ste(is)=tempe_xyz(x,y,z,1)
c     salphac(is)=0.0d0
c     salphal(is)=2.d0*ckvipol
c Smirnov970105 beg
c BH991017   sdpwr(is)=2.d0*ckvipl_e   ! electron damping coefficient
c BH991017   salphal(is)=2.d0*ckvipl_i  ! ion damping coefficient
      salphac(is)=2.d0*ckvipl_cl  ! collisional damping coefficient
c Smirnov970105 end
cHarvey991017 beg
      sdpwr(is)=  2.d0*ckvipl_i     ! ion damping coefficient
      salphal(is)=2.d0*ckvipl_e     ! electron damping coefficient
cBH041009
      do kk=1,nbulk  !2,nbulk YuP[2016] include e why not?
         salphas(is,kk)=2.d0*ckvipl_s(kk)
         ! Find the thermal speed for each species:
         ! vth=sqrt(2T/m),    vth in (cm/sec),T(keV)
         ! vth= 1.87d9*dsqrt(T(i)/dmas(i))
         wvthermal(is,kk)=1.87d9*dsqrt(tempe_xyz(x,y,z,kk)/dmas(kk))
         ! cm/sec !
      enddo
cHarvey991017 end
c------------------------------------------------------------
c     for xdraw plotter
      xarr(is)=x*r00
      yarr(is)=y*r00
      rez(is)=cdabs(cez)
c     if (istart.eq.2) then
c        start point is inside the plasma (for lh and fw)
c        using cold plasma dielectric tensor
c        rez(is)=ezcold
c     end if
c------------------------------------------------------------
c-----------------------------------------------------------
c       data for onetwo
c-----------------------------------------------------------
      if(ionetwo.eq.1) then
cyup        if(is.gt.1) then       
cyup           call check_monotonic_radius_in_ray_element(is,
cyup     &     i_non_monotonic,z_center,r_center,rho_center)
           ! YuP: output: i_non_monotonic Not used?
cyup        endif
        call p_c_prof_xyz(is,rhoold,rho,cnpar,cnper,cex,cey,cez)
      !Note: rho is found at beginning of prep3d_xyz, by dense_xyz()
      endif

      rhoold=rho

      return
      end ! prep3d_xyz
      
      
      
c======================================================================
c======================================================================

      real*8 function psilim_psif_xyz(p)
c---------------------------------------------------------------
c     It is used to find the itersection of the sraight line
c     with the plasma boundary. 
c     It calculates the difference: psilim-psi(x,y,z)
c     Here x=x(p),y=y(p),z=z(p) are the coordinates of the point  
c     along the line:
c     xp= xma + (xst(i_cone)-xma)*p
c     yp= yma + (yst(i_cone)-yma)*p
c     zp= zma + (zst(i_cone)-zma)*p
c
c     i_cone is the index of the EC cone
c           It is in antenna.i common block
c---------------------------------------------------------------       
c      implicit real*8 (a-h,o-z)
      implicit none
c-----input 
      real*8 p ! the parameter which determines the point
             ! on the straight line            
      include 'param.i'     
      include 'three.i'   ! rma,zma,psilim
      include 'cone.i'    ! xst,yst,zst (ncone) 
      include 'antenna.i' ! gives i_cone the number
c-----locals 
      real*8 rp, xp,yp,zp
c-----external
      real*8  psif_xyz 

      xp= xma + (xst(i_cone)-xma)*p
      yp= yma + (yst(i_cone)-yma)*p
      zp= zma + (zst(i_cone)-zma)*p
      psilim_psif_xyz= psilim-psif_xyz(xp,yp,zp) 
      return
      end


c======================================================================
c======================================================================

c     The subroutines for the optimaL OXB launch.
      subroutine antenna_surface_xyz(x,y,z, x_ant,y_ant,z_ant)
c     determine the possible antenna position 

c     x_ant(theta), y_ant(theta), z_ant(theta), r_ant(theta) [m]
  
c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none

      include 'param.i'                !for i_ant=2 case
      include 'three.i' ! rma,zma,xma,yma, etc.
c----------------------------------------
      include 'five.i' !zmax,rzmax     !
      include 'cone.i' !xst(),yst(),zst(),rst()    !
c------------------------------------------
      include 'antenna.i'
  
      logical first !used in i_ant=1 case
      data first / .true./
      save first

c-----input
      real*8 x,y,z
c-----output
      real*8 x_ant, y_ant, z_ant
c-----local
      real*8 p,p1,p2,det, det_a,det_b,a_sq,b_sq,
     & accuracy,p_intersection,
     +r_intersection, x_intersection, y_intersection, z_intersection,
     &distance,theta_l,pi, r, r_ant, theta, ela, rax

c-----externals
      real*8 rtbis,psilim_psif_xyz
      external psilim_psif_xyz
 
      pi=4.d0*datan(1.d0)
      
      if((i_ant.lt.1).or.(i_ant.gt.2)) then
        write(*,*)'antenna_surface_xyz'
        write(*,*)'(i_ant.lt.1).or.(i_ant.gt.2)' 
        write(*,*)'i_ant=',i_ant
        write(*,*)'Please check i_ant'
        stop 'in antenna_surface_xyz'
      endif

      r= dsqrt(x*x+y*y)
      xma= rma*x/r
      yma= rma*y/r !ok for now? Assume symm.in phi. Otherwise needs work
      theta= atan2(z-zma,r-rma) ! pol.angle in range [-pi,pi] ! ok for now?
      
      if(i_ant.eq.2) then  ! revised from original genray !
        ela=1.d0 
        ! Specify elongation ela==b/a of an ellipse (r/a)^2 + (z/b)^2 =1
        ! at which the starting point rst(icone),zst(icone) is located.
        ! The point r_ant,z_ant will be looked at this ellipse.
        ! No need for accuracy; 
        ! just to match approximately the flux surface shape.
c------------------------------------------------------------------------
c       EC cone vertex is at the ellipse curve: 
c
c       (r_ant-rma)**2     (z_ant-zma)**2 
c       --------------  +   -------------  =1
c            a**2               b**2
c
c       Assume the following points are at the ellipse curve:
c       (rst,zst) the given EC cone vertex,  and  (r_ant,z_ant) 
c
c       Find the value of 'a' from 
c
c       (rst-rma)**2         (zst-zma)**2 
c       --------------  +   -------------  =1
c            a**2              (a*ela)**2
c
c        a= sqrt( (rst-rma)**2  + (zst-zma)**2/ela**2 )
c
c        Then,
c        r_ant==r(theta)= rma + a*cos(theta)
c        z_ant==z(theta)= zma + b*sin(theta), where b==a*ela
c-------------------------------------------------------------------------
c        i_cone is the number of the EC cone i_cone=1,...,nconea
        rax= dsqrt( (rst(i_cone)-rma)**2  +(zst(i_cone)-zma)**2/ela**2 )
        r_ant= rma + rax*dcos(theta)
        z_ant= zma + rax*ela*dsin(theta)
        x_ant= r_ant*x/r ! Assume symm.in phi. Otherwise needs work
        y_ant= r_ant*y/r ! Assume symm.in phi. Otherwise needs work
c        write(*,*)'antenna_surface: rst,zst=',rst(i_cone),zst(i_cone)
c        write(*,*)'antenna_surface: r_ant,z_ant=',r_ant,z_ant
      endif  !i_ant.eq.2

      if (i_ant.eq.1) then
c---------------------------------------------------------------------
c        The antenna surface will be at the same given distance from
c        the plasma edge. This distances is along the straight line 
c        from the magnetic axis. The distance is equal to the
c        distance of the given EC cone vertex from the plasma edge.
c-------------------------------------------------------------------
c         r-rma      z-zma
c        ------- = --------=p is the line through the magnetic axis (rma,zma)
c        rst-rma    zst-zma   and the EC cone vertex (rst,zst)
c
c        r=rma + (rst-rma)p 
c        z=zma + (zst-zma)p
c
c        or   z=zma + (r-rma)*(zst-zma)/(rst-rma) 
c
c        Plasma boundary: rmin < r < rmax
c
c        The point of this line intersection with the plasma boundary
c
c        psilim-psif_xyz(x,y,z_intersecsion) =0 
c
c        Here z_intersection=rma+(rst-rma)p_intersection
c             r_intersection=zma+(zst-zma)p_intersection
c---------------------------------------------------------------------
c         write(*,*)'antenna_surface first=',first
        
         if (first) then
c-----------create the antenna surface          
            do iray=1,ncone
c--------------the loop by all EC cone vertexs
                write(*,*)'antenna_surface iray=',iray
c-----------------------------------------------------
c              calculate intersection [r_intersection,z_intersection]
c              of plasma boundary with the stright line
c              from the magnetic axis to the given EC cone vertex
c------------------------------------------------------
c               write(*,*)'antenna surface before p_intersection=rtbis'
               i_cone=iray !set this parameter to antenna.i
               accuracy=1.d-7
               p1=0.d0
               p2=1.d0
               p_intersection= rtbis(psilim_psif_xyz,p1,p2,accuracy)
               x_intersection= xma+(xst(i_cone)-xma)*p_intersection
               y_intersection= yma+(yst(i_cone)-yma)*p_intersection
               z_intersection= zma+(zst(i_cone)-zma)*p_intersection
               r_intersection=dsqrt(x_intersection**2+y_intersection**2)
               distance=dsqrt((xst(i_cone)-x_intersection)**2+
     +                        (yst(i_cone)-y_intersection)**2+
     +                        (zst(i_cone)-z_intersection)**2 )
c------------------------------------------------------------------
c              calculate spline coefficients for the antenna points
c              at the surface which is at the equal distance from
c              the plasma {z_an(theta),r_ant(theta)}
c------------------------------------------------------------------  
c               write(*,*)'antenna surface before antenna_equal_dist'
c               write(*,*)'first,i_cone,distance,theta,z_ant,r_ant',
c     &         first,i_cone,distance,theta,z_ant,r_ant
cyup               theta_l=theta
               call antenna_equal_dist_xyz(first,i_cone,distance,theta,
     &                                 z_ant,r_ant)
               x_ant= r_ant*x/r ! needs work
               y_ant= r_ant*y/r ! needs work
               write(*,*)'********************************'
               write(*,*)'antenna surface was created antenna surface'
               write(*,*)'********************************'
            enddo !iray   
             
            first=.false.
         endif !first
c--------calculate the antenna vertex coordinates z_ant,r_ant
c        at the antenna surface for the given poloidal angle theta
c        using spline

c         write(*,*)'antenna surface before antenna_equal_dist'
c         write(*,*)'******************************************'

c         write(*,*)'first,i_cone,distance,theta',
c     &         first,i_cone,distance,theta 
cyup         theta_l=theta
         call antenna_equal_dist_xyz(first,i_cone,distance,theta,
     &                           z_ant,r_ant)
         x_ant= r_ant*x/r ! needs work
         y_ant= r_ant*y/r ! needs work
c         write(*,*)'antenna surface after antenna_equal_dist'
c         write(*,*)'theta,z_ant,r_ant',theta,z_ant,r_ant
      endif  !i_ant.eq.1

      return      
      end


c======================================================================
c======================================================================



      subroutine antenna_equal_dist_xyz(first,i_cone,distance,theta,
     & z_ant,r_ant)
c-----------------------------------------------------------
      implicit none 
      include 'param.i' 
      include 'three.i'
      include 'cone.i'
      
c-----input      
      logical first !
      integer i_cone !the index of a given EC cone 
      real*8 theta, ! poloidal angle (radians)
     &distance      ! the distance from the plasma to antenna
                    ! along the sraight line 
                    ! xp=xma + (xst(icone)-xma)*p 
                    ! yp=yma + (yst(icone)-yma)*p 
                    ! zp=zma + (zst(icone)-zma)*p      
c-----output
      real*8 r_ant,z_ant ! coordinates of the antenna vertex at
                         ! the given poloidal angle theta
                         ! for the given i_cone 
 
c-----local
      real*8 theta_ar(nteta),
     &z_ant_ar(nconea,nteta),r_ant_ar(nconea,nteta),
     &z_bound(nteta),r_bound(nteta),tabl(3),
     &temp_function(nteta),temp_deriv(nteta),pi,dth_pol,rho_bound,p,
     &workk(3*nteta+1),two_pi,theta_loc
      integer iop(2),itabl(3),i

      real*8 work(3*nteta+1),d2z_ant(nconea,nteta),d2r_ant(nconea,nteta)
c-----externals period_argument
      real*8 period_argument
      
      save d2z_ant,d2r_ant,theta_ar,r_ant_ar,z_ant_ar
    
      pi=4.d0*datan(1.d0) 

c      write(*,*)'antenna_equal_dist first,i_cone,distance,theta',
c     &first,i_cone,distance,theta

      dth_pol=2.d0*pi/(nteta-1.d0)
c      write(*,*)'antenna_equal_dist nteta,i_cone',nteta,i_cone
      if (first) then
c-------------------------------------------------------------------
c       creation of spline coefficients for the antenna surface
c-------------------------------------------------------------------
        do i=1,nteta-1
          theta_ar(i)=(i-1)*dth_pol
c          write(*,*)'antenna_equal_dist i,theta_ar(i)',i,theta_ar(i)
c---------creation mesh [z_bound,r_bound] along plasma boundary 
          !-> Get z(psilim,theta_ar),r(psilim,theta_ar)
          call zr_psith(psilim,theta_ar(i),z_bound(i),r_bound(i))
          rho_bound=dsqrt((z_bound(i)-zma)**2+(r_bound(i)-rma)**2)
          p=1.d0+distance/rho_bound
c---------creation mesh [z_ant_ar,r_ant_ar] along antenna surface
          z_ant_ar(i_cone,i)=z_bound(i)+(p-1.d0)*(z_bound(i)-zma)
          r_ant_ar(i_cone,i)=r_bound(i)+(p-1.d0)*(r_bound(i)-rma)
c          write(*,*)'antenna_equal_dist z_ant_ar(i_cone,i),
c     &    r_ant_ar(i_cone,i)',z_ant_ar(i_cone,i),r_ant_ar(i_cone,i)
        enddo !i

        theta_ar(nteta)= 2.d0*pi
        z_ant_ar(i_cone,nteta)= z_ant_ar(i_cone,1)
        r_ant_ar(i_cone,nteta)= r_ant_ar(i_cone,1)

c        do i=1,nteta
c        write(*,*)'i,theta_ar(i),z_ant_ar(i_cone,i),r_ant_ar(i_cone,i)',
c     &  i,theta_ar(i),z_ant_ar(i_cone,i),r_ant_ar(i_cone,i)
c        enddo
c----------------------------------------------------------------
c     
        iop(1)=3 ! periodic spline boundary conditions
        iop(2)=3

c-------creation spline coefficients [d2z_ant,d2r_ant] for the antenna surface
        call take_1d_array_from_2d_array(temp_function,z_ant_ar,
     &  nconea,ncone,nteta,nteta,i_cone)
        call coeff1(nteta,theta_ar,temp_function,temp_deriv,
     &              iop,1,workk)
        call put_1d_array_to_2d_array(temp_deriv,d2z_ant,
     &  nconea,ncone,nteta,nteta,i_cone)

        call take_1d_array_from_2d_array(temp_function,r_ant_ar,
     &  nconea,ncone,nteta,nteta,i_cone)
        call coeff1(nteta,theta_ar,r_ant_ar(i_cone,1),d2r_ant(i_cone,1),
     &              iop,1,workk)
        call put_1d_array_to_2d_array(temp_deriv,d2r_ant,
     &  nconea,ncone,nteta,nteta,i_cone)
     
        write(*,*)'antenna_equal_dist after spline coeff1' 
       
      endif !first
     
      itabl(1)=1
      itabl(2)=0
      itabl(3)=0

      two_pi=8.d0*atan(1.d0)         
      theta_loc=period_argument(theta,two_pi)
     
      call take_1d_array_from_2d_array(temp_function,z_ant_ar,
     &nconea,ncone,nteta,nteta,i_cone)
      call take_1d_array_from_2d_array(temp_deriv,d2z_ant,
     &nconea,ncone,nteta,nteta,i_cone)
      call terp1(nteta,theta_ar,temp_function,temp_deriv,
     &           theta_loc,1,tabl,itabl)
      z_ant=tabl(1)

      call take_1d_array_from_2d_array(temp_function,r_ant_ar,
     &nconea,ncone,nteta,nteta,i_cone)
      call take_1d_array_from_2d_array(temp_deriv,d2r_ant,
     &nconea,ncone,nteta,nteta,i_cone)
      call terp1(nteta,theta_ar,temp_function,temp_deriv,
     &           theta_loc,1,tabl,itabl)
      r_ant=tabl(1)
c      write(*,*)'antenna_equal_dist theta,r_ant,z_ant',
c     &theta,r_ant,z_ant
      return
      end


c======================================================================
c======================================================================

    
      real*8 function f_antenna_min_xyz(p)
c----------------------------------------------------------------
c     Calculate the function:
c     f_antenna_min=       (r_ray(p)-r_ant(theta(p)))**2+
c                         +(z_ray(p)-z_ant(theta(p)))**2)
c
c     Here
c     x_ant, y_ant, z_ant are the coordinates of the antenna surface
c
c     x_ray(p), y_ray(p), z_ray(p) are the coordinates along 
c     straight line starting at the plasma edge x_0,y_0,z_0 with the
c     direction n_0_x, n_0_y, n_0_z
c     
c      implicit real*8 (a-h,o-z)
      implicit none

c-----input
      real*8 p !parameter along the ray

      real*8         x_0,y_0,z_0, n_vac_x,n_vac_y,n_vac_z
      common /ox_xyz/x_0,y_0,z_0, n_vac_x,n_vac_y,n_vac_z
     
      include 'param.i'
      include 'three.i' !rma,zma
      
      real*8 x_ray, y_ray, z_ray, x_ant, y_ant, z_ant

      x_ray=  x_0 + p * n_vac_x
      y_ray=  y_0 + p * n_vac_y
      z_ray=  z_0 + p * n_vac_z
c      write(*,*)'x_0,    z_0    =', x_0,     z_0
c      write(*,*)'x_ray,  z_ray  =', x_ray,   z_ray
      !-> get x_ant, y_ant, z_ant for this point (p) on the line:
      call antenna_surface_xyz(x_ray,y_ray,z_ray, x_ant,y_ant,z_ant)
      
      f_antenna_min_xyz= (x_ray-x_ant)**2 + 
     +                   (y_ray-y_ant)**2 +
     +                   (z_ray-z_ant)**2
      return
      end

c======================================================================
c======================================================================

      subroutine edg_vac_refr_index_xyz(n_x, n_y, n_z, 
     + dpsi_dx, dpsi_dy, dpsi_dz, !   input
     & n_vac_x, n_vac_y, n_vac_z) !-> output
c-----calculate vacuum refractive index at the plasma edge
c     n_vac_x, n_vac_y, n_vac_z
c     using the refractive index in plasma at the plasma edge
c     n_x, n_y, n_z 
c     and the poloidal flux gradient
c     dpsi_dx, dpsi_dy, dpsi_dz

      implicit none
c-----input
      real*8 n_x, n_y, n_z,  !refractive index at the plasma boundary 
                             !in plasma side
     &dpsi_dx, dpsi_dy, dpsi_dz  !gradiend of the poloidal flux        
                             !at the plasma boundary
c-----output
      real*8 n_vac_x, n_vac_y, n_vac_z  ! refractive index at the plasma
                                        ! boundary in vacuum side

c-----local
      real*8 p,
     &gradpsi_mod,
     +n_gradpsi_pl_x, n_gradpsi_pl_y,
     &n_gradpsi_pl_z,           ! refractive index along the vector
     &n_gradpsi_pl_mod,         ! grad_psi at plasma side

     &n_psi_x, n_psi_y, n_psi_z,          ! refractive index along 
     &n_psi_mod_s,                        ! the magnetic surface

     &n_gradpsi_vac_x,n_gradpsi_vac_y,n_gradpsi_vac_z, ! refr index along
     &n_gradpsi_vac_mod                   ! grad_psi at vacuum side

c-----refractive index along grad_psi at plasma side
      gradpsi_mod=  dsqrt(dpsi_dx**2 + dpsi_dy**2  + dpsi_dz**2)
      n_gradpsi_pl_mod= (n_x*dpsi_dx + n_y*dpsi_dy + n_z*dpsi_dz)
     /                              /gradpsi_mod
      n_gradpsi_pl_x= n_gradpsi_pl_mod*dpsi_dx/gradpsi_mod
      n_gradpsi_pl_y= n_gradpsi_pl_mod*dpsi_dy/gradpsi_mod
      n_gradpsi_pl_z= n_gradpsi_pl_mod*dpsi_dz/gradpsi_mod

c-----refractive index along the magnetic surface 
c     at plasma side and at vacuum side
      n_psi_x= n_x - n_gradpsi_pl_x
      n_psi_y= n_y - n_gradpsi_pl_y
      n_psi_z= n_z - n_gradpsi_pl_z
      n_psi_mod_s= n_psi_x**2 + n_psi_y**2 + n_psi_z**2 ! should be .le.1

c-----refractive index along grad_psi at vacuum side
      p=1.d0-n_psi_mod_s
      if (p.gt.0.d0)then
          n_gradpsi_vac_mod=dsqrt(p)
      else
          write(*,*)'*****WARNING****'
          write(*,*)'oxb  edg_vac_refr_index np_psi_mod_s>1'
          write(*,*)'n_psi_mod_s',n_psi_mod_s
          n_gradpsi_vac_mod=0.d0
c          stop 'oxb'
      endif

      n_gradpsi_vac_x= n_gradpsi_vac_mod*dpsi_dx/gradpsi_mod
      n_gradpsi_vac_y= n_gradpsi_vac_mod*dpsi_dy/gradpsi_mod
      n_gradpsi_vac_z= n_gradpsi_vac_mod*dpsi_dz/gradpsi_mod

c-----vacuum refractive index
      n_vac_x= n_psi_x + n_gradpsi_vac_x
      n_vac_y= n_psi_y + n_gradpsi_vac_y
      n_vac_z= n_psi_z + n_gradpsi_vac_z
      return
      end

c======================================================================
c======================================================================

      subroutine antenna_vertex_xyz(x_0,y_0,z_0, n0_x,n0_y,n0_z,
     & b_x0,b_y0,b_z0, dpsi_dx,dpsi_dy,dpsi_dz,
     & x_st,y_st,z_st, alpha_st,beta_st,icone)
c-----calculate O mode cone vertex
c     x_st,y_st,z_st
c     and the central ray direction
c     alpha_st,beta_st
c
c     At the intersection of the vacuum ray (straight line)
c     starting in the point (x_0,y_0,z_0) in the direction    
c     determined by vacuum refractive index
c
      implicit none

c-----input
      real*8 x_0,y_0,z_0,   !cartesian coords. of location of the O-mode
                            !at the plasma boundary [m]
     &n0_x,n0_y,n0_z,       !o-mode refractive index at the plasma boundary
                            ! in plasma side
     &b_x0,b_y0,b_z0,       !the magnetic field in  M_0(r_0,phi_0,z_0)
     &dpsi_dx,dpsi_dy,dpsi_dz    !gradiend of the poloidal flux        
                            !at the plasma boundary
      integer icone         !the number of EC cone vertex =1,...,ncone
                            !ncone =< nconea

c-----output
      real*8 x_st,y_st,z_st, !EC cone vertex [m, degree]
     &alpha_st,!central ray toroidal angle [degree]
     &beta_st !central ray poloidal angle
     
c-----externals
      real*8 x_bin_min,f_antenna_min_xyz
      external f_antenna_min_xyz

c-----locals
      real*8 n_vac_x, n_vac_y, n_vac_z,
     & pi,gamma,delta,trans, r_0, r_st, phi_0, phi_st
      double precision
     &pacc,p_antenna,p1,p2, p_range, delta2,cos_gamma
     
c-----------------------------------------------------------------
      real*8 r_min,r_max
c-------------------------------------------------------------------
      include 'param.i'    
      include 'antenna.i'   
      include 'three.i' ! rma,zma,xma,yma, etc.
      include 'five.i' !zmax,rzmax     !
      include 'cone.i' !xst(),yst(),zst(),rst()    !
c--------------------------------------------------------------------   
      real*8         x_00,y_00,z_00, nvac_x,nvac_y,nvac_z
      common /ox_xyz/x_00,y_00,z_00, nvac_x,nvac_y,nvac_z

      write(*,*)'ENTER antenna_vertex_xyz'
c      i_cone=1      ! the number of EC cone =1,...,nconea
      i_cone=icone  
      n_vac_x= n0_x ! initialize ! 
      n_vac_y= n0_y ! initialize ! 
      n_vac_z= n0_z ! initialize ! 
      
      if(i_ant.eq.1 .or. i_ant.eq.2) then
c        calculate (x,y,z) components of the vacuum refractive index
c        n_vac_x,n_vac_y,n_vac_z   at the plasma boundary: 
         call edg_vac_refr_index_xyz(n0_x,    n0_y,    n0_z,
     &                         dpsi_dx, dpsi_dy, dpsi_dz,  
     &                         n_vac_x, n_vac_y, n_vac_z) !-> out
         write(*,*)'n_vac_x,n_vac_z=', n_vac_x, n_vac_z
         x_00= x_0 ! put these values into common /ox/
         y_00= y_0 
         z_00= z_0 
         nvac_x= n_vac_x
         nvac_y= n_vac_y
         nvac_z= n_vac_z
c        calculate the vertex coordinates
         pacc=1.d-10 
         ! distance range for searching the surface with antenna,
         ! in proximity of the surface containing (x_0,y_0,z_0) point
         p_range= sqrt(rst(icone)**2+zst(icone)**2) 
     +          - sqrt(x_0**2+y_0**2+z_0**2)  
         p1= -p_range*2.5 ! extended lower limit
         p2=  p_range*2.5 ! extended upper limit
         write(*,*)'antenna_vertex before x_bin_min p1,p2,pacc',
     &   p1,p2,pacc
         !-------------------------------------------------
         p_antenna= x_bin_min(f_antenna_min_xyz,p1,p2,pacc)
         !-------------------------------------------------
         write(*,*)'after x_bin_min p_antenna',p_antenna
         x_st= x_0 + p_antenna*n_vac_x
         y_st= y_0 + p_antenna*n_vac_y 
         z_st= z_0 + p_antenna*n_vac_z
      endif ! i_ant.eq.1 .or. i_ant.eq.2

      if(i_ant.eq.3) then ! the ray was stopped at R=rst(icone) surface
         x_st= x_0 
         y_st= y_0  
         z_st= z_0 
      endif ! i_ant.eq.3

      r_st=dsqrt(x_st**2+y_st**2)
      r_0=dsqrt(x_0**2+y_0**2)

      phi_0 = atan2(y_0,x_0) ! [-pi; pi]
      phi_st= dacos(x_st/r_st)
      if (y_st.lt.0.d0) phi_st=-phi_st

c-----calculation the central ray direction

c--------poloidal angle beta_st calculation
      pi=4.d0*datan(1.d0)
      if(z_0.le.z_st) then
         if(r_0.le.r_st) then
            !  0=<beta_st<pi/2
            beta_st=dasin(n_vac_z)
         else  
            ! pi/2<beta_st<pi
            beta_st=pi-dasin(n_vac_z)
         endif
      else
         if(r_0.le.r_st) then
            !  -pi/2<beta_st=<0
            beta_st=dasin(n_vac_z)
         else  
            ! -pi<beta_st<-pi/2
            beta_st=-pi-dasin(n_vac_z)
         endif
      endif
      beta_st=-beta_st

c-------toroidal angle calculation alpha_st
      delta2= r_st**2+r_0**2-2.d0*r_st*r_0*dcos(phi_st-phi_0)
      if(delta2.gt.0.d0) then
         delta=dsqrt(delta2)
         cos_gamma= (r_st**2+delta2-r_0**2)/(2.d0*r_st*delta)
      else
         delta=0.d0
         cos_gamma= 1.d0 ! ok?
      endif
      if(abs(cos_gamma).lt.1d0)then
         gamma=dacos(cos_gamma)
      else
         if(cos_gamma.ge. 1.d0) gamma=0.d0
         if(cos_gamma.le.-1.d0) gamma=pi
      endif

      if (phi_st.ge.phi_0) then
        alpha_st=pi+gamma
      else
        alpha_st=pi-gamma
      endif 

c-----transformation from radians to degrees
      trans=180.d0/pi
      beta_st=  beta_st*trans
      alpha_st= alpha_st*trans
      phi_st=   phi_st*trans
      write(*,*)'EXIT antenna_vertex_xyz: r_st=',r_st
      return
      end


c======================================================================
c======================================================================

      subroutine find_maximal_OX_transmission_xyz(xe0,
     & x_in,y_in,z_in, nx_in,ny_in,nz_in,
     & i_ox_conversion)
c-----if((xe.gt.(xe0-eps_xe)).and.(xe.le.xe0)) then it
c     calculates OX transmission coefficient transm_ox in u() point
c     and compares this coefficient with the transmission coefficient
c     transm_ox_old at the previous point u_old()
c     Finds the point where the transmission coefficient has the maximal
c     value and put i_ox_conversion=1
      implicit none
c-----input
      real*8 xe0, x_in,y_in,z_in, nx_in,ny_in,nz_in
     
      include 'param.i' 
      include 'one.i' 
      include 'output.i'
      include 'oxb.i'
c-----output 
c     r_in,z_in,phi_in,nr_in,nz_in,m_in the ray coordinates at
c             the point with the maxima transmission coefficient
      integer i_ox_conversion !=0 the ray is far from 
                              !   the OX conversion point
                              !=1 the ray is near the OX
                              !   conversion point
c-----externals
      real*8 wpw_2,bxyz

c-----local
      real*8 xe,transm_ox,u(6),u_old(6),cnpar,cnpar2,cn2,cnper2
      integer i

      save  u_old

c      write(*,*)'--> find_maximal_OX_transmission_xyz x_in=',x_in

      u(1)= x_in
      u(2)= y_in
      u(3)= z_in
      u(4)= nx_in
      u(5)= ny_in 
      u(6)= nz_in

      bmod= bxyz(x_in,y_in,z_in)  
      cnpar=(bx*nx_in+by*ny_in+bz*nz_in)*o_bmod
      cnpar2=cnpar*cnpar
      cn2= nx_in*nx_in + ny_in*ny_in + nz_in*nz_in
      cnper2=dabs(cn2-cnpar2)
      
      xe= wpw_2(u(1),u(2),u(3),1)
      i_ox_conversion=0

cyup      write(*,*)'was_in_ox_vicinity',was_in_ox_vicinity

      !YuP[11-2016] some adjustments      
      if( (xe.gt.(xe0-eps_xe) .or. (cn2-cnpar2.le.0.01)) 
     +  .and. (xe.le.xe0) )then
cyup      if((xe.gt.(xe0-eps_xe)) .and. (xe.le.xe0))then
         !Case 2. Ray has approached cutoff [Xe>(Xe0-eps_xe)]
         !Find transm_ox, save it as transm_ox_old
         ! Continue moving along ray while transm_ox is improving.
         ! But if transm_ox starts to deteriorate, 
         ! change i_ox_conversion=0->1, which will initiate
         ! a jump across cutoff / search for X-mode after cutoff
      
         write(*,*)'Xe is close to xe0. BEFORE transmit_coef_ox_xyz:'
         write(*,'(a,3e13.4)')' xe, xe0-eps_xe, cn2-cnpar2', 
     +                          xe, xe0-eps_xe, cn2-cnpar2
         !pause !!!
         
         was_in_ox_vicinity=.true.

         call transmit_coef_ox_xyz(u(1),u(2),u(3),u(4),u(5),u(6),
     &                           transm_ox)
         write(*,*)'AFTER: transm_ox, transm_ox_old, Nper^2=',
     &    transm_ox, transm_ox_old, cn2-cnpar2

         if( ((transm_ox-0.001).le.transm_ox_old)
     +      .or. (cn2-cnpar2.le.0.d0) )then 
            !YuP[11-2016] added -0.001 accuracy
            !YuP[11-2016] added (cn2-cnpar2.le.0.d0) condition.
            !If transm_ox starts to deteriorate (with accuracy of 0.001), 
            !or Nperp^2 <= 0, then initiate a jump over cutoff
            !(usually happens in few steps after entering Case 2)
            write(*,*)'transm_ox<transm_ox_old, initiate OX jump'
            !pause
            i_ox_conversion=1 ! means: will be stepping across cutoff, 
                      ! until Xmode is detected (using find_rho_X_xyz)
            ! use values from the previous time step,
            ! redefine *_in point to the point where transm.coef is highest
            !(i.e. at the previous time step)
            x_in=  u_old(1)
            y_in=  u_old(2)
            z_in=  u_old(3)
            nx_in= u_old(4)
            ny_in= u_old(5) 
            nz_in= u_old(6)
            !pause !!!
            goto 10 !-> finish
         endif
         do i=1,6
            u_old(i)=u(i)
         enddo
         transm_ox_old=transm_ox ! transm_ox_old is in oxb.i
        write(*,*)'Continue moving to Xe=1 while transm_ox is improving'
         
      else
      
         !Case 1. Ray is still far from cutoff [Xe<Xe0 and Xe<(Xe0-eps_xe)]
         !Do nothing. i_ox_conversion remains 0 (meaning: do not initiate stepping
      
         if((xe.gt.xe0).and.was_in_ox_vicinity) then
            !Case 3. Ray has just passed the point Xe=Xe0.
            !It can happen if in Case 2 above 
            !the condition transm_ox<transm_ox_old was never reached,
            !meaning that the transm.coef. was only improving 
            !(increasing or flat) during approach to Xe=1. 
            !So, there was no max point for transm.coef.
            ! In this case, force the OX jump as soon as Xe>Xe0, here:
            i_ox_conversion=1 ! means: will be stepping across cutoff, 
                      ! until Xmode is detected (using find_rho_X_xyz)
            transm_ox=transm_ox_old 
            ! Use values from the previous time step 
            !(when Xe was just below Xe0==1, or maybe Xe was equal Xe0).
            ! Redefine *_in point to the point where transm.coef is highest
            !(i.e. at the previous time step):
            x_in=  u_old(1)
            y_in=  u_old(2)
            z_in=  u_old(3)
            nx_in= u_old(4)
            ny_in= u_old(5) 
            nz_in= u_old(6)
            write(*,*) 'find_maximal_ox: just past Xe=1.  Xe, Nper^2=',
     +        xe,cn2-cnpar2
            !pause !!!
            goto 10 !-> finish
          endif
          
      endif

 10   continue
cyup      write(*,*)'after 10 i_ox_conversion',i_ox_conversion
      return
      end     


c======================================================================
c======================================================================

      subroutine find_rho_X_xyz(
     & x_ini,y_ini,z_ini, cnx_ini,cny_ini,cnz_ini,
     & x_x,y_x,z_x, cnx_x,cny_x,cnz_x,cnpar_x,cnper_x,
     + i_ox_conversion)

c-----finds the point (x_x, y_x, z_x) where
c     X mode with (nx, ny, nz) has the cutoff.
      implicit none
      include 'param.i'
      include 'one.i'
      include 'output_no_nml.i' ! was_not_ox_conversion indicator for OX jump
      include 'three.i'
      include 'five.i'
      include 'oxb.i'
c-----input
      real*8 x_ini,y_ini,z_ini,r_ini ! coordinates before O cutoff
      real*8 cnx_ini             ! Nx before O cutoff
      real*8 cny_ini             ! Ny 
      real*8 cnz_ini             ! Nz 
c     real*8 rho is the small radius at (x_ini,y_ini,z_ini) point .
c     It is in common/one/ in one.i
c-----output
      real*8 x_x,y_x,z_x
      real*8 cnpar_x,cnper_x    ! N_parallel N_perp in X_cutoff point
      real*8 cnx_x,cny_x,cnz_x  ! refractive index coordinates of
                                ! X cutoff point      
c-----externals
      real*8 bxyz, wpw_2, dense_xyz
c-----local
      real*8 psi,hstep,ppp,gradpsi,den
      integer iraystop,id_loc,ioxm_loc, iter, i_ox_conversion
      real*8 x,y,z,r,cnx,cny,cnz, cnpar
      real*8 cnrho_ini, cn2_ini, cnpar_ini, cntan2_ini 
      real*8 rho_ini,rho_loc,
     + cnper2p,cnper2m,cnperp_loc,cnperm_loc,cnper2,cn2,
     & cntan2,cnrho2,cnrho,cnperp,cnper_tang, xe, xe_old,xe_new
      real*8 uvtx, uvty, uvtz,  uvnx, uvny, uvnz, cnper_norm
      real*8 dwpw2dx,dwpw2dy,dwpw2dz,dndx,dndy,dndz,gradne,
     + pp,anorm_x,anorm_y,anorm_z , vgroup
      logical was_not_o_cutoff
      character*8 step_dir !Choice of direction: 'gradne' or 'vgroup'

      ! A step size for changing (x,y,z) until X-mode is found:
      hstep=1.d-4*rmax
      
      step_dir= ox_step_dir  !YuP[11-2016] now in namelist /ox/ section
      ! ox_step_dir defines the type of stepping in O-X conversion case,
      !         when looking for X-mode across evanescent layer.
      !         'gradne' - step along grad(ne) direction, or
      !         'vgroup' - along Vgroup velocity direction
      !                  taken just before O-X jump is initiated.
      !         The original (only) option was effectively 'gradne'.
      !         The stepping is done in subroutine find_rho_X_xyz()
      !         until the X mode is found. 
      !         Sometimes such stepping along grad(ne) direction
      !         results in an angled (rather than smooth) ray jump
      !         across the evanescent layer. The second option, 
      !         'vgroup' usually gives a smooth transition,
      !         however, as the ray tends to slow down 
      !         in normal-to-flux-surface direction, 
      !         while the tangential group velocity may remain 
      !         almost unchanged, the ray may "slide"
      !         along the surface, and the O-X transition can be missed.
      if(step_dir.eq.'gradne' .or. step_dir.eq.'vgroup')then
        ! ok
      else
        write(*,*)
     +  'find_rho_X_xyz: ox_step_dir can only be "gradne" or "vgroup" '
        stop
      endif

      pi=4.d0*datan(1.d0)
      cnper2p=0.d0 ! to initialize
      cnper2m=0.d0 ! to initialize

      x= x_ini  !initialization
      y= y_ini  
      z= z_ini 
      cnx= cnx_ini   
      cny= cny_ini   
      cnz= cnz_ini  
      ! N^2 total:
      cn2= cnx**2 + cny**2 + cnz**2
      xe= wpw_2(x,y,z,1) ! to get Xe; it will also get rho value
      rho_ini= rho

      if(step_dir.eq.'gradne')then
        ! Determine the direction of grad(ne).
        ! The stepping/iterations towards X-mode cutoff
        ! will be performed along grad(ne) direction,
        ! which itself can change during stepping, 
        ! so it will be updated at each step
        call dwpw_2(x,y,z,1, dwpw2dx,dwpw2dy,dwpw2dz) 
        dndx= dwpw2dx ! strictly, it is d(wpe/w)^2/dx,
        dndy= dwpw2dy !     but we'll normalize
        dndz= dwpw2dz !     to make a unit vector along grad(ne)
        gradne= dsqrt(dndx**2 + dndy**2 + dndz**2)
        if(gradne.eq.0.d0)then
          write(*,*)'  find_rho_X_xyz ENTER: gradne=0'
          stop
        endif
        pp= 1.0d0/gradne
        ! Unit vector opposite to grad(el.density)  :
        anorm_x= -dndx*pp ! Note: At plasma edge, grad(ne)<0 (inward)
        anorm_y= -dndy*pp
        anorm_z= -dndz*pp !(anorm_x, anorm_y, anorm_z) is pointing outward
      endif

      if(step_dir.eq.'vgroup')then
        !----------------------
        ! Alternative method for stepping_inward (for searching Xmode):
        ! Along group velocity. 
        !In this case, the direction remains unchanged during stepping.
        vgroup= sqrt(vgr_x**2+vgr_y**2+vgr_z**2)
        if(vgroup.eq.0.d0)then
          write(*,*)'  find_rho_X_xyz ENTER: vgroup=0.  xyz=',x,y,z
          stop
        endif
        anorm_x= -vgr_x/vgroup !(anorm_x, anorm_y, anorm_z) is pointing outward
        anorm_y= -vgr_y/vgroup
        anorm_z= -vgr_z/vgroup
        !----------------------
      endif
       
      was_not_o_cutoff= .true.
      id_loc=id
      ioxm_loc=ioxm
       
      write(*,*)'=======> find_rho_X_xyz START:'
      write(*,'(a,5e13.4)')'   x,y,z, Xe,rho=', x,y,z,Xe,rho
      write(*,'(a,3e13.4)')'   anorm_x,anorm_y,anorm_z=',
     +                         anorm_x,anorm_y,anorm_z
      write(*,'(a,3e13.4)')'   vgr_x,vgr_y,vgr_z=',
     +                         vgr_x,vgr_y,vgr_z
      
      iter=0 ! to count iterations -------------------------------------
      xe_old=xe-1.d-3 !starting Xe (should be ~1.0), just before OX-jump
                      !(and reduced a bit, to avoid accuracy issues)
                        
 10   continue ! handle for new iteration with new x,y,z
      iter=iter+1
      write(*,*)'  LOOP-10:iter,was_not_o_cutoff,cnz',
     +    iter,was_not_o_cutoff,cnz
      xe= wpw_2(x,y,z,1) ! to get Xe; it will also get rho value
      xe_new=xe
      if(step_dir.eq.'gradne')then ! step along grad(ne)
        call dwpw_2(x,y,z,1, dwpw2dx,dwpw2dy,dwpw2dz) 
        dndx= dwpw2dx ! strictly, it is d(wpe/w)^2/dx,
        dndy= dwpw2dy !     but we'll normalize
        dndz= dwpw2dz !     to make a unit vector along grad(ne)
        gradne= dsqrt(dndx**2 + dndy**2 + dndz**2)
        if(gradne.eq.0.d0)then
          write(*,*)'  find_rho_X_xyz ENTER: gradne=0'
          stop
        endif
        pp= 1.0d0/gradne
        ! Unit vector opposite to grad(el.density)  :
        anorm_x= -dndx*pp ! Note: At plasma edge, grad(ne)<0 (inward)
        anorm_y= -dndy*pp
        anorm_z= -dndz*pp !(anorm_x, anorm_y, anorm_z) is pointing outward
      endif
      write(*,'(a,i5,3e13.4,f8.5)') '   ioxm,x,y,z,rho=',ioxm,x,y,z,rho
      !pause !!!
 
      if (abs(rho-rho_ini).gt.0.2  .or. xe_new.lt.xe_old) then 
         write(*,*)'  find_rho_X_xyz-10: |rho-rho_ini|>0.2   rho=',rho
         write(*,*)'  Or Xe decreasing  Xe_old,Xe_new=',xe_old,xe_new
         write(*,*)'=================================================='
         write(*,*)'  No roots (Try smaller hstep in find_rho_X_xyz?)'
         write(*,*)'=================================================='
         write(*,*)'  RETURN/Continue from the initial point as O-mode.'
         !stop
         ! For output:
         i_ox_conversion=0 ! means: No conversion took place
         was_not_ox_conversion=.false. !To prevent calling ox_conversion_xyz
         i_call_prep3d_in_output_at_i_ox_conversion_eq_1=0 ! reset
         x_x=x_ini
         y_x=y_ini
         z_x=z_ini
         cnx_x=cnx_ini
         cny_x=cny_ini
         cnz_x=cnz_ini
         bmod=bxyz(x_x,y_x,z_x) ! to find bx,by,bz,bmod -> /one.i/
         cnpar_x= (cnx_x*bx + cny_x*by + cnz_x*bz)/bmod !for the local (x,y,z)
         cnper_x= dsqrt(cnx_x**2 + cny_x**2 + cnz_x**2 - cnpar_x**2)
         id=id_loc ! restore
         ioxm=ioxm_loc
            !pause !!!
         return
      endif
       
      ! Find cnper_tang and cnper_norm from known cnx,cny,cnz.
      ! First, Calculate the components of unit vectors perp to b-field;
      ! one is TANGENTIAL to psi-flux-surface, another is NORMAL to it:
      call unit_vectors_perp_b(x,y,z, 
     +                     uvtx, uvty, uvtz,  uvnx, uvny, uvnz) !->out
      ! Nperp to b-field and tangential to psi-flux-surface:
      cnper_tang= cnx*uvtx + cny*uvty + cnz*uvtz 
      ! Nperp to b-field and normal to psi-flux-surface:
      cnper_norm= cnx*uvnx + cny*uvny + cnz*uvnz 
      ! Npar along b-field:
      cnpar= (cnx*bx + cny*by + cnz*bz) / bmod ! for the local (x,y,z)
      !den=dense_xyz(x,y,z,1) !-> get rho (stored in one.i)
      !rho_loc=rho !from one.i;  
      ! N^2 total:
      cn2=    cnx**2 + cny**2 + cnz**2
      ! N^2 tangential to flux surface:
      cntan2= cnper_tang**2 + cnpar**2  ! or cn2 - cnper_norm**2

      id=2 ! cold plasma approximation

      if (was_not_o_cutoff) then
         ioxm=1 ! search as O-mode   (ioxm is stored in one.i)
      else
         ioxm=-1 ! search as X-mode
      endif
      
      ! Now find new Nperp^2
      call npernpar_xyz(x,y,z,cnpar,cnper2p,cnper2m) !-> cnper2p,cnper2m

      write(*,'(a,3e13.4)')
     +  '  find_rho_X npernpar: Xe,cnper2p,cnper2m',Xe,cnper2p,cnper2m

      ! skip this part when was_not_o_cutoff=True 
      ! (meaning: still an O-mode, before cutoff)
      if (.not.was_not_o_cutoff) then 
        ! search for X-mode when was_not_o_cutoff=False
        if ((cnper2p.le.0.d0).and.(cnper2m.le.0.d0)) then ! both roots<0
           write(*,*)'  find_rho_X_xyz: Nper2p<0 and Nper2m<0. New iter'
           x= x - hstep*anorm_x ! step inward, along -grad(ne) direction
           y= y - hstep*anorm_y
           z= z - hstep*anorm_z
           !xe_old= xe_new ! update for the next iteration; 
           !Do not "update" xe_old: 
           !Rather, compare Xe_new with Xe before OX-jump (xe_old, not updated).
           !Why: because in one step, Xe_new and Xe at previous step
           !can be so close that within accuracy, Xe_new may look 
           !like decreasing.
           goto 10
        else  !-> at least one root is positive.
           ! Select the smaller (but positive) root:
           cnper2=min(cnper2p,cnper2m)
           if (cnper2.le.0.d0) then
             cnper2=max(cnper2p,cnper2m)
           endif 
           cn2= cnper2+cnpar**2
           if(cn2 .lt. cntan2)then
             write(*,*)'  find_rho_X_xyz: cn2<cntan2. New iter.'
             x= x - hstep*anorm_x ! step inward
             y= y - hstep*anorm_y
             z= z - hstep*anorm_z
             !xe_old= xe_new ! for the next iteration
             goto 10
           endif

           ! Adjust Nrho (keeping Npar and Nper_tang unchanged)
           cn2= cnper2+cnpar**2
           cnrho2= cn2-cntan2
           cnrho=  dsqrt(cnrho2)
           cnper_norm= cnrho
           call cnxyz(x,y,z, cnpar, cnper_tang, cnper_norm,
     +                 cnx,cny,cnz) !->out
           write(*,*)'  find_rho_X: found X-mode. cnper2=',cnper2
           write(*,*)'  find_rho_X: found X-mode. cnpar =',cnpar
           write(*,*)'  find_rho_X: found X-mode. cn2   =',cn2
           write(*,*)'  find_rho_X: found X-mode. cnx=',cnx
           write(*,*)'  find_rho_X: found X-mode. cny=',cny
           write(*,*)'  find_rho_X: found X-mode. cnz=',cnz
           iraystop=0 
           ! Almost done...
           goto 30
        endif
      endif
      
      ! solve the dispersion relation N=N(n_par) for cold plasma:
      write(*,*)'  find_rho_X_xyz before cninit12_xyz: cnpar,cnz=',
     + cnpar,cnz
      call cninit12_xyz(x,y,z, cnpar, cnper_tang,
     &                  cnx,cny,cnz,iraystop) !-> cnx,cny,cnz
      write(*,*)'  find_rho_X_xyz after  cninit12_xyz: cnpar,cnz=',
     + cnpar,cnz

 30   continue
      write(*,*)'  find_rho_X_xyz-30: ioxm,iraystop===',ioxm,iraystop
      !if(iraystop.eq.1) pause !!!

      id=id_loc
      ioxm=ioxm_loc

      if (was_not_o_cutoff) then ! true, meaning - looking for O-mode
      
        if(iraystop.eq.0) then
c---------ioxm=1 O-mode exists in this point
          write(*,*)'  find_rho_X_xyz: Omode still exist here. New iter'
          x= x - hstep*anorm_x ! step inward
          y= y - hstep*anorm_y
          z= z - hstep*anorm_z
          !xe_old= xe_new ! for the next iteration
          go to 10                   
        else ! iraystop=1: Reached a point where O-mode vanishes (Nper^2<0)
          was_not_o_cutoff=.false.   ! Switch to a search of X-mode
          ioxm=-1 ! search as X-mode
          write(*,*)'  was_not_o_cutoff=F. No O-mode; searching X-mode'
          !xe_old=xe
 20       continue ! stepping inward until Xe=1.0 (or larger)
          xe= wpw_2(x,y,z,1)
          xe_new=xe
          write(*,*)'   LOOP-20: was_not_o_cutoff=', was_not_o_cutoff,
     +              '   Xe=', xe
          if(step_dir.eq.'gradne')then
            call dwpw_2(x,y,z,1, dwpw2dx,dwpw2dy,dwpw2dz) 
            dndx= dwpw2dx ! strictly, it is d(wpe/w)^2/dx,
            dndy= dwpw2dy !     but we'll normalize
            dndz= dwpw2dz !     to make a unit vector along grad(ne)
            gradne= dsqrt(dndx**2 + dndy**2 + dndz**2)
            if(gradne.eq.0.d0)then
              write(*,*)'  find_rho_X_xyz ENTER: gradne=0'
              stop
            endif
            pp= 1.0d0/gradne
            ! Unit vector opposite to grad(el.density)  :
            anorm_x= -dndx*pp ! Note: At plasma edge, grad(ne)<0(inward)
            anorm_y= -dndy*pp
            anorm_z= -dndz*pp !(anorm_x, anorm_y, anorm_z) is outward
          endif
          if (abs(rho-rho_ini).gt.0.2  .or. xe_new.lt.xe_old) then 
          write(*,*)'  find_rho_X_xyz-20: |rho-rho_ini|>0.2   rho=',rho
          write(*,*)'  Or Xe decreasing  xe_old,xe_new=',xe_old,xe_new
          write(*,*)'=================================================='
          write(*,*)'  No roots (Try smaller hstep in find_rho_X_xyz?)'
          write(*,*)'=================================================='
          write(*,*)'  RETURN/Continue from the initial point as O-mode'
             !stop
             ! For output:
             i_ox_conversion=0 ! means: No conversion took place
             was_not_ox_conversion=.false. !To prevent calling ox_conversion_xyz
             i_call_prep3d_in_output_at_i_ox_conversion_eq_1=0 ! reset
             x_x=x_ini
             y_x=y_ini
             z_x=z_ini
             cnx_x=cnx_ini
             cny_x=cny_ini
             cnz_x=cnz_ini
             bmod=bxyz(x_x,y_x,z_x) ! to find bx,by,bz,bmod -> /one.i/
             cnpar_x= (cnx_x*bx+cny_x*by+cnz_x*bz)/bmod !for the local (x,y,z)
             cnper_x= dsqrt(cnx_x**2 + cny_x**2 + cnz_x**2 - cnpar_x**2)
             id=id_loc ! restore
             ioxm=ioxm_loc
                !pause !!!
             return
          endif
          if (xe.le. 1.d0) then
             write(*,'(a,e13.4, f8.5)')
     +         '    find_rho_X_xyz-20: reducing (x,y,z). Xe,rho=',xe,rho
             x= x - hstep*anorm_x ! step inward until Xe>1
             y= y - hstep*anorm_y
             z= z - hstep*anorm_z
             !xe_old= xe_new ! for the next iteration
             goto 20  ! stepping inward until Xe=1.0 (or larger)
          endif
          write(*,'(a,f8.3)')'    find_rho_X_xyz:goto10 SearchX; Xe=',Xe
          !pause !!!
          !xe_old= xe_new ! for the next iteration
          call cninit12_xyz(x,y,z, cnpar, cnper_tang,
     &                  cnx,cny,cnz,iraystop) !-> cnx,cny,cnz
c      write(*,*)'  find_rho_X_xyz-cninit12_xyz/X-mode/Xe>1: cnpar,cnz=',
c     + cnpar,cnz
          goto 10 ! Xe>1.0 is achieved. Search for X-mode 
        endif
        
      else ! was_not_o_cutoff=false (was Looking for X-mode)
          
        if(iraystop.eq.1) then
          ! Never happens?
          ! ioxm=-1 X-mode does not exist at this point
          write(*,*)'  find_rho_X_xyz: X-mode not here yet. Step inward'
          x= x - hstep*anorm_x ! step inward
          y= y - hstep*anorm_y
          z= z - hstep*anorm_z
          !xe_old= xe_new ! for the next iteration
          !pause !!!
          goto 10   
        else ! iraystop=0   Final result for X-mode: 
          x_x=x
          y_x=y
          z_x=z
          cnx_x=cnx
          cny_x=cny
          cnz_x=cnz
          cnpar_x=cnpar
          cnper_x=dsqrt(cnx**2 + cny**2 + cnz**2 - cnpar**2)
          id=id_loc
          ioxm=ioxm_loc
          ! Finished
        endif
        
      endif

      xe=wpw_2(x,y,z,1)
      write(*,*)'=======> find_rho_X_xyz END/RETURN  Xe,rho=',xe,rho
        !pause !!!
      return
      end


c======================================================================
c======================================================================

      subroutine ox_conversion_xyz(x_in,y_in,z_in, nx_in,ny_in,nz_in,
     &eps_xe,  ! The parameter which sets the vicinitity
               ! of the O-mode cutoff surface
               ! If xe > (1-eps_xe) then this subroutine  
               ! makes the ray jump in small radius direction
               ! and finds X mode.
c-------the point is near the OX mode conversion area
     &x_x,y_x,z_x, nx_x,ny_x,nz_x, i_ox_conversion)
c     Creates the jump of the ray point through the OX mode conversion
c     area where X_e=1 (V_perp=0)

      implicit none               

c-----input
      real*8 x_in,y_in,z_in,nx_in,ny_in,nz_in !the coordinates before
                                               !O mode cutoff
      real*8 eps_xe

c-----output
      real*8 x_x,y_x,z_x,nx_x,ny_x,nz_x !the coordinates of
                                         !X mode cutoff point
      integer i_ox_conversion !=0 was not OX conversion
                              !=1 was the jump in the radial direction
c-----locals
      real*8 xe_0,xe_in,bmod,
     &cnpar,
     &cnpar_x,cnper_x,m_in_loc,
     &transm_ox_loc !for test

      integer iraystop

c-----externals
      real*8 wpw_2,bxyz
 
      i_ox_conversion=0
      xe_0=1.d0

      xe_in= wpw_2(x_in,y_in,z_in,1)

      call find_maximal_OX_transmission_xyz(xe_0,
     & x_in,y_in,z_in, nx_in,ny_in,nz_in,
     & i_ox_conversion) ! redefines x_in,y_in,z_in, nx_in,ny_in,nz_in
                        ! to the point along ray 
                        ! where transm.coeff is highest
      
      ! YuP[Nov-2014] Why needed? The transm.coef is calculated later,
      ! in outpt_xyz.  Commenting this call:
cyup      call transmit_coef_ox_xyz(x_in,y_in,z_in, nx_in,ny_in,nz_in,
cyup     &                           transm_ox_loc)

      if(i_ox_conversion.eq.1) then
        !find the point where Xmode with nx_x,ny_x,nz_x has the cutoff.
        call find_rho_X_xyz(x_in, y_in, z_in, nx_in, ny_in, nz_in,
     &    x_x, y_x, z_x, nx_x,ny_x,nz_x, cnpar_x,cnper_x,
     +    i_ox_conversion) !-> last 2 lines: out
      endif
      return
      end
      

c======================================================================
c======================================================================

c        ********************* reflect_xyz ********************
c        *                        -                           *
c        * This subroutine reflects                           *
c        * the ray from plasma boundary	                      *
c        ******************************************************
c
c------------------------------------------------------------------
c								   
c      Input:					   
c      								   
c      x,y,z - cartesian coords. at point of reflection
c      cnx,cny,cnz						  
c------------------------------------------------------------------!
c      Output:					   
c	       	cnxref,cnyref,cnzref= N_x,N_y,N_z after reflection
c----------------------------------------------------------------------
      subroutine reflect_xyz(x,y,z,cnx,cny,cnz, cnxref,cnyref,cnzref)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
      double precision
     1 ias1r,ias2r,ias2r_Sm, bxyz

      nr4=nx+4
      nz4=ny+4
      ncx=nr4
      ncy=nz4
      rabs=sqrt(x*x+y*y+z*z) ! r to the given point
      
      if(model_rho_dens.eq.2 .or. model_rho_dens.eq.3 .or. 
     +   model_rho_dens.eq.4) then !reflect at wall
        epsbnd=2.d-3 ![m]! Should be exactly same as in bound_xyz()
        r=sqrt(x*x+y*y)
        if(r.ge. wall_rmax-epsbnd) then
          ! Reflection from cylinder wall.
          ! Unit vector perp to wall  :
          anorm_x= x/r !  (outward)
          anorm_y= y/r
          anorm_z= 0.d0
          ! Magnitude of N-component along the unit vector :
          cnpsi= (anorm_x*cnx + anorm_y*cny + anorm_z*cnz)
          ! Components of Npsi (vector):              
          cnpsi_x= anorm_x*cnpsi
          cnpsi_y= anorm_y*cnpsi
          cnpsi_z= anorm_z*cnpsi
          ! Reflect N by subtracting 2*Npsi as vector
          cnxref= cnx -2.d0*cnpsi_x
          cnyref= cny -2.d0*cnpsi_y
          cnzref= cnz -2.d0*cnpsi_z
          return
        endif
      endif
      
cyup      bmod= bxyz(x,y,z) !-> get components of grad(psi) in one.i
cyup      gradpsi= dsqrt(dpdxd*dpdxd+dpdyd*dpdyd+dpdzd*dpdzd)
c----> YuP: New method: use grad(el.density) instead of grad(psi)
c----> They could be the same (in tokamaks).
      call dwpw_2(x,y,z,1, dwpw2dx,dwpw2dy,dwpw2dz) ! derivatives
      wpw2= wpw_2(x,y,z,1) ! = (wpe/w)^2
      
      dndx= dwpw2dx ! strictly, it is d(wpe/w)^2/dx,
      dndy= dwpw2dy !     but we'll normalize
      dndz= dwpw2dz !     to make a unit vector along grad(ne)
      gradne= dsqrt(dndx**2 + dndy**2 + dndz**2)
      dndr=gradne
      grpsi=  dsqrt(dpdxd**2 +dpdyd**2 +dpdzd**2) ! |grad(psi)|
     
      !YuP[12-2016] added more logical branches for reflection procedure:
      !if grad(ne)~0, use grad(psi); If both ~0, use local r vector.
      if( dndr*rabs .gt. wpw2*1.e-8 )then
        pp= 1.0d0/dndr
        ! Unit vector opposite to grad(el.density)  :
        anorm_x= -dndx*pp ! Note: grad(ne) <0 (INward)
        anorm_y= -dndy*pp
        anorm_z= -dndz*pp
      elseif( grpsi*rabs .gt. psid*1.e-8 )then
        !grad(n)~0; Try using grad(psi):
        pp= 1.0d0/grpsi
        ! Unit vector along grad(psi)  :
        anorm_x= dpdxd*pp ! Note: grad(psi) >0 (OUTward)
        anorm_y= dpdyd*pp
        anorm_z= dpdzd*pp
      else ! both grad(n) and grad(psi) are almost zero;
        ! then - use (x,y,z) direction
        pp= 1.0d0/rabs
        ! Unit vector along r vector  :
        anorm_x= x*pp 
        anorm_y= y*pp
        anorm_z= z*pp
      endif

      ! Magnitude of N-component along grad() direction:
      cnpsi= (anorm_x*cnx + anorm_y*cny + anorm_z*cnz)
      ! Components of Npsi (vector):              
      cnpsi_x= anorm_x*cnpsi
      cnpsi_y= anorm_y*cnpsi
      cnpsi_z= anorm_z*cnpsi
      ! Reflect N by subtracting 2*Npsi as vector
      cnxref= cnx -2.d0*cnpsi_x
      cnyref= cny -2.d0*cnpsi_y
      cnzref= cnz -2.d0*cnpsi_z
c      write(*,*)'reflect_xyz: cnx,cny,cnz',cnx,cny,cnz
c      write(*,*)'cnpsi_x,cnpsi_y,cnpsi_z',cnpsi_x,cnpsi_y,cnpsi_z
c      write(*,*)'reflect_xyz: cnxref,cnyref,cnzref',cnxref,cnyref,cnzref
      return
      end



c======================================================================
c======================================================================


      subroutine wall_limiter_N_reflection_xyz(x,y,z, cnx,cny,cnz, 
     & n_rho, i_wall_ar_p, i_wall_ar_m,
     + cnxrefl,cnyrefl,cnzrefl) !-> out
c---------------------------------------------------------------------------
c     calculate wall refractive index after reflection:
c     cnxrefl,cnyrefl,cnzrefl
c--------------------------------------------------------------------
c     the wall element M_p(z_p,r_p),M_m(z_m,r_m) near the reflection 
c     point is a straight line
c
c     The unit vector parallel to the wall element lw(lw_z.lw_r):
c
c     lw_z=(z_p-z_m)/sqrt((z_p-z_m)**2+(r_p-r_m)**2)
c     lw_r=(r_p-r_m)/sqrt((z_p-z_m)**2+(r_p-r_m)**2)
c     
c     The unit vector perpendicular to the wall element Nw(nw_z,nw_r)
c     is the solution of the following equation: 
c     lw_z*nm_z+lw_r*nm_r=0
c
c     If (lw_z.ne.0) then nm_z=-(lw_r/lw_z)*nm_r
c     else (lw_z.ne.0) then nm_r=-(lw_z/lw_r)*nm_z
c
c     Unit Nw: nm_z/sqrt(mn_z**2+mnr**2),nm_r/sqrt(mn_z**2+mnr**2)
c
c     The refractive vector perpendicular to the wall element :
c     N_perp^=(cnz_perp,cnr_perp)
c
c     (N^.Nw^)=cnz*nm_z+cnr*nm_r
c
c     cnr_perp= (N^.Nw^)*nm_r,  cnr_perp= (N^.Nw^)*nm_r
c
c     The reflected refractive index N_refl^=N^-2*N_perp^

      implicit none
      include 'param.i'
      include 'one_nml.i'
      include 'fourb.i'
c-----input 
      integer 
     &n_rho,                      !index of point at which wall small radius
                                  !is close to ray small radius 
                                  !at the reflection point.
     &i_wall_ar_m(n_rho_wall_a),  !indices of wall points on a straight 
     &i_wall_ar_p(n_rho_wall_a)   !line in poloidal direction.
      real*8 x,y,z       ! cartesian coords. at refl.point
      real*8 cnx,cny,cnz !refractive index before reflection
c-----output
      real*8
     &cnxrefl,         ! refractive index after reflection            
     &cnyrefl,         ! at wall element 
     &cnzrefl          ! 
c-----locals
      real*8
     &z_p,r_p,z_m,r_m, !coordinates of points M_p,M_m
     &lw_z,lw_r,       !unit vector along wall element (M_m,M_p)
     &nw_r,nw_x,nw_y,nw_z, !unit vector normal to wall element (M_m,M_p)
     &cnx_perp_w,      !The refractive vector perpendicular 
     &cny_perp_w,      !to the wall element 
     &cnz_perp_w,      !
     &p, r
c-------------------------------------------------------------------
c     the wall element M_p(z_p,r_p),M_m(z_m,r_m) near the reflection 
c     point 
c------------------------------------------------------------------
      z_p=z_wall(i_wall_ar_p(n_rho)) ! For now, assume the wall is 
      r_p=r_wall(i_wall_ar_p(n_rho)) ! axially symmetric => no phi-dep.
      z_m=z_wall(i_wall_ar_m(n_rho))
      r_m=r_wall(i_wall_ar_m(n_rho))
c------------------------------------------------------------------
c     The unit vector parallel to the wall element lw(lw_z.lw_r):
c     lw_z=(z_p-z_m)/sqrt((z_p-z_m)**2+(r_p-r_m)**2)
c     lw_r=(r_p-r_m)/sqrt((z_p-z_m)**2+(r_p-r_m)**2)
c------------------------------------------------------------------
      lw_z= z_p-z_m
      lw_r= r_p-r_m
      p= 1.d0/dsqrt(lw_z**2+lw_r**2)
      lw_z= lw_z*p
      lw_r= lw_r*p
      write(*,*)'lw_z,lw_r',lw_z,lw_r
c-----------------------------------------------------------------------
c     The unit vector perpendicular to the wall element Nw(nw_z,nw_r)
c     is the solution of the following equation: 
c     lw_z*nw_z+lw_r*nw_r=0
c-----------------------------------------------------------------------
      if (lw_z.ne.0) then
         nw_r=1.d0
         nw_z=-(lw_r/lw_z)*nw_r
         p=1.d0/dsqrt(nw_z**2+nw_r**2)
         nw_z=  nw_z*p
         nw_r=  nw_r*p
      else
c--------lw_r.ne.0) 
         nw_z=1.d0
         nw_r=-(lw_z/lw_r)*nw_z
         p=1.d0/dsqrt(nw_z**2+nw_r**2)
         nw_z=  nw_z*p
         nw_r=  nw_r*p
      endif
      r= dsqrt(x*x+y*y)
      nw_x= nw_r*x/r ! For now, assume the wall is axially symmetric
      nw_y= nw_r*y/r ! so that the nw vector has no phi-component
      write(*,*)'nw_z,nw_r',nw_z,nw_r
      write(*,*)'lw_z*nw_z+lw_r*nc',lw_z*nw_z+lw_r*nw_r
c-----------------------------------------------------------------------
c     The refractive vector perpendicular to the wall element :
c     N_perp^= (cnx_perp_w, cny_perp_w, cnz_perp_w)
c
c     (N^.Nw^)= cnx*nw_x + cny*nw_y+ cnz*nw_z
c
c     cnx_perp_w= (N^.Nw^)*nw_x,  etc.
c-----------------------------------------------------------------------
      p= cnx*nw_x + cny*nw_y+ cnz*nw_z !(N^.Nw^)
      write(*,*)'(N^.Nw^)',p
      cnx_perp_w= nw_x*p
      cny_perp_w= nw_y*p
      cnz_perp_w= nw_z*p
c------------------------------------------------------------------------
c     The refactive index after reflection N_refl^= N^-2*N_perp^
c-----------------------------------------------------------------------
      cnxrefl= cnx -2.d0*cnx_perp_w
      cnyrefl= cny -2.d0*cny_perp_w
      cnzrefl= cnz -2.d0*cnz_perp_w
      write(*,*)'cnz,cnz_perp_w,cnzrefl',cnz,cnz_perp_w,cnzrefl
      return
      end

c======================================================================
c======================================================================

      subroutine callc_eps_ham_xyz(u,ham)
c-------------------------------------
c     calculates hamiltonian
c     to check rk accuracy
c-----------------------------
      implicit none
      include 'param.i'
      include 'one.i'
c-----input
      real*8 u(6)    !ray cartesian coordinates and N-components
c-----output
      real*8 ham
c-----externals
      real*8 bxyz,gamma1_xyz,hamilt_xyz
c-----local
      real*8 cnt2      

      bmod= bxyz(u(1),u(2),u(3)) !-> get b and derivs of b
      gam=  gamma1_xyz(u(1),u(2),u(3),u(4),u(5),u(6))
      cnt2= u(4)**2 + u(5)**2 + u(6)**2
      ham=  hamilt_xyz(u(1),u(2),u(3),cnt2)

      return
      end

c======================================================================
c======================================================================

      subroutine rho_wall_at_theta_xyz(x_q,y_q,z_q,
     + rho_wall_ar, n_rho_wall, !-> out
     &i_wall_ar_m, i_wall_ar_p) !-> out
c-----calculates wall small radius at poloidal angle theta [radian]
c     of the given ray point (x_q,y_q,z_q)
c     It is assumed that between two neighboring wall points
c     the wall small radius is a linear function on poloidal angle.
c     inside this subroutine 0=<theta=thetapol<2*pi
      implicit none
      include 'param.i'
      include 'one_nml.i'
      include 'three.i'
      include 'fourb.i'
c-----input
      real*8
     &x_q,y_q,z_q  !cartesian coordinates of the ray point
c-----output
      integer  n_rho_wall,        !number of small radii at the given poloidal angle
     &i_wall_ar_m(n_rho_wall_a),  !numbers of wall points around the straight 
     &i_wall_ar_p(n_rho_wall_a)   !line from the magnetic directed at the 
                                  !poloidal angle theta
      real*8 rho_wall_ar(n_rho_wall_a) ! small radii at given poloidal angle
    
c-----locals
      integer i,j
      real*8     
     &theta,       !poloidal angle [radians] of the ray point(z_q,r_q)
     & pi, r_p,z_p, r_m,z_m,
     &a,b,a_q,b_q, r_q,
     &r_b, x_b,y_b,z_b,  l_z_b_r_b,  r_k,  z_k,
     &r_l, z_l         

c      write(*,*)'In rho_wall_at_theta_xyz'

      pi=4*datan(1.d0)

      r_q= dsqrt(x_q**2 + y_q**2)
      r_l= r_q-rma 
      z_l= z_q-zma
      call theta_rz(r_l,z_l, theta)  ! get pol.angle 0=< theta <2*pi
      
      n_rho_wall=0
c------------------------------------------------------------------------
c     determine 
c     numbers of wall points:
c     i_wall_ar_m(1:n_rho__wall),i_wall_ar_p(1:n_rho__wall),
c
c     thetapol_wall(i-1) =< theta < thetapol_wall(i)
c     or
c     thetapol_wall(i) < theta =< thetapol_wall(i-1)
c     
c-------------------------------------------------------------------------
      do i=1,n_wall
        if (i.eq.1) then
           if((theta.le.thetapol_wall(1)).and.(theta.ge.0.d0)) then
              n_rho_wall=n_rho_wall+1
              i_wall_ar_m(n_rho_wall)=n_wall-1
              i_wall_ar_p(n_rho_wall)=1
           endif
           goto 10
        endif

        if (i.eq.n_wall) then
          if((theta.gt.thetapol_wall(n_wall-1)).and.
     &      (theta.le.2.d0*pi)) then
          n_rho_wall=n_rho_wall+1
          i_wall_ar_m(n_rho_wall)=n_wall-1
          i_wall_ar_p(n_rho_wall)=n_wall
          endif
          goto 10
        endif

        if(((theta.gt.thetapol_wall(i-1)).and.
     &      (theta.le.thetapol_wall(i))).or.
     &     ((theta.gt.thetapol_wall(i)).and.
     &      (theta.le.thetapol_wall(i-1)))) then
          n_rho_wall=n_rho_wall+1
          i_wall_ar_m(n_rho_wall)=i-1
          i_wall_ar_p(n_rho_wall)=i
          write(*,*)'n_rho_wall',n_rho_wall
        endif

 10   continue
      enddo ! i=1,n_wall
c--------------------------------------------------------------------------
c
c     equation for the wall small radius versus poloidal angle theta
c
c     rho=rho_wall_m+(theta-thetapol_wall_m)/(thetapol_wall_p-thetapol_wall_m)    
      
      do j=1, n_rho_wall

         r_p=r_wall(i_wall_ar_p(j))
         z_p=z_wall(i_wall_ar_p(j))
         r_m=r_wall(i_wall_ar_m(j))
         z_m=z_wall(i_wall_ar_m(j))

         if (r_p.ne.r_m) then
            b=(z_p-z_m)/(r_p-r_m)
            a=z_m-b*r_m

            if (r_q.ne.rma) then
              b_q=(z_q-zma)/(r_q-rma)
              a_q=z_q-b_q*r_q
c-------------intersection point K
              r_k=-(a-a_q)/(b-b_q)
              z_k=a+b*r_k
            else
c-------------r_q.eq.rma
              r_k=rma
              z_k=a+b*r_k
            endif 

         else
c-----------r_m=r_p
            if (r_q.ne.rma)then 
               b_q=(z_q-zma)/(r_q-rma)
               a_q=z_q-b_q*r_q
               r_k=r_m
               z_k=a_q+b_q*r_k
            else
c-------------r_q.eq.rma
              r_k=r_m
              write(*,*)'r_m.eq.r_p, r_1.eq_r_m'
              stop ' rho_wall_at_theta_2'
            endif

         endif
  
c-----calculate coordinates r_b,z_b at the LCFS surface psi=psilim
c     and the poloidal angle theta at the K point
         r_l=r_k-rma 
         z_l=z_k-zma
         call theta_rz(r_l,z_l,theta)
         call zr_psith(psilim,theta,z_b,r_b)   
         l_z_b_r_b=dsqrt((r_b-rma)**2+(z_b-zma)**2)  ! the distance between
                                                     ! the magnetic axis and
                                                     ! LCFS at theta
         rho_wall_ar(j)=dsqrt((r_k-rma)**2+(z_k-zma)**2)/l_z_b_r_b
      enddo
      ! needs work? assumed axial symmetry for walls
      return
      end

c======================================================================
c======================================================================

       subroutine wall_limiter_reflection_point_xyz(x_n,y_n,z_n,
     &cnx_n,cny_n,cnz_n, i_reflection,
     &cnx_refl,cny_refl,cnz_refl, x_refl,y_refl,z_refl)
c---------------------------------------------------------------------------
c     calculate wall or limiter reflection point coordinates:
c     x_refl,y_refl,z_refl,
c     and the refractive index after reflection
c     cnx_refl,cny_refl,cnz_refl
c     output: 
c     index i_reflection=1 if the given point x_n,y_n,z_n
c                          and the previous ray 'old' point 
c                          are at the different sides of the wall.
c                       =0 x_n,y_n,z_n point does not give the reflection
c                          
c     x_refl,y_refl,z_refl  - reflection point coordinates:
c     cnx_old,cny_old,cnz_old - refractive index at
c                               the previous ray point
c 
c--------------------------------------------------------------------
c     the wall element M_p(z_p,r_p),M_m(z_m,r_m) near the reflection 
c     point is a straight line
c
c     The unit vector parallel to the wall element lw(lw_z.lw_r):
c
c     lw_z=(z_p-z_m)/sqrt((z_p-z_m)**2+(r_p-r_m)**2)
c     lw_r=(r_p-r_m)/sqrt((z_p-z_m)**2+(r_p-r_m)**2)
c     
c     The unit vector perpendicular to the wall element Nw(nw_z,nw_r)
c     is the solution of the following equation: 
c     lw_z*nm_z+lw_r*nm_r=0
c
c     If (lw_z.ne.0) then nm_z=-(lw_r/lw_z)*nm_r
c     else (lw_z.ne.0) then nm_r=-(lw_z/lw_r)*nm_z
c
c     Unit Nw: nm_z/sqrt(mn_z**2+mnr**2),nm_r/sqrt(mn_z**2+mnr**2)
c
c     The refractive vector perpendicular to the wall element :
c     N_perp^=(cnz_perp,cnr_perp)
c
c     (N^.Nw^)=cnz*nm_z+cnr*nm_r
c
c     cnr_perp= (N^.Nw^)*nm_r,  cnr_perp= (N^.Nw^)*nm_r
c
c     The reflected refactive index N_refl^=N^-2*N_perp^
c-------------------------------------------------------------------
      implicit none
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'fourb.i'

c-----input
      real*8 x_n,y_n,z_n, !ray coordinates
     &cnx_n,cny_n,cnz_n   !ray refractive index 
c-----output
      integer i_reflection !=1 at reflection point
                           !=0 no reflection
      real*8 x_refl,y_refl,z_refl, !reflection point coordinates
     &cnx_refl,cny_refl,cnz_refl     !refractive index in the ray point
                                     !after reflection
c-----locals
      real*8
     &x_old,y_old,z_old,   !ray coordinates at previous point
     &thetapol_n,          !poloidal angle at given (z_n,r_n) point
     &x_l,y_l,z_b,r_b,l_zr,l_z_b_r_b,rho_n,rho_old,
     &rho_wall_ar(n_rho_wall_a), !wall small radii at given thetapol_n
     &dif_rho_min,               !min(abs(rho_n-rho_wall_arr))   
     &dif_rho,                   !rho_n-rho_wall_ar at minimal
                                 ! abs(rho_n-rho_wall_arr))   
     &dif_rho_old, r_n, r_l, z_l,
     &cnx_old,cny_old,cnz_old,   !refractive index at the previous point
cfor test
     &u_l(6),hamilt                                   
cend for test 
      integer i,
     &n_rho_wall,                 !number of wall small radii at 
                                  !given thetapol_n
     &n_rho,                      !number of point at which wall small radius
                                  !is close to ray small radius 
                                  !at given polidal angle in reflection point
     &i_wall_ar_m(n_rho_wall_a),  !numbers of wall points around the straight 
     &i_wall_ar_p(n_rho_wall_a)   !line from the magnetic directed at the 
                                  !poloidal angle theta
     
      logical first

      data first/.true./
      save x_old,y_old,z_old, rho_old,dif_rho_old,
     &cnx_old,cny_old,cnz_old 

      i_reflection=0
      
      if(rho.le.rhowall) goto 10 !nothing to do
cyup      write(*,*)'wall_limiter_reflection: rho,rhowall=',rho,rhowall

c      write(*,*)'n_wall',n_wall
      if (n_wall.eq.0) then
c--------There is no wall 
         goto 10
      endif

c      write(*,*)'wall_limiter_reflection_point rho',rho
      
      r_n= dsqrt(x_n*x_n + y_n*y_n)
      r_l=r_n-rma 
      z_l=z_n-zma
c
      call theta_rz(r_l,z_l,thetapol_n) !poloidal angle at given (z_n,r_n) point

c      write(*,*)'z_n,r_n,thetapol_n', z_n,r_n,thetapol_n
c
c-----calculate coordinates r_b,z_b at the LCFS surface psi=psilim
c     and the poloidal angle thetapol_wall
      call zr_psith(psilim,thetapol_n,z_b,r_b)

      l_zr=dsqrt(r_l**2+z_l**2)                   !distance between the magnetic axis 
                                                  !and the ray point
      l_z_b_r_b=dsqrt((r_b-rma)**2+(z_b-zma)**2)  !distance between the magnetic axis and
                                                  !LCFS at thetapol_wall_n 
      rho_n = l_zr / l_z_b_r_b                    !small radius at given z_n,r_n

c      write(*,*)'rho_n',rho_n
c-----------------------------------------------------------------------------
c     calculates wall small radius at the given point
c     It is assumed that between two wall neighboring points
c     the wall small radius is a linear function on poloidal angle
c-----------------------------------------------------------------------------
      call rho_wall_at_theta_xyz(x_n,y_n,z_n,
     +rho_wall_ar, n_rho_wall,  !-> out
     &i_wall_ar_m, i_wall_ar_p) !-> out

      if (n_rho_wall.eq.0) then
c--------There is no ray-wall intersection 
         goto 10
      endif
c      write(*,*)'after rho_wall_at_theta_2 n_rho',n_rho_wall
c      write(*,*)'rho_wall_ar(i)',(rho_wall_ar(i),i=1,n_rho_wall)
c-------------------------------------------------------------------
c     choose small wall radius closer to rho_n
c--------------------------------------------------------------------
c      write(*,*)'rho_n',rho_n
      n_rho=0
      dif_rho_min=1.d9
      do i=1,n_rho_wall
        if ( dif_rho_min .gt. dabs(rho_n-rho_wall_ar(i))) then
           dif_rho_min=dabs(rho_n-rho_wall_ar(i))
           dif_rho=rho_n-rho_wall_ar(i)
           n_rho=i
        endif 
      enddo
c      write(*,*)'n_rho',n_rho
c      write(*,*)'first',first       
      if (first) then
          x_old= x_n
          y_old= y_n
          z_old= z_n
          cnx_old= cnx_n
          cny_old= cny_n 
          cnz_old= cnz_n 
          dif_rho_old=-1.d0 
          first=.false.
      endif
c      write(*,*)'dif_rho_old,dif_rho',dif_rho_old,dif_rho
c      write(*,*)'cnz_n,cnr_n,cm_n',cnz_n,cnr_n,cm_n
c      write(*,*)'cnz_old,cnr_old,cm_old',
c     &           cnz_old,cnr_old,cm_old

      if (dif_rho_old*dif_rho.lt.0.d0) then
c--------ray-wall intersection point is between old and new points
c        The simple solution is to use the old point as a reflection point. 
         i_reflection=1
         x_refl= x_old
         y_refl= y_old
         z_refl= z_old  
         first=.true.
c--------check hamiltonian value before reflection
c         u_l(1)=x_old
c         u_l(2)=y_old 
c         u_l(3)=z_old
c         u_l(4)=cnx_old
c         u_l(5)=cny_old 
c         u_l(6)=cnz_old
c         call callc_eps_ham_xyz(u_l,hamilt)
c         write(*,*)'before reflection: hamilt=',hamilt
c         write(*,*)'before reflection: u_l=',u_l
c--------calculate refractive index after reflection
         call wall_limiter_N_reflection_xyz(x_old,y_old,z_old,
     +       cnx_old,cny_old,cnz_old, n_rho, i_wall_ar_p, i_wall_ar_m,
     &       cnx_refl,cny_refl,cnz_refl) !-> out

         cnx_old= cnx_refl
         cny_old= cny_refl
         cnz_old= cnz_refl
c---------check hamiltonian value after reflection
c          u_l(1)= x_refl
c          u_l(2)= y_refl 
c          u_l(3)= z_refl
c          u_l(4)= cnx_refl
c          u_l(5)= cny_refl
c          u_l(6)= cnz_refl
c          call callc_eps_ham_xyz(u_l,hamilt)
c          write(*,*)'after reflection: hamilt=',hamilt
c          write(*,*)'after reflection: u_l=',u_l
          goto 10
      endif

c-----set old point equal to given point
      x_old= x_n
      y_old= y_n     
      z_old= z_n
      rho_old=rho_n
      dif_rho_old=dif_rho 
      cnx_old= cnx_n
      cny_old= cny_n
      cnz_old= cnz_n

 10   continue
      return
      end

c======================================================================
c======================================================================

c        *********************BOUND_xyz ***********************
c        *                        -                           *
c        * This function gives the ray position in            *
c        * the plasma (If the point of the ray is inside the  *
c        * grid limits, then BOUND=-1 , else BOUND=1),	      *
c        * and it reflects the ray from plasma boundary.      *
c        * It calculates the number of the reflections =irefl *
c        * It gives the command to stop ray calculation	      *
c        * if irefl.ge.ireflm				      *
c        ******************************************************
c
c------------------------------------------------------------------
c								   !
c        Input parameters					   !
c      								   !
c      x,y,z - cartesian point on the ray
c
c     the input parameter epsbnd(the small distance from the boundary)
c     is set inside this subroutine
c------------------------------------------------------------------!
c     output parameters: iflref=1 after reflection =-1 before refl.
c			cnxref,cnyref,cnzref= Nx,Ny,Nz after reflection point
c                       ibound=1 point is outside the plasma =-1 in plasma
c----------------------------------------------------------------------
      subroutine bound_xyz(x,y,z, cnx,cny,cnz, iflref,
     &x_ref,y_ref,z_ref, cnxref,cnyref,cnzref,
     &ibound,dxdt,dydt,dzdt)
c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none

c-----input
      real*8 
     &x,y,z,         !cartesian space coordinates
     &cnx,cny,cnz,   !refractive index coordinates
     &dxdt,dydt,dzdt !wave group velocity in RHS of the ray-tracing equations
c-----output
      integer
     &iflref,     !iflref=1 after reflection =-1 before refl.
     &ibound      !ibound=1 point is outside the plasma =-1 in plasma

      real*8 
     & x_ref,y_ref,z_ref,  ! reflection point coordinates
     & cnxref,cnyref,cnzref  ! Nx,Ny,Nz after reflection point
     
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
      include 'rkutta.i'
      include 'antenna.i' ! contains i_cone, i_ant
      include 'cone.i' ! contains rst(),zst(),...
      include 'onetwo.i' ! rho_bin()
c-----externals
      real*8
     & ias1r,ias2r,ias2r_Sm,bxyz,rhopsi, dense_xyz
c-----locals
      real*8
     & epsbnd,   !accuracy of point near the plasma edle LCFS
     & zp,zm,zzrm,zzrp,res,rrr, dndx,dndy,dndz,dndr,rabs,
     & edg,  r, den,wpw2,wpw_2,dwpw2dx,dwpw2dy,dwpw2dz, gradne, vvdot,
     + grpsi

      integer nr4,nz4,ipx,ipx4,idx,
     & i_reflection  !=1 it was wall/limiter reflection
                     !=0 it was not wall/limiter reflection
      nr4=nx+4
      nz4=ny+4
      ncx=nr4
      ncy=nz4

      cnxref=cnx
      cnyref=cny
      cnzref=cnz
      x_ref=x
      y_ref=y
      z_ref=z
      r= dsqrt(x*x+y*y)
      ! get rho based on model for density (model_rho_dens=1,2,4), 
      ! 2D spline of data (model_rho_dens=3),
      ! or based on magnetic flux (model_rho_dens=0,5):
      den=dense_xyz(x,y,z,1) !-> get rho (stored in one.i) 

      iflref=-1
      ibound=-1
      epsbnd=2.d-3 ![m] ! should be same as in reflect_xyz
      
c-----------------------
      ! Check: if the ray is out of xyz-grid:
      if ( (x.lt.xeqmin+epsbnd).or.(x.ge.xeqmax-epsbnd) .or.
     +     (y.lt.yeqmin+epsbnd).or.(y.ge.yeqmax-epsbnd) .or.
     +     (z.lt.zeqmin+epsbnd).or.(z.ge.zeqmax-epsbnd)  ) then
         ibound=1
         irefl=ireflm !  to stop the ray
         write(*,*)'bound_xyz: out of grid x,y,z,r',x,y,z,r
         !pause
         goto 10 !-> procedure to make a reflection/stop
      end if
c--------------------------------------------------------------YuP added
      !check that the ray reached the wall (cylinder with r=wall_rmax)
      if(r.ge. wall_rmax-epsbnd) then
         ibound=1
         !irefl=ireflm !  uncomment to stop the ray
         write(*,*)'bound_xyz:  r = wall_rmax'
         !pause
         goto 10 !-> procedure to make a reflection/stop
      endif

c--------------------------------------------------------------YuP added
      !YuP [03-2016] Added: check that rho>rho_reflect
      if((rho.ge.rho_reflect-0.01) .and. (rho.le.rho_reflect*1.5)) then
         write(*,'(a,3e13.3)')
     +      'bound_xyz: rho~rho_reflect;  x,y,z=',x,y,z
         !Note: in a mirror machine, 
         ! definition of rho=1 is quite arbitrary;
         ! That is why it is useful to define rho_reflect;
         ! it could be quite different from 1.0.
         !--------------------------------------
         !Not always the reflection should be triggered:
         ! Only when Vgroup vector is directed outward of plasma.
         ! Check the vector dot-product: ( Vgr.grad(ne) )
         wpw2= wpw_2(x,y,z,1) ! = (wpe/w)^2
         call dwpw_2(x,y,z,1, dwpw2dx,dwpw2dy,dwpw2dz)
         dndx= dwpw2dx !strictly, it is d(wpe/w)^2/dx,
         dndy= dwpw2dy !not a density itself, but it will be norm-ed
         dndz= dwpw2dz ! to the unit vector
         ! Define inverse scale length for grad((wpe/w)^2)
         dndr= sqrt(dndx*dndx+dndy*dndy+dndz*dndz)
         !one_rgrad= dndr/wpw2 ! 1/r scale of n, i.e. abs(grad(n)/n)
         ! Check that it is not zero (it can be ~0 if |grad(n)|~0)
         ! If it is ~zero, consider grad(psi) instead:
         grpsi= dsqrt(dpdxd**2 +dpdyd**2 +dpdzd**2) ! |grad(psi)|
         ! Similarly, 1/r scale of PSI is |grpsi/psid|
         rabs=sqrt(x*x+y*y+z*z)
         if( dndr*rabs .gt. wpw2*1.e-8 )then
           dndx= dwpw2dx/dndr  !grad(n) is supposed to be directed INward
           dndy= dwpw2dy/dndr 
           dndz= dwpw2dz/dndr 
         elseif( grpsi*rabs .gt. psid*1.e-8 )then
           write(*,'(a,4e13.3)')
     +      'bound_xyz: grad(n)~0;  USE psid, dpdxd,dpdyd,dpdzd=',
     +       psid, dpdxd,dpdyd,dpdzd
           !grad(n)~0; Try using grad(psi):
           dndx=-dpdxd/grpsi !grad(psi) is supposed to be directed OUTward
           dndy=-dpdyd/grpsi 
           dndz=-dpdzd/grpsi 
         else ! both grad(n) and grad(psi) are almost zero;
           ! then - use (x,y,z) direction
           write(*,'(a,2e13.3)')
     +      'bound_xyz: both grad(n) and grad(psi)~0; dndr,grpsi=',
     +       dndr,grpsi
           dndx=-x/rabs
           dndy=-y/rabs
           dndz=-z/rabs
         endif
         ! Note: We assume that grad(ne) <0 (directed inward)
         vvdot= dxdt*dndx +dydt*dndy +dzdt*dndz
         if(vvdot.lt.-1.d-33)then ! Negative ( Vgr.grad(ne) )
           ! grad(ne) is inward, Vgroup is outward : make reflection
           ibound=1
           !irefl=ireflm !  uncomment to stop the ray
           write(*,'(a,6e13.3)')
     +      'bound_xyz:  rho=rho_reflect. Vgr,grad(ne)=',
     +       dxdt,dydt,dzdt,dndx,dndy,dndz
           goto 10 !-> procedure to make a reflection/stop
         endif
      endif
      
      ! YuP Added: skip checking psi>psi_LCFS and rho>rhowall
      goto 30 ! done  
      
      if (no_reflection.eq.1) then 
          if(rho.le.rho_bin(NR))then
             i_reflection=0 ! Nothing to do
          else
          call wall_limiter_reflection_point_xyz(x,y,z,cnx,cny,cnz,
     &         i_reflection,
     &         cnxref,cnyref,cnzref,x_ref,y_ref,z_ref) !-> out
          endif
          if (i_reflection.eq.1) then  
             write(*,*) 'bound_xyz: i_reflection=1'
             !pause !!!
             irefl=irefl+1
             ibound=1
             iflref=1
          endif
          goto 30 !-> return/end
      endif
c------------------------------------------------------------------
      if(i_ox.eq.1) then ! check either rho>rhowall, or r>rst(icone)
         if ( (i_ant.le.2 .and. rho.ge.rho_bin(NR)-0.01) .or.
     +        (i_ant.eq.3 .and. r.ge.rst(i_cone))  ) then
            ibound=1
            write(*,*)'bound_xyz: rho>rhowall or r>rst.  rho,r=',rho,r
            goto 10
         else
            goto 30 !-> return/end
         endif
      endif ! i_ox.eq.1

      if(model_b.ne.0) goto 30 
      if(model_rho_dens.gt.0 .and. ixyz.eq.1) goto 30
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
      ! Check that the ray is inside flux surface
      if ((z.ge.zp-epsbnd).or.(z.le.zm+epsbnd)) then
         write(*,*)'in bound_xyz r,zm,z,zp,epsbnd',r,zm,z,zp,epsbnd
         !pause
         ibound=1
         goto 10
      endif
c-----------------------
c      res=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,0,0,r,z)
      res=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,0,0,r,z,nr4a)
      if (res.ge.(psilim*(1.d0-epsbnd))) then
         ibound=1
         write(*,*)'in bound_xyz res,psilim,epsbnd',res,psilim,epsbnd
         !pause
         goto 10
      end if
c-----------------------
  5   continue ! checking: Is rho>rhowall ?
      ! get rho based on model for density (model_rho_dens=1,2,4), 
      ! 2D spline of data (model_rho_dens=3),
      ! or based on magnetic flux (model_rho_dens=0,5):
      den=dense_xyz(x,y,z,1) !-> get rho (stored in one.i) 
      if (rho.ge.rho_bin(NR)-0.01) then ! Here: bound_xyz, for model_b=0
         ibound=1
         write(*,*)'bound_xyz: rho>rhowall  rho==',rho
         goto 10
      else
         !write(*,*)'bound_xyz: rho<1.1  Exiting...'
         !pause
         goto 30 !-> return/end
      end if
c-----------------------
 10   continue 

c      write(*,*)'in bound no_reflection,i',no_reflection,'ibound',ibound
      
      if ((no_reflection.eq.1).and.(ibound.eq.1)) then 
          ibound=-1           !initialization
          ! get rho based on model for density (model_rho_dens=1,2,4), 
          ! 2D spline of data (model_rho_dens=3),
          ! or based on magnetic flux (model_rho_dens=0,5):
          den=dense_xyz(x,y,z,1) !-> get rho (stored in one.i) 
          if(rho.le.rho_bin(NR))then
             i_reflection=0 ! Nothing to do
          else
          call wall_limiter_reflection_point_xyz(x,y,z,cnx,cny,cnz,
     &           i_reflection,
     &           cnxref,cnyref,cnzref,x_ref,y_ref,z_ref)
          endif

c          write(*,*)'after wall_limiter_reflection_point i_reflection=',
c     &    i_reflection
          !pause

          if (i_reflection.eq.1) then  
             irefl=irefl+1
             ibound=1
             iflref=1
          endif
          goto 30 !-> return/end
      endif

c-------------------------------------------------------
 100  continue
      ncx=nr4
      ncy=nz4

c----> YuP: New method: use grad(el.density) instead of grad(psi)
c----> They could be the same (in tokamaks).
      call dwpw_2(x,y,z,1, dwpw2dx,dwpw2dy,dwpw2dz)
      dndx= dwpw2dx ! strictly, it is d(wpe/w)^2/dx,
      dndy= dwpw2dy !     but we'll normalize
      dndz= dwpw2dz !     to make a unit vector along grad(ne)
      gradne= dsqrt(dndx**2 + dndy**2 + dndz**2)  ! |grad(n)|
      grpsi=  dsqrt(dpdxd**2 +dpdyd**2 +dpdzd**2) ! |grad(psi)|
      if(gradne.ne.0.d0)then
         edg=-(dndx*dxdt+dndy*dydt+dndz*dzdt)/gradne !grad(ne).er<0 (inward)
      else ! in case of flat ne(rho) profile:
         edg=(dpdxd*dxdt+dpdyd*dydt+dpdzd*dzdt)/grpsi !(inward)
      endif
      edg= edg*prmt(3)/dabs(prmt(3)) ! prmt(3)>0
      write(*,*)'bound_xyz: rho===',rho, '   edg=',edg
      write(*,'(a,3e12.4)')'x, y, z= ',x,y,z
      write(*,'(a,3e12.4)')'dpdxd,dpdyd,dpdzd',dpdxd,dpdyd,dpdzd
      write(*,'(a,3e12.4)')'dndx,dndy,dndz   ',dndx,dndy,dndz
      write(*,'(a,3e12.4)')'dxdt,dydt,dzdt   ',dxdt,dydt,dzdt
      if(edg.gt.0.d0) then !group velocity is directed outward of plasma
         irefl=irefl+1
         write(*,*)'bound_xyz->reflect_xyz:  irefl,rho===',irefl,rho
         call reflect_xyz(x,y,z,cnx,cny,cnz,cnxref,cnyref,cnzref)
         iflref=1
         !pause
      elseif(edg.lt.0.d0) then !group velocity is directed inward
         write(*,*)'bound_xyz: edg<0, Vgroup=inward.  rho===',rho
      else ! edg=0.
         irefl=ireflm ! stop the ray
         write(*,*)'bound_xyz->edg=0/stopping:  irefl,rho===',irefl,rho
      endif
 30   continue

      return
      end

c======================================================================
c======================================================================


c        *********************BOUNDC_xyz ****************
c        *                        -                           *
c        * This function gives  the ray position in           *
c        * the plasma. If the point of the ray is inside the  *
c        * limiting flux surface, leave iBOUNDc=-1, 
c        * else set iBOUNDc=1                                 *
c        ******************************************************
c
c------------------------------------------------------------------
c								   !
c        Input parameters					   !
c      								   !
c      x,y,z - cartesian coord of the ray point
c
c     the value of epsbnd(the small distance from the boundary)
c     is set inside this subroutine
c------------------------------------------------------------------!
c        Output parameter:iboundc=1 point is outside the  plasma   !				   !
c                         iboundc=-1 point is inside the  plasma   !				   !
c------------------------------------------------------------------!
      subroutine boundc_xyz(x,y,z,dxdt,dydt,dzdt, iboundc )
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
      include 'antenna.i' ! contains i_cone, i_ant
      include 'cone.i' ! contains rst(),zst(),...
      include 'onetwo.i' ! rho_bin()

      double precision
     1 ias1r,ias2r,ias2r_Sm,dxdt,dydt,dzdt,vvdot
     
      iboundc=-1
      epsbnd= 2.d-3 ![m]! 

      r= dsqrt(x*x+y*y)
      ! get rho based on model for density (model_rho_dens=1,2,4), 
      ! 2D spline of data (model_rho_dens=3),
      ! or based on magnetic flux (model_rho_dens=0,5):
      den=dense_xyz(x,y,z,1) !-> get rho (stored in one.i) 
      
c------------------------------------------------------------------
      ! Check: if the ray is out of xyz-grid:
      if ( (x.lt.xeqmin+epsbnd).or.(x.ge.xeqmax-epsbnd) .or.
     +     (y.lt.yeqmin+epsbnd).or.(y.ge.yeqmax-epsbnd) .or.
     +     (z.lt.zeqmin+epsbnd).or.(z.ge.zeqmax-epsbnd) ) then
         write(*,*)'boundc_xyz: out of grid. x,y,z,epsbnd',x,y,z,epsbnd
         iboundc=1
         !pause
         goto 10
      end if
c--------------------------------------------------------------YuP added
c      !check that the ray reached the wall (cylinder with r=wall_rmax)
      if(r.ge. wall_rmax-epsbnd) then
         iboundc=1
         write(*,*)'boundc_xyz:  r = wall_rmax'
         goto 10 !-> return/end
      endif
c------------------------------------------------------------------
      !YuP [03-2016] Added: check that rho>rho_reflect
      if((rho.ge.rho_reflect-0.01) .and. (rho.le.rho_reflect*1.5)) then
           write(*,'(a,3e13.3)')
     +      'boundc_xyz: rho~rho_reflect;  x,y,z=',x,y,z
         !Note: in a mirror machine, 
         ! definition of rho=1 is quite arbitrary;
         ! That is why it is useful to define rho_reflect;
         ! it could be quite different from 1.0.
         !--------------------------------------
         !Not always the reflection should be triggered:
         ! Only when Vgroup vector is directed outward of plasma.
         ! Check the vector dot-product: ( Vgr.grad(ne) )
         wpw2= wpw_2(x,y,z,1) ! = (wpe/w)^2
         call dwpw_2(x,y,z,1, dwpw2dx,dwpw2dy,dwpw2dz)
         dndx= dwpw2dx !strictly, it is d(wpe/w)^2/dx,
         dndy= dwpw2dy !not a density itself, but it will be norm-ed
         dndz= dwpw2dz ! to the unit vector
         ! Define inverse scale length for grad((wpe/w)^2)
         dndr= dsqrt(dndx*dndx+dndy*dndy+dndz*dndz)
         !one_rgrad= dndr/wpw2 ! 1/r scale of n, i.e. abs(grad(n)/n)
         ! Check that it is not zero (it can be ~0 if |grad(n)|~0)
         ! If it is ~zero, consider grad(psi) instead:
         grpsi= dsqrt(dpdxd**2 +dpdyd**2 +dpdzd**2) ! |grad(psi)|
         ! Similarly, 1/r scale of PSI is |grpsi/psid|
         rabs=sqrt(x*x+y*y+z*z)
         if( dndr*rabs .gt. wpw2*1.e-8 )then
           dndx= dwpw2dx/dndr  !grad(n) is supposed to be directed INward
           dndy= dwpw2dy/dndr 
           dndz= dwpw2dz/dndr 
         elseif( grpsi*rabs .gt. psid*1.e-8 )then
           write(*,'(a,4e13.3)')
     +      'boundc_xyz: grad(n)~0;  USE psid, dpdxd,dpdyd,dpdzd=',
     +       psid, dpdxd,dpdyd,dpdzd
           !grad(n)~0; Try using grad(psi):
           dndx=-dpdxd/grpsi !grad(psi) is supposed to be directed OUTward
           dndy=-dpdyd/grpsi 
           dndz=-dpdzd/grpsi 
         else ! both grad(n) and grad(psi) are almost zero;
           ! then - use (x,y,z) direction
           write(*,'(a,2e13.3)')
     +      'boundc_xyz: both grad(n) and grad(psi)~0; dndr,grpsi=',
     +       dndr,grpsi
           dndx=-x/rabs
           dndy=-y/rabs
           dndz=-z/rabs
         endif
         ! Note: We assume that grad(ne) <0 (directed inward)
         vvdot= dxdt*dndx +dydt*dndy +dzdt*dndz
         if(vvdot.lt.-1.d-33)then ! Negative ( Vgr.grad(ne) )
           ! grad(ne) is inward, Vgroup is outward : make reflection
           ibound=1
           !irefl=ireflm !  uncomment to stop the ray
           write(*,'(a,6e13.3)')
     +      'boundc_xyz:  rho=rho_reflect. Vgr,grad(ne)=',
     +       dxdt,dydt,dzdt,dndx,dndy,dndz
           goto 10 !-> return/end
         endif
      endif

c------------------------------------------------------------------
      if (no_reflection.eq.1) goto 10 !-> return/end
c------------------------------------------------------------------
      ! YuP Added: skip checking psi>psi_LCFS and rho>rhowall
      goto 10 ! done  

      if(i_ox.eq.1) then ! check either rho>rhowall, or r>rst(icone)
         if ( (i_ant.le.2 .and. rho.ge.rho_bin(NR)-0.01) .or.
     +        (i_ant.eq.3 .and. r.ge.rst(i_cone))   ) then
            iboundc=1
            write(*,*)'boundc_xyz: rho>rhowall or r>rst.  rho,r=',rho,r
            goto 10 !-> return/end
         endif
      endif ! i_ox.eq.1
      
      if(model_b.ne.0) goto 10 !-> return/end
      if(model_rho_dens.gt.0 .and. ixyz.eq.1) goto 10 
c------------------------------------------------------------------      
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
c      write(*,*)'in boundc rrr,zm,z,zp',rrr,zm,z,zp
      if ((z.gt.zp-epsbnd).or.(z.lt.zm+epsbnd)) then
         write(*,*)'in boundc_xyz rrr,zm,z,zp',rrr,zm,z,zp
         write(*,*)'in boundc_xyz zm+epsbnd,z,zp-epsbnd',
     1	 zm+epsbnd,z,zp-epsbnd
         iboundc=1
         goto 10  !-> return/end
      end if
c------------------------------------------------------------------
      nr4=nx+4
      nz4=ny+4
      ncx=nr4
      ncy=nz4
c      res=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,0,0,r,z)
      res=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,0,0,r,z,nr4a)
c      rho=rhopsi(res)
      if (res.gt.(psilim*(1.d0-epsbnd))) then
      write(*,*)'boundc_xyz res,psilim*(1-epsbnd)',res,psilim*(1-epsbnd)
         iboundc=1
         goto 10 !-> return/end
      end if
c------------------------------------------------------------------
  5   continue ! checking rho>rhowall
      ! get rho based on model for density (model_rho_dens=1,2,4), 
      ! 2D spline of data (model_rho_dens=3),
      ! or based on magnetic flux (model_rho_dens=0,5):
      den=dense_xyz(x,y,z,1) !-> get rho (stored in one.i) 
      if (rho.ge.rho_bin(NR)-0.01) then !Here: boundc_xyz, for model_b=0
         write(*,*)'boundc_xyz:  rho>rhowall;  rho===', rho
         iboundc=1
         goto 10 !-> return/end
      end if
c------------------------------------------------------------------
 10   continue
c------------------------------------------------------------------
      return
      end


c======================================================================
c======================================================================

c        ********************** outpt_xyz *******************
c        *                      -----                       *
c        *   output is used by drkgs2_xyz as an	            *
c        *   output subroutine. it has not to change the    *
c        *   values of its formal input  parameters.        *
c        *   will be transferred to the main program.       *
c        *       its formal parameters are:                 *
c        *         x,y,z, dery,dery,derz, ihlf, ndim, prmt  *
c        *                                                  *
c        ****************************************************
c        !   this code prints: values of u(ndim);	    !
c        !   creates the data for 3D code; 		    !
c        !   it writes array u(6) in file( the name of      !
c        !   file is given as  the first parameter in	    !
c        !   genray.in file),                   	    !
c        !   eps- value of Hamiltonian. 		    !
c        !   It controls the conditions: if the ray point   !
c        !   is inside the plasma or not .		    !
c        !   It controls the conditions: if the refractive  !
c        !   index less the cnmax or not .		    !
c        !   It writes the output date in mnemonic.txt file	    !
c        !   It corrects the trajectory for Hamiltonian     !
c        !   conservation 
c        !   It scatters the perpendicular refractive index
c----------------------------------------------------------------------
c         input parameters: us,u,deru,ihlf,ndim,prmt,iflagh
c         output parameters:
c            if iraystop=1  then end of the j_ray calculation
c            if iraystop=0  then continuation of the j_ray calculations
c            iflagh=(2 rays outside  the plasma after the correction;
c                    1 ray is near the plasma boundary after Runge-Kutta
c                      procedure (it is after reflection)
c                    3 ordinary situation ray is inside the plasma
c                      after correction procedure)
c            u(i) after correction and reflection
c            it prepares and writes  parameters for 3d and onetwo
c            codes
c            it changes the u():
c            1)due to correction procedure that conserves Hamiltonian D=0
c            2)due to boundary reflection in reflection points
c            3)due to n_perp scattering after reflection
c
c           For i_ox.eq.1 it calculates antenna vertex coordinates
c           x_st_ox, y_st_ox, z_st_ox, alpha_st_ox, beta_st_ox
c           and puts these coordinates into  cone_ec
c           For i_ox.eq.2 it creates OX mode conversion jump
c           in the small radial direction and calculates the 
c           OX transmission coefficient: transm_ox
c-----------------------------------------------------------------*
c           this program uses the following functions and
c           subroutines bxyz, gamma1_xyz, hamilt_xyz,  
c           bound_xyz, prep3d_xyz, write3d 
c-----------------------------------------------------------------*
      subroutine outpt_xyz(us,u,deru,ihlf,ndim,prmt,iflagh,iraystop)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'write.i'
      include 'cefield.i'
      include 'cone.i'
c-----the following variable is for N_perp scatterin procedure   
      !!! include 'scatnper.i'
c-----the following to get rma,zma for ox_conversion
      include 'three.i'

      dimension u(6),deru(6),prmt(*),pu(6)
c-----the following include is for hot dispersion with electron+ion, 
c     it has the arrays for the correction subroutines 
      include 'nperpcom.i'
      include 'ions.i'
      
      include 'output.i'
      include 'oxb.i'
      include 'antenna.i'
    
      save x_old,y_old,z_old      

c------------------------------------
      !===> For i_ox=1 case:
      !-> Specify the shape of surface which contains antenna.
      !-> The ray is traced from Xe=1 layer to this surface.
      !i_ant=1 ! Equal distance from antenna to plasma.
      !i_ant=2 ! Ellipse curve for EC cone vertex positions
               ! (requires the value of ela in antenna_surface_xyz)
      i_ant=3 ! Cylindrical surface with radius r0=rst(icone).
              ! The ray is stopped at r0=rst(icone).
c------------------------------------

      if (first) then
        first=.false.
        x_old=u(1)
        y_old=u(2)
        z_old=u(3)
      endif

      iraystop=0  ! outpt_xyz: initialize
      iflagh=3    ! outpt_xyz: initialize
c-----------------------------------------------------------

      x1=u(1)
      y1=u(2)
      z1=u(3)
      cnx1=u(4)
      cny1=u(5)
      cnz1=u(6)
      iter=0

      bmod= bxyz(x1,y1,z1) ! for wcw

      cnpar=(bx*cnx1+by*cny1+bz*cnz1)*o_bmod
     
      cnpar2=cnpar*cnpar
      cn2= cnx1*cnx1 + cny1*cny1 + cnz1*cnz1
      cnper=dsqrt(dabs(cn2-cnpar2))
c--------------------------------------------------------------     

c------------------------------------------------------------
c     if the number of the time step nstep_rk is bigger that maxsteps_rk
c     then stop ray calculations
      if(nstep_rk.gt.maxsteps_rk) then
        write(*,*)'**********nstep_rk.gt.maxsteps_rk****************'
        write(*,*)'  INCREASE nstep_rk', '   iray=',iray
        write(*,*)'  iraystop->1'
        iraystop=1
        !pause
        return
      end if
      nstep_rk=nstep_rk+1
c--------------------------------------------------------------------
c     if ray is close to resonance point then stop ray calculations
      cnmode= dsqrt(cnx1*cnx1+cny1*cny1+cnz1*cnz1)
      !cnmax=10000.d0
      cnmax=1.d+6

      if(cnmode.gt.cnmax) then
	write(*,*)'***************nmode.gt.cnmax***********'
      write(*,*)'  iraystop->1'
	iraystop=1
	return
      end if

c------------------------------------------------------------------
c      check that cold plasma X mode is close to UH resonance
c      and can be EBW conversion 
c-------------------------------------------------------------
      xe= wpw_2(u(1),u(2),u(3),1) ! rho is calc. inside
      ye= wcw(u(1),u(2),u(3),1)
      uh=dsqrt(xe+ye*ye)
      uh=xe+ye*ye
      del_uh=1.d-2
      temp_e=tempe_xyz(u(1),u(2),u(3),1) 
      rme=9.1094d-28 
      vt_e=dsqrt(2.d0*temp_e*1.6022d-9/rme) !(cm/sec)
                                             ! thermal velocity
      clight=2.99792458d10
      cnper_max_ebw=ye*(clight/vt_e)
      
      ! print-out: number of steps along ray when it crosses rho=1
      if( abs(rho-1.d0).lt.0.01 .and. irefl.eq.0) then
      write(*,'(a,f6.3,i9)')
     +  ' CROSSING rho=1. NUMBER OF STEPS(from launch): rho, nstep_rk=',
     +    rho,nstep_rk
        !pause
      endif
      
      if((((id.eq.1).or.(id.eq.2)).and.(uh.gt.1.d0).and.
     & ((uh-1.d0).lt.del_uh))
     &.and.(cnper.gt.cnper_max_ebw)) then
        write(*,*)'********************************************'
        write(*,*)'For cold plasma dispersion relation id=',id
        write(*,*)'the ray is close to upper hybrid resonance'
        write(*,*)' (uh-1.d0).lt.del_uh, uh=',uh,'del_uh=',del_uh
        write(*,*)'N>1, N=',dsqrt(cn2),'cnper=',cnper
        write(*,*)'EBW condition: cnper*(vt_e/clight)/(wce/w)=',
     &  cnper*(vt_e/clight)/ye
        write(*,*)'It can be Xmode to close to the UH resonance'
        write(*,*)'Xmode can be transformed to EBW'
        write(*,*)'For Xmode- EBW convertion hot plasma dispersion'
        write(*,*)'can be used id=6,9'
        write(*,*)'The ray calculation stopped' 
        write(*,*)'********************************************'
        write(*,*)'  iraystop->1'
        iraystop=1
	return
      endif

      if (us.gt.poldist_mx) then
	 write(*,*)'***************us.gt.poldist_mx***********'
       write(*,*)'  iraystop->1'
       iraystop=1
	 return
      end if

c     if Group Velocity, normalized to c, is greater than 1.1
c     (give 10 percent grace) then stop the ray.  [RWH: 030427].
      vgrmods= deru(1)**2+deru(2)**2+deru(3)**2

      if (vgrmods.gt.1.1  .and. id.eq.14) then
         write(*,*)
         write(*,*) '*************************************************'
         write(*,*) 'outpt_xyz: Stopping ray (for id=14 only).'
         write(*,*)' vgroup>1.1,   abs(vgroup) = ',dsqrt(vgrmods)  
         write(*,*) 'rho=',rho
         write(*,*) '*************************************************'
         write(*,*)
CENM 1Sep05 -- best to really stop the ray when using the relativistic
C    dispersion relation, otherwise it goes on and on without much progress.
C    vgrmods.gt.1.1 usually when it has nearly reached full depletion of power
C    in the ray anyway.
      	 iraystop=1
      	 !pause !!!
      end if
      
      if (vgrmods.gt.1.1) then ! print when it's above 1.1
         write(*,*)
         write(*,*) '*************************************************'
         write(*,*) 'outpt_xyz: WARNING  (Vgroup/c)^2=',vgrmods
         write(*,'(a,3e12.3)') ' dx/dt/c, dy/dt/c, dz/dt/c=', deru(1:3)
         write(*,'(a,4e12.3)') ' x, Nx,Ny,Nz=', x,cnx1,cny1,cnz1
         write(*,*)
      	 if (vgrmods.gt.2.0) then ! but stop when it's above 2.0
          write(*,*)
          write(*,*) 'outpt_xyz: STOPPING RAY  (Vgroup/c)^2=',vgrmods
          write(*,*) '*************************************************'
          write(*,*)
            iraystop=1
      	    !pause !!!
      	 endif
      	 return
      end if

c---------------------------------------------------------------------
c     change of the dispersion function near the cyclotron resonance 
      
      yj= wcw(x1,y1,z1,jy_d)

      if(iswitch.eq.1) then
c-------  change of the dispersion function and the absorption subroutines
c         near the cyclotron resonance points       
          call switch_da(yj,del_y,id,iabsorp,idswitch,iabswitch)
      endif
      
c---------------------------------------------------------
c     correction
c---------------------------------------------------------
c     the  switch off the Hamiltonian correction procedure
      if(icorrect.eq.0) goto 11
c---------------------------------------------------------
      epscor=prmt(4)
      bmod= bxyz(x1,y1,z1) !-> get b and derivs of b
      gam=  gamma1_xyz(x1,y1,z1,cnx1,cny1,cnz1)
      cnt2= u(4)**2 + u(5)**2 + u(6)**2
      eps=  hamilt_xyz(u(1),u(2),u(3),cnt2)
      eps1=eps
  
      cnpar=(bx*cnx1+by*cny1+bz*cnz1)*o_bmod
     
      cnpar2=cnpar*cnpar
      cn2= cnx1*cnx1 + cny1*cny1 + cnz1*cnz1
      cnper=dsqrt(dabs(cn2-cnpar2))

      if(dabs(eps).lt.epscor) goto 53 ! Done with correction

      if ((cnper.ge.0.d0).and.(id.eq.6)) then
c--------correction from the solution n_perp=n_perp(n_parallel)
c        Now it is for id=4,5,6,7
         
        do j=1,nbulk
          massc(j)=dmas(j)
          xc(j)= wpw_2(x1,y1,z1,j)
          yc(j)= wcw(x1,y1,z1,j)
          if(j.eq.1) yc(1)=-yc(1)
          tec(j)=tempe_xyz(x1,y1,z1,j)*1.d+3 !(eV) averaged temperature
          tpopc(j)=tpoprho(rho,j)
          vflowc(j)=vflowrho(rho,j)         
        enddo
 
        cnpar=(bx*cnx1+by*cny1+bz*cnz1)*o_bmod
        cn2= cnx1*cnx1 + cny1*cny1 + cnz1*cnz1
        accurcy=epscor 
        naccurc=5
        cnpero=dsqrt(dabs((cn2-cnpar**2))) 
        write(*,*)'output before solvnperp cnpar,cnpero',cnpar,cnpero
        ihermloc=iherm
        
        call solvnperp(nbulk,massc,xc,yc,tec,tpopc,vflowc,cnpar,id,
     *  ihermloc,accurcy,naccurc,cnprim,cnper) !-> get cnper
        write(*,*)'outpt after solvnperp cnpar,cnper',cnpar,cnper
        !write(*,*)'outpt before correct2 cnx1,cny1,cnz1',cnx1,cny1,cnz1
        
        call correct2_xyz(cnpar,cnper,cnpern,
     *  cnx1,cny1,cnz1, bx,by,bz,  cnxnew,cnynew,cnznew) 
        goto 21 ! YuP: skipping correct3()?
cYuP?        eps=1.d-8
cYuP?        itermax=20
c        write(*,*)'output before correct3 z1,r1,phi1',z1,r1,phi1
cYuP?        call correct3(cnpar,cnper,cnz1,cnr1,cm1,z1,r1,phi1,
cYuP?     .  eps,itermax,cnznew,cnrnew,rnew)
c        write(*,*)'output after correct3 rnew,cnznew,cnrnew',
c     .  rnew,cnznew,cnrnew        
cYuP?         u(2)=rnew
 21     continue
c        write(*,*)'outpt after correct2 cnx-y-znew',cnxnew,cnynew,cnznew
        u(4)=cnxnew
        u(5)=cnynew
        u(6)=cnznew
        goto 11 ! same as 53 Done with correction
        
      endif ! (cnper.ge.0.d0).and.(id.eq.6)
cSAP080829
 60   continue

c--------correction from the problem
c        sum{delt(x_i)**2}=min
c        disp_func(x+delt(x))=0 
	 iter=1
	 icor=0 
c----------------------------------------------------------------------  
ctest
         bmod=bxyz(x1,y1,z1) !-> get b and derivs of b
         gam= gamma1_xyz(x1,y1,z1,cnx1,cny1,cnz1)
         cnt2= u(4)**2 + u(5)**2 + u(6)**2
         eps=  hamilt_xyz(u(1),u(2),u(3),cnt2)

 52	 continue ! handle to repeat iteration
         dip=1.d0
         call dddrz1_xyz(u,deru) ! Here, deru(1:3) is dHamiltonian/dN
         sum= deru(1)**2 +deru(2)**2 +deru(3)**2
 12      dlambd=2.d0*eps/sum*dip
c-----------------------------------------------------------
	 do i=1,3
          pu(i)=u(i)
       enddo
 
	 do 51 i=4,6 ! cnx,cny,cnz are adjusted => Nphi will not conserve
 51      pu(i)=u(i)-0.5*dlambd*deru(i-3) ! needs work? In r-phi version,
                                         !only cnr and cnz were adjusted
cyup Another method - but no advantage?
c       call unit_vectors_perp_b(u(1),u(2),u(3), 
c     +                     uvtx, uvty, uvtz,  uvnx, uvny, uvnz) !->out
c       cnpar= (u(4)*bx + u(5)*by + u(6)*bz)*o_bmod
c       cnper_tang= u(4)*uvtx + u(5)*uvty + u(6)*uvtz !Nperp-b_tang-psi
c       cnper_norm= u(4)*uvnx + u(5)*uvny + u(6)*uvnz !Nperp-b_norm-psi
c       ! Adjust the component of N normal to psi flux surface:
c       cnper_norm= cnper_norm 
c     -          -0.5*dlambd*(deru(1)*uvnx + deru(2)*uvny + deru(3)*uvnz)
c       ! Find (Nx,Ny,Nz) from (Npar,Nper_tang,Nper_norm)
c       call nparper_to_nxyz(cnpar, cnper_tang, cnper_norm,
c     +                           uvtx, uvty, uvtz,  uvnx, uvny, uvnz,
c     +                           cnx,cny,cnz) !->out
c       pu(4)=cnx
c       pu(5)=cny
c       pu(6)=cnz
c-----------------------------------------------------------
         x1=pu(1)
         y1=pu(2)
         z1=pu(3)
         cnx1=pu(4)
         cny1=pu(5)
         cnz1=pu(6)
        ! check that the ray is within xyz-grid and within r-grid, 
        dxdt=deru(1)
        dydt=deru(2)
        dzdt=deru(3)
        call boundc_xyz(x1,y1,z1,dxdt,dydt,dzdt, iboundc) 

	 if (iboundc.eq.1) then
c           write(*,*)'in output iboundc=1 lflagh=2'
c	   write(*,*)'outside the plasma after the correction'
c---------------------------------------------------------------------
c          corrected ray point is outside the plasma
c          it is necessary to reduce the time step in the Runge -Kutta
c          procedure  and to recalculate the array u(i)
c----------------------------------------------------------------------
           !iflagh=2 !! outpt_xyz:  Not used now?
c           goto 120

c          after the correction step ray is out of plasma in up()
c          stop the correction procedure at the last correction step u()
           iflagh=3 ! outpt_xyz: 
           goto 11
	 endif

         bmod=bxyz(x1,y1,z1) !-> get b and derivs of b
         gam=gamma1_xyz(x1,y1,z1,cnx1,cny1,cnz1)

         cnt2= pu(4)**2 + pu(5)**2 + pu(6)**2
         eps=  hamilt_xyz(pu(1),pu(2),pu(3),cnt2)

       if(eps.gt.1.d120) then
           write(*,*)'in output eps.gt.1.d20 iter=',iter
           write(*,*)'Hamiltonian correction procedure stoped'
           write(*,*)'Hamiltonian=',eps
           goto 53 ! Done with correction
	 end if
  
	 if (dabs(eps).gt.dabs(eps1)) then
	     write(*,*)'in output eps.gt.eps1,dip',eps,eps1,dip             
	     dip=0.5d0*dip
	     goto 12
	 end if
         
c-------------------------------------------------------------
	 eps1=eps
c-------------------------------------------------------------
	 do 14 i=1,ndim
 14      u(i)=pu(i)

	 iter=iter+1
	 if(dabs(eps).lt.epscor) goto 53 ! Done with correction
         
	 if(iter.gt.20)then
	   write(*,*)'Hamiltonian correction procedure
     1	   made 20 iterations and stopped , Hamiltonian=',eps
	   goto 53 ! Done with correction
	 end if

	 goto 52
 53    continue

c--------------------------------------------------------------
c     end of correction
c--------------------------------------------------------------
 11   continue
c--------------------------------------------------------------
c     measure error in the dispersion realation
c     If D/(N|gradD|) > toll_hamilt stop ray calculation
c--------------------------------------------------------------

c      write(*,*)'output.f before refractive_index_relative_error'

      call refractive_index_relative_error_xyz(
     + u(1),u(2),u(3),u(4),u(5),u(6),iraystop)

c      write(*,*)'output.f after refractive_index_relative_error iraystop
c     &',iraystop

      if(iraystop.eq.1) then
         write(*,*)
         write(*,*) '*************************************************'
         write(*,*) 'outpt: Stopping ray.'
         write(*,*) ' D/(N|gradD|) > toll_hamilt'  
         write(*,*) '*************************************************'
         write(*,*)
         !pause !!!
         return
      endif
c------------------------------------------------------------------
c     Creates the jump of the ray point through the OX mode conversion
c     area where X_e=1 (N_perp=0)
c-------------------------------------------------------------------
      if (i_ox.eq.2) then

        i_call_prep3d_in_output_at_i_ox_conversion_eq_1=0 ! initialize
        
        !Vgroup can be used as 
        ! the direction for stepping across OX cutoff layer
        ! until X mode is found.
        !Save Vgroup at some point (last or nearly last) before OX jump
        ! Check the vector dot-product: ( Vgr(t).Vgr(t=0) )
        vvdot= deru(1)*wvgr_x(1) +deru(2)*wvgr_y(1) +deru(3)*wvgr_z(1)
        ! If it's negative, it means that
        ! at the present time step (just before OX jump)
        ! the vector of Vgroup(t_OX) is approximately opposite
        ! to the initial direction Vgroup(t=0).
        ! It can happen because |Vgroup| is quickly decreasing
        ! as the ray approaches the cutoff layer, and because of
        ! accuracy, it might even become slightly negative,
        ! comparing to the initial direction.
        if(was_not_ox_conversion .and. vvdot.gt.0.d0)then 
           !-> Save Vgroup only before OX jump, and only if vvdot>0
           vgr_x= deru(1) ! Vgroup/c ! store in one.i
           vgr_y= deru(2)
           vgr_z= deru(3)
        else
           vgr_x= wvgr_x(1) ! Vgroup/c at t=0
           vgr_y= wvgr_y(1)
           vgr_z= wvgr_z(1)
        endif ! for saving Vgroup
        !Alternatively: Take Vgroup values at the starting point:
        !vgr_x= wvgr_x(1) ! Vgroup/c at t=0
        !vgr_y= wvgr_y(1)
        !vgr_z= wvgr_z(1)
       
        
        if (was_not_ox_conversion) then
          eps_xe_loc=eps_xe
          
c---------calculate the coordinates of X mode cutoff point
c---------x_x,y_x,z_x,cnx_x,cny_x,cnz_x 
cyup          write(*,*)'output'
cyup          write(*,*)'before ox_conversion u(1),u(2),u(3),u(4),u(5),u(6)'
cyup     &    ,u(1),u(2),u(3),u(4),u(5),u(6)

          call ox_conversion_xyz(u(1),u(2),u(3),u(4),u(5),u(6),
     &    eps_xe_loc,
     &    x_x,y_x,z_x, cnx_x,cny_x,cnz_x, i_ox_conversion_loc)
     
          i_ox_conversion=i_ox_conversion_loc

cyup          write(*,*)'i_ox_conversion=',i_ox_conversion
cyup          write(*,*)'after ox_conversion u(1),u(2),u(3),u(4),u(5),u(6)'
cyup     &    ,u(1),u(2),u(3),u(4),u(5),u(6)
          
          if (i_ox_conversion.eq.1) then
c           calculate transmission coefficient for OX mode
c           conversion: transm_ox -> save into write.i

c--------------------------------------------------------
c           The following data are for write.i to prepare 
c           the output data for mnemonic.nc file            
              
            Y_abs=dabs(wpw_2(u(1),u(2),u(3),1))
            cn_par_optimal=dsqrt(Y_abs/(Y_abs+1)) !optimal N parallel
                                                  !for OX conversion
            bmod=bxyz(u(1),u(2),u(3)) ! get b and grad(psi)
            cnpar_ox=(bx*u(4)+by*u(5)+bz*u(6))*o_bmod ! N_parallel
                                                    ! before OX
                                                    ! conversion
            ! components of [b x grad(psi)]:                              
            b_gradpsi_x= by*dpdzd - bz*dpdyd
            b_gradpsi_y= bz*dpdxd - bx*dpdzd
            b_gradpsi_z= bx*dpdyd - by*dpdxd
            ! |[b x grad(psi)]| 
            b_grad_psi= dsqrt(b_gradpsi_x**2 
     +                      + b_gradpsi_y**2 
     +                      + b_gradpsi_z**2)
            ! N*[b x gradpsi]/|[b x gradpsi]|
            cn_b_gradpsi= ( u(4)*b_gradpsi_x 
     +                    + u(5)*b_gradpsi_y 
     +                    + u(6)*b_gradpsi_z ) / b_grad_psi
c-------------------------------------------------------
   
c-----------transmission coefficient in the point before the jump 
            call transmit_coef_ox_xyz(u(1),u(2),u(3),u(4),u(5),u(6),
     &                           transm_ox) !transm_ox -> save into write.i
c--------------------------------------------------------------------
c           write the data fo O mode  before jump
            i_call_prep3d_in_output_at_i_ox_conversion_eq_1=1 !for O-mode
            call prep3d_xyz(u,deru,iraystop)
c--------------------------------------------------------------------
            u(1)=  x_x
            u(2)=  y_x
            u(3)=  z_x
            u(4)=cnx_x
            u(5)=cny_x
            u(6)=cnz_x
            !write the data for X mode  after jump (set delpwr(X-mode))
            i_call_prep3d_in_output_at_i_ox_conversion_eq_1=2 !for X-mode
            call prep3d_xyz(u,deru,iraystop)
c--------------------------------------------------------------------
            was_not_ox_conversion=.false. !To prevent calling ox_conversion_xyz
            i_call_prep3d_in_output_at_i_ox_conversion_eq_1=0 ! reset
          endif ! i_ox_conversion.eq.1
        endif   ! was_not_ox_conversion 
      endif     ! i_ox.eq.2
c-------------------------------------------------------------------
    
 40   continue
     
      
 41   continue

      i_output_data=0 !this time step is not the output step
     
      if(us.lt.prmt(7))goto 20 ! us<prmt(7) not the output step
      
      ! Here: us>prmt(7) (or =prmt(7)): 
      i_output_data=1 !this time step is the output step
      !Prepare to save data, and also increase prmt(7): add prmt(6), 
      !so the new value of prmt(7) will be used 
      !in the above "if(..) goto 20"  during next call.
c      prmt(7)=prmt(7)+prmt(6)
      prmt(7)=(dint(us/prmt(6))+1)*prmt(6) ! here: in outpt_xyz
      
c-------------------------------------------------------
c     The work on the data preparation for the output
c     at given time steps
c---------------------------------------------------------
c     denormalization of ray point coordinates
c---------------------------------------------------------
      uout1=r0x*u(1)
      uout2=r0x*u(2)
      uout3=r0x*u(3)
      uout4=u(4)
      uout5=u(5)
      uout6=u(6)
c-----------------------------------------------------------
c      write(*,130)uout2,uout1,uout3,uout4,uout5,uout6
 130  format(3x,6(' ',e16.9))
c***************************************************
      bmod= bxyz(x1,y1,z1) ! get b and derivs of b
      gam= gamma1_xyz(u(1),u(2),u(3),u(4),u(5),u(6))
      ds2=dsin(gam)**2
      dc2=dcos(gam)**2
      call abc_xyz(u(1),u(2),u(3),ds2,dc2,ad,bd,cd)

      cnt2= u(4)**2 + u(5)**2 + u(6)**2
      eps=  hamilt_xyz(u(1),u(2),u(3),cnt2)

      cn2= cnx1*cnx1+ cny1*cny1+ cnz1*cnz1
      cnpar=(bx*cnx1+by*cny1+bz*cnz1)*o_bmod
      cnper2=cn2*ds2
      cnper=dsqrt(cnper2)
            
      xe=wpw_2(u(1),u(2),u(3),1)
      ye=wcw(u(1),u(2),u(3),1)
      uh=dsqrt(xe+ye*ye)
      uh=xe+ye*ye

      if (i_uh_switch.eq.1) then
         !change the output step prmt(6) for prmt6_uh_switch 
         if(uh.lt.uh_switch) then
           prmt(6)=prmt6_uh_switch
         else
         ! regular output step
           prmt(6)=prmt6
         endif
      endif


      if (i_power_switch_resonance.eq.1) then
         !change the output step prmt(6) for prmt6_power_switch_resonance
         do k=1,n_power_switch_resonance 
           if(dabs(ye-y_power_switch_resonance(k)).lt.
     &       del_y_power_switch_resonance) then
             prmt(6)=prmt6_power_switch_resonance
           else
           ! regular output step
             prmt(6)=prmt6
           endif
         enddo
      endif  

      do i=2,nbulk
         xi=wpw_2(u(1),u(2),u(3),i)
         yi=wcw(u(1),u(2),u(3),i)
c         write(*,*)'output i,xi,yi',i,xi,yi
      enddo      
c end test
 568  format(3x,'eps=',d13.6)
     
      call prep3d_xyz(u,deru,iraystop)

      if (iraystop.eq.1) then
	 return
      end if

      if (nrayelt.eq.nrelt) then
	 write(*,*)'***************nrayelt=nrelt***********'
        write(*,*)'  iraystop->1'
      	 iraystop=1
	 return
      end if

  20  continue
c---------------------------------------------------
c     end of output data preparation at the given time steps
c--------------------------------------------------
  30  format(1x,6(1x,e11.4))
c-----------------------------------------------------------
c     control of the reflection moment
c-----------------------------------------------------------
c     this call is for the calculation of deru().
c     subroutine bound_xyz uses deru() to detemine
c     if the ray point goes into or out the plasma.
      call rside_xyz(u,deru)
c-----------------------------------------------------------
c      write(*,*)'20 outpt before bound x,y,z,rho=',u(1),u(2),u(3),rho
      
      call bound_xyz(u(1),u(2),u(3),u(4),u(5),u(6),iflref,
     &      x_ref,y_ref,z_ref, cnxref,cnyref,cnzref,
     &      ibound,deru(1),deru(2),deru(3))

c      write(*,*)'20 outpt after bound ibound,iflref=',ibound,iflref

      if (iflref.eq.1) then
        write(*,*)'the data in reflection point before reflection'
        write(*,'(a,3e13.5)')'  x,y,z=',u(1:3)
        write(*,'(a,3e13.5)')'  cnx,cny,cnz=',u(4:6)
        bmod=bxyz(u(1),u(2),u(3))
        cnpar=(u(4)*bx+u(5)*by+u(6)*bz)/bmod
        cn=sqrt(u(4)**2 +u(5)**2 +u(6)**2)
        write(*,'(a,2e13.5)')'  Npar,N=',cnpar,cn
        write(*,*)'after reflection:'
        u(1)= x_ref
        u(2)= y_ref
        u(3)= z_ref
        write(*,'(a,3e13.5)')'  cnx,cny,cnz=',cnxref,cnyref,cnzref
        bmod=bxyz(x_ref,y_ref,z_ref)
        cnpar_ref=(cnxref*bx+cnyref*by+cnzref*bz)/bmod
        cn_ref=sqrt(cnxref*cnxref+cnyref*cnyref+cnzref*cnzref)
        write(*,'(a,2e13.5)')'  Npar,N=',cnpar_ref,cn_ref
        if( abs(cnpar_ref-cnpar).gt.1.d-8*abs(cnpar) )then
          write(*,*)'WARNING: after reflection Npar is changed!'
          !!pause
        endif
      endif

      if ((i_ox.eq.1).and.(iflref.eq.1)) then
c--------calculate the EC cone vertex coordinates
c        for the optimal OX mode conversion
           bmod=bxyz(u(1),u(2),u(3))  !-> get b and grad(psi)  
           b_x0=bx
           b_y0=by
           b_z0=bz
           dpsi_dx=dpdxd
           dpsi_dy=dpdyd
           dpsi_dz=dpdzd
           x0=u(1)
           y0=u(2)
           z0=u(3)
           cnx0=u(4)
           cny0=u(5)
           cnz0=u(6)
           icone=icone_loc
           write(*,*)'before antenna_vertex: x0,z0=',x0,z0
           write(*,*)'before antenna_vertex: cnx0,cnz0=',cnx0,cnz0
           call antenna_vertex_xyz(u(1),u(2),u(3), -u(4),-u(5),-u(6),
     &          b_x0,b_y0,b_z0, dpsi_dx,dpsi_dy,dpsi_dz,
     &          x_st,y_st,z_st, alpha_st,beta_st,icone)
c----------put the data into cone.i
           x_st_ox= x_st
           y_st_ox= y_st 
           z_st_ox= z_st
           r_st_ox= dsqrt(x_st**2 +y_st**2) 
c           phi_st_ox=phi_st
           alpha_st_ox=alpha_st
           beta_st_ox=beta_st
           write(*,*)'r_st_ox is found:',r_st_ox
           write(*,*)'x_st_ox=',x_st_ox
           write(*,*)'y_st_ox=',y_st_ox
           write(*,*)'z_st_ox=',z_st_ox
      endif ! ((i_ox.eq.1).and.(iflref.eq.1)  

      u(4)=cnxref
      u(5)=cnyref
      u(6)=cnzref

      if (iflref.eq.1) call prep3d_xyz(u,deru,iraystop)
c----------------------------------------------------------
c     call b() to calculate the small radius rho (inside b())
c     the resulting will be in common block  one.i WHY NEEDED?
c-----------------------------------------------------------
      bmod=bxyz(u(1),u(2),u(3)) ! for what?
      if (iflref.eq.1) then         
         iflagh=1 ! outpt_xyz: iflref.eq.1
      else
c---------------------------------------------------------------
      end if

      ! get rho based on model for density (model_rho_dens=1,2,4), 
      ! 2D spline of data (model_rho_dens=3),
      ! or based on magnetic flux (model_rho_dens=0,5):
      den=dense_xyz(u(1),u(2),u(3),1) !-> get rho (stored in one.i) 
      rhooldsc=rho

      if (irefl.ge.ireflm) then
         write(*,*)'irefl.ge.ireflm  iray,rho,Z=',iray,rho,u(3)
         !write(*,*)'i_output_data ',i_output_data 
         write(*,*)'  iraystop->1'
         !pause
         iraystop=1            
         if (i_output_data.eq.0) then
             call prep3d_xyz(u,deru,iraystop)
         endif
      else
         iraystop=0
      end if

 120  continue
c----------------------------------------------------
cyup      write(*,*)'end outpt_xyz'
      return
      end  ! outp_xyz

c======================================================================
c======================================================================

      double precision function b_average_xyz(psi_in)
c---------------------------------------------------------------------
c     calculates the afveraged magnetic field 
c     at the flux surface psi_in
c     It uses the function bxyz that calcultes the magnetic field.
c     It should be used after first call of subroutine rhospl
c     that calculates the cubic spline coefficients for small radius
c--------------------------------------------------------------------
c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none
      include 'param.i'
      include 'one.i'
      include 'gr.i'
      include 'rho.i'

      integer n_work
      parameter (n_work=3*npsi+1)
c-----input
      double precision  psi_in !poloidal flux
c-----external
       double precision bxyz
c-----local
      double precision pollength(npsi),bmod_av(npsi),
     &work(3*npsi+1),d2bmod_psi(npsi),tabl(3),x,y,z,r,pl,psi_l

      integer i1p(2),itabl(3),i_first,i,j

      data i_first /0/
      save i_first,d2bmod_psi,bmod_av

      if (i_first.eq.0) then
c-------------------------------------------------------------------
c        calculate the array bmod_av(npsi) of the averaged magnetic field
c        along the magnetic surfaces
c        and calculate cubic spline coefficients
c--------------------------------------------------------------------
         pollength(1)=0.d0 ! poloidal length of the magnetic surface
         z=zpsi(1,1)
         r=rpsi(1,1)
         x=r
         y=0.d0 ! ok for now; needs work
         bmod_av(1)=bxyz(x,y,z) ! get bmod
c         write(*,*)'bmod_av(1)',bmod_av(1)
         do 10 j=2,npsi             
            pollength(j)=0.d0 ! initialization of the poloidal length  
            psi_l=arpsi(j)
            bmod_av(j)=0
	    do 20 i=1,nteta
               z=0.5d0*(zpsi(j,i)+zpsi(j,i+1))
               r=0.5d0*(rpsi(j,i)+rpsi(j,i+1))
               x=r
               y=0.d0 ! ok for now; needs work
               bmod=bxyz(x,y,z) ! get bmod
               pl=dsqrt((rpsi(j,i+1)-rpsi(j,i))**2+
     1                 (zpsi(j,i+1)-zpsi(j,i))**2)
               pollength(j)=pollength(j)+pl
               bmod_av(j)=bmod_av(j)+bmod*pl
c               write(*,*)'i,pl,bmod,bmod_av(j)',i,pl,bmod,bmod_av(j)
 20         continue
            bmod_av(j)=bmod_av(j)/pollength(j) !normalize by poloidal length
c            write(*,*)'bmod_av(j),pollength(j)',bmod_av(j),pollength(j)
 10      continue
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
      itabl(1)=1
      itabl(2)=0
      itabl(3)=0
      call terp1(npsi,arpsi,bmod_av,d2bmod_psi,psi_in,1,tabl,itabl)
        b_average_xyz=tabl(1)
      return
      end


c======================================================================
c======================================================================


c        ********************** outinit_xyz *****************
c        *                      -----                       *
c        *   outinit is used in dinit to prepare output	    *
c        *   data for 3D  F-P code       		    *
c        *   outinit subroutine.                            *
c        *                                                  *
c        ****************************************************
      subroutine outinit_xyz(u)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'write.i'
      dimension u(*)
      x1=  u(1)
      y1=  u(2)
      z1=  u(3)
      cnx1=u(4)
      cny1=u(5)
      cnz1=u(6)
        nrayelt=0
      xold=u(1)
      yold=u(2)
      zold=u(3)
      return
      end



c======================================================================
c======================================================================

c        ********************cnpermuz_xyz *********************
c        * 						      *
c        * this subroutine finds the perpendicular component  *
c        * of the refractive index from the solution of the   *
c        * dispersion relation D(cnpar,cnper)=0	   using      *
c        * Mazzucato code, using full dielectric tensor with  *
c        * the hermitian and anti-hermitian parts 	      *		
c        *****************************************************
c-------------------------------------------------------------
c       input:				     !
c       cnpar - paralel (to magnetic field) component	     !
c               of the refractive index			     !
c       ioxm =1 O mode, =-1 X mode			     !
c       ihermloc=1 Hermit dielectric tenzor,	 	     !
c               =2 full	dielectric tenzor	 	     !
c       Attention! iherm is in common one.i		     !
c     	x,y,z coordinates of ray point !
c       rho      the small radius (from common one.i) 	     !
c       ioptmaz  is the option for the estimation
c                of the perpendicular refractive index       !
c       ioptmaz=1 np2  = estimation  of perp ref. index      ! 
c                (from coldm)                                !
c       ioptmaz=2   np2  = estimateion of perp ref. index    !
c                   (from input parameter cnper)             !
c       cnper - real part of the perpendicular component of  !
c               the refractive index           		     !
c       output: !
c       cnper - real part of the perpendicular component of  !
c               the refractive index           		     !
c       cnprim- imaginal part of the perpendicular component !
c               of the refractive index        		     !
c       reps(3,3)-complex dielectric tensor		     !
c                 (with antihermitian part) in commmon/eps/  !
c------------------------------------------------------------
      subroutine cnpermuz_xyz(cnpar,ihermloc,x,y,z,
     +  cnper,cnprim,ioptmaz)
        implicit integer (i-n), real*8 (a-h,o-z)
        double precision mode,mu,nz,np2,npara
        double precision ncld(2)
        include 'param.i'
        include 'one.i'
        include 'eps.i'
        dimension an2(2),icuto(2)
        double complex ceps(3,3)
        double complex cpp,sol,chamilmz

        mode=dfloat(ioxm)
        xe= wpw_2(x,y,z,1)         
        bmod=bxyz(x,y,z) !-> get b (needed for wcw)
        ye= wcw(x,y,z,1)


        te= tempe_xyz(x,y,z,1)
        mu=512.44d00/te

        npara=cnpar
cSm060225
        call n_par_min_maz(cnpar,npara)

        az2=npara*npara
        accrcy=1.d-06
        naccrcy=500

        if(ioptmaz.eq.1) then 
c---------calculation of nperp for X and O modes
c         estimation of nperp from the cold plasma 
          call coldm(npara,az2,xe,ye,an2,ncld,icuto)
          if (ioxm.eq.1) then
             np2=an2(2) !omode
          end if
          if (ioxm.eq.-1) then
             np2=an2(1) !xmode
          end if
c-------------------------------------------------------------------
c         it calcultes the solution of dispersion relation 
c         d (n_perp_real,d_nperp_im)=0 using
c         if ihermloc=2 the total Mazzucato tensor
c         if ihermloc=1 the hermitian part of the total Mazzucato tensor       
c----------------------------------------------------------------------

          call complx2(mode,xe,ye,mu,npara,np2,accrcy,naccrcy,ceps,sol,
     +    ihermloc)

          !write(*,*)'cold zero iteration cnpermuz sol',sol
          a1=dreal(sol)
          cnper=a1
          a2=dimag(sol)
          cnprim=a2
        endif !ioptmaz=1
c------------------------------------------------------------------- 
        if(ioptmaz.eq.2) then 
c---------calculation of the imaginary part cnprim of N_perpendicular
c         from the total dispersion function 
c         using the given real part cnper 
          
c---------calculation of the complex value of the dispersion function
c         from the total hermitian and anti-hermitian parts of the
c         dielectric tensor, using the real N_perp refractive index
c         cnper

          ihermlc1=2
          cpp=dcmplx(cnper,0.0d0)
          np2=cnper*cnper
          sol=cpp*cpp
          call complx1(mode,xe,ye,mu,npara,np2,sol,ihermlc1,chamilmz)
c          write(*,*)'cnpermuz 0 ihermloc=2 cpp=cnper',cpp
c          write(*,10)dreal(chamilmz),dimag(chamilmz)
 10       format('total chamilmz iherm=2 initial cpp=cnper'/
     *    2(' ',1pe15.7))
          hamimag=dimag(chamilmz)
          hamr=dreal(chamilmz)

c---------calculation of the numerical derivative from the
c         real part of the hermitian dispersion function
c         d(dispersion_function)/d(N_perp_real)
          ihermlc1=1
          step=1.0d-7  ! cnpermuz_xyz
          if(dabs(cnper).gt.1.d0) step=cnper*step
            
          cnperp=cnper+step
          cpp=dcmplx(cnperp,0.0d0)
          np2=cnperp*cnperp
          sol=cpp*cpp
          call complx1(mode,xe,ye,mu,npara,np2,sol,ihermlc1,chamilmz)
          hamrp=dreal(chamilmz)
          
          cnperm=cnper-step
          cpp=dcmplx(cnperm,0.0d0)
          np2=cnperm*cnperm
          sol=cpp*cpp
          call complx1(mode,xe,ye,mu,npara,np2,sol,ihermlc1,chamilmz)
          hamrm=dreal(chamilmz)          

          dhamrdnr=(hamrp-hamrm)/(2.d0*step)
          cnprim=-hamimag/dhamrdnr !imaginary part of N_perp

          write(*,*)'cnpermuz.f ipotmaz=2 hamimag,dhamrdnr,cnprim',
     .    hamimag,dhamrdnr,cnprim

ctest
 11       format ('chamilmz',2(' ',1pe15.7))
 12       format ('hermitian chamilmz ihermloc=1 for initial cpp=cnper'/
     *    2(' ',1pe15.7))

          ihermloc=2
          cpp=dcmplx(cnper,cnprim)
          np2=cnper*cnper
          sol=cpp*cpp
c---------this call will create dielectric tensor reps and will put this tensor into eps.i 
          call complx1(mode,xe,ye,mu,npara,np2,sol,ihermlc1,chamilmz)
c          write(*,*)'cnpermuz 2 iherml=2 cpp=(cnper,cnprim)',cpp
c          write(*,13)dreal(chamilmz),dimag(chamilmz)
 13       format('total chamilmz iherm=2 cpp=(initial_cnper,cnprim)'/
     *    2(' ',1pe15.7))
 21       format('solutition of Mazzucato solver for full dispesion'/
     *    'using the initial complexNperp=ReN_perp+i*ImN_perp'/
     *    'Im_Nperp=Im(D_full(N_par,ReN_perp)/dReD_full/dReN_perp'/
     *    'sol=',2(' ',1pe15.7))         
           write(*,21)dreal(sol),dimag(sol) 
           goto 300           

c-------- calculation cnper,cnprim from the solution of the system
c         dhamr/dnr*delnr+dhamr/dni*delni=-hamr
c         dhami/dnr*delnr+dhami/dni*delni=-hami
          ihermloc=2
          step=1.0d-7 ! cnpermuz_xyz
          cnprim=0.d0
          epsnprim=1.0d-2
          inprmax=1
          inprim=0

 100      continue

          cnperp=cnper+step*cnper
          cpp=dcmplx(cnperp,cnprim)
          np2=cnperp*cnperp
          sol=cpp*cpp
          call complx1(mode,xe,ye,mu,npara,np2,sol,ihermloc,chamilmz)
          hamrp=dreal(chamilmz)
          hamip=dimag(chamilmz)
c          write(*,*)'chamilmz',chamilmz
c          write(*,*)'hamrp,hamip',hamrp,hamip
          cnperm=cnper-step*cnper
          cpp=dcmplx(cnperm,cnprim)
          np2=cnperm*cnperm
          sol=cpp*cpp
          call complx1(mode,xe,ye,mu,npara,np2,sol,ihermloc,chamilmz)
          hamrm=dreal(chamilmz)
          hamim=dimag(chamilmz) 
c          write(*,*)'chamilmz,2*step*cnper',chamilmz,2*step*cnper
c          write(*,*)'hamrm,hamim',hamrm,hamim 
          dhamrdnr=(hamrp-hamrm)/(2.d0*step*cnper)
          dhamidnr=(hamip-hamim)/(2.d0*step*cnper)
c          write(*,*)'dhamrdnr,dhamidnr',dhamrdnr,dhamidnr
          
          cnprimp=cnprim+step
          cpp=dcmplx(cnper,cnprimp)
          np2=cnper*cnper
          sol=cpp*cpp
          call complx1(mode,xe,ye,mu,npara,np2,sol,ihermloc,chamilmz)
          hamrp=dreal(chamilmz)
          hamip=dimag(chamilmz)
c          write(*,*)'chamilmz',chamilmz
c          write(*,*)'hamrp,hamip',hamrp,hamip
          cnprimm=cnprim-step
          cpp=dcmplx(cnper,cnprimm)
          np2=cnper*cnper
          sol=cpp*cpp
          call complx1(mode,xe,ye,mu,npara,np2,sol,ihermloc,chamilmz)
          hamrm=dreal(chamilmz)
          hamim=dimag(chamilmz)
c          write(*,*)'chamilmz',chamilmz
c          write(*,*)'hamrm,hamim',hamrm,hamim 
          dhamrdni=(hamrp-hamrm)/(2.d0*step)
          dhamidni=(hamip-hamim)/(2.d0*step)
c          write(*,*)'dhamrdni,dhamidni',dhamrdni,dhamidni
          delt=dhamrdnr*dhamidni-dhamrdni*dhamidnr
          deltr=-(hamr*dhamidni-hamimag*dhamrdni)
          delti=-(dhamrdnr*hamimag-dhamidnr*hamr)
c         write(*,*)'delt,deltr,delti',delt,deltr,delti
          dcnper=deltr/delt
          cnprim=delti/delt
          cnper=cnper+dcnper
c          write(*,*)'3 dcnper,cnper,cnprim',dcnper,cnper,cnprim

          cpp=dcmplx(cnper,cnprim)
          np2=cnper*cnper
          sol=cpp*cpp
          call complx1(mode,xe,ye,mu,npara,np2,sol,ihermloc,chamilmz)
c          write(*,*)'cnpermuz 3 ihermloc=2 cpp=(cnper,cnprim)',cpp
c          write(*,14)dreal(chamilmz),dimag(chamilmz)
 14       format('total chamilmz iherm=2 cpp=(cnper,cnprim)'/
     *    2(' ',1pe15.7))
          inprim=inprim+1
          if(inprim.gt.inprmax) then
c             write(*,*)'Warning: cnpermuz inprim > inprmax'
             goto 200
          endif

          if ( dsqrt(dreal(chamilmz)**2+dimag(chamilmz)**2)
     *        .gt.epsnprim) goto 100
 200      continue
         
cendtest
        endif !iopmuz=2
c-------------------------------------------------------------------
 300    continue
	return
	end



c======================================================================
c======================================================================

c======================================================================
c======================================================================




      double complex function relativistic_nperp_xyz(x,y,z,
     & nll,cnper_tang,
     + K,iraystop) !-> out
c     calculates n_perp( and the relativistic dielectric tensor K)
c     from the relativistic electron
c     dispersion function (combined Eric Nelson-Melby - Abhay Ram)
c     if ibw=0 it calculates the O and X modes, using the cold plasma roots
c              as the initial approximation
c     if ibw=1 it calculates BW root N_perp_bw> N_perp_O and N_perp_X
c--------------------------------------------------------------------
c     INPUTS:    
c     x,y,z - space cordinates
c     nll     - parallel index of refraction n
c     cnper_tang - refractive index tangential to psi flux surface
c                  and perp to b-field (contributes to Nperp) 		       	  
c    
c     OUTPUTS:
c       function hotnperp_xyz=N_perpendiculae (double complex)
c       K double complex dielectric tensor
c       iraystop=0 the root have been found,
c       iraystop=1 the root was not found
c-----------------------------------------------------    
c      implicit integer (i-n), real*8 (a-h,o-z) 
      implicit none
      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'nperpcom.i'
c-----input
c      integer nbulk 
      double precision x,y,z
      double precision nll ! N_parallel
      double precision cnper_tang
c-----external 
      double complex hotnp
      double precision bxyz,wpw_2,wcw,tempe_xyz,tpoprho,vflowrho
      double precision ddwrap_relativistic, rtbis

cBH090202:  After pathscale compiler error:
      external ddwrap_relativistic

      external bxyz,wpw_2,wcw,tempe_xyz,tpoprho,vflowrho
 
c-----local  
      double precision mass_ar(nbulka),x_ar(nbulka),y_ar(nbulka),
     .t_av_ar(nbulka),tpop_ar(nbulka),vflow_ar(nbulka), cnx,cny,cnz
      integer j,id_old

      double precision cnper,cnprim,
     &cnper_o,cnper_x,
     &cnpar2,cntang2,cnper_02,cnper_0,xacc,cnperp_o,cnperp_x,
     &cnpermax,x1,x2
  
      integer ioxm_old,ihermloc,ioptmaz,
cSm060725
     &ntry
 
      double precision te_kev
cSm060725
      double complex compl_nper
      double complex d
      double complex eps_weiss(3,3) !eps_weiss is in Weiss system
                                    !k is in z-y plane
      
c-----output      
      double complex K(3,3)
      integer iraystop

      !write(*,*)'relativistic_nperp_xyz'
      !pause
      
      cnpar2= nll*nll
      cntang2= cnper_tang*cnper_tang + cnpar2 !cnteta*cnteta+cnphi*cnphi
      bmod=bxyz(x,y,z) !-> get b (needed for wcw)
      
      if(ibw.eq.0) then ! X and O mode cnper calculation 
c-------calculation of perpendicular refractive index cnper
c       for the given ioxm=(+/-)1 mode  using Mazzucato solver
c       ioxm is in common /one/

        ihermloc=1
        ioptmaz=1
        call cnpermuz_xyz(nll,ihermloc,x,y,z,cnper,cnprim,ioptmaz)
        write(*,*)'in function relativistic_nperp cnper =',cnper
        if(cnper.ge.0.d0) then
           iraystop=0
        else
           write(*,*)'  iraystop->1'
           iraystop=1
           write(*,*)'relativistic_nperp_xyz: cnpermuz could not
     &     find the root for mode ioxm=',ioxm
           id_old=id
           id=1            ! to use the cold plasma dispersion
   
           call cninit12_xyz(x,y,z,nll,cnper_tang,
     &                       cnx,cny,cnz,iraystop)
           id=id_old 
           cnper_02= cnx*cnx + cny*cny + cnz*cnz -nll*nll
           write(*,*)'cold plasma cnper_02=',cnper_02
           cnper_0=dsqrt(cnper_02)
           write(*,*)'cold plasma cnper_0=',cnper_0  
           stop 'in function relativistic_nperp'
        endif
      else
        ! calculation of the EBW
c-------initialization common/nperpcom/, it is in nperpcom.i
        nllc=nll
        nbulkc=nbulk
        if(nbulka.lt.nbulk) then
          write(*,*)'relativistic_nperp_xyz: nbulka.lt.nbulk'          
          write(*,*)'nbulka=',nbulka,'nbulk=',nbulk
          write(*,*)'change parameter nbulka in file param.i'  
          stop
        endif

        do j=1,nbulk
          massc(j)=dmas(j)
          xc(j)=wpw_2(x,y,z,j)
          yc(j)=wcw(x,y,z,j)
          if(j.eq.1) yc(1)=-yc(1) ! negative Y=(omega_ce/omega)
                                  ! for electrons
          tec(j)=tempe_xyz(x,y,z,j)*1.d+3 !(eV) averaged temperature
          tpopc(j)=tpoprho(rho,j)
          vflowc(j)=vflowrho(rho,j)         
        enddo

        xacc=1.e-5 !1.d-14  ! accuracy of the root calculations
        ntry=200


c-------calculate X and O modes using Mazucato solver 
c       to find the boundary for EBW root calculations 
        ihermloc=1
        ioptmaz=1
        ioxm_old=ioxm

c-------O mode calculations 
        ioxm=1
        call cnpermuz_xyz(nll,ihermloc,x,y,z,cnper,cnprim,ioptmaz)
        write(*,*)'in function relativistic_nperp ioxm,cnper',ioxm,cnper
         
        if(cnper.gt.0.d0) then
           cnperp_o=cnper
        else !cnper =<0
           cnperp_o=0.d0
c----------check the cnper root using cold plasma: cnper_0
           write(*,*)'in function relativistic_nperp cnpermuz could not
     &     find the root for mode ioxm=',ioxm
           id_old=id
           id=1            ! to use the cold plasma dispersion        
           call cninit12_xyz(x,y,z,nll,cnper_tang,
     &                   cnx,cny,cnz,iraystop)
           id=id_old 
           cnper_02= cnx*cnx + cny*cny + cnz*cnz -nll*nll
           write(*,*)'cold plasma cnper_02=',cnper_02
           cnper_0=dsqrt(cnper_02)
           write(*,*)'cold plasma cnper_0=',cnper_0  
        endif

c-------X mode calculations   
        ioxm=-1
        call cnpermuz_xyz(nll,ihermloc,x,y,z,cnper,cnprim,ioptmaz)
        write(*,*)'in function relativistic_nperp ioxm,cnper',ioxm,cnper
         
        if(cnper.gt.0.d0) then
           cnperp_x=cnper
        else !cnper =<0
           cnperp_x=0.d0
c----------check the cnper root using cold plasma: cnper_0
           write(*,*)'in function relativistic_nperp cnpermuz could not
     &     find the root for mode ioxm=',ioxm
           id_old=id
           id=1            ! to use the cold plasma dispersion        
           call cninit12_xyz(x,y,z,nll,cnper_tang,
     &                   cnx,cny,cnz,iraystop)
           id=id_old 
           cnper_02= cnx*cnx + cny*cny + cnz*cnz -nll*nll
           write(*,*)'cold plasma cnper_02=',cnper_02
           cnper_0=dsqrt(cnper_02)
           write(*,*)'cold plasma cnper_0=',cnper_0  
        endif  

        ioxm=ioxm_old

        cnpermax=dmax1(cnperp_x,cnperp_o)
        x1 = dlog(cnpermax+1.d-2)  !nperp
        x2 = x1+dlog(1.d0)         !nperp

c       expand the range of x2 until the first root is found
        iraystop=0
        do j = 1, ntry
          if (j.eq.ntry) then
            write(*,*)'relativistic_nperp could not find ebw root'
            write(*,*)'  iraystop->1'
            iraystop=1
            return
          endif
          if(ddwrap_relativistic(x2)*ddwrap_relativistic(x1).lt.0) then
            go to 12
          endif
          x1 = x2
          if(ddwrap_relativistic(x2)*ddwrap_relativistic(x1).ge.0) then
            x2 = x2 +dlog(2.0D0)
          endif   
        enddo
 12     continue

        cnper=dexp(rtbis(ddwrap_relativistic,x1,x2,xacc))
        write(*,*)'in relativistic_nperp cnper=',cnper  

        compl_nper=dcmplx(cnper,0.d0)
        iraystop=0
        write(*,*)'in relativistic_nperp compl_nper=',compl_nper
        
        te_kev = tec(1)*1.d-3 !electron temperature [kev]
        call Disp_Ram(te_kev,nllc,xc(1),dabs(yc(1)),
     +                compl_nper,K,d) !K is in z-y plane 

      endif   !ibw=1 BW

      relativistic_nperp_xyz= cnper

      return
      end
c======================================================================
c======================================================================

      double complex function hotnperp_xyz(x,y,z,nll,cnper_tang,
     +                                     K,iraystop) !->out
c     calculates n_perp( and the hot dielectric tensor K)
c     from the hot non-relativistic electron+ions
c     dispresion function.
c     if ibw=0 it calculates the O and X modes , using the cold plasma roots
c              as the initial approximation
c     if ibw=1 it calculates BW root N_perp_bw> N_perp_O and N_perp_X
c--------------------------------------------------------------------
c     INPUTS:    
c     x,y,z - space cordinates
c     nll     - parallel index of refraction n
c     cnper_tang - refractive index tangential to psi flux surface
c                  and perp to b-field (contributes to Nperp) 		       	  
c    
c     OUTPUTS:
c       function hotnperp_xyz=N_perpendicular (double complex)
c       K double complex dielectric tensor
c       iraystop=0 the root have been found,
c       iraystop=1 the root was not found
c-----------------------------------------------------    
      implicit integer (i-n), real*8 (a-h,o-z) 
      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'nperpcom.i'

c-----input
      double precision x,y,z
      double precision nll ! N_parallel

c-----external 
      double complex hotnp
      double precision bxyz,wpw_2,wcw,tempe_xyz,tpoprho,vflowrho
      external hotnp,bxyz,wpw_2,wcw,tempe_xyz,tpoprho,vflowrho
c-----local  
     
      double precision mass_ar(nbulka),x_ar(nbulka),y_ar(nbulka),
     .t_av_ar(nbulka),tpop_ar(nbulka),vflow_ar(nbulka)
      integer j,id_old

c-----output      
      double complex K(3,3)
      integer iraystop

      !write(*,*)'hotnperp_xyz'
      !pause

      cnpar2=nll*nll
      cntang2= cnper_tang*cnper_tang + cnpar2 !cnteta*cnteta+cnphi*cnphi
      bmod=bxyz(x,y,z) !-> get b (needed for wcw)

c-----calculation of two roots N_perp**2(N_parallel)
c     from the cold electron+ions plasma dispersion function
      id_old=id
      id=1            ! to use the cold plasma dispersion
      cnper2p=-1.d0
      cnper2m=-1.d0 
      call npernpar_xyz(x,y,z,nll, cnper2p,cnper2m)  !-> cnper2p,cnper2m
      id=id_old      
      write(*,*)'hotnperp_xyz: nll,cnper2p,cnper2m',
     .nll,cnper2p,cnper2m
      write(*,*)'hotnperp_xyz: ibw=',ibw
      
      iraystop=0

c-----choosing of the cold mode N_perp=cnper_0 accoding to ioxm parameter      
      if(ibw.eq.0) then   
        id_old=id
        id=1            ! to use the cold plasma dispersion
        
        call cninit12_xyz(x,y,z,nll,cnper_tang,
     .               cnx,cny,cnz,iraystop) !->out

        id=id_old  

        if(iraystop.eq.0) then
          cnper_0=dsqrt(cnx*cnx+cny*cny+cnz*cnz -nll*nll)
          write(*,*)'hotnperp_xyz: after cninit12 cnper_0,cnper_0**2',
     *    cnper_0,cnper_0**2
        else
          write(*,*)'hotnperp_xyz: iraystop=1, could not find 
     .    the root for the cold plasma mode'
        endif
      endif

c-----initialization common/nperpcom/, it is in nperpcom.i
      nllc=nll
      nbulkc=nbulk
      if(nbulka.lt.nbulk) then
        write(*,*)'hotnperp_xyz: nbulka.lt.nbulk'          
        write(*,*)'nbulka=',nbulka,'nbulk=',nbulk
        write(*,*)'change parameter nbulka in file param.i'  
        stop
      endif

      do j=1,nbulk
        massc(j)=dmas(j)
        xc(j)=wpw_2(x,y,z,j)
        yc(j)=wcw(x,y,z,j)
        if(j.eq.1) yc(1)=-yc(1) ! negative Y=(omega_ce/omega) for electrons
        tec(j)=tempe_xyz(x,y,z,j)*1.d+3 !(eV) averaged temperature
        tpopc(j)=tpoprho(rho,j)
        vflowc(j)=vflowrho(rho,j)         
      enddo
      
      hotnperp_xyz= hotnp(nbulk,ibw,cnper_0,cnper2p,cnper2m,K,iraystop)
      
      return
      end

c======================================================================
c======================================================================

      subroutine calculate_hot_nperp_roots_xyz(x,y,z,cnpar,
     &n_hot_roots,N_perp_root_ar)
c---------------------------------------------------------------------
c     calculates hot plasma dispersion function roots
c     N_perp and polarization at given point z,r,phi
c     and cnpar=N_parallel
c----------------------------------------------------------
      implicit none
      include 'param.i'
      include 'one.i'      
c-----input
      double precision
     &x,y,z,           !space coordinates
     &cnpar              !N_parallel
         
c-----output
      integer n_hot_roots              !number of hot roots
      double precision 
     &N_perp_root_ar(n_hot_roots_a)    !hot plasma roots
      double complex e_xyz_nperp_root_ar(3,n_hot_roots_a) ! electric filed 
                                                      ! for ecach root
                               ! e_nperp_root_ar(1,k) E_x compomnent
                               ! e_nperp_root_ar(2,k) E_y compomnent
                               ! e_nperp_root_ar(3,k) E_z compomnent

      double precision  e_mode_root_ar(3,n_hot_roots_a)  !electric field
                                                         !polarization 
                                                         ! E+,E-,E_parallel
                                                         !at each N_perp_root(i)]
                               !e_mode_root_ar(1,k)      ! E+
                               !e_mode_root_ar(2,k)      ! E-
                               !e_mode_root_ar(3,k)      ! E||


c-----externals
      double precision bxyz,wpw_2,wcw,tempe_xyz,tpoprho,vflowrho
c-----locals
      integer j
      double precision x_ar(nbulka),y_ar(nbulka),T_ar_ev(nbulka),     
     &tpop_ar(nbulka),vflow_ar(nbulka)![sm/sec]

      character*60  text
      character*6 format
      character*2 format1

      bmod=bxyz(x,y,z) ! for wcw
      do j=1,nbulk
        x_ar(j)=wpw_2(x,y,z,j)
        y_ar(j)=wcw(x,y,z,j)
        if(j.eq.1) y_ar(1)=-y_ar(1) ! negative Y=(omega_ce/omega)
                                    ! for electrons
        T_ar_ev(j)=tempe_xyz(x,y,z,j)*1.d+3 !(eV) averaged temperature
        tpop_ar(j)=tpoprho(rho,j)
        vflow_ar(j)=vflowrho(rho,j)         
      enddo

      call hot_roots_solver(nbulk,T_ar_ev,tpop_ar,vflow_ar,
     &x_ar,y_ar,cnpar,
     &cN_perp_root_max,n_points_root,  
     &N_perp_root_ar,n_hot_roots,
     &e_xyz_nperp_root_ar,e_mode_root_ar)

      write(*,*)'Number of roots:n_hot_roots',n_hot_roots 
      open(10,file='hot_roots_at_given_point.dat')
      format='d21.15'
      format1='i2'
      text='(1X,"x=",'//format//')'
      write(10,fmt=text)x
      text='(1X,"y=",'//format//')'
      write(10,fmt=text)y
      text='(1X,"z=",'//format//')'
      write(10,fmt=text)z
      text='(1X,"cnpar=",'//format//')'
      write(10,fmt=text)cnpar
      text='(1X,"The code searched roots at the interval")'
      write(10,fmt=text)
      text='(1X,"0 < N_perp < cN_perp_root_max ")'
      write(10,fmt=text)
      text='(1X,"using n_points_root uniform N_perp mesh")'
      write(10,fmt=text)
      text='(1X,"cN_perp_root_max=",'//format//')'
      write(10,fmt=text)cN_perp_root_max
      text='(1X," n_points_root=",'//format1//')'
      write(10,fmt=text)n_points_root

      do j=1,n_hot_roots
        write(*,*)'===============hot root================='
        write(*,*)'j,N_perp_root_ar(j)',j,N_perp_root_ar(j)
        write(*,*)'EX e_xyz_nperp_root_ar(1,j)',e_xyz_nperp_root_ar(1,j)
        write(*,*)'EY e_xyz_nperp_root_ar(2,j)',e_xyz_nperp_root_ar(2,j)
        write(*,*)'EZ e_xyz_nperp_root_ar(3,j)',e_xyz_nperp_root_ar(3,j)
        write(*,*)'E+ e_mode_root_ar(1,j)',e_mode_root_ar(1,j)
        write(*,*)'E- e_mode_root_ar(2,j)',e_mode_root_ar(2,j)
        write(*,*)'E|| e_mode_root_ar(3,j)',e_mode_root_ar(3,j)
        write(*,*)'========================================'
c-----------------------------------------------------------------------
c       write roots and polarization to the file:
c       hot_roots_at_given_point.dat
c------------------------------------------------------------------
        text='(1X,"========================================")'
        write(10,fmt=text)
        text='(1X,"root number=",'//format1//')'
        write(10,fmt=text)j
        text='(1X,"root N_perpendicular=",'//format//')'
        write(10,fmt=text)N_perp_root_ar(j)
        text='(1X,"ReEx=",'//format//')'
        write(10,fmt=text) dreal(e_xyz_nperp_root_ar(1,j))
        text='(1X,"ImEx=",'//format//')'
        write(10,fmt=text) dimag(e_xyz_nperp_root_ar(1,j))
        text='(1X,"ReEy=",'//format//')'
        write(10,fmt=text) dreal(e_xyz_nperp_root_ar(2,j))
        text='(1X,"ImEy=",'//format//')'
        write(10,fmt=text) dimag(e_xyz_nperp_root_ar(2,j))
        text='(1X,"ReEz=",'//format//')'
        write(10,fmt=text) dreal(e_xyz_nperp_root_ar(3,j))
        text='(1X,"ImEz=",'//format//')'
        write(10,fmt=text) dimag(e_xyz_nperp_root_ar(3,j))
        text='(1X,"|E+|=",'//format//')'
        write(10,fmt=text)e_mode_root_ar(1,j)
        text='(1X,"|E-|=",'//format//')'
        write(10,fmt=text)e_mode_root_ar(2,j)
        text='(1X,"|E|||=",'//format//')'
        write(10,fmt=text)e_mode_root_ar(3,j)

      enddo

      close(10)

      return
      end

c======================================================================
c======================================================================
      
      
c        **********************cninit3_xyz ********************
c        *                        -                           *
c        * It solves the dispersion relation N=N(n_par)       *
c        * for Appleton-Hartree disperstion relation
c        * Then subroutine calculates the initial components  *
c        * of the refractive index  cnx,cny,cnz        	      *
c        ******************************************************
c
c------------------------------------------------------------------
c					                          !
c        input parameters				          !
c        x,y,z,cnpar,
c        cnper_tang - refractive index tangential to psi flux surface
c                     and perp to b-field (contributes to Nperp) 		       	  
c        istart from common block/ one/                           !
c        output parameters				          !
c        u(4)=cnx,u(5)=cny,u(6)=cnz 			          !
c        iraystop=1 end ray calculation                           !
c------------------------------------------------------------------
c        it uses the following functions and subroutines          !
c        ias1r,bxyz,wcw,wpw_2,gamma1_xyz,s_xyz,abc_xyz,hamilt_xyz,                     !
c------------------------------------------------------------------
      subroutine cninit3_xyz(x,y,z,cnpar,cnper_tang,
     1                  cnx,cny,cnz,iraystop)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'

      save
      double complex cmplnper      

      iraystop=0
      pi=4*datan(1.d0)
      cnpar2=cnpar*cnpar
      cntang2= cnper_tang*cnper_tang + cnpar2 !cnteta*cnteta+cnphi*cnphi
      bmod=bxyz(x,y,z) !-> get components of b and derivs of b

c-----------------------------------------------------------------
      iroot=0
c------------------------------------------------------------------
c     epsmode is the accuracy  for selection mode (O or X) on
c             the plasma boundary
c------------------------------------------------------------------
      epsmode=1.d-4 !0.0000001d0
c--------------------------------------------------------------------
c     initial condition in plasma boundary
c     for cn2 from cnpar by using dispersion relation
c     f*n**4 +g*n**2+w=0
c     in the following form n**2=(-g+ioxm*sqrt(g**2-4*f*w))/(2*f)
c------------------------------------------------------------------
c------------------------------------------------------------------
c     Appleton - Hartree dispersion relation
c
      xi=wpw_2(x,y,z,1)
      yi=wcw(x,y,z,1)
      py2=yi*yi
      py4=py2*py2
      px=1.d0-xi
      px2=px*px
c------------------------------------------------------------------
      pyp=xi/(1.d0+yi)
      f1e=1.d0
      f0e=-pyp
      g1e=cnpar2*(-xi)+pyp+xi-2.d0
      g0e=cnpar2*pyp+xi*(1.-pyp)+pyp*(1.d0-xi)
      w1e=cnpar2*(xi-pyp)+(1.d0-pyp)*(1.d0-xi)
      w0e=cnpar2*(pyp*(1.d0-xi)-xi*(1.d0-pyp))-(1.d0-xi)*(1-pyp)*xi
      dele=1.d0-yi
      fd=f1e*dele+f0e
      gd=g1e*dele+g0e
      wd=w1e*dele+w0e
      detin=gd**2-4.d0*fd*wd
      if (detin.lt.0d0) then ! negative discriminant
         !YuP[11-2016] adjustment for a very small neg.discriminant
         if(abs(detin).lt. 1.d-8*(gd**2+4.d0*fd*wd) ) then
           ! Could be a very small negative value, from rounding.
           ! Then, set it to zero:
           detin=0.d0
         else ! large negative discriminant: no wave here.
           write(*,*)'cninit3_xyz: detin,gd**2,Xe=',
     +                             detin,gd**2,xe
           write(*,*)' 3 in cninit3_xyz detin<0'
           write(*,*)'  iraystop->1'
           iraystop=1
           return
         endif
      end if

c     cn2=(-gd+ioxm*dsqrt(detin))/(2.d0*fd)
      cn2p=(-gd+dsqrt(detin))/(2.d0*fd)
      cn2m=(-gd-dsqrt(detin))/(2.d0*fd)
      WRITE(*,*)'apl cninit cn2p,cn2m',cn2p,cn2m
      if((cn2m.lt.0.d0).and.(cn2p.lt.0d0)) then
            write(*,*)'in cninit2 two roots of the dispersion < 0'
            write(*,*)'cn2m,cn2p',cn2m,cn2p
            write(*,*)'the given wave can not exist in plasma'
            write(*,*)'  iraystop->1'
	    iraystop=1
	    return
      endif
10    iroot=iroot+1
      if(iroot.eq.1) then
        cn2=cn2p
	if(cn2.lt.cntang2)then
           write(*,*)'in cninit3_xyz cn2p.lt.cntang2'
           go to 10
        end if
      else
        cn2=cn2m
        if(cn2.lt.cntang2)then
           write(*,*)'in cninit3_xyz cn2m.lt.cntang2'
           write(*,*)'the given mode can not exist in plasma'
           write(*,*)'  iraystop->1'
           iraystop=1
           !pause
           return
        end if
      endif
c     write(*,*)'in cninit cn2=',cn2
      cnrho2=cn2-cntang2
      cnrho=dsqrt(cnrho2)
c------------------------------------------------------------------
      cnper_norm=cnrho ! found above
      call cnxyz(x,y,z, cnpar, cnper_tang, cnper_norm,
     +                 cnx,cny,cnz) !->out
c------------------------------------------------------------------YuP
      gam=gamma1_xyz(x,y,z,cnx,cny,cnz)
      ds=dsin(gam)
      dc=dcos(gam)
      ds2=ds*ds
      dc2=dc*dc
      ds4=ds2*ds2
c--------------------------------------------------------------------
c     control that cn2 and gam are the solution of the dispersion
c     relation  n**2=(-b+ioxm*sqrt(b**2-4*a*c))/(2*a)
c--------------------------------------------------------------------
      sqrdet=dsqrt(py4*ds4+4.d0*py2*px2*dc2)
      pz=2.d0*px-py2*ds2+ioxm*sqrdet
      cn2new=1.d0-2.d0*xi*px/pz
c     write(*,*)'cn2new=',cn2new
      if (iroot.eq.1) then
	dnp=dabs(cn2-cn2new)/dabs(cn2+cn2new)
	if (dnp.gt.epsmode)then
	   goto 10
	end if
      else
	dnm=dabs(cn2-cn2new)/dabs(cn2+cn2new)
	if (dnm.gt.epsmode)then
	          write(*,*)'cninit3_xyz:'
           write(*,*)'the given mode can not exist in plasma'
           !pause
         write(*,*)'  iraystop->1'
	   iraystop=1
	   return
	end if
      end if
      goto 111
      
  111 continue
      return
      end
      

c======================================================================
c======================================================================
c        **********************cninit_xyz *********************
c        *                        -                           *
c        * It solves the dispersion relation N=N(n_par)       *
c        * Then subroutine calculates the initial components  *
c        * of the refractive index N: cnx=Nx,cny=Ny,cnz=Nz    *
c        ******************************************************
c
c------------------------------------------------------------------
c					                          !
c        input :				          !
c        x,y,z,cnpar,
c        cnper_tang - refractive index tangential to psi flux surface
c                     and perp to b-field (contributes to Nperp) 		       	  
c        if i_n_poloidal=4 cnpar=N_parallel will be calculated
c        inside this subroutine using cnteta,cnphi
c        istart from common block/ one/                           !
c        output :				          !
c        u(4)=cnx,u(5)=cny,u(6)=cnz 			          !
c        iraystop=1 is switch to stop the ray calculation         !
c------------------------------------------------------------------
c        it uses the following functions and subroutines          !
c        ias1r,  bxyz ,wcw,wpw_2,gamma1_xyz,s_xyz,abc_xyz,hamilt_xyz                     !
c------------------------------------------------------------------
      subroutine cninit_xyz(x,y,z, cnpar,cnper_tang,
     1                  cnx,cny,cnz, iraystop) !->out
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'      
      include 'grill.i'
      
      save
      double complex hotnperp_xyz,cmplnper,relativistic_nperp_xyz
      double precision bxyz
      external hotnperp_xyz, relativistic_nperp_xyz, bxyz      

c-----for hot plasma roots
      double precision 
     &N_perp_root_ar(n_hot_roots_a)     !hot plasma roots
     
      double complex K(3,3)

      !write(*,*)'cninit_xyz'
      !pause

      iraystop=0
      pi=4*datan(1.d0)
      cnpar2=cnpar*cnpar
      cntang2=cnper_tang*cnper_tang + cnpar2 !cnteta*cnteta+cnphi*cnphi
      bmod=bxyz(x,y,z) !-> get components of b and derivs of b

      if(istart.eq.2) then !grill conditions
      
        if (i_n_poloidal.eq.3) then
c----------refractive index is specified by N_parallel 
c          and ksi_nperp the angle between grad(psi) and_ N_perpendicular
c
c          set zero values for following variables
c          to get the solution of the dispersion function N_perp(N_parallel)
c          In this case N_perp=N
           cntang2=0.d0
        endif

        if (i_n_poloidal.eq.4) then ! needs work?
          stop 'i_n_poloidal=4   Not setup yet'
c---------refractive index is specified by N_toroidal and N_poloidal
c         calculation of the parallel refractive index component
c         N_parallel=cnpar
cyup          gradpsi= dsqrt(dpdxd*dpdxd+dpdyd*dpdyd+dpdzd*dpdzd)
cyup          b_teta=(bz*dpdrd-br*dpdzd)/gradpsi !poloidal magnetic field
cyup          cnpar=(cnphi*bphi+cnteta*b_teta)*o_bmod !N_parallel
cyup          write(*,*)'cninit.f cnteta,cnphi,cnpar',cnteta,cnphi,cnpar
        endif

      endif !istart.eq.2
c-----------------------------------------------------------------
      iroot=0
c--------------------------------------------------------------------
c     the initial condition in plasma boundary
c     for cn2 from cnpar by using dispersion relation
c     f*n**4 +g*n**2+w=0
c     in the following form n**2=(-g+ioxm*sqrt(g**2-4*f*w))/(2*f)
c------------------------------------------------------------------
c     Appleton - Hartree dispersion relation
c
      if (id.eq.3) then
         call cninit3_xyz(x,y,z,cnpar,cnper_tang,
     1                  cnx,cny,cnz,iraystop)
         cnper=dsqrt(cnx**2 +cny**2 +cnz**2 -cnpar**2)
         goto 111
      end if !id=3
      
c-----------------------------------------------------------------
c     initial condition in plasma boundary
c     for cn2 from cnpar by using dispersin relation
c     f*n**4 +g*n**2+w=0
c     in the following form n**2=(-g+ioxm*sqrt(g**2-4*f*w))/(2*f)
c------------------------------------------------------------------
c     cold plasma dispresion relation
c------------------------------------------------------------------
      if ((id.eq.1).or.(id.eq.2)) then
         call cninit12_xyz(x,y,z,cnpar,cnper_tang,
     1                  cnx,cny,cnz,iraystop)
         if (iraystop.eq.1)then
           write(*,*)'the given conditions did not give root' 
           return
         endif
         cnper=dsqrt(cnx**2 +cny**2 +cnz**2 -cnpar**2)
         goto 111
      end if !id=1 or id=2

      if(id.eq.6) then
c-----------------------------------------------------------------
c       initial condition on plasma boundary
c       for cnper2 from cnpar by using the hot dispersion relation
c-----------------------------------------------------------------
        if (i_look_roots.eq.2) then     
          call calculate_hot_nperp_roots_xyz(x,y,z,cnpar,
     &    n_hot_roots,N_perp_root_ar)
          cnper=N_perp_root_ar(k_hot_root)
          cnprim=0.d0
        else
          ihermloc=iherm
          write(*,*)'cninit_xyz before hotnperp_xyz. cnpar,cnper_tang=',
     +                                               cnpar,cnper_tang
          cmplnper=hotnperp_xyz(x,y,z,cnpar,cnper_tang,K,iraystop)
          write(*,*)'cninit_xyz after hotnperp_xyz cmplnper',cmplnper
          if (iraystop.eq.1)then
            write(*,*)'the given conditions did not give hot root' 
            return
          endif
          cnper=dreal(cmplnper)
          cnprim=dimag(cmplnper)            
        endif
      endif !id=6

      if (id.eq.14) then
c--------------------------------------------------------------------
c        id=14 Abhay Ram relativistic electron dispersion function
c----------------------------------------------------------------------
         cmplnper=relativistic_nperp_xyz(x,y,z,cnpar,cnper_tang,K,
     &                              iraystop)
         cnper=dreal(cmplnper)
         cnprim=dimag(cmplnper)
      endif ! (id.eq.14)

c--------------------------------------------------------------------
c-----calculations of initial values cnx,cny,cnz from cnper for id=6    
      cn2=cnper**2+cnpar**2   
      cnrho2=cn2-cntang2
      if (cnrho2.lt.0.d0) cnrho2=1.d-16
      cnrho=dsqrt(cnrho2)

c---------------------------------------------------------------------
      cnper_norm=cnrho ! found above
      call cnxyz(x,y,z, cnpar, cnper_tang, cnper_norm,
     +                 cnx,cny,cnz) !->out
c------------------------------------------------------------------YuP
      gam=gamma1_xyz(x,y,z,cnx,cny,cnz)
      cnt2= cnx**2 + cny**2 + cnz**2
      ddd=hamilt_xyz(x,y,z,cnt2)

  111 continue

c----------------------------------------------------------------
c     for the grill conditions it will use the different grill types
c---------------------------------------------------------------
      if(istart.eq.2) then !grill conditions
      
         if (i_n_poloidal.eq.2) then ! known: Npar, N_theta_pol, Nperp
             ! Calculate Nper_tang, Nper_norm (==cnrho)
             gradpsi=dsqrt(dpdxd*dpdxd+dpdyd*dpdyd+dpdzd*dpdzd)
             r= dsqrt(x*x+y*y)
             ! Bpoloid = (b.[grad(psi) X e_phi]) / |gradpsi| 
             ! e_phi= {-y/r      ,   x/r     ,     0    }
             ! Find poloidal and toroidal magnetic field
             b_teta=(bz*(dpdxd*x+dpdyd*y)-(bx*x+by*y)*dpdzd)/(r*gradpsi)
             bphi=  (-bx*y +by*x)/r
             if (bphi .ne. 0.d0) then
                cnper_tang_full= (n_theta_pol*bmod -cnpar*b_teta)/bphi
             else ! bphi=0.  In this case Npar== +/- N_theta_pol
                cnper_tang_full= cnper_tang ! take from input value
             endif
             cnrho2= cnper**2 - cnper_tang_full**2
             if (cnrho2 .ge. 0.d0) then
                cnrho_full= dsqrt(cnrho2)
             else
                write(*,*) 'cninit_xyz  i_n_poloidal=2 case
     &          cnper<cnper_tang_full.   Impossible to find N_rho'
                write(*,*)'  iraystop->1'
                iraystop=1
                return                
             endif
         endif !i_n_poloidal=2

         if (i_n_poloidal.eq.4) then ! known: Nper_tang (?was? expressed
             ! through given N_toroidal, N_poloidal);  also known: Nper
             ! Calculate N_rho
             cnper_tang_full= cnper_tang ! from input ! needs work/check?
             cnrho2= cnper**2 - cnper_tang_full**2
             if (cnrho2 .ge. 0.d0) then
                cnrho_full= dsqrt(cnrho2)*i_vgr_ini
             else
                write(*,*) 'cninit_xyz  i_n_poloidal=4 case
     &          cnper<cnper_tang_full.   Impossible to find N_rho'
                write(*,*)'  iraystop->1'
                iraystop=1
                return                
             endif
         endif ! i_n_pol=4

         if((i_n_poloidal.ne.1).and.(i_n_poloidal.ne.3)) then 
            ! Find (Nx,Ny,Nz) from known (Npar, Nper_tang, Nper_norm)
             cnper_norm=cnrho_full ! found above
             cnper_tang=cnper_tang_full
             call cnxyz(x,y,z, cnpar, cnper_tang, cnper_norm,
     +                 cnx,cny,cnz) !->out
         endif

      endif !(istart.eq.2) grill conditions


      return
      end
      

c======================================================================
c======================================================================
      subroutine rho_ini_LHFW_xyz(x,y,z, cnpar,
     & i_n_poloidal, n_theta_pol, n_toroidal,
     & rho_ini, x_ini, y_ini, z_ini, cnx_ini, cny_ini, cnz_ini, !-> out
     & i_rho_ini_LHFW_found)                          !-> out
c-----finds the small radius rho_ini at the vector rho^ where
c     the wave specified by ioxm (LH or  FW) has the cutoff.
c     The vector rho^ starting at the magnetic axis O(rma,zma).
c     The vector points at x,y,z.
c
c     If (this subroutine could find rho_ini) then
c        it will set i_rho_ini_LHFW_found=1
c     else
c        it will set i_rho_ini_LHWF_found=0

      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
c-----input
      real*8 x,y,z
      double precision cnpar !N_parallel
      integer i_n_poloidal ! the grill case (see inigrill) 
      double precision n_theta_pol, !N  poloidal for i_n_poloidal=4 case
     &n_toroidal                    !N  toroidal for i_n_poloidal=4 case
c-----output
      double precision  rho_ini !normalized small radius in cutoff point
      double precision x_ini,y_ini,z_ini
      double precision cnx_ini,cny_ini,cnz_ini !N at rho_ini point
      integer i_rho_ini_LHFW_found
c-----external 
      double precision psi_rho,bxyz,dense_xyz
      double precision rho_loc
c-----local
      real*8 xl,yl,zl,r

      write(*,*)'grill_lh.f rho_ini_LHFW i_n_poloidal,n_theta_pol',
     &i_n_poloidal,n_theta_pol

      write(*,*)'rho_ini_LHFW_xyz'

      hstep= rho_step_find_LHFW_cutoff 

      pi=4.d0*datan(1.d0)
       
cyup      rho_loc= rho_initial_find_LHFW_cutoff     !initialization
      xl=x
      yl=y
      zl=z
      den=dense_xyz(xl,yl,zl,1) !-> get rho
      rho_loc=rho

 10   continue  ! handle for iterations

      r= dsqrt(xl*xl+yl*yl)
cSm061107
cyup      psi=psi_rho(rho_loc)
      psi= psif_xyz(xl,yl,zl) ! YuP New version ! spline or model
      write(*,*)'grill_lh in rho_ini_LHFW rho_loc= ',rho_loc
cSm050408   
      if(i_n_poloidal.eq.1) then
        cnper_tang= 0.d0 ! Nperp is purely along grad(psi) in this case
        call cninit12_xyz(xl,yl,zl,cnpar,cnper_tang,
     +                    cnx,cny,cnz,iraystop)
      else    
        bmod= bxyz(xl,yl,zl) !-> get b and grad(psi)
        gradpsi=dsqrt(dpdxd*dpdxd+dpdyd*dpdyd+dpdzd*dpdzd)
        cs=xl/r
        sn=yl/r
        ! Bpoloid = (b.[grad(psi) X e_phi]) / |gradpsi| 
        ! e_phi= {-y/r      ,   x/r     ,     0    }
        ! poloidal magnetic field
        b_theta=(bz*(dpdxd*cs+dpdyd*sn)-(bx*cs+by*sn)*dpdzd)/gradpsi
        bphi= -bx*sn + by*cs  
        cnper_tang= (n_theta_pol*bphi - n_toroidal*b_theta)/bmod
        if(i_n_poloidal.eq.2) then !find Nper_tang from Npar and Ntheta
           if(bphi .ne. 0.d0) then
             cnper_tang= (n_theta_pol*bmod - cnpar*b_theta)/bphi
           else ! bphi=0 ->  In this case Npar== +/- N_theta_pol 
             ! Note: n_toroidal is not known for this case. Set to 0?
             cnper_tang= -n_toroidal*b_theta/bmod !b_theta/b gives +or-
             ! stop 'rho_ini_LHFW_xyz: bphi=0; cnper_tang cannot be found'
           endif
        endif  !i_n_poloidal.eq.2
        call cninit_xyz(xl,yl,zl,cnpar,cnper_tang,
     &              cnx,cny,cnz,iraystop) !-> Nx,Ny,Nz found
      endif

c---------------------------------------------------------------
c     create the plot:
c     cold plasma dispersion function on N_perp
c     at the initial ray point
c---------------------------------------------------------------
      if(iraystop.eq.1) then
c-------the given wave does not exist in this point
        xl= xl*(1.d0-hstep)
        yl= yl*(1.d0-hstep)
        zl= zl*(1.d0-hstep)
        den=dense_xyz(xl,yl,zl,1) !-> get new rho
        rho_loc=rho
        if (rho_loc .lt. 0.1d0) then 
           write(*,*)'rho_ini_LHFW did not find the rho point'
           i_rho_ini_LHFW_found=0
           return
c           stop
        endif
        go to 10                   
      else
        cnper2= cnx**2+cny**2+cnz**2 - cnpar**2
        if (i_n_poloidal.eq.2) then 
          ! check the condition
          p= cnper2 - cnper_tang**2 
          if (p.lt.0.d0) then
             xl= xl*(1.d0-hstep)
             yl= yl*(1.d0-hstep)
             zl= zl*(1.d0-hstep)
             den=dense_xyz(xl,yl,zl,1) !-> get new rho
             rho_loc=rho
             goto 10
          endif
        endif
        i_rho_ini_LHFW_found=1        
        write(*,*)'in sub rho_ini_LHFW i_rho_ini_LHFW_found=',
     &            i_rho_ini_LHFW_found
        rho_ini=rho_loc
        x_ini=xl
        y_ini=yl
        z_ini=zl
        cnx_ini=cnx
        cny_ini=cny
        cnz_ini=cnz
      endif

      return
      end


c======================================================================
c======================================================================


      subroutine rho_ini_hot_nperp_roots_xyz(x_edge,y_edge,z_edge,
     + cnpar)     
c-----finds the small radius rho_ini > rho_min_find_hot_nperp_roots
c     at the vector rho^ where
c     hot plasma dispersdion function D_hot(npar) has three roots.
c     The vector rho^ is starting at the edge point (r_edge,z_edge,phi_edge),
c     and directed to the magnetic axis O(rma,zma,phi_edge)
c      
      implicit none
c      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
c-----input
      double precision
     &x_edge,y_edge,z_edge,          ! edge point coordinates
     &cnpar                            ! N_parallel
     
c-----output
      double precision  rho_ini(n_hot_roots_a), ! normalized small radius of M point 
                                                ! where D_hot(Nperp)=0 has 'i' roots
     &x_ini(n_hot_roots_a),y_ini(n_hot_roots_a),z_ini(n_hot_roots_a),
     +r_ini(n_hot_roots_a),
      !space coordinates of M point
 
     &N_perp_roots_diff_points(n_hot_roots_a,n_hot_roots_a) 
      ! hot plasma roots in n_hot_roots_a points
      
      double complex e_xyz_diff_points(3,n_hot_roots_a,n_hot_roots_a) 
                     ! electric filed Ex,Ey,Ez
                     ! for each root at i=1,... n_hot_roots_a points
                     ! e_xyz_diff_points(1,k,i) = E_x compomnent
                     ! e_xyz_diff_points(2,k,i) = E_y compomnent
                     ! e_xyz_diff_points(3,k,i) = E_z compomnent

      double precision e_mode_diff_points(3,n_hot_roots_a,n_hot_roots_a)  
                      !electric field polarization E+,E-,E_parallel
                      !at each N_perp_root(k)  an i=1,... n_hot_roots_a points
                      !e_mode_diff_points(1,k,i) = E+
                      !e_mode_diff_points(2,k,i) = E-
                      !e_mode_diff_points(3,k,i) = E||

c-----external
      double precision psi_rho,bxyz,wpw_2,wcw,tempe_xyz,
     + tpoprho,vflowrho,psif_xyz

c-----local
      double precision x_ar(nbulka),y_ar(nbulka),T_ar_ev(nbulka),     
     &tpop_ar(nbulka),vflow_ar(nbulka)![sm/sec]
      double precision psi,rho_loc,zaxis,raxis,theta,
     &costheta,sintheta,step,x,y,z,r,rho_edge, r_edge, xl,yl,zl

      integer i,j,k,i_found_n_hot_roots(n_hot_roots_a),n_hot_roots
 
      double precision N_perp_root_ar(n_hot_roots_a)    !hot plasma roots
      
      double complex e_xyz_nperp_root_ar(3,n_hot_roots_a) ! electric filed 
                                                          ! for each root
                               ! e_nperp_root_ar(1,k) E_x compomnent
                               ! e_nperp_root_ar(2,k) E_y compomnent
                               ! e_nperp_root_ar(3,k) E_z compomnent

      double precision  e_mode_root_ar(3,n_hot_roots_a)  !electric field
                                                         !polarization 
                                                         ! E+,E-,E_parallel
                                                         !at each N_perp_root(i)]
                               !e_mode_root_ar(1,k)      ! E+
                               !e_mode_root_ar(2,k)      ! E-
                               !e_mode_root_ar(3,k)      ! E||

      character*60  text
      character*6 format

      pi=4.d0*datan(1.d0)
      
      do i=1,n_hot_roots_a 
        i_found_n_hot_roots(i)=0
        do k=1,n_hot_roots_a 
          N_perp_roots_diff_points(k,i)=0.d0
        enddo
      enddo 

      step=rho_step_find_hot_nperp_roots 

      xl=x_edge
      yl=y_edge
      zl=z_edge
      r_edge = dsqrt(x_edge**2 + y_edge**2)

 10   continue


      bmod= bxyz(xl,yl,zl) ! for wcw
      do j=1,nbulk
        x_ar(j)=wpw_2(xl,yl,zl,j) ! finds rho  for density
        y_ar(j)=wcw(xl,yl,zl,j)
        if(j.eq.1) y_ar(1)=-y_ar(1) ! negative Y=(omega_ce/omega)
                                    ! for electrons
        T_ar_ev(j)=tempe_xyz(xl,yl,zl,j)*1.d+3 !(eV) averaged temperature
        tpop_ar(j)=tpoprho(rho,j)
        vflow_ar(j)=vflowrho(rho,j)         
      enddo

      call hot_roots_solver(nbulk,T_ar_ev,tpop_ar,vflow_ar,
     &x_ar,y_ar,cnpar,
     &cN_perp_root_max,n_points_root,  
     &N_perp_root_ar,n_hot_roots,
     &e_xyz_nperp_root_ar,e_mode_root_ar)

      write(*,*)' rho_ini_hot_nperp_roots after hot_roots_solver'
      write(*,*)'n_hot_roots',n_hot_roots
      write(*,*)'n_hot_roots_a',n_hot_roots_a

      do i=1,n_hot_roots_a

         write(*,*)'i,i_found_n_hot_roots(i)',i,i_found_n_hot_roots(i)

         if ((n_hot_roots.eq.i).and.(i_found_n_hot_roots(i).eq.0))then
c---------------------------------------------------------
c          coordinates of the first point along the line
c          where the solution D_hot(N_perp)=0 had 'i' roots
c-------------------------------------------------------           
           x_ini(i)=xl
           y_ini(i)=yl
           z_ini(i)=zl
           r= dsqrt(xl**2+yl**2)
           r_ini(i)= r 
           rho_ini(i)=rho ! found above in wpw_2
           i_found_n_hot_roots(i)=1
           do k=1,n_hot_roots_a 
             N_perp_roots_diff_points(k,i)=N_perp_root_ar(k)

             do j=1,3
                e_xyz_diff_points(j,k,i)= e_xyz_nperp_root_ar(j,k)
                e_mode_diff_points(j,k,i)=e_mode_root_ar(j,k)
             enddo
 
           enddo
         endif
      enddo

      if  (i_found_n_hot_roots(3).eq.1) goto 20 ! point with three roots
                                              ! was found 
      xl= xl*(1.d0-step)                                              
      yl= yl*(1.d0-step)                                              
      zl= zl*(1.d0-step)                                              
cyup      rho_loc= rho_loc*(1.d0-step)
      rho_loc=rho ! found above in wpw_2
        
      if (rho_loc.lt.rho_min_find_hot_nperp_roots) then 
         write(*,*)'forest.f rho_ini_hot_nperp_roots'
         write(*,*)'did not find 3 roots at interval'
         write(*,*)'(rho_min_find_hot_nperp_roots,rho_edge)'
         goto 20 
c         stop
      endif

      go to 10                   
    
 20   continue

      write(*,*)'open.f before open10' 
      open(10,file='find_hot_roots.dat')
      format='d21.15'
      text='(1X,"r_edge=",'//format//')'
      write(10,fmt=text)r_edge
      text='(1X,"z_edge=",'//format//')'
      write(10,fmt=text)z_edge
      text='(1X,"cnpar=",'//format//')'
      write(10,fmt=text)cnpar

      do i=1,n_hot_roots_a !loop over space points along rho

         write(*,*)'i,i_found_n_hot_roots(i)',i,i_found_n_hot_roots(i)

         if ((i_found_n_hot_roots(i).eq.1).and.(i.eq.1))then
            write(*,*)'i_found_n_hot_roots(i).eq.1)'
          write(*,*)'i,i_found_n_hot_roots(i)',i,i_found_n_hot_roots(i)

       text='(1X,"one root was found along small radius in point with")'
            write(10,fmt=text)
            text='(1X,"rho=",'//format//')'
            write(10,fmt=text)rho_ini(i)
            text='(1X,"r=",'//format//')'
            write(10,fmt=text)r_ini(i)
            text='(1X,"z=",'//format//')'
            write(10,fmt=text)z_ini(i)
c-------------------------------------------------------------------------------            
            text='(1X,"root N_perp=",'//format//')'
            write(10,fmt=text) N_perp_roots_diff_points(1,i)
            text='(1X,"ReEx=",'//format//')'
            write(10,fmt=text) dreal(e_xyz_diff_points(1,1,i))
            text='(1X,"ImEx=",'//format//')'
            write(10,fmt=text) dimag(e_xyz_diff_points(1,1,i))
            text='(1X,"ReEy=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(2,1,i))
            text='(1X,"ImEy=",'//format//')'
            write(10,fmt=text)dimag(e_xyz_diff_points(2,1,i))
            text='(1X,"ReEz=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(3,1,i))
            text='(1X,"ImEz=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(3,1,i))
            text='(1X,"   E+=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(1,1,i)
            text='(1X,"   E-=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(2,1,i)
            text='(1X,"   E||=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(3,1,i)
c-----------------------------------------------------------------------------

         endif

         if ((i_found_n_hot_roots(i).eq.1).and.(i.eq.2))then
           write(*,*)'i_found_n_hot_roots(i) i.eq.2)'
         write(*,*)'i,i_found_n_hot_roots(i)',i,i_found_n_hot_roots(i)

          text='(1X,"two roots was found along small radius in point")'
            write(10,fmt=text)
            text='(1X,"rho=",'//format//')'
            write(10,fmt=text)rho_ini(i)
            text='(1X,"r=",'//format//')'
            write(10,fmt=text)r_ini(i)
            text='(1X,"z=",'//format//')'
            write(10,fmt=text)z_ini(i)
c------------------------------------------------------------------------
            text='(1X,"root N_perp(1)=",'//format//')'
            write(10,fmt=text) N_perp_roots_diff_points(1,i)
            text='(1X,"ReEx=",'//format//')'
            write(10,fmt=text) dreal(e_xyz_diff_points(1,1,i))
            text='(1X,"ImEx=",'//format//')'
            write(10,fmt=text) dimag(e_xyz_diff_points(1,1,i))
            text='(1X,"ReEy=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(2,1,i))
            text='(1X,"ImEy=",'//format//')'
            write(10,fmt=text)dimag(e_xyz_diff_points(2,1,i))
            text='(1X,"ReEz=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(3,1,i))
            text='(1X,"ImEz=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(3,1,i))
            text='(1X,"   E+=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(1,1,i)
            text='(1X,"   E-=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(2,1,i)
            text='(1X,"   E||=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(3,1,i)
c------------------------------------------------------------------------
            text='(1X,"root N_perp(2)=",'//format//')'
            write(10,fmt=text) N_perp_roots_diff_points(2,i)
            text='(1X,"ReEx=",'//format//')'
            write(10,fmt=text) dreal(e_xyz_diff_points(1,2,i))
            text='(1X,"ImEx=",'//format//')'
            write(10,fmt=text) dimag(e_xyz_diff_points(1,2,i))
            text='(1X,"ReEy=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(2,2,i))
            text='(1X,"ImEy=",'//format//')'
            write(10,fmt=text)dimag(e_xyz_diff_points(2,2,i))
            text='(1X,"ReEz=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(3,2,i))
            text='(1X,"ImEz=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(3,2,i))
            text='(1X,"   E+=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(1,2,i)
            text='(1X,"   E-=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(2,2,i)
            text='(1X,"   E||=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(3,2,i)
c--------------------------------------------------------------------------
         endif

         if ((i_found_n_hot_roots(i).eq.1).and.(i.eq.3))then
          write(*,*)'i_found_n_hot_roots(i) i.eq.3)'
        write(*,*)'i,i_found_n_hot_roots(i)',i,i_found_n_hot_roots(i)

        text='(1X,"three roots was found along small radius in point")'
            write(10,fmt=text)
            text='(1X,"rho=",'//format//')'
            write(10,fmt=text)rho_ini(i)
            text='(1X,"r=",'//format//')'
            write(10,fmt=text)r_ini(i)
            text='(1X,"z=",'//format//')'
            write(10,fmt=text)z_ini(i)
c-------------------------------------------------------------------------
            text='(1X,"root N_perp(1)=",'//format//')'
            write(10,fmt=text) N_perp_roots_diff_points(1,i)
            text='(1X,"ReEx=",'//format//')'
            write(10,fmt=text) dreal(e_xyz_diff_points(1,1,i))
            text='(1X,"ImEx=",'//format//')'
            write(10,fmt=text) dimag(e_xyz_diff_points(1,1,i))
            text='(1X,"ReEy=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(2,1,i))
            text='(1X,"ImEy=",'//format//')'
            write(10,fmt=text)dimag(e_xyz_diff_points(2,1,i))
            text='(1X,"ReEz=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(3,1,i))
            text='(1X,"ImEz=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(3,1,i))
            text='(1X,"   E+=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(1,1,i)
            text='(1X,"   E-=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(2,1,i)
            text='(1X,"   E||=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(3,1,i)
c---------------------------------------------------------------------------
            text='(1X,"root N_perp(2)=",'//format//')'
            write(10,fmt=text) N_perp_roots_diff_points(2,i)
             text='(1X,"ReEx=",'//format//')'
            write(10,fmt=text) dreal(e_xyz_diff_points(1,2,i))
            text='(1X,"ImEx=",'//format//')'
            write(10,fmt=text) dimag(e_xyz_diff_points(1,2,i))
            text='(1X,"ReEy=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(2,2,i))
            text='(1X,"ImEy=",'//format//')'
            write(10,fmt=text)dimag(e_xyz_diff_points(2,2,i))
            text='(1X,"ReEz=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(3,2,i))
            text='(1X,"ImEz=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(3,2,i))
            text='(1X,"   E+=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(1,2,i)
            text='(1X,"   E-=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(2,2,i)
            text='(1X,"   E||=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(3,2,i)
c---------------------------------------------------------------------------
            text='(1X,"root N_perp(3)=",'//format//')'
            write(10,fmt=text) N_perp_roots_diff_points(3,i)
            text='(1X,"ReEx=",'//format//')'
            write(10,fmt=text) dreal(e_xyz_diff_points(1,3,i))
            text='(1X,"ImEx=",'//format//')'
            write(10,fmt=text) dimag(e_xyz_diff_points(1,3,i))
            text='(1X,"ReEy=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(2,3,i))
            text='(1X,"ImEy=",'//format//')'
            write(10,fmt=text)dimag(e_xyz_diff_points(2,3,i))
            text='(1X,"ReEz=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(3,3,i))
            text='(1X,"ImEz=",'//format//')'
            write(10,fmt=text)dreal(e_xyz_diff_points(3,3,i))
            text='(1X,"   E+=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(1,3,i)
            text='(1X,"   E-=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(2,3,i)
            text='(1X,"   E||=",'//format//')'
            write(10,fmt=text)e_mode_diff_points(3,3,i)
         endif
      enddo !i=1,n_hot_roots_a !
      close(10)
  
      return
      end



c======================================================================
c======================================================================

c======================================================================
c======================================================================





c        **********************dinit_1ray_xyz *****************
c        * this subroutine                                    *
c        * 1)if ray was launched outside the plasma           *
c        *   (istart=1, ECR wave)                             *
c        *   it detemines the point where the ray inter-      *
c        *   sects the plasma boundary(xu0,yu0,zu0),and the   *
c        *   tangent components of the refractive index	      *
c        *   cnteta,cnphi}                                    *
c        *   \hat{theta}=grad(psi) x \hat{phi}/abs(grad(psi)).*
c        *   Here \hat{theta} is the poloidal unit vector.    *
c        *   Here x is the vector product.		      *
c        *   If istart.ne.1 it launch the ray                 *
c        *   in the given point inside the plasma .	      *
c        * 2)	Then it calculates the parallel to magnetic   *
c        *   field component of the refractive index  cnpar.  *
c        *      Then it determines the normal   to   magnetic *
c        *    surface component of refractive index cnrho. *
c        *   It is directed inside the plasma		      *
c        *      Then it determines the initial components   *
c        *   of the refractive index:			      *
c        *              cnx,cny,cnz	      *
c        *   Then it calculate and prints the value of        *
c        *   dispersion function d=eps                        *
c        *   The result is  the initial condition  for the    *
c        *   ray equations				      *
c        ******************************************************
c
c------------------------------------------------------------------
c								   !
c        input parameters					   !
c        xst,yst,zst                                            !
c	 EC wave:						   !
c        alfast- toroidal angle(radian) of the ray at antenna      !
c        betast- angle(radian) between the  horizontal             !
c                          plane and the ray at antenna            !
c        LH and FW waves:					   !
c        cnx,cny,cnz						   !
c        output parameters					   !
c        u(1)=xu0,u(2)=yu0,u(3)=zu0,u(4)=cnx,u(5)=cny,u(6)=cnz	   !
c        iraystop- index to stop the ray determination(=1)	   !
c       	   or make the ray determination (=0)		   !
c------------------------------------------------------------------
c        it uses the following functions and subroutines           !
c        ias1r,bxyz,wcw,wpw_2,gamma1_xyz,s_xyz,abc_xyz,hamilt_xyz
c------------------------------------------------------------------
      subroutine dinit_1ray_xyz(xst,yst,zst,alfast,betast,cnteta,cnphi,
     1                       u,iraystop) !-> out
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'three.i'
      include 'five.i'
      include 'grill.i'
      include 'write.i'
      include 'rkutta.i'

      dimension u(6),deru(6)

cfor test_to plot D(ReN_perp,ImN_perp)
      integer n_param
      parameter (n_param=3)
      character*(15) name_param(n_param) !names of the input parameters
      real param(n_param)               !values of the input parameters     
cendtest

c-----for hot plasma roots
      double precision 
     &N_perp_root_ar(n_hot_roots_a)     !hot plasma roots
c---------------------------------------------------------------
c------for  rho_ini_LHFW
      integer  i_rho_ini_LHWH_found

cyup      write(*,*)' in dinit_1ray istart=',istart

      iraystop=0
      pi=4*datan(1.d0)
      trnspi=pi/180.d0
c---------------------------------------------------------------
c     if the ray starting point is outside the plasma (ECR -case)
c     then:
c          determine:    1) where the vacuum ray
c                           intersects the plasma boundary
c                        2) cnpar,cnper_tang- parallel and 
c                           perp-to-b_tang-to-psi
c                           components of the vacuum
c                           refractive index
c---------------------------------------------------------------
      call set_output ! initialize output.i
                      ! for each new ray

      xu0=xst ! to initialize
      yu0=yst
      zu0=zst
            
c      if (i_ox.eq.1) then
      call set_oxb         ! initialize oxb.i for each new ray
c      endif

cyup      write(*,*)'dinit_1ray_xyz: xst,yst,zst=',xst,yst,zst
      
      if (istart.eq.1) then
c--------EC wave 
         bmod= bxyz(xst,yst,zst) 
         !-----------------------------------------------!
c         call plasmray_xyz(xst,yst,zst,alfast,betast,    !
c     1                     xu0,yu0,zu0,iraystop)         !
         !-----------------------------------------------!
	   !write(*,*)'after plasmray_xyz: x,y=',x,y
	   !pause
         if (iraystop.eq.1) then
	      return
         end if
c	 -------------------------------
c        shift the initial point inside the plasma from the boundary
         x=xu0
         y=yu0
	   z=zu0
         !-YuP: do not shift:    call edgcor_xyz(x,y,z, xu0,yu0,zu0) 
c        end of the shift
c	 -------------------------------
c        nx,ny,nz at starting point:
         cnxst=dcos(betast)*dcos(alfast+phist)
         cnyst=dcos(betast)*dsin(alfast+phist)
         cnzst=dsin(betast)
c----------------------------------------------------------------
c        Calculate the components of unit vector  which is
c        perpendicular to b-field and TANGENTIAL to psi-flux-surface:
         call unit_vectors_perp_b(xu0,yu0,zu0, 
     +                     uvtx, uvty, uvtz,  uvnx, uvny, uvnz) !->out
         ! Refractive index perp. to b and tang. to psi:
         cnper_tang= cnxst*uvtx + cnyst*uvty + cnzst*uvtz
         ! Refractive index parallel to b:
         cnpar1=    (cnxst*bx   + cnyst*by   + cnzst*bz)/bmod
      endif !(istart.eq.1)


      if ((istart.eq.2).or.(istart.eq.3)) then
c--------LH and FW wave, OX-conversion pt.
         ! Find cnper_tang from known (Nteta, Nphi):
         x=xst
         y=yst
         z=zst
         r=dsqrt(x*x+y*y)
         cs=x/r
         sn=y/r
         bmod= bxyz(x,y,z) !-> get components of b and grad(psi)
         grpsi= dsqrt(dpdxd**2 + dpdyd**2 + dpdzd**2) ! |grad(psi)| 
         ! poloidal and toroidal magnetic field
         bpol=(bz*(dpdxd*cs+dpdyd*sn)-(bx*cs+by*sn)*dpdzd)/grpsi
         bphi= -bx*sn + by*cs  
         ! Refractive index perp. to b and tang. to psi:
         cnper_tang= (cnteta*bphi-cnphi*bpol)*o_bmod ! needs work ?
         ! Refractive index parallel to b:
         cnpar1=(cnphi*bphi+cnteta*bpol)*o_bmod ! needs work ?
      end if
      write(*,'(a,3e11.3)')'dinit_1ray_xyz(entry): x,y,z=',x,y,z
c--------------------------------------------------------------
      cnpar2= cnpar1**2
c-------------------------------------------------------------------
      write(*,*)'dinit_1ray_xyz: cnpar1,cnper_tang=',cnpar1,cnper_tang
      write(*,*)'i_rho_find_hot_nperp_roots',i_rho_find_hot_nperp_roots
c---------------------------------------------------------------
      if (i_rho_find_hot_nperp_roots.eq.1) then
c-----------------------------------------------------------------
c       finds the small radius rho_ini > rho_min_find_hot_nperp_roots
c       at the vector rho^ where
c       hot plasma dispersdion function D_hot(npar) has three roots.
c       The vector rho^ is starting at the edge point (x_edge,y_edge,z_edge),
c       and directed to the magnetic axis 
c-------------------------------------------------------------------
        call rho_ini_hot_nperp_roots_xyz(x,y,z,cnpar1)     
        stop 'dinit_1ray_xyz after call rho_ini_hot_nperp_roots'
      endif
c------------------------------------------------------------------    
c-------------------------------------------------------------
c     fit the initial value of rho_ini for LH or FW cutoff
c--------------------------------------------------------------
      write(*,*)'dinit_1ray_xyz i_rho_cutoff ',i_rho_cutoff
      write(*,*)'dinit_1ray_xyz before i_rho_cutoff=1 z,r',z,r
      if (i_rho_cutoff.eq.1) then
         call rho_ini_LHFW_xyz(x,y,z,cnpar1,
     &   i_n_poloidal, n_theta_pol, cnphi,
     &   rho_ini, x_ini,y_ini,z_ini, cnx_ini,cny_ini,cnz_ini,
     &   i_rho_ini_LHFW_found)

         write(*,*)
     +    'dinit_1ray_xyz after rho_ini_LHFW i_rho_ini_LHFW_found=',
     &              i_rho_ini_LHFW_found

         if (i_rho_ini_LHFW_found.eq.1) then
c-----------cutoff point with new x,y,z coordinates was found
            x=x_ini 
            y=y_ini 
            z=z_ini 
            r= dsqrt(x*x+y*y)
            !write(*,*)'dinit_1ray_xyz after  rho_ini_LHFW z,r',z,r
         else
c-----------cutoff point was not found
            write(*,*)'dinit_1ray_xyz: cutoff point was not found'
            write(*,*)'dinit_1ray_xyz:  iraystop->1'
            iraystop=1
            return
         endif
      endif

      if (i_rho_cutoff.eq.0) then
         cnper_tang_ini= cnper_tang
      endif

c--------------------------------------------------------------------
c     cninit solves the dispersion relation N=N(n_par)
c     Then subroutine calculates the initial components
c     of the refractive index  cnx,cny,cnz
c---------------------------------------------------------
cBH070123 start
      
cyup      bmod=bxyz(x,y,z) ! find rho

cSAP090518
      x_e=wpw_2(x,y,z,1)
      y_e=wcw(x,y,z,1)
      !write(*,*)'x_e,y_e',x_e,y_e

      if(nbulk.ge.2) then
        x_2=wpw_2(x,y,z,2)
        y_2=wcw(x,y,z,2)
        write(*,*)'x_2,y_2',x_2,y_2
      endif

      if (nbulk.ge.3)then
        x_3=wpw_2(x,y,z,3) 
        y_3=wcw(x,y,z,3)
        write(*,*)'x_3,y_3',x_3,y_3
      endif

      if(nbulk.eq.2) then
         w_lh_d_w=dsqrt(x_2)
         write(*,*)'w_lh_d_w',w_lh_d_w
       endif
       if(nbulk.eq.3) then
         w_lh_d_w=dsqrt(x_2+x_3)
         write(*,*)'w_lh_d_w',w_lh_d_w
       endif
cBH070123 end

      write(*,'(a,3i3)')'dinit_1ray_xyz: original id,ioxm_n_npar,ioxm=',
     + id, ioxm_n_npar, ioxm
      
      ioxm_n_npar_loc=ioxm_n_npar ! save the original value
      id_loc=id  ! save the original value
      ioxm_loc=ioxm ! save the original
      
      if(i_ox.eq.1) id=2
      
      write(*,*)
      write(*,*)'Dinit_1ray_xyz-> nper_npar_ioxm_n_npar_xyz'
      write(*,*)'-------------------------------------------------'
      write(*,*)'xe==(omega_pe/omega)^2 = ',x_e
      write(*,*)'ye==(omega_ce/omega)   = ',y_e
      write(*,*)'For Whistler need: Npar^2= 1- xe/(1-ye)=',1-x_e/(1-y_e)
      write(*,*)'For Whistler launch, make sure xe>1, ye>1, '
      write(*,*)'and Npar = sqrt[1- xe/(1-ye)] or a little larger.'
      write(*,*)'Npar^2 from input data  ==', cnpar2
      write(*,*)'-------------------------------------------------'
      
      cnper=0.d0 ! to initialize
      
      ioxm_n_npar=1 ! Try +1 
      iraystop=0
      write(*,*)
      write(*,*)'Dinit calling nper_npar_ioxm_n_npar_xyz ioxm_n_npar=+1'
      call nper_npar_ioxm_n_npar_xyz(2,x,y,z,cnpar1,
     & cnper_p,iraystop) ! ioxm_n_npar is in one.i    !-> Nper(Npar)
      write(*,*)'iraystop,ioxm_n_npar,Nper_p=',iraystop,1,cnper_p

      if(iraystop.eq.0) then ! root found
         wn_perp_ioxm_p(ifreq_write)=cnper_p
         cnper= cnper_p
      else ! no root
         wn_perp_ioxm_p(ifreq_write)=-1.d0 
      endif

      ioxm_n_npar=-1 ! Try -1
      iraystop=0
      write(*,*)
      write(*,*)'Dinit calling nper_npar_ioxm_n_npar_xyz ioxm_n_npar=-1'
      call nper_npar_ioxm_n_npar_xyz(2,x,y,z,cnpar1,
     & cnper_m,iraystop) ! ioxm_n_npar is in one.i    !-> Nper(Npar)
      write(*,*)'iraystop,ioxm_n_npar,Nper_m',iraystop,-1,cnper_m

      if(iraystop.eq.0) then ! root found
         wn_perp_ioxm_m(ifreq_write)=cnper_m
         cnper= cnper_m
         if(wn_perp_ioxm_p(ifreq_write).gt.0.d0) then
            write(*,*)'Dinit after nper_npar_ioxm_n_npar_xyz:'
            write(*,*)'Two roots found:', cnper_m,cnper_p
            if(ioxm.eq. 1) then ! O-mode
               cnper= cnper_p !min(cnper_m,cnper_p)
               write(*,*)'ioxm=+1: Select O-mode.  Nper=',cnper
            endif
            if(ioxm.eq.-1) then ! X-mode
               cnper= cnper_m !max(cnper_m,cnper_p)
               write(*,*)'ioxm=-1: Select X-mode.  Nper=',cnper
            endif
         endif
      else ! no root
         wn_perp_ioxm_m(ifreq_write)=-1.d0
      endif
            
      ioxm_n_npar= ioxm_n_npar_loc ! restore the original value

c      if (i_look_roots.eq.1)then   
c-----------------------------------------------------------------
c       plot ReD_hot(nperp) at given npar
c-----------------------------------------------------------------
c        name_param(1)='xe'
c        name_param(2)='ye'
c        name_param(3)='npar'
c        param(1)=xe
c        param(2)=ye
c        param(3)=cnpar1
c        write(*,*)'dinit.f before map_d_hot n_param', n_param
c        write(*,*)'name_param ',name_param 
c        write(*,*)'param ', param
c        call map_d_cold(z,r,phi,cnpar1,n_nperp_plot,cnperp_plot_min,
c     &  cnperp_plot_max,
c     &  name_param,param,n_param)   
c        call map_d_hot(z,r,phi,cnpar1,n_nperp_plot,cnperp_plot_min,
c     &  cnperp_plot_max,
c     &  name_param,param,n_param)
c        call map_d_hot(z,r,phi,cnpar1,n_nperp_plot,cnperp_plot_min,
c     &  cnperp_plot_max,
c     &  name_param,param,n_param)
c        write(*,*)'dinit.f after map_d_hot'
c        call rho_ini_hot_nperp_roots(r,z,phi,cnpar1)  
c--------------------------------------------------------------
c       calculate all roots N_perpendicular of the hot plasma
c       dispersion function D_hot(N_perp=0) at the interval
c       0 < N_perpendicular < cN_perp_root_max
c----------------------------------------------------------
c        call calculate_hot_nperp_roots(z,r,phi,cnpar1,
c     &  n_hot_roots,N_perp_root_ar)
c        iraystop=1
c        return
c      endif !   i_look_roots=1    
c-----------------------------------------------------------------
      if ((id.eq.1).or.(id.eq.2)) then
c-------calculate two cold plasma roots for different ioxm
        ioxm_loc=ioxm ! save the original ! Careful! hamilt_xyz() depends on ioxm !!!
        !--------------------------
        !1. Try the original ioxm
        iraystop=0
        write(*,*)
        write(*,*)'dinit_1ray-> cninit_xyz trying ioxm====', ioxm
        call cninit_xyz(x,y,z,cnpar1,cnper_tang_ini,
     1              cnx,cny,cnz,iraystop) !-> out
        cnper=-1.d0 ! to initialize
        cnper_p=-1.d0 ! to initialize
        if (iraystop.eq.0) then
          cnper=dsqrt(cnx**2+cny**2+cnz**2 -cnpar1**2)   
          cnper_p=cnper
          ioxm_p=ioxm ! save    ! Careful! hamilt_xyz() depends on ioxm !!!
          wn_perp_ioxm_p(ifreq_write)=cnper_p
        else
          write(*,*)
     +      'dinit_1ray:WARNING COULD NOT FIND ROOT FOR ORIGINAL ioxm'
          !2. Try the opposite ioxm
          ioxm=-ioxm 
          iraystop=0
          write(*,*)
          write(*,*)'dinit_1ray-> cninit_xyz trying ioxm====', ioxm
          call cninit_xyz(x,y,z,cnpar1,cnper_tang_ini,
     1              cnx,cny,cnz,iraystop) !-> out
          cnper_m=-1.d0 ! to initialize
          if (iraystop.eq.0) then
            cnper=dsqrt(cnx**2+cny**2+cnz**2 -cnpar1**2)   
            cnper_m=cnper
            ioxm_m=ioxm ! save    ! Careful! hamilt_xyz() depends on ioxm !!!
            wn_perp_ioxm_m(ifreq_write)=cnper_m
            ioxm_loc=ioxm !OVER-WRITE THE ORIGINAL
            write(*,*)
     +       'dinit_1ray:WARNING FOUND ROOT FOR opposite ioxm=',ioxm
          else
            write(*,*)
     +      'dinit_1ray:WARNING COULD NOT FIND ROOT FOR opposite ioxm=',
     +        ioxm
            stop
            !if(cnper_m.gt.0.d0) then
            !  ! if two roots, cnper_p and cnper_m were found, 
            !  ! choose the smaller Nperp (for whistler mode)
            !  cnper=min(cnper_p,cnper_m)
            !  if(cnper.eq.cnper_p) ioxm_loc=ioxm_p ! save    
            !  if(cnper.eq.cnper_m) ioxm_loc=ioxm_m ! save    
            !  ! Careful! hamilt_xyz() depends on ioxm !!!
            !endif
          endif
        endif
        wye_0(ifreq_write)=wcw(x,y,z,1)
        wxe_0(ifreq_write)=wpw_2(x,y,z,1)
      endif !id=1,2

c Run 3rd time, (for id=1 or 2) with proper ioxm:
      write(*,'(a,3i3)')'dinit_1ray RERUN for  id, ioxm_n_npar, ioxm=',
     +  id, ioxm_n_npar, ioxm
      call cninit_xyz(x,y,z,cnpar1,cnper_tang_ini,
     1                cnx,cny,cnz,iraystop) !-> out

      id=id_loc ! restore to original value
      ioxm=ioxm_loc
      ioxm_n_npar=ioxm_n_npar_loc
      write(*,'(a,3i3)')'dinit_1ray_xyz: restore id, ioxm_n_npar,ioxm=',
     +  id, ioxm_n_npar, ioxm
      
      if (iraystop.eq.1) then
         return
      end if

      cn2=   cnx**2+cny**2+cnz**2
      cnper= dsqrt(cn2-cnpar1**2)

      bmod=bxyz(x,y,z) !-> get b and derivs of b
      gam= gamma1_xyz(x,y,z,cnx,cny,cnz)
      cnt2= cnx**2 + cny**2 + cnz**2
      dh=hamilt_xyz(x,y,z,cnt2)
c-------------------------------------------------------------------
c      write(*,*)'before outinit'
      call outinit_xyz(u)
c      write(*,*)'in dinit_1ray after call outini nrayelt= ',nrayelt
      irefl=0
      u(1)=x
      u(2)=y
      u(3)=z
      u(4)=cnx
      u(5)=cny
      u(6)=cnz
c      write(*,*)'dinit_1ray before prep3d powini',powini
c      write(*,*)'dinit_1ray before prep3d u',u

      call prep3d_xyz(u,deru,iraystop)

      write(*,*)'-----------------------------------------------------'
      write(*,*)'end dinit_1ray_xyz         INITIAL DATA for this ray:'
      write(*,'(a,6e11.3)')' x,y,z; Vx,Vy,Vz=', x,y,z, deru(1:3)
      write(*,'(a,3e11.3)')' Nx,Ny,Nz=',cnx,cny,cnz
      write(*,'(a,3e11.3)')' N^2, Nper, Npar =',cn2,cnper,cnpar1
      write(*,*)'rho=',rho
      write(*,*)'(wpe/w)^2=',wpw_2(x,y,z,1),'wce/w=',wcw(x,y,z,1)
      write(*,*)'-----------------------------------------------------'
      return
      end


c======================================================================
c======================================================================

c                               grill_lh_xyz
c     This subroutine calculates the initial conditions for      
c     the grill type wave launch
c     ***************************************************************
c       input parameters:				               
c
c       igrillpw: option specifying the N_parallel power spectra       
c                =1 > power=powers/nnkpar	(default)
c                =2 > power=sin**2(x)/x**2,			      
c                     x=2pi(n_par-0.5(anmax-anmin))/(anmax-anmin)
c                =3 > power=exp-((n_par-a1)/a2)**2,
c                           and normalized to powers().
c                           a1=anmin(1:ngrill), a2=anmax(1:ngrill)     
c       ngrill   is a number of the poloidal grill angles
c       ngrilla    is a maximal number of the poloidal grill angles
c       thgrill(ngrilla)    poloidal angle of grill,measured counter   
c                          clockwise from horizontal through the        
c                          magnetic axis (degrees)		       
c       height(ngrilla)     is a poloidal length (m) of grill 	       
c                          (giving poloidal power distribution	       
c                          of each grill).    			       
c       nthin(ngrilla)      is a number of ray near the each poloidal
c                          center, simulating a grill                  
c       phigrill(1;ngrilla) is a toroidal grill angle of grill (degrees)
c                n_parallel (ngrill=nspect)
c
c       igrilltw: option specifying poloidal distribution of power
c                =1 > spread with equal weight over height 
c                =2 > {cos(pi(theta_pol-thgrill(i))/(height/radius))}**2
c                     (default).
c                     -0.5height/radius<theta_pol-thgrill(i)<0.5height/radius
c
c       i_grill_pol_mesh: option specifying the poloidal mesh wtheta(j)
c                         near the central grill angle thgrill(i)
c                         =1 equispaced mesh 
c                            wtheta(j)-wtheta(j-1)=zdth=Const (default)
c                         =2 poloidal mesh will be chosen to get the equal
c			     power fpwth(j) for all rays near the central 
c                            grill angle fpwth(j)=1/nthini
c
c       i_grill_npar_ntor_npol_mesh: option specifying the refractive
c                         index meshes.
c
c                         For  i_n_poloidal=1,2,3 it is specifying
c                         n_parallel mesh anzin(n) for the power
c                         spectrum pwcpl(n) n=1,...,nnkpari
c                         =1 equispaced mesh 
c                            anzin(n)-anzin(n-1)=hnpar=Const (default)
c                         =2 n_parallel mesh will be chosen to get the equal
c			     power pwcpl(n) for all rays in the given power
c                            spectrum  pwcpl(n)=1.d0/nnkpari 
c                            pwcpl(n)=power_spectrum(anzin(n))*
c                                     delta_npar_bin(n)= 1.d0/nnkpari
c
c                            For  i_n_poloidal=4 it is specifying two meshes:
c                            a) n_toroidal mesh anztorin(ntor) and             
c                            b) n_poloidal mesh anzpolin(npolmesh) 
c                            for the power spectrum
c                            pwcpl_tp(1:nnktori,1:nnkpoli)=pwcpl_t*pwcpl_t 
c                         =1 equispaced meshes (default)
c                            anztorin(ntor)- anztorin(ntor-1)=hntor=Const 
c                            anzpolin(npol)- anzpolin(npol-1)=hnpol=Const
c                         =2 the meshes anztorin(1:nntori) anzpolin(1:nnkpoli)
c                            will be chosen to get the equal
c			     power pwcpl_tp(ntor,npol) for all rays in 
c                            the given power spectrum 
c                            pwcpl_tp(ntor,npol)=1.d0/(nnktori*nnkpoli) 
c                         
c
c       anmin(1:ngrilla)    position of the left bound of	       
c                          power spectrum P(n_parallel) (Can be neg).
c                          Or, as specified for igrillpw=3 above.
c                          It needs for i_n_poloidal=1,2,3  
c       anmax(1:ngrilla)    position of the right bounds  	       
c                          of power spectrum P(n_parallel)	       
c                          Or, as specified for igrillpw=3 above.
c                          It needs for i_n_poloidal=1,2,3   
c       nnkpar(1:ngrilla)   number of points  of power spectrum	       
c                          P(n_parallel)		
c                          It needs for i_n_poloidal=1,2,3 	       
c       powers(1:ngrilla)   power in one grill (MWatts)	    	       
c                          (total power of grill(in MWatts) will be    
c                          powtotlh=sum{powers}			       
c       rhopsi0(1:ngrilla)  initial psi for wave front (0<rhopsi0<1)
c       rma (m)	           major radius of magnetic axis	       
c       zma (m)            vertical shift of magnetic axis
c       psimag	       
c       nnkprmax=max{i=1,ngrill}nnkpar(i)
c
c       i_n_poloidal gives the type of the grill launch	
c                    =1,2,3 N_parallel and other variables will be given
c                    =4     N_toroidal and N_poloidal will be given
c
c       The following variables work for i_n_poloidal=4
c
c       antormin(1:ngrilla)    position of the left bound of	       
c                          power spectrum P(n_parallel) (Can be neg).
c                          Or, as specified for igrillpw=3 above.  
c       antormax(1:ngrilla)    position of the right bounds  	       
c                          of power spectrum P(n_parallel)	       
c                          Or, as specified for igrillpw=3 above.  
c       nnktor(1:ngrilla)   number of points  of power spectrum	       
c                          P(n_parallel)	
c       anpolmin(1:ngrilla)    position of the left bound of	       
c                          power spectrum P(n_parallel) (Can be neg).
c                          Or, as specified for igrillpw=3 above.  
c       anpolmax(1:ngrilla)    position of the right bounds  	       
c                          of power spectrum P(n_parallel)	       
c                          Or, as specified for igrillpw=3 above.  
c       nnkpol(1:ngrilla)   number of points  of power spectrum	       
c                          P(n_parallel)	
c       nnktormax=max{i=1,ngrill}nnktor(i)
c       nnkpolmax=max{i=1,ngrill}nnkpol(i)
c
c       n_theta_pol   the poloidal refractive index
c                     used in i_n_poloidal=2 case 
c--------------------------------------------------------------------
c       output parameters:					       
c       nray               total number of the rays 		       
c                          nray=Sum{i=1,ngrill}Sum_{j=1,nthin(i)}      
c                          Sum{n=1,nnkpar(i)}nnkpar(n)    
c------------------------- YuP:             
c       arxu0(nray)   (m)  initial cartesian coords. of	                       
c       aryu0(ntay)   (m)                  ray_iray                    
c       arzu0(nray)   (m)             starting points		       
c       arnpar(nray)       initial parallel-to-b refractive index	       
c       arnper_tang(nray)  initial refr. index perp. to b, tang. to psi
c       arntheta(nray)     initial poloidal refractive index	       
c       arnphi(nray)       initial toroidal refractive index	       
c------------------------- YuP            
c       powinilh(nray)     initial power flowing in the iray wave channel
c		           normalized Sum_{i=1,nray}powinilh(i)=powtotlh
c                          (erg/c)				       
c       fpwth(nthinmax)    spectrum poloidal distributions             
c                          normalized 1=Sum{j=1,nthin}fpwth(j)
c       wtheta(nthinmax)   poloidal angle mesh wtheta(j) radian        
c                          near the grill poloidal angles thgrill(i)   
c                          j=1,nthin(i)
c       wtheta_edge(nthinmax+1): wtheta_edge(j) j=1,,nthin=1
c                          wtheta(j)=0.5*
c                          (wtheta_edge(j)+wtheta_edge(j+1)),j=1,nthin
c		       
c       anzin(nnkprmax)    are the points of mesh n_parallel	       
c                          P(n_parallel)  anzin(j), j=1,nnkpar(i)
c
c       anzin_edge(nnkprmax+1) : anzin_edge(j), j=1,nnkpar(i)+1
c                         anzin(j)=0.5(anzin_edge(j)+anzin_edge(j+1))
c          	          j=1,nnkpar(i)
c
c       pwcpl(nnkprmax)    power spectrum pwcpl(n) on   	       
c                          n_parallel, normalization		       
c                          Sum{n=1,nnkpar(i)}pwspl(n)=powers(i)        
c       wdnpar0(nray)	   initial values of the n_parallel width of   
c                          the rays
c  
c       pwcpl_tp(nnktormax,nnkpolmax)    power spectrum pwcpl on
c                          (N_toroidal,N_poloidal) mesh parallel,
c                          It is for i_n_poloidal=4
c       normalization		       
c       Sum{ntor=1,nnktor(i),npol=1,nnkpol(i)}pwspl_tp(ntpr,npol)
c               =powers(i)  
c
c       powtotlh
c----------------------------------------------------------------
c       for i_n_poloidal=4 case
c       anztorin(nnktormax)   are the points of mesh n_toroidal       
c                          P(n_toroidal,n_ploidal), j=1:nnktori
c       anztpolin(nnkpolmax)   are the points of mesh n_poloidal
c                              j=1:nnkpoli
c       anztorin_edge(nnktormax+1) N_toroidal mesh, j=1:nnktori+1 
c                anztorin(j)=0.5(anztorin_edge(j+1)+anztorin_edge(j))
c
c       anzpolin_edge(nnkpolmax+1) N_poloidal mesh, j=1:nnkpoli+1
c                anztpolin(j)=0.5(anzpolin_edge(j+1)+anzpolin_edge(j))
c
c       pwcpl_t_ar(nnktormax), power spectrum on n_toroidal mesh 
c       pwcpl_p_ar(nnkpolmax), power spectrum on n_poloidal mesh 
c       pwcpl_tp(nnktormax,nnkpolmax)    power spectrum on   	       
c                          n_toroidal,n_poloidal mesh, 
c       normalization
c       Sum{ntor=1,nnktor(i),n_pol=1,nnkpol(i)}(pwspl_tp(ntor,npol)
c                                 =powers(i)
c     	       
c       ilaunch=1 gives explicit r0launch,phi0launch,z0launch launch
c                 location (meters and radians). 
c
c-----------------------------------------------------------------------
c  uses: double precision functions r(psi,theta),z(psi,theta),psiff(rho)
c-----------------------------------------------------------------------

      subroutine grill_lh_xyz(rhopsi0,ngrill,thgrill,phigrill,
     1 height,nthin,
     1 anmin,anmax,nnkpar,powers,powtotlh,
     1 antormin,antormax,nnktor,anpolmin,anpolmax,nnkpol,
     1 n_theta_pol,
     1 rma,zma,psimag,
     1 fpwth,wtheta,wtheta_edge,anzin,anzin_edge,pwcpl,
     1 anztorin,anzpolin,pwcpl_tp,
     1 anztorin_edge,anzpolin_edge,pwcpl_t_ar,pwcpl_p_ar,
     1 nray, arxu0, aryu0, arzu0,  !!!! YuP
     1 arnpar, arnper_tang, arntheta, arnphi, !!!! YuP
     1 powinilh,  
     1 wdnpar0,igrillpw,igrilltw,i_n_poloidal,
     1 ilaunch,r0launch,phi0launch,z0launch,i_grill_pol_mesh,
     1 i_grill_npar_ntor_npol_mesh)
     
      implicit none
      include 'param.i'
      include 'one.i'
      
      integer ngrill,nray,iray,nthini,nnkpari,
     .igrillpw,igrilltw,i_n_poloidal,ilaunch,
     .i_grill_pol_mesh,i_grill_npar_ntor_npol_mesh
      double precision thgrill(ngrilla),height(ngrilla),
     1          phigrill(ngrilla),rhopsi0(ngrilla),
     1          anmin(ngrilla),anmax(ngrilla),powers(ngrilla),
     1          antormin(ngrilla),antormax(ngrilla),
     1          anpolmin(ngrilla),anpolmax(ngrilla),
     1          n_theta_pol
      double precision r0launch,phi0launch,z0launch
      integer   nnkpar(ngrilla),nthin(ngrilla),
     & nnktor(ngrilla),nnkpol(ngrilla)
      double precision arxu0(*),aryu0(*),arzu0(*),
     1       arnpar(*),arnper_tang(*),  arntheta(*),
     1       arnphi(*),  powinilh(*)
      double precision fpwth(nthinmax),wtheta(nthinmax),
     1       wtheta_edge(nthinmax+1)
      double precision anzin(nnkprmax),anzin_edge(nnkprmax+1),
     &pwcpl(nnkprmax)

      double precision anztorin(nnktormax),anzpolin(nnkpolmax),
     &pwcpl_tp(nnktormax,nnkpolmax),
     &anztorin_edge(nnktormax+1),anzpolin_edge(nnkpolmax+1),
     &pwcpl_t_ar(nnktormax),pwcpl_p_ar(nnkpolmax)

      double precision rma,zma,psimag,powtotlh,tet,psi0,z,r,
     .rhomag,thaper,zdth,delth,zcos,znorm,an0,difrnpar,andwdth,
     .hnpar,anorm,delnpar,x1,xx,phi0,theta0,x0,y0,z0,r0,
     1 cnpar,cnper_tang, cnteta,cnphi
     
      double precision psi_rho,powert
c-----for i_n_poloidal=4 case, input N_toroidal and N_poloidal
      double precision antor0,anpol0,difrntor,difrnpol,
     &andwdth_tor,andwdth_pol,hntor,hnpol,hnwidth,delntor,delnpol,
     &pwcpl_t,pwcpl_p,y1,yy,anpar
      integer nnktori,nnkpoli
      integer i,j,n,ntor,npol,k

      double precision wdnpar0(*)
      double precision f_pow_poloid_1,f_pow_poloid_0, bxyz
      external f_pow_poloid_1,f_pow_poloid_0, bxyz
      
      real*8 cs,sn,b_theta,grpsi,x,y

      double precision aperture
      common /aperture/ aperture
      double precision f_pow_npar_2,f_pow_npar_3
      external f_pow_npar_2,f_pow_npar_3
      
      if (ilaunch.eq.1) then
         write(*,*)
      write(*,*)'grill_lh_xyz: ilaunch=1. All rays start from one point'
         write(*,*)
      endif

      write(*,*)'in grill_lh ngrill,i_n_poloidal,ilaunch',
     &                       ngrill,i_n_poloidal,ilaunch

      if (ngrill.lt.1) then 
        write(*,*)'in ngrill_lh ngrill<1 but it should be .ge.1'
        write(*,*)'Please change ngrill in genray.in file' 
        stop
      endif

      if (ngrill.gt.ngrilla) then
        write(*,*)'in ngrill_lh ngrilla<ngrill but it should be .ge.'
        write(*,*)'please change ngrilla in param.i and recompile'
        stop
      endif
   
      do i=1,ngrill
         write(*,*)'i,thgrill(i),nthin(i)',i,thgrill(i),nthin(i)
         write(*,*)'height(i),phigrill(i),rhopsi0(i),powers(i)',
     1              height(i),phigrill(i),rhopsi0(i),powers(i)
         if ((i_n_poloidal.ge.1).and.(i_n_poloidal.le.3)) then
           write(*,*)'anmin(i),anmax(i),nnkpar(i)',
     1                anmin(i),anmax(i),nnkpar(i)
         else !i_n_poloidal=4
           write(*,*)'antormin(i),antormax(i),nnktor(i)',
     1               antormin(i),antormax(i),nnktor(i)
           write(*,*)'anpolmin(i),anpolmax(i),nnkpol(i)',
     1                anpolmin(i),anpolmax(i),nnkpol(i)
         endif
      enddo !i
c-----------------------------------------------------------
c     powtotlh is in MWatts
c     Transformation of the total antenna power from MW to erg/sec
c-----------------------------------------------------------------
      powtotlh=0.d0
      do i=1,ngrill
         powers(i)=powers(i)*1.d13   !erg/sec
         powtotlh=powtotlh+powers(i)
      enddo

      write(*,*)'powtotlh',powtotlh
c----------------------------------------------------------------
      pi=4.d0*datan(1.d0)
      tet=pi/180.d0
c-------------------------------------------------------------------
      do i=1,ngrill
         thgrill(i)=thgrill(i)*tet    ! in rad. now
         phigrill(i)=phigrill(i)*tet
      end do
c----------------------------------------------------------------
c     For each i, calculation of
c     wtheta (j)  -poloidal mesh, j=1,nthin(i)
c     fpwth  (j)  -poloidal power distributions normalized to
c                    1=Sum{j=1,nthin(i)}fpwth(j)
c----------------------------------------------------------------

      iray=0      !Ray counter

      write(*,*)'grill_lh iray=0, ngrill',ngrill

      do i=1,ngrill
         nthini=nthin(i)
         write(*,*)'nthini',nthini

         if(nthini.eq.1) then
            wtheta(1)=thgrill(i)
            fpwth(1)=1.d0
         else
c------------------------------------------------------------
c           calculations: r=r(psi0,thgrill(i)),
c                         z=z(psi0,thgrill(i))
c------------------------------------------------------------
            if(model_rho_dens.eq.0) then
               psi0=psi_rho(rhopsi0(i))
               write(*,*)'grill_lh_xyz->zr_psith  #1'
               call zr_psith(psi0,thgrill(i),z,r)
            elseif(model_rho_dens.eq.1 .or. model_rho_dens.eq.2 .or.
     +             model_rho_dens.eq.3 .or. model_rho_dens.eq.4) then
               call rho_theta_phi_xyz(rhopsi0(i),thgrill(i),phigrill(i),
     +                             x,y,z) !->out
               r=dsqrt(x*x+y*y)
            else
               stop 'grill_lh_xyz: only for model_rho_dens=0,1,2,3,4'
            endif
            
c------------------------------------------------------------
	    rhomag=dsqrt((r-rma)**2+(z-zma)**2)
	    thaper=height(i)/rhomag
cSm050311
            aperture=thaper !put the data to common /aperture/

cSm050309
            if((i_grill_pol_mesh.eq.1).or.(igrilltw.eq.1))then
c------------------------------------------------------------
c             equispaced poloidal mesh 
c             wtheta(j)-wtheta(j-1)=zdth=Const(default)
c------------------------------------------------------------
              zdth=thaper/dfloat(nthini)

              do j=1,nthini
	         delth=-.5d0*thaper+zdth*(j-0.5)
	         wtheta(j)=thgrill(i)+delth
	         zcos=dcos(pi*delth/thaper)
                
c ----------------------------------------------------
c                fpwth(j) is the poloidal power distribution on
c                       poloidal angle mesh wtheta (j)
c                       near the grill poloidal angle thgrill(i)
c----------------------------------------------------               
                 if (igrilltw.eq.1) then
                    fpwth(j)=1.0
                 else
                    fpwth(j)=zcos*zcos
                 endif
c                 write(*,*)'fpwth(j)',fpwth(j)
              enddo !j
            endif !i_grill_pol_mesh.eq.1
cSm050309
            if((i_grill_pol_mesh.eq.2).and.(igrilltw.ne.1))then
c-------------------------------------------------------------
c             poloidal mesh will be chosen to get the equal
c             power fpwth(j) for all rays near the central 
c             grill angle fpwth(j)=1/nthini
c-------------------------------------------------------------
        
c_old_method using   f_pow_poloid_1,        
c              call create_equal_mesh(-0.5*thaper,0.5*thaper,
c     &        f_pow_poloid_1,nthini,
c     &        wtheta_edge,wtheta,fpwth)

c              do k=1,nthini
c                wtheta(k)= wtheta(k)+thgrill(i)
c              enddo ! k

c              do k=1,nthini+1
c                wtheta_edge(k)= wtheta_edge(k)+thgrill(i)
c              enddo ! k

c              do k=1,nthini
c              write(*,*)'!1 k,wtheta(k),wtheta_edge(k),fpwth(k)',
c     &                    k,wtheta(k),wtheta_edge(k),fpwth(k)
c              enddo
c              write(*,*)'1 wtheta_edge(nthini+1)',wtheta_edge(nthini+1)

c_new method using f_pow_poloid_0
              call create_equal_mesh(-0.5*pi,0.5*pi,
     &        f_pow_poloid_0,nthini,
     &        wtheta_edge,wtheta,fpwth)

              do k=1,nthini
                wtheta(k)= wtheta(k)*aperture/pi+thgrill(i)
              enddo ! k
            
              do k=1,nthini+1
                wtheta_edge(k)= wtheta_edge(k)*aperture/pi+thgrill(i)
              enddo ! k

           endif !i_grill_pol_mesh.eq.2
	 endif ! (nthini.eq.1)
	 
c--------------------------------------------
c        normalization of fpwth
	 znorm=0.d0
	 do j=1,nthini
	    znorm=znorm+fpwth(j)
	 enddo !j

	 do j=1,nthini
	    fpwth(j)=fpwth(j)/znorm
	    write(*,*)'j,fpwth(j),wtheta(j)',j,fpwth(j),wtheta(j)
	 enddo !j
        
         if((i_n_poloidal.ge.1).and.(i_n_poloidal.le.3)) then
c-----------------------------------------------------------------------
c          calculation of anzin(n), the n_parallel mesh,
c          and the corresponding 
c          power spectrum pwcpl(n).
c      	   Normalization:  Sum{n=1,nnkpar(i)}pwspl(n)=powers(i)
c          BH020826:  Added igrillpw.eq.3 option, and slightly
c                     modified n_par grid of points.
c-----------------------------------------------------------------------
          
           if(igrillpw.eq.1 .or. igrillpw.eq.2) then
             an0=0.5d0*(anmin(i)+anmax(i)) !central value of N_parallel
	     difrnpar=anmax(i)-anmin(i)
           elseif (igrillpw.eq.3) then
             an0=anmin(i)
	     difrnpar=2.*anmax(i)  ! power out to exp(-4.)
           endif
            
	   nnkpari=nnkpar(i)
	   write(*,*)'nnkpari',nnkpari
	   if(nnkpari.eq.1)then
	     anzin(1)=an0
	     pwcpl(1)=powers(i)
c	     write(*,*)'nnkpari',nnkpari,'anzin(1),pwcpl(1)',
c     1                                   anzin(1),pwcpl(1)
	     goto 10
	   else
ccBH020826    difrnpar=anmax(i)-anmin(i)
	     andwdth=2.d0*pi/difrnpar
            
             if(i_grill_npar_ntor_npol_mesh.eq.1) then
c---------------------------------------------------------
c              equispaced mesh 
c              anzin(n)-anzin(n-1)=hnpar=Const (default)
c--------------------------------------------------------
cBH020826      hnpar=difrnpar/dfloat(nnkpari+1)
	       hnpar=difrnpar/dfloat(nnkpari)

	       do n=1,nnkpari
cBH020826        anzin(n)=anmin(i)+n*hnpar
                 anzin(n)=an0-0.5*difrnpar+(n-0.5)*hnpar
	       enddo

	       anorm=0.d0
	       do n=1,nnkpari
               delnpar=anzin(n)-an0
               if(igrillpw.eq.1) pwcpl(n)=1.d0/nnkpari ! YuP[08-2016] 
               !YuP: moved it outside of if(delnpar.ne.0.d0) condition (below)
               !For igrillpw.eq.1 All rays should have same power, 
               !no matter what delnpar is 
	         if(igrillpw.ne.1)then
	         if(delnpar.ne.0.d0) then
	           x1=andwdth*delnpar
	           xx=x1*x1
    	           if(igrillpw.eq.2) pwcpl(n)=dsin(x1)*dsin(x1)/xx
                 if(igrillpw.eq.3) pwcpl(n)=exp(-(delnpar/anmax(i))**2)
	         else
	           pwcpl(n)=1.d0
	         endif
	         endif
	         anorm=anorm+pwcpl(n)
	       enddo !n
c              write(*,*)'anorm='anorm ! should be 1.0 for igrillpw.eq.1

	       do n=1,nnkpari
	         pwcpl(n)= (pwcpl(n)/anorm)*powers(i)
             write(*,*)'grill_lh_xyz: POWERS in RAY#n: n,pwcpl(n)',
     +       n,pwcpl(n)
	       enddo !n
             endif !i_grill_npar_ntor_npol_mesh.eq.1
            
             if(i_grill_npar_ntor_npol_mesh.eq.2) then
c--------------------------------------------------------------------
c               n_parallel mesh will be chosen to get the equal
c		power pwcpl(n) for all rays in the given power
c               spectrum  pwcpl(n)=1.d0/nnkpari 
c               pwcpl(n)=power_spectrum(anzin(n))*
c                        delta_npar_bin(n)= 1.d0/nnkpari
c---------------------------------------------------------------------
                if (igrillpw.eq.2) then
                  andwdth=2.d0*pi/difrnpar
                  call create_equal_mesh(-pi,pi,f_pow_npar_2,nnkpari,
     &            anzin_edge,anzin,pwcpl)

                  do k=1,nnkpari 
                    anzin(k)=anzin(k)/andwdth+an0 
                    write(*,*)'k, anzin(k)',k, anzin(k)
                  enddo 

                  do k=1,nnkpari+1                   
                    anzin_edge(k)=anzin_edge(k)/andwdth+an0
                    write(*,*)'k, anzin_edge(k)',k,anzin_edge(k)
                  enddo
                endif !igrillpw.eq.2

                if(igrillpw.eq.3) then 
                  call create_equal_mesh(-1.d0,1.d0,f_pow_npar_3,
     &            nnkpari,anzin_edge,anzin,pwcpl)

                  do k=1,nnkpari 
                    anzin(k)=anzin(k)*0.5d0*difrnpar+an0 
                    write(*,*)'k,anzin(k)',k,anzin(k)
                  enddo
 
                  do k=1,nnkpari+1 
                    anzin_edge(k)=anzin_edge(k)*0.5d0*difrnpar+an0 
                    write(*,*)'k,anzin_edge(k)',k,anzin_edge(k)
                  enddo
                endif !igrillpw.eq.3

c------------------------------------------------------------------------
c               power normalization
                anorm=0.d0
	        do n=1,nnkpari
	           anorm=anorm+pwcpl(n)
	        enddo !n
c               write(*,*)'anorm='anorm

	        do n=1,nnkpari
	         pwcpl(n)=pwcpl(n)/anorm*powers(i)
	         write(*,*)'n,pwcpl(n)',n,pwcpl(n)
	        enddo !n
 
            endif !i_grill_npar_ntor_npol_mesh.eq.2

	   endif ! nnkpari.eq.1
10         continue
c--------------------------------------------------------------------
         endif !i_n_poloidal=1,2,3
         
         if(i_n_poloidal.eq.4) then
c-----------------------------------------------------------------------
c          calculation of anztorin(ntor), the N_toroidal mesh,
c          calculation of anzpolin(ntor), the N_poloidal mesh,
c
c          and the corresponding 
c          power spectrum pwcpl_tp(ntor,npol).
c      	   Normalization: 
c           Sum{ntor=1,nnktor(i),npol=1,nnkpol(i)}pwspl_tp(ntor,npol)
c                   =powers(i)
c          BH020826:  Added igrillpw.eq.3 option, and slightly
c                     modified n_par grid of points.
c-----------------------------------------------------------------------

           write(*,*)'i_n_poloidal.eq.4 igrillpw',igrillpw

           if((igrillpw.eq.1).or.(igrillpw.eq.2)) then
             antor0=0.5d0*(antormin(i)+antormax(i)) 
                                           !central value of N_toroidal
             anpol0=0.5d0*(anpolmin(i)+anpolmax(i)) 
                                           !central value of N_poloidal
	     difrntor=antormax(i)-antormin(i)
             difrnpol=anpolmax(i)-anpolmin(i)

             write(*,*)'i,antormax(i),antormin(i),antor0,difrntor',
     &                  i,antormax(i),antormin(i),antor0,difrntor
             write(*,*)'i,anpolmax(i),anpolmin(i),anpol0,difrnpol',
     &                  i,anpolmax(i),anpolmin(i),anpol0,difrnpol
             

           elseif (igrillpw.eq.3) then
             antor0=antormin(i)
             anpol0=anpolmin(i)
	     difrntor=2.*antormax(i)  ! power out to exp(-4.)
             difrnpol=2.*anpolmax(i)  ! power out to exp(-4.)
           endif

	   nnktori=nnktor(i)
           nnkpoli=nnkpol(i)       

	   write(*,*)'nnktori,nnkpoli',nnktori,nnkpoli
	   if((nnktori.eq.1).and.(nnkpoli.eq.1))then
	     anztorin(1)=antor0
             anzpolin(1)=anpol0
 	     pwcpl_tp(1,1)=powers(i)

	     write(*,*)'nnktori,nnkpoli',nnktori,nnkpoli

             write(*,*)'anztorin(1),anzpolin(1),pwcpl_tp(1,1)',
     &       anztorin(1),anzpolin(1),pwcpl_tp(1,1)

c            Save composite width for nominal npar_width, wdnpar:
             hnwidth=0.05*sqrt(antor0**2+anpol0**2)

	     goto 20
	   else

	     andwdth_tor=2.d0*pi/difrntor
             andwdth_pol=2.d0*pi/difrnpol
 
             if(i_grill_npar_ntor_npol_mesh.eq.1) then
c---------------------------------------------------------
c              equispaced meshs 
c              anztorin(n)-anztorin(n-1)=hntor=Const (default)
c              anzpolin(n)-anzpolin(n-1)=hnpol=Const (default)
c---------------------------------------------------------
	       hntor=difrntor/dfloat(nnktori)
               hnpol=difrnpol/dfloat(nnkpoli)
             
	       do ntor=1,nnktori
                 anztorin(ntor)=antor0-0.5*difrntor+(ntor-0.5)*hntor
	       enddo

	       do npol=1,nnkpoli
                 anzpolin(npol)=anpol0-0.5*difrnpol+(npol-0.5)*hnpol
	       enddo

c              Save composite width for nominal npar_width, wdnpar:
               hnwidth=sqrt(hntor**2+hnpol**2)

	       anorm=0.d0

	       do ntor=1,nnktori
                 delntor=anztorin(ntor)-antor0
	         x1=andwdth_tor*delntor
	         xx=x1*x1
                 do npol=1,nnkpoli
                   delnpol=anzpolin(npol)-anpol0
	           y1=andwdth_tor*delntor
	           yy=y1*y1
	           if(delntor.ne.0.d0) then
                     if(igrillpw.eq.1) pwcpl_t=1.d0/nnktori
    	             if(igrillpw.eq.2) pwcpl_t=dsin(x1)*dsin(x1)/xx
               if(igrillpw.eq.3) pwcpl_t=exp(-(delntor/antormax(i))**2)
	           else
	              pwcpl_t=1.d0
	           endif

	           if(delnpol.ne.0.d0) then
                     if(igrillpw.eq.1) pwcpl_p=1.d0/nnktori
    	             if(igrillpw.eq.2) pwcpl_p=dsin(y1)*dsin(y1)/yy
               if(igrillpw.eq.3) pwcpl_p=exp(-(delnpol/anpolmax(i))**2)
	           else
	             pwcpl_t=1.d0
	           endif

                   pwcpl_tp(ntor,npol)=pwcpl_t*pwcpl_t

	           anorm=anorm+pwcpl_tp(ntor,npol)
	         enddo !npol
               enddo ! ntor
c              write(*,*)'anorm='anorm

	       do ntor=1,nnktori
                 do npol=1,nnkpoli
	         pwcpl_tp(ntor,npol)=pwcpl_tp(ntor,npol)/anorm*powers(i)
	           write(*,*)'ntor,npol,pwcpl_tp(ntor,npol)',
     &                    ntor,npol,pwcpl_tp(ntor,npol)
                  enddo !npol
	       enddo !ntor
             endif ! i_grill_npar_ntor_npol_mesh.eq.1

             if(i_grill_npar_ntor_npol_mesh.eq.2) then
c------------------------------------------------------------------
c               n_toroidal mesh anztorin() will be chosen to get the equal
c		power pwcpl_t_ar(n) for all rays in the given power
c
c               n_poloidal mesh anzpolin() will be chosen to get the equal
c		power pwcpl_p_ar(n) for all rays in the given power
c
c               two dimentional array pwcpl_tp(nnktormax,nnkpolmax)
c               for power distribution on (N_toroidal,N_poloidal)
c               In this case it will pwcpl_tp(i,j)=Const for all rays
c------------------------------------------------------------------
                if(igrillpw.eq.2) then

c-----------------N_toroidal mesh calculations

                  call create_equal_mesh(-pi,pi,f_pow_npar_2,nnktori,
     &            anztorin_edge,anztorin,pwcpl_t_ar)

                  write(*,*)'andwdth_tor,antor0',andwdth_tor,antor0
                  do k=1,nnktori 
                    anztorin(k)=anztorin(k)/andwdth_tor+antor0
                    write(*,*)'k,anztorin(k)',k,anztorin(k)                  
                  enddo 

                  do k=1,nnktori+1                   
                     anztorin_edge(k)=anztorin_edge(k)/andwdth_tor+
     &                                 antor0
                    write(*,*)'k,anztorin_edge(k)',k,anztorin_edge(k)
                  enddo

c-----------------N_poloidal mesh calculations

                  call create_equal_mesh(-pi,pi,f_pow_npar_2,nnkpoli,
     &            anzpolin_edge,anzpolin,pwcpl_p_ar)

                  write(*,*)'andwdth_pol,anpol0',andwdth_pol,anpol0
                  do k=1,nnkpoli 
                    anzpolin(k)=anzpolin(k)/andwdth_pol+anpol0
                    write(*,*)'k,anzpolin(k)',k,anzpolin(k)                  
                  enddo 

                  do k=1,nnkpoli+1                   
                     anzpolin_edge(k)=anzpolin_edge(k)/andwdth_pol+
     &                                anpol0
                    write(*,*)'k,anzpolin_edge(k)',k,anzpolin_edge(k)
                  enddo

                endif ! igrillpw.eq.2

                if(igrillpw.eq.3) then

c-----------------N_toroidal mesh calculations

                  call create_equal_mesh(-1.d0,1.d0,f_pow_npar_3,nnktori
     &            ,anztorin_edge,anztorin,pwcpl_t_ar)

                  do k=1,nnktori 
                    anztorin(k)=anztorin(k)*0.5d0*difrntor+antor0
                    write(*,*)'k,anztorin(k)',k,anztorin(k)                  
                  enddo 

                  do k=1,nnktori+1                   
                     anztorin_edge(k)=anztorin_edge(k)*0.5d0*difrntor+
     &                                antor0
                    write(*,*)'k,anztorin_edge(k)',k,anztorin_edge(k)
                  enddo

c-----------------N_poloidal mesh calculations

                  call create_equal_mesh(-1.d0,1.d0,f_pow_npar_3,nnkpoli
     &            ,anzpolin_edge,anzpolin,pwcpl_p_ar)

                  do k=1,nnkpoli 
                    anzpolin(k)=anzpolin(k)*0.5d0*difrnpol+anpol0
                    write(*,*)'k,anzpolin(k)',k,anzpolin(k)                  
                  enddo 

                  do k=1,nnkpoli+1                   
                     anzpolin_edge(k)=anzpolin_edge(k)*0.5d0*difrnpol+
     &                                anpol0
                    write(*,*)'k,anzpolin_edge(k)',k,anzpolin_edge(k)
                  enddo

                endif ! igrillpw.eq.3

                do ntor=1,nnktori
                   pwcpl_t_ar(ntor)=1.d0/nnktori
                enddo

                do npol=1,nnkpoli
                   pwcpl_p_ar(npol)=1.d0/nnkpoli
                enddo

                do ntor=1,nnktori
                  do npol=1,nnkpoli
                    pwcpl_tp(ntor,npol)=pwcpl_t_ar(ntor)*
     &                                  pwcpl_p_ar(npol)
                  enddo
                enddo

                anorm=0.d0
 
                do ntor=1,nnktori
                  do npol=1,nnkpoli
                    anorm=anorm+pwcpl_tp(ntor,npol)                    
                  enddo
                enddo
c               write(*,*)'anorm='anorm

                do ntor=1,nnktori
                  do npol=1,nnkpoli
                 pwcpl_tp(ntor,npol)=pwcpl_tp(ntor,npol)/anorm*powers(i)
	           write(*,*)'ntor,npol,pwcpl_tp(ntor,npol)',
     &                        ntor,npol,pwcpl_tp(ntor,npol)
                  enddo !npol
	        enddo !ntor

             endif ! i_grill_npar_ntor_npol_mesh.eq.2

	   endif
20         continue
         endif !i_n_poloidal=4

c        calculation of the initial data for rays_iray
	 phi0=phigrill(i)
	 if(istart.eq.3) then ! already found as xconv,yconv,zconv
         x0=arxu0(i)
         y0=aryu0(i)
         z0=arzu0(i)
         r0=sqrt(x0**2+y0**2)
	 endif
c	 write(*,*)'phi0=',phi0
	 do j=1,nthini
            if (ilaunch.ne.1) then
               theta0=wtheta(j)
               if(istart.ne.3) then 
                 ! find  z0=z(psi0,theta0) and r0
                 if(model_rho_dens.eq.0) then
                    psi0=psi_rho(rhopsi0(i))
                    write(*,*)'grill_lh_xyz->zr_psith  #2'
                    call zr_psith(psi0,theta0,z0,r0)
                    x0=r0*cos(phi0)
                    y0=r0*sin(phi0)
                 elseif(model_rho_dens.eq.1 .or.model_rho_dens.eq.2 .or.
     +                  model_rho_dens.eq.3 .or.model_rho_dens.eq.4)then
                    call rho_theta_phi_xyz(rhopsi0(i),theta0,phi0,
     +                               x0,y0,z0) !->out
                    r0=dsqrt(x0*x0+y0*y0)
                 else
                   stop 'grill_lh_xyz: only for model_rho_dens=0,1,2,3'
                 endif
               endif ! istart.ne.3
            else ! ilaunch=1
               r0=r0launch
               z0=z0launch
               phigrill(1)=phi0launch*tet
               x0=r0*cos(phigrill(1))
               y0=r0*sin(phigrill(1))
            endif

            ! Get components of b and grad(psi):   [calls bxyz()]
            !call b_grpsi_xyz(x0,y0,z0, bx,by,bz,dpdxd,dpdyd,dpdzd) 
            bmod= bxyz(x0,y0,z0) !dsqrt(bx*bx + by*by + bz*bz)
            r0= dsqrt(x0**2+y0**2)
            grpsi=dsqrt(dpdxd*dpdxd+dpdyd*dpdyd+dpdzd*dpdzd)
            cs=x0/r0
            sn=y0/r0
            ! poloidal and toroidal magnetic field
            b_theta=(bz*(dpdxd*cs+dpdyd*sn)-(bx*cs+by*sn)*dpdzd)/grpsi
            bphi= -bx*sn + by*cs  

            if((i_n_poloidal.ge.1).and.(i_n_poloidal.le.3)) then 
c-------------i_n_poloidal=1,2,3
              do n=1,nnkpari
                write(*,*)'grill_lh: 1 n=',n
                iray=iray+1
                write(*,*)'grill_lh: iray',iray                
                if (iray.gt.nraymaxl) then
                 write(*,*)'in grill_lh iray=',iray,'nraymaxl',nraymaxl
                 write(*,*)'in grill_lh the total number of the rays is'
                 write(*,*)'greater than number of elements in arrays'
                 write(*,*)'iray.gt.nraymaxl'
                 write(*,*)'please change parameter nraymaxl in param.i'
                 stop
                endif

                arxu0(iray)=x0
                aryu0(iray)=y0
                arzu0(iray)=z0

c-----------------------------------------------------
                cnpar= anzin(n)

              if(i_n_poloidal.eq.1)then
                cnper_tang= 0.d0 !Nperp is purely along grad(psi) in this case
                cnteta= cnpar*(b_theta/bmod)
                cnphi=  cnpar*(bphi/bmod)
              endif

              if(i_n_poloidal.eq.3)then !Find Nper_tang from Nperp and angle
                cnper_tang= 0.d0 ! Nperp is unknown yet; 
                cnteta= cnpar*(b_theta/bmod)
                cnphi=  cnpar*(bphi/bmod)
              endif              ! cnper_tang will be found later
               
              if(i_n_poloidal.eq.2)then !Find Nper_tang from Npar and Ntheta
                cnteta= n_theta_pol ! given
                if(bphi .ne. 0.d0) then
                  cnper_tang= (cnteta*bmod - cnpar*b_theta)/bphi
                  cnphi=  (cnpar*bmod - cnteta*b_theta)/bphi
                else ! bphi=0 ->  In this case Npar== +/- N_theta_pol 
                  ! Note: n_toroidal is not known for this case. Set to 0?
                  cnper_tang= 0.d0 
                  cnphi= 0.d0
                  ! stop 'grill_lh_xyz: bphi=0; cnper_tang cannot be found'
                endif
              endif 

	        arnpar(iray)=cnpar
	        arnper_tang(iray)=cnper_tang
	        arnphi(iray)=cnphi
	        arntheta(iray)=cnteta
c-----------------------------------------------------
	        powinilh(iray)=fpwth(j)*pwcpl(n)

c calculation of the initial values of the n_parallel width
	        if (nnkpari.eq.1) then
                 write(*,*)'grill_lh nnkpari.eq.1'
                 wdnpar0(iray)=dabs(anzin(n))*0.05d0
	        else
                 wdnpar0(iray)=dabs(difrnpar)/dfloat(nnkpari)
	        endif
                write(*,*)'grill_lh  i_n_poloidal,iray,wdnpar0(iray)',
     +                               i_n_poloidal,iray,wdnpar0(iray)

	      end do !n
            endif !i_n_poloidal=1,2,3

            if(i_n_poloidal.eq.4)then ! Find Nper_tang from Ntor, Ntheta
c-------------i_n_poloidal=4
              do ntor=1,nnktori
              do npol=1,nnkpoli
                 iray=iray+1
                 if (iray.gt.nraymaxl) then
                 write(*,*)'in grill_lh iray=',iray,'nraymaxl',nraymaxl
                 write(*,*)'in grill_lh the total number of the rays is'
                 write(*,*)'greater than number of elements in arrays'
                 write(*,*)'iray.gt.nraymaxl'
                 write(*,*)'please change parameter nraymaxl in param.i'
                 stop
	           endif

                 cnpar= 
     =              (anzpolin(npol)*b_theta +anztorin(ntor)*bphi)/bmod
                 cnper_tang= 
     =              (anzpolin(npol)*bphi -anztorin(ntor)*b_theta)/bmod

                 arxu0(iray)=x0
                 aryu0(iray)=y0
                 arzu0(iray)=z0
                 arnpar(iray)=cnpar
                 arnper_tang(iray)=cnper_tang
                 arnphi(iray)=anztorin(ntor)
                 arntheta(iray)=anzpolin(npol)
	           powinilh(iray)=fpwth(j)*pwcpl_tp(ntor,npol)
c Setting the initial values of the n_parallel width
                 wdnpar0(iray)=hnwidth

                 write(*,*)'grill_lh i_n_poloidal,iray,wdnpar0(iray)',
     +                               i_n_poloidal,iray,wdnpar0(iray)
	         end do !npol
               end do !ntor
            endif !i_n_poloidal=4
              
	 end do !j
      end do !i=1,ngrill
      nray=iray		   

      powert=0.d0
      do iray=1,nray
         powert=powert+powinilh(iray)
      enddo
      write(*,*)'grill_lh_xyz nray,powert',nray,powert
      
c--------------------------------------------------------------------
      return
      end



c======================================================================
c======================================================================
      subroutine b_grpsi_xyz(x,y,z, 
     +           b_x,b_y,b_z,dpdx,dpdy,dpdz) !->out
c     This subroutine is called by grill_lh_xyz 
c     to get components of b and grad(psi).
c     A direct call of bxyz function from grill_lh_xyz
c     is prevented by conflict between arguments of grill_lh_xyz
c     and arrays in 'one.i'.
      implicit none
      include 'param.i'
      include 'one.i' ! contains components of b and grad(psi)
      include 'three.i'
      include 'five.i'
      include 'six.i'
c-----Input:
      real*8 x,y,z ! cartesian coords
c-----Output: components of b and grad(psi)
      real*8 b_x,b_y,b_z, dpdx,dpdy,dpdz
      external bxyz
      real*8 bxyz
      bmod= bxyz(x,y,z) !-> get bx,by,bz
      b_x=bx
      b_y=by
      b_z=bz
      dpdx=dpdxd
      dpdy=dpdyd
      dpdz=dpdzd
      return
      end
c======================================================================
c======================================================================

      subroutine owconvr_xyz(theta,phi,rmn,rmx,wpw2in,
     + xconv,yconv,zconv) !->out
c     It determines the point (xconv,yconv,zconv) 
c     where wpw2()=wpw2in == (omega_pe/omega)**2.
c     For wpw2in=1, it gives the point of O_X conversion. 
c     The search of the point is made along the line originating from
c     (x,y,z)=0, at a given (fixed) theta and given phi angle;
c     poloidal angle theta (degree); theta is from thgrill() input. 
c     if theta=0  point is on the outer part of the magnetic surface
c     if theta=90 point is at the top of the magnetic surface
c     The search is done among r=sqrt(x^2+y^2) such that rmn < r < rmx.
c     input data:
c          wpw2(r=rmx) < wpw2in < wpw2(r=rmn)  
c          a given value wpw2in == (omega_pe/omega)**2
c          theta is a poloidal angle (degree)
c          phi is a toroidal angle (degree)
c     output:
c          the point (xconv,yconv,zconv)  
c          where wpw2()=wpw2in
c----------------------------------------------------------------------
c     It solves the equation 
c     wpw_2(xconv(r),yconv(r),zconv(r),1) - xe0 = 0
      implicit none
      double precision theta,phi,rmn,rmx,wpw2in ! Input
      double precision racc , r
      double precision xconv,yconv,zconv ! out
      double precision thetax,phix,xe0 
      double precision rtbis,xe_eq_xyz
      common /convers/ thetax,phix,xe0
      external xe_eq_xyz

      thetax=theta*datan(1.d0)*4.d0/180.d0 ! poloidal angle in radian
      phix=phi*datan(1.d0)*4.d0/180.d0 ! toroidal angle in radian
      
      xe0= wpw2in ! used by xe_eq
      racc = 1.d-7
      r= rtbis(xe_eq_xyz,rmn,rmx,racc) !-> find such r that xe_eq(r)~0
      xconv= r*cos(phix)*cos(thetax)
      yconv= r*sin(phix)*cos(thetax)
      zconv= r*sin(thetax)
      return
      end

c======================================================================
c======================================================================
      double precision function	xe_eq_xyz(r)
c-----The function xe_eq= wpw_2(x,y,z,1)-xe0
c     input:
c     r==sqrt(x^2+y^2)  
c     thetax poloidal angle 
c     phix   toroidal angle
c     xe0     the given xe value
c     output:
c     xe_eq= wpw_2(x,y,z,1)-xe0
c--------------------------------------
      implicit none
      double precision r
      double precision thetax,phix,xe0
      common /convers/thetax,phix,xe0  ! as parameters for xe_eq_xyz
      double precision xx,yx,zx
      double precision wpw_2, xp
      xx= r*cos(phix)*cos(thetax)
      yx= r*sin(phix)*cos(thetax)
      zx= r*sin(thetax)
      xp=wpw_2(xx,yx,zx,1)
      xe_eq_xyz= xp-xe0
      return
      end

c======================================================================
c======================================================================

c        **********************_dinit_mr_xyz ************************
c        *                        -                           *
c        * this subroutine reads the data from genray.in     *
c        ******************************************************
c
c-------------------------dinit_mr_xyz---------------------------------
c        It creates data for multiple ray case. 
c        for spline approximations                       	   
c        density,temperature and Z_effective profiles.             
c        Calculates initial data for all rays.                     
c        It uses input data from genray.in or genray.dat 
c        file for multiple ray case.                               
c        Data from genray.in or genray.dat file were read          
c        previously in genray.f file   
c        It reads or creates the non-maxwellian electron distribution
c
c        it uses the following functions and subroutines           
c        bmod,spldens1                                             

c      	 This program reads data from genray.in file for	   !
c        multiple ray case.                                        !
c        It creates the files for spline approximations 	   !
c        density,temperature and Z_effective profiles.             !
c        it uses the following functions and subroutines           !
c        bmod,spldens1                                             !
c------------------------------------------------------------------
c        output parameters:					   !
c                          ndim-number of the ray-tracing equations!
c                          nray-number of the rays at antenna      !
c------------------------------------------------------------------

      subroutine dinit_mr_xyz(ndim,nray)

      implicit none

      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'three.i'
      include 'five.i'
      include 'cone.i'
      include 'grill.i'
      include 'rkutta.i'
      include 'six.i'
      include 'write.i'
      include 'writencdf.i'
      include 'onetwo.i'
      include 'output.i'
      !!! include 'scatnper.i'
      !!! include 'emissa.i'
      double precision thetax,phix,xe0
      common /convers/ thetax,phix,xe0

c-----output
      integer
     &ndim,  !number of the ray-tracing equations
     &nray   !number of the rays at antenna
c..............................................................
c     these two arrays are for namelist work only
c      dimension prof2(nbulka,ndensa)
c      dimension prof(nbulka*ndensa)
c     Following are working arrays for the namelist input of density 
c     and temperature,... at nonuniform radii grids
      real*8 prof2_nonuniform(nbulka,ndensa),prof_radii(nbulka,ndensa)
c      integer nj_tab(nbulka) !the number of profile points for each species
c..............................................................
      include 'dinit_nml.i'
c..............................................................
      integer nbulk1,i,j,k,iraystop,i1,j1,ifreq,ii,initial,nray1,icone,
     &imax

      real*8 trnspi,h,v0,w0,hfreq,delt,zefftest,zion,psi,denstot,
     &szini,szi2ni,stini,pressure,den_test,prestest,tem_test,
     &energy,pitch,fdist,dfdx,dfdpitch,dfdp,psi_mag,xi,yi,tetan,wpw2in,
     &zconv,rconv,theta,rhoconv,phi,xconv,yconv,cnparopt,
     &dens,temp,psi_loc,
     &tpop,
cSAP090311 for test
     &zeff_loc,vflow_loc

c-----externals
      real*8 bxyz,psi_rho,prespsi,densrho,temperho,
     + wpw_2,wcw,tpoprho,
     & zeffrho,vflowrho, psif_xyz
c-----locals
      real*8 wpw2_conv, wcw_conv, rmn,rmx

c-----genray input file: genray.in or genray.dat
      character*10 genray_in_dat

      iraystop=0
      pi=4*datan(1.d0)
      trnspi=pi/180.d0

      write(*,*)'dinit_mr_xyz: rst,zst=',rst(1),zst(1)
      !pause

c---------------------------------------------------------------
c     read all namelists from input genray.in or genray.dat file
c---------------------------------------------------------------
cSAP080731
c      call read_all_namelists(genray_in_dat,ndim,nray)
c      write(*,*)'in dinit_mr, genray_in_dat=', genray_in_dat
c-----If input file is genray.in, then it will change
c     input data to genray.dat file format 
c      if (genray_in_dat.eq.'genray.in')  then
c         write(*,*)'genray.f before transform_genray_in_to_dat'
c         call transform_genray_in_to_dat
c      endif
c---------------------------------------------------------------------
      if (n_wall.gt.0) then    
c-------calculate poloidal angles [at radian] and small radius
c       of wall and limiter points
c       thetapol_wall(i=1,..,n_wall)
c       thetapol_limiter(i=1,...,n_limiter(j),j=1,max_limiters)
c       rho_wall(i=1,..,n_wall)
c       rho_limiter(i=1,..,n_limiter(j),j=1,max_limiters))
c  
c       These arrays will be in common /fourb/ in file fourb.i
        call wall_limiter_theta_pol_rho
      endif
c------------------------------------------------------------------------
      write(*,*)'dinit_mr_xyz: Absorption'
      if(iabsorp.le.4  .or.
     +   iabsorp.eq.6  .or.  iabsorp.eq.7  .or.  iabsorp.eq.12) then
      ! continue
      else
      stop 'Not setup for this value of iabsorp. Choose from 1-4,6,7,12'
      endif     
c     iabsorp=2 for LH waves
c     iabsorp=3 for FW waves
c-----------------------------------------------------
      write(*,*)'dinit_mr_xyz: n_relt_harm1',n_relt_harm1
      write(*,*)'dinit_mr_xyz: n_relt_harm',n_relt_harm
      write(*,*)'dinit_mr_xyz: n_relt_harma',n_relt_harma
      if (n_relt_harm1.eq.9999)then
        n_relt_harm1=-n_relt_harm
        n_relt_harm2= n_relt_harm
      else
        n_relt_harm2=n_relt_harm1+n_relt_harm-1
      endif
      print*,'dinit_mr_xyz: n_relt_harm1,2',n_relt_harm1,n_relt_harm2
      print*,'n=omega/omega_ce along rays should be within range above'
      !pause !!!
       
      if(n_relt_harm1.lt.n_relt_harm1a) then
         write(*,*)'dinit_mr_xyz: n_relt_harm1<n_relt_harm1a'
         write(*,*)'it should be n_relt_harm1=>n_relt_harm1a'
         write(*,*)'please change n_relt_harm1a'        
         write(*,*)'in param.i and recompile the code'
         write(*,*)'n_relt_harm1,n_relt_harm1a',
     &              n_relt_harm1,n_relt_harm1a
         stop 'in dinit_mr_xyz:' 
      endif
 
      if(n_relt_harm2.gt.n_relt_harma) then
         write(*,*)'dinit_mr_xyz: n_relt_harm2>n_relt_harma'
         write(*,*)'it should be n_relt_harm2=<n_relt_harma'
         write(*,*)'please change n_relt_harma'
         write(*,*)'in param.i and recompile the code'
         write(*,*)'n_relt_harm2,n_relt_harma',
     &              n_relt_harm2,n_relt_harma
         stop 'in dinit_mr_xyz:' 
      endif

      do i=1,nbulk 
        if ((i_salphal(i).ne.1).and.(i_salphal(i).ne.0)) then
          write(*,*)'(i_salphal(i).ne.1).or.(i_salphal(i).ne.0)'
          write(*,*)'dinit_mr_xyz: It should i_salphal(i) =0 or =1'
          write(*,*)'i=',i,'i_salphal(i)',i_salphal(i)
          write(*,*)'Change i_saplhal(i)'
          write(*,*)'in input genray.in or genray.dat file'
          stop 'in dinit_mr_xyz: i_salphal problem'
        endif
      enddo 
           
      write(*,*)'dinit_mr_xyz: Plasma parameters'
      write(*,*)'model_rho_dens=',model_rho_dens,'    idens=', idens
      if(model_rho_dens.ne.0  .and. idens.ne.0) then
        stop 'dinit_mr_xyz: model_rho_dens>0 only works with idens=0'
      endif

      if (model_rho_dens.eq.1) then
        write(*,'(a,3e12.4)')
     +       'Ellipsoid Center:   elx0,ely0,elz0=',elx0,ely0,elz0
        write(*,'(a,3e12.4)')
     +       'Semi-axes [meters]: elax,elay,elaz=',elax,elay,elaz
        if(elax*elay*elaz .eq. 0.d0) then ! all should be >0.
          stop 'Set elax,elay,elaz to positive real values'
        endif
      endif
            
      if (model_rho_dens.eq.2) then
        write(*,'(a,2e12.4)')
     +       'Ellipse Center [m]:   elx0,ely0=',elx0,ely0
        write(*,'(a,2e12.4)')
     +       'Semi-axes [meters]: elax,elay=',elax,elay
        if(elax*elay .eq. 0.d0) then ! all should be >0.
          stop 'Set elax,elay to positive real values'
        endif
      endif

      if (model_rho_dens.eq.4) then
        write(*,'(a,e12.4)')
     +       'Radius of wall [m]:   wall_rmax=',wall_rmax
        write(*,'(a,e12.4)')
     +       'Separatrix radius [meters]: rs_frc=',rs_frc
        if(rs_frc .ge. wall_rmax) then ! 
          stop 'Set rs_frc to be less than wall_rmax'
        endif
      endif

      if (model_rho_dens.eq.5) then
         if(eqdsktype.ne.'TAE')then
           stop 'model_rho_dens=5 can only be used for eqdsktype="TAE" '
         endif
         if(model_b.ne.0)then
           stop 'model_rho_dens=5 can only be used for model_b=0(eqdsk)'
         endif
         if(nbulk.gt.1)then
           stop 'model_rho_dens=5 is only checked for nbulk=1 '
         endif
         if(idens.eq.1)then
           write(*,*)'model_rho_dens=5 is only checked for idens=0 :'
           write(*,*)'ne is read from eqdsk, but Te - analytically.'
           stop 
         endif
      endif
      
      write(*,*)'dinit_mr_xyz: izeff=',izeff
      do i=1,nbulk
         write(*,*)'i,temp_scale(i),den_scale(i)',
     .   i,temp_scale(i),den_scale(i)
      enddo 

      do i=1,nbulk
         te0(i)=ate0(i)
         teb(i)=ateb(i)
      enddo
c-----------------------------------------------------
c      /species/
c   plasma component charges charge(i)=mod(charge(i)/charge_electron)
c----------------------------------------------------

c-----------------------------------------------------
c  plasma components mass dmas(i)=Mass(i)/Mass_electron
c-----------------------------------------------------
      do i=1,nbulk
         write(*,*)'dinit_mr_xyz: i, dmas(i)',i,dmas(i)
      enddo
c--------------------------------------------------------------

c     calculation of nbulk1
      if(((izeff.eq.0).or.(izeff.eq.2)).or.(izeff.eq.3)) then
c        izeff=0, zeff will be calculated using the given ions;
c                 the electron density will be calculated using ion's densities;
c             =1  ion's densities(i) i=nbulk and i=nbulk-1 will be calculated  using
c                 Zeff, the electon density and ion's densities(i), i=2,nbulk-1;
c        izeff=2, zeff will not coincide with the plasma components
c             =3  it uses eqdsk pres (pressure) and ions densities_i
c                 for i=2,... nbulk
c                 Let temperature T_E=T_i
c                 pres=dens1(k,1)*temp1(k,1)+
c                      Sum(i=2,nbulk)(dens1(k,i)*temp1(k,i)
c                 In this case we will calculate Zeff(rho),
c                 dens_electron(rho) and T_e(rho)=T_i(rho)
c             =4  it uses eqdsk pres (pressure), zeff,ions densities
c                 for i=2,... nbulk-2 (nbulk>2) and T_e(rho),T_i(rho)
c                 pres=dens1(k,1)*temp1(k,1)+
c                      Sum(i=2,nbulk)(dens1(k,i)*temp1(k,i)
c                 In this case we will calculate dens_electron(rho) and
c                 ion densities for i=nbulk and i=nbulk-1
         nbulk1=nbulk         
      else
c        (izeff=1 or izeff=4), zeff is given, the ions component will be calculated
         if (nbulk.le.2) nbulk1=nbulk
         if (nbulk.eq.2) then
         write(*,*)'dinit_mr_xyz: nbulk=2, Zeff must be equal charge(2)'
         write(*,*)'Please check it or use the option izeff=0'
         endif
         if (nbulk.gt.2) nbulk1=nbulk-2
      endif !izeff
      write(*,*)'dinit_mr_xyz: nbulk1=',nbulk1
c------------------------------------------------------------------
      h=1.d0/(ndens-1)
      do i=1,ndens
        !rhom(i)=h*(i-1) ! should be defined by now
        write(*,*) 'dinit_mr_xyz: i, rhom(i)=', i, rhom(i)
      enddo
c------------------------------------------------------------------

c------------------------------------------------------------------
      if(idens.eq.0) then        
         if(model_rho_dens.eq.0 .or. model_rho_dens.eq.5) then
c-dense_(i)=(dense0(i)-denseb(i))*(1-rho**rn1de(i))**rn2de(i)+denseb(i)
c------------------------------------------------------------------\   
            if (model_rho_dens.eq.0) i1=1 ! analytical for all species  
            if (model_rho_dens.eq.5) i1=2 ! analytical for all except e
            if(i1.le.nbulk1)then 
            write(*,*)'dinit_mr_xyz: Analyt.radial profiles for dens'
            do i=i1,nbulk1
               write(*,*)'i, dense0(i)',i,dense0(i)
               write(*,*)'i, denseb(i)',i,denseb(i)
               write(*,*)'i, rn1de(i)',i,rn1de(i)
               write(*,*)'i, rn2de(i)',i,rn2de(i)
               write(*,*)'i, dense0(i)',i,dense0(i)
            enddo
c--------creation the array dens1(ndensa,nbulka)
            do i=i1,nbulk1
            do k=1,ndens
               rho=rhom(k)
               dens1(k,i)=(dense0(i)-denseb(i))*
     1	                  (1-rho**rn1de(i))**rn2de(i)+denseb(i)
            enddo
            enddo

            do i=nbulk1,i1,-1
            do k=1,ndens
               rho=rhom(k)
               if (((izeff.eq.0).or.(izeff.eq.3)).and.(i.eq.1)) then
                  dens1(k,1)=0.d0
                  do j=2,nbulk
                    dens1(k,1)=dens1(k,1)+charge(j)*dens1(k,j)
                  enddo
               else
                  dens1(k,i)=(dense0(i)-denseb(i))*
     1	                     (1-rho**rn1de(i))**rn2de(i)+denseb(i)
c-----------------multiply density profiles by den_scale
                  dens1(k,i)=dens1(k,i)*den_scale(i)
                  write(*,*)'dinit i,k,dens1(k,i)',i,k,dens1(k,i)
                  write(*,*)'dense0(i),denseb(i),rho,rn1de(i),rn2de(i)',
     *            dense0(i),denseb(i),rho,rn1de(i),rn2de(i)
               endif
            enddo
            enddo
            write(*,*)'end of analytical density profiles input'
            endif ! i1.le.nbulk1
         endif ! model_rho_dens=0,5

c------------------------------------------------------------------         
c         /tpopprof/
c         /vflprof/
c--------------------------------------------------------------------       
c        creation the arrays for analytical profiles
c        tpop1(ndensa,nbulka), vflow1(ndensa,nbulka)
	 do i=1,nbulk
	    do k=1,ndens
	       rho=rhom(k)	
	       tpop1(k,i)=(tp0(i)-tpb(i))*
     1                    (1-rho**rn1tp(i))**rn2tp(i)+tpb(i)

               vflow1(k,i)=(vfl0(i)-vflb(i))*
     1                    (1-rho**rn1vfl(i))**rn2vfl(i)+vflb(i)
	    enddo
	 enddo

c---------------------------------------------------------------------
c        /zprof/          
	 if (izeff.eq.3) goto 10
c        /tprof/
c----------------------------------------------------------
cTemperature
ctempe_(i)=(te0(i)-teb(i))*(1-rho**rn1te(i))**rn2te(i)+teb(i)
c-----------------------------------------------------------
	 do i=1,nbulk
           write(*,*)'i, te0(i)',i,te0(i)	
           write(*,*)'i, teb(i)',i,teb(i)
           write(*,*)'i, rn1te(i)',i,rn1te(i)
           write(*,*)'i, rn2te(i)',i,rn2te(i)
           write(*,*)'i,tp0(i),tpb(i)',i,tp0(i),tpb(i)
           write(*,*)'rn1tp(i),rn2tp(i)',rn1tp(i),rn2tp(i)
	 enddo
c------- creation of array temp1(ndensa,nbulka)
	 do i=1,nbulk
          do k=1,ndens
             rho=rhom(k)
             temp1(k,i)=(te0(i)-teb(i))*
     1                    (1-rho**rn1te(i))**rn2te(i)+teb(i)
cSm070426
c--------------multiply temperature profiles by temp_scale
             temp1(k,i)=temp1(k,i)*temp_scale(i)
          enddo
	 enddo
10	 continue         


         write(*,*)'zeff0,zeffb,rn1zeff,rn2zeff'
         write(*,*)zeff0,zeffb,rn1zeff,rn2zeff
	 if(((izeff.eq.1).or.(izeff.eq.2)).or.(izeff.eq.4)) then
c           the given analytical Zeff profile
c-------------------------------------------
c           zeff=(zeff0-zeffb)*(1-rho**rn1zeff)**rn2zeff+zeffb
c-------------------------------------------
c           the creation of array zeff1(ndens)
	    do k=1,ndens
               rho=rhom(k)
               zeff1(k)=(zeff0-zeffb)*(1-rho**rn1zeff)**rn2zeff+zeffb
	    enddo
	 endif

	 if((izeff.eq.0).or.(izeff.eq.3)) then
	    write(*,*)'izeff=1 zeff will be calculated using the given
     1	    ions densities'
	 endif
      endif ! idens=0  analytical
c--------------------------------------------------------------------

c--------------------------------------------------------------
      if(idens.eq.1) then
c        spline approximation of the density, temperature, zeff,
c        tpop=T_perp/T_parallel, vflow
c        radial profiles
c        input of the arrays on the radial mesh	from dtzprof.dat
       write(*,*)'dinit_mr_xyz: idens=1 Spline radial profiles (dentab)'
       write(*,*)'nbulka,ndensa,nbulk,ndens',nbulka,ndensa,nbulk,ndens
         write(*,*)'dinit_mr_xyz: before mult by den_scale'
         do i=1,nbulk
            write(*,*)'i',i,'dens1(k,i)'
            write(*,*)(dens1(k,i),k=1,ndens)
         enddo
c---------multiply density profiles by den_scale(i)
cBh080102          do i=i1,nbulk          
         do i=1,nbulk
         do k=1,ndens
            dens1(k,i)=den_scale(i)*dens1(k,i)
         enddo
         enddo
        
         write(*,*)'dinit_mr_xyz: nbulk1',nbulk1,'ndens',ndens
         if ((izeff.eq.0).or.(izeff.eq.3)) then
            i1=2
         else
            i1=1
         endif

         write(*,*)'dinit_mr_xyz: after mult by den_scale' 
         do i=i1,nbulk1
            write(*,*)'i',i,'dens1(k,i)'
            write(*,*)(dens1(k,i),k=1,ndens)
         enddo
         
c----------------------------------------------------------
c        calculation of the electron density from
c        the charge neutrality 
c----------------------------------------------------------
         if ((izeff.eq.0).or.(izeff.eq.3)) then
           do k=1,ndens
              dens1(k,1)=0.d0
           do j=2,nbulk
              dens1(k,1)=dens1(k,1)+charge(j)*dens1(k,j)
           enddo
           enddo
         endif

         write(*,*)'dinit_mr_xyz: dens1(k,1)'
         write(*,*)(dens1(k,1),k=1,ndens)

c--------multiply temperature profiles by temp_scale
         do i=1,nbulk
         do k=1,ndens
c               write(*,*)'i,k,temp1(k,i),temp_scale(i)',
c     &                    i,k,temp1(k,i),temp_scale(i)
               temp1(k,i)=temp1(k,i)*temp_scale(i)
c               write(*,*)'temp1(k,i)',temp1(k,i)
         enddo
         enddo

         do i=1,nbulk
            write(*,*)'i',i,'temp1(k,i)'
            write(*,*) (temp1(k,i),k=1,ndens)
         enddo

c---------tpoptab
         do i=1,nbulk
            write(*,*)'i',i,'tpop1(k,i)'
            write(*,*) (tpop1(k,i),k=1,ndens)
         enddo

c--------vflowtab
         do i=1,nbulk
            write(*,*)'i',i,'vflow1(k,i)'
            write(*,*) (vflow1(k,i),k=1,ndens)
         enddo

         write(*,*)'dinit_mr_xyz: izeff=',izeff

         if(((izeff.eq.1).or.(izeff.eq.2)).or.(izeff.eq.4)) then
c           the given Zeff profile is in the table form
            write(*,*)'zeff1(k)'
         else
            write(*,*)'uniform zeff1',zeff1
         endif !izeff              
      endif ! idens=1
c-----------------------------------------------------------


 20   continue 
      write(*,*)'dinit_mr_xyz after 20_continue finished idens=',idens
      !pause !!!
c------------------------------------------------------------

c---------------------------------------------------------
c     read the data for for EC cone vertex coordinates calculations
c     This case is used for the optimal OX mode conversion.
      
c      /ox/


      if(i_ox.eq.1) then
        istart=3
        prmt(3)=-prmt3 !to create the negative time step of integration
        i_vgr_ini=+1
        ireflm=1   
        write(*,*)'dinit_mr_xyz: i_ox.eq.1 prmt(3)',prmt(3)
      endif

      if(((i_ox.ne.0).and.(i_ox.ne.1)).and.(i_ox.ne.2)) then
         write(*,*)'i_ox can  =0 or =1 or =2'
         write(*,*)'in namelist /ox/ i_ox=',i_ox
         write(*,*)'please change i_ox in input file'
         stop 'in dinit_mr_xyz: /ox/'
      endif  
         
    

c------------------------------------------------------------
c       no emission calculations, set one frequency
        nfreq=1
c-----------------------------------------------
c     for the emission multi frequency case
c     calculate the electron gyro-frequency freqncy0 at the plasma center
c---------------------------------------------
      bmod=bxyz(rma,0.d0,zma) !TL  ! YuP Needs work? Dependence on phi?
      freqncy0=28.0*b0*bmod !GHZ

      if (nfreq.eq.1) then
c--------no emission or only one frequency in the emission calculations   
         v0=806.2/(frqncy*frqncy)
         w0=28.0*b0/frqncy 
c-----------------------------------------------------
c        set the arrays for v and w
         i=1 ! electrons
         v(i)=v0
         w(i)=w0      
         write(*,*)'dinit_mr_xyz: i charge(i),dmas(i),v(i),w(i)'
     +      ,i,charge(i),dmas(i),v(i),w(i)
         do i=2,nbulk
            v(i)=v0*charge(i)**2/dmas(i)
            w(i)=w0*charge(i)/dmas(i)
            write(*,*)'dinit_mr_xyz: i charge(i),dmas(i),v(i),w(i)'
     +      ,i,charge(i),dmas(i),v(i),w(i)
         enddo
      else
      endif
     
c-------------------------------------------------------------
      if ((izeff.eq.0).or.(izeff.eq.3))then
c---------------------------------------------------------
c        calculation of the table for the radial profile zeff1(ndens)
         call zeffcalc
c        zeff1(ndens) is in common six.i
	 do i=1,ndens
	   write(*,*)'i',i,'zeff1(i)',zeff1(i)
         enddo
      endif !izeff=0
c---------------------------------------------------------
      if(izeff.eq.1) then
c---------------------------------------------------------
c        calculation of the table for the ion densities profiles
c---------------------------------------------------------
         if (nbulk.lt.3) then
            write(*,*)'dinit_mr_xyz: nbulk.lt.3, 
     1       Zeff must be equal charge(2),
     1       control it and use the option izeff=0'
            stop
         else
c           nbulk.ge.3
c           calculation of the tables for the radial profile
c           dens1(ndens,nbulk) and dens1(ndens,nbulk-1)
            if( charge(nbulk).eq.charge(nbulk-1)) then
              write(*,*)'Warning in dinit: nbulk(.ge.3)=',nbulk
              write(*,*)'in dinit: charge(nbulk)=charge(nbulk-1)'
              write(*,*)'it is impossible to find the ions densities'
              write(*,*)'change charge(nulk) or charge(nbulk-1)'
              write(*,*)'it should be charge(nulk)>charge(nbulk-1)'
              write(*,*)'or use the option izeff=0'
              stop
            endif

            call denscalc
c------------------------------------------------
c for test
            do i1=1,nbulk
               write(*,*)'dinit_mr_xyz: after call denscalc i1=',i1
               do j1=1,ndens
	          write(*,*)'j1=',j1,'dens1(j1,i1)',dens1(j1,i1)
               enddo
	    enddo
	    do j1=1,ndens
	       zefftest=0.d0
	       zion=0.d0
	       do i1=2,nbulk
	          if(dens1(j1,1).ne.0.d0) then
	             zefftest=zefftest+(dens1(j1,i1)/dens1(j1,1))*
     1                        charge(i1)*charge(i1)
		     zion=zion+charge(i1)*dens1(j1,i1)/dens1(j1,1)
		  else
		     write(*,*)'dinit_mr_xyz: dens1(j1,1)=0'
		     zefftest=zefftest+1.d0*charge(i1)*charge(i1)
		  endif
	       enddo
	       write(*,*)'j1',j1,'zefftest',zefftest,'zion',zion
	    enddo
c end test
c------------------------------------------------

	 endif ! nbulk
      endif ! izeff=1
      
      if (izeff.eq.3) then
	 do i=1,nbulk
	    do k=1,ndens
	       rho=rhom(k)
	       psi=psi_rho(rho)
	       denstot=dens1(k,1)
	       do ii=2,nbulk
                denstot=denstot+dens1(k,ii)
	       enddo
	       temp1(k,i)=prespsi(psi)/denstot/(1.6d3)
	       write(*,*)'dinit_mr_xyz: izeff=3,k,rho',k,rho
	       write(*,*)'prespsi(psi),denstot',prespsi(psi),denstot
	       write(*,*)'dinit_mr_xyz: i,temp1(k,i)',i,temp1(k,i)
	       write(*,*)'dinit_mr_xyz: i,dens1(k,i)',i,dens1(k,i)
	       write(*,*)'dinit_mr_xyz: i,dens1(k,1)',i,dens1(k,1)
	    enddo
	 enddo
      endif !izeff=3

      if(izeff.eq.4) then
cSAP090801 
c        if(nbulk.lt.3) then
c	  write(*,*)'in dinit: izeff=4 nbulk=',nbulk
c	  write(*,*)'in dinit: it should be nbulk>2, use another '
c	  write(*,*)'in dinit: izeff option '
c	  stop
c	else

c         nbulk.ge.3
c         calculation of the tables for the radial profiles
c         dens1(ndens,nbulk),dens1(ndens,nbulk-1) and
c         dens1(ndens,1)
	  call denscalp
c test izeff=4	beg
cSAP090801
        if (nbulk.ge.3) then
cyup	  do j=1,ndens
cyup	    szini=0.d0	   !sum(i=2,nbulk){charge(i)*dens1(j,i)}
cyup	    szi2ni=0.d0	   !sum(i=2,nbulk){charge(i)**2*dens1(j,i)}
cyup	    stini=0.d0	   !sum(i=2,nbulk){temp1(j,i)*dens1(j,i})
cyup	    rho=rhom(j)
cyup	    psi=psi_rho(rho)
cyup	    pressure=prespsi(psi)/1.6d3
cyup	    do i=2,nbulk
cyup	       szini=szini+charge(i)*dens1(j,i)
cyup	       szi2ni=szi2ni+charge(i)*charge(i)*dens1(j,i)
cyup	       stini=stini+temp1(j,i)*dens1(j,i)
cyup	    enddo
cyup	    prestest=temp1(j,1)*dens1(j,1)+stini
cyup	    zefftest=szi2ni/dens1(j,1)
cyup	    write(*,*)'in dinit izeff=4,j,rho',j,rho
cyup	    write(*,*)'pressure,prestest',pressure,prestest
cyup	    write(*,*)'zeff,zefftest',zeff1(j),zefftest
cyup	    write(*,*)'dens1(j,1),szini',dens1(j,1),szini
cyup	  enddo !j
c test izeff=4	end
	endif !nbulk.ge.3
      endif !izeff=4

 30   continue ! if (partner.eq. ....) goto 20

c     creation of the density,temperature,zeff and
c     tpop, vflow
c     spline coefficients
      call spldens1
      write(*,*)'dinit_mr_xyz: after spldens1'
c-----test printing density,temperature,zeff
cyup      do j=1,ndens
cyup         write(*,*)'j,rhom(j)',j,rhom(j)
cyup         do i=1,nbulk
cyup           den_test=densrho(rhom(j),i)
cyup           tem_test=temperho(rhom(j),i)
cyup           write(*,*)'i,dens(i,rhom(j)),temp(i,rhom(j))',
cyup     &     i,den_test,tem_test
cyup         enddo
cyup      enddo

      do i=1,ncone ! from genray.in
         phist(i)=phist(i)*trnspi ! degree to radians
         betast(i)=betast(i)*trnspi
         alfast(i)=alfast(i)*trnspi
      enddo
c--------------------------------------------------------------
c  normalization of the start parameters
      do i=1,ncone
         zst(i)=zst(i)/r0x ! from genray.in
         rst(i)=rst(i)/r0x ! from genray.in
         xst(i)=rst(i)*dcos(phist(i))
         yst(i)=rst(i)*dsin(phist(i))
      enddo
c  normaliszation 'diskdisk' parameters  
      d_disk=d_disk/r0x
      rho_launching_disk=rho_launching_disk/r0x
      rho_focus_disk=rho_focus_disk/r0x
      sigma_launching_disk=sigma_launching_disk/r0x
c-----------------------------------------------------------
      psi_mag=psif_xyz(rma,0.d0,zma)
      write(*,*)'dinit_mr_xyz: rma,zma,psi_mag,psilim=',
     +    rma,zma,psi_mag,psilim
ctest
      if (nfreq.eq.1) then
c--------no emission or only one frequency in the emission calculations   
         ! YuP: is xi,yi here only for print-out?
         write(*,*)'dinit_mr_xyz: 1'
         xi=wpw_2(rma,0.d0,zma,1) ! Needs work? Dependence on phi?
         write(*,*)'dinit_mr_xyz: 2'
         yi=  wcw(rma,0.d0,zma,1) ! Needs work? Dependence on phi?
         write(*,*)'dinit_mr_xyz: at magnetic axis Xe,Ye ',xi,yi
c---------------------------------------------------------------
c     the creation the data for the contours X_e=const,Y_e=const
c     B_tot,B_tor, B_pol on the plate (rho,theta) 
c     These data will have the form of the tables in 
c     xybrhoth.dat: rho(i),theta(j),xe,ye,(xe+ye*ye),bmod,bphi,
c     *              dsqrt(bz**2+br**2)
c      nrhomap=30
c      nthetmap=100
c      call mapxyb(nrhomap,nthetmap)
c      write(*,*)'mapxyb was created'
c      stop
      endif
c---------------------------------------------------------------
c     if ray start point is outside the plasma (ECR -case)
c     then:
c     determination of arrays :1)for ray coordinates on the ECR cone
c     alphaj(nray),betaj(nray)-spherical and cylindrical angles
c                              2)for wave power	angle distribution
c     powj(nray) -power flowing in the ray chanel at antenna
c     with normalized  Sum(i=1,nray)delpw0(i)=1
c---------------------------------------------------------------
      if (istart.eq.1) then
c---------EC waves---------------------------------
         write(*,*)'dinit_mr_xyz: istart=',istart
         write(*,*)'raypatt=',raypatt
         write(*,*)'ncone=',ncone

         nray1=1                ! Counter for position in ray arrays
         do icone=1,ncone
            if (raypatt.ne.'toray') then
            else
               tetan=betast(icone)
            endif
            write(*,*)'icone=',icone
            write(*,*)'alpha1,na1,na2,alpha2,phist,alfast,tetan',
     1           alpha1(icone),na1,na2,alpha2(icone),phist(icone),
     1           alfast(icone),tetan
c            if (raypatt.ne.'toray') then
            if (raypatt.eq.'genray') then
               alpha1(icone)=alpha1(icone)*trnspi
               alpha2(icone)=alpha2(icone)*trnspi
               tetan=0.5d0*pi-betast(icone) !Polar angle (radians)
               write(*,*)'in dinit_mr before cone_ec'
               !-----------
               call cone_ec(alpha1(icone),na1,na2,alpha2(icone),
     1              phist(icone),alfast(icone),tetan,powtot(icone),nray,
     1              alphaj(nray1),betaj(nray1),powj(nray1)) !->out
               !-----------
               do i=nray1,(nray1-1)+nray
                  zstj(i)=zst(icone)
                  rstj(i)=rst(icone)
                  phistj(i)=phist(icone)
                  xstj(i)=rstj(i)*dcos(phistj(i))
                  ystj(i)=rstj(i)*dsin(phistj(i))
               enddo
               write(*,*)'in dinit_mr after cone_ec:'
               do i=nray1,(nray1-1)+nray
                  write(*,*)' i,powj(i)',i,powj(i)
                  write(*,*)' i,betaj(i),alphaj(i)',i,betaj(i),alphaj(i)
               enddo
               !pause
            endif
            if(raypatt.eq.'toray') then
               alfast(icone)=alfast(icone)/trnspi
               tetan=betast(icone)/trnspi
               write(*,*)'dinit bef raypat: tetan,alfast,gzone,nray_in',
     1              tetan,alfast(icone),gzone,nray_in
               call raypat(tetan,alfast(icone),alpha1(icone),cr,nray_in,
     1              gzone,mray,betaj(nray1),alphaj(nray1))
               
               write(*,*)'dinit after raypat: gzone,nray_in',
     1              gzone,nray_in
               nray=nray_in
               do i=nray1,(nray1-1)+nray
                  zstj(i)=zst(icone)
                  rstj(i)=rst(icone)
                  phistj(i)=phist(icone)
                  xstj(i)=rstj(i)*dcos(phistj(i))
                  ystj(i)=rstj(i)*dsin(phistj(i))
               enddo
               
               do i=nray1,(nray1-1)+nray
                  write(*,*)'Raypat ray starting angles (degrees):'
                  write(*,*)' i,betaj(i),alphaj(i)',i,betaj(i),alphaj(i)
               enddo
               do i=nray1,(nray1-1)+nray
                  powj(i)=(powtot(icone)/nray)*1.e13
                  alphaj(i)=alphaj(i)*trnspi
                  betaj(i)=(90.-betaj(i))*trnspi
               enddo
               write(*,*)'in dinit_mr after raypat powj(i)'
               do i=nray1,(nray1-1)+nray
                  write(*,*)' i,powj(i)',i,powj(i)
                  write(*,*)' i,betaj(i),alphaj(i)',i,betaj(i),alphaj(i)
               enddo
            endif               !(On raypatt)
             
            if (raypatt.eq.'diskdisk') then
               call disk_to_disk_rays_initial_launching_data(nray)
            endif !diskdisk

            if (raypatt.eq.'diskbeam') then
               write(*,*)'dinit.f before'
               write(*,*)'disk_beam_rays_initial_launching_data'
                 call disk_beam_rays_initial_launching_data(nray)
               write(*,*)'dinit.f after'
               write(*,*)'disk_beam_rays_initial_launching_data'
            endif !diskbeam
            
            powtott=powtott+powtot(icone)*1.e13   !(erg/sec)
            nray1=nray1+nray
         enddo                  !(On icone)
         nray=nray1-1
        write(*,*)'dinit_mr_xyz: after EC start conditions, total nray',
     1        nray
      endif                     !(On istart.eq.1, EC)

      if (istart.eq.2) then
c---------LH or WF----------
          write(*,*)'dinit_mr_xyz: before grill_lh '
          write(*,*)'ngrilla,ngrill,rma,zma',ngrilla,ngrill,rma,zma
c--------------------------------------------------------------
          call grill_lh_xyz(rhopsi0,ngrill,thgrill,phigrill,
     1    height,nthin,
     1    anmin,anmax,nnkpar,powers,powtott,
     1    antormin,antormax,nnktor,anpolmin,anpolmax,nnkpol,
     1    n_theta_pol,
     1    rma,zma,psimag,
     1    fpwth,wtheta,wtheta_edge,anzin,anzin_edge,pwcpl,
     1    anztorin,anzpolin,pwcpl_tp,
     1    anztorin_edge,anzpolin_edge,pwcpl_t_ar,pwcpl_p_ar,
     1    nray,
     1    arxu0, aryu0, arzu0,  !!! YuP  output
     1    arnpar, arnper_tang, arntheta, arnphi, !!! YuP  output
     1    powinilh,
     1    wdnpar0,igrillpw,igrilltw,i_n_poloidal,
     1    ilaunch,r0launch,phi0launch,z0launch,
     1    i_grill_pol_mesh,i_grill_npar_ntor_npol_mesh)
             write(*,*)'dinit_mr_xyz: after grill nray=',nray
          do iray=1,nray         
             write(*,*)'iray,arxu0(iray),aryu0(iray),arzu0(iray)',
     &                  iray,arxu0(iray),aryu0(iray),arzu0(iray)
          enddo
      endif !LH or FW

      if (istart.eq.3) then
c---------ECR O_X_EBW mode conversion case----------
c         It uses the wave input data from the grill form
c         It sets i_n_poloidal=1
c	  It calculates: 
c	  1) the space point M_conv=(rhoconv,theta)=(zconv,rconv)
c            at which x_e=wpw2in=1. The value of the poloidal
c            angle theta (degree) is given bellow
c            if theta=0  point is on the outer part of the magnetic surface
c            if theta=90 point is on the top of the magnetic surface
c         2) the optimal value N_parallel_optimal=cnparopt
c            for O_X mode conversion
c         3) gives
c            anmin(1)=cnparopt-0.01d0 
c            anmax(1)=cnparopt+0.01d0 
c         4) call grill_lh to calculate the initial data
c         5) gives the coordinates of the initial O_X mode conversion poin 
c            arxu0(1)=xconv
c            aryu0(1)=yconv
c            arzu0(1)=zconv
c-------------------------------------------------------------
          i_n_poloidal=1
c         these lines are to create the start point inside the 
c         plasma for the given point xe=wpw2in
cSAP050510	  wpw2in=1.d0
          wpw2in=1.002d0 !<-yup  ! was:  1.002d0
c          theta=0.d0   !poloidal angle  (degree)
c          theta=-30.d0
          theta=thgrill(1)
          phi=0.d0 ! YuP Added. Needs work? To be specified in genray.in?
          rmn=rmin+1.d-5 ! Limits for searching O-X conversion point
          rmx=rmax  ! Specify limits for searching O-X conversion point

          write(*,*)'dinit_mr_xyz: before owconvr_xyz theta=',theta
          
          call owconvr_xyz(theta,phi,rmn,rmx,wpw2in,xconv,yconv,zconv)

     	    write(*,*)'dinit_mr_xyz xconv,yconv,zconv',xconv,yconv,zconv
          write(*,*)'dinit_mr_xyz old anmin,anmax(1)',anmin(1),anmax(1)

          bmod=       bxyz(xconv,yconv,zconv) !-> get bmod for wcw
          wcw_conv=    wcw(xconv,yconv,zconv,1)
          wpw2_conv= wpw_2(xconv,yconv,zconv,1) !-> rho is in one.i
          rhopsi0(1)=rho
          rconv= sqrt(xconv**2+yconv**2)
          
c---------calculation of the optimal value N_parallel_optimal
c         for O_X mode conversion
          cnparopt=dsqrt(wcw_conv/(1.d0+wcw_conv))
          anmin(1)=cnparopt-0.01d0 
          anmax(1)=cnparopt+0.01d0 
          write(*,*)'dinit new anmin(1),anmax(1)',anmin(1),anmax(1)
          write(*,*)'dinit new rhopsi0(1),Xe=',rhopsi0(1), wpw2_conv
c---------------------------------------------------------------
          write(*,*)'dinit_mr_xyz before grill_lh ngrilla,ngrill',
     1    ngrilla,ngrill
     
          arxu0(1)=xconv
          aryu0(1)=yconv
          arru0(1)=rconv
          arzu0(1)=zconv

          call grill_lh_xyz(rhopsi0,ngrill,thgrill,phigrill,
     1    height,nthin,
     1    anmin,anmax,nnkpar,powers,powtott,
     1    antormin,antormax,nnktor,anpolmin,anpolmax,nnkpol,
     1    n_theta_pol,
     1    rma,zma,psimag,
     1    fpwth,wtheta,wtheta_edge,anzin,anzin_edge,pwcpl,
     1    anztorin,anzpolin,pwcpl_tp,
     1    anztorin_edge,anzpolin_edge,pwcpl_t_ar,pwcpl_p_ar,
     1    nray,
     1    arxu0, aryu0, arzu0,  !!! YuP  output
     1    arnpar, arnper_tang, arntheta, arnphi, !!! YuP  output
     1    powinilh,
     1    wdnpar0,igrillpw,igrilltw,i_n_poloidal,
     1    ilaunch,r0launch,phi0launch,z0launch,
     1    i_grill_pol_mesh,i_grill_npar_ntor_npol_mesh)

          write(*,*)'dinit_mr_xyz aft grill_lh: arxu0,arzu0',
     +               arxu0(1),arzu0(1)
      endif ! O_X_EBW  mode conversion case
           
c-----allocate pointers at writencdf.i and write_i
      write(*,*)'dinit_mr_xyz before ainalloc_writencdf nray',nray
      call ainalloc_writencdf_i(nray)
      call ainalloc_write_i(nray)

      return
      end ! dinit_mr_xyz


c======================================================================
c======================================================================


      subroutine ox_conversion_grill_in_poloidal_point_xyz(theta_pol,
     & i_n_optimal)   !  output data are in 'grill.i'

c    It calculates: 
c	  1) the space point M_conv=(rhoconv,theta)=(zconv,rconv)
c            at which x_e=wpw2in=1. The value of the poloidal
c            angle theta (degree)=theta_pol is given bellow
c            if theta=0  point is on the outer part of the magnetic surface
c            if theta=90 point is on the top of the magnetic surface
c         2) the optimal value N_parallel_optimal=cnparopt
c            for O_X mode conversion
c            if i_n_optimal=1 cnparopt >0
c            if i_n_optimal=2 cnparopt <0
c            
c         3) gives
c            anmin(1)=cnparopt-0.01d0 
c            anmax(1)=cnparopt+0.01d0 
c         4) call grill_lh to calculate the initial data
c         5) gives the coordinates of the initial O_X mode conversion poin 
c            arxu0(1)=xconv
c            aryu0(1)=yconv
c            arzu0(1)=zconv
c 
c     output data are in 'grill.i'
c-------------------------------------
c     that gives the optimal N_parallel at X_e=1 surface 
     
      implicit integer (i-n), real*8 (a-h,o-z)

c-----input 
      real*8 theta_pol
      integer i_n_optimal
     
      include 'param.i'
      include 'one.i'
      include 'ions.i'
      include 'three.i'
      include 'five.i'
      include 'cone.i'
      include 'grill.i'
      include 'rkutta.i'
      include 'six.i'
      !!! include 'scatnper.i'
      include 'write.i'

      istart=3
      
c-----ECR O_X_EBW mode conversion case----------
c         It uses the wave input data from the grill form
c         It sets i_n_poloidal=1
c	  It calculates: 
c	  1) the space point M_conv=(rhoconv,theta)=(zconv,rconv)
c            at which x_e=wpw2in=1. The value of the poloidal
c            angle theta (degree) is given bellow
c            if theta=0  point is on the outer part of the magnetic surface
c            if theta=90 point is on the top of the magnetic surface
c         2) the optimal value N_parallel_optimal=cnparopt
c            for O_X mode conversion
c            if i_n_optimal=1 cnparopt >0
c            if i_n_optimal=2 cnparopt <0
c            
c         3) gives
c            anmin(1)=cnparopt-0.01d0 
c            anmax(1)=cnparopt+0.01d0 
c         4) call grill_lh to calculate the initial data
c         5) gives the coordinates of the initial O_X mode conversion poin 
c            arxu0(1)=xconv
c            aryu0(1)=yconv
c            arzu0(1)=zconv
c-------------------------------------------------------------
      i_n_poloidal=1
c     the following lines are to create the start point inside the 
c     plasma for the given point xe=wpw2in
      wpw2in=0.8d0  !<-yup  !was 0.999d0
      theta=theta_pol
      thgrill(1)=theta
      phi=0.d0 ! YuP Added. Needs work? To be specified in genray.in?
      rmn=1.d-5 ! Specify limits for searching O-X conversion point
      rmx=min(xeqmax,yeqmax,abs(xeqmin),abs(yeqmin)) ! Specify limits for searching O-X conversion point

     
      write(*,*)'ox_conversion before owconvr theta,Xein='
     &,theta,wpw2in
     
      call owconvr_xyz(theta,phi,rmn,rmx,wpw2in,xconv,yconv,zconv)

      write(*,*)'ox_conversion xconv,yconv,zconv',
     &           xconv,yconv,zconv

      bmod=       bxyz(xconv,yconv,zconv) !-> get bmod for wcw
      wcw_conv=    wcw(xconv,yconv,zconv,1)
      wpw2_conv= wpw_2(xconv,yconv,zconv,1) !-> rho is in one.i
      rhopsi0(1)=rho
      rconv= sqrt(xconv**2+yconv**2)
c-----calculation of the optimal value N_parallel_optimal
c     for O_X mode conversion
      cnparopt=dsqrt(wpw2_conv/(1.d0+wpw2_conv))
      if  (i_n_optimal.eq.2) cnparopt=-cnparopt
      anmin(1)=cnparopt-0.01d0 
      anmax(1)=cnparopt+0.01d0
      
      write(*,*)'ox_conversion rhopsi0(1),Xe=',rhopsi0(1),wpw2_conv
      write(*,*)'ox_conversion anmin(1),anmax(1)',anmin(1),anmax(1)

c---------------------------------------------------------------
      arxu0(1)=xconv
      aryu0(1)=yconv
      arzu0(1)=zconv
      arru0(1)=rconv

      call grill_lh_xyz(rhopsi0,ngrill,thgrill,phigrill,
     1height,nthin,
     1anmin,anmax,nnkpar,powers,powtott,
     1antormin,antormax,nnktor,anpolmin,anpolmax,nnkpol,
     1n_theta_pol,
     1rma,zma,psimag,
     1fpwth,wtheta,wtheta_edge,anzin,anzin_edge,pwcpl,
     1anztorin,anzpolin,pwcpl_tp,
     1anztorin_edge,anzpolin_edge,pwcpl_t_ar,pwcpl_p_ar,
     1 nray, arxu0, aryu0, arzu0,  !!!! YuP
     1 arnpar, arnper_tang, arntheta, arnphi, !!!! YuP
     1powinilh,
     1wdnpar0,igrillpw,igrilltw,i_n_poloidal,
     1ilaunch,r0launch,phi0launch,z0launch,
     1i_grill_pol_mesh,i_grill_npar_ntor_npol_mesh)

	write(*,*)'ox_conversion: arxu0,arzu0=', arxu0(1),arzu0(1)
      !pause
      
      return
      end
      

c======================================================================
c======================================================================


      subroutine gr_OX_optimal_direction_xyz(ndim)
c-----calculate the OX optimal direction of the EC cone central ray
      
      implicit integer (i-n), real*8 (a-h,o-z)
  
      include 'param.i'
      include 'one.i'
      include 'cone.i'
      include 'rkutta.i'
      include 'grill.i'
      include 'three.i' 
      dimension u(6),deru(6),aux(8,6)
      external rside_xyz, outpt_xyz, bxyz

c-----input integer ndim ! the number of ODE equations     

      write(*,*)'in gr_OX_optimal_direction_xyz ncone,i_ox',ncone,i_ox
      

      do icone=1,ncone

         icone_loc=icone ! for cone.i file
         
         if(i_ox_poloidal_max.gt.1)then
           dtheta= 
     +       (theta_top(icone)-theta_bot(icone))/(i_ox_poloidal_max-1)
         else
           dtheta= 0.d0
         endif
      
         do i_n_optimal=1,2 ! the loop for pos. and neg. N_parallel_optimal

           i_ox_poloidal=0 ! the counter of the poloidal rays for OX
                           ! optimal direction calculations

40         continue ! begin of the loop for the calculations 
                    ! of the optimal direction in the given EC cone vertex
                    ! for OX mode conversion.
                    ! i_ox_poloidal is incremented up to i_ox_poloidal_max
                    
c-----------------------------------------------------------------1beg
           if(i_ox.eq.1) then
             i_ox_poloidal=i_ox_poloidal+1 ! increment 
             write(*,*)'i_ox_poloidal=== ',i_ox_poloidal
             if(i_ox_poloidal.gt.i_ox_poloidal_max)then
                if(i_n_optimal.eq.1) then 
                   goto 20 ! Try Npar with other sign
                else
                write(*,*)'All angles between  theta_bot(icone) '
                write(*,*)'and theta_top(icone) are tried.'
                write(*,*)'If file ECcone_optimal.dat is not saved,'
                write(*,*)'none of rays came close to antenna.'
                write(*,*)'Try to increase i_ox_poloidal_max, '
                write(*,*)'or increase [theta_bot theta_top] range,'
                write(*,*)'or the vertical size of antenna eps_antenna.'
                stop 'i_ox_poloidal > i_ox_poloidal_max'
                endif
             endif

             theta_pol_0= theta_bot(icone) + (i_ox_poloidal-1)*dtheta
             if(i_ox_poloidal_max.eq.1)then
               theta_pol_0= 0.5*(theta_bot(icone) + theta_top(icone))
             endif
             !theta_bot and theta_top are given in genray.in file 
             !They determine the range of poloidal angles at Xe=1 surface
             !where the rays are started and traced back to the edge.
             !Those rays that come close to antenna (determined by 
             !rst(icone and zst(icone)) are selected for initial conditions
             !which are saved into ECcone_optimal.dat
             theta_pol=theta_pol_0              
c--------------------------------------------------------
c            for the given poloidal angle theta_pol
c            it will calculate
c            the coordinates of the initial O_X mode conversion point 
c            the optimal value N_parallel_optimal=(+,-) cnparopt
c            If i_n_optimal.eq.1 then  N_parallel_optimal=+cnparopt
c            If i_n_optimal.eq.2 then  N_parallel_optimal=-cnparopt
c            set the n_parallel grill spectrum
c            anmin(1)=cnparopt-0.01d0 
c            anmax(1)=cnparopt+0.01d0 
c--------------------------------------------------------
             write(*,*)'gr_OX_optimal_direction_xyz'
             write(*,*)'before ox_conversion_grill_in_poloidal_point'
             write(*,*)'theta_pol =',theta_pol
             call ox_conversion_grill_in_poloidal_point_xyz(theta_pol,
     &       i_n_optimal)
           endif !i_ox.eq.1
c-----------------------------------------------------------------1end

           iray=1

           write(*,*)'gr_OX_optimal_direction_xyz bef dinit_1ray'

           call dinit_1ray_xyz(arxu0(iray),aryu0(iray),arzu0(iray),
     1	   alfast1,betast1,arntheta(iray),arnphi(iray),u,iraystop)
          
           if(iraystop.eq.1) then
              write(*,*)'icone=',icone,'iray=',iray,'iraystop=1
     &                bad initial conditions'
              nrayelt=0
              goto 20
           endif 
  
           prmt(7)=prmt(6)  !here: gr_OX_optimal_direction_xyz
c----------------------------------------------------------
c          call bxyz() to calculate the small radius rho 
c          the result will be in common block  one.i
c-----------------------------------------------------------
cyup           bmod= bxyz(u(1),u(2),u(3)) !-> get rho_magnetic (no more)
           den=dense_xyz(u(1),u(2),u(3),1) !-> rho_magnetic or rho_dens

CMPIINSERTPOSITION STARTRUNGEKUTTA

           nstep_rk=1 ! initialize the number of Runge-Kutta time step 
            
c--------------------------------------------------------------beg2
           if(isolv.eq.1) then
c--------------------------------------------------------------
c            The Runge-Kutta solution of 6 hamiltonian equations
c            with correction which gives 
c            hamiltonian conservation with
c            accuracy epscor=prmt(4)
c--------------------------------------------------------------
             if(irkmeth.eq.0) then
c              4_th order Runge-Kutta method with constant time step
ccc               call drkgs(prmt,u,deru,ndim,ihlf,rside1,outpt,aux)
                 stop 'set irkmeth=2 or 3; other values not supported'
             endif
             if(irkmeth.eq.1) then
c              5-th order Runge-Kutta method with variable time step,
ccc               call drkgs1(prmt,u,deru,ndim,ihlf,rside1,outpt,aux)
                 stop 'set irkmeth=2 or 3; other values not supported'
             endif
             if(irkmeth.eq.2) then
c              4th order Runge-Kutta with variable time step,
c              time step can be reduce or enlarge
               write(*,*)'prmt',prmt
               write(*,*)'u===',u
               call drkgs2_xyz(prmt,u,deru,ndim,ihlf,
     +                         rside_xyz, outpt_xyz, aux, i_output)
               write(*,*)'u===',u
               write(*,*)'deru=',deru
             endif
             if(irkmeth.eq.3) then
c              4th order Runge-Kutta with variable time step,
c              automatically set, based on dL_step and dN_step
               call drkgs_auto(prmt, u, deru, ndim, ihlf,
     &                      rside_xyz, outpt_xyz, aux, dL_step,dN_step)
             endif
           endif ! isolv.eq.1
c-------------------------------------------------------------------end2


           if(isolv.eq.2) then
c-------------------------------------------------------------------
c            The Runge-Kutta solution of 5 hamiltonian equations
c            and the solution of the dispersion relation for the sixth
c            ray variable'
c-------------------------------------------------------------------------
ccc              call rkb1(prmt,u,deru,ndim,ihlf,rsideb1,outptb1,aux)
c------------------------------------------------------------------------
                 stop 'set isolv=1; other values not supported'
           endif ! isov.eq.2

c---------------------------------------------------------------------beg5
           if(i_ox.eq.1) then
c------------set the boundaries for the calculations of
c            OX optimal direction in EC cone vertex           
             write(*,*)'i_n_optimal===',i_n_optimal
             write(*,*)'i_ox_poloidal=',i_ox_poloidal
c---------------calculate the poloidal angle theta_st of the given
c               cone vertex
                rho_st=dsqrt((rst(icone)-rma)**2+(zst(icone)-zma)**2)
                cos_theta_st=(rst(icone)-rma)/rho_st
                sin_theta_st=(zst(icone)-zma)/rho_st
                if (cos_theta_st .gt.  1.d0)  cos_theta_st= 1.d0
                if (cos_theta_st .lt. -1.d0)  cos_theta_st=-1.d0
                if (sin_theta_st.ge.0.d0) then
                   theta_st=+dacos(cos_theta_st)
                else  
                   theta_st=-dacos(cos_theta_st)
                endif
               
                write(*,*)'***************************************'
                write(*,*)'rst(icone),zst(icone)',rst(icone),zst(icone)
                write(*,*)'rma,zma',rma,zma
                write(*,*)'rho_st',rho_st
                write(*,*)'rst(icone)-rma,cos_theta_st',
     &                    rst(icone)-rma,cos_theta_st
                write(*,*)'zst(icone)-zma,sin_theta_st',
     &                     zst(icone)-zma,sin_theta_st
                write(*,*)'theta_st[degree]',theta_st*180.d0/pi
                write(*,*)'r_st_ox,z_st_ox=', r_st_ox,z_st_ox
                write(*,*)'***************************************'

                rho_st_ox=dsqrt((r_st_ox-rma)**2+(z_st_ox-zma)**2)
                cos_theta=(r_st_ox-rma)/rho_st_ox
                sin_theta=(z_st_ox-zma)/rho_st_ox

                if (cos_theta .gt.  1.d0)  cos_theta= 1.d0
                if (cos_theta .lt. -1.d0)  cos_theta=-1.d0
                if (sin_theta.ge.0.d0) then
                  theta_st_ox=+dacos(cos_theta)
                else  
                  theta_st_ox=-dacos(cos_theta)
                endif
                                                
                write(*,*)'z_st_ox, zst(icone), eps_antenna=',
     &          z_st_ox, zst(icone), eps_antenna
                !pause

                if(abs(z_st_ox-zst(icone)) .le. eps_antenna*0.5)then
                   theta_ox=theta_pol_0
                   call write_optimal_O_cone_central_ray_multivertex(
     &                                              icone,i_n_optimal)
                   write(*,*)'ray is close to antenna;'
                   write(*,*)'data is saved into ECcone_optimal.dat'
                   goto 20 !-> next i_n_optimal (other Npar sign)
                else
                   goto 40 
                   !-> Try next theta 
                   !   in the the range [theta_bot theta_top]
                endif

           endif !i_ox.eq.1
c---------------------------------------------------------------------end5

 20      continue  
         enddo ! i_n_optimal=1,2 
      enddo   ! icone 
      
      stop 'oxb'                    
      return
      end


c======================================================================
c======================================================================

c        ********************hamiltmuz_ful_xyz ***************
c        * 						      *
c        * this subroutine  calculates: hamiltmz -real part of*
c        *                              Hamiltonian, and      *
c        * reps(3,3) -components of complex dielectric tensor *
c        *            for Mazzucato code                      *
c        * for complex n_perp
c        *****************************************************
c-------------------------------------------------------------
c       input parameters				     !
c       cnpar - parallel (to magnetic field) component	     !
c               of the refractive index			     !
c       cnper - real part of the perpendicular component of
c                the refractive index
c       cnperim -Imaginary part N_perp                        !
c       ihermloc=1 Hermitian dielectric tensor       	     !
c               =2 full      dielectric tensor               !
c       Attention: iherm is in common one.i		     !
c     	x,y,z  cartesian coordinates  of ray point		     !
c       .....						     !
c       rho-small radius was calculated in subroutine: b     !
c            rho is in common one    			     !
c       mode = +1., O-mode   is in common one		     !
c              -1., X-mode				     !
c       output:				     !
c       hamiltmz-Hamiltonian				     !
c       reps(3,3)-complex dielectric tensor in common /eps/  !
c------------------------------------------------------------
      subroutine hamiltmuz_ful_xyz(cnpar,ihermloc,x,y,z,cnper,cnperim,
     .  hamiltmz) !->out (also reps(3,3) in erps.i)
      implicit integer (i-n), real*8 (a-h,o-z)
      double precision mode,mu,npar,np2
      double precision ncld(2)
      include 'param.i'
      include 'one.i'
      dimension an2(2),icuto(2)
      double complex sol,cpp
      double complex chamilmz
      !----------------
      mode=dfloat(ioxm)
c	 write(*,*)'in hamilmuz mode=',mode
c         write(*,*)' hamiltmuz_ful cnpar,ihermloc,x,y,z,cnper,cnprim',
c     &   cnpar,ihermloc,x,y,z,cnper,cnprim
       xe=wpw_2(x,y,z,1)
       bmod=bxyz(x,y,z) !-> get b (needed for wcw)
       ye=wcw(x,y,z,1)
       te=tempe_xyz(x,y,z,1)
       mu=512.44d00/te
c	 write(*,*)'in hamiltmuz_ful xe,ye,te',xe,ye,te
         if(dabs(cnpar).lt.1.d-4) then
           if (cnpar.ge.0.d0) then
             sign=1.d0
           else 
             sign=-1.d0
           endif
           cnpar=sign*1.d-4
         endif
      npar=cnpar
c	 write(*,*)'in hamilmuz cnpar,cnper,ihermloc',
c     *   cnpar,cnper,ihermloc
c	 cpp=dcmplx(cnper,0.0d0)
      cpp=dcmplx(cnper,cnperim)
      np2=cnper*cnper
      xe=wpw_2(x,y,z,1)
      ye=wcw(x,y,z,1)
      sol=cpp*cpp
c         write(*,*)'hamiltmuz_ful cnper,cnperim,cpp',
c     &   cnper,cnperim,cpp
      call complx1(mode,xe,ye,mu,npar,np2,sol,ihermloc,chamilmz)
      hamiltmz=dreal(chamilmz)
      return
      end
c---------------------------------------------------------------------
      


c======================================================================
c======================================================================
      real*8 function rho_dens_xyz(x,y,z)
c===> Yu.P. 2011
c     Define rho for density, temperature, Tpop, Vflow, Zeff profiles.
c     Also define derivatives drho/dx, drho/dy, drho/dz !store in one.i
! -----------------------------------------------------
! Yu.P. Added in 2011                    model_rho_dens
! If model_rho_dens=0 or 5, the definition of rho for plasma profiles 
! is based on psi - magnetic flux, as before.
!.................
!.................
! model_rho_dens=1 or 2.   
!      Only works with idens=0.
!      The profile of density (also Temperature, Tpop, Vflow, Zeff)
!      is determined on analytically defined surfaces, 
!      different from magnetic flux surfaces.  
!      The models are adequate for a weakly-ionized FRC.
!.................
!.................
! model_rho_dens=1 
! Defines rho-coordinate on ellipsoid surfaces.
! See functions rho_dens_xyz(x,y,z) and dense_xyz(x,y,z,i) for details.
! rho is defined as:
!      rho2=   (( (x-elx0)*costt + (y-ely0)*sintt )/elax)**2
!             +((-(x-elx0)*sintt + (y-ely0)*costt )/elay)**2 
!             +((z-elz0)/elaz)**2  
!      rho_dens_xyz= sqrt(rho2)
! where
!      elthet=eltheta*pi/180 ! [rad] inclination of the ellipse in x-y
!      sintt= sin(elthet) !-> to one.i
!      costt= cos(elthet) !-> to one.i
! where elax, elay, elaz are the semi-axes of the ellipsoid 
! in cartesian coordinates, 
! and elx0, ely0, elz0 are the coordinates of its center.
! Axis z is in the same direction as the vertical coordinate
! in tokamaks; x and y are in r-phi plane.
! The profile of density is then set as (idens=0 option)
!   dens1(k,i)=(dense0(i)-denseb(i))*
!              (1-rhom(k)**rn1de(i))**rn2de(i) + denseb(i)
! Same definition of rho is applied to setting profiles of
! temp1, tpop1, vflow1, zeff1.
!.................
!.................
! model_rho_dens=2 
! The profile of density is set as the sum of three types:
! Rigid Rotor ("rr") profile, 
! Ellipsoidal Spindle ("es") profile,
! Uniform Background ("ub") profile.
! See functions rho_dens_xyz(x,y,z) and dense_xyz(x,y,z,i) for details.
! rho is defined similar to model_rho_dens=1:
!      rho2=   (( (x-elx0)*costt + (y-ely0)*sintt )/elax)**2
!             +((-(x-elx0)*sintt + (y-ely0)*costt )/elay)**2 
! (uniform in z, for now)
!      rho_dens_xyz= sqrt(rho2)
! where
!      elthet=eltheta*pi/180.d0 ! [rad] inclination of the ellipse in x-y
!      sintt= sin(elthet) !-> to one.i
!      costt= cos(elthet) !-> to one.i
! where elax, elay are the semi-axes of the ellipse 
! and   elx0, ely0 are the coordinates of its center.
!-1-> Rigid Rotor profile:
!       Rs0rr= max(elax,elay)
!       rrk=(Rm0rr/Rs0rr)**2  ! K in Eq.(38)
!       Rs=Rs0rr ! Assume no dependence in z, for now
!       dens_rr= dens0rr*( sech(rho2-rrk)/sech(rrk) )**2
!       ! Note: At rho=0 -> nrr= n0rr
!       ! Note: At rho=1 -> nrr= n0rr*[sech(1-K)/sech(K)]^2
!-2-> Ellipsoidal Spindle profile:
!       dens_es= dens0es* exp( -(Rs*rho-Rm0es)**2 / (2*rtau**2) )
!-3-> Uniform background density profile:
!       dens_ub= dens0ub 
!       r= sqrt(x*x+y*y)
!       if (r .gt. r_ub_edge) then 
!          ! Linear drop in region  r_ub_edge < r < wall_rmax
!          dens_ub= dens0ub*(1.d0 -(r-r_ub_edge)/(wall_rmax-r_ub_edge))
!       endif
!===> NEED TO SPECIFY:
! Peak densities for the three types of profiles:
! dens0rr= ! [m^-3] Peak density for Rigid Rotor profile
! dens0es= ! [m^-3] Peak density for Ellipsoidal Spindle profile
! dens0ub= ! [m^-3] Density for Uniform Background profile
!-1-> Input for the Rigid Rotor profile:
! Position of plasma center:
! elx0= ! [m]
! ely0= ! [m]
! Semi-axes of the ellipse:
! elax= ! [m] 
! elay= ! [m]
! eltheta= ! [degrees] inclination of the ellipse-profile in x-y plane
! (Note: Rs0rr (Separatrix radius at z=z0) is set to max(elax,elay) )
! Rm0rr= ! [m] Radius of peak power emission at z=z0, for Rigid Rotor
!-2-> Input for the Ellipsoidal Spindle profile:
! Rm0es= ! [m] Radius of peak power emission at z=z0, for Ell.Spindle
! rtau=  ! [m] Radial decay distance (==tau_R in Eq.(44))
!-3-> Input for the uniform background profile:
! r_ub_edge= ! [m] Radius for the uniform background density; 
             !     drops linearly to zero from r_ub_edge to wall_rmax
!.................
! (Other models can be added later).   
!------------------------------------------------------

      implicit none
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
c-----Input:
      real*8 x,y,z
c-----Externals:
      real*8 psif_xyz, rhopsi, drhodr, drhodz
c     Local:
      real*8 psi,rho2,r,r2,rs_frc2,phi,drho_dr
      !------------      
      rho_dens_xyz= 1.d0 ! to initialize
      
      goto (10,11,12,13,14) model_rho_dens+1 ! 

  10  continue ! model_rho_dens=0 or model_rho_dens=5
        ! rho definition is based on magnetic flux psi
        r2= x*x+y*y
        r=  dsqrt(r2)
        psi= psif_xyz(x,y,z) ! spline 
        ! YuP[2016] Spline may give a value of psi 
        ! a little smaller than psimag (at r=0) => 
        ! may cause error in rhopsi (fixed).
        rho_dens_xyz= rhopsi(psi)  ! spline
        
        phi=0.d0 ! not important
        drho_dr= drhodr(z,r,phi) ! not actually using phi
        drho_dz= drhodz(z,r,phi) ! not actually using phi

        if(r.gt.0.d0) then
           drho_dx= drho_dr*(x/r)
           drho_dy= drho_dr*(y/r)
        else
           drho_dx= 0.d0
           drho_dy= 0.d0
        endif
      return
      
  11  continue ! model_rho_dens=1
        ! Ellipsoid surface for n,T,Zeff,...
        rho2=   (( (x-elx0)*costt + (y-ely0)*sintt )/elax)**2
     +         +((-(x-elx0)*sintt + (y-ely0)*costt )/elay)**2 
     +         +((z-elz0)/elaz)**2  
        rho_dens_xyz= sqrt(rho2)
        drho_dx= (   
     +    (( (x-elx0)*costt + (y-ely0)*sintt )/elax)*( costt/elax)
     +   +((-(x-elx0)*sintt + (y-ely0)*costt )/elay)*(-sintt/elay)
     +            )/rho_dens_xyz
        drho_dy= (   
     +    (( (x-elx0)*costt + (y-ely0)*sintt )/elax)*( sintt/elax)
     +   +((-(x-elx0)*sintt + (y-ely0)*costt )/elay)*( costt/elay)
     +            )/rho_dens_xyz
        drho_dz= ((z-elz0)/elaz**2)/rho_dens_xyz
      return


  12  continue ! model_rho_dens=2
        rho2=   (( (x-elx0)*costt + (y-ely0)*sintt )/elax)**2
     +         +((-(x-elx0)*sintt + (y-ely0)*costt )/elay)**2 
               ! uniform in z, for now
        rho_dens_xyz= sqrt(rho2)
        drho_dx= (   
     +    (( (x-elx0)*costt + (y-ely0)*sintt )/elax)*( costt/elax)
     +   +((-(x-elx0)*sintt + (y-ely0)*costt )/elay)*(-sintt/elay)
     +            )/rho_dens_xyz
        drho_dy= (   
     +    (( (x-elx0)*costt + (y-ely0)*sintt )/elax)*( sintt/elax)
     +   +((-(x-elx0)*sintt + (y-ely0)*costt )/elay)*( costt/elay)
     +            )/rho_dens_xyz
        drho_dz= 0.d0
      return
      
  13  continue ! model_rho_dens=3 
        stop 'rho_dens_xyz() is not setup for model_rho_dens=3'
        ! Note: For model_rho_dens=3, rho is defined in dense_xyz
        ! The derivatives of rho are not ready yet.
      return
      
  14  continue ! model_rho_dens=4
        r2= x*x+y*y
        r=  dsqrt(r2)
        rs_frc2= rs_frc**2 
        rho_dens_xyz=abs(2.*r2/rs_frc2 -1.d0) ! rs_frc is in one.i
        ! With such definition, rho=0 at r=rs_frs/sqrt(2), 
        ! and rho=1 at r=0 or r=rs_frc .
        if(2.*r2>rs_frc2)then ! r>r_null
           drho_dx= 4.*x/rs_frc2
           drho_dy= 4.*y/rs_frc2
        else
           drho_dx=-4.*x/rs_frc2
           drho_dy=-4.*y/rs_frc2
        endif
        drho_dz= 0.d0
      return

      ! other values of model_rho_dens can be added later.
      PRINT*,'Set model_rho_dens to 0, 1, 2, 3, or 4 or 5'
      PRINT*,'model_rho_dens=1-5 works with ixyz=1 only (cartesian)'
      PRINT*,'model_rho_dens=0 works with ixyz=0 or 1'
      stop 'rho_dens_xyz: can be called for model_rho_dens=0-5'
      end ! rho_dens_xyz(x,y,z)
      
      
c======================================================================
c======================================================================

      subroutine rho_theta_phi_xyz(rhon,theta,phi,
     +                             x,y,z) !-> out
c     Yu.P. 2011
c     Calculate x,y,z, cartesian coords. 
c     from known rho_density, theta_poloidal, phi_toroidal
      implicit none
      include 'param.i'
      include 'one.i'
c-----Input:
      real*8 rhon,theta,phi
c-----Output:
      real*8 x,y,z
c-----Local
      real*8 cos_phi,sin_phi,cos_theta,sin_theta
! model_rho_dens=1 defines rho_dens as 
! sqrt( (x-elx0)^2/elax^2 + (y-ely0)^2/elay^2 + (z-elz0)^2/elaz^2 )
! where elax, elay, elaz are the semi-axes of the ellipsoid 
! in cartesian coordinates, 
! and elx0, ely0, elz0 are the coordinates of its center.
      
      if(model_rho_dens.eq.1) then ! ellipsoid surface
         cos_phi=dcos(phi)
         sin_phi=dsin(phi)
         cos_theta=dcos(theta)
         sin_theta=dsin(theta)
         x= elx0 +rhon*cos_theta*(elax*cos_phi*costt-elay*sin_phi*sintt)
         y= ely0 +rhon*cos_theta*(elax*cos_phi*sintt+elay*sin_phi*costt)
         z= elz0 +rhon*sin_theta*elaz
         return
      endif
      
      if(model_rho_dens.eq.2) then ! cylinder surface at ellipse in X-Y
         cos_phi=dcos(phi)
         sin_phi=dsin(phi)
         x= elx0 +rhon*(elax*cos_phi*costt-elay*sin_phi*sintt)
         y= ely0 +rhon*(elax*cos_phi*sintt+elay*sin_phi*costt)
         z= sqrt(x*x+y*y)*dtan(theta)
         return
      endif
            
      stop 'rho_theta_phi_xyz: Only works for model_rho_dens=1 or 2'
      return
      end      

c======================================================================
c======================================================================

c*************************rbmax_xyz **********************************
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
c------------------------------------------------------------------
c     output: arrays ar_min(npsi),ar_max(npsi),ab_max(npsi) in gr.i
c------------------------------------------------------------------
      subroutine rbmax_xyz ! for y=0 only ! (needed for graphics?)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'gr.i'

      y=0.d0
      do j=1,npsi
        x=rpsi(j,1)  ! x=r for y=0, for now
        z=zpsi(j,1)
        bmod=bxyz(x,y,z) !-> get bmod
        ar_max(j)=x
        ar_min(j)=x
        ab_max(j)=bmod
        ab_min(j)=bmod
        do i=1,nteta
           x=rpsi(j,i)  ! x=r for y=0, for now
           z=zpsi(j,i)
           bmod=bxyz(x,y,z) !-> get bmod
           if (ar_min(j).gt.x) ar_min(j)=x
           if (ar_max(j).lt.x) ar_max(j)=x
           if (ab_max(j).lt.bmod) ab_max(j)=bmod
           if (ab_min(j).gt.bmod) ab_min(j)=bmod
        enddo ! i=1,nteta
      enddo ! j=1,npsi
c     arrays ar_min(npsi),ar_min(npsi),ab_max(npsi),ab_min(npsi)
c     have been created. They are in gr.i
      return
      end



c======================================================================
c======================================================================
c-----------------------------------------------------------------
c     this  subroutine calculates current drive efficiency
c     for the toroidal  plasma using CURBA code
c     ATTENTION:all parameters and variables input to curba
c       are real(no double precision)
c-----------------------------------------------------------------
      subroutine effcurb_xyz(x,y,z,zma,rma,r0x,z_eff,temp,den,
     1                   jwave,cnpar,ye,
     +                   efficient) !->out
c     input: x,y,z - cartesian coordinates of the ray point(cm) !!!
c                       rma,zma -coordinares of magnetic axis(cm) !!!
c                       r0x character length (m)
c                       z_eff-effective charge
c                       temp-temperature (keV)
c                       den-density(10**13 cm**-3)
c                       jwave -wave harmonic number
c                       cnpar -parallel to magnetic field refractive
c                              index
c                     	ye-omega_Be/omega_wave
c     output: J/P=efficient  in (A/m**2)/(Wt/m**3)
c------------------------------------------------------------------
c     these double precision parameters are to call function: 
c        psif_xyz 
c-------------------------------------------------------------------
      implicit none
c-----input
      real*8 x,y,z,r,zma,rma,z_eff,temp,den,cnpar,ye
      real*8 r0x
      integer jwave
c-----output
      real*8 efficient

c-----locals
      real*8 psid,zd,xd,yd,rd,rmaxpsid,rminpsid,pid,
     &zbmin,zbmax,zb
      real*8 tc,tol,elomom,prho,ctheta,stheta,pi,theta,thtc,enpar,denom,
     *aspct,zacosarg,rjpd,rjpd0,ratjpd,efficien
      integer ig,n0,model,lh
c-----external 
      real*8 bmin_psi,bmax_psi,rmax_psi,rmin_psi,bxyz,zfac_f_xyz,
     + psif_xyz

      
c-------------------------------------------------------------------
c     for curba_GA
      real*8 rjpd_d,rjpd0_d,ratjpd_d,denom_d,aspct_d,enpar_d, 
     &tc_d, thtc_d, theta_d, elomom_d,zeff_d,tol_d
c-----external
c-----for efficiency transformation like in currn in toray
      real*8 zfac
c-------------------------------------------------------------------

c      write(*,*)'effcurb_xyz: x,y,z, rma,zma',x,y,z,rma,zma

c     ig=+1 selects new	Green's fcn in curba.-1 selects old
c        +2 or +3 are older models
c---------------------------------------------------------------------
      ig=1
c--------------------------------------------------------------------
c     n0 is number of mesh points for gaussian integration of Green's
c        function; n0=2,4,6,8,10,12,16,20,24,32,48 or 64. n0=128 gives
c        64 points over resonanse if resonanse has only single passing
c        particle segment; 64 points on each if two segments.20,24,32,
c        48 or 64.
c---------------------------------------------------------------------
c      n0=4
       n0=32
c      n0=20

c--------------------------------------------------------------------
c    tc is a bulk electron temperature in keV
c--------------------------------------------------------------------
      tc=temp
      
c--------------------------------------------------------------------
c     tol is absolute toleranse for D02GBF integrator; starting in
c       version	1.1 the variables integrated by D02GBF are normalized
c       to be O(1), so tol can be viewed as relative error tolerense.
c--------------------------------------------------------------------
c      tol=2.e-3
      tol=1.d-2
c--------------------------------------------------------------------
c     model Absolute value of model selectes collisional model: 1 for
c     full bounce av, 2 for square well(numerical solution), 3 for
c     analytic solution to square well.negative model does parallel
c     heating (lower hybrid, fast wave)
c-------------------------------------------------------------------
      model=3
c-------------------------------------------------------------------
c     z_eff is ion effective charge
c-------------------------------------------------------------------
c     lh ABS(lh) is power of e_perp in diffusion coeff; sign governs
c       power of p_parallel in diffusion coeff: + gives p_par^0;
c       - gives p_par^2. For ECRH, |lh|	is harmonic number (lh=1
c       for fundamental);+ for E_- contribution, - for E_parallel;
c        + with lh-->lh+2 for E_+ component of electric field E.
c       For now ,there is noprovision for the p_par^2 option with lh=0,
c       as the compiler doesn't know the differaence between +0 and -0.
c       If there is any interest, sent a message to use 313 and
c       revision including this option will be created.
c--------------------------------------------------------------------
      lh=jwave
c      lh=2
c      lh=0
c--------------------------------------------------------------------
c     elomom is ABS (lh)*cyclotron frec over wave frec for Y<0,
c       interpret as evaluated at poloidal angle where electrons are
c       first trapped in bucket.
c--------------------------------------------------------------------
      elomom=lh*dabs(ye)
c--------------------------------------------------------------------
c     theta is a poloidal angle at which power is absorbed (measured
c       from outside); for thtc<0,theta is poloidal angle at which
c       electrons become trapped in bucket; used to calculate resonant
c       energy given elomom and enpar. Note if calling programm fixes
c       minimum resonant energy and calculates enpar, then theta has no
c       significance to the physics; rjpd depends only on
c       B_min*Y/B(theta) ,not on theta or Y alone. But it is useful to
c       be able to speciffy theta and Y separately in order to keep
c       track of what higher harmonic are doing.
c--------------------------------------------------------------------
      r=dsqrt(x*x+y*y)
      prho=dsqrt((z-zma)**2+(r-rma)**2)
      ctheta=(r-rma)/prho
      stheta=(z-zma)/prho
      pi=4.d0*datan(1.d0)
c      write(*,*)' effcurba ctheta_geom,stheta',ctheta,stheta
      if(stheta.ge.0.0d0) then
         theta=dacos(ctheta)
      else
         theta=2.0d0*pi-dacos(ctheta)
      end if
c      write(*,*)' effcurba acos(ctheta),theta',acos(ctheta),theta
c--------------------------------------------------------------------
c     thtc is ratio of temperature along characteristic to that of
c       bulk, except for thtc <0, -thtc is energy (in keV) of bucket
c       rise;
c--------------------------------------------------------------------
      thtc=1.0d0
c--------------------------------------------------------------------
c     enpar is k_parallel*c/wave frec. Note enpar**2 <1-elomom**2
c       implies no resonant particles; rjpd set to zero and rjpd0
c       determined for vparallel given by nonrelativistic resonance
c       condition
c--------------------------------------------------------------------
      enpar=cnpar
c--------------------------------------------------------------------
c     aspct is inverse aspct ratio
c--------------------------------------------------------------------
c bmax,bmin  in toray
c     if(dabs(bmax-bmin).lt.1.0d-8) bmax=bmin+1.0d-8
c 	aspct=(bmax-bmin)/(bmax+bmin)
c--------------------------------------------------------------------
c     Here :aspct is calculaed using rmax_psi and rmin_psi
c  conversion from real to double precision
      pid=4.d0*datan(1.d0)
c conversion from non-dimensional to m
      xd=dble(x)*0.01d0/r0x
      yd=dble(y)*0.01d0/r0x
      rd=dble(r)*0.01d0/r0x
      zd=dble(z)*0.01d0/r0x
      psid=psif_xyz(xd,yd,zd)
       zbmin=bmin_psi(psid)    
       zbmax=bmax_psi(psid)
       aspct=((zbmax-zbmin)/(zbmax+zbmin))
cSm060201
       if (aspct.lt.1.d-8) aspct=1.d-8 

c     function bxyz() uses the arguments x,y,z in [m]
      zb=bxyz(.01d0*x/r0x, .01d0*y/r0x, .01d0*z/r0x)!it changes 
                                         !bx,by,bz in one.i
      zacosarg=(zbmax + zbmin - 2.d0*zbmax*zbmin/zb)/
     &         (zbmin - zbmax)

c      write(*,*)'zbmax,zb,zbmin,zacosarg',zbmax,zb,zbmin,zacosarg
      if (dABS (zacosarg) .gt. 1.0d0) zacosarg = 
     1                               DSIGN (1.0d0, dble(zacosarg))
      theta =dACOS (zacosarg) 
      denom=1.d0
c      write(*,*)'effcurb befre curba denom,aspct,enpar,tc,thtc',
c     &denom,aspct,enpar,tc,thtc
c      write(*,*) 'theta,elomom,lh,z_eff,model,tol,n0,ig',
c     &            theta,elomom,lh,z_eff,model,tol,n0,ig

c-----using  TorGAcurba with real*8 arguments
      denom_d= denom
      aspct_d=aspct
      enpar_d=enpar 
      tc_d=tc
      thtc_d=thtc
      theta_d=theta
      elomom_d=elomom
      zeff_d=z_eff
      tol_d=tol
cSAP080617
c      write(*,*)'before TorGA_curba '

      call TorGA_curba (rjpd_d,rjpd0_d,ratjpd_d,denom_d,aspct_d,enpar_d, 
     &tc_d, thtc_d, theta_d, elomom_d, lh, zeff_d, model, tol_d, n0, ig)
      rjpd=rjpd_d
      rjpd0=rjpd0_d
      ratjpd=ratjpd_d
c      write(*,*)'after TorGA_curba rjpd,rjpd0,ratjpd,denom_d',
c     &rjpd,rjpd0,ratjpd,denom_d

c----------------------------------------------------------
c     efficiency  in (A/cm**2)/(erg/(sec*cm**3))
c     temperature temp in kev
c     density     den in 10**19/m**3
c-----------------------------------------------------------
c-----genray transformation of CD efficiency 
c     cln -Coulomb Ln
c     it is necessary to write the formula to calculate cln

c      cln=17.0
c      arg1=1.d3/temp*dsqrt(10.d0*den) !temp KeV, den 10**13/cm**3
c      cln=24.d0-alog(arg1)
       efficien=rjpd 
c      efficient=efficien*temp/den*(17.0/cln)*4.5*1.e-6
c      write(*,*)'efficiecy coefficient=temp/den*(17.0/cln)*4.5*1.e-6',
c     1temp/den*(17.0/cln)*4.5*1.e-6
c      write(*,*)'in effic rjpd,rjpd0,ratjpd,denom',
c     1           rjpd,rjpd0,ratjpd,denom
c      write(*,*)'end effic curba efficient',efficient

c-----the efficiency transformation from torayGA
c      write(*,*)'before zfac_f temp',temp
      zfac=zfac_f_xyz(x*1.d-2,y*1.d-2,z*1.d-2,dble(temp)) 
      zfac=zfac*1.d-7 ! watt => egr/sec
      zfac=zfac*1.d-4 ! m**2 => cm**2
c      write(*,*)'zfac',zfac
      efficient=efficien*zfac
cSm050923
      efficient=-efficient 

      return
      end

c======================================================================
c======================================================================

c     this  subroutine calculates current drive efficiency
c     for the toroidal  plasma using TorGA_curgap codo written by Lin-Liu
c     ATTENTION:all parameters and variables input to curba
c       are real(no double precision)
c-----------------------------------------------------------------
      subroutine eff_Lin_Liu_xyz(x,y,z,zma,rma,r0x,z_eff,temp,den,jwave,
     &                   cnpar,cnper,ioxm,ye,
     &                   cefldx,cefldy,cefldz,
     &                   efficient) !-> out

c     input parameters: x,y,z -cartesian coordinates of the ray point(cm)

c                       rma,zma -coordinares of magnetic axis(cm)
c                       r0x character length (m)
c                       z_eff-effective charge
c                       temp-temperature (keV)
c                       den-density(10**13 cm**-3)
c                       jwave -wave harmonic number
c                       cnpar -parallel to magnetic field refractive
c                              index
c                       cnper Re(N_perp)    
c                       ioxm= +1 O mode,
c                             -1 X mode
c                     	ye-omega_Be/omega_wave
c     cefldx = the x-component of the wave electric field (COMPLEX)
c     cefldy = the y-component                            (COMPLEX)
c     cefldy = the z-component                            (COMPLEX)
c     output parameter: J/P=efficient  in (A/m**2)/(Wt/m**3)
c------------------------------------------------------------------
c     these double precision parameters are to call function: 
c        psif_xyz 
c-------------------------------------------------------------------
      implicit none
c-----input
      real*8 x,y,z,r,zma,rma,z_eff,temp,den,cnpar,ye,cnper
      real*8 r0x
      complex*16 cefldx,cefldy,cefldz
      integer jwave,ioxm
c-----output
      real*8 efficient

c-----locals
      real*8 psid,zd,xd,yd,rd,rmaxpsid,rminpsid,pid,
     &zbmin,zbmax,zb
      real*8 tc,tol,elomom,prho,ctheta,stheta,pi,theta,thtc,enpar,denom,
     *aspct,zacosarg,rjpd,rjpd0,ratjpd,efficien
      integer ig,n0,model,lh
c-----external 
      real*8 bmin_psi,bmax_psi,rmax_psi,rmin_psi,bxyz,zfac_f_xyz,
     + psif_xyz

      
c-------------------------------------------------------------------
c     for TorGA_curgap
      real*8 rjpd_d,rjpd0_d,ratjpd_d,denom_d,aspct_d,enpar_d, 
     &tc_d, thtc_d, theta_d, elomom_d,zeff_d,tol_d, omode
c-----external
c-----for efficiency transformation like in currn in toray
      real*8 zfac
c-------------------------------------------------------------------
c      write(*,*)'eff_lin_liu z,r,zma,rma,r0x,z_eff',
c     &z,r,zma,rma,r0x,z_eff
c----------------------------------------------------------------
c     ig     = 1 : the relativistic Green's function
c            = 0 : the non-relativistic approximation
c---------------------------------------------------------------------
      ig=1
c--------------------------------------------------------------------
c     n0 number of points in the Gaussian quadrature (ngauss = 64
c              is recommended.)
c---------------------------------------------------------------------
      n0=64
c--------------------------------------------------------------------
c     tc is a bulk electron temperature in keV
c--------------------------------------------------------------------
      tc=temp
c--------------------------------------------------------------------
c     tol = The relative error tolerence in numerical integration is set
c           to be MAX (tol, 1.0E-6).
c--------------------------------------------------------------------
c     tol=2.e-3
      tol=1.d-2
c--------------------------------------------------------------------
c     model Absolute value of model selectes collisional model: 1 for
c     full bounce av, 2 for square well(numerical solution), 3 for
c     analytic solution to square well.negative model does parallel
c     heating (lower hybrid, fast wave)
c
c     model < 5  gives rjpd in CURGAC
c           = 5  gives rjpd using the exact polarization-dependent
c                rf diffusion operator
c           > 5  gives rjpd using the polarization-dependent rf diffusion
c                operator with but small gyro-radius expansion 
c-------------------------------------------------------------------
      model=3
      model=5
c-------------------------------------------------------------------
c     z_eff is ion effective charge
c-------------------------------------------------------------------
c     lh = the cyclotron harmonic number
c--------------------------------------------------------------------
      lh=jwave 
c      lh=2
c      lh=0
c--------------------------------------------------------------------
c     elomom is ABS (lh)*cyclotron frec over wave frec for Y<0,
c       interpret as evaluated at poloidal angle where electrons are
c       first trapped in bucket.
c     elomom=yy = lh*omega_c/omega (y in Refs [1] and [2])
c--------------------------------------------------------------------
      elomom=lh*dabs(ye)
c--------------------------------------------------------------------
c     theta  = poloidal angle at which power is absorbed in radians
c              0: outborad; pi: inboard
c--------------------------------------------------------------------
      r=dsqrt(x*x+y*y)  ! needs work?
      prho=dsqrt((z-zma)**2+(r-rma)**2)
      ctheta=(r-rma)/prho
      stheta=(z-zma)/prho
      pi=4.d0*datan(1.d0)
c      write(*,*)' effcurba ctheta_geom,stheta',ctheta,stheta
      if(stheta.ge.0.0d0) then
         theta=dacos(ctheta)
      else
         theta=2.0d0*pi-dacos(ctheta)
      end if
c      write(*,*)' effcurba acos(ctheta),theta',acos(ctheta),theta
c--------------------------------------------------------------------
c     thtc is ratio of temperature along characteristic to that of
c          bulk, except for thtc <0, -thtc is energy (in keV) of bucket
c          rise; 
c          an obsolete variable
c--------------------------------------------------------------------
      thtc=1.0d0
c--------------------------------------------------------------------
c     enpar is k_parallel*c/wave frec. Note enpar**2 <1-elomom**2
c       implies no resonant particles; rjpd set to zero and rjpd0
c       determined for vparallel given by nonrelativistic resonance
c       condition
c--------------------------------------------------------------------
      enpar=cnpar
c--------------------------------------------------------------------
c     aspct is inverse aspct ratio
c--------------------------------------------------------------------
c bmax,bmin  in toray
c     if(dabs(bmax-bmin).lt.1.0d-8) bmax=bmin+1.0d-8
c 	aspct=(bmax-bmin)/(bmax+bmin)
c--------------------------------------------------------------------
c     Here :aspct is calculaed using rmax_psi and rmin_psi
c  conversion from real to double precision
      pid=4.d0*datan(1.d0)
c conversion from non-dimensional to m
      xd=dble(x)*0.01d0/r0x
      yd=dble(y)*0.01d0/r0x
      rd=dble(r)*0.01d0/r0x
      zd=dble(z)*0.01d0/r0x
      psid=psif_xyz(xd,yd,zd)
       zbmin=bmin_psi(psid)    
       zbmax=bmax_psi(psid)
       aspct=((zbmax-zbmin)/(zbmax+zbmin))
cSm060201
       if (aspct.lt.1.d-8) aspct=1.d-8 

c     function bxyz() uses the arguments x,y,z in [m]; 
      zb=bxyz(.01d0*x/r0x,.01d0*y/r0x,.01d0*z/r0x)!it changes 
                                         !bx,by,bz in one.i
      zacosarg=(zbmax + zbmin - 2.d0*zbmax*zbmin/zb)/
     &         (zbmin - zbmax)

c      write(*,*)'zbmax,zb,zbmin,zacosarg',zbmax,zb,zbmin,zacosarg
      if (dABS (zacosarg) .gt. 1.0d0) zacosarg = 
     1                               DSIGN (1.0d0, dble(zacosarg))
      theta =dACOS (zacosarg) 
      denom=1.d0
c      write(*,*)'effcurb befre curba denom,aspct,enpar,tc,thtc',
c     &denom,aspct,enpar,tc,thtc
c      write(*,*) 'theta,elomom,lh,z_eff,model,tol,n0,ig',
c     &            theta,elomom,lh,z_eff,model,tol,n0,ig

c-----using  TorGAcurba with real*8 arguments
      denom_d= denom
      aspct_d=aspct
      enpar_d=enpar 
      tc_d=tc
      thtc_d=thtc
      theta_d=theta
      elomom_d=elomom
      zeff_d=z_eff
      tol_d=tol
      omode=dble(float(ioxm)) ! YuP Arg.#8 should be real*8

      call TorGA_curgap(rjpd_d,rjpd0_d,ratjpd_d,denom_d,
     &aspct_d,enpar_d,
     &cnper, omode,cefldx,cefldy,cefldz,
     &tc_d,thtc_d,theta_d,elomom_d,lh,zeff_d,model,tol_d,n0,ig)

      rjpd=rjpd_d
      rjpd0=rjpd0_d
      ratjpd=ratjpd_d
      write(*,*)'prep3d after TorGA_curgap' 
      write(*,*)'rjpd,rjpd0,ratjpd',rjpd,rjpd0,ratjpd

c----------------------------------------------------------
c     efficiency  in (A/cm**2)/(erg/(sec*cm**3))
c     temperature temp in kev
c     density     den in 10**19/m**3
c-----------------------------------------------------------
c-----genray transformation of CD efficiency 
c     cln -Coulomb Ln
c     it is necessary to write the formula to calculate cln

c      cln=17.0
c      arg1=1.d3/temp*dsqrt(10.d0*den) !temp KeV, den 10**13/cm**3
c      cln=24.d0-alog(arg1)
       efficien=rjpd 
c      efficient=efficien*temp/den*(17.0/cln)*4.5*1.e-6
c      write(*,*)'efficiecy coefficient=temp/den*(17.0/cln)*4.5*1.e-6',
c     1temp/den*(17.0/cln)*4.5*1.e-6
c      write(*,*)'in effic rjpd,rjpd0,ratjpd,denom',
c     1           rjpd,rjpd0,ratjpd,denom
c      write(*,*)'end effic curba efficient',efficient

c-----the efficiency transformation from torayGA
c      write(*,*)'before zfac_f temp',temp
      zfac=zfac_f_xyz((x*1.d-2),(y*1.d-2),(z*1.d-2),dble(temp))
      zfac=zfac*1.d-7 ! watt => egr/sec
      zfac=zfac*1.d-4 ! m**2 => cm**2
c      write(*,*)'zfac',zfac
      efficient=efficien*zfac
cSm080629
c      efficient=-efficient 

      return
      end



c======================================================================
c======================================================================


      double precision function zfac_f_xyz(x,y,z,temp_kev) 
c-----It calculates the coefficient
c     that transforms the CD efficiency from
c     curba variables to A/cm**2/erg/sec*cm**3))
c     It uses the formula from currn in toray.
c-----input  x(m),y(m),z(m)
c            temp_kev is the electron temperature in keV    
c     It uses function wpw_2() that uses the small radius rho.

      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      double precision ld
c-----externals: wpw_2
      data zelchg/4.8032d-10/, zconst/10.83259d0/,zconvr/33.33333d0/
      data cvac/2.9979d10/
c       write(*,*)'in zfac_f z,r,phi,temp_kev', z,r,phi,temp_kev
      alfa=wpw_2(x,y,z,1) !(omega_pe/omega)**2
      zralfa=dsqrt(alfa)

c     Subroutine to interface between TORAY and R. Cohen's calculation
c     of current.
c
c     zelchg is the electron charge in statcoulombs.
c     zconst is the ln of Boltzmann's constant (1.3807e-16) times the
c            factor to convert temperature from eV to degrees Kelvin
c            (1.1605e4) divided by the product of cvac times Planck's
c            constant (1.0546e-27).
c     cvac   is the light speed sm/sec
c     ld     is cvac divided by omega (2*pi*f).
c     zconvr is the constant that converts from (statamps/cm**2) divided
c            by (ergs/sec) to (amps/m**2) divided by watts.
      
      pi=4.d0*datan(1.d0)
      ld=cvac/(2.d0*pi*frqncy*1.d9) !c/omega cm
      zlds=ld*ld
      zfac1=zconst+log(ld)
      zte=temp_kev*1.d3 !eV       
      alfa=wpw_2(x,y,z,1) !(omega_pe/omega)**2
      zralfa=sqrt(alfa)
      zfac2=(zlds /(zelchg*alfa))*zte*2.d0/511.0d3
      zfac3=zfac1+ log (zte/zralfa)
      zfac_f_xyz=zconvr*zfac2/zfac3

      return
      end
      
      

c======================================================================
c======================================================================
      double precision function psif_xyz(x,y,z) 
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i' ! contains bz0
      include 'five.i'
      include 'three.i' 
      include 'fourb.i' ! contains req,zeq,xeq,yeq grids
      double precision ias2r,ias2r_Sm
         r=dsqrt(x*x+y*y) ! ok for now, change to xyz-spline later
         r=min(r,req(nreqd)) ! psi is not defined outside of req grid
         ! nr4a , nz4a and nrya  are given in param.i 
         nr4=nx+4
         nz4=ny+4
         ncx=nr4
         ncy=nz4
         psif_xyz=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,0,0,r,z,nr4a)
      return
      end


c======================================================================
c======================================================================
c        *********************edgcor_xyz **********************
c        *                        -                           *
c        * this subroutine  shifts the point(where	      *
c        * the EC ray intersects the plasma)inside the plasma * 
c        ******************************************************
c
c------------------------------------------------------------------
c								   !
c        input:					   !
c        x,y,z   coordinates before the shift                   !
c        output:				   !
c        xu0,yu0,zu0              	   !
c------------------------------------------------------------------
      subroutine edgcor_xyz(x,y,z, 
     + xu0,yu0,zu0) !-> out
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
      double precision ias1r,ias2r,ias2r_Sm

      if (model_rho_dens.ne.0) then
        xu0=x
        yu0=y
        zu0=z
        return
      end if

     
      epsbnd=2.d-3 ![m]! 
      delcor=1.d0-epsbnd

      rho_start=1.d0
      ! YuP: For direct launch of O-mode, no need to step inside rho=1.
      ! YuP: Specify rho_start here
      if(jwave.eq.1 .and. ibw.eq.0) rho_start=1.43d0
      
      rho_start =min(rho_start,rhowall)

1     continue ! handle for iterations
      iedg=0
      
      r=dsqrt(x*x+y*y)
      xma= rma*x/r
      yma= rma*y/r  ! ok for now?
      
      if ( (x.lt.xeqmin+epsbnd).or.(x.ge.xeqmax-epsbnd) .or.
     +     (y.lt.yeqmin+epsbnd).or.(y.ge.yeqmax-epsbnd) .or.
     +     (z.lt.zeqmin+epsbnd).or.(z.ge.zeqmax-epsbnd)  ) then 
         ! Outside of grid
        write(*,*)'edgecor_xyz: r,rmin,rmax===',r,rmin,rmax
        write(*,*)'edgecor_xyz: y,yeqmin,yeqmax===',y,yeqmin,yeqmax
        write(*,*)'edgecor_xyz: z,zeqmin,zeqmax===',z,zeqmin,zeqmax
        !pause
        iedg=1
      endif

c-----------------------------------------------------
      ! Yu.P. Added: skip checking on LCFS   
      goto 5   

      if(model_rho_dens.eq.0) then
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
        if ((z.gt.zp-epsbnd).or.(z.lt.zm+epsbnd)) then
          iedg=1
        end if
      endif
      
  5   continue
  
      ! get rho based on model for density (model_rho_dens=1,2,4), 
      ! 2D spline of data (model_rho_dens=3),
      ! or based on magnetic flux (model_rho_dens=0,5):
      den=dense_xyz(x,y,z,1) !-> get rho (stored in one.i) 
      !write(*,*)'edgecor_xyz: x,y,z,rho===',x,y,z,rho
      if(rho.gt.rho_start) then
        iedg=1
      endif
        
      if (iedg.eq.1) then ! another iteration
c       shift (x,y,z) towards (0,0,0) point
        delx=x-xma
        dely=y-yma
        delz=z-zma
        x= xma+delcor*delx
        y= yma+delcor*dely
        z= zma+delcor*delz
        write(*,*)'edgecor_xyz: x,y,z,rho===',x,y,z,rho
        goto 1
      endif
      xu0=x
      yu0=y
      zu0=z
      return
      end


c======================================================================
c======================================================================

      subroutine vacray_xyz(x0,y0,z0,p,cnx,cny,cnz, 
     +                      xp,yp,zp) !-> out
c----------------------------------------------------------------------
c       this subroutine determines the vacuum ray point(xp,yp,zp)
c      	input: x0,y0,z0,p,cnx,cny,cnz
c       output:xp,yp,zp
c---------------------------------------------------------------------
      implicit integer (i-n), real*8 (a-h,o-z)
      xp= x0+p*cnx
      yp= y0+p*cny
      zp= z0+p*cnz
      ! YuP: needs work / checking?
      return
      end
	

c======================================================================
c======================================================================

c        *********************plasmray_xyz *******************
c        *                        -                           *
c        * this subroutine  determines the point where	      *
c        * the EC ray intersects the plasma                   *
c        * boundary(zu0,ru0,phiu0)                            *
c        ******************************************************
c
c------------------------------------------------------------------
c								   !
c        input parameters					   !
c        xst,yst,zst [m] -start point cartesian coordinates	   !
c        alfast,betast-tor. and poloid.angl. of the start	   !
c        refractive index in the vacuum(radian)                    !
c        output parameters					   !
c        xu0,yu0,zu0                                  	   !
c        iraystop -parameter for stoppage the ray calculation 	   !
c                  if iraystop =1                                  !    
c------------------------------------------------------------------
c        this program uses the following functions and subroutines !
c        double precision function proot(p2,p1,p0)                 !
c        subroutine vacray_xyz()	   !
c        double precision function psif_xyz()  		           !
c------------------------------------------------------------------
      subroutine plasmray_xyz(xst,yst,zst,alfast,betast,
     1                        xu0,yu0,zu0,iraystop) !-> out
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
      !---------
      
      if (model_rho_dens.ne.0) then
        xu0=xst
        yu0=yst
        zu0=zst
        return
      end if

      
      epsbnd=2.d-3 ![m] 
      rho_start=1.d0
      ! YuP: For direct launch of O-mode, no need to step inside rho=1.
      ! YuP: Specify rho_start here
      if(jwave.eq.1 .and. ibw.eq.0) rho_start=1.43d0
      
      rho_start =min(rho_start,rhowall)
        
      iraystop=0 ! initialize      
c-------------------------------------------------------------
c     check: Is the starting point is inside the plasma?
      ! get rho based on model for density (model_rho_dens=1,2,4), 
      ! 2D spline of data (model_rho_dens=3),
      ! or based on magnetic flux (model_rho_dens=0,5):
      den=dense_xyz(xst,yst,zst,1) !-> get rho (stored in one.i) 
      if (rho.le.rho_start) then
        write(*,*)'starting point is inside rho_start=',rho_start
        xu0=xst
        yu0=yst
        zu0=zst
        return
      end if
c---------------------------------------------------------------
      cosbet=dcos(betast)
      sinbet=dsin(betast)
      cosalf=dcos(alfast)
      sinalf=dsin(alfast)
      cnx=cosbet*dcos(alfast+phist)
      cny=cosbet*dsin(alfast+phist)
      cnz=sinbet
c-----------------------------------------------------------
c     vacuum ray trajectory
c     X=X_st+p*cnx
c     Y=Y_st+p*cny
c     Z=Z_st+p*cnz
c     p>0 -parameter
c-----------------------------------------------------------
      xu0=xst
      yu0=yst
      zu0=zst

      write(*,*)
     +'plasmray_xyz/start iterations: shift starting point to rho_start'
      
10    continue ! handle for iterations

c----> YuP: use grad(el.density) instead of grad(psi)
c----> to find the unit vector pointing into plasma.
      if (model_rho_dens.ne.0) then
        xu0=xst
        yu0=yst
        zu0=zst
        return
      end if

      call dwpw_2(xu0,yu0,zu0,1, dwpw2dx,dwpw2dy,dwpw2dz)
      write(*,'(a,4e12.3)')'plasmray_xyz: x,y,z,rho',xu0,yu0,zu0,rho
      if (rho.le.rho_start) then
         write(*,*)'starting point is inside rho_start=',rho_start
         return
      end if
      ! Vector along grad(el.density)  :
      dndx= dwpw2dx ! strictly, it is d(wpe/w)^2/dx,
      dndy= dwpw2dy ! but we only need a sign/direction of grad(ne)
      dndz= dwpw2dz ! Note: grad(ne) <0 (pointing inward)   
      grne= dsqrt(dndx**2 + dndy**2 + dndz**2)
      pp=1.d0/grne
      ! Unit vector along grad(ne):
      dndx= dndx*pp
      dndy= dndy*pp
      dndz= dndz*pp
      ! Dot product (N.grad(ne)) :
      cn_grne= (dndx*cnx + dndy*cny + dndz*cnz)
c      write(*,'(a,3e12.3)')'plasmray_xyz:dndx,dndy,dndz=',dndx,dndy,dndz
c      write(*,'(a,3e12.3)')'plasmray_xyz: Nx,Ny,Nz=',cnx,cny,cnz
      ! If cn_grne>0, N-vector is in the same dir. as grad(ne), inward
      if(cn_grne.le.0.d0) then
         write(*,*)'plasmray_xyz: ray cannot go into plasma'
         write(*,*)'  iraystop->1'
         iraystop=1
         return
      endif
      ! cn_grne>0 (refractive index vector pointing inwards):
      ! shift the point inwards:
      rrr0=sqrt(xu0**2+yu0**2+zu0**2)
      xu0= xu0+ rrr0*dndx*epsbnd  ! shifting inwards
      yu0= yu0+ rrr0*dndy*epsbnd  ! shifting inwards
      zu0= zu0+ rrr0*dndz*epsbnd  ! shifting inwards
      go to 10
      
      write(*,*)'end of plasmray_xyz'
      return
      end


c======================================================================
c======================================================================

      subroutine bfield_coils(x,y,z,                        ! model_b=2
     +                        b_x, b_y, b_z) !-> out
c     Yu.P. 2011
c     Calculate magnetic field components 
c     produced by coils (loops) with currents.
c
c     Input: 
c     x,y,z - coordinates [m] where the field is calculated.
c
c     Input from one.i (Specify in namelist &tokamak):
c     ncoils - number of coils.
c     radc(ncoils)  - radii of coils [m].
c     zcoil(ncoils) - z-coordinates of coils' centers;
c                     the centers are on the z-axis,
c                     with loop's plane perpendicular to the z-axis.
c     curc(ncoils)  - current in each coil [A].
c
c     Output:
c     b_x,b_y,b_z - magnetic field components [Tesla] at (x,y,z).
c
c     Uses COMPLETE ELLIPTICAL INTEGRALS K() & E() 
c     expressed through Carlson integrals DRF(),DRD()

      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'

      upi2 = 2.d-7  ! == mu0/2pi == 12.56637e-07/(2*3.1415926535898)
      R2   = x*x+y*y
      R    = dsqrt(R2)
      o3   = 1.d0/3.d0

      b_x= 0.d0
      b_y= 0.d0
      b_r= 0.d0
      b_z= 0.d0
      
      r_eps= 0.001*minval(radc(1:ncoils)) !YuP[04-2016]
      !If the (R,Z) point happen to get onto the (Rcoil,Zcoil) position
      !it may result in singularity.
      !To avoid this situation, the value of R-radc(i) 
      !is adjusted to avoid the exact equality R=radc(i),
      !(when the distance from given point to the coil (wire) 
      !is smaller than r_eps).
      !The effect from this correction seems to be negligible,
      !but could allow avoiding 1/0 singularity.
       
      do 10 i = 1,ncoils !-----------------> sum over coils
         R_Rc = R - radc(i)
         z_zc = z - zcoil(i)
         dist2= R_Rc*R_Rc +z_zc*z_zc
         dist= sqrt(dist2) !distance from given point (x,y,z) to the coil
         if(dist.lt.r_eps)then
           R_Rc = r_eps*0.5 ! adjusted
           !z_zc = r_eps*0.5 ! alternatively, could also be adjusted
           dist2= R_Rc*R_Rc +z_zc*z_zc ! re-evaluated
           dist= sqrt(dist2) !distance from given point (x,y,z) to the coil
         endif
         rrad = 2.0*radc(i)*R
         radc2= radc(i)*radc(i)
         r2zc = R2 + z_zc*z_zc ! = radc^2 at coil position
         ar1  = radc2 + r2zc + rrad ! = 4*radc^2 at coil position
         sq1  = curc(i)/dsqrt(ar1)
         ez1  = (radc2-r2zc)/dist2 ! 0/0 at (R=radc;Z=zcoil)
         ak2  = 2.d0*rrad/ar1 != k^2
         !if(1.d0-ak2.le.0.d0)then
         !  write(*,*)' bfield_coils: 1.d0-ak2=',1.d0-ak2
         !  !!pause
         !endif
      !write(*,'(a,i4,2e12.3)')' bfield_coils: i,ar1,ak2',i,ar1,ak2
         elk1 = DRF(0.d0,1.d0-ak2,1.d0,ifail) != K(k^2)      IMSL: ELK(ak2)
      ! write(*,'(a,2i4,e12.3)')' bfield_coils: i,ifail,elk1',i,ifail,elk1
         ele1 = elk1- o3*ak2*DRD(0.d0,1.d0-ak2,1.d0,ifail) !E(k^2) ELE(ak2)
      ! write(*,'(a,2i4,e12.3)')' bfield_coils: i,ifail,ele1',i,ifail,ele1
         ex1  = 0.5d0*(2.d0-ak2)/(1.d0-ak2)
         elz1 =   elk1 + ele1*ez1
         elx1 = (-elk1 + ele1*ex1)*z_zc
         b_z  = b_z + elz1*sq1
         b_r  = b_r + elx1*sq1
   10 continue
     
      !write(*,*)' bfield_coils: b_r,b_z',b_r,b_z

      b_z= b_z*upi2
      if(r.ne.0.d0) then
         one_r = 1.d0/r
         b_r = b_r*upi2*one_r
         b_x = b_r*x*one_r
         b_y = b_r*y*one_r
      endif

      return
      end

c======================================================================
c======================================================================

      subroutine eq_mirror1(x,y,z,             ! model_b=3 YuP[04-2016]
     +                       PSI, b_x, b_y, b_z) !-> out)  
c Calculate values of pol. flux PSI (norm-ed by 2pi), 
c and components of B field for eqsource="mirror1" model:    
c     A model for magnetic field (Br,Bz) in a mirror machine, 
c     set over (R,Z) plane,
c     assuming symmetry in toroidal (azimuthal) angle
c     and min value of |B| being at Z=0 point
c     for each field line PSI=const (constant poloidal flux value).
c     The model is defined by 
c        Bz(R,Z)=  B00*J0(R/glb)*cosh(Z/glb)
c        Br(R,Z)= -B00*J1(R/glb)*sinh(Z/glb)
c     where B00 is the value of |B| at (Z=0,R=0),
c     glb is the scale length of magnetic field increase.
c     J0 and J1 are the Bessel functions of argument (R/glb).
c     The "toroidal" (azimuthal in cylindr.geometry) component
c     of magnetic field is assumed to be zero for now, BPHI=0. 
c     In addition, the poloidal flux function can be found as 
c        PSI= INTEGRAL[0_R]{Bz*R*dR)= B00*R*glb*J1(R/glb)*cosh(Z/glb)
c     It can be verified that for the model equations above,
c     div.B =0 in cylindrical coords,
c     which is (using d/dphi =0): (1/R)d(R*Br)/dR + d(Bz)/dZ = 0.
c     Also it can be shown that curl(B) = dBr/dZ - dBz/dR == 0,
c     for the specified model equations.
c     Note that the mirror ratio is Bz(Zmax,R=0)/B00 = cosh(Zmax/glb).
c     From here, the value of glb can be found 
c     if the mirror ratio is known:
c       glb_mirror= 0.5*zbox_mirror/acosh(rmirror) 
c     ! Note: zbox_mirror=2*Zmax
c     (set zbox_mirror [m] and rmirror in genray.in).
c     (Also set b00_mirror [T] in genray.in, for the value of B00).
c
c     Input: 
c        x,y,z - coordinates [m] where the field is calculated;
c           also b00_mirror [T] and glb_mirror [m] - access through one.i
c     Output: poloidal flux (norm. by 2pi) PSI
c        and  b_x,b_y,b_z - magnetic field components [Tesla] at (x,y,z).

      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i' ! contains b00_mirror and glb_mirror
      
      ! local working arrays:
      real*8  xarg_bess ! argument of Bessel functions
      integer nb_bess
      real*8  alpha_bess
      real*8  b_bess(2) ! J0=b_bess(1) and J1=b_bess(2) 
                        ! SHOULD MATCH nb_bess !

      r2= x*x+y*y
      r=  dsqrt(r2)
      Rc= r ! radius in cyl.coord.

      b_x= 0.d0
      b_y= 0.d0
      b_r= 0.d0
      b_z= 0.d0

      nb_bess=2  ! 1+highest order to be calculated (=2 for J0 and J1)  
      alpha_bess=0.d0 ! for J_{alpha+k-1} functions (k=1:nb_bess)
      xarg_bess= Rc/glb_mirror
      ncalc= nb_bess ! if an error happens inside DBESJ (or zzbeslri),
                     ! the value of ncalc will be changed.
C     ncalc - number of components of b_bess(i) set to zero due to
C                    underflow,
C                    ncalc=0   , normal return, computation completed
C                    ncalc.NE.0, last ncalc components of b_bess(i) 
C                       set to 0., b_bess(i)=0.0D0, 
C                                         i=nb_bess-ncalc+1,...,nb_bess.

cYuP: this version is used in CQL3D distribution, for the Bessel functions:
c      call zzbeslri(xarg_bess,nb_bess,ize_bess,b_bess,ncalc) !in zcunix.f
c     Local version:
      call DBESJ(xarg_bess,alpha_bess,nb_bess, b_bess,ncalc) 
      
      if (ncalc.lt.0) then
        WRITE(*,*)'eq_mirror1/bessel: an argument is out of range.'
        stop
      endif
      
c      if(z.eq.0.d0)then
c      write(*,'(a,3e13.3)')'eq_mirror1: Rc/glb_mirror,J0,J1=',
c     +                 xarg_bess, b_bess(1), b_bess(2)
c      endif
     
      Z_glb= z/glb_mirror
      b_z=  b00_mirror*b_bess(1)*cosh(Z_glb) ! ~J0(R/glb)
      b_r= -b00_mirror*b_bess(2)*sinh(Z_glb) ! ~J1(R/glb)
      !bphi= 0.d0 ! Could be added later through rbphi0 ?
      if(r.ne.0.d0) then
         one_r = 1.d0/r
         b_x = b_r*x*one_r
         b_y = b_r*y*one_r
      endif
      PSI= b00_mirror*Rc*glb_mirror*b_bess(2)*cosh(Z_glb) ! ~J1(R/glb)
      ! Note that PSI=0 at R=0 (at the axis), and grows with R 
      ! if b00_mirror>0
      ! (grows up to the first null point of J1 Bessel function;
      !  the range in R is limited by the R_b0, see equilib)
      ! Note: In cql3d, PSI must be descending from magnetic axis 
      ! to the edge, so the value of b00_mirror is reset in cql3d
      ! to a negative value to satisfy this condition.
      ! In Genray, b00_mirror can be positive or negative;
      ! The pol.flux PSI will also be consistent,
      ! But for the code purposes, array peqd()is always ascending
      ! func. of rho. The physical sign is tracked by value of dpsimax.
      
      return
      end


c======================================================================
c====================================================================== 
      subroutine eq_mirror1_field_line(PSI,Rline,Zline,nline) !YuP[03-2016]
c Define arrays Rline(nline) and Zline(nline) that describe 
c the magnetic field line in mirror1 model, 
c for a given value of pol.flux PSI.
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i' !contains b00_mirror, glb_mirror, zbox_mirror
c In mirror1 model, the poloidal flux function (norm-ed by 2pi)
c is given by 
c PSI= INTEGRAL[0_R]{Bz*R*dR)= B00*R*glb*J1(R/glb)*cosh(Z/glb)
c The field line is determined by such set of points
c Rline(nline) and Zline(nline) that
c (Rline/glb)*J1(Rline/glb)*cosh(Zline/glb) = PSI/(B00*glb**2)
c where PSI is fixed for a given field line.

c INPUT: PSI= pol.flux (norm-ed by 2pi) 
c        nline= number of points along line.
c OUTPUT: Rline(nline) and Zline(nline)

      real*8 Rline(nline) , Zline(nline)
      
      ! All lines in mirror machine pass through 
      ! Zmin=-0.5*zbox_mirror and Zmax=+0.5*zbox_mirror,
      ! which are the 1st and last points in ez grid.
      ! Define Zline points equidistantly in that range.
      if(nline.le.1)then ! should not normally happen
        dzline=0.d0
      else ! nline>1  normal case
        dzline= zbox_mirror/(nline-1)
      endif
      
      glb=glb_mirror

      ! For solving 
      !(Rline/glb)*J1(Rline/glb) = PSI/(B00*glb**2*cosh(Zline/glb))
      ! define
      psi_b00= PSI/(b00_mirror*glb**2)
      do iline=1,nline
         Zline(iline)= -0.5*zbox_mirror + dzline*(iline-1)
         ! For a given value of Zline, we need to solve
         ! (Rline/glb)*J1(Rline/glb) = PSI/(cosh(Zline/glb)*B00*glb**2)
         ! to find the unknown value of Rline.
         ! Find R from searching the root
         ! (R/glb)*J1(R/glb) = cz0 , where
         cz0= psi_b00/cosh(Zline(iline)/glb)
         ! Approximate analyt. formula (by expansion of J1):
         Rline(iline)= glb*sqrt(cz0*(2.0+0.5*cz0+0.25*cz0*cz0))
         ! It is quite accurate for R/glb <1 :
         ! 0.1% error at R/glb=1.0, 1% at R/glb=1.35,
         ! but not so good near 1st null of J0(R/glb), i.e. R/glb=2.4048 :
         ! 5% at R/glb=1.7, and 12% error at R/glb=1.9.
         ! Ok for graphics.
      enddo
      
      return
      end  
c======================================================================
c======================================================================

      subroutine density_profile_read  ! Only used for model_rho_dens=3

c     Yu.P. 2011
c     Read the density profile(x,y) data from file
c-------------------------------------------------------------------!
! Read 2D density profile from file.
! Assumed: uniform (x,y)-grid,
! Density profile is a function of (x,y) only.
!===> NEED TO SPECIFY THE NAME OF INPUT FILE in genray.in:
! dendsk="19846_0_25densf.dat"   (example)
! Check that the data in file is in these units: 10^19 [m^-3]
c-------------------------------------------------------------------!
c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none 

      include 'param.i'
      include 'one.i'
      include 'fourb.i'
      include 'five.i'

c---- local: 
      real*8 den_unit,dstep, cnorm,den_t ,dn,dn_min,dn_max,dndx,dndy  
      integer i,j , lx,ly,mxa,mxb,mya,myb,nfx,nfy,nxy

      double precision
     1     arfxa(nzeqda),arfxb(nydena),arfya(nxdena),arfyb(nxdena),
     1	   xwk(nxdena),ywk(nydena),
     2	   denwk(nxdena,nydena)
c---- external:     
      double precision	ias1r,ias2r,ias2r_Sm
      
c     Data from /plasma/  namelist were read
c     in genray.f using read_all_namelists
      
      ! den_unit= 1.0d19 ! [m^-3] 
      ! Check that the data in file is in these units !
      
c Read the dendsk file
      write(*,'(//,a)') ' density_profile_read: dendsk=',dendsk
c---------------------------------------------------------
      open(30,file=dendsk)
      
c Read the dimensions of x-grid and y-grid:
      read(30,2) nxden,nyden
2     format(2i4)
      write(*,*) 'nxden,nyden=',nxden,nyden

      if(nxden.gt.nxdena)then
        write(*,*) 'density_profile_read: nxden>nxdena'
        stop 'Increase nxdena in param.i'
      endif
      if(nyden.gt.nydena)then
        write(*,*) 'density_profile_read: nyden>nydena'
        stop 'Increase nydena in param.i'
      endif

c Read the min/max values for the x-grid and y-grid:
c (We assume that the grids are uniform)
      read(30,3) xdenmin,xdenmax, ydenmin,ydenmax ! Units:[m]
3     format(4e16.9)
      write(*,*) 'xdenmin, xdenmax [m]=',xdenmin,xdenmax
      write(*,*) 'ydenmin, ydenmax [m]=',ydenmin,ydenmax

c Read the profile of density on the (x,y)-grid:
      read(30,4) ((dengrid(i,j),i=1,nxden),j=1,nyden)
4     format(5e16.9)
c Note: data is recorded in 5 columns !
      denmin=MINVAL(dengrid, MASK=dengrid.gt.0.d0)
      denmax=MAXVAL(dengrid)
      write(*,*) 'dengrid min/max [file]:',denmin,denmax
      write(*,*) 'The units in file should be [10^19 m^-3]'

      close(30)

c----- GENERATE GRIDS -----------------------------------------
      dstep= (xdenmax-xdenmin)/(nxden-1)
      do i=1,nxden
        xden(i)= xdenmin+dstep*(i-1) ! [xdenmin; xdenmax]
        xwk(i)=xden(i)
      enddo
 
      dstep= (ydenmax-ydenmin)/(nyden-1)
      do j=1,nyden
        yden(j)= ydenmin+dstep*(j-1) ! [ydenmin; ydenmax]
        ywk(j)=yden(j)
      enddo

c----- GENERATE SPLINE COEFFS for dengrid -----------------------
      nx_den=nxden
      ny_den=nyden
      lx=1
      ly=1
      mxa=0
      mxb=0
      mya=0
      myb=0  ! yup: presumably 3 means periodic conditions (?)
      ! Working array for spline coeffs.
      do i=1,nxden
      do j=1,nyden
         denwk(i,j)=max(dengrid(i,j),denmin) 
         ! Note the Lower Limit for denwk; modify if needed.
      enddo
      enddo
      do i=1,nyden
         arfxa(i)=0.
         arfxb(i)=0.
      enddo
      do i=1,nxden
         arfya(i)=0.
         arfyb(i)=0.
      enddo
      nfx=nxden
      nfy=nyden
      ncx_den=nxden+4
      ncy_den=nyden+4
      nxy=max(nxden,nyden)+4
      call iac2r_Sm(xwk,nx_den,ywk,ny_den,
     +              denwk,nfx,nfy,lx,ly,
     &              mxa,arfxa,mxb,arfxb,
     &              mya,arfya,myb,arfyb,
     +              tx_den,ty_den,cxy_den,ncx_den,ncy_den,
     +              nxy,cy_den,
     &              nxdena,nx4a) 
      write(*,*) 'spline coefficients for dengrid were created'

ctest     
      cnorm=0.d0
      dn_max=0.d0     
      dn_min=100*(denmax-denmin)/xdenmax
      do i=1,nx_den
      do j=1,ny_den
         den_t=ias2r_Sm(tx_den,nx_den,ty_den,ny_den,
     +                  cxy_den,ncx_den,ncy_den,0,0,xwk(i),ywk(j),nx4a)
         cnorm=cnorm+(dengrid(i,j)-den_t)**2
         if(den_t.ge.denmin) then
           ! dn/dx:
           dndx=ias2r_Sm(tx_den,nx_den,ty_den,ny_den,
     +                  cxy_den,ncx_den,ncy_den,1,0,xwk(i),ywk(j),nx4a)
           ! dn/dy:
           dndy=ias2r_Sm(tx_den,nx_den,ty_den,ny_den,
     +                  cxy_den,ncx_den,ncy_den,0,1,xwk(i),ywk(j),nx4a)
           dn=sqrt(dndx**2+dndy**2)
           dn_min=min(dn_min,dn) ! min{|grad(n)|}
           dn_max=max(dn_max,dn) ! max{|grad(n)|}
         endif
      enddo
      enddo
      cnorm=dsqrt(cnorm)/(nx_den*ny_den*denmax)
      write(*,*)'cnorm/denmax of dengrid',cnorm
      write(*,*)'|grad(n)| min/max=',dn_min,dn_max
      write(*,*)'(denmax-denmin)/xdenmax=',(denmax-denmin)/xdenmax
cendtest

      return
      end


c======================================================================
c======================================================================
      double precision function func_kappa(arg)    ! Only for model_b=4
c     Yu.P. 2012
c     Function beta_frc*arg-tanh(arg)
c     The root of this function, i.e. the value of arg
c     that makes func_kappa(arg)=0, corresponds to "Kappa"
c     which is used in definition of density profile and Bz profile
c     in FRC-like plasma.
c-------------------------------------------------------------------!
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i' ! stores beta_frc
      func_kappa= beta_frc*arg - tanh(arg)
      end

c======================================================================
c======================================================================

      subroutine bfield_frc(x,y,z,                          ! model_b=4
     +                        b_x, b_y, b_z) !-> out
c     Yu.P. 2012
c     Calculate magnetic field components 
c     that correspond to an FRC plasma.
c
c     Input: 
c     x,y,z - coordinates [m] where the field is calculated.
c
c     Input from one.i (Specify in namelist &tokamak):
c     rs_frc  - separatrix radius (specify in genray.in) [m]
c     bz0     - magnetic field at r=+INF, approximately at wall [T]
c     Input from one.i (calculated in equilib.f):
c     akappa  - "Kappa" value such that beta_frc*Kappa-tanh(Kappa)=0
c
c     Output:
c     b_x,b_y,b_z - magnetic field components [Tesla] at (x,y,z).
c

      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'

      r2= x*x+y*y
      if(model_b.eq.4) then
         ! In the present model, the magnetic field has only one 
         ! component (Bz) which is reversed at r=rs_frc/sqrt(2) radius.
         b_x= 0.d0
         b_y= 0.d0
         b_z= bz0*tanh( akappa*(2.*r2/rs_frc**2 -1.d0) )
         ! B_z=0 at r= r0= rs_frc/sqrt(2)
         ! B_z -> +bz0 at very large r
         ! B_z ~~ -bz0 at r=0
      endif
      ! Later: Could add Solov'ev model ? (maybe as model_b.eq.3)
      return
      end

c======================================================================
c======================================================================

      subroutine dnd_xyz(x,y,z, cnx,cny,cnz,
     .dnpdx,dnpdy,dnpdz,    dnpdcnx,dnpdcny,dnpdcnz,
     .dnlldx,dnlldy,dnlldz, dnlldcnx,dnlldcny,dnlldcnz,
     .dnpdw,dnlldw)
c----------------------------------------------------------
c     Yu.P. 2012
c
c     calculates the derivatives of N_perp=np and N_parallel=nll
c     over x,y,z, Nx,Ny,Nz
c
c     it uses bx,by,bz from common one.i
c
c     bxyz(x,y,z) should be run before this subroutine.
c     bxyz calculates magnetic field components
c     and derivatives of B by x,y,z  
c     (stored one.i)
c--------------------------------------------------------- 
      implicit none
c      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
c     input:
      double precision x,y,z, cnx,cny,cnz
c     output:
      double precision cnp,nll,
     .dnpdx,dnpdy,dnpdz,    dnpdcnx,dnpdcny,dnpdcnz,    !derivatives N_perp
     .dnlldx,dnlldy,dnlldz, dnlldcnx,dnlldcny,dnlldcnz, !derivatives N_par
     &dnpdw,dnlldw
c     local:
      real*8 o_cnp
      
      nll= (cnx*bx+cny*by+cnz*bz)*o_bmod ! o_bmod is 1/|B|
      cnp= dsqrt(cnx**2+cny**2+cnz**2 -nll**2)
c      write(*,*)'dnd_xyz: rho, b=',rho,bx,by,bz
c      write(*,*)'dnd_xyz: dbz/d =',dbzdx,  dbzdy,  dbzdz
c      write(*,*)'dnd_xyz: nll,cnp',nll,cnp
 
c-----derivatives of N_parallel    
      dnlldx= (cnx*dbxdx+cny*dbydx+cnz*dbzdx -nll*dbmdx)*o_bmod     
      dnlldy= (cnx*dbxdy+cny*dbydy+cnz*dbzdy -nll*dbmdy)*o_bmod     
      dnlldz= (cnx*dbxdz+cny*dbydz+cnz*dbzdz -nll*dbmdz)*o_bmod     
      dnlldcnx= bx*o_bmod
      dnlldcny= by*o_bmod
      dnlldcnz= bz*o_bmod
      dnlldw=  -nll/frqncy  !dN_parallel/d_omega

c-----derivatives of N_perpendicular  
      o_cnp= 1.d0/cnp
      dnpdx= -nll*dnlldx*o_cnp
      dnpdy= -nll*dnlldy*o_cnp
      dnpdz= -nll*dnlldz*o_cnp
      dnpdcnx= (cnx-nll*dnlldcnx)*o_cnp
      dnpdcny= (cny-nll*dnlldcny)*o_cnp
      dnpdcnz= (cnz-nll*dnlldcnz)*o_cnp
      dnpdw=  -cnp/frqncy   !dN_perpendicular/d_omega
c	write(*,*)'dnd dnlldcnz,dnpdcnz',dnlldcnz,dnpdcnz
      return
      end


c======================================================================
c======================================================================
	subroutine dtemp_xyz(x,y,z,i, dtemp_dx, dtemp_dy, dtemp_dz)
c----------------------------------------------------------
c     Yu.P. 2012
c
c     calculates the derivatives of temperature of species 'i'
c     over x,y,z
c----------------------------------------------------------
      !implicit integer (i-n), real*8 (a-h,o-z)
      implicit none
      include 'param.i'
      include 'one.i' ! stores rho, drho_dx,...
      real*8 rho_dens_xyz, dtemdrho ! external
      real*8 x,y,z ! Input
      integer i ! Input
      real*8 rholoc, dtemp_drho
      real*8 dtemp_dx,dtemp_dy,dtemp_dz ! Output
      !--------------
	rholoc= rho_dens_xyz(x,y,z) ! get rho and derivs of rho
	dtemp_drho= dtemdrho(rholoc,i) ! spline
	dtemp_dx= dtemp_drho*drho_dx
	dtemp_dy= dtemp_drho*drho_dy
	dtemp_dz= dtemp_drho*drho_dz
      return
      end


c======================================================================
c======================================================================
	subroutine dtpop_xyz(x,y,z,i, dtpop_dx, dtpop_dy, dtpop_dz)
c----------------------------------------------------------
c     Yu.P. 2012
c
c     calculates the derivatives 
c     of tpop=T_perp/T_parallel 
c     of species 'i'
c     over x,y,z
c----------------------------------------------------------
      !implicit integer (i-n), real*8 (a-h,o-z)
      implicit none
      include 'param.i'
      include 'one.i' ! stores rho
      real*8 rho_dens_xyz, dtpoprho ! external
      real*8 x,y,z ! Input
      integer i ! Input
      real*8 rholoc, dtpop_drho
      real*8 dtpop_dx,dtpop_dy,dtpop_dz ! Output
      !--------------
	rholoc= rho_dens_xyz(x,y,z) ! get rho and derivs of rho
	dtpop_drho= dtpoprho(rholoc,i)
	dtpop_dx= dtpop_drho*drho_dx
	dtpop_dy= dtpop_drho*drho_dy
	dtpop_dz= dtpop_drho*drho_dz
      return
      end


c======================================================================
c======================================================================
	subroutine dvflow_xyz(x,y,z,i, dvflow_dx, dvflow_dy, dvflow_dz)
c----------------------------------------------------------
c     Yu.P. 2012
c
c     calculates the derivatives of vflow 
c     d(cm/sec)/d_rho 
c     of species 'i'
c     over x,y,z
c----------------------------------------------------------
      !implicit integer (i-n), real*8 (a-h,o-z)
      implicit none
      include 'param.i'
      include 'one.i' ! stores rho
      real*8 rho_dens_xyz, dvflwrho ! external
      real*8 x,y,z ! Input
      integer i ! Input
      real*8 rholoc, dvflow_drho
      real*8 dvflow_dx,dvflow_dy,dvflow_dz ! Output
      !--------------
	rholoc= rho_dens_xyz(x,y,z) ! get rho and derivs of rho
	dvflow_drho= dvflwrho(rholoc,i)
	dvflow_dx= dvflow_drho*drho_dx
	dvflow_dy= dvflow_drho*drho_dy
	dvflow_dz= dvflow_drho*drho_dz
      return
      end
      
c======================================================================
c======================================================================
      subroutine hotdervs_xyz(u,wf,dddcnx,dddcny,dddcnz,
     .dddx,dddy,dddz,dddw)
c---------------------------------------------------------
c     analytical calculation of the derivatives for the ray-tracing
c     equations from the non-relativistic hot plasma
c     with ions and electrons
c--------------------------------------------------------- 
      implicit none
cSm030226
      include 'param.i'
      include 'one.i'
c-----input
      double precision u(6)  ! x,y,z,Nx,Ny,Nz
      double precision wf    ! wave friquency
c-----output
      double precision dddcnz,dddcnx,dddcny,
     .dddz,dddx,dddy,dddw   !derivatives from D
c-----uses
      external dhot_sum,bxyz,rho_dens_xyz,tempe_xyz,
     .temperho,tpoprho,vflowrho,wcw,wpw_2
      double complex dhot_sum
      double precision bxyz,rho_dens_xyz,tempe_xyz,
     .temperho,tpoprho,vflowrho,wcw,wpw_2

c-----locals
      double complex dd(5,nbulka) !dD/d(X_s,Y_s,Tav_s,tpop_as,V_s)
      double precision t_kev,nll,nperp,step,cnp

      double complex K_sum(3,3),dK(3,3,7),ddnp_h,ddnll_h,ddnp,d
      double precision x,y,z,cnz,cnx,cny,rholoc
      integer i
      double precision dwpw2dx,dwpw2dy,dwpw2dz,dwcwdz,dwcwdx,dwcwdy, 
     .dtemp_dz,dtemp_dx,dtpop_dz,dtpop_dx,dvflow_dz,dvflow_dx,
     .dtemp_dy,dtpop_dy,dvflow_dy,
     .dxdwe,dydwe,
     .dnpdz,dnpdx,dnpdy,dnpdcnz,dnpdcnx,dnpdcny,
     .dnlldz,dnlldx,dnlldy,dnlldcnz,dnlldcnx,dnlldcny,
     .dnpdw,dnlldw

      double precision mass_ar(nbulka),x_ar(nbulka),y_ar(nbulka),
     .t_av_ar(nbulka),tpop_ar(nbulka),vflow_ar(nbulka)

      integer k      
            
      x=u(1)
      y=u(2)
      z=u(3)
      cnx=u(4)
      cny=u(5)
      cnz=u(6)

      bmod=bxyz(x,y,z) ! also finds components of B !->one.i
      rholoc=rho_dens_xyz(x,y,z)

      if(nbulk.gt.nbulka) then
          write(*,*)'in hotdervs_xyz: nbulk.gt.nbulka'
          write(*,*)'Set the value of parameter nbulka.ge.nbulk'
          write(*,*)'parameter nbulka is given in param.i in hotdervs'
          write(*,*)'nbulka,nbulk',nbulka,nbulk
          stop
      endif

c-----The initilization mass_ar
      call put_mass(mass_ar,nbulk)
              
      nll= (cnx*bx+cny*by+cnz*bz)*o_bmod ! o_bmod is 1/|B|  !== Npar
      cnp= dsqrt(cnx**2+cny**2+cnz**2 -nll**2)  !== Nperp
      nperp=cnp
      
	do i=1,nbulk                  
          x_ar(i)=wpw_2(x,y,z,i)
          y_ar(i)=wcw(x,y,z,i)
          if(i.eq.1) y_ar(1)=-y_ar(1)
          t_keV=temperho(rholoc,i)   !keV averaged temperature
          t_av_ar(i)=t_keV*1.d3    !eV
          tpop_ar(i)=tpoprho(rholoc,i)
          vflow_ar(i)=vflowrho(rholoc,i) !cm/sec
	enddo

c      write(*,*)'in hotdervs_xyz-1: x,y,z, Y=',x,y,z, y_ar(1)
	d=dhot_sum(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,
     .vflow_ar,nll,nperp,1,K_sum)

      call Ddhot(nbulk,mass_ar,x_ar,y_ar,t_av_ar,tpop_ar,vflow_ar
     .,nll,nperp,K_sum,dd,ddnp_h,ddnll_h,ddnp)
     
      call dnd_xyz(x,y,z, cnx,cny,cnz,
     .dnpdx,dnpdy,dnpdz,    dnpdcnx,dnpdcny,dnpdcnz,
     .dnlldx,dnlldy,dnlldz, dnlldcnx,dnlldcny,dnlldcnz,
     .dnpdw,dnlldw)
     
      dddx = ddnp_h*dnpdx+ddnll_h*dnlldx
      dddy = ddnp_h*dnpdy+ddnll_h*dnlldy
      dddz = ddnp_h*dnpdz+ddnll_h*dnlldz
      dddcnx=ddnp_h*dnpdcnx+ddnll_h*dnlldcnx
      dddcny=ddnp_h*dnpdcny+ddnll_h*dnlldcny
      dddcnz=ddnp_h*dnpdcnz+ddnll_h*dnlldcnz
      
      dddw=ddnp_h*dnpdw+ddnll_h*dnlldw

      do i=1,nbulk
        call dwpw_2(x,y,z,i,dwpw2dx,dwpw2dy,dwpw2dz) !derivs of (wp/w)^2
        call dwcw(x,y,z,i,  dwcwdx, dwcwdy, dwcwdz) ! derivs of wc/w
        if(i.eq.1) then ! negative wc/w for electrons
           dwcwdz=-dwcwdz
           dwcwdx=-dwcwdx
           dwcwdy=-dwcwdy
        endif 
        call dtemp_xyz(x,y,z,i,  dtemp_dx,  dtemp_dy,  dtemp_dz)
        call dtpop_xyz(x,y,z,i,  dtpop_dx,  dtpop_dy,  dtpop_dz)
        call dvflow_xyz(x,y,z,i, dvflow_dx, dvflow_dy, dvflow_dz)
        dddx= dddx+
     .       dd(1,i)*dwpw2dx + dd(2,i)*dwcwdx 
     +     + dd(3,i)*dtemp_dx*1.d3 
     .     + dd(4,i)*dtpop_dx + dd(5,i)*dvflow_dx
        dddy= dddy+
     .       dd(1,i)*dwpw2dy + dd(2,i)*dwcwdy 
     +     + dd(3,i)*dtemp_dy*1.d3 
     .     + dd(4,i)*dtpop_dy + dd(5,i)*dvflow_dy
        dddz= dddz+
     .       dd(1,i)*dwpw2dz + dd(2,i)*dwcwdz 
     +     + dd(3,i)*dtemp_dz*1.d3 
     .     + dd(4,i)*dtpop_dz + dd(5,i)*dvflow_dz
        ! Note: The terms like dtemp_dx may cause problem if T(rho) 
        ! changes too fast at the edge 
        ! (edge of flat T profile, for example).
        ! Try to set profiles to be as smooth as possible.
        dxdwe=-2.d0*x_ar(i)/wf
        dydwe=-y_ar(i)/wf
        dddw=dddw+dd(1,i)*dxdwe+dd(2,i)*dydwe
      enddo ! i

      return
      end  ! hotdervs_xyz end


c======================================================================
c======================================================================


c======================================================================
c======================================================================


c======================================================================
c======================================================================
      