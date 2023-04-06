c        ************************ b ***************************     
c        this function calculates the components  of  the   
c        magnetic field bz,br,bphi, 
c        the absolute value of the total magnetic field bmod
c        and the derivatives of bz,br,bphi,bmod by rz,r,phi 
c        It calculates the rho,rho2,a for the programms     
c        which use the density profile :dense,dxdr,dxdz     
c	 						      
c        It controls if the ray point is inside the plasma  
c        if the ray point is outside the plasma then the    
c        calculations will be finished 		      
c        ******************************************************
c
c------------------------------------------------------------------
c								   !
c        input parameters					   !
c      								   !
c      z, r, phi - point which the components of magnetic field	   !
c                  are calculated in.				   !
c------------------------------------------------------------------
      double precision
     1function b(z,r,phi)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
      include 'six.i'
      real*8
     1 ias1r,ias2r_Sm
cSm070324
      real*8
     & thetapol, l_zr, l_z_b_r_b
c-----local
      real*8 sin_phi,cos_phi

      data ifirstc/1/
      save ifirstc

      write(*,*) 'b called'
      pause
      
      sin_phi=dsin(phi)
      cos_phi=dcos(phi)

c--------------------------------------------------------------
c     control that the point is inside the plasma
      epsbnd=1.d-10
      iboundb=-1
c-----------------------
      x_l=r-rma
      y_l=z-zma
      if ((x_l**2+y_l**2).lt.1.d-10) then ! close to mag.axis
         goto 10
      endif

c-----calculate poloidal angle  thetapol
      call theta_rz(x_l,y_l,theta_l) !0 =< theta_l < 2*pi

c-----calculate coordinates x_b,z_b at the limiter surface psi=psilim
c     and the poloidal angle thetapol
      call zr_psith(psilim,theta_l,z_b,r_b)
      
      l_zr=dsqrt(x_l**2+y_l**2)
      l_z_b_r_b=dsqrt((r_b-rma)**2+(z_b-zma)**2)

      if (l_zr.ge.l_z_b_r_b) then
         rho = l_zr / l_z_b_r_b         
         iboundb=3
      endif
      goto 10
c----------------------------------------------------
      if ((r.le.rmin+epsbnd).or.(r.ge.rmax-epsbnd)) then
         iboundb=2
         theta_l=thetapol(z,r) 
         call zr_psith(psilim,theta_l,z_b,r_b)
         rho=dsqrt(((r-rma)**2+(z-zma)**2)/((r_b-rma)**2+(z_b-zma)**2))
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
cSAP080816
      if ((z.ge.zp-epsbnd).or.(z.le.zm+epsbnd)) then 
         iboundb=1
cSm070324 
         theta_l=thetapol(z,r) 
         call zr_psith(psilim,theta_l,z_b,r_b)
         rho=dsqrt(((r-rma)**2+(z-zma)**2)/((r_b-rma)**2+(z_b-zma)**2))
         drhodzb=(z-zma)/(((r_b-rma)**2+(z_b-zma)**2)*rho)
         goto 10
      end if
c---------------------------------------------------------------------
 10     continue
c end of the control that the point is inside the plasma
c--------------------------------------------------------------
c
c nr4a , nz4a and nrya  are given in common five as parameters
c nreqda,nzeqda,nlimit are given in common four as parameters
c
      nr4=nx+4
      nz4=ny+4
      ncx=nr4
      ncy=nz4
      res=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,0,0,r,z,nr4a)
      dpsidr=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,1,0,r,z,nr4a)
      dpsidz=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,0,1,r,z,nr4a)
   
c--------------------------------------------------------------------
c     idx derivativs order 0.ge.idx.le.3
c--------------------------------------------------------------------
      nr4=nx+4  
      if (iboundb.ge.1) then
c        the point is outside the plasma
         idx=0
         ffr=ias1r(txf,nx,nr4,cx,idx,psilim)
         ffd=ffr
         dffr=0.d0
         dres=dffr
      else
c        the point is inside the plasma
         idx=0
         ffr=ias1r(txf,nx,nr4,cx,idx,res)
         ffd=ffr
         idx=1
         dffr=ias1r(txf,nx,nr4,cx,idx,res)
         dres=dffr
      endif
      psid=res
      dpdzd=dpsidz
      dpdrd=dpsidr
      pp=1.d0/r
      dpdph= 0.d0 ! dpsi/dphi ! Zero for now, generalize later
      dpdxd= cos_phi*dpdrd - (sin_phi*pp)*dpdph ! dpsi/dx ! Added for now
      dpdyd= sin_phi*dpdrd + (cos_phi*pp)*dpdph ! dpsi/dy ! Added for now
c------------------------------------------------------------------
c this  part of the programm calculates the rho,rho2 and a
c for the functions dense,dxdr,dxdz

       a=1.d0       
       if (iboundb.eq.-1) rho=rhopsi(res)
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
      bz=-dpdrd*pp
      br=dpdzd*pp
      bphi=ffd*pp
 
      if(toteqd.ge.0.d0) spol=-1.d0
      if(toteqd.lt.0.d0) spol=1.d0
      spol=-spol
      if(toteqd.eq.0.d0) then
           if (ifirstc.eq.1) then
              write(*,*)'current is not given in eqdsk (toteqd=0)' 
              ifirstc=0
           endif
c          The total current is not given. The sign of the B_pol will be
c          determined from psi, using the original poloidal flux from eqdsk.
c          In this case the poloidal flux had max on the magnetic axis
c          (and dpsimax=-1d0), or it had min (and dpsimax=1.d0).
           spol=dpsimax
      endif
      
      bz=bz*spol
      br=br*spol
      bx= br*cos_phi - bphi*sin_phi ! YuP added for now
      by= br*sin_phi + bphi*cos_phi ! YuP added for now

      dpdzr=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,1,1,r,z,nr4a)
      dpdrr=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,2,0,r,z,nr4a)
      dpdzz=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,0,2,r,z,nr4a)
      dpdzrd=dpdzr
      dpdrrd=dpdrr
      dpdzzd=dpdzz
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
      b=dsqrt(bphi*bphi+bz*bz+br*br)
c------------------------------------------------------------------
      ob=1.d0/b
ccc      dbmdz= ob*(bz*dbzdz + br*dbrdz + bphi*dbphdz) ! YuP: Commented-out; defined below
      dbmdr= ob*(bz*dbzdr + br*dbrdr + bphi*dbphdr)
      dbmdph=ob*(bz*dbzdph+ br*dbrdph+ bphi*dbphdph)
      
      dbmdx= ob*(bx*dbxdx + by*dbydx + bz*dbzdx) ! YuP added for now
      dbmdy= ob*(bx*dbxdy + by*dbydy + bz*dbzdy) ! YuP added for now
      dbmdz= ob*(bx*dbxdz + by*dbydz + bz*dbzdz) ! YuP added for now
c---------------------------------------------------------------------
      
      return
      end



