c        *********************BOUND****************************
c        *                        -                           *
c        * This function gives the ray position in            *
c        * the plasma (If the point of the ray is inside the  *
c        * limitter the BOUND=-1 , else BOUND=1),	      *
c        * and it reflects the ray from plasma boundary.      *
c        * It calculates the number of the reflections-irefl  *
c        * It gives the command to stop ray calculation	      *
c        * if irefl.ge.ireflm				      *
c        ******************************************************
c
c------------------------------------------------------------------
c								   !
c        Input parameters					   !
c      								   !
c      Z, R, PHI - point which the components of magnetic field	   !
c                  are calculated in.
c
c     the input parameter epsbnd(the small distance from the boundary)
c     is set inside this subroutine
c------------------------------------------------------------------!
c     output parameters:iflref=1 after reflection =-1 before refl.
c			cnzref,cnrref,cmref= N_z,N_r,M after reflectin point,or
c				     =N_z,N_r,M before reflection point
c                       ibound=1 point is outside the plasma =-1 in plasma
c----------------------------------------------------------------------
      subroutine bound(z,r,phi,cnz,cnr,cm,iflref,
     &z_ref,r_ref,phi_ref,cnzref,cnrref,cmref,
     &ibound,dzdt,drdt)
c      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none

c-----input
      real*8 
     &z,r,phi,     !space coordinates
     &cnz,cnr,cm,  !refractive index coordinates
     &dzdt,drdt    !wave group velocity in RHS of the ray-tracing equations
c-----output
      integer
     &iflref,     !iflref=1 after reflection =-1 before refl.
     &ibound      !ibound=1 point is outside the plasma =-1 in plasma

      real*8 
     &z_ref,r_ref,phi_ref, ! reflection point coordinates
     &cnzref,cnrref,cmref  ! N_z,N_r,M after reflectin point,or
                           ! =N_z,N_r,M before reflection point
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
cSm041006
      include 'rkutta.i'
c-----externals
      real*8
     &ias1r,ias2r,ias2r_Sm,b,rhopsi

c-----locals
      real*8
     &epsbnd,   !accuracy of point near the plasma edle LCFS
     &zp,zm,zzrm,zzrp,res,rrr,dmrdt,dpsidr,dpsidz,
     &edg,edag,grad

      integer nr4,nz4,ipx,ipx4,idx,
     &i_reflection  !=1 it was wall limiter reflection
                    !=0 it was not walllimiter reflection
      cnzref=cnz
      cnrref=cnr
      cmref=cm
      z_ref=z
      r_ref=r
      phi_ref=phi

      iflref=-1
      epsbnd=1.d-8
cSm070121
      epsbnd=1.d-7
cSm070613
      epsbnd=1.d-3
cSm070613
      epsbnd=1.d-5

      ibound=-1

cSAP080727
c      if (no_reflection.eq.1) goto 30
cSAP081028
      write(*,*)'in bound no_reflection',no_reflection
      
      if (no_reflection.eq.1) then 
          
c          write(*,*)'before wall_limiter_reflection_point'
c          write(*,*)'z,r,phi,cnz,cnr,cm',
c     &               z,r,phi,cnz,cnr,cm
          
          call wall_limiter_reflection_point(z,r,phi,cnz,cnr,cm,
     &    i_reflection,
     &    cnzref,cnrref,cmref,z_ref,r_ref,phi_ref)

          write(*,*)'after wall_limiter_reflection_point i_reflection=',
     &    i_reflection

          if (i_reflection.eq.1) then  

c             write(*,*)'z,r,phi,cnz,cnr,cm',
c     &                  z,r,phi,cnz,cnr,cm
c             write(*,*)'z_ref,r_ref,phi_ref',z_ref,r_ref,phi_ref
c             write(*,*)'cnzref,cnrref,cmref',
c     &                  cnzref,cnrref,cmref

             irefl=irefl+1
             ibound=1.d0
             iflref=1
          endif
          goto 30
      endif
c-----------------------
      if ((r.le.rmin+epsbnd).or.(r.ge.rmax-epsbnd)) then
	 write(*,*)'in bound rmin,r,rmax,epsbnd',rmin,r,rmax,epsbnd
         ibound=1.d0
         goto 10
      end if
      rrr=r
c---------------- idx derivativs order 0.ge.idx.le.3---------------
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
         write(*,*)'in bound r,zm,z,zp,epsbnd',r,zm,z,zp,epsbnd
         ibound=1.d0
         goto 10
      endif
c-----------------------
cSm030224
      nr4=nx+4
      nz4=ny+4

      ncx=nr4
      ncy=nz4
cSm030224
c      res=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,0,0,r,z)
      res=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,0,0,r,z,nr4a)

      if (res.ge.(psilim*(1.d0-epsbnd))) then
         ibound=1.d0
         write(*,*)'in bound res,psilim,epsbnd',res,psilim,epsbnd
         goto 10
      end if
      rho=rhopsi(res)
c-----------------------
      bmod=b(z,r,phi)
      if (rho.ge.1.d0-epsbnd) then
         ibound=1.d0
         write(*,*)'in bound rho,epsbnd',rho,epsbnd
         goto 10
      else
         goto 30
      end if
c-----------------------
 10   continue 

      write(*,*)'in bound no_reflection,i',no_reflection,'ibound',ibound
      
      if ((no_reflection.eq.1).and.(ibound.eq.1)) then 
          ibound=-1           !initialization
          write(*,*)'before wall_limiter_reflection_point'
          write(*,*)'z,r,phi,cnz,cnr,cm',
     &               z,r,phi,cnz,cnr,cm

          call wall_limiter_reflection_point(z,r,phi,cnz,cnr,cm,
     &    i_reflection,
     &    cnzref,cnrref,cmref,z_ref,r_ref,phi_ref)

          write(*,*)'after wall_limiter_reflection_point i_reflection=',
     &    i_reflection

          if (i_reflection.eq.1) then  

             write(*,*)'z,r,phi,cnz,cnr,cm',
     &                  z,r,phi,cnz,cnr,cm
             write(*,*)'z_ref,r_ref,phi_ref',z_ref,r_ref,phi_ref
             write(*,*)'cnzref,cnrref,cmref',
     &                  cnzref,cnrref,cmref

             irefl=irefl+1

             write(*,*)'###irefl=',irefl

             ibound=1.d0
             iflref=1
          endif
          goto 30
      endif

c-------------------------------------------------------
 100  continue
      ncx=nr4
      ncy=nz4
cSm030224
c      dpsidr=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,1,0,r,z)
c      dpsidz=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,0,1,r,z)
      dpsidr=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,1,0,r,z,nr4a)
      dpsidz=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,0,1,r,z,nr4a)

      edg=dpsidz*dzdt+dpsidr*drdt

cSm041006
      edg=edg*prmt(3)/dabs(prmt(3))

c      write(*,*)'in bound dpsidr,dpsidz',dpsidr,dpsidz
c      write(*,*)'in bound drdt,dzdt',drdt,dzdt,'edg',edg
      grad=dsqrt(dpsidz*dpsidz+dpsidr*dpsidr)
      dmrdt=dsqrt(drdt*drdt+dzdt*dzdt)
      edag=edg/(grad*dmrdt)
      if (edg.gt.0d0) then
         irefl=irefl+1
         write(*,*)'in bound irefl=',irefl
         call reflect(z,r,phi,cnz,cnr,cnzref,cnrref)
	 iflref=1
      else
	 grad=dsqrt(dpsidz*dpsidz+dpsidr*dpsidr)
	 dmrdt=dsqrt(drdt*drdt+dzdt*dzdt)
         write(*,*)'in bound edg.le.0,edg,grad,dmrdt',edg,grad,dmrdt
      endif
 30   continue

      write(*,*)'in bound after 30'

      return
      end

c        *********************BOUNDC correction****************
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
c      Z, R, - point which the components of magnetic field	   !
c                  are calculated in.
c
c     the input parameter epsbnd(the small distance from the boundary)
c     is set inside this subroutine
c------------------------------------------------------------------!
c        Output parameter:iboundc=1 point is outside the  plasma   !				   !
c                         iboundc=-1 point is inside the  plasma   !				   !
c------------------------------------------------------------------!
      subroutine boundc(z,r,iboundc)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
      double precision
     1 ias1r,ias2r,ias2r_Sm
     
      iboundc=-1
      epsbnd=  1.d-3 !1.d-10 ! YuP: changed to 1e-3

      if (no_reflection.eq.1) goto 10 !-> return/end
c-----------------------
      ! Check that the ray is within rz-grid
      if ( (z.le.zeqmin+epsbnd).or.(z.ge.zeqmax-epsbnd) .or.
     +     (r.le.rmin+epsbnd)  .or.(r.ge.rmax-epsbnd)   ) then
         write(*,*)'in boundc r,z,epsbnd',r,z,epsbnd
         iboundc=1
         goto 10
      end if
      rrr=r
c---------------- idx derivativs order 0.ge.idx.le.3---------------
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
         write(*,*)'in boundc rrr,zm,z,zp',rrr,zm,z,zp
         write(*,*)'in boundc zm+epsbnd,z,zp-epsbnd',
     1	 zm+epsbnd,z,zp-epsbnd
         iboundc=1
         goto 10  !-> return/end
      end if
c---------------------------------------------------------------------
      nr4=nx+4
      nz4=ny+4
      ncx=nr4
      ncy=nz4
c      res=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,0,0,r,z)
      res=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,0,0,r,z,nr4a)
c      rho=rhopsi(res)
      if (res.gt.(psilim*(1.d0-epsbnd))) then
         write(*,*)'in boundc res,psilim*(1-epsbnd)',
     1	 res,psilim*(1.d0-epsbnd)
         iboundc=1
         goto 10 !-> return/end
      end if
c-----------------------
      rho=rhopsi(res)
      if (rho.gt.1.d0-epsbnd) then
         write(*,*)'in boundc rho,1.d0-epsbnd',rho,1.d0-epsbnd
         iboundc=1
         goto 10 !-> return/end
      end if
c---------------------------------------------------------------------
 10     continue

      return
      end
      
      
      
      

c        *********************REFLECT**************************
c        *                        -                           *
c        * This subroutine reflects                           *
c        * the ray from plasma boundary	                      *
c        ******************************************************
c
c------------------------------------------------------------------
c								   !
c        Input parameters					   !
c      								   !
c      Z, R, PHI - point which the components of magnetic field	   !
c                  are calculated in.				   !
c      cnz,cnr						   !
c------------------------------------------------------------------!
c     output parameters:					   !
c	       	cnzref,cnrref= N_z,N_r after reflectin point	   !
c----------------------------------------------------------------------
      subroutine reflect(z,r,phi,cnz,cnr,cnzref,cnrref)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'three.i'
      include 'five.i'
      double precision
     1 ias1r,ias2r,ias2r_Sm

      nr4=nx+4
      nz4=ny+4
      ncx=nr4
      ncy=nz4
c       dpsidz=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,0,1,r,z)
c       dpsidr=ias2r(tx,nx,ty,ny,cxy,ncx,ncy,1,0,r,z)
      dpsidz=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,0,1,r,z,nr4a)
      dpsidr=ias2r_Sm(tx,nx,ty,ny,cxy,ncx,ncy,1,0,r,z,nr4a)

      dpsimd=dsqrt(dpsidz**2+dpsidr**2)
      pp=1.0d0/dpsimd
      cnpsi=-pp*(dpsidz*cnz+dpsidr*cnr)
c      cnteta= pp*(-dpsidz*cnr+dpsidr*cnz)
c      cnzref=pp*(cnpsi*dpsidz+cnteta*dpsidr)
c      cnrref=pp*(cnpsi*dpsidr-cnteta*dpsidz)
              write(*,*)'reflect1 cnzref,cnrref',cnzref,cnrref
      ! Unit vector along grad(psi):
      grad_psi_z= dpsidz*pp
      grad_psi_r= dpsidr*pp
              write(*,*)'grad_psi_z,grad_psi_r',grad_psi_z,grad_psi_r
              write(*,*)'cnpsi,cnz,cnr',cnpsi,cnz,cnr
      ! Components of Npsi, along grad(psi)              
      cnz_psi=-grad_psi_z*cnpsi
      cnr_psi=-grad_psi_r*cnpsi
              write(*,*)'cnz_psi,cnr_psi',cnz_psi,cnr_psi
      cnzref=cnz-2.d0*cnz_psi
      cnrref=cnr-2.d0*cnr_psi
              write(*,*)'reflect2 cnzref,cnrref',cnzref,cnrref
c	      cnpsi=-pp*(dpsidz*cnzref+dpsidr*cnrref)
c	      cnteta=pp*(-dpsidz*cnrref+dpsidr*cnzref)
              write(*,*)'reflection point'
              !pause
      return
      end


