c        ********************** sdr1  ***********************
c        *                      -----                       *
c        * this subroutine calculates                       *
c        * value u(i0) as  a solution of the        	    *
c        * dispersion relation for isolv=2	            *
c        * The hamiltonian derivatives  are calculated	    *
c        * numerically or analytically   		    *
c        ****************************************************
c
c-------------------------------------------------------------------
c								    !
c        input parameters					    !
c      								    !
c                                  		         	    !
c      u - solution of geometrical optics equations at the point t  !
c	u(1) = z                                                    !
c	u(2) = r                                                    !
c	u(3) = phi                                                  !
c	u(4) = n_z                                                  !
c	u(5) = n_r                                                  !
c	u(6) = r*n_phi                                              !
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
c      								    !
c        output parameters					    !
c	u(1) = z                                                    !
c	u(2) = r                                                    !
c	u(3) = phi                                                  !
c	u(4) = n_z                                                  !
c	u(5) = n_r                                                  !
c	u(6) = r*n_phi                                              !
c                                                                   !
c-------------------------------------------------------------------
c     this program uses following functions and subroutines	    !
c     b,gamma1,dxdz,dxdr,dxdphi,dydz,dydr,dydphi,x,y,s,		    !
c-------------------------------------------------------------------
	double precision
     1function sdr1(u,deru,i0,epsten)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      dimension u(*),deru(*)
      write(*,*)'in sdr1 u',u(1),u(2)
         z=u(1)
         r=u(2)
         phi=u(3)
         cnz=u(4)
         cnr=u(5)
         cm=u(6)
c---------------------------------------------------------------
c u(i0) is the value from previous time step
c---------------------------------------------------------------
c ind=1 newton method is used for the solution of the dispersion
c relation
c---------------------------------------------------------------
c the solution of the dispersion relation
c by newton iterative process
c x^n+1=x^n-d(x^n)/dddx(x^n)
c---------------------------------------------------------------
 1       continue
c-------------------for analytical derivatives-----------
c	 idif=1
c-------------------for numerical derivatives------------
c	 idif=2
c------------------------------------------------------------------
      if (idif.eq.1) goto 1955

c----------------------------------------------------------------
c numerical hamiltonian derivatives calculations
c----------------------------------------------------------------
           write(*,*)'in sdr1 numerical deriv'
           write(*,*)'in sdr1 z,r,phi',z,r,phi
	   bmod=b(z,r,phi)
	   gam=gamma1(z,r,phi,cnz,cnr,cm)
c          write(*,*)'in sdr1 z=',z,'r=',r,'phi=',phi
c          write(*,*)'in sdr1 cnz=',cnz,'cnr=',cnr,'cm=',cm
	   ddd=hamilt1(z,r,phi,cnz,cnr,cm)
           write(*,*)'in sdr1 ddd=',ddd,'epsten=',epsten
             if (dabs(ddd).lt.epsten) then
	       goto 60
	     endif
           step=0.0000001d0
	   hz=step
	   hr=step
	   hphi=step
	   hnz=step
	   hnr=step
	   hm=step
c---------------------------------------------------------------
       if (i0.eq.4) then         
	   cnzplus=cnz+hnz
	   cnzmins=cnz-hnz
           write(*,*)'in 1 z,r',z,r
	   bmod=b(z,r,phi)
	   gam=gamma1(z,r,phi,cnzplus,cnr,cm)
	   hp=hamilt1(z,r,phi,cnzplus,cnr,cm)
	   gam=gamma1(z,r,phi,cnzmins,cnr,cm)
	   hm=hamilt1(z,r,phi,cnzmins,cnr,cm)
	   dddcnz=(hp-
     1  	   hm)/(2.*hnz)
c           deru(1)=dddcnz
            goto 1953
       endif ! i0=4

       if (i0.eq.5) then         
	   cnrplus=cnr+hnr
	   cnrmins=cnr-hnr
	   gam=gamma1(z,r,phi,cnz,cnrplus,cm)
	   hp=hamilt1(z,r,phi,cnz,cnrplus,cm)
	   gam=gamma1(z,r,phi,cnz,cnrmins,cm)
	   hm=hamilt1(z,r,phi,cnz,cnrmins,cm)
	   dddcnr=(hp-
     2  	   hm)/(2.*hnr)
c           deru(2)=dddcnr
            goto 1953
       endif ! i0=5

       if (i0.eq.6) then         
	   cmplus=cm+hm
	   cmminus=cm-hm
	   gam=gamma1(z,r,phi,cnz,cnr,cmplus)
	   hp=hamilt1(z,r,phi,cnz,cnr,cmplus)
	   gam=gamma1(z,r,phi,cnz,cnr,cmminus)
	   hmin=hamilt1(z,r,phi,cnz,cnr,cmminus)
	   dddcm=(hp-
     3  	   hmin)/(2.*hm)
c           deru(3)=dddcm
            goto 1953
       endif ! i0=6

       if (i0.eq.1) then         
	   zplus=z+hz
	   zminus=z-hz
           write(*,*)'in 2 zplus,r',z,r
	   bmod=b(zplus,r,phi)
	   gam=gamma1(zplus,r,phi,cnz,cnr,cm)
	   hp=hamilt1(zplus,r,phi,cnz,cnr,cm)
           write(*,*)'in 2 zminus,r',zminus,r
	   bmod=b(zminus,r,phi)
	   gam=gamma1(zminus,r,phi,cnz,cnr,cm)
	   hm=hamilt1(zminus,r,phi,cnz,cnr,cm)
	   dddz=(hp-
     4  	 hm)/(2.*hz)
c           deru(4)=-dddz
            goto 1953
       endif ! i0=1

       if (i0.eq.2) then         
	   rplus=r+hr
	   rminus=r-hr
           write(*,*)'in 3 z,rplus',z,rplus
	   bmod=b(z,rplus,phi)
	   gam=gamma1(z,rplus,phi,cnz,cnr,cm)
	   hp=hamilt1(z,rplus,phi,cnz,cnr,cm)
           write(*,*)'in 3 z,rminus',z,rminus
	   bmod=b(z,rminus,phi)
	   gam=gamma1(z,rminus,phi,cnz,cnr,cm)
	   hm=hamilt1(z,rminus,phi,cnz,cnr,cm)
	   dddr=(hp-
     5  	   hm)/(2.*hr)
c           deru(5)=-dddr
            goto 1953
       endif ! i0=2

       if (i0.eq.3) then         
	   phipls=phi+hphi
	   phimin=phi-hphi
           write(*,*)'in 4 z,r',z,r
	   bmod=b(z,r,phipls)
	   gam=gamma1(z,r,phipls,cnz,cnr,cm)
	   hp=hamilt1(z,r,phipls,cnz,cnr,cm)
	   bmod=b(z,r,phimin)
	   gam=gamma1(z,r,phimin,cnz,cnr,cm)
	   hm=hamilt1(z,r,phimin,cnz,cnr,cm)
	   dddphi=(hp-
     6  	   hm)/(2.*hphi)
c           deru(6)=-dddphi
            goto 1953
       endif ! i0=3
          write(*,*)'numerical rigth side calculations deru(i)='
	   write(*,*)deru(1),deru(2),deru(3),deru(4),deru(5),deru(6)
	   go to 1953
c----------------------------------------------------------------
c  end of numerical calculation of hamiltonian derivatives
c------------------------------------------------------------------
 1955	  continue
c-----------------------------------------------------------------
c analytical calculations of hamiltonian derivatives
c-----------------------------------------------------------------
	   bmod=b(z,r,phi)
	   gam=gamma1(z,r,phi,cnz,cnr,cm)
           ds=dsin(gam)
           dc=dcos(gam)
             ds2=ds*ds
             dc2=dc*dc
                cnt=cn(r,cnz,cnr,cm)
                  cn2=cnt*cnt
c-------------------------------------------------------------
          if (id.eq.3) then
c
c         Appltoh-Hartry dispersion relation
c
          dxidz=dxdz(z,r,phi,1)
          dxidr=dxdr(z,r,phi,1)
          dxidph=dxdphi(z,r,phi,1)
          dyidz=dydz(z,r,phi,1)
          dyidr=dydr(z,r,phi,1)
          dyidph=dydphi(z,r,phi,1)

          xi=x(z,r,phi,1)
          yi=y(z,r,phi,1)

            py2=yi*yi
            py3=py2*yi
            py4=py2*py2
            px=1.-xi
            px2=px*px
            ds4=ds2*ds2

              det=py4*ds4+4.*py2*px2*dc2
              zer=0.d0
              if (det.lt.zer) then
               write(*,*)'det in rside less then 0 det=',det
	       stop
              end if
              sqrdet=dsqrt(det)
              p1=0.5/sqrdet
                ddetdx=-8.*py2*px*dc2
                ddetdy=4.*py3*ds4+8.*yi*px2*dc2
                pz=2.*px-py2*ds2+ioxm*sqrdet
c---------------------------------------------------------------------
c                cnt=cn(r,cnz,cnr,cm)
c                  cn2=cnt*cnt
c                  cn2od=1.-2.*xi*px/(2.*px-py2*ds2+sqrdet)
c                  cn2ex=1.-2.*xi*px/(2.*px-py2*ds2-sqrdet)
c                  write(*,*)'cn2od=',cn2od,'cn2ex=',cn2ex
c                  write(*,*)'cn2=',cn2
             ddd=cn2-(1.-2.*xi*px/pz)
             if (dabs(ddd).lt.epsten) then
	       goto 60
	     endif
c------------------------------------------------------------------
                pz2=1./(pz*pz)
                  dpzdx=-2.+ioxm*p1*ddetdx
                  dpzdy=-2.*yi*ds2+ioxm*p1*ddetdy
                   dddx=-2.*((2.*xi-1.)*pz+xi*px*dpzdx)*pz2
                   dddy=-2.*xi*px*dpzdy*pz2

                    dddc2=-2.*xi*px*pz2*(py2+ioxm*(-2.*py4*(1.-dc2)+
     1                    4.*py2*px2)*p1)
c-----------------------------------------------------------------
                      dddz=dddx*dxidz+dddy*dyidz+dddc2*dc2dz
                      dddr=dddx*dxidr+dddy*dyidr+dddc2*dc2dr-
     1                     2*cm*cm/(r**3)
                      dddphi=dddx*dxidph+dddy*dyidph+dddc2*dc2dph

                      dddcnz=2.*cnz+dddc2*dc2dnz
                      dddcnr=2.*cnr+dddc2*dc2dnr
                      dddcm=2.*cm/(r*r)+dddc2*dc2dm
c--------------------------------------------------------------------
            goto 50
           end if
c   end if Appelton - Hartree conditions
c-----------------------------------------------------------------
c  ---------------------------------------------------------------
c  cold plasma dispersion relation with electrons and ions
c  x(i=1),y(i=1) electrons component
c  x(i may be=2,nbulk),y(i may be=2,nbulk) ions components
c  ib number of component for which delib=1-yib may be equal
c  zero inside the plasma,ib  may be=from 1 to nbulk
c  dispersion relation is multiplied by delib
c
	   pc=1.+dc2
         call s(z,r,phi,s1,s2,s3,s4,s6,s7)
         xib=x(z,r,phi,ib)
         yib=y(z,r,phi,ib)
           delib=1.-yib

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
         peym=xe/(1.-ye)
	 peyp=xe/(1.+ye)
	   a0e=-peyp*ds2
	   a1e=s7*ds2+s4*dc2
	   b0e=s4*peyp*pc+xe*(s6-peyp)*ds2
	   b1e=-s4*s7*pc-s3*(s6-peyp)*ds2
	   c0e=-xe*s4*(s6-peyp)
	   c1e=s4*s3*(s6-peyp)
	   dele=1.-ye
	     ad=dele*a1e+a0e
	     bd=dele*b1e+b0e
	     cd=dele*c1e+c0e
c----------------------------------------------------------------------
	     da1edc=-s7+s4
	     da0edc=peyp
	     db1edc=-s7*s4+s3*(s6-peyp)
	     db0edc=s4*peyp-xe*(s6-peyp)
	     dc1edc=0.d0
	     dc0edc=0.d0

	     dadc2=dele*da1edc+da0edc
	     dbdc2=dele*db1edc+db0edc
	     dcdc2=dele*dc1edc+dc0edc
c-------------------------------------------------------------------
      end if
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
       xe=x(z,r,phi,1)
       ye=y(z,r,phi,1)
         ppe=1./(1.-ye)
         peym=xe/(1.-ye)
	 peyp=xe/(1.+ye)
	 peym2=peyp*ppe

         ppb=1./(1.-yib)
         pbym=xib/(1.-yib)
	 pbyp=xib/(1.+yib)
	 pbym2=pbyp*ppb

	   a0b=-pbyp*ds2
	   a1b=(s1-peym2)*ds2+s4*dc2
	   b0b=s4*pbyp*pc+xib*(s3-peym)*ds2
	   b1b=-s4*(s1-peym2)*pc-(s2-peyp)*(s3-peym)*ds2
	   c0b=-xib*s4*(s3-peym)
	   c1b=s4*(s2-peyp)*(s3-peym)
	   delib=1.-yib
	     ad=delib*a1b+a0b
	     bd=delib*b1b+b0b
	     cd=delib*c1b+c0b
c--------------------------------------------------------------------
	     da1bdc=-(s1-peym2)+s4
	     da0bdc=pbyp
	     db1bdc=-s4*(s1-peym2)+(s2-peyp)*(s3-peym)
	     db0bdc=s4*pbyp-xib*(s3-peym)
	     dc1bdc=0.d0
	     dc0bdc=0.d0

	     dadc2=delib*da1bdc+da0bdc
	     dbdc2=delib*db1bdc+db0bdc
 	     dcdc2=delib*dc1bdc+dc0bdc
c-------------------------------------------- -----------0-----------
      end if
c  if 3 end
c-----------------------------------------------------------------

        if (id.eq.1) then
c
c  dispersion relation a*n**4+b*n**2+c=0
c
          ddd=ad*cn2*cn2+bd*cn2+cd
             if (dabs(ddd).lt.epsten) then
	       goto 60
	     end if
	end if

        if (id.eq.2) then
c
c  dispersion relation n**2-(-b+ioxm*sqrt(b**2-4*a*c))/(2*a)=0
c
         ddd=cn2-(-bd+ioxm*dsqrt(bd*bd-4.*ad*cd))/(2.*ad)
             if (dabs(ddd).lt.epsten) then
	       goto 60
	     end if
	end if
c------------------------------------------------------------------
c  derivative calculation
c   dadz,dadr,dadphi
c   dbdz,dbdr,dbdphi
c   dcdz,dcdr,dcdphi
c----------------------------------------------------------------
       dadz=0.d0
       dadr=0.d0
       dadphi=0.d0
       dbdz=0.d0
       dbdr=0.d0
       dbdphi=0.d0
       dcdz=0.d0
       dcdr=0.d0
       dcdphi=0.d0

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

         dxidz=dxdz(z,r,phi,i)
         dxidr=dxdr(z,r,phi,i)
         dxidph=dxdphi(z,r,phi,i)
         dyidz=dydz(z,r,phi,i)
         dyidr=dydr(z,r,phi,i)
         dyidph=dydphi(z,r,phi,i)

         if ( i.eq.ib) goto 20
c
c  i.ne.ib
c
	 xi=x(z,r,phi,i)
         yi=y(z,r,phi,i)
	   pyim=1./(1.-yi)
	   pyip=1./(1.+yi)
	   pyim2=pyim*pyim

c if 5 begin
	     if (ib.gt.1) then
	       ds1dxi=pyim*pyip
	       ds2dxi=pyim
	       ds1dyi=2.*xi*yi*ds1dxi*ds1dxi
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
	       ds7dyi=2.*xi*yi*ds7dxi*ds7dxi
	     end if
c if 6 end

c if 7 begin
	     if (ib.eq.1) then

               da1exi=-ds7dxi*ds2-ds4dxi*dc2
	       da0exi=0.d0
	       db1exi=(ds4dxi*s7+ds7dxi*s4)*pc+
     1		      (ds3dxi*(s6-xe/(1.+ye))+ds6dxi*s3)*ds2
               db0exi=-ds4dxi*xe/(1.+ye)*pc-
     1                 xe*ds6dxi*ds2
               dc1exi=-ds4dxi*s3*(s6-xe/(1.+ye))-
     1                 s4*ds3dxi*(s6-xe/(1.+ye))-s4*s3*ds6dxi
               dc0exi=xe*ds4dxi*(s6-xe/(1.+ye))+xe*s4*ds6dxi

	       da1eyi=-ds7dyi*ds2-ds4dyi*dc2
	       da0eyi=0.d0
	       db1eyi=(ds4dyi*s7+ds7dyi*s4)*pc+
     1                (ds3dyi*(s6-xe/(1.+ye))+ds6dyi*s3)*ds2
               db0eyi=-ds4dyi*xe/(1.+ye)*pc-xe*ds6dyi*ds2
	       dc1eyi=-ds4dyi*s3*(s6-xe/(1.+ye))-
     1    	       ds3dyi*s4*(s6-xe/(1.+ye))-ds6dyi*s4*s3
               dc0eyi=xe*ds4dyi*(s6-xe/(1.+ye))+ds6dyi*xe*s4

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
	       db1bxi=(ds4dxi*(s1-xe/(1.-ye*ye))+ds1dxi*s4)*pc+
     1		      (ds2dxi*(s3-xe/(1.-ye))+ds3dxi*(s2-xe/(1.+ye)))*
     2                 ds2
               db0bxi=-ds4dxi*xib/(1.+yib)*pc-
     1                 xib*ds3dxi*ds2
               dc1bxi=-ds4dxi*(s2-xe/(1.+ye))*(s3-xe/(1.-ye))-
     1                 ds2dxi*s4*(s3-xe/(1.-ye))-
     2                 ds3dxi*s4*(s2-xe/(1.+ye))
               dc0bxi=xib*ds4dxi*(s3-xe/(1.-ye))+xib*s4*ds3dxi

	       da1byi=-ds1dyi*ds2-ds4dyi*dc2
	       da0byi=0.d0
	       db1byi=(ds4dyi*(s1-xe/(1.-ye*ye))+s4*ds1dyi)*pc+
     1                (ds2dyi*(s3-xe/(1.-ye))+ds3dyi*(s2-xe/(1.+ye)))*
     2                 ds2
               db0byi=-ds4dyi*xib/(1.+yib)*pc-xib*ds3dyi*ds2
	       dc1byi=s4*(-ds2dyi*(s3-xe/(1.-ye))-
     1    	       ds3dyi*(s2-xe/(1.+ye)))
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
 20     continue
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
             ds3dxb=1./(1.+yib)
	     ds4dxb=1.d0
             ds3dyb=-xib/(1.+yib)**2
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
               da0bxb=-1./(1.+yib)*ds2
               db1bxb=(ds4dxb*(s1-xe/(1.-ye*ye))+ds1dxb*s4)*pc+
     1                (ds2dxb*(s3-xe/(1.-ye))+ds3dxb*(s2-xe/(1.+ye)))*
     2                 ds2
               db0bxb=(-ds4dxb*xib/(1.+yib)+s4/(1.+yib))*pc+
     1                ((s3-xe/(1.-ye))-xib*ds3dxb)*ds2
               dc1bxb=-ds4dxb*(s2-xe/(1.+ye))*(s3-xe/(1.-ye))-
     1                 ds2dxb*s4*(s3-xe/(1.-ye))-
     2                 ds3dxb*s4*(s2-xe/(1.+ye))
               dc0bxb=-s4*(s3-xe/(1.-ye))+xib*ds4dxb*(s3-xe/(1.-ye))+
     +                xib*s4*ds3dxb

               da1byb=0.d0
               da0byb=xib/(1.+yib)**2*ds2
               db1byb=(ds4dyb*(s1-xe/(1.-ye*ye))+s4*ds1dyb)*pc+
     1                (ds2dyb*(s3-xe/(1.-ye))+ds3dyb*(s2-xe/(1.+ye)))*
     2                 ds2
               db0byb=(-ds4dyb*xib/(1.+yib)-s4*xib/(1.+yib)**2)*pc-
     -                xib*ds3dyb*ds2
               dc1byb=s4*(-ds2dyb*(s3-xe/(1.-ye))-
     1                 ds3dyb*(s2-xe/(1.+ye)))
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

          dadz=dadz+dadxi*dxidz+dadyi*dyidz
          dadr=dadr+dadxi*dxidr+dadyi*dyidr
          dadphi=dadphi+dadxi*dxidph+dadyi*dyidph

          dbdz=dbdz+dbdxi*dxidz+dbdyi*dyidz
          dbdr=dbdr+dbdxi*dxidr+dbdyi*dyidr
          dbdphi=dbdphi+dbdxi*dxidph+dbdyi*dyidph

          dcdz=dcdz+dcdxi*dxidz+dcdyi*dyidz
          dcdr=dcdr+dcdxi*dxidr+dcdyi*dyidr
          dcdphi=dcdphi+dcdxi*dxidph+dcdyi*dyidph

 10       continue
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
               da1bxe=-1./(1.-ye*ye)*ds2-ds4dxe*dc2
               da0bxe=0.d0
               db1bxe=(ds4dxe*(s1-xe/(1.-ye*ye))+s4/(1.-ye*ye))*pc+
     1                (1./(1.+ye)*(s3-xe/(1.-ye))+
     +                 1./(1.-ye)*(s2-xe/(1.+ye)))*ds2
               db0bxe=-ds4dxe*xib/(1.+yib)*pc-xib/(1.-ye)*ds2
               dc1bxe=s4*(-1./(1.+ye)*(s3-xe/(1.-ye))-
     -                     1./(1.-ye)*(s2-xe/(1.+ye)))-
     -                 ds4dxe*((s2-xe/(1.+ye))*(s3-xe/(1.-ye)))
               dc0bxe=xib*s4/(1.-ye)+xib*ds4dxe*(s3-xe/(1.-ye))

               da1bye=-2.*xe*ye/(1.-ye*ye)**2*ds2-ds4dye*dc2
               da0bye=0.d0
               db1bye=s4*2.*xe*ye/(1.-ye*ye)**2*pc+
     +                (-xe/(1.+ye)**2*(s3-xe/(1.-ye))+
     +                  xe/(1.-ye)**2*(s2-xe/(1.+ye)))*ds2
               db0bye=-xib*xe/(1.-ye)**2*ds2
               dc1bye=s4*(xe/(1.+ye)**2*(s3-xe/(1.-ye))-
     -                    xe/(1.-ye)**2*(s2-xe/(1.+ye)))
               dc0bye=xib*s4*xe/(1.-ye)**2

                 dadxi=delib*da1bxe+da0bxe
                 dbdxi=delib*db1bxe+db0bxe
                 dcdxi=delib*dc1bxe+dc0bxe
                 dadyi=delib*da1bye+da0bye
                 dbdyi=delib*db1bye+db0bye
                 dcdyi=delib*dc1bye+dc0bye

               goto 40
	     end if
c if 10 end

c if 11 begin
c
c   ib = 1, i = 1
c
             if (ib.eq.1) then

               da1exe=-ds4dxe*dc2
               da0exe=-1./(1.+ye)*ds2
               db1exe=ds4dxe*s7*pc+s3/(1.+ye)*ds2
               db0exe=(-ds4dxe*xe/(1.+ye)+s4/(1.+ye))*pc+
     +                 (s6-2.*xe/(1.+ye))*ds2
               dc1exe=-ds4dxe*s3*(s6-xe/(1.+ye))-
     -                s4*s3/(1.+ye)
               dc0exe=-s4*(s6-2.*xe/(1.+ye))+xe*(s6-xe/(1.+ye))

               da1eye=0.d0
               da0eye=xe/(1.+ye)**2*ds2
               db1eye=-s3*xe/(1.+ye)**2*ds2
               db0eye=-s4*xe/(1.+ye)**2*pc+xe*xe/(1.+ye)**2*ds2
               dc1eye=-ds4dye*s3*(s6-xe/(1.+ye))+s4*s3*xe/(1.+ye)**2
               dc0eye=-xe*xe*s4/(1.+ye)**2

                 dadxi=dele*da1exe+da0exe
                 dbdxi=dele*db1exe+db0exe
                 dcdxi=dele*dc1exe+dc0exe
                 dadyi=dele*da1eye+da0eye-a1e
                 dbdyi=dele*db1eye+db0eye-b1e
                 dcdyi=dele*dc1eye+dc0eye-c1e


	     end if
c if 11 end

 40     continue
         dxidz=dxdz(z,r,phi,1)
         dxidr=dxdr(z,r,phi,1)
         dxidph=dxdphi(z,r,phi,1)
         dyidz=dydz(z,r,phi,1)
         dyidr=dydr(z,r,phi,1)
         dyidph=dydphi(z,r,phi,1)
c------------------------------------------------------------------
        dadz=dadz+dadxi*dxidz+dadyi*dyidz+dadc2*dc2dz
        dadr=dadr+dadxi*dxidr+dadyi*dyidr+dadc2*dc2dr
        dadphi=dadphi+dadxi*dxidph+dadyi*dyidph+dadc2*dc2dph

        dbdz=dbdz+dbdxi*dxidz+dbdyi*dyidz+dbdc2*dc2dz
        dbdr=dbdr+dbdxi*dxidr+dbdyi*dyidr+dbdc2*dc2dr
        dbdphi=dbdphi+dbdxi*dxidph+dbdyi*dyidph+dbdc2*dc2dph

        dcdz=dcdz+dcdxi*dxidz+dcdyi*dyidz+dcdc2*dc2dz
        dcdr=dcdr+dcdxi*dxidr+dcdyi*dyidr+dcdc2*dc2dr
        dcdphi=dcdphi+dcdxi*dxidph+dcdyi*dyidph+dcdc2*dc2dph

c---------------------------------------------------------------------
        dcn=cn(r,cnz,cnr,cm)
        dcn2=dcn*dcn
        dcn4=dcn2*dcn2

	dadnz=dadc2*dc2dnz
	dbdnz=dbdc2*dc2dnz
	dcdnz=dcdc2*dc2dnz

	dadnr=dadc2*dc2dnr
	dbdnr=dbdc2*dc2dnr
	dcdnr=dcdc2*dc2dnr

	dadm=dadc2*dc2dm
	dbdm=dbdc2*dc2dm
	dcdm=dcdc2*dc2dm
c------------------------------------------------------------------

        if (id.eq.1) then

c
c  dispersion relation a*n**4+b*n**2+c=0
c
          dddz=dcn4*dadz+dcn2*dbdz+dcdz
          dddr=dcn4*dadr+dcn2*dbdr+dcdr-
     1	       cm*cm*(4.*dcn2*ad+2.*bd)/(r**3)
          dddphi=dcn4*dadphi+dcn2*dbdphi+dcdphi

         dddn=4.*ad*dcn2+2.*bd
          dddcnz=dddn*cnz+dcn4*dadnz+dcn2*dbdnz+dcdnz
          dddcnr=dddn*cnr+dcn4*dadnr+dcn2*dbdnr+dcdnr
          dddcm=dddn*cm/(r*r)+dcn4*dadm+dcn2*dbdm+dcdm
        end if

        if (id.eq.2) then

c
c  dispersion relation n**2-(-b+ioxm*sqrt(b**2-4*a*c))/(2*a)=0
c
         det=dsqrt(bd*bd-4.*ad*cd)
         p4=0.5/(ad*ad)
         p5=1./det
         p6=-bd+ioxm*det
c------------------------------------------------------------------
         p7=ioxm*(bd*dbdz-2.*dadz*cd-2.*ad*dcdz)*p5
         dddz=-((-dbdz+p7)*ad-dadz*p6)*p4
         p7=ioxm*(bd*dbdr-2.*dadr*cd-2.*ad*dcdr)*p5
         dddr=-((-dbdr+p7)*ad-dadr*p6)*p4-
     1		2.*cm**2/(r**3)
         p7=ioxm*(bd*dbdphi-2.*dadphi*cd-2.*ad*dcdphi)*p5
         dddphi=-((-dbdphi+p7)*ad-dadphi*p6)*p4

         p7=ioxm*(bd*dbdnz-2.*dadnz*cd-2.*ad*dcdnz)*p5
          dddcnz=2.*cnz-
     1      	 ((-dbdnz+p7)*ad-dadnz*p6)*p4
         p7=ioxm*(bd*dbdnr-2.*dadnr*cd-2.*ad*dcdnr)*p5
          dddcnr=2.*cnr-
     1      	 ((-dbdnr+p7)*ad-dadnr*p6)*p4
         p7=ioxm*(bd*dbdm-2.*dadm*cd-2.*ad*dcdm)*p5
          dddcm=2.*cm/(r*r)-
     1      	 ((-dbdm+p7)*ad-dadm*p6)*p4
c--------------------------------------------------------------------
        end if
 50        continue
 1953   continue
           deru(4)=dddcnz
           deru(5)=dddcnr
           deru(6)=dddcm
           deru(1)=dddz
           deru(2)=dddr
           deru(3)=dddphi
	   write(*,*)'in sdr1 i0=',i0
	   write(*,*)'in sdr1 u(i0)=',u(i0),'ddd=',ddd
     	   write(*,*)'sdr1 deru(i0)=',deru(i0)
	   u(i0)=u(i0)-ddd/deru(i0)
	   write(*,*)'in sdr1 u(i0)=',u(i0)
	 if(i0.eq.1) then
           z=u(1)
	 endif
	 if(i0.eq.2) then
           r=u(2)
	 endif
         phi=u(3)
	 if(i0.eq.4) then
           cnz=u(4)
	 endif
	 if(i0.eq.5) then
           cnr=u(5)
	 endif
         cm=u(6)
	   goto 1
 60       continue
 	  sdr1=u(i0)
          return
          end
