
c*****************************************************
c*----------------subroutine GRAPHnpr----------------*
c* Prepaires file drawgenr.in for xdraw              *
c* INPUT: from common blocks'gr.i''n_parb.i' and   *
c         'write'                                    *
c*****************************************************
      SUBROUTINE GRAPHnpr
      implicit integer (i-n), real*8 (a-h,o-z)
      INTEGER I
      include 'param.i'
      include 'gr.i'
      include 'write.i'
      include 'n_parb.i'
      include 'three.i'

      OPEN(11,file='npar.bin',form='unformatted')

      DO 10 I=1,nrayelt
        WRITE(11) REAL(wthet(I)),REAL(wr(I)),REAL(wz(I)),
     1  REAL(ws(I)),REAL(rez(I)),REAL(spsi(I)),
     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(wmdevr(I)),
     3  REAL(wnpar(I)-wmdevr(I))
10    CONTINUE

      WRITE(11)

      CLOSE(11)
13	Format(16g15.5)
14	Format(14e11.4)
15	Format(12e12.4)

      OPEN(21,file='npar.sap')
      OPEN(23,file='nparadd.sap')
 1    format(' ',i3)
      write(21,*)'wthet wr wz ws psi npar nper npol
     1 bpol btor btot rex imey rez flux wmtor '
      write(23,*)'delpwr phi wx wy rmaga zmaga wnrho gnpar wxe wp 
     1wnparplb wnparmnb '
      DO 120 I=1,nrayelt
       WRITE(21,13) wthet(I),wr(I),wz(I),ws(I),spsi(I),
     2 wnpar(I),wnper(I),wnpol(I),
     3 bpoloid(I),btor(I),btot(I),
     4 dreal(cwexde(I)),dimag(cweyde(I)),dreal(cwezde(I)),fluxn(I),
     5 wmtor(I)
       WRITE(23,15) delpwr(I),wphi(I),wr(I)*cos(wphi(I)),
     1 wr(I)*sin(wphi(I)),100.d0*rma*dcos(wphi(I)),
     2 100.d0*rma*dsin(wphi(I)),
     3 wnrho(I),gnpar(I),wxe(I),wp(I),wnparplb(I),wnparmnb(I)
120    CONTINUE
      CLOSE(21)
      CLOSE(23)

c-------------------------------------------- 
c     5 contours psi=const
      OPEN(24,file='section.sap')
      write(24,*)
16    Format(10e12.4)
      do i=1,NP+1
        WRITE(24,16)AR(1,i),AZ(1,i),AR(2,i),AZ(2,i),AR(3,i),AZ(3,i),
     1              AR(4,i),AZ(4,i),AR(5,i),AZ(5,i)
      enddo
      CLOSE(24)
c-------------------------------------------- 
      RETURN
      END



c        ********************** prepebw***************************
c        *                      -----                             *
c        *  prepebw -subroutine to prepar the  data for studying   *
C        *  ebw                                                   *
c        **********************************************************
      subroutine prepebw(t,u,is)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'write.i'
      include 'three.i'
      include 'n_parb.i'
      double precision modr,modrold
      dimension u(6)
c---------------------------------------------
c     z,r (m)
      z=u(1)
      r=u(2)
      phi=u(3)
      cnz=u(4)
      cnr=u(5)
      cm=u(6)

c---------------------------------------
c     the creation of the array wthet(i): poloidal angle
c     along the ray trajectory (ir radians)
      zcomp=z-zma
      rcomp=r-rma
      modr=dsqrt(zcomp*zcomp+rcomp*rcomp)  !=|r|

      btor(is)=bphi
      btot(is)=bmod
      gradpsi(is)=dsqrt(dpdzd**2+dpdrd**2)
      gradpdr(is)=gradpsi(is)/r
c--------------------------------------------------------------------
c     e_theta=e_psi x e_phi		   (x is a vector production)
c     e_psi=(e_z*dpsidz+e_r*dpsidr)/sqrt(dpsidz**2+dpsidr**2)
c     e_theta=(e_z*dpsidr-e_r*dpsidz)/sqrt(dpsidz**2+dpsidr**2)
c     Npol=N * e_theta
c--------------------------------------------------------------------
      bpoloid(is)=(bz*dpdrd-br*dpdzd)/gradpsi(is)
      wnpol(is)=(cnz*dpdrd-cnr*dpdzd)/gradpsi(is)
c--------------------------------------------------------------------
      if (is.eq.1) then
         if(zcomp.ge.0.d0) then
            wthet(1)=dacos(rcomp/modr)
         else
            wthet(1)=-dacos(rcomp/modr)
         endif
      else

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

      wmdevr(is)=cm/r*bphi/bmod

c--------------------------------------------------------------------
      

      xe=x(z,r,phi,1)
      wxe(is)=dsqrt(xe)
c     wnrho is the radial component of the refructive index.
c     It is positive if it is directed along the gradient 
c     of the flux surface. 
      wnrho(is)=(cnr*dpdrd+cnz*dpdzd)/gradpsi(is)
      gnpar(is)=bpoloid(is)/btot(is)
      wp(is)=gnpar(is)*wxe(is)
     
      cnphi=cm/r
     
 10   continue
      return
      END

c*****************************************************
c*----------------subroutine GRAPHebw----------------*
c* Prepares for file drawgenr.in for xdraw           *
c* INPUT: from common blocks'gr.i''n_parb.i' and     *
c         'write'                                    *
c*****************************************************
      SUBROUTINE GRAPHebw
      implicit integer (i-n), real*8 (a-h,o-z)
      INTEGER I
      include 'param.i'
      include 'gr.i'
      include 'write.i'
      include 'n_parb.i'
      include 'three.i'

      OPEN(11,file='npar.bin',form='unformatted')

      DO 10 I=1,nrayelt
        WRITE(11) REAL(wthet(I)),REAL(wr(I)),REAL(wz(I)),
     1  REAL(ws(I)),REAL(rez(I)),REAL(spsi(I)),
     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(wmdevr(I)),
     3  REAL(wnpar(I)-wmdevr(I))
10    CONTINUE

      WRITE(11)

      CLOSE(11)
13	Format(16g15.5)
14	Format(14e11.4)
15	Format(12e12.4)

      OPEN(21,file='npar.sap')
      OPEN(23,file='nparadd.sap')
 1    format(' ',i3)
      write(21,*)'wthet wr wz ws psi npar nper npol
     1 bpol btor btot rex imey rez flux wmtor '
      write(23,*)'delpwr phi wx wy rmaga zmaga wnrho gnpar wxe wp 
     1wnparplb wnparmnb '
      DO 120 I=1,nrayelt
       WRITE(21,13) wthet(I),wr(I),wz(I),ws(I),spsi(I),
     2 wnpar(I),wnper(I),wnpol(I),
     3 bpoloid(I),btor(I),btot(I),
     4 dreal(cwexde(I)),dimag(cweyde(I)),dreal(cwezde(I)),fluxn(I),
     5 wmtor(I)
       WRITE(23,15) delpwr(I),wphi(I),wr(I)*cos(wphi(I)),
     1 wr(I)*sin(wphi(I)),100.d0*rma*dcos(wphi(I)),
     2 100.d0*rma*dsin(wphi(I)),
     3 wnrho(I),gnpar(I),wxe(I),wp(I),wnparplb(I),wnparmnb(I)
120    CONTINUE
      CLOSE(21)
      CLOSE(23)

c-------------------------------------------- 
c     5 contours psi=const
      OPEN(24,file='section.sap')
      write(24,*)
16    Format(10e12.4)
      do i=1,NP+1
        WRITE(24,16)AR(1,i),AZ(1,i),AR(2,i),AZ(2,i),AR(3,i),AZ(3,i),
     1              AR(4,i),AZ(4,i),AR(5,i),AZ(5,i)
      enddo
      CLOSE(24)
c-------------------------------------------- 
      RETURN
      END



