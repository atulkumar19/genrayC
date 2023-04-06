
c the data for the spline coefficients for psi(r,z),feqd(psi),pres(psi),
c limiters :z_plus(r) and z_minus(r)
c They are created in equilib.f

c      parameter (nr4a=nreqda+4,nz4a=nzeqda+4,nrya=nz4a)
c      parameter (nlim4=nlimit+4)
c      parameter nrya=max(nreqda,nzeqda)+4

      real*8      rmax,zrmax,rmin,zrmin,zmax,rzmax,zmin,rzmin,
     1            tx,ty,txf,cx,cy,
     4            trlimp,trlimm,
     5            cxlimp,cxlimm,cxy,
     +		  tpres,cpres
      integer     nx,ny,ncx,ncy,ip,im
      common/five/rmax,zrmax,rmin,zrmin,zmax,rzmax,zmin,rzmin,
     1            tx(nr4a),ty(nz4a),txf(nr4a),cx(nr4a),cy(nrya),
     4            trlimp(nlim4),trlimm(nlim4),
     5            cxlimp(nlim4),cxlimm(nlim4),cxy(nr4a,nz4a),
     +		  tpres(nr4a),cpres(nr4a),
     6            nx,ny,ncx,ncy,ip,im
c
c  nreqda must be >= nreqd	number of the mech points along r in eqdsk.file
c  nzeqda must be >= nzeqd	number of the mech points along z in eqdsk.file
c

c   coordinates of Lackner rectangle:
c   rmax,zrmax,rmin,zrmin,zmax,rzmax,zmin,rzmin

c   coefficients for the spline approximation of
c   psi(r,z): tx(nr4a),ty(nz4a),cxy (cy is a working array)
c   feqd(psi):txf,cx
c   pres(psi):tpres(nr4a),cpres(nr4a)
c   z_limmiter plus(r) or (z_above): trlimp,cxlimp
c   z_limmiter minus(r) or (z_under): trlimm,cxlimm


c==============================================================
c Arrays for the spline coefficients for dengrid(x,y)
c   parameter (nx4a=nxdena+4,ny4a=nydena+4,nxya=ny4a)
c   coefficients for the spline approximation of dengrid(i,j):
c   tx_den,ty_den,cxy_den (cy_den is a working array)

      real*8   tx_den, ty_den, cy_den, cxy_den
      integer  nx_den, ny_den, ncx_den, ncy_den
      
      common/den_spline/
     1         tx_den(nx4a),ty_den(ny4a),cy_den(nxya),
     1         cxy_den(nx4a,ny4a),
     1         nx_den, ny_den, ncx_den, ncy_den
     
