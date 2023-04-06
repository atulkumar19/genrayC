
      subroutine netcdf_create
C===> YuP 120505 
c     Create a netCDF file. Write only title, version, mnemonic.
c     This file will be opened several times to add data.

      implicit none
      include 'param.i'
      include 'writencdf.i'
      include 'one.i'
      include 'ions.i'
      include 'fourb.i'
      include 'five.i'
      include 'cone_nml.i'     !nccone
      include 'grill_nml.i'    !ngrill

c --- include file for netCDF declarations 
c --- (obtained from NetCDF distribution)
      include 'netcdf.inc'

c --- some stuff for netCDF file ---
      integer ncid,vid,istatus,nccre2,ncvdef2,ncdid2,ncddef2
      integer length_char
      integer nraysid,neltmaxid,twoid,char64id,char128id,char512id,
     &char8id

      character ltitle*512

      ncid=nccre2(filenc,NCCLOB,istatus)
      call check_err(istatus)
      
      
c     Brief description added to file:
      ltitle='netCDF file of ray data from GENRAY version: '//version
      if( length_char(ltitle).gt.512 ) stop 'Adjust ltitle in netcdfrw2'
      call ncaptc2(ncid,NCGLOBAL,'title',NCCHAR,length_char(ltitle),
     +     ltitle,istatus)

c--------------------------
c     Genray version
c--------------------------
      char64id=ncddef2(ncid,'char64dim',64,istatus)
      vid=ncvdef2(ncid,'version',NCCHAR,1,char64id,istatus)
      write(*,*)'netcdf_create: after version:'
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,21,
     +          'GENRAY version number',istatus)

c--------------------------
c     Mnemonic for the run
c--------------------------
      char128id=ncddef2(ncid,'char128dim',128,istatus)
      vid=ncvdef2(ncid,'mnemonic',NCCHAR,1,char128id,istatus)
      write(*,*)'netcdf_create: after mnemonic:'
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23,
     +          'Mnemonic run identifier',istatus)

      char512id=ncddef2(ncid,'char512dim',512,istatus)
      vid=ncvdef2(ncid,'eqdskin',NCCHAR,1,char512id,istatus)
      write(*,*)'netcdf_create: after eqdskin:'
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,39,
     +          'Name of input eqdsk, for eqsource=eqdsk',istatus)
    
c     Close netCDF file
      call ncclos3(ncid,istatus) 
      write(*,*)'netcdf_create: after ncclos3:'
      call check_err(istatus)
      ! it will be re-opened later for writing data
      
      write(*,*)'netcdf file is created'
     
      return
      end



C=======================================================================
C=======================================================================

      subroutine wrtnetcdf_plasma_prof(netcdfnml)
C===> YuP 120505
C     Write plasma profiles into existing netcdf file.     
C     These data are not related to ray-tracing results. 
c     Called by genray.f (line~424), just once.
      !implicit integer (i-n), real*8 (a-h,o-z)
      implicit none
      include 'param.i'     
      include 'one.i'
      include 'three.i'
      include 'fourb.i'
      include 'five.i' ! contains rmax
      include 'ions.i'     
      include 'onetwo.i'  !For profiles and total current.
      include 'rho.i'     !For areatot,voltot,torftot,totlength
      include 'writencdf.i'     !
      include 'eps.i'  ! reps(3,3)
      
      
c --- include file for netCDF declarations 
c --- (obtained from NetCDF distribution)
      include 'netcdf.inc'

      include 'nperpcom.i'
                                 
c --- some stuff for netCDF file ---
      integer ncid,vid,istatus,istat,ncvdef2,ncdid2,ncddef2
      integer nrhoid,nrhomid,nbulkid,nbulkmid
      integer NRgrid_id, NZgrid_id
      integer nxeqd_id,nyeqd_id,nzeqd_id
      integer i,rdims(2),rdimss(2)
      integer start(2),counts(2)
      integer xzdims(2),countxz(2)
      integer yzdims(2),countyz(2)
      integer xydims(2),countxy(2)
      integer nscanid,nparid,nperid
      integer npdims(2),countnp(2)
      integer dddims(3),startddd(3),countddd(3)
      
      character(*) netcdfnml ! input filename
      
      real*8 dense_xyz,zeffrho,temperho,wpw_2,wcw,bxyz,hamilt_xyz,
     +       tempe_xyz,tpoprho,vflowrho
      double complex hotnp
      real*8 tem1(NR*nbulk), r_max,x,y,z,den,dstep,xe,ye
      real*8 xmin,xmax
      integer kk,j, ixeq,iyeq,izeq, id_loc,ibw_loc
            
      real*8 dnpar,dnper,cnper2m,cnper2p,cnper_0,
     + cnt2,arg, cnpar, 
     + cnparp,cnperp,te,cnprim,dDcold
c-----local
      real*8 x_ar(nbulka),y_ar(nbulka),t_av_ar(nbulka),
     & tpop_ar(nbulka),vflow_ar(nbulka)
      real*8 tempiar(nbulka), cnper, dens_e,temp_e,z_eff,
     + cnprim_cl,cnprim_e,cnprim_i,cnprim_s(nbulka),ckvipl_s(nbulka)
      double complex dK(3,3,7),dd(5),d,ddn,aK(3,3)
      double complex dd5(5,nbulka),ddnp_h,ddnll_h,ddnp
      double complex disp_func,dcold_rlt, dhot, dhot_rlt

      double precision nll ! N_parallel
      double complex K(3,3), hotnpc
      integer iraystop,is,i_fkin
      
      real*8 f000, rho_larm_max0
      
      ! For i_save_disp=1 option:
      ! (saving Nper^2(xscan,Npar) for cold plasma and hot plasma)
      ! Use allocatable arrays 
      !(could take large memory space, if xscan-grid is refined,
      ! and/or inper_save_disp is large)
      ! These are real*4 arrays, to save memory:
      REAL(4),ALLOCATABLE :: xscan(:)    ! (ixscan_save_disp)
      REAL(4),ALLOCATABLE :: rhoscan(:)  ! (ixscan_save_disp)
      REAL(4),ALLOCATABLE :: wcwscan(:)  ! (ixscan_save_disp)
      REAL(4),ALLOCATABLE :: wpwscan(:)  ! (ixscan_save_disp)
      REAL(4),ALLOCATABLE :: wuwscan(:)  ! (ixscan_save_disp)
      REAL(4),ALLOCATABLE :: Npara(:)   ! (inpar_save_disp)
      REAL(4),ALLOCATABLE :: Npera(:)   ! (inper_save_disp)
      REAL(4),ALLOCATABLE :: Nperp2m(:,:)   !(ixscan_save_disp, inpar_save_disp)
      REAL(4),ALLOCATABLE :: Nperp2p(:,:)   !(ixscan_save_disp, inpar_save_disp)
      REAL(4),ALLOCATABLE :: Nperp_im(:,:)  !(ixscan_save_disp, inpar_save_disp)
      REAL(4),ALLOCATABLE :: Nperp2hot(:,:) !(ixscan_save_disp, inpar_save_disp)
      REAL(4),ALLOCATABLE :: ddd(:,:,:) 
      !(ixscan_save_disp, inper_save_disp, inpar_save_disp)
      

      data start/1,1/, startddd/1,1,1/
      
      f000=frqncy

c     1.1 open existing netCDF file: netcdfnml, 
c         and define additional;
c         dimensions,variables and attributes
c
c.......................................................................
c     1.1.1 open existing netCDF file: netcdfnml
      write(*,*)'wrtnetcdf_plasma_prof: Opening file'
      call ncopn2(netcdfnml,NCWRITE,ncid,istatus) ! Open existing file
      call check_err(istatus)
c.......................................................................
c     1.1.2 
c     Put Open NetCDF file into Define Mode
      write(*,*)'wrtnetcdf_plasma_prof: Define mode'
      call ncredf3(ncid,istatus)
      call check_err(istatus)

c-------------------------------------------------------------------
c     1.1.3 define dimensions
c     integer function ncddef(ncid,character(*) dim_name,
c                             integer cdim_siz, integer error_code)
c     returns dimension id.

c-----define dimensions      
      nbulkid=ncddef2(ncid,'nbulk',nbulk,istatus) 
      call check_err(istatus)     
      nrhoid=ncddef2(ncid,'nrho',NR,istatus)
      call check_err(istatus)     
      nrhomid=ncddef2(ncid,'nrhom',NR-1,istatus)
      call check_err(istatus)     
      rdimss(1)=nrhoid
      rdimss(2)=nbulkid
      counts(1)=NR
      counts(2)=nbulk
      
      ! YuP[Nov-2014] Added: Saving Local absorbed power over (R,Z) grid.
      NRgrid_id=ncddef2(ncid,'NRgrid',NRgrid,istatus)
      call check_err(istatus)     
      NZgrid_id=ncddef2(ncid,'NZgrid',NZgrid,istatus)
      call check_err(istatus)
      

      ! Get dimension ID from dimension name
      nxeqd_id=ncdid2(ncid,'nxeqd',istatus) !defined in netcdf_eqdsk_data
      nyeqd_id=ncdid2(ncid,'nyeqd',istatus) !defined in netcdf_eqdsk_data
      nzeqd_id=ncdid2(ncid,'nzeqd',istatus) !defined in netcdf_eqdsk_data
      xzdims(1)=nxeqd_id 
      xzdims(2)=nzeqd_id 
      countxz(1)=nxeqd
      countxz(2)=nzeqd
      yzdims(1)=nyeqd_id 
      yzdims(2)=nzeqd_id 
      countyz(1)=nyeqd
      countyz(2)=nzeqd
      xydims(1)=nxeqd_id 
      xydims(2)=nyeqd_id 
      countxy(1)=nxeqd
      countxy(2)=nyeqd
      
      
      ! For i_save_disp=1 option:
      nscanid=ncddef2(ncid,'nscan',ixscan_save_disp,istatus)
      call check_err(istatus)     
      nparid=ncddef2(ncid,'inpar_save_disp',inpar_save_disp,istatus)
      call check_err(istatus)     
      nperid=ncddef2(ncid,'inper_save_disp',inper_save_disp,istatus)
      call check_err(istatus)     
      npdims(1)=nscanid
      npdims(2)=nparid
      countnp(1)=ixscan_save_disp !size(xscan)
      countnp(2)=inpar_save_disp
      dddims(1)=nscanid
      dddims(2)=nperid
      dddims(3)=nparid
      countddd(1)=ixscan_save_disp !size(xscan)
      countddd(2)=inper_save_disp
      countddd(3)=inpar_save_disp
      
      vid=ncvdef2(ncid,'nbulk',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,43,
     +          'Number of Maxl plasma cmpts, electrons+ions',istatus)
      call check_err(istatus)  
         
      vid=ncvdef2(ncid,'freqcy',NCDOUBLE,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,14,
     +            'Wave frequency',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,2,
     +           'Hz',istatus)

      vid=ncvdef2(ncid,'rs_frc',NCDOUBLE,0,0,istatus) ! for model_b=4
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,39,
     +            'For model_b=4: rs_frc Separatrix radius',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,1,
     +           'm',istatus)

      vid=ncvdef2(ncid,'model_b',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,7,
     +            'model_b',istatus)

      vid=ncvdef2(ncid,'model_rho_dens',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,14,
     +            'model_rho_dens',istatus)

      vid=ncvdef2(ncid,'y_save_disp',NCDOUBLE,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +            'Y-coord for scans of wc/w vs x',istatus)
     
      vid=ncvdef2(ncid,'z_save_disp',NCDOUBLE,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +            'Z-coord for scans of wc/w vs x',istatus)

      vid=ncvdef2(ncid,'dmas',NCDOUBLE,1,nbulkid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,41,
     +           'plasma species mass: electrons, then ions',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,27,
     +           'Normalized to electron mass',istatus)

      vid=ncvdef2(ncid,'charge',NCDOUBLE,1,nbulkid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,44,
     +           'plasma species charge: electrons, then ions',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,31,
     +           'Normalized to electronic charge',istatus)

      vid=ncvdef2(ncid,'rho_bin',NCDOUBLE,1,nrhoid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,38,
     +            'normalized small radius bin boundaries',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'rho_bin_center',NCDOUBLE,1,nrhomid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,35,
     +            'normalized small radius bin centers',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'densprof',NCDOUBLE,2,rdimss,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,44,
     +     'plasma density at bin boundaries, e and ions',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,14,
     +     'particles/cm^3',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'temprof',NCDOUBLE,2,rdimss,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,49,
     +     'plasma temperatures at bin boundaries, e and ions',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +     'keV',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'zefprof',NCDOUBLE,1,rdimss(1),istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,41,
     +     'plasma Zeff at bin boundaries, e and ions',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,8,
     +     'unitless',istatus)
      call check_err(istatus)

      !--------------------------------------------------

      vid=ncvdef2(ncid,'bmodprofxz',NCDOUBLE,2,xzdims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,24,
     +     'Total B on xz-grid (y~0)',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +     'Tesla',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'densprofxz',NCDOUBLE,2,xzdims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,42,
     +     'plasma density on xz-grid (y~0), electrons',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,14,
     +     'particles/cm^3',istatus)
      call check_err(istatus)
      
      vid=ncvdef2(ncid,'densprofyz',NCDOUBLE,2,yzdims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,42,
     +     'plasma density on yz-grid (x~0), electrons',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,14,
     +     'particles/cm^3',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'densprofxy',NCDOUBLE,2,xydims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,42,
     +     'plasma density on xy-grid (z~0), electrons',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,14,
     +     'particles/cm^3',istatus)
      call check_err(istatus)


      !--------------------------------------------------
      ! plasma profiles VS x [m], at y=0,z=0
      !--------------------------------------------------
      vid=ncvdef2(ncid,'w_x_densprof_nc',NCDOUBLE,1,nrhoid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,26,
     +            'x-grid for plasma profiles',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,2,
     +     'm ',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'w_bmod_vs_x_nc',NCDOUBLE,1,nrhoid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,32,
     +            'B profile vs x, at fixed y and z',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,2,
     +     'T ',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'w_dens_vs_x_nc',NCDOUBLE,2,rdimss,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,46,
     +     'plasma density at bin bndries, e and ions vs x',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,13,
     +     'particles/cm^3',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'w_temp_vs_x_nc',NCDOUBLE,2,rdimss,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,51,
     +    'plasma temperatures at bin bndries, e and ions vs x',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +     'keV',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'w_zeff_vs_x_nc',NCDOUBLE,1,nrhoid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,31,
     +     'plasma Zeff at bin bndries vs x',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,8,
     +     'unitless',istatus)
      call check_err(istatus)

      !--------------------------------------------------
      ! plasma profiles VS y [m], at x=0,z=0
      !--------------------------------------------------
      vid=ncvdef2(ncid,'w_y_densprof_nc',NCDOUBLE,1,nrhoid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,26,
     +            'y-grid for plasma profiles',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,2,
     +     'm ',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'w_bmod_vs_y_nc',NCDOUBLE,1,nrhoid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,32,
     +            'B profile vs y, at fixed x and z',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,2,
     +     'T ',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'w_dens_vs_y_nc',NCDOUBLE,2,rdimss,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,46,
     +     'plasma density at bin bndries, e and ions vs y',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,13,
     +     'particles/cm^3',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'w_temp_vs_y_nc',NCDOUBLE,2,rdimss,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,51,
     +    'plasma temperatures at bin bndries, e and ions vs y',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +     'keV',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'w_zeff_vs_y_nc',NCDOUBLE,1,nrhoid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,31,
     +     'plasma Zeff at bin bndries vs y',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,8,
     +     'unitless',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'xscan',NCFLOAT,1,nscanid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,42,
     +            'refined x-grid for Nperp2(x,Npar) profiles',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,2,
     +     'm ',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'rhoscan',NCFLOAT,1,nscanid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23,
     +            'rho vs x (at fixed y,z)',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'wcwscan',NCFLOAT,1,nscanid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23,
     +            'omega_ce/omega profiles',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'wpwscan',NCFLOAT,1,nscanid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23,
     +            'omega_pe/omega profiles',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'wuwscan',NCFLOAT,1,nscanid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,27,
     +            'Upper_Hybrid/omega profiles',istatus)
      call check_err(istatus)

      !--------------------------------------------------
      ! For i_save_disp=1 option:
      if(i_save_disp.eq.1)then
      vid=ncvdef2(ncid,'Npara',NCFLOAT,1,nparid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,36,
     +            'Npar grid for Nperp2(x,Npar) profiles',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'Npera',NCFLOAT,1,nperid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,41,
     +            'Nper grid for ddd(x,Nper,Npar) dispersion',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'Nperp2m',NCFLOAT,2,npdims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,51,
     +    'Nperp2(x,Npar) profile vs x; Cold root with ioxm=-1',istatus)
      call check_err(istatus)
      
      vid=ncvdef2(ncid,'Nperp2p',NCFLOAT,2,npdims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,51,
     +    'Nperp2(x,Npar) profile vs x; Cold root with ioxm=+1',istatus)
      call check_err(istatus)
      
      vid=ncvdef2(ncid,'Nperp_im',NCFLOAT,2,npdims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,39,
     +    'imaginary Nperp (damping) profile vs x',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'Nperp2hot',NCFLOAT,2,npdims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,44,
     +   'Nperp2(x,Npar) profile vs x; Hot plasma root',istatus)
      call check_err(istatus)
      
      vid=ncvdef2(ncid,'ddd',NCFLOAT,3,dddims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,52,
     +   'ddd(x,Nper,Npar) dispers.function over (x,Nper) grid',istatus)
      call check_err(istatus)
      endif ! i_save_disp=1
c--------------------------------------------------------------------
      
cl    1.1.5 end the define-mode and start the data-mode
      call ncendf3(ncid,istatus)
      call check_err(istatus)
      write(*,*)'wrtnetcdf_plasma_prof: End of Define mode'



      ! Calculate data

      do i=1,NR
         zefprof(i)=zeffrho(rho_bin(i))
      enddo

c-----evaluate radial profiles of density, temperature, zeff
c-----to put with radial profile data.
c-----(units: 10**13/cm**3, keV)
      ! Get density profile as a function of x;
      ! Scan along x view line;  y=0, z=0 plane.
      ! (modify if needed; used for netcdf file only)
      !r_max= xeqmax ! Limits: from equilibr.grid
c      if(rlim>0.) r_max=rlim
      r_max=rmax ! rmax is set in dinitr(), for all cases of model_b
      z= 0.
      y= 0.
      dstep= (rmax-rmin)/(NR-1)
      do kk=1,nbulk
         do i=1,NR
            x= rmin+dstep*(i-1) ! [rmin; rmax]
            ! This is only used for saving data into nc-file:
            !(Not plotted anymore in genray_plot.py, so can be removed?)
            densprof(i,kk)=dense_xyz(x,y,z,kk)*1.e13 !-> get rho
            !The above rho should be same as rho_bin(i).
            !But if x coord. is not following rho_bin(0):rho_bin(NR)
            !then they are different. In this case do not save densprof:
            if (abs(rho-rho_bin(i)) .gt. 1.d-3) then
               densprof(i,kk)=0.d0
               ! (Better nothing than wrong data.)
            endif
            ! But temprof is a function of rho_bin, not x :
            temprof(i,kk)=temperho(rho_bin(i),kk) !can be plotted vs rho_bin
            !write(*,*)'x,rho,n=',x,rho,densprof(i,kk)
         enddo
      enddo
      write(*,*)'wrtnetcdf_plasma_prof: after 1D prof. densprof,temprof'

      !--------------------------------------------------
      ! Total B on the grid [Tesla]
      allocate(  bmodprofxz(nxeqd,nzeqd) ,STAT=istatus )
      if(istatus.eq.0)then
        call bcast( bmodprofxz, 0.d0, SIZE(bmodprofxz) )
      else
        stop 'Cannot allocate bmodprofxz'
      endif

      ! Density on xz-grid and yz-grid (for electrons only, for now)
      ! This is just for netcdf output file.
c     for output of electron density profile in xz-plane to .nc file
      allocate(  densprofxz(nxeqd,nzeqd) ,STAT=istatus )
      if(istatus.eq.0)then
        call bcast( densprofxz, 0.d0, SIZE(densprofxz) )
      else
        stop 'Cannot allocate densprofxz'
      endif
      
c     for output of electron density profile in yz-plane to .nc file
      allocate(  densprofyz(nyeqd,nzeqd) ,STAT=istatus )
      if(istatus.eq.0)then
        call bcast( densprofyz, 0.d0, SIZE(densprofyz) )
      else
        stop 'Cannot allocate densprofyz'
      endif

c     for output of electron density profile in xy-plane to .nc file
      allocate(  densprofxy(nxeqd,nyeqd) ,STAT=istatus )
      if(istatus.eq.0)then
        call bcast( densprofxy, 0.d0, SIZE(densprofxy) )
      else
        stop 'Cannot allocate densprofxy'
      endif

      write(*,*)'wrtnetcdf_plasma_prof: working on 2D prof. for plots..'
      kk=1 ! save for electrons !!!
      iyeq=nyeqd/2 
      y=yeq(iyeq)   ! approximately y~0 (phi~0)
      do izeq=1,nzeqd
         z=zeq(izeq)
      do ixeq=1,nxeqd
         x=xeq(ixeq)
         !uses rho which is found in dense_xyz
         densprofxz(ixeq,izeq)=dense_xyz(x,y,z,kk)*1.e13 
         bmodprofxz(ixeq,izeq)=bxyz(x,y,z)
      enddo
      enddo
      write(*,*)'wrtnetcdf_plasma_prof: after 2D prof. densprofxz'

      kk=1 ! save for electrons !!!      
      dstep= (zeqmax-zeqmin)/(nzeqd-1)
      izeq= (zmideqd-zeqmin)/dstep +1
      z= zeq(izeq)   ! approximately z~midplane 
      do ixeq=1,nxeqd
         x=xeq(ixeq)
      do iyeq=1,nyeqd
         y=yeq(iyeq)
         !dense_xyz is either analytical or 2D-(x,y)-spline
         densprofxy(ixeq,iyeq)=dense_xyz(x,y,z,kk)*1.e13 
         !densprofxy(ixeq,iyeq)=dengrid(ixeq,iyeq)*1.e13 
      enddo
      enddo
      write(*,*)'wrtnetcdf_plasma_prof: after 2D prof. densprofxy'
      
      kk=nbulk ! save for last species: Can be fast ions !!!      
      ixeq=nxeqd/2 
      x=xeq(ixeq)   ! approximately x~0 (phi~pi/2)
      do izeq=1,nzeqd
         z=zeq(izeq)
      do iyeq=1,nyeqd
         y=yeq(iyeq)
         !uses rho which is found in dense_xyz
         densprofyz(iyeq,izeq)=dense_xyz(x,y,z,kk)*1.e13 
      enddo
      enddo
      write(*,*)'wrtnetcdf_plasma_prof: after 2D prof. densprofyz'
      
c     Calculate 1D arrays of density,temperature,zeff
c     vs x-coord. and vs y-coord., at given fixed z
      call plasma_profiles_vs_xy(0.d0)  ! YuP: maybe from zma?
      write(*,*)'wrtnetcdf_plasma_prof: after 1D w_dens_vs_x_nc,...'

      ! This is used For i_save_disp=1 option,
      ! for calc. of dispersion function along (xscan,y=fixed,z=fixed).
      ! But find/save omega_ce/omega as a func. of xscan anyway,
      ! even when i_save_disp=0.  (Does not take much time or memory) 
      allocate(xscan(ixscan_save_disp),   STAT=istat)
      if(istat.ne.0) then
         print*,'Allocation problem for xscan()',istat
      endif
      xscan=0.d0 ! Initialize
      allocate(wcwscan(ixscan_save_disp), STAT=istat)
      if(istat.ne.0) then
         print*,'Allocation problem for wcwscan()',istat
      endif
      wcwscan=0.d0 ! Initialize omega_ce/omega
      allocate(wpwscan(ixscan_save_disp), STAT=istat)
      if(istat.ne.0) then
         print*,'Allocation problem for wpwscan()',istat
      endif
      wcwscan=0.d0 ! Initialize Upper_hybrid/omega
      allocate(wuwscan(ixscan_save_disp), STAT=istat)
      if(istat.ne.0) then
         print*,'Allocation problem for wuwscan()',istat
      endif
      wuwscan=0.d0 ! Initialize Upper_hybrid/omega
      allocate(rhoscan(ixscan_save_disp), STAT=istat)
      if(istat.ne.0) then
         print*,'Allocation problem for rhoscan()',istat
      endif
      rhoscan=0.d0 ! Initialize rhoscan( )
      if(i_save_disp.eq.1) then
         xmax= min(xeqmax,xmax_save_disp) ! To stay within equilib grid
         xmin= max(xeqmin,xmin_save_disp) ! To stay within equilib grid
      else ! i_save_disp=0 (no saving of Nperp(x) data)
         xmax= xeqmax
         xmin= xeqmin
      endif
      write(*,*) 'wrtnetcdf_plasma_prof: xmin,xmax=',xmin,xmax

      dstep= (xmax-xmin)/(countnp(1)-1)
      y= y_save_disp ! use default value, which is 0.
      z= z_save_disp ! use default value, which is 0.
      do i=1,countnp(1) ! size of xscan grid
         x= xmax -dstep*(i-1) ! scan from xmax to xmin
         xscan(i)= x ! scan along x-axis
         bmod=bxyz(x,y,z) !-> get b (needed for wcw)
         ye= wcw(x,y,z,1)   ! Ye== omega_ce/omega
         xe= wpw_2(x,y,z,1) ! Xe== (wpe/w)^2    ! also gets rho
         rhoscan(i)= rho ! as a function of xscan, for a fixed (y,z)
         wcwscan(i)= ye
         wpwscan(i)= sqrt(xe) ! wpe/w = omega_pe/omega
         wuwscan(i)= sqrt(xe+ye*ye) ! (Upper hybrid)/omega
               write(*,'(a,2e12.3,a,e12.3,a,e12.3,a,e12.3)') 
     +         'x,z[m]=', x,z,
     +         '  (w/wpe)^2=',  1.0/xe,
     +         '   w/wce=',     1.0/ye,
     +         '   bmod=',      bmod
      enddo
      !pause


      !--------------------------------------------------
      ! For i_save_disp=1 option:
      if(i_save_disp.eq.1) then
      ! scan x, find Nper^2(x,Npar) for cold plasma and hot plasma
      ! First, allocate arrays 
      !(could take large memory space, if xscan-grid is refined,
      ! and/or inper_save_disp is large)
      ! These are real*4 arrays, to save memory:
      allocate(Npara(inpar_save_disp),  STAT=istat)
      if(istat.ne.0) then
         print*,'Allocation problem for Npara()',istat
      endif
      Npara=0.d0
      allocate(Npera(inper_save_disp),  STAT=istat)
      if(istat.ne.0) then
         print*,'Allocation problem for Npera()',istat
      endif
      Npera=0.d0
      allocate(Nperp2m(ixscan_save_disp,inpar_save_disp),  STAT=istat)
      if(istat.ne.0) then
         print*,'Allocation problem for Nperp2m()',istat
      endif
      Nperp2m=0.d0
      allocate(Nperp2p(ixscan_save_disp,inpar_save_disp),  STAT=istat)
      if(istat.ne.0) then
         print*,'Allocation problem for Nperp2p()',istat
      endif
      Nperp2p=0.d0
      
      allocate(Nperp_im(ixscan_save_disp,inpar_save_disp),  STAT=istat)
      if(istat.ne.0) then
         print*,'Allocation problem for Nperp_im()',istat
      endif
      Nperp_im=0.d0

      allocate(Nperp2hot(ixscan_save_disp,inpar_save_disp),STAT=istat)
      if(istat.ne.0) then
         print*,'Allocation problem for Nperp2hot()',istat
      endif
      Nperp2hot=0.d0
      allocate(ddd(ixscan_save_disp,inper_save_disp,inpar_save_disp),
     +         STAT=istat)
      if(istat.ne.0) then
         print*,'Allocation problem for ddd()',istat
      endif
      ddd=0.d0
      !---

      id_loc=id ! save
      ibw_loc=ibw ! save
      id=6 ! used by hamilt_xyz ! HOT ROOTS
      print*,'wrtnetcdf_plasma_prof: Dispersion Function D(R,Nper,Npar)'
      dnpar= (Npar_mx_save_disp-Npar_mn_save_disp)/(inpar_save_disp-1)
      ! Linear scale:
      !dnper= (Nper_mx_save_disp-Nper_mn_save_disp)/(inper_save_disp-1)
      ! Parabolic scale (for better resolution at small Nperp):
      dnper=(Nper_mx_save_disp-Nper_mn_save_disp)/(inper_save_disp**2-1)
      if(dnper.le.0.d0)then
        WRITE(*,*)'wrtnetcdf_plasma_prof: check Npar_mn_save_disp'
        WRITE(*,*)'wrtnetcdf_plasma_prof: check Npar_mx_save_disp'
        WRITE(*,*)'set Npar_mx_save_disp>Npar_mn_save_disp'
        stop
      endif
      do kk=1,inpar_save_disp ! several values of Npar
         Npara(kk)= Npar_mn_save_disp +(kk-1)*dnpar ! min:max of Npar
         print*,'wrtnetcdf_plasma_prof: Dhot() for Npar=',Npara(kk)
         do j=1,inper_save_disp ! Nperp grid
            ! Linear scale:
            !Npera(j)= Nper_mn_save_disp +(j-1)*dnper ! min:max of Nper
            ! Parabolic scale (for better resolution at small Nperp):
            Npera(j)= Nper_mn_save_disp +(j**2 -1)*dnper 
            cnt2= Npara(kk)**2 + Npera(j)**2 ! N^2
            arg= Npara(kk)/sqrt(cnt2)  ! == N.b/(|N||b|) == Npar/N
            gam=dacos(arg) !-> to one.i ; used by hamilt_xyz
            print*,'wrtnetcdf_plasma_prof: Dhot() for Nperp=',Npera(j)
            do i=1,countnp(1) ! size of xscan
               x= xscan(i) ! scan from xmax to xmin (along x-axis)
               bmod=bxyz(x,y,z) !-> get b (needed for wcw)
c               write(*,'(a,2e12.3,a,e12.3,a,e12.3)') 'x,z[m]=', x,z,
c     +         '  (w/wpe)^2=',  1.0/wpw_2(x,y,z,1),
c     +         '   w/wce=',     1.0/wcw(x,y,z,1)
               ddd(i,j,kk)=hamilt_xyz(x,y,z,cnt2)
               ! Here, we simply save the value of dispersion function 
               ! over 2D grid (Nper,Npar). 
               ! Then, in genray_plot.py, 
               ! we plot the contour levels of ddd().
               ! Level=0 marks the roots of ddd()=0.
               ! (In fact, we only plot level=0.)
            enddo ! i=1,countnp(1) !  xscan
         enddo ! j=1,inper_save_disp
      enddo ! kk=1,inpar_save_disp

      cnper_0=1.d0 ! initialize
      do kk=1,inpar_save_disp ! several values of Npar
         do i=1,countnp(1) ! size of xscan
            x=xscan(i) ! scan along x-axis ! scan from xmax to xmin
            cnpar=dble(Npara(kk))
            id=2 ! used by npernpar_xyz  ! COLD ROOTS
            call npernpar_xyz(x,y,z,cnpar,cnper2p,cnper2m)
            Nperp2m(i,kk)=cnper2m ! X-mode !could be neg., if non-propagating
            Nperp2p(i,kk)=cnper2p ! O-mode !could be neg., if non-propagating
            !write(*,'(a,4e12.4)')'x,Npara(kk),Nperp2m,Nperp2p=',
     +      !  x,Npara,cnper2m,cnper2p

            cnper=0.d0 ! to initialize: cold root for Nperp
            if(ioxm.eq.-1 .and. cnper2m.gt.0.d0) cnper= sqrt(cnper2m)
            if(ioxm.eq. 1 .and. cnper2p.gt.0.d0) cnper= sqrt(cnper2p)
            bmod=bxyz(x,y,z) !-> get b
            dens_e= dense_xyz(x,y,z,1) !-> get rho and dens_e
            temp_e=tempe_xyz(x,y,z,1)

            !----------- YuP[04-2016] Add absorption -> Nper_im array
            
            !if ((iabsorp.eq.3).or.(iabsorp.eq.2)) then
	      !absorption for lh and fw
            !electric field using the cold plasma dielectric tensor
            call tensrcld_xyz(x,y,z)
            do is=2,nbulk
               tempiar(is)=tempe_xyz(x,y,z,is)
            enddo
            z_eff=zeffrho(rho) 

            if (iabsorp.eq.3 .and. cnper.gt.0.d0) then 
              !FW absorption, using cnper from id=2
              cnprim_cl=0.d0
            !absorpfd uses complex function modified bessel zfunc(argz,zf,zfp) 
            !call absorpf1(z,r,phi,cnpar,cnper,temp_e,dens_e,tempiar,
            !absorpfd uses double complex function modified bessel 
            !czeta(argz,zf,zfp,ierror) and calculates the dielectric tensor
            !reps )(in common eps.i)
            !for the electron plasma with the hot correction using
            !Chiu et al, Nucl.Fus Vol. 29, No.12(1989) p.2175
            !formula (2),(3),(4),(5),(6) and (7)
              rho_larm_max0=rho_larm_max ! YuP[11-2016]
              call absorpfd_xyz(x,y,z,f000,rho_larm_max0,cnpar,cnper,
     1             temp_e,dens_e,tempiar,nbulk,bmod,
     1             cnprim_e,cnprim_i,cnprim_s)
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
     +          'netcdfr3d: cnprim<0. cnprim_e,cnprim_i,cnpar,cnper=',
     +           cnprim_e,cnprim_i, cnpar,cnper,
     +          '  Resetting to abs(cnprim)'
                !if(cnprim_i.lt.0.d0) cnprim_i=0.d0 ! another version
                !if(cnprim_e.lt.0.d0) cnprim_e=0.d0 ! another version
                cnprim_e=dabs(cnprim_e)
                cnprim_i=dabs(cnprim_i)
                do is=1,nbulk
                   cnprim_s(is)=dabs(cnprim_s(is))
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
              Nperp_im(i,kk)=cnprim
            endif ! iabsorp=3

            if(iabsorp.eq.7) then
                 !EC wave case.The complex electric field calculations
                 !using Cold plasma tensor +antihermition relativistic tensor
                 !EC relativistic absorption 
                 !dielectric tensor=hermitian part(cold plasma)+
                 !anti-hermitian part(full relativistic) 
                 cnparp=cnpar
                 !if(cnper2p.gt.0.d0)then
                 !  cnperp=sqrt(cnper2p)
                 !elseif(cnper2m.gt.0.d0)then
                 !  cnperp=sqrt(cnper2m)
                 !else
                   cnperp=10. ! as a guess for hot plasma root.
                 !endif
                 !calculate Hermitian cold plasma complex tensor reps. 
                 !It will be in eps.i        
                 call tensrcld_xyz(x,y,z)
                 do is=1,nbulk
                    x_ar(is)=wpw_2(x,y,z,is)
                    y_ar(is)=wcw(x,y,z,is)
                    if(is.eq.1) y_ar(1)=-y_ar(1)
                    te=tempe_xyz(x,y,z,is) ! kev
                    t_av_ar(is)=te*1000.d0      ! ev 
                    tpop_ar(is)=tpoprho(rho,is)
                    vflow_ar(is)=vflowrho(rho,is)
                 enddo
                 if ((i_diskf.eq.0).or.(i_diskf.eq.5)) then
                  !usage of the analytical relativistic function and its derivatives 
                    i_fkin=0
                 else 
                  !usage of the mech relativistic function and its derivatives
                    i_fkin=1
                 endif
                 call anth_rlt_xyz(x_ar(1),y_ar(1), t_av_ar(1)*1.d-3,
     +                      cnparp,cnperp,
     +                      n_relt_harm1,n_relt_harm2,n_relt_intgr,
     +                      i_resonance_curve_integration_method,epsi,
     +                      i_fkin,x,y,z,
     +                      aK)
                 !complex dispersion function calculated from the sum of
                 !of the cold electron plasma dielectric tensor eps_h
                 !and the relativistic electron anti-hermition dielectric tensor eps_a
                 disp_func=dcold_rlt(reps,aK,cnparp,cnperp)
                 !calculate the derivative d(D_hermitian)/d(ReN_perp)
                 !from the electron cold plasma dispersion function D
                 ddnp=dDcold(reps,cnparp,cnperp)
                 cnprim = dabs(DIMAG(disp_func) / DREAL(ddnp))
                 !cnprim_cl=0.d0
                 !cnprim_e=cnprim
                 !cnprim_i=0.d0
                 Nperp_im(i,kk)=cnprim
            endif ! iabsorp.eq.7
c------------------------------------------------------------------
            if(iabsorp.eq.6) then
c       EC wave case.The complex electric field calculations
c       using Forest tensor (with antihermition part)
c       reps(3,3), and N_per=N_per(N_par) from forest solver
c-------EC relativistic absorption 
c       dielectric tensor=hermitian part(from Forest code)+
c                         anti-hermitian part(full relativistic) 
                 !-------Hermitian non-relativistic tensor reps        
                 cnparp=cnpar
                 !if(cnper2p.gt.0.d0)then
                 !  cnperp=sqrt(cnper2p)
                 !elseif(cnper2m.gt.0.d0)then
                 !  cnperp=sqrt(cnper2m)
                 !else
                   cnperp=10. ! as a guess for hot plasma root.
                 !endif
                 do is=1,nbulk
                    x_ar(is)=wpw_2(x,y,z,is)
                    y_ar(is)=wcw(x,y,z,is)
                    if(is.eq.1) y_ar(1)=-y_ar(1)
                    te=tempe_xyz(x,y,z,is) ! kev
                    t_av_ar(is)=te*1000.d0      ! ev 
                    tpop_ar(is)=tpoprho(rho,is)
                    vflow_ar(is)=vflowrho(rho,is)
                 enddo
                 if ((i_diskf.eq.0).or.(i_diskf.eq.5)) then
                  !usage of the analytical relativistic function and its derivatives 
                    i_fkin=0
                 else 
                  !usage of the mech relativistic function and its derivatives
                    i_fkin=1
                 endif
                 call anth_rlt_xyz(x_ar(1),y_ar(1), t_av_ar(1)*1.d-3,
     +                      cnparp,cnperp,
     +                      n_relt_harm1,n_relt_harm2,n_relt_intgr,
     +                      i_resonance_curve_integration_method,epsi,
     +                      i_fkin,x,y,z,
     +                      aK)
                 cnprim=0.d0
                 D=dhot_rlt(reps,aK,cnparp,cnperp,cnprim)
                 !dham=dreal(d)
                 call Ddhot(nbulk,dmas,x_ar,y_ar,
     +                t_av_ar,tpop_ar,vflow_ar,
     +                cnparp,cnperp,reps,dd5,ddnp_h,ddnll_h, ddnp)
                 cnprim = dabs(DIMAG(D) / DREAL(ddnp))
                 !cnprim_cl=0.d0
                 !cnprim_e=cnprim
                 !cnprim_i=0.d0
                 Nperp_im(i,kk)=cnprim
            endif ! iabsorp.eq.6

            !-> Try to find a Hot plasma root.
            Nperp2hot(i,kk)= cnper2m ! X-mode ! initialize
            goto 10 !---> skip this part - takes too long !comment if needed
            !initialization for common/nperpcom/,  in nperpcom.i
            nllc= Npara(kk)
            nbulkc=nbulk
            do j=1,nbulk
               massc(j)=dmas(j)
               xc(j)=wpw_2(x,y,z,j)
               bmod=bxyz(x,y,z) !-> get b (needed for wcw)
               yc(j)=wcw(x,y,z,j)
               if(j.eq.1) yc(1)=-yc(1) ! negative Y=(omega_ce/omega) for electrons
               tec(j)=tempe_xyz(x,y,z,j)*1.d+3 !(eV) averaged temperature
               tpopc(j)=tpoprho(rho,j)
               vflowc(j)=vflowrho(rho,j)         
            enddo
            iraystop=0 ! initialize
            if(cnper2m.gt.0.d0 .or. cnper2p.gt.0.d0)then
               ! Try X-mode or O-mode as an initial guess (the larger)
               ibw=0
               cnper_0= sqrt(max(cnper2m,cnper2p)) ! the larger root
               hotnpc=
     +         hotnp(nbulk,ibw,cnper_0,cnper2p,cnper2m,K,iraystop)
               cnper_0= real(hotnpc) ! save for the next r(i) step
            else ! No O or X mode
               !if(cnper_0.ne.1.d0) then ! use cnper_0 from previous r(i)
               !   ibw=1
               !   hotnpc=
     +         !   hotnp(nbulk,ibw,cnper_0,cnper_0,cnper_0,K,iraystop)
               !else
                  iraystop=1
               !endif
            endif
            if (iraystop.eq.0) then ! the root was found
               Nperp2hot(i,kk)=hotnpc*CONJG(hotnpc)
            else ! could not find the root
               Nperp2hot(i,kk)=1.d-6
            endif
   10       continue ! handle to skip hot root calculations

         enddo ! i=1,countnp(1) !  xscan
      enddo ! kk=1,inpar_save_disp
      id=id_loc   ! restore
      ibw=ibw_loc ! restore
      endif ! i_save_disp=1
c-------------------------------------------------------

      ! Write data
      write(*,*)'wrtnetcdf_plasma_prof: Start writing data'
      
      call ncvid2(vid,ncid,'nbulk',istatus)
      call ncvpt_int2(ncid,vid,0,0,nbulk,istatus)

      call ncvid2(vid,ncid,'model_b',istatus)
      call ncvpt_int2(ncid,vid,0,0,model_b,istatus)

      call ncvid2(vid,ncid,'model_rho_dens',istatus)
      call ncvpt_int2(ncid,vid,0,0,model_rho_dens,istatus)

      call ncvid2(vid,ncid,'freqcy',istatus)
      call ncvpt_doubl2(ncid,vid,0,0,freqcy,istatus)

      call ncvid2(vid,ncid,'rs_frc',istatus) ! for model_b=4
      call ncvpt_doubl2(ncid,vid,0,0,rs_frc,istatus)

      call ncvid2(vid,ncid,'y_save_disp',istatus)
      call ncvpt_doubl2(ncid,vid,0,0,y_save_disp,istatus)

      call ncvid2(vid,ncid,'z_save_disp',istatus)
      call ncvpt_doubl2(ncid,vid,0,0,z_save_disp,istatus)

      call ncvid2(vid,ncid,'dmas',istatus)
      call ncvpt_doubl2(ncid,vid,1,nbulk,dmas,istatus)

      call ncvid2(vid,ncid,'charge',istatus)
      call ncvpt_doubl2(ncid,vid,1,nbulk,charge,istatus)

      call ncvid2(vid,ncid,'rho_bin',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR,rho_bin,istatus)
      call check_err(istatus)
      
      call ncvid2(vid,ncid,'rho_bin_center',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR-1,rho_bin_center,istatus)
      call check_err(istatus)
            
      call ncvid2(vid,ncid,'zefprof',istatus)
      call ncvpt_doubl2(ncid,vid,start(1),counts(1),zefprof,istatus)
      call check_err(istatus)
      
      call ncvid2(vid,ncid,'temprof',istatus)
      call ncvpt_doubl2(ncid,vid,start,counts,temprof,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'densprof',istatus)
      call ncvpt_doubl2(ncid,vid,start,counts,densprof,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'bmodprofxz',istatus) ! big increase in file size
      call ncvpt_doubl2(ncid,vid,start,countxz,bmodprofxz,istatus)
      call check_err(istatus)
      
      call ncvid2(vid,ncid,'densprofxz',istatus) ! big increase in file size
      call ncvpt_doubl2(ncid,vid,start,countxz,densprofxz,istatus)
      call check_err(istatus)
      
      call ncvid2(vid,ncid,'densprofyz',istatus) ! big increase in file size
      call ncvpt_doubl2(ncid,vid,start,countyz,densprofyz,istatus)
      call check_err(istatus)
      
      call ncvid2(vid,ncid,'densprofxy',istatus) ! big increase in file size
      call ncvpt_doubl2(ncid,vid,start,countxy,densprofxy,istatus)
      call check_err(istatus)
      
      deallocate(bmodprofxz, densprofxz, densprofyz, densprofxy)  
      
      call ncvid2(vid,ncid,'w_x_densprof_nc',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR,w_x_densprof_nc,istatus)
      call check_err(istatus)
 
      call ncvid2(vid,ncid,'w_bmod_vs_x_nc',istatus)
      call ncvpt_doubl2(ncid,vid,1,NR,w_bmod_vs_x_nc,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'w_dens_vs_x_nc',istatus)
      call pack21(w_dens_vs_x_nc,1,NR,1,nbulk,tem1,NR,nbulk)
      call ncvpt_doubl2(ncid,vid,start,counts,tem1,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'w_temp_vs_x_nc',istatus)
      call pack21(w_temp_vs_x_nc,1,NR,1,nbulk,tem1,NR,nbulk)
      call ncvpt_doubl2(ncid,vid,start,counts,tem1,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'w_zeff_vs_x_nc',istatus)
      call ncvpt_doubl2(ncid,vid,1,NR,w_zeff_vs_x_nc,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'w_y_densprof_nc',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR,w_y_densprof_nc,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'w_bmod_vs_y_nc',istatus)
      call ncvpt_doubl2(ncid,vid,1,NR,w_bmod_vs_y_nc,istatus)
      call check_err(istatus)
 
      call ncvid2(vid,ncid,'w_dens_vs_y_nc',istatus)
      call pack21(w_dens_vs_y_nc,1,NR,1,nbulk,tem1,NR,nbulk)
      call ncvpt_doubl2(ncid,vid,start,counts,tem1,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'w_temp_vs_y_nc',istatus)
      call pack21(w_temp_vs_y_nc,1,NR,1,nbulk,tem1,NR,nbulk)
      call ncvpt_doubl2(ncid,vid,start,counts,tem1,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'w_zeff_vs_y_nc',istatus)
      call ncvpt_doubl2(ncid,vid,1,NR,w_zeff_vs_y_nc,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'xscan',istatus) ! x-grid for Nperp^2(x,Npar)
      call ncvpt_real(ncid,vid,1,countnp(1),xscan,istatus)
      call check_err(istatus)
      
      call ncvid2(vid,ncid,'wcwscan',istatus) ! omega_ce/omega
      call ncvpt_real(ncid,vid,1,countnp(1),wcwscan,istatus)
      call check_err(istatus)
      
      call ncvid2(vid,ncid,'wpwscan',istatus) ! omega_pe/omega
      call ncvpt_real(ncid,vid,1,countnp(1),wpwscan,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'wuwscan',istatus) ! Upper_hybrid/omega
      call ncvpt_real(ncid,vid,1,countnp(1),wuwscan,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'rhoscan',istatus) ! rho vs x
      call ncvpt_real(ncid,vid,1,countnp(1),rhoscan,istatus)
      call check_err(istatus)
      
      ! For i_save_disp=1 option:
      if(i_save_disp.eq.1)then !----------------------------------------
      call ncvid2(vid,ncid,'Npara',istatus) !Npar-grid for Nperp^2(x,Npar)
      call ncvpt_real(ncid,vid,1,countnp(2),Npara,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'Npera',istatus) !Nper-grid for ddd(x,Nper,Npar)
      call ncvpt_real(ncid,vid,1,countddd(2),Npera,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'Nperp2m',istatus) ! Nperp^2(x,Npar) Cold'-1'
      call ncvpt_real(ncid,vid,start,countnp,Nperp2m,istatus)
      call check_err(istatus)
      
      call ncvid2(vid,ncid,'Nperp2p',istatus) ! Nperp^2(x,Npar) Cold'+1'
      call ncvpt_real(ncid,vid,start,countnp,Nperp2p,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'Nperp_im',istatus) !Nperp_im(x,Npar) damping
      call ncvpt_real(ncid,vid,start,countnp,Nperp_im,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'Nperp2hot',istatus) !Nperp^2(x,Npar) Hot root
      call ncvpt_real(ncid,vid,start,countnp,Nperp2hot,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'ddd',istatus) !ddd(x,Nper,Npar) dispers.function
      call ncvpt_real(ncid,vid,startddd,countddd,ddd,istatus)
      call check_err(istatus)
      endif ! i_save_disp=1 !-------------------------------------------
      
      write(*,*)'wrtnetcdf_plasma_prof: Finished writing data'

c     Close netCDF file
      call ncclos3(ncid,istatus)
      call check_err(istatus)
      write(*,*)'wrtnetcdf_plasma_prof: File is Closed'

      return
      end


C=======================================================================
C=======================================================================
     
     
      subroutine wrtnetcdf(kopt)
c      implicit integer (i-n), double precision (a-h,o-z)
      implicit none
c
c     Write ray tracing data (as in mnemonic.txt) into a netCDF file.

cSm0940727
c     If the parameter ionetwo.eq.1 it will write the current
c     and power profiles into a netcdf file

      include 'param.i'
      include 'writencdf.i'
      include 'one.i'
      include 'ions.i'
      include 'fourb.i'
      include 'five.i'
      include 'cone_nml.i'     !nccone
      include 'grill_nml.i'    !ngrill
c--------------------------
cSm040727 to write the current and power profiles into netcdf file
c     Done in subroutine, wrtnetcdf_prof
c      include 'onetwo.i'
c--------------------------
c-----input
      integer
     & kopt  !1 define dimensions,variables
c            !and attributes
             !0 Write data into created netcdf file
     
c-----local     
      integer n1n2,istat
c     Storage tem1 is used in netcdf writes, including complex numbers.
cSAP090905
c      parameter (n1n2=2*nrelta*nraymax)
       real*8, pointer :: tem1(:)
       real*4, pointer :: tem1_single(:)
    
c --- include file for netCDF declarations 
c --- (obtained from NetCDF distribution)
      include 'netcdf.inc'

c --- some stuff for netCDF file ---
      integer ncid,vid,istatus,nccre2,ncvdef2,ncdid2,ncddef2
      integer nraysid,neltmaxid,twoid,char512id,
cSAP
     &char8id
    
      integer nbulkmid  !ion species dimension
      integer nbulkid
      integer nrho_adjid

      integer ray_dims(3),start(3),ray_count(3)
      integer ray_dimss(3),starts(3),ray_counts(3)
      integer ray_dimsp(3)
      ! Note: ray_dimsp includes ALL species, while ray_dimss - ions only

      double complex cei

      character ltitle*512

      integer neltmax,iray,i,j,ii,ll
      integer length_char

      data start/1,1,1/,starts/1,1,1/
       
      save

cSAP090903
      n1n2=2*nrelta*nrayl
c      write(*,*)'nrelta,nrayl,n1n2',nrelta,nrayl,n1n2
c------------------------------------------
c     allocate pointers tem1
c-------------------------------------------
      allocate( tem1(1:n1n2),STAT=istat)
      allocate( tem1_single(1:n1n2),STAT=istat)
      tem1=0.d0
      tem1_single=0.0 
c      write(*,*)'wrtnetcdf after allocate tem1 istat=',istat
c------------------------------------------
      cei=(0.d0,1.d0)

c     Maximum number of ray elements per ray:
      neltmax=0

c      write(*,*)'wrtnetcdf nray=',nrayl


      do iray=1,nrayl
         neltmax=max(neltmax,nrayelt_nc(iray))
c         write(*,*)'wrtnetcdf iray,nrayelt_nc(iray),neltmax',
c     &                        iray,nrayelt_nc(iray),neltmax
      enddo

cSm05038
      if (neltmax.eq.0) neltmax=1

      ray_count(1)=neltmax
      ray_count(2)=nrayl

      ray_count(3)=2
      ray_counts(1)=neltmax
      ray_counts(2)=nrayl
      ray_counts(3)=1

c      write(*,*)'ray_count',ray_count

c.......................................................................
cl    1. Initialize part
c

c --- begin if ---
      if ( kopt.eq.1 ) then
      write(*,*)'wrtnetcdf Open existing file, define dims,  kopt=1'

c.......................................................................
c     1.1.1 open existing netCDF file: netcdfnml,
c     create netCDF filename (Entering define mode.)
c     integer function ncopn(filename,write,error_code)
c     Ref to page 40 of NetCDF-2 manual.
c     CLOBber old file, if it exists.
c     istatus is 0, if no errors.
      
      call ncopn2(filenc,NCWRITE,ncid,istatus) ! Open existing netCDF file
      call check_err(istatus)
      
c.......................................................................
c     1.1.2 
c     Put Open NetCDF file into Define Mode
c     subroutine ncredf3(integer ncid,integer error_code)
c     p. 41 of netcdf manual
c
c     call ncredf3(integer ncid,integer error_code)
     
      call ncredf3(ncid,istatus)
      call check_err(istatus)
        
c.......................................................................
cl    1.1.2 define dimensions
c     integer function ncddef(ncid,dim_name,dim_siz,error_code)
c       returns dimension id.
c     p. 67 of netcdf-2 manual

c     For ray data:
c      write(*,*)'wrtnetcdf before ncddef(neltmax) neltmax=',neltmax

      neltmaxid=ncddef2(ncid,'neltmax',neltmax,istatus)         
     
c      write(*,*)'wrtnetcdf before ncddef(ncid,nrays) nrayl=',nrayl

      nraysid=ncddef2(ncid,'nrays',nrayl,istatus)

c      write(*,*)'wrtnetcdf before ncddef(ncid,two,2,istatus)'
      
      twoid=ncddef2(ncid,'two',2,istatus)
    
c      write(*,*)'wrtnetcdf before ncddef(ncid,nbulk)'

      !nbulkid=ncddef2(ncid,'nbulk',nbulk,istatus) 
      ! Get dimension ID from dimension name
      nbulkid=ncdid2(ncid,'nbulk',istatus)

cSAP080303
c      write(*,*)'wrtnetcdf before ncddef(char8dim)'
      char8id=ncddef2(ncid,'char8dim',8,istatus)

c      write(*,*)'wrtnetcdf before ncddef(char64dim)'


c      write(*,*)'wrtnetcdf before ncddef(char128dim)'


c      write(*,*)'wrtnetcdf before ncddef(char5212dim)'

      char512id=ncddef2(ncid,'char512dim',512,istatus)

c      write(*,*)'neltmaxid',neltmaxid
c      write(*,*)'nraysid',nraysid
c      write(*,*)'twoid',twoid

      ray_dims(1)=neltmaxid
      ray_dims(2)=nraysid
      ray_dims(3)=twoid
    
      if (nbulk.gt.1) then
         nbulkmid=ncddef2(ncid,'nbulkm',nbulk-1,istatus)
         ray_dimss(1)=ray_dims(1)
         ray_dimss(2)=ray_dims(2)
         ray_dimss(3)=nbulkmid
      endif
      ! Note: ray_dimsp includes ALL species, while ray_dimss - ions only
      ray_dimsp(1)=ray_dims(1)
      ray_dimsp(2)=ray_dims(2)
      ray_dimsp(3)=nbulkid
      
c.......................................................................
cl    1.1.3 define variables

c     Note, the variable IDs (denoted below as vid) are
c     not saved here in this subroutine; rather, the IDs
c     are retrieved from the netCDF data file, as needed,
c     by calling the netCDF routine ncvid.

c     netCDF variable_define:
c     integer function ncvdef(ncid,variable_name,variable_type,
c                number_of_dimensions,
c                vector_for_length_per_dimension,error_code)
c       returns varid.
c     Note: Unlimited dimension must be last
c     Refer to p. 77, netcdf-2 manual

c     NCDOUBLE for REAL*8.  This is independent of the
c       internal representation of data on a particular
c       machine, e.g., 32- or 64-bit based calculations.



c-------------------------------
c     Run Descriptive Parameters
c-------------------------------
c      write(*,*)'before ncvdef2(vid)'
      vid=ncvdef2(ncid,'id',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,24,
     +            'Disp relation identifier',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(iabsorp)'
      vid=ncvdef2(ncid,'iabsorp',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,22,
     +            'Absorp calc identifier',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(ieffic)'
      vid=ncvdef2(ncid,'ieffic',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,29,
     +            'Current drive calc identifier',istatus)
      call check_err(istatus)

cSAP080303
c      write(*,*)'before ncvdef2(ion_absorption)'
      vid=ncvdef2(ncid,'ion_absorption',NCCHAR,1,char8id,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,49,
     +'Switch on/off ion absorption at iabsorp=3,9,91,92',istatus)
      call check_err(istatus)


c      write(*,*)'before ncvdef2(refl_loss)'
      vid=ncvdef2(ncid,'refl_loss',NCDOUBLE,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,41,
     +    'fraction of power loss at each reflection',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(iflux)'
      vid=ncvdef2(ncid,'iflux',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,34,
     +            'Flux calc, non-Westerhof-Tokman id',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(ioxm)'
      vid=ncvdef2(ncid,'ioxm',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,38,
     +            'wave mode indicator (1 - om, -1 - xm )',istatus)
      call check_err(istatus)


c      write(*,*)'before ncvdef2(ioxm_n_npar)'
      vid=ncvdef2(ncid,'ioxm_n_npar',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,68,
     +'wave mode indicator: sign before square root to find '//
     +' N(N_parallel)',
     +istatus)
      call check_err(istatus)


c      write(*,*)'before ncvdef2(jwave)'
      vid=ncvdef2(ncid,'jwave',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,37,
     +            'Wave harmonic, for CD efficiency calc',istatus)
      call check_err(istatus)

      if (iabsorp.eq.4) then
      vid=ncvdef2(ncid,'i_im_nperp',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,57,
     +       'iabsorp=4  nperp:1,ImD_full/dD/dnp; 2,Cmplx soln',istatus)
      call check_err(istatus)
      endif

c      write(*,*)'before ncvdef2(i_geom_optic)'
      vid=ncvdef2(ncid,'i_geom_optic',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,34,
     +            'Integrate rays wrt (1)time,(2)dist',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(istart)'
      vid=ncvdef2(ncid,'istart',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,50,
     +     'Ray launch type: 1,eccone; 2,grill, 3,OX in plasma',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(ncone)'
      vid=ncvdef2(ncid,'ncone',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name0',NCCHAR,39,
     +     'Number of rf source ray cones, istart=1',istatus)
      call ncaptc2(ncid,vid,'long_name1',NCCHAR,40,
     +     'Each cone has (nray/ncone) launched rays, istart=1',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(ngrill)'
      vid=ncvdef2(ncid,'ngrill',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name0',NCCHAR,41,
     +     'Number of rf source ray grills, istart=2,3',istatus)
      call ncaptc2(ncid,vid,'long_name0',NCCHAR,47,
     +     ',nray/ngrill rays launched per grill istart=2,3',istatus)
      call ncaptc2(ncid,vid,'long_name0',NCCHAR,55,
     +     ',Each grill has (nray/ngrill) launched rays, istart=2,3',
     +     istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(ionetwo)'
      vid=ncvdef2(ncid,'ionetwo',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +            'if ionetwo=1 then calculate CD',istatus)
      call check_err(istatus)
c--------------------------
c     Ray data
c--------------------------
c      write(*,*)'before ncvdef2(nray)'
      vid=ncvdef2(ncid,'nray',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,14,
     +            'Number of rays',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(nharm)'
      vid=ncvdef2(ncid,'nharm',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,21,
     +            'First harmonic number',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(nrayelt)'

      vid=ncvdef2(ncid,'nrayelt',NCLONG,1,nraysid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,35,
     +            'Number of ray elements for each ray',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(ws)'

      vid=ncvdef2(ncid,'ws',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,20,
     +           'distance along a ray',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'cms',istatus)

c      write(*,*)'before ncvdef2(seikon)'

      vid=ncvdef2(ncid,'seikon',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,7,
     +           'eikonal',istatus)


c      write(*,*)'before ncvdef2(spsi)'

      vid=ncvdef2(ncid,'spsi',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,45,
     +         'normalized small radius=rho given by indexrho',istatus)

c      write(*,*)'before ncvdef2(wr)'

      vid=ncvdef2(ncid,'wr',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,12,
     +           'major radius',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'cms',istatus)

c      write(*,*)'before ncvdef2(wphi)'

      vid=ncvdef2(ncid,'wphi',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,14,
     +           'toroidal angle',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'radians',istatus)

c      write(*,*)'before ncvdef2(wz)'
      vid=ncvdef2(ncid,'wx',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,12,
     +           'x-cartesian',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'cms',istatus)

      vid=ncvdef2(ncid,'wy',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,12,
     +           'y-cartesian',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'cms',istatus)

      vid=ncvdef2(ncid,'wz',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,15,
     +           'vertical height',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'cms',istatus)

c      write(*,*)'before ncvdef2(w_theta_pol)'

      vid=ncvdef2(ncid,'w_theta_pol',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,14,
     +           'poloidal angle',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,6,
     +           'degree',istatus)

c      write(*,*)'before ncvdef2(wnpar)'

      vid=ncvdef2(ncid,'wnpar',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +           'parallel refractive index',istatus)

c      write(*,*)'before ncvdef2(wnper)'

      vid=ncvdef2(ncid,'wnper',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +           'perpendicular refractive index',istatus)

c      write(*,*)'before ncvdef2(delpwr)'

      vid=ncvdef2(ncid,'delpwr',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,20,
     +           'power in ray channel',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +           'erg/sec',istatus)

      vid=ncvdef2(ncid,'sdpwr',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,48,
     +      'Ion collisionless absorption coeff (all species)',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           '1/cms',istatus)

      if (iabsorp.eq.3 .and. nbulk.gt.1) then !salphas 
      vid=ncvdef2(ncid,'salphas',NCFLOAT,3,ray_dimss,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,49,
     +     'Ion collisionless absorption coeff (each species)',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           '1/cms',istatus)
      endif

      vid=ncvdef2(ncid,'vthermal',NCFLOAT,3,ray_dimsp,istatus)
      call check_err(istatus)
      ! Note: ray_dimsp includes ALL species, while ray_dimss - ions only
      call ncaptc2(ncid,vid,'long_name',NCCHAR,39,
     +     'Thermal speed along ray (each species)',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,4,
     +           'cm/s',istatus)

      vid=ncvdef2(ncid,'wdnpar',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,6,
     +           'wdnpar',istatus)

c     Added 3rd dimension equal to 2 accomodates complex data.
c      write(*,*)'ncid=',ncid
c      write(*,*)'ray_dims',ray_dims
      vid=ncvdef2(ncid,'cwexde',NCDOUBLE,3,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +           'Complex Ex/E Polarization',istatus)

      vid=ncvdef2(ncid,'cweyde',NCDOUBLE,3,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +           'Complex Ey/E Polarization',istatus)

      vid=ncvdef2(ncid,'cwezde',NCDOUBLE,3,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +           'Complex Ez/E Polarization',istatus)

      vid=ncvdef2(ncid,'fluxn',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23,
     +           'fluxn, Stix norm, |E|=1',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,13,
     +           'ergs/sec/cm^2',istatus)

      vid=ncvdef2(ncid,'sbtot',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23,
     +           'Magnetic field strength',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           'Gauss',istatus)

c      write(*,*)'before ncvdef2(cene)'

      vid=ncvdef2(ncid,'sene',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,17,
     +           'Density along ray',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,14,
     +           'particles/cm^3',istatus)
     
c      write(*,*)'before ncvdef2(ste)'

      vid=ncvdef2(ncid,'ste',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,21,
     +           'Temperature along ray',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'keV',istatus)

c      write(*,*)'before ncvdef2(salphac)'

      vid=ncvdef2(ncid,'salphac',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +           'Collisional damping wavenumber',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           '1/cms',istatus)

c      write(*,*)'before ncvdef2(salphal)'

      vid=ncvdef2(ncid,'salphal',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +           'Linear damping wavenumber',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           '1/cms',istatus)

c      write(*,*)'before ncvdef2(sb_r)'

      vid=ncvdef2(ncid,'sb_r',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,28,
     +           'B_r magnetic field component',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           'Gauss',istatus)

c      write(*,*)'before ncvdef2(sb_z)'

      vid=ncvdef2(ncid,'sb_x',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,28,
     +           'B_x magnetic field component',istatus) 
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           'Gauss',istatus)

      vid=ncvdef2(ncid,'sb_y',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,28,
     +           'B_y magnetic field component',istatus) 
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           'Gauss',istatus)

      vid=ncvdef2(ncid,'sb_z',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,28,
     +           'B_z magnetic field component',istatus) 
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           'Gauss',istatus)

c      write(*,*)'before ncvdef2(sb_phi)'

      vid=ncvdef2(ncid,'sb_phi',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +           'B_phi magnetic field component',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           'Gauss',istatus)

c      write(*,*)'before ncvdef2(wn_r)'

      vid=ncvdef2(ncid,'wn_r',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +           'N_r refractive index component',istatus)

c      write(*,*)'before ncvdef2(wn_z)'

      vid=ncvdef2(ncid,'wn_x',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +           'N_x refractive index component',istatus) 

      vid=ncvdef2(ncid,'wn_y',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +           'N_y refractive index component',istatus) 

      vid=ncvdef2(ncid,'wn_z',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +           'N_z refractive index component',istatus) 

c      write(*,*)'before ncvdef2(wn_phi)'

      vid=ncvdef2(ncid,'wn_phi',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,32,
     +           'N_phi refractive index component',istatus)
     
c      write(*,*)'before ncvdef2(wgr_r)'

      vid=ncvdef2(ncid,'vgr_r',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,24,
     +           'vgroup_r normalized to c',istatus)
     
c      write(*,*)'before ncvdef2(wgr_z)'

      vid=ncvdef2(ncid,'vgr_x',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,24,
     +           'vgroup_x normalized to c',istatus)
     
      vid=ncvdef2(ncid,'vgr_y',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,24,
     +           'vgroup_y normalized to c',istatus)
     
      vid=ncvdef2(ncid,'vgr_z',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,24,
     +           'vgroup_z normalized to c',istatus)
     
      vid=ncvdef2(ncid,'vgr_phi',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,26,
     +           'vgroup_phi normalized to c',istatus)
     
      vid=ncvdef2(ncid,'flux_z',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,22,
     +           'flux_z normalized to c',istatus)
       
      vid=ncvdef2(ncid,'flux_r',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,22,
     +           'flux_r normalized to c',istatus)
            
      vid=ncvdef2(ncid,'flux_phi',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,24,
     +           'flux_phi normalized to c',istatus)

c--------------------------------------------------------
c      CD efficiency along rays
c--------------------------------------------------------
      if (ionetwo.eq.1)then
c         write(*,*)'before ncvdef2(w_eff_nc)'
         vid=ncvdef2(ncid,'w_eff_nc',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +           'CD efficiency along a ray',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,27,
     +           '(A/cm**2)/(erg/(sec*cm**3))',istatus)
       endif
c--------------------------------------------------------
c      DC electric field conductivity from adj 
c--------------------------------------------------------
c--------------------------
c     determine OX conversion data for i_ox=2 case
c--------------------------
      if (i_ox.eq.2) then

         vid=ncvdef2(ncid,'i_ox_conversion',NCLONG,1,nraysid,istatus)
         call check_err(istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,35,
     &            'Option equals 1 after OX conversion',istatus)
         call check_err(istatus)

         vid=ncvdef2(ncid,'transm_ox',NCDOUBLE,1,nraysid,istatus)
         call check_err(istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,57,
     &   'OX transmission coefficient Preinhaelter and Kopecky 1973',
     &              istatus)
         call check_err(istatus)

         vid=ncvdef2(ncid,'cn_par_optimal',NCDOUBLE,1,nraysid,istatus)
         call check_err(istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,38,
     &              'optimal N parallel for OX transmission',
     &              istatus)
         call check_err(istatus)

         vid=ncvdef2(ncid,'cnpar_ox',NCDOUBLE,1,nraysid,istatus)
         call check_err(istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,33,
     &              'N parallel before OX transmission',
     &              istatus)
         call check_err(istatus)

         vid=ncvdef2(ncid,'cn_b_gradpsi',NCDOUBLE,1,nraysid,istatus)
         call check_err(istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,45,
     &              'N along [B^*grad(psi)] before OX transmission',
     &              istatus)
         call check_err(istatus)
         
      endif !i_ox=2

c--------------------------
c     Plasma data
c--------------------------
            
      ! Get dimension ID from dimension name
      nbulkid=ncdid2(ncid,'nbulk',istatus)



c------------------------------------------------------------------
c     for total power [erg/sec] absorbed at all reflections at all rays
c------------------------------------------------------------------
c      write(*,*)'before ncvdef2(w_tot_pow_absorb_at_refl_nc)'
      vid=ncvdef2(ncid,'w_tot_pow_absorb_at_refl_nc',NCDOUBLE,
     &0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,47,
     +    'Total power absorbed at reflections of all rays',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +           'erg/sec',istatus)
c.................................
cl    1.1.4 end the define-mode and start the data-mode
c     p. 51-2 of manual
c      write(*,*)'ncid',ncid
      call ncendf3(ncid,istatus)
      call check_err(istatus)

c     Close netCDF file
      call ncclos3(ncid,istatus) 
      call check_err(istatus)
      ! it will be re-opened later for writing data
      
      endif               ! End initialize




c.......................................................................
cl    1. Writing data
c

      if ( kopt.eq.0 ) then

      call ncopn2(filenc,NCWRITE,ncid,istatus) !Open existing netCDF file
      call check_err(istatus)

c      write(*,*)'wrtnetcdf,  kopt=',kopt
cl    1.2 Write data
c

c     First get variable_id:
c             nf_inq_varid(netcdf_id,variable_name,integer_info)
c     Then write data with nc variable_put:
c             nf_put_var... ()
c
c     function nf_put_var...(ncid,variable_id,index,val)
c     p. 52,, 54-65 of netcdf-3 manual
c
c--------------------------

      call ncvid2(vid,ncid,'version',istatus)
      ll=length_char(version)
      call ncvptc2(ncid,vid,1,ll,version,ll,istatus)
c      write(*,*)'wrtnetcdf after ncptc version istatus',istatus

      call ncvid2(vid,ncid,'mnemonic',istatus)
      ll=length_char(mnemonic)
      call ncvptc2(ncid,vid,1,ll,mnemonic,ll,istatus)
c      write(*,*)'wrtnetcdf after ncptc mnemonic istatus',istatus

      call ncvid2(vid,ncid,'eqdskin',istatus)
      ll=length_char(eqdskin)
      call ncvptc2(ncid,vid,1,ll,eqdskin,ll,istatus)
c      write(*,*)'wrtnetcdf after ncptc eqdskin istatus',istatus

      call ncvid2(vid,ncid,'istart',istatus)
      call ncvpt_int2(ncid,vid,0,0,istart,istatus)
c      write(*,*)'wrtnetcdf after ncptc istart istatus',istatus
  
      call ncvid2(vid,ncid,'ncone',istatus)
      call ncvpt_int2(ncid,vid,0,0,ncone,istatus)
c      write(*,*)'wrtnetcdf after ncptc nccone istatus',istatus

      call ncvid2(vid,ncid,'ngrill',istatus)
      call ncvpt_int2(ncid,vid,0,0,ngrill,istatus)
c      write(*,*)'wrtnetcdf after ncptc ngrill istatus',istatus

      call ncvid2(vid,ncid,'ionetwo',istatus)
      call ncvpt_int2(ncid,vid,0,0,ionetwo,istatus)
c      write(*,*)'wrtnetcdf after ncptc ionetwo istatus',istatus

c--------------------------
c     Run specs
c--------------------------
      call ncvid2(vid,ncid,'id',istatus)
      call ncvpt_int2(ncid,vid,0,0,id,istatus)

c      write(*,*)'netcdfr3d.f refl_loss',refl_loss

      call ncvid2(vid,ncid,'refl_loss',istatus)
      call ncvpt_doubl2(ncid,vid,0,0,refl_loss,istatus)
c      write(*,*)'wrtnetcdf after ncvpt refl_loss istatus',istatus

      call ncvid2(vid,ncid,'iabsorp',istatus)
      call ncvpt_int2(ncid,vid,0,0,iabsorp,istatus)
c      write(*,*)'wrtnetcdf after ncptc iabsorp istatus',istatus


      call ncvid2(vid,ncid,'ieffic',istatus)
      call ncvpt_int2(ncid,vid,0,0,ieffic,istatus)
c      write(*,*)'wrtnetcdf after ncptc ieffic istatus',istatus


cSAP080303
      call ncvid2(vid,ncid,'ion_absorption',istatus)
      ll=length_char(ion_absorption)
      call ncvptc2(ncid,vid,1,ll,ion_absorption,ll,istatus)
c      write(*,*)'wrtnetcdf after ncptc ion_absorption istatus',istatus

      call ncvid2(vid,ncid,'iflux',istatus)
      call ncvpt_int2(ncid,vid,0,0,iflux,istatus)
c      write(*,*)'wrtnetcdf after ncpvt iflux istatus',istatus

      call ncvid2(vid,ncid,'ioxm',istatus)
      call ncvpt_int2(ncid,vid,0,0,ioxm,istatus)
c      write(*,*)'wrtnetcdf after ncvpt ioxm istatus',istatus

      call ncvid2(vid,ncid,'ioxm_n_npar',istatus)
      call ncvpt_int2(ncid,vid,0,0,ioxm_n_npar,istatus)
c      write(*,*)'wrtnetcdf after ncvpt ioxm_n_npar istatus',istatus


      call ncvid2(vid,ncid,'jwave',istatus)
      call ncvpt_int2(ncid,vid,0,0,jwave,istatus)
c      write(*,*)'wrtnetcdf after ncvpt jwave istatus',istatus

      call ncvid2(vid,ncid,'i_geom_optic',istatus)
      call ncvpt_int2(ncid,vid,0,0,i_geom_optic,istatus)
c      write(*,*)'wrtnetcdf after ncvpt i_geom_optic istatus',istatus

      
      if (iabsorp.eq.4) then
         call ncvid2(vid,ncid,'i_im_nperp',istatus)
         call ncvpt_int2(ncid,vid,0,0,i_im_nperp,istatus)
      endif


c--------------------------
c     Ray data
c--------------------------
      call ncvid2(vid,ncid,'nray',istatus)
      call ncvpt_int2(ncid,vid,0,0,nrayl,istatus)

      call ncvid2(vid,ncid,'nharm',istatus)
      call ncvpt_int2(ncid,vid,0,0,nharm,istatus)

      call ncvid2(vid,ncid,'nrayelt',istatus)
      call ncvpt_int2(ncid,vid,1,nrayl,nrayelt_nc,istatus)

      call pack21(ws_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'ws',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(seikon_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'seikon',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(spsi_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'spsi',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(wr_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'wr',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(wphi_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'wphi',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(wx_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'wx',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(wy_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'wy',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(wz_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'wz',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(w_theta_pol_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'w_theta_pol',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(wnpar_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'wnpar',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(wnper_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'wnper',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(delpwr_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'delpwr',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(sdpwr_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'sdpwr',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      if (iabsorp.eq.3) then !salphas_nc
      do i=2,nbulk
      call pack21_single(salphas_nc(1,1,i),1,nrelta,1,nrayl,
     +                   tem1_single,neltmax,nrayl)
      starts(1)=start(1)
      starts(2)=start(2)
      starts(3)=i-1  ! nbulk-1 ion species salphas are put in .nc file
      call ncvid2(vid,ncid,'salphas',istatus)
      call ncvpt_real(ncid,vid,starts,ray_counts,tem1_single,istatus)
      enddo
      endif

      do i=1,nbulk
      call pack21_single(wvthermal_nc(1,1,i),1,nrelta,1,nrayl,
     +                   tem1_single,neltmax,nrayl)
      starts(1)=start(1)
      starts(2)=start(2)
      starts(3)=i  ! nbulk species are put in .nc file
      call ncvid2(vid,ncid,'vthermal',istatus)
      call ncvpt_real(ncid,vid,starts,ray_counts,tem1_single,istatus)
      enddo

      call pack21(wdnpar_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'wdnpar',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      ii=0
      do j=1,nrayl
         do i=1,neltmax
            ii=ii+1
            tem1(ii)=0.5d0*(cwexde_nc(i,j)+dconjg(cwexde_nc(i,j)))
         enddo
      enddo
      do j=1,nrayl
         do i=1,neltmax
            ii=ii+1
            tem1(ii)=-cei*0.5d0*(cwexde_nc(i,j)-dconjg(cwexde_nc(i,j)))
         enddo
      enddo
      call ncvid2(vid,ncid,'cwexde',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      ii=0
      do j=1,nrayl
         do i=1,neltmax
            ii=ii+1
            tem1(ii)=0.5d0*(cweyde_nc(i,j)+dconjg(cweyde_nc(i,j)))
         enddo
      enddo
      do j=1,nrayl
         do i=1,neltmax
            ii=ii+1
            tem1(ii)=-cei*0.5d0*(cweyde_nc(i,j)-dconjg(cweyde_nc(i,j)))
         enddo
      enddo
      call ncvid2(vid,ncid,'cweyde',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      ii=0
      do j=1,nrayl
         do i=1,neltmax
            ii=ii+1
            tem1(ii)=0.5d0*(cwezde_nc(i,j)+dconjg(cwezde_nc(i,j)))
         enddo
      enddo
      do j=1,nrayl
         do i=1,neltmax
            ii=ii+1
            tem1(ii)=-cei*0.5d0*(cwezde_nc(i,j)-dconjg(cwezde_nc(i,j)))
         enddo
      enddo
      call ncvid2(vid,ncid,'cwezde',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(fluxn_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'fluxn',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(sbtot_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'sbtot',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(sene_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'sene',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(ste_nc,1,nrelt,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'ste',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(salphac_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'salphac',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(salphal_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'salphal',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(sb_r_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'sb_r',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(sb_x_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'sb_x',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(sb_y_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'sb_y',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(sb_z_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'sb_z',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(sb_phi_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'sb_phi',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(wn_r_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'wn_r',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(wn_x_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'wn_x',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(wn_y_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'wn_y',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(wn_z_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'wn_z',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(wn_phi_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'wn_phi',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(vgr_r_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'vgr_r',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(vgr_x_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'vgr_x',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(vgr_y_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'vgr_y',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(vgr_z_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'vgr_z',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(vgr_phi_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'vgr_phi',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(flux_r_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'flux_r',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(flux_z_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'flux_z',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call pack21(flux_phi_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
      call ncvid2(vid,ncid,'flux_phi',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)
c--------------------------------------------------------------
c      CD efficiency along rays
c--------------------------------------------------------
c      write(*,*)'netcdfr3d ionetwo',ionetwo
       if (ionetwo.eq.1) then            
         call pack21(w_eff_nc,1,nrelta,1,nrayl,tem1,neltmax,nrayl)
          call ncvid2(vid,ncid,'w_eff_nc',istatus)
          call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)
       endif
c--------------------------
c     write OX conversion data for i_ox=2 case
c--------------------------
      if (i_ox.eq.2) then

         call ncvid2(vid,ncid,'i_ox_conversion',istatus)
         call ncvpt_int2(ncid,vid,1,nrayl,i_ox_conversion_nc,istatus)

         call ncvid2(vid,ncid,'transm_ox',istatus)
         call ncvpt_doubl2(ncid,vid,1,nrayl,transm_ox_nc,istatus)

         call ncvid2(vid,ncid,'cn_par_optimal',istatus)
         call ncvpt_doubl2(ncid,vid,1,nrayl,cn_par_optimal_nc,istatus)

c         call ncvid2(vid,ncid,'cn_par_optimal',istatus)
c         call ncvpt_doubl2(ncid,vid,1,nrayl,cn_par_optimal_nc,istatus)

         call ncvid2(vid,ncid,'cnpar_ox',istatus)
         call ncvpt_doubl2(ncid,vid,1,nrayl,cnpar_ox_nc,istatus)

         call ncvid2(vid,ncid,'cn_b_gradpsi',istatus)
         call ncvpt_doubl2(ncid,vid,1,nrayl,cn_b_gradpsi_nc,istatus)

      endif !i_ox=2 

c------------------------------------------------------------------
c     for total power absorbed at all reflections at all rays
c------------------------------------------------------------------
      call ncvid2(vid,ncid,'w_tot_pow_absorb_at_refl_nc',istatus)
      call ncvpt_doubl2(ncid,vid,0,0, w_tot_pow_absorb_at_refl_nc,
     +                  istatus)
      
C-----------------------------------------------------------------------
c
cl    4. Close netCDF file
c        Page 53-54 netCDF-2 Manual
c
      call ncclos3(ncid,istatus)
      call check_err(istatus)

      endif                    ! End write netcdf data

cSAP090903
      deallocate (tem1,STAT=istat)

      return
      end
c
c




C=======================================================================
C=======================================================================



      subroutine wrtnetcdf_prof(netcdfnml,kopt)
      implicit integer (i-n), real*8 (a-h,o-z)

c     If the namelist ionetwo.eq.1 it will write the current
c     and power profiles into existing netcdf file: netcdfnml 

c-----profiles description:

c     spower(i)   i=1,NR-1, power [erg/sec] in small radius bin(i)
c     powden(i)  i=1,NR-1, power density [erg/(cm**3*sec)]
c                           in small radius bin(i)
c
c     s_cur_den_parallel(i) i=1,NR-1 parallel averaged current density
c                           [A/cm**2] 
c               <j_parallel>=Integral{dl_poloidal*j_parallel/B_poloidal}
c     s_cur_den_onetwo(i) i=1,NR-1 ONETWO current density 
c                        <j^.B^>/B_0=<>j_parallel*modB>/B_0=
c                        =<j_parallel><B**2>/(B_0*<B>)
c     s_cur_den_toroidal(i) i=1,NR-1 toroidal current density
c                         j_phi=<j_parallel>f<1/R**2>/(<B><1/R>)
c     s_cur_den_poloidal(i) i=1,NR-1   poloidal current density
c                         i_poloidal=<j_parallel>B_poloidal/<B>
c                         B_poloidal is taken at theta_poloidal=0
c

c-----total: absorbed power and toroidal current
c     total power = power_total [erg/sec]
c     total toroidal current = tor_curr_total [A]
c     do i=1,NR-1
c        power_total=power_total+spower(i)
c        tor_curr_total=tor_curr_total+scurrent(i)
c     enndo

c    YuP[Nov-2014] Added: Saving Local absorbed power over (R,Z) grid.
c   spwr_rz_e(NRgrid,NZgrid)=  coll-less damping for e
c   spwr_rz_i(NRgrid,NZgrid)        !coll-less damping for i
c   spwr_rz_cl(NRgrid,NZgrid)       !collisional damping
c    (R,Z) grid for power deposition profiles:
c   Rgrid(NRgrid)
c   Zgrid(NZgrid)

      include 'param.i'     
      include 'one.i'
      include 'three.i'
      include 'fourb.i'
      include 'five.i' ! contains rmax
      include 'onetwo.i'  !For profiles and total current.
      include 'rho.i'     !For areatot,voltot,torftot,totlength
cSm070201
      include 'writencdf.i'     !
c --- include file for netCDF declarations 
c --- (obtained from NetCDF distribution)
      include 'netcdf.inc'

c --- some stuff for netCDF file ---
      integer ncid,vid,istatus,ncvdef2,ncdid2,ncddef2
      integer nrhoid,nrhomid,nbulkid,nbulkmid
      integer i,rdims(2)
      integer start(2),starts(2),count(2),counts(2)
      
      integer NRgrid_id, NZgrid_id, rz_dims(2), count_rz(2)

      character(*) netcdfnml ! input filename
      
      real*8 tem1(NR*nbulk) 

      data start/1,1/, starts/1,1/

      save

c      write(*,*)'netcdfr3d.f in wrtnetcdf_prof start=',start

c      write(*,*)'wrtnetcdf_prof,  ionetwo=',ionetwo

c      write(*,*)'wrtnetcdf_prof: powtot_e,powtot_i',powtot_e,powtot_i

      if (ionetwo.ne.1) return    !nothing to do
      
c.......................................................................
cl    1. Initialize part, open netcdf file, define dims
c

c --- begin if ---
      if ( kopt.eq.1 ) then
c      write(*,*)'wrtnetcdf_prof,  kopt=',kopt

C-----------------------------------------------------------------------
c
c     1.1 open existing netCDF file: netcdfnml, 
c         and define additional;
c         dimensions,variables and attributes
c

c.......................................................................
c     1.1.1 open existing netCDF file: netcdfnml,
c     create netCDF filename (Entering define mode.)
c     integer function ncopn(filename,write,error_code)
c     Ref to page 40 of NetCDF-2 manual.
c     CLOBber old file, if it exists.
c     istatus is 0, if no errors.

      call ncopn2(netcdfnml,NCWRITE,ncid,istatus)    ! Open existing netCDF file
      call check_err(istatus)
      
c.......................................................................
c     1.1.2 
c     Put Open NetCDF file into Define Mode
c     subroutine ncredf3(integer ncid,integer error_code)
c     p. 41 of netcdf manual
c
c     call ncredf3(integer ncid,integer error_code)

      call ncredf3(ncid,error_code)
      call check_err(istatus)

c-------------------------------------------------------------------
c     1.1.3 define dimensions
c     integer function ncddef(ncid,character(*) dim_name,
c                             integer cdim_siz, integer error_code)
c     returns dimension id.
c     p. 52 of netcdf manual

      
c-----For radial profiles data:

c      write(*,*)'netcdfr3d.f in wrtnetcdf_prof 2'



      ! Get dimension ID from dimension name
      nbulkid=ncdid2(ncid,'nbulk',istatus)
      call check_err(istatus)
      nrhoid=ncdid2(ncid,'nrho',istatus)
      call check_err(istatus)
      nrhomid=ncdid2(ncid,'nrhom',istatus)
      call check_err(istatus)
      
      if (nbulk.gt.1) then
         nbulkmid=ncdid2(ncid,'nbulkm',istatus)
         call check_err(istatus)
         rdims(1)=nrhomid  ! (NR-1)    ! only used for 'powden_s'
         rdims(2)=nbulkmid ! (nbulk-1) ! only used for 'powden_s'
         !write(*,*)'netcdfr3d: rdims=',rdims
      endif
      
      count(1)=NR-1    ! only used for 'powden_s'
      count(2)=nbulk-1 ! only used for 'powden_s'
      counts(1)=NR
      counts(2)=nbulk
      
      !---------------------- 
      ! YuP[Nov-2014] Added: Saving Local absorbed power over (R,Z) grid.
      ! Get dimension ID from dimension name
      NRgrid_id=ncdid2(ncid,'NRgrid',istatus)
      call check_err(istatus)
      NZgrid_id=ncdid2(ncid,'NZgrid',istatus)
      call check_err(istatus)
      rz_dims(1)= NRgrid_id     
      rz_dims(2)= NZgrid_id     
      count_rz(1)=NRgrid
      count_rz(2)=NZgrid
      !----------------------
      vid=ncvdef2(ncid,'Rgrid',NCDOUBLE,1,rz_dims(1),istatus) !--- Rgrid
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,38,
     +     'R grid for Power deposition over (R,Z)',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,1,'m',istatus)
      call check_err(istatus)
      vid=ncvdef2(ncid,'Zgrid',NCDOUBLE,1,rz_dims(2),istatus) !--- Zgrid
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,38,
     +     'Z grid for Power deposition over (R,Z)',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,1,'m',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'spwr_rz_e',NCDOUBLE,2,rz_dims,istatus)   !--- e
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,51,
     +    'Power deposited to electrons, locally at (R,Z) grid',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +     'erg/sec',istatus)
      call check_err(istatus)
      vid=ncvdef2(ncid,'spwr_rz_i',NCDOUBLE,2,rz_dims,istatus)   !--- i
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,46,
     +    'Power deposited to ions, locally at (R,Z) grid',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +     'erg/sec',istatus)
      call check_err(istatus)
      vid=ncvdef2(ncid,'spwr_rz_cl',NCDOUBLE,2,rz_dims,istatus)  !--- cl
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,50,
     +    'Collisional Power, deposited locally at (R,Z) grid',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +     'erg/sec',istatus)
      call check_err(istatus)
      !----------------------
      
      
c      write(*,*)'netcdfr3d.f in wrtnetcdf_prof 4'
      
c    
c     For toroidal current and power:
c-----define variables and attributes
      vid=ncvdef2(ncid,'parallel_cur_total',NCDOUBLE,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,22,
     +            'Total parallel current',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,1,
     +           'A',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'toroidal_cur_total',NCDOUBLE,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,22,
     +            'Total toroidal current',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,1,
     +           'A',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'poloidal_cur_total',NCDOUBLE,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,22,
     +            'Total poloidal current',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,1,
     +           'A',istatus)
      call check_err(istatus)
      
      vid=ncvdef2(ncid,'power_inj_total',NCDOUBLE,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,20,
     +            'Total injected power',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +           'erg/sec',istatus)
      call check_err(istatus)
      
      vid=ncvdef2(ncid,'power_total',NCDOUBLE,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,20,
     +            'Total absorbed power',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +           'erg/sec',istatus)
      call check_err(istatus)
      
      vid=ncvdef2(ncid,'powtot_e',NCDOUBLE,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,24,
     +            'Total power to electrons',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +           'erg/sec',istatus)
      call check_err(istatus)
      
      vid=ncvdef2(ncid,'powtot_i',NCDOUBLE,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,19,
     +            'Total power to ions',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +           'erg/sec',istatus)
      call check_err(istatus)
      
      vid=ncvdef2(ncid,'powtot_cl',NCDOUBLE,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,26,
     +            'Collisional power absorbed',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +           'erg/sec',istatus)
      call check_err(istatus)

      if (iabsorp.eq.3. and. nbulk.gt.1) then !YuP[Nov-2014] why not for any iabsorp?
         vid=ncvdef2(ncid,'powtot_s',NCDOUBLE,1,nbulkmid,istatus)
         call check_err(istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +        'Power to individual ion species',istatus)
         call check_err(istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +        'erg/sec',istatus)
         call check_err(istatus)
      endif         

c.......................................................................
cl    1.1.4 define variables

c     Note, the variable IDs (denoted below as vid) are
c     not saved here in this subroutine; rather, the IDs
c     are retrieved from the netCDF data file, as needed,
c     by calling the netCDF routine ncvid.

c     netCDF variable_define:
c     integer function ncvdef(ncid,variable_name,variable_type,
c                number_of_dimensions,
c                vector_for_length_per_dimension,error_code)
c       returns varid.
c     Note: Unlimited dimension must be last
c     Refer to p. 77, netcdf-2 manual

c     NCDOUBLE for REAL*8.  This is independent of the
c       internal representation of data on a particular
c       machine, e.g., 32- or 64-bit based calculations.


c--------------------------
c     current and power profiles data
c--------------------------
      vid=ncvdef2(ncid,'NR',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,29,
     +            'Number of small radius points',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'voltot',NCDOUBLE,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,40,
     +            'Total volume in last closed flux surface',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,4,
     +           'cm^3',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'areatot',NCDOUBLE,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,34,
     +            'Total cross-sectional area of LCFS',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,4,
     +           'cm^2',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'pollentot',NCDOUBLE,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23,
     +            'Poloidal length of LCFS',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,2,
     +           'cm',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'torftot',NCDOUBLE,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,26,
     +            'Toroidal flux through LCFS',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,9,
     +           'Tesla*m^2',istatus)
      call check_err(istatus)
 
      vid=ncvdef2(ncid,'indexrho',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,40,
     +           'Radial coord type: 2 gives sqrt(tor flx)',istatus)
      call check_err(istatus)


      vid=ncvdef2(ncid,'psifactr',NCDOUBLE,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +            'Reduces Psi-Value of LCFS',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'Should be .le.1',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'binvol',NCDOUBLE,1,nrhomid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,22,
     +            'Volumes of radial bins',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,4,
     +           'cm^3',istatus)
      call check_err(istatus)
 
      vid=ncvdef2(ncid,'binarea',NCDOUBLE,1,nrhomid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,20,
     +            'Areas of radial bins',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,4,
     +           'cm^2',istatus)
      call check_err(istatus)
 
      vid=ncvdef2(ncid,'pollen',NCDOUBLE,1,nrhomid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,31,
     +            'Poloidal lengths of radial bins',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,2,
     +           'cm',istatus)
      call check_err(istatus)
 
      vid=ncvdef2(ncid,'spower',NCDOUBLE,1,nrhomid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,21,
     +            'power per bin profile',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +           'erg/sec',istatus)
      call check_err(istatus)

c----------------------------------------------------------------
c     averaged current densities according to GA-memo
c----------------------------------------------------------------
      vid=ncvdef2(ncid,'GA_tor_cur_total',NCDOUBLE,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +            'GA memo total toroidal current',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,1,
     +           'A',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'s_cur_den_parallel',NCDOUBLE,1,nrhomid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,41,
     &            'averaged parallel current density profile',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     &           'A/cm**2',istatus)
      call check_err(istatus)


      vid=ncvdef2(ncid,'s_cur_den_onetwo',NCDOUBLE,1,nrhomid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,40,
     &            'ONETWO current <j.B>/B_0 density profile',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     &           'A/cm**2',istatus)
      call check_err(istatus)


      vid=ncvdef2(ncid,'s_cur_den_toroidal',NCDOUBLE,1,nrhomid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,60,
     &'toroidal current <j_par>f<1/r**2>/(<B><1/r>) density profile',
     &istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +           'A/cm**2',istatus)
      call check_err(istatus)


      vid=ncvdef2(ncid,'s_cur_den_poloidal',NCDOUBLE,1,nrhomid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,62,
     &'poloidal current <j_par>B_pol(theta_pol=0)/<B> density profile',
     &istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +           'A/cm**2',istatus)
      call check_err(istatus)

c-------------------------------------------------------

      vid=ncvdef2(ncid,'powden',NCDOUBLE,1,nrhomid,istatus)
      call check_err(istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,21,
     +            'power density profile',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'erg/(cm**3*sec)',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'powden_e',NCDOUBLE,1,nrhomid,istatus)
      call check_err(istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,48,
     +     'power density profile to electrons, bin centered',istatus)  
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'erg/(cm**3*sec)',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'powden_i',NCDOUBLE,1,nrhomid,istatus)
      call check_err(istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,43,
     +     'power density profile to ions, bin centered',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'erg/(cm**3*sec)',istatus)
      call check_err(istatus)

      if (iabsorp.eq.3) then !YuP[Nov-2014] why not for any iabsorp?
         vid=ncvdef2(ncid,'powden_s',NCDOUBLE,2,rdims,istatus)
         write(*,*)'netcdfr3d: istatus after powden_s def. rdims=',rdims
         call check_err(istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,57,
     +      'power density profile to individual ion species, bin cent',
     +      istatus)
         call check_err(istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +        'erg/(cm**3*sec)',istatus)
         call check_err(istatus)
      endif

      vid=ncvdef2(ncid,'powden_cl',NCDOUBLE,1,nrhomid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,31,
     +     'collisional power, bin centered',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'erg/(cm**3*sec)',istatus)
      call check_err(istatus)
    
c.......................................................................
cl    1.1.5 end the define-mode and start the data-mode
c     p. 51-2 of manual

      call ncendf3(ncid,istatus)
      call check_err(istatus)

      endif               ! End initialize, kopt=1


c.......................................................................
cl    1. Writing data
c
      

      if ( kopt.eq.0 ) then
c      write(*,*)'wrtnetcdf_prof,  kopt,NR=',kopt,NR
cl    1.2 Write data
c

c     First get variable_id:
c             nf_inq_varid(netcdf_id,variable_name,integer_info)
c     Then write data with nc variable_put:
c             nf_put_var... ()
c
c     function nf_put_var...(ncid,variable_id,index,val)
c     p. 52,, 54-65 of netcdf-3 manual
c
c--------------------------
         
c-----calculate total toroidal current and total absorbed power
      power_total=0.d0
      tor_curr_total=0.d0
      do i=1,NR-1
         power_total=power_total+spower(i)              
      enddo


c--------------------------
c     write
c     the number of radial points NR for bin boundaries  and 
c     total toroidal current and total power data
c--------------------------     
c      write(*,*)'in wrtnetcdf_prof kopt===0  pt.1'
      call ncvid2(vid,ncid,'NR',istatus)
      call ncvpt_int2(ncid,vid,0,0,NR,istatus)
c      write(*,*)'in wrtnetcdf_prof kopt===0  pt.2'
      call ncvid2(vid,ncid,'voltot',istatus)
      call ncvpt_doubl2(ncid,vid,0,0,1.d6*voltot,istatus)

      call ncvid2(vid,ncid,'areatot',istatus)
      call ncvpt_doubl2(ncid,vid,0,0,1.d4*areatot,istatus)

      call ncvid2(vid,ncid,'pollentot',istatus)
      call ncvpt_doubl2(ncid,vid,0,0,1.d2*totlength,istatus)

      call ncvid2(vid,ncid,'torftot',istatus)
      call ncvpt_doubl2(ncid,vid,0,0,torftot,istatus)

      call ncvid2(vid,ncid,'indexrho',istatus)
      call ncvpt_int2(ncid,vid,0,0,indexrho,istatus)

      call ncvid2(vid,ncid,'psifactr',istatus)
      call ncvpt_doubl2(ncid,vid,0,0,psifactr,istatus)

      call ncvid2(vid,ncid,'binvol',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR-1,binvol,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'binarea',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR-1,binarea,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'pollen',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR-1,pollen,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'parallel_cur_total',istatus)
      call ncvpt_doubl2(ncid,vid,0,0,parallel_cur_total,istatus)

      call ncvid2(vid,ncid,'toroidal_cur_total',istatus)
      call ncvpt_doubl2(ncid,vid,0,0,toroidal_cur_total,istatus)
   
      call ncvid2(vid,ncid,'poloidal_cur_total',istatus)
      call ncvpt_doubl2(ncid,vid,0,0,poloidal_cur_total,istatus)

      call ncvid2(vid,ncid,'power_inj_total',istatus)
      call ncvpt_doubl2(ncid,vid,0,0,powtott,istatus)

      call ncvid2(vid,ncid,'power_total',istatus)
      call ncvpt_doubl2(ncid,vid,0,0,power_total,istatus)

      call ncvid2(vid,ncid,'powtot_e',istatus)
      call ncvpt_doubl2(ncid,vid,0,0,powtot_e,istatus)

      call ncvid2(vid,ncid,'powtot_i',istatus)
      call ncvpt_doubl2(ncid,vid,0,0,powtot_i,istatus)

      call ncvid2(vid,ncid,'powtot_cl',istatus)
      call ncvpt_doubl2(ncid,vid,0,0,powtot_cl,istatus)

      if (iabsorp.eq.3) then !YuP[Nov-2014] why not for any iabsorp?
         call ncvid2(vid,ncid,'powtot_s',istatus)
         call ncvpt_doubl2(ncid,vid,1,nbulk-1,powtot_s(2),istatus)
         call check_err(istatus)
      endif
c      write(*,*)'in wrtnetcdf_prof kopt===0  pt.3'



c --------------------------
c     current and power profiles data
c--------------------------  
      call ncvid2(vid,ncid,'spower',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR-1,spower,istatus)
      call check_err(istatus)
      
c----------------------------------------------------------------
c     averaged current densities according to GA-memo
c----------------------------------------------------------------
c-----calculate total GA memo toroidal current 
      GA_tor_cur_total=0.d0
c      write(*,*)'in wrtnetcdf_prof'
      do i=1,NR-1
         GA_tor_cur_total=GA_tor_cur_total+s_cur_den_toroidal(i)*
     &                     binarea(i)
c         GA_tor_cur_total=GA_tor_cur_total+s_cur_den_onetwo(i)*
c     &                     binarea(i)

      enddo

c      write(*,*)'GA_tor_cur_total ',GA_tor_cur_total
cSAP090306
cyup      write(*,1002)GA_tor_cur_total     
 1002 format('total toroidal current GA_tor_cur_total (A)=',
     &   1pe14.6)

      call ncvid2(vid,ncid,'GA_tor_cur_total',istatus)
      call ncvpt_doubl2(ncid,vid,0,0,GA_tor_cur_total,istatus)

      call ncvid2(vid,ncid,'s_cur_den_parallel',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR-1,s_cur_den_parallel,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'s_cur_den_onetwo',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR-1,s_cur_den_onetwo,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'s_cur_den_toroidal',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR-1,s_cur_den_toroidal,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'s_cur_den_poloidal',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR-1,s_cur_den_poloidal,istatus)
      call check_err(istatus)


c----------------------------------------------------------------
      call ncvid2(vid,ncid,'powden',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR-1,powden,istatus) 
      ! NR-1 because profiles are defined at rho_bin_center  [rhomid]
      call check_err(istatus)
 
      call ncvid2(vid,ncid,'powden_e',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR-1,powden_e,istatus)
      call check_err(istatus)
 
      call ncvid2(vid,ncid,'powden_i',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR-1,powden_i,istatus)
      call check_err(istatus)
     
      call ncvid2(vid,ncid,'powden_cl',istatus)  
      call ncvpt_doubl2(ncid,vid,1,NR-1,powden_cl,istatus)
      call check_err(istatus)

      if (iabsorp.eq.3 .and. nbulk.gt.1) then ! powden_s
      call pack21(powden_s(1,2),1,NR,1,nbulk-1,tem1,NR-1,nbulk-1)      
      call ncvid2(vid,ncid,'powden_s',istatus)
      call ncvpt_doubl2(ncid,vid,start,count,tem1,istatus)
      write(*,*)'netcdfr3d: istatus after powden_s rec.:', istatus 
      call check_err(istatus)    
      endif
 
      !---------------------- YuP[Nov-2014] Added: power over (R,Z) grid.
      call ncvid2(vid,ncid,'Rgrid',istatus)  !--- Rgrid
      call ncvpt_doubl2(ncid,vid,1,NRgrid,Rgrid,istatus)
      call check_err(istatus)
      call ncvid2(vid,ncid,'Zgrid',istatus)  !--- Zgrid
      call ncvpt_doubl2(ncid,vid,1,NZgrid,Zgrid,istatus)
      call check_err(istatus)
      call ncvid2(vid,ncid,'spwr_rz_e',istatus)  !--- e
      call ncvpt_doubl2(ncid,vid,start,count_rz,spwr_rz_e,istatus)
      call check_err(istatus)
      call ncvid2(vid,ncid,'spwr_rz_i',istatus)  !--- i
      call ncvpt_doubl2(ncid,vid,start,count_rz,spwr_rz_i,istatus)
      call check_err(istatus)
      call ncvid2(vid,ncid,'spwr_rz_cl',istatus) !--- cl
      call ncvpt_doubl2(ncid,vid,start,count_rz,spwr_rz_cl,istatus)
      call check_err(istatus)
      !----------------------


c      do i=1,NR-1
c         write(*,*)'i,scurrent(i),currden(i)',i,scurrent(i),currden(i)
c         write(*,*)'i,spower(i),powden(i)',i,spower(i),powden(i)
c      enddo     

c      write(*,1000)'i rho_bin_center scurrent    currden     spower    '
c     1          //'powden'
c      do i=1,NR-1
c         write(*,1001)i,rho_bin_center(i),scurrent(i),currden(i),
c     1        spower(i),powden(i)
c      enddo
 1000    format(/,1x,a)
 1001    format(i3,5(1pe12.4))



C-----------------------------------------------------------------------
c
cl    4. Close netCDF file
c        Page 53-54 netCDF-2 Manual
c
      call ncclos3(ncid,istatus)
      call check_err(istatus)

      endif                    ! End write netcdf data, kopt=0

      return
      end
c
c

      subroutine wrtnetcdf_eps(netcdfnml)
      implicit integer (i-n), real*8 (a-h,o-z)

c     If the namelist variable dielectric_op="enabled", this routine
c     will write the complex dielectric tensor elements along the rays
c     into into existing netcdf file: netcdfnml 


      include 'param.i'     
      include 'one.i'
      include 'writencdf.i'

      integer n1n2,istat
c     Storage tem1 is used in netcdf writes, including complex numbers.
      real*8, pointer :: tem1(:)

c --- include file for netCDF declarations 
c --- (obtained from NetCDF distribution)
      include 'netcdf.inc'

c --- some stuff for netCDF file ---
      integer ncid,vid,istatus,ncvdef2,ncdid2,ncddef2
      integer nraysid,neltmaxid,twoid
      integer ray_dims(3),start(3),ray_count(3)

c----input
      character(*) netcdfnml ! input filename

      data start/1,1,1/

cSAP090903
      n1n2=2*nrelta*nrayl
c------------------------------------------
c     allocate pointers tem1
c-------------------------------------------
      allocate( tem1(1:n1n2),STAT=istat)
      call bcast(tem1,0.d0,SIZE(tem1))
c----------------------------------------------
 
c     Maximum number of ray elements per ray:
      neltmax=0
      do iray=1,nrayl
         neltmax=max(neltmax,nrayelt_nc(iray))
      enddo
      
      ray_count(1)=neltmax
      ray_count(2)=nrayl
      ray_count(3)=2

c.......................................................................
cl    1. Initialize part
c

C-----------------------------------------------------------------------
c
c     1.1 open existing netCDF file: netcdfnml, 
c         and define additional;
c         dimensions,variables and attributes
c

c.......................................................................
c     1.1.1 open existing netCDF file: netcdfnml,
c     create netCDF filename (Entering define mode.)
c     integer function ncopn(filename,write,error_code)
c     Ref to page 40 of NetCDF-2 manual.
c     CLOBber old file, if it exists.
c     istatus is 0, if no errors.
      
      call ncopn2(netcdfnml,NCWRITE,ncid,istatus) ! Open existing netCDF file
      call check_err(istatus)
      
c.......................................................................
c     1.1.2 
c     Put Open NetCDF file into Define Mode
c     subroutine ncredf3(integer ncid,integer error_code)
c     p. 41 of netcdf manual
c
c     call ncredf3(integer ncid,integer error_code)
     
      call ncredf3(ncid,error_code)
      call check_err(istatus)

c-------------------------------------------------------------------
c     1.1.3 determine dimension ID of 
c     dimension named 'lat'
c     assumed to have been defined previosly in an existing netCDF file
c     with ID ncid 
c
c     integer function ncdid(ncid,character(*) lat,integer error_code)
c     returns dimension id.
c     p. 54 of netcdf manual



      neltmaxid=ncdid2(ncid,'neltmax',istatus)
      nraysid=ncdid2(ncid,'nrays',istatus)
      twoid=ncdid2(ncid,'two',istatus)
     
      ray_dims(1)=neltmaxid
      ray_dims(2)=nraysid
      ray_dims(3)=twoid

c     For ray data:
c-----define variables 

c-----Added 3rd dimension equal to 2 accomodates complex data.
  
      vid=ncvdef2(ncid,'cweps11',NCDOUBLE,3,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,13,
     +           'Complex eps11',istatus)

      vid=ncvdef2(ncid,'cweps12',NCDOUBLE,3,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,13,
     +           'Complex eps12',istatus)  

      vid=ncvdef2(ncid,'cweps13',NCDOUBLE,3,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,13,
     +           'Complex eps13',istatus)

      vid=ncvdef2(ncid,'cweps21',NCDOUBLE,3,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,13,
     +           'Complex eps21',istatus)

      vid=ncvdef2(ncid,'cweps22',NCDOUBLE,3,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,13,
     +           'Complex eps22',istatus)  

      vid=ncvdef2(ncid,'cweps23',NCDOUBLE,3,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,13,
     +           'Complex eps23',istatus)

      vid=ncvdef2(ncid,'cweps31',NCDOUBLE,3,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,13,
     +           'Complex eps31',istatus)

      vid=ncvdef2(ncid,'cweps32',NCDOUBLE,3,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,13,
     +           'Complex eps32',istatus)  

      vid=ncvdef2(ncid,'cweps33',NCDOUBLE,3,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,13,
     +           'Complex eps33',istatus)
     
c.......................................................................
cl    1.1.5 end the define-mode and start the data-mode
c     p. 51-2 of manual

      call ncendf3(ncid,istatus)
      call check_err(istatus)
   
c     End initialize

c.......................................................................
cl    1. Writing data
c
 
cl    1.2 Write data
c

c     First get variable_id:
c             nf_inq_varid(netcdf_id,variable_name,integer_info)
c     Then write data with nc variable_put:
c             nf_put_var... ()
c
c     function nf_put_var...(ncid,variable_id,index,val)
c     p. 52,, 54-65 of netcdf-3 manual
c 
c--------------------------
c     write eps data
c--------------------------     
      call storage_compl_2d(nrelta,nrayl,nrayl,neltmax,
     &cweps11_nc,tem1)    
      call ncvid2(vid,ncid,'cweps11',istatus)     
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call storage_compl_2d(nrelta,nrayl,nrayl,neltmax,
     &cweps12_nc,tem1)
      call ncvid2(vid,ncid,'cweps12',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call storage_compl_2d(nrelta,nrayl,nrayl,neltmax,
     &cweps13_nc,tem1)
      call ncvid2(vid,ncid,'cweps13',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call storage_compl_2d(nrelta,nrayl,nrayl,neltmax,
     &cweps21_nc,tem1)
      call ncvid2(vid,ncid,'cweps21',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call storage_compl_2d(nrelta,nrayl,nrayl,neltmax,
     &cweps22_nc,tem1)
      call ncvid2(vid,ncid,'cweps22',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call storage_compl_2d(nrelta,nrayl,nrayl,neltmax,
     &cweps23_nc,tem1)
      call ncvid2(vid,ncid,'cweps23',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call storage_compl_2d(nrelta,nrayl,nrayl,neltmax,
     &cweps31_nc,tem1)
      call ncvid2(vid,ncid,'cweps31',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call storage_compl_2d(nrelta,nrayl,nrayl,neltmax,
     &cweps32_nc,tem1)
      call ncvid2(vid,ncid,'cweps32',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      call storage_compl_2d(nrelta,nrayl,nrayl,neltmax,
     &cweps33_nc,tem1)
      call ncvid2(vid,ncid,'cweps33',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cl    4. Close netCDF file
c        Page 53-54 netCDF-2 Manual
c
      call ncclos3(ncid,istatus)
      call check_err(istatus)

c     End write netcdf data

      deallocate( tem1,STAT=istat)

      return
      end
c
c


      subroutine storage_compl_2d(nrelta,nraya,nrayl,neltmax,
     &compl_2d,tem1)
c     put complex*16 2D array: compl_2d(neltmax,nrayl)
c     to real*8 1D array: tem1
c
c     It will be used in  netcdf subroutines

      
      implicit none

c-----input
      integer nrelta,nraya,nrayl,neltmax
      complex*16 compl_2d(nrelta,nraya)

c-----output
c     Storage tem1 is used in netcdf writes, including complex numbers.
      real*8 tem1(*) 

c-----locals    
      integer j,i,ii
      complex*16 cei

      cei=(0.d0,1.d0)

      ii=0
      do j=1,nrayl
         do i=1,neltmax
            ii=ii+1
            tem1(ii)=0.5d0*(compl_2d(i,j)+dconjg(compl_2d(i,j)))
         enddo
      enddo

      do j=1,nrayl
         do i=1,neltmax
            ii=ii+1
            tem1(ii)=-cei*0.5d0*(compl_2d(i,j)-dconjg(compl_2d(i,j)))
         enddo
      enddo

      return
      end



      subroutine wrtnetcdf_grill_launch(netcdfnml)
      implicit integer (i-n), real*8 (a-h,o-z)

c     If the namelist variable istart=2 or =3  this routine
c     will write the ray starting coordinates 
c     into existing netcdf file: netcdfnml 



      include 'param.i'     
      include 'one.i'
      include 'writencdf.i'
      include 'grill.i'
     
c --- include file for netCDF declarations 
c --- (obtained from NetCDF distribution)
      include 'netcdf.inc'

c --- some stuff for netCDF file ---
      integer ncid,vid,istatus,ncvdef2,ncdid2,ncddef2
      integer nraysid,neltmaxid
      
c----input
      character(*) netcdfnml ! input filename

c.......................................................................
c     1. Initialize part
c
C-----------------------------------------------------------------------
c
c     1.1 open existing netCDF file: netcdfnml, 
c         and define additional;
c         dimensions,variables and attributes
c
c.......................................................................
c     1.1.1 open existing netCDF file: netcdfnml,
c     create netCDF filename (Entering define mode.)
c     integer function ncopn(filename,write,error_code)
c     Ref to page 40 of NetCDF-2 manual.
c     CLOBber old file, if it exists.
c     istatus is 0, if no errors.
      
c      write(*,*)'netcdfr3d.f in wrtnetcdf_grill_launch'
c      write(*,*)'before ncid=ncopn(netcdfnml,NCWRITE,istatus)'
c      write(*,*)'netcdfnml = ',netcdfnml

      call ncopn2(netcdfnml,NCWRITE,ncid,istatus) ! Open existing netCDF file
      call check_err(istatus)

c      write(*,*)'after ncid=ncopn2(netcdfnml,NCWRITE,) ncid= ',ncid
c.......................................................................
c     1.1.2 
c     Put Open NetCDF file into Define Mode
c     subroutine ncredf3(integer ncid,integer error_code)
c     p. 41 of netcdf manual
c
c     call ncredf3(integer ncid,integer error_code)
     
      call ncredf3(ncid,error_code)
      call check_err(istatus)

c-------------------------------------------------------------------
c     1.1.3 determine dimension ID of 
c     dimension named 'lat'
c     assumed to have been defined previosly in an existing netCDF file
c     with ID ncid 
c
c     integer function ncdid(ncid,character(*) lat,integer error_code)
c     returns dimension id.
c     p. 54 of netcdf manual
     
      nraysid=ncdid2(ncid,'nrays',istatus)


c     For ray data:
c-----define variables 

c-----starting coordinate.
c      write(*,*)'before ncvdef2(z_starting) ncid,nraysid',ncid,nraysid
      vid=ncvdef2(ncid,'z_starting',NCDOUBLE,1,nraysid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,10,
     +           'z_starting',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,1,
     +           'm',istatus)

c      write(*,*)'before ncvdef2(r_starting) ncid,nraysid',ncid,nraysid
      vid=ncvdef2(ncid,'r_starting',NCDOUBLE,1,nraysid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,10,
     +           'r_starting',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,1,
     +           'm',istatus)

c      write(*,*)'before ncvdef2(phi_strting) ncid,nraysid',ncid,nraysid
      vid=ncvdef2(ncid,'phi_starting',NCDOUBLE,1,nraysid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,12,
     +           'phi_starting',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,6,
     +           'radians',istatus)

c      write(*,*)'before ncvdef2(N_toroidal_starting) ncid,nraysid',
c     +           ncid,nraysid
      vid=ncvdef2(ncid,'N_toroidal_starting',NCDOUBLE,1,nraysid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,19,
     +           'N_toroidal_starting',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'non-dimensional',istatus)

c      write(*,*)'before ncvdef(N_pol_starting)ncid,nraysid',ncid,nraysid
      vid=ncvdef2(ncid,'N_poloidal_starting',NCDOUBLE,1,nraysid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,19,
     +           'N_poloidal_starting',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'non-dimensional',istatus)

c.......................................................................
cl    1.1.5 end the define-mode and start the data-mode
c     p. 51-2 of manual

      call ncendf3(ncid,istatus)
      call check_err(istatus)

c     End initialize

c.......................................................................
cl    1. Writing data
c
 
cl    1.2 Write data
c

c     First get variable_id:
c             nf_inq_varid(netcdf_id,variable_name,integer_info)
c     Then write data with nc variable_put:
c             nf_put_var... ()
c
c     function nf_put_var...(ncid,variable_id,index,val)
c     p. 52,, 54-65 of netcdf-3 manual
c 
c--------------------------
c     write ray starting coordinates data for all rays
c--------------------------     
c      write(*,*)'before ncvid( z_starting ) ncid',ncid
      call ncvid2(vid,ncid,'z_starting',istatus)     
c      write(*,*)'after ncvid( z_starting vid',vid
      call ncvpt_doubl2(ncid,vid,1,nrayl,arzu0,istatus)

c      write(*,*)'before ncvid( r_starting ) ncid',ncid
      call ncvid2(vid,ncid,'r_starting',istatus)     
      call ncvpt_doubl2(ncid,vid,1,nrayl,arru0,istatus)
           
c      write(*,*)'before ncvid( phi_starting ) ncid',ncid
      call ncvid2(vid,ncid,'phi_starting',istatus)     
      call ncvpt_doubl2(ncid,vid,1,nrayl,arphiu0,istatus)
   
c      write(*,*)'before ncvid( N_toroidal_starting ) ncid',ncid
      call ncvid2(vid,ncid,'N_toroidal_starting',istatus)     
      call ncvpt_doubl2(ncid,vid,1,nrayl,arnphi,istatus)

c      write(*,*)'before ncvid( N_poloidal_starting ncid',ncid
      call ncvid2(vid,ncid,'N_poloidal_starting',istatus)     
      call ncvpt_doubl2(ncid,vid,1,nrayl,arntheta,istatus)

cl    4. Close netCDF file
c        Page 53-54 netCDF-2 Manual
c
      call ncclos3(ncid,istatus)
      call check_err(istatus)

c     End write netcdf data

      return
      end



      subroutine wrtnetcdf_EC_launch(netcdfnml)
      implicit integer (i-n), real*8 (a-h,o-z)

c     If the namelist variable istart=1  this routine
c     will write the ray starting coordinates 
c     into existing netcdf file: netcdfnml 

      include 'param.i'     
      include 'one.i'
      include 'writencdf.i'
      include 'cone.i'
     
c --- include file for netCDF declarations 
c --- (obtained from NetCDF distribution)
      include 'netcdf.inc'

c --- some stuff for netCDF file ---
      integer ncid,vid,istatus,ncvdef2,ncdid2,ncddef2
      integer nraysid,neltmaxid
      
c----input
      character(*) netcdfnml ! input filename
      
      real*8 wk(nraymax)

c.......................................................................
c     1. Initialize part
c
C-----------------------------------------------------------------------
c
c     1.1 open existing netCDF file: netcdfnml, 
c         and define additional;
c         dimensions,variables and attributes
c
c.......................................................................
c     1.1.1 open existing netCDF file: netcdfnml,
c     create netCDF filename (Entering define mode.)
c     integer function ncopn(filename,write,error_code)
c     Ref to page 40 of NetCDF-2 manual.
c     CLOBber old file, if it exists.
c     istatus is 0, if no errors.
      
c      write(*,*)'netcdfr3d.f in wrtnetcdf_grill_launch'
c      write(*,*)'before ncid=ncopn2(netcdfnml,NCWRITE,istatus)'
c      write(*,*)'netcdfnml = ',netcdfnml

      call ncopn2(netcdfnml,NCWRITE,ncid,istatus) ! Open existing netCDF file
      call check_err(istatus)

c      write(*,*)'after ncid=ncopn2(netcdfnml,NCWRITE,) ncid= ',ncid
c.......................................................................
c     1.1.2 
c     Put Open NetCDF file into Define Mode
c     subroutine ncredf3(integer ncid,integer error_code)
c     p. 41 of netcdf manual
c
c     call ncredf3(integer ncid,integer error_code)
     
      call ncredf3(ncid,error_code)
      call check_err(istatus)

c-------------------------------------------------------------------
c     1.1.3 determine dimension ID of 
c     dimension named 'lat'
c     assumed to have been defined previosly in an existing netCDF file
c     with ID ncid 
c
c     integer function ncdid(ncid,character(*) lat,integer error_code)
c     returns dimension id.
c     p. 54 of netcdf manual

c      write(*,*)'before nraysid=ncdid(ncid,nrays,istatus)'

      nraysid=ncdid2(ncid,'nrays',istatus)

c      write(*,*)'after nraysid=ncdid(ncid,nrays,) nraysid ',nraysid

c     For ray data:
c-----define variables 

c-----starting coordinates.
      vid=ncvdef2(ncid,'z_starting',NCDOUBLE,1,nraysid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,10,
     +           'z_starting',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,1,
     +           'm',istatus)

      vid=ncvdef2(ncid,'r_starting',NCDOUBLE,1,nraysid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,10,
     +           'r_starting',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,1,
     +           'm',istatus)

      vid=ncvdef2(ncid,'phi_starting',NCDOUBLE,1,nraysid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,12,
     +           'phi_starting',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +           'degree',istatus)

      vid=ncvdef2(ncid,'alphast',NCDOUBLE,1,nraysid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,52,
     &'toroidal angle measured from R-vector through source',istatus)     
      call ncaptc2(ncid,vid,'units',NCCHAR,6,'degree',istatus)

      vid=ncvdef2(ncid,'betast',NCDOUBLE,1,nraysid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,73,
     &'poloidal angle measured from z=constant plane,'//
     & ' pos above plane, neg below',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,6,
     +           'degree',istatus)

c.......................................................................
cl    1.1.5 end the define-mode and start the data-mode
c     p. 51-2 of manual

      call ncendf3(ncid,istatus)
      call check_err(istatus)

c     End initialize

c.......................................................................
cl    1. Writing data
c
 
cl    1.2 Write data
c

c     First get variable_id:
c             nf_inq_varid(netcdf_id,variable_name,integer_info)
c     Then write data with nc variable_put:
c             nf_put_var... ()
c
c     function nf_put_var...(ncid,variable_id,index,val)
c     p. 52,, 54-65 of netcdf-3 manual
c
c--------------------------
c     write ray starting coordinates data for all rays
c--------------------------     
      call ncvid2(vid,ncid,'z_starting',istatus)     
      call ncvpt_doubl2(ncid,vid,1,nrayl,zstj,istatus)

      call ncvid2(vid,ncid,'r_starting',istatus)     
      call ncvpt_doubl2(ncid,vid,1,nrayl,rstj,istatus)

      call ncvid2(vid,ncid,'phi_starting',istatus)     
      wk=180.d0/pi*phistj  ! was phist before 03-2016
      call ncvpt_doubl2(ncid,vid,1,nrayl,wk,istatus)
   
      call ncvid2(vid,ncid,'alphast',istatus) 
      wk=180.d0/pi*alphaj
      call ncvpt_doubl2(ncid,vid,1,nrayl,wk,istatus)

      call ncvid2(vid,ncid,'betast',istatus)
      wk=180.d0/pi*betaj
      call ncvpt_doubl2(ncid,vid,1,nrayl,wk,istatus)

cl    4. Close netCDF file
c        Page 53-54 netCDF-2 Manual
c
      call ncclos3(ncid,istatus)
      call check_err(istatus)

c     End write netcdf data

      return
      end




   
      subroutine read_nc(netcdfnm_out)
c-----read output *.nc file to check the resonance curves

c     implicit none
      implicit none
c      implicit integer (i-n), real*8 (a-h,o-z)
      character(*) netcdfnm_out ! input filename
c      character*9 netcdfnm ! input filename

      include 'netcdf.inc'
      include 'param.i' 
      include 'one.i'
      include 'writencdf.i'
      !!! include 'emissa.i'

c      integer n1n2
c      parameter (n1n2=2*nrelta*nfreqa)
c      real*8 tem1(n1n2) 
      integer n1n2
      real*8, pointer :: tem1(:)      !(n1n2)
 
      integer ncid,vid,istatus,ncvdef2,ncdid2,ncddef2

      integer nfreqdim,neltmaxid,neltmax_initial_mesh_dim

      integer ray_dims(2),ray_dims_initial_mesh(2)

      integer start(2),ray_count(2),ray_count_initial_mesh(2)

      character*128 name

      integer i,j,
     &istat,ifreq,neltmax,neltmax_id,nfreq_id,
     & neltmax_emis_initial_mesh
       
      data start/1,1/
       
c      write(*,*)'********************************************'
c      write(*,*)'in read_nc'
c      write(*,*)'********************************************'
c      write(*,*)'netcdfnm_out= ',netcdfnm_out
c-------------------------------------------------
c     Open previous netCDF file
c-------------------------------------------------
      call ncopn2(netcdfnm_out,NCNOWRIT,ncid,istatus)
      call check_err(istatus)
c      write(*,*)'ncid',ncid
c------------------------------------------------
c     Put Open NetCDF file into Define Mode
c------------------------------------------------
c      call ncredf3(ncid,error_code)
c      call check_err(istatus)
c----------------------------------------------------
c     determine dimension ID of  dimension named '...'
c     assumed to have been defined previosly in an existing netCDF file
c     with ID ncid 
c     read ID name and sizes
c-----------------------------------------------------
      neltmax_id =ncdid2(ncid,'neltmax',istatus)
      call check_err(istatus)      
      call ncdinq3(ncid,neltmax_id,name,neltmax,istatus)
      call check_err(istatus)   
c      write(*,*)'neltmax_id,neltmax',neltmax_id,neltmax

      nfreq_id =ncdid2(ncid,'nfreq',istatus)
      call check_err(istatus)   
      call ncdinq3(ncid,nfreq_id,name,nfreq,istatus)
      call check_err(istatus)  
c      write(*,*)'nfreq_id,nfreq',nfreq_id,nfreq
     
c-----allocate tem1
      n1n2=nrelta*nfreq
      allocate( tem1(1:n1n2),STAT=istat)
      call bcast(tem1,0.d0,SIZE(tem1))
c--------------------------------------
      neltmax_initial_mesh_dim=ncdid2(ncid,'neltmax_emis_initial_mesh',
     &istatus)
      call check_err(istatus)
      call ncdinq3(ncid,neltmax_initial_mesh_dim,name,
     &neltmax_emis_initial_mesh,istatus)
      call check_err(istatus)   
c      write(*,*)'neltmax_emis_initial_mesh',neltmax_emis_initial_mesh

      ray_count_initial_mesh(1)=neltmax_emis_initial_mesh
      ray_count_initial_mesh(2)=nfreq

      ray_dims(1)=neltmax
      ray_dims(2)=nfreq

      ray_dims_initial_mesh(1)=neltmax_emis_initial_mesh
      ray_dims_initial_mesh(2)=nfreq
   
c-----------------------------------------------------------------------
c     read numbers of non-refindable mesh points
c     for all frequencies
c     nrayelt_emis_initial_mesh_nc(ifreq), ifreq=1,nfreq
c-----------------------------------------------------------------
      call ncvid2(vid,ncid,'nrayelt_emis_initial_mesh_nc',istatus)
      call ncvgt_int2(ncid,vid,1,nfreq, nrayelt_emis_initial_mesh_nc,
     +                istatus)

c-----------------------------------------------------------------------
c     read the major radus initial (non-refindable) mesh
c     along rays for all frequencies 
c     wr_emis_initial_mesh_nc(i,2,ifreq)
c     Here
c     i=1,nrayelt_emis_initial_mesh_nc(ifreq)
c     ifreq=1,nfreq
c----------------------------------------------------------------------
      call ncvid2(vid,ncid,'wr_emis_initial_mesh_nc',istatus)
      call check_err(istatus)  
      call ncvgt_doubl2(ncid,vid,start, ray_count_initial_mesh,tem1,
     +                  istatus)
      call check_err(istatus)  
      
      j=0
      do ifreq=1,nfreq       
         do i=1,neltmax_emis_initial_mesh
            j=j+1            
            wr_emis_initial_mesh_nc(i,ifreq)=tem1(j)
         enddo      
      enddo      
c-----------------------------------------------------------------------
c     read p_perpmax/me/clight at the initial (non-refindable) mesh
c     along rays for all frequencies for the second harmonic
c     wp_perpmax_dmc_nc(i,2,ifreq)
c     Here
c     i=1,nrayelt_emis_initial_mesh_nc(ifreq)
c     ifreq=1,nfreq
c----------------------------------------------------------------------
      call ncvid2(vid,ncid,'wp_perpmax_dmc_nc',istatus)
      call check_err(istatus)  
      call ncvgt_doubl2(ncid,vid,start, ray_count_initial_mesh,tem1,
     +                  istatus)
      call check_err(istatus)  
       
      j=0
      do ifreq=1,nfreq        
         do i=1,neltmax_emis_initial_mesh
            j=j+1
            wp_perpmax_dmc_nc(i,2,ifreq)=tem1(j)
         enddo      
      enddo      

c-----------------------------------------------------------------------
c     read p_parmax/me/clight and p_parmin/me/cligh
c     at the initial (non-refindable) mesh
c     along rays for all frequencies for the second harmonic
c     wp_parpmax_dmc_nc(i,2,ifreq), wp_parpm_dmc_nc(i,2,ifreq)
c     Here
c     i=1,nrayelt_emis_initial_mesh_nc(ifreq)
c     ifreq=1,nfreq
c----------------------------------------------------------------------
      call ncvid2(vid,ncid,'wp_parmax_dmc_nc',istatus)
      call check_err(istatus)  
      call ncvgt_doubl2(ncid,vid,start,ray_count_initial_mesh,tem1,
     +                  istatus)
      call check_err(istatus)  

      j=0
      do ifreq=1,nfreq
         do i=1,neltmax_emis_initial_mesh
            j=j+1
            wp_parmax_dmc_nc(i,2,ifreq)=tem1(j)
         enddo      
      enddo      

      call ncvid2(vid,ncid,'wp_parmin_dmc_nc',istatus)
      call check_err(istatus)  
      call ncvgt_doubl2(ncid,vid,start,ray_count_initial_mesh,tem1,
     +                  istatus)
      call check_err(istatus)  

      j=0
      do ifreq=1,nfreq
         do i=1,neltmax_emis_initial_mesh
            j=j+1
            wp_parmin_dmc_nc(i,2,ifreq)=tem1(j)
         enddo      
      enddo      
      call ncclos3(ncid,istatus)
      call check_err(istatus)

c--------------------------------------
      deallocate( tem1,STAT=istat)
c---------------------------------
      return
      end



     


      subroutine plasma_profiles_vs_xy(z)
c------------------------------------------------------
c     Creates 1D arrays of density,temperature,zeff
c     versus x and versus y,  at given z
c-------------------------------------------------------
      implicit none
      include 'param.i'
      include 'one.i' ! contains/stores rho
      include 'onetwo_nml.i'
      include 'three.i'
      include 'five.i'
      include 'writencdf.i'
      
c-----input
      double precision z    ! vertical coordinate         
c-----externals
      double precision rhopsi,zeffrho,densrho,temperho
     + ,dense_xyz, bxyz
c-----locals
      double precision 
     &accuracy,  ! accuracy of the solutions
                 ! r=r_min_z ans r=r_max_z
                 ! of the equation psi_limitter-psi(r,z)=0   
                 ! at the given z
     &psi_loc,h_r,rho_loc,z_loc,
     + r, x,y, den

      integer i,k

      accuracy=1.d-4
      z_loc=z
      
      !-1-> Versus x, at y=0
      h_r=(xeqmax-xeqmin)/dfloat(nr-1)
      y= 0.d0 
      do i=1,nr
         x=xeqmin+h_r*(i-1)
         w_x_densprof_nc(i)=x
         w_bmod_vs_x_nc(i)= bxyz(x,y,z) ! to find bmod -> /one.i/
         do k=1,nbulk    
           !-------------- Get value of density and value of rho:
           ! get rho based on model for density (model_rho_dens=1,2,4),
           ! 2D spline of data (model_rho_dens=3),
           ! or based on magnetic flux (model_rho_dens=0,5):
           den=dense_xyz(x,y,z,k) !-> get rho (stored in one.i) 
           rho_loc=rho
           !--------------
           w_dens_vs_x_nc(i,k)=den*1.e13 ! [cm^-3]
           w_temp_vs_x_nc(i,k)=temperho(rho_loc,k)
         enddo
         w_zeff_vs_x_nc(i)=zeffrho(rho_loc)
      enddo

      !-2-> Versus y, at x=0     
      h_r=(yeqmax-yeqmin)/dfloat(nr-1)
      x= 0.d0 
      do i=1,nr
         y=yeqmin+h_r*(i-1)
         w_y_densprof_nc(i)=y
         w_bmod_vs_y_nc(i)= bxyz(x,y,z) ! to find bmod -> /one.i/
         do k=1,nbulk    
           !-------------- Get value of density and value of rho:
           ! get rho based on model for density (model_rho_dens=1,2,4),
           ! 2D spline of data (model_rho_dens=3),
           ! or based on magnetic flux (model_rho_dens=0,5):
           den=dense_xyz(x,y,z,k) !-> get rho (stored in one.i) 
           rho_loc=rho
           !--------------
           w_dens_vs_y_nc(i,k)=den*1.e13 ! [cm^-3]
           w_temp_vs_y_nc(i,k)=temperho(rho_loc,k)
         enddo
         w_zeff_vs_y_nc(i)=zeffrho(rho_loc)
      enddo

      return
      end

c-----------------------------------------------------------------------                           
      subroutine wrtnetcdf_one_ray_point(kopt,is0,nray0)            
      implicit none
c
c     Write ray tracing data into a netCDF file genray_one_ray_point.nc
c     in one ray point:
c     at the ray with number nray0 
c     at the pointwith number is_.

cSm0940727
c     If the parameter ionetwo.eq.1 it will write the current
c     and power profiles into a netcdf file

      include 'param.i'
      include 'writencdf.i'
      include 'one.i'
      include 'ions.i'
      !!! include 'adj.i'
      include 'cone_nml.i'     !nccone
      include 'grill_nml.i'    !ngrill
c-----input
      integer kopt
      integer nray0, !the number of ray
     &is0            !the number of ray point 
 
c-----locals
      character filenc_one_ray_point* 128
cSAP090903
      integer n1n2,istat
c     Storage tem1 is used in netcdf writes, including complex numbers.
      real*8, pointer :: tem1(:)      !(n1n2)
      real*4, pointer :: tem1_single(:)      !(n1n2)
 
cSAP090314
c     Storage to put array(nraya) to one point array(nraya1=1)
c     at one ray 
      integer nrelta1,nraya1
      parameter(nrelta1=1)
      parameter(nraya1=1)
      real*8 point_2D(nrelta1,nraya1)
      real*4 point_3D(nrelta1,nraya1,nbulka) !single prec to save storage
      integer  ipoint_1D(nraya1)
c --- include file for netCDF declarations 
c --- (obtained from NetCDF distribution)
      include 'netcdf.inc'

c --- some stuff for netCDF file ---
      integer ncid,vid,istatus,nccre2,ncvdef2,ncdid2,ncddef2
      integer nraysid,neltmaxid,twoid,char64id,char128id,char512id,
cSAP
     &char8id

    
      integer nbulkmid  !ion species dimension
      integer nbulkid
      integer nrho_adjid

      integer ray_dims(3),start(3),ray_count(3)
      integer ray_dimss(3),starts(3),ray_counts(3)
      integer ray_dimsp(3)
      ! Note: ray_dimsp includes ALL species, while ray_dimss - ions only

      double complex cei

      character ltitle*512

      integer neltmax,iray,i,j,ii,ll,nrayl0
      integer length_char

      data start/1,1,1/,starts/1,1,1/


      save

c      write(*,*)'@@@ wrtnetcdf_one_ray_point begin'

c----------------------------------------------
c     allocate tem1
c----------------------------------------------
      n1n2=2*nrelta*nrayl
      allocate( tem1(1:n1n2),  STAT=istat)
      allocate( tem1_single(1:n1n2), STAT=istat)
      tem1=0.d0
      tem1_single=0.0
c----------------------------------------------

cSAP090309
      filenc_one_ray_point='genray_one_ray_point.nc'

      cei=(0.d0,1.d0)

c     Maximum number of ray elements per ray:
      neltmax=0
cSAP090309 for one point of the ray
      neltmax=1
cSAP090314 number of rays nrayl0
      nrayl0=1  
c      write(*,*)'wrtnetcdf nray=',nrayl

cSm05038
      if (neltmax.eq.0) neltmax=1

      ray_count(1)=neltmax
cSAP090314
c     ray_count(2)=nrayl
      ray_count(2)=nrayl0
      ray_count(3)=2
      ray_counts(1)=neltmax
cSAP090314
c      ray_counts(2)=nrayl
      ray_counts(2)=nrayl0
      ray_counts(3)=1

c      write(*,*)'ray_count',ray_count

c.......................................................................
cl    1. Initialize part, creating new netcdf file
c

c --- begin if ---
      if ( kopt.eq.1 ) then
c      write(*,*)'wrtnetcdf,  kopt=',kopt

C-----------------------------------------------------------------------
c
cl     1.1 create netCDF file and define dimensions,variables
c          and attributes
c

c.......................................................................
cl    1.1.1 create netCDF filename (Entering define mode.)
c     integer function nccre(filename,overwrite?,error_code)
c     Ref to page 46 of NetCDF-2 manual.
c     CLOBber old file, if it exists.
c     istatus is 0, if no errors.

c      write(*,*)'wrtnetcdf before 
c     &ncid=nccre(filenc_one_ray_point,NCCLOB,istatus)'

      ncid=nccre2(filenc_one_ray_point,NCCLOB,istatus)
      call check_err(istatus)
c      write(*,*)'ncid',ncid      
c     Brief description added to file:
      ltitle='netCDF file of ray data from GENRAY version: '//version
      if( length_char(ltitle).gt.512 ) stop 'Adjust ltitle in netcdfrw2'
      call ncaptc2(ncid,NCGLOBAL,'title',NCCHAR,length_char(ltitle),
     +     ltitle,istatus)
        
c.......................................................................
cl    1.1.2 define dimensions
c     integer function ncddef(ncid,dim_name,dim_siz,error_code)
c       returns dimension id.
c     p. 67 of netcdf-2 manual

c     For ray data:
c      write(*,*)'wrtnetcdf before ncddef(neltmax) neltmax=',neltmax

      neltmaxid=ncddef2(ncid,'neltmax',neltmax,istatus)         
     
c      write(*,*)'wrtnetcdf before ncddef(ncid,nrays) nrayl0=',nrayl0
cSAP090314
c      nraysid=ncddef2(ncid,'nrays',nrayl,istatus)
      nraysid=ncddef2(ncid,'nrays',nrayl0,istatus)

c      write(*,*)'wrtnetcdf before ncddef(ncid,two,2,istatus)'
      
      twoid=ncddef2(ncid,'two',2,istatus)
    
c      write(*,*)'wrtnetcdf before ncddef(ncid,nbulk)'

      nbulkid=ncddef2(ncid,'nbulk',nbulk,istatus) 

cSAP080303
c      write(*,*)'wrtnetcdf before ncddef(char8dim)'
      char8id=ncddef2(ncid,'char8dim',8,istatus)

c      write(*,*)'wrtnetcdf before ncddef(char64dim)'

      char64id=ncddef2(ncid,'char64dim',64,istatus)

c      write(*,*)'wrtnetcdf before ncddef(char128dim)'

      char128id=ncddef2(ncid,'char128dim',128,istatus)

c      write(*,*)'wrtnetcdf before ncddef(char5212dim)'

      char512id=ncddef2(ncid,'char512dim',512,istatus)

c      write(*,*)'neltmaxid',neltmaxid
c      write(*,*)'nraysid',nraysid
c      write(*,*)'twoid',twoid

      ray_dims(1)=neltmaxid
      ray_dims(2)=nraysid
      ray_dims(3)=twoid
    
      if (nbulk.gt.1) then
         nbulkmid=ncddef2(ncid,'nbulkm',nbulk-1,istatus)
         ray_dimss(1)=ray_dims(1)
         ray_dimss(2)=ray_dims(2)
         ray_dimss(3)=nbulkmid
      endif
      ! Note: ray_dimsp includes ALL species, while ray_dimss - ions only
      ray_dimsp(1)=ray_dims(1)
      ray_dimsp(2)=ray_dims(2)
      ray_dimsp(3)=nbulkid
      
c.......................................................................
cl    1.1.3 define variables

c     Note, the variable IDs (denoted below as vid) are
c     not saved here in this subroutine; rather, the IDs
c     are retrieved from the netCDF data file, as needed,
c     by calling the netCDF routine ncvid.

c     netCDF variable_define:
c     integer function ncvdef2(ncid,variable_name,variable_type,
c                number_of_dimensions,
c                vector_for_length_per_dimension,error_code)
c       returns varid.
c     Note: Unlimited dimension must be last
c     Refer to p. 77, netcdf-2 manual

c     NCDOUBLE for REAL*8.  This is independent of the
c       internal representation of data on a particular
c       machine, e.g., 32- or 64-bit based calculations.


c--------------------------
c     Genray version
c--------------------------

c      write(*,*)'before ncvdef2(version)'
c      write(*,*)'version',version
      vid=ncvdef2(ncid,'version',NCCHAR,1,char64id,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,21,
     +          'GENRAY version number',istatus)

c--------------------------
c     Mnemonic for the run
c--------------------------
c      write(*,*)'before ncvdef2(mnemonic)'
      vid=ncvdef2(ncid,'mnemonic',NCCHAR,1,char128id,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23,
     +          'Mnemonic run identifier',istatus)

c-------------------------------
c     Run Descriptive Parameters
c-------------------------------
c      write(*,*)'before ncvdef2(vid)'
      vid=ncvdef2(ncid,'id',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,24,
     +            'Disp relation identifier',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(iabsorp)'
      vid=ncvdef2(ncid,'iabsorp',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,22,
     +            'Absorp calc identifier',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(ieffic)'
      vid=ncvdef2(ncid,'ieffic',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,29,
     +            'Current drive calc identifier',istatus)
      call check_err(istatus)

cSAP080303
c      write(*,*)'before ncvdef2(ion_absorption)'
      vid=ncvdef2(ncid,'ion_absorption',NCCHAR,1,char8id,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,49,
     +'Switch on/off ion absorption at iabsorp=3,9,91,92',istatus)
      call check_err(istatus)


c      write(*,*)'before ncvdef2(refl_loss)'
      vid=ncvdef2(ncid,'refl_loss',NCDOUBLE,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,41,
     +    'fraction of power loss at each reflection',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(iflux)'
      vid=ncvdef2(ncid,'iflux',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,34,
     +            'Flux calc, non-Westerhof-Tokman id',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(ioxm)'
      vid=ncvdef2(ncid,'ioxm',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,38,
     +            'wave mode indicator (1 - om, -1 - xm )',istatus)
      call check_err(istatus)


c      write(*,*)'before ncvdef2(ioxm_n_npar)'
      vid=ncvdef2(ncid,'ioxm_n_npar',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,68,
     +'wave mode indicator: sign before square root to find '//
     +' N(N_parallel)',
     +istatus)
      call check_err(istatus)


c      write(*,*)'before ncvdef2(jwave)'
      vid=ncvdef2(ncid,'jwave',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,37,
     +            'Wave harmonic, for CD efficiency calc',istatus)
      call check_err(istatus)

      if (iabsorp.eq.4) then
      vid=ncvdef2(ncid,'i_im_nperp',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,57,
     +       'iabsorp=4  nperp:1,ImD_full/dD/dnp; 2,Cmplx soln',istatus)
      call check_err(istatus)
      endif

c      write(*,*)'before ncvdef2(i_geom_optic)'
      vid=ncvdef2(ncid,'i_geom_optic',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,34,
     +            'Integrate rays wrt (1)time,(2)dist',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(istart)'
      vid=ncvdef2(ncid,'istart',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,50,
     +     'Ray launch type: 1,eccone; 2,grill, 3,OX in plasma',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(ncone)'
      vid=ncvdef2(ncid,'ncone',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name0',NCCHAR,39,
     +     'Number of rf source ray cones, istart=1',istatus)
      call ncaptc2(ncid,vid,'long_name1',NCCHAR,40,
     +     'Each cone has (nray/ncone) launched rays, istart=1',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(ngrill)'
      vid=ncvdef2(ncid,'ngrill',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name0',NCCHAR,41,
     +     'Number of rf source ray grills, istart=2,3',istatus)
      call ncaptc2(ncid,vid,'long_name0',NCCHAR,47,
     +     ',nray/ngrill rays launched per grill istart=2,3',istatus)
      call ncaptc2(ncid,vid,'long_name0',NCCHAR,55,
     +     ',Each grill has (nray/ngrill) launched rays, istart=2,3',
     +     istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(ionetwo)'
      vid=ncvdef2(ncid,'ionetwo',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +            'if ionetwo=1 then calculate CD',istatus)
      call check_err(istatus)
c--------------------------
c     Ray data
c--------------------------
c      write(*,*)'before ncvdef2(nray)'
      vid=ncvdef2(ncid,'nray',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,14,
     +            'Number of rays',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(nharm)'
      vid=ncvdef2(ncid,'nharm',NCLONG,0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,21,
     +            'First harmonic number',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(nrayelt)'

      vid=ncvdef2(ncid,'nrayelt',NCLONG,1,nraysid,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,35,
     +            'Number of ray elements for each ray',istatus)
      call check_err(istatus)

c      write(*,*)'before ncvdef2(ws)'

      vid=ncvdef2(ncid,'ws',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,29,
     +           'poloidal distance along a ray',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'cms',istatus)

c      write(*,*)'before ncvdef2(seikon)'

      vid=ncvdef2(ncid,'seikon',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,7,
     +           'eikonal',istatus)


c      write(*,*)'before ncvdef2(spsi)'

      vid=ncvdef2(ncid,'spsi',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,45,
     +         'normalized small radius=rho given by indexrho',istatus)

c      write(*,*)'before ncvdef2(wr)'

      vid=ncvdef2(ncid,'wr',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,12,
     +           'major radius',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'cms',istatus)


c      write(*,*)'before ncvdef2(wphi)'

      vid=ncvdef2(ncid,'wphi',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,14,
     +           'toroidal angle',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'radians',istatus)

c      write(*,*)'before ncvdef2(wz)'

      vid=ncvdef2(ncid,'wx',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,15,
     +           'x-cartesian',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'cms',istatus)

      vid=ncvdef2(ncid,'wy',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,15,
     +           'y-cartesian',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'cms',istatus)

      vid=ncvdef2(ncid,'wz',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,15,
     +           'vertical height',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'cms',istatus)

c      write(*,*)'before ncvdef2(w_theta_pol)'

      vid=ncvdef2(ncid,'w_theta_pol',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,14,
     +           'poloidal angle',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,6,
     +           'degree',istatus)

c      write(*,*)'before ncvdef2(wnpar)'

      vid=ncvdef2(ncid,'wnpar',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +           'parallel refractive index',istatus)

c      write(*,*)'before ncvdef2(wnper)'

      vid=ncvdef2(ncid,'wnper',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +           'perpendicular refractive index',istatus)

c      write(*,*)'before ncvdef2(delpwr)'

      vid=ncvdef2(ncid,'delpwr',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,20,
     +           'power in ray channel',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +           'erg/sec',istatus)

c      write(*,*)'before ncvdef2(sdpwr)'


      vid=ncvdef2(ncid,'sdpwr',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,48,
     +      'Ion collisionless absorption coeff (all species)',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           '1/cms',istatus)

      if (iabsorp.eq.3.and.nbulk.gt.1) then ! salphas
      vid=ncvdef2(ncid,'salphas',NCFLOAT,3,ray_dimss,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,49,
     +     'Ion collisionless absorption coeff (each species)',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           '1/cms',istatus)
      endif

      vid=ncvdef2(ncid,'vthermal',NCFLOAT,3,ray_dimsp,istatus)
      call check_err(istatus)
      ! Note: ray_dimsp includes ALL species, while ray_dimss - ions only
      call ncaptc2(ncid,vid,'long_name',NCCHAR,39,
     +     'Thermal speed along ray (each species)',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,4,
     +           'cm/s',istatus)

c      write(*,*)'before ncvdef2(wdnpar)'

      vid=ncvdef2(ncid,'wdnpar',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,6,
     +           'wdnpar',istatus)

c      write(*,*)'before ncvdef2(cwexde)'

c     Added 3rd dimension equal to 2 accomodates complex data.
c      write(*,*)'ncid=',ncid
c      write(*,*)'ray_dims',ray_dims
      vid=ncvdef2(ncid,'cwexde',NCDOUBLE,3,ray_dims,istatus)
      call check_err(istatus)
c      write(*,*)'istatus',istatus
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +           'Complex Ex/E Polarization',istatus)

c      write(*,*)'before ncvdef2(cweyde)'

      vid=ncvdef2(ncid,'cweyde',NCDOUBLE,3,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +           'Complex Ey/E Polarization',istatus)

c      write(*,*)'before ncvdef2(cwezde)'

      vid=ncvdef2(ncid,'cwezde',NCDOUBLE,3,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +           'Complex Ez/E Polarization',istatus)

c      write(*,*)'before ncvdef2(fluxn)'

      vid=ncvdef2(ncid,'fluxn',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23,
     +           'fluxn, Stix norm, |E|=1',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,13,
     +           'ergs/sec/cm^2',istatus)

c      write(*,*)'before ncvdef2(sbtot)'

      vid=ncvdef2(ncid,'sbtot',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23,
     +           'Magnetic field strength',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           'Gauss',istatus)

c      write(*,*)'before ncvdef2(cene)'

      vid=ncvdef2(ncid,'sene',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,17,
     +           'Density along ray',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,14,
     +           'particles/cm^3',istatus)
     
c      write(*,*)'before ncvdef2(ste)'

      vid=ncvdef2(ncid,'ste',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,21,
     +           'Temperature along ray',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,
     +           'keV',istatus)

c      write(*,*)'before ncvdef2(salphac)'

      vid=ncvdef2(ncid,'salphac',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +           'Collisional damping wavenumber',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           '1/cms',istatus)

c      write(*,*)'before ncvdef2(salphal)'

      vid=ncvdef2(ncid,'salphal',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +           'Linear damping wavenumber',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           '1/cms',istatus)

c      write(*,*)'before ncvdef2(sb_r)'

      vid=ncvdef2(ncid,'sb_r',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,28,
     +           'B_r magnetic field component',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           'Gauss',istatus)

c      write(*,*)'before ncvdef2(sb_z)'

      vid=ncvdef2(ncid,'sb_x',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,28,
     +           'B_x magnetic field component',istatus) 
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           'Gauss',istatus)

      vid=ncvdef2(ncid,'sb_y',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,28,
     +           'B_y magnetic field component',istatus) 
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           'Gauss',istatus)

      vid=ncvdef2(ncid,'sb_z',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,28,
     +           'B_z magnetic field component',istatus) 
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           'Gauss',istatus)

c      write(*,*)'before ncvdef2(sb_phi)'

      vid=ncvdef2(ncid,'sb_phi',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +           'B_phi magnetic field component',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5,
     +           'Gauss',istatus)

c      write(*,*)'before ncvdef2(wn_r)'

      vid=ncvdef2(ncid,'wn_r',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +           'N_r refractive index component',istatus)

c      write(*,*)'before ncvdef2(wn_z)'

      vid=ncvdef2(ncid,'wn_x',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +           'N_x refractive index component',istatus) 

      vid=ncvdef2(ncid,'wn_y',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +           'N_y refractive index component',istatus) 

      vid=ncvdef2(ncid,'wn_z',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30,
     +           'N_z refractive index component',istatus) 


      vid=ncvdef2(ncid,'wn_phi',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,32,
     +           'N_phi refractive index component',istatus)
     
c      write(*,*)'before ncvdef2(wgr_r)'

      vid=ncvdef2(ncid,'vgr_r',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,24,
     +           'vgroup_r normalized to c',istatus)
     
c      write(*,*)'before ncvdef2(wgr_z)'

      vid=ncvdef2(ncid,'vgr_x',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,24,
     +           'vgroup_x normalized to c',istatus)
     
      vid=ncvdef2(ncid,'vgr_y',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,24,
     +           'vgroup_y normalized to c',istatus)
     
      vid=ncvdef2(ncid,'vgr_z',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,24,
     +           'vgroup_z normalized to c',istatus)
     
c      write(*,*)'before ncvdef2(wgr_hpi)'

      vid=ncvdef2(ncid,'vgr_phi',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,26,
     +           'vgroup_phi normalized to c',istatus)
     
c      write(*,*)'before ncvdef2(flux_z)'

       
      vid=ncvdef2(ncid,'flux_z',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,22,
     +           'flux_z normalized to c',istatus)
     
c      write(*,*)'before ncvdef2(flux_r)'

       
      vid=ncvdef2(ncid,'flux_r',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,22,
     +           'flux_r normalized to c',istatus)
     
c      write(*,*)'before ncvdef2(flux_phi)'

       
      vid=ncvdef2(ncid,'flux_phi',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,24,
     +           'flux_phi normalized to c',istatus)

c--------------------------------------------------------
c      CD efficiency along rays
c--------------------------------------------------------
      if (ionetwo.eq.1)then
c         write(*,*)'before ncvdef2(w_eff_nc)'
         vid=ncvdef2(ncid,'w_eff_nc',NCDOUBLE,2,ray_dims,istatus)
      call check_err(istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +           'CD efficiency along a ray',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,27,
     +           '(A/cm**2)/(erg/(sec*cm**3))',istatus)
       endif
c--------------------------
c     determine OX conversion data for i_ox=2 case
c--------------------------
      goto 10
      if (i_ox.eq.2) then

         vid=ncvdef2(ncid,'i_ox_conversion',NCLONG,1,nraysid,istatus)
      call check_err(istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,35,
     &            'Option equals 1 after OX conversion',istatus)
         call check_err(istatus)

         vid=ncvdef2(ncid,'transm_ox',NCDOUBLE,1,nraysid,istatus)
      call check_err(istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,57,
     &   'OX transmission coefficient Preinhaelter and Kopecky 1973',
     &              istatus)
         call check_err(istatus)

         vid=ncvdef2(ncid,'cn_par_optimal',NCDOUBLE,1,nraysid,istatus)
      call check_err(istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,38,
     &              'optimal N parallel for OX transmission',
     &              istatus)
         call check_err(istatus)

         vid=ncvdef2(ncid,'cnpar_ox',NCDOUBLE,1,nraysid,istatus)
      call check_err(istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,33,
     &              'N parallel before OX transmission',
     &              istatus)
         call check_err(istatus)

         vid=ncvdef2(ncid,'cn_b_gradpsi',NCDOUBLE,1,nraysid,istatus)
      call check_err(istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,45,
     &              'N along [B^*grad(psi)] before OX transmission',
     &              istatus)
         call check_err(istatus)
         
      endif !i_ox=2
 10   continue

c------------------------------------------------------------------
c     for total power [erg/sec] absorbed at all reflections at all rays
c------------------------------------------------------------------
c      write(*,*)'before ncvdef2(w_tot_pow_absorb_at_refl_nc)'
      vid=ncvdef2(ncid,'w_tot_pow_absorb_at_refl_nc',NCDOUBLE,
     &0,0,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,47,
     +    'Total power absorbed at reflections of all rays',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +           'erg/sec',istatus)
c.................................
cl    1.1.4 end the define-mode and start the data-mode
c     p. 51-2 of manual
c      write(*,*)'ncid',ncid
      call ncendf3(ncid,istatus)
      call check_err(istatus)

c      write(*,*)'end initialization'

      endif               ! End initialize




c.......................................................................
cl    1. Writing data
c

      if ( kopt.eq.0 ) then
c      write(*,*)'wrtnetcdf,  kopt=',kopt
cl    1.2 Write data
c

c     First get variable_id:
c             nf_inq_varid(netcdf_id,variable_name,integer_info)
c     Then write data with nc variable_put:
c             nf_put_var... ()
c
c     function nf_put_var...(ncid,variable_id,index,val)
c     p. 52,, 54-65 of netcdf-3 manual
c
c--------------------------
c      write(*,*)'one_ray_point version',version
c      write(*,*)'ncid',ncid
      call ncvid2(vid,ncid,'version',istatus)
c      write(*,*)'vid',vid 
      ll=length_char(version)
c      write(*,*)'ll',ll
      call ncvptc2(ncid,vid,1,ll,version,ll,istatus)
      call check_err(istatus)

c      write(*,*)'befor vid=ncvid(ncid,mnemonic'
      call ncvid2(vid,ncid,'mnemonic',istatus)
      ll=length_char(mnemonic)
      call ncvptc2(ncid,vid,1,ll,mnemonic,ll,istatus)

c      write(*,*)'befor vid=ncvid(ncid,eqdskin'
      call ncvid2(vid,ncid,'eqdskin',istatus)
      ll=length_char(eqdskin)
      call ncvptc2(ncid,vid,1,ll,eqdskin,ll,istatus)

c      write(*,*)'befor vid=ncvid(ncid,dmas'
      call ncvid2(vid,ncid,'dmas',istatus)
      call ncvpt_doubl2(ncid,vid,1,nbulk,dmas,istatus)

c      write(*,*)'befor vid=ncvid(ncid,charge'
      call ncvid2(vid,ncid,'charge',istatus)
      call ncvpt_doubl2(ncid,vid,1,nbulk,charge,istatus)

c      write(*,*)'befor vid=ncvid(ncid,istart'
      call ncvid2(vid,ncid,'istart',istatus)
      call ncvpt_int2(ncid,vid,0,0,istart,istatus)

c      write(*,*)'befor vid=ncvid(ncid,ncone'
      call ncvid2(vid,ncid,'ncone',istatus)
      call ncvpt_int2(ncid,vid,0,0,ncone,istatus)

c      write(*,*)'befor vid=ncvid(ncid,ngrill'
      call ncvid2(vid,ncid,'ngrill',istatus)
      call ncvpt_int2(ncid,vid,0,0,ngrill,istatus)

c      write(*,*)'befor vid=ncvid(ncid,ionetwo'
      call ncvid2(vid,ncid,'ionetwo',istatus)
      call ncvpt_int2(ncid,vid,0,0,ionetwo,istatus)

c--------------------------
c     Run specs
c--------------------------  
c      write(*,*)'befor vid=ncvid(ncid,id'
      call ncvid2(vid,ncid,'id',istatus)
      call ncvpt_int2(ncid,vid,0,0,id,istatus)

c      write(*,*)'netcdfr3d.f refl_loss',refl_loss
      call ncvid2(vid,ncid,'refl_loss',istatus)
      call ncvpt_doubl2(ncid,vid,0,0,refl_loss,istatus)

c      write(*,*)'befor vid=ncvid(ncid,iabsorp'
      call ncvid2(vid,ncid,'iabsorp',istatus)
      call ncvpt_int2(ncid,vid,0,0,iabsorp,istatus)

c      write(*,*)'befor vid=ncvid(ncid,ieffic'
      call ncvid2(vid,ncid,'ieffic',istatus)
      call ncvpt_int2(ncid,vid,0,0,ieffic,istatus)

cSAP080303
c      write(*,*)'befor vid=ncvid(ncid,ion_absoption'
      call ncvid2(vid,ncid,'ion_absorption',istatus)
      ll=length_char(ion_absorption)
      call ncvptc2(ncid,vid,1,ll,ion_absorption,ll,istatus)

c      write(*,*)'befor vid=ncvid(ncid,iflux'
      call ncvid2(vid,ncid,'iflux',istatus)
      call ncvpt_int2(ncid,vid,0,0,iflux,istatus)

      call ncvid2(vid,ncid,'ioxm',istatus)
      call ncvpt_int2(ncid,vid,0,0,ioxm,istatus)

      call ncvid2(vid,ncid,'ioxm_n_npar',istatus)
      call ncvpt_int2(ncid,vid,0,0,ioxm_n_npar,istatus)

      call ncvid2(vid,ncid,'jwave',istatus)
      call ncvpt_int2(ncid,vid,0,0,jwave,istatus)

      call ncvid2(vid,ncid,'i_geom_optic',istatus)
      call ncvpt_int2(ncid,vid,0,0,i_geom_optic,istatus)

      
      if (iabsorp.eq.4) then
         call ncvid2(vid,ncid,'i_im_nperp',istatus)
         call ncvpt_int2(ncid,vid,0,0,i_im_nperp,istatus)
      endif

      call ncvid2(vid,ncid,'nbulk',istatus)
      call ncvpt_int2(ncid,vid,0,0,nbulk,istatus)

c--------------------------
c     Ray data
c--------------------------
      call ncvid2(vid,ncid,'nray',istatus)
      call ncvpt_int2(ncid,vid,0,0,nrayl0,istatus)

      call ncvid2(vid,ncid,'nharm',istatus)
      call ncvpt_int2(ncid,vid,0,0,nharm,istatus)

      call ncvid2(vid,ncid,'nrayelt',istatus)
      ipoint_1D(1)=1
      call ncvpt_int2(ncid,vid,1,nrayl0,ipoint_1D,istatus)

      point_2D(1,1)=ws_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'ws',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      point_2D(1,1)=seikon_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'seikon',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=spsi_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'spsi',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=wr_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'wr',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314      
      point_2D(1,1)=wphi_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'wphi',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=wx_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'wx',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      point_2D(1,1)=wy_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'wy',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      point_2D(1,1)=wz_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'wz',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=w_theta_pol_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'w_theta_pol',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=wnpar_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'wnpar',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=wnper_nc(is0,nray0)
      
cSAP090903
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'wnper',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=delpwr_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'delpwr',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

cSAP090314
      point_2D(1,1)=sdpwr_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'sdpwr',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      if (iabsorp.eq.3) then ! salphas_nc
      do j=1,nbulka
         point_3D(1,1,j)=salphas_nc(is0,nray0,j) ! in single prec.
      enddo
      do i=2,nbulk
        call pack21_single(point_3D(1,1,i),1,nrelta,1,nrayl,
     &  tem1_single,neltmax,nrayl0)
        starts(1)=start(1)
        starts(2)=start(2)
        starts(3)=i-1  ! nbulk-1 ion species salphas are put in .nc file
        call ncvid2(vid,ncid,'salphas',istatus)
        call ncvpt_real(ncid,vid,starts,ray_counts,tem1_single,istatus)
      enddo
      endif

      do i=1,nbulka
         point_3D(1,1,i)= wvthermal_nc(is0,nray0,i) ! in single prec.
      enddo
      do i=1,nbulk
        call pack21_single(point_3D(1,1,i),1,nrelta,1,nrayl,
     &  tem1_single,neltmax,nrayl0)
        starts(1)=start(1)
        starts(2)=start(2)
        starts(3)=i  ! nbulk species are put in .nc file
        call ncvid2(vid,ncid,'vthermal',istatus)
        call ncvpt_real(ncid,vid,starts,ray_counts,tem1_single,istatus)
      enddo

      point_2D(1,1)=wdnpar_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'wdnpar',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)
c-----------------------------------------------------------------
      ii=0
      j=nray0
      i=is0
      ii=ii+1
      tem1(ii)=0.5d0*(cwexde_nc(i,j)+dconjg(cwexde_nc(i,j)))
      
      ii=ii+1
      tem1(ii)=-cei*0.5d0*(cwexde_nc(i,j)-dconjg(cwexde_nc(i,j)))
      call ncvid2(vid,ncid,'cwexde',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)
c-------------------------------------------------------------------
      ii=0
      j=nray0
      i=is0
      ii=ii+1
      tem1(ii)=0.5d0*(cweyde_nc(i,j)+dconjg(cweyde_nc(i,j)))

      ii=ii+1
      tem1(ii)=-cei*0.5d0*(cweyde_nc(i,j)-dconjg(cweyde_nc(i,j)))
      call ncvid2(vid,ncid,'cweyde',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

c-------------------------------------------------------------
      ii=0
      j=nray0
      i=is0
      ii=ii+1
      tem1(ii)=0.5d0*(cwezde_nc(i,j)+dconjg(cwezde_nc(i,j)))

      ii=ii+1
      tem1(ii)=-cei*0.5d0*(cwezde_nc(i,j)-dconjg(cwezde_nc(i,j)))
      call ncvid2(vid,ncid,'cwezde',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)
c---------------------------------------------------------------

      point_2D(1,1)=fluxn_nc(is0,nray0)
      call pack21(point_2D,1,nrelt,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'fluxn',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      point_2D(1,1)=sbtot_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'sbtot',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      point_2D(1,1)=sene_nc(is0,nray0)
      call pack21(point_2D,1,nrelt,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'sene',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      point_2D(1,1)=ste_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'ste',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      point_2D(1,1)=salphac_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'salphac',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      point_2D(1,1)=salphal_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'salphal',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      point_2D(1,1)=sb_r_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'sb_r',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      point_2D(1,1)=sb_x_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'sb_x',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      point_2D(1,1)=sb_y_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'sb_y',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      point_2D(1,1)=sb_z_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'sb_z',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      point_2D(1,1)=sb_phi_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'sb_phi',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      point_2D(1,1)=wn_r_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'wn_r',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      point_2D(1,1)=wn_x_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'wn_x',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      point_2D(1,1)=wn_y_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'wn_y',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      point_2D(1,1)=wn_z_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'wn_z',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      point_2D(1,1)=wn_phi_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'wn_phi',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      point_2D(1,1)=vgr_r_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'vgr_r',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      point_2D(1,1)=vgr_x_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'vgr_x',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      point_2D(1,1)=vgr_y_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'vgr_y',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      point_2D(1,1)=vgr_z_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'vgr_z',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      point_2D(1,1)=vgr_phi_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'vgr_phi',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      point_2D(1,1)=flux_r_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'flux_r',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      point_2D(1,1)=flux_z_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'flux_z',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)

      point_2D(1,1)=flux_phi_nc(is0,nray0)
      call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
      call ncvid2(vid,ncid,'flux_phi',istatus)
      call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)
c--------------------------------------------------------------
c      CD efficiency along rays
c--------------------------------------------------------
      point_2D(1,1)=w_eff_nc(is0,nray0)
       if (ionetwo.eq.1) then            
          call pack21(point_2D,1,nrelta,1,nrayl,tem1,neltmax,nrayl0)
          call ncvid2(vid,ncid,'w_eff_nc',istatus)
          call ncvpt_doubl2(ncid,vid,start,ray_count,tem1,istatus)
       endif
c--------------------------
c     write OX conversion data for i_ox=2 case
c--------------------------
      goto 20
      if (i_ox.eq.2) then

c         write(*,*)'before vid=ncvid(ncid,i_ox_conversion'
         call ncvid2(vid,ncid,'i_ox_conversion',istatus)
         call ncvpt_int2(ncid,vid,1,nrayl,i_ox_conversion_nc,istatus)

c          write(*,*)'before vid=ncvid(ncid,transm_ox'
         call ncvid2(vid,ncid,'transm_ox',istatus)
         call ncvpt_doubl2(ncid,vid,1,nrayl,transm_ox_nc,istatus)

c         write(*,*)'before vid=ncvid(ncid,cn_par_optimal'
         call ncvid2(vid,ncid,'cn_par_optimal',istatus)
         call ncvpt_doubl2(ncid,vid,1,nrayl,cn_par_optimal_nc,istatus)

c         write(*,*)'before vid=ncvid(ncid,w_tot_pow_absorb_at_refl_nc'
c         call ncvid2(vid,ncid,'cn_par_optimal',istatus)
c         call ncvpt_doubl2(ncid,vid,1,nrayl,cn_par_optimal_nc,istatus)

c         write(*,*)'before vid=ncvid(ncid,cnpar_ox'
         call ncvid2(vid,ncid,'cnpar_ox',istatus)
         call ncvpt_doubl2(ncid,vid,1,nrayl,cnpar_ox_nc,istatus)

c         write(*,*)'before vid=ncvid(ncid,cn_b_gradpsi'
         call ncvid2(vid,ncid,'cn_b_gradpsi',istatus)
         call ncvpt_doubl2(ncid,vid,1,nrayl,cn_b_gradpsi_nc,istatus)

      endif !i_ox=2 
 20   continue
c------------------------------------------------------------------
c     for total power absorbed at all reflections at all rays
c------------------------------------------------------------------
c      write(*,*)'befor vid=ncvid(ncid,w_tot_pow_absorb_at_refl_nc'
      call ncvid2(vid,ncid,'w_tot_pow_absorb_at_refl_nc',istatus)
      call ncvpt_doubl2(ncid,vid,0,0, w_tot_pow_absorb_at_refl_nc,
     +                  istatus)
      
C-----------------------------------------------------------------------
c
cl    4. Close netCDF file
c        Page 53-54 netCDF-2 Manual
c
      call ncclos3(ncid,istatus)
      call check_err(istatus)

      endif                    ! End write netcdf data

      deallocate (tem1,STAT=istat)
      return
      end
c
c


C=YuP=> ADDED: conversion from Netcdf-2 to Netcdf-3 or higher ==========
C These routines/function convert old routines:

      integer function ncvdef2(NCID,VARNAM,VARTYP,NDIMS,VDIMS,istatus)
      ! vid=ncvdef2() is renamed to vid=ncvdef2() in *.f files
      INCLUDE 'netcdf.inc'
      CHARACTER*(*) VARNAM
      INTEGER istatus,vid,NCID, VARTYP,XTYPE, NDIMS,VDIMS(*)
      if (VARTYP.eq.NCFLOAT)  XTYPE=NF_FLOAT    ! 32 BITS
      if (VARTYP.eq.NCDOUBLE) XTYPE=NF_DOUBLE   ! 64 BITS
      if (VARTYP.eq.NCCHAR)   XTYPE=NF_CHAR
      if (VARTYP.eq.NCBYTE)   XTYPE=NF_BYTE
      if (VARTYP.eq.NCSHORT)  XTYPE=NF_SHORT    ! 16 BITS
      if (VARTYP.eq.NCLONG)   XTYPE=NF_INT      ! 32 BITS
      istatus = NF_DEF_VAR(NCID,VARNAM,XTYPE,NDIMS,VDIMS,vid)
      ncvdef2 = vid
      end

      integer function ncdid2(ncid,VARNAM,istatus)
      INCLUDE 'netcdf.inc'
      CHARACTER*(*) VARNAM
      INTEGER istatus,ncid,ndim
c     Get dimension ID from dimension name
c     integer function ncdid(ncid,character(*) lat,integer error_code)
c-YuP:      neltdim=ncdid(ncid,'neltmax',istatus)
      istatus= NF_INQ_DIMID(ncid,VARNAM,ndim) 
      ncdid2 = ndim
      end

      integer function ncddef2(ncid,VARNAM,LEN,istatus)
      INCLUDE 'netcdf.inc'
      CHARACTER*(*) VARNAM
      INTEGER istatus,ncid,LEN,RECID
c     returns dimension id.
      if(LEN.eq.NCUNLIM) then
        ! unlimited dimension for time, dimension name= 'time'
        istatus= NF_DEF_DIM(ncid, VARNAM, NF_UNLIMITED, RECID)
      else
        istatus= NF_DEF_DIM(ncid, VARNAM, LEN, RECID)
      endif
      ncddef2 = RECID
      end

      integer function nccre2(filename,MODE,istatus)
      INCLUDE 'netcdf.inc'
      CHARACTER*(*) filename
      INTEGER istatus,ncid,MODE
cl    create netCDF filename (Entering define mode.)
      if(MODE.eq.NCCLOB)then
        ! over-write existing file
        istatus= NF_CREATE(filename, NF_CLOBBER, ncid) 
      else
        ! Do not over-write existing file
        istatus= NF_CREATE(filename, NF_NOCLOBBER, ncid) 
      endif
      nccre2 = ncid
      end

      subroutine ncvgtc3(NCID, vid, START, COUNTS, TEXT, LEN, istatus)
      INCLUDE 'netcdf.inc'
      CHARACTER*(*) TEXT
      INTEGER istatus,vid,NCID, LEN, START(*), COUNTS(*)
      istatus= NF_GET_VARA_TEXT (NCID, vid, START, COUNTS, TEXT)
      return 
      end

      subroutine ncvgt_int2(NCID, vid, START, COUNTS,  IVALS, istatus)
      INCLUDE 'netcdf.inc'
      INTEGER istatus,vid,NCID, START(*), COUNTS(*)
      INTEGER IVALS(*)
      istatus= NF_GET_VARA_INT(NCID, vid, START, COUNTS, IVALS)
      return 
      end

      subroutine ncvgt_doubl2(NCID, vid, START, COUNTS,  DVALS, istatus)
      INCLUDE 'netcdf.inc'
      INTEGER istatus,vid,NCID, START(*), COUNTS(*)
      REAL*8 DVALS(*)
      istatus= NF_GET_VARA_DOUBLE(NCID, vid, START, COUNTS, DVALS)
      return 
      end
      
      
      subroutine ncaptc2(NCID, vid, NAME, ATTYPE, LEN, TEXT, istatus)
      INCLUDE 'netcdf.inc'
      CHARACTER*(*) NAME
      CHARACTER*(*) TEXT
      INTEGER istatus,vid,NCID, LEN, ATTYPE
      istatus= NF_PUT_ATT_TEXT(NCID, vid, NAME, LEN, TEXT)
      return 
      end
        
      subroutine ncvptc2(NCID, vid, START, COUNTS, TEXT, LEN, istatus)
      INCLUDE 'netcdf.inc'
      CHARACTER*(*) TEXT
      INTEGER istatus,vid,NCID, LEN, START(*), COUNTS(*)
      istatus= NF_PUT_VARA_TEXT(NCID, vid, START, COUNTS, TEXT)
      return 
      end

      subroutine ncvpt_doubl2(NCID, vid, START, COUNTS,  DVALS, istatus)
      INCLUDE 'netcdf.inc'
      INTEGER istatus,vid,NCID, START(*), COUNTS(*)
      REAL*8 DVALS(*)
      istatus= NF_PUT_VARA_DOUBLE(NCID, vid, START, COUNTS, DVALS)
      return 
      end

      subroutine ncvpt_real(NCID, vid, START, COUNTS,  DVALS, istatus)
      INCLUDE 'netcdf.inc'
      INTEGER istatus,vid,NCID, START(*), COUNTS(*)
      REAL*4 DVALS(*)
      istatus= NF_PUT_VARA_REAL(NCID, vid, START, COUNTS, DVALS)
      return 
      end

      subroutine ncvpt_int2(NCID, vid, START, COUNTS,  IVALS, istatus)
      INCLUDE 'netcdf.inc'
      INTEGER istatus,vid,NCID, START(*), COUNTS(*)
      INTEGER IVALS(*)
      istatus= NF_PUT_VARA_INT(NCID, vid, START, COUNTS, IVALS)
      return 
      end
      
      subroutine ncopn2(filename,MODE,ncid,istatus)
      INCLUDE 'netcdf.inc'
      INTEGER istatus,MODE,ncid
      CHARACTER*(*) filename
      !integer function ncopn(filename,write,error_code) ! NetCDF-2
      !ncid=ncopn(netcdfnml,NCWRITE,istatus) !Open existing netCDF file
      if(MODE.eq.0) then
         istatus= NF_OPEN(filename, NF_NOWRITE, ncid) 
         else
         istatus= NF_OPEN(filename, NF_WRITE, ncid) 
      endif
      if (istatus .NE. NF_NOERR) then         
         write(*,*)'   ***   Problem opening .nc data file   ***'
         Stop
      endif  
      return                                
      end

      subroutine ncclos3(ncid,istatus)
      INCLUDE 'netcdf.inc'
      INTEGER istatus,ncid
      istatus= NF_CLOSE(ncid)
      return 
      end

      subroutine ncvid2(vid,ncid,NAME,istatus)
      ! vid=ncvid(ncid,NAME,istatus) ! NetCDF-2
      INCLUDE 'netcdf.inc'
      CHARACTER*(*) NAME
      INTEGER istatus,ncid,vid
      istatus= NF_INQ_VARID(ncid,NAME,vid)
      return
      end

      subroutine ncdinq3(ncid,DIMID,NAME,LEN,istatus) 
c     Query netcdf file for dimensions:
      INCLUDE 'netcdf.inc'
      CHARACTER*(*) NAME
      INTEGER istatus,ncid,DIMID,LEN
      istatus= NF_INQ_DIM(ncid,DIMID,NAME,LEN)
      return
      end

      subroutine ncredf3(ncid,istatus)
c     Put Open NetCDF file into Define Mode
      INCLUDE 'netcdf.inc'
      INTEGER istatus,ncid
      istatus= NF_REDEF(ncid) ! put into define mode
      return
      end

      subroutine ncendf3(ncid,istatus)
c     end the define-mode and start the data-mode
      INCLUDE 'netcdf.inc'
      INTEGER istatus,ncid
      istatus= NF_ENDDEF(ncid) 
      return
      end

         