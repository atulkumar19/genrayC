
c************************************************************
      character version*64
      parameter(version="Genray-c_160826.1")
c************************************************************

c************************************************************
c number of rays in common/cone/ EC-cone
       integer nraymax
       parameter(nraymax=200)
c nraymax must be greater than or equal to nray 
c (nray is calculated)

c maximum number of EC cones
       integer nconea
       parameter(nconea=20)
c maximal number of radial bins at the launching disk
       integer n_mesh_disk_radial_bin_a
       parameter (n_mesh_disk_radial_bin_a=5)
       
c************************************************************
c for reading dendsk file
      integer nxdena,nydena, nx4a,ny4a,nxya
      parameter (nxdena=256,nydena=256)
c for the 2D-spline of dengrid:
      parameter (nx4a=nxdena+4,ny4a=nydena+4,nxya=ny4a)
       
c************************************************************
c     for common/five/
      integer nxeqda,nyeqda,nreqda,nzeqda,nlimit,nr4a,nz4a,nrya,nlim4

      parameter (nreqda=530,nzeqda=530,nlimit=100)
      parameter (nxeqda=2*nreqda, nyeqda=2*nreqda)
      
      parameter (nr4a=nreqda+4,nz4a=nzeqda+4,nrya=nz4a)
      parameter (nlim4=nlimit+4)
c It should be nrya=max(nreqda,nzeqda)+4
c************************************************************
c     for common/fourb/
      integer nves
      parameter (nves=62)
c************************************************************
c     for common gr.i
      INTEGER NL,NP,nteta,npsi,nteta1
      real*8 epspsi
      parameter (nteta=100,npsi=50,epspsi=0.0001d0,nteta1=nteta+1)
      PARAMETER (NL=5,NP=nteta)
c npsi is a number of contours
c nteta  is a number of the points along the each contours
c                       (in the poloidal direction)
c epspsi is the accuracy for the determination of the contours
c         poins coordinates	(zpsi,rpsi)
c NL is a number of contours for the ploting of the tokamak
c              cross section
c NP is a number of points along each contours in poloidal direction
c     for the plotting of the tokamak cross section
c************************************************************
c     for common/grill/
      integer ngrilla,nnkprmax,nraymaxl,nthinmax
      parameter (ngrilla=10,nnkprmax=60,nraymaxl=1000,nthinmax=64)
c     for N_toroidal N_poloidal initial condition i_n_poloidal=4 
      integer nnktormax,nnkpolmax
      parameter (nnktormax=60,nnkpolmax=60)

c------------------------------------
c nnkprmax =max{i=1,ngril1}nnkpar(i)
c nnktormax=max{i=1,ngril1}nnktor(i)
c nnkpolmax=max{i=1,ngril1}nnkpol(i)
c------------------------------------
c ngrilla    is a maximal number of the N_parallel spectra
c ngrill     is a number of the N_parallel spectra
c nnkprmax=max{i=1,ngrill}nnkpar(i)
c nthinmax=max{i=1,ngrill}nthin(i)
c nraymaxl is the max number of rays from the grill
c************************************************************
c     for common/ions/
      integer nbulka,nbulkma
      parameter (nbulka=4)
      parameter (nbulkma=nbulka-1)
c ncomp=nbulk -the number of plasma species.ge.1 & .le.nbulka
c************************************************************
c     for common/onetwo/
      integer NRA
      PARAMETER (NRA=401)
c NRA is the maximal number of bin boundaries in the small radius direction 
c for the calculation of the power and current drive radial profiles.
c Power and current is tabulated at (NR-1) bin centers.
c************************************************************
c     for common/rho/
      integer npsi4,nteta14,nzy
      parameter (npsi4=npsi+4,nteta14=nteta1+4)
c npsi  and nteta1 are in common/gr/
c nzy=max(npsi,nteta1)+4
      parameter(nzy=nteta1+4)
c************************************************************
c     for common/six/
      integer ndensa,ndens4a
      parameter (ndensa=201)
      parameter (ndens4a=ndensa+4)
c ndensa is the max number of points in arrays with the 
c plasma density, temperature, zeff. tpop, and vflow
c************************************************************
c     for common/write/
      integer nrelta
      integer n_relt_harma,n_relt_harm1a,n_relt_harm2a
      parameter (nrelta=32000) 
     
c nrelta is maximum value for nrelt
c nrelt is the max number of the ray elements along every ray
c nraymax is the max value for nray (the number of the rays)

c nfreqa is the max value for nfreq 
c (nfreq is the number of the emission frequencies)
c n_relt_harma  is the number of EC harmonics used in anti-hermitian
c               dielectric tensor calculations  n=<n_relt_harm 
      parameter (n_relt_harma=500)
      parameter (n_relt_harm1a=-5)
      parameter (n_relt_harm2a=500)
c************************************************************
c     for common /output/
      integer n_plot_dispa,n_plot_disp_colda,
     &m_r_nperp_a,m_i_nperp_a,n_contour_plot_disp_a
      parameter (n_plot_dispa=10)
      parameter (n_plot_disp_colda=10)
      parameter(m_r_nperp_a=20,m_i_nperp_a=20,n_contour_plot_disp_a=20)
c*************************************************************
c-----for hot plasma roots n_hot_roots_a max number of hot plasma root
      integer n_hot_roots_a
      parameter (n_hot_roots_a=4)
c******************************************************************
c     for small output step near EC resonance points for power calculation 
c     It will be used for data at namelist /numercl/ and common /one/
      integer n_power_switch_resonance_a
      parameter (n_power_switch_resonance_a=3)
c******************************************************************
c     for toray EC launch
      integer gzonemax
      parameter (gzonemax=20)
c********************************************************************
c     to read wall and limiter coordinates (r,z)
c     for writencdf.i:
      integer n_wall_a,n_limiter_a,max_limiters_a
c     n_wall_a,n_limiter_a are maximal numbes of wall and limiter points
c     max_limiters_a is a maximal number of limiters
      parameter (n_wall_a=200)
      parameter (n_limiter_a=200)
      parameter (max_limiters_a=1)

cYuP  Max. number of coils (current loops) for calculating magnetic field
c     (model_b=2)
      integer ncoilsa
      parameter (ncoilsa=10) ! 

c     maximal number of wall points having the given poloidal angle theta
      integer n_rho_wall_a
      parameter (n_rho_wall_a=10)
c***********************************************************************
c     to read normalized exponential density falls at poloidal mesh
c     for edge_prof_nml,i
c     n_pol_edge_dens_a is a maximal number of mesh points
      integer  n_pol_edge_dens_a
      parameter (n_pol_edge_dens_a=100)
c********************************************************************
c     to create fine meshes of wall coordinates with additional points
      integer n_wall_add_a
      parameter (n_wall_add_a=15000) 
c*********************************************************************
c     for (R,Z) meshes rr_add, zz_add at the poloidal plane
c     for the density fall near the wall
      integer
     &nxeqd_add_a, !max number of x mesh points at the poloidal plane
     &nyeqd_add_a, !max number of y mesh points at the poloidal plane
     &nreqd_add_a, !max number of r mesh points at the poloidal plane
     &nzeqd_add_a  !max number of z mesh points at the poloidal plane

      parameter (nreqd_add_a=1000) 
      parameter (nzeqd_add_a=1000)
      parameter (nxeqd_add_a= 2*nreqd_add_a) 
      parameter (nyeqd_add_a= 2*nreqd_add_a) 
