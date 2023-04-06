c      all input namelists

      include 'name_genr.i'               !/namelist /genr/
     
      include 'name_tokamak.i'            !/namelist/ tokamak/     

      namelist /wave/ frqncy,ioxm,ireflm,jwave,istart,delpwrmn,ibw,
     *i_vgr_ini,poldist_mx,ioxm_n_npar,
     &cnperp_plot_min,cnperp_plot_max,n_nperp_plot,
     &cN_perp_root_max,n_points_root,
     &i_look_roots,k_hot_root,
     &i_rho_find_hot_nperp_roots,
     &rho_step_find_hot_nperp_roots,rho_min_find_hot_nperp_roots,
     &no_reflection, rho_reflect

      namelist /dispers/ ib,id,iherm,iabsorp,iswitch,del_y,jy_d,
     *idswitch,iabswitch,n_relt_harm,n_relt_intgr,iflux,
     &i_im_nperp,i_geom_optic,ray_direction,errabs0,errrel0,navg,
     &diff_err,relres,iabsorp_collisional,coll_mult,refl_loss,
     &n_relt_harm1,i_salphal,ion_absorption,
     &iabsorp_ql,  rho_larm_max

      namelist /numercl/ irkmeth,ndim1,isolv,idif,nrelt,
     * prmt1,prmt2,prmt3,prmt4,prmt6,prmt9,icorrect,
     * maxsteps_rk,i_output,
     & i_uh_switch,uh_switch,prmt6_uh_switch,
     &toll_hamilt,  
     + dL_step, dN_step, der_r, der_n, der_f,
     &i_power_switch_resonance, 
     &prmt6_power_switch_resonance,
     &n_power_switch_resonance,
     &y_power_switch_resonance,
     &del_y_power_switch_resonance,
     &i_resonance_curve_integration_method,epsi

      namelist /output/ iwcntr,iwopen,iwj,itools,i_plot_b,
     &n_plot_disp,r_plot_disp,id_plot_disp,
     &s_poloid_plot_disp,point_plot_disp,
     &i_plot_disp_cold,
     &n_plot_disp_cold,s_poloid_plot_disp_cold,r_plot_disp_cold,
     &point_plot_disp_cold,     
     &number_map_points_real_nperp,number_map_points_image_nperp,
     &ratio_min_r_nperp,ratio_max_r_nperp,
     &ratio_min_i_nperp,ratio_max_i_nperp,
     &n_contour_plot_disp,
     &r_freq,z_freq,alpha_freq,beta_freq,dist_freq,max_plot_freq,
     &nsteps_freq,n_ec_harmonics_freq,npar_freq,
     + i_save_disp,
     + ixscan_save_disp,
     + inpar_save_disp,
     + inper_save_disp,
     + xmin_save_disp,
     + xmax_save_disp, 
     + y_save_disp, 
     + z_save_disp, 
     + Npar_mn_save_disp,
     + Npar_mx_save_disp,
     + Nper_mn_save_disp,
     + Nper_mx_save_disp
     

      namelist /plasma/ nbulk,izeff,idens,model_rho_dens,
     +  temp_scale,den_scale,
     +  elx0,ely0,elz0, elax,elay,elaz,  ! YuP added !
     +  dens0rr,dens0es,dens0ub,eltheta,sintt,costt, ! YuP added !
     +  Rm0rr,Rm0es,rtau,r_ub_edge,      ! YuP added !
     +  ndens,
     &  nonuniform_profile_mesh,
     +  dendsk  ! YuP added !

      namelist /species/ charge,dmas

      namelist /denprof/ dense0,denseb,rn1de,rn2de,
     + zbegin_den_dropoff, zlength_den_dropoff ! YuP[April,2014] added   

      namelist /tpopprof/ tp0,tpb,rn1tp,rn2tp

      namelist /vflprof/vfl0,vflb,rn1vfl,rn2vfl

      namelist /tprof/ ate0,ateb,rn1te,rn2te

      namelist /zprof/ zeff0,zeffb,rn1zeff,rn2zeff

c----------------------------------------------------------
c     namelists for plasma profiles at uniform radial mesh:
c----------------------------------------------------------
      include 'name_uniform_mesh_profiles.i'     
c      namelist /dentab/    prof  
c      namelist /temtab/    prof  
c      namelist /tpoptab/   prof  
c      namelist /vflowtab/  prof  
c      namelist /zeftab/    zeff1
c---------------------------------------------------------
c     namelists for plasma profiles at non-uniform radial mesh
c     written by lines
c------------------------------------------------------------
      include 'name_non_uniform_mesh_profiles_line.i'  
c      namelist /dentab_nonuniform_line/   nj_tab,prof_2d,radii_2d
c      namelist /temtab_nonuniform_line/   nj_tab,prof_2d,radii_2d
c      namelist /tpoptab_nonuniform_line/  nj_tab,prof_2d,radii_2d
c      namelist /vflowtab_nonuniform_line / nj_tab,prof_2d,radii_2d
c      namelist /zeftab_nonuniform_line /   nj_tab,prof_2d,radii_2d
c-------------------------------------------------------------------

      namelist /read_diskf/ i_diskf,
     . netcdfnm,
     . rtem0,
     . rtail,r1t,r2t,ttail,rhot,r1h,r2h,thotpar,thotper,
     . hotexp,hotmnpar,hotmxpar,hotmnper,hotmxper,
     . rbeam,r1b,r2b,tbeam,ebeam,thbeam,
     . jx,iym,lrz,ngen,
CENM 31Aug05 Added (optional) parameters at the end if dispers namelist
C    to be used in the relativistic dispersion relation in abhay_disp.f
     & rvtovte1,rvtovte2,rtemp1,rtemp2,rtemp3

      include 'name_grill.i'              !namelist /grill/
     
      namelist /ox/ i_ox,
     &theta_bot,theta_top,i_ox_poloidal_max,eps_antenna,
     &eps_xe, ox_step_dir


      include 'name_eccone.i'               !namelist /eccone/
     
      include 'name_edge_prof_nml.i'       !namelist /edge_prof_nml/
