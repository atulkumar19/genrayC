
c the arrays for the text ray data for FP (CQL3D) code
c      nrelta is maximum value of nrelt
c      in all array nfreqa -> nfreq
       double complex cwexde,cweyde,cwezde
       double complex w_ceps
       real*8
     1              ws,seikon,spsi,
     2              wx,wy,wz, 
     2              wr,wphi,
     3              wnpar,wnper,delpwr,sdpwr,
     4              wdnpar,fluxn,sbtot,
     4	            sb_x,sb_y,sb_z, 
     4              sb_r,sb_phi,
     5              sene,ste,salphac,
     5              salphal,
     5              wye,wyi,wyi2,
     5              wxi,wxi2,
     6              xarr,yarr,
     7              rez,
     9              eff,wmtor,
     &              wn_x,wn_y,wn_z, 
     +              wn_r,wn_phi,
     &              wvgr_x,wvgr_y,wvgr_z,
     &              wvgr_r,wvgr_phi, 
cSAP080906
     &              wtheta_pol,
     &              salphas, 
c YuP[2016] Thermal velocity along ray trajectory, for each species:
     +              wvthermal,

c------------------------------------------------------
c    the data for the emission calculations           
     +              wcnz_em,wcnr_em,wcm_em,
     +              wz_em,wr_em,wphi_em,  
     +              wal_emis,wj_emis,            
     +              wnray,wsn,     
     +              win_sn,win_0,
c     +              wi_0,
     &wi_0sn,
     +              wtemp_em,wtemp_rad_em,
c     +              wtemp_rad_fr_wall,
c     +              wtemp_rad_fr,
c     +              wtemp_pl_fr,
c     +              wr0_em,wz0_em,wrho0_em,
     +              wtaun_em,
c     &wfreq,
c     +              wtau_em,
c     +              wi_0t,
c     +              wr_2nd_harm,
c     +              wtemp_2nd_harm,
     +              freqncy0,
c     +              wn_perp_ioxm_p,wn_perp_ioxm_m,
c     +              wn_perp_ioxm_n_npar_p,wn_perp_ioxm_n_npar_m,
c     +              wye_0,wxe_0,
c-------------------------------------------------------
c    for emission at O-X ebw jump case
c-------------------for O mode ray part 
     +              wsn_o,wz_o,wr_o,wphi_o,
     +              wcnz_em_o,wcnr_em_o,wcm_em_o,
c-------------------for X_EBW mode ray part 
     +              wsn_x,wz_x,wr_x,wphi_x,
     +              wcnz_em_x,wcnr_em_x,wcm_em_x,
c-------------------------------------------------------
c   the data for the determination of the N_par boundaries
     9		    wnrho,gnpar,wxe,wp,
     9		    wnparplb,wnparmnb,
c-------------------------------------------------------
     1              phiold,
     8              xold,yold,zold,rold,rhoold,cld,cvac,
c--------------------------------------------------------
c     These data are for OX transmission coefficient       
     &              transm_ox,cn_par_optimal,cnpar_ox,cn_b_gradpsi,
c-------------------------------------------------------
c     These data are to plot the resonance ellipse boundary
c     for several harmonics    
     &              wp_perpmax_dmvt,wp_parmin_dmvt,wp_parmax_dmvt,
     &              wp_perpmax_dmc,wp_parmin_dmc,wp_parmax_dmc,
cSm060303
     &              wp_0_dmvt,wdel_p_dmvt,
c-------------------------------------------------------------------
c     data for power absorprion at reflections
     &w_tot_pow_absorb_at_refl,
cSAP080902
c-----data to plot delta power_e and delta current along the ray
     &             delpow_e_ar,delcur_par_ar 

      integer       w_ires  

      integer       nrayelt,
     +              nrayelt_emis,ifreq0,i_ox_conversion,
     +              ifreq_write 

      real*8, pointer ::
c-----data for the emission calculations     
     & wi_0(:,:),                    !(nraya,nfreqa)
     & wtemp_rad_fr_wall(:),         !(nfreqa)
     & wtemp_rad_fr(:),              !(nfreqa)
     & wtemp_pl_fr(:),               !(nfreqa)
     & wr0_em(:),                    !(nfreqa)
     & wz0_em(:),                    !(nfreqa)
     & wrho0_em(:),                  !(nfreqa)
     & wfreq(:),                     !(nfreqa)
     & wtau_em(:,:),                 !(nraya,nfreqa)
     & wi_0t(:,:),                   !(nraya,nfreqa)
     & wr_2nd_harm(:),               !(nfreqa)
     & wtemp_2nd_harm(:),            !(nfreqa)
c-------------------to plot cold N_perp in initial ray point for test
     & wn_perp_ioxm_p(:),            !(nfreqa)
     &wn_perp_ioxm_m(:),             !(nfreqa)
     &wn_perp_ioxm_n_npar_p(:),      !(nfreqa)
     &wn_perp_ioxm_n_npar_m(:),      !(nfreqa)
     &wye_0(:),                      !(nfreqa)
     &wxe_0(:)                       !(nfreqa)

      common/write/
     &cwexde(nrelta),
     &cweyde(nrelta),
     &cwezde(nrelta),
     &ws(nrelta),
     &seikon(nrelta),
     &spsi(nrelta),
     &wx(nrelta),wy(nrelta),wz(nrelta),
     &wr(nrelta),wphi(nrelta),
     3              wnpar(nrelta),wnper(nrelta),
     3              delpwr(nrelta),sdpwr(nrelta),
     4              wdnpar(nrelta),fluxn(nrelta),sbtot(nrelta),
     4	            sb_x(nrelta),sb_y(nrelta),sb_z(nrelta),
     +              sb_r(nrelta),sb_phi(nrelta),
     5              sene(nrelta),ste(nrelta),salphac(nrelta),
     5              salphal(nrelta),
     5              wye(nrelta),wyi(nrelta),wyi2(nrelta),
     5              wxi(nrelta),wxi2(nrelta),
     6              xarr(nrelta),yarr(nrelta),
     7              rez(nrelta),
     9              eff(nrelta),wmtor(nrelta),
     &              wn_x(nrelta),wn_y(nrelta),wn_z(nrelta),
     +              wn_r(nrelta),wn_phi(nrelta),
     &              wvgr_x(nrelta),wvgr_y(nrelta),wvgr_z(nrelta),
     +              wvgr_r(nrelta),
     &              wvgr_phi(nrelta),
     &              wtheta_pol(nrelta),
c   for eps checking
     &              w_ceps(3,3,nrelta), 
c    tokman flux
     &              salphas(nrelta,nbulka), 
c YuP[2016] Thermal velocity along ray trajectory, for each species:
     +              wvthermal(nrelta,nbulka),

c------------------------------------------------------
c    the data for the emission calculations           
     +              wcnz_em(nrelta),wcnr_em(nrelta),wcm_em(nrelta),
     +              wz_em(nrelta),wr_em(nrelta),wphi_em(nrelta),  
     +              wal_emis(nrelta),wj_emis(nrelta),            
     +              wnray(nrelta),wsn(nrelta),     
     +              win_sn(nrelta),win_0(nrelta),
     &wi_0,         ! pointer (nraya,nfreqa),
     +              wi_0sn(nrelta),
     +              wtemp_em(nrelta),wtemp_rad_em(nrelta),
     &wtemp_rad_fr_wall, !pointer (nfreqa),
     &wtemp_rad_fr,      !pointer (nfreqa),
     &wtemp_pl_fr,       !pointer (nfreqa),
     &wr0_em,            !pointer (nfreqa),
     &wz0_em,            !pointer (nfreqa),
     &wrho0_em,          !pointer (nfreqa),
     +              wtaun_em(nrelta),
     &wfreq,             !pointer (nfreqa),
     &wtau_em,           !pointer (nraya,nfreqa),
     &wi_0t,             !pointer (nraya,nfreqa),
     &wr_2nd_harm,       !pointer (nfreqa),
     &wtemp_2nd_harm,    !poiunter (nfreqa),
     +              freqncy0,
c-------------------to plot cold N_perp in initial ray point for test
     &wn_perp_ioxm_p,        !pointer (nfreqa),
     &wn_perp_ioxm_m,        !pointer (nfreqa),
     &wn_perp_ioxm_n_npar_p, !pointer (nfreqa),
     &wn_perp_ioxm_n_npar_m, !pointer (nfreqa),
     &wye_0,                 !pointer (nfreqa),
     &wxe_0,                 !pointer (nfreqa),
c-------------------------------------------------------
c    for emission at O-X ebw jump case
c-------------------for O mode ray part 
     +              wsn_o(nrelta),
     +              wz_o(nrelta),wr_o(nrelta),wphi_o(nrelta),
     +              wcnz_em_o(nrelta),wcnr_em_o(nrelta),
     +              wcm_em_o(nrelta),
c-------------------for X_EBW mode ray part 
     +              wsn_x(nrelta),
     +              wz_x(nrelta),wr_x(nrelta),wphi_x(nrelta),
     +              wcnz_em_x(nrelta),wcnr_em_x(nrelta),
     +              wcm_em_x(nrelta),
c-------------------------------------------------------
c   the data for the determination of the N_par boundaries
     9		    wnrho(nrelta),gnpar(nrelta),wxe(nrelta),wp(nrelta),
     9		    wnparplb(nrelta),wnparmnb(nrelta),
c-------------------------------------------------------
     1              phiold,
     8              xold,yold,zold,rold,rhoold,cld,cvac,nrayelt,
     +              nrayelt_emis,ifreq0,
c--------------------------------------------------------
c     These data are for OX trasmition coefficient       
     &              transm_ox,cn_par_optimal,cnpar_ox,cn_b_gradpsi,
     &              i_ox_conversion,
c-------------------------------------------------------
c     These data are to plot the resonance ellipse boundary
c     for the several harmonics    
     &              wp_perpmax_dmvt(nrelta,n_relt_harm1a:n_relt_harm2a),
     &              wp_parmin_dmvt(nrelta,n_relt_harm1a:n_relt_harm2a),
     &              wp_parmax_dmvt(nrelta,n_relt_harm1a:n_relt_harm2a),
     &              wp_perpmax_dmc(nrelta,n_relt_harm1a:n_relt_harm2a),
     &              wp_parmin_dmc(nrelta,n_relt_harm1a:n_relt_harm2a),
     &              wp_parmax_dmc(nrelta,n_relt_harm1a:n_relt_harm2a),
     &              w_ires(nrelta,n_relt_harm1a:n_relt_harma),
cSm060303
     &              wp_0_dmvt(nrelta,n_relt_harm1a:n_relt_harm2a),
     &              wdel_p_dmvt(nrelta,n_relt_harm1a:n_relt_harm2a),
c-------------------------------------------------------------------
c     data for power absorprion at reflections 
     &              w_tot_pow_absorb_at_refl,
c-----data to plot delta power_e and delta current along the ray
     &             delpow_e_ar(nrelta),delcur_par_ar(nrelta), 
c-------------------to plot data in initial ray point 
     &              ifreq_write


c the data for mnemonic.txt file
c nrelta is the max number of the ray elements along any ray
c---------------------------------------------------------
c     ws poloidal ray distance (cm)
c----------------------------------------------------------
c   wal_emis is the absorption(1/cm)           along the ray
c   wj_emis  is the emmisivity (egr*sec/cm**3) along the ray
c   wcnz_em,wcnr_em,wcm_em   the refractive index coordinates along the ray
c   wnray_em           the ray refractive index N_ray along the ray  
c   win_sn_em          emission I_n at the detector side of nth bin at s=s_n
c   win_0_em           emission at the plasma boundary s=s_1 from one nth bin s_n
c   wi_0_em            emission I_0 at the plasma boundary from the ray at s=s(1)=0
c   wsn_em total ray distance (cm)
c   wtemp_em           temperature along the ray
c   wtemp_rad_em       radiation temperature along the ray from the ratio
c                      wtepm_rad_em=(2.d0*pi)**3*(clight)**2/(omega*cnray)**2
c                     *j_emis/al_emis/1.6022d-9
c
c            
c   nrayelt_emis       the number of points along the ray with the additional
c                      emission points
c   wi_0sn(n)          sum{k=1,n}[in_0(k)]! emission at the plasma boundary
c                      from the part 0<s<sn(n) of the ray
c  wtaun_em(n)         tau_n from n-th bin (for ploting only)
c  wfreq(nfreqa)       the array for the frequencies.It is used for the
c                      emission calculations
c  ifreq0              is the number of the frequency in wfreq(ifreq) more
c                      closed to the second elecron gyro-frequency
c                      at plasma center
c                      wfreq(ifreq0) gives min{dabs(wfreq(ifreq)-2*freqncy0)}
c
c  wtemp_rad_fr()       radiation temperature from multi-frequence case
c                      =2.d0*pi*clight**2/1.6022d-9*wi_0_em()/omega**2       
c
c  wtau_em(nraya,nfreqa) is the tau(optical length) from one ray pass
c
c  wi_0t               is the flux wi_0 divided by the coefficient with
c                      wallr wall reflection coefficient
c  wtemp_rad_fr_wall   radiation temperature from malti freq. case
c                      with wallr coefficient
c  wtemp_pl_fr(nfreqa) the plasma temperature at the point along the ray
c                      where In_0/ds has the maximal value 
c
c  wr0_em(nfreqa)      the space coordinates asnd the small radius 
c  wz0_em(nfreqa)      of the ray point where In_0/ds has the maximal value 
c  wrho0_em(nfreqa)
c
c  wr_2nd_harm(nfreqa) the major radius of EC 2nd harmonic resonace points
c                       at Z=0  
c  wtemp_2nd_harm(nfreqa) the bulk plasma temperature at EC 2nd harmonic
c                         points T(z=0,wr_2nd_harm,phi=0)
c
c  w_ceps(3,3,nrelta)   complex dielectric tensor along the ray, for checking 
c
c  wn_z(nrelta),wn_r(nrelta),wn_phi(nrelta) are the refructive index
c                                           coordinates N_r,N_z,N_phi
c
c  wvgr_z(nrelta),wvgr_r(nrelta),wvgr_phi(nrelta) are the group velocity
c  components normalized to c

c YuP[2016] Thermal velocity along ray trajectory, for each species:
c                   wvthermal(nrelta,nbulka)

c  For iabsorp=3, individual ion absorption damping coeffs.
c                 [Could be extended to some additional iabsorp]:
c                  salphas(nrelta,nbulka)
c-----------------------------------------------------------
c          transm_ox OX transmission coefficient
c          i_ox_conversion=1 was the jump in the radial direction
c                         =0 was not OX conversion
c----------------------------------------------------------------
c   cn_par_optimal=dsqrt(Y_abs/(Y_abs+1)) is the optimal N parallel
c                  for OX conversion. It is used for mnemonic.nc file
c-----------------------------------------------------------------
c   cnpar_ox is N parallel before OX conversion procedure.
c            It is used to prepare the data for mnemonic.nc file
c-----------------------------------------------------------------
c   cn_b_gradpsi=N^*[b^*gradpsi^]/|[b^*gradpsi^]| is the refractive index
c                component before OX conversion.
c                It is perpendicular to the magnetic field b^ and
c                to the gradiend(psi)
c----------------------------------------------------------------
c     These data are to plot the resonance ellipse boundary
c     for the several harmonics    
c
c wp_perpmax_dmvt(nrelta,n_relt_harm1a:n_relt_harm2a), maximal value of
c                                                     p_perp/ (m_e*V_thermal)
c                                                     at resonace allopse
c
c wp_parmin_dmvt(nrelta,n_relt_harm1a:n_relt_harm2a)  minimal and the maximal 
c wp_parmax_dmvt(nrelta,n_relt_harm1a:n_relt_harm2a) p_parallel/(m_e*V_thermal)
c                                        at resonance ellipse      
c
c w_ires(nrelta,n_relt_harm1a:n_relt_harm2a) =0 no resonance
c                              1 ellipse
c                              2 parabola
c                              3 hyperbole
c ellipse center for ires=1 case
c wp_0_dmvt(nrelata,n_relt_harm1a:n_relt_harm2a),
c shift from the ellipse center for  ires=1 case
c wdel_p_dmvt(is,n)nrelata,n_relt_harm2a:n_relt_harm2a)
c
c-------------------to plot cold N_perp in initial ray point for test
c                   wn_perp_ioxm_p(nfreqa) is nperp at ioxm=+1
c                   wn_perp_ioxm_m(nfreqa) is nperp at ioxm=-1
c                   wn_perp_ioxm_n_npar_p(freqa) is nperp at ioxm_n_npar=+1
c                   wn_perp_ioxm_n_npar_m(nfreqa) is nperp at ioxm_n_npar=-1
c                   wye_0(nfreqa) ye at initial poin
c                   wxe_0(nfreqa) xe at initial point
c                   ifreq_write
c------------------------------------------------------------------------
c w_tot_pow_absorb_at_refl total power [MWatt] absorbed at all reflections
c                        at all rays
c------------------------------------------------------------------------
c-----data to plot delta power_e and delta current along the ray
c                  delpow_e_ar(nrelta),delcur_par_ar(nrelta)
c---------------------------------------------------------------------------
c     wtheta_pol(nrelta) poloidal angle [degree] -180< wtheta_pol<+180
c---------------------------------------------------------------------------
