


c        ********************** prep3d***********************
c        *                      -----                       *
c        *  prep3d -subroutine to prepare the parameters    *
c        *  for	 output files for FP code  (e.g., CQL3D)    *
c        ****************************************************
c         input parameters: iray -number of the ray from antenna  
c                                 is in common/cone/ 
c output parameter: iraystop (if power in the ray channel 
c   delpwr(is).lt. delpwrmn*delpwr(1), then iraystop=1). 
c Also, stop ray if fluxn(is).lt.0.0, iraystop=1, and set fluxn 
c  to previous value fluxn(is-1).   [RWH:030428]
c-----------------------------------------------------------------------
      subroutine prep3d(t,u,deru,iraystop)
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
      double complex cflown,cnx,fr_func_noeps,dfrfunc,roots(2)
      integer info(2),ier,nsig
      double complex dhot,dhot_rlt,dhot_sum

      double precision x_ar(nbulka),y_ar(nbulka),t_av_ar(nbulka),
     .tpop_ar(nbulka),vflow_ar(nbulka),image_d
      
      double complex hotnp
      double complex dhot_sumt,dhot_sum_e,dhot_sum_i
cfor_test
      double complex k_aherm(3,3)
      double complex eps1(3,3),cez1,eps2(3,3),cez2
      double complex cd 
CENM For the iabsorp=12 damping option to work (find full complex solution of
C    D(nperp)=0), the muller root-finding subroutine is called, and either
C    Ram's relativistic function is used (for id=14) or Nelson-Melby's (for id=11)
      external fr_func_noeps  ! for dispersion function root finder (id=14)     

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

cSAP081111
      save ckvipold

c      write(*,*)'prep3d t',t
c      write(*,*)'prep3d (u(i),i=1,6)', (u(i),i=1,6)
c      write(*,*)'prep3d begin rma,zma',rma,zma
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

c----------------------------------------
c     cvac (cm/sec)
      cvac=2.997930D+10
c----------------------------------------
c     cld (cm),frgncy(GHz)
      cld=cvac/(2.d0*pi*frqncy*1.0d+09)
c----------------------------------------
c     now proposed that r0x=1 m
      r00=100.d0*r0x
      t00=cvac/(2.d0*pi*frqncy*r00)
c-----------------------------------------
c     z,r (m)
      z=u(1)
      r=u(2)
      phi=u(3)
      cnz=u(4)
      cnr=u(5)
      cm=u(6)
cSAP080906
c----------------------------------------------------
c     poloidal angle 0< theta_pol [degree] <360
c-----------------------------------------------------
c      write(*,*)'prerp3d r,z,rma,zma', r,z,rma,zma
      call theta_rz((r-rma),(z-zma),wtheta_pol(is))
      wtheta_pol(is)=(wtheta_pol(is))*180d0/pi
      
      if (is.gt.1) then
c         write(*,*)'wtheta_pol(is-1),wtheta_pol(is)',
c     &              wtheta_pol(is-1),wtheta_pol(is)
         if ((wtheta_pol(is-1).le.90d0).and.
     &       (wtheta_pol(is).ge.180.d0)) then
            wtheta_pol(is)=wtheta_pol(is)-2*180.d0
         endif
      endif
c      write(*,*)'prep3d is,wtheta_pol(is)',is,wtheta_pol(is)
c---------------------------------------
c     bmod,bz,br,bphi (Tl)
      bmod=b(z,r,phi)
      bf(1)=bz
      bf(2)=br
      bf(3)=bphi
      gam=gamma1(z,r,phi,cnz,cnr,cm)
      ds=dsin(gam)
      dc=dcos(gam)
      cnt=cn(r,cnz,cnr,cm)
      cnpar=cnt*dc
      cnper=cnt*ds

c-----waves absorption and electric field calculations----
      if ((iabsorp.eq.3).or.(iabsorp.eq.2)) then
c	absorption for lh and fw
c------------------------------------------------------------
c       electric field using the cold plasma dielectric tensor
        call tensrcld(u(1),u(2),u(3))
        cnx=dcmplx(cnper,0.d0)
        call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
c       electric field parallel to wave
        enp=(ex*cnper+ez*cnpar)/dsqrt(cnper**2+cnpar**2)
        
c-------put cold plasma tensor to w_ceps array
        do i=1,3
           do j=1,3
              w_ceps(i,j,is)=reps(i,j) !from  eps.i
           enddo
        enddo
           
c-------------------------------------------------------------
        temp_e=tempe(z,r,phi,1)
        do i=2,nbulk
          tempiar(i)=tempe(z,r,phi,i)
	enddo

        dens_e=dense(z,r,phi,1)
        z_eff=zeffrho(rho) !zeffi(z,r,phi)

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

           call absorpfd(z,r,phi,cnpar,cnper,temp_e,dens_e,tempiar,
     1                   nbulk,bmod,cnprim_e,cnprim_i,cnprim_s)

           if (ion_absorption.eq.'enabled') then 
              cnprim=cnprim_e+cnprim_i
           else
              cnprim=cnprim_e !only electron absorption
           endif 

            if (cnprim.lt.0.d0) then
             write(*,*)'Warning cnprim <0 cnprim=',cnprim
             write(*,*)'The code wil use abs(cnprim)'
             cnprim=dabs(cnprim) 
             cnprim_e=dabs(cnprim_e)
             cnprim_i=dabs(cnprim_i)
             do i=1,nbulk
                cnprim_s(i)=dabs(cnprim_s(i))
             enddo
           endif
c----------electric field calculations using the dielectric tensor
c          from absorpfd (the electron plasma with the thermal correction)
c
c          electric field using the cold plasma dielectric tensor
           cnx=dcmplx(cnper,(cnprim_e+cnprim_i))

           do i=1,3
              do j=1,3
                 w_ceps(i,j,is)=reps(i,j) !cold plasma with thermal correction
              enddo
           enddo
         
           call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
c          electric field parallel to wave
           enp=(ex*cnper+ez*cnpar)/dsqrt(cnper**2+cnpar**2)

           if (is.eq.1000000) then
c------------ test N_perp and polarization in one given ray point WF case
              write(*,*)'is=81'
              write(*,*)'cnpar,cnper',cnpar,cnper
              write(*,*)'cex,cey,cez',cex,cey,cez

              call npernpar(z,r,phi,cnpar,cnper2p,cnper2m)

              if (cnper2p.ge.0.d0) then
                 cnper_p=dsqrt(cnper2p)
                 write(*,*)'cnper_p,cnper',cnper_p,cnper
              endif

              if (cnper2m.ge.0.d0) then
                 cnper_m=dsqrt(cnper2m)
                 write(*,*)'cnper_m,cnper',cnper_m,cnper
              endif

              psi_loc=psi_rho(rho)
              r_m=r
              z_m=z  
              rho_loc=dsqrt((r_m-rma)**2+(z_m-zma)**2)
              cos_theta_pol=(r_m-rma)/rho_loc
              sin_theta_pol=(z_m-zma)/rho_loc
              if (cos_theta_pol .gt. 1.d0)  cos_theta_pol= 1.d0
              if (cos_theta_pol .lt. -1.d0)  cos_theta_pol=-1.d0
              if (sin_theta_pol.ge.0.d0) then
                 theta_pol=+dacos(cos_theta_pol)
              else  
                 theta_pol=-dacos(cos_theta_pol)
              endif

              if (theta_pol.lt.0.d0) then
                 theta_pol=theta_pol+2.d0*pi
              endif

              if (theta_pol.gt.2.d0*pi) then
                 n_theta_pol=theta_pol/(2.d0*pi)
                 theta_pol=theta_pol-2.d0*pi*n_theta_pol
              endif

c-------------calculate the hot plasma full dielectric tensor and the 
c             electric field
                 x_ar(i)=x(z,r,phi,i)
            endif !is.eq.81)
          
c-----------end test

        endif ! iabsorp=3

	if (iabsorp.eq.2) then
c----------LH absorption
           call absorplh(u,cnpar,cnper,temp_e,dens_e,tempiar
     1                  ,bz,br,bphi,nbulk,bmod,frqncy,z_eff,
     1                   cnprim_e,cnprim_i,cnprim_cl)
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
        ihermloc=iherm
c       ihermloc=1
        ihermloc=2
c-------EC absorption (from Mazzucato code)
        ioptmaz=1 ! estimation of cnper1 from coldm
        cnper1=cnper
c        write(*,*)'prep3d  bef cnperm ioptmaz=1 cnper1,cnper',
c     *  cnper1,cnper
c        write(*,*)'perp3d cnper will calculate cnper1,cnprim' 
c        write(*,*)'using the estimation of nperp from cold plasma'
        ihermloc=2
        call cnpermuz(cnpar,ihermloc,z,r,phi,cnper1,cnprim,ioptmaz)
        write(*,*)'ioptmaz=1 ihermloc,new cnper1,old cnper,new cnprim'
     *, ihermloc,cnper1,cnper,cnprim

        ioptmaz=2 ! estimation of cnper1 from input parameter cnper
c        cnper1=cnper
c        ihermloc=2
c        call cnpermuz(cnpar,ihermloc,z,r,phi,cnper1,cnprim,ioptmaz)
c        write(*,*)'ioptmaz=2 cnper1,cnper,cnprim',cnper1,cnper,cnprim

        cnprim_cl=0.d0
        cnprim_e=cnprim
        cnprim_i=0.d0

c-------electric field for mazzucato tensor

        ihermloc=iherm
c        ihermloc=1
        ihermloc=2 !full hermition + antihermition Mazzucato tens.

c        call hamiltmuz(cnpar,ihermloc,z,r,phi,cnper1,hamiltmz)
        
        call hamiltmuz_ful(cnpar,ihermloc,z,r,phi,cnper1,
     .  cnprim,hamiltmz)

c-------put Mazzucato tensor from complex N perpendicular to w_ceps array
        do i=1,3
           do j=1,3
              w_ceps(i,j,is)=reps(i,j) !from  eps.i
           enddo
        enddo

        cnx=dcmplx(cnper1,cnprim)
        call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
        enp=(ex*cnper+ez*cnpar)/dsqrt(cnper**2+cnpar**2)
c       write(*,*)'in prep3d efield mazz ex,ey,ez,enp,cnper1,cnpar',
c     1	ex,ey,ez,enp,cnper1,cnpar
      endif


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

c           write(*,*)'prep3d cnpar,cnper,cnper2p,cnper2m',
c     +     cnpar,cnper,cnper2p,cnper2m
cc----------set the data to common npercom.i for hotnp function
            call set_nperpcom(cnpar,nbulk,z,r,phi,dmas)
cc          calculates Nper(Npar) from hot plasma with the initial 
cc          iteration cnper from cold plasma

c           hotnperp_xyz=hotnp(nbulk,ibw,cnperp,cnper2p,cnper2m,K,iraystop)
c           write(*,*)'prep3d cnper,hotnperp_xyz',cnper,hotnperp_xyz
c           cnperp=hotnperp_xyz
        endif      
 
        do i=1,nbulk
           x_ar(i)=x(z,r,phi,i)
           y_ar(i)=y(z,r,phi,i)
           if(i.eq.1) y_ar(1)=-y_ar(1)
           te=tempe(z,r,phi,i)        ! kev
           t_av_ar(i)=te*1000.d0      ! ev 
           tpop_ar(i)=tpoprho(rho,i)
           vflow_ar(i)=vflowrho(rho,i)
        enddo

        if(i_im_nperp.eq.1) then
c---------calculate ImN_perp using the formula
c         ImN_perp=abs(ImD_full/dD_hermitian/dReN_perp))
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

c-------------------------------
          cnper_new=cnperp !it will be used in the electric field calculations
        endif  ! i_im_nperp.eq.1  

        if(i_im_nperp.eq.2) then 
c---------find (Im_N_perp,ReN_perp) the root of the complex dispersion relation
c         using the Newton method with numerical derivatives (the chord method)
          iter_max=100 !max number of the iterations
          iter_max=3 !max number of the iterations
c---------initial values of Im_N_perp=cnprim Re_N_perp=cnperp
          cnprim=0.d0   
c          cnper_new=cnperp
c          cnprim=cnprim_old

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

        endif !i_im_nperp.eq.2

        cnprim_cl=0.d0
        cnprim_e=dabs(cnprim)
        cnprim_i=0.d0

 30     continue
        cnprim=dabs(cnprim)

c-------electric field for Forest tensor
c       cnx=dcmplx(cnper,cnprim)
        cnx=dcmplx(cnper,0.d0)

c-------full tensor with new n_perp calculated by solver
        d=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     .  vflow_ar,cnparp,cnper,2,reps)      
        cnx=dcmplx(cnper_new,cnprim)
        cnx=dcmplx(cnper,0.d0)
 
        call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
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
           x_ar(i)=x(z,r,phi,i)
	   y_ar(i)=y(z,r,phi,i)
           if(i.eq.1) y_ar(1)=-y_ar(1)
	   te=tempe(z,r,phi,i) ! kev
	   t_av_ar(i)=te*1000.d0      ! ev 
           tpop_ar(i)=tpoprho(rho,i)
           vflow_ar(i)=vflowrho(rho,i)
        enddo
 
        d=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     .   vflow_ar,cnparp,cnperp,1,reps)
        if ((i_diskf.eq.0).or.(i_diskf.eq.5)) then
c----------usage of the analytical relativistic function and its derivatives 
           i_fkin=0
        else 
c----------usage of the mech relativistic function and its derivatives
           i_fkin=1
        endif

        call anth_rlt(x_ar(1),y_ar(1),t_av_ar(1)*1.d-3,cnparp,cnperp,
     +  n_relt_harm1,n_relt_harm2,n_relt_intgr,
     +  i_resonance_curve_integration_method,epsi,
     +  i_fkin,r,z,phi,
     +  aK)
         
        cnprimp=0.d0
        d=dhot_rlt(reps,aK,cnparp,cnperp,cnprimp)
        dham=dreal(d)
        write(*,*)'reps',reps
        write(*,*)'aK',aK
        write(*,*)'d',d
        call Ddhot(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,vflow_ar
     . ,cnparp,cnperp,reps,dd5,ddnp_h,ddnll_h,ddnp)
        
        cnprim = dabs(DIMAG(D) / DREAL(ddnp))

c	write(*,*)'prep3d forest cnper, relativistic cnprim',cnper,cnprim
        write(*,*)'dimag(d)',dimag(d)
        write(*,*)'dreal(ddnp),cnprim',dreal(ddnp),cnprim

        cnprim_cl=0.d0
        cnprim_e=cnprim
        cnprim_i=0.d0

c-------electric field for Forest tensor
        cnx=dcmplx(cnper,cnprim)
        cnx=dcmplx(cnper,0.d0)

        d=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     .  vflow_ar,cnparp,cnperp,2,reps)

        call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
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
        call tensrcld(z,r,phi)

        do i=1,nbulk
           x_ar(i)=x(z,r,phi,i)
	   y_ar(i)=y(z,r,phi,i)
           if(i.eq.1) y_ar(1)=-y_ar(1)
	   te=tempe(z,r,phi,i) ! kev
	   t_av_ar(i)=te*1000.d0      ! ev 
           tpop_ar(i)=tpoprho(rho,i)
           vflow_ar(i)=vflowrho(rho,i)
        enddo
         
        if ((i_diskf.eq.0).or.(i_diskf.eq.5)) then
c----------use of the analytical relativistic function and its derivatives 
           i_fkin=0
        else 
c----------use of the numerical relativistic function and its derivatives
           i_fkin=1
        endif
         call anth_rlt(x_ar(1),y_ar(1),t_av_ar(1)*1.d-3,cnparp,cnperp,
     +   n_relt_harm1,n_relt_harm2,n_relt_intgr,
     +   i_resonance_curve_integration_method,epsi,
     +   i_fkin,r,z,phi,
     +  aK)
        
c------ complex dispersion function calculated from the sum of
c       of the cold electron plasma dielectric tensor eps_h
c       and the relativistic electron anti-hermition dielectric tensor eps_a

        disp_func=dcold_rlt(reps,aK,cnparp,cnperp)

c-------calculate the derivative d(D_hermitian)/d(ReN_perp)
c       from the electron cold plasma dispersion function D
        ddnp=dDcold(reps,cnparp,cnperp)
        
        cnprim = dabs(DIMAG(disp_func) / DREAL(ddnp))
	
        cnprim_cl=0.d0
        cnprim_e=cnprim
        cnprim_i=0.d0
     
c-------electric field for cold plasma
        cnx=dcmplx(cnper,cnprim)
        cnx=dcmplx(cnper,0.d0)
        
        call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
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
         call absorp_relativist_disp_combined(z,r,phi,cnpar,cnper,
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
        cnx=dcmplx(cnper,cnprim)
        call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
         
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
C************ won't be much of an issue.
         errabs=1.d-6
         nsig=6
         nknown=0
         nrts=1                 ! just look for one root
         nguess=nrts 
         nnew=nrts
         itmax=50
C******* For initial guess, just search with the real part from before and just
C******* 0 imaginary part. Usually, if the ray is propagating, the damping
C******* is low anyway, so it shouldn't be too large imaginary part.
         roots(1)=dcmplx(cnper,0.0d0)
      write(*,*)'()()(()()(()()initial data'
      write(*,*)'z=',z,'r=',r,'phi=',phi
      write(*,*)'cnz=',cnz,'cnr=',cnr,'cm=',cm
      write(*,*)'rho=',rho
      print *,'====== is=',is
      print *,'cnprim: ',cnprim,' cnpar: ',cnpar,' cnper: ',cnper
      cn2=cnz**2+cnr**2+(cm/r)**2

      print *,'cn2=',cn2,' cnper(calculated) ',dsqrt(cn2-cnpar**2)
         if (is.lt.1) then
            print *,'****** is=',is
         else
            print *,'nper before: ',wnper(is-1)
         endif

C id.eq.14, call fr_func_noeps, using Ram's dielectric function
c-----------print relativistic tensor for testing
            X_e=x(z,r,phi,1)
            Y_e=y(z,r,phi,1)
            T_e=tempe(z,r,phi,1)  
            cnx=dcmplx(cnper,0.d0)
            write(*,*)'Xe,Y_e ',Xe,Y_e
            write(*,*)'cnpar',cnpar
            write(*,*)'cnx',cnx
            call Disp_Ram(T_e,cnpar,X_e,Y_e,cnx,K,d) !K is in z-y plane
                               !in Stix coordinates
            write(*,*)'prep3d K',K 
            call herm(K,K_herm)
            write(*,*)'prep3d K_herm',K_herm
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

         cnx=dcmplx(cnper,cnprim)

c-------------------------------------------------------------
cSm060719
c--------calculate dielectric tensor reps for electric field calculations
         call Disp_combined(T_e,cnpar,X_e,Y_e,cnx,reps,d)

         call efield1(cnpar,cnx,ex,ey,ez,eplus,eminus,epar)
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
c        write(*,*)'in prep3d is',is
      wye(is)=y(z,r,phi,1)
      if(nbulk.gt.1) wyi(is)=y(z,r,phi,2)

      wxe(is)=x(z,r,phi,1)
      if(nbulk.gt.1) wxi(is)=x(z,r,phi,2)

               if (nbulk.gt.3) then
                 wyi2(is)=y(z,r,phi,4)
                 wxi2(is)=x(z,r,phi,4)
               endif
c---------------------------------------------------------
c     ws (cm)
      if (is.eq.1) then
         ws(1)=0.d0  !poloidal ray distance 
         rhoold=rho
         zold=0.d0
         rold=0.d0 
         phiold=0.d0
         i_ox_conversion=0
      else
         ws(is)=ws(is-1)+dsqrt((z-zold)**2+(r-rold)**2)*r0x*100.d0
         delws=dsqrt((z-zold)**2+(r-rold)**2)*r0x*100.d0
      end if

      psi_s=psi_rho(rho)
      q_s=qsafety_psi(psi_s)
      wphi(is)=u(3)
      
      call prepebw(t,u,is)
           
      zold=z
      rold=r
      phiold=phi
c--------------------------------------------------------
c     cflown - dimensionless for E_x/E,E_y/E,E_z/E
c  !!!!now flown is calculated using Mazzucato dielectric tensor
c  !!!!and electric field was calculated using Mazzucato tensor
c  !!!!it is only for EC wave .For LH and FW it is necessery
c !!!!!to create new subroutine flown ?what tensor?
           
      if (iflux.eq.1) then
          call flown(u(1),u(2),u(3),u(4),u(5),u(6),cflown)
      endif

      if (iflux.eq.2) then
c------- flux from the cold plasma    
         xe=x(z,r,phi,1)
         ye=y(z,r,phi,1)
         rnpar=cnpar
         rnper=cnper
         romega=1/ye
         romegpe=dsqrt(xe)/ye
        write(*,*)'prep3d before grpde xe,ye',xe,ye
         nsigma=ioxm
c        nsigma=1
c        nsigma=-1
         call grpde2 (rnpar, rnper, romega, romegpe, nsigma,
     .                   rvgrpdc, redenfac, exde, eyde, ezde)
c        write(*,*)'+ redenfac',redenfac
         cflown=2.d0*redenfac
      endif !cold electron plasma flux
    
c-------------------------------------------------------

c-----calculate the group velocity
      id_old=id

      i_geom_optic_loc=i_geom_optic
      i_geom_optic=1 !to get group velocity in deru 
      call rside1(u,deru)
      i_geom_optic=i_geom_optic_loc

      vgr(1)=deru(1)
      vgr(2)=deru(2)
      vgr(3)=r*deru(3)
      
      vgrmods=vgr(1)**2+vgr(2)**2+vgr(3)**2
      if(vgrmods.gt.1.d0) then
         write(*,*)
         write(*,*) '*************************************************'
         write(*,*)'WARNING vgroup>1,   abs(vgroup) = ',dsqrt(vgrmods)  
         write(*,*) '*************************************************'
         write(*,*)
      endif

c-----the data for mnemonic.nc output file
c     the group velocity normalized to c
      wvgr_z(is)=vgr(1)
      wvgr_r(is)=vgr(2)
      wvgr_phi(is)=vgr(3)
c-----refractive index
      wn_z(is)=cnz
      wn_r(is)=cnr
      wn_phi(is)=cm/r
      
c-----------------------------------------------------------
c     vdotb -projection of the group velocity  on the
c            magnetic field multiplited by bmod
c-----------------------------------------------------------
      vdotb=vgr(1)*bz+vgr(2)*br+vgr(3)*bphi
c-----------------------------------------------------------
c     vgperps -perpendicular (to magnetic field)
c              component of the group velocity
c              in the second degree
c-----------------------------------------------------------
      vgperps=0.0d0
      do i=1,3
        vgperps=vgperps+(vgr(i)-vdotb*bf(i)/bmod**2)**2
      enddo
c----------------------------------------------------------
c     collisional damping 
c----------------------------------------------------------
      if( iabsorp_collisional.eq.1) then
         temp_e=tempe(z,r,phi,1)
         dens_e=dense(z,r,phi,1)
         z_eff=zeffrho(rho) !zeffi(z,r,phi)
         frqncy_l=frqncy
         v_gr_perp= dsqrt(vgperps)
         call absorp_collisional(temp_e,dens_e,frqncy_l,z_eff,
     &   v_gr_perp,coll_mult,
     &   cnprim_cl)
         write(*,*)'prep3d cnprim,cnprim_cl',cnprim,cnprim_cl
 2       cnprim=cnprim+cnprim_cl
      endif !iabsorp_collisional=1 
c----------------------------------------------------------
      vgrs=vgr(1)**2+vgr(2)**2+vgr(3)**2
      ckvi=dsqrt(vgperps/vgrmods)*cnprim/cld
      vgrpls=vgr(1)**2+vgr(2)**2
      vgrpol=dsqrt(vgrpls)
      wf=frqncy
c----------------------------------------------------------
c     ckvipol  (1/cm)
      vratio=dsqrt(vgperps/vgrpls)
      ckvipol=vratio*cnprim/cld
      ckvipl_e=vratio*cnprim_e/cld
      ckvipl_i=vratio*cnprim_i/cld
      ckvipl_cl=vratio*cnprim_cl/cld

      write(*,*)'ckvipl_e,vratio,cnprim_e,cld',
     +           ckvipl_e,vratio,cnprim_e,cld
 
cBH041009  Only germaine if iabsorp.eq.3:
      do kk=2,nbulk
         ckvipl_s(kk)=vratio*cnprim_s(kk)/cld
      enddo
c--------------------------------------------------------
c----------------------------------------------------------
      seikon(is)=0.d0
      spsi(is)=rho/a
c---------------------------------------------------------
c     wr and wz (cm),wphi (radian)
      wr(is)=r*r00
      wphi(is)=phi
      wz(is)=z*r00

      wnpar(is)=cnpar
c      write(*,*)'prep3d wnpar(is)',wnpar(is)
      wnper(is)=cnper
c      write(*,*)'prep3d is,cnper,wnper(is)',is,cnper,wnper(is)
      wmtor(is)=u(6)

c-------------------------------------------------------------------
c     Here delpwr(is) is the power(erg/c) in the ray
c     channel.It is equal powini_iray at antenna.
      if (is.eq.1) then
c        powj(iray) (erg/c) was calculated in cone_ec
c        powinilh(iray) (erg/c) was calculated in grill_lh
c        powini=powj or powinilh
         p=powini
         write(*,*)'prep3d is=1,powini',powini
cSAP081202
c        delpwr(is)=dexp(-2.d0*ckvipol*(ws(is)))*p
         delpwr(is)=p
         write(*,*)'is=1 delpwr(1)',delpwr(1)
         optical_depth=0.d0
cSm021101
         t_old=0.d0
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
           write(*,*)'ckvipold',ckvipold
           write(*,*)'0.5d0*(sdpwr(is-1)+salphal(is-1)+salphac(is-1))',
     &               0.5d0*(sdpwr(is-1)+salphal(is-1)+salphac(is-1))

           argexp=(-(ckvipol+ckvipold)*(ws(is)-ws(is-1)))
cSAPO81111
           ckvipold=ckvipol
cSmirnov970105
c          write(*,*)'in prep3d ckvipol,argexp',ckvipol,argexp
	   if (dabs(argexp).gt.90.d0) then
	     write(*,*)'in prep3d argexp.gt.90',argexp
             delpwr(is)=0.d0
	   else 
             pwexp=dexp(argexp)
	     if (pwexp.lt.1.d-50) then
	         write(*,*)'in prep3d pwexp.lt.1.d-50',pwexp
                 delpwr(is)=0.d0
	      endif
	   endif
           optical_depth=optical_depth+dabs(argexp)
c          write(*,*)'optical_depth',optical_depth
                  
cSm050225
c-YuP-130605: changed ray-stopping criterion from argexp>0.d0 to this: 
           if(argexp.gt.1.d-30)then
              write(*,*)'******************************************'
              write(*,*)'WARNING in prep3d.f argexp>0' 
              write(*,*)'It will give the growing ray power'
              write(*,*)'******************************************'
              argexp=0.d0
              iraystop=1
           endif
c-YuP-130605: From print-out: 
c Even though sometimes ckvipol and ckvipold both are zero,
c yet the value of  argexp= -(ckvipol+ckvipold)*(ws(is)-ws(is-1))
c is not zero (but rather a small number ~ 1e-321).
c Because of this seemingly insignificant rounding error,
c the rays were stopped prematurely.
c It only happens on Hopper/PGI, not on IntelVisualFortran. 
cSm021101
c          write(*,*)'delpwr(is-1),argexp',delpwr(is-1),argexp
           delpwr(is)=delpwr(is-1)*dexp(argexp)
c          delpwr(is)=delpwr(is-1)*dexp(-2.d0*dabs(wi)*
c     &               (ws(is)-ws(is-1))/(dsqrt(vgrs)*cvac)*    
c     &              (2.d0*pi*frqncy*1.d+9))
         else
c----------to use QL flux for absorption calculations
c          iabsorp_ql=1

           call absorbed_power_using_ql_flux(wnpar(is-1),wnper(is-1),
     &     wz(is-1),wr(is-1),wphi(is-1),
     &     fluxn(is-1),delpwr(is-1),(ws(is)-ws(is-1)),
     &     absorbed_power_ql)
           delpwr(is)=delpwr(is-1)-absorbed_power_ql

           write(*,*)'QL absorption'
           write(*,*)'delpwr(is-1),absorbed_power_ql,delpwr(is)',
     &                delpwr(is-1),absorbed_power_ql,delpwr(is)

         endif !iabsorp_ql  
c-------------------------------------------------------------------
c        reflection lost at the plasma edge
c------------------------------------------------------------------
         write(*,*)'prep3d refl_loss,irefl,irefl_old',
     &             refl_loss,irefl,irefl_old
         tot_pow_absorb_at_refl=tot_pow_absorb_at_refl+
     &           delpwr(is)*refl_loss*(irefl-irefl_old)

         write(*,*)'prep3d before refl_looss delpwr(i)',delpwr(is)

         delpwr(is)=delpwr(is)*(1.d0-refl_loss*(irefl-irefl_old))
         irefl_old=irefl 
         w_tot_pow_absorb_at_refl=tot_pow_absorb_at_refl
c-------------------------------------------------------------------    
          
cSAP090603
         write(*,*)'prep3d delpwr(is)',delpwr(is)
         write(*,*)'delpwr(1)',delpwr(1)

cSAP090603
         if((delpwr(is).gt.1.d-200).and.(delpwr(1).gt.1.d-200))then
           if(i_ox.ne.1) write(*,*)'-dlog(delpwr(is)/delpwr(1))',
     &             -dlog(delpwr(is)/delpwr(1))
         endif
 
c         if(argexp.gt.0.d0)then
c           write(*,*)'******************************************'
c           write(*,*)'WARNING in prep3d.f argexp>0' 
c           write(*,*)'It will give the growing ray power'
c          write(*,*)'******************************************'
c         endif
           
         if(delpwr(is).lt.delpwrmn*delpwr(1))then
            write(*,*)'***in prep3d delpwr(is).lt.delpwrmn*delpwr(1)**'
c           stop ray_iray calculations
            iraystop=1
         endif
      end if

      if(i_ox.eq.2) then
c---------------------------------------------------------------
c       It works for i_ox=2 case after OX mode conversion point,
c       where i_ox_conversion=1
c       It will reduce the power from O mode to X mode using 
c       transmission coefficient transm_ox

        write(*,*)'i_call_prep3d_in_output_at_i_ox_conversion_eq_1',
     &  i_call_prep3d_in_output_at_i_ox_conversion_eq_1

        if (i_call_prep3d_in_output_at_i_ox_conversion_eq_1.eq.0)goto 20

c        write(*,*)'prep3d is,delpwr(is),delpwr_o',
c     &  is,delpwr(is),delpwr_o

        if (i_call_prep3d_in_output_at_i_ox_conversion_eq_1.eq.1) then
          nrayelt=nrayelt-1
          is=is-1
          nrayelt_o_cutoff=is  !the number of ray point where O
                               ! cutoff was found
          write(*,*)'prep3d: nrayelt_o_cutoff',nrayelt_o_cutoff
          goto 20 !the first call of prep3d in output after OX conversion 
        endif

        if (i_call_prep3d_in_output_at_i_ox_conversion_eq_1.eq.2) then        
          transm_ox_loc=transm_ox
c          write(*,*)'prep3d OX transm i_ox_conversion,delpwr(is),
c     &    transm_ox',i_ox_conversion,delpwr(is),transm_ox
          delpwr_o=delpwr(is-1) !o-mode before jump
c          call OX_power_transmission(is,i_ox,i_ox_conversion,
c     &    delpwr_o,transm_ox_loc,delpwr_x)
          delpwr_x=delpwr_o*transm_ox
c          write(*,*)'prep3d delpwr_o,transm_ox_loc,delpwr_x',
c     &                      delpwr_o,transm_ox_loc,delpwr_x
          delpwr(is)=delpwr_x   !x-mode after jump
        endif

 20   continue
  
      write(*,*)'prep3d after transm delpwr(is)',delpwr(is)
      endif ! i_ox.eq.2

c     sdpwr(is)=0.d0
c----------------------------------------------------------------------
cSmirnov961122
      if(istart.eq.1) then
c        electron cyclotron waves
         wdnpar(is)=dabs(0.05d0*wnpar(1))
c        write(*,*)'in prep3d old wdnpar(is)',wdnpar(is)
      else
cSmirnov961205
c        grill conditions for the waves (LH or FW)
         wdnpar(is)=wdnpar0(iray)
c        write(*,*)'in prep3d new wdnpar(is)',wdnpar(is)
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
            write(*,*) '***in prep3d fluxn(is).lt.0.0. Set iraystop=1'
            write(*,*) '***in prep3d Set fluxn(is)=fluxn(is-1)'
            write(*,*) '***********************************************'
            write(*,*)
            if(is.gt.1) then
               fluxn(is)=fluxn(is-1)
            else
               fluxn(is)=1.
            endif
c           stop ray_iray calculations
            iraystop=1
         endif
     
c-----------------------
cBH040816      sbtot(is)=bmod*10000.d0*b0
      one=1d0
      sbtot(is)=bmod*10000.d0*b0*sign(one,feqd(1))
cBH040915:  Magnetic field components
      sb_z(is)=1.e4*bz
      sb_r(is)=1.e4*br
      sb_phi(is)=1.e4*bphi
c-----------------------------------------------------------
c     dense - dimensionless
      sene(is)=dense(z,r,phi,1)*1.0d+13
      ste(is)=tempe(z,r,phi,1)
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
      do kk=2,nbulk
         salphas(is,kk)=2.d0*ckvipl_s(kk)
      enddo
cHarvey991017 end
c------------------------------------------------------------
c     for xdraw plotter
      xarr(is)=r*dcos(phi)*r00
      yarr(is)=r*dsin(phi)*r00
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
        if(is.gt.1) then       
cyup           call check_monotonic_radius_in_ray_element(is,
cyup     &     i_non_monotonic,z_center,r_center,rho_center)
           ! YuP: output: i_non_monotonic Not used?
        endif
c        write(*,*)'rhoold,rho======',rhoold,rho
c        pause
        call p_c_prof(is,rhoold,rho,cnpar,cnper,cex,cey,cez)
      endif

      rhoold=rho

      call cpu_time(time_prep3d_2)
      time_prep3d=time_prep3d+(time_prep3d_2-time_prep3d_1)
      time_prep3d_emis=time_prep3d_emis+(time_prep3d_emis_2-
     &                                    time_prep3d_emis_1)

      return
      END ! prep3d
      
      
      

c======================================================================
c======================================================================
	  double precision
     1    FUNCTION u_res(jwave,cnpar,temp,ye)
c resonanse velosity (for nonrelativistic case)
c----------------------------------------------------------------------
c      input parameters: jwave - wave harmonic number
c                        cnpar - refractive index parallel to 
c                                magnetic field
c                        temp -  temperature in kev
c			 ye -    (omega_Be/omega)
c      output parameter: u_res in thermal velosity v/ve, ve=sqrt(2Te/me)
c-----------------------------------------------------------------------

	  implicit integer (i-n), real*8 (a-h,o-z)

c				        ve in (m/sec)
	  ve=1.87d7*dsqrt(temp)
c                               c - light velocity in (m/sec)
	  c=3.0d8
	  u_res=c/ve*(1.0d0-dble(jwave)*ye)/cnpar
c	  write(*,*)'in u_res ye,cnpar,u_res',ye,cnpar,u_res
	  return
	  END


      double precision
     1FUNCTION efficien(z_eff,u1,jwave,temp,den)
c------------------------------------------------------------------
c     RF current drive efficiency(asymptotic formulas
c     ,nonrelativistic case)
c------------------------------------------------------------------
c     input parameters: z_eff-effective charge
c               	u1-resonanse velosity  v/ve,ve=sqrt(2T_e/m_e)
c                       jwave -wave harmonic number
c                       temp - temperature in kev
c                       den - density in 10**19/m**3
c     output parameter:  efficieny in (A/cm**2)/(erg/(sec*cm**3))
c------------------------------------------------------------------
      implicit integer (i-n), real*8 (a-h,o-z)
      if(jwave.eq.1) then
c        EC wave first harmonic
         efficien=1.5d0/(5.d0+z_eff)*
     1   (u1*u1+(2.0d0+3.0d0/2.0d0/(3.d0+z_eff)))
      endif
      if(jwave.eq.0) then
c        LH wave
         efficien=2.d0/(5.d0+z_eff)*
     1   (u1*u1+(7.0d0/4.0d0+9.0d0/4.0d0/(3.d0+z_eff)))
      endif
c      write(*,*)'nondim asimptotic efficiency',efficien
c-----------------------------------------------------------
c     efficiency  in (A/cm**2)/(erg/(sec*cm**3))
c     temperature temp in kev
c     density     den in 10**19/m**3
c-----------------------------------------------------------
c     cln -Coulomb Ln
      cln=17.d0
      arg1=1.d3/temp*dsqrt(10.d0*den)	!temp KeV, den 10**13/cm**3
      cln=24.d0-dlog(arg1)
      efficien=efficien*temp/den*(17.d0/cln)*9.d0*1.d-6
      if (u1.ge.0.d0) then
         s=1.d0
      else
         s=-1.d0
      endif
cSm050329
c      efficien=efficien*s
      efficien=-efficien*s

c      write(*,*)'dim asimptotic efficien',efficien
      return
      END

c-----------------------------------------------------------------
c     this  subroutine calculates current drive efficiency
c     for the toroidal  plasma using CURBA code
c     ATTENTION:all parameters and variables input to curba
c       are real(no double precision)
c-----------------------------------------------------------------
      subroutine effcurb(z,r,zma,rma,r0x,z_eff,temp,den,jwave,cnpar,ye,
     1                   efficient)
c     input parameters: z,r -coordinates of the ray point(cm)

c                       rma,zma -coordinares of magnetic axis(cm)
c                       r0x character length (m)
c                       z_eff-effective charge
c                       temp-temperature (keV)
c                       den-density(10**13 cm**-3)
c                       jwave -wave harmonic number
c                       cnpar -parallel to magnetic field refractive
c                              index
c                     	ye-omega_Be/omega_wave
c     output parameter: J/P=efficient  in (A/m**2)/(Wt/m**3)
c------------------------------------------------------------------
c     these double precision parameters are to call function: 
c        psif(z,r) and  subroutine: zr_psith(psi,theta,z,r)
c-------------------------------------------------------------------
      implicit none
c-----input
      real*8 z,r,zma,rma,z_eff,temp,den,cnpar,ye
      real*8 r0x
      integer jwave
c-----output
      real*8 efficient

c-----locals
      real*8 psid,zd,rd,rmaxpsid,rminpsid,pid,
     &zbmin,zbmax,zb
      real*8 tc,tol,elomom,prho,ctheta,stheta,pi,theta,thtc,enpar,denom,
     *aspct,zacosarg,rjpd,rjpd0,ratjpd,efficien
      integer ig,n0,model,lh
c-----external 
      real*8 psif,bmin_psi,bmax_psi,rmax_psi,rmin_psi,b,zfac_f

      
c-------------------------------------------------------------------
c     for curba_GA
      real*8 rjpd_d,rjpd0_d,ratjpd_d,denom_d,aspct_d,enpar_d, 
     &tc_d, thtc_d, theta_d, elomom_d,zeff_d,tol_d
c-----external
c-----for efficiency transformation like in currn in toray
      real*8 phi,zfac
c-------------------------------------------------------------------
      write(*,*)'effcurb z,r,zma,rma,r0x,z_eff',
     &z,r,zma,rma,r0x,z_eff
cSAP080617
      write(*,*)'temp,den,jwave,cnpar,ye',temp,den,jwave,cnpar,ye

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
c      write(*,*)'z,r',z,r
      zd=dble(z)*0.01d0/r0x
      rd=dble(r)*0.01d0/r0x
      psid=psif(zd,rd)
c      write(*,*)'zd,rd,psid',zd,rd,psid
cSm040426
cc     rmax_psid and rmin_psid are the largest and and the least
cc     values of the radius (r) on the given flux surface
cc     psid=psid(zd,rd)
c------------------------------
c      call zr_psith(psid,0.d0,zd,rmaxpsid)
c      call zr_psith(psid,pid,zd,rminpsid)
      
c      if(abs(rmaxpsi-rminpsi).lt.1.0e-8) rminpsi=rmaxpsi-1.0e-8
c      aspct=(rmaxpsi-rminpsi)/(rmaxpsi+rminpsi)
      
c----------------------------------
c     Here :aspct is calculating using double prcision functions
c           rmax_psi(psid) and rmin_psi(psid)
c----------------------------------------
c     rmaxpsi and rminpsi are the largest and and the least
c     values of the radius (rmaxpsi and rminpsi in (m))
c     on the given flux surface psi=psi(zd,rd)
c--------------------------------------------
c      rmaxpsid=rmax_psi(psid)
c      rminpsid=rmin_psi(psid)
c      if (rmaxpsid.lt.rd)rmaxpsid=rd
c      if (rminpsid.gt.rd)rminpsid=rd
c      if(dabs(rmaxpsid-rminpsid).lt.1.0d-8) rminpsid=rmaxpsid-1.0d-8
c      aspct=(rmaxpsid-rminpsid)/(rmaxpsid+rminpsid)
cSm040426
       zbmin=bmin_psi(psid)    
       zbmax=bmax_psi(psid)
       aspct=((zbmax-zbmin)/(zbmax+zbmin))
cSm060201
       if (aspct.lt.1.d-8) aspct=1.d-8 

c       write(*,*)'psid,zbmax,zbmin,aspct',psid,zbmax,zbmin,aspct
c      write(*,*)'aspct from functionsrmin_psi, rmax_psi',aspct
c      write(*,*)'aspct',aspct,'enpar',enpar,'tc',tc,'thtc',thtc
c      write(*,*)'theta',theta,'elomom',elomom,'lh',lh
c      write(*,*)'z_eff',z_eff,'model',model,'n0',n0,'ig',ig

cSm040503
c     function b() uses the arguments z and r in [m]
      zb=b(.01d0*dble(z)/r0x,.01d0*dble(r)/r0x,0.d0)!it changes 
                                         !bz,br,bphi and rho in one.i
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
      zfac=zfac_f(dble(z*1.d-2),dble(r*1.d-2),0.d0,dble(temp))
      zfac=zfac*1.d-7 ! watt => egr/sec
      zfac=zfac*1.d-4 ! m**2 => cm**2
c      write(*,*)'zfac',zfac
      efficient=efficien*zfac
cSm050923
      efficient=-efficient 

      return
      end

      double precision function zfac_f(z,r,phi,temp_kev)
c-----It calculates the coefficient
c     that transforms the CD efficiency from
c     curba variables to A/cm**2/erg/sec*cm**3))
c     It uses the formula from currn in toray.
c-----input  z(m),r(m),phi(radians)
c            temp_kev is the electron temperature in keV    
c     It uses function x() that uses the small radius rho.
c     rho is calculated inside function b() and b() puts rho to common/one/
c     So, wee need to call function b() before using zfac_f 


      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      double precision ld
c-----externals: x
      data zelchg/4.8032d-10/, zconst/10.83259d0/,zconvr/33.33333d0/
      data cvac/2.9979d10/
c       write(*,*)'in zfac_f z,r,phi,temp_kev', z,r,phi,temp_kev
      alfa=x(z,r,phi,1) !(omega_pe/omega)**2
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
      alfa=x(z,r,phi,1) !(omega_pe/omega)**2
      zralfa=sqrt(alfa)
      zfac2=(zlds /(zelchg*alfa))*zte*2.d0/511.0d3
      zfac3=zfac1+ log (zte/zralfa)
      zfac_f=zconvr*zfac2/zfac3

      return
      end
      
      


      subroutine p_c_prof(is,rhobegin,rhoend,cnpar,cnper,
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
     
c-----input
      real*8 rhobegin,rhoend,cnpar,cnper
      integer is
      complex*16 cefldx,cefldy,cefldz !polarization

c-----locals
      real*8 delpow_s,ppow_s
      dimension delpow_s(nbulka) !Added for indiv ion contrib[BH041009]
      dimension ppow_s(nbulka)  !Added for indiv ion contrib[BH041009]

      real*8 z_r,r_r,zmar,rmar,z_effr,tempr,denr,cnparr,yer,effic_r

      real*8 r_m,z_m,temp,ye,u1,den,z_eff,hro,rholeft,rhoright,
     &rhomax,rhomin,eff_rho_max,eff_rho_min,hrho,delpower,delpow_i,
     &delpow_e,delpow_cl,del,r0,z0,r0m,z0m,psiloc,zfacgeom,aspct,
     &delcurr,rho0,poloidlen,rho0_pol,delrho,
     &ppow,pcur,ppow_e,ppow_i,ppow_cl,
     &r_bmin_right,r_bmin_left,z_bmin_right,z_bmin_left,theta,psi

      integer i,kk,jbinmin,jbinmax,j
cSAP070831
      real*8 psi_loc,cos_theta_pol,sin_theta_pol,rho_loc,theta_pol
      integer n_theta_pol,ir
c-----external
      real*8  tempe,y,u_res,dense,zeff,psi_rho,b,efficien,rhov,rhos,
     &qsafety_psi,bmin_psi,bmax_psi,dvol_dpsi,rho_lrho,
     &b_average, zeffi, zeffrho

      pi=4.d0*datan(1.d0)

c     wr and wz are in (cm) 
      r_m=wr(is)*0.01d0/r0x ! normalization
      z_m=wz(is)*0.01d0/r0x ! normalization

c       write(*,*)'p_c_prof NR,is,ieffic',NR,is,ieffic

c  radius rho is in common/one/,it was calculated in subroutine: b
c  rho is used in functions:tempe,dense,z_eff
c  The magnetic field bmode is in common/one/ .bmode is used
c  in function y(z,r,phi)


c begin if is=1      

      if(is.eq.1) then

         write(*,*)'prep3d.f p_c_prof is=1 rho=',rho
cSAP080731
        if (rho.ge.1.d0) then
c---------------------------------------------------------------
c          zero CD efficiency outside the plasma
c----------------------------------------------------------------
           eff(is)=0.d0
           goto 100
        endif
  


c        calculation of efficiency on the first step
         temp=tempe(z_m,r_m,wphi(is),1)
         ye=y(z_m,r_m,wphi(is),1)
         u1=u_res(jwave,cnpar,temp,ye)
         den=dense(z_m,r_m,wphi(is),1)
         z_eff=zeffrho(rho) !zeffi(z_m,r_m,wphi(is))
         if (ieffic.eq.1) then
c -------------------------------------------------------------------
c     calculation of current drive efficiency using asymptotic formula
c--------------------------------------------------------------------
          eff(is)=efficien(z_eff,u1,jwave,temp,den) 

c         write(*,*)'prep3d asymptotic eff',eff(is)
         endif

         if (ieffic.eq.2) then
c -------------------------------------------------------------------
c     calculation of current drive efficiency using Ehst-Karney formula
c--------------------------------------------------------------------
           write(*,*)'prep3d z_eff',z_eff
           call efKarney(z_m*100.d0*r0x,r_m*100.0d0*r0x,r0x,
     +                 z_eff,temp,den,jwave,
     1                 cnpar,eff(is))
           write(*,*)'prep3d Karney is,cnpar,eff(is)',is,cnpar,eff(is)
c----------------------------------------------------
c          ef_Karney using Bonoli subroutine:
c          For the FW NSTX case the Bonoli subroutine gave efficiency
c          which coincided exactly with the efficiency obtained 
c          calculated bu subroutine call efKarney
c-------------------------------------------------------
c          call efKarney_Bonoli(z_m*100.d0*r0x,r_m*100.0d0*r0x,r0x,
c     +                 z_eff,temp,den,jwave,
c     1                 cnpar,eff(is))
c          write(*,*)'prep3d Karney_Bonoli is,cnpar,eff(is)',
c     &                is,cnpar,eff(is)
c------------------------------------------------------------
	 endif   ! ieffic.eq.2

         if (ieffic.eq.3) then
c -------------------------------------------------------------------
c          calculation of current drive efficiency using curba
c--------------------------------------------------------------------
           z_r=(wz(is))
           r_r=(wr(is))
           zmar=(zma*100.d0*r0x) !cm
           rmar=(rma*100.d0*r0x) !cm
           z_effr=(z_eff)
           tempr=(temp)
           denr=(den)
           cnparr=(cnpar)
           yer=(ye)

      call effcurb(z_r,r_r,zmar,rmar,r0x,z_effr,tempr,denr,jwave,cnparr,
     +             yer,effic_r)
           write(*,*)'in p_c_prof after effcurb efffic_r',effic_r
           eff(is)=(effic_r)
         endif     ! ieffic.eq.3

         if (ieffic.eq.4) then
c -------------------------------------------------------------------
c          calculation of current drive efficiency using Lin_liu
c          TorGA_curgap
c--------------------------------------------------------------------
           z_r=(wz(is))
           r_r=(wr(is))
           zmar=(zma*100.d0*r0x) !cm
           rmar=(rma*100.d0*r0x) !cm
           z_effr=(z_eff)
           tempr=(temp)
           denr=(den)
           cnparr=(cnpar)
           yer=(ye)

           call eff_Lin_Liu(z_r,r_r,zmar,rmar,r0x,z_effr,tempr,denr,
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
         bmod=b(z_m,r_m,0.d0)        
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

c--------binvol(NR) calculations
         theta=0.d0 
         hro=1.d0/dble(NR-1)
         do i=1,NR-1 !for onetwo.i in [cm**3]
            rholeft=hro*(i-1)
            rhoright=rholeft+hro
            binvol(i)=voltot*(rhov(rhoright)**2-rhov(rholeft)**2)*1.d6

c            write(*,*)'NR,i,binvol(i)',NR,i,binvol(i)

           binarea(i)=areatot*(rhos(rhoright)**2-rhos(rholeft)**2)*1.d4

            psi=psi_rho(rholeft)
            call zr_psith(psi,theta,z_bmin_left,r_bmin_left)
            psi=psi_rho(rhoright)
            call zr_psith(psi,theta,z_bmin_left,r_bmin_right)

            binarea_pol(i)=pi*(r_bmin_right-r_bmin_left)*
     &                        (r_bmin_right+r_bmin_left)*1.d4

         enddo
 
         if (iabsorp.eq.3) then
            do kk=2,nbulk
               do i=1,NR
                  power_s(i,kk)=0.0d0
               enddo
            enddo
         endif

	 goto 30
      endif
c end if is=1

c      write(*,*)'rhobegin,rhoend',rhobegin,rhoend

      if(rhoend.gt.rhobegin) then
         rhomax=rhoend
         rhomin=rhobegin
      else
	 rhomin=rhoend
	 rhomax=rhobegin
      endif
cSAP080831
      if(rhomax.gt.1.d0) rhomax=1.d0
      if(rhomin.gt.1.d0) rhomin=1.d0-1.d-10    

      hrho=1.0d0/(NR-1)
      jbinmin=1
      jbinmax=NR-1

c      write(*,*)'NR,hrho,rhomin,rhomax',NR,hrho,rhomin,rhomax

      do j=1,NR-1
         if(rhomin.lt.(hrho*j)) then
           jbinmin=j
	    goto 10
	 endif
      enddo
10    continue
      do j=jbinmin,NR-1
         if(rhomax.lt.(hrho*j)) then
            jbinmax=j
	    goto 20
	 endif
      enddo
20    continue

      if (jbinmin.eq.NR) jbinmin=NR-1

c-----------------------------------------------------------
c     here delpower and allpower are in (erg/sec)
c------------------------------------------------------------
      delpower=delpwr(is-1)-delpwr(is)
      delpow_i=0.5d0*(ws(is)-ws(is-1))*
     1(delpwr(is-1)*sdpwr(is-1)+delpwr(is)*sdpwr(is))
      if (iabsorp.eq.3) then
         do kk=2,nbulk
            delpow_s(kk)=0.5d0*(ws(is)-ws(is-1))*
     1        (delpwr(is-1)*salphas(is-1,kk)+delpwr(is)*salphas(is,kk))
         enddo
      endif
      
      delpow_e=0.5d0*(ws(is)-ws(is-1))*
     1(delpwr(is-1)*salphal(is-1)+delpwr(is)*salphal(is))
cHarvey991017 end
      delpow_cl=0.5d0*(ws(is)-ws(is-1))*
     1(delpwr(is-1)*salphac(is-1)+delpwr(is)*salphac(is))

      del=delpow_i+delpow_e+delpow_cl
cHarvey970111 beg
cSmirnov050215      if (del.ne.0.d0) then
cSmirnov050215      Change 0.d0 to a small number
      if (del.gt.1.d-100) then
         delpow_e=delpow_e*delpower/del
         delpow_i=delpow_i*delpower/del
         if (iabsorp.eq.3) then
            do kk=2,nbulk
               delpow_s(kk)=delpow_s(kk)*delpower/del
            enddo
         endif
         delpow_cl=delpow_cl*delpower/del

      endif
      allpower=allpower+delpower
      allpw_e=allpw_e+delpow_e
      allpw_i=allpw_i+delpow_i
      allpw_cl=allpw_cl+delpow_cl
      
c-----------------------------------------------------------
c     r0 (cm)
      r0=0.5d0*(wr(is)+wr(is-1))

      if (rho.ge.1.d0) then
c---------------------------------------------------------------
c        zero CD efficiency outside the plasma
c----------------------------------------------------------------
         eff(is)=0.d0
         goto 110
      endif
c---- calculation of the efficiency 
      temp=tempe(z_m,r_m,wphi(is),1)
      ye=y(z_m,r_m,wphi(is),1)
      u1=u_res(jwave,cnpar,temp,ye)
      den=dense(z_m,r_m,wphi(is),1) 
      z_eff=zeffrho(rho) !zeffi(z_m,r_m,wphi(is))
   
      if (ieffic.eq.1) then
c -------------------------------------------------------------------
c       calculation of current drive efficiency using asymptotic formula
c--------------------------------------------------------------------
        eff(is)=efficien(z_eff,u1,jwave,temp,den)
c       write(*,*)'prep3d asymptotic eff',eff(is)
      endif

      if (ieffic.eq.2) then
c -------------------------------------------------------------------
c       calculation of current drive efficiency using Ehst-Karney formula
c--------------------------------------------------------------------
      call efKarney(z_m*100.0d0*r0x,r_m*100.0d0*r0x,r0x,
     & z_eff,temp,den,jwave,cnpar,
     &              eff(is))
      write(*,*)'prep3d Karney is,cnpar,eff(is)',is,cnpar,eff(is)
c----------------------------------------------------
c          ef_Karney using Bonoli subroutine:
c          For the FW NSTX case the Bonoli subroutine gave efficiency
c          which coinsided exectly wiht the efficiency obtained 
c          calculated bu subroutine call efKarney
c-------------------------------------------------------
c      call efKarney_Bonoli(z_m*100.0d0*r0x,r_m*100.0d0*r0x,r0x,
c     & z_eff,temp,den,jwave,cnpar,
c     &              eff(is))
c      write(*,*)'prep3d Karney_Bonoli is,cnpar,eff(is)',is,cnpar,eff(is)

      endif  !ieffic.eq.2

      if (ieffic.eq.3) then
c -------------------------------------------------------------------
c       calculation of current drive efficiency using curba
c-------------------------------------------------------------------
        z_r=(wz(is))
        r_r=(wr(is))
        zmar=(zma*100.d0*r0x)
        rmar=(rma*100.d0*r0x)
        z_effr=(z_eff)
        tempr=(temp)
        denr=(den)
        cnparr=(cnpar)
        yer=real(ye)

ctest
        u1=u_res(jwave,cnparr,tempr,yer)
        write(*,*)'jwave,cnpar,tempr,yer,u1',jwave,cnparr,tempr,yer,u1
        effic_r=efficien(z_effr,u1,jwave,tempr,denr)
        write(*,*)'asimptotic: z_effr,denr,effic_r',z_effr,denr,effic_r
cendtest
        call effcurb(z_r,r_r,zmar,rmar,r0x,z_effr,tempr,denr,jwave,
     +  cnparr,yer,effic_r)

        write(*,*)'in p_c_prof after effcurb is,effic_r',is,effic_r
 
        eff(is)=effic_r
      endif !ieffic.eq.3

      if (ieffic.eq.4) then
c -------------------------------------------------------------------
c       calculation of current drive efficiency using eff_Lin_Liu
c-------------------------------------------------------------------
        z_r=(wz(is))
        r_r=(wr(is))
        zmar=(zma*100.d0*r0x)
        rmar=(rma*100.d0*r0x)
        z_effr=(z_eff)
        tempr=(temp)
        denr=(den)
        cnparr=(cnpar)
        yer=real(ye)

        call effcurb(z_r,r_r,zmar,rmar,r0x,z_effr,tempr,denr,jwave,
     +  cnparr,yer,effic_r)

        write(*,*)'in p_c_prof after effcurb is,effic_r',is,effic_r

           call eff_Lin_Liu(z_r,r_r,zmar,rmar,r0x,z_effr,tempr,denr,
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
      z0=0.5d0*(wz(is)+wz(is-1))
      r0m=r0*1.d-2
      z0m=z0*1.d-2      
      bmod=b(z0m,r0m,0.d0)
c      write(*,*)'z0m,r0m,bmod,rho',z0m,r0m,bmod,rho
      psiloc=psi_rho(rho)
c-----geometric factor ~ 1/b_averaged 
      zfacgeom = -2.d0 * pi * qsafety_psi(psiloc)/
     &     (dvol_dpsi(psiloc)*dpsimax*b_average(psiloc))/100.d0


      delcurr=delpow_e*0.5d0*(eff(is-1)+eff(is))*zfacgeom !toroidal current
                                                          !created by delpow_e
c
      write(*,*)' prep3d p_c_prof delpow_e,eff(is-1),eff(is),delcurr',
     +delpow_e,eff(is-1),eff(is),delcurr
cSAP080902 to plot delta power and delta current along the ray
      delpow_e_ar(is)=delpow_e
      delcur_par_ar(is)=delpow_e*0.5d0*(eff(is-1)+eff(is))
      
c-----toroidal and poloidal CD from old genray version      
      rho0=0.5*(spsi(is)+spsi(is-1)) ! the small radius
      rho0_pol=rho_lrho(rho0)
      poloidlen=rho0_pol*totlength*100.d0     ! poloidal length cm
            
cSmirnov970105 end
      allcur=allcur+delcurr                     !total toroidal current
      
      if(rhoend.gt.rhobegin) then
          eff_rho_max=eff(is)
          eff_rho_min=eff(is-1)
      else
          eff_rho_max=eff(is-1)
          eff_rho_min=eff(is)
      endif

c      write(*,*)'jbinmin,jbinmax',jbinmin,jbinmax     
       
      if (jbinmin.eq.jbinmax) then

c         write(*,*)'jbinmon=jbinmax,delpower,power(jbinmin)',
c     .              jbinmin,delpower,power(jbinmin)

         power(jbinmin)=power(jbinmin)+delpower

cSAP080731
        write(*,*)'power(jbinmin)',power(jbinmin)
 
         cur_den_parallel(jbinmin)=cur_den_parallel(jbinmin)+
     &   delpow_e*0.5d0*(eff_rho_min+eff_rho_max)/binvol(jbinmin)

cSAP080731
        write(*,*)'delpow_e,eff_rho_min,eff_rho_max,binvol(jbinmin)',
     &             delpow_e,eff_rho_min,eff_rho_max,binvol(jbinmin)

        write(*,*)'cur_den_parallel(jbinmin)',cur_den_parallel(jbinmin)

         current(jbinmin)=current(jbinmin)+delcurr   !toroidal current from bin

         power_e(jbinmin)=power_e(jbinmin)+delpow_e
         power_i(jbinmin)=power_i(jbinmin)+delpow_i
         if (iabsorp.eq.3) then
            do kk=2,nbulk
               power_s(jbinmin,kk)=power_s(jbinmin,kk)+delpow_s(kk)
            enddo
         endif
       
         power_cl(jbinmin)=power_cl(jbinmin)+delpow_cl     

         goto 30
      endif

      delrho=rhomax-rhomin

      ppow=delpower/delrho 

      pcur=delcurr/delrho
      
cSmirnov970106 beg
      ppow_e=delpow_e/delrho
      ppow_i=delpow_i/delrho
      if (iabsorp.eq.3) then
         do kk=2,nbulk
            ppow_s(kk)=delpow_s(kk)/delrho
         enddo
      endif
      
      ppow_cl=delpow_cl/delrho
cSmirnov970106 end

c------------------------------------------------------------
c     power (erg/sec), current(A)
c------------------------------------------------------------

c      write(*,*)'jbinmin,power(jbinmin),ppow,(hrho*jbinmin-rhomin)',
c     &jbinmin,power(jbinmin),ppow,(hrho*jbinmin-rhomin)

      power(jbinmin)=power(jbinmin)+ppow*(hrho*jbinmin-rhomin)     

c      write(*,*)'power(jbinmin)',power(jbinmin)

cSmirnov970106 beg
      power_e(jbinmin)=power_e(jbinmin)+ppow_e*(hrho*jbinmin-rhomin)
      power_i(jbinmin)=power_i(jbinmin)+ppow_i*(hrho*jbinmin-rhomin)
      
      if(ppow_i.lt.0.d0)then
         write(*,*)'p_c_prof: jbinmin,ppow_i=',
     +                        jbinmin,ppow_i,delpow_i/delrho
         pause
      endif
      if((hrho*jbinmin-rhomin).lt.0.d0)then
         write(*,*)'p_c_prof: jbinmin,(hrho*jbinmin-rhomin)=',
     +                        jbinmin,(hrho*jbinmin-rhomin)
         pause
      endif
      if((rhomax-hrho*(jbinmax-1)).lt.0.d0)then
         write(*,*)'p_c_prof: jbinmin,(rhomax-hrho*(jbinmax-1))=',
     +                        jbinmin,(rhomax-hrho*(jbinmax-1))
         pause
      endif
      
      power_cl(jbinmin)=power_cl(jbinmin)+ppow_cl*(hrho*jbinmin-rhomin)
      if (iabsorp.eq.3) then
         do kk=2,nbulk
            power_s(jbinmin,kk)=
     &           power_s(jbinmin,kk)+ppow_s(kk)*(hrho*jbinmin-rhomin)
         enddo
      endif

cSmirnov970106 end
c      write(*,*)'p_c_prof delpow_e',delpow_e

      cur_den_parallel(jbinmin)=cur_den_parallel(jbinmin)+
     &(delpow_e/delrho)*(hrho*jbinmin-rhomin)/binvol(jbinmin)*
     &0.5d0*(eff_rho_min+
     &      (eff_rho_min+(eff_rho_max-eff_rho_min)*
     &                   (hrho*jbinmin-rhomin)/(delrho)))

      current(jbinmin)=
     1 current(jbinmin)+pcur*(hrho*jbinmin-rhomin)  !toroidal current from bin
 
      power(jbinmax)=power(jbinmax)+ppow*(rhomax-hrho*(jbinmax-1))

cSmirnov970106 beg
      power_e(jbinmax)=power_e(jbinmax)+ppow_e*(rhomax-hrho*(jbinmax-1))
      power_i(jbinmax)=power_i(jbinmax)+ppow_i*(rhomax-hrho*(jbinmax-1))
      power_cl(jbinmax)=power_cl(jbinmax)+
     1                  ppow_cl*(rhomax-hrho*(jbinmax-1))
      if (iabsorp.eq.3) then
         do kk=2,nbulk
            power_s(jbinmax,kk)=
     &          power_s(jbinmax,kk)+ppow_s(kk)*(rhomax-hrho*(jbinmax-1))
         enddo
      endif
cSmirnov970106 end

      cur_den_parallel(jbinmax)=cur_den_parallel(jbinmax)+
     &(delpow_e/delrho)*(rhomax-hrho*(jbinmax-1))/binvol(jbinmax)*
     &0.5d0*(eff_rho_max+
     &   (eff_rho_min+(eff_rho_max-eff_rho_min)*
     &   (rhomax-hrho*(jbinmax-1))/(delrho)))
     
      current(jbinmax)=
     1	  current(jbinmax)+pcur*(rhomax-hrho*(jbinmax-1)) !toroidal current from bin

      if(jbinmax.gt.(jbinmin+1)) then
         do j=(jbinmin+1),(jbinmax-1)
c            write(*,*)'j,power(j),ppow,hrho',j,power(j),ppow,hrho
            power(j)=power(j)+ppow*hrho
c            write(*,*)'power(j)',power(j)
cSmirnov970106 beg
            power_e(j)=power_e(j)+ppow_e*hrho
            power_i(j)=power_i(j)+ppow_i*hrho
            power_cl(j)=power_cl(j)+ppow_cl*hrho
            if (iabsorp.eq.3) then
               do kk=2,nbulk
                  power_s(j,kk)=power_s(j,kk)+ppow_s(kk)*hrho
               enddo
            endif

cSmirnov970106 end
            cur_den_parallel(j)=cur_den_parallel(j)+
     &      (delpow_e/delrho)*hrho/binvol(j)*
     &      (eff_rho_min+(eff_rho_max-eff_rho_min)*
     &      (hrho*(j-0.5d0)-rhomin)/delrho)
     
            current(j)=current(j)+pcur*hrho             !toroidal current from bins

         enddo
      endif

30    continue
      
c99    return
      END




c-----------------------------------------------------------------------
c     this subroutine SONETWO is called after each ray finished
c     it calculates power (spower(i)) and current
c     (scurrent(i)) profiles
c     as sum the same profiles for all rays.
c-----------------------------------------------------------------------
c     input parameter:iray -number of the ray is in common/cone/
c-----------------------------------------------------------------------
      subroutine sonetwo
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'onetwo.i'
      INCLUDE 'cone.i'
c for test
      INCLUDE 'write.i'
      INCLUDE 'one.i'     
      INCLUDE 'three.i'    
      include 'five.i' ! contains rmax
      
c      write(*,*)'*******************in sonetwo********************start'
c      write(*,*)'Absorbed power (one ray) as a func. of rho_bin_center:'
c      write(*,*)
c     + 'i, rho_bin_center,  power_e,  power_i,  power_cl,  power(total)'
c      do i=1,NR-1
c         write(*,'(i4,5e12.3)')
c     1	i,rho_bin_center(i),power_e(i),power_i(i),power_cl(i),power(i)
c      enddo
c      write(*,*)
c     + 'i, rho_bin_center,  power_e,  power_i,  power_cl,  power(total)'
c      write(*,*)'*******************in sonetwo*************************'
c-----------------------------------------------------------------------
c     spower (erg/sec), scurrent(A)
c-----------------------------------------------------------------------
      do i=1,NR-1
         spower(i)=spower(i)+power(i)
         spower_e(i)=spower_e(i)+power_e(i)
         spower_i(i)=spower_i(i)+power_i(i)
         spower_cl(i)=spower_cl(i)+power_cl(i)
         scurrent(i)=scurrent(i)+current(i)      !toroidal current
         s_cur_den_parallel(i)=s_cur_den_parallel(i)+cur_den_parallel(i)
      enddo
      if (iabsorp.eq.3) then
         do kk=1,nbulk
            do i=1,NR-1
               spower_s(i,kk)=spower_s(i,kk)+power_s(i,kk)
            enddo
         enddo
      endif

      do j=1,NZgrid ! YuP[Nov-2014] Sum power over (R,Z) grid:
      do i=1,NRgrid ! add power from the given ray to total spwr(R,Z)
         spwr_rz_e(i,j)= spwr_rz_e(i,j) +pwr_rz_e(i,j) !coll-less.damp e
         spwr_rz_i(i,j)= spwr_rz_i(i,j) +pwr_rz_i(i,j) !coll-less.damp i
         spwr_rz_cl(i,j)= spwr_rz_cl(i,j) +pwr_rz_cl(i,j) ! coll.damp
         do kk=1,nbulk ! Collisionless, for each plasma species:
            spwr_rz_s(i,j,kk)= spwr_rz_s(i,j,kk) +pwr_rz_s(i,j,kk)
         enddo
      enddo
      enddo

c----------------------------------------------------------------------
c   test
c   absorbed power on the ray abspwer
      abspwer=0.d0
      abspw_e=0.d0
      abspw_i=0.d0
      abspw_cl=0.d0
      curntray=0.d0
      curtrray=0.d0     !toroidal current
      curplray=0.d0     !poloidal current
      do i=1,NR
         abspwer=abspwer+power(i) ! Total
         abspw_e=abspw_e+power_e(i) ! e
         abspw_i=abspw_i+power_i(i) ! i
         abspw_cl=abspw_cl+power_cl(i) ! collisional
         curntray=curntray+current(i) !toroidal current from one ray
      end do
      write(*,*)
      write(*,*)'   Absorbed power for the ray (sum over all rho):'
      write(*,*)
     +'   abspw_e,    abspw_i,    abspw_cl,   abspwer,
     +    delpwr(1)-delpwr(nrayelt)'
      write(*,'(5e12.3)')
     + abspw_e, abspw_i, abspw_cl, abspwer, delpwr(1)-delpwr(nrayelt)
c      absdpwr=0.0d0
c      write(*,*)'in sonetwo nrayelt',nrayelt
c      do i=1,nrayelt
c        write(*,*)'i,delpwr(i)',i,delpwr(i)
c      end do
c      absdpwr=delpwr(1)-delpwr(nrayelt)
c      write(*,*)'absdpwr',absdpwr
      write(*,*)'*******************in sonetwo**********************end'
      !pause
      return
      end


c-----------------------------------------------------------------------
c     this subroutine DNONETWO calculates power and current density
c     profiles on rho. Arrays: powden(NR) (erg/(cm**3*c)
c                           and currden(NR).
c-----------------------------------------------------------------------
      subroutine dnonetwo
!      implicit integer (i-n), real*8 (a-h,o-z)
      implicit none
      include 'param.i'
      include 'onetwo.i'
      include 'gr.i'
      include 'rho.i'
      include 'three.i'
      include 'five.i'
      include 'one.i'

c-----externals
      real*8 rhov,rhos,rho_lrho,psi_rho,ias1r,b

c-----locals
      real*8 pwtot_s
      dimension pwtot_s(nbulka)  !Added for indiv ion contrib[BH041009]
      real*8 hro,rholeft,rhoright,binplarea,h_rho,psi,f_eqd,
     &b_av,bs_av,r_av,dr_av,drs_av,theta,r_bmin,z_bmin,b_pol_bmin,
     & powertot,rho0,rho0_pol,poloidlen,
     &powrtot1,currtot1,pwtot_e,pwtot_i,pwtot_cl,
     &rho_l,psi_l  
      integer i,kk,k,idx,nr4,j
     
      pi=4.0d0*datan(1.0d0)
c-----------------------------------------------------------------------
cSAP080731
c     write(*,*)'in dnonetwo'
c     write(*,*)'spower'
c     write(*,*)(spower(i), i=1,NR-1)
c     write(*,*)'spower_e'
c     write(*,*)(spower_e(i), i=1,NR-1)
c     write(*,*)'spower_i'
c     write(*,*)(spower_i(i), i=1,NR-1)
c     write(*,*)'spower_cl'
c     write(*,*)(spower_cl(i), i=1,NR-1)
c      write(*,*)'spower_s'
c      write(*,*)((spower_s(i,kk), i=1,NR-1),kk=1,nbulk)
c      write(*,*)'scurrent'
c      write(*,*)(scurrent(i), i=1,NR-1)      
c     write(*,*)'volume and area total:',voltot,areatot

c----------------------------------------------------------------------
c     spower (erg/sec), scurrent(A), binvol(cm**3),binarea(cm**2)
c     powden(erg/(sec*cm**3)),currden(A/cm**2)
c     voltot (m**3), areatot (m**2)
c----------------------------------------------------------------------
      hro=1.d0/dble(NR-1)
      do i=1,NR-1
         rholeft=hro*(i-1)
         rhoright=rholeft+hro

        write(*,*)'NR,i,rholeft,rhoright',NR,i,rholeft,rhoright

         binvol(i)=voltot*(rhov(rhoright)**2-rhov(rholeft)**2)
     1          *1000000.d0
c         write(*,*)'rhov(rhoright)',rhov(rhoright)
c         write(*,*)'rhov(rholeft)',rhov(rholeft)
         write(*,*)'voltot(m**3),binvol(cm**3',voltot,binvol(i)

         binarea(i)=areatot*(rhos(rhoright)**2-rhos(rholeft)**2)
     1            *10000.d0
c         write(*,*)'rhos(rhoright)',rhos(rhoright)
c         write(*,*)'rhos(rholeft)',rhos(rholeft)
c         write(*,*)'areatot(m**2),binarea(cm**2)',areatot,binarea(i)

         pollen(i)=
     1        totlength*50.d0*(rho_lrho(rhoright)+rho_lrho(rholeft))
c         write(*,*)'rho_lrho(rhoright)',rho_lrho(rhoright)
c         write(*,*)'prep3d totlength(m),pollen(cm)',
c     1   totlength,,pollen(i)
         binplarea=binvol(i)/pollen(i)
c	 write(*,*)'i',i, 'binvol!!(cm**3)',binvol(i)
c         write(*,*)'binarea cm**2',binarea(i)
c         write(*,*)'poloidal pollen (cm)',pollen(i)
c         write(*,*)'binplarea cm**2',binplarea

	 powden(i)=spower(i)/binvol(i)
cSmirnov970106 beg
	 powden_e(i)=spower_e(i)/binvol(i)
	 powden_i(i)=spower_i(i)/binvol(i)


         if(iabsorp.eq.3) then
            do kk=2,nbulk
               powden_s(i,kk)=spower_s(i,kk)/binvol(i)
            enddo
         endif

	 powden_cl(i)=spower_cl(i)/binvol(i)

cSmirnov970106 end
	 currden(i)=scurrent(i)/binarea(i)    !toroidal current density      
      enddo  ! i=1,NR-1

c      do kk=2, nbulk
c         write(*,*)'prep3d: first powden_s(i,kk)=', 
c     1                   (powden_s(i,kk),i=1,NR-1)
c      enddo
         
c       write(*,*)'spower'
c       write(*,*)(spower(i), i=1,NR-1)
c       write(*,*)'scurrent'
c       write(*,*)(scurrent(i), i=1,NR-1)

c$$$      do i=2,NR-1
c$$$cSm040426
c$$$         powden_s(i)=0.5d0*(powden(i-1)+powden(i))
c$$$
c$$$         powden_e_s(i)=0.5d0*(powden_e(i-1)+powden_e(i))
c$$$         powden_i_s(i)=0.5d0*(powden_i(i-1)+powden_i(i))
c$$$         powden_cl_s(i)=0.5d0*(powden_cl(i-1)+powden_cl(i))
c$$$
c$$$         currden_s(i)=0.5d0*(currden(i-1)+currden(i))
c$$$      end do
c$$$
c$$$      powden_s(1)=powden(1)
c$$$      powden_e_s(1)=powden_e(1)
c$$$      powden_i_s(1)=powden_i(1)
c$$$      powden_cl_s(1)=powden_cl(1)
c$$$      currden_s(1)=currden(1)
c$$$
c$$$      powden_s(NR)=powden(NR-1)
c$$$      powden_e_s(NR)=powden_e(NR-1)
c$$$      powden_i_s(NR)=powden_i(NR-1)
c$$$      powden_cl_s(NR)=powden_cl(NR-1)
c$$$      currden_s(NR)=currden(NR-1)

      write(*,*)' powden_e',(powden_e(i), i=1,NR)

      do k=2,nbulk
        write(*,*)'k=',k,'powden_s'
        write(*,*)(powden_s(i,k), i=1,NR)
      enddo

c      write(*,*)'curdens'
c      write(*,*)(currden_s(i), i=1,NR)

c-----------------------------------------------------------
c     CD calculation using GA memo
c-----------------------------------------------------------
c     cur_den_onetwo=<j_parallel.B>/B_0=<j_parallel><B**2>/<B>/B_0
c----------------------------------------------------------
      h_rho=1.d0/(NR-1)
      idx=0 
      nr4=nx+4
 
      write(*,1010)'cur_den_tor=p_tor*cur_den_par, '
     &   // 'p_tor = drs_av*f_eqd/(b_av*dr_av)'
      write(*,1010)'cur_den_pol=p_pol*cur_den_par, '
     &   // 'p_pol = b_pol_bmin/b_av'

      write(*,1010)'i rho_bin_center powden_e   cur_den_par  p_tor'
     & // '      cur_den_tor'
     & // '   p_pol    cur_den_pol'
    
 1010    format(/,1x,a)
 1011    format(i3,7(1pe12.4))

      do i=1,NR-1
        rho_l=h_rho*(i-0.5d0)    
        psi_l=psi_rho(rho_l)
        write(*,*)'i,rho_l,psi_l',i,rho_l,psi_l
cSAP080321 argument psi in f_eqd was changed from psi to psi_l
        f_eqd=ias1r(txf,nx,nr4,cx,idx,psi_l)

        call average_variables(psi_l,b_av,bs_av,r_av,dr_av,drs_av)

c        write(*,*)'after average_variables psi_l',psi_l
        
c        write(*,*)'b_av,bs_av,r_av,dr_av,drs_av',
c     &             b_av,bs_av,r_av,dr_av,drs_av

c        write(*,*)'s_cur_den_parallel(i)',s_cur_den_parallel(i)

        s_cur_den_onetwo(i)   = s_cur_den_parallel(i)*bs_av/(b_av*beqd)       
        s_cur_den_toroidal(i) = s_cur_den_parallel(i)*drs_av*f_eqd/
     &                          (b_av*dr_av)

c        write(*,*)'i,rho_l,s_cur_den_parallel(i)',
c     &  i,rho_l,s_cur_den_parallel(i)
c        write(*,*)'drs_av,f_eqd,b_av,dr_av,s_cur_den_toroidal(i)', 
c     &             drs_av,f_eqd,b_av,dr_av,s_cur_den_toroidal(i)

        theta=0.d0

c        write(*,*)'before zr_psith psi_l,theta ',psi_l,theta
 
        call zr_psith(psi_l,theta,z_bmin,r_bmin)

c        write(*,*)'psi_l,theta,z_bmin,r_bmin ',
c     &  psi_l,theta,z_bmin,r_bmin

        bmod=b(z_bmin,r_bmin,0.d0)
        b_pol_bmin=dsqrt(bz**2+br**2) ! poloidal B at the point
                                      ! with minimal B
                                      ! with theta poloidal =0

c        write(*,*)'bz,br,b_pol_bmin ',bz,br,b_pol_bmin 

        s_cur_den_poloidal(i) = s_cur_den_parallel(i)*b_pol_bmin/
     &                          b_av
c        write(*,*)'s_cur_den_onetwo(i)', s_cur_den_onetwo(i)
c        write(*,*)'s_cur_den_toroidal(i)',s_cur_den_toroidal(i)
c        write(*,*)'s_cur_den_poloidal(i)',s_cur_den_poloidal(i) 

c        write(*,*)'b_pol_bmin,b_av,b_pol_bmin/b_av ',
c     &             b_pol_bmin,b_av,b_pol_bmin/b_av
c        write(*,*)'drs_av,f_eqd,b_av,dr_av ',drs_av,f_eqd,b_av,dr_av 
c        write(*,*)'drs_av*f_eqd/(b_av*dr_av)',drs_av*f_eqd/(b_av*dr_av)
c----------------------------------------------------------------
         write(*,1011)i,rho_l,powden_e(i),s_cur_den_parallel(i),
     &   drs_av*f_eqd/(b_av*dr_av),s_cur_den_toroidal(i),
     &   b_pol_bmin/b_av,s_cur_den_poloidal(i)
      enddo 
c-----------------------------------------------------------
c     powertot (erg/sec), currtot(A)
c-----------------------------------------------------------
     
cSm040426begin
      powertot=0.0d0
      powtot_e=0.0
      powtot_i=0.0d0
      if (iabsorp.eq.3) then
         do kk=2,nbulk
            powtot_s(kk)=0.0d0
         enddo
      endif
      powtot_cl=0.0d0
      currtot=0.0d0
   
      do i=1,NR-1
c     integration formulas ****INT1**
         powertot=powertot+powden(i)*
     1            (rhov(hro*i)**2-rhov(hro*(i-1))**2)*
     2	          voltot*1000000.0d0
         powtot_e=powtot_e+powden_e(i)*
     1            (rhov(hro*i)**2-rhov(hro*(i-1))**2)*
     2	          voltot*1000000.0d0
         powtot_i=powtot_i+powden_i(i)*
     1            (rhov(hro*i)**2-rhov(hro*(i-1))**2)*
     2	          voltot*1000000.0d0
         if (iabsorp.eq.3) then
            do kk=2,nbulk
               powtot_s(kk)=powtot_s(kk)+powden_s(i,kk)*
     1              (rhov(hro*i)**2-rhov(hro*(i-1))**2)*
     2              voltot*1000000.0d0
            enddo
         endif
         powtot_cl=powtot_cl+powden_cl(i)*
     1            (rhov(hro*i)**2-rhov(hro*(i-1))**2)*
     2	          voltot*1000000.0d0

         currtot=currtot+currden(i)*
cSm030428     1	         (rhov(hro*i)**2-rhov(hro*(i-1))**2)*
     1	         (rhos(hro*i)**2-rhos(hro*(i-1))**2)*
     2	          areatot*10000.0d0

c--------poloidal total current
         rho0=hro*(i+0.5d0)                 !  small radius
         rho0_pol=rho_lrho(rho0)
c        poloidlen=rho0_pol*totlength*100.d0*r0x! poloidal length cm
cSm040728         poloidlen=rho0_pol*totlength*100.d0 ! poloidal length cm
         poloidlen=totlength*50.d0*(rho_lrho(hro*i)+rho_lrho(hro*(i-1)))
         binvol(i)=voltot*(rhov(hro*i)**2-rhov(hro*(i-1))**2)
     1          *1000000.d0
      enddo
     
      parallel_cur_total=0.d0
      toroidal_cur_total=0.d0
      poloidal_cur_total=0.d0
      do j=1,NR-1
        parallel_cur_total=parallel_cur_total+
     &                     s_cur_den_parallel(j)*binarea(j)
        toroidal_cur_total=toroidal_cur_total+
     &                     s_cur_den_toroidal(j)*binarea(j)
        poloidal_cur_total=poloidal_cur_total+
     &                     s_cur_den_poloidal(j)*binarea_pol(j)
      enddo


      write(*,*)'parallel_cur_total, toroidal_cur_total ',
     &parallel_cur_total, toroidal_cur_total 
      write(*,*)'poloidal_cur_total ',poloidal_cur_total

      write(*,*)'INT1 testing DNONETWO. powertot=erg/sec',powertot,
     1' currtot=A',currtot
 
      write(*,*)'testing 1 DNONETWO powtot_e,powtot_i,powtot_cl',
     &     powtot_e,powtot_i,powtot_cl
      write(*,*)'powtot_s(1:nbulk)=',(powtot_s(kk),kk=1,nbulk)
cSm040426end


      powrtot1=0.d0
      currtot1=0.d0
      
cSmirnov970106 beg
      pwtot_e=0.d0
      pwtot_i=0.d0
      if (iabsorp.eq.3) then
         do kk=2,nbulk
            pwtot_s(kk)=0.d0
         enddo
      endif
      pwtot_cl=0.d0
cSmirnov970106 end
      do i=1,NR-1
         powrtot1=powrtot1+spower(i)
cSmirnov970106 beg
         pwtot_e=pwtot_e+spower_e(i)
         pwtot_i=pwtot_i+spower_i(i)
         if (iabsorp.eq.3) then
            do kk=2,nbulk
               pwtot_s(kk)=pwtot_s(kk)+spower_s(i,kk)
            enddo
         endif
         pwtot_cl=pwtot_cl+spower_cl(i)
cSmirnov970106 end
         currtot1=currtot1+scurrent(i)             !totaql toroidal current        
      enddo
      write(*,*)'INT2 testing DNONETWO. powrtot1=erg/sec',powrtot1,
     1' currtot1=A',currtot1
      
      write(*,*)'testing 2 DNONETWO pwtot_e,pwtot_i,pwtot_cl',
     &     pwtot_e,pwtot_i,pwtot_cl
      write(*,*)'pwtot_s(1:nbulk)=',(pwtot_s(kk),kk=1,nbulk)

      write(*,*)'powden: ',powden
      write(*,*)'powden_e: ',powden_e
      write(*,*)'powden_i: ',powden_i
      do kk=2,nbulk
         write(*,*)'kk,powden_s(,kk): ',kk,(powden_s(i,kk),i=1,NR-1)
      enddo
         

c -----------------------------------------------------------------
c     normalization of density profiles
c     it dependents on numerical integration formulas
c     In the case  ***INT**** we have the following normalization:

      write(*,*)'NR',NR

      do i=1,NR
cHarvey970112
         if (powertot.ne.0.d0) then
             powden(i)=powden(i)*powrtot1/powertot
         endif
cSmirnov970106 beg
         if (powtot_e.ne.0.d0) then
             powden_e(i)=powden_e(i)*pwtot_e/powtot_e
         endif
         if (powtot_i.ne.0.d0) then
           powden_i(i)=powden_i(i)*pwtot_i/powtot_i
         endif
         if (iabsorp.eq.3) then
         do kk=2,nbulk
            if (powtot_s(kk).ne.0.d0) then
               powden_s(i,kk)=powden_s(i,kk)*pwtot_s(kk)/powtot_s(kk)
            endif
         enddo
         endif
         if (powtot_cl.ne.0.d0) then
           powden_cl(i)=powden_cl(i)*pwtot_cl/powtot_cl
         endif
cSmirnov970106 end
         if (currtot.ne.0.d0) then
            currden(i)=currden(i)*currtot1/currtot
         endif        
      enddo
c------------------------------------------------------------------

      do kk=2, nbulk
         write(*,*)'prep3d: second powden_s(i,kk)=', 
     1        (powden_s(i,kk),i=1,NR-1)
      enddo


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
      if (iabsorp.eq.3) then
cSAP090306
c         write(*,10021)(kk-1,pwtot_s(kk),kk=2,nbulk)
c10021    format(' total absrbd power_s(erg/sec) for ion species',
c     1       i2,' =',1pe14.6)
        write(*,10021)(kk-1,pwtot_s(kk)*1.d-7,kk=2,nbulk)
10021    format(' total absrbd power_s(power) for ion species',
     1       i2,' =',1pe14.6)
      endif
cSAP090306
c      write(*,1003)pwtot_cl
c 1003 format(' total absorbed power_cl(erg/sec)=',1pe14.6)
      write(*,1003)pwtot_cl*1.d-7
 1003 format(' total absorbed power_cl(power)=',1pe14.6)

cSAP080928
c      write(*,1004)currtot1
 1004 format(' total tor curr drive per LL memo[Cohen/(1+eps)] (A)=',
     1     1pe14.6)
 1005 format(' number of rays with bad initial conditions ='1i6)
cSAP090306
c 1006 format(' total launched power with good initial conditions',/,
c     1 '                        (egs/sec)=',1pe14.6)
 1006 format(' total launched power with good initial conditions',/,
     1 '                        (watt)=',1pe14.6)
ctest      
c      write(*,*)' testing DNONETWO. allpower=',allpower

cSAP090306
      write(*,1007)parallel_cur_total
      write(*,1008)toroidal_cur_total
      write(*,1009)poloidal_cur_total
 1007 format('total parallel current parallel_cur_total (A)=',
     &    1pe14.6)
 1008 format('total toroidal current toroidal_cur_total (A)=',
     &    1pe14.6)
 1009 format('total poloidal current poloidal_cur_total (A)=',
     &    1pe14.6)



      return ! dnonetwo
      END





      subroutine grpde2 (npar, nper, omega, omegpe, nsigma,
     .                   vgrpdc, edenfac, exde, eyde, ezde)
c
      complex          rootm1,exde,eyde,ezde,d11,d12,d13,d22,d33,denom
      real             npar,npar2,nper,nper2,n2,nperol
****  double precision deps1,deps2,deps3,domega,domegpe,
****  double precision dquad
      real             dquad
      real deps1,deps2,deps3,domega,domegpe,aquad,bquad,cquad
c
c     calculate group velocity, wave polarizations, and energy density factor,
c     using cold plasma relations. (M O'Brien's version)
c
      nperol = nper
      rootm1 = (0.0, 1.0)
c
c     DOUBLE PRECISION is used since a problem was encountered in the
c                 real is used since a problem was encountered in the
c     following subtraction within the formation of rootqd.
c
c --- double precision mentioned above was disabled 21 Aug 94 by Joe Freeman
c --- not needed since the WHOLE PROGRAM is done with 64-bit real arithmetic
c --- on CRAY this is automatic, elsewhere a compilation switch ensures this
c
      domega  = omega
      domegpe = omegpe
c      write(*,*)'in grpde2 npar,domega,domegpe,nsigma',
c     1                     npar,domega,domegpe,nsigma
      deps1   = 1.0 - domegpe*domegpe/        (domega*domega-1.0)
      deps2   =      -domegpe*domegpe/(domega*(domega*domega-1.0))
      deps3   = 1.0 - domegpe*domegpe/        (domega*domega)
      eps1    = deps1
      eps2    = deps2
      eps3    = deps3
c
      npar2 = npar*npar
      aquad = deps1
      bquad = - ( (deps1-npar2)*(deps1+deps3) - deps2*deps2 )
      cquad = deps3 * ( (deps1-npar2)*(deps1-npar2) - deps2*deps2 )
      dquad = bquad*bquad - 4.0*aquad*cquad
      
      if (dquad .gt. 0.0) then
        rootqd = SQRT (dquad)
      else
        rootqd = 0.0
      end if
      nper2=(-bquad+nsigma*SIGN (1.0,omega-1.0)*rootqd)/(2.0*aquad)

c
c     Let nper2 become (only) slightly negative:
c
      if (nper2 .le. 0.0) then
        if (nper2 .le. -1.0e-2) then
          STOP 'subroutine GRPDE2: nper2 is too negative'
        else
          nper2=1.0e-10
          nper=1.0e-5
        end if
      else
        nper = SQRT (nper2)
      end if
c
      n2=nper2+npar2

      d11=eps1-npar2
      d12=-rootm1*eps2
      d13=npar*nper
      d22=eps1-n2
      d33=eps3-nper2
c 
      
      exde=d22/d12
      ezde=-d13*exde/d33
      denom = SQRT (exde*conjg(exde)+1.0+ezde*conjg(ezde))
      exde=exde/denom
      eyde=1.0/denom
      ezde=ezde/denom
c      write(*,*)'in grpde2  exde,eyde,ezde',exde,eyde,ezde

c
c calculate derivatives of eps's wrt omega
c
      dep1do=omegpe*omegpe*2.0*omega/
     .        ((omega*omega-1.0)*(omega*omega-1.0))
      dep2do=-eps2/omega-eps2*2.0*omega/(omega*omega-1.0)
      dep3do=omegpe*omegpe*2.0/(omega**3)
c
c find derivatives of D the dispersion reln wrt omega, kper and kpar
c
      dddnpp=-(2.0*nper)*( (eps1-npar2)*(eps1-n2+eps3-nper2)
     .                   - eps2*eps2 + npar2*(eps1-n2-nper2) )
      dddnpl=-(2.0*npar)*( (eps3-nper2)*(eps1-npar2+eps1-n2)
     .                   + nper2*(eps1-n2-npar2) )
c
      dddome=dddnpp*(-nper/omega) + dddnpl*(-npar/omega)
     .     + dep1do*( (eps3-nper2)*(eps1-n2+eps1-npar2) - nper2*npar2 )
     .     + dep2do*( -2.0*eps2*(eps3-nper2) )
     .     + dep3do*( (eps1-npar2)*(eps1-n2) - eps2*eps2 )
c
c find vgroup/c
c
      vgrpdc = SQRT (dddnpp*dddnpp+dddnpl*dddnpl) / ABS (dddome*omega)
      
c     energy density factor
c
c      1  | ~ | 2    1  ~*      ~       ~     ~
c   =  -  | B |    + -  E .tens.E  with B and E normalized to |E|
c      2  |   |      2       _      _
c                         d |     h  |    h
c    and tens the tensor  - |  w e   |   e  is the hermitian component
c                         dw|_      _!
c
c    of the dielectric tensor
c
      bmod2 = REAL ( n2*eyde*conjg(eyde)
     .    + (npar*exde-nper*ezde)*conjg(npar*exde-nper*ezde) )
c      write(*,*)'in grpde2  bmod2',bmod2
      tens11=eps1+omega*dep1do
      tens12=-(eps2+omega*dep2do)
      tens21=-tens12
      tens22=eps1+omega*dep1do
      tens33=eps3+omega*dep3do
      emod2 = REAL (conjg(exde)*exde*tens11
     .            + conjg(exde)*eyde*rootm1*tens12
     .            + conjg(eyde)*exde*rootm1*tens21
     .            + conjg(eyde)*eyde*tens22
     .            + conjg(ezde)*ezde*tens33)
c      write(*,*)'in grpde2  emod2',emod2
c
      edenfac=0.5*(bmod2+emod2)
c      write(*,*)'in grpde2  edenfac',edenfac

c
c     nper=nperol
      return
c
      end


c       *******************testel*******************
c       test of polarization and fluxn calculations
c----------------------------------------------------
c       input: cnpar,cnpert -refractive index components
c              xe,ye
c       output:ex,ey,ez,flux
c----------------------------------------------------
      subroutine testel(cnpar,cnper,xe,ye,ex,ey,ez,flux)
      implicit integer (i-n), real*8 (a-h,o-z)
c------------------------------------------------------
c			| s , -id,  0|
c     dielectric tensor=| id,	s,  0|
c			| 0 ,	0,  p|
c-----------------------------------------------------
      cn2=cnper*cnper+cnpar*cnpar
      cs=cnpar/dsqrt(cn2)
      cos2=cs*cs
      sn=cnper/dsqrt(cn2)
      sin2=sn*sn
      cn2=cnper*cnper+cnpar*cnpar
      ye2=ye*ye
      s=1.d0-xe/(1.d0-ye2)
      d=-xe*ye/(1.d0-ye2)
      p=1.d0-xe
      write(*,*)'in testel xe,ye,cos2,sin2,cn2'
      write(*,*)xe,ye,cos2,sin2,cn2
      write(*,*)'in testel s,d,p'
      write(*,*)s,d,p
c------------------------------------------------------
      ez=(s-cn2)/d*cn2*sn*cs/(p-cn2*sin2)
      ex=-(s-cn2)/d
      emod2=1.d0+ez*ez+ex*ex
      ey2=1.d0/emod2
c----------let .ge.0
      if(ez.lt.0.d0) then
        ey=-dsqrt(ey2)
      else
        ey=-dsqrt(ey2)
      endif
c-----------------------------
      ex=ex*ey
      ez=ez*ey
      write(*,*)'in testel real_ex imag_ey real_ez'
      write(*,*)ex,ey,ez
c-----------------------------
      dwsdw=s+2.d0*xe/((1.d0-ye2)*(1.d0-ye2))
      dwddw=d+xe*ye*(3.d0-ye2)/((1.d0-ye2)*(1.d0-ye2))
      write(*,*)'1 dwddw',dwddw,'dwsdw',dwsdw
      dwddw=2.d0*xe*ye/((1.d0-ye2)*(1.d0-ye2))
      write(*,*)'2 dwddw',dwddw
      dwpdw=1.d0+xe
c------------------------------------------------------
      fluxel1=(ex*ex+ey*ey)*dwsdw
      fluxel2=ez*ez*dwpdw
      fluxel3=2.d0*ex*ey*dwddw
      write(*,*)'ex,ey,ex*ex+ey*ey',ex,ey,ex*ex+ey*ey
      write(*,*)'fluxel1',fluxel1
      write(*,*)'fluxel21',fluxel2
      write(*,*)'fluxel3',fluxel3
      fluxel=fluxel1+fluxel2+fluxel3
      write(*,*)'fluxel',fluxel
c------------------------------------------------------
      cbx=-cnpar*ey
      cby=-cnper*ez+cnpar*ex
      cbz=cnper*ey
      write(*,*)'cbx,cby,cbz',cbx,cby,cbz
      bmodb=cbx*cbx+cby*cby+cbz*cbz
      write(*,*)'in testel bmodb',bmodb
      flux=bmodb+fluxel
      write(*,*)'in testel flux',flux
      return
      end
      
      
      
      subroutine efKarney(z,r,r0x,z_eff,temp,den,jwave,cnpar,effKarn)
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
c     double precision function: psif(z,r) and
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
c----------------------------------------
c      write(*,*)'in efKarney zd,rd',zd,rd
      psid=psif(zd,rd)
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
      phi=0.d0 !   ????
      bmod=b(zd,rd,phi)
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



      double complex function dcold_rlt(K,aK,nll,np)
c     calculates complex dispersion function dcold_rlt
c     using the cold plasma hermitian tensor K and 
c     relativistic anti-hermition dielectric tensor aK
c     for electron plasma
c     input
c       K(3,3)   complex hermition tensor
c       aK(3,3)  complec anti-hermition relativistic tensor
c       nll      N_parallel
c       np       N_perpendicular
      implicit none
c     input
      double complex K(3,3), aK(3,3)
      double precision nll,np
c     local
      double complex Kxx,Kxy,Kxz,Kyx,Kyy,Kyz,Kzx,Kzy,Kzz
      double precision nlls,nps
      double complex i

c      write(*,*)'dcold_rlt nll,np',nll,np		
c      write(*,*)'dcold_rlt K',K
c      write(*,*)'dcold_rlt aK',aK
      i=dcmplx(0.d0,1.d0)

      Kxx=K(1,1)+aK(1,1)*i
      Kxy=K(1,2)+aK(1,2)*i
      Kxz=K(1,3)+aK(1,3)*i
      
      Kyx=K(2,1)+aK(2,1)*i
      Kyy=K(2,2)+aK(2,2)*i
      Kyz=K(2,3)+aK(2,3)*i
      
      Kzx=K(3,1)+aK(3,1)*i
      Kzy=K(3,2)+aK(3,2)*i
      Kzz=K(3,3)+aK(3,3)*i

      nlls=nll*nll
      nps=np*np

c      write(*,*)'nlls,nps',nlls,nps 
c      write(*,*)'Kxx,Kxy,Kxz',Kxx,Kxy,Kxz 
c      write(*,*)'Kyx,Kyy,Kyz',Kyx,Kyy,Kyz 
c      write(*,*)'Kzx,Kzy,KXz',Kzx,Kzy,Kzz 
      dcold_rlt=(Kxx-nlls) * (Kyy-nlls-nps) * (Kzz-nps)
     .+  Kxy * Kyz * (Kzx+np*nll) 
     .+ (Kxz+np*nll) * Kyx * Kzy
     .- (Kzx+np*nll) * (Kyy-nlls-nps) * (Kxz+np*nll)
     .-  Kzy * Kyz * (Kxx-nlls)
     .- (Kzz - nps) * Kyx * Kxy

c      write(*,*)'dcold_rlt',dcold_rlt

      return
      end

      double precision function dDcold(K,nll,np)
c     calculates the derivative from dD/d(N_perp) the electron cold plasma
c     dispersion function D
c
c     input
c       K(3,3)   complex cold plasma tensor
c       nll      N_parallel
c       np       N_perpendicular
      implicit none
c     input
      double complex K(3,3)
      double precision nll,np
c     local
      double complex Kxx,Kxy,Kxz,Kyx,Kyy,Kyz,Kzx,Kzy,Kzz
      double precision nlls,nps

      nlls=nll*nll
      nps=np*np

c-----cold plasma dispersion function 
c      dcold_rlt=(Kxx-nlls) * (Kyy-nlls-nps) * (Kzz-nps) 
c     .- (Kzx+np*nll) * (Kyy-nlls-nps) * (Kxz+np*nll)
c     .- (Kzz - nps) * Kyx * Kxy
c     .+  Kxy * Kyz * (Kzx+np*nll)   !YuP[11-2016] added this line
c     .+ (Kxz+np*nll) * Kyx * Kzy    !YuP[11-2016] added this line


      Kxx=K(1,1)
      Kxy=K(1,2)
      Kxz=K(1,3)
      !
      Kyx=K(2,1)
      Kyy=K(2,2)
      Kyz=K(2,3) !dDcold()  YuP[11-2016] : was K(3,3)
                 !however, Kyz was not used (now added, see below)
      Kzx=K(3,1)
      Kzy=K(3,2)
      Kzz=K(3,3)
 
      dDcold=(Kxx-nlls) * (-2.d0*np) * ((Kzz-nps)+(Kyy-nlls-nps)) 
     .- (nll) * (Kyy-nlls-nps) * (Kxz+np*nll)
     .- (Kzx+np*nll) * (-2.d0*np) * (Kxz+np*nll)
     .- (Kzx+np*nll) * (Kyy-nlls-nps) * (nll)
     .- (-2.d0*np) * Kyx * Kxy
     ++ (Kxy*Kyz+Kyx*Kzy)*nll ! YuP[11-2016] added this line: zero?
      !usually Kyz=0 and Kzy=0
      !write(*,*)'dDcold',dDcold,Kxy,Kyz,Kyx,Kzy
      return
      end


      subroutine check_monotonic_radius_in_ray_element(is,
     &i_non_monotonic,z_center,r_center,rho_center)
c YuP: Not needed? All calls are commented out.
c-----------------------------------------------------------
c     check if the small radius is a monotonic function along 
c     the ray element [M_(is-1),M_is]
c     Here: M_is={z(is),r(is),phi(is)}
c-----------------------------------------------------------
      implicit none
      include 'param.i'
      include 'write.i'
c-----input
      integer is  !  >1 the number of ray element
c-----output
      integer i_non_monotonic   ! =0 the maximal and minimal rho are
                                !    at the boundaries of the ray element  
                                ! =1 the maximal or minimal rho are
                                !    inside the ray element
      real*8  z_center,r_center ! The coordinates of the extremum rho
                                ! at  i_non_monotonic=1 case
      real*8 rho_center         !extremal rho at i_non_monotonic=1 case
c-----locals
      real*8 psi_left,psi_right,rho_left,rho_right,h_z,h_r,
     &z_i,r_i,psi_i,rho_i,rho_min,rho_max

      integer n_steps,i_max,i_min,i
c-----externals
      real*8 psif, rhopsi

      psi_left=psif(wz(is-1),wr(is-1))
      psi_right=psif(wz(is),wr(is))

      rho_left=rhopsi(psi_left)
      rho_right=rhopsi(psi_right)

      n_steps=100
      h_z=(wz(is)-wz(is-1))/(n_steps-1)
      h_r=(wr(is)-wr(is-1))/(n_steps-1)

      rho_min = rho_left
      rho_max = rho_left
      i_max=0
      i_min=0
      do i=1,n_steps
        z_i=wz(is-1)+h_z*i
        r_i=wr(is-1)+h_r*i
        psi_i=psif(z_i,r_i)
        rho_i=rhopsi(psi_i)
       
        if (rho_i.gt.rho_max) then
           i_max=i
           rho_max=rho_i
        endif

        if (rho_i.lt.rho_max) then
           i_min=i
           rho_min=rho_i
        endif
      enddo

      if((i_min.gt.1).and.(i_min.lt.n_steps)) then
c--------the minimal small radius value is inside the ray element
         i_non_monotonic=1
         rho_center=rho_min
         z_center=wz(is-1)+h_z*i_min
         r_center=wr(is-1)+h_r*i_min
         goto 10
      endif

      if((i_max.gt.1).and.(i_max.lt.n_steps)) then
c--------the maximal small radius value is inside the ray element
         i_non_monotonic=1
         rho_center=rho_max
         z_center=wz(is-1)+h_z*i_max
         r_center=wr(is-1)+h_r*i_max
         goto 10
      endif

      i_non_monotonic=0

 10   continue

      return
      end     


c     this  subroutine calculates current drive efficiency
c     for the toroidal  plasma using TorGA_curgap codo written by Lin-Liu
c     ATTENTION:all parameters and variables input to curba
c       are real(no double precision)
c-----------------------------------------------------------------
      subroutine eff_Lin_Liu(z,r,zma,rma,r0x,z_eff,temp,den,jwave,cnpar,
     &                   cnper,ioxm,ye,
     &                   cefldx,cefldy,cefldz,
     &                   efficient)

c     input parameters: z,r -coordinates of the ray point(cm)

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
c        psif(z,r) and  subroutine: zr_psith(psi,theta,z,r)
c-------------------------------------------------------------------
      implicit none
c-----input
      real*8 z,r,zma,rma,z_eff,temp,den,cnpar,ye,cnper
      real*8 r0x
      complex*16 cefldx,cefldy,cefldz
      integer jwave,ioxm
c-----output
      real*8 efficient

c-----locals
      real*8 psid,zd,rd,rmaxpsid,rminpsid,pid,
     &zbmin,zbmax,zb
      real*8 tc,tol,elomom,prho,ctheta,stheta,pi,theta,thtc,enpar,denom,
     *aspct,zacosarg,rjpd,rjpd0,ratjpd,efficien
      integer ig,n0,model,lh
c-----external 
      real*8 psif,bmin_psi,bmax_psi,rmax_psi,rmin_psi,b,zfac_f

      
c-------------------------------------------------------------------
c     for TorGA_curgap
      real*8 rjpd_d,rjpd0_d,ratjpd_d,denom_d,aspct_d,enpar_d, 
     &tc_d, thtc_d, theta_d, elomom_d,zeff_d,tol_d, omode
c-----external
c-----for efficiency transformation like in currn in toray
      real*8 phi,zfac
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
c      write(*,*)'z,r',z,r
      zd=dble(z)*0.01d0/r0x
      rd=dble(r)*0.01d0/r0x
      psid=psif(zd,rd)
c      write(*,*)'zd,rd,psid',zd,rd,psid
cSm040426
cc     rmax_psid and rmin_psid are the largest and and the least
cc     values of the radius (r) on the given flux surface
cc     psid=psid(zd,rd)
c------------------------------
c      call zr_psith(psid,0.d0,zd,rmaxpsid)
c      call zr_psith(psid,pid,zd,rminpsid)
      
c      if(abs(rmaxpsi-rminpsi).lt.1.0e-8) rminpsi=rmaxpsi-1.0e-8
c      aspct=(rmaxpsi-rminpsi)/(rmaxpsi+rminpsi)
      
c----------------------------------
c     Here :aspct is calculating using double precision functions
c           rmax_psi(psid) and rmin_psi(psid)
c----------------------------------------
c     rmaxpsi and rminpsi are the largest and and the least
c     values of the radius (rmaxpsi and rminpsi in (m))
c     on the given flux surface psi=psi(zd,rd)
c--------------------------------------------
c      rmaxpsid=rmax_psi(psid)
c      rminpsid=rmin_psi(psid)
c      if (rmaxpsid.lt.rd)rmaxpsid=rd
c      if (rminpsid.gt.rd)rminpsid=rd
c      if(dabs(rmaxpsid-rminpsid).lt.1.0d-8) rminpsid=rmaxpsid-1.0d-8
c      aspct=(rmaxpsid-rminpsid)/(rmaxpsid+rminpsid)
cSm040426
       zbmin=bmin_psi(psid)    
       zbmax=bmax_psi(psid)
       aspct=((zbmax-zbmin)/(zbmax+zbmin))
cSm060201
       if (aspct.lt.1.d-8) aspct=1.d-8 

c       write(*,*)'psid,zbmax,zbmin,aspct',psid,zbmax,zbmin,aspct
c      write(*,*)'aspct from functionsrmin_psi, rmax_psi',aspct
c      write(*,*)'aspct',aspct,'enpar',enpar,'tc',tc,'thtc',thtc
c      write(*,*)'theta',theta,'elomom',elomom,'lh',lh
c      write(*,*)'z_eff',z_eff,'model',model,'n0',n0,'ig',ig

cSm040503
c     function b() uses the arguments z and r in [m]
      zb=b(.01d0*dble(z)/r0x,.01d0*dble(r)/r0x,0.d0)!it changes 
                                         !bz,br,bphi and rho in one.i
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
      zfac=zfac_f(dble(z*1.d-2),dble(r*1.d-2),0.d0,dble(temp))
      zfac=zfac*1.d-7 ! watt => egr/sec
      zfac=zfac*1.d-4 ! m**2 => cm**2
c      write(*,*)'zfac',zfac
      efficient=efficien*zfac
cSm080629
c      efficient=-efficient 

      return
      end



       subroutine etajrf(z,eps,w,rto,eta,eta0,irfq) 
! +-------------------------------------------------------------------+
! | Last modified: July 3,1997                                        |
! | Written by P. Bonoli & J. Wright                                  |
! +-------------------------------------------------------------------+
! |  Evaluates the normalized current-drive efficiency using          |
! |  the parametrization of Ehst and Karney                           |
! |  (Argonne National Laboratory Report ANL/FPP/TM-247 (1990)).      |
! |  IRFQ = 1 indicates alfven wave type damping                      |
! |  IRFQ = 2 indicates landau wave type damping                      |
! +-------------------------------------------------------------------+
!   w=wte
!   rto=bmod/bmax
!   irfq=jwave+2 = For Alfven +1
!                = for Ladaw damp +2
      implicit none

      integer, intent(in) :: irfq

      real*8, intent(in) :: z, eps, w, rto
      real*8, intent(out) :: eta, eta0
      real*8 :: xr, en, ek, a, fk, fc, ec, em, r, flamt, fm,
     &            yt, ytt, gc, rtp

c      write(*,*)'in etajrf z,eps,w,rto,irfq',z,eps,w,rto,irfq 

      xr = 3.5d0
      en = 0.77d0
      ek = 3.d0

      if(irfq.eq.1)  then
         a = 12.3d0
         fk = 11.91d0/(0.678d0 + z)
         fc = 4.13d0/z**0.707d0
         ec = 0.0987d0
         em = 2.48d0
      endif

      if(irfq .eq. 2) then
         a = 0.d0
         fk = 3.d0/z
         fc = 3.83d0/z**0.707d0
         ec = 0.389d0
         em = 1.38d0
      endif

c      write(*,*)'in etajrf a,fk,fc,ec ,em',a,fk,fc,ec ,em

      eta0 = (fk/w + fc + 4.d0*w*w/(5.d0 + z))
      
c      write(*,*)'in etajrf eta0',eta0
     
      r = 1.d0 - eps**en*sqrt(xr*xr + w*w)/(eps**en*xr + w)
c      R=r=rprof

       rtp = min(rto,0.9999999d0)
       flamt = sqrt(1.d0 - rtp)
c------flamt=alt1=sqrt(1-bmod/bmx) => rtp=bmod/bmax 
       fm = 1.d0 + a*(flamt/w)**ek

       yt = (1.d0 - flamt*flamt)/(flamt*flamt)*w*w
c-----      Y_t=yt=ytt
c-----      flamt=alt1
       ytt = (ec*yt)**em
       ytt = min(500.d0,ytt)
       gc = 1.d0 - exp(-ytt)
       eta = gc*fm*eta0*r !normalized eta

       return
       end subroutine etajrf


       Subroutine efKarney_Bonoli(z,r,r0x,z_eff,temp,den,jwave,cnpar,
     & effKarn)
c-------RF current drive efficiency(asymptotic formula, nonrelativistic case)  
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
c
c It uses Bonoli subroutine etajrf
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
c     double precision function: psif(z,r) and
c     subroutine: zr_psith(psi,theta,z,r)
c     double precision functions from zr_psith.f: rmax_psi(psi),
c     rmin_psi(psi), bmax_psi(psi)

      implicit none
c-----input
      real*8 z,r, !cm
     &r0x,        !normalization lemgth
     &z_eff,      !effective charge
     &temp,       !temperatura KeV
     &den,        !density 10**13
     &cnpar       !parallel to magnetic field refractive index

      integer jwave 
c-----output
      real*8 ffKarn    !efficiency  in (A/cm**2)/(erg/(c*cm**3))


c-----locals
      real*8
     &w,           !  w=wte=u_res=v_par/v_e=cvac/(ve*cnpar)
     &cvac,        !  light speed [cm/sec]
     &ve,          !  ve=sqrt(T_e/m_e)
     &zd,rd,       !  coordinates in [m]
     &phi,psid,rmaxpsi,rminpsi,aspct,eps,bmod,bmax,
     &rto,         !  b/bmax
     &eta0,eta,    !  normalized efficiency from Bonoli subrourine
     &arg1,cln,s
      integer irfq !=jwave+2
       
c-----output
      real*8 effKarn !efficiency  in (A/cm**2)/(erg/(c*cm**3)

c-----externals
      real*8 psif,b,bmax_psi,rmin_psi,rmax_psi

      cvac=2.997930d+10
      ve=1.32d9*dsqrt(temp)	 
      w=cvac/(ve*cnpar)
      w=dabs(w)     !w=|wte|
      write(*,*)'Bon w',w
c-------------------------------------------
c     conversion from cm to m
      zd=z*0.01d0/r0x
      rd=r*0.01d0/r0x
      phi=0.d0       !arbitrary vatoroidal angle
      write(*,*)'Bon z,r,cnpar',z,r,cnpar 
c----------------------------------------------
c     calculate epsilon
c-------------------------------------------
      psid=psif(zd,rd) !poloidal flux at z,r
      rmaxpsi=rmax_psi(psid)
      rminpsi=rmin_psi(psid)
      if (rmaxpsi.lt.rd)rmaxpsi=rd
      if (rminpsi.gt.rd)rminpsi=rd 
      if(dabs(rmaxpsi-rminpsi).lt.1.0d-8) rminpsi=rmaxpsi-1.0d-8
      aspct=(rmaxpsi-rminpsi)/(rmaxpsi+rminpsi)
      eps=aspct 
      write(*,*)'Bon eps',eps
c------------------------------------------------
c     calculate rto=b/bmax
c--------------------------------------------------
      bmod=b(zd,rd,phi)
      bmax=bmax_psi(psid)
      if (bmod.gt.bmax) bmod=bmax 
      rto=bmod/bmax
      write(*,*)'Bon rto',rto
c--------------------------------------------------
      irfq=jwave+2
      write(*,*)'Bon irfq',irfq
c----------------------------------------------------
c     calculate normalized efficiency eta 
c---------------------------------------------------
      write(*,*)'Bon before etajrf z_eff',z_eff,eps,w,rto,irfq
      call etajrf(z_eff,eps,w,rto,eta,eta0,irfq)
      write(*,*)'Bon after etajr rta0,eta',eta0,eta
c-------------------------------------------------- 
c     efficiency  in (A/m**2)/(joule/(c*m**3))
c     temperature temp in kev
c     density     dens in 10**19 /m**3
c-----------------------------------------------------------
c     cln -Coulomb Ln
      arg1 = 1.d3/temp*dsqrt(10.d0*den)	!temp KeV, den 10**13/cm**3
      cln=24.d0-dlog(arg1)

      effKarn=eta*3.84d0*temp/(cln*den)	!(a/m**2)/joule/(c*m**3)
      write(*,*)'Bon temp,den,cln,eta,effKarn',temp,den,cln,eta,effKarn
c--------------------------------------------------------------  
c     determination the of the current drive sign
c------------------------------------------------------
      if (cnpar.ge.0.d0) then
         s=1.d0
      else
         s=-1.d0
      endif
c------------------------------------------------------------
      effKarn=-effKarn*s*1.d-5  !  A/cm**2/(erg/(sec*cm**3)

      return
      end

      subroutine plot_1_Karney
c----------------------------------------------------
c     creates data (arrays) for plot like Fig. 1
c     at Karney article
c     Nuclear Fusion 1991 p. 1934
c--------------------------------------------------- 
      implicit none
      include 'param.i'
      include 'three.i'

      integer n_points_w,  !number of plot points in w direction
     &n_epsilon,  !number of epsilons
     &n_theta_pol !number of polidal angles
      parameter (n_points_w=11)
c      parameter (n_points_w=3)
      parameter (n_epsilon=3)
      parameter (n_theta_pol=4)

      integer i,j,k

      real*8 
     &epsilon,        !inverse aspect ratio
     &epsilon_ar(n_epsilon),
     &theta_pol,      !poloidal angle [radians]
     &theta_pol_ar(n_theta_pol),
     &w_ar(n_points_w),w,    !=clight/(N_parallel*V_te)
                      !here m_e*V_te**2=T_e
     &w_min,w_max,step,
     &pi,
     &psi,            !poloidal flux(epsilon)
     &z_eff,
     &eta,eta0,       !CD efficiency
     &eta_ar_fw(n_points_w,n_theta_pol,n_epsilon),
     &eta0_ar_fw(n_points_w,n_theta_pol,n_epsilon),
     &eta_ar_lh(n_points_w,n_theta_pol,n_epsilon),
     &eta0_ar_lh(n_points_w,n_theta_pol,n_epsilon),
     &z,r,phi,        !space coordinates
     &rto,            !rto=b(z,r,phi)/bmax
     &bmod,bmax,
     &accuracy       !accuracy to solve the equation
                     ! epsilon_psi(psi)=epsilon
      real*8 rmaxpsi,rminpsi,rmax_psi,rmin_psi

      integer irfq !=1 indicates alfven wave type damping                      
                    != 2 indicates landau wave type damping       

      character*16 file_nm ! name of output netcdf file
                   
c-----externals
      real*8 psi_epsilon,bmax_psi,b,rhopsi


      pi=4.d0*datan(1.d0)

      theta_pol_ar(1)=0.d0
      theta_pol_ar(2)=pi/2.d0
      theta_pol_ar(3)=3.d0*pi/4.d0
      theta_pol_ar(4)=pi

      epsilon_ar(1)=0.d0
      epsilon_ar(2)=0.03
      epsilon_ar(3)=0.1d0

      z_eff=1.d0
      irfq=1
      phi=0.d0
      w_min=0.1d0
      w_max=10.d0
c      step=(w_max-w_min)/(n_points_w-1)
      step=(dlog(w_max/w_min)/dlog(10.d0))/(n_points_w-1 )    
      write(*,*)'in prep3d subroutine plot_1_Karney'
      accuracy=1.d-7
      do i=1,n_epsilon
         epsilon=0.d0
         epsilon=epsilon_ar(i)
         write(*,*)'i,epsilon',i,epsilon

         psi=psi_epsilon(epsilon,psimag,psilim, accuracy)
         write(*,*)'after psi_epsilon psi=',psi
         write(*,*)'rhopsi(psi)',rhopsi(psi)
c--------check epsilon
         rmaxpsi=rmax_psi(psi)
         rminpsi=rmin_psi(psi)
     
         if(dabs(rmaxpsi-rminpsi).lt.1.0d-8) rminpsi=rmaxpsi-1.0d-8
         epsilon=(rmaxpsi-rminpsi)/(rmaxpsi+rminpsi)
         write(*,*)'check epsilon',epsilon

         do j=1,n_theta_pol
            theta_pol=theta_pol_ar(j)
            call zr_psith(psi,theta_pol,z,r)
            bmod=b(z,r,phi)
            bmax=bmax_psi(psi)
            if (bmod.gt.bmax) bmod=bmax 
            rto=bmod/bmax
            write(*,*)'i,j,rto',i,j,rto
            do k=1,n_points_w
               w_ar(k)=w_min+step*(k-1)
c               x=dlog(w_min)/dlog(10.d0)+step*(k-1)
               w_ar(k)=dexp(dlog(10.d0)*
     &                      (dlog(w_min)/dlog(10.d0)+step*(k-1)))
               w=w_ar(k)
               write(*,*)'k, w_min,w',k, w_min,w
               irfq=1 !fw
               call etajrf(z_eff,epsilon,w,rto,eta,eta0,irfq)
               eta_ar_fw(k,j,i)=eta 
               eta0_ar_fw(k,j,i)=eta0 
               write(*,*)'fw eta0,eta',eta0,eta
               irfq=2 !lh
               call etajrf(z_eff,epsilon,w,rto,eta,eta0,irfq)
               eta_ar_lh(k,j,i)=eta 
               eta0_ar_lh(k,j,i)=eta0
               write(*,*)'lh eta0,eta',eta0,eta
            enddo
        enddo
      enddo

c--------------------------------------------------------------
c     write CD efficiency eta and w for fig. 1 to netcdf file
c------------------------------------------------------------
      file_nm='Ehst-Karney_plot' ! set name of output netcdf file

      write(*,*)'before netcdf_Karne !sqrt(T/m)y file_nm=',file_nm

cSAP081110
c      call netcdf_Karney(file_nm,n_points_w,w_ar,
c     &n_theta_pol,theta_pol_ar,n_epsilon,epsilon_ar,
c     &eta0_ar_fw,eta_ar_fw,eta0_ar_lh,eta_ar_lh)

      write(*,*)'after netcdf_Karney '

      return
      end

      real*8 function epsilon_psi(psi)
c-------------------------------------------------------------
c     calculates inverse aspect ratio 
c     epsilon=(rmaxpsi-rminpsi)/(rmaxpsi+rminpsi)
c     at given ploidal flux psi
c---------------------------------------------------------------
      implicit none
     
c-----input
      real*8 psi ! poloidal flux

c-----externals rmax_psi,rmin_psi
      real*8 rmax_psi,rmin_psi


c-----locals
      real*8 rmaxpsi,rminpsi
      
      stop 'epsilon_psi Not setup for xyz?'

      rmaxpsi=rmax_psi(psi)
      rminpsi=rmin_psi(psi)

      if(dabs(rmaxpsi-rminpsi).lt.1.0d-8) rminpsi=rmaxpsi-1.0d-8

      epsilon_psi=(rmaxpsi-rminpsi)/(rmaxpsi+rminpsi)

      return
      end

      real*8 function epsilon_psi_minus_epsilon(psi)
c------------------------------------------------------------------
c     Calculates (epsilon_psi(psi)-epsilon)
c     It used to find the root of the equation
c     epsilon_psi(psi)-epsilon=0
c-----------------------------------------------------------------
      implicit none
c-----input
      real*8 psi !poloidal flux

      real*8 epsilon_in
      common /epsilon_psi_minos_epsilon/ epsilon_in     
c-----externals
      real*8 epsilon_psi

      epsilon_psi_minus_epsilon=epsilon_psi(psi)-epsilon_in
       
      return
      end

      real*8 function psi_epsilon(epsilon,psimag,psilim, accuracy)
c-----------------------------------------------------------
c     calculates poloidal flux at given inverse aspect ratio epsioln
c     epsilon=(rmaxpsi-rminpsi)/(rmaxpsi+rminpsi)
c
c     It uses bisection method at psimag < psi < psilim
c-------------------------------------------------------------
      implicit none
     
c-----input
      real*8 epsilon, !inverse aspect ratio
     &psimag,         !polidal flux at the magnetic axis 
     &psilim,         !poloidal flux at the last closed flux surface
     &accuracy        !bisection method accuracy 
                      !It uses the condition |difference of root| < accuracy
      real*8 epsilon_in
      common /epsilon_psi_minos_epsilon/ epsilon_in ! will be set equal to
                                                    ! the input epsilon
c-----externals
      real*8  epsilon_psi_minus_epsilon,rtbis
      external  epsilon_psi_minus_epsilon

c-----initialization of the common block
      epsilon_in =epsilon

      if(epsilon.lt.1.d-10) then
        psi_epsilon=psimag
      else

        write(*,*)'in psi_epsilon psimag,psilim,accuracy',
     &                            psimag,psilim,accuracy

        psi_epsilon=RTBIS(epsilon_psi_minus_epsilon,
     &psimag,psilim,accuracy)

      endif

      return
      end
     
      subroutine netcdf_Karney(file_nm,n_points_w,w_ar,
     &n_theta_pol,theta_pol_ar,n_epsilon,epsilon_ar,
     &eta0_ar_fw,eta_ar_fw,eta0_ar_lh,eta_ar_lh)
c-------------------------------------------------------------
c     writes data for fig1 Ehst-Karney  to the netcdf file
c     file_nm.nc
c--------------------------------------------------------------
      implicit none
      include 'netcdf.inc'
c-----input
      character*(*) file_nm        !name of output nc file
      integer n_points_w,          !number of points
     &n_theta_pol,                 !number of poloidal angle values
     &n_epsilon                    !number of epsilon values
      real*8,  dimension(1:n_points_w)  :: w_ar         ! c/(N_par*v_te)
      real*8,  dimension(1:n_theta_pol) :: theta_pol_ar !ploidal angle [radian]
      real*8,  dimension(1:n_epsilon)   :: epsilon_ar   !inverse aspect ratio
      real*8,  dimension(1:n_points_w, 1:n_theta_pol, 1:n_epsilon) 
     &                                  :: eta_ar_fw,eta0_ar_fw, !FW efficiency
     &                                     eta_ar_lh,eta0_ar_lh  !LH efficiency
c-----locals----------------------------------------  
      integer i,j,k
c --- some stuff for netCDF file ---
      integer ncid,vid,istatus,n_points_w_id,n_theta_pol_id,
     &n_epsilon_id,
     &start(3),starts(3),eta_dims(3),eta_count(3),
     +nccre2,ncvdef2,ncdid2,ncddef2

      character ltitle*512,filenc*128

c-----Storage tem1 is used in netcdf write.   
      real*8, dimension (1: n_points_w*n_theta_pol) :: tem1
      real*8, dimension (1: n_points_w,1:n_theta_pol) :: tem2
c-----externals
      integer length_char

      data start/1,1,1/

      write(*,*)'in netcdf_Karney begin'

      eta_count(1)=n_points_w
      eta_count(2)=n_theta_pol
      eta_count(3)=1

      write(*,*)'in netcdf_Karney n_points_w,n_theta_pol,n_epsilon',
     &                            n_points_w,n_theta_pol,n_epsilon
      write(*,*)'in netcdf_Karney file_nm=',file_nm

      if( length_char(file_nm).gt.124)
     &   stop 'Adjust file_nm in netcdf_Karney'

      write(filenc,1002) file_nm
 1002 format(a,".nc")
      write(*,*)'in netcdf_Karney after length_char(file_nm)'
c-------------------------------------------------------------
c      create net CDF file  define dimensions,variables
c          and attributes
c------------------------------------------------------------
      ncid=nccre2(filenc,NCCLOB,istatus)
      call check_err(istatus)

c     Brief description added to file:
      ltitle='netCDF file of data for fig.1 WEhst-Karney article'
      if( length_char(ltitle).gt.512 ) 
     &   stop 'Adjust ltitle in netcdf_Karney'
      call ncaptc2(ncid,NCGLOBAL,'title',NCCHAR,length_char(ltitle),
     &     ltitle,istatus)

c.......................................................................
cl    1.1.2 define dimensions
c.......................................................................
      n_epsilon_id=ncddef2(ncid,'n_epsilon',n_epsilon,istatus)
      n_theta_pol_id=ncddef2(ncid,'n_theta_pol',n_theta_pol,istatus)
      n_points_w_id=ncddef2(ncid,'n_points_w',n_points_w,istatus)
     
c.......................................................................
cl    1.1.3 define variables
c.......................................................................
      vid=ncvdef2(ncid,'n_points_w',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,27,
     +            'Number of argument points_w',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'n_epsilon',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,18,
     +            'Number of epsilons',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'n_theta_pol',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25,
     +            'Number of poloidal angles',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'w_ar',NCDOUBLE,1,n_points_w_id,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,41,
     +     'array of phase velocities w=c/(N_par*v_e)',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'non-dimensional',istatus)
      call check_err(istatus)


      vid=ncvdef2(ncid,'theta_pol_ar',NCDOUBLE,1,n_theta_pol_id,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,24,
     +     'array of poloidal angles',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,
     +           'radians',istatus)
      call check_err(istatus)  

      vid=ncvdef2(ncid,'epsilon_ar',NCDOUBLE,1,n_epsilon_id,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,16,
     +     'array of epsilon',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,15,
     +           'non-dimensional',istatus)
      call check_err(istatus)

      eta_dims(1)=n_points_w_id
      eta_dims(2)=n_theta_pol_id
      eta_dims(3)=n_epsilon_id

      write(*,*)'ncid=',ncid
      write(*,*)'eta_dims',eta_dims
      vid=ncvdef2(ncid,'eta_ar_fw',NCDOUBLE,3,eta_dims,istatus)
      write(*,*)'istatus',istatus
      call ncaptc2(ncid,vid,'long_name',NCCHAR,20,
     +           'FW CD efficiency eta',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,10,
     +           'normalized',istatus)

      write(*,*)'ncid=',ncid
      write(*,*)'eta_dims',eta_dims
      vid=ncvdef2(ncid,'eta0_ar_fw',NCDOUBLE,3,eta_dims,istatus)
      write(*,*)'istatus',istatus
      call ncaptc2(ncid,vid,'long_name',NCCHAR,21,
     +           'FW CD efficiency eta0',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,10,
     +           'normalized',istatus)

      write(*,*)'ncid=',ncid
      write(*,*)'eta_dims',eta_dims
      vid=ncvdef2(ncid,'eta_ar_lh',NCDOUBLE,3,eta_dims,istatus)
      write(*,*)'istatus',istatus
      call ncaptc2(ncid,vid,'long_name',NCCHAR,20,
     +           'LH CD efficiency eta',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,10,
     +           'normalized',istatus)

      write(*,*)'ncid=',ncid
      write(*,*)'eta_dims',eta_dims
      vid=ncvdef2(ncid,'eta0_ar_lh',NCDOUBLE,3,eta_dims,istatus)
      write(*,*)'istatus',istatus
      call ncaptc2(ncid,vid,'long_name',NCCHAR,21,
     +           'LH CD efficiency eta0',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,10,
     +           'normalized',istatus)

c------------------------------------------------------------
cl    1.1.4 end the define-mode and start the data-mode
c     p. 51-2 of manual
c----------------------------------------------------------
      call ncendf3(ncid,istatus)
      call check_err(istatus)

      write(*,*)'end initialization'
c----------------------------------------------------------------
c     write data
c-----------------------------------------------------------------
      write(*,*)'before vid=ncvid(ncid,n_epsilon'
      call ncvid2(vid,ncid,'n_epsilon',istatus)
      call ncvpt_int2(ncid,vid,0,0,n_epsilon,istatus)

      call ncvid2(vid,ncid,'n_theta_pol',istatus)
      call ncvpt_int2(ncid,vid,0,0,n_theta_pol,istatus)

      call ncvid2(vid,ncid,'n_points_w',istatus)
      call ncvpt_int2(ncid,vid,0,0,n_points_w,istatus)

      call ncvid2(vid,ncid,'epsilon_ar',istatus)  
      call ncvpt_doubl2(ncid,vid,1,n_epsilon,epsilon_ar,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'theta_pol_ar',istatus)  
      call ncvpt_doubl2(ncid,vid,1,n_theta_pol,theta_pol_ar,istatus)
      call check_err(istatus)

      call ncvid2(vid,ncid,'w_ar',istatus)  
      call ncvpt_doubl2(ncid,vid,1,n_points_w,w_ar,istatus)
      call check_err(istatus)
 
c-----FW    
      do i=1,n_epsilon
         do j=1,n_theta_pol
            do k=1,n_points_w
               tem2(k,j)=eta_ar_fw(k,j,i) 
            enddo
         enddo
          
         call pack21(tem2,1,n_points_w,1,n_theta_pol,
     &               tem1,n_points_w,n_theta_pol)

         starts(1)=start(1)
         starts(2)=start(2)
         starts(3)=i  ! number of epsilon put in .nc file

         call ncvid2(vid,ncid,'eta_ar_fw',istatus)
         call ncvpt_doubl2(ncid,vid,starts,eta_count,tem1,istatus)
      enddo

      do i=1,n_epsilon
         do j=1,n_theta_pol
            do k=1,n_points_w
               tem2(k,j)=eta0_ar_fw(k,j,i)  
            enddo
         enddo
          
         call pack21(tem2,1,n_points_w,1,n_theta_pol,
     &               tem1,n_points_w,n_theta_pol)

         starts(1)=start(1)
         starts(2)=start(2)
         starts(3)=i  ! number of epsilon put in .nc file

         call ncvid2(vid,ncid,'eta0_ar_fw',istatus)
         call ncvpt_doubl2(ncid,vid,starts,eta_count,tem1,istatus)
      enddo

c-----LH
      do i=1,n_epsilon
         do j=1,n_theta_pol
            do k=1,n_points_w             
                tem2(k,j)=eta_ar_lh(k,j,i)
            enddo
         enddo
          
         call pack21(tem2,1,n_points_w,1,n_theta_pol,
     &               tem1,n_points_w,n_theta_pol)

         starts(1)=start(1)
         starts(2)=start(2)
         starts(3)=i  ! number of epsilon put in .nc file

         call ncvid2(vid,ncid,'eta_ar_lh',istatus)
         call ncvpt_doubl2(ncid,vid,starts,eta_count,tem1,istatus)
      enddo

      do i=1,n_epsilon
         do j=1,n_theta_pol
            do k=1,n_points_w
               tem2(k,j)=eta0_ar_lh(k,j,i)  
            enddo
         enddo
          
         call pack21(tem2,1,n_points_w,1,n_theta_pol,
     &               tem1,n_points_w,n_theta_pol)

         starts(1)=start(1)
         starts(2)=start(2)
         starts(3)=i  ! number of epsilon put in .nc file

         call ncvid2(vid,ncid,'eta0_ar_lh',istatus)
         call ncvpt_doubl2(ncid,vid,starts,eta_count,tem1,istatus)
      enddo
C-----------------------------------------------------------------------
c
cl    4. Close netCDF file
c        Page 53-54 netCDF-2 Manual
c
      call ncclos3(ncid,istatus)
      call check_err(istatus)
      return
      end

      

 
