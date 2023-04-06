*****************************************************
c----------------subroutine MK_GRAP----------------*
c Prepares for file drawgenr.in for xdraw              *
c INPUT: from common blocks 'gr' and 'write'        *
c        iray -number of the ray at antenna         *
c*****************************************************
      subroutine MK_GRAP(iray)
      implicit integer (i-n), real*8 (a-h,o-z)
      INTEGER I,J
      include 'param.i'
      include 'gr.i'
      include 'write.i'

      write(*,*)'mk_grap iray=',iray
      if(iray.eq.1) then
c         CALL ASSIGN("assign -F f77 -N ieee u:11",ier)
         OPEN(11,file='genray.bin',form='unformatted')
      end if
      write(*,*)'mk_grap nrayelt=',nrayelt
      if (nrayelt.eq.0)then
         goto 70
      endif
      if (iray.gt.1) then
c        data for noncentral antenna rays
         goto 60
      endif

c-------------------------------------------------------------------
c     data for the central antenna ray
c-------------------------------------------------------------------
      DO 10 I=1,nrayelt
       WRITE(11) REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I))
10    CONTINUE
      WRITE(11)

      DO 30 J=1,3
       DO 20 I=1,NP+1
        WRITE(11) REAL(AR(J,I)),REAL(AZ(J,I)),REAL(XT(J,I)),
     1   REAL(YT(J,I)),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt))
20     CONTINUE
       WRITE(11)
30    CONTINUE

      DO 50 J=4,NL
       DO 40 I=1,NP+1
        WRITE(11) REAL(AR(J,I)),REAL(AZ(J,I)),
     1   REAL(XT(3,NP+1)),REAL(YT(3,NP+1)),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt))
40     CONTINUE
       WRITE(11)
50    CONTINUE
      goto 70
c  end data for the central ray
c---------------------------------------------------------------------

60    continue
c--------------------------------------------------------------
c        data for noncentral antenna rays
c--------------------------------------------------------------
c      p=0.0
      DO 61 I=1,nrayelt
       WRITE(11) REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
c     1  REAL(p),REAL(p),REAL(p),REAL(p),
c     2  REAL(p),REAL(p),REAL(p)
     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I))
61    CONTINUE
      WRITE(11)
c  end data for the noncentral rays
c---------------------------------------------------------------------

70    continue
c      CLOSE(11)
      RETURN
      END


c----------------subroutine MK_GRAPH----------------
c Prepares for file drawgenr.in for xdraw              
c INPUT: from common blocks 'gr' and 'write'        
c        iray -number of the ray at antenna        
c        n_wall is a number of points in arrays
c                which set the wall form
c*****************************************************
cSAP090313
c      subroutine MK_GRAPH(iray,nray,ifreq,nfreq,nbulk)
      subroutine MK_GRAPH(iray,nray,ifreq)
      implicit integer (i-n), real*8 (a-h,o-z)
      INTEGER I,J

c-----input
      integer iray,nray,ifreq

      integer isave

cSm061209
      integer isave_freq

      save isave

cSm061209
      save isave_freq

      include 'param.i'
      include 'gr.i'
      include 'write.i'
cSAP090313
      include 'one.i'
      include 'fourb.i'
      include 'three.i'

      data isave /0/
      data isave_freq/0/

      if (isave.lt.1) then
        isave=iray
      endif

cSm061229
      if(isave_freq.lt.1) then
        isave_freq=ifreq
      endif

cyup      write(*,*)'!!!!!!in mk_graph iray,nrayelt,isave,isave_freq',
cyup     + iray,nrayelt,isave,isave_freq

       if((iray.eq.1).or.(isave.eq.iray)) then
c         CALL ASSIGN("assign -F f77 -N ieee u:82",ier)
         write(*,*)'mk_grph before open 82 and 81'
         OPEN(82,file='genray.bin',form='unformatted')
         open(81,file='genray.doc')
         open(85,file='absorp.bin',form='unformatted')
         open(86,file='absorp.doc')
         open(75,file='efield.bin',form='unformatted')
         open(76,file='eps_r.bin',form='unformatted')
         open(77,file='eps_i.bin',form='unformatted')
         open(78,file='em_res.bin',form='unformatted')
         open(79,file='cd.bin',form='unformatted')
      end if
cSm030508
c      if (nrayelt.eq.0)then
      write(*,*)'nrayelt,iray,isave,ifreq',nrayelt,iray,isave,ifreq
      if (nrayelt.lt.2)then
        if((iray.gt.isave).or.(ifreq.gt.1)) then
          write(*,*)'before goto 70'
          goto 70  
        endif
        
      endif
     
c      write(*,*)'iray,isave,ifreq,isave_freq',
c     &            iray,isave,ifreq,isave_freq

      if (((iray.gt.isave).or.(ifreq.gt.1)).and.
     &    (ifreq.gt.isave_freq))then
c        data for noncentral antenna rays      
c         write(*,*)'before goto 60' 
         goto 60
      endif

c-------------------------------------------------------------------
c     data for the central antenna ray
c-------------------------------------------------------------------
 1    format(15(1pe10.3))
 7    format(7(1pe10.3))
 2    format(/)
      write(*,*)'data for central ray nrayelt',nrayelt

      DO 10 I=1,nrayelt
cSm070720
        if (nrayelt.eq.1) goto 11

c        write(*,*)'2before 82 I',I
c        write(*,*)'REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
c     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
c     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I)),wye(I),wyi(I)'
c        write(*,*)REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
c     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
c     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I)),real(wye(I)),
c     3  REAL(wyi(i))

c        write(*,*)'I,wal_emis(I)',I,wal_emis(I)
        WRITE(82) REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I)),REAL(wye(I)),
     3  REAL(wyi(I)),real(wal_emis(I)),real(wj_emis(I))

        WRITE(81,1)wr(I),wz(I),xarr(I),yarr(I),
     1  ws(I),delpwr(I),rez(I),spsi(I),
     2  wnpar(I),wnper(I),salphal(I),wye(I),wyi(I),
     3  wal_emis(I),wj_emis(I)
         
c       for ion ion resonance for nbulk >2
c        write(*,*)'nbulk,I,wyi(I),wxi2(I),wyi2(I),wxi(I)',
c     &             nbulk,I,wyi(I),wxi2(I),wyi2(I),wxi(I)

c        if (nbulk.gt.2) then
        if (nbulk.gt.3) then
           p12=wyi(I)*wxi2(I)/(wyi2(I)*wxi(I))

c           write(*,*)'p12',p12

           del=1.d0/(1.d0+p12)
           yii=wyi(I)*wyi2(I)*((1.d0-del)*wyi(I)+del*(wyi2(I)))/
     .          ((1.d0-del)*wyi2(I)+del*(wyi(I)))   
        else
           p12=0.d0
           del=1.d0
           yii=0.d0
        endif

c        write(*,*)'mk_graph bef 85'
                       
        WRITE(85)REAL(ws(I)),REAL(delpwr(I)),REAL(spsi(I)),REAL(wye(I)),
     3  REAL(wyi(I)),REAL(wyi2(I)),REAL(dsqrt(yii))

c        write(*,*)'mk_graph bef 86'

        WRITE(86,7)ws(I),delpwr(I),spsi(I),wye(I),
     3  wyi(I),wyi2(I),dsqrt(dabs(yii))

c        write(*,*)'mk_graph bef 75'

        write(75)real(ws(I)),
     +   real(dreal(cwexde(I))),real(dimag(cwexde(I))),
     +   real(dreal(cweyde(I))),real(dimag(cweyde(I))),
     +   real(dreal(cwezde(I))),real(dimag(cwezde(I))),
     +   real(fluxn(I)),REAL(spsi(I))

c        write(*,*)'mk_graph bef 76'

        write(76)real(ws(I)),
     &  real(dreal(w_ceps(1,1,I))),real(dreal(w_ceps(1,2,I))),
     &  real(dreal(w_ceps(1,3,I))),
     &  real(dreal(w_ceps(2,1,I))),real(dreal(w_ceps(2,2,I))),
     &  real(dreal(w_ceps(2,3,I))),
     &  real(dreal(w_ceps(3,1,I))),real(dreal(w_ceps(3,2,I))),
     &  real(dreal(w_ceps(3,3,I)))


c        write(*,*)'mk_graph bef 77'

        write(77)real(ws(I)),
     &  real(dimag(w_ceps(1,1,I))),real(dimag(w_ceps(1,2,I))),
     &  real(dimag(w_ceps(1,3,I))),
     &  real(dimag(w_ceps(2,1,I))),real(dimag(w_ceps(2,2,I))),
     &  real(dimag(w_ceps(2,3,I))),
     &  real(dimag(w_ceps(3,1,I))),real(dimag(w_ceps(3,2,I))),
     &  real(dimag(w_ceps(3,3,I)))

c        write(*,*)'mk_graph bef 78'

        WRITE(78) REAL(wr(I)),REAL(wye(I)),
     &  real(wp_perpmax_dmvt(I,2)),
     &  real(wp_parmin_dmvt(I,2)),real(wp_parmax_dmvt(I,2)),
     &  real(wp_0_dmvt(I,2)),real(wdel_p_dmvt(I,2)),
     &  REAL(wnpar(I))

c      write(*,*)'mk_graph i,wp_0_dmvt(I,2),wdel_p_dmvt(I,2)',
c     &i,wp_0_dmvt(I,2),wdel_p_dmvt(I,2)

c       write(*,*)'mk_graph bef 79'

       WRITE(79) REAL(ws(I)),REAL(wr(I)),REAL(spsi(I)),REAL(wye(I)),
     &  real(eff(I)),REAL(wnpar(I)),real(delpow_e_ar(I)),
     &  real(delcur_par_ar(I)),real(wtheta_pol(I))
     

10    CONTINUE
c      write(*,*)'after 10'
      WRITE(82)
      write(85)
      write(75)
      write(78)
      write(79)
      WRITE(81,2)

cSm070720
 11   continue
      DO 30 J=1,3 
       DO 20 I=1,NP+1
         WRITE(82) REAL(AR(J,I)),REAL(AZ(J,I)),REAL(XT(J,I)),
     1   REAL(YT(J,I)),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt)),real(wye(nrayelt)),real(wyi(nrayelt)),
     5   real(wal_emis(nrayelt)),real(wj_emis(nrayelt))

        WRITE(81,1)AR(J,I),AZ(J,I),XT(J,I),
     1   YT(J,I),
     2   ws(nrayelt),delpwr(nrayelt),rez(nrayelt),
     3   spsi(nrayelt),wnpar(nrayelt),wnper(nrayelt),
     4   salphal(nrayelt),wye(nrayelt),wyi(nrayelt),
     5   wal_emis(nrayelt),wj_emis(nrayelt)
20     CONTINUE
       WRITE(82)
       WRITE(81,2)
30    CONTINUE

      DO 50 J=4,NL
       DO 40 I=1,NP+1
         WRITE(82) REAL(AR(J,I)),REAL(AZ(J,I)),
     1   REAL(XT(3,NP+1)),REAL(YT(3,NP+1)),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt)),real(wye(nrayelt)),real(wyi(nrayelt)),
     5   real(wal_emis(nrayelt)),real(wj_emis(nrayelt))
    
         WRITE(81,1)AR(J,I),AZ(J,I),
     1   XT(3,NP+1),YT(3,NP+1),
     2   ws(nrayelt),delpwr(nrayelt),rez(nrayelt),
     3   spsi(nrayelt),wnpar(nrayelt),wnper(nrayelt),
     4   salphal(nrayelt),wye(nrayelt),wyi(nrayelt),
     5   wal_emis(nrayelt),wj_emis(nrayelt)
40     CONTINUE
       WRITE(82)
       WRITE(81,2)
50    CONTINUE
     
cSAP090313 write wall coordinates

      if (n_wall.gt.0) then
         do i=1,n_wall  
         WRITE(82) 
     1    REAL(r_wall(i)*100.0),REAL(z_wall(i)*100.0),
     &    REAL(xarr(nrayelt)),REAL(yarr(nrayelt)),
     2    REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3    REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4    REAL(salphal(nrayelt)),real(wye(nrayelt)),real(wyi(nrayelt)),
     5    real(wal_emis(nrayelt)),real(wj_emis(nrayelt))
         enddo

         WRITE(82)

      endif

      if (n_wall.gt.1) then        
         do m=0,max_limiters          
            do i=1,n_wall_add(m)           
               WRITE(82) 
     &         REAL(r_wall_add(i,m)*100.0),
     &         REAL(z_wall_add(i,m)*100.0),
     &         REAL(xarr(nrayelt)),REAL(yarr(nrayelt)),
     2    REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4    REAL(salphal(nrayelt)),real(wye(nrayelt)),real(wyi(nrayelt)),
     5         real(wal_emis(nrayelt)),real(wj_emis(nrayelt))  
           enddo ! i=1,n_wall_add(m)       
 
           WRITE(82)
         enddo !m
      endif
         
cSAP090413     
      p=pi/180.d0
      do m=1,max_limiters
         r_min=1.d0
         do i=1,n_limiter(m)
            if(r_min.gt.r_limiter(i,m)) r_min=r_limiter(i,m)
         enddo
         r_max=req(nreqd)*100.d0
         r_min=r_min*100.d0
         write(*,*)'degree phi_limiter(1,m),phi_limiter(2,m)',
     &              phi_limiter(1,m),phi_limiter(2,m)
         write(*,*)'radian phi_limiter(1,m)*p,phi_limiter(2,m)*p',
     &              phi_limiter(1,m)*p,phi_limiter(2,m)*p

         x_lim_min=r_min*dcos(phi_limiter(1,m)*p)
         x_lim_max=r_max*dcos(phi_limiter(1,m)*p)
         y_lim_min=r_min*dsin(phi_limiter(1,m)*p)
         y_lim_max=r_max*dsin(phi_limiter(1,m)*p)
 
         write(*,*)'mk_graph r_min,r_max',r_min,r_max
         write(*,*)'x_lim_min,y_lim_min',x_lim_min,y_lim_min
         write(*,*)'x_lim_max,y_lim_max',x_lim_max,y_lim_max

         WRITE(82) REAL(AR(1,1)),REAL(AZ(1,1)),
     1   REAL(x_lim_min),REAL(y_lim_min),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt)),real(wye(nrayelt)),real(wyi(nrayelt)),
     5   real(wal_emis(nrayelt)),real(wj_emis(nrayelt))

         WRITE(82) REAL(AR(1,1)),REAL(AZ(1,1)),
     1   REAL(x_lim_max),REAL(y_lim_max),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt)),real(wye(nrayelt)),real(wyi(nrayelt)),
     5   real(wal_emis(nrayelt)),real(wj_emis(nrayelt))
         WRITE(82)

         x_lim_min=r_min*dcos(phi_limiter(2,m)*p)
         x_lim_max=r_max*dcos(phi_limiter(2,m)*p)
         y_lim_min=r_min*dsin(phi_limiter(2,m)*p)
         y_lim_max=r_max*dsin(phi_limiter(2,m)*p)

         write(*,*)'x_lim_min,y_lim_min',x_lim_min,y_lim_min
         write(*,*)'x_lim_max,y_lim_max',x_lim_max,y_lim_max

         WRITE(82) REAL(AR(1,1)),REAL(AZ(1,1)),
     1   REAL(x_lim_min),REAL(y_lim_min),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt)),real(wye(nrayelt)),real(wyi(nrayelt)),
     5   real(wal_emis(nrayelt)),real(wj_emis(nrayelt))

         WRITE(82) REAL(AR(1,1)),REAL(AZ(1,1)),
     1   REAL(x_lim_max),REAL(y_lim_max),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt)),real(wye(nrayelt)),real(wyi(nrayelt)),
     5   real(wal_emis(nrayelt)),real(wj_emis(nrayelt))
         WRITE(82)
      enddo !m
c-------------------------------------------------------------------          
      goto 70
c  end data for the central ray
c---------------------------------------------------------------------

60    continue
cyup      write(*,*)'!!!!!!in mk_graph 60 iray,nrayelt',iray,nrayelt
c--------------------------------------------------------------
c        data for noncentral antenna rays
c--------------------------------------------------------------
cSm070729
      if (nrayelt.eq.1) goto 62

      DO 61 I=1,nrayelt

c       write(*,*)'befor 82 i',i
c       write(*,*)'REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
c     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
c     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I)),real(wye(I))'
c       write(*,*)REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
c     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
c     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I)),real(wye(I))


       WRITE(82) REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I)),real(wye(I)),
     3  REAL(wyi(I)),real(wal_emis(I)),real(wj_emis(I))

       WRITE(81,1)wr(I),wz(I),xarr(I),yarr(I),
     1  ws(I),delpwr(I),rez(I),spsi(I),
     2  wnpar(I),wnper(I),salphal(I),wye(I),wyi(I),
     3  wal_emis(I),wj_emis(I)

        write(75)real(ws(I)),
     +   real(dreal(cwexde(I))),real(dimag(cwexde(I))),
     +   real(dreal(cweyde(I))),real(dimag(cweyde(I))),
     +   real(dreal(cwezde(I))),real(dimag(cwezde(I))),
     +   real(fluxn(I)),REAL(spsi(I))

        WRITE(78) REAL(wr(I)),REAL(wye(I)),
     &  real(wp_perpmax_dmvt(I,2)),
     &  real(wp_parmin_dmvt(I,2)),real(wp_parmax_dmvt(I,2)),
     &  real(wp_0_dmvt(I,2)),real(wdel_p_dmvt(I,2)),
     &  REAL(wnpar(I))

c      write(*,*)'mk_graph i,wp_0_dmvt(I,2),wdel_p_dmvt(I,2)',
c     &i,wp_0_dmvt(I,2),wdel_p_dmvt(I,2)

       WRITE(79) REAL(ws(I)),REAL(wr(I)),REAL(spsi(I)),REAL(wye(I)),
     &  real(eff(I)),REAL(wnpar(I)),real(delpow_e_ar(I)),
     &  real(delcur_par_ar(I)),real(wtheta_pol(I))
61    CONTINUE
c      write(*,*)'after 61'
      WRITE(82)
      WRITE(81,2)
      write(75)
      WRITE(78)
      write(79) 

cSm070720 
 62   CONTINUE

c  end data for the noncentral rays
c---------------------------------------------------------------------

70    continue
c      if (iray.eq.nray) then
       
c      write(*,*)'mk_graph iray,nraym,ifreq,nfreq',
c     & iray,nray,ifreq,nfreq
      if ((iray.eq.nray).and.(ifreq.eq.nfreq)) then
c       write(*,*)'mk_graph before close 82'
       CLOSE(82)
       CLOSE(81)
       close(85)
       close(86)
       close(75)
       close(76)
       close(77)
       close(78)
       close(79)
c       write(*,*)'mk_graph after close 82 and 81'
      endif
      RETURN
      END

				    

c----------------subroutine MK_GRAPc----------------
c Prepares for file drawgenc.in for xdraw              
c  +  plots of contours 1/Yc_i=n
c INPUT: from common blocks 'gr' and 'write'        
c         iray -number of the ray at antenna         
c    	  iwopen and iwj are in common one.i         
c iwopen =1 mk_grapc will calculate open contours wb_c=n (using contrb1)
c         2 mk_grapc will calculate close contours wb_c=n (using contrb2)
c iwj    =mk_grapc will calculate contours wb_cj=n.
c         Here j gives the plasma  component, 
c        must be.le.nbulk, j=1 for the electron gyrofrequency
c*****************************************************
      subroutine MK_GRAPc(iray,iwopen,iwj)
      implicit integer (i-n), real*8 (a-h,o-z)
      INTEGER I,J
      include 'param.i'
      include 'gr.i'
      include 'write.i'
c      write(*,*)'!!!!!!in mk_graphc iray,nrayelt',iray,nrayelt
      if(iray.eq.1) then
c         CALL ASSIGN("assign -F f77 -N ieee u:84",ier)
         OPEN(84,file='genrac.bin',form='unformatted')
         open(83,file='genrac.doc')

      end if
      if (nrayelt.eq.0)then
         goto 70
      endif
      if (iray.gt.1) then
c        data for noncentral antenna rays
         goto 60
      endif

c-------------------------------------------------------------------
c     data for the central antenna ray
c-------------------------------------------------------------------
 1    format(11(1pe10.3))
 2    format(/)
      
      do i=1,nrayelt
       WRITE(83,1)wr(I),wz(I),xarr(I),yarr(I),
     1 ws(I),delpwr(I),rez(I),spsi(I),
     2 wnpar(I),wnper(I),salphal(I)
      enddo
      DO 10 I=1,nrayelt
cSm070720
        if (nrayelt.eq.1) goto 11

       WRITE(84) REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I))
10    CONTINUE
      WRITE(84)
      WRITE(83,2)
cc      write(*,*)'in mk_graph before 30'
cc       read(*,*)

cSm070720
 11     continue

      DO 30 J=1,3
       DO 20 I=1,NP+1
        WRITE(84) REAL(AR(J,I)),REAL(AZ(J,I)),REAL(XT(J,I)),
     1   REAL(YT(J,I)),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt))
        WRITE(83,1)AR(J,I),AZ(J,I),XT(J,I),
     1   YT(J,I),
     2   ws(nrayelt),delpwr(nrayelt),rez(nrayelt),
     3   spsi(nrayelt),wnpar(nrayelt),wnper(nrayelt),
     4   salphal(nrayelt)
20     CONTINUE
       WRITE(84)
       WRITE(83,2)
30    CONTINUE
cc      write(*,*)'in mk_grapc after 30'
cc      read(*,*)
cc      write(*,*)'in mk_grapc before 50'
      DO 50 J=4,NL
       DO 40 I=1,NP+1
        WRITE(84) REAL(AR(J,I)),REAL(AZ(J,I)),
     1   REAL(XT(3,NP+1)),REAL(YT(3,NP+1)),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt))
        WRITE(83,1)AR(J,I),AZ(J,I),
     1   XT(3,NP+1),YT(3,NP+1),
     2   ws(nrayelt),delpwr(nrayelt),rez(nrayelt),
     3   spsi(nrayelt),wnpar(nrayelt),wnper(nrayelt),
     4   salphal(nrayelt)
40     CONTINUE
       WRITE(84)
       WRITE(83,2)
50    CONTINUE
c      write(*,*)'in mk_graphc after 50'
c	write(*,*)'in mk_graphc before call contrb1 iwopen,iwj',
c     *  iwopen,iwj 
        if (iwopen.eq.1) then
c	  open contours R=R(z), iwj are in common one.i
          call contrb1
	else 
          if(iwopen.eq.2) then
c           close contours rho=rho(theta), iwj are in common one
cSm030514
            write(*,*)'in mk_graph iwopen=',iwopen
            call contrb2
            write(*,*)'after contrb2'
          endif
 	endif
c	write(*,*)'in mk_graphc after call contrb1'
      goto 70
c  end data for the central ray
c---------------------------------------------------------------------

60    continue
cc      write(*,*)'!!!!!!in mk_graph 60 iray,nrayelt',iray,nrayelt
c--------------------------------------------------------------
c        data for noncentral antenna rays
c--------------------------------------------------------------
c      p=0.0
cc      write(*,*)'in mk_graph before 61'
cSm070720
        if (nrayelt.eq.1) goto 62

      DO 61 I=1,nrayelt
       WRITE(84) REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
c     1  REAL(p),REAL(p),REAL(p),REAL(p),
c     2  REAL(p),REAL(p),REAL(p)
     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I))
cc       write(*,*)'in mk_graph 61 i',i
cc       WRITE(*,1)wr(I),wz(I),xarr(I),yarr(I),
cc     1  ws(I),delpwr(I),rez(I),spsi(I),
cc     2  wnpar(I),wnper(I),salphal(I)
       WRITE(83,1)wr(I),wz(I),xarr(I),yarr(I),
     1  ws(I),delpwr(I),rez(I),spsi(I),
     2  wnpar(I),wnper(I),salphal(I)
61    CONTINUE
cc      write(*,*)'in mk_graph after 61 30'
      WRITE(84)
      WRITE(83,2)

cSm070720 
 62   CONTINUE


c  end data for the noncentral rays
c---------------------------------------------------------------------

70    continue
c      CLOSE(84)
c      CLOSE(83)
      RETURN
      END


c----------------subroutine mk_gronetwo-------------*
c Prepares for file drawonet.in for xdraw            *
c INPUT: from common block  'onetwo.i'              *
c****************************************************
      subroutine mk_gronetwo
      implicit integer (i-n), real*8 (a-h,o-z)
      INTEGER I
      include 'param.i'
      include 'onetwo.i'
cSAP090404
      include 'rho.i'
      include 'one.i'
  
c         CALL ASSIGN("assign -F f77 -N ieee u:11",ier)
      OPEN(11,file='onetwo.bin',form='unformatted')
      OPEN(21,file='onetwo.dat')
20    format(8e12.4)
      WRITE(21,*)'rho spower powden powden_e
     1powden_i powden_cl'
     
      hrho= rho_bin(2)-rho_bin(1) ! YuP: was: hrho=1.d0/(NR-1)
      p=1.d0/(hrho*dsqrt(1.d4*areatot/pi)) !1/[cm]

      DO 10 I=1,NR-1
         rho=rho_bin_center(i)  ! YuP: was: hrho*(I-1)
         WRITE(11) REAL(rho),REAL(spower(I)),REAL(powden(I)),
     1   REAL(powden_e(I)),REAL(powden_i(I)),REAL(powden_cl(I)),
     &   real(spower(I)*1.d-7*p),real(rho*dsqrt(1.d4*areatot/pi))   
         WRITE(21,20)rho,spower(I),powden(I),
     1   powden_e(I),powden_i(I),powden_cl(I)
10    CONTINUE
      WRITE(11)
     
      CLOSE(11)
      CLOSE(21)

      write(*,*)'mk_graph.f: mk_gronetwo created onewto.bin onetwo.dat'
    
      RETURN
      END


c----------------subroutine MK_GRAPT----------------*
c Prepares for file drawgenr.in for xdraw toray_data  *
c INPUT: from common blocks 'gr' and 'write'        *
c        iray -number of the ray at antenna         *
c****************************************************
      subroutine MK_GRAPT (iray,nray)
      implicit integer (i-n), real*8 (a-h,o-z)
      INTEGER I,J

      integer isave
      save isave

      include 'param.i'
      include 'gr.i'
      include 'write.i'

      if (isave.lt.1) then
        isave=iray
      endif

c      write(*,*)'!!!!!!in mk_graph iray,nrayelt,isave',
c     + iray,nrayelt,isave

c      if(iray.eq.1) then
       if((iray.eq.1).or.(isave.eq.iray)) then
c         CALL ASSIGN("assign -F f77 -N ieee u:52",ier)
c         write(*,*)'mk_grph before open 52 and 51'
         OPEN(52,file='genrayt.bin',form='unformatted')
         open(51,file='genrayt.doc')
         open(55,file='absorpt.bin',form='unformatted')
         open(56,file='absorpt.doc')
      end if
      if (nrayelt.eq.0)then
         goto 70
      endif
c      if (iray.gt.1) then
      if (iray.gt.isave) then
c        data for noncentral antenna rays
         goto 60
      endif

c-------------------------------------------------------------------
c     data for the central antenna ray
c-------------------------------------------------------------------
 1    format(11(1pe10.3))
 7    format(7(1pe10.3))
 2    format(/)
c      write(*,*)'in mk_graph nrayelt=',nrayelt,'iray',iray

      DO 10 I=1,nrayelt

c        write(*,*)'before 82 I',I
c        write(*,*)'REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
c     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
c     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I)),wye(I),wyi(I)'
c        write(*,*)REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
c     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
c     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I)),real(wye(I)),
c     3  REAL(wyi(i))

cSm070720
        if (nrayelt.eq.1) goto 11

        WRITE(52) REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I)),REAL(wye(I)),
     3  REAL(wyi(I))
        WRITE(51,1)wr(I),wz(I),xarr(I),yarr(I),
     1  ws(I),delpwr(I),rez(I),spsi(I),
     2  wnpar(I),wnper(I),salphal(I),wye(I),wyi(I)
         
c       for ion resonance for nbulk >2
        p12=wyi(I)*wxi2(I)/(wyi2(I)*wxi(I))
        del=1.d0/(1.d0+p12)
        yii=wyi(I)*wyi2(I)*((1.d0-del)*wyi(I)+del*(wyi2(I)))/
     .   ((1.d0-del)*wyi2(I)+del*(wyi(I)))   
                   
        WRITE(55)REAL(ws(I)),REAL(delpwr(I)),REAL(spsi(I)),REAL(wye(I)),
     3  REAL(wyi(I)),REAL(wyi2(I)),REAL(dsqrt(yii))
        WRITE(56,7)ws(I),delpwr(I),spsi(I),wye(I),
     3  wyi(I),wyi2(I),dsqrt(dabs(yii))
10    CONTINUE
      WRITE(52)
      write(55)
      WRITE(51,2)
c      write(*,*)'in mk_graph before 30'
cc      read(*,*)

cSm070720
 11   continue

      DO 30 J=1,3
       DO 20 I=1,NP+1
        WRITE(52) REAL(AR(J,I)),REAL(AZ(J,I)),REAL(XT(J,I)),
     1   REAL(YT(J,I)),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt)),real(wye(nrayelt)),real(wyi(nrayelt))
        WRITE(51,1)AR(J,I),AZ(J,I),XT(J,I),
     1   YT(J,I),
     2   ws(nrayelt),delpwr(nrayelt),rez(nrayelt),
     3   spsi(nrayelt),wnpar(nrayelt),wnper(nrayelt),
     4   salphal(nrayelt),wye(nrayelt),wyi(nrayelt)
20     CONTINUE
       WRITE(52)
       WRITE(51,2)
30    CONTINUE
c      write(*,*)'in mk_graph before 50'
      DO 50 J=4,NL
       DO 40 I=1,NP+1
        WRITE(52) REAL(AR(J,I)),REAL(AZ(J,I)),
     1   REAL(XT(3,NP+1)),REAL(YT(3,NP+1)),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt)),real(wye(nrayelt)),real(wyi(nrayelt))
        WRITE(51,1)AR(J,I),AZ(J,I),
     1   XT(3,NP+1),YT(3,NP+1),
     2   ws(nrayelt),delpwr(nrayelt),rez(nrayelt),
     3   spsi(nrayelt),wnpar(nrayelt),wnper(nrayelt),
     4   salphal(nrayelt),wye(nrayelt),wyi(nrayelt)
40     CONTINUE
       WRITE(52)
       WRITE(51,2)
50    CONTINUE
c      write(*,*)'in mk_graph after 50'
      goto 70
c  end data for the central ray
c---------------------------------------------------------------------

60    continue
cc      write(*,*)'!!!!!!in mk_graph 60 iray,nrayelt',iray,nrayelt
c--------------------------------------------------------------
c        data for noncentral antenna rays
c--------------------------------------------------------------
cc      write(*,*)'in mk_graph before 61'

cSm070720
      if (nrayelt.eq.1) goto 62

      DO 61 I=1,nrayelt

       WRITE(52) REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I)),real(wye(I)),
     3  REAL(wyi(I))
       WRITE(51,1)wr(I),wz(I),xarr(I),yarr(I),
     1  ws(I),delpwr(I),rez(I),spsi(I),
     2  wnpar(I),wnper(I),salphal(I),wye(I),wyi(I)
61    CONTINUE
      WRITE(52)
      WRITE(51,2)
cSm070720 
 62   continue
c  end data for the noncentral rays
c---------------------------------------------------------------------

70    continue
      if (iray.eq.nray) then
       CLOSE(52)
       CLOSE(51)
       close(55)
       close(56)


c       write(*,*)'mk_grapht  after close 52 and 51'
      endif
      RETURN
      END



c----------------subroutine mkgrtool-------------*
c Prepares for file tool.bin for xdraw              *
c INPUT: from common blocks: onetwo, one, 'five' *
c************************************************
      subroutine mkgrtool
      implicit integer (i-n), real*8 (a-h,o-z)
      INTEGER I
      include 'param.i'
      include 'onetwo.i'
      include 'one.i'
      include 'five.i'

c     CALL ASSIGN("assign -F f77 -N ieee u:11",ier)
      write(*,*)'in mkgrtool'
      OPEN(11,file='tool.bin',form='unformatted')
      OPEN(21,file='tool.dat')
20    format(8e12.4)
      
      hrho= rho_bin(2)-rho_bin(1) ! YuP: was: hrho=1.d0/(NR-1)
      WRITE(21,*) 'rho,dens_e,temp_e,bmod,ye,
     1 xe,ur,cutoffp,cutoffm,rx'
      thetax=0.d0

      write(*,*)' in mkgrtool before do 30'

c-----write data to tool.bin and tool.dat files

      DO 30 I=1,NR
         rho=rho_bin(i) ! YuP: was: hrho*(I-1)

         psix=psi_rho(rho)
         call zr_psith(psix,thetax,zx,rx)
         bmod=b(zx,rx,0.d0)
         ye=y(zx,rx,0.d0,1)
         xe=x(zx,rx,0.d0,1)
         ur=xe+ye**2
         cnpar=0.d0
         cutoffp=xe/(1.d0+ye)+cnpar**2
         if (ye.ne.0.d0) then
          cutoffm=xe/(1.d0-ye)+cnpar**2
         else
          cutoffm=0.d0
         endif
         
         WRITE(11) REAL(rho),REAL(densrho(rho,1)),
     1   REAL(temperho(rho,1)),REAL(bmod),REAL(ye),
     1   REAL(xe),REAL(ur),REAL(cutoffp),REAL(cutoffm),
     1   REAL(rx)
    
         WRITE(21,20)rho,densrho(rho,1),
     1   temperho(rho,1),bmod,ye,
     1   xe,ur,cutoffp,cutoffm,
     1   rx
 30   CONTINUE  

      write(*,*)' in mkgrtool after do 30'

      WRITE(11)
     
      CLOSE(11)
      CLOSE(21)
     
      OPEN(11,file='tool1.bin',form='unformatted')
      OPEN(21,file='tool1.dat')

c-----data for x-mode optimal parameters
      OPEN(31,file='tool2.bin',form='unformatted')
      OPEN(41,file='tool2.dat')
c---------------------------------------


      N_r=50
c      write(*,*)'rmax,rmin',rmax,rmin
      hr=(rmax-rmin)/(N_r-1)
      WRITE(21,*) 'r,dens_e,temp_e bmod,y_e,x_e,ur,cutoofp,
     1cutoffm,omega_p/omega_b'
 
      nmax=0
      mmax=0

      write(*,*)' in mkgrtool before do 10'

      DO 10 I=1,N_r
         r=rmin+hr*(I-1)
c         write(*,*)'i,r',i,r
         bmod=b(0.d0,r,0.d0)
c         write(*,*)'bmod',bmod
         ye=y(0.d0,r,0.d0,1)
         xe=x(0.d0,r,0.d0,1)
         ur=xe+ye**2
c         write(*,*)'ye,xe,ur',ye,xe,ur

         cnpar=0.d0
         cutoffp=xe/(1.d0+ye)+cnpar**2
         if (ye.ne.0.d0) then
          cutoffm=xe/(1.d0-ye)+cnpar**2
         else
          cutoffm=0.d0
         endif

c--------sqrt(x(r))=omega_p(r)/omega=eta o-mode cutoff for eta=1. 
c        y(r)=omega_b(r)/omega=1/n  'n' harmonics EC resonance condition
         eta=0.95d0  


         WRITE(11) REAL(r),REAL(densrho(rho,1)),
     1   REAL(temperho(rho,1)),REAL(bmod),REAL(ye),
     1   REAL(xe),REAL(ur),REAL(cutoffp),REAL(cutoffm),
     1   REAL(dsqrt(xe)/(eta*ye))
        
         WRITE(21,20)r,densrho(rho,1),
     1   temperho(rho,1),bmod,ye,xe,ur,cutoffp,cutoffm,
     1   dsqrt(xe)/(eta*ye)
         
         n=dsqrt(xe)/(eta*ye)
         if(n.gt.nmax) then
           nmax=n
           r_l=r
         endif        
         
         WRITE(31) REAL(r),REAL(densrho(rho,1)),
     1   REAL(temperho(rho,1)),REAL(bmod),REAL(ye),
     1   REAL(xe),REAL(ur),REAL(cutoffp-1.d0),REAL(cutoffm-1.d0),
     1   REAL(1.d0/(xe/(1.d0-cnpar**2)-1.d0))
        
         WRITE(41,20)r,densrho(rho,1),
     1   temperho(rho,1),bmod,ye,xe,ur,cutoffp-1.,cutoffm-1.,
     1   (1.d0/(xe/(1.d0-cnpar**2)-1.d0))

         m=1.d0/(xe/(1.d0-cnpar**2)-1.d0)
         
         if(iabs(m).gt.mmax) then
           mmax=iabs(m)

           if(m.lt.0)then
             i_mmax=-1
           else
             i_mmax=1
           endif

           r_lm=r
c           write(*,*)'mmax,r_lm',mmax,r_lm
         endif
c         write(*,*)'0 r_lm',r_lm
     
10    CONTINUE

      write(*,*)' in mkgrtool after do 10'

c      write(*,*)'1 r_lm',r_lm     
      WRITE(11)
      WRITE(31)

      CLOSE(11)
      CLOSE(21)

      CLOSE(31)
      CLOSE(41)
c      write(*,*) 'mk_graph.f mkgrtool nmax',nmax
c      write(*,*) 'mk_graph.f mkgrtool mmax,i_mmax',mmax,i_mmax
      z=0.d0
c      write(*,*)'2 r_lm',r_lm    

c      write(*,*)' in mkgrtool before do i=1,nmax nmax=',nmax

c      write(*,*)'eta',eta
 
      goto 100

      do i=1,nmax
       write(*,*)'i',i
       call solvropt(i,z,eta,r_l,rmax,r_opt_o)
       bmod=b(z,r_opt_o,0.d0)
c       write(*,*)'mk_graph.f mkgrtool i,r_opt_o,rho',i,r_opt_o,rho
       
       den_opt=dense(z,r_otp_o,0.d0,1) ! the optimal density
c------omega_optimal=omega_pe/eta
       friq_o=dsqrt(806.2*den_opt/eta**2)
c       write(*,*)'den_opt,friq_o',den_opt,friq_o

      enddo

      write(*,*)' in mkgrtool after do i=1,nmax'


c      write(*,*)'3 r_lm',r_lm     
      mmax=mmax*i_mmax
      k1=i_mmax
      k2=i_mmax
      cnpar=0.d0

      write(*,*)' in mkgrtool before do i=k1,nmax,k2'

      do i=k1,mmax,k2
        write(*,*)'i',i
c       write(*,*)'mk_graph before solvropx i=',i
c       write(*,*)'r_lm,rmax',r_lm,rmax
       call solvropx(i,z,r_lm,rmax,cnpar,r_opt_x)
       bmod=b(z,r_opt_x,0.d0)
c       write(*,*)'mk_graph.f mkgrtool i,r_opt_x,rho',i,r_opt_x,rho
       
       den_opt=dense(z,r_otp_x,0.d0,1)
c------omega_optimal=omega_pe/eta
       friq_x=dsqrt(806.2*den_opt/((1-cnpar**2)*(1.d0+1.d0/i)))
c       write(*,*)'den_opt,friq_x',den_opt,friq_x

      enddo

c      write(*,*)' in mkgrtool after do i=k1,nmax,k2'
 100  continue
      write(*,*) 'end of mkgrtool'

      RETURN
      END


      subroutine solvropt(n,z,eta,r_left,r_right,r_opt_o)
c     calculates the root r_opt_o from the equation( for o-mode)
c     sqrt(x_e(r,z)/(eta*y_e(r,z))=n
c     using the binary method
c     on (r_left<r<r_right) interval
c     input:
c     n is a number of harmonics,
c     z is vericle variable
c     eta is a parameter ~<1,
c     r_left,r_right are the boundaries of the interval
c
c     output:
c     r_opt_o is the root of the equation 
      
      implicit integer (i-n), real*8 (a-h,o-z)
      
      include 'param.i'
      include 'one.i'
cSm030514      
c      common /cf_r_o/z_loc,n_loc,eta_loc
      common /cf_r_o/z_loc,eta_loc,n_loc

      external f_r_o 

      z_loc=z
      n_loc=n
      eta_loc=eta

      racc=1.d-6
      
      r_opt_o=rtbis(f_r_o,r_left,r_right,racc)
           
      return
      end 

 
      double precision function f_r_o(r)
      implicit integer (i-n), real*8 (a-h,o-z)      
      include 'param.i'
      include 'one.i'
cSm030514
c      common /cf_r_o/z_loc,n_loc,eta_loc
      common /cf_r_o/z_loc,eta_loc,n_loc

c-----input
c     r is a major radius
c     z_loc is a vertical distance: from common /cf_r_o/
c     n_loc is the EC hatmonic number: from common /cf_r_o/
c     eta_loc ia a parameter xe=eta**2: from common /cf_r_o/

c      write(*,*)'f_r_o ,z_loc,n_loc,eta_loc,r',z_loc,n_loc,eta_loc,r
      z_loc1=z_loc
      bmod=b(z_loc1,r,0.d0)      
      ye=y(z_loc1,r,0.d0,1)
      xe=x(z_loc1,r,0.d0,1)
      
      n_loc1=n_loc
      eta_loc1=eta_loc
      if (ye.ne.0.d0) then
        f_r_o=dsqrt(xe)/(eta_loc1*ye)-n_loc1
      else
        f_r_o=1000.d0
      endif

      return
      end

      subroutine solvropx(n,z,r_left,r_right,cnpar,r_opt_x)
c     calculates the root r_opt_o from the equation( for x-mode)
c     1/(x/(1-N_par**2)-1)=n
c     using the binary method
c     on (r_left<r<r_right) interval
c     input:
c     n is a number og harmonics,
c     z is vericle variable
c     r_left,r_right are the boundaries of the interval
c     cnpar is N_par
c     output:
c     r_opt_o is the root of the equation 
      
      implicit integer (i-n), real*8 (a-h,o-z)
      
      include 'param.i'
      include 'one.i'
      
c      common /cf_r_x/z_loc,n_loc,cnpar_l
      common /cf_r_x/z_loc,cnpar_l,n_loc


      external f_r_x 

      z_loc=z
      n_loc=n
      cnpar_l=cnpar 

      racc=1.d-6

      r_opt_x=rtbis(f_r_x,r_left,r_right,racc)

           
      return
      end

      double precision function f_r_x(r)
      implicit integer (i-n), real*8 (a-h,o-z)      
      include 'param.i'
      include 'one.i'
c      common /cf_r_x/z_loc,n_loc,cnpar_l
      common /cf_r_x/z_loc,cnpar_l,n_loc


c-----input
c     r is a major radius
c     z_loc is avertical distance: from common /cf_r_x/
c     n_loc is the EC hatmonic number: from common /cf_r_x/
c     cnpar_l is N_parrael from common/cf_r_x/

c      write(*,*)'f_r_x ,z_loc,n_loc,cnpar_l,r',
c     *z_loc,n_loc,cnpar_l,r
      z_loc1=z_loc
      bmod=b(z_loc1,r,0.d0)      
      ye=y(z_loc1,r,0.d0,1)
      xe=x(z_loc1,r,0.d0,1)
      
      n_loc1=n_loc
      cnpar_l1=cnpar_l

      if((1.d0-cnpar_l1**2.ne.0.d0).and.
     *((xe/(1.d0-cnpar_l1**2)-1).ne.0.d0))then
        f_r_x=1/(xe/(1.d0-cnpar_l1**2)-1)-n_loc1
      else
        f_r_x=1000.d0
      endif

      return
      end




c----------------subroutine MK_GR3----------------*
c Prepares for file drawgr3d.in for xdraw              *
c INPUT: from common blocks 'gr' and 'write'        *
c         iray -number of the ray at antenna         *
c*****************************************************
      subroutine MK_GR3d(iray,nray)
      implicit integer (i-n), real*8 (a-h,o-z)
      INTEGER I,J

      integer isave
      save isave

      include 'param.i'
      include 'gr.i'
      include 'write.i'

      if (isave.lt.1) then
        isave=iray
      endif

c      write(*,*)'!!!!!!in mk_gr3d iray,nrayelt,isave',
c     + iray,nrayelt,isave

       if((iray.eq.1).or.(isave.eq.iray)) then
c         CALL ASSIGN("assign -F f77 -N ieee u:82",ier)
         write(*,*)'mk_gr3d before open 87 '
         OPEN(71,file='gr3d_1.bin',form='unformatted')
         open(73,file='gr3d_2.bin',form='unformatted')         
         open(72,file='gr3d_1,doc')
         open(74,file='gr3d_2,doc')
      end if
 11   format(10(1pe10.3))
 9    format(9(1pe10.3))
 1    format(/)

      if (nrayelt.eq.0)then
         goto 70
      endif

      write(*,*)'mk_gr3d iray,isave,nrayelt',iray,isave,nrayelt
      if (iray.gt.isave) then
c        data for noncentral antenna rays
         goto 60
      endif

c-------------------------------------------------------------------
c     data for the central antenna ray
c-------------------------------------------------------------------
cSm070720
      if (nrayelt.eq.1) goto 12

      DO 10 I=1,nrayelt

        WRITE(71) real(ws(I)),REAL(seikon(I)),real(spsi(I)),
     +  REAL(wr(I)),real(wphi(I)),REAL(wz(I)),
     2  REAL(wnpar(I)),REAL(wnper(I)),real(delpwr(I)),real(sdpwr(I))
  
        write(72,11)real(ws(I)),REAL(seikon(I)),real(spsi(I)),
     +  REAL(wr(I)),real(wphi(I)),REAL(wz(I)),
     2  REAL(wnpar(I)),REAL(wnper(I)),real(delpwr(I)),real(sdpwr(I))

        write(73)real(ws(I)),real(wdnpar(I)),real(dreal(cwexde(I))),
     +  real(dreal(cweyde(I))),real(dreal(cwezde(I))),
     +  real(fluxn(I)),real(sbtot(I)),real(sene(I)),
     +  REAL(salphac(I)),REAL(salphal(I))

        write(74,9)real(ws(I)),real(wdnpar(I)),real(dreal(cwexde(I))),
     +  real(dreal(cweyde(I))),real(dreal(cwezde(I))),
     +  real(fluxn(I)),real(sbtot(I)),real(sene(I)),
     +  REAL(salphac(I)),REAL(salphal(I))

10    CONTINUE
      WRITE(71)
      write(73)
      write(72,1)
      write(74,1)

cSm070720
 12   continue

      goto 70
c  end data for the central ray
c--------------------------------------------------------------
60    continue
c--------------------------------------------------------------
c        data for noncentral antenna rays
c--------------------------------------------------------------
cSm070720
      if (nrayelt.eq.1) goto 62

      DO 61 I=1,nrayelt
        WRITE(71) real(ws(I)),REAL(seikon(I)),real(spsi(I)),
     +  REAL(wr(I)),real(wphi(I)),REAL(wz(I)),
     2  REAL(wnpar(I)),REAL(wnper(I)),real(delpwr(I)),real(sdpwr(I))
  
        write(72,11)real(ws(I)),REAL(seikon(I)),real(spsi(I)),
     +  REAL(wr(I)),real(wphi(I)),REAL(wz(I)),
     2  REAL(wnpar(I)),REAL(wnper(I)),real(delpwr(I)),real(sdpwr(I))

        write(73)real(ws(I)),real(wdnpar(I)),real(dreal(cwexde(I))),
     +  real(dreal(cweyde(I))),real(dreal(cwezde(I))),
     +  real(fluxn(I)),real(sbtot(I)),real(sene(I)),
     +  REAL(salphac(I)),REAL(salphal(I))

        write(74,9)real(ws(I)),real(wdnpar(I)),real(dreal(cwexde(I))),
     +  real(dreal(cweyde(I))),real(dreal(cwezde(I))),
     +  real(fluxn(I)),real(sbtot(I)),real(sene(I)),
     +  REAL(salphac(I)),REAL(salphal(I))




61    CONTINUE
      
      WRITE(71)
      write(73)
      write(72,1)
      write(74,1)

cSm070720
 62   continue

c  end data for the noncentral rays
c---------------------------------------------------------------------
70    continue
      if (iray.eq.nray) then
       CLOSE(71)
       CLOSE(72)
       CLOSE(73)
       CLOSE(74)      
      endif

      write(*,*)'end MK_GR3d' 

      RETURN
      END

	


c----------------subroutine MK_gremis----------------
c Prepares for file drawem.in for xdraw              
c INPUT: from common blocks 'gr' and 'write'        
c         iray -number of the ray at antenna         
c*****************************************************
      subroutine mk_gremis(iray,nray,ifreq,nfreq)
      implicit integer (i-n), real*8 (a-h,o-z)
      INTEGER I,J

      integer isave
cSm070108
      integer isave_freq

      save isave  
      save isave_freq

      include 'param.i'
      include 'gr.i'
      include 'write.i'

      data isave /0/
      data isave_freq/0/

      if (isave.lt.1) then
        isave=iray
      endif

cSm061229
      if(isave_freq.lt.1) then
        isave_freq=ifreq
      endif

       if((iray.eq.1).or.(isave.eq.iray)) then
         OPEN(90,file='emis.bin',form='unformatted')
         open(91,file='emis.doc') 
      end if

      if (nrayelt_emis.lt.3)then
        if((iray.gt.isave).or.(ifreq.gt.1)) then
         goto 70
        endif
      endif
 

       if (((iray.gt.isave).or.(ifreq.gt.1)).and.
     &    (ifreq.gt.isave_freq))then
c        data for noncentral antenna rays
c         write(*,*)'before goto 60' 
         goto 60
      endif

c-------------------------------------------------------------------
c     data for the central antenna ray
c-------------------------------------------------------------------
 1    format(13(1pe10.3))
 2    format(/)
    
 11   continue

      goto 70
c  end data for the central ray
c---------------------------------------------------------------------

60    continue

c--------------------------------------------------------------
c        data for noncentral antenna rays
c--------------------------------------------------------------

      if (nrayelt_emis.lt.3) goto 62
 
 62   continue

c  end data for the noncentral rays
c---------------------------------------------------------------------

70    continue
c      if (iray.eq.nray) then
c      write(*,*)'mk_graph iray,nraym,ifreq,nfreq',
c     & iray,nray,ifreq,nfreq
      if ((iray.eq.nray).and.(ifreq.eq.nfreq)) then
c       write(*,*)'mk_graph before close 90'
       CLOSE(90)
       CLOSE(91)     

c       write(*,*)'mk_gremis after close 90 and 91'
      endif
      RETURN
      END


c----------------subroutine mk_gremfr----------------
c Prepares for file drawemfr.in for xdraw              
c INPUT: from common blocks 'gr' and 'write'        
c         iray -number of the ray at antenna         
c*****************************************************
      subroutine mk_gremfr(iray,nray)
      implicit integer (i-n), real*8 (a-h,o-z)
      INTEGER I,J

      integer isave
      save isave

      include 'param.i'
      include 'one.i'
      include 'write.i'

      if (isave.lt.1) then
        isave=iray
      endif

c      write(*,*)'!!!!!!in mk_gremfr iray,isave',
c     + iray,isave

      if((iray.eq.1).or.(isave.eq.iray)) then
c         CALL ASSIGN("assign -F f77 -N ieee u:94",ier)
c         write(*,*)'mk_gremfr before open 93 and 94'
         OPEN(95,file='emfr.bin',form='unformatted')
         open(93,file='emfr.doc') 
      end if
    
      if (iray.gt.isave) then
c        data for noncentral antenna rays
         goto 60
      endif

c-------------------------------------------------------------------
c     data for the central antenna ray
c-------------------------------------------------------------------
 1    format(10(1pe10.3))
 2    format(/)
      write(*,*)'in mk_gremfr iray',iray,'nfreq',nfreq
      write(*,*)'mk_gremf freqncy0',freqncy0


      DO 10 I=1,nfreq

c        write(*,*)'in mk_gremfr I,wfreq(I),wtemp_rad_fr(I),
c     +  wi_0(iray,I),wtemp_rad_fr_wall(I)',
c     +  I,wfreq(I),wtemp_rad_fr(I),wi_0(iray,I),
c     +  wtemp_rad_fr_wall(I)

c        write(*,*)'in mk_gremfr wtemp_pl_fr(I),wr0_em(I),wz0_em(I),
c     +  wrho0_em(I)',
c     +  wtemp_pl_fr(I),
c     +  wr0_em(I),wz0_em(I),wrho0_em(I)

        WRITE(95) REAL(wfreq(I)/freqncy0),Real(wtemp_rad_fr(I)),
     +  REAL(wi_0(iray,I)),real(wtemp_rad_fr_wall(I)),
     +  real(wtemp_pl_fr(I)),real(wr0_em(I)),real(wz0_em(I)),
     +  real(wrho0_em(I)),real(wr_2nd_harm(I)),
     +  real(wtemp_2nd_harm(I)),real(wtau_em(iray,I))      

        WRITE(93,1) wfreq(I)/freqncy0,wtemp_rad_fr(I),wi_0(iray,I),
     +  wtemp_rad_fr_wall(I),wtemp_pl_fr(I),wr0_em(I),wz0_em(I),
     +  wrho0_em(I),wr_2nd_harm(I),wtemp_2nd_harm(I),
     +  wtau_em(iray,I)    


10    CONTINUE
      WRITE(95)
      WRITE(93,2)
     
      goto 70
c  end data for the central ray
c---------------------------------------------------------------------

60    continue

c--------------------------------------------------------------
c        data for noncentral antenna rays
c--------------------------------------------------------------
      
      DO 61 I=1,nfreq
        
        WRITE(95) REAL(wfreq(I)/freqncy0),REAL(wtemp_rad_fr(I)),
     +  REAL(wi_0(iray,I)),real(wtemp_rad_fr_wall(I)),
     +  real(wtemp_pl_fr(I)),real(wr0_em(I)),real(wz0_em(I)),
     +  real(wrho0_em(I)),real(wr_2nd_harm(I)),
     +  real(wtemp_2nd_harm(I)),real(wtau_em(iray,I))    
        WRITE(93,1) wfreq(I),wtemp_rad_fr(I),wi_0(iray,I),
     +  wtemp_rad_fr_wall(I),wtemp_pl_fr(I),wr0_em(I),wz0_em(I),
     +  wrho0_em(I),wr_2nd_harm(I),wtemp_2nd_harm(I),
     +  wtau_em(iray,I)
    
61    CONTINUE
      WRITE(95)
      WRITE(93,2)
c  end data for the noncentral rays
c---------------------------------------------------------------------

70    continue

      if (iray.eq.nray) then
         CLOSE(95)
         CLOSE(93)
         write(*,*)'mk_gremfr after close 95 and 93'
      endif

      RETURN
      END


      
      subroutine read_emfr_bin(iray)
c-----reads emfr.bin file        
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.i'
      include 'one.i'
      include 'write.i'
      real a1,a2,a3,a4,a5,a6,a7,a8,a9,a10

      OPEN(95,file='emfr.bin',form='unformatted',status='old')
      write(*,*)' in read_emfr_bin nfreq=',nfreq
      DO I=1,nfreq
      write(*,*)'I',I
      READ(95) a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
      enddo  
      close(95)

      RETURN
      END
       

      subroutine map_d_cold(z,r,phi,npar,n_nperp,nperp_min,nperp_max,
     &name_param,param,n_param)
c-----creates----------------------------------------------------------
c     1D array D(nperp) for plotting of the cold plasma dispersion function
c     D(nperp) versus perpendicular refractive index N_perpendicular=nperp
c     at the given point (r,z,phi) and given N_parallel=npar
c
c     Plots D_cold(ReN_perp) to plot.ps file
c----------------------------------------------------------------------------
      implicit none
c-----input:
      integer n_nperp                         !number of points in nperp mesh
      double precision nperp_min,nperp_max,   !minimal and maximal nperp values
     &z,r,phi, ! space point coordinates
     &npar     ! parallel refractive index
      integer n_param                    !number of input parameters
      character*(*) name_param(*)        !names  of the input parameters
      real param(*)                      !values of the input parameters    
c-----locals
      integer n_nperp_a !maximal number of n_nperp
c      parameter(n_nperp_a=10000)
      parameter(n_nperp_a=100000)
      integer i

      double precision d_ar(n_nperp_a), !array of dispersion function values
     &nperp_ar(n_nperp_a)   !mesh of nperp

c      real d_ar(n_nperp_a), !array of dispersion function values
c     &nperp_ar(n_nperp_a)   !mesh of nperp

      double precision step,nperp
      double complex disp_func

      write(*,*)'in subroutine map_d_cold z,r,phi',z,r,phi
      write(*,*)'npar,n_nperp,nperp_min,nperp_max',
     &npar,n_nperp,nperp_min,nperp_max


      if(n_nperp.gt.n_nperp_a)then
        write(*,*)'in map_d_cold n_nperp>n_nperp_a'
        write(*,*)'n_nperp,n_nperp_a',n_nperp,n_nperp_a
        write(*,*)'Please increase n_nperp_a or decrease n_nperp'
        write(*,*)'and recompile the code'
        stop 'in map_d_cold'
      endif

      step=(nperp_max-nperp_min)/dfloat(n_nperp-1)

c      write(*,*)'mk_graph.f step',step
   
      do i=1,n_nperp
        nperp=nperp_min+(i-1)*step           
c        write(*,*)'i,nperp',i,nperp    
        call d_cold(z,r,phi,npar,nperp,disp_func)
        nperp_ar(i)=nperp
        d_ar(i)=dreal(disp_func)
c        write(*,*)'nperp_ar(i),d_ar(i)',nperp_ar(i),d_ar(i)
      enddo
      return
      end

      subroutine d_cold(z,r,phi,cnpar,cnper,disp_func)
c-----calculates cold plasma dispersion function disp_func
c     in space point (z,r,phi)
c      with refractive index components: cnpar,cnper
      implicit none
      include 'param.i'
      include 'one.i'
      include 'eps.i'
      double complex aK(3,3)
      double precision z,r,phi,cnpar,cnper
      double complex disp_func
c-----external
      double complex dcold_rlt
      integer i,j

      do i=1,3
        do j=1,3
           aK(3,3)=dcmplx(0.d0,0.d0)
        enddo
      enddo 
      call tensrcld(z,r,phi) 
      disp_func=dcold_rlt(reps,aK,cnpar,cnper)
        
      return
      end


      subroutine map_d_hot(z,r,phi,npar,n_nperp,nperp_min,nperp_max,
     &name_param,param,n_param)
c-----creates----------------------------------------------------------
c     1D array D(nperp) for plotting hot plasma dispersion function
c     Re(D(nperp)) on perpendicular refractive index N_perpendicular=nperp
c     in given the point (r,z,phi) and given N_parallel=npar
c
c     Plots ReD_hot(ReN_perp) to plot.ps file
c----------------------------------------------------------------------------
      implicit none
      include 'param.i'
      include 'one.i'
      include 'ions.i'
cfor test
      include 'eps.i'
c-----input:
      integer n_nperp                         !number of points in nperp mesh
      double precision nperp_min,nperp_max,   !minimal and maximal nperp values
     &z,r,phi, ! space point coordinates
     &npar     ! parallel refractive index
      integer n_param                         !number of input parameters
      character*(*) name_param(n_param)        !names  of the input parameters
      real param(n_param)                      !values of the input parameters    
c-----locals
      integer n_nperp_a !maximal number of n_nperp
c      parameter(n_nperp_a=10000)
      parameter(n_nperp_a=1000000)
      integer i

      double precision d_ar(n_nperp_a), !array of dispersion function values
     &nperp_ar(n_nperp_a)   !mesh of nperp

c      real d_ar(n_nperp_a), !array of dispersion function values
c     &nperp_ar(n_nperp_a)   !mesh of nperp

      double precision step,nperp
      double complex disp_func

      double precision x_ar(nbulka),y_ar(nbulka),
     &t_av_ar(nbulka),tpop_ar(nbulka),vflow_ar(nbulka),te

      double complex K_sum(3,3),d_complex

      integer iherm_loc,k1,k2


c-----externals
      double precision b,x,y,tempe,tpoprho,vflowrho
      double complex dhot_sum

      write(*,*)'in subroutine map_d_hot z,r,phi',z,r,phi
      write(*,*)'npar,n_nperp,nperp_min,nperp_max',
     &npar,n_nperp,nperp_min,nperp_max


      if(n_nperp.gt.n_nperp_a)then
        write(*,*)'in map_d_hot n_nperp>n_nperp_a'
        write(*,*)'n_nperp,n_nperp_a',n_nperp,n_nperp_a
        write(*,*)'Please increase n_nperp_a or decrease n_nperp'
        write(*,*)'and recompile the code'
        stop 'in map_d_hot'
      endif

      bmod=b(z,r,phi)

      do i=1,nbulk
        x_ar(i)=x(z,r,phi,i)
	y_ar(i)=y(z,r,phi,i)
        if(i.eq.1) y_ar(1)=-y_ar(1)
	te=tempe(z,r,phi,i) ! kev
	t_av_ar(i)=te*1000.d0      ! ev 
        tpop_ar(i)=tpoprho(rho,i)
        vflow_ar(i)=vflowrho(rho,i)
      enddo

c      write(*,*)'dmas',dmas
ctest      
      call tensrcld(z,r,phi) !tensor reps is in one.i  
cendtest
      step=(nperp_max-nperp_min)/dfloat(n_nperp-1)

      write(*,*)'n_nperp, step',n_nperp,step
      iherm_loc=1

      write(*,*)' nbulk,npar',nbulk,npar
      do i=1,nbulk
        write(*,*)'=i,dmas(i),x_ar(i),y_ar(i)',i,dmas(i),x_ar(i),y_ar(i)
        write(*,*)'t_av_ar(i),tpop_ar(i),vflow_ar(i)',
     &             t_av_ar(i),tpop_ar(i),vflow_ar(i)
      enddo
      do i=1,n_nperp
        nperp=nperp_min+(i-1)*step           
        nperp_ar(i)=nperp
        d_complex=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     &   vflow_ar,npar,nperp,iherm_loc,K_sum)
        d_ar(i)=dreal(d_complex)

      enddo
      return
      end

     
      double precision function freq_p(z,r,phi,i)
c--------------------------------------------------------
c     calculates plasma frequency [GHZ]
c-----------------------------------------------------
      implicit none
      include 'param.i'
      include 'one.i'
      include 'ions.i'
c--------------------------------------------
c     input
      double precision z,r,phi
      integer i !number of plasma component
c-----externals
      double precision dense

c-----locals
      double precision den,mass_e,charge_electron

      den=dense(z,r,phi,i) !10**13 cm**-3
      if(den.lt.0.d0)then
         den=0.d0
      endif

      mass_e=9.1094d-28         !electron rest mass (g)
      charge_electron=4.8032d-10 !electron charge (statcoulomb)  

c     greq_pe=sqrt(4*pi*n_e*e**2/m_e)/2pi      HZ
c     freq_pi=sqrt(4*pi*n_i*Z**2*e**2/m_i)/2pi      

      freq_p=dsqrt(4.d0*pi*den*1.d13*(charge_electron*charge(i))**2/
     &             (mass_e*dmas(i)))/(2.d0*pi)*1.d-9          !GHZ
     
    
      return
      end

       
      double precision function freq_c(z,r,phi,i)
c--------------------------------------------------------
c     calculates gyrofrequency [GHZ]
c-----------------------------------------------------
      implicit none
      include 'param.i'
      include 'one.i'
      include 'ions.i'
c--------------------------------------------
c     input
      double precision z,r,phi
      integer i !number of plasma component 

c-----locals
      double precision den,mass_e,charge_electron,clight

      mass_e=9.1094d-28          !electron rest mass (g)
      charge_electron=4.8032d-10 !electron charge (statcoulomb) 
      clight=2.99792458d10       !light speed (cm/sec)

c     freq_ce=(e*B/m_e*clight)/2pi      HZ
c     freq_ci=(Z*e*B/m_i*clihgt)/2pi      

      freq_c=charge_electron*charge(i)*bmod*1.d4/
     &       (mass_e*dmas(i)*clight)/(2.d0*pi)*1.d-9   

      return
      end
 
    
cyup not used
      subroutine plot_fcefuh(z_freq,r_freq,alpha_freq,
     &beta_freq,dist_freq,
     &nsteps_freq,n_ec_harmonics_freq,npar_freq,
     &max_plot_freq)
 
C     straightpath_plot_fcefuh to plot f_ce*j, f_uh, f_pe
C     with a path through the plasma (midplane for example)  

c      implicit none
      implicit integer (i-n), real*8 (a-h,o-z)

      include 'param.i' ! specifies code parameters 
c      include 'commons.i'
      include 'one.i'
      include 'three.i'
      include 'ions.i'
c-----input
      double precision
     & r_freq,z_freq, !straight line edge point [m]
     & dist_freq,     !straight line length [m]
     & alpha_freq,    !toroidal angle [degree] of straight line
     & beta_freq,     !angle between Z direction and straigth
                      !line direction [degree]
     &npar_freq,      !N_parallel to plot X mode cutoff
     &max_plot_freq   ! maximal frquency at plot GHZ

      integer nsteps_nfreq, !the number of points along the line
     &n_ec_harmonics_nfreq  !number of EC harmonic
c-----external
      double precision b,freq_c,freq_p
     
      integer maxnstep
      parameter (maxnstep=1000) !max value of nsteps_freq

      double precision raxis,zaxis,
     &fpi,fci,flh,fx1,fx2
      real ss(maxnstep),fuhs(maxnstep)
      real fces(maxnstep),rs(maxnstep)
      real fpes(maxnstep),rhos(maxnstep),den_e_s(maxnstep)
      real xmin,xmax,ymin,ymax,newx(maxnstep),newy(maxnstep)
     
      integer xtitlelen,titlelen,nchoice,j,k
      CHARACTER filename*40, line*1000, xtitle*100, titlestr*200
C Color index:
      INTEGER BLACK, WHITE, RED, GREEN, BLUE, CYAN, MAGENT, YELLOW
      PARAMETER (BLACK=0)
      PARAMETER (WHITE=1)
      PARAMETER (RED=2)
      PARAMETER (GREEN=3)
      PARAMETER (BLUE=4)
      PARAMETER (CYAN=5)
      PARAMETER (MAGENT=6)
      PARAMETER (YELLOW=7)
C Line style:
      INTEGER FULL, DASHED, DOTDSH, DOTTED, FANCY
      PARAMETER (FULL=1)
      PARAMETER (DASHED=2)
      PARAMETER (DOTDSH=3)
      PARAMETER (DOTTED=4)
      PARAMETER (FANCY=5)
C Character font:
      INTEGER NORMAL, ROMAN, ITALIC, SCRIPT
      PARAMETER (NORMAL=1)
      PARAMETER (ROMAN=2)
      PARAMETER (ITALIC=3)
      PARAMETER (SCRIPT=4)
C Fill-area style:
      INTEGER SOLID, HOLLOW
      PARAMETER (SOLID=1)
      PARAMETER (HOLLOW=2)

      if (nsteps_freq.gt.maxnstep) then
         write(*,*)'mk_graph.f in plot_fcefuh'
         write(*,*)'nsteps_freq.gt.maxnstep'
         write(*,*)'it should be nsteps_freq.le.maxnstep'
         write(*,*)'nsteps_freq,maxnstep',nsteps_freq,maxnstep
         write(*,*)'Please reduse nsteps_freq in gebray.dat'
         stop 'in plot_fcefuh'
      endif


C-------- find magnetic axis (psi=0, theta=0)
      call zr_psith(psimag,0.0,zaxis,raxis)
C**** FOR ARIESST: Magnetic axis at: R= 4.678357129999998  Z= -4.316799816857955E-18

c-----plot magnetic field contours to plot.ps file

C****** convert alpha, beta into radians from degrees
      alpha=alpha_freq*PI/180.0
      beta=beta_freq*PI/180.0
c      print *,'PI= ',PI
C****** starting coordinates in x,y,z
      xc=r_freq
      yc=0
      zc=z_freq
      sc=0
      ds=dist_freq/nsteps_freq
      dx=ds*sin(beta)*cos(alpha)
      dy=ds*sin(beta)*sin(alpha)
      dz=ds*cos(beta)
    
      do j=1,nsteps_freq 
C------- keep track of coordinates x,y,z, then convert to r,z,phi for
C------- calling bmod, dense, tempe
         r=sqrt(xc**2+yc**2)
         z=zc
         phi=atan2(yc,xc)
         bmod=b(z,r,phi)         
c         print *,'---r=',r,' z=',z,' phi=',phi,' rho=',rho,'--'
C----------- now calculate cold plasma stuff
         fces(j)=freq_c(z,r,phi,1)
         fpes(j)=freq_p(z,r,phi,1)
         fuhs(j)=sqrt(fces(j)**2+fpes(j)**2)
         ss(j)=sc
         rs(j)=r
         rhos(j)=rho
         den_e_s(j)=densrho(rho,1)
c         print *,'sc=',sc,'fces(j)=',fces(j),'fpes(j)=',fpes(j)
C**** step forward
         sc=sc+ds
         xc=xc+dx
         yc=yc+dy
         zc=zc+dz
      enddo

c      print *,'psimag=',psimag
c      print *,'Magnetic axis at: R=',raxis,' Z=',zaxis
    
c------------------------------------------------------------------
c     calulate minimal and maximal frequencies: ymin,ymax [GHZ]
c------------------------------------------------------------------
      ymin=fces(1) 
      ymax=fces(1)

      do j=1,nsteps_freq 
         if (fpes(j).lt.ymin) ymin=fpes(j)
         if (fces(j).lt.ymin) ymin=fces(j)
         if (fuhs(j).lt.ymin) ymin=fuhs(j)

         if (fpes(j).gt.ymax) ymax=fpes(j)
         if (fces(j)*n_ec_harmonics_freq.gt.ymax) then
             ymax=fces(j)*n_ec_harmonics_freq
         endif
         if (fuhs(j).gt.ymax) ymax=fuhs(j)
      enddo

C**** find limits for plot of real parts
C**** was ss before
        xmin=rs(nsteps_freq)
        xmax=rs(1)
C******* (reverse)
c        print *,'xmin=',xmin,' xmax=',xmax
C**** machine axis is at reqd (not magnetic axis with shafranov shift)
c        print *,'reqd: ',reqd
C***** Plot freq0, and dimensionless scale, and eta
        ymin=0.0      ! GHz for y axis
        ymax=max_plot_freq

C*** Sets limit for plot.
        xtitle='R along midplane (m)'
        xtitlelen=LEN_TRIM(xtitle)
        titlestr=''
        titlelen=LEN_TRIM(titlestr)

      OPEN(10,file='freqelec.bin',form='unformatted') 

      do j=1,nsteps_freq 
        write(10) real(rs(j)),real(fpes(j)),real(rhos(j)),
     &            real(den_e_s(j))
      enddo
      WRITE(10)
      do j=1,nsteps_freq 
        write(10) real(rs(j)),real(fces(j)),real(rhos(j)),
     &            real(den_e_s(j))
      enddo
      WRITE(10)
      do j=1,nsteps_freq 
        write(10) real(rs(j)),real(fuhs(j)),real(rhos(j)),
     &            real(den_e_s(j))
      enddo
      WRITE(10)
      do j=1,nsteps_freq 
        write(10) real(rs(j)),real(2*fces(j)),real(rhos(j)),
     &            real(den_e_s(j))
      enddo
      WRITE(10)
       do j=1,nsteps_freq 
        write(10) real(rs(j)),real(3*fces(j)),real(rhos(j)),
     &            real(den_e_s(j))
      enddo
      WRITE(10)
      do j=1,nsteps_freq 
        write(10) real(rs(j)),real(4*fces(j)),real(rhos(j)),  
     &            real(den_e_s(j))
      enddo   
      WRITE(10)
      do j=1,nsteps_freq 
        write(10) real(rs(j)),real(5*fces(j)),real(rhos(j)),
     &            real(den_e_s(j))
      enddo
      WRITE(10)
      do j=1,nsteps_freq 
        write(10) real(rs(j)),real(6*fces(j)),real(rhos(j)),
     &            real(den_e_s(j))
      enddo    
      WRITE(10)
      do j=1,nsteps_freq
        fx1=0.5d0*(-fces(j)+dsqrt(fces(j)**2+4.d0*fpes(j)**2/
     &                           (1.d0-npar_freq**2))) 
        write(10) real(rs(j)),real(fx1),real(rhos(j)),
     &            real(den_e_s(j))
c        write(*,*)'xe-(1+ye)*(1-npar_freq**2)'
c      write(*,*)'j,rs(j),fx1',j,rs(j),fx1
c      write(*,*)'fpes(j)/fx1)**2-(1.d0+fces(j)/fx1)*(1.d0-npar_freq**2)'
c      write(*,*)(fpes(j)/fx1)**2-(1.d0+fces(j)/fx1)*(1.d0-npar_freq**2)
      enddo 
      WRITE(10)
      do j=1,nsteps_freq
        fx2=0.5d0*(fces(j)+dsqrt(fces(j)**2+4.d0*fpes(j)**2/
     &                           (1.d0-npar_freq**2))) 
        write(10) real(rs(j)),real(fx2),real(rhos(j)),
     &            real(den_e_s(j))
c      write(*,*)'j,rs(j),fx2',j,rs(j),fx2
c      write(*,*)'fpes(j)/fx2)**2-(1.d0-fces(j)/fx2)*(1.d0-npar_freq**2)'
c      write(*,*)(fpes(j)/fx2)**2-(1.d0-fces(j)/fx2)*(1.d0-npar_freq**2)
      enddo 

      CLOSE(10)


      OPEN(10,file='freqion.bin',form='unformatted') 
c-----fpi for the first ion component
      do j=1,nsteps_freq 
        fpi=charge(2)*fpes(j)/dsqrt(dmas(2))
        write(10) real(rs(j)),real(fpi),real(rhos(j))
      enddo
      WRITE(10)
c-----fci for the first ion component
      do j=1,nsteps_freq
        fci= charge(2)*fces(j)/dmas(2)
        write(10) real(rs(j)),real(fci),real(rhos(j))
      enddo
      WRITE(10)
c-----lh
      do j=1,nsteps_freq 
        fpi=charge(2)*fpes(j)/dsqrt(dmas(2))
        fci= charge(2)*fces(j)/dmas(2)
        flh=dsqrt(fci**2+fpi**2*fces(j)**2/(fces(j)**2+fpes(j)**2))
        write(10) real(rs(j)),real(flh),real(rhos(j)) 
      enddo
      WRITE(10)

c-----2*fci for the first ion component
      do j=1,nsteps_freq
        fci=2*charge(2)*fces(j)/dmas(2)
        write(10) real(rs(j)),real(fci),real(rhos(j))  
      enddo

      WRITE(10)
c-----3*fci for the first ion component
      do j=1,nsteps_freq
        fci=3*charge(2)*fces(j)/dmas(2)
        write(10) real(rs(j)),real(fci),real(rhos(j))  
      enddo
      WRITE(10)

c-----4*fci for the first ion component
      do j=1,nsteps_freq
        fci=4*charge(2)*fces(j)/dmas(2)
        write(10) real(rs(j)),real(fci),real(rhos(j))  
      enddo
      WRITE(10)

      CLOSE(10)

      return
      end
c----------------subroutine mk_gronetwo1-------------*
c Prepares for file drawonet.in for xdraw            *
c INPUT: from common block  'onetwo.i'              *
c****************************************************
      subroutine mk_gronetwo_1
      implicit integer (i-n), real*8 (a-h,o-z)
      INTEGER I
      include 'param.i'
      include 'onetwo.i'

      OPEN(11,file='onetwo1.bin',form='unformatted')
      OPEN(21,file='onetwo1.dat')
20    format(11e12.3)

c      write(*,*)'mk_gronetwo_1 NR=',NR

      WRITE(21,*)'i rho powden spower powden_e powden_i powden_cl
     *curden_par curden_onetwo curden_tor cur_den_pol currden '

      hrho= rho_bin(2)-rho_bin(1) ! YuP: was: hrho=1.d0/(NR-1)
 
      DO 10 I=1,NR-1
         rho= rho_bin_center(i) ! YuP: was: hrho*(I-1+0.5d0)
         WRITE(11) REAL(rho),REAL(powden(I)),
     &   REAL(spower(I)),
     &   REAL(powden_e(I)),REAL(powden_i(I)),REAL(powden_cl(I)),   
     &   REAL(s_cur_den_parallel(I)),REAL(s_cur_den_onetwo(I)),
     &   REAL(s_cur_den_toroidal(I)),REAL(s_cur_den_poloidal(I)),
     &   REAL(currden(I))
 
         WRITE(21,20)rho,powden(I),
     &   spower(I),
     &   powden_e(I),powden_i(I),powden_cl(I),   
     &   s_cur_den_parallel(I),s_cur_den_onetwo(I),
     &   s_cur_den_toroidal(I),s_cur_den_poloidal(I),
     &   currden(I)
        
c         WRITE(*,*)i,rho,powden(I),spower(I),
c     &   powden_e(I),powden_i(I),powden_cl(I),   
c     &   s_cur_den_parallel(I),s_cur_den_onetwo(I),
c     &   s_cur_den_toroidal(I),s_cur_den_poloidal(I),
c     &   currden(I)

10    CONTINUE
      WRITE(11)
     
      CLOSE(11)
      CLOSE(21)

      write(*,*)'mk_graph.f: mk_gronetwo1 created onewto1.bin'
    
      RETURN
      END



      subroutine MK_graph_chamber_wall
c---------------------------------------------------
c     writes file wall.bin
c     with chamber wall coordinates
c----------------------------------------------------
      implicit none
      include 'param.i'
      include 'one_nml.i'
      include 'gr.i'
      include 'write.i'
      include 'fourb.i'
       
      integer i,j,m

      open(1,file='wall.bin',form='unformatted')
        
      DO J=1,3 
       DO I=1,NP+1
         WRITE(1) REAL(AR(J,I)),REAL(AZ(J,I))
       enddo
      enddo
      WRITE(1) 

      DO J=4,NL
       DO I=1,NP+1
         WRITE(1) REAL(AR(J,I)),REAL(AZ(J,I))
       enddo
      enddo
      WRITE(1) 

c      write(*,*)'mk_graph.f n_wall_add',n_wall_add

      if (n_wall.gt.1) then
         do i=1,n_wall  
           WRITE(1) REAL(r_wall(i)*100.0),REAL(z_wall(i)*100.0)
         enddo
         WRITE(1)
         do m=0,max_limiters          
            do i=1,n_wall_add(m)  
              WRITE(1) REAL(r_wall_add(i,m)*100.0),
     &                 REAL(z_wall_add(i,m)*100.0)
            enddo !i
            WRITE(1)
         enddo !m
      endif

      CLOSE(1)
      
      return
      end

      subroutine plot_cold_n_perp_omega_npar(z,r,phi,cnpar,
     &ratio_freq_min,ratio_freq_max,n_plot_freq)
c     
      implicit none

      include 'param.i'
      include 'one.i'
c-----input
      real*8 
     &z,r,phi,         ! space coordinates
     & ratio_freq_min, ! frequency_min/frqncy
     & ratio_freq_max,  ! frequency_max/frqncy
     & cnpar           ! N_parallel

      integer n_plot_freq !number of points of frequency mesh
c-----locals
      real*8 freq_loc,freq_min,freq_max,step,df, cnpar_loc,
     &eps,g,eta,det,yi,xi,a2,a0,n_z,n_x

      complex*16 im_one,det_x,det_y,det_z,e_x,e_y,e_z
      real*8     det_xy,det_yz,det_xz,det_max,emod,emod2

      integer j

      real*8, dimension(1:nbulk) :: v_loc,w_loc
      real*8, dimension(1:n_plot_freq) :: freq_ar,n_perp2_p,n_perp2_m
      real*8, dimension(1:n_plot_freq) :: cnpar_ar
      complex*16, dimension(1:n_plot_freq) :: Ex_p,Ey_p,Ez_p,E_long_p
      complex*16, dimension(1:n_plot_freq) :: Ex_m,Ey_m,Ez_m,E_long_m
      real*8, dimension(1:n_plot_freq) :: temp

      integer n,i,n_rung
c-----external
      real*8 b,x,y

      im_one=dcmplx(0.d9,1.d0)

      freq_min=ratio_freq_min*frqncy
      freq_max=ratio_freq_max*frqncy

      step=(freq_max-freq_min)/(n_plot_freq-1)

      bmod=b(z,r,phi)
      do i=1,nbulk
	  v_loc(i)=v(i)
	  w_loc(i)=w(i)
      enddo

      do n=1,n_plot_freq

         freq_loc=  freq_min+step*(n-1)  
         freq_ar(n)= freq_loc
         df=frqncy/freq_loc 

         do i=1,nbulk
       	    v(i)=v_loc(i)*df* df
	    w(i)=w_loc(i)*df
         enddo !i
       
c         cnpar_loc=cnpar*df
         cnpar_loc=cnpar
         cnpar_ar(n)=cnpar_loc

         eps=1.d0
         g=0.d0
         eta=1.d0

         do i=1,nbulk
            xi=x(z,r,phi,i)
            yi=y(z,r,phi,i)

       	    eps=eps-xi/(1.d0-yi**2)
            eta=eta-xi

            if (i.eq.1) then
               g=g+xi*yi/(1.d0-yi**2)
            else
	       g=g-xi*yi/(1.d0-yi**2)
            endif  
         enddo !i

c--------dispersion relation
c        eps*N_perp**4-a2*N_perp**2+a0=0

         a2=(eps+eta)*(eps-cnpar_loc**2)-g*g
         a0=eta*((eps-cnpar_loc**2)**2-g*g)

         det=a2*a2-4.d0*eps*a0
  
         if (det.ge.0) then
           n_perp2_p(n)=(a2+dsqrt(det))/(2.d0*eps)
           n_perp2_m(n)=(a2-dsqrt(det))/(2.d0*eps)
         else
           n_perp2_p(n)=0.d0
           n_perp2_m(n)=0.d0
         endif

c        write(*,*)'n,eps*n_perp2_p(n)**2-a2*n_perp2_p(n)+a0',
c     &              n,eps*n_perp2_p(n)**2-a2*n_perp2_p(n)+a0
c        write(*,*)'n_perp2_p(n)',n_perp2_p(n)

cc         write(*,*)'n,eps*n_perp2_m(n)**2-a2*n_perp2_m(n)+a0',
cc     &              n,eps*n_perp2_m(n)**2-a2*n_perp2_m(n)+a0

c------------------------------------------------------------
c        polarizatiion
c
c        (eps-N_z^2)E_x + igE_y                + N_xN_zE_z      =0
c
c         - igE_x       + (eps-N_x^2-N_z^2)E_y                  =0
c
c         N_zN_xE_x     +                        (eta-N_x^2)E_z =0
c
c-------------------------------------------------------------------
         
          if (n_perp2_p(n).ge.0) then
             n_z=cnpar_loc
             n_x=dsqrt(n_perp2_p(n))

             det_xy=(eps-n_z**2)*(eps-n_x**2-n_z**2)-g**2
             det_yz=(eps-n_x**2-n_z**2)*(eta-n_x**2)    
             det_xz=(eps-n_z**2)*(eta-n_x**2)-n_x**2*n_z**2    
            
             det_max=dabs(det_yz)
             j=1
             if(dabs(det_yz).gt.det_max) then
                det_max=dabs(det_yz)
                j=2
             endif
             if(dabs(det_xz).gt.det_max) then
                det_max=dabs(det_xz)
                j=3
             endif

             if (det_max.gt.0) then
                n_rung=2
             else
                n_rung=1
             endif  

             if(j.eq.1) then
               e_z=1.d0
               det_x=-n_z*n_x*(eps-n_x**2-n_z**2)
               det_y= im_one*g*n_z*n_x
               e_x=det_x/det_xy
               e_y=det_y/det_xy
             endif

             if(j.eq.2) then
               e_x=1.d0
               det_y=im_one*g*(eta-n_x**2)
               det_z=-(eps-n_x**2-n_z**2)*n_x*n_z 
               e_y=det_y/det_yz
               e_z=det_z/det_yz
             endif

             if(j.eq.3) then
               e_y=1.d0
               det_x=im_one*g*(eta-n_x**2)
               det_z=im_one*g*n_x*n_z 
               e_x=det_x/det_xz
               e_z=det_z/det_xz
             endif
             emod2=dreal(e_x*dconjg(e_x)+e_y*dconjg(e_y)+
     &                   e_z*dconjg(e_z))
             emod=dsqrt(emod2)
             Ex_p(n)=e_x/emod
             Ey_p(n)=e_y/emod
             Ez_p(n)=e_z/emod
             E_long_p(n)=(Ez_p(n)*cnpar_loc+Ex_p(n)*n_perp2_p(n))/
     &                   dsqrt(cnpar_loc**2+n_perp2_p(n)**2)
          else
             Ex_p(n)=(0.0,0.0)
             Ey_p(n)=(0.0,0.0)
             Ez_p(n)=(0.0,0.0)
             E_long_p(n)=(0.0,0.0)
          endif
c          write(*,*)'Ex N_p', Ex_p(n)
c          write(*,*)'Ey N_p', Ey_p(n)
c          write(*,*)'Ez N_p', Ez_p(n)
c          write(*,*)'E_long N_p', E_long_p(n)
c-----------------------------------------      
c          write(*,*)'n_perp2_m(n)',n_perp2_m(n)

          if (n_perp2_m(n).ge.0) then
             n_z=cnpar_loc
             n_x=dsqrt(n_perp2_m(n))

             det_xy=(eps-n_z**2)*(eps-n_x**2-n_z**2)-g**2
             det_yz=(eps-n_x**2-n_z**2)*(eta-n_x**2)    
             det_xz=(eps-n_z**2)*(eta-n_x**2)-n_x**2*n_z**2    
            
             det_max=dabs(det_yz)
           
             j=1
             if(dabs(det_yz).gt.det_max) then
                det_max=dabs(det_yz)
                j=2
             endif
             if(dabs(det_xz).gt.det_max) then
                det_max=dabs(det_xz)
                j=3
             endif

c             write(*,*)'det_max,j',det_max,j

             if (det_max.gt.0) then
                n_rung=2
             else
                n_rung=1
             endif  

             if(j.eq.1) then
               e_z=1.d0
               det_x=-n_z*n_x*(eps-n_x**2-n_z**2)
               det_y= im_one*g*n_z*n_x
               e_x=det_x/det_xy
               e_y=det_y/det_xy
c               write(*,*)'j,det_xy,e_z,e_x,e_y',j,det_xy,e_z,e_x,e_y
             endif

             if(j.eq.2) then
               e_x=1.d0
               det_y=im_one*g*(eta-n_x**2)
               det_z=-(eps-n_x**2-n_z**2)*n_x*n_z 
               e_y=det_y/det_yz
               e_z=det_z/det_yz
c               write(*,*)'j,det_yz,e_z,e_x,e_y',j,det_yz,e_z,e_x,e_y
             endif

             if(j.eq.3) then
               e_y=1.d0
               det_x=im_one*g*(eta-n_x**2)
               det_z=im_one*g*n_x*n_z 
               e_x=det_x/det_xz
               e_z=det_z/det_xz
c               write(*,*)'j,det_xz,e_z,e_x,e_y',j,det_xz,e_z,e_x,e_y
             endif
             emod2=dreal(e_x*dconjg(e_x)+e_y*dconjg(e_y)+
     &                   e_z*dconjg(e_z))
             emod=dsqrt(emod2)
             Ex_m(n)=e_x/emod
             Ey_m(n)=e_y/emod
             Ez_m(n)=e_z/emod
             E_long_m(n)=(Ez_m(n)*cnpar_loc+Ex_m(n)*n_perp2_m(n))/
     &                   dsqrt(cnpar_loc**2+n_perp2_m(n)**2)
          else
             Ex_m(n)=(0.0,0.0)
             Ey_m(n)=(0.0,0.0)
             Ez_m(n)=(0.0,0.0)
             E_long_m(n)=(0.0,0.0)
          endif
c          write(*,*)'Ex N_m', Ex_m(n)
c          write(*,*)'Ey N_m', Ey_m(n)
c          write(*,*)'Ez N_m', Ez_m(n)
c          write(*,*)'E_long N_m', E_long_m(n)
c---------------------------------------------------     
      enddo !n

      do i=1,nbulk
       	 v(i)=v_loc(i)
	 w(i)=w_loc(i)
      enddo !i   




      return
      end

      

