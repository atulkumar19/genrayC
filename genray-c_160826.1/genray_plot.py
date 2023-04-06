# genray_plot.py'
# Plots genray.nc (output data file produced by GENRAY)
# Yuri Petrov   CompX   2011

from numpy import *
from mpl_toolkits.mplot3d import Axes3D

from pylab import *
from matplotlib import rc 
from matplotlib.pyplot import cm, figure, axes, plot, xlabel, ylabel,  \
     title, savefig, show

from netCDF4 import Dataset

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import time
import pylab as pylab
#matplotlib.interactive(True) # no plots on screen
matplotlib.interactive(False) # with plots on screen

# Specify for plots:
rhomax=5.0  # In the code, rhomax is defined as max of rho 
            # over all equilibrium grid - it could be very large.
            # If the power deposition happened within plasma (rho<0)
            # the power profiles may look "compressed".
            # For clarity, set the limit of rho for plots.
            # If the whole rho range is desired, set rhomax to 0.
            # it will be determined automatically.             
            
            
fnt  = 10     # font size for axis numbers (see 'param=' below) 
linw = 1.0    # LineWidth for contour plots
Ncont= 100 #50     # Number of contour levels for PSI (pol.flux)

# For plots of Dispersion Function vs x, showing Hot plasma roots
# and also plots of cold Nperp vs x
plot_Nperp_prof=1 # 0- do not plot; 1-plot (takes some time)

pi=3.14159265358979323846264338327950288419716939937510

arr_len=3.  # Specify arrow length for plots of (Nr,Nz) at start point

e0 = time.time()  # elapsed time since the epoch
c0 = time.clock() # total cpu time spent in the script so far


filenm='genray.nc'
# Open the genray netcdf file, read only
dat= Dataset(filenm, 'r', format='NETCDF4')


print dat.file_format # print which format was used for the genray.nc file

print 'The genray file, ',filenm,', contains:'
print '========================================'
print "The global attributes: ",dat.dimensions.keys()        
print "File contains variables: ",dat.variables.keys()
print '========================================'

print '---------------------------------------GENERAL PROFILES'
# some eqdsk data:
eqdsk_x=dat.variables['eqdsk_x']
print 'eqdsk_x is: ', eqdsk_x.long_name, eqdsk_x.shape
eqdsk_y=dat.variables['eqdsk_y']
print 'eqdsk_y is: ', eqdsk_y.long_name, eqdsk_y.shape
eqdsk_r=dat.variables['eqdsk_r']
print 'eqdsk_r is: ', eqdsk_r.long_name, eqdsk_r.shape
eqdsk_z=dat.variables['eqdsk_z']
print 'eqdsk_z is: ', eqdsk_z.long_name, eqdsk_z.shape
eqdsk_psi=dat.variables['eqdsk_psi']
print 'eqdsk_psi is: ', eqdsk_psi.long_name, eqdsk_psi.shape
nxeqd=np.size(eqdsk_x)
nzeqd=np.size(eqdsk_z)
nreqd=np.size(eqdsk_r)
print 'nxeqd,nzeqd,nreqd=',nxeqd,nzeqd,nreqd

model_b=dat.variables['model_b'].getValue()  #getValue() for scalar
model_b=np.asscalar(model_b)
model_rho_dens=dat.variables['model_rho_dens'].getValue()  #getValue() for scalar
model_rho_dens=np.asscalar(model_rho_dens)
if model_b==4: # FRC
    rs_frc=dat.variables['rs_frc'].getValue()  #getValue() for scalar
    rs_frc= np.asscalar(rs_frc)*100 # cm FRC separatrix
    r0_frc= rs_frc/sqrt(2)  #   cm
        
rscan=dat.variables['xscan'] # it was 'rscan' before 09-30-2014
rscan=rscan[:]*100  # m->cm
rhoscan=dat.variables['rhoscan']  #  rho(xscan)
wcw=dat.variables['wcwscan']  # wce/w as a func of rscan
wpw=dat.variables['wpwscan']  # wpe/w as a func of rscan
wuw=dat.variables['wuwscan']  # wuh/w as a func of rscan (Upper_hybrid/omega)
# Values of Y and Z at which the scan is done:
yscan=dat.variables['y_save_disp'] # a single value in the code
zscan=dat.variables['z_save_disp']
yscan=np.asscalar(yscan[:])*100 # cm
zscan=np.asscalar(zscan[:])*100 # cm
rscan_size=np.size(rscan)

# Note: r_wall array may not exist in genray.nc:
try:
    try:
        r_wall=dat.variables['r_wall']
    except:
        print('No data on r_wall, z_wall in genray.nc')
        n_wall=0
    else:
        r_wall=dat.variables['r_wall']
        z_wall=dat.variables['z_wall']
        n_wall=r_wall.size
        print 'r_wall is: ', r_wall.long_name, r_wall.shape
        print 'z_wall is: ', z_wall.long_name, z_wall.shape  
        print 'n_wall=',n_wall
        r_wall=np.asarray(r_wall)*100
        z_wall=np.asarray(z_wall)*100
finally:
    print '----------------------------------------'

if plot_Nperp_prof==1:
    try:
        try:
            Nperp2m=dat.variables['Nperp2m']
        except:
            print('No data on Nperp2(r) in genray.nc')
            plot_Nperp_prof=0 # reset to 0, meaning: no plots of Disp.Func.
        else:
            Npara=dat.variables['Npara']
            Npera=dat.variables['Npera']
            Nperp2m=dat.variables['Nperp2m']
            Nperp2p=dat.variables['Nperp2p']
            Nperp2hot=dat.variables['Nperp2hot']
            ddd=dat.variables['ddd']
            inpar=np.size(Nperp2m,0)
            inper=np.size(ddd,1)
            print 'inpar=',inpar
            print 'rscan_size=',rscan_size
            print 'inper=',inper
    finally:
        print '----------------------------------------'    
    

freqcy=dat.variables['freqcy']
print 'freqcy =',freqcy.long_name, freqcy[:], freqcy.units
freqcy=freqcy.getValue()  #getValue() for scalar

mass=dat.variables['dmas']
print 'mass=', mass.long_name, mass[:], mass.units

charge=dat.variables['charge']
print 'charge=', charge.long_name, charge[:], charge.units

rho_bin_center=dat.variables['rho_bin_center']
print 'rho_bin_center is: ', rho_bin_center.long_name, rho_bin_center.shape
#print rho_bin_center[:]

rho_bin=dat.variables['rho_bin']
print 'rho_bin is: ', rho_bin.long_name, rho_bin.shape

binvol=dat.variables['binvol']
print 'binvol is: ', binvol.long_name, binvol.shape

densprof=dat.variables['densprof']
print 'densprof is: ', densprof.long_name, densprof.shape

temprof=dat.variables['temprof']
print 'temprof is: ', temprof.long_name, temprof.shape

Nsp=temprof[:,0].size  # number of species
print 'Number of species: Nsp=',Nsp


w_dens_vs_x_nc=dat.variables['w_dens_vs_x_nc']
print 'w_dens_vs_x_nc is: ', w_dens_vs_x_nc.long_name, w_dens_vs_x_nc.shape
w_temp_vs_x_nc=dat.variables['w_temp_vs_x_nc']
print 'w_temp_vs_x_nc is: ', w_temp_vs_x_nc.long_name, w_temp_vs_x_nc.shape
w_x_densprof_nc=dat.variables['w_x_densprof_nc']
print 'w_x_densprof_nc is: ', w_x_densprof_nc.long_name, w_x_densprof_nc.shape

w_dens_vs_y_nc=dat.variables['w_dens_vs_y_nc']
print 'w_dens_vs_y_nc is: ', w_dens_vs_y_nc.long_name, w_dens_vs_y_nc.shape
w_temp_vs_y_nc=dat.variables['w_temp_vs_y_nc']
print 'w_temp_vs_y_nc is: ', w_temp_vs_y_nc.long_name, w_temp_vs_y_nc.shape
w_y_densprof_nc=dat.variables['w_y_densprof_nc']
print 'w_y_densprof_nc is: ', w_y_densprof_nc.long_name, w_y_densprof_nc.shape

# Added 01-12-2016
bmodprofxz=dat.variables['bmodprofxz']
print 'bmodprofxz is: ', bmodprofxz.long_name, bmodprofxz.shape
# Can be plotted vs eqdsk_x and eqdsk_z grids generated by Genray-C

# Density profiles on xy, xz and yz-grids
densprofxy=dat.variables['densprofxy']
print 'densprofxy is: ', densprofxy.long_name, densprofxy.shape
densprofxz=dat.variables['densprofxz']
print 'densprofxz is: ', densprofxz.long_name, densprofxz.shape
densprofyz=dat.variables['densprofyz']
print 'densprofyz is: ', densprofyz.long_name, densprofyz.shape


# O-X modes Transmission data: Only present for i_ox=2:
try:
    try:
        i_ox_conversion=dat.variables['i_ox_conversion']
    except:
        print('No data on O-X conversion in genray.nc')
        i_ox=0
    else:
        i_ox=2
        i_ox_conversion=dat.variables['i_ox_conversion'] #
        print 'i_ox_conversion is:',i_ox_conversion.long_name
        transm_ox=dat.variables['transm_ox'] #
        print 'transm_ox is: ', transm_ox.long_name, transm_ox.shape
        cnpar_ox=dat.variables['cnpar_ox'] #
        print 'cnpar_ox is: ', cnpar_ox.long_name, cnpar_ox.shape
        cn_b_gradpsi=dat.variables['cn_b_gradpsi'] #
        print 'cn_b_gradpsi is: ', cn_b_gradpsi.long_name, cn_b_gradpsi.shape
finally:
    print '----------------------------------------'
    


#[Nov-2014] Added: power deposition profiles over (R,Z) rectangular grid.
#  [erg/sec]
Rgrid=dat.variables['Rgrid']
print 'Rgrid is: ', Rgrid.long_name, Rgrid.shape
Zgrid=dat.variables['Zgrid']
print 'Rgrid is: ', Zgrid.long_name, Zgrid.shape
NRgrid=Rgrid[:].size  
NZgrid=Zgrid[:].size  
spwr_rz_e=dat.variables['spwr_rz_e']
print 'spwr_rz_e is: ', spwr_rz_e.long_name #, spwr_rz_e.shape
spwr_rz_i=dat.variables['spwr_rz_i']
print 'spwr_rz_i is: ', spwr_rz_i.long_name #, spwr_rz_i.shape
spwr_rz_cl=dat.variables['spwr_rz_cl']
print 'spwr_rz_cl is: ', spwr_rz_cl.long_name #, spwr_rz_cl.shape

Rgrid= Rgrid[:]*100  # m->cm
Zgrid= Zgrid[:]*100  # m->cm
spwr_rz_e=  spwr_rz_e[:]/1.0e10  # kW
spwr_rz_i=  spwr_rz_i[:]/1.0e10  # kW
spwr_rz_cl= spwr_rz_cl[:]/1.0e10  # kW
sum_spwr_rz_e= sum(spwr_rz_e)
sum_spwr_rz_i= sum(spwr_rz_i)
sum_spwr_rz_cl=sum(spwr_rz_cl)
print 'sum(spwr_rz_e)=    ',sum_spwr_rz_e,   ' kW'
print 'sum(spwr_rz_i)=    ',sum_spwr_rz_i,   ' kW'
print 'sum(spwr_rz_cl)=   ',sum_spwr_rz_cl,  ' kW'

# Define min/max for major radius range:
if n_wall==0:  # wall is not defined; get limits from eqdsk:
    Rmax=np.amax(eqdsk_r)*100
    Rmin=np.amin(eqdsk_r)*100 # meters -> cm
else:
    Rmax=np.amax(r_wall)
    Rmin=np.amin(r_wall)
    Zmin=np.amin(z_wall)
    Zmax=np.amax(z_wall)
    Rmax=max(Rmax,abs(Rmin))
    Rmin=-Rmax
    Zmax=min(Zmax,Zgrid[NZgrid-1])
    Zmin=max(Zmin,Zgrid[0])
    
print 'Rmin,Rmax=', Rmin,Rmax , 'Zmin,Zmax=', Zmin,Zmax    


    
# Define boundary for top-view (R-phi) plots:
phi=arange(301.)*2*pi/300.  # tor.angle in rad.
bndry_in_X= Rmin*cos(phi)
bndry_in_Y= Rmin*sin(phi)
bndry_out_X= Rmax*cos(phi)
bndry_out_Y= Rmax*sin(phi)


#===================  PLOTS =============================================
#set fonts and line thicknesses
params = {
    'axes.linewidth': linw,
    'lines.linewidth': linw,
    'axes.labelsize': fnt+4,
    'font.size': fnt+4,
    'legend.fontsize': fnt,
    'xtick.labelsize':fnt,
    'ytick.labelsize':fnt,
    'font.weight'  : 'regular'
}

mpl.rcParams.update(params)
#rc.defaults() #to restore defaults

mpl.rcParams['font.size']=fnt+2  # set font size for text in mesh-plots



#--------------------------------------------------------------------------
fig0=plt.figure()  # Te and ne profiles vs x or y coord.
plt.subplot(221)
plt.hold(True)
plt.ylabel('$ T $   $ (keV) $')
plt.grid(True)
T_max= 1.05*np.amax(w_temp_vs_x_nc[0,:])+0.001  # upper limit for plots
for i in range(0,Nsp,1):
    if remainder(i,4)==0: col='k'
    if remainder(i,4)==1: col='b'
    if remainder(i,4)==2: col='g'
    if remainder(i,4)==3: col='r'    
    #plt.plot(rho_bin[:],temprof[i,:],color=col,linewidth=linw*(Nsp-i) )
    plt.plot(w_x_densprof_nc[:]*100,w_temp_vs_x_nc[i,:],color=col,linewidth=linw*(Nsp-i) )
if model_b==4: # FRC
    plt.plot([+rs_frc, +rs_frc],[0, T_max],'r--')
    plt.plot([+r0_frc, +r0_frc],[0, T_max],'r--')
    plt.plot([-rs_frc, -rs_frc],[0, T_max],'r--')
    plt.plot([-r0_frc, -r0_frc],[0, T_max],'r--')
    plt.title('$Dashed:$ $B=0$ $and$ $Separatrix$')
axis([-Rmax,Rmax,0.,T_max])  

plt.subplot(223)
plt.hold(True)
plt.xlabel('$ x $   $ (cm) $')
#plt.xlabel(r'$\rho$')
plt.ylabel('$ n $   $ (cm^{-3}) $')
plt.grid(True)
n_max= 1.05*np.amax(w_dens_vs_x_nc[0,:])+1 # upper limit for plots
for i in range(0,Nsp,1):
    if remainder(i,4)==0: col='k'
    if remainder(i,4)==1: col='b'
    if remainder(i,4)==2: col='g'
    if remainder(i,4)==3: col='r'    
    #plt.plot(rho_bin[:],densprof[i,:],color=col,linewidth=linw*(Nsp-i) )   
    plt.plot(w_x_densprof_nc[:]*100,w_dens_vs_x_nc[i,:],color=col,linewidth=linw*(Nsp-i) )
if model_b==4: # FRC
    plt.plot([+rs_frc, +rs_frc],[0, n_max],'r--')
    plt.plot([+r0_frc, +r0_frc],[0, n_max],'r--')
    plt.plot([-rs_frc, -rs_frc],[0, n_max],'r--')
    plt.plot([-r0_frc, -r0_frc],[0, n_max],'r--')
axis([-Rmax,Rmax,0.,n_max])  

plt.subplot(222)
plt.hold(True)
#plt.ylabel('$ T $   $ (keV) $')
plt.grid(True)
for i in range(0,Nsp,1):
    if remainder(i,4)==0: col='k'
    if remainder(i,4)==1: col='b'
    if remainder(i,4)==2: col='g'
    if remainder(i,4)==3: col='r'    
    plt.plot(w_y_densprof_nc[:]*100,w_temp_vs_y_nc[i,:],color=col,linewidth=linw*(Nsp-i) )
axis([-Rmax,Rmax,0.,1.05*np.amax(w_temp_vs_y_nc[0,:])+0.001])  

plt.subplot(224)
plt.hold(True)
plt.xlabel('$ y $   $ (cm) $')
#plt.ylabel('$ n $   $ (cm^{-3}) $')
plt.grid(True)
for i in range(0,Nsp,1):
    if remainder(i,4)==0: col='k'
    if remainder(i,4)==1: col='b'
    if remainder(i,4)==2: col='g'
    if remainder(i,4)==3: col='r'    
    plt.plot(w_y_densprof_nc[:]*100,w_dens_vs_y_nc[i,:],color=col,linewidth=linw*(Nsp-i))   
axis([-Rmax,Rmax,0.,1.05*np.amax(w_dens_vs_y_nc[0,:])+1])    
savefig('genray_profiles_T-n.png',format='png') # try pdf,eps,ps
plt.show() 


#-----------------------------------------------------------------------------
fig0=plt.figure()  # B vs Z coord. at x=0,y=0 axis
ix0= int(nxeqd/2) #This is the index for x=0       ix0=530
x00= eqdsk_x[ix0] #value of X will be printed in plot (may be not exactly 0)
Baxial=bmodprofxz[:,ix0]
f=freqcy
plt.subplot(211) #------------------------- B(Z)
plt.hold(True)
plt.grid(True)
plt.title('$Z-scan$  $at$ $y=0,$ $x=$' +"%5.3f" %(x00) +'$cm$' )
plt.ylabel('$B(Z)$  $(T)$' )
#plt.xlabel('$ Z $   $ (cm) $')
plt.plot(eqdsk_z[:]*100,Baxial[:],'b',linewidth=1.5)
#--------------
plt.subplot(212) #------------------------- omega/omega_c vs Z
plt.hold(True)
plt.grid(True)
if Nsp>1:
    isp=1 # ion species #1
    mime= mass[isp]/mass[0]      # mi/me ratio
    Zi= charge[isp]/charge[0]    # Zi/e
    plt.ylabel('$\omega/\omega_{ci}$  $for$ $m/m_e=$'+r"$%1.0f$" %(mime) +' $Z=$'+r"$%1.0f$" %(Zi))
else:  # electrons only: 
    isp=0 # electrons
    mime=1.0
    Zi=1.0
    plt.ylabel('$\omega/\omega_{ce}$')
plt.xlabel('$Z$   $ (cm) $')
plt.plot(eqdsk_z[:]*100,f/(28e9*Baxial[:]*Zi/mime),'b',linewidth=1.5)       
#--------------
savefig('genray_B_vs_Z.png',format='png')
plt.show() 

#-----------------------------------------------------------------------------
fig0=plt.figure()  # B vs X coord. at y=0 z=0 plane
iz0= int(nzeqd/2) #This is the index for z=0 (approximately)  
z00= eqdsk_z[iz0] # The exact value of Z (will be printed in plot)
BvsX=bmodprofxz[iz0,:]
BvsX_min= np.amin(BvsX[:])
BvsX_max= np.amax(BvsX[:])
BvsX_lim= max(BvsX_min,BvsX_max)
f=freqcy
plt.subplot(221) #------------------------- B(X)
plt.hold(True)
plt.grid(True)
plt.title('$X-scan$  $at$ $y=0,$ $z=$' +"%5.3f" %(z00) +'$cm$' )
plt.ylabel('$B(X)$  $(T)$' )
#plt.xlabel('$ X $   $ (cm) $')
if model_b !=4:  # not equal to 4
    #Not FRC: Plot the magnitude of B (as saved in bmodprofxz)
    plt.plot(eqdsk_x[:]*100,BvsX[:],'b',linewidth=1.5)
if model_b==4: # FRC
    # Reverse sign of B at x=r0_frc and plot the position of separatrix and B=0
    # set the vertical limits for this plot:
    #ylim( (-BvsX_lim, BvsX_lim) ) 
    # Find indexes corresponding to |X|<r0_frc (null point):
    ir0= np.where( (eqdsk_x<r0_frc/100.) & (eqdsk_x>-r0_frc/100.) )
    #print 'ir0=',ir0
    # Find indexes corresponding to |X|>r0_frc :
    irn= np.where( (eqdsk_x<-r0_frc/100.) ) # where X<-r0
    irp= np.where( (eqdsk_x> r0_frc/100.) ) # where X>+r0
    #print 'irs=',irs
    # Plot B inside null points as negative:
    plt.plot(eqdsk_x[ir0]*100,-BvsX[ir0],'b',linewidth=1.5)
    # and outside of null points as positive:
    plt.plot(eqdsk_x[irn]*100,+BvsX[irn],'b',linewidth=1.5)
    plt.plot(eqdsk_x[irp]*100,+BvsX[irp],'b',linewidth=1.5)  
    ## Plot B inside null points as positive:
    # plt.plot(eqdsk_x[ir0]*100,+BvsX[ir0],'b',linewidth=1.5)
    ## and outside of null points as negative:
    # plt.plot(eqdsk_x[irn]*100,-BvsX[irn],'b',linewidth=1.5)
    # plt.plot(eqdsk_x[irp]*100,-BvsX[irp],'b',linewidth=1.5)    
    # Designate null point and separatrix point by dashed lines:
    plt.plot([+rs_frc, +rs_frc],[-BvsX_lim, BvsX_lim],'r--')
    plt.plot([+r0_frc, +r0_frc],[-BvsX_lim, BvsX_lim],'r--')
    plt.plot([-rs_frc, -rs_frc],[-BvsX_lim, BvsX_lim],'r--')
    plt.plot([-r0_frc, -r0_frc],[-BvsX_lim, BvsX_lim],'r--')
    # Also plot the horizontal line for B=0 level:
    plt.plot(eqdsk_x[:]*100,BvsX[:]*0,'k',linewidth=1.)
    #--------------
plt.subplot(223) #------------------------- omega/omega_c vs X
plt.hold(True)
plt.grid(True)
if Nsp>1:
    isp=1 # ion species #1
    mime= mass[isp]/mass[0]      # mi/me ratio
    Zi= charge[isp]/charge[0]    # Zi/e
    plt.ylabel('$\omega/\omega_{ci}$  $for$ $m/m_e=$'+r"$%1.0f$" %(mime) +' $Z=$'+r"$%1.0f$" %(Zi))
else:  # electrons only: 
    isp=0 # electrons
    mime=1.0
    Zi=1.0
    plt.ylabel('$\omega/\omega_{ce}$')
plt.xlabel('$X$   $ (cm) $')
plt.yscale('log') # log 10 scale !!! Comment it, if you want a linear scale
plt.plot(eqdsk_x[:]*100,f/(28e9*BvsX[:]*Zi/mime),'b',linewidth=1.5)       
#--------------
savefig('genray_B_vs_X.png',format='png')
plt.show() 


#--------------------------------------------------------------------------

fig0=plt.figure()  # omega_ce/omega profiles vs x coord.
f=freqcy  #[0]
plt.subplot(221) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title( r"$f (GHz) =$"+r"$%1.2f$" %(f*1e-9) )
plt.ylabel('$\omega_{ce}/\omega$ $,$ $\omega_{pe}/\omega$ $,$ $\omega_{UH}/\omega$')
#plt.xlabel('$ x $   $ (cm) $')
plt.grid(True)
plt.plot(rscan[:],wcw[:],'b',linewidth=1.5) # fce/f
plt.plot(rscan[:],wpw[:],'r',linewidth=1.5) # fpe/f
plt.plot(rscan[:],wuw[:],'k',linewidth=1.5) # fuh/f
wcw_max=np.amax(wcw[:])
wpw_max=np.amax(wpw[:])
wpw_min=np.amin(wpw[:])
print 'max of wpe/w =',wpw_max, '  (along Z[cm]=',zscan,')'
print 'max of wce/w =',wcw_max, '  (along Z[cm]=',zscan,')'
www_max=max(3*wcw_max,3*wpw_min)
plt.plot([0, Rmax],[0, 0],'k')
if model_b==4: # FRC
    plt.plot([rs_frc, rs_frc],[0, www_max],'r--')
    plt.plot([r0_frc, r0_frc],[0, www_max],'r--')
xlim( (0., Rmax) )  
# set the upper limit for this plot:
ylim( (0., www_max) ) 
#--------------
plt.subplot(222) #-------------------------
plt.hold(True)
plt.grid(True)
plt.ylabel('$\omega/\omega_{ce}$ $,$ $\omega/\omega_{pe}$ $,$ $\omega/\omega_{UH}$')
#plt.xlabel('$ x $   $ (cm) $')
plt.grid(True)
plt.title('$Scan$  $at$ $(y,z)=($' +"%5.0f" %(yscan) +',' +\
"%5.0f" %(zscan) +'$)cm$' )
plt.plot(rscan[:],1./wcw[:],'b',linewidth=1.5) # w/wce
plt.plot(rscan[:],1./wpw[:],'r',linewidth=1.5) # w/wpe
plt.plot(rscan[:],1./wuw[:],'k',linewidth=1.5) # w/wuh
plt.plot([0, Rmax],[0, 0],'k')
owcw_edge=np.amax(1./wcw[rscan_size-1])   
owcw_max= ceil(owcw_edge*3) # set upper limit to 3*(omega/omega_c)_edge
if model_b==4: # FRC
    plt.plot([rs_frc, rs_frc],[0, owcw_max],'r--')
    plt.plot([r0_frc, r0_frc],[0, owcw_max],'r--')
xlim( (0., Rmax) )  
# To avoid large peaks where B->0,  set the upper limit for this plot:
ylim( (0., owcw_max) ) 
#--------------
plt.subplot(223) #-------------------------
plt.hold(True)
plt.grid(True)
plt.ylabel(r'$\rho$ $coordinate$')
plt.xlabel('$ x $   $ (cm) $')
plt.grid(True)
plt.plot(rscan[:],rhoscan[:],'b',linewidth=1.5)  
plt.plot([0, Rmax],[0, 0],'k')
rho_max= rhoscan[rscan_size-1]
if model_b==4: # FRC
    plt.plot([rs_frc, rs_frc],[0, rho_max],'r--')
    plt.plot([r0_frc, r0_frc],[0, rho_max],'r--')
xlim( (0., Rmax) )  
#--------------
plt.subplot(224) #-------------------------
plt.hold(True)
plt.grid(True)
rho_mx= np.max(rho_bin_center[:])
binvol_mx= np.max(binvol[:])
plt.xlabel(r'$\rho$ $coordinate$')
text(0.5,binvol_mx,'$ binvol $   $ (cm^3) $')
plt.grid(True)
plt.plot(rho_bin_center[:],binvol[:],'b',linewidth=1.5)  
plt.plot([0, rho_mx],[0, 0],'k')
xlim( (0., rho_mx) )
ylim( (0., binvol_mx*1.1) )  

savefig('genray_wcw_vs_R.png',format='png')
plt.show() 



fig0=plt.figure()  # Characteristic frequencies vs R coord.
f9=freqcy/1e9  #[0] GHz units
plt.subplot(221) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title( '$-f_{pe}$    $--f_{UH}$    $-.-f_{ce}$' )
plt.ylabel('$f/GHz$')
plt.xlabel('$R$   $(cm)$')
plt.yscale('log') 
plt.grid(True)
plt.plot(rscan[:],wpw[:]*f9,'r',linewidth=1.5) # fpe
plt.plot(rscan[:],wuw[:]*f9,'k--',linewidth=1.5) # fuh
plt.plot(rscan[:],wcw[:]*f9,'b-.',linewidth=1.5) # fce  GHz
fce_max=np.amax(wcw[:]*f9)
fpe_max=np.amax(wpw[:]*f9)
fpe_min=np.amin(wpw[:]*f9)
print 'max of fpe[GHz] =',fpe_max, '  (along Z[cm]=',zscan,')'
print 'max of fce[GHz] =',fce_max, '  (along Z[cm]=',zscan,')'
fff_max=max(3*fce_max,3*fpe_min)
xlim( (0., Rmax) )  # limits in R axis in cm
# set the upper limit for this plot:
ylim( (0.1, 100) )    # in GHz
#--------------
plt.subplot(222) #-------------------------
txt1= ' $Scan$  $at$ $(y,z)=($' +"%5.0f" %(yscan) +',' +\
"%5.0f" %(zscan) +'$)cm$' 
txt2=  r"$f (GHz) =$"+r"$%1.2f$" %(f*1e-9)
dx=0.7/3.0
plt.text(0.03, 0.7-dx*0, txt1 , va='center',fontsize=fnt+4)
plt.text(0.03, 0.7-dx*1, txt2 , va='center',fontsize=fnt+4) 
plt.axis([0., 1., 0., 1.])
plt.axis('off')
#--------------
savefig('genray_fce_vs_R.png',format='png')
plt.show() 



if plot_Nperp_prof==1:     # Nperp2 profiles vs x coord.
#--------------------------------------------------------------------------
    for i in range(0,inpar,1):
        fig0=plt.figure()  # Nperp2 profiles vs x coord.
        plt.hold(True)
        plt.title('$ioxm=-1 (blue)$  $and$  $ioxm=+1(green)$     $For$ $N_{||}=$'\
        +"%4.3f" %(Npara[i]))
        plt.ylabel('$N_{\perp} ^2 $     $Cold$ $plasma$ $roots$ $(id=2)$')
        plt.xlabel('$x$ $(cm)$   $Scan$ $at$ $(y,z)=($' +"%5.0f" %(yscan) +',' +\
"%5.0f" %(zscan) +'$)cm$' )
        plt.grid(True)
        plt.plot(rscan[:],Nperp2m[i,:],'b',linewidth=1.5)  
        plt.plot(rscan[:],Nperp2p[i,:],'g',linewidth=1.0)
        ###plt.plot(rscan[:],Nperp2hot[i,:],'r--',linewidth=1.0)
        plt.plot([0, Rmax],[0, 0],'k')
        xlim( (0., Rmax) )
        # Comment next line if you want automatic limits in Nperp (cold)
        ylim((-10.,+10.))  # Un-comment if you want specific limits in Nperp
        savefig('genray_Nperp_vs_R_'+"%3.3f" %(Npara[i])+'.png',format='png')
        print 'file saved for Npara=',Npara[i]
        #plt.show() 
    #stop


if plot_Nperp_prof==1:   # ddd(r,Nper,Npar) profiles vs x,Nperp coord.
#--------------------------------------------------------------------------
    X,Y = np.meshgrid(rscan, Npera)
    #D=ddd[0,:,:]
    #print shape(D)
    #print shape(X)
    for i in range(0,inpar,1):
        fig0=plt.figure()  # ddd(r,Nper,Npar) profiles vs x coord.
        plt.hold(True)
        plt.title('$Dispersion$ $function$ $D(x,N_{\perp})=0$ $levels$  $For$ $N_{||}=$'\
        +"%4.3f" %(Npara[i]))
        plt.ylabel('$ N_{\perp}$     $Hot$ $roots$ $(id=6)$')
        plt.xlabel('$x$ $(cm)$   $Scan$ $at$ $(y,z)=($' +"%5.0f" %(yscan) +',' +\
"%5.0f" %(zscan) +'$)cm$' )
        plt.grid(True)
        #plt.plot([0, Rmax],[0, 0],'k')
        D=ddd[i,:,:]
        #CS=plt.contour(X,Y,D[:],Ncont*5,linewidth=linw,cmap=plt.cm.jet)
        #CB=plt.colorbar(orientation='vertical', shrink=0.5, format='%.2e')
        levels=[0]
        CS=plt.contour(X,Y,D[:],levels,linewidth=linw,colors='r')
        axis([0,Rmax,0.,np.max(Npera)])    
        savefig('genray_ddd_vs_R-Nperp_'+"%3.3f" %(Npara[i])+'.png',format='png')
        print 'file saved for Npara=',Npara[i]
        #plt.show() 
    #stop






#============== MORE DATA, FROM RAY-TRACING ====================================
print '---------------------------------------- RAYS DATA and PLOTS'
# Averaged current densities
s_cur_den_parallel=dat.variables['s_cur_den_parallel']
print 's_cur_den_parallel is: ', s_cur_den_parallel.long_name, s_cur_den_parallel.shape
s_cur_den_onetwo=dat.variables['s_cur_den_onetwo']
print 's_cur_den_onetwo is: ', s_cur_den_onetwo.long_name, s_cur_den_onetwo.shape
s_cur_den_toroidal=dat.variables['s_cur_den_toroidal']
print 's_cur_den_toroidal is: ', s_cur_den_toroidal.long_name, s_cur_den_toroidal.shape
s_cur_den_poloidal=dat.variables['s_cur_den_poloidal']
print 's_cur_den_poloidal is: ', s_cur_den_poloidal.long_name, s_cur_den_poloidal.shape

# Power profiles [erg/sec/cm^3]
powden=dat.variables['powden']
print 'powden is: ', powden.long_name, powden.shape
powden_e=dat.variables['powden_e']
print 'powden_e is: ', powden_e.long_name, powden_e.shape
powden_i=dat.variables['powden_i']
print 'powden_i is: ', powden_i.long_name, powden_i.shape

# Note: powden_cl array did not exist in early genray.nc:
try:
    try:
        powden_cl=dat.variables['powden_cl']
    except:
        print('No data on powden_cl in genray.nc')
        n_powden_cl=0
    else:
        n_powden_cl=1
        powden_cl=dat.variables['powden_cl']
        print 'powden_cl is: ', powden_cl.long_name, powden_cl.shape
finally:
    print '----------------------------------------'
    
    

# Total power [erg/sec]
print '---------------------------------------- [erg/sec]'
power_total=dat.variables['power_total'].getValue()  #getValue() for scalar
print 'power_total is: ', power_total
power_inj_total=dat.variables['power_inj_total'].getValue()  #getValue() for scalar
print 'power_inj_total is: ', power_inj_total
powtot_e=dat.variables['powtot_e'].getValue()  #getValue() for scalar
print 'powtot_e is: ', powtot_e
powtot_i=dat.variables['powtot_i'].getValue()  #getValue() for scalar
print 'powtot_i is: ', powtot_i
powtot_cl=dat.variables['powtot_cl'].getValue()  #getValue() for scalar
print 'powtot_cl is: ', powtot_cl
# For iabsorp.eq.3 only:
#powtot_s=dat.variables['powtot_s'].getValue()  #getValue() for scalar
#print 'powtot_s is: ', powtot_s
print '---------------------------------------- [erg/sec]'

ws=dat.variables['ws']
print 'ws is: ', ws.long_name, ws.shape

delpwr=dat.variables['delpwr']
print 'delpwr is: ', delpwr.long_name, delpwr.shape
#sdpwr=dat.variables['sdpwr']
#print 'sdpwr is: ', sdpwr.long_name, sdpwr.shape
spower=dat.variables['spower']
print 'spower is: ', spower.long_name, spower.shape


# RAY trajectories
wx=dat.variables['wx']
print 'wx is: ', wx.long_name, wx.shape
wy=dat.variables['wy']
print 'wy is: ', wy.long_name, wy.shape
wr=dat.variables['wr']
print 'wr is: ', wr.long_name, wr.shape
wz=dat.variables['wz']
print 'wz is: ', wz.long_name, wz.shape
wphi=dat.variables['wphi']
print 'wphi is: ', wphi.long_name, wphi.shape
# Number of elements for each ray:
nrayelt=dat.variables['nrayelt']
print 'nrayelt is: ', nrayelt.long_name, nrayelt.shape
# Refractive indices along rays:
wnper=dat.variables['wnper']
print 'wnper is: ', wnper.long_name, wnper.shape
wnpar=dat.variables['wnpar']
print 'wnpar is: ', wnpar.long_name, wnpar.shape
wn_x=dat.variables['wn_x']
print 'wn_x is: ', wn_x.long_name, wn_x.shape
wn_y=dat.variables['wn_y']
print 'wn_y is: ', wn_y.long_name, wn_y.shape
wn_r=dat.variables['wn_r']
print 'wn_r is: ', wn_r.long_name, wn_r.shape
wn_z=dat.variables['wn_z']
print 'wn_z is: ', wn_z.long_name, wn_z.shape
wn_phi=dat.variables['wn_phi']
print 'wn_phi is: ', wn_phi.long_name, wn_phi.shape
# E-wave-field along rays
cwexde=dat.variables['cwexde']
print 'cwexde is: ', cwexde.long_name, cwexde.shape
cweyde=dat.variables['cweyde']
print 'cweyde is: ', cweyde.long_name, cweyde.shape
cwezde=dat.variables['cwezde']
print 'cwezde is: ', cwezde.long_name, cwezde.shape
# fluxn = Power Flux along rays:
fluxn=dat.variables['fluxn']
print 'fluxn is: ', fluxn.long_name, fluxn.shape
# Total magnetic field along rays:
sbtot=dat.variables['sbtot']
print 'sbtot is: ',sbtot.long_name, sbtot.shape, sbtot.units
# density along rays [1/cm^3]
sene=dat.variables['sene']
print 'sene is: ', sene.long_name, sene.shape, sene.units
# Ki along rays
salphal=dat.variables['salphal']
print 'salphal is: ', salphal.long_name, salphal.shape

# Vgroup normalized to c:
vgr_x=dat.variables['vgr_x']
print 'vgr_x is: ', vgr_x.long_name, vgr_x.shape
vgr_y=dat.variables['vgr_y']
print 'vgr_y is: ', vgr_y.long_name, vgr_y.shape
vgr_z=dat.variables['vgr_z']
print 'vgr_z is: ', vgr_z.long_name, vgr_z.shape

# Number of rays
Nrays=wz[:,0].size  
print 'Number of rays: Nrays=',Nrays

print '----------------------------------------'



#--------------------------------------------------------------------------
fig4=plt.figure()  # RAYS in top-view (R-phi or X-Y)
#ax = plt.subplot(1, 1, 1)
#ax.set_aspect(1.0)
plt.hold(True)
plt.grid(True)
plt.title('$Rays$  $x(t),y(t)$  $and$  $Electron$ $Density$ $at$ $z=0$  $(cm^{-3})$')
plt.xlabel('$X$  $(cm)$')
plt.ylabel('$Y$  $(cm)$')
X,Y = np.meshgrid(eqdsk_x, eqdsk_y)
X=X*100
Y=Y*100
D=densprofxy
CS=plt.contour(X,Y,D[:],Ncont,linewidth=linw,cmap=cm.jet)
CB=plt.colorbar(orientation='vertical', shrink=0.5, format='%.2e') 
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    print i,Nm
    X=wx[i,0:Nm]
    Y=wy[i,0:Nm]
    plt.plot(X,Y,color=col,linewidth=linw)  
    # plot arrow for refractive vector (Nx,Nz) at starting point: 
    Nx0=wn_x[i,0]
    Ny0=wn_y[i,0]
    N0=sqrt(Nx0*Nx0+Ny0*Ny0)
    #plt.arrow(wx[i,0],wy[i,0],(Nx0/N0)*arr_len,(Ny0/N0)*arr_len,'->') 
plt.plot(wx[0,0],wy[0,0],'ko') # small circle at launching point
plt.plot(bndry_in_X,  bndry_in_Y,  'k',linewidth=linw*2)
plt.plot(bndry_out_X, bndry_out_Y, 'k',linewidth=linw*2)
plt.axis('equal')    
plt.savefig('genray_rays_in_R-phi.png') 
plt.show()

#stop


#--------------------------------------------------------------------------
fig0=plt.figure() # RAYS and Density profile in cross-sectional view X-Z
ax = plt.subplot(111)
ax.set_aspect(1.0)
ax.axis([Rmin,Rmax,Zmin,Zmax])
plt.hold(True)
xmin=amin(eqdsk_x)*100*1.1  # Limits for plots
xmax=amax(eqdsk_x)*100*1.1  # m->cm, and add abit
ymin=amin(eqdsk_y)*100*1.1  # Limits for plots
ymax=amax(eqdsk_y)*100*1.1  # m->cm, and add abit
zmin=amin(eqdsk_z)*100*1.1  # Limits for plots
zmax=amax(eqdsk_z)*100*1.1  # m->cm, and add abit
#ax.axis([xmin,xmax,zmin,zmax])
plt.title('$Rays$  $x(t),z(t)$  $and$  $Electron$ $Density$  $(cm^{-3})$')
plt.xlabel('$X$  $(cm)$')
plt.ylabel('$Z$  $(cm)$')
X,Z = np.meshgrid(eqdsk_x, eqdsk_z)
X=X*100
Z=Z*100
D=densprofxz
CS=plt.contour(X,Z,D[:],Ncont,linewidth=linw,cmap=cm.jet)
CB=plt.colorbar(orientation='vertical', shrink=0.5, format='%.2e') 
#plt.plot(eqdsk_x, densprofxz[32,:])
plt.grid(True)
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    print i,Nm
    plt.plot(wx[i,0:Nm],wz[i,0:Nm],color=col,linewidth=linw)
    # plot arrow for refractive vector (Nx,Nz) at starting point: 
    Nx0=wn_x[i,0]
    Nz0=wn_z[i,0]
    N0=sqrt(Nx0*Nx0+Nz0*Nz0)
    #plt.arrow(wx[i,0],wz[i,0],(Nx0/N0)*arr_len,(Nz0/N0)*arr_len,'->') 
# plot small circle at the launching point:
plt.plot(wx[0,0],wz[0,0],'ko')
# plot walls, if any:
if n_wall>>0:
    plt.plot(r_wall,z_wall,'k',linewidth=linw*2)
plt.savefig('genray_rays_Ne_inXZ.png') 
plt.show()

#stop

#--------------------------------------------------------------------------
fig0=plt.figure() # RAYS and Density profile in cross-sectional view Y-Z
ax = plt.subplot(111)
ax.set_aspect(1.0)
ax.axis([Rmin,Rmax,Zmin,Zmax])
#ax.axis([ymin,ymax,zmin,zmax])
plt.hold(True)
plt.title('$Rays$  $y(t),z(t)$  $and$  $Electron$ $Density$  $(cm^{-3})$')
plt.xlabel('$Y$  $(cm)$')
plt.ylabel('$Z$  $(cm)$')
Y,Z = np.meshgrid(eqdsk_y, eqdsk_z)
Y=Y*100
Z=Z*100
D=densprofyz
CS=plt.contour(Y,Z,D[:],Ncont,linewidth=linw,cmap=cm.jet)
CB=plt.colorbar(orientation='vertical', shrink=0.5, format='%.2e') 
#plt.plot(eqdsk_y, densprofyz[32,:])
plt.grid(True)
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    plt.plot(wy[i,0:Nm],wz[i,0:Nm],color=col,linewidth=linw)
    # plot arrow for refractive vector (Ny,Nz) at starting point: 
    Ny0=wn_y[i,0]
    Nz0=wn_z[i,0]
    N0=sqrt(Ny0*Ny0+Nz0*Nz0)
    #plt.arrow(wy[i,0],wz[i,0],(Ny0/N0)*arr_len,(Nz0/N0)*arr_len,'->') 
# plot small circle at the launching point:
plt.plot(wy[0,0],wz[0,0],'ko')
# plot walls, if any:
if n_wall>>0:
    plt.plot(r_wall,z_wall,'k',linewidth=linw*2)
plt.savefig('genray_rays_Ne_inYZ.png') 
plt.show()



#--------------------------------------------------------------------------
fig3=plt.figure()  # RAYS and Pol.Flux in cross-section view X-Z
ax = plt.subplot(1, 1, 1)
ax.set_aspect(1.0)
ax.axis([Rmin,Rmax,Zmin,Zmax])
plt.hold(True)
plt.title('$Rays$  $x(t),z(t)$  $and$  $\Psi_p/2\pi$  $(Tesla*m^2)$')
plt.xlabel('$X$  $(cm)$')
plt.ylabel('$Z$  $(cm)$')
R,Z = np.meshgrid(eqdsk_r, eqdsk_z)
R=R*100
Z=Z*100
PSI = eqdsk_psi
PSImin=np.min(PSI)
PSImax=np.max(PSI)
if PSImin<0.:
    PSImax=-PSImin # Note: for FRC: PSI=0 at r=0 and r=rs
levels=np.arange(PSImin,PSImax,(PSImax-PSImin)/(Ncont-1))
CS=plt.contour( R,Z,PSI[:],levels,linewidth=linw,cmap=cm.jet)
CS=plt.contour(-R,Z,PSI[:],levels,linewidth=linw,cmap=cm.jet)
CB=plt.colorbar(orientation='vertical', shrink=0.5, format='%.2e') 
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    plt.plot(wx[i,0:Nm],wz[i,0:Nm],color=col,linewidth=linw)
    # plot arrow for refractive vector (Nx,Nz) at starting point: 
    Nx0=wn_x[i,0]
    Nz0=wn_z[i,0]
    N0=sqrt(Nx0*Nx0+Nz0*Nz0)
    #plt.arrow(wx[i,0],wz[i,0],(Nx0/N0)*arr_len,(Nz0/N0)*arr_len,'->',color=col) 
# plot small circle at the launching point:
plt.plot(wx[0,0],wz[0,0],'ko')
# plot walls, if any:
if n_wall>>0:
    plt.plot(r_wall,z_wall,'k',linewidth=linw*2)
plt.savefig('genray_rays_PolFlux_inXZ.png') 
plt.show()

#stop

#--------------------------------------------------------------------------
fig8=plt.figure()  # E-wave-field along RAYS vs pol.distance(t)
plt.subplot(231) #-------------------------
plt.hold(True)
plt.title('$|E_X/E|$ $along$ $ray$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    ereal=cwexde[0,i,0:Nm].copy()
    eimag=cwexde[1,i,0:Nm].copy()
    ea= sqrt(ereal**2 + eimag**2)
    plt.plot(ws[i,0:Nm],ea,color=col,linewidth=linw)  
plt.subplot(232) #-------------------------
plt.hold(True)
plt.title('$|E_Y/E|$ $along$ $ray$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    ereal=copy(cweyde[0,i,0:Nm])
    eimag=cweyde[1,i,0:Nm].copy()
    ea= sqrt(ereal**2 + eimag**2)
    plt.plot(ws[i,0:Nm],ea,color=col,linewidth=linw)  
plt.subplot(233) #-------------------------
plt.hold(True)
plt.title('$|E_Z/E|$ $along$ $ray$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    ereal=cwezde[0,i,0:Nm].copy()
    eimag=cwezde[1,i,0:Nm].copy()
    ea= sqrt(ereal**2 + eimag**2)
    plt.plot(ws[i,0:Nm],ea,color=col,linewidth=linw)
plt.subplot(234) #-------------------------
plt.hold(True)
plt.title('$Power Flux$ $along$ $ray$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    plt.plot(ws[i,0:Nm],fluxn[i,0:Nm],color=col,linewidth=linw)
plt.xlabel('$distance$ $along$ $ray$  $(cm)$')
plt.subplot(236) #-------------------------
plt.hold(True)
plt.title('  $Lin.Damp.$'+' $K_i$'+' $(cm^{-1})$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    plt.plot(ws[i,0:Nm],salphal[i,0:Nm],color=col,linewidth=linw)
plt.xlabel('$distance$ $along$ $ray$  $(cm)$')
plt.savefig('genray_rays_Ewave_p.png') 
plt.show()


#--------------------------------------------------------------------------
fig81=plt.figure()  # Vgroup/c along RAYS vs pol.distance(t)
plt.subplot(231) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title('$V_{group,x}/c$  $along$ $ray$')
Vgr_mn= np.amin(vgr_x[:])
Vgr_mx= np.amax(vgr_x[:])
Vgr_mn=max(Vgr_mn,-1.0)
Vgr_mx=min(Vgr_mx, 1.0)
ylim((Vgr_mn-1e-4,Vgr_mx+1e-4))
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    plt.plot(ws[i,0:Nm],vgr_x[i,0:Nm],color=col,linewidth=linw)  
plt.subplot(234) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title('$V_{group,y}/c$  $along$ $ray$')
Vgr_mn= np.amin(vgr_y[:])
Vgr_mx= np.amax(vgr_y[:])
Vgr_mn=max(Vgr_mn,-1.0)
Vgr_mx=min(Vgr_mx, 1.0)
ylim((Vgr_mn-1e-4,Vgr_mx+1e-4))
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    plt.plot(ws[i,0:Nm],vgr_y[i,0:Nm],color=col,linewidth=linw)  
plt.xlabel('$distance$ $along$ $ray$  $(cm)$')
plt.subplot(233) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title('$V_{group,z}/c$  $along$ $ray$')
Vgr_mn= np.amin(vgr_z[:])
Vgr_mx= np.amax(vgr_z[:])
Vgr_mn=max(Vgr_mn,-1.0)
Vgr_mx=min(Vgr_mx, 1.0)
ylim((Vgr_mn-1e-4,Vgr_mx+1e-4))
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    plt.plot(ws[i,0:Nm],vgr_z[i,0:Nm],color=col,linewidth=linw)
plt.subplot(236) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title('$|V_{group}/c|$')
vgr=sqrt(vgr_x[:]**2 +vgr_y[:]**2 +vgr_z[:]**2)
#print 'vgr shape:', vgr.shape
Vgr_mn= np.amin(vgr[:])
Vgr_mx= np.amax(vgr[:])
Vgr_mn=max(Vgr_mn, 0.0)
Vgr_mx=min(Vgr_mx, 1.0)
ylim((Vgr_mn,Vgr_mx))
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    plt.plot(ws[i,0:Nm],vgr[i,0:Nm],color=col,linewidth=linw)
plt.xlabel('$distance$ $along$ $ray$  $(cm)$')
plt.savefig('genray_rays_Vgroup_p.png') 
plt.show()

#stop

#--------------------------------------------------------------------------
fig5=plt.figure()  # wce/w and (wpe/w)^2 along RAYS vs distance
f=freqcy
q2m= charge[0]**2/mass[0]
plt.subplot(231) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title('$\omega_{ce}/\omega$  $along$ $rays$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    plt.plot(ws[i,0:Nm],28e5*sbtot[i,0:Nm]/f,color=col,linewidth=linw)
plt.xlabel('$distance$ $along$ $ray$  $(cm)$')

plt.subplot(232) #-------------------------
txt=r"$f (GHz) =$"+r"$%1.2f$" %(f*1e-9)
plt.title(txt)   # 
plt.axis('off')

plt.subplot(233) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title('$(\omega_{pe}/\omega)^2$  $along$ $rays$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    #print 806.2e5*sene[0,0]*q2m/f**2, sene[0,0], q2m
    plt.plot(ws[i,0:Nm],806.2e5*sene[i,0:Nm]*q2m/f**2,color=col,linewidth=linw)
plt.xlabel('$distance$ $along$ $ray$  $(cm)$')
#--------------------------------------------
plt.subplot(234) #-------------------------
plt.hold(True)
plt.grid(True)
plt.ylabel('$\omega/\omega_{ce}$  $along$ $rays$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    plt.plot(ws[i,0:Nm],f/(28e5*sbtot[i,0:Nm]),color=col,linewidth=linw)
plt.xlabel('$distance$ $along$ $ray$  $(cm)$')

plt.savefig('genray_rays_wc_wp_p.png') 
plt.show()

#stop

#--------------------------------------------------------------------------
fig5=plt.figure()  # Refractive indices along RAYS vs R(t)
ax=plt.subplot(231) #-------------------------
plt.hold(True)
plt.title('$N_R$ $along$ $ray$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    plt.plot(wr[i,0:Nm],wn_r[i,0:Nm],color=col,linewidth=linw)  
plt.subplot(232) #-------------------------
plt.hold(True)
plt.title('$N_Z$ $along$ $ray$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    plt.plot(wr[i,0:Nm],wn_z[i,0:Nm],color=col,linewidth=linw)  
plt.subplot(233) #-------------------------
plt.hold(True)
plt.title(r'      $N_\phi$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    plt.plot(wr[i,0:Nm],wn_phi[i,0:Nm],color=col,linewidth=linw)
plt.subplot(234) #-------------------------
plt.hold(True)
plt.title(r'$N_{\perp}$'+' $along$ $ray$')
plt.xlabel('$R$ $along$ $ray$  $(cm)$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    plt.plot(wr[i,0:Nm],wnper[i,0:Nm],color=col,linewidth=linw)
plt.subplot(235) #-------------------------
plt.hold(True)
plt.title(r'$N_{||}$'+' $along$ $ray$')
plt.xlabel('$R$ $along$ $ray$  $(cm)$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    plt.plot(wr[i,0:Nm],wnpar[i,0:Nm],color=col,linewidth=linw)
plt.subplot(236) #-------------------------
plt.hold(True)
plt.title('$   Lin.Damping$'+' $K_i$'+' $(cm^{-1})$')
plt.xlabel('$R$ $along$ $ray$  $(cm)$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    plt.plot(wr[i,0:Nm],salphal[i,0:Nm],color=col,linewidth=linw)
plt.savefig('genray_rays_refr-index_R.png') 
plt.show()

#stop




#--------------------------------------------------------------------------
fig6=plt.figure()  # Refractive indices along RAYS vs pol.distance(t)
plt.subplot(231) #-------------------------
plt.hold(True)
plt.title('$N_R$ $along$ $ray$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    plt.plot(ws[i,0:Nm],wn_r[i,0:Nm],color=col,linewidth=linw)  
plt.subplot(232) #-------------------------
plt.hold(True)
plt.title('$N_Z$ $along$ $ray$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    plt.plot(ws[i,0:Nm],wn_z[i,0:Nm],color=col,linewidth=linw)  
plt.subplot(233) #-------------------------
plt.hold(True)
plt.title(r'$N_\phi$'+' $along$ $ray$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    plt.plot(ws[i,0:Nm],wn_phi[i,0:Nm],color=col,linewidth=linw)
plt.subplot(234) #-------------------------
plt.hold(True)
plt.title(r'$N_{\perp}$'+' $along$ $ray$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    plt.plot(ws[i,0:Nm],wnper[i,0:Nm],color=col,linewidth=linw)
plt.subplot(235) #-------------------------
plt.hold(True)
plt.title(r'$N_{||}$'+' $along$ $ray$')
plt.xlabel('$distance$ $along$ $ray$  $(cm)$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    plt.plot(ws[i,0:Nm],wnpar[i,0:Nm],color=col,linewidth=linw)
plt.subplot(236) #-------------------------
plt.hold(True)
plt.title(r'$|N|$'+' $along$ $ray$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    plt.plot(ws[i,0:Nm],sqrt(wnpar[i,0:Nm]**2+wnper[i,0:Nm]**2),color=col,linewidth=linw)
plt.savefig('genray_rays_refr-index_p.png') 
plt.show()


#--------------------------------------------------------------------------
fig2=plt.figure()
plt.hold(True)
plt.title(delpwr.long_name)
plt.xlabel('$distance$ $along$ $ray$' + '  $ws$  $(cm)$')
plt.ylabel('$ delpwr $  $ (ergs/sec) $')
ws_max=0.
plt.grid(True)
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    ws_max1=np.amax(ws[:,Nm-1])
    ws_max=np.amax([ws_max1,ws_max])
    plt.plot(ws[i,0:Nm],delpwr[i,0:Nm],color=col,linewidth=6*linw-0.6*i)
axis([0.,1.05*ws_max,0.,1.05*np.amax(delpwr[:,0])])
plt.savefig('genray_rays_delpwr.png') 
plt.show()



#--------------------------------------------------------------------------
if i_ox==2:  # only for O-X transmission case
    fig8=plt.figure()
    plt.subplot(221) #-------------------------
    plt.hold(True)
    plt.grid(True)
    #plt.xlabel('$N_{||}$ $at$ $starting$ $point$')
    plt.ylabel('$1=converted;$   $0-not$')
    plt.title('Indicator of OX conversion')
    ylim((0.,1.1))
    for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
        if remainder(i,6)==0: col='b'
        if remainder(i,6)==1: col='g'
        if remainder(i,6)==2: col='r'
        if remainder(i,6)==3: col='c'    
        if remainder(i,6)==4: col='m' 
        if remainder(i,6)==5: col='k'  
        plt.plot(wnpar[i,0],i_ox_conversion[i],'o',color=col)
    plt.subplot(222) #-------------------------
    plt.hold(True)
    plt.grid(True)
    plt.xlabel('$N_{||}$ $at$ $starting$ $point$')
    #text(wnpar[0,0],1.02,'   $pwr(Xmode)/pwr(Omode)$')
    plt.title('Transm.coef=P(X)/P(O)')
    
    ylim((0.,1.1))
    for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
        if remainder(i,6)==0: col='b'
        if remainder(i,6)==1: col='g'
        if remainder(i,6)==2: col='r'
        if remainder(i,6)==3: col='c'    
        if remainder(i,6)==4: col='m' 
        if remainder(i,6)==5: col='k'  
        plt.plot(wnpar[i,0],transm_ox[i],'o',color=col)
    plt.subplot(223) #-------------------------
    plt.hold(True)
    plt.grid(True)
    plt.xlabel('$N_{||}$ $at$ $starting$ $point$')
    plt.ylabel('$N_{||}$ $at$ $OX$ $transm.$')
    y_mn= min(0.,1.1*np.amin(cnpar_ox[:]))
    y_mx= max(0.,1.1*np.amax(cnpar_ox[:]))
    ylim((y_mn,y_mx))
    for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
        if remainder(i,6)==0: col='b'
        if remainder(i,6)==1: col='g'
        if remainder(i,6)==2: col='r'
        if remainder(i,6)==3: col='c'    
        if remainder(i,6)==4: col='m' 
        if remainder(i,6)==5: col='k'  
        plt.plot(wnpar[i,0],cnpar_ox[i],'o',color=col)
#    plt.subplot(224) #-------------------------
#    plt.hold(True)
#    plt.xlabel('$N_{||}$ $at$ $starting$ $point$')
#    plt.ylabel('$N$ $along$ $B$x$grad(\Psi)$')
#    for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
#        if remainder(i,6)==0: col='b'
#        if remainder(i,6)==1: col='g'
#        if remainder(i,6)==2: col='r'
#        if remainder(i,6)==3: col='c'    
#        if remainder(i,6)==4: col='m' 
#        if remainder(i,6)==5: col='k'  
#        plt.plot(wnpar[i,0],cn_b_gradpsi[i],'o',color=col)
    
    plt.savefig('genray_rays_OX_transmission.png') 
    plt.show()



#--------------------------------------------------------------------------
fig7=plt.figure()   # Power profiles

powden_e= powden_e[:]/1.0e10  # kW/cm^3 (unless binvol=1.0, then kW)
powden_i= powden_i[:]/1.0e10  # kW/cm^3 (unless binvol=1.0, then kW)
powden_cl=powden_cl[:]/1.0e10 # kW/cm^3 (unless binvol=1.0, then kW)
powden=   powden[:]/1.0e10    # kW/cm^3 (unless binvol=1.0, then kW)
#------------
powtot_e=    powtot_e/1.0e10  # kW
powtot_i=    powtot_i/1.0e10  # kW
powtot_cl=   powtot_cl/1.0e10 # kW
power_total= power_total/1.0e10 # kW
power_inj_total= power_inj_total/1.0e10 # kW

NR=rho_bin_center[:].size

if rhomax==0:
    rhomax=1.05*np.amax(rho_bin_center)
rhomax=min(rhomax, rho_bin_center[NR-1])
print 'rho_bin_center: min/max',np.amin(rho_bin_center),np.amax(rho_bin_center)
print 'rho_bin: min/max', np.amin(rho_bin),np.amax(rho_bin)
drho= rho_bin_center[1]-rho_bin_center[0]  # step in rho grid
irmx= int(rhomax/drho) # index corresponding to rhomax
irmx=min(irmx,NR-1)

plt.subplot(221) #-------------------------
txt="$p_{e}$"+"$=$%3.3f" %(powtot_e) +" $kW$"
plt.title(txt)
plt.hold(True)
plt.grid(True)
axis([0.,rhomax,0.,1.05*np.amax(powden_e[0:irmx])+0.01])
plt.plot(rho_bin_center[0:irmx],powden_e[0:irmx],'r.',linewidth=linw)
plt.plot(rho_bin_center[0:irmx],powden_e[0:irmx],'k',linewidth=linw)
plt.subplot(222) #-------------------------
plt.hold(True)
plt.grid(True)
txt="$p_{i}$"+"$=$%3.3f" %(powtot_i) +" $kW$"
plt.title(txt)
axis([0.,rhomax,0.,1.05*np.amax(powden_i[0:irmx])+0.01])
plt.plot(rho_bin_center[0:irmx],powden_i[0:irmx],'r.',linewidth=linw)
plt.plot(rho_bin_center[0:irmx],powden_i[0:irmx],'k',linewidth=linw)
plt.subplot(223) #-------------------------
plt.hold(True)
plt.grid(True)
plt.xlabel(r'$\rho$')
txt="$p_{cl}$"+"$=$%3.3f" %(powtot_cl) +" $kW$"
plt.title(txt,verticalalignment='center')
if n_powden_cl>0:
    axis([0.,rhomax,0.,1.05*np.amax(powden_cl[0:irmx])+0.01])
    plt.plot(rho_bin_center[0:irmx],powden_cl[0:irmx],'r.',linewidth=linw)
    plt.plot(rho_bin_center[0:irmx],powden_cl[0:irmx],'k',linewidth=linw)
plt.subplot(224) #-------------------------
plt.hold(True)
plt.grid(True)
plt.xlabel(r'$\rho$')
txt="$p_{total}$"+"$=$%3.3f" %(power_total) +" $kW$"
plt.title(txt,verticalalignment='center')
axis([0.,rhomax,0.,1.05*np.amax(powden[0:irmx])+0.01])
plt.plot(rho_bin_center[0:irmx],powden[0:irmx],'r.',linewidth=linw)
plt.plot(rho_bin_center[0:irmx],powden[0:irmx],'k',linewidth=linw)
savefig('genray_profiles_power.png')
plt.show() 


print 'Sum over rays:'
print 'integral(powden*binvol)=    ', sum(powden*binvol),   ' kW'
print 'integral(powden_e*binvol)=  ', sum(powden_e*binvol), ' kW'
print 'integral(powden_i*binvol)=  ', sum(powden_i*binvol), ' kW'
print 'integral(powden_cl*binvol)= ', sum(powden_cl*binvol),' kW'
print 'Note: check plots of binvol: Could be set to 1.0 if no flux surfaces'
print 'ABSORBED/LOST-in-plasma:  power_total =   ',power_total,    ' kW'
print 'INJECTED (STARTED power): power_inj_total=',power_inj_total,' kW'
print '==========================================================='
#--------------------------------------------------------------------------
#fig1=plt.figure()   # Current profiles
#plt.subplot(221) #-------------------------
#title('Current Profiles  '+'$(A/cm^2)$')
#plt.hold(True)
#plt.grid(True)
#plt.ylabel('$<j_{||}>$')
#plt.plot(rho_bin_center,s_cur_den_parallel,'r',linewidth=linw)
#plt.subplot(222) #-------------------------
#plt.hold(True)
#plt.grid(True)
#plt.title('<j.B>'+'$/B_0$')
#plt.plot(rho_bin_center,s_cur_den_onetwo,'r',linewidth=linw)
#plt.subplot(223) #-------------------------
#plt.hold(True)
#plt.grid(True)
#plt.xlabel(r'$\rho$')
#plt.title('$<j_{||}>f<1/r^2>/(<B><1/r>)$',verticalalignment='center')
#plt.plot(rho_bin_center,s_cur_den_toroidal,'r',linewidth=linw)
#plt.subplot(224) #-------------------------
#plt.hold(True)
#plt.grid(True)
#plt.xlabel(r'$\rho$')
#plt.title('$<j_{||}>B_{pol}($'+r'$\theta$'+'$=0)/<B>$',verticalalignment='center')
#plt.plot(rho_bin_center,s_cur_den_poloidal,'r',linewidth=linw)
#savefig('genray_profiles_J.png')
#plt.show() 



#--------------------------------------------------------------------------
fig10=plt.figure()   # Power profiles over R,Z grid

Rgr,Zgr = np.meshgrid(Rgrid,Zgrid) # 2D grids [cm]
Nc_pwr=10 # number of contour levels in power(R,Z) plots

ax=plt.subplot(131) #-------------------------
ax.set_aspect(1.0)
ax.axis([0,Rmax,Zmin,Zmax])
txt="$p_{e}$"+"$=$%3.3f" %(sum_spwr_rz_e) +" $kW$"
plt.title(txt)
plt.xlabel('$R$  $(cm)$')
plt.ylabel('$Z$  $(cm)$')
if sum_spwr_rz_e>0:
    CSP=plt.contour(Rgr,Zgr,spwr_rz_e[:],Nc_pwr,cmap=plt.cm.jet)
#CBP=plt.colorbar(orientation='vertical', shrink=0.4,format='%.1e') # colorbar
levels=np.arange(PSImin,PSImax,(PSImax-PSImin)/20)
CS=plt.contour( R,Z,PSI[:],levels,linewidth=linw,cmap=cm.jet)
# plot small circle at the launching point:
plt.plot(wr[0,0],wz[0,0],'ko')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    # Rays projected to (R,Z) plane:
    plt.plot(wr[i,0:Nm],wz[i,0:Nm],color=col,linewidth=linw)
    
ax=plt.subplot(132) #-------------------------
ax.set_aspect(1.0)
ax.axis([0,Rmax,Zmin,Zmax])
txt="$p_{i}$"+"$=$%3.3f" %(sum_spwr_rz_i) +" $kW$"
plt.title(txt)
plt.xlabel('$R$  $(cm)$')
#plt.ylabel('$Z$  $(cm)$')
if sum_spwr_rz_i>0:
    CSP=plt.contour(Rgr,Zgr,spwr_rz_i[:],Nc_pwr,cmap=plt.cm.jet)
#CBP=plt.colorbar(orientation='vertical', shrink=0.4,format='%.1e') # colorbar
levels=np.arange(PSImin,PSImax,(PSImax-PSImin)/20)
CS=plt.contour( R,Z,PSI[:],levels,linewidth=linw,cmap=cm.jet)
# plot small circle at the launching point:
plt.plot(wr[0,0],wz[0,0],'ko')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    # Rays projected to (R,Z) plane:
    plt.plot(wr[i,0:Nm],wz[i,0:Nm],color=col,linewidth=linw)

ax=plt.subplot(133) #-------------------------
ax.set_aspect(1.0)
ax.axis([0,Rmax,Zmin,Zmax])
txt="$p_{cl}$"+"$=$%3.3f" %(sum_spwr_rz_cl) +" $kW$"
plt.title(txt)
plt.xlabel('$R$  $(cm)$')
#plt.ylabel('$Z$  $(cm)$')
if sum_spwr_rz_cl>0:
    CSP=plt.contour(Rgr,Zgr,spwr_rz_cl[:],Nc_pwr,cmap=plt.cm.jet)
#CBP=plt.colorbar(orientation='vertical', shrink=0.4,format='%.1e') # colorbar
levels=np.arange(PSImin,PSImax,(PSImax-PSImin)/20)
CS=plt.contour( R,Z,PSI[:],levels,linewidth=linw,cmap=cm.jet)
# plot small circle at the launching point:
plt.plot(wr[0,0],wz[0,0],'ko')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    # Rays projected to (R,Z) plane:
    plt.plot(wr[i,0:Nm],wz[i,0:Nm],color=col,linewidth=linw)  
    
savefig('genray_profiles_power_RZ.png')
plt.show() 



dat.close() # close genray.nc
         
elapsed_time = time.time() - e0
cpu_time = time.clock() - c0
print 'elapsed and cpu time since start (sec.) =', elapsed_time, cpu_time
print 'FINISHED'
