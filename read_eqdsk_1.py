#!/usr/bin/ipython
# If starting in interactive mode,  %autoindent toggles autoindent in ipython.
# To continue the interactive session, do an
# from read_eqdsk import *
# This will leave all variables available at the prompt.


#BH, 110618    
#YuP 110620   Added plots of B-field

from numpy import *
from mpl_toolkits.mplot3d import Axes3D

from pylab import *
from matplotlib import rc 
from matplotlib.pyplot import cm,figure,axes,plot,xlabel,ylabel,title,savefig,show

import os
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pylab as pylab
#matplotlib.interactive(True) # no plots on screen
matplotlib.interactive(False) # with plots on screen



#eqdskin='eqdsk_NSTX' #'eqdsk'  #'g128739.00295'
eqdskin='eqdsk_C2W.out'  #Renamed eqdsk.out (obtained from TAE on 12/30/2016)
eqdsktype='TAE' # Type of eqdsk data file.
#          For TAE (Tri Alpha Energy) FRC case, the eqdsk file 
#          does not contain qpsi array, but instead it contains 
#          necut and tmcut arrays (electron density [1/m^3] and
#          el.temperature [eV]). So, one extra read block is implemented.
#          Options so far: "TAE" (default) or "tokamak" (or "mirror")
# [2017-01-01] New option: 'TAEi' - additional data on ions ???


#set fonts and line thicknesses  ------------------------------
fnt  = 12 #14     # font size for axis numbers (see 'params=' below) 
linw = 1.0    # LineWidth for plots
Ncont= 50     # Number of contour levels
Rmn_plots=0.
Rmx_plots=0. #40. #160. #. # limits for plots.  Set to 0.0 for auto-setup
Zmn_plots=-150 #-170. #-130. # limits for plots.  Set to 0.0 for auto-setup
Zmx_plots=+150 #+170. #130. # limits for plots.  Set to 0.0 for auto-setup
Bmn_plots=0.0 # [Tesla] Lower Limit for contour levels
Bmx_plots=0. #2.5 #0.06 #9e-3 #0.06 # [Tesla] Upper Limit for contour plots
PSImn_plots=0.0 # [Web/rad]  Lower Limit for contour levels
PSImx_plots=0.03 #1.e-3 # [Web/rad]  Upper Limit for contour plots
#--------------------------------------------------------------


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




#==================================================================
def read_vector(flnm,ndim,nums_per_line):
    global nlines  #Can be used to return number of lines.
    """
    Reads delimited items from an open file, flnm.
    Returns them in a vector of input length ndim.
    nums_per_line is (assumed constant) number of
    items per line (except last line may be shorter).
    BH2009.
    """
    a=np.ones(ndim)
    nlines=ndim/nums_per_line
    if nlines*nums_per_line != ndim: nlines=nlines+1
    #print nlines
    for i in range(nlines):
        ibegin=i*nums_per_line
        iend=min(ibegin+nums_per_line,ndim)
        #print ibegin,iend
        a[ibegin:iend]=np.array(flnm.readline().split(),float)
    return a
#==================================================================



#******************************************************************
#******************************************************************

eqdsk=open(eqdskin,'r')
lines=eqdsk.readlines()
#fortran format 110  format(6a8,4i4), picking off 2nd and 3rd slots of 4i4.
nr=lines[0][52:56]; nz=lines[0][56:60]
nr=int(nr); nz=int(nz)
try:
    nv=int(lines[0][60:64])
except:   #line[60:64] not characters representing integers
    nv=nr

print 'nr=',nr, '   nz=',nz, '    nv=',nv
#A difficulty making reading of eqdsk files difficult is that
#when a number is procedded by a "-" sign, then there is no
#blank space separating adjacent numbers.
#For remaining lines after 1st down to possible occurance of
#an "&", fortran format is format(5e16.9) [Except possible 
#sequence of integers extending less that 1st 16 columns.
#Therefore, Add blanks in columns 17,33,49,65, and use
#above read_vector() function.

#Determine line number of "&", if exists:
#(lines[iamp] will be line before "&".)
iamp=0
for line in lines[1:]:
    if line.strip().find('&') == -1:
        iamp=iamp+1
    else:
        break

#Introduce spaces in lines[1:iamp], skipping 1st and & lines.
#This should be ok for standard eqdsks through reading xlimiter,ylimiter.
#BH noticed a nonstandard format (i5,e16.9) further on.
cols=range(16,65,16); cols.reverse()  #cols=[64, 48, 32, 16
b=' '
for i in range(1,iamp+1):
    for j in range(4):
        lines[i]=lines[i][0:cols[j]]+b+lines[i][cols[j]:]

# Open tmp file, put adjusted lines in it, 
# then read using previously constructed read_vector.
tmp=open('./tmp','w')
tmp.writelines(lines)
tmp.close()
tmp=open('./tmp','r')
tmp.readline()  #space down one line
#---
rbox,zbox,radmaj,rboxdst,ymideqd=read_vector(tmp,5,5)
raxis,zaxis,psimag,psilim,btor=read_vector(tmp,5,5)
toteqd,psimx1,psimx2,xax1,xax2=read_vector(tmp,5,5)
zax1,zax2,psisep,xsep,zsep=read_vector(tmp,5,5)
fpsiar=read_vector(tmp,nv,5)
print 'min/max of R*B_phi array fpsiar=', min(fpsiar), max(fpsiar)
prar=read_vector(tmp,nv,5)
print 'min/max of pressure array prar=', min(prar), max(prar)
ffpar=read_vector(tmp,nv,5)
print 'min/max of ffpar=', min(ffpar), max(ffpar)
ppar=read_vector(tmp,nv,5)
print 'min/max of ppar=', min(ppar), max(ppar)
epsi=read_vector(tmp,nr*nz,5)
print 'min/max of epsi(i,j)=', min(epsi), max(epsi), epsi.shape
if eqdsktype=='TAE':
    necut=read_vector(tmp,nr*nz,5)
    print 'min/max of necut=', min(necut), max(necut)
    tmcut=read_vector(tmp,nr*nz,5)
    print 'min/max of tmcut=', min(tmcut), max(tmcut), tmcut.shape
elif eqdsktype=='TAEi':
    # electrons:
    necut=read_vector(tmp,nr*nz,5)
    print 'min/max of necut=', min(necut), max(necut)
    tmcut=read_vector(tmp,nr*nz,5)
    print 'min/max of tmcut=', min(tmcut), max(tmcut), tmcut.shape  
    # ions:    
    ni_thermal=read_vector(tmp,nr*nz,5)
    print 'min/max of ni_thermal=', min(ni_thermal), max(ni_thermal)
    ti_thermal=read_vector(tmp,nr*nz,5)
    print 'min/max of ti_thermal=', min(ti_thermal), max(ti_thermal), ti_thermal.shape        
    ## ions/fast:    
    #ni_fast=read_vector(tmp,nr*nz,5)
    #print 'min/max of ni_fast=', min(ni_fast), max(ni_fast)
    #ti_fast=read_vector(tmp,nr*nz,5)
    #print 'min/max of ti_fast=', min(ti_fast), max(ti_fast), ti_fast.shape        
else:   # tokamak type of eqdsk
    qar=read_vector(tmp,nv,5)
    print 'min/max of qar=', min(qar), max(qar)
    
ncontr, nlimiter=read_vector(tmp,2,2)
ncontr=int(ncontr); nlimiter=int(nlimiter)

#---
# raxis,zaxis      are the major and vertical height of magnetic axis.
# psimag,psilim    are the poloidal flux function values at the
#                  magnetic axis and the last closed flux surface
#                   (touching the limiter or the separatrix).
# btor             is the vacuum toroidal magnetic field at radmaj.
# toteqd           is the toroidal current.
#
# feqd(nnr),pres(nnr), usually, although dimension is
#                      specified by nnv if it is present.
# fpsiar = r * B_phi  at nnr equispaced points in psi
#          from psimag to psilim.
# prar   = are the pressure values at the same points.
# ffpar  = fpsiar_prime, [i.e., d f/d (psi)] at the same points.
# ppeqd  =  p_prime [i.e., d prar/d (psi)] at the same points.
# epsi(nnr,nnz) are the psi values on the nnr * nnz
#               equispaced grid.
# qar(nnr)  gives safety factor q on the equispaced psi grid.
#
# The following quantities are given in some eqdsks, but are not
# presently used in cql3d:
# nlimit,nves  are the numbers of point at limiters and
#              vacuum vessel wall.
# rlimit(nlimit),zlimit(nlimit):
# rlimit,zlimit      is the r,z location of the limiters wall.
#       rves(nves),zves(nves):
# rves,zves          is the r,z location of the limiters
#                    vacuum vessel wall.

if ncontr > 0:
    rzcontr=read_vector(tmp,2*ncontr,5)
    rzcontr.resize((ncontr,2))
    rcontr=rzcontr[:,0]
    zcontr=rzcontr[:,1]

if nlimiter > 0:
    rzlimiter=read_vector(tmp,2*nlimiter,5)
    rzlimiter.resize((nlimiter,2))
    rlimiter=rzlimiter[:,0]
    zlimiter=rzlimiter[:,1]
    
eqdsk.close()
tmp.close()

epsi.resize((nz,nr)) 
if eqdsktype=='TAE':
    tmcut.resize((nz,nr)) 
    necut.resize((nz,nr)) 
    
if eqdsktype=='TAEi':  # electrons and ions (thermal-ion and fast-ion species)
    # electrons:
    tmcut.resize((nz,nr)) 
    necut.resize((nz,nr)) 
    # ions:
    ti_thermal.resize((nz,nr)) 
    ni_thermal.resize((nz,nr)) 
    #ti_fast.resize((nz,nr)) 
    #ni_fast.resize((nz,nr)) 
    
    
print '========================================'
print 'nr=',nr
print 'nz=',nz
print 'nv=',nv
print 'psimag=', psimag
print 'psilim=', psilim
print 'zaxis [m]=', zaxis
print 'raxis [m]=', raxis
print 'radmaj [m]=', radmaj
print 'btor at radmaj [Tesla]=',btor
print 'fpsiar shape =', fpsiar.shape
print 'epsi shape =', epsi.shape
print 'read_eqdsk done'
print '========================================'

#This is all that is read in cql3d: equilib.f.  Additional data
#has been added to the eqdsk over the years, some of which may
#not be standard from one implementation to the next.


# Formats used in reading eqdsk, in cql3d:equilib.f
# 110  format(6a8,4i4)
# 120  format(5e16.9)
# 8200    format(5e16.9)
# 8210    format (2i5)

# convert to cgs:
rbox=rbox*1.e+2 
zbox=zbox*1.e+2
rboxdst=rboxdst*1.e+2
radmaj=radmaj*1.e+2
toteqd=toteqd*3.e+9       
R_axis  = raxis*1.e+2 # cgs
Z_axis  = zaxis*1.e+2 # cgs
fpsiar=fpsiar*1.e+6   # r * B_phi  [cgs now]
epsi=  epsi*1.e+8     # cgs now
psilim= psilim*1.e8     
psimag= psimag*1.e8 
 
psisgn=+1. # for ascending psi(rho)
if psilim<psimag:
    psisgn=-1.  # descending psi(rho)

btor=btor*1e4 # vacuum toroidal magnetic field at radmaj [Gauss]

# Form R,Z grids
dzz=zbox/(nz-1) # cgs
drr=rbox/(nr-1) # cgs
er=np.zeros(nr)
ez=np.zeros(nz)
er[0]=rboxdst    # cgs
ez[0]= zaxis -zbox*.5
for nn in range(1,nr,1):  # nn goes from 1 to nr-1
    er[nn]=er[nn-1]+drr   # R-grid [cm] 

for nn in range(1,nz,1):  # nn goes from 1 to nz-1
    ez[nn]=ez[nn-1]+dzz   # Z-grid [cm]

R,Z = np.meshgrid(er,ez) # 2D grids [cm]

ez_mn= min(ez)
ez_mx= max(ez)
er_mn= min(er)
er_mx= max(er)
print 'drr,dzz [cm] =',drr,dzz
print ' min/max of er grid:', er_mn, er_mx, ' [cm]'
print ' min/max of ez grid:', ez_mn, ez_mx, ' [cm]'

#........................................................
# Form the equally spaced psi array/grid.
# psimag < psilim;  epsi has min. at m.axis
delpsi= (psilim-psimag)/(nv-1)  
psiar=np.zeros(nv)
for ix in range(0,nv,1):  # ix goes from 0 to nv-1
    psiar[ix]= psimag + ix*delpsi # [psimag; psilim]
#.........................................................

# Grad(psi) components:
grpsi_z, grpsi_r = np.gradient(epsi,dzz,drr)
grpsi= sqrt(grpsi_r*grpsi_r + grpsi_z*grpsi_z)
Bpol= zeros((nz,nr))
Btor= zeros((nz,nr))
Bz=   zeros((nz,nr))
Bpmn=1.e10 # will be found below
Bpmx=0.    # will be found below
# Define |Bpol|  and  Btor  [Gauss]
for ir in range(0,nr,1):  # ir goes from 0 to nr-1
    if er[ir] > 1.e-5: # R>0 area (R=0 is calculated below)
        for iz in range(0,nz,1):  # iz goes from 0 to nz-1
            Bpol[iz,ir]= grpsi[iz,ir]/er[ir]
            Bz[iz,ir]= grpsi_r[iz,ir]/er[ir]
            if epsi[iz,ir]*psisgn<psilim*psisgn:   # inside LCFS (plasma)
                iv= int((epsi[iz,ir]-psimag)/delpsi)
                iv= max(iv,0) # make sure iv is not negative
                #print 'ir,iz,iv=',ir,iz,iv
                #if iz==64: print epsi[iz,ir], iv, psiar[iv],psiar[iv+1]             
                Btor[iz,ir]=fpsiar[iv]/er[ir]
                Bpmn=min(Bpol[iz,ir],Bpmn)
                Bpmx=max(Bpol[iz,ir],Bpmx)
            else:   # outside LCFS (nearly vacuum)
                Btor[iz,ir]=btor*radmaj/er[ir]
            
        
# Now consider R=0, if present
if er[0]<=1.e-5:  # R-grid starts at R=0. (This is a typical FRC case)
    for iz in range(0,nz,1):   # iz goes from 0 to nz-1
        Btor[iz,0]= 0. # Set Btor(R=0) to 0.
        # Bpol at R=0 is calculated as 2*d(psi)/d(R^2)
        Bz[iz,0]= 2*(epsi[iz,1]-epsi[iz,0])/(drr*drr)
        Bpol[iz,0]= abs(Bz[iz,0])
        Bpmn=min(Bpol[iz,0],Bpmn)
        Bpmx=max(Bpol[iz,0],Bpmx)

print ' min/max of Bpol:', Bpol.min(), Bpol.max(), ' [Gauss]'
print ' min/max of Btor:', Btor.min(), Btor.max(), ' [Gauss]'
if Rmx_plots==0:
    Rmx_plots=er[nr-1]*1.2
if Zmx_plots==0:
    Zmx_plots=ez[nz-1]
if Zmn_plots==0:
    Zmn_plots=ez[0]
print ' Rmx_plots:', Rmx_plots, ' [cm]'


#--------------------------------------------------------------------------
fig0=plt.figure() #   PSI(R) and Bz(R) along Z=Zmidplane
nzmid= int(nz/2)  # index for the midplane
Zmid= ez[nzmid]
if PSImn_plots==0:
    PSImin=psimag/1e8
else:
    PSImin=PSImn_plots
    
if PSImx_plots==0:
    PSImax=psilim/1e8
else:
    PSImax=PSImx_plots
    
#--------------------------
ax = plt.subplot(222)
#plt.xlabel('$R$  $(cm)$')
plt.title('$B_{z}$  $(Tesla)$' '  $along$ $Z=$' + r"$%1.2f$"%(Zmid) + '$cm$')
plt.grid(True)
plt.hold(True)
plot(er,Bz[nzmid,:]/1e4)
xlim( (0.0, Rmx_plots) )
#-------------------------
ax = plt.subplot(224)
plt.xlabel('$R$  $(cm)$')
plt.title('$\Psi_p/2\pi$   $(Wb/rad)$')
plt.grid(True)
plt.hold(True)
plot(er,epsi[nzmid,:]/1e8)
xlim( (0.0, Rmx_plots) )
ylim( (PSImin, PSImax) )
#-------------------------
if eqdsktype=='TAE': # data on electrons only
    ax= plt.subplot(223)
    plt.xlabel('$R$  $(cm)$')
    plt.ylabel('$T_{e}$  $(eV)$' '  $along$ $Z=$'+r"$%1.2f$"%(Zmid)+'$cm$')
    plt.grid(True)
    plt.hold(True)
    plot(er,tmcut[nzmid,:],'b')
    xlim( (0.0, Rmx_plots) )
    #-----------------------------
    ax= plt.subplot(221)
    #plt.xlabel('$R$  $(cm)$')
    plt.ylabel('$n_{e}$  $(m^{-3})$' '  $along$ $Z=$'+r"$%1.2f$"%(Zmid)+'$cm$')
    plt.grid(True)
    plt.hold(True)
    plot(er,necut[nzmid,:],'b')
    xlim( (0.0, Rmx_plots) )

if eqdsktype=='TAEi':  # electrons, and additional data on ions
    ax= plt.subplot(223)
    plt.xlabel('$R$  $(cm)$')
    plt.ylabel('$T$  $(eV)$' '  $along$ $Z=$'+r"$%1.2f$"%(Zmid)+'$cm$')
    plt.grid(True)
    plt.hold(True)
    plot(er,tmcut[nzmid,:],     'b') # Te
    plot(er,ti_thermal[nzmid,:],'g') # Ti
    #plot(er,ti_fast[nzmid,:],   'r') # Ti fast/hot
    xlim( (0.0, Rmx_plots) )
    #-----------------------------
    ax= plt.subplot(221)
    #plt.xlabel('$R$  $(cm)$')
    plt.ylabel('$n$  $(m^{-3})$' '  $along$ $Z=$'+r"$%1.2f$"%(Zmid)+'$cm$')
    plt.grid(True)
    plt.hold(True)
    plot(er,necut[nzmid,:],     'b') # ne
    plot(er,ni_thermal[nzmid,:],'g') # ni
    #plot(er,ni_fast[nzmid,:],   'r') # ni fast/hot
    xlim( (0.0, Rmx_plots) )
    
savefig('eqdsk'+'_psi_B_along_midplane.png')
show() 
#--------------------------------------------------------------------------


#--------------------------------------------------------------------------
fig1=plt.figure() #  ne, Te, B, Psi, in R-Z view

#-------------------------
if eqdsktype=='TAE':
    ax = plt.subplot(141)
    ax.set_aspect(1.0)
    #ax.axis([40,60,-60,55]) # cm     Specify limits for the plot  
    ylim( (Zmn_plots, Zmx_plots) )   
    xlim( (0.0, Rmx_plots) )
    plt.xlabel('$R$  $(cm)$')
    plt.title('$n_e$   $(m^{-3})$')
    plt.grid(True)
    plt.hold(True)
    
    if ncontr > 0:
        plt.plot(rcontr*100,zcontr*100,'b',linewidth=2) # LCFS
        
    if nlimiter > 0:
        plt.plot(rlimiter*100,zlimiter*100,'k',linewidth=2) # limiter/chamber
    
    plot(R_axis,Z_axis,'k+') # plot '+' at magnetic axis
    
    CS=plt.contour(R,Z,necut,Ncont,linewidths=linw,cmap=plt.cm.jet)
    CB=plt.colorbar(orientation='vertical', shrink=0.4,format='%.1e') # colorbar
    #------------------------
    ax = plt.subplot(142)
    ax.set_aspect(1.0)
    #ax.axis([40,60,-60,55]) # cm     Specify limits for the plot  
    ylim( (Zmn_plots, Zmx_plots) ) 
    xlim( (0.0, Rmx_plots) )
    plt.xlabel('$R$  $(cm)$')
    plt.title('$T_e$   $(eV)$')
    plt.grid(True)
    plt.hold(True)
    
    if ncontr > 0:
        plt.plot(rcontr*100,zcontr*100,'b',linewidth=2) # LCFS
        
    if nlimiter > 0:
        plt.plot(rlimiter*100,zlimiter*100,'k',linewidth=2) # limiter/chamber
    
    plot(R_axis,Z_axis,'k+') # plot '+' at magnetic axis
    
    CS=plt.contour(R,Z,tmcut,Ncont,linewidths=linw,cmap=plt.cm.jet)
    CB=plt.colorbar(orientation='vertical', shrink=0.4,format='%.1e') # colorbar
    #------------------------    

#.....................
ax = plt.subplot(143)
ax.set_aspect(1.0)
#ax.axis([40,60,-60,55]) # cm     Specify limits for the plot
ylim( (Zmn_plots, Zmx_plots) )  
xlim( (0.0, Rmx_plots) )    
plt.xlabel('$R$  $(cm)$')
plt.title('$|B_{total}|$  $(Tesla)$')
plt.grid(True)
plt.hold(True)

if ncontr > 0:
    plt.plot(rcontr*100,zcontr*100,'b',linewidth=2) # LCFS
    
if nlimiter > 0:
    plt.plot(rlimiter*100,zlimiter*100,'k',linewidth=2) # limiter/chamber

plot(R_axis,Z_axis,'k+') # plot '+' at magnetic axis

B=sqrt(Btor*Btor + Bpol*Bpol) # Gauss
D=B/1e4 # to [Tesla]
if Bmn_plots==0:
    Dmin=np.min(D)
else:
    Dmin=Bmn_plots
    
if Bmx_plots==0:
    Dmax=np.max(D)
else:
    Dmax=Bmx_plots
    
print 'Bmin[T]=',Dmin, '  Bmax[T]=',Dmax
if ((Dmax>1000*Dmin) & (Dmin>1e-8)): 
    # if Bmin is at least 1 Gauss, check range in |B|, 
    # reduce to avoid peaks near wires (of magnetic coils, if present)
    Dmin= np.amin(abs(Bz[nzmid,:]))/1e4
    Dmax=Dmin*300
    print 'For plots of rays over |B| levels, Bmin and Bmax are adjusted to:'
    print 'Bmin[T]=',Dmin, '  Bmax[T]=',Dmax

levels_B=np.arange(Dmin,Dmax,(Dmax-Dmin)/(Ncont-1))
CS=plt.contour(R,Z,D,levels_B,linewidths=linw,cmap=plt.cm.jet)
CB=plt.colorbar(orientation='vertical', shrink=0.4,format='%.1e') # colorbar

#-----------------------
ax = plt.subplot(144)
ax.set_aspect(1.0)
#ax.axis([40,60,-60,55]) # cm     Specify limits for the plot  
ylim( (Zmn_plots, Zmx_plots) )   
xlim( (0.0, Rmx_plots) )
plt.xlabel('$R$  $(cm)$')
plt.title('$\Psi_p/2\pi$   $(Wb/rad)$')
plt.grid(True)
plt.hold(True)

if ncontr > 0:
    plt.plot(rcontr*100,zcontr*100,'b',linewidth=2) # LCFS
    
if nlimiter > 0:
    plt.plot(rlimiter*100,zlimiter*100,'k',linewidth=2) # limiter/chamber

plot(R_axis,Z_axis,'k+') # plot '+' at magnetic axis
if PSImn_plots==0:
    PSImin=psimag/1e8
else:
    PSImin=PSImn_plots
    
if PSImx_plots==0:
    PSImax=psilim/1e8
else:
    PSImax=PSImx_plots

levels_B=np.arange(Dmin,Dmax,(Dmax-Dmin)/(Ncont-1))
CS=plt.contour(R,Z,D,levels_B,linewidths=linw/2,color='k')
    
levels=np.arange(PSImin,PSImax,(PSImax-PSImin)/(Ncont-1))
CS=plt.contour(R,Z,epsi/1e8,levels,linewidths=linw,cmap=plt.cm.jet)
CB=plt.colorbar(orientation='vertical', shrink=0.4,format='%.1e') # colorbar




savefig('eqdsk_RZ'+'_psi_B.png')
show() 
#--------------------------------------------------------------------------


#--------------------------------------------------------------------------
fig2=plt.figure() # Btor and Bpol in R-Z view

ax = plt.subplot(122)
ax.set_aspect(1.0)
xlim( (0.0, Rmx_plots) )
ylim( (Zmn_plots, Zmx_plots) )        
plt.xlabel('$R$  $(cm)$')
plt.title('$B_{tor}$  $(Tesla)$')
plt.grid(True)
plt.hold(True)

if ncontr > 0:
    plt.plot(rcontr*100,zcontr*100,'b',linewidth=2) # LCFS
    
if nlimiter > 0:
    plt.plot(rlimiter*100,zlimiter*100,'k',linewidth=2) # limiter/chamber

plot(R_axis,Z_axis,'k+') # plot '+' at magnetic axis

if btor>0:
    CS=plt.contour(R,Z,Btor/1e4,Ncont,linewidths=linw,cmap=plt.cm.jet)
    CB=plt.colorbar(orientation='vertical', shrink=0.4,format='%.2e') # colorbar

#........................
ax = plt.subplot(121)
ax.set_aspect(1.0)
xlim( (0.0, Rmx_plots) )
ylim( (Zmn_plots, Zmx_plots) )      
plt.xlabel('$R$  $(cm)$')
plt.title('$B_{pol}$  $(Tesla)$')
plt.grid(True)
plt.hold(True)

if ncontr > 0:
    plt.plot(rcontr*100,zcontr*100,'b',linewidth=2) # LCFS
    
if nlimiter > 0:
    plt.plot(rlimiter*100,zlimiter*100,'k',linewidth=2) # limiter/chamber

plot(R_axis,Z_axis,'k+') # plot '+' at magnetic axis

dB=(Bpmx*0.99-Bpmn)/Ncont
levels=np.arange(Bpmn,Bpmx*0.99,dB)
CS=plt.contour(R,Z,Bpol/1e4,levels/1e4,linewidths=linw,cmap=plt.cm.jet)
CB=plt.colorbar(orientation='vertical', shrink=0.4,format='%.2e') # colorbar
savefig('eqdsk_RZ'+'_Bpol_Btor.png')
show() #--------------------------------------------------------------------------


