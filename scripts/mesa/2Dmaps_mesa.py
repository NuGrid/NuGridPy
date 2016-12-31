#!/usr/bin/env python
###############################################################################
## INFO
## example: ./2Dmaps_mesa.py ../mesa_run/LOGS
##
## Script developed by R. Hirschi, Samuel Jones and Jacqueline den Hartogh
###############################################################################

import sys
import matplotlib.pyplot as mpl
import matplotlib.colors
import os 
import numpy as np 
from matplotlib.font_manager import fontManager, FontProperties
from scipy import interpolate
#import mesa as ms
from NuGridPy import mesa as ms

format = '.png'
##----------------------------------plots--------------------------------------#
## set plotting parameters
params = {'backend': 'Qt4Agg',
          'axes.labelsize':  10,
          'text.fontsize':   10,
          'legend.fontsize': 9,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'text.usetex': True,
	  'figure.dpi': 100}
mpl.rcParams.update(params)

#          'text.usetex': True,
###############################################################################
## options
if len(sys.argv)<1:
  print "****************************************************************"
  print "Not enough arguments! ex: ./2Dmaps_mesa.py ../mesa_run/LOGS"
  print "****************************************************************"
  sys.exit()

path = str(sys.argv[1])
ptitle = 'mesa'

modelnr_start = 35000 #norot
modelnr_end = 41000 #norot
#modelnr_start = 39350 #norot
#modelnr_end = 43000 #norot
#modelnr_start = 33100 #rot
#modelnr_end = 38500 #rot 
#modelnr_start = 48000 #43000 # 
#modelnr_end =  50000 #45000 #
modelnr_step = 20



## note: the line below requires that (modelnr_end-modelnr_start) can be divided by modelnr_step
ifile=(modelnr_end-modelnr_start)/modelnr_step


last_profile=ms.mesa_profile(path,modelnr_end,num_type='nearest_model')
mini = last_profile.header_attr.get('initial_mass')
#print 'mini= ', mini

outfile = '2Dmaps'+'mesa'+format
fig1 = mpl.figure()
ax=mpl.subplot(111)
mpl.suptitle(ptitle)

## Choose the x-axis: age, log_time_left, model_number
ixaxis = 2
if ixaxis == 0:
  print 'plot versus age'
  xxaxis = 'star_age'
  xxlabel = 'Age [yr]'
elif ixaxis == 1:
  print 'plot versus time left'
#  current_Xvar="ageadv"
  xxaxis = 'star_age'
  agemin = 1.e-5
  xxlabel = "log$_{10}$(Time left until last model) [yr]"
elif ixaxis == 2:
  print 'plot versus model number'
  xxaxis = 'model_number'
  xxlabel = "Model number"
elif ixaxis == 3:
  print 'plot versus (age - age_reset) [yr]'
  xxaxis = "star_age"
  xxlabel = 'Age - age$_{\mathrm{ref}}$ [yr]'
  fileagereset=ms.mesa_profile(path,41000,num_type='nearest_model')

## Choose the y-axis quantity: Mass fraction, M_r/Mo, log(r_cm), shell number
iyaxis = 1
if iyaxis == 0:
  print 'plot versus Mr/Mtot'
  yyaxis = 'mass'
  yylabel = "$M_r/M_{\mathrm{total}}$"
elif iyaxis == 1:
  print 'plot versus Mr [Mo]'
  yyaxis = 'mass'
  yylabel = "$M_r$ [M$_\odot$]"
elif iyaxis == 2:
  print 'plot versuslog$_{10}$(radius) [cm]'
  yyaxis = 'logR'
  yylabel = 'log$_{10}$(Radius) [cm]'
elif iyaxis == 3:
  print 'plot versus shell number'
  yyaxis = 'zone'
  yylabel = "Shell number"
#  ax.invert_yaxis()
elif iyaxis == 4:
  print 'plot versus radius [cm]'
  yyaxis = "radius"
  yylabel = 'Radius [$R_{\odot}$]'
  yminlim = 0
  ymaxlim = 1 #500.

## Choose the z-axis quantity (convection plotted by default): energy, abundances, ...
## zmin is minimum for plots of relevant quantity
izaxis = 2
if izaxis == 0:
  print 'plot energy generation'
  zzaxis = 'eps_nuc'
  zmin = 1.0e0
  title = 'Energy generation contours (red/contours: heating, blue: cooling, grey: convection)\n' + path
elif izaxis == 1:
  print 'plot abundance'
  zzaxis = 'temperature' #'eps_nuc' #'o16' #'gradL_composition_term' #'am_log_D_GSF' # 'c12' #'gradL_composition_term' #'log_D_ovr' #'cno' #'log_D_ovr' #'c13' #'am_log_D_GSF' #'omega' #'log_D_mix' # 'c13' #radius' #'gradr' #'am_log_nu_rot' #'log_D_mix' #'gradr' 
  zmin = 1.e-1
  title = zzaxis +' ' ' abundance contours (reds; grey: convection)\n' + path
elif izaxis == 2:
  print 'plot abundance or qtty varying linearly or over short range'
  zzaxis = 'temperature' #'tri_alfa' #'c12' #'am_log_D_GSF' #'gradL_composition_term' #'j_rot' #'mass' #'log_D_mix' #'j_rot' #'omega' #'logR' #'omega' #'gradr' #'grada' #'he4' #'mixing_type'
  zmin = 100
  title = zzaxis + ' contours (reds; grey: convection)\n' + path
elif izaxis == 3: #(1 minus 2)
  print 'plot abundance'
  zzaxis1 = 'c13' #'gradr' 'log_D_mix' #
  zzaxis2 = 'n14' #'grada' 'log_D_ovr' #
  zmin = 1.e-5
  title = zzaxis1 + '-' + zzaxis2 + ' abundance contours (reds; grey: convection)\n' + path
elif izaxis == 4: #(1 minus 2 minus 3)
  print 'plot abundance'
  zzaxis1 = 'log_D_mix' #'c13' #'gradr' 
  zzaxis2 = 'log_D_ovr' #'n14' #'grada' 
  zzaxis3 = 'log_D_conv' #'n14' #'grada' 
  zmin = 1.e-8
  title = zzaxis1 + '-' + zzaxis2 + '-' + zzaxis3 +' abundance contours (reds; grey: convection)\n' + path
elif izaxis == 5: #(1 + 2 + 3 + 4)
  print 'plot abundance'
  zzaxis1 = 'am_log_D_ES' #'c13' #'gradr' 
  zzaxis2 = 'am_log_D_GSF' #'n14' #'grada' 
  zmin = 1.e-8
  title = zzaxis1 + '+' + zzaxis2 +' abundance contours (reds; grey: convection)\n' + path



  ### where we get('neut')if zzaxis== 'neut': log(Nn)-8

	
#levels2 = [0.05,0.06,0.07,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.175,0.2]
#clabel2 = "$\Omega$ [Hz]"
clabel0 = "$\epsilon_{\mathrm{nuc}}-\epsilon_{\mathrm{nu}}$ [erg \,g$^{-1}$\,s$^{-1}$]"

last_profile=ms.mesa_profile(path,modelnr_end,num_type='nearest_model')
first_profile=ms.mesa_profile(path,modelnr_start,num_type='nearest_model')

## age/modnb of star at moment of last modelnr included in plot:
xxmax = last_profile.header_attr.get(xxaxis)
#print 'xxmax= ', xxmax, 'mini= ',mini

## setting x & y axes limits:
if ixaxis == 1:
  xmax = np.log10(agemin/2.)
  xmin = np.log10(xxmax*2.)
elif ixaxis==3:
  age_reset = fileagereset.header_attr.get(xxaxis)
  print 'age_reset= ',age_reset
  xmax = (xxmax-age_reset)#*24.*365.*1.05
  xmin = -xmax
  #xmin = -10.
else:
  xmax = last_profile.header_attr.get(xxaxis)+100
  xmin = first_profile.header_attr.get(xxaxis)-100

if iyaxis == 2: 
  ## yy= log radius in cm
  yy = last_profile.get('logR')+np.log10(6.955e10)
elif iyaxis == 0: 
  yy = last_profile.get(yyaxis)/last_profile.header_attr.get('star_mass')
  print 'yy', yyaxis, yy
else: 
  yy = 1.*last_profile.get(yyaxis)

if iyaxis ==2:
  ymin = np.min(yy)-1.
  ymax = np.max(yy)+1.
elif iyaxis ==1:  #### mass limits
  ymin = 0.52 #0.58 #
  ymax = 0.57 #0.60 #
  mpl.axis([xmin, xmax, 0., mini])
elif iyaxis == 3:
  ymin = 500.
  ymax = 700.
elif iyaxis == 4:
  ymin = np.min(yy)*.95
  ymax = np.max(yy)*1.05
  if yminlim > ymin: 
    ymin = yminlim
  if ymaxlim < ymax:
    ymax = ymaxlim
else:
  ymin = 0.
  ymax = 1.
print ymin,ymax
  
mpl.axis([xmin, xmax, ymin, ymax])
mpl.xlabel(xxlabel)
mpl.ylabel(yylabel)

## colour-map info
## y-axis resolution
ny=2000
dy=(ymax-ymin)/float(ny)

## x-axis resolution
## chosen by profile sampling
## nx must be equal to the number of profiles to be plotted
nx = ifile
dx=1

y = np.arange(ymin, ymax, dy)
x = np.arange(0, ifile, dx)

#print y,ny,x,nx

xcontours=np.zeros(len(x),float)
mtot=np.zeros(len(x),float)
mtot[0]=mini
mrplot=np.zeros(len(x),float)
ccont=np.zeros([len(y),len(x)],float)

zcont=np.zeros([len(y),len(x)],float)
zcont2=np.zeros([len(y),len(x)],float)
zcontm=np.zeros([len(y),len(x)],float)
zcontr=np.zeros([len(y),len(x)],float)


ifile=1
for i in range(modelnr_start,modelnr_end,modelnr_step):
## Now load profile using mesa.py
            profile=ms.mesa_profile(path,i,num_type='nearest_model')

            ## prepare grid x-y coordinates for colour map
            if iyaxis == 2: 
              yy = profile.get(yyaxis)+np.log10(6.955e10)
            elif iyaxis == 0: 
              yy = profile.get(yyaxis)/profile.header_attr.get('star_mass')
	      print 'yy', yyaxis, yy
            else: 
              yy = 1.*profile.get(yyaxis)
            if ixaxis == 1:
              xx = np.log10(xxmax -profile.header_attr.get(xxaxis)+agemin)
            elif ixaxis == 3:
              xx = (profile.header_attr.get(xxaxis)-age_reset)#*24.*365.
            else:
              xx = profile.header_attr.get(xxaxis)
            #print 'xx', xx
            xcontours[ifile-1]=xx
            #print ifile-1, xcontours[ifile-1],len(yy), ny
            
            ## load quantity to plot into zcontv (zcontv2 contains negative values of zcontv in positive array)
            if izaxis == 0: 
              zcontv = profile.get('eps_nuc')-profile.get('non_nuc_neu')
              zcontv2 = np.minimum(zcontv,-zmin)
	    elif izaxis == 3:
	      zcontv = profile.get(zzaxis1)-profile.get(zzaxis2)
	    elif izaxis == 4:
	      zcontv = profile.get(zzaxis1)-profile.get(zzaxis2)-profile.get(zzaxis3)
	    elif izaxis == 5:
	      zcontv = profile.get(zzaxis1)+profile.get(zzaxis2)
            else:
              zcontv = profile.get(zzaxis)

            ## use zmin to avoid range of values too large for colour map:
            zcontv= np.maximum(zcontv,zmin)


            ## interpolate quantity to map zcontv into grid xx-yy: results in zcont
            f = interpolate.interp1d(yy, zcontv, bounds_error=False, fill_value=0.)
            zcont[:,ifile-1] = f(y)
            ## zcont2 to plot negative energy colour map
            if izaxis == 0: 
              f2 = interpolate.interp1d(yy, zcontv2, bounds_error=False, fill_value=0.)
              zcont2[:,ifile-1] = f2(y)
            
            ## Plot grey area for convective zones
            zz = profile.get('mixing_type')
            #print 'ifile= ', ifile ,'zz: ',zz
            ii=0
            while ii < len(zz):
               ## ignore non-convective mixing types for now:
               ##
               if (zz[ii]>1.):
                 zz[ii]=0.
               #if (zz[ii]!= 0. and zz[ii]!= 1.):
            	 #print len(zz),ii, zz[ii], zz[ii-1], zz[ii+1], xx, yy[ii]
               ii+=1
            ## now zz=1. in convective zones, zz=0. outside

            ## interpolate grey area into grid xx-yy: results in ccont
            g = interpolate.interp1d(yy, np.minimum(zz,1.), bounds_error=False, fill_value=0.)
            ccont[:,ifile-1] = g(y)

            ## interpolate mass/radius to map zcontm/r into grid xx-yy: results in zcontm/r
            mcontv = profile.get('mass')
            fm = interpolate.interp1d(yy,mcontv, bounds_error=False)
            zcontm[:,ifile-1] = fm(y)
            mcontr = profile.get('logR')
            fr = interpolate.interp1d(yy, mcontr, bounds_error=False)
            zcontr[:,ifile-1] = fr(y)

                 
            ii=0
            while ii < len(zz):
               ## add a line for the total mass
               if (ii == 0):
                 mtot[ifile-1]=yy[ii]
               ## add dots for the convective zones edges
               if (zz[ii]==1.):
            	 if (ii > 0 and ii< len(zz)-1):
            	   ## checks if a zone is convective while a neighbour is not
                   if (zz[ii-1] != 1. or zz[ii+1] != 1.):
            	     #print len(zz),ii, zz[ii], zz[ii-1], zz[ii+1], xx, yy[ii]
            	     mpl.plot(xx,yy[ii], 'k.',alpha=0.2)            	   
               ii+=1
            ifile+=1

ifile -=1
print 'number of files used= ',ifile, 'now creating plot'

mpl.plot(xcontours,mtot, 'k')

## B&W colour map:
cmapbw = matplotlib.colors.ListedColormap(['w', 'k'])
## create grey areas for convective zones (alpha is transparency factor):
cont_c = ax.contourf(xcontours,y,ccont, cmap=mpl.get_cmap('Greys'), alpha=0.3,levels=[0.9,1.5])



## create colour maps for chosen quantity:
print max(zcontv),min(zcontv)
print zcontv
print xcontours
print y
if izaxis == 2 or izaxis == 4 or izaxis == 5:
##  cont_num = ax.contour(xcontours,y,zcont,levels2,locator=matplotlib.ticker.LinearLocator())
  cont_num = ax.contour(xcontours,y,zcont)#,locator=matplotlib.ticker.LinearLocator())
##  CB = mpl.colorbar(cont_num, ticks=levels2, orientation='horizontal', extend='both')
  CB = mpl.colorbar(cont_num, orientation='horizontal', extend='both')
##  CB.set_label(clabel2)
elif izaxis >=0 :
  cont_num = ax.contour(xcontours,y,zcont,locator=matplotlib.ticker.LogLocator(base=10.))
  CB = mpl.colorbar(cont_num, orientation='horizontal', extend='both')
## optional colour bar for the contour lines:
#  CB = mpl.colorbar(cont_num, orientation='horizontal', extend='both')
#  CB.set_label(clabel0)
## optional numerical values on contour lines:
##  mpl.clabel(cont_num, inline=0, fmt='%1.0e',fontsize=10)
## mpl.clabel(cont_num, inline=1, fmt="%r %%",fontsize=10)

## energy generation (or any quantity which contains both positive and negative values)
if izaxis == 0: 
## positive values
  cont_z = ax.contourf(xcontours,y,zcont, cmap=mpl.get_cmap('Reds'), alpha=0.3,locator=matplotlib.ticker.LogLocator())
## negative values
  cont_z2 = ax.contourf(xcontours,y,-zcont2, cmap=mpl.get_cmap('Blues'), alpha=0.3,locator=matplotlib.ticker.LogLocator())
## optional colour bars for these colour maps:
##  CBARBURN1 = mpl.colorbar(cont_z)
##  CBARBURN2 = mpl.colorbar(cont_z2)

## for quantities, which contain only positive values (or zeros)
elif izaxis == 1 or izaxis == 3 or izaxis == 4 or izaxis == 5:
  cont_z = ax.contourf(xcontours,y,zcont, cmap=mpl.get_cmap('Reds'), alpha=0.3,locator=matplotlib.ticker.LogLocator())
  CBARBURN1 = mpl.colorbar(cont_z)
elif izaxis == 2:
  cont_z = ax.contourf(xcontours,y,zcont, cmap=mpl.get_cmap('Reds'), alpha=0.3,locator=matplotlib.ticker.LinearLocator())
##  cont_z = ax.contourf(xcontours,y,zcont,levels2, cmap=mpl.get_cmap('Reds'), alpha=0.3,locator=matplotlib.ticker.LinearLocator())
  CBARBURN1 = mpl.colorbar(cont_z)

## mass/radius contours:
levelsm = [mini*0.001,mini*0.01,mini*0.1,mini*0.2,mini*0.3,mini*0.4,mini*0.5,mini*0.6,mini*0.7,mini*0.8,mini*0.9,mini]

#levelsr  = [7.,8.,9.,10.,11.,12.,13.,14.]
levelsr  = [-2,-1.88,-1.85,-1.75,-1.5,1.25,-1,-0.5,0,0.5,0.75,1.]#,1.25,1.5,1.75,2.,2.5] 

if iyaxis == 0 or iyaxis == 1:
  print "zontr"
  cont_r = ax.contour(xcontours,y,zcontr,levelsr,colors='k')#,locator=matplotlib.ticker.LogLocator(base=10.))
  #CB = mpl.colorbar(cont_r, orientation='horizontal', extend='both')
  mpl.clabel(cont_r,ticks=levelsr, fmt='logR = %3.3f', fontsize=9, inline=1)
elif  iyaxis == 2 or iyaxis == 4:
  print "zcontm"
  cont_m = ax.contour(xcontours,y,zcontm,levelsm,colors='k',locator=matplotlib.ticker.LinearLocator())
#  CB = mpl.colorbar(cont_m, orientation='horizontal', extend='both')
  mpl.clabel(cont_m,ticks=levelsm,fmt='$M_r$ = %3.2f $M_{\odot}$', fontsize=9, inline=1)


mpl.suptitle(title)

fig1.savefig(outfile)
print outfile,' is done'
mpl.show()
######################################################################


  

