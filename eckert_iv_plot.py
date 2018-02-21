# This script will produce:
#
# 1. Colorfill of desired oceanic data in Eckert IV projection (Equal-Area, but with some distortion; easy to look at, however).
# 2. Masked land, so only ocean data shows up.
#

from mpl_toolkits.basemap import Basemap,shiftgrid,interp
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
from scipy.interpolate import griddata as griddata2
import numpy.ma as ma


# Which file and what vertical level do you wish to plot from?
lvl = 55
file = '../data/onexone/mini/MAA_B1850C4CN_f19_g16_cret_4x_sewall.pop.h.0500.nc.1x1.nc'
f = netCDF4.Dataset(file)

#print(f.variables)
u = f.variables['UVEL'][0,lvl,:,:]
v = f.variables['VVEL'][0,lvl,:,:]
w = f.variables['WVEL'][0,lvl,:,:]
pt = f.variables['TEMP'][0,lvl,:,:]
s = f.variables['SALT'][:][0,lvl,:,:]
pd = f.variables['PD'][:]
iage = f.variables['IAGE'][0,lvl,:,:]
print(u.shape)

lat = f.variables['lat'][:]
lon = f.variables['lon'][:]
z_t = f.variables['z_t'][:]
x,y = np.meshgrid(lon,lat)

print(z_t[lvl]/100.0) # To display depth level in m
fig = plt.figure(figsize=[18,12])
ax = fig.add_subplot(111)

#m = Basemap(projection='mbtfpq',lon_0=0,resolution=None)
m = Basemap(projection='eck4',lon_0=0,resolution='c')
xxold,yyold = m(x,y) 

# generate a grid that is equally spaced in a plot with the current pojection
lons,lats,xxnew,yynew = m.makegrid(500,500,returnxy=True)
#print(lon.shape)
#print(xx.flatten().shape,yy.flatten().shape,u.flatten().shape)

# project the data onto the new grid
iagenew = griddata2((xxold.ravel(),yyold.ravel()),pt.ravel(),(xxnew,yynew), method = 'linear')
iagenew = ma.masked_where((iagenew > 100000), iagenew)
unew = griddata2((xxold.ravel(),yyold.ravel()),u.ravel(),(xxnew,yynew), method = 'linear')
unew = ma.masked_where((unew > 100000), unew)
vnew = griddata2((xxold.ravel(),yyold.ravel()),v.ravel(),(xxnew,yynew), method = 'linear')
vnew = ma.masked_where((vnew > 100000), vnew)

sc1 = m.scatter(xxnew,yynew,c=iagenew,edgecolor='None',s=5,cmap='jet')
#m.streamplot(xxnew,yynew,unew,vnew,density=3,arrowsize=2,arrowstyle='-|>',color='black')

m.drawmeridians(np.arange(0,360,30))
m.drawparallels(np.arange(-90,90,30))
layer = '%3.0fm'%(z_t[lvl]/100)
cb = plt.colorbar(orientation='horizontal',extend='both')
cb.set_label('$Potential\/\/Temperature\/\/(K)$',size=20)
cb.ax.tick_params(labelsize=16) 
ax.set_title('$MAA\/\/4x\/\/PI\/\/CO_2:\/\/\/\/%s$'%layer,size=24)

plt.show()
