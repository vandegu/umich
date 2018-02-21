# This script produces a zoomed-in look at a portion of the CESM output. It produces the following:
#
# 1. Color fill of chosen geospatial oceanic data.
# 2. Streamlines of horizontal oceanic motion.
#

from mpl_toolkits.basemap import Basemap,shiftgrid,interp
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
from scipy.interpolate import griddata as griddata2
import numpy.ma as ma

# Which file to read, and from what vertical level?
lvl = 50
file = '../data/onexone/mini/MAA_B1850C4CN_f19_g16_cret_2x_sewall.pop.h.0850.nc.1x1.nc'
f = netCDF4.Dataset(file)

u = f.variables['UVEL'][0,lvl,:,:]
v = f.variables['VVEL'][0,lvl,:,:]
w = f.variables['WVEL'][0,lvl,:,:]
pt = f.variables['TEMP'][:]
s = f.variables['SALT'][:]
pd = f.variables['PD'][:]
iage = f.variables['IAGE'][:][0,lvl,:,:]
print(u.shape)

lat = f.variables['lat'][:]
lon = f.variables['lon'][:]
z_t = f.variables['z_t'][:]
x,y = np.meshgrid(lon,lat)

print(z_t[lvl]/100.0) # To display depth level in m
fig = plt.figure(figsize=[14,12])
ax = fig.add_subplot(111)

# Create basemap instance (choose your projection here):
mp = Basemap(llcrnrlon=-65.0,llcrnrlat=-10.0,urcrnrlon=20.0,urcrnrlat=60.0,
             resolution='i',projection='cass',lon_0=-35,lat_0=45.0)
#mp = Basemap(width=10000000,height=8000000,
#            resolution='l',projection='stere',
#            lat_ts=50,lat_0=-60.0,lon_0=-60.0)
#mp = Basemap(projection='nplaea',boundinglat=30,lon_0=270,resolution='l')


xxold,yyold = mp(x,y) # Projecting the old data locations. Use these to interpolate from. 

# Build new grid to interpolate to...
lons,lats,xxnew,yynew = mp.makegrid(350,350,returnxy=True)

# Project the data onto the new grid from the old grid.
unew = griddata2((xxold.ravel(),yyold.ravel()),u.ravel(),(xxnew,yynew), method = 'linear')
unew2 = ma.masked_where((unew > 100000), unew)
vnew = griddata2((xxold.ravel(),yyold.ravel()),v.ravel(),(xxnew,yynew), method = 'linear')
vnew2 = ma.masked_where((vnew > 100000), vnew)
iagenew = griddata2((xxold.ravel(),yyold.ravel()),iage.ravel(),(xxnew,yynew), method = 'linear')
iagenew2 = ma.masked_where((iagenew > 100000), iagenew)
landmask = iagenew > 100000

# Create streamplot:
sc1 = mp.scatter(xxnew,yynew,c=iagenew2,edgecolor='None',s=10,vmin=0,vmax=1800,cmap='jet')
mp.streamplot(xxnew,yynew,unew2,vnew2,density=3,arrowsize=2,arrowstyle='-|>',color='black')
mp.scatter(xxnew[landmask],yynew[landmask],c='white',s=10,edgecolor='None',zorder=2)
mp.drawmeridians(np.arange(0,360,10),labels=[False,False,False,True],zorder=3)
mp.drawparallels(np.arange(-90,90,10),labels=[False,True,False,False],zorder=3)

layer = '%3.0fm'%(z_t[lvl]/100)
cb = plt.colorbar(sc1,orientation='horizontal',extend='both')
cb.set_label('$Ideal\/\/Age\/\/(yr)$',size=20)
cb.ax.tick_params(labelsize=16) 
ax.set_title('$Control:\/\/\/\/%s$'%layer,size=24)

plt.show()
