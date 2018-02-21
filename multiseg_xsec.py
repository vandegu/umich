from mpl_toolkits.basemap import Basemap,shiftgrid,interp
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
from scipy.interpolate import griddata as griddata2
import numpy.ma as ma

class multi_segment_xsection(object):
    '''This class produces a continuously-interpolated cross-section through the ocean of CESM of any
       continuous line drawn across the globe, but one that does not need to be confined to any 
       particular latitude or longitude line, nor does it need to be a continuously-differentiable
       line.
       
       To produce a plot, proceed with the following:
           
        1. Create an instance of multi_segment_xsection with the desired arguments (see __init__ below)
        2. Call interpolate_grid().
        3. Call makeplot(vmin,vmax), where the arguments are the minimum and maximum fill values,
            respectively.
        4. (optional) Call insetplot() to place a reference plot in the corner of the image. This 
            reference plot will show exactly where the cross section is, and where the alphabet 
            labels are located. It is useful to show your readers where your data was pulled from
            in the dataset.
            
        Andrew Vande Guchte -- February 2, 2018
    '''

    
    def __init__(self,lats,lons,data,dlats,dlons,ddepth):
        '''The arguments passed are a list of latitudes and respective longitudes of the vertices of
           the multi-segmented cross-section. Then, pass the data desired to be plotted, and the 
           lats, lons, and depths of its location (must be 1d arrays for where the data is currently 
           located).'''
        
        self.lats = np.deg2rad(np.array(lats))
        self.lons = np.deg2rad(np.array(lons))
        self.create_plane()
        self.alph = ['A','B','C','D','E','F']
        
        self.data = data
        self.dlats = dlats
        self.dlons = dlons
        self.ddepth = ddepth
        
    
    def surface_distance(self,lats,lons):
        '''This method calculates the haversine distance between two points.
           Must be given in radians.'''
        
        if lats.size <= 1:
            raise SystemExit('surface_distance requires at least (2) latitudes and longitudes')
            
        if lats.size != lons.size:
            raise SystemExit('surface_distance requires that there be the same number of lons and lats')
        
        surface_dist = []
        
        for i,qqq in enumerate(lats[:-1]):
            a = (np.sin((lats[i+1]-lats[i])/2)**2 + np.cos(lats[i])*np.cos(lats[i+1]) * 
                 np.sin((lons[i+1]-lons[i])/2)**2)
            c = 2.0 * np.arctan2(np.sqrt(a),np.sqrt(1.0-a))
            d = 6371.0*c
            surface_dist.append(d)
        
        return np.array(surface_dist)
    

    def create_plane(self,):
        
        scaling_factor = 1000.0
        td = self.surface_distance(self.lats,self.lons)#*scaling_factor
        td = td.astype(int)
        self.td = td
        
        self.plats = np.empty(0)
        self.plons = np.empty(0)
        
        for i,qqq in enumerate(self.lats[:-1]):
            self.plats = np.hstack((self.plats,np.linspace(self.lats[i],self.lats[i+1],td[i])))
            self.plons = np.hstack((self.plons,np.linspace(self.lons[i],self.lons[i+1],td[i])))
        
        # convert the new lat and lon arrays back to degrees:
        
        self.plats = np.rad2deg(self.plats)
        self.plons = np.rad2deg(self.plons)
        
        #print(self.plats)
        #print(self.plons)
    

    def interpolate_grid(self,):#data,dlats,dlons,ddepth):
        '''This method will interpolate to the planar grid from each level of the model.
           The ddepth should be in meters.'''
        
        # first create the lats and lons of the multi-segmented panel that we want.
        self.create_plane()
        
        # create meshgrid of locations of the data on the surface of the earth (lons/lats) from which 
        # to interpolate once it is projected.
        x,y = np.meshgrid(self.dlons,self.dlats)
        
        # create instance of the basemap projection from which to interpolate meaningfully and efficiently.
        m = Basemap(projection='eck4',lon_0=0,resolution='c')
        
        # create the projected lon/lat of the actual data. xxold and yyold are where the 
        # ACTUAL data exists and where it will be interpolated FROM.
        xxold,yyold = m(x,y) 
        
        # create the projected lon/lat of the new planar data. This will be 1 dimensional because it only
        # represents the projection of our plane through the ocean at a horizontal slice.
        xxnew,yynew = m(self.plons,self.plats)
        print(xxnew.shape)
        
        plotdata = np.empty(0)
        
        # interpolate at each vertical level in the model to the desired plane.
        for i,z in enumerate(self.ddepth):
            
            print(i)
            datanew = griddata2((xxold.ravel(),yyold.ravel()),self.data[i,:,:].ravel(),(xxnew,yynew),method='linear')
            datanew = ma.masked_where((datanew>100000),datanew)
            plotdata = np.hstack((plotdata,datanew))
        
        # reshape the plotdata to be in the correct array shape.
        plotdata = plotdata.reshape(plotdata.size/self.plons.size,self.plons.size)
        
        # mask the plotdata where there is land. 
        plotdata = ma.masked_where((plotdata>100000),plotdata)
        
        self.plotdata = plotdata
        
        # create a more realistic grid from which to interpolate to. The xxplot and zzplot (actually a z-axis
        # for the data) are created here and used to interpolate the data to and then to plot with.
        self.xx = np.arange(0,self.plons.size,10)
        self.zz = np.arange(0,5.500,.01)
        xxplotold,zzplotold = np.meshgrid(np.arange(self.plons.size),self.ddepth/100000.0) # make sure both are in km. 
        xxplotnew,zzplotnew = np.meshgrid(self.xx,self.zz) # both are in km.
        
        #print(xxplotold,zzplotold,xxplotnew,zzplotnew)
        
        # interpolate to finer resolution grid within the desired plane.
        paneldata = griddata2((xxplotold.ravel(),zzplotold.ravel()),
                              plotdata.ravel(),(xxplotnew,zzplotnew),
                              method='linear')
        
        paneldata = ma.masked_where((paneldata>100000),paneldata)
        self.paneldata = paneldata
        print(paneldata.shape)    
        
    
    def makeplot(self,vmin=0,vmax=1800):
        
        fig = plt.figure(figsize=[16,13])
        ax = fig.add_subplot(111)
        ax.set_xlim(0,self.xx[-1])
        ax.set_ylim(self.zz[-1],0)
        pc = ax.pcolor(self.xx,self.zz,self.paneldata,vmin=vmin,vmax=vmax)
        
        ax.set_ylabel('Depth (km)',size=20)
        #ax.set_xlabel('')
        
        ax.set_xticks([np.sum(self.td[:i]) for i in range(self.td.size+1)])
        ax.set_xticklabels([self.alph[i] for i in range(self.td.size+1)])
        ax.tick_params(labelsize=18)
        
        cb = plt.colorbar(pc,orientation='horizontal',extend='max')
        cb.set_label('Ideal Age (yr)',size=22)
        cb.ax.tick_params(labelsize=16)
        
        self.ax = ax
        

    def insetplot(self,pt,lonflip=False):
        '''Create an inset plot that shows where the cross-section is. Be sure to include the top level
           of potential temperature (pt) here as the only argument passed. Technically, you could use any
           level of pt, but whichever level you choose will be the one that is used to find the ocean and
           continents in the inset plot.'''
        
        if lonflip:
            pt = np.hstack((pt[:,(pt.shape[1]/2):],pt[:,:(pt.shape[1]/2)]))
        
        axin = plt.axes([.14, .35, .15, .12])
        xinset,yinset = np.meshgrid(self.dlons,self.dlats)
        m = Basemap(ax=axin,projection='eck4',lon_0=0,resolution='c')
        xxinsetold,yyinsetold = m(xinset,yinset) 
        lons,lats,xxinsetnew,yyinsetnew = m.makegrid(500,500,returnxy=True)
        ptnew = griddata2((xxinsetold.ravel(),yyinsetold.ravel()),pt.ravel(),(xxinsetnew,yyinsetnew),method='linear')
        l = ptnew>100000 # find where land exists
        s = ptnew<100000 # find where sea exists
            
        sc1 = m.scatter(xxinsetnew[l],yyinsetnew[l],c='coral',edgecolor='None',s=.6)
        sc1 = m.scatter(xxinsetnew[s],yyinsetnew[s],c='aqua',edgecolor='None',s=.6)

        # project the coordinates of the vertices of the multi-segmented panel:
        if lonflip:
            linex,liney = m(np.rad2deg(self.lons)+180.0,np.rad2deg(self.lats))
        else:
            linex,liney = m(np.rad2deg(self.lons),np.rad2deg(self.lats))
        
        m.plot(linex,liney,'r-',lw=2)

        for i in range(len(linex)):
            plt.annotate(self.alph[i],xy=(linex[i]+1000000,liney[i]),size=14)


# Example of how to produce a plot:

if __name__=='__main__':

    file = '../data/onexone/mini/MAA_B1850C4CN_f19_g16_cret_4x_sewall.pop.h.0500.nc.1x1.nc'
    f = netCDF4.Dataset(file)
    iage = f.variables['IAGE'][0,:,:,:]
    lat = f.variables['lat'][:]
    lon = f.variables['lon'][:]
    z_t = f.variables['z_t'][:]
    a = multi_segment_xsection([35,-15,-65],[-350,-240,-100],iage,lat,lon,z_t) # Tethys
    a.interpolate_grid()

    a.makeplot(vmin=0,vmax=1800)
    pt = f.variables['TEMP'][0,0,:,:]
    a.insetplot(pt,lonflip=1)

    a.ax.set_ylim(5.0,0)

    a.ax.set_title('Control: Water Age',size=28)

    plt.show()
