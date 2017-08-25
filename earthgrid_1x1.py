# earthgrid_1x1.py

# Create rough area meshgrid for a latitude x longitude grid on the earth. Can be for any sphere, if R is specified.

class Earth(object):
    
    def __init__(self,R=6371.0):
        self.R = R # Radius of Earth (km): 6371.0 km
        self.margins = self.latlon_grid_boundaries()
        A = self.area_on_sphere(self.margins[:,0],self.margins[:,1],self.margins[:,2],self.margins[:,3])
        self.A = np.reshape(A,(180,360))
    
    def latlon_grid_boundaries(self,):
        '''This method produces the four boundary points that serve as input for the area on a sphere method.'''

        margins = np.empty((360*180,4))
        for x,lat in enumerate(np.arange(-90,90)):
            for y,lon in enumerate(np.arange(0,360)):
                margins[x*360+y,0] = lat
                margins[x*360+y,1] = lat+1
                margins[x*360+y,2] = lon
                margins[x*360+y,3] = lon+1

        return margins
        
    def area_on_sphere(self,lat1,lat2,lon1,lon2):
        '''Returns the area on the Earth bounded by a set of lons and lats (entered in degrees).'''

        lat1,lat2,lon1,lon2 = np.deg2rad((lat1,lat2,lon1,lon2)) # Convert the lats,lons into radians.
        area = self.R**2*(lon2-lon1)*(np.sin(lat2)-np.sin(lat1))

        return area
    
    def earthmask(self,lolat,hilat,lolon,hilon):
        '''Defines a mask that can be applied to a 2D spherical Earth 1x1 dataset to isolate specific areas.'''

        lat = np.arange(-89.5,90.5,1.0)
        lon = np.arange(0.5,360.5,1.0)
        Lon,Lat = np.meshgrid(lon,lat)
        Latlow = Lat > lolat
        Lathigh = Lat < hilat
        goodLat = Latlow*Lathigh
        Lonlow = Lon > lolon
        Lonhigh = Lon < hilon
        goodLon = Lonlow*Lonhigh
        meshmask = goodLat*goodLon

        return Lon[meshmask],Lat[meshmask]
