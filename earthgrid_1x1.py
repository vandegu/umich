# earthgrid_1x1.py

# Create rough area meshgrid for a latitude x longitude grid on the earth. Can be for any sphere, if R is specified.

class Earth(object):
    '''This class is used to analyze the output of CESM on a 1x1 latxlon grid.'''
    
    import numpy as np
    
    def __init__(self,):
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
        
    def area_on_sphere(self,lat1,lat2,lon1,lon2,R=6371.0):
        '''Returns the area on the Earth bounded by a set of lons and lats (entered in degrees).'''

        #R = 6371.0 # Radius of Earth; km.
        lat1,lat2,lon1,lon2 = np.deg2rad((lat1,lat2,lon1,lon2)) # Convert the lats,lons into radians.
        area = R**2*(lon2-lon1)*(np.sin(lat2)-np.sin(lat1))

        return area
    
    def quadrangle(self,lat0,lat1,lon0,lon1):

        lat = np.arange(-89.5,90.5,1.0)
        lon = np.arange(0.5,360.5,1.0)
        Lon,Lat = np.meshgrid(lon,lat)
        
        Latlow = Lat > lat0
        Lathigh = Lat < lat1
        goodLat = Latlow*Lathigh
        Lonlow = Lon > lon0
        Lonhigh = Lon < lon1
        goodLon = Lonlow*Lonhigh
        meshmask = goodLat*goodLon
    
        return meshmask
    
    def lonflip(self,lolon,hilon):
        
        lonbounds = [lolon,hilon]
        for i,x in enumerate(lonbounds):
            if x >= 0.0 and x < 180.0:
                lonbounds[i] += 180.0
            elif x >= 180.0 and x <= 360.0:
                lonbounds[i] -= 180.0
        newlolon = lonbounds[0]
        newhilon = lonbounds[1]
        
        return newlolon,newhilon

    
    def earthmask(self,lolat,hilat,lolon,hilon,pcen):
        '''Defines a mask that can be applied to a 2D spherical Earth 1x1 dataset to isolate specific areas.
           Input: Bottom, top, left, and right bounds of the quadrangle, plus pcen (is it pacific-centered?)
           Output: Meshgrid mask 180x360 to fit a 1x1 latxlon grid and only allow the quadrangle specified 
                   through.
        '''

        bounds = np.array([lolat,hilat,lolon,hilon])
        if bounds.any() < 0.0:
            raise ValueError('* * * Cannot submit negative coordinates to this function! Replace with 0-360 degrees.')

        if pcen: # With the PM at the edges of the map, in line with 180 longitude of the meshgrid.

            lolon,hilon = self.lonflip(lolon,hilon)

            if lolon > hilon:

                lowbox = self.quadrangle(lolat,hilat,lolon,360.0)
                highbox = self.quadrangle(lolat,hilat,0.0,hilon)
                box = np.maximum(lowbox,highbox) # Returns the element-wise maximum (in this case, 1 == True).

                return box

            else:

                return self.quadrangle(lolat,hilat,lolon,hilon)

        else: # With the PM in the center of the map, in line with the 0 longitude of the meshgrid. 

            if lolon > hilon:

                lowbox = self.quadrangle(lolat,hilat,lolon,360.0)
                highbox = self.quadrangle(lolat,hilat,0.0,hilon)
                box = np.maximum(lowbox,highbox) # Returns the element-wise maximum (in this case, 1 == True).

                return box

            else:

                return self.quadrangle(lolat,hilat,lolon,hilon)

