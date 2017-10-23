# This script is to replace the missing values in a modified paleogeography.
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import os
import pyproj as proj4
import scipy.interpolate as si

# The CoordinateSystem and GeographicSystem classes are hereby used courtesy of Dr. Eric Bruning; his
# repository of excellent cartography tools (among other things) can be found at https://github.com/deeplycloudy.

# Developer notes: This code could be more elegant, and because it involves several loops through data, can be
# somewhat slow. The 'nhn' interpolation method is the slowest (~35 minutes), the 'nn' method follows, and the
# 'nvn' method is the quickest. A global averaging method may be very fast in the future when it is implemented.
# Finally, this code could be rewritten to calculate exactly which locations within the grids need to be modified
# and from where for each nearest neighbor interpolation scheme, and then that matrix could be passed to all of
# the necessary components of the code to speed it up. However, this improvement will have to wait, as I am too
# busy to implement it now.

# Last few IMPORTANT notes about running this script:
#
# 1) To change the input, output, and other model variables go to the very bottom of the script, under the
#    'if name = main' conditional. These should be the only thing you need to change, unless you wish to 
#    use a different interpolation scheme, which is determined when the initialize_new_paleobath instance 
#    is created with the interp_method argument.
#
# 2) This script assumes that there are no 'overhanging' ledges in the bathymetry (i.e. land cells above 
#    water cells). If there are--and I'm not even sure that would be allowed in CESM--this code will break
#    in all kinds of places. You should probably fix those overhangs.
#
# 3) So far, this script has been built to interpolate the restart data to NEW OCEAN CELLS, and is not 
#    equipped to remove ocean cells (turning them to land cells). This addition will likely come soon, it
#    was simply not necessary so far in this research project. It may not even be necessary, so long as 
#    KMT file of the new paleobathymetry is accurate. Not sure exactly how CESM would handle that.
#

# Published Oct. 20, 2017 - Andrew Vande Guchte

class CoordinateSystem(object):
    """The abstract coordinate system handling provided here works as follows.

    Each coordinate system must be able to convert data to a common coordinate system, which is chosen to be ECEF cartesian.
    data -> common system
    common system -> dislpay coordinates
    This is implemented by the fromECEF and toECEF methods in each coordinate system object.
    User code is responsible for taking data in its native coord system,
        transforming it using to/fromECEF using the a coord system appropriate to the data, and then
        transforming that data to the final coordinate system using another coord system.

    This class maintains an attribute WGS84xyz that can be used in
        transformations to/from the WGS84 ECEF cartesian system, e.g.
        >>> WGS84lla = proj4.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
        >>> projectedData = proj4.transform(WGS84lla, coordinateSystem.WGS84xyz, lat, lon, alt )
    The ECEF system has its origin at the center of the earth, with the +Z toward the north pole,
        +X toward (lat=0, lon=0), and +Y right-handed orthogonal to +X, +Z

    Depends on pyproj, http://code.google.com/p/pyproj/ to handle the ugly details of
    various map projections, geodetic transforms, etc.

    "You can think of a coordinate system as being something like character encodings,
    but messier, and without an obvious winner like UTF-8." - Django OSCON tutorial, 2007
    http://toys.jacobian.org/presentations/2007/oscon/tutorial/
    """

    WGS84xyz = proj4.Proj(proj='geocent',  ellps='WGS84', datum='WGS84')

    def coordinates():
        """Return a tuple of standarized coordinate names"""
        raise NotImplemented

    def fromECEF(self, x, y, z):
        """Take ECEF x, y, z values and return x, y, z in the coordinate system defined by the object subclass"""
        raise NotImplemented

    def toECEF(self, x, y, z):
        """Take x, y, z in the coordinate system defined by the object subclass and return ECEF x, y, z"""
        raise NotImplemented


class GeographicSystem(CoordinateSystem):
    """
    Coordinate system defined on the surface of the earth using latitude, longitide, and altitude, referenced to WGS84 ellipse
    """

    WGS84lla = proj4.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

    def toECEF(self, lon, lat, alt):
        projectedData = np.array(proj4.transform(GeographicSystem.WGS84lla, CoordinateSystem.WGS84xyz, lon, lat, alt ))
        if len(projectedData.shape) == 1:
            return projectedData[0], projectedData[1], projectedData[2]
        else:
            return projectedData[0,:], projectedData[1,:], projectedData[2,:]

    def fromECEF(self, x, y, z):
        projectedData = np.array(proj4.transform(CoordinateSystem.WGS84xyz, GeographicSystem.WGS84lla, x, y, z ))
        if len(projectedData.shape) == 1:
            return projectedData[0], projectedData[1], projectedData[2]
        else:
            return projectedData[0,:], projectedData[1,:], projectedData[2,:]



class initialize_new_paleobath(GeographicSystem):

    def __init__(self,oldtemp,newkmt,tlat,tlon,tdepth,interp_method='nhn'):
        '''The first argument is the variable from the restart file.
           The second argument is the shape of the new paleobath.
           The third argument
           The third argument is the interpolation method:
           '''

        self.tlat = tlat # Please make sure this is the T-grid, where scalar variables are defined.
        self.tlon = tlon # See above.
        self.tdepth = tdepth # See above.
        self.im = interp_method
        self.ns = newkmt # 2D new shape (new kmt data), not pythonic (ie unlike self.oldkmt).
        self.ov = oldtemp # 3D variable data from restart file (preferably temperature).
        self.oldkmt = self.getNearestPointAbove() # 2D kmt from oldvariable. Assumes no overhang.
                                                  # NOTE: This is kmt-1, so that it is a pythonic index.
        self.new = self.find_new_cells()

    def find_new_cells(self,):
        '''This function identifies which new ocean cells are being introduced with the new 
           kmt. It does not at this time identify where the ocean is being removed in favor
           of land.'''

        # Define a list of new cells to be filled for the new rst file.
        # Also define a list that will gather data for the interpolation scheme, if necessary.
        new = []
        interp_info = []

        print('\n\nFinding new seafloor...')

        for j in range(self.ov.shape[1]): # y-coord
            for i in range(self.ov.shape[2]): # x-coord
                for k in range(self.ov.shape[0]): # z-coord of new data

                    # Check to see if cursor is below old seafloor and above new seafloor...
                    if k > self.oldkmt[j,i] and k < self.ns[j,i]:

                    # ...and if it is, append the location to the 'new cell' array as [depth,lon,lat] indices.
                        new.append([k,j,i])
                    
                    # ...and if it is above old seafloor but below the new one (i.e. user added land cells)...
                    #elif k <= self.oldkmt[j,i] and k > self.ns[j,i]:
                    #    
                    #    pass # will be implemented when necessary

                    # ...and if it is above old seafloor, create an organized array of d,lon,lat,value for
                    # use in the interpolation methods of this code.
                    elif k <= self.oldkmt[j,i]:
                            if self.im == 'nn' or self.im == 'nhn':

                                interp_info.append([self.tdepth[k],self.tlon[j,i],self.tlat[j,i],k,j,i])

        # Convert new list into an array...size [x,3], where x is the number of cells to be filled,
        # and the other dimension is 0=k, 1=j, and 2=i. Also convert interp_info into an array.
        new = np.array(new)
        if self.im == 'nn' or self.im == 'nhn':
            self.interp_info = np.array(interp_info)

        return new

    def create_new_rst(self,rstin,rstout):
        '''Fill the empty cells in self.new with the determined values which can be created via any of
           the following methods:

           Nearest vertical neighbor (nvn)
           Nearest horizontal neighbor (nhn)
           Nearest neighbor (nn)

           Reads in rstin, writes out to rstout
           '''

        if self.im == 'nn': # Nearest-Neighbor (nn)

            # Convert lat, lon, depth grid to x, y, z on data grid:
            print('\n\nConverting from LLA to XYZ coordinate system...')
            self.xyz_grid = self.geodetic2geocentric(self.interp_info[:,1],self.interp_info[:,2],self.interp_info[:,0])

        # Read in old restart file and make a global attribute to be referenced in the 'writeout' function.
        self.f0 = nc.Dataset(rstin)
        keys = list(self.f0.variables.keys())

        # Initialize new data dictionary to be written to the output netcdf file.
        newrst = dict()

        for var in keys:#[-4:-3]: # Commented part is for testing.

            print(var)
            olddata = f0.variables['%s'%var][:]
            newdata = olddata

            # Test to see if a 3D variable...if it is, get the new values.
            if len(olddata.shape) == 3:

                newdata = self.fill_new_cells(newdata,olddata)

                #print(newdata[:,218,166]) # THIS LINE IS FOR TESTING THE VIABILITY OF FILLING.

            # Save the new data to the dictionary containing all the data to be written out as the new rst file:
            newrst[var] = newdata

        # Write new netcdf rst file with the new variables.
        self.writeout(rstout,newrst)

    def fill_new_cells(self,newdata,olddata):
        '''Gets and fills new cells with whichever method of interpolation was chosen.'''

        # Nearest Vertical Neighbor
        if self.im == 'nvn':

            for cell in self.new:

                newdata[cell[0],cell[1],cell[2]] = olddata[self.oldkmt[cell[1],cell[2]],cell[1],cell[2]]

            return newdata

        # Nearest Neighbor (based off of a Earth-centered, Earth-fixed (ECEF) coordinate system);
        # requires the lat, lon, and depth grids to convert.
        elif self.im == 'nn':

            # Create data grid that corresponds to the xyz_grid (in other words, only the defined (above-seafloor)
            # points that were in the olddata).
            print('\n\nIdentifying data from which to interpolate from in the old restart file...')
            defined_olddata = []
            for x in range(len(self.interp_info[:,3])):
                defined_olddata.append(olddata[self.interp_info[x,3],self.interp_info[x,4],self.interp_info[x,5]])
            defined_olddata = np.array(defined_olddata)

            # Create interpolation class instance:
            print('\n\nCreating nearest-neighbor interpolation matrix...')
            interpEngine = si.NearestNDInterpolator(self.xyz_grid,defined_olddata)

            for cell in self.new:

                xyz_cell = self.geodetic2geocentric(self.tlon[cell[1],cell[2]],self.tlat[cell[1],cell[2]],self.tdepth[cell[0]])
                newdata[cell[0],cell[1],cell[2]] = interpEngine.__call__(xyz_cell)

            return newdata

        # Nearest Horizontal Neighbor
        elif self.im == 'nhn':

            # Defined data and their locations (defined_olddata,xyz_grid) will be dictionaries, where the key is the
            # klvl (depth level). In other words, each level of data will be interpolated seperately, and based off of what
            # depth level they are at.
            xyz_grid = dict()
            defined_olddata = dict()

            print('\n\nFinding data and converting from LLA to XYZ coordinate systems on grid level: ')

            for klvl in range(self.ov.shape[0]):

                print('k = %d'%klvl)

                # Initialize the lists that will capture the data on each level.
                gridlla_at_this_level = []
                values_at_this_level = []

                for k,j,i in self.interp_info[:,3:6]:

                    if k == klvl:

                        # Append any data and locations in LLA that are on this level and are defined in the old rst file.
                        gridlla_at_this_level.append([self.tdepth[k],self.tlon[j,i],self.tlat[j,i]])
                        values_at_this_level.append(olddata[k,j,i])

                gridlla_at_this_level = np.array(gridlla_at_this_level)

                xyz_grid[klvl] = self.geodetic2geocentric(gridlla_at_this_level[:,1],gridlla_at_this_level[:,2],gridlla_at_this_level[:,0])
                defined_olddata[klvl] = np.array(values_at_this_level)

            print('\n\nBuilding interpolation matrix and interpolating to new grid cells on grid level: ')

            for klvl in range(self.ov.shape[0]):

                print('k = %d'%klvl)

                interpEngine = si.NearestNDInterpolator(xyz_grid[klvl],defined_olddata[klvl])

                for cell in self.new:

                    if cell[0] == klvl:

                        xyz_cell = self.geodetic2geocentric(self.tlon[cell[1],cell[2]],self.tlat[cell[1],cell[2]],self.tdepth[cell[0]])
                        newdata[cell[0],cell[1],cell[2]] = interpEngine.__call__(xyz_cell)

        return newdata

        # other interpolation methods will be built-in later...

    def geodetic2geocentric(self,llon,llat,ddepth):
        '''The latitude, longitude, and altitude of the gridded points are converted to earth-centered,
           earth-fixed cartesian grid (all longitude will be east longitude, as is convention for geodesy).
           Formulae from http://clynchg3c.com/Technote/geodesy/coordcvt.pdf.'''

        g = GeographicSystem()

        if llon.size > 1:

            x,y,z = g.toECEF(llon,llat,ddepth)
            xyz = [x,y,z]

            # Convert to [x,3] array with [x,y,z] in the second dimension.
            return np.array(xyz).T

        else:

            x,y,z = g.toECEF(llon,llat,ddepth)

            return np.array([x,y,z])


    def global_avg(self,tarea):
        '''Not implemented at this time; in the future this will likely be the global averege for each level'''

        pass

    def getNearestPointAbove(self,):
        '''This function creates an array that locates the last true value in a downwards-facing column search.
           This is a rough approximation, even less delicate than a nearest-neighbor approach. Effectively, it
           calculates the KMT value of the rst file.'''

        print('\n\nGetting the nearest vertical points...')

        last_k = np.empty(self.ns.shape)

        for j in range(self.ov.shape[1]): # y-coord
            for i in range(self.ov.shape[2]): # x-coord

                # Count the number of nonzero entries in the temperature data...subtract one to convert to an index.
                # This assumes there are no overhangs...which I believe is the case.
                last_k[j,i] = (np.count_nonzero(self.ov[:,j,i]))-1

        return last_k

    def writeout(self,outfilename,outdata):

        print('\n\nCreating output file...')

        # Determine the size of each dimension of data.
        numk = self.ov.shape[0]
        numj = self.ov.shape[1]
        numi = self.ov.shape[2]

        if os.path.isfile(outfilename):
            clobber = input('\n\nWARNING: %s already exists. Overwrite? 1=yes; 0=no. '%outfilename)

            if clobber:
                os.remove(outfilename)
            elif not clobber:
                outfilename = input('\n\nPlease input new output file name. To exit without saving, type "q". ')

            if outfilename == 'q':
                raise SystemExit('Exiting...')

        fout = nc.Dataset(outfilename,'w')
        fout.createDimension('k',numk)
        fout.createDimension('j',numj)
        fout.createDimension('i',numi)

        for v in outdata.keys():

            if len(outdata[v].shape) == 2:
                # MAY NEED TO DOUBLECHECK IF THIS IS THE CORRECT DATA TYPE...
                variable_instance = fout.createVariable(v,np.float64,('j','i',))
                variable_instance[:,:] = outdata[v]
            elif len(outdata[v].shape) == 3:
                variable_instance = fout.createVariable(v,np.float64,('k','j','i',))
                variable_instance[:,:,:] = outdata[v]
            else:
                raise SystemExit('Unsupported number of dimensions for variable %s on writeout.'%v)

        # Copy old restart file global attributes to new restart file:
        for name in self.f0.ncattrs():
            fout.setncattr(name,self.f0.getncattr(name))

        # Critical finalizing step to write attributes and variables to netcdf4 file:
        fout.close()

        print('\n\nFinished! New restart file: %s \n\n'%outfilename)

if __name__=='__main__':

    # Load test dataset and a test variable for the class to operate upon.
    filename0 = '../MAA_B1850C4CN_f19_g16_cret_4x_sewall.pop.r.0847-01-01-00000.nc'
    f0 = nc.Dataset(filename0)
    t = f0.variables['TEMP_CUR'][:]
    print(t.shape)

    filename1 = 'kmt.deeppanama1.ieeei4' #
    f1 = open(filename1,'r')
    kmt = np.fromfile(f1,dtype='>i4',count=-1,sep='')
    kmt = kmt.reshape((384,320))
    print(kmt.shape)

    # Technically not even needed...the code calculates the old kmt from the restart file temperature grid...
    #f = open('kmt.m.150519.ieeei4','r')
    #kmtc = np.fromfile(f,dtype='>i4',count=-1,sep='') # Big-endian (most significant byte first) 32-bit (4-byte) integers.
    #kmtc = kmtc.reshape((384,320))
    #print(kmtc.shape)

    f = nc.Dataset('../../data/pop2_clim/MAA_B1850C4CN_f19_g16_cret_4x_sewall.pop.h.851-850.nc')
    tlat = f.variables['TLAT'][:]
    tlon = f.variables['TLONG'][:]
    depth = f.variables['z_t'][:]/100.0*-1
    zzz = initialize_new_paleobath(t,kmt,tlat,tlon,depth)
    zzz.create_new_rst('../MAA_B1850C4CN_f19_g16_cret_4x_sewall.pop.r.0851-01-01-00000.nc','../testout_rst.nc')
