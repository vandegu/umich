# read_cm1.py

class readcm1(object): ####################################################################################

    '''
       The binary output option in cm1 returns files in a very unusable and inconvenient
       output structure, unlike some nicer output options such as netCDF. This class is
       designed to easily and nicely unpack the cm1 binary datafiles into a usable 
       dictionary format, designed to emulate netCDF output.
       
       Currently, this class is equipped with methods to handle:

            Scalar output files (_s.dat, _w.dat)
            Parcel output files
            Statistics output files

            from cm1r18

       ^ While the vector output files are arranged in a very similar fashion to the 
       scalar output, there is not a reading method for them implemented at this time.

       Initial publishing: Oct 20,2015

       Update April 5, 2016: Implemented a method 'mesh' which produces a meshgrid of
                             any 2 dimensions from the 3 dimensional cm1r18 output,
                             as long as there exists a '_s.ctl' file.

       Update July 7, 2016: Implemented a method 'w_read' which reads '_w.dat' output 
                            in the same fashion as the '_s.dat' output files.
    '''

###########################################################################################################

    def __init__(self,filetag):
        
        '''
           Init only requires the filetag that is common to all of your model's output.
           eg: '/lustre/scratch/jsmith/testrun'
        '''
 
        self.filetag    = filetag 
        self.statsctl   = filetag+'_stats.ctl'
        self.statsfile  = filetag+'_stats.dat'
        self.scalarctl  = filetag+'_s.ctl'
        self.wctl       = filetag+'_w.ctl'
        # Scalar .dat files are defined in 'scalar read' below.
        self.pfile      = filetag+'_pdata.dat'
        self.pctl       = filetag+'_pdata.ctl'

###########################################################################################################

    def scalar_read(self,filenumber,):
        
        '''
           Reads the _s.dat files from cm1. Note that this is specifically for cm1r18
           or earlier. The only input necessary is the integer value of the scalar 
           data file. Output is a dictionary of all of the data stored in the scalar
           file.

           NOTE: This code requires the scalar data to be output at EACH TIMESTEP in 
           in a SEPERATE FILE (namelist: output_format = 1, output_filetype = 2).
        '''

        import numpy as np

        number_svars = 0
        number_twodvars = 0
        number_threedvars = 0
        skeys =[]
        start_var_key_list = 99999

        # Determine the size of the model domain.
        f = open(self.scalarctl,'r')
        for i,line in enumerate(f):
            linewords = line.split()
            if linewords[0] == 'xdef':
                nx = int(linewords[1])
            elif linewords[0] == 'ydef':
                ny = int(linewords[1])
            elif linewords[0] == 'zdef':
                nz = int(linewords[1])
            elif linewords[0] == 'vars':
                number_svars = int(linewords[1])
                start_var_key_list = i
            elif linewords[0] == 'tdef':
                number_stimes = int(linewords[1])
            if i > start_var_key_list and i < start_var_key_list+number_svars+1:
                skeys.append(linewords[0])
                if linewords[1] == '0':
                    number_twodvars = number_twodvars+1
                if linewords[0] == 'endvars': # Head off a common error here.
                    raise SystemExit('\n\nError reading scalar file....reached EOF.\n\n')
        f.close()        

        # Read the data file.

        sdata = []
        
        number_threedvars = number_svars-number_twodvars # This will differentiate between 2d and 3d vars.
        scalarfile = self.filetag+'_%06d_s.dat' %filenumber
        f = open(scalarfile,'r')

        # Read the 2D variables at the start of the file.
        twodvar_data = np.fromfile(f,dtype='float32',count=nx*ny*number_twodvars,sep='')
        twodvar_data = twodvar_data.reshape(number_twodvars,ny,nx)
        for i in range(number_twodvars):
            sdata.append(twodvar_data[i,:,:])
        
        # Read the 3D variables afterward.
        threedvar_data = np.fromfile(f,dtype='float32',count=nx*ny*nz*number_threedvars,sep='')
        threedvar_data = threedvar_data.reshape(number_threedvars,nz,ny,nx)
        for i in range(number_threedvars):
            sdata.append(threedvar_data[i,:,:,:])

        print('\n'+scalarfile+'-> 2D vars: '+str(number_twodvars)+' 3D vars: '+str(number_threedvars)+'\n')
         
        return dict(zip(skeys,sdata))

###########################################################################################################

    def w_read(self,filenumber,):
        
        '''
           Reads the _w.dat files from cm1. Note that this is specifically for cm1r18
           or earlier. The only input necessary is the integer value of the scalar 
           data file. Output is a dictionary of all of the data stored in the scalar
           file.

           NOTE: This code requires the scalar data to be output at EACH TIMESTEP in 
           in a SEPERATE FILE (namelist: output_format = 1, output_filetype = 2).
        '''

        import numpy as np

        number_svars = 0
        number_twodvars = 0
        number_threedvars = 0
        skeys =[]
        start_var_key_list = 99999

        # Determine the size of the model domain.
        f = open(self.wctl,'r')
        for i,line in enumerate(f):
            linewords = line.split()
            if linewords[0] == 'xdef':
                nx = int(linewords[1])
            elif linewords[0] == 'ydef':
                ny = int(linewords[1])
            elif linewords[0] == 'zdef':
                nz = int(linewords[1])
            elif linewords[0] == 'vars':
                number_svars = int(linewords[1])
                start_var_key_list = i
            elif linewords[0] == 'tdef':
                number_stimes = int(linewords[1])
            if i > start_var_key_list and i < start_var_key_list+number_svars+1:
                skeys.append(linewords[0])
                if linewords[1] == '0':
                    number_twodvars = number_twodvars+1
                if linewords[0] == 'endvars': # Head off a common error here.
                    raise SystemExit('\n\nError reading scalar file....reached EOF.\n\n')
        f.close()        

        # Read the data file.

        sdata = []
        
        number_threedvars = number_svars-number_twodvars # This will differentiate between 2d and 3d vars.
        wfile = self.filetag+'_%06d_w.dat' %filenumber
        f = open(wfile,'r')

        # Read the 2D variables at the start of the file.
        twodvar_data = np.fromfile(f,dtype='float32',count=nx*ny*number_twodvars,sep='')
        twodvar_data = twodvar_data.reshape(number_twodvars,ny,nx)
        for i in range(number_twodvars):
            sdata.append(twodvar_data[i,:,:])
        
        # Read the 3D variables afterward.
        threedvar_data = np.fromfile(f,dtype='float32',count=nx*ny*nz*number_threedvars,sep='')
        threedvar_data = threedvar_data.reshape(number_threedvars,nz,ny,nx)
        for i in range(number_threedvars):
            sdata.append(threedvar_data[i,:,:,:])

        print('\n'+wfile+'-> 2D vars: '+str(number_twodvars)+' 3D vars: '+str(number_threedvars)+'\n')
         
        return dict(zip(skeys,sdata))

###########################################################################################################

    def stats_read(self,):

        '''
           Read the binary statistics file output from cm1, and return a dictionary
           of the data contained therein, much like netCDF output. No arguments
           are necessary so long as an instance of readcm1 is extant.
        '''

        import os
        import numpy as np
        
        # Read the CTL file and determine keys, timesteps, and number of variables:

        keys = [] # Initialize the key list.
        f = open(self.statsctl,'r')
        for i,line in enumerate(f):
            linewords = line.split()
            if i == 7:
                statsvars = int(linewords[1])
            elif i > 7 and linewords[0] != 'endvars':
                keys.append(linewords[0])
        f.close()

        # NEW METHOD TO DETERMINE TIMES INCLUDED THAT ISN'T AS KLUDGY.
        # Determine the times included in a new method with the filesize.

        statstimes = os.path.getsize(self.statsfile)/(4*statsvars)
        print('\n'+self.statsfile+': '+str(statsvars)+' stats variables with '+str(statstimes)+' timesteps.\n')

        # Read the statistics outfile.        

        stats_data = [] # Initialize the data list.
        total = statsvars*statstimes
        read_stats_data = np.fromfile(self.statsfile,dtype='float32',count=total,sep='')
        #
        # A note how the stats are organized:
        #
        #           32 bit (4 byte) per value.
        #           Printed time sequentially, such that
        #           all the variables for a time are
        #           printed before the next time starts.
        #
        read_stats_data = read_stats_data.reshape(statstimes,statsvars)
        for i in range(statsvars):
            stats_data.append(read_stats_data[:,i])

        return dict(zip(keys,stats_data))

###########################################################################################################

    def parcel_read(self,parcel_desired,):

        '''
           Reads the timeseries data along a single parcel trajectory from a cm1 _pdata.dat
           file. Please be aware that this method reads the pdata CTL file to obtain info
           on the contents of the pdata.dat file, so if there are modifications to the cm1
           source code to output different variables such that they are not recorded in the
           pdata CTL file, you will need to modify the pdata CTL file accordingly in order 
           for this method to work properly.

           Input: integer callvalue of parcel.
           Output: dictionary of timeseries data along that parcel's trajectory.
        '''

        import os
        import numpy as np

        # Determine the number of parcels, timesteps, and variables written to the pdata file.

        keys =[]
        f = open(self.pctl,'r')
        for i,line in enumerate(f):
            linewords = line.split()
            if i == 3:
                number_parcels = int(linewords[1])-1 # on xdef line, for some reason.
            if i == 6:
                ptimes = int(linewords[1])
            elif i == 7:
                pvars = int(linewords[1])
            elif i > 7 and linewords[0] != 'endvars':
                keys.append(linewords[0])
        f.close()
        
        if parcel_desired >= number_parcels:
            raise SystemExit('\n\nError: Desired parcel must be between 0 and '+str(number_parcels-1)+'.\n\n')
        
        print('\n'+self.pfile+': Reading '+str(pvars)+' variables from parcel # '+str(parcel_desired)+'.\n')

        # Read the parcel data. NOTE: cm1r18 outputs a pdata file much larger than just 
        # the parcel data here, however, this still works because the actual parcel 
        # variable data is contained at the front of the pdata file.

        parcel_data = []
        pdata = np.zeros((ptimes,pvars))
        h = open(self.pfile,'r')
        for var in range(pvars):
            for time in range(ptimes):
                h.seek(time*number_parcels*pvars*4+number_parcels*var*4+parcel_desired*4)
                a = np.fromfile(h,dtype='float32',count=1,sep='')
                #print(a,time,var)
                pdata[time,var] = a
        for i in range(pvars):
            parcel_data.append(pdata[:,i])

        return dict(zip(keys,parcel_data))

###########################################################################################################

    def mesh(self,dimensions):

        '''
           Creates and returns a meshgrid for plotting any 2 dimensions of the 3 dimensional
           model output.

           Input: 2 dimensions (x,y, or z) with which to create the meshgrid (NOTE:
                  THE ARGUMENTS MUST BE IN HORIZONTAL, VERTICAL AXES ORDER, AND 
                  LOWER CASE, LIKE ['x','y']).
           Output: 2D meshgrid of input dimensions in the form X,Y (two outputs to unpack).
        '''

        import numpy as np

        # Assert that the proper type of argument was received by the 'mesh' method.
        assert isinstance(dimensions, list), '\n\nERROR: Wrong argument...readcm1.mesh requires a list.\n\n'

        # Determine the size of the model domain and define the model horizontal resolution.
        f = open(self.scalarctl,'r')
        for i,line in enumerate(f):
            linewords = line.split()
            if linewords[0] == 'xdef':
                nx = int(linewords[1])
                x0 = float(linewords[3])
                dx = float(linewords[4])
            elif linewords[0] == 'ydef':
                ny = int(linewords[1])
                y0 = float(linewords[3])
                dy = float(linewords[4])
            elif linewords[0] == 'zdef':
                nz = int(linewords[1])
                start_z_levels = i + 1
                end_z_levels = start_z_levels + nz

            mesh_dims = dict()
        
        # Create the meshgrid as specified by the 'dimensions' input argument.
        for i,dimension in enumerate(dimensions):
            if dimension == 'x':
                mesh_dims[i] = np.arange(x0,x0+nx*dx,dx)
            elif dimension == 'y':
                mesh_dims[i] = np.arange(y0,y0+ny*dy,dy)
            elif dimension == 'z':
                zs = []
                f.seek(0)
                for j,line in enumerate(f):
                    linewords = line.split()
                    if j >= start_z_levels and j < end_z_levels:
                        zs.append(float(linewords[0]))
                mesh_dims[i] = np.array(zs) # Not necessary to make array, but will for consistency's sake.
            else:
                raise SystemExit('\n\nERROR: readcm1.mesh requires dimensions x,y,or z in a 2 member list\n\n')
                
        return(np.meshgrid(mesh_dims[0],mesh_dims[1]))

###########################################################################################################

if __name__=='__main__': # For testing purposes:

    test = readcm1('/lustre/scratch/avandegu/extrapLBC/idd1')
    pdata = test.parcel_read(240000)
    x = pdata['x']
    y = pdata['y']
    z = pdata['z']
    print('X: %s, Y: %s, Z: %s' % (x[0],y[0],z[0]))
    

