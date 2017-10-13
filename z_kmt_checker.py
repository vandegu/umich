import os
import numpy as np
import netCDF4 as nc

class avoid_mask_edits(object):

    def __init__(self,controlkmt,testkmt,show=1,netcdf=False):
        '''This module is intended to replace disperate points in the testkmt
           file from the controlkmt file. Specifically, it will isolate where
           there are differences in the land/sea distribution between the two
           files. The user will be asked to edit these disperate points.
        '''

        #import os
        #import numpy as np
        #import netCDF4 as nc
        self.netcdf = netcdf
        file1 = controlkmt
        file2 = testkmt

        if self.netcdf:
            self.fc = nc.Dataset(file1,'r')
            self.f1 = nc.Dataset(file2,'r')
            self.kmtc = self.fc.variables['kmt'][:]
            self.kmt1 = self.f1.variables['kmt'][:]

        else:
            infile1 = open(file1)
            # Read binary big-endian 32-bit integer array of KMT values:
            kmtc = np.fromfile(infile1,dtype='>i4',count=-1,sep='')
            self.kmtc = kmtc.reshape((384,320))
            infile2 = open(file2)
            kmt1 = np.fromfile(infile2,dtype='>i4',count=-1,sep='')
            self.kmt1= kmt1.reshape((384,320))

        self.maskdiff,self.count = self.mask_check(self.kmtc,self.kmt1)
        self.show = show

        #if show:
        #    print(self.fc.variables)

    def mask_check(self,c,t):
        self.maskc = np.copy(c)
        self.mask1 = np.copy(t)
        seac = np.invert(c == 0.)
        sea1 = np.invert(t == 0.)
        self.maskc[seac] = 1.
        self.mask1[sea1] = 1.
        maskdiff = self.maskc-self.mask1
        count = np.count_nonzero(maskdiff)

        return maskdiff,count

    def replace(self,newname):

        points_to_change = np.array(np.nonzero(self.maskdiff)).T

        for pointx in points_to_change:
            self.edit_point(pointx)

        self.writeout(newname)


    def edit_point(self,point):

        if self.show:
            print(point)

        print('''Control:
                 %4.f%4.f%4.f%4.f%4.f
                 %4.f%4.f%4.f%4.f%4.f
                 %4.f%4.f%4.f*%3.f%4.f
                 %4.f%4.f%4.f%4.f%4.f
                 %4.f%4.f%4.f%4.f%4.f\n'''%(
                self.kmtc[point[0]+2,point[1]-2],
                self.kmtc[point[0]+2,point[1]-1],
                self.kmtc[point[0]+2,point[1]],
                self.kmtc[point[0]+2,point[1]+1],
                self.kmtc[point[0]+2,point[1]+2],
                self.kmtc[point[0]+1,point[1]-2],
                self.kmtc[point[0]+1,point[1]-1],
                self.kmtc[point[0]+1,point[1]],
                self.kmtc[point[0]+1,point[1]+1],
                self.kmtc[point[0]+1,point[1]+2],
                self.kmtc[point[0],point[1]-2],
                self.kmtc[point[0],point[1]-1],
                self.kmtc[point[0],point[1]],
                self.kmtc[point[0],point[1]+1],
                self.kmtc[point[0],point[1]+2],
                self.kmtc[point[0]-1,point[1]-2],
                self.kmtc[point[0]-1,point[1]-1],
                self.kmtc[point[0]-1,point[1]],
                self.kmtc[point[0]-1,point[1]+1],
                self.kmtc[point[0]-1,point[1]+2],
                self.kmtc[point[0]-2,point[1]-2],
                self.kmtc[point[0]-2,point[1]-1],
                self.kmtc[point[0]-2,point[1]],
                self.kmtc[point[0]-2,point[1]+1],
                self.kmtc[point[0]-2,point[1]+2]))
        print('''Test:
                 %4.f%4.f%4.f%4.f%4.f
                 %4.f%4.f%4.f%4.f%4.f
                 %4.f%4.f%4.f*%3.f%4.f
                 %4.f%4.f%4.f%4.f%4.f
                 %4.f%4.f%4.f%4.f%4.f\n'''%(
                self.kmt1[point[0]+2,point[1]-2],
                self.kmt1[point[0]+2,point[1]-1],
                self.kmt1[point[0]+2,point[1]],
                self.kmt1[point[0]+2,point[1]+1],
                self.kmt1[point[0]+2,point[1]+2],
                self.kmt1[point[0]+1,point[1]-2],
                self.kmt1[point[0]+1,point[1]-1],
                self.kmt1[point[0]+1,point[1]],
                self.kmt1[point[0]+1,point[1]+1],
                self.kmt1[point[0]+1,point[1]+2],
                self.kmt1[point[0],point[1]-2],
                self.kmt1[point[0],point[1]-1],
                self.kmt1[point[0],point[1]],
                self.kmt1[point[0],point[1]+1],
                self.kmt1[point[0],point[1]+2],
                self.kmt1[point[0]-1,point[1]-2],
                self.kmt1[point[0]-1,point[1]-1],
                self.kmt1[point[0]-1,point[1]],
                self.kmt1[point[0]-1,point[1]+1],
                self.kmt1[point[0]-1,point[1]+2],
                self.kmt1[point[0]-2,point[1]-2],
                self.kmt1[point[0]-2,point[1]-1],
                self.kmt1[point[0]-2,point[1]],
                self.kmt1[point[0]-2,point[1]+1],
                self.kmt1[point[0]-2,point[1]+2]))
        change = input('\nWhat should the middle be changed into?')
        self.kmt1[point[0],point[1]] = change
        print('''Result:
                 %4.f%4.f%4.f%4.f%4.f
                 %4.f%4.f%4.f%4.f%4.f
                 %4.f%4.f%4.f*%3.f%4.f
                 %4.f%4.f%4.f%4.f%4.f
                 %4.f%4.f%4.f%4.f%4.f\n'''%(
                self.kmt1[point[0]+2,point[1]-2],
                self.kmt1[point[0]+2,point[1]-1],
                self.kmt1[point[0]+2,point[1]],
                self.kmt1[point[0]+2,point[1]+1],
                self.kmt1[point[0]+2,point[1]+2],
                self.kmt1[point[0]+1,point[1]-2],
                self.kmt1[point[0]+1,point[1]-1],
                self.kmt1[point[0]+1,point[1]],
                self.kmt1[point[0]+1,point[1]+1],
                self.kmt1[point[0]+1,point[1]+2],
                self.kmt1[point[0],point[1]-2],
                self.kmt1[point[0],point[1]-1],
                self.kmt1[point[0],point[1]],
                self.kmt1[point[0],point[1]+1],
                self.kmt1[point[0],point[1]+2],
                self.kmt1[point[0]-1,point[1]-2],
                self.kmt1[point[0]-1,point[1]-1],
                self.kmt1[point[0]-1,point[1]],
                self.kmt1[point[0]-1,point[1]+1],
                self.kmt1[point[0]-1,point[1]+2],
                self.kmt1[point[0]-2,point[1]-2],
                self.kmt1[point[0]-2,point[1]-1],
                self.kmt1[point[0]-2,point[1]],
                self.kmt1[point[0]-2,point[1]+1],
                self.kmt1[point[0]-2,point[1]+2]))

    def writeout(self,newfilename):

        if os.path.isfile(newfilename):
            overwrite = input('\n* * * The file %s already exists! Overwrite? (1=yes,0=no)\n'%newfilename)
            if overwrite:
                os.remove(newfilename)
            else:
                newfilename = input("\nEnter a new file name. Be sure it isn't taken, or an error will be thrown.\n")

        if self.netcdf:

            orig_variables = self.f1.variables.keys()
            f3 = nc.Dataset(newfilename,'w',format='NETCDF4')
            f3.createDimension('latitude',self.kmt1.shape[0])
            f3.createDimension('longitude',self.kmt1.shape[1])

            for ov in orig_variables:
                data = self.fc.variables[ov][:]
                if ov == 'kmt' or ov == 'KMT':
                    variable_instance = f3.createVariable('kmt',np.float32,('latitude','longitude',))
                    variable_instance[:,:] = self.kmt1
                else:
                    variable_instance = f3.createVariable(ov,np.float32,('latitude','longitude',))
                    variable_instance[:,:] = data

            print('\n\nFinished! New data saved as "kmt" variable in %s.\n\n'%newfilename)

        else:

            # Write out to a binary big-endian 32-bit integer file.
            self.kmt1.astype('>i4').tofile(newfilename,sep='')

if __name__ == "__main__":

    filec = input('\nControl kmt file (the one that you wish to emulate its kmt mask structure)\n')
    if not os.path.isfile(filec):
        raise SystemExit('\n\n* * * CONTROL FILE DOES NOT EXIST * * *\n\n')
    file1 = input('\nThe kmt file that you wish to edit:\n')
    if not os.path.isfile(file1):
        raise SystemExit('\n\n* * * TEST FILE DOES NOT EXIST * * *\n\n')
    file2 = input('''\nOutput kmt file. If you wish to overwrite, simply place the same filename.
Later, you will be asked to confirm that you want to overwrite this file.\n''')

    ins = avoid_mask_edits(filec,file1)
    ins.replace(file2)
