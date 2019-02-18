'''
This file convert output modes to vtk file for SPOD visualization
'''

import sys
import numpy as np
from pyevtk.hl import gridToVTK 
import h5py

# params
input_file = 'pod_done.h5'
freq_id = 3  # index of frequency to plot
mode_id = 0  # index of mode to plot
n_xyz = (121, 51, 51)  # number of points in each dir
fields = ['pressure']  # fields to plot
D_n=0.0508 #nozzle diameter

# open output hdf5 file
houtput = h5py.File(input_file, 'r')
#read attribution
in_field = houtput.attrs["fields"].tolist()
for i in range(0,len(in_field)):
    in_field[i]=in_field[i].decode('ascii')

df = houtput.attrs["df"]
n_probe = houtput.attrs["n_probe"]

#find out field id
field_id = list()
for i in fields:
    if i in in_field:
        field_id.append(in_field.index(i))  # save index to field_id
    else:
        print("field not found")
        sys.exit(1)

# read coordinate, store seprately
coord = houtput["coord"][...]  # probe*dimension
coord=coord/D_n
# load corresponding field
data = {}

for i in range(0, len(fields)):
    data[fields[i]+"_real"] = houtput["modes_real"][freq_id, mode_id,
                                            field_id[i]*n_probe:(field_id[i]+1)*n_probe].reshape(n_xyz,order="F")
    data[fields[i]+"_imag"] = houtput["modes_imag"][freq_id, mode_id,
                                            field_id[i]*n_probe:(field_id[i]+1)*n_probe].reshape(n_xyz,order="F")
    data[fields[i]+"_mag"]=np.abs(data[fields[i]+"_real"]+1j*data[fields[i]+"_imag"])                                     

x = coord[:, 0].reshape(n_xyz,order="F").copy()
y = coord[:, 1].reshape(n_xyz,order="F").copy()
z = coord[:, 2].reshape(n_xyz,order="F").copy()

gridToVTK("./st_{0:04d}_{1:04d}".format(freq_id,mode_id), x, y, z, pointData=data)

houtput.close()
