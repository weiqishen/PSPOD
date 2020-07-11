'''
This file convert output modes to vtk file for SPOD visualization
'''

import sys,os
import numpy as np
from pyevtk.hl import gridToVTK 
import h5py

# params
def plot_vtk(freq_id,mode_id,input_file):
    n_xyz = (130, 61, 61)  # number of points in each dir
    D_n=0.0508 #nozzle diameter
    Uj = 762.389735 #jet velocity
    # open output hdf5 file
    houtput = h5py.File(input_file, 'r')
    #read attribution
    in_field = houtput.attrs["fields"].tolist()
    for i in range(0,len(in_field)):
        in_field[i]=in_field[i].decode('ascii') # bytes to string

    df = houtput.attrs["df"]
    print('St={}'.format(df*freq_id*D_n/Uj))
    n_probe = int(houtput["modes_real"].shape[2]/len(in_field))

    # read coordinate, store seprately
    coord = houtput["coord"][...]  # probe*dimension
    coord=coord/D_n
    # load corresponding field
    data = {}

    for i in range(0, len(in_field)):
        data[in_field[i]+"_real"] = houtput["modes_real"][freq_id, mode_id,
                                               i*n_probe:(i+1)*n_probe].reshape(n_xyz,order="F") #fortran layout
        data[in_field[i]+"_imag"] = houtput["modes_imag"][freq_id, mode_id,
                                                i*n_probe:(i+1)*n_probe].reshape(n_xyz,order="F")
        data[in_field[i]+"_mag"]=np.abs(data[in_field[i]+"_real"]+1j*data[in_field[i]+"_imag"])                                     

    x = coord[:, 0].copy().reshape(n_xyz,order="F")#fortran layout
    y = coord[:, 1].copy().reshape(n_xyz,order="F")#fortran layout
    z = coord[:, 2].copy().reshape(n_xyz,order="F")#fortran layout

    gridToVTK("./st_{0:04d}_{1:04d}".format(freq_id,mode_id), x, y, z, pointData=data)

    houtput.close()

for i in range(0,5):
    input_fname='pod_result.h5'
    plot_vtk(5,i,input_fname)