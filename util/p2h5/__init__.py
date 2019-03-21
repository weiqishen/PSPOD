from p2h5.probereader import ProbeReader
import h5py, sys
import numpy as np
from mpi4py import MPI

def convert(ppath,args):

    # initialize MPI
    rank = MPI.COMM_WORLD.Get_rank()
    nproc=MPI.COMM_WORLD.Get_size()
    # initialize probe reader object
    if rank==0:
        print("Reading probe info...")

    pr = ProbeReader(ppath,args)

    if rank == 0:
        print( 'Number of probes: {0:d}\n\
                Is surface: {1}\n\
                Sampling interval: {2} sec\n\
                Sampling fields: {3}\n\
                Simulation dimension: {4:d}\n\
                Number of snapshots: {5:d}'.format(
                pr.num_probe, pr.surf_flag, pr.dt, pr.fields, pr.n_dim, pr.n_snaps))

    # create h5file and dataset
    hfile=h5py.File(pr.prefix+'.h5','w',driver='mpio',comm=MPI.COMM_WORLD)
    hdata=hfile.create_dataset("data",shape=(len(pr.fields),pr.num_probe,pr.n_snaps),
                                maxshape=(len(pr.fields),pr.num_probe,None),chunks=(len(pr.fields),pr.num_probe,1)
                                ,dtype=float) #field*probe*snap
    hcoord = hfile.create_dataset("coord", shape=(pr.num_probe, pr.n_dim), dtype=float)  # probe*dim
    if pr.surf_flag:
        hnormal=hfile.create_dataset("normal", shape=(pr.num_probe, pr.n_dim), dtype=float)  # probe*dim
        harea=hfile.create_dataset("area", shape=(pr.num_probe,), dtype=float)  # probe

    # load one probe file at a time and fill into slice of dataset
    if rank==0:
        print("loading probe file...")

    #set range of probe to read locally
    n_read = int(pr.num_probe/nproc)
    start_id = n_read*rank
    if rank == nproc-1:
        n_read += pr.num_probe-nproc*n_read
    
    #declare array to hold local data
    data = np.empty(shape=(len(pr.fields), n_read, pr.n_snaps), dtype=np.float)

    #read and write data
    for i in range(start_id, start_id+n_read):
        if rank == 0:
            print("Progress:{0}%".format(
                round((i+1 - start_id) * 100 / n_read, 2)), end="\r")
        data[:, i-start_id, :], coord, normal, area = pr.load_probe(i)
        hcoord[i, :] = coord
        if pr.surf_flag:
            hnormal[i,:]=normal
            harea[i]=area

    if rank == 0:
        print("\nWritinng to hdf5 file")
    for i in range(0, pr.n_snaps):
        if rank == 0:
            print("Progress:{0}%".format(
                round((i+1) * 100 / pr.n_snaps, 2)), end="\r")
        hdata[:, start_id:start_id+n_read, i] = data[:, :, i]
    # all processors write attribute to the dataset
    hfile.attrs["dt"]=pr.dt
    hfile.attrs.create("fields", pr.fields, dtype=h5py.special_dtype(vlen=str))
    hfile.attrs["fnl_time"]=pr.fnl_time
    hfile.close()

    if rank==0:
        print("\nDone.")
