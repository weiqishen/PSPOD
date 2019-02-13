from probe2snap.probereader import ProbeReader
import h5py, sys
from err import Fatal_Error
from mpi4py import MPI

def p2s(ppath,args):

    # initialize MPI
    comm=MPI.COMM_WORLD
    rank = comm.Get_rank()
    nproc=comm.Get_size()
    # initialize probe reader object
    if rank==0:
        print("Reading probe info...")

    pr = ProbeReader(ppath, args)
    if rank == 0:
        print('Number of probes:{0:d}\nSampling interval: {1} sec\nSampling fields: {2}\nNumber of snapshots: {3:d}'.format(
            pr.num_probe, pr.dt, pr.fields, pr.n_snaps))

    # create h5file and dataset
    hfile=h5py.File(pr.prefix+'.h5','w',driver='mpio',comm=comm)
    hdata=hfile.create_dataset("snapshots",shape=(len(pr.fields_id),pr.num_probe,pr.n_snaps),dtype=float) #field*probe*snap

    # load one probe file at a time and fill into slice of dataset
    if rank==0:
        print("Converting probe file...")

    #set range of probe to read locally
    n_read = int(pr.num_probe/nproc)
    start_id = n_read*rank
    if rank == nproc-1:
        n_read = n_read+pr.num_probe-nproc*n_read
    
    #read corresponding file
    for i in range(start_id, start_id+n_read):
        if rank == 0:
            print("Progress:{0}%".format(
                round((i+1 - start_id) * 100 / n_read, 2)), end="\r")
        data = pr.load_probe(i)
        hdata[:, i, :] = data
    
    # all processors write attribute to the dataset
    hdata.attrs.create("n_probe", pr.num_probe, dtype="i4")
    hdata.attrs.create("dt", pr.dt, dtype=float)
    hdata.attrs.create("n_snaps", pr.n_snaps, dtype="i4")
    hdata.attrs.create("fields", pr.fields, dtype=h5py.special_dtype(vlen=str))
    hfile.close()

    if rank==0:
        print("\nDone.")
