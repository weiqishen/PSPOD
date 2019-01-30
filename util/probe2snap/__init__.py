from probe2snap.probereader import ProbeReader
import h5py, os
from err import Fatal_Error

def p2s(ppath,args):
    
    # initialize probe reader object
    print("Reading probe info...")
    pr=ProbeReader(ppath,args)

    # create h5file and dataset
    hfile=h5py.File(pr.prefix+'.h5','w')
    hdata=hfile.create_dataset("snapshots",shape=(len(pr.fields_id),pr.num_probe,pr.n_snaps),dtype=float) #field*probe*snap

    # load one probe file at a time and fill into slice of dataset
    print("Converting probe file...")
    for i in range(0,pr.num_probe):
        print("Progress:{0}%".format(round((i + 1) * 100 / pr.num_probe,2)),end="\r")
        data=pr.load_probe(i)
        hdata[:,i,:]=data
    # write attribute to the dataset
    hdata.attrs.create("num_probe", pr.num_probe, dtype="i4")
    hdata.attrs.create("dt", pr.dt, dtype=float)
    hdata.attrs.create("n_snaps", pr.n_snaps, dtype="i4")
    hdata.attrs.create("n_fields",len(pr.fields),dtype="i4")
    hdata.attrs.create("fields", pr.fields, dtype=h5py.special_dtype(vlen=str))
    hfile.close()
    
    print("\nDone.")
