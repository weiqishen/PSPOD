'''
ProbeReader reads data from the probe files.

Attributes:
    prefix           - folder name of the probe
    ppath            - absolute path of probe files
    surf_flag        - is surface or not
    num_probe        - number of probe files
    fields           - sampled fields
    dt               - sampling interval
    fnl_time         - final time
    n_snaps          - number of snapshots
    skip             - number of lines before data
    fields_id        - index of selected fields
    n_dim            - number of coordinate dimension
'''
import sys,os
import numpy as np
from mpi4py import MPI

class ProbeReader():

    def __init__(self, ppath, args):     
        self.rank = MPI.COMM_WORLD.Get_rank()
        self.ppath = ppath
        self.surf_flag=False
        self.setup()
        self.get_info(args)

    def setup(self):
        #  get abs dir and count number of probes
        if os.path.isdir(self.ppath):
            self.ppath = os.path.abspath(self.ppath)
        else:
            if self.rank==0:
                print("Input should be a path to probe files")
            sys.exit(1)
        self.prefix = os.path.basename(self.ppath)
        self.num_probe = len(os.listdir(self.ppath))

    def get_info(self, args):
        # open first file
        fname = self.ppath+'/{0:s}_{1:06d}.dat'.format(self.prefix, 0)
        pfile = open(fname, 'rt')

        self.skip = 0
        while True:  # get variable fields
            self.fields = pfile.readline().split()
            if not self.fields:  # if eof
                if self.rank==0:
                    print("Cant find fields")
                sys.exit(1)
            else:  # there's a line
                self.skip += 1
                if self.skip == 3: # coordinate line
                    self.n_dim = len(self.fields)
                if self.fields[0] == "time":  # if fields
                    break
                if self.skip == 4: # is a surface
                    self.surf_flag=True

        self.fields_id = []
        
        #initialize list of fields index
        if not args:# all fields is needed
            self.fields_id=[i for i in range(1,len(self.fields))]
            self.fields=self.fields[1:] # new fields without time
        else:
            for i in args:
                if i in self.fields:
                    self.fields_id.append(self.fields.index(i))
                else:
                    if self.rank==0:
                        print("{0:s} is not sampled".format(i))
                    sys.exit(1)
            self.fields = args # replace fields with args
        
        # close file
        pfile.close()

        # get sampling interval and number of samples, check interval uniformity
        time = np.loadtxt(fname, skiprows=self.skip,
                          usecols=(0), dtype=np.float)
        self.n_snaps = len(time)  # set number of snaps
        self.dt = np.round(time[1]-time[0], 10)
        self.fnl_time = time[-1]

    def load_probe(self, probe_id):
        fname = self.ppath+'/{0:s}_{1:06d}.dat'.format(self.prefix, probe_id)
        data = np.loadtxt(fname, skiprows=self.skip,
                          usecols=(self.fields_id), dtype=np.float)
        coord = np.loadtxt(fname, skiprows=2, max_rows=1, dtype=np.float)

        if self.surf_flag:
            normal = np.loadtxt(fname, skiprows=4, max_rows=1, dtype=np.float)
            area = np.loadtxt(fname, skiprows=6, max_rows=1,dtype=np.float)
        else:
            normal = None
            area = None
        # validate file
        if data.shape[0] != self.n_snaps:
            if self.rank == 0:
                print("Number of snapshots not equal to {0:d} in {1:s}".format(
                    self.n_snaps, fname))
                sys.exit(1)
        return data.T, coord, normal, area
