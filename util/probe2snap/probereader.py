'''
ProbeReader reads data from the probe files.

Attributes:
    prefix           - folder name of the probe
    ppath            - absolute path of probe files
    num_probe        - number of probe files
    fields           - sampled fields
    dt               - sampling interval
    n_snaps          - number of snapshots
    skip             - number of lines before data
    fields_id         - index of selected fields
'''
import os
from err import Fatal_Error
import numpy as np


class ProbeReader():

    def __init__(self, ppath, args):
        if not args:
            Fatal_Error("Please specify fields to extract")
        self.ppath = ppath
        self.n_snaps = 0

        self.setup()
        self.get_info(args)
        print('Number of probes:{0:d}\nSampling interval: {1} sec\nSampling fields: {2}\nNumber of snapshots: {3:d}'.format(
            self.num_probe, self.dt, self.fields, self.n_snaps))

    def setup(self):
        #  get abs dir and count number of probes
        if os.path.isdir(self.ppath):
            self.ppath = os.path.abspath(self.ppath)
        else:
            Fatal_Error("Input should be a path to probe files")
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
                Fatal_Error("Can't find fields")
            else:  # there's a line
                self.skip += 1
                if self.fields[0] == "time":  # if fields
                    break

        self.fields_id = []
        
        #initialize list of fields index
        for i in args:
            if i in self.fields:
                self.fields_id.append(self.fields.index(i))
            else:
                Fatal_Error("{0:s} is not sampled".format(i))

        self.fields = args  # replace fields with args
        pfile.close()

        # get sampling interval and number of samples, check interval uniformity
        time = np.loadtxt(fname, skiprows=self.skip,
                          usecols=(0), dtype=np.float)
        self.n_snaps = len(time)  # set number of snaps

        for i in range(0, self.n_snaps):
            if i == 0:  # first snapshot
                t1 = time[i]
            elif i == 1:  # second snapshot
                # round to 10 digit after dot
                self.dt = np.round(time[i]-t1, 10)
                t1 = time[i]
            else:
                if np.round(time[i]-t1, 10) != self.dt:
                    Fatal_Error("non-uniform sampling interval")
                else:
                    t1 = time[i]

    def load_probe(self, probe_id):
        fname = self.ppath+'/{0:s}_{1:06d}.dat'.format(self.prefix, probe_id)
        data = np.loadtxt(fname, skiprows=self.skip, usecols=(self.fields_id), dtype=np.float)
        if data.shape[0] != self.n_snaps:
            Fatal_Error("Number of snapshots not equal to {0:d} in {1:s}".format(self.n_snaps,fname))
        return data.T
