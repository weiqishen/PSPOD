''' 
This program serves as a utility for the PSPOD main program.

Features:
    1. Convert HiFiLES surface/volume probe data into PSPOD snapshot file
    2. Convert PSPOD snapshot file or PSPOD output files into visualization data

Options:
    -h, --help                    : print this help
    -p, --probe <probe_dir>       : specify directory to the probe files and convert to snapshots 

For more information goto https://github.com/weiqishen/HiFiLES-solver/wiki
'''
import sys
from probe2snap import p2s
from err import Fatal_Error
import getopt
from mpi4py import MPI

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    # parse command line options
    try:
        opts, args = getopt.getopt(argv, "hp:", ["help", "probe="])
    except:
        Fatal_Error("for help use --help")

    # read command line options
    if not opts:
        Fatal_Error("Missing options. For help use -h, --help")
    for opt, value in opts:
        if opt in ("-h", "--help"):  # load help document
            if MPI.COMM_WORLD.Get_rank()==0:
                print(__doc__)
            return 0
        elif opt in ("-p", "--probe"): # work mode: probe2snap
            p2s(value,args)
            return 0
    else:
        Fatal_Error("for help use --help")

# execute main funtion
if __name__ == "__main__":
    sys.exit(main())
