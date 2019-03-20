''' 
This program serves as a format converter for the PSPOD main program.

Features:
    1. Convert ASCII probe files into hdf5 file format

Options:
    -h, --help                                      : print this help
    -p, --probe= <ascii_dir>                         : specify directory to the ascii probe files and convert to hdf5 

For more information goto https://github.com/weiqishen/HiFiLES-solver/wiki
'''
import sys
import p2h5
import getopt
from mpi4py import MPI

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    # parse command line options
    try:
        opts, args = getopt.getopt(argv, "hp:", ["help", "probe="])
    except:
        if MPI.COMM_WORLD.Get_rank()==0:
            print("for help use --help")
        return 0

    # read command line options
    if not opts:
        if MPI.COMM_WORLD.Get_rank()==0:
            print("Missing options. For help use -h, --help")
        return 1
    for opt, value in opts:
        if opt in ("-h", "--help"):  # load help document
            if MPI.COMM_WORLD.Get_rank()==0:
                print(__doc__)
            return 0
        elif opt in ("-p", "--probe"): # work mode: probe2snap
            return p2h5.convert(value,args)
    else:
        if MPI.COMM_WORLD.Get_rank()==0:
            print("Options not recognized. For help use h or --help")
        return 1

# execute main funtion
if __name__ == "__main__":
    sys.exit(main())
