import sys
from mpi4py import MPI


def Fatal_Error(err_str, err_code=1):
    if MPI.COMM_WORLD.Ger_rank == 0:
        print(err_str)
    sys.exit(err_code)
