# MPI wrapper module
# Use this to load the real mpi module, which might live under different names
try :
#    print("Importing boost.mpi module...")
    from boost.mpi import *
except ImportError :
#    print("Caught exception, trying to load mpi module...")
    from mpi import *
