# MPI wrapper module
# Use this to load the real mpi module, which might live under different names
import logging
import esutil

log=logging.getLogger(__name__)

# Try to load mpi module
try :
    from boost.mpi import *
except ImportError :
    log.debug("Caught import exception when importing boost.mpi, loading mpi module...")
    from mpi import *

# extend the MPI Communicator class by the collective operations
class ExtendMPI(Communicator) :
    __metaclass__ = esutil.ExtendBaseClass

    def broadcast(self, value=None, root=0) :
        """Calls __module__.broadcast(self, value, root)."""
        return broadcast(self, value, root)

    def gather(self, value=None, root=0) :
        """Calls __module__.gather(self, value, root)."""
        return gather(self, value, root)

    def reduce(self, value, op, root) :
        """Calls __module__.reduce(self, value, op, root)."""
        return reduce(self, *args)
