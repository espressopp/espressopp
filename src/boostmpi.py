# Adapted from boostmpi (src/python/__init__.py)
import sys

if sys.platform == 'linux2':
    import DLFCN as dl
    flags = sys.getdlopenflags()
    sys.setdlopenflags(dl.RTLD_NOW|dl.RTLD_GLOBAL)
    from _boostmpi import *
    sys.setdlopenflags(flags)
else:
    from _boostmpi import *

import esutil

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
        return reduce(self, value, op, root)
