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

import _boostmpi
from _boostmpi import __doc__, __author__, __copyright__, __license__

class _ExtendBaseClass (type) :
    def __new__(self, name, bases, dict):
        del dict['__metaclass__']
        del dict['__module__']

        theClass = bases[0]
        # loop over all items in the class and replace it
        for k,v in dict.iteritems():
            setattr(theClass, k, v)
        return theClass

# extend the MPI Communicator class by the collective operations
class ExtendMPI(Communicator) :
    __metaclass__ = _ExtendBaseClass

    def broadcast(self, value=None, root=0) :
        """Calls broadcast(self, value, root).
        """
        return broadcast(self, value, root)

    def gather(self, value=None, root=0) :
        """Calls gather(self, value, root).
        """
        return gather(self, value, root)

    def reduce(self, value, op, root) :
        """Calls reduce(self, value, op, root).
        """
        return reduce(self, value, op, root)

    def all_gather(self, value) :
        """Calls all_gather(self, value).
        """
        return all_gather(self, value)

    def all_reduce(self, value, op) :
        """Calls all_reduce(self, value, op).
        """
        return all_reduce(self, value, op)

__author__ = __author__ + ', Olaf Lenz <lenzo@mpip-mainz.mpg.de>'
__copyright__ = __copyright__ + ', 2009 Olaf Lenz'

# export all symbols from _boostmpi
__all__ = [n for n in dir(_boostmpi) if not n.startswith("_")]
