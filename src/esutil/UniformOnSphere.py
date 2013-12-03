"""
***********************************
**espresso.esutil.UniformOnSphere**
***********************************

"""
from espresso import pmi

from _espresso import esutil_UniformOnSphere

class UniformOnSphereLocal(esutil_UniformOnSphere):
    pass

if pmi.isController:
    class UniformOnSphere(object):
        __metaclass__ = pmi.Proxy
        """A random variate that generates 3D vectors that are uniformly 
        distributed on a sphere."""
        pmiproxydefs = dict(
            cls = 'espresso.esutil.UniformOnSphereLocal',
            localcall = [ '__call__' ],
            )

