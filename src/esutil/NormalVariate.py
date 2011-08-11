from espresso import pmi

from _espresso import esutil_NormalVariate

class NormalVariateLocal(esutil_NormalVariate):
    def __init__(self, mean=0.0, sigma=1.0):
        cxxinit(self, esutil_NormalVariate, mean, sigma)

if pmi.isController:
    class NormalVariate(object):
        __metaclass__ = pmi.Proxy
        """A random normal variate."""
        pmiproxydefs = dict(
            cls = 'espresso.esutil.NormalVariateLocal',
            localcall = [ '__call__' ],
            )

