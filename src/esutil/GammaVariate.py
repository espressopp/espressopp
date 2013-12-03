"""
********************************
**espresso.esutil.GammaVariate**
********************************

"""
from espresso import pmi

from _espresso import esutil_GammaVariate

class GammaVariateLocal(esutil_GammaVariate):
    def __init__(self, alpha, beta):
        cxxinit(self, esutil_GammaVariate, alpha, beta)

if pmi.isController:
    class GammaVariate(object):
        __metaclass__ = pmi.Proxy
        """A random gamma variate."""
        pmiproxydefs = dict(
            cls = 'espresso.esutil.GammaVariateLocal',
            localcall = [ '__call__' ],
            )

