from espresso import pmi

from _espresso import esutil_RNG

class RNGLocal(esutil_RNG):
    pass

if pmi.isController:
    class RNG(object):
        __metaclass__ = pmi.Proxy
        'Random number generator.'
        pmiproxydefs = dict(
            cls = 'espresso.esutil.RNGLocal',
            localcall = [ '__call__', 'normal', 'uniformOnSphere' ],
            pmicall = [ 'seed' ]
            )
    
