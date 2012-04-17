from espresso import pmi

from _espresso import esutil_RNG

class RNGLocal(esutil_RNG):
  pass

#    def gamma(self, a=None):
#          if pmi._PMIComm and pmi._PMIComm.isActive():
#            if pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
#                if a==None:
#                    return self.cxxclass.gammaArg(self, 1)
#                else:
#                    return self.cxxclass.gammaArg(self, a)
#            else :
#                pass
    
if pmi.isController:
    class RNG(object):
        __metaclass__ = pmi.Proxy
        'Random number generator.'
        pmiproxydefs = dict(
            cls = 'espresso.esutil.RNGLocal',
            localcall = [ '__call__', 'normal', 'gamma', 'uniformOnSphere' ],
            pmicall = [ 'seed' ]
            )
    
