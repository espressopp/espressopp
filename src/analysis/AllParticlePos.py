from espresso import pmi
from _espresso import analysis_AllParticlePos

class AllParticlePosLocal(object):
    """Abstract local base class for observables."""
    def gatherAllPositions(self):
      return self.cxxclass.gatherAllPositions(self)

if pmi.isController :
    class AllParticlePos(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            pmicall = [ "gatherAllPositions" ]
        )
