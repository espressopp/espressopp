from espresso import pmi
from _espresso import analysis_AllParticlePos

class AllParticlePosLocal(object):
    pass

if pmi.isController :
    class AllParticlePos(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict()
