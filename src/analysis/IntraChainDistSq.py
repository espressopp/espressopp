from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.AllParticlePos import *
from _espresso import analysis_IntraChainDistSq

class IntraChainDistSqLocal(AllParticlePosLocal, analysis_IntraChainDistSq):
    'The (local) IntraChainDistSq object'
    def __init__(self, system, fpl):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
          cxxinit(self, analysis_IntraChainDistSq, system, fpl)
    def compute(self):
        return self.cxxclass.compute(self)

if pmi.isController :
    class IntraChainDistSq(AllParticlePos):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.analysis.IntraChainDistSqLocal',
            pmicall = [ "compute" ]
            )
