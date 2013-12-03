"""
************************************************
**espresso.analysis.ParticleRadiusDistribution**
************************************************

"""
from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.AnalysisBase import *
from _espresso import analysis_ParticleRadiusDistribution

class ParticleRadiusDistributionLocal(AnalysisBaseLocal, analysis_ParticleRadiusDistribution):
    'The (local) compute of the particle radius distribution.'
    def __init__(self, system):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, analysis_ParticleRadiusDistribution, system)

if pmi.isController :
    class ParticleRadiusDistribution(AnalysisBase):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.analysis.ParticleRadiusDistributionLocal'
            )
