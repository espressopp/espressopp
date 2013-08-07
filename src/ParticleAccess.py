"""
***************************************************************
**ParticleAccess** - abstract base class for analysis/measurement/io
***************************************************************

"""

from espresso import pmi
from _espresso import ParticleAccess

class ParticleAccessLocal(ParticleAccess):
    """Abstract local base class"""
    def perform_action(self):
      if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        self.cxxclass.perform_action(self)

if pmi.isController :
    class ParticleAccess(object):
        """Abstract base class"""
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            pmicall = [ 'perform_action' ]
        )
