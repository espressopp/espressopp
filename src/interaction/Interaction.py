from espresso import pmi
from _espresso import interaction_Interaction

class InteractionLocal(object):
    """Abstract local base class for interactions."""
    def computeEnergy(self):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.computeEnergy(self)

if pmi.isController :
    class Interaction(object):
        """Abstract base class for interaction."""
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            pmicall = [ "computeEnergy" ]
            )
