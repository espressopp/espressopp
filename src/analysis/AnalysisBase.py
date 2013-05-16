from espresso import pmi
from _espresso import analysis_AnalysisBase

class AnalysisBaseLocal(analysis_AnalysisBase):
    """Abstract local base class for observables."""
    def __init__(self, system):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, analysis_AnalysisBase, system)

    def compute(self):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.compute(self)

    def reset(self):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.reset(self)

    def getInstantValue(self):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            res = self.cxxclass.getInstantValue(self)
            if len(res) > 1:
                return res
            else:
                return res[0]

    def getAverageValue(self):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getAverageValue(self)

    def getNumberOfMeasurements(self):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getNumberOfMeasurements(self)

if pmi.isController :
    class AnalysisBase(object):
        """Abstract base class for observable."""
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            pmicall = [ "compute", "reset", "getInstantValue", "getAverageValue", "getNumberOfMeasurements" ]
            )
