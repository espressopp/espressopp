"""
******************************
**espresso.analysis.LBOutput**
******************************

"""
from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.AnalysisBase import *
from _espresso import analysis_LBOutput

class LBOutputLocal(AnalysisBaseLocal, analysis_LBOutput):
    'The (local) compute of LBOutput.'
    def writeOutput(self):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.writeOutput(self)
        
if pmi.isController :
    class LBOutput(AnalysisBase):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
#            cls =  'espresso.analysis.LBOutputLocal',
#            pmicall = ["writeOutput"]
            )
        