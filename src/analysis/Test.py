"""
**************************
**espresso.analysis.Test**
**************************

"""
from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.AnalysisBase import *
from _espresso import analysis_Test

class TestLocal(AnalysisBaseLocal, analysis_Test):
    'The (local) test of analysis.'
    def __init__(self, system):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, analysis_Test, system)

if pmi.isController :
    class Test(AnalysisBase):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.analysis.TestLocal'
            )
