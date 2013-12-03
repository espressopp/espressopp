"""
******************************************
**espresso.analysis.LBOutputProfileVzOfX**
******************************************

"""
from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.LBOutput import *
from _espresso import analysis_LBOutputProfile_VzOfX

class LBOutputProfileVzOfXLocal(LBOutputLocal, analysis_LBOutputProfile_VzOfX):
    'The (local) compute of LBOutputProfileVzOfX.'
    def __init__(self, system, latticeboltzmann):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, analysis_LBOutputProfile_VzOfX, system, latticeboltzmann)

if pmi.isController :
    class LBOutputProfileVzOfX(LBOutput):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.analysis.LBOutputProfileVzOfXLocal',
            pmicall = ["writeOutput"]
            )