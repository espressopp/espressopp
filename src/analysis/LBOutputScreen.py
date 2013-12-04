from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.LBOutput import *
from _espresso import analysis_LBOutput_Screen

class LBOutputScreenLocal(LBOutputLocal, analysis_LBOutput_Screen):
    'The (local) compute of LBOutputScreen.'
    def __init__(self, system, latticeboltzmann):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, analysis_LBOutput_Screen, system, latticeboltzmann)

if pmi.isController :
    class LBOutputScreen(LBOutput):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.analysis.LBOutputScreenLocal',
            pmicall = ["writeOutput"]
            )