from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.LBOutput import *
from _espresso import analysis_LBOutput_VzInTime

class LBOutputVzInTimeLocal(LBOutputLocal, analysis_LBOutput_VzInTime):
    'The (local) compute of LBOutputVzInTime.'
    def __init__(self, system, latticeboltzmann):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, analysis_LBOutput_VzInTime, system, latticeboltzmann)

if pmi.isController :
    class LBOutputVzInTime(LBOutput):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.analysis.LBOutputVzInTimeLocal',
            pmicall = ["writeOutput"]
            )