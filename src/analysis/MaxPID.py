"""
****************************
**espresso.analysis.MaxPID**
****************************

"""
from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.Observable import *
from _espresso import analysis_MaxPID

class MaxPIDLocal(ObservableLocal, analysis_MaxPID):
    'The (local) compute of the maximum pid number of the system.'
    def __init__(self, system):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, analysis_MaxPID, system)

if pmi.isController :
    class MaxPID(Observable):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.analysis.MaxPIDLocal'
            )
