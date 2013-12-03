"""
***************************
**espresso.analysis.NPart**
***************************

"""
from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.Observable import *
from _espresso import analysis_NPart

class NPartLocal(ObservableLocal, analysis_NPart):
    'The (local) compute of the number of particles of the system.'
    def __init__(self, system):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, analysis_NPart, system)

if pmi.isController :
    class NPart(Observable):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.analysis.NPartLocal'
        )
