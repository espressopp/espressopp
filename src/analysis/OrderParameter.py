from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.AnalysisBase import *
from _espresso import analysis_OrderParameter

class OrderParameterLocal(AnalysisBaseLocal, analysis_OrderParameter):
    'The (local) compute of temperature.'
    def __init__(self, system, cutoff=0.32, angular_momentum=6):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, analysis_OrderParameter, system, cutoff, angular_momentum)

if pmi.isController :
    class OrderParameter(AnalysisBase):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espresso.analysis.OrderParameterLocal',
          pmiproperty = [ 'cutoff', 'l' ]
        )
