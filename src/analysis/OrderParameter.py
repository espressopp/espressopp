"""
************************************
**espresso.analysis.OrderParameter**
************************************

"""
from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.AnalysisBase import *
from _espresso import analysis_OrderParameter

class OrderParameterLocal(AnalysisBaseLocal, analysis_OrderParameter):
    'The (local) compute of temperature.'
    def __init__(self, system, cutoff, angular_momentum=6,
                      do_cluster_analysis=False, include_surface_particles=False,
                      ql_low=-1.0, ql_high=1.0):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, analysis_OrderParameter, system, cutoff, angular_momentum,
                      do_cluster_analysis, include_surface_particles,
                      ql_low, ql_high)

if pmi.isController :
    class OrderParameter(AnalysisBase):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espresso.analysis.OrderParameterLocal',
          pmiproperty = [ 'cutoff', 'l' ]
        )
