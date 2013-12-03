"""
*********************************
**espresso.analysis.Temperature**
*********************************

"""
from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.AnalysisBase import *
from _espresso import analysis_Temperature

class TemperatureLocal(AnalysisBaseLocal, analysis_Temperature):
    'The (local) compute of temperature.'
    def __init__(self, system):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, analysis_Temperature, system)

if pmi.isController :
    class Temperature(AnalysisBase):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.analysis.TemperatureLocal'
            )
