from espresso.esutil import cxxinit
from espresso import pmi
from espresso.integrator.Extension import *
from _espresso import integrator_ExtAnalyze 

class ExtAnalyzeLocal(ExtensionLocal, integrator_ExtAnalyze):
    'The (local) extension analyze.'
    def __init__(self, observable, interval=1):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
           cxxinit(self, integrator_ExtAnalyze, observable, interval)

if pmi.isController :
    class ExtAnalyze(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.ExtAnalyzeLocal',
            pmicall = ['getAverage', 'getVariance']
            )
