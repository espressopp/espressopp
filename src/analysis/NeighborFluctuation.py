"""
*****************************************
**espresso.analysis.NeighborFluctuation**
*****************************************

"""
from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.Observable import *
from _espresso import analysis_NeighborFluctuation

class NeighborFluctuationLocal(ObservableLocal, analysis_NeighborFluctuation):
    'The (local) compute of the neighbor fluctuations (<n^2>-<n>^2) in the number of particles found in a sphere of radius d around particle i.'
    def __init__(self, system, radius):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, analysis_NeighborFluctuation, system, radius)

if pmi.isController :
    class NeighborFluctuation(Observable):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.analysis.NeighborFluctuationLocal'
            )
