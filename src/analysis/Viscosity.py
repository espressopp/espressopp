from espresso import pmi
from espresso.esutil import cxxinit

from _espresso import analysis_Viscosity

from espresso.analysis.Autocorrelation import *

class ViscosityLocal(AutocorrelationLocal, analysis_Viscosity):  
    'The (local) storage of configurations.'
    def __init__(self, system):
      if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        cxxinit(self, analysis_Viscosity, system)
      
    def gather(self):
      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        return self.cxxclass.gather(self)
      
    def compute(self, t0, dt, T):
      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        return self.cxxclass.compute(self, t0, dt, T)
    
if pmi.isController:
  class Viscosity(Autocorrelation):
    """Class for parallel analysis"""
    #__metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espresso.analysis.ViscosityLocal',
      pmicall = [ 'gather', 'compute' ]
    )
    def __init__(self, system):
      self.pmiinit(system)
