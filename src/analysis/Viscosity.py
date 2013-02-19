from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.Autocorrelation import *
from _espresso import analysis_Viscosity

class ViscosityLocal(AutocorrelationLocal, analysis_Viscosity):
    'The (local) storage of configurations.'
    def __init__(self, system):
      cxxinit(self, analysis_Viscosity, system)
      
    def gather(self):
      return self.cxxclass.gather(self)
      
    def compute(self, t0, dt):
      return self.cxxclass.compute(self, t0, dt)
    
if pmi.isController:
  class Viscosity(Autocorrelation):
    """Class for parallel analysis"""
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espresso.analysis.ViscosityLocal',
      pmicall = [ "gather", "compute" ]
    )
