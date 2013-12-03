"""
*************************************
**espresso.analysis.Autocorrelation**
*************************************

"""
from espresso.esutil import cxxinit
from espresso import pmi

from _espresso import analysis_Autocorrelation

class AutocorrelationLocal(analysis_Autocorrelation):
    'The (local) storage of configurations.'
    def __init__(self, system):
      cxxinit(self, analysis_Autocorrelation, system)
    def gather(self, value):
      return self.cxxclass.gather(self, value)
    def clear(self):
      return self.cxxclass.clear(self)
      
    def compute(self):
      return self.cxxclass.compute(self)
    
if pmi.isController:
  class Autocorrelation(object):
    """Class for parallel analysis"""
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espresso.analysis.AutocorrelationLocal',
      pmicall = [ "gather", "clear", "compute" ],
      localcall = ["__getitem__", "all"],
      pmiproperty = ["size"]
    )
