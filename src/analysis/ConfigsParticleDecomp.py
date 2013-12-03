"""
*******************************************
**espresso.analysis.ConfigsParticleDecomp**
*******************************************

"""
#from espresso.esutil import cxxinit
from espresso import pmi

from _espresso import analysis_ConfigsParticleDecomp

class ConfigsParticleDecompLocal(analysis_ConfigsParticleDecomp):
    'The (local) storage of configurations.'
    def __init__(self, system):
      cxxinit(self, analysis_ConfigsParticleDecomp, system)
    def gather(self):
      return self.cxxclass.gather(self)
    def clear(self):
      return self.cxxclass.clear(self)
    def __iter__(self):
      return self.cxxclass.all(self).__iter__()
      
    def compute(self):
      return self.cxxclass.compute(self)

if pmi.isController:
  class ConfigsParticleDecomp(object):
    """Abstract base class for parallel analysis based on particle decomposition."""
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      #cls =  'espresso.analysis.ConfigsParticleDecompLocal',
      pmicall = [ "gather", "clear", "compute" ],
      localcall = ["__getitem__", "all"],
      pmiproperty = ["size"]
    )
