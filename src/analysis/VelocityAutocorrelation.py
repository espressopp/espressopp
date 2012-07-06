from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.ConfigsParticleDecomp import *
from _espresso import analysis_VelocityAutocorrelation

class VelocityAutocorrelationLocal(ConfigsParticleDecompLocal, analysis_VelocityAutocorrelation):
    'The (local) compute autocorrelation f.'
    def __init__(self, system):
      cxxinit(self, analysis_VelocityAutocorrelation, system)
      
if pmi.isController:
  class VelocityAutocorrelation(ConfigsParticleDecomp):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espresso.analysis.VelocityAutocorrelationLocal'
    )
