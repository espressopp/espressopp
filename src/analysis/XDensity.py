"""
******************************
**espresso.analysis.XDensity**
******************************

"""
from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.Observable import *
from _espresso import analysis_XDensity

class XDensityLocal(ObservableLocal, analysis_XDensity):
  'The (local) compute the density profile in x direction.'
  def __init__(self, system):
    cxxinit(self, analysis_XDensity, system)
    
  def compute(self, rdfN):
    return self.cxxclass.compute(self, rdfN)
    
if pmi.isController :
  class XDensity(Observable):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      pmicall = [ "compute" ],
      cls = 'espresso.analysis.XDensityLocal'
    )
