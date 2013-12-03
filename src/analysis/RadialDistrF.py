"""
**********************************
**espresso.analysis.RadialDistrF**
**********************************

"""
from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.Observable import *
from _espresso import analysis_RadialDistrF

class RadialDistrFLocal(ObservableLocal, analysis_RadialDistrF):
  'The (local) compute the radial distr function.'
  def __init__(self, system):
    cxxinit(self, analysis_RadialDistrF, system)
    
  def compute(self, rdfN):
    return self.cxxclass.compute(self, rdfN)
    
if pmi.isController :
  class RadialDistrF(Observable):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      pmiproperty = [ 'print_progress' ],
      pmicall = [ "compute" ],
      cls = 'espresso.analysis.RadialDistrFLocal'
    )
