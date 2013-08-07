from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.Observable import *
from _espresso import analysis_StaticStructF

class StaticStructFLocal(ObservableLocal, analysis_StaticStructF):
  'The (local) compute the radial distr function.'
  def __init__(self, system):
    cxxinit(self, analysis_StaticStructF, system)
    
  def compute(self, rdfN):
    return self.cxxclass.compute(self, rdfN)
    
if pmi.isController :
  class StaticStructF(Observable):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      pmicall = [ "compute" ],
      cls = 'espresso.analysis.StaticStructFLocal'
    )
