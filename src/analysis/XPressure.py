from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.Observable import *
from _espresso import analysis_XPressure

class XPressureLocal(ObservableLocal, analysis_XPressure):
  'The (local) compute the pressure profile in x direction.'
  def __init__(self, system):
    cxxinit(self, analysis_XPressure, system)
    
  def compute(self, N):
    return self.cxxclass.compute(self, N)
    
if pmi.isController :
  class XPressure(Observable):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      pmicall = [ "compute" ],
      cls = 'espresso.analysis.XPressureLocal'
    )
