from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.Observable import *
from _espresso import analysis_PressureTensor

class PressureTensorLocal(ObservableLocal, analysis_PressureTensor):
    'The (local) compute of pressure tensor.'
    def __init__(self, system):
        cxxinit(self, analysis_PressureTensor, system)

    def compute(self, xmin = None, xmax = None, ymin = None, ymax = None, zmin = None, zmax = None):
      if (xmin == None) or (xmax == None) or (ymin == None) or (ymax == None) or (zmin == None) or (zmax == None):
        return self.cxxclass.compute1(self)
      else:
        return self.cxxclass.compute2(self, xmin, xmax, ymin, ymax, zmin, zmax)

if pmi.isController:
  class PressureTensor(Observable):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espresso.analysis.PressureTensorLocal',
      #pmiinvoke = ['computeLocal']
      pmicall = ['compute']
    )
