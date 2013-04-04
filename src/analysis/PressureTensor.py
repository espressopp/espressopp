from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.Observable import *
from _espresso import analysis_PressureTensor

class PressureTensorLocal(ObservableLocal, analysis_PressureTensor):
    'The (local) compute of pressure tensor.'
    def __init__(self, system):
        cxxinit(self, analysis_PressureTensor, system)

    '''
    def compute(self, z = None, dz = None):
      if (z == None) or (dz == None):
        ret = self.cxxclass.compute1(self)
        print 'pressure tensor return:', ret
        exit(0)
        return ret
      elif(isinstance(z, int) and z!=0):
        return self.cxxclass.compute2(self, z, dz)
      else:
        return self.cxxclass.compute3(self, z, dz)
    '''
        

if pmi.isController:
  class PressureTensor(Observable):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espresso.analysis.PressureTensorLocal',
      #pmiinvoke = ['computeLocal']
      #pmicall = ['compute']
    )
