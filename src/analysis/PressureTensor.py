from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.Observable import *
from _espresso import analysis_PressureTensor

class PressureTensorLocal(ObservableLocal, analysis_PressureTensor):
    'The (local) compute of pressure tensor.'
    def __init__(self, system):
        cxxinit(self, analysis_PressureTensor, system)

if pmi.isController :
    class PressureTensor(Observable):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.analysis.PressureTensorLocal'
            )
