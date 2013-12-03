"""
******************************
**espresso.analysis.Pressure**
******************************

"""
from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.Observable import *
from _espresso import analysis_Pressure

class PressureLocal(ObservableLocal, analysis_Pressure):
    'The (local) compute of pressure.'
    def __init__(self, system):
        cxxinit(self, analysis_Pressure, system)

if pmi.isController :
    class Pressure(Observable):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.analysis.PressureLocal'
            )
