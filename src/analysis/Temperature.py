from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.Observable import *
from _espresso import analysis_Temperature

class TemperatureLocal(ObservableLocal, analysis_Temperature):
    'The (local) compute of temperature.'
    def __init__(self, system):
        cxxinit(self, analysis_Temperature, system)

if pmi.isController :
    class Temperature(Observable):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.analysis.TemperatureLocal'
            )
