from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.Observable import *
from _espresso import analysis_Configurations

class ConfigurationsLocal(ObservableLocal, analysis_Configurations):
    'The (local) storage of configurations.'
    def __init__(self, system):
        cxxinit(self, analysis_Configurations, system)
    def getNParticles(self, stackpos):
        return self.cxxclass.getNParticles(self, stackpos)
    def getCoordinates(self, index, stackpos):
        return self.cxxclass.getCoordinates(self, index, stackpos)
    def push(self):
        return self.cxxclass.push(self)

if pmi.isController :
    class Configurations(Observable):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.analysis.ConfigurationsLocal',
            pmicall = [ "push" ],
            localcall = ["getNParticles", "getCoordinates"],
            pmiproperty = ["capacity", "size"]
            )
