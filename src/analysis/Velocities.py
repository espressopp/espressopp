from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.Observable import *
from _espresso import analysis_Velocities

class VelocitiesLocal(ObservableLocal, analysis_Velocities):
    'The (local) storage of configurations.'
    def __init__(self, system):
        cxxinit(self, analysis_Velocities, system)
    def gather(self):
        return self.cxxclass.gather(self)
    def clear(self):
        return self.cxxclass.clear(self)
    def __iter__(self):
        return self.cxxclass.all(self).__iter__()

if pmi.isController :
    class Velocities(Observable):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.analysis.VelocitiesLocal',
            pmicall = [ "gather", "clear" ],
            localcall = ["getNParticles", "getCoordinates", 
                         "__getitem__", "__iter__", "all"],
            pmiproperty = ["capacity", "size"]
            )
