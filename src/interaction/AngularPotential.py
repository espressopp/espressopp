from espresso import pmi
from espresso import toReal3DFromVector
from _espresso import interaction_AngularPotential

# Python base class for angular potentials
class AngularPotentialLocal(object):
    def computeEnergy(self, *args):
        if len(args) == 1:
            arg0 = args[0]
            if isinstance(arg0, float) or isinstance(arg0, int):
                return self.cxxclass.computeEnergy(self, arg0)
        return self.cxxclass.computeEnergy(self, toReal3DFromVector(*args))

    def computeForce(self, *args):
        return self.cxxclass.computeForce(self, toReal3DFromVector(*args))

if pmi.isController:
    class AngularPotential(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            localcall = [ 'computeForce', 'computeEnergy' ],
            pmiproperty = [ 'cutoff' ]
            )
