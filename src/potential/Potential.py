from espresso import pmi
from espresso import toReal3DFromVector

from _espresso import potential_PythonPotential

class PythonPotentialLocal(potential_PythonPotential):
    def getCutoffSqr(self):
        pass

    def computeForce(self, *args):
        """Override this method to compute the force for a given distance.
        
        It should at least be able to handle a Real3D distance input.
        """
        pass

    def computeEnergy(self, *args):
        """Override this method to compute the energy at a given distance.
        
        It should at least be able to handle a Real3D distance input.
        """
        pass

# Python base class for potentials
class PotentialLocal(object):
    def computeForce(self, *args):
        return \
            self.cxxclass.computeForce(self,
                toReal3DFromVector(*args))
    def computeEnergy(self, *args):
        return \
            self.cxxclass.computeEnergy(self,
                toReal3DFromVector(*args))

# # this class is the base of all potentials implemented in Python
# from _espresso import potential_PythonPotential
# class PythonPotentialLocal(PotentialLocal, potential_PythonPotential):
#     pass

if pmi.IS_CONTROLLER:
    class Potential:
        __metaclass__ = pmi.Proxy
        pmiproxydefs = {
            'localcall' : [ 'computeForce', 'computeEnergy' ]
            }

        
