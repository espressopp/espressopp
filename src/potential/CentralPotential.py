from espresso import pmi
from espresso.potential.Potential import *

from _espresso import potential_PythonCentralPotential

class PythonCentralPotentialLocal(PotentialLocal, potential_PythonCentralPotential):
    def getCutoffSqr(self):
        pass

    def computeForce(self, *args):
        """Override this method to compute the force for a given distance.
        
        It should at least be able to handle a Real3D distance input.
        """
        pass

    def computeEnergySqr(self, distSqr):
        """Override this method to compute the energy at a given distance.
        
        It should at least be able to handle a Real3D distance input.
        """
        pass

    
class CentralPotentialLocal(PotentialLocal):
    def computeEnergy(self, *args):
        if len(args) == 1:
            arg0 = args[0]
            if type(arg0) is float or type(arg0) is int:
                return self.cxxclass.computeEnergy(self, arg0)
        return PotentialLocal.computeEnergy(self, *args)

CentralPotentialLocal = PythonCentralPotentialLocal

# from _espresso import potential_PythonCentralPotential
# class PythonCentralPotentialLocal(CentralPotentialLocal, potential_PythonCentralPotential):
#     pass

if pmi.IS_CONTROLLER:
    class CentralPotential(Potential):
        pmiproxydefs = Potential.pmiproxydefs
