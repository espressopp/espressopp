from espresso import pmi
from espresso.potential.Potential import *

from _espresso import potential_PythonCentralPotential
class PythonCentralPotentialLocal(PotentialLocal, potential_PythonCentralPotential):
    pass

CentralPotentialLocal = PythonCentralPotentialLocal

# from _espresso import potential_PythonCentralPotential
# class PythonCentralPotentialLocal(CentralPotentialLocal, potential_PythonCentralPotential):
#     pass

if pmi.IS_CONTROLLER:
    class CentralPotential(Potential):
        pmiproxydefs = Potential.pmiproxydefs
