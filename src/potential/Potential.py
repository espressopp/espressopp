from espresso import pmi
from espresso import toReal3DFromVector

from _espresso import potential_PythonPotential
class PythonPotentialLocal(object):
    def computeForce(self, *args):
        return \
            self.cxxclass.computeForce(self,
                toReal3DFromVector(*args))
    def computeEnergy(self, *args):
        return \
            self.cxxclass.computeEnergy(self,
                toReal3DFromVector(*args))

PotentialLocal = PythonPotentialLocal

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

        
