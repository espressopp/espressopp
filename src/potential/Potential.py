from espresso import pmi

from _espresso import potential_PythonPotential
class PythonPotentialLocal(potential_PythonPotential):
    pass

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

        
