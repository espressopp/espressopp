from espresso import pmi
from espresso import toReal3DFromVector

from _espresso import interaction_PotentialVSpherePair

# Python base class for potentials
class PotentialVSpherePairLocal(object):
    def computeEnergy(self, *args):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if len(args) == 1:
                arg0 = args[0]
                if isinstance(arg0, float) or isinstance(arg0, int):
                    return self.cxxclass.computeEnergy(self, arg0)
            return self.cxxclass.computeEnergy(self, toReal3DFromVector(*args))

    def computeForce(self, *args):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if len(args) == 1:
                arg0 = args[0]
                if isinstance(arg0, float) or isinstance(arg0, int):
                    newargs = [arg0, 0, 0]
                    return self.cxxclass.computeForce(self, toReal3DFromVector(*newargs))[0]
            return self.cxxclass.computeForce(self, toReal3DFromVector(*args))

    def _setShift(self, shift="auto"):
        if (shift == "auto"):
            if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
                self.cxxclass.setAutoShift(self)
        else:
            if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
                self.cxxclass.shift.fset(self, shift)

    def _getShift(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.shift.fget(self)

    shift = property(_getShift, _setShift)

if pmi.isController:
    class PotentialVSpherePair(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            localcall = [ 'computeForce', 'computeEnergy' ],
            pmiproperty = ['cutoff', 'shift']
            )
        

        
# class PythonPotentialLocal(potential_PythonPotential):
#     def getCutoffSqr(self):
#         pass

#     def computeForce(self, *args):
#         """Override this method to compute the force for a given distance.
        
#         It should at least be able to handle a Real3D distance input.
#         """
#         pass

#     def computeEnergy(self, *args):
#         """Override this method to compute the energy at a given distance.
        
#         It should at least be able to handle a Real3D distance input.
#         """
#         pass
